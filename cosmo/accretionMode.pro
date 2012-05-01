; accretionMode.pro
; gas accretion project - past substructure history of gas elements
; dnelson may.2012

; -----------------------------------------------------------------------------------------------------
; accretionMode(): for eaching gas particle/tracer with a recorded accretion time, starting at some 
;                  redshift, track backwards in time with respect to the tracked parent halos (using 
;                  mergerTree) and classify the mode of accretion by considering its group membership
;                  at the outer bracketing snapshot just prior to accretion (r>~r_vir) as:
;                  
;  1. smooth  = in the (primary) parent of the original halo or not in any subgroup
;  2. sclumpy = in any subgroup other than the (primary) parent halo that is a secondary subgroup ("small")
;  3. bclumpy = in any subgroup other than the (primary) parent halo that is a primary subgroup ("big")
; -----------------------------------------------------------------------------------------------------

function accretionMode, sP=sP

  forward_function cosmoTracerChildren, cosmoTracerVelParents, accretionTimes
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; first, walk back through the merger tree and find primary subhalos with good parent histories
  mt = mergerTreeSubset(sP=sP,/verbose)
  at = accretionTimes(sP=sP)

  ; set saveFilename and check for existence
  saveTag = ''
  if sP.trMCPerCell eq -1 then saveTag = '.trVel'
  if sP.trMCPerCell gt 0  then saveTag = '.trMC'
  if sP.trMCPerCell eq 0  then saveTag = '.SPH'
  
  saveFilename = sP.derivPath + 'accMode'+saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(mt.maxSnap)+'-'+str(mt.minSnap)+'.sav'  
  
  if file_test(saveFilename) then begin
    restore, saveFilename
    return, r
  endif
  
  ; select those particles/tracers with recorded accretion times
  gal_w_at  = where(at.AccTime_gal[0,*] ne -1,count_gal)
  gmem_w_at = where(at.AccTime_gmem[0,*] ne -1,count_gmem)
  
  ; load galaxy/group member catalogs at zMin for gas ids to search for
  galcat = galaxyCat(sP=sP)  
  
  count_bracketed = 0UL
  
  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    print,'Calculating new accretion mode using ( SPH Particles ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
      
    ; only find accretion mode for those gas particles with accretion times
    galcatSub = { gal  : mt.galcatSub.gal[gal_w_at]    ,$
                  gmem : mt.galcatSub.gmem[gmem_w_at]   }
      
    ; convert the accretion time (scale factor) into the outer bracketing snapshot number for each member
    bracketSnap = { gal  : intarr(n_elements(galcatSub.gal))   ,$
                    gmem : intarr(n_elements(galcatSub.gmem))   }
                    
    snapTimes = snapNumToRedshift(sP=sP,/all,/time)
    snapTimes = snapTimes[where(snapTimes ne -1)]

    bracketSnap.gal  = value_locate(snapTimes,at.accTime_gal[0,gal_w_at])
    bracketSnap.gmem = value_locate(snapTimes,at.accTime_gmem[0,gmem_w_at])
    
    w = where(bracketSnap.gal eq -1 or bracketSnap.gal eq n_elements(snapTimes)-1,count)
    if count then message,'Error: Bad gal bracketing.'
    w = where(bracketSnap.gmem eq -1 or bracketSnap.gmem eq n_elements(snapTimes)-1,count)
    if count then message,'Error: Bad gmem bracketing.'
    
    ; list of parent subgroup IDs at starting snapshot for each gas particle
    gcIndOrig = galCatRepParentIDs(galcat=galcat,gcIDList=mt.galcatIDList)
    
    origParIDs = { gal  : gcIndOrig.gal[gal_w_at]   ,$
                   gmem : gcIndOrig.gmem[gmem_w_at]  }
    
    ; store the main arrays as a structure so we can write them directly
    r = {accMode_gal       : intarr(n_elements(galcatSub.gal))-1  ,$
         accMode_gmem      : intarr(n_elements(galcatSub.gmem))-1  }
         
    ; debugging
    rMask = { gal  : intarr(n_elements(galcatSub.gal))    ,$
              gmem : intarr(n_elements(galcatSub.gmem))    }
    
    for m=mt.maxSnap,mt.minSnap,-1 do begin
      sP.snap = m
      ; load local group catalog
      h = loadSnapshotHeader(sP=sP)
      gc = loadGroupCat(sP=sP,/readIDs)
      
      ; select those gas particles whose accretion modes are determined at this snapshot
      w_gal  = where(bracketSnap.gal eq sP.snap,count_gal)
      w_gmem = where(bracketSnap.gmem eq sP.snap,count_gmem)
      
      ;print,'  ['+str(m)+'] bracketed gal: '+str(count_gal)+' gmem: '+str(count_gmem)
      
      count_smooth  = 0UL
      count_sclumpy = 0UL
      count_bclumpy = 0UL
      
      if count_gal gt 0 then begin
        ; global match against subgroup ID list (this includes fuzz and unbound gas)
        match,galcat.galaxyIDs[galcatSub.gal[w_gal]],gc.IDs,galcat_ind,gc_gal_ind,count=countGal,/sort
        
        ; those that don't match, assign to category 1. smooth
        if countGal lt count_gal then begin
          ; if none matched, assign all the bracketed particles to smooth
          if countGal eq 0 then begin
            r.accMode_gal[w_gal] = 1
            rMask.gal[w_gal] += 1
            count_smooth += count_gal
          endif else begin
            ; if some matched, assign the complement of those to smooth
            all = bytarr(count_gal)
            all[galcat_ind] = 1B
            wNotInSG = where(all eq 0B, ncomp)
            
            ;if ncomp gt 0 then begin
              r.accMode_gal[w_gal[wNotInSG]] = 1
              rMask.gal[w_gal[wNotInSG]] += 1 ;debug
              count_smooth += ncomp
            ;endif
          endelse
        endif
        
        ; those that match, calculate subgroup ID they belong to
        if countGal gt 0 then begin
          ; calculate parent subgroup ID each gas particle belongs to
          parIDs = value_locate(gc.subgroupOffsetType[partTypeNum('gas'),*],gc_gal_ind)
          
          ; exclude those that aren't actually within the offsetLenType range
          diff = gc_gal_ind - gc.subgroupOffsetType[partTypeNum('gas'),parIDs]
          wInSG = where(diff lt gc.subgroupLenType[partTypeNum('gas'),parIDs],countInSG,comp=wc,ncomp=ncomp)
          
          ;print,'  ['+str(m)+'] sg gal bounds in: '+str(count)+' out: '+str(ncomp)
          
          ; those failing the subgroup bound are smooth
          if ncomp gt 0 then begin
            rInds = w_gal[galcat_ind[wc]]
            r.accMode_gal[rInds] = 1
            rMask.gal[rInds] += 1 ;debug
            count_smooth += ncomp
          endif
          
          if countInSG gt 0 then begin
            ; if their subgroup ID is the traced parent, assign to category 1. smooth
            rInds = w_gal[galcat_ind[wInSG]]
            w0 = where(parIDs[wInSG] eq origParIDs.gal[galcat_ind[wInSG]],count,comp=wc,ncomp=ncomp)
            if count gt 0 then begin
              r.accMode_gal[rInds[w0]] = 1
              count_smooth += count
              rMask.gal[rInds[w0]] += 1 ;debug
            endif
            
            ; only separate into sclumpy/bclumpy for the non-smooth remainder
            if ncomp gt 0 then begin
              wInSG = wInSG[wc]
              rInds = rInds[wc]
            endif
            
            ; make ID lists of primary,secondary subgroups
            priSGids = gcIDList(gc=gc,select='pri')
            secSGids = gcIDList(gc=gc,select='sec')
            
            ; if their subgroup ID is not the parent but is a primary sg, assign to category 3. bclumpy
            qInds = value_locate(priSGids,parIDs[wInSG])
            w1 = where(priSGids[qInds] eq parIDs[wInSG],count,comp=wc,ncomp=ncomp)
            if count gt 0 then begin
              r.accMode_gal[rInds[w1]] = 3
              count_bclumpy += count
              rMask.gal[rInds[w1]] += 1 ; debug
            endif
            
            ; only search for secondary parents for the non-primary remainder (this is mainly a doublecheck
            ; since we have already excluded every other possibility)
            if ncomp gt 0 then begin
              wInSG = wInSG[wc]
              rInds = rInds[wc]
            
              qInds = value_locate(secSGids,parIDs[wInSG])
              w2 = where(secSGids[qInds] eq parIDs[wInSG],count)
              if count gt 0 then begin
                r.accMode_gal[rInds[w2]] = 2
                count_sclumpy += count
                rMask.gal[rInds[w2]] += 1 ; debug
              endif
            endif
          endif ; countInSG
        endif ; countGal
        
        ; check that all bracketed gal particles found accretion modes
        if min(r.accMode_gal[w_gal]) eq 0 then message,'Error: Not all gal done.'
        if (count_smooth+count_sclumpy+count_bclumpy) ne count_gal then message,'Error: Gal totals.'
        if max(rMask.gal) gt 1 then message,'Error: gal mask.'
      endif ; count_gal
        
      count_smooth  = 0UL
      count_sclumpy = 0UL
      count_bclumpy = 0UL
        
      if count_gmem gt 0 then begin
        ; global match against subgroup ID list (this includes fuzz and unbound gas)
        match,galcat.groupmemIDs[galcatSub.gmem[w_gmem]],gc.IDs,galcat_ind,gc_gmem_ind,count=countGmem,/sort
        
        ; those that don't match, assign to category 1. smooth
        if countGmem lt count_gmem then begin
          ; if none matched, assign all the bracketed particles to smooth
          if countGmem eq 0 then begin
            r.accMode_gmem[w_gmem] = 1
            rMask.gmem[w_gmem] += 1
            count_smooth += count_gmem
          endif else begin
            ; if some matched, assign the complement of those to smooth
            all = bytarr(count_gmem)
            all[galcat_ind] = 1B
            wNotInSG = where(all eq 0B, ncomp)
    
            ;if ncomp gt 0 then begin
              r.accMode_gmem[w_gmem[wNotInSG]] = 1
              rMask.gmem[w_gmem[wNotInSG]] += 1 ;debug
              count_smooth += ncomp
            ;endif
          endelse
        endif
        
        ; those that match, calculate subgroup ID they belong to
        if countGmem gt 0 then begin
          ; calculate parent subgroup ID each gas particle belongs to
          parIDs = value_locate(gc.subgroupOffsetType[partTypeNum('gas'),*],gc_gmem_ind)
          
          ; exclude those that aren't actually within the offsetLenType range
          diff = gc_gmem_ind - gc.subgroupOffsetType[partTypeNum('gas'),parIDs]
          wInSG = where(diff lt gc.subgroupLenType[partTypeNum('gas'),parIDs],countInSG,comp=wc,ncomp=ncomp)

          ;print,'  ['+str(m)+'] sg gmem bounds in: '+str(count)+' out: '+str(ncomp)
          
          ; those failing the subgroup bound are smooth
          if ncomp gt 0 then begin
            rInds = w_gmem[galcat_ind[wc]]
            r.accMode_gmem[rInds] = 1
            rMask.gmem[rInds] += 1 ;debug
            count_smooth += ncomp
          endif
          
          if countInSG gt 0 then begin
            ; if their subgroup ID is the traced parent, assign to category 1. smooth
            rInds = w_gmem[galcat_ind[wInSG]]
            w0 = where(parIDs[wInSG] eq origParIDs.gmem[galcat_ind[wInSG]],count,comp=wc,ncomp=ncomp)
            if count gt 0 then begin
              r.accMode_gmem[rInds[w0]] = 1
              count_smooth += count
              rMask.gmem[rInds[w0]] += 1 ;debug
            endif
            
            ; only separate into sclumpy/bclumpy for the non-smooth remainder
            if ncomp gt 0 then begin
              wInSG = wInSG[wc]
              rInds = rInds[wc]
            endif
            
            ; make ID lists of primary,secondary subgroups
            priSGids = gcIDList(gc=gc,select='pri')
            secSGids = gcIDList(gc=gc,select='sec')
            
            ; if their subgroup ID is not the parent but is a primary sg, assign to category 3. bclumpy
            qInds = value_locate(priSGids,parIDs[wInSG])
            w1 = where(priSGids[qInds] eq parIDs[wInSG],count,comp=wc,ncomp=ncomp)
            if count gt 0 then begin
              r.accMode_gmem[rInds[w1]] = 3
              count_bclumpy += count
              rMask.gmem[rInds[w1]] += 1 ; debug
            endif
            
            ; only search for secondary parents for the non-primary remainder (this is mainly a doublecheck
            ; since we have already excluded every other possibility)
            if ncomp gt 0 then begin
              wInSG = wInSG[wc]
              rInds = rInds[wc]
  
              qInds = value_locate(secSGids,parIDs[wInSG])
              w2 = where(secSGids[qInds] eq parIDs[wInSG],count)
              if count gt 0 then begin
                r.accMode_gmem[rInds[w2]] = 2
                count_sclumpy += count
                rMask.gmem[rInds[w2]] += 1 ; debug
              endif
            endif
          endif ; countInSG
        endif ; countGmem
        
        ; check that all bracketed gmem particles found accretion modes
        if min(r.accMode_gmem[w_gmem]) eq 0 then message,'Error: Not all gmem done.'
        if (count_smooth+count_sclumpy+count_bclumpy) ne count_gmem then message,'Error: Gmem totals.'
        if max(rMask.gmem) gt 1 then message,'Error: gmem mask.'
      endif ; count_gal
          
      if m ne mt.minSnap then begin
        ; load mergerTree and move to Parent
        Parent = mergerTree(sP=sP)
        
        origParIDs.gal  = Parent[origParIDs.gal] ; change to parent IDs
        origParIDs.gmem = Parent[origParIDs.gmem]
        
        ; sanity check no IDs are -1 (we should only be processing the mergerTreeSubset)
        if min(origParIDs.gal)  lt 0 then message,'Error: Bad gal parent.'
        if min(origParIDs.gmem) lt 0 then message,'Error: Bad gmem parent.'
      endif
      
      count_bracketed += count_gal+count_gmem
      cur_frac_found = float(count_bracketed) / (n_elements(gal_w_at)+n_elements(gmem_w_at))
      print,' ['+string(m,format='(i3)')+'] counts smooth: '+string(count_smooth,format='(i5)')+$
        ' sclumpy: '+string(count_sclumpy,format='(i5)')+' bclumpy: '+$
        string(count_bclumpy,format='(i5)')+'  ('+string(cur_frac_found*100,format='(f4.1)')+'%)'
      
      ; free some memory for next load
      Parent = !NULL
      all    = !NULL
      w_gal  = !NULL
      w_gmem = !NULL
    endfor
    
    ; verify we found an accretion mode for every gas particle
    if min(rMask.gal)  lt 1 then message,'Error: Not all gal found.'
    if min(rMask.gmem) lt 1 then message,'Error: Not all gmem found.'
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
  endif  
  
  ; MONTE CARLO TRACERS CASE - for all gas cells, track back all child tracers
  ; ------------------------
  if sP.trMCPerCell gt 0 then begin
    print,'Calculating new accretion mode using ( TracerMC ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
      
    ; make gcIndOrigTr and ids of sorted tracer children
    gcIndOrigTr = mergerTreeRepParentIDs(mt=mt,galcat=galcat,sP=sP,$
                    trids_gal=galcat_gal_trids, trids_gmem=galcat_gmem_trids)
    
    origParIDs = { gal  : gcIndOrigTr.gal[gal_w_at]   ,$
                   gmem : gcIndOrigTr.gmem[gmem_w_at]  }
                   
    gcIndOrigTr = !NULL
    galcat = !NULL
    
    ; only find accretion mode for those tracers with accretion times
    galcatSub = { gal  : galcat_gal_trids[gal_w_at]    ,$
                  gmem : galcat_gmem_trids[gmem_w_at]   }
                  
    galcat_gal_trids  = !NULL
    galcat_gmem_trids = !NULL
      
    ; convert the accretion time (scale factor) into the outer bracketing snapshot number for each member
    bracketSnap = { gal  : intarr(n_elements(galcatSub.gal))   ,$
                    gmem : intarr(n_elements(galcatSub.gmem))   }
                    
    snapTimes = snapNumToRedshift(sP=sP,/all,/time)
    snapTimes = snapTimes[where(snapTimes ne -1)]

    bracketSnap.gal  = value_locate(snapTimes,at.accTime_gal[0,gal_w_at])
    bracketSnap.gmem = value_locate(snapTimes,at.accTime_gmem[0,gmem_w_at])
    
    w = where(bracketSnap.gal eq -1 or bracketSnap.gal eq n_elements(snapTimes)-1,count)
    if count then message,'Error: Bad gal bracketing.'
    w = where(bracketSnap.gmem eq -1 or bracketSnap.gmem eq n_elements(snapTimes)-1,count)
    if count then message,'Error: Bad gmem bracketing.'
    
    ; store the main arrays as a structure so we can write them directly
    r = {accMode_gal       : intarr(n_elements(galcatSub.gal))-1  ,$
         accMode_gmem      : intarr(n_elements(galcatSub.gmem))-1  }
         
    ; debugging
    rMask = { gal  : intarr(n_elements(galcatSub.gal))    ,$
              gmem : intarr(n_elements(galcatSub.gmem))    }
    
    for m=mt.maxSnap,mt.minSnap,-1 do begin
      sP.snap = m
      ; load local group catalog
      h = loadSnapshotHeader(sP=sP)
      gc = loadGroupCat(sP=sP,/readIDs)
      
      ; differing from the sph case, we need to map gc.IDs (all) to a child tracerMC list for each snapshot
      groupcat_trIDs = cosmoTracerChildren(sP=sP, /getIDs, gasIDs=gc.IDs, child_counts=groupcat_cc)
    
      ; and construct the subgroupLen and subgroupOffset arrays for the tracerMCs
      subgroupOffsetTr = ulonarr(gc.nSubgroupsTot)
      subgroupLenTr    = ulonarr(gc.nSubgroupsTot)
      
      globalTrCount = 0UL
      parPT = partTypeNum('gas')
      
      for SubNr=0L, gc.nSubgroupsTot-1 do begin
        ; set offset to current marker position
        subgroupOffsetTr[SubNr] = globalTrCount
        
        if gc.subgroupLenType[parPT,SubNr] gt 0 then begin
          ; indices into gc.IDs for this subgroup
          id_inds = lindgen(gc.subgroupLenType[parPT,SubNr]) + gc.subgroupOffsetType[parPT,SubNr]
          
          ; total number of tracers in this subgroup
          subgroupLenTr[SubNr] = total(groupcat_cc[id_inds],/int)

          globalTrCount += subgroupLenTr[SubNr]
        endif
      endfor
      
      ; select those gas particles whose accretion modes are determined at this snapshot
      w_gal  = where(bracketSnap.gal eq sP.snap,count_gal)
      w_gmem = where(bracketSnap.gmem eq sP.snap,count_gmem)
      
      ;print,'  ['+str(m)+'] bracketed gal: '+str(count_gal)+' gmem: '+str(count_gmem)
      
      count_smooth  = 0UL
      count_sclumpy = 0UL
      count_bclumpy = 0UL
      
      if count_gal gt 0 then begin
        ; global match against subgroup ID list (this includes fuzz and unbound gas)
        match,galcatSub.gal[w_gal],groupcat_trIDs,galcat_ind,gc_gal_ind,count=countGal,/sort
        
        ; those that don't match, assign to category 1. smooth
        if countGal lt count_gal then begin
          ; if none matched, assign all the bracketed particles to smooth
          if countGal eq 0 then begin
            r.accMode_gal[w_gal] = 1
            rMask.gal[w_gal] += 1
            count_smooth += count_gal
          endif else begin
            ; if some matched, assign the complement of those to smooth
            all = bytarr(count_gal)
            all[galcat_ind] = 1B
            wNotInSG = where(all eq 0B, ncomp)
            
            ;if ncomp gt 0 then begin
              r.accMode_gal[w_gal[wNotInSG]] = 1
              rMask.gal[w_gal[wNotInSG]] += 1 ;debug
              count_smooth += ncomp
            ;endif
          endelse
        endif
        
        ; those that match, calculate subgroup ID they belong to
        if countGal gt 0 then begin
          ; calculate parent subgroup ID each gas particle belongs to
          parIDs = value_locate(subgroupOffsetTr,gc_gal_ind)
          
          ; exclude those that aren't actually within the offsetLenType range
          diff = gc_gal_ind - subgroupOffsetTr[parIDs]
          wInSG = where(diff lt subgroupLenTr[parIDs],countInSG,comp=wc,ncomp=ncomp)
          
          ;print,'  ['+str(m)+'] sg gal bounds in: '+str(count)+' out: '+str(ncomp)
          
          ; those failing the subgroup bound are smooth
          if ncomp gt 0 then begin
            rInds = w_gal[galcat_ind[wc]]
            r.accMode_gal[rInds] = 1
            rMask.gal[rInds] += 1 ;debug
            count_smooth += ncomp
          endif
          
          if countInSG gt 0 then begin
            ; if their subgroup ID is the traced parent, assign to category 1. smooth
            rInds = w_gal[galcat_ind[wInSG]]
            w0 = where(parIDs[wInSG] eq origParIDs.gal[galcat_ind[wInSG]],count,comp=wc,ncomp=ncomp)
            if count gt 0 then begin
              r.accMode_gal[rInds[w0]] = 1
              count_smooth += count
              rMask.gal[rInds[w0]] += 1 ;debug
            endif
            
            ; only separate into sclumpy/bclumpy for the non-smooth remainder
            if ncomp gt 0 then begin
              wInSG = wInSG[wc]
              rInds = rInds[wc]
            endif
            
            ; make ID lists of primary,secondary subgroups
            priSGids = gcIDList(gc=gc,select='pri')
            secSGids = gcIDList(gc=gc,select='sec')
            
            ; if their subgroup ID is not the parent but is a primary sg, assign to category 3. bclumpy
            qInds = value_locate(priSGids,parIDs[wInSG])
            w1 = where(priSGids[qInds] eq parIDs[wInSG],count,comp=wc,ncomp=ncomp)
            if count gt 0 then begin
              r.accMode_gal[rInds[w1]] = 3
              count_bclumpy += count
              rMask.gal[rInds[w1]] += 1 ; debug
            endif
            
            ; only search for secondary parents for the non-primary remainder (this is mainly a doublecheck
            ; since we have already excluded every other possibility)
            if ncomp gt 0 then begin
              wInSG = wInSG[wc]
              rInds = rInds[wc]
            
              qInds = value_locate(secSGids,parIDs[wInSG])
              w2 = where(secSGids[qInds] eq parIDs[wInSG],count)
              if count gt 0 then begin
                r.accMode_gal[rInds[w2]] = 2
                count_sclumpy += count
                rMask.gal[rInds[w2]] += 1 ; debug
              endif
            endif
          endif ; countInSG
        endif ; countGal
        
        ; check that all bracketed gal particles found accretion modes
        if min(r.accMode_gal[w_gal]) eq 0 then message,'Error: Not all gal done.'
        if (count_smooth+count_sclumpy+count_bclumpy) ne count_gal then message,'Error: Gal totals.'
        if max(rMask.gal) gt 1 then message,'Error: gal mask.'
      endif ; count_gal
    
      count_smooth  = 0UL
      count_sclumpy = 0UL
      count_bclumpy = 0UL
        
      if count_gmem gt 0 then begin
        ; global match against subgroup ID list (this includes fuzz and unbound gas)
        match,galcatSub.gmem[w_gmem],groupcat_trIDs,galcat_ind,gc_gmem_ind,count=countGmem,/sort
        
        ; those that don't match, assign to category 1. smooth
        if countGmem lt count_gmem then begin
          ; if none matched, assign all the bracketed particles to smooth
          if countGmem eq 0 then begin
            r.accMode_gmem[w_gmem] = 1
            rMask.gmem[w_gmem] += 1
            count_smooth += count_gmem
          endif else begin
            ; if some matched, assign the complement of those to smooth
            all = bytarr(count_gmem)
            all[galcat_ind] = 1B
            wNotInSG = where(all eq 0B, ncomp)
    
            ;if ncomp gt 0 then begin
              r.accMode_gmem[w_gmem[wNotInSG]] = 1
              rMask.gmem[w_gmem[wNotInSG]] += 1 ;debug
              count_smooth += ncomp
            ;endif
          endelse
        endif
        
        ; those that match, calculate subgroup ID they belong to
        if countGmem gt 0 then begin
          ; calculate parent subgroup ID each gas particle belongs to
          parIDs = value_locate(subgroupOffsetTr,gc_gmem_ind)
          
          ; exclude those that aren't actually within the offsetLenType range
          diff = gc_gmem_ind - subgroupOffsetTr[parIDs]
          wInSG = where(diff lt subgroupLenTr[parIDs],countInSG,comp=wc,ncomp=ncomp)

          ;print,'  ['+str(m)+'] sg gmem bounds in: '+str(count)+' out: '+str(ncomp)
          
          ; those failing the subgroup bound are smooth
          if ncomp gt 0 then begin
            rInds = w_gmem[galcat_ind[wc]]
            r.accMode_gmem[rInds] = 1
            rMask.gmem[rInds] += 1 ;debug
            count_smooth += ncomp
          endif
          
          if countInSG gt 0 then begin
            ; if their subgroup ID is the traced parent, assign to category 1. smooth
            rInds = w_gmem[galcat_ind[wInSG]]
            w0 = where(parIDs[wInSG] eq origParIDs.gmem[galcat_ind[wInSG]],count,comp=wc,ncomp=ncomp)
            if count gt 0 then begin
              r.accMode_gmem[rInds[w0]] = 1
              count_smooth += count
              rMask.gmem[rInds[w0]] += 1 ;debug
            endif
            
            ; only separate into sclumpy/bclumpy for the non-smooth remainder
            if ncomp gt 0 then begin
              wInSG = wInSG[wc]
              rInds = rInds[wc]
            endif
            
            ; make ID lists of primary,secondary subgroups
            priSGids = gcIDList(gc=gc,select='pri')
            secSGids = gcIDList(gc=gc,select='sec')
            
            ; if their subgroup ID is not the parent but is a primary sg, assign to category 3. bclumpy
            qInds = value_locate(priSGids,parIDs[wInSG])
            w1 = where(priSGids[qInds] eq parIDs[wInSG],count,comp=wc,ncomp=ncomp)
            if count gt 0 then begin
              r.accMode_gmem[rInds[w1]] = 3
              count_bclumpy += count
              rMask.gmem[rInds[w1]] += 1 ; debug
            endif
            
            ; only search for secondary parents for the non-primary remainder (this is mainly a doublecheck
            ; since we have already excluded every other possibility)
            if ncomp gt 0 then begin
              wInSG = wInSG[wc]
              rInds = rInds[wc]
  
              qInds = value_locate(secSGids,parIDs[wInSG])
              w2 = where(secSGids[qInds] eq parIDs[wInSG],count)
              if count gt 0 then begin
                r.accMode_gmem[rInds[w2]] = 2
                count_sclumpy += count
                rMask.gmem[rInds[w2]] += 1 ; debug
              endif
            endif
          endif ; countInSG
        endif ; countGmem
        
        ; check that all bracketed gmem particles found accretion modes
        if min(r.accMode_gmem[w_gmem]) eq 0 then message,'Error: Not all gmem done.'
        if (count_smooth+count_sclumpy+count_bclumpy) ne count_gmem then message,'Error: Gmem totals.'
        if max(rMask.gmem) gt 1 then message,'Error: gmem mask.'
      endif ; count_gal
          
      if m ne mt.minSnap then begin
        ; load mergerTree and move to Parent
        Parent = mergerTree(sP=sP)
        
        origParIDs.gal  = Parent[origParIDs.gal] ; change to parent IDs
        origParIDs.gmem = Parent[origParIDs.gmem]
        
        ; sanity check no IDs are -1 (we should only be processing the mergerTreeSubset)
        if min(origParIDs.gal)  lt 0 then message,'Error: Bad gal parent.'
        if min(origParIDs.gmem) lt 0 then message,'Error: Bad gmem parent.'
      endif
      
      count_bracketed += count_gal+count_gmem
      cur_frac_found = float(count_bracketed) / (n_elements(gal_w_at)+n_elements(gmem_w_at))
      print,' ['+string(m,format='(i3)')+'] counts smooth: '+string(count_smooth,format='(i5)')+$
        ' sclumpy: '+string(count_sclumpy,format='(i5)')+' bclumpy: '+$
        string(count_bclumpy,format='(i5)')+'  ('+string(cur_frac_found*100,format='(f4.1)')+'%)'
      
      ; free some memory for next load
      Parent = !NULL
      all    = !NULL
      w_gal  = !NULL
      w_gmem = !NULL
    endfor ;m
    
  endif
  
  ; VELOCITY TRACERS case - will be similar to above since there could be multiple
  ; ---------------------
  if sP.trMCPerCell eq -1 then begin
    message,'todo'
  endif  
  
end

; accModeInds(): subselect in the mtS/atS/traj subsets for a particular accretion mode

function accModeInds, at=at, sP=sP, accMode=accMode, mask=mask

  if n_elements(at) eq 0 or n_elements(accMode) eq 0 then message,'Error: Inputs'
  
  gal_w  = where(at.AccTime_gal[0,*] ne -1,count_gal)
  gmem_w = where(at.AccTime_gmem[0,*] ne -1,count_gmem)
  
  ; select on accretion mode by modifying gal_w and gmem_w
  if accMode ne 'all' then begin
    if sP.trMCPerCell ne 0 then message,'Error: accMode for tracers not yet.'
    am = accretionMode(sP=sP)
    
    if accMode eq 'smooth' then begin
      gal_w  = gal_w[where(am.accMode_gal eq 1,count_gal)]
      gmem_w = gmem_w[where(am.accMode_gmem eq 1,count_gmem)]
    endif
    if accMode eq 'bclumpy' then begin
      gal_w  = gal_w[where(am.accMode_gal eq 3,count_gal)]
      gmem_w = gmem_w[where(am.accMode_gmem eq 3,count_gmem)]
    endif
    if accMode eq 'sclumpy' then begin
      gal_w  = gal_w[where(am.accMode_gal eq 2,count_gal)]
      gmem_w = gmem_w[where(am.accMode_gmem eq 2,count_gmem)]
    endif
    if accMode eq 'clumpy' then begin
      gal_w  = gal_w[where(am.accMode_gal eq 2 or am.accMode_gal eq 3,count_gal)]
      gmem_w = gmem_w[where(am.accMode_gmem eq 2 or am.accMode_gmem eq 3,count_gmem)]
    endif
    
    am = !NULL
  endif
  
  ; make mask to combine selection with additional selections, if requested
  if keyword_set(mask) then begin
    galMask  = bytarr(n_elements(at.accTime_gal[0,*]))
    gmemMask = bytarr(n_elements(at.accTime_gmem[0,*]))
    galMask[gal_w] = 1B
    gmemMask[gmem_w] = 1B
    r = {gal:gal_w,gmem:gmem_w,galMask:galmask,gmemMask:gmemMask}
  endif else begin
    r = {gal:gal_w,gmem:gmem_w}
  endelse
  
  return,r
end
