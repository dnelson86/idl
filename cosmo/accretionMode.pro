; accretionMode.pro
; gas accretion project - past substructure history of gas elements
; dnelson jun.2014

; -----------------------------------------------------------------------------------------------------
; accretionMode(): for eaching gas particle/tracer with a recorded accretion time, starting at some 
;                  redshift, track backwards in time with respect to the tracked parent halos (using 
;                  mergerTree) and classify the mode of accretion by considering its group membership
;                  at the outer bracketing snapshot just prior to accretion (r>~r_vir) as:
;                  
;  1. smooth   = in the (primary) parent of the original halo or not in any subgroup, and likewise back
;                to the beginning of the group catalogs
;  2. sclumpy  = in any subgroup other than the (primary) parent halo that is a secondary subgroup ("small")
;  3. bclumpy  = in any subgroup other than the (primary) parent halo that is a primary subgroup ("big")
;  4. stripped = in the (primary) parent of the original halo or not in any subgroup, but in some
;                non-parent subgroup at some earlier time back to
;  >10. recycled = for GFM_WINDS runs, in a wind at any previous time (in which case 10-accMode gives the 
;                  additional information of the accMode at the time of the -first- rvir crossing, e.g.
;                  accMode=11 indicates recycled but smooth at its first entrance to the halo
; -----------------------------------------------------------------------------------------------------

function accretionMode, sP=sP

  forward_function cosmoTracerChildren, cosmoTracerVelParents, accretionTimes
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; first, walk back through the merger tree and find primary subhalos with good parent histories
  mt = mergerTreeSubset(sP=sP,/verbose)
  at = accretionTimes(sP=sP)

  ; set saveFilename and check for existence  
  saveFilename = sP.derivPath + 'accMode.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(mt.maxSnap)+'-'+str(mt.minSnap)+'.sav'  
  
  if file_test(saveFilename) then begin
    restore, saveFilename
    return, accMode
  endif
  
  ; select those particles/tracers with recorded accretion times
  w_at = where(at.AccTime[sP.atIndMode,*] ne -1,nTotSearch,comp=w_at_comp)
  
  ; load galaxy/group member catalogs at zMin for gas ids to search for
  origSnap = sP.snap
  galcat = galaxyCat(sP=sP)
  
  count_bracketed = 0UL
  
  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    print,'Calculating new accretion mode using ( SPH Particles ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
      
    ; only find accretion mode for those gas particles with accretion times
    galcatSub = galcat.ids[w_at]

    ; convert the accretion time (scale factor) into the outer bracketing snapshot number for each member
    bracketSnap = intarr(n_elements(galcatSub))
                    
    snapTimes = snapNumToRedshift(sP=sP,/all,/time)
    snapTimes = [ snapTimes[where(snapTimes ne -1)], 1.00001 ] ; fail if any are past 1.00001

    bracketSnap = value_locate(snapTimes, reform(at.accTime[sP.atIndMode,w_at]) )
    
    w = where(bracketSnap eq -1 or bracketSnap eq n_elements(snapTimes)-1,count)
    if count gt 0 then message,'Error: Bad gal bracketing.'
    
    ; list of parent subgroup IDs at starting snapshot for each gas particle
    origParIDs = mt.gcIndOrig[w_at]
    
    ; replicate hMinSnap for each child gas element
    gasMinSnap = mt.hMinSnap[mt.gcIndOrig[w_at]]
    
    ; store the main arrays as a structure so we can write them directly
    accMode = intarr(n_elements(galcat.ids))-1
         
    ; debugging
    accMask = intarr(n_elements(galcat.ids))
    
    for m=mt.maxSnap,mt.minSnap,-1 do begin
    
      ; load mergerTree and change to parent IDs
      if m ne mt.maxSnap then begin
        Parent = mergerTree(sP=sP)

        w = where(origParIDs ne -1,count)
        if count eq 0 then message,'error'
        origParIDs[w] = Parent[origParIDs[w]] ; change to parent IDs
          
        ; sanity check no IDs are -1 (we should only be processing the mergerTreeSubset)
        w = where(gasMinSnap lt m-1,count)
        if count gt 0 then if min(origParIDs[w]) lt 0 then message,'Error: Bad parent.'
      endif
      
      sP.snap = m
      ; load local group catalog
      h = loadSnapshotHeader(sP=sP)
      gc = loadGroupCat(sP=sP,/readIDs)
      
      ; keep track of how many of each accretion mode we find at this snapshot
      counts  = { smooth:0UL, sclumpy:0UL, bclumpy:0UL, stripped:0UL, recycled:-1}
      
      ; select those gas particles whose accretion modes are determined at this snapshot
      w_now = where(bracketSnap eq sP.snap,count_now)
      w_now_at = w_at[w_now] ; stamp indices into accMode/accMask

      ; debug: make sure we aren't searching for a tracer that is now untracked
      if count_now gt 0 and max(gasMinSnap[w_now]) gt sP.snap+1 then $
        message,'Error: Searching for untracked'
      
      if count_now eq 0 then continue
      
      ; global match against subgroup ID list (this includes fuzz and unbound gas)
      calcMatch,galcatSub[w_now],gc.IDs,galcat_ind,gc_ind,count=countNow
        
      ; those that don't match, assign to category 1. smooth
      if countNow lt count_now then begin
        ; if none matched, assign all the bracketed particles to smooth
        if countNow eq 0 then begin
          accMode[w_now_at] = 1
          accMask[w_now_at] += 1
          counts.smooth += count_now
        endif else begin
          ; if some matched, assign the complement of those to smooth
          all = bytarr(count_now)
          all[galcat_ind] = 1B
          wNotInSG = where(all eq 0B, ncomp)
          
          accMode[w_now_at[wNotInSG]] = 1
          accMask[w_now_at[wNotInSG]] += 1 ;debug
          counts.smooth += ncomp
        endelse
      endif
      
      ; those that match, calculate subgroup ID they belong to
      if countNow gt 0 then begin
        ; calculate parent subgroup ID each gas particle belongs to
        parIDs = value_locate(gc.subgroupOffset,gc_ind)
        
        ; exclude those that aren't actually within the offsetLenType range
        diff = gc_ind - gc.subgroupOffset[parIDs]
        wInSG = where(diff lt gc.subgroupLen[parIDs],countInSG,comp=wc,ncomp=ncomp)
        
        ; those failing the subgroup bound are smooth
        if ncomp gt 0 then begin
          rInds = w_now_at[galcat_ind[wc]]
          accMode[rInds] = 1
          accMask[rInds] += 1 ;debug
          counts.smooth += ncomp
        endif
          
        if countInSG gt 0 then begin
          ; if their subgroup ID is the traced parent, assign to category 1. smooth
          rInds = w_now_at[galcat_ind[wInSG]]
          w0 = where(parIDs[wInSG] eq origParIDs[galcat_ind[wInSG]],count,comp=wc,ncomp=ncomp)
          if count gt 0 then begin
            accMode[rInds[w0]] = 1
            counts.smooth += count
            accMask[rInds[w0]] += 1 ;debug
          endif
          
          ; only separate into sclumpy/bclumpy for the non-smooth remainder
          if ncomp gt 0 then begin
            wInSG = wInSG[wc]
            rInds = rInds[wc]
          
            ; make ID lists of primary,secondary subgroups
            priSGids = gcIDList(gc=gc,select='pri')
            secSGids = gcIDList(gc=gc,select='sec')
            
            ; if their subgroup ID is not the parent but is a primary sg, assign to category 3. bclumpy
            qInds = value_locate(priSGids,parIDs[wInSG])
            w1 = where(priSGids[qInds] eq parIDs[wInSG],count,comp=wc,ncomp=ncomp)
            if count gt 0 then begin
              accMode[rInds[w1]] = 3
              counts.bclumpy += count
              accMask[rInds[w1]] += 1 ; debug
            endif
            
            ; only search for secondary parents for the non-primary remainder (this is mainly a doublecheck
            ; since we have already excluded every other possibility)
            wInSG = wInSG[wc]
            rInds = rInds[wc]
          
            qInds = value_locate(secSGids,parIDs[wInSG])
            w2 = where(secSGids[qInds] eq parIDs[wInSG],count)
            if count gt 0 then begin
              accMode[rInds[w2]] = 2
              counts.sclumpy += count
              accMask[rInds[w2]] += 1 ; debug
            endif
          endif ; ncomp
        endif ; countInSG
      endif ; countNow

      ; check that all bracketed particles found accretion modes
      if min(accMode[w_now_at]) eq 0 then message,'Error: Not all done.'
      if (counts.smooth+counts.sclumpy+counts.bclumpy) ne count_now then $
        message,'Error: Totals fail to add up.'
      
      ; any gas previously marked smooth, categorize as stripped if it is now in a non-parent halo
      ; note: require we are still tracking its progenitor branch to avoid mischaracterizing
      ; smooth gas as stripped if it is just within its true, untracked parent at some earlier time
      wSmooth = where(accMode[w_at] eq 1 and gasMinSnap le sP.snap,countSmooth)
      
      if countSmooth eq 0 then continue
      
      ; global match against subgroup ID list (this includes fuzz and unbound gas)
      calcMatch,galcatSub[wSmooth],gc.IDs,galcat_ind,gc_ind,count=countNow
      if countNow eq 0 then continue
      
      ; calculate parent subgroup ID each gas particle belongs to
      parIDs = value_locate(gc.subgroupOffset,gc_ind)
          
      ; exclude those that aren't actually within the offsetLenType range
      diff = gc_ind - gc.subgroupOffset[parIDs]
      wInSG = where(diff lt gc.subgroupLen[parIDs],countInSG)
      
      if countInSG gt 0 then begin
        ; if their subgroup ID is -not- the traced parent, assign to category 4. stripped
        
        rInds = w_at[ wSmooth[galcat_ind[wInSG]] ]
        w0 = where(parIDs[wInSG] ne origParIDs[galcat_ind[wInSG]],count)
        if count gt 0 then begin
          accMode[rInds[w0]] = 4
          counts.stripped += count
        endif
      endif ; countInSG
      
      ; count stripped/(smooth+stipped) and clumpy/total fraction
      w = where(accMode eq 1,totSmooth)
      w = where(accMode eq 4,totStripped)
      w = where(accMode eq 2 or accMode eq 3,totClumpy)
      
      totAllModes  = totSmooth + totStripped + totClumpy
      fracSmooth   = float(totSmooth)   / totAllModes * 100
      fracStripped = float(totStripped) / totAllModes * 100
      fracClumpy   = float(totClumpy)   / totAllModes * 100
      
      cur_frac_found = float(count_now) / nTotSearch * 100
      
      print,' ['+string(m,format='(i3)')+'] counts'+$
        ' smooth: ' +string(counts.smooth,format='(i5)')+$
        ' sclumpy: '+string(counts.sclumpy,format='(i5)')+$
        ' bclumpy: '+string(counts.bclumpy,format='(i5)')+$
        ' stripped: '+string(counts.stripped,format='(i5)')
      print, '  curFracFound: ['+string(cur_frac_found,format='(f4.1)')+'%]'+$
        ' (smooth: '+string(fracSmooth,format='(f4.1)')+'%)'+$
        ' (stripped: '+string(fracStripped,format='(f4.1)')+'%)'+$
        ' (clumpy: '+string(fracClumpy,format='(f4.1)')+'%)'
      
      ; free some memory for next load
      Parent  = !NULL
      all     = !NULL
      w_now   = !NULL

    endfor ;m
    
    ; verify we found an accretion mode for every galaxyCat member with an accretion time
    if min(accMask[w_at]) lt 1 then message,'Error: Not all found.'
    if max(accMask[w_at]) gt 1 then message,'Error: Overassignment.'
    w = where(accMask[w_at_comp] ne 0,count)
    if count ne 0 then message,'Error: Assigned to negative one accretion time member.'
    
    ; save
    save,accMode,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
  endif  
  
  ; MONTE CARLO TRACERS CASE - for all gas cells, track back all child tracers
  ; ------------------------
  if sP.trMCPerCell gt 0 then begin
    print,'Calculating new accretion mode using ( TracerMC ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
      
    ; make gcIndOrigTr and ids of sorted tracer children
    origParIDs = mt.gcIndOrigTrMC[w_at]
                   
    if min(origParIDs) lt 0 then message,'Error: Bad starting parents.'
             
    ; replicate hMinSnap for each child tracer
    trMinSnap = mt.hMinSnap[mt.gcIndOrigTrMC[w_at]]
                  
    ; convert the accretion time (scale factor) into the outer bracketing snapshot number for each member
    bracketSnap = intarr(n_elements(galcat.trMC_ids))
                    
    snapTimes = snapNumToRedshift(sP=sP,/all,/time)
    snapTimes = [ snapTimes[where(snapTimes ne -1)], 1.00001 ] ; fail if any are past 1.00001

    bracketSnap = value_locate(snapTimes, reform(at.accTime[sP.atIndMode,w_at]) )
    
    ; main save array so we can write them directly
    accMode = intarr(n_elements(galcat.trMC_ids))
    accMask = intarr(n_elements(galcat.trMC_ids))
              
    ; for recycled mode, we can decide immediately based on the tracer windcounter
    if sP.gfmWinds ne 0 then begin
      ; load wind counters and tracer parents
      tr_windcounter = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_windcounter')
      tr_parids      = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentids')
      
      ; gal: find recycled and override bracketSnap using earliest rvir crossing
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
      
      idIndexMap = getIDIndexMap(tr_ids,minid=minid)
      tr_ids = !NULL
      
      galcat_trinds = idIndexMap[galcat.trMC_ids - minid]

      w = where(tr_windcounter[galcat_trinds[w_at]] gt 0,count)
      if count gt 0 then accMode[w_at[w]] = 10
      
      tr_parids = !NULL
      galcat_trinds = !NULL
      tr_windcounter = !NULL
      
      print,'Found recycled: ['+str(count)+'] ('+string(float(count)/nTotSearch*100,format='(f4.1)')+'%)'
    endif
    
    ; verify accretion time -> snapshot mapping
    w = where(bracketSnap eq -1 or bracketSnap eq n_elements(snapTimes)-1,count)
    if count then message,'Error: Bad gal bracketing.'
    
    for m=mt.maxSnap,mt.minSnap,-1 do begin
      ; load mergerTree and move to Parent
      if m ne mt.maxSnap then begin
        Parent = mergerTree(sP=sP)
        
        ; change to parent IDs for the next snapshot
        w = where(origParIDs ne -1,count)
        if count eq 0 then message,'error'
        origParIDs[w] = Parent[origParIDs[w]]
          
        ; sanity check no IDs are -1 (we should only be processing the mergerTreeSubset)
        w = where(trMinSnap lt m-1,count)
        if count gt 0 then $
          if min(origParIDs[w]) lt 0 then message,'Error: Bad parent.'
      endif
      
      sP.snap = m
      ; load local group catalog
      h = loadSnapshotHeader(sP=sP)
      gc = loadGroupCat(sP=sP,/readIDs)
      
      ; differing from the sph case, we need to map gc.IDs (gas -only-) to a child tracerMC list for each snapshot
      gasIDs = gcPIDList(gc=gc,select='all',partType='gas')
      groupcat_trIDs = cosmoTracerChildren(sP=sP, /getIDs, gasIDs=gasIDs, child_counts=groupcat_cc)

      ; and construct the subgroupLen and subgroupOffset arrays for the tracerMCs in groupcat_trIDs
      subgroupOffsetTr = ulonarr(gc.nSubgroupsTot)
      subgroupLenTr    = ulonarr(gc.nSubgroupsTot)
      
      globalTrCount = 0L
      globalGasCount = 0L
      parPT = partTypeNum('gas')
      
      mask = intarr(n_elements(gasIDs)) ; debugging
      
      for SubNr=0L, gc.nSubgroupsTot-1 do begin
        ; set offset to current marker position
        subgroupOffsetTr[SubNr] = globalTrCount
        
        if gc.subgroupLenType[parPT,SubNr] gt 0 then begin
          ; indices into gc.IDs for this subgroup
          id_inds = lindgen(gc.subgroupLenType[parPT,SubNr]) + globalGasCount ;gc.subgroupOffsetType[parPT,SubNr]
          mask[id_inds] += 1
          ; total number of tracers in this subgroup
          subgroupLenTr[SubNr] = total(groupcat_cc[id_inds],/int)

          globalTrCount  += subgroupLenTr[SubNr]
          globalGasCount += gc.subgroupLenType[parPT,SubNr]
        endif
      endfor
      
      w = where(mask ne 1,countFM)
      if countFM gt 0 then message,'Error: Failed creating new offset array.'
      mask = !NULL
      
      ; keep track of how many of each accretion mode we find at this snapshot
      counts = { smooth : 0UL, sclumpy: 0UL, bclumpy: 0UL, stripped: 0UL, recycled: 0UL }
      
      ; select those gas particles whose accretion modes are determined at this snapshot
      w_now    = where(bracketSnap eq sP.snap,count_now)
      w_now_at = w_at[w_now] ; stamp indices into accMode/accMask
      
      ; debug: make sure we aren't searching for a tracer that is now untracked
      if count_now gt 0 and max(trMinSnap[w_now]) gt sP.snap+1 then message,'Error: Searching for untracked'
      
      if count_now eq 0 then continue
        
      ; global match against subgroup tracers ID list (this includes -only- bound gas)
      calcMatch,galcat.trMC_ids[w_now_at],groupcat_trIDs,galcat_ind,gc_ind,count=countNow
        
      ; those that don't match, assign to category 1. smooth
      if countNow lt count_now then begin
        ; if none matched, assign all the bracketed particles to smooth
        if countNow eq 0 then begin
          accMode[w_now_at] += 1
          accMask[w_now_at] += 1
          counts.smooth += count_now
        endif else begin
          ; if some matched, assign the complement of those to smooth
          all = bytarr(count_now)
          all[galcat_ind] = 1B
          wNotInSG = where(all eq 0B, ncomp)
            
          accMode[w_now_at[wNotInSG]] += 1
          accMask[w_now_at[wNotInSG]] += 1 ;debug
          counts.smooth += ncomp
        endelse
      endif
        
      ; those that match, calculate subgroup ID they belong to
      if countNow gt 0 then begin
        ; calculate parent subgroup ID each gas particle belongs to
        parIDs = value_locate(subgroupOffsetTr,gc_ind)
          
        ; exclude those that aren't actually within the offsetLenType range
        diff = gc_ind - subgroupOffsetTr[parIDs]
        wInSG = where(diff lt subgroupLenTr[parIDs],countInSG,comp=wc,ncomp=ncomp)
          
        ; those failing the subgroup bound are smooth
        if ncomp gt 0 then message,'This should not happen for tracerMC.'
          
        if countInSG gt 0 then begin
          ; if their subgroup ID is the traced parent, assign to category 1. smooth
          rInds = w_now_at[galcat_ind[wInSG]]
          w0 = where(parIDs[wInSG] eq origParIDs[galcat_ind[wInSG]],count,comp=wc,ncomp=ncomp)
          if count gt 0 then begin
            accMode[rInds[w0]] += 1
            counts.smooth += count
            accMask[rInds[w0]] += 1 ;debug
          endif
            
          ; only separate into sclumpy/bclumpy for the non-smooth remainder
          if ncomp gt 0 then begin
            wInSG = wInSG[wc]
            rInds = rInds[wc]
            
            ; make ID lists of primary,secondary subgroups
            priSGids = gcIDList(gc=gc,select='pri')
            secSGids = gcIDList(gc=gc,select='sec')
              
            ; if their subgroup ID is not the parent but is a primary sg, assign to category 3. bclumpy
            qInds = value_locate(priSGids,parIDs[wInSG])
            w1 = where(priSGids[qInds] eq parIDs[wInSG],count,comp=wc,ncomp=ncomp)
            if count gt 0 then begin
              accMode[rInds[w1]] += 3
              counts.bclumpy += count
              accMask[rInds[w1]] += 1 ; debug
            endif
              
            ; only search for secondary parents for the non-primary remainder (this is mainly a doublecheck
            ; since we have already excluded every other possibility)
            wInSG = wInSG[wc]
            rInds = rInds[wc]
            
            qInds = value_locate(secSGids,parIDs[wInSG])
            w2 = where(secSGids[qInds] eq parIDs[wInSG],count)
            if count gt 0 then begin
              accMode[rInds[w2]] += 2
              counts.sclumpy += count
              accMask[rInds[w2]] += 1 ; debug
            endif
          endif ; ncomp
        endif ; countInSG
      endif ; countNow
        
      ; check that all bracketed particles found accretion modes
      if min(accMode[w_now_at]) eq 0 then message,'Error: Not all done.'
      if (counts.smooth+counts.sclumpy+counts.bclumpy) ne count_now then message,'Error: Totals fail to add up.'
        
      ; any gas previously marked smooth, categorize as stripped if it is now in a non-parent halo
      smoothAccModes = [1,11]
        
      foreach smoothAccMode,smoothAccModes do begin       
        ; find either smooth only or smooth+recycled
        wSmooth = where(accMode[w_at] eq smoothAccMode and trMinSnap le sP.snap,countSmooth)
        
        if countSmooth eq 0 then continue
          
        ; global match against subgroup ID list (this includes fuzz and unbound gas)
        calcMatch,galcat.trMC_ids[w_at[wSmooth]],groupcat_trIDs,galcat_ind,gc_ind,count=countNow
        
        if countNow eq 0 then continue
          
        ; calculate parent subgroup ID each gas particle belongs to
        parIDs = value_locate(subgroupOffsetTr,gc_ind)
          
        ; exclude those that aren't actually within the offsetLenType range
        diff = gc_ind - subgroupOffsetTr[parIDs]
        wInSG = where(diff lt subgroupLenTr[parIDs],countInSG)

        if countInSG ne countNow then message,'Error: Should not happen (gal).'
          
        ; if their subgroup ID is -not- the traced parent, assign to category 4. stripped
        rInds = w_at[ wSmooth[galcat_ind[wInSG]] ]
        w0 = where(parIDs[wInSG] ne origParIDs[galcat_ind[wInSG]],count)
        if count gt 0 then begin
          accMode[rInds[w0]] = 4 + (smoothAccMode-1) ; addition is zero for accMode=1, 10 for recycled
          counts.stripped += count
        endif
          
      endforeach ; smoothAccModes
      
      ; count stripped/(smooth+stipped) and clumpy/total fraction
      w = where(accMode eq 1,totSmooth)
      w = where(accMode eq 4,totStripped)
      w = where(accMode eq 2 or accMode eq 3,totClumpy)
      
      fracStripped = float(totStripped) / (totSmooth+totStripped)
      fracClumpy   = float(totClumpy) / (totSmooth+totStripped+totClumpy)
      
      cur_frac_found = float(count_now) / nTotSearch
      
      ; count fractions again for recycled
      w = where(accMode eq 11,totSmooth)
      w = where(accMode eq 14,totStripped)
      w = where(accMode eq 12 or accMode eq 13,totClumpy)
      
      fracStripped2 = float(totStripped) / (totSmooth+totStripped)
      fracClumpy2   = float(totClumpy) / (totSmooth+totStripped+totClumpy)
      
      print,' ['+string(m,format='(i3)')+'] counts'+$
        ' smooth: ' +string(counts.smooth,format='(i5)')   + $
        ' sclumpy: '+string(counts.sclumpy,format='(i5)')  + $
        ' bclumpy: '+string(counts.bclumpy,format='(i5)') + $
        ' stripped: '+string(counts.stripped,format='(i5)')
        
      print,' ['+string(m,format='(i3)')+'] fractions tot ('+string(cur_frac_found*100,format='(f4.1)')+'%)'+$
        ' stripped ('+string(fracStripped*100,format='(f4.1)')+'%)'+$
        ' clumpy ('+string(fracClumpy*100,format='(f4.1)')+'%)'+$
        ' strippedRec ('+string(fracStripped2*100,format='(f4.1)')+'%)'+$
        ' clumpyRec ('+string(fracClumpy2*100,format='(f4.1)')+'%)'
      
      ; free some memory for next load
      Parent  = !NULL
      all     = !NULL
      w_now   = !NULL
    endfor ;m
    
    ; verify we found an accretion mode for every gas particle
    if min(accMask[w_at]) lt 1 then message,'Error: Not all found.'
    if max(accMask[w_at]) gt 1 then message,'Error: Overassignment.'
    w = where(accMask[w_at_comp] ne 0,count)
    if count ne 0 then message,'Error: Assigned to negative one accretion time member.'
    
    ; verify no strange accMode numbers
    w = where((accMode gt 4 and accMode lt 10) or accMode gt 14,count)
    if count gt 0 then message,'Error: Bad accMode.'
    
    ; save
    save,accMode,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
  endif
  
  ; VELOCITY TRACERS case - will be similar to above since there could be multiple
  ; ---------------------
  if sP.trMCPerCell eq -1 then begin
    message,'todo'
  endif
  
  sP.snap = origSnap ; restore sP.snap
  return, accMode
end

; accModeInds(): subselect for a particular accretion mode (indices are for galcat.ids or galcat.trMC_ids)

function accModeInds, at=at, mt=mt, sP=sP, accMode=accMode, maskAndInds=maskAndInds

  if (n_elements(at) eq 0 and n_elements(mt) eq 0) or n_elements(accMode) eq 0 then message,'Error: Inputs'
  if ~sP.gfmWinds and accMode eq 'recycled' then message,'Error: Request recycled on non-winds run.'
  
  if n_elements(at) gt 0 then nElem = n_elements(at.accTime[0,*])
  if n_elements(mt) gt 0 then begin
    if sP.trMCPerCell gt 0 then nElem = n_elements(mt.gcIndOrigTrMC)
    if sP.trMCPerCell eq 0 then nElem = n_elements(mt.gcIndOrig)
    if sP.trMCPerCell lt 0 then nElem = n_elements(mt.gcIndOrigTrVel)
  endif
  
  ; select on accretion mode
  if sP.trMCPerCell lt 0 then message,'Error: accMode for tracerVel not yet.'
  am = accretionMode(sP=sP)
  if n_elements(am) ne nElem then message,'Error: Inconsistency.'
    
  if accMode eq 'smooth'       then at_w = where(am eq 1,count)
  if accMode eq 'smooth_rec'   then at_w = where(am eq 1 or am eq 11,count)
  if accMode eq 'bclumpy'      then at_w = where(am eq 3,count)
  if accMode eq 'sclumpy'      then at_w = where(am eq 2,count)
  if accMode eq 'clumpy'       then at_w = where(am eq 2 or am eq 3,count)
  if accMode eq 'clumpy_rec'   then at_w = where(am eq 2 or am eq 3 or am eq 12 or am eq 13,count)
  if accMode eq 'stripped'     then at_w = where(am eq 4,count)
  if accMode eq 'stripped_rec' then at_w = where(am eq 4 or am eq 14,count)
  if accMode eq 'recycled'     then at_w = where(am ge 10,count) ; recycled+any other mode
  if accMode eq 'all'          then at_w = where(am gt 0,count) ; full atS
  
  if n_elements(count) eq 0 then message,'Error: Unrecognized accretion mode.'
  if count eq 0 then message,'Error: Bad accModeInds selection.'
  
  am = !NULL

  ; make mask
  mask = bytarr(nElem)
  mask[at_w] = 1B
  
  ; load galaxyCat to replicate type for each tracer
  galcat = galaxyCat(sP=sP)
  
  if sP.trMCPerCell gt 0 then type = ( galcat.type[ replicate_var(galcat.trMC_cc) ] )
  if sP.trMCPerCell eq 0 then type = ( galcat.type )
  if sP.trMCPerCell lt 0 then type = ( galcat.type[ replicate_var(galcat.trVel_cc) ] )
  if n_elements(type) ne n_elements(mask) then message,'Error: Inconsistency.'
  
  ; split at_w into types
  rr = {}  
  totCount = 0L

  for i=0,n_tags(galcat.types)-1 do begin
    ; select from atS (note: these indices give the gal/gmem/etc subsets from the GLOBAL mtS)
    w_type = where(type eq galcat.types.(i) and mask eq 1B,typeCount)
    
    if typeCount gt 0 then $
      rr = mod_struct(rr, (tag_names(galcat.types))[i], w_type) $
    else $
      rr = mod_struct(rr, (tag_names(galcat.types))[i], -1)
    totCount += typeCount
  endfor
  
  if totCount ne n_elements(at_w) then message,'Error in splitting.'
  
  ; include mask and unsplit indices at_w in return
  if keyword_set(maskAndInds) then begin
    rr = mod_struct( rr, 'mask', mask )
    rr = mod_struct( rr, 'at_w', at_w )
  endif
  
  return,rr
end
