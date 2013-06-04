; accretionMode.pro
; gas accretion project - past substructure history of gas elements
; dnelson jun.2013

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
    return, r
  endif
  
  ; select those particles/tracers with recorded accretion times
  gal_w_at   = where(at.AccTime_gal[0,*] ne -1,count_gal)
  gmem_w_at  = where(at.AccTime_gmem[0,*] ne -1,count_gmem)
  stars_w_at = where(at.AccTime_stars[0,*] ne -1,count_stars)
  
  nTotSearch = n_elements(gal_w_at) + n_elements(gmem_w_at) + n_elements(stars_w_at)
  
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
    galcatSub = { gal   : galcat.galaxyIDs[mt.galcatSub.gal[gal_w_at]]    ,$
                  gmem  : galcat.groupmemIDs[mt.galcatSub.gmem[gmem_w_at]]  ,$
                  stars : galcat.stellarIDs[mt.galcatSub.stars[stars_w_at]] }

    ; convert the accretion time (scale factor) into the outer bracketing snapshot number for each member
    bracketSnap = { gal   : intarr(n_elements(galcatSub.gal))   ,$
                    gmem  : intarr(n_elements(galcatSub.gmem))  ,$
                    stars : intarr(n_elements(galcatSub.stars))  }
                    
    snapTimes = snapNumToRedshift(sP=sP,/all,/time)
    snapTimes = snapTimes[where(snapTimes ne -1)]

    bracketSnap.gal   = value_locate(snapTimes,at.accTime_gal[0,gal_w_at])
    bracketSnap.gmem  = value_locate(snapTimes,at.accTime_gmem[0,gmem_w_at])
    bracketSnap.stars = value_locate(snapTimes,at.accTime_stars[0,stars_w_at])
    
    w = where(bracketSnap.gal eq -1 or bracketSnap.gal eq n_elements(snapTimes)-1,count)
    if count gt 0 then message,'Error: Bad gal bracketing.'
    w = where(bracketSnap.gmem eq -1 or bracketSnap.gmem eq n_elements(snapTimes)-1,count)
    if count gt 0 then message,'Error: Bad gmem bracketing.'
    w = where(bracketSnap.stars eq -1 or bracketSnap.stars eq n_elements(snapTimes)-1,count)
    if count gt 0 then message,'Error: Bad stars bracketing.'
    
    ; list of parent subgroup IDs at starting snapshot for each gas particle
    gcIndOrig = galCatRepParentIDs(galcat=galcat,gcIDList=mt.galcatIDList)
    
    origParIDs = { gal   : gcIndOrig.gal[gal_w_at]    ,$
                   gmem  : gcIndOrig.gmem[gmem_w_at]  ,$
                   stars : gcIndOrig.stars[stars_w_at] }
                   
    gcIndOrig = !NULL ; includes primaries+secondaries for Parent (not compacted)
    
    ; replicate hMinSnap for each child gas element
    gasMinSnap = { gal   : mt.hMinSnap[mt.gcIndOrig.gal[gal_w_at]]     ,$
                   gmem  : mt.hMinSnap[mt.gcIndOrig.gmem[gmem_w_at]]   ,$
                   stars : mt.hMinSnap[mt.gcIndOrig.stars[stars_w_at]]  }     
    
    gal_w_at   = !NULL
    gmem_w_at  = !NULL
    stars_w_at = !NULL
    
    ; store the main arrays as a structure so we can write them directly
    r = {accMode_gal       : intarr(n_elements(galcatSub.gal))-1  ,$
         accMode_gmem      : intarr(n_elements(galcatSub.gmem))-1 ,$
         accMode_stars     : intarr(n_elements(galcatSub.stars))-1 }
         
    ; debugging
    rMask = { gal   : intarr(n_elements(galcatSub.gal))    ,$
              gmem  : intarr(n_elements(galcatSub.gmem))   ,$
              stars : intarr(n_elements(galcatSub.stars))   }
    
    for m=mt.maxSnap,mt.minSnap,-1 do begin
      sP.snap = m
      ; load local group catalog
      h = loadSnapshotHeader(sP=sP)
      gc = loadGroupCat(sP=sP,/readIDs)
      
      ; we only want to match galcat.stellarIDs against gas IDs in groupcat, not against any star ids
      ; just set all star blocks in gc.IDs to negative to avoid any such matches later (keep unique)
      star_ids = loadSnapshotSubset(sP=sP,partType='stars',field='ids')
      calcMatch,star_ids,gc.IDs,star_ids_ind,gc_ids_ind,count=countMatch
      if countMatch gt 0 then gc.IDs[gc_ids_ind] *= -1
      
      star_ids     = !NULL
      star_ids_ind = !NULL
      gc_ids_ind   = !NULL
      
      ; keep track of how many of each accretion mode we find at this snapshot
      count_smooth   = { gal: 0UL, gmem: 0UL, stars: 0UL }
      count_sclumpy  = { gal: 0UL, gmem: 0UL, stars: 0UL }
      count_bclumpy  = { gal: 0UL, gmem: 0UL, stars: 0UL }
      count_stripped = { gal: 0UL, gmem: 0UL, stars: 0UL }
      
      ; loop over tracer parent types (gal,gmem,stars)
      ; ----------------------------------------------
      for k=0,n_tags(bracketSnap)-1 do begin
        ; select those gas particles whose accretion modes are determined at this snapshot
        w_now   = where(bracketSnap.(k) eq sP.snap,count_now)
      
        ; debug: make sure we aren't searching for a tracer that is now untracked
        if count_now gt 0 and max(gasMinSnap.(k)[w_now]) gt sP.snap+1 then message,'Error: Searching for untracked'
      
        if count_now eq 0 then continue
      
        ; global match against subgroup ID list (this includes fuzz and unbound gas)
        calcMatch,galcatSub.(k)[w_now],gc.IDs,galcat_ind,gc_ind,count=countNow
        
        ; those that don't match, assign to category 1. smooth
        if countNow lt count_now then begin
          ; if none matched, assign all the bracketed particles to smooth
          if countNow eq 0 then begin
            r.(k)[w_now] = 1
            rMask.(k)[w_now] += 1
            count_smooth.(k) += count_now
          endif else begin
            ; if some matched, assign the complement of those to smooth
            all = bytarr(count_now)
            all[galcat_ind] = 1B
            wNotInSG = where(all eq 0B, ncomp)
            
            r.(k)[w_now[wNotInSG]] = 1
            rMask.(k)[w_now[wNotInSG]] += 1 ;debug
            count_smooth.(k) += ncomp
          endelse
        endif
        
        ; those that match, calculate subgroup ID they belong to
        if countNow gt 0 then begin
          ; calculate parent subgroup ID each gas particle belongs to
          parIDs = value_locate(gc.subgroupOffsetType[partTypeNum('gas'),*],gc_ind)
          
          ; exclude those that aren't actually within the offsetLenType range
          diff = gc_ind - gc.subgroupOffsetType[partTypeNum('gas'),parIDs]
          wInSG = where(diff lt gc.subgroupLenType[partTypeNum('gas'),parIDs],countInSG,comp=wc,ncomp=ncomp)
          
          ; those failing the subgroup bound are smooth
          if ncomp gt 0 then begin
            rInds = w_now[galcat_ind[wc]]
            r.(k)[rInds] = 1
            rMask.(k)[rInds] += 1 ;debug
            count_smooth.(k) += ncomp
          endif
          
          if countInSG gt 0 then begin
            ; if their subgroup ID is the traced parent, assign to category 1. smooth
            rInds = w_now[galcat_ind[wInSG]]
            w0 = where(parIDs[wInSG] eq origParIDs.(k)[galcat_ind[wInSG]],count,comp=wc,ncomp=ncomp)
            if count gt 0 then begin
              r.(k)[rInds[w0]] = 1
              count_smooth.(k) += count
              rMask.(k)[rInds[w0]] += 1 ;debug
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
                r.(k)[rInds[w1]] = 3
                count_bclumpy.(k) += count
                rMask.(k)[rInds[w1]] += 1 ; debug
              endif
              
              ; only search for secondary parents for the non-primary remainder (this is mainly a doublecheck
              ; since we have already excluded every other possibility)
              wInSG = wInSG[wc]
              rInds = rInds[wc]
            
              qInds = value_locate(secSGids,parIDs[wInSG])
              w2 = where(secSGids[qInds] eq parIDs[wInSG],count)
              if count gt 0 then begin
                r.(k)[rInds[w2]] = 2
                count_sclumpy.(k) += count
                rMask.(k)[rInds[w2]] += 1 ; debug
              endif
            endif ; ncomp
          endif ; countInSG
        endif ; countNow
        
        ; check that all bracketed particles found accretion modes
        if min(r.(k)[w_now]) eq 0 then message,'Error: Not all done.'
        if (count_smooth.(k)+count_sclumpy.(k)+count_bclumpy.(k)) ne count_now then $
          message,'Error: Totals fail to add up.'
        if max(rMask.(k)) gt 1 then message,'Error: Mask exceeds one.'
        
        ; any gas previously marked smooth, categorize as stripped if it is now in a non-parent halo
        ; note: require we are still tracking its progenitor branch to avoid mischaracterizing
        ; smooth gas as stripped if it is just within its true, untracked parent at some earlier time
        wSmooth = where(r.(k) eq 1 and gasMinSnap.(k) le sP.snap,countSmooth)
        
        if countSmooth eq 0 then continue
        
        ; global match against subgroup ID list (this includes fuzz and unbound gas)
        calcMatch,galcatSub.(k)[wSmooth],gc.IDs,galcat_ind,gc_ind,count=countNow

        if countNow eq 0 then continue
        
        ; calculate parent subgroup ID each gas particle belongs to
        parIDs = value_locate(gc.subgroupOffsetType[partTypeNum('gas'),*],gc_ind)
          
        ; exclude those that aren't actually within the offsetLenType range
        diff = gc_ind - gc.subgroupOffsetType[partTypeNum('gas'),parIDs]
        wInSG = where(diff lt gc.subgroupLenType[partTypeNum('gas'),parIDs],countInSG)

        if countInSG gt 0 then begin
          ; if their subgroup ID is -not- the traced parent, assign to category 4. stripped
          rInds = wSmooth[galcat_ind[wInSG]]
          w0 = where(parIDs[wInSG] ne origParIDs.(k)[galcat_ind[wInSG]],count)
          if count gt 0 then begin
            r.(k)[rInds[w0]] = 4
            count_stripped.(k) += count
          endif
        endif ; countInSG
        
      endfor ; k
          
      if m ne mt.minSnap then begin
        ; load mergerTree and change to parent IDs
        Parent = mergerTree(sP=sP)
        
        for k=0,n_tags(bracketSnap)-1 do begin
        
          w = where(origParIDs.(k) ne -1,count)
          if count eq 0 then message,'error'
          origParIDs.(k)[w] = Parent[origParIDs.(k)[w]] ; change to parent IDs
        
          ; sanity check no IDs are -1 (we should only be processing the mergerTreeAdaptiveSubset)
          w = where(gasMinSnap.(k) lt m-1,count)
          if count gt 0 then $
            if min(origParIDs.(k)[w]) lt 0 then message,'Error: Bad parent.'
          
        endfor ;k
      endif
      
      ; count stripped/(smooth+stipped) and clumpy/total fraction
      w = where(r.accMode_gal eq 1,totSmoothGal)
      w = where(r.accMode_gal eq 4,totStrippedGal)
      w = where(r.accMode_gal eq 2 or r.accMode_gal eq 3,totClumpyGal)
      
      w = where(r.accMode_gmem eq 1,totSmoothGmem)
      w = where(r.accMode_gmem eq 4,totStrippedGmem)
      w = where(r.accMode_gmem eq 2 or r.accMode_gmem eq 3,totClumpyGmem)
      
      w = where(r.accMode_stars eq 1,totSmoothStars)
      w = where(r.accMode_stars eq 4,totStrippedStars)
      w = where(r.accMode_stars eq 2 or r.accMode_stars eq 3,totClumpyStars)
      
      fracStrippedGal   = float(totStrippedGal)   / (totSmoothGal+totStrippedGal)
      fracStrippedGmem  = float(totStrippedGmem)  / (totSmoothGmem+totStrippedGmem)
      fracStrippedStars = float(totStrippedStars) / (totSmoothStars+totStrippedStars)
      
      fracClumpyGal    = float(totClumpyGal)   / (totSmoothGal+totStrippedGal+totClumpyGal)
      fracClumpyGmem   = float(totClumpyGmem)  / (totSmoothGmem+totStrippedGmem+totClumpyGmem)
      fracClumpyStars  = float(totClumpyStars) / (totSmoothStars+totStrippedStars+totClumpyStars)
      
      count_bracketed += count_gal + count_gmem + count_stars
      cur_frac_found = float(count_bracketed) / nTotSearch
      
      print,' ['+string(m,format='(i3)')+'] counts'+$
        ' smooth: ' +string(count_smooth.gal,format='(i5)')   +'.'+string(count_smooth.gmem,format='(i5)')+'.'+string(count_smooth.stars,format='(i5)')+$
        ' sclumpy: '+string(count_sclumpy.gal,format='(i5)')  +'.'+string(count_sclumpy.gmem,format='(i5)')+'.'+string(count_sclumpy.stars,format='(i5)')+$
        ' bclumpy: '+ string(count_bclumpy.gal,format='(i5)') +'.'+string(count_bclumpy.gmem,format='(i5)')+'.'+string(count_bclumpy.stars,format='(i5)')+$
        ' stripped: '+string(count_stripped.gal,format='(i5)')+'.'+string(count_stripped.gmem,format='(i5)')+'.'+string(count_stripped.stars,format='(i5)')+$
        '  ('+string(cur_frac_found*100,format='(f4.1)')+'%)'+$
        ' ('+string(fracStrippedGal*100,format='(f4.1)')+'% '+string(fracStrippedGmem*100,format='(f4.1)')+'% '+string(fracStrippedStars*100,format='(f4.1)')+'%)'+$
        ' ('+string(fracClumpyGal*100,format='(f4.1)')+'% '+string(fracClumpyGmem*100,format='(f4.1)')+'% '+string(fracClumpyStars*100,format='(f4.1)')+'%)'
      
      ; free some memory for next load
      Parent  = !NULL
      all     = !NULL
      w_gal   = !NULL
      w_gmem  = !NULL
      w_stars = !NULL
    endfor
    
    ; verify we found an accretion mode for every gas particle
    if min(rMask.gal)   lt 1 then message,'Error: Not all gal found.'
    if min(rMask.gmem)  lt 1 then message,'Error: Not all gmem found.'
    if min(rMask.stars) lt 1 then message,'Error: Not all stars found.'
    
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
                    trids_gal=galcat_gal_trids, trids_gmem=galcat_gmem_trids, trids_stars=galcat_stars_trids)
    
    origParIDs = { gal   : gcIndOrigTr.gal[gal_w_at]    ,$
                   gmem  : gcIndOrigTr.gmem[gmem_w_at]  ,$
                   stars : gcIndOrigTr.stars[stars_w_at] }
                   
    if min(origParIDs.gal) lt 0 or min(origParIDs.gmem) lt 0 or $
       min(origParIDs.stars) lt 0 then message,'Error: Bad starting parents.'
             
    ; replicate hMinSnap for each child tracer (do compactMtS on gcIndOrigTr first)
    placeMap = getIDIndexMap(mt.galcatIDList,minid=minid)
    gcIndOrigTr.gal   = placeMap[gcIndOrigTr.gal-minid]
    gcIndOrigTr.gmem  = placeMap[gcIndOrigTr.gmem-minid]
    gcIndOrigTr.stars = placeMap[gcIndOrigTr.stars-minid]
    placeMap = !NULL
    
    trMinSnap = { gal   : mt.hMinSnap[gcIndOrigTr.gal[gal_w_at]]  ,$
                  gmem  : mt.hMinSnap[gcIndOrigTr.gmem[gmem_w_at]] ,$
                  stars : mt.hMinSnap[gcIndOrigTr.stars[stars_w_at]] } 
         
    gcIndOrigTr = !NULL
    
    ; only find accretion mode for those tracers with accretion times
    galcatSub = { gal   : galcat_gal_trids[gal_w_at]    ,$
                  gmem  : galcat_gmem_trids[gmem_w_at]  ,$
                  stars : galcat_stars_trids[stars_w_at] }
                  
    galcat_gal_trids   = !NULL
    galcat_gmem_trids  = !NULL
    galcat_stars_trids = !NULL
                  
    ; convert the accretion time (scale factor) into the outer bracketing snapshot number for each member
    bracketSnap = { gal   : intarr(n_elements(galcatSub.gal))   ,$
                    gmem  : intarr(n_elements(galcatSub.gmem))  ,$
                    stars : intarr(n_elements(galcatSub.stars))  }
                    
    snapTimes = snapNumToRedshift(sP=sP,/all,/time)
    snapTimes = snapTimes[where(snapTimes ne -1)]

    bracketSnap.gal   = value_locate(snapTimes,at.accTime_gal[0,gal_w_at])
    bracketSnap.gmem  = value_locate(snapTimes,at.accTime_gmem[0,gmem_w_at])
    bracketSnap.stars = value_locate(snapTimes,at.accTime_stars[0,stars_w_at])
    
    ; store the main arrays as a structure so we can write them directly
    r = {accMode_gal       : intarr(n_elements(galcatSub.gal))  ,$
         accMode_gmem      : intarr(n_elements(galcatSub.gmem)) ,$
         accMode_stars     : intarr(n_elements(galcatSub.stars)) }
         
    ; debugging
    rMask = { gal   : intarr(n_elements(galcatSub.gal))   ,$
              gmem  : intarr(n_elements(galcatSub.gmem))  ,$
              stars : intarr(n_elements(galcatSub.stars))  }
              
    ; for recycled mode, we can decide immediately based on the tracer windcounter
    if sP.gfmWinds ne 0 then begin
      ; load wind counters and tracer parents
      tr_windcounter = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_windcounter')
      tr_parids      = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentids')
      
      ; gal: find recycled and override bracketSnap using earliest rvir crossing
      ; changed june 2013: so that 'secondary' modes for recycled tracers are consistent, 
      ; make the determination at the latest rvir crossing as normal
      galcat_trinds = cosmoTracerChildren(sP=sP, /getInds, tr_parids=tr_parids, gasIDs=galcat.galaxyIDs[mt.galcatSub.gal])
      w = where(tr_windcounter[galcat_trinds[gal_w_at]] gt 0,count1)
      if count1 gt 0 then begin
        r.accMode_gal[w] = 10
        ;bracketSnap.gal[w] = value_locate(snapTimes,at.accTime_gal[-1,gal_w_at[w]])
      endif
      
      ; gmem
      galcat_trinds = cosmoTracerChildren(sP=sP, /getInds, tr_parids=tr_parids, gasIDs=galcat.groupmemIDs[mt.galcatSub.gmem])
      w = where(tr_windcounter[galcat_trinds[gmem_w_at]] gt 0,count2)
      if count2 gt 0 then begin
        r.accMode_gmem[w] = 10
        ;bracketSnap.gmem[w] = value_locate(snapTimes,at.accTime_gmem[-1,gmem_w_at[w]])
      endif
      
      ; stars
      galcat_trinds = cosmoTracerChildren(sP=sP, /getInds, tr_parids=tr_parids, starIDs=galcat.stellarIDs[mt.galcatSub.stars])
      w = where(tr_windcounter[galcat_trinds[stars_w_at]] gt 0,count3)
      if count3 gt 0 then begin
        r.accMode_stars[w] = 10
        ;bracketSnap.stars[w] = value_locate(snapTimes,at.accTime_stars[-1,stars_w_at[w]])
      endif
      
      tr_parids = !NULL
      galcat_trinds = !NULL
      tr_windcounter = !NULL
      print,'Found recycled: ['+str(count1)+'] gal, ['+str(count2)+'] gmem, ['+str(count3)+'] stars ('+$
        string(float(count1)/count_gal*100,format='(f4.1)')+'% '+$
        string(float(count2)/count_gmem*100,format='(f4.1)')+'% '+$
        string(float(count3)/count_stars*100,format='(f4.1)')+'%)'
    endif
    
    ; verify accretion time -> snapshot mapping
    w = where(bracketSnap.gal eq -1 or bracketSnap.gal eq n_elements(snapTimes)-1,count)
    if count then message,'Error: Bad gal bracketing.'
    w = where(bracketSnap.gmem eq -1 or bracketSnap.gmem eq n_elements(snapTimes)-1,count)
    if count then message,'Error: Bad gmem bracketing.'
    w = where(bracketSnap.stars eq -1 or bracketSnap.stars eq n_elements(snapTimes)-1,count)
    if count then message,'Error: Bad stars bracketing.'
    
    gal_w_at   = !NULL
    gmem_w_at  = !NULL
    stars_w_at = !NULL
    galcat = !NULL
    
    for m=mt.maxSnap,mt.minSnap,-1 do begin
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
      count_smooth   = { gal: 0UL, gmem: 0UL, stars: 0UL }
      count_sclumpy  = { gal: 0UL, gmem: 0UL, stars: 0UL }
      count_bclumpy  = { gal: 0UL, gmem: 0UL, stars: 0UL }
      count_stripped = { gal: 0UL, gmem: 0UL, stars: 0UL }
      
      ; loop over tracer parent types (gal,gmem,stars)
      ; ----------------------------------------------
      for k=0,n_tags(bracketSnap)-1 do begin
 
        ; select those gas particles whose accretion modes are determined at this snapshot
        w_now   = where(bracketSnap.(k) eq sP.snap,count_now)
      
        ; debug: make sure we aren't searching for a tracer that is now untracked
        if count_now gt 0 and max(trMinSnap.(k)[w_now]) gt sP.snap+1 then message,'Error: Searching for untracked'
      
        if count_now eq 0 then continue
        
        ; global match against subgroup tracers ID list (this includes -only- bound gas)
        calcMatch,galcatSub.(k)[w_now],groupcat_trIDs,galcat_ind,gc_ind,count=countNow
        
        ; those that don't match, assign to category 1. smooth
        if countNow lt count_now then begin
          ; if none matched, assign all the bracketed particles to smooth
          if countNow eq 0 then begin
            r.(k)[w_now] += 1
            rMask.(k)[w_now] += 1
            count_smooth.(k) += count_now
          endif else begin
            ; if some matched, assign the complement of those to smooth
            all = bytarr(count_now)
            all[galcat_ind] = 1B
            wNotInSG = where(all eq 0B, ncomp)
            
            r.(k)[w_now[wNotInSG]] += 1
            rMask.(k)[w_now[wNotInSG]] += 1 ;debug
            count_smooth.(k) += ncomp
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
            rInds = w_now[galcat_ind[wInSG]]
            w0 = where(parIDs[wInSG] eq origParIDs.(k)[galcat_ind[wInSG]],count,comp=wc,ncomp=ncomp)
            if count gt 0 then begin
              r.(k)[rInds[w0]] += 1
              count_smooth.(k) += count
              rMask.(k)[rInds[w0]] += 1 ;debug
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
                r.(k)[rInds[w1]] += 3
                count_bclumpy.(k) += count
                rMask.(k)[rInds[w1]] += 1 ; debug
              endif
              
              ; only search for secondary parents for the non-primary remainder (this is mainly a doublecheck
              ; since we have already excluded every other possibility)
              wInSG = wInSG[wc]
              rInds = rInds[wc]
            
              qInds = value_locate(secSGids,parIDs[wInSG])
              w2 = where(secSGids[qInds] eq parIDs[wInSG],count)
              if count gt 0 then begin
                r.(k)[rInds[w2]] += 2
                count_sclumpy.(k) += count
                rMask.(k)[rInds[w2]] += 1 ; debug
              endif
            endif ; ncomp
          endif ; countInSG
        endif ; countNow
        
        ; check that all bracketed particles found accretion modes
        if min(r.(k)[w_now]) eq 0 then message,'Error: Not all done.'
        if (count_smooth.(k)+count_sclumpy.(k)+count_bclumpy.(k)) ne count_now then message,'Error: Totals fail to add up.'
        if max(rMask.(k)) gt 1 then message,'Error: Mask exceeds one.'
        
        ; any gas previously marked smooth, categorize as stripped if it is now in a non-parent halo
        smoothAccModes = [1,11]
        
        foreach smoothAccMode,smoothAccModes do begin       
          ; find either smooth only or smooth+recycled
          wSmooth = where(r.(k) eq smoothAccMode and trMinSnap.(k) le sP.snap,countSmooth)
        
          if countSmooth eq 0 then continue
          
          ; global match against subgroup ID list (this includes fuzz and unbound gas)
          calcMatch,galcatSub.(k)[wSmooth],groupcat_trIDs,galcat_ind,gc_ind,count=countNow
          
          if countNow eq 0 then continue
          
          ; calculate parent subgroup ID each gas particle belongs to
          parIDs = value_locate(subgroupOffsetTr,gc_ind)
          
          ; exclude those that aren't actually within the offsetLenType range
          diff = gc_ind - subgroupOffsetTr[parIDs]
          wInSG = where(diff lt subgroupLenTr[parIDs],countInSG)

          if countInSG ne countNow then message,'Error: Should not happen (gal).'
          
          ; if their subgroup ID is -not- the traced parent, assign to category 4. stripped
          rInds = wSmooth[galcat_ind[wInSG]]
          w0 = where(parIDs[wInSG] ne origParIDs.(k)[galcat_ind[wInSG]],count)
          if count gt 0 then begin
            r.(k)[rInds[w0]] = 4 + (smoothAccMode-1) ; addition is zero for accMode=1, 10 for recycled
            count_stripped.(k) += count
          endif
          
        endforeach ; smoothAccModes
        
      endfor ; k
          
      if m ne mt.minSnap then begin
        ; load mergerTree and move to Parent
        Parent = mergerTree(sP=sP)
        
        for k=0,n_tags(bracketSnap)-1 do begin
          ; change to parent IDs for the next snapshot
          w = where(origParIDs.(k) ne -1,count)
          if count eq 0 then message,'error'
          origParIDs.(k)[w] = Parent[origParIDs.(k)[w]]
          
          ; sanity check no IDs are -1 (we should only be processing the mergerTreeAdaptiveSubset)
          w = where(trMinSnap.(k) lt m-1,count)
          if count gt 0 then $
            if min(origParIDs.(k)[w]) lt 0 then message,'Error: Bad parent.'
        endfor
      endif
      
      ; count stripped/(smooth+stipped) and clumpy/total fraction
      w = where(r.accMode_gal eq 1,totSmoothGal)
      w = where(r.accMode_gal eq 4,totStrippedGal)
      w = where(r.accMode_gal eq 2 or r.accMode_gal eq 3,totClumpyGal)
      
      w = where(r.accMode_gmem eq 1,totSmoothGmem)
      w = where(r.accMode_gmem eq 4,totStrippedGmem)
      w = where(r.accMode_gmem eq 2 or r.accMode_gmem eq 3,totClumpyGmem)
      
      w = where(r.accMode_stars eq 1,totSmoothStars)
      w = where(r.accMode_stars eq 4,totStrippedStars)
      w = where(r.accMode_stars eq 2 or r.accMode_stars eq 3,totClumpyStars)
      
      fracStrippedGal   = float(totStrippedGal)   / (totSmoothGal+totStrippedGal)
      fracStrippedGmem  = float(totStrippedGmem)  / (totSmoothGmem+totStrippedGmem)
      fracStrippedStars = float(totStrippedStars) / (totSmoothStars+totStrippedStars)
      
      fracClumpyGal    = float(totClumpyGal)   / (totSmoothGal+totStrippedGal+totClumpyGal)
      fracClumpyGmem   = float(totClumpyGmem)  / (totSmoothGmem+totStrippedGmem+totClumpyGmem)
      fracClumpyStars  = float(totClumpyStars) / (totSmoothStars+totStrippedStars+totClumpyStars)
      
      count_bracketed += count_gal + count_gmem + count_stars
      cur_frac_found = float(count_bracketed) / nTotSearch
      
      ; count fractions again for recycled
      w = where(r.accMode_gal eq 11,totSmoothGal)
      w = where(r.accMode_gal eq 14,totStrippedGal)
      w = where(r.accMode_gal eq 12 or r.accMode_gal eq 13,totClumpyGal)
      
      w = where(r.accMode_gmem eq 11,totSmoothGmem)
      w = where(r.accMode_gmem eq 14,totStrippedGmem)
      w = where(r.accMode_gmem eq 12 or r.accMode_gmem eq 13,totClumpyGmem)
      
      w = where(r.accMode_stars eq 11,totSmoothStars)
      w = where(r.accMode_stars eq 14,totStrippedStars)
      w = where(r.accMode_stars eq 12 or r.accMode_stars eq 13,totClumpyStars) 
      
      fracStrippedGal2   = float(totStrippedGal)   / (totSmoothGal+totStrippedGal)
      fracStrippedGmem2  = float(totStrippedGmem)  / (totSmoothGmem+totStrippedGmem)
      fracStrippedStars2 = float(totStrippedStars) / (totSmoothStars+totStrippedStars)
      
      fracClumpyGal2    = float(totClumpyGal)   / (totSmoothGal+totStrippedGal+totClumpyGal)
      fracClumpyGmem2   = float(totClumpyGmem)  / (totSmoothGmem+totStrippedGmem+totClumpyGmem)
      fracClumpyStars2  = float(totClumpyStars) / (totSmoothStars+totStrippedStars+totClumpyStars)
      
      print,' ['+string(m,format='(i3)')+'] counts'+$
        ' smooth: ' +string(count_smooth.gal,format='(i5)')   +'.'+string(count_smooth.gmem,format='(i5)')+'.'+string(count_smooth.stars,format='(i5)')+$
        ' sclumpy: '+string(count_sclumpy.gal,format='(i5)')  +'.'+string(count_sclumpy.gmem,format='(i5)')+'.'+string(count_sclumpy.stars,format='(i5)')+$
        ' bclumpy: '+ string(count_bclumpy.gal,format='(i5)') +'.'+string(count_bclumpy.gmem,format='(i5)')+'.'+string(count_bclumpy.stars,format='(i5)')+$
        ' stripped: '+string(count_stripped.gal,format='(i5)')+'.'+string(count_stripped.gmem,format='(i5)')+'.'+string(count_stripped.stars,format='(i5)')
      print,' ['+string(m,format='(i3)')+'] fractions tot ('+string(cur_frac_found*100,format='(f4.1)')+'%)'+$
        ' stripped ('+string(fracStrippedGal*100,format='(f4.1)')+'% '+string(fracStrippedGmem*100,format='(f4.1)')+'% '+string(fracStrippedStars*100,format='(f4.1)')+'%)'+$
        ' clumpy ('+string(fracClumpyGal*100,format='(f4.1)')+'% '+string(fracClumpyGmem*100,format='(f4.1)')+'% '+string(fracClumpyStars*100,format='(f4.1)')+'%)'+$
        ' strippedRec ('+string(fracStrippedGal2*100,format='(f4.1)')+'% '+string(fracStrippedGmem2*100,format='(f4.1)')+'% '+string(fracStrippedStars2*100,format='(f4.1)')+'%)'+$
        ' clumpyRec ('+string(fracClumpyGal2*100,format='(f4.1)')+'% '+string(fracClumpyGmem2*100,format='(f4.1)')+'% '+string(fracClumpyStars2*100,format='(f4.1)')+'%)'
      
      ; free some memory for next load
      Parent  = !NULL
      all     = !NULL
      w_gal   = !NULL
      w_gmem  = !NULL
      w_stars = !NULL
    endfor ;m
    
    for k=0,n_tags(bracketSnap)-1 do begin
      ; verify we found an accretion mode for every gas particle
      if min(rMask.(k))   lt 1 then message,'Error: Not all found.'
      if min(r.(k))   lt 1 then message,'Error: Not all accMode set.'
    
      ; verify no strange accMode numbers
      w = where((r.(k) gt 4 and r.(k) lt 10) or r.(k) gt 14,count)
      if count gt 0 then message,'Error: Bad accMode.'
    endfor
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
  endif
  
  ; VELOCITY TRACERS case - will be similar to above since there could be multiple
  ; ---------------------
  if sP.trMCPerCell eq -1 then begin
    message,'todo'
  endif
  
  sP.snap = origSnap ; restore sP.snap
  return, r
end

; accModeInds(): subselect in the mtS/atS/traj subsets for a particular accretion mode

function accModeInds, at=at, sP=sP, accMode=accMode, mask=mask

  if n_elements(at) eq 0 or n_elements(accMode) eq 0 then message,'Error: Inputs'
  if ~sP.gfmWinds and accMode eq 'recycled' then message,'Error: Request recycled on non-winds run.'
  
  gal_w   = where(at.AccTime_gal[0,*] ne -1,count_gal)
  gmem_w  = where(at.AccTime_gmem[0,*] ne -1,count_gmem)
  stars_w = where(at.AccTime_stars[0,*] ne -1,count_stars)
  
  ; select on accretion mode by modifying gal_w, gmem_w, and stars_w
  if accMode ne 'all' then begin
    if sP.trMCPerCell lt 0 then message,'Error: accMode for tracerVel not yet.'
    am = accretionMode(sP=sP)
    
    if accMode eq 'smooth' then begin
      gal_w   = gal_w[where(am.accMode_gal eq 1,count_gal)]
      gmem_w  = gmem_w[where(am.accMode_gmem eq 1,count_gmem)]
      stars_w = stars_w[where(am.accMode_stars eq 1,count_stars)]
    endif
    if accMode eq 'smooth_rec' then begin
      gal_w   = gal_w[where(am.accMode_gal eq 1 or am.accMode_gal eq 11,count_gal)]
      gmem_w  = gmem_w[where(am.accMode_gmem eq 1 or am.accMode_gmem eq 11,count_gmem)]
      stars_w = stars_w[where(am.accMode_stars eq 1 or am.accMode_stars eq 11,count_stars)]
    endif
    
    if accMode eq 'bclumpy' then begin
      gal_w   = gal_w[where(am.accMode_gal eq 3,count_gal)]
      gmem_w  = gmem_w[where(am.accMode_gmem eq 3,count_gmem)]
      stars_w = stars_w[where(am.accMode_stars eq 3,count_stars)]
    endif
    if accMode eq 'sclumpy' then begin
      gal_w   = gal_w[where(am.accMode_gal eq 2,count_gal)]
      gmem_w  = gmem_w[where(am.accMode_gmem eq 2,count_gmem)]
      stars_w = stars_w[where(am.accMode_stars eq 2,count_stars)]
    endif
    
    if accMode eq 'clumpy' then begin
      gal_w   = gal_w[where(am.accMode_gal eq 2 or am.accMode_gal eq 3,count_gal)]
      gmem_w  = gmem_w[where(am.accMode_gmem eq 2 or am.accMode_gmem eq 3,count_gmem)]
      stars_w = stars_w[where(am.accMode_stars eq 2 or am.accMode_stars eq 3,count_stars)]
    endif
    if accMode eq 'clumpy_rec' then begin
      gal_w   = gal_w[where(am.accMode_gal eq 2 or am.accMode_gal eq 3 or $
                            am.accMode_gal eq 12 or am.accMode_gal eq 13,count_gal)]
      gmem_w  = gmem_w[where(am.accMode_gmem eq 2 or am.accMode_gmem eq 3 or $
                             am.accMode_gmem eq 12 or am.accMode_gmem eq 13,count_gmem)]
      stars_w = stars_w[where(am.accMode_stars eq 2 or am.accMode_stars eq 3 or $
                              am.accMode_stars eq 12 or am.accMode_stars eq 13,count_stars)]
    endif
    
    if accMode eq 'stripped' then begin
      gal_w   = gal_w[where(am.accMode_gal eq 4,count_gal)]
      gmem_w  = gmem_w[where(am.accMode_gmem eq 4,count_gmem)]
      stars_w = stars_w[where(am.accMode_stars eq 4,count_stars)]
    endif
    if accMode eq 'stripped_rec' then begin
      gal_w   = gal_w[where(am.accMode_gal eq 4 or am.accMode_gal eq 14,count_gal)]
      gmem_w  = gmem_w[where(am.accMode_gmem eq 4 or am.accMode_gmem eq 14,count_gmem)]
      stars_w = stars_w[where(am.accMode_stars eq 4 or am.accMode_stars eq 14,count_stars)]
    endif
    
    if accMode eq 'recycled' then begin ; recycled+any other mode
      gal_w   = gal_w[where(am.accMode_gal ge 10,count_gal)]
      gmem_w  = gmem_w[where(am.accMode_gmem ge 10,count_gmem)]
      stars_w = stars_w[where(am.accMode_stars ge 10,count_stars)]
    endif
    
    am = !NULL
  endif
  
  ; make mask to combine selection with additional selections, if requested
  if keyword_set(mask) then begin
    galMask   = bytarr(n_elements(at.accTime_gal[0,*]))
    gmemMask  = bytarr(n_elements(at.accTime_gmem[0,*]))
    starsMask = bytarr(n_elements(at.accTime_stars[0,*]))
    
    galMask[gal_w]     = 1B
    gmemMask[gmem_w]   = 1B
    starsMask[stars_w] = 1B
    r = {gal:gal_w,gmem:gmem_w,stars:stars_w,galMask:galmask,gmemMask:gmemMask,starsMask:starsMask}
  endif else begin
    r = {gal:gal_w,gmem:gmem_w,stars:stars_w}
  endelse
  
  return,r
end
