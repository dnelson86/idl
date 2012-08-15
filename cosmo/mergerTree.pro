; mergerTree.pro
; cosmological halo tracking / merger tree through time
; dnelson aug.2012

; mergerTree(): construct simplified merger tree for tracking halos/subhalos through time across snaps
; 
; makeNum : if specified, make the catalogs for the specified number of snaps back from sP.snap
;           (otherwise just make it for sP.snap and/or return existing catalog for sP.snap)

function mergerTree, sP=sP, makeNum=makeNum

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  partMatchFracTol = 0.6   ; 60% minimum match between particle members of the specified type
  massDiffFracTol  = 0.4   ; 40% agreement in total mass or better
  positionDiffTol  = 100.0 ; 200kpc maximum separation
  minSubgroupLen   = 100   ; do not try to match subgroups with less than N total particles (for speed)

  ptNum = partTypeNum('dm') ; use dark matter particles (only) for matching

  ; set minimum/ending snapshot (maxmimum/ending redshift)
  minSnap = sP.snap
  if keyword_set(makeNum) then minSnap = sP.snap - makeNum + 1
  
  ; earliest possible snapshot to create Parents for is the one after the first group catalog
  if minSnap le sP.groupCatRange[0] then minSnap = sP.groupCatRange[0] + 1

  ; set maximum/starting snapshot (minimum/starting redshift)
  maxSnap = sP.snap

  ; check for existing catalog
  if ~keyword_set(makeNum) then begin
    saveFilename = sP.derivPath+'mergerTree/mergerTree.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.sav'
    
    restore, saveFilename
    return, Parent
  endif
  
  ; if creating new catalog, start loop over snapshots
  for m=maxSnap,minSnap,-1 do begin
    sP.snap = m    
    ; if at maxSnap, load current group catalog
    if sP.snap eq maxSnap then gcCur = loadGroupCat(sP=sP,/readIDs)
    
    ; or, if not at maxSnap, move "previous" group catalog to current
    if sP.snap lt maxSnap then gcCur = gcPrev
    gcPrev = !NULL
    
    ; load "previous" group catalog (one snapshot back in time)
    sP.snap -= 1
    gcPrev = loadGroupCat(sP=sP,/readIDs)

    time=systime(/sec) ; start timer
    
    ; create/zero Parent to hold sgID of matching subgroup in previous snapshot
    Parent = lonarr(gcCur.nSubgroupsTot) - 1 ;-1 indicates no Parent found

    ; if at maxSnap, construct particle ID list for all subgroups in gcCur
    if m eq maxSnap then partIDs_cur = gcPIDList(gc=gcCur,select='all',partType=ptNum)
    
    ; or, if not at maxSnap, move "previous" particle ID list to current
    if m lt maxSnap then partIDs_cur = partIDs_prev
    partIDs_prev = !NULL
    
    ; construct particle ID list for all subgroups in gcPrev
    partIDs_prev = gcPIDList(gc=gcPrev,select='all',partType=ptNum)
    
    ; do global match between current and previous particle IDs
    match,partIDs_cur,partIDs_prev,cur_ind,prev_ind,count=matchCount,/sort
   
    ; IMPORTANT! rearrange cur_ind to be in the order of partIDs_cur (needed for while walk)
    prev_ind    = prev_ind[sort(cur_ind)]    
    cur_ind     = cur_ind[sort(cur_ind)]

    partIDs_cur = !NULL ; no longer used
    
    curSGEndCum = 0LL
    indCount = 0LL
    
    ; find parent of each current subgroup in "previous" catalog: loop over each current subgroup
    for i=0LL,gcCur.nSubgroupsTot-1 do begin
      
      ; find subset of matched indices of this current subgroup in all "previous" subgroups
      if gcCur.subgroupLenType[ptNum,i] eq 0 then continue
      
      ; METHOD A. one touch per element walk
      countTot = 0LL
      while (indCount+countTot) lt matchCount do begin
        if cur_ind[indCount+countTot] ge curSGEndCum+gcCur.subgroupLenType[ptNum,i] then break
        countTot += 1
      endwhile
      
      if countTot gt 0 then begin
        pInds = prev_ind[indCount:indCount+countTot-1]
        indCount += countTot ; increment progress counter for -matched- partIDs_cur
      endif else begin
        pInds = [0]
      endelse
      ; END (a)
      
      ; METHOD B. this works but probably slower on bigger group cats:
      ;pInds = where(cur_ind ge curSGEndCum and $
      ;              cur_ind lt curSGEndCum+gcCur.subgroupLenType[ptNum,i],countTot)
      ;pInds = prev_ind[pInds] ; want the indices of partIDs_prev, not of cur/prev_ind or partIDs_cur
      ; END (b)
      
      pInds = pInds[sort(pInds)] ; sort ascending for the j walk
      
      curSGEndCum += gcCur.subgroupLenType[ptNum,i] ; increment progress counter for -all- partIDs_cur

      if countTot eq 0 then continue
      if gcCur.subgroupLen[i] lt minSubgroupLen then continue
      
      ; create an array to store the number of particles found in each "previous" subgroup
      prevSGIndex  = 0LL
      prevSGEndCum = ulong(gcPrev.subgroupLenType[ptNum,0])
      partCount    = fltarr(gcPrev.nSubgroupsTot)
      
      ; walk through matched indices in cur_ind (and so prev_ind) and assign counts to each 
      ; subgroup owner in the "previous" snapshot
      for j=0LL,countTot-1 do begin
        if pInds[j] lt prevSGEndCum then begin
          ; this particle is within the subgroupOffset for this prevSGIndex
          partCount[prevSGIndex] += 1.0
        endif else begin
          ; this particle is not in this prevSGIndex subgroup, move on to the next and repeat this j
          prevSGIndex += 1
          prevSGEndCum += gcPrev.subgroupLenType[ptNum,prevSGIndex]
          j -= 1
        endelse
      endfor
      
      ; convert particle counts to fraction of current subgroup's particles in each previous subgroup
      partCount /= gcCur.subgroupLenType[ptNum,i]
      
      ; find "previous" subgroups with first and second highest fractions
      maxFrac = max(partCount,max_index)

      ; enforce tolerances and save parent
      massDiffFrac = abs(gcCur.subgroupMass[i] - gcPrev.subgroupMass[max_index]) / gcCur.subgroupMass[i]

      xyzDist = gcCur.subgroupPos[*,i]-gcPrev.subgroupPos[*,max_index]
      correctPeriodicDistVecs, xyzDist, sP=sP
      positionDiff = (sqrt( xyzDist[0]*xyzDist[0] + xyzDist[1]*xyzDist[1] + xyzDist[2]*xyzDist[2] ) )[0]

      if massDiffFrac lt massDiffFracTol and positionDiff lt positionDiffTol and $
         maxFrac gt partMatchFracTol then Parent[i] = max_index
         
      ; DEBUG: what matches are we loosing?
      ;if massDiffFrac gt massDiffFracTol or positionDiff gt positionDiffTol or $
      ;  maxFrac lt partMatchFracTol and gcCur.subgroupMass[i] ge 1.0 then $
      ;  print,codeMassToLogMsun(gcCur.subgroupMass[i]),maxFrac,massDiffFrac,positionDiff
         
      ; DEBUG: verify
      ;debug_curIDs  = gcCur.IDs [gcCur.subgroupOffsetType[ptNum,i]:$
      ;                           gcCur.subgroupOffsetType[ptNum,i]+gcCur.subgroupLenType[ptNum,i]-1]
      ;debug_prevIDs = gcPrev.IDs[gcPrev.subgroupOffsetType[ptNum,max_index]:$
      ;                           gcPrev.subgroupOffsetType[ptNum,max_index]+$
      ;                           gcPrev.subgroupLenType[ptNum,max_index]-1]
      ;match,debug_curIDs,debug_prevIDs,debug_ind1,debug_ind2,count=debug_count
      ;debug_frac = float(debug_count)/n_elements(debug_curIDs) ;TODO: cur or prev
      ;if abs(debug_frac-maxFrac) gt 1e-6 then message,'DEBUG FAIL'
      ;partCount[max_index] = 0.0
      ;maxFrac2 = max(partCount,max_index2)
      ;print,gcCur.subgroupGrnr[i],gcCur.subgroupLen[i],gcCur.subgroupLenType[ptNum,i],$
      ;      massDiffFrac,positionDiff,maxFrac;,maxFrac2,debug_frac
    endfor
    
    ; count number of parents found and calculate delta(age of universe) back so far
    w = where(Parent ne -1,count)
    wMS = where(gcCur.subgroupLen ge minSubgroupLen,countMS)
    
    cur_z     = snapnumToRedshift(sP=sP)
    delta_z   = cur_z - sP.redshift
    delta_age = (redshiftToAgeFlat(sP.redshift) - redshiftToAgeFlat(cur_z))*1000 ;Myr
    
    ; save Parent ("merger tree") for this snapshot
    saveFilename = sP.derivPath+'mergerTree/mergerTree.'+sP.savPrefix+str(sP.res)+'.'+str(m)+'.sav'
    if file_test(saveFilename) then begin
      print,'SKIP : '+strmid(saveFilename,strlen(sp.derivPath))
    endif else begin
      save,Parent,filename=saveFilename
      print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))+' (Matched: '+$
        string(float(count)/countMS*100,format='(f5.1)')+'%, '+$
        string(systime(/sec)-time,format='(f5.1)')+' sec, Dz = '+$
        string(delta_z,format='(f5.3)')+' Dt = '+string(delta_age,format='(f6.1)')+' Myr)'
    endelse
  endfor

end

; mergerTreeAdaptiveSubset():
;                     walk the merger tree from sP.snap back to start of the group catalogs and calculate
;                     the snapshot to which each halo can be tracked in time. also calculate
;                     a few properties for these subgroups back in time (r_vir,T_vir,position) and 
;                     the subset of the galaxycat at maxSnap that we can track back using this selection

function mergerTreeAdaptiveSubset, sP=sP, verbose=verbose

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  maxSnap = sP.snap
  minSnap = sP.groupCatRange[0] + 1 ;z=6
  
  smoothKer = 3 ; 1,3,5 number of snapshots
  units = getUnits()
  
  ; set saveFilename and check for existence
  saveFilename = sP.derivPath + 'mTreeAdaSub.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(maxSnap)+'-'+str(minSnap)+'.sK'+str(smoothKer)+'.sav'  

  if file_test(saveFilename) then begin
    restore, saveFilename
    return, r
  endif

  print,'Walking merger tree for adaptive halo selection...'

  numBack = maxSnap - minSnap
  
  for i=0,numBack+floor(smoothKer/2.0)-1 do begin
    ; load
    h = loadSnapshotHeader(sP=sP)
    gc = loadGroupCat(sP=sP,/skipIDs)
    subgroupCen = subgroupPosByMostBoundID(sP=sP)

    ; on first snapshot select primary halos and allocate arrays
    if i eq 0 then begin
      gcIDs = gcIDList(gc=gc,select='pri')

      gcIDPlace = lindgen(n_elements(gcIDs)) ; initial placement in ascending order of ID
      
      times     = fltarr(numBack+floor(smoothKer/2.0))
      hPos      = fltarr(numBack+floor(smoothKer/2.0),3,n_elements(gcIDPlace))
      hVel      = fltarr(numBack+floor(smoothKer/2.0),3,n_elements(gcIDPlace))
      hMass     = fltarr(numBack+floor(smoothKer/2.0),n_elements(gcIDPlace))
      hVirRad   = fltarr(numBack+floor(smoothKer/2.0),n_elements(gcIDPlace))
      hVirTemp  = fltarr(numBack+floor(smoothKer/2.0),n_elements(gcIDPlace))
      noParMask = bytarr(n_elements(gcIDPlace))
      hMinSnap  = intarr(n_elements(gcIDPlace)) - 1
    endif
  
    ; which gcIDs still valid? only new "no parents"
    ;wPar = where(gcIDs ne -1,ncomp=count,comp=wNoPar)
    wNoPar = where(gcIDs eq -1 and noParMask eq 0B)
    
    noParMask[wNoPar] = 1B
    hMinSnap[wNoPar] = sP.snap + 1 ; mark end of tracking for halos with no valid parent at this snapshot
    gcIDs[wNoPar] = -1 ; from this point on, will contain incorrect values, must use noParMask to select
    
    wPar = where(hMinSnap eq -1,nleft)
    wPar10 = where(hMinSnap eq -1 and hMass[0,*] ge 1.0,nleft10)
    tot10 = where(hMass[0,*] ge 1.0,ntot10)
    
    if keyword_set(verbose) then $
      print,' ['+string(sP.snap,format='(i3)')+'] remaining: '+string(nleft,format='(i5)')+' ('+$
        string(float(nleft)/n_elements(gcIDs)*100,format='(f4.1)')+'%) massive: '+string(nleft10,format='(i5)')+' ('+$
        string(float(nleft10)/ntot10*100,format='(f4.1)')+'%) '
        
    ; store only currently valid gcIDs
    w = where(gcIDs ne -1,count)
    if count eq 0 then message,'error'
    
    times[i] = h.time
    hPos[i,*,gcIDPlace[w]]   = subgroupCen[*,gcIDs[w]]
    hVel[i,*,gcIDPlace[w]]   = gc.subgroupVel[*,gcIDs[w]]
    hMass[i,gcIDPlace[w]]    = gc.subgroupMass[gcIDs[w]]
    hVirRad[i,gcIDPlace[w]]  = gc.group_r_crit200[gc.subgroupGrNr[gcIDs[w]]]
    hVirTemp[i,gcIDPlace[w]] = alog10(codeMassToVirTemp(gc.subgroupMass[gcIDs[w]],sP=sP))

    ; load mergerTree and change subgroup IDs to parents at the prior snapshot
    Parent = mergerTree(sP=sP)
    
    gcIDs[w] = Parent[gcIDs[w]] ; change to parent IDs

    sP.snap -= 1
  endfor
  
  Parent = !NULL
  gcIDs  = !NULL
  gcIDPlace = !NULL
  noParMask = !NULL
  subgroupCen = !NULL
  
  hPosSm = hPos
  hVelSm = hVel
  
  ; created smoothed estimates of pos(t), vel(t), mass(t), r_vir(t) and T_vir(t) for halos
  for i=0,n_elements(hMass[0,*])-1 do begin
    endInd = maxSnap - hMinSnap[i]
    if endInd eq 0 or endInd lt smoothKer or hMinSnap[i] eq -1 then continue
    
    hPosSm[0:endInd,0,i]   = smooth(reform(hPos[0:endInd,0,i]),smoothKer)
    hPosSm[0:endInd,1,i]   = smooth(reform(hPos[0:endInd,1,i]),smoothKer)
    hPosSm[0:endInd,2,i]   = smooth(reform(hPos[0:endInd,2,i]),smoothKer)
    hVelSm[0:endInd,0,i]   = smooth(reform(hVel[0:endInd,0,i]),smoothKer)
    hVelSm[0:endInd,1,i]   = smooth(reform(hVel[0:endInd,1,i]),smoothKer)
    hVelSm[0:endInd,2,i]   = smooth(reform(hVel[0:endInd,2,i]),smoothKer)
    
    hMass[0:endInd,i]    = smooth(reform(hMass[0:endInd,i]),smoothKer)
    hVirRad[0:endInd,i]  = smooth(reform(hVirRad[0:endInd,i]),smoothKer)
    hVirTemp[0:endInd,i] = smooth(reform(hVirTemp[0:endInd,i]),smoothKer)
  endfor
  
  ; load galaxy/group member catalogs at zMin for gas ids to search for
  sP.snap = maxSnap
  galcat = galaxyCat(sP=sP)
  
  ; replicate parent IDs (of PRIMARY/parent)
  gc = loadGroupCat(sP=sP,/skipIDs)
  galcatIDList = gcIDList(gc=gc,select='pri') ; this is the starting gcIDs from above
  
  ; replicated galcat parent IDs for each galcatSub member
  gcIndOrig = galCatRepParentIDs(galcat=galcat,gcIDList=galcatIDList)
  
  ; want to use these parent IDs to access hVirRad,etc so compact the same way (ascending ID->index)
  placeMap = getIDIndexMap(galcatIDList,minid=minid)
  gcIndOrig.gal   = placeMap[gcIndOrig.gal-minid]
  gcIndOrig.gmem  = placeMap[gcIndOrig.gmem-minid]
  gcIndOrig.stars = placeMap[gcIndOrig.stars-minid]
  placeMap = !NULL
  
  ; get the subset of the galcat indices in sgIDList
  galcatSub = galcatINDList(galcat=galcat,gcIDList=galcatIDList)
  
  r = {galcatIDList:galcatIDList,gcIndOrig:gcIndOrig,galcatSub:galcatSub,$
       hPos:hPos,hVel:hVel,hPosSm:hPosSm,hVelSm:hVelSm,$
       hVirRad:hVirRad,hVirTemp:hVirTemp,hMass:hMass,hMinSnap:hMinSnap,$
       times:times,minSnap:minSnap,maxSnap:maxSnap,smoothKer:smoothKer}

  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))

  return,r
end

; mergerTreeSubset(): walk the merger tree from sP.snap back to some min snap and return the subset of the 
;                     subgroups at sP.snap that are well-connected all the way back. also calculate
;                     a few properties for these subgroups back in time (r_vir,T_vir,position) and 
;                     the subset of the galaxycat at maxSnap that we can track back using this selection
;
; NOTE: minSnap is decided here (change will require recalculation of all accretionTimes/accretionModes)

function mergerTreeSubset, sP=sP, verbose=verbose

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; --------------------------------------------------------------------------------------------------
  ; switch to new adaptive mtS
  ;print,'Note: mergerTreeAdaptiveSubset'
  mt = mergerTreeAdaptiveSubset(sP=sP,verbose=verbose)
  return,mt
  ; --------------------------------------------------------------------------------------------------
  
  ; config
  maxSnap = sP.snap
  minSnap = redshiftToSnapnum(4.0,sP=sP)
  
  smoothKer = 3 ; 1,3,5 number of snapshots
  units = getUnits()
  
  ; set saveFilename and check for existence
  saveFilename = sP.derivPath + 'mTreeSub.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(maxSnap)+'-'+str(minSnap)+'.sK'+str(smoothKer)+'.sav'  

  if file_test(saveFilename) then begin
    restore, saveFilename
    return, r
  endif

  print,'Walking merger tree for halo selection...'

  numBack = maxSnap - minSnap
  
  for i=0,numBack+floor(smoothKer/2.0)-1 do begin
    ; load
    h = loadSnapshotHeader(sP=sP)
    gc = loadGroupCat(sP=sP,/skipIDs)
    subgroupCen = subgroupPosByMostBoundID(sP=sP)

    ; on first snapshot select primary halos and allocate arrays
    if i eq 0 then begin
      gcIDs = gcIDList(gc=gc,select='pri')

      gcIDPlace = lindgen(n_elements(gcIDs)) ; initial placement in ascending order of ID
      
      times     = fltarr(numBack+floor(smoothKer/2.0))
      hPos      = fltarr(numBack+floor(smoothKer/2.0),3,n_elements(gcIDPlace))
      hVel      = fltarr(numBack+floor(smoothKer/2.0),3,n_elements(gcIDPlace))
      hMass     = fltarr(numBack+floor(smoothKer/2.0),n_elements(gcIDPlace))
      hVirRad   = fltarr(numBack+floor(smoothKer/2.0),n_elements(gcIDPlace))
      hVirTemp  = fltarr(numBack+floor(smoothKer/2.0),n_elements(gcIDPlace))
      noParMask = bytarr(n_elements(gcIDPlace))
    endif
  
    ; which gcIDs still valid?
    wPar = where(gcIDs ne -1,ncomp=count,comp=wNoPar)
    
    noParMask[wNoPar] = 1B ; mark halos with no parents
    gcIDs[wNoPar] = 0 ; from this point on, will contain incorrect values, must use noParMask to select
    
    wPar = where(noParMask eq 0B,nleft)
    wPar10 = where(noParMask eq 0B and hMass[0,*] ge 1.0,nleft10)
    tot10 = where(hMass[0,*] ge 1.0,ntot10)
    
    if keyword_set(verbose) then $
      print,' ['+string(sP.snap,format='(i3)')+'] remaining: '+string(nleft,format='(i5)')+' ('+$
        string(float(nleft)/n_elements(gcIDs)*100,format='(f4.1)')+'%) massive: '+string(nleft10,format='(i5)')+' ('+$
        string(float(nleft10)/ntot10*100,format='(f4.1)')+'%)
     
    ; store
    times[i] = h.time
    hPos[i,*,gcIDPlace]   = subgroupCen[*,gcIDs]
    hVel[i,*,gcIDPlace]   = gc.subgroupVel[*,gcIDs]
    hMass[i,gcIDPlace]    = gc.subgroupMass[gcIDs]
    hVirRad[i,gcIDPlace]  = gc.group_r_crit200[gc.subgroupGrNr[gcIDs]]
    hVirTemp[i,gcIDPlace] = alog10(codeMassToVirTemp(gc.subgroupMass[gcIDs],sP=sP))

    ; load mergerTree and change subgroup IDs to parents at the prior snapshot
    Parent = mergerTree(sP=sP)
    
    gcIDs = Parent[gcIDs] ; change to parent IDs
    
    sP.snap -= 1
  endfor
  
  ; use subgroup ID list comprised only of those with good parent histories (indices in gcIDs)
  galcatIDList = where(noParMask eq 0B)
  
  Parent = !NULL
  gcIDs  = !NULL
  gcIDPlace = !NULL
  noParMask = !NULL
  subgroupCen = !NULL
  
  ; created smoothed estimates of pos(t), vel(t), mass(t),  r_vir(t) and T_vir(t) for only those kept halos
  hPos     = hPos[*,*,galcatIDList] ; don't smooth positions or velocities here
  hVel     = hVel[*,*,galcatIDList]
  hMass    = smooth(hMass[*,galcatIDList],[1,smoothKer])
  hVirRad  = smooth(hVirRad[*,galcatIDList],[1,smoothKer])
  hVirTemp = smooth(hVirTemp[*,galcatIDList],[1,smoothKer])
  print,'smooth ker warning'
  
  ; load galaxy/group member catalogs at zMin for gas ids to search for
  sP.snap = maxSnap
  galcat = galaxyCat(sP=sP)
  
  ; replicate parent IDs (of PRIMARY/parent)
  gc = loadGroupCat(sP=sP,/skipIDs)
  priParentIDs = gcIDList(gc=gc,select='pri') ; this is the starting gcIDs from above
  
  ; transform galcatIDList from indices to IDs
  galcatIDList = priParentIDs[galcatIDList]
  
  ; replicated galcat parent IDs for each galcatSub member
  gcIndOrig = galCatRepParentIDs(galcat=galcat,priParentIDs=priParentIDs,gcIDList=galcatIDList)
  
  ; want to use these parent IDs to access hVirRad,etc so compact the same way (ascending ID->index)
  placeMap = getIDIndexMap(galcatIDList,minid=minid)
  gcIndOrig.gal = placeMap[gcIndOrig.gal-minid]
  gcIndOrig.gmem = placeMap[gcIndOrig.gmem-minid]
  placeMap = !NULL
  
  ; get the subset of the galcat indices in sgIDList
  galcatSub = galcatINDList(galcat=galcat,gcIDList=galcatIDList)
  
  r = {galcatIDList:galcatIDList,gcIndOrig:gcIndOrig,galcatSub:galcatSub,$
       hPos:hPos,hVel:hVel,hVirRad:hVirRad,hVirTemp:hVirTemp,hMass:hMass,$
       times:times,minSnap:minSnap,maxSnap:maxSnap,smoothKer:smoothKer}
       
  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
       
  return,r
end

; mergerTreeRepParentIDs(): wrapper which returns mt.gcIndOrig for SPH simulations and the similarly
;                           replicated galcat parent ID list for tracer simulations

function mergerTreeRepParentIDs, mt=mt, galcat=galcat, sP=sP, compactMtS=compactMtS, $
                                 trids_gal=galcat_gal_trids, trids_gmem=galcat_gmem_trids, $ ; outputs
                                 trids_stars=galcat_stars_trids, gc_gal_cc=galcat_gal_cc, $ ; outputs
                                 gc_gmem_cc=galcat_gmem_cc, gc_stars_cc=galcat_stars_cc ; outputs

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function cosmoTracerChildren, galaxyCat

  if n_elements(sP) eq 0 then message,'Error: sP required.'
  if n_elements(mt) eq 0 then mt = mergerTreeSubset(sP=sP)
  if n_elements(galcat) eq 0 then galcat = galaxyCat(sP=sP)

  if sP.trMCPerCell eq 0 then begin
    gcIndOrig = galCatRepParentIDs(galcat=galcat,gcIDList=mt.galcatIDList)
  
    ; want to use these parent IDs to access hVirRad,etc so compact the same way (ascending ID->index)
    if keyword_set(compactMtS) then begin
      placeMap = getIDIndexMap(mt.galcatIDList,minid=minid)
      gcIndOrig.gal   = placeMap[gcIndOrig.gal-minid]
      gcIndOrig.gmem  = placeMap[gcIndOrig.gmem-minid]
      gcIndOrig.stars = placeMap[gcIndOrig.stars-minid]
      placeMap = !NULL
    endif
    
    return,gcIndOrig
  endif

  if sP.trMCPerCell gt 0 then begin
  ; note: all this ID loading seems unnecessary, since ids_gal=galcat.galaxyIDs[mt.galcatSub.gal] etc
  ;  ; load gas ids
  ;  gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')

  ;  ; match galcat IDs to gas_ids
  ;  match,galcat.galaxyIDs[mt.galcatSub.gal],gas_ids,galcat_ind,ids_gal_ind,count=countGal,/sort
  ;  ids_gal = gas_ids[ids_gal_ind[sort(galcat_ind)]]

  ;  match,galcat.groupmemIDs[mt.galcatSub.gmem],gas_ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
  ;  ids_gmem = gas_ids[ids_gmem_ind[sort(galcat_ind)]]
    
  ;  gas_ids = !NULL
    
  ;  ; load star ids and match
  ;  star_ids = loadSnapshotSubset(sP=sP,partType='stars',field='ids')
    
  ;  match,galcat.stellarIDs[mt.galcatSub.stars],star_ids,galcat_ind,ids_stars_ind,count=countStars,/sort
  ;  ids_stars = star_ids[ids_stars_ind[sort(galcat_ind)]]
    
  ;  star_ids = !NULL
 ;   if countGal ne n_elements(mt.galcatSub.gal) or countGmem ne n_elements(mt.galcatSub.gmem) then $
 ;     message,'Error: Check.'
 
    ids_gal   = galcat.galaxyIDs[mt.galcatSub.gal] ; new
    ids_gmem  = galcat.groupmemIDs[mt.galcatSub.gmem]
    ids_stars = galcat.stellarIDs[mt.galcatSub.stars]
    
    ; locate tracer children (indices) of gas id subsets
    galcat_gal_trids   = cosmoTracerChildren(sP=sP, /getInds, gasIDs=ids_gal, child_counts=galcat_gal_cc)
    galcat_gmem_trids  = cosmoTracerChildren(sP=sP, /getInds, gasIDs=ids_gmem, child_counts=galcat_gmem_cc)
    galcat_stars_trids = cosmoTracerChildren(sP=sP, /getInds, starIDs=ids_stars, child_counts=galcat_stars_cc)
     
    ids_gal   = !NULL
    ids_gmem  = !NULL
    ids_stars = !NULL
    
    ; convert tracer children indices to tracer IDs at this zMin if we are returning them
    if keyword_set(galcat_gal_trids) then begin
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
      galcat_gal_trids   = tr_ids[galcat_gal_trids]
      galcat_gmem_trids  = tr_ids[galcat_gmem_trids]
      galcat_stars_trids = tr_ids[galcat_stars_trids]
      tr_ids   = !NULL
    endif
    
    ; create a gcIndOrig for the tracers
    gcIndOrigTr = galCatRepParentIDs(galcat=galcat,gcIDList=mt.galcatIDList,$
                                     child_counts={gal:galcat_gal_cc,gmem:galcat_gmem_cc,stars:galcat_stars_cc}) 
                  
    ; want to use these parent IDs to access hVirRad,etc so compact the same way (ascending ID->index)
    if keyword_set(compactMtS) then begin
      placeMap = getIDIndexMap(mt.galcatIDList,minid=minid)
      gcIndOrigTr.gal   = placeMap[gcIndOrigTr.gal-minid]
      gcIndOrigTr.gmem  = placeMap[gcIndOrigTr.gmem-minid]
      gcIndOrigTr.stars = placeMap[gcIndOrigTr.stars-minid]
      placeMap = !NULL
    endif
    
    return,gcIndOrigTr
  endif
  
  if sP.trMCPerCell eq -1 then begin
    ; load gas ids
    gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')

    ; match galcat IDs to gas_ids
    match,galcat.galaxyIDs[mt.galcatSub.gal],gas_ids,galcat_ind,ids_gal_ind,count=countGal,/sort
    inds_gal = ids_gal_ind[sort(galcat_ind)]
    
    match,galcat.groupmemIDs[mt.galcatSub.gmem],gas_ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
    inds_gmem = ids_gmem_ind[sort(galcat_ind)]
    
    ; no stars
    
    if countGal ne n_elements(mt.galcatSub.gal) or countGmem ne n_elements(mt.galcatSub.gmem) then $
      message,'Error: Failed to locate all of galcat in gas_ids (overflow64?).'

    gas_ids = !NULL
    
    ; locate tracer children (indices) of gas id subsets
    galcat_gal_trids  = cosmoTracerVelChildren(sP=sP,/getInds,gasInds=inds_gal,child_counts=galcat_gal_cc)
    galcat_gmem_trids = cosmoTracerVelChildren(sP=sP,/getInds,gasInds=inds_gmem,child_counts=galcat_gmem_cc)
    
    inds_gal  = !NULL
    inds_gmem = !NULL
    
    ; convert tracer children indices to tracer IDs at this zMin if we are returning them
    if keyword_set(galcat_gal_trids) then begin
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerVel',field='ids')
      galcat_gal_trids  = tr_ids[galcat_gal_trids]
      galcat_gmem_trids = tr_ids[galcat_gmem_trids]
      tr_ids   = !NULL
    endif
    
    ; create a gcIndOrig for the tracers
    gcIndOrigTr = galCatRepParentIDs(galcat=galcat,gcIDList=mt.galcatIDList,$
                                     child_counts={gal:galcat_gal_cc,gmem:galcat_gmem_cc}) 
                  
    ; want to use these parent IDs to access hVirRad,etc so compact the same way (ascending ID->index)
    if keyword_set(compactMtS) then begin
      placeMap = getIDIndexMap(mt.galcatIDList,minid=minid)
      gcIndOrigTr.gal = placeMap[gcIndOrigTr.gal-minid]
      gcIndOrigTr.gmem = placeMap[gcIndOrigTr.gmem-minid]
      placeMap = !NULL
    endif
    
    return,gcIndOrigTr
  endif

end

; mergerTreeINDList(): return a list of indices into the mergerTreeSubset/accretionTimes/Traj catalog 
;                      for a subset of the members defined by the subgroup ID list gcIDList

function mergerTreeINDList, sP=sP, galcat=galcat, mt=mt, gcIDList=gcIDList

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function galaxyCat

  ; load galaxy cat and mergerTreeSubset if necessary
  if ~keyword_set(galcat) then galcat = galaxyCat(sP=sP)
  if ~keyword_set(mt) then mt = mergerTreeSubset(sP=sP)
  
  ; verify that all the requested subgroup IDs are part of the mergerTreeSubset
  match,mt.galcatIDList,gcIDList,ind1,ind2,count=count
  if count ne n_elements(gcIDList) then message,'Error: Not all gcIDs found in mergerTreeSubset.'
  
  ; make mask for requested subgroup IDs
  gcIDMask = bytarr(n_elements(galcat.galaxyLen))
  if keyword_set(gcIDList) then gcIDMask[gcIDList] = 1B  
  if ~keyword_set(gcIDList) then gcIDMask[*] = 1B
  
  ; normal indices return (handle zeros which will happen for individual halo ID requests at low masses)
  ; note: must handle the -1 return inside the calling function
  tots = { gal : total(galcat.galaxyLen[gcIDList],/int), gmem : total(galcat.groupmemLen[gcIDList],/int) }
  
  if tots.gal gt 0 and tots.gmem gt 0 then r = {gal  : ulonarr(tots.gal) , gmem : ulonarr(tots.gmem) }
  if tots.gal gt 0 and tots.gmem eq 0 then r = {gal  : ulonarr(tots.gal) , gmem : -1 }
  if tots.gal eq 0 and tots.gmem gt 0 then r = {gal  : -1 , gmem : ulonarr(tots.gmem) }
  
  offsetGal  = 0L
  offsetGmem = 0L
  
  offsetGal_mt  = 0L
  offsetGmem_mt = 0L
  
  ; (1) make list for gas cells/particles
  foreach gcID,mt.galcatIDList do begin
    ; galaxy
    if galcat.galaxyLen[gcID] gt 0 and gcIDMask[gcID] eq 1B then begin
      galInds    = ulindgen(galcat.galaxyLen[gcID]) + offsetGal_mt
      r.gal[offsetGal:offsetGal+galcat.galaxyLen[gcID]-1] = galInds
      offsetGal += galcat.galaxyLen[gcID]
    endif
    
    ; group member
    if galcat.groupmemLen[gcID] gt 0 and gcIDMask[gcID] eq 1B then begin
      gmemInds    = ulindgen(galcat.groupmemLen[gcID]) + offsetGmem_mt
      r.gmem[offsetGmem:offsetGmem+galcat.groupmemLen[gcID]-1] = gmemInds
      offsetGmem += galcat.groupmemLen[gcID]
    endif
    
    offsetGal_mt  += galcat.galaxyLen[gcID]
    offsetGmem_mt += galcat.groupmemLen[gcID]
  endforeach
  
  ; (2) make list including child counts if simulation has a tracer type
  if sP.trMCPerCell eq 0 then return,r
  
  ; load child counts
  maxt_gal = maxTemps(sP=sP,/loadAllTRGal)
  child_counts_gal = maxt_gal.child_counts
  maxt_gal = !NULL
  
  maxt_gmem = maxTemps(sP=sP,/loadAllTRGmem)
  child_counts_gmem   = maxt_gmem.child_counts
  maxt_gmem = !NULL
    
  if n_elements(child_counts_gal) ne total(galcat.galaxyLen,/int) or $
     n_elements(child_counts_gmem) ne total(galcat.groupmemLen,/int) then $
     message,'Error: Child_counts gal/gmem should have same size as full galcat subset.'

  if total(child_counts_gal,/int) gt 2e9 then stop ; consider lon64/removing /int
  if total(child_counts_gmem,/int) gt 2e9 then stop

  rcc = { gal  : ulonarr(total(child_counts_gal[r.gal],/int))  ,$
          gmem : ulonarr(total(child_counts_gmem[r.gmem],/int)) }
       
  offsetGal  = 0L
  offsetGmem = 0L
  
  offsetGal_all  = 0UL
  offsetGmem_all = 0UL
  
  for gcID=0UL,n_elements(galcat.galaxyLen)-1 do begin
    ; galaxy
    if galcat.galaxyLen[gcID] gt 0 then begin
      ; calculate total number of children in this subgroup
      tot_children_gal  = total(child_counts_gal[galcat.galaxyOff[gcID]:$
                                                 galcat.galaxyOff[gcID]+galcat.galaxyLen[gcID]-1],/int)
      ; add indices only for specified galaxy IDs
      if gcIDMask[gcID] then begin
        ; calculate place and store indices
        galInds = ulindgen(tot_children_gal) + offsetGal_all
        rcc.gal[offsetGal:offsetGal+tot_children_gal-1] = galInds
  
        offsetGal += tot_children_gal
      endif
      
      ; add child counts to indice array offset
      offsetGal_all  += tot_children_gal
    endif      
                  
    ; group member
    if galcat.groupmemLen[gcID] gt 0 then begin
      ; calculate total number of children in this subgroup             
      tot_children_gmem = total(child_counts_gmem[galcat.groupmemOff[gcID]:$
                                                  galcat.groupmemOff[gcID]+galcat.groupmemLen[gcID]-1],/int)
            
      ; add indices only for specified group member IDs
      if gcIDMask[gcID] then begin
        ; calculate place and store indices
        gmemInds     = ulindgen(tot_children_gmem) + offsetGmem_all
        rcc.gmem[offsetGmem:offsetGmem+tot_children_gmem-1] = gmemInds
        
        offsetGmem += tot_children_gmem
      endif
      
      offsetGmem_all += tot_children_gmem
    endif
    
  endfor
  
  return,rcc
  
end

; trackHaloPosition(): return a new position for some halo gcID at an earlier snapshot using mergerTree
    
function trackHaloPosition, sP=sP, gcID=gcID, endSnap=endSnap
    
    ; try to trace back to this snapshot for new halo center position
    failSnap = -1
    sgCenters = fltarr(3,sP.snap-endSnap+1) ; xyz ckpc
    times     = fltarr(sP.snap-endSnap+1) ; Gyr age
    
    gcIDcur = gcID ; starting halo ID
    startSnap = sP.snap ; starting snapshot number
    
    for m=startSnap,endSnap,-1 do begin
      ;print,m,gcIDcur,failSnap
      sP.snap = m
      
      ; haven't failed yet, keep recording positions
      if failSnap eq -1 then begin
        sgcen = subgroupPosByMostBoundID(sP=sP)
        sgCenters[*,startSnap-m] = sgcen[*,gcIDcur]
        
        ; record times for possible extrapolation
        h = loadSnapshotHeader(sP=sP)
        times[startSnap-m] = redshiftToAgeFlat(1/h.time-1)
      
        ; load most massive progenitor and move index
        Parent = mergerTree(sP=sP)
        
        if Parent[gcIDcur] eq -1 then failSnap = sP.snap
        gcIDcur = Parent[gcIDcur]
      endif
    endfor
    
    if failSnap eq -1 then begin
      ; successfully tracked back, load new subgroup centers and use
      haloPos = sgCenters[*,-1]
    endif else begin
      ; failed to find a good parent history
      if failSnap eq startSnap or failSnap eq startSnap-1 then begin
        ; zero or one parents (one or two centers), just use starting position
        haloPos = sgCenters[*,0]
      endif else begin
        ; at least two parents (three positions), do an extrapolation to the ending snapshot
        print,'untested'
        sP.snap = endSnap
        h = loadSnapshotHeader(sP=sP)
        timeEnd = redshiftToAgeFlat(1/h.time-1)
        
        inds = where(times ne 0.0,count)
        if count lt 3 then message,'error'
        
        xPos = interpol(sgCenters[0,w],times[w],timeEnd)
        yPos = interpol(sgCenters[1,w],times[w],timeEnd)
        zPos = interpol(sgCenters[2,w],times[w],timeEnd)
        haloPos = [xPos,yPos,zPos]
      endelse
    endelse
    
    return,haloPos
end

; plotHaloTracking(): plot fraction of successful halo adaptive tracking vs redshift

pro plotHaloTracking, sP=sP

  minMasses = [10.0,10.5,11.0,11.5] ; z=2 masses in log(msun)
  
  mt = mergerTreeAdaptiveSubset(sP=sP,/verbose)
  numSnaps = mt.maxSnap - mt.minSnap
  
  ; arrays
  hMasses = codeMassToLogMsun(reform(mt.hMass[0,*]))
  fracs = fltarr(n_elements(minMasses),numSnaps)
  
  redshifts = reverse(snapNumToRedshift(sP=sP,snap=indgen(numSnaps)+mt.minSnap))
  
  ; loop over each snapshot, calculate number above each mass threshold remaining
  for i=0,numSnaps-1 do begin
    curSnap = mt.maxSnap - i
    
    for j=0,n_elements(minMasses)-1 do begin
      w = where(mt.hMinSnap lt curSnap and hMasses ge minMasses[j],count)
      fracs[j,i] = float(count)
    endfor
    print,i,curSnap,redshifts[i],fracs[0,i]
  endfor
  
  ; convert totals to fractions
  for j=0,n_elements(minMasses)-1 do begin
    w = where(hMasses ge minMasses[j],count)
    fracs[j,*] /= float(count)
  endfor

  ; plot
  start_PS, sP.plotPath + 'halotracking.eps'
    
      xrange = [2.0,6.0]
      yrange = [0.0,max(fracs)*1.1]
      
      fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
               xtitle="Redshift",ytitle="Fraction Tracked"
      
      ; plot by mass threshold
      for j=0,n_elements(minMasses)-1 do begin
        fsc_plot,redshifts,fracs[j,*],/overplot,color=getColor(j)
      endfor
      
      legend,['M > '+string(minMasses,format='(f4.1)')],$
        textcolors=getColor(indgen(n_elements(minMasses)),/name),box=0,/bottom,/left
    
  end_PS
  
  stop

end

; plotHaloEvo(): plot evolution of halo masses, virial temperatures, positions, r200 using mergerTree()

pro plotHaloEvo, sP=sP

  minLogM = 10.0
  
  mt = mergerTreeSubset(sP=sP,/verbose)
  
  ; convert times into delta_times from starting redshift
  redshifts = 1/mt.times-1
  ages = redshiftToAgeFlat(redshifts)
  ages = (ages-ages[0])*1000.0 ;Myr
  
  ; convert masses into log
  mt.hMass = codeMassToLogMsun(mt.hMass)
  
  ; convert positions into relative to the starting position
  for i=0,n_elements(mt.hPos[0,0,*])-1 do begin
    startPos = reform(mt.hPos[0,*,i])
    
    dists_xyz = [[mt.hPos[*,0,i] - startPos[0]],$
                 [mt.hPos[*,1,i] - startPos[1]],$
                 [mt.hPos[*,2,i] - startPos[2]]]
                 
    correctPeriodicDistVecs, dists_xyz, sP=sP
   
    mt.hPos[*,*,i] = dists_xyz
  endfor
  
  ;massBins = [11.6,11.7,11.8,12.0,12.4]
  ;massBins = [10.5,10.55,10.6,10.65,10.7]
  massBins = [10.0,10.5,11.0,11.5,12.0]

  ires = 20
  respts = findgen(ires)/(ires-1) * min(ages)
  
  ; plot
  start_PS, sP.plotPath + 'haloevo_mass.eps',xs=10,ys=8
    !p.multi = [0,2,2]
    !p.charsize -= 0.4
    
    for i=0,n_elements(massBins)-2 do begin
      w = where(mt.hMass[0,*] ge massBins[i] and mt.hMass[0,*] lt massBins[i+1],count)
  
      xrange = [min(ages),0.0]
      yrange = [min(mt.hMass[*,w])*0.99,max(mt.hMass[*,w])*1.01]
      
      fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
               xtitle="Time Back [Myr]",ytitle="log ( halo mass [M"+textoidl("_{sun}")+"] )",$
               title=textoidl("z_{start} = "+string(sP.redshift,format='(f3.1)'))+" "+$
               string(massBins[i],format='(f4.1)')+" < logM < "+string(massBins[i+1],format='(f4.1)')
      
      ; smooth/interpolate testing
      for j=0,count-1 do begin
        fsc_plot,ages,smooth(mt.hMass[*,w[j]],5),line=0,color=getColor(j),/overplot
        ;res=interpol(mt.hMass[*,w[j]],ages,respts,/lsquadratic)
        ;fsc_plot,respts,res,line=0,color=getColor(j),/overplot
      endfor
      ; original data
      for j=0,count-1 do fsc_plot,ages,mt.hMass[*,w[j]],line=1,color=getColor(j),/overplot
    endfor
    
  end_PS
  
  start_PS, sP.plotPath + 'haloevo_virtemp.eps',xs=10,ys=8
    !p.multi = [0,2,2]
    !p.charsize -= 0.4
    
    for i=0,n_elements(massBins)-2 do begin
      w = where(mt.hMass[0,*] ge massBins[i] and mt.hMass[0,*] lt massBins[i+1],count)
      
      xrange = [min(ages),0.0]
      yrange = [min(mt.hVirTemp[*,w])*0.99,max(mt.hVirTemp[*,w])*1.01]
      
      fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
               xtitle="Time Back [Myr]",ytitle="halo T"+textoidl("_{vir}")+" [K]",$
               title=textoidl("z_{start} = "+string(sP.redshift,format='(f3.1)'))+" "+$
               string(massBins[i],format='(f4.1)')+" < logM < "+string(massBins[i+1],format='(f4.1)')
               
      for j=0,count-1 do fsc_plot,ages,smooth(mt.hVirTemp[*,w[j]],5),line=0,color=getColor(j),/overplot
      for j=0,count-1 do fsc_plot,ages,mt.hVirTemp[*,w[j]],line=1,color=getColor(j),/overplot
    endfor
    
  end_PS

  start_PS, sP.plotPath + 'haloevo_virrad.eps',xs=10,ys=8
    !p.multi = [0,2,2]
    !p.charsize -= 0.4
    
    for i=0,n_elements(massBins)-2 do begin
      w = where(mt.hMass[0,*] ge massBins[i] and mt.hMass[0,*] lt massBins[i+1],count)
      
      xrange = [min(ages),0.0]
      yrange = [min(mt.hVirRad[*,w])*0.94,max(mt.hVirRad[*,w])*1.06]
      
      fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
               xtitle="Time Back [Myr]",ytitle="log ( halo R"+textoidl("_{vir}")+" [kpc] )",$
               title=textoidl("z_{start} = "+string(sP.redshift,format='(f3.1)'))+" "+$
               string(massBins[i],format='(f4.1)')+" < logM < "+string(massBins[i+1],format='(f4.1)')
               
      for j=0,count-1 do fsc_plot,ages,smooth(mt.hVirRad[*,w[j]],5),line=0,color=getColor(j),/overplot
      for j=0,count-1 do fsc_plot,ages,mt.hVirRad[*,w[j]],line=1,color=getColor(j),/overplot
    endfor
    
  end_PS
  
  start_PS, sP.plotPath + 'haloevo_pos0.eps',xs=10,ys=8
    !p.multi = [0,2,2]
    !p.charsize -= 0.4
    
    for i=0,n_elements(massBins)-2 do begin
      w = where(mt.hMass[0,*] ge massBins[i] and mt.hMass[0,*] lt massBins[i+1],count)
      
      xrange = [min(ages),0.0]
      yrange = [-300,300] ;kpc offset from initial
      
      fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
               xtitle="Time Back [Myr]",ytitle="Delta Position X [kpc]",$
               title=textoidl("z_{start} = "+string(sP.redshift,format='(f3.1)'))+" "+$
               string(massBins[i],format='(f4.1)')+" < logM < "+string(massBins[i+1],format='(f4.1)')
               
      for j=0,count-1 do fsc_plot,ages,smooth(mt.hPos[*,0,w[j]],5),$
        line=0,color=getColor(j),/overplot
      for j=0,count-1 do fsc_plot,ages,mt.hPos[*,0,w[j]],$
        line=1,color=getColor(j),/overplot
    endfor
    
  end_PS
  stop
end
