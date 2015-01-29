; mergerTree.pro
; cosmological halo tracking / merger tree through time
; dnelson jan.2015

; mergerTree(): construct simplified merger tree for tracking halos/subhalos through time across snaps
; 
; makeNum : if specified, make the catalogs for the specified number of snaps back from sP.snap
;           (otherwise just make it for sP.snap and/or return existing catalog for sP.snap)

function mergerTree, sP=sP, makeNum=makeNum
  compile_opt idl2, hidden, strictarr, strictarrsubs
  if sP.snap eq 0 then message,'Error'
  
  ; config
  partMatchFracTol = 0.1   ; 55% minimum match between particle members of the specified type
  ;massDiffFracTol  = 0.1   ; 40% agreement in total mass or better
  ;positionDiffTol  = 200.0 ; 200kpc maximum separation
  minSubgroupLen   = 60    ; do not try to match subgroups with less than N total particles (for speed)

  ptNum = partTypeNum('dm') ; use dark matter particles (only) for matching

  ; set minimum/ending snapshot (maxmimum/ending redshift)
  minSnap = sP.snap
  if keyword_set(makeNum) then minSnap = sP.snap - makeNum + 1
  
  ; earliest possible snapshot to create Parents for is the one after the first group catalog
  if minSnap le sP.groupCatRange[0] then minSnap = 1 ;sP.groupCatRange[0] + 1

  ; set maximum/starting snapshot (minimum/starting redshift)
  maxSnap = sP.snap

  ; check for existing catalog
  saveFilename = sP.derivPath+'mergerTree/mergerTree.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.sav'
  
  if ~keyword_set(makeNum) then begin
    if file_test(saveFilename) then begin
      restore, saveFilename
      return, Parent
    endif
    
    ; not found, just make the one
    minSnap = sP.snap
  endif
  
  ; if creating new catalog, start loop over snapshots
  for m=maxSnap,minSnap,-1 do begin
    sP.snap = m    
    ; if at maxSnap, load current group catalog
    if sP.snap eq maxSnap then gcCur = loadGroupCat(sP=sP,/readIDs)
    
    ; or, if not at maxSnap, move "previous" group catalog to current
    if sP.snap lt maxSnap then gcCur = gcPrev
    
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
    calcMatch,partIDs_cur,partIDs_prev,cur_ind,prev_ind,count=matchCount
   
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

      ;if massDiffFrac lt massDiffFracTol and positionDiff lt positionDiffTol and $
      ;   maxFrac gt partMatchFracTol then Parent[i] = max_index
         
      ;if positionDiff lt positionDiffTol and maxFrac gt partMatchFracTol then Parent[i] = max_index
      
      ; NEW: no pos or mass tolerances, always pick max index (for zooms)
      if maxFrac gt partMatchFracTol then Parent[i] = max_index
         
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
  
  sP.snap = maxSnap ; restore sP.snap

  if ~keyword_set(makeNum) then return,Parent
end

; mergerTreeSubset():
;                     walk the merger tree from sP.snap back to start of the group catalogs and calculate
;                     the snapshot to which each halo can be tracked in time. also calculate
;                     a few properties for these subgroups back in time (r_vir,T_vir,position) and 
;                     the subset of the galaxycat at maxSnap that we can track back using this selection

function mergerTreeSubset, sP=sP, verbose=verbose
  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  maxSnap = sP.snap
  minSnap = sP.groupCatRange[0] + 1 ;z=6
  
  smoothKer = 3 ; 1,3,5 number of snapshots
  units = getUnits()
  
  ; set saveFilename and check for existence
  saveFilename = sP.derivPath + 'mTreeSubset.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(maxSnap)+'-'+str(minSnap)+'.sav'  

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
      gcIDs = gcIDList(gc=gc,select='all')

      gcIDs_pri = gcIDList(gc=gc,select='pri')
      
      ; mark secondary subgroups as invalid at the start
      gcIDs[ gcIDList(gc=gc,select='sec') ] = -1
      
      times     = fltarr(numBack+floor(smoothKer/2.0))
      snaps     = lonarr(numBack+floor(smoothKer/2.0))
      hPos      = fltarr(numBack+floor(smoothKer/2.0),3,n_elements(gcIDs))
      hVel      = fltarr(numBack+floor(smoothKer/2.0),3,n_elements(gcIDs))
      hMass     = fltarr(numBack+floor(smoothKer/2.0),n_elements(gcIDs))
      hVirRad   = fltarr(numBack+floor(smoothKer/2.0),n_elements(gcIDs))
      hVirTemp  = fltarr(numBack+floor(smoothKer/2.0),n_elements(gcIDs))
      hInd      = lonarr(numBack+floor(smoothKer/2.0),n_elements(gcIDs)) - 1
      noParMask = bytarr(n_elements(gcIDs))
      hMinSnap  = intarr(n_elements(gcIDs)) - 1
    endif
  
    ; which gcIDs still valid? only new "no parents"
    wNoPar = where(gcIDs eq -1 and noParMask eq 0B,wNoParCount)
    
    if wNoParCount gt 0 then begin
      noParMask[wNoPar] = 1B
      hMinSnap[wNoPar] = sP.snap + 1 ; mark end of tracking for halos with no valid parent at this snapshot
    endif
    
    wPar = where(hMinSnap[gcIDs_pri] eq -1,nleft)
    wPar10 = where(hMinSnap[gcIDs_pri] eq -1 and hMass[0,gcIDs_pri] ge 1.0,nleft10)
    tot10 = where(hMass[0,gcIDs_pri] ge 1.0,ntot10)
    
    if keyword_set(verbose) then $
      print,' ['+string(sP.snap,format='(i3)')+'] remaining: '+string(nleft,format='(i5)')+' ('+$
        string(float(nleft)/n_elements(gcIDs)*100,format='(f4.1)')+'%) above 10^10: '+$
        string(nleft10,format='(i5)')+' ('+string(float(nleft10)/ntot10*100,format='(f4.1)')+'%) '
        
    ; store only currently valid gcIDs
    w = where(gcIDs ne -1,count)
    if count eq 0 then message,'error'
    
    times[i] = h.time
    hPos[i,*,w]   = subgroupCen[*,gcIDs[w]]
    hVel[i,*,w]   = gc.subgroupVel[*,gcIDs[w]]
    hMass[i,w]    = gc.subgroupMass[gcIDs[w]]
    hVirRad[i,w]  = gc.group_r_crit200[gc.subgroupGrNr[gcIDs[w]]]
    hVirTemp[i,w] = alog10(codeMassToVirTemp(gc.subgroupMass[gcIDs[w]],sP=sP))
    hInd[i,w]     = gcIDs[w]
    snaps[i]      = sP.snap
    
    ; load mergerTree and change subgroup IDs to parents at the prior snapshot
    Parent = mergerTree(sP=sP)

    gcIDs[w] = Parent[gcIDs[w]] ; change to parent IDs

    sP.snap -= 1
  endfor
  
  Parent = !NULL
  gcIDs  = !NULL
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
  
  ; replicated galcat parent IDs for each galcatSub member
  gcIndOrig = galCatRepParentIDs(galcat=galcat)
  
  r = {gcIndOrig:gcIndOrig,gcIndOrigTrMC:-1,gcIndOrigTrVel:-1,$
       hPos:hPos,hVel:hVel,hPosSm:hPosSm,hVelSm:hVelSm,$
       hVirRad:hVirRad,hVirTemp:hVirTemp,hMass:hMass,hMinSnap:hMinSnap,$
       hInd:hInd,snaps:snaps,times:times,minSnap:minSnap,maxSnap:maxSnap,smoothKer:smoothKer}

  if sP.trMCPerCell gt 0 then begin
    gcIndOrigTrMC = galCatRepParentIDs(galcat=galcat,child_counts=galcat.trMC_cc)
    r = mod_struct( r, 'gcIndOrigTrMC', gcIndOrigTrMC )
  endif
  
  if sP.trVelPerCell gt 0 then begin
    if n_elements(galcat.trVel_cc) gt 1 then begin ; is just -1 if not included
      gcIndOrigTrVel = galCatRepParentIDs(galcat=galcat,child_counts=galcat.trVel_cc)
      r = mod_struct( r, 'gcIndOrigTrVel', gcIndOrigTrVel )
    endif
  endif
       
  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
  
  return,r
end

; trackHaloPosition(): return a new position for some halo gcID at an earlier snapshot using mergerTree
    
function trackHaloPosition, sP=sP, gcID=gcID, endSnap=endSnap
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
    step = -1 ; back in snapshots, back in time
    if endSnap gt sP.snap then step = 1 ; forward in snapshots
    
    ; try to trace back to this snapshot for new halo center position
    failSnap = -1
    sgCenters = fltarr( 3,abs( sP.snap-endSnap-step) ) ; xyz ckpc
    gcIDs     = lonarr( abs(sP.snap-endSnap-step) )
    snaps     = intarr( abs(sP.snap-endSnap-step) )
    times     = fltarr( abs(sP.snap-endSnap-step) ) ; Gyr age
    
    gcIDcur = gcID ; starting halo ID
    startSnap = sP.snap ; starting snapshot number

    for m=startSnap,endSnap,step do begin
      ind = abs(startSnap-m) ; index to store this snapshot result in
      
      gcIDs[ind] = gcIDcur
      print,m,ind,gcIDcur,failSnap
      
      ; record times for possible extrapolation
      h = loadSnapshotHeader(sP=sP)
      times[ind] = redshiftToAgeFlat(1/h.time-1)
      snaps[ind] = m
      
      sP.snap = m + 1 * (step eq 1) ; want to load one snapshot ahead if moving forward
      
      ; haven't failed yet, keep recording positions
      if failSnap eq -1 and sP.snap gt 0 then begin
        sgcen = subgroupPosByMostBoundID(sP=sP)
        sgCenters[*,ind] = sgcen[*,gcIDcur]
              
        ; load most massive progenitor
        Parent = mergerTree(sP=sP)
        
        ; if going backwards, directly move index
        if step eq -1 then begin
          gcIDcur = Parent[gcIDcur]
          
          if gcIDcur eq -1 then failSnap = sP.snap
        endif
        
        ; if going forwards, find which child has this parent
        if step eq 1 then begin
          gcIDcur = ( where(Parent eq gcIDcur,count) )[0]
          if count eq 0 then failSnap = sP.snap
        endif
                
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
        sP.snap = endSnap
        h = loadSnapshotHeader(sP=sP)
        timeEnd = redshiftToAgeFlat(1/h.time-1)
        
        inds = where(times ne 0.0,count)
        if count lt 3 then message,'error'
        
        xPos = interpol(sgCenters[0,inds],times[inds],timeEnd)
        yPos = interpol(sgCenters[1,inds],times[inds],timeEnd)
        zPos = interpol(sgCenters[2,inds],times[inds],timeEnd)
        haloPos = [xPos,yPos,zPos]
      endelse
    endelse
    
    sP.snap = startSnap
    return,{sgCenters:sgCenters,haloPosEnd:haloPos,gcIDs:gcIDs,snaps:snaps,times:times}
end

; plotHaloTracking(): plot fraction of successful halo adaptive tracking vs redshift

pro plotHaloTracking, sP=sP
  compile_opt idl2, hidden, strictarr, strictarrsubs
  minMasses = [10.0,10.5,11.0,11.5] ; z=2 masses in log(msun)
  
  mt = mergerTreeSubset(sP=sP,/verbose)
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

; plotHaloEvoSingle(): plot evolution of a single (zoom) halo properties

pro plotHaloEvoSingle
  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function zoomTargetHalo

  ; config
  sP = simParams(run='zoom_20Mpc',res=9,hInd=0,redshift=2.0)
  colors = ['red','blue','orange']
  
  ; load
  hInd = zoomTargetHalo(sP=sP)
  mt   = mergerTreeSubset(sP=sP,/verbose)
  
  ; calculate
  redshifts = 1/mt.times-1
  ages = redshiftToAgeFlat(redshifts)
  ages = (ages-ages[0])*1000 ; Myr since sP.redshift
  mt.hMass = codeMassToLogMsun(mt.hMass)
  
  print,'Target halo: ['+str(hInd)+'] tracked to minSnap: ['+str(mt.minSnap[hInd])+'] (z='+$
    string(min(redshifts),format='(f4.2)')+' to z='+string(max(redshifts),format='(f4.2)')+')'
  
  startPos = reform(mt.hPos[0,*,hInd])
  
  xyz_delta    = fltarr(3,n_elements(ages))
  xyz_delta_sm = fltarr(3,n_elements(ages))
  
  for i=0,2 do begin
    xyz_delta[i,*]    = reform(mt.hPos[*,i,hInd]) - startPos[i]
    xyz_delta_sm[i,*] = reform(mt.hPosSm[*,i,hInd]) - startPos[i]
  endfor
  
  ; plot
  pos = plot_pos(rows=2,cols=2,/gap)
  start_PS,sP.plotPath + 'haloEvoSingle_' + sP.saveTag + '.eps', /big
    cgPlot,[0],[0],/nodata,xrange=minmax(ages),yrange=minmax(xyz_delta),$
      xtitle="Delta Time Back [Myr]",ytitle="Delta Position [ckpc]",pos=pos[0]
      
    for i=0,2 do begin
      cgPlot,ages,xyz_delta[i,*],psym=-4,line=0,color=cgColor(colors[i]),/overplot
      ;cgPlot,ages,xyz_delta_sm[i,*],psym='filled circle',color=cgColor(colors[i]),/overplot
    endfor
    
    legend,['x','y','z'],textcolors=colors,/bottom,/right
    
    cgPlot,[0],[0],/nodata,xrange=minmax(ages),yrange=minmax(mt.hMass[*,hInd]),$
      xtitle="Delta Time Back [Myr]",ytitle="Halo Mass [log Msun]",pos=pos[1],/noerase
      
    cgPlot,ages,mt.hMass[*,hInd],psym=-4,color=cgColor('black'),/overplot
    ;cgPlot,ages,1.95*mt.hVirTemp[*,hInd],psym=-4,color=cgColor('forest green'),/overplot
    ;legend,['mass','2*virTemp'],textcolors=['black','forest green'],/top,/left
    
    cgPlot,[0],[0],/nodata,xrange=minmax(ages),yrange=minmax(mt.hVirRad[*,hInd]),$
      xtitle="Delta Time Back [Myr]",ytitle="Halo VirRad [ckpc]",pos=pos[2],/noerase
      
    cgPlot,ages,mt.hVirRad[*,hInd],psym=-4,color=cgColor('black'),/overplot
    
    cgPlot,[0],[0],/nodata,xrange=minmax(ages),yrange=minmax(mt.hVel[*,*,hInd]),$
      xtitle="Delta Time Back [Myr]",ytitle="Halo Velocity [km/s]",pos=pos[3],/noerase
      
    cgPlot,minmax(ages),[0,0],line=1,color=cgColor('gray'),/overplot
    for i=0,2 do $
      cgPlot,ages,mt.hVel[*,i,hInd],psym=-4,line=0,color=cgColor(colors[i]),/overplot
    legend,['x','y','z'],textcolors=colors,/bottom,/right
  end_PS
  
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
  
  massBins = [11.6,11.7,11.8,12.0,12.4]
  ;massBins = [10.5,10.55,10.6,10.65,10.7]
  ;massBins = [10.0,10.5,11.0,11.5,12.0]

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
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
             xtitle="Time Back [Myr]",ytitle="log ( halo mass [M"+textoidl("_{sun}")+"] )",$
             title=textoidl("z_{start} = "+string(sP.redshift,format='(f3.1)'))+" "+$
             string(massBins[i],format='(f4.1)')+" < logM < "+string(massBins[i+1],format='(f4.1)')
      
      ; smooth/interpolate testing
      for j=0,count-1 do begin
        cgPlot,ages,smooth(mt.hMass[*,w[j]],5),line=0,color=getColor(j),/overplot
        ;res=interpol(mt.hMass[*,w[j]],ages,respts,/lsquadratic)
        ;cgPlot,respts,res,line=0,color=getColor(j),/overplot
      endfor
      ; original data
      for j=0,count-1 do cgPlot,ages,mt.hMass[*,w[j]],line=1,color=getColor(j),/overplot
    endfor
    
  end_PS
  
  start_PS, sP.plotPath + 'haloevo_virtemp.eps',xs=10,ys=8
    !p.multi = [0,2,2]
    !p.charsize -= 0.4
    
    for i=0,n_elements(massBins)-2 do begin
      w = where(mt.hMass[0,*] ge massBins[i] and mt.hMass[0,*] lt massBins[i+1],count)
      
      xrange = [min(ages),0.0]
      yrange = [min(mt.hVirTemp[*,w])*0.99,max(mt.hVirTemp[*,w])*1.01]
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
             xtitle="Time Back [Myr]",ytitle="halo T"+textoidl("_{vir}")+" [K]",$
             title=textoidl("z_{start} = "+string(sP.redshift,format='(f3.1)'))+" "+$
             string(massBins[i],format='(f4.1)')+" < logM < "+string(massBins[i+1],format='(f4.1)')
               
      for j=0,count-1 do cgPlot,ages,smooth(mt.hVirTemp[*,w[j]],5),line=0,color=getColor(j),/overplot
      for j=0,count-1 do cgPlot,ages,mt.hVirTemp[*,w[j]],line=1,color=getColor(j),/overplot
    endfor
    
  end_PS

  start_PS, sP.plotPath + 'haloevo_virrad.eps',xs=10,ys=8
    !p.multi = [0,2,2]
    !p.charsize -= 0.4
    
    for i=0,n_elements(massBins)-2 do begin
      w = where(mt.hMass[0,*] ge massBins[i] and mt.hMass[0,*] lt massBins[i+1],count)
      
      xrange = [min(ages),0.0]
      yrange = [min(mt.hVirRad[*,w])*0.94,max(mt.hVirRad[*,w])*1.06]
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
             xtitle="Time Back [Myr]",ytitle="log ( halo R"+textoidl("_{vir}")+" [kpc] )",$
             title=textoidl("z_{start} = "+string(sP.redshift,format='(f3.1)'))+" "+$
             string(massBins[i],format='(f4.1)')+" < logM < "+string(massBins[i+1],format='(f4.1)')
               
      for j=0,count-1 do cgPlot,ages,smooth(mt.hVirRad[*,w[j]],5),line=0,color=getColor(j),/overplot
      for j=0,count-1 do cgPlot,ages,mt.hVirRad[*,w[j]],line=1,color=getColor(j),/overplot
    endfor
    
  end_PS
  
  start_PS, sP.plotPath + 'haloevo_pos0.eps',xs=10,ys=8
    !p.multi = [0,2,2]
    !p.charsize -= 0.4
    
    for i=0,n_elements(massBins)-2 do begin
      w = where(mt.hMass[0,*] ge massBins[i] and mt.hMass[0,*] lt massBins[i+1],count)
      
      xrange = [min(ages),0.0]
      yrange = [-300,300] ;kpc offset from initial
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
             xtitle="Time Back [Myr]",ytitle="Delta Position X [kpc]",$
             title=textoidl("z_{start} = "+string(sP.redshift,format='(f3.1)'))+" "+$
             string(massBins[i],format='(f4.1)')+" < logM < "+string(massBins[i+1],format='(f4.1)')
               
      for j=0,count-1 do cgPlot,ages,smooth(mt.hPos[*,0,w[j]],5),$
        line=0,color=getColor(j),/overplot
      for j=0,count-1 do cgPlot,ages,mt.hPos[*,0,w[j]],$
        line=1,color=getColor(j),/overplot
    endfor
    
  end_PS
  stop
end

; plotTrackedDists(): plot mass/tvir at z=2 and at time when each halo stopped tracking

pro plotTrackedDists

  ; config
  sP = simParams(res=512,run='tracer',redshift=2.0)
  mt = mergerTreeSubset(sP=sP)
  binsize = 0.15
  
  ; "IGM temperature" model (constant)
  T_IGM = 4e4 ; K
  massRanges = list([10.0,12.5],[10.0,10.75],[10.75,11.5],[11.5,12.5]) ; log Msun
  
  start_PS, sP.plotPath + 'td.endz.'+sP.plotPrefix+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps',/big
  
    !p.multi = [0,2,2]
    !p.charsize = 1.1
  
    foreach massRange,massRanges,k do begin
    
      ; choose halos
      w = where(mt.hMinSnap ne -1 and codeMassToLogMsun(mt.hMass[0,*]) ge massRange[0] and $
                                      codeMassToLogMsun(mt.hMass[0,*]) lt massRange[1],count)
      ms_ind = mt.maxSnap - (mt.hMinSnap[w])
  
      ; pick ending snapshots based on TIGM
      temp_ratio = 10.0^mt.hVirTemp[*,w] / T_IGM
  
      end1 = lonarr(count) ; snapshot number
      end2 = lonarr(count)
  
      for i=0,count-1 do begin
        loc_temp_ratio = temp_ratio[*,i]
        ww = where(loc_temp_ratio ge 1.0,count_loc)
        if count_loc eq 0 then print,'error',i
        end1[i] = mt.maxSnap - max(ww)
    
        ww = where(loc_temp_ratio ge 2.0,count_loc)
        if count_loc eq 0 then print,'error',i
        end2[i] = mt.maxSnap - max(ww)   
      endfor    
    
      ; plot
      title = "masses "+string(massRange[0],format='(f4.1)')+" - "+string(massRange[1],format='(f4.1)')
      if k eq 0 then title = "all masses"
      
      cgPlot,[0],[0],/nodata,xtitle="Ending Redshift",ytitle="Fraction Lost",$
        yrange=[0.0,1.1],/ys,xrange=[6,2],/xs,title=title
    
      hist = histogram(snapNumToRedshift(snap=mt.hMinSnap[w],sP=sP),binsize=2*binsize,loc=loc)
      cgPlot,loc+binsize,total(hist,/cum)/float(count),color=cgColor('blue'),/overplot
    
      hist = histogram(snapNumToRedshift(snap=end1,sP=sP),binsize=2*binsize,loc=loc)
      cgPlot,loc+binsize,total(hist,/cum)/float(count),color=cgColor('red'),/overplot    
    
      hist = histogram(snapNumToRedshift(snap=end2,sP=sP),binsize=2*binsize,loc=loc)
      cgPlot,loc+binsize,total(hist,/cum)/float(count),color=cgColor('green'),/overplot        
    
      legend,['end track','end Tvir>TIGM','end Tvir>2 TIGM'],textcolors=['blue','red','green'],$
        box=0,/bottom,/left
    endforeach
  end_PS
  
  stop

end
