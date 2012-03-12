; cosmoUtil.pro
; cosmological simulations - utility functions
; dnelson mar.2012

; redshiftToSnapNum(): convert redshift to the nearest snapshot number

function redshiftToSnapNum, redshiftList, sP=sP, verbose=verbose

  if not keyword_set(verbose) then verbose = 0
  
  saveFileName = sP.derivPath + sP.savPrefix + '_snapnum.redshift.sav'

  if not (file_test(saveFileName)) then begin

    nSnaps = 400 ; maximum
    
    redshifts = fltarr(nSnaps) - 1
    times     = fltarr(nSnaps) - 1
    
    for m=0,nSnaps-1 do begin
      ; format filename
      ext = str(string(m,format='(I3.3)'))
      f = sP.simPath + 'snapdir_' + ext + '/snap_' + ext + '.0.hdf5'
      
      ; single file per group catalog
      if (not file_test(f)) then $
        f = sP.simPath + 'snap_' + ext + '.hdf5'
        
      ; single groupordered file per group catalog
      if (not file_test(f)) then $
        f = sP.simPath + 'snap-groupordered_' + ext + '.hdf5'
      
      ; if file doesn't exist yet, skip (e.g. every other file deleted)
      if (not file_test(f)) then continue
    
      ; load hdf5 header and save time+redshift
      fileID   = h5f_open(f)
      headerID = h5g_open(fileID,"Header")
      
      redshifts[m] = h5a_read(h5a_open_name(headerID,"Redshift"))
      times[m]     = h5a_read(h5a_open_name(headerID,"Time"))
      
      h5g_close, headerID
      h5f_close, fileID
    endfor
  
    ; save/restore
    save,nSnaps,redshifts,times,filename=saveFileName
  endif else begin
    restore,saveFileName
  endelse

  ; for -1 return !NULL
  if (total(redshiftList) eq -1) then return,!NULL
  
  ; return array
  snapNum = intarr(n_elements(redshiftList))
  
  foreach redshift,redshiftList,i do begin
    if (redshift eq 0.0) then begin
      ; for redshift zero hard select #314 (ComparisonProject)
      snapNum[i] = sP.snapRange[1]
    endif else begin
      w = where(abs(redshifts - redshift) eq min(abs(redshifts - redshift)),count)
      if (count eq 2) then w = w[0]
      if (count eq 0) then return,-1
      snapNum[i] = w[0]
    endelse
  
    if (verbose) then $
      print,'Found nearest snapshot to z = ' + str(redshift) + ' (num = ' + str(snapNum) + ')'
  endforeach

  if (n_elements(snapNum) eq 1) then snapNum = snapNum[0]

  return,snapNum
end

; sgIDList(): return a sub-list of subgroup IDs from a given group catalog sg
; 
; select: one of the following
;   'pri' : members of the first subgroup of each group only ("background"/"main subhalo"/"halo")
;   'sec' : members of the non-first subgroups of each group only ("satellites"/"subhalos")
;   'all' : all subgroups
; 
; minNumPart=1 : require a minimum number of particles in (subfind) group to include

function sgIDList, sP=sP, sg=sg, select=select, minNumPart=minNumPart

  forward_function loadSubhaloGroups

  if (not keyword_set(select)) then begin
    print,'Error: Must specify select.'
    stop
  endif
    
  ; load galaxy cat if necessary
  if not keyword_set(sg) then begin
    if not keyword_set(sP) then begin
      print,'Error: Must specify sg or sP.'
      stop
    endif
    sg = loadSubhaloGroups(sP.simPath,sP.snap)
  endif

  minNumPartVal = 100 ; in subgroup (total, dm+gas+stars)
  ; TODO: reconcile with the other minNumGasPart check on gas only

  prevGrNr   = -1
  valSGids   = []

  if (select eq 'pri') then begin
  
    ; "background"/"main subhalos"/"halos" only
    for i=0,n_elements(sg.subgroupLen)-1 do begin
      if (sg.subgroupGrnr[i] eq prevGrNr) then begin
        prevGrNr = sg.subgroupGrnr[i]
      endif else begin
      
        if (keyword_set(minNumPart)) then begin
          if (sg.subgroupLen[i] ge minNumPartVal) then $
            valSGids = [valSGids,i]
        endif else begin
          valSGids = [valSGids,i]
        endelse ;minNumPart
        
        prevGrNr = sg.subgroupGrnr[i]
      endelse
    endfor
    
  endif
  
  if (select eq 'sec') then begin
  
    ; "satellites"/"subhalos" only
    for i=0,n_elements(sg.subgroupLen)-1 do begin
      if (sg.subgroupGrnr[i] ne prevGrNr) then begin
        prevGrNr = sg.subgroupGrnr[i]
      endif else begin
      
        if (keyword_set(minNumPart)) then begin
          if (sg.subgroupLen[i] ge minNumPartVal) then $
            valSGids = [valSGids,i]
        endif else begin
          valSGids = [valSGids,i]
        endelse ;minNumPart
        
      endelse
    endfor
    
  endif
  
  if (select eq 'all') then begin
  
    ; both primary and secondary
    for i=0,n_elements(sg.subgroupLen)-1 do begin
    
      if (keyword_set(minNumPart)) then begin
        if (sg.subgroupLen[i] ge minNumPartVal) then $
          valSGids = [valSGids,i]
      endif else begin
        valSGids = [valSGids,i]
      endelse ;minNumPart
      
    endfor
  
  endif
  
  return, valSGids

end

; sgPriChildInd(): get the subgroup index of the primary subgroup of a given group, or -1 if none

function sgPriChildInd, gc=gc, haloID=haloID

  if (n_elements(gc) eq 0 or n_elements(haloID) eq 0) then stop

  sgInd = min(where(gc.subgroupGrNr eq haloID,count))
  fsInd = gc.groupFirstSub[haloID]
  
  ; skip if group has no subgroup to obtain this value from
  if (count eq 0) then return,-1
  
  ; double-check (failing either should indicate a problem in the group catalogs, right?)
  if (sgInd ne fsInd) then return,-1 ;stop
  if (gc.subgroupGrNr[sgInd] ne haloID) then return,-1 ;stop

  return,sgInd

end

; sgPIDList(): return a list of member particle IDs from a given group catalog sg
;
; select: one of the following
;   'pri' : members of the first subgroup of each group only ("background"/"main subhalo"/"halo")
;   'sec' : members of the non-first subgroups of each group only ("satellites"/"subhalos")
;   'all' : all subgroups
;
; gasOnly=1 : load gas ids from snapshot and match to restrict return to only gas particle ids
; dmOnly=1  : load dm ids from snapshot and match to restrict return to only dark matter ids

function sgPIDList, sg=sg, select=select, gasOnly=gasOnly, dmOnly=dmOnly

  if (not keyword_set(sg) or not keyword_set(select)) then begin
    print,'Error: sgPIDList: Bad inputs.'
    return,0
  endif

  ; get list of appropriate group ids
  valSGids = sgIDList(sg=sg,select=select)
  
  ; make list of particle ids in these groups
  subgroupPIDs = []
  
  foreach sgID, valSGids do begin
    ; select subgroup
    sgIDs = sg.subGroupIDs[sg.subGroupOffset[sgID] : sg.subGroupOffset[sgID] + $
                           sg.subGroupLen[sgID] - 1]
    
    subgroupPIDs = [subgroupPIDs, sgIDs]
  endforeach
    
  ; particle type restriction
  if (keyword_set(gasOnly) or keyword_set(dmOnly)) then begin
    if (keyword_set(gasOnly)) then partName='gas'
    if (keyword_set(dmOnly))  then partName='dm'
    
    ; load gas/dm ids from targetSnap and restrict
    ;pIDs = loadSnapshotSubset(gadgetPath,snapNum=snap,partType=partName,field='ids')
    
    ;match, pIDs, sgEnd_subgroupIDs, pIDs_ind, sgIDs_ind, count=count_pID
    
    ;pIDs = !NULL
    ;pIDs_ind = !NULL
    
    ; keep only sgIDs[sgIDs_ind]
    ;if (count_pID gt 0) then $
    ;    sgEnd_subgroupIDs = sgEnd_subgroupIDs[sgIDs_ind]
    print,'ERROR need to fix this part'
    exit
  endif

  return, subgroupPIDs
end

; galCatRepParentIDs(): for the galaxy catalog, replicate the list of ordered parent IDs such that
;                       the return array is the same size as the number of gas particle ids with
;                       each element the id of its parent
;
; sgIDListPri=1 : should be a list of primary subhalo IDs, e.g. sgIDListPri=sgIDList(sg=sg,/pri)
;                 if set, for secondary groups, return instead of the ID of the primary parent
;                 (for primary groups the return is unchanged)

function galCatRepParentIDs, gc=gc, sgIDListPri=sgIDListPri

    ; arrays to hold parent IDs for each gas particle
    sg_ind_gal  = lonarr(n_elements(gc.galaxyIDs))
    sg_ind_gmem = lonarr(n_elements(gc.groupmemIDs))
    
    if not keyword_set(sgIDListPri) then $
      sgIDListPri = indgen(n_elements(gc.galaxyLen)) ; valid id list set to all
    
    for sgID=0,n_elements(gc.galaxyLen)-1 do begin
    
        ; parent ID (pri or sec)
        if (total(sgID eq sgIDListPri) eq 0) then begin
          dists = sgID - sgIDListPri
          dists = dists[where(dists gt 0)]
          w = where(dists eq min(dists),count)
          parID = sgIDListPri[w[0]]
        endif else begin
          parID = sgID
        endelse
        
        ; galaxies
        groupStart = gc.galaxyOff[sgID]
        groupEnd   = groupStart + gc.galaxyLen[sgID]
        
        if (gc.galaxyLen[sgID] gt 0) then $
          sg_ind_gal[groupStart:groupEnd-1] = cmreplicate(parID,groupEnd-groupStart)
          
        ; group members
        groupStart = gc.groupmemOff[sgID]
        groupEnd   = groupStart + gc.groupmemLen[sgID]
        
        if (gc.groupmemLen[sgID] gt 0) then $
          sg_ind_gmem[groupStart:groupEnd-1] = cmreplicate(parID,groupEnd-groupStart)
    endfor
    
    r = {gal:sg_ind_gal,gmem:sg_ind_gmem}
    return,r
end

; galCatParentProperties: calculate some property of the parent galaxy/group for every gas particle
;                         in the galaxy catalog at some snapshot
; virTemp=1 : virial temperature
; mass=1    : total mass (from catalog, dm+baryon)
; rVir=1    : virial radius (r_200 critical)

function galCatParentProperties, sP=sP, virTemp=virTemp, mass=mass, rVir=rVir

  forward_function galaxyCat, snapNumToRedshift, codeMassToLogMsun

  ; load group catalog for masses
  sg = loadSubhaloGroups(sP.simPath,sP.snap)

  ; load galaxy catalog
  gc = galaxyCat(res=sP.res,run=sP.run,snap=sP.snap)

  ; replicate parent IDs
  sgInd = galCatRepParentIDs(gc=gc)
  
  ; arrays
  gal  = fltarr(n_elements(gc.galaxyIDs))
  gmem = fltarr(n_elements(gc.groupmemIDs))
  
  ; masses (log msun)
  if keyword_set(mass) then begin
    gal  = sg.subgroupMass[sgInd.gal]
    gmem = sg.subgroupMass[sgInd.gmem]
    
    gal  = codeMassToLogMsun(gal)
    gmem = codeMassToLogMsun(gmem)
  endif

  ; calculate virial temperatures (K)
  if keyword_set(virTemp) then begin
    gal  = sg.subgroupMass[sgInd.gal]
    gmem = sg.subgroupMass[sgInd.gmem]
    
    redshift = snapNumToRedshift(snap=sP.snap)
  
    gal  = codeMassToVirTemp(gal,redshift)
    gmem = codeMassToVirTemp(gmem,redshift)
  endif
  
  if keyword_set(rVir) then begin
    gal  = sg.subgroupGrnr[sgInd.gal]
    gmem = sg.subgroupGrnr[sgInd.gmem]
    
    gal  = sg.group_r_crit200[gal]
    gmem = sg.group_r_crit200[gmem]
  endif

  r = {gal:gal,gmem:gmem}
  return,r
end

; gcINDList(): return a list of member particle indices from a given list of subgroup IDs

function gcINDList, sP=sP, gc=gc, sgIDList=sgIDList

  forward_function galaxyCat

  ; load galaxy cat if necessary
  if not keyword_set(gc) then begin
    if not keyword_set(sP) then begin
      print,'Error: Must specific gc or sP.'
      return,0
    endif
    gc = galaxyCat(res=sP.res,run=sP.run,snap=sP.snap)
  endif

  galaxyInds = []
  groupmemInds = []
  
  foreach sgID, sgIDList do begin
    ; galaxy
    if (gc.galaxyLen[sgID] gt 0) then begin
      galInds  = lindgen(gc.galaxyLen[sgID]) + gc.galaxyOff[sgID] ; int overflow on indgen
      galaxyInds   = [galaxyInds, galInds]
    endif
    
    ; group member
    if (gc.groupmemLen[sgID] gt 0) then begin
      gmemInds = lindgen(gc.groupmemLen[sgID]) + gc.groupmemOff[sgID]
      groupmemInds = [groupmemInds, gmemInds]
    endif    
    
  endforeach
  
  r = {gal:galaxyInds,gmem:groupmemInds}
  return,r
  
end

; subgroupPosByMostBoundID(): compute a "best" center position in space for all subfind groups by 
;                             using the position of the most bound particles, whose IDs are stored in 
;                             the group catalog but without knowing the particle types we have to 
;                             load all gas+dm+stars particle positions

function subgroupPosByMostBoundID, sP=sP

  forward_function loadSnapshotSubset, loadGroupCat
  
  ; save/restore
  saveFilename = sP.derivPath + sP.savPrefix + str(sP.res) + '.gCen.mbID.'+str(sP.snap)+'.sav'
                 
  if file_test(saveFilename) then begin
    restore,saveFilename
  endif else begin 
  
    sg = loadGroupCat(sP.simPath,sP.snap)
  
    groupCen = fltarr(3,sg.nSubgroupsTot)
    
    ; load gas ids and pos, find matches
    ids = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='ids')
    match,sg.subgroupIdMostBound,ids,sg_ind,ids_ind,count=count1,/sort
    
    ids_ind = ids_ind[sort(sg_ind)]
    sg_ind  = sg_ind[sort(sg_ind)]
    
    pos = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='pos')
    
    groupCen[*,sg_ind] = pos[*,ids_ind]
    
    ; load stars ids and pos, find matches
    ids = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='stars',field='ids')
    
    ; skip stars for non-radiative runs or early times
    if (n_elements(ids) ne 1) then begin
      match,sg.subgroupIdMostBound,ids,sg_ind,ids_ind,count=count2,/sort
      
      ids_ind = ids_ind[sort(sg_ind)]
      sg_ind  = sg_ind[sort(sg_ind)]
      
      pos = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='stars',field='pos')
      
      groupCen[*,sg_ind] = pos[*,ids_ind]
    endif else begin
      count2 = 0
    endelse
    
    ; load dm ids and pos, find matches
    ids = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='dm',field='ids')
    match,sg.subgroupIdMostBound,ids,sg_ind,ids_ind,count=count3,/sort
    
    ids_ind = ids_ind[sort(sg_ind)]
    sg_ind  = sg_ind[sort(sg_ind)]
    
    pos = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='dm',field='pos')
    
    groupCen[*,sg_ind] = pos[*,ids_ind]
  
    if ((count1 + count2 + count3) ne sg.nSubgroupsTot) then begin
      print,'ERROR'
      stop
    endif
  
    ; save
    save,groupCen,filename=saveFilename
  endelse
  
  return, groupCen
end

; correctPeriodicDistVecs(): enforce periodic B.C. for distance vecotrs (effectively component by 
;                            component), input vecs in format fltarr[3,n]

pro correctPeriodicDistVecs, vecs, sP=sP

  w = where(vecs gt sP.boxSize*0.5,count)
  if (count ne 0) then $
    vecs[w] = vecs[w] - sP.boxSize
    
  w = where(vecs lt -sP.boxSize*0.5,count)
  if (count ne 0) then $
    vecs[w] = sP.boxSize + vecs[w]

end

; periodicDists(): calculate distances correctly taking into account periodic B.C.
; 
; if pt is one point: distance from pt to all vecs
; if pt is several points: distance from each pt to each vec (must have same number of points)

function periodicDists, pt, vecs, sP=sP

  nDimsPt = (size(pt))[0]

  if ( ((size(vecs))[0] ne 1 and (size(vecs))[0] ne 2) or $
       (size(vecs))[1] ne 3) then stop ; vecs not in expected shape
  if (nDimsPt ne 1 and nDimsPt ne 2) then stop ; something is wrong

  ; distances from one point to a vector of other points [3,n]
  if (nDimsPt eq 1) then begin
    xDist = vecs[0,*]-pt[0]
    yDist = vecs[1,*]-pt[1]
    zDist = vecs[2,*]-pt[2]
  endif
  
  ; distances from a vector of points [3,n] to another vector of other points [3,n]
  if (nDimsPt eq 2) then begin
    xDist = vecs[0,*]-pt[0,*]
    yDist = vecs[1,*]-pt[1,*]
    zDist = vecs[2,*]-pt[2,*]
  endif
  
  correctPeriodicDistVecs, xDist, sP=sP
  correctPeriodicDistVecs, yDist, sP=sP
  correctPeriodicDistVecs, zDist, sP=sP
  
  dists = reform( sqrt( xDist*xDist + yDist*yDist + zDist*zDist ) )
  
  return, dists

end

; periodicPairwiseDists(): calculate pairwise distances between all 3D points, correctly
;                          taking into account periodic B.C.

function periodicPairwiseDists, pts, sP=sP

  ; construct indices for vector computation of pairwise distances
  dims = size(pts,/dimensions)
  
  num = dims[1]*(dims[1]-1)/2
  
  ii = 0L
  index0 = lindgen(dims[1]-1) + 1
  index1 = lonarr(num,/nozero)
  index2 = lonarr(num,/nozero)
  
  for i=0,dims[1]-2 do begin
    n1 = dims[1] - (i+1)
    index1[ii:ii+n1-1] = i
    index2[ii] = index0[0:n1-1] + i
    ii += n1
  endfor
  
  ; component wise difference
  xDist = pts[0,index1] - pts[0,index2]
  yDist = pts[1,index1] - pts[1,index2]
  zDist = pts[2,index1] - pts[2,index2]
  
  ; account for periodic distance function
  correctPeriodicDistVecs, xDist, sP=sP
  correctPeriodicDistVecs, yDist, sP=sP
  correctPeriodicDistVecs, zDist, sP=sP
  
  dists = reform( sqrt( xDist*xDist + yDist*yDist + zDist*zDist ) )

  return,dists

end

; groupCenterPosByIterativeCM(): compute a "better" center position in space for all FoF groups by
;                                iteratively searching for the center of mass

function groupCenterPosByIterativeCM, sP=sP, gc=gc, haloIDs=haloIDs

  forward_function loadSnapshotSubset, loadSnapshotHeader

  if (n_elements(haloIDs) eq 0) then stop

  ; config
  radStartStop = [1.0,0.2] ; fractions of r_vir
  radStep      = 0.01 ; decrease search radius by this fraction each step

  nSteps = fix((radStartStop[0]-radStartStop[1]) / radStep + 1)
  radSteps = linspace(radStartStop[0],radStartStop[1],nSteps)
  
  iterDM  = fltarr(3,n_elements(haloIDs))
  iterGAS = fltarr(3,n_elements(haloIDs))

  ; load
  h  = loadSnapshotHeader(sP.simPath, snapNum=sP.snap)
  
  pos_gas   = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='pos')
  pos_dm    = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='dm',field='pos')
  
  gas_mass   = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='mass')
  dm_mass    = h.massTable[1]  
  
  ; locate halo
  foreach haloID,haloIDs,j do begin
    haloPos = gc.groupPos[*,haloID]

    ; NaNs in gc?
    if (not finite(haloPos[0])) then begin
      iterDM[*,j]  = [-1,0,0]
      iterGAS[*,j] = [-1,0,0]
      continue
    endif

    ;haloRad = 100.0
    ;if tag_exist(gc,'group_r_crit200') then stop ; catch
    haloRad = gc.group_r_crit200[haloID]

    ; calculate radii and make radial cut
    rad = periodicDists(haloPos,pos_gas,sP=sP)
    
    gas_ind = where(rad le haloRad*radStartStop[0],gas_count)
    pos_gas_sub = pos_gas[*,gas_ind]

    rad = periodicDists(haloPos,pos_dm,sP=sP)

    dm_ind = where(rad le haloRad*radStartStop[0],dm_count)
    pos_dm_sub = pos_dm[*,dm_ind]

    rad = !NULL
    
    ;print,'('+str(haloID)+') Found ['+str(gas_count)+'] gas and ['+str(dm_count)+'] dm  stars inside cut.'
                     
    ; subselect gas masses
    gas_mass_sub = gas_mass[gas_ind]
  
    cm_gas = fltarr(3,nSteps) ;xyz
    cm_dm  = fltarr(3,nSteps) ;xyz
    
    cm_diffs = fltarr(2,nSteps) ;gas,dm
    cm_drift = fltarr(2,nSteps) ;gas,dm
    
    ; compute a center of mass (gas) and compare to haloPos
    foreach radCut,radSteps,k do begin

      ; calculate radii based on last CM
      rad_gas = periodicDists(cm_gas[*,k-1>0],pos_gas_sub,sP=sP)
      rad_dm  = periodicDists(cm_dm[*,k-1>0],pos_dm_sub,sP=sP)

      ; make selection inside this radial cut
      w_gas = where(rad_gas le radCut*haloRad,count_gas)
      w_dm  = where(rad_dm  le radCut*haloRad,count_dm)
      
      ;if (count_dm lt 100) then continue
      
      ; re-calculate center of mass and distance from FoF center
      for i=0,2 do cm_gas[i,k] = total(gas_mass_sub[w_gas] * pos_gas_sub[i,w_gas]) / $
                                 total(gas_mass_sub[w_gas])
      for i=0,2 do cm_dm[i,k]  = mean(pos_dm_sub[i,w_dm])
      
      cm_diffs[*,k] = [sqrt(total((haloPos-cm_gas[*,k])^2.0)),sqrt(total((haloPos-cm_dm[*,k])^2.0))]
      cm_drift[0,k] = sqrt(total((cm_gas[*,k-1>0]-cm_gas[*,k])^2.0))
      cm_drift[1,k] = sqrt(total((cm_dm[*,k-1>0]-cm_dm[*,k])^2.0))
      
      ;debug:
      ;print,k,radCut,count_gas,count_dm
      ;print,'haloPos: ',haloPos
      ;print,'gasCM:   ',cm_gas[*,k]
      ;print,'dmCM:    ',cm_dm[*,k]
      ;print,'cmDiffs:  ',cm_diffs[*,k]/haloRad
    endforeach ;radSteps

    ; save final
    iterDM[*,j]  = cm_dm[*,nSteps-1]
    iterGAS[*,j] = cm_gas[*,nSteps-1]
  endforeach ;haloIDs
  
  ; debug
  ;start_PS,'itercm_'+str(haloID)+'.eps'
  ;  plotsym,0,/fill
  ;  fsc_plot,[0],[0],/nodata,xrange=radStartStop,yrange=[0.0002,4.0],/xs,/ys,/ylog,$
  ;    xtitle="radial search cut [r/"+textoidl("r_{vir}")+"]",$
  ;    ytitle="difference [r/"+textoidl("r_{vir}")+"]",$
  ;    title="iterative CoM convergence (haloID="+str(haloID)+")"
  ;  
  ;  fsc_plot,radSteps,cm_diffs[0,*]/haloRad,line=0,color=getColor(1),/overplot
  ;  fsc_plot,radSteps,cm_diffs[1,*]/haloRad,line=0,color=getColor(2),/overplot
  ;  fsc_plot,radSteps[1:*],cm_drift[0,1:*]/haloRad,line=0,color=getColor(3),/overplot
  ;  fsc_plot,radSteps[1:*],cm_drift[1,1:*]/haloRad,line=0,color=getColor(7),/overplot
  ;  
  ;  legend,['gas FoF delta','dm FoF delta','gas drift','dm drift'],textcolors=getColor([1,2,3,7],/name),$
  ;    box=0,/left,/top
  ;end_PS
  
  r = {iterDM:iterDM,iterGAS:iterGAS}
  return,r
end

; findMatchedHalos(): cross-match halos between two group catalogs

function findMatchedHalos, sP1=sP1, sP2=sP2

  forward_function loadGroupCat

  ; config
  distTol = 10.0 ;kpc
  massTol = 0.1  ;10%
  
  ; load group catalogs
  gc1 = loadGroupCat(sP1.simPath, sP1.snap)
  gc2 = loadGroupCat(sP2.simPath, sP2.snap)

  ; matched Ind in groupcat2 for each FoF halo in groupcat1
  matchedInds = lonarr(gc1.nGroupsTot) - 1
  massDiffs   = fltarr(gc1.nGroupsTot)
  posDiffs    = fltarr(gc1.nGroupsTot)
  
  targetPos   = fltarr(3,gc2.nGroupsTot)
 
  ; load most bound particle ID positions for subfind groups
  idMBCM1 = subgroupPosByMostBoundID(sP=sP1)
  idMBCM2 = subgroupPosByMostBoundID(sP=sP2)
  
  ; fill target positions with mbID positions for each first subgroup of each FoF group in groupcat2
  for i=0,gc2.nGroupsTot-1 do begin
    sgInd = sgPriChildInd(gc=gc2,haloID=i)
    targetPos[*,i] = idMBCM2[*,sgInd]
  endfor

  ; loop over all FoF halos in groupcat1
  for i=0,gc1.nGroupsTot-1 do begin
    ; find primary subgroup child and so position, skip if none
    sgInd = sgPriChildInd(gc=gc1,haloID=i)
    if (sgInd eq -1) then continue
    
    ; compute distances from subfind CM center and cut at maximum
    hPos  = idMBCM1[*,sgInd]
    hMass = gc1.groupMass[i]
    
    dists = periodicDists(hPos,targetPos,sP=sP1)
    minDist = min(dists)
    
    if (minDist le distTol) then begin
      ; select closest halo
      w = where(dists eq minDist,count)
      if (count eq 1) then begin 
        ; check mass agreement
        targInd = w[0]
        massDiff = abs(hMass-gc2.groupMass[targInd]) / hMass
        
        massDiffs[i]   = massDiff
        posDiffs[i]    = minDist        
        
        if (massDiff lt massTol) then begin
          ; select halo
          ;print,'['+str(i)+'] Matched ID ['+str(targInd)+'] at distance = '+string(minDist,format='(f5.2)')+$
          ;      ' Mass Diff = '+string(massDiff,format='(f5.2)')+' '+str(massTol)
          matchedInds[i] = targInd
        endif else begin;mass
          ;print,massDiff,massTol
        endelse
      endif ;count
    endif ;dist
  endfor
  
  wMatch = where(matchedInds ne -1,nMatched)

  r = {matchedInds:matchedInds,massDiffs:massDiffs,posDiffs:posDiffs,wMatch:wMatch,nMatched:nMatched}
  return,r
end

; cosmoCompareHaloCenters(): compare relative differences between different calcuations of the centers
;                            of halos, and verify group catalog consistency (group-subgroup mappings)

pro cosmoCompareHaloCenters

  res = 128
  run = 'dev.tracer.nograd'
  redshift = 3.0
  
  sP = simParams(res=res,run=run,redshift=redshift)
  
  gc    = loadGroupCat(sP.simPath,sP.snap)
  sgCen = subgroupPosByMostBoundID(sP=sP)
  
  binSize = 1.0
  min = 0.0
  max = 50.0
  
  num = 3500
  
  ; group catalog consistency checks
  print,'nGroups nSubgroups',gc.nGroupsTot,gc.nSubgroupsTot
  ; should pass this
  for i=0,gc.nGroupsTot-1 do begin
    sgInd = sgPriChildInd(gc=gc,haloID=i)
    if (sgInd eq -1) then continue
    grNr = gc.subgroupGrNr[sgInd]
    if (grNr ne i) then stop
  endfor
  
  ; centers
  cen_fof  = gc.groupPos[*,0:num]
  
  cen_sfcm = fltarr(3,num)
  cen_mb   = fltarr(3,num)
  
  w = []
  
  for i=0,num-1 do begin
    sgInd = sgPriChildInd(gc=gc,haloID=i)
    if (sgInd ne -1) then begin
      w = [w,i]
      cen_sfcm[*,i] = gc.subgroupCM[*,sgInd]
      cen_mb[*,i] = sgCen[*,sgInd]
    endif
  endfor
  
  ; distances
  dist_fof_sfcm = periodicDists(cen_fof,cen_sfcm,sP=sP)
  dist_fof_mb   = periodicDists(cen_fof,cen_mb,sP=sP)
  dist_sfcm_mb  = periodicDists(cen_sfcm,cen_mb,sP=sP)

  ; plot histos
  hist_fof_sfcm = histogram(dist_fof_sfcm[w],binsize=binSize,min=min,max=max,loc=loc1)
  hist_fof_mb   = histogram(dist_fof_mb[w],binsize=binSize,min=min,max=max,loc=loc2)
  hist_sfcm_mb  = histogram(dist_sfcm_mb[w],binsize=binSize,min=min,max=max,loc=loc3)
  
  start_PS,'hist_fof_sfcm.eps'
    fsc_plot,loc1,hist_fof_sfcm,xtitle="dist [kpc]",ytitle="N"
  end_PS
  
  start_PS,'hist_fof_mb.eps'
    fsc_plot,loc2,hist_fof_mb,xtitle="dist [kpc]",ytitle="N"
  end_PS
  
  start_PS,'hist_sfcm_mb.eps'
    fsc_plot,loc3,hist_sfcm_mb,xtitle="dist [kpc]",ytitle="N"
  end_PS
  
  stop

end

; snapNumToRedshift(): convert snapshot number to redshift or time

function snapNumToRedshift, time=time, all=all, sP=sP, snap=snap

  saveFileName = sP.derivPath + sP.savPrefix + '_snapnum.redshift.sav'

  if not file_test(saveFileName) then stop
  if not keyword_set(snap) then snap = -1
  
  ; restore
  restore,saveFilename

  if (not keyword_set(time)) then begin
    if (snap ge 0 and snap lt n_elements(redshifts)) then $
      return,redshifts[snap]
      
    if (keyword_set(all)) then $
      return,redshifts
  endif else begin
    if (snap ge 0 and snap lt n_elements(redshifts)) then $
      return,times[snap]
      
    if (keyword_set(all)) then $
      return,times
  endelse

end

; codeMassToVirTemp(): convert halo mass (in code units) to virial temperature at specified redshift

function codeMassToVirTemp, mass, redshift, meanmolwt=meanmolwt

  units = getUnits()

  ; mean molecular weight default (valid for ComparisonProject)
  if not keyword_set(meanmolwt) then meanmolwt = 0.6

  ; mass to msun
  mass_msun = mass * (units.UnitMass_in_g / units.Msun_in_g)
  
  ; cosmo
  omega_m   = 0.27
  omega_L   = 0.73
  omega_k   = 0.0
  little_h  = 0.7
  
  omega_m_z = omega_m * (1+redshift)^3.0 / $
              ( omega_m*(1+redshift)^3.0 + omega_L + omega_k*(1+redshift)^2.0 )
  
  Delta_c = 18*!pi^2 + 82*(omega_m_z-1.0) - 39*(omega_m_z-1.0)^2.0

  Tvir = 1.98e4 * (meanmolwt/0.6) * (mass_msun/1e8*little_h)^(2.0/3.0) * $
         (omega_m/omega_m_z * Delta_c / 18.0 / !pi^2.0)^(1.0/3.0) * $
         (1.0 + redshift)/10.0 ;K

  return, Tvir
  
end

; codeMassToLogMsun(): convert mass in code units to log(msun)

function codeMassToLogMsun, mass

  units = getUnits()
  
  mass_msun = mass * (units.UnitMass_in_g / units.Msun_in_g)
  
  ; log of nonzero
  w = where(mass_msun eq 0.0,count)
  if (count ne 0) then $
    mass_msun[w] = 1.0
  
  return,alog10(mass_msun)
end

; convertUtoTemp():

function convertUtoTemp, u, nelec, gamma=gamma, hmassfrac=hmassfrac

  units = getUnits()
  
  ; adiabatic index and hydrogen mass fraction defaults (valid for ComparisonProject)
  if not keyword_set(gamma)     then gamma = 5.0/3.0
  if not keyword_set(hmassfrac) then hmassfrac = 0.76
  
  ; calculate mean molecular weight
  meanmolwt = 4.0/(1.0 + 3.0 * hmassfrac + 4.0* hmassfrac * nelec) * units.mass_proton

  ; calculate temperature
  temp = (gamma-1.0) * u / units.boltzmann * units.UnitEnergy_in_cgs / units.UnitMass_in_g * meanmolwt
  
  return, temp
end

; convertCoolingRatetoCGS():

function convertCoolingRatetoCGS, coolrate, h=h

  units = getUnits()
  
  ; default little h
  if not keyword_set(h) then h = 0.7
  
  ; convert code units (du/dt) to erg/s/g (cgs)
  coolrate_cgs = coolrate * units.UnitEnergy_in_cgs * units.UnitTime_in_s^(-1.0) * $
                 units.UnitMass_in_g^(-1.0) * h
                 
  return, coolrate_cgs
end

; calcEntropy():

function calcEntropy, u, dens, gamma=gamma

  ; adiabatic index default (valid for ComparisonProject)
  if not keyword_set(gamma)     then gamma = 5.0/3.0
  
  pressure = (gamma-1.0) * u * dens
  entropy  = pressure / (dens^gamma)
  
  return, entropy
end

; rhoRatioToCrit():

function rhoRatioToCrit, rho, omega_b=omega_b, redshift=redshift

  units = getUnits()
  
  ; cosmo
  omega_m   = 0.27
  omega_L   = 0.73
  omega_k   = 0.0
  
  ; default omega_b (valid for ComparisonProject)
  if not keyword_set(omega_b) then omega_b = 0.044
  
  rho_b = omega_b * units.rhoCrit
  
  ; scale for redshift other than zero
  if keyword_set(redshift) then begin
    H_z_fact = ( omega_m*(1+redshift)^3.0 + omega_L + omega_k*(1+redshift)^2.0 )
    rho_b *= H_z_fact
  endif
  
  return, rho/rho_b

end

; redshiftToAge(): convert redshift to age of the universe (approximate)
function dtdz, z, lambda0 = lambda0, q0 = q0
  term1 = (1.0d + z)
  term2 = 2.0d * (q0 + lambda0) * z + 1.0d - lambda0
  term3 = (1.0d + z) * (1.0d +z)
  return, 1.0 / (term1 * sqrt(term2 * term3 + lambda0))
end
   
function redshiftToAge, z

  units = getUnits()
  
  ; config
  zform = 1000.0
  H0 = 70.0
  k = 0.0
  Omega_m = 0.27
  Lambda0 = 0.73
  q0 = -0.55
  
  ; arrays
  nz  = N_elements(z)
  age = z * 0.0
  
  ; integrate with qsimp
  for i= 0L, nz-1 do begin
    if (z[i] ge zform) then age_z = 0 else $
        qsimp,'dtdz', z[i], zform, age_z, q0 = q0, lambda0 = lambda0
    age[i] = age_z
  endfor

  return, age * 3.085678e+19 / 3.15567e+7 / H0 / 1e9 ;Gyr

end

function snapNumToAge, snap
  z = snapNumToRedshift(snapNum=snap)
  return, redshiftToAge(z)
end
  
; rhoTHisto(): make mass-weighted density-temperature 2d histogram

function rhoTHisto, dens_in, temp_in, mass=mass, nbins=nbins, plot=plot

  units = getUnits()

  ; config
  if not keyword_set(nbins) then nbins = 20
  dens = alog10(rhoRatioToCrit(dens_in))
  temp = alog10(temp_in)
  rMinMax = [-2.0,8.0]
  tMinMax = [3.0,7.0]
  
  ; calculate bin sizes
  binSizeRho  = (rMinMax[1]-rMinMax[0]) / nbins
  binSizeTemp = (tMinMax[1]-tMinMax[0]) / nbins
  
  if not keyword_set(mass) then begin
    ; hist_2d (no weighting)
    h2rt = hist_2d(dens,temp,bin1=binSizeRho,bin2=binSizeTemp,$
                   min1=rMinMax[0],min2=tMinMax[0],max1=rMinMax[1]+binSizeRho,max2=tMinMax[1]+binSizeTemp)
    ; h2rt /= float(max(h2rt)) ;norm s.t. colorbar gives fraction wrt max cell
  endif else begin
    ; hist2d (weighted)
    h2rt = hist2d(dens,temp,mass,$
                  binsize1=binsizeRho,binsize2=binSizeTemp,$
                  min1=rMinMax[0],min2=tMinMax[0],max1=rMinMax[1]+binSizeRho,max2=tMinMax[1]+binSizeTemp)
  endelse

  ; plot
  if keyword_set(plot) then begin
    ; color table
    loadct, 2, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
      
    tvim,h2rt,pcharsize=!p.charsize-1.0,scale=1,clip=[10,100],$;,/c_map
         xtitle="log ("+textoidl("\rho / \rho_{crit}")+")",ytitle="log (T [K])",$
         stitle="Total Mass (Msun)",barwidth=0.5,lcharsize=!p.charsize-1.5,$
         xrange=[-2.0,8.0],yrange=[3.0,7.0],$;xrange=rMinMax,yrange=tMinMax,$
         /rct;,nodata=0,rgb_nodata=[1.0,1.0,1.0] ;display zeros as white not black
  endif
  
  return,h2rt

end

; redshift_axis(): draw redshift axis

pro redshift_axis, xRange, yRange, ylog=ylog, snapnum=snapnum, dotted=dotted, zTicknames=zTicknames

  logFac = 1.0
  if keyword_set(ylog) then logFac = 10.0
  
  yTickLen = (yRange[1]-yRange[0]) / 40.0 * logFac
  yTextOff = (yRange[1]-yRange[0]) / 50.0 * logFac

  if (not keyword_set(zTicknames)) then $
    zTicknames = ['30','6','4','3','2','1','0.5','0.25','0']
  nZ = n_elements(zTicknames)
  
  zXPos = fltarr(nZ)
  
  ; plot "maximum" redshift label
  if (xRange[0] le 0.0) then $
    fsc_text,0.0,yRange[1]+yTextOff,zTicknames[0],alignment=0.5
  
  ; skip z=30 (highest) at t=0
  for i=1,nZ-1 do begin
    if keyword_set(snapnum) then begin ;x-axis in snapshot number
      zXPos[i] = redshiftToSnapNum(float(zTicknames[i]))
    endif else begin ;x-axis in time [gyr]
      zXPos[i] = redshiftToAge(float(zTicknames[i]))
    endelse
    
    ; plot tick mark and label if inside plotrange
    if (zXPos[i] ge xRange[0] and zXPos[i] le xRange[1]) then begin
      fsc_plot,[zXPos[i],zXPos[i]],[yRange[1],yRange[1]-yTickLen],/overplot
      fsc_text,zXPos[i],yRange[1]+yTextOff,zTicknames[i],alignment=0.5
    endif
    
    ; plot vertical dotted line at each redshift mark if requested
    if keyword_set(dotted) then $
      fsc_plot,[zXPos[i],zXPos[i]],yRange,line=dotted,/overplot,thick=!p.thick-0.5
  endfor
  
  fsc_plot,xRange,[yRange[1],yRange[1]],/overplot
  fsc_text,0.5,0.94,"Redshift",/normal

end

; universeage_axis(): draw age of universe (elapsed) axis

pro universeage_axis, xRange, yRange, ylog=ylog, dotted=dotted

  logFac = 1.0
  if keyword_set(ylog) then logFac = 10.0
  
  yTickLen = (yRange[1]-yRange[0]) / 40.0 * logFac
  yTextOff = (yRange[1]-yRange[0]) / 50.0 * logFac

  if (not keyword_set(zTicknames)) then begin
    zTicknames = ['0.1','1','1.5','2','3','4','5','6','8','13.8'] ;10=0.33
    zRedshifts = [30,5.7,4.21,3.25,2.23,1.65,1.265,0.98,0.595,0.0]
  endif
  nZ = n_elements(zTicknames)
  
  zXPos = fltarr(nZ)
  
  ; plot "maximum" redshift label
  if (xRange[0] le 0.0) then $
    fsc_text,0.0,yRange[1]+yTextOff,zTicknames[0],alignment=0.5
  
  ; skip t=0 (highest) at z=inf (?)
  for i=1,nZ-1 do begin

    zXPos[i] = zRedshifts[i]
    
    ; plot tick mark and label if inside plotrange
    if (zXPos[i] le xRange[0] and zXPos[i] ge xRange[1]) then begin
      fsc_plot,[zXPos[i],zXPos[i]],[yRange[1],yRange[1]-yTickLen],/overplot
      fsc_text,zXPos[i],yRange[1]+yTextOff,zTicknames[i],alignment=0.5
    endif
    
    ; plot vertical dotted line at each redshift mark if requested
    if keyword_set(dotted) then $
      fsc_plot,[zXPos[i],zXPos[i]],yRange,line=dotted,/overplot,thick=!p.thick-0.5
  endfor
  
  fsc_plot,xRange,[yRange[1],yRange[1]],/overplot
  fsc_text,(xRange[0]-xRange[1])/2.0,yRange[1]+yTextOff*5.0,"Time [Gyr]",alignment=0.5

end
