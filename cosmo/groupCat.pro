; groupCat.pro
; cosmological simulations - group (fof/subfind) catalog utilities
; dnelson sep.2013

; getMatchedIDs(): return matched subgroup IDs given a "haloID" from the HaloComparisonProject
;                  account for resolution, redshift, and run. try to prevent IDs from floating
;                  around in various functions

function getMatchedIDs, sPa=sPa, sPg=sPg, simParams=simParams, haloID=haloID, mosaicIDs=mosaicIDs
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if keyword_set(sPa) and keyword_set(sPg) then begin
    if sPa.res ne sPg.res then message,'Error: Resolutions differ.'
    if abs(sPa.redshift - sPg.redshift) gt 0.02 then message,'Error: Redshifts differ.'
  endif
  
  ; hash container
  sgIDs = hash()
  
  ; return 8 matched IDs for mosaic?
  if keyword_set(mosaicIDs) then begin
    ; 512^3 z=2
    sgIDs['z2_res512_tracer']   = [5611,4518,2874,2389,2058,1037,816,252]
    sgIDs['z2_res512_gadget']   = [6369,5151,3338,2824,2538,1117,981,266]
    sgIDs['z2_res512_feedback'] = [5751,4678,2910,2400,2233,1040,639,252]
    sgIDs['z2_res512_axes']     = list([0,2],[0,1],[0,2],[0,2],[1,2],[0,2],[1,2],[0,1])
    
    ; 256^3 z=2 (matched from 512)
    sgIDs['z2_res256_tracer']        = [1527,1229,733,605,539,203,169,28]
    sgIDs['z2_res256_gadget']        = [1888,1565,983,845,758,330,278,74]
    sgIDs['z2_res256_feedback']      = [1502,1107,694,559,517,235,160,51]
    sgIDs['z2_res256_feedback_noz']  = [1545,1162,770,647,560,219,200,60]
    sgIDs['z2_res256_feedback_nofb'] = [1604,1257,778,637,543,196,162,58]
    sgIDs['z2_res256_axes']          = list([0,2],[0,1],[0,2],[0,2],[1,2],[0,2],[1,2],[0,1])
    
    ; 512^3 z=3 (~same mass spacing as z2, but randomly selected at z=3 for run=tracer)
    sgIDs['z3_res512_tracer']   = [3664,2920,1812,1459,1210,620,496,137]
    sgIDs['z3_res512_gadget']   = [4398,3194,2061,1593,1275,742,573,154] ; matched from tracer
    sgIDs['z3_res512_feedback'] = [4434,2962,1735,1409,1238,651,447,134] ; matched from tracer
    sgIDs['z3_res512_axes']     = list([0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1])
    
    if keyword_set(simParams) then begin
      key  = 'z' + str(fix(simParams.redshift)) + '_res' + str(simParams.res) + '_' + simParams.run
      keyAxes = 'z' + str(fix(simParams.redshift)) + '_res' + str(simParams.res) + '_axes'
      if ~sgIDs.HasKey(key) then message,'Error: Unknown request.'
      
      return, { gcIDs : sgIDs[key], axes : sgIDs[keyAxes] }
    endif
    
    keyA = 'z' + str(fix(sPa.redshift)) + '_res' + str(sPa.res) + '_' + sPa.run
    keyG = 'z' + str(fix(sPg.redshift)) + '_res' + str(sPg.res) + '_' + sPg.run
    keyAxes = 'z' + str(fix(sPa.redshift)) + '_res' + str(sPa.res) + '_axes'
    
    if ~sgIDs.HasKey(keyA) or ~sgIDs.HasKey(keyG) then message,'Error: Unknown request.'
    
    return, { gcIDsA : sgIDs[keyA], gcIDsG : sgIDs[keyG], axes : sgIDs[keyAxes] }
  endif
  
  ; z2.314 (disappearing fan along bottom, all-sky views for paper)
  sgIDs['z2_h314_res512_tracer']   = 816
  sgIDs['z2_h314_res512_gadget']   = 981
  sgIDs['z2_h314_res512_arepo']    = 879 ;had as 927 in comparison project?
  sgIDs['z2_h314_res512_feedback'] = 639
  
  sgIDs['z2_h314_res256_tracer']   = 169
  sgIDs['z2_h314_res256_gadget']   = 278
  sgIDs['z2_h314_res256_feedback'] = 160
  sgIDs['z2_h314_res256_feedback_noz'] = 200
  sgIDs['z2_h314_res256_feedback_nofb'] = 162
  
  sgIDs['z2_h314_res128_tracer']   = 50
  sgIDs['z2_h314_res128_gadget']   = 71
  sgIDs['z2_h314_res128_feedback'] = 56
  
  ; z3 h314 (tracked back)
  sgIDs['z3_h314_res512_tracer']   = 183
  sgIDs['z3_h314_res512_gadget']   = 211
  sgIDs['z3_h314_res512_feedback'] = 184
  
  sgIDs['z3_h314_res256_tracer']   = 46
  sgIDs['z3_h314_res256_gadget']   = 51
  sgIDs['z3_h314_res256_feedback'] = 45
  
  sgIDs['z3_h314_res128_tracer']   = 16
  sgIDs['z3_h314_res128_gadget']   = 19
  sgIDs['z3_h314_res128_feedback'] = 10
    
  ; z2.304 (four filaments aligned in a plus pattern, high mass example for paper)
  sgIDs['z2_h304_res512_tracer']   = 2004
  sgIDs['z2_h304_res512_gadget']   = 2342
  sgIDs['z2_h304_res512_arepo']    = 2132
  sgIDs['z2_h304_res512_feedback'] = 2020
  
  sgIDs['z2_h304_res256_tracer']   = 510
  sgIDs['z2_h304_res256_gadget']   = 673
  sgIDs['z2_h304_res256_feedback'] = 458
  sgIDs['z2_h304_res256_feedback_noz'] = 518
  sgIDs['z2_h304_res256_feedback_nofb'] = 524
  
  sgIDs['z2_h304_res128_tracer']   = 150
  sgIDs['z2_h304_res128_gadget']   = 217
  sgIDs['z2_h304_res128_feedback'] = 151
  
  ; z3 h304 (tracked back from z=2)
  sgIDs['z3_h304_res512_tracer']   = 1669
  sgIDs['z3_h304_res512_gadget']   = 1902
  sgIDs['z3_h304_res512_arepo']    = 1657
  sgIDs['z3_h304_res512_feedback'] = 1621
  
  sgIDs['z3_h304_res256_tracer']   = 271
  sgIDs['z3_h304_res256_gadget']   = 533
  sgIDs['z3_h304_res256_feedback'] = 385
  
  sgIDs['z3_h304_res128_tracer']   = 91
  sgIDs['z3_h304_res128_gadget']   = 162
  sgIDs['z3_h304_res128_feedback'] = 100
  
  ; z0 h304 (tracked forward from z=2)
  sgIDs['z0_h304_res256_tracer']   = 1310
  sgIDs['z0_h304_res256_gadget']   = 0
  sgIDs['z0_h304_res256_feedback'] = 1379
  
  ; z2.301
  sgIDs['z2_h301_res512_tracer']   = 1948
  sgIDs['z2_h301_res512_gadget']   = 2289
  sgIDs['z2_h301_res512_arepo']    = 2034
  sgIDs['z2_h301_res512_feedback'] = 2076
  
  ; z2.130 (low mass example for paper)
  sgIDs['z2_h130_res512_tracer']   = 5611
  sgIDs['z2_h130_res512_gadget']   = 6369
  sgIDs['z2_h130_res512_arepo']    = 5966
  sgIDs['z2_h130_res512_feedback'] = 5755
  
  sgIDs['z2_h130_res256_tracer']   = 1527
  sgIDs['z2_h130_res256_gadget']   = 1888
  sgIDs['z2_h130_res256_feedback'] = 1502
  
  sgIDs['z2_h130_res128_tracer']   = 619
  sgIDs['z2_h130_res128_gadget']   = 802
  sgIDs['z2_h130_res128_feedback'] = 604
  
  ; z2.64
  sgIDs['z2_h64_res512_gadget']    = 5498
  sgIDs['z2_h64_res512_arepo']     = 5097
  sgIDs['z2_h64_res512_feedback']  = 5363
  
  ; zoom.0 (target halo)
  sgIDs['z2_h0_res9_zoom_20mpc_hind0']  = 0
  sgIDs['z2_h0_res10_zoom_20mpc_hind0'] = 0
  sgIDs['z2_h0_res9_zoom_20mpc_derefgal_hind0'] = 0
  sgIDs['z2_h0_res9_zoom_20mpc_derefgal_nomod_hind0'] = 0

  ; illustris.0 (most massive halo in box at z=0)
  sgIDs['snap123_h0_res455_illustris']  = 0
  sgIDs['snap123_h0_res910_illustris']  = 0
  sgIDs['snap123_h0_res1820_illustris'] = 0
  
  sgIDs['z0_h0_res128_feedback'] = 0 ; fake it for testing
  sgIDs['z0_h0_res256_feedback'] = 0 ; fake it for testing
  
  ; illustris.1000 (~10^12 example at z=0)
  sgIDs['snap123_h1000_res455_illustris']  = 15912
  sgIDs['snap123_h1000_res910_illustris']  = 70342
  sgIDs['snap123_h1000_res1820_illustris'] = 356960
  
  ; single simulation requested?
  if keyword_set(simParams) then begin
    if simParams.redshift ge 0.0 then $ ; z*
      key = 'z' + str(fix(simParams.redshift)) + '_h' + str(haloID) + '_res' + str(simParams.res) + $
            '_' + simParams.run
          
    if simParams.redshift lt 0.0 then $ ;snap*
      key = 'snap' + str(simParams.snap) + '_h' + str(haloID) + '_res' + str(simParams.res) + $
            '_' + simParams.run
            
    if simParams.zoomLevel ne 0 then key += '_hind'+str(simParams.hInd)
    
    if ~sgIDs.HasKey(key) then message,'Error: Unknown request.'
    
    return, sgIDs[key]
  endif
  
  ; make key and pull the requested subgroup ID out
  keyA = 'z' + str(fix(sPa.redshift)) + '_h' + str(haloID) + '_res' + str(sPa.res) + '_' + sPa.run
  keyG = 'z' + str(fix(sPg.redshift)) + '_h' + str(haloID) + '_res' + str(sPg.res) + '_' + sPg.run
    
  if sPa.zoomLevel ne 0 then keyA += '_hind'+str(sPa.hInd)
  if sPg.zoomLevel ne 0 then keyG += '_hind'+str(sPg.hInd)
  
  if (sPa.zoomLevel eq 0 and sPg.zoomLevel ne 0) or (sPa.zoomLevel ne 0 and sPg.zoomLevel eq 0) then $
    message,'Error: One zoom and one non zoom? No comparable halos.'
  if (sPa.zoomLevel eq 0 or sPg.zoomLevel eq 0) and haloID eq 0 then message,'Error: haloID=0 for zooms.'
  if ~sgIDs.HasKey(keyA) and ~sgIDs.HasKey(keyG) then message,'Error: Unknown request.'
    
  ; make return
  r = { a:9999999, g:9999999, axes:[0,1] }
    
  if sgIDs.HasKey(keyA) then r.a = sgIDs[keyA]
  if sgIDs.HasKey(keyG) then r.g = sgIDs[keyG]
  
  return,r

end

; rematch(): crossmatch a few subhalos between different sims

pro rematch

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function loadGroupCat

  sP1 = simParams(res=455,run='illustris',snap=123)
  sP2 = simParams(res=1820,run='illustris',snap=123)

  gc1 = loadGroupCat(sP=sP1,/skipIDs,/skipOffsets)
  gc2 = loadGroupCat(sP=sP2,/skipIDs,/skipOffsets)

  known1gcIDs = [15912]
  
  foreach gcID,known1gcIDs do begin
    ; find closest spatial match
    dists = periodicDists(gc1.subgroupPos[*,gcID],gc2.subgroupPos,sP=sP1)
    w = where(dists eq min(dists))
    
    ; sort dists to find second closest
    dists_sort = sort(dists)
    dists = dists[dists_sort]
    
    ; report closest
    print,'Matched known ['+str(gcID)+'] to new ['+str(w[0])+'].'
    print,'Masses ',gc1.subgroupMass[gcID],gc2.subgroupMass[w[0]]
    print,'Distance of match: ['+string(min(dists),format='(f5.1)')+'] and second closest: ['+$
      string(dists[1],format='(f5.1)')+'] with mass: '+$
      string(gc2.subgroupMass[dists_sort[1]],format='(f6.2)')+' ind: ['+str(dists_sort[1])+']'
    print,''
  endforeach
  
  stop

end

; gcPriChildInd(): get the subgroup index of the primary subgroup of a given group, or -1 if none

function gcPriChildInd, gc=gc, haloID=haloID

  if (n_elements(gc) eq 0 or n_elements(haloID) eq 0) then stop

  gcInd = min(where(gc.subgroupGrNr eq haloID,count))
  fsInd = gc.groupFirstSub[haloID]
  
  ; skip if group has no subgroup to obtain this value from
  if (count eq 0) then return,-1
  
  ; double-check (failing either should indicate a problem in the group catalogs, right?)
  if (gcInd ne fsInd) then return,-1 ;stop
  if (gc.subgroupGrNr[gcInd] ne haloID) then return,-1 ;stop

  return,gcInd

end

; gcIDList(): return a sub-list of subgroup IDs from a given group catalog gc
;             note these subgroup IDs are the same used to index a galaxy catalog
; 
; select: one of the following
;   'pri' : members of the first subgroup of each group only ("background"/"main subhalo"/"halo")
;   'sec' : members of the non-first subgroups of each group only ("satellites"/"subhalos")
;   'all' : all subgroups
;
; massRange=[min,max] : if input, restrict to halos in this log(msun) mass range

function gcIDList, sP=sP, gc=gc, select=select, massRange=massRange

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function loadGroupCat

  if ~keyword_set(select) then message,'Error: Must specify select.'
  if select ne 'pri' and select ne 'sec' and select ne 'all' then message,'Error: Unknown select.'
    
  ; load galaxy cat if necessary
  if not keyword_set(gc) then begin
    if not keyword_set(sP) then begin
      message,'Error: Must specify gc or sP.'
    endif
    gc = loadGroupCat(sP=sP,/skipIDs)
  endif

  prevGrNr   = -1
  valGCids   = []

  if (select eq 'pri') then begin
  
    ; "background"/"main subhalos"/"halos" only
    for i=0L,n_elements(gc.subgroupLen)-1 do begin
      if (gc.subgroupGrnr[i] eq prevGrNr) then begin
        prevGrNr = gc.subgroupGrnr[i]
      endif else begin
      
        if (keyword_set(minNumPart)) then begin
          if (gc.subgroupLen[i] ge minNumPartVal) then $
            valGCids = [valGCids,i]
        endif else begin
          valGCids = [valGCids,i]
        endelse ;minNumPart
        
        prevGrNr = gc.subgroupGrnr[i]
      endelse
    endfor
    
  endif
  
  if (select eq 'sec') then begin
  
    ; "satellites"/"subhalos" only
    for i=0L,n_elements(gc.subgroupLen)-1 do begin
      if (gc.subgroupGrnr[i] ne prevGrNr) then begin
        prevGrNr = gc.subgroupGrnr[i]
      endif else begin
      
        if (keyword_set(minNumPart)) then begin
          if (gc.subgroupLen[i] ge minNumPartVal) then $
            valGCids = [valGCids,i]
        endif else begin
          valGCids = [valGCids,i]
        endelse ;minNumPart
        
      endelse
    endfor
    
  endif
  
  if (select eq 'all') then begin
  
    ; both primary and secondary
    for i=0L,n_elements(gc.subgroupLen)-1 do begin
    
      if (keyword_set(minNumPart)) then begin
        if (gc.subgroupLen[i] ge minNumPartVal) then $
          valGCids = [valGCids,i]
      endif else begin
        valGCids = [valGCids,i]
      endelse ;minNumPart
      
    endfor
    
  endif
  
  ; restrict by mass
  if n_elements(massRange) gt 0 then begin
    gcMasses = codeMassToLogMsun(gc.subgroupMass)
    w = where(gcMasses[valGCids] ge massRange[0] and gcMasses[valGCids] le massRange[1],count)
    if count eq 0 then message,'Error'
    valGCids = valGCids[w]
  endif
  
  return,valGCids

end

; gcPIDList(): return a list of member particle IDs from a given group catalog gc
;
; select: one of the following
;   'pri' : members of the first subgroup of each group only ("background"/"main subhalo"/"halo")
;   'sec' : members of the non-first subgroups of each group only ("satellites"/"subhalos")
;   'all' : all subgroups
; valGCids: return member ids from only these groups
;
; partType=N     : use subgroup type offset table to return only IDs of a specific particle type
; partType='all' : return all particle IDs in that subgroup regardless of type

function gcPIDList, gc=gc, select=select, valGCids=valGCids, partType=PT

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if (n_elements(gc) eq 0 or (n_elements(select) eq 0 and n_elements(valGCids) eq 0) or $
      n_elements(PT) eq 0) then message,'Error: gcPIDList: Bad inputs.'

  ; get list of appropriate group ids
  if n_elements(valGCids) eq 0 then valGCids = gcIDList(gc=gc,select=select)
  
  ; make list of particle ids in these groups
  start = 0L
  
  partType = strlowcase(string(PT)) ; so we don't change the input
  
  if strcmp(partType,'all') then begin
    ; allocate return
    if size(gc.IDs,/tname) eq 'LONG'   then subgroupPIDs  = lonarr(total(gc.subGroupLen[valGCids],/int))
    if size(gc.IDs,/tname) eq 'LONG64' then subgroupPIDs  = lon64arr(total(gc.subGroupLen[valGCids],/int))
    if size(gc.IDs,/tname) eq 'ULONG64' then subgroupPIDs = ulon64arr(total(gc.subGroupLen[valGCids],/int))
    
    foreach gcID, valGCids do begin
      ; select particle IDs in subgroup
      subgroupPIDs[start:start+gc.subGroupLen[gcID]-1] = $
        gc.IDs[gc.subGroupOffset[gcID] : gc.subGroupOffset[gcID] + gc.subGroupLen[gcID] - 1]
      start += gc.subGroupLen[gcID]
    endforeach
    
    if start ne n_elements(subgroupPIDs) then message,'Error: Failed to locate all partTypes.'
    if min(subgroupPIDs) lt 0 then stop ; check 32 bit long overflow
    return, subgroupPIDs
  endif
  
  ; convert partType to number if it is a name
  partType = partTypeNum(partType)

  ; check if this particle type is present in the subgroup selection
  if total(gc.subgroupLenType[partType,valGCids] gt 0) then begin
    ; allocate return
    if size(gc.IDs,/tname) eq 'LONG'    then subgroupPIDs = lonarr(total(gc.subGroupLenType[partType,valGCids],/int))
    if size(gc.IDs,/tname) eq 'LONG64'  then subgroupPIDs = lon64arr(total(gc.subGroupLenType[partType,valGCids],/int))
    if size(gc.IDs,/tname) eq 'ULONG64' then subgroupPIDs = ulon64arr(total(gc.subGroupLenType[partType,valGCids],/int))
    
    ; store particle IDs of this type from each subgroup
    foreach gcID, valGCids do begin
      ; may be none
      if gc.subgroupLenType[partType,gcID] gt 0 then begin
        subgroupPIDs[start:start+gc.subGroupLenType[partType,gcID]-1] = $
          gc.IDs[gc.subGroupOffsetType[partType,gcID] : gc.subGroupOffsetType[partType,gcID] + $
          gc.subGroupLenType[partType,gcID] - 1]
        
        start += gc.subGroupLenType[partType,gcID]
      endif
    endforeach

    if start ne n_elements(subgroupPIDs) then message,'Error: Failed to locate all of this partType.'
    if min(subgroupPIDs) lt 0 then stop ; check 32 bit long overflow
    return, subgroupPIDs    
  endif
  
  print,'Warning! Empty gcPIDList return.'
  return,[]

end

; massTargetToHaloID(): return the halo ID (subgroupInd) nearest to the target mass in log(Msun)

function massTargetToHaloID, hMassTargets, sP=sP, verbose=verbose

  if ~keyword_set(sP) then message,'Error'
  
  ; load group catalog and calculate log(M) masses of all subgroups of the requested type
  gc = loadGroupCat(sP=sP,/skipIDs)
  priSGIDs = gcIDList(gc=gc,select='pri')
  hMasses = codeMassToLogMsun(gc.subgroupMass[priSGIDs])

  ; locate nearest masses to requested masses in a consistent way
  hInds = value_locate(hMasses,hMassTargets)
  w = where(hInds eq -1,count)
  if count gt 0 then hInds[w] = 0 ; largest mass halo available
 
  ; convert indices in hMasses to sgIDs
  subgroupIDs = priSGIDs[hInds] 

  if keyword_set(verbose) then $
      print,'selected halo ind ['+str(hInd)+'] sgID ['+str(subgroupID)+'] mass = '+$
        string(hMasses[hInd],format='(f5.2)')

  return,subgroupIDs
end

; subgroupPosByMostBoundID(): compute a "best" center position in space for all subfind groups by 
;                             using the position of the most bound particles, whose IDs are stored in 
;                             the group catalog but without knowing the particle types we have to 
;                             load all gas+dm+stars+bh particle positions

function subgroupPosByMostBoundID, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function loadSnapshotSubset, loadGroupCat
  
  ; save/restore
  saveFilename = sP.derivPath+'sgCen.mbID/sgCen.mbID.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.sav'
                 
  if file_test(saveFilename) then begin
    restore,saveFilename
    return, groupCen
  endif
  
  h = loadSnapshotHeader(sP=sP)
  gc = loadGroupCat(sP=sP,/skipIDs)
  
  groupCen = fltarr(3,gc.nSubgroupsTot)
  
  partTypes = ['gas','dm','star','bh','tracervel']
  
  numFound = 0L
  
  foreach partType,partTypes do begin
    ; skip if no particles of this type, or if we have all our matches already
    if h.nPartTot[partTypeNum(partType)] eq 0 then continue
    if numFound eq gc.nSubgroupsTot then continue
    
    ; load ids and find matches for this type
    ids = loadSnapshotSubset(sP=sP,partType=partType,field='ids')
    calcMatch,gc.subgroupIdMostBound,ids,gc_ind,ids_ind,count=count
    
    numFound += count
    
    if count gt 0 then begin
      ; load positions for matched most bound particles
      sort_gc_ind = calcSort(gc_ind)
      ids_ind = ids_ind[sort_gc_ind]
      gc_ind  = gc_ind[sort_gc_ind]

      pos = loadSnapshotSubset(sP=sP,partType=partType,field='pos',inds=ids_ind)
      groupCen[*,gc_ind] = pos
    endif
  endforeach

  if numFound ne gc.nSubgroupsTot then message,'Error: Failed to find all most bound IDs.'
  
  ; save
  save,groupCen,filename=saveFilename
  
  return, groupCen
end

; groupCenterPosByIterativeCM(): compute a "better" center position in space for all FoF groups by
;                                iteratively searching for the center of mass

function groupCenterPosByIterativeCM, sP=sP, gc=gc, haloIDs=haloIDs

  compile_opt idl2, hidden, strictarr, strictarrsubs
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
  h  = loadSnapshotHeader(sP=sP)
  
  pos_gas   = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  pos_dm    = loadSnapshotSubset(sP=sP,partType='dm',field='pos')
  
  gas_mass   = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
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
      
    endforeach ;radSteps

    ; save final
    iterDM[*,j]  = cm_dm[*,nSteps-1]
    iterGAS[*,j] = cm_gas[*,nSteps-1]
  endforeach ;haloIDs
  
  r = {iterDM:iterDM,iterGAS:iterGAS}
  return,r
end

; findMatchedHalos(): cross-match halos between two group catalogs
; (different runs or resolutions of the same ICs where the same halos will be in nearly the same spots)

function findMatchedHalos, sP1=sP1, sP2=sP2

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function loadGroupCat

  ; config
  distTol = 40.0 ; kpc
  massTol = 0.4  ; 40%
  
  ; save/restore
  saveFilename = sP1.derivPath + 'matchCat_'+sP2.run+'_'+str(sP2.res)+'-'+str(sP1.snap)+'_m'+$
                 str(fix(massTol*10))+'_d'+str(fix(distTol))+'.sav'
    
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif  
  
  ; load group catalogs
  gc1 = loadGroupCat(sP=sP1,/skipIDs)
  gc2 = loadGroupCat(sP=sP2,/skipIDs)

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
    gcInd = gcPriChildInd(gc=gc2,haloID=i)
    targetPos[*,i] = idMBCM2[*,gcInd]
  endfor

  ; loop over all FoF halos in groupcat1
  for i=0,gc1.nGroupsTot-1 do begin
    ; find primary subgroup child and so position, skip if none
    gcInd = gcPriChildInd(gc=gc1,haloID=i)
    if (gcInd eq -1) then continue
    
    ; compute distances from subfind CM center and cut at maximum
    hPos  = idMBCM1[*,gcInd]
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
  
  ; save
  save,r,filename=saveFilename
  
  return,r
end

; cosmoVerifySubfindIntegrity(): verify ordering etc within a subfind catalog

pro cosmoVerifySubfindIntegrity

  sP = simParams(res=512,run='tracer',redshift=2.0)
  
m=189
  ;for m=50,80 do begin
    ; load
    sP.snap = m
    print,sP.snap
    gc = loadGroupCat(sP=sP,/readIDs)
    
    gas_ids  = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    dm_ids   = loadSnapshotSubset(sP=sP,partType='dm',field='ids')
    star_ids = loadSnapshotSubset(sP=sP,partType='stars',field='ids')
    
    ; general
    if gc.nGroupsTot ne n_elements(gc.groupLen) then message,'Fail nG'
    if gc.nSubgroupsTot ne n_elements(gc.subgroupLen) then message,'Fail nSG'
    if gc.nIDsTot ne n_elements(gc.IDs) then message,'Fail IDs'
    
    if total(gc.groupNSubs,/int) ne gc.nSubgroupsTot then message,'Fail NS'
    if max(gc.groupFirstSub) gt gc.nSubgroupsTot then message,'Fail gFS'
    if max(gc.subgroupGrNr) gt gc.nGroupsTot then message,'Fail gGN'
    if max(gc.subgroupParent) gt gc.nSubgroupsTot then message,'Fail sP'
 
    if max(gc.groupOffset) gt gc.nIDsTot then message,'Fail O1'
    if max(gc.groupOffsetType) gt gc.nIDsTot then message,'Fail O2'
    if max(gc.subgroupOffset) gt gc.nIDsTot then message,'Fail 03'
    if max(gc.subgroupOffsetType) gt gc.nIDsTot then message,'Fail 04'

    if total(gc.groupLen,/int) ne gc.nIDsTot then message,'Fail tot1'
    if total(gc.groupLenType,/int) ne gc.nIDsTot then message,'Fail tot2'
    if total(gc.subgroupLen,/int) ne total(gc.subgroupLenType,/int) then message,'Fail tot3'   

    ; verify LenType (OffsetType)
    testIDs = gcPIDList(gc=gc,select='all',partType='gas')
    match,testIDs,gas_ids,ind1,ind2,count=count1
    if count1 ne n_elements(testIDs) then message,'Fail 1'
    
    testIDs = gcPIDList(gc=gc,select='all',partType='dm')
    match,testIDs,dm_ids,ind1,ind2,count=count2
    if count2 ne n_elements(testIDs) then message,'Fail 21'
    
    testIDs = gcPIDList(gc=gc,select='all',partType='stars')
    match,testIDs,star_ids,ind1,ind2,count=count3
    if count3 ne n_elements(testIDs) then message,'Fail 3'
    
    testIDs = gcPIDList(gc=gc,select='all',partType='all')
    if n_elements(testIDs) ne (count1+count2+count3) then message,'Fail 4'
    
    ; verify FirstSub, NSubs
    for i=0L,gc.nGroupsTot-1 do begin
      if gc.groupNSubs[i] gt 0 then begin
        FirstSub = gc.groupFirstSub[i]
        LastSub  = gc.groupFirstSub[i] + gc.groupNSubs[i] - 1
        ParIDs = gc.subgroupGrNr[FirstSub:LastSub]
        if nuniq(ParIDs) gt 1 or ParIDs[0] ne i then message,'Fail 5'
      endif
      
      ; check LenType totals
      if total(gc.groupLenType[*,i],/int) ne gc.groupLen[i] then message,'Fail gLT'
    endfor
    
    ; verify GrNr ordering and SubLen ordering
    prevGrNr = 0L
    prevSGLen = 1000000000LL
    
    for i=0L,gc.nSubgroupsTot-1 do begin
      if prevGrNr ne gc.subgroupGrNr[i] then prevSGLen = 1000000000LL
      
      if gc.subgroupGrNr[i] lt prevGrNr then message,'Fail GrNrOrder'
      if gc.subgroupLen[i] gt prevSGLen then message,'Fail SubLenorder'
      
      prevSGLen = gc.subgroupLen[i]
      prevGrNr = gc.subgroupGrNr[i]
      
      ; check LenType totals
      if total(gc.subgroupLenType[*,i],/int) ne gc.subgroupLen[i] then message,'Fail sgLT'
    endfor
  
  ;endfor ;m
  
  print,'Pass.'
  
end

; cosmoCompareHaloCenters(): compare relative differences between different calcuations of the centers
;                            of halos, and verify group catalog consistency (group-subgroup mappings)

pro cosmoCompareHaloCenters

  compile_opt idl2, hidden, strictarr, strictarrsubs

  res = 128
  run = 'dev.tracer.nograd'
  redshift = 3.0
  
  sP = simParams(res=res,run=run,redshift=redshift)
  
  gc    = loadGroupCat(sP=sP,/skipIDs)
  gcCen = subgroupPosByMostBoundID(sP=sP)
  
  binSize = 1.0
  min = 0.0
  max = 50.0
  
  num = 3500
  
  ; group catalog consistency checks
  print,'nGroups nSubgroups',gc.nGroupsTot,gc.nSubgroupsTot
  ; should pass this
  for i=0,gc.nGroupsTot-1 do begin
    gcInd = gcPriChildInd(gc=gc,haloID=i)
    if (gcInd eq -1) then continue
    grNr = gc.subgroupGrNr[gcInd]
    if (grNr ne i) then stop
  endfor
  
  ; centers
  cen_fof  = gc.groupPos[*,0:num]
  
  cen_sfcm = fltarr(3,num)
  cen_mb   = fltarr(3,num)
  
  w = []
  
  for i=0,num-1 do begin
    gcInd = gcPriChildInd(gc=gc,haloID=i)
    if (gcInd ne -1) then begin
      w = [w,i]
      cen_sfcm[*,i] = gc.subgroupCM[*,gcInd]
      cen_mb[*,i] = gcCen[*,gcInd]
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
