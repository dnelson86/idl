; cosmoVisTry.pro
; cosmological boxes - 2d visualization (unused in papers / experimental)
; dnelson jul.2013

; scatterPlotBox():

pro scatterPlotBox

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  res      = 128
  run      = 'tracer'
  redshift = 0.0
  partType = 'dm' ;tracervel, tracermc, dm

  ; plot config
  zoomFac   = 1     ; only in axes, not along projection direction
  sliceFac  = 1     ; along projection direction (1=whole box, 5=20% width)
  axes      = [0,1] ; x,y
  sliceAxis = 2     ; z

  ; paths and render config
  sP = simParams(res=res,run=run,redshift=redshift)
  h = loadSnapshotHeader(sP=sP)
  
  ; view parameters
  boxSizeImg = [h.boxSize,h.boxSize,h.boxSize]              ;kpc
  boxCen     = [h.boxSize/2.0,h.boxSize/2.0,h.boxSize/2.0]  ;kpc
  
  foreach k,axes do boxSizeImg[k] /= zoomFac
  
  plotFilename = 'scatter.box.'+sP.saveTag+'.'+str(sP.res)+'.'+str(sP.snap)+'.zF'+str(zoomFac)+$
    '.sF'+str(sliceFac)+'.axisX.'+str(axes[0])+'.axisY.'+str(axes[0])+'.'+partType+'.eps'

  ; load positions and color quantity
  pos   = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
    
  ; make slice in projection direction
  sliceMin = boxCen[sliceAxis] - (h.boxSize/2.0)/sliceFac
  sliceMax = boxCen[sliceAxis] + (h.boxSize/2.0)/sliceFac
  
  wSlice = where(pos[sliceAxis,*] ge sliceMin and pos[sliceAxis,*] le sliceMax,countSlice)
  if countSlice eq 0 then message,'Error: Slice empty.'
    
  ; scale color value
  ; todo

  ; plot
  xMinMax = [boxCen[0]-boxSizeImg[0]/2.0,boxCen[0]+boxSizeImg[0]/2.0]
  yMinMax = [boxCen[1]-boxSizeImg[1]/2.0,boxCen[1]+boxSizeImg[1]/2.0]
  
  start_PS, sP.plotPath + plotFilename, xs=8, ys=8
    ;loadColorTable,'helix'
    xrange = [boxCen[axes[0]] - boxSizeImg[axes[0]]*0.5, boxCen[axes[0]] + boxSizeImg[axes[0]]*0.5]
    yrange = [boxCen[axes[1]] - boxSizeImg[axes[1]]*0.5, boxCen[axes[1]] + boxSizeImg[axes[1]]*0.5]
    
    cgPlot,[0],[0],xrange=xrange,yrange=yrange,/xs,/ys
    plots,pos[axes[0],wSlice],pos[axes[1],wSlice],psym=3
  end_PS, pngResize=50, /deletePS

  stop
end

; makeTrackingCutouts(): smooth a halo position in time then make spatial cutouts and save
; the idea was to plot individual gas elements, interpolating all the quantities in time, but 
; dealing with all the new/disappearing cells especially across 4 snaps (for hermite) was impossible

pro makeTrackingCutouts

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  sP = simParams(res=512,run='feedback',redshift=2.0)
  
  sizeFac   = 3.5
  velVecFac = 0.01 ; must match cosmoVisCutout()!
  haloID    = 304  ; z2.304 z2.301 z2.130 z2.64

  ; movie config
  redshiftStart = 2.5
  redshiftEnd   = 2.0
  axes          = [0,1,2] ; rotation about first axis
  nFrames       = 60      ; number of interpolation points along each trajectory
  nFramesOrbit  = 60      ; do one complete orbit at fixed time at the end (0=disable)

  ; calculate frame stepping in Gyr
  timeStart  = redshiftToAgeFlat(redshiftStart) + 1e-4
  timeEnd    = redshiftToAgeFlat(redshiftEnd) - 1e-4
  timeStep   = (timeEnd - timeStart) / (nFrames-1)
  
  frameTimes    = timeStep * findgen(nFrames) + timeStart
  frameTimesRev = reverse(frameTimes)
  
  ; get snapshot times
  cutoutTimes = snapNumToRedshift(sP=sP,/all)
  cutoutTimes = redshiftToAgeFlat(cutoutTimes[where(cutoutTimes ge 0.0)])
  
  ; get list of tracked halo indices
  gcID = ( getMatchedIDs(sPa=sP,sPg=sP,haloID=haloID) ).a
  gcIDList = trackedHaloInds(sP=sP,gcID=gcID)
  
  ; load mergerTreeSubset and get this halo, make smoothed position
  mt = mergerTreeSubset(sP=sP)  
  
  mtInd = where(mt.galcatIDList eq gcID,count)
  
  smoothPos = fltarr(n_elements(cutoutTimes),3)
  smoothPos[mt.minSnap : mt.maxSnap, *] = reverse(smoothHaloPos(mt=mt,hInd=mtInd,sP=sP))
  
  ; sanity checks
  if count eq 0 then message,'Error: Failed to find halo in mtS.'
  if mt.hMinSnap[mtInd] ne -1 then message,'Error: Halo is not tracked all the way.'
  
  if redshiftStart lt min(1/mt.times-1) then message,'Error: redshiftStart beyond snapshot range.'
  if redshiftEnd gt max(1/mt.times-1) then message,'Error: redshiftEnd beyond snapshot range.'
  if redshiftEnd ne sP.redshift then message,'Are you sure?'  
  
  ; decide the initial bracketing
  if min(frameTimes) lt min(cutoutTimes) or max(frameTimes) gt max(cutoutTimes) then $
    message,'Error: Requested frame times out of bounds.'
    
  firstSnap = value_locate(cutoutTimes, frameTimes)
  
  curBracket = [firstSnap[0]-1,firstSnap[0],firstSnap[0]+1,firstSnap[0]+2]
  gcIDCur = gcIDList[curBracket]
  
  ; load the initial cutouts
  print,'[0]    start bracketing, '+str(curBracket[0])+' - '+str(curBracket[3])
  sP.snap = curBracket[0]
  cutout_0 = cosmoVisCutout(sP=sP,gcInd=gcIDCur[0],sizeFac=sizeFac)
  sP.snap = curBracket[1]
  cutout_1 = cosmoVisCutout(sP=sP,gcInd=gcIDCur[1],sizeFac=sizeFac)
  sP.snap = curBracket[2]
  cutout_2 = cosmoVisCutout(sP=sP,gcInd=gcIDCur[2],sizeFac=sizeFac)
  sP.snap = curBracket[3]
  cutout_3 = cosmoVisCutout(sP=sP,gcInd=gcIDCur[3],sizeFac=sizeFac)
  
  ; correct for the smoothed halo position by adding the delta with the unsmoothed pos
  pos_delta = cutout_0.boxCen - smoothPos[curBracket[0],*]
  for j=0,2 do cutout_0.loc_pos[j,*] += pos_delta[j]
  pos_delta = cutout_1.boxCen - smoothPos[curBracket[1],*]
  for j=0,2 do cutout_1.loc_pos[j,*] += pos_delta[j]
  pos_delta = cutout_2.boxCen - smoothPos[curBracket[2],*]
  for j=0,2 do cutout_2.loc_pos[j,*] += pos_delta[j]
  pos_delta = cutout_3.boxCen - smoothPos[curBracket[3],*]
  for j=0,2 do cutout_3.loc_pos[j,*] += pos_delta[j]
  
  ; for hermite we need gas cells to exist at all 4 snaps, do match between extremes
  calcMatch,cutout_0.loc_ids,cutout_3.loc_ids,ind0,ind3,count=countMatch
  calcMatch,cutout_1.loc_ids,cutout_3.loc_ids,ind1,ind3b,count=countMatch2
  calcMatch,cutout_2.loc_ids,cutout_3.loc_ids,ind2,ind3c,count=countMatch3
  
  if countMatch2 ne countMatch or countMatch3 ne countMatch3 then message,'a bit unexpected'
  
  ; take intersection of ind3,ind3b,ind3c
  calcMatch,ind3,ind3b,inda,indb
  ind3 = ind3[inda]
  
  ; TODO
  
  stop
  
  ; grab subsets
  cutout_1_vel  = (cutout_1.loc_pos2[*,ind1] - cutout_1.loc_pos[*,ind1]) / velVecFac
  cutout_2_vel  = (cutout_2.loc_pos2[*,ind2] - cutout_2.loc_pos[*,ind2]) / velVecFac
  cutout_1_pos  = cutout_1.loc_pos[*,ind1]
  cutout_2_pos  = cutout_2.loc_pos[*,ind2]
  cutout_1_temp = cutout_1.loc_temp[ind1]
  cutout_2_temp = cutout_2.loc_temp[ind2]
  
  ; start render frame loop
  for i=0,nFrames-1 do begin
    print,i
    ; at the requested time of this frame, we need the 2 bracketing cutouts lower and higher
    ; so if necessary, shift them back one and load the next
    if firstSnap[i] ne curBracket[1] then begin
      curBracket = [firstSnap[i]-1,firstSnap[i],firstSnap[i]+1,firstSnap[i]+2]
      
      ; truncate at last snapshot if necessary
      if curBracket[3] gt sP.snap then curBracket = curBracket[0:2]
      
      gcIDCur = gcIDList[curBracket]
      print,'['+str(i)+'] progress bracketing, '+str(curBracket[0])+' - '+str(curBracket[-1])
      
      ; shift cutouts and load next
      cutout_0 = cutout_1
      cutout_1 = cutout_2
      cutout_2 = cutout_3
      sP.snap = curBracket[-1]
      cutout_3 = cosmoVisCutout(sP=sP,gcInd=gcIDCur[-1],sizeFac=sizeFac)
      
      ; correct for smoothed halo position
      pos_delta = cutout_3.boxCen - smoothPos[curBracket[-1],*]
      for j=0,2 do cutout_3.loc_pos[j,*] += pos_delta[j]
      
      ; since we are between cutout_1 and cutout_2 (always), which gas IDs match?
      calcMatch,cutout_1.loc_ids,cutout_2.loc_ids,ind1,ind2,count=countMatch
      
      ; those in cutout_1 but not in cutout_2 have disappeared/derefined
      cutout_1_pos = cutout_1.loc_pos[*,ind1]
      cutout_2_pos = cutout_2.loc_pos[*,ind2]
  
      ; those in cutout_2 but not in cutout_1 have appeared/refined from other cells
    endif
    
    ; spline interpolation of all axes coordinates and current temperatures for this subset
    gasPos  = fltarr(3,countMatch)
    gasVel  = fltarr(3,countMatch)
    gasTemp = fltarr(countMatch)

    ; interpolate without using known derivatives
    times = [cutoutTimes[curBracket[1]], cutoutTimes[curBracket[2]]]
    for j=0,countMatch-1 do begin
      gasPos[0,j] = interpol([cutout_1_pos[0,j],cutout_2_pos[0,j]],times,frameTimes[i])
      gasPos[1,j] = interpol([cutout_1_pos[1,j],cutout_2_pos[1,j]],times,frameTimes[i])
      gasPos[2,j] = interpol([cutout_1_pos[2,j],cutout_2_pos[2,j]],times,frameTimes[i])
      gasVel[0,j] = interpol([cutout_1_vel[0,j],cutout_2_vel[0,j]],times,frameTimes[i])
      gasVel[1,j] = interpol([cutout_1_vel[1,j],cutout_2_vel[1,j]],times,frameTimes[i])
      gasVel[2,j] = interpol([cutout_1_vel[2,j],cutout_2_vel[2,j]],times,frameTimes[i])
      gasTemp[j]  = interpol([cutout_1_temp[j],cutout_2_temp[j]],times,frameTimes[i])
    endfor
    
    stop
  
  endfor
  
  stop
end

; scatterMapPastPosComp(): plot past positions of gas segregated by hot/cold (compare arepo/gadget top/bottom)

pro scatterMapPastPosComp

  compile_opt idl2, hidden, strictarr, strictarrsubs

  sPg = simParams(res=512,run='gadget',redshift=2.0)
  sPa = simParams(res=128,run='tracer',redshift=2.0)
  print,'change tracer 128 to 512'
  message,'better test this all, maybe never finished it'
  
  ; OLD512 z2.304 g2342 a2132 -- z2.301 g2289 a2034 -- z2.130 g6369 a5966 axes02 -- z2.64 g5498 a5097
  ; NEW512 z2.304 g2342 a2004 -- 
  ; NEW256 z2.304 g673 a510   -- 
  ; NEW128 z2.304 g217 a150   -- 
  gcIDg = 2342
  gcIDa = 150

  ; config
  sizeFac     = 3.5       ; times rvir
  tempMinMax  = [4.0,7.0] ; log(K)
  coldTempCut = 5.0       ; log(K)
  timeBack    = 500.0     ; Myr
  velVecFac   = 0.01      ; times velocity (km/s) in plotted kpc
  accMode     = 'all'
  axisPair    = [0,1]     ; xy
  
  ; GADGET
  ; ------
  gc = loadGroupCat(sP=sPg,/skipIDs,/verbose)

  saveFilename = sPg.derivPath + 'cutout.' + sPg.savPrefix + str(sPg.res) + '.' + str(sPg.snap) + $
                 '.h' + str(gcIDg) + '.tb' + str(fix(timeBack)) + '.sf' + str(fix(sizeFac*10)) + '.sav'
  
  if file_test(saveFilename) then begin
    restore,saveFilename
  endif else begin     
    ; load ids at starting redshift and make a hot/cold galaxy selection
    print,'Locating...'
    at = accretionTimes(sP=sPg)
    mt = mergerTreeSubset(sP=sPg)
    
    wAm = accModeInds(at=at,accMode=accMode,sP=sPg,/mask)
    
    ; find group members
    gcIndOrig = mergerTreeRepParentIDs(mt=mt,sP=sPg)
    ww = where(gcIndOrig.gal[wAm.gal] eq gcIDg,count)
    
    gcIndOrig = !NULL & at = !NULL & mt = !NULL
    
    ; load max temps, current tvir, tvir at accretion
    accTvir = gcSubsetProp(sP=sPg,select='pri',/accTvir,/accretionTimeSubset,accMode=accMode)
    maxTemp = gcSubsetProp(sP=sPg,select='pri',/maxPastTemp,/accretionTimeSubset,accMode=accMode)
    
    accTvir = accTvir.gal[ww]
    maxTemp = maxTemp.gal[ww]
    ratio   = 10.0^maxTemp / 10.0^accTvir
    
    ; load galcat and make gas ID list
    galIDs = gcSubsetProp(sP=sPg,select='pri',/elemIDs,/accretionTimeSubset,accMode=accMode)
    galIDs = galIDs.gal[ww]
    wAm = !NULL
    
    galIDs = { hot : galIDs[where(ratio ge 1.0)], cold : galIDs[where(ratio lt 1.0)] }
    if n_elements(galIDs.hot) + n_elements(galIDs.cold) ne n_elements(ratio) then message,'error'
    
    ; calculate previous snapshot required and load
    snapTimes = snapNumToRedshift(sP=sPg,/all)
    snapTimes = snapTimes[where(snapTimes ne -1)]
    snapTimes = redshiftToAgeFlat(snapTimes)
    
    curAge = redshiftToAgeFlat(snapNumToRedshift(sP=sPg))
    
    newSnap = value_locate(curAge-snapTimes,timeBack/1000.0)
    
    ; track halo back in time (if possible) to get an earlier center position
    boxCen = trackHaloPosition(sP=sPg,gcID=gcIDg,endSnap=newSnap)
    
    ; calculate bounds
    boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcIDg]] / 10.0) * 10.0
    boxSizeImg = [boxSize,boxSize,boxSize] ; cube   
    
    sPg.snap = newSnap

    ; load ids and match for indices
    print,'Loading...'
    ids = loadSnapshotSubset(sP=sPg,partType='gas',field='ids')
    calcMatch,ids,galIDs.hot,ids_ind_hot,inds2,count=countHot
    calcMatch,ids,galIDs.cold,ids_ind_cold,inds2,count=countCold
    if countHot ne n_elements(galIDs.hot) or countCold ne n_elements(galIDs.cold) then message,'error2'
    ids = !NULL
    
    ; load u,nelec and calculate temperature
    u     = loadSnapshotSubset(sP=sPg,partType='gas',field='u')
    nelec = loadSnapshotSubset(sP=sPg,partType='gas',field='nelec')
    temp  = alog10(convertUtoTemp(u,nelec))
    u     = !NULL
    nelec = !NULL
    temp = { hot : temp[ids_ind_hot], cold : temp[ids_ind_cold] }
  
    ; load gas positions and velocities
    pos = loadSnapshotSubset(sP=sPg,partType='gas',field='pos')
    pos = { hot : pos[*,ids_ind_hot], cold : pos[*,ids_ind_cold] }
    vel = loadSnapshotSubset(sP=sPg,partType='gas',field='vel')
    vel = { hot : vel[*,ids_ind_hot], cold : vel[*,ids_ind_cold] }
    
    ; adjust positions periodic relative
    xDist = pos.hot[0,*] - boxCen[0]
    yDist = pos.hot[1,*] - boxCen[1]
    zDist = pos.hot[2,*] - boxCen[2]
    
    correctPeriodicDistVecs, xDist, sP=sPg
    correctPeriodicDistVecs, yDist, sP=sPg
    correctPeriodicDistVecs, zDist, sP=sPg

    pos.hot[0,*] = xDist & pos.hot[1,*] = yDist & pos.hot[2,*] = zDist
    
    xDist = pos.cold[0,*] - boxCen[0]
    yDist = pos.cold[1,*] - boxCen[1]
    zDist = pos.cold[2,*] - boxCen[2]
    
    correctPeriodicDistVecs, xDist, sP=sPg
    correctPeriodicDistVecs, yDist, sP=sPg
    correctPeriodicDistVecs, zDist, sP=sPg
    
    pos.cold[0,*] = xDist & pos.cold[1,*] = yDist & pos.cold[2,*] = zDist
        
    xDist = !NULL & yDist = !NULL & zDist = !NULL
    
    ; create endpoint for each position point for the velocity vector line
    pos2 = { hot : pos.hot, cold : pos.cold }
    pos2.hot[0,*] = pos.hot[0,*] + vel.hot[0,*]*velVecFac
    pos2.hot[1,*] = pos.hot[1,*] + vel.hot[1,*]*velVecFac
    pos2.hot[2,*] = pos.hot[2,*] + vel.hot[2,*]*velVecFac
    pos2.cold[0,*] = pos.cold[0,*] + vel.cold[0,*]*velVecFac
    pos2.cold[1,*] = pos.cold[1,*] + vel.cold[1,*]*velVecFac
    pos2.cold[2,*] = pos.cold[2,*] + vel.cold[2,*]*velVecFac

    ; save
    save,pos,temp,pos2,sPg,galIDs,gcIDg,sizeFac,boxCen,boxSizeImg,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sPg.derivPath))
  
  endelse

  ; create color index mapping
  colorinds_hot = (temp.hot-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
  colorinds_hot = fix(colorinds_hot + 50.0) > 0 < 255 ;50-255  
  colorinds_cold = (temp.cold-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
  colorinds_cold = fix(colorinds_cold + 50.0) > 0 < 255 ;50-255  
  
  print,'rendering gadget...'
  ; get box center (in terms of specified axes)
  haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcIDg]] ;ckpc
  haloMass = codeMassToLogMsun(gc.subgroupMass[gcIDg])
  
  plotFilename = 'scatter.'+sPg.savPrefix+str(sPg.res)+'.h'+str(gcIDg)+'-'+sPa.savPrefix+str(sPa.res)+$
                 '.h'+str(gcIDa)+'.tb-'+str(fix(timeBack))+'.'+str(sPg.snap)+'.axes'+str(axisPair[0])+str(axisPair[1])+'.eps'

  config = {boxSizeImg:boxSizeImg,plotFilename:plotFilename,haloVirRad:haloVirRad,haloMass:haloMass,$
            axisPair:axisPair,sP:sPg,barMM:tempMinMax,barType:'1temp'}
  
  ; plot
  start_PS, sPg.plotPath + config.plotFilename, xs=8, ys=8
  plotScatterComp,pos.hot,pos2.hot,pos.cold,pos2.cold,colorinds_hot,colorinds_cold,config=config,/top           

  ; AREPO
  ; -----
  gc = loadGroupCat(sP=sPa,/skipIDs,/verbose)
  
  saveFilename = sPa.derivPath + 'cutout.' + sPa.savPrefix + str(sPa.res) + '.' + str(sPa.snap) + $
                 '.h' + str(gcIDa) + '.tb' + str(fix(timeBack)) + '.sf' + str(fix(sizeFac*10)) + '.sav'
  
  if file_test(saveFilename) then begin
    restore,saveFilename
  endif else begin
    ; load ids at starting redshift and make a hot/cold galaxy selection
    print,'Locating...'
    at = accretionTimes(sP=sPa)
    mt = mergerTreeSubset(sP=sPa)
    
    wAm = accModeInds(at=at,accMode=accMode,sP=sPa)
    at = !NULL
    
    ; find group members (and get tracer ID list)
    gcIndOrig = mergerTreeRepParentIDs(mt=mt,sP=sPa,trids_gal=galcat_gal_trids)
    ww = where(gcIndOrig.gal[wAm.gal] eq gcIDa,count)
    
    galIDs = galcat_gal_trids[wAm.gal[ww]]
    
    gcIndOrig = !NULL & mt = !NULL & wAm = !NULL
    
    ; load max temps, current tvir, tvir at accretion
    accTvir = gcSubsetProp(sP=sPa,select='pri',/accTvir,/accretionTimeSubset,accMode=accMode)
    maxTemp = gcSubsetProp(sP=sPa,select='pri',/maxPastTemp,/accretionTimeSubset,accMode=accMode)

    accTvir = accTvir.gal[ww]
    maxTemp = maxTemp.gal[ww]
    ratio   = 10.0^maxTemp / 10.0^accTvir
    
    galIDs = { hot : galIDs[where(ratio ge 1.0)], cold : galIDs[where(ratio lt 1.0)] }
    if n_elements(galIDs.hot) + n_elements(galIDs.cold) ne n_elements(ratio) then message,'error'
    
    ; calculate previous snapshot required and load
    snapTimes = snapNumToRedshift(sP=sPa,/all)
    snapTimes = snapTimes[where(snapTimes ne -1)]
    snapTimes = redshiftToAgeFlat(snapTimes)
    
    curAge = redshiftToAgeFlat(snapNumToRedshift(sP=sPa))
    
    newSnap = value_locate(curAge-snapTimes,timeBack/1000.0)
    
    ; track halo back in time (if possible) to get an earlier center position
    boxCen = trackHaloPosition(sP=sPa,gcID=gcIDa,endSnap=newSnap)
    
    ; calculate bounds
    boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcIDa]] / 10.0) * 10.0
    boxSizeImg = [boxSize,boxSize,boxSize] ; cube   
    
    sPa.snap = newSnap

    ; load ids and match for tracer indices
    print,'Loading...'
    ids = loadSnapshotSubset(sP=sPa,partType='tracerMC',field='tracerIDs')
    calcMatch,ids,galIDs.hot,ids_ind_hot,inds2,count=countHot
    calcMatch,ids,galIDs.cold,ids_ind_cold,inds2,count=countCold
    if countHot ne n_elements(galIDs.hot) or countCold ne n_elements(galIDs.cold) then message,'error2'
    ids = !NULL
    
    ; load parent IDs, gas ids, and convert to parent indices
    parIDs = loadSnapshotSubset(sP=sPa,partType='tracerMC',field='parentIDs')
    parIDs = { hot : parIDs[ids_ind_hot], cold : parIDs[ids_ind_cold] }
    
    ids = loadSnapshotSubset(sP=sPa,partType='gas',field='ids')
    placeMap = getIDIndexMap(ids,minid=minid)
    ids = !NULL
    
    ids_ind_hot  = placeMap[parIDs.hot-minid]  ; replace ids_ind_hot with gas indices
    ids_ind_cold = placeMap[parIDs.cold-minid] ; same
    placeMap = !NULL
  
    ; load u,nelec and parent gas cells calculate temperature
    u     = loadSnapshotSubset(sP=sPa,partType='gas',field='u')
    nelec = loadSnapshotSubset(sP=sPa,partType='gas',field='nelec')
    temp  = alog10(convertUtoTemp(u,nelec))
    u     = !NULL
    nelec = !NULL
    temp = { hot : temp[ids_ind_hot], cold : temp[ids_ind_cold] }
  
    ; load gas positions and velocities
    pos = loadSnapshotSubset(sP=sPa,partType='gas',field='pos')
    pos = { hot : pos[*,ids_ind_hot], cold : pos[*,ids_ind_cold] }
    vel = loadSnapshotSubset(sP=sPa,partType='gas',field='vel')
    vel = { hot : vel[*,ids_ind_hot], cold : vel[*,ids_ind_cold] }
    
    ; adjust positions periodic relative
    xDist = pos.hot[0,*] - boxCen[0]
    yDist = pos.hot[1,*] - boxCen[1]
    zDist = pos.hot[2,*] - boxCen[2]
    
    correctPeriodicDistVecs, xDist, sP=sPa
    correctPeriodicDistVecs, yDist, sP=sPa
    correctPeriodicDistVecs, zDist, sP=sPa

    pos.hot[0,*] = xDist & pos.hot[1,*] = yDist & pos.hot[2,*] = zDist
    
    xDist = pos.cold[0,*] - boxCen[0]
    yDist = pos.cold[1,*] - boxCen[1]
    zDist = pos.cold[2,*] - boxCen[2]
    
    correctPeriodicDistVecs, xDist, sP=sPa
    correctPeriodicDistVecs, yDist, sP=sPa
    correctPeriodicDistVecs, zDist, sP=sPa
    
    pos.cold[0,*] = xDist & pos.cold[1,*] = yDist & pos.cold[2,*] = zDist
        
    xDist = !NULL & yDist = !NULL & zDist = !NULL
    
    ; create endpoint for each position point for the velocity vector line
    pos2 = { hot : pos.hot, cold : pos.cold }
    pos2.hot[0,*] = pos.hot[0,*] + vel.hot[0,*]*velVecFac
    pos2.hot[1,*] = pos.hot[1,*] + vel.hot[1,*]*velVecFac
    pos2.hot[2,*] = pos.hot[2,*] + vel.hot[2,*]*velVecFac
    pos2.cold[0,*] = pos.cold[0,*] + vel.cold[0,*]*velVecFac
    pos2.cold[1,*] = pos.cold[1,*] + vel.cold[1,*]*velVecFac
    pos2.cold[2,*] = pos.cold[2,*] + vel.cold[2,*]*velVecFac

    ; save
    save,pos,temp,pos2,sPa,galIDs,gcIDa,sizeFac,boxCen,boxSizeImg,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sPa.derivPath))
  
  endelse

  ; create color index mapping
  colorinds_hot = (temp.hot-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
  colorinds_hot = fix(colorinds_hot + 50.0) > 0 < 255 ;50-255  
  colorinds_cold = (temp.cold-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
  colorinds_cold = fix(colorinds_cold + 50.0) > 0 < 255 ;50-255  

  print,'rendering arepo...'
  ; get box center (in terms of specified axes)
  haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcIDa]] ;ckpc
  haloMass = codeMassToLogMsun(gc.subgroupMass[gcIDa])

  config = {boxSizeImg:boxSizeImg,plotFilename:plotFilename,haloVirRad:haloVirRad,haloMass:haloMass,$
            axisPair:axisPair,sP:sPa,barMM:tempMinMax,barType:'1temp'}
  
  ; plot
  plotScatterComp,pos.hot,pos2.hot,pos.cold,pos2.cold,colorinds_hot,colorinds_cold,config=config,/bottom
  end_PS, pngResize=60
  
  stop
end

; scatterMapHalosGasDM: plot temperature colored scatter plots with velocity vectors on boxes centered on halos
;                       gas on left and DM on right (usually more zoomed out)

pro scatterMapHalosGasDM;, sP=sP, gcIDs=gcIDs

  sP = simParams(res=512,run='gadgetold',redshift=2.0)
  gcIDs = [2342] ;z2.304 g2342 a2132 ;z2.301 g2289 a2034

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  sizeFac         = 10.5      ; times rvir
  tempMinMax      = [4.0,7.0] ; log(K)
  velDispMaxDisp  = 301.0     ; km/s
  velVecFac       = 0.01      ; times velocity (km/s) in plotted kpc
  
  axes = list([0,1],[0,2],[1,2]) ;xy,xz,yz
  
  ; target list
  h     = loadSnapshotHeader(sP=sP)
  gc    = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sP) 

  ; load u,nelec and calculate temperature
  u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
  nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
  temp  = convertUtoTemp(u,nelec,/log)
  u     = !NULL
  nelec = !NULL

  ; load gas positions and velocities
  pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  vel = loadSnapshotSubset(sP=sP,partType='gas',field='vel')
  
  ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
  print,'shuffling gas...'
  iseed = 424242L
  sort_inds = sort(randomu(iseed,n_elements(temp)))
  
  temp = temp[sort_inds]
  pos  = pos[*,sort_inds]
  vel  = vel[*,sort_inds]
  
  sort_inds = !NULL

  ; find the DM positions in the hsmldir for the veldisps (have to load all particle IDs)
  veldisp = loadHsmlDir(sP=sP,partType='dm',/readVelDisp,/verbose)
  
  ; load dm positions and velocities
  print,'loading dm...'
  pos_dm  = loadSnapshotSubset(sP=sP,partType='dm',field='pos')
  vel_dm  = loadSnapshotSubset(sP=sP,partType='dm',field='vel')
  
  ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
  print,'shuffling dm...'
  iseed = 434343L
  sort_inds = sort(randomu(iseed,n_elements(veldisp)))
  
  veldisp = veldisp[sort_inds]
  pos_dm  = pos_dm[*,sort_inds]
  vel_dm  = vel_dm[*,sort_inds]
  
  sort_inds = !NULL
  
  print,'rendering...'
  ; loop over all requested halos and image
  foreach gcID, gcIDs do begin
  
    ; get subhalo position and size of imaging box
    boxCen     = sgcen[*,gcID]
    boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcID]] / 10.0) * 10.0
    boxSizeImg = [boxSize,boxSize,boxSize] ; cube
  
    ; make conservative cutout greater than boxsize accounting for periodic (do cube not sphere)
    xDist = pos[0,*] - boxCen[0]
    yDist = pos[1,*] - boxCen[1]
    zDist = pos[2,*] - boxCen[2]
    
    correctPeriodicDistVecs, xDist, sP=sP
    correctPeriodicDistVecs, yDist, sP=sP
    correctPeriodicDistVecs, zDist, sP=sP
  
    ; local (cube) cutout
    wCut = where(abs(xDist) le 0.5*boxSize and abs(yDist) le 0.5*boxSize and $
                 abs(zDist) le 0.5*boxSize,nCutout)
    
    loc_temp = temp[wCut]
    loc_pos  = fltarr(3,nCutout)
    loc_pos[0,*] = xDist[wCut] ; delta
    loc_pos[1,*] = yDist[wCut]
    loc_pos[2,*] = zDist[wCut]
    
    xDist = !NULL
    yDist = !NULL
    zDist = !NULL
    
    loc_vel = vel[*,wCut]
    
    ; create endpoint for each position point for the velocity vector line
    loc_pos2 = fltarr(3,nCutout)
    loc_pos2[0,*] = loc_pos[0,*] + loc_vel[0,*]*velVecFac
    loc_pos2[1,*] = loc_pos[1,*] + loc_vel[1,*]*velVecFac
    loc_pos2[2,*] = loc_pos[2,*] + loc_vel[2,*]*velVecFac
    loc_vel = !NULL
  
    ; create color index mapping
    colorinds = (loc_temp-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
    colorinds = fix(colorinds + 50.0) > 0 < 255 ;50-255
    loc_temp = !NULL
  
    ; DM: make cutout
    xDist = pos_dm[0,*] - boxCen[0]
    yDist = pos_dm[1,*] - boxCen[1]
    zDist = pos_dm[2,*] - boxCen[2]
    
    correctPeriodicDistVecs, xDist, sP=sP
    correctPeriodicDistVecs, yDist, sP=sP
    correctPeriodicDistVecs, zDist, sP=sP
    
    ; local (cube) cutout
    wCut = where(abs(xDist) le 0.5*boxSize and abs(yDist) le 0.5*boxSize and $
                 abs(zDist) le 0.5*boxSize and $
                 veldisp lt velDispMaxDisp,nCutoutDM)
    
    loc_veldisp = veldisp[wCut]
    loc_pos_dm  = fltarr(3,nCutoutDM)
    
    loc_pos_dm[0,*] = xDist[wCut] ; delta
    loc_pos_dm[1,*] = yDist[wCut]
    loc_pos_dm[2,*] = zDist[wCut]
    
    xDist = !NULL
    yDist = !NULL
    zDist = !NULL
    
    loc_vel_dm = vel_dm[*,wCut]
    
    ; create endpoint for each position point for the velocity vector line
    loc_pos2_dm = fltarr(3,nCutoutDM)
    loc_pos2_dm[0,*] = loc_pos_dm[0,*] + loc_vel_dm[0,*]*velVecFac
    loc_pos2_dm[1,*] = loc_pos_dm[1,*] + loc_vel_dm[1,*]*velVecFac
    loc_pos2_dm[2,*] = loc_pos_dm[2,*] + loc_vel_dm[2,*]*velVecFac
    loc_vel_dm = !NULL
  
    ; create color index mapping
    veldispMM = [0.0,floor(max(loc_veldisp)/100.0)*100.0 > 100 < 400]
    colorinds_dm = (loc_veldisp-veldispMM[0])*205.0 / (veldispMM[1]-veldispMM[0]) ;0-205
    colorinds_dm = fix(colorinds_dm + 50.0) > 0 < 255 ;50-255
    loc_veldisp = !NULL
    
    ; make a plot for each requested projection direction
    foreach axisPair, axes do begin
           
      ; get box center (in terms of specified axes)
      boxCenImg  = [sgcen[axisPair[0],gcID],sgcen[axisPair[1],gcID],sgcen[3-axisPair[0]-axisPair[1],gcID]]
      haloMass = codeMassToLogMsun(gc.subgroupMass[gcID])
      haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcID]] ;ckpc
      
      ; fix halo mass if we're using old (x2 bug) catalogs
      if sP.run eq 'gadgetold' then haloMass = codeMassToLogMsun(0.5*gc.subgroupMass[gcID])
      
      print,'['+string(gcID,format='(i4)')+'] Mapping ['+str(axisPair[0])+' '+$
            str(axisPair[1])+'] with '+str(boxSize[0])+$
            ' kpc box around subhalo center ['+str(boxCen[0])+' '+str(boxCen[1])+' '+str(boxCen[2])+']'

      plotFilename = 'gasdm.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+$
                     '.axes'+str(axisPair[0])+str(axisPair[1])+'.eps'

      config = {boxSizeImg:boxSizeImg,plotFilename:plotFilename,haloVirRad:haloVirRad,haloMass:haloMass,$
                axisPair:axisPair,sP:sP,barMM_left:tempMinMax,barMM_right:veldispMM,barType:'2tempvdisp'}
      
      ; plot
      plotScatterComp,loc_pos,loc_pos2,loc_pos_dm,loc_pos2_dm,colorinds,colorinds_dm,config=config

    endforeach ;axisPair

  endforeach ;gcIDs

end

; scatterMapHalosDM: plot veldisp colored scatter plots with velocity vectors on boxes centered on halos
;                    all DM on left and DM in some overdensity range on right

pro scatterMapHalosDM, sP=sP, gcIDs=gcIDs

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  sizeFac         = 3.5      ; times rvir
  cutFac          = 1.0       ; times boxSize
  overdensMinMax  = [0.0,1.0] ; log(rho/mean rho)
  velVecFac       = 0.01      ; times velocity (km/s) in plotted kpc
  
  axes = list([0,1],[0,2],[1,2]) ;xy,xz,yz
  
  ; target list
  h     = loadSnapshotHeader(sP=sP)
  gc    = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sP) 
  
  ; load dm positions and velocities
  print,'loading dm...'
  pos_dm  = loadSnapshotSubset(sP=sP,partType='dm',field='pos')
  vel_dm  = loadSnapshotSubset(sP=sP,partType='dm',field='vel')
  
  ; find the DM positions in the hsmldir for the veldisps (have to load all particle IDs)
  dens_dm = loadHsmlDir(sP=sP,partType='dm',/readDens,/verbose)
  
  ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
  print,'shuffling dm...'
  iseed = 434343L
  sort_inds = sort(randomu(iseed,n_elements(dens_dm)))
  
  dens_dm = dens_dm[sort_inds]
  pos_dm  = pos_dm[*,sort_inds]
  vel_dm  = vel_dm[*,sort_inds]
  
  sort_inds = !NULL
  
  ; convert densities to log(overdensities)
  meanDensityBox = h.nPartTot[partTypeNum('dm')] * h.masstable[partTypeNum('dm')] / (sP.boxSize)^3.0
  dens_dm = alog10( dens_dm / meanDensityBox )
  
  print,'rendering...'
  ; loop over all requested halos and image
  foreach gcID, gcIDs do begin
  
    ; get subhalo position and size of imaging box
    boxCen     = sgcen[*,gcID]
    boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcID]] / 10.0) * 10.0
    boxSizeImg = [boxSize,boxSize,boxSize] ; cube
  
    ; DM: make cutout
    xDist = pos_dm[0,*] - boxCen[0]
    yDist = pos_dm[1,*] - boxCen[1]
    zDist = pos_dm[2,*] - boxCen[2]
    
    correctPeriodicDistVecs, xDist, sP=sP
    correctPeriodicDistVecs, yDist, sP=sP
    correctPeriodicDistVecs, zDist, sP=sP
    
    ; local (cube) cutout
    wCut = where(abs(xDist) le 0.5*cutFac*boxSize and abs(yDist) le 0.5*cutFac*boxSize and $
                 abs(zDist) le 0.5*cutFac*boxSize,nCutoutDM)
    
    loc_dens_dm = dens_dm[wCut]
    loc_pos_dm  = fltarr(3,nCutoutDM)
    
    loc_pos_dm[0,*] = xDist[wCut] ; delta
    loc_pos_dm[1,*] = yDist[wCut]
    loc_pos_dm[2,*] = zDist[wCut]
    
    xDist = !NULL
    yDist = !NULL
    zDist = !NULL
    
    loc_vel_dm = vel_dm[*,wCut]
    
    ; create endpoint for each position point for the velocity vector line
    loc_pos2_dm = fltarr(3,nCutoutDM)
    loc_pos2_dm[0,*] = loc_pos_dm[0,*] + loc_vel_dm[0,*]*velVecFac
    loc_pos2_dm[1,*] = loc_pos_dm[1,*] + loc_vel_dm[1,*]*velVecFac
    loc_pos2_dm[2,*] = loc_pos_dm[2,*] + loc_vel_dm[2,*]*velVecFac
    loc_vel_dm = !NULL
  
    ; create color index mapping
    overdensMM = [ceil(min(loc_dens_dm)) > (-2.0),floor(max(loc_dens_dm)) > 1.0 < 8.0]
    colorinds_dm = (loc_dens_dm-overdensMM[0])*205.0 / (overdensMM[1]-overdensMM[0]) ;0-205
    colorinds_dm = fix(colorinds_dm + 50.0) > 0 < 255 ;50-255

    ; OVERDENSE DM: make cutout
    wOD = where(loc_dens_dm ge overDensMinMax[0] and loc_dens_dm lt overDensMinMax[1],nOverDens)
    print,nCutoutDM,nOverDens
    
    loc_pos_od   = loc_pos_dm[*,wOD]
    loc_pos2_od  = loc_pos2_dm[*,wOD]
    
    ; create restricted color index mapping
    colorinds_od = (loc_dens_dm[wOD]-overdensMinMax[0])*205.0 / (overdensMinMax[1]-overdensMinMax[0]) ;0-205
    colorinds_od = fix(colorinds_dm + 50.0) > 0 < 255 ;50-255
  
    ; make a plot for each requested projection direction
    foreach axisPair, axes do begin
           
      ; get box center (in terms of specified axes)
      boxCenImg  = [sgcen[axisPair[0],gcID],sgcen[axisPair[1],gcID],sgcen[3-axisPair[0]-axisPair[1],gcID]]
      haloMass = codeMassToLogMsun(gc.subgroupMass[gcID])
      haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcID]] ;ckpc
      
      print,'['+string(gcID,format='(i4)')+'] Mapping ['+str(axisPair[0])+' '+$
            str(axisPair[1])+'] with '+str(boxSize[0])+$
            ' kpc box around subhalo center ['+str(boxCen[0])+' '+str(boxCen[1])+' '+str(boxCen[2])+']'

      plotFilename = 'gasdm.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+$
                     '.axes'+str(axisPair[0])+str(axisPair[1])+'.eps'

      config = {boxSizeImg:boxSizeImg,plotFilename:plotFilename,haloVirRad:haloVirRad,haloMass:haloMass,$
                axisPair:axisPair,sP:sP,barMM_left:overdensMM,barMM_right:overDensMinMax,barType:'2overdens'}
      
      ; plot
      plotScatterComp,loc_pos_dm,loc_pos2_dm,loc_pos_od,loc_pos2_od,colorinds_dm,colorinds_od,config=config
stop
    endforeach ;axisPair

  endforeach ;gcIDs

end
