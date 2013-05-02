; cosmoVis.pro
; cosmological boxes - animated 2d visualizations, usually tracking halos in time
; dnelson apr.2013

; trackedHaloInds():

function trackedHaloInds, sP=sP, gcID=gcID, maxZ=maxZ, minZ=minZ

  gcIDList = lonarr(sP.groupCatRange[1]+1) - 1 ; progressing through snapshots
  gcIDList[sP.snap] = gcID
  
  origSnap = sP.snap
  
  ; go backwards in time
  for m=sP.snap,sP.groupCatRange[0],-1 do begin
    sP.snap = m
    
    ; past maximum redshift requested?
    if snapNumToRedshift(sP=sP) gt maxZ+0.01 then break
    
    ; load most massive progenitor and move index
    if m gt sP.groupCatRange[0] then begin
      Parent = mergerTree(sP=sP)
      if Parent[gcIDList[m]] eq -1 then message,'not tracked this far'
      gcIDList[m-1] = Parent[gcIDList[m]]
    endif
    
  endfor
  
  sP.snap = origSnap
  
  ; go forwards in time
  for m=sP.snap+1,sP.groupCatRange[1],1 do begin
    sP.snap = m
    
    ; out of snapshots or passed minimum redshift requested?
    if (getSnapFilelist(sP.simPath,snapNum=sP.snap,/checkExists))[0] eq '-1' then break
    if snapNumToRedshift(sP=sP) lt minZ-0.01 then break
    
    ; load most massive progenitor and find child, move index
    Parent = mergerTree(sP=sP)
    w = where(Parent eq gcIDList[m-1],count)
    if count eq 0 then message,'not tracked this far'
    gcIDList[m] = w[0]
  endfor
  
  sP.snap = origSnap
  
  return, gcIDList

end

; makeTrackingCutouts(): smooth a halo position in time then make spatial cutouts and save

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
  mt = mergerTreeAdaptiveSubset(sP=sP)  
  
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

; makeProgressionPanels()

pro makeProgressionPanels

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  sP  = simParams(res=512,run='feedback',redshift=2.0)
  
  ; best match in sP2 is found for each sP snap, comment out sP2 to not make a comparison
  sP2 = simParams(res=512,run='tracer',redshift=2.0)

  haloID        = 304 ; z2.314 z2.304 z2.301 z2.130 z2.64
  redshiftStart = 2.0 ; when to start (can differ from sP.redshift, which is where the halo is targeted)
  redshiftEnd   = 4.0 ; how far back to go
  
  singleColorScale = 0 ; 1=use same color scale for right panel, 0=rescale
  secondGt         = 0 ; 1=show greater than cut, 0=show less than cut
  axes             = [0,1]
  nbottom          = 50
  sizeFac          = 3.5 ; times rvir
  ctName           = 'helix'
  
  ; use which field and minmax for color mapping? cut value for right panel?
  ;colorField = 'temp'     & fieldMinMax  = [4.0,7.0] & secondCutVal = 5.0
  ;colorField = 'entropy'  & fieldMinMax  = [5.0,9.0] & secondCutVal = 7.5 ; log(CGS)
  ;colorField = 'metal'     & fieldMinMax = [-4.0,-1.0] & secondCutVal = -2.5 ; log(Z/Zsun)
  ;colorField = 'vrad'      & fieldMinMax = [-400,400] & secondCutVal = -200.0 ; km/s
  colorField = 'vradnorm'  & fieldMinMax = [-3.0,3.0] & secondCutVal = -1.0 ; vrad/v200
  
  subtitles = ['all gas','inflow'] ;['low entr'], ['cold gas']
  
  ; get snapshot times
  cutoutTimes = snapNumToRedshift(sP=sP,/all)
  minSnap = min(where(cutoutTimes le redshiftEnd,numSnaps))
  maxSnap = min( [sP.snap, max(where(cutoutTimes ge redshiftStart,wCheck))] )
  if numSnaps eq 0 then message,'error'
  if wCheck eq 0 then message,'error'
  numSnaps = sP.snap - minSnap + 1
  
  ; get list of tracked halo indices
  gcID     = ( getMatchedIDs(sPa=sP,sPg=sP,haloID=haloID) ).a
  gcIDList = trackedHaloInds(sP=sP,gcID=gcID,minZ=redshiftStart,maxZ=redshiftEnd)

  sP2tag = ''
  if n_elements(sP2) gt 0 then begin
    gcID2     = ( getMatchedIDs(sPa=sP2,sPg=sP2,haloID=haloID) ).a
    gcIDList2 = trackedHaloInds(sP=sP2,gcID=gcID2,minZ=redshiftStart,maxZ=redshiftEnd)
    sP2tag    = '.'+sP2.saveTag+'.'+str(sP2.res)
  endif

  ; start plot
  pFilename = 'scatter.'+sP.saveTag+'.'+str(sP.res)+'.'+str(sP.snap)+'-'+str(minSnap)+$
              '.h'+str(haloID)+'.axes'+str(axes[0])+str(axes[1])+'-'+$
              colorField+'-'+str(secondGt)+'sCS'+str(singleColorScale)+sP2tag+'.eps'
              
  pWidth  = 8.0 * (2 * n_elements(sP2)) ; change 8 to 16 if comparing runs
  pHeight = 4.0 * numSnaps
  
  start_PS, sP.plotPath + pFilename, xs=pWidth, ys=pHeight
                     
  ; loop over snapshots
  for m=sP.snap,minSnap,-1 do begin
    sP.snap = m & print,m,gcIDList[sP.snap]
    sP.redshift = snapNumToRedshift(sP=sP)
    
    ; load cutouts
    cutout = cosmoVisCutout(sP=sP,gcInd=gcIDList[sP.snap],sizeFac=sizeFac)
    
    config = {boxSizeImg:cutout.boxSizeImg,plotFilename:'',haloVirRad:cutout.haloVirRad,$
              haloMass:cutout.haloMass,axisPair:axes,sP:sP,singleColorScale:singleColorScale,$
              colorField:colorField,fieldMinMax:fieldMinMax,secondCutVal:secondCutVal,$
              secondGt:secondGt,nbottom:nbottom,barMM:fieldMinMax,ctName:ctName,barType:'2bar'}
    
    if m eq minSnap then config = mod_struct(config, 'subtitle', subtitles)
    
    sub = cosmoVisCutoutSub(cutout=cutout,config=config)

    ; making a comparison to a sP2 halo?
    if n_elements(sP2) gt 0 then begin
      ; find closest snapshot
      sP2.snap = redshiftToSnapNum(sP.redshift,sP=sP2)
      redshift_err = sP.redshift-snapNumToRedshift(sP=sP2)
      print,sP2.snap,gcIDList2[sP2.snap],redshift_err
      
      ; load cutouts
      cutout2 = cosmoVisCutout(sP=sP2,gcInd=gcIDList2[sP2.snap],sizeFac=sizeFac)
      config = mod_struct(config, 'sP2', sP2)
      sub2 = cosmoVisCutoutSub(cutout=cutout2,config=config)
      
      ; rearrange subs, want the two all gas together, and the two cold gas together
      sub_temp = sub
      
      sub = {cinds_left:sub.cinds_left,pos_left:sub.pos_left,pos_left2:sub.pos_left2,$
             cinds_right:sub2.cinds_left,pos_right:sub2.pos_left,pos_right2:sub2.pos_left2}
      sub2 = {cinds_left:sub_temp.cinds_right,pos_left:sub_temp.pos_right,pos_left2:sub_temp.pos_right2,$
              cinds_right:sub2.cinds_right,pos_right:sub2.pos_right,pos_right2:sub2.pos_right2}
    endif    
    
    ; add row to plot
    plotScatterComp,sub=sub,config=config,row=[m-minSnap,numSnaps],first=(n_elements(sP2) gt 0)
    
    if n_elements(sP2) gt 0 then $
      plotScatterComp,sub=sub2,config=mod_struct(config,'subtitle',/delete),$
        row=[m-minSnap,numSnaps],/second
  
  endfor
  
  ; finish plot
  end_PS, pngResize=40
  stop
end
