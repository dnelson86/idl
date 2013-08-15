; cosmoVis.pro
; cosmological boxes - animated 2d visualizations, usually tracking halos in time
; dnelson jul.2013

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

; makeSubboxFrames() : generate frames for a movie based on a subbox (either whole subbox or 
;   following one halo through time within that subbox) (single run or comparison of two runs)

pro makeSubboxFrames

  ; config
  zStart = 6.0
  zEnd   = 2.0
  hInd   = -1 ; negative one to show whole subbox
  
  sP1 = simParams(res=512,run='tracer',redshift=zStart)
  sP2 = simParams(res=512,run='feedback',redshift=zStart)
  
  ; load all subbox times, find ranges within [zStart,zEnd]
  
  
  nFrames = 0
  
  ; if subbox times don't match, take the number of frames as the larger, decide for each frame
  ; the corresponding subbox number for the lesser set (duplicate frames) (use value_locate)
  snaps = { run1 : intarr(nFrames), run2 : intarr(nFrames) }
  
  ; TODO
  
  
  ; if following one halo, calculate tracked position vs time from main snapshots, interpolate to subbox times
  
  ; loop over all frames
  for fN=0,nFrames-1 do begin
    sP1.snap = snaps.run1[fN]
    sP2.snap = snaps.run2[fN]
    
    ; load
    x = loadSnapshotSubset(sP=sP1,partType='gas',field='pos',/subBox)
  
    ; calculate colors
    
    ; render (call existing function)
  
  endfor
  
  stop
  
end

; makeProgressionPanels()

pro makeProgressionPanels

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  sP  = simParams(res=256,run='feedback',redshift=0.0)
  
  ; best match in sP2 is found for each sP snap, comment out sP2 to not make a comparison
  sP2 = simParams(res=256,run='tracer',redshift=0.0)

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
  colorField = 'temp'     & fieldMinMax  = [4.0,7.0] & secondCutVal = 5.0
  ;colorField = 'entropy'  & fieldMinMax  = [5.0,9.0] & secondCutVal = 7.5 ; log(CGS)
  ;colorField = 'metal'     & fieldMinMax = [-4.0,-1.0] & secondCutVal = -2.5 ; log(Z/Zsun)
  ;colorField = 'vrad'      & fieldMinMax = [-400,400] & secondCutVal = -200.0 ; km/s
  ;colorField = 'vradnorm'  & fieldMinMax = [-3.0,3.0] & secondCutVal = -1.0 ; vrad/v200
  
  subtitles = ['all gas','cold gas'] ;['low entr'], ['cold gas'], ['inflow']
  
  ; get snapshot times
  cutoutTimes = snapNumToRedshift(sP=sP,/all)
  minSnap = min(where(cutoutTimes le redshiftEnd,numSnaps))
  maxSnap = min( [sP.snap, max(where(cutoutTimes ge redshiftStart,wCheck))] )
  if numSnaps eq 0 then message,'error'
  if wCheck eq 0 then message,'error'
  numSnaps = sP.snap - minSnap + 1
  
  ; get list of tracked halo indices
  gcID     = getMatchedIDs(simParams=sP,haloID=haloID)
  gcIDList = trackedHaloInds(sP=sP,gcID=gcID,minZ=redshiftStart,maxZ=redshiftEnd)

  sP2tag = ''
  if n_elements(sP2) gt 0 then begin
    gcID2     = getMatchedIDs(simParams=sP2,haloID=haloID)
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
