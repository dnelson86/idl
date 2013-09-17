; cosmoVis.pro
; cosmological boxes - animated 2d visualizations, usually tracking halos in time
; dnelson sep.2013

; trackedHaloInds(): starting from a halo at one time, return a list of tracked halo indices
;   both earlier and/or later using the mergerTree

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

; subBoxVisCutout(): make a spatial cutout around a halo within a subbox, or whole subbox

function subBoxVisCutout, sP=sP, gcInd=gcInd, sizeFac=sizeFac

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  velVecFac = 0.01 ; times velocity (km/s) in plotted kpc
  hsmlFac   = 1.75 ; increase arepo 'hsml' to decrease visualization noise  

  if ~keyword_set(sP) or n_elements(gcInd) eq 0 then message,'Error'
  
  ; check existence of requested saves if more than one halo
  saveFilename = sP.derivPath + 'cutouts/cutout.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.h' + str(gcInd) + '_subbox0.sav'
        
  ; if single halo requested and save exists, load it
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif
    
  ; proceed with load

  ; have metallicities to include in cutout?
  metalsFlag = 0
  if snapshotFieldExists(sP=sP,field='Metallicity',/subBox) or $
     snapshotFieldExists(sP=sP,field='GFM_Metallicity',/subBox) then metalsFlag = 1
  
  ; load u,nelec,dens and calculate temperature,entropy
  u     = loadSnapshotSubset(sP=sP,partType='gas',field='u',/subBox)
  nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec',/subBox)
  temp  = alog10(convertUtoTemp(u,nelec))
  nelec = !NULL
  dens  = loadSnapshotSubset(sP=sP,partType='gas',field='dens',/subBox)
  ent   = calcEntropyCGS(u,dens,/log,sP=sP)
  u     = !NULL
  
  if metalsFlag then begin
    metal = loadSnapshotSubset(sP=sP,partType='gas',field='metallicity',/subBox)
    ; convert to log(metallicity) for positive values, otherwise set GFM_MIN_METAL = -20 (log)
    w = where(metal gt 0.0,count,comp=wc,ncomp=ncomp)
    if count gt 0 then metal[w] = alog10(metal[w])
    if ncomp gt 0 then metal[wc] = -20.0
  endif
    
  ; load HSMLs or volumes (convert to sizes)
  if sP.trMCPerCell eq 0 then begin
    hsml = loadSnapshotSubset(sP=sP,partType='gas',field='hsml',/subBox)
    hsml = 1.0 * temporary(hsml); increase hsml to decrease visualization noise
  endif else begin
    hsml = loadSnapshotSubset(sP=sP,partType='gas',field='vol',/subBox)
    hsml = (temporary(hsml) * 3.0 / (4*!pi))^(1.0/3.0) ;cellrad [ckpc]
    hsml = hsmlFac * temporary(hsml) ; increase hsml to decrease visualization noise
  endelse  
  
  ; load gas ids, positions, velocities and masses
  ids  = loadSnapshotSubset(sP=sP,partType='gas',field='ids',/subBox)
  pos  = loadSnapshotSubset(sP=sP,partType='gas',field='pos',/subBox)
  vel  = loadSnapshotSubset(sP=sP,partType='gas',field='vel',/subBox)
  mass = loadSnapshotSubset(sP=sP,partType='gas',field='mass',/subBox)
  
  ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
  iseed = 424242L
  sort_inds = sort(randomu(iseed,n_elements(temp)))
  
  ids   = ids[sort_inds]
  temp  = temp[sort_inds]
  ent   = ent[sort_inds]
  dens  = dens[sort_inds]
  hsml  = hsml[sort_inds]
  mass  = mass[sort_inds]
  pos   = pos[*,sort_inds]
  vel   = vel[*,sort_inds]
  
  if metalsFlag then metal = metal[sort_inds]
  
  sort_inds = !NULL
  
  if gcInd ne -1 then begin
    ; cutout around individual halo
    gc    = loadGroupCat(sP=sP,/skipIDs,/verbose)
    sgcen = subgroupPosByMostBoundID(sP=sP)
  
    ; halo properties
    haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcInd]] ;ckpc
    haloMass = codeMassToLogMsun(gc.subgroupMass[gcInd])
    haloM200 = codeMassToLogMsun(gc.group_m_crit200[gc.subgroupGrNr[gcInd]])
    haloV200 = sqrt(units.G * gc.subgroupMass[gcInd] / haloVirRad )
    
    ; get subhalo position and size of imaging box
    boxCen     = sgcen[*,gcInd]
    boxSize    = replicate( ceil(sizeFac * haloVirRad / 10.0) * 10.0, 3)
    sgcen_vel  = gc.subgroupVel[*,gcInd]
  endif else begin
    ; whole subbox
    haloVirRad = -1
    haloMass   = -1
    haloM200   = -1
    haloV200   = -1
    
    boxCen     = sP.subboxCen
    boxSize    = sP.subboxSize
    sgcen_vel  = [0.0,0.0,0.0]
  endelse
  
  ; make conservative cutout greater than boxsize accounting for periodic (do cube not sphere)
  xDist = reform( pos[0,*] - boxCen[0] )
  yDist = reform( pos[1,*] - boxCen[1] )
  zDist = reform( pos[2,*] - boxCen[2] )
    
  correctPeriodicDistVecs, xDist, sP=sP
  correctPeriodicDistVecs, yDist, sP=sP
  correctPeriodicDistVecs, zDist, sP=sP
    
  ; local (cube) cutout
  wCut = where(abs(xDist) le 0.5*boxSize[0] and abs(yDist) le 0.5*boxSize[1] and $
               abs(zDist) le 0.5*boxSize[2],nCutout)
                    
  ; take selection of fields
  loc_ids   = ids[wCut]
  loc_temp  = temp[wCut]
  loc_ent   = ent[wCut]
  loc_dens  = dens[wCut]
  loc_hsml  = hsml[wCut]
  loc_mass  = mass[wCut]
  loc_pos   = fltarr(3,nCutout)
  loc_pos[0,*] = xDist[wCut] ; delta
  loc_pos[1,*] = yDist[wCut]
  loc_pos[2,*] = zDist[wCut]
    
  xDist = !NULL
  yDist = !NULL
  zDist = !NULL
    
  if metalsFlag then loc_metal = metal[wCut]
  if ~metalsFlag then loc_metal = -1    
    
  loc_vel = vel[*,wCut]
    
  ; only available for halo(cat) selection
  loc_dynTime = -1
  loc_coolTime = -1
    
  ; calculate norm of radial velocity vector
  rad = reform(loc_pos[0,*]^2.0 + loc_pos[1,*]^2.0 + loc_pos[2,*]^2.0)

  loc_vrad = reform( (loc_vel[0,*]-sgcen_vel[0]) * loc_pos[0,*] + $
                     (loc_vel[1,*]-sgcen_vel[1]) * loc_pos[1,*] + $
                     (loc_vel[2,*]-sgcen_vel[2]) * loc_pos[2,*]) / sqrt(rad)
                       
  rad = !NULL
    
  ; create endpoint for each position point for the velocity vector line
  loc_pos2 = fltarr(3,nCutout)
  loc_pos2[0,*] = loc_pos[0,*] + loc_vel[0,*]*velVecFac
  loc_pos2[1,*] = loc_pos[1,*] + loc_vel[1,*]*velVecFac
  loc_pos2[2,*] = loc_pos[2,*] + loc_vel[2,*]*velVecFac
  loc_vel = !NULL

  ; save
  r = {loc_pos:loc_pos,loc_temp:loc_temp,loc_ent:loc_ent,loc_vrad:loc_vrad,loc_hsml:loc_hsml,$
       loc_pos2:loc_pos2,loc_mass:loc_mass,loc_metal:loc_metal,loc_ids:loc_ids,$
       loc_dynTime:loc_dynTime,loc_coolTime:loc_coolTime,loc_dens:loc_dens,$
       sP:sP,gcID:gcInd,boxCen:boxCen,boxSize:boxSize,sizeFac:sizeFac,$
       haloVirRad:haloVirRad,haloMass:haloMass,haloM200:haloM200,haloV200:haloV200}
         
  ;saveFilename = sP.derivPath + 'cutouts/cutout.' + sP.savPrefix + str(sP.res) + '.' + $
  ;str(sP.snap) + '.h' + str(gcInd) + '_subbox0.sav'
         
  ;save,r,filename=saveFilename
  ;print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  return,r
  
end

; makeSubboxFrames() : generate frames for a movie based on a subbox (either whole subbox or 
;   following one halo through time within that subbox) (single run or comparison of two runs)

pro makeSubboxFrames

  ; config
  zStart  = 9.0
  zEnd    = 2.0
  hInd    = -1 ; negative one to show whole subbox
  
  sP1 = simParams(res=512,run='tracer',redshift=zStart)
  sP2 = simParams(res=512,run='feedback',redshift=zStart)
  
  savePath = '/n/home07/dnelson/data3/frames/'
  
  ; vis config
  sizeFac          = 1.0 ; fraction of rvir if hInd>=0, or fraction of subbox if hInd=-1
  axes             = [0,1] ;list([0,1],[0,2],[1,2]) ;xy,xz,yz
  nPixels          = [800,800] ; px for sphMap
  singleColorScale = 0 ; 1=use same color scale for right panel, 0=rescale
  secondGt         = 0 ; 1=show greater than cut, 0=show less than cut
  nbottom          = 50
  ctName           = 'helix'
  
  ; use which field and minmax for color mapping? cut value for right panel?
  colorField = 'temp'     & fieldMinMax  = [4.0,7.0] & secondCutVal = 5.0
  ;colorField = 'entropy'  & fieldMinMax  = [6.5,8.5] & secondCutVal = 7.5 ; log(CGS)
  ;colorField = 'metal'     & fieldMinMax = [-4.0,-1.0] & secondCutVal = -2.5 ; log(Z/Zsun)
  ;colorField = 'vrad'      & fieldMinMax = [-400,400] & secondCutVal = -200.0 ; km/s
  ;colorField = 'vradnorm'  & fieldMinMax = [0.0,4.0] & secondCutVal = 2.0 ; vrad/v200
  ;colorField  = 'coolTime' & fieldMinMax = [0.0,5.0] & secondCutVal = 1.0 ; halo only
  ;colorField = 'dynTime'   & fieldMinMax = [0.0,0.8] & secondCutVal = 0.4 ; halo only
  ;colorField = 'timeRatio' & fieldMinMax = [0.0,4.0] & secondCutVal = 1.0 ; halo only
  
  ; load all subbox times, find ranges within [zStart,zEnd]
  snapZ1 = snapNumToRedshift(/all,/subBox,sP=sP1)
  w_snapZ1 = where(snapZ1 le zStart and snapZ1 ge zEnd,count1)
  
  snapZ2 = snapNumToRedshift(/all,/subBox,sP=sP2)
  w_snapZ2 = where(snapZ2 le zStart and snapZ2 ge zEnd,count2)
  
  print,'For run = '+sP1.run+' found: ['+str(count1)+'] snapshots.'
  print,'For run = '+sP2.run+' found: ['+str(count2)+'] snapshots.'
  
  nFrames = max([count1,count2])
  
  ; if subbox times don't match, take the number of frames as the larger, decide for each frame
  ; the corresponding subbox number for the lesser set (duplicate frames) (use value_locate)
  snaps = { run1 : intarr(nFrames), run2 : intarr(nFrames) }
  
  if count1 ge count2 then snaps.run1 = w_snapZ1
  if count1 le count2 then snaps.run2 = w_snapZ2
  
  if count1 lt count2 then snaps.run1 = value_locate(snapZ1,snapZ2[w_snapZ2]) > 0
  if count2 lt count1 then snaps.run2 = value_locate(snapZ2,snapZ1[w_snapZ1]) > 0
    
  ; if following one halo, calculate tracked position vs time from main snapshots, interpolate to subbox times
  if hInd ne -1 then begin
    gcInd = getMatchedIDs(sPa=sP1,sPg=sP2,haloID=haloID)
    boxSizeImg = [0,0]
    scaleBarLen = 200.0
    message,'todo'
  endif else begin
    boxSizeImg = sizeFac * sP1.subboxSize[axes]
    gcInd = {a:-1,g:-1}
    scaleBarLen = 800.0
  endelse
  
  ; loop over all frames
  for fN=0,nFrames-1 do begin
    print,'['+string(fN,format='(I4.4)')+'] progress: '+string(float(fN)/nFrames*100,format='(f5.1)')+$
      '%'+' (snap = '+string( (snaps.(0))[fN],format='(I4.4)' )+' / '+$
      string( (snaps.(1))[fN],format='(I4.4)' )+') (redshift = '+string(sP1.redshift,format='(f5.2)')+')'
    
    saveFilename = 'subbox_frame_' + string(fN,format='(I4.4)') + '.eps'
                 
    start_PS, savePath + saveFilename, xs=6, ys=9
  
    ; loop over the two runs
    for i=0,1 do begin
    
      if i eq 0 then sP = sP1
      if i eq 1 then sP = sP2
    
      sP.snap = (snaps.(i))[fN]
      sP.redshift = snapNumToRedshift(sP=sP,/subBox)
    
      ; load
      cutout = subBoxVisCutout(sP=sP,gcInd=gcInd.(i),sizeFac=sizeFac)
  
      ; calculate sph density weighted temperature projection map
      sphmap = calcSphMap(cutout.loc_pos,cutout.loc_hsml,cutout.loc_mass,cutout.loc_temp,$
                          boxSizeImg=cutout.boxSize,boxSizeSim=0,boxCen=[0,0,0],$
                          nPixels=nPixels,axes=axes,ndims=3)
  
      ; calculate colors
      config = {boxSizeImg:boxSizeImg,plotFilename:'',haloVirRad:cutout.haloVirRad,$
                haloMass:cutout.haloMass,axisPair:[0,0],sP:sP,singleColorScale:singleColorScale,$
                colorField:colorField,fieldMinMax:fieldMinMax,secondCutVal:secondCutVal,$
                secondGt:secondGt,nbottom:nbottom,barMM:fieldMinMax,ctName:ctName,barType:'2bar'}
    
      sub = cosmoVisCutoutSub(cutout=cutout,config=config)
    
      ; render scatterplot and sphmap
      scatter = {pos_left:sub.pos_left,  pos2_left:sub.pos_left2,  cinds_left:sub.cinds_left,$
                 pos_right:sub.pos_right,pos2_right:sub.pos_right2,cinds_right:sub.cinds_right}  
    
      config = {saveFilename:'',nPixels:nPixels,axes:axes,fieldMinMax:fieldMinMax,$
                gcID:gcInd.(i),haloMass:cutout.haloMass,haloVirRad:cutout.haloVirRad,$
                boxCen:[0,0,0],boxSizeImg:cutout.boxSize,boxSizeScat:cutout.boxSize,$
                ctNameScat:'helix',ctNameMap:'blue-red2',sP:sP,bartype:'none',scaleBarLen:scaleBarLen}

      plotScatterAndSphmap, map=sphmap, scatter=scatter, config=config, left=(i eq 0), right=(i eq 1)
    
    endfor ;i
    
    end_PS, pngResize=60, /deletePS
  
  endfor
  
  print,'Done.'
  
end

; makeProgressionPanels(): for one run, or comparing two runs, do panels or individual frames of 
;  plotScatterComp for a halo tracked through time

pro makeProgressionPanels

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  ;sP = simParams(res=256,run='feedback',redshift=0.0)
  sP = simParams(res=512,run='zoom_20Mpc',redshift=2.0,hInd=0)
  
  ; best match in sP2 is found for each sP snap, comment out sP2 to not make a comparison
  ;sP2 = simParams(res=256,run='tracer',redshift=0.0)

  haloID        = 0   ; zoom.0 z2.314 z2.304 z2.301 z2.130 z2.64
  redshiftStart = 2.0 ; when to start (can differ from sP.redshift, which is where the halo is targeted)
  redshiftEnd   = 4.0 ; how far back to go
  indivFrames   = 1   ; 1=write individual frames (for a movie) instead of all panels in one big image
  
  singleColorScale = 0 ; 1=use same color scale for right panel, 0=rescale
  secondGt         = 0 ; 1=show greater than cut, 0=show less than cut
  axes             = [0,1]
  nbottom          = 50
  sizeFac          = 3.5 ; times rvir
  ctName           = 'helix'
  pngResize        = 80 ; 40=900px wide for one sim, 1800px for two sims
  deletePS         = 1  ; 0=false, 1=true
  
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
  if indivFrames then numSnaps = 1 ; override
  
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
  
  if ~indivFrames then start_PS, sP.plotPath + pFilename, xs=pWidth, ys=pHeight

  ; loop over snapshots
  for m=sP.snap,minSnap,-1 do begin
    sP.snap = m & print,m,m-minSnap,gcIDList[sP.snap]
    sP.redshift = snapNumToRedshift(sP=sP)
    
    ; individual frames? set filename and open plot
    if indivFrames then begin
      pFilename = 'scatter_frame_'+string(m-minSnap,format='(I4.4)')+'.eps'
      start_PS, sP.plotPath + pFilename, xs=pWidth, ys=pHeight
    endif
    
    ; load cutouts
    cutout = cosmoVisCutout(sP=sP,gcInd=gcIDList[sP.snap],sizeFac=sizeFac)
    
    config = {boxSizeImg:cutout.boxSizeImg,plotFilename:'',haloVirRad:cutout.haloVirRad,$
              haloMass:cutout.haloMass,axisPair:axes,sP:sP,singleColorScale:singleColorScale,$
              colorField:colorField,fieldMinMax:fieldMinMax,secondCutVal:secondCutVal,$
              secondGt:secondGt,nbottom:nbottom,barMM:fieldMinMax,ctName:ctName,barType:'2bar',$
              lineThick:3.0}
    
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
    row = [m-minSnap,numSnaps]
    if indivFrames then row = []
    
    plotScatterComp,sub=sub,config=config,row=row,first=(n_elements(sP2) gt 0)
    
    if n_elements(sP2) gt 0 then $
      plotScatterComp,sub=sub2,config=mod_struct(config,'subtitle',/delete),row=row,/second
        
    if indivFrames then end_PS, pngResize=pngResize, deletePS=deletePS
  
  endfor
  
  ; finish plot
  if ~indivFrames then end_PS, pngResize=pngResize, deletePS=deletePS
  stop
end
