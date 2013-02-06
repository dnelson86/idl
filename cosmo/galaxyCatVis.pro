; galaxyCatVis.pro
; 2d visualization (using cosmoVis.pro) specialized for the galaxyCat
; dnelson dec.2012

; galaxyCatVisCutout(): make a cutout of a halo (gmem only) based on the galaxy catalog
;                       call with multiple gcInd's for one load and save cutouts
;                       call with one gcInd to return results

function galaxyCatVisCutout, sP=sP, gcInd=gcInd

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  velVecFac = 0.01 ; times velocity (km/s) in plotted kpc
  hsmlFac = 1.75   ; increase arepo 'hsml' to decrease visualization noise  
  
  h = loadSnapshotHeader(sP=sP)

  ; check existence of requested saves if more than one halo
  saveFilenames = sP.derivPath + 'cutouts/gmemCut.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + $
                  '.h' + str(gcInd) + '.sav'
                  
  readFlag = 0
  foreach saveFilename,saveFilenames do if ~file_test(saveFilename) then readFlag = 1
  
  ; if single halo requested and save exists, load it
  if n_elements(saveFilenames) eq 1 then $
    if file_test(saveFilenames) then restore,saveFilename

  ; proceed with cutouts if at least one save is missing
  if readFlag then begin
  
  gc     = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen  = subgroupPosByMostBoundID(sP=sP)
  
  ; load u,nelec,dens and calculate temperature,entropy
  u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
  nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
  temp  = alog10(convertUtoTemp(u,nelec))
  nelec = !NULL
  dens  = loadSnapshotSubset(sP=sP,partType='gas',field='dens')
  ent   = calcEntropyCGS(u,dens,/log,sP=sP)
  u     = !NULL
  dens  = !NULL
  
  ; load HSMLs or volumes (convert to sizes)
  if sP.trMCPerCell eq 0 then begin
    hsml = loadSnapshotSubset(sP=sP,partType='gas',field='hsml')
    hsml = 1.0 * temporary(hsml); increase hsml to decrease visualization noise
  endif else begin
    hsml = loadSnapshotSubset(sP=sP,partType='gas',field='vol')
    hsml = (temporary(hsml) * 3.0 / (4*!pi))^(1.0/3.0) ;cellrad [ckpc]
    hsml = hsmlFac * temporary(hsml) ; increase hsml to decrease visualization noise
  endelse  
  
  ; load gas ids, positions, velocities and masses
  ids  = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
  pos  = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  vel  = loadSnapshotSubset(sP=sP,partType='gas',field='vel')
  mass = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
  dens = loadSnapshotSubset(sP=sP,partType='gas',field='dens')
  
  ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
  print,'shuffling...'
  iseed = 424242L
  sort_inds = sort(randomu(iseed,n_elements(temp)))
  
  ids  = ids[sort_inds]
  temp = temp[sort_inds]
  ent  = ent[sort_inds]
  hsml = hsml[sort_inds]
  mass = mass[sort_inds]
  dens = dens[sort_inds]
  pos  = pos[*,sort_inds]
  vel  = vel[*,sort_inds]
  
  sort_inds = !NULL
  
  ; now restrict all these quantities to gmem only
  galcat = galaxyCat(sP=sP)
  match,ids,galcat.groupmemIDs,ids_ind,gmem_ind,count=countGmem
  ids = !NULL
  if countGmem ne n_elements(galcat.groupmemIDs) then message,'Error: Failed to find all gmem in gas ids.'
  ids_ind = ids_ind[sort(gmem_ind)]
  
  temp = temp[ids_ind]
  ent  = ent[ids_ind]
  hsml = hsml[ids_ind]
  mass = mass[ids_ind]
  dens = dens[ids_ind]
  pos  = pos[*,ids_ind]
  vel  = vel[*,ids_ind]
  
  ; load cooling and dynamical timescales (already gmem only)
  encMass = enclosedMass(sP=sP) ; code units

  gasRadii = galaxyCatRadii(sP=sP)
  gasRadii = gasRadii.gmem_sec
  meanDensEnc = 3*encMass / (4 * !pi * gasRadii^3.0) / (h.time)^3.0 ; code units (physical)
  dynTime = sqrt( 3*!pi / (32 * float(units.G) * meanDensEnc * units.HubbleParam) ) ; code units (Gyr)
  encMass = !NULL
  gasRadii = !NULL
  
  ct = coolingTime(sP=sP)
  coolTime = ct.coolTime
  ct = !NULL
  
  print,'cutout...'
  foreach gcIndCur,gcInd,k do begin
    ; get indices into gmem array for this halo
    loc_inds = galcatINDList(sP=sP, galcat=galcat, gcIDList=[gcIndCur])
    loc_inds = loc_inds.gmem
    
    nCutout  = n_elements(loc_inds)
    boxCen   = sgcen[*,gcIndCur]
    
    ; make positions relative
    xDist = pos[0,loc_inds] - boxCen[0]
    yDist = pos[1,loc_inds] - boxCen[1]
    zDist = pos[2,loc_inds] - boxCen[2]
    
    correctPeriodicDistVecs, xDist, sP=sP
    correctPeriodicDistVecs, yDist, sP=sP
    correctPeriodicDistVecs, zDist, sP=sP
    
    loc_pos  = fltarr(3,nCutout)
    loc_pos[0,*] = xDist
    loc_pos[1,*] = yDist
    loc_pos[2,*] = zDist
    
    xDist = !NULL
    yDist = !NULL
    zDist = !NULL
    
    ; take subsets of other gmem arrays
    loc_temp = temp[loc_inds]
    loc_ent  = ent[loc_inds]
    loc_hsml = hsml[loc_inds]
    loc_mass = mass[loc_inds]
    loc_dens = dens[loc_inds]
    loc_vel  = vel[*,loc_inds]
    
    loc_dynTime  = dynTime[loc_inds]
    loc_coolTime = coolTime[loc_inds]
    
    ; determine a bounding box
    boxSize    = ceil(max(loc_pos) / 10.0) * 10.0
    boxSizeImg = [boxSize,boxSize,boxSize] ; cube
    
    ; calculate norm of radial velocity vector
    rad = reform(loc_pos[0,*]^2.0 + loc_pos[1,*]^2.0 + loc_pos[2,*]^2.0)
    const_vrad = (loc_vel[0,*] * loc_pos[0,*] + $
                  loc_vel[1,*] * loc_pos[1,*] + $
                  loc_vel[2,*] * loc_pos[2,*]) / rad
                
    vrad = fltarr(3,nCutout)
    for i=0,2 do vrad[i,*] = const_vrad * loc_pos[i,*]
    const_vrad = !NULL
    rad = !NULL
    
    loc_vrad = sqrt(vrad[0,*]^2.0 + vrad[1,*]^2.0 + vrad[2,*]^2.0)
    
    ; create endpoint for each position point for the velocity vector line
    loc_pos2 = fltarr(3,nCutout)
    loc_pos2[0,*] = loc_pos[0,*] + loc_vel[0,*]*velVecFac
    loc_pos2[1,*] = loc_pos[1,*] + loc_vel[1,*]*velVecFac
    loc_pos2[2,*] = loc_pos[2,*] + loc_vel[2,*]*velVecFac
    loc_vel = !NULL
    
    haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcIndCur]] ;ckpc
    haloMass = codeMassToLogMsun(gc.subgroupMass[gcIndCur])
    haloM200 = codeMassToLogMsun(gc.group_m_crit200[gc.subgroupGrNr[gcIndCur]])
    
    ; fix halo mass if we're using old (x2 bug) catalogs
    if sP.run eq 'gadgetold' or sP.run eq 'arepo' then $
      haloMass = codeMassToLogMsun(0.5*gc.subgroupMass[gcIndCur])  
  
    ; save
    r = {loc_pos:loc_pos,loc_temp:loc_temp,loc_ent:loc_ent,loc_vrad:loc_vrad,loc_hsml:loc_hsml,$
	   loc_pos2:loc_pos2,loc_mass:loc_mass,sP:sP,gcID:gcIndCur,boxCen:boxCen,$
         loc_dynTime:loc_dynTime,loc_coolTime:loc_coolTime,loc_dens:loc_dens,$
         boxSizeImg:boxSizeImg,haloVirRad:haloVirRad,haloMass:haloMass}
         
    saveFilename = sP.derivPath + 'cutouts/gmemCut.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + $
                   '.h' + str(gcIndCur) + '.sav'  
         
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  endforeach
  
  endif ; readFlag
  
  if n_elements(gcInd) eq 1 then return,r
  return,1
  
end

; galCatScatterMapHalos: plot temperature colored scatter plots with velocity vectors on boxes centered on halos

pro galCatScatterMapHalos;, sP=sP, gcIDs=gcIDs

  sP = simParams(res=512,run='tracer',redshift=2.0)
  units = getUnits()

  haloID = 304 ;z2.304 z2.301 z2.130 z2.64
  gcID = getMatchedIDs(sPa=sP,sPg=sP,haloID=haloID)
  gcIDs = [gcID.a]

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  singleColorScale = 0 ; 1=use same color scale for right panel, 0=rescale
  secondGt         = 0 ; 1=show greater than cut, 0=show less than cut
  
  ; use which field and cut value for right panel?
  ;secondField = 'temp'     & secondCutVal = 5.0
  ;secondField = 'entropy'  & secondCutVal = 7.5
  ;secondField = 'vrad'     & secondCutVal = 200.0
  ;secondField = 'vradnorm' & secondCutVal = 2.0
  ;secondField = 'coolTime'  & secondCutVal = 1.0
  ;secondField = 'dynTime'   & secondCutVal = 0.0
  secondField = 'timeRatio' & secondCutVal = 1.0
  
  ; use which field and minmax for color mapping?
  ;colorField = 'temp'     & fieldMinMax  = [4.0,7.0]
  ;colorField = 'entropy'  & fieldMinMax  = [6.5,8.5] ; log(CGS)
  ;colorField = 'vrad'     & fieldMinMax = [0.0,400.0] ; km/s
  ;colorField = 'vradnorm' & fieldMinMax  = [0.0,4.0] ; vrad/v200
  ;colorField  = 'coolTime' & fieldMinMax = [0.0,5.0]
  ;colorField = 'dynTime'  & fieldMinMax = [0.0,5.0]
  colorField = 'timeRatio' & fieldMinMax = [0.0,4.0]
  
  axes = list([0,1],[0,2],[1,2]) ;xy,xz,yz
  
  ; make cutouts (multiple or single)
  cutout = galaxyCatVisCutout(sP=sP,gcInd=gcIDs)

  ; target list
  gc    = loadGroupCat(sP=sP,/skipIDs)
  sgcen = subgroupPosByMostBoundID(sP=sP) 
  
  print,'rendering...'
  ; loop over all requested halos and image
  foreach gcID, gcIDs do begin
  
    ; load cutout
    cutout = galaxyCatVisCutout(sP=sP,gcInd=gcID)

    v200 = sqrt(units.G * (10.0^cutout.haloMass / units.UnitMass_in_Msun) / cutout.haloVirRad )

    ; create color index mapping
    if colorField eq 'temp'      then fieldVal = cutout.loc_temp
    if colorField eq 'entropy'   then fieldVal = cutout.loc_ent
    if colorField eq 'vrad'      then fieldVal = reform(cutout.loc_vrad)
    if colorField eq 'vradnorm'  then fieldVal = reform(cutout.loc_vrad)/v200
    if colorField eq 'coolTime'  then fieldVal = cutout.loc_coolTime
    if colorField eq 'dynTime'   then fieldVal = cutout.loc_dynTime
    if colorField eq 'timeRatio' then fieldVal = cutout.loc_coolTime / cutout.loc_dynTime
    
    colorinds = (fieldVal-fieldMinMax[0])*205.0 / (fieldMinMax[1]-fieldMinMax[0]) ;0-205
    colorinds = fix(colorinds + 50.0) > 0 < 255 ;50-255  
  
    ; second/right panel cutout
    wSecond = where(fieldVal le secondCutVal,nCutoutSecond,comp=wComp)
      
    ; show gas above this cut value (instead of below)?
    if secondGt then wSecond = wComp
    
    loc_pos_second   = cutout.loc_pos[*,wSecond]
    loc_pos2_second  = cutout.loc_pos2[*,wSecond]
    colorinds_second = colorinds[wSecond]
    
    ; use instead a differently scaled color mapping for the second panel?
    if singleColorScale eq 0 then begin
      if secondGt eq 1 then $
        colorinds_second = (fieldVal-secondCutVal)*205.0 / (fieldMinMax[1]-secondCutVal)
      if secondGt eq 0 then $
        colorinds_second = (fieldVal-fieldMinMax[0])*205.0 / (secondCutVal-fieldMinMax[0])
      
      colorinds_second = colorinds_second[wSecond]
    endif

    ; make a plot for each requested projection direction
    foreach axisPair, axes do begin
           
      ; get box center (in terms of specified axes)
      boxCenImg  = [sgcen[axisPair[0],gcID],sgcen[axisPair[1],gcID],sgcen[3-axisPair[0]-axisPair[1],gcID]]
      
      print,'['+string(gcID,format='(i4)')+'] Mapping ['+str(axisPair[0])+' '+$
            str(axisPair[1])+'] with '+str(cutout.boxSizeImg[0])+$
            ' kpc box around subhalo center ['+str(boxCenImg[0])+' '+str(boxCenImg[1])+' '+str(boxCenImg[2])+']'

      plotFilename = 'scatter.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+$
                     '.axes'+str(axisPair[0])+str(axisPair[1])+'-'+$
                     colorField+'-'+secondField+'-'+str(secondGt)+'sCS'+str(singleColorScale)+'.eps'

      config = {boxSizeImg:cutout.boxSizeImg*1.5,plotFilename:plotFilename,haloVirRad:cutout.haloVirRad,$
                haloMass:cutout.haloMass,axisPair:axisPair,sP:sP,$
                colorField:colorField,fieldMinMax:fieldMinMax,secondCutVal:secondCutVal,secondGt:secondGt,$
                barMM:fieldMinMax,barType:'2bar'}
      
      ; plot
      plotScatterComp,cutout.loc_pos,cutout.loc_pos2,loc_pos_second,loc_pos2_second,$
        colorinds,colorinds_second,config=config            

    endforeach ;axisPair

  endforeach ;gcIDs
  
  stop

end
