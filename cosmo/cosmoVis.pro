; cosmoVis.pro
; cosmological boxes - 2d visualization
; dnelson apr.2013

; cosmoVisCutout(): make a spatial cutout around a halo
;                   call with multiple gcInd's for one load and save cutouts
;                   call with one gcInd to return results
;                   call with selectGmem

function cosmoVisCutout, sP=sP, gcInd=gcInd, sizeFac=sizeFac, selectGmem=selectGmem

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  velVecFac = 0.01 ; times velocity (km/s) in plotted kpc
  hsmlFac   = 1.75 ; increase arepo 'hsml' to decrease visualization noise  

  if keyword_set(sizeFac) and keyword_set(selectGmem) then message,'The two are exclusive.'
  if ~keyword_set(sP) or n_elements(gcInd) eq 0 then message,'Error'
  
  ; filename tag for selectGmem or not
  cutTag = 'cutout'
  if keyword_set(selectGmem) then begin
    cutTag = 'gmemCut'
    sizeFac = 0
  endif
  
  ; check existence of requested saves if more than one halo
  saveFilenames = sP.derivPath + 'cutouts/' + cutTag + '.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.h' + str(gcInd) + '.sf' + str(fix(sizeFac*10)) + '.sav'
        
  ; if single halo requested and save exists, load it
  if n_elements(saveFilenames) eq 1 then $
    if file_test(saveFilenames) then restore,saveFilenames
        
  ; more than one halo requested, see if any need to be made
  readFlag = 0
  foreach saveFilename,saveFilenames do if ~file_test(saveFilename) then readFlag = 1
    
  ; proceed with cutouts if at least one save is missing
  if readFlag then begin
    
  gc    = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sP)   
  
  ; have metallicities to include in cutout?
  metalsFlag = 0
  if snapshotFieldExists(sP=sP,field='Metallicity') or $
     snapshotFieldExists(sP=sP,field='GFM_Metallicity') then metalsFlag = 1
  
  ; load u,nelec,dens and calculate temperature,entropy
  u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
  nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
  temp  = alog10(convertUtoTemp(u,nelec))
  nelec = !NULL
  dens  = loadSnapshotSubset(sP=sP,partType='gas',field='dens')
  ent   = calcEntropyCGS(u,dens,/log,sP=sP)
  u     = !NULL
  
  if metalsFlag then begin
    metal = loadSnapshotSubset(sP=sP,partType='gas',field='metallicity')
    ; convert to log(metallicity) for positive values, otherwise set GFM_MIN_METAL = -20 (log)
    w = where(metal gt 0.0,count,comp=wc,ncomp=ncomp)
    if count gt 0 then metal[w] = alog10(metal[w])
    if ncomp gt 0 then metal[wc] = -20.0
  endif
    
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
  
  ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
  print,'shuffling...'
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
  
  ; now restrict all these quantities to gmem only if requested
  if keyword_set(selectGmem) then begin
    h = loadSnapshotHeader(sP=sP)
    galcat = galaxyCat(sP=sP)
    
    idIndexMap = getIDIndexMap(ids,minid=minid)
    ids_ind = idIndexMap[galcat.groupmemIDs-minid]
  
    ids  = ids[ids_ind]
    temp = temp[ids_ind]
    ent  = ent[ids_ind]
    dens = dens[ids_ind]
    hsml = hsml[ids_ind]
    mass = mass[ids_ind]
    dens = dens[ids_ind]
    pos  = pos[*,ids_ind]
    vel  = vel[*,ids_ind]
    
    if metalsFlag then metal = metal[ids_ind]
  
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
  endif
  
  print,'cutout...'
  foreach gcIndCur,gcInd,k do begin
    ; halo properties
    haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcIndCur]] ;ckpc
    haloMass = codeMassToLogMsun(gc.subgroupMass[gcIndCur])
    haloM200 = codeMassToLogMsun(gc.group_m_crit200[gc.subgroupGrNr[gcIndCur]])
    haloV200 = sqrt(units.G * gc.subgroupMass[gcIndCur] / haloVirRad )
    
    ; decide cutout
    if keyword_set(sizeFac) then begin
      ; get subhalo position and size of imaging box
      boxCen     = sgcen[*,gcIndCur]
      boxSize    = ceil(sizeFac * haloVirRad / 10.0) * 10.0
    
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
    endif else begin
      ; gmem selection (all quantities are in gmem catalog order)
      wCut = galcatINDList(sP=sP, galcat=galcat, gcIDList=[gcIndCur])
      wCut = wCut.gmem
    
      nCutout  = n_elements(wCut)
    endelse
                 
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
    
    if keyword_set(selectGmem) then begin
      loc_dynTime  = dynTime[wCut]
      loc_coolTime = coolTime[wCut]
    endif else begin
      loc_dynTime = -1
      loc_coolTime = -1
    endelse
    
    ; calculate norm of radial velocity vector
    rad = reform(loc_pos[0,*]^2.0 + loc_pos[1,*]^2.0 + loc_pos[2,*]^2.0)

    sgcen_vel = gc.subgroupVel[*,gcIndCur]
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

    ; fix halo mass if we're using old (x2 bug) catalogs
    if sP.run eq 'gadgetold' or sP.run eq 'arepo' then $
      haloMass = codeMassToLogMsun(0.5*gc.subgroupMass[gcIndCur])  
  
    boxSizeImg = [boxSize,boxSize,boxSize] ; cube
    
    ; save
    r = {loc_pos:loc_pos,loc_temp:loc_temp,loc_ent:loc_ent,loc_vrad:loc_vrad,loc_hsml:loc_hsml,$
	   loc_pos2:loc_pos2,loc_mass:loc_mass,loc_metal:loc_metal,loc_ids:loc_ids,$
         loc_dynTime:loc_dynTime,loc_coolTime:loc_coolTime,loc_dens:loc_dens,$
         sP:sP,gcID:gcIndCur,boxCen:boxCen,boxSizeImg:boxSizeImg,sizeFac:sizeFac,$
         haloVirRad:haloVirRad,haloMass:haloMass,haloM200:haloM200,haloV200:haloV200}
         
    saveFilename = sP.derivPath + 'cutouts/' + cutTag + '.' + sP.savPrefix + str(sP.res) + '.' + $
      str(sP.snap) + '.h' + str(gcIndCur) + '.sf' + str(fix(sizeFac*10)) + '.sav'  
         
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  endforeach
  
  endif ; readFlag
  
  if n_elements(gcInd) eq 1 then return,r
  return,1
  
end

; plotScatterComp(): plot side by side colored/vectorized scatter plots

pro plotScatterComp, pos_left, pos2_left, pos_right, pos2_right, cinds_left, cinds_right, $
                     config=config, top=top, bottom=bottom, subtitle=subtitle

      xMinMax = [-config.boxSizeImg[0]/2.0,config.boxSizeImg[0]/2.0]
      yMinMax = [-config.boxSizeImg[1]/2.0,config.boxSizeImg[1]/2.0]
      
      plotPath = '/n/home07/dnelson/coldflows/'
      
      posLeft = [0.0,0.0,0.5,1.0]
      if keyword_set(top) then posLeft = [0.0,0.5,0.5,1.0]
      if keyword_set(bottom) then posLeft = [0.0,0.0,0.5,0.5]
      posRight = posLeft + [0.5,0.0,0.5,0.0]
      
      if ~keyword_set(top) and ~keyword_set(bottom) then $
        start_PS, plotPath + config.plotFilename, xs=8, ys=4
      
        !p.thick = 3.0
        !p.charsize = 0.8
      
        ; fill with black background
        if ~keyword_set(bottom) then $
          cgColorfill,[1,1,0,0,1],[1,0,0,1,1],/normal,color=cgColor('black')
      
        ; color table and establish temperature -> color mapping
        loadColorTable,config.ctName
        
        ; (left panel)
        ; ------------
        cgPlot, /nodata, xMinMax, yMinMax, pos=posLeft, xs=5, ys=5, /noerase
        
        ; circle at virial radius
        tvcircle,config.haloVirRad,0,0,cgColor('light gray'),thick=0.6,/data

        ; particle loop for velocity vector plotting
        nCutoutLeft = n_elements(pos_left[0,*])
        for i=0L,nCutoutLeft-1 do $
          oplot,[pos_left[config.axisPair[0],i],pos2_left[config.axisPair[0],i]],$
                 [pos_left[config.axisPair[1],i],pos2_left[config.axisPair[1],i]],$
                 line=0,color=cinds_left[i]
        
        ; scale bar
        if ~keyword_set(bottom) then begin
          len = 100.0 ;ckpc
          cgText,mean([-config.boxSizeImg[0]/2.2,-config.boxSizeImg[0]/2.2+len]),config.boxSizeImg[0]/2.3,$
                 string(len,format='(i3)')+' kpc',alignment=0.5,color=cgColor('light gray')
          cgPlot,[-config.boxSizeImg[0]/2.2,-config.boxSizeImg[0]/2.2+len],$
                 [config.boxSizeImg[1]/2.1,config.boxSizeImg[1]/2.1],$
                 color=cgColor('light gray'),thick=4.0,/overplot
        endif
        
        if keyword_set(bottom) then $
          cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=1.0,color=cgColor('dark gray'),/overplot
               
        ; (right panel)
        ; -------------
        loadColorTable,config.ctName
        
        cgPlot, /nodata, xMinMax, yMinMax, pos=posRight, xs=5, ys=5, /noerase

        tvcircle,config.haloVirRad,0,0,cgColor('light gray'),thick=0.6,/data
        
        ; particle loop for velocity vector plotting (cold gas only)
        nCutoutRight = n_elements(pos_right[0,*])
        for i=0L,nCutoutRight-1 do $
          oplot,[pos_right[config.axisPair[0],i],pos2_right[config.axisPair[0],i]],$
                 [pos_right[config.axisPair[1],i],pos2_right[config.axisPair[1],i]],$
                 line=0,color=cinds_right[i]
        
        ; redshift and halo mass
        if ~keyword_set(bottom) then begin
          cgText,config.boxSizeImg[0]/2.1,config.boxSizeImg[1]/2.4,alignment=1.0,$
            "z = "+string(config.sP.redshift,format='(f3.1)'),color=cgColor('light gray')
          cgText,config.boxSizeImg[0]/2.1,config.boxSizeImg[1]/2.2,alignment=1.0,$
            "M = "+string(config.haloMass,format='(f4.1)'),color=cgColor('light gray')
        endif
        
        ; dividing line
        cgPlot,[xMinMax[0],xMinMax[0]],yMinMax,line=0,thick=1.0,color=cgColor('dark gray'),/overplot
        
        if keyword_set(bottom) then $ ; horizontal dividing line
          cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=1.0,color=cgColor('dark gray'),/overplot
        
        ; colorbar(s) on bottom
        loadColorTable,config.ctName
        
        !x.thick = 1.0
        !y.thick = 1.0
        
        if ~keyword_set(top) then begin
        
        ; choose center colorbar label
        if config.colorField eq 'vrad'      then labelText = "v_{rad} [km/s]"
        if config.colorField eq 'vradnorm'  then labelText = "v_{rad} / v_{200}"
        if config.colorField eq 'temp'      then labelText = "log T_{gas} [K]"
        if config.colorField eq 'entropy'   then labelText = "log (S) [cgs]"
        if config.colorField eq 'density'   then labelText = "\rho_{gas}"
        if config.colorField eq 'metal'     then labelText = "log Z"
        if config.colorField eq 'overdens'  then labelText = "log \rho_{DM} / <\rho_{DM}>"
        if config.colorField eq 'coolTime'  then labelText = "t_{cool} [Gyr]"
        if config.colorField eq 'dynTime'   then labelText = "t_{dyn} [Gyr]"
        if config.colorField eq 'timeRatio' then labelText = "t_{cool} / t_{dyn}"

        ; first text: integer or has a decimal?
        if round(config.fieldMinMax[1]) eq config.fieldMinMax[1] then begin
          firstText = str(fix(config.fieldMinMax))
          extraPad = 0.0
            
          ; more than one digit?
          if alog10(abs(config.fieldMinMax[1])) ge 1.0 then $
            extraPad += floor(alog10(abs(config.fieldMinMax[1]))) * 0.007
        endif else begin
          firstText = string(config.fieldMinMax,format='(f3.1)')
          extraPad = 0.007
        endelse
        
        if config.barType eq '1bar' then begin
          ; one colorbar (centered)
          pos = [0.35,0.02,0.65,0.078]
          if keyword_set(bottom) then pos *= [1.0,0.5,1.0,0.5]
          cgColorbar,position=pos,divisions=0,charsize=0.000001,ticklen=0.00001,$
            bottom=config.nbottom,ncolors=(255-3-config.nbottom) ; 3 to remove white band at rightside
          
          ; colorbar labels
          cgText,0.5,0.021,textoidl(labelText),alignment=0.5,color=cgColor('black'),/normal
          cgText,pos[0]+0.015+extraPad,0.0195,firstText[0],alignment=0.5,color=cgColor('black'),/normal
          cgText,pos[2]-0.015-extraPad,0.0195,firstText[1],alignment=0.5,color=cgColor('black'),/normal
        endif
        
        if config.barType eq '2bar' then begin
          ; two colorsbars (separate ranges)
          cgColorbar,position=[0.02,0.1,0.076,0.4],divisions=0,charsize=0.000001,$
            bottom=config.nbottom,ncolor=(255-config.nbottom),ticklen=0.00001
          
          ; second (right panel)
          loc_mm = [config.fieldMinMax[0],config.secondCutVal]
          if config.secondGt then loc_mm = [config.secondCutVal,config.fieldMinMax[1]]
          
          loadColorTable,config.ctName
          
          cgColorbar,position=[0.02,0.6,0.076,0.9],divisions=0,charsize=0.000001,$
            bottom=config.nbottom,ncolor=(255-config.nbottom),ticklen=0.00001
          
          ; first text
          cgText,0.25,0.0375,textoidl(labelText),alignment=0.5,color=cgColor('black'),/normal
          cgText,0.115+extraPad,0.036,firstText[0],alignment=0.5,color=cgColor('black'),/normal
          cgText,0.385-extraPad,0.036,firstText[1],alignment=0.5,color=cgColor('black'),/normal
          
          ; second text: integer or has a decimal?
          if round(loc_mm[1]) eq loc_mm[1] and round(loc_mm[0]) eq loc_mm[0] then begin
            secondText = str(fix(loc_mm))
            extraPad = 0.0
            
            ; more than one digit?
            if alog10(abs(loc_mm[1])) ge 1.0 then $
              extraPad += floor(alog10(abs(loc_mm[1]))) * 0.007
          endif else begin
            secondText = str(string(loc_mm,format='(f4.1)'))
            extraPad = 0.007
          endelse
            
          cgText,0.75,0.0375,textoidl(labelText),alignment=0.5,color=cgColor('black'),/normal
          cgText,0.615+extraPad,0.036,secondText[0],alignment=0.5,color=cgColor('black'),/normal
          cgText,0.885-extraPad,0.036,secondText[1],alignment=0.5,color=cgColor('black'),/normal
        endif
        
        endif ;top
        
        ; simulation name
        if keyword_set(top) and ~keyword_set(subtitle) then cgText,0.5,0.96,"GADGET",charsize=!p.charsize+0.5,$
          alignment=0.5,/normal,color=cgColor('white')
        if keyword_set(bottom) and ~keyword_set(subtitle) then cgText,0.5,0.46,"AREPO",charsize=!p.charsize+0.5,$
          alignment=0.5,/normal,color=cgColor('white')
          
        ; or subtitle for each panel
        if keyword_set(subtitle) then begin
          if keyword_set(top) then ypos = 0.52
          if keyword_set(bottom) then ypos = 0.462
          cgText,0.48,ypos,subtitle[0],alignment=1.0,/normal,charsize=!p.charsize+0.5,color=cgColor('white')
          cgText,0.52,ypos,subtitle[1],alignment=0.0,/normal,charsize=!p.charsize+0.5,color=cgColor('white')
        endif
        
      if ~keyword_set(top) and ~keyword_set(bottom) then $        
        end_PS, pngResize=60;, /deletePS
end

; scatterMapHalos: plot colored scatter plots with velocity vectors on boxes centered on halos
; showGmem/showAll (use galaxy catalog and display gmem only, or all in spatial subset)

pro scatterMapHalos, showGmem=showGmem, showAll=showAll

  compile_opt idl2, hidden, strictarr, strictarrsubs
  if ~keyword_set(showGmem) and ~keyword_set(showAll) then message,'Choose one.'
  
  ; config
  sP = simParams(res=512,run='tracer',redshift=2.0)
  units = getUnits()

  haloID = 304 ;z2.304 z2.301 z2.130 z2.64
  gcID = getMatchedIDs(sPa=sP,sPg=sP,haloID=haloID)
  gcIDs = [gcID.a]

  if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; compare to a second run (2x2 panels instead of 2x1)?
  ; TODO
  ; xs=8, ys=8, /top, /bottom
  
  ; plot config
  singleColorScale = 0 ; 1=use same color scale for right panel, 0=rescale
  secondGt         = 1 ; 1=show greater than cut, 0=show less than cut
  axes             = list([0,1]) ;list([0,1],[0,2],[1,2]) ;xy,xz,yz
  nbottom          = 50
  sizeFac          = 3.5 ; times rvir
  ctName           = 'helix'
  
  ; use which field and cut value for right panel?
  ;secondField = 'temp'     & secondCutVal = 5.0
  ;secondField = 'entropy'  & secondCutVal = 7.5
  secondField = 'metal'     & secondCutVal = -2.5
  ;secondField = 'vrad'      & secondCutVal = -200.0
  ;secondField = 'vradnorm'  & secondCutVal = 2.0
  ;secondField = 'coolTime'  & secondCutVal = 1.0 ; gmem only
  ;secondField = 'dynTime'   & secondCutVal = 0.4 ; gmem only
  ;secondField = 'timeRatio' & secondCutVal = 1.0 ; gmem only
  
  ; use which field and minmax for color mapping?
  ;colorField = 'temp'     & fieldMinMax  = [4.0,7.0]
  ;colorField = 'entropy'  & fieldMinMax  = [6.5,8.5] ; log(CGS)
  colorField = 'metal'     & fieldMinMax = [-4.0,-1.0] ; log(Z/Zsun)
  ;colorField = 'vrad'      & fieldMinMax = [-400,400] ; km/s
  ;colorField = 'vradnorm'  & fieldMinMax = [0.0,4.0] ; vrad/v200
  ;colorField  = 'coolTime' & fieldMinMax = [0.0,5.0] ; gmem only
  ;colorField = 'dynTime'   & fieldMinMax = [0.0,0.8] ; gmem only
  ;colorField = 'timeRatio' & fieldMinMax = [0.0,4.0] ; gmem only

  ; pre-make cutouts (multiple or single)
  ;if keyword_set(showGmem) then cutout = cosmoVisCutout(sP=sP,gcInd=gcIDs,/selectGmem)
  ;if keyword_set(showAll)  then cutout = cosmoVisCutout(sP=sP,gcInd=gcIDs,sizeFac=sizeFac)

  ; target list
  gc    = loadGroupCat(sP=sP,/skipIDs)
  sgcen = subgroupPosByMostBoundID(sP=sP) 
  
  print,'rendering...'
  ; loop over all requested halos and image
  foreach gcID, gcIDs do begin
    
    ; load cutout
    if keyword_set(showGmem) then cutout = cosmoVisCutout(sP=sP,gcInd=gcID,/selectGmem)
    if keyword_set(showAll)  then cutout = cosmoVisCutout(sP=sP,gcInd=gcID,sizeFac=sizeFac)
    
    ; create color index mapping
    if colorField eq 'temp'      then fieldVal = cutout.loc_temp
    if colorField eq 'entropy'   then fieldVal = cutout.loc_ent
    if colorField eq 'density'   then fieldVal = cutout.loc_dens
    if colorField eq 'metal'     then fieldVal = cutout.loc_metal
    if colorField eq 'vrad'      then fieldVal = reform(cutout.loc_vrad)
    if colorField eq 'vradnorm'  then fieldVal = reform(cutout.loc_vrad)/cutout.haloV200
    if colorField eq 'coolTime'  then fieldVal = cutout.loc_coolTime
    if colorField eq 'dynTime'   then fieldVal = cutout.loc_dynTime
    if colorField eq 'timeRatio' then fieldVal = cutout.loc_coolTime / cutout.loc_dynTime
    
    colorinds = (fieldVal-fieldMinMax[0])*(255.0-nbottom) / (fieldMinMax[1]-fieldMinMax[0])
    colorinds = fix(colorinds + nbottom) > 0 < 255 ;nbottom-255  
  
    ; second/right panel cutout
    wSecond = where(fieldVal le secondCutVal,nCutoutSecond,comp=wComp)
    if nCutoutSecond eq 0 then message,'None in cutout.'
      
    ; show gas above this cut value (instead of below)?
    if secondGt then wSecond = wComp
    
    loc_pos_second   = cutout.loc_pos[*,wSecond]
    loc_pos2_second  = cutout.loc_pos2[*,wSecond]
    colorinds_second = colorinds[wSecond]
    
    ; use instead a differently scaled color mapping for the second panel?
    if singleColorScale eq 0 then begin
      if secondGt eq 1 then $
        colorinds_second = (fieldVal-secondCutVal)*(255.0-nbottom) / (fieldMinMax[1]-secondCutVal)
      if secondGt eq 0 then $
        colorinds_second = (fieldVal-fieldMinMax[0])*(255.0-nbottom) / (secondCutVal-fieldMinMax[0])
      
      colorinds_second = fix(colorinds_second[wSecond] + nbottom) > 0 < 255 ; nbottom-255
    endif

    ; make a plot for each requested projection direction
    foreach axisPair, axes do begin
           
      ; get box center (in terms of specified axes)
      boxCenImg  = [sgcen[axisPair[0],gcID],sgcen[axisPair[1],gcID],sgcen[3-axisPair[0]-axisPair[1],gcID]]
      
      print,'['+string(gcID,format='(i4)')+'] Mapping ['+str(axisPair[0])+' '+$
            str(axisPair[1])+'] with '+str(cutout.boxSizeImg[0])+$
            ' kpc box around center ['+str(boxCenImg[0])+' '+str(boxCenImg[1])+' '+str(boxCenImg[2])+']'

      pFilename = 'scatter.'+sP.saveTag+'.'+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+$
                     '.axes'+str(axisPair[0])+str(axisPair[1])+'-'+$
                     colorField+'-'+secondField+'-'+str(secondGt)+'sCS'+str(singleColorScale)+'.eps'

      config = {boxSizeImg:cutout.boxSizeImg,plotFilename:pFilename,haloVirRad:cutout.haloVirRad,$
                haloMass:cutout.haloMass,axisPair:axisPair,sP:sP,$
                colorField:colorField,fieldMinMax:fieldMinMax,secondCutVal:secondCutVal,$
                secondGt:secondGt,nbottom:nbottom,barMM:fieldMinMax,ctName:ctName,barType:'2bar'}
      
      ; plot
      plotScatterComp,cutout.loc_pos,cutout.loc_pos2,loc_pos_second,loc_pos2_second,$
        colorinds,colorinds_second,config=config            

    endforeach ;axisPair

  endforeach ;gcIDs
  
  stop

end

; scatterMap4Panels(): four slices of some field for one halo
; showGmem/showAll (use galaxy catalog and display gmem only, or all in spatial subset)

pro scatterMap4Panels, showGmem=showGmem, showAll=showAll

  compile_opt idl2, hidden, strictarr, strictarrsubs

  sP = simParams(res=512,run='tracer',redshift=2.0)
  
  haloID = 304 ;z2.304 z2.301 z2.130 z2.64
  gcID = getMatchedIDs(sPa=sP,sPg=sP,haloID=haloID)
  gcIDs = [gcID.a]

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  singleColorScale = 1 ; 1=use same color scale for all panels, 0=rescale
  secondGt         = 0 ; 1=show greater than cut, 0=show less than cut
  axes             = list([0,1]) ;xy,xz,yz
  sizeFac          = 1.5 ; times rvir
  
  ; use which field and minmax for color mapping?
  ; temp [log K]
  ;colorField  = 'temp'
  ;fieldMinMax = [4.0,7.0]
  ;fieldRanges = list([4.5,5.0],[5.0,5.5],[5.5,6.0],[6.0,7.0])
  ;subtitles   = ['< 5.0','5.0 - 5.5','5.5 - 6.0','> 6.0']
  ;nbottom     = 50
  ;ctName      = 'helix' 
  
  ; vrad [km/s]
  colorField  = 'vrad'
  fieldMinMax = [-400,400]
  fieldRanges = list([-400,400],[-50,50],[-400,-300],[100,400])
  subtitles   = ['all','zero','inflow','outflow']
  nbottom     = 0
  ctName      = 'brewer-redblue'
  
  ; vrad norm [vrad/v200]
  ;colorField  = 'vradnorm'
  ;fieldMinMax = [0.0,4.0]
  ;fieldRanges = list([-5,3],[-1,1],[-5,-2],[2,3])
  ;subtitles   = ['all','zero','inflow,'outflow']
  ;nbottom     = 0
  ;ctName      = 'brewer-redblue'
  
  ; make cutouts (multiple or single)
  ;if keyword_set(showAll)  then cutout = cosmoVisCutout(sP=sP,gcInd=gcIDs,/selectGmem)
  ;if keyword_set(showAll)  then cutout = cosmoVisCutout(sP=sP,gcInd=gcIDs,sizeFac=sizeFac)

  ; target list
  gc    = loadGroupCat(sP=sP,/skipIDs)
  sgcen = subgroupPosByMostBoundID(sP=sP) 
  
  print,'rendering...'
  ; loop over all requested halos and image
  foreach gcID, gcIDs do begin
  
    ; load cutout
    if keyword_set(showAll)  then cutout = cosmoVisCutout(sP=sP,gcInd=gcID,/selectGmem)
    if keyword_set(showAll)  then cutout = cosmoVisCutout(sP=sP,gcInd=gcID,sizeFac=sizeFac)

    ; create color index mapping
    if colorField eq 'temp'      then fieldVal = cutout.loc_temp
    if colorField eq 'vrad'      then fieldVal = reform(cutout.loc_vrad)
    if colorField eq 'vradnorm'  then fieldVal = reform(cutout.loc_vrad)/cutout.haloV200
    
    ; make a plot for each requested projection direction
    foreach axisPair, axes do begin    
    
      pFilename = 'scatter4.'+sP.saveTag+'.'+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+$
                     '.axes'+str(axisPair[0])+str(axisPair[1])+'-'+$
                     colorField+'-'+str(secondGt)+'sCS'+str(singleColorScale)+'.eps'
      
      config = {boxSizeImg:cutout.boxSizeImg*sizeFac,plotFilename:pFilename,haloVirRad:cutout.haloVirRad,$
                haloMass:cutout.haloMass,axisPair:axisPair,sP:sP,barMM:fieldMinMax,$
                colorField:colorField,fieldMinMax:fieldMinMax,$
                ctName:ctName,nbottom:nbottom,barType:'1bar'}

      if ~singleColorScale then config.barType = ''
                
      start_PS, sP.plotPath + config.plotFilename, xs=8, ys=8
                
      ; cutouts and plot (k=0 top, k=1 bottom)
      for k=0,1 do begin
        wCutout = where(fieldVal ge fieldRanges[k*2+0,0] and fieldVal lt fieldRanges[k*2+0,1],nCutout)
          
        loc_pos_left = cutout.loc_pos[*,wCutout]
        loc_pos2_left = cutout.loc_pos2[*,wCutout]
        
        if singleColorScale then $
          colorinds = (fieldVal[wCutout]-fieldMinMax[0])*(255.0-nbottom) / $
                      (fieldMinMax[1]-fieldMinMax[0])
        if ~singleColorScale then $
          colorinds = (fieldVal[wCutout]-fieldRanges[k*2+0,0])*(255.0-nbottom) / $
                      (fieldRanges[k*2+0,1]-fieldRanges[k*2+0,0])
          
        colorinds_left = fix(colorinds + nbottom) > 0 < 255 ;nbottom-255  
        
        wCutout = where(fieldVal ge fieldRanges[k*2+1,0] and fieldVal lt fieldRanges[k*2+1,1],nCutout)
          
        loc_pos_right = cutout.loc_pos[*,wCutout]
        loc_pos2_right = cutout.loc_pos2[*,wCutout]
        
        if singleColorScale then $
          colorinds = (fieldVal[wCutout]-fieldMinMax[0])*(255.0-nbottom) / $
                      (fieldMinMax[1]-fieldMinMax[0])
        if ~singleColorScale then $
          colorinds = (fieldVal[wCutout]-fieldRanges[k*2+1,0])*(255.0-nbottom) / $
                      (fieldRanges[k*2+1,1]-fieldRanges[k*2+1,0])
          
        colorinds_right = fix(colorinds + nbottom) > 0 < 255 ;nbottom-255
        
        
        plotScatterComp,loc_pos_left,loc_pos2_left,loc_pos_right,loc_pos2_right,$
          colorinds_left,colorinds_right,config=config,$
          top=(k eq 0),bottom=(k eq 1),subtitle=subtitles[(k*2+0):(k*2+1)]
        
      endfor
      
      end_PS, pngResize=60  
  
    endforeach ;axisPair  
  
  endforeach ;gcIDs  
  stop
end

; mosaicHalosComp(): mosaic 4x2 comparison between arepo/gadget (cold only)

pro mosaicHalosComp, redshift=redshift

  compile_opt idl2, hidden, strictarr, strictarrsubs

  sPg = simParams(res=256,run='gadget',redshift=redshift)
  sPa = simParams(res=256,run='tracer',redshift=redshift)
  
  ; get list of matched IDs for good comparison
  mag = getMatchedIDs(sPa=sPa,sPg=sPg,/mosaicIDs)
  
  ; config
  sizeFac     = 2.5       ; times rvir
  tempMinMax  = [4.0,5.0] ; log(K)
  coldTempCut = 5.0       ; log(K)
  
  ; pre-cutout
  cutout = cosmoVisCutout(sP=sPg,gcInd=mag.gcIDsG,sizeFac=sizeFac)
  cutout = cosmoVisCutout(sP=sPa,gcInd=mag.gcIDsA,sizeFac=sizeFac)
  
  ; GADGET
  ; ------
  gc    = loadGroupCat(sP=sPg,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sPg) 

  ; start plot
  plotFilename = 'scatter.mosaic.'+sPg.saveTag+'.'+str(sPg.res)+'-'+sPa.saveTag+str(sPa.res)+$
                 '.'+str(sPg.snap)+'.eps'
                 
  start_PS, sPg.plotPath + plotFilename, xs=8, ys=8
  
  ; color table and establish temperature -> color mapping
  cgColorfill,[1,1,0,0,1],[1,0,0,1,1],/normal,color=cgColor('black')
  loadColorTable,'helix'
  
  ; initial position: upper left corner
  pos = [0.0,0.75,0.25,1.0]                 
                 
  gaHaloMasses = []
                 
  foreach gcIDg,mag.gcIDsG,k do begin
    ; load
    cutout = cosmoVisCutout(sP=sPg,gcInd=gcIDg,sizeFac=sizeFac)
    
    ; create color index mapping
    colorinds = (cutout.loc_temp-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
    colorinds = fix(colorinds + 50.0) > 0 < 255 ;50-255  
  
    ; local (cold) cutout
    wCold = where(cutout.loc_temp le coldTempCut,nCutoutCold)
    
    loc_pos_cold   = cutout.loc_pos[*,wCold]
    loc_pos2_cold  = cutout.loc_pos2[*,wCold]
    colorinds_cold = colorinds[wCold]
    
    ; get box center (in terms of specified axes)
    axisPair   = (mag.axes)[k]
    boxCenImg  = [sgcen[axisPair[0],gcIDg],sgcen[axisPair[1],gcIDg],sgcen[3-axisPair[0]-axisPair[1],gcIDg]]

    gaHaloMasses = [gaHaloMasses,cutout.haloMass]
  
    config = {boxSizeImg:cutout.boxSizeImg,plotFilename:plotFilename,haloVirRad:cutout.haloVirRad,$
              colorField:'na',fieldMinMax:tempMinMax,secondCutVal:coldTempCut,secondGt:0,$
              haloMass:cutout.haloMass,axisPair:axisPair,sP:sPg,barMM:tempMinMax,barType:'na'}
    
    ; plot
    xMinMax = [-config.boxSizeImg[0]/2.0,config.boxSizeImg[0]/2.0]
    yMinMax = [-config.boxSizeImg[1]/2.0,config.boxSizeImg[1]/2.0]
    
    !p.thick = 1.0
    !p.charsize = 0.8

    ; cold gas / dark matter (right panel)
    cgPlot, /nodata, xMinMax, yMinMax, pos=pos, xs=5, ys=5, /noerase

    tvcircle,config.haloVirRad,0,0,cgColor('light gray'),thick=0.6,/data
    
    ; particle loop for velocity vector plotting (cold gas only)
    nCutoutRight = n_elements(loc_pos_cold[0,*])
    for i=0L,nCutoutRight-1 do $
      oplot,[loc_pos_cold[config.axisPair[0],i],loc_pos2_cold[config.axisPair[0],i]],$
             [loc_pos_cold[config.axisPair[1],i],loc_pos2_cold[config.axisPair[1],i]],$
             line=0,color=colorinds_cold[i]
    
    ; scale bar
    if ~keyword_set(bottom) then begin
      len = 100.0 ;ckpc
      cgText,mean([-config.boxSizeImg[0]/2.2,-config.boxSizeImg[0]/2.2+len]),config.boxSizeImg[0]/2.4,$
             string(len,format='(i3)')+' kpc',alignment=0.5,color=cgColor('light gray'),$
             charsize=!p.charsize-0.2
      cgPlot,[-config.boxSizeImg[0]/2.2,-config.boxSizeImg[0]/2.2+len],$
             [config.boxSizeImg[1]/2.1,config.boxSizeImg[1]/2.1],$
             color=cgColor('light gray'),thick=4.0,/overplot
    endif
    
    if keyword_set(bottom) then $
      cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=1.0,color=cgColor('dark gray'),/overplot
    
    ; redshift and halo mass
    cgText,config.boxSizeImg[0]/2.1,config.boxSizeImg[1]/2.3,alignment=1.0,$
      "M = "+string(config.haloMass,format='(f4.1)'),color=cgColor('light gray'),charsize=!p.charsize-0.2
    
    ; dividing lines
    cgPlot,[xMinMax[0],xMinMax[0]],yMinMax,line=0,thick=0.5,color=cgColor('dark gray'),/overplot
    if k gt 3 then $
      cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=0.5,color=cgColor('dark gray'),/overplot
        
    ; update position
    pos += [0.25,0.0,0.25,0.0] ; move one space to the right
    if k eq 3 then pos -= [1.0,0.25,1.0,0.25] ; move to start of second row
        
  endforeach

  ; AREPO
  ; -----
  gc    = loadGroupCat(sP=sPa,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sPa) 
  
  ; initial position: 3rd row left space
  pos = [0.0,0.25,0.25,0.5]                 
                 
  foreach gcIDa,mag.gcIDsA,k do begin
    ; load
    cutout = cosmoVisCutout(sP=sPa,gcInd=gcIDa,sizeFac=sizeFac)
    
    ; create color index mapping
    colorinds = (cutout.loc_temp-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
    colorinds = fix(colorinds + 50.0) > 0 < 255 ;50-255  
  
    ; local (cold) cutout
    wCold = where(cutout.loc_temp le coldTempCut,nCutoutCold)
    
    loc_pos_cold   = cutout.loc_pos[*,wCold]
    loc_pos2_cold  = cutout.loc_pos2[*,wCold]
    colorinds_cold = colorinds[wCold]
  
    ; get box center (in terms of specified axes)
    axisPair   = (mag.axes)[k]
    boxCenImg  = [sgcen[axisPair[0],gcIDa],sgcen[axisPair[1],gcIDa],sgcen[3-axisPair[0]-axisPair[1],gcIDa]]
  
    config = {boxSizeImg:cutout.boxSizeImg,plotFilename:plotFilename,haloVirRad:cutout.haloVirRad,$
              colorField:'na',fieldMinMax:tempMinMax,secondCutVal:coldTempCut,secondGt:0,$
              haloMass:cutout.haloMass,axisPair:axisPair,sP:sPa,barMM:tempMinMax,barType:'na'}
    
    ; plot
    xMinMax = [-config.boxSizeImg[0]/2.0,config.boxSizeImg[0]/2.0]
    yMinMax = [-config.boxSizeImg[1]/2.0,config.boxSizeImg[1]/2.0]
    
    !p.thick = 1.0
    !p.charsize = 0.8

    ; cold gas / dark matter (right panel)
    cgPlot, /nodata, xMinMax, yMinMax, pos=pos, xs=5, ys=5, /noerase

    tvcircle,config.haloVirRad,0,0,cgColor('light gray'),thick=0.6,/data
    
    ; particle loop for velocity vector plotting (cold gas only)
    nCutoutRight = n_elements(loc_pos_cold[0,*])
    for i=0L,nCutoutRight-1 do $
      oplot,[loc_pos_cold[config.axisPair[0],i],loc_pos2_cold[config.axisPair[0],i]],$
             [loc_pos_cold[config.axisPair[1],i],loc_pos2_cold[config.axisPair[1],i]],$
             line=0,color=colorinds_cold[i]
    
    ; scale bar
    if ~keyword_set(bottom) then begin
      len = 100.0 ;ckpc
      cgText,mean([-config.boxSizeImg[0]/2.2,-config.boxSizeImg[0]/2.2+len]),config.boxSizeImg[0]/2.4,$
             string(len,format='(i3)')+' kpc',alignment=0.5,color=cgColor('light gray'),$
             charsize=!p.charsize-0.2
      cgPlot,[-config.boxSizeImg[0]/2.2,-config.boxSizeImg[0]/2.2+len],$
             [config.boxSizeImg[1]/2.1,config.boxSizeImg[1]/2.1],$
             color=cgColor('light gray'),thick=4.0,/overplot
    endif
    
    if keyword_set(bottom) then $
      cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=1.0,color=cgColor('dark gray'),/overplot
    
    ; redshift and halo mass
    cgText,config.boxSizeImg[0]/2.1,config.boxSizeImg[1]/2.3,alignment=1.0,$
      "M = "+string(gaHaloMasses[k],format='(f4.1)'),color=cgColor('light gray'),charsize=!p.charsize-0.2
    
    ; dividing lines
    cgPlot,[xMinMax[0],xMinMax[0]],yMinMax,line=0,thick=0.5,color=cgColor('dark gray'),/overplot
    if k le 3 then $
      cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=0.5,color=cgColor('light gray'),/overplot
    if k gt 3 then $
      cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=0.5,color=cgColor('dark gray'),/overplot
        
    ; update position
    pos += [0.25,0.0,0.25,0.0] ; move one space to the right
    if k eq 3 then pos -= [1.0,0.25,1.0,0.25] ; move to start of second row
        
  endforeach
  
  ; colorbar(s) on bottom
  !x.thick = 1.0
  !y.thick = 1.0

  ; gas temperature one colorbar (centered)
  pos = [0.02,0.35,0.076,0.65]
  fac = 0.5
  pos *= [fac,1,fac,1]
  
  loadColorTable,'helix'
  cgColorbar,position=pos,divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
  cgText,0.5,0.0375*fac,textoidl("log T_{gas} [K]"),alignment=0.5,color=cgColor('black'),/normal
  
  ; left/right colorbar labels
  ;cgText,0.365,0.036*fac,str(fix(tempMinMax[0])),alignment=0.5,color=cgColor('black'),/normal
  ;cgText,0.635,0.036*fac,str(fix(tempMinMax[1])),alignment=0.5,color=cgColor('black'),/normal
  ;cgText,0.365,0.036*fac,'4',alignment=0.5,color=cgColor('black'),/normal
  ;cgText,0.635,0.036*fac,'5',alignment=0.5,color=cgColor('black'),/normal
  cgText,0.370,0.036*fac,string(tempMinMax[0],format='(f3.1)'),alignment=0.5,color=cgColor('black'),/normal
  cgText,0.627,0.036*fac,string(tempMinMax[1],format='(f3.1)'),alignment=0.5,color=cgColor('black'),/normal
  
  ; simname labels
  cgText,0.5,0.95,sPg.simName,charsize=!p.charsize+0.4,alignment=0.5,/normal,color=cgColor('white')
  cgText,0.5,0.45,sPa.simName,charsize=!p.charsize+0.4,alignment=0.5,/normal,color=cgColor('white')
          
  end_PS, pngResize=60
  
  stop
end
