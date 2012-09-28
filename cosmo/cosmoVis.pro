; cosmoVis.pro
; cosmological boxes - 2d visualization
; dnelson sep.2012

; cosmoVisCutout(): make a spatial cutout around a halo
;                   call with multiple gcInd's for one load and save cutouts
;                   call with one gcInd to return results

function cosmoVisCutout, sP=sP, gcInd=gcInd, sizeFac=sizeFac

  velVecFac = 0.01 ; times velocity (km/s) in plotted kpc

  ; check existence of requested saves if more than one halo
  saveFilenames = sP.derivPath + 'cutout.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + $
                  '.h' + str(gcInd) + '.sf' + str(fix(sizeFac*10)) + '.sav'
                  
  readFlag = 0
  foreach saveFilename,saveFilenames do if ~file_test(saveFilename) then readFlag = 1
  
  ; if single halo requested and save exists, load it
  if n_elements(saveFilenames) eq 1 then $
    if file_test(saveFilenames) then restore,saveFilename

  ; proceed with cutouts if at least one save is missing
  if readFlag then begin
  
  gc    = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sP)   
  
  ; load u,nelec and calculate temperature
  u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
  nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
  temp  = alog10(convertUtoTemp(u,nelec))
  u     = !NULL
  nelec = !NULL

  ; load gas positions and velocities
  pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  vel = loadSnapshotSubset(sP=sP,partType='gas',field='vel')
  
  ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
  print,'shuffling...'
  iseed = 424242L
  sort_inds = sort(randomu(iseed,n_elements(temp)))
  
  temp = temp[sort_inds]
  pos  = pos[*,sort_inds]
  vel  = vel[*,sort_inds]
  
  sort_inds = !NULL
  
  print,'cutout...'
  foreach gcIndCur,gcInd,k do begin
    ; get subhalo position and size of imaging box
    boxCen     = sgcen[*,gcIndCur]
    boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcIndCur]] / 10.0) * 10.0
    boxSizeImg = [boxSize,boxSize,boxSize] ; cube
  
    ; make conservative cutout greater than boxsize accounting for periodic (do cube not sphere)
    xDist = pos[0,*] - boxCen[0]
    yDist = pos[1,*] - boxCen[1]
    zDist = pos[2,*] - boxCen[2]
    
    correctPeriodicDistVecs, xDist, sP=sP
    correctPeriodicDistVecs, yDist, sP=sP
    correctPeriodicDistVecs, zDist, sP=sP
    
    rvir = gc.group_r_crit200[gc.subgroupGrNr[gcIndCur]]
    rad = reform(sqrt(xDist*xDist + yDist*yDist + zDist*zDist)) / rvir[0]
  
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
    
    haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcIndCur]] ;ckpc
    haloMass = codeMassToLogMsun(gc.subgroupMass[gcIndCur])
    
    ; fix halo mass if we're using old (x2 bug) catalogs
    if sP.run eq 'gadgetold' or sP.run eq 'arepo' then $
      haloMass = codeMassToLogMsun(0.5*gc.subgroupMass[gcIndCur])  
  
    ; save
    r = {loc_pos:loc_pos,loc_temp:loc_temp,loc_pos2:loc_pos2,sP:sP,gcID:gcIndCur,$
         sizeFac:sizeFac,boxSizeImg:boxSizeImg,haloVirRad:haloVirRad,haloMass:haloMass}
         
    saveFilename = sP.derivPath + 'cutout.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + $
                   '.h' + str(gcIndCur) + '.sf' + str(fix(sizeFac*10)) + '.sav'  
         
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
      
        !p.thick = 1.0
        !p.charsize = 0.8
      
        ; fill with black background
        if ~keyword_set(bottom) then $
          cgColorfill,[1,1,0,0,1],[1,0,0,1,1],/normal,color=cgColor('black')
      
        ; color table and establish temperature -> color mapping
        loadColorTable,'helix';,/reverse
        
        TVLCT, rr, gg, bb, /GET

        newcolors_left = getColor24([[rr[cinds_left]], [gg[cinds_left]], [bb[cinds_left]]])

        ; all gas (left panel)
        cgPlot, /nodata, xMinMax, yMinMax, pos=posLeft, xs=5, ys=5, /noerase
        
        ; circle at virial radius
        tvcircle,config.haloVirRad,0,0,cgColor('light gray'),thick=0.6,/data
        
        ; particle loop for velocity vector plotting
        nCutoutLeft = n_elements(pos_left[0,*])
        for i=0L,nCutoutLeft-1 do $
          oplot,[pos_left[config.axisPair[0],i],pos2_left[config.axisPair[0],i]],$
                 [pos_left[config.axisPair[1],i],pos2_left[config.axisPair[1],i]],$
                 line=0,color=newcolors_left[i]
                 
        newcolors_left = !NULL
        
        ; scale bar
        if ~keyword_set(bottom) then begin
          len = 250.0 ;ckpc
          cgText,mean([-config.boxSizeImg[0]/2.2,-config.boxSizeImg[0]/2.2+len]),config.boxSizeImg[0]/2.3,$
                 string(len,format='(i3)')+' ckpc',alignment=0.5,color=cgColor('light gray')
          cgPlot,[-config.boxSizeImg[0]/2.2,-config.boxSizeImg[0]/2.2+len],$
                 [config.boxSizeImg[1]/2.1,config.boxSizeImg[1]/2.1],$
                 color=cgColor('light gray'),thick=4.0,/overplot
        endif
        
        if keyword_set(bottom) then $
          cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=1.0,color=cgColor('dark gray'),/overplot
               
        ; cold gas / dark matter (right panel)
        cgPlot, /nodata, xMinMax, yMinMax, pos=posRight, xs=5, ys=5, /noerase

        tvcircle,config.haloVirRad,0,0,cgColor('light gray'),thick=0.6,/data
        
        newcolors_right = getColor24([[rr[cinds_right]], [gg[cinds_right]], [bb[cinds_right]]])
        
        ; particle loop for velocity vector plotting (cold gas only)
        nCutoutRight = n_elements(pos_right[0,*])
        for i=0L,nCutoutRight-1 do $
          oplot,[pos_right[config.axisPair[0],i],pos2_right[config.axisPair[0],i]],$
                 [pos_right[config.axisPair[1],i],pos2_right[config.axisPair[1],i]],$
                 line=0,color=newcolors_right[i]
                 
        newcolors_right = !NULL
        
        ; redshift and halo mass
        if ~keyword_set(bottom) then begin
          cgText,config.boxSizeImg[0]/2.1,config.boxSizeImg[1]/2.4,alignment=1.0,$
            "z = "+string(config.sP.redshift,format='(f3.1)'),color=cgColor('light gray')
          cgText,config.boxSizeImg[0]/2.1,config.boxSizeImg[1]/2.2,alignment=1.0,$
            "M = "+string(config.haloMass,format='(f4.1)'),color=cgColor('light gray')
        endif
        
        ; dividing line
        cgPlot,[xMinMax[0],xMinMax[0]],yMinMax,line=0,thick=1.0,color=cgColor('dark gray'),/overplot
        
        if keyword_set(bottom) then $
          cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=1.0,color=cgColor('dark gray'),/overplot
        
        ; colorbar(s) on bottom
        !x.thick = 1.0
        !y.thick = 1.0
        
        if ~keyword_set(top) then begin
        if config.barType eq '2tempvdisp' then begin
          ; temp and veldisp separate colorbars
          colorbar,position=[0.02,0.1,0.076,0.4],divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
          cgText,0.25,0.0375,textoidl("log T_{gas} [K]"),alignment=0.5,color=cgColor('black'),/normal
          cgText,0.115,0.036,'4',alignment=0.5,color=cgColor('black'),/normal
          cgText,0.385,0.036,'7',alignment=0.5,color=cgColor('black'),/normal
          
          colorbar,position=[0.02,0.6,0.076,0.9],divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
          cgText,0.75,0.0375,textoidl("\sigma_{vel} [km/s]"),alignment=0.5,color=cgColor('black'),/normal
          cgText,0.615,0.036,'0',alignment=0.5,color=cgColor('black'),/normal
          cgText,0.875,0.036,string(config.barMM_right[1],format='(i3)'),$
            alignment=0.5,color=cgColor('black'),/normal
        endif
        
        if config.barType eq '1temp' then begin
          ; gas temperature one colorbar (centered)
          pos = [0.02,0.35,0.076,0.65]
          fac = 1.0
          if keyword_set(bottom) then fac = 0.5
          pos *= [fac,1,fac,1]
          colorbar,position=pos,divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
          cgText,0.5,0.0375*fac,textoidl("log T_{gas} [K]"),alignment=0.5,color=cgColor('black'),/normal
          ;cgText,0.365,0.036*fac,'4',alignment=0.5,color=cgColor('black'),/normal
          ;cgText,0.635,0.036*fac,'7',alignment=0.5,color=cgColor('black'),/normal
        endif
        
        if config.barType eq '1overdens' then begin
          ; local overdensity one colorbar (centered)
          colorbar,position=[0.02,0.35,0.076,0.65],divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
          cgText,0.5,0.0375,textoidl("log \rho_{DM} / <\rho_{DM}>"),alignment=0.5,color=cgColor('black'),/normal
          cgText,0.365,0.036,string(config.barMM[0],format='(i3)'),alignment=0.5,color=cgColor('black'),/normal
          cgText,0.635,0.036,string(config.barMM[1],format='(i3)'),alignment=0.5,color=cgColor('black'),/normal
        endif
        
        if config.barType eq '2overdens' then begin
          ; local overdensity two colorsbars (separate ranges)
          colorbar,position=[0.02,0.1,0.076,0.4],divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
          cgText,0.25,0.0375,textoidl("log \rho_{DM} / <\rho_{DM}>"),alignment=0.5,color=cgColor('black'),/normal
          cgText,0.115,0.036,string(config.barMM_left[0],format='(i3)'),alignment=0.5,color=cgColor('black'),/normal
          cgText,0.385,0.036,string(config.barMM_left[1],format='(i3)'),alignment=0.5,color=cgColor('black'),/normal

          colorbar,position=[0.02,0.6,0.076,0.9],divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
          cgText,0.75,0.0375,textoidl("log \rho_{DM} / <\rho_{DM}>"),alignment=0.5,color=cgColor('black'),/normal
          cgText,0.615,0.036,string(config.barMM_right[0],format='(i3)'),alignment=0.5,color=cgColor('black'),/normal
          cgText,0.875,0.036,string(config.barMM_right[1],format='(i3)'),alignment=0.5,color=cgColor('black'),/normal
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
        end_PS, pngResize=60;, im_options='-negate';, /deletePS
end

; scatterMapHalos: plot temperature colored scatter plots with velocity vectors on boxes centered on halos

pro scatterMapHalos, sP=sP, gcIDs=gcIDs

  sP = simParams(res=512,run='gadget',redshift=2.0)
  gcIDs = [5498] ;z2.304 g2342 a2132 ;z2.301 g2289 a2034 ;z2.64 g5498 a5097

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  sizeFac     = 3.5       ; times rvir
  tempMinMax  = [4.0,7.0] ; log(K)
  coldTempCut = 5.0       ; log(K)
  
  axes = list([0,1],[0,2],[1,2]) ;xy,xz,yz
  
  ; make cutouts (multiple or single)
  cutout = cosmoVisCutout(sP=sP,gcInd=gcIDs,sizeFac=sizeFac)
  
  ; target list
  gc    = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sP) 
  
  print,'rendering...'
  ; loop over all requested halos and image
  foreach gcID, gcIDs do begin
  
    ; load cutout
    cutout = cosmoVisCutout(sP=sP,gcInd=gcID,sizeFac=sizeFac)
  
    ; create color index mapping
    colorinds = (cutout.loc_temp-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
    colorinds = fix(colorinds + 50.0) > 0 < 255 ;50-255  
  
    ; local (cold) cutout
    wCold = where(cutout.loc_temp le coldTempCut,nCutoutCold)
    print,nCutout,nCutoutCold
    
    loc_pos_cold   = cutout.loc_pos[*,wCold]
    loc_pos2_cold  = cutout.loc_pos2[*,wCold]
    colorinds_cold = colorinds[wCold]
  
    ; make a plot for each requested projection direction
    foreach axisPair, axes do begin
           
      ; get box center (in terms of specified axes)
      boxCenImg  = [sgcen[axisPair[0],gcID],sgcen[axisPair[1],gcID],sgcen[3-axisPair[0]-axisPair[1],gcID]]
      
      print,'['+string(gcID,format='(i4)')+'] Mapping ['+str(axisPair[0])+' '+$
            str(axisPair[1])+'] with '+str(boxSize[0])+$
            ' kpc box around subhalo center ['+str(boxCen[0])+' '+str(boxCen[1])+' '+str(boxCen[2])+']'

      plotFilename = 'scatter.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+$
                     '.axes'+str(axisPair[0])+str(axisPair[1])+'.eps'

      config = {boxSizeImg:boxSizeImg,plotFilename:plotFilename,haloVirRad:cutout.haloVirRad,$
                haloMass:cutout.haloMass,axisPair:axisPair,sP:sP,barMM:tempMinMax,barType:'none'} ;'1temp'
      
      ; plot
      plotScatterComp,cutout.loc_pos,cutout.loc_pos2,loc_pos_cold,loc_pos2_cold,$
        colorinds,colorinds_cold,config=config            

    endforeach ;axisPair

  endforeach ;gcIDs

end

; scatterMapTemp4Panels(): four slices of temperature for one halo

pro scatterMapTemp4Panels

  compile_opt idl2, hidden, strictarr, strictarrsubs

  sPg = simParams(res=512,run='gadget',redshift=2.0)
  sPa = simParams(res=512,run='tracer',redshift=2.0)
  
  haloID = 304 ;z2.304 z2.301 z2.130 z2.64
  gcID = getMatchedIDs(sPa=sPa,sPg=sPg,haloID=haloID)

  ; config
  sizeFac     = 3.5       ; times rvir
  tempMinMax  = [4.0,7.0] ; log(K)
  tempRanges  = list([4.0,4.5],[4.5,5.0],[5.0,5.5],[5.5,7.0])       ; log(K)
  singleColorScale = 1
  
  ; select gadget or arepo
  gcInd = gcID.a
  sP = sPa
  
  ; load
  gc    = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sP)

  cutout = cosmoVisCutout(sP=sP,gcInd=gcInd,sizeFac=sizeFac)

  ; get box center (in terms of specified axes)
  boxCenImg  = [sgcen[gcID.axes[0],gcInd],$
                sgcen[gcID.axes[1],gcInd],$
                sgcen[3-gcID.axes[0]-gcID.axes[1],gcInd]]

  plotFilename = 'scatter4temp.'+sP.savPrefix+str(sP.res)+'.h'+str(gcInd)+$
                 '.'+str(sP.snap)+'.axes'+str(gcID.axes[0])+str(gcID.axes[1])+'.eps'

  config = {boxSizeImg:cutout.boxSizeImg,plotFilename:plotFilename,haloVirRad:cutout.haloVirRad,$
            haloMass:cutout.haloMass,axisPair:gcID.axes,sP:sP,barMM:tempMinMax,barType:'1temp'}

  ; temperature range cutouts and plot (top)
  wCutout = where(cutout.loc_temp ge tempRanges[0,0] and cutout.loc_temp lt tempRanges[0,1],nCutout)
    
  loc_pos_left = cutout.loc_pos[*,wCutout]
  loc_pos2_left = cutout.loc_pos2[*,wCutout]
  
  if singleColorScale then $
    colorinds = (cutout.loc_temp[wCutout]-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
  if ~singleColorScale then $
    colorinds = (cutout.loc_temp[wCutout]-tempRanges[0,0])*205.0 / (tempRanges[0,1]-tempRanges[0,0]) ;0-205
    
  colorinds_left = fix(colorinds + 50.0) > 0 < 255 ;50-255  
  
  wCutout = where(cutout.loc_temp ge tempRanges[1,0] and cutout.loc_temp lt tempRanges[1,1],nCutout)
    
  loc_pos_right = cutout.loc_pos[*,wCutout]
  loc_pos2_right = cutout.loc_pos2[*,wCutout]
  
  if singleColorScale then $
    colorinds = (cutout.loc_temp[wCutout]-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
  if ~singleColorScale then $
    colorinds = (cutout.loc_temp[wCutout]-tempRanges[1,0])*205.0 / (tempRanges[1,1]-tempRanges[1,0]) ;0-205
    
  colorinds_right = fix(colorinds + 50.0) > 0 < 255 ;50-255  
  
  start_PS, sP.plotPath + config.plotFilename, xs=8, ys=8
  
  plotScatterComp,loc_pos_left,loc_pos2_left,loc_pos_right,loc_pos2_right,$
    colorinds_left,colorinds_right,config=config,/top,subtitle=['4.0 - 4.5','4.5 - 5.0']           
  
  ; temperature range cutouts and plot (bottom)
  wCutout = where(cutout.loc_temp ge tempRanges[2,0] and cutout.loc_temp lt tempRanges[2,1],nCutout)
    
  loc_pos_left = cutout.loc_pos[*,wCutout]
  loc_pos2_left = cutout.loc_pos2[*,wCutout]
  
  if singleColorScale then $
    colorinds = (cutout.loc_temp[wCutout]-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
  if ~singleColorScale then $
    colorinds = (cutout.loc_temp[wCutout]-tempRanges[2,0])*205.0 / (tempRanges[2,1]-tempRanges[2,0]) ;0-205
    
  colorinds_left = fix(colorinds + 50.0) > 0 < 255 ;50-255  
  
  wCutout = where(cutout.loc_temp ge tempRanges[3,0] and cutout.loc_temp lt tempRanges[3,1],nCutout)
    
  loc_pos_right = cutout.loc_pos[*,wCutout]
  loc_pos2_right = cutout.loc_pos2[*,wCutout]
  
  if singleColorScale then $
    colorinds = (cutout.loc_temp[wCutout]-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
  if ~singleColorScale then $
    colorinds = (cutout.loc_temp[wCutout]-tempRanges[3,0])*205.0 / (tempRanges[3,1]-tempRanges[3,0]) ;0-205
    
  colorinds_right = fix(colorinds + 50.0) > 0 < 255 ;50-255  
  
  plotScatterComp,loc_pos_left,loc_pos2_left,loc_pos_right,loc_pos2_right,$
    colorinds_left,colorinds_right,config=config,/bottom,subtitle=['5.0 - 5.5','5.5 - 7.0'] 
  
  end_PS, pngResize=60  
  
end

; scatterMapHalosComp(): same as above but compare arepo/gadget top/bottom

pro scatterMapHalosComp

  compile_opt idl2, hidden, strictarr, strictarrsubs

  sPg = simParams(res=512,run='gadget',redshift=2.0)
  sPa = simParams(res=512,run='tracer',redshift=2.0)
  
  haloID = 130 ;z2.304 z2.301 z2.130 z2.64
  gcID = getMatchedIDs(sPa=sPa,sPg=sPg,haloID=haloID)

  ; config
  sizeFac     = 3.5       ; times rvir
  tempMinMax  = [4.0,7.0] ; log(K)
  coldTempCut = 5.0       ; log(K)
  
  ; GADGET
  ; ------
  gc    = loadGroupCat(sP=sPg,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sPg)

  cutout = cosmoVisCutout(sP=sPg,gcInd=gcID.g,sizeFac=sizeFac)

  ; create color index mapping
  colorinds = (cutout.loc_temp-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
  colorinds = fix(colorinds + 50.0) > 0 < 255 ;50-255  

  ; local (cold) cutout
  wCold = where(cutout.loc_temp le coldTempCut,nCutoutCold)
  
  loc_pos_cold   = cutout.loc_pos[*,wCold]
  loc_pos2_cold  = cutout.loc_pos2[*,wCold]
  colorinds_cold = colorinds[wCold]

  print,'rendering gadget...'
  
  ; get box center (in terms of specified axes)
  boxCenImg  = [sgcen[gcID.axes[0],gcID.g],$
                sgcen[gcID.axes[1],gcID.g],$
                sgcen[3-gcID.axes[0]-gcID.axes[1],gcID.g]]

  plotFilename = 'scatter.'+sPg.savPrefix+str(sPg.res)+'.h'+str(gcID.g)+'-'+sPa.savPrefix+str(sPa.res)+$
                 '.h'+str(gcID.a)+'.'+str(sPg.snap)+'.axes'+str(gcID.axes[0])+str(gcID.axes[1])+'.eps'

  config = {boxSizeImg:cutout.boxSizeImg,plotFilename:plotFilename,haloVirRad:cutout.haloVirRad,$
            haloMass:cutout.haloMass,axisPair:gcID.axes,sP:sPg,barMM:tempMinMax,barType:'1temp'}
  
  ; plot
  start_PS, sPg.plotPath + config.plotFilename, xs=8, ys=8
  plotScatterComp,cutout.loc_pos,cutout.loc_pos2,loc_pos_cold,loc_pos2_cold,$
  colorinds,colorinds_cold,config=config,/top           

  ; AREPO
  ; -----
  gc    = loadGroupCat(sP=sPa,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sPa) 
  
  cutout = cosmoVisCutout(sP=sPa,gcInd=gcID.a,sizeFac=sizeFac)
  
  ; create color index mapping
  colorinds = (cutout.loc_temp-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
  colorinds = fix(colorinds + 50.0) > 0 < 255 ;50-255  

  ; local (cold) cutout
  wCold = where(cutout.loc_temp le coldTempCut,nCutoutCold)
  
  loc_pos_cold   = cutout.loc_pos[*,wCold]
  loc_pos2_cold  = cutout.loc_pos2[*,wCold]
  colorinds_cold = colorinds[wCold]

  print,'rendering arepo...'
  
  ; get box center (in terms of specified axes)
  boxCenImg  = [sgcen[gcID.axes[0],gcID.a],$
                sgcen[gcID.axes[1],gcID.a],$
                sgcen[3-gcID.axes[0]-gcID.axes[1],gcID.a]]

  config = {boxSizeImg:cutout.boxSizeImg,plotFilename:plotFilename,haloVirRad:cutout.haloVirRad,$
            haloMass:cutout.haloMass,axisPair:gcID.axes,sP:sPa,barMM:tempMinMax,barType:'1temp'}
  
  ; plot
  plotScatterComp,cutout.loc_pos,cutout.loc_pos2,loc_pos_cold,loc_pos2_cold,$
    colorinds,colorinds_cold,config=config,/bottom
  end_PS, pngResize=60
  
end

; scatterMapPastPosComp(): plot past positions of gas segregated by hot/cold (compare arepo/gadget top/bottom)

pro scatterMapPastPosComp

  compile_opt idl2, hidden, strictarr, strictarrsubs

  sPg = simParams(res=512,run='gadget',redshift=2.0)
  sPa = simParams(res=128,run='tracer',redshift=2.0)
  print,'change tracer 256 to 512'
  print,'128!'
  
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
    accTvir = gcSubsetProp(sP=sPg,select='pri',/accTvir,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
    maxTemp = gcSubsetProp(sP=sPg,select='pri',/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
    
    accTvir = accTvir.gal[ww]
    maxTemp = maxTemp.gal[ww]
    ratio   = 10.0^maxTemp / 10.0^accTvir
    
    ; load galcat and make gas ID list
    galIDs = gcSubsetProp(sP=sPg,select='pri',/elemIDs,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
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
    match,ids,galIDs.hot,ids_ind_hot,inds2,count=countHot
    match,ids,galIDs.cold,ids_ind_cold,inds2,count=countCold
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
    accTvir = gcSubsetProp(sP=sPa,select='pri',/accTvir,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
    maxTemp = gcSubsetProp(sP=sPa,select='pri',/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)

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
    match,ids,galIDs.hot,ids_ind_hot,inds2,count=countHot
    match,ids,galIDs.cold,ids_ind_cold,inds2,count=countCold
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

; mosaicHalosComp(): mosaic 4x2 comparison between arepo/gadget (cold only)

pro mosaicHalosComp

  compile_opt idl2, hidden, strictarr, strictarrsubs

  sPg = simParams(res=512,run='gadgetold',redshift=2.0)
  sPa = simParams(res=512,run='arepo',redshift=2.0)
  
  ; good candidates:
  ; z2.10  g4467 a4103 axes12 -- z2.40  g5151 z4789 axes01 -- z2.130 g6369 a5966 axes02
  ; z2.167 g3338 a3022 axes02 -- z2.251 g2538 a2180 axes12 -- z2.262 g2824 a2513 axes02
  ; z2.314 g981  a927  axes12 -- z2.315 g1117 a966  axes02 -- z2.291 g1478 a1341 axes02
  ; z2.322 g266  a257  axes01
  
  gcIDsG = [6369,5151,3338,2824,2538,1117,981,266]
  gcIDsA = [5966,4789,3022,2513,2180,966,927,257]
  axes   = list([0,2],[0,1],[0,2],[0,2],[1,2],[0,2],[1,2],[0,1])
  
  ; config
  sizeFac     = 2.5       ; times rvir
  tempMinMax  = [4.0,7.0] ; log(K)
  coldTempCut = 5.0       ; log(K)
  
  ; GADGET
  ; ------
  gc    = loadGroupCat(sP=sPg,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sPg) 

  ; start plot
  plotFilename = 'scatter.mosaic.'+sPg.savPrefix+str(sPg.res)+'-'+sPa.savPrefix+str(sPa.res)+$
                 '.'+str(sPg.snap)+'.eps'
                 
  start_PS, sPg.plotPath + plotFilename, xs=8, ys=8
  
  ; color table and establish temperature -> color mapping
  cgColorfill,[1,1,0,0,1],[1,0,0,1,1],/normal,color=cgColor('black')
  loadColorTable,'helix'
  
  TVLCT, rr, gg, bb, /GET
  
  textColor = getColor24([180,180,180]) ; light gray
  
  ; initial position: upper left corner
  pos = [0.0,0.75,0.25,1.0]                 
                 
  gaHaloMasses = []
                 
  foreach gcIDg,gcIDsG,k do begin
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
    axisPair   = axes[k]
    boxCenImg  = [sgcen[axisPair[0],gcIDg],sgcen[axisPair[1],gcIDg],sgcen[3-axisPair[0]-axisPair[1],gcIDg]]

    gaHaloMasses = [gaHaloMasses,cutout.haloMass]
  
    config = {boxSizeImg:boxSizeImg,plotFilename:plotFilename,haloVirRad:cutout.haloVirRad,$
              haloMass:cutout.haloMass,axisPair:axisPair,sP:sPg,barMM:tempMinMax,barType:'1temp'}
    
    ; plot
    xMinMax = [-config.boxSizeImg[0]/2.0,config.boxSizeImg[0]/2.0]
    yMinMax = [-config.boxSizeImg[1]/2.0,config.boxSizeImg[1]/2.0]
    
    !p.thick = 1.0
    !p.charsize = 0.8

    ; cold gas / dark matter (right panel)
    cgPlot, /nodata, xMinMax, yMinMax, pos=pos, xs=5, ys=5, /noerase

    tvcircle,config.haloVirRad,0,0,cgColor('light gray'),thick=0.6,/data
    
    newcolors_right = getColor24([[rr[colorinds_cold]], [gg[colorinds_cold]], [bb[colorinds_cold]]])
           
    ; particle loop for velocity vector plotting (cold gas only)
    nCutoutRight = n_elements(loc_pos_cold[0,*])
    for i=0L,nCutoutRight-1 do $
      oplot,[loc_pos_cold[config.axisPair[0],i],loc_pos2_cold[config.axisPair[0],i]],$
             [loc_pos_cold[config.axisPair[1],i],loc_pos2_cold[config.axisPair[1],i]],$
             line=0,color=newcolors_right[i]
             
    newcolors_right = !NULL
    
    ; scale bar
    if ~keyword_set(bottom) then begin
      len = 100.0 ;ckpc
      cgText,mean([-config.boxSizeImg[0]/2.2,-config.boxSizeImg[0]/2.2+len]),config.boxSizeImg[0]/2.4,$
             string(len,format='(i3)')+' ckpc',alignment=0.5,color=textColor,charsize=!p.charsize-0.2
      cgPlot,[-config.boxSizeImg[0]/2.2,-config.boxSizeImg[0]/2.2+len],$
             [config.boxSizeImg[1]/2.1,config.boxSizeImg[1]/2.1],$
             color=textColor,thick=4.0,/overplot
    endif
    
    if keyword_set(bottom) then $
      cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=1.0,color=cgColor('dark gray'),/overplot
    
    ; redshift and halo mass
    cgText,config.boxSizeImg[0]/2.1,config.boxSizeImg[1]/2.3,alignment=1.0,$
      "M = "+string(config.haloMass,format='(f4.1)'),color=textColor,charsize=!p.charsize-0.2
    
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
                 
  foreach gcIDa,gcIDsA,k do begin
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
    axisPair   = axes[k]
    boxCenImg  = [sgcen[axisPair[0],gcIDa],sgcen[axisPair[1],gcIDa],sgcen[3-axisPair[0]-axisPair[1],gcIDa]]
  
    config = {boxSizeImg:boxSizeImg,plotFilename:plotFilename,haloVirRad:cutout.haloVirRad,$
              haloMass:cutout.haloMass,axisPair:axisPair,sP:sPa,barMM:tempMinMax,barType:'1temp'}
    
    ; plot
    xMinMax = [-config.boxSizeImg[0]/2.0,config.boxSizeImg[0]/2.0]
    yMinMax = [-config.boxSizeImg[1]/2.0,config.boxSizeImg[1]/2.0]
    
    !p.thick = 1.0
    !p.charsize = 0.8

    ; cold gas / dark matter (right panel)
    cgPlot, /nodata, xMinMax, yMinMax, pos=pos, xs=5, ys=5, /noerase

    tvcircle,config.haloVirRad,0,0,cgColor('light gray'),thick=0.6,/data
    
    newcolors_right = getColor24([[rr[colorinds_cold]], [gg[colorinds_cold]], [bb[colorinds_cold]]])
           
    ; particle loop for velocity vector plotting (cold gas only)
    nCutoutRight = n_elements(loc_pos_cold[0,*])
    for i=0L,nCutoutRight-1 do $
      oplot,[loc_pos_cold[config.axisPair[0],i],loc_pos2_cold[config.axisPair[0],i]],$
             [loc_pos_cold[config.axisPair[1],i],loc_pos2_cold[config.axisPair[1],i]],$
             line=0,color=newcolors_right[i]
             
    newcolors_right = !NULL
    
    ; scale bar
    if ~keyword_set(bottom) then begin
      len = 100.0 ;ckpc
      cgText,mean([-config.boxSizeImg[0]/2.2,-config.boxSizeImg[0]/2.2+len]),config.boxSizeImg[0]/2.4,$
             string(len,format='(i3)')+' ckpc',alignment=0.5,color=textColor,charsize=!p.charsize-0.2
      cgPlot,[-config.boxSizeImg[0]/2.2,-config.boxSizeImg[0]/2.2+len],$
             [config.boxSizeImg[1]/2.1,config.boxSizeImg[1]/2.1],$
             color=textColor,thick=4.0,/overplot
    endif
    
    if keyword_set(bottom) then $
      cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=1.0,color=cgColor('dark gray'),/overplot
    
    ; redshift and halo mass
    cgText,config.boxSizeImg[0]/2.1,config.boxSizeImg[1]/2.3,alignment=1.0,$
      "M = "+string(gaHaloMasses[k],format='(f4.1)'),color=textColor,charsize=!p.charsize-0.2
    
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
  colorbar,position=pos,divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
  cgText,0.5,0.0375*fac,textoidl("log T_{gas} [K]"),alignment=0.5,color=cgColor('black'),/normal
  cgText,0.365,0.036*fac,'4',alignment=0.5,color=cgColor('black'),/normal
  cgText,0.635,0.036*fac,'7',alignment=0.5,color=cgColor('black'),/normal
  
  cgText,0.5,0.95,"GADGET",charsize=!p.charsize+0.4,alignment=0.5,/normal,color=cgColor('white')
  cgText,0.5,0.45,"AREPO",charsize=!p.charsize+0.4,alignment=0.5,/normal,color=cgColor('white')
          
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
