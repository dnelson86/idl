; cosmoVis.pro
; cosmological boxes - 2d visualization
; dnelson apr.2012

; sphMapBox: run sph kernel density projection on whole box

pro sphMapBox, res=res, run=run, partType=partType

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  ;res = 128
  ;run = 'dev.tracer.nocomov'
  ;partType = 'gas'

  ;redshift = 3.0 ;5.0
  ;snap     = redshiftToSnapNum(redshift,sP=sP)
  snap = 19
  
  nPixels = [800,800] ;px

  zoomFac = 1    ; only in axes, not along projection direction
  nNGB    = 64   ; use CalcHSML for HSML with nNGB
  axes    = [0,1] ; x,y

  ; paths and render config
  sP = simParams(res=res,run=run)
  h = loadSnapshotHeader(sP.simPath,snapNum=snap)
  
  boxSize = [h.boxSize,h.boxSize,h.boxSize]              ;kpc
  boxCen  = [h.boxSize/2.0,h.boxSize/2.0,h.boxSize/2.0]  ;kpc
  
  foreach k,axes do boxSize[k] /= zoomFac
  
  outFilename = 'sphmap.box_'+str(zoomFac)+'.nNGB='+str(nNGB)+'.snap='+str(snap)+$
                '.box.axis0='+str(axes[0])+'.axis1='+str(axes[0])+'.'+partType

  ; save/restore
  if (file_test(sP.derivPath + outFilename + '.sav')) then begin
    restore,sP.derivPath + outFilename + '.sav',/verbose
  endif else begin

    if (partType eq 'gas') then begin
      mass = loadSnapshotSubset(sP.simPath,snapNum=snap,partType='gas',field='mass',/verbose)
    endif
    
    if (partType eq 'tracer') then begin
       mass_gas = loadSnapshotSubset(sP.simPath,snapNum=snap,partType='gas',field='mass',/verbose)
       mass = replicate(total(mass_gas) / h.nPartTot[3], h.nPartTot[3])
    endif
    
    ; load positions from snapshot
    pos  = loadSnapshotSubset(sP.simPath,snapNum=snap,partType=partType,field='pos',/verbose) 
    
    ; load HSML from snapshot (only stored for gas)
    if partType eq 'gas' then $
      hsml = loadSnapshotSubset(sP.simPath,snapNum=snap,partType='gas',field='hsml',/verbose) $
    ; if different particle type, calculate HSMLs
    else hsml = calcHSML(pos,ndims=3,nNGB=nNGB,boxSize=boxSize[0])
    
    colMassMap = calcSphMap(pos,hsml,mass,boxSize=boxSize,boxCen=boxCen,nPixels=nPixels,$
                            axes=axes,ndims=3)
              
    save,colMassMap,hsml,filename=sP.derivPath + outFilename + '.sav'
  endelse

  ; rescale
  ;maxVal = max(colMassMap)/2.0
  maxVal = 0.5
  minVal = maxVal / 1e4
  
  print,'min max val: ',minVal,maxVal
  
  w = where(colMassMap eq 0, count, complement=ww)
  if (count ne 0) then colMassMap[w] = min(colMassMap[ww])
  
  colMassMap = colMassMap > minVal < maxVal
  
  colMassMap = alog10(colMassMap)
  colMassMap = (colMassMap-min(colMassMap))*254.0 / (max(colMassMap)-min(colMassMap)) ;0-254
  ;h2d = filter_image(h2d,FWHM=[1.1,1.1],/ALL) ;gaussian kernel convolution
  colMassMap += 1.0 ;1-255

  ; plot
  xMinMax = [boxCen[0]-boxSize[0]/2.0,boxCen[0]+boxSize[0]/2.0]
  yMinMax = [boxCen[1]-boxSize[1]/2.0,boxCen[1]+boxSize[1]/2.0]
  
  start_PS, sP.plotPath + outFilename + '.eps'
    loadct, 4, bottom=1, /silent
    tvim,colMassMap,xrange=xMinMax,yrange=yMinMax,pcharsize=0.0001
  end_PS, pngResize=68, /deletePS
  
end

; plotScatterComp(): plot side by side colored/vectorized scatter plots

pro plotScatterComp, pos_left, pos2_left, pos_right, pos2_right, cinds_left, cinds_right, config=config

      xMinMax = [-config.boxSizeImg[0]/2.0,config.boxSizeImg[0]/2.0]
      yMinMax = [-config.boxSizeImg[1]/2.0,config.boxSizeImg[1]/2.0]
      
      plotPath = '/n/home07/dnelson/data3/HaloComp/'
      
      start_PS, plotPath + config.plotFilename, xs=8, ys=4
      
        !p.thick = 1.0
        !p.charsize = 0.8
      
        ; color table and establish temperature -> color mapping
        loadColorTable,'helix',/reverse
        
        TVLCT, rr, gg, bb, /GET

        newcolors_left = getColor24([[rr[cinds_left]], [gg[cinds_left]], [bb[cinds_left]]])

        ; all gas
        cgPlot, /nodata, xMinMax, yMinMax, pos=[0.0,0.0,0.5,1.0], xs=5, ys=5
        
        ; circle at virial radius
        tvcircle,config.haloVirRad,0,0,cgColor('dark gray'),thick=0.6,/data
        
        ; particle loop for velocity vector plotting
        nCutoutLeft = n_elements(pos_left[0,*])
        for i=0L,nCutoutLeft-1 do $
          oplot,[pos_left[config.axisPair[0],i],pos2_left[config.axisPair[0],i]],$
                 [pos_left[config.axisPair[1],i],pos2_left[config.axisPair[1],i]],$
                 line=0,color=newcolors_left[i]
                 
        newcolors_left = !NULL
        
        ; scale bar
        len = 250.0 ;ckpc
        cgText,mean([-config.boxSizeImg[0]/2.2,-config.boxSizeImg[0]/2.2+len]),config.boxSizeImg[0]/2.3,$
               string(len,format='(i3)')+' ckpc',alignment=0.5,color=cgColor('dark gray')
        cgPlot,[-config.boxSizeImg[0]/2.2,-config.boxSizeImg[0]/2.2+len],$
               [config.boxSizeImg[1]/2.1,config.boxSizeImg[1]/2.1],$
               color=cgColor('dark gray'),thick=4.0,/overplot
               
        ; dark matter
        cgPlot, /nodata, xMinMax, yMinMax, pos=[0.5,0.0,1.0,1.0], xs=5, ys=5, /noerase

        tvcircle,config.haloVirRad,0,0,cgColor('dark gray'),thick=0.6,/data
        
        newcolors_right = getColor24([[rr[cinds_right]], [gg[cinds_right]], [bb[cinds_right]]])
        
        ; particle loop for velocity vector plotting (cold gas only)
        nCutoutRight = n_elements(pos_right[0,*])
        for i=0L,nCutoutRight-1 do $
          oplot,[pos_right[config.axisPair[0],i],pos2_right[config.axisPair[0],i]],$
                 [pos_right[config.axisPair[1],i],pos2_right[config.axisPair[1],i]],$
                 line=0,color=newcolors_right[i]
                 
        newcolors_right = !NULL
        
        ; redshift and halo mass
        cgText,0.99,0.96,"z = "+string(config.sP.redshift,format='(f3.1)'),alignment=1.0,color=cgColor('dark gray'),/normal
        cgText,0.99,0.92,"M = "+string(config.haloMass,format='(f4.1)'),alignment=1.0,color=cgColor('dark gray'),/normal
        
        ; dividing line
        cgPlot,[xMinMax[0],xMinMax[0]],yMinMax,line=0,thick=1.0,color=cgColor('light gray'),/overplot
        
        ; colorbar(s) on bottom
        !x.thick = 1.0
        !y.thick = 1.0
        
        if config.barType eq '2tempvdisp' then begin
          ; temp and veldisp separate colorbars
          colorbar,position=[0.02,0.1,0.076,0.4],divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
          cgText,0.25,0.0375,textoidl("log T_{gas} [K]"),alignment=0.5,color=cgColor('white'),/normal
          cgText,0.115,0.036,'4',alignment=0.5,color=cgColor('white'),/normal
          cgText,0.385,0.036,'7',alignment=0.5,color=cgColor('white'),/normal
          
          colorbar,position=[0.02,0.6,0.076,0.9],divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
          cgText,0.75,0.0375,textoidl("\sigma_{vel} [km/s]"),alignment=0.5,color=cgColor('white'),/normal
          cgText,0.615,0.036,'0',alignment=0.5,color=cgColor('white'),/normal
          cgText,0.875,0.036,string(config.barMM_right[1],format='(i3)'),$
            alignment=0.5,color=cgColor('white'),/normal
        endif
        
        if config.barType eq '1temp' then begin
          ; gas temperature one colorbar (centered)
          colorbar,position=[0.02,0.35,0.076,0.65],divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
          cgText,0.5,0.0375,textoidl("log T_{gas} [K]"),alignment=0.5,color=cgColor('white'),/normal
          cgText,0.365,0.036,'4',alignment=0.5,color=cgColor('white'),/normal
          cgText,0.635,0.036,'7',alignment=0.5,color=cgColor('white'),/normal
        endif
        
        if config.barType eq '1overdens' then begin
          ; local overdensity one colorbar (centered)
          colorbar,position=[0.02,0.35,0.076,0.65],divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
          cgText,0.5,0.0375,textoidl("log \rho_{DM} / <\rho_{DM}>"),alignment=0.5,color=cgColor('white'),/normal
          cgText,0.365,0.036,string(config.barMM[0],format='(i3)'),alignment=0.5,color=cgColor('white'),/normal
          cgText,0.635,0.036,string(config.barMM[1],format='(i3)'),alignment=0.5,color=cgColor('white'),/normal
        endif
        
      end_PS, pngResize=60, im_options='-negate', /deletePS
end

; scatterMapHalos: plot temperature colored scatter plots with velocity vectors on boxes centered on halos

pro scatterMapHalos, sP=sP, gcIDs=gcIDs

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  sizeFac = 3.5       ; times rvir
  cutFac  = 1.0       ; times boxSize
  nPixels = [800,800] ; px
  
  tempMinMax  = [4.0,7.0] ; log(K)
  coldTempCut = 5.0       ; log(K)
  velVecFac   = 0.01     ; times velocity (km/s) in plotted kpc
  
  axes = list([0,1],[0,2],[1,2]) ;xy,xz,yz
  
  ; target list
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
    wCut = where(abs(xDist) le 0.5*cutFac*boxSize and abs(yDist) le 0.5*cutFac*boxSize and $
                 abs(zDist) le 0.5*cutFac*boxSize,nCutout)
    
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
  
    ; local (cold) cutout
    wCold = where(loc_temp le coldTempCut,nCutoutCold)
    loc_temp = !NULL
    print,nCutout,nCutoutCold
    
    loc_pos_cold   = loc_pos[*,wCold]
    loc_pos2_cold  = loc_pos2[*,wCold]
    colorinds_cold = colorinds[wCold]
  
    ; make a plot for each requested projection direction
    foreach axisPair, axes do begin
           
      ; get box center (in terms of specified axes)
      boxCenImg  = [sgcen[axisPair[0],gcID],sgcen[axisPair[1],gcID],sgcen[3-axisPair[0]-axisPair[1],gcID]]
      haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcID]] ;ckpc
      haloMass = codeMassToLogMsun(gc.subgroupMass[gcID])
      
      print,'['+string(gcID,format='(i4)')+'] Mapping ['+str(axisPair[0])+' '+$
            str(axisPair[1])+'] with '+str(boxSize[0])+$
            ' kpc box around subhalo center ['+str(boxCen[0])+' '+str(boxCen[1])+' '+str(boxCen[2])+']'

      plotFilename = 'scatter.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+$
                     '.axes'+str(axisPair[0])+str(axisPair[1])+'.eps'

      config = {boxSizeImg:boxSizeImg,plotFilename:plotFilename,haloVirRad:haloVirRad,haloMass:haloMass,$
                axisPair:axisPair,sP:sP,barMM:tempMinMax,barType:'1temp'}
      
      ; plot
      plotScatterComp,loc_pos,loc_pos2,loc_pos_cold,loc_pos2_cold,colorinds,colorinds_cold,config=config            

    endforeach ;axisPair

  endforeach ;gcIDs

end

; scatterMapHalosGasDM: plot temperature colored scatter plots with velocity vectors on boxes centered on halos
;                       gas on left and DM on right (usually more zoomed out)

pro scatterMapHalosGasDM, sP=sP, gcIDs=gcIDs

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  sizeFac = 10.0      ; times rvir
  cutFac  = 1.0       ; times boxSize
  nPixels = [800,800] ; px
  
  tempMinMax  = [4.0,7.0] ; log(K)
  velVecFac   = 0.01      ; times velocity (km/s) in plotted kpc
  
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
    wCut = where(abs(xDist) le 0.5*cutFac*boxSize and abs(yDist) le 0.5*cutFac*boxSize and $
                 abs(zDist) le 0.5*cutFac*boxSize,nCutout)
    
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
    wCut = where(abs(xDist) le 0.5*cutFac*boxSize and abs(yDist) le 0.5*cutFac*boxSize and $
                 abs(zDist) le 0.5*cutFac*boxSize,nCutoutDM)
    
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
  sizeFac = 10.0      ; times rvir
  cutFac  = 1.0       ; times boxSize
  nPixels = [800,800] ; px
  
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
    ;loc_dens_dm = !NULL
    print,nCutoutDM,nOverDens
    
    loc_pos_od   = loc_pos_dm[*,wOD]
    loc_pos2_od  = loc_pos2_dm[*,wOD]
    colorinds_od = colorinds_dm[wOD]
  
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
                axisPair:axisPair,sP:sP,barMM:overdensMM,barType:'1overdens'}
      
      ; plot
      plotScatterComp,loc_pos_dm,loc_pos2_dm,loc_pos_od,loc_pos2_od,colorinds_dm,colorinds_od,config=config
stop
    endforeach ;axisPair

  endforeach ;gcIDs

end

; plotSphmapDensQuant(): plot side by side projection results from CalcSphMap

pro plotSphmapDensQuant, map=sphmap, config=config

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; color map and rescale
  w = where(sphmap.dens_out eq 0.0,count,comp=wc)
  if count gt 0 then sphmap.dens_out[w] = min(sphmap.dens_out[wc])

  sphmap.dens_out = alog10(sphmap.dens_out)
  sphmap.dens_out = (sphmap.dens_out-min(sphmap.dens_out))*254.0 / $
                    (max(sphmap.dens_out)-min(sphmap.dens_out)) ;0-254
  sphmap.dens_out += 1.0 ;1-255
  
  w = where(sphmap.quant_out eq 0.0,count,comp=wc)
  if count gt 0 then sphmap.quant_out[w] = min(sphmap.quant_out[wc])
  
  sphmap.quant_out = alog10(sphmap.quant_out)
  sphmap.quant_out = (sphmap.quant_out-min(sphmap.quant_out))*254.0 / $
                    (max(sphmap.quant_out)-min(sphmap.quant_out)) ;0-254
  sphmap.quant_out += 1.0 ;1-255

  xMinMax = [config.boxCen[0]-config.boxSizeImg[0]/2.0,config.boxCen[0]+config.boxSizeImg[0]/2.0]
  yMinMax = [config.boxCen[1]-config.boxSizeImg[1]/2.0,config.boxCen[1]+config.boxSizeImg[1]/2.0]
  
  plotPath = '/n/home07/dnelson/data3/HaloComp/'
  
  start_PS, plotPath + strmid(config.saveFilename,0,strlen(config.saveFilename)-4)+'.eps', xs=8, ys=4
  
    !p.charsize = 0.8
    loadColorTable,'helix'
    
    ; density
    cgPlot, /nodata, xMinMax, yMinMax, pos=[0.0,0.0,0.5,1.0], xs=5, ys=5
    tv, sphmap.dens_out,0.0,0.0,/normal,xsize=0.5
    
    ; circle at virial radius
    ;tvcircle,config.haloVirRad,0,0,cgColor('dark gray'),thick=0.6,/data
    
    ; scale bar
    len = 100.0 ;ckpc
    cgText,mean([config.boxCen[0]-config.boxSizeImg[0]/2.2,$
                 config.boxCen[0]-config.boxSizeImg[0]/2.2+len]),$
           config.boxCen[1]+config.boxSizeImg[0]/2.3,$
           string(len,format='(i3)')+' ckpc',alignment=0.5,color=cgColor('black')
    cgPlot,[config.boxCen[0]-config.boxSizeImg[0]/2.2,config.boxCen[0]-config.boxSizeImg[0]/2.2+len],$
           [config.boxCen[1]+config.boxSizeImg[1]/2.1,config.boxCen[1]+config.boxSizeImg[1]/2.1],$
           color=cgColor('black'),/overplot
           
    ; mass weighted temperature
    cgPlot, /nodata, xMinMax, yMinMax, pos=[0.0,0.0,0.5,1.0], xs=5, ys=5, /noerase
    tv, sphmap.quant_out,0.5,0.0,/normal,xsize=0.5
    
    ; redshift and halo mass
    cgText,0.99,0.96,"z = "+string(config.sP.redshift,format='(f3.1)'),alignment=1.0,$
           color=cgColor('black'),/normal
    cgText,0.99,0.92,"M = "+string(config.haloMass,format='(f4.1)'),alignment=1.0,$
           color=cgColor('black'),/normal
             
  end_PS, pngResize=60, /deletePS

end

; sphMapHalos: run sph kernel density projection on gas particles/cells with boxes centered on halos

pro sphMapHalos, sP=sP, gcIDs=gcIDs, coldOnly=coldOnly

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  sizeFac = 3.5       ; times rvir
  nPixels = [800,800] ; px
  
  axes = list([0,1],[0,2],[1,2]) ;xy,xz,yz
  
  if keyword_set(coldOnly) then coldFlag = '.cold' else coldFlag = ''
  
  ; target list
  gc    = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sP) 

  ; load u,nelec and calculate temperature
  u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
  nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
  temp  = convertUtoTemp(u,nelec)
  u     = !NULL
  nelec = !NULL

  ; load gas positions and masses from snapshot
  pos  = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  mass = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
  
  ; load HSMLs or volumes (convert to sizes)
  if sP.trMCPerCell eq 0 then begin
    hsml = loadSnapshotSubset(sP=sP,partType='gas',field='hsml')
    hsml = 1.0 * temporary(hsml); increase hsml to decrease visualization noise
  endif else begin
    hsml = loadSnapshotSubset(sP=sP,partType='gas',field='vol')
    hsml = (temporary(hsml) * 3.0 / (4*!pi))^(1.0/3.0) ;cellrad [ckpc]
    hsml = 1.75 * temporary(hsml) ; increase hsml to decrease visualization noise
  endelse
  
  ; cold only
  if keyword_set(coldOnly) then begin
    wCold = where(alog10(temp) le 5.0,count)
    print,' Cold cut have remaining ['+str(count)+'] of ['+str(n_elements(temp))+'] particles.'
    
    temp = temp[wCold]
    mass = mass[wCold]
    hsml = hsml[wCold]
    pos  = pos[*,wCold]
  endif

  ; loop over all non-background subhalos and image
  foreach gcID, gcIDs do begin
  
    ; get subhalo position and size of imaging box
    haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcID]] ;ckpc
    haloMass   = codeMassToLogMsun(gc.subgroupMass[gcID])
    
    boxCen     = sgcen[*,gcID]
    boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcID]])
    boxSizeImg = [boxSize,boxSize,boxSize] ; cube
  
    foreach axisPair, axes do begin
  
      saveFilename = 'map.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+$
                     '.axes'+str(axisPair[0])+str(axisPair[1])+coldFlag+'.sav'
                    
      if ~file_test(sP.derivPath+'sphMaps/'+saveFilename) then begin
      
        ; get box center (in terms of specified axes)
        boxCenImg  = [sgcen[axisPair[0],gcID],sgcen[axisPair[1],gcID],sgcen[3-axisPair[0]-axisPair[1],gcID]]
        
        print,'['+string(gcID,format='(i4)')+'] Mapping ['+str(axisPair[0])+' '+$
              str(axisPair[1])+'] with '+str(boxSize[0])+$
              ' kpc box around subhalo center ['+str(boxCen[0])+' '+str(boxCen[1])+' '+str(boxCen[2])+']'

        ; calculate projection using sph kernel
        sphmap = calcSphMap(pos,hsml,mass,temp,$
                            boxSizeImg=boxSizeImg,boxSizeSim=sP.boxSize,boxCen=boxCenImg,$
                            nPixels=nPixels,axes=axisPair,ndims=3)
      
        ; save
        config = {saveFilename:saveFilename,sizeFac:sizeFac,nPixels:nPixels,axes:axes,boxCen:boxCen,$
                  gcID:gcID,haloMass:haloMass,haloVirRad:haloVirRad,$
                  boxCenImg:boxCenImg,boxSize:boxSize,boxSizeImg:boxSizeImg,sP:sP}
        ;save,sphmap,config,filename=sP.derivPath+'sphMaps/'+saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
      endif else begin
        restore,sP.derivPath +'sphMaps/'+ saveFilename
      endelse
    
      ; plot
      plotSphmapDensQuant, map=sphmap, config=config

    endforeach ;axisPair
  
  endforeach ;gcIDs
  
end

; sphMapHalosDM: run sph kernel density projection on DM particles with boxes centered on halos

pro sphMapHalosDM, sP=sP, gcIDs=gcIDs

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  sizeFac = 3.5      ; times rvir
  nPixels = [800,800] ; px
  
  ;axes = list([0,1],[0,2],[1,2]) ;xy,xz,yz
  axes = list([0,1])
  
  ; target list
  gc    = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sP) 

  ; load dm positions from snapshot and replicate masses from header
  pos  = loadSnapshotSubset(sP=sP,partType='dm',field='pos')
  
  h = loadSnapshotHeader(sP=sP)
  dmPartMass = float(h.massTable[partTypeNum('dm')])
  mass = replicate(dmPartMass,h.nPartTot[1])
  
  ; load HSMLs and veldisps (~temp), assume ID sorted with no gaps (though DM need not be at the end)
  hsmldir = loadHsmlDir(sP=sP,/readHsml,/readVelDisp,/verbose)
  ids     = loadSnapshotSubset(sP=sP,partType='dm',field='ids')
  hsml    = hsmldir.hsml[ids-1]
  veldisp = hsmldir.veldisp[ids-1]
  
  hsmldir = !NULL
  ids = !NULL
  
  ; loop over all non-background subhalos and image
  foreach gcID, gcIDs do begin
  
    ; get subhalo position and size of imaging box
    haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcID]] ;ckpc
    haloMass   = codeMassToLogMsun(gc.subgroupMass[gcID])
    
    boxCen     = sgcen[*,gcID]
    boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcID]])
    boxSizeImg = [boxSize,boxSize,boxSize] ; cube
  
    foreach axisPair, axes do begin
  
      saveFilename = 'dmMap.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+$
                     '.axes'+str(axisPair[0])+str(axisPair[1])+'.sav'
                    
      if ~file_test(sP.derivPath+'sphMaps/'+saveFilename) then begin
      
        ; get box center (in terms of specified axes)
        boxCenImg  = [sgcen[axisPair[0],gcID],sgcen[axisPair[1],gcID],sgcen[3-axisPair[0]-axisPair[1],gcID]]
        
        print,'['+string(gcID,format='(i4)')+'] Mapping ['+str(axisPair[0])+' '+$
              str(axisPair[1])+'] with '+str(boxSize[0])+$
              ' kpc box around subhalo center ['+str(boxCen[0])+' '+str(boxCen[1])+' '+str(boxCen[2])+']'

        ; calculate projection using sph kernel
        dmMap = calcSphMap(pos,hsml,mass,veldisp,$
                           boxSizeImg=boxSizeImg,boxSizeSim=sP.boxSize,boxCen=boxCenImg,$
                           nPixels=nPixels,axes=axisPair,ndims=3)
                            
        ; vertical+horizontal flip the C image arrays so they agree with the (x,y) order in scatterMapHalos
        dmMap.dens_out  = reverse(dmMap.dens_out,2)
        dmMap.dens_out  = reverse(dmMap.dens_out,1)
        dmMap.quant_out = reverse(dmMap.quant_out,2)
        dmMap.quant_out = reverse(dmMap.quant_out,1)
        
        ; save
        config = {saveFilename:saveFilename,sizeFac:sizeFac,nPixels:nPixels,axes:axes,boxCen:boxCen,$
                  gcID:gcID,haloMass:haloMass,haloVirRad:haloVirRad,$
                  boxCenImg:boxCenImg,boxSize:boxSize,boxSizeImg:boxSizeImg,sP:sP}
        ;save,dmMap,config,filename=sP.derivPath+'sphMaps/'+saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
      endif else begin
        restore,sP.derivPath +'sphMaps/'+ saveFilename
      endelse
    
      ; plot
      plotSphmapDensQuant, map=dmMap, config=config

    endforeach ;axisPair
  
  endforeach ;gcIDs
  
end

; makeHaloComparisonImages(): create a mosaic of halo comparison images between gadget and arepo

pro makeHaloComparisonImages, select=select, redshift=redshift

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  res = 512  
  minLogMsun = 10.75
  maxLogMsun = 13.60
  
  sPa = simParams(res=res,run='arepo',redshift=redshift)
  sPg = simParams(res=res,run='gadget',redshift=redshift)
  
  ; gadget group catalog
  if select eq 'gadget' then begin
    gcg = loadGroupCat(sP=sPg,/skipIDs)
    priGIDs = gcIDList(gc=gcg,select='pri')
    priGMasses = codeMassToLogMsun(gcg.subgroupMass[priGIDs])
    w = where(priGMasses ge minLogMsun and priGMasses le maxLogMsun,countG)
    gcIDs_gadget = priGIDs[w]
    
    print,'Mapping ['+str(countG)+'] gadget halos above minLogMsun.'
    sphMapHalos,sP=sPg,gcIDs=gcIDs_gadget;,/coldOnly
    ;scatterMapHalos,sP=sPg,gcIDs=gcIDs_gadget
    ;scatterMapHalosGasDM,sP=sPg,gcIDs=gcIDs_gadget
  endif
  
  ; arepo group catalog
  if select eq 'arepo' then begin
    gca = loadGroupCat(sP=sPa,/skipIDs)
    priAIDs = gcIDList(gc=gca,select='pri')
    priAMasses = codeMassToLogMsun(gca.subgroupMass[priAIDs])
    w = where(priAMasses ge minLogMsun and priAMasses le maxLogMsun,countA)
    gcIDs_arepo = priAIDs[w]
    
    print,'Mapping ['+str(countA)+'] arepo halos above minLogMsun.'
    sphMapHalos,sP=sPa,gcIDs=gcIDs_arepo;,/coldOnly
    ;scatterMapHalos,sP=sPa,gcIDs=gcIDs_arepo
    ;scatterMapHalosGasDM,sP=sPa,gcIDs=gcIDs_arepo
  endif
  
  print,'done.'

end

pro makeHaloComparisonPage
  
  compile_opt idl2, hidden, strictarr, strictarrsubs  
  
  ; config
  res = 512
  redshifts = [0.0,1.0,2.0,3.0]
  
  massBins = [10.75,11.00,11.25,11.50,11.75,12.00,13.60]
  
  axesPairs = ['01','02','12']
  axesNames = ['xy','xz','yz']
  
  redshift = 0.0
  
  haloCounter = 0UL
  
  ; load group catalogs and find halo IDs
  sPa = simParams(res=res,run='arepo',redshift=redshift)
  sPg = simParams(res=res,run='gadget',redshift=redshift)
  
  gcg = loadGroupCat(sP=sPg,/skipIDs)
  gca = loadGroupCat(sP=sPa,/skipIDs)
  
  sgceng = subgroupPosByMostBoundID(sP=sPg)
  sgcena = subgroupPosByMostBoundID(sP=sPa)
  
  priGIDs = gcIDList(gc=gcg,select='pri')
  priGMasses = codeMassToLogMsun(gcg.subgroupMass[priGIDs])   
  priAIDs = gcIDList(gc=gca,select='pri')
  priAMasses = codeMassToLogMsun(gca.subgroupMass[priAIDs])
  
  ; match only for webpage output
  match = findMatchedHalos(sP1=sPa,sP2=sPg)  

  ; loop over mass bins
  for j=0,n_elements(massBins)-2 do begin
    print,massBins[j],massBins[j+1]
    
    ; open file and write header with links
    outputName = "HaloComp_z"+string(redshift,format='(i1)')+"mb"+string(massBins[j],format='(f5.2)')+".htm"
    openw,lun,outputName,/get_lun
    
    printf,lun,"<div class='compHead'><table>"
    printf,lun,"  <tr><th>Redshifts</th><th colspan='"+str(n_elements(massBins)-1)+"'>Log(Msun) Bins</th></tr>"
    foreach red,redshifts,k do begin
      printf,lun,"  <tr><td>"+string(red,format='(f3.1)')+"</td>"
      for mb=0,n_elements(massBins)-2 do begin
        printf,lun,"    <td><a href='HaloComp_z"+string(red,format='(i1)')+"mb"+$
                   string(massBins[mb],format='(f5.2)')+".htm'>"+string(massBins[mb],format='(f5.2)')+$
                   " - "+string(massBins[mb+1],format='(f5.2)')+"</a></td>"
      endfor
      printf,lun,"  </tr>"
    endforeach
    printf,lun,"</table></div>"
    
    printf,lun,""
    printf,lun,"<div class='compBody'>"
    
    ; redshift header
    printf,lun,""
    printf,lun,"<h2>z = "+string(redshift,format='(f3.1)')+"</h2>"
    printf,lun,""
  
    ; halo IDs in this mass bin
    wG = where(priGMasses ge massBins[j] and priGMasses lt massBins[j+1],countG)
    if countG gt 0 then gcIDs_gadget = priGIDs[wG]
    wA = where(priAMasses ge massBins[j] and priAMasses lt massBins[j+1],countA)
    if countA gt 0 then gcIDs_arepo = priAIDs[wA]
    
    ; decide which arepo halos have matching gadget halos
    if countA gt 0 then begin
      matchedInds = match.matchedInds[wA]
      wMatch = where(matchedInds ne -1,countMatch,comp=wNoMatch,ncomp=countNoMatchA)
      
      ; find subset of gadget halos that were not matched
      match,match.matchedInds,wG,inds_match,inds_wG,count=countNoMatchG,/sort
      if countNoMatchG gt 0 then begin
        wGMatched = wG[inds_wG]
        wGNoMatch = removeIntersectionFromB(wGMatched,wG)
        countNoMatchG = n_elements(wGNoMatch)
      endif
    endif else begin
      ; if arepo has no halos in this massbin, there can be no matches
      countMatch = 0
      countNoMatchA = countA
      countNoMatchG = countG
      wGNoMatch = wG
    endelse

    ; mass bin header
    printf,lun,""
    printf,lun,"<a name='z"+str(k)+"_massbin"+str(j)+"'></a>"
    printf,lun,"<h3>"+string(massBins[j],format='(f4.1)')+" < log(M) < "+$
                      string(massBins[j+1],format='(f4.1)')+"</h3>"
    printf,lun,""    
    printf,lun,"<table>"
    printf,lun," <tr><th>Info</th><th>Gadget</th><th>Arepo</th></tr>"
    
    ; matched
    ; -------
    for i=0,countMatch-1 do begin
      haloCounter += 1
      
      ; make stats string
      gcIndA = priAIDs[wA[wMatch[i]]]
      gcIndG = priGIDs[matchedInds[wMatch[i]]]
      
      haloMassA = codeMassToLogMsun(gca.subgroupMass[gcIndA])
      haloMassG = codeMassToLogMsun(gcg.subgroupMass[gcIndG])
      
      virRadA  = gca.group_r_crit200[gca.subgroupGrNr[gcIndA]]
      virRadG  = gcg.group_r_crit200[gcg.subgroupGrNr[gcIndG]]
      virTempA = alog10(codeMassToVirTemp(gca.subgroupMass[gcIndA],sP=sPa))
      virTempG = alog10(codeMassToVirTemp(gcg.subgroupMass[gcIndG],sP=sPg))
      
      statsString = "<span class='haloid'>ID: z"+string(redshift,format='(i1)')+"."+$
                    str(haloCounter)+"</span><br><br> "+$
                    "Gadget:<br> log(M) = "+string(haloMassG,format='(f5.2)')+" <br> "+$
                    " xyz = "+string(sgceng[0,gcIndG],format='(f7.1)')+" "+$
                              string(sgceng[1,gcIndG],format='(f7.1)')+" "+$
                              string(sgceng[2,gcIndG],format='(f7.1)')+" <br> "+$
                    " virRad = "+string(virRadG,format='(i3)')+" ckpc<br> "+$
                    " virTemp = "+string(virTempG,format='(f3.1)')+" <br><br> "+$
                    "Arepo:<br> log(M) = "+string(haloMassA,format='(f5.2)')+" <br> "+$
                    " xyz = "+string(sgcena[0,gcIndA],format='(f7.1)')+" "+$
                              string(sgcena[1,gcIndA],format='(f7.1)')+" "+$
                              string(sgcena[2,gcIndA],format='(f7.1)')+" <br> "+$
                    " virRad = "+string(virRadA,format='(i3)')+" ckpc<br> "+$
                    " virTemp = "+string(virTempA,format='(f3.1)')+" <br> "
                              
      ; make gadget image thumbnail/link string
      gadgetString = ''
      
      foreach axisPair,axesPairs,m do begin        
        ; sphMap
        fpath = "sphMaps.gadget.z"+string(redshift,format='(f3.1)')+"/"
        fname = "map.G"+str(res)+"."+str(sPg.snap)+".h"+str(gcIndG)+".axes"+axisPair+".png"
        gadgetString += axesNames[m]+": <a title='Gadget' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a>"
                        
        ; scatterMap
        fpath = "scatterMaps.gadget.z"+string(redshift,format='(f3.1)')+"/"
        fname = "scatter.G"+str(res)+"."+str(sPg.snap)+".h"+str(gcIndG)+".axes"+axisPair+".png"
        gadgetString += " <a title='Gadget' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a>"
        
        ; gasdm
        fpath = "gasdm.gadget.z"+string(redshift,format='(f3.1)')+"/"
        fname = "gasdm.G"+str(res)+"."+str(sPg.snap)+".h"+str(gcIndG)+".axes"+axisPair+".png"
        gadgetString += " <a title='Gadget' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a><br>"
      endforeach
      
      gadgetString = strmid(gadgetString,0,strlen(gadgetString)-4)
      
      ; make arepo image thumbnail/link string
      arepoString = ''
      
      foreach axisPair,axesPairs,m do begin        
        ; sphMap
        fpath = "sphMaps.arepo.z"+string(redshift,format='(f3.1)')+"/"
        fname = "map.A"+str(res)+"."+str(sPa.snap)+".h"+str(gcIndA)+".axes"+axisPair+".png"
        arepoString += axesNames[m]+": <a title='Arepo' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a>"
                        
        ; scatterMap
        fpath = "scatterMaps.arepo.z"+string(redshift,format='(f3.1)')+"/"
        fname = "scatter.A"+str(res)+"."+str(sPa.snap)+".h"+str(gcIndA)+".axes"+axisPair+".png"
        arepoString += " <a title='Arepo' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a>"
       
        ; gasdm
        fpath = "gasdm.arepo.z"+string(redshift,format='(f3.1)')+"/"
        fname = "gasdm.A"+str(res)+"."+str(sPa.snap)+".h"+str(gcIndA)+".axes"+axisPair+".png"
        arepoString += " <a title='Arepo' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a><br>"
      endforeach
      
      arepoString = strmid(arepoString,0,strlen(arepoString)-4)
      
      printf,lun," <tr class='row'>"
      printf,lun,"  <td class='stcell'>"+statsString+"</td>"
      printf,lun,"  <td class='gacell'>"+gadgetString+"</td>"
      printf,lun,"  <td class='arcell'>"+arepoString+"</td>"
      printf,lun," </tr>"
    endfor
    
    ; arepo (unmatched)
    ; -----------------
    for i=0,countNoMatchA-1 do begin
      haloCounter += 1
      
      ; make stats string
      gcIndA = priAIDs[wA[wNoMatch[i]]]
      
      haloMassA = codeMassToLogMsun(gca.subgroupMass[gcIndA])
      virRadA  = gca.group_r_crit200[gca.subgroupGrNr[gcIndA]]
      virTempA = alog10(codeMassToVirTemp(gca.subgroupMass[gcIndA],sP=sPa))
      
      statsString = "<span class='haloid'>ID: z"+string(redshift,format='(i1)')+"."+$
                    str(haloCounter)+"</span><br><br> "+$
                    "Gadget:<br> <i>Unmatched.</i> <br><br> "+$
                    "Arepo:<br> log(M) = "+string(haloMassA,format='(f5.2)')+" <br> "+$
                    " xyz = "+string(sgcena[0,gcIndA],format='(f7.1)')+" "+$
                              string(sgcena[1,gcIndA],format='(f7.1)')+" "+$
                              string(sgcena[2,gcIndA],format='(f7.1)')+" <br> "+$
                    " virRad = "+string(virRadA,format='(i3)')+" ckpc<br> "+$
                    " virTemp = "+string(virTempA,format='(f3.1)')+" <br> "
                       
      ; make arepo image thumbnail/link string
      arepoString = ''
      
      foreach axisPair,axesPairs,m do begin        
        ; sphMap
        fpath = "sphMaps.arepo.z"+string(redshift,format='(f3.1)')+"/"
        fname = "map.A"+str(res)+"."+str(sPa.snap)+".h"+str(gcIndA)+".axes"+axisPair+".png"
        arepoString += axesNames[m]+": <a title='Arepo' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a>"
                        
        ; scatterMap
        fpath = "scatterMaps.arepo.z"+string(redshift,format='(f3.1)')+"/"
        fname = "scatter.A"+str(res)+"."+str(sPa.snap)+".h"+str(gcIndA)+".axes"+axisPair+".png"
        arepoString += " <a title='Arepo' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a>"
        
        ; gasdm
        fpath = "gasdm.arepo.z"+string(redshift,format='(f3.1)')+"/"
        fname = "gasdm.A"+str(res)+"."+str(sPa.snap)+".h"+str(gcIndA)+".axes"+axisPair+".png"
        arepoString += " <a title='Arepo' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a><br>"
      endforeach
      
      arepoString = strmid(arepoString,0,strlen(arepoString)-4)
      
      printf,lun," <tr class='row'>"
      printf,lun,"  <td class='stcell'>"+statsString+"</td>"
      printf,lun,"  <td class='gacell'><i>(Unmatched)</i></td>"
      printf,lun,"  <td class='arcell'>"+arepoString+"</td>"
      printf,lun," </tr>"
    endfor
    
    ; gadget (unmatched)
    ; -----------------
    for i=0,countNoMatchG-1 do begin
      ; check: was this gadget halo matched by an arepo halo in some other mass bin? if so just skip
      ;gcIndG = priGIDs[wGNoMatch[i]]
      ;w = where(match.matchedInds eq gcIndG,countOMM)
      ;if countOMM gt 0 then continue
      
      haloCounter += 1
      
      ; make stats string
      haloMassG = codeMassToLogMsun(gcg.subgroupMass[gcIndG])
      virRadG   = gcg.group_r_crit200[gcg.subgroupGrNr[gcIndG]]
      virTempG  = alog10(codeMassToVirTemp(gcg.subgroupMass[gcIndG],sP=sPg))
      
      statsString = "<span class='haloid'>ID: z"+string(redshift,format='(i1)')+"."+$
                    str(haloCounter)+"</span><br><br> "+$
                    "Gadget:<br> log(M) = "+string(haloMassG,format='(f5.2)')+" <br> "+$
                    " xyz = "+string(sgceng[0,gcIndG],format='(f7.1)')+" "+$
                              string(sgceng[1,gcIndG],format='(f7.1)')+" "+$
                              string(sgceng[2,gcIndG],format='(f7.1)')+" <br> "+$
                    " virRad = "+string(virRadG,format='(i3)')+" ckpc<br> "+$
                    " virTemp = "+string(virTempG,format='(f3.1)')+" <br><br> "+$
                    "Arepo:<br> <i>Unmatched.</i>  <br> "
                       
      ; make arepo image thumbnail/link string
      gadgetString = ''
      
      foreach axisPair,axesPairs,m do begin        
        ; sphMap
        fpath = "sphMaps.gadget.z"+string(redshift,format='(f3.1)')+"/"
        fname = "map.G"+str(res)+"."+str(sPg.snap)+".h"+str(gcIndG)+".axes"+axisPair+".png"
        gadgetString += axesNames[m]+": <a title='Gadget' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a>"
                        
        ; scatterMap
        fpath = "scatterMaps.gadget.z"+string(redshift,format='(f3.1)')+"/"
        fname = "scatter.G"+str(res)+"."+str(sPg.snap)+".h"+str(gcIndG)+".axes"+axisPair+".png"
        gadgetString += " <a title='Gadget' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a>"
        
        ; gasdm
        fpath = "gasdm.gadget.z"+string(redshift,format='(f3.1)')+"/"
        fname = "gasdm.G"+str(res)+"."+str(sPg.snap)+".h"+str(gcIndG)+".axes"+axisPair+".png"
        gadgetString += " <a title='Gadget' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a><br>"
      endforeach
      
      gadgetString = strmid(gadgetString,0,strlen(gadgetString)-4)
      
      printf,lun," <tr class='row'>"
      printf,lun,"  <td class='stcell'>"+statsString+"</td>"
      printf,lun,"  <td class='gacell'>"+gadgetString+"</td>"
      printf,lun,"  <td class='arcell'><i>(Unmatched)</i></td>"
      printf,lun," </tr>"
    endfor
    
    printf,lun,"</table>"
    printf,lun,""
    printf,lun,"</div>"
    printf,lun,""
    printf,lun,"<p id='footer'>The end.</p>"
    printf,lun,"</body>"
    printf,lun,"</html>"
     
    ; close file
    free_lun,lun
    
    ; if header.txt exists, cat the two together
    if file_test("header.txt") then begin
      print,' added header'
      spawn,'cat header.txt '+outputName+' > newpage.htm' ;make new
      spawn,'mv newpage.htm '+outputName ;overwrite original
    endif
        
  endfor ;massBins
  
end
