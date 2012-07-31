; cosmoVis.pro
; cosmological boxes - 2d visualization
; dnelson jul.2012

; webglCutouts(): make halo centered cutouts for the webGL app

pro webglCutouts;, sP=sP, gcIDs=gcIDs

  compile_opt idl2, hidden, strictarr, strictarrsubs

  sP = simParams(res=512,run='gadget',redshift=2.0)
  gcIDs = [2342] ;g2342 a2132 (z2.304)
  ;if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  sizeFac = 5.0 ; times rvir for the bounding box of each cutout
  
  ; target list
  gc    = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sP) 

  ; load gas properties
  u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
  nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
  pos   = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  dens  = loadSnapshotSubset(sP=sP,partType='gas',field='dens')
  hsml  = loadSnapshotSubset(sP=sP,partType='gas',field='hsml')
  vel   = loadSnapshotSubset(sP=sP,partType='gas',field='vel')

  ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
  print,'shuffling...'
  iseed = 424242L
  sort_inds = sort(randomu(iseed,n_elements(u)))
  
  u     = u[sort_inds]
  nelec = nelec[sort_inds]
  dens  = dens[sort_inds]
  hsml  = hsml[sort_inds]
  pos   = pos[*,sort_inds]
  vel   = vel[*,sort_inds]
  
  sort_inds = !NULL
  
  if sP.trMCPerCell eq 0 then simType = 1 ;gadget
  if sP.trMCPerCell ne 0 then simType = 2 ;arepo
  
  print,'saving...'
  ; loop over all requested halos
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
    
    rvir = gc.group_r_crit200[gc.subgroupGrNr[gcID]]
    virtemp = codeMassToVirTemp(gc.subgroupmass[gcID],sP=sP)
    rad = reform(sqrt(xDist*xDist + yDist*yDist + zDist*zDist)) / rvir[0]
  
    ; local (cube) cutout
    wCut = where(abs(xDist) le 0.5*boxSize and abs(yDist) le 0.5*boxSize and $
                 abs(zDist) le 0.5*boxSize,nCutout)
    print,gcID,nCutout
           
    ; field subsets   
    loc_u     = u[wCut]
    loc_nelec = nelec[wCut]
    loc_dens  = alog10(dens[wCut] * 1e10) ; log msun/ckpc^3 (i.e. comoving density)
    loc_hsml  = hsml[wCut]
    
    ; derived
    loc_temp = alog10(convertUtoTemp(loc_u,loc_nelec)) ; log T (K)
    loc_pres = alog10(calcPressureCGS(loc_u,loc_dens,sP=sP)) ; log P/k (K/cm^3)
    loc_entr = alog10(calcEntropyCGS(loc_u,loc_dens,sP=sP)) ; log S (K cm^2)
    
    ; relative positions
    loc_pos  = fltarr(3,nCutout)
    loc_pos[0,*] = xDist[wCut] ; delta
    loc_pos[1,*] = yDist[wCut]
    loc_pos[2,*] = zDist[wCut]
    
    xDist = !NULL
    yDist = !NULL
    zDist = !NULL
    
    ; velocities: calculate radial velocity relative to bulk halo motion
    loc_vel = vel[*,wCut]
    
    gVel = gc.subgroupVel[*,gcID]
    loc_vel[0,*] = reform(loc_vel[0,*] - gVel[0])
    loc_vel[1,*] = reform(loc_vel[1,*] - gVel[1])
    loc_vel[2,*] = reform(loc_vel[2,*] - gVel[2])
    
    ; make normalized position vector wrt halo center = vec(r) / ||r|| where r from particle to center
    ; means that radvel<0 is inflow and radvel>0 is outflow
    rnorm0 = reform(loc_pos[0,*] - boxCen[0])
    rnorm1 = reform(loc_pos[1,*] - boxCen[1])
    rnorm2 = reform(loc_pos[2,*] - boxCen[2])
    
    correctPeriodicDistVecs, rnorm0, sP=sP
    correctPeriodicDistVecs, rnorm1, sP=sP
    correctPeriodicDistVecs, rnorm2, sP=sP
    
    ; denominator and do divide
    rnorm = sqrt(rnorm0*rnorm0 + rnorm1*rnorm1 + rnorm2*rnorm2)

    rnorm0 /= rnorm
    rnorm1 /= rnorm
    rnorm2 /= rnorm
    
    ; dot(vel,rnorm) gives the magnitude of the projection of vel onto vec(r)
    loc_radvel = loc_vel[0,*]*rnorm0 + loc_vel[1,*]*rnorm1 + loc_vel[2,*]*rnorm2 ; 1xN

    ; velocities: create endpoint for each position point for the velocity vector line
    ;loc_pos2 = fltarr(3,nCutout)
    ;loc_pos2[0,*] = loc_pos[0,*] + loc_vel[0,*]*velVecFac
    ;loc_pos2[1,*] = loc_pos[1,*] + loc_vel[1,*]*velVecFac
    ;loc_pos2[2,*] = loc_pos[2,*] + loc_vel[2,*]*velVecFac
  
    loc_hsml = reform(loc_hsml,[1,nCutout]) ; make 1xN vector
    loc_temp = reform(loc_temp,[1,nCutout])
    loc_dens = reform(loc_dens,[1,nCutout])
    loc_pres = reform(loc_pres,[1,nCutout])
    loc_entr = reform(loc_entr,[1,nCutout])

    ; prepare output data
    dataout = [loc_pos,loc_hsml,loc_temp,loc_dens,loc_pres,loc_entr,loc_radvel]
    
    ; write binary file
    fileName = 'cutout.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+'.dat'
    
    outStruct = { fileType : fix(1)              ,$ ; 0=gas[x,y,z,hsml,temp,dens]
                                                  $ ; 1=gas[x,y,z,hsml,temp,dens,pres,entr,radvel]
                  simType  : fix(simType)        ,$ ; 1=gadget, 2=arepo
                  nPts     : long(nCutout)       ,$
                  redshift : float(sP.redshift)  ,$
                  hInd     : long(gcID)          ,$ ; subgroup id
                  hVirRad  : float(rvir)         ,$ ; ckpc
                  hVirTemp : float(virtemp)      ,$ ; log K
                  sizeFac  : float(sizeFac)      ,$ ; bounding box size
                  data     : reform(dataout,n_elements(dataout)) $
                }
    
    ; write file
    openw,lun,fileName,/get_lun
    writeu,lun,outStruct
    close,lun
  endforeach
end

; sphMapBox: run sph kernel density projection on whole box

pro sphMapBox, partType=partType

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  res = 128
  run = 'tracernew'
  ;partType = 'tracervel'

  redshift = 0.0
  
  nPixels = [800,800] ;px

  zoomFac = 1    ; only in axes, not along projection direction
  nNGB    = 64   ; use CalcHSML for HSML with nNGB
  axes    = [0,1] ; x,y

  ; paths and render config
  sP = simParams(res=res,run=run,redshift=redshift)
  h = loadSnapshotHeader(sP=sP)
  
  boxSizeImg = [h.boxSize,h.boxSize,h.boxSize]              ;kpc
  boxCen     = [h.boxSize/2.0,h.boxSize/2.0,h.boxSize/2.0]  ;kpc
  
  foreach k,axes do boxSizeImg[k] /= zoomFac
  
  outFilename = 'sphmap.box_'+str(zoomFac)+'.nNGB='+str(nNGB)+'.snap='+str(sP.snap)+$
                '.box.axis0='+str(axes[0])+'.axis1='+str(axes[0])+'.'+partType

  ; save/restore
  if (file_test(sP.derivPath + outFilename + '.sav')) then begin
    restore,sP.derivPath + outFilename + '.sav',/verbose
  endif else begin

    ; gas cells or sph particles (pos,hsml stored in snapshots)
    if partType eq 'gas' then begin
      mass = loadSnapshotSubset(sP=sP,partType=partType,field='mass',/verbose)
      pos  = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
      hsml = loadSnapshotSubset(sP=sP,partType=partType,field='hsml')
      quant = replicate(1.0,h.nPartTot[0]) ; dummy quant for now
    endif
    
    ; velocity tracers (pos stored in snapshots, calculate hsmls and use constant eff mass)
    if partType eq 'tracervel' then begin
      mass = replicate(sP.targetGasMass, h.nPartTot[2])
      pos  = loadSnapshotSubset(sP=sP,partType=partType,field='pos',/verbose)
      hsml = calcHSML(pos,ndims=3,nNGB=nNGB,boxSize=h.boxSize)
      quant = replicate(1.0,h.nPartTot[2]) ; dummy quant for now
    endif
    
    ; monte carlo tracers (use gas pos,hsml and replace mass by num_child_tr*tr_mass_eff)
    if partType eq 'tracermc' then begin
      mass = loadSnapshotSubset(sP=sP,partType='gas',field='numtr',/verbose)
      mass *= sP.trMassConst
      pos  = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
      hsml = loadSnapshotSubset(sP=sP,partType='gas',field='hsml')
      quant = replicate(1.0,h.nPartTot[0]) ; dummy quant for now
    endif
    
    sphMap = calcSphMap(pos,hsml,mass,quant,boxSizeImg=boxSizeImg,boxSizeSim=h.boxSize,boxCen=boxCen,$
                        nPixels=nPixels,axes=axes,ndims=3)
    colMassMap = sphMap.dens_out
    
    save,colMassMap,filename=sP.derivPath + outFilename + '.sav'
  endelse

  ; rescale
  maxVal = max(colMassMap)/2.0
  ;maxVal = 0.5
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
  xMinMax = [boxCen[0]-boxSizeImg[0]/2.0,boxCen[0]+boxSizeImg[0]/2.0]
  yMinMax = [boxCen[1]-boxSizeImg[1]/2.0,boxCen[1]+boxSizeImg[1]/2.0]
  
  start_PS, sP.plotPath + outFilename + '.eps'
    loadct, 4, bottom=1, /silent
    tvim,colMassMap,xrange=xMinMax,yrange=yMinMax,pcharsize=0.0001
  end_PS, pngResize=68, /deletePS
  
end

; plotScatterComp(): plot side by side colored/vectorized scatter plots

pro plotScatterComp, pos_left, pos2_left, pos_right, pos2_right, cinds_left, cinds_right, $
                     config=config, top=top, bottom=bottom

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
      
        ; color table and establish temperature -> color mapping
        cgColorfill,[1,1,0,0,1],[1,0,0,1,1],/normal,color=cgColor('black')
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
          cgText,0.365,0.036*fac,'4',alignment=0.5,color=cgColor('black'),/normal
          cgText,0.635,0.036*fac,'7',alignment=0.5,color=cgColor('black'),/normal
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
        if keyword_set(top) then cgText,0.5,0.96,"GADGET",charsize=!p.charsize+0.5,$
          alignment=0.5,/normal,color=cgColor('white')
        if keyword_set(bottom) then cgText,0.5,0.46,"AREPO",charsize=!p.charsize+0.5,$
          alignment=0.5,/normal,color=cgColor('white')
        
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
  sizeFac = 3.5       ; times rvir
  cutFac  = 1.0       ; times boxSize
  
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
    
    rvir = gc.group_r_crit200[gc.subgroupGrNr[gcID]]
    rad = reform(sqrt(xDist*xDist + yDist*yDist + zDist*zDist)) / rvir[0]
  
    ; local (cube) cutout
    wCut = where(abs(xDist) le 0.5*cutFac*boxSize and abs(yDist) le 0.5*cutFac*boxSize and $
                 abs(zDist) le 0.5*cutFac*boxSize,nCutout)
                 
    ; local (cube) cutout - only blobs or only non-blobs
    ;wCut = where(abs(xDist) le 0.5*cutFac*boxSize and abs(yDist) le 0.5*cutFac*boxSize and $
    ;             abs(zDist) le 0.5*cutFac*boxSize and alog10(dens*1e10) gt 5.0 and rad gt 0.125,nCutout)
    ;wCut = where(abs(xDist) le 0.5*cutFac*boxSize and abs(yDist) le 0.5*cutFac*boxSize and $
    ;             abs(zDist) le 0.5*cutFac*boxSize and (alog10(dens*1e10) le 5.0 or rad le 0.125),nCutout)
    
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
                axisPair:axisPair,sP:sP,barMM:tempMinMax,barType:'none'} ;'1temp'
      
      ; plot
      plotScatterComp,loc_pos,loc_pos2,loc_pos_cold,loc_pos2_cold,colorinds,colorinds_cold,config=config            

    endforeach ;axisPair

  endforeach ;gcIDs

end

; scatterMapHalosComp(): same as above but compare arepo/gadget top/bottom

pro scatterMapHalosComp

  compile_opt idl2, hidden, strictarr, strictarrsubs

  sPg = simParams(res=512,run='gadgetold',redshift=2.0)
  sPa = simParams(res=512,run='arepo',redshift=2.0)
  
  ;z2.304 g2342 a2132 -- z2.301 g2289 a2034 -- z2.130 g6369 a5966 axes02 -- z2.64 g5498 a5097
  gcIDg = 6369
  gcIDa = 5966

  ; config
  sizeFac     = 3.5       ; times rvir
  tempMinMax  = [4.0,7.0] ; log(K)
  coldTempCut = 5.0       ; log(K)
  velVecFac   = 0.01      ; times velocity (km/s) in plotted kpc
  axisPair    = [0,2]     ; xy
  
  ; GADGET
  ; ------
  gc    = loadGroupCat(sP=sPg,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sPg) 
  
  saveFilename = sPg.derivPath + 'cutout.' + sPg.savPrefix + str(sPg.res) + '.' + str(sPg.snap) + $
                 '.h' + str(gcIDg) + '.sf' + str(fix(sizeFac*10)) + '.sav'
  
  if file_test(saveFilename) then begin
    restore,saveFilename
  endif else begin
    ; load u,nelec and calculate temperature
    u     = loadSnapshotSubset(sP=sPg,partType='gas',field='u')
    nelec = loadSnapshotSubset(sP=sPg,partType='gas',field='nelec')
    temp  = alog10(convertUtoTemp(u,nelec))
    u     = !NULL
    nelec = !NULL
  
    ; load gas positions and velocities
    pos = loadSnapshotSubset(sP=sPg,partType='gas',field='pos')
    vel = loadSnapshotSubset(sP=sPg,partType='gas',field='vel')
    
    ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
    print,'shuffling...'
    iseed = 424242L
    sort_inds = sort(randomu(iseed,n_elements(temp)))
    
    temp = temp[sort_inds]
    pos  = pos[*,sort_inds]
    vel  = vel[*,sort_inds]
    
    sort_inds = !NULL
    
    print,'cutout...'
    ; get subhalo position and size of imaging box
    boxCen     = sgcen[*,gcIDg]
    boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcIDg]] / 10.0) * 10.0
    boxSizeImg = [boxSize,boxSize,boxSize] ; cube
  
    ; make conservative cutout greater than boxsize accounting for periodic (do cube not sphere)
    xDist = pos[0,*] - boxCen[0]
    yDist = pos[1,*] - boxCen[1]
    zDist = pos[2,*] - boxCen[2]
    
    correctPeriodicDistVecs, xDist, sP=sPg
    correctPeriodicDistVecs, yDist, sP=sPg
    correctPeriodicDistVecs, zDist, sP=sPg
    
    rvir = gc.group_r_crit200[gc.subgroupGrNr[gcIDg]]
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
    
    ; save
    save,loc_pos,loc_temp,loc_pos2,sPg,gcIDg,sizeFac,boxSizeImg,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sPg.derivPath))
  
  endelse

  ; create color index mapping
  colorinds = (loc_temp-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
  colorinds = fix(colorinds + 50.0) > 0 < 255 ;50-255  

  ; local (cold) cutout
  wCold = where(loc_temp le coldTempCut,nCutoutCold)
  loc_temp = !NULL
  
  loc_pos_cold   = loc_pos[*,wCold]
  loc_pos2_cold  = loc_pos2[*,wCold]
  colorinds_cold = colorinds[wCold]

  print,'rendering gadget...'
  ; get box center (in terms of specified axes)
  boxCenImg  = [sgcen[axisPair[0],gcIDg],sgcen[axisPair[1],gcIDg],sgcen[3-axisPair[0]-axisPair[1],gcIDg]]
  haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcIDg]] ;ckpc
  haloMass = codeMassToLogMsun(gc.subgroupMass[gcIDg])
  
  ; fix halo mass if we're using old (x2 bug) catalogs
  if sPg.run eq 'gadgetold' then haloMass = codeMassToLogMsun(0.5*gc.subgroupMass[gcIDg])
  
  plotFilename = 'scatter.'+sPg.savPrefix+str(sPg.res)+'.h'+str(gcIDg)+'-'+sPa.savPrefix+str(sPa.res)+$
                 '.h'+str(gcIDa)+'.'+str(sPg.snap)+'.axes'+str(axisPair[0])+str(axisPair[1])+'.eps'

  config = {boxSizeImg:boxSizeImg,plotFilename:plotFilename,haloVirRad:haloVirRad,haloMass:haloMass,$
            axisPair:axisPair,sP:sPg,barMM:tempMinMax,barType:'1temp'}
  
  ; plot
  start_PS, sPg.plotPath + config.plotFilename, xs=8, ys=8
  plotScatterComp,loc_pos,loc_pos2,loc_pos_cold,loc_pos2_cold,colorinds,colorinds_cold,config=config,/top           

  ; AREPO
  ; -----
  gc    = loadGroupCat(sP=sPa,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sPa) 
  
  saveFilename = sPa.derivPath + 'cutout.' + sPa.savPrefix + str(sPa.res) + '.' + str(sPa.snap) + $
                 '.h' + str(gcIDa) + '.sf' + str(fix(sizeFac*10)) + '.sav'
  
  if file_test(saveFilename) then begin
    restore,saveFilename
  endif else begin
    ; load u,nelec and calculate temperature
    u     = loadSnapshotSubset(sP=sPa,partType='gas',field='u')
    nelec = loadSnapshotSubset(sP=sPa,partType='gas',field='nelec')
    temp  = alog10(convertUtoTemp(u,nelec))
    u     = !NULL
    nelec = !NULL
  
    ; load gas positions and velocities
    pos = loadSnapshotSubset(sP=sPa,partType='gas',field='pos')
    vel = loadSnapshotSubset(sP=sPa,partType='gas',field='vel')
    
    ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
    print,'shuffling...'
    iseed = 424242L
    sort_inds = sort(randomu(iseed,n_elements(temp)))
    
    temp = temp[sort_inds]
    pos  = pos[*,sort_inds]
    vel  = vel[*,sort_inds]
    
    sort_inds = !NULL
    
    print,'cutout...'
    ; get subhalo position and size of imaging box
    boxCen     = sgcen[*,gcIDa]
    boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcIDa]] / 10.0) * 10.0
    boxSizeImg = [boxSize,boxSize,boxSize] ; cube
  
    ; make conservative cutout greater than boxsize accounting for periodic (do cube not sphere)
    xDist = pos[0,*] - boxCen[0]
    yDist = pos[1,*] - boxCen[1]
    zDist = pos[2,*] - boxCen[2]
    
    correctPeriodicDistVecs, xDist, sP=sPa
    correctPeriodicDistVecs, yDist, sP=sPa
    correctPeriodicDistVecs, zDist, sP=sPa
    
    rvir = gc.group_r_crit200[gc.subgroupGrNr[gcIDa]]
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
    
    ; save
    save,loc_pos,loc_temp,loc_pos2,sPa,gcIDa,boxSizeImg,sizeFac,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sPa.derivPath))
  
  endelse
  
  ; create color index mapping
  colorinds = (loc_temp-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
  colorinds = fix(colorinds + 50.0) > 0 < 255 ;50-255  

  ; local (cold) cutout
  wCold = where(loc_temp le coldTempCut,nCutoutCold)
  loc_temp = !NULL
  
  loc_pos_cold   = loc_pos[*,wCold]
  loc_pos2_cold  = loc_pos2[*,wCold]
  colorinds_cold = colorinds[wCold]

  print,'rendering arepo...'
  ; get box center (in terms of specified axes)
  boxCenImg  = [sgcen[axisPair[0],gcIDa],sgcen[axisPair[1],gcIDa],sgcen[3-axisPair[0]-axisPair[1],gcIDa]]
  haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcIDa]] ;ckpc
  haloMass = codeMassToLogMsun(gc.subgroupMass[gcIDa])

  ; fix halo mass if we're using old (x2 bug) catalogs
  if sPa.run eq 'arepo' then haloMass = codeMassToLogMsun(0.5*gc.subgroupMass[gcIDa])

  config = {boxSizeImg:boxSizeImg,plotFilename:plotFilename,haloVirRad:haloVirRad,haloMass:haloMass,$
            axisPair:axisPair,sP:sPa,barMM:tempMinMax,barType:'1temp'}
  
  ; plot
  plotScatterComp,loc_pos,loc_pos2,loc_pos_cold,loc_pos2_cold,colorinds,colorinds_cold,config=config,/bottom
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
  velVecFac   = 0.01      ; times velocity (km/s) in plotted kpc
  
  ; GADGET
  ; ------
  gc    = loadGroupCat(sP=sPg,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sPg) 
  
  saveFilenames = sPg.derivPath + 'cutout.' + sPg.savPrefix + str(sPg.res) + '.' + str(sPg.snap) + $
                  '.h' + str(gcIDsG) + '.sf' + str(fix(sizeFac*10)) + '.sav'
  readFlag = 0
  foreach saveFilename,saveFilenames do if ~file_test(saveFilename) then readFlag = 1
  
  if readFlag then begin
    ; load u,nelec and calculate temperature
    u     = loadSnapshotSubset(sP=sPg,partType='gas',field='u')
    nelec = loadSnapshotSubset(sP=sPg,partType='gas',field='nelec')
    temp  = alog10(convertUtoTemp(u,nelec))
    u     = !NULL
    nelec = !NULL
  
    ; load gas positions and velocities
    pos = loadSnapshotSubset(sP=sPg,partType='gas',field='pos')
    vel = loadSnapshotSubset(sP=sPg,partType='gas',field='vel')
    
    ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
    print,'shuffling...'
    iseed = 424242L
    sort_inds = sort(randomu(iseed,n_elements(temp)))
    
    temp = temp[sort_inds]
    pos  = pos[*,sort_inds]
    vel  = vel[*,sort_inds]
    
    sort_inds = !NULL
    
    print,'cutout...'
    foreach gcIDg,gcIDsG,k do begin
      ; get subhalo position and size of imaging box
      boxCen     = sgcen[*,gcIDg]
      boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcIDg]] / 10.0) * 10.0
      boxSizeImg = [boxSize,boxSize,boxSize] ; cube
    
      ; make conservative cutout greater than boxsize accounting for periodic (do cube not sphere)
      xDist = pos[0,*] - boxCen[0]
      yDist = pos[1,*] - boxCen[1]
      zDist = pos[2,*] - boxCen[2]
      
      correctPeriodicDistVecs, xDist, sP=sPg
      correctPeriodicDistVecs, yDist, sP=sPg
      correctPeriodicDistVecs, zDist, sP=sPg
      
      rvir = gc.group_r_crit200[gc.subgroupGrNr[gcIDg]]
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
      
      ; save
      save,loc_pos,loc_temp,loc_pos2,sPg,gcIDg,sizeFac,boxSizeImg,filename=saveFilenames[k]
      print,'Saved: '+strmid(saveFilenames[k],strlen(sPg.derivPath))
    endforeach
  
  endif ;readFlag

  print,'rendering gadget...'
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
    restore,saveFilenames[k]
    
    ; create color index mapping
    colorinds = (loc_temp-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
    colorinds = fix(colorinds + 50.0) > 0 < 255 ;50-255  
  
    ; local (cold) cutout
    wCold = where(loc_temp le coldTempCut,nCutoutCold)
    loc_temp = !NULL
    
    loc_pos_cold   = loc_pos[*,wCold]
    loc_pos2_cold  = loc_pos2[*,wCold]
    colorinds_cold = colorinds[wCold]
  
    ; get box center (in terms of specified axes)
    axisPair   = axes[k]
    boxCenImg  = [sgcen[axisPair[0],gcIDg],sgcen[axisPair[1],gcIDg],sgcen[3-axisPair[0]-axisPair[1],gcIDg]]
    haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcIDg]] ;ckpc
    haloMass = codeMassToLogMsun(gc.subgroupMass[gcIDg])
    
    ; fix halo mass if we're using old (x2 bug) catalogs
    if sPg.run eq 'gadgetold' then haloMass = codeMassToLogMsun(0.5*gc.subgroupMass[gcIDg])
  
    gaHaloMasses = [gaHaloMasses,haloMass]
  
    config = {boxSizeImg:boxSizeImg,plotFilename:plotFilename,haloVirRad:haloVirRad,haloMass:haloMass,$
              axisPair:axisPair,sP:sPg,barMM:tempMinMax,barType:'1temp'}
    
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
  
  saveFilenames = sPa.derivPath + 'cutout.' + sPa.savPrefix + str(sPa.res) + '.' + str(sPa.snap) + $
                  '.h' + str(gcIDsA) + '.sf' + str(fix(sizeFac*10)) + '.sav'
  readFlag = 0
  foreach saveFilename,saveFilenames do if ~file_test(saveFilename) then readFlag = 1
  
  if readFlag then begin
    ; load u,nelec and calculate temperature
    u     = loadSnapshotSubset(sP=sPa,partType='gas',field='u')
    nelec = loadSnapshotSubset(sP=sPa,partType='gas',field='nelec')
    temp  = alog10(convertUtoTemp(u,nelec))
    u     = !NULL
    nelec = !NULL
  
    ; load gas positions and velocities
    pos = loadSnapshotSubset(sP=sPa,partType='gas',field='pos')
    vel = loadSnapshotSubset(sP=sPa,partType='gas',field='vel')
    
    ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
    print,'shuffling...'
    iseed = 424242L
    sort_inds = sort(randomu(iseed,n_elements(temp)))
    
    temp = temp[sort_inds]
    pos  = pos[*,sort_inds]
    vel  = vel[*,sort_inds]
    
    sort_inds = !NULL
    
    print,'cutout...'
    foreach gcIDa,gcIDsA,k do begin
      ; get subhalo position and size of imaging box
      boxCen     = sgcen[*,gcIDa]
      boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcIDa]] / 10.0) * 10.0
      boxSizeImg = [boxSize,boxSize,boxSize] ; cube
    
      ; make conservative cutout greater than boxsize accounting for periodic (do cube not sphere)
      xDist = pos[0,*] - boxCen[0]
      yDist = pos[1,*] - boxCen[1]
      zDist = pos[2,*] - boxCen[2]
      
      correctPeriodicDistVecs, xDist, sP=sPa
      correctPeriodicDistVecs, yDist, sP=sPa
      correctPeriodicDistVecs, zDist, sP=sPa
      
      rvir = gc.group_r_crit200[gc.subgroupGrNr[gcIDg]]
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
      
      ; save
      save,loc_pos,loc_temp,loc_pos2,sPa,gcIDa,sizeFac,boxSizeImg,filename=saveFilenames[k]
      print,'Saved: '+strmid(saveFilenames[k],strlen(sPa.derivPath))
    endforeach
  
  endif ;readFlag

  print,'rendering arepo...'
  ; initial position: 3rd row left space
  pos = [0.0,0.25,0.25,0.5]                 
                 
  foreach gcIDa,gcIDsA,k do begin
    restore,saveFilenames[k]
    
    ; create color index mapping
    colorinds = (loc_temp-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
    colorinds = fix(colorinds + 50.0) > 0 < 255 ;50-255  
  
    ; local (cold) cutout
    wCold = where(loc_temp le coldTempCut,nCutoutCold)
    loc_temp = !NULL
    
    loc_pos_cold   = loc_pos[*,wCold]
    loc_pos2_cold  = loc_pos2[*,wCold]
    colorinds_cold = colorinds[wCold]
  
    ; get box center (in terms of specified axes)
    axisPair   = axes[k]
    boxCenImg  = [sgcen[axisPair[0],gcIDa],sgcen[axisPair[1],gcIDa],sgcen[3-axisPair[0]-axisPair[1],gcIDa]]
    haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcIDa]] ;ckpc
    haloMass = codeMassToLogMsun(gc.subgroupMass[gcIDa])
    
    ; fix halo mass if we're using old (x2 bug) catalogs
    if sPa.run eq 'arepo' then haloMass = codeMassToLogMsun(0.5*gc.subgroupMass[gcIDa])
  
    config = {boxSizeImg:boxSizeImg,plotFilename:plotFilename,haloVirRad:haloVirRad,haloMass:haloMass,$
              axisPair:axisPair,sP:sPa,barMM:tempMinMax,barType:'1temp'}
    
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
    tvcircle,config.haloVirRad,0,0,cgColor('dark gray'),thick=0.6,/data
    
    ; scale bar
 ;   len = 100.0 ;ckpc
 ;   cgText,mean([config.boxCen[0]-config.boxSizeImg[0]/2.2,$
 ;                config.boxCen[0]-config.boxSizeImg[0]/2.2+len]),$
 ;          config.boxCen[1]+config.boxSizeImg[0]/2.3,$
 ;          string(len,format='(i3)')+' ckpc',alignment=0.5,color=cgColor('black')
 ;   cgPlot,[config.boxCen[0]-config.boxSizeImg[0]/2.2,config.boxCen[0]-config.boxSizeImg[0]/2.2+len],$
 ;          [config.boxCen[1]+config.boxSizeImg[1]/2.1,config.boxCen[1]+config.boxSizeImg[1]/2.1],$
 ;          color=cgColor('black'),/overplot
           
    ; mass weighted temperature
    cgPlot, /nodata, xMinMax, yMinMax, pos=[0.0,0.0,0.5,1.0], xs=5, ys=5, /noerase
    tv, sphmap.quant_out,0.5,0.0,/normal,xsize=0.5
    
    ; redshift and halo mass
 ;   cgText,0.99,0.96,"z = "+string(config.sP.redshift,format='(f3.1)'),alignment=1.0,$
 ;          color=cgColor('black'),/normal
 ;   cgText,0.99,0.92,"M = "+string(config.haloMass,format='(f4.1)'),alignment=1.0,$
 ;          color=cgColor('black'),/normal
             
  end_PS, pngResize=60;, /deletePS

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
  minLogMsun = 11.25
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
    ;sphMapHalos,sP=sPg,gcIDs=gcIDs_gadget;,/coldOnly
    ;scatterMapHalos,sP=sPg,gcIDs=gcIDs_gadget
    ;scatterMapHalosGasDM,sP=sPg,gcIDs=gcIDs_gadget

    valNames  = ['temp','density','pressure','radvel','entropy','metallicity','angmom']
    foreach valName,valNames do $
      hsv = haloShellValue(sP=sPg,partType='gas',valName=valName,subgroupIDs=gcIDs_gadget,/cutSubS)

  endif
  
  ; arepo group catalog
  if select eq 'arepo' then begin
    gca = loadGroupCat(sP=sPa,/skipIDs)
    priAIDs = gcIDList(gc=gca,select='pri')
    priAMasses = codeMassToLogMsun(gca.subgroupMass[priAIDs])
    w = where(priAMasses ge minLogMsun and priAMasses le maxLogMsun,countA)
    gcIDs_arepo = priAIDs[w]
    
    print,'Mapping ['+str(countA)+'] arepo halos above minLogMsun.'
    ;sphMapHalos,sP=sPa,gcIDs=gcIDs_arepo;,/coldOnly
    ;scatterMapHalos,sP=sPa,gcIDs=gcIDs_arepo
    ;scatterMapHalosGasDM,sP=sPa,gcIDs=gcIDs_arepo

    valNames  = ['temp','density','pressure','radvel','entropy','metallicity','angmom']
    foreach valName,valNames do $
      hsv = haloShellValue(sP=sPa,partType='gas',valName=valName,subgroupIDs=gcIDs_arepo,/cutSubS)

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
  
  redshift = 2.0
  
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
