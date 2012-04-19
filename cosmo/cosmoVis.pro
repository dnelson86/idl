; cosmoVis.pro
; cosmological boxes - 2d visualization
; dnelson apr.2012

; sphMapBox: run sph kernel density projection on whole box

pro sphMapBox, res=res, run=run, partType=partType

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
    
    hsml = calcHSML(pos,ndims=3,nNGB=nNGB,boxSize=boxSize[0])

    ; OR: load HSML from snapshot (only stored for gas)
    ;hsml = loadSnapshotSubset(sP.simPath,snapNum=snap,partType='gas',field='hsml',/verbose)
      
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

; sphDensityProjection(): (OLD) make density projection using SPH kernel (inspired by Mark's sphMap)
;                         NOTE: kernel coeffs only valid for 3D!

function sphDensityProjection, pos, hsml, mass, quantity=quantity, imgSize=imgSize, boxSize=boxSize,$
                               boxCen=boxCen, axis0=axis0, axis1=axis1, mode=mode, periodic=periodic,$
                               verbose=verbose

  print,'You should switch this to the calcSphMap C-routine.'
  stop

  ; config
  if not keyword_set(axis0) then axis0 = 0
  if not keyword_set(axis1) then axis1 = 1
  if not keyword_set(verbose) then verbose = 0
  
  if keyword_set(periodic) then begin
    print,'ERROR: PERIODIC not supported.'
    return,0
  endif
  
  if (mode ne 1 and mode ne 2 and mode ne 3) then begin
    print,'ERROR: Unsupported mode='+str(mode)+' parameter.'
    return,0
  endif
  
  ; storage
  p    = dblarr(3)
  pos0 = double(0.0)
  pos1 = double(0.0)
  binnedParticles = 0UL
  
  ; init
  npart = n_elements(hsml)

  grid = fltarr(imgSize[0],imgSize[1])
  
  if keyword_set(quantity) then $
    gridQuantity = fltarr(imgSize[0],imgSize[1])
  
  pxSize = [float(boxSize[0]) / imgSize[0], float(boxSize[1]) / imgSize[1]]
  pxArea = pxSize[0] * pxSize[1]

  if (pxSize[0] lt pxSize[1]) then $
    hMin = 1.001 * pxSize[0] / 2.0
  if (pxSize[0] ge pxSize[1]) then $
    hMin = 1.001 * pxSize[1] / 2.0
    
  hMax = pxSize[0] * 50.0
  
  for part=0, npart-1, 1 do begin
    ; progress report
    if (part mod round(npart/10.0) eq 0 and verbose) then $
      print,'Progress: '+string(100.0*part/npart,format='(I3)')+'%'
      
    ; get particle data
    p[0] = pos[0,part]
    p[1] = pos[1,part]
    p[2] = pos[2,part]
    h    = double(hsml[part])
    v    = double(mass[part])
    
    if keyword_set(quantity) then $
      w    = double(quantity[part])
    
    ; early exit if out of z-bounds
    if (abs(p[3-axis0-axis1] - boxCen[2]) gt boxSize[2] / 2.0) then $
      continue
      
    pos0 = p[axis0] - (boxCen[0] - boxSize[0] / 2.0)
    pos1 = p[axis1] - (boxCen[1] - boxSize[1] / 2.0)
    
    ; clamp hsml
    if (h lt hMin) then h = hMin;
    if (h gt hMax) then h = hMax;
    
    ; early exit if ...
    if (pos0 - 0.0 lt -h or pos1 - 0.0 lt -h or pos0 - boxSize[0] gt h or pos1 - boxSize[1] gt h) then $
      continue
      
    binnedParticles += 1
    
    h2 = h * h;
    
    ; number of pixels covered by particle
    nx = h / pxSize[0] + 1;
    ny = h / pxSize[1] + 1;
    
    ; coordinates of pixel center of particle
    x = (floor(pos0 / pxSize[0]) + 0.5) * pxSize[0]
    y = (floor(pos1 / pxSize[1]) + 0.5) * pxSize[1]
    
    ; normalization constant
    sum = 0.0
    
    for dx = -nx, nx, 1 do begin
      for dy = -ny, ny, 1 do begin
        ; dist of covered pixel from actual position
        xx = x + dx * pxSize[0] - pos0
        yy = y + dy * pxSize[1] - pos1
        r2 = xx*xx + yy*yy
        
        if (r2 < h2) then begin
          ; sph kernel (inlined): sum += _getkernel(h,r2);
          hinv = double(1.0) / h
          u    = sqrt(r2) * hinv
          
          if (u lt 0.5) then begin
            sum += (2.546479089470 + 15.278874536822 * (u - 1.0) * u * u)
          endif else begin
            sum += (5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u))
          endelse
        endif ;r2 < h2
      endfor
    endfor
    
    ; exit if negligible
    if (sum lt 1.0e-10) then $
      continue
      
    ; add contribution to image
    for dx = -nx, nx, 1 do begin
      for dy = -ny, ny, 1 do begin
        ; coordinates of pixel center of covering pixels
        xxx = x + dx * pxSize[0]
        yyy = y + dy * pxSize[1]
        
        ; pixel array indices
        i = floor(xxx / pxSize[0]) ;implicit C cast to int
        j = floor(yyy / pxSize[1]) ;same
        
        if (i ge 0 and i lt imgSize[0] and j ge 0 and j lt imgSize[1]) then begin
          xx = x + dx * pxSize[0] - pos0
          yy = y + dy * pxSize[1] - pos1
          r2 = xx*xx + yy*yy
          
          if (r2 lt h2) then begin
            ; divide by sum for normalization
            ; divide by pixelarea to get column density (optional: /pxArea)
            ; sph kernel (inlined): grid[] += _getkernel(h,r2) * v / sum
            hinv = double(1.0) / h
            u    = sqrt(r2) * hinv
            
            if (u lt 0.5) then begin
              grid[i * imgSize[1] + j] += $
                (2.546479089470 + 15.278874536822 * (u - 1.0) * u * u) * v / sum
              if keyword_set(quantity) then $
                gridQuantity[i * imgSize[1] + j] += $
                  (2.546479089470 + 15.278874536822 * (u - 1.0) * u * u) * v * w / sum
            endif else begin
              grid[i * imgSize[1] + j] += $
                (5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u)) * v / sum
                  if keyword_set(quantity) then $
                  gridQuantity[i * imgSize[1] + j] += $
                  (5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u)) * v * w / sum
            endelse
          
          endif ;r2 < h2
        endif ;i,j
      
      endfor
    endfor

  endfor ;part
  
  if (verbose) then print,'Number of binned particles: ',binnedParticles
  
  if (mode eq 1) then begin
    if (verbose) then print,'Returning: Column Mass Map'
    return,grid
  endif
  if (mode eq 2) then begin
    if (verbose) then print,'Returning: Quantity Mass-Weighted Map'
    return,gridQuantity
  endif
  if (mode eq 3) then begin
    if (verbose) then print,'Returning: Column Density Map'
    for i=0,i lt imgSize[0] do begin
      for j=0,j lt imgSize[1] do begin
        grid[i + imgSize[1] * j] /= pxArea
      endfor
    endfor
    
    return,grid
  endif

end

; scatterMapHalos: plot temperature colored scatter plots with velocity vectors on boxes centered on halos

pro scatterMapHalos, sP=sP, gcIDs=gcIDs

  if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  sizeFac = 3.5       ; times rvir
  cutFac  = 1.0       ; times boxSize
  nPixels = [800,800] ; px
  
  tempMinMax  = [4.0,7.0] ; log(K)
  coldTempCut = 5.0       ; log(K)
  velVecFac   = 0.001     ; times velocity (km/s) in plotted kpc
  
  axes = list([0,1],[0,2],[1,2]) ;xy,xz,yz
  ;axes = list([0,1])
  
  ; target list
  print,'loading...'
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
  
  ; randomly shuffle the points (break the peano orderin to avoid "square" visualization artifacts)
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
    xDist = boxCen[0] - pos[0,*]
    yDist = boxCen[1] - pos[1,*]
    zDist = boxCen[2] - pos[2,*]
    
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
    colorinds = fix(colorinds + 50.0) ;50-255  
  
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
      
      print,'['+string(gcID,format='(i4)')+'] Mapping ['+str(axisPair[0])+' '+$
            str(axisPair[1])+'] with '+str(boxSize[0])+$
            ' kpc box around subhalo center ['+str(boxCen[0])+' '+str(boxCen[1])+' '+str(boxCen[2])+']'

      xMinMax = [-boxSizeImg[0]/2.0,boxSizeImg[0]/2.0]
      yMinMax = [-boxSizeImg[1]/2.0,boxSizeImg[1]/2.0]
      
      plotFilename = 'scatter.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+$
                     '.axes'+str(axisPair[0])+str(axisPair[1])+'.eps'
      
      start_PS, sP.plotPath + plotFilename, xs=8, ys=4
      
        ; color table and establish temperature -> color mapping
        loadColorTable,'helix',/reverse
        
        TVLCT, rr, gg, bb, /GET

        newcolors = getColor24([[rr[colorinds]], [gg[colorinds]], [bb[colorinds]]])
        newcolors_cold = getColor24([[rr[colorinds_cold]], [gg[colorinds_cold]], [bb[colorinds_cold]]])

        ; all gas
        cgPlot, /nodata, xMinMax, yMinMax, pos=[0.0,0.0,0.5,1.0], xs=5, ys=5
        
        haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcID]] ;ckpc
        tvcircle,haloVirRad,0,0,cgColor('dark gray'),thick=0.6,/data
        
        ; particle loop for velocity vector plotting
        for i=0L,nCutout-1 do $
          plots,reform([loc_pos[axisPair[0],i],loc_pos2[axisPair[0],i]]),$
                 reform([loc_pos[axisPair[1],i],loc_pos2[axisPair[1],i]]),$
                 line=0,color=newcolors[i]

        ; scale bar
        len = 100.0 ;ckpc
        cgText,mean([-boxSizeImg[0]/2.2,-boxSizeImg[0]/2.2+len]),boxSizeImg[0]/2.4,$
               string(len,format='(i3)')+' ckpc',alignment=0.5,charsize=!p.charsize-0.6,color=cgColor('dark gray')
        cgPlot,[-boxSizeImg[0]/2.2,-boxSizeImg[0]/2.2+len],[boxSizeImg[1]/2.2,boxSizeImg[1]/2.2],$
               color=cgColor('dark gray'),/overplot
               
        ; cold only
        cgPlot, /nodata, xMinMax, yMinMax, pos=[0.5,0.0,1.0,1.0], xs=5, ys=5, /noerase

        tvcircle,haloVirRad,0,0,cgColor('dark gray'),thick=0.6,/data
        
        ; particle loop for velocity vector plotting (cold gas only)
        for i=0L,nCutoutCold-1 do $
          plots,reform([loc_pos_cold[axisPair[0],i],loc_pos2_cold[axisPair[0],i]]),$
                 reform([loc_pos_cold[axisPair[1],i],loc_pos2_cold[axisPair[1],i]]),$
                 line=0,color=newcolors_cold[i]
        
        ; redshift and halo mass
        !p.charsize = 1.0
        haloMass = codeMassToLogMsun(gc.subgroupMass[gcID])
        cgText,0.98,0.94,"z = "+string(sP.redshift,format='(f3.1)'),alignment=1.0,color=cgColor('dark gray'),/normal
        cgText,0.98,0.89,"M = "+string(haloMass,format='(f4.1)'),alignment=1.0,color=cgColor('dark gray'),/normal
        
        ; dividing line
        cgPlot,[xMinMax[0],xMinMax[0]],yMinMax,line=0,thick=1.0,color=cgColor('light gray'),/overplot
        
        ; temp colorbar on right
        !x.thick = 1.0
        !y.thick = 1.0
        colorbar,position=[0.02,0.35,0.08,0.65],divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
        cgText,0.5,0.04,textoidl("log T_{gas} [K]"),alignment=0.5,charsize=!p.charsize-0.2,color=cgColor('white'),/normal
        cgText,0.365,0.04,'4',alignment=0.5,charsize=!p.charsize-0.2,color=cgColor('white'),/normal
        cgText,0.635,0.04,'7',alignment=0.5,charsize=!p.charsize-0.2,color=cgColor('white'),/normal
        
      end_PS, pngResize=60, im_options='-negate', /deletePS

    endforeach ;axisPair

  endforeach ;gcIDs

end

; sphMapHalos: run sph kernel density projection on boxes centered on halos

pro sphMapHalos, sP=sP, gcIDs=gcIDs, coldOnly=coldOnly

  if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  sizeFac = 2.0       ; times rvir
  cutFac  = 1.1       ; times boxSize
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
    hsml = 1.25 * temporary(hsml); increase hsml to decrease visualization noise
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
    boxCen = sgcen[*,gcID]
    boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcID]] / 10.0) * 10.0
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
        config = {sizeFac:sizeFac,cutFac:cutFac,nPixels:nPixels,axes:axes,gcID:gcID,boxCen:boxCen,$
                  boxCenImg:boxCenImg,boxSize:boxSize,boxSizeImg:boxSizeImg,sP:sP}
        save,sphmap,config,filename=sP.derivPath+'sphMaps/'+saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
      endif else begin
        restore,sP.derivPath +'sphMaps/'+ saveFilename
      endelse
    
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
      
      start_PS, sP.plotPath + strmid(saveFilename,0,strlen(saveFilename)-4)+'.eps', xs=8, ys=4
      
        ;loadct, 4, bottom=1, /silent
        loadColorTable,'helix'
        ; density
        cgPlot, /nodata, xMinMax, yMinMax, pos=[0.0,0.0,0.5,1.0], xs=5, ys=5
        tv, sphmap.dens_out,0.0,0.0,/normal,xsize=0.5
        
        ; scale bar
        len = 100.0 ;ckpc
        cgText,mean([config.boxCen[0]-config.boxSizeImg[0]/2.2,$
                     config.boxCen[0]-config.boxSizeImg[0]/2.2+len]),$
               config.boxCen[1]+config.boxSizeImg[0]/2.4,$
               string(len,format='(i3)')+' ckpc',alignment=0.5,charsize=!p.charsize-0.6,color=cgColor('black')
        cgPlot,[config.boxCen[0]-config.boxSizeImg[0]/2.2,config.boxCen[0]-config.boxSizeImg[0]/2.2+len],$
               [config.boxCen[1]+config.boxSizeImg[1]/2.2,config.boxCen[1]+config.boxSizeImg[1]/2.2],$
               color=cgColor('black'),/overplot
               
        ; mass weighted temperature
        cgPlot, /nodata, xMinMax, yMinMax, pos=[0.0,0.0,0.5,1.0], xs=5, ys=5, /noerase
        tv, sphmap.quant_out,0.5,0.0,/normal,xsize=0.5
        
        ; redshift and halo mass
        !p.charsize = 1.0
        haloMass = codeMassToLogMsun(gc.subgroupMass[gcID])
        cgText,0.98,0.94,"z = "+string(sP.redshift,format='(f3.1)'),alignment=1.0,color=cgColor('black'),/normal
        cgText,0.98,0.89,"M = "+string(haloMass,format='(f4.1)'),alignment=1.0,color=cgColor('black'),/normal
                 
      end_PS, pngResize=60, /deletePS

    endforeach ;axisPair
  
  endforeach ;gcIDs
  
end

; makeHaloComparisonImages(): create a mosaic of halo comparison images between gadget and arepo

pro makeHaloComparisonImages, select=select

  ; config
  res      = 512
  redshift = 0.0
  
  minLogMsun = 11.0
  
  sPa = simParams(res=res,run='arepo',redshift=redshift)
  sPg = simParams(res=res,run='gadget',redshift=redshift)
  
  ; gadget group catalog
  if select eq 'gadget' then begin
    gcg = loadGroupCat(sP=sPg,/skipIDs)
    priGIDs = gcIDList(gc=gcg,select='pri')
    priGMasses = codeMassToLogMsun(gcg.subgroupMass[priGIDs])
    w = where(priGMasses ge minLogMsun,countG)
    gcIDs_gadget = priGIDs[w]
    
    gcIDs_gadget = [10233,10264,10288,10311,10348,10384,10403,10439,10472,10487,10503,10526,$
                    10544,10577,10615,10649,10666,10686,10706,10727,10757,10786,10800,10826,$
                    10865,10881,10889,10936,10953,10980,10994,11020,11041,11061,11074,11088,$
                    11109,11122,11155,11185]
    
    print,'Mapping ['+str(countG)+'] gadget halos above minLogMsun.'
    sphMapHalos,sP=sPg,gcIDs=gcIDs_gadget;,/coldOnly
    ;scatterMapHalos,sP=sPg,gcIDs=gcIDs_gadget
  endif
  
  ; arepo group catalog
  if select eq 'arepo' then begin
    gca = loadGroupCat(sP=sPa,/skipIDs)
    priAIDs = gcIDList(gc=gca,select='pri')
    priAMasses = codeMassToLogMsun(gca.subgroupMass[priAIDs])
    w = where(priAMasses ge minLogMsun,countA)
    ;gcIDs_arepo = priAIDs[w]
    
    gcIDs_arepo = [9249,9270,9275,9291,9312,9327,9350,9386,9400,9417,9427,9442,9466,9480,9481,9490,$
                   9504,9522,9523,9531,9543,9559,9588,9604,9620,9635,9662,9679,9694,9700,9716,9739,$
                   9759,9772,9782,9787,9798,9820,9831,9846,9868,9889,9907,9919,9934,9946,9961]
    
    print,'Mapping ['+str(countA)+'] arepo halos above minLogMsun.'
    ;sphMapHalos,sP=sPa,gcIDs=gcIDs_arepo;,/coldOnly
    scatterMapHalos,sP=sPa,gcIDs=gcIDs_arepo
  endif
  
  print,'done.'

end

pro makeHaloComparisonPage
  
  ; config
  res = 512
  redshifts = [0.0,0.5,1.0,2.0]
  
  massBins = [10.5,11.0,11.5,12.0,14.0]
  
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
  
  ; match only for webpage output (maybe redo with higher tols)
  match = findMatchedHalos(sP1=sPa,sP2=sPg)  
  
  ; loop over mass bins
  for j=0,n_elements(massBins)-2 do begin
    print,massBins[j],massBins[j+1]
    
    ; open file and write header with links
    outputName = "HaloComp_z"+string(redshift,format='(f3.1)')+"mb"+string(massBins[j],format='(f4.1)')+".htm"
    openw,lun,outputName,/get_lun
    
    printf,lun,"<div class='compHead'><table>"
    printf,lun,"  <tr><th>Redshifts</th><th colspan='4'>Log(Msun) Bins</th></tr>"
    foreach red,redshifts,k do begin
      printf,lun,"  <tr><td>"+string(red,format='(f3.1)')+"</td>"
      for mb=0,n_elements(massBins)-2 do begin
        printf,lun,"    <td><a href='HaloComp_z"+string(red,format='(f3.1)')+"mb"+$
                   string(massBins[mb],format='(f4.1)')+".htm'>"+string(massBins[mb],format='(f4.1)')+$
                   " - "+string(massBins[mb+1],format='(f4.1)')+"</a></td>"
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
    matchedInds = match.matchedInds[wA]
    wMatch = where(matchedInds ne -1,countMatch,comp=wNoMatch,ncomp=countNoMatchA)
    
    ; find subset of gadget halos that were not matched
    match,match.matchedInds,wG,inds_match,inds_wG,count=countNoMatchG,/sort
    if countNoMatchG gt 0 then wGNoMatch = wG[inds_wG]
    
    ; mass bin header
    printf,lun,""
    printf,lun,"<a name='z"+str(k)+"_massbin"+str(j)+"'></a>"
    printf,lun,"<h3>"+string(massBins[j],format='(f4.1)')+" < log(M) < "+$
                      string(massBins[j+1],format='(f4.1)')+"</h3>"
    printf,lun,""
    
    ; TEMPORARY: skip lowest mass bin (not run)
    if j eq 0 then begin
      printf,lun,"(Not yet rendered)."
      printf,lun,"</div>"
      printf,lun,""
      printf,lun,"</body>"
      printf,lun,"</html>"
       
      ; close file
      free_lun,lun
      
      ; if header.txt exists, cat the two together
      if file_test("header.txt") then begin
        print,' added header (skipped lowest massbin)'
        spawn,'cat header.txt '+outputName+' > newpage.htm' ;make new
        spawn,'mv newpage.htm '+outputName ;overwrite original
      endif
      
      continue
    endif
    
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
      
      statsString = "<span class='haloid'>ID: "+str(haloCounter)+"</span><br><br> "+$
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
        
        gadgetString += axesNames[m]+": <a href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+"'></a>"
                        
        ; scatterMap
        fpath = "scatterMaps.gadget.z"+string(redshift,format='(f3.1)')+"/"
        fname = "scatter.G"+str(res)+"."+str(sPg.snap)+".h"+str(gcIndG)+".axes"+axisPair+".png"
        
        gadgetString += " <a href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+"'></a><br>"
      endforeach
      
      gadgetString = strmid(gadgetString,0,strlen(gadgetString)-4)
      
      ; make arepo image thumbnail/link string
      arepoString = ''
      
      foreach axisPair,axesPairs,m do begin        
        ; sphMap
        fpath = "sphMaps.arepo.z"+string(redshift,format='(f3.1)')+"/"
        fname = "map.A"+str(res)+"."+str(sPa.snap)+".h"+str(gcIndA)+".axes"+axisPair+".png"
        
        arepoString += axesNames[m]+": <a href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+"'></a>"
                        
        ; scatterMap
        fpath = "scatterMaps.arepo.z"+string(redshift,format='(f3.1)')+"/"
        fname = "scatter.A"+str(res)+"."+str(sPa.snap)+".h"+str(gcIndA)+".axes"+axisPair+".png"
        
        arepoString += " <a href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+"'></a><br>"
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
      
      statsString = "<span class='haloid'>ID: "+str(haloCounter)+"</span><br><br> "+$
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
        
        arepoString += axesNames[m]+": <a href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+"'></a>"
                        
        ; scatterMap
        fpath = "scatterMaps.arepo.z"+string(redshift,format='(f3.1)')+"/"
        fname = "scatter.A"+str(res)+"."+str(sPa.snap)+".h"+str(gcIndA)+".axes"+axisPair+".png"
        
        arepoString += " <a href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+"'></a><br>"
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
      haloCounter += 1
      ; make stats string
      gcIndG = priGIDs[wGNoMatch[i]]
      
      haloMassG = codeMassToLogMsun(gcg.subgroupMass[gcIndG])
      virRadG   = gcg.group_r_crit200[gcg.subgroupGrNr[gcIndG]]
      virTempG  = alog10(codeMassToVirTemp(gcg.subgroupMass[gcIndG],sP=sPg))
      
      statsString = "<span class='haloid'>ID: "+str(haloCounter)+"</span><br><br> "+$
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
        
        gadgetString += axesNames[m]+": <a href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+"'></a>"
                        
        ; scatterMap
        fpath = "scatterMaps.gadget.z"+string(redshift,format='(f3.1)')+"/"
        fname = "scatter.G"+str(res)+"."+str(sPg.snap)+".h"+str(gcIndG)+".axes"+axisPair+".png"
        
        gadgetString += " <a href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+"'></a><br>"
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
