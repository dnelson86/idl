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

; sphMapSubhalos: run sph kernel density projection on boxes centered on halos/subhalos

pro sphMapHalos, sP=sP, gcIDs=gcIDs

  ; config
  boxSize = [200,200,200] ;kpc
  imgSize = [800,800]     ;px
  
  axes = list([0,1]) ;x,y
  mode = 1 ;1=col mass, 2=mass-weighted quantity, 3=col density
  
  ; target list
  gcTarget = loadGroupCat(sP=sP,/verbose)
    
  if (not keyword_set(gcIDs)) then begin
    ; make list of primary halo IDs (map them all)
    valGCids = gcIDList(gc=gcTarget,select='pri')
  endif else begin
    ; make specified IDs
    valGCids = gcIDs
  endelse

  ; load properties from snapshot
  pos  = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  hsml = loadSnapshotSubset(sP=sP,partType='gas',field='hsml')
  mass = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
  
  ; loop over all non-background subhalos and image
  foreach gcID, valGCids do begin
  
    foreach axisPair, axes do begin
  
      imgFilename = 'sphmap.G.snap='+str(sP.snap)+'.gcID='+str(gcID)+'.axis0='+$
                    str(axisPair[0])+'.axis1='+str(axisPair[1])+'.res='+str(sP.res)+'.box='+str(boxSize[0])
                    
      if (file_test(sP.plotPath+imgFilename+'.png') or file_test(sP.plotPath+imgFilename+'.eps')) then begin
        print,'Skipping: ' + imgFilename
        continue
      endif     
  
      ; get subhalo position
      boxCen = gcTarget.subgroupPos[*,gcID]
      
      print,'['+string(gcID,format='(I04)')+'] Mapping ['+str(axisPair[0])+' '+$
            str(axisPair[1])+'] with '+str(boxSize[0])+$
            ' kpc box around subhalo center ['+str(boxCen[0])+' '+str(boxCen[1])+' '+str(boxCen[2])+']'
    
      colMassMap = sphDensityProjection(pos, hsml, mass, imgSize=imgSize, boxSize=boxSize,$
                                        boxCen=boxCen, axis0=axisPair[0], axis1=axisPair[1], $
                                        mode=mode, periodic=0, /verbose)
    
      ; rescale
      w = where(colMassMap eq 0, count, complement=ww)
      if (count ne 0) then colMassMap[w] = min(colMassMap[ww])
      
      colMassMap = alog10(colMassMap)
      colMassMap = (colMassMap-min(colMassMap))*254.0 / (max(colMassMap)-min(colMassMap)) ;0-254
      ;h2d = filter_image(h2d,FWHM=[1.1,1.1],/ALL) ;gaussian kernel convolution
      colMassMap += 1.0 ;1-254
    
      ; plot PS and PNG
      xMinMax = [boxCen[0]-boxSize[0]/2.0,boxCen[0]+boxSize[0]/2.0]
      yMinMax = [boxCen[1]-boxSize[1]/2.0,boxCen[1]+boxSize[1]/2.0]
      
      start_PS, sP.plotPath+imgFilename+'.eps'
        loadct, 4, bottom=1, /silent
        tvim,colMassMap,xrange=xMinMax,yrange=yMinMax,pcharsize=0.0001
        fsc_text,0.72,0.05,"z=3 id="+string(gcID,format='(I04)'),alignment=0.5,$
                 color=fsc_color('white'),/normal
      end_PS, pngResize=68, /deletePS
    
    endforeach ;axisPair
  
  endforeach ;valGCids
  
end

; cosmoTrajMovieMatchedHalo(): load positions(time) and make interpolated images of tracer/sph positions
;                              from a cosmological box via accretionTraj() comparing a matched halo

pro cosmoTrajMovieMatchedHalo

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  sPa = simParams(res=128,run='tracer',redshift=3.0)
  sPa.trMCPerCell = -1 ; velocity tracers
  sPg = simParams(res=128,run='gadget',redshift=3.0)

  ; get gas/tracer positions with time
  at_ga  = accretionTraj(sP=sPg)
  mt_ga  = mergerTreeSubset(sP=sPg)
  at_vel = accretionTraj(sP=sPa)
  mt_vel = mergerTreeSubset(sP=sPa)

  ; plot a single halo? 
  matchID = 1
  
  if matchID ne -1 then begin
    print,'matching halos...'
    ; get the matched group catalog
    match = findMatchedHalos(sP1=sPa,sP2=sPg)
  
    haloIDs_vel = []
    haloIDs_ga  = []
  
    ; find intersection of matched halos with tracked halos in accretionTimes()
    for i=0,match.nMatched-1 do begin
      haloID_vel = match.wMatch[i]
      haloID_ga = match.matchedInds[haloID_vel]
      
      ; this arepo halo is tracked?
      w = where(mt_vel.galcatIDList eq haloID_vel,count)
      if count eq 0 then continue
      
      ; the matched gadget halo is tracked?
      w = where(mt_ga.galcatIDList eq haloID_ga,count)
      if count eq 0 then continue
      
      haloIDs_vel = [haloIDs_vel,haloID_vel]
      haloIDs_ga  = [haloIDs_ga,haloID_ga]
    endfor
      
    ; load group catalog to find masses of these matched halos
    gc = loadGroupCat(sP=sPa,/skipIDs)
    haloMasses = codeMassToLogMsun(gc.subgroupMass[haloIDs_vel])
    gc = !NULL
    
    ; choose halo
    haloID_vel = haloIDs_vel[matchID]
    haloID_ga  = haloIDs_ga[matchID]
    haloMass   = haloMasses[matchID]

    ; load child counts for all members of galcat
    maxt_vel_gal = maxTemps(sP=sPa,/loadAllTRGal)
    cc_vel_gal   = maxt_vel_gal.child_counts
    maxt_vel_gal = !NULL
    
    maxt_vel_gmem = maxTemps(sP=sPa,/loadAllTRGmem)
    cc_vel_gmem   = maxt_vel_gmem.child_counts
    maxt_vel_gmem = !NULL
    
    ;maxt_ga_gal = maxTemps(sP=sPg,/loadByGas)
    ;cc_ga_gal   = maxt_ga_gal.child_counts
    ;maxt_ga_gal = !NULL
    
    ; get indices for the subset of the galaxycat corresponding to the matched halo in each sim
    galcatSub_vel = galCatINDList(sP=sPa,gcIDList=[haloID_vel],$
                               child_counts={gal:cc_vel_gal,gmem:cc_vel_gmem})          
    galcatSub_ga  = galCatINDList(sP=sPg,gcIDList=[haloID_ga])
  endif else begin
    ; make pseudo galcatSub which accesses all members of galaxycat
    galcatSub_vel = { gal  : lindgen(n_elements(at_vel.relPos_gal[0,0,*]))  ,$
                      gmem : lindgen(n_elements(at_vel.relPos_gmem[0,0,*]))  }
    galcatSub_ga  = { gal  : lindgen(n_elements(at_ga.relPos_gal[0,0,*]))  ,$
                      gmem : lindgen(n_elements(at_ga.relPos_gmem[0,0,*]))  }        
  endelse

  ; count
  nTrMC  = 0
  nTrVel = n_elements(at_vel.relPos_gal[0,0,galcatSub_vel.gal])
  nSph   = n_elements(at_ga.relPos_gal[0,0,galcatSub_ga.gal])

  ; view configuration
  boxCen   = [0,0]
  zoomSize = 400.0 ;kpc
  axes     = [0,2]
  
  ; movie configuration
  nFrames = 401
  timeStart = 4.0 ;redshift
  timeEnd   = 3.0 ;redshift
  
  ; convert scale factors to redshifts (or TODO: ages)
  mt_ga.times  = 1/mt_ga.times-1
  mt_vel.times = 1/mt_vel.times-1
  
  if timeStart lt min([mt_ga.times,mt_vel.times]) then $
    message,'Error: Requested timeStart beyond snapshot range.'
  if timeEnd gt max([mt_ga.times,mt_vel.times]) then $
    message,'Error: Requested timeEnd beyond snapshot range.'
  
  ; trail length and color sequence
  trailLenTime = 0.05 ;Gyr

  ; deriv
  timeStep = (timeEnd - timeStart) / (nFrames-1)
  frameTimes = timeStep * findgen(nFrames) + timeStart
  frameTimesRev = reverse(frameTimes)
  
  trailLenFrame = abs(fix(trailLenTime / timeStep))
  
  ; setup fading trail
  rgb = fix(findgen(trailLenFrame)/(trailLenFrame-1) * 254)
  trailColors = lonarr(trailLenFrame)
  for i=0,trailLenFrame-1 do $
    trailColors[i] = rgb[i] + rgb[i]*(256L) + rgb[i]*(256L)^2

  trailColors = reverse(trailColors)

  ; spline interpolation of halo properties
  if matchID ne -1 then begin
    print,'interp halo...'
    
    ; reverse halo properties so they progress forward in time
    velHaloMass = reverse(codeMassToLogMsun(spline(mt_vel.times,mt_vel.hMass[*,haloID_vel],frameTimesRev)))
    gaHaloMass  = reverse(codeMassToLogMsun(spline(mt_ga.times,mt_ga.hMass[*,haloID_ga],frameTimesRev)))

    velHaloRvir = reverse(spline(mt_vel.times,mt_vel.hVirRad[*,haloID_vel],frameTimesRev))
    gaHaloRvir  = reverse(spline(mt_ga.times,mt_ga.hVirRad[*,haloID_ga],frameTimesRev))
  endif

  ; spline interpolation of axes[0],axes[1] coordinates
  velPos = fltarr(2,nTrVel,nFrames)
  gaPos  = fltarr(2,nSph,nFrames)

  print,'interp vel...'
  for i=0,nTrVel-1 do begin
    velPos[0,i,*] = spline(mt_vel.times,at_vel.relPos_gal[*,axes[0],galcatSub_vel.gal[i]],frameTimesRev)
    velPos[1,i,*] = spline(mt_vel.times,at_vel.relPos_gal[*,axes[1],galcatSub_vel.gal[i]],frameTimesRev)
  ;  ;velPos[0,i,*] = interpol(pos.vel[*,axes[0],i],pos.times,frameTimesRev)
  ;  ;velPos[1,i,*] = interpol(pos.vel[*,axes[1],i],pos.times,frameTimesRev)
  endfor
  
  print,'interp sph...'
  for i=0,nSph-1 do begin
    ; spline interpolate to reversed frameTimes (so they progress back in time like the data)
    gaPos[0,i,*] = spline(mt_ga.times,at_ga.relPos_gal[*,axes[0],galcatSub_ga.gal[i]],frameTimesRev)
    gaPos[1,i,*] = spline(mt_ga.times,at_ga.relPos_gal[*,axes[1],galcatSub_ga.gal[i]],frameTimesRev)
    ;gaPos[0,i,*] = interpol(pos.ga[*,axes[0],i],pos.times,frameTimesRev)
    ;gaPos[1,i,*] = interpol(pos.ga[*,axes[1],i],pos.times,frameTimesRev)
  endfor
  
  ; reverse positions so they progress forward in time
  gaPos  = reverse(gaPos,3,/overwrite)
  velPos = reverse(velPos,3,/overwrite)

  ; clip tracers at box dimensions
  w = where(gaPos lt boxCen[0]-zoomSize/2.0,count)
  if count gt 0 then gaPos[w] = !values.f_nan
  
  w = where(velPos gt boxCen[0]+zoomSize/2.0,count)
  if count gt 0 then velPos[w] = !values.f_nan
  w = !NULL

  print,'rendering frames...'
  for fn=0,nFrames-1 do begin
  ;fn = 100
    ; determine time of this frame and bracketing snapshots
    time = timeStart + timeStep * fn
    
    ; determine earliest frame number to draw trail to
    fnTrail = fn - trailLenFrame + 1 < fn > 0
    numTrail = fn - fnTrail + 1

    print,' ['+string(fn,format='(i3)')+'] z = '+string(time,format='(f5.3)')+' ( '+str(fn)+' - '+str(fnTrail)+' )'
    
    ; plot
    start_PS, sPa.plotPath + 'frames256/frame_'+str(fn)+'.eps', xs=8, ys=4
    
      !p.thick = 1.0
      
      ; determine halo rvir at this time
      haloRvir_ga  = 200.0
      haloRvir_vel = 200.0
      
      if matchID ne -1 then haloRvir_ga  = gaHaloRvir[fn]
      if matchID ne -1 then haloRvir_vel = velHaloRvir[fn]
      
      ; fixed viewport
      xrange = [boxCen[0]-zoomSize/2.0,boxCen[0]+zoomSize/2.0]
      yrange = [boxCen[1]-zoomSize/2.0,boxCen[1]+zoomSize/2.0] ; top+bottom halves
    
      ; left: arepo
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=5,ystyle=5,$
             xtitle="",ytitle="",title="",position=[0,0,0.5,1],xtickname=replicate(' ',10),$
             ytickname=replicate(' ',10)
           
      ; left: draw rvir and rvir/10
      tvcircle,haloRvir_vel/10,0,0,cgColor('green'),thick=0.5,/data
      tvcircle,haloRvir_vel/100,0,0,cgColor('green'),thick=0.5,/data
      tvcircle,haloRvir_vel,0,0,cgColor('green'),thick=0.5,/data
        
      ; to plot the trail, plots each particle path separately
      for i=0,nTrVel-1 do $
        for j=fnTrail,fn-1 do $
          plots,velPos[0,i,j:(j+1)],velPos[1,i,j:(j+1)],color=trailColors[j-fnTrail]
        
      ; scale bar
      len = 100.0 ;ckpc
      cgPlot,[boxCen[0]-zoomSize/2.2,boxCen[0]-zoomSize/2.2+len],[boxCen[1]+zoomSize/2.2,boxCen[1]+zoomSize/2.2],$
        line=0,color=cgColor('dark gray'),/overplot
      cgText,mean([boxCen[0]-zoomSize/2.2,boxCen[0]-zoomSize/2.2+len]),boxCen[1]+zoomSize/2.4,$
        string(len,format='(i3)')+' ckpc',alignment=0.5,charsize=!p.charsize-0.6,color=cgColor('dark gray')
        
      ; right: gadget (tile x2 in the left-right direction)
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=5,ystyle=5,$
             xtitle="",ytitle="",title="",position=[0.5,0,1,1],xtickname=replicate(' ',10),$
             ytickname=replicate(' ',10),/noerase
        
      tvcircle,haloRvir_ga/10,0,0,cgColor('green'),thick=0.5,/data
      tvcircle,haloRvir_ga/100,0,0,cgColor('green'),thick=0.5,/data
      tvcircle,haloRvir_ga,0,0,cgColor('green'),thick=0.5,/data 

      for i=0,nSph-1 do $
        for j=fnTrail,fn-1 do $
          plots,gaPos[0,i,j:(j+1)],gaPos[1,i,j:(j+1)],color=trailColors[j-fnTrail]
    
      ; draw time ticker
      strTime = "z = "+string(time,format='(f5.3)')
      cgText,0.5,0.03,strTime,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
      
      ; halo mass
      if matchID ne -1 then begin
        strHalo = "log(M) = "+string(velHaloMass[fn],format='(f4.1)')
        cgText,0.5,0.08,strHalo,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
        ;strHalo = textoidl("R_{vir} = ")+string(velHaloRvir[fn],format='(i3)')
        ;cgText,0.5,0.13,strHalo,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
      endif
      
      ; simulation names
      cgText,0.05,0.03,"Arepo",alignment=0.5,color=cgColor('dark orchid'),/normal,charsize=!p.charsize-0.4
      cgText,0.95,0.03,"Gadget",alignment=0.5,color=cgColor('dark orchid'),/normal,charsize=!p.charsize-0.4
      
      ; dividing line
      cgPlot,[boxCen[0]-zoomSize/2.0,boxCen[0]-zoomSize/2.0],[boxCen[1]-zoomSize/2.4,boxCen[1]+zoomSize/2.0],$
        line=0,color=cgColor('light gray'),/overplot
      
    end_PS, pngResize=60, im_options='-negate', /deletePS

  endfor

end

; cosmoTrajMovieHalo(): load positions(time) and make interpolated images of tracer/sph positions
;                       from a cosmological box via accretionTraj() for one halo or a mass bin

pro cosmoTrajMovieHalo

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  sP = simParams(res=256,run='gadget',redshift=3.0)

  ; choose one: plot a single halo? or a mass bin of halos?
  massBinLog = [] ;[12.0,12.5]
  mtHaloInd  = 0
  
  ; view configuration
  boxCen   = [0,0]
  zoomSize = 400.0 ;kpc
  axes     = [0,2]
  
  ; movie configuration
  nFrames   = 401
  timeStart = 4.0 ;redshift
  timeEnd   = 3.0 ;redshift
  trailLenTime = 0.05 ;Gyr
  
  ; color config
  ggTintFac    = 1.5   ; color skew for gal,gmem,hot
  tempLogFac   = -0.05 ; within ~90% of Tvir or above
  
  ; get gas/tracer positions with time
  at  = accretionTraj(sP=sP)
  mt  = mergerTreeSubset(sP=sP)

  haloIDs = []

  ; load group catalog to find masses of tracked halos
  gc = loadGroupCat(sP=sP,/skipIDs)
  haloMasses = codeMassToLogMsun(gc.subgroupMass[mt.galcatIDList])
  gc = !NULL
  
  ; choose halo(s)
  if n_elements(massBinLog) gt 0 then begin
    w = where(haloMasses ge massBinLog[0] and haloMasses le massBinLog[1],count)
    if count eq 0 then message,'Error: No halos in mass bin.'
    print,'Found ['+str(count)+'] halos in mass bin.'
    
    haloIDs    = mt.galcatIDList[w]
    haloMasses = halomasses[w]
  endif else begin
    haloIDs    = [mt.galcatIDList[mtHaloInd]]
    haloMasses = haloMasses[mtHaloInd]
  endelse

  ; get indices for the subset of the galaxycat corresponding to the selected halo(s)
  if sP.trMCPerCell ne 0 then begin
    maxt_vel_gal = maxTemps(sP=sPa,/loadAllTRGal)
    cc_vel_gal   = maxt_vel_gal.child_counts
    maxt_vel_gal = !NULL
    
    maxt_vel_gmem = maxTemps(sP=sPa,/loadAllTRGmem)
    cc_vel_gmem   = maxt_vel_gmem.child_counts
    maxt_vel_gmem = !NULL
    
    galcatSub = galCatINDList(sP=sP,gcIDList=haloIDs,child_counts={gal:cc_vel_gal,gmem:cc_vel_gmem})
  endif else begin    
    galcatSub = galCatINDList(sP=sP,gcIDList=haloIDs)
  endelse
  
  ; count
  nPart = { gal  : n_elements(at.relPos_gal[0,0,galcatSub.gal])   ,$
            gmem : n_elements(at.relPos_gmem[0,0,galcatSub.gmem])  }

  ; convert scale factor to redshift
  mt.times  = 1/mt.times-1
  
  if timeStart lt min(mt.times) then message,'Error: Requested timeStart beyond snapshot range.'
  if timeEnd gt max(mt.times) then message,'Error: Requested timeEnd beyond snapshot range.'
  
  ; deriv
  timeStep = (timeEnd - timeStart) / (nFrames-1)
  frameTimes = timeStep * findgen(nFrames) + timeStart
  frameTimesRev = reverse(frameTimes)
  
  trailLenFrame = abs(fix(trailLenTime / timeStep))
  
  ; setup fading trail color sequence
  rgb = fix(findgen(trailLenFrame)/(trailLenFrame-1) * 254)
  
  trailColorsGal  = lonarr(trailLenFrame)
  trailColorsGmem = lonarr(trailLenFrame)
  
  for i=0,trailLenFrame-1 do begin
    trailColorsGal[i]  = rgb[i] + rgb[i]*(256L) + fix(rgb[i]/ggTintFac)*(256L)^2 ; blue tint
    trailColorsGmem[i] = rgb[i] + fix(rgb[i]/ggTintFac)*(256L) + rgb[i]*(256L)^2 ; green tint
  endfor
  
  trailColorsGal = reverse(trailColorsGal)
  trailColorsGmem = reverse(trailColorsGmem)
  
  ; setup high temperature color trail
  trailColorModGal  = lonarr(trailLenFrame)
  trailColorModGMem = lonarr(trailLenFrame)
  
  for i=0,trailLenFrame-1 do begin
    trailColorModGal[i]  = fix(rgb[i]/ggTintFac*0.25) + rgb[i]*(256L) + fix(rgb[i]/ggTintFac)*(256L)^2
    trailColorModGmem[i] = fix(rgb[i]/ggTintFac*0.25) + fix(rgb[i]/ggTintFac)*(256L) + rgb[i]*(256L)^2
  endfor
  
  trailColorModGal  = reverse(trailColorModGal)
  trailColorModGmem = reverse(trailColorModGmem)

  ; spline interpolation of halo properties (single or mean)
  print,'interp halo...'
  interpHaloMass = fltarr(n_elements(haloIDs),nFrames)
  interpHaloRvir = fltarr(n_elements(haloIDs),nFrames)
  interpHaloTvir = fltarr(n_elements(haloIDs),nFrames)
  
  for i=0,n_elements(haloIDs)-1 do begin
    ; reverse halo properties so they progress forward in time
    interpHaloMass[i,*] = reverse(codeMassToLogMsun(spline(mt.times,mt.hMass[*,haloIDs[i]],frameTimesRev)))
    interpHaloRvir[i,*] = reverse(spline(mt.times,mt.hVirRad[*,haloIDs[i]],frameTimesRev))
    interpHaloTvir[i,*] = reverse(spline(mt.times,mt.hVirTemp[*,haloIDs[i]],frameTimesRev))
  endfor
  
  ; spline interpolation of axes[0],axes[1] coordinates and current temperatures
  galPos   = fltarr(2,nPart.gal,nFrames)
  gmemPos  = fltarr(2,nPart.gmem,nFrames)
  galTemp  = fltarr(nPart.gal,nFrames)
  gmemTemp = fltarr(nPart.gmem,nFrames) 

  print,'interp gal...'
  for i=0,nPart.gal-1 do begin
    galPos[0,i,*] = spline(mt.times,at.relPos_gal[*,axes[0],galcatSub.gal[i]],frameTimesRev)
    galPos[1,i,*] = spline(mt.times,at.relPos_gal[*,axes[1],galcatSub.gal[i]],frameTimesRev)
    galTemp[i,*]  = spline(mt.times,at.curTemp_gal[*,galcatSub.gal[i]],frameTimesRev)
  endfor
  
  print,'interp gmem...'
  for i=0,nPart.gmem-1 do begin
    gmemPos[0,i,*] = spline(mt.times,at.relPos_gmem[*,axes[0],galcatSub.gmem[i]],frameTimesRev)
    gmemPos[1,i,*] = spline(mt.times,at.relPos_gmem[*,axes[1],galcatSub.gmem[i]],frameTimesRev)
    gmemTemp[i,*]  = spline(mt.times,at.curTemp_gmem[*,galcatSub.gmem[i]],frameTimesRev)
  endfor
  
  ; reverse positions/temperatures so they progress forward in time
  galPos   = reverse(galPos,3,/overwrite)
  gmemPos  = reverse(gmemPos,3,/overwrite)
  galTemp  = reverse(galTemp,2,/overwrite)
  gmemTemp = reverse(gmemTemp,2,/overwrite)

  print,'rendering frames...'
  for fn=0,nFrames-1 do begin
  ;fn = 100
    ; determine time of this frame and bracketing snapshots
    time = timeStart + timeStep * fn
    
    ; determine earliest frame number to draw trail to
    fnTrail = fn - trailLenFrame + 1 < fn > 0
    numTrail = fn - fnTrail + 1

    print,' ['+string(fn,format='(i3)')+'] z = '+string(time,format='(f5.3)')+' ( '+str(fn)+' - '+str(fnTrail)+' )'
    
    ; plot
    start_PS, sP.plotPath + 'frames/frame_'+str(fn)+'.eps', xs=8, ys=8
    
      !p.thick = 1.0
      
      ; fixed viewport
      xrange = [boxCen[0]-zoomSize/2.0,boxCen[0]+zoomSize/2.0]
      yrange = [boxCen[1]-zoomSize/2.0,boxCen[1]+zoomSize/2.0]
    
      ; start plot
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=5,ystyle=5,$
             xtitle="",ytitle="",title="",position=[0,0,1.0,1.0],xtickname=replicate(' ',10),$
             ytickname=replicate(' ',10)
           
      ; draw rvir, 10% rvir and 1% rvir
      if n_elements(haloIDs) eq 1 then begin
        haloRvir = interpHaloRvir[fn]
        
        tvcircle,haloRvir/10,0,0,cgColor('green'),thick=0.5,/data
        tvcircle,haloRvir/100,0,0,cgColor('green'),thick=0.5,/data
        tvcircle,haloRvir,0,0,cgColor('green'),thick=0.5,/data
      endif else begin
        ; draw min and max for population of halos
        haloRvir_min = min(interpHaloRvir[*,fn])
        haloRvir_max = max(interpHaloRvir[*,fn])
        
        tvcircle,haloRvir_min/10,0,0,cgColor('green'),thick=0.5,/data
        tvcircle,haloRvir_min,0,0,cgColor('green'),thick=0.5,/data
        
        tvcircle,haloRvir_max/10,0,0,cgColor('green'),thick=0.5,/data
        tvcircle,haloRvir_max,0,0,cgColor('green'),thick=0.5,/data
      endelse
      
      ; decide a color modifier value based on current temperature of each particle vs. mean Tvir
      haloTvir_mean = mean(interpHaloTvir[*,fn])

      ; gal: to plot the trail, plots each particle path separately
      w = where(alog10(10.0^galTemp[*,fn]/10.0^haloTvir_mean) ge tempLogFac,count,comp=wc,ncomp=countc)
      if count gt 0 then $
        foreach i,w do for j=fnTrail,fn-1 do $
          plots,galPos[0,i,j:(j+1)],galPos[1,i,j:(j+1)],color=trailColorModGal[j-fnTrail]
      if countc gt 0 then $
        foreach i,wc do for j=fnTrail,fn-1 do $
          plots,galPos[0,i,j:(j+1)],galPos[1,i,j:(j+1)],color=trailColorsGal[j-fnTrail]
                
      ; gmem
      w = where(alog10(10.0^gmemTemp[*,fn]/10.0^haloTvir_mean) ge tempLogFac,count,comp=wc,ncomp=countc)
      if count gt 0 then $
        foreach i,w do for j=fnTrail,fn-1 do $
          plots,gmemPos[0,i,j:(j+1)],gmemPos[1,i,j:(j+1)],color=trailColorModGmem[j-fnTrail]
      if countc gt 0 then $
        foreach i,wc do for j=fnTrail,fn-1 do $
          plots,gmemPos[0,i,j:(j+1)],gmemPos[1,i,j:(j+1)],color=trailColorsGmem[j-fnTrail]
          
      ; scale bar
      len = 100.0 ;ckpc
      cgPlot,[boxCen[0]-zoomSize/2.1,boxCen[0]-zoomSize/2.1+len],[boxCen[1]+zoomSize/2.1,boxCen[1]+zoomSize/2.1],$
        line=0,color=cgColor('dark gray'),/overplot
      cgText,mean([boxCen[0]-zoomSize/2.1,boxCen[0]-zoomSize/2.1+len]),boxCen[1]+zoomSize/2.2,$
        string(len,format='(i3)')+' ckpc',alignment=0.5,charsize=!p.charsize-0.6,color=cgColor('dark gray')
        
      ; write time ticker
      strTime = "z = "+string(time,format='(f5.3)')
      cgText,0.5,0.03,strTime,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
      
      ; write halo mass or mean
      if n_elements(haloIDs) eq 1 then $
        strHalo = "log(M) = "+string(interpHaloMass[fn],format='(f4.1)')
      if n_elements(haloIDs) gt 1 then $
        strHalo = string(min(interpHaloMass[*,fn]),format='(f4.1)') + " < log(M) < " + $
                  string(max(interpHaloMass[*,fn]),format='(f4.1)')
                  
      cgText,0.5,0.06,strHalo,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
      ;strHalo = textoidl("R_{vir} = ")+string(interpHaloRvir[fn],format='(i3)')
      ;cgText,0.5,0.13,strHalo,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
      
    end_PS, pngResize=60, im_options='-negate', /deletePS

  endfor

end
