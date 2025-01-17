; accretionTrajVis.pro
; particle trajectories with time (visualization routines)
; dnelson jul.2013

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
  
    haloIDs_vel = [match.wMatch]
    haloIDs_ga  = [match.matchedInds[match.wMatch]]
    print,'CHECK the above line (subset matchedInds?)'
        
    ; load group catalog to find masses of these matched halos
    gc = loadGroupCat(sP=sPa,/skipIDs)
    haloMasses = codeMassToLogMsun(gc.subgroupMass[haloIDs_vel])
    gc = !NULL
    
    ; choose halo
    haloID_vel = haloIDs_vel[matchID]
    haloID_ga  = haloIDs_ga[matchID]
    haloMass   = haloMasses[matchID]

    ; get indices for the subset of the galaxycat corresponding to the matched halo in each sim
    galcatSub_vel = galCatINDList(sP=sPa,gcIDList=[haloID_vel],/trRep)          
    galcatSub_ga  = galCatINDList(sP=sPg,gcIDList=[haloID_ga])
  endif else begin
    ; make pseudo galcatSub which accesses all members of galaxycat
    galcatSub_vel = lindgen(n_elements(at_vel.relPos[0,0,*]))
    galcatSub_ga  = lindgen(n_elements(at_ga.relPos[0,0,*]))   
  endelse

  ; count
  nTrMC  = 0
  nTrVel = n_elements(at_vel.relPos[0,0,galcatSub_vel])
  nSph   = n_elements(at_ga.relPos[0,0,galcatSub_ga])

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
    velPos[0,i,*] = spline(mt_vel.times,at_vel.relPos[*,axes[0],galcatSub_vel[i]],frameTimesRev)
    velPos[1,i,*] = spline(mt_vel.times,at_vel.relPos[*,axes[1],galcatSub_vel[i]],frameTimesRev)
  endfor
  
  print,'interp sph...'
  for i=0,nSph-1 do begin
    ; spline interpolate to reversed frameTimes (so they progress back in time like the data)
    gaPos[0,i,*] = spline(mt_ga.times,at_ga.relPos[*,axes[0],galcatSub_ga[i]],frameTimesRev)
    gaPos[1,i,*] = spline(mt_ga.times,at_ga.relPos[*,axes[1],galcatSub_ga[i]],frameTimesRev)
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
      cgText,0.05,0.03,"AREPO",alignment=0.5,color=cgColor('dark orchid'),/normal,charsize=!p.charsize-0.4
      cgText,0.95,0.03,"GADGET",alignment=0.5,color=cgColor('dark orchid'),/normal,charsize=!p.charsize-0.4
      
      ; dividing line
      cgPlot,[boxCen[0]-zoomSize/2.0,boxCen[0]-zoomSize/2.0],[boxCen[1]-zoomSize/2.4,boxCen[1]+zoomSize/2.0],$
        line=0,color=cgColor('light gray'),/overplot
      
    end_PS, pngResize=60, im_options='-negate', /deletePS

  endfor

end

; cosmoTrajMovieHalo(): load positions(time) and make interpolated images of tracer/sph positions
;                           from a cosmological box via accretionTraj() for one halo or a mass bin
; note: uses cubic hermite spline interpolation with optionally known derivatives (velocities)

pro cosmoTrajMovieHalo

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  sP = simParams(res=256,run='gadget',redshift=0.0)

  ; choose one: plot a single halo? or a mass bin of halos?
  massBinLog = [] ;[12.0,12.5]
  mtHaloInd  = 0
  
  ; view configuration
  useVels  = 0 ; use velocities as known derivatives
  boxCen   = [0,0]
  zoomSize = 800.0 ;kpc
  axes     = [0,2]
  
  ; movie configuration
  nFrames       = 1501
  redshiftStart = 4.0
  redshiftEnd   = 0.0
  trailLenTime  = 0.06 ;Gyr
  
  ; color config
  ggTintFac    = 1.5   ; color skew for gal,gmem,hot
  tempLogFac   = -0.05 ; within ~90% of Tvir or above
  
  ; smooth halocenter tracks?
  smoothKer = 5 ; 0 to disable
  wrapTol   = 100.0 ;kpc  
  
  ; get gas/tracer positions with time
  at  = accretionTraj(sP=sP)
  mt  = mergerTreeSubset(sP=sP)

  ; convert velocities to kpc/Gyr
  if useVels then at.relVel *= (units.s_in_Gyr / units.kpc_in_km)
  
  haloIDs = []

  ; load group catalog to find masses of tracked halos
  gc = loadGroupCat(sP=sP,/skipIDs)
  haloMasses = codeMassToLogMsun(gc.subgroupMass[mt.galcatIDList])
  gc = !NULL
  
  ; choose halo(s)
  if n_elements(massBinLog) gt 0 then begin
    haloIDs = where(haloMasses ge massBinLog[0] and haloMasses le massBinLog[1],count)
    if count eq 0 then message,'Error: No halos in mass bin.'
    print,'Found ['+str(count)+'] halos in mass bin.'
    
    haloMasses = halomasses[mtHaloInd]
  endif else begin
    haloIDs    = [mtHaloInd]
    haloMasses = haloMasses[mtHaloInd]
  endelse
  
  ; count
  nPart = n_elements(at.relPos[0,0,*])

  ; convert scale factor to redshift
  mt.times = 1/mt.times-1
  
  if redshiftStart lt min(mt.times) then message,'Error: Requested redshiftStart beyond snapshot range.'
  if redshiftEnd gt max(mt.times) then message,'Error: Requested redshiftEnd beyond snapshot range.'
  if redshiftEnd ne sP.redshift then message,'Are you sure?'
  
  ; stepping in redshift
  redshiftStep   = (redshiftEnd - redshiftStart) / (nFrames-1)
  frameRedshifts = redshiftStep * findgen(nFrames) + redshiftStart
  
  ; stepping in Gyr
  timeStart  = redshiftToAgeFlat(redshiftStart) + 1e-4
  timeEnd    = redshiftToAgeFlat(redshiftEnd) - 1e-4
  timeStep   = (timeEnd - timeStart) / (nFrames-1)
  
  frameTimes = timeStep * findgen(nFrames) + timeStart
  frameTimesRev = reverse(frameTimes)
  
  trailLenFrame = abs(fix(trailLenTime / timeStep))

  ; convert mt.times to Gyr
  mt.times = redshiftToAgeFlat(mt.times)

  ; setup fading trail color sequence
  rgb = fix(findgen(trailLenFrame)/(trailLenFrame-1) * 254)
  
  trailColorsGal  = lonarr(trailLenFrame)
  
  for i=0,trailLenFrame-1 do begin
    ;trailColors[i]  = rgb[i] + rgb[i]*(256L) + fix(rgb[i]/ggTintFac)*(256L)^2 ; blue tint
    trailColors[i]  = rgb[i] + rgb[i]*(256L) + rgb[i]*(256L)^2 ; blue tint
  endfor
  
  trailColors = reverse(trailColors)
  
  ; setup high temperature color trail
  trailColorMod  = lonarr(trailLenFrame)
  
  for i=0,trailLenFrame-1 do $
    trailColorMod[i] = fix(rgb[i]/ggTintFac*0.25) + rgb[i]*(256L) + fix(rgb[i]/ggTintFac)*(256L)^2
  
  trailColorMod = reverse(trailColorMod)

  ; spline interpolation of halo properties (single or mean)
  print,'interp halo...'
  interpHaloMass = fltarr(n_elements(haloIDs),nFrames)
  interpHaloRvir = fltarr(n_elements(haloIDs),nFrames)
  interpHaloTvir = fltarr(n_elements(haloIDs),nFrames)
  
  for i=0,n_elements(haloIDs)-1 do begin
    ; reverse halo properties so they progress forward in time
    interpHaloMass[i,*] = reverse(codeMassToLogMsun(hermite(mt.times,mt.hMass[*,mtHaloInd[i]],frameTimesRev)))
    interpHaloRvir[i,*] = reverse(hermite(mt.times,mt.hVirRad[*,mtHaloInd[i]],frameTimesRev))
    interpHaloTvir[i,*] = reverse(hermite(mt.times,mt.hVirTemp[*,mtHaloInd[i]],frameTimesRev))
  endfor

  ; spline interpolation of axes[0],axes[1] coordinates and current temperatures
  pos  = fltarr(2,nPart,nFrames)
  temp = fltarr(nPart,nFrames)

  print,'interp...'
  
  if useVels then begin
    ; interpolate using velocities are known derivatives
    for i=0L,nPart-1 do begin
      pos[0,i,*] = hermite(mt.times,at.relPos[*,axes[0],i],frameTimesRev,$
                              fderiv=at.relVel[*,axes[0],i])
      pos[1,i,*] = hermite(mt.times,at.relPos[*,axes[1],i],frameTimesRev,$
                              fderiv=at.relVel[*,axes[1],i])
      temp[i,*]  = hermite(mt.times,at.curTemp[*,i],frameTimesRev)
    endfor
  endif else begin
    ; interpolate without using known derivatives
    for i=0L,nPart-1 do begin
      pos[0,i,*] = hermite(mt.times,at.relPos[*,axes[0],i],frameTimesRev)
      pos[1,i,*] = hermite(mt.times,at.relPos[*,axes[1],i],frameTimesRev)
      temp[i,*]  = hermite(mt.times,at.curTemp[*,i],frameTimesRev)
    endfor
  endelse
  
  ; reverse positions/temperatures so they progress forward in time
  pos  = reverse(pos,3,/overwrite)
  temp = reverse(temp,2,/overwrite)

  print,'rendering frames...'
  ;frames = [0,300,600,900,1200]
  ;foreach fn,frames do begin
  for fn=0,nFrames-1 do begin
  ;fn = 400
    ; determine time of this frame and bracketing snapshots
    deltatime = timeStep * fn ;Gyr
    redshift  = frameRedshifts[fn]
    
    ; determine earliest frame number to draw trail to
    fnTrail = fn - trailLenFrame + 1 < fn > 0
    numTrail = fn - fnTrail + 1

    print,' ['+string(fn,format='(i4)')+'] z = '+string(redshift,format='(f5.3)')+' ( '+$
      string(fn,format='(i4)')+' - '+string(fnTrail,format='(i4)')+' )'
    
    ; plot
    start_PS, sP.plotPath + 'frames3/frame_'+str(fn)+'.eps', xs=8, ys=8
    
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
      haloTvir_mean = interpHaloTvir[fn]
      if n_elements(haloIDs) gt 1 then haloTvir_mean = mean(interpHaloTvir[*,fn])

      ; to plot the trail, plots each particle path separately
      w = where(alog10(10.0^temp[*,fn]/10.0^haloTvir_mean) ge tempLogFac,count,comp=wc,ncomp=countc)
      
      if count gt 0 then $
        foreach i,w do for j=fnTrail,fn-1 do $
          plots,pos[0,i,j:(j+1)],pos[1,i,j:(j+1)],color=trailColorMod[j-fnTrail]
      if countc gt 0 then $
        foreach i,wc do for j=fnTrail,fn-1 do $
          plots,pos[0,i,j:(j+1)],pos[1,i,j:(j+1)],color=trailColors[j-fnTrail]
                
      ; scale bar
      len = 100.0 ;ckpc
      cgPlot,[boxCen[0]-zoomSize/2.1,boxCen[0]-zoomSize/2.1+len],[boxCen[1]+zoomSize/2.1,boxCen[1]+zoomSize/2.1],$
        line=0,color=cgColor('dark gray'),/overplot
      cgText,mean([boxCen[0]-zoomSize/2.1,boxCen[0]-zoomSize/2.1+len]),boxCen[1]+zoomSize/2.2,$
        string(len,format='(i3)')+' ckpc',alignment=0.5,charsize=!p.charsize-0.6,color=cgColor('dark gray')
        
      ; write time ticker
      strTime = textoidl("\Delta")+"t = "+string(deltatime,format='(f4.1)')+" Gyr"
      cgText,0.5,0.06,strTime,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
      strTime = "z = "+string(redshift,format='(f5.3)')
      cgText,0.5,0.03,strTime,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
      
      ; write halo mass or mean
      if n_elements(haloIDs) eq 1 then $
        strHalo = "log(M) = "+string(interpHaloMass[fn],format='(f4.1)')
      if n_elements(haloIDs) gt 1 then $
        strHalo = string(min(interpHaloMass[*,fn]),format='(f4.1)') + " < log(M) < " + $
                  string(max(interpHaloMass[*,fn]),format='(f4.1)')
                  
      cgText,0.5,0.09,strHalo,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
      
    end_PS, pngResize=60, im_options='-negate', /deletePS

  endfor

end

; cosmoTrajMovieGasDM(): comparison between gas and DM for a single halo
; note: uses hermite interpolation -without- velocities and uses accretionTrajSingle()
;
; split = [nProcs,curProc] to split work into disjoint sets for parallelization (no good with hotMask)

pro cosmoTrajMovieGasDM, split=split

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  sP = simParams(res=256,run='gadget',redshift=0.0)
  plotPath = sP.plotPath + 'frames.256.h0/'
  
  ; view configuration
  boxCen   = [0,0]
  zoomSize = 900.0   ; kpc
  axes     = [2,0,1] ; rotation about first axis
  numRots  = 0.0     ; total number of rotations throughout duration of movie (0=disable)
  
  ; movie configuration
  nFrames       = 1501 ; sets time spacing, etc.
  runFrames     = pSplit(indgen(nFrames),split=split)
  
  nFramesEndOrbit  = 300 ; do one complete orbit at fixed time at the end (0=disable)
  midOrbitRedshift = 2 ; do one complete orbit at some intermediate redshift (99=disable)
  
  redshiftStart = 4.0
  redshiftEnd   = 0.0
  trailLenTime  = 0.04 ;Gyr
  pssize = [8,4] ; inches
  
  ; color config
  ggTintFac    = 2.0   ; color skew for gal,gmem,hot
  tempLogFac   = -0.05 ; within ~90% of Tvir or above
  
  ; load mergerTreeSubset to find masses of tracked halos and for time sequence
  mt = mergerTreeSubset(sP=sP)

  hInd = 0 ; 0,2,15,74,200,430,640 (13.1 12.5 12.0 11.5 11.0 10.5 10.0 for 256z=0) or select e.g. by mass
           ; also 1,3 are ~12.5 massive examples

  ; load trajectories of chosen halo
  at = accretionTrajSingle(sP=sP,hInd=hInd)
  nPart = {gas : n_elements(at.relPos_gas[0,0,*]), dm : n_elements(at.relPos_dm[0,0,*])}

  mt.times  = 1/mt.times-1 ; convert scale factor to redshift
  
  if redshiftStart lt min(mt.times) then message,'Error: Requested redshiftStart beyond snapshot range.'
  if redshiftEnd gt max(mt.times) then message,'Error: Requested redshiftEnd beyond snapshot range.'
  if redshiftEnd ne sP.redshift then message,'Are you sure?'
  
  ; stepping in redshift
  redshiftStep   = (redshiftEnd - redshiftStart) / (nFrames-1)
  frameRedshifts = redshiftStep * findgen(nFrames) + redshiftStart
  
  ; stepping in Gyr
  timeStart  = redshiftToAgeFlat(redshiftStart) + 1e-4
  timeEnd    = redshiftToAgeFlat(redshiftEnd) - 1e-4
  timeStep   = (timeEnd - timeStart) / (nFrames-1)
  
  frameTimes    = timeStep * findgen(nFrames) + timeStart
  frameTimesRev = reverse(frameTimes)
  trailLenFrame = abs(fix(trailLenTime / timeStep))
  
  mt.times = redshiftToAgeFlat(mt.times) ; convert mt.times to Gyr

  ; setup fading trail color sequence (blue tint) and high temperature color trail (red tint)
  rgb = fix(findgen(trailLenFrame)/(trailLenFrame-1) * 254)
  
  trailColors   = lonarr(trailLenFrame)
  trailColorMod = lonarr(trailLenFrame)
  
  for i=0,trailLenFrame-1 do trailColors[i]   = rgb[i] + rgb[i]*(256L) + fix(rgb[i]/ggTintFac*1.5)*(256L)^2
  for i=0,trailLenFrame-1 do trailColorMod[i] = fix(rgb[i]/ggTintFac) + rgb[i]*(256L) + rgb[i]*(256L)^2
  
  trailColors   = reverse(trailColors)
  trailColorMod = reverse(trailColorMod)

  ; rotations
  frameRotations = findgen(nFrames)/(nFrames-1) * 360.0 * numRots ; deg

  ; make mask for gas particles that have gone "hot"
  hotMask = bytarr(nPart.gas)

  ; spline interpolation of halo properties (reverse so they progress forward in time)
  print,'interp halo...'
  interpHaloMass = reverse(codeMassToLogMsun(hermite(mt.times,mt.hMass[*,hInd],frameTimesRev)))
  interpHaloRvir = reverse(hermite(mt.times,mt.hVirRad[*,hInd],frameTimesRev))
  interpHaloTvir = reverse(hermite(mt.times,mt.hVirTemp[*,hInd],frameTimesRev))

  ; spline interpolation of axes[0],axes[1] coordinates and current temperatures
  gasPos  = fltarr(3,nPart.gas,nFrames)
  dmPos   = fltarr(3,nPart.dm,nFrames)
  gasTemp = fltarr(nPart.gas,nFrames)

  ; interpolate without using known derivatives
  print,'interp gas...'
  for i=0,nPart.gas-1 do begin
    gasPos[0,i,*] = hermite(mt.times,at.relPos_gas[*,0,i],frameTimesRev)
    gasPos[1,i,*] = hermite(mt.times,at.relPos_gas[*,1,i],frameTimesRev)
    gasPos[2,i,*] = hermite(mt.times,at.relPos_gas[*,2,i],frameTimesRev)
    gasTemp[i,*]  = hermite(mt.times,at.curTemp_gas[*,i],frameTimesRev)
  endfor
  
  print,'interp dm...'
  for i=0,nPart.dm-1 do begin
    dmPos[0,i,*] = hermite(mt.times,at.relPos_dm[*,0,i],frameTimesRev)
    dmPos[1,i,*] = hermite(mt.times,at.relPos_dm[*,1,i],frameTimesRev)
    dmPos[2,i,*] = hermite(mt.times,at.relPos_dm[*,2,i],frameTimesRev)
  endfor
  
  ; reverse positions/temperatures so they progress forward in time
  gasPos  = reverse(gasPos,3,/overwrite)
  dmPos   = reverse(dmPos,3,/overwrite)
  gasTemp = reverse(gasTemp,2,/overwrite)

  ; find frame closest to z=2 for orbit in the middle of the movie
  midOrbitFrame = value_locate(frameRedshifts,midOrbitRedshift) + 1
  fnOff = 0
  
  print,'rendering frames...'
  foreach fn,runFrames do begin
  ;fn = 200
    ; determine time of this frame and bracketing snapshots
    deltatime = timeStep * fn ;Gyr
    redshift  = frameRedshifts[fn]
    
    ; do an orbit in the middle near z=2
    if fn eq midOrbitFrame then begin
      endingRot = curRot
      
      for fnOrbit=0,nFramesEndOrbit-1 do begin        
        ; determine current rotation and create v1,v2 (screen plane projections)
        curRot = endingRot + (360.0/nFramesEndOrbit) * fnOrbit
        
        v2_gas = reform(gasPos[axes[0],*,fnTrail:fn]) ; keep axes[0] always constant
        v1_gas = reform(gasPos[axes[1],*,fnTrail:fn]) * cos(curRot*!dtor) + $
                 reform(gasPos[axes[2],*,fnTrail:fn]) * sin(curRot*!dtor)
                        
        v2_dm  = reform(dmPos[axes[0],*,fnTrail:fn])
        v1_dm  = reform(dmPos[axes[1],*,fnTrail:fn]) * cos(curRot*!dtor) + $
                 reform(dmPos[axes[2],*,fnTrail:fn]) * sin(curRot*!dtor)
    
        ; clip positions at box dimensions (left,right)
        w = where(v1_dm lt boxCen[0]-zoomSize/2.0,count)
        if count gt 0 then v1_dm[w] = !values.f_nan
        
        w = where(v1_gas gt boxCen[0]+zoomSize/2.0,count)
        if count gt 0 then v1_gas[w] = !values.f_nan
        w = !NULL
    
        print,' ['+string(fn+fnOff,format='(i4)')+'] orbit ['+$
          str(fnOrbit)+' of '+str(nFramesEndOrbit)+']'
        
        ; plot
        start_PS, plotPath + 'frame_'+str(fn+fnOff)+'.eps', xs=pssize[0], ys=pssize[1]
        
          !p.thick = 1.0
    
          ; fixed viewport
          xrange = [boxCen[0]-zoomSize/2.0,boxCen[0]+zoomSize/2.0]
          yrange = [boxCen[1]-zoomSize/2.0,boxCen[1]+zoomSize/2.0] ; top+bottom halves
        
          ; left: gas
          cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=5,ystyle=5,$
                 xtitle="",ytitle="",title="",position=[0,0,0.5,1],xtickname=replicate(' ',10),$
                 ytickname=replicate(' ',10)
               
          ; left: draw rvir and rvir/10
          tvcircle,interpHaloRvir[fn]/10,0,0,cgColor('green'),thick=0.5,/data
          tvcircle,interpHaloRvir[fn],0,0,cgColor('green'),thick=0.5,/data
            
          ; plot gas with constant hotMask
          for i=0,nPart.gas-1 do begin
            ; to plot the trail, plots each particle path separately
            if hotMask[i] eq 1B then begin
              for j=0,numTrail-2 do $
                plots,v1_gas[i,j:(j+1)],v2_gas[i,j:(j+1)],color=trailColorMod[j]
            endif else begin
              for j=0,numTrail-2 do $
                plots,v1_gas[i,j:(j+1)],v2_gas[i,j:(j+1)],color=trailColors[j]
            endelse
          endfor
            
          ; scale bar
          len = 100.0 ;ckpc
          cgPlot,[boxCen[0]-zoomSize/2.2,boxCen[0]-zoomSize/2.2+len],[boxCen[1]+zoomSize/2.2,boxCen[1]+zoomSize/2.2],$
            line=0,color=cgColor('dark gray'),/overplot
          cgText,mean([boxCen[0]-zoomSize/2.2,boxCen[0]-zoomSize/2.2+len]),boxCen[1]+zoomSize/2.4,$
            string(len,format='(i3)')+' ckpc',alignment=0.5,charsize=!p.charsize-0.6,color=cgColor('dark gray')
            
          ; right: dm (tile x2 in the left-right direction)
          cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=5,ystyle=5,$
                 xtitle="",ytitle="",title="",position=[0.5,0,1,1],xtickname=replicate(' ',10),$
                 ytickname=replicate(' ',10),/noerase
            
          tvcircle,interpHaloRvir[fn]/10,0,0,cgColor('green'),thick=0.5,/data
          tvcircle,interpHaloRvir[fn],0,0,cgColor('green'),thick=0.5,/data
    
          for i=0,nPart.dm-1 do begin
            ; to plot the trail, plots each particle path separately
            for j=0,numTrail-2 do $
              plots,v1_dm[i,j:(j+1)],v2_dm[i,j:(j+1)],color=trailColors[j]
          endfor
        
          ; draw time ticker (constant)
          strTime = textoidl("\Delta")+"t = "+string(deltatime,format='(f4.1)')+" Gyr"
          cgText,0.5,0.06,strTime,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
          strTime = "z = "+string(redshift,format='(f5.3)')
          cgText,0.5,0.02,strTime,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
          
          ; halo mass (constant)
          strHalo = "log(M) = "+string(interpHaloMass[fn],format='(f4.1)')
          cgText,0.5,0.10,strHalo,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
          
          ; particle type names
          cgText,0.03,0.03,"Gas",alignment=0.5,color=cgColor('dark orchid'),/normal,charsize=!p.charsize-0.4
          cgText,0.97,0.03,"DM",alignment=0.5,color=cgColor('dark orchid'),/normal,charsize=!p.charsize-0.4
          
          ; dividing line
          cgPlot,[boxCen[0]-zoomSize/2.0,boxCen[0]-zoomSize/2.0],[boxCen[1]-zoomSize/2.4,boxCen[1]+zoomSize/2.0],$
            line=0,color=cgColor('light gray'),/overplot
    
        end_PS, pngResize=60, im_options='-negate', /deletePS
        
        ; increment counter for frame number
        fnOff += 1
      endfor
    
    endif
    
    ; determine earliest frame number to draw trail to
    fnTrail = fn - trailLenFrame + 1 < fn > 0
    numTrail = fn - fnTrail + 1
    
    ; determine current rotation and create v1,v2 (screen plane projections)
    curRot = frameRotations[fn]
    
    v2_gas = reform(gasPos[axes[0],*,fnTrail:fn]) ; keep axes[0] always constant
    v1_gas = reform(gasPos[axes[1],*,fnTrail:fn]) * cos(curRot*!dtor) + $
             reform(gasPos[axes[2],*,fnTrail:fn]) * sin(curRot*!dtor)
                    
    v2_dm  = reform(dmPos[axes[0],*,fnTrail:fn])
    v1_dm  = reform(dmPos[axes[1],*,fnTrail:fn]) * cos(curRot*!dtor) + $
             reform(dmPos[axes[2],*,fnTrail:fn]) * sin(curRot*!dtor)

    ; clip positions at box dimensions (left,right)
    w = where(v1_dm lt boxCen[0]-zoomSize/2.0,count)
    if count gt 0 then v1_dm[w] = !values.f_nan
    
    w = where(v1_gas gt boxCen[0]+zoomSize/2.0,count)
    if count gt 0 then v1_gas[w] = !values.f_nan
    w = !NULL

    print,' ['+string(fn+fnOff,format='(i4)')+'] z = '+string(redshift,format='(f5.3)')+' ( '+$
      string(fn,format='(i4)')+' - '+string(fnTrail,format='(i4)')+' )'
    
    ; plot
    start_PS, plotPath + 'frame_'+str(fn+fnOff)+'.eps', xs=pssize[0], ys=pssize[1]
    
      !p.thick = 1.0

      ; fixed viewport
      xrange = [boxCen[0]-zoomSize/2.0,boxCen[0]+zoomSize/2.0]
      yrange = [boxCen[1]-zoomSize/2.0,boxCen[1]+zoomSize/2.0] ; top+bottom halves
    
      ; left: gas
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=5,ystyle=5,$
             xtitle="",ytitle="",title="",position=[0,0,0.5,1],xtickname=replicate(' ',10),$
             ytickname=replicate(' ',10)
           
      ; left: draw rvir and rvir/10
      tvcircle,interpHaloRvir[fn]/10,0,0,cgColor('green'),thick=0.5,/data
      tvcircle,interpHaloRvir[fn],0,0,cgColor('green'),thick=0.5,/data
        
      ; gas: update color modifier mask based on current temperature of each particle vs. Tvir
      w = where(alog10(10.0^gasTemp[*,fn]/10.0^interpHaloTvir[fn]) ge tempLogFac,count)
      if count gt 0 then hotMask[w] = 1B
      
      for i=0,nPart.gas-1 do begin
        ; to plot the trail, plots each particle path separately
        if hotMask[i] eq 1B then begin
          for j=0,numTrail-2 do $
            plots,v1_gas[i,j:(j+1)],v2_gas[i,j:(j+1)],color=trailColorMod[j]
        endif else begin
          for j=0,numTrail-2 do $
            plots,v1_gas[i,j:(j+1)],v2_gas[i,j:(j+1)],color=trailColors[j]
        endelse
      endfor
        
      ; scale bar
      len = 100.0 ;ckpc
      cgPlot,[boxCen[0]-zoomSize/2.2,boxCen[0]-zoomSize/2.2+len],[boxCen[1]+zoomSize/2.2,boxCen[1]+zoomSize/2.2],$
        line=0,color=cgColor('dark gray'),/overplot
      cgText,mean([boxCen[0]-zoomSize/2.2,boxCen[0]-zoomSize/2.2+len]),boxCen[1]+zoomSize/2.4,$
        string(len,format='(i3)')+' ckpc',alignment=0.5,charsize=!p.charsize-0.6,color=cgColor('dark gray')
        
      ; right: dm (tile x2 in the left-right direction)
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=5,ystyle=5,$
             xtitle="",ytitle="",title="",position=[0.5,0,1,1],xtickname=replicate(' ',10),$
             ytickname=replicate(' ',10),/noerase
        
      tvcircle,interpHaloRvir[fn]/10,0,0,cgColor('green'),thick=0.5,/data
      tvcircle,interpHaloRvir[fn],0,0,cgColor('green'),thick=0.5,/data

      for i=0,nPart.dm-1 do begin
        ; to plot the trail, plots each particle path separately
        for j=0,numTrail-2 do $
          plots,v1_dm[i,j:(j+1)],v2_dm[i,j:(j+1)],color=trailColors[j]
      endfor
    
      ; draw time ticker
      strTime = textoidl("\Delta")+"t = "+string(deltatime,format='(f4.1)')+" Gyr"
      cgText,0.5,0.06,strTime,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
      strTime = "z = "+string(redshift,format='(f5.3)')
      cgText,0.5,0.02,strTime,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
      
      ; halo mass
      strHalo = "log(M) = "+string(interpHaloMass[fn],format='(f4.1)')
      cgText,0.5,0.10,strHalo,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
      
      ; particle type names
      cgText,0.03,0.03,"Gas",alignment=0.5,color=cgColor('dark orchid'),/normal,charsize=!p.charsize-0.4
      cgText,0.97,0.03,"DM",alignment=0.5,color=cgColor('dark orchid'),/normal,charsize=!p.charsize-0.4
      
      ; dividing line
      cgPlot,[boxCen[0]-zoomSize/2.0,boxCen[0]-zoomSize/2.0],[boxCen[1]-zoomSize/2.4,boxCen[1]+zoomSize/2.0],$
        line=0,color=cgColor('light gray'),/overplot

    end_PS, pngResize=60, im_options='-negate', /deletePS

  endforeach
  
  ; do ending orbit
  endingRot = curRot
  
  for fnOrbit=0,nFramesEndOrbit-1 do begin    
    ; determine current rotation and create v1,v2 (screen plane projections)
    curRot = endingRot + (360.0/nFramesEndOrbit) * fnOrbit
    
    v2_gas = reform(gasPos[axes[0],*,fnTrail:fn]) ; keep axes[0] always constant
    v1_gas = reform(gasPos[axes[1],*,fnTrail:fn]) * cos(curRot*!dtor) + $
             reform(gasPos[axes[2],*,fnTrail:fn]) * sin(curRot*!dtor)
                    
    v2_dm  = reform(dmPos[axes[0],*,fnTrail:fn])
    v1_dm  = reform(dmPos[axes[1],*,fnTrail:fn]) * cos(curRot*!dtor) + $
             reform(dmPos[axes[2],*,fnTrail:fn]) * sin(curRot*!dtor)

    ; clip positions at box dimensions (left,right)
    w = where(v1_dm lt boxCen[0]-zoomSize/2.0,count)
    if count gt 0 then v1_dm[w] = !values.f_nan
    
    w = where(v1_gas gt boxCen[0]+zoomSize/2.0,count)
    if count gt 0 then v1_gas[w] = !values.f_nan
    w = !NULL

    print,' ['+string(fn+fnOff+fnOrbit+1,format='(i4)')+'] ending orbit ['+$
      str(fnOrbit)+' of '+str(nFramesEndOrbit)+']'
    
    ; plot
    start_PS, plotPath + 'frame_'+str(fn+fnOff+fnOrbit+1)+'.eps', xs=pssize[0], ys=pssize[1]
    
      !p.thick = 1.0

      ; fixed viewport
      xrange = [boxCen[0]-zoomSize/2.0,boxCen[0]+zoomSize/2.0]
      yrange = [boxCen[1]-zoomSize/2.0,boxCen[1]+zoomSize/2.0] ; top+bottom halves
    
      ; left: gas
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=5,ystyle=5,$
             xtitle="",ytitle="",title="",position=[0,0,0.5,1],xtickname=replicate(' ',10),$
             ytickname=replicate(' ',10)
           
      ; left: draw rvir and rvir/10
      tvcircle,interpHaloRvir[fn]/10,0,0,cgColor('green'),thick=0.5,/data
      tvcircle,interpHaloRvir[fn],0,0,cgColor('green'),thick=0.5,/data
        
      ; plot gas with constant hotMask
      for i=0,nPart.gas-1 do begin
        ; to plot the trail, plots each particle path separately
        if hotMask[i] eq 1B then begin
          for j=0,numTrail-2 do $
            plots,v1_gas[i,j:(j+1)],v2_gas[i,j:(j+1)],color=trailColorMod[j]
        endif else begin
          for j=0,numTrail-2 do $
            plots,v1_gas[i,j:(j+1)],v2_gas[i,j:(j+1)],color=trailColors[j]
        endelse
      endfor
        
      ; scale bar
      len = 100.0 ;ckpc
      cgPlot,[boxCen[0]-zoomSize/2.2,boxCen[0]-zoomSize/2.2+len],[boxCen[1]+zoomSize/2.2,boxCen[1]+zoomSize/2.2],$
        line=0,color=cgColor('dark gray'),/overplot
      cgText,mean([boxCen[0]-zoomSize/2.2,boxCen[0]-zoomSize/2.2+len]),boxCen[1]+zoomSize/2.4,$
        string(len,format='(i3)')+' ckpc',alignment=0.5,charsize=!p.charsize-0.6,color=cgColor('dark gray')
        
      ; right: dm (tile x2 in the left-right direction)
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=5,ystyle=5,$
             xtitle="",ytitle="",title="",position=[0.5,0,1,1],xtickname=replicate(' ',10),$
             ytickname=replicate(' ',10),/noerase
        
      tvcircle,interpHaloRvir[fn]/10,0,0,cgColor('green'),thick=0.5,/data
      tvcircle,interpHaloRvir[fn],0,0,cgColor('green'),thick=0.5,/data

      for i=0,nPart.dm-1 do begin
        ; to plot the trail, plots each particle path separately
        for j=0,numTrail-2 do $
          plots,v1_dm[i,j:(j+1)],v2_dm[i,j:(j+1)],color=trailColors[j]
      endfor
    
      ; draw time ticker (constant)
      strTime = textoidl("\Delta")+"t = "+string(deltatime,format='(f4.1)')+" Gyr"
      cgText,0.5,0.06,strTime,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
      strTime = "z = "+string(redshift,format='(f5.3)')
      cgText,0.5,0.02,strTime,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
      
      ; halo mass (constant)
      strHalo = "log(M) = "+string(interpHaloMass[fn],format='(f4.1)')
      cgText,0.5,0.10,strHalo,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
      
      ; particle type names
      cgText,0.03,0.03,"Gas",alignment=0.5,color=cgColor('dark orchid'),/normal,charsize=!p.charsize-0.4
      cgText,0.97,0.03,"DM",alignment=0.5,color=cgColor('dark orchid'),/normal,charsize=!p.charsize-0.4
      
      ; dividing line
      cgPlot,[boxCen[0]-zoomSize/2.0,boxCen[0]-zoomSize/2.0],[boxCen[1]-zoomSize/2.4,boxCen[1]+zoomSize/2.0],$
        line=0,color=cgColor('light gray'),/overplot

    end_PS, pngResize=60, im_options='-negate', /deletePS
  endfor

end

; cosmoTrajectoryGasDM(): comparison between gas and DM trajectories for a single halo
; note: uses hermite interpolation -without- velocities and uses accretionTrajSingle()
; (plots a single image, or a rotation movie, of the full trajectory streamline)

pro cosmoTrajectoryGasDM

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  sP = simParams(res=512,run='gadget',redshift=2.0)
  plotPath = sP.plotPath + ''
  
  ; view configuration
  boxCen   = [0,0]
  zoomSize = 3.5     ; times r_vir
  axes     = [2,0,1] ; rotation about first axis
  nFrames  = 600     ; number of interpolation points along each trajectory
  nFramesOrbit = 60  ; do one complete orbit at fixed time at the end (0=disable)
  
  virRadCut = 0.15   ; only plot particles ending within this sphere  
  maxTempCut = 5.5   ; only plot particles below this temp threshold for all time
  
  redshiftStart = 4.0
  redshiftEnd   = 0.0
  pssize = [6,6] ; inches
  
  haloID = 304 ;z2.304 z2.301 z2.130 z2.64
  hInd = getMatchedIDs(simParams=sP,haloID=haloID)
  
  ; load mergerTreeSubset to find masses of tracked halos and for time sequence
  mt = mergerTreeSubset(sP=sP)

  ; load trajectories of chosen halo
  at = accretionTrajSingle(sP=sP,hInd=hInd)

  mt.times  = 1/mt.times-1 ; convert scale factor to redshift
  
  if redshiftStart lt min(mt.times) then message,'Error: Requested redshiftStart beyond snapshot range.'
  if redshiftEnd gt max(mt.times) then message,'Error: Requested redshiftEnd beyond snapshot range.'
  if redshiftEnd ne sP.redshift then message,'Are you sure?'
  
  ; stepping in redshift
  redshiftStep   = (redshiftEnd - redshiftStart) / (nFrames-1)
  frameRedshifts = redshiftStep * findgen(nFrames) + redshiftStart
  
  ; stepping in Gyr
  timeStart  = redshiftToAgeFlat(redshiftStart) + 1e-4
  timeEnd    = redshiftToAgeFlat(redshiftEnd) - 1e-4
  timeStep   = (timeEnd - timeStart) / (nFrames-1)
  
  frameTimes    = timeStep * findgen(nFrames) + timeStart
  frameTimesRev = reverse(frameTimes)
  
  mt.times = redshiftToAgeFlat(mt.times) ; convert mt.times to Gyr

  ; make a temperature subselection for all time
  

  ; make a spatial subselection at the ending time
  print,'interp gas...'
  gasRad = sqrt(at.relPos_gas[0,0,*]*at.relPos_gas[0,0,*] + $
                at.relPos_gas[0,1,*]*at.relPos_gas[0,1,*] + $
                at.relPos_gas[0,2,*]*at.relPos_gas[0,2,*])
                
  w = where(gasRad/mt.hVirRad[0,hInd] lt virRadCut,count_gas)
  gasRad = !NULL

  ; spline interpolation of all axes coordinates and current temperatures for this subset
  gasPos  = fltarr(3,count_gas,nFrames)
  gasTemp = fltarr(count_gas,nFrames)

  ; interpolate without using known derivatives
  for i=0,count_gas-1 do begin
    gasPos[0,i,*] = hermite(mt.times,at.relPos_gas[*,0,w[i]],frameTimesRev)
    gasPos[1,i,*] = hermite(mt.times,at.relPos_gas[*,1,w[i]],frameTimesRev)
    gasPos[2,i,*] = hermite(mt.times,at.relPos_gas[*,2,w[i]],frameTimesRev)
    gasTemp[i,*]  = hermite(mt.times,at.curTemp_gas[*,w[i]],frameTimesRev)
  endfor
  
  print,'inside cut: gas '+string(100*float(count_gas)/n_elements(at.relPos_gas[0,0,*]),format='(f4.1)')+'%'
  
  w = !NULL
  at = !NULL
  
  ; reverse positions/temperatures so they progress forward in time
  gasPos  = reverse(gasPos,3,/overwrite)
  gasTemp = reverse(gasTemp,2,/overwrite)

  ; calculate viewbox
  zoomSize = zoomSize * mt.hVirRad[0,hInd]

  ; colormap
  tempMinMax = [4.0,6.0]
  gasCM = (gasTemp-tempMinMax[0])*205.0 / (tempMinMax[1]-tempMinMax[0]) ;0-205
  gasCM = fix(gasCM + 50.0) > 0 < 255 ;50-255
  
  ; start plot
  start_PS, plotPath + 'frame_0.eps', xs=pssize[0], ys=pssize[1]
  
    !p.thick = 0.6
    
    ; clip positions at box dimensions (left,right)
    v1 = reform(gasPos[axes[1],*,*])
    v2 = reform(gasPos[axes[0],*,*])
    
    ;w = where(v1 gt boxCen[0]+zoomSize/2.0,count)
    ;if count gt 0 then v1[w] = !values.f_nan
    ;w = !NULL
    
    ; color table
    loadColorTable,'helix',/reverse
    TVLCT, rr, gg, bb, /GET

    gasCM = reform(gasCM,[n_elements(gasCM)])
    gasCM = getColor24([[rr[gasCM]], [gg[gasCM]], [bb[gasCM]]])
    gasCM = reform(gasCM,[count_gas,nFrames])
    
    ; fixed viewport
    xrange = [boxCen[0]-zoomSize/2.0,boxCen[0]+zoomSize/2.0]
    yrange = [boxCen[1]-zoomSize/2.0,boxCen[1]+zoomSize/2.0] ; top+bottom halves
  
    ; gas
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=5,ystyle=5,$
           xtitle="",ytitle="",title="",position=[0,0,1,1],xtickname=replicate(' ',10),$
           ytickname=replicate(' ',10)
         
    ; draw rvir and rvir/10
    tvcircle,mt.hVirRad[0,hInd]/10,0,0,cgColor('black'),thick=0.8,/data
    tvcircle,mt.hVirRad[0,hInd],0,0,cgColor('black'),thick=0.8,/data  
  
    ; draw gas trajectories
    print,'rendering gas trajectories...'
    
    for i=0,count_gas-1 do plots,v1[i,*],v2[i,*],color=gasCM[i,*]
  
    ; scale bar
    len = 100.0 ;ckpc
    cgPlot,[boxCen[0]-zoomSize/2.2,boxCen[0]-zoomSize/2.2+len],[boxCen[1]+zoomSize/2.2,boxCen[1]+zoomSize/2.2],$
      line=0,color=cgColor('dark gray'),/overplot
    cgText,mean([boxCen[0]-zoomSize/2.2,boxCen[0]-zoomSize/2.2+len]),boxCen[1]+zoomSize/2.4,$
      string(len,format='(i3)')+' ckpc',alignment=0.5,charsize=!p.charsize-0.6,color=cgColor('dark gray')
      
    ; halo mass
    strHalo = "log(M) = "+string(codeMassToLogMsun(mt.hMass[0,hInd]),format='(f4.1)')
    cgText,0.75,0.04,strHalo,alignment=0.5,color=cgColor('dark gray'),/normal,charsize=!p.charsize-0.4
        
    ; temperature colorbar
    !p.thick = 1.0
    !p.charsize = 0.8
    !x.thick = 1.0
    !y.thick = 1.0
    colorbar,position=[0.02,0.1,0.076,0.4],divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
    cgText,0.25,0.0375,textoidl("log T_{gas} [K]"),alignment=0.5,color=cgColor('white'),/normal
    cgText,0.115,0.036,'4',alignment=0.5,color=cgColor('white'),/normal
    cgText,0.385,0.036,'7',alignment=0.5,color=cgColor('white'),/normal
    
  end_PS, pngResize=60, im_options='-negate', /deletePS

  ; make a movie of an orbit if requested
  for fnOrbit=0,nFramesOrbit-1 do begin    
    ; determine current rotation and create v1,v2 (screen plane projections)
    curRot = (360.0/nFramesOrbit) * fnOrbit

    print,' ['+string(fnOrbit+1,format='(i4)')+'] orbit frame ['+$
      str(fnOrbit)+' of '+str(nFramesOrbit)+']'
    
    ; plot
    start_PS, plotPath + 'frame_'+str(fnOrbit+1)+'.eps', xs=pssize[0], ys=pssize[1]
    
      !p.thick = 0.6
      
      ; clip positions at box dimensions (left,right)
      v2 = reform(gasPos[axes[0],*,*]) ; keep axes[0] always constant
      v1 = reform(gasPos[axes[1],*,*]) * cos(curRot*!dtor) + $
           reform(gasPos[axes[2],*,*]) * sin(curRot*!dtor)
          
      ;w = where(v1 gt boxCen[0]+zoomSize/2.0,count)
      ;if count gt 0 then v1[w] = !values.f_nan
      ;w = !NULL
      
      ; color table
      loadColorTable,'helix',/reverse
      ; gasCM already made
      
      ; fixed viewport
      xrange = [boxCen[0]-zoomSize/2.0,boxCen[0]+zoomSize/2.0]
      yrange = [boxCen[1]-zoomSize/2.0,boxCen[1]+zoomSize/2.0] ; top+bottom halves
    
      ; left: gas
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=5,ystyle=5,$
             xtitle="",ytitle="",title="",position=[0,0,0.5,1],xtickname=replicate(' ',10),$
             ytickname=replicate(' ',10)
           
      ; left: draw rvir and rvir/10
      tvcircle,mt.hVirRad[0,hInd]/10,0,0,cgColor('black'),thick=0.8,/data
      tvcircle,mt.hVirRad[0,hInd],0,0,cgColor('black'),thick=0.8,/data  
    
      ; draw gas trajectories
      for i=0,count_gas-1 do plots,v1[i,*],v2[i,*],color=gasCM[i,*]
    
      ; scale bar
      len = 100.0 ;ckpc
      cgPlot,[boxCen[0]-zoomSize/2.2,boxCen[0]-zoomSize/2.2+len],[boxCen[1]+zoomSize/2.2,boxCen[1]+zoomSize/2.2],$
        line=0,color=cgColor('dark gray'),/overplot
      cgText,mean([boxCen[0]-zoomSize/2.2,boxCen[0]-zoomSize/2.2+len]),boxCen[1]+zoomSize/2.4,$
        string(len,format='(i3)')+' ckpc',alignment=0.5,charsize=!p.charsize-0.6,color=cgColor('dark gray')
        
      ; halo mass
      strHalo = "log(M) = "+string(codeMassToLogMsun(mt.hMass[0,hInd]),format='(f4.1)')
      cgText,0.75,0.04,strHalo,alignment=0.5,color=cgColor('dark gray'),/normal,charsize=!p.charsize-0.4
          
      ; particle type names
      ;cgText,0.04,0.04,"Gas",alignment=0.5,color=cgColor('dark gray'),/normal,charsize=!p.charsize-0.4
      ;cgText,0.96,0.04,"DM",alignment=0.5,color=cgColor('dark gray'),/normal,charsize=!p.charsize-0.4
      
      ; temperature colorbar
      !p.thick = 1.0
      !p.charsize = 0.8
      !x.thick = 1.0
      !y.thick = 1.0
      colorbar,position=[0.02,0.1,0.076,0.4],divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
      cgText,0.25,0.0375,textoidl("log T_{gas} [K]"),alignment=0.5,color=cgColor('white'),/normal
      cgText,0.115,0.036,'4',alignment=0.5,color=cgColor('white'),/normal
      cgText,0.385,0.036,'7',alignment=0.5,color=cgColor('white'),/normal

    end_PS, pngResize=60, im_options='-negate', /deletePS
  endfor

end
