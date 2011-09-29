; spiralsPatSpeed.pro
; fourier pattern speed routines
; dnelson may.2011

; barSpeed()

function barSpeed, simName, old=old

  units = getUnits()

  ; sim config
  basePath = '/n/home07/dnelson/spirals/lambdacrit/sims/'

  filePath    = basePath + simName + "/"
  workingPath = filePath + "/output3/"
  
  nSnaps      = max([n_elements(file_search(filePath+"output/Stars_*_0")),$
                     n_elements(file_search(filePath+"output/Stars_*.bin")),$
                     n_elements(file_search(filePath+"output/snap_*"))])

  h2 = loadSimParams(simName)

  saveFileName = filePath + 'barSpeed.radii.sav'

  ; config
  errTol         = 0.03
  radii          = [0.25,0.5,0.75,1.0,1.5,2.0,2.5,3.0,4.0,5.0] * h2.h
  nAngularBins   = 180 ;360

  ; arrays
  meanPhase = fltarr(nSnaps,n_elements(radii))
  errPhase  = fltarr(nSnaps,n_elements(radii))
  priPhase  = fltarr(nSnaps,n_elements(radii))
  secPhase  = fltarr(nSnaps,n_elements(radii))
  times     = fltarr(nSnaps)

  dTheta = 2.0 * !pi / nAngularBins 
  angularMidBins  = findgen(nAngularBins) / nAngularBins * 2.0 * !pi ;[0,+2pi]
  angularMidBins += dTheta/2.0

  if not (file_test(saveFileName)) then begin
  
    ; find error optimized annulus widths
    binWidths = fltarr(n_elements(radii))
    
    if not keyword_set(old) then begin
      xyz = loadStars(filePath,1)
      r     = reform(sqrt(xyz.x^2.0 + xyz.y^2.0))
    endif else begin
      old = loadSnapshot(filePath+"output/snap_",1,s_pos,s_vel,s_id,c_pos,c_vel,c_id)
      r   = reform(sqrt(s_pos[0,*]^2.0 + s_pos[1,*]^2.0))
    endelse
    
    for i=0,n_elements(radii)-1 do begin
      binWidths[i] = 0.1 ;kpc
      w = where(abs(r-radii[i]) lt binWidths[i],count)
      estError = 1.0 / sqrt(count/nAngularBins)
    
      while (abs(estError-errTol) gt errTol/100) do begin
        estError = 1.0 / sqrt(count/nAngularBins)
        binWidths[i] *= (1.0 + (estError-Errtol)/2)
      
        w = where(abs(r-radii[i]) lt binWidths[i],count)
      endwhile
      print,'i binwidth N errEst',i,binWidths[i],count,estError
    endfor
    
    ; loop over all time
    for m=0,nSnaps-1,1 do begin ;bar_a=stepsize of 10, bar_b/LC=stepsize of 1
    ;m=1600
  
    ; azimuthal histogram
    h1d = fltarr(nAngularBins,n_elements(radii))
    
    ; load
    if not keyword_set(old) then begin
      xyz      = loadStars(filePath,(m+1))
      mcs      = loadMCs(filePath,(m+1),h)
      times[m] = h[1] * 1000.0 ;Myr
    
      r     = reform(sqrt(xyz.x^2.0 + xyz.y^2.0))
      arctan,xyz.x,xyz.y,theta_rad,theta_deg ;[0,2pi]
    endif else begin
      old = loadSnapshot(filePath+"output/snap_",m,s_pos,s_vel,s_id,c_pos,c_vel,c_id)
      times[m] = old.time * 1000.0 ;Myr
      
      r   = reform(sqrt(s_pos[0,*]^2.0 + s_pos[1,*]^2.0))
      arctan,s_pos[0,*],s_pos[1,*],theta_rad,theta_deg ;[0,2pi]
    endelse
  
    ; loop over each desired radius
    for i=0,n_elements(radii)-1 do begin
      
      ; bin surface density wrt phase
      w = where(abs(r-radii[i]) lt binWidths[i],count)
  
      h1d[*,i] = (histogram(theta_rad[w],binsize=dTheta))[0:nAngularBins-1]
      h1d[*,i] /= mean(h1d[*,i]) ;overdensity
    
      ; duplicate profile once on each side
      gx = fltarr(nAngularBins*3)
      gy = fltarr(nAngularBins*3)
      
      for j=0,2 do begin
        gx[nAngularBins*j:nAngularBins*(j+1)-1] = angularMidBins + 2*!pi*(j-1)
        gy[nAngularBins*j:nAngularBins*(j+1)-1] = h1d[*,i]
      endfor
    
      ; fit maximum with gaussian
      w = where(gy eq max(gy),count)
      maxInd = w[count/2]
      
      gxFit = gx[maxInd-nAngularBins/4:maxInd+nAngularBins/4] ;[max-pi/2,max+pi/2] subset
      gyFit = gy[maxInd-nAngularBins/4:maxInd+nAngularBins/4]
    
      res = gaussfit(gxFit,gyFit,A)
      priAmp        = A[0]
      priPhase[m,i] = A[1]
    
      ; fit secondary maximum with gaussian
      if (priPhase[m,i] lt !pi) then $
        gyFit = gy[maxInd+(nAngularBins/4):maxInd+(3*nAngularBins/4)-1]
      if (priPhase[m,i] ge !pi) then $
        gyFit = gy[maxInd-(3*nAngularBins/4):maxInd-(nAngularBins/4)-1]
      
      w = where(gyFit eq max(gyFit),count)
      
      if (priPhase[m,i] lt !pi) then $
        minInd = maxInd + (nAngularBins/4) + w[0] ; convert back to global gy[] index
      if (priPhase[m,i] ge !pi) then $
        minInd = maxInd - (3*nAngularBins/4) + w[0]
      
      gxFit = gx[minInd-nAngularBins/4:minInd+nAngularBins/4] ;[max-pi/2,max+pi/2] subset
      gyFit = gy[minInd-nAngularBins/4:minInd+nAngularBins/4]  
      
      res = gaussfit(gxFit,gyFit,A)
      secAmp        = A[0]
      secPhase[m,i] = A[1]
      
      ; avg two maxima -> angle of "normal" of bar in disk plane
      meanPhase[m,i] = (priPhase[m,i]+secPhase[m,i])/2.0
      errPhase[m,i]  = abs(abs(priPhase[m,i]-secPhase[m,i])-!pi)
      print,'m i pri sec mean err',m,i,priPhase[m,i],secPhase[m,i],meanPhase[m,i],errPhase[m,i]
    
    endfor ;i
    
    ; plot
    if (m mod 5 eq 0) then begin ;bar=mod100, spiral=mod40, bar_b=mod5
    start_PS, workingPath+"barSpeed."+str(m)+".eps"    
      fsc_plot,[0],[0],/nodata,xrange=[0,2.5*!pi],yrange=[0,max(h1d)*1.1],$
               xtitle="Azimuthal Phase",ytitle=textoidl("\Sigma(\theta) / \Sigma_0"),/xs,/ys,$
               xticks=5,xtickv=[0,!pi/2,!pi,3*!pi/2,2*!pi,5*!pi/2],$
               xtickname=['0',textoidl('\pi/2'),textoidl('\pi'),textoidl('3\pi/2'),$
                          textoidl('2\pi'),textoidl('5\pi/2')]
      
      for i=0,n_elements(radii)-1,1 do begin
        fsc_plot,angularMidBins,h1d[*,i],color=fsc_color(units.colors[i]),/overplot
        fsc_plot,[2*!pi-dTheta/2.0,angularMidBins+2*!pi],[h1d[nAngularBins-1,i],h1d[*,i]],$
                 color=fsc_color(units.colors[i]),line=1,/overplot ;wrap
      
        ; gaussian fits
        fsc_plot,[priPhase[m,i],priPhase[m,i]],[max(h1d)*0.05,max(h1d)*0.15],$
                 line=0,color=fsc_color(units.colors[i+1]),/overplot
        fsc_plot,[secPhase[m,i],secPhase[m,i]],[max(h1d)*0.05,max(h1d)*0.15],$
                 line=0,color=fsc_color(units.colors[i+1]),/overplot
      endfor
    end_PS
    
    endif ;m mod X
    
    endfor ;m

    ; save/restore
    save,priPhase,secPhase,meanPhase,errPhase,times,radii,filename=saveFileName
  endif else begin
    restore,saveFileName
  endelse 
  
  nSnaps = n_elements(times) ;override nSnaps in case we have added more snaps but not in save

  ; calculate bar angle
  SA = fltarr(n_elements(times))

  for m=0,nSnaps-1,1 do begin
    SA[m] = (meanPhase[m,0] + meanPhase[m,1]) / 2.0 ;avg two inner measurements
    SA[m] -= !pi/2.0   ;normal to parallel
    SA[m] *= (180/!pi) ;rad to deg
  endfor

  ; calculate bar angular frequency with radius 
  angMean    = fltarr(n_elements(radii))
  angMeanErr = fltarr(n_elements(radii))
  frac       = fltarr(n_elements(radii))
  
  for j=0,n_elements(radii)-1 do begin
  
    priPhaseRad  = priPhase[*,j] ;nz,j (bar)
    secPhaseRad  = secPhase[*,j]
    meanPhaseRad = meanPhase[*,j]
    errPhaseRad  = errPhase[*,j]
    
    ; maxima
    plotPhase = fltarr(n_elements(meanPhaseRad))
    flip = 0
    
    ; unwrap
    for i=0,n_elements(meanPhaseRad)-2 do begin
      plotPhase[i] = meanPhaseRad[i] + flip*!pi
      
      if (meanPhaseRad[i] ge meanPhaseRad[i+1]) then begin
        flip += 1
      endif
    endfor
    
    plotPhase2 = plotPhase mod (2*!pi)

    ; calculate derivative
    angFreq            = deriv(times,plotPhase) ;plotPhase is monotonic increasing
    angFreq_err        = derivsig(times,plotPhase,0.0,errPhase)

    angFreq_kmskpc     = angFreq * units.kpc_in_km/units.s_in_Myr
    angFreq_err_kmskpc = angFreq_err * units.kpc_in_km/units.s_in_Myr

    ; calculate mean angFreq and error
    w = where(angFreq_kmskpc gt 0 and angFreq_kmskpc lt 60,count)
    
    angMean[j]    = mean(angFreq_kmskpc[w])
    angMeanErr[j] = stddev(angFreq_kmskpc[w])
    frac[j]       = float(count)/n_elements(angFreq_kmskpc)
    
    ;print,'j angMean angMeanErr frac',j, angMean[j], angMeanErr[j], frac[j]

    ; plot theta(t)
    ;if 0 then begin ;for patSpeedSpiralMult suppress plotting
    start_PS, workingPath+"barSpeed.theta_t."+str(j)+".eps"    
      fsc_plot,[0],[0],/nodata,xrange=[0,round(max(times)/100)*100],yrange=[-1.05,1.05],$
               xtitle="Time [Myr]",ytitle="Sin (Bar Normal Angle) [rad]",/xs,/ys
     
      ;fsc_plot,times,priPhase,line=0,color=fsc_color('red'),/overplot
      ;fsc_plot,times,secPhase,line=0,color=fsc_color('blue'),/overplot
      ;fsc_plot,times,meanPhase,line=0,/overplot,color=fsc_color('orange')
      ;fsc_plot,times,plotPhase,line=0,thick=2.0,/overplot
      fsc_plot,times,sin(plotPhase2),line=0,thick=2.0,/overplot
    end_PS
    
    ; plot dtheta/dt(t) = angular velocity / radius = angular frequency
    start_PS, workingPath+"barSpeed.angfreq."+str(j)+".eps"    
      fsc_plot,[0],[0],/nodata,xrange=[0,round(max(times)/100)*100],yrange=[0,120.0],$
               xtitle="Time [Myr]",ytitle="Angular Frequency [km/s/kpc]",/xs,/ys
      
      fsc_plot,times,angFreq_kmskpc,psym=16,symsize=0.6,/overplot
      
      ;errors
      ;for i=0,n_elements(errPhase)-2 do begin
      ;  fsc_plot,[times[i],times[i]],$
      ;           [angFreq_kmskpc[i]-angFreq_err_kmskpc[i],angFreq_kmskpc[i]+angFreq_err_kmskpc[i]],$
      ;           line=0,/overplot
      ;endfor
      
      ;t>500Myr avg line (bar)
      ;angMean = mean(angFreq_kmskpc[where(where(angFreq_kmskpc gt 0) ge min(where(times ge 500.0)))])
      ;fsc_plot,[0,round(max(times)/100)*100],[angMean,angMean],line=1,/overplot
      
      ; 0<angFreq<60 avg line (LC old algorithm)
      fsc_plot,[0,round(max(times)/100)*100],[angMean[j],angMean[j]],line=1,/overplot
      fsc_plot,[0,round(max(times)/100)*100],[angMean[j]+angMeanErr[j],angMean[j]+angMeanErr[j]],$
               line=1,color=fsc_color('gray'),/overplot
      fsc_plot,[0,round(max(times)/100)*100],[angMean[j]-angMeanErr[j],angMean[j]-angMeanErr[j]],$
               line=1,color=fsc_color('gray'),/overplot
    end_PS
    
    ; plot dtheta/dt(t) full range
    start_PS, workingPath+"barSpeed.angfreq.full."+str(j)+".eps"    
      fsc_plot,[0],[0],/nodata,xrange=[0,round(max(times)/100)*100],yrange=[-300,2000],$
               xtitle="Time [Myr]",ytitle="Angular Frequency [km/s/kpc]",/xs,/ys
      
      fsc_plot,times,angFreq_kmskpc,psym=16,symsize=0.6,/overplot
    end_PS
    ;endif ;0
  
  endfor ;j
  
  ; plot angMean vs radius
  fracThreshold = 0.05
  
  w = where(frac ge fracThreshold)
  
  ;if 0 then begin ;patSpeedSeriesMult
  start_PS, workingPath+"barSpeed.angMean.vs.rad.eps"    
    fsc_plot,[0],[0],/nodata,xrange=[0,round(max(radii/h2.h)/2)*2],yrange=[0,60],$
             xtitle="Radius / Disk Scale Length [kpc]",ytitle="Angular Frequency [km/s/kpc]",/xs,/ys
    
    fsc_plot,radii[w]/h2.h,angMean[w],psym=16,symsize=0.8,/overplot
    
    ; errors
    for i=0,n_elements(w)-1 do begin
      fsc_plot,[radii[w[i]],radii[w[i]]]/h2.h,$
               [angMean[w[i]]-angMeanErr[w[i]],angMean[w[i]]+angMeanErr[w[i]]],$
               line=0,/overplot
    endfor
  end_PS
  ;endif ;0
  
  r = {radii:radii,SA:SA,angMean:angMean,angMeanErr:angMeanErr,frac:frac}
  
  return,r
end

; barRadProfile(): use bar normal angle with time to find bar radial extent

pro barRadProfile, simName, old=old

  units = getUnits()

  ; sim config  
  basePath    = '/n/home07/dnelson/spirals/lambdacrit/sims/'
  filePath    = basePath + simName + "/"
  workingPath = filePath + "output3/"
  
  h2 = loadSimParams(simName)

  saveFileName = filePath + 'barRadProfile.sav'

  ; restore barSpeed save
  fileName = filePath + "barSpeed.radii.sav"
  restore,fileName
  
  ; config
  nRad   = 100
  maxRad = 6.0 ;kpc
  yWidth = 0.2  ;kpc
  startTime = 600.0 ;Myr, of "stable" bar period
  
  break_radius = 3.6  ;kpc
  frac_amp_ext = 0.04 ;fraction of r=0 density to call extent of bar
  
  nSnaps    = n_elements(times)  
  
  if not (file_test(saveFileName)) then begin
    ; arrays
    radSlices = fltarr(nRad,nSnaps)
    
    bar_ext_derived   = fltarr(nSnaps)
    bar_ext_dererr    = fltarr(nSnaps)
    bar_gaussfit_fwhm = fltarr(nSnaps)
    bar_gaussfit_cent = fltarr(nSnaps)
  
    ; loop over snapshots
    for m=0,nSnaps-1,1 do begin ;bar_a=stepsize of 10, bar_b=stepsize of 1
  
      ; load
      if not keyword_set(old) then begin
        xyz = loadStars(filePath,1)
        x = xyz.x
        y = xyz.y
      endif else begin
        old = loadSnapshot(filePath+"output/snap_",1,s_pos,s_vel,s_id,c_pos,c_vel,c_id)
        x = reform(s_pos[0,*])
        y = reform(s_pos[1,*])
      endelse
      
      r = reform(sqrt(x^2.0 + y^2.0))
      arctan,x,y,theta_rad,theta_deg ;[0,2pi]
        
      ; calc rotation angle
      SA = (meanPhase[m,0] + meanPhase[m,1]) / 2.0 ;avg two inner measurements
      SA += !pi/2.0   ;normal to parallel
      SA *= (180/!pi) ;rad to deg
      
      ; rotate to align bar with x-axis
      xr = x * cos(-1.0*SA*!dtor) - y * sin(-1.0*SA*!dtor)
      yr = y * cos(-1.0*SA*!dtor) + x * sin(-1.0*SA*!dtor)
      
      ; slice and histogram
      w = where(abs(yr) lt yWidth,count)
      hist1d = histogram(xr[w],nbins=nRad,min=-maxRad,max=maxRad,locations=radPts)
      
      radPts += (2*maxRad)/nRad/2.0  ;adjust to midbins
      hist1d /= float(max(hist1d))          ;norm to r=0 value
  
      ; fit amplitude level to calculate extent
      w = where(abs(radPts-break_radius) eq min(abs(radPts-break_radius)))
      
      w = where(radPts ge 0.0)
      w1 = where(abs(hist1d[w]-frac_amp_ext) eq min(abs(hist1d[w]-frac_amp_ext)))
      rad1 = radPts[w[w1]]
      
      w = where(radPts lt 0.0)
      w2 = where(abs(hist1d[w]-frac_amp_ext) eq min(abs(hist1d[w]-frac_amp_ext)))
      rad2 = radPts[w[w2]]
      
      bar_ext_derived[m] = (rad1 + abs(rad2)) / 2.0
      bar_ext_dererr[m]  = abs(rad1 + rad2)
      
      ;gaussian+constant fit
      res = gaussfit(radPts,hist1d,A,nterms=4)
      bar_gaussfit_fwhm[m] = 2.0 * sqrt(2.0 * alog(2.0)) * A[2]
      bar_gaussfit_cent[m] = A[1]
  
      ; plot pos/neg slices
      if (m mod 5 eq 0) then begin ;bar_a=mod100, bar_b=mod5
      start_PS, workingPath+"barRadPro."+str(m)+".eps"    
        fsc_plot,[0],[0],/nodata,xrange=[0,maxRad],yrange=[0,1.0],$
                 xtitle="Radius [kpc]",ytitle=textoidl("\Sigma(r) / \Sigma_0"),/xs,/ys
        
          ;r>0
          w = where(radPts ge 0.0)
          fsc_plot,radPts[w],res[w],/overplot,color=fsc_color(units.colors[3]),thick=1.5 ;gaussian fit  
          fsc_plot,radPts[w],hist1d[w],thick=3.0,/overplot
          
          ;r<0
          w = where(radPts lt 0.0)
          fsc_plot,-1.0*radPts[w],res[w],/overplot,color=fsc_color(units.colors[3]),thick=1.5 ;gaussian fit 
          fsc_plot,-1.0*radPts[w],hist1d[w],thick=3.0,/overplot
          
          ;extents
          fsc_plot,[break_radius,break_radius],[0.0,1.0],/overplot,$
                   line=1,color=fsc_color(units.colors[1])
          fsc_plot,[bar_ext_derived[m],bar_ext_derived[m]],[0.0,1.0],/overplot,$
                   line=1,color=fsc_color(units.colors[2])
    
      end_PS
      
      endif ;m mod X
  
      print,m,bar_ext_derived[m],bar_ext_dererr[m],bar_gaussfit_fwhm[m],bar_gaussfit_cent[m]
    endfor

    ; save/restore
    save,bar_ext_derived,bar_ext_dererr,bar_gaussfit_fwhm,bar_gaussfit_cent,filename=saveFileName
  endif else begin
    restore,saveFileName
  endelse

  ; plot 
  w1 = where(bar_ext_derived ne 0.0 and times ge startTime)
   
  start_PS, workingPath+"barRadPro.gauss.vs.amp.eps"    
    fsc_plot,[0],[0],/nodata,xrange=[min(bar_gaussfit_fwhm[w1])*0.9,max(bar_gaussfit_fwhm[w1])*1.1],$
             yrange=[min(bar_ext_derived[w1])*0.9,max(bar_ext_derived[w1])*1.1],$
             xtitle="Gaussian Fit FWHM [kpc]",ytitle="Amplitude Fit Extent [kpc]",/xs,/ys
             
    fsc_plot,bar_gaussfit_fwhm[w1],bar_ext_derived[w1],psym=16,/overplot
      
  end_PS
  
  w2 = where(bar_ext_derived ne 0.0)
  
  start_PS, workingPath+"barRadPro.extents.wtime.eps"    
    fsc_plot,[0],[0],/nodata,xrange=[0,round(max(times)/100)*100],$
             yrange=[-1.5,6.3],$
             xtitle="Time [Myr]",ytitle="Bar Extent / Center [kpc]",/xs,/ys
   
    ;break_radius 
    fsc_plot,[startTime,2000.0],[break_radius,break_radius],line=1,/overplot  
    
    ;bar extents
    fsc_plot,times[w2],bar_ext_derived[w2],psym=-16,symsize=0.8,color=fsc_color(units.colors[2]),/overplot
    fsc_plot,times[w2],bar_gaussfit_fwhm[w2],psym=-16,symsize=0.8,color=fsc_color(units.colors[3]),/overplot
    
    ;t>600 Myr fits
    fsc_plot,[startTime,2000.0],[mean(bar_ext_derived[w1]),mean(bar_ext_derived[w1])],line=1,/overplot
    fsc_plot,[startTime,2000.0],[mean(bar_gaussfit_fwhm[w1]),mean(bar_gaussfit_fwhm[w1])],line=1,/overplot  
    
    ;gaussfit center
    fsc_plot,times[w2],bar_gaussfit_cent[w2],psym=-16,symsize=0.8,color=fsc_color(units.colors[4]),/overplot
      
  end_PS
  
  ; plot rotation angle
  SA = fltarr(nSnaps)

  for m=0,nSnaps-1,1 do begin
    SA[m] = (meanPhase[m,0] + meanPhase[m,1]) / 2.0 ;avg two inner measurements
    SA[m] -= !pi/2.0   ;normal to parallel
    SA[m] *= (180/!pi) ;rad to deg
  endfor
  
  w3 = where(SA ne 0.0)
  
  start_PS, workingPath+"barRadPro.SA.wtime.eps"    
    fsc_plot,[0],[0],/nodata,xrange=[startTime,round(max(times)/100)*100],$
             yrange=[0.0,180.0],$
             xtitle="Time [Myr]",ytitle="Bar Angle [deg]",xs=9,/ys
   
    ;bar extents
    fsc_plot,times[w3],SA[w3],psym=-16,symsize=0.5,color=fsc_color(units.colors[0]),/overplot
    fsc_axis,0.0,180.0,0.0,xaxis=1,/xs,xrange=minmax(w1)
  end_PS

  print,interpol(SA,times,[600.0,650.0,700.0,750.0,800.0,850.0,900.0,950.0,1000.0])
  print,SA[60],SA[61]
            
end

; mpModel(): omega_p=const (1 parameter)
function mpModelConst, X, P, _EXTRA=fargs

  y = fltarr(n_elements(X))
  
  for i=0,n_elements(X)-1 do begin
    ;start integrand
    int = fargs.xMid * fargs.surfDensDY[*,i] - X[i] * fargs.surfDensDX[*,i]
    
    ;omega_p model
    omega_p = P[0]
    
    ;finish integrand
    int *= omega_p
    
    y[i] = int_tabulated(fargs.xMid,int)
  endfor
  
  return,y
end

; mpModel(): omega_p = 2 segment broken constant at a break radius (3 parameter, 1 fixed)
function mpModelBar, X, P, _EXTRA=fargs

  y       = fltarr(n_elements(X))
  omega_p = fltarr(n_elements(X))
  
  ; P[0]=const_1 (bar), P[1]=const_2 (outer), P[2]=break_radius
  
  for i=0,n_elements(X)-1 do begin
    ;start integrand
    int = fargs.xMid * fargs.surfDensDY[*,i] - X[i] * fargs.surfDensDX[*,i]
    
    r       = sqrt(fargs.xMid^2.0 + X[i]^2.0)
    omega_p = fltarr(n_elements(fargs.xMid))
    
    ;omega_p model
    w = where(r lt P[2],count)
    if (count ne 0) then $
      omega_p[w] = P[0]
    w = where(r ge P[2],count)
    if (count ne 0) then $
      omega_p[w] = P[1];/r[w]^2.0
      
    ;finish integrand
    int *= omega_p
    
    y[i] = int_tabulated(fargs.xMid,int)
  endfor
  
  return,y
end

; mpModel(): omega_p(r) = a/r^b (2 parameter)
function mpModelSpiralPowerlaw, X, P, _EXTRA=fargs

  y       = fltarr(n_elements(X))
  omega_p = fltarr(n_elements(X))
  
  ; P[0]=a (amplitude), P[1]=b (plaw index), P[2]=c (const offset)
  
  for i=0,n_elements(X)-1 do begin
    ;start integrand
    int = fargs.xMid * fargs.surfDensDY[*,i] - X[i] * fargs.surfDensDX[*,i]
    
    r       = sqrt(fargs.xMid^2.0 + X[i]^2.0)
    omega_p = fltarr(n_elements(fargs.xMid))
    
    ;omega_p model
    omega_p = P[0] / (r / fargs.h)^(P[1]) ; + P[2]
      
    ;finish integrand
    int *= omega_p
    
    y[i] = int_tabulated(fargs.xMid,int)
  endfor
  
  return,y
end

; doMpfit()

function doMpfit, y_r_max, yWidth, break_radius, xBins, sgN, h2, s_pos, s_vel, $
                  SA=SA, inc=inc, bar=bar

  ; derivative config
  N = floor(y_r_max / yWidth) + 1
  
  yStarts = findgen(N)/(N-1) * y_r_max - y_r_max/2.0 - yWidth/2.0 ;y_j
  yMids   = yStarts + yWidth/2.0

  xMinMax = [-10.0,10.0]*h2.h ; whole distribution
  
  sgOrder = 6
  
  dx = (xMinMax[1] - xMinMax[0]) / xBins
  xMid = findgen(xBins)/(xBins-1) * (xMinMax[1] - xMinMax[0]) + xMinMax[0] + dx/2.0
  
  ; precompute surface density and derivative fields
  surfDens    = fltarr(xBins,N)
  surfDensDX  = fltarr(xBins,N)
  surfDensDY  = fltarr(xBins,N)
  vyMean      = fltarr(xBins,N)
  vyDY        = fltarr(xBins,N)
  
  surfDensErr   = fltarr(xBins,N)
  surfDensDXErr = fltarr(xBins,N)
  surfDensDYErr = fltarr(xBins,N)
  vyMeanErr     = fltarr(xBins,N)
  vyDYErr       = fltarr(xBins,N)
  
  ; separate components
  x = reform(s_pos[0,*])
  y = reform(s_pos[1,*])
  z = reform(s_pos[2,*])
  
  vx = reform(s_vel[0,*])
  vy = reform(s_vel[1,*])
  vz = reform(s_vel[2,*])
  
  if not keyword_set(SA) then SA = 0.0
  if not keyword_set(inc) then inc = 0.0

  ; slice angle rotation (CCW, so rot x,y coords and [vx,vy] vector CW)
  x1 = x * cos(-1.0*SA*!dtor) - y * sin(-1.0*SA*!dtor)
  y1 = y * cos(-1.0*SA*!dtor) + x * sin(-1.0*SA*!dtor)
  
  vx1 = vx * cos(-1.0*SA*!dtor) - vy * sin(-1.0*SA*!dtor)
  vy1 = vy * cos(-1.0*SA*!dtor) + vx * sin(-1.0*SA*!dtor)

  ; sky plane inclination (about x-axis, with y>0 away from observer)
  xr = x1
  yr = y1 * cos(-1.0*inc*!dtor) + z  * sin(-1.0*inc*!dtor)
  yr /= cos(-1.0*inc*!dtor) ; correct for known inclination (zero error)

  if (inc ne 0.0) then begin
    vyr = vy1 * sin(-1.0*inc*!dtor) + vz * cos(-1.0*inc*!dtor) ;contour LOS/radial vel
    vyr /= sin(-1.0*inc*!dtor) ; de-project for disk plane component
  endif else begin
    vyr = vy1
  endelse

  ;loop over y_j
  for i=0,N-1 do begin

    yMin = yStarts[i]
    yMax = yMin + yWidth
    
    ;select y subset
    w1 = where(yr ge yMin and yr lt yMax, count1)
    
    if (count1 eq 0) then continue
    
    x1  = xr[w1]
    vy1 = vyr[w1]
    
    ;loop over x bins
    for j=0,xBins-1 do begin
      xMin = xMid[j] - dx/2.0
      xMax = xMid[j] + dx/2.0
      
      ;find stars in bin
      w = where(x1 ge xMin and x1 lt xMax, count)
      
      if (count ne 0) then begin
        surfDens[j,i]    = count ;divide by area later
        vyMean[j,i]      = mean(vy1[w])
        ;vyMeanErr[j,i]   = stddev(vy1[w])
      endif
    endfor
    
  endfor

  ; smoothing
  surfDens = convol(surfDens,savgol(sgN,sgN,0,sgOrder),/edge_truncate)
  vyMean   = convol(vyMean,  savgol(sgN,sgN,0,sgOrder),/edge_truncate)

  ; sanity check
  w = where(surfDens lt 0,count)
  if (count ne 0) then $
    surfDens[w] = 0.0

  w = where(finite(vyMeanErr) eq 0,count)
  if (count ne 0) then $
    vyMeanErr[w] = 0.0
  vyMeanErr = convol(vyMeanErr,  savgol(sgN,sgN,0,sgOrder),/edge_truncate)

  ; convert binned counts to surface density
  surfDens    /= (dx*yWidth)
  
  ; poisson errors on smoothed surface density
  ;surfDensErr = 1.0/sqrt(surfDens)
  surfDensErr = sqrt(surfDens)
  
  w = where(finite(surfDensErr) eq 0,count)
  if (count ne 0) then $
    surfDensErr[w] = 0.0

  ; numerical derivatives
  for i=0,N-1 do begin
    surfDensDX[*,i]    = deriv(xMid,surfDens[*,i])
    surfDensDXErr[*,i] = derivsig(xMid,surfDens[*,i],surfDensErr[*,i],0.0)
  endfor
  
  for j=0,xBins-1 do begin
    surfDensDY[j,*]    = deriv(yStarts,surfDens[j,*])
    surfDensDYErr[j,*] = derivsig(yStarts,surfDens[j,*],surfDensErr[j,*],0.0)
    vyDY[j,*]          = deriv(yStarts,vyMean[j,*])
    vyDYErr[j,*]       = derivsig(yStarts,vyMean[j,*],vyMeanErr[j,*],0.0)
  endfor

  ; calculate f(y_j)
  f    = fltarr(N)
  
  for i=0,N-1 do begin
    fVals = surfDensDY[*,i] * vyMean[*,i] + surfDens[*,i] * vyDY[*,i] ;product rule
    f[i]  = int_tabulated(xMid,fVals)
  endfor
  
  ; calculate ferr(y_j)
  if 0 then begin
  ferr = fltarr(N)
  nR = 1000
  
  for i=0,N-1 do begin
    ; randomly sample integral nR times
    ints = fltarr(nR)
    
    for j=0,nR-1 do begin
      R = randomn(seed,4)         ;normal distribution, mean=0 stddev=1
      ;R = randomu(seed,4)*2.0-1.0 ;uniform distribution [-1,1]
      
      fValsErr = (surfDensDY[*,i]+R[0]*surfDensDYErr[*,i]) * (vyMean[*,i]+R[1]*vyMeanErr[*,i]) + $
                  (surfDens[*,i]+R[2]*surfDensErr[*,i]) * (vyDY[*,i]+R[3]*vyDYErr[*,i])
  
      ints[j] = int_tabulated(xMid,fValsErr)
    endfor
    
    ; report error as stddev of realizations
    ferr[i] = stddev(ints)
  endfor
  endif ;0
    
  ; functargs
  fargs = {surfDensDX:surfDensDX,surfDensDY:surfDensDY,xMid:xMid,h:h2.h}

  if keyword_set(bar) then begin
    ; mpfit model (bar)
    ; -----------------
    nArgs = 3
    err   = 1.0 ;uniform weighting
    
    pi = replicate({parname:'', value:0.0, relstep: 0.0, fixed:0, limited:[0,0], limits:[0.D,0.D]},nArgs)
    pi[0].parname = 'OMEGA_P_BAR'
    pi[1].parname = 'OMEGA_P_OUTER'
    pi[2].parname = 'BREAK_RADIUS'
    
    pi[2].fixed = 1
    
    pi[*].value = [10.0,10.0,break_radius] ;starting, y_r_max/4.0
    
    ;pi[0].relstep = 0.5
    ;pi[1].relstep = 0.5
    ;pi[2].relstep = 0.20
    
    ; do fit
    res = mpfitfun('mpModelBar',yMids,f,err,functargs=fargs,parinfo=pi,$
                   bestnorm=bestnorm,dof=dof,perror=perror,/quiet)
                   
  endif else begin
    ; mpfit model (spiral)
    ; --------------------
    nArgs = 2
    err   = 1.0
    ;err = ferr
    
    pi = replicate({parname:'', value:0.0, relstep: 0.0, fixed:0, limited:[0,0], limits:[0.D,0.D]},nArgs)
    pi[0].parname = 'OMEGA_P_AMP'
    pi[1].parname = 'OMEGA_P_SLOPE'
    ;pi[2].parname = 'OMEGA_P_CONST'
    
    pi[*].value = [10.0,1.0] ;starting
   
    ;pi[0].relstep = 0.1
    ;pi[1].relstep = 0.01
    
    ; do fit
    res = mpfitfun('mpModelSpiralPowerlaw',yMids,f,err,functargs=fargs,parinfo=pi,$
                   bestnorm=bestnorm,dof=dof,perror=perror,/quiet)
                 
  endelse

  perrors = perror*sqrt(bestnorm/dof)
  
  r = {res:res,perrors:perrors,yMids:yMids,f:f,fargs:fargs}
  
  return,r
end

; patSpeedRadialMpfit(): radial TW using MPFIT model

function patSpeedRadialMpfit, simName, m, bar=bar

  basePath    = '/n/home07/dnelson/spirals/lambdacrit/sims/'
  filePath    = basePath + simName + "/"
  workingPath = filePath + "/output3/"

  units = getUnits()

  snapPath = filePath + "output/snap_"
  
  ;load snapshot
  h  = loadSnapshot(snapPath,m,s_pos,s_vel,s_id,c_pos,c_vel,c_id)
  h2 = loadSimParams(filePath)

  maxY = 85.0 ;km/s/kpc for plotting
  minY = 35.0

  ; fiducial parameters
  if keyword_set(bar) then begin
    y_r_max = 20.0 ;kpc, bar
  endif else begin
    y_r_max = 30.0 ;kpc, spiral
  endelse
  
  yWidth  = 0.5
  xBins   = 500
  SA      = 0.0 ;slice angle
  
  ; model config
  break_radius = 3.6 ; 1 scale length (bar model only)
  sgN = 15 ;savgol filter width
  
  ; run single
  r = doMpfit(y_r_max,yWidth,break_radius,xBins,sgN,h2,s_pos,s_vel,SA=SA,bar=bar)

  ; visualize raw fit
  PS_Start, FILENAME=workingPath+"mpfit.rawfit."+str(m)+".eps", /nomatch, /quiet, font=0, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=7.5, ys=5.0, /inches, tt_font='Times';3/2  
    
    fsc_plot,[0],[0],/nodata,xrange=[min(r.yMids)*1.1,max(r.yMids)*1.1],$
             yrange=[min(r.f)/1e5*1.1,max(r.f)/1e5*1.1],$
             xtitle="",ytitle="f / 1e5 (RHS)",/xs,/ys,position=[0.15,0.35,0.95,0.95],$
             xtickname=replicate(' ',10)
    
    fsc_plot,r.yMids,r.f/1e5,psym=16,symsize=0.6,/overplot
    
    if keyword_set(bar) then begin
      fsc_plot,r.yMids,mpModelBar(r.yMids,r.res,_EXTRA=r.fargs)/1e5,psym=15,symsize=0.6,$
               color=fsc_color('orange'),/overplot
      resids = (r.f-mpModelBar(r.yMids,r.res,_EXTRA=r.fargs))/1e5
    endif else begin
      fsc_plot,r.yMids,mpModelSpiralPowerlaw(r.yMids,r.res,_EXTRA=r.fargs)/1e5,psym=15,symsize=0.6,$
               color=fsc_color('orange'),/overplot
      resids = (r.f-mpModelSpiralPowerlaw(r.yMids,r.res,_EXTRA=r.fargs))/1e5
    endelse
    
    fsc_plot,[0],[0],ytitle="Resids",/nodata,xrange=[min(r.yMids)*1.1,max(r.yMids)*1.1],$
             yrange=[min(resids)*1.1,max(resids)*1.2],/xs,/ys,$
             xtitle="<y> [kpc]",position=[0.15,0.15,0.95,0.35],/noerase
             
    fsc_plot,r.yMids,resids,psym=16,symsize=0.6,/overplot
    fsc_plot,[min(r.yMids)*1.1,max(r.yMids)*1.1],[0,0],line=1,/overplot
  PS_End
  
  ; calc theory IC curves
  M_star = h2.f_d * h2.m200
  rPts = findgen(300)/300.0 * (y_r_max)/2.0 + 0.02
  
  omega2 = func_omega2_ExpHern(units.G, rPts, h2.m200, M_star, h2.h, h2.a)
  kappa2 = func_kappa2_ExpHern(units.G, rPts, h2.m200, M_star, h2.h, h2.a)
  
  ; visualize omega_p fit
  PS_Start, FILENAME=workingPath+"mpfit.omega_p."+str(m)+".eps", /nomatch, /quiet, font=0, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=7.5, ys=5.0, /inches, tt_font='Times';3/2  
    
    fsc_plot,[0],[0],/nodata,xrange=[0,max(rPts/h2.h)],yrange=[0,120],$
             xtitle="Radius / Disk Scale Length [kpc]",ytitle="Pattern Speed [km/s/kpc]",/xs,/ys
    
    if keyword_set(bar) then begin
      ; double constant with break model
      fsc_plot,[r.res[2]/2],[r.res[0]],psym=16,symsize=1.0,/overplot
      fsc_plot,[0.02,r.res[2]],[r.res[0],r.res[0]],line=0,thick=3.0,/overplot
      fsc_plot,[(y_r_max/2.0+r.res[2])/2],[r.res[1]],psym=16,symsize=1.0,/overplot
      fsc_plot,[r.res[2],y_r_max/2.0],[r.res[1],r.res[1]],line=0,thick=3.0,/overplot
      fsc_plot,[r.res[2],r.res[2]],[r.res[0],r.res[1]],line=1,thick=3.0,/overplot
    endif else begin
      ; spiral powerlaw model
      omegaPts = r.res[0] / (rPts/h2.h)^(r.res[1])
      fsc_plot,rPts/h2.h,omegaPts,line=0,thick=4.0,/overplot
    endelse
    
    ;overplot omega, omega+-kappa/2, omega+-kappa/4
    fsc_plot,rPts/h2.h,sqrt(omega2),color=fsc_color('orange'),/overplot
    fsc_plot,rPts/h2.h,sqrt(omega2)-sqrt(kappa2)/2.0,color=fsc_color('orange'),line=2,/overplot
    fsc_plot,rPts/h2.h,sqrt(omega2)+sqrt(kappa2)/2.0,color=fsc_color('orange'),line=2,/overplot
    fsc_plot,rPts/h2.h,sqrt(omega2)-sqrt(kappa2)/4.0,color=fsc_color('orange'),line=1,/overplot
    fsc_plot,rPts/h2.h,sqrt(omega2)+sqrt(kappa2)/4.0,color=fsc_color('orange'),line=1,/overplot  
  PS_End

  return,r
end

; sliceAngleStability(): check spread on variable slice angle

pro sliceAngleStability, simName, bar=bar
  
  basePath    = '/n/home07/dnelson/spirals/lambdacrit/sims/'
  filePath    = basePath + simName + "/"
  workingPath = filePath + "output3/"
  snapPath    = filePath + "output/snap_"
  
  units = getUnits()

  for m=0,20,1 do begin ;bar_a=12,20,1, bar_b=60,98,2, LC2=1,20,1
  
    saveFileName = workingPath+"mpfit.stability.SA."+str(m)+".sav"
  
    if not (file_test(saveFileName)) then begin  
  
      ;load snapshot
      h  = loadSnapshot(snapPath,m,s_pos,s_vel,s_id,c_pos,c_vel,c_id)
      h2 = loadSimParams(filePath)
    
      ; fiducial parameters
      if keyword_set(bar) then begin
        y_r_max = 20.0 ;kpc, bar
      endif else begin
        y_r_max = 30.0 ;kpc, spiral
      endelse
      
      yWidth  = 0.5
      xBins   = 500
      SA      = 0.0 ;slice angle, CCW from x-axis SA=0
      
      ; model config
      break_radius = 3.6 ; 1 scale length (bar model only)
      sgN = 15 ;savgol filter width
      
      ; slice angle stability (bar)
      param_min  = 0.0
      param_max  = 380.0
      param_step = 2.0
      param_num  = (param_max-param_min)/param_step + 1
      
      param1    = findgen(param_num)/(param_num-1) * (param_max-param_min) + param_min
      omegaP1   = fltarr(2,param_num)
    
      for i=0,param_num-1 do begin
        print,m,i
        r = doMpfit(y_r_max,yWidth,break_radius,xBins,sgN,h2,s_pos,s_vel,SA=param1[i],bar=bar)
        omegaP1[0,i] = r.res[0]
        omegaP1[1,i] = r.res[1]
      endfor 
    
      ; save/restore
      save,param1,omegaP1,simName,m,filename=saveFileName
    endif else begin
      restore,saveFileName
    endelse
    
    ; load barSpeed results
    r = barSpeed(simName)
  
    ; plot
    start_PS, workingPath+"mpfit.SA."+str(m)+".eps"
      fsc_plot,[0],[0],/nodata,yrange=[0.0,80.0],xrange=[min(param1),max(param1)],$
               xtitle="Slice Angle [deg]",ytitle="Pattern Speed [km/s/kpc]",/xs,/ys
      
      if keyword_set(bar) then begin
        ; double constant with break model
        fsc_plot,param1,omegaP1[0,*],psym=-16,symsize=0.4,color=fsc_color('orange'),/overplot
        fsc_plot,param1,omegaP1[1,*],psym=-16,symsize=0.4,color=fsc_color('blue'),/overplot
     
        ; overplot bar angle
        fsc_plot,[r.SA[m],r.SA[m]],[10.0,70.0],line=0,color=fsc_color('red'),thick=1.5,/overplot
        fsc_plot,[r.SA[m],r.SA[m]]-90.0,[10.0,70.0],line=0,color=fsc_color('red'),thick=1.0,/overplot
        fsc_plot,[r.SA[m],r.SA[m]]+90.0,[10.0,70.0],line=0,color=fsc_color('red'),thick=1.0,/overplot 
        fsc_plot,[r.SA[m],r.SA[m]]-45.0,[10.0,70.0],line=0,color=fsc_color('green'),thick=1.0,/overplot
        fsc_plot,[r.SA[m],r.SA[m]]+45.0,[10.0,70.0],line=0,color=fsc_color('green'),thick=1.0,/overplot        
        fsc_plot,[r.SA[m],r.SA[m]]+135.0,[10.0,70.0],line=0,color=fsc_color('green'),thick=1.0,/overplot
        
        ; overplots/legend
        fsc_plot,[min(param1),max(param1)],[45.0,45.0],line=1,/overplot
        fsc_text,max(param1)*0.93,70,"Inner (Bar)",color=fsc_color('orange'),alignment=1.0,charsize=1.5
        fsc_text,max(param1)*0.93,64,"Outer",color=fsc_color('blue'),alignment=1.00,charsize=1.5
      endif else begin
        ; powerlaw slope/amp model
        fsc_plot,param1,omegaP1[0,*],psym=-16,symsize=0.4,color=fsc_color('orange'),/overplot
        fsc_plot,param1,omegaP1[1,*]*100,psym=-16,symsize=0.4,color=fsc_color('blue'),/overplot
        
        ; overplot averages
        fsc_plot,[0.0,max(param1)],[mean(omegaP1[0,*]),mean(omegaP1[0,*])],line=1,/overplot
        fsc_plot,[0.0,max(param1)],[mean(omegaP1[1,*]),mean(omegaP1[1,*])]*100,line=1,/overplot
        
        ; legend
        fsc_text,max(param1)*0.93,10,"Amplitude (at h)",color=fsc_color('orange'),alignment=1.0,charsize=1.5
        fsc_text,max(param1)*0.93,4,"Slope x100",color=fsc_color('blue'),alignment=1.00,charsize=1.5
        
        ; principal axes
        fsc_plot,[90.0,90.0],[0.0,80.0],line=1,/overplot
        fsc_plot,[180.0,180.0],[0.0,80.0],line=1,/overplot
        fsc_plot,[360.0,360.0],[0.0,80.0],line=1,/overplot
      endelse
      
    end_PS

  endfor ;m

end

; mpfitStability(): run numerical parameter stability tests

pro mpfitStability, simName, m, bar=bar

  basePath    = '/n/home07/dnelson/spirals/lambdacrit/sims/'
  filePath    = basePath + simName + "/"
  workingPath = filePath + "/output3/"

  units = getUnits()

  snapPath = filePath + "output/snap_"
  
  ;load snapshot
  h  = loadSnapshot(snapPath,m,s_pos,s_vel,s_id,c_pos,c_vel,c_id)
  h2 = loadSimParams(filePath)

  maxY = 85.0 ;km/s/kpc for plotting
  minY = 35.0

  ; fiducial parameters
  if keyword_set(bar) then begin
    y_r_max = 20.0 ;kpc, bar
  endif else begin
    y_r_max = 30.0 ;kpc, LC
  endelse
  
  yWidth  = 0.5
  xBins   = 500
  SA      = 0.0 ;slice angle, CCW from x-axis SA=0
  
  ; model config
  break_radius = 3.6 ; 1 scale length (bar model only)
  sgN = 15 ;savgol filter width

  ; stability: break_radius (bar)
  ; -----------------------------
  if 0 then begin
  param_min  = 1.0
  param_max  = 6.0
  param_step = 0.05
  param_num  = (param_max-param_min)/param_step + 1
  
  param1    = findgen(param_num)/(param_num-1) * (param_max-param_min) + param_min
  omegaP1   = fltarr(2,param_num)

  for i=0,param_num-1 do begin
    print,i
    r = doMpfit(y_r_max,yWidth,param1[i],xBins,sgN,h2,s_pos,s_vel,SA=SA,bar=bar)
    omegaP1[0,i] = r.res[0]
    omegaP1[1,i] = r.res[1]
  endfor 
  
  ; plot
  start_PS, workingPath+"mpfit.break_radius."+str(m)+".eps"
    fsc_plot,[0],[0],/nodata,yrange=[0,maxY],xrange=[min(param1),max(param1)],$
             xtitle="Break Radius [kpc]",ytitle="Pattern Speed [km/s/kpc]",/xs,/ys
    
    ; double constant with break model
    fsc_plot,param1,omegaP1[0,*],psym=-16,symsize=0.4,color=fsc_color('orange'),/overplot
    fsc_plot,param1,omegaP1[1,*],psym=-16,symsize=0.4,color=fsc_color('blue'),/overplot
 
    ; overplots/legend
    fsc_plot,[min(param1),max(param1)],[45.0,45.0],line=1,/overplot
    fsc_text,max(param1)*0.95,10,"Inner (Bar)",color=fsc_color('orange'),alignment=1.0,charsize=1.5
    fsc_text,max(param1)*0.95,6,"Outer",color=fsc_color('blue'),alignment=1.00,charsize=1.5
  end_PS
  endif ;0
  
  ; stability: savgolN (spiral)
  ; ---------------------------
  if 0 then begin
  param_min  = 3
  param_max  = 25
  param_step = 1
  param_num  = (param_max-param_min)/param_step + 1
  
  param1    = findgen(param_num)/(param_num-1) * (param_max-param_min) + param_min
  omegaP1   = fltarr(2,param_num)

  for i=0,param_num-1 do begin
    print,i
    r = doMpfit(y_r_max,yWidth,break_radius,xBins,param1[i],h2,s_pos,s_vel,SA=SA,bar=bar)
    omegaP1[0,i] = r.res[0]
    omegaP1[1,i] = r.res[1]
  endfor 
  
  ; plot
  start_PS, workingPath+"mpfit.sgN."+str(m)+".eps"
    fsc_plot,[0],[0],/nodata,yrange=[0,maxY],xrange=[min(param1),max(param1)],$
             xtitle="Savgol Filter Width [#]",ytitle="Amplitude [km/s/kpc] / 100x Slope",/xs,/ys
    
    ; double constant with break model
    fsc_plot,param1,omegaP1[0,*],psym=-16,symsize=0.4,color=fsc_color('orange'),/overplot
    fsc_plot,param1,omegaP1[1,*]*100,psym=-16,symsize=0.4,color=fsc_color('blue'),/overplot
 
    ; overplots/legend
    fsc_plot,[min(param1),max(param1)],[45.0,45.0],line=1,/overplot
    fsc_text,max(param1)*0.95,10,"Amplitude",color=fsc_color('orange'),alignment=1.0,charsize=1.5
    fsc_text,max(param1)*0.95,6,"Slope",color=fsc_color('blue'),alignment=1.00,charsize=1.5
  end_PS

  ; stability - y_r_max
  ; -------------------
  param_min  = 8.0
  param_max  = 40.0
  param_step = 0.4
  param_num  = (param_max-param_min)/param_step + 1
  
  param2    = findgen(param_num)/(param_num-1) * (param_max-param_min) + param_min
  omegaP2   = fltarr(2,param_num)

  for i=0,param_num-1 do begin
    print,i
    r = doMpfit(param2[i],yWidth,break_radius,xBins,sgN,h2,s_pos,s_vel,SA=SA,bar=bar)
    omegaP2[0,i] = r.res[0]
    omegaP2[1,i] = r.res[1]
  endfor 
  
  ; plot
  start_PS, workingPath+"mpfit.y_r_max."+str(m)+".eps"
    fsc_plot,[0],[0],/nodata,yrange=[0,maxY],xrange=[min(param2),max(param2)],$
             xtitle="Maximum Y Extent [kpc]",ytitle="Amplitude [km/s/kpc] / 100x Slope",/xs,/ys
    
    ; double constant with break model
    fsc_plot,param2,omegaP2[0,*],psym=-16,symsize=0.4,color=fsc_color('orange'),/overplot
    fsc_plot,param2,omegaP2[1,*]*100,psym=-16,symsize=0.4,color=fsc_color('blue'),/overplot
 
    ; overplots/legend
    fsc_plot,[min(param2),max(param2)],[45.0,45.0],line=1,/overplot
    fsc_text,max(param2)*0.95,10,"Amplitude",color=fsc_color('orange'),alignment=1.0,charsize=1.5
    fsc_text,max(param2)*0.95,6,"Slope",color=fsc_color('blue'),alignment=1.00,charsize=1.5
  end_PS
  
  ; stability - yWidth
  ; -------------------
  param_min  = 0.05 ;50pc
  param_max  = 1.0  ;1kpc
  param_step = 0.0125
  param_num  = (param_max-param_min)/param_step + 1
  
  param3    = findgen(param_num)/(param_num-1) * (param_max-param_min) + param_min
  omegaP3   = fltarr(2,param_num)

  for i=0,param_num-1 do begin
    print,i
    r = doMpfit(y_r_max,param3[i],break_radius,xBins,sgN,h2,s_pos,s_vel,SA=SA,bar=bar)
    omegaP3[0,i] = r.res[0]
    omegaP3[1,i] = r.res[1]
  endfor 
  
  ; plot
  start_PS, workingPath+"mpfit.yWidth."+str(m)+".eps"
    fsc_plot,[0],[0],/nodata,yrange=[0,maxY],xrange=[min(param3),max(param3)],$
             xtitle="Y Slice Width [kpc]",ytitle="Amplitude [km/s/kpc] / 100x Slope",/xs,/ys
    
    ; double constant with break model
    fsc_plot,param3,omegaP3[0,*],psym=-16,symsize=0.4,color=fsc_color('orange'),/overplot
    fsc_plot,param3,omegaP3[1,*]*100,psym=-16,symsize=0.4,color=fsc_color('blue'),/overplot
 
    ; overplots/legend
    fsc_plot,[min(param3),max(param3)],[45.0,45.0],line=1,/overplot
    fsc_text,max(param3)*0.95,10,"Amplitude",color=fsc_color('orange'),alignment=1.0,charsize=1.5
    fsc_text,max(param3)*0.95,6,"Slope",color=fsc_color('blue'),alignment=1.00,charsize=1.5
  end_PS
  
  ; stability - xBins
  ; -------------------
  param_min  = 50 ;50pc
  param_max  = 1000  ;1kpc
  param_step = 10
  param_num  = (param_max-param_min)/param_step + 1
  
  param4    = findgen(param_num)/(param_num-1) * (param_max-param_min) + param_min
  omegaP4   = fltarr(2,param_num)

  for i=0,param_num-1 do begin
    print,i
    r = doMpfit(y_r_max,yWidth,break_radius,param4[i],sgN,h2,s_pos,s_vel,SA=SA,bar=bar)
    omegaP4[0,i] = r.res[0]
    omegaP4[1,i] = r.res[1]
  endfor 
  
  ; plot
  start_PS, workingPath+"mpfit.xBins."+str(m)+".eps"
    fsc_plot,[0],[0],/nodata,yrange=[0,maxY],xrange=[min(param4),max(param4)],$
             xtitle="xBins (Line Integral Binning) [#]",ytitle="Amplitude [km/s/kpc] / 100x Slope",/xs,/ys
    
    ; double constant with break model
    fsc_plot,param4,omegaP4[0,*],psym=-16,symsize=0.4,color=fsc_color('orange'),/overplot
    fsc_plot,param4,omegaP4[1,*]*100,psym=-16,symsize=0.4,color=fsc_color('blue'),/overplot
 
    ; overplots/legend
    fsc_plot,[min(param4),max(param4)],[45.0,45.0],line=1,/overplot
    fsc_text,max(param4)*0.95,10,"Amplitude",color=fsc_color('orange'),alignment=1.0,charsize=1.5
    fsc_text,max(param4)*0.95,6,"Slope",color=fsc_color('blue'),alignment=1.00,charsize=1.5
  end_PS
  
  ; save/restore
  endif ;0
  ;save,param1,param2,param3,param4,omegaP1,omegaP2,omegaP3,omegaP4,simName,m,$
  ;     filename=workingPath+"mpfit.stability."+str(m)+".sav"
  restore,workingPath+"mpfit.stability."+str(m)+".sav"

  ; fix bad pt
  w = where(OmegaP1[1,*] lt 0,count)
  if (count ne 0) then begin
    omegaP1[*,w] = omegaP1[*,w-1]
  endif

  ; plot composite
  start_PS, workingPath+"mpfit.stability."+str(m)+".eps"
    xm = !x.margin
    ym = !y.margin
    !x.margin = [0,0]
    !y.margin = [0,0]
      
    ; break radius / sgN (upper left)
    fsc_plot,[0],[0],/nodata,yrange=[minY,maxY],xrange=[min(param1),max(param1)*0.99],$
             xtitle="",ytitle="",xs=9,/ys,position=[0.15,0.5,0.5,0.85],xtickname=replicate(' ',10)
    fsc_plot,param1,omegaP1[0,*],psym=-16,symsize=0.4,color=fsc_color('orange'),/overplot
    fsc_plot,param1,omegaP1[1,*]*100,psym=-16,symsize=0.4,color=fsc_color('blue'),/overplot 
    ;fsc_plot,[min(param1),max(param1)],[45.0,45.0],line=1,/overplot 
    fsc_axis,0.0,maxY,0.0,/xaxis,/xs,xrange=[min(param1),max(param1)*0.99]
    fsc_text,(0.5+0.15)/2.0,0.93,"SG Filter Size [#]",alignment=0.5,/normal

    ; y_r_max (upper right)
    fsc_plot,[0],[0],/nodata,yrange=[minY,maxY],xrange=[min(param2),max(param2)],$
             xtitle="",ytitle="",xs=9,ys=9,position=[0.5,0.5,0.85,0.85],/noerase,$
             xtickname=replicate(' ',10),ytickname=replicate(' ',10)
    fsc_plot,param2,omegaP2[0,*],psym=-16,symsize=0.4,color=fsc_color('orange'),/overplot
    fsc_plot,param2,omegaP2[1,*]*100,psym=-16,symsize=0.4,color=fsc_color('blue'),/overplot
    ;fsc_plot,[min(param2),max(param2)],[45.0,45.0],line=1,/overplot 
    fsc_axis,0.0,maxY,0.0,/xaxis,/xs,xrange=[min(param2),max(param2)]
    fsc_axis,max(param2),0.0,0.0,/yaxis,/ys,yrange=[minY,maxY]
    fsc_text,(0.85+0.5)/2.0,0.93,"Maximum Y [kpc]",alignment=0.5,/normal
    
    ; yWidth (lower left)
    fsc_plot,[0],[0],/nodata,yrange=[minY,maxY],xrange=[min(param3),max(param3)*0.99],$
             xtitle="yWidth [kpc]",ytitle="",/xs,/ys,position=[0.15,0.15,0.5,0.5],/noerase
    fsc_plot,param3,omegaP3[0,*],psym=-16,symsize=0.4,color=fsc_color('orange'),/overplot
    fsc_plot,param3,omegaP3[1,*]*100,psym=-16,symsize=0.4,color=fsc_color('blue'),/overplot 
    ;fsc_plot,[min(param3),max(param3)],[45.0,45.0],line=1,/overplot 
    
    ; xBins (lower right)
    fsc_plot,[0],[0],/nodata,yrange=[minY,maxY],xrange=[min(param4),max(param4)],$
             xtitle="xBins [#]",ytitle="",/xs,ys=9,position=[0.5,0.15,0.85,0.5],/noerase,ytickname=replicate(' ',10)
    fsc_plot,param4,omegaP4[0,*],psym=-16,symsize=0.4,color=fsc_color('orange'),/overplot
    fsc_plot,param4,omegaP4[1,*]*100,psym=-16,symsize=0.4,color=fsc_color('blue'),/overplot 
    ;fsc_plot,[min(param4),max(param4)],[45.0,45.0],line=1,/overplot 
    fsc_axis,max(param4),0.0,0.0,/yaxis,/ys,yrange=[minY,maxY]
    
    ; legend
    fsc_text,max(param4)*0.95,58,"Amplitude",color=fsc_color('orange'),alignment=1.0,charsize=1.5
    fsc_text,max(param4)*0.95,50,"Slope",color=fsc_color('blue'),alignment=1.00,charsize=1.5
    
    ; y axes
    fsc_text,0.07,0.5,"Amplitude [km/s/kpc] / 100x Slope",alignment=0.5,orientation=90.0,/normal
    
    !x.margin = xm
    !y.margin = ym
  end_PS
  
end

; patSpeedRadial(): RTW method unregularized backsubstitution (UNUSED)

function patSpeedRadial;, filePath, workingPath, m, time, extent

  m = 14

  basePath = '/n/home07/dnelson/spirals/lambdacrit/sims/'
  simName  = 'bar_1m_a'

  filePath    = basePath + simName + "/"
  workingPath = filePath + "/output3/"

  snapPath = filePath + "output/snap_"
  
  h2 = loadSimParams(filePath)
  
  ; y/r config
  y_r_max = 25.0 ;kpc
  N       = 71
  yWidth  = 0.1
  
  tiltang = 45 ;?
  
  deltaR  = y_r_max / N  
  
  yStarts = findgen(N)/(N-1) * y_r_max + deltaR/2.0 ;y_j
  midR    = findgen(N)/(N-1) * y_r_max + deltaR/2.0 ;r_i
  
  ; x-int config
  xBins   = 500
  xMinMax = [-10.0,10.0]*h2.h ; whole distribution
  
  dx = (xMinMax[1] - xMinMax[0]) / xBins
  xMid = findgen(xBins)/(xBins-1) * (xMinMax[1] - xMinMax[0]) + xMinMax[0] + dx/2.0
  
  ; arrays
  K       = fltarr(N,N)
  v_y     = fltarr(N)
  omega_p = fltarr(N)

  ;load snapshot
  h = loadSnapshot(snapPath,m,s_pos,s_vel,s_id,c_pos,c_vel,c_id)
  time = h.time  
    
  x = reform(s_pos[0,*])
  y = reform(s_pos[1,*])
  
  vx = reform(s_vel[0,*])
  vy = reform(s_vel[1,*])
  
  ; calculate kernel
  for i=0,N-1 do begin
    for j=0,N-1 do begin
      ;upper triangle only
      if (j gt i) then continue
      
      ;coordinates (r_i, y_j -> K_ji)
      r = midR[i]
      yp = yStarts[j] * cos(tiltang*!dtor)
      
      xp = sqrt(r^2 - yp^2)
      
      ;surface density
      w = where(x gt xp-dx/2.0 and x le xp+dx/2.0 and $
                y gt yp and y le yp+yWidth, count1)
      w = where(x gt -xp-dx/2.0 and x le -xp+dx/2.0 and $
                y gt yp and y le yp+yWidth, count2)               
                
      if (count1 ne 0) then $
        sigma_pos = count1 / (dx*yWidth)
      if (count2 ne 0) then $
        sigma_neg = count2 / (dx*yWidth)
      
      ;kernel
      K[i,j] = (sigma_pos - sigma_neg) * r

    endfor
  endfor

  ; calculate f(y)
  for i=0,N-1 do begin
    ;select y subset
    yMin = yStarts[i] * cos(tiltang*!dtor)
    yMax = yMin + yWidth
    w1 = where(y ge yMin and y lt yMax, count1)
    
    x1  = x[w1]
    vy1 = vy[w1] ;* sin(tiltang*!dtor)
    
    vAvgTemp = fltarr(xBins)
    
    ;for a given y, calculate <v_y> along dx
    for j=0,xBins-1 do begin
      xMin = xMid[j] - dx/2.0
      xMax = xMid[j] + dx/2.0
      
      ;find stars in bin
      w = where(x1 ge xMin and x1 lt xMax, count)
      
      if (count ne 0) then begin
        surfaceDens = count / dx / yWidth
        vyMean      = mean(vy1[w])
        
        vAvgTemp[j] = dx * surfaceDens * vyMean
      endif
    endfor
  
    v_y[i] = total(vAvgTemp)
  endfor
  
  ; solution via back-substitution
  omega_p[N-1] = v_y[N-1] / K[N-1,N-1]
  
  for i=N-2,0,-1 do begin
    omega_p[i] = v_y[i]
    omega_p[i] -= total(K[i+1:*,i+1]*omega_p[i+1:*])
    omega_p[i] /= K[i,i]
    print,i,total(K[i+1:*,i+1]*omega_p[i+1:*]),total(K[i+1,i+1:*]*omega_p[i+1:*])
  endfor
  
  ; calculate a chi-squared
  chi2 = 0
  for i=0,N-1 do begin
    stddev = 1.0
    chi2 += total( (K[i,*]*omega_p - v_y)^2.0 / stddev^2.0 )
  endfor
  print,'chi2 ',chi2
  
  ; plot
  start_PS, workingPath+"patSpeedRad_"+str(m)+"_"+str(N)+".eps"
  
    fsc_plot,[0],[0],psym=4, $ ;1=cross, 3=dot, 4=diamond
         xtitle="radius [kpc]",ytitle=""+textoidl("\Omega_p")+" [km/s/kpc]",/nodata, $
         xrange=[0,max(midR)*1.1],yrange=[min(omega_p)*1.1,max(omega_p)*1.1],/xs,/ys
         
    fsc_plot,midR,omega_p,psym=16,/overplot ;filled circle
  
    ;draw some bin widths
    for i=0,N-1 do begin
    
      xMin = midR[i]-deltaR/2.0
      xMax = midR[i]+deltaR/2.0
      fsc_plot,[xMin,xMax],[omega_p[i],omega_p[i]],line=0,/overplot
    endfor
  
  end_PS
  
  print,'K ',K
  print,'v_y ',v_y
  print,'omega_p ',omega_p
  stop
end

; patSpeedConst(): calculate constant pattern speed and radially bin

function patSpeedConst, filePath, workingPath, m, extent, nRadBins

  units = getUnits()

  snapPath = filePath + "output/snap_"
  
  omega_p = fltarr(n_elements(extent))
  
  ;load snapshot
  h  = loadSnapshot(snapPath,m,s_pos,s_vel,s_id,c_pos,c_vel,c_id)
  h2 = loadSimParams(filePath)
  time = h.time

  x = reform(s_pos[0,*])
  y = reform(s_pos[1,*])
  
  vx = reform(s_vel[0,*])
  vy = reform(s_vel[1,*])  
  
  ;loop over extents
  for k=0,n_elements(extent)-1 do begin
  
  ;config(x)
  xBins   = 500
  xMinMax = [-10.0,10.0]*h2.h ; whole distribution
  
  dx = (xMinMax[1] - xMinMax[0]) / xBins
  xMid = findgen(xBins)/(xBins-1) * (xMinMax[1] - xMinMax[0]) + xMinMax[0] + dx/2.0

  ;config(y)
  yWidth  = 0.5
  yMax = extent[k]
  
  yRes = floor(yMax / yWidth)
  yStarts = findgen(yRes)/(yRes-1) * yMax - yMax/2
  
  ;arrays
  yBins = n_elements(yStarts)
  vAvg = fltarr(yBins)
  xAvg = fltarr(yBins)
  
  ;loop over annuli
  for i=0,yBins-1 do begin

    yMin = yStarts[i]
    yMax = yMin + yWidth
    
    surfaceDens = fltarr(xBins)
    vyMean      = fltarr(xBins)
    
    ;select y subset
    w1 = where(y ge yMin and y lt yMax, count1)
    
    x1  = x[w1]
    vy1 = vy[w1]
    
    ;for a given annulus, calculate pattern speed by integrating along dtheta
    for j=0,xBins-1 do begin
      xMin = xMid[j] - dx/2.0
      xMax = xMid[j] + dx/2.0
      
      ;find stars in bin
      w = where(x1 ge xMin and x1 lt xMax, count)
      
      if (count ne 0) then begin
        surfaceDens[j] = count / dx / yWidth
        vyMean[j]      = mean(vy1[w])
      endif
    endfor
    
    ; INT_TAB method
    f_xAvg = xMid * surfaceDens
    f_vAvg = vyMean * surfaceDens
    
    xAvgt = int_tabulated(xMid,f_xAvg)
    vAvgt = int_tabulated(xMid,f_vAvg)
    norm = int_tabulated(xMid,surfaceDens)
    
    xAvg[i] = total(xAvgt) / total(norm)
    vAvg[i] = total(vAvgt) / total(norm)

  endfor
  
  ; global fiting
  res = linfit(xAvg,vAvg, yfit=yfit)
  print,'fit',time,res

  omega_p[k] = res[1]
  
  ; crude radial fitting
  radBinWidth = extent[k] / nRadBins
  
  avgOmegaP = fltarr(nRadBins)
  avgYPts   = fltarr(nRadBins)
  
  for i=0,nRadBins-1 do begin
    yMin = -extent[k]/2.0 + radBinWidth*i
    yMax = -extent[k]/2.0 + radBinWidth*(i+1)
    w = where(yStarts ge yMin and yStarts lt yMax,count)
    if (count ne 0) then begin
      res = linfit(xAvg[w],vAvg[w])
      avgOmegaP[i] = res[1]
      avgYPts[i]   = round(mean([yMin,yMax])*100)/100.0
    endif
  endfor

  ; IC-theory circular velocity lines
  M_star = h2.f_d * h2.m200
  rPts = findgen(300)/300.0 * extent[k] + 0.1
  
  omega2 = func_omega2_ExpHern(units.G, rPts, h2.m200, M_star, h2.h, h2.a)
  kappa2 = func_kappa2_ExpHern(units.G, rPts, h2.m200, M_star, h2.h, h2.a)

  ; plot fit
  start_PS, workingPath+"patSpeedConst_"+str(m)+"_"+str(k)+".eps"

    fsc_plot,xAvg,vAvg,psym=4, $ ;1=cross, 3=dot, 4=diamond
         xtitle="<x> [kpc]",ytitle="<v> [km/s]",/nodata
         
    ; already sorted by yStarts
    fsc_plot,xAvg,vAvg,line=1,/overplot
    for i=0,n_elements(xAvg)-1,1 do begin
      fsc_text,xAvg[i],vAvg[i],string(round(yStarts[i]),format='(i3)'),alignment=0.5,$
               charsize=!p.charsize*0.5,color=fsc_color('blue')
    endfor

    fsc_plot,xAvg,yfit,line=0,/overplot
    fsc_plot,[0,0],minmax(vAvg),line=1,/overplot
    fsc_plot,minmax(xAvg),[0,0],line=1,/overplot
    fsc_text,0.80,0.25,textoidl("\Omega_p=")+str(omega_p[k])+" km/s/kpc",$
             /normal,charsize=!p.charsize*0.75,alignment=0.5
    fsc_text,0.80,0.20,"extent="+str(extent[k])+" kpc",$
             /normal,charsize=!p.charsize*0.75,alignment=0.5

  end_PS
  
  ; plot crude radial
  start_PS, workingPath+"patSpeedConst_rad_"+str(m)+"_"+str(k)+".eps"
  
    fsc_plot,[0],[0],psym=4, $ ;1=cross, 3=dot, 4=diamond
         xtitle="<y> [kpc]",ytitle="<"+textoidl("\Omega_p")+"> [km/s/kpc]",/nodata, $
         xrange=[0,max(avgYPts)*1.1],yrange=[0,max(avgOmegaP)*1.1],/xs,/ys
         
    w1 = where(avgYPts gt 0)
    fsc_plot,avgYPts[w1],avgOmegaP[w1],psym=15,/overplot ;filled square
    w2 = where(avgYPts lt 0)
    fsc_plot,-avgYPts[w2],avgOmegaP[w2],psym=6,/overplot ;open square
  
    ;draw some bin widths
    for i=0,nRadBins-1 do begin
      if (avgYPts[i] le 0) then continue
      
      xMin = avgYPts[i]-radBinWidth/2.0
      xMax = avgYPts[i]+radBinWidth/2.0
      fsc_plot,[xMin,xMax],[avgOmegaP[i],avgOmegaP[i]],line=0,/overplot
    endfor
  
    ;overplot omega, omega-kappa/2, omega+kappa/2
    fsc_plot,rPts,sqrt(omega2),color=fsc_color('orange'),/overplot
    fsc_plot,rPts,sqrt(omega2)-sqrt(kappa2)/2.0,color=fsc_color('orange'),line=2,/overplot
    fsc_plot,rPts,sqrt(omega2)+sqrt(kappa2)/2.0,color=fsc_color('orange'),line=2,/overplot
    fsc_plot,rPts,sqrt(omega2)-sqrt(kappa2)/4.0,color=fsc_color('orange'),line=1,/overplot
    fsc_plot,rPts,sqrt(omega2)+sqrt(kappa2)/4.0,color=fsc_color('orange'),line=1,/overplot     
  
  end_PS
  
  endfor ;k
  
  r = {omega_p:omega_p,avgOmegaP:avgOmegaP,avgYPts:avgYPts,time:time}

  return,r
end

; patSpeedSeriesBar(): run many

pro patSpeedSeriesBar, simName

  ;config
  basePath    = '/n/home07/dnelson/spirals/lambdacrit/sims/'
  filePath    = basePath + simName + "/"
  workingPath = filePath + "/output3/"
  
  saveFileName = filePath + 'pSS.72.sav'  
  
  nSnaps = n_elements(file_search(filePath+"output/snap_*"))
  
  h2 = loadSimParams(filePath)  
  
  extent   = [7.2]
  nRadBins = 4
  sizeRes  = 3
  
  start = 60 ;bar_a=12, bar_b=60 (t > 600 Myr)
  
  colors = ['pink','skyblue','forest green','purple','cyan']
  
  if not (file_test(saveFileName)) then begin
  
    omega_p_const      = fltarr(nSnaps)
    omega_p_const_rad  = fltarr(nRadBins,nSnaps)
    omega_p_const_ypts = fltarr(nRadBins,nSnaps)
    omega_p_rad_mpfit  = fltarr(sizeRes,nSnaps)
    times              = fltarr(nSnaps)
    
    for m=0,nSnaps-1,1 do begin 
      r = patSpeedConst(filePath, workingPath, m, extent, nRadBins) ;extent ~ bar length
      
      omega_p_const[m]        = r.omega_p   ;extents
      times[m]                = r.time
      
      r = patSpeedConst(filePath, workingPath, m, extent*2.0, nRadBins) ;extent ~ 2x bar length
      
      omega_p_const_rad[*,m]  = r.avgOmegaP ;radially binned
      omega_p_const_ypts[*,m] = r.avgYPts   ;radial midbins
      
      r = patSpeedRadialMpfit(simName, m, /bar)
      
      omega_p_rad_mpfit[*,m]  = r.res       ;omega_bar,omega_outer,break_radius
    endfor
  
    ; save/restore
    save,omega_p_const,omega_p_const_rad,omega_p_const_ypts,omega_p_rad_mpfit,$
         times,filename=saveFileName
  endif else begin
    restore,saveFileName
  endelse

  ; reform radially binned
  const_rad = fltarr(nRadBins/2,nSnaps)
  ypts      = fltarr(nRadBins/2,nSnaps)
  
  for m=0,nSnaps-1,1 do begin
    for i=0,nRadBins/2-1 do begin
      val1 = omega_p_const_rad[i,m]
      val2 = omega_p_const_rad[nRadBins-i-1,m]
      ypts[i,m] = omega_p_const_ypts[nRadBins-i-1,m]
      
      if (val1 gt 0 and val2 gt 0) then $
        const_rad[i,m] = (val1+val2)/2.0
      if (val1 gt 0 and val2 lt 0) then $
        const_rad[i,m] = val1
      if (val1 lt 0 and val2 gt 0) then $
        const_rad[i,m] = val2
    endfor
  endfor
  
  ; load barSpeed results
  barR = barSpeed(simName)
  
  ; print mean+stddev
  print,'global constant ', mean(omega_p_const[start:*]),stddev(omega_p_const[start:*])
  print,'binned tw (4-8) ', mean(const_rad[0,start:*]),stddev(const_rad[0,start:*])
  print,'binned tw (0-4) ', mean(const_rad[1,start:*]),stddev(const_rad[1,start:*])
  print,'rtw inner (<3.6) ',mean(omega_p_rad_mpfit[0,start:*]),stddev(omega_p_rad_mpfit[0,start:*])
  print,'rtw outer (>3.6) ',mean(omega_p_rad_mpfit[1,start:*]),stddev(omega_p_rad_mpfit[1,start:*])
  
  ; plot
  start_PS, workingPath+"pSS_OmegaP_histos_72.eps"

    !p.multi = [0,2,3]

    histoplot,omega_p_const[start:*],title="Global Constant",/fillpolygon,nbins=8
    fsc_plot,[mean(omega_p_const[start:*]),mean(omega_p_const[start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[median(omega_p_const[start:*]),median(omega_p_const[start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[2]),/overplot
             
    histoplot,const_rad[0,start:*],title="Binned TW (3.6-7.2 kpc)",/fillpolygon,nbins=8
    fsc_plot,[mean(const_rad[0,start:*]),mean(const_rad[0,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[median(const_rad[0,start:*]),median(const_rad[0,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[2]),/overplot
    histoplot,const_rad[1,start:*],title="Binned TW (0-3.6 kpc)",/fillpolygon,nbins=8
    fsc_plot,[mean(const_rad[1,start:*]),mean(const_rad[1,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[median(const_rad[1,start:*]),median(const_rad[1,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[2]),/overplot
             
    histoplot,omega_p_rad_mpfit[0,start:*],title="RTW Inner (<3.6)",/fillpolygon,nbins=8
    fsc_plot,[mean(omega_p_rad_mpfit[0,start:*]),mean(omega_p_rad_mpfit[0,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[median(omega_p_rad_mpfit[0,start:*]),median(omega_p_rad_mpfit[0,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[2]),/overplot
    histoplot,omega_p_rad_mpfit[1,start:*],title="RTW Outer (>3.6)",/fillpolygon,nbins=8
    fsc_plot,[mean(omega_p_rad_mpfit[1,start:*]),mean(omega_p_rad_mpfit[1,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[median(omega_p_rad_mpfit[1,start:*]),median(omega_p_rad_mpfit[1,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[2]),/overplot
    
    fsc_text,0.8,0.15,"Mean",color=fsc_color(colors[1]),charsize=1.2,/normal
    fsc_text,0.8,0.10,"Median",color=fsc_color(colors[2]),charsize=1.2,/normal
    
    !p.multi = [0,0,0]

  end_PS

  start_PS, workingPath+"pSS_OmegaP_wtime_72.eps"

    !p.multi = [0,2,3]
  
    xMinMax = [0.5,2.0]

    fsc_plot,times[start:*],omega_p_const[start:*],title="Global Constant",xrange=xMinMax,/xs,/ys
    fsc_plot,[0.5,1.0],[mean(omega_p_const[start:*]),mean(omega_p_const[start:*])],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[0.5,1.0],[median(omega_p_const[start:*]),median(omega_p_const[start:*])],$
             line=1,color=fsc_color(colors[2]),/overplot  
                        
    fsc_plot,times[start:*],const_rad[0,start:*],title="Binned TW (3.6-7.2 kpc)",xrange=xMinMax,/xs,/ys
    fsc_plot,[0.5,1.0],[mean(const_rad[0,start:*]),mean(const_rad[0,start:*])],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[0.5,1.0],[median(const_rad[0,start:*]),median(const_rad[0,start:*])],$
             line=1,color=fsc_color(colors[2]),/overplot 
    fsc_plot,times[start:*],const_rad[1,start:*],title="Binned TW (0-3.6 kpc)",xrange=xMinMax,/xs,/ys
    fsc_plot,[0.5,1.0],[mean(const_rad[1,start:*]),mean(const_rad[1,start:*])],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[0.5,1.0],[median(const_rad[1,start:*]),median(const_rad[1,start:*])],$
             line=1,color=fsc_color(colors[1]),/overplot

    fsc_plot,times[start:*],omega_p_rad_mpfit[0,start:*],title="RTW Inner (<3.6)",xrange=xMinMax,/xs,/ys,$
             xtitle="Time [Gyr]"
    fsc_plot,[0.5,1.0],[mean(omega_p_rad_mpfit[0,start:*]),mean(omega_p_rad_mpfit[0,start:*])],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[0.5,1.0],[median(omega_p_rad_mpfit[0,start:*]),median(omega_p_rad_mpfit[0,start:*])],$
             line=1,color=fsc_color(colors[2]),/overplot
    fsc_plot,times[start:*],omega_p_rad_mpfit[1,start:*],title="RTW Outer (>3.6)",xrange=xMinMax,/xs,/ys
    fsc_plot,[0.5,1.0],[mean(omega_p_rad_mpfit[1,start:*]),mean(omega_p_rad_mpfit[1,start:*])],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[0.5,1.0],[median(omega_p_rad_mpfit[1,start:*]),median(omega_p_rad_mpfit[1,start:*])],$
             line=1,color=fsc_color(colors[2]),/overplot
                 
    fsc_text,0.8,0.15,"Mean",color=fsc_color(colors[1]),charsize=1.2,/normal
    fsc_text,0.8,0.10,"Median",color=fsc_color(colors[2]),charsize=1.2,/normal
    
    !p.multi = [0,0,0]

  end_PS
  
  ; correlations
  start_PS, workingPath+"pSS_OmegaP_corrs.eps"

    !p.multi = [0,2,3]

    fsc_plot,barR.SA[start:*],omega_p_rad_mpfit[0,start:*],psym=16,symsize=0.5,$
             xrange=[-10,190],/xs,ytitle="RTW Inner (<3.6)"
    
    fsc_plot,[0,0],[0,100],line=1,thick=1.0,/overplot  
    fsc_plot,[90,90],[0,100],line=1,thick=1.0,/overplot
    fsc_plot,[180,180],[0,100],line=1,thick=1.0,/overplot
    fsc_plot,[-10,190],[45,45],line=1,thick=1.0,/overplot
    
    fsc_plot,barR.SA[start:*],omega_p_rad_mpfit[1,start:*],psym=16,symsize=0.5,$
             xrange=[-10,190],/xs,ytitle="RTW Outer (>3.6)"
    
    fsc_plot,[0,0],[0,100],line=1,thick=1.0,/overplot  
    fsc_plot,[90,90],[0,100],line=1,thick=1.0,/overplot
    fsc_plot,[180,180],[0,100],line=1,thick=1.0,/overplot
    fsc_plot,[-10,190],[45,45],line=1,thick=1.0,/overplot
    
    fsc_plot,barR.SA[start:*],omega_p_const[start:*],psym=16,symsize=0.5,$
             xrange=[-10,190],/xs,ytitle="Global Constant"
    
    fsc_plot,[0,0],[0,100],line=1,thick=1.0,/overplot  
    fsc_plot,[90,90],[0,100],line=1,thick=1.0,/overplot
    fsc_plot,[180,180],[0,100],line=1,thick=1.0,/overplot
    fsc_plot,[-10,190],[45,45],line=1,thick=1.0,/overplot
    
    fsc_plot,barR.SA[start:*],const_rad[1,start:*],psym=16,symsize=0.5,$
             xrange=[-10,190],/xs,ytitle="Binned TW (0-3.6)"
    
    fsc_plot,[0,0],[0,100],line=1,thick=1.0,/overplot  
    fsc_plot,[90,90],[0,100],line=1,thick=1.0,/overplot
    fsc_plot,[180,180],[0,100],line=1,thick=1.0,/overplot
    fsc_plot,[-10,190],[45,45],line=1,thick=1.0,/overplot
    
    fsc_plot,barR.SA[start:*],const_rad[0,start:*],psym=16,symsize=0.5,$
             xrange=[-10,190],/xs,ytitle="Binned TW (3.6-7.2)"
    
    fsc_plot,[0,0],[0,100],line=1,thick=1.0,/overplot  
    fsc_plot,[90,90],[0,100],line=1,thick=1.0,/overplot
    fsc_plot,[180,180],[0,100],line=1,thick=1.0,/overplot
    fsc_plot,[-10,190],[45,45],line=1,thick=1.0,/overplot
    
    !p.multi = [0,0,0]

  end_PS
  
  ; RTW inner corr (folded)
  start_PS, workingPath+"pSS_OmegaP_RTWinner_corr.eps"
  
    angles = barR.SA[start:*] mod 90.0 ;first quadrant
    
    w = where(angles ge 45.0,count)
    angles[w] = 90.0 - angles[w] ;angular distance from closest axis
  
    fsc_plot,angles,omega_p_rad_mpfit[0,start:*],psym=16,symsize=0.5,xrange=[0,45],/xs,$
             ytitle="Bar "+textoidl("\Omega_p")+" [km/s/kpc]",yrange=[-90,100],/ys,$
             xtitle="",xtickname=replicate(' ',10),position=[0.2,0.4,0.9,0.9]
    
    fsc_plot,[0,0],[0,100],line=1,thick=1.0,/overplot  
    fsc_plot,[90,90],[0,100],line=1,thick=1.0,/overplot
    fsc_plot,[180,180],[0,100],line=1,thick=1.0,/overplot
    fsc_plot,[-10,190],[45,45],line=1,thick=1.0,/overplot
    
    nBins    = 9
    maxAngle = ceil(max(angles)/5.0)*5.0
    binSize  = maxAngle / nBins
    midBins  = findgen(nBins)/nBins * maxAngle + binSize/2.0
    sigma_omega_p = fltarr(nBins)
    
    for i=0,nBins-1 do begin
      w = where(angles ge binSize*i and angles lt binSize*(i+1),count)
      if (count ne 0) then begin
        sigma_omega_p[i] = stddev(omega_p_rad_mpfit[0,w+start])
        print,count,mean(omega_p_rad_mpfit[0,w+start]),stddev(omega_p_rad_mpfit[0,w+start])
      endif
    endfor
    
    fsc_plot,midBins,sigma_omega_p,psym=-16,symsize=0.6,xrange=[0,45],/xs,$
             ytitle=textoidl("\sigma_{\Omega_p}"),xtitle="Bar Angular Offset from Axes [deg]",$
             position=[0.2,0.2,0.9,0.4],/noerase
  
  end_PS
  
  stop
  
  ; IC-theory circular velocity lines
  units = getUnits()
  M_star = h2.f_d * h2.m200
  rPts = findgen(300)/300.0 * extent[0] + 0.02
  
  omega2 = func_omega2_ExpHern(units.G, rPts, h2.m200, M_star, h2.h, h2.a)
  kappa2 = func_kappa2_ExpHern(units.G, rPts, h2.m200, M_star, h2.h, h2.a)

  start_PS, workingPath+"pSS_OmegaP_R.eps"

    pS = 0.1
    colors = ['pink','green','skyblue']

    fsc_plot,[0],[0],/nodata, $
         xtitle="Radius / Disk Scale Length",ytitle="Pattern Speed [km/s/kpc]", $
         xrange=[0,max(extent)/h2.h],yrange=[0,80.0],/xs,/ys

    ;global constant
    polyfill,[0,max(extent),max(extent),0]/h2.h,$
             [mean(omega_p_const[start:*])-stddev(omega_p_const[start:*]),$
              mean(omega_p_const[start:*])-stddev(omega_p_const[start:*]),$
              mean(omega_p_const[start:*])+stddev(omega_p_const[start:*]),$
              mean(omega_p_const[start:*])+stddev(omega_p_const[start:*])],$
             color=fsc_color(colors[0]),/line_fill,spacing=pS,orientation=45.0
    fsc_plot,[0,max(extent)/h2.h],[mean(omega_p_const[start:*])-stddev(omega_p_const[start:*]),$
                              mean(omega_p_const[start:*])-stddev(omega_p_const[start:*])],$
                             color=fsc_color(colors[0]),/overplot
    fsc_plot,[0,max(extent)/h2.h],[mean(omega_p_const[start:*])+stddev(omega_p_const[start:*]),$
                              mean(omega_p_const[start:*])+stddev(omega_p_const[start:*])],$
                             color=fsc_color(colors[0]),/overplot                             
    
    ;binned TW
    polyfill,[0,max(extent)/2.0,max(extent)/2.0,0]/h2.h,$
             [mean(const_rad[1,start:*])-stddev(const_rad[1,start:*]),$
              mean(const_rad[1,start:*])-stddev(const_rad[1,start:*]),$
              mean(const_rad[1,start:*])+stddev(const_rad[1,start:*]),$
              mean(const_rad[1,start:*])+stddev(const_rad[1,start:*])],$
             color=fsc_color(colors[1]),/line_fill,spacing=pS,orientation=-45.0
    fsc_plot,[0,max(extent)/2.0/h2.h],[mean(const_rad[1,start:*])-stddev(const_rad[1,start:*]),$
                                  mean(const_rad[1,start:*])-stddev(const_rad[1,start:*])],$
                                  color=fsc_color(colors[1]),/overplot
    fsc_plot,[0,max(extent)/2.0/h2.h],[mean(const_rad[1,start:*])+stddev(const_rad[1,start:*]),$
                                  mean(const_rad[1,start:*])+stddev(const_rad[1,start:*])],$
                                  color=fsc_color(colors[1]),/overplot    
    polyfill,[max(extent)/2.0,max(extent),max(extent),max(extent)/2.0]/h2.h,$
             [mean(const_rad[0,start:*])-stddev(const_rad[0,start:*]),$
              mean(const_rad[0,start:*])-stddev(const_rad[0,start:*]),$
              mean(const_rad[0,start:*])+stddev(const_rad[0,start:*]),$
              mean(const_rad[0,start:*])+stddev(const_rad[0,start:*])],$
             color=fsc_color(colors[1]),/line_fill,spacing=pS,orientation=-45.0             
    fsc_plot,[max(extent)/2.0,max(extent)]/h2.h,[mean(const_rad[0,start:*])-stddev(const_rad[0,start:*]),$
                                            mean(const_rad[0,start:*])-stddev(const_rad[0,start:*])],$
                                            color=fsc_color(colors[1]),/overplot
    fsc_plot,[max(extent)/2.0,max(extent)]/h2.h,[mean(const_rad[0,start:*])+stddev(const_rad[0,start:*]),$
                                            mean(const_rad[0,start:*])+stddev(const_rad[0,start:*])],$
                                            color=fsc_color(colors[1]),/overplot    
    
    ;RTW broken const model
    polyfill,[0,max(extent)/2.0,max(extent)/2.0,0]/h2.h,$
             [mean(omega_p_rad_mpfit[0,start:*])-stddev(omega_p_rad_mpfit[0,start:*]),$
              mean(omega_p_rad_mpfit[0,start:*])-stddev(omega_p_rad_mpfit[0,start:*]),$
              mean(omega_p_rad_mpfit[0,start:*])+stddev(omega_p_rad_mpfit[0,start:*]),$
              mean(omega_p_rad_mpfit[0,start:*])+stddev(omega_p_rad_mpfit[0,start:*])],$
             color=fsc_color(colors[2]),/line_fill,spacing=pS,orientation=45.0
    fsc_plot,[0,max(extent)/2.0/h2.h],[mean(omega_p_rad_mpfit[0,start:*])-stddev(omega_p_rad_mpfit[0,start:*]),$
                                  mean(omega_p_rad_mpfit[0,start:*])-stddev(omega_p_rad_mpfit[0,start:*])],$
                                  color=fsc_color(colors[2]),/overplot
    fsc_plot,[0,max(extent)/2.0/h2.h],[mean(omega_p_rad_mpfit[0,start:*])+stddev(omega_p_rad_mpfit[0,start:*]),$
                                  mean(omega_p_rad_mpfit[0,start:*])+stddev(omega_p_rad_mpfit[0,start:*])],$
                                  color=fsc_color(colors[2]),/overplot
    polyfill,[max(extent)/2.0,max(extent),max(extent),max(extent)/2.0]/h2.h,$
             [mean(omega_p_rad_mpfit[1,start:*])-stddev(omega_p_rad_mpfit[1,start:*]),$
              mean(omega_p_rad_mpfit[1,start:*])-stddev(omega_p_rad_mpfit[1,start:*]),$
              mean(omega_p_rad_mpfit[1,start:*])+stddev(omega_p_rad_mpfit[1,start:*]),$
              mean(omega_p_rad_mpfit[1,start:*])+stddev(omega_p_rad_mpfit[1,start:*])],$
             color=fsc_color(colors[2]),/line_fill,spacing=pS,orientation=45.0
    fsc_plot,[max(extent)/2.0,max(extent)]/h2.h,[mean(omega_p_rad_mpfit[1,start:*])-stddev(omega_p_rad_mpfit[1,start:*]),$
                                            mean(omega_p_rad_mpfit[1,start:*])-stddev(omega_p_rad_mpfit[1,start:*])],$
                                            color=fsc_color(colors[2]),/overplot
    fsc_plot,[max(extent)/2.0,max(extent)]/h2.h,[mean(omega_p_rad_mpfit[1,start:*])+stddev(omega_p_rad_mpfit[1,start:*]),$
                                            mean(omega_p_rad_mpfit[1,start:*])+stddev(omega_p_rad_mpfit[1,start:*])],$
                                            color=fsc_color(colors[2]),/overplot                                  
                                  
    ;45 (fourier/barspeed)
    fsc_plot,[0,max(extent)/2.0/h2.h],[45.0,45.0],line=2,thick=4.0,/overplot

    ; barSpeed results                         
    fracThreshold = 0.05
    w = where(barR.frac ge fracThreshold)
  
    fsc_plot,barR.radii[w]/h2.h,barR.angMean[w],psym=16,symsize=0.8,/overplot
    
    ; errors
    for i=0,n_elements(w)-1 do begin
      fsc_plot,[barR.radii[w[i]],barR.radii[w[i]]]/h2.h,$
               [barR.angMean[w[i]]-barR.angMeanErr[w[i]],barR.angMean[w[i]]+barR.angMeanErr[w[i]]],$
               line=0,/overplot
    endfor

    ;overplot omega, omega-kappa/2, omega+kappa/2
    fsc_plot,rPts/h2.h,sqrt(omega2),color=fsc_color('orange'),/overplot
    fsc_plot,rPts/h2.h,sqrt(omega2)-sqrt(kappa2)/2.0,color=fsc_color('orange'),line=2,/overplot
    fsc_plot,rPts/h2.h,sqrt(omega2)+sqrt(kappa2)/2.0,color=fsc_color('orange'),line=2,/overplot
    fsc_plot,rPts/h2.h,sqrt(omega2)-sqrt(kappa2)/4.0,color=fsc_color('orange'),line=1,/overplot
    fsc_plot,rPts/h2.h,sqrt(omega2)+sqrt(kappa2)/4.0,color=fsc_color('orange'),line=1,/overplot   

  end_PS
  
end

; patSpeedSeriesSpiral(): run many

pro patSpeedSeriesSpiral, simName, barePlot=barePlot

  basePath = '/n/home07/dnelson/spirals/lambdacrit/sims/'

  filePath    = basePath + simName + "/"
  workingPath = filePath + "/output3/"
  
  saveFileName = filePath+"pSS.30.sav"
  
  nSnaps = n_elements(file_search(filePath+"output/snap_*"))
  
  h2 = loadSimParams(filePath)  
  
  extent   = [10.0]*h2.h
  nRadBins = 8
  sizeRes  = 2
  
  ; plot config
  pS = 0.1
  colors = ['pink','skyblue','forest green','purple','cyan','orange']
  ;colors = ['black','gray','dark gray','black','black','black'] ;grayscale
  
  if not (file_test(saveFileName)) then begin
  
    omega_p_const          = fltarr(nSnaps)
    omega_p_const_rad      = fltarr(nRadBins,nSnaps)
    omega_p_const_ypts     = fltarr(nRadBins,nSnaps)
    omega_p_rad_mpfit      = fltarr(sizeRes,nSnaps)
    omega_p_rad_mpfit_errs = fltarr(sizeRes,nSnaps)
    times                  = fltarr(nSnaps)
    
    for m=0,nSnaps-1,1 do begin 
      r = patSpeedConst(filePath, workingPath, m, extent, nRadBins)
      
      omega_p_const[m]        = r.omega_p     ;extents
      omega_p_const_rad[*,m]  = r.avgOmegaP   ;radially binned
      omega_p_const_ypts[*,m] = r.avgYPts     ;radial midbins
      
      times[m]                = r.time
      
      r = patSpeedRadialMpfit(simName, m)
      
      omega_p_rad_mpfit[*,m]      = r.res     ;omega_amp,omega_slope
      omega_p_rad_mpfit_errs[*,m] = r.perrors
    endfor
  
    ; save/restore
    save,omega_p_const,omega_p_const_rad,omega_p_const_ypts,omega_p_rad_mpfit,$
       omega_p_rad_mpfit_errs,times,filename=saveFileName
  endif else begin
    restore,saveFileName
  endelse
  
  nSnaps = n_elements(times) ;override nSnaps in the case we have added more snaps but not in save

  ; reform radially binned
  const_rad = fltarr(nRadBins/2,nSnaps)
  ypts      = fltarr(nRadBins/2)
  
  for m=0,nSnaps-1,1 do begin
    for i=0,nRadBins/2-1 do begin
      val1 = omega_p_const_rad[i,m]
      val2 = omega_p_const_rad[nRadBins-i-1,m]
      ypts[i] = omega_p_const_ypts[nRadBins-i-1,m]
      
      if (val1 gt 0 and val2 gt 0) then $
        const_rad[i,m] = (val1+val2)/2.0
      if (val1 gt 0 and val2 lt 0) then $
        const_rad[i,m] = val1
      if (val1 lt 0 and val2 gt 0) then $
        const_rad[i,m] = val2 
    endfor
  endfor
  
  ; load barSpeed results
  barR = barSpeed(simName)
  
  ; t > 200 Myr = 4 (LC spirals)
  start = 4
  
  ; print mean+stddev
  print,'global constant ', mean(omega_p_const[start:*]),stddev(omega_p_const[start:*])
  print,'binned tw (0-7.5) ',   mean(const_rad[3,start:*]),stddev(const_rad[3,start:*])
  print,'binned tw (7.5-15) ',  mean(const_rad[2,start:*]),stddev(const_rad[2,start:*])
  print,'binned tw (15-22.5) ', mean(const_rad[1,start:*]),stddev(const_rad[1,start:*]) 
  print,'binned tw (22.5-30) ', mean(const_rad[0,start:*]),stddev(const_rad[0,start:*])
  print,'rtw Amp (at h) ',  mean(omega_p_rad_mpfit[0,start:*]),stddev(omega_p_rad_mpfit[0,start:*])
  print,'rtw slope ',       mean(omega_p_rad_mpfit[1,start:*]),stddev(omega_p_rad_mpfit[1,start:*])
 
  ; plot
  if not keyword_set(barePlot) then begin
  start_PS, workingPath+"pSS_OmegaP_histos_30.eps"

    !p.multi = [0,2,4]

    histoplot,omega_p_const[start:*],title="Global Constant",/fillpolygon,nbins=8
    fsc_plot,[mean(omega_p_const[start:*]),mean(omega_p_const[start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[median(omega_p_const[start:*]),median(omega_p_const[start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[2]),/overplot
             
    histoplot,const_rad[3,start:*],title="Binned TW (0-7.5 kpc)",/fillpolygon,nbins=8
    fsc_plot,[mean(const_rad[3,start:*]),mean(const_rad[3,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[median(const_rad[3,start:*]),median(const_rad[3,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[2]),/overplot
    histoplot,const_rad[2,start:*],title="Binned TW (7.5-15 kpc)",/fillpolygon,nbins=8
    fsc_plot,[mean(const_rad[2,start:*]),mean(const_rad[2,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[median(const_rad[2,start:*]),median(const_rad[2,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[2]),/overplot
    histoplot,const_rad[1,start:*],title="Binned TW (15-22.5 kpc)",/fillpolygon,nbins=8
    fsc_plot,[mean(const_rad[1,start:*]),mean(const_rad[1,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[median(const_rad[1,start:*]),median(const_rad[1,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[2]),/overplot
    histoplot,const_rad[0,start:*],title="Binned TW (22.5-30 kpc)",/fillpolygon,nbins=8 
    fsc_plot,[mean(const_rad[0,start:*]),mean(const_rad[0,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[median(const_rad[0,start:*]),median(const_rad[0,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[2]),/overplot  
             
    histoplot,omega_p_rad_mpfit[0,start:*],title="RTW Amp (at h)",/fillpolygon,nbins=8
    fsc_plot,[mean(omega_p_rad_mpfit[0,start:*]),mean(omega_p_rad_mpfit[0,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[median(omega_p_rad_mpfit[0,start:*]),median(omega_p_rad_mpfit[0,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[2]),/overplot
    histoplot,omega_p_rad_mpfit[1,start:*],title="RTW Slope",/fillpolygon,nbins=8
    fsc_plot,[mean(omega_p_rad_mpfit[1,start:*]),mean(omega_p_rad_mpfit[1,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[median(omega_p_rad_mpfit[1,start:*]),median(omega_p_rad_mpfit[1,start:*])],[0.0,10.0],$
             line=1,color=fsc_color(colors[2]),/overplot
    
    fsc_text,0.8,0.15,"Mean",color=fsc_color(colors[1]),charsize=1.2,/normal
    fsc_text,0.8,0.10,"Median",color=fsc_color(colors[2]),charsize=1.2,/normal
    
    !p.multi = [0,0,0]

  end_PS
  
  start_PS, workingPath+"pSS_OmegaP_wtime_30.eps"

    !p.multi = [0,2,4]

    fsc_plot,times[start:*],omega_p_const[start:*],title="Global Constant",xrange=[0.0,1.0],/xs,/ys
    fsc_plot,[0.0,1.0],[mean(omega_p_const[start:*]),mean(omega_p_const[start:*])],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[0.0,1.0],[median(omega_p_const[start:*]),median(omega_p_const[start:*])],$
             line=1,color=fsc_color(colors[2]),/overplot  
                        
    fsc_plot,times[start:*],const_rad[3,start:*],title="Binned TW (0-7.5 kpc)",xrange=[0.0,1.0],/xs,/ys
    fsc_plot,[0.0,1.0],[mean(const_rad[3,start:*]),mean(const_rad[3,start:*])],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[0.0,1.0],[median(const_rad[3,start:*]),median(const_rad[3,start:*])],$
             line=1,color=fsc_color(colors[2]),/overplot 
    fsc_plot,times[start:*],const_rad[2,start:*],title="Binned TW (7.5-15 kpc)",xrange=[0.0,1.0],/xs,/ys
    fsc_plot,[0.0,1.0],[mean(const_rad[2,start:*]),mean(const_rad[2,start:*])],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[0.0,1.0],[median(const_rad[2,start:*]),median(const_rad[2,start:*])],$
             line=1,color=fsc_color(colors[2]),/overplot 
    fsc_plot,times[start:*],const_rad[1,start:*],title="Binned TW (15-22.5 kpc)",xrange=[0.0,1.0],/xs,/ys
    fsc_plot,[0.0,1.0],[mean(const_rad[1,start:*]),mean(const_rad[1,start:*])],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[0.0,1.0],[median(const_rad[1,start:*]),median(const_rad[1,start:*])],$
             line=1,color=fsc_color(colors[2]),/overplot 
    fsc_plot,times[start:*],const_rad[0,start:*],title="Binned TW (22.5-30 kpc)",xrange=[0.0,1.0],/xs,/ys  
    fsc_plot,[0.0,1.0],[mean(const_rad[0,start:*]),mean(const_rad[0,start:*])],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[0.0,1.0],[median(const_rad[0,start:*]),median(const_rad[0,start:*])],$
             line=1,color=fsc_color(colors[2]),/overplot 

    fsc_plot,times[start:*],omega_p_rad_mpfit[0,start:*],title="RTW Amp (at h)",xrange=[0.0,1.0],/xs,/ys,$
             xtitle="Time [Gyr]"
    fsc_plot,[0.0,1.0],[mean(omega_p_rad_mpfit[0,start:*]),mean(omega_p_rad_mpfit[0,start:*])],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[0.0,1.0],[median(omega_p_rad_mpfit[0,start:*]),median(omega_p_rad_mpfit[0,start:*])],$
             line=1,color=fsc_color(colors[2]),/overplot
    fsc_plot,times[start:*],omega_p_rad_mpfit[1,start:*],title="RTW Slope",xrange=[0.0,1.0],/xs,/ys
    fsc_plot,[0.0,1.0],[mean(omega_p_rad_mpfit[1,start:*]),mean(omega_p_rad_mpfit[1,start:*])],$
             line=1,color=fsc_color(colors[1]),/overplot
    fsc_plot,[0.0,1.0],[median(omega_p_rad_mpfit[1,start:*]),median(omega_p_rad_mpfit[1,start:*])],$
             line=1,color=fsc_color(colors[2]),/overplot
                 
    fsc_text,0.8,0.15,"Mean",color=fsc_color(colors[1]),charsize=1.2,/normal
    fsc_text,0.8,0.10,"Median",color=fsc_color(colors[2]),charsize=1.2,/normal
    
    !p.multi = [0,0,0]

  end_PS
  endif ;barePlot

  ; IC-theory circular velocity lines
  units = getUnits()
  M_star = h2.f_d * h2.m200
  rPts = findgen(300)/300.0 * extent[0] + 0.01
  
  omega2 = func_omega2_ExpHern(units.G, rPts, h2.m200, M_star, h2.h, h2.a)
  kappa2 = func_kappa2_ExpHern(units.G, rPts, h2.m200, M_star, h2.h, h2.a)

  if not keyword_set(barePlot) then begin
  start_PS, workingPath+"pSS_OmegaP_R.eps"

    fsc_plot,[0],[0],/nodata, $
         xtitle="Radius / Disk Scale Length",ytitle="Pattern Speed [km/s/kpc]", $
         xrange=[0,max(extent/2)/h2.h],yrange=[0,80.0],/xs,/ys
  endif ;barePlot
    
    ;binned TW
    nRadBins = float(nRadBins)
    
    for i=0.0,nRadBins/2-1 do begin
      ;print,i,max(extent/2)*(i/(nRadBins/2)),max(extent/2)*(i+1)/(nRadBins/2)
      polyfill,[max(extent/2)*(i/(nRadBins/2)),max(extent/2)*(i+1)/(nRadBins/2),$
                max(extent/2)*(i+1)/(nRadBins/2),max(extent/2)*(i/(nRadBins/2))]/h2.h,$
               [mean(const_rad[(nRadBins/2-1)-i,start:*])-stddev(const_rad[(nRadBins/2-1)-i,start:*]),$
                mean(const_rad[(nRadBins/2-1)-i,start:*])-stddev(const_rad[(nRadBins/2-1)-i,start:*]),$
                mean(const_rad[(nRadBins/2-1)-i,start:*])+stddev(const_rad[(nRadBins/2-1)-i,start:*]),$
                mean(const_rad[(nRadBins/2-1)-i,start:*])+stddev(const_rad[(nRadBins/2-1)-i,start:*])],$
               color=fsc_color(colors[1]),/line_fill,spacing=pS,orientation=-45.0
      fsc_plot,[max(extent/2)*(i/(nRadBins/2)),max(extent/2)*(i+1)/(nRadBins/2)]/h2.h,$
               [mean(const_rad[(nRadBins/2-1)-i,start:*])-stddev(const_rad[(nRadBins/2-1)-i,start:*]),$
                mean(const_rad[(nRadBins/2-1)-i,start:*])-stddev(const_rad[(nRadBins/2-1)-i,start:*])],$
                color=fsc_color(colors[1]),/overplot
      fsc_plot,[max(extent/2)*(i/(nRadBins/2)),max(extent/2)*(i+1)/(nRadBins/2)]/h2.h,$
               [mean(const_rad[(nRadBins/2-1)-i,start:*])+stddev(const_rad[(nRadBins/2-1)-i,start:*]),$
                mean(const_rad[(nRadBins/2-1)-i,start:*])+stddev(const_rad[(nRadBins/2-1)-i,start:*])],$
                color=fsc_color(colors[1]),/overplot    
    endfor

    ;RTW powerlaw model
    oPts    = fltarr(n_elements(rPts))
    oPtsMax = fltarr(n_elements(rPts))
    oPtsMin = fltarr(n_elements(rPts))
    
    for i=0,n_elements(rPts)-1 do begin
      op = omega_p_rad_mpfit[0,start:*] / (rPts[i]/h2.h)^(omega_p_rad_mpfit[1,start:*])
      
      oPts[i]    = mean(op)
      oPtsMax[i] = mean(op)+stddev(op)
      oPtsMin[i] = mean(op)-stddev(op)
    endfor
    
    ; restrict range
    w = where(rPts ge 0.0*h2.h and rPts lt max(extent/2))
    
    ; visually fix first two points
    oPtsMin[0] = (linfit(rPts[2:5],oPtsMin[2:5]))[0]
    oPtsMin[1] = (oPtsMin[0]+oPtsMin[2])/2.0
    
    polyfill,[rPts[w],reverse(rPts[w]),rPts[min(w)]]/h2.h,$
             [oPtsMax[w],reverse(oPtsMin[w]),oPtsMax[min(w)]]<80.0,$
             color=fsc_color(colors[2]),/line_fill,spacing=pS,orientation=45.0
           
    ;fsc_plot,rPts/h2.h,oPts,line=0,thick=4.0,color=fsc_color(colors[2]),/overplot
    fsc_plot,rPts/h2.h,oPtsMax,line=0,thick=2.0,color=fsc_color(colors[2]),/overplot
    fsc_plot,rPts/h2.h,oPtsMin,line=0,thick=2.0,color=fsc_color(colors[2]),/overplot
             
    ; barSpeed results                         
    ;fracThreshold = 0.05
    ;w = where(barR.frac ge fracThreshold)
    w = where(barR.radii gt 0.5*h2.h and barR.radii lt 5.0*h2.h)
  
    fsc_plot,barR.radii[w]/h2.h,barR.angMean[w],psym=16,symsize=1.2,/overplot
             
    ; errors
    for i=0,n_elements(w)-1 do begin
      fsc_plot,[barR.radii[w[i]],barR.radii[w[i]]]/h2.h,$
               [barR.angMean[w[i]]-barR.angMeanErr[w[i]],barR.angMean[w[i]]+barR.angMeanErr[w[i]]],$
               line=0,thick=!p.thick+1,/overplot
    endfor
             
    ;overplot omega, omega-kappa/2, omega+kappa/2
    fsc_plot,rPts/h2.h,sqrt(omega2),color=fsc_color(colors[5]),/overplot,thick=!p.thick+2,line=2
    fsc_plot,rPts/h2.h,sqrt(omega2)-sqrt(kappa2)/2.0,color=fsc_color(colors[5]),line=2,/overplot
    fsc_plot,rPts/h2.h,sqrt(omega2)+sqrt(kappa2)/2.0,color=fsc_color(colors[5]),line=2,/overplot
    fsc_plot,rPts/h2.h,sqrt(omega2)-sqrt(kappa2)/4.0,color=fsc_color(colors[5]),line=1,/overplot
    fsc_plot,rPts/h2.h,sqrt(omega2)+sqrt(kappa2)/4.0,color=fsc_color(colors[5]),line=1,/overplot   

  if not keyword_set(barePlot) then begin
  end_PS
  endif ;barePlot

end

; patSpeedSpiralMult():

pro patSpeedSpiralMult

  workingPath = '/n/home07/dnelson/spirals/lambdacrit/'
  
  simNames = ['LC_10m_disk_2','LC_10m_disk_4','LC_10m_disk_6','LC_10m_disk_8']
  
  minMaxY = [0.0,80.0]
  
  start_PS, workingPath+"pSS_OmegaP_R_mult.eps"
  
    xm = !x.margin
    ym = !y.margin
    !x.margin = [0,0]
    !y.margin = [0,0]

    ; upper left (LC2)
    h2       = loadSimParams(simNames[0])  
    extent   = [10.0]*h2.h
    
    fsc_plot,[0],[0],/nodata, $
         xtitle="",ytitle="",xtickname=replicate(' ',10), $
         xrange=[0,max(extent/2)/h2.h],yrange=minMaxY+0.01,xs=1,ys=1,position=[0.1,0.525,0.525,0.95]

    patSpeedSeriesSpiral,simNames[0],/barePlot
    fsc_text,4.3,66,"LC-2",alignment=0.5,/data

    ; upper right (LC4)
    h2       = loadSimParams(simNames[1])  
    extent   = [10.0]*h2.h
    
    fsc_plot,[0],[0],/nodata,/noerase, $
         xtitle="",ytitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10), $
         xrange=[0,max(extent/2)/h2.h],yrange=minMaxY,xs=1,ys=1,position=[0.525,0.525,0.95,0.95]
      
    patSpeedSeriesSpiral,simNames[1],/barePlot
    fsc_text,4.3,66,"LC-4",alignment=0.5,/data
    
    ; lower left (LC6)
    h2       = loadSimParams(simNames[2])  
    extent   = [10.0]*h2.h
    
    fsc_plot,[0],[0],/nodata,/noerase, $
         xtitle="",ytitle="", $
         xrange=[0,max(extent/2)/h2.h],yrange=minMaxY,xs=1,ys=1,position=[0.1,0.15,0.525,0.525]
    
    patSpeedSeriesSpiral,simNames[2],/barePlot
    fsc_text,4.3,66,"LC-6",alignment=0.5,/data
    
    ; lower right (LC8)
    h2       = loadSimParams(simNames[3])  
    extent   = [10.0]*h2.h
    
    fsc_plot,[0],[0],/nodata,/noerase, $
         xtitle="",ytitle="",ytickname=replicate(' ',10), $
         xrange=[0.01,max(extent/2)/h2.h],yrange=minMaxY,xs=1,ys=1,position=[0.525,0.15,0.95,0.525]
      
    patSpeedSeriesSpiral,simNames[3],/barePlot
    fsc_text,4.3,66,"LC-8",alignment=0.5,/data
    
    ; axes labels
    fsc_text,0.03,0.5,"Pattern Speed [km/s/kpc]",alignment=0.5,orientation=90.0,/normal
    fsc_text,0.5,0.02,"Radius / Disk Scale Length",alignment=0.5,orientation=0.0,/normal
    
    if 0 then begin
    ; break radius / sgN (upper left)
    fsc_plot,[0],[0],/nodata,yrange=[minY,maxY],xrange=[min(param1),max(param1)*0.99],$
             xtitle="",ytitle="",xs=9,/ys,position=[0.15,0.5,0.5,0.85],xtickname=replicate(' ',10)
    fsc_plot,param1,omegaP1[0,*],psym=-16,symsize=0.4,color=fsc_color('orange'),/overplot
    fsc_plot,param1,omegaP1[1,*]*100,psym=-16,symsize=0.4,color=fsc_color('blue'),/overplot 
    ;fsc_plot,[min(param1),max(param1)],[45.0,45.0],line=1,/overplot 
    fsc_axis,0.0,maxY,0.0,/xaxis,/xs,xrange=[min(param1),max(param1)*0.99]
    fsc_text,(0.5+0.15)/2.0,0.93,"SG Filter Size [#]",alignment=0.5,/normal

    ; y_r_max (upper right)
    fsc_plot,[0],[0],/nodata,yrange=[minY,maxY],xrange=[min(param2),max(param2)],$
             xtitle="",ytitle="",xs=9,ys=9,position=[0.5,0.5,0.85,0.85],/noerase,$
             xtickname=replicate(' ',10),ytickname=replicate(' ',10)
    fsc_plot,param2,omegaP2[0,*],psym=-16,symsize=0.4,color=fsc_color('orange'),/overplot
    fsc_plot,param2,omegaP2[1,*]*100,psym=-16,symsize=0.4,color=fsc_color('blue'),/overplot
    ;fsc_plot,[min(param2),max(param2)],[45.0,45.0],line=1,/overplot 
    fsc_axis,0.0,maxY,0.0,/xaxis,/xs,xrange=[min(param2),max(param2)]
    fsc_axis,max(param2),0.0,0.0,/yaxis,/ys,yrange=[minY,maxY]
    fsc_text,(0.85+0.5)/2.0,0.93,"Maximum Y [kpc]",alignment=0.5,/normal
    
    ; yWidth (lower left)
    fsc_plot,[0],[0],/nodata,yrange=[minY,maxY],xrange=[min(param3),max(param3)*0.99],$
             xtitle="yWidth [kpc]",ytitle="",/xs,/ys,position=[0.15,0.15,0.5,0.5],/noerase
    fsc_plot,param3,omegaP3[0,*],psym=-16,symsize=0.4,color=fsc_color('orange'),/overplot
    fsc_plot,param3,omegaP3[1,*]*100,psym=-16,symsize=0.4,color=fsc_color('blue'),/overplot 
    ;fsc_plot,[min(param3),max(param3)],[45.0,45.0],line=1,/overplot 
    
    ; xBins (lower right)
    fsc_plot,[0],[0],/nodata,yrange=[minY,maxY],xrange=[min(param4),max(param4)],$
             xtitle="xBins [#]",ytitle="",/xs,ys=9,position=[0.5,0.15,0.85,0.5],/noerase,ytickname=replicate(' ',10)
    fsc_plot,param4,omegaP4[0,*],psym=-16,symsize=0.4,color=fsc_color('orange'),/overplot
    fsc_plot,param4,omegaP4[1,*]*100,psym=-16,symsize=0.4,color=fsc_color('blue'),/overplot 
    ;fsc_plot,[min(param4),max(param4)],[45.0,45.0],line=1,/overplot 
    fsc_axis,max(param4),0.0,0.0,/yaxis,/ys,yrange=[minY,maxY]
    
    ; legend
    fsc_text,max(param4)*0.95,58,"Amplitude",color=fsc_color('orange'),alignment=1.0,charsize=1.5
    fsc_text,max(param4)*0.95,50,"Slope",color=fsc_color('blue'),alignment=1.00,charsize=1.5
    
    ; y axes
    fsc_text,0.07,0.5,"Amplitude [km/s/kpc] / 100x Slope",alignment=0.5,orientation=90.0,/normal
    endif ;0
    
    !x.margin = xm
    !y.margin = ym
  end_PS

end