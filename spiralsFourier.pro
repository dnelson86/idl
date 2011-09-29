; spiralsFourier.pro
; fourier analysis routines
; dnelson may.2011

; fourierContour(): load C(r,t) and compute C(r,omega)

pro fourierContour, simName

  units  = getUnits()

  ; config
  basePath = "/n/home07/dnelson/spirals/lambdacrit/sims/"
  
  filePath = basePath + simName + "/"
  fourierPath = basePath + simName + "/output2/"
  workingPath = basePath + simName + "/output3/"

  modes = [1,2,3,4,5,6,7,8]
  
  h2 = loadSimParams(simName)
  
  nSnaps = n_elements(file_search(fourierPath+"fourierSaved_*.sav"))
  
  ; restrict range ;bar_1m_a: [600,1000] ;LC_10m_disk_2: [200,1000]
  startTime = 200.0
  endTime   = 1000.0
 
  snapStart = 0
  for j=0,nSnaps-1 do begin
    savePath = fourierPath+"fourierSaved_"+str(j+1)+".sav"
    restore,savePath
 
    if (time le startTime) then $
      snapStart = j
    if (time le endTime) then $
      snapEnd = j
  endfor
  
  nSnaps = snapEnd - snapStart
  
  print,'snap start end number startTime endTime',snapStart,snapEnd,nSnaps,startTime,endTime
  
  ; load sizing info
  savePath = fourierPath+"fourierSaved_1.sav"
  restore,savePath
  
  nRadialBins   = n_elements(f.R)
  radialMidBins = f.R
  
  ; 2h slice arrays
  sliceRadii = [0.5,1.0,2.0,3.0] * h2.h
  nSliceFreqs = 100 ;m101=18
  sliceAmps  = fltarr(nSliceFreqs,n_elements(modes),n_elements(sliceRadii))
  
  ; loop over requested modes
  for k=0,n_elements(modes)-1 do begin
  
    m=modes[k]
    print,' m=',str(m)

    ; result arrays  
    A   = fltarr(nRadialBins,nSnaps)
    phi = fltarr(nRadialBins,nSnaps)
    aj  = fltarr(nRadialBins,nSnaps)
    bj  = fltarr(nRadialBins,nSnaps)
    
    ; c(r,t) arrays
    c_a = fltarr(nRadialBins,nSnaps)
    c_b = fltarr(nRadialBins,nSnaps)
    c = complexarr(nRadialBins,nSnaps)
    c_interp = complexarr(nRadialBins,nSnaps)
    times = fltarr(nSnaps)
    
    ; load
    for j=0,nSnaps-1 do begin
      savePath = fourierPath+"fourierSaved_"+str(j+1+snapStart)+".sav"
      restore,savePath
      
      c_a[*,j] = f.aj[*,m]
      c_b[*,j] = f.bj[*,m]
      c[*,j] = complex(f.aj[*,m],f.bj[*,m])
      
      times[j] = time
    endfor
    
    times -= min(times)
    
    ; evenly spaced frequencies
    period = max(times)-min(times)
    dt = period / nSnaps
    
    tarTimes = findgen(nSnaps)/nSnaps * period
    freqs = findgen(nSnaps)/(nSnaps*dt)
    
    ; interpolate c onto regularly gridded tarTimes
    for i=0,nRadialBins-1 do begin
      aj = interpol(reform(c_a[i,*]),times,tarTimes,/spline)
      bj = interpol(reform(c_b[i,*]),times,tarTimes,/spline)
      c_interp[i,*] = complex(aj,bj)
    endfor
    
    ; IDL FFT (OF INTERP!)
    fft_idl = fft(c_interp,DIMENSION=2)
    A_idl   = fltarr(nRadialBins,nSnaps)
    
    for i=0,nRadialBins-1 do begin
      A_idl[i,*] = sqrt( real_part(fft_idl[i,*])^2.0 + imaginary(fft_idl[i,*])^2.0 )
    endfor
   
    ; Manual FFT (disabled)
    if 0 then begin
    for i=0,nRadialBins-1 do begin
      ca = reform(c_a[i,*])
      cb = reform(c_b[i,*])
  
      A[i,0] = 1.0/nSnaps * total(sqrt(ca^2.0 + cb^2.0)) ;??? (dt/period)=(1/N)
      
      for j=0,nSnaps-1 do begin
        aj[i,j] = 1.0/nSnaps * total(ca * cos(-2*!pi*freqs[j]*tarTimes) - $
                                     cb * sin(-2*!pi*freqs[j]*tarTimes))
        bj[i,j] = 1.0/nSnaps * total(cb * cos(-2*!pi*freqs[j]*tarTimes) + $
                                     ca * sin(-2*!pi*freqs[j]*tarTimes))
  
        arctan,bj[i,j],aj[i,j],aTanVal,aTanValDeg
        
        A[i,j] = sqrt(aj[i,j]^2.0 + bj[i,j]^2.0)
        phi[i,j] = aTanVal
        ;print,j,aj[i,j],bj[i,j],A[i,j],phi[i,j]
        ;print,total(ca * cos(-2*!pi*freqs[j]*tarTimes) - cb * sin(-2*!pi*freqs[j]*tarTimes)),$
        ;      total(cb * cos(-2*!pi*freqs[j]*tarTimes) + ca * sin(-2*!pi*freqs[j]*tarTimes)),$
        ;      total( complex(cos(-2*!pi*freqs[j]*tarTimes),sin(-2*!pi*freqs[j]*tarTimes))*complex(ca,cb))
      endfor
    endfor
    endif
    
    ; make IC-theory circular frequency curves
    M_star = h2.f_d * h2.m200
    rPts = findgen(300)/300.0 * max(radialMidBins) + 0.1
    
    omega2 = func_omega2_ExpHern(units.G, rPts, h2.m200, M_star, h2.h, h2.a)
    kappa2 = func_kappa2_ExpHern(units.G, rPts, h2.m200, M_star, h2.h, h2.a)
    
    ; contour plot
    freqs_plot = 2*!pi * freqs * units.kpc_in_km/units.s_in_Myr ;angular km/s/kpc
    
    start_PS, workingPath+"fCont.m="+str(m)+".eps"
      nl      = 15.0
      nSurf   = 200
      radMax  = 5.0
      freqMax = 220.0 ;bar=120.0 ;spiral=220
      
      w = where(freqs_plot le freqMax*1.2)
      
      levels = ( findgen(nl)/(nl+1) + (1/nl) ) * max(A_idl[*,w])
      levels = levels[ceil(nl/10):*]
      
      radSpacing  = radMax / (nSurf-1)
      freqSpacing = freqMax / (nSurf-1)
      
      conSurf = tri_surf(A_idl[*,w],xvalues=radialMidBins/h2.h,yvalues=freqs_plot[w],$
                               gs=[radSpacing,freqSpacing],bounds=[0.0,0.0,radMax,freqMax])
                               
      conSurfX = findgen(nSurf)/nSurf * radMax
      conSurfY = findgen(nSurf)/nSurf * freqMax
      
      fsc_contour,conSurf,conSurfX,consurfY,$
      ;fsc_contour,A_idl,radialMidBins/h2.h,freqs_plot,charsize=1.2,$
        xtitle="Radius / Disk Scale Length",ytitle="Circular Frequency [km/s/kpc]",$
        levels=levels,c_labels=fltarr(nl),yrange=[0,freqMax],xrange=[0,4],/xs,/ys
  
      ;overplot omega, omega+-kappa/2, omega+-kappa/4
      fsc_plot,rPts/h2.h,sqrt(omega2),color=fsc_color('orange'),/overplot
      fsc_plot,rPts/h2.h,sqrt(omega2)-sqrt(kappa2)/2.0,color=fsc_color('orange'),line=2,/overplot
      fsc_plot,rPts/h2.h,sqrt(omega2)+sqrt(kappa2)/2.0,color=fsc_color('orange'),line=2,/overplot
      fsc_plot,rPts/h2.h,sqrt(omega2)-sqrt(kappa2)/4.0,color=fsc_color('orange'),line=1,/overplot
      fsc_plot,rPts/h2.h,sqrt(omega2)+sqrt(kappa2)/4.0,color=fsc_color('orange'),line=1,/overplot    
    end_PS
  
    ; store slices
    for i=0,n_elements(sliceRadii)-1 do begin
      w = where(abs(radialMidBins-sliceRadii[i]) eq min(abs(radialMidBins-sliceRadii[i])))
      sliceAmps[*,k,i] = A_idl[w,0:nSliceFreqs-1]
    endfor
  
  endfor ;m
  
  ; plot mode amplitude as a function of freq at fixed radius for all modes
  for i=0,n_elements(sliceRadii)-1 do begin
  
    start_PS,workingPath+"fCont.slice."+str(i)+".eps"
      
      fsc_plot,[0],[0],/nodata,xrange=[0,freqMax],yrange=[0,max(sliceAmps)*1.1],$
               xtitle="Angular Frequency [km/s/kpc]",ytitle="Mode Amplitude",/xs,/ys    
      
      for k=0,n_elements(modes)-1 do begin
        fsc_plot,freqs_plot,sliceAmps[*,k,i],color=fsc_color(units.colors[k+1]),thick=4.0,/overplot
  
        ; legend
        legendYStart = max(sliceAmps)*1.03
        legendYSpace = ( max(sliceAmps)-min(sliceAmps) )/22.0
              
        fsc_plot,[freqMax*0.05,freqMax*0.10],[legendYStart-legendYSpace*k,legendYStart-legendYSpace*k],$
              color=fsc_color(units.colors[k+1]),thick=4.0,/overplot
        fsc_text,freqMax*0.15,legendYStart-legendYSpace/5-legendYSpace*k,"m = "+str(modes[k]),charsize=1.2,alignment=0.5,/data
 
        fsc_text,freqMax*0.90,legendYStart-legendYSpace/2,"r = "+string(sliceRadii[i]/h2.h,format='(f3.1)'),$
                 charsize=1.2,alignment=0.5
      endfor    
    end_PS  
  
  endfor ;i
  
end

; fourierContourNonInt(): contour power in radius and non-integer frequency
;                         (unused)

pro fourierContourNonInt

  basePath = '/n/home07/dnelson/spirals/lambdacrit/sims/LC_10m_disk_6/'
  workingPath = '/n/home07/dnelson/spirals/lambdacrit/sims/LC_10m_disk_6/output3/'
  snapPath = basePath + "output/snap_"
  m=5

  ; config
  nRadialBins  = 120
  nAngularBins = 360
  
  minMaxRadius = [0.0,30.0]
  
  nModes = 7
  modeMax = 8.0
  
  modes = findgen(nModes)/nModes * (modeMax-1) + 1.0
  
  units = getUnits()

  ; arrays
  A   = fltarr(nRadialBins,nModes+1)
  phi = fltarr(nRadialBins,nModes+1)

  ; bins  
  dTheta = 2.0 * !pi / nAngularBins
  dR     = (minMaxRadius[1]-minMaxRadius[0]) / nRadialBins
  
  angularMidBins  = findgen(nAngularBins) / nAngularBins * 2.0 * !pi - !pi ;[-pi,+pi]
  angularMidBins += dTheta/2.0
  
  radialMidBins = findgen(nRadialBins) / nRadialBins * minMaxRadius[1] + minMaxRadius[0]
  radialMidBins += dR/2.0

  ; load snapshot
  h = loadSnapshot(snapPath,m,s_pos,s_vel,s_id,c_pos,c_vel,c_id)

  ; convert to polar coordinates
  r     = reform(sqrt(s_pos[0,*]^2.0 + s_pos[1,*]^2.0))
  arctan,s_pos[0,*],s_pos[1,*],theta_rad,theta_deg
  
  theta_rad -= !pi ;[-pi,+pi]

  ; bin data
  h2d = hist_2d(r,theta_rad,bin1=dR,bin2=dTheta,min1=minMaxRadius[0],min2=-!pi, $
                                                max1=minMaxRadius[1],max2=!pi)
  
  h2d = h2d[0:nRadialBins-1,0:nAngularBins-1]

  ; make Fourier transform
  for i=0,nRadialBins-1 do begin
    dd = reform(h2d[i,*])
    
    A[i,0] = dTheta/!pi * total(dd)
    
    for j=1,nModes do begin
      aj = dTheta/!pi * total(dd * cos(modes[j-1]*angularMidBins))
      bj = dTheta/!pi * total(dd * sin(modes[j-1]*angularMidBins))
      
      arctan,bj,aj,aTanVal,aTanValDeg
      
      A[i,j] = sqrt(aj^2.0 + bj^2.0)
      phi[i,j] = aTanVal
      ;print,j,aj,bj,A[i,j],phi[i,j]
    endfor
  endfor
  
  ; plot contour
  start_PS, workingPath+"fCont.1.eps"  
    fsc_plot,[0],[0],/nodata,xrange=minMaxRadius,yrange=[0,max(A[*,1:*])*1.1],charsize=1.2,$
      xtitle="radius [kpc]",ytitle="mode amplitude",/xs,/ys
    
    ; legend
    legendYStart = max(A[*,1:*])*1.05
    legendYSpace = max(A[*,1:*])*0.035
            
    for j=1,nModes do begin
      oplot,radialMidBins,reform(A[*,j]),color=fsc_color(colors[j mod n_elements(units.colors)]),thick=2.0
      oplot,[15,17],[legendYStart-legendYSpace*j,legendYStart-legendYSpace*j],$
            color=fsc_color(units.colors[j mod n_elements(units.colors)]),thick=2.0
      xyouts,18,legendYStart-0.2-legendYSpace*j,"m = "+str(modes[j-1]),charsize=1.2,alignment=0.5,/data
    endfor
  end_PS
  
  start_PS, workingPath+"fCont.2.eps"
    fsc_plot,[0],[0],/nodata,xrange=[0,max(modes)],yrange=[0,max(A[*,1:*])*1.1],$
      xtitle="mode",ytitle="mode amplitude at 2r_s",/xs,/ys
    
    targetRad = 8.0
    w = where(abs(radialMidBins - targetRad) eq min(abs(radialMidBins - targetRad)))
    
    fsc_plot,modes,A[w,*],/overplot
  end_PS
  
  start_PS, workingPath+"fCont.3.eps"
    Ac = A[*,1:*]
    fsc_contour,Ac,radialMidBins,modes,xtitle="radius [kpc]",ytitle="mode",$
      levels=max(Ac)*[0.3,0.5,0.7,0.9],c_charsize=0.4
  end_PS
  
end

; fourierModeAnalysis(): do binned Fourier mode analysis, based on Mark's code

function fourierModeAnalysis, snapPath, iSnap, epsOut=epsOut, plot=plot, old=old

  ; config
  nRadialBins  = 120 ;120
  nAngularBins = 180 ;360
  
  minMaxRadius = [0.0,30.0]
  
  ; fourier modes
  nModes = 8

  ; arrays
  A   = fltarr(nRadialBins,nModes+1)
  phi = fltarr(nRadialBins,nModes+1)
  aj  = fltarr(nRadialBins,nModes+1)
  bj  = fltarr(nRadialBins,nModes+1)
  
  ; bins  
  dTheta = 2.0 * !pi / nAngularBins
  dR     = (minMaxRadius[1]-minMaxRadius[0]) / nRadialBins
  
  angularMidBins  = findgen(nAngularBins) / nAngularBins * 2.0 * !pi - !pi ;[-pi,+pi]
  angularMidBins += dTheta/2.0
  
  radialMidBins = findgen(nRadialBins) / nRadialBins * minMaxRadius[1] + minMaxRadius[0]
  radialMidBins += dR/2.0

  ; load snapshot and convert to polar coordinates
  if not keyword_set(old) then begin
    xyz = loadStars(snapPath, iSnap)
    r     = reform(sqrt(xyz.x^2.0 + xyz.y^2.0))
    arctan,xyz.x,xyz.y,theta_rad,theta_deg
  endif else begin
    old = loadSnapshot(snapPath+"output/snap_",iSnap,s_pos,s_vel,s_id,c_pos,c_vel,c_id)
    r   = reform(sqrt(s_pos[0,*]^2.0 + s_pos[1,*]^2.0))
    arctan,s_pos[0,*],s_pos[1,*],theta_rad,theta_deg
  endelse
  
  theta_rad -= !pi ;[-pi,+pi]
  
  ; bin data
  h2d = hist_2d(r,theta_rad,bin1=dR,bin2=dTheta,min1=minMaxRadius[0],min2=-!pi, $
                                                max1=minMaxRadius[1],max2=!pi)
  
  h2d = h2d[0:nRadialBins-1,0:nAngularBins-1]
  
  ; make Fourier transform
  for i=0,nRadialBins-1 do begin
    dd = reform(h2d[i,*])
    
    A[i,0] = dTheta/!pi * total(dd)
    
    for j=1,nModes do begin
      aj[i,j] = dTheta/!pi * total(dd * cos(j*angularMidBins))
      bj[i,j] = dTheta/!pi * total(dd * sin(j*angularMidBins))
      
      arctan,bj[i,j],aj[i,j],aTanVal,aTanValDeg
      
      A[i,j] = sqrt(aj[i,j]^2.0 + bj[i,j]^2.0)
      phi[i,j] = aTanVal
      ;print,j,aj[i,j],bj[i,j],A[i,j],phi[i,j]
    endfor
  endfor
  
  ; plot
  if keyword_set(plot) then begin
    plotFourierModes, iSnap, radialMidBins, minMaxRadius, nModes, A, epsOut=epsOut
  endif
  
  ; set returns
  f = {r:radialMidBins,nModes:nModes,aj:aj,bj:bj,A:A,phi:phi}
  return,f
  
end

; fourierSequence():

pro fourierSequence, simName, old=old

  ; config
  basePath = '/n/home07/dnelson/spirals/lambdacrit/sims/
  filePath = basePath + simName + "/"

  nSnaps      = max([n_elements(file_search(filePath+"output/Stars_*_0")),$
                     n_elements(file_search(filePath+"output/Stars_*.bin")),$
                     n_elements(file_search(filePath+"output/snap_*"))])
  snapSpacing = 1

  ; calc fourier modes
  for j=0,nSnaps-1,1 do begin
    ; if pre-calculated fourier file already exists, skip
    if (file_test(filePath+"output2/fourierSaved_"+str(j+1)+".sav")) then begin
      print,'skipped ',j+1
      continue
    endif
  
    ; calculate power spectrum
    f = fourierModeAnalysis(filePath, (j+1)*snapSpacing, old=old)
 
    ; load MCs for time
    if not keyword_set(old) then begin
      mcs  = loadMCs(filePath,(j+1)*snapSpacing,h)
      time = h[1] * 1000.0 ;Myr
    endif else begin
      time = old.time * 1000.0 ;Myr
    endelse
    
    ; save fourier file
    fileName = filePath + "output2/fourierSaved_" + str(j+1) + ".sav"
    save,f,time,filename=fileName  
    
    print,j+1
  endfor

end

; plotFourierModes()

pro plotFourierModes, iSnap, radialMidBins, minMaxRadius, nModes, A, epsOut=epsOut

  units = getUnits()

  ; plot mode coefficients
  if keyword_set(epsOut) then begin
    start_PS, workingPath + 'fmodes_'+str(iSnap)+'.eps'
  endif

    fsc_plot,[0],[0],/nodata,xrange=minMaxRadius,yrange=[0,15], $
             xtickname=replicate(' ',10),ytickname=replicate(' ',10)
    
    ; legend
    legendYStart = 14
    legendYSpace = 0.7
    ; xyouts,17.5,legendYStart-legendYSpace*(nModes+4),"t = "+string(h.time*1000,format='(f5.1)')+" Myr",alignment=0.5
    
    for j=1,nModes do begin
      fsc_plot,radialMidBins,reform(A[*,j]),color=fsc_color(units.colors[j]),thick=2.0,/overplot
      fsc_plot,[13,15],[legendYStart-legendYSpace*j,legendYStart-legendYSpace*j],$
            color=fsc_color(units.colors[j]),thick=2.0,/overplot
      fsc_text,17,legendYStart-0.2-legendYSpace*j,"m = "+str(j),charsize=1.2,alignment=0.5,/data
    endfor

  if keyword_set(epsOut) then begin
    end_PS
  endif

end

