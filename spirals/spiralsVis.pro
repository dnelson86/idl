; spiralsVis.pro
; visualization (mostly 2d or derived quantities with time)
; dnelson july.2011

; plotZZDot()

pro plotZZDot

  ; sim select
  filePath = "/n/home07/dnelson/spirals/lambdacrit/sims/LC_10m_disk_6/"
  
  snapPath    = filePath + "output/snap_"
  workingPath = filePath + "output3/"
  
  m = 4
  
  ; load snapshot
  h  = loadSnapshot(snapPath,m,s_pos,s_vel,s_id,c_pos,c_vel,c_id)

  ; subset selection
  frac = 0.01
  inds = lindgen(n_elements(s_id)*frac) / frac

  z    = reform(s_pos[2,inds])
  zdot = reform(s_vel[2,inds])
  
  ; all points
  ;z    = reform(s_pos[2,*])
  ;zdot = reform(s_vel[2,*])

  ; plot
  PS_Start, FILENAME=workingPath+"zzdot2."+str(m)+".eps", /nomatch, /quiet, font=1, $
            bits_per_pixel=8, color=1, /encapsulated, decomposed=0, xs=7.5, ys=5.0, /inches ;3/2  
    
    zm  = max([min(z),max(z)])*1.1
    zdm = max([min(zdot),max(zdot)])*1.1
    
    fsc_plot,[0],[0],/nodata,xrange=[-zm,zm],yrange=[-zdm,zdm],$
             charsize=1.2,xtitle="z [kpc]",ytitle="dz/dt [km/s]",/xs,/ys    
    
    fsc_plot,z,zdot,psym=3,/overplot
    
  PS_End, /PNG

end

; plotClouds(): 2D plot MCs (barebones, no axis labels or eps)

pro plotClouds, filePath, i, h

  units = getUnits()
  
  xyRange = [-20,20]
  lifeCycleWidth = 0.0035

  ; load
  mcs = loadMCs(filePath,i,h)
  
  xc1 = mcs.x * units.UnitLength_in_cm / (units.pc_in_cm*1000) ;kpc
  yc1 = mcs.y * units.UnitLength_in_cm / (units.pc_in_cm*1000) ;kpc
  
  ; plot
  fsc_plot,xc1,yc1,charsize=1.0,psym=9,symsize=0.5,xrange=xyRange,yrange=xyRange,/xs,/ys, $
       xtickname=replicate(' ',10),ytickname=replicate(' ',10)
       
  curTime = float(h[1])     
  w_d = where((mcs.timeborn + mcs.lifetime lt curTime+lifeCycleWidth) and $
              (mcs.timeborn lt curTime), count_d)
  w_n = where(mcs.timeborn gt curTime-lifeCycleWidth, count_n)

  if (count_n ne 0) then $
    fsc_plot,xc1[w_n],yc1[w_n],psym=9,symsize=0.5,color=fsc_color('green'),/overplot
  if (count_d ne 0) then $
    fsc_plot,xc1[w_d],yc1[w_d],psym=9,symsize=0.5,color=fsc_color('red'),/overplot
 
end

; visStarsCloudsCC: two plots, stellar density overlaid with clouds, and cross correlation

pro visStarsCloudsCC, snapPath, workingPath, i, title

  pMinMax = [-10,10] ;plot range

  ; load stellar density field
  histoStars, snapPath, i, h2d
  
  ; load crosscor contours
  contourCC2D, workingPath, i, nNeighbors, xc, yc, B, xps, yps
  
  ; ps
  set_plot,'PS'
  epsFilename = workingPath+'SCCC_'+string(i,format='(i3.3)')+'.eps'
  device,filename=epsFilename,bits_per_pixel=8,color=1,/encapsulated,decomposed=0,$
         xs=11,ys=8.5,/inches  

    ; set levels (for cc contour) and color table
    levels = [0.0,0.2,0.3,0.5,1.0,1.5,2.0]
    colors = indgen(n_elements(levels))*(252/n_elements(levels))+3
    loadct, 33, ncolors=252, bottom=3  

    ; plot stellar density byte image
    
    ;!p.position = [0.05,0.30,0.45,0.70]
    tvim,h2d,xrange=pMinMax,yrange=pMinMax,pcharsize=1.0,pos=[0.04,0.24,0.46,0.76]
         ;xtitle="x [kpc]",ytitle="y [kpc]"
         ;title="RUN_mc_pro log stellar density + clouds (snap "+str(i)+")"
                     
    ; overplot cloud points
    plotsym,0,/fill ;0=circle, 3=star, 4=tri, 8=sq
    oplot,xc,yc,psym=8,symsize=0.5
    
    ; title
    xyouts,0.48,0.83,title+": log stellar density + clouds | cross correlation (snap "+str(i)+")", $
           alignment=0.5,/normal,charsize=1.2    
    
    ; plot cross correlation contour
    ;!p.position=[0.55,0.30,0.95,0.70]
    ;contour, nNeighbors, xc, yc, /irregular, /fill, /isotropic, levels=levels, c_colors=colors, $
    ;         xrange=pMinMax, yrange=pMinMax, /xs, /ys, charsize=1.0, $
    ;         /noerase, pos=[0.46,0.24,0.98,0.76], ytickname=replicate(' ',10)
    contour, B, xps, yps, /fill, levels=levels, c_colors=colors, $
             xrange=pMinMax, yrange=pMinMax, /xs, /ys, $
             charsize=1.0, /noerase, pos=[0.52,0.24,0.94,0.76], ytickname=replicate(' ',10)
             ;xtitle="x [kpc]", ytitle="y [kpc]", title="RUN_mc_pro cc (snap "+str(i)+")"

  device,/close_file
  set_plot,'X'

end

; ccVis2D(): unsmoothed & smoothed 2D visualization of the cloud/star cross correlation

pro ccVis2D, snapPath, workingPath, i
  
  ; load crosscor contours
  contourCC2D, workingPath, i, nNeighbors, xc, yc, B, xps, yps
  
  ; contour plot
  set_plot,'PS'
  epsFilename = workingPath+'cc2d_'+string(i,format='(i3.3)')+'.eps'
  device,filename=epsFilename,bits_per_pixel=8,color=1,/landscape,/encapsulated,decomposed=0,$
         xs=10,ys=6,/inches

    ; set levels and color table
    levels = [0.0,0.2,0.3,0.5,1.0,1.5,2.0]
    colors = indgen(n_elements(levels))+3
    loadct, 33, ncolors=n_elements(levels), bottom=3

    ; contour
    !p.multi = [0,2,1]
    contour, nNeighbors, xc, yc, /irregular, /fill, /isotropic, levels=levels, c_colors=colors, $
             xrange=[-10,10], yrange=[-10,10], /xs, /ys, charsize=1.2, $
             xtitle="x [kpc]", ytitle="y [kpc]", title="RUN_mc_pro cc (snap "+str(i)+")"
    contour, B, xps, yps, /fill, /isotropic, levels=levels, c_colors=colors, $
             xrange=[-10,10], yrange=[-10,10], /xs, /ys, charsize=1.2, $
             xtitle="x [kpc]", ytitle="y [kpc]", title="RUN_mc_pro cc (snap "+str(i)+")"
    !p.multi = 0
    
  device,/close_file
  set_plot,'X'
  
end

; ccPlotRadial(): plot cross correlation as a function of radius for one snapshot

pro ccPlotRadial, snapPath, workingPath, i, title

  ; bin cc radially
  numBins = 40
  radMinMax = [0.0,20.0] ;kpc
  
  ccBinRadial, snapPath, workingPath, i, numBins, radMinMax, radBinC, corrCoeffRad
  
  ; exponential profile
  zH = 3.13                 ;stellar scale height (kpc)
  cc0 = 2.2                 ;central value (guess)
  rPts = findgen(200)/10.0
  ePts = cc0 * exp(-rPts / zH)
  
  ; plot corrCoeff vs. radius
  set_plot,'PS'
  epsFilename = workingPath+'ccrad_'+string(i,format='(i3.3)')+'.eps'
  device,filename=epsFilename,bits_per_pixel=8,color=1,/landscape,/encapsulated,/decomposed

    plot,radBinC,corrCoeffRad,charsize=1.5,psym=4, $ ;1=+, 2=star, 3=dot, 4=dia, 5=tri, 6=sq, 7=X
         xrange=radMinMax,yrange=[0,max(corrCoeffRad)*1.2],/xs,/ys, $
         xtitle="radius [kpc]",ytitle="Correlation Coefficient ("+textoidl("<\rho>/<\rho_0>")+")", $
         title=title+" Cross Correlation vs Radius (snap "+str(i)+")"

    ;oplot,rPts,ePts,line=2

  device,/close_file
  set_plot,'X'

end

; ccPlotRadialSeries(): plot cross correlation with radius over whole run, normed to t=0

pro ccPlotRadialSeries, snapPath, workingPath

  nSnaps = n_elements(file_search(snapPath+"*"))

  ; bin config
  numBins = 20
  radMinMax = [0.0,20.0] ;kpc
  
  ; arrays
  corrCoeffRad = fltarr(nSnaps, numBins)
  snapTimes    = fltarr(nSnaps)

  for i=0,nSnaps-1,1 do begin
    ; bin crosscorr radially
    ccBinRadial, snapPath, workingPath, i, radBinC, corrCoeffRad[i,*]
    
    print,str(i)+' done.'
  endfor
  
  ; normalize to t=0 snapshot (remove disk profile)
  for i=1,nSnaps-1,1 do begin
    corrCoeffRad[i,*] /= corrCoeffRad[0,*]
  endfor

  ; overplot all
  set_plot,'PS'
  device,filename=workingPath+'ccrad_series.eps',bits_per_pixel=8,color=1,/landscape,/encapsulated,/decomposed

    plot,radBinC,corrCoeffRad[0,*],charsize=1.5,psym=-4,/nodata, $ ;1=cross, 3=dot, 4=diamond
         xrange=[minRad,maxRad],yrange=[0,max(corrCoeffRad[1:*,*])*1.2],/xs,/ys, $
         xtitle="radius [kpc]",ytitle="Correlation Coefficient ("+textoidl("<\rho>/<\rho_0>")+")", $
         title="RUN_mc_pro Cross Correlation vs Radius (series)"
         
    for i=1,nSnaps-1,1 do begin
      oplot,radBinC,corrCoeffRad[i,*],line=1 ;0=solid, 1=dotted, 2=dashed
      ;xyouts
    endfor

  device,/close_file
  set_plot,'X'

end

; ccPlotTime(): plot cross correlation with time for a simulation run

pro ccPlotTime, snapPath, workingPath

  units = getUnits()

  ; config
  nSnaps = n_elements(file_search(snapPath+"*"))
  
  ; load results and plot over time
  corrCoeff = fltarr(nSnaps)
  snapTimes = fltarr(nSnaps)
  
  for i=0,nSnaps-1,1 do begin
    saveFilename = workingPath+'cc2d_'+string(i,format='(i3.3)')+'.sav'
    restore,saveFilename
    
    snapTimes[i] = h.time * units.UnitTime_in_s / units.s_in_Myr ;Myr
    
    ; compute single statistic
    corrCoeff[i] = mean(nNeighbors)
    
    ; convert <N> to <rho>
    sph_rad = crossCorrRadius * units.UnitLength_in_cm / units.pc_in_cm ;pc
    mass_per_star = (h.massTable[2]) * (units.UnitMass_in_g/units.Msun_in_g) ;Msun
    
    corrCoeff[i] *= mass_per_star / (4.0/3.0*!pi*(sph_rad)^3.0) ;Msun/pc^3
  endfor
  
  ; normalize corrCoeff by t=0 (~random) value
  corrCoeff /= corrCoeff[0]
  
  ; plot
  start_PS, workingPath+'cc_time.eps'
    fsc_plot,snapTimes,corrCoeff,xrange=minmax(snapTimes), yrange=[min(corrCoeff)*0.9,max(corrCoeff)*1.1], $
             xtitle="Time [Myr]",ytitle="Correlation Coefficient ("+textoidl("<\rho>/<\rho_0>")+")", $
             title="RUN_mc_pro"
         
    fsc_plot,snapTimes,corrCoeff,psym=4,/overplot
    fsc_plot,[0,1000],[corrCoeff[0],corrCoeff[0]],line=2,/overplot
  end_PS

end

; compareTwoRows()

pro compareTwoRows, filePath1, filePath2, workingPath, i

  ; config
  fileBase = workingPath+'tworows_'+string(i,format='(I04)')
  units = getUnits()
  
  pMinMax = [-10,10]
  minMaxRadius = [0.0,20.0]
  nModes  = 8
  
  ; histo stars
  histoStars, filePath1, i, h2d1, /norm
  histoStars, filePath2, i, h2d2, /norm

  ; plot
  PS_Start, FILENAME=fileBase+".eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=7.5, ys=5.0, /inches ;3/2

      !p.multi = [0,3,2]
      xm = !x.margin
      ym = !y.margin
      !x.margin = [0,0]
      !y.margin = [0,0]
      
      ; color table (TODO)
      loadct, 33, ncolors=252, bottom=3, /silent

      ; upper left (stars1)
      tvim,h2d1,xrange=pMinMax,yrange=pMinMax,pcharsize=1.0,/noaxis

      ; upper middle (mcs1)
      plotClouds, filePath1, i, header

      ; upper right (fourier1)
      ;f = fourierModeAnalysis(filePath1, i, plot=1)

      ; lower left (stars2)
      tvim,h2d1,xrange=pMinMax,yrange=pMinMax,pcharsize=1.0,/noaxis
      
      ; lower middle (mcs2)
      plotClouds, filePath2, i, header

      ; lower right (fourier2)
      ;f = fourierModeAnalysis(filePath2, i, plot=1)
      
        ; time ticker
        fsc_text,15,7,"t = " + header[1],alignment=0.5,charsize=1.2,/data
  
      !x.margin = xm
      !y.margin = ym
      !p.multi = 0
  
  PS_End, /PNG, /Delete_PS, Resize=68 ;PNG size=[xs,ys]*300*(resize/100)

end

; tripleRow()

pro tripleRow, filePath1, filePath2, filePath3, off3, workingPath, m, old=old

  ;skip if frame is already made
   if (file_test(workingPath+"tr_"+string(m,format='(I04)')+".png")) then begin
     print,'skipped ',m
     return
   endif

  ; config
  fileBase = workingPath+'tr_'+string(m,format='(I04)')
  units = getUnits()
  
  sp1      = loadSimParams(filePath1)
  sp2      = loadSimParams(filePath2)
  sp3      = loadSimParams(filePath3)
  
  pMinMax1 = [-5.0*sp1.h,5.0*sp1.h]
  pMinMax2 = [-5.0*sp2.h,5.0*sp2.h]
  pMinMax3 = [-5.0*sp3.h,5.0*sp3.h]

  ; histo stars
  histoStars, filePath1, m, h2d1, sp1.h, /norm, old=old
  print,'1 ',old.time
  histoStars, filePath2, m, h2d2, sp2.h, /norm, old=old
  print,'2 ',old.time
  histoStars, filePath3, m-off3, h2d3, sp3.h, /norm, old=old
  print,'3 ',old.time
  
  ; get header for this time
  if not keyword_set(old) then begin
    mcs  = loadMCs(filePath1,m,h)
    time = h[1]
  endif else begin
    time = old.time
  endelse  
  
  ; plot
  PS_Start, FILENAME=fileBase+".eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=7.5, ys=2.5, /inches ;3/2

      !p.multi = [0,3,1]
      xm = !x.margin
      ym = !y.margin
      !x.margin = [0,0]
      !y.margin = [0,0]
      
      ; color table
      loadct, 33, bottom=1, /silent ;ncolors=252

      ; left (stars1)
      tvim,h2d1,xrange=pMinMax1,yrange=pMinMax1,pcharsize=1.0,/noaxis
      ;fsc_plot,mcs1.x,mcs1.y,psym=16,symsize=0.2,color=fsc_color('purple'),/overplot

      ; center (stars2)
      tvim,h2d2,xrange=pMinMax2,yrange=pMinMax2,pcharsize=1.0,/noaxis
      ;fsc_plot,mcs2.x,mcs2.y,psym=16,symsize=0.2,color=fsc_color('purple'),/overplot
      
      ; right (stars3)
      tvim,h2d3,xrange=pMinMax3,yrange=pMinMax3,pcharsize=1.0,/noaxis
      ;fsc_plot,mcs3.x,mcs3.y,psym=16,symsize=0.2,color=fsc_color('purple'),/overplot      
      
      ; time ticker
      fsc_text,(pMinMax3[1]-pMinMax3[1]/4),(pMinMax3[0]-pMinMax3[0]/20),$
               "t = " + string(time,format='(f4.2)'),$
               alignment=0.5,color=fsc_color('white'),charsize=1.2,/data
  
      !x.margin = xm
      !y.margin = ym
      !p.multi = 0
  
  PS_End, /PNG, /Delete_PS, Resize=68 ;PNG size=[xs,ys]*300*(resize/100)

end

; sideBySide()

pro sideBySide, filePath1, filePath2, workingPath, i

  ;skip if frame is already made
   if (file_test(workingPath+"sbs_"+string(i,format='(I04)')+".png")) then begin
     print,'skipped ',i
     return
   endif

  ; config
  fileBase = workingPath+'sbs_'+string(i,format='(I04)')
  units = getUnits()
  
  sp1      = loadSimParams(filePath1)
  sp2      = loadSimParams(filePath2)
  
  pMinMax1 = [-5.0*sp1.h,5.0*sp1.h]
  pMinMax2 = [-5.0*sp2.h,5.0*sp2.h]

  ; histo stars
  histoStars, filePath1, i, h2d1, sp1.h, /norm
  histoStars, filePath2, i, h2d2, sp2.h, /norm
  
  ; get header for this time
  mcs1 = loadMCs(filePath1,i,h1)
  ;mcs2 = loadMCs(filePath2,i,h2)

  ; plot
  PS_Start, FILENAME=fileBase+".eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=4.0, /inches

      !p.multi = [0,2,1]
      xm = !x.margin
      ym = !y.margin
      !x.margin = [0,0]
      !y.margin = [0,0]
      
      ; color table
      loadct, 33, bottom=1, /silent ;ncolors=252

      ; left (stars1)
      tvim,h2d1,xrange=pMinMax1,yrange=pMinMax1,pcharsize=0.0001
      ;fsc_plot,mcs1.x,mcs1.y,psym=16,symsize=0.2,color=fsc_color('purple'),/overplot

      ; right (stars2)
      tvim,h2d2,xrange=pMinMax2,yrange=pMinMax2,pcharsize=0.0001
      ;fsc_plot,mcs2.x,mcs2.y,psym=16,symsize=0.2,color=fsc_color('purple'),/overplot
      
      ; time ticker
      fsc_text,(pMinMax2[1]-pMinMax2[1]/5),(pMinMax2[0]-pMinMax2[0]/20),"t = " + string(h1[1],format='(f6.4)'),alignment=0.5,$
               color=fsc_color('white'),charsize=1.2,/data
  
      !x.margin = xm
      !y.margin = ym
      !p.multi = 0
  
  PS_End, /PNG, /Delete_PS, Resize=68 ;PNG size=[xs,ys]*300*(resize/100)

end

; singleMovie()

pro singleMovie, filePath, workingPath, m, inc=inc, old=old

  ;skip if frame is already made
   if (file_test(workingPath+"sm_"+string(m,format='(I04)')+".png")) then begin
     print,'skipped ',m
     return
   endif

  ; config
  fileBase = workingPath+'sm_'+string(m,format='(I04)')
  units = getUnits()
  
  h2 = loadSimParams(filePath)
  
  pMinMax = [-5.0*h2.h,5.0*h2.h]

  ; histo stars
  histoStars, filePath, m, h2d, h2.h, /norm, inc=inc, old=old
  
  ; get header for this time
  if not keyword_set(old) then begin
    mcs  = loadMCs(filePath,m,h)
    time = h[1]
  endif else begin
    time = old.time
  endelse

  ; plot
  PS_Start, FILENAME=fileBase+".eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=4.0, ys=4.0, /inches

      xm = !x.margin
      ym = !y.margin
      !x.margin = [0,0]
      !y.margin = [0,0]
      
      ; color table
      loadct, 33, bottom=1, /silent ;ncolors=252 ;blue-red rainbow
      ;loadct, 9, ncolors=252, bottom=1, /silent ;green-white exponential
      ;loadct, 0, bottom=1, /silent ;bw linear

      ; (stars1)
      tvim,h2d,xrange=pMinMax,yrange=pMinMax,pcharsize=0.0001
      ;fsc_plot,mcs.x,mcs.y,psym=16,symsize=0.2,color=fsc_color('white'),/overplot
      
      ; important radii markers
      ;tvcircle,h2.h,0,0,color=fsc_color('white'),line=1,/data
      ;tvcircle,2*h2.h,0,0,color=fsc_color('white'),line=1,/data

      ; time ticker
      fsc_text,(pMinMax[1]-pMinMax[1]/5),(pMinMax[0]-pMinMax[0]/20),"t = " + string(time,format='(f6.4)'),$
               alignment=0.5,color=fsc_color('white'),charsize=1.2,/data
  
      !x.margin = xm
      !y.margin = ym
  
  PS_End;, /PNG, /Delete_PS, Resize=68 ;PNG size=[xs,ys]*300*(resize/100)

end

; singleXYVR(): single frame of XY projected surf dens and VR observed radial velocity, with inc i

pro singleXYVR, filePath, workingPath, m, inc=inc, SA=SA

  ;skip if frame is already made
   if (file_test(workingPath+"xyvr_"+string(m,format='(I04)')+".png")) then begin
     print,'skipped ',m
     return
   endif

  units = getUnits()
  h2 = loadSimParams(filePath)

  ; config
  fileBase = workingPath+'xyvr_'+string(m,format='(I04)')
  
  nBinsXY  = 512
  nBinsVR  = 128
  
  pMinMax  = [-5.0*h2.h,5.0*h2.h]

  binSizeXY = (pMinMax[1]-pMinMax[0]) / nBinsXY
  binSizeVR = (pMinMax[1]-pMinMax[0]) / nBinsVR
  
  vrMean = fltarr(nBinsVR,nBinsVR)
  vrDisp = fltarr(nBinsVR,nBinsVR)

  ; load full snapshot (vel field required)
  h  = loadSnapshot(filePath+"output/snap_",m,s_pos,s_vel,s_id,c_pos,c_vel,c_id)
  time = round(h.time*1000)

  x = reform(s_pos[0,*])
  y = reform(s_pos[1,*])
  z = reform(s_pos[2,*])
  
  vx = reform(s_vel[0,*])
  vy = reform(s_vel[1,*])
  vz = reform(s_vel[2,*])

  if not keyword_set(SA) then SA = [0.0]
  if not keyword_set(inc) then inc = [0.0]
  
  for k=0,n_elements(SA)-1 do begin

    ; rotate
    v1 = x  * cos(-1.0*SA[k]*!dtor) - y  * sin(-1.0*SA[k]*!dtor)
    v2 = y  * cos(-1.0*SA[k]*!dtor) + x  * sin(-1.0*SA[k]*!dtor)
    v3 = vy * cos(-1.0*SA[k]*!dtor) + vx * sin(-1.0*SA[k]*!dtor)
    
    for n=0,n_elements(inc)-1 do begin
    
      print,m,k,SA[k],n,inc[n]
    
      fileName = fileBase + '_sa_'  + string(k,format='(I03)') + $
                            '_inc_' + string(n,format='(I03)') + '.eps'
      
      ; incline
      v1 = v1
      v2 = v2 * cos(-1.0*inc[n]*!dtor) + z  * sin(-1.0*inc[n]*!dtor)
      ;v2 /= cos(-1.0*inc[n]*!dtor) ; correct for known inclination (zero error)
      
      if (inc[n] ne 0.0) then begin
        v3 = v3 * sin(-1.0*inc[n]*!dtor) + vz * cos(-1.0*inc[n]*!dtor) ;contour LOS/radial vel
        v3 /= sin(-1.0*inc[n]*!dtor) ; de-project for disk plane component
      endif else begin
        v3 = v3
      endelse
    
      ; histogram (x,y)
      h2xy = hist_2d(v1,v2,bin1=binSizeXY,bin2=binSizeXY,$
                     min1=pMinMax[0],min2=pMinMax[0],max1=pMinMax[1]+binSizeXY,max2=pMinMax[1]+binSizeXY)
      h2xy = h2xy[0:nBinsXY-1,0:nBinsXY-1]  
      
      ; rescale h2xy
      w = where(h2xy eq 0, count, complement=ww)
      if (count ne 0) then h2xy[w] = min(h2xy[ww])
      
      h2xy = alog10(h2xy)
      h2xy = filter_image(h2xy,FWHM=[1.1,1.1],/ALL) ;gaussian kernel convolution
      h2xy = (h2xy-min(h2xy))*252.0 / (max(h2xy)-min(h2xy)) + 1.0 ;1-253
      
      ; bin (v_y)
      for i=0,nBinsVR-1 do begin
    
        yMin = pMinMax[0] + i * binSizeVR
        yMax = yMin + binSizeVR
        
        ;select y subset
        w1 = where(v2 ge yMin and v2 lt yMax, count1)
        
        if (count1 eq 0) then continue
        
        x1 = v1[w1]
        vr = v3[w1]
        
        for j=0,nBinsVR-1 do begin
          xMin = pMinMax[0] + j * binSizeVR
          xMax = xMin + binSizeVR
          
          ;find stars in bin
          w = where(x1 ge xMin and x1 lt xMax, count)
          
          if (count ne 0) then begin
            vrMean[j,i]   = mean(vr[w])
            vrDisp[j,i]   = stddev(vr[w])
            
            if (not finite(vrDisp[j,i])) then vrDisp[j,i] = 0.0
          endif
        endfor ;j
      endfor ;i
      
      ; rescale velocity moment fields
      smoothFac = 64.0 ;10
      vrMean2 = filter_image(vrMean, FWHM=[nBinsVR/smoothFac,nBinsVR/smoothFac],/ALL)
      vrDisp2 = filter_image(vrDisp, FWHM=[nBinsVR/smoothFac,nBinsVR/smoothFac],/ALL)  
    
      ; plot
      start_PS, fileName, xs=8.0, ys=4.0
    
          !p.multi  = [0,2,1]
          xm        = !x.margin
          ym        = !y.margin
          !x.margin = [0,0]
          !y.margin = [0,0]
          
          ; color table
          loadct, 33, bottom=1, /silent ;bw linear=0 (33=blue-red)
    
          ; (stars1)
          tvim,h2xy,xrange=pMinMax,yrange=pMinMax,pcharsize=0.001,$
               position=[0.05,0.05,0.5,0.95],/c_map
          
          ; important radii markers
          tvellipse,h2.h,h2.h*cos(inc*!dtor),0,0,$
            color=fsc_color('white'),thick=1.0,line=2,npoints=240,/data
          tvellipse,2*h2.h,2*h2.h*cos(inc*!dtor),0,0,$
            color=fsc_color('white'),thick=1.0,line=2,npoints=240,/data
                
          ; time ticker
          fsc_text,(pMinMax[0]-pMinMax[0]/4),(pMinMax[1]-pMinMax[1]/7),string(time,format='(i4)')+" Myr",$
                   alignment=0.5,color=fsc_color('white'),charsize=1.4,/data
           
          ; v_r contour            
          fsc_contour,vrDisp2,xtickname=replicate(' ',10),ytickname=replicate(' ',10),c_charsize=0.0001,$
                xtitle="",ytitle="",/xs,/ys,nlevel=6,position=[0.5,0.05,0.95,0.95],$
                color=fsc_color('orange'),thick=1.0
                
          fsc_contour,vrMean2,xtickname=replicate(' ',10),ytickname=replicate(' ',10),c_charsize=0.0001,$
                xtitle="",ytitle="",/xs,/ys,nlevel=10,position=[0.5,0.05,0.95,0.95],/overplot
               
          !x.margin = xm
          !y.margin = ym
          !p.multi  = 0
      
      end_PS, pngResize=68
      
    endfor ;n,inc

  endfor ;k,SA

end

; singleXYRT()

pro singleXYRT, filePath, workingPath, m, old=old, radPro=radPro

  ;skip if frame is already made
   if (file_test(workingPath+"xyrt_"+string(m,format='(I04)')+".png")) then begin
     print,'skipped ',m
     return
   endif

  units = getUnits()
  h2 = loadSimParams(filePath)

  mass_per_star = 8.34e-07 * 1e10 ;LC_8_1m in msun = 8.34e-6, 10m=8.34e-7

  ; config
  fileBase = workingPath+'xyrt_'+string(m,format='(I04)')
  
  nBinsXY  = 512
  nBinsR   = 300 ;256
  nBinsAng = 300 ;180
  pMinMax  = [-5.0*h2.h,5.0*h2.h]
  
  binSizeXY  = (pMinMax[1]-pMinMax[0]) / nBinsXY
  binSizeR   = (pMinMax[1]) / nBinsR
  binSizeAng = 2 * !pi / nBinsAng

  ; load
  if not keyword_set(old) then begin
    ; new: load outputStarPositions dump
    xyz = loadStars(filePath, m)
    
    mcs  = loadMCs(filePath,m,h)
    time = round(h[1]*1000.0)
  
    v1 = reform(xyz.x) ;x, kpc
    v2 = reform(xyz.y) ;y, kpc  
  endif else begin
    ; old: load full snapshot
    h  = loadSnapshot(filePath+"output/snap_",m,s_pos,s_vel,s_id,c_pos,c_vel,c_id)
    time = round(h.time*1000.0)
  
    v1 = reform(s_pos[0,*]) ;x, kpc
    v2 = reform(s_pos[1,*]) ;y, kpc 
  endelse
  
  print,m,time

  r     = reform(sqrt(v1^2.0 + v2^2.0))
  arctan,v1,v2,theta_rad,theta_deg
  
  ; histogram (x,y)
  h2xy = hist_2d(v1,v2,bin1=binSizeXY,bin2=binSizeXY,$
                 min1=pMinMax[0],min2=pMinMax[0],max1=pMinMax[1]+binSizeXY,max2=pMinMax[1]+binSizeXY)
  h2xy = h2xy[0:nBinsXY-1,0:nBinsXY-1]  
  
  ; histogram (r)
  radPro = histogram(r,binsize=binSizeR,min=0.0,max=pMinMax[1]+binSizeR,locations=radPts)
  
  ; histogram (r,theta)
  h2rt = hist_2d(r,theta_rad,bin1=binSizeR,bin2=binSizeAng,$
                 min1=0.0,min2=0.0,max1=pMinMax[1]+binSizeR,max2=2*!pi)
  h2rt = h2rt[0:nBinsR-1,0:nBinsAng-1]  

  ; fit exp profile and subtract from h2rt
  meanCount = fltarr(nBinsR)
  rMid      = findgen(nBinsR)*binSizeR + binSizeR/2.0
  tMid      = findgen(nBinsAng)*binSizeAng + binSizeAng/2.0
  
  for j=0,nBinsR-1 do begin
    w = where(r gt binSizeR*j and r le binSizeR*(j+1),count)
    if (count ne 0) then begin
      meanCount[j] = float(count)/nBinsAng
    endif
  endfor
  
  densSmooth = rebin(meanCount,nBinsR,nBinsAng)
  
  h2rt -= densSmooth
  
  ; rescale h2xy
  w = where(h2xy eq 0, count, complement=ww)
  if (count ne 0) then h2xy[w] = min(h2xy[ww])
  
  h2xy = alog10(h2xy)
  h2xy = filter_image(h2xy,FWHM=[1.1,1.1],/ALL) ;gaussian kernel convolution
  h2xy = (h2xy-min(h2xy))*252.0 / (max(h2xy)-min(h2xy)) + 1.0 ;1-253
  
  ; rescale h2rt
  h2rt = (h2rt-min(h2rt))*252.0 / (max(h2rt)-min(h2rt)) ;0-252
  w = where(h2rt eq 0, count, complement=ww)
  if (count ne 0) then h2rt[w] = min(h2rt[ww])

  h2rt = alog10(h2rt)

  h2rt = filter_image(h2rt,FWHM=[1.1,1.1],/ALL) ;gaussian kernel convolution
  h2rt = (h2rt-min(h2rt))*252.0 / (max(h2rt)-min(h2rt)) + 1.0 ;1-253

  ; plot
  start_PS, fileBase+".eps", xs=8.0, ys=4.0

      !p.multi  = [0,2,1]
      xm        = !x.margin
      ym        = !y.margin
      !x.margin = [0,0]
      !y.margin = [0,0]
      
      ; color table
      loadct, 33, bottom=1, /silent ;bw linear=0 (33=blue-red)

      ; (stars1)
      tvim,h2xy,xrange=pMinMax,yrange=pMinMax,pcharsize=0.001,$
           position=[0.05,0.05,0.5,0.95],/c_map
      
      ; important radii markers
      tvcircle,h2.h,0,0,color=fsc_color('white'),line=0,thick=0.5,/data
      tvcircle,2*h2.h,0,0,color=fsc_color('white'),line=0,thick=0.5,/data
        
      ; time ticker
      fsc_text,(pMinMax[0]-pMinMax[0]/4),(pMinMax[0]-pMinMax[0]/10),string(time,format='(i4)')+" Myr",$
               alignment=0.5,color=fsc_color('white'),charsize=1.4,/data

      ; (r,theta) or radial profile
      if keyword_set(radPro) then begin
        annuli_areas = !pi * (((rMid+binSizeR/2.0)*1000)^2.0 - ((rMid-binSizeR/2.0)*1000)^2.0)
        meanCount = meanCount * nBinsAng * mass_per_star / annuli_areas
        radPro = radPro * mass_per_star / annuli_areas
        fsc_plot,rMid/h2.h,meanCount,/noerase,position=[0.50,0.05,0.95,0.95],charsize=0.5,ys=8,$
                 ytickname=replicate(' ',10),/xlog,/ylog,xrange=[0.1,5.0],/xs;,yrange=[1.0,1000.0]
        ;fsc_plot,(radPts+binSizeR/2.0)/h2.h,radPro,color=fsc_color('red'),/overplot
        fsc_axis,5.0,0.0,0.0,/yaxis,/ys,charsize=0.5
      endif else begin
        tvim,h2rt,xrange=[0,pMinMax[1]],yrange=[0,2*!pi],pcharsize=1.0,$
             position=[0.50,0.05,0.95,0.95],/noaxis,/c_map
      endelse
           
      !x.margin = xm
      !y.margin = ym
      !p.multi  = 0
  
  end_PS, pngResize=68  

end

; surfDensXYRT: (x,y) and (r,theta) side by side, overdensity tracking

pro surfDensXYRT

  basePath = "/n/home07/dnelson/spirals/lambdacrit/sims/"

  simName     = "LC_10m_disk_2"
  filePath    = basePath + simName + "/"
  workingPath = basePath + simName + "/output2/"
  
  ; 200 Myr start
  i = 4

  units = getUnits()
  h2 = loadSimParams(filePath)

  ; config
  fileBase = workingPath+'xyrt_'+string(i,format='(I04)')
  snapPath = filePath + "output/snap_"
  
  nBinsXY  = 512
  nBinsR   = 300 ;256
  nBinsAng = 300 ;180
  pMinMax  = [-5.0*h2.h,5.0*h2.h]
  
  binSizeXY  = (pMinMax[1]-pMinMax[0]) / nBinsXY
  binSizeR   = (pMinMax[1]) / nBinsR
  binSizeAng = 2 * !pi / nBinsAng

  ; load full snapshot (need star IDs)
  h  = loadSnapshot(snapPath,i,s_pos,s_vel,s_id,c_pos,c_vel,c_id)
  time = round(h.time*1000)
  
  v1 = reform(s_pos[0,*]) * units.UnitLength_in_cm / (units.pc_in_cm*1000) ;x, kpc
  v2 = reform(s_pos[1,*]) * units.UnitLength_in_cm / (units.pc_in_cm*1000) ;y, kpc
  
  r     = reform(sqrt(v1^2.0 + v2^2.0))
  arctan,v1,v2,theta_rad,theta_deg
  
  ; histogram (x,y)
  h2xy = hist_2d(v1,v2,bin1=binSizeXY,bin2=binSizeXY,$
                 min1=pMinMax[0],min2=pMinMax[0],max1=pMinMax[1]+binSizeXY,max2=pMinMax[1]+binSizeXY)
  h2xy = h2xy[0:nBinsXY-1,0:nBinsXY-1]  
  
  ; histogram (r,theta)
  h2rt = hist_2d(r,theta_rad,bin1=binSizeR,bin2=binSizeAng,$
                 min1=0.0,min2=0.0,max1=pMinMax[1]+binSizeR,max2=2*!pi)
  h2rt = h2rt[0:nBinsR-1,0:nBinsAng-1]  

  ; fit exp profile and subtract from h2rt
  meanCount = fltarr(nBinsR)
  rMid      = findgen(nBinsR)*binSizeR + binSizeR/2.0
  tMid      = findgen(nBinsAng)*binSizeAng + binSizeAng/2.0
  
  for j=0,nBinsR-1 do begin
    w = where(r gt binSizeR*j and r le binSizeR*(j+1),count)
    if (count ne 0) then begin
      meanCount[j] = float(count)/nBinsAng
    endif
  endfor
  
  densSmooth = rebin(meanCount,nBinsR,nBinsAng)
  
  h2rt -= densSmooth
  
  ; make overdensity selection
  thresh = max(h2rt)*0.35 ;0.35 ;0.8
  wC = where(h2rt ge thresh,count)
  print,count,float(count)/n_elements(h2rt)
  
  tM_max = -!pi/12.0
  tB_max = !pi*0.6
  
  tM_min = -!pi/12.0
  tB_min = !pi*0.4
  
  tR_min = 3.0 ;kpc
  
  wC2D = array_indices(h2rt,wC)
  wCxr = reform(rMid[wC2D[0,*]])
  wCyt = reform(tMid[wC2D[1,*]])
  
  w = where(wCyt le (wCxr*tM_max+tB_max) and wCyt ge (wCxr*tM_min+tB_min) and wCxr ge tR_min,count)
  print,count,float(count)/n_elements(h2rt)
  
  if (count ne 0) then begin
    wC   = wC[w]
    wC2D = wC2D[*,w]
    wCxr = wCxr[w]
    wCyt = wCyt[w]
  endif
    
  ; find s_id in selected pixels
  s_id_mask = intarr(n_elements(s_id))
  
  for i=0,n_elements(wC)-1 do begin
    r_min = wCxr[i]-binSizeR/2.0
    r_max = wCxr[i]+binSizeR/2.0
    t_min = wCyt[i]-binSizeAng/2.0
    t_max = wCyt[i]+binSizeAng/2.0
    
    w = where(r ge r_min and r lt r_max and theta_rad ge t_min and theta_rad lt t_max,count)
    ;print,r_min,r_max,t_min,t_max,count
    if (count ne 0) then $
      s_id_mask[w] = 1

  endfor
  
  w = where(s_id_mask eq 1,count)
  sel_id = s_id[w]
  print,'sel ',n_elements(sel_id),float(n_elements(sel_id))/n_elements(s_id)

  ; rescale h2xy
  w = where(h2xy eq 0, count, complement=ww)
  if (count ne 0) then h2xy[w] = min(h2xy[ww])
  
  h2xy = alog10(h2xy)
  h2xy = filter_image(h2xy,FWHM=[1.1,1.1],/ALL) ;gaussian kernel convolution
  h2xy = (h2xy-min(h2xy))*252.0 / (max(h2xy)-min(h2xy)) + 1.0 ;1-253
  
  ; rescale h2rt
  h2rt = (h2rt-min(h2rt))*252.0 / (max(h2rt)-min(h2rt)) ;0-252
  w = where(h2rt eq 0, count, complement=ww)
  if (count ne 0) then h2rt[w] = min(h2rt[ww])

  h2rt = alog10(h2rt)

  h2rt = filter_image(h2rt,FWHM=[1.1,1.1],/ALL) ;gaussian kernel convolution
  h2rt = (h2rt-min(h2rt))*252.0 / (max(h2rt)-min(h2rt)) + 1.0 ;1-253

  ; replace selected pixels with custom color table data value
  h2rt[wC] = 255.0

  ; plot
  start_PS, fileBase+".eps", xs=8.0, ys=4.0

      !p.multi  = [0,2,1]
      xm        = !x.margin
      ym        = !y.margin
      !x.margin = [0,0]
      !y.margin = [0,0]
      
      ; color table
      loadct, 33, bottom=1, /silent ;bw linear=0 (33=blue-red)

      ; (stars1)
      tvim,h2xy,xrange=pMinMax,yrange=pMinMax,pcharsize=0.001,$
           position=[0.05,0.05,0.5,0.95],/c_map
      
      ; important radii markers
      tvcircle,h2.h,0,0,color=fsc_color('white'),line=2,/data
      tvcircle,2*h2.h,0,0,color=fsc_color('white'),line=2,/data
            
      ; time ticker
      fsc_text,(pMinMax[0]-pMinMax[0]/4),(pMinMax[1]-pMinMax[1]/7),string(time,format='(i4)')+" Myr",$
               alignment=0.5,color=fsc_color('white'),charsize=1.4,/data
           
      ; custom color table
      tvlct,r,g,b,/get
      cS = 255 ;val=255
 
      r[cS] = 255.0 ;white
      g[cS] = 255.0
      b[cS] = 255.0
      
      tvlct,r,g,b
           
      tvim,h2rt,xrange=[0,pMinMax[1]],yrange=[0,2*!pi],pcharsize=1.0,$
           position=[0.50,0.05,0.95,0.95],/noaxis,/c_map
           
      ; important radii markers
      ;fsc_plot,[h2.h,h2.h],[0,2*!pi],line=2,color=255,/overplot
      ;fsc_plot,2*[h2.h,h2.h],[0,2*!pi],line=2,color=255,/overplot      

      ; debug
      ;rpts = [0.0,10.0]
      ;fsc_plot,rpts,tM_max*rpts + tB_max,line=1,/overplot
      ;fsc_plot,rpts,tM_min*rpts + tB_min,line=1,/overplot

      !x.margin = xm
      !y.margin = ym
      !p.multi  = 0
  
  end_PS, pngResize=68
stop
  ; i+1 snapshot (250 Myr), follow s_id
  i = 5
  
  fileBase = workingPath+'xyrt_'+string(i,format='(I04)')

  h  = loadSnapshot(snapPath,i,s_pos,s_vel,s_id,c_pos,c_vel,c_id)
  time = round(h.time*1000)
  
  v1 = reform(s_pos[0,*]) * units.UnitLength_in_cm / (units.pc_in_cm*1000) ;x, kpc
  v2 = reform(s_pos[1,*]) * units.UnitLength_in_cm / (units.pc_in_cm*1000) ;y, kpc
  
  r     = reform(sqrt(v1^2.0 + v2^2.0))
  arctan,v1,v2,theta_rad,theta_deg
  
  ; histogram (x,y)
  h2xy = hist_2d(v1,v2,bin1=binSizeXY,bin2=binSizeXY,$
                 min1=pMinMax[0],min2=pMinMax[0],max1=pMinMax[1]+binSizeXY,max2=pMinMax[1]+binSizeXY)
  h2xy = h2xy[0:nBinsXY-1,0:nBinsXY-1]  
  
  ; histogram (r,theta)
  h2rt = hist_2d(r,theta_rad,bin1=binSizeR,bin2=binSizeAng,$
                 min1=0.0,min2=0.0,max1=pMinMax[1]+binSizeR,max2=2*!pi)
  h2rt = h2rt[0:nBinsR-1,0:nBinsAng-1]  

  ; fit exp profile and subtract from h2rt
  meanCount = fltarr(nBinsR)
  rMid      = findgen(nBinsR)*binSizeR + binSizeR/2.0
  tMid      = findgen(nBinsAng)*binSizeAng + binSizeAng/2.0
  
  for j=0,nBinsR-1 do begin
    w = where(r gt binSizeR*j and r le binSizeR*(j+1),count)
    if (count ne 0) then begin
      meanCount[j] = float(count)/nBinsAng
    endif
  endfor
  
  densSmooth = rebin(meanCount,nBinsR,nBinsAng)
  
  h2rt -= densSmooth

  ; take sel_id and choose wC
  inds = find_elements(s_id,sel_id)
  
  h2rt_sel = hist_2d(r[inds],theta_rad[inds],bin1=binSizeR,bin2=binSizeAng,$
                 min1=0.0,min2=0.0,max1=pMinMax[1]+binSizeR,max2=2*!pi)
  h2rt_sel = h2rt_sel[0:nBinsR-1,0:nBinsAng-1] 
  
  w = where(h2rt_sel ne 0,count)
  if (count ne 0) then $
    wC = w
  
  ; rescale h2xy
  w = where(h2xy eq 0, count, complement=ww)
  if (count ne 0) then h2xy[w] = min(h2xy[ww])
  
  h2xy = alog10(h2xy)
  h2xy = filter_image(h2xy,FWHM=[1.1,1.1],/ALL) ;gaussian kernel convolution
  h2xy = (h2xy-min(h2xy))*252.0 / (max(h2xy)-min(h2xy)) + 1.0 ;1-253
  
  ; rescale h2rt
  h2rt = (h2rt-min(h2rt))*252.0 / (max(h2rt)-min(h2rt)) ;0-252
  w = where(h2rt eq 0, count, complement=ww)
  if (count ne 0) then h2rt[w] = min(h2rt[ww])

  h2rt = alog10(h2rt)

  h2rt = filter_image(h2rt,FWHM=[1.1,1.1],/ALL) ;gaussian kernel convolution
  h2rt = (h2rt-min(h2rt))*252.0 / (max(h2rt)-min(h2rt)) + 1.0 ;1-253

  ; replace selected pixels with custom color table data value
  h2rt[wC] = 255.0

  ; plot
  start_PS, fileBase+".eps", xs=8.0, ys=4.0

      !p.multi  = [0,2,1]
      xm        = !x.margin
      ym        = !y.margin
      !x.margin = [0,0]
      !y.margin = [0,0]
      
      ; color table
      loadct, 33, bottom=1, /silent ;bw linear=0 (33=blue-red)

      ; (stars1)
      tvim,h2xy,xrange=pMinMax,yrange=pMinMax,pcharsize=0.001,$
           position=[0.05,0.05,0.5,0.95],/c_map
      
      ; important radii markers
      tvcircle,h2.h,0,0,color=fsc_color('white'),line=2,/data
      tvcircle,2*h2.h,0,0,color=fsc_color('white'),line=2,/data
            
      ; time ticker
      fsc_text,(pMinMax[0]-pMinMax[0]/4),(pMinMax[1]-pMinMax[1]/7),string(time,format='(i4)')+" Myr",$
               alignment=0.5,color=fsc_color('white'),charsize=1.4,/data
           
      ; custom color table
      tvlct,r,g,b,/get
      cS = 255 ;val=255
 
      r[cS] = 255.0 ;white
      g[cS] = 255.0
      b[cS] = 255.0
      
      tvlct,r,g,b
           
      tvim,h2rt,xrange=[0,pMinMax[1]],yrange=[0,2*!pi],pcharsize=1.0,$
           position=[0.50,0.05,0.95,0.95],/noaxis,/c_map
           
      !x.margin = xm
      !y.margin = ym
      !p.multi  = 0
  
  end_PS, pngResize=68

  stop
end

; rowZoomIn()

pro rowZoomIn, filePath, workingPath, m

  ; config
  fileBase = workingPath+'rowzoomin_'+string(m,format='(I04)')
  units = getUnits()
  
  pMinMax = [-10,10]
  
  zoom1 = [2,7]
  zoom2 = [4,5]

  ; histo stars
  histoStars, filePath, m, h2d, /norm
  histoStars, filePath, m, i=90.0, h2dEdgeOn
  
  ; load stars
  xyz = loadStars(filePath, m)
  
  ; load mcs (for time)
  mcs = loadMCs(filePath,m,header)

  ; plot
  PS_Start, FILENAME=fileBase+".eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=3.5, /inches ;3/1 2.667

      ;!p.multi = [0,3,1]
      xm = !x.margin
      ym = !y.margin
      !x.margin = [0,0]
      !y.margin = [0,0]
      
      ; color table (TODO)
      loadct, 33, ncolors=252, bottom=3, /silent ;3,252

      ; left (stars1 no zoom)
      tvim,h2d,xrange=pMinMax,yrange=pMinMax,pcharsize=1.0,/noaxis,position=[0.0,0.25,0.33,1.0]

        ; zoom boxes
        ;fsc_plot,[zoom1[0],zoom1[1],zoom1[1],zoom1[0],zoom1[0]], $
        ;         [zoom1[0],zoom1[0],zoom1[1],zoom1[1],zoom1[0]], $
        ;         line=0, color=fsc_color('black'), /overplot

      ; middle (stars zoom-middle)
      fsc_plot,xyz.x,xyz.y,psym=9,symsize=0.25,xrange=zoom1,yrange=zoom1,/xs,/ys, $
        xtickname=replicate(' ',10),ytickname=replicate(' ',10),xticklen=0,yticklen=0, $
        position=[0.33,0.25,0.66,1.0],/noerase
        
        ;fsc_plot,[zoom2[0],zoom2[1],zoom2[1],zoom2[0],zoom2[0]], $
        ;         [zoom2[0],zoom2[0],zoom2[1],zoom2[1],zoom2[0]], $
        ;         line=0, color=fsc_color('black'), /overplot

      ; right (stars zoom-max)
      fsc_plot,xyz.x,xyz.y,psym=9,symsize=0.25,xrange=zoom1,yrange=zoom2,/xs,/ys, $
        xtickname=replicate(' ',10),ytickname=replicate(' ',10),xticklen=0,yticklen=0, $
        position=[0.66,0.25,1.0,1.0],/noerase
      
      ; bottom (edge on star histo)
      tvimage,h2dEdgeOn,position=[0.0,0.0,1.0,0.25];,/keep_aspect_ratio
      ;fsc_plot,[0],[0],xrange=[-10,10],yrange=[-10,10],/nodata,xtickname=replicate(' ',10),ytickname=replicate(' ',10), $
      ;  position=[0.0,0.0,1.0,0.25],/noerase
      
        ; time ticker
        fsc_text,0.95,0.02,"t = " + string(header[1],format='(f6.4)'),$
                 alignment=0.5,charsize=1.2,color=fsc_color('white'),/normal
        
      !x.margin = xm
      !y.margin = ym
      !p.multi = 0
  
  PS_End, /PNG, Resize=80;,/Delete_PS ;PNG size=[xs,ys]*300*(resize/100)

end
