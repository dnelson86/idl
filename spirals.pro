; spirals.pro
; spirals project
; dnelson may.2011

@ helper
@ spiralsLoad
@ spiralsVis
@ spiralsMCs
@ spiralsICs
@ spiralsFourier
@ spiralsPatSpeed
@ spiralsLC

; fitExpProfile():

function fitExpProfile, r

  nRadBins = 180
  rMax     = 20.0
  dr       = rMax / nRadBins
  
  meanProfile = fltarr(nRadBins)
  rMid = findgen(nRadBins)*dr + dr/2.0
  
  for j=0,nRadBins-1 do begin
    w = where(r gt dr*j and r le dr*(j+1),count)
    if (count ne 0) then begin
      ;rMid = ( dr*j + dr*(j+1) ) / 2.0
      meanProfile[j] = count/(rMid[j]*dr*2*!pi)
    endif
  endfor
  
  meanFit = linfit(rMid,alog(meanProfile))
  ;fitRS = -1.0 / meanFit[1]
  
  return, meanFit
end

; stackDensity():

pro stackDensity, filePath, filePathMCs, workingPath, i

  ; config
  fileBase = workingPath+'sD_1sL_'+string(i,format='(I04)')
  units = getUnits()
  
  simParams = loadSimParams(filePath)
  
  nBins = 128
  xypMin = -0.5 * simParams.h
  xypMax = 0.5 * simParams.h
  binSize = (xypMax - xypMin) / nBins
  bins = (xypMin + binSize*findgen(nBins)) / simParams.h
  
  ; rs strip
  mc_r_min = 0.9*simParams.h ;0.8*2*sL
  mc_r_max = 1.1*simParams.h ;1.2*2*sL
  
  ;if (0) then begin
  h2d = fltarr(nBins,nBins)

  ; load MC positions
  mcs = loadMCs(filePathMCs,i,h)

  r_mcs     = reform(sqrt(mcs.x^2.0 + mcs.y^2.0))
  arctan,mcs.x,mcs.y,theta_rad_mcs,theta_deg_mcs

  ; histo stars
  xyz = loadStars(filePath, i)
  
  r     = reform(sqrt(xyz.x^2.0 + xyz.y^2.0))
  arctan,xyz.x,xyz.y,theta_rad,theta_deg
  
  ; fit axisymmetric surface density to remove exp profile
  meanFit = fitExpProfile(r)
  
  for j=0,n_elements(mcs)-1 do begin
  ;for j=0,100 do begin
  
    if (r_mcs[j] lt xypMax*1.5) then begin
      ;print,'skip smallr ',j
      continue
    endif
    if (r_mcs[j]*theta_rad_mcs[j] lt xypMax) then begin
      ;print,'warning1 ',j
      continue
    endif
    if (r_mcs[j]*theta_rad_mcs[j] ge (2*!pi*r_mcs[j] - xypMax)) then begin
      ;print,'warning2 ',j
      continue
    endif
    if (r_mcs[j] lt mc_r_min or r_mcs[j] gt mc_r_max) then begin
      ;print,'rs strip ',j
      continue
    endif
  
    v1 = r - r_mcs[j]
    v2 = r * (theta_rad - theta_rad_mcs[j])
    ;v3 = r_mcs[j] * (theta_rad - theta_rad_mcs[j])
  
    ; histogram local stars
    h2dt = hist_2d(v1,v2,bin1=binSize,bin2=binSize,$
                  min1=xypMin,min2=xypMin,max1=xypMax+binSize,max2=xypMax+binSize)
    h2dt = h2dt[0:nBins-1,0:nBins-1]
    
    ; subtract off exp (radial) profile
    for k=0,nBins-1 do begin
      rMid = r_mcs[j] + xypMin + binSize*k
      h2dt[k,*] -= exp(meanFit[0] + meanFit[1]*rMid)*binSize*binSize
    endfor

    h2d += h2dt
    ;print,j
  endfor

  save,h2d,filename=fileBase+".sav"
  ;endif
  ;restore,fileBase+".sav"
  h2draw = h2d
  
  ; rescale
  if (min(h2d) lt 0) then begin
    h2d += abs(min(h2d))+1
  endif
  w = where(h2d eq 0, count, complement=ww)
  if (count ne 0) then h2d[w] = min(h2d[ww])
  
  h2d = alog10(h2d)
  h2d = (h2d-min(h2d))*250.0 / (max(h2d)-min(h2d)) ;0-250
  ;h2d = filter_image(h2d,FWHM=[1.1,1.1],/ALL) ;gaussian kernel convolution
  h2d += 1.0 ;1-251
  
  PS_Start, FILENAME=fileBase+".4.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=4.0, ys=4.0, /inches

      xm = !x.margin
      ym = !y.margin
      !x.margin = [4,0]
      !y.margin = [0,0]
      
    ; color table
    loadct, 33, bottom=1, /silent ;ncolors=252            
            
    pMinMax = [xypMin,xypMax]
    tvim,h2d,pcharsize=1.0,xrange=pMinMax,yrange=pMinMax
      
    fsc_contour,filter_image(h2draw,FWHM=[nBins/40.0,nBins/40.0],/ALL),bins,bins,charsize=1.2,c_charsize=1.0,$
            xtitle="radial / rs",ytitle="azimuthal / rs",/xs,/ys,nlevel=5
      
    xp = findgen(nBins)/nBins * (xypMax-xypMin) + xypMin
    yp = fltarr(nBins)
    for j=0,nBins-1 do begin
      yp[j] = mean(h2draw[j,*])
    endfor
    fsc_plot,xp,yp,charsize=1.2,xtitle="x [kpc]",ytitle="radial profile"  
      
      !x.margin = xm
      !y.margin = ym
      
  PS_End

  ; plot
  PS_Start, FILENAME=fileBase+".rp2.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=4.0, ys=4.0, /inches
            
    xp = findgen(nBins)/nBins * (xypMax-xypMin) + xypMin
    yp = fltarr(nBins)
    for j=0,nBins-1 do begin
      yp[j] = mean(h2draw[j,*])
    endfor
    fsc_plot,xp,yp,charsize=1.2,xtitle="x [kpc]",ytitle="radial profile"  
                  
  PS_End, /PNG, /Delete_PS
  
  PS_Start, FILENAME=fileBase+".2.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=4.0, ys=4.0, /inches
            
    xm = !x.margin
    ym = !y.margin
    !x.margin = [0,0]
    !y.margin = [0,0]
      
    ; color table (TODO)
    loadct, 33, bottom=1, /silent ;ncolors=252            
            
    pMinMax = [xypMin,xypMax]
    tvim,h2d,xrange=pMinMax,yrange=pMinMax,pcharsize=1.0,/noaxis   
            
    !x.margin = xm
    !y.margin = ym            
            
  PS_End, /PNG, /Delete_PS

  PS_Start, FILENAME=fileBase+".c.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=6.0, ys=6.0, /inches
            
    contour,filter_image(h2draw,FWHM=[nBins/40.0,nBins/40.0],/ALL),bins,bins,charsize=1.2,c_charsize=1.0,$
            xtitle="radial / rs",ytitle="azimuthal / rs",/xs,/ys,nlevel=5
            
  PS_End, /PNG, /Delete_PS

end

; histoStars(): 2d histogram stars
;   inputs: snapPath, m
;   outputs: h2d

pro histoStars, snapPath, m, h2d, scaleLength, inc=inc, norm=norm, old=old

  units = getUnits()
  
  if not keyword_set(old) then begin
    ; new: load outputStarPositions dump
    xyz = loadStars(snapPath, m)
  
    x = reform(xyz.x) * units.UnitLength_in_cm / (units.pc_in_cm*1000) ;x, kpc
    y = reform(xyz.y) * units.UnitLength_in_cm / (units.pc_in_cm*1000) ;y, kpc
    z = reform(xyz.z) * units.UnitLength_in_cm / (units.pc_in_cm*1000) ;z, kpc
  endif else begin
    ; old: load full snapshot
    old = loadSnapshot(snapPath+"output/snap_",m,s_pos,s_vel,s_id,c_pos,c_vel,c_id)
    
    x = reform(s_pos[0,*]) * units.UnitLength_in_cm / (units.pc_in_cm*1000) ;x, kpc
    y = reform(s_pos[1,*]) * units.UnitLength_in_cm / (units.pc_in_cm*1000) ;y, kpc
    z = reform(s_pos[2,*]) * units.UnitLength_in_cm / (units.pc_in_cm*1000) ;y, kpc
  endelse

  ; histo config
  nBins = 512
  xypMin = -5.0 * scaleLength
  xypMax = 5.0 * scaleLength
  binSize = (xypMax - xypMin) / nBins
  
  ;xps = xypMin + binSize * findgen((xypMax-xypMin)/binSize + 1)
  ;yps = xypMin + binSize * findgen((xypMax-xypMin)/binSize + 1)
  
  ; incline by angle i (deg) about the x-axis (in y,z plane)
  if keyword_set(inc) then begin
    v1 = x
    v2 = y * cos(-1.0*inc*!dtor) + z * sin(-1.0*inc*!dtor)
  endif else begin
    v1 = x
    v2 = y
  endelse

  ; histogram stars
  h2d = hist_2d(v1,v2,bin1=binSize,bin2=binSize,$
                min1=xypMin,min2=xypMin,max1=xypMax+binSize,max2=xypMax+binSize)
  h2d = h2d[0:nBins-1,0:nBins-1]

  ; rescale
  w = where(h2d eq 0, count, complement=ww)
  if (count ne 0) then h2d[w] = min(h2d[ww])
  
  ;todo convert h2d to physical density
  
  h2d = alog10(h2d)
  h2d = (h2d-min(h2d))*250.0 / (max(h2d)-min(h2d)) ;0-250
	h2d = filter_image(h2d,FWHM=[1.1,1.1],/ALL) ;gaussian kernel convolution
  h2d += 1.0 ;1-251
  
end

; crossCorrRun3D(): do cross correlation for snapshots and save results

pro crossCorrRun3D, snapPath, workingPath

  ; config
  nSnaps = n_elements(file_search(snapPath+"*"))
  crossCorrRadius  = 0.1 ; code units, 0.1 = 100pc
  
  ; loop over each snapshot
  for i=0,nSnaps-1,1 do begin 

    ; load snapshot
    print,'loading snap ['+str(i)+']'
    h = loadSnapshot(snapPath,i,s_pos,s_vel,s_id,c_pos,c_vel,c_id)
  
    ; get number of cloud / star particles
    nStars  = h.nPartTot[2]
    nClouds = h.nPartTot[3]
    
    nNeighbors = intarr(nClouds)
  
    ; loop over each GMC particle
    print,format='($, " computing cross correlation:")'
    for j=0, nClouds-1 do begin
    
      ; find number of nearby neighbors
      distsSQ = (s_pos[0,*] - c_pos[0,j])^2.0 + $
                (s_pos[1,*] - c_pos[1,j])^2.0 + $
                (s_pos[2,*] - c_pos[2,j])^2.0
                
      w = where(distsSQ le crossCorrRadius^2.0, count)
                 
      nNeighbors[j] = count
      
      ; statusbar
      if (j mod 100 eq 0) then print,format='($, " ")' 
      if (j mod 10 eq 0) then print,format='($, ".")' 
  
    endfor
    print,' done'
  
    ; save cross corr
    saveFilename = workingPath+'cc3d_'+string(i,format='(i3.3)')+'.sav'
    save,h,nNeighbors,c_pos,crossCorrRadius,filename = saveFilename
  
    ; calculate single statistic
    corrCoeff = mean(nNeighbors)
    
  endfor
  
end

; crossCorrRun2D(): do cross correlation for snapshots and save results

pro crossCorrRun2D, snapPath, workingPath

  ; config
  nSnaps = n_elements(file_search(snapPath+"*"))
  crossCorrRadius  = 0.1 ; code units, 0.1 = 100pc
  
  ; loop over each snapshot
  for i=0,nSnaps-1,1 do begin 

    ; load snapshot
    print,'loading snap ['+str(i)+']'
    h = loadSnapshot(snapPath,i,s_pos,s_vel,s_id,c_pos,c_vel,c_id)
  
    ; get number of cloud / star particles
    nStars  = h.nPartTot[2]
    nClouds = h.nPartTot[3]
    
    nNeighbors = intarr(nClouds)
  
    ; loop over each GMC particle
    print,format='($, " computing cross correlation:")'
    for j=0, nClouds-1 do begin
    
      ; find number of nearby neighbors
      distsSQ = (s_pos[0,*] - c_pos[0,j])^2.0 + $
                (s_pos[1,*] - c_pos[1,j])^2.0 ;neglect z-component (thin disk assumption)
                
      w = where(distsSQ le crossCorrRadius^2.0, count)
                 
      nNeighbors[j] = count
      
      ; statusbar
      if (j mod 100 eq 0) then print,format='($, " ")' 
      if (j mod 10 eq 0) then print,format='($, ".")' 
  
    endfor
    print,' done'
  
    ; save cross corr
    saveFilename = workingPath+'cc2d_'+string(i,format='(i3.3)')+'.sav'
    save,h,nNeighbors,c_pos,crossCorrRadius,filename = saveFilename
  
    ; calculate single statistic
    corrCoeff = mean(nNeighbors)
    
  endfor
  
end

; contourCC2D(): contour crosscor in 2D
;   inputs:  workingPath, i
;   outputs: nNeighbors, xc, yc, B, xps, yps

pro contourCC2D, workingPath, i, nNeighbors, xc, yc, B, xps, yps, norm=norm

  units = getUnits()

  ; load t=0 crosscor for normalization
;  saveFilename = workingPath+'cc2d_'+string(0,format='(i3.3)')+'.sav'
;  restore,saveFilename
;  rc = sqrt(xc^2.0 + yc^2.0)

  ; load crosscor results (now includes c_pos, crossCorRadius)
  saveFilename = workingPath+'cc2d_'+string(i,format='(i3.3)')+'.sav'
  restore,saveFilename

  ; arrays
  xc = c_pos[0,*] * units.UnitLength_in_cm / (units.pc_in_cm*1000) ;kpc
  yc = c_pos[1,*] * units.UnitLength_in_cm / (units.pc_in_cm*1000) ;kpc
  rc = sqrt(xc^2.0 + yc^2.0)
  
  snapTime = h.time * units.UnitTime_in_s / units.s_in_Myr ;Myr  

  ; convert <N> to <rho>
  sph_rad = crossCorrRadius * units.UnitLength_in_cm / units.pc_in_cm ;pc
  mass_per_star = (h.massTable[2]) * (units.UnitMass_in_g/units.Msun_in_g) ;Msun
  
  nNeighbors *= mass_per_star / (4.0/3.0*!pi*(sph_rad)^3.0) ;Msun/pc^3
  
  ; normalize? - exponential profile
  if keyword_set(norm) then begin
    zH  = 3.13 ;stellar scale height (kpc)
    cc0 = 2.2  ;central value (guess)
    nNeighbors /= ( cc0 * exp(-rc / zH) )
  endif
  
  ; reform data for contour
  w = where(nNeighbors gt 0)
  nNeighbors = nNeighbors[w]
  xc         = xc[w]
  yc         = yc[w]
  
  ; fit minimum curvature spline surface
  nps = 128
  xypMin = -11.0
  xypMax = 11.0
  xps = findgen(nps)/nps * (xypMax-xypMin) + xypMin
  yps = findgen(nps)/nps * (xypMax-xypMin) + xypMin
  B = min_curve_surf(nNeighbors, xc, yc, xout=xps, yout=yps)
  
  ; keep surface positive
  w = where(B lt 0, count)
  if (count ne 0) then B[w] = 0.0
  
  ;rescale (loose physical units)
  ;B = alog10(B)
  B = (B-min(B))*2.5 / (max(B)-min(B)) ;0-2.5 (similar to pre-norm)
  
end

; ccBinRadial(): bin cross correlation in radial/annular bins

pro ccBinRadial, snapPath, workingPath, i, numBins, radMinMax, radBinC, corrCoeffRad

  units = getUnits()

  ; load crosscor results
  saveFilename = workingPath+'cc2d_'+string(i,format='(i3.3)')+'.sav'
  restore,saveFilename
  
  ; arrays
  xc = c_pos[0,*] * units.UnitLength_in_cm / (units.pc_in_cm*1000) ;kpc
  yc = c_pos[1,*] * units.UnitLength_in_cm / (units.pc_in_cm*1000) ;kpc
  rc = sqrt(xc^2.0 + yc^2.0)
  
  snapTime = h.time * units.UnitTime_in_s / units.s_in_Myr ;Myr

  ; radial binning config
  radBins = findgen(numBins+1)/numBins * radMinMax[1] + radMinMax[0]
  radBinC = radBins + 0.5 * (radMinMax[1]-radMinMax[0])/numBins
  
  ; bin radially
  corrCoeffRad  = fltarr(numBins)
  corrCoeffNorm = fltarr(numBins)
  
  for j=0, numBins-1,1 do begin
    w = where(rc ge radBins[j] and rc lt radBins[j+1], count)
  
    if (count ne 0.0) then begin
      corrCoeffRad[j] = mean(nNeighbors[w])
    endif
    
  endfor
  
  ; convert <N> to <rho>
  sph_rad = crossCorrRadius * units.UnitLength_in_cm / units.pc_in_cm ;pc
  mass_per_star = (h.massTable[2]) * (units.UnitMass_in_g/units.Msun_in_g) ;Msun
  
  corrCoeffRad *= mass_per_star / (4.0/3.0*!pi*(sph_rad)^3.0) ;Msun/pc^3
  
  ; normalize? exponential profile
  zH = 3.13                 ;stellar scale height (kpc)
  cc0 = 2.2                 ;central value (guess)
  
  corrCoeffRad /= ( cc0 * exp(-radBinC / zH) )
  
end

; batch
; -----

pro newFullAnalysis, simName, old=old

  smFrames,simName

  r = barSpeed(simName, old=old)
  
  patSpeedSeriesBar,simName
  ;patSpeedSeriesSpiral,simName
  
  barRadProfile,simName,old=old
  sliceAngleStability,simName
  
  fourierSequence,simName

end

pro smFrames, simName, old=old

  basePath = "/n/home07/dnelson/spirals/lambdacrit/sims/"

  filePath    = basePath + simName + "/"
  workingPath = basePath + simName + "/output2/"

  nSnaps      = max([n_elements(file_search(filePath+"output/Stars_*_0")),$
                     n_elements(file_search(filePath+"output/Stars_*.bin")),$
                     n_elements(file_search(filePath+"output/snap_*"))])

  ;inc = [0.0,30.0,45.0,60.0,75.0] ;inclination
  ;SA  = findgen(6)*20.0

  ;for m=1,nSnaps-1,1 do begin
  m=197
    singleMovie, filePath, workingPath, m, old=old
    ;singleXYVR,  filePath, workingPath, m, inc=inc, SA=SA
    ;singleXYRT, filePath, workingPath, m, old=old, /radPro
  ;endfor  

end

pro sbsFrames, simName1, simName2

  basePath = "/n/home07/dnelson/spirals/lambdacrit/sims/"

  filePath1    = basePath + simName1 + "/"
  filePath2    = basePath + simName2 + "/"
  workingPath = basePath + simName1 + "/output2/"
  
  nSnaps      = max([n_elements(file_search(filePath1+"output/Stars_*_0")),$
                     n_elements(file_search(filePath1+"output/Stars_*.bin")),$
                     n_elements(file_search(filePath1+"output/snap_*"))])

  for i=1,nSnaps,1 do begin
    ;compareTwoRows, filePath1, filePath2, workingPath, i
    ;rowZoomIn, filePath1, workingPath, i
    sideBySide, filePath1, filePath2, workingPath, i
    print,i
  endfor  
  
 ;ffmpeg command:
 ;module load viz/ffmpeg-r21113
 ;ffmpeg -f image2 -i MCs_%03d.png -vcodec libx264 -vpre default -crf 22 -threads 0 test6.mkv
 ;adjust crf to change quality, lower numbers = higher quality (sane range 18-28)
 ;use -vpre slow if possible  
 ;-r 15 should set 15fps instead of 30 (default?)
  
end

pro trFrames

  simName1 = 'LC_10m_disk_2'
  simName2 = 'LC_10m_disk_2_nomc'
  simName3 = 'LC_10m_disk_2_mcsoff'
  off3 = 20 ;snapshot numbering offset

  basePath = "/n/home07/dnelson/spirals/lambdacrit/sims/"

  filePath1    = basePath + simName1 + "/"
  filePath2    = basePath + simName2 + "/"
  filePath3    = basePath + simName3 + "/"
  workingPath  = basePath + simName1 + "/output2/"
  
  nSnaps      = n_elements(file_search(filePath1+"output/snap_*"))
                     
;  for i=off3+1,nSnaps-1,1 do begin
i=off3
    tripleRow, filePath1, filePath2, filePath3, off3, workingPath, i, /old
    print,i
;  endfor 
      
end

pro stackSeries, simName

  basePath    = "/n/home07/dnelson/spirals/lambdacrit/sims/"

  filePath    = basePath + simName + "/"
  filePathMCs = basePath + simName + "/"
  workingPath = basePath + simName + "/output3/"

  nSnaps      = max([n_elements(file_search(filePath+"output/Stars_*_0")),$
                     n_elements(file_search(filePath+"output/Stars_*.bin"))])

  ;for i=1,nSnaps,20 do begin
  for i=10,200,10 do begin
    stackDensity, filePath, filePathMCs, workingPath, i
    print,i
  endfor  

end
