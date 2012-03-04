; spiralsMCs.pro
; MC related analysis
; dnelson feb.2011

function makenlog,xlo,xhi,n
  xlhi = alog(xhi)
  xllo = alog(xlo)
  xst = (xlhi-xllo)/float(n-1)    ; step size.
 
  return, exp(xllo+xst*findgen(n))
end

; plotMCs()
; plot output from outputMolClouds() for one timestep

pro plotMCs, filePath, workingPath, i, iOut

  pts = loadMCs(filePath, i, header)
  
  fileBase = workingPath+'MCs_'+string(iOut,format='(I04)')
  
  PS_Start, FILENAME=fileBase+".eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,/decomposed, xs=8, ys=8, /inches

      fsc_plot,pts.x,pts.y,charsize=1.5,psym=9, xrange=[-20,20],yrange=[-20,20],/xs,/ys, $
           xtitle="x [kpc]",ytitle="y [kpc]",title="MCs t="+string(header[1],format='(f6.4)'), $
           position=[0.1,0.1,0.95,0.95]
  
  PS_End, /PNG, /NoFix, /Delete_PS

end

; findPoissonBins()

function findPoissonBins, ICsPath, nRadialBins, nAngularBins, minMaxRadius

  ; load ICs
  h = loadSnapshot(ICsPath,'none',s_pos,s_vel,s_id,c_pos,c_vel,c_id)
  
  r     = reform(sqrt(s_pos[0,*]^2.0 + s_pos[1,*]^2.0))
  rSort = r[sort(r)]
  
  print,'minmax radius: ',minmax(r)
  
  nStarsTot = n_elements(s_id)
  ;errMax    = 0.02 ;5%
  ;numPerBin = 1/errMax^2.0
  numPerBin = nStarsTot / (nRadialBins+1)
 
  radialBins = fltarr(nRadialBins+1)
  
  for i=1,nRadialBins-1 do begin
    radialBins[i] = round(rSort[numPerBin*(i+1)-1]*100)/100.0
  endfor
  
  radialBins[nRadialBins] = minMaxRadius[1]
  
  return,radialBins
end

; findBins()
; misc code not used

pro findBins

  ; create bins  
  
  ;LINEAR
  ;------
  
  ;for (k=0; k <= NUM_RADIAL_BINS; k++)
  ;  radialBins[k] = MIN_RADIUS + (MAX_RADIUS - MIN_RADIUS) / NUM_RADIAL_BINS * k;
  ;for (m=0; m <= NUM_ANGULAR_BINS; m++)
  ;  angularBins[m] = 360.0 / NUM_ANGULAR_BINS * m;
  ;  

  radialBins = fltarr(nRadialBins+1)
  angularBins = fltarr(nAngularBins+1)
  
  r0 = 5
  p0 = 2.0
  radialBins[0] = 0.0
  
  for k=1,nRadialBins do begin
    radialBins[k] = minMaxRadius[0] + (minMaxRadius[1]-minMaxRadius[0]) / nRadialBins * k
    radialBins[k] *= (radialBins[k]/r0)^p0 / (minMaxRadius[1])^p0
  endfor
  for m=0,nAngularBins do begin
    angularBins[m] = 360.0 / nAngularBins * m;
  endfor
  
  radialBins /= max(radialBins)/minMaxRadius[1]
  
  ;LOGARITHMIC
  ;-----------
  xst = (alog10(50.0) - alog10(1.0))/float(10-1)    ; step size.
 
  ;radialBins = exp(1.0+xst*findgen(10))

  ;for k=0,nRadialBins do begin
  ;  radialBins[k] = 10.0^(alog10(1.0) + k*xst)
  ;endfor
  
end


; checkPoissonBins()

pro checkPoissonBins, ICsPath, radialBins, angularBins

  ; load ICs
  h = loadSnapshot(ICsPath,'none',s_pos,s_vel,s_id,c_pos,c_vel,c_id)
  
  r     = reform(sqrt(s_pos[0,*]^2.0 + s_pos[1,*]^2.0))
  arctan,s_pos[0,*],s_pos[1,*],theta_rad,theta_deg
  
  print,'minmax theta:  ',minmax(theta_deg)
  
  starHisto = fltarr(n_elements(radialBins)-1,n_elements(angularBins)-1)
  
  ;stop
  
  ; do binning
  for i=0,n_elements(radialBins)-2 do begin
    for j=0,n_elements(angularBins)-2,10 do begin
      radMin = radialBins[i]
      radMax = radialBins[i+1]
      
      thetaMin = angularBins[j]
      thetaMax = angularBins[j+1]
      
      ;print,i,j,radMin,radMax,thetaMin,thetaMax
      
      w = where(r gt radMin and r le radMax and theta_deg gt thetaMin and theta_deg le thetaMax, count)
      starHisto[i,j] = count
      
      ; check Poisson level
      print,radMin,radMax,thetaMin,thetaMax,count,1.0/sqrt(count)
  
    endfor
  endfor
  
  rbString = ''
  for i=0,n_elements(radialBins)-1 do begin
    rbString = rbString + string(radialBins[i],format='(F5.2, $)') + ","
  endfor
  
  print,'radialBins: ', rbString
  print,'angularBins: ',angularBins
  
  ;verify counting
  print,'total stars: ',n_elements(s_id)
  print,'total count: ',total(starHisto)
  
  ;errors
  print,'minmax mean stddev errors:',minmax(1.0/sqrt(starHisto)),mean(1.0/sqrt(starHisto)),stddev(1.0/sqrt(starHisto))

  stop
end