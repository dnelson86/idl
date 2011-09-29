; coldflowsUtil.pro
; cold flows - utility functions
; dnelson sep.2011

; redshiftToSnapNum(): convert redshift to the nearest snapshot number

function redshiftToSnapNum, redshift, verbose=verbose

  if not keyword_set(verbose) then verbose = 0
  
  gadgetPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/128_20Mpc/Gadget/output/'
  saveFileName = '/n/home07/dnelson/coldflows/snapnum.redshift.sav'

  if not (file_test(saveFileName)) then begin

    nSnaps = n_elements(file_search(gadgetPath+'snapdir_*'))
    
    redshifts = fltarr(nSnaps)
    times     = fltarr(nSnaps)
    
    for m=0,nSnaps-1 do begin
      ; format filename
      ext = str(string(m,format='(I3.3)'))
      f = gadgetPath + 'snapdir_' + ext + '/snap_' + ext + '.0.hdf5'
    
      ; load hdf5 header and save time+redshift
      fileID   = h5f_open(f)
      headerID = h5g_open(fileID,"Header")
      
      redshifts[m] = h5a_read(h5a_open_name(headerID,"Redshift"))
      times[m]     = h5a_read(h5a_open_name(headerID,"Time"))
      
      h5g_close, headerID
      h5f_close, fileID
    endfor
  
    ; save/restore
    save,gadgetPath,nSnaps,redshifts,times,filename=saveFileName
  endif else begin
    restore,saveFileName
  endelse
  
  ; for redshift zero hard select #314
  if (redshift eq 0.0) then begin
    snapNum = 314
  endif else begin
    w = where(abs(redshifts - redshift) eq min(abs(redshifts - redshift)),count)
    if (count eq 2) then w = w[0]
    if (count eq 0) then return,-1
    snapNum = w[0]
  endelse

  if (verbose) then $
    print,'Found nearest snapshot to z = ' + str(redshift) + ' (num = ' + str(snapNum) + ')'

  return,snapNum
end

; snapNumToRedshift(): convert snapshot number to redshift or time

function snapNumToRedshift, snapNum=snapNum, time=time, all=all

  saveFileName = '/n/home07/dnelson/coldflows/snapnum.redshift.sav'
  
  if not (file_test(saveFileName)) then begin
    print,'ERROR: Save missing.'
    return,0
  endif
  
  if not keyword_set(snapNum) then snapNum = -1
  
  ; restore
  restore,saveFilename

  if (not keyword_set(time)) then begin
    if (snapNum ge 0 and snapNum lt n_elements(redshifts)) then $
      return,redshifts[snapNum]
      
    if (keyword_set(all)) then $
      return,redshifts
  endif else begin
    if (snapNum ge 0 and snapNum lt n_elements(redshifts)) then $
      return,times[snapNum]
      
    if (keyword_set(all)) then $
      return,times
  endelse

end

; convertUtoTemp():

function convertUtoTemp, u, nelec, gamma=gamma, hmassfrac=hmassfrac

  units = getUnits()
  
  ; adiabatic index and hydrogen mass fraction defaults (valid for ComparisonProject)
  if not keyword_set(gamma)     then gamma = 5.0/3.0
  if not keyword_set(hmassfrac) then hmassfrac = 0.76
  
  ; calculate mean molecular weight
  meanmolwt = 4.0/(1.0 + 3.0 * hmassfrac + 4.0* hmassfrac * nelec) * units.mass_proton
  
  ; calculate temperature
  temp = (gamma-1.0) * u / units.boltzmann * units.UnitEnergy_in_cgs / units.UnitMass_in_g * meanmolwt
  
  return, temp
end

; convertCoolingRatetoCGS():

function convertCoolingRatetoCGS, coolrate, h=h

  units = getUnits()
  
  ; default little h
  if not keyword_set(h) then h = 0.7
  
  ; convert code units (du/dt) to erg/s/g (cgs)
  coolrate_cgs = coolrate * units.UnitEnergy_in_cgs * units.UnitTime_in_s^(-1.0) * $
                 units.UnitMass_in_g^(-1.0) * h
                 
  return, coolrate_cgs
end

; calcEntropy():

function calcEntropy, u, dens, gamma=gamma

  ; adiabatic index default (valid for ComparisonProject)
  if not keyword_set(gamma)     then gamma = 5.0/3.0
  
  pressure = (gamma-1.0) * u * dens
  entropy  = pressure / (dens^gamma)
  
  return, entropy
end

; rhoRatioToCrit():

function rhoRatioToCrit, rho, omega_b=omega_b

  units = getUnits()
  
  ; default omega_b (valid for ComparisonProject)
  if not keyword_set(omega_b) then omega_b = 0.044
  
  rho_b = omega_b * units.rhoCrit
  
  return, rho/rho_b

end

; redshiftToAge(): convert redshift to age of the universe (approximate)
function dtdz, z, lambda0 = lambda0, q0 = q0
  term1 = (1.0d + z)
  term2 = 2.0d * (q0 + lambda0) * z + 1.0d - lambda0
  term3 = (1.0d + z) * (1.0d +z)
  return, 1.0 / (term1 * sqrt(term2 * term3 + lambda0))
end
   
function redshiftToAge, z

  units = getUnits()
  
  ; config
  zform = 1000.0
  H0 = 70.0
  k = 0.0
  Omega_m = 0.27
  Lambda0 = 0.73
  q0 = -0.55
  
  ; arrays
  nz  = N_elements(z)
  age = z * 0.0
  
  ; integrate with qsimp
  for i= 0L, nz-1 do begin
    if (z[i] ge zform) then age_z = 0 else $
        qsimp,'dtdz', z[i], zform, age_z, q0 = q0, lambda0 = lambda0
    age[i] = age_z
  endfor

  return, age * 3.085678e+19 / 3.15567e+7 / H0 / 1e9 ;Gyr

end
  
; rhoTHisto(): make mass-weighted density-temperature 2d histogram

function rhoTHisto, dens_in, temp_in, mass=mass, nbins=nbins, plot=plot

  units = getUnits()

  ; config
  if not keyword_set(nbins) then nbins = 20
  dens = alog10(rhoRatioToCrit(dens_in))
  temp = alog10(temp_in)
  rMinMax = [-2.0,8.0]
  tMinMax = [3.0,7.0]
  
  ; calculate bin sizes
  binSizeRho  = (rMinMax[1]-rMinMax[0]) / nbins
  binSizeTemp = (tMinMax[1]-tMinMax[0]) / nbins
  
  if not keyword_set(mass) then begin
    ; hist_2d (no weighting)
    h2rt = hist_2d(dens,temp,bin1=binSizeRho,bin2=binSizeTemp,$
                   min1=rMinMax[0],min2=tMinMax[0],max1=rMinMax[1]+binSizeRho,max2=tMinMax[1]+binSizeTemp)
    h2rt /= float(max(h2rt)) ;norm s.t. colorbar gives fraction wrt max cell
  endif else begin
    ; hist2d (weighted)
    h2rt = hist2d(dens,temp,mass*(units.UnitMass_in_g/units.Msun_in_g),$
                  binsize1=binsizeRho,binsize2=binSizeTemp,$
                  min1=rMinMax[0],min2=tMinMax[0],max1=rMinMax[1]+binSizeRho,max2=tMinMax[1]+binSizeTemp)
  endelse

  ; plot
  if keyword_set(plot) then begin
    ; color table
    loadct, 2, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
      
    tvim,h2rt,pcharsize=!p.charsize-1.0,scale=1,clip=[10,100],$;,/c_map
         xtitle="log ("+textoidl("\rho / \rho_{crit}")+")",ytitle="log (T [K])",$
         stitle="Total Mass (Msun)",barwidth=0.5,lcharsize=!p.charsize-1.5,$
         xrange=[-2.0,8.0],yrange=[3.0,7.0],$;xrange=rMinMax,yrange=tMinMax,$
         /rct;,nodata=0,rgb_nodata=[1.0,1.0,1.0] ;display zeros as white not black
  endif
  
  return,h2rt

end

; sphDensityProjection(): make density projection using SPH kernel (inspired by Mark's sphMap)

function sphDensityProjection, pos, hsml, mass, quantity=quantity, imgSize=imgSize, boxSize=boxSize,$
                               boxCen=boxCen, axis0=axis0, axis1=axis1, mode=mode, periodic=periodic,$
                               verbose=verbose

  ; config
  if not keyword_set(axis0) then axis0 = 0
  if not keyword_set(axis1) then axis1 = 1
  if not keyword_set(verboses) then verbose = 0
  
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
      print,'Progress: '+str(100.0*part/npart)+'%'
      
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
