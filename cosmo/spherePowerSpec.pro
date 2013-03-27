; spherePowerSpec.pro
; healpix power spectrum, isotropy parameter
; dnelson mar.2013

; read_alm(): parse FITS output of anafast_cxx healpix for a_lm coefficients

function read_alm, l_vals, fitsFileName
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  alm_in = mrdfits(fitsFileName,1,/silent)
  spawn, 'rm '+fitsFileName, ret
        
  ;alm.index ; l*l+l+m+1 (m <= l) (which means l+1 m for each l)
  totCoeffs = total(l_vals+1,/int)
      
  if totCoeffs ne n_elements(alm_in) then message,'Error: Bad size of fits input.'
      
  ; return structure
  alm = {l : lonarr(totCoeffs) ,$
         m : lonarr(totCoeffs) ,$
         a : fltarr(totCoeffs)  }
      
  ; separate out l,m indices for each alm
  l = 0L & m = 0L
      
  for i=0L,n_elements(alm_in)-1 do begin
    alm.l[i] = l
    alm.m[i] = m
        
    m += 1
    if m gt l then begin
      l += 1
      m = 0
    endif
        
  endfor
      
  if ~array_equal(alm_in.index, alm.l*alm.l+alm.l+alm.m+1) then message,'Error: Bad l,m decomposition.'
    
  ; overwrite alm structure with alm amplitudes
  alm.a = sqrt(alm_in.real * alm_in.real + alm_in.imag * alm_in.imag)
      
  return,alm
end

; isotropyParam(): calculate "isotropy parameter" in various ways

function isotropyParam, l_vals, cl=cl, l_split=l_split, $
                                alm=alm, l_window=l_window
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  r = {cl1      : !values.f_nan ,$
       cl2      : !values.f_nan ,$
       cl_split : !values.f_nan ,$
       alm1     : !values.f_nan ,$
       alm2     : !values.f_nan  }

  ; power spectrum
  if keyword_set(cl) then begin
    pSpec = l_vals*(l_vals+1)*cl/2/!pi
    if pSpec[0] ne 0.0 then message,'strange, this is mean subtracted'
    
    r.cl1 = total(pSpec[2:*]) / pSpec[1]
    r.cl2 = total(pSpec)
    
    if keyword_set(l_split) then $
      r.cl_split = total(pSpec[0:l_split]) / total(pSpec[l_split+1:*])
  endif
  
  ; alm coefficients
  if keyword_set(alm) then begin
    if ~keyword_set(l_window) then message,'Error: Must specify l window for alm IP.'
    
    window_vals = lindgen(l_window[1]-l_window[0]) + l_window[0]
    
    r.alm1 = 0.0
    r.alm2 = 0.0
    
    foreach l_val,window_vals do begin
      w = where(alm.l eq l_val,count)
      if count eq 0 then message,'Error'
    
      ; (1)
      a = alm.a[w] ; / mean(alm.a[w])
      a -= mean(a)
    
      r.alm1 += total(abs(a))
    
      ; (2)
      a = alm.a[w]
      r.alm2 += stddev(a)
    endforeach
    
  endif
  
  return,r
end

; haloShellAngPSpecTheory(): make synthetic healpix maps and test the power spectra output

pro haloShellAngPSpecTheory

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  nSide    = 128
  angSizes = [10,30,50,70] ; deg, radius (180=all sky)
  l_modes  = [2,7,10,20]
  
  ; choose maximum spherical harmonic order l_max
  ;l_max = fix(3.0 * nSide - 1)
  l_max = fix(2.0 * nSide)

  ; start plot (1)
  sP = simParams(res=128,run='tracer',redshift=2.0)
  start_PS, sP.plotPath + 'powerspec.theory.eps'
      
    l_vals = findgen(l_max+1)
    unit_lambda = 2*!pi / (l_vals + 0.5) ; wavelength on unit sphere
    ang_size = unit_lambda * 180.0 / !pi
      
    first_l = 1
    xrange = [first_l,max(l_vals)]
    yrange = [1e-8,1.0]
      
    !y.margin[1] += 1.0
      
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,/xlog,/ylog,yminor=0,xminor=0,$
      ytitle=textoidl("l(l+1)C_l/2\pi"),xtitle="l (Spherical Wavenumber)",title=""
        
    ; draw wavelength axis (linear map)
    cgAxis,xs=1,xaxis=1,xminor=0,xtitle="Angular Size [deg]",xrange=[ang_size[first_l],min(ang_size)]
      
    ; loop over angular sizes
    foreach angSize,angSizes,k do begin
      ; empty healpix map
      healpix_data = fltarr(nside2npix(nSide))
        
      ; add a feature of this angular size to the map
      angSizeRad = angSize * !dtor
      query_disc, nSide, [0,0,1], angSizeRad, pxInds, /nested
      
      pxFrac = float(n_elements(pxInds)) / n_elements(healpix_data)
      healpix_data[pxInds] = 1.0
      healpix_data -= mean(healpix_data)
      healpix_data /= (max(healpix_data)-min(healpix_data))
        
      ; calculate power spectrum and overplot
      ianafast,healpix_data,cl,nlmax=l_max,/cxx,/nested,tmpdir='/n/home07/dnelson/',$
        /silent;,/won,iter_order=2
      
      cgPlot,l_vals,l_vals*(l_vals+1)*cl/2/!pi,line=0,color=cgColor(units.colors[k]),/overplot

      ; isotropy parameter
      ip = isotropyParam(l_vals,cl=cl)
        
      print,angSize,pxFrac,ip.cl1,ip.cl2
    endforeach
    
    ; legends
    legend,string(angSizes),textcolors=units.colors[0:k],box=0,/bottom,/left
    
  end_PS
  
  ; start plot (2)
  rot_ang = [0,0]
  
  start_PS, sP.plotPath + 'powerspec.theory.vis_' + str(min(l_modes)) + '_' + str(max(l_modes)) + '.eps'
  
    pos = ['ul_nb','ur_nb','ll_nb','lr_nb']
    xtpos = [0.06,0.56,0.06,0.56]
    ytpos = [0.55,0.55,0.13,0.13]
        
    for i=0,3 do begin
      ; make c_l spectrum
      cl = fltarr(n_elements(l_vals))
        
      ; fill single mode
      cl[l_modes[i]] = 1.0
        
      ; generate synthetic healpix map
      isynfast,cl,healpix_data,simul_type=1,nside=nSide,tmpdir='/n/home07/dnelson/',/silent
      print,'l = ' + str(l_modes[i])
        
      ; convert healpix map to nested system
      healpix_data = reorder(healpix_data, /r2n)
  
      if i gt 0 then bartitle = "" else bartitle = "Relative Intensity"
      plotMollweideProj,healpix_data,rot_ang=rot_ang,title="",bartitle=bartitle,$
        pos=pos[i],noerase=(i gt 0),/bigbar
  
      cgText,xtpos[i]-0.01,ytpos[i]-0.01,'l = '+str(l_modes[i]),/normal,alignment=0.5
    endfor
    
  end_PS, pngResize=60, /deletePS
  
  ; start plot (3)
  start_PS, sP.plotPath + 'powerspec.theory.vis2_' + str(min(l_modes)) + '_' + str(max(l_modes)) + '.eps'
  
    yrange = [1e-13,1e3]
    !y.margin[1] += 1.0
      
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,/xlog,/ylog,yminor=0,xminor=0,$
      ytitle=textoidl("l(l+1)C_l/2\pi"),xtitle="l (Spherical Wavenumber)",title=""
        
    ; draw wavelength axis (linear map)
    cgAxis,xs=1,xaxis=1,xminor=0,xtitle="Angular Size [deg]",xrange=[ang_size[first_l],min(ang_size)]
        
    for i=0,3 do begin
      ; make c_l spectrum
      cl = fltarr(n_elements(l_vals))
        
      ; fill single mode
      cl[l_modes[i]:l_modes[i]+2] = 1.0
        
      ; generate synthetic healpix map
      isynfast,cl,healpix_data,simul_type=1,nside=nSide,tmpdir='/n/home07/dnelson/',/silent
        
      ; convert healpix map to nested system
      healpix_data = reorder(healpix_data, /r2n)
      healpix_data -= mean(healpix_data)
      healpix_data /= (max(healpix_data)-min(healpix_data))
      
      ; run power spectra
      ianafast,healpix_data,cl,nlmax=l_max,/cxx,/nested,tmpdir='/n/home07/dnelson/',$
        /silent,alm1_out='/n/home07/dnelson/alm.temp.fits';,/won,iter_order=2
        
      cgPlot,l_vals,l_vals*(l_vals+1)*cl/2/!pi,line=0,color=cgColor(units.colors[i]),/overplot

      ; recover a_lm coefficients and calculate isotropy parameters
      alm = read_alm(l_vals,'/n/home07/dnelson/alm.temp.fits')
      
      l_window = [20,23]
      ip = isotropyParam(l_vals,cl=cl,alm=alm,l_window=l_window)
        
      print,i,l_modes[i],ip.cl1,ip.cl2,ip.alm1,ip.alm2
    endfor
    
    ; legends
    legend,'l = ' + str(l_modes),textcolors=units.colors[0:i],box=0,/top,/right
  end_PS
  
  stop

end

; haloShellAngPowerSpec(): compute and plot the angular power spectrum
;                          e.g. the variation of the spherical harmonic coefficients

pro haloShellAngPowerSpec
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  radInds      = [4,5,7]  ; pre-saved radFacs
  hMassTargets = [12.0,11.5,11.0,10.5]
  cutSubS      = 1 ; cut satellite substructures out from halo
  fieldName    = 'density'
  l_split      = 10 ; take isotropy as ratio of power below to power above this wavenumber    
  
  sP = simParams(res=512,run='tracer',redshift=2.0)  
 
  ; locate primary subgroupID for mass target
  subgroupIDs = massTargetToHaloID(hMassTargets,sP=sP)
 
  foreach subgroupID,subgroupIDs,m do begin
    ; interpolate onto shells at a set of fixed radii
    hsd_gas = haloShellValue(sP=sP,partType='gas',valName=fieldName,subgroupID=subgroupID,cutSubS=cutSubS)
    hsd_dm  = haloShellValue(sP=sP,partType='dm',valName=fieldName,subgroupID=subgroupID,cutSubS=cutSubS)
    
    ; choose maximum spherical harmonic order l_max
    ;l_max = fix(3.0 * hsd_gas.nSide - 1)
    l_max = fix(2.0 * hsd_gas.nSide)
  
    l_vals = findgen(l_max+1)
    unit_lambda = 2*!pi / (l_vals + 0.5) ; wavelength on unit sphere
    ang_size = unit_lambda * 180.0 / !pi
      
    first_l = 1
    xrange = [first_l,max(l_vals)]
    yrange =[1e-7,1e-1]
  
    ; plot
    start_PS, sP.plotPath + 'powerspec.shell.' + sP.savPrefix + str(sP.res) + '.' + $
      str(sP.snap) + '.h' + str(subgroupID) + '.eps', /big
      !y.margin[1] += 1.0
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,/ylog,yminor=0,/xlog,xminor=0,$
        ytitle=textoidl("l(l+1)C_l/2\pi"),xtitle="l (Spherical Wavenumber)",title=""
        
      ; draw wavelength axis (linear map)
      cgAxis,xs=1,xaxis=1,xminor=0,xtitle="Angular Size [deg]",xrange=[ang_size[first_l],min(ang_size)]
      
      ; for this radius, calculate angular power spectrum of gas and dm
      foreach radInd,radInds,k do begin
        ;lambda = unit_lambda * hsd_gas.radFacs[radInd] * hsd_gas.rVir ; wavelength in ckpc
        
        healpix_data = hsd_gas.value[*,radInd] - mean(hsd_gas.value[*,radInd])
        healpix_data = reform(healpix_data)
        healpix_data /= (max(healpix_data)-min(healpix_data))
        ianafast,healpix_data,cl_gas,nlmax=l_max,/cxx,/nested,tmpdir='/n/home07/dnelson/',$
          /silent,alm1_out='/n/home07/dnelson/alm.temp.fits';,/won,iter_order=2
      
        alm_gas = read_alm(l_vals,'/n/home07/dnelson/alm.temp.fits')
      
        healpix_data = hsd_dm.value[*,radInd] - mean(hsd_dm.value[*,radInd])
        healpix_data = reform(healpix_data)
        healpix_data /= (max(healpix_data)-min(healpix_data))
        ianafast,healpix_data,cl_dm,nlmax=l_max,/cxx,/nested,tmpdir='/n/home07/dnelson/',$
          /silent,alm1_out='/n/home07/dnelson/alm.temp.fits';,/won,iter_order=2        
        
        alm_dm = read_alm(l_vals,'/n/home07/dnelson/alm.temp.fits')
        
        ; plot power spectra
        cgPlot,l_vals,l_vals*(l_vals+1)*cl_gas/2/!pi,line=k,color=cgColor(units.colors[0]),/overplot
        cgPlot,l_vals,l_vals*(l_vals+1)*cl_dm/2/!pi,line=k,color=cgColor(units.colors[1]),/overplot

        ; isotropy parameters
        l_window = [10,50]
        
        ip_gas = isotropyParam(l_vals,cl=cl_gas,l_split=l_split,alm=alm_gas,l_window=l_window)
        ip_dm  = isotropyParam(l_vals,cl=cl_dm,l_split=l_split,alm=alm_dm,l_window=l_window)
                
        print,hMassTargets[m],hsd_gas.radFacs[radInd],$
              ip_gas.cl1,ip_dm.cl1,ip_gas.cl_split,ip_dm.cl_split,$
              ip_gas.alm1,ip_gas.alm2,ip_dm.alm1,ip_dm.alm2
      endforeach
      
      ; legends
      legend,[textoidl("r/r_{vir} = ")+string(hsd_gas.radFacs[radInds],format='(f4.2)')],$
        linestyle=indgen(n_elements(radInds)),linesize=0.25,box=0,/top,/right
      legend,['gas','dm'],textcolor=units.colors[0:1],box=0,/bottom,/left
      legend,[sP.run+" "+str(sP.res)+textoidl("^3")+" log(M) = "+$
        string(hMassTargets[m],format='(f4.1)')],box=0,/top,/left,textcolor=['light gray']
    end_PS
    
  endforeach ;subgroupIDs
  stop
end

; subsetIsotropy(): measure the isotropy of the particle distribution using a healpix powerspectrum 
;  for different subsets of the particle population

function subsetIsotropy, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; config
  sgSelect        = 'pri'
  minNumGasInHalo = 3000 ; in any case this must be >>nNGB to make any sense
  subsetProp      = 'vradnorm'
  subsetRanges    = list([-5.0,5.0],[-5.0,-3.0],[-0.5,0.5],[1.0,5.0])
  
  nSide    = 64 ; for healpix map
  nNGB     = 20 ; in ThVal search
  l_split  = 10 ; CL: take isotropy as ratio of power below to power above this wavenumber
  l_window = [20,50] ; ALM: take isotropy as a_lm coefficient width within l window
  nIPs     = 4 ; number of isotropy parameters to store per subset
  
  nSubsets = n_elements(subsetRanges)
  nSphPx   = nSide2nPix(nSide) ; number of healpix pixels per map
  l_max    = fix(2.0 * nSide) ; maximum spherical harmonic order
  l_vals   = findgen(l_max+1)
  
  ; check if save exists
  saveFilename = sP.derivPath + 'binnedVals/binIsoNew.' + sP.saveTag + '.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.' + subsetProp + '_' + str(n_elements(subsetRanges)) + '.' + sgSelect + '.sav'
  
  ; results exist, return
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif 
  
  ; halo list and centers
  gc     = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen  = subgroupPosByMostBoundID(sP=sP)
  
  gcIDs  = gcIDList(gc=gc,select=sgSelect)
  
  ; restrict to reasonable sized halos
  ww = where(gc.subgroupLenType[partTypeNum('gas'),gcIDs] ge minNumGasInHalo,count)
  gcIDs = gcIDs[ww]
  haloMasses = codeMassToLogMsun(gc.subgroupMass[gcIDs])

  print,'Processing ['+str(n_elements(ww))+'] halos, with masses ['+string(min(haloMasses),format='(f4.1)')+' to '+$
    string(max(haloMasses),format='(f4.1)')+']'
  
  ; return array
  r = { powerSpecs   : fltarr(l_max+1,n_elements(gcIDs),nSubsets) ,$
        isoIndex     : fltarr(nIPs,n_elements(gcIDs),nSubsets) + !values.f_nan ,$
        sgSelect     : sgSelect      ,$
        subsetProp   : subsetProp    ,$
        subsetRanges : subsetRanges  ,$
        nSide        : nSide         ,$
        nNGB         : nNGB          ,$
        l_vals       : l_vals        ,$
        l_split      : l_split       ,$
        minNumGasInHalo : minNumGasInHalo ,$
        haloMasses      : haloMasses      ,$
        gcIDs           : gcIDs            }
  
  ; load gas ids, positions, velocities and masses
  ids  = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
  pos  = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  vel  = loadSnapshotSubset(sP=sP,partType='gas',field='vel')
  mass = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
  
  ; now restrict all these quantities to gmem only
  galcat = galaxyCat(sP=sP)
  calcMatch,ids,galcat.groupmemIDs,ids_ind,gmem_ind,count=countGmem
  ids = !NULL
  if countGmem ne n_elements(galcat.groupmemIDs) then message,'Error: Failed to find all gmem in gas ids.'
  ids_ind = ids_ind[calcSort(gmem_ind)]
  
  mass = mass[ids_ind]
  pos  = pos[*,ids_ind]
  vel  = vel[*,ids_ind]
  
  ; loop over each halo
  foreach gcIndCur,gcIDs,k do begin
    if (k mod 50) eq 0 then print,'['+str(k)+'] ' + string(float(k)/n_elements(gcIDs),format='(f4.1)')+'%'
    
    ; halo properties
    haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcIndCur]] ;ckpc
    haloV200   = sqrt(units.G * gc.subgroupMass[gcIndCur] / haloVirRad )
    
    ; get indices into gmem array for this halo
    loc_inds = galcatINDList(sP=sP, galcat=galcat, gcIDList=[gcIndCur])
    loc_inds = loc_inds.gmem
    
    ; make positions relative
    loc_pos  = fltarr(3,n_elements(loc_inds))
    
    for i=0,2 do begin
      cDist = pos[i,loc_inds] - sgcen[i,gcIndCur]
      correctPeriodicDistVecs, cDist, sP=sP
      loc_pos[i,*] = cDist
    endfor
    
    cDist = !NULL
    rad = sqrt(reform(loc_pos[0,*]^2.0 + loc_pos[1,*]^2.0 + loc_pos[2,*]^2.0))
    
    ; take subsets of other gmem arrays
    loc_mass = mass[loc_inds]
    
    if subsetProp eq 'vradnorm' then begin
      loc_vel  = vel[*,loc_inds]
      
      loc_subProp = reform( (loc_vel[0,*]-gc.subgroupVel[0,gcIndCur]) * loc_pos[0,*] + $
                            (loc_vel[1,*]-gc.subgroupVel[1,gcIndCur]) * loc_pos[1,*] + $
                            (loc_vel[2,*]-gc.subgroupVel[2,gcIndCur]) * loc_pos[2,*])
      loc_subProp /= rad ; vrad [km/s]
      loc_subProp /= haloV200 ; vradnorm
    endif
    
    ; move all points to surface of virial radius sphere
    for i=0,2 do loc_pos[i,*] *= (haloVirRad/rad)

    rad = !NULL
    
    ; loop over subsets
    foreach subsetRange,subsetRanges,j do begin
      ; take subset
      wSubsetCut = where(loc_subProp ge subsetRange[0] and loc_subProp le subsetRange[1],count)
      
      if count le nNGB then continue ; CalcTHVal requires nNGB points or more
      
      ; ThVal config (from valName eq 'density' in haloShellValue)
      loc_subset_value = reform(loc_mass[wSubsetCut],[1,n_elements(wSubsetCut)]) ; 1xN
      loc_subset_pos   = loc_pos[*,wSubsetCut]
      
      posval = [loc_subset_pos,loc_subset_value]
      thMode = 3 ; total/volume
      
      ; calculate healpix
      sphereXYZ = sphereXYZCoords(Nside=nSide,radius=haloVirRad,center=[0,0,0])
      
      if sP.trMCPerCell eq 0 then begin
        healpix_data = CalcTHVal(posval,sphereXYZ,ndims=3,nNGB=nNGB,thMode=thMode,boxSize=sP.boxSize)
      endif else begin
        ; if this is an Arepo run, make the mass subset for weighting and do tophat estimate
        weights = reform(loc_mass[wSubsetCut],[1,n_elements(wSubsetCut)])
        posvalwt = [posval,weights]
        healpix_data = CalcTHVal(posvalwt,sphereXYZ,ndims=3,nNGB=nNGB,thMode=thMode,boxSize=sP.boxSize,/weighted)
      endelse
      
      ; calculate powerspectrum
      healpix_data -= mean(healpix_data)
      healpix_data /= (max(healpix_data)-min(healpix_data))
      ianafast,healpix_data,cl_gas,nlmax=l_max,/cxx,/nested,tmpdir='/n/home07/dnelson/',$
        /silent,alm1_out='/n/home07/dnelson/alm.temp.fits';,/won,iter_order=2
          
      alm_gas = read_alm(l_vals,'/n/home07/dnelson/alm.temp.fits')    
          
      ; save power spectrum and isotropy index
      pSpec = l_vals*(l_vals+1)*cl_gas/2/!pi
      ip    = isotropyParam(l_vals,cl=cl_gas,l_split=l_split,alm=alm_gas,l_window=l_window)
      
      r.powerSpecs[*,k,j] = pSpec
      r.isoIndex[0,k,j]   = ip.cl_split
      r.isoIndex[1,k,j]   = ip.cl1
      r.isoIndex[2,k,j]   = ip.alm1
      r.isoIndex[3,k,j]   = ip.alm2
    endforeach
    
  endforeach

  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))

  return, r

end

; plotIsotropyVsHaloMass()

pro plotIsotropyVsHaloMass

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; config
  sP = simParams(res=512,run='tracer',redshift=2.0)

  ; load
  si = subsetIsotropy(sP=sP)
  gc = loadGroupCat(sP=sP,/skipIDs)
  h  = loadSnapshotHeader(sP=sP)
  
  ; plot
  colors = ['red','blue','forest green','orange']
  sK     = 1 ; smoothing kernel
  xrange = [10.0,12.0]
  yrange = [0.01,3.0]
  frange = [10.26,11.97]

  start_PS, sP.plotPath+'isotropy_vs_halomass.'+str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/ylog,/xs,/ys,yminor=0,$
      xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
      ytitle=textoidl("\zeta_{gas} (Isotropy Index)")
         
    ; median lines    
    foreach subsetRange,si.subsetRanges,j do begin
      ; individual elements
      cgPlot,si.haloMasses,si.isoIndex[*,j],psym=4+j,color=cgColor('gray'),thick=1.0,/overplot
      
      ; radial fit and plot median
      radFit = fitRadProfile(radii=si.haloMasses,vals=si.isoIndex[*,j],range=frange,radBins=8)
      cgPlot,radFit.binCen,smooth(radFit.radMedian,sK),psym=-16,$
        color=cgColor(colors[j]),/overplot ; filled circle
      
    endforeach
    
    ; minimum halo mass measured
    cgPlot,[frange[0],frange[0]],[yrange[0]*2,yrange[1]*0.7],$
      line=1,color=cgColor('dark gray'),/overplot
    
    ; legend, redo plot borders
    strings = ['all','-5 < v_{rad} / v_{circ} < -3','-1/2 < v_{rad} / v_{circ} < 1/2',$
               '1 < v_{rad} / v_{circ} < 5']
    strings = ['all','inflow','zero','outflow']
    legend,textoidl(strings),textcolor=colors,psym=4+indgen(n_elements(si.subsetRanges)),box=0,/top,/right
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/ylog,yminor=0,/xs,/ys,/noerase
    
  end_PS
  stop

end
