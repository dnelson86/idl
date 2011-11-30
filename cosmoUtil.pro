; cosmoUtil.pro
; cosmological simulations - utility functions
; dnelson nov.2011

; simParams(): return structure of simulation parameters

function simParams, res=res, run=run

  if not keyword_set(res) or not keyword_set(run) then begin
     print,'Error: simParams: arguments not specified.'
     exit
  endif

  r = {simPath:      '',$    ; root path to simulation snapshots and group catalogs
       savPrefix:    '',$    ; save prefix for simulation (make unique, e.g. 'G')
       boxSize:      0.0,$   ; boxsize of simulation, kpc
       snapRange:    [0,0],$ ; snapshot range of simulation
       groupCatRange:[0,0],$ ; snapshot range of fof/subfind catalogs (subset of above)
       plotPath:     '',$    ; path to put plots
       galCatPath:   '',$    ; path to galaxy catalog files
       thistPath:    '',$    ; path to thermal history files
       derivPath:    '',$    ; path to put derivative (data) files
       snap:         0,$     ; convenience for passing between functions
       res:          0,$     ; copied from input
       run:          '',$    ; copied from input
       minNumGasPart: 0}     ; minimum number of gas particles in a galaxy to include in e.g. mass func
       
  r.boxSize       = 20000.0
  r.snapRange     = [0,314]
  r.groupCatRange = [50,314]
  r.minNumGasPart = 100
  
  if (isnumeric(res)) then $
    r.res = res
  r.run = run
  
  r.plotPath   = '/n/home07/dnelson/coldflows/'
  r.galCatPath = '/n/home07/dnelson/coldflows/galaxycat/'
  r.thistPath  = '/n/home07/dnelson/coldflows/thermhist/'
  r.derivPath  = '/n/home07/dnelson/coldflows/data.files/'
       
  if (run eq 'gadget') then begin
    r.simPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'
    r.savPrefix = 'G'
    return,r
  endif

  print,'simParams: ERROR.'
  exit
end

; sgIDList(): return a sub-list of subgroup IDs from a given group catalog sg
; 
; pri=1 : members of the first subgroup of each group only ("background"/"main subhalo"/"halo")
; sec=1 : members of the non-first subgroups of each group only ("satellites"/"subhalos")
; all=1 : all subgroups
; 
; minNumPart=1 : require a minimum number of particles in (subfind) group to include

function sgIDList, sP=sP, sg=sg, pri=pri, sec=sec, all=all, minNumPart=minNumPart

    if (not keyword_set(pri) and not keyword_set(sec) and not keyword_set(all)) then begin
      print,'Error: Must specify sgID parameter.'
      return,0
    endif
    
    ; load galaxy cat if necessary
    if not keyword_set(sg) then begin
      if not keyword_set(sP) then begin
        print,'Error: Must specify sg or sP.'
        return,0
      endif
      sg = loadSubhaloGroups(sP.simPath,sP.snap)
    endif

    minNumPartVal = 100 ; in subgroup (total, dm+gas+stars)

    prevGrNr   = -1
    valSGids   = []

    if keyword_set(pri) then begin
    
      ; "background"/"main subhalos"/"halos" only
      for i=0,n_elements(sg.subgroupLen)-1 do begin
        if (sg.subgroupGrnr[i] eq prevGrNr) then begin
          prevGrNr = sg.subgroupGrnr[i]
        endif else begin
        
          if (keyword_set(minNumPart)) then begin
            if (sg.subgroupLen[i] ge minNumPartVal) then $
              valSGids = [valSGids,i]
          endif else begin
            valSGids = [valSGids,i]
          endelse ;minNumPart
          
          prevGrNr = sg.subgroupGrnr[i]
        endelse
      endfor
      
    endif
    
    if keyword_set(sec) then begin
    
      ; "satellites"/"subhalos" only
      for i=0,n_elements(sg.subgroupLen)-1 do begin
        if (sg.subgroupGrnr[i] ne prevGrNr) then begin
          prevGrNr = sg.subgroupGrnr[i]
        endif else begin
        
          if (keyword_set(minNumPart)) then begin
            if (sg.subgroupLen[i] ge minNumPartVal) then $
              valSGids = [valSGids,i]
          endif else begin
            valSGids = [valSGids,i]
          endelse ;minNumPart
          
        endelse
      endfor
      
    endif
    
    if keyword_set(all) then begin
    
      ; both primary and secondary
      for i=0,n_elements(sg.subgroupLen)-1 do begin
      
        if (keyword_set(minNumPart)) then begin
          if (sg.subgroupLen[i] ge minNumPartVal) then $
            valSGids = [valSGids,i]
        endif else begin
          valSGids = [valSGids,i]
        endelse ;minNumPart
        
      endfor
    
    endif
    
    return, valSGids

end

; sgPIDList(): return a list of member particle IDs from a given group catalog sg
;
; pri=1 : members of the first subgroup of each group only ("background"/"main subhalo"/"halo")
; sec=1 : members of the non-first subgroups of each group only ("satellites"/"subhalos")
; all=1 : the union of the above (one of these three options must be set)
;
; gasOnly=1 : load gas ids from snapshot and match to restrict return to only gas particle ids
; dmOnly=1  : load dm ids from snapshot and match to restrict return to only dark matter ids

function sgPIDList, sg=sg, pri=pri, sec=sec, all=all, gasOnly=gasOnly, dmOnly=dmOnly

  if (not keyword_set(sg) or (not keyword_set(pri) and not keyword_set(sec) and $
                              not keyword_set(all))) then begin
    print,'Error: sgPIDList: Bad inputs.'
    return,0
  endif

  ; get list of appropriate group ids
  valSGids = sgIDList(sg=sg,pri=pri,sec=sec,all=all)
  
  ; make list of particle ids in these groups
  subgroupPIDs = []
  
  foreach sgID, valSGids do begin
    ; select subgroup
    sgIDs = sg.subGroupIDs[sg.subGroupOffset[sgID] : sg.subGroupOffset[sgID] + $
                           sg.subGroupLen[sgID] - 1]
    
    subgroupPIDs = [subgroupPIDs, sgIDs]
  endforeach
    
  ; particle type restriction
  if (keyword_set(gasOnly) or keyword_set(dmOnly)) then begin
    if (keyword_set(gasOnly)) then partName='gas'
    if (keyword_set(dmOnly))  then partName='dm'
    
    ; load gas/dm ids from targetSnap and restrict
    ;pIDs = loadSnapshotSubset(gadgetPath,snapNum=snap,partType=partName,field='ids')
    
    ;match, pIDs, sgEnd_subgroupIDs, pIDs_ind, sgIDs_ind, count=count_pID
    
    ;pIDs = !NULL
    ;pIDs_ind = !NULL
    
    ; keep only sgIDs[sgIDs_ind]
    ;if (count_pID gt 0) then $
    ;    sgEnd_subgroupIDs = sgEnd_subgroupIDs[sgIDs_ind]
    print,'ERROR need to fix this part'
    exit
  endif

  return, subgroupPIDs
end

; galCatRepParentIDs(): for the galaxy catalog, replicate the list of ordered parent IDs such that
;                       the return array is the same size as the number of gas particle ids with
;                       each element the id of its parent
;
; sgIDListPri=1 : should be a list of primary subhalo IDs, e.g. sgIDListPri=sgIDList(sg=sg,/pri)
;                 if set, for secondary groups, return instead of the ID of the primary parent
;                 (for primary groups the return is unchanged)

function galCatRepParentIDs, gc=gc, sgIDListPri=sgIDListPri

    ; arrays to hold parent IDs for each gas particle
    sg_ind_gal  = lonarr(n_elements(gc.galaxyIDs))
    sg_ind_gmem = lonarr(n_elements(gc.groupmemIDs))
    
    if not keyword_set(sgIDListPri) then $
      sgIDListPri = indgen(n_elements(gc.galaxyLen)) ; valid id list set to all
    
    for sgID=0,n_elements(gc.galaxyLen)-1 do begin
    
        ; parent ID (pri or sec)
        if (total(sgID eq sgIDListPri) eq 0) then begin
          dists = sgID - sgIDListPri
          dists = dists[where(dists gt 0)]
          w = where(dists eq min(dists),count)
          parID = sgIDListPri[w[0]]
        endif else begin
          parID = sgID
        endelse
        
        ; galaxies
        groupStart = gc.galaxyOff[sgID]
        groupEnd   = groupStart + gc.galaxyLen[sgID]
        
        if (gc.galaxyLen[sgID] gt 0) then $
          sg_ind_gal[groupStart:groupEnd-1] = cmreplicate(parID,groupEnd-groupStart)
          
        ; group members
        groupStart = gc.groupmemOff[sgID]
        groupEnd   = groupStart + gc.groupmemLen[sgID]
        
        if (gc.groupmemLen[sgID] gt 0) then $
          sg_ind_gmem[groupStart:groupEnd-1] = cmreplicate(parID,groupEnd-groupStart)
    endfor
    
    r = {gal:sg_ind_gal,gmem:sg_ind_gmem}
    return,r
end

; correctPeriodicDistVecs(): enforce periodic B.C. for distance vecotrs (effectively component by 
;                            component), input vecs in format fltarr[3,n]

pro correctPeriodicDistVecs, vecs

  boxSize = 20000.0 ;kpc, valid for ComparisonProject
  
  w = where(vecs gt boxSize*0.5,count)
  if (count ne 0) then $
    vecs[w] = vecs[w] - boxSize
    
  w = where(vecs lt -boxSize*0.5,count)
  if (count ne 0) then $
    vecs[w] = boxSize + vecs[w]

end

; redshiftToSnapNum(): convert redshift to the nearest snapshot number

function redshiftToSnapNum, redshiftList, verbose=verbose

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
  
  ; for -1 return !NULL
  if (total(redshiftList) eq -1) then return,!NULL
  
  ; for redshift zero hard select #314
  snapNum = intarr(n_elements(redshiftList))
  
  foreach redshift,redshiftList,i do begin
    if (redshift eq 0.0) then begin
      snapNum[i] = 314
    endif else begin
      w = where(abs(redshifts - redshift) eq min(abs(redshifts - redshift)),count)
      if (count eq 2) then w = w[0]
      if (count eq 0) then return,-1
      snapNum[i] = w[0]
    endelse
  
    if (verbose) then $
      print,'Found nearest snapshot to z = ' + str(redshift) + ' (num = ' + str(snapNum) + ')'
  endforeach

  if (n_elements(snapNum) eq 1) then snapNum = snapNum[0]

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

; codeMassToVirTemp(): convert halo mass (in code units) to virial temperature at specified redshift

function codeMassToVirTemp, mass, redshift, meanmolwt=meanmolwt

  units = getUnits()

  ; mean molecular weight default (valid for ComparisonProject)
  if not keyword_set(meanmolwt) then meanmolwt = 0.6

  ; mass to msun
  mass_msun = mass * (units.UnitMass_in_g / units.Msun_in_g)
  
  ; cosmo
  omega_m   = 0.27
  omega_L   = 0.73
  omega_k   = 0.0
  little_h  = 0.7
  
  omega_m_z = omega_m * (1+redshift)^3.0 / $
              ( omega_m*(1+redshift)^3.0 + omega_L + omega_k*(1+redshift)^2.0 )
  
  Delta_c = 18*!pi^2 + 82*(omega_m_z-1.0) - 39*(omega_m_z-1.0)^2.0

  Tvir = 1.98e4 * (meanmolwt/0.6) * (mass_msun/1e8*little_h)^(2.0/3.0) * $
         (omega_m/omega_m_z * Delta_c / 18.0 / !pi^2.0)^(1.0/3.0) * $
         (1.0 + redshift)/10.0 ;K

  return, Tvir
  
end

; codeMassToLogMsun(): convert mass in code units to log(msun)

function codeMassToLogMsun, mass

  units = getUnits()
  
  mass_msun = mass * (units.UnitMass_in_g / units.Msun_in_g)
  
  ; log of nonzero
  w = where(mass_msun eq 0.0,count)
  if (count ne 0) then $
    mass_msun[w] = 1.0
  
  return,alog10(mass_msun)
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

function snapNumToAge, snap
  z = snapNumToRedshift(snapNum=snap)
  return, redshiftToAge(z)
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
    ; h2rt /= float(max(h2rt)) ;norm s.t. colorbar gives fraction wrt max cell
  endif else begin
    ; hist2d (weighted)
    h2rt = hist2d(dens,temp,mass,$
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

; redshift_axis(): draw redshift axis

pro redshift_axis, xRange, yRange, ylog=ylog, snapnum=snapnum, dotted=dotted, zTicknames=zTicknames

  logFac = 1.0
  if keyword_set(ylog) then logFac = 10.0
  
  yTickLen = (yRange[1]-yRange[0]) / 40.0 * logFac
  yTextOff = (yRange[1]-yRange[0]) / 50.0 * logFac

  if (not keyword_set(zTicknames)) then $
    zTicknames = ['30','6','4','3','2','1','0.5','0.25','0']
  nZ = n_elements(zTicknames)
  
  zXPos = fltarr(nZ)
  
  ; plot "maximum" redshift label
  if (xRange[0] le 0.0) then $
    fsc_text,0.0,yRange[1]+yTextOff,zTicknames[0],alignment=0.5
  
  ; skip z=30 (highest) at t=0
  for i=1,nZ-1 do begin
    if keyword_set(snapnum) then begin ;x-axis in snapshot number
      zXPos[i] = redshiftToSnapNum(float(zTicknames[i]))
    endif else begin ;x-axis in time [gyr]
      zXPos[i] = redshiftToAge(float(zTicknames[i]))
    endelse
    
    ; plot tick mark and label if inside plotrange
    if (zXPos[i] ge xRange[0] and zXPos[i] le xRange[1]) then begin
      fsc_plot,[zXPos[i],zXPos[i]],[yRange[1],yRange[1]-yTickLen],/overplot
      fsc_text,zXPos[i],yRange[1]+yTextOff,zTicknames[i],alignment=0.5
    endif
    
    ; plot vertical dotted line at each redshift mark if requested
    if keyword_set(dotted) then $
      fsc_plot,[zXPos[i],zXPos[i]],yRange,line=dotted,/overplot,thick=!p.thick-0.5
  endfor
  
  fsc_plot,xRange,[yRange[1],yRange[1]],/overplot
  fsc_text,0.5,0.94,"Redshift",/normal

end

; universeage_axis(): draw age of universe (elapsed) axis

pro universeage_axis, xRange, yRange, ylog=ylog, dotted=dotted

  logFac = 1.0
  if keyword_set(ylog) then logFac = 10.0
  
  yTickLen = (yRange[1]-yRange[0]) / 40.0 * logFac
  yTextOff = (yRange[1]-yRange[0]) / 50.0 * logFac

  if (not keyword_set(zTicknames)) then begin
    zTicknames = ['0.1','1','1.5','2','3','4','5','6','8','13.8'] ;10=0.33
    zRedshifts = [30,5.7,4.21,3.25,2.23,1.65,1.265,0.98,0.595,0.0]
  endif
  nZ = n_elements(zTicknames)
  
  zXPos = fltarr(nZ)
  
  ; plot "maximum" redshift label
  if (xRange[0] le 0.0) then $
    fsc_text,0.0,yRange[1]+yTextOff,zTicknames[0],alignment=0.5
  
  ; skip t=0 (highest) at z=inf (?)
  for i=1,nZ-1 do begin

    zXPos[i] = zRedshifts[i]
    
    ; plot tick mark and label if inside plotrange
    if (zXPos[i] le xRange[0] and zXPos[i] ge xRange[1]) then begin
      fsc_plot,[zXPos[i],zXPos[i]],[yRange[1],yRange[1]-yTickLen],/overplot
      fsc_text,zXPos[i],yRange[1]+yTextOff,zTicknames[i],alignment=0.5
    endif
    
    ; plot vertical dotted line at each redshift mark if requested
    if keyword_set(dotted) then $
      fsc_plot,[zXPos[i],zXPos[i]],yRange,line=dotted,/overplot,thick=!p.thick-0.5
  endfor
  
  fsc_plot,xRange,[yRange[1],yRange[1]],/overplot
  fsc_text,(xRange[0]-xRange[1])/2.0,yRange[1]+yTextOff*5.0,"Time [Gyr]",alignment=0.5

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
