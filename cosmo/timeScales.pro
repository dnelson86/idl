; timeScales.pro
; cooling times of halo gas vs. dynamical/hubble timescales
; dnelson feb.2013

; coolingTime(): calculate primordial network cooling times for all halo gas

function coolingTime, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; set saveFilename and check for existence
  saveFilename = sP.derivPath + 'coolTime.gmem' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.sav'
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif
  
  ; load galaxy/group member catalogs for gas ids to search for
  h = loadSnapshotHeader(sP=sP)
  galcat = galaxyCat(sP=sP)
  
  ; load gas ids and match to catalog
  ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')

  match,galcat.groupmemIDs,ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
  ids_gmem_ind = ids_gmem_ind[sort(galcat_ind)]
  
  ids = !NULL
  galcat_ind = !NULL

  ; save structure
  r = { temp : fltarr(countGmem), dens : fltarr(countGmem), coolTime : fltarr(countGmem) }
  
  ; load u,dens,nelec
  u = loadSnapshotSubset(sP=sP,partType='gas',field='u')
  u = u[ids_gmem_ind]

  dens = loadSnapshotSubset(sP=sP,partType='gas',field='dens')
  dens = dens[ids_gmem_ind]

  nelec = loadSnapshotSubset(sP=sP,partType='gas',field='ne')
  nelec = nelec[ids_gmem_ind]
  
  ; calculate a temperature to save
  r.temp = convertUtoTemp(u,nelec,/log)
  r.dens = dens
  
  ; convert code units -> (physical) cgs units
  dens = codeDensToPhys(dens, scalefac=h.time, /cgs)
  u *= units.UnitPressure_in_cgs / units.UnitDensity_in_cgs ; same as UnitEnergy/UnitMass (e.g. km^2/s^2)
  
  ; calculate cooling rates using primordial network
  r.cooltime = CalcCoolTime(u,dens,nelec,scalefac=h.time) ;scalefac=h.time
  
  ; little h factor is because u is per mass
  r.cooltime *= units.HubbleParam / units.UnitTime_in_s ; convert cgs -> code units (Gyr)
  
  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
  
  return, r
  
end

; enclosedMass(): calculate the total enclosed mass (all particle types) for each gas parent in galaxy/halo catalog

function enclosedMass, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  partTypes = ['dm','gas','star','bh'] ; count mass in each of these particle types
  
  ; set saveFilename and check for existence
  saveFilename = sP.derivPath + 'massEnc.gmem.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.sav'
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,massEnc
  endif
  
  ; load galaxy/group member catalogs for gas ids to search for
  h = loadSnapshotHeader(sP=sP)
  galcat = galaxyCat(sP=sP)
  
  ; load subhalo catalog for mostBoundParticleID and for priParentIDs
  gc = loadGroupCat(sP=sP,/readIDs)

  ; load group center positions using most bound particles for each primary subgroup
  subgroupCen = subgroupPosByMostBoundID(sP=sP)
    
  ; save array
  massEnc = fltarr(n_elements(galcat.groupmemIDs))
  
  ; load radii of target gas parents
  galcatRadii = galaxyCatRadii(sP=sP)
  
  ; loop over each particle type to consider
  foreach partType,partTypes,j do begin
  
    ptNum = partTypeNum(partType)
    print,partType
    
    ; load ids and positions, restrict to fof members
    ids = loadSnapshotSubset(sP=sP,partType=partType,field='ids')
    if n_elements(ids) eq 0 then continue ; no particles of this type present in this snapshot
    
    match,ids,gc.IDs,ids_ind,gc_IDs_ind,count=countMatch,/sort
    ids_ind = ids_ind[sort(gc_IDs_ind)] ; ids_ind indexes position in snapshot in the order of gc.IDs
    gc_IDs_ind = gc_IDs_ind[sort(gc_IDs_ind)]
    
    if countMatch ne total(gc.groupLenType[ptNum,*],/int) then message,'counts'
    if countMatch eq 0 then message,'error'
    
    pos = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
    pos = pos[*,ids_ind]
    
    ; load masses
    if partType ne 'dm' then begin
      mass = loadSnapshotSubset(sP=sP,partType=partType,field='mass')
      mass = mass[ids_ind]
    endif
    
    ; now process each halo
    wMask = intarr(n_elements(gc_IDs_ind))
    
    for i=0L,gc.nGroupsTot-1 do begin
      if i mod fix(gc.nGroupsTot/10) eq 0 then print,string(100*float(i)/gc.nGroupsTot,format='(f4.1)')+'%'
      
      w = where(gc_IDs_ind ge gc.groupOffset[i] and gc_IDs_ind lt gc.groupOffset[i]+gc.groupLen[i],count)
      if count eq 0 then continue ; halo has no particles of this type

      wMask[w] += 1
      
      ; last few groups may have no subgroup
      if gc.groupFirstSub[i] lt 0 or gc.groupFirstSub[i] ge gc.nSubgroupsTot then continue
      
      ; calculate distance to parent halo centers
      ptRad = periodicDists(subgroupCen[*,gc.groupFirstSub[i]],pos[*,w],sP=sP)
    
      ; radial sort ascending
      ptRadSort = sort(ptRad)
      ptRad = ptRad[ptRadSort]
    
      ; sum cumulative masses
      if partType ne 'dm' then begin
        ptMass = (mass[w])[ptRadSort]
      endif else begin
        ptMass = replicate(h.massTable[ptNum],count)
      endelse
      
      ptMass = total(ptMass,/cum)
      
      ; handle single particle of this type in this halo (for value_locate)
      if n_elements(ptRad) eq 1 then begin
        ptRad  = [0,ptRad]
        ptMass = [0,ptMass]
      endif
      
      ; gmem
      if galcat.groupMemLen[gc.groupFirstSub[i]] gt 0 then begin
        ; match target gas parent radii and this particle type radii
        target_inds = lindgen(galcat.groupMemLen[gc.groupFirstSub[i]]) + galcat.groupmemOff[gc.groupFirstSub[i]]
        ptMatch = value_locate(ptRad,galcatRadii.gmem_sec[target_inds])
    
        ; add enclosed mass contribution
        w = where(ptMatch ge 0,count)
        if count gt 0 then massEnc[target_inds[w]] += ptMass[ptMatch[w]]
      endif
      
    endfor
    
    if min(wMask) lt 1 or max(wMask) gt 1 then message,'error wMask'

  endforeach
  
  ; save
  save,massEnc,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))

  return, massEnc
  
end

; fitRadProfile(): helper function to fit median radial profile

function fitRadProfile, radii=radii, vals=vals, range=range, radBins=radBins

  r = { binSize    : (range[1]-range[0])/radBins          ,$
        binCen     : linspace(range[0],range[1],radBins)  ,$
        radMean    : fltarr(radBins)                      ,$
        radMedian  : fltarr(radBins)                      ,$
        radStddev  : fltarr(radBins)                      ,$
        radNum     : lonarr(radBins)                       }
        
  for i=0,radBins-1 do begin
    binEdges = range[0] + [i,i+1]*r.binSize
    w = where(radii ge binEdges[0] and radii lt binEdges[1],count)
    if count gt 0 then begin
      r.radMean[i]   = mean(vals[w])
      r.radMedian[i] = median(vals[w])
      r.radStddev[i] = stddev(vals[w])
      r.radNum[i]    = count
    endif
  endfor
         
  return,r

end

; binHisto2D()

function binHisto2D, xx=xx, yy=yy, wt=wt, xmm=xmm, ymm=ymm, xbs=xbs, ybs=ybs

  compile_opt idl2, hidden, strictarr, strictarrsubs
  if ~keyword_set(xx) or ~keyword_set(yy) or ~keyword_set(xmm) or ~keyword_set(ymm) or $
     ~keyword_set(xbs) or ~keyword_set(ybs) then message,'error'
     
  if ~keyword_set(wt) then wt = replicate(1.0,n_elements(xx))
  
  nXBins = ceil((xmm[1]-xmm[0])/xbs)
  nYBins = ceil((ymm[1]-ymm[0])/ybs)
  
  binCenX = linspace(xmm[0],xmm[1]-xbs,nXBins) + xbs/2
  binCenY = linspace(ymm[0],ymm[1]-ybs,nYBins) + ybs/2
  
  ; return array
  r = { nXBins:nXBins, nYBins:nYBins, binCenX:binCenX, binCenY:binCenY, $
        binSizeX:xbs, binSizeY:ybs, $
        h2:fltarr(nXBins,nYBins)  }
        
  ; bin manually for mass weighting
  for i=0,nXBins-1 do begin
    xBin = [binCenX[i]-xbs/2,binCenX[i]+xbs/2]
  
    for j=0,nYBins-1 do begin
      yBin = [binCenY[j]-ybs/2,binCenY[j]+ybs/2]
      
      w = where(xx ge xBin[0] and xx lt xBin[1] and yy ge yBin[0] and yy lt yBin[1],count)
      if count gt 0 then r.h2[i,j] = total(wt[w])
    
    endfor
  endfor

  return, r

end

; loadFitTimescales(): load gas properties->timescales for a list of halos and also make best model fits

function loadFitTimescales, sP=sP, gcIDList=gcIDList

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  if ~keyword_set(sP) or ~n_elements(gcIDList) then message,'error'
  
  gc = loadGroupCat(sP=sP,/skipIDs)
  h  = loadSnapshotHeader(sP=sP)
  
  ; get indices for gas in halo(s)
  galcat = galaxyCat(sP=sP)
  inds = galcatINDList(sP=sP, galcat=galcat, gcIDList=gcIDList)

  ; load cooling time (Gyr) and radius (ckpc) for each gas cell
  ct = coolingTime(sP=sP)

  coolTime = ct.coolTime[inds.gmem]
  
  curTemp = ct.temp[inds.gmem]
  curDens = ct.dens[inds.gmem]
  
  gasRadii = galaxyCatRadii(sP=sP)
  
  gasVRad  = gasRadii.gmem_vrad_pri[inds.gmem]
  gasRadii = gasRadii.gmem_sec[inds.gmem]

  ct = !NULL
  
  ; for normalizing radii and vel
  gasRvir  = galcatParentProperties(sP=sP,/rVir,parNorm='pri')
  gasRvir  = gasRvir.gmem[inds.gmem]
  
  gasVcirc = galcatParentProperties(sP=sP,/vCirc,parNorm='pri')
  gasVcirc = gasVcirc.gmem[inds.gmem]
  
  ; estimate dynamical timescale using total enclosed mass at each gas cell
  encMass = enclosedMass(sP=sP)
  encMass = encMass[inds.gmem] ; code units
  
  meanDensEnc = codeDensToPhys( 3*encMass / (4 * !pi * gasRadii^3.0), scalefac=h.time ) ; code units (physical)
  
  dynTime = sqrt( 3*!pi / (32 * float(units.G) * meanDensEnc) ) ; code units (Gyr)
  
  ; age of universe
  hubbleTime = snapNumToAgeFlat(sP=sP)
  
  ; cooling and dynamical timescales for the halo as a whole (Tvir,subgroup Mtot at rvir)
  grNr = gc.subgroupGrNr[gcIDList[0]]
  
  virTemp_halo = codeMassToVirTemp(gc.subgroupMass[gcIDList[0]],sP=sP,/log)
  meanDensEnc_halo = 3*gc.subgroupMassType[partTypeNum('gas'),gcIDList[0]] / (4*!pi*gc.group_r_crit200[grNr]^3.0)
  meanDensEnc_halo = codeDensToPhys( meanDensEnc_halo, scalefac=h.time )
    
  meanDens_cgs = meanDensEnc_halo * units.UnitDensity_in_cgs * units.HubbleParam * units.HubbleParam
  coolTime_halo = CalcCoolTime(virTemp_halo,meanDens_cgs,1.0,scalefac=h.time,flag=1)
  coolTime_halo *= units.HubbleParam / units.UnitTime_in_s
  
  meanDensEnc_halo = 3*gc.subgroupMass[gcIDList[0]] / (4*!pi*gc.group_r_crit200[grNr]^3.0) ; gas+DM+stars
  meanDensEnc_halo = codeDensToPhys( meanDensEnc_halo, scalefac=h.time )
  dynTime_halo = sqrt( 3*!pi / (32 * float(units.G) * meanDensEnc_halo) )

  ; alternative definitions for dynamical time
  dynTime_halo2 = gc.group_r_crit200[grNr] / sqrt(units.G * gc.subgroupMass[gcIDList[0]] / gc.group_r_crit200[grNr] )
  dynTime_halo2 *= units.hubbleParam * h.time
  
  dynTime2 = (gasRadii * units.HubbleParam * h.time) / sqrt(units.G*encMass/gasRadii/h.time)
  
  ; radial fitting
  radFitRange = [0.15,1.5]
  
  radCt = fitRadProfile(radii=gasRadii/gasRvir,vals=coolTime,range=radFitRange,radBins=15)
  radDt = fitRadProfile(radii=gasRadii/gasRvir,vals=dynTime,range=radFitRange,radBins=40)
  radDens = fitRadProfile(radii=gasRadii/gasRvir,vals=curDens,range=radFitRange,radBins=15)
  radTemp = fitRadProfile(radii=gasRadii/gasRvir,vals=curTemp,range=radFitRange,radBins=15)
  radVrad = fitRadProfile(radii=gasRadii/gasRvir,vals=gasVRad,range=radFitRange,radBins=15)
  
  ; calculate average hot halo gas mass
  if sP.trMCPerCell eq 0 then begin
    ; SPH case: all particles have constant mass
    massPerPart = sP.targetGasMass
    masses = replicate(massPerPart,h.nPartTot[partTypeNum('gas')])
  endif else begin
    ids_gmem = galcat.groupmemIDs[inds.gmem]
  
    ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    match,ids,ids_gmem,ids_ind,ids_gmem_ind,count=countMatch
    if countMatch ne n_elements(inds.gmem) then message,'error'
    ids = !NULL
    ids_ind = ids_ind[sort(ids_gmem_ind)]
  
    masses = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
    masses = masses[ids_ind]
  endelse
  
    ; load gas ids to match to gmem ids
    ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    match,ids,galcat.groupmemIDs,ids_ind,galcat_ind,count=countGmem
    if countGmem ne n_elements(galcat.groupmemIDs) then message,'Error: Failed to find all gmem ids.'
    
    ids_ind = ids_ind[sort(galcat_ind)]
    ids = !NULL
    galcat_ind = !NULL
    
    ; load gas masses and take gmem subset
    mass = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
    mass = mass[ids_ind]
    ids_ind = !NULL
  
  mass_hot = total(masses) / n_elements(gcIDList) * units.HubbleParam

  ; profile fitting
  x = linspace(0.13,1.53,50)
  meanHaloDMMass  = mean(gc.subgroupMassType[partTypeNum('dm'),gcIDList]) * units.HubbleParam
  meanHaloGasMass = mean(gc.subgroupMassType[partTypeNum('gas'),gcIDList]) * units.HubbleParam
  
  ; SIS
  sis_dm  = sis_profile(x, mass=meanHaloDMMass, redshift=sP.redshift)
  sis_gas = sis_gas_profile(mass_hot=mass_hot, sis_dm=sis_dm, tables=tables)
  sis_fit = sis_gas_fit(rad_frac=gasRadii/gasRvir, dens=curDens, temp=curTemp, x=x)
  
  ; NFW
  nfw_dm  = nfw_profile(x, mass=meanHaloDMMass, redshift=sP.redshift)
  nfw_gas = nfw_gas_suto(mass_hot=mass_hot, nfw_dm=nfw_dm, tables=tables)
  nfw_fit = nfw_gas_fit(rad_frac=gasRadii/gasRvir, dens=curDens, temp=curTemp, x=x, nfw_dm=nfw_dm, tables=tables)

  ; convert to same units as simulation (comoving/h)
  sis_gas.rho_gas *= (h.time)^3.0/units.HubbleParam
  nfw_gas.rho_gas *= (h.time)^3.0/units.HubbleParam
  nfw_gas.rho_gas_iso *= (h.time)^3.0/units.HubbleParam

  r = {coolTime:coolTime, dynTime:dynTime, hubbleTime:hubbleTime, $
       gasRadii:gasRadii, gasVRad:gasVRad, gasRvir:gasRvir, gasVcirc:gasVcirc, $
       curTemp:curTemp, curDens:curDens, $
       virTemp_halo:virTemp_halo, coolTime_halo:coolTime_halo, dynTime_halo:dynTime_halo, $
       radCt:radCt, radDt:radDt, radDens:radDens, radTemp:radTemp, radVRad:radVRad, $
       mass_hot:mass_hot, masses:masses, x:x, $
       sis_dm:sis_dm, sis_gas:sis_gas, sis_fit:sis_fit, $
       nfw_dm:nfw_dm, nfw_gas:nfw_gas, nfw_fit:nfw_fit}
  
  return, r
end

; subsetIsotropy(): measure the isotropy of the particle distribution using a healpix powerspectrum 
;  for different subsets of the particle population

function subsetIsotropy, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; config
  sgSelect        = 'pri'
  minNumGasInHalo = 400 ; in any case this must be >>nNGB to make any sense
  subsetProp      = 'vradnorm'
  subsetRanges    = list([-5.0,5.0],[-5.0,-3.0],[-0.5,0.5],[1.0,5.0])
  
  nSide   = 64 ; for healpix map
  nNGB    = 20 ; in ThVal search
  l_split = 10 ; take isotropy as ratio of power below to power above this wavenumber
  
  nSubsets = n_elements(subsetRanges)
  nSphPx   = nSide2nPix(nSide) ; number of healpix pixels per map
  l_max    = fix(2.0 * nSide) ; maximum spherical harmonic order
  l_vals   = findgen(l_max+1)
  
  ; check if save exists
  saveFilename = sP.derivPath + 'binnedVals/binIso.' + sP.saveTag + '.' + sP.savPrefix + str(sP.res) + '.' + $
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
  print,n_elements(ww)
  
  haloMasses = codeMassToLogMsun(gc.subgroupMass[gcIDs])

  ; return array
  r = { powerSpecs   : fltarr(l_max+1,n_elements(gcIDs),nSubsets) ,$
        isoIndex     : fltarr(n_elements(gcIDs),nSubsets) + !values.f_nan ,$
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
  match,ids,galcat.groupmemIDs,ids_ind,gmem_ind,count=countGmem
  ids = !NULL
  if countGmem ne n_elements(galcat.groupmemIDs) then message,'Error: Failed to find all gmem in gas ids.'
  ids_ind = ids_ind[sort(gmem_ind)]
  
  mass = mass[ids_ind]
  pos  = pos[*,ids_ind]
  vel  = vel[*,ids_ind]
  
  ; loop over each halo
  foreach gcIndCur,gcIDs,k do begin
    print,k
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
        /silent;,/won,iter_order=2
          
      ; save power spectrum and isotropy index
      pSpec    = l_vals*(l_vals+1)*cl_gas/2/!pi
      isoIndex = total(pSpec[0:l_split]) / total(pSpec[l_split+1:*])
      
      r.powerSpecs[*,k,j] = pSpec
      r.isoIndex[k,j] = isoIndex
    endforeach
    
  endforeach

end

; timescaleFracsVsHaloMass(): gas mass fractions with tcool<tdyn, tcool<tage as a function of halo mass

function timescaleFracsVsHaloMass, sP=sP, sgSelect=sgSelect

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; config
  tsRatioVals = [0.5,1.0,1.5]
  minNumGasInHalo = 10
  nCuts = n_elements(tsRatioVals)
  
  ; check if save exists
  saveFilename = sP.derivPath + 'binnedVals/binTS.' + sP.saveTag + '.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.cut' + str(nCuts) + '.' + sgSelect + '.sav'
  
  ; results exist, return
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif 
  
  ; load catalogs
  gc = loadGroupCat(sP=sP,/skipIDs)
  gcIDList = gcIDList(gc=gc,select=sgSelect)
  h  = loadSnapshotHeader(sP=sP)
  
  galcat = galaxyCat(sP=sP)
  gcMasses = codeMassToLogMsun(gc.subgroupMass[gcIDList])
  
  gmem_hotmasses = fltarr(n_elements(gcIDList))                     
  
  ; gas mass
  if sP.trMCPerCell ne 0 then begin
    ; SPH case: all particles have constant mass
    massPerPart = sP.targetGasMass
    mass = replicate(massPerPart,h.nPartTot[partTypeNum('gas')])
  endif else begin
    ; load gas ids to match to gmem ids
    ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    match,ids,galcat.groupmemIDs,ids_ind,galcat_ind,count=countGmem
    if countGmem ne n_elements(galcat.groupmemIDs) then message,'Error: Failed to find all gmem ids.'
    
    ids_ind = ids_ind[sort(galcat_ind)]
    ids = !NULL
    galcat_ind = !NULL
    
    ; load gas masses and take gmem subset
    mass = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
    mass = mass[ids_ind]
    ids_ind = !NULL
  endelse
  
  ; load cooling time (Gyr) and radius (ckpc) for each gas cell
  coolTime = coolingTime(sP=sP)
  coolTime = coolTime.coolTime
  
  gasRadii = galaxyCatRadii(sP=sP)
  gasRadii = gasRadii.gmem_sec
  
  ; estimate dynamical timescale using total enclosed mass at each gas cell
  encMass = enclosedMass(sP=sP)
  meanDensEnc = codeDensToPhys( 3*encMass / (4 * !pi * gasRadii^3.0), scalefac=h.time ) 
  dynTime = sqrt( 3*!pi / (32 * float(units.G) * meanDensEnc) ) ; Gyr
  
  ; age of universe
  hubbleTime = snapNumToAgeFlat(sP=sP)
  
  ; bin fractions into halo mass bins and make median lines
  logMassBins   = [8.5,linspace(9.0,10.0,21),linspace(10.1,11.0,10),11.1,11.25,11.5,11.75,11.9,13.1]
  logMassNBins  = n_elements(logMassBins)-1
  logMassBinCen = 0.5 * (logMassBins + shift(logMassBins,-1))
  logMassBinCen = logMassBinCen[0:-2]
  
  ; save arrays (per halo)
  tsFracs = { gmem_tcool_tdyn    : fltarr(nCuts,n_elements(gcIDList)) + !values.f_nan  ,$
              gmem_tcool_tage    : fltarr(nCuts,n_elements(gcIDList)) + !values.f_nan   }
                  
  ; save arrays (binned in halo mass)
  tsMedian = { gmem_tcool_tdyn    : fltarr(nCuts,logMassNBins) + !values.f_nan  ,$
               gmem_tcool_tage    : fltarr(nCuts,logMassNBins) + !values.f_nan   }
  
  ; loop over each halo
  for i=0L,n_elements(gcIDList)-1 do begin
    ; indices for this halo
    gcInd = gcIDList[i]
    if galcat.groupmemLen[gcInd] eq 0 then continue
       
    inds = lindgen(galcat.groupmemLen[gcInd]) + galcat.groupmemOff[gcInd]
    
    ; get gas in this halo
    loc_coolTime = coolTime[inds]
    loc_dynTime  = dynTime[inds]
    loc_mass     = mass[inds]
    
    loc_massTot  = total(loc_mass)
    
    ; record total "hot mass" (halo mass)
    gmem_hotmasses[i] = loc_massTot
    
    ; enforce a minimum number of gas elements in gmem (leave as NaN, skip in median)
    if galcat.groupmemLen[gcInd] lt minNumGasInHalo then continue
    
    ; count
    for j=0,nCuts-1 do begin
      ; tcool/tdyn < ratio
      w = where(loc_coolTime/loc_dynTime le tsRatioVals[j],count)
      if count gt 0 then tsFracs.gmem_tcool_tdyn[j,i] = total(loc_mass[w]) / loc_massTot
      
      ; tcool/tage < ratio
      w = where(loc_coolTime/hubbleTime le tsRatioVals[j],count)
      if count gt 0 then tsFracs.gmem_tcool_tage[j,i] = total(loc_mass[w]) / loc_massTot
    endfor
    
  endfor
  
  ; bin medians vs. halo mass
  for i=0,logMassNbins-1 do begin
  
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1],count)
    
    if count gt 0 then begin
      for j=0,nCuts-1 do begin
        tsMedian.gmem_tcool_tdyn[j,i] = median(tsFracs.gmem_tcool_tdyn[j,w])
        tsMedian.gmem_tcool_tage[j,i] = median(tsFracs.gmem_tcool_tage[j,w])
      endfor
    endif
    
  endfor
  
  r = { tsFracs:tsFracs, tsMedian:tsMedian, tsRatioVals:tsRatioVals, gcMasses:gcMasses, $
        logMassBins:logMassBins, logMassBinCen:logMassBinCen, sP:sP, sgSelect:sgSelect, $
        minNumGasInHalo:minNumGasInHalo, gmem_hotmasses:gmem_hotmasses }
  
  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  return, r

end

; modelMassFracs(): calculate timescale related gas mass fractions for many SIS/NFW models
; ts: return of timescaleFracsVsHaloMass(sP=sP,sgSelect=sgSelect)
; simHotMasses: 0=analytical (f_b*halo_mass), 1=derive hot halo masses from sim

; calcMMF(): helper (timescale is either t_dyn or t_H)

function calcMMF, coolTime, timescale, tsRatio=tsRatio, sis_dm=sis_dm, $
  nfw_dm=nfw_dm, nfw_gas=nfw_gas, nfwPoly=nfwPoly, nfwIso=nfwIso

  if keyword_set(nfw_dm) and keyword_set(sis_dm) then message,'error'
  if keyword_set(nfw_dm) and ~keyword_set(nfw_gas) then message,'error'
  if keyword_set(nfw_dm) and (~keyword_set(nfwPoly) and ~keyword_set(nfwIso)) then message,'error'
  
  if keyword_set(sis_dm) then begin
    zz = coolTime ge tsRatio * timescale
    modelRadRes = n_elements(sis_dm.rad)
  endif
  if keyword_set(nfw_dm) then begin
    zz = coolTime ge tsRatio * timescale
    modelRadRes = n_elements(nfw_dm.rad)
  endif
  
  zz_tot = total(zz)
  
  ; tcool is nowhere greater than tsRatio*tdyn (all halo gas has tcool<tdyn)
  if zz_tot eq 0 then return, 1.0
    
  ; all halo gas has tcool > tsRatio*tdyn
  if zz_tot eq n_elements(modelRadRes) then return, !values.f_nan
      
  ; otherwise, find the radius where tcool = tsRatio*tdyn
  w = where(zz eq 1,count)
  if count eq 0 then message,'error'
  w = min(w)
    
  if w eq 0 then return, !values.f_nan ; left (smallest radius) boundary
    
  if keyword_set(nfwIso) then begin
    ; calculate gas mass inside e_rad and take ratio
    e_mass   = int_tabulated(nfw_dm.rad[0:w],4*!pi*nfw_dm.rad[0:w]^2.0*nfw_gas.rho_gas_iso[0:w])
    ; slight inconsistency for ISO total mass: equals hotHaloMasses[i] only for poly
    tot_mass = int_tabulated(nfw_dm.rad,4*!pi*nfw_dm.rad^2.0*nfw_gas.rho_gas_iso) 

    mmf = (e_mass/tot_mass)
  endif
      
  if keyword_set(nfwPoly) then begin
    e_mass   = int_tabulated(nfw_dm.rad[0:w],4*!pi*nfw_dm.rad[0:w]^2.0*nfw_gas.rho_gas[0:w])
    tot_mass = int_tabulated(nfw_dm.rad,4*!pi*nfw_dm.rad^2.0*nfw_gas.rho_gas) ; equals hotHaloMasses[i]
          
    mmf = (e_mass/tot_mass)
  endif
    
  if keyword_set(sis_dm) then begin
    e_rad = sis_dm.rad[w]
    mmf = (e_rad/sis_dm.r200)
  endif

  return, mmf
end

function modelMassFracs, sP=sP, ts=ts, simHotMasses=simHotMasses
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  tables = interpLambdaSD93()  

  ; config
  modelMassRange = [8.8,12.2]
  modelMassRes   = 250
  
  modelRadRange = [0.02,1.0]
  modelRadRes   = 1000
  
  modelMasses = [10.0^linspace(modelMassRange[0],modelMassRange[1],modelMassRes)/1e10] ; code units

  ; check if save exists
  simHotMassesTag = '.mAna' ; analytical (f_b * M_halo)
  if keyword_set(simHotMasses) then simHotMassesTag = '.mSim' ; derive hot halo mass from sim
  
  saveFilename = sP.derivPath + 'binnedVals/binMMF.' + sP.saveTag + '.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.res' + str(modelMassRes) + '.rad' + str(modelRadRes) + simHotMassesTag + '.sav'
  
  ; results exist, return
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif
  
  ; hubble time
  t_H = redshiftToAgeFlat(sP.redshift)
  
  ; derive hot halo masses from simulation?
  if keyword_set(simHotMasses) then begin
    polyDeg = 2
    radBins = 6
    
    ; coarse median binning
    x_fit = ts.gcMasses
    y_fit = ts.gmem_hotMasses / logMsunToCodeMass(ts.gcMasses)
    fit = fitRadProfile(radii=x_fit, vals=y_fit, range=[10.0,modelMassRange[1]], radBins=radBins)
    
    w = where(fit.radMedian eq 0,count)
    if count gt 0 then message,'Error: Cut out the zeros first.'
    
    ; spline interp the median binned points
    fracs = interpol(fit.radMedian,fit.binCen,codeMassToLogMsun(modelMasses),/spline)
    hotHaloMasses = fracs * modelMasses
    
    if min(hotHaloMasses) lt 0.0 then message,'Error: Gone negative interpolating the fraction.'

  endif else begin
    hotHaloMasses = modelMasses * units.f_b
    print,'Warning: Analytical is overestimate since we neglect what has already entered the galaxy.'
  endelse
  
  r = { SIS1_tDyn     : fltarr(n_elements(ts.tsRatioVals),modelMassRes) ,$
        SIS2_tDyn     : fltarr(n_elements(ts.tsRatioVals),modelMassRes) ,$
        NFW_iso_tDyn  : fltarr(n_elements(ts.tsRatioVals),modelMassRes) ,$
        NFW_poly_tDyn : fltarr(n_elements(ts.tsRatioVals),modelMassRes) ,$
        SIS1_tH       : fltarr(n_elements(ts.tsRatioVals),modelMassRes) ,$
        SIS2_tH       : fltarr(n_elements(ts.tsRatioVals),modelMassRes) ,$
        NFW_iso_tH    : fltarr(n_elements(ts.tsRatioVals),modelMassRes) ,$
        NFW_poly_tH   : fltarr(n_elements(ts.tsRatioVals),modelMassRes) ,$
        modelMasses            : modelMasses    ,$
        modelMassRange         : modelMassRange ,$
        modelRadRes            : modelRadRes    ,$
        modelRadRange          : modelRadRange  ,$
        hotHaloMasses          : hotHaloMasses   }
  
  x = linspace(modelRadRange[0],modelRadRange[1],modelRadRes)
  
  for i=0,modelMassRes-1 do begin
    ; MODEL 1: SIS uniform halo
    sis_dm_rVir  = sis_profile(1.0, mass=modelMasses[i], redshift=sP.redshift)
    sis_gas_rVir = sis_gas_profile(mass_hot=hotHaloMasses[i], sis_dm=sis_dm_rVir, tables=tables)
    
    ; MODEL 2: SIS with tcool(r)
    sis_dm  = sis_profile(x, mass=modelMasses[i], redshift=sP.redshift)
    sis_gas = sis_gas_profile(mass_hot=hotHaloMasses[i], sis_dm=sis_dm, tables=tables)
    
    ; MODEL 3: NFW
    nfw_dm  = nfw_profile(x, mass=modelMasses[i], redshift=sP.redshift)
    nfw_gas = nfw_gas_suto(mass_hot=hotHaloMasses[i], nfw_dm=nfw_dm, tables=tables)    
    
    foreach tsRatio,ts.tsRatioVals,j do begin
    
      ; MODEL 1: SIS uniform halo (all halo gas either has tcool<ts or tcool>ts)
      if sis_gas_rVir.coolTime lt tsRatio * sis_gas_rVir.dynTime then r.SIS1_tDyn[j,i] = 1.0
      if sis_gas_rVir.coolTime ge tsRatio * sis_gas_rVir.dynTime then r.SIS1_tDyn[j,i] = 0.0
      
      if sis_gas_rVir.coolTime lt tsRatio * t_H then r.SIS1_tH[j,i] = 1.0
      if sis_gas_rVir.coolTime ge tsRatio * t_H then r.SIS1_tH[j,i] = 0.0
    
      ; MODEL 2: SIS with tcool(r)
      r.SIS2_tDyn[j,i] = calcMMF(sis_gas.coolTime, sis_gas.dynTime, tsRatio=tsRatio, sis_dm=sis_dm)
      r.SIS2_tH[j,i]   = calcMMF(sis_gas.coolTime, t_H, tsRatio=tsRatio, sis_dm=sis_dm)
    
      ; MODEL 3: NFW ISO
      r.NFW_iso_tDyn[j,i] = calcMMF(nfw_gas.coolTime_iso, nfw_gas.dynTime, $
                                    tsRatio=tsRatio, nfw_dm=nfw_dm, nfw_gas=nfw_gas, /nfwIso)
      r.NFW_iso_tH[j,i]   = calcMMF(nfw_gas.coolTime_iso, t_H, $
                                    tsRatio=tsRatio, nfw_dm=nfw_dm, nfw_gas=nfw_gas, /nfwIso)

      ; MODEL 3: NFW POLY
      r.NFW_poly_tDyn[j,i] = calcMMF(nfw_gas.coolTime, nfw_gas.dynTime, $
                                     tsRatio=tsRatio, nfw_dm=nfw_dm, nfw_gas=nfw_gas, /nfwPoly)
      r.NFW_poly_tH[j,i]   = calcMMF(nfw_gas.coolTime, t_H, $
                                     tsRatio=tsRatio, nfw_dm=nfw_dm, nfw_gas=nfw_gas, /nfwPoly)
      
    endforeach
    
  endfor ; modelMassRes
  
  ; make extremely low values not plot
  w = where(r.NFW_poly_tDyn lt 0.0001,count)
  if count gt 0 then r.NFW_poly_tDyn[w] = !values.f_nan
  w = where(r.NFW_iso_tDyn lt 0.0001,count)
  if count gt 0 then r.NFW_iso_tDyn[w] = !values.f_nan
  
  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  return,r
  
end

; verifyGadgetCoolingTimes(): check calculated cooling times vs. output in snapshots (where available)

pro verifyGadgetCoolingTimes

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  sP = simParams(res=128,run='gadget',redshift=2.0)
  
  ; load galaxy/group member catalogs for gas ids to search for
  h = loadSnapshotHeader(sP=sP)
  galcat = galaxyCat(sP=sP)
  gc = loadGroupCat(sP=sP,/skipIDs)
  
  ; calculate cooling times
  ct = coolingTime(sP=sP)
  
  ; load cooling rates from snapshot and calculate cooling times
  cr  = loadSnapshotSubset(sP=sP,partType='gas',field='coolingRate')
  u   = loadSnapshotSubset(sP=sP,partType='gas',field='u')
  
  ; older (Gadget code): cr = u / tcool (all code units = Gyr)
  ct_snap = u / cr
  
  ; load snapshot ids for gmem subset
  ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')

  match,galcat.groupmemIDs,ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
  ids_gmem_ind = ids_gmem_ind[sort(galcat_ind)]
  
  ct_snap_gmem = ct_snap[ids_gmem_ind]
    
  start_PS,'test.eps'
    cgPlot,ct.cooltime,ct_snap_gmem,psym=3,xtitle="calc cooltime",ytitle="snap cooltime",$
      xrange=[0,50],yrange=[0,50],/xs,/ys
  end_PS
  
  start_PS,'test2.eps'
    w = where(ct.cooltime ne 0 and finite(ct_snap_gmem))
    xx = ct.cooltime[w]-ct_snap_gmem[w]
    w = where(xx le 10.0)
    hist = histogram(xx[w],bin=0.1,loc=loc)
    cgPlot,loc+0.05,hist,xtitle="delta cooltime",ytitle="N",/ylog,yrange=[1,1e5],yminor=0
  end_PS
  
  start_PS,'test3.eps'
    w = where(ct.cooltime ne 0 and finite(ct_snap_gmem))
    xx = alog10(abs(ct.cooltime[w]-ct_snap_gmem[w]))
    w = where(finite(xx))
    hist = histogram(xx[w],bin=0.1,loc=loc)
    cgPlot,loc+0.05,hist,xtitle="log abs delta cooltime",ytitle="N"
  end_PS
  
  stop

end

; crotonTest

pro crotonTest
 
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
 
  ; config
  zBounds = [0.0,7.0]
  zSteps  = 100
  
  M0 = 100.0 ; 10^12 Msun
  
  redshifts = linspace(zBounds[1],zBounds[0],zSteps)
  times = redshiftToAgeFlat(redshifts) * 1e9 ; yr
  
  ; load Sutherland & Dopita (1993) cooling tables
  tables = interpLambdaSD93()
  
  ; arrays: halo properties
  M_DM    = fltarr(zSteps) ; Msun, total halo mass
  M_b     = fltarr(zSteps) ; Msun, total baryonic halo mass
  R_vir   = fltarr(zSteps) ; km
  V_c     = fltarr(zSteps) ; km/s
  T_vir   = fltarr(zSteps) ; K
  dM_b    = fltarr(zSteps) ; Msun/yr
  
  ; Croton model
  M_cold     = fltarr(zSteps) ; Msun, e.g. the only destination for hot gas with no SF/ejecta
  M_hot      = fltarr(zSteps) ; Msun, e.g. Mgas0 (hot gas available to cool)
  dM_cool    = fltarr(zSteps) ; Msun/yr
  r_cool     = fltarr(zSteps) ; km
  rapid_mode = bytarr(zSteps) ; 0=no, 1=yes
  
  ; Kang model
  M_cold_kang     = fltarr(zSteps) ; Msun, e.g. the only destination for hot gas with no SF/ejecta
  M_hot_kang      = fltarr(zSteps) ; Msun, e.g. Mgas0 (hot gas available to cool)
  dM_cool_kang    = fltarr(zSteps) ; Msun/yr
  r_cool_kang     = fltarr(zSteps) ; km
  rapid_mode_kang = bytarr(zSteps) ; 0=no, 1=yes
  
  ; timestep zero
  M_b[0] = units.f_b * haloMAH(redshifts[0],M0=M0,z0=zBounds[0])
  
  ; evolve in time
  ; --------------
  for i=1,zSteps-1 do begin
    ; halo properties
    delta_t = ( times[i] - times[i-1] )
    
    ; new halo baryonic mass
    M_DM[i]  = haloMAH(redshifts[i],M0=M0,z0=zBounds[0])
    M_b[i]   = M_DM[i] * units.f_b
    
    ; rate of growth of halo baryons
    dM_b[i] = (M_b[i] - M_b[i-1]) / delta_t ; 1e10 Msun/yr
    
    ; SIS DM halo
    sis_dm   = sis_profile(1.0, mass=M_DM[i], redshift=redshifts[i])
    R_vir[i] = sis_dm.r200
    
    ; Croton model
    ; ------------
      M_hot[i] = M_DM[i] * units.f_b - M_cold[i-1]
      
      if M_hot[i] gt 0.0 then begin
      
        sis_gas = sis_gas_profile(mass_hot=M_hot[i], sis_dm=sis_dm, tables=tables)
        r_cool[i] = sis_gas.r_cool_croton
    
        ; compare rcool to rvir, "rapid mode" if rcool>rvir
        if sis_gas.r_cool_croton ge sis_dm.r200 then rapid_mode[i] = 1B
      
        ; set cold gas accretion rate based on mode (1e10 Msun/yr)
        if rapid_mode[i] then begin
          ; rapid mode
          dM_cool[i] = M_hot[i] / delta_t ;(sis_gas.dynTime*1e9)
        endif else begin
          ; slow mode
          dM_cool[i] = 0.5 * M_hot[i] * sis_gas.r_cool_croton / (sis_dm.r200 * sis_gas.dynTime_halo*1e9)
        endelse
      
        if dM_cool[i] * delta_t gt M_hot[i] then begin
          print,'dM_cool * dt > M_hot',i
          dM_cool[i] = M_hot[i] / delta_t
        endif
      
        ; update cold gas
        M_cold[i] = M_cold[i-1] + dM_cool[i] * delta_t
      
      endif else begin
        print,'warning M_hot zero',i
      endelse
      
    ; Kang model
    ; ----------
      M_hot_kang[i] = M_DM[i] * units.f_b - M_cold_kang[i-1]
      
      if M_hot_kang[i] le 0.0 then message,'error'
      
      sis_gas = sis_gas_profile(mass_hot=M_hot_kang[i], sis_dm=sis_dm, tables=tables)
      r_cool_kang[i] = sis_gas.r_cool_hubble
      
      ; compare rcool_h to rvir, "rapid mode" if rcool_h>rvir
      if sis_gas.r_cool_hubble ge sis_dm.r200 then rapid_mode_kang[i] = 1B
      
      ; set cold gas accretion rate based on mode
      kang_timenorm = times[i] ;sis_gas.hubbleTime*1e9 ;
      
      if rapid_mode_kang[i] then begin
        dM_cool_kang[i] = M_hot_kang[i] / kang_timenorm
      endif else begin
        dM_cool_kang[i] = 0.5 * M_hot_kang[i] * sis_gas.r_cool_hubble / (sis_dm.r200 * kang_timenorm)
      endelse
      
      if dM_cool_kang[i] * delta_t gt M_hot_kang[i] then message,'error'
      
      ; update cold gas
      M_cold_kang[i] = M_cold_kang[i-1] + dM_cool_kang[i] * delta_t
      
  endfor
  
  ; plot
  xrange = [1,7]
  
  start_PS,'halo_SAM_vs_redshift.eps', ys=8, xs=6
    ; halo baryons
    cgPlot,1+redshifts,dM_b*1e10,xtitle='',ytitle=textoidl('dM/dt [M_{sun}/year]'),$
      /xlog,xrange=xrange,/xs,xminor=0,/ylog,yrange=[0.06,60.0],/ys,yminor=0,$
      xtickv=[1,2,3,4,5,6,7],xtickname=replicate(' ',10),xticks=6,position=(sP.pos_3x1)[0]
      
    ; croton model, dM/dt cold split into rapid and slow
    cgPlot,1+redshifts,dM_cool*1e10,/overplot,color=cgColor('red')
    
    w_rapid = where(rapid_mode eq 1B,count)
    cgPlot,1+redshifts[w_rapid],dM_cool[w_rapid]*1e10,psym=4,/overplot,color=cgColor('red')
    
    ; kang model
    cgPlot,1+redshifts,dM_cool_kang*1e10,/overplot,color=cgColor('blue')
    
    w_rapid = where(rapid_mode_kang eq 1B,count)
    cgPlot,1+redshifts[w_rapid],dM_cool_kang[w_rapid]*1e10,psym=4,/overplot,color=cgColor('blue')

    ; halo baryons
    cgPlot,1+redshifts,M_b*1e10,xtitle='',ytitle=textoidl('M_{cold} [M_{sun}]'),$
      /xlog,xrange=[1,7],/xs,xminor=0,/ylog,yrange=[1.5e8,2e11],/ys,yminor=0,$
      xtickv=[1,2,3,4,5,6,7],xtickname=replicate(' ',10),xticks=6,/noerase,position=(sP.pos_3x1)[1]
      
    ; croton model, dM/dt cold split into rapid and slow
    cgPlot,1+redshifts,M_cold*1e10,/overplot,color=cgColor('red')
    
    w_rapid = where(rapid_mode eq 1B,count)
    cgPlot,1+redshifts[w_rapid],M_cold[w_rapid]*1e10,psym=4,/overplot,color=cgColor('red')
    
    ; kang model
    cgPlot,1+redshifts,M_cold_kang*1e10,/overplot,color=cgColor('blue')
    
    w_rapid = where(rapid_mode_kang eq 1B,count)
    cgPlot,1+redshifts[w_rapid],M_cold_kang[w_rapid]*1e10,psym=4,/overplot,color=cgColor('blue')
    
    legend,['croton','kang'],textcolors=['red','blue'],box=0,/top,/right

    cgPlot,[0],[0],/nodata,xtitle='1+z',ytitle='Radius [kpc]',$
      /xlog,xrange=[1,7],/xs,xminor=0,yrange=[0,150],/ys,$
      xtickv=[1,2,3,4,5,6,7],xtickname=['1','2','3','4','5','6','7'],xticks=6,/noerase,position=(sP.pos_3x1)[2]
      
    cgPlot,1+redshifts,R_vir,/overplot,color=cgColor('black')
    cgPlot,1+redshifts,R_cool,/overplot,color=cgColor('red')
    cgPlot,1+redshifts,r_cool_kang,/overplot,color=cgColor('blue')
    legend,textoidl(['R_{vir}','R_{cool,croton}','R_{cool,kang}']),textcolors=['black','red','blue'],box=0,$
      position=[2.5,140]
  end_PS
  
end
  
; compVsMass(): at one redshift, compute models as a function of halo mass
; calculate critical halo mass crossover for rcool<>rvir
  
pro compVsMass
  
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()  
  
  ; config
  redshift = 2.0
  M0_range = [10.0,13.0] ; log(msun) at z=0
  
  hRes = 200
  hMasses = 10.0^(linspace(M0_range[0],M0_range[1],hRes))/1e10
  
  rPts = linspace(0.01,10.0,5000)
  
  ; load Sutherland & Dopita (1993) cooling tables
  tables = interpLambdaSD93()
  
  ; arrays
  hMass_z = fltarr(hRes)
  virRad  = fltarr(hRes)
  virTemp = fltarr(hRes)
  
  sis_coolRad = fltarr(hRes)
  nfw_coolRad = fltarr(hRes)
  nfw_coolRad_iso = fltarr(hRes)
  
  for i=0,hRes-1 do begin
    ; calculate halo masses / hot gas mass, at target redshift
    hMass_z[i] = haloMAH(redshift,M0=hMasses[i],z0=0.0)
    m_gas0 = hMass_z[i] * units.f_b
    print,codeMassToLogMsun(hMass_z[i])
    
    ; SIS DM halo and gas profiles
    sis_dm  = sis_profile(rPts, mass=hMass_z[i], redshift=redshift)
    sis_gas = sis_gas_profile(mass_hot=m_gas0, sis_dm=sis_dm, tables=tables)

    nfw_dm  = nfw_profile(rPts, mass=hMass_z[i], redshift=redshift)
    nfw_gas = nfw_gas_suto(mass_hot=m_gas0, nfw_dm=nfw_dm, tables=tables)
    
    ; radii, code units (kpc)
    virRad[i]  = sis_dm.r200
    virTemp[i] = sis_dm.Tvir_rVir
    
    sis_coolRad[i] = sis_gas.r_cool
    nfw_coolRad[i] = nfw_gas.r_cool
    nfw_coolRad_iso[i] = nfw_gas.r_cool_iso
  endfor
  
  ; plot
  start_PS,'halo_rads_vs_mass.eps'
    cgPlot,[0],[0],/nodata,xtitle='Halo Mass at z=2',ytitle='Radius [kpc]',$
      xrange=[8.5,12],/xs,yrange=[0,220],/ys
      
    ; mark tvir=10^4 (end of cooling tables)
    w = min(where(alog10(virTemp) ge 4.0))
    cgPlot,codeMassToLogMsun([hMass_z[w],hMass_z[w]]),[50,200],/overplot,line=1
      
    ; rvir and rcool
    cgPlot,codeMassToLogMsun(hMass_z),virRad,/overplot,color=cgColor('black'),line=2
    
    cgPlot,codeMassToLogMsun(hMass_z),sis_coolRad,/overplot,color=cgColor('orange'),line=0
    cgPlot,codeMassToLogMsun(hMass_z),nfw_coolRad,/overplot,color=cgColor('red'),line=0
    cgPlot,codeMassToLogMsun(hMass_z),nfw_coolRad_iso,/overplot,color=cgColor('blue'),line=0
    
    legend,textoidl(['R_{vir}','SIS R_{cool}','NFW poly R_{cool}','NFW iso R_{cool}']),$
      textcolors=['black','orange','red','blue'],linestyle=[2,0,0,0],box=0,/top,/right,linesize=0.25
  end_PS
  
  stop
  
end
  
; compProfiles(): for one redshift and halo mass, compare different theoretical profiles
;                 density,temp,timescales vs radius
  
pro compProfiles
  
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()    
  
  ; config
  redshift = 2.0
  hMass    = 10.0 ; code units at redshift
  
  ; load Sutherland & Dopita (1993) cooling tables
  tables = interpLambdaSD93()
  
  ; arrays
  rRes = 2000
  rPts = linspace(0.01,5.0,rRes) ; r/rvir
  
  ; DM and gas profiles
  sis_dm  = sis_profile(rPts, mass=hMass, redshift=redshift)
  sis_gas = sis_gas_profile(mass_hot=hMass*units.f_b, sis_dm=sis_dm, tables=tables)
  
  nfw_dm  = nfw_profile(rPts, mass=hMass, redshift=redshift)
  nfw_gas = nfw_gas_suto(mass_hot=hMass*units.f_b, nfw_dm=nfw_dm, tables=tables)

  ; plot
  xrange = [0.01,1.5]
  pos = list([0.18,0.67,0.95,0.95] ,$
             [0.18,0.39,0.95,0.67] ,$
             [0.18,0.11,0.95,0.39])
  
  start_PS,'halo_denstempts_vs_rad.eps', ys=8, xs=6
    cgPlot,[0],[0],/nodata,xtitle='',ytitle=textoidl('\rho(r) / \rho_{crit,z=2}'),$
      xrange=xrange,/xs,/xlog,xminor=0,yrange=10.0^[1.1,8],/ys,/ylog,yminor=0,position = pos[0],xtickname=replicate(' ',10)
      
    ; DM profile
    cgPlot,rPts,rhoRatioToCrit(sis_dm.rho_DM,redshift=redshift),/overplot,color=cgColor('gray')
    cgPlot,rPts,rhoRatioToCrit(nfw_dm.rho_DM,redshift=redshift),/overplot,color=cgColor('black')
    
    ; r_s line for NFW
    cgPlot,[nfw_dm.r_s/nfw_dm.r200,nfw_dm.r_s/nfw_dm.r200],[1.5,3],line=2,/overplot
    
    ; gas profile
    cgPlot,rPts,rhoRatioToCrit(sis_gas.rho_gas,redshift=redshift),/overplot,line=2,color=cgColor('orange')
    cgPlot,rPts,rhoRatioToCrit(nfw_gas.rho_gas,redshift=redshift),/overplot,line=2,color=cgColor('red')
    cgPlot,rPts,rhoRatioToCrit(nfw_gas.rho_gas_iso,redshift=redshift),/overplot,line=2,color=cgColor('blue')
    
    legend,['SIS DM','NFW DM'],textcolors=['gray','black'],box=0,/right,/top
    
    ; temp
    cgPlot,[0],[0],/nodata,xtitle="",xtickname=replicate(' ',10),ytitle=textoidl('T_{gas} [log K]'),$
       xrange=xrange,/xs,/xlog,xminor=0,yrange=[5,7],/ys,/noerase,position = pos[1]
      
    cgPlot,rPts,alog10(sis_gas.temp_gas),/overplot,color=cgColor('orange')
    cgPlot,rPts,alog10(nfw_gas.temp_gas),/overplot,color=cgColor('red')
    cgPlot,rPts,alog10(replicate(nfw_gas.T_0,n_elements(rPts))),/overplot,color=cgColor('blue')
    
    legend,['SIS','NFW poly','NFW iso'],textcolors=['orange','red','blue'],box=0,/right,/top
    
    ; timescales
    cgPlot,[0],[0],/nodata,xtitle=textoidl(' r / r_{vir} '),ytitle='Timescale [Gyr]',$
      xrange=xrange,/xs,yrange=[0.001,5.0],/ys,/ylog,/xlog,xminor=0,yminor=0,/noerase,position = pos[2]
      
    cgPlot,rPts,sis_gas.coolTime,/overplot,color=cgColor('orange'),line=0
    cgPlot,rPts,sis_gas.dynTime,/overplot,color=cgColor('orange'),line=2
    cgPlot,rPts,nfw_gas.coolTime,/overplot,color=cgColor('red'),line=0
    cgPlot,rPts,nfw_gas.dynTime,/overplot,color=cgColor('red'),line=2
    cgPlot,rPts,nfw_gas.coolTime_iso,/overplot,color=cgColor('blue'),line=0
    cgPlot,rPts,nfw_gas.dynTime*0.97,/overplot,color=cgColor('blue'),line=2 ; visual offset
    
    ; verify r_cool
    y_rcool = [0.003,0.03]
    cgPlot,[sis_gas.r_cool,sis_gas.r_cool]/nfw_dm.r200,y_rcool,line=1,color=cgColor('orange'),/overplot
    cgPlot,[nfw_gas.r_cool,nfw_gas.r_cool]/nfw_dm.r200,y_rcool,line=1,color=cgColor('red'),/overplot
    cgPlot,[nfw_gas.r_cool_iso,nfw_gas.r_cool_iso]/nfw_dm.r200,y_rcool,line=1,color=cgColor('blue'),/overplot
    
    legend,textoidl(['\tau_{cool}','\tau_{dyn}']),linestyle=[0,2],box=0,/top,/left,linesize=0.25
    
  end_PS
  
  stop
end

