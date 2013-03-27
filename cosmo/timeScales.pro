; timeScales.pro
; cooling times of halo gas vs. dynamical/hubble timescales
; dnelson mar.2013

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

  calcMatch,galcat.groupmemIDs,ids,galcat_ind,ids_gmem_ind,count=countGmem
  ids_gmem_ind = ids_gmem_ind[calcSort(galcat_ind)]
  
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
    
    calcMatch,ids,gc.IDs,ids_ind,gc_IDs_ind,count=countMatch
    ids_ind = ids_ind[calcSort(gc_IDs_ind)] ; ids_ind indexes position in snapshot in the order of gc.IDs
    gc_IDs_ind = gc_IDs_ind[calcSort(gc_IDs_ind)]
    
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

; loadFitTimescales(): load gas properties->timescales for a list of halos and also make best model fits
; accTimesRepTR=1 : also return gmem accretion times (for each tracer) and so also need to replicate 
;                   all other values by the child tracer counts of each gmem gas parent at sP.snap

function loadFitTimescales, sP=sP, gcIDList=gcIDList, accTimesRepTR=accTimesRepTR

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
  
  dynTime = sqrt( 3*!pi / (32 * float(units.G * meanDensEnc)) ) ; code units (Gyr)
  
  ; age of universe
  hubbleTime = snapNumToAgeFlat(sP=sP)
  
  ; cooling and dynamical timescales for the halo as a whole (Tvir,subgroup Mtot at rvir)
  grNr = gc.subgroupGrNr[gcIDList[0]]
  
  virTemp_halo = codeMassToVirTemp(gc.subgroupMass[gcIDList[0]],sP=sP,/log)
  meanDensEnc_halo = 3*gc.subgroupMassType[partTypeNum('gas'),gcIDList[0]] / (4*!pi*gc.group_r_crit200[grNr]^3.0)
  meanDensEnc_halo = codeDensToPhys( meanDensEnc_halo, scalefac=h.time )
    
  meanDens_cgs = meanDensEnc_halo * units.UnitDensity_in_cgs * units.HubbleParam * units.HubbleParam
  coolTime_halo = CalcCoolTime(virTemp_halo,meanDens_cgs,1.0,scalefac=h.time,flag=1)
  coolTime_halo = coolTime_halo[0] * float(units.HubbleParam / units.UnitTime_in_s)
  
  meanDensEnc_halo = 3*gc.subgroupMass[gcIDList[0]] / (4*!pi*gc.group_r_crit200[grNr]^3.0) ; gas+DM+stars
  meanDensEnc_halo = codeDensToPhys( meanDensEnc_halo, scalefac=h.time )
  dynTime_halo = float(sqrt( 3*!pi / (32 * float(units.G) * meanDensEnc_halo) ))

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
    calcMatch,ids,ids_gmem,ids_ind,ids_gmem_ind,count=countMatch
    if countMatch ne n_elements(inds.gmem) then message,'error'
    ids = !NULL
    ids_ind = ids_ind[calcSort(ids_gmem_ind)]
  
    masses = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
    masses = masses[ids_ind]
  endelse
  
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
  
  ; finally, should we load the gmemAccretionTimes and therefor replicate all other values by 
  ; the child tracer counts? (the fits and the models are unaffected)
  if keyword_set(accTimesRepTR) then begin
    ; load accretionRates
    ar = accretionRates(sP=sP)
    
    ; replicate parent indices (duplicates) of all tracers of requested halos only
    inds_tr = replicate_var(ar.galcat_gmem_cc,subset_inds=inds.gmem)

    ; convert accSnap -> accSnapOffset -> Gyr elapsed for hot gas to galaxy accretion
    redshifts = snapNumToRedshift(sP=sP,/all)
    ages      = redshiftToAgeFlat(redshifts)
    
    w = where(ar.accSnap_gmem[inds_tr.child_inds] ne -1,count)
    accTime = ages[ar.accSnap_gmem[inds_tr.child_inds[w]]] - (ages[sP.snap])[0]
    ar = !NULL
    
    ; replicate all properties which are currently one per gas cell, only for those with known accTime

    ; cannot here take e.g. coolTime[parInds] since parInds are global (all gmem), but coolTime is
    ; already the inds.gmem subset for a single halo, a mass range, etc, so do a value_locate approach
    parInds = value_locate(inds.gmem,inds_tr.parent_inds[w])
    
    if ~array_equal(inds.gmem[parInds],inds_tr.parent_inds[w]) then message,'Error: Not 1-to-1 in locate.'

    coolTime = coolTime[parInds]
    dynTime  = dynTime [parInds]
    gasRadii = gasRadii[parInds]
    gasVRad  = gasVRad [parInds]
    gasRvir  = gasRvir [parInds]
    gasVcirc = gasVcirc[parInds]
    curTemp  = curTemp [parInds]
    curDens  = curDens [parInds]
    
    ; override masses (weights) with constant mass
    if sP.trMCPerCell le 0 then masses = replicate(sP.targetGasMass,n_elements(accTime)) ; SPH or vel tracer
    if sP.trMCPerCell gt 0 then masses = replicate(sP.trMassConst,n_elements(accTime)) ; MC tracer  
    
    inds_tr = !NULL
    parInds = !NULL
  endif else begin
    accTime = -1
  endelse

  r = {coolTime:coolTime, dynTime:dynTime, hubbleTime:hubbleTime, accTime:accTime, $
       gasRadii:gasRadii, gasVRad:gasVRad, gasRvir:gasRvir, gasVcirc:gasVcirc, $
       curTemp:curTemp, curDens:curDens, nGas:n_elements(inds.gmem), $
       virTemp_halo:virTemp_halo, coolTime_halo:coolTime_halo, dynTime_halo:dynTime_halo, $
       radCt:radCt, radDt:radDt, radDens:radDens, radTemp:radTemp, radVRad:radVRad, $
       mass_hot:mass_hot, masses:masses, x:x, $
       sis_dm:sis_dm, sis_gas:sis_gas, sis_fit:sis_fit, $
       nfw_dm:nfw_dm, nfw_gas:nfw_gas, nfw_fit:nfw_fit}
  
  return, r
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
    calcMatch,ids,galcat.groupmemIDs,ids_ind,galcat_ind,count=countGmem
    if countGmem ne n_elements(galcat.groupmemIDs) then message,'Error: Failed to find all gmem ids.'
    
    ids_ind = ids_ind[calcSort(galcat_ind)]
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
