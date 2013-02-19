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
  gasRadii = gasRadii.gmem_sec[inds.gmem]
  
  gasRvir  = galcatParentProperties(sP=sP,/rVir,parNorm='pri')
  gasRvir  = gasRvir.gmem[inds.gmem]
  
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
       gasRadii:gasRadii, gasRvir:gasRvir, curTemp:curTemp, curDens:curDens, $
       virTemp_halo:virTemp_halo, coolTime_halo:coolTime_halo, dynTime_halo:dynTime_halo, $
       radCt:radCt, radDt:radDt, radDens:radDens, radTemp:radTemp, $
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
  if sP.trMCPerCell eq 0 then begin
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
  logMassBins   = [9.5,10.0,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11.0,$
                   11.1,11.25,11.5,11.75,11.9,13.1]
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

function modelMassFracs, sP=sP, ts=ts
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  tables = interpLambdaSD93()  

  ; config
  modelMassRange = [10.8,12.2]
  modelMassRes   = 150
  
  modelRadRange = [0.04,1.0]
  modelRadRes   = 400
  
  modelMasses = [1.0,2.0,5.0,$
                 10.0^linspace(modelMassRange[0],modelMassRange[1],modelMassRes)/1e10] ; code units

  ; check if save exists
  saveFilename = 'temp.mmf.sav'
  
  ; results exist, return
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif 
  
  r = { modelMassFrac_SIS1     : fltarr(modelMassRes)                            ,$
        modelMassFrac_SIS2     : fltarr(n_elements(ts.tsRatioVals),modelMassRes) ,$
        modelMassFrac_NFW_iso  : fltarr(n_elements(ts.tsRatioVals),modelMassRes) ,$
        modelMassFrac_NFW_poly : fltarr(n_elements(ts.tsRatioVals),modelMassRes) ,$
        modelMasses            : modelMasses    ,$
        modelMassRange         : modelMassRange ,$
        modelRadRes            : modelRadRes    ,$
        modelRadRange          : modelRadRange   }
  
  x = linspace(modelRadRange[0],modelRadRange[1],modelRadRes)
  
  for i=0,modelMassRes-1 do begin
    m_gas0 = modelMasses[i] * units.f_b
    
    ; MODEL 1: SIS uniform halo
    sis_dm  = sis_profile(1.0, mass=modelMasses[i], redshift=sP.redshift)
    sis_gas = sis_gas_profile(mass_hot=m_gas0, sis_dm=sis_dm, tables=tables)
    
    if sis_gas.coolTime lt sis_gas.dynTime then r.modelMassFrac_SIS1[i] = 1.0 ; all halo gas has tcool<tdyn
    if sis_gas.coolTime ge sis_gas.dynTime then r.modelMassFrac_SIS1[i] = 0.0 ; all halo gas has tcool>tdyn
    
    ; MODEL 2: SIS with tcool(r)
    sis_dm  = sis_profile(x, mass=modelMasses[i], redshift=sP.redshift)
    sis_gas = sis_gas_profile(mass_hot=m_gas0, sis_dm=sis_dm, tables=tables)

    foreach tsRatio,ts.tsRatioVals,j do begin
    
      if max(sis_gas.coolTime) lt tsRatio * sis_gas.dynTime then begin
        ; all halo gas has tcool<tdyn
        r.modelMassFrac_SIS2[j,i] = 1.0
      endif else begin
        if min(sis_gas.coolTime) gt tsRatio * sis_gas.dynTime then begin
          ; all halo gas has tcool>tdyn
          r.modelMassFrac_SIS2[j,i] = 0.0
        endif else begin
          ; find the radius where tcool=frac*tdyn
          w = where(sis_gas.coolTime ge tsRatio * sis_gas.dynTime,count)
          if count eq 0 then message,'error'
          w = min(w)
          e_rad = sis_dm.rad[w]
        
          ; calculate gas mass inside e_rad and take ratio
          r.modelMassFrac_SIS2[j,i] = (e_rad/sis_dm.r200)
        endelse
      endelse
    
    endforeach
    
    ; MODEL 3: NFW
    nfw_dm  = nfw_profile(x, mass=modelMasses[i], redshift=sP.redshift)
    nfw_gas = nfw_gas_suto(mass_hot=m_gas0, nfw_dm=nfw_dm, tables=tables)
    
    foreach tsRatio,ts.tsRatioVals,j do begin
    
      ; ISO
      zz = nfw_gas.coolTime_iso ge tsRatio * nfw_gas.dynTime
      zz_tot = total(zz)
      
      if zz_tot eq 0 then begin
        ; tcool is nowhere greater than tsRatio*tdyn (all halo gas has tcool<tdyn)
        r.modelMassFrac_NFW_iso[j,i] = 1.0
      endif else begin
        if zz_tot eq n_elements(modelRadRes) then begin
          ; all halo gas has tcool > tsRatio*tdyn
          r.modelMassFrac_NFW_iso[j,i] = 0.0
        endif else begin
          ; find the radius where tcool = tsRatio*tdyn
          w = where(nfw_gas.coolTime_iso ge tsRatio * nfw_gas.dynTime,count)
          if count eq 0 then message,'error'
          w = min(w)
          
          if w eq 0 then continue
        
          ; calculate gas mass inside e_rad and take ratio
          e_mass   = int_tabulated(nfw_dm.rad[0:w],4*!pi*nfw_dm.rad[0:w]^2.0*nfw_gas.rho_gas_iso[0:w])
          tot_mass = int_tabulated(nfw_dm.rad,4*!pi*nfw_dm.rad^2.0*nfw_gas.rho_gas_iso) ; equals m_gas0 only for poly

          r.modelMassFrac_NFW_iso[j,i] = (e_mass/tot_mass)
        endelse
      endelse
      
      ; POLY
      zz = nfw_gas.coolTime ge tsRatio * nfw_gas.dynTime
      zz_tot = total(zz)
      
      if zz_tot eq 0 then begin
        ; tcool is nowhere greater than tsRatio*tdyn (all halo gas has tcool<tdyn)
        r.modelMassFrac_NFW_poly[j,i] = 1.0
      endif else begin
        if zz_tot eq n_elements(modelRadRes) then begin
          ; all halo gas has tcool > tsRatio*tdyn
          r.modelMassFrac_NFW_poly[j,i] = 0.0
        endif else begin
          ; find the radius where tcool = tsRatio*tdyn
          w = where(nfw_gas.coolTime ge tsRatio * nfw_gas.dynTime,count)
          if count eq 0 then message,'error'
          w = min(w)
        
          ; calculate gas mass inside e_rad and take ratio
          e_mass   = int_tabulated(nfw_dm.rad[0:w],4*!pi*nfw_dm.rad[0:w]^2.0*nfw_gas.rho_gas[0:w])
          tot_mass = int_tabulated(nfw_dm.rad,4*!pi*nfw_dm.rad^2.0*nfw_gas.rho_gas_iso) ; equals m_gas0

          r.modelMassFrac_NFW_poly[j,i] = (e_mass/tot_mass)
        endelse
      endelse
      
    endforeach
    
  endfor ; modelMassRes
  
  ; make extremely low values not plot
  w = where(r.modelMassFrac_NFW_poly lt 0.0001,count)
  if count gt 0 then r.modelMassFrac_NFW_poly[w] = !values.f_nan
  w = where(r.modelMassFrac_NFW_iso lt 0.0001,count)
  if count gt 0 then r.modelMassFrac_NFW_iso[w] = !values.f_nan
  
  ; save
  ;save,r,filename=saveFilename
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
        r_cool[i] = sis_gas.r_cool
    
        ; compare rcool to rvir, "rapid mode" if rcool>rvir
        if sis_gas.r_cool ge sis_dm.r200 then rapid_mode[i] = 1B
      
        ; LU CODE
        lu_Hubble_h = 0.7
        lu_H0kpc = 0.1
        
        lu_rvir = sis_dm.r200
        lu_mvir = M_DM[i]
        lu_vvir = sqrt(4.30071e-6 * lu_mvir / lu_rvir)
        lu_tgas = 35.9 * lu_vvir^2.0
        
        lu_rho_0 = M_hot[i] / 4 / !pi / lu_rvir
        
        lu_cof = 1.0/0.28086 * units.mass_proton * units.boltzmann * lu_tgas / sis_gas.lambdaNet $ 
                 / units.UnitMass_in_g * units.UnitLength_in_cm^3.0
        lu_cof *= lu_H0kpc * lu_Hubble_h * units.UnitVelocity_in_cm_per_s / units.UnitLength_in_cm
        
        lu_dyn_time = lu_H0kpc * lu_rvir / lu_vvir
        lu_cool_time = lu_dyn_time
        lu_rcool = lu_rho_0 * lu_cool_time / lu_cof
        lu_rcool = sqrt(lu_rcool) * lu_Hubble_h
        ; END LU CODE
        
        ;print,lu_dyn_time/1e4,sis_gas.dynTime
        ;print,lu_rcool/sis_gas.r_cool
      
        ; set cold gas accretion rate based on mode (1e10 Msun/yr)
        if rapid_mode[i] then begin
          ; rapid mode
          dM_cool[i] = M_hot[i] / delta_t ;(sis_gas.dynTime*1e9)
        endif else begin
          ; slow mode
          ;dM_cool[i] = 0.5 * M_hot[i] * sis_gas.r_cool * (sis_dm.Vcirc_DM*units.s_in_yr/units.kpc_in_km) / sis_dm.r200^2.0 ;equivalent
          dM_cool[i] = 0.5 * M_hot[i] * sis_gas.r_cool / (sis_dm.r200 * sis_gas.dynTime*1e9)
        endelse
      
        if dM_cool[i] * delta_t gt M_hot[i] then begin
          print,'warning',i
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
      r_cool_kang[i] = sis_gas.r_cool_h
      
      ; compare rcool_h to rvir, "rapid mode" if rcool_h>rvir
      if sis_gas.r_cool_h ge sis_dm.r200 then rapid_mode_kang[i] = 1B
      
      ; set cold gas accretion rate based on mode
      kang_timenorm = times[i] ;sis_gas.hubbleTime*1e9 ;
      
      if rapid_mode_kang[i] then begin
        dM_cool_kang[i] = M_hot_kang[i] / kang_timenorm
      endif else begin
        dM_cool_kang[i] = 0.5 * M_hot_kang[i] * sis_gas.r_cool_h / (sis_dm.r200 * kang_timenorm)
      endelse
      
      if dM_cool_kang[i] * delta_t gt M_hot_kang[i] then message,'error'
      
      ; update cold gas
      M_cold_kang[i] = M_cold_kang[i-1] + dM_cool_kang[i] * delta_t
      
  endfor
  
  ; plot
  start_PS,'test_dMdt_vs_redshift.eps'
    ; halo baryons
    cgPlot,1+redshifts,dM_b*1e10,xtitle='1+z',ytitle='dM/dt [Msun/year]',$
      /xlog,xrange=[1,7],/xs,xminor=0,/ylog,yrange=[0.06,60.0],/ys,yminor=0,$
      xtickv=[1,2,3,4,5,6,7],xtickname=['1','2','3','4','5','6','7'],xticks=6,color=cgColor('black')
      
    ; croton model, dM/dt cold split into rapid and slow
    cgPlot,1+redshifts,dM_cool*1e10,/overplot,color=cgColor('red')
    
    w_rapid = where(rapid_mode eq 1B,count)
    cgPlot,1+redshifts[w_rapid],dM_cool[w_rapid]*1e10,psym=4,/overplot,color=cgColor('red')
    
    ; kang model
    cgPlot,1+redshifts,dM_cool_kang*1e10,/overplot,color=cgColor('blue')
    
    w_rapid = where(rapid_mode_kang eq 1B,count)
    cgPlot,1+redshifts[w_rapid],dM_cool_kang[w_rapid]*1e10,psym=4,/overplot,color=cgColor('blue')
    
    legend,['croton','kang'],textcolors=['red','blue'],box=0,/top,/right
  end_PS
  
  start_PS,'test_Mcold_vs_redshift.eps'
    ; halo baryons
    cgPlot,1+redshifts,M_b*1e10,xtitle='1+z',ytitle='Mcold [Msun]',$
      /xlog,xrange=[1,7],/xs,xminor=0,/ylog,yrange=[1.5e8,2e11],/ys,yminor=0,$
      xtickv=[1,2,3,4,5,6,7],xtickname=['1','2','3','4','5','6','7'],xticks=6,color=cgColor('black')
      
    ; croton model, dM/dt cold split into rapid and slow
    cgPlot,1+redshifts,M_cold*1e10,/overplot,color=cgColor('red')
    
    w_rapid = where(rapid_mode eq 1B,count)
    cgPlot,1+redshifts[w_rapid],M_cold[w_rapid]*1e10,psym=4,/overplot,color=cgColor('red')
    
    ; kang model
    cgPlot,1+redshifts,M_cold_kang*1e10,/overplot,color=cgColor('blue')
    
    w_rapid = where(rapid_mode_kang eq 1B,count)
    cgPlot,1+redshifts[w_rapid],M_cold_kang[w_rapid]*1e10,psym=4,/overplot,color=cgColor('blue')
    
    legend,['croton','kang'],textcolors=['red','blue'],box=0,/top,/right
  end_PS
  
  start_PS,'test_rads_vs_redshift.eps'
    cgPlot,[0],[0],/nodata,xtitle='1+z',ytitle='Radius [kpc]',$
      /xlog,xrange=[1,7],/xs,xminor=0,yrange=[0,200],/ys,$
      xtickv=[1,2,3,4,5,6,7],xtickname=['1','2','3','4','5','6','7'],xticks=6
      
    cgPlot,1+redshifts,R_vir,/overplot,color=cgColor('red')
    cgPlot,1+redshifts,R_cool,/overplot,color=cgColor('blue')
    legend,['$R_{vir}$','$R_{cool}$'],textcolors=['red','blue'],box=0,/top,/right
  end_PS    
  stop
end
  
; compVsMass(): at one redshift, compute models as a function of halo mass
; calculate critical halo mass crossover for rcool<>rvir
  
pro compVsMass
  
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()  
  
  ; config
  redshift = 2.0
  M0_range = [9.0,13.0] ; log(msun) at z=0
  
  hRes = 200
  hMasses = 10.0^(linspace(M0_range[0],M0_range[1],hRes))/1e10
  
  ; load Sutherland & Dopita (1993) cooling tables
  tables = interpLambdaSD93()
  
  ; arrays
  hMass_z = fltarr(hRes)
  virRad  = fltarr(hRes)
  virTemp = fltarr(hRes)
  
  sis_coolRad = fltarr(hRes)
  nfw_coolRad = fltarr(hRes)
  
  for i=0,hRes-1 do begin
    ; calculate halo masses / hot gas mass, at target redshift
    hMass_z[i] = haloMAH(redshift,M0=hMasses[i],z0=0.0)
    m_gas0 = hMass_z[i] * units.f_b
    
    ; SIS DM halo and gas profiles
    sis_dm  = sis_profile(1.0, mass=hMass_z[i], redshift=redshift)
    sis_gas = sis_gas_profile(mass_hot=m_gas0, sis_dm=sis_dm, tables=tables)

    nfw_dm  = nfw_profile(linspace(0.01,10.0,1000), mass=hMass_z[i], redshift=redshift)
    nfw_gas = nfw_gas_profile(mass_hot=m_gas0, nfw_dm=nfw_dm, tables=tables)
    
    ; radii, code units (kpc)
    virRad[i]  = sis_dm.r200
    virTemp[i] = sis_dm.Tvir_DM
    
    sis_coolRad[i] = sis_gas.r_cool
    nfw_coolRad[i] = nfw_gas.r_cool
  endfor
  
  ; plot
  start_PS,'test_rads_vs_mass.eps'
    cgPlot,[0],[0],/nodata,xtitle='Halo Mass at z=2',ytitle='Radius [kpc]',$
      xrange=[8.5,12],/xs,yrange=[0,220],/ys
      
    ; mark tvir=10^4 (end of cooling tables)
    w = min(where(alog10(virTemp) ge 4.0))
    cgPlot,codeMassToLogMsun([hMass_z[w],hMass_z[w]]),[50,200],/overplot,line=1
      
    ; rvir and rcool
    cgPlot,codeMassToLogMsun(hMass_z),virRad,/overplot,color=cgColor('black'),line=2
    
    cgPlot,codeMassToLogMsun(hMass_z),sis_coolRad,/overplot,color=cgColor('red'),line=0
    cgPlot,codeMassToLogMsun(hMass_z),nfw_coolRad,/overplot,color=cgColor('blue'),line=0
    
    legend,textoidl(['R_{vir}','SIS R_{cool}','NFW R_{cool}']),$
      textcolors=['black','red','blue'],linestyle=[2,0,0],box=0,/top,/right,linesize=0.25
  end_PS
  
end
  
; compProfiles(): for one redshift and halo mass, compare different theoretical profiles
;                 density,temp,timescales vs radius
  
pro compProfiles
  
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()    
  
  ; config
  redshift = 2.0
  hMass    = 100.0 ; code units at redshift
  
  ; load Sutherland & Dopita (1993) cooling tables
  tables = interpLambdaSD93()
  
  ; arrays
  rRes = 200
  rPts = linspace(0.01,1.5,rRes) ; r/rvir
  
  ; DM and gas profiles
  sis_dm  = sis_profile(rPts, mass=hMass, redshift=redshift)
  sis_gas = sis_gas_profile(mass_hot=hMass*units.f_b, sis_dm=sis_dm, tables=tables)
  
  nfw_dm  = nfw_profile(rPts, mass=hMass, redshift=redshift)
  nfw_gas = nfw_gas_profile(mass_hot=hMass*units.f_b, nfw_dm=nfw_dm, tables=tables)
  
  nfw_gas_poly = nfw_gas_poly(mass_hot=hMass*units.f_b, nfw_dm=nfw_dm, tables=tables)
  
  ; suto model
  nfw_gas_suto2 = nfw_gas_suto(mass_hot=hMass*units.f_b, nfw_dm=nfw_dm, tables=tables)
  
  T_0   = 1.2 * nfw_dm.Tvir_rs
  rho_0 = 150 * nfw_dm.mass_codeunits / (4/3*!pi*nfw_dm.r200^3.0)
  n     = 20.0
  
  nfw_gas_suto = suto_model(T_0,rho_0,n,nfw_dm)
  
  print,'final',alog10(T_0),alog10(rho_0),nfw_gas_suto.B_p
  
  ; plot
  start_PS,'halo_ts_vs_rad.eps'
    cgPlot,[0],[0],/nodata,xtitle=textoidl(' r / r_{vir} '),ytitle='Timescale [Gyr]',$
      xrange=[0.0,1.5],/xs,yrange=[0,5],/ys
      
    cgPlot,rPts,sis_gas.coolTime,/overplot,color=cgColor('red'),line=0
    cgPlot,rPts,sis_gas.dynTime,/overplot,color=cgColor('red'),line=2
    cgPlot,rPts,nfw_gas.coolTime,/overplot,color=cgColor('blue'),line=0
    cgPlot,rPts,nfw_gas.dynTime,/overplot,color=cgColor('blue'),line=2
    
    legend,textoidl(['SIS \tau_{cool}','SIS \tau_{dyn}','NFW \tau_{cool}','NFW \tau_{dyn}']),$
      textcolors=['red','red','blue','blue'],linestyle=[0,2,0,2],box=0,/top,/left,linesize=0.25
  end_PS
  
  ; plot density/temp profiles
  start_PS,'halo_denstemp_vs_rad.eps', ys=8, xs=6
    cgPlot,[0],[0],/nodata,xtitle='',ytitle='Density (ratioToCrit)',$
      xrange=[0.01,1.5],/xs,/xlog,xminor=0,yrange=[1.1,8],/ys,$
      position = [0.15,0.55,0.95,0.95],xtickname=replicate(' ',10)
      
    ; DM profile
    cgPlot,rPts,alog10(rhoRatioToCrit(sis_dm.rho_DM,redshift=redshift)),/overplot,color=cgColor('red')
    cgPlot,rPts,alog10(rhoRatioToCrit(nfw_dm.rho_DM,redshift=redshift)),/overplot,color=cgColor('black')
    
    ; r_s line for NFW
    cgPlot,[nfw_dm.r_s/nfw_dm.r200,nfw_dm.r_s/nfw_dm.r200],[1.5,3],line=2,/overplot
    
    ; gas profile
    cgPlot,rPts,alog10(rhoRatioToCrit(sis_gas.rho_gas,redshift=redshift)),/overplot,line=2,color=cgColor('red')
    cgPlot,rPts,alog10(rhoRatioToCrit(nfw_gas.rho_gas,redshift=redshift)),/overplot,line=2,color=cgColor('blue')
    cgPlot,rPts,alog10(rhoRatioToCrit(nfw_gas_suto.rho_gas,redshift=redshift)),/overplot,line=2,color=cgColor('orange')
    cgPlot,rPts,alog10(rhoRatioToCrit(nfw_gas_suto.rho_gas_iso,redshift=redshift)),/overplot,line=2,color=cgColor('magenta')
    cgPlot,rPts,alog10(rhoRatioToCrit(nfw_gas_suto2.rho_gas,redshift=redshift)),/overplot,line=2,color=cgColor('green')
    
    legend,['SIS','NFW iso','suto poly','suto iso','autoPoly'],textcolors=['red','blue','orange','magenta','green'],box=0,/right,/top
    
    cgPlot,[0],[0],/nodata,xtitle=textoidl('r / r_{vir}'),ytitle='Temp (logK)',$
       xrange=[0.01,1.5],/xs,/xlog,xminor=0,yrange=[5,7],/ys,/noerase,position = [0.15,0.15,0.95,0.55]
      
    cgPlot,rPts,alog10(sis_gas.temp_gas),/overplot,color=cgColor('red')
    cgPlot,rPts,alog10(nfw_gas.temp_gas),/overplot,color=cgColor('blue')
    cgPlot,rPts,alog10(nfw_gas_suto.temp_gas),/overplot,color=cgColor('orange')
    cgPlot,rPts,alog10(nfw_gas_suto2.temp_gas),/overplot,color=cgColor('green')
    
  end_PS
  
  stop
end

