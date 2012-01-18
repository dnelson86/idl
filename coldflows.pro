; coldflows.pro
; cold flows - main
; dnelson nov.2011

; ---------------------------------------------------------------------------------
; ---------------------------------------------------------------------------------
; NOTE: Nothing in this file is still in use and may be out of date and/or buggy.
;       Only still around for ideas when working on cosmo* codebase.
; ---------------------------------------------------------------------------------
; ---------------------------------------------------------------------------------

@helper
@cosmoUtil
@cosmoLoad
@coldflowsVis

; findAccretionRate(): calculate smooth accretion rate at every snapshot split by hot/cold mode
;

function findAccretionRate, res=res, bSH=bSH, halos=halos

  if not keyword_set(res) then begin
     print,'Error: findAccretionRate() arguments not specified.'
     return,0
  endif

  ; config
  gadgetPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'  
  dataPath     = '/n/home07/dnelson/coldflows/thermhist/'
  workingPath  = '/n/home07/dnelson/coldflows/thermhist.deriv/'  
  
  snapRange = [51,314]

  critLogTemp = 5.5  
  
  saveFilename = workingPath + 'accretion.rate.' + str(res) + '.' + str(snapRange[0]) + '-' + $
                 str(snapRange[1]) + '.sav'
  
  if (keyword_set(halos)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.halos.sav'
  if (keyword_set(bSH)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.bSH.sav'  
  
  if (file_test(saveFilename)) then begin
    restore,saveFilename
    r = {smoothAccHot:smoothAccHot,smoothAccCold:smoothAccCold}
    return,r
  endif

  ; hot mode (temp > critLogTemp) mask and smooth accretion  mask
  h = loadSnapshotHeader(gadgetPath,snapNum=0)
  nTot = h.nPartTot[0] + h.nPartTot[1] ; gas + dm (no stars at snap 0)
  
  hmMask = bytarr(nTot+1) ; 0 (default) = temp has not exceeded critLogTemp at an point in the past
                          ; 1 = temp has exceeded critLogTemp at some point in the past
   
  sgMask = bytarr(nTot+1) ; 0 (default) = not found at any previous time in any bound substructure
                          ; 1 = found at some previous time in a subfind group, not "smooth"   
   
  ; result arrays
  smoothAccHot  = fltarr(snapRange[1]+1)
  smoothAccCold = fltarr(snapRange[1]+1)
   
  ; loop over starting snapshots to find maxtemps during this period
  for m=0,snapRange[0]-1 do begin
    maxTempsFilename = dataPath + 'thermhist.gas.'+str(res)+'_' + str(m) + '.sav'
    
    ; load maximum temperatures
    restore,maxTempsFilename
    density = !NULL
    
    ; take log
    w = where(temp eq 0,count)
    temp[w] = 1.0
    
    temp = alog10(temp)
    
    ; update hot mode mask
    w = where(temp ge critLogTemp,count)
    temp = !NULL
    
    if (count ge 1) then $
      hmMask[w] = 1B
    
    w = !NULL
  endfor
  
  ; constant masses (GADGET ONLY)
  masses = loadSnapshotSubset(gadgetPath,snapNum=snapRange[0],partType='gas',field='mass')
  masses = masses[0]
  
  ; loop over all for accretion rate calculation
  for m=snapRange[0],snapRange[1] do begin
    ; find smoothly accreted gas
    sgIDs_All = getSubgroupIDList(res=res,snap=m,bSH=bSH,halos=halos,/gasOnly)
    
    sgCur_pIDs = where(sgMask eq 1B,countMask)
    sgIDs_All = removeIntersectionFromB(sgCur_pIDs,sgIDs_All)    
    sgCur_pIDs = !NULL
   
    ; select cold as those not in mask
    mask_pIDs  = where(hmMask eq 1B,countMask)
    sgIDs_Hot  = 1 ; set as union return
    sgIDs_Cold = removeIntersectionFromB(mask_pIDs,sgIDs_All,union=sgIDs_Hot)
    
    sgIDs_All  = !NULL
    mask_pIDs  = !NULL
    
    ; load masses and pids
    ; mass+id load for arepo
    ; ids    = loadSnapshotSubset(gadgetPath,snapNum=m,partType='gas',field='ids')
    
    ; store (arepo w/ masses)
    ;smoothAccCold[m] = total(masses[ids[sgIDs_Cold]])
    
    ;if (sgIDs_Hot[0] ne 1) then $
    ;  smoothAccHot[m]  = total(masses[ids[sgIDs_Hot]])
      
    ; store (GADGET ONLY)
    smoothAccCold[m] = masses * n_elements(sgIDs_Cold)
    smoothAccHot[m]  = masses * n_elements(sgIDs_Hot)
    
    sgIDs_Cold = !NULL
    sgIDs_Hot  = !NULL
    
    print,'snap ['+str(m)+'] cold = ' + str(smoothAccCold[m]) + ' hot = ' + str(smoothAccHot[m])
    
    ; masses = !NULL ; for arepo
    ids = !NULL
    
    ; load maximum temperatures
    maxTempsFilename = dataPath + 'thermhist.gas.'+str(res)+'_' + str(m) + '.sav'
    
    ; load maximum temperatures
    restore,maxTempsFilename
    density = !NULL
    
    ; take log
    w = where(temp eq 0,count)
    temp[w] = 1.0
    
    temp = alog10(temp)
    
    ; update hot mode mask
    w = where(temp ge critLogTemp,count)
    temp = !NULL
    
    hmMask[w] = 1B
    
    w = !NULL
    
    ; load group catalog for current snapshot m
    sgCur_pIDs = getSubgroupIDList(res=res,snap=m,bSH=bSH) ; do not include /halos
    
    ; add particles to exclusion mask (modulo bSH and halos)
    sgMask[sgCur_pIDs] = 1B
    sgCur_pIDs = !NULL
  endfor
  
  ;save
  save,smoothAccHot,smoothAccCold,filename=saveFilename
  print,'Saved: ' + saveFilename
  
  r = {smoothAccHot:smoothAccHot,smoothAccCold:smoothAccCold}
  return,r
                          
end

; findAccretedGas(): find gas particle ids from subhalos in targetSnap that were either:
;                    (a) not members of any subhalo at one specified earlier snapshot (smoothCutSnap)
;                    (b) not members of any subhalo for all previous snapshots (smoothCutSnap not set)
;                    [background subhalos optionally included or not through bSH flag]
;
; zWidth=0 (accretion only from immediately preceeding snap if redshiftsCut not set, or if
;           redshiftsCut is set, then accretion allowed with the implicit redshift width
;           back until the cut, with no checks in the middle to insure that gas particles
;           which are not a member of any subhalo at the cut redshift are also not members
;           of any subhalo at any redshift prior to that when they were found in a subhalo
;           at the target redshift)
; zWidth=1 : allow a number of snapshots counting back from the target through which
;            smooth accretion is allowed (resolved in any one snapshot but not in the previous).
;            After the width is ended unresolved in all previous snapshots is enforced
;
; halos=1: find cosmological smooth accretion onto main subhalos ("background" subfind groups)
;  - should not set bSH, smoothCutSnap, zWidth for most recent approach
; DM=1 : dark matter instead of gas

function findAccretedGas, res=res, bSH=bSH, targetSnap=targetSnap, sgTarget=sgTarget, $
                          smoothCutSnap=smoothCutSnap, zWidth=zWidth, halos=halos, DM=DM

  if not keyword_set(res) or not keyword_set(targetSnap) then begin
     print,'Error: findAccretedGas() arguments not specified.'
     return,0
  endif
  if keyword_set(zWidth) and not keyword_set(smoothCutSnap) then begin
      print,'Error: Cannot use zWidth without a smoothCutSnap.'
      return,0
  endif
  
  ; config
  gadgetPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'  
  workingPath  = '/n/home07/dnelson/coldflows/thermhist.deriv/'  
  
  ; set range for enforcing smooth accretion
  if keyword_set(smoothCutSnap) then begin
      smoothStart = smoothCutSnap
      smoothEnd   = smoothCutSnap
  endif else begin
      smoothStart = targetSnap - 1
      smoothEnd   = 50 ; first group catalogs at z=6
  endelse
  
  ; set range for allowing a width of smooth accretion redshifts
  if keyword_set(zWidth) then begin
    widthStart = targetSnap
    widthEnd   = smoothCutSnap
    smoothStart = widthEnd - 1
    smoothEnd   = 50
  endif else begin
    widthStart = targetSnap
    widthEnd   = targetSnap
  endelse
  
  if keyword_set(DM) then partName = 'dm'
  if not keyword_set(DM) then partName = 'gas'
  
  ; load group catalog at targetSnap if not passed in
  if not keyword_set(sgTarget) then $
    sgTarget  = loadSubhaloGroups(gadgetPath,targetSnap)
  
  ; set saveFilename and check existence
  saveFilename = workingPath + 'accreted.'+partName+'.'+str(res)+'.target='+str(targetSnap)+'.'+$
                  str(smoothStart)+'-'+str(smoothEnd)+'.sav'
                  
  if (keyword_set(halos)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.halos.sav'
  if (keyword_set(zWidth)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.zW.sav'
  if (keyword_set(bSH)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.bSH.sav'
                  
  if (file_test(saveFilename)) then begin
    restore,saveFilename
    return,sgIDs_All
  endif

  ; arrays to hold all valid particle ids
  sgIDs_All = []
  
  for m=widthStart,widthEnd,-1 do begin

    ; get list of gas ids in subhalos (or halos) at this snap
    if not keyword_set(DM) then $
      sgIDs_Cur = getSubgroupIDList(res=res,snap=m,bSH=bSH,halos=halos,/gasOnly)
    if keyword_set(DM) then $
      sgIDs_Cur = getSubgroupIDList(res=res,snap=m,bSH=bSH,halos=halos,/dmOnly)
    
    ; enforce smooth accretion on immediately prior if using zWidth:
    if keyword_set(zWidth) then begin
      ; make (exclusion) list of particles belonging to subhalos
      sgCur_pIDs = getSubgroupIDList(res=res,snap=m-1,bSH=bSH) ; do not include /halos
      sgIDs_Cur  = removeIntersectionFromB(sgCur_pIDs,sgIDs_Cur)
      sgCur_pIDs = !NULL
      print,' After zWidth m='+str(m)+' smooth accretion cut have ['+str(n_elements(sgIDs_Cur))+'] left.'
    endif
    
    ; add this part of zWidth (or single selection at targetSnap) to keeper array
    sgIDs_All = [sgIDs_All,sgIDs_Cur]
    sgIDs_Cur = !NULL
    
  endfor ;m
  
  ; enforce smooth accretion to beginning of simulation
  print,'Running smooth accretion cut ('+str(n_elements(sgIDs_All))+' particles) over snapshots ['+$
        str(smoothStart)+'-'+str(smoothEnd)+'].'
        
  for m=smoothStart,smoothEnd,-1 do begin
    ; make (exclusion) list of particles belonging to subhalos
    sgCur_pIDs = getSubgroupIDList(res=res,snap=m,bSH=bSH) ; do not include /halos
    sgIDs_All  = removeIntersectionFromB(sgCur_pIDs,sgIDs_All)
  endfor ;m
    
  print,'After smooth accretion cut have ['+str(n_elements(sgIDs_All))+'] left.'
  
  save,sgIDs_All,filename=saveFilename
  
  return, sgIDs_All
  
end

; saveThermalHistories(): record the density+temperature for all gas particles at all snapshots

pro saveThermalHistories, res=res
  
  if not keyword_set(res) then begin
    print,'Error: Must specify resolution.'
    return
  endif
  
  ; config
  gadgetPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'
  workingPath  = '/n/hernquistfs1/dnelson/coldflows/thermhist/'
  saveFilebase = workingPath + 'thermhist.gas.'+str(res)+'_'
  
  nSnaps = n_elements(file_search(gadgetPath+'snapdir_*')) - 1 ;exclude final duplicate
  
  ; load first snapshot header
  h = loadSnapshotHeader(gadgetPath,snapNum=0)
  nGas = h.nPartTot[0]
  
  print,'Saving ['+str(nGas)+'] gas particles (rho,T) over ['+str(nSnaps)+'] snapshots.'
  
  ; loop from target snapshot through all prior snapshots
  for m=0,nSnaps-1,1 do begin
  
    ; check we haven't already done this
    saveFilename=saveFilebase + str(m) + '.sav'
    
    if (file_test(saveFilename)) then begin
      print,' Skipping: ' + saveFilename
      continue
    endif
  
    ; arrays
    temp    = fltarr(nGas+1)
    density = fltarr(nGas+1)
    
    ; setup such that first array index = particle ID, so index=0 is undefined
    temp[0]    = -1.0
    density[0] = -1.0
  
    if (m mod 10 eq 0) then print, ' '+str(m)
    ; load u and nelec and calculate temperature
    u     = loadSnapshotSubset(gadgetPath,snapNum=m,partType='gas',field='u')
    nelec = loadSnapshotSubset(gadgetPath,snapNum=m,partType='gas',field='ne')
    
    t = convertUtoTemp(u,nelec)
    
    u     = !NULL
    nelec = !NULL
    
    ; load density
    rho   = loadSnapshotSubset(gadgetPath,snapNum=m,partType='gas',field='density')
    
    ; load ids for sorted insert
    ids = loadSnapshotSubset(gadgetPath,snapNum=m,partType='gas',field='ids')
    
    ; save properties
    temp[ids]    = t
    density[ids] = rho

    rho = !NULL
    t   = !NULL
    ids = !NULL
    
    ; save
    save,temp,density,nGas,nSnaps,filename=saveFileName
    
  endfor

end

; calcRhoTemp2DHisto(): calculate/save mass-weighted 2d histogram of thermal evolution tracks 
;                       in (rho,temp) plane in a redshift range [zMin,zMax], separated into
;                       all tracks, hot mode tracks, and cold mode tracks
;
; timeWeight=1 : convert snapshot spacing to Gyr and use to weight mass contributions to each bin
;                such that the contribution of each gas cell to each bin scales with the amount of
;                time it spends in that bin. if set to zero, the "total mass" units of each bin
;                are biased by the snapshot time sampling (uneven) of the (rho,temp) plane

function calcRhoTemp2DHisto, bSH=bSH, halos=halos, zMin=zMin, zMax=zMax, res=res, nbins=nbins, tW=tW

  units = getUnits()

  if not keyword_set(res) or not keyword_set(nbins) then begin
    print,'Error: Must specific resolution set.'
    return,0
  endif
  
  critLogTemp   = 5.5  ; cold/hot mode separator    
  
  if not keyword_set(zMin) then zMin = 0.0
  if not keyword_set(zMax) then zMax = 30.0  
  
  ; restrict arrays to redshift range
  maxSnap = redShiftToSnapnum(zMin)
  minSnap = redShiftToSnapnum(zMax)

  ; config
  gadgetPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'
  dataPath     = '/n/hernquistfs1/dnelson/coldflows/thermhist/'
  workingPath  = '/n/home07/dnelson/coldflows/thermhist.deriv/'
  
  ; set saveFilename
  saveFilename = workingPath + 'rhot.2dhisto.nB='+str(nbins)+'.'+str(res)+'.'+str(minSnap)+'-'+str(maxSnap)+'.sav'
  
  if (keyword_set(halos)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.halos.sav'  
  if (keyword_set(bSH)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.bSH.sav'
  if (keyword_set(tW)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.tW.sav'  

  if (file_test(saveFilename)) then begin
    restore,saveFilename
    r = {h2rt:h2rt,h2rt_hot:h2rt_hot,h2rt_cold:h2rt_cold}
    return,r
  endif
  
  print,'Calculating new (rho,temp) 2d histo for res = '+str(res)+' in range ['+str(minSnap)+'-'+str(maxSnap)+'].'
  
  ; get list of smoothly accreted gas particle ids
  sgIDs_Acc = findAccretedGas(res=res,bSH=bSH,targetSnap=maxSnap,halos=halos)  
  
  ; load masses (assume all are equal and constant in time)
  masses = loadSnapshotSubset(gadgetPath,snapNum=maxSnap,partType='gas',field='mass')
  
  if (min(masses) ne max(masses)) then begin
    print,'ERROR' ;gadget check
    return,0
  end
  
  ;masses = masses[sgIDs_Acc] ;Arepo
  masses = masses[0] ;Gadget (assume all are equal)
  
  h2rt      = fltarr(nbins,nbins)
  h2rt_hot  = fltarr(nbins,nbins)
  h2rt_cold = fltarr(nbins,nbins)
  
  ; time spacing
  redshifts = snapNumToRedshift(/all)
  
  ; load maximum temperatures
  maxTemps = maxTemperatures(res=res, zMin=zMin)
  maxTemps = maxTemps[sgIDs_Acc]
  
  wH = where(maxTemps ge critLogTemp,countH,complement=wC,ncomplement=countC)
  
  ; load thermal history
  for m=maxSnap,minSnap,-1 do begin
    thFilename = dataPath + 'thermhist.gas.'+str(res)+'_'+str(m)+'.sav'
    
    restore,thFilename

    ; restrict to desired gas particles
    density = density[sgIDs_Acc]
    temp    = temp[sgIDs_Acc]

    ; avoid log(0)
    w = where(temp eq 0,count)
    temp[w] = 1.0
    
    ; calculate time spacing to previous snapshot, use as weight if requested
    if keyword_set(tW) then begin
      curTime  = redshiftToAge(redshifts[m])
      prevTime = redshiftToAge(redshifts[(m-1) + ((m-1) lt 0)]) ; zero weight for snap=0
      timeWeight = curTime - prevTime
    endif else begin
      timeWeight = 1.0
    endelse
    
    ; calculate total histogram
    h2rt_cur = rhoTHisto(density,temp,nbins=nbins-1) ;mass=masses for Arepo
    h2rt += (h2rt_cur[0:nbins-1,0:nbins-1] * timeWeight)
    
    ; calculate hot mode histogram
    if (countH ne 0) then begin
      h2rt_cur = rhoTHisto(density[wH],temp[wH],nbins=nbins-1) ;mass=masses for Arepo
      h2rt_hot += (h2rt_cur[0:nbins-1,0:nbins-1] * timeWeight)
    endif
    
    ; calculate cold mode histogram
    if (countC ne 0) then begin
      h2rt_cur = rhoTHisto(density[wC],temp[wC],nbins=nbins-1) ;mass=masses for Arepo
      h2rt_cold += (h2rt_cur[0:nbins-1,0:nbins-1] * timeWeight)
    endif
  endfor ;m
  
  ; for GADGET, multiply by constant mass
  h2rt      *= masses
  h2rt_hot  *= masses
  h2rt_cold *= masses
  
  ; save results
  save,h2rt,h2rt_hot,h2rt_cold,filename=saveFilename
  print,'Saved: '+saveFilename

  r = {h2rt:h2rt,h2rt_hot:h2rt_hot,h2rt_cold:h2rt_cold} ;total mass (code units) per bin
  return, r
end

; calcMaxTempVsMassHisto

function calcMaxTempVsMassHisto, res=res, halos=halos, targetSnap=targetSnap, smoothCutSnap=smoothCutSnap

  units = getUnits()

  if (not keyword_set(res) or not keyword_set(targetSnap)) then begin
    print,'Error: calcMaxTempVsMassHisto: Inputs not set.'
    return,0
  endif
  
  ; config
  gadgetPath  = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'
  workingPath  = '/n/home07/dnelson/coldflows/thermhist.deriv/'

  massMM = [9.0,13.0] ; minmax for mass
  tempMM = [4.0,7.0]  ; minmax for temp
  massB  = 0.025
  tempB  = 0.025
  
  ; check savefile for existence
  if (keyword_set(smoothCutSnap)) then $
    endingSnap = smoothCutSnap
  if (not keyword_set(smoothCutSnap)) then $
    endingSnap = 'all'
    
  saveFilename = workingPath + 'tmaxmass.'+str(res)+'.'+str(targetSnap)+'-'+str(endingSnap)+'.sav'
  
  if (keyword_set(halos)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.halos.sav' 

  if (file_test(saveFilename)) then begin
    restore,saveFilename
    r = {shmass:shmass,tm2d:tm2d}
    return,r
  endif
  
  ; calc max temp
  redshift = snapNumToRedshift(snapNum=targetSnap)
  maxTemps = maxTemperatures(res=res, zMin=redshift)
  
  ; get smoothly accreted gas particle id list
  sgTarget  = loadSubhaloGroups(gadgetPath,targetSnap)
  sgIDs_Acc = findAccretedGas(res=res,bSH=bSH,targetSnap=targetSnap,sgTarget=sgTarget,$
                              smoothCutSnap=smoothCutSnap,halos=halos)  

  maxTemps = maxTemps[sgIDs_Acc]

  ; make list of sgIDs of either subhalos or halos
  valSGids = getPrimarySubhaloList(sgTarget,halos=halos)
  shmass   = sgTarget.subgroupMass[valSGids]
  ;shmass   = sgTarget.groupMass[sgTarget.subgroupGrnr[valSGids]] ;use/record parent fof masses
  
  ; reform subgroupIDs to valid & acc only
  subgroupIDs = lonarr(total(sgTarget.subgroupLen[valSGids]))
  subgroupOff = lonarr(n_elements(valSGids))
  subgroupLen = lonarr(n_elements(valSGids))
  
  match,sgTarget.subgroupIDs,sgIDs_Acc,sgtarg_ind,sgIDs_ind,count=matchAcc
  
  nextOff = 0
  
  for i=0,n_elements(valSGids)-1 do begin
    origOff = sgTarget.subGroupOffset[valSGids[i]]
    origLen = sgTarget.subGroupLen[valSGids[i]]
    
    w = where(sgtarg_ind ge origOff and sgtarg_ind lt (origOff+origLen),count)
    
    subgroupLen[i] = count
    subgroupOff[i] = nextOff
    nextOff = subgroupOff[i] + subgroupLen[i]
    ;print,i,count,subgroupOff[i],nextOff-1
    if (count ne 0) then $
      subgroupIDs[subgroupOff[i]:nextOff-1] = sgIDs_Acc[sgIDs_ind[w]]
    
  endfor
  
  sgTarget = !NULL
  
  ; 2d histogram
  massNB = (massMM[1]-massMM[0])/massB
  tempNB = (tempMM[1]-tempMM[0])/tempB
  
  tm2d = fltarr(massNB,tempNB)
  
  print,'Processing ['+str(n_elements(valSGids))+'] halos.'

  for i=0,n_elements(valSGids)-1 do begin
    if (i mod 1000 eq 0) then print,' '+str(i)
    
    ; method with reformed arrays
    if (subgroupLen[i] ne 0) then begin
      pIDs = subgroupIDs[subgroupOff[i] : subgroupOff[i] + subgroupLen[i] - 1]
      match,sgIDs_Acc,pIDs,sgIDs_ind,pIDs_ind,count=count
    endif

    if (count eq 0) then $
      continue
    
    ; 2d histogram
    hist = histogram(maxTemps[sgIDs_ind],binsize=tempB,min=tempMM[0],max=tempMM[1])
    hist = hist[0:tempNB-1]
    mass_ind = alog10( shmass[i] * (units.UnitMass_in_g / units.Msun_in_g) )
    mass_ind = floor((mass_ind - massMM[0]) / massB)
    
    ; out of mass range we're saving
    if (mass_ind lt 0 or mass_ind ge massNB) then $
      continue
    
    tm2d[mass_ind,*] += hist
    
  endfor
  
  ; subhalo masses (total baryonic+DM in log(msun))
  shmass = alog10( (units.UnitMass_in_g / units.Msun_in_g) * shmass )  
  
  ; save
  save,shmass,tm2d,valSGids,filename=saveFilename
  
  r = {shmass:shmass,tm2d:tm2d}
  return,r
  
end

; calcColdFracVsSubhaloMass()

function calcColdFracVsSubhaloMass, res=res, halos=halos, bSH=bSH, critLogTemp=critLogTemp, $
                                    targetSnap=targetSnap, smoothCutSnap=smoothCutSnap

  units = getUnits()

  if (not keyword_set(res) or not keyword_set(targetSnap) or not keyword_set(critLogTemp)) then begin
    print,'Error: calcColdFracVsSubhaloMass: Inputs not set.'
    return,0
  endif
  
  if (str(critLogTemp) eq 'vTC') then begin
    critLogTempStr = 'vTC'
  endif else begin
    critLogTempStr = string(critLogTemp,format='(f3.1)')
  endelse
  
  ; config
  gadgetPath  = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'
  workingPath  = '/n/home07/dnelson/coldflows/thermhist.deriv/'
  
  minNumGasPart = 6    ; minimum smoothly accreted gas particles per subhalo
  
  ; check savefile for existence
  if (keyword_set(smoothCutSnap)) then $
    endingSnap = smoothCutSnap
  if (not keyword_set(smoothCutSnap)) then $
    endingSnap = 'all'
    
  saveFilename = workingPath + 'coldfracmass.'+str(res)+'.CT='+critLogTempStr+$
                 '.'+str(targetSnap)+'-'+str(endingSnap)+'.sav'
  
  if (keyword_set(halos)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.halos.sav'  
  if (keyword_set(bSH)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.bSH.sav'

  if (file_test(saveFilename)) then begin
    restore,saveFilename
    r = {shmass:shmass,coldfrac:coldfrac,valSGids:valSGids,counts:counts}
    return,r
  endif
  
  ; calc max temp
  redshift = snapNumToRedshift(snapNum=targetSnap)
  maxTemps = maxTemperatures(res=res, zMin=redshift)
  
  ; get smoothly accreted gas particle id list
  sgTarget  = loadSubhaloGroups(gadgetPath,targetSnap)
  sgIDs_Acc = findAccretedGas(res=res,bSH=bSH,targetSnap=targetSnap,sgTarget=sgTarget,$
                              smoothCutSnap=smoothCutSnap,halos=halos)  

  maxTemps = maxTemps[sgIDs_Acc]

  ; make list of sgIDs of either subhalos or halos
  valSGids = getPrimarySubhaloList(sgTarget,halos=halos)
  shmass   = sgTarget.subgroupMass[valSGids]
  ;shmass   = sgTarget.groupMass[sgTarget.subgroupGrnr[valSGids]] ;use/record parent fof masses
  
  ; reform subgroupIDs to valid & acc only
  subgroupIDs = lonarr(total(sgTarget.subgroupLen[valSGids]))
  subgroupOff = lonarr(n_elements(valSGids))
  subgroupLen = lonarr(n_elements(valSGids))
  
  match,sgTarget.subgroupIDs,sgIDs_Acc,sgtarg_ind,sgIDs_ind,count=matchAcc
  
  nextOff = 0
  
  for i=0,n_elements(valSGids)-1 do begin
    origOff = sgTarget.subGroupOffset[valSGids[i]]
    origLen = sgTarget.subGroupLen[valSGids[i]]
    
    w = where(sgtarg_ind ge origOff and sgtarg_ind lt (origOff+origLen),count)
    
    if (count gt 0) then begin
      subgroupLen[i] = count
      subgroupOff[i] = nextOff
  
      nextOff = subgroupOff[i] + subgroupLen[i]
      ;print,i,count,subgroupOff[i],nextOff-1  

      subgroupIDs[subgroupOff[i]:nextOff-1] = sgIDs_Acc[sgIDs_ind[w]]
    endif
    
  endfor
  
  sgTarget = !NULL
  
  ; subhalo cold fractions (y-axis)
  coldfrac = fltarr(n_elements(valSGids))
  counts   = lonarr(n_elements(valSGids))
  
  print,'Processing ['+str(n_elements(valSGids))+'] halos.'

  for i=0,n_elements(valSGids)-1 do begin
    if (i mod 1000 eq 0) then print,' '+str(i)
    
    ; method without reformed subgroupIDs,subgroupOff,subgroupLen (much slower)
    ;sgIDs = sgTarget.subGroupIDs[sgTarget.subGroupOffset[valSGids[i]] : $
    ;                             sgTarget.subGroupOffset[valSGids[i]] + sgTarget.subGroupLen[valSGids[i]] - 1]
    ;match,sgIDs,sgIDs_Acc,sh_ind,sgIDs_ind,count=matchCount
    ;counts[i] = matchCount
    
    ; method with reformed arrays
    if (subgroupLen[i] ne 0) then begin
      pIDs = subgroupIDs[subgroupOff[i] : subgroupOff[i] + subgroupLen[i] - 1]
      match,sgIDs_Acc,pIDs,sgIDs_ind,pIDs_ind,count=count
      counts[i] = subgroupLen[i]
    endif

    ; common
    if (counts[i] eq 0 or counts[i] lt minNumGasPart) then begin
      coldfrac[i] = -1
      continue
    endif
    
    ; select count based on logCritTemp or halo virial temp
    if (str(critLogTemp) eq 'vTC') then begin
      critT = alog10( codeMassToVirTemp(shmass[i],redshift) )
      wCold = where(maxTemps[sgIDs_ind] lt critT, count_cold)
    endif else begin
      wCold = where(maxTemps[sgIDs_ind] lt critLogTemp, count_cold)
    endelse
      
    ; either "cold fraction" or "never heated above Tvir fraction"
    if (count_cold ne 0) then begin
      coldfrac[i] = float(count_cold) / counts[i]
    endif
  endfor
  
  ; subhalo masses (total baryonic+DM in log(msun))
  shmass = alog10( (units.UnitMass_in_g / units.Msun_in_g) * shmass )  
  
  ; save
  save,shmass,coldfrac,valSGids,counts,filename=saveFilename
  
  r = {shmass:shmass,coldfrac:coldfrac,valSGids:valSGids,counts:counts}
  return,r

end

; calcRadialPosAndVel(): calculate radial position distribution (with respect to r_200 of each parent halo)
;                        and decomposed radial/tangential velocity vectors of individual gas particles
;                        prior to accretion

function calcRadialPosAndVel, res=res, redshift=redshift, halos=halos, bSH=bSH, DM=DM

  ; config
  gadgetPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'
  workingPath  = '/n/home07/dnelson/coldflows/thermhist.deriv/'

  targetSnap = redshiftToSnapNum(redshift)
  endingSnap = targetSnap - 1 ; immediately preceeding
  
  if keyword_set(DM) then partName = 'dm'
  if not keyword_set(DM) then partName = 'gas'
  
  ; set saveFilename
  saveFilename = workingPath + 'radposvel.'+str(res)+'.'+str(targetSnap)+'-'+str(endingSnap)+'.sav'
  
  if (keyword_set(halos)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.halos.sav'
  if (keyword_set(bSH)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.bSH.sav'
  if (keyword_set(DM)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.DM.sav'    
    
  if (file_test(saveFilename)) then begin
    restore,saveFilename
    r = {radPos:radPos,radVel:radVel,tanVel:tanVel,parentIDs:parentIDs}
    return,r
  endif  
  
  ; load subhalo group catalogs
  sgTarget = loadSubhaloGroups(gadgetPath,targetSnap) 
  
  ; get list of smoothly accreted gas/dm particle ids
  sgIDs_Acc = findAccretedGas(res=res,bSH=bSH,targetSnap=targetSnap,sgTarget=sgTarget,halos=halos,DM=DM)
  
  ; arrays
  radPos = list()
  radVel = list()
  tanVel = list()
  
  radialVectors  = fltarr(3,n_elements(sgIDs_Acc))
  
  ; find parent subhalo IDs 
  parentIDs = findParentIDs(sgTarget,sgIDs_Acc)
  
  ; load particle ids
  ids = loadSnapshotSubset(gadgetPath,snapNum=endingSnap,partType=partName,field='ids')
  
  match,ids,sgIDs_Acc,ids_ind,sg_ind,count=count
  ids = !NULL
  
  if (count ne n_elements(sgIDs_Acc)) then begin
    print,'ERROR2'
    return,0
  endif
  
  ; load particle positions
  pos = loadSnapshotSubset(gadgetPath,snapNum=endingSnap,partType=partName,field='pos')  

  radialVectors = pos[*,ids_ind] - sgTarget.subgroupPos[*,parentIDs]
  pos = !NULL
  sgTarget = !NULL
  
  ; periodic B.C. for distances (effectively component by component)
  correctPeriodicDistVecs, radialVectors
    
  ; load particle velocities
  vel = loadSnapshotSubset(gadgetPath,snapNum=endingSnap,partType=partName,field='vel')  
  vel = vel[*,ids_ind]

  ; find unique parent IDs
  uniqParentIDs = parentIDs[uniq(parentIDs,sort(parentIDs))]
  
  print,'Processing [' + str(n_elements(uniqParentIDs)) + '] unique parent subhalos (of '+$
         str(n_elements(parentIDs))+' total):'

  ; loop over each unique parent ID
  foreach parentID,uniqParentIDs,k do begin
  
    if (k mod 100 eq 0) then print,' '+str(k)
    
      ; find all accreted gas particles of this parent subhalo
      w = where(parentIDs eq parentID,count)
 
      if (count gt 1) then begin
        bufRadPos = fltarr( ulong64(count) )
        bufRadVel = fltarr( ulong64(count) )
        bufTanVel = fltarr( ulong64(count) )
        ind = 0ULL
        
        ; save properties of each particle into buffer
        for i=0,count-1 do begin
          radNorm = sqrt(radialVectors[*,w[i]] ## transpose(radialVectors[*,w[i]]))
          
          radVelNorm = vel[*,w[i]] ## transpose(-radialVectors[*,w[i]]) / radNorm ; flipped sign for radVec
          radVelNorm = radVelNorm[0]
          tanVelVec  = vel[*,w[i]] - radVelNorm * radialVectors[*,w[i]] ; flipped sign for radVec
          tanVelNorm = sqrt(tanVelVec ## transpose(tanVelVec))
          stop ;TODO FINISH
          bufRadPos[ind] = radNorm
          bufRadVel[ind] = radVelNorm
          bufTanVel[ind] = tanVelNorm
          ind++
        endfor
        
        ; save buffer once per parentID
        radPos->add, bufRadPos
        radVel->add, bufRadVel
        tanVel->add, buftanVel
        parentIDs->add, parentID
      endif
  
  endforeach

  ; save and return  
  save,radPos,radVel,tanVel,parentIDs,filename=saveFilename
  
  print,'Saved: '+saveFilename

  r = {radPos:radPos,radVel:radVel,tanVel:tanVel,parentIDs:parentIDs}
  return,r  
  
end

; calcNormScalarProd(): calculate normalized scalar products of velocity vectors of individual
;                       gas particles prior to accretion, towards the subhalo center
;
; bSH=1 : include background subhalos
; halos=1 : find accretion onto main halos (do not set bSH,zWidth,redshiftsCut)
; DM=1 : calculate for dark matter particles instead of gas

function calcNormScalarProd, res=res, redshift=redshift, bSH=bSH, halos=halos, DM=DM

  ; config
  gadgetPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'
  workingPath  = '/n/home07/dnelson/coldflows/thermhist.deriv/'
  
  critLogTemp = 5.5

  targetSnap = redshiftToSnapNum(redshift)
  endingSnap = targetSnap - 1 ; immediately preceeding

  if keyword_set(DM) then partName = 'dm'
  if not keyword_set(DM) then partName = 'gas'

  ; set saveFilename
  saveFilename = workingPath + 'normscalar.'+str(res)+'.'+str(targetSnap)+'-'+str(endingSnap)+'.sav'
  
  if (keyword_set(halos)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.halos.sav'  
  if (keyword_set(bSH)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.bSH.sav'
  if (keyword_set(DM)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.DM.sav'
    
  if (file_test(saveFilename)) then begin
    restore,saveFilename
    r = {normProd_Hot:normProd_Hot,normProd_Cold:normProd_Cold,$
         parentIDs_Hot:parentIDs_Hot,parentIDs_Cold:parentIDs_Cold}
    return,r
  endif

  ; load subhalo group catalogs
  sgTarget = loadSubhaloGroups(gadgetPath,targetSnap) 
  
  ; get list of smoothly accreted gas/dm particle ids
  sgIDs_Acc = findAccretedGas(res=res,bSH=bSH,targetSnap=targetSnap,sgTarget=sgTarget,halos=halos,DM=DM)

  ; get maximum temperatures
  if not keyword_set(DM) then begin
    maxTemps = maxTemperatures(res=res,zMin=redshift)
    maxTemps = maxTemps[sgIDs_Acc]
  endif

  ; arrays
  normProd_Cold = list() ; or for dm
  normProd_Hot  = list() ; empty for dm
  parentIDs_Cold = list() ; or dm
  parentIDs_Hot  = list() ; empty for dm
  
  radialVectors  = fltarr(3,n_elements(sgIDs_Acc))
  
  ; find parent subhalo IDs 
  parentIDs = findParentIDs(sgTarget,sgIDs_Acc)
  
  ; load particle ids
  ids = loadSnapshotSubset(gadgetPath,snapNum=endingSnap,partType=partName,field='ids')
  
  match,ids,sgIDs_Acc,ids_ind,sg_ind,count=count
  ids = !NULL
  
  if (count ne n_elements(sgIDs_Acc)) then begin
    print,'ERROR2'
    return,0
  endif
  
  ; load particle positions
  pos = loadSnapshotSubset(gadgetPath,snapNum=endingSnap,partType=partName,field='pos')  

  radialVectors = pos[*,ids_ind] - sgTarget.subgroupPos[*,parentIDs]
  pos = !NULL
  sgTarget = !NULL
  
  ; periodic B.C. for distances (effectively component by component)
  boxSize = 20000.0 ;kpc
  w = where(radialVectors gt boxSize*0.5,count)
  if (count ne 0) then $
    radialVectors[w] = radialVectors[w] - boxSize
    
  w = where(radialVectors lt -boxSize*0.5,count)
  if (count ne 0) then $
    radialVectors[w] = boxSize + radialVectors[w]
    
  ;debug
  ;uniqind = n_elements(ids_ind[uniq(ids_ind,sort(ids_ind))])
  ;if (uniqind ne n_elements(ids_ind)) then begin
  ;  print,'ERROR3'
  ;  return,0
  ;endif

  ; find unique parent IDs
  uniqParentIDs = parentIDs[uniq(parentIDs,sort(parentIDs))]
  
  print,'Processing [' + str(n_elements(uniqParentIDs)) + '] unique parent subhalos (of '+$
         str(n_elements(parentIDs))+' total):'

  ; loop over each unique parent ID
  foreach parentID,uniqParentIDs,k do begin
  
    if (k mod 100 eq 0) then print,' '+str(k)
    
    if not keyword_set(DM) then begin
      ; find all COLD accreted gas particles of this parent subhalo
      w = where(parentIDs eq parentID and maxTemps lt critLogTemp,count)
  
      if (count gt 5e4) then begin
        print,'EXITING. Not a good idea, gas count='+str(count)+' and pairwise size would be ='+$
              str(ulong64(count)*(count-1)/2)
        return,0
      endif
  
      if (count gt 1) then begin
        buf = fltarr( ulong64(count)*(count-1)/2 ) ; long overflow on count
        ind = 0ULL
        
        ; COLD: pairwise (double loop) compute: normalized scalar product and save
        for i=0,count-1 do begin
          iip = sqrt(radialVectors[*,w[i]] ## transpose(radialVectors[*,w[i]]))
          for j=i+1,count-1 do begin
            dotp = radialVectors[*,w[i]] ## transpose(radialVectors[*,w[j]])
            norm = iip * sqrt(radialVectors[*,w[j]] ## transpose(radialVectors[*,w[j]]))
            buf[ind] = dotp/norm
            ind++
          endfor
        endfor
        
        normProd_Cold->add, buf
        parentIDs_Cold->add, parentID
      endif
  
      ; find all HOT accreted gas particles of this parent subhalo
      w = where(parentIDs eq parentID and maxTemps ge critLogTemp,countH)
      
      if (countH gt 5e4) then begin
        print,'EXITING. Not a good idea, countH='+str(count)+' and pairwise size would be ='+$
              str(ulong64(countH)*(countH-1)/2)
        return,0
      endif      
      
      if (countH gt 1) then begin
        buf = fltarr( ulong64(countH)*(countH-1)/2 )
        ind = 0ULL ;ugh int overflow from ~count^2
        
        ; HOT: pairwise product
        for i=0,countH-1 do begin
          iip = sqrt(radialVectors[*,w[i]] ## transpose(radialVectors[*,w[i]]))
          for j=i+1,countH-1 do begin
            dotp = radialVectors[*,w[i]] ## transpose(radialVectors[*,w[j]])
            norm = iip * sqrt(radialVectors[*,w[j]] ## transpose(radialVectors[*,w[j]]))
            buf[ind] = dotp/norm
            ind++
          endfor
        endfor
        normProd_Hot->add, buf
        parentIDs_Hot->add, parentID
      endif
    endif else begin ;DM
      w = where(parentIDs eq parentID,count)
      
      if (count gt 5e4) then begin
        print,'EXITING. Not a good idea, DM count='+str(count)+' and pairwise size would be ='+$
              str(ulong64(count)*(count-1)/2)
        return,0
      endif
      
      if (count gt 1) then begin
        buf = fltarr( ulong64(count)*(count-1)/2 ) ; long overflow on count
        ind = 0ULL
        
        ; COLD: pairwise (double loop) compute: normalized scalar product and save
        for i=0,count-1 do begin
          iip = sqrt(radialVectors[*,w[i]] ## transpose(radialVectors[*,w[i]]))
          for j=i+1,count-1 do begin
            dotp = radialVectors[*,w[i]] ## transpose(radialVectors[*,w[j]])
            norm = iip * sqrt(radialVectors[*,w[j]] ## transpose(radialVectors[*,w[j]]))
            buf[ind] = dotp/norm
            ind++
          endfor
        endfor
        
        normProd_Cold->add, buf
        parentIDs_Cold->add, parentID
      endif
      
    endelse
    
  endforeach
  
  save,normProd_Hot,normProd_Cold,parentIDs_Cold,parentIDs_Hot,filename=saveFilename
  
  print,'Saved: '+saveFilename

  r = {normProd_Hot:normProd_Hot,normProd_Cold:normProd_Cold,$
       parentIDs_Hot:parentIDs_Hot,parentIDs_Cold:parentIDs_Cold}
  return,r
  
end

; selectColdFilHalos(): select halos or subhalos by making cuts in the normalized scalar dot product
;                       (prefer accreting gas has preferential direction) and the thermal history
;                       of smoothly accreted gas (prefer cold mode dominated)

function selectColdFilHalos, res=res, halos=halos, redshift=redshift

  if (not keyword_set(res) or (halos ne 0 and halos ne 1)) then begin
    print,'Error: selectColdFilHalos: Bad inputs.'
    return,0
  endif

  ; config
  gadgetPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'
  workingPath  = '/n/home07/dnelson/coldflows/thermhist.deriv/'
  
  critLogTemp  = 5.5
  normProdMin  = 0.8 ; minimum mean normScalarProd (cold mode) to include a halo
  coldFracMin  = 0.5 ; minimum fraction of smoothly accreted gas that was in the cold mode
  minNumPart   = 100 ; minimum number of accreting gas particles
  
  bSH = 0
  
  targetSnap    = redshiftToSnapNum(redshift)
  smoothCutSnap = !NULL

  ; set saveFilename
  saveFilename = workingPath + 'selcoldfil.'+str(res)+'.'+str(targetSnap)+'-all.sav'
  
  if (keyword_set(halos)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.halos.sav'  
  if (keyword_set(bSH)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.bSH.sav'

  if (file_test(saveFilename)) then begin
    restore,saveFilename
    return,valSGids
  endif

  ; load catalog
  sgTarget = loadSubhaloGroups(gadgetPath,targetSnap)
  valSGids = getPrimarySubhaloList(sgTarget,halos=halos)

  print,'Searching ['+str(n_elements(valSGids))+'] halos for candidates:'

  ; load norm scalar dot product
  nsp = calcNormScalarProd(res=res, redshift=redshift, bSH=bSH, halos=halos)
  
  normProd_Cold  = nsp.normProd_Cold
  parentIDs_Cold = nsp.parentIDs_Cold
  nsp = !NULL
  
  npInclude = []
  
  ; save parentIDs meeting normProd cut  
  for i=0UL,n_elements(normProd_Cold)-1 do begin
    if (mean(normProd_Cold[i]) ge normProdMin and $
        n_elements(normProd_Cold[i]) ge minNumPart) then begin
      npInclude = [npInclude,parentIDs_Cold[i]]
      print,median(normProd_Cold[i]),mean(normProd_Cold[i]),n_elements(normProd_Cold[i])
    endif
  endfor
  
  ; keep only intersection
  match,valSGids,npInclude,val_ind,np_ind,count=count

  if (count eq 0) then begin
    print,'Out of candidates after normProd cut, exiting.'
    return,0
  endif
  
  print,' After normProd cut have ['+str(count)+'] candidates (lost '+$
        str(n_elements(valSGids)-count)+').'
        
  valSGids = valSGids[val_ind]
  
  ; load coldmode smoothly accreted gas fraction
  print,' Skipping coldFrac cut.'
  if 0 then begin
    cf = calcColdFracVsSubhaloMass(res=res,halos=halos,bSH=bSH,$
                                   targetSnap=targetSnap,smoothCutSnap=smoothCutSnap)
  
    cfSGind = where(cf.coldfrac ge coldFracMin,count)
    cfSGids = cf.valSGids[cfSGind]
  
    if (count eq 0) then begin
      print,'Out of candidates after colddFrac cut, exiting.'
      return,0
    endif
    
    ; keep only intersection
    match,valSGids,cfSGids,val_ind,cf_ind,count=countCF
    
    print,' After coldFrac cut have ['+str(countCF)+'] candidates (lost '+$
          str(n_elements(valSGids)-countCF)+').'
    
    if (countCF ne 0) then $
      valSGids = valSGids[val_ind]
    if (countCF eq 0) then $
      valSGids = -1
  endif
    
  ; save and return
  save,valSGids,filename=saveFilename
  print,'Saved: '+saveFilename
  
  return,valSGids

end

; selectFilament(): beginning of an algorithm to select a filamentary structure incident on a 
;                   specified subhalo [x,y,z] position based on a column-like spatial structure
;                   for all properties falling inside selection, plot various quantities

pro selectFilament

  units = getUnits()

  ; config
  targetSnap  = 189
  subhaloPath = '/n/hernquistfs1/mvogelsberger/ComparisonProject/512_20Mpc/Arepo_ENERGY/output/'
  workingPath = '/n/home07/dnelson/coldflows/vis/cutout_512A/'
  
  ; filament description
  x1 = [1023.2,7468.8,16090.0] ;start
  x2 = [1110.0,7540.0,16130.0] ;end  
  
  filRad   = 10.0      ; kpc
  filBound = [0.0,1.0] ; t
  
  ; 200kpc and 400kpc bboxes
  bbox    = [1023.20,1223.20,7468.80,7668.80,16044.2,16244.2]/1e5 ;Mpc
  bboxBig = [923.2,1323.2,7368.8,7768.8,15944.2,16344.2]/1e5 ;Mpc
  
  ; load subhalo cat and select in bbox
  sh = loadSubhaloGroups(subhaloPath,targetSnap,/verbose)
  
  wGR = where(sh.groupPos[0,*] gt bboxBig[0]*1e5 and sh.groupPos[0,*] lt bboxBig[1]*1e5 and $
              sh.groupPos[1,*] gt bboxBig[2]*1e5 and sh.groupPos[1,*] lt bboxBig[3]*1e5 and $
              sh.groupPos[2,*] gt bboxBig[4]*1e5 and sh.groupPos[2,*] lt bboxBig[5]*1e5, countGR)
  wSH = where(sh.subgroupPos[0,*] gt bboxBig[0]*1e5 and sh.subgroupPos[0,*] lt bboxBig[1]*1e5 and $
              sh.subgroupPos[1,*] gt bboxBig[2]*1e5 and sh.subgroupPos[1,*] lt bboxBig[3]*1e5 and $
              sh.subgroupPos[2,*] gt bboxBig[4]*1e5 and sh.subgroupPos[2,*] lt bboxBig[5]*1e5, countSH)
  print,'Found ['+str(countGR)+'] groups and ['+str(countSH)+'] subgroups in big bbox.'

  ; temp: 3d bbox pts for plotting
  bboxPts = [[bbox[0],bbox[2],bbox[4]],$
             [bbox[0],bbox[3],bbox[4]],$
             [bbox[1],bbox[3],bbox[4]],$
             [bbox[1],bbox[2],bbox[4]],$
             [bbox[0],bbox[2],bbox[4]],$
             [bbox[0],bbox[2],bbox[5]],$
             [bbox[0],bbox[3],bbox[5]],$
             [bbox[1],bbox[3],bbox[5]],$
             [bbox[1],bbox[2],bbox[5]],$
             [bbox[0],bbox[2],bbox[5]]]
  
  ; load particle positions
  x0_gas   = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='gas',field='pos',/verbose)
  x0_stars = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='stars',field='pos',/verbose) 
  x0_dm    = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='dm',field='pos',/verbose)
  
  ; parametric solution and minimum distance to line (gas)  
  n21  = (x2[0]-x1[0])^2.0 + (x2[1]-x1[1])^2.0 + (x2[2]-x1[2])^2.0 ;n21=transpose(x1-x2)#(x2-x1)
  n10  = reform( (x1[0]-x0_gas[0,*])^2.0 + (x1[1]-x0_gas[1,*])^2.0 + (x1[2]-x0_gas[2,*])^2.0 )
  dotp = reform( (x1[0]-x0_gas[0,*])*(x2[0]-x1[0]) + (x1[1]-x0_gas[1,*])*(x2[1]-x1[1]) + $
                 (x1[2]-x0_gas[2,*])*(x2[2]-x1[2]) )

  t_gas = -1.0 * dotp / n21
  d2    = ( n10 * n21 - dotp^2.0 ) / n21
  d_gas = sqrt(d2)
  
  ; (stars)
  n21  = (x2[0]-x1[0])^2.0 + (x2[1]-x1[1])^2.0 + (x2[2]-x1[2])^2.0 ;n21=transpose(x1-x2)#(x2-x1)
  n10  = reform( (x1[0]-x0_stars[0,*])^2.0 + (x1[1]-x0_stars[1,*])^2.0 + (x1[2]-x0_stars[2,*])^2.0 )
  dotp = reform( (x1[0]-x0_stars[0,*])*(x2[0]-x1[0]) + (x1[1]-x0_stars[1,*])*(x2[1]-x1[1]) + $
                 (x1[2]-x0_stars[2,*])*(x2[2]-x1[2]) )

  t_stars = -1.0 * dotp / n21
  d2      = ( n10 * n21 - dotp^2.0 ) / n21
  d_stars = sqrt(d2)
  
  ; (dm)
  n21  = (x2[0]-x1[0])^2.0 + (x2[1]-x1[1])^2.0 + (x2[2]-x1[2])^2.0 ;n21=transpose(x1-x2)#(x2-x1)
  n10  = reform( (x1[0]-x0_dm[0,*])^2.0 + (x1[1]-x0_dm[1,*])^2.0 + (x1[2]-x0_dm[2,*])^2.0 )
  dotp = reform( (x1[0]-x0_dm[0,*])*(x2[0]-x1[0]) + (x1[1]-x0_dm[1,*])*(x2[1]-x1[1]) + $
                 (x1[2]-x0_dm[2,*])*(x2[2]-x1[2]) )

  t_dm = -1.0 * dotp / n21
  d2   = ( n10 * n21 - dotp^2.0 ) / n21
  d_dm = sqrt(d2)
  
  ; make selection
  wGas   = where(t_gas   gt filBound[0] and d_gas   lt filRad and t_gas   lt filBound[1], countGas)
  wStars = where(t_stars gt filBound[0] and d_stars lt filRad and t_stars lt filBound[1], countStars)
  wDM    = where(t_dm    gt filBound[0] and d_dm    lt filRad and t_dm    lt filBound[1], countDM)
  
  print,'Filament selection: '+str(countGas)+ ' gas, '+str(countStars)+' stars, '+str(countDM)+' DM.'
  
  ; plot dist from line vs t
  start_PS,workingPath+'dist_vs_t.eps',xs=7,ys=10
    ; gas
    fsc_plot,[0],[0],/nodata,xrange=[0,40],yrange=[-1.4,2.0],$
             xtitle="Distance from Line",ytitle="Parameter t (Dist along Line)",/xs,/ys,$
             title="Gas",position=[0.1,0.1,0.3,0.9],charsize=!p.charsize-1
    fsc_plot,d_gas,t_gas,psym=3,symsize=1.0,/overplot
    fsc_plot,d_gas[wGas],t_gas[wGas],psym=3,symsize=2.0,color=fsc_color('green'),/overplot

    ; stars
    fsc_plot,[0],[0],/nodata,xrange=[0,40],yrange=[-1.4,2.0],$
             xtitle="Distance from Line",ytitle="Parameter t (Dist along Line)",/xs,/ys,$
             title="Stars",position=[0.4,0.1,0.6,0.9],charsize=!p.charsize-1,/noerase
    fsc_plot,d_stars,t_stars,psym=3,symsize=1.0,/overplot
    fsc_plot,d_stars[wStars],t_stars[wStars],psym=3,symsize=2.0,color=fsc_color('green'),/overplot
    
    ; dm
    fsc_plot,[0],[0],/nodata,xrange=[0,40],yrange=[-1.4,2.0],$
             xtitle="Distance from Line",ytitle="Parameter t (Dist along Line)",/xs,/ys,$
             title="DM",position=[0.7,0.1,0.9,0.9],charsize=!p.charsize-1,/noerase
    fsc_plot,d_dm,t_dm,psym=3,symsize=1.0,/overplot
    fsc_plot,d_dm[wDM],t_dm[wDM],psym=3,symsize=2.0,color=fsc_color('green'),/overplot
  end_PS
  
  ; plot selection (gas)
  start_PS, workingPath+'gas_select_xyz.eps', xs=8.0, ys=8.0
    !p.multi = [0,2,2]
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[0],bboxBig[1]],yrange=[bboxBig[2],bboxBig[3]],$
               xtitle="x [Mpc]",ytitle="y [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_gas[0,*]/1e5,x0_gas[1,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_gas[0,wGas]/1e5,x0_gas[1,wGas]/1e5,psym=3,symsize=3.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[0],x2[0]]/1e5,[x1[1],x2[1]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      ;for i=0,countGR-1 do begin
      ;  print,sh.grouppos[0,wGR[i]]/1e5,sh.groupPos[1,wGR[i]]/1e5,sh.group_R_Mean200[wGR[i]]/1e5
      ;  tvcircle,sh.group_r_mean200[wGR[i]]/1e5,sh.grouppos[0,wGR[i]]/1e5,sh.grouppos[1,wGR[i]]/1e5,$
      ;           color=fsc_color('blue'),/data,thick=!p.thick-1.0
      ;endfor
      for i=0,countSH-1 do begin
        ;print,sh.subgroupMass[wSH[i]]/5e4,sh.subgroupPos[0,wSH[i]]/1e5,sh.subgroupPos[1,wSH[i]]/1e5
        tvcircle,sh.subgroupMass[wSH[i]]/5e4,sh.subgroupPos[0,wSH[i]]/1e5,sh.subgroupPos[1,wSH[i]]/1e5,$
                 color=fsc_color('purple'),/data,thick=!p.thick-1.0
      endfor
      
      fsc_plot,[bbox[0],bbox[1],bbox[1],bbox[0],bbox[0]],$
               [bbox[2],bbox[2],bbox[3],bbox[3],bbox[2]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
      
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[2],bboxBig[3]],yrange=[bboxBig[4],bboxBig[5]],$
               xtitle="y [Mpc]",ytitle="z [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_gas[1,*]/1e5,x0_gas[2,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_gas[1,wGas]/1e5,x0_gas[2,wGas]/1e5,psym=3,symsize=1.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[1],x2[1]]/1e5,[x1[2],x2[2]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      for i=0,countSH-1 do begin
        tvcircle,sh.subgroupMass[wSH[i]]/5e4,sh.subgroupPos[1,wSH[i]]/1e5,sh.subgroupPos[2,wSH[i]]/1e5,$
                 color=fsc_color('purple'),/data,thick=!p.thick-1.0
      endfor
      
      fsc_plot,[bbox[2],bbox[3],bbox[3],bbox[2],bbox[2]],$
               [bbox[4],bbox[4],bbox[5],bbox[5],bbox[4]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
      
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[0],bboxBig[1]],yrange=[bboxBig[4],bboxBig[5]],$
               xtitle="x [Mpc]",ytitle="z [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_gas[0,*]/1e5,x0_gas[2,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_gas[0,wGas]/1e5,x0_gas[2,wGas]/1e5,psym=3,symsize=1.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[0],x2[0]]/1e5,[x1[2],x2[2]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      for i=0,countSH-1 do begin
        tvcircle,sh.subgroupMass[wSH[i]]/5e4,sh.subgroupPos[0,wSH[i]]/1e5,sh.subgroupPos[2,wSH[i]]/1e5,$
                 color=fsc_color('purple'),/data,thick=!p.thick-1.0
      endfor
      
      fsc_plot,[bbox[0],bbox[1],bbox[1],bbox[0],bbox[0]],$
               [bbox[4],bbox[4],bbox[5],bbox[5],bbox[4]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
               
    !p.multi = 0
  end_PS
  
  ; (stars)
  start_PS, workingPath+'stars_select_xyz.eps', xs=8.0, ys=8.0
    !p.multi = [0,2,2]
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[0],bboxBig[1]],yrange=[bboxBig[2],bboxBig[3]],$
               xtitle="x [Mpc]",ytitle="y [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_stars[0,*]/1e5,x0_stars[1,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_stars[0,wStars]/1e5,x0_stars[1,wStars]/1e5,psym=3,symsize=3.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[0],x2[0]]/1e5,[x1[1],x2[1]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      fsc_plot,[bbox[0],bbox[1],bbox[1],bbox[0],bbox[0]],$
               [bbox[2],bbox[2],bbox[3],bbox[3],bbox[2]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
      
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[2],bboxBig[3]],yrange=[bboxBig[4],bboxBig[5]],$
               xtitle="y [Mpc]",ytitle="z [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_stars[1,*]/1e5,x0_stars[2,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_stars[1,wStars]/1e5,x0_stars[2,wStars]/1e5,psym=3,symsize=1.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[1],x2[1]]/1e5,[x1[2],x2[2]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      fsc_plot,[bbox[2],bbox[3],bbox[3],bbox[2],bbox[2]],$
               [bbox[4],bbox[4],bbox[5],bbox[5],bbox[4]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
      
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[0],bboxBig[1]],yrange=[bboxBig[4],bboxBig[5]],$
               xtitle="x [Mpc]",ytitle="z [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_stars[0,*]/1e5,x0_stars[2,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_stars[0,wStars]/1e5,x0_stars[2,wStars]/1e5,psym=3,symsize=1.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[0],x2[0]]/1e5,[x1[2],x2[2]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      fsc_plot,[bbox[0],bbox[1],bbox[1],bbox[0],bbox[0]],$
               [bbox[4],bbox[4],bbox[5],bbox[5],bbox[4]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
               
    !p.multi = 0
  end_PS
  
  ; (dm)
  start_PS, workingPath+'dm_select_xyz.eps', xs=8.0, ys=8.0
    !p.multi = [0,2,2]
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[0],bboxBig[1]],yrange=[bboxBig[2],bboxBig[3]],$
               xtitle="x [Mpc]",ytitle="y [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_dm[0,*]/1e5,x0_dm[1,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_dm[0,wDM]/1e5,x0_dm[1,wDM]/1e5,psym=3,symsize=3.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[0],x2[0]]/1e5,[x1[1],x2[1]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      fsc_plot,[bbox[0],bbox[1],bbox[1],bbox[0],bbox[0]],$
               [bbox[2],bbox[2],bbox[3],bbox[3],bbox[2]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
      
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[2],bboxBig[3]],yrange=[bboxBig[4],bboxBig[5]],$
               xtitle="y [Mpc]",ytitle="z [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_dm[1,*]/1e5,x0_dm[2,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_dm[1,wDM]/1e5,x0_dm[2,wDM]/1e5,psym=3,symsize=1.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[1],x2[1]]/1e5,[x1[2],x2[2]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      fsc_plot,[bbox[2],bbox[3],bbox[3],bbox[2],bbox[2]],$
               [bbox[4],bbox[4],bbox[5],bbox[5],bbox[4]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
      
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[0],bboxBig[1]],yrange=[bboxBig[4],bboxBig[5]],$
               xtitle="x [Mpc]",ytitle="z [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_dm[0,*]/1e5,x0_dm[2,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_dm[0,wDM]/1e5,x0_dm[2,wDM]/1e5,psym=3,symsize=1.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[0],x2[0]]/1e5,[x1[2],x2[2]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      fsc_plot,[bbox[0],bbox[1],bbox[1],bbox[0],bbox[0]],$
               [bbox[4],bbox[4],bbox[5],bbox[5],bbox[4]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
               
    !p.multi = 0
  end_PS
  
  ; load gas properties
  dens = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='gas',field='density')
  u     = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='gas',field='u')
  nelec = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='gas',field='ne')
  mass = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='gas',field='masses')
  
  dens = dens[wGas]
  u     = u[wGas]
  nelec = nelec[wGas]
  mass = mass[wGas]

  temp = convertUtoTemp(u,nelec)
  entropy = calcEntropy(u,dens)
  
  ; plot gas properties
  start_PS, workingPath+'gas_properties.eps'
    !p.multi = [0,2,2]
       plothist,alog10(rhoRatioToCrit(dens)),/auto,$
                xtitle="log ("+textoidl("\rho / \rho_{crit}")+")",ytitle="N",$
                title="Gas Density",charsize=!p.charsize-1
       
       plothist,alog10(temp),/auto,xtitle="log (T [K])",ytitle="N",$
                title="Gas Temperature",charsize=!p.charsize-1
                
       h2rt = rhoTHisto(dens,temp,mass=mass,/plot)
       
       plothist,alog10(entropy),/auto,xtitle="log (Entropy) [Code]",ytitle="N",$
                title="Gas Entropy",charsize=!p.charsize-1
    !p.multi = 0
  end_PS
  
  ; velocities - find distribution of cos(theta) between line and individual vel vectors
  vels = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='gas',field='vel')
  vels = vels[*,wGas]

  line = x2-x1

  dotp  = reform( line[0]*vels[0,*] + line[1]*vels[1,*] + line[2]*vels[2,*] )
  norm1 = line[0]*line[0] + line[1]*line[1] + line[2]*line[2]
  norm2 = reform( (vels[0,*]*vels[0,*] + vels[1,*]*vels[1,*] + vels[2,*]*vels[2,*]) )

  cos_theta = dotp / sqrt(norm1*norm2)
  theta_deg_gas = acos(cos_theta) * 180.0/!pi
  
  ; (dm)
  vels = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='dm',field='vel')
  vels = vels[*,wDM]

  line = x2-x1

  dotp  = reform( line[0]*vels[0,*] + line[1]*vels[1,*] + line[2]*vels[2,*] )
  norm1 = line[0]*line[0] + line[1]*line[1] + line[2]*line[2]
  norm2 = reform( (vels[0,*]*vels[0,*] + vels[1,*]*vels[1,*] + vels[2,*]*vels[2,*]) )

  cos_theta = dotp / sqrt(norm1*norm2)
  theta_deg_dm = acos(cos_theta) * 180.0/!pi

  start_PS, workingPath+'select_vel_dir.eps'
      plothist,theta_deg_gas,/auto,xtitle="Relative Angle between Part. Vel and Filament [deg]",ytitle="N",$
               xrange=[-2.0,180.0],yrange=[0,1.05],/peak,/ys,/xs
      plothist,theta_deg_dm,/auto,/overplot,/peak,color=fsc_color('orange')
      fsc_text,150.0,0.9,"Gas",color=fsc_color('black')
      fsc_text,150.0,0.84,"DM",color=fsc_color('orange')
  end_PS
  
  stop
end

@coldflowsPlot
