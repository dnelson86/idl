; timeScales.pro
; cooling times of halo gas vs. dynamical/hubble timescales
; dnelson nov.2012

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
  dens *= units.UnitDensity_in_cgs * units.HubbleParam * units.HubbleParam / (h.time)^3.0
  u *= units.UnitPressure_in_cgs / units.UnitDensity_in_cgs
  
  ; calculate cooling rates using primordial network
  r.cooltime = CalcCoolTime(u,dens,nelec,scalefac=h.time)
  
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
  
  ; cum mass plot
  start_PS, sP.plotPath + 'cumMass.'+sP.plotPrefix+'.'+str(sP.res) + '.' + str(sP.snap)+'.eps'
    cgPlot,galcatRadii.gmem_sec,massEnc,psym=3,xtitle="Radius [ckpc]",ytitle="Enclosed Mass [ Code ]"
  end_PS
  
  return, massEnc
  
end

; compareTimescalesHalo(): explore for a single halo or a mass range stacked

pro compareTimescalesHalo, redshift=redshift

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  if ~keyword_set(redshift) then message,'set redshift'
  sP = simParams(res=512,run='tracer',redshift=redshift)
  gc = loadGroupCat(sP=sP,/skipIDs)
  h  = loadSnapshotHeader(sP=sP)
  
  ; 1. single halo by haloID
  ;haloID = 304 ;z2.304 z2.301 z2.130 z2.64
  ;gcID = getMatchedIDs(sPa=sP,sPg=sP,haloID=haloID)
  ;gcIDList = [gcID.a]
  ;hTag = 'h'+str(haloID)+'_m='+string(codeMassToLogMsun(gc.subgroupMass[gcID.a]),format='(f4.1)')

  ; 2. single/several halo(s) closest to given mass(es)
  ;massTarget = [11.0] ; log Msun
  ;gcIDList = massTargetToHaloID(massTarget,sP=sP)
  ;hTag = 'mt'+str(n_elements(massTarget))
  
  ; 3. mass range stacked
  massRange = [9.98,10.0] ; log Msun
  gcIDList = where(codeMassToLogMsun(gc.subgroupMass) ge massRange[0] and $
                   codeMassToLogMsun(gc.subgroupMass) lt massRange[1])
  hTag = 'mbin.'+string(massRange[0],format='(f5.2)')+'-'+string(massRange[1],format='(f5.2)')

  ; get indices for gas in halo(s)
  galcat = galaxyCat(sP=sP)
  inds = galcatINDList(sP=sP, galcat=galcat, gcIDList=gcIDList)
  print,n_elements(gcIDList),n_elements(inds.gmem)

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
  
  meanDensEnc = 3*encMass / (4 * !pi * gasRadii^3.0) / (h.time)^3.0 ; code units (physical)
  
  dynTime = sqrt( 3*!pi / (32 * float(units.G) * meanDensEnc * units.HubbleParam) ) ; code units (Gyr)
  
  ; age of universe
  hubbleTime = snapNumToAgeFlat(sP=sP)
  
  ; cooling and dynamical timescales for the halo as a whole (Tvir,subgroup Mtot at rvir)
  grNr = gc.subgroupGrNr[gcIDList[0]]
  
  virTemp_halo = codeMassToVirTemp(gc.subgroupMass[gcIDList[0]],sP=sP,/log)
  meanDensEnc_halo = 3*gc.subgroupMassType[partTypeNum('gas'),gcIDList[0]] / (4*!pi*gc.group_r_crit200[grNr]^3.0) / (h.time)^3.0                
    
  meanDens_cgs = meanDensEnc_halo * units.UnitDensity_in_cgs * units.HubbleParam * units.HubbleParam
  coolTime_halo = CalcCoolTime(virTemp_halo,meanDens_cgs,1.0,scalefac=(h.time+1.0))
  coolTime_halo *= units.HubbleParam / units.UnitTime_in_s
  
  meanDensEnc_halo = 3*gc.subgroupMass[gcIDList[0]] / (4*!pi*gc.group_r_crit200[grNr]^3.0) / (h.time)^3.0 ; gas+DM+stars
  dynTime_halo = sqrt( 3*!pi / (32 * float(units.G) * meanDensEnc_halo * units.HubbleParam) )

  ; plots
  coolingRange = [0.01,30]
  tempRange = [4.0,7.0]
  densRange = 10.0^[-9.0,-6.0]
  radRange = [0.0,1.75]
  binsize = 0.1 / (sP.res/128)
  
  psym = 4
  if n_elements(inds.gmem) gt 2000 then psym = 3
  
  ; plot (1) - histograms
  start_PS, sP.plotPath + 'timescales_histo.' + hTag + '.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps'
    
    hist = histogram(alog10(coolTime),binsize=binsize,loc=loc)
    cgPlot,10.0^loc,hist/float(total(hist)),color=cgColor('red'),$
      xtitle="Timescale [Gyr]",ytitle="Halo Gas Mass Fraction",xrange=coolingRange,yrange=[0.0,0.7]/(sP.res/128),/xlog
    
    hist = histogram(alog10(dynTime[where(finite(dynTime))]),binsize=binsize,loc=loc)
    cgPlot,10.0^loc,hist/float(total(hist)),color=cgColor('blue'),/overplot
    
    cgPlot,[hubbleTime,hubbleTime],[0,1],color=cgColor('green'),/overplot
    
    legend,textoidl(['t_{cool}','t_{dyn}','t_{age}']),textcolors=['red','blue','green'],/top,/right,box=0
    
  end_PS
  
  ; plot (2) - cooling/dynamical scatter
  start_PS, sP.plotPath + 'timescales_scat.' + hTag + '.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps'
    cgPlot,coolTime,dynTime,psym=psym,xtitle="Cooling Time [Gyr]",ytitle="Dynamical Time [Gyr]",xrange=coolingRange
  end_PS

  ; plot (3) - cooling time vs current temperature scatter
  start_PS, sP.plotPath + 'timescales_cooltemp_scat.' + hTag + '.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps'
    cgPlot,coolTime,curTemp,psym=psym,xtitle="Cooling Time [Gyr]",ytitle="Gas Temperature [log K]",$
      xrange=coolingRange,yrange=tempRange
  end_PS
  
  ; plot (4) - cooling time vs current density scatter
  start_PS, sP.plotPath + 'timescales_cooldens_scat.' + hTag + '.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps'
    cgPlot,coolTime,curDens,psym=psym,xtitle="Cooling Time [Gyr]",ytitle="Density [Code]",$
      xrange=coolingRange,yrange=densRange,/ylog
  end_PS
  
  ; plot (5) - cooling/dynamical vs radius (temp colormap)
  start_PS, sP.plotPath + 'timescales_vsrad1.' + hTag + '.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps'
    cgPlot,[0],[0],/nodata,xtitle=textoidl("r / r_{vir}"),ytitle="Timescale [Gyr]",$
      xrange=radRange,yrange=coolingRange,/ylog,/xs,/ys,yminor=0,$
      title=hTag + ' (z='+string(sP.redshift,format='(f3.1)')+')',position=[0.15,0.13,0.85,0.9]
    
    ; individual gas elements: cooling time color mapped by temperature
    loadColorTable,'blue-red';,/reverse
    TVLCT, rr, gg, bb, /GET
    
    fieldMinMax = [5.0,6.5]
    colorinds = (curTemp-fieldMinMax[0])*235.0 / (fieldMinMax[1]-fieldMinMax[0]) ;0-235
    colorinds = fix(colorinds + 20.0) > 0 < 255 ;20-255
    colorinds = getColor24([[rr[colorinds]], [gg[colorinds]], [bb[colorinds]]])
    
    for i=0L,n_elements(coolTime)-1 do $
      oplot,[gasRadii[i]/gasRvir[i]],[coolTime[i]],psym=psym,color=colorinds[i]
    
    ;cgPlot,gasRadii/gasRvir,coolTime,psym=psym,color=cgColor('red'),/overplot
    
    ; individual gas elements
    cgPlot,gasRadii/gasRvir,dynTime,psym=psym,color=cgColor('blue'),/overplot
    cgPlot,[0.05,1.5],[hubbleTime,hubbleTime],color=cgColor('green'),/overplot
    
    ; mean halo
    cgPlot,[0.05,1.0],[coolTime_halo,coolTime_halo],color=cgColor('magenta'),line=0,/overplot
    cgPlot,[0.05,1.0],[dynTime_halo,dynTime_halo],color=cgColor('orange'),line=0,/overplot
     
    ; legend
    legend,textoidl(['t_{cool}','t_{dyn}','t_{age}','t_{cool,halo}','t_{dyn,halo}']),$
      textcolors=['red','blue','green','magenta','orange'],/bottom,/right,box=0
      
    ; colorbar
    colorbar,/right,/vertical,position=[0.88,0.13,0.93,0.9],range=fieldMinMax,title="log Temp"
  end_PS
  
  ; plot (6) - cooling/dynamical vs radius (dens colormap)
  start_PS, sP.plotPath + 'timescales_vsrad2.' + hTag + '.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps'
    cgPlot,[0],[0],/nodata,xtitle=textoidl("r / r_{vir}"),ytitle="Timescale [Gyr]",$
      xrange=radRange,yrange=coolingRange,/ylog,/xs,/ys,yminor=0,$
      title=hTag + ' (z='+string(sP.redshift,format='(f3.1)')+')',position=[0.15,0.13,0.85,0.9]
    
    ; individual gas elements: cooling time color mapped by temperature
    loadColorTable,'blue-red';,/reverse
    TVLCT, rr, gg, bb, /GET
    
    fieldMinMax = [-8.0,-6.0]
    colorinds = (alog10(curDens)-fieldMinMax[0])*235.0 / (fieldMinMax[1]-fieldMinMax[0]) ;0-235
    colorinds = fix(colorinds + 20.0) > 0 < 255 ;20-255
    colorinds = getColor24([[rr[colorinds]], [gg[colorinds]], [bb[colorinds]]])
    
    for i=0L,n_elements(coolTime)-1 do $
      oplot,[gasRadii[i]/gasRvir[i]],[coolTime[i]],psym=psym,color=colorinds[i]
    
    ;cgPlot,gasRadii/gasRvir,coolTime,psym=psym,color=cgColor('red'),/overplot
    
    ; individual gas elements
    cgPlot,gasRadii/gasRvir,dynTime,psym=psym,color=cgColor('blue'),/overplot
    cgPlot,[0.05,1.5],[hubbleTime,hubbleTime],color=cgColor('green'),/overplot
    
    ; mean halo
    cgPlot,[0.05,1.0],[coolTime_halo,coolTime_halo],color=cgColor('magenta'),line=0,/overplot
    cgPlot,[0.05,1.0],[dynTime_halo,dynTime_halo],color=cgColor('orange'),line=0,/overplot
     
    ; legend
    legend,textoidl(['t_{cool}','t_{dyn}','t_{age}','t_{cool,halo}','t_{dyn,halo}']),$
      textcolors=['red','blue','green','magenta','orange'],/bottom,/right,box=0
      
    ; colorbar
    colorbar,/right,/vertical,position=[0.88,0.13,0.92,0.9],range=fieldMinMax,title="log Dens"
  end_PS

  stop
  
end

; reesOstrikerFig1(): recreate Fig.1 of Rees & Ostriker (1977) to verify our cooling time calculations

pro reesOstrikerFig1

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  tempMinMax = [2.0,9.0] ; log(K)
  densMinMaxPlot = [-7,3] ; log(cm^-3)
  densMinMaxCalc = [-5,3] ; log(cm^-3)
  nTemps = 200
  nDens  = 200
  
  redshifts = [0.0,1.0,2.0]
  colors = ['red','blue','green']
  
  ; start plot
  start_PS,'fig1.eps'
    cgPlot,[0],[0],/nodata,xrange=10.0^densMinMaxPlot,yrange=10.0^tempMinMax,xtitle=textoidl("n [cm^{-3}]"),$
      ytitle="Temp [K]",/xlog,/ylog,/xs,/ys,xminor=1,yminor=1
  
  ; generate arrays
  temps = 10.0^(findgen(nTemps)/(nTemps-1) * (tempMinMax[1]-tempMinMax[0]) + tempMinMax[0]) ; K
  dens  = 10.0^(findgen(nDens)/(nDens-1) * (densMinMaxCalc[1]-densMinMaxCalc[0]) + densMinMaxCalc[0]) ; 1/cm^3
  dens_gcm3 = dens * units.mass_proton * units.HubbleParam * units.HubbleParam ; g/cm^3
  dens_code = float(dens_gcm3 / units.UnitDensity_in_cgs) ; code units (10^10 msun/kpc^3)
  nelec = fltarr(nDens) + 1.0

  ; constant jeans mass lines (M_jeans = 1e8 * (temps/1e4)^(1.5) * dens^(-0.5) ; Msun)
  M_jeans_targets = [1e12,1e11,1e9]
  temp_jeans0 = ( (M_jeans_targets[0]/1e8) / dens^(-0.5) ) ^ (2.0/3.0) * 1e4
  temp_jeans1 = ( (M_jeans_targets[1]/1e8) / dens^(-0.5) ) ^ (2.0/3.0) * 1e4
  temp_jeans2 = ( (M_jeans_targets[2]/1e8) / dens^(-0.5) ) ^ (2.0/3.0) * 1e4
  
  ; plot: constant jeans mass lines
  w = where(dens ge 10.0^(-5.5))
  cgPlot,dens[w],temp_jeans0[w],line=2,/overplot
  cgPlot,dens[w],temp_jeans1[w],line=2,/overplot
  cgPlot,dens[w],temp_jeans2[w],line=2,/overplot
  cgText,10^(-5.9),10^(5.2),textoidl('10^{12} M_{sun}'),/data,alignment=0.5
  cgText,10^(-5.5),10^(3.7),textoidl('10^{11} M_{sun}'),/data,alignment=0.5
  cgText,10^(-4.0),10^(2.5),textoidl('M_J=10^{9} M_{sun}'),/data,alignment=0.5
  
  ; ABC regions
  cgText,10^(-4.0),1e8,"A",/data,alignment=0.5
  cgText,10^(-0.6),1e8,"B",/data,alignment=0.5
  cgText,10^(1.0),10^(5.5),"C",/data,alignment=0.5
  
  foreach redshift,redshifts,m do begin
  
    ; age of universe
    scalefac = 1.0/(1+redshift)
    print,scalefac
  
    tage = redshiftToAgeFlat(redshift)

    ; tcool = tage and tcool = tdyn lines
    dens_tcool_tage = fltarr(nTemps)
    dens_tcool_tdyn = fltarr(nTemps)
  
    ; estimate the dynamical (free-fall) timescale for these densities
    tdyn_eval = sqrt( 3*!pi / (32 * float(units.G) * dens_code/scalefac^3.0) ) ; code units (Gyr)
  
    for i=0,nTemps-1 do begin
      ; for each temperature, calculate the tcool across all densities
      temps_loc = replicate(alog10(temps[i]),nDens)
    
      tcool_eval = calcCoolTime(temps_loc,dens_gcm3,nelec,scalefac=scalefac+1.0)
      tcool_eval *= float(units.HubbleParam / units.UnitTime_in_s) ; convert cgs -> code units (Gyr)
  
      ; interpolate the evaluated tcool to tage and save 
      dens_tcool_tage[i] = interpol(dens,tcool_eval-tage,0.0)

      ; interpolate the evaluate tcool to tdyn and save
      dens_tcool_tdyn[i] = interpol(dens,tcool_eval-tdyn_eval,0.0)
      
      ;if i eq 150 then stop
    endfor
      
    ; plot tcool=tage/tdyn lines
    w = where(temps ge 10.0^(5.1))
    cgPlot,dens_tcool_tage[w],temps[w],line=1,/overplot,color=cgColor(colors[m])
    cgPlot,dens_tcool_tdyn[w],temps[w],line=0,/overplot,color=cgColor(colors[m])
  
  endforeach ;redshifts

  legend,'z='+string(redshifts,format='(f3.1)'),textcolors=colors,box=0,/top,/left
  
  ; end plot
  end_PS
  
  ; for a fixed temperature, plot the cooling time vs density at different redshifts
  redshifts = [0.0,0.5,1.0,1.5,2.0]
  colors = ['red','blue','green','magenta','orange']
  log_temp = [7.0,8.0]
  
  start_PS,'fig1_b.eps'
    cgPlot,[0],[0],/nodata,xrange=10.0^densMinMaxPlot,yrange=[0.01,1e4],xtitle=textoidl("n [cm^{-3}]"),$
      ytitle="Cooling/Dynamical Time [Gyr]",/xlog,/ylog,/xs,/ys,xminor=1,yminor=1
      
  foreach redshift,redshifts,m do begin
    scalefac = 1.0/(1+redshift)
    
    ; temp 1
    temps_loc = replicate(log_temp[0],nDens)
    tcool_eval = calcCoolTime(temps_loc,dens_gcm3,nelec,scalefac=scalefac+1.0)
     tcool_eval *= float(units.HubbleParam / units.UnitTime_in_s)
     
    cgPlot,dens,tcool_eval,color=cgColor(colors[m]),/overplot
    
    ; temp 2
    temps_loc = replicate(log_temp[1],nDens)
    tcool_eval = calcCoolTime(temps_loc,dens_gcm3,nelec,scalefac=scalefac+1.0)
     tcool_eval *= float(units.HubbleParam / units.UnitTime_in_s)
     
    cgPlot,dens,tcool_eval,color=cgColor(colors[m]),line=2,/overplot
    
    ; dynamical
    tdyn_eval = sqrt( 3*!pi / (32 * float(units.G) * dens_code/scalefac^3.0) ) * float(units.HubbleParam)
    cgPlot,dens,tdyn_eval,color=cgColor(colors[m]),line=1,/overplot
  endforeach
  
  legend,'z='+string(redshifts,format='(f3.1)'),textcolors=colors,box=0,/top,/right
  
  end_PS
  
  stop

end