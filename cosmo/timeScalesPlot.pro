; timeScalesPlot.pro
; plotting: cooling times of halo gas vs. dynamical/hubble timescales
; dnelson dec.2012

; plotTSFracsVsHaloMass(): plot gas mass fractions for timescale ratios vs halo mass

pro plotTSFracsVsHaloMass

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  sP = simParams(res=512,run='tracer',redshift=2.0)
  sgSelect = 'pri'
  
  ts = timescaleFracsVsHaloMass(sP=sP,sgSelect=sgSelect)
  
  xrange = [10.0,12.0]
  yrange = [0.0,1.05]
  cInd   = 1 ; colorindex
  
  ; plot
  start_PS, sP.plotPath + 'tsfrac_tdyn_vsmass.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),ytitle="Gas Mass Fraction"
      
    j=1
    cgPlot,ts.gcMasses,ts.tsFracs.gmem_tcool_tdyn[j,*],psym=4,color=sP.colorsA[cInd],/overplot
      
    for j=0,n_elements(ts.tsRatioVals)-1 do $
      cgPlot,ts.logMassBinCen,ts.tsMedian.gmem_tcool_tdyn[j,*],color=cgColor('black'),line=j,/overplot
      
    ; legend
    strings = textoidl("t_{cool} / t_{dyn} < "+string(ts.tsRatioVals,format='(f3.1)'))
    legend,strings,linestyle=indgen(n_elements(ts.tsRatioVals)),box=0,linesize=0.25,$
      position=[11.25,0.97],charsize=!p.charsize-0.1
  end_PS
  
  start_PS, sP.plotPath + 'tsfrac_tage_vsmass.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),ytitle="Gas Mass Fraction"
      
    j=1
    cgPlot,ts.gcMasses,ts.tsFracs.gmem_tcool_tage[j,*],psym=4,color=sP.colorsA[cInd],/overplot
      
    for j=0,n_elements(ts.tsRatioVals)-1 do $
      cgPlot,ts.logMassBinCen,ts.tsMedian.gmem_tcool_tage[j,*],color=cgColor('black'),line=j,/overplot
      
    ; legend
    strings = textoidl("t_{cool} / t_{age} < "+string(ts.tsRatioVals,format='(f3.1)'))
    legend,strings,linestyle=indgen(n_elements(ts.tsRatioVals)),box=0,linesize=0.25,$
      position=[11.25,0.97],charsize=!p.charsize-0.1
  end_PS

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
  haloID = 130 ;z2.304 z2.301 z2.130 z2.64
  gcID = getMatchedIDs(sPa=sP,sPg=sP,haloID=haloID)
  gcIDList = [gcID.a]
  hTag = 'h'+str(haloID)+'_m='+string(codeMassToLogMsun(gc.subgroupMass[gcID.a]),format='(f4.1)')

  ; 2. single/several halo(s) closest to given mass(es)
  ;massTarget = [11.0] ; log Msun
  ;gcIDList = massTargetToHaloID(massTarget,sP=sP)
  ;hTag = 'mt'+str(n_elements(massTarget))
  
  ; 3. mass range stacked
  ;massRange = [10.6,10.7] ; log Msun
  ;gcIDList = where(codeMassToLogMsun(gc.subgroupMass) ge massRange[0] and $
  ;                 codeMassToLogMsun(gc.subgroupMass) lt massRange[1])
  ;hTag = 'mbin.'+string(massRange[0],format='(f5.2)')+'-'+string(massRange[1],format='(f5.2)')

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
  coolTime_halo = CalcCoolTime(virTemp_halo,meanDens_cgs,1.0,scalefac=(h.time+1.0))
  coolTime_halo *= units.HubbleParam / units.UnitTime_in_s
  
  meanDensEnc_halo = 3*gc.subgroupMass[gcIDList[0]] / (4*!pi*gc.group_r_crit200[grNr]^3.0) ; gas+DM+stars
  meanDensEnc_halo = codeDensToPhys( meanDensEnc_halo, scalefac=h.time )
  dynTime_halo = sqrt( 3*!pi / (32 * float(units.G) * meanDensEnc_halo) )

  ; alternative definitions for dynamical time
  dynTime_halo2 = gc.group_r_crit200[grNr] / sqrt(units.G * gc.subgroupMass[gcIDList[0]] / gc.group_r_crit200[grNr] )
  dynTime_halo2 *= units.hubbleParam * h.time
  
  dynTime2 = (gasRadii * units.HubbleParam * h.time) / sqrt(units.G*encMass/gasRadii/h.time)
  
  ; plots
  coolingRange = [0.05,20]
  tempRange = [4.0,7.0]
  densRange = 10.0^[-9.0,-5.5]
  radRange = [0.1,1.75]
  binsize = 0.1 / (sP.res/128)
  
  psym = 4
  if n_elements(inds.gmem) gt 10000 then psym = 3
  
  ; radial fitting
  radFitRange = [0.15,1.5]
  
  radCt = fitRadProfile(radii=gasRadii/gasRvir,vals=coolTime,range=radFitRange,radBins=15)
  radDt = fitRadProfile(radii=gasRadii/gasRvir,vals=dynTime,range=radFitRange,radBins=40)
  radDens = fitRadProfile(radii=gasRadii/gasRvir,vals=curDens,range=radFitRange,radBins=15)
  radTemp = fitRadProfile(radii=gasRadii/gasRvir,vals=curTemp,range=radFitRange,radBins=15)
  
  ; calculate hot halo gas mass
  ids_gmem = galcat.groupmemIDs[inds.gmem]
  
  ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
  match,ids,ids_gmem,ids_ind,ids_gmem_ind,count=countMatch
  if countMatch ne n_elements(inds.gmem) then message,'error'
  ids = !NULL
  ids_ind = ids_ind[sort(ids_gmem_ind)]
  
  masses = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
  masses = masses[ids_ind]
  mass_hot = total(masses) * units.HubbleParam
  
  ; profile fitting
  x = linspace(0.13,1.6,50)
  meanHaloDMMass  = mean(gc.subgroupMassType[partTypeNum('dm'),gcIDList]) * units.HubbleParam
  meanHaloGasMass = mean(gc.subgroupMassType[partTypeNum('gas'),gcIDList]) * units.HubbleParam
  
  ; SIS
  sis_dm  = sis_profile(x, mass=meanHaloDMMass, redshift=sP.redshift)
  sis_gas = sis_gas_profile(mass_hot=mass_hot, sis_dm=sis_dm, tables=tables)
  sis_fit = sis_gas_fit(rad_frac=gasRadii/gasRvir, dens=curDens, temp=curTemp, x=x)
  
  ; NFW
  nfw_dm  = nfw_profile(x, mass=meanHaloDMMass, redshift=sP.redshift)
  nfw_gas = nfw_gas_suto(mass_hot=mass_hot, nfw_dm=nfw_dm, tables=tables)
  nfw_fit = nfw_gas_fit(rad_frac=gasRadii/gasRvir, dens=curDens, temp=curTemp, x=x, nfw_dm=nfw_dm)
  
  ; convert to same units as simulation (comoving/h)
  sis_gas.rho_gas *= (h.time)^3.0/units.HubbleParam
  nfw_gas.rho_gas *= (h.time)^3.0/units.HubbleParam
  nfw_gas.rho_gas_iso *= (h.time)^3.0/units.HubbleParam
  
  ; plot (1) - histograms
  start_PS, sP.plotPath + 'timescales_histo.' + hTag + '.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps'
    
    w = where(coolTime eq 0.0,count)
    logCt = alog10(coolTime)
    logCt[w] = 0.0
    hist = histogram(logCt,binsize=binsize,loc=loc)
    cgPlot,10.0^loc,hist/float(total(hist)),color=cgColor('red'),/xs,$
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
  
  ; plot (5) - cooling/dynamical vs radius (temp/dens colormap)
  start_PS, sP.plotPath + 'timescales_vsrad1.' + hTag + '.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps'
    cgPlot,[0],[0],/nodata,xtitle=textoidl("r / r_{vir}"),ytitle="Timescale [Gyr]",$
      xrange=radRange,yrange=coolingRange,/ylog,/xs,/ys,yminor=0,$
      title=hTag + ' (z='+string(sP.redshift,format='(f3.1)')+')',position=[0.15,0.13,0.85,0.9]
    
    ; individual gas elements: cooling time color mapped by temperature
    loadColorTable,'blue-red';,/reverse
    TVLCT, rr, gg, bb, /GET
    
    ; temp colormap:
    fieldMinMax = [5.0,6.5]
    colorinds = (curTemp-fieldMinMax[0])*235.0 / (fieldMinMax[1]-fieldMinMax[0]) ;0-235
    
    ; dens colormap:
    ;fieldMinMax = [-8.0,-6.0]
    ;colorinds = (alog10(curDens)-fieldMinMax[0])*235.0 / (fieldMinMax[1]-fieldMinMax[0]) ;0-235
    
    ; individual gas elements
    colorinds = fix(colorinds + 20.0) > 0 < 255 ;20-255
    colorinds = getColor24([[rr[colorinds]], [gg[colorinds]], [bb[colorinds]]])
    
    for i=0L,n_elements(coolTime)-1 do $
      oplot,[gasRadii[i]/gasRvir[i]],[coolTime[i]],psym=psym,color=colorinds[i]
    
    ;cgPlot,gasRadii/gasRvir,coolTime,psym=psym,color=cgColor('red'),/overplot
    
    ; dynTime(r) and hubbleTime
    cgPlot,gasRadii/gasRvir,dynTime,psym=psym,color=cgColor('blue'),/overplot
    cgPlot,[0.05,1.5],[hubbleTime,hubbleTime],color=cgColor('green'),/overplot
    
    ; mean halo
    cgPlot,[0.05,1.0],[coolTime_halo,coolTime_halo],color=cgColor('magenta'),line=0,/overplot
    cgPlot,[0.05,1.0],[dynTime_halo,dynTime_halo],color=cgColor('orange'),line=0,/overplot
     
    ; radial fits
    cgPlot,radCt.binCen,radCt.radMedian,line=0,color=cgColor('magenta'),/overplot
    cgPlot,radDt.binCen,radDt.radMedian,line=0,color=cgColor('orange'),/overplot
     
    ; SIS/NFW models
    cgPlot,x,sis_gas.coolTime,line=0,color=cgColor('brown'),/overplot
    cgPlot,[0.15,1.5],[sis_gas.dynTime,sis_gas.dynTime],line=0,color=cgColor('pink'),/overplot
    
    cgPlot,x,nfw_gas.dynTime,line=0,color=cgColor('dark gray'),/overplot
    cgPlot,x,nfw_gas.coolTime,line=0,color=cgColor('forest green'),/overplot
     
    ; mark various r_cool
    cgPlot,[sis_gas.r_cool_h/sis_dm.r200,sis_gas.r_cool_h/sis_dm.r200],[0.06,0.1],line=2,color=cgColor('brown'),/overplot
    cgPlot,[sis_gas.r_cool/sis_dm.r200,sis_gas.r_cool/sis_dm.r200],[0.06,0.1],line=0,color=cgColor('brown'),/overplot
     
    cgPlot,[nfw_gas.r_cool/nfw_dm.r200,nfw_gas.r_cool/nfw_dm.r200],[0.06,0.1],line=0,color=cgColor('dark gray'), /overplot
     
    ; legend
    legend,textoidl(['t_{cool}','t_{dyn}','t_{age}','t_{cool,halo}','t_{dyn,halo}',$
                     't_{cool,sis}','t_{dyn,sis}','t_{cool,nfw}','t_{dyn,nfw}']),$
      textcolors=['red','blue','green','magenta','orange','brown','pink','forest green','dark gray'],$
      /bottom,/right,box=0
      
    ; colorbar
    cgColorbar,/right,/vertical,position=[0.88,0.13,0.93,0.9],range=fieldMinMax,title="log Temp"
  end_PS
  
  ; plot (6) - gas density/temperature vs radius
  start_PS, sP.plotPath + 'timescales_denstemp_rad.' + hTag + '.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps',xs=6,ys=8
    ; temp
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="log(T) [K]",$
      xrange=radRange,yrange=tempRange+[0.1,0.0],/xs,/ys,/xlog,xminor=0,$
      position=[0.15,0.55,0.95,0.95],xtickname=replicate(' ',10)
      
    ; individual gas elements and median binned
    cgPlot,gasRadii/gasRvir,curTemp,psym=psym,color=cgColor('blue'),/overplot
    cgPlot,radTemp.binCen,radTemp.radMedian,psym=-4,color=cgColor('black'),/overplot
    
    ; SIS model and SIS fit
    cgPlot,x,alog10(sis_gas.temp_gas),line=0,color=cgColor('forest green'),/overplot
    cgPlot,x,sis_fit.f_temp,line=2,color=cgColor('green'),/overplot
    
    ; NFW iso
    cgPlot,x,alog10(replicate(nfw_gas.T_0,n_elements(x))),line=0,color=cgColor('red'),/overplot
    
    ; NFW model and NFW fit (poly)
    cgPlot,x,alog10(nfw_gas.temp_gas),line=0,color=cgColor('orange'),/overplot    
    cgPlot,x,nfw_fit.f_temp,line=2,color=cgColor('orange'),/overplot
    
    ; dens
    cgPlot,[0],[0],/nodata,xtitle=textoidl("r / r_{vir}"),ytitle="Density [Code]",$
      xrange=radRange,yrange=densRange,/ylog,/xs,/ys,/xlog,yminor=0,xminor=0,$
      /noerase,position = [0.15,0.15,0.95,0.55]
      
    ; individual gas elements and median binned
    cgPlot,gasRadii/gasRvir,curDens,psym=psym,color=cgColor('blue'),/overplot
    cgPlot,radDens.binCen,radDens.radMedian,psym=-4,color=cgColor('black'),/overplot
    
    ; SIS model and SIS fit
    cgPlot,x,sis_gas.rho_gas,line=0,color=cgColor('forest green'),/overplot
    cgPlot,x,sis_fit.f_dens,line=2,color=cgColor('forest green'),/overplot
    
    ; NFW iso dens
    cgPlot,x,nfw_gas.rho_gas_iso,line=0,color=cgColor('red'),/overplot
    
    ; NFW model and NFW fit (poly)
    cgPlot,x,nfw_gas.rho_gas,line=0,color=cgColor('orange'),/overplot
    cgPlot,x,nfw_fit.f_dens,line=2,color=cgColor('orange'),/overplot
    
    legend,['median','SIS','NFW iso','NFW poly'],textcolors=['black','forest green','red','orange'],$
      box=0,/bottom,/left
    legend,['model','fit'],linestyle=[0,2],linesize=0.25,box=0,/top,/right
    
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
