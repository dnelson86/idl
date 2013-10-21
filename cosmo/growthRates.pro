; growthRates.pro
; feedback project - calculate actual growth rates of galaxies/halos, compare to net inflow rates
; dnelson sep.2013

pro growthRates;, sP=sP, timeWindow=timeWindow
  forward_function cosmoTracerChildren, cosmoTracerVelParents
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; config
  sP=simParams(res=512,run='gadget',redshift=2.0)
  timeWindow = 500.0 ; Myr or 'all'
  
  ; load net accretion rates for these same galaxies over the same timewindow (must be accMode=all)
  mbv = haloMassBinValues(sP=sP,accMode='all',timeWindow=timeWindow)
  
  if sP.trMCPerCell le 0 then massPerPart = sP.targetGasMass ; SPH or vel tracer
  if sP.trMCPerCell gt 0 then massPerPart = sP.trMassConst ; MC tracer
  
  ; which snapshot corresponds to the beginning of the timeWindow?
  if str(timeWindow) eq 'all' then begin
    targetSnap = redshiftToSnapNum(6.0,sP=sP)
    timeWindow = redshiftToAgeFlat(sP=sP) - redshiftToAgeFlat(6.0)
  endif else begin
    ages = redshiftToAgeFlat(snapNumToRedshift(sP=sP,/all))
    curAge = redshiftToAgeFlat(sP.redshift)
    targetSnap = closest( curAge-ages, timeWindow/1e3 )
  endelse
  
  actualTW = snapNumToAgeFlat(snap=sP.snap, sP=sP) - snapNumToAgeFlat(snap=targetSnap, sP=sP)
  actualTW *= 1e3 ; Myr

  print,'Snapshot range: ['+str(sP.snap)+'] to ['+str(targetSnap)+']'
  print,'Actual timeWindow: ['+str(actualTW)+'] vs requested ['+str(timeWindow)+']'
  
  ; load galaxyCatalog at this end of the timeWindow
  galcatCur = galaxyCat(sP=sP)
  
  countCur = total(galcatCur.len gt 0,/int)
  
  ; save arrays
  curMasses  = { galaxy : fltarr(galcatCur.nGroups), halo : fltarr(galcatCur.nGroups) }
  targMasses = { galaxy : fltarr(galcatCur.nGroups), halo : fltarr(galcatCur.nGroups) }
  
  ; calculate current masses: SPH/trMC/trVel
  if sP.trMCPerCell eq 0 then childCounts = intarr(galcatCur.countTot) + 1
  if sP.trMCPerCell gt 0 then childCounts = galcatCur.trMC_cc
  if sP.trMCPerCell lt 0 then childCounts = galcatCur.trVel_cc  
    
  for i=0L,galcatCur.nGroups-1 do begin
    if galcatCur.len[i] eq 0 then continue
    ; total mass = number of child tracers/SPH particles * const mass
    locInds = lindgen(galcatCur.len[i]) + galcatCur.off[i]
      
    galaxyInds = where(galcatCur.type[locInds] eq 1 or $
                       galcatCur.type[locInds] eq 4 or $
                       galcatCur.type[locInds] eq 5,countGal)
    haloInds   = where(galcatCur.type[locInds] eq 2,countHalo) ; type 3 = inter
      
    if countGal gt 0 then $
      curMasses.galaxy[i] = total(childCounts[locInds[galaxyInds]],/int) * massPerPart
    if countHalo gt 0 then $
      curMasses.halo[i] = total(childCounts[locInds[haloInds]],/int) * massPerPart
  endfor
  
  ; walk back through snapshots, determine indices of main progenitors at targetSnap
  origParIDs = lindgen(galcatCur.nGroups)
  
  for m=sP.snap,targetSnap+1,-1 do begin
    sP.snap = m
    
    ; how many remaining galaxies (with original length>0) remain?
    w = where(origParIDs ne -1 and galcatCur.len gt 0,countRemain)
    
    print,'['+str(m)+'] Remain: '+str(countRemain)+' '+$
      string(float(countRemain)/countCur*100,format='(f4.1)')
    
    ; load mergerTree parents
    Parent = mergerTree(sP=sP)
    
    ; shift indices to parent in previous snapshot
    w = where(origParIDs ne -1,count)
    if count eq 0 then message,'error'
    origParIDs[w] = Parent[origParIDs[w]] ; change to parent IDs    
    
  endfor
  
  ; load galaxyCatalog at this end of the timeWindow
  sP.snap -= 1
  galcatTarget = galaxyCat(sP=sP)
  
  ; calculate target masses: SPH/trMC/trVel
  if sP.trMCPerCell eq 0 then childCounts = intarr(galcatTarget.countTot) + 1
  if sP.trMCPerCell gt 0 then childCounts = galcatTarget.trMC_cc
  if sP.trMCPerCell lt 0 then childCounts = galcatTarget.trVel_cc  
    
  for i=0L,galcatCur.nGroups-1 do begin
    tInd = origParIDs[i]
    if tInd eq -1 or galcatTarget.len[tInd] eq 0 then continue
    
    ; total mass = number of child tracers/SPH particles * const mass
    locInds = lindgen(galcatTarget.len[tInd]) + galcatTarget.off[tInd]
      
    galaxyInds = where(galcatTarget.type[locInds] eq 1 or $
                       galcatTarget.type[locInds] eq 4 or $
                       galcatTarget.type[locInds] eq 5,countGal)
    haloInds   = where(galcatTarget.type[locInds] eq 2,countHalo) ; type 3 = inter
      
    if countGal gt 0 then $
      targMasses.galaxy[i] = total(childCounts[locInds[galaxyInds]],/int) * massPerPart
    if countHalo gt 0 then $
      targMasses.halo[i] = total(childCounts[locInds[haloInds]],/int) * massPerPart
  endfor
  
  ; calculate mass differences of galaxy(+stars) and halo
  massDiffs   = { galaxy : fltarr(galcatCur.nGroups) + !values.f_nan ,$
                  halo   : fltarr(galcatCur.nGroups) + !values.f_nan }
  growthRates = { galaxy : fltarr(galcatCur.nGroups) + !values.f_nan ,$
                  halo   : fltarr(galcatCur.nGroups) + !values.f_nan }
                
  wValidGal  = where(curMasses.galaxy gt 0 and targMasses.galaxy gt 0,countGal)
  wValidHalo = where(curMasses.halo gt 0 and targMasses.halo gt 0,countHalo)
  print,'Found: '+str(countGal)+' galaxies, '+str(countHalo)+' halos, of '+str(countCur)+' total.'
  
  massDiffs.galaxy[wValidGal] = curMasses.galaxy[wValidGal] - targMasses.galaxy[wValidGal]
  growthRates.galaxy[wValidGal] = (massDiffs.galaxy[wValidGal] * units.UnitMass_in_Msun) / (actualTW*1e6) ; msun/yr
  
  massDiffs.halo[wValidHalo] = curMasses.halo[wValidHalo] - targMasses.halo[wValidHalo]
  growthRates.halo[wValidHalo] = (massDiffs.halo[wValidHalo] * units.UnitMass_in_Msun) / (actualTW*1e6) ; msun/yr
  
  ; plots
  xyrange = [0.05,10.0]
  linexy  = [0.07,8.0]
  
  start_PS,sP.plotPath + 'mass_vs_mass_gal_'+sP.plotPrefix+str(sP.res)+'.eps'
    cgPlot,curMasses.galaxy[wValidGal],targMasses.galaxy[wValidGal],psym=4,$
      xtitle="Gal Mass Current",ytitle="Gal Mass Target",$
      /xlog,/ylog,xminor=0,yminor=0,xrange=xyrange,yrange=xyrange
    cgPlot,linexy,linexy,line=0,color=cgColor('orange'),/overplot
  end_PS
  
  start_PS,sP.plotPath + 'mass_vs_mass_halo_'+sP.plotPrefix+str(sP.res)+'.eps'
    cgPlot,curMasses.halo[wValidHalo],targMasses.halo[wValidHalo],psym=4,$
      xtitle="Halo Mass Current",ytitle="Halo Mass Target",$
      /xlog,/ylog,xminor=0,yminor=0,xrange=xyrange,yrange=xyrange
    cgPlot,linexy,linexy,line=0,color=cgColor('orange'),/overplot
  end_PS
  
  ; how do the actual growthRates compare to the net inflow rates over the same tW?
  ; 0 index = 1.0 tvir (total is independent of this choice)
  mbvGal  = reform( mbv.galaxy.netRate.total.hot.tVirAcc[0,*] ) + $
            reform( mbv.galaxy.netRate.total.cold.tVirAcc[0,*] )
  mbvHalo = reform( mbv.halo.netRate.total.hot.tVirAcc[0,*] ) + $
            reform( mbv.halo.netRate.total.cold.tVirAcc[0,*] )
  
  xyrange = [0.05,100]
  linexy  = [0.07,70]
  
  start_PS,sP.plotPath + 'growthRate_vs_netRate_gal_'+sP.plotPrefix+str(sP.res)+'.eps'
    cgPlot,growthRates.galaxy[wValidGal],mbvGal[wValidGal],psym=4,$
      xrange=xyrange,yrange=xyrange,/xs,/ys,/xlog,/ylog,xminor=0,yminor=0,$
      xtitle="Galaxy Growth Rate [Msun/yr]",ytitle="Net Gal Inflow Rate [Msun/yr]",$
      title="tW: "+str(timeWindow)+" "+sP.run+" "+str(sP.res)
      
    cgPlot,linexy,linexy,line=0,color=cgColor('orange'),/overplot
    
    offset = mean(mbvGal[wValidGal] / growthRates.galaxy[wValidGal],/nan)
    cgPlot,linexy,linexy*offset,line=0,color=cgColor('blue'),/overplot
    
    legend,['1-to-1','offset = '+str(offset)],textcolor=['orange','blue'],/bottom,/right
  end_PS
  
  start_PS,sP.plotPath + 'growthRate_vs_netRate_halo_'+sP.plotPrefix+str(sP.res)+'.eps'
    cgPlot,growthRates.halo[wValidHalo],mbvHalo[wValidHalo],psym=4,$
      xrange=xyrange,yrange=xyrange,/xs,/ys,/xlog,/ylog,xminor=0,yminor=0,$
      xtitle="Halo Growth Rate [Msun/yr]",ytitle="Net Halo Inflow Rate [Msun/yr]",$
      title="tW: "+str(timeWindow)+" "+sP.run+" "+str(sP.res)
    
    cgPlot,linexy,linexy,line=0,color=cgColor('orange'),/overplot
    
    offset = mean(mbvHalo[wValidHalo] / growthRates.halo[wValidHalo],/nan)
    cgPlot,linexy,linexy*offset,line=0,color=cgColor('blue'),/overplot
    
    legend,['1-to-1','offset = '+str(offset)],textcolor=['orange','blue'],/bottom,/right
    
  end_PS
  
  stop
  
end