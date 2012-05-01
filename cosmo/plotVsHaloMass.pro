; plotVsHaloMass.pro
; gas accretion project - plots as a function of halo mass
; dnelson may.2012

; haloMassBinDeltaAccTime(): bin mean time for accreting gas as a function of parent halo mass

function haloMassBinDeltaAccTime, sP=sP, sgSelect=sgSelect, accMode=accMode

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  xrange = [9.5,12.5]
  logMassBinSize = 0.25

  ; load
  at = accretionTimes(sP=sP)
  mt = mergerTreeSubset(sP=sP)
  
  ; get indices and mask for either smooth,clumpy,all
  wAm = accModeInds(at=at,accMode=accMode,sP=sP,/mask)
  
  gcIndOrig = mergerTreeRepParentIDs(mt=mt,sP=sP)
  
  ; get parent mass of each gas element (all of mtS since accretion time enforced later)
  parMasses = gcSubsetProp(sP=sP,select=sgSelect,/parMass,/mergerTreeSubset,accMode=accMode)

  ; convert (k=0) r=rvir accretion times to nearest snap
  snapTimes = snapNumToRedshift(snap=0,/time,/all,sP=sP)
  snapTimes = snapTimes[where(snapTimes ge 0)] ; remove -1
    
  accSnap = {gal  : intarr(n_elements(at.accTime_gal[0,*]))   ,$
             gmem : intarr(n_elements(at.accTime_gmem[0,*]))   }
             
  accSnap.gal[wAm.gal]   = value_locate(snapTimes,at.accTime_gal[0,wAm.gal]) 
  accSnap.gmem[wAm.gmem] = value_locate(snapTimes,at.accTime_gmem[0,wAm.gmem])
  
  w = where(accSnap.gal[wAm.gal] le 0 or accSnap.gal[wAm.gal] gt sP.snapRange[1],count1)
  w = where(accSnap.gmem[wAm.gmem] le 0 or accSnap.gmem[wAm.gmem] gt sP.snapRange[1],count2)
  if count1 gt 0 or count2 gt 0 then message,'Error: Bad snapshot time mapping.'
    
  ; each accDt now stores the scale factor when that particle crossed the radius of interest
  ; convert to delta time (Myr) since crossing rvir
  for k=1,n_elements(at.rVirFacs)-1 do begin
    gal_w  = where(wAm.galMask eq 1B and at.accTime_gal[k,*] ne -1,count_gal)
    gmem_w = where(wAm.gmemMask eq 1B and at.accTime_gmem[k,*] ne -1,count_gmem)

    if count_gal gt 0 then begin
      ; for r=rvir, convert scale factor to redshift to age of universe in Gyr
      galAgeRvir  = redshiftToAgeFlat(1/reform(at.accTime_gal[0,gal_w])-1)
      
      ; for r=rVirFac, convert scale factor to redshift to age of universe in Gyr
      galAgeRfac  = redshiftToAgeFlat(1/reform(at.accTime_gal[k,gal_w])-1)
      
      ; delta time in Myr
      dt_gal  = (galAgeRfac - galAgeRvir) ; Gyr    
      
      w = where(dt_gal le 0.0,count)
      if count gt 0 then message,'Error: Negative gal times.'
      
      ; double check we have good accretion snapshots
      w = where(accSnap.gal[gal_w] le 0 or accSnap.gal[gal_w] gt sP.snapRange[1],count)
      if count gt 0 then message,'Error: Bad snapshot time mapping gal two.'    
      
      ; normalize by t_circ for parent halos given their mass and rvir crossing times
      parSnapMasses = mt.hMass[accSnap.gal[gal_w]-mt.minSnap,gcIndOrig.gal[gal_w]]
      parVirRad = mt.hVirRad[accSnap.gal[gal_w]-mt.minSnap,gcIndOrig.gal[gal_w]]
      circTime  = reform(2*!pi * sqrt(parVirRad^3.0 / units.G / parSnapMasses)) ; Gyr
      
      at.accTime_gal[k,gal_w] = dt_gal / circTime    
    endif
    
    if count_gmem gt 0 then begin
      ; repeat for gmem
      gmemAgeRvir = redshiftToAgeFlat(1/reform(at.accTime_gmem[0,gmem_w])-1)
      gmemAgeRfac = redshiftToAgeFlat(1/reform(at.accTime_gmem[k,gmem_w])-1)
  
      dt_gmem = (gmemAgeRfac - gmemAgeRvir)
  
      w = where(dt_gmem le 0.0,count)
      if count gt 0 then message,'Error: Negative gmem times.'
      
      w = where(accSnap.gmem[gmem_w] le 0 or accSnap.gmem[gmem_w] gt sP.snapRange[1],count)
      if count gt 0 then message,'Error: Bad snapshot time mapping gmem two.'
  
      parSnapMasses = mt.hMass[accSnap.gmem[gmem_w]-mt.minSnap,gcIndOrig.gmem[gmem_w]]
      parVirRad = mt.hVirRad[accSnap.gmem[gmem_w]-mt.minSnap,gcIndOrig.gmem[gmem_w]]
      circTime  = reform(2*!pi * sqrt(parVirRad^3.0 / units.G / parSnapMasses))
      
      at.accTime_gmem[k,gmem_w] = dt_gmem / circTime
    endif
  endfor

  ; bin fractions into halo mass bins and make median lines
  logMassNbins  = floor((xrange[1]-xrange[0]) / logMassBinSize)
  logMassBins   = linspace(xrange[0],xrange[1],logMassNbins+1) ; edges
  logMassBinCen = linspace(xrange[0],xrange[1],logMassNbins+1) + logMassBinSize/2.0
  logMassBinCen = logMassBinCen[0:-2] ; remove last

  binnedVals = { mean_gal    : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                 mean_gmem   : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                 mean_both   : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                 median_gal  : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                 median_gmem : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                 median_both : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                 stddev_gal  : fltarr(n_elements(at.rVirFacs),logMassNBins) + !values.f_nan ,$
                 stddev_gmem : fltarr(n_elements(at.rVirFacs),logMassNBins) + !values.f_nan ,$
                 stddev_both : fltarr(n_elements(at.rVirFacs),logMassNBins) + !values.f_nan ,$
                 logMassBinCen : logMassBinCen ,$
                 rVirFacs      : at.rVirFacs       }
  
  ; calculate median in bins (Tmax) and enforce minimum particle numbers
  for i=0,logMassNbins-2 do begin
    for k=1,n_elements(at.rVirFacs)-1 do begin
      ; gal
      w_gal = where(parMasses.gal gt logMassBins[i] and parMasses.gal le logMassBins[i+1] and $
                    at.accTime_gal[0,*] ne -1 and at.accTime_gal[k,*] ne -1,count_gal)
      if count_gal gt 0 then begin
        binnedVals.mean_gal[k,i]   = mean(at.accTime_gal[k,w_gal])
        binnedVals.median_gal[k,i] = median(at.accTime_gal[k,w_gal])
        binnedVals.stddev_gal[k,i] = stddev(at.accTime_gal[k,w_gal])
      endif
      
      ; gmem
      w_gmem = where(parMasses.gmem gt logMassBins[i] and parMasses.gmem le logMassBins[i+1] and $
                     at.accTime_gmem[0,*] ne -1 and at.accTime_gmem[k,*] ne -1,count_gmem)
      if count_gmem gt 0 then begin
        binnedVals.mean_gmem[k,i]   = mean(at.accTime_gmem[k,w_gmem])
        binnedVals.median_gmem[k,i] = median(at.accTime_gmem[k,w_gmem])
        binnedVals.stddev_gmem[k,i] = stddev(at.accTime_gmem[k,w_gmem])
      endif
      
      ; both
      vals = []
      if count_gal gt 0  then vals = [vals,reform(at.accTime_gal[k,w_gal])]
      if count_gmem gt 0 then vals = [vals,reform(at.accTime_gmem[k,w_gmem])]
      
      if (count_gal+count_gmem) gt 0 then begin
        binnedVals.mean_both[k,i]   = mean(vals)
        binnedVals.median_both[k,i] = median(vals)
        binnedVals.stddev_both[k,i] = stddev(vals)
      endif
    endfor ;k
  endfor ;i
  
  ; debug plot
  start_PS,sP.plotPath + 'accdt.histos.'+accMode+'.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    binsize = 0.1 & strings = []
    cgPlot,[0],[0],/nodata,xrange=[0.0,0.5],yrange=[0.01,0.6],/ylog,$
      xtitle=textoidl("\Delta t_{acc} / \tau_{circ}"),ytitle="Fraction"
      
    for k=1,n_elements(at.rVirFacs)-1 do begin
      w = where(at.accTime_gal[k,*] ne -1,count)
      if count gt 1 then begin
        hist = histogram(at.accTime_gal[k,w],binsize=binsize,loc=loc)
        cgplot,loc+binsize*0.5,hist/float(total(hist)),line=0,color=getColor(k),/overplot
      endif
      w = where(at.accTime_gmem[k,*] ne -1,count)
      if count gt 1 then begin
        hist = histogram(at.accTime_gmem[k,w],binsize=binsize,loc=loc)
        cgplot,loc+binsize*0.5,hist/float(total(hist)),line=1,color=getColor(k),/overplot
      endif
      strings = [strings,string(at.rVirFacs[k],format='(f4.2)')]
    endfor
    
    legend,strings,textcolors=getColor([1,2,3,4,5,6],/name),box=0,/top,/right
  end_PS
  
  return,binnedVals
  
end

; haloMassBinModeMasses(): bin total mass in hot/cold modes as a function of halo mass

function haloMassBinModeMasses, sP=sP, sgSelect=sgSelect, accMode=accMode, KDE=KDE
  
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  TcutVals = [5.3,5.4,5.5,5.6,5.7] ; for constant threshold
  xrange   = [10.0,12.5] ; mass bin range
  
  logMassBinSize = 0.5 / (sP.res/128)
  nCuts = n_elements(TcutVals)
  
  ; make a uniform gas selection at the start
  at = accretionTimes(sP=sP)
  mt = mergerTreeSubset(sP=sP)
  
  wAm = accModeInds(at=at,accMode=accMode,sP=sP)
    
  if sP.trMCPerCell eq 0 then massGasPart = sP.targetGasMass
  if sP.trMCPerCell ne 0 then massGasPart = sP.trMassConst

  ; load max temps, current tvir, tvir at accretion
  accTvir = gcSubsetProp(sP=sP,select=sgSelect,/accTvir,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  curTvir = gcSubsetProp(sP=sP,select=sgSelect,/virTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  maxTemp = gcSubsetProp(sP=sP,select=sgSelect,/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)

  ; load parent mass for each gas element
  parMasses = gcSubsetProp(sP=sP,select=sgSelect,/parMass,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  
  if n_elements(accTvir.gal) ne n_elements(curTvir.gal) or $
     n_elements(accTvir.gal) ne n_elements(maxTemp.gal) or $
     n_elements(accTvir.gal) ne n_elements(parMasses.gal) or $
     n_elements(accTvir.gal) ne n_elements(wAm.gal) then message,'Data error'
  
  ; bin fractions into halo mass bins and make median lines
  logMassNbins  = floor((xrange[1]-xrange[0]) / logMassBinSize)
  logMassBins   = reverse(linspace(xrange[0],xrange[1],logMassNbins+1)) ; edges, descending
  logMassBinCen = linspace(xrange[0],xrange[1],logMassNbins+1) + logMassBinSize/2.0
  logMassBinCen = reverse(logMassBinCen[0:-2]) ; remove last, descending
  
  coldMass = { const_gal    : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
               const_gmem   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
               tVircur_gal  : fltarr(logMassNbins)       + !values.f_nan ,$
               tVircur_gmem : fltarr(logMassNbins)       + !values.f_nan ,$
               tViracc_gal  : fltarr(logMassNbins)       + !values.f_nan ,$
               tViracc_gmem : fltarr(logMassNbins)       + !values.f_nan  }
  
  hotMass = { const_gal    : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
              const_gmem   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
              tVircur_gal  : fltarr(logMassNbins)       + !values.f_nan ,$
              tVircur_gmem : fltarr(logMassNbins)       + !values.f_nan ,$
              tViracc_gal  : fltarr(logMassNbins)       + !values.f_nan ,$
              tViracc_gmem : fltarr(logMassNbins)       + !values.f_nan  }
  
  totalMass = { const_gal    : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                const_gmem   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                tVircur_gal  : fltarr(logMassNbins)       + !values.f_nan ,$
                tVircur_gmem : fltarr(logMassNbins)       + !values.f_nan ,$
                tViracc_gal  : fltarr(logMassNbins)       + !values.f_nan ,$
                tViracc_gmem : fltarr(logMassNbins)       + !values.f_nan  }
                   
  ; calculate total in each halo mass bin
  for i=0,logMassNbins-2 do begin
    ; gal (total mass)
    loc_inds = where(parMasses.gal lt logMassBins[i] and parMasses.gal ge logMassBins[i+1],count)
    if count gt 0 then begin
      loc_maxt = maxTemp.gal[loc_inds]
      
      ; count fraction Tmax below each constant temperature threshold
      for j=0,nCuts-1 do begin
        w = where(loc_maxt le TcutVals[j],count_below,ncomp=count_above)
        coldMass.const_gal[j,i] = count_below * massGasPart
        hotMass.const_gal[j,i]    = count_above * massGasPart
        totalMass.const_gal[j,i]  = (count_below+count_above) * massGasPart
        if (count_below+count_above) ne count then stop
      endfor
      
      ; count fraction Tmax below Tvir at current time
      w = where(loc_maxt le curTvir.gal[loc_inds],count_below,ncomp=count_above)
      coldMass.tVircur_gal[i]  = count_below * massGasPart
      hotMass.tVircur_gal[i]   = count_above * massGasPart
      totalMass.tVircur_gal[i] = (count_below+count_above) * massGasPart
      if (count_below+count_above) ne count then stop
  
      ; count fraction Tmax below Tvir at accretion time
      w = where(loc_maxt le accTvir.gal[loc_inds],count_below,ncomp=count_above)
      coldMass.tViracc_gal[i]  = count_below * massGasPart
      hotMass.tViracc_gal[i]   = count_above * massGasPart
      totalMass.tViracc_gal[i] = (count_below+count_above) * massGasPart
      if (count_below+count_above) ne count then stop
    endif
    
    ; gmem (total mass)
    loc_inds = where(parMasses.gmem lt logMassBins[i] and parMasses.gmem ge logMassBins[i+1],count)
    if count gt 0 then begin
      loc_maxt = maxTemp.gmem[loc_inds]
      
      ; count fraction Tmax below each constant temperature threshold
      for j=0,nCuts-1 do begin
        w = where(loc_maxt le TcutVals[j],count_below,ncomp=count_above)
        coldMass.const_gmem[j,i]  = count_below * massGasPart
        hotMass.const_gmem[j,i]   = count_above * massGasPart
        totalMass.const_gmem[j,i] = (count_below+count_above) * massGasPart
        if (count_below+count_above) ne count then stop
      endfor
      
      ; count fraction Tmax below Tvir at current time
      w = where(loc_maxt le curTvir.gmem[loc_inds],count_below,ncomp=count_above)
      coldMass.tVircur_gmem[i]  = count_below * massGasPart
      hotMass.tVircur_gmem[i]   = count_above * massGasPart
      totalMass.tVircur_gmem[i] = (count_below+count_above) * massGasPart
      if (count_below+count_above) ne count then stop
  
      ; count fraction Tmax below Tvir at accretion time
      w = where(loc_maxt le accTvir.gmem[loc_inds],count_below,ncomp=count_above)
      coldMass.tViracc_gmem[i]  = count_below * massGasPart
      hotMass.tViracc_gmem[i]   = count_above * massGasPart
      totalMass.tViracc_gmem[i] = (count_below+count_above) * massGasPart
      if (count_below+count_above) ne count then stop
    endif
  endfor
  
  ; normalize histograms by binsize and convert to log msun
  coldMass.tViracc_gal[*]  = coldMass.tVircur_gal[*] / logMassBinSize
  hotMass.tViracc_gal[*]   = hotMass.tVircur_gal[*] / logMassBinSize
  totalMass.tViracc_gal[*] = totalMass.tVircur_gal[*] / logMassBinSize
  
  coldMass.tViracc_gmem[*]  = coldMass.tViracc_gmem[*] / logMassBinSize
  hotMass.tViracc_gmem[*]   = hotMass.tViracc_gmem[*] / logMassBinSize
  totalMass.tViracc_gmem[*] = totalMass.tViracc_gmem[*] / logMassBinSize
  
  ; test: kernel density estimator
  if keyword_set(KDE) then begin
    xpts = findgen(50)/50.0*2.0 + 10.5
    spacing = 2.0/100
    ypts = parMasses.gal[where(maxTemp.gal lt accTvir.gal)]
    tt_gal = akde(ypts,xpts)
    tt_gal = tt_gal * spacing * n_elements(ypts) * logMassBinSize / 2.0
    
    ypts = parMasses.gmem[where(maxTemp.gmem lt accTvir.gmem)]
    tt_gmem = akde(ypts,xpts)
    tt_gmem = tt_gmem * spacing * n_elements(ypts) * logMassBinSize / 2.0
    print,'Warning: This de-normalization is wrong.'
  endif else begin
    tt_gal = 0 & tt_gmem = 0
  endelse
  
  r = {coldMass:coldMass,hotMass:hotMass,totalMass:totalMass,logMassBinCen:logMassBinCen,xrange:xrange,$
       tt_gal:tt_gal,tt_gmem:tt_gmem}
  return,r
end

; plotDeltaAccTimeVsHaloMass(): plot the mean time for accreting gas to reach various radii from the 
;                               virial radius (normalized by t_circ)

pro plotDeltaAccTimeVsHaloMass

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect = 'pri'
  accMode = 'all'
  
  res = 128
  sPg = simParams(res=res,run='gadget',redshift=2.0)
  ;sPa = simParams(res=res,run='tracer',redshift=2.0)
  
  binnedGadget = haloMassBinDeltaAccTime(sP=sPg,sgSelect=sgSelect,accMode=accMode)
  ;binnedArepo  = haloMassBinDeltaAccTime(sP=sPa,sgSelect=sgSelect,accMode=accMode)

  ; plot (1) - lines for all rVirFacs vs halo mass
  start_PS, sPg.plotPath + 'accdt.vshalo.comp.'+accMode+'.'+str(sPg.res)+'_'+str(sPg.snap)+'.eps', /big
    
    xrange = [9.5,12.5]
    yrange = [0.0,0.28]
    sigFac = 2.0
    
    xtickv = [10,11,12]
    
    x0 = 0.15 & x1 = 0.42 & x2 = 0.69 & x3 = 0.96
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; uc
                [x2,y1,x3,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1] ,$ ; lc
                [x2,y0,x3,y1] )  ; lr
    
    for j=1,n_elements(binnedGadget.rVirFacs)-1 do begin
      
      if j eq 1 or j eq 4 then ytickname = '' else ytickname = replicate(' ',10)
      if j ge 4 then xtickname = '' else xtickname = replicate(' ',10)
      if j gt 1 then noerase = 1 else noerase = 0
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,pos=pos[j-1],$
        ytitle="",xtitle="",xticks=n_elements(xtickv)-1,xtickv=xtickv,xminor=2,$
        xtickname=xtickname,ytickname=ytickname,noerase=noerase       
      
      cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
      
      w = where(binnedGadget.logMassBinCen gt 8.0)
      
      ; gadget both
      cgPlot,binnedGadget.logMassBinCen[w],binnedGadget.mean_both[j,w],color=getColor(1),line=1,/overplot
      cgPlot,binnedGadget.logMassBinCen[w],binnedGadget.median_both[j,w],color=getColor(1),line=2,/overplot
      
      ; sigma envelope
      cgPlot,binnedGadget.logMassBinCen[w],$
             binnedGadget.mean_both[j,w]+sigFac*binnedGadget.stddev_both[j,w],$
             color=getColor(1),line=0,/overplot
      cgPlot,binnedGadget.logMassBinCen[w],$
             binnedGadget.mean_both[j,w]-sigFac*binnedGadget.stddev_both[j,w],$
             color=getColor(1),line=0,/overplot
      
      ; arepo both
      ;cgPlot,binnedArepo.logMassBinCen[w],binnedArepo.mean_both[j,w],color=getColor(3),line=1,/overplot
      ;cgPlot,binnedArepo.logMassBinCen[w],binnedArepo.median_both[j,w],color=getColor(3),line=2,/overplot
      
      ; sigma envelope
      ;cgPlot,binnedArepo.logMassBinCen[w],$
      ;       binnedArepo.mean_both[j,w]+sigFac*binnedArepo.stddev_both[j,w],$
      ;       color=getColor(3),line=0,/overplot
      ;cgPlot,binnedArepo.logMassBinCen[w],$
      ;       binnedArepo.mean_both[j,w]-sigFac*binnedArepo.stddev_both[j,w],$
      ;       color=getColor(3),line=0,/overplot
      
      ; legend
      massBinStr = textoidl('r/r_{vir}')+' = ' + string(binnedGadget.rVirFacs[j],format='(f4.2)')
      cgText,mean(xrange),yrange[1]*0.9,massBinStr,charsize=!p.charsize-0.2,alignment=0.5,$
        color=cgColor('forest green')
    
    endfor
    
    ; axis labels
    cgText,0.05,0.5,textoidl("\Delta t_{acc} / \tau_{circ}"),alignment=0.5,orientation=90.0,/normal
    cgText,0.5,0.05,textoidl("log ( M_{halo} ) [_{ }M_{sun }]"),alignment=0.5,/normal
    
    ; annotations
    cgText,0.95,0.06,"gadget",alignment=0.5,charsize=!p.charsize-0.2,color=getColor(1),/normal
    cgText,0.95,0.02,"arepo",alignment=0.5,charsize=!p.charsize-0.2,color=getColor(3),/normal
    
  end_PS
end

; plotModeMassesVsHaloMass(): plot the "cold mass" and "hot mass" (not fraction) vs halo mass

pro plotModeMassesVsHaloMass

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect = 'pri'
  accMode  = 'all' ; accretion mode: all, smooth, bclumpy, sclumpy
  
  res = 256
  sPg = simParams(res=res,run='gadget',redshift=2.0)
  sPa = simParams(res=res,run='tracer',redshift=2.0)
  
  mmG = haloMassBinModeMasses(sP=sPg,sgSelect=sgSelect,accMode=accMode)
  mmA = haloMassBinModeMasses(sP=sPa,sgSelect=sgSelect,accMode=accMode)
  
  ; plot - cold,hot,total masses
  xrange = [10.5,mmG.xrange[1]]
  w = where(mmG.logMassBinCen gt 10.5)
    
  ; plot (1) - tviracc only, cold+hot+total
  start_PS, sPg.plotPath + 'massBudget.'+accMode+'.'+sPg.run+'.'+str(sPg.res)+'_'+str(sPg.snap)+'.eps',$
    xs = 7.5, ys = 10
    
    x0 = 0.15 & x1 = 0.96
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; uc
                [x0,y0,x1,y1]  ) ; lc
   
    ; gal
    yrange = [0.0,max(mmA.totalMass.tviracc_gal,/nan)*1.02]
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0],xminor=2,yminor=2
    
    ;xpts = findgen(50)/50.0*2.0 + 10.5
    ;cgPlot,xpts,mmG.tt_gal,line=3,/overplot
    
    cgPlot,mmG.logMassBinCen[w],mmG.coldMass.tViracc_gal[w],color=getColor(1),line=1,/overplot
    cgPlot,mmG.logMassBinCen[w],mmG.hotMass.tViracc_gal[w],color=getColor(1),line=2,/overplot
    cgPlot,mmG.logMassBinCen[w],mmG.totalMass.tViracc_gal[w],color=getColor(1),line=0,/overplot
    
    cgPlot,mmA.logMassBinCen[w],mmA.coldMass.tViracc_gal[w],color=getColor(3),line=1,/overplot
    cgPlot,mmA.logMassBinCen[w],mmA.hotMass.tViracc_gal[w],color=getColor(3),line=2,/overplot
    cgPlot,mmA.logMassBinCen[w],mmA.totalMass.tViracc_gal[w],color=getColor(3),line=0,/overplot
    
    ; legend
    legend,['cold','hot','total'],linestyle=[1,2,0],box=0,linesize=0.25,/top,/right
    
    ; gmem
    yrange = [0.0,max(mmG.totalMass.tviracc_gmem,/nan)*1.02]
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle=textoidl("log( M_{halo} ) [_{ }M_{sun }]"),xminor=2,yminor=2,pos=pos[1],/noerase
    
    ;xpts = findgen(50)/50.0*2.0 + 10.5
    ;cgPlot,xpts,mmG.tt_gmem,line=3,/overplot
    
    cgPlot,mmG.logMassBinCen[w],mmG.coldMass.tViracc_gmem[w],color=getColor(1),line=1,/overplot
    cgPlot,mmG.logMassBinCen[w],mmG.hotMass.tViracc_gmem[w],color=getColor(1),line=2,/overplot
    cgPlot,mmG.logMassBinCen[w],mmG.totalMass.tViracc_gmem[w],color=getColor(1),line=0,/overplot
    
    cgPlot,mmA.logMassBinCen[w],mmA.coldMass.tViracc_gmem[w],color=getColor(3),line=1,/overplot
    cgPlot,mmA.logMassBinCen[w],mmA.hotMass.tViracc_gmem[w],color=getColor(3),line=2,/overplot
    cgPlot,mmA.logMassBinCen[w],mmA.totalMass.tViracc_gmem[w],color=getColor(3),line=0,/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=getColor([1,3],/name),box=0,/top,/right
    
    ; labels
    cgText,0.05,0.5,"Accreted Gas Mass "+textoidl("[_{ }M_{sun } / log( M_{halo} )]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,mean([x0,x1]),0.17,"Halo Atmosphere",alignment=0.5,color=cgColor('forest green'),/normal
    cgText,mean([x0,x1]),0.91,"Central Galaxy",alignment=0.5,color=cgColor('forest green'),/normal
    
  end_PS
  
  stop
end

; plotColdFracVsHaloMass(): plot the "cold mode fraction" vs halo mass in a few different ways
;                           using the mergerTreeSubset with accretionTimes and accretionModes

pro plotColdFracVsHaloMass, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect = 'pri'
  TcutVals = [5.3,5.4,5.5,5.6,5.7] ; for constant threshold
  nCuts = n_elements(TcutVals)
  
  accMode = 'all' ; accretion mode: all, smooth, bclumpy, sclumpy
  minNum = 6
  xrange = [9.5,12.5]
  yrange = [0.0,1.1]
  
  logMassBinSize = 0.15 / (sP.res/128)
  
  ; make a uniform gas selection at the start
  at = accretionTimes(sP=sP)
  mt = mergerTreeSubset(sP=sP)
  
  wAm = accModeInds(at=at,accMode=accMode,sP=sP)
    
  ; reverse histogram parent IDs of all particles/tracers in this selection
  gcIndOrig = mergerTreeRepParentIDs(mt=mt,sP=sP)
  
  hist_gal  = histogram(gcIndOrig.gal[wAm.gal],min=0,loc=loc_gal,rev=rev_gal)
  hist_gmem = histogram(gcIndOrig.gmem[wAm.gmem],min=0,loc=loc_gmem,rev=rev_gmem)
  
  ; load max temps, current tvir, tvir at accretion
  accTvir = gcSubsetProp(sP=sP,select=sgSelect,/accTvir,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  curTvir = gcSubsetProp(sP=sP,select=sgSelect,/virTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  maxTemp = gcSubsetProp(sP=sP,select=sgSelect,/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  
  ; load current temps and current SFR
  curTemp = gcSubsetProp(sP=sP,select=sgSelect,/curTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  curSFR  = gcSubsetProp(sP=sP,select=sgSelect,/curSingleVal,singleValField='sfr',$
                         /mergerTreeSubset,/accretionTimeSubset,accMode=accMode)

  ; load group cat for subgroup masses
  gc = loadGroupCat(sP=sP,/skipIDs)
  gcMasses = codeMassToLogMsun(gc.subgroupMass[mt.galcatIDList])
  gc = !NULL

  ; structures to store results (Tmax)
  coldFrac = { gal_const    : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
               gmem_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
               gal_tvircur  : fltarr(n_elements(mt.galcatIDList))         ,$
               gmem_tvircur : fltarr(n_elements(mt.galcatIDList))         ,$
               gal_tviracc  : fltarr(n_elements(mt.galcatIDList))         ,$
               gmem_tviracc : fltarr(n_elements(mt.galcatIDList))         ,$
               both_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
               both_tvircur : fltarr(n_elements(mt.galcatIDList))         ,$
               both_tviracc : fltarr(n_elements(mt.galcatIDList))         ,$
               gal_num      : lonarr(n_elements(mt.galcatIDList))         ,$
               gmem_num     : lonarr(n_elements(mt.galcatIDList))          }
  
  ; loop over all tracked subgroups (galaxy, Tmax)
  for i=0,n_elements(hist_gal)-1 do begin
    if hist_gal[i] gt 0 then begin
      ; list of indices of galaxy gas particles in this subgroup
      loc_inds_gal = rev_gal[rev_gal[i]:rev_gal[i+1]-1]
      ;loc_inds_gal2 = where(gcIndOrigTr.gal[gal_w] eq i)
      ;if ~array_equal(loc_inds_gal,loc_inds_gal2) then print,'error'
      
      loc_maxt_gal = maxTemp.gal[loc_inds_gal]
      nloc = n_elements(loc_maxt_gal)
      
      coldFrac.gal_num[i] = nloc
      
      ; count fraction Tmax below each constant temperature threshold
      for j=0,nCuts-1 do begin
        w = where(loc_maxt_gal le TcutVals[j],count_below)
        coldFrac.gal_const[j,i] = float(count_below) / nloc
      endfor
      
      ; count fraction Tmax below Tvir at current time
      w = where(loc_maxt_gal le curTvir.gal[loc_inds_gal],count_below)
      coldFrac.gal_tvircur[i] = float(count_below) / nloc

      ; count fraction Tmax below Tvir at accretion time
      w = where(loc_maxt_gal le accTvir.gal[loc_inds_gal],count_below)
      coldFrac.gal_tviracc[i] = float(count_below) / nloc
    endif
  endfor
  
  ; loop over all tracked subgroups (groupmem, Tmax)
  for i=0,n_elements(hist_gmem)-1 do begin
    if hist_gmem[i] gt 0 then begin
      ; list of indices of groupmem gas particles in this subgroup
      loc_inds_gmem = rev_gmem[rev_gmem[i]:rev_gmem[i+1]-1]
      loc_maxt_gmem = maxTemp.gmem[loc_inds_gmem]
      nloc = n_elements(loc_maxt_gmem)
      
      coldFrac.gmem_num[i] = nloc
      
      ; count fraction Tmax below each constant temperature threshold
      for j=0,nCuts-1 do begin
        w = where(loc_maxt_gmem le TcutVals[j],count_below)
        coldFrac.gmem_const[j,i] = float(count_below) / nloc
      endfor
      
      ; count fraction Tmax below Tvir at current time
      w = where(loc_maxt_gmem le curTvir.gmem[loc_inds_gmem],count_below)
      coldFrac.gmem_tvircur[i] = float(count_below) / nloc
      
      ; count fraction Tmax below Tvir at accretion time
      w = where(loc_maxt_gmem le accTvir.gmem[loc_inds_gmem],count_below)
      coldFrac.gmem_tviracc[i] = float(count_below) / nloc
    endif
  endfor
  
  ; create composite gal+gmem (Tmax)
  coldFrac.both_tvircur = (coldFrac.gal_tvircur * coldFrac.gal_num + $
                           coldFrac.gmem_tvircur * coldFrac.gmem_num) / $
                          (coldFrac.gal_num + coldFrac.gmem_num)
  coldFrac.both_tviracc = (coldFrac.gal_tviracc * coldFrac.gal_num + $
                           coldFrac.gmem_tviracc * coldFrac.gmem_num) / $
                          (coldFrac.gal_num + coldFrac.gmem_num)
  for j=0,nCuts-1 do $
    coldFrac.both_const[j,*] = (coldFrac.gal_const[j,*] * coldFrac.gal_num + $
                                coldFrac.gmem_const[j,*] * coldFrac.gmem_num) / $
                               (coldFrac.gal_num + coldFrac.gmem_num)  
  
  ; structures to store results (Tcur)
  coldFrac_cur = { gal_const    : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                   gmem_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                   gal_tvircur  : fltarr(n_elements(mt.galcatIDList))         ,$
                   gmem_tvircur : fltarr(n_elements(mt.galcatIDList))         ,$
                   gal_tviracc  : fltarr(n_elements(mt.galcatIDList))         ,$
                   gmem_tviracc : fltarr(n_elements(mt.galcatIDList))         ,$
                   both_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                   both_tvircur : fltarr(n_elements(mt.galcatIDList))         ,$
                   both_tviracc : fltarr(n_elements(mt.galcatIDList))         ,$
                   gal_num      : lonarr(n_elements(mt.galcatIDList))         ,$
                   gmem_num     : lonarr(n_elements(mt.galcatIDList))          }
  
  ; loop over all tracked subgroups (galaxy, Tcur)
  for i=0,n_elements(hist_gal)-1 do begin
    if hist_gal[i] gt 0 then begin
      ; list of indices of galaxy gas particles in this subgroup
      loc_inds_gal   = rev_gal[rev_gal[i]:rev_gal[i+1]-1]
      loc_curt_gal   = curTemp.gal[loc_inds_gal]
      loc_cursfr_gal = curSFR.gal[loc_inds_gal]
      
      ; select only non-eff eos gas
      w_sfr = where(loc_cursfr_gal ne 0.0,count_nosfr)
      if count_nosfr gt 0 then begin
        loc_curt_gal = loc_curt_gal[w_sfr] ; select
        loc_inds_gal = loc_inds_gal[w_sfr]
        
        coldFrac_cur.gal_num[i] = count_nosfr
      
        ; count fraction Tcur below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_curt_gal le TcutVals[j],count_below)
          coldFrac_cur.gal_const[j,i] = float(count_below) / count_nosfr
        endfor
        
        ; count fraction Tmax below Tvir at current time
        w = where(loc_curt_gal le curTvir.gal[loc_inds_gal],count_below)
        coldFrac_cur.gal_tvircur[i] = float(count_below) / count_nosfr
        
        ; count fraction Tmax below Tvir at accretion time
        w = where(loc_curt_gal le accTvir.gal[loc_inds_gal],count_below)
        coldFrac_cur.gal_tviracc[i] = float(count_below) / count_nosfr
      endif ; cursfr!=0
    endif ; hist_gal>0
  endfor  
  
  ; loop over all tracked subgroups (groupmem, Tcur)
  for i=0,n_elements(hist_gmem)-1 do begin
    if hist_gmem[i] gt 0 then begin
      ; list of indices of groupmem gas particles in this subgroup
      loc_inds_gmem   = rev_gmem[rev_gmem[i]:rev_gmem[i+1]-1]
      loc_curt_gmem   = curTemp.gmem[loc_inds_gmem]
      loc_cursfr_gmem = curSFR.gmem[loc_inds_gmem]
      
      ; select only non-eff eos gas
      w_sfr = where(loc_cursfr_gmem ne 0.0,count_nosfr)
      if count_nosfr gt 0 then begin
        loc_curt_gmem = loc_curt_gmem[w_sfr] ; select
        loc_inds_gmem = loc_inds_gmem[w_sfr]
        
        coldFrac_cur.gmem_num[i] = count_nosfr
      
        ; count fraction Tcur below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_curt_gmem le TcutVals[j],count_below)
          coldFrac_cur.gmem_const[j,i] = float(count_below) / count_nosfr
        endfor
        
        ; count fraction Tmax below Tvir at current time
        w = where(loc_curt_gmem le curTvir.gmem[loc_inds_gmem],count_below)
        coldFrac_cur.gmem_tvircur[i] = float(count_below) / count_nosfr
        
        ; count fraction Tmax below Tvir at accretion time
        w = where(loc_curt_gmem le accTvir.gmem[loc_inds_gmem],count_below)
        coldFrac_cur.gmem_tviracc[i] = float(count_below) / count_nosfr
      endif ; cursfr!=0
    endif ; hist_gal>0
  endfor
                        
  ; create composite gal+gmem (Tcur)
  coldFrac_cur.both_tvircur = (coldFrac_cur.gal_tvircur * coldFrac_cur.gal_num + $
                               coldFrac_cur.gmem_tvircur * coldFrac_cur.gmem_num) / $
                              (coldFrac_cur.gal_num + coldFrac_cur.gmem_num)
  coldFrac_cur.both_tviracc = (coldFrac_cur.gal_tviracc * coldFrac_cur.gal_num + $
                               coldFrac_cur.gmem_tviracc * coldFrac_cur.gmem_num) / $
                              (coldFrac_cur.gal_num + coldFrac_cur.gmem_num)
  for j=0,nCuts-1 do $
    coldFrac_cur.both_const[j,*] = (coldFrac_cur.gal_const[j,*] * coldFrac_cur.gal_num + $
                                    coldFrac_cur.gmem_const[j,*] * coldFrac_cur.gmem_num) / $
                                   (coldFrac_cur.gal_num + coldFrac_cur.gmem_num)                       
  
  ; bin fractions into halo mass bins and make median lines
  logMassNbins  = floor((xrange[1]-xrange[0]) / logMassBinSize)
  logMassBins   = linspace(xrange[0],xrange[1],logMassNbins+1) ; edges
  logMassBinCen = linspace(xrange[0],xrange[1],logMassNbins+1) + logMassBinSize/2.0
  logMassBinCen = logMassBinCen[0:-2] ; remove last

  ; combine last two bins (always poor statistics)
  ;logMassNbins  -= 1
  ;logMassBins   = [logMassBins[0:-3],logMassBins[-1]]
  ;logMassBinCen = [logMassBinCen[0:-3],mean(logMassBinCen[-2:-1])]

  medianVals = { const_gal    : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 const_gmem   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 const_both   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 tVircur_gal  : fltarr(logMassNbins)       + !values.f_nan ,$
                 tVircur_gmem : fltarr(logMassNbins)       + !values.f_nan ,$
                 tVircur_both : fltarr(logMassNbins)       + !values.f_nan ,$
                 tViracc_gal  : fltarr(logMassNbins)       + !values.f_nan ,$
                 tViracc_gmem : fltarr(logMassNbins)       + !values.f_nan ,$
                 tViracc_both : fltarr(logMassNbins)       + !values.f_nan  }
                 
  medianVals_cur = { const_gal    : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                     const_gmem   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                     const_both   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                     tVircur_gal  : fltarr(logMassNbins)       + !values.f_nan ,$
                     tVircur_gmem : fltarr(logMassNbins)       + !values.f_nan ,$
                     tVircur_both : fltarr(logMassNbins)       + !values.f_nan ,$
                     tViracc_gal  : fltarr(logMassNbins)       + !values.f_nan ,$
                     tViracc_gmem : fltarr(logMassNbins)       + !values.f_nan ,$
                     tViracc_both : fltarr(logMassNbins)       + !values.f_nan  }
  
  ; calculate median in bins (Tmax) and enforce minimum particle numbers
  for i=0,logMassNbins-2 do begin
    ; gal (Tmax)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac.gal_num ge minNum,count)
    if count gt 0 then begin
      medianVals.tVircur_gal[i] = median(coldFrac.gal_tvircur[w])
      medianVals.tViracc_gal[i] = median(coldFrac.gal_tviracc[w])
      for j=0,nCuts-1 do medianVals.const_gal[j,i] = median(coldFrac.gal_const[j,w])
    endif
    
    ; gmem (Tmax)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac.gmem_num ge minNum,count)
    if count gt 0 then begin
      medianVals.tVircur_gmem[i] = median(coldFrac.gmem_tvircur[w])
      medianVals.tViracc_gmem[i] = median(coldFrac.gmem_tviracc[w])
      for j=0,nCuts-1 do medianVals.const_gmem[j,i] = median(coldFrac.gmem_const[j,w])
    endif
    
    ; both (Tmax)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac.gal_num+coldFrac.gmem_num ge minNum,count)
    if count gt 0 then begin
      medianVals.tVircur_both[i] = median(coldFrac.both_tvircur[w])
      medianVals.tViracc_both[i] = median(coldFrac.both_tviracc[w])
      for j=0,nCuts-1 do medianVals.const_both[j,i] = median(coldFrac.both_const[j,w])
    endif
    
    ; gal (Tcur)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac_cur.gal_num ge minNum,count)
    if count gt 0 then begin
      medianVals_cur.tVircur_gal[i] = median(coldFrac_cur.gal_tvircur[w])
      medianVals_cur.tViracc_gal[i] = median(coldFrac_cur.gal_tviracc[w])
      for j=0,nCuts-1 do medianVals_cur.const_gal[j,i] = median(coldFrac_cur.gal_const[j,w])
    endif
    
    ; gmem (Tcur)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac_cur.gmem_num ge minNum,count)
    if count gt 0 then begin
      medianVals_cur.tVircur_gmem[i] = median(coldFrac_cur.gmem_tvircur[w])
      medianVals_cur.tViracc_gmem[i] = median(coldFrac_cur.gmem_tviracc[w])
      for j=0,nCuts-1 do medianVals_cur.const_gmem[j,i] = median(coldFrac_cur.gmem_const[j,w])
    endif
    
    ; both (Tcur)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac_cur.gal_num+coldFrac_cur.gmem_num ge minNum,count)
    if count gt 0 then begin
      medianVals_cur.tVircur_both[i] = median(coldFrac_cur.both_tvircur[w])
      medianVals_cur.tViracc_both[i] = median(coldFrac_cur.both_tviracc[w])
      for j=0,nCuts-1 do medianVals_cur.const_both[j,i] = median(coldFrac_cur.both_const[j,w])
    endif    
  endfor
  
  ; plot (1) - all data points and median lines (Tmax)
  start_PS, sP.plotPath + 'coldFrac.'+accMode+'.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("log( M_{halo} ) [M_{sun}]"),$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+") "+accMode
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; plot raw points
    psym = 4
    symsize = 0.5
    
    w = where(coldFrac.gal_num+coldFrac.gmem_num ge minNum,count)
    print,count,n_elements(coldFrac.gal_num)
    
    cgPlot,gcMasses,coldFrac.both_tVircur[w],psym=psym,symsize=symsize,thick=1.0,$
      color=getColor(1),/overplot
    cgPlot,gcMasses,coldFrac.both_tViracc[w],psym=psym,symsize=symsize,thick=1.0,$
      color=getColor(2),/overplot
    
    for j=0,nCuts-1 do $
      cgPlot,gcMasses,coldFrac.both_const[j,w],psym=psym,symsize=symsize,thick=1.0,$
        color=getColor(3+j),/overplot
        
    ; plot median lines
    cgPlot,logMassBinCen,medianVals.tVircur_both,color=getColor(1),/overplot
    cgPlot,logMassBinCen,medianVals.tViracc_both,color=getColor(2),/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen,medianVals.const_both[j,*],color=getColor(3+j),/overplot
    
    ; legend
    strings = [textoidl("T_{vir,cur}"),textoidl("T_{vir,acc}"),$
               textoidl("T_{c} = ")+string(TcutVals,format='(f4.1)')]
    colors  = getColor([[1,2],indgen(nCuts)+3],/name)
    legend,strings,textcolor=colors,box=0,/bottom,/left
    
  end_PS
  
  ; plot (2) - just median lines (Tmax)
  start_PS, sP.plotPath + 'coldFrac2.'+accMode+'.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=[9.5,12.25],yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("log( M_{halo} ) [M_{sun}]")
      ;title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+") "+accMode
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; plot raw points
    psym = 4
    symsize = 0.5
    
    w = where(logMassBinCen gt 10.1)
    
    ; gal
    cgPlot,logMassBinCen[w],medianVals.tVircur_gal[w],color=getColor(1),line=0,/overplot
    cgPlot,logMassBinCen[w],medianVals.tViracc_gal[w],color=getColor(2),line=0,/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen[w],medianVals.const_gal[j,w],color=getColor(3+j),line=0,/overplot
      
    ; gmem
    cgPlot,logMassBinCen[w],medianVals.tVircur_gmem[w],color=getColor(1),line=1,/overplot
    cgPlot,logMassBinCen[w],medianVals.tViracc_gmem[w],color=getColor(2),line=1,/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen[w],medianVals.const_gmem[j,w],color=getColor(3+j),line=1,/overplot 
    
    ; both
    ;cgPlot,logMassBinCen,medianVals.tVircur_both,color=getColor(1),/overplot
    ;cgPlot,logMassBinCen,medianVals.tViracc_both,color=getColor(2),/overplot
    ;for j=0,nCuts-1 do $
    ;  cgPlot,logMassBinCen,medianVals.const_both[j,*],color=getColor(3+j),/overplot
    
    ; legend
    strings = [textoidl("T_{vir,cur}"),textoidl("T_{vir,acc}"),$
               textoidl("T_{c} = ")+string(TcutVals,format='(f4.1)')]
    colors  = getColor([[1,2],indgen(nCuts)+3],/name)
    legend,strings,textcolor=colors,box=0,/top,/left
    legend,['central galaxy','halo atmosphere'],linestyle=[0,1],box=0,linesize=0.25,/bottom,/left
    
  end_PS
  
  ; plot (3) - all data points and median lines (Tcur)
  start_PS, sP.plotPath + 'coldFrac_tcur.'+accMode+'.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("log( M_{halo} ) [M_{sun}]")+"",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+") "+accMode
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; plot raw points
    psym = 4
    symsize = 0.5
    
    w = where(coldFrac_cur.gal_num+coldFrac_cur.gmem_num ge minNum,count)
    print,count,n_elements(coldFrac_cur.gal_num)
    
    cgPlot,gcMasses,coldFrac_cur.both_tVircur[w],psym=psym,symsize=symsize,thick=1.0,$
      color=getColor(1),/overplot
    cgPlot,gcMasses,coldFrac_cur.both_tViracc[w],psym=psym,symsize=symsize,thick=1.0,$
      color=getColor(2),/overplot
    
    for j=0,nCuts-1 do $
      cgPlot,gcMasses,coldFrac_cur.both_const[j,w],psym=psym,symsize=symsize,thick=1.0,$
        color=getColor(3+j),/overplot
        
    ; plot median lines
    cgPlot,logMassBinCen,medianVals_cur.tVircur_both,color=getColor(1),/overplot
    cgPlot,logMassBinCen,medianVals_cur.tViracc_both,color=getColor(2),/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen,medianVals_cur.const_both[j,*],color=getColor(3+j),/overplot
    
    ; legend
    strings = [textoidl("T_{vir,cur}"),textoidl("T_{vir,acc}"),$
               textoidl("T_{c} = ")+string(TcutVals,format='(f4.1)')]
    colors  = getColor([[1,2],indgen(nCuts)+3],/name)
    legend,strings,textcolor=colors,box=0,/bottom,/left
    
  end_PS
  
  ; plot (4) - just median lines (Tcur)
  start_PS, sP.plotPath + 'coldFrac_tcur2.'+accMode+'.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("log(M_{halo})")+"",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+") "+accMode
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; plot raw points
    psym = 4
    symsize = 0.5
    
    ; gal
    cgPlot,logMassBinCen,medianVals_cur.tVircur_gal,color=getColor(1),line=1,/overplot
    cgPlot,logMassBinCen,medianVals_cur.tViracc_gal,color=getColor(2),line=1,/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen,medianVals_cur.const_gal[j,*],color=getColor(3+j),line=1,/overplot
      
    ; gmem
    cgPlot,logMassBinCen,medianVals_cur.tVircur_gmem,color=getColor(1),line=2,/overplot
    cgPlot,logMassBinCen,medianVals_cur.tViracc_gmem,color=getColor(2),line=2,/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen,medianVals_cur.const_gmem[j,*],color=getColor(3+j),line=2,/overplot 
    
    ; both
    cgPlot,logMassBinCen,medianVals_cur.tVircur_both,color=getColor(1),/overplot
    cgPlot,logMassBinCen,medianVals_cur.tViracc_both,color=getColor(2),/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen,medianVals_cur.const_both[j,*],color=getColor(3+j),/overplot
    
    ; legend
    strings = [textoidl("T_{vir,cur}"),textoidl("T_{vir,acc}"),$
               textoidl("T_{cut} = ")+string(TcutVals,format='(f4.1)')]
    colors  = getColor([[1,2],indgen(nCuts)+3],/name)
    legend,strings,textcolor=colors,box=0,/bottom,/left
    legend,['gal','gmem','both'],linestyle=[1,2,0],box=0,linesize=0.25,/top,/left
    
  end_PS
end

; plotColdFracVsHaloMassAll(): plot the "cold mode fraction" vs halo mass in a few different ways
;                              without making the mtS/atS cuts

pro plotColdFracVsHaloMassAll, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect = 'pri'
  TcutVals = [5.3,5.4,5.5,5.6,5.7] ; for constant threshold
  nCuts = n_elements(TcutVals)
  
  minNum = 32
  xrange = [9.0,12.5]
  yrange = [0.0,1.1]
  
  logMassBinSize = 0.1
    
  ; load galaxy catalog
  galcat = galaxyCat(sP=sP)
  gcIDList = gcIDList(sP=sP,select=sgSelect)    
    
  ; for all the child gas particles/tracers in the halo selection, replicate parent IDs
  if sP.trMCPerCell eq 0 then begin
    gcIndOrig = galcatRepParentIDs(galcat=galcat,gcIDList=gcIDList)
  endif else begin
    gcIndOrig = galCatRepParentIDs(galcat=galcat,gcIDList=gcIDList,$
                  child_counts={gal:at.child_counts_gal,gmem:at.child_counts_gmem}) 
  endelse
  
  ; compact parents (ascending ID->index)
  placeMap = getIDIndexMap(gcIDList,minid=minid)
  gcIndOrig.gal = placeMap[gcIndOrig.gal-minid]
  gcIndOrig.gmem = placeMap[gcIndOrig.gmem-minid]
  placeMap = !NULL  
  
  ; reverse histogram parent IDs of all particles/tracers in this selection
  hist_gal  = histogram(gcIndOrig.gal,min=0,loc=loc_gal,rev=rev_gal)
  hist_gmem = histogram(gcIndOrig.gmem,min=0,loc=loc_gmem,rev=rev_gmem)

  gcIndOrig = !NULL
  galcat = !NULL
  
  ; load max temps, current tvir, tvir at accretion
  curTvir = gcSubsetProp(sP=sP,select=sgSelect,/virTemp)
  maxTemp = gcSubsetProp(sP=sP,select=sgSelect,/maxPastTemp)

  ; load group cat for subgroup masses
  gc = loadGroupCat(sP=sP,/skipIDs)
  gcMasses = codeMassToLogMsun(gc.subgroupMass[gcIDList])
  gc = !NULL

  ; structures to store results (Tmax)
  coldFrac = { gal_const    : fltarr(nCuts,n_elements(gcIDList))   ,$
               gmem_const   : fltarr(nCuts,n_elements(gcIDList))   ,$
               gal_tvircur  : fltarr(n_elements(gcIDList))         ,$
               gmem_tvircur : fltarr(n_elements(gcIDList))         ,$
               gal_tviracc  : fltarr(n_elements(gcIDList))         ,$
               gmem_tviracc : fltarr(n_elements(gcIDList))         ,$
               both_const   : fltarr(nCuts,n_elements(gcIDList))   ,$
               both_tvircur : fltarr(n_elements(gcIDList))         ,$
               both_tviracc : fltarr(n_elements(gcIDList))         ,$
               gal_num      : lonarr(n_elements(gcIDList))         ,$
               gmem_num     : lonarr(n_elements(gcIDList))          }
  
  ; loop over all tracked subgroups (galaxy, Tmax)
  for i=0,n_elements(hist_gal)-1 do begin
    if hist_gal[i] gt 0 then begin
      ; list of indices of galaxy gas particles in this subgroup
      loc_inds_gal = rev_gal[rev_gal[i]:rev_gal[i+1]-1]
      loc_maxt_gal = maxTemp.gal[loc_inds_gal]
      nloc = n_elements(loc_maxt_gal)
      
      coldFrac.gal_num[i] = nloc
      
      ; count fraction Tmax below each constant temperature threshold
      for j=0,nCuts-1 do begin
        w = where(loc_maxt_gal le TcutVals[j],count_below)
        coldFrac.gal_const[j,i] = float(count_below) / nloc
      endfor
      
      ; count fraction Tmax below Tvir at current time
      w = where(loc_maxt_gal le curTvir.gal[loc_inds_gal],count_below)
      coldFrac.gal_tvircur[i] = float(count_below) / nloc
    endif
  endfor
  
  ; loop over all tracked subgroups (groupmem, Tmax)
  for i=0,n_elements(hist_gmem)-1 do begin
    if hist_gmem[i] gt 0 then begin
      ; list of indices of groupmem gas particles in this subgroup
      loc_inds_gmem = rev_gmem[rev_gmem[i]:rev_gmem[i+1]-1]
      loc_maxt_gmem = maxTemp.gmem[loc_inds_gmem]
      nloc = n_elements(loc_maxt_gmem)
      
      coldFrac.gmem_num[i] = nloc
      
      ; count fraction Tmax below each constant temperature threshold
      for j=0,nCuts-1 do begin
        w = where(loc_maxt_gmem le TcutVals[j],count_below)
        coldFrac.gmem_const[j,i] = float(count_below) / nloc
      endfor
      
      ; count fraction Tmax below Tvir at current time
      w = where(loc_maxt_gmem le curTvir.gmem[loc_inds_gmem],count_below)
      coldFrac.gmem_tvircur[i] = float(count_below) / nloc
    endif
  endfor
  
  ; create composite gal+gmem (Tmax)
  coldFrac.both_tvircur = (coldFrac.gal_tvircur * coldFrac.gal_num + $
                           coldFrac.gmem_tvircur * coldFrac.gmem_num) / $
                          (coldFrac.gal_num + coldFrac.gmem_num)
  coldFrac.both_tviracc = (coldFrac.gal_tviracc * coldFrac.gal_num + $
                           coldFrac.gmem_tviracc * coldFrac.gmem_num) / $
                          (coldFrac.gal_num + coldFrac.gmem_num)
  for j=0,nCuts-1 do $
    coldFrac.both_const[j,*] = (coldFrac.gal_const[j,*] * coldFrac.gal_num + $
                                coldFrac.gmem_const[j,*] * coldFrac.gmem_num) / $
                               (coldFrac.gal_num + coldFrac.gmem_num)             
  
  ; bin fractions into halo mass bins and make median lines
  logMassNbins  = (xrange[1]-xrange[0]) / logMassBinSize
  logMassBins   = linspace(xrange[0],xrange[1],logMassNbins+1) ; edges
  logMassBinCen = linspace(xrange[0],xrange[1],logMassNbins+1) + logMassBinSize/2.0
  
  medianVals = { const_gal    : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 const_gmem   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 const_both   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 tVircur_gal  : fltarr(logMassNbins)       + !values.f_nan ,$
                 tVircur_gmem : fltarr(logMassNbins)       + !values.f_nan ,$
                 tVircur_both : fltarr(logMassNbins)       + !values.f_nan ,$
                 tViracc_gal  : fltarr(logMassNbins)       + !values.f_nan ,$
                 tViracc_gmem : fltarr(logMassNbins)       + !values.f_nan ,$
                 tViracc_both : fltarr(logMassNbins)       + !values.f_nan  }
                 
  ; calculate median in bins (Tmax) and enforce minimum particle numbers
  for i=0,logMassNbins-2 do begin
    ; gal (Tmax)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac.gal_num ge minNum,count)
    if count gt 0 then begin
      medianVals.tVircur_gal[i] = median(coldFrac.gal_tvircur[w])
      medianVals.tViracc_gal[i] = median(coldFrac.gal_tviracc[w])
      for j=0,nCuts-1 do medianVals.const_gal[j,i] = median(coldFrac.gal_const[j,w])
    endif
    
    ; gmem (Tmax)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac.gmem_num ge minNum,count)
    if count gt 0 then begin
      medianVals.tVircur_gmem[i] = median(coldFrac.gmem_tvircur[w])
      medianVals.tViracc_gmem[i] = median(coldFrac.gmem_tviracc[w])
      for j=0,nCuts-1 do medianVals.const_gmem[j,i] = median(coldFrac.gmem_const[j,w])
    endif
    
    ; both (Tmax)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac.gal_num+coldFrac.gmem_num ge minNum,count)
    if count gt 0 then begin
      medianVals.tVircur_both[i] = median(coldFrac.both_tvircur[w])
      medianVals.tViracc_both[i] = median(coldFrac.both_tviracc[w])
      for j=0,nCuts-1 do medianVals.const_both[j,i] = median(coldFrac.both_const[j,w])
    endif
  endfor
  
  ; plot (1) - all data points and median lines (Tmax)
  start_PS, sP.plotPath + 'coldFrac_noMTs.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("log(M_{halo})")+"",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+") noMTs"
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; plot raw points
    psym = 4
    symsize = 0.5
    
    w = where(coldFrac.gal_num+coldFrac.gmem_num ge minNum,count)
    print,count,n_elements(coldFrac.gal_num)
    
    cgPlot,gcMasses,coldFrac.both_tVircur[w],psym=psym,symsize=symsize,thick=1.0,$
      color=getColor(1),/overplot
    cgPlot,gcMasses,coldFrac.both_tViracc[w],psym=psym,symsize=symsize,thick=1.0,$
      color=getColor(2),/overplot
    
    for j=0,nCuts-1 do $
      cgPlot,gcMasses,coldFrac.both_const[j,w],psym=psym,symsize=symsize,thick=1.0,$
        color=getColor(3+j),/overplot
        
    ; plot median lines
    cgPlot,logMassBinCen,medianVals.tVircur_both,color=getColor(1),/overplot
    cgPlot,logMassBinCen,medianVals.tViracc_both,color=getColor(2),/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen,medianVals.const_both[j,*],color=getColor(3+j),/overplot
    
    ; legend
    strings = [textoidl("T_{vir,cur}"),textoidl("T_{vir,acc}"),$
               textoidl("T_{cut} = ")+string(TcutVals,format='(f4.1)')]
    colors  = getColor([[1,2],indgen(nCuts)+3],/name)
    legend,strings,textcolor=colors,box=0,/bottom,/left
    
  end_PS
  
  ; plot (2) - just median lines (Tmax)
  start_PS, sP.plotPath + 'coldFrac_noMTs2.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("log(M_{halo})")+"",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+") noMTs"
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; plot raw points
    psym = 4
    symsize = 0.5
    
    ; gal
    cgPlot,logMassBinCen,medianVals.tVircur_gal,color=getColor(1),line=1,/overplot
    cgPlot,logMassBinCen,medianVals.tViracc_gal,color=getColor(2),line=1,/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen,medianVals.const_gal[j,*],color=getColor(3+j),line=1,/overplot
      
    ; gmem
    cgPlot,logMassBinCen,medianVals.tVircur_gmem,color=getColor(1),line=2,/overplot
    cgPlot,logMassBinCen,medianVals.tViracc_gmem,color=getColor(2),line=2,/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen,medianVals.const_gmem[j,*],color=getColor(3+j),line=2,/overplot 
    
    ; both
    cgPlot,logMassBinCen,medianVals.tVircur_both,color=getColor(1),/overplot
    cgPlot,logMassBinCen,medianVals.tViracc_both,color=getColor(2),/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen,medianVals.const_both[j,*],color=getColor(3+j),/overplot
    
    ; legend
    strings = [textoidl("T_{vir,cur}"),textoidl("T_{vir,acc}"),$
               textoidl("T_{cut} = ")+string(TcutVals,format='(f4.1)')]
    colors  = getColor([[1,2],indgen(nCuts)+3],/name)
    legend,strings,textcolor=colors,box=0,/bottom,/left
    legend,['gal','gmem','both'],linestyle=[1,2,0],box=0,linesize=0.25,/top,/left
    
  end_PS
end