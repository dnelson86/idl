; plotVsHaloMass.pro
; gas accretion project - bin quantities as a function of halo mass
; dnelson may.2012

; haloMassBinAngMom(): bin angular momentum for accreting gas as a function of parent halo mass

function haloMassBinAngMom, sP=sP, sgSelect=sgSelect, accMode=accMode

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  xrange = [9.5,12.5]
  logMassBinSize = 0.25

  ; load
  at = accretionTimes(sP=sP)
  mt = mergerTreeSubset(sP=sP)
  
  ; get indices and mask for either smooth,clumpy,all
  wAm = accModeInds(at=at,accMode=accMode,sP=sP,/mask)
  
  ;gcIndOrig = mergerTreeRepParentIDs(mt=mt,sP=sP,/compactMtS)
  
  ; get maxtemp and tviracc for hot/cold separation
  accTvirAtS = gcSubsetProp(sP=sP,select=sgSelect,/accTvir,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  maxTempAtS = gcSubsetProp(sP=sP,select=sgSelect,/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)

  ; place these in mtS sized arrays
  accTvir = { gal  : fltarr(n_elements(wAm.galMask))  + !values.f_nan ,$
              gmem : fltarr(n_elements(wAm.gmemMask)) + !values.f_nan  }
  maxTemp = { gal  : fltarr(n_elements(wAm.galMask))  + !values.f_nan ,$
              gmem : fltarr(n_elements(wAm.gmemMask)) + !values.f_nan  }
  
  accTvir.gal[wAm.gal] = accTvirAtS.gal
  accTvir.gmem[wAm.gmem] = accTvirAts.gmem
  maxTemp.gal[wAm.gal] = maxTempAtS.gal
  maxTemp.gmem[wAm.gmem] = maxTempAtS.gmem
  
  accTvirAtS = !NULL
  maxTempAtS = !NULL

  ; get parent mass of each gas element (all of mtS since accretion time enforced later)
  parMasses = gcSubsetProp(sP=sP,select=sgSelect,/parMass,/mergerTreeSubset)

  ; load accretionTraj and accretionTrajVel (also all of mtS)
  atj = accretionTraj(sP=sP)
  atv = accretionTraj(sP=sP,/getVel)

  ; convert snapshot numbers to scale factors
  snapTimes = snapNumToRedshift(snap=0,/time,/all,sP=sP)
  snapTimes = snapTimes[where(snapTimes ge 0)] ; remove -1
    
  if n_elements(parMasses.gal) ne n_elements(wAm.galMask) or $
     n_elements(atj.relpos_gal[0,0,*]) ne n_elements(wAm.galMask) or $
     n_elements(atv.vel_gal[0,0,*]) ne n_elements(wAm.galMask) then message,'Error: Data error.'

  ; for each radial point, find the crossing snapshot for each gas element and calculate
  ; specific angular momentum at that time using accretionTraj()
  for k=0,n_elements(at.rVirFacs)-1 do begin
  
    gal_w  = where(wAm.galMask eq 1B and at.accTime_gal[k,*] ne -1,count_gal)
    gmem_w = where(wAm.gmemMask eq 1B and at.accTime_gmem[k,*] ne -1,count_gmem)

    if count_gal gt 0 then begin
      ; find accretion snapshots and double check
      accSnap = value_locate(snapTimes,at.accTime_gal[k,gal_w])
      w = where(accSnap le 0 or accSnap gt sP.snapRange[1],count)
      if count gt 0 then message,'Error: Bad snapshot time mapping gal two.'    
      
      ; calculate angular momentum
      radxyz = fltarr(3,count_gal)
      velxyz = fltarr(3,count_gal)
      
      for i=0UL,count_gal-1 do begin
        radxyz[*,i] = atj.relPos_gal[mt.maxSnap-accSnap[i],*,gal_w[i]]
        velxyz[*,i] = atv.vel_gal[mt.maxSnap-accSnap[i],*,gal_w[i]]
      endfor
      print,mean(radxyz),mean(velxyz)
      
      jvec = fltarr(3,count_gal)
      jvec[0,*] = radxyz[1,*] * velxyz[2,*] - radxyz[2,*] * velxyz[1,*]
      jvec[1,*] = radxyz[2,*] * velxyz[0,*] - radxyz[0,*] * velxyz[2,*]
      jvec[2,*] = radxyz[0,*] * velxyz[1,*] - radxyz[1,*] * velxyz[0,*]
      
      jnorm = reform(sqrt(jvec[0,*]*jvec[0,*] + jvec[1,*]*jvec[1,*] + jvec[2,*]*jvec[2,*]))
      
      ; normalize by j_circ for parent halos given their mass and rvir crossing times
      ;radnorm = reform(sqrt(radxyz[0,*]*radxyz[0,*] + radxyz[1,*]*radxyz[1,*] + radxyz[2,*]*radxyz[2,*]))
      ;velnorm = reform(sqrt(velxyz[0,*]*velxyz[0,*] + velxyz[1,*]*velxyz[1,*] + velxyz[2,*]*velxyz[2,*]))
      ;jcirc = radnorm * velnorm
      
      at.accTime_gal[k,gal_w] = jnorm ;/ jcirc
    endif
    
    if count_gmem gt 0 then begin
      ; repeat for gmem
      accSnap = value_locate(snapTimes,at.accTime_gmem[k,gmem_w])
      w = where(accSnap le 0 or accSnap gt sP.snapRange[1],count)
      if count gt 0 then message,'Error: Bad snapshot time mapping gal two.'    
      
      ; calculate angular momentum
      radxyz = fltarr(3,count_gmem)
      velxyz = fltarr(3,count_gmem)
      
      for i=0UL,count_gmem-1 do begin
        radxyz[*,i] = atj.relPos_gmem[mt.maxSnap-accSnap[i],*,gmem_w[i]]
        velxyz[*,i] = atv.vel_gmem[mt.maxSnap-accSnap[i],*,gmem_w[i]]
      endfor
      
      jvec = fltarr(3,count_gmem)
      jvec[0,*] = radxyz[1,*] * velxyz[2,*] - radxyz[2,*] * velxyz[1,*]
      jvec[1,*] = radxyz[2,*] * velxyz[0,*] - radxyz[0,*] * velxyz[2,*]
      jvec[2,*] = radxyz[0,*] * velxyz[1,*] - radxyz[1,*] * velxyz[0,*]
      
      jnorm = reform(sqrt(jvec[0,*]*jvec[0,*] + jvec[1,*]*jvec[1,*] + jvec[2,*]*jvec[2,*]))
      
      ; normalize by j_circ for parent halos given their mass and rvir crossing times
      ;radnorm = reform(sqrt(radxyz[0,*]*radxyz[0,*] + radxyz[1,*]*radxyz[1,*] + radxyz[2,*]*radxyz[2,*]))
      ;velnorm = reform(sqrt(velxyz[0,*]*velxyz[0,*] + velxyz[1,*]*velxyz[1,*] + velxyz[2,*]*velxyz[2,*]))
      ;jcirc = radnorm * velnorm
      
      at.accTime_gmem[k,gmem_w] = jnorm ;/ jcirc
    endif
  endfor

  ; bin fractions into halo mass bins and make median lines
  ;logMassNbins  = floor((xrange[1]-xrange[0]) / logMassBinSize)
  ;logMassBins   = linspace(xrange[0],xrange[1],logMassNbins+1) ; edges
  ;logMassBinCen = linspace(xrange[0],xrange[1],logMassNbins+1) + logMassBinSize/2.0
  ;logMassBinCen = logMassBinCen[0:-2] ; remove last
  
  ; manual
  logMassBins=[9.5,10.0,10.15,10.3,10.45,10.6,10.75,10.9,11.0,$
               11.1,11.25,11.5,11.75,12.0,12.25,13.0]
  logMassNBins = n_elements(logMassBins)-1
  logMassBinCen = 0.5 * (logMassBins + shift(logMassBins,-1))
  logMassBinCen = logMassBinCen[0:-2]

  binnedVals = { hotMode : { $
                   mean_gal    : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   mean_gmem   : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   mean_both   : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   median_gal  : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   median_gmem : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   median_both : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   stddev_gal  : fltarr(n_elements(at.rVirFacs),logMassNBins) + !values.f_nan ,$
                   stddev_gmem : fltarr(n_elements(at.rVirFacs),logMassNBins) + !values.f_nan ,$
                   stddev_both : fltarr(n_elements(at.rVirFacs),logMassNBins) + !values.f_nan $
                 } ,$
                 coldMode : { $
                   mean_gal    : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   mean_gmem   : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   mean_both   : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   median_gal  : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   median_gmem : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   median_both : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   stddev_gal  : fltarr(n_elements(at.rVirFacs),logMassNBins) + !values.f_nan ,$
                   stddev_gmem : fltarr(n_elements(at.rVirFacs),logMassNBins) + !values.f_nan ,$
                   stddev_both : fltarr(n_elements(at.rVirFacs),logMassNBins) + !values.f_nan $
                 } ,$
                 logMassBinCen : logMassBinCen ,$
                 rVirFacs      : at.rVirFacs       }
  
  ; calculate median in bins (Tmax) and enforce minimum particle numbers
  for i=0,logMassNbins-1 do begin
    for k=0,n_elements(at.rVirFacs)-1 do begin
      ; hot gal
      w_gal = where(parMasses.gal gt logMassBins[i] and parMasses.gal le logMassBins[i+1] and $
                    at.accTime_gal[k,*] ne -1 and $
                    maxTemp.gal gt accTvir.gal and wAm.galMask eq 1B,count_gal)   
      if count_gal gt 0 then begin
        binnedVals.hotMode.mean_gal[k,i]   = mean(at.accTime_gal[k,w_gal])
        binnedVals.hotMode.median_gal[k,i] = median(at.accTime_gal[k,w_gal])
        binnedVals.hotMode.stddev_gal[k,i] = stddev(at.accTime_gal[k,w_gal])
      endif
      
      ; hot gmem
      w_gmem = where(parMasses.gmem gt logMassBins[i] and parMasses.gmem le logMassBins[i+1] and $
                     at.accTime_gmem[k,*] ne -1 and $
                     maxTemp.gmem gt accTvir.gmem and wAm.gmemMask eq 1B,count_gmem)
      if count_gmem gt 0 then begin
        binnedVals.hotMode.mean_gmem[k,i]   = mean(at.accTime_gmem[k,w_gmem])
        binnedVals.hotMode.median_gmem[k,i] = median(at.accTime_gmem[k,w_gmem])
        binnedVals.hotMode.stddev_gmem[k,i] = stddev(at.accTime_gmem[k,w_gmem])
      endif
      
      ; hot both
      vals = []
      if count_gal gt 0  then vals = [vals,reform(at.accTime_gal[k,w_gal])]
      if count_gmem gt 0 then vals = [vals,reform(at.accTime_gmem[k,w_gmem])]
      
      if (count_gal+count_gmem) gt 0 then begin
        binnedVals.hotMode.mean_both[k,i]   = mean(vals)
        binnedVals.hotMode.median_both[k,i] = median(vals)
        binnedVals.hotMode.stddev_both[k,i] = stddev(vals)
      endif
      
      ; cold gal
      w_gal = where(parMasses.gal gt logMassBins[i] and parMasses.gal le logMassBins[i+1] and $
                    at.accTime_gal[k,*] ne -1 and $
                    maxTemp.gal le accTvir.gal and wAm.galMask eq 1B,count_gal)   
      if count_gal gt 0 then begin
        binnedVals.coldMode.mean_gal[k,i]   = mean(at.accTime_gal[k,w_gal])
        binnedVals.coldMode.median_gal[k,i] = median(at.accTime_gal[k,w_gal])
        binnedVals.coldMode.stddev_gal[k,i] = stddev(at.accTime_gal[k,w_gal])
      endif
      
      ; cold gmem
      w_gmem = where(parMasses.gmem gt logMassBins[i] and parMasses.gmem le logMassBins[i+1] and $
                     at.accTime_gmem[k,*] ne -1 and $
                     maxTemp.gmem le accTvir.gmem and wAm.gmemMask eq 1B,count_gmem)
      if count_gmem gt 0 then begin
        binnedVals.coldMode.mean_gmem[k,i]   = mean(at.accTime_gmem[k,w_gmem])
        binnedVals.coldMode.median_gmem[k,i] = median(at.accTime_gmem[k,w_gmem])
        binnedVals.coldMode.stddev_gmem[k,i] = stddev(at.accTime_gmem[k,w_gmem])
      endif
      
      ; cold both
      vals = []
      if count_gal gt 0  then vals = [vals,reform(at.accTime_gal[k,w_gal])]
      if count_gmem gt 0 then vals = [vals,reform(at.accTime_gmem[k,w_gmem])]
      
      if (count_gal+count_gmem) gt 0 then begin
        binnedVals.coldMode.mean_both[k,i]   = mean(vals)
        binnedVals.coldMode.median_both[k,i] = median(vals)
        binnedVals.coldMode.stddev_both[k,i] = stddev(vals)
      endif
    endfor ;k
  endfor ;i

  ; debug plot
  start_PS,sP.plotPath + 'angmom.histos.all.'+accMode+'.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    binsize = 0.1 & strings = []
    cgPlot,[0],[0],/nodata,xrange=[1,6],yrange=[0.001,1.0],/ylog,/ys,/xs,$
      xtitle=textoidl("log ( j ) [kpc km/s]"),ytitle="Fraction"
      
    for k=0,n_elements(at.rVirFacs)-1 do begin
      w = where(at.accTime_gal[k,*] ne -1 and wAm.galMask eq 1B,count)
      at.accTime_gal[k,w] = alog10(at.accTime_gal[k,w])
      if count gt 1 then begin
        hist = histogram(at.accTime_gal[k,w],binsize=binsize,loc=loc)
        cgplot,loc+binsize*0.5,hist/float(total(hist)),line=0,color=getColor(k),/overplot
        cgplot,[mean(at.accTime_gal[k,w]),mean(at.accTime_gal[k,w])],[0.8,1.0],line=0,color=getColor(k),/overplot
        cgplot,[median(at.accTime_gal[k,w]),median(at.accTime_gal[k,w])],[0.5,0.7],line=0,color=getColor(k),/overplot
      endif
      at.accTime_gal[k,w] = 10.0^(at.accTime_gal[k,w])
      
      w = where(at.accTime_gmem[k,*] ne -1 and wAm.gmemMask eq 1B,count)
      at.accTime_gmem[k,w] = alog10(at.accTime_gmem[k,w])
      if count gt 1 then begin
        hist = histogram(at.accTime_gmem[k,w],binsize=binsize,loc=loc)
        cgplot,loc+binsize*0.5,hist/float(total(hist)),line=1,color=getColor(k),/overplot
        cgplot,[mean(at.accTime_gmem[k,w]),mean(at.accTime_gmem[k,w])],[0.8,1.0],line=1,color=getColor(k),/overplot
        cgplot,[median(at.accTime_gmem[k,w]),median(at.accTime_gmem[k,w])],[0.5,0.7],line=1,color=getColor(k),/overplot
      endif
      strings = [strings,string(at.rVirFacs[k],format='(f4.2)')]
      at.accTime_gmem[k,w] = 10.0^(at.accTime_gmem[k,w])
    endfor
    
    legend,strings,textcolors=getColor([0,1,2,3,4,5,6],/name),box=0,/top,/left
  end_PS
  
  start_PS,sP.plotPath + 'angmom.histos.cold.'+accMode+'.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    binsize = 0.1 & strings = []
    cgPlot,[0],[0],/nodata,xrange=[1,6],yrange=[0.001,1.0],/ylog,/ys,/xs,$
      xtitle=textoidl("log ( j ) [kpc km/s]"),ytitle="Fraction"
      
    for k=0,n_elements(at.rVirFacs)-1 do begin
      w = where(at.accTime_gal[k,*] ne -1 and wAm.galMask eq 1B and maxTemp.gal le accTvir.gal,count)
      at.accTime_gal[k,w] = alog10(at.accTime_gal[k,w])
      if count gt 1 then begin
        hist = histogram(at.accTime_gal[k,w],binsize=binsize,loc=loc)
        cgplot,loc+binsize*0.5,hist/float(total(hist)),line=0,color=getColor(k),/overplot
        cgplot,[mean(at.accTime_gal[k,w]),mean(at.accTime_gal[k,w])],[0.8,1.0],line=0,color=getColor(k),/overplot
        cgplot,[median(at.accTime_gal[k,w]),median(at.accTime_gal[k,w])],[0.5,0.7],line=0,color=getColor(k),/overplot
      endif
      at.accTime_gal[k,w] = 10.0^(at.accTime_gal[k,w])
      
      w = where(at.accTime_gmem[k,*] ne -1 and wAm.gmemMask eq 1B and maxTemp.gmem le accTvir.gmem,count)
      at.accTime_gmem[k,w] = alog10(at.accTime_gmem[k,w])
      if count gt 1 then begin
        hist = histogram(at.accTime_gmem[k,w],binsize=binsize,loc=loc)
        cgplot,loc+binsize*0.5,hist/float(total(hist)),line=1,color=getColor(k),/overplot
        cgplot,[mean(at.accTime_gmem[k,w]),mean(at.accTime_gmem[k,w])],[0.8,1.0],line=1,color=getColor(k),/overplot
        cgplot,[median(at.accTime_gmem[k,w]),median(at.accTime_gmem[k,w])],[0.5,0.7],line=1,color=getColor(k),/overplot
      endif
      strings = [strings,string(at.rVirFacs[k],format='(f4.2)')]
      at.accTime_gmem[k,w] = 10.0^(at.accTime_gmem[k,w])
    endfor
    
    legend,strings,textcolors=getColor([0,1,2,3,4,5,6],/name),box=0,/top,/left
  end_PS
  
  start_PS,sP.plotPath + 'angmom.histos.hot.'+accMode+'.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    binsize = 0.1 & strings = []
    cgPlot,[0],[0],/nodata,xrange=[1,6],yrange=[0.001,1.0],/ylog,/ys,/xs,$
      xtitle=textoidl("log ( j ) [kpc km/s]"),ytitle="Fraction"
      
    for k=0,n_elements(at.rVirFacs)-1 do begin
      w = where(at.accTime_gal[k,*] ne -1 and wAm.galMask eq 1B and maxTemp.gal gt accTvir.gal,count)
      at.accTime_gal[k,w] = alog10(at.accTime_gal[k,w])
      if count gt 1 then begin
        hist = histogram(at.accTime_gal[k,w],binsize=binsize,loc=loc)
        cgplot,loc+binsize*0.5,hist/float(total(hist)),line=0,color=getColor(k),/overplot
        cgplot,[mean(at.accTime_gal[k,w]),mean(at.accTime_gal[k,w])],[0.8,1.0],line=0,color=getColor(k),/overplot
        cgplot,[median(at.accTime_gal[k,w]),median(at.accTime_gal[k,w])],[0.5,0.7],line=0,color=getColor(k),/overplot
      endif
      at.accTime_gal[k,w] = 10.0^(at.accTime_gal[k,w])
      
      w = where(at.accTime_gmem[k,*] ne -1 and wAm.gmemMask eq 1B and maxTemp.gmem gt accTvir.gmem,count)
      at.accTime_gmem[k,w] = alog10(at.accTime_gmem[k,w])
      if count gt 1 then begin
        hist = histogram(at.accTime_gmem[k,w],binsize=binsize,loc=loc)
        cgplot,loc+binsize*0.5,hist/float(total(hist)),line=1,color=getColor(k),/overplot
        cgplot,[mean(at.accTime_gmem[k,w]),mean(at.accTime_gmem[k,w])],[0.8,1.0],line=1,color=getColor(k),/overplot
        cgplot,[median(at.accTime_gmem[k,w]),median(at.accTime_gmem[k,w])],[0.5,0.7],line=1,color=getColor(k),/overplot
      endif
      strings = [strings,string(at.rVirFacs[k],format='(f4.2)')]
      at.accTime_gmem[k,w] = 10.0^(at.accTime_gmem[k,w])
    endfor
    
    legend,strings,textcolors=getColor([0,1,2,3,4,5,6],/name),box=0,/top,/left
  end_PS
  
  return,binnedVals
end

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
  
  gcIndOrig = mergerTreeRepParentIDs(mt=mt,sP=sP,/compactMtS)
  
  ; get maxtemp and tviracc for hot/cold separation
  accTvirAtS = gcSubsetProp(sP=sP,select=sgSelect,/accTvir,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  maxTempAtS = gcSubsetProp(sP=sP,select=sgSelect,/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)

  ; place these in mtS sized arrays
  accTvir = { gal  : fltarr(n_elements(gcIndOrig.gal))  + !values.f_nan ,$
              gmem : fltarr(n_elements(gcIndOrig.gmem)) + !values.f_nan  }
  maxTemp = { gal  : fltarr(n_elements(gcIndOrig.gal))  + !values.f_nan ,$
              gmem : fltarr(n_elements(gcIndOrig.gmem)) + !values.f_nan  }
  
  accTvir.gal[wAm.gal] = accTvirAtS.gal
  accTvir.gmem[wAm.gmem] = accTvirAts.gmem
  maxTemp.gal[wAm.gal] = maxTempAtS.gal
  maxTemp.gmem[wAm.gmem] = maxTempAtS.gmem
  
  accTvirAtS = !NULL
  maxTempAtS = !NULL

  ; get parent mass of each gas element (all of mtS since accretion time enforced later)
  parMasses = gcSubsetProp(sP=sP,select=sgSelect,/parMass,/mergerTreeSubset)

  ; convert (k=0) r=rvir accretion times to nearest snap
  snapTimes = snapNumToRedshift(snap=0,/time,/all,sP=sP)
  snapTimes = snapTimes[where(snapTimes ge 0)] ; remove -1
    
  accSnap = {gal  : intarr(n_elements(at.accTime_gal[0,*]))   ,$
             gmem : intarr(n_elements(at.accTime_gmem[0,*]))   }
             
  accSnap.gal[wAm.gal]   = value_locate(snapTimes,at.accTime_gal[0,wAm.gal]) 
  accSnap.gmem[wAm.gmem] = value_locate(snapTimes,at.accTime_gmem[0,wAm.gmem])

  if n_elements(gcIndOrig.gal) ne n_elements(parMasses.gal) or $
     n_elements(gcIndOrig.gal) ne n_elements(wAm.galMask) then message,'Error: Data error.'
  w = where(accSnap.gal[wAm.gal] le 0 or accSnap.gal[wAm.gal] gt sP.snapRange[1],count1)
  w = where(accSnap.gmem[wAm.gmem] le 0 or accSnap.gmem[wAm.gmem] gt sP.snapRange[1],count2)
  if count1 gt 0 or count2 gt 0 then message,'Error: Bad snapshot time mapping.'

  ; each accDt now stores the scale factor when that particle crossed the radius of interest
  ; convert to delta time (Myr) since crossing rvir
  for k=1,n_elements(at.rVirFacs)-1 do begin
    gal_w  = where(wAm.galMask eq 1B and at.accTime_gal[0,*] ne -1 and at.accTime_gal[k,*] ne -1,$
                   count_gal,comp=gal_wc,ncomp=countc_gal)
    gmem_w = where(wAm.gmemMask eq 1B and at.accTime_gmem[0,*] ne -1 and at.accTime_gmem[k,*] ne -1,$
                   count_gmem,comp=gmem_wc,ncomp=countc_gmem)

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
  ;logMassNbins  = floor((xrange[1]-xrange[0]) / logMassBinSize)
  ;logMassBins   = linspace(xrange[0],xrange[1],logMassNbins+1) ; edges
  ;logMassBinCen = linspace(xrange[0],xrange[1],logMassNbins+1) + logMassBinSize/2.0
  ;logMassBinCen = logMassBinCen[0:-2] ; remove last
  
  ; manual
  logMassBins=[9.5,10.0,10.15,10.3,10.45,10.6,10.75,10.9,11.0,$
               11.1,11.25,11.5,11.75,12.0,12.25,13.0]
  logMassNBins = n_elements(logMassBins)-1
  logMassBinCen = 0.5 * (logMassBins + shift(logMassBins,-1))
  logMassBinCen = logMassBinCen[0:-2]

  binnedVals = { hotMode : { $
                   mean_gal    : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   mean_gmem   : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   mean_both   : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   median_gal  : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   median_gmem : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   median_both : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   stddev_gal  : fltarr(n_elements(at.rVirFacs),logMassNBins) + !values.f_nan ,$
                   stddev_gmem : fltarr(n_elements(at.rVirFacs),logMassNBins) + !values.f_nan ,$
                   stddev_both : fltarr(n_elements(at.rVirFacs),logMassNBins) + !values.f_nan $
                 } ,$
                 coldMode : { $
                   mean_gal    : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   mean_gmem   : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   mean_both   : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   median_gal  : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   median_gmem : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   median_both : fltarr(n_elements(at.rVirFacs),logMassNbins) + !values.f_nan ,$
                   stddev_gal  : fltarr(n_elements(at.rVirFacs),logMassNBins) + !values.f_nan ,$
                   stddev_gmem : fltarr(n_elements(at.rVirFacs),logMassNBins) + !values.f_nan ,$
                   stddev_both : fltarr(n_elements(at.rVirFacs),logMassNBins) + !values.f_nan $
                 } ,$
                 logMassBinCen : logMassBinCen ,$
                 rVirFacs      : at.rVirFacs       }
  
  ; calculate median in bins (Tmax) and enforce minimum particle numbers
  for i=0,logMassNbins-1 do begin
    for k=1,n_elements(at.rVirFacs)-1 do begin
      ; hot gal
      w_gal = where(parMasses.gal gt logMassBins[i] and parMasses.gal le logMassBins[i+1] and $
                    at.accTime_gal[0,*] ne -1 and at.accTime_gal[k,*] ne -1 and $
                    maxTemp.gal gt accTvir.gal and wAm.galMask eq 1B,count_gal)   
      if count_gal gt 0 then begin
        binnedVals.hotMode.mean_gal[k,i]   = mean(at.accTime_gal[k,w_gal])
        binnedVals.hotMode.median_gal[k,i] = median(at.accTime_gal[k,w_gal])
        binnedVals.hotMode.stddev_gal[k,i] = stddev(at.accTime_gal[k,w_gal])
      endif
      
      ; hot gmem
      w_gmem = where(parMasses.gmem gt logMassBins[i] and parMasses.gmem le logMassBins[i+1] and $
                     at.accTime_gmem[0,*] ne -1 and at.accTime_gmem[k,*] ne -1 and $
                     maxTemp.gmem gt accTvir.gmem and wAm.gmemMask eq 1B,count_gmem)
      if count_gmem gt 0 then begin
        binnedVals.hotMode.mean_gmem[k,i]   = mean(at.accTime_gmem[k,w_gmem])
        binnedVals.hotMode.median_gmem[k,i] = median(at.accTime_gmem[k,w_gmem])
        binnedVals.hotMode.stddev_gmem[k,i] = stddev(at.accTime_gmem[k,w_gmem])
      endif
      
      ; hot both
      vals = []
      if count_gal gt 0  then vals = [vals,reform(at.accTime_gal[k,w_gal])]
      if count_gmem gt 0 then vals = [vals,reform(at.accTime_gmem[k,w_gmem])]
      
      if (count_gal+count_gmem) gt 0 then begin
        binnedVals.hotMode.mean_both[k,i]   = mean(vals)
        binnedVals.hotMode.median_both[k,i] = median(vals)
        binnedVals.hotMode.stddev_both[k,i] = stddev(vals)
      endif
      
      ; cold gal
      w_gal = where(parMasses.gal gt logMassBins[i] and parMasses.gal le logMassBins[i+1] and $
                    at.accTime_gal[0,*] ne -1 and at.accTime_gal[k,*] ne -1 and $
                    maxTemp.gal le accTvir.gal and wAm.galMask eq 1B,count_gal)   
      if count_gal gt 0 then begin
        binnedVals.coldMode.mean_gal[k,i]   = mean(at.accTime_gal[k,w_gal])
        binnedVals.coldMode.median_gal[k,i] = median(at.accTime_gal[k,w_gal])
        binnedVals.coldMode.stddev_gal[k,i] = stddev(at.accTime_gal[k,w_gal])
      endif
      
      ; cold gmem
      w_gmem = where(parMasses.gmem gt logMassBins[i] and parMasses.gmem le logMassBins[i+1] and $
                     at.accTime_gmem[0,*] ne -1 and at.accTime_gmem[k,*] ne -1 and $
                     maxTemp.gmem le accTvir.gmem and wAm.gmemMask eq 1B,count_gmem)
      if count_gmem gt 0 then begin
        binnedVals.coldMode.mean_gmem[k,i]   = mean(at.accTime_gmem[k,w_gmem])
        binnedVals.coldMode.median_gmem[k,i] = median(at.accTime_gmem[k,w_gmem])
        binnedVals.coldMode.stddev_gmem[k,i] = stddev(at.accTime_gmem[k,w_gmem])
      endif
      
      ; cold both
      vals = []
      if count_gal gt 0  then vals = [vals,reform(at.accTime_gal[k,w_gal])]
      if count_gmem gt 0 then vals = [vals,reform(at.accTime_gmem[k,w_gmem])]
      
      if (count_gal+count_gmem) gt 0 then begin
        binnedVals.coldMode.mean_both[k,i]   = mean(vals)
        binnedVals.coldMode.median_both[k,i] = median(vals)
        binnedVals.coldMode.stddev_both[k,i] = stddev(vals)
      endif
    endfor ;k
  endfor ;i

  ; debug plot
  start_PS,sP.plotPath + 'accdt.histos.all.'+accMode+'.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    binsize = 0.01 & strings = []
    cgPlot,[0],[0],/nodata,xrange=[0.0,0.25],yrange=[0.0001,1.1],/ylog,/ys,/xs,$
      xtitle=textoidl("\Delta t_{acc} / \tau_{circ}"),ytitle="Fraction"
      
    for k=1,n_elements(at.rVirFacs)-1 do begin
      w = where(at.accTime_gal[0,*] ne -1 and at.accTime_gal[k,*] ne -1 and wAm.galMask eq 1B,count)
      if count gt 1 then begin
        hist = histogram(at.accTime_gal[k,w],binsize=binsize,loc=loc)
        cgplot,loc+binsize*0.5,hist/float(total(hist)),line=0,color=getColor(k),/overplot
        cgplot,[mean(at.accTime_gal[k,w]),mean(at.accTime_gal[k,w])],[0.8,1.0],line=0,color=getColor(k),/overplot
        cgplot,[median(at.accTime_gal[k,w]),median(at.accTime_gal[k,w])],[0.5,0.7],line=0,color=getColor(k),/overplot
      endif
      w = where(at.accTime_gmem[0,*] ne -1 and at.accTime_gmem[k,*] ne -1 and wAm.gmemMask eq 1B,count)
      if count gt 1 then begin
        hist = histogram(at.accTime_gmem[k,w],binsize=binsize,loc=loc)
        cgplot,loc+binsize*0.5,hist/float(total(hist)),line=1,color=getColor(k),/overplot
        cgplot,[mean(at.accTime_gmem[k,w]),mean(at.accTime_gmem[k,w])],[0.8,1.0],line=1,color=getColor(k),/overplot
        cgplot,[median(at.accTime_gmem[k,w]),median(at.accTime_gmem[k,w])],[0.5,0.7],line=1,color=getColor(k),/overplot
      endif
      strings = [strings,string(at.rVirFacs[k],format='(f4.2)')]
    endfor
    
    legend,strings,textcolors=getColor([1,2,3,4,5,6],/name),box=0,/top,/right
  end_PS
  
  start_PS,sP.plotPath + 'accdt.histos.cold.'+accMode+'.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    binsize = 0.01 & strings = []
    cgPlot,[0],[0],/nodata,xrange=[0.0,0.25],yrange=[0.0001,1.1],/ylog,/ys,/xs,$
      xtitle=textoidl("\Delta t_{acc} / \tau_{circ}"),ytitle="Fraction"
      
    for k=1,n_elements(at.rVirFacs)-1 do begin
      w = where(at.accTime_gal[0,*] ne -1 and at.accTime_gal[k,*] ne -1 and wAm.galMask eq 1B and $
                maxTemp.gal le accTvir.gal,count)
      if count gt 1 then begin
        hist = histogram(at.accTime_gal[k,w],binsize=binsize,loc=loc)
        cgplot,loc+binsize*0.5,hist/float(total(hist)),line=0,color=getColor(k),/overplot
        cgplot,[mean(at.accTime_gal[k,w]),mean(at.accTime_gal[k,w])],[0.8,1.0],line=0,color=getColor(k),/overplot
        cgplot,[median(at.accTime_gal[k,w]),median(at.accTime_gal[k,w])],[0.5,0.7],line=0,color=getColor(k),/overplot
      endif
      w = where(at.accTime_gmem[0,*] ne -1 and at.accTime_gmem[k,*] ne -1 and wAm.gmemMask eq 1B,count)
      if count gt 1 then begin
        hist = histogram(at.accTime_gmem[k,w],binsize=binsize,loc=loc)
        cgplot,loc+binsize*0.5,hist/float(total(hist)),line=1,color=getColor(k),/overplot
        cgplot,[mean(at.accTime_gmem[k,w]),mean(at.accTime_gmem[k,w])],[0.8,1.0],line=1,color=getColor(k),/overplot
        cgplot,[median(at.accTime_gmem[k,w]),median(at.accTime_gmem[k,w])],[0.5,0.7],line=1,color=getColor(k),/overplot
      endif
      strings = [strings,string(at.rVirFacs[k],format='(f4.2)')]
    endfor
    
    legend,strings,textcolors=getColor([1,2,3,4,5,6],/name),box=0,/top,/right
  end_PS
  
  start_PS,sP.plotPath + 'accdt.histos.hot.'+accMode+'.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    binsize = 0.01 & strings = []
    cgPlot,[0],[0],/nodata,xrange=[0.0,0.25],yrange=[0.0001,1.1],/ylog,/ys,/xs,$
      xtitle=textoidl("\Delta t_{acc} / \tau_{circ}"),ytitle="Fraction"
      
    for k=1,n_elements(at.rVirFacs)-1 do begin
      w = where(at.accTime_gal[0,*] ne -1 and at.accTime_gal[k,*] ne -1 and wAm.galMask eq 1B and $
                maxTemp.gal gt accTvir.gal,count)
      if count gt 1 then begin
        hist = histogram(at.accTime_gal[k,w],binsize=binsize,loc=loc)
        cgplot,loc+binsize*0.5,hist/float(total(hist)),line=0,color=getColor(k),/overplot
        cgplot,[mean(at.accTime_gal[k,w]),mean(at.accTime_gal[k,w])],[0.8,1.0],line=0,color=getColor(k),/overplot
        cgplot,[median(at.accTime_gal[k,w]),median(at.accTime_gal[k,w])],[0.5,0.7],line=0,color=getColor(k),/overplot
      endif
      w = where(at.accTime_gmem[0,*] ne -1 and at.accTime_gmem[k,*] ne -1 and wAm.gmemMask eq 1B,count)
      if count gt 1 then begin
        hist = histogram(at.accTime_gmem[k,w],binsize=binsize,loc=loc)
        cgplot,loc+binsize*0.5,hist/float(total(hist)),line=1,color=getColor(k),/overplot
        cgplot,[mean(at.accTime_gmem[k,w]),mean(at.accTime_gmem[k,w])],[0.8,1.0],line=1,color=getColor(k),/overplot
        cgplot,[median(at.accTime_gmem[k,w]),median(at.accTime_gmem[k,w])],[0.5,0.7],line=1,color=getColor(k),/overplot
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
  
  logMassBinSize = 0.2 ;/ (sP.res/128)
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
  ;logMassNbins  = floor((xrange[1]-xrange[0]) / logMassBinSize)
  ;logMassBins   = linspace(xrange[0],xrange[1],logMassNbins+1) ; edges
  ;logMassBinCen = linspace(xrange[0],xrange[1],logMassNbins+1) + logMassBinSize/2.0
  ;logMassBinCen = logMassBinCen[0:-2] ; remove last
  
  ; manual
  logMassBins=[9.5,10.0,10.15,10.3,10.45,10.6,10.75,10.9,11.0,$
               11.1,11.25,11.5,11.75,12.0,12.25,13.0]
  logMassNBins = n_elements(logMassBins)-1
  logMassBinCen = 0.5 * (logMassBins + shift(logMassBins,-1))
  logMassBinCen = logMassBinCen[0:-2]
  
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
        coldMass.const_gal[j,i]   = count_below * massGasPart
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
  coldMass.tViracc_gal[*]  = codeMassToLogMsun(coldMass.tVircur_gal[*])  ;/ logMassBinSize
  hotMass.tViracc_gal[*]   = codeMassToLogMsun(hotMass.tVircur_gal[*])   ;/ logMassBinSize
  totalMass.tViracc_gal[*] = codeMassToLogMsun(totalMass.tVircur_gal[*]) ;/ logMassBinSize
  
  coldMass.tViracc_gmem[*]  = codeMassToLogMsun(coldMass.tViracc_gmem[*])  ;/ logMassBinSize
  hotMass.tViracc_gmem[*]   = codeMassToLogMsun(hotMass.tViracc_gmem[*])   ;/ logMassBinSize
  totalMass.tViracc_gmem[*] = codeMassToLogMsun(totalMass.tViracc_gmem[*]) ;/ logMassBinSize
  
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

; haloMassBinAccRate(): bin accretion rate (in cold/hot) as a function of halo mass

function haloMassBinAccRate, sP=sP, sgSelect=sgSelect, accMode=accMode, radInd=radInd

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; current time
  h = loadSnapshotHeader(sP=sP)
  curtime = 1/h.time - 1 ; redshift
  curtime = redshiftToAgeFlat(curtime)*1e9 ; yr  
  
  timeWindow = curtime - redshiftToAgeFlat(4.0)*1e9 ; go back to z=4
  print,timeWindow/1e6,' Myr'

  ; config
  TcutVals = [5.3,5.4,5.5,5.6,5.7] ; for constant threshold
  TvirVals = [1.1,1.0,0.9,0.75,0.5] ; for Tviracc threshold
  
  nCuts = n_elements(TcutVals)
  nVirs = n_elements(TvirVals)
  
  xrange = [9.5,12.6]
  
  logMassBinSize = 0.3 / (sP.res/128)
  
  ; make a uniform gas selection at the start
  at = accretionTimes(sP=sP)
  mt = mergerTreeSubset(sP=sP)
  
  wAm = accModeInds(at=at,accMode=accMode,sP=sP,/mask)
    
  ; reverse histogram parent IDs of all particles/tracers in this selection
  gcIndOrig = mergerTreeRepParentIDs(mt=mt,sP=sP,/compactMtS)
  
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
  coldAccRate = { gal_const    : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                  gmem_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                  gal_tvircur  : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                  gmem_tvircur : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                  gal_tviracc  : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                  gmem_tviracc : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                  both_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                  both_tvircur : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                  both_tviracc : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                  gal_num      : lonarr(n_elements(mt.galcatIDList))         ,$
                  gmem_num     : lonarr(n_elements(mt.galcatIDList))          }
                  
  hotAccRate = { gal_const    : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                 gmem_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                 gal_tvircur  : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                 gmem_tvircur : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                 gal_tviracc  : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                 gmem_tviracc : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                 both_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                 both_tvircur : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                 both_tviracc : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                 gal_num      : lonarr(n_elements(mt.galcatIDList))         ,$
                 gmem_num     : lonarr(n_elements(mt.galcatIDList))          }
               
  ; loop over all tracked subgroups (galaxy)
  for i=0,n_elements(hist_gal)-1 do begin
    if hist_gal[i] gt 0 then begin
      ; list of indices of galaxy gas particles in this subgroup
      loc_inds_gal = rev_gal[rev_gal[i]:rev_gal[i+1]-1]
      
      ; corresponding accretion times for these particles
      loc_atime_gal = reform(at.accTime_gal[radInd,wAm.gal[loc_inds_gal]])
      loc_atime_gal = 1/loc_atime_gal - 1 ; redshift
      loc_atime_gal = redshiftToAgeFlat(loc_atime_gal)*1e9 ; yr
      
      ; make a count of those falling in the time window
      w = where(curtime - loc_atime_gal le timeWindow,nloc)
      
      coldAccRate.gal_num[i] = nloc
      
      ; maximum past temps, cur and acc tvirs for only those particles in the time window
      if nloc gt 0 then begin
        loc_maxt_gal = maxTemp.gal[loc_inds_gal[w]]
        loc_curtvir_gal = curTvir.gal[loc_inds_gal[w]]
        loc_acctvir_gal = accTvir.gal[loc_inds_gal[w]]
  
        ; count mass elements with Tmax below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_maxt_gal le TcutVals[j],count_cold,ncomp=count_hot)
          coldAccRate.gal_const[j,i] = count_cold
          hotAccRate.gal_const[j,i]  = count_hot
        endfor
        
        for j=0,nVirs-1 do begin
          ; count mass elements with Tmax below Tvir at current time
          w = where(10.0^loc_maxt_gal / 10.0^loc_curtvir_gal le TvirVals[j],count_cold,ncomp=count_hot)
          coldAccRate.gal_tvircur[j,i] = count_cold
          hotAccRate.gal_tvircur[j,i]  = count_hot
    
          ; count mass elements with Tmax below Tvir at accretion time
          w = where(10.0^loc_maxt_gal / 10.0^loc_acctvir_gal le TvirVals[j],count_cold,ncomp=count_hot)
          coldAccRate.gal_tviracc[j,i] = count_cold
          hotAccRate.gal_tviracc[j,i]  = count_hot
        endfor
      endif ;nloc>0
    endif
  endfor
  
  ; loop over all tracked subgroups (groupmem)
  for i=0,n_elements(hist_gmem)-1 do begin
    if hist_gmem[i] gt 0 then begin
      ; list of indices of galaxy gas particles in this subgroup
      loc_inds_gmem = rev_gmem[rev_gmem[i]:rev_gmem[i+1]-1]
      
      ; corresponding accretion times for these particles
      loc_atime_gmem = reform(at.accTime_gmem[radInd,wAm.gmem[loc_inds_gmem]])
      loc_atime_gmem = 1/loc_atime_gmem - 1 ; redshift
      loc_atime_gmem = redshiftToAgeFlat(loc_atime_gmem)*1e9 ; yr
      
      ; make a count of those falling in the time window
      w = where(curtime - loc_atime_gmem le timeWindow,nloc)
      
      coldAccRate.gmem_num[i] = nloc
      
      ; maximum past temps, cur and acc tvirs for only those particles in the time window
      loc_maxt_gmem = maxTemp.gmem[loc_inds_gmem[w]]
      loc_curtvir_gmem = curTvir.gmem[loc_inds_gmem[w]]
      loc_acctvir_gmem = accTvir.gmem[loc_inds_gmem[w]]

      ; count mass elements with Tmax below each constant temperature threshold
      for j=0,nCuts-1 do begin
        w = where(loc_maxt_gmem le TcutVals[j],count_cold,ncomp=count_hot)
        coldAccRate.gmem_const[j,i] = count_cold
        hotAccRate.gmem_const[j,i]  = count_hot
      endfor
      
      for j=0,nVirs-1 do begin
        ; count mass elements with Tmax below Tvir at current time
        w = where(10.0^loc_maxt_gmem / 10.0^loc_curtvir_gmem le TvirVals[j],count_cold,ncomp=count_hot)
        coldAccRate.gmem_tvircur[j,i] = count_cold
        hotAccRate.gmem_tvircur[j,i]  = count_hot
  
        ; count mass elements with Tmax below Tvir at accretion time
        w = where(10.0^loc_maxt_gmem / 10.0^loc_acctvir_gmem le TvirVals[j],count_cold,ncomp=count_hot)
        coldAccRate.gmem_tviracc[j,i] = count_cold
        hotAccRate.gmem_tviracc[j,i]  = count_hot
      endfor
    endif
  endfor
  
  ; convert total(counts) to msun/year
  if sP.trMCPerCell le 0 then massPerPart = sP.targetGasMass ; SPH or vel tracer
  if sP.trMCPerCell gt 0 then massPerPart = sP.trMassConst ; MC tracer
  
  coldAccRate.gal_const    *= massPerPart * units.UnitMass_in_Msun / timeWindow
  coldAccRate.gmem_const   *= massPerPart * units.UnitMass_in_Msun / timeWindow
  coldAccRate.gal_tvircur  *= massPerPart * units.UnitMass_in_Msun / timeWindow
  coldAccRate.gmem_tvircur *= massPerPart * units.UnitMass_in_Msun / timeWindow
  coldAccRate.gal_tviracc  *= massPerPart * units.UnitMass_in_Msun / timeWindow
  coldAccRate.gmem_tviracc *= massPerPart * units.UnitMass_in_Msun / timeWindow
  
  hotAccRate.gal_const    *= massPerPart * units.UnitMass_in_Msun / timeWindow
  hotAccRate.gmem_const   *= massPerPart * units.UnitMass_in_Msun / timeWindow
  hotAccRate.gal_tvircur  *= massPerPart * units.UnitMass_in_Msun / timeWindow
  hotAccRate.gmem_tvircur *= massPerPart * units.UnitMass_in_Msun / timeWindow
  hotAccRate.gal_tviracc  *= massPerPart * units.UnitMass_in_Msun / timeWindow
  hotAccRate.gmem_tviracc *= massPerPart * units.UnitMass_in_Msun / timeWindow
  
  ; create composite gal+gmem
  for j=0,nVirs-1 do begin
    coldAccRate.both_tvircur[j,*] = coldAccRate.gal_tvircur[j,*] + coldAccRate.gmem_tvircur[j,*]
    coldAccRate.both_tviracc[j,*] = coldAccRate.gal_tviracc[j,*] + coldAccRate.gmem_tviracc[j,*]
    
    hotAccRate.both_tvircur[j,*] = hotAccRate.gal_tvircur[j,*] + hotAccRate.gmem_tvircur[j,*]
    hotAccRate.both_tviracc[j,*] = hotAccRate.gal_tviracc[j,*] + hotAccRate.gmem_tviracc[j,*]
  endfor
  
  for j=0,nCuts-1 do begin
    coldAccRate.both_const[j,*] = coldAccRate.gal_const[j,*] + coldAccRate.gmem_const[j,*]
    hotAccRate.both_const[j,*] = hotAccRate.gal_const[j,*] + hotAccRate.gmem_const[j,*]
  endfor
  
  ; bin fractions into halo mass bins and make median lines
  ;logMassNbins  = floor((xrange[1]-xrange[0]) / logMassBinSize)
  ;logMassBins   = linspace(xrange[0],xrange[1],logMassNbins+1) ; edges
  ;logMassBinCen = linspace(xrange[0],xrange[1],logMassNbins+1) + logMassBinSize/2.0
  ;logMassBinCen = logMassBinCen[0:-2] ; remove last
  
  ; manual
  logMassBins=[9.5,10.0,10.15,10.3,10.45,10.6,10.75,10.9,11.0,$
               11.1,11.25,11.5,11.75,12.0,12.25,13.0]
  logMassNBins = n_elements(logMassBins)-1
  logMassBinCen = 0.5 * (logMassBins + shift(logMassBins,-1))
  logMassBinCen = logMassBinCen[0:-2]
  
  ; structures to store the binned values
  coldMedian = { gal_const    : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 gmem_const   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 both_const   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 gal_tVircur  : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 gmem_tVircur : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 both_tVircur : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 gal_tViracc  : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 gmem_tViracc : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 both_tViracc : fltarr(nVirs,logMassNbins) + !values.f_nan  }
                 
  hotMedian = { gal_const    : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                gmem_const   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                both_const   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                gal_tVircur  : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                gmem_tVircur : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                both_tVircur : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                gal_tViracc  : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                gmem_tViracc : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                both_tViracc : fltarr(nVirs,logMassNbins) + !values.f_nan  }
                
  totalMedian = { gal    : fltarr(logMassNbins) + !values.f_nan ,$
                  gmem   : fltarr(logMassNbins) + !values.f_nan ,$
                  both   : fltarr(logMassNbins) + !values.f_nan  }         
                 
  ; calculate median accretion rate in bins of halo mass
  for i=0,logMassNbins-1 do begin

    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1],count)
    
    if count gt 0 then begin
      for j=0,nVirs-1 do begin
        ; gal (hot+cold)
        coldMedian.gal_tVirCur[j,i]  = median(coldAccRate.gal_tvircur[j,w])
        coldMedian.gal_tVirAcc[j,i]  = median(coldAccRate.gal_tviracc[j,w])
        hotMedian.gal_tVirCur[j,i]   = median(hotAccRate.gal_tvircur[j,w])
        hotMedian.gal_tVirAcc[j,i]   = median(hotAccRate.gal_tviracc[j,w])
        
        ; gmem (hot+cold)
        coldMedian.gmem_tVirCur[j,i]  = median(coldAccRate.gmem_tvircur[j,w])
        coldMedian.gmem_tVirAcc[j,i]  = median(coldAccRate.gmem_tviracc[j,w])
        hotMedian.gmem_tVirCur[j,i]   = median(hotAccRate.gmem_tvircur[j,w])
        hotMedian.gmem_tVirAcc[j,i]   = median(hotAccRate.gmem_tviracc[j,w]) 
        
        ; both (hot+cold)
        coldMedian.both_tVirCur[j,i]  = median(coldAccRate.both_tvircur[j,w])
        coldMedian.both_tVirAcc[j,i]  = median(coldAccRate.both_tviracc[j,w])
        hotMedian.both_tVirCur[j,i]   = median(hotAccRate.both_tvircur[j,w])
        hotMedian.both_tVirAcc[j,i]   = median(hotAccRate.both_tviracc[j,w])
      endfor
      
      for j=0,nCuts-1 do begin
        coldMedian.gal_const[j,i]    = median(coldAccRate.gal_const[j,w])
        hotMedian.gal_const[j,i]     = median(hotAccRate.gal_const[j,w])
        
        coldMedian.gmem_const[j,i]   = median(coldAccRate.gmem_const[j,w])
        hotMedian.gmem_const[j,i]    = median(hotAccRate.gmem_const[j,w])
        
        coldMedian.both_const[j,i]   = median(coldAccRate.both_const[j,w])
        hotMedian.both_const[j,i]    = median(hotAccRate.both_const[j,w])
      endfor
      
      ; totals (same under any cold/hot definition)
      totalMedian.gal[i]   = median(coldAccRate.gal_const[0,w]+hotAccRate.gal_const[0,w])
      totalMedian.gmem[i]  = median(coldAccRate.gmem_const[0,w]+hotAccRate.gmem_const[0,w])
      totalMedian.both[i]  = median(coldAccRate.both_const[0,w]+hotAccRate.both_const[0,w])
    endif    
    
  endfor  
  
  ; debug: plot individual points
  xrange = [10.0,12.5]
  yrange = [0.1,50.0]
  
  constInd = 2 ; log(T)=5.5
  tVirInd  = 1 ; Tmax/Tvir=1

  hMasses = codeMassToLogMsun(mt.hMass[0,*])  
  
  start_PS, sP.plotPath + 'accRateRaw.const.'+accMode+'.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+$
    '_tw'+string(timeWindow/1e6,format='(i4)')+'.eps'
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="Const Accretion Rate [Msun/yr]",xtitle=textoidl("log ( M_{halo} ) [_{ }M_{sun }]")

    cgPlot,hMasses,coldAccRate.gal_const[constInd,*],color=getColor(1),psym=4,/overplot
    cgPlot,hMasses,hotAccRate.gal_const[constInd,*],color=getColor(3),psym=4,/overplot
    cgPlot,logMassBinCen,coldMedian.gal_const[constInd,*],color=getColor(1),line=0,/overplot
    cgPlot,logMassBinCen,hotMedian.gal_const[constInd,*],color=getColor(3),line=0,/overplot
    
    cgPlot,hMasses,coldAccRate.gmem_const[constInd,*],color=getColor(4),psym=4,/overplot
    cgPlot,hMasses,hotAccRate.gmem_const[constInd,*],color=getColor(5),psym=4,/overplot
    cgPlot,logMassBinCen,coldMedian.gmem_const[constInd,*],color=getColor(4),line=0,/overplot
    cgPlot,logMassBinCen,hotMedian.gmem_const[constInd,*],color=getColor(5),line=0,/overplot
    legend,['cold gal','hot gal','cold gmem','hot gmem'],textcolor=getColor([1,3,4,5],/name),box=0,/top,/left
  end_PS
  
  start_PS, sP.plotPath + 'accRateRaw.tvircur.'+accMode+'.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+$
    '_tw'+string(timeWindow/1e6,format='(i4)')+'.eps'
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="Tvircur Accretion Rate [Msun/yr]",xtitle=textoidl("log ( M_{halo} ) [_{ }M_{sun }]")
    
    cgPlot,hMasses,coldAccRate.gal_tvircur[tVirInd,*],color=getColor(1),psym=4,/overplot
    cgPlot,hMasses,hotAccRate.gal_tvircur[tVirInd,*],color=getColor(3),psym=4,/overplot
    cgPlot,logMassBinCen,coldMedian.gal_tvircur[tVirInd,*],color=getColor(1),line=0,/overplot
    cgPlot,logMassBinCen,hotMedian.gal_tvircur[tVirInd,*],color=getColor(3),line=0,/overplot
    
    cgPlot,hMasses,coldAccRate.gmem_tvircur[tVirInd,*],color=getColor(4),psym=4,/overplot
    cgPlot,hMasses,hotAccRate.gmem_tvircur[tVirInd,*],color=getColor(5),psym=4,/overplot
    cgPlot,logMassBinCen,coldMedian.gmem_tvircur[tVirInd,*],color=getColor(4),line=0,/overplot
    cgPlot,logMassBinCen,hotMedian.gmem_tvircur[tVirInd,*],color=getColor(5),line=0,/overplot
    legend,['cold gal','hot gal','cold gmem','hot gmem'],textcolor=getColor([1,3,4,5],/name),box=0,/top,/left
  end_PS
  
  start_PS, sP.plotPath + 'accRateRaw.tviracc.'+accMode+'.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+$
    '_tw'+string(timeWindow/1e6,format='(i4)')+'.eps'
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="Tviracc Accretion Rate [Msun/yr]",xtitle=textoidl("log ( M_{halo} ) [_{ }M_{sun }]")
    
    cgPlot,hMasses,coldAccRate.gal_tviracc[tVirInd,*],color=getColor(1),psym=4,/overplot
    cgPlot,hMasses,hotAccRate.gal_tviracc[tVirInd,*],color=getColor(3),psym=4,/overplot
    cgPlot,logMassBinCen,coldMedian.gal_tviracc[tVirInd,*],color=getColor(1),line=0,/overplot
    cgPlot,logMassBinCen,hotMedian.gal_tviracc[tVirInd,*],color=getColor(3),line=0,/overplot
    
    cgPlot,hMasses,coldAccRate.gmem_tviracc[tVirInd,*],color=getColor(4),psym=4,/overplot
    cgPlot,hMasses,hotAccRate.gmem_tviracc[tVirInd,*],color=getColor(5),psym=4,/overplot
    cgPlot,logMassBinCen,coldMedian.gmem_tviracc[tVirInd,*],color=getColor(4),line=0,/overplot
    cgPlot,logMassBinCen,hotMedian.gmem_tviracc[tVirInd,*],color=getColor(5),line=0,/overplot
    legend,['cold gal','hot gal','cold gmem','hot gmem'],textcolor=getColor([1,3,4,5],/name),box=0,/top,/left
  end_PS
  
  r = {coldMedian:coldMedian,hotMedian:hotMedian,totalMedian:totalMedian,$
       xrange:xrange,radInd:radInd,logMassBins:logMassBins,logMassBinCen:logMassBinCen,$
       TcutVals:TcutVals,TvirVals:TvirVals}
  return,r
end

; haloMassBinColdFracs(): bin cold fraction as a function of halo mass (using different definitions)

function haloMassBinColdFracs, sP=sP, sgSelect=sgSelect, accMode=accMode

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  TcutVals = [5.3,5.4,5.5,5.6,5.7] ; for constant threshold
  TvirVals = [1.1,1.0,0.9,0.75,0.5] ; for Tviracc threshold
  
  nCuts = n_elements(TcutVals)
  nVirs = n_elements(TvirVals)
  
  minNum = 6
  xrange = [9.5,12.5]
  
  logMassBinSize = 0.45 / (sP.res/128)
  
  ; make a uniform gas selection at the start
  at = accretionTimes(sP=sP)
  mt = mergerTreeSubset(sP=sP)
  
  wAm = accModeInds(at=at,accMode=accMode,sP=sP)
    
  ; reverse histogram parent IDs of all particles/tracers in this selection
  gcIndOrig = mergerTreeRepParentIDs(mt=mt,sP=sP,/compactMtS)
  
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
               gal_tvircur  : fltarr(nVirs,n_elements(mt.galcatIDList))         ,$
               gmem_tvircur : fltarr(nVirs,n_elements(mt.galcatIDList))         ,$
               gal_tviracc  : fltarr(nVirs,n_elements(mt.galcatIDList))         ,$
               gmem_tviracc : fltarr(nVirs,n_elements(mt.galcatIDList))         ,$
               both_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
               both_tvircur : fltarr(nVirs,n_elements(mt.galcatIDList))         ,$
               both_tviracc : fltarr(nVirs,n_elements(mt.galcatIDList))         ,$
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
      
      for j=0,nVirs-1 do begin
        ; count fraction Tmax below Tvir at current time
        w = where(10.0^loc_maxt_gal / 10.0^curTvir.gal[loc_inds_gal] le TvirVals[j],count_below)
        coldFrac.gal_tvircur[j,i] = float(count_below) / nloc
  
        ; count fraction Tmax below Tvir at accretion time
        w = where(10.0^loc_maxt_gal / 10.0^accTvir.gal[loc_inds_gal] le TvirVals[j],count_below)
        coldFrac.gal_tviracc[j,i] = float(count_below) / nloc
      endfor
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
      
      for j=0,nVirs-1 do begin
        ; count fraction Tmax below Tvir at current time
        w = where(10.0^loc_maxt_gmem / 10.0^curTvir.gmem[loc_inds_gmem] le TvirVals[j],count_below)
        coldFrac.gmem_tvircur[j,i] = float(count_below) / nloc
  
        ; count fraction Tmax below Tvir at accretion time
        w = where(10.0^loc_maxt_gmem / 10.0^accTvir.gmem[loc_inds_gmem] le TvirVals[j],count_below)
        coldFrac.gmem_tviracc[j,i] = float(count_below) / nloc
      endfor
    endif
  endfor
  
  ; create composite gal+gmem (Tmax)
  for j=0,nVirs-1 do $
  coldFrac.both_tvircur[j,*] = (coldFrac.gal_tvircur[j,*] * coldFrac.gal_num + $
                           coldFrac.gmem_tvircur[j,*] * coldFrac.gmem_num) / $
                          (coldFrac.gal_num + coldFrac.gmem_num)
  for j=0,nVirs-1 do $
  coldFrac.both_tviracc[j,*] = (coldFrac.gal_tviracc[j,*] * coldFrac.gal_num + $
                           coldFrac.gmem_tviracc[j,*] * coldFrac.gmem_num) / $
                          (coldFrac.gal_num + coldFrac.gmem_num)
  for j=0,nCuts-1 do $
    coldFrac.both_const[j,*] = (coldFrac.gal_const[j,*] * coldFrac.gal_num + $
                                coldFrac.gmem_const[j,*] * coldFrac.gmem_num) / $
                               (coldFrac.gal_num + coldFrac.gmem_num)  
  
  ; structures to store results (Tcur)
  coldFrac_cur = { gal_const    : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                   gmem_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                   gal_tvircur  : fltarr(nVirs,n_elements(mt.galcatIDList))         ,$
                   gmem_tvircur : fltarr(nVirs,n_elements(mt.galcatIDList))         ,$
                   gal_tviracc  : fltarr(nVirs,n_elements(mt.galcatIDList))         ,$
                   gmem_tviracc : fltarr(nVirs,n_elements(mt.galcatIDList))         ,$
                   both_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                   both_tvircur : fltarr(nVirs,n_elements(mt.galcatIDList))         ,$
                   both_tviracc : fltarr(nVirs,n_elements(mt.galcatIDList))         ,$
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
        
        for j=0,nVirs-1 do begin
          ; count fraction Tmax below Tvir at current time
          w = where(10.0^loc_curt_gal / 10.0^curTvir.gal[loc_inds_gal] le TvirVals[j],count_below)
          coldFrac_cur.gal_tvircur[j,i] = float(count_below) / count_nosfr
    
          ; count fraction Tmax below Tvir at accretion time
          w = where(10.0^loc_curt_gal / 10.0^accTvir.gal[loc_inds_gal] le TvirVals[j],count_below)
          coldFrac_cur.gal_tviracc[j,i] = float(count_below) / count_nosfr
        endfor
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
        
        for j=0,nVirs-1 do begin
          ; count fraction Tmax below Tvir at current time
          w = where(10.0^loc_curt_gmem / 10.0^curTvir.gmem[loc_inds_gmem] le TvirVals[j],count_below)
          coldFrac_cur.gmem_tvircur[j,i] = float(count_below) / count_nosfr
    
          ; count fraction Tmax below Tvir at accretion time
          w = where(10.0^loc_curt_gmem / 10.0^accTvir.gmem[loc_inds_gmem] le TvirVals[j],count_below)
          coldFrac_cur.gmem_tviracc[j,i] = float(count_below) / count_nosfr
        endfor
      endif ; cursfr!=0
    endif ; hist_gal>0
  endfor
                        
  ; create composite gal+gmem (Tcur)
  for j=0,nVirs-1 do $
  coldFrac_cur.both_tvircur[j,*] = (coldFrac_cur.gal_tvircur[j,*] * coldFrac_cur.gal_num + $
                                    coldFrac_cur.gmem_tvircur[j,*] * coldFrac_cur.gmem_num) / $
                                   (coldFrac_cur.gal_num + coldFrac_cur.gmem_num)
  for j=0,nVirs-1 do $
  coldFrac_cur.both_tviracc[j,*] = (coldFrac_cur.gal_tviracc[j,*] * coldFrac_cur.gal_num + $
                                    coldFrac_cur.gmem_tviracc[j,*] * coldFrac_cur.gmem_num) / $
                                   (coldFrac_cur.gal_num + coldFrac_cur.gmem_num)
  for j=0,nCuts-1 do $
    coldFrac_cur.both_const[j,*] = (coldFrac_cur.gal_const[j,*] * coldFrac_cur.gal_num + $
                                    coldFrac_cur.gmem_const[j,*] * coldFrac_cur.gmem_num) / $
                                   (coldFrac_cur.gal_num + coldFrac_cur.gmem_num)                       
  
  ; bin fractions into halo mass bins and make median lines
  ;logMassNbins  = floor((xrange[1]-xrange[0]) / logMassBinSize)
  ;logMassBins   = linspace(xrange[0],xrange[1],logMassNbins+1) ; edges
  ;logMassBinCen = linspace(xrange[0],xrange[1],logMassNbins+1) + logMassBinSize/2.0
  ;logMassBinCen = logMassBinCen[0:-2] ; remove last
  
  ; manual
  logMassBins=[9.5,10.0,10.15,10.3,10.45,10.6,10.75,10.9,11.0,$
               11.1,11.25,11.5,11.75,12.0,12.25,13.0]
  logMassNBins = n_elements(logMassBins)-1
  logMassBinCen = 0.5 * (logMassBins + shift(logMassBins,-1))
  logMassBinCen = logMassBinCen[0:-2]
 
  medianVals = { const_gal    : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 const_gmem   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 const_both   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 tVircur_gal  : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 tVircur_gmem : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 tVircur_both : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 tViracc_gal  : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 tViracc_gmem : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 tViracc_both : fltarr(nVirs,logMassNbins) + !values.f_nan  }
                 
  medianVals_cur = { const_gal    : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                     const_gmem   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                     const_both   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                     tVircur_gal  : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                     tVircur_gmem : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                     tVircur_both : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                     tViracc_gal  : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                     tViracc_gmem : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                     tViracc_both : fltarr(nVirs,logMassNbins) + !values.f_nan  }
  
  ; calculate median in bins (Tmax) and enforce minimum particle numbers
  for i=0,logMassNbins-1 do begin
    ; gal (Tmax)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac.gal_num ge minNum,count)
    if count gt 0 then begin
      for j=0,nVirs-1 do medianVals.tVircur_gal[j,i] = median(coldFrac.gal_tvircur[j,w])
      for j=0,nVirs-1 do medianVals.tViracc_gal[j,i] = median(coldFrac.gal_tviracc[j,w])
      for j=0,nCuts-1 do medianVals.const_gal[j,i]   = median(coldFrac.gal_const[j,w])
    endif
    
    ; gmem (Tmax)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac.gmem_num ge minNum,count)
    if count gt 0 then begin
      for j=0,nVirs-1 do medianVals.tVircur_gmem[j,i] = median(coldFrac.gmem_tvircur[j,w])
      for j=0,nVirs-1 do medianVals.tViracc_gmem[j,i] = median(coldFrac.gmem_tviracc[j,w])
      for j=0,nCuts-1 do medianVals.const_gmem[j,i]   = median(coldFrac.gmem_const[j,w])
    endif
    
    ; both (Tmax)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac.gal_num+coldFrac.gmem_num ge minNum,count)
    if count gt 0 then begin
      for j=0,nVirs-1 do medianVals.tVircur_both[j,i] = median(coldFrac.both_tvircur[j,w])
      for j=0,nVirs-1 do medianVals.tViracc_both[j,i] = median(coldFrac.both_tviracc[j,w])
      for j=0,nCuts-1 do medianVals.const_both[j,i]   = median(coldFrac.both_const[j,w])
    endif
    
    ; gal (Tcur)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac_cur.gal_num ge minNum,count)
    if count gt 0 then begin
      for j=0,nVirs-1 do medianVals_cur.tVircur_gal[j,i] = median(coldFrac_cur.gal_tvircur[j,w])
      for j=0,nVirs-1 do medianVals_cur.tViracc_gal[j,i] = median(coldFrac_cur.gal_tviracc[j,w])
      for j=0,nCuts-1 do medianVals_cur.const_gal[j,i]   = median(coldFrac_cur.gal_const[j,w])
    endif
    
    ; gmem (Tcur)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac_cur.gmem_num ge minNum,count)
    if count gt 0 then begin
      for j=0,nVirs-1 do medianVals_cur.tVircur_gmem[j,i] = median(coldFrac_cur.gmem_tvircur[j,w])
      for j=0,nVirs-1 do medianVals_cur.tViracc_gmem[j,i] = median(coldFrac_cur.gmem_tviracc[j,w])
      for j=0,nCuts-1 do medianVals_cur.const_gmem[j,i]   = median(coldFrac_cur.gmem_const[j,w])
    endif
    
    ; both (Tcur)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac_cur.gal_num+coldFrac_cur.gmem_num ge minNum,count)
    if count gt 0 then begin
      for j=0,nVirs-1 do medianVals_cur.tVircur_both[j,i] = median(coldFrac_cur.both_tvircur[j,w])
      for j=0,nVirs-1 do medianVals_cur.tViracc_both[j,i] = median(coldFrac_cur.both_tviracc[j,w])
      for j=0,nCuts-1 do medianVals_cur.const_both[j,i]   = median(coldFrac_cur.both_const[j,w])
    endif    
  endfor
  
  r = {coldFrac:coldFrac,coldFrac_cur:coldFrac_cur,medianVals:medianVals,medianVals_cur:medianVals_cur,$
       logMassBinCen:logMassBinCen,xrange:xrange,TcutVals:TcutVals,TvirVals:TvirVals}
  return, r
  
end