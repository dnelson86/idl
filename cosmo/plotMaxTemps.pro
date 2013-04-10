; plotMaxTemps.pro
; gas accretion project - plots related to maximum past temperature of the gas
; dnelson oct.2012

; binTmaxHistos()

function binTmaxHistos, sP=sP, sgSelect=sgSelect, accMode=accMode, timeWindow=TW

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  binSizeLog = 0.15 / (sP.res/128)
  massBins   = [9.0,9.5,10.0,10.5,11.0,11.5,12.0] ; log(M)  
  
  mmRatio = [-2.2,1.2]-[0.0,0.0001] ; log(tmax/tvir)
  mmTemp  = [3.8,7.2]-[0.0,0.0001]  ; log(tmax [K])
  
  nMassBins  = n_elements(massBins)-1
  nRatioBins = ceil((mmRatio[1]-mmRatio[0])/binSizeLog)
  nTempBins  = ceil((mmTemp[1]-mmTemp[0])/binSizeLog)
  
  ; current time
  h = loadSnapshotHeader(sP=sP)
  curtime = 1/h.time - 1 ; redshift
  curtime = redshiftToAgeFlat(curtime)*1e9 ; yr
  
  snapTimes = redshiftToAgeFlat(snapNumToRedshift(/all,sP=sP))*1e9 ; yr
  
  ; time window to consider accretion over
  if ~keyword_set(TW) then message,'time window required (in Myr)'
  
  if str(TW) eq 'all' then begin
    timeWindow = curtime - redshiftToAgeFlat(6.0)*1e9 ; go back to z=6 (in yr)
  endif else begin
    timeWindow = TW * 1e6 ; convert input Myr to yr
  endelse

  ; check if save exists
  saveFilename = sP.derivPath + 'binnedVals/binTemp.' + sP.saveTag + '.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.mb' + str(n_elements(massBins)) + '.' + sgSelect + '.' + accMode + $
    '_tw' + str(TW) + '.sav'
  
  ; results exist, return
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif  

  ; load accretion times (need wAm for accMode selection) for timewindow
  at = accretionTimes(sP=sP)
  wAm = accModeInds(at=at,accMode=accMode,sP=sP)

  ; do the timewindow restriction immediately
  ; ---
  
    ; GALAXY: accretion defined as (rho,temp) joining time or 0.15rvir crossing time (most recent)
    loc_atime_gal = reform(at.accTimeRT_gal[wAm.gal])
    
    r_crossing_time = reform(at.accTime_gal[sP.radIndGalAcc,wAm.gal])
    w = where(r_crossing_time gt loc_atime_gal,count)
    if count gt 0 then loc_atime_gal[w] = r_crossing_time[w]
    
    loc_atime_gal = 1/loc_atime_gal - 1 ; redshift
    loc_atime_gal = redshiftToAgeFlat(loc_atime_gal)*1e9 ; yr
    
    ; make a count of those falling in the time window
    w_gal = where(curtime - loc_atime_gal le timeWindow,nloc)
    loc_atime_gal = !NULL
    
    ; STARS:
    loc_atime_stars = reform(at.accTimeRT_stars[wAm.stars])
    
    r_crossing_time = reform(at.accTime_stars[sP.radIndGalAcc,wAm.stars])
    w = where(r_crossing_time gt loc_atime_stars,count)
    if count gt 0 then loc_atime_stars[w] = r_crossing_time[w]
    
    ; convert from scale factor to age of the universe
    loc_atime_stars = 1/loc_atime_stars - 1 ; redshift
    loc_atime_stars = redshiftToAgeFlat(loc_atime_stars)*1e9 ; yr
    
    ; note: if no 0.15rvir crossing time exists (crossed as a star particle, not as a gas cell) then
    ;       set the accretion time to t=-1 which moves it outside any timeWindow (do not consider)
    w = where(r_crossing_time eq -1,count)
    if count gt 0 then loc_atime_stars[w] = -1.0
    
    ; make a count of those falling in the time window
    w_stars = where(curtime - loc_atime_stars le timeWindow,nloc)
    loc_atime_stars = !NULL
    
    ; GMEM:
    loc_atime_gmem = reform(at.accTime_gmem[sP.radIndHaloAcc,wAm.gmem])
    loc_atime_gmem = 1/loc_atime_gmem - 1 ; redshift
    loc_atime_gmem = redshiftToAgeFlat(loc_atime_gmem)*1e9 ; yr
    
    ; make a count of those falling in the time window
    w_gmem = where(curtime - loc_atime_gmem le timeWindow,nloc)
    loc_atime_gmem = !NULL
    
  ;print,'gal ',n_elements(w_gal),n_elements(wAm.gal)
  ;print,'gmem ',n_elements(w_gmem),n_elements(wAm.gmem)
  ;print,'stars ',n_elements(w_stars),n_elements(wAm.stars) 
  
  w   = !NULL
  at  = !NULL
  wAm = !NULL
  r_crossing_time = !NULL
  
  ; load temps and do TW subsets
  accTvir = gcSubsetProp(sP=sP,select=sgSelect,/accTvir,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  accTvir = { gal : accTvir.gal[w_gal], gmem : accTvir.gmem[w_gmem], stars : accTvir.stars[w_stars] }

  maxTemp = gcSubsetProp(sP=sP,select=sgSelect,/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)  
  maxTemp = { gal : maxTemp.gal[w_gal], gmem : maxTemp.gmem[w_gmem], stars : maxTemp.stars[w_stars] }

  ; load parent halo masses so we can make halo massbins
  parentMass = gcSubsetProp(sP=sP,select=sgSelect,/parMass,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  parentMass = { gal : parentMass.gal[w_gal], gmem : parentMass.gmem[w_gmem], stars : parentMass.stars[w_stars] }
  
  w_gal   = !NULL
  w_gmem  = !NULL
  w_stars = !NULL
  
  ; return array
  r = { massBins           : massBins                     ,$
        binSizeLog         : binSizeLog                   ,$
        binLocRatio        : fltarr(nRatioBins)           ,$
        binLocTemp         : fltarr(nTempBins)            ,$
                                                           $
        hGalTmaxTviracc    : fltarr(nMassBins,nRatioBins) ,$
        hGalTmax           : fltarr(nMassBins,nTempBins)  ,$
        hStarsTmaxTviracc  : fltarr(nMassBins,nRatioBins) ,$
        hStarsTmax         : fltarr(nMassBins,nTempBins)  ,$
        hBothTmaxTviracc   : fltarr(nMassBins,nRatioBins) ,$
        hBothTmax          : fltarr(nMassBins,nTempBins)  ,$
        hGmemTmaxTviracc   : fltarr(nMassBins,nRatioBins) ,$
        hGmemTmax          : fltarr(nMassBins,nTempBins)  ,$
                                                           $
        hGalGlobalTmaxTviracc   : fltarr(nRatioBins)      ,$
        hStarsGlobalTmaxTviracc : fltarr(nRatioBins)      ,$
        hBothGlobalTmaxTviracc  : fltarr(nRatioBins)      ,$
        hGmemGlobalTmaxTviracc  : fltarr(nRatioBins)      ,$
        hGalGlobalTmax          : fltarr(nTempBins)       ,$
        hStarsGlobalTmax        : fltarr(nTempBins)       ,$
        hBothGlobalTmax         : fltarr(nTempBins)       ,$
        hGmemGlobalTmax         : fltarr(nTempBins)        }
        
  ; bin
  for j=0,n_elements(massBins)-2 do begin
  
    w_gal   = where(parentMass.gal gt massBins[j] and parentMass.gal le massBins[j+1],count1)
    w_gmem  = where(parentMass.gmem gt massBins[j] and parentMass.gmem le massBins[j+1],count2)
    w_stars = where(parentMass.stars gt massBins[j] and parentMass.stars le massBins[j+1],count3)
    
    if count1 eq 0 or count2 eq 0 or count3 eq 0 then begin
      print,'WARNING' ; empty
      continue
    endif
    
    ; log(tmax/tviracc)
    vals = [10.0^maxTemp.gal[w_gal]/10.0^accTvir.gal[w_gal]]
    hist = histogram(alog10(vals),binsize=binsizeLog,min=mmRatio[0],max=mmRatio[1],loc=loc)    
    if n_elements(loc) ne nRatioBins then message,'error'
    r.hGalTmaxTviracc[j,*] = hist
    
    vals = [10.0^maxTemp.stars[w_stars]/10.0^accTvir.stars[w_stars]]
    hist = histogram(alog10(vals),binsize=binsizeLog,min=mmRatio[0],max=mmRatio[1],loc=loc)  
    if n_elements(loc) ne nRatioBins then message,'error'
    r.hStarsTmaxTviracc[j,*] = hist
    
    r.hBothTmaxTviracc[j,*] = r.hGalTmaxTviracc[j,*] + r.hStarsTmaxTviracc[j,*]
    
    vals = [10.0^maxTemp.gmem[w_gmem]/10.0^accTvir.gmem[w_gmem]]
    hist = histogram(alog10(vals),binsize=binsizeLog,min=mmRatio[0],max=mmRatio[1],loc=loc)    
    if n_elements(loc) ne nRatioBins then message,'error'
    r.hGmemTmaxTviracc[j,*] = hist
  
    ; log(tmax)
    hist = histogram(maxTemp.gal[w_gal],binsize=binsizeLog,min=mmTemp[0],max=mmTemp[1],loc=loc)
    if n_elements(loc) ne nTempBins then message,'error'
    r.hGalTmax[j,*] = hist
  
    hist = histogram(maxTemp.stars[w_stars],binsize=binsizeLog,min=mmTemp[0],max=mmTemp[1],loc=loc)
    if n_elements(loc) ne nTempBins then message,'error'
    r.hStarsTmax[j,*] = hist
    
    r.hBothTmax[j,*] = r.hGalTmax[j,*] + r.hStarsTmax[j,*]
    
    hist = histogram(maxTemp.gmem[w_gmem],binsize=binsizeLog,min=mmTemp[0],max=mmTemp[1],loc=loc)
    if n_elements(loc) ne nTempBins then message,'error'
    r.hGmemTmax[j,*] = hist
  endfor
  
  ; bin global tmax/tviracc
  vals = [10.0^maxTemp.gal/10.0^accTvir.gal]
  hist = histogram(alog10(vals),binsize=binsizeLog,min=mmRatio[0],max=mmRatio[1],loc=loc) 
  if n_elements(loc) ne nRatioBins then message,'error'
  r.hGalGlobalTmaxTviracc[*] = hist
  
  vals = [10.0^maxTemp.stars/10.0^accTvir.stars]
  hist = histogram(alog10(vals),binsize=binsizeLog,min=mmRatio[0],max=mmRatio[1],loc=loc)  
  if n_elements(loc) ne nRatioBins then message,'error'
  r.hStarsGlobalTmaxTviracc[*] = hist
  
  r.hBothGlobalTmaxTviracc[*] = r.hGalGlobalTmaxTviracc[*] + r.hStarsGlobalTmaxTviracc[*]
  
  vals = [10.0^maxTemp.gmem/10.0^accTvir.gmem]
  hist = histogram(alog10(vals),binsize=binsizeLog,min=mmRatio[0],max=mmRatio[1],loc=loc)   
  if n_elements(loc) ne nRatioBins then message,'error'
  r.hGmemGlobalTmaxTviracc[*] = hist
  
  r.binLocRatio = loc + binSizeLog*0.5 ; save ratio midbins
  
  ; bin global tmax
  hist = histogram(maxTemp.gal,binsize=binsizeLog,min=mmTemp[0],max=mmTemp[1],loc=loc)
  if n_elements(loc) ne nTempBins then message,'error'
  r.hGalGlobalTmax[*] = hist
  
  hist = histogram(maxTemp.stars,binsize=binsizeLog,min=mmTemp[0],max=mmTemp[1],loc=loc)
  if n_elements(loc) ne nTempBins then message,'error'
  r.hStarsGlobalTmax[*] = hist
  
  r.hBothGlobalTmax[*] = r.hGalGlobalTmax[*] + r.hStarsGlobalTmax[*]
  
  hist = histogram(maxTemp.gmem,binsize=binsizeLog,min=mmTemp[0],max=mmTemp[1],loc=loc)
  if n_elements(loc) ne nTempBins then message,'error'
  r.hGmemGlobalTmax[*] = hist
  
  r.binLocTemp = loc + binSizeLog*0.5 ; save temp [K] midbins

  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  return, r

end

; binTmaxHisto2D()

function binTmaxHisto2D, sP=sP, sgSelect=sgSelect, accMode=accMode, timeWindow=TW

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  binSizeRatio = 0.15 / (sP.res/128)
  binSizeTemp  = 0.10 / (sP.res/128)
  binSizeMass  = 0.15 / (sP.res/128)
  
  mmRatio = [-2.0,1.2]-[0.0,0.0001] ; log(tmax/tvir)
  mmTemp  = [4.0,7.0]-[0.0,0.0001]  ; log(tmax [K])
  mmMass  = [9.0,12.0]-[0.0,0.0001] ; log(mhalo [msun])
  
  nRatioBins = ceil((mmRatio[1]-mmRatio[0])/binSizeRatio)
  nTempBins  = ceil((mmTemp[1]-mmTemp[0])/binSizeTemp)
  nMassBins  = ceil((mmMass[1]-mmMass[0])/binSizeMass)
  
  ; current time
  h = loadSnapshotHeader(sP=sP)
  curtime = 1/h.time - 1 ; redshift
  curtime = redshiftToAgeFlat(curtime)*1e9 ; yr
  
  snapTimes = redshiftToAgeFlat(snapNumToRedshift(/all,sP=sP))*1e9 ; yr
  
  ; time window to consider accretion over
  if ~keyword_set(TW) then message,'time window required (in Myr)'
  
  if str(TW) eq 'all' then begin
    timeWindow = curtime - redshiftToAgeFlat(6.0)*1e9 ; go back to z=6 (in yr)
  endif else begin
    timeWindow = TW * 1e6 ; convert input Myr to yr
  endelse

  ; check if save exists
  saveFilename = sP.derivPath + 'binnedVals/binTemp2D.' + sP.saveTag + '.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.mb' + str(nMassBins) + '.' + sgSelect + '.' + accMode + $
    '_tw' + str(TW) + '.sav'
  
  ; results exist, return
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif  

  ; load accretion times (need wAm for accMode selection) for timewindow
  at = accretionTimes(sP=sP)
  wAm = accModeInds(at=at,accMode=accMode,sP=sP)

  ; do the timewindow restriction immediately
  ; ---
  
    ; GALAXY: accretion defined as (rho,temp) joining time or 0.15rvir crossing time (most recent)
    loc_atime_gal = reform(at.accTimeRT_gal[wAm.gal])
    
    r_crossing_time = reform(at.accTime_gal[sP.radIndGalAcc,wAm.gal])
    w = where(r_crossing_time gt loc_atime_gal,count)
    if count gt 0 then loc_atime_gal[w] = r_crossing_time[w]
    
    loc_atime_gal = 1/loc_atime_gal - 1 ; redshift
    loc_atime_gal = redshiftToAgeFlat(loc_atime_gal)*1e9 ; yr
    
    ; make a count of those falling in the time window
    w_gal = where(curtime - loc_atime_gal le timeWindow,nloc)
    loc_atime_gal = !NULL
    
    ; STARS:
    loc_atime_stars = reform(at.accTimeRT_stars[wAm.stars])
    
    r_crossing_time = reform(at.accTime_stars[sP.radIndGalAcc,wAm.stars])
    w = where(r_crossing_time gt loc_atime_stars,count)
    if count gt 0 then loc_atime_stars[w] = r_crossing_time[w]
    
    ; convert from scale factor to age of the universe
    loc_atime_stars = 1/loc_atime_stars - 1 ; redshift
    loc_atime_stars = redshiftToAgeFlat(loc_atime_stars)*1e9 ; yr
    
    ; note: if no 0.15rvir crossing time exists (crossed as a star particle, not as a gas cell) then
    ;       set the accretion time to t=-1 which moves it outside any timeWindow (do not consider)
    w = where(r_crossing_time eq -1,count)
    if count gt 0 then loc_atime_stars[w] = -1.0
    
    ; make a count of those falling in the time window
    w_stars = where(curtime - loc_atime_stars le timeWindow,nloc)
    loc_atime_stars = !NULL
    
    ; GMEM:
    loc_atime_gmem = reform(at.accTime_gmem[sP.radIndHaloAcc,wAm.gmem])
    loc_atime_gmem = 1/loc_atime_gmem - 1 ; redshift
    loc_atime_gmem = redshiftToAgeFlat(loc_atime_gmem)*1e9 ; yr
    
    ; make a count of those falling in the time window
    w_gmem = where(curtime - loc_atime_gmem le timeWindow,nloc)
    loc_atime_gmem = !NULL
  
  w   = !NULL
  at  = !NULL
  wAm = !NULL
  r_crossing_time = !NULL
  
  ; load temps and do TW subsets
  accTvir = gcSubsetProp(sP=sP,select=sgSelect,/accTvir,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  accTvir = { gal : accTvir.gal[w_gal], gmem : accTvir.gmem[w_gmem], stars : accTvir.stars[w_stars] }

  curTvir = gcSubsetProp(sP=sP,select=sgSelect,/virTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  curTvir = { gal : curTvir.gal[w_gal], gmem : curTvir.gmem[w_gmem], stars : curTvir.stars[w_stars] }

  maxTemp = gcSubsetProp(sP=sP,select=sgSelect,/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)  
  maxTemp = { gal : maxTemp.gal[w_gal], gmem : maxTemp.gmem[w_gmem], stars : maxTemp.stars[w_stars] }

  ; load parent halo masses so we can make halo massbins
  parentMass = gcSubsetProp(sP=sP,select=sgSelect,/parMass,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  parentMass = { gal : parentMass.gal[w_gal], gmem : parentMass.gmem[w_gmem], stars : parentMass.stars[w_stars] }
  
  w_gal   = !NULL
  w_gmem  = !NULL
  w_stars = !NULL
  
  ; return array
  r = { binSizeRatio       : binSizeRatio                 ,$
        binSizeTemp        : binSizeTemp                  ,$
        binSizeMass        : binSizeMass                  ,$
        mmRatio            : mmRatio                      ,$
        mmTemp             : mmTemp                       ,$
        mmMass             : mmMass                       ,$
                                                           $
        h2_tmax_gal      : fltarr(nMassBins,nTempBins)  ,$
        h2_tmax_stars    : fltarr(nMassBins,nTempBins)  ,$
        h2_tmax_both     : fltarr(nMassBins,nTempBins)  ,$
        h2_tmax_gmem     : fltarr(nMassBins,nTempBins)  ,$
        h2_tvircur_gal   : fltarr(nMassBins,nRatioBins) ,$
        h2_tvircur_stars : fltarr(nMassBins,nRatioBins) ,$
        h2_tvircur_both  : fltarr(nMassBins,nRatioBins) ,$
        h2_tvircur_gmem  : fltarr(nMassBins,nRatioBins) ,$
        h2_tviracc_gal   : fltarr(nMassBins,nRatioBins) ,$
        h2_tviracc_stars : fltarr(nMassBins,nRatioBins) ,$
        h2_tviracc_both  : fltarr(nMassBins,nRatioBins) ,$
        h2_tviracc_gmem  : fltarr(nMassBins,nRatioBins)  }
        
  ; bin Tmax
  r.h2_tmax_gal = hist_nd(transpose([[parentMass.gal],[maxTemp.gal]]),[binSizeMass,binSizeTemp],$
                          min=[mmMass[0],mmTemp[0]],max=[mmMass[1],mmTemp[1]])
  r.h2_tmax_stars = hist_nd(transpose([[parentMass.stars],[maxTemp.stars]]),[binSizeMass,binSizeTemp],$
                          min=[mmMass[0],mmTemp[0]],max=[mmMass[1],mmTemp[1]])
  r.h2_tmax_gmem = hist_nd(transpose([[parentMass.gmem],[maxTemp.gmem]]),[binSizeMass,binSizeTemp],$
                          min=[mmMass[0],mmTemp[0]],max=[mmMass[1],mmTemp[1]])
        
  ; bin Tmax/Tvircur
  vals = alog10(10.0^maxTemp.gal/10.0^curTvir.gal)
  r.h2_tvircur_gal = hist_nd(transpose([[parentMass.gal],[vals]]),[binSizeMass,binSizeRatio],$
                             min=[mmMass[0],mmRatio[0]],max=[mmMass[1],mmRatio[1]])
  vals = alog10(10.0^maxTemp.stars/10.0^curTvir.stars)                      
  r.h2_tvircur_stars = hist_nd(transpose([[parentMass.stars],[vals]]),[binSizeMass,binSizeRatio],$
                             min=[mmMass[0],mmRatio[0]],max=[mmMass[1],mmRatio[1]])
  vals = alog10(10.0^maxTemp.gmem/10.0^curTvir.gmem)                         
  r.h2_tvircur_gmem = hist_nd(transpose([[parentMass.gmem],[vals]]),[binSizeMass,binSizeRatio],$
                             min=[mmMass[0],mmRatio[0]],max=[mmMass[1],mmRatio[1]])
  
  ; bin Tmax/Tviracc
  vals = alog10(10.0^maxTemp.gal/10.0^accTvir.gal)
  r.h2_tviracc_gal = hist_nd(transpose([[parentMass.gal],[vals]]),[binSizeMass,binSizeRatio],$
                             min=[mmMass[0],mmRatio[0]],max=[mmMass[1],mmRatio[1]])
  vals = alog10(10.0^maxTemp.stars/10.0^accTvir.stars)                         
  r.h2_tviracc_stars = hist_nd(transpose([[parentMass.stars],[vals]]),[binSizeMass,binSizeRatio],$
                             min=[mmMass[0],mmRatio[0]],max=[mmMass[1],mmRatio[1]])
  vals = alog10(10.0^maxTemp.gmem/10.0^accTvir.gmem)                         
  r.h2_tviracc_gmem = hist_nd(transpose([[parentMass.gmem],[vals]]),[binSizeMass,binSizeRatio],$
                             min=[mmMass[0],mmRatio[0]],max=[mmMass[1],mmRatio[1]])
          
  ; uniform weighting by mass
  if sP.trMCPerCell gt 0 then massWt = sP.trMassConst * units.UnitMass_in_Msun
  if sP.trMCPerCell eq 0 then massWt = sP.targetGasMass * units.UnitMass_in_Msun
  if sP.trMCPerCell eq -1 then message,'error'
  
  r.h2_tmax_gal *= massWt & r.h2_tmax_stars *= massWt & r.h2_tmax_gmem *= massWt
  r.h2_tvircur_gal *= massWt & r.h2_tvircur_stars *= massWt & r.h2_tvircur_gmem *= massWt
  r.h2_tviracc_gal *= massWt & r.h2_tviracc_stars *= massWt & r.h2_tviracc_gmem *= massWt
  
  ; form composite histograms
  r.h2_tmax_both = r.h2_tmax_gal + r.h2_tmax_stars
  r.h2_tvircur_both = r.h2_tvircur_gal + r.h2_tvircur_stars
  r.h2_tviracc_both = r.h2_tviracc_gal + r.h2_tviracc_stars
  
  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  return, r

end

; plotTmaxHistos(); plot (1) the previous max temp normalized by tviracc for arepo vs. gadget, gal vs. 
;                   gmem, (2) same but unnormalized by tviracc, (3) global not binned by halo mass but
;                   normalized by each parent tviracc, (4) same global without normalization

pro plotTmaxHistos

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  runs       = ['gadget','tracer','feedback']
  redshift   = 2.0
  res        = 512
  sgSelect   = 'pri'
  timeWindow = 1000.0 ; Myr
  accModes   = ['smooth'] ;['all','smooth','clumpy','stripped']
  
  ; plot config
  lines   = [0,1] ; gal,gmem
  cInd    = 1 ; color index
  
  if isnumeric(timeWindow) then twStr = str(fix(timeWindow))
  if ~isnumeric(timeWindow) then twStr = "-"+timeWindow
  
  ; loop over requested accretion modes
  foreach accMode,accModes do begin
   
    ; load
    foreach run,runs,i do begin
      sP  = mod_struct(sP, 'sP'+str(i), simParams(res=res,run=run,redshift=redshift))
      bth = mod_struct(bth, 'bth'+str(i), $
        binTmaxHistos(sP=sP.(i),sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow))
    endforeach
    
    ; strings
    plotStr   = ''
    simNames  = []
    simColors = []
  
    foreach run,runs,i do begin
      plotStr   = plotStr + sP.(i).plotPrefix + '.'
      simNames  = [simNames, sP.(i).simName]
      simColors = [simColors, sP.(i).colors[cInd]]
    endforeach

    plotStr += str(res) + '_' + str(sP.(0).snap) + '_tw' + twStr + '_am-' + accMode
    
    ; plot (1) - 3x2 mass bins separated out and each panel with gadget+arepo, GAL vs. gmem (VIR)
    start_PS, sP.(0).plotPath + 'tmax_3x2_gal_tviracc.'+plotStr+'.eps', /big
      !p.thick += 1
      xrange = [-2.2,1.2]
      yrange = [6e-4,1.0]
      
      xtickv = [-2.0,-1.0,0.0,1.0]
      
      x0 = 0.15 & x1 = 0.42 & x2 = 0.69 & x3 = 0.96
      y0 = 0.15 & y1 = 0.55 & y2 = 0.95
      pos = list( [x0,y1,x1,y2] ,$ ; ul
                  [x1,y1,x2,y2] ,$ ; uc
                  [x2,y1,x3,y2] ,$ ; ur
                  [x0,y0,x1,y1] ,$ ; ll
                  [x1,y0,x2,y1] ,$ ; lc
                  [x2,y0,x3,y1] )  ; lr
      
      for j=0,n_elements(bth.(0).massBins)-2 do begin
        
        if j eq 0 or j eq 3 then ytickname = '' else ytickname = replicate(' ',10)
        if j ge 3 then xtickname = '' else xtickname = replicate(' ',10)
        if j gt 0 then noerase = 1 else noerase = 0
        
        cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,pos=pos[j],$
          ytitle="",xtitle="",xticks=3,xtickv=xtickv,xtickname=xtickname,ytickname=ytickname,noerase=noerase       
        
        cgPlot,[0,0],[8e-4,0.25],line=2,color=cgColor('black'),thick=!p.thick-0.0,/overplot
        
        for i=0,n_tags(sP)-1 do begin
          ; gal
          hist = bth.(i).hGalTmaxTviracc[j,*]
          cgPlot,bth.(i).binLocRatio,float(hist)/total(hist),line=lines[0],color=sP.(i).colors[cInd],/overplot

          ; gmem
          hist = bth.(i).hGmemTmaxTviracc[j,*]
          cgPlot,bth.(i).binLocRatio,float(hist)/total(hist),line=lines[1],color=sP.(i).colors[cInd],/overplot
        endfor
        
        ; legends
        massBinStr = string(bth.(0).massBins[j],format='(f4.1)') + ' < log(M) < ' + $
                     string(bth.(0).massBins[j+1],format='(f4.1)')
  
        cgText,mean(xrange),yrange[1]*0.4,massBinStr,charsize=!p.charsize-0.0,alignment=0.5
            
        if j eq 0 then $
          legend,simNames,textcolors=simColors,box=0,position=[-2.2,0.2],charsize=!p.charsize-0.27
      
      endfor
      
      legend,['central galaxy','halo atmosphere'],linestyle=[0,1],$
        color=cgColor('dark gray'),textcolors=['dark gray','dark gray'],$
        linesize=0.25,box=0,position=[0.08,0.09],/normal,charsize=!p.charsize-0.2
      
      ; axis labels
      cgText,0.05,mean([y0,y2]),"Gas Mass Fraction",alignment=0.5,orientation=90.0,/normal
      cgText,mean([x0,x3]),0.05,textoidl("log ( T_{max} / T_{vir,acc} )"),alignment=0.5,/normal
      
    end_PS
    
    ; plot (2) - 3x2 mass bins separated out and each panel with gadget+arepo, BOTH vs. gmem (VIR)
    start_PS, sP.(0).plotPath + 'tmax_3x2_both_tviracc.'+plotStr+'.eps', /big
      !p.thick += 1
      xrange = [-2.2,1.2]
      yrange = [6e-4,1.0]
      
      xtickv = [-2.0,-1.0,0.0,1.0]
      
      x0 = 0.15 & x1 = 0.42 & x2 = 0.69 & x3 = 0.96
      y0 = 0.15 & y1 = 0.55 & y2 = 0.95
      pos = list( [x0,y1,x1,y2] ,$ ; ul
                  [x1,y1,x2,y2] ,$ ; uc
                  [x2,y1,x3,y2] ,$ ; ur
                  [x0,y0,x1,y1] ,$ ; ll
                  [x1,y0,x2,y1] ,$ ; lc
                  [x2,y0,x3,y1] )  ; lr
      
      for j=0,n_elements(bth.(0).massBins)-2 do begin
        
        if j eq 0 or j eq 3 then ytickname = '' else ytickname = replicate(' ',10)
        if j ge 3 then xtickname = '' else xtickname = replicate(' ',10)
        if j gt 0 then noerase = 1 else noerase = 0
        
        cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,pos=pos[j],$
          ytitle="",xtitle="",xticks=3,xtickv=xtickv,xtickname=xtickname,ytickname=ytickname,noerase=noerase       
        
        cgPlot,[0,0],[8e-4,0.25],line=2,color=cgColor('black'),thick=!p.thick-0.0,/overplot
        
        for i=0,n_tags(sP)-1 do begin
          ; gal
          hist = bth.(i).hBothTmaxTviracc[j,*]
          cgPlot,bth.(i).binLocRatio,float(hist)/total(hist),line=lines[0],color=sP.(i).colors[cInd],/overplot

          ; gmem
          hist = bth.(i).hGmemTmaxTviracc[j,*]
          cgPlot,bth.(i).binLocRatio,float(hist)/total(hist),line=lines[1],color=sP.(i).colors[cInd],/overplot
        endfor
                
        ; legends
        massBinStr = string(bth.(0).massBins[j],format='(f4.1)') + ' < log(M) < ' + $
                     string(bth.(0).massBins[j+1],format='(f4.1)')
  
        cgText,mean(xrange),yrange[1]*0.4,massBinStr,charsize=!p.charsize-0.0,alignment=0.5
            
        if j eq 0 then $
          legend,simNames,textcolors=simColors,box=0,position=[-2.2,0.2],charsize=!p.charsize-0.27
      
      endfor
      
      legend,['central galaxy','halo atmosphere'],linestyle=[0,1],$
        color=cgColor('dark gray'),textcolors=['dark gray','dark gray'],$
        linesize=0.25,box=0,position=[0.08,0.09],/normal,charsize=!p.charsize-0.2
      
      ; axis labels
      cgText,0.05,mean([y0,y2]),"Gas Mass Fraction",alignment=0.5,orientation=90.0,/normal
      cgText,mean([x0,x3]),0.05,textoidl("log ( T_{max} / T_{vir,acc} )"),alignment=0.5,/normal
      
    end_PS
    
    ; plot (3) - 3x2 mass bins separated out and each panel with gadget+arepo, GAL vs. gmem (CONST)
    start_PS, sP.(0).plotPath + 'tmax_3x2.'+plotStr+'.eps', /big
      !p.thick += 1
      xrange = [3.8,7.2]
      yrange = [6e-4,1.0]
      
      xtickv = [4.0,5.0,6.0,7.0]
      
      for j=0,n_elements(bth.(0).massBins)-2 do begin
        
        if j eq 0 or j eq 3 then ytickname = '' else ytickname = replicate(' ',10)
        if j ge 3 then xtickname = '' else xtickname = replicate(' ',10)
        if j gt 0 then noerase = 1 else noerase = 0
        
        cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,pos=pos[j],$
          ytitle="",xtitle="",xticks=3,xtickv=xtickv,xtickname=xtickname,ytickname=ytickname,noerase=noerase       
        
        cgPlot,[0,0],[8e-4,0.25],line=2,color=cgColor('black'),thick=!p.thick-0.0,/overplot
        
        for i=0,n_tags(sP)-1 do begin
          ; gal
          hist = bth.(i).hGalTmax[j,*]
          cgPlot,bth.(i).binLocTemp,float(hist)/total(hist),line=lines[0],color=sP.(i).colors[cInd],/overplot

          ; gmem
          hist = bth.(i).hGmemTmax[j,*]
          cgPlot,bth.(i).binLocTemp,float(hist)/total(hist),line=lines[1],color=sP.(i).colors[cInd],/overplot
        endfor
        
        ; legends
        massBinStr = string(bth.(0).massBins[j],format='(f4.1)') + ' < log(M) < ' + $
                     string(bth.(0).massBins[j+1],format='(f4.1)')
  
        cgText,mean(xrange),yrange[1]*0.4,massBinStr,charsize=!p.charsize-0.0,alignment=0.5
            
        if j eq 0 then $
          legend,simNames,textcolors=simColors,box=0,position=[-2.2,0.2],charsize=!p.charsize-0.27
      
      endfor
      
      legend,['central galaxy','halo atmosphere'],linestyle=[0,1],$
        color=cgColor('dark gray'),textcolors=['dark gray','dark gray'],$
        linesize=0.25,box=0,position=[0.08,0.09],/normal,charsize=!p.charsize-0.2
      
      ; axis labels
      cgText,0.05,mean([y0,y2]),"Gas Mass Fraction",alignment=0.5,orientation=90.0,/normal
      cgText,mean([x0,x3]),0.05,textoidl("log ( T_{max} )"),alignment=0.5,/normal
      
    end_PS
  
  endforeach ; accModes
  stop
end

; plotTmaxHisto2D(): plot 2D histogram of Tmax vs. halo mass (e.g. fig 8 of vdv11a)

pro plotTmaxHisto2D

  compile_opt idl2, hidden, strictarr, strictarrsubs
  !except = 0 ;suppress floating point underflow/overflow errors
  units = getUnits()
  
  ; config
  runs       = ['feedback','tracer']
  sgSelect   = 'pri'
  timeWindow = 1000.0 ; Myr
  redshift   = 2.0
  res        = 512
  accModes   = ['all','smooth'] ;['all','smooth','clumpy','stripped']  
  
  ; plot config
  exp    = 0.5   ; gamma exponent for non-linear color scaling
  ndivs  = 5     ; number of divisions on colorbar   
  Tc_val = 5.5   ; log(K) for constant temp line
  lines  = [2,0] ; Tc,Tvir line styles
  colors = ['dark gray','dark gray'] ; Tc,Tvir line colors
  cInd   = 1     ; color index for simName labels
  
  ; loop over requested accretion modes
  foreach accMode,accModes do begin
   
    ; load
    sP1 = simParams(res=res,run=runs[0],redshift=redshift)
    sP2 = simParams(res=res,run=runs[1],redshift=redshift)
    hh1 = binTmaxHisto2D(sP=sP1,sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow)
    hh2 = binTmaxHisto2D(sP=sP2,sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow)
  
    plotStr = accMode+'.'+sP1.plotPrefix+'.'+sP2.plotPrefix+'.'+str(res)+'_'+str(sP1.snap)
  
    ; plot (1) - two runs side by side (gal tmax)
    start_PS, sP1.plotPath + 'temp2d_comp_tmax_both.'+plotStr+'.eps', xs=9.0, ys=4.0
          
      loadColorTable, 'helix', /reverse      
      ytitle = textoidl('log ( T_{max} ) [K]')
      
      ; first run histogram
      h2mt = hh1.h2_tmax_gal
      yrange = [hh1.mmTemp[0],(round(hh1.mmTemp[1])*10)/10]
      
      tvim,h2mt^exp,pcharsize=!p.charsize,scale=0,clip=-1,$
           xtitle="",ytitle=ytitle,$
           barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=hh1.mmMass,yrange=yrange,$
           xticks=2,xtickv=[9.0,10.0,11.0],position=[0.1,0.18,0.47,0.94]
           
      ; temp lines and legend
      tvir_vals = alog10(codeMassToVirTemp(10.0^hh1.mmMass/units.UnitMass_in_Msun,sP=sP1))
      
      cgPlot,hh1.mmMass,[Tc_val,Tc_val],line=lines[0],color=cgColor(colors[0]),/overplot
      cgPlot,hh1.mmMass,tvir_vals,line=lines[1],color=cgColor(colors[1]),/overplot
      cgText,hh1.mmMass[0]+0.15,Tc_val+0.15,sP1.simName,color=sP1.colors[cInd],charsize=!p.charsize-0.2
           
      legend,textoidl(['T_c = 5.5','T_{vir }(z='+str(fix(redshift))+')']),linestyle=lines,$
        linesize=0.3,position=[9.05,6.85],box=0,spacing=!p.charsize+0.5
        
      ; second run histogram
      h2mt = hh2.h2_tmax_gal
      yrange = [hh2.mmTemp[0],(round(hh2.mmTemp[1])*10)/10]
      
      loadColorTable, 'helix', /reverse
      tvim,h2mt^exp,pcharsize=!p.charsize,scale=0,clip=-1,$
           xtitle="",ytitle="",$
           barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=hh2.mmMass,yrange=yrange,$
           xticks=3,xtickv=[9.0,10.0,11.0,12.0],position=[0.47,0.18,0.84,0.94],$
           ytickname=replicate(' ',10),/noerase
           
      ; colorbar
      barvals = findgen(ndivs+1)/ndivs*(max(h2mt^exp)-min(h2mt^exp)) + min(h2mt^exp)
      ticknames = textoidl(str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)')))
      ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
      
      loadColorTable, 'helix', /reverse
      cgColorbar,bottom=1,range=minmax(h2mt),position=[0.86,0.18,0.9,0.94],$
         /vertical,/right,title="",$
         divisions=ndivs,ticknames=ticknames,ncolors=255,charsize=!p.charsize-0.2
           
      ; temp lines
      tvir_vals = alog10(codeMassToVirTemp(10.0^hh2.mmMass/units.UnitMass_in_Msun,sP=sP2))
      
      cgPlot,hh2.mmMass,[Tc_val,Tc_val],line=lines[0],color=cgColor(colors[0]),/overplot
      cgPlot,hh2.mmMass,tvir_vals,line=lines[1],color=cgColor(colors[1]),/overplot
      cgText,hh2.mmMass[0]+0.15,Tc_val+0.15,sP2.simName,color=sP2.colors[cInd],charsize=!p.charsize-0.2
           
      ; labels
      cgText,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),(0.1+0.8)/2,0.04,alignment=0.5,/normal
      cgText,textoidl("M_{gas,tr} [_{ }log h^{-1} M_{sun }]"),0.98,(0.18+0.94)/2,alignment=0.5,orientation=90,/normal
                
    end_PS
    
    ; plot (2) - two runs side by side (gal tviracc)
    start_PS, sP1.plotPath + 'temp2d_comp_tviracc_both.'+plotStr+'.eps', xs=9.0, ys=4.0
          
      loadColorTable, 'helix', /reverse      
      ytitle = textoidl('log ( T_{max} / T_{vir,acc} )')
      
      ; first run histogram
      h2mt = hh1.h2_tviracc_gal
      yrange = [hh1.mmRatio[0],(round(hh1.mmRatio[1])*10)/10]
      
      tvim,h2mt^exp,pcharsize=!p.charsize,scale=0,clip=-1,$
           xtitle="",ytitle=ytitle,$
           barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=hh1.mmMass,yrange=yrange,$
           xticks=2,xtickv=[9.0,10.0,11.0],position=[0.1,0.18,0.47,0.94]
           
      ; temp lines and legend
      tvir_vals = alog10(10.0^Tc_val / codeMassToVirTemp(10.0^hh1.mmMass/units.UnitMass_in_Msun,sP=sP1))
      
      cgPlot,hh1.mmMass,tvir_vals,line=lines[0],color=cgColor(colors[1]),/overplot
      cgText,hh1.mmMass[0]+0.15,yrange[0]+0.15,sP1.simName,color=sP1.colors[cInd],charsize=!p.charsize-0.2
           
      ; second run histogram
      h2mt = hh2.h2_tviracc_gal
      yrange = [hh2.mmRatio[0],(round(hh2.mmRatio[1])*10)/10]
      
      loadColorTable, 'helix', /reverse
      tvim,h2mt^exp,pcharsize=!p.charsize,scale=0,clip=-1,$
           xtitle="",ytitle="",$
           barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=hh1.mmMass,yrange=yrange,$
           xticks=3,xtickv=[9.0,10.0,11.0,12.0],position=[0.47,0.18,0.84,0.94],$
           ytickname=replicate(' ',10),/noerase
           
      ; colorbar
      barvals = findgen(ndivs+1)/ndivs*(max(h2mt^exp)-min(h2mt^exp)) + min(h2mt^exp)
      ticknames = textoidl(str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)')))
      ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
      
      loadColorTable, 'helix', /reverse
      cgColorbar,bottom=1,range=minmax(h2mt),position=[0.86,0.18,0.9,0.94],$
         /vertical,/right,title="",$
         divisions=ndivs,ticknames=ticknames,ncolors=255,charsize=!p.charsize-0.2
           
      ; temp lines
      cgPlot,hh2.mmMass,tvir_vals,line=lines[0],color=cgColor(colors[1]),/overplot
      cgText,hh2.mmMass[0]+0.15,yrange[0]+0.15,sP2.simName,color=sP2.colors[cInd],charsize=!p.charsize-0.2
           
      ; labels
      cgText,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),(0.1+0.8)/2,0.04,alignment=0.5,/normal
      cgText,textoidl("M_{gas,tr} [_{ }log h^{-1} M_{sun }]"),0.98,(0.18+0.94)/2,alignment=0.5,orientation=90,/normal
                
    end_PS
    
    ; plot (3) - two runs side by side (gal+gmem 2x2 tmax)
    start_PS, sP1.plotPath + 'temp2d_comp_tmax_gal_gmem.'+plotStr+'.eps', xs=9.0, ys=8.0
          
      loadColorTable, 'helix', /reverse      
      ytitle = textoidl('log ( T_{max} ) [K]')
      
      pos = list([0.1,0.56,0.47,0.94],$
                 [0.47,0.56,0.84,0.94],$
                 [0.1,0.1,0.47,0.48],$
                 [0.47,0.1,0.84,0.48])
      
      ; UL - first run gal
      h2mt = hh1.h2_tmax_gal
      yrange = [hh1.mmTemp[0],(round(hh1.mmTemp[1])*10)/10]
      
      tvim,h2mt^exp,pcharsize=!p.charsize,scale=0,clip=-1,$
           xtitle="",ytitle=ytitle,barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=hh1.mmMass,yrange=yrange,$
           xticks=2,xtickv=[9.0,10.0,11.0],xtickname=replicate(' ',10),position=pos[0]
           
      ; temp lines and legend
      tvir_vals = alog10(codeMassToVirTemp(10.0^hh1.mmMass/units.UnitMass_in_Msun,sP=sP1))
      
      cgPlot,hh1.mmMass,[Tc_val,Tc_val],line=lines[0],color=cgColor(colors[0]),/overplot
      cgPlot,hh1.mmMass,tvir_vals,line=lines[1],color=cgColor(colors[1]),/overplot
      cgText,hh1.mmMass[0]+0.15,Tc_val+0.15,sP1.simName,color=sP1.colors[cInd],charsize=!p.charsize-0.2
           
      legend,textoidl(['T_c = 5.5','T_{vir }(z='+str(fix(redshift))+')']),linestyle=lines,$
        linesize=0.3,position=[9.05,6.85],box=0,spacing=!p.charsize+0.5
        
      ; UR - second run gal
      h2mt = hh2.h2_tmax_gal
      yrange = [hh2.mmTemp[0],(round(hh2.mmTemp[1])*10)/10]
      
      loadColorTable, 'helix', /reverse
      tvim,h2mt^exp,pcharsize=!p.charsize,scale=0,clip=-1,$
           xtitle="",ytitle="",$
           barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=hh1.mmMass,yrange=yrange,$
           xticks=3,xtickv=[9.0,10.0,11.0,12.0],xtickname=replicate(' ',10),position=pos[1],$
           ytickname=replicate(' ',10),/noerase
           
      ; temp lines
      tvir_vals = alog10(codeMassToVirTemp(10.0^hh2.mmMass/units.UnitMass_in_Msun,sP=sP2))
      
      cgPlot,hh2.mmMass,[Tc_val,Tc_val],line=lines[0],color=cgColor(colors[0]),/overplot
      cgPlot,hh2.mmMass,tvir_vals,line=lines[1],color=cgColor(colors[1]),/overplot
      cgText,hh2.mmMass[0]+0.15,Tc_val+0.15,sP2.simName,color=sP2.colors[cInd],charsize=!p.charsize-0.2
      
      ; colorbar
      barvals = findgen(ndivs+1)/ndivs*(max(h2mt^exp)-min(h2mt^exp)) + min(h2mt^exp)
      ticknames = textoidl(str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)')))
      ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
      
      loadColorTable, 'helix', /reverse
      cgColorbar,bottom=1,range=minmax(h2mt),position=[pos[1,2]+0.02,pos[1,1],pos[1,2]+0.06,pos[1,3]],$
         /vertical,/right,title="",$
         divisions=ndivs,ticknames=ticknames,ncolors=255,charsize=!p.charsize-0.2
      
      ; LL - first run gmem
      h2mt = hh1.h2_tmax_gmem
      yrange = [hh1.mmTemp[0],(round(hh1.mmTemp[1])*10)/10]
      
      loadColorTable, 'helix', /reverse
      tvim,h2mt^exp,pcharsize=!p.charsize,scale=0,clip=-1,$
           xtitle="",ytitle=ytitle,barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=hh1.mmMass,yrange=yrange,$
           xticks=2,xtickv=[9.0,10.0,11.0],position=pos[2],/noerase
      
      ; temp lines and legend
      tvir_vals = alog10(codeMassToVirTemp(10.0^hh1.mmMass/units.UnitMass_in_Msun,sP=sP1))
      
      cgPlot,hh1.mmMass,[Tc_val,Tc_val],line=lines[0],color=cgColor(colors[0]),/overplot
      cgPlot,hh1.mmMass,tvir_vals,line=lines[1],color=cgColor(colors[1]),/overplot
      cgText,hh1.mmMass[0]+0.15,Tc_val+0.15,sP1.simName,color=sP1.colors[cInd],charsize=!p.charsize-0.2
           
      legend,textoidl(['T_c = 5.5','T_{vir }(z='+str(fix(redshift))+')']),linestyle=lines,$
        linesize=0.3,position=[9.05,6.85],box=0,spacing=!p.charsize+0.5
      
      ; LR - second run gmem
      h2mt = hh2.h2_tmax_gmem
      yrange = [hh2.mmTemp[0],(round(hh2.mmTemp[1])*10)/10]
      
      loadColorTable, 'helix', /reverse
      tvim,h2mt^exp,pcharsize=!p.charsize,scale=0,clip=-1,$
           xtitle="",ytitle="",$
           barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=hh1.mmMass,yrange=yrange,$
           xticks=3,xtickv=[9.0,10.0,11.0,12.0],position=pos[3],$
           ytickname=replicate(' ',10),/noerase
      
      ; temp lines
      tvir_vals = alog10(codeMassToVirTemp(10.0^hh2.mmMass/units.UnitMass_in_Msun,sP=sP2))
      
      cgPlot,hh2.mmMass,[Tc_val,Tc_val],line=lines[0],color=cgColor(colors[0]),/overplot
      cgPlot,hh2.mmMass,tvir_vals,line=lines[1],color=cgColor(colors[1]),/overplot
      cgText,hh2.mmMass[0]+0.15,Tc_val+0.15,sP2.simName,color=sP2.colors[cInd],charsize=!p.charsize-0.2

      ; colorbar
      barvals = findgen(ndivs+1)/ndivs*(max(h2mt^exp)-min(h2mt^exp)) + min(h2mt^exp)
      ticknames = textoidl(str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)')))
      ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
      
      loadColorTable, 'helix', /reverse
      cgColorbar,bottom=1,range=minmax(h2mt),position=[pos[3,2]+0.02,pos[3,1],pos[3,2]+0.06,pos[3,3]],$
         /vertical,/right,title="",$
         divisions=ndivs,ticknames=ticknames,ncolors=255,charsize=!p.charsize-0.2
           
      ; labels
      cgText,"Central Galaxy",pos[0,2],pos[1,3]+0.02,/normal,alignment=0.5
      cgText,"Halo Atmosphere",pos[0,2],pos[2,3]+0.02,/normal,alignment=0.5
      cgText,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),(0.1+0.8)/2,0.02,alignment=0.5,/normal
      cgText,textoidl("M_{gas} [_{ }log h^{-1} M_{sun }]"),0.98,mean([pos[1,1],pos[1,3]]),alignment=0.5,orientation=90,/normal
      cgText,textoidl("M_{gas} [_{ }log h^{-1} M_{sun }]"),0.98,mean([pos[3,1],pos[3,3]]),alignment=0.5,orientation=90,/normal
                
    end_PS
    
  endforeach
  
  !except = 1
end
