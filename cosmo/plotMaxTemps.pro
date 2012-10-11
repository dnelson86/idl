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
  sPg = simParams(res=512,run='gadget',redshift=2.0)
  sPa = simParams(res=512,run='tracer',redshift=2.0) ; f=-1 use velocity tracers

  sgSelect   = 'pri'
  timeWindow = 1000.0 ; Myr
  accModes   = ['all','smooth','clumpy','stripped']
  
  lines   = [1,0,2] ; tvircur,tviracc,const5.5
  thicks  = [4,6,8] ; 128,256,512
  cInd    = 1 ; color index
  
  foreach accMode,accModes do begin
   
    ; load
    thA = binTmaxHistos(sP=sPa,sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow)
    thG = binTmaxHistos(sP=sPg,sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow)
  
    ; plot (1) - 3x2 mass bins separated out and each panel with gadget+arepo, gal vs. gmem
    start_PS, sPg.plotPath + 'tmax_3x2_tviracc.'+accMode+'.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
              str(sPg.res)+'_'+str(sPg.snap)+'.eps', /big
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
      
      for j=0,n_elements(thG.massBins)-2 do begin
        
        if j eq 0 or j eq 3 then ytickname = '' else ytickname = replicate(' ',10)
        if j ge 3 then xtickname = '' else xtickname = replicate(' ',10)
        if j gt 0 then noerase = 1 else noerase = 0
        
        cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,pos=pos[j],$
          ytitle="",xtitle="",xticks=3,xtickv=xtickv,xtickname=xtickname,ytickname=ytickname,noerase=noerase       
        
        cgPlot,[0,0],[8e-4,0.25],line=2,color=fsc_color('black'),thick=!p.thick-0.0,/overplot
        
        ; (gal)
        hist = thG.hGalTmaxTviracc[j,*]
        cgPlot,thG.binLocRatio,float(hist)/total(hist),line=lines[1],color=sPg.colorsG[cInd],/overplot
    
        hist = thA.hGalTmaxTviracc[j,*]
        cgPlot,thA.binLocRatio,float(hist)/total(hist),line=lines[1],color=sPa.colorsA[cInd],/overplot
      
        ; (gmem)
        hist = thG.hGmemTmaxTviracc[j,*]
        cgPlot,thG.binLocRatio,float(hist)/total(hist),line=lines[0],color=sPg.colorsG[cInd],/overplot
    
        hist = thA.hGmemTmaxTviracc[j,*]
        cgPlot,thA.binLocRatio,float(hist)/total(hist),line=lines[0],color=sPa.colorsA[cInd],/overplot
      
        ; legends
        massBinStr = string(thG.massBins[j],format='(f4.1)') + ' < log(M) < ' + $
                     string(thG.massBins[j+1],format='(f4.1)')
  
        cgText,mean(xrange),yrange[1]*0.4,massBinStr,charsize=!p.charsize-0.0,alignment=0.5
            
        if j eq 0 then $
          legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,$
            position=[-2.0,0.1],charsize=!p.charsize-0.1,spacing=!p.charsize+0.5
      
      endfor
      
      legend,['central galaxy','halo atmosphere'],linestyle=[0,1],$
        color=cgColor('dark gray'),textcolors=['dark gray','dark gray'],$
        linesize=0.25,box=0,position=[0.08,0.09],/normal,charsize=!p.charsize-0.2
      
      ; axis labels
      cgText,0.05,mean([y0,y2]),"Gas Mass Fraction",alignment=0.5,orientation=90.0,/normal
      cgText,mean([x0,x3]),0.05,textoidl("log ( T_{max} / T_{vir,acc} )"),alignment=0.5,/normal
      
    end_PS
    
    ; plot (1b) - 3x2 mass bins separated out and each panel with gadget+arepo, gal vs. gmem
    start_PS, sPg.plotPath + 'tmax_3x2b_tviracc.'+accMode+'.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
              str(sPg.res)+'_'+str(sPg.snap)+'.eps', /big
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
      
      for j=0,n_elements(thG.massBins)-2 do begin
        
        if j eq 0 or j eq 3 then ytickname = '' else ytickname = replicate(' ',10)
        if j ge 3 then xtickname = '' else xtickname = replicate(' ',10)
        if j gt 0 then noerase = 1 else noerase = 0
        
        cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,pos=pos[j],$
          ytitle="",xtitle="",xticks=3,xtickv=xtickv,xtickname=xtickname,ytickname=ytickname,noerase=noerase       
        
        cgPlot,[0,0],[8e-4,0.25],line=2,color=fsc_color('black'),thick=!p.thick-0.0,/overplot
        
        ; (gal)
        hist = thG.hBothTmaxTviracc[j,*]
        cgPlot,thG.binLocRatio,float(hist)/total(hist),line=lines[1],color=sPg.colorsG[cInd],/overplot
    
        hist = thA.hBothTmaxTviracc[j,*]
        cgPlot,thA.binLocRatio,float(hist)/total(hist),line=lines[1],color=sPa.colorsA[cInd],/overplot
      
        ; (gmem)
        hist = thG.hGmemTmaxTviracc[j,*]
        cgPlot,thG.binLocRatio,float(hist)/total(hist),line=lines[0],color=sPg.colorsG[cInd],/overplot
    
        hist = thA.hGmemTmaxTviracc[j,*]
        cgPlot,thA.binLocRatio,float(hist)/total(hist),line=lines[0],color=sPa.colorsA[cInd],/overplot
      
        ; legends
        massBinStr = string(thG.massBins[j],format='(f4.1)') + ' < log(M) < ' + $
                     string(thG.massBins[j+1],format='(f4.1)')
  
        cgText,mean(xrange),yrange[1]*0.4,massBinStr,charsize=!p.charsize-0.0,alignment=0.5
            
        if j eq 0 then $
          legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,$
            position=[-2.0,0.1],charsize=!p.charsize-0.1,spacing=!p.charsize+0.5
      
      endfor
      
      legend,['central galaxy','halo atmosphere'],linestyle=[0,1],$
        color=cgColor('dark gray'),textcolors=['dark gray','dark gray'],$
        linesize=0.25,box=0,position=[0.08,0.09],/normal,charsize=!p.charsize-0.2
      
      ; axis labels
      cgText,0.05,mean([y0,y2]),"Gas Mass Fraction",alignment=0.5,orientation=90.0,/normal
      cgText,mean([x0,x3]),0.05,textoidl("log ( T_{max} / T_{vir,acc} )"),alignment=0.5,/normal
      
    end_PS
    
    ; plot (2) - 3x2 mass bins separated out and each panel with gadget+arepo, gal vs. gmem
    if 0 then begin
    start_PS, sPg.plotPath + 'tmax_3x2.'+accMode+'.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
              str(sPg.res)+'_'+str(sPg.snap)+'.eps', /big
      !p.thick += 1
      xrange = [3.8,7.2]
      yrange = [6e-4,1.0]
      
      xtickv = [4.0,5.0,6.0,7.0]
      
      for j=0,n_elements(thG.massBins)-2 do begin
        
        if j eq 0 or j eq 3 then ytickname = '' else ytickname = replicate(' ',10)
        if j ge 3 then xtickname = '' else xtickname = replicate(' ',10)
        if j gt 0 then noerase = 1 else noerase = 0
        
        cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,pos=pos[j],$
          ytitle="",xtitle="",xticks=3,xtickv=xtickv,xtickname=xtickname,ytickname=ytickname,noerase=noerase       
        
        cgPlot,[0,0],[8e-4,0.25],line=2,color=fsc_color('black'),thick=!p.thick-0.0,/overplot
        
        ; (gal)
        hist = thG.hGalTmax[j,*]
        cgPlot,thG.binLocRatio,float(hist)/total(hist),line=lines[1],color=sPg.colorsG[cInd],/overplot
    
        hist = thA.hGalTmax[j,*]
        cgPlot,thA.binLocRatio,float(hist)/total(hist),line=lines[1],color=sPa.colorsA[cInd],/overplot
      
        ; (gmem)
        hist = thG.hGmemTmax[j,*]
        cgPlot,thG.binLocRatio,float(hist)/total(hist),line=lines[0],color=sPg.colorsG[cInd],/overplot
    
        hist = thA.hGmemTmax[j,*]
        cgPlot,thA.binLocRatio,float(hist)/total(hist),line=lines[0],color=sPa.colorsA[cInd],/overplot
      
        ; legends
        massBinStr = string(thG.massBins[j],format='(f4.1)') + ' < log(M) < ' + $
                     string(thG.massBins[j+1],format='(f4.1)')
  
        cgText,mean(xrange),yrange[1]*0.4,massBinStr,charsize=!p.charsize-0.0,alignment=0.5
            
        if j eq 0 then $
          legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,$
            position=[-2.0,0.1],charsize=!p.charsize-0.1,spacing=!p.charsize+0.5
      
      endfor
      
      legend,['central galaxy','halo atmosphere'],linestyle=[0,1],$
        color=cgColor('dark gray'),textcolors=['dark gray','dark gray'],$
        linesize=0.25,box=0,position=[0.08,0.09],/normal,charsize=!p.charsize-0.2
      
      ; axis labels
      cgText,0.05,mean([y0,y2]),"Gas Mass Fraction",alignment=0.5,orientation=90.0,/normal
      cgText,mean([x0,x3]),0.05,textoidl("log ( T_{max} )"),alignment=0.5,/normal
      
    end_PS  
  
    ; plot (3) - global histos of tmax/tviracc
    start_PS, sPg.plotPath + 'tmax_global_tviracc.'+accMode+'.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
              str(sPg.res)+'_'+str(sPg.snap)+'.eps', /big
      
      xrange = [-2.0,1.5]
      yrange = [1e-3,1e3]
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
        ytitle=textoidl("M_{tot}"),xtitle=textoidl("log ( T_{max} / T_{vir,acc} )")
      
      ; histogram gadget (gal) differences
      hist = thG.hGalGlobalTmaxTviracc[*]
      cgPlot,thG.binLocRatio,float(hist)*sPg.targetGasMass,line=lines[1],color=sPg.colorsG[cInd],/overplot
  
      ; histogram tracer (gal) differences
      hist = thA.hGalGlobalTmaxTviracc[*]
      cgPlot,thA.binLocRatio,float(hist)*sPa.trMassConst,line=lines[1],color=sPa.colorsA[cInd],/overplot
    
      ; histogram gadget (gmem) differences
      hist = thG.hGmemGlobalTmaxTviracc[*]
      cgPlot,thG.binLocRatio,float(hist)*sPg.targetGasMass,line=lines[0],color=sPg.colorsG[cInd],/overplot
  
      ; histogram tracer (gmem) differences
      hist = thA.hGmemGlobalTmaxTviracc[*]
      cgPlot,thA.binLocRatio,float(hist)*sPa.trMassConst,line=lines[0],color=sPa.colorsA[cInd],/overplot
    
      legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,$
        /top,/right,spacing=!p.charsize+0.5
      
      legend,['central galaxy','halo atmosphere'],linestyle=[0,1],$
        color=cgColor('dark gray'),textcolors=['dark gray','dark gray'],$
        linesize=0.25,box=0,/top,/left
      
    end_PS
    
    ; plot (4) - global histos of unnormalized tmax
    start_PS, sPg.plotPath + 'tmax_global.'+accMode+'.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
              str(sPg.res)+'_'+str(sPg.snap)+'.eps', /big
      ;!p.thick += 1
      
      xrange = [3.5,7.5]
      yrange = [1e-3,1e3]
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
        ytitle=textoidl("M_{tot}"),xtitle=textoidl("log ( T_{max} )")
      
      ; histogram gadget (gal) differences
      hist = thG.hGalGlobalTmax[*]
      cgPlot,thG.binLocRatio,float(hist)*sPg.targetGasMass,line=lines[1],color=sPg.colorsG[cInd],/overplot
  
      ; histogram tracer (gal) differences
      hist = thA.hGalGlobalTmax[*]
      cgPlot,thA.binLocRatio,float(hist)*sPa.trMassConst,line=lines[1],color=sPa.colorsA[cInd],/overplot
    
      ; histogram gadget (gmem) differences
      hist = thG.hGmemGlobalTmax[*]
      cgPlot,thG.binLocRatio,float(hist)*sPg.targetGasMass,line=lines[0],color=sPg.colorsG[cInd],/overplot
  
      ; histogram tracer (gmem) differences
      hist = thA.hGmemGlobalTmax[*]
      cgPlot,thA.binLocRatio,float(hist)*sPa.trMassConst,line=lines[0],color=sPa.colorsA[cInd],/overplot
    
      legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,$
        /top,/right,spacing=!p.charsize+0.5
      
      legend,['central galaxy','halo atmosphere'],linestyle=[0,1],$
        color=cgColor('dark gray'),textcolors=['dark gray','dark gray'],$
        linesize=0.25,box=0,/top,/left
      
    end_PS
    endif ;0
  
  endforeach ; accModes
  
end

; plotTmaxHisto2D(): plot 2D histogram of Tmax vs. halo mass (e.g. fig 8 of vdv11a)

pro plotTmaxHisto2D

  compile_opt idl2, hidden, strictarr, strictarrsubs
  !except = 0 ;suppress floating point underflow/overflow errors
  units = getUnits()
  
  ; config
  sPg = simParams(res=512,run='gadget',redshift=2.0)
  sPa = simParams(res=512,run='tracer',redshift=2.0) ; f=-1 use velocity tracers

  sgSelect   = 'pri'
  timeWindow = 1000.0 ; Myr
  accModes   = ['all','smooth'] ;['all','smooth','clumpy','stripped']  
  
  exp    = 0.5   ; gamma exponent for non-linear color scaling
  ndivs  = 5     ; number of divisions on colorbar   
  Tc_val = 5.5   ; log(K) for constant temp line
  lines  = [2,0] ; Tc,Tvir line styles
  colors = ['dark gray','dark gray'] ; Tc,Tvir line colors
  cInd   = 1     ; color index for gadget/arepo labels
  
  foreach accMode,accModes do begin
   
    ; load
    h2A = binTmaxHisto2D(sP=sPa,sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow)
    h2G = binTmaxHisto2D(sP=sPg,sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow)
  
    ; plot (1) - arepo tmax (debug)
    start_PS, sPg.plotPath + 'temp2d_test_tmax_both.'+accMode+'.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
              str(sPg.res)+'_'+str(sPg.snap)+'.eps', /big
          
      loadColorTable, 'helix', /reverse      
            
      h2mt = h2A.h2_tmax_both
      yrange = h2A.mmTemp
      ytitle = textoidl('log ( T_{max} ) [K]')
      
      ; arepo histogram
      tvim,h2mt^exp,pcharsize=!p.charsize,scale=0,clip=-1,$
           xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),ytitle=ytitle,$
           barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=h2A.mmMass,yrange=yrange,$
           xticks=3,xtickv=[9.0,10.0,11.0,12.0],position=[0.1,0.14,0.8,0.94]
           
      ; colorbar
      barvals = findgen(ndivs+1)/ndivs*(max(h2mt^exp)-min(h2mt^exp)) + min(h2mt^exp)
      ticknames = textoidl(str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)')))
      ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
      colorbar,bottom=1,range=minmax(h2mt),position=[0.83,0.14,0.87,0.94],$
         /vertical,/right,title=textoidl("M_{gas,tr} [_{ }log h^{-1} M_{sun }]"),divisions=ndivs,ticknames=ticknames,ncolors=255
           
      ; temp lines and legend
      tvir_vals = alog10(codeMassToVirTemp(10.0^h2A.mmMass/units.UnitMass_in_Msun,sP=sPa))
      
      cgPlot,h2A.mmMass,[Tc_val,Tc_val],line=lines[0],color=cgColor(colors[0]),/overplot
      cgPlot,h2A.mmMass,tvir_vals,line=lines[1],color=cgColor(colors[1]),/overplot
           
      legend,textoidl(['T_c = 5.5','T_{vir }(M)']),linestyle=lines,$;textcolor=colors,$
        linesize=0.3,/top,/left,box=0,spacing=!p.charsize+0.5
                
    end_PS
    
    ; plot (2) - arepo and gadget side by side (tmax)
    start_PS, sPg.plotPath + 'temp2d_comp_tmax_both.'+accMode+'.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
              str(sPg.res)+'_'+str(sPg.snap)+'.eps', xs=9.0, ys=4.0
          
      loadColorTable, 'helix', /reverse      
      ytitle = textoidl('log ( T_{max} ) [K]')
      
      ; arepo histogram
      h2mt = h2A.h2_tmax_gal
      yrange = [h2A.mmTemp[0],(round(h2A.mmTemp[1])*10)/10]
      
      tvim,h2mt^exp,pcharsize=!p.charsize,scale=0,clip=-1,$
           xtitle="",ytitle=ytitle,$
           barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=h2A.mmMass,yrange=yrange,$
           xticks=2,xtickv=[9.0,10.0,11.0],position=[0.1,0.18,0.47,0.94]
           
      ; temp lines and legend
      tvir_vals = alog10(codeMassToVirTemp(10.0^h2A.mmMass/units.UnitMass_in_Msun,sP=sPa))
      
      cgPlot,h2A.mmMass,[Tc_val,Tc_val],line=lines[0],color=cgColor(colors[0]),/overplot
      cgPlot,h2A.mmMass,tvir_vals,line=lines[1],color=cgColor(colors[1]),/overplot
      cgText,h2A.mmMass[0]+0.3,Tc_val+0.2,"arepo",color=cgColor(sPa.colorsA[cInd])
           
      legend,textoidl(['T_c = 5.5','T_{vir }(z=2)']),linestyle=lines,$;textcolor=colors,$
        linesize=0.3,position=[9.05,6.85],box=0,spacing=!p.charsize+0.5
        
      ; gadget histogram
      h2mt = h2G.h2_tmax_gal
      yrange = [h2G.mmTemp[0],(round(h2G.mmTemp[1])*10)/10]
      
      loadColorTable, 'helix', /reverse
      tvim,h2mt^exp,pcharsize=!p.charsize,scale=0,clip=-1,$
           xtitle="",ytitle="",$
           barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=h2A.mmMass,yrange=yrange,$
           xticks=3,xtickv=[9.0,10.0,11.0,12.0],position=[0.47,0.18,0.84,0.94],$
           ytickname=replicate(' ',10),/noerase
           
      ; colorbar
      barvals = findgen(ndivs+1)/ndivs*(max(h2mt^exp)-min(h2mt^exp)) + min(h2mt^exp)
      ticknames = textoidl(str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)')))
      ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
      colorbar,bottom=1,range=minmax(h2mt),position=[0.86,0.18,0.9,0.94],$
         /vertical,/right,title="",$
         divisions=ndivs,ticknames=ticknames,ncolors=255,charsize=!p.charsize-0.2
           
      ; temp lines
      tvir_vals = alog10(codeMassToVirTemp(10.0^h2G.mmMass/units.UnitMass_in_Msun,sP=sPa))
      
      cgPlot,h2G.mmMass,[Tc_val,Tc_val],line=lines[0],color=cgColor(colors[0]),/overplot
      cgPlot,h2G.mmMass,tvir_vals,line=lines[1],color=cgColor(colors[1]),/overplot
      cgText,h2G.mmMass[0]+0.3,Tc_val+0.2,"gadget",color=cgColor(sPg.colorsG[cInd])
           
      ; labels
      cgText,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),(0.1+0.8)/2,0.04,alignment=0.5,/normal
      cgText,textoidl("M_{gas,tr} [_{ }log h^{-1} M_{sun }]"),0.98,(0.18+0.94)/2,alignment=0.5,orientation=90,/normal
                
    end_PS
    
    ; plot (3) - arepo and gadget side by side (tvircur)
    start_PS, sPg.plotPath + 'temp2d_comp_tvircur_both.'+accMode+'.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
              str(sPg.res)+'_'+str(sPg.snap)+'.eps', xs=9.0, ys=4.0
          
      loadColorTable, 'helix', /reverse      
      ytitle = textoidl('log ( T_{max} / T_{vir,cur} )')
      
      ; arepo histogram
      h2mt = h2A.h2_tvircur_gal
      yrange = [h2A.mmRatio[0],(round(h2A.mmRatio[1])*10)/10]
      
      tvim,h2mt^exp,pcharsize=!p.charsize,scale=0,clip=-1,$
           xtitle="",ytitle=ytitle,$
           barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=h2A.mmMass,yrange=yrange,$
           xticks=2,xtickv=[9.0,10.0,11.0],position=[0.1,0.18,0.47,0.94]
           
      ; temp lines and legend
      tvir_vals = alog10(10.0^Tc_val / codeMassToVirTemp(10.0^h2A.mmMass/units.UnitMass_in_Msun,sP=sPa))
      
      cgPlot,h2A.mmMass,tvir_vals,line=lines[0],color=cgColor(colors[1]),/overplot
      cgText,h2A.mmMass[0]+0.3,yrange[0]+0.3,"arepo",color=cgColor(sPa.colorsA[cInd])
           
      ; gadget histogram
      h2mt = h2G.h2_tvircur_gal
      yrange = [h2G.mmRatio[0],(round(h2G.mmRatio[1])*10)/10]
      
      loadColorTable, 'helix', /reverse
      tvim,h2mt^exp,pcharsize=!p.charsize,scale=0,clip=-1,$
           xtitle="",ytitle="",$
           barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=h2A.mmMass,yrange=yrange,$
           xticks=3,xtickv=[9.0,10.0,11.0,12.0],position=[0.47,0.18,0.84,0.94],$
           ytickname=replicate(' ',10),/noerase
           
      ; colorbar
      barvals = findgen(ndivs+1)/ndivs*(max(h2mt^exp)-min(h2mt^exp)) + min(h2mt^exp)
      ticknames = textoidl(str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)')))
      ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
      colorbar,bottom=1,range=minmax(h2mt),position=[0.86,0.18,0.9,0.94],$
         /vertical,/right,title="",$
         divisions=ndivs,ticknames=ticknames,ncolors=255,charsize=!p.charsize-0.2
           
      ; temp lines
      cgPlot,h2G.mmMass,tvir_vals,line=lines[0],color=cgColor(colors[1]),/overplot
      cgText,h2G.mmMass[0]+0.3,yrange[0]+0.3,"gadget",color=cgColor(sPg.colorsG[cInd])
           
      ; labels
      cgText,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),(0.1+0.8)/2,0.04,alignment=0.5,/normal
      cgText,textoidl("M_{gas,tr} [_{ }log h^{-1} M_{sun }]"),0.98,(0.18+0.94)/2,alignment=0.5,orientation=90,/normal
                
    end_PS
    
    ; plot (4) - arepo and gadget side by side (tviracc)
    start_PS, sPg.plotPath + 'temp2d_comp_tviracc_both.'+accMode+'.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
              str(sPg.res)+'_'+str(sPg.snap)+'.eps', xs=9.0, ys=4.0
          
      loadColorTable, 'helix', /reverse      
      ytitle = textoidl('log ( T_{max} / T_{vir,acc} )')
      
      ; arepo histogram
      h2mt = h2A.h2_tviracc_gal
      yrange = [h2A.mmRatio[0],(round(h2A.mmRatio[1])*10)/10]
      
      tvim,h2mt^exp,pcharsize=!p.charsize,scale=0,clip=-1,$
           xtitle="",ytitle=ytitle,$
           barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=h2A.mmMass,yrange=yrange,$
           xticks=2,xtickv=[9.0,10.0,11.0],position=[0.1,0.18,0.47,0.94]
           
      ; temp lines and legend
      tvir_vals = alog10(10.0^Tc_val / codeMassToVirTemp(10.0^h2A.mmMass/units.UnitMass_in_Msun,sP=sPa))
      
      cgPlot,h2A.mmMass,tvir_vals,line=lines[0],color=cgColor(colors[1]),/overplot
      cgText,h2A.mmMass[0]+0.3,yrange[0]+0.3,"arepo",color=cgColor(sPa.colorsA[cInd])
           
      ; gadget histogram
      h2mt = h2G.h2_tviracc_gal
      yrange = [h2G.mmRatio[0],(round(h2G.mmRatio[1])*10)/10]
      
      loadColorTable, 'helix', /reverse
      tvim,h2mt^exp,pcharsize=!p.charsize,scale=0,clip=-1,$
           xtitle="",ytitle="",$
           barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=h2A.mmMass,yrange=yrange,$
           xticks=3,xtickv=[9.0,10.0,11.0,12.0],position=[0.47,0.18,0.84,0.94],$
           ytickname=replicate(' ',10),/noerase
           
      ; colorbar
      barvals = findgen(ndivs+1)/ndivs*(max(h2mt^exp)-min(h2mt^exp)) + min(h2mt^exp)
      ticknames = textoidl(str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)')))
      ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
      colorbar,bottom=1,range=minmax(h2mt),position=[0.86,0.18,0.9,0.94],$
         /vertical,/right,title="",$
         divisions=ndivs,ticknames=ticknames,ncolors=255,charsize=!p.charsize-0.2
           
      ; temp lines
      cgPlot,h2G.mmMass,tvir_vals,line=lines[0],color=cgColor(colors[1]),/overplot
      cgText,h2G.mmMass[0]+0.3,yrange[0]+0.3,"gadget",color=cgColor(sPg.colorsG[cInd])
           
      ; labels
      cgText,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),(0.1+0.8)/2,0.04,alignment=0.5,/normal
      cgText,textoidl("M_{gas,tr} [_{ }log h^{-1} M_{sun }]"),0.98,(0.18+0.94)/2,alignment=0.5,orientation=90,/normal
                
    end_PS
    
  endforeach
  
  !except = 1
              stop
end

; -------------------------------------------------------------------------------------------------

pro plotTest

  resolutions = [128,256,512]
  lines       = [1,2,0]
  sgSelect    = 'pri'
  accMode     = 'all'

  start_PS, '/n/home07/dnelson/test.eps'
    foreach res,resolutions,j do begin
      print,res
      
      ; load
      sPg = simParams(res=res,run='gadget',redshift=2.0)
      parentMass_ga = gcSubsetProp(sP=sPg,select=sgSelect,/parMass,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
      
      ; histogram
      bin = 0.1
      hist = histogram([parentMass_ga.gal,parentMass_ga.stars],bin=bin,loc=loc)
      hist *= sPg.targetGasMass
      
      ; plot
      overplot = 0
      if j gt 0 then overplot = 1
      cgPlot,loc+bin/2.0,hist,psym=10,linestyle=lines[j],overplot=overplot,$
        yrange=[1e-1,1e3],/ylog,xrange=[8.0,12.5],xtitle="Parent Halo Mass",ytitle="Total Galaxy Mass per 0.1dex"
      
    endforeach ;res
    
    ; legend
    legend,textoidl(['128^3','256^3','512^3']),linestyle=[lines[0],lines[1],lines[2]],box=0,linesize=0.25,/top,/right

  end_PS

  stop
end

pro plotTest2

  ; config
  sPg = simParams(res=256,run='gadget',redshift=2.0)
  sgSelect = 'pri'
  accMode  = 'all'
  binRad   = 0.005
  binDens  = 0.05
  virInd   = 0

  ; to decide hot/cold
  accTvir = gcSubsetProp(sP=sPg,select=sgSelect,/accTvir,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  maxTemp = gcSubsetProp(sP=sPg,select=sgSelect,/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  
  ; radial distribution at z=2 and density
  radPri = gcSubsetProp(sP=sPg,select=sgSelect,/rVirNorm,parNorm='pri',/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  dens   = gcSubsetProp(sP=sPg,select=sgSelect,/curSingleVal,singleValField='dens',/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)

  ; select cold/hot
  wCold = where(10.0^maxTemp.gal / 10.0^accTvir.gal le sPg.TvirVals[virInd],countCold,comp=wHot,ncomp=countHot)

  start_PS, sPg.plotPath + 'raddist.'+accMode+'.'+sPg.plotPrefix+'.'+$
            str(sPg.res)+'_'+str(sPg.snap)+'.eps'
            
    cgPlot,[0],[0],/nodata,xrange=[0.0,0.2],yrange=[0.0,1.1],xtitle="r/rvir",ytitle="N"
            
    hist = histogram(radPri.gal[wCold],bin=binRad,min=0.0,loc=loc)
    cgPlot,loc+binRad*0.5,hist/float(max(hist)),psym=10,color=cgColor('blue'),/overplot
    
    hist = histogram(radPri.gal[wHot],bin=binRad,min=0.0,loc=loc)
    cgPlot,loc+binRad*0.5,hist/float(max(hist)),psym=10,color=cgColor('red'),/overplot
    
    hist = histogram(radPri.gal,bin=binRad,min=0.0,loc=loc)        
    cgPlot,loc+binRad*0.5,hist/float(max(hist)),psym=10,linestyle=1,color=cgColor('gray'),/overplot
    
  end_PS
  
  start_PS, sPg.plotPath + 'densdist.'+accMode+'.'+sPg.plotPrefix+'.'+$
            str(sPg.res)+'_'+str(sPg.snap)+'.eps'
            
    cgPlot,[0],[0],/nodata,xrange=[-6.5,-2.5],yrange=[0.0,1.1],xtitle="r/rvir",ytitle="N"
            
    dens.gal = alog10(dens.gal)
    
    hist = histogram(dens.gal[wCold],bin=binDens,loc=loc)
    cgPlot,loc+binDens*0.5,hist/float(max(hist)),psym=10,color=cgColor('blue'),/overplot
    
    hist = histogram(dens.gal[wHot],bin=binDens,loc=loc)
    cgPlot,loc+binDens*0.5,hist/float(max(hist)),psym=10,color=cgColor('red'),/overplot
    
    hist = histogram(dens.gal,bin=binDens,loc=loc)        
    cgPlot,loc+binDens*0.5,hist/float(max(hist)),psym=10,linestyle=1,color=cgColor('gray'),/overplot
    
  end_PS
  
  stop
  
end