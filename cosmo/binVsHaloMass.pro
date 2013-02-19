; plotVsHaloMass.pro
; gas accretion project - bin quantities as a function of halo mass
; dnelson nov.2012

; haloMassBinValues(): bin accretion rate (in cold/hot) and cold fraction as a function of halo mass

function haloMassBinValues, sP=sP, sgSelect=sgSelect, accMode=accMode, timeWindow=TW

  compile_opt idl2, hidden, strictarr, strictarrsubs
  if ~sP.gfmWinds and accMode eq 'recycled' then message,'Error: Request recycled on non-winds run.'
  units = getUnits()
  
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
    if isnumeric(TW) then begin
      timeWindow = TW * 1e6 ; convert input Myr to yr
    endif else begin
      timeWindow = TW ; "tVir_tIGM" or "tVir_tIGM_bin"
    endelse
  endelse

  ; config
  nCuts = n_elements(sP.TcutVals)
  nVirs = n_elements(sP.TvirVals)
  minNumFrac = 6
  
  if sP.trMCPerCell le 0 then massPerPart = sP.targetGasMass ; SPH or vel tracer
  if sP.trMCPerCell gt 0 then massPerPart = sP.trMassConst ; MC tracer
  
  ; check if save exists
  saveFilename = sP.derivPath + 'binnedVals/binVals.' + sP.saveTag + '.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.cut' + str(nCuts) + '.vir' + str(nVirs) + '.' + sgSelect + '.' + accMode + $
    '.r' + str(sP.radIndHaloAcc) + '.r' + str(sP.radIndGalAcc) + '_tw' + str(TW) + '.sav'
  
  ; results exist, return
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif   
  
  ; make a uniform gas selection at the start
  at = accretionTimes(sP=sP)
  mt = mergerTreeSubset(sP=sP)
  
  wAm = accModeInds(at=at,accMode=accMode,sP=sP,/mask)
    
  ; reverse histogram parent IDs of all particles/tracers in this selection
  gcIndOrig = mergerTreeRepParentIDs(mt=mt,sP=sP,/compactMtS)
  
  maxHist = max([gcIndOrig.gal[wAm.gal],gcIndOrig.gmem[wAm.gmem],gcIndOrig.stars[wAm.stars]])
  hist_gal   = histogram(gcIndOrig.gal[wAm.gal],min=0,max=maxHist,loc=loc_gal,rev=rev_gal)
  hist_gmem  = histogram(gcIndOrig.gmem[wAm.gmem],min=0,max=maxHist,loc=loc_gmem,rev=rev_gmem)
  hist_stars = histogram(gcIndOrig.stars[wAm.stars],min=0,max=maxHist,loc=loc_stars,rev=rev_stars)
  gcIndOrig = !NULL
  
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
  coldAccRate = { gal_const     : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                  gmem_const    : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                  stars_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                  gal_tvircur   : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                  gmem_tvircur  : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                  stars_tvircur : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                  gal_tviracc   : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                  gmem_tviracc  : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                  stars_tviracc : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                  both_const    : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$ ; both=gal+stars
                  both_tvircur  : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                  both_tviracc  : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                  total_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$ ; total=gal+stars+gmem
                  total_tvircur : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                  total_tviracc : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                  gal_num       : lonarr(n_elements(mt.galcatIDList))         ,$
                  gmem_num      : lonarr(n_elements(mt.galcatIDList))         ,$
                  stars_num     : lonarr(n_elements(mt.galcatIDList))          }
                  
  hotAccRate = coldAccRate
  coldFrac   = coldAccRate
  totalColdMass = coldAccRate
  totalHotMass  = coldAccRate
              
  coldFrac_cur = { gal_const     : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                   gmem_const    : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                   stars_const   : -1   ,$  ; leave empty for Tcur so we don't use it by accident
                   gal_tvircur   : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                   gmem_tvircur  : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                   stars_tvircur : -1   ,$
                   gal_tviracc   : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                   gmem_tviracc  : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                   stars_tviracc : -1   ,$
                   both_const    : -1   ,$
                   both_tvircur  : -1   ,$
                   both_tviracc  : -1   ,$
                   total_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                   total_tvircur : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                   total_tviracc : fltarr(nVirs,n_elements(mt.galcatIDList))   ,$
                   gal_num       : lonarr(n_elements(mt.galcatIDList))         ,$
                   gmem_num      : lonarr(n_elements(mt.galcatIDList))         ,$
                   stars_num     : -1                                           }              
               
  ; loop over all tracked subgroups
  for i=0L,maxHist do begin
  
    ; calculate min(timeWindow,length of time this halo was tracked back) to use to normalize rates
    loc_timeWindow = timeWindow

    ; special types of timeWindows
    if ~isnumeric(loc_timeWindow) then begin
      temp_IGM     = 4e4 ; K (based on FG11)
      massRanges   = [10.0,10.7,11.5] ; log Msun
      massRangeTWs = [500.0,1500.0,2000.0,2500.0]*1e6 ; yr
      
      ; 1. special timewindow: until tvir drops below tIGM
      if loc_timeWindow eq 'tVir_tIGM' then begin
        temp_ratios = 10.0^mt.hVirTemp[*,i] / temp_IGM
	ww = where(temp_ratios ge 1.0,count_loc)
	
	if count_loc eq 0 then begin
	  ; low mass
	  loc_timeWindow = massRangeTWs[0] ; 500Myr
	endif else begin
	  ; set snapshot -> timeWindow
	  endSnapLoc = mt.maxSnap - max(ww)
	  endAgeLoc = snapNumToAgeFlat(snap=endSnapLoc,sP=sP)
	  loc_timeWindow = curtime - endAgeLoc*1e9
	endelse
      endif else begin
    
        ; 2. special timewindow: until tvir of mass bin drops below tIGM
        if loc_timeWindow eq 'tVir_tIGM_bin' then begin
          codeMassLoc = codeMassToLogMsun(mt.hMass[0,i])
          if codeMassLoc lt massRanges[0] then $
	    loc_timeWindow = massRangeTWs[0]
	  if codeMassLoc ge massRanges[0] and codeMassLoc lt massRanges[1] then $
	    loc_timeWindow = massRangeTWs[1]
	  if codeMassLoc ge massRanges[1] and codeMassLoc lt massRanges[2] then $
	    loc_timeWindow = massRangeTWs[2]
	  if codeMassLoc ge massRanges[2] then $
	    loc_timeWindow = massRangeTWs[3]
        endif ;tVir_tIGM_bin
      endelse ;tVir_tIGM
    endif ;isnumeric
    
    ; if current halo is tracked less than timewindow, shorten timewindow to tracked period
    if mt.hMinSnap[i] ne -1 then begin
      loc_minSnapTime = snapTimes[mt.hMinSnap[i]]
      if curtime-loc_minSnapTime lt loc_timeWindow then loc_timeWindow = curtime-loc_minSnapTime
    endif
  
    if hist_gal[i] gt 0 then begin
      ; list of indices of galaxy gas particles in this subgroup
      loc_inds_gal = rev_gal[rev_gal[i]:rev_gal[i+1]-1]
      
      ; debug
      ;loc_inds_gal2 = where(gcIndOrigTr.gal[gal_w] eq i)
      ;if ~array_equal(loc_inds_gal,loc_inds_gal2) then print,'error'
      
      ; corresponding accretion times for these particles
      ; OLD: galaxy accretion defined as rvir (or any other radius) crossing time
      ;loc_atime_gal = reform(at.accTime_gal[sP.radIndGalAcc,wAm.gal[loc_inds_gal]])
      
      ; NEW: galaxy accretion defined as (rho,temp) joining time or 0.15rvir crossing time (most recent)
      loc_atime_gal = reform(at.accTimeRT_gal[wAm.gal[loc_inds_gal]])
      
      r_crossing_time = reform(at.accTime_gal[sP.radIndGalAcc,wAm.gal[loc_inds_gal]])
      w = where(r_crossing_time gt loc_atime_gal,count)
      if count gt 0 then loc_atime_gal[w] = r_crossing_time[w]
      
      loc_atime_gal = 1/loc_atime_gal - 1 ; redshift
      loc_atime_gal = redshiftToAgeFlat(loc_atime_gal)*1e9 ; yr
      
      ; make a count of those falling in the time window
      w = where(curtime - loc_atime_gal le loc_timeWindow,nloc)
      
      coldFrac.gal_num[i]    = nloc
      coldAccRate.gal_num[i] = nloc
      
      ; maximum past temps, cur and acc tvirs for only those particles in the time window
      if nloc gt 0 then begin
        loc_maxt_gal = maxTemp.gal[loc_inds_gal[w]]
        loc_curtvir_gal = curTvir.gal[loc_inds_gal[w]]
        loc_acctvir_gal = accTvir.gal[loc_inds_gal[w]]
  
        ; count mass elements with Tmax below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_maxt_gal le sP.TcutVals[j],count_cold,ncomp=count_hot)
          coldFrac.gal_const[j,i]    = float(count_cold) / nloc
          coldAccRate.gal_const[j,i] = count_cold
          hotAccRate.gal_const[j,i]  = count_hot
        endfor
        
        for j=0,nVirs-1 do begin
          ; count mass elements with Tmax below Tvir at current time
          w = where(10.0^loc_maxt_gal / 10.0^loc_curtvir_gal le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldFrac.gal_tvircur[j,i]    = float(count_cold) / nloc
          coldAccRate.gal_tvircur[j,i] = count_cold
          hotAccRate.gal_tvircur[j,i]  = count_hot
    
          ; count mass elements with Tmax below Tvir at accretion time
          w = where(10.0^loc_maxt_gal / 10.0^loc_acctvir_gal le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldFrac.gal_tviracc[j,i]    = float(count_cold) / nloc
          coldAccRate.gal_tviracc[j,i] = count_cold
          hotAccRate.gal_tviracc[j,i]  = count_hot
        endfor
      endif ;nloc>0
      
      ; current temp block (for cold fractions)
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
          w = where(loc_curt_gal le sP.TcutVals[j],count_below)
          coldFrac_cur.gal_const[j,i] = float(count_below) / count_nosfr
        endfor
        
        for j=0,nVirs-1 do begin
          ; count fraction Tmax below Tvir at current time
          w = where(10.0^loc_curt_gal / 10.0^curTvir.gal[loc_inds_gal] le sP.TvirVals[j],count_below)
          coldFrac_cur.gal_tvircur[j,i] = float(count_below) / count_nosfr
    
          ; count fraction Tmax below Tvir at accretion time
          w = where(10.0^loc_curt_gal / 10.0^accTvir.gal[loc_inds_gal] le sP.TvirVals[j],count_below)
          coldFrac_cur.gal_tviracc[j,i] = float(count_below) / count_nosfr
        endfor
      endif ; cursfr!=0
      
    endif ; hist_gal[i]

    if hist_gmem[i] gt 0 then begin
      ; list of indices of group member gas particles in this subgroup
      loc_inds_gmem = rev_gmem[rev_gmem[i]:rev_gmem[i+1]-1]
      
      ; corresponding accretion times for these particles
      loc_atime_gmem = reform(at.accTime_gmem[sP.radIndHaloAcc,wAm.gmem[loc_inds_gmem]])
      loc_atime_gmem = 1/loc_atime_gmem - 1 ; redshift
      loc_atime_gmem = redshiftToAgeFlat(loc_atime_gmem)*1e9 ; yr
      
      ; make a count of those falling in the time window
      w = where(curtime - loc_atime_gmem le loc_timeWindow,nloc)
      
      coldFrac.gmem_num[i]    = nloc
      coldAccRate.gmem_num[i] = nloc
      
      ; maximum past temps, cur and acc tvirs for only those particles in the time window
      if nloc gt 0 then begin
        loc_maxt_gmem = maxTemp.gmem[loc_inds_gmem[w]]
        loc_curtvir_gmem = curTvir.gmem[loc_inds_gmem[w]]
        loc_acctvir_gmem = accTvir.gmem[loc_inds_gmem[w]]
  
        ; count mass elements with Tmax below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_maxt_gmem le sP.TcutVals[j],count_cold,ncomp=count_hot)
          coldFrac.gmem_const[j,i]    = float(count_cold) / nloc
          coldAccRate.gmem_const[j,i] = count_cold
          hotAccRate.gmem_const[j,i]  = count_hot
        endfor
        
        for j=0,nVirs-1 do begin
          ; count mass elements with Tmax below Tvir at current time
          w = where(10.0^loc_maxt_gmem / 10.0^loc_curtvir_gmem le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldFrac.gmem_tvircur[j,i]    = float(count_cold) / nloc
          coldAccRate.gmem_tvircur[j,i] = count_cold
          hotAccRate.gmem_tvircur[j,i]  = count_hot
    
          ; count mass elements with Tmax below Tvir at accretion time
          w = where(10.0^loc_maxt_gmem / 10.0^loc_acctvir_gmem le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldFrac.gmem_tviracc[j,i]    = float(count_cold) / nloc
          coldAccRate.gmem_tviracc[j,i] = count_cold
          hotAccRate.gmem_tviracc[j,i]  = count_hot
        endfor
      endif ; nloc>0
      
      ; current temp (for cold fractions)
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
          w = where(loc_curt_gmem le sP.TcutVals[j],count_below)
          coldFrac_cur.gmem_const[j,i] = float(count_below) / count_nosfr
        endfor
        
        for j=0,nVirs-1 do begin
          ; count fraction Tmax below Tvir at current time
          w = where(10.0^loc_curt_gmem / 10.0^curTvir.gmem[loc_inds_gmem] le sP.TvirVals[j],count_below)
          coldFrac_cur.gmem_tvircur[j,i] = float(count_below) / count_nosfr
    
          ; count fraction Tmax below Tvir at accretion time
          w = where(10.0^loc_curt_gmem / 10.0^accTvir.gmem[loc_inds_gmem] le sP.TvirVals[j],count_below)
          coldFrac_cur.gmem_tviracc[j,i] = float(count_below) / count_nosfr
        endfor
      endif ; cursfr!=0
      
    endif ; hist_gmem[i]

    if hist_stars[i] gt 0 then begin
      ; list of indices of star gas particles in this subgroup
      loc_inds_stars = rev_stars[rev_stars[i]:rev_stars[i+1]-1]
      
      ; corresponding accretion times for these particles:
      ; stellar accretion defined as (rho,temp) joining time or 0.15rvir crossing time (most recent)
      loc_atime_stars = reform(at.accTimeRT_stars[wAm.stars[loc_inds_stars]])
      
      r_crossing_time = reform(at.accTime_stars[sP.radIndGalAcc,wAm.stars[loc_inds_stars]])
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
      w = where(curtime - loc_atime_stars le loc_timeWindow,nloc)
      
      coldFrac.stars_num[i]    = nloc
      coldAccRate.stars_num[i] = nloc
      
      ; maximum past temps, cur and acc tvirs for only those particles in the time window
      if nloc gt 0 then begin
        loc_maxt_stars = maxTemp.stars[loc_inds_stars[w]]
        loc_curtvir_stars = curTvir.stars[loc_inds_stars[w]]
        loc_acctvir_stars = accTvir.stars[loc_inds_stars[w]]
  
        ; count mass elements with Tmax below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_maxt_stars le sP.TcutVals[j],count_cold,ncomp=count_hot)
          coldFrac.stars_const[j,i]    = float(count_cold) / nloc
          coldAccRate.stars_const[j,i] = count_cold
          hotAccRate.stars_const[j,i]  = count_hot
        endfor
        
        for j=0,nVirs-1 do begin
          ; count mass elements with Tmax below Tvir at current time
          w = where(10.0^loc_maxt_stars / 10.0^loc_curtvir_stars le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldFrac.stars_tvircur[j,i]    = float(count_cold) / nloc
          coldAccRate.stars_tvircur[j,i] = count_cold
          hotAccRate.stars_tvircur[j,i]  = count_hot
    
          ; count mass elements with Tmax below Tvir at accretion time
          w = where(10.0^loc_maxt_stars / 10.0^loc_acctvir_stars le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldFrac.stars_tviracc[j,i]    = float(count_cold) / nloc
          coldAccRate.stars_tviracc[j,i] = count_cold
          hotAccRate.stars_tviracc[j,i]  = count_hot
        endfor
      endif ;nloc>0
    endif ; hist_stars[i]
    
    ; convert accretion total mass (counts) to msun
    totalColdMass.gal_const[*,i]     = coldAccRate.gal_const[*,i]     * massPerPart * units.UnitMass_in_Msun
    totalColdMass.gmem_const[*,i]    = coldAccRate.gmem_const[*,i]    * massPerPart * units.UnitMass_in_Msun
    totalColdMass.stars_const[*,i]   = coldAccRate.stars_const[*,i]   * massPerPart * units.UnitMass_in_Msun
    totalColdMass.gal_tvircur[*,i]   = coldAccRate.gal_tvircur[*,i]   * massPerPart * units.UnitMass_in_Msun
    totalColdMass.gmem_tvircur[*,i]  = coldAccRate.gmem_tvircur[*,i]  * massPerPart * units.UnitMass_in_Msun
    totalColdMass.stars_tvircur[*,i] = coldAccRate.stars_tvircur[*,i] * massPerPart * units.UnitMass_in_Msun
    totalColdMass.gal_tviracc[*,i]   = coldAccRate.gal_tviracc[*,i]   * massPerPart * units.UnitMass_in_Msun
    totalColdMass.gmem_tviracc[*,i]  = coldAccRate.gmem_tviracc[*,i]  * massPerPart * units.UnitMass_in_Msun
    totalColdMass.stars_tviracc[*,i] = coldAccRate.stars_tviracc[*,i] * massPerPart * units.UnitMass_in_Msun
    
    totalHotMass.gal_const[*,i]     = hotAccRate.gal_const[*,i]     * massPerPart * units.UnitMass_in_Msun
    totalHotMass.gmem_const[*,i]    = hotAccRate.gmem_const[*,i]    * massPerPart * units.UnitMass_in_Msun
    totalHotMass.stars_const[*,i]   = hotAccRate.stars_const[*,i]   * massPerPart * units.UnitMass_in_Msun
    totalHotMass.gal_tvircur[*,i]   = hotAccRate.gal_tvircur[*,i]   * massPerPart * units.UnitMass_in_Msun
    totalHotMass.gmem_tvircur[*,i]  = hotAccRate.gmem_tvircur[*,i]  * massPerPart * units.UnitMass_in_Msun
    totalHotMass.stars_tvircur[*,i] = hotAccRate.stars_tvircur[*,i] * massPerPart * units.UnitMass_in_Msun
    totalHotMass.gal_tviracc[*,i]   = hotAccRate.gal_tviracc[*,i]   * massPerPart * units.UnitMass_in_Msun
    totalHotMass.gmem_tviracc[*,i]  = hotAccRate.gmem_tviracc[*,i]  * massPerPart * units.UnitMass_in_Msun
    totalHotMass.stars_tviracc[*,i] = hotAccRate.stars_tviracc[*,i] * massPerPart * units.UnitMass_in_Msun
    
    ; convert accretion total(counts) to msun/year
    coldAccRate.gal_const[*,i]     *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    coldAccRate.gmem_const[*,i]    *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    coldAccRate.stars_const[*,i]   *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    coldAccRate.gal_tvircur[*,i]   *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    coldAccRate.gmem_tvircur[*,i]  *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    coldAccRate.stars_tvircur[*,i] *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    coldAccRate.gal_tviracc[*,i]   *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    coldAccRate.gmem_tviracc[*,i]  *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    coldAccRate.stars_tviracc[*,i] *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    
    hotAccRate.gal_const[*,i]     *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    hotAccRate.gmem_const[*,i]    *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    hotAccRate.stars_const[*,i]   *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    hotAccRate.gal_tvircur[*,i]   *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    hotAccRate.gmem_tvircur[*,i]  *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    hotAccRate.stars_tvircur[*,i] *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    hotAccRate.gal_tviracc[*,i]   *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    hotAccRate.gmem_tviracc[*,i]  *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    hotAccRate.stars_tviracc[*,i] *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    
  endfor ; i
  
  ; create totals of both=gal+stars, total=both+gmem
  for j=0,nVirs-1 do begin
    coldAccRate.both_tvircur[j,*]  = coldAccRate.gal_tvircur[j,*] + coldAccRate.stars_tvircur[j,*]
    coldAccRate.both_tviracc[j,*]  = coldAccRate.gal_tviracc[j,*] + coldAccRate.stars_tviracc[j,*]
    coldAccRate.total_tvircur[j,*] = coldAccRate.both_tvircur[j,*] + coldAccRate.gmem_tvircur[j,*]
    coldAccRate.total_tviracc[j,*] = coldAccRate.both_tviracc[j,*] + coldAccRate.gmem_tviracc[j,*]
    
    hotAccRate.both_tvircur[j,*]  = hotAccRate.gal_tvircur[j,*] + hotAccRate.stars_tvircur[j,*]
    hotAccRate.both_tviracc[j,*]  = hotAccRate.gal_tviracc[j,*] + hotAccRate.stars_tviracc[j,*]
    hotAccRate.total_tvircur[j,*] = hotAccRate.both_tvircur[j,*] + hotAccRate.gmem_tvircur[j,*]
    hotAccRate.total_tviracc[j,*] = hotAccRate.both_tviracc[j,*] + hotAccRate.gmem_tviracc[j,*]
    
    totalColdMass.both_tvircur[j,*]  = totalColdMass.gal_tvircur[j,*] + totalColdMass.stars_tvircur[j,*]
    totalColdMass.both_tviracc[j,*]  = totalColdMass.gal_tviracc[j,*] + totalColdMass.stars_tviracc[j,*]
    totalColdMass.total_tvircur[j,*] = totalColdMass.both_tvircur[j,*] + totalColdMass.gmem_tvircur[j,*]
    totalColdMass.total_tviracc[j,*] = totalColdMass.both_tviracc[j,*] + totalColdMass.gmem_tviracc[j,*]
    
    totalHotMass.both_tvircur[j,*]  = totalHotMass.gal_tvircur[j,*] + totalHotMass.stars_tvircur[j,*]
    totalHotMass.both_tviracc[j,*]  = totalHotMass.gal_tviracc[j,*] + totalHotMass.stars_tviracc[j,*]
    totalHotMass.total_tvircur[j,*] = totalHotMass.both_tvircur[j,*] + totalHotMass.gmem_tvircur[j,*]
    totalHotMass.total_tviracc[j,*] = totalHotMass.both_tviracc[j,*] + totalHotMass.gmem_tviracc[j,*]
    
    coldFrac.both_tvircur[j,*]  = (coldFrac.gal_tvircur[j,*] * coldFrac.gal_num + coldFrac.stars_tvircur[j,*] * coldFrac.stars_num) / (coldFrac.gal_num + coldFrac.stars_num)
    coldFrac.total_tvircur[j,*] = (coldFrac.gal_tvircur[j,*] * coldFrac.gal_num + coldFrac.gmem_tvircur[j,*] * coldFrac.gmem_num + coldFrac.stars_tvircur[j,*] * coldFrac.stars_num) / (coldFrac.gal_num + coldFrac.stars_num + coldFrac.gmem_num)                
    coldFrac.both_tviracc[j,*]  = (coldFrac.gal_tviracc[j,*] * coldFrac.gal_num + coldFrac.stars_tviracc[j,*] * coldFrac.stars_num) / (coldFrac.gal_num + coldFrac.stars_num)
    coldFrac.total_tviracc[j,*] = (coldFrac.gal_tviracc[j,*] * coldFrac.gal_num + coldFrac.gmem_tviracc[j,*] * coldFrac.gmem_num + coldFrac.stars_tviracc[j,*] * coldFrac.stars_num) / (coldFrac.gal_num + coldFrac.stars_num + coldFrac.gmem_num)  
    
    coldFrac_cur.total_tvircur[j,*] = (coldFrac_cur.gal_tvircur[j,*] * coldFrac_cur.gal_num + coldFrac_cur.gmem_tvircur[j,*] * coldFrac_cur.gmem_num) / (coldFrac_cur.gal_num + coldFrac_cur.gmem_num)              
    coldFrac_cur.total_tviracc[j,*] = (coldFrac_cur.gal_tviracc[j,*] * coldFrac_cur.gal_num + coldFrac_cur.gmem_tviracc[j,*] * coldFrac_cur.gmem_num) / (coldFrac_cur.gal_num + coldFrac_cur.gmem_num)  
  endfor
  
  for j=0,nCuts-1 do begin
    coldAccRate.both_const[j,*]  = coldAccRate.gal_const[j,*] + coldAccRate.stars_const[j,*]
    coldAccRate.total_const[j,*] = coldAccRate.both_const[j,*] + coldAccRate.gmem_const[j,*]
    hotAccRate.both_const[j,*]   = hotAccRate.gal_const[j,*] + hotAccRate.stars_const[j,*]
    hotAccRate.total_const[j,*]  = hotAccRate.both_const[j,*] + hotAccRate.gmem_const[j,*]
    
    totalColdMass.both_const[j,*]  = totalColdMass.gal_const[j,*] + totalColdMass.stars_const[j,*]
    totalColdMass.total_const[j,*] = totalColdMass.both_const[j,*] + totalColdMass.gmem_const[j,*]
    totalHotMass.both_const[j,*]   = totalHotMass.gal_const[j,*] + totalHotMass.stars_const[j,*]
    totalHotMass.total_const[j,*]  = totalHotMass.both_const[j,*] + totalHotMass.gmem_const[j,*]
    
    coldFrac.both_const[j,*]  = (coldFrac.gal_const[j,*] * coldFrac.gal_num + coldFrac.stars_const[j,*] * coldFrac.stars_num) / (coldFrac.gal_num + coldFrac.stars_num)  
    coldFrac.total_const[j,*] = (coldFrac.gal_const[j,*] * coldFrac.gal_num + coldFrac.gmem_const[j,*] * coldFrac.gmem_num + coldFrac.stars_const[j,*] * coldFrac.stars_num) / (coldFrac.gal_num + coldFrac.stars_num + coldFrac.gmem_num)      
    coldFrac_cur.total_const[j,*] = (coldFrac_cur.gal_const[j,*] * coldFrac_cur.gal_num + coldFrac_cur.gmem_const[j,*] * coldFrac_cur.gmem_num) / (coldFrac_cur.gal_num + coldFrac_cur.gmem_num)  
  endfor
  
  ; bin fractions into halo mass bins and make median lines
  logMassBins   = [9.5,10.0,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11.0,$
                   11.1,11.25,11.5,11.75,11.9,13.1]
  logMassNBins  = n_elements(logMassBins)-1
  logMassBinCen = 0.5 * (logMassBins + shift(logMassBins,-1))
  logMassBinCen = logMassBinCen[0:-2]
  
  ; structures to store the binned values
  coldMedian = { gal_const     : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 stars_const   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 gmem_const    : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 both_const    : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 total_const   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 gal_tVircur   : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 stars_tVircur : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 gmem_tVircur  : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 both_tVircur  : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 total_tVircur : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 gal_tViracc   : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 stars_tViracc : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 gmem_tViracc  : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 both_tViracc  : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
                 total_tViracc : fltarr(nVirs,logMassNbins) + !values.f_nan  }
                 
  hotMedian      = coldMedian
  fracMedian     = coldMedian
  fracMedian_cur = coldMedian
  coldTotal      = coldMedian ; total mass (not rates)
  hotTotal       = coldMedian ; total mass (not rates)
  
  ; total of hot+cold (separated by gal, gmem, stars, both=gal+stars, tot=gal+stars+gmem)
  totalHCMedian = { gal    : fltarr(logMassNbins) + !values.f_nan ,$
                    gmem   : fltarr(logMassNbins) + !values.f_nan ,$
                    stars  : fltarr(logMassNbins) + !values.f_nan ,$
                    both   : fltarr(logMassNbins) + !values.f_nan ,$
                    tot    : fltarr(logMassNbins) + !values.f_nan  }
                    
  totalMassHC = { gal    : fltarr(logMassNbins) + !values.f_nan ,$
                  gmem   : fltarr(logMassNbins) + !values.f_nan ,$
                  stars  : fltarr(logMassNbins) + !values.f_nan ,$
                  both   : fltarr(logMassNbins) + !values.f_nan ,$
                  tot    : fltarr(logMassNbins) + !values.f_nan  }
                 
  ; calculate median values in bins of halo mass
  for i=0,logMassNbins-1 do begin

    ; --- accretion rate ---
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1],count)
    
    if count gt 0 then begin
      for j=0,nVirs-1 do begin
        ; accretion rates:
        ; gal (hot+cold)
        coldMedian.gal_tVirCur[j,i]   = median(coldAccRate.gal_tvircur[j,w])
        coldMedian.gal_tVirAcc[j,i]   = median(coldAccRate.gal_tviracc[j,w])
        hotMedian.gal_tVirCur[j,i]    = median(hotAccRate.gal_tvircur[j,w])
        hotMedian.gal_tVirAcc[j,i]    = median(hotAccRate.gal_tviracc[j,w])
        
        ; gmem (hot+cold)
        coldMedian.gmem_tVirCur[j,i]  = median(coldAccRate.gmem_tvircur[j,w])
        coldMedian.gmem_tVirAcc[j,i]  = median(coldAccRate.gmem_tviracc[j,w])
        hotMedian.gmem_tVirCur[j,i]   = median(hotAccRate.gmem_tvircur[j,w])
        hotMedian.gmem_tVirAcc[j,i]   = median(hotAccRate.gmem_tviracc[j,w]) 
        
        ; stars (hot+cold)
        coldMedian.stars_tVirCur[j,i] = median(coldAccRate.stars_tvircur[j,w])
        coldMedian.stars_tVirAcc[j,i] = median(coldAccRate.stars_tviracc[j,w])
        hotMedian.stars_tVirCur[j,i]  = median(hotAccRate.stars_tvircur[j,w])
        hotMedian.stars_tVirAcc[j,i]  = median(hotAccRate.stars_tviracc[j,w]) 
        
        ; both=gal+stars (hot+cold)
        coldMedian.both_tVirCur[j,i]  = median(coldAccRate.both_tvircur[j,w])
        coldMedian.both_tVirAcc[j,i]  = median(coldAccRate.both_tviracc[j,w])
        hotMedian.both_tVirCur[j,i]   = median(hotAccRate.both_tvircur[j,w])
        hotMedian.both_tVirAcc[j,i]   = median(hotAccRate.both_tviracc[j,w])
        
        ; total=gal+stars+gmem (hot+cold)
        coldMedian.total_tVirCur[j,i] = median(coldAccRate.total_tvircur[j,w])
        coldMedian.total_tVirAcc[j,i] = median(coldAccRate.total_tviracc[j,w])
        hotMedian.total_tVirCur[j,i]  = median(hotAccRate.total_tvircur[j,w])
        hotMedian.total_tVirAcc[j,i]  = median(hotAccRate.total_tviracc[j,w])
        
        ; total masses:
        ; gal (hot+cold)
        coldTotal.gal_tVirCur[j,i]   = median(totalColdMass.gal_tvircur[j,w])
        coldTotal.gal_tVirAcc[j,i]   = median(totalColdMass.gal_tviracc[j,w])
        hotTotal.gal_tVirCur[j,i]    = median(totalHotMass.gal_tvircur[j,w])
        hotTotal.gal_tVirAcc[j,i]    = median(totalHotMass.gal_tviracc[j,w])
        
        ; gmem (hot+cold)
        coldTotal.gmem_tVirCur[j,i]  = median(totalColdMass.gmem_tvircur[j,w])
        coldTotal.gmem_tVirAcc[j,i]  = median(totalColdMass.gmem_tviracc[j,w])
        hotTotal.gmem_tVirCur[j,i]   = median(totalHotMass.gmem_tvircur[j,w])
        hotTotal.gmem_tVirAcc[j,i]   = median(totalHotMass.gmem_tviracc[j,w]) 
        
        ; stars (hot+cold)
        coldTotal.stars_tVirCur[j,i] = median(totalColdMass.stars_tvircur[j,w])
        coldTotal.stars_tVirAcc[j,i] = median(totalColdMass.stars_tviracc[j,w])
        hotTotal.stars_tVirCur[j,i]  = median(totalHotMass.stars_tvircur[j,w])
        hotTotal.stars_tVirAcc[j,i]  = median(totalHotMass.stars_tviracc[j,w]) 
        
        ; both=gal+stars (hot+cold)
        coldTotal.both_tVirCur[j,i]  = median(totalColdMass.both_tvircur[j,w])
        coldTotal.both_tVirAcc[j,i]  = median(totalColdMass.both_tviracc[j,w])
        hotTotal.both_tVirCur[j,i]   = median(totalHotMass.both_tvircur[j,w])
        hotTotal.both_tVirAcc[j,i]   = median(totalHotMass.both_tviracc[j,w])
        
        ; total=gal+stars+gmem (hot+cold)
        coldTotal.total_tVirCur[j,i] = median(totalColdMass.total_tvircur[j,w])
        coldTotal.total_tVirAcc[j,i] = median(totalColdMass.total_tviracc[j,w])
        hotTotal.total_tVirCur[j,i]  = median(totalHotMass.total_tvircur[j,w])
        hotTotal.total_tVirAcc[j,i]  = median(totalHotMass.total_tviracc[j,w])
      endfor
      
      for j=0,nCuts-1 do begin
        ; rates:
        coldMedian.gal_const[j,i]   = median(coldAccRate.gal_const[j,w])
        hotMedian.gal_const[j,i]    = median(hotAccRate.gal_const[j,w])
        coldMedian.gmem_const[j,i]  = median(coldAccRate.gmem_const[j,w])
        hotMedian.gmem_const[j,i]   = median(hotAccRate.gmem_const[j,w])
        coldMedian.stars_const[j,i] = median(coldAccRate.stars_const[j,w])
        hotMedian.stars_const[j,i]  = median(hotAccRate.stars_const[j,w])
        
        coldMedian.both_const[j,i]  = median(coldAccRate.both_const[j,w])
        hotMedian.both_const[j,i]   = median(hotAccRate.both_const[j,w])
        coldMedian.total_const[j,i] = median(coldAccRate.total_const[j,w])
        hotMedian.total_const[j,i]  = median(hotAccRate.total_const[j,w])
        
        ; total mass:
        coldTotal.gal_const[j,i]   = median(totalColdMass.gal_const[j,w])
        hotTotal.gal_const[j,i]    = median(totalHotMass.gal_const[j,w])
        coldTotal.gmem_const[j,i]  = median(totalColdMass.gmem_const[j,w])
        hotTotal.gmem_const[j,i]   = median(totalHotMass.gmem_const[j,w])
        coldTotal.stars_const[j,i] = median(totalColdMass.stars_const[j,w])
        hotTotal.stars_const[j,i]  = median(totalHotMass.stars_const[j,w])
        
        coldTotal.both_const[j,i]  = median(totalColdMass.both_const[j,w])
        hotTotal.both_const[j,i]   = median(totalHotMass.both_const[j,w])
        coldTotal.total_const[j,i] = median(totalColdMass.total_const[j,w])
        hotTotal.total_const[j,i]  = median(totalHotMass.total_const[j,w])
      endfor
      
      ; rate totals (same under any cold/hot definition)
      totalHCMedian.gal[i]   = median(coldAccRate.gal_const[0,w]+hotAccRate.gal_const[0,w])
      totalHCMedian.gmem[i]  = median(coldAccRate.gmem_const[0,w]+hotAccRate.gmem_const[0,w])
      totalHCMedian.stars[i] = median(coldAccRate.stars_const[0,w]+hotAccRate.stars_const[0,w])
      totalHCMedian.both[i]  = median(coldAccRate.both_const[0,w]+hotAccRate.both_const[0,w])
      totalHCMedian.tot[i]   = median(coldAccRate.total_const[0,w]+hotAccRate.total_const[0,w])
      
      ; total mass totals
      totalMassHC.gal[i]   = median(totalColdMass.gal_const[0,w]+totalHotMass.gal_const[0,w])
      totalMassHC.gmem[i]  = median(totalColdMass.gmem_const[0,w]+totalHotMass.gmem_const[0,w])
      totalMassHC.stars[i] = median(totalColdMass.stars_const[0,w]+totalHotMass.stars_const[0,w])
      totalMassHC.both[i]  = median(totalColdMass.both_const[0,w]+totalHotMass.both_const[0,w])
      totalMassHC.tot[i]   = median(totalColdMass.total_const[0,w]+totalHotMass.total_const[0,w])
    endif  
    
    ; --- cold fraction ---
    ; gal (Tmax)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac.gal_num ge minNumFrac,count)
    if count gt 0 then begin
      for j=0,nVirs-1 do fracMedian.gal_tVircur[j,i] = median(coldFrac.gal_tvircur[j,w])
      for j=0,nVirs-1 do fracMedian.gal_tViracc[j,i] = median(coldFrac.gal_tviracc[j,w])
      for j=0,nCuts-1 do fracMedian.gal_const[j,i]   = median(coldFrac.gal_const[j,w])
    endif
    
    ; gmem (Tmax)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac.gmem_num ge minNumFrac,count)
    if count gt 0 then begin
      for j=0,nVirs-1 do fracMedian.gmem_tVircur[j,i] = median(coldFrac.gmem_tvircur[j,w])
      for j=0,nVirs-1 do fracMedian.gmem_tViracc[j,i] = median(coldFrac.gmem_tviracc[j,w])
      for j=0,nCuts-1 do fracMedian.gmem_const[j,i]   = median(coldFrac.gmem_const[j,w])
    endif
    
    ; stars (Tmax)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac.stars_num ge minNumFrac,count)
    if count gt 0 then begin
      for j=0,nVirs-1 do fracMedian.stars_tVircur[j,i] = median(coldFrac.stars_tvircur[j,w])
      for j=0,nVirs-1 do fracMedian.stars_tViracc[j,i] = median(coldFrac.stars_tviracc[j,w])
      for j=0,nCuts-1 do fracMedian.stars_const[j,i]   = median(coldFrac.stars_const[j,w])
    endif
    
    ; both (Tmax)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac.gal_num+coldFrac.stars_num ge minNumFrac,count)
    if count gt 0 then begin
      for j=0,nVirs-1 do fracMedian.both_tVircur[j,i] = median(coldFrac.both_tvircur[j,w])
      for j=0,nVirs-1 do fracMedian.both_tViracc[j,i] = median(coldFrac.both_tviracc[j,w])
      for j=0,nCuts-1 do fracMedian.both_const[j,i]   = median(coldFrac.both_const[j,w])
    endif
    
    ; total (Tmax)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac.gal_num+coldFrac.gmem_num+coldFrac.stars_num ge minNumFrac,count)
    if count gt 0 then begin
      for j=0,nVirs-1 do fracMedian.total_tVircur[j,i] = median(coldFrac.total_tvircur[j,w])
      for j=0,nVirs-1 do fracMedian.total_tViracc[j,i] = median(coldFrac.total_tviracc[j,w])
      for j=0,nCuts-1 do fracMedian.total_const[j,i]   = median(coldFrac.total_const[j,w])
    endif
    
    ; gal (Tcur)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac_cur.gal_num ge minNumFrac,count)
    if count gt 0 then begin
      for j=0,nVirs-1 do fracMedian_cur.gal_tVircur[j,i] = median(coldFrac_cur.gal_tvircur[j,w])
      for j=0,nVirs-1 do fracMedian_cur.gal_tViracc[j,i] = median(coldFrac_cur.gal_tviracc[j,w])
      for j=0,nCuts-1 do fracMedian_cur.gal_const[j,i]   = median(coldFrac_cur.gal_const[j,w])
    endif
    
    ; gmem (Tcur)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac_cur.gmem_num ge minNumFrac,count)
    if count gt 0 then begin
      for j=0,nVirs-1 do fracMedian_cur.gmem_tVircur[j,i] = median(coldFrac_cur.gmem_tvircur[j,w])
      for j=0,nVirs-1 do fracMedian_cur.gmem_tViracc[j,i] = median(coldFrac_cur.gmem_tviracc[j,w])
      for j=0,nCuts-1 do fracMedian_cur.gmem_const[j,i]   = median(coldFrac_cur.gmem_const[j,w])
    endif
    
    ; total (Tcur)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac_cur.gal_num+coldFrac_cur.gmem_num ge minNumFrac,count)
    if count gt 0 then begin
      for j=0,nVirs-1 do fracMedian_cur.total_tVircur[j,i] = median(coldFrac_cur.total_tvircur[j,w])
      for j=0,nVirs-1 do fracMedian_cur.total_tViracc[j,i] = median(coldFrac_cur.total_tviracc[j,w])
      for j=0,nCuts-1 do fracMedian_cur.total_const[j,i]   = median(coldFrac_cur.total_const[j,w])
    endif     
    
  endfor
  
  r = {coldAccRate:coldAccRate,hotAccRate:hotAccRate,coldFrac:coldFrac,coldFrac_cur:coldFrac_cur,$
       fracMedian:fracMedian,fracMedian_cur:fracMedian_cur,$
       coldMedian:coldMedian,hotMedian:hotMedian,totalHCMedian:totalHCMedian,$
       totalColdMass:totalColdMass,totalHotMass:totalHotMass,coldTotal:coldTotal,hotTotal:hotTotal,totalMassHC:totalMassHC,$
       radIndGalAcc:sP.radIndGalAcc,radIndHaloAcc:sP.radIndHaloAcc,$
       logMassBins:logMassBins,logMassBinCen:logMassBinCen,TcutVals:sP.TcutVals,TvirVals:sP.TvirVals}

  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  return, r
  
end

; ------------------------------------------------------------------------------------------------------

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
  if 0 then begin
  start_PS,sP.plotPath + 'angmom.histos.all.'+accMode+'.'+sP.plotPrefix+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
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
  
  start_PS,sP.plotPath + 'angmom.histos.cold.'+accMode+'.'+sP.plotPrefix+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
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
  
  start_PS,sP.plotPath + 'angmom.histos.hot.'+accMode+'.'+sP.plotPrefix+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
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
  endif ;0
  
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
  if 0 then begin
  start_PS,sP.plotPath + 'accdt.histos.all.'+accMode+'.'+sP.plotPrefix+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
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
  
  start_PS,sP.plotPath + 'accdt.histos.cold.'+accMode+'.'+sP.plotPrefix+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
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
  
  start_PS,sP.plotPath + 'accdt.histos.hot.'+accMode+'.'+sP.plotPrefix+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
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
  endif ;0

  return,binnedVals
  
end
