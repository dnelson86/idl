; plotVsHaloMass.pro
; gas accretion project - bin quantities as a function of halo mass
; dnelson jun.2013

; haloMassBinValues(): bin accretion rate (in cold/hot) and cold fraction as a function of halo mass

function haloMassBinValues, sP=sP, accMode=accMode, timeWindow=TW

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  sgSelect = 'pri' ; only option for atS
  
  if ~sP.gfmWinds and accMode eq 'recycled' then message,'Error: Request recycled on non-winds run.'
  
  if sP.gfmWinds eq 1 and (accMode eq 'smooth' or accMode eq 'stripped' or accMode eq 'clumpy') then begin
    print,'Switching [' + accMode + '] to [' + accMode + '_rec] for GFM run.'
    accMode = accMode + '_rec'
  endif
  
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
  ;if file_test(saveFilename) then begin
  ;  restore,saveFilename
  ;  return,r
  ;endif   
  
  ; make a uniform gas selection at the start
  at = accretionTimes(sP=sP)
  mt = mergerTreeSubset(sP=sP)
  
  wAm = accModeInds(at=at,accMode=accMode,sP=sP)
    
  ; reverse histogram parent IDs of all particles/tracers in this selection
  gcIndOrig = mergerTreeRepParentIDs(mt=mt,sP=sP,/compactMtS)
  
  maxHist = max(gcIndOrig)
  
  hist_gal   = histogram(gcIndOrig[wAm.gal],min=0,max=maxHist,loc=loc_gal,rev=rev_gal)
  hist_gmem  = histogram(gcIndOrig[wAm.gmem],min=0,max=maxHist,loc=loc_gmem,rev=rev_gmem)
  hist_stars = histogram(gcIndOrig[wAm.stars],min=0,max=maxHist,loc=loc_stars,rev=rev_stars)
  hist_inter = histogram(gcIndOrig[wAm.inter],min=0,max=maxHist,loc=loc_inter,rev=rev_inter)
  if sP.gfMBHs ne 0 then $
    hist_bhs   = histogram(gcIndOrig[wAm.bhs],min=0,max=maxHist,loc=loc_bhs,rev=rev_bhs)
    
  gcIndOrig = !NULL
  
  ; load max temps, current tvir, tvir at accretion
  accTvir = gcSubsetProp(sP=sP,select=sgSelect,/accTvir,/accretionTimeSubset,accMode=accMode)
  curTvir = gcSubsetProp(sP=sP,select=sgSelect,/virTemp,/accretionTimeSubset,accMode=accMode)
  maxTemp = gcSubsetProp(sP=sP,select=sgSelect,/maxPastTemp,/accretionTimeSubset,accMode=accMode)
  
  ; load group cat for subgroup masses
  gc = loadGroupCat(sP=sP,/skipIDs)
  gcMasses = codeMassToLogMsun(gc.subgroupMass[mt.galcatIDList])
  gc = !NULL

  ; structure config
  constInds = [0] ; indices into template for constant
  tVirInds  = [1,2] ; indices into template for tvir
  
  ; which galaxyCat types contribute to the "allGal" and "total"?
  allgalInds = [0,2,4] ; gal,stars,bhs
  totalInds  = [0,1,2,3,4] ; gal,gmem,stars,bhs  
  
  typeLabels = ['gal','gmem','stars','inter','bhs','allgal','total']
  
  allgalInd = ( where( typeLabels eq 'allgal' ) )[0]
  totalInd  = ( where( typeLabels eq 'total' ) )[0]
  
  ; structures to store results (Tmax)
  template = { tConst  : fltarr(nCuts,n_elements(mt.galcatIDList))  ,$
               tVirCur : fltarr(nVirs,n_elements(mt.galcatIDList))  ,$
               tVirAcc : fltarr(nVirs,n_elements(mt.galcatIDList))  ,$
               num     : lonarr(n_elements(mt.galcatIDList))         }

  coldAccRate = {}

  foreach typeLabel,typeLabels do $  
    coldAccRate = mod_struct(coldAccRate, typeLabel, template)
  
  hotAccRate    = coldAccRate
  coldFrac      = coldAccRate
  totalColdMass = coldAccRate ; total mass
  totalHotMass  = coldAccRate
                   
  ; outflow rates, filled only for .gal (which searches in gmem and inter) and .gmem
  coldOutRate = { gal : template, gmem : template }
  hotOutRate  = coldOutRate
  coldNetRate = coldOutRate ; net (inflow-outflow) rates
  hotNetRate  = coldOutRate

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
      
      ; galaxy accretion defined as most recent of (rho,temp) joining time or 0.15rvir crossing time (headed in)
      loc_atime_gal = reform(at.accTimeRT[wAm.gal[loc_inds_gal]])
      
      r_crossing_time = reform(at.accTime[sP.radIndGalAcc,wAm.gal[loc_inds_gal]])
      w = where(r_crossing_time gt loc_atime_gal,count)
      if count gt 0 then loc_atime_gal[w] = r_crossing_time[w]
      
      loc_atime_gal = 1/loc_atime_gal - 1 ; redshift
      loc_atime_gal = redshiftToAgeFlat(loc_atime_gal)*1e9 ; yr
      
      ; make a count of those falling in the time window
      w = where(curtime - loc_atime_gal le loc_timeWindow,nloc)
      
      coldFrac.gal.num[i]    = nloc
      coldAccRate.gal.num[i] = nloc
      
      ; maximum past temps, cur and acc tvirs for only those particles in the time window
      if nloc gt 0 then begin
        loc_maxt_gal = maxTemp.gal[loc_inds_gal[w]]
        loc_curtvir_gal = curTvir.gal[loc_inds_gal[w]]
        loc_acctvir_gal = accTvir.gal[loc_inds_gal[w]]
  
        ; count mass elements with Tmax below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_maxt_gal le sP.TcutVals[j],count_cold,ncomp=count_hot)
          coldFrac.gal.tConst[j,i]    = float(count_cold) / nloc
          coldAccRate.gal.tConst[j,i] = count_cold
          hotAccRate.gal.tConst[j,i]  = count_hot
        endfor
        
        for j=0,nVirs-1 do begin
          ; count mass elements with Tmax below Tvir at current time
          w = where(10.0^loc_maxt_gal / 10.0^loc_curtvir_gal le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldFrac.gal.tvircur[j,i]    = float(count_cold) / nloc
          coldAccRate.gal.tvircur[j,i] = count_cold
          hotAccRate.gal.tvircur[j,i]  = count_hot
    
          ; count mass elements with Tmax below Tvir at accretion time
          w = where(10.0^loc_maxt_gal / 10.0^loc_acctvir_gal le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldFrac.gal.tviracc[j,i]    = float(count_cold) / nloc
          coldAccRate.gal.tviracc[j,i] = count_cold
          hotAccRate.gal.tviracc[j,i]  = count_hot
        endfor
      endif ;nloc>0
      
    endif ; hist_gal[i]

    if hist_gmem[i] gt 0 then begin
      ; list of indices of group member gas particles in this subgroup
      loc_inds_gmem = rev_gmem[rev_gmem[i]:rev_gmem[i+1]-1]
      
      ; corresponding accretion times for these particles
      ; -------------------------------------------------
      loc_atime_gmem = reform(at.accTime[sP.radIndHaloAcc,wAm.gmem[loc_inds_gmem]])
      loc_atime_gmem = 1/loc_atime_gmem - 1 ; redshift
      loc_atime_gmem = redshiftToAgeFlat(loc_atime_gmem)*1e9 ; yr
      
      ; make a count of those falling in the time window
      w = where(curtime - loc_atime_gmem le loc_timeWindow,nloc)
      
      coldFrac.gmem.num[i]    = nloc
      coldAccRate.gmem.num[i] = nloc
      
      ; maximum past temps, cur and acc tvirs for only those particles in the time window
      if nloc gt 0 then begin
        loc_maxt_gmem = maxTemp.gmem[loc_inds_gmem[w]]
        loc_curtvir_gmem = curTvir.gmem[loc_inds_gmem[w]]
        loc_acctvir_gmem = accTvir.gmem[loc_inds_gmem[w]]
  
        ; count mass elements with Tmax below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_maxt_gmem le sP.TcutVals[j],count_cold,ncomp=count_hot)
          coldFrac.gmem.tConst[j,i]    = float(count_cold) / nloc
          coldAccRate.gmem.tConst[j,i] = count_cold
          hotAccRate.gmem.tConst[j,i]  = count_hot
        endfor
        
        for j=0,nVirs-1 do begin
          ; count mass elements with Tmax below Tvir at current time
          w = where(10.0^loc_maxt_gmem / 10.0^loc_curtvir_gmem le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldFrac.gmem.tvircur[j,i]    = float(count_cold) / nloc
          coldAccRate.gmem.tvircur[j,i] = count_cold
          hotAccRate.gmem.tvircur[j,i]  = count_hot
    
          ; count mass elements with Tmax below Tvir at accretion time
          w = where(10.0^loc_maxt_gmem / 10.0^loc_acctvir_gmem le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldFrac.gmem.tviracc[j,i]    = float(count_cold) / nloc
          coldAccRate.gmem.tviracc[j,i] = count_cold
          hotAccRate.gmem.tviracc[j,i]  = count_hot
        endfor
      endif ; nloc>0
      
      ; corresponding --OUTFLOW-- times from the GALAXY for these particles
      ; defined as most recent of (rho,temp) leaving time or 0.15rvir crossing time (headed out)
      ; -------------------------------------------------------------------
      loc_otime = reform(at.outTimeRT[wAm.gmem[loc_inds_gmem]])
      
      r_crossing_time = reform(at.outTime[sP.radIndGalAcc,wAm.gmem[loc_inds_gmem]])
      w = where(r_crossing_time gt loc_otime,count)
      if count gt 0 then loc_otime[w] = r_crossing_time[w]
      
      loc_otime = 1/loc_otime - 1 ; redshift
      loc_otime = redshiftToAgeFlat(loc_otime)*1e9 ; yr
      
      ; make a count of those falling in the time window
      w = where(curtime - loc_otime le loc_timeWindow,nloc)
      
      coldOutRate.gal.num[i] += nloc
      
      ; maximum past temps, cur and acc tvirs for only those particles in the time window
      if nloc gt 0 then begin
        loc_maxt = maxTemp.gmem[loc_inds_gmem[w]]
        loc_curtvir = curTvir.gmem[loc_inds_gmem[w]]
        loc_acctvir = accTvir.gmem[loc_inds_gmem[w]]
  
        ; count mass elements with Tmax below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_maxt le sP.TcutVals[j],count_cold,ncomp=count_hot)
          coldOutRate.gal.tConst[j,i] += count_cold
          hotOutRate.gal.tConst[j,i]  += count_hot
        endfor
        
        for j=0,nVirs-1 do begin
          ; count mass elements with Tmax below Tvir at current time
          w = where(10.0^loc_maxt / 10.0^loc_curtvir le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldOutRate.gal.tvircur[j,i] += count_cold
          hotOutRate.gal.tvircur[j,i]  += count_hot
    
          ; count mass elements with Tmax below Tvir at accretion time
          w = where(10.0^loc_maxt / 10.0^loc_acctvir le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldOutRate.gal.tviracc[j,i] += count_cold
          hotOutRate.gal.tviracc[j,i]  += count_hot
        endfor
      endif ; nloc>0
      
      ; corresponding --OUTFLOW-- times from the HALO for these particles
      ; -------------------------------------------------------------------
      loc_otime = reform(at.outTime[sP.radIndHaloAcc,wAm.gmem[loc_inds_gmem]])
      loc_otime = 1/loc_otime - 1 ; redshift
      loc_otime = redshiftToAgeFlat(loc_otime)*1e9 ; yr
      
      ; make a count of those falling in the time window
      w = where(curtime - loc_otime le loc_timeWindow,nloc)
      
      coldOutRate.gmem.num[i] += nloc
      
      ; maximum past temps, cur and acc tvirs for only those particles in the time window
      if nloc gt 0 then begin
        loc_maxt = maxTemp.gmem[loc_inds_gmem[w]]
        loc_curtvir = curTvir.gmem[loc_inds_gmem[w]]
        loc_acctvir = accTvir.gmem[loc_inds_gmem[w]]
  
        ; count mass elements with Tmax below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_maxt le sP.TcutVals[j],count_cold,ncomp=count_hot)
          coldOutRate.gmem.tConst[j,i] += count_cold
          hotOutRate.gmem.tConst[j,i]  += count_hot
        endfor
        
        for j=0,nVirs-1 do begin
          ; count mass elements with Tmax below Tvir at current time
          w = where(10.0^loc_maxt / 10.0^loc_curtvir le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldOutRate.gmem.tvircur[j,i] += count_cold
          hotOutRate.gmem.tvircur[j,i]  += count_hot
    
          ; count mass elements with Tmax below Tvir at accretion time
          w = where(10.0^loc_maxt / 10.0^loc_acctvir le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldOutRate.gmem.tviracc[j,i] += count_cold
          hotOutRate.gmem.tviracc[j,i]  += count_hot
        endfor
      endif ; nloc>0
      
    endif ; hist_gmem[i]

    if hist_stars[i] gt 0 then begin
      ; list of indices of star gas particles in this subgroup
      loc_inds_stars = rev_stars[rev_stars[i]:rev_stars[i+1]-1]
      
      ; corresponding accretion times for these particles:
      ; stellar accretion defined as (rho,temp) joining time or 0.15rvir crossing time (most recent)
      loc_atime_stars = reform(at.accTimeRT[wAm.stars[loc_inds_stars]])
      
      r_crossing_time = reform(at.accTime[sP.radIndGalAcc,wAm.stars[loc_inds_stars]])
      w = where(r_crossing_time gt loc_atime_stars,count)
      if count gt 0 then loc_atime_stars[w] = r_crossing_time[w]
      
      ; convert from scale factor to age of the universe
      loc_atime_stars = 1/loc_atime_stars - 1 ; redshift
      loc_atime_stars = redshiftToAgeFlat(loc_atime_stars)*1e9 ; yr
      
      ; note: if no 0.15rvir crossing time exists (crossed as a star particle, not as a gas cell) then
      ;       set the accretion time to t=-1 which moves it outside any timeWindow (do not consider) (only SPH)
      w = where(r_crossing_time eq -1,count)
      if count gt 0 then loc_atime_stars[w] = -1.0
      
      ; make a count of those falling in the time window
      w = where(curtime - loc_atime_stars le loc_timeWindow,nloc)
      
      coldFrac.stars.num[i]    = nloc
      coldAccRate.stars.num[i] = nloc
      
      ; maximum past temps, cur and acc tvirs for only those particles in the time window
      if nloc gt 0 then begin
        loc_maxt_stars = maxTemp.stars[loc_inds_stars[w]]
        loc_curtvir_stars = curTvir.stars[loc_inds_stars[w]]
        loc_acctvir_stars = accTvir.stars[loc_inds_stars[w]]
  
        ; count mass elements with Tmax below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_maxt_stars le sP.TcutVals[j],count_cold,ncomp=count_hot)
          coldFrac.stars.tConst[j,i]    = float(count_cold) / nloc
          coldAccRate.stars.tConst[j,i] = count_cold
          hotAccRate.stars.tConst[j,i]  = count_hot
        endfor
        
        for j=0,nVirs-1 do begin
          ; count mass elements with Tmax below Tvir at current time
          w = where(10.0^loc_maxt_stars / 10.0^loc_curtvir_stars le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldFrac.stars.tvircur[j,i]    = float(count_cold) / nloc
          coldAccRate.stars.tvircur[j,i] = count_cold
          hotAccRate.stars.tvircur[j,i]  = count_hot
    
          ; count mass elements with Tmax below Tvir at accretion time
          w = where(10.0^loc_maxt_stars / 10.0^loc_acctvir_stars le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldFrac.stars.tviracc[j,i]    = float(count_cold) / nloc
          coldAccRate.stars.tviracc[j,i] = count_cold
          hotAccRate.stars.tviracc[j,i]  = count_hot
        endfor
      endif ;nloc>0
    endif ; hist_stars[i]
    
    if hist_inter[i] gt 0 then begin
    
      ; list of indices of inter gas particles in this subgroup
      loc_inds = rev_inter[rev_inter[i]:rev_inter[i+1]-1]
      
      ; galaxy --OUTFLOW-- defined as most recent of (rho,temp) leaving time or 0.15rvir crossing time (headed out)
      loc_otime = reform(at.outTimeRT[wAm.inter[loc_inds]])
      
      r_crossing_time = reform(at.outTime[sP.radIndGalAcc,wAm.inter[loc_inds]])
      w = where(r_crossing_time gt loc_otime,count)
      if count gt 0 then loc_otime[w] = r_crossing_time[w]
      
      loc_otime = 1/loc_otime - 1 ; redshift
      loc_otime = redshiftToAgeFlat(loc_otime)*1e9 ; yr
      
      ; make a count of those falling in the time window
      w = where(curtime - loc_otime le loc_timeWindow,nloc)
      
      coldOutRate.gal.num[i] += nloc
      
      ; maximum past temps, cur and acc tvirs for only those particles in the time window
      if nloc gt 0 then begin
        loc_maxt    = maxTemp.inter[loc_inds[w]]
        loc_curtvir = curTvir.inter[loc_inds[w]]
        loc_acctvir = accTvir.inter[loc_inds[w]]
  
        ; count mass elements with Tmax below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_maxt le sP.TcutVals[j],count_cold,ncomp=count_hot)
          coldOutRate.gal.tConst[j,i] += count_cold
          hotOutRate.gal.tConst[j,i]  += count_hot
        endfor
        
        for j=0,nVirs-1 do begin
          ; count mass elements with Tmax below Tvir at current time
          w = where(10.0^loc_maxt / 10.0^loc_curtvir le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldOutRate.gal.tvircur[j,i] += count_cold
          hotOutRate.gal.tvircur[j,i]  += count_hot
    
          ; count mass elements with Tmax below Tvir at accretion time
          w = where(10.0^loc_maxt / 10.0^loc_acctvir le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldOutRate.gal.tviracc[j,i] += count_cold
          hotOutRate.gal.tviracc[j,i]  += count_hot
        endfor
      endif ;nloc>0
    
    endif ; hist_inter[i]
    
    if sP.gfmBHs ne 0 then if hist_bhs[i] gt 0 then begin
    
      ; list of indices of galaxy gas particles in this subgroup
      loc_inds = rev_bhs[rev_bhs[i]:rev_bhs[i+1]-1]
      
      ; galaxy accretion defined as most recent of (rho,temp) joining time or 0.15rvir crossing time (headed in)
      loc_atime = reform(at.accTimeRT[wAm.bhs[loc_inds]])
      
      r_crossing_time = reform(at.accTime[sP.radIndGalAcc,wAm.bhs[loc_inds]])
      w = where(r_crossing_time gt loc_atime,count)
      if count gt 0 then loc_atime[w] = r_crossing_time[w]
      
      loc_atime = 1/loc_atime - 1 ; redshift
      loc_atime = redshiftToAgeFlat(loc_atime)*1e9 ; yr
      
      ; make a count of those falling in the time window
      w = where(curtime - loc_atime le loc_timeWindow,nloc)
      
      coldFrac.bhs.num[i]    = nloc
      coldAccRate.bhs.num[i] = nloc
      
      ; maximum past temps, cur and acc tvirs for only those particles in the time window
      if nloc gt 0 then begin
        loc_maxt = maxTemp.bhs[loc_inds[w]]
        loc_curtvir = curTvir.bhs[loc_inds[w]]
        loc_acctvir = accTvir.bhs[loc_inds[w]]
  
        ; count mass elements with Tmax below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_maxt le sP.TcutVals[j],count_cold,ncomp=count_hot)
          coldFrac.bhs.tConst[j,i]    = float(count_cold) / nloc
          coldAccRate.bhs.tConst[j,i] = count_cold
          hotAccRate.bhs.tConst[j,i]  = count_hot
        endfor
        
        for j=0,nVirs-1 do begin
          ; count mass elements with Tmax below Tvir at current time
          w = where(10.0^loc_maxt / 10.0^loc_curtvir le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldFrac.bhs.tvircur[j,i]    = float(count_cold) / nloc
          coldAccRate.bhs.tvircur[j,i] = count_cold
          hotAccRate.bhs.tvircur[j,i]  = count_hot
    
          ; count mass elements with Tmax below Tvir at accretion time
          w = where(10.0^loc_maxt / 10.0^loc_acctvir le sP.TvirVals[j],count_cold,ncomp=count_hot)
          coldFrac.bhs.tviracc[j,i]    = float(count_cold) / nloc
          coldAccRate.bhs.tviracc[j,i] = count_cold
          hotAccRate.bhs.tviracc[j,i]  = count_hot
        endfor
      endif ;nloc>0
    
    endif ; hist_bhs[i]
    
    ; right now, cold/hotAccRate holds just particle/tracer counts, convert to total masses and rates
    foreach k,[constInds,tVirInds] do begin ; 0=const, 1=tvircur, 2=tviracc (skip last one which is num)
      for j=0,n_tags(totalColdMass)-1 do begin
        ; convert accretion total mass (counts) to msun
        totalColdMass.(j).(k)[*,i] = coldAccRate.(j).(k)[*,i] * massPerPart * units.UnitMass_in_Msun
        totalHotMass.(j).(k)[*,i]  = coldAccRate.(j).(k)[*,i] * massPerPart * units.UnitMass_in_Msun
        
        ; convert accretion total(counts) to msun/year
        coldAccRate.(j).(k)[*,i] *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
        hotAccRate.(j).(k)[*,i]  *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
      endfor
      
      ; convert outflow total(counts) to msun/year
      coldOutRate.gal.(k)[*,i]  *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
      coldOutRate.gmem.(k)[*,i] *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
      hotOutRate.gal.(k)[*,i]   *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
      hotOutRate.gmem.(k)[*,i]  *= massPerPart * units.UnitMass_in_Msun / loc_timeWindow
    endforeach
  endfor ; i
  
  ; create totals of allgal=gal+stars+bhs, total=both+gmem+inter 
  for j=0,nVirs-1 do begin
    foreach k,tVirInds do begin
      coldFracNorm = 0
      
      ; loop over each galaxyCat type contributing to allGal
      foreach q,allgalInds do begin
        coldAccRate.(allgalInd).(k)[j,*] += coldAccRate.(q).(k)[j,*]
        hotAccRate.(allgalInd).(k)[j,*]  += hotAccRate.(q).(k)[j,*]
        
        totalColdMass.(allgalInd).(k)[j,*] += totalColdMass.(q).(k)[j,*]
        totalHotMass.(allgalInd).(k)[j,*]  += totalHotMass.(q).(k)[j,*]
        
        coldFrac.(allgalInd).(k)[j,*] += coldFrac.(q).(k)[j,*] * coldFrac.(q).num
        coldFracNorm += coldFrac.(q).num
      endforeach
      
      coldFrac.(allgalInd).(k)[j,*] /= coldFracNorm
      coldFrac.(allgalInd).num = coldFracNorm
      coldFracNorm = 0
      
      ; loop over each galaxyCat type contributing to total
      foreach q,totalInds do begin
        coldAccRate.(totalInd).(k)[j,*] += coldAccRate.(q).(k)[j,*]
        hotAccRate.(totalInd).(k)[j,*]  += hotAccRate.(q).(k)[j,*]
        
        totalColdMass.(totalInd).(k)[j,*] += totalColdMass.(q).(k)[j,*]
        totalHotMass.(totalInd).(k)[j,*]  += totalHotMass.(q).(k)[j,*]
        
        coldFrac.(totalInd).(k)[j,*] += coldFrac.(q).(k)[j,*] * coldFrac.(q).num
        coldFracNorm += coldFrac.(q).num
      endforeach
      
      coldFrac.(totalInd).(k)[j,*] /= coldFracNorm
      coldFrac.(totalInd).num = coldFracNorm
      
      ; calculate net rates
      ; TODO: add +stars+BHs terms
      coldNetRate.gal.(k)[j,*]  = coldAccRate.gal.(k)[j,*] - coldOutRate.gal.(k)[j,*]
      coldNetRate.gmem.(k)[j,*] = coldAccRate.gmem.(k)[j,*] - coldOutRate.gmem.(k)[j,*]
      hotNetRate.gal.(k)[j,*]   = hotAccRate.gal.(k)[j,*] - hotOutRate.gal.(k)[j,*]
      hotNetRate.gmem.(k)[j,*]  = hotAccRate.gmem.(k)[j,*] - hotOutRate.gmem.(k)[j,*]
      
    endforeach
  endfor ; nVirs
  
  for j=0,nCuts-1 do begin
    foreach k,constInds do begin
      coldFracNorm = 0
      
      ; loop over each galaxyCat type contributing to allGal
      foreach q,allgalInds do begin
        coldAccRate.(allgalInd).(k)[j,*] += coldAccRate.(q).(k)[j,*]
        hotAccRate.(allgalInd).(k)[j,*]  += hotAccRate.(q).(k)[j,*]
        
        totalColdMass.(allgalInd).(k)[j,*] += totalColdMass.(q).(k)[j,*]
        totalHotMass.(allgalInd).(k)[j,*]  += totalHotMass.(q).(k)[j,*]
        
        coldFrac.(allgalInd).(k)[j,*] += coldFrac.(q).(k)[j,*] * coldFrac.(q).num
        coldFracNorm += coldFrac.(q).num
      endforeach
      
      coldFrac.(allgalInd).(k)[j,*] /= coldFracNorm
      coldFrac.(allgalInd).num = coldFracNorm ; redundant for clarity
      coldFracNorm = 0
      
      ; loop over each galaxyCat type contributing to total
      foreach q,totalInds do begin
        coldAccRate.(totalInd).(k)[j,*] += coldAccRate.(q).(k)[j,*]
        hotAccRate.(totalInd).(k)[j,*]  += hotAccRate.(q).(k)[j,*]
        
        totalColdMass.(totalInd).(k)[j,*] += totalColdMass.(q).(k)[j,*]
        totalHotMass.(totalInd).(k)[j,*]  += totalHotMass.(q).(k)[j,*]
        
        coldFrac.(totalInd).(k)[j,*] += coldFrac.(q).(k)[j,*] * coldFrac.(q).num
        coldFracNorm += coldFrac.(q).num
      endforeach
      
      coldFrac.(totalInd).(k)[j,*] /= coldFracNorm
      coldFrac.(totalInd).num = coldFracNorm ; redundant for clarity
      
      ; calculate net rates
      ; TODO: add +stars+BHs terms
      coldNetRate.gal.(k)[j,*]  = coldAccRate.gal.(k)[j,*] - coldOutRate.gal.(k)[j,*]
      coldNetRate.gmem.(k)[j,*] = coldAccRate.gmem.(k)[j,*] - coldOutRate.gmem.(k)[j,*]
      hotNetRate.gal.(k)[j,*]   = hotAccRate.gal.(k)[j,*] - hotOutRate.gal.(k)[j,*]
      hotNetRate.gmem.(k)[j,*]  = hotAccRate.gmem.(k)[j,*] - hotOutRate.gmem.(k)[j,*]
      
    endforeach
  endfor ; nCuts
  
  ; bin fractions into halo mass bins and make median lines
  logMassBins   = [findgen(21)/10 + 9.0,11.1,11.25,11.5,11.75,11.9,13.1]
  logMassNBins  = n_elements(logMassBins)-1
  logMassBinCen = 0.5 * (logMassBins + shift(logMassBins,-1))
  logMassBinCen = logMassBinCen[0:-2]
  
  ; structures to store the binned values
  template = { tConst  : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
               tVirCur : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
               tVirAcc : fltarr(nVirs,logMassNbins) + !values.f_nan  }
  
  coldMedian = {}
  
  foreach typeLabel,typeLabels do $  
    coldMedian = mod_struct(coldMedian, typeLabel, template)
  
  hotMedian   = coldMedian
  fracMedian  = coldMedian
  coldTotal   = coldMedian ; total mass (not rates)
  hotTotal    = coldMedian ; total mass (not rates)
  
  ; outflow and net rates
  coldOutMedian = { gal : template, gmem : template }
  hotOutMedian  = coldOutMedian
  coldNetMedian = coldOutMedian
  hotNetMedian  = coldOutMedian
  
  ; total of hot+cold (separated by gal, gmem, stars, both=gal+stars, tot=gal+stars+gmem)
  template = fltarr(logMassNbins) + !values.f_nan
  
  totalHCMedian = {}
                    
  foreach typeLabel,typeLabels do $  
    totalHCMedian = mod_struct(totalHCMedian, typeLabel, template)
                    
  totalMassHC = totalHCMedian
                 
  ; calculate median values in bins of halo mass
  for i=0,logMassNbins-1 do begin

    ; --- accretion rate ---
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1],count)
    
    if count gt 0 then begin
    
      for j=0,nVirs-1 do begin
        ; loop over each tVir type
        foreach q,tVirInds do begin
          ; loop over each galaxyCat type (gal,gmem,stars,...)
          for k=0,n_tags(coldMedian)-1 do begin
            ; median accretion rates
            coldMedian.(k).(q)[j,i] = median(coldAccRate.(k).(q)[j,w])
            hotMedian.(k).(q)[j,i]  = median(hotAccRate.(k).(q)[j,w])
            
            ; median total masses
            coldTotal.(k).(q)[j,i] = median(totalColdMass.(k).(q)[j,w])
            hotTotal.(k).(q)[j,i]  = median(totalHotMass.(k).(q)[j,w])
          endfor
          
          ; outflow and net (gal,gmem)
          for k=0,n_tags(coldOutMedian)-1 do begin
            coldOutMedian.(k).(q)[j,i] = median(coldOutRate.(k).(q)[j,w])
            hotOutMedian.(k).(q)[j,i]  = median(hotOutRate.(k).(q)[j,w])
            coldNetMedian.(k).(q)[j,i] = median(coldNetRate.(k).(q)[j,w])
            hotNetMedian.(k).(q)[j,i]  = median(hotNetRate.(k).(q)[j,w])
          endfor
        endforeach
      endfor
      
      for j=0,nCuts-1 do begin
        ; loop over each tVir type
        foreach q,constInds do begin
          ; loop over each galaxyCat type (gal,gmem,stars,...)
          for k=0,n_tags(coldMedian)-1 do begin
            ; median accretion rates
            coldMedian.(k).(q)[j,i] = median(coldAccRate.(k).(q)[j,w])
            hotMedian.(k).(q)[j,i]  = median(hotAccRate.(k).(q)[j,w])
            
            ; median total masses
            coldTotal.(k).(q)[j,i] = median(totalColdMass.(k).(q)[j,w])
            hotTotal.(k).(q)[j,i]  = median(totalHotMass.(k).(q)[j,w])
          endfor
          
          ; outflow and net (gal,gmem)
          for k=0,n_tags(coldOutMedian)-1 do begin
            coldOutMedian.(k).(q)[j,i] = median(coldOutRate.(k).(q)[j,w])
            hotOutMedian.(k).(q)[j,i]  = median(hotOutRate.(k).(q)[j,w])
            coldNetMedian.(k).(q)[j,i] = median(coldNetRate.(k).(q)[j,w])
            hotNetMedian.(k).(q)[j,i]  = median(hotNetRate.(k).(q)[j,w])
          endfor
        endforeach
      endfor
      
      ; rate totals (same under any cold/hot definition)
      for j=0,n_tags(totalHCMedian)-1 do begin
        totalHCMedian.(j)[i] = median( coldAccRate.(j).(0)[0,w] + hotAccRate.(j).(0)[0,w] )
        totalMassHC.(j)[i]   = median( totalColdMass.(j).(0)[0,w] + totalHotMass.(j).(0)[0,w] )
      endfor
      
    endif 
    
    ; --- cold fraction ---
    for k=0,n_tags(coldFrac)-1 do begin

      w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and coldFrac.(k).num ge minNumFrac,count)
      
      if count gt 0 then begin
        for j=0,nVirs-1 do fracMedian.(k).tVircur[j,i] = median(coldFrac.(k).tVirCur[j,w])
        for j=0,nVirs-1 do fracMedian.(k).tViracc[j,i] = median(coldFrac.(k).tVirAcc[j,w])
        for j=0,nCuts-1 do fracMedian.(k).tConst[j,i]  = median(coldFrac.(k).tConst[j,w])
      endif
    endfor
    
  endfor ; logMassNBins

  r = {coldAccRate:coldAccRate,hotAccRate:hotAccRate,coldFrac:coldFrac,$
       coldOutRate:coldOutRate,hotOutRate:hotOutRate,coldNetRate:coldNetRate,hotNetRate:hotNetRate,$
       fracMedian:fracMedian,coldMedian:coldMedian,hotMedian:hotMedian,totalHCMedian:totalHCMedian,$
       totalColdMass:totalColdMass,totalHotMass:totalHotMass,coldTotal:coldTotal,hotTotal:hotTotal,totalMassHC:totalMassHC,$
       coldOutMedian:coldOutMedian,hotOutMedian:hotOutMedian,coldNetMedian:coldNetMedian,hotNetMedian:hotNetMedian,$
       allgalInd:allgalInd,totalInd:totalInd,typeLabels:typeLabels,allgalInds:allgalInds,totalInds:totalInds,$
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
  accTvirAtS = gcSubsetProp(sP=sP,select=sgSelect,/accTvir,/accretionTimeSubset,accMode=accMode)
  maxTempAtS = gcSubsetProp(sP=sP,select=sgSelect,/maxPastTemp,/accretionTimeSubset,accMode=accMode)

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
  parMasses = gcSubsetProp(sP=sP,select=sgSelect,/parMass,/accretionTimeSubset)

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
  accTvirAtS = gcSubsetProp(sP=sP,select=sgSelect,/accTvir,/accretionTimeSubset,accMode=accMode)
  maxTempAtS = gcSubsetProp(sP=sP,select=sgSelect,/maxPastTemp,/accretionTimeSubset,accMode=accMode)

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
  parMasses = gcSubsetProp(sP=sP,select=sgSelect,/parMass,/accretionTimeSubset)

  ; convert (k=0) r=rvir accretion times to nearest snap
  snapTimes = snapNumToRedshift(snap=0,/time,/all,sP=sP)
  snapTimes = snapTimes[where(snapTimes ge 0)] ; remove -1
    
  accSnap = {gal  : intarr(n_elements(at.accTime_gal[sP.atIndMode,*]))   ,$
             gmem : intarr(n_elements(at.accTime_gmem[sP.atIndMode,*]))   }
             
  accSnap.gal[wAm.gal]   = value_locate(snapTimes,at.accTime_gal[sP.atIndMode,wAm.gal]) 
  accSnap.gmem[wAm.gmem] = value_locate(snapTimes,at.accTime_gmem[sP.atIndMode,wAm.gmem])

  if n_elements(gcIndOrig.gal) ne n_elements(parMasses.gal) or $
     n_elements(gcIndOrig.gal) ne n_elements(wAm.galMask) then message,'Error: Data error.'
  w = where(accSnap.gal[wAm.gal] le 0 or accSnap.gal[wAm.gal] gt sP.snapRange[1],count1)
  w = where(accSnap.gmem[wAm.gmem] le 0 or accSnap.gmem[wAm.gmem] gt sP.snapRange[1],count2)
  if count1 gt 0 or count2 gt 0 then message,'Error: Bad snapshot time mapping.'

  ; each accDt now stores the scale factor when that particle crossed the radius of interest
  ; convert to delta time (Myr) since crossing rvir
  for k=1,n_elements(at.rVirFacs)-1 do begin
    gal_w  = where(wAm.galMask eq 1B and at.accTime_gal[sP.atIndMode,*] ne -1 and at.accTime_gal[k,*] ne -1,$
                   count_gal,comp=gal_wc,ncomp=countc_gal)
    gmem_w = where(wAm.gmemMask eq 1B and at.accTime_gmem[sP.atIndMode,*] ne -1 and at.accTime_gmem[k,*] ne -1,$
                   count_gmem,comp=gmem_wc,ncomp=countc_gmem)

    if count_gal gt 0 then begin
      ; for r=rvir, convert scale factor to redshift to age of universe in Gyr
      galAgeRvir  = redshiftToAgeFlat(1/reform(at.accTime_gal[sP.atIndMode,gal_w])-1)
      
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
      gmemAgeRvir = redshiftToAgeFlat(1/reform(at.accTime_gmem[sP.atIndMode,gmem_w])-1)
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
                    at.accTime_gal[sP.atIndMode,*] ne -1 and at.accTime_gal[k,*] ne -1 and $
                    maxTemp.gal gt accTvir.gal and wAm.galMask eq 1B,count_gal)   
      if count_gal gt 0 then begin
        binnedVals.hotMode.mean_gal[k,i]   = mean(at.accTime_gal[k,w_gal])
        binnedVals.hotMode.median_gal[k,i] = median(at.accTime_gal[k,w_gal])
        binnedVals.hotMode.stddev_gal[k,i] = stddev(at.accTime_gal[k,w_gal])
      endif
      
      ; hot gmem
      w_gmem = where(parMasses.gmem gt logMassBins[i] and parMasses.gmem le logMassBins[i+1] and $
                     at.accTime_gmem[sP.atIndMode,*] ne -1 and at.accTime_gmem[k,*] ne -1 and $
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
                    at.accTime_gal[sP.atIndMode,*] ne -1 and at.accTime_gal[k,*] ne -1 and $
                    maxTemp.gal le accTvir.gal and wAm.galMask eq 1B,count_gal)   
      if count_gal gt 0 then begin
        binnedVals.coldMode.mean_gal[k,i]   = mean(at.accTime_gal[k,w_gal])
        binnedVals.coldMode.median_gal[k,i] = median(at.accTime_gal[k,w_gal])
        binnedVals.coldMode.stddev_gal[k,i] = stddev(at.accTime_gal[k,w_gal])
      endif
      
      ; cold gmem
      w_gmem = where(parMasses.gmem gt logMassBins[i] and parMasses.gmem le logMassBins[i+1] and $
                     at.accTime_gmem[sP.atIndMode,*] ne -1 and at.accTime_gmem[k,*] ne -1 and $
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
      w = where(at.accTime_gal[sP.atIndMode,*] ne -1 and at.accTime_gal[k,*] ne -1 and wAm.galMask eq 1B,count)
      if count gt 1 then begin
        hist = histogram(at.accTime_gal[k,w],binsize=binsize,loc=loc)
        cgplot,loc+binsize*0.5,hist/float(total(hist)),line=0,color=getColor(k),/overplot
        cgplot,[mean(at.accTime_gal[k,w]),mean(at.accTime_gal[k,w])],[0.8,1.0],line=0,color=getColor(k),/overplot
        cgplot,[median(at.accTime_gal[k,w]),median(at.accTime_gal[k,w])],[0.5,0.7],line=0,color=getColor(k),/overplot
      endif
      w = where(at.accTime_gmem[sP.atIndMode,*] ne -1 and at.accTime_gmem[k,*] ne -1 and wAm.gmemMask eq 1B,count)
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
      w = where(at.accTime_gal[sP.atIndMode,*] ne -1 and at.accTime_gal[k,*] ne -1 and wAm.galMask eq 1B and $
                maxTemp.gal le accTvir.gal,count)
      if count gt 1 then begin
        hist = histogram(at.accTime_gal[k,w],binsize=binsize,loc=loc)
        cgplot,loc+binsize*0.5,hist/float(total(hist)),line=0,color=getColor(k),/overplot
        cgplot,[mean(at.accTime_gal[k,w]),mean(at.accTime_gal[k,w])],[0.8,1.0],line=0,color=getColor(k),/overplot
        cgplot,[median(at.accTime_gal[k,w]),median(at.accTime_gal[k,w])],[0.5,0.7],line=0,color=getColor(k),/overplot
      endif
      w = where(at.accTime_gmem[sP.atIndMode,*] ne -1 and at.accTime_gmem[k,*] ne -1 and wAm.gmemMask eq 1B,count)
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
      w = where(at.accTime_gal[sP.atIndMode,*] ne -1 and at.accTime_gal[k,*] ne -1 and wAm.galMask eq 1B and $
                maxTemp.gal gt accTvir.gal,count)
      if count gt 1 then begin
        hist = histogram(at.accTime_gal[k,w],binsize=binsize,loc=loc)
        cgplot,loc+binsize*0.5,hist/float(total(hist)),line=0,color=getColor(k),/overplot
        cgplot,[mean(at.accTime_gal[k,w]),mean(at.accTime_gal[k,w])],[0.8,1.0],line=0,color=getColor(k),/overplot
        cgplot,[median(at.accTime_gal[k,w]),median(at.accTime_gal[k,w])],[0.5,0.7],line=0,color=getColor(k),/overplot
      endif
      w = where(at.accTime_gmem[sP.atIndMode,*] ne -1 and at.accTime_gmem[k,*] ne -1 and wAm.gmemMask eq 1B,count)
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
