; binVsHaloMass.pro
; gas accretion project - bin quantities as a function of halo mass
; dnelson dec.2013

; haloMassBinValues(): bin accretion rate (in cold/hot) and cold fraction as a function of halo mass
; approach controlled by sP.newAccRates value:
;   model (0):
;    accretion time defined as: switching side of (rho,temp) cut or 0.15rvir crossing (headed in)
;    outflow time defined as: switching side of (rho,temp) cut or 0.15rvir crossing (headed out)
;   model (1): 
;    accretion time defined as: first 0.15rvir crossing
;    outflow time defined as: no outflow
;   model (3): (see accretionFlags.pro)
;    accretion time defined as: membership mismatch between origSnap and targetSnap
;    outflow time defined as: membership mismatch between targetSnap and origSnap
      
function haloMassBinValues, sP=sP, accMode=accMode, timeWindow=TW

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
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
  if sP.newAccRates gt 0 then arTag = '_new' + str(sP.newAccRates) else arTag = ''
  saveFilename = sP.derivPath + 'binnedVals/binVals.' + sP.saveTag + '.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.cut' + str(nCuts) + '.vir' + str(nVirs) + '.' + accMode + $
    '.r' + str(sP.radIndHaloAcc) + '.r' + str(sP.radIndGalAcc) + '_tw' + str(TW) + arTag + '.sav'

  if file_test(saveFilename) then begin
    restore, saveFilename
    return, r
  endif
  
  ; make a uniform gas selection at the start
  mt  = mergerTreeSubset(sP=sP)
  wAm = accModeInds(mt=mt,accMode=accMode,sP=sP)
  galcat = galaxyCat(sP=sP)
  
  ; reverse histogram parent IDs of all particles/tracers in this selection
  if sP.trMCPerCell gt 0 then gcIndOrig = mt.gcIndOrigTrMC
  if sP.trMCPerCell eq 0 then gcIndOrig = mt.gcIndOrig
  if sP.trMCPerCell lt 0 then gcIndOrig = mt.gcIndOrigTrVel
  
  maxHist = max(gcIndOrig)
  
  ; load group cat for subgroup masses
  gc = loadGroupCat(sP=sP,/skipIDs)
  gcMasses = codeMassToLogMsun(gc.subgroupMass)
  
  ; for zoom runs, only process the single target halo (should be the most massive)
  if sP.zoomLevel ne 0 then begin
    hInd = zoomTargetHalo(sP=sP,gc=gc)
    if hInd ne 0 then message,'Handle.'
    maxHist = 0
  endif
  
  gc = !NULL

  for i=0,n_tags(wAm)-1 do begin
    hist_type = histogram(gcIndOrig[wAm.(i)],min=0,max=maxHist,loc=loc,rev=rev_type)
    hist = mod_struct( hist, (tag_names(wAm))[i], hist_type )
    rev  = mod_struct( rev,  (tag_names(wAm))[i], rev_type )
  endfor

  hist_type = !NULL
  gcIndOrig = !NULL
  
  ; load max temps, current tvir, tvir at accretion
  accTvir = gcSubsetProp(sP=sP,/accTvir,/accretionTimeSubset,accMode=accMode)
  curTvir = gcSubsetProp(sP=sP,/virTemp,/accretionTimeSubset,accMode=accMode)
  maxTemp = gcSubsetProp(sP=sP,/maxPastTemp,/accretionTimeSubset,accMode=accMode)

  ; structure config
  constInds = [0] ; indices into template for constant
  tVirInds  = [1,2] ; indices into template for tvir
  
  ; which galaxyCat types contribute to "total"?
  typeLabels = [ tag_names(galcat.types), 'total' ]
  
  totalInds = where( typeLabels ne 'total' ) ; where to take total from
  totalInd  = ( where( typeLabels eq 'total' ) )[0] ; where to store total

  galcat = !NULL
  
  ; make templates: structures to store results (Tmax)
  tempHC   = { tConst  : fltarr(nCuts,maxHist+1)  ,$
               tVirCur : fltarr(nVirs,maxHist+1)  ,$
               tVirAcc : fltarr(nVirs,maxHist+1)  ,$
               num     : lonarr(maxHist+1)         }

  template = { hot : tempHC, cold : tempHC }

  foreach typeLabel,typeLabels do $  
    templateTypes = mod_struct(templateTypes, typeLabel, template)
    
  foreach typeLabel,typeLabels do $
    templateCF = mod_struct(templateCF, typeLabel, tempHC)
  
  ; for (1) galaxy and (2) halo, duplicate template for: accretion (inflow), outflow, net, and coldfrac
  galaxy = { accRate : templateTypes, outRate : templateTypes, $
             netRate : templateTypes, coldFrac : templateCF }

  halo   = { accRate : templateTypes, outRate : templateTypes, $
             netRate : templateTypes, coldFrac : templateCF }

  if sP.newAccRates gt 0 then begin
    if sP.newAccRates eq 1 then begin
      ; (1) accretionRates just based on first 0.15rvir crossing
      print,'NEW ACCRETION RATES (1):'
      at = accretionTimes(sP=sP)
      
      sz = size(at.accTime)
      if sz[0] ne 2 or sz[1] ne n_elements(at.rVirFacs)+2 then message,'Error'
      
    endif else begin
      ; (3) accretionRates based on membership at previous snapshot
      print,'NEW ACCRETION RATES (3):'
      targetSnap = closest( snapTimes, curtime-timeWindow )
      twError = (curtime-snapTimes[targetSnap]) / timeWindow - 1.0
    
      print,'Selected targetSnap = '+str(targetSnap)+' with timeWindow error ['+$
        string(twError*100,format='(f4.1)')+'%]'

      af = accretionFlags(sP=sP,targetSnap=targetSnap)
      
      ; now we also need quite a few values for the galcatTarget members
      sP.snap = targetSnap
      
      mtTarget  = mergerTreeSubset(sP=sP)
      wAmTarget = accModeInds(mt=mtTarget,accMode=accMode,sP=sP)
  
      ; reverse histogram parent IDs (at origSnap!) of all particles/tracers in this targetSnap selection
      for i=0,n_tags(wAmTarget)-1 do begin
        hist_type = histogram(af.origParIDsTarg[wAmTarget.(i)],min=0,max=maxHist,loc=loc,rev=rev_type)
        histTarget = mod_struct( histTarget, (tag_names(wAmTarget))[i], hist_type )
        revTarget  = mod_struct( revTarget,  (tag_names(wAmTarget))[i], rev_type )
      endfor
      
      mtTarget  = !NULL
      rev_type  = !NULL
      hist_type = !NULL
      gcIndOrig = !NULL
   
      ; load max temps, current tvir, tvir at accretion
      accTvirTarget = gcSubsetProp(sP=sP,/accTvir,/accretionTimeSubset,accMode=accMode)
      curTvirTarget = gcSubsetProp(sP=sP,/virTemp,/accretionTimeSubset,accMode=accMode)
      maxTempTarget = gcSubsetProp(sP=sP,/maxPastTemp,/accretionTimeSubset,accMode=accMode)
      
    endelse
  endif else begin
    ; (0) accretionRates based on accretionTimes (crossing of certain radii in certain directions)
    at = accretionTimes(sP=sP)
    maxAccTimeRT = max(at.accTimeRT)
    maxOutTimeRT = max(at.outTimeRT)
  endelse  
        
  ; loop over all tracked subgroups
  for i=0L,maxHist do begin

    if i mod (maxHist/10) eq 0 then print,string(float(i)/(maxHist+1)*100,format='(f5.2)')+'%'
  
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
    
    ; loop over types
    for k=0,n_tags(wAm)-1 do begin
      if (tag_names(wAm))[k] ne (tag_names(rev))[k] then message,'Check'
      
      ; local indices of this type for this halo
      loc_rev  = rev.(k)
	  
	if loc_rev[i] eq loc_rev[i+1] then continue ; no particles of type in this group, or whole group empty (secondary)
	  
      loc_inds = loc_rev[ loc_rev[i] : loc_rev[i+1]-1 ]
      
      ; targetSnap: local indices of this type for this halo (the tracked progenitor of this halo)
      if sP.newAccRates eq 3 then begin
        loc_rev_target = revTarget.(k)
        targetSkipFlag = 1
        
        if loc_rev_target[i] ne loc_rev_target[i+1] then begin
          targetSkipFlag = 0 ; process below
          loc_inds_target = loc_rev_target[ loc_rev_target[i] : loc_rev_target[i+1]-1 ]
        endif
      endif
      
      ; --- GALAXY ACCRETION ---

      if sP.newAccRates gt 0 then begin
        if sP.newAccRates eq 1 then begin
          loc_atime = reform( at.accTime[ sP.accRateInd, (wAm.(k))[ loc_inds ] ] )
          
          loc_atime = 1/loc_atime - 1 ; redshift
          loc_atime = redshiftToAgeFlat(loc_atime) * 1e9 ; yr
        
          w = where(curtime - loc_atime le loc_timeWindow, nloc)
        endif else begin
          loc_accFlag = af.accFlag[ (wAm.(k))[loc_inds] ]
          
          w = where(loc_accFlag eq 1, nloc) ; 1 = galaxy tracer accretion
        endelse
      endif else begin
        ; decide between using outTimeRT and accTimeRT (satisfying/failing the rho,temp cut, respectively)
        ; based on status at sP.snap (e.g. gmem by definition fails rho,temp cut)
        loc_atime = reform( at.accTimeRT[ (wAm.(k))[loc_inds] ] )
      
        ; those that fail at sP.snap (e.g. gmem), use satisfying the cut as the accretion time
        w = where( at.accTimeRT[ (wAm.(k))[loc_inds] ] eq maxAccTimeRT, count)
        if count gt 0 then loc_atime[w] = at.outTimeRT[ (wAm.(k))[loc_inds[w]] ]
      
        ; for those that neither satisfy nor fail the cut at sP.snap (all stars, stars in inter) we
        ; cannot use the (rho,temp) cut, since don't know where they will fall when switching back
        ; to a gas parent, so just use the radial cut
        w_noCut = where( at.outTimeRT[ (wAm.(k))[loc_inds] ] ne maxOutTimeRT and $
                         at.accTimeRT[ (wAm.(k))[loc_inds] ] ne maxAccTimeRT, count_noCut)
        if count_noCut gt 0 then loc_atime[w_noCut] = -1
		
        ; modify based on 0.15rvir crossing time (the latest / lowest redshift of the two criteria)
        r_crossing_time = reform( at.accTime[ sP.radIndGalAcc, (wAm.(k))[ loc_inds ] ] )
        w = where(r_crossing_time gt loc_atime, count)
        if count gt 0 then loc_atime[w] = r_crossing_time[w]
	  
        loc_atime = 1/loc_atime - 1 ; redshift
        loc_atime = redshiftToAgeFlat(loc_atime) * 1e9 ; yr
      
        w = where(curtime - loc_atime le loc_timeWindow, nloc)
      endelse
      
      galaxy.coldFrac.(k).num[i] = nloc
      galaxy.coldFrac.total.num[i] += nloc
      
      if nloc gt 0 then begin
        ; maximum past temps, cur and acc tvirs for only those particles in the time window
        loc_maxTemp = maxTemp.(k)[ loc_inds[w] ]
        loc_curTvir = curTvir.(k)[ loc_inds[w] ]
        loc_accTvir = accTvir.(k)[ loc_inds[w] ]
        
        ; count mass elements with Tmax below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_maxTemp le sP.TcutVals[j],count_cold,ncomp=count_hot)
          galaxy.coldFrac.(k).tConst[j,i]     = float(count_cold) / nloc
          galaxy.coldFrac.total.tConst[j,i]   += count_cold
          galaxy.accRate.(k).cold.tConst[j,i] = count_cold
          galaxy.accRate.(k).hot.tConst[j,i]  = count_hot
        endfor
        
        for j=0,nVirs-1 do begin
          ; count mass elements with Tmax below Tvir at current time
          w = where(10.0^loc_maxTemp / 10.0^loc_curTvir le sP.TvirVals[j],count_cold,ncomp=count_hot)
          galaxy.coldFrac.(k).tVirCur[j,i]     = float(count_cold) / nloc
          galaxy.coldFrac.total.tVirCur[j,i]   += count_cold
          galaxy.accRate.(k).cold.tVirCur[j,i] = count_cold
          galaxy.accRate.(k).hot.tVirCur[j,i]  = count_hot
    
          ; count mass elements with Tmax below Tvir at accretion time
          w = where(10.0^loc_maxTemp / 10.0^loc_accTvir le sP.TvirVals[j],count_cold,ncomp=count_hot)
          galaxy.coldFrac.(k).tVirAcc[j,i]     = float(count_cold) / nloc
          galaxy.coldFrac.total.tVirAcc[j,i]   += count_cold
          galaxy.accRate.(k).cold.tVirAcc[j,i] = count_cold
          galaxy.accRate.(k).hot.tVirAcc[j,i]  = count_hot
        endfor
      endif ;nloc>0
      
      ; --- GALAXY OUTFLOW ---
      
      if sP.newAccRates gt 0 then begin
        if sP.newAccRates eq 1 then begin
          ; (1) no outflow under this model
          nloc = 0
        endif else begin
          ; (3) use the values at targetSnap for outflow calculation
          nloc = 0
          
          if targetSkipFlag eq 0 then begin
            loc_accFlagTarg = af.accFlagTarg[ (wAmTarget.(k))[loc_inds_target] ]
          
            w = where(loc_accFlagTarg eq 1, nloc) ; 1 = galaxy tracer outflow from targetSnap
          
            if nloc gt 0 then begin
              loc_maxTemp = maxTempTarget.(k)[ loc_inds_target[w] ]
              loc_curTvir = curTvirTarget.(k)[ loc_inds_target[w] ]
              loc_accTvir = accTvirTarget.(k)[ loc_inds_target[w] ]
            endif
          endif
          
        endelse ; (nAR: 1/3)
      endif else begin
        loc_atime = reform( at.outTimeRT[ (wAm.(k))[loc_inds] ] )
      
        ; those that satisfy at sP.snap (e.g. gal gas), use failing the cut as the outflow time
        w = where( at.outTimeRT[ (wAm.(k))[loc_inds] ] eq maxOutTimeRT, count)
        if count gt 0 then loc_atime[w] = at.accTimeRT[ (wAm.(k))[loc_inds[w]] ]
      
        if count_noCut gt 0 then loc_atime[w_noCut] = -1 ; see above
      
        ; modify by 0.15rvir crossing time
        r_crossing_time = reform( at.outTime[ sP.radIndGalAcc, (wAm.(k))[ loc_inds ] ] )
        w = where(r_crossing_time gt loc_atime, count)
        if count gt 0 then loc_atime[w] = r_crossing_time[w]
      
        loc_atime = 1/loc_atime - 1 ; redshift
        loc_atime = redshiftToAgeFlat(loc_atime) * 1e9 ; yr
      
        w = where(curtime - loc_atime le loc_timeWindow, nloc)
        
        if nloc gt 0 then begin
          ; maximum past temps, cur and acc tvirs for only those particles in the time window
          loc_maxTemp = maxTemp.(k)[ loc_inds[w] ]
          loc_curTvir = curTvir.(k)[ loc_inds[w] ]
          loc_accTvir = accTvir.(k)[ loc_inds[w] ]
        endif
        
      endelse
      
      ; process outflow members (common)
      galaxy.outRate.(k).cold.num[i]  = nloc
      
      if nloc gt 0 then begin
        ; count mass elements with Tmax below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_maxTemp le sP.TcutVals[j],count_cold,ncomp=count_hot)
          galaxy.outRate.(k).cold.tConst[j,i] = count_cold
          galaxy.outRate.(k).hot.tConst[j,i]  = count_hot
        endfor
        
        for j=0,nVirs-1 do begin
          ; count mass elements with Tmax below Tvir at current time
          w = where(10.0^loc_maxTemp / 10.0^loc_curTvir le sP.TvirVals[j],count_cold,ncomp=count_hot)
          galaxy.outRate.(k).cold.tVirCur[j,i] = count_cold
          galaxy.outRate.(k).hot.tVirCur[j,i]  = count_hot
      
          ; count mass elements with Tmax below Tvir at accretion time
          w = where(10.0^loc_maxTemp / 10.0^loc_accTvir le sP.TvirVals[j],count_cold,ncomp=count_hot)
          galaxy.outRate.(k).cold.tVirAcc[j,i] = count_cold
          galaxy.outRate.(k).hot.tVirAcc[j,i]  = count_hot
        endfor
      endif ;nloc>0
      
      ; --- HALO ACCRETION ---
      
      if sP.newAccRates gt 0 then begin
        if sP.newAccRates eq 1 then begin
          loc_atime = reform( at.accTime[ sP.atIndMode, (wAm.(k))[ loc_inds ] ] )
          
          loc_atime = 1/loc_atime - 1 ; redshift
          loc_atime = redshiftToAgeFlat(loc_atime) * 1e9 ; yr
          
          w = where(curtime - loc_atime le loc_timeWindow, nloc)
        endif else begin
          loc_accFlag = af.accFlag[ (wAm.(k))[loc_inds] ]
          
          w = where(loc_accFlag eq 2, nloc) ; 2 = halo tracer accretion
        endelse
      endif else begin
        ; accretion time defined as 1.0rvir crossing (headed in)
        ; outflow time defined as 1.0rvir crossing (headed out)
        loc_atime = reform( at.accTime[ sP.radIndHaloAcc, (wAm.(k))[ loc_inds ] ] )      
        loc_atime = 1/loc_atime - 1 ; redshift
        loc_atime = redshiftToAgeFlat(loc_atime) * 1e9 ; yr
      
        w = where(curtime - loc_atime le loc_timeWindow, nloc)
      endelse
      
      halo.coldFrac.(k).num[i] = nloc
	halo.coldFrac.total.num[i] += nloc
      
      if nloc gt 0 then begin
        ; maximum past temps, cur and acc tvirs for only those particles in the time window
        loc_maxTemp = maxTemp.(k)[ loc_inds[w] ]
        loc_curTvir = curTvir.(k)[ loc_inds[w] ]
        loc_accTvir = accTvir.(k)[ loc_inds[w] ]
        
        ; count mass elements with Tmax below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_maxTemp le sP.TcutVals[j],count_cold,ncomp=count_hot)
          halo.coldFrac.(k).tConst[j,i]     = float(count_cold) / nloc
	    halo.coldFrac.total.tConst[j,i]   += count_cold
          halo.accRate.(k).cold.tConst[j,i] = count_cold
          halo.accRate.(k).hot.tConst[j,i]  = count_hot
        endfor
        
        for j=0,nVirs-1 do begin
          ; count mass elements with Tmax below Tvir at current time
          w = where(10.0^loc_maxTemp / 10.0^loc_curTvir le sP.TvirVals[j],count_cold,ncomp=count_hot)
          halo.coldFrac.(k).tVirCur[j,i]     = float(count_cold) / nloc
		  halo.coldFrac.total.tVirCur[j,i]   += count_cold
          halo.accRate.(k).cold.tVirCur[j,i] = count_cold
          halo.accRate.(k).hot.tVirCur[j,i]  = count_hot
    
          ; count mass elements with Tmax below Tvir at accretion time
          w = where(10.0^loc_maxTemp / 10.0^loc_accTvir le sP.TvirVals[j],count_cold,ncomp=count_hot)
          halo.coldFrac.(k).tVirAcc[j,i]     = float(count_cold) / nloc
	    halo.coldFrac.total.tVirAcc[j,i]   += count_cold
          halo.accRate.(k).cold.tVirAcc[j,i] = count_cold
          halo.accRate.(k).hot.tVirAcc[j,i]  = count_hot
        endfor
      endif ;nloc>0
      
      ; --- HALO OUTFLOW ---
      
      if sP.newAccRates gt 0 then begin
        if sP.newAccRates eq 1 then begin
          ; no outflow under this model
          nloc = 0
        endif else begin
          ; (3) use the values at targetSnap for outflow calculation
          nloc = 0
          
          if targetSkipFlag eq 0 then begin
            loc_accFlagTarg = af.accFlagTarg[ (wAmTarget.(k))[loc_inds_target] ]
          
            w = where(loc_accFlagTarg eq 2, nloc) ; 2 = halo tracer outflow from targetSnap
          
            if nloc gt 0 then begin
              loc_maxTemp = maxTempTarget.(k)[ loc_inds_target[w] ]
              loc_curTvir = curTvirTarget.(k)[ loc_inds_target[w] ]
              loc_accTvir = accTvirTarget.(k)[ loc_inds_target[w] ]
            endif
          endif
          
        endelse
      endif else begin
        loc_atime = reform( at.outTime[ sP.radIndHaloAcc, (wAm.(k))[ loc_inds ] ] )      
        loc_atime = 1/loc_atime - 1 ; redshift
        loc_atime = redshiftToAgeFlat(loc_atime) * 1e9 ; yr
      
        w = where(curtime - loc_atime le loc_timeWindow, nloc)
        
        if nloc gt 0 then begin
          ; maximum past temps, cur and acc tvirs for only those particles in the time window
          loc_maxTemp = maxTemp.(k)[ loc_inds[w] ]
          loc_curTvir = curTvir.(k)[ loc_inds[w] ]
          loc_accTvir = accTvir.(k)[ loc_inds[w] ]
        endif
      endelse
      
      ; process outflow members (common)
      halo.outRate.(k).cold.num[i]  = nloc
      
      if nloc gt 0 then begin        
        ; count mass elements with Tmax below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_maxTemp le sP.TcutVals[j],count_cold,ncomp=count_hot)
          halo.outRate.(k).cold.tConst[j,i] = count_cold
          halo.outRate.(k).hot.tConst[j,i]  = count_hot
        endfor
        
        for j=0,nVirs-1 do begin
          ; count mass elements with Tmax below Tvir at current time
          w = where(10.0^loc_maxTemp / 10.0^loc_curTvir le sP.TvirVals[j],count_cold,ncomp=count_hot)
          halo.outRate.(k).cold.tVirCur[j,i] = count_cold
          halo.outRate.(k).hot.tVirCur[j,i]  = count_hot
    
          ; count mass elements with Tmax below Tvir at accretion time
          w = where(10.0^loc_maxTemp / 10.0^loc_accTvir le sP.TvirVals[j],count_cold,ncomp=count_hot)
          halo.outRate.(k).cold.tVirAcc[j,i] = count_cold
          halo.outRate.(k).hot.tVirAcc[j,i]  = count_hot
        endfor
      endif ;nloc>0
	  
	; loop over tConst, tVirCur, tVirAcc
	foreach ind,[constInds,tVirInds] do begin
	  ; right now, accRate/outRate hold just particle/tracer counts, convert to mass rates
	  factor = massPerPart * units.UnitMass_in_Msun / loc_timeWindow
		
	  for hotCold=0,1 do begin
	    galaxy.accRate.(k).(hotCold).(ind)[*,i] *= factor
	    galaxy.outRate.(k).(hotCold).(ind)[*,i] *= factor
	    halo.accRate.(k).(hotCold).(ind)[*,i] *= factor
	    halo.outRate.(k).(hotCold).(ind)[*,i] *= factor
		
	    ; compute net=(inflow-outflow) rates by type
	    galaxy.netRate.(k).(hotCold).(ind)[*,i] = $
	    galaxy.accRate.(k).(hotCold).(ind)[*,i] - galaxy.outRate.(k).(hotCold).(ind)[*,i]
	    halo.netRate.(k).(hotCold).(ind)[*,i] = $
	    halo.accRate.(k).(hotCold).(ind)[*,i] - halo.outRate.(k).(hotCold).(ind)[*,i]
		
	    ; calculate total rates (not by type)
	    galaxy.accRate.total.(hotCold).(ind)[*,i] += galaxy.accRate.(k).(hotCold).(ind)[*,i]
	    galaxy.outRate.total.(hotCold).(ind)[*,i] += galaxy.outRate.(k).(hotCold).(ind)[*,i]
	    galaxy.netRate.total.(hotCold).(ind)[*,i] += galaxy.netRate.(k).(hotCold).(ind)[*,i]
	    halo.accRate.total.(hotCold).(ind)[*,i] += halo.accRate.(k).(hotCold).(ind)[*,i]
	    halo.outRate.total.(hotCold).(ind)[*,i] += halo.outRate.(k).(hotCold).(ind)[*,i]
	    halo.netRate.total.(hotCold).(ind)[*,i] += halo.netRate.(k).(hotCold).(ind)[*,i]
	  endfor
      endforeach

    endfor ; n_tags(wAm)
	
    ; loop over tConst, tVirCur, tVirAcc: normalize total cold fractions
    foreach ind,[constInds,tVirInds] do begin
      galaxy.coldFrac.total.(ind)[*,i] /= galaxy.coldFrac.total.num[i]
      halo.coldFrac.total.(ind)[*,i] /= halo.coldFrac.total.num[i]
    endforeach
    
  endfor ; i
  
  if sP.zoomLevel eq 0 then begin
    ; bin fractions into halo mass bins and make median lines
    logMassBins   = [findgen(21)/10 + 9.0,11.1,11.25,11.5,11.75,11.9,13.1]
    logMassNBins  = n_elements(logMassBins)-1
    logMassBinCen = 0.5 * (logMassBins + shift(logMassBins,-1))
    logMassBinCen = logMassBinCen[0:-2]
  
    ; structures to store the binned values
    tempHC = { tConst  : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
               tVirCur : fltarr(nVirs,logMassNbins) + !values.f_nan ,$
               tVirAcc : fltarr(nVirs,logMassNbins) + !values.f_nan  }
  
    template = { hot : tempHC, cold : tempHC }

    foreach typeLabel,typeLabels do $  
      templateTypes = mod_struct(templateTypes, typeLabel, template)
    
    foreach typeLabel,typeLabels do $
      templateCF = mod_struct(templateCF, typeLabel, tempHC)
  
    ; for (1) galaxy and (2) halo, duplicate template for: accretion (inflow), outflow, net, and coldfrac
    galaxyMedian = { accRate : templateTypes, outRate : templateTypes, $
                     netRate : templateTypes, coldFrac : templateCF }

    haloMedian   = { accRate : templateTypes, outRate : templateTypes, $
                     netRate : templateTypes, coldFrac : templateCF }
 				   
    ; calculate median values in bins of halo mass
    for i=0,logMassNbins-1 do begin

    ; since minNumFrac > 0, effectively restrict medians to primary only
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              (galaxy.coldFrac.total.num+halo.coldFrac.total.num) gt minNumFrac,count)
	
    if count eq 0 then continue
	
    print,i,logMassBins[i],count
	
    ; loop over tConst, tVirCur, tVirAcc
    foreach ind,[constInds,tVirInds] do begin
	  
	  ; how many cuts?
	  nCuts = n_elements( galaxy.accRate.gal.hot.(ind)[*,0] )
		
	  for j=0,nCuts-1 do begin
          ; loop over each galaxyCat type (gal,gmem,stars,inter,bhs,total)
          for k=0,n_tags(galaxy.accRate)-1 do begin
          
		; loop  over both hot and cold modes
		for hotCold=0,1 do begin
              ; median accretion, outflow and net rates
              galaxyMedian.accRate.(k).(hotCold).(ind)[j,i] = median(galaxy.accRate.(k).(hotCold).(ind)[j,w])
		  galaxyMedian.outRate.(k).(hotCold).(ind)[j,i] = median(galaxy.outRate.(k).(hotCold).(ind)[j,w])
		  galaxyMedian.netRate.(k).(hotCold).(ind)[j,i] = median(galaxy.netRate.(k).(hotCold).(ind)[j,w])
		  
              haloMedian.accRate.(k).(hotCold).(ind)[j,i] = median(halo.accRate.(k).(hotCold).(ind)[j,w])
		  haloMedian.outRate.(k).(hotCold).(ind)[j,i] = median(halo.outRate.(k).(hotCold).(ind)[j,w])
		  haloMedian.netRate.(k).(hotCold).(ind)[j,i] = median(halo.netRate.(k).(hotCold).(ind)[j,w])
            endfor
		  
		; cold fractions
		galaxyMedian.coldFrac.(k).(ind)[j,i] = median(galaxy.coldFrac.(k).(ind)[j,w])
		haloMedian.coldFrac.(k).(ind)[j,i]   = median(halo.coldFrac.(k).(ind)[j,w])
        endfor ;k
      endfor ;j
    endforeach ; ind
    
    endfor ; logMassNBins
  
  endif else begin ; zoomLevel > 0
    galaxyMedian  = -1
    haloMedian    = -1
    logMassBins   = -1
    logMassBinCen = gcMasses[hInd]
  endelse

  r = {galaxy:galaxy, halo:halo, galaxyMedian:galaxyMedian, haloMedian:haloMedian, $
       radIndGalAcc:sP.radIndGalAcc, radIndHaloAcc:sP.radIndHaloAcc,$
       logMassBins:logMassBins, logMassBinCen:logMassBinCen, TcutVals:sP.TcutVals, TvirVals:sP.TvirVals}

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
