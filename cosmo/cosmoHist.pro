; cosmoHist.pro
; gas accretion project - past history of gas (all routines that step through multiple snapshots)
; dnelson mar.2012

; gasOrigins(): from sP.snap load the galaxy catalog and consider the evolution of all the gas 
;               elements snapshot by snapshot backwards in time until either endingSnap (if specified) or
;               back five snapshots (default). at each step save:
;                 1. the radial distance of each from its primary/secondary parent
;                 2. temperature and entropy
;                 3. implicitly, whether the gas is in a primary or secondary subhalo (corresponding to 
;                    the two radial distances being equal or different)              

function gasOrigins, sP=sP, endingSnap=endingSnap

  compile_opt idl2, hidden, strictarr, strictarrsubs
  print,'think more about this function (gadget only anyways?)' & stop
  ; TODO: make this use the tracers if running on an Arepo SP?

  ; if snap specified, run only one snapshot (and/or just return previous results)
  if n_elements(endingSnap) gt 0 then begin
    saveFilename = sP.derivPath + 'gas.origins.'+str(sP.res)+'.'+str(endingSnap)+'-'+str(sP.snap)+'.sav'
    
    ; results exist, return
    if (file_test(saveFilename)) then begin
      restore,saveFilename
      r = {temp_gal:temp_gal,temp_gmem:temp_gmem,entropy_gal:entropy_gal,entropy_gmem:entropy_gmem,$
           indMatch:indMatch}
      return,r
    endif
    
    ; need to compute, set restricted range of snapshots to process
    snapRange = [sP.snap,endingSnap]
  endif else begin
    ; default config
    numSnapsBack = 5
    
    snapRange = [sP.snap,sP.snap-numSnapsBack]
  endelse  
  
  ; load galaxy catalog at starting redshift
  galcat = galaxyCat(sP=sP)
  
  print,'Loaded  ['+str(n_elements(gc.galaxyIDs))+'] ['+str(n_elements(gc.groupmemIDs))+'] from galCat.'
  
  for m=snapRange[0],snapRange[1],-1 do begin
  
    ; set save filename and check for existence
    sP.snap = m
    saveFilename = sP.derivPath + 'gas.origins.'+str(sP.res)+'.'+str(endingSnap)+'-'+str(sP.snap)+'.sav'
    
    if (file_test(saveFilename)) then begin
      print,'Skipping: '+strmid(saveFilename,strlen(sP.derivPath))
      continue
    endif  
  
    ; load gas IDs and match
    ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    
    ; IMPORTANT! rearrange ids_ind to be in the order of galcat.IDs, need this if we want ids[ids_ind], 
    ; temp[ids_ind], etc to be in the same order as the galaxy catalog id list     
    match,galcat.galaxyIDs,ids,galcat_ind,ids_gal_ind,count=countGal,/sort
    ids_gal_ind = ids_gal_ind[sort(galcat_ind)]
    
    match,galcat.groupmemIDs,ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
    ids_gmem_ind = ids_gmem_ind[sort(galcat_ind)]
    
    ids    = !NULL
    galcat_ind = !NULL
    
    ; load u,nelec and calculate temp of gas
    u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
    nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')

    temp_gal  = convertUtoTemp(u[ids_gal_ind], nelec[ids_gal_ind])
    temp_gmem = convertUtoTemp(u[ids_gmem_ind],nelec[ids_gmem_ind])
    
    nelec = !NULL
    
    ; load gas density to calculate entropy
    dens = loadSnapshotSubset(sP=sP,partType='gas',field='density')
    
    entropy_gal  = calcEntropy(u[ids_gal_ind], dens[ids_gal_ind])
    entropy_gmem = calcEntropy(u[ids_gmem_ind],dens[ids_gmem_ind])
    
    u    = !NULL
    dens = !NULL
    
    ; load the galaxy catalog at this redshift and match IDs from the galCat at the target redshift
    galcatCur = galaxyCat(sP=sP)
    
    match,galcat.galaxyIDs,galcatCur.galaxyIDs,galcat_ind_gal,galcatCur_ind_gal,count=countGal,/sort
    match,galcat.galaxyIDs,galcatCur.groupmemIDs,galcat_ind_gal2,galcatCur_ind_gal2,count=countGal2,/sort
    match,galcat.groupmemIDs,galcatCur.groupmemIDs,galcat_ind_gmem,galcatCur_ind_gmem,count=countGmem,/sort
    match,galcat.groupmemIDs,galcatCur.galaxyIDs,galcat_ind_gmem2,galcatCur_ind_gmem2,count=countGmem2,/sort
    
    galcatCur = !NULL
    
    print,'['+str(endingSnap-m)+'] Matched ['+str(countGal)+' + '+str(countGal2)+'] ['+$
          str(countGmem)+' + '+str(countGmem2)+'].'

    ; keep the match indices relating gas elements in the two galaxy catalogs
    ; NOTE: -1 for empty here    
    indMatch = {galcat_ind_gal:galcat_ind_gal,         galcat_ind_gal2:galcat_ind_gal2,         $
                galcatCur_ind_gal:galcatCur_ind_gal,   galcatCur_ind_gal2:galcatCur_ind_gal2,   $
                galcat_ind_gmem:galcat_ind_gmem,       galcat_ind_gmem2:galcat_ind_gmem2,       $
                galcatCur_ind_gmem:galcatCur_ind_gmem, galcatCur_ind_gmem2:galcatCur_ind_gmem2}
    
    ; save
    save,temp_gal,temp_gmem,entropy_gal,entropy_gmem,indMatch,filename=saveFilename
    print,'    Saved: '+strmid(saveFilename,strlen(sP.derivPath))

  endfor

end

; maxTemps(): find maximum temperature for gas particles in galaxy/group member catalogs at redshift
;             through the redshift range (redshift,zStart] where zStart is typically the start of 
;             the simulation
;
; NOTE: currently temps are only saved for gas in groups at the end of the interval (not all gas)
; saveRedshifts: if not set, then the history trace runs down to sP.snap and is saved
;                if set, save at each redshift as it is reached, running to the lowest

function maxTemps, sP=sP, zStart=zStart, saveRedshifts=saveRedshifts, $
                   loadByGas=loadByGas, loadAllTrGal=loadAllTrGal, loadAllTrGmem=loadAllTrGmem ; load opt

  forward_function cosmoTracerChildren, cosmoTracerVelParents
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; set minimum snapshot (maxmimum redshift)
  if not keyword_set(zStart) then zStart = 30.0
  minSnap = redshiftToSnapnum(zStart,sP=sP)

  ; set maximum snapshot (minimum redshift) (-1 do not include final snapshot in Tmax search)
  if not keyword_set(saveRedshifts) then $
    saveSnaps = [sP.snap - 1]
  if keyword_set(saveRedshifts) then $
    saveSnaps = redshiftToSnapnum(saveRedshifts,sP=sP) - 1

  if keyword_set(saveRedshifts) then stop ; DO NOT USE (in order to know which tracers to follow the
                                          ; properties of we have to first know the galaxy catalog of
                                          ; interest, but for multiple save times the only realistic
                                          ; approach would be to process all tracers and at each save
                                          ; time extract the subset in that galaxy catalog (TODO)
  maxSnap = max(saveSnaps)
  
  ; set saveFilename and check for existence
  if keyword_set(loadByGas) or keyword_set(loadAllTrGal) or keyword_set(loadAllTrGmem) then begin
    if ((keyword_set(loadAllTrGas) or keyword_set(loadAllTrGmem)) and sP.trMCPerCell eq 0) then $
      message,'Cannot load all tracers for SPH type.'
    
    saveTag = ''
    if sP.trMCPerCell eq -1 then saveTag = '.trVel'
    if sP.trMCPerCell gt 0 then  saveTag = '.trMC'
    if sP.trMCPerCell eq 0 then  saveTag = '.SPH'
    
    ; load the condensed results (one per gas particle in galcat)
    if keyword_set(loadByGas) then begin
      saveFilename = sP.derivPath+'maxtemp'+saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                     str(minSnap)+'-'+str(maxSnap)+'.sav'
      
      if not file_test(saveFilename) then message,'Error: Specified maxTemps not found!'
      restore, saveFilename
      return, r
    endif
    
    ; load the results for all tracers (galaxy)
    if keyword_set(loadAllTrGal) then begin
      saveFilename = sP.derivPath + 'maxtemp'+saveTag+'.gal.'+sP.savPrefix+str(sP.res)+'.'+$
                     str(minSnap)+'-'+str(maxSnap)+'.sav'
                     
      if not file_test(saveFilename) then message,'Error: Specified maxTemps gal not found!'
      restore, saveFilename
      return, rtr_gal
    endif
    
    ; load the results for all tracers (groupmem)
    if keyword_set(loadAllTrGmem) then begin
        saveFilename = sP.derivPath + 'maxtemp'+saveTag+'.gmem.'+sP.savPrefix+str(sP.res)+'.'+$
                        str(minSnap)+'-'+str(maxSnap)+'.sav'
                     
      if not file_test(saveFilename) then message,'Error: Specified maxTemps gmem not found!'
      restore, saveFilename
      return, rtr_gmem
    endif
  endif
  
  ; load galaxy/group member catalogs at zMin for gas ids to search for
  galcat = galaxyCat(sP=sP) ;sP.snap is still at zMin
  
  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    print,'Calculating new maxtemp using ( SPH Particles ) res = '+str(sP.res)+$
      ' in range ['+str(minSnap)+'-'+str(maxSnap)+'].'
      
    ; store the main arrays for all tracers as structures so we can write them directly
      r = {maxTemps_gal       : fltarr(n_elements(galcat.galaxyIDs))   ,$
           maxTemps_gmem      : fltarr(n_elements(galcat.groupmemIDs)) ,$
           maxTempTime_gal    : fltarr(n_elements(galcat.galaxyIDs))   ,$
           maxTempTime_gmem   : fltarr(n_elements(galcat.groupmemIDs)) ,$
           maxTempDisp_gal    : 0 ,$
           maxTempDisp_gmem   : 0 ,$ ; no dispersions for sph
           maxTemps_min_gal   : 0 ,$
           maxTemps_min_gmem  : 0 ,$ ; no minimum values for sph
           maxTemps_mean_gal  : 0 ,$
           maxTemps_mean_gmem : 0}   ; no mean values for sph
    
    for m=minSnap,maxSnap,1 do begin
      sP.snap = m
      print,m
      
      ; load gas ids and match to catalog
      ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs, need this if we want ids[ids_ind], 
      ; temp[ids_ind], etc to be in the same order as the group catalog id list    
      match,galcat.galaxyIDs,ids,galcat_ind,ids_gal_ind,count=countGal,/sort
      ids_gal_ind = ids_gal_ind[sort(galcat_ind)]
      
      match,galcat.groupmemIDs,ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
      ids_gmem_ind = ids_gmem_ind[sort(galcat_ind)]
      
      ids        = !NULL
      galcat_ind = !NULL
      
      ; load u,nelec to calculate temperatures
      u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='ne')
      
      u_gal  = u[ids_gal_ind]
      u_gmem = u[ids_gmem_ind]
      u = !NULL
      
      nelec_gal  = nelec[ids_gal_ind]
      nelec_gmem = nelec[ids_gmem_ind]
      nelec = !NULL
      
      temp_gal  = convertUtoTemp(u_gal,nelec_gal)
      temp_gmem = convertUtoTemp(u_gal,nelec_gal)
      
      ; load gas SFR and select off the effective EOS (SFR must be zero)
      sfr = loadSnapshotSubset(sP=sP,partType='gas',field='sfr')
      
      sfr_gal  = sfr[ids_gal_ind]
      sfr_gmem = sfr[ids_gmem_ind]
      
      sfr = !NULL
      
      ; take log of temperatures
      w = where(temp_gal le 0,count_gal)
      if (count_gal ne 0) then temp_gal[w] = 1.0
      temp_gal = alog10(temp_gal)
      
      w = where(temp_gmem le 0,count_gmem)
      if (count_gmem ne 0) then temp_gmem[w] = 1.0
      temp_gmem = alog10(temp_gmem)
      
      ; replace existing values if current snapshot has higher temps (enforce off effective EOS)
      w1 = where(temp_gal gt r.maxTemps_gal and sfr_gal eq 0.0,count1)
      if (count1 gt 0) then begin
        r.maxTemps_gal[w1]    = temp_gal[w1]
        r.maxTempTime_gal[w1] = snapNumToRedshift(sP=sP,/time)
      endif
      
      ; replace existing values if current snapshot has higher temps (group members)
      w2 = where(temp_gmem gt r.maxTemps_gmem and sfr_gmem eq 0.0,count2)
      if (count2 gt 0) then begin
        r.maxTemps_gmem[w2]    = temp_gmem[w2]
        r.maxTempTime_gmem[w2] = snapNumToRedshift(sP=sP,/time)
      endif
      
      ; SAVE?
      if total(sP.snap eq saveSnaps) gt 0 then begin
        ; set savefilename
        w = where(sP.snap eq saveSnaps)
        saveFilename = sP.derivPath + 'maxtemp.SPH.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(saveSnaps[w[0]])+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
          continue
        endif
      
        save,r,filename=saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
      endif ; save
    endfor ;m
  endif
  
  ; MONTE CARLO TRACERS CASE - for each original gas cell, determine some statistics of its
  ; population of tracers and an estimate for the dispersion in those statistics
  ; ------------------------
  if sP.trMCPerCell gt 0 then begin
  
    print,'Calculating new maxtemp using ( TracerMC ) res = '+str(sP.res)+$
      ' in range ['+str(minSnap)+'-'+str(maxSnap)+'].'
      
    ; load gas ids
    gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')

    ; match galcat IDs to gas_ids
    match,galcat.galaxyIDs,gas_ids,galcat_ind,ids_gal_ind,count=countGal,/sort
    ids_gal = gas_ids[ids_gal_ind[sort(galcat_ind)]]
    
    match,galcat.groupmemIDs,gas_ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
    ids_gmem = gas_ids[ids_gmem_ind[sort(galcat_ind)]]
    
    gas_ids = !NULL
    
    ; locate tracer children (indices) of gas id subsets
    galcat_gal_trids  = cosmoTracerChildren(sP=sP, /getInds, gasIDs=ids_gal, child_counts=galcat_gal_cc)
    galcat_gmem_trids = cosmoTracerChildren(sP=sP, /getInds, gasIDs=ids_gmem, child_counts=galcat_gmem_cc)
    
    ; convert tracer children indices to tracer IDs at this zMin
    tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
    
    galcat_gal_trids  = tr_ids[galcat_gal_trids]
    galcat_gmem_trids = tr_ids[galcat_gmem_trids]

    tr_ids   = !NULL
    ids_gal  = !NULL
    ids_gmem = !NULL
    galcat   = !NULL ; not used past this point
    
    ; store the main arrays for all tracers as structures so we can write them directly
    rtr_gal  = { maxTemps      : fltarr(n_elements(galcat_gal_trids))  ,$
                 maxTempTime   : fltarr(n_elements(galcat_gal_trids))  ,$
                 child_counts  : galcat_gal_cc }
    galcat_gal_cc  = !NULL
           
    rtr_gmem = { maxTempTime  : fltarr(n_elements(galcat_gmem_trids)) ,$
                 maxTemps     : fltarr(n_elements(galcat_gmem_trids)) ,$
                 child_counts : galcat_gmem_cc }
    galcat_gmem_cc = !NULL

    for m=minSnap,maxSnap,1 do begin
      sP.snap = m
      print,m  
  
      ; load tracer ids and match to child ids from zMin
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs, need this if we want ids[ids_ind], 
      ; temp[ids_ind], etc to be in the same order as the group catalog id list    
      match,galcat_gal_trids,tr_ids,galcat_ind,trids_gal_ind,count=countGal,/sort
      trids_gal_ind = trids_gal_ind[sort(galcat_ind)]
      
      match,galcat_gmem_trids,tr_ids,galcat_ind,trids_gmem_ind,count=countGmem,/sort
      trids_gmem_ind = trids_gmem_ind[sort(galcat_ind)]
      
      tr_ids     = !NULL
      galcat_ind = !NULL
      
      ; load tracer maximum temperature and sub-snapshot time at this snapshot
      tr_maxtemp = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxtemp')
      
      temp_gal  = tr_maxtemp[trids_gal_ind]  * units.UnitTemp_in_cgs ; tracer output still in unit system
      temp_gmem = tr_maxtemp[trids_gmem_ind] * units.UnitTemp_in_cgs
      tr_maxtemp = !NULL
      
      tr_maxtemp_time = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxtemp_time')
      
      temp_time_gal  = tr_maxtemp_time[trids_gal_ind]
      temp_time_gmem = tr_maxtemp_time[trids_gmem_ind]
      tr_maxtemp_time = !NULL
      
      ; take log of temperatures
      w = where(temp_gal le 0,count_gal)
      if (count_gal ne 0) then temp_gal[w] = 1.0
      temp_gal = alog10(temp_gal)
      
      w = where(temp_gmem le 0,count_gmem)
      if (count_gmem ne 0) then temp_gmem[w] = 1.0
      temp_gmem = alog10(temp_gmem)
      
      ; replace existing values if current snapshot has higher temps (galaxy members)
      w1 = where(temp_gal gt rtr_gal.maxTemps,count1)
      if (count1 gt 0) then begin
        rtr_gal.maxTemps[w1]    = temp_gal[w1]
        rtr_gal.maxTempTime[w1] = temp_time_gal[w1] ; sub-snapshot timing
      endif
      
      ; replace existing values if current snapshot has higher temps (group members)
      w2 = where(temp_gmem gt rtr_gmem.maxTemps,count2)
      if (count2 gt 0) then begin
        rtr_gmem.maxTemps[w2]    = temp_gmem[w2]
        rtr_gmem.maxTempTime[w2] = temp_time_gmem[w2] ; sub-snapshot timing
      endif
      
      temp_gal  = !NULL
      temp_gmem = !NULL
      temp_time_gal  = !NULL
      temp_time_gmem = !NULL
      
      ; SAVE?
      if total(sP.snap eq saveSnaps) gt 0 then begin
        ; (1) full tracer information (galaxy members) - set savefilename
        w = where(sP.snap eq saveSnaps)
        saveFilename = sP.derivPath + 'maxtemp.trMC.gal.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(saveSnaps[w[0]])+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
        endif else begin
          ; output all gal/gmem tracer values at this point and the child counts so they can be 
          ; matched back up to their parents later
          save,rtr_gal,filename=saveFilename
          print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        endelse
        
        ; (2) full tracer information (group members) - set savefilename
        w = where(sP.snap eq saveSnaps)
        saveFilename = sP.derivPath + 'maxtemp.trMC.gmem.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(saveSnaps[w[0]])+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
        endif else begin
          ; output all gal/gmem tracer values at this point and the child counts so they can be 
          ; matched back up to their parents later
          save,rtr_gmem,filename=saveFilename
          print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        endelse
        
        ; (3) values condensed to gas parents - set savefilename
        w = where(sP.snap eq saveSnaps)
        saveFilename = sP.derivPath + 'maxtemp.trMC.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(saveSnaps[w[0]])+'.sav'
                       
        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
          continue
        endif
        
        ; prepare output
        offset = 0L
        
        r = {maxTemps_gal       : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; max(maxTemps)
             maxTempTime_gal    : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; time of max(maxTemps)
             maxTempDisp_gal    : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; stddev(maxTemps)
             maxTemps_gmem      : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTempTime_gmem   : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTempDisp_gmem   : fltarr(n_elements(rtr_gmem.child_counts)) ,$
                                                                      $
             maxTemps_min_gal   : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; min(maxTemps)
             maxTemps_mean_gal  : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; mean(maxTemps)
             maxTemps_min_gmem  : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTemps_mean_gmem : fltarr(n_elements(rtr_gmem.child_counts)) ,$
                                                                      $
             maxTempTime_min_gal   : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; time of min(maxTemps)
             maxTempTime_mean_gal  : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; mean(times maxTemps)
             maxTempTime_min_gmem  : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTempTime_mean_gmem : fltarr(n_elements(rtr_gmem.child_counts))  }
        
        for i=0,n_elements(rtr_gal.child_counts)-1 do begin
          if rtr_gal.child_counts[i] gt 0 then begin
            ; for each gas cell, collect its child tracers
            childInds = lindgen(rtr_gal.child_counts[i]) + offset
            
            locMaxTemps    = rtr_gal.maxTemps[childInds]
            locMaxTempTime = rtr_gal.maxTempTime[childInds] ; this is really time (scale factor) for tracers
            
            ; save the maximum value reached by any child tracer and the stddev of the population maxima
            r.maxTemps_gal[i]    = max(locMaxTemps,indmax)
            r.maxTempTime_gal[i] = locMaxTempTime[indmax]
            r.maxTempDisp_gal[i] = stddev(locMaxTemps)
            
            ; save minimum and means (and associated times)
            r.maxTemps_min_gal[i]     = min(locMaxTemps,indmin)
            r.maxTempTime_min_gal[i]  = locMaxTempTime[indmin]
            r.maxTemps_mean_gal[i]    = mean(locMaxTemps)
            r.maxTempTime_mean_gal[i] = mean(locMaxTempTime)
            
            offset += rtr_gal.child_counts[i]
          endif
        endfor
        
        offset = 0L
        
        for i=0,n_elements(rtr_gmem.child_counts)-1 do begin
          if rtr_gmem.child_counts[i] gt 0 then begin
            ; for each gas cell, collect its child tracers
            childInds = lindgen(rtr_gmem.child_counts[i]) + offset
            
            locMaxTemps    = rtr_gmem.maxTemps[childInds]
            locMaxTempTime = rtr_gmem.maxTempTime[childInds]
            
            ; save the maximum value reached by any child tracer and the stddev of the population maxima
            r.maxTemps_gmem[i]    = max(locMaxTemps,indmax)
            r.maxTempTime_gmem[i] = locMaxTempTime[indmax]
            r.maxTempDisp_gmem[i] = stddev(locMaxTemps)
            
            ; save minimum and means (and associated times)
            r.maxTemps_min_gmem[i]     = min(locMaxTemps,indmin)
            r.maxTempTime_min_gmem[i]  = locMaxTempTime[indmin]
            r.maxTemps_mean_gmem[i]    = mean(locMaxTemps)
            r.maxTempTime_mean_gmem[i] = mean(locMaxTempTime)
            
            offset += rtr_gmem.child_counts[i]
          endif
        endfor
        
        ; save for future lookups        
        save,r,filename=saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        r = !NULL
      endif ;save
    endfor ;m  
  endif
  
  ; VELOCITY TRACERS case - will be similar to above since there could be multiple
  ; ---------------------
  if sP.trMCPerCell eq -1 then begin
  
    print,'Calculating new maxtemp using ( TracerVEL ) res = '+str(sP.res)+$
      ' in range ['+str(minSnap)+'-'+str(maxSnap)+'].'

    ; load gas ids
    gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')

    ; match galcat IDs to gas_ids
    match,galcat.galaxyIDs,gas_ids,galcat_ind,ids_gal_ind,count=countGal,/sort
    inds_gal = ids_gal_ind[sort(galcat_ind)]
    
    match,galcat.groupmemIDs,gas_ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
    inds_gmem = ids_gmem_ind[sort(galcat_ind)]

    gas_ids = !NULL
    
    ; locate tracer children (indices) of gas id subsets
    galcat_gal_trids  = cosmoTracerVelChildren(sP=sP,/getInds,gasInds=inds_gal,child_counts=galcat_gal_cc)
    galcat_gmem_trids = cosmoTracerVelChildren(sP=sP,/getInds,gasInds=inds_gmem,child_counts=galcat_gmem_cc)
    
    ; convert tracer children indices to tracer IDs at this zMin
    tr_ids = loadSnapshotSubset(sP=sP,partType='tracerVel',field='ids')
    
    galcat_gal_trids  = tr_ids[galcat_gal_trids]
    galcat_gmem_trids = tr_ids[galcat_gmem_trids]

    tr_ids   = !NULL
    inds_gal  = !NULL
    inds_gmem = !NULL
    galcat   = !NULL ; not used past this point

    ; store the main arrays for all tracers as structures so we can write them directly
    rtr_gal  = { maxTemps      : fltarr(n_elements(galcat_gal_trids))  ,$
                 maxTempTime   : fltarr(n_elements(galcat_gal_trids))  ,$
                 child_counts  : galcat_gal_cc }
    galcat_gal_cc  = !NULL
           
    rtr_gmem = { maxTempTime  : fltarr(n_elements(galcat_gmem_trids)) ,$
                 maxTemps     : fltarr(n_elements(galcat_gmem_trids)) ,$
                 child_counts : galcat_gmem_cc }
    galcat_gmem_cc = !NULL
  
    for m=minSnap,maxSnap,1 do begin
      sP.snap = m
      print,m  
  
      ; load tracer ids and match to child ids from zMin
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerVel',field='ids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs, need this if we want ids[ids_ind], 
      ; temp[ids_ind], etc to be in the same order as the group catalog id list    
      match,galcat_gal_trids,tr_ids,galcat_ind,trids_gal_ind,count=countGal,/sort
      trids_gal_ind = trids_gal_ind[sort(galcat_ind)]
      
      match,galcat_gmem_trids,tr_ids,galcat_ind,trids_gmem_ind,count=countGmem,/sort
      trids_gmem_ind = trids_gmem_ind[sort(galcat_ind)]
      
      tr_ids     = !NULL
      galcat_ind = !NULL
      
      ; load tracer maximum temperature and sub-snapshot time at this snapshot
      tr_maxtemp = loadSnapshotSubset(sP=sP,partType='tracerVel',field='tracer_maxtemp')
      
      temp_gal  = tr_maxtemp[trids_gal_ind]  * units.UnitTemp_in_cgs ; tracer output still in unit system
      temp_gmem = tr_maxtemp[trids_gmem_ind] * units.UnitTemp_in_cgs
      tr_maxtemp = !NULL
      
      tr_maxtemp_time = loadSnapshotSubset(sP=sP,partType='tracerVel',field='tracer_maxtemp_time')
      
      temp_time_gal  = tr_maxtemp_time[trids_gal_ind]
      temp_time_gmem = tr_maxtemp_time[trids_gmem_ind]
      tr_maxtemp_time = !NULL
      
      ; take log of temperatures
      w = where(temp_gal le 0,count_gal)
      if (count_gal ne 0) then temp_gal[w] = 1.0
      temp_gal = alog10(temp_gal)
      
      w = where(temp_gmem le 0,count_gmem)
      if (count_gmem ne 0) then temp_gmem[w] = 1.0
      temp_gmem = alog10(temp_gmem)
      
      ; replace existing values if current snapshot has higher temps (galaxy members)
      w1 = where(temp_gal gt rtr_gal.maxTemps,count1)
      if (count1 gt 0) then begin
        rtr_gal.maxTemps[w1]    = temp_gal[w1]
        rtr_gal.maxTempTime[w1] = temp_time_gal[w1] ; sub-snapshot timing
      endif
      
      ; replace existing values if current snapshot has higher temps (group members)
      w2 = where(temp_gmem gt rtr_gmem.maxTemps,count2)
      if (count2 gt 0) then begin
        rtr_gmem.maxTemps[w2]    = temp_gmem[w2]
        rtr_gmem.maxTempTime[w2] = temp_time_gmem[w2] ; sub-snapshot timing
      endif
      
      temp_gal  = !NULL
      temp_gmem = !NULL
      temp_time_gal  = !NULL
      temp_time_gmem = !NULL
     
      ; SAVE?
      if total(sP.snap eq saveSnaps) gt 0 then begin
        ; (1) full tracer information (galaxy members) - set savefilename
        w = where(sP.snap eq saveSnaps)
        saveFilename = sP.derivPath + 'maxtemp.trVel.gal.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(saveSnaps[w[0]])+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
        endif else begin
          ; output all gal/gmem tracer values at this point and the child counts so they can be 
          ; matched back up to their parents later
          save,rtr_gal,filename=saveFilename
          print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        endelse
        
        ; (2) full tracer information (group members) - set savefilename
        w = where(sP.snap eq saveSnaps)
        saveFilename = sP.derivPath + 'maxtemp.trVel.gmem.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(saveSnaps[w[0]])+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
        endif else begin
          ; output all gal/gmem tracer values at this point and the child counts so they can be 
          ; matched back up to their parents later
          save,rtr_gmem,filename=saveFilename
          print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        endelse
        
        ; (3) values condensed to gas parents - set savefilename
        w = where(sP.snap eq saveSnaps)
        saveFilename = sP.derivPath + 'maxtemp.trVel.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(saveSnaps[w[0]])+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
          continue
        endif
        
        ; prepare output
        offset = 0L
        
        r = {maxTemps_gal       : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; max(maxTemps)
             maxTempTime_gal    : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; time of max(maxTemps)
             maxTempDisp_gal    : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; stddev(maxTemps)
             maxTemps_gmem      : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTempTime_gmem   : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTempDisp_gmem   : fltarr(n_elements(rtr_gmem.child_counts)) ,$
                                                                      $
             maxTemps_min_gal   : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; min(maxTemps)
             maxTemps_mean_gal  : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; mean(maxTemps)
             maxTemps_min_gmem  : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTemps_mean_gmem : fltarr(n_elements(rtr_gmem.child_counts)) ,$
                                                                      $
             maxTempTime_min_gal   : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; time of min(maxTemps)
             maxTempTime_mean_gal  : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; mean(times maxTemps)
             maxTempTime_min_gmem  : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTempTime_mean_gmem : fltarr(n_elements(rtr_gmem.child_counts))  }
        
        for i=0,n_elements(rtr_gal.child_counts)-1 do begin
          if rtr_gal.child_counts[i] gt 0 then begin
            ; for each gas cell, collect its child tracers
            childInds = lindgen(rtr_gal.child_counts[i]) + offset
            
            locMaxTemps    = rtr_gal.maxTemps[childInds]
            locMaxTempTime = rtr_gal.maxTempTime[childInds] ; this is really time (scale factor) for tracers
            
            ; save the maximum value reached by any child tracer and the stddev of the population maxima
            r.maxTemps_gal[i]    = max(locMaxTemps,indmax)
            r.maxTempTime_gal[i] = locMaxTempTime[indmax]
            r.maxTempDisp_gal[i] = stddev(locMaxTemps)
            
            ; save minimum and means (and associated times)
            r.maxTemps_min_gal[i]     = min(locMaxTemps,indmin)
            r.maxTempTime_min_gal[i]  = locMaxTempTime[indmin]
            r.maxTemps_mean_gal[i]    = mean(locMaxTemps)
            r.maxTempTime_mean_gal[i] = mean(locMaxTempTime)
            
            offset += rtr_gal.child_counts[i]
          endif
        endfor
        
        offset = 0L
        
        for i=0,n_elements(rtr_gmem.child_counts)-1 do begin
          if rtr_gmem.child_counts[i] gt 0 then begin
            ; for each gas cell, collect its child tracers
            childInds = lindgen(rtr_gmem.child_counts[i]) + offset
            
            locMaxTemps    = rtr_gmem.maxTemps[childInds]
            locMaxTempTime = rtr_gmem.maxTempTime[childInds]
            
            ; save the maximum value reached by any child tracer and the stddev of the population maxima
            r.maxTemps_gmem[i]    = max(locMaxTemps,indmax)
            r.maxTempTime_gmem[i] = locMaxTempTime[indmax]
            r.maxTempDisp_gmem[i] = stddev(locMaxTemps)
            
            ; save minimum and means (and associated times)
            r.maxTemps_min_gmem[i]     = min(locMaxTemps,indmin)
            r.maxTempTime_min_gmem[i]  = locMaxTempTime[indmin]
            r.maxTemps_mean_gmem[i]    = mean(locMaxTemps)
            r.maxTempTime_mean_gmem[i] = mean(locMaxTempTime)
            
            offset += rtr_gmem.child_counts[i]
          endif
        endfor
        
        ; save for future lookups        
        save,r,filename=saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        r = !NULL
      endif ;save
    endfor ;m
  endif

end