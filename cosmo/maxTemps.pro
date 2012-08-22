; maxTemps.pro
; gas accretion project - past temperature history of gas
; dnelson aug.2012

; -----------------------------------------------------------------------------------------------------
; maxTemps(): find maximum temperature for gas particles in galaxy/group member catalogs at redshift
;             through the redshift range (redshift,zStart] where zStart is typically the start of 
;             the simulation
;
; NOTE: currently temps are only saved for gas in groups at the end of the interval (not all gas)
; -----------------------------------------------------------------------------------------------------

function maxTemps, sP=sP, zStart=zStart, restart=restart, $
                   loadByGas=loadByGas, $ ; load options
                   loadAllTrGal=loadAllTrGal, loadAllTrGmem=loadAllTrGmem, loadAllTrStars=loadAllTrStars

  forward_function cosmoTracerChildren, cosmoTracerVelParents
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; set minimum snapshot (maxmimum redshift)
  if not keyword_set(zStart) then zStart = snapNumToRedshift(snap=0,sP=sP)
  minSnap = redshiftToSnapnum(zStart,sP=sP)

  ; set maximum snapshot (minimum redshift) (-1 do not include final snapshot in Tmax search)
  maxSnap = sP.snap - 1
  
  snapRange = [minSnap,maxSnap]
  
  ; set saveFilename and check for existence
  if (keyword_set(loadByGas) or keyword_set(loadAllTrGal) or keyword_set(loadAllTrGmem) or $
      keyword_set(loadAllTrStars)) then begin
    if ((keyword_set(loadAllTrGas) or keyword_set(loadAllTrGmem) or keyword_set(loadAllTrStars)) and $
      sP.trMCPerCell eq 0) then message,'Cannot load all tracers for SPH type.'
    
    ; load the condensed results (one per gas particle in galcat)
    if keyword_set(loadByGas) then begin
      saveFilename = sP.derivPath+'maxtemp.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                     str(minSnap)+'-'+str(maxSnap)+'.sav'
      
      if not file_test(saveFilename) then message,'Error: Specified maxTemps not found!'
      restore, saveFilename
      return, r
    endif
    
    ; load the results for all tracers (galaxy)
    if keyword_set(loadAllTrGal) then begin
      saveFilename = sP.derivPath + 'maxtemp.'+sP.saveTag+'.gal.'+sP.savPrefix+str(sP.res)+'.'+$
                     str(minSnap)+'-'+str(maxSnap)+'.sav'
                     
      if not file_test(saveFilename) then message,'Error: Specified maxTemps gal not found!'
      restore, saveFilename
      return, rtr_gal
    endif
    
    ; load the results for all tracers (groupmem)
    if keyword_set(loadAllTrGmem) then begin
        saveFilename = sP.derivPath + 'maxtemp.'+sP.saveTag+'.gmem.'+sP.savPrefix+str(sP.res)+'.'+$
                        str(minSnap)+'-'+str(maxSnap)+'.sav'
                     
      if not file_test(saveFilename) then message,'Error: Specified maxTemps gmem not found!'
      restore, saveFilename
      return, rtr_gmem
    endif
    
    ; load the results for all tracers (stars)
    if keyword_set(loadAllTrStars) then begin
        saveFilename = sP.derivPath + 'maxtemp.'+sP.saveTag+'.stars.'+sP.savPrefix+str(sP.res)+'.'+$
                        str(minSnap)+'-'+str(maxSnap)+'.sav'
                     
      if not file_test(saveFilename) then message,'Error: Specified maxTemps stars not found!'
      restore, saveFilename
      return, rtr_stars
    endif
  endif
  
  resFilename = sP.derivPath + 'maxtemp.'+sP.saveTag+'.restart.'+sP.savPrefix+str(sP.res)+'.'+$
                        str(minSnap)+'-'+str(maxSnap)+'.sav'
                        
  ; check for maxTempsAll save (in this case we don't have to loop over any snapshots, but load directly)
  maxTempsAllFlag = 0
  maxTempsAllSaveFilename = sP.derivPath + 'maxTempAll.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+str(maxSnap)+'.sav'

  if file_test(maxTempsAllSaveFilename) then maxTempsAllFlag = 1
  
  ; load galaxy/group member catalogs at zMin for gas ids to search for
  galcat = galaxyCat(sP=sP) ;sP.snap is still at zMin
  
  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
  
    print,'Calculating new maxtemp using ( SPH Particles ) res = '+str(sP.res)+$
      ' in range ['+str(minSnap)+'-'+str(maxSnap)+'].'
      
    if ~file_test(resFilename) then begin ; no restart  
        ; store the main arrays for all tracers as structures so we can write them directly
        r = {maxTemps_gal       : fltarr(n_elements(galcat.galaxyIDs))   ,$
             maxTemps_gmem      : fltarr(n_elements(galcat.groupmemIDs)) ,$
             maxTemps_stars     : fltarr(n_elements(galcat.stellarIDs))  ,$
             maxTempTime_gal    : fltarr(n_elements(galcat.galaxyIDs))   ,$
             maxTempTime_gmem   : fltarr(n_elements(galcat.groupmemIDs)) ,$
             maxTempTime_stars  : fltarr(n_elements(galcat.stellarIDs))  ,$
             maxTempDisp_gal    : 0 ,$
             maxTempDisp_gmem   : 0 ,$ ; no dispersions for sph
             maxTempDisp_stars  : 0 ,$
             maxTemps_min_gal   : 0 ,$
             maxTemps_min_gmem  : 0 ,$ ; no minimum values for sph
             maxTemps_min_stars : 0 ,$
             maxTemps_mean_gal  : 0 ,$ ; no mean values for sph
             maxTemps_mean_gmem : 0 ,$
             maxTemps_mean_stars : 0 }
    endif else begin
      ; restart
      if ~keyword_set(restart) then message,'Error: Restart file exists but restart not requested.'
      restore,resFilename,/verbose
      snapRange[0] = m
    endelse    
    
    for m=snapRange[0],snapRange[1],1 do begin
      sP.snap = m
      print,m
      
      ; save restart?
      if m mod 10 eq 0 and m gt minSnap and keyword_set(restart) then begin
        print,' --- Writing restart! ---'
        save,r,m,filename=resFilename
        print,' --- Done! ---'
      endif
      
      ; load gas ids and match to catalog
      ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs, need this if we want ids[ids_ind], 
      ; temp[ids_ind], etc to be in the same order as the group catalog id list    
      match,galcat.galaxyIDs,ids,galcat_ind,ids_gal_ind,count=countGal,/sort
      ids_gal_ind = ids_gal_ind[sort(galcat_ind)]
      
      match,galcat.groupmemIDs,ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
      ids_gmem_ind = ids_gmem_ind[sort(galcat_ind)]
      
      ; match galcat star ids to current gas ids
      match,galcat.stellarIDs,ids,galcat_ind,ids_stars_ind,count=countStars,/sort
      ids_stars_ind = ids_stars_ind[sort(galcat_ind)]
      
      ids        = !NULL
      galcat_ind = !NULL
      
      ; load u,nelec to calculate temperatures
      u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='ne')
      
      u_gal  = u[ids_gal_ind]
      u_gmem = u[ids_gmem_ind]
      if countStars gt 0 then u_stars = u[ids_stars_ind]
      u = !NULL
      
      nelec_gal  = nelec[ids_gal_ind]
      nelec_gmem = nelec[ids_gmem_ind]
      if countStars gt 0 then nelec_stars = nelec[ids_stars_ind]
      nelec = !NULL
      
      temp_gal  = convertUtoTemp(u_gal,nelec_gal,/log)
      temp_gmem = convertUtoTemp(u_gmem,nelec_gmem,/log)
      if countStars gt 0 then temp_stars = convertUtoTemp(u_stars,nelec_stars,/log)
      
      ; load gas SFR and select off the effective EOS (SFR must be zero)
      sfr = loadSnapshotSubset(sP=sP,partType='gas',field='sfr')
      
      sfr_gal  = sfr[ids_gal_ind]
      sfr_gmem = sfr[ids_gmem_ind]
      if countStars gt 0 then sfr_stars = sfr[ids_stars_ind]
      
      sfr = !NULL
      
      ; replace existing values if current snapshot has higher temps (enforce off effective EOS)
      w1 = where(temp_gal gt r.maxTemps_gal and sfr_gal eq 0.0,count1)
      if count1 gt 0 then begin
        r.maxTemps_gal[w1]    = temp_gal[w1]
        r.maxTempTime_gal[w1] = snapNumToRedshift(sP=sP,/time)
      endif
      
      ; replace existing values if current snapshot has higher temps (group members)
      w2 = where(temp_gmem gt r.maxTemps_gmem and sfr_gmem eq 0.0,count2)
      if count2 gt 0 then begin
        r.maxTemps_gmem[w2]    = temp_gmem[w2]
        r.maxTempTime_gmem[w2] = snapNumToRedshift(sP=sP,/time)
      endif
      
      ; if galcat stars are gas at this snapshot, and off eEOS, replace with higher temps
      if countStars gt 0 then begin
        w3 = where(temp_stars gt r.maxTemps_stars and sfr_stars eq 0.0,count3)
        if count3 gt 0 then begin
          r.maxTemps_stars[w3] = temp_stars[w3]
          r.maxTempTime_stars[w3] = snapNumToRedshift(sP=sP,/time)
        endif
      endif
      
      ; SAVE?
      if sP.snap eq maxSnap then begin
        ; set savefilename
        saveFilename = sP.derivPath + 'maxtemp.SPH.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(maxSnap)+'.sav'

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
      
    if ~file_test(resFilename) then begin ; no restart   
      ; load gas ids
      gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
  
      ; match galcat IDs to gas_ids
      match,galcat.galaxyIDs,gas_ids,galcat_ind,ids_gal_ind,count=countGal,/sort
      ids_gal = gas_ids[ids_gal_ind[sort(galcat_ind)]]
      
      match,galcat.groupmemIDs,gas_ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
      ids_gmem = gas_ids[ids_gmem_ind[sort(galcat_ind)]]
      
      gas_ids = !NULL
      
      ; load star ids and match
      star_ids = loadSnapshotSubset(sP=sP,partType='stars',field='ids')
      
      match,galcat.stellarIDs,star_ids,galcat_ind,ids_stars_ind,count=countStars,/sort
      ids_stars = star_ids[ids_stars_ind[sort(galcat_ind)]]
      
      star_ids = !NULL
      
      ; locate tracer children (indices) of gas id subsets
      galcat_gal_trids   = cosmoTracerChildren(sP=sP, /getInds, gasIDs=ids_gal, child_counts=galcat_gal_cc)
      galcat_gmem_trids  = cosmoTracerChildren(sP=sP, /getInds, gasIDs=ids_gmem, child_counts=galcat_gmem_cc)
      galcat_stars_trids = cosmoTracerChildren(sP=sP, /getInds, starIDs=ids_stars, child_counts=galcat_stars_cc)
      
      ; convert tracer children indices to tracer IDs at this zMin
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
      
      galcat_gal_trids   = tr_ids[galcat_gal_trids]
      galcat_gmem_trids  = tr_ids[galcat_gmem_trids]
      galcat_stars_trids = tr_ids[galcat_stars_trids]
  
      tr_ids    = !NULL
      ids_gal   = !NULL
      ids_gmem  = !NULL
      ids_stars = !NULL
      galcat    = !NULL ; not used past this point
      
      ; store the main arrays for all tracers as structures so we can write them directly
      rtr_gal  = { maxTemps      : fltarr(n_elements(galcat_gal_trids))  ,$
                   maxTempTime   : fltarr(n_elements(galcat_gal_trids))  ,$
                   child_counts  : galcat_gal_cc }
      galcat_gal_cc  = !NULL
             
      rtr_gmem = { maxTempTime  : fltarr(n_elements(galcat_gmem_trids)) ,$
                   maxTemps     : fltarr(n_elements(galcat_gmem_trids)) ,$
                   child_counts : galcat_gmem_cc }
      galcat_gmem_cc = !NULL
      
      rtr_stars = { maxTempTime  : fltarr(n_elements(galcat_stars_trids)) ,$
                    maxTemps     : fltarr(n_elements(galcat_stars_trids)) ,$
                    child_counts : galcat_stars_cc }
      galcat_stars_cc = !NULL
    endif else begin
      ; restart
      if ~keyword_set(restart) then message,'Error: Restart file exists but restart not requested.'
      restore,resFilename,/verbose
      snapRange[0] = m
    endelse 

    for m=snapRange[0],snapRange[1],1 do begin
      sP.snap = m
      print,m  
  
      ; save restart?
      if m mod 10 eq 0 and m gt minSnap and keyword_set(restart) then begin
        print,' --- Writing restart! ---'
        save,rtr_gal,rtr_gmem,rtr_stars,galcat_gal_trids,galcat_gmem_trids,galcat_stars_trids,m,filename=resFilename
        print,' --- Done! ---'
      endif
      
      if maxTempsAllFlag eq 0 then begin ; normal
  
        ; load tracer ids and match to child ids from zMin
        tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
        
        ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs, need this if we want ids[ids_ind], 
        ; temp[ids_ind], etc to be in the same order as the group catalog id list    
        match,galcat_gal_trids,tr_ids,galcat_ind,trids_gal_ind,count=countGal,/sort
        trids_gal_ind = trids_gal_ind[sort(galcat_ind)]
        
        match,galcat_gmem_trids,tr_ids,galcat_ind,trids_gmem_ind,count=countGmem,/sort
        trids_gmem_ind = trids_gmem_ind[sort(galcat_ind)]
        
        match,galcat_stars_trids,tr_ids,galcat_ind,trids_stars_ind,count=countStars,/sort
        trids_stars_ind = trids_stars_ind[sort(galcat_ind)]
        
        tr_ids     = !NULL
        galcat_ind = !NULL
        
        ; load tracer maximum temperature and sub-snapshot time at this snapshot
        tr_maxtemp = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxtemp')
        
        ; tracerMC maxtemp field changed to Kelvin in recent code, tracerVEL still in unit system
        w = where(tr_maxtemp eq 0,count)
        if count gt 0 then tr_maxtemp[w] = 1.0
        tr_maxtemp = alog10(tr_maxtemp)
        w = !NULL
        
        temp_gal   = tr_maxtemp[trids_gal_ind]
        temp_gmem  = tr_maxtemp[trids_gmem_ind]
        temp_stars = tr_maxtemp[trids_stars_ind]
        tr_maxtemp = !NULL
        
        tr_maxtemp_time = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxtemp_time')
        
        temp_time_gal   = tr_maxtemp_time[trids_gal_ind]
        temp_time_gmem  = tr_maxtemp_time[trids_gmem_ind]
        temp_time_stars = tr_maxtemp_time[trids_stars_ind]
        tr_maxtemp_time = !NULL
        
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
        
        ; replace existing values if current snapshot has higher temps (stars)
        ; note: if tracers are still inside a star particle, their maxtemp entry will be zero
        w3 = where(temp_stars gt rtr_stars.maxTemps,count3)
        if (count3 gt 0) then begin
          rtr_stars.maxTemps[w3]    = temp_stars[w3]
          rtr_stars.maxTempTime[w3] = temp_time_stars[w3] ; sub-snapshot timing
        endif
        
        temp_gal   = !NULL
        temp_gmem  = !NULL
        temp_stars = !NULL
        temp_time_gal   = !NULL
        temp_time_gmem  = !NULL
        temp_time_stars = !NULL
        
      endif else begin ; we have a maxTempsAll save, just load and move immediately to save in correct format
        
        ; load tracer ids and match to maxTempsAll save order (sorted ascending)
        tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
        tr_ids = tr_ids[sort(tr_ids)]
        
        ; IMPORTANT! rearrange trids_x_ind to be in the order of galcat_x_trids
        match,galcat_gal_trids,tr_ids,galcat_ind,trids_gal_ind,count=countGal,/sort
        trids_gal_ind = trids_gal_ind[sort(galcat_ind)]
        
        match,galcat_gmem_trids,tr_ids,galcat_ind,trids_gmem_ind,count=countGmem,/sort
        trids_gmem_ind = trids_gmem_ind[sort(galcat_ind)]
        
        match,galcat_stars_trids,tr_ids,galcat_ind,trids_stars_ind,count=countStars,/sort
        trids_stars_ind = trids_stars_ind[sort(galcat_ind)]
        
        tr_ids = !NULL
        galcat_ind = !NULL
        
        if countGal ne n_elements(galcat_gal_trids) or $
           countGmem ne n_elements(galcat_gmem_trids) or $
           countStars ne n_elements(galcat_stars_trids) then message,'Error: maxTempsAll recovery failed.'
        
        ; fill rtr_gal, rtr_gmem, rtr_stars
        restore,maxTempsAllSaveFilename,/verbose
        
        rtr_gal.maxTemps      = rtr_all.maxTemps[trids_gal_ind]
        rtr_gal.maxTempTime   = rtr_all.maxTempTime[trids_gal_ind]
        rtr_gmem.maxTemps     = rtr_all.maxTemps[trids_gmem_ind]
        rtr_gmem.maxTempTime  = rtr_all.maxTempTime[trids_gmem_ind]
        rtr_stars.maxTemps    = rtr_all.maxTemps[trids_stars_ind]
        rtr_stars.maxTempTime = rtr_all.maxTempTime[trids_stars_ind]
        
        rtr_all = !NULL
        
        ; modify current snapshot
        sP.snap = maxSnap
        m = maxSnap
        
      endelse
      
      ; SAVE?
      if sP.snap eq maxSnap then begin
        ; (1) full tracer information (galaxy members) - set savefilename
        saveFilename = sP.derivPath + 'maxtemp.trMC.gal.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(maxSnap)+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
        endif else begin
          ; output all gal/gmem tracer values at this point and the child counts so they can be 
          ; matched back up to their parents later
          save,rtr_gal,filename=saveFilename
          print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        endelse
        
        ; (2) full tracer information (group members) - set savefilename
        saveFilename = sP.derivPath + 'maxtemp.trMC.gmem.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(maxSnap)+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
        endif else begin
          ; output all gal/gmem tracer values at this point and the child counts so they can be 
          ; matched back up to their parents later
          save,rtr_gmem,filename=saveFilename
          print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        endelse
        
        ; (2) full tracer information (group members) - set savefilename
        saveFilename = sP.derivPath + 'maxtemp.trMC.stars.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(maxSnap)+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
        endif else begin
          ; output all star tracer values at this point and the child counts so they can be 
          ; matched back up to their parents later
          save,rtr_stars,filename=saveFilename
          print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        endelse
        
        ; (3) values condensed to gas parents - set savefilename
        saveFilename = sP.derivPath + 'maxtemp.trMC.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(maxSnap)+'.sav'
                       
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
             maxTemps_stars     : fltarr(n_elements(rtr_stars.child_counts)) ,$
             maxTempTime_stars  : fltarr(n_elements(rtr_stars.child_counts)) ,$
             maxTempDisp_stars  : fltarr(n_elements(rtr_stars.child_counts)) ,$
                                                                      $
             maxTemps_min_gal    : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; min(maxTemps)
             maxTemps_mean_gal   : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; mean(maxTemps)
             maxTemps_min_gmem   : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTemps_mean_gmem  : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTemps_min_stars  : fltarr(n_elements(rtr_stars.child_counts)) ,$
             maxTemps_mean_stars : fltarr(n_elements(rtr_stars.child_counts)) ,$
                                                                      $
             maxTempTime_min_gal    : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; time of min(maxTemps)
             maxTempTime_mean_gal   : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; mean(times maxTemps)
             maxTempTime_min_gmem   : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTempTime_mean_gmem  : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTempTime_min_stars  : fltarr(n_elements(rtr_stars.child_counts)) ,$
             maxTempTime_mean_stars : fltarr(n_elements(rtr_stars.child_counts))  }
        
        for i=0L,n_elements(rtr_gal.child_counts)-1 do begin
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
        
        for i=0L,n_elements(rtr_gmem.child_counts)-1 do begin
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
        
        offset = 0L
        
        for i=0L,n_elements(rtr_stars.child_counts)-1 do begin
          if rtr_stars.child_counts[i] gt 0 then begin
            ; for each star particle, collect its child tracers
            childInds = lindgen(rtr_stars.child_counts[i]) + offset
            
            locMaxTemps    = rtr_stars.maxTemps[childInds]
            locMaxTempTime = rtr_stars.maxTempTime[childInds]
            
            ; save the maximum value reached by any child tracer and the stddev of the population maxima
            r.maxTemps_stars[i]    = max(locMaxTemps,indmax)
            r.maxTempTime_stars[i] = locMaxTempTime[indmax]
            r.maxTempDisp_stars[i] = stddev(locMaxTemps)
            
            ; save minimum and means (and associated times)
            r.maxTemps_min_stars[i]     = min(locMaxTemps,indmin)
            r.maxTempTime_min_stars[i]  = locMaxTempTime[indmin]
            r.maxTemps_mean_stars[i]    = mean(locMaxTemps)
            r.maxTempTime_mean_stars[i] = mean(locMaxTempTime)
            
            offset += rtr_stars.child_counts[i]
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

    print,'NOTE: no stars for tracerVEL'

    if ~file_test(resFilename) then begin ; no restart   
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
    endif else begin
      ; restart
      if ~keyword_set(restart) then message,'Error: Restart file exists but restart not requested.'
      restore,resFilename,/verbose
      snapRange[0] = m
    endelse 
    
    for m=snapRange[0],snapRange[1],1 do begin
      sP.snap = m
      print,m  
  
      ; save restart?
      if m mod 10 eq 0 and m gt minSnap and keyword_set(restart) then begin
        print,' --- Writing restart! ---'
        save,rtr_gal,rtr_gmem,galcat_gal_trids,galcat_gmem_trids,m,filename=resFilename
        print,' --- Done! ---'
      endif
  
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
      
      temp_gal  = codeTempToLogK(tr_maxtemp[trids_gal_ind]) ; tracer output still in unit system
      temp_gmem = codeTempToLogK(tr_maxtemp[trids_gmem_ind])
      tr_maxtemp = !NULL
      
      tr_maxtemp_time = loadSnapshotSubset(sP=sP,partType='tracerVel',field='tracer_maxtemp_time')
      
      temp_time_gal  = tr_maxtemp_time[trids_gal_ind]
      temp_time_gmem = tr_maxtemp_time[trids_gmem_ind]
      tr_maxtemp_time = !NULL
      
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
      if sP.snap eq maxSnap then begin
        ; (1) full tracer information (galaxy members) - set savefilename
        saveFilename = sP.derivPath + 'maxtemp.trVel.gal.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(maxSnap)+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
        endif else begin
          ; output all gal/gmem tracer values at this point and the child counts so they can be 
          ; matched back up to their parents later
          save,rtr_gal,filename=saveFilename
          print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        endelse
        
        ; (2) full tracer information (group members) - set savefilename
        saveFilename = sP.derivPath + 'maxtemp.trVel.gmem.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(maxSnap)+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
        endif else begin
          ; output all gal/gmem tracer values at this point and the child counts so they can be 
          ; matched back up to their parents later
          save,rtr_gmem,filename=saveFilename
          print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        endelse
        
        ; (3) values condensed to gas parents - set savefilename
        saveFilename = sP.derivPath + 'maxtemp.trVel.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(maxSnap)+'.sav'

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

; maxTempsAll(): record all max temps (to save time, can start from snap>0 given an existing save)

pro maxTempsAll, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; set minimum snapshot (maxmimum redshift) and maximum snapshot (minimum redshift)
  minSnap = 0
  maxSnap = sP.snap - 1
  
  ; search for pre-existing maxTempAll save
  results = file_search(sP.derivPath+'maxTempAll.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.*.sav')
  
  if results[0] ne '' then begin
    ; use maximum previous save that is less than our maxSnap
    prevSaveSnap = -1
    
    foreach res,results do begin
      ; get saved snapshot number from filename
      sSplit = strsplit(res,'.',/extract)
      curSaveSnap = fix(sSplit[-2])
      help,curSaveSnap
      if curSaveSnap gt prevSaveSnap and curSaveSnap lt maxSnap then prevSaveSnap = curSaveSnap
    endforeach
    
    ; modify minSnap if we have a valid pre-existing save
    if prevSaveSnap ne -1 then minSnap = prevSaveSnap+1
  endif  
  
  snapRange = [minSnap,maxSnap]
  
  ; set save filename and check existence
  saveFilename = sP.derivPath + 'maxTempAll.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+str(maxSnap)+'.sav'

  if file_test(saveFilename) then message,'Error: Save file already exists.'

  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    message,'todo'
  endif
  
  ; MONTE CARLO TRACERS CASE - track all tracers (constant number) by ID
  ; ------------------------
  if sP.trMCPerCell gt 0 then begin
  
    ; for debugging:
    ;  ; load tracer IDs at minimum snap (this is really unneccessary, but gives us a debug check)
    ;  sP.snap = minSnap
    ;  tr_ids_orig = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
    ;   
    ;   ; sort and verify
    ;   tr_ids_orig_sort = sort(tr_ids_orig)
    ;   tr_ids_orig = tr_ids_orig[tr_ids_orig_sort] ; sorted ascending
    ;   nTracers = n_elements(tr_ids_orig)
      
    h = loadSnapshotHeader(sP=sP)
    nTracers = h.nPartTot[partTypeNum('tracerMC')]
    
    ; store the main arrays for all tracers as structures so we can write them directly
    if minSnap eq 0 then begin
      rtr_all  = { maxTemps      : fltarr(nTracers)  ,$
                   maxTempTime   : fltarr(nTracers)   }
    endif else begin
      ; starting with a previous save, just load rtr_all
      loadFilename = sP.derivPath + 'maxTempAll.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+str(prevSaveSnap)+'.sav'
      restore,loadFilename,/verbose
    endelse
    
    for m=snapRange[0],snapRange[1],1 do begin
      sP.snap = m
      print,m
      
      ; load tracer ids and match to original list
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
      
      ; for debugging:
      ;match,tr_ids,tr_ids_orig,ind_tr_ids,ind_tr_ids_orig,count=countOrig
      ;if countOrig ne n_elements(tr_ids_orig) or $
      ;   n_elements(tr_ids) ne n_elements(tr_ids_orig) or $
      ;   array_equal(sort(tr_ids),ind_tr_ids) ne 1 or $
      ;   array_equal(tr_ids[ind_tr_ids],tr_ids_orig) ne 1 then message,'Error: Bad tracer IDs.'
      ;ind_tr_ids_orig = !NULL
      
      ind_tr_ids = sort(tr_ids) ; quicker than match
      tr_ids = !NULL
      
      ; load tracer maximum temperature and sub-snapshot time at this snapshot
      tr_maxtemp = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxtemp')
      
      ; tracerMC maxtemp field in Kelvin
      w = where(tr_maxtemp eq 0,count)
      if count gt 0 then tr_maxtemp[w] = 1.0
      tr_maxtemp = alog10(tr_maxtemp)
      w = !NULL
      
      tr_maxtemp = tr_maxtemp[ind_tr_ids] ; ind_tr_ids same as sort(tr_ids)
      
      ; replace existing values if current snapshot has higher temps (galaxy members)
      w = where(tr_maxtemp gt rtr_all.maxTemps,count)
      if (count gt 0) then begin
        rtr_all.maxTemps[w] = tr_maxtemp[w]
        tr_maxtemp = !NULL
        
        tr_maxtemp_time = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxtemp_time')
        tr_maxtemp_time = tr_maxtemp_time[ind_tr_ids]
        
        rtr_all.maxTempTime[w] = tr_maxtemp_time[w] ; sub-snapshot timing
      endif
      
      ind_tr_ids = !NULL
      tr_maxtemp_time = !NULL
      
      ; SAVE?
      if sP.snap eq maxSnap then begin
        save,rtr_all,filename=saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
      endif ; save
      
    endfor
    
  endif

  ; VELOCITY TRACERS case - track all tracers (constant number) by ID
  ; ---------------------
  if sP.trMCPerCell eq -1 then begin
    message,'todo'
  endif  

end

; evolIGMTemp(): evolution of IGM gas (outside halos) with time

function evolIGMTemp, sP=sP
  
  ; config
  redshiftEnd = 2.0
  
  snapRange = [sP.groupCatRange[0],redshiftToSnapNum(redshiftEnd,sP=sP)] ;z=6-2
  nSnaps = snapRange[1]-snapRange[0]+1
  
  ; compare full distributions at these snapshots
  distRedshifts = [6.0,5.0,4.0,3.0,2.0]
  distSnaps     = redshiftToSnapNum(distRedshifts,sP=sP)
  
  tempBinSize = 0.05 / (sP.res/128)
  tempMinMax = [1.0,7.0]
  nTempBins = fix((tempMinMax[1]-tempMinMax[0]) / tempBinSize + 1)
  
  ; check if save exists
  saveFilename = sP.derivPath + 'igmtemp.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snapRange[0]) + '-' + str(sP.snapRange[1]) + '.sav'
  
  ; results exist, return
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif
    
  ; arrays
  r = { meanTemp      : fltarr(nSnaps)                           ,$
        medianTemp    : fltarr(nSnaps)                           ,$
        minTemp       : fltarr(nSnaps)                           ,$
        redshifts     : fltarr(nSnaps)                           ,$
        tempDist      : fltarr(n_elements(distSnaps),nTempBins)  ,$
        tempBinCen    : fltarr(nTempBins)                        ,$
        distRedshifts : distRedshifts                            ,$
        nTempBins     : nTempBins                                ,$
        tempMinMax    : tempMinMax                               ,$
        tempBinSize   : tempBinSize                               }
  
  distCount = 0
  
  for m=snapRange[0],snapRange[1] do begin

    ; load group catalog and gas ids
    sP.snap = m & print,m
    h = loadSnapshotHeader(sP=sP)
    gc = loadGroupCat(sP=sP,/readIDs)
    
    ; exclude all gas in groups from consideration
    gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    
    match,gas_ids,gc.IDs,ind1,ind2,count=count
    gas_igm_inds = bytarr(n_elements(gas_ids))
    if count gt 0 then gas_igm_inds[ind1] = 1B
    gas_igm_inds = where(gas_igm_inds eq 0B)
    gas_ids = !NULL
    
    ; load gas temperatures and store mean and median of IGM gas
    u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
    nelec = loadSnapshotSubset(sP=sP,partType='gas',field='ne')
    temp  = convertUtoTemp(u,nelec,/log)
    temp  = temp[gas_igm_inds]
    
    ; store
    r.meanTemp[m-snapRange[0]]   = mean(temp)
    r.medianTemp[m-snapRange[0]] = median(temp)
    r.minTemp[m-snapRange[0]]    = min(temp)
    r.redshifts[m-snapRange[0]]  = 1/h.time-1
    
    ; histogram distribution if requested
    if total(distSnaps eq m) gt 0 then begin
      hist = histogram(temp,min=tempMinMax[0],max=tempMinMax[1],binsize=tempBinSize,loc=loc)
      r.tempDist[distCount,*] = hist
      r.tempBinCen = loc + tempBinSize*0.5
      distCount += 1
    endif
    
  endfor
  
  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
   
  return,r
end

pro plotEvolIGMTemp

  sP_UV   = simParams(res=256,run='tracernew')
  sP_noUV = simParams(res=256,run='tracer_nouv')
  
  igmEvol_UV   = evolIGMTemp(sP=sP_UV)
  igmEvol_noUV = evolIGMTemp(sP=sP_noUV)
  
  ; plot (1) - mean/median evolution of IGM gas temperature with redshift
  start_PS, sP_UV.plotPath + 'igmTemp.comp.'+sP_UV.savPrefix+str(sP_UV.res)+'.eps'
    
    xrange = [6.0,2.0]
    yrange = [1.0,4.5]
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="IGM Gas Temperature [ log K ]",xtitle="Redshift"
    
    ; with UV background
    cgPlot,igmEvol_UV.redshifts,igmEvol_UV.meanTemp,color=getColor(3),line=1,/overplot
    cgPlot,igmEvol_UV.redshifts,igmEvol_UV.medianTemp,color=getColor(3),line=2,/overplot
    cgPlot,igmEvol_UV.redshifts,igmEvol_UV.minTemp,color=getColor(3),line=0,/overplot

    ; without UV background
    cgPlot,igmEvol_noUV.redshifts,igmEvol_noUV.meanTemp,color=getColor(1),line=1,/overplot
    cgPlot,igmEvol_noUV.redshifts,igmEvol_noUV.medianTemp,color=getColor(1),line=2,/overplot
    cgPlot,igmEvol_noUV.redshifts,igmEvol_noUV.minTemp,color=getColor(1),line=0,/overplot
    
    ; legend
    strings = ["w/ UV (mean)","w/ UV (median)","w/ UV (min)",$
               "no UV (mean)","no UV (median)","no UV (min)"]
    legend,strings,linestyle=[1,2,0,1,2,0],textcolors=getColor([3,3,3,1,1,1],/name),$
      box=0,linesize=0.5,charsize=!p.charsize-0.4,position=[5.9,3.3]
    
  end_PS
  
  ; plot (2) - distributions
  start_PS, sP_UV.plotPath + 'igmTemp.compdist.'+sP_UV.savPrefix+str(sP_UV.res)+'.eps'
  
  x0 = 0.15 & x1 = 0.55 & x2 = 0.95
  y0 = 0.15 & y1 = 0.55 & y2 = 0.95
  
  pos = list( [x0,y1,x1,y2] ,$ ; ul
              [x1,y1,x2,y2] ,$ ; ur
              [x0,y0,x1,y1] ,$ ; ll
              [x1,y0,x2,y1]  ) ; lr    

  inds = [0,2,3,4] ;z=6,4,3,2
  
  xrange = [1.1,6.6]
  yrange = [50,max(igmEvol_UV.tempDist)*2.0]    
  
  ; ul
  cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,$
    ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0];ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4,
  
  cgPlot,igmEvol_UV.tempBinCen,igmEvol_UV.tempDist[inds[0],*],line=0,color=getColor(3),/overplot
  cgPlot,igmEvol_noUV.tempBinCen,igmEvol_noUV.tempDist[inds[0],*],line=0,color=getColor(1),/overplot  
  cgText,xrange[1]*0.86,yrange[1]*0.2,"z="+string(igmEvol_UV.distRedshifts[inds[0]],format='(f3.1)'),alignment=0.5
  
  legend,["w/ UV","no UV"],textcolors=getColor([3,1],/name),box=0,linesize=0.25,position=[4.7,2e3]
  
  ; ur
  cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,$
    ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
  
  cgPlot,igmEvol_UV.tempBinCen,igmEvol_UV.tempDist[inds[1],*],line=0,color=getColor(3),/overplot
  cgPlot,igmEvol_noUV.tempBinCen,igmEvol_noUV.tempDist[inds[1],*],line=0,color=getColor(1),/overplot 
  cgText,xrange[1]*0.86,yrange[1]*0.2,"z="+string(igmEvol_UV.distRedshifts[inds[1]],format='(f3.1)'),alignment=0.5
  
  ; ll
  cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
    ytitle="",xtitle="",pos=pos[2],/noerase
  
  cgPlot,igmEvol_UV.tempBinCen,igmEvol_UV.tempDist[inds[2],*],line=0,color=getColor(3),/overplot
  cgPlot,igmEvol_noUV.tempBinCen,igmEvol_noUV.tempDist[inds[2],*],line=0,color=getColor(1),/overplot 
  cgText,xrange[1]*0.86,yrange[1]*0.2,"z="+string(igmEvol_UV.distRedshifts[inds[2]],format='(f3.1)'),alignment=0.5
  
  ; lr
  cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
    ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],/noerase;xticks=4,xtickv=[10.5,11.0,11.5,12.0,12.5]
  
  cgPlot,igmEvol_UV.tempBinCen,igmEvol_UV.tempDist[inds[3],*],line=0,color=getColor(3),/overplot
  cgPlot,igmEvol_noUV.tempBinCen,igmEvol_noUV.tempDist[inds[3],*],line=0,color=getColor(1),/overplot 
  cgText,xrange[1]*0.86,yrange[1]*0.2,"z="+string(igmEvol_UV.distRedshifts[inds[3]],format='(f3.1)'),alignment=0.5
   
  ; labels
  cgText,0.05,y1,"N",alignment=0.5,orientation=90.0,/normal
  cgText,x1,0.01,textoidl("IGM Gas Temperature [ log K ]"),alignment=0.5,/normal    
  
  end_PS

end

; checkStarIDs(): make sure SPH star ids turn once into gas IDs moving backwards in time
; note: in Arepo runs spawned (not converted) stars will have new IDs with no progenitor gas cell info

pro checkStarIDs

  sP = simParams(res=128,run='gadget',redshift=2.0)
  
  ; load all star particle IDs
  ids = loadSnapshotSubset(sP=sP,partType='stars',field='ids')
  ids_sort = sort(ids)
  mask = intarr(n_elements(ids))
   
  for m=sP.snap,0,-1 do begin
    sP.snap = m
    
    ; load gas ids and match
    gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    match,gas_ids,ids,gas_ind,star_ind,count=count1
    star_ind = star_ind[ids_sort]
    
    if count1 gt 0 then mask[star_ind] += 1
    
    ; load star ids and match
    star_ids = loadSnapshotSubset(sP=sP,partType='stars',field='ids')
    match,star_ids,ids,star_ind_cur,star_ind_orig,count=count2
    
    if count1+count2 ne n_elements(ids) then message,'Did not find all original star IDs.'
    
    ; for those stars that are still stars, make sure we have never seen them as gas
    star_ind_orig = star_ind_orig[ids_sort]
    if max(mask[star_ind_orig]) gt 0 then message,'Error: Flip gas back to star.'
    
    print,m,float(count2)/(count1+count2),float(count1)/(count1+count2)
  endfor

end
