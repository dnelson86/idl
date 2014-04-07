; maxVals.pro
; gas accretion project - past temperature/entropy/density history of gas
; dnelson oct.2013

; -----------------------------------------------------------------------------------------------------
; maxVals(): find maximum temperature/ent/dens for gas particles in galaxyCat from the start of the
;   simulation up to sP.redshift (maxTemp only while not on eEOS)
;
; NOTE: vals are only saved for gas in the galaxyCat at the end of the interval (not all gas)
; NOTE: vals are only recorded until gas crosses 0.15rvir of main progenitor branch for the first time
; -----------------------------------------------------------------------------------------------------

function maxVals, sP=sP, restart=restart

  forward_function cosmoTracerChildren, cosmoTracerVelParents
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; set minimum snapshot (maxmimum redshift)
  zStart = snapNumToRedshift(snap=0,sP=sP)
  minSnap = redshiftToSnapnum(zStart,sP=sP)

  ; set maximum snapshot (minimum redshift) (-1 do not include final snapshot in Tmax search)
  maxSnap = sP.snap - 1

  snapStep = 1 ; process every snapshot
  ;if sP.snapRange[1] ge 300 then snapStep = 2 ; every other for run=tracer
  maxSnap = maxSnap - (maxSnap mod snapStep) ; make sure we save at the end
  
  snapRange = [minSnap,maxSnap]
  
  ; set saveFilename and check for existence
  saveFilename = sP.derivPath+'maxVals.' + sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(minSnap)+'-'+str(maxSnap)+'.sav'

  if file_test(saveFilename) then begin
    restore, saveFilename
    return, rtr
  endif
  
  ; load first 0.15rvir crossing times
  at = accretionTimes(sP=sP)
  if (size(at.accTime))[1] ne n_elements(at.rVirFacs)+2 then message,'Error: Old at, update.'

  at = reform( at.AccTime[-2,*] ) ; last two are special (first) 0.15 and 1.0 crossings
  w = where(at lt 0.0,count)
  if count gt 0 then at[w] = 1.0 ; do not restrict maxvals recording if no recorded 0.15 accTime
  
  ; set filename for restart file
  resFilename = sP.derivPath + 'maxVals.' + sP.saveTag+'.restart.'+sP.savPrefix+str(sP.res)+'.'+$
                        str(minSnap)+'-'+str(maxSnap)+'.sav'
                        
  ; check for maxTempsAll save (in this case we don't have to loop over any snapshots, but load directly)
  maxTempsAllFlag = 0
  maxTempsAllSaveFilename = sP.derivPath + 'maxVals.' + 'All.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+str(maxSnap)+'.sav'

  if file_test(maxTempsAllSaveFilename) then maxTempsAllFlag = 1
  
  ; load galaxy/group member catalogs at zMin for gas ids to search for
  galcat = galaxyCat(sP=sP) ;sP.snap is still at zMin
  
  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
  
    print,'Calculating new maxVals using ( SPH Particles ) res = '+str(sP.res)+$
      ' in range ['+str(minSnap)+'-'+str(maxSnap)+'].'
      
    if ~file_test(resFilename) then begin ; no restart  
        accMask = bytarr(galcat.countTot)
        
        ; store the main arrays for all tracers as structures so we can write them directly
        rtr = {maxTemps      : fltarr(galcat.countTot)   ,$
               maxTempTime   : fltarr(galcat.countTot)   ,$
               maxEnt        : fltarr(galcat.countTot)   ,$
               maxDens       : fltarr(galcat.countTot)   ,$
               maxMachNum    : 0 ,$ ; no machNum for sph
               maxTempDisp   : 0 ,$
               maxTemps_min  : 0 ,$
               maxTemps_mean : 0  } ; no disp,min,mean values for sph
             
    endif else begin
      ; restart
      if ~keyword_set(restart) then message,'Error: Restart file exists but restart not requested.'
      restore,resFilename,/verbose
      snapRange[0] = m
    endelse    
    
    for m=snapRange[0],snapRange[1],snapStep do begin
      sP.snap = m
      
      ; save restart?
      if m mod 20 eq 0 and m gt snapRange[0] and keyword_set(restart) then begin
        print,' --- Writing restart! ---'
        save,rtr,m,accMask,filename=resFilename
        print,' --- Done! ---'
        ;exit,status=33 ; requeue
      endif
      
      ; load gas ids and match to catalog
      ids_gas = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs, need this if we want ids[ids_ind], 
      ; temp[ids_ind], etc to be in the same order as the group catalog id list
      calcMatch,galcat.ids,ids_gas,galcat_ind,ids_ind,count=countMatch
      ids_ind = ids_ind[sort(galcat_ind)]
      galcat_ind = galcat_ind[sort(galcat_ind)]
      
      ids_gas = !NULL
      
      ; update mask for any tracers passing 0.15rvir point (stop recording maxvals)
      h = loadSnapshotHeader(sP=sP)
      w = where(at lt h.time,countAccMask)
      if countAccMask gt 0 then accMask[w] = 1B
      if sP.gfmWinds then message,'A bit strange.' ; update mask if windcounter>0
      
      ; temp
      ; ----
      u = loadSnapshotSubset(sP=sP,partType='gas',field='u')
      u = u[ids_ind]
        
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='ne')
      nelec = nelec[ids_ind]
       
      temp = convertUtoTemp(u,nelec,/log)
      
      nelec = !NULL
      
      ; entropy, density
      ; ----------------
      dens = loadSnapshotSubset(sP=sP,partType='gas',field='dens')
      dens = dens[ids_ind]
        
      ent = calcEntropyCGS(u,dens,sP=sP,/log)
      
      u = !NULL
      
      ; load gas SFR and select off the effective EOS (SFR must be zero)
      sfr = loadSnapshotSubset(sP=sP,partType='gas',field='sfr')
      sfr = sfr[ids_ind]
      
      ; replace existing values if current snapshot has higher temps (enforce off effective EOS)
      w1 = where(temp gt rtr.maxTemps[galcat_ind] and sfr eq 0.0 and accMask eq 0B,countTemp)
      if countTemp gt 0 then begin
        rtr.maxTemps[galcat_ind[w1]]    = temp[w1]
        rtr.maxTempTime[galcat_ind[w1]] = snapNumToRedshift(sP=sP,/time)
      endif
            
      w1 = where(ent gt rtr.maxEnt[galcat_ind] and accMask eq 0B,countEnt)
      if countEnt gt 0 then rtr.maxEnt[galcat_ind[w1]] = ent[w1]
      
      w1 = where(dens gt rtr.maxDens[galcat_ind] and accMask eq 0B,countDens)
      if countDens gt 0 then rtr.maxDens[galcat_ind[w1]] = dens[w1]
      
      ; output
      fracInMask = string( float(countAccMask)/galcat.countTot * 100, format='(f5.1)')
      fracInWind = string( 0.0,                                       format='(f5.1)')
      fracTemp   = string( float(countTemp)/galcat.countTot * 100,    format='(f5.1)')
      fracEnt    = string( float(countEnt)/galcat.countTot * 100,     format='(f5.1)')
      fracDens   = string( float(countDens)/galcat.countTot * 100,    format='(f5.1)')
      fracMach   = string( 0.0,                                       format='(f5.1)')
        
      print,'['+string(m,format='(I3.3)')+'] z='+string(snapNumToRedshift(sP=sP),format='(f5.1)')+$
        ' fracInMask: '+fracInMask+' fracInWind: '+fracInWind+' fracTemp: '+fracTemp+$
        ' fracEnt: '+fracEnt+' fracDens: '+fracDens+' fracMach: '+fracMach
      
      ; SAVE?
      if sP.snap eq maxSnap then begin
        ; set savefilename
        saveFilename = sP.derivPath+'maxVals.' + sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(maxSnap)+'.sav'
      
        save,rtr,filename=saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
      endif ; save
    endfor ;m
  endif
  
  ; MONTE CARLO TRACERS CASE - for each original gas cell, determine some statistics of its
  ; population of tracers and an estimate for the dispersion in those statistics
  ; ------------------------
  if sP.trMCPerCell gt 0 then begin
  
    print,'Calculating new maxVals using ( TracerMC ) res = '+str(sP.res)+$
      ' in range ['+str(minSnap)+'-'+str(maxSnap)+'].'
    if maxTempsAllFlag eq 1 then print,'Using ['+maxTempsAllSaveFilename+']'
      
    if ~file_test(resFilename) then begin ; no restart   
      
      galcat_trids   = galcat.trMC_ids
      num_galcat_trs = n_elements(galcat_trids)
      accMask = bytarr(num_galcat_trs)
  
      ; store the main arrays for all tracers as structures so we can write them directly
      rtr  = { maxTemps      : fltarr(num_galcat_trs)  ,$
               maxTempTime   : fltarr(num_galcat_trs)  ,$
               maxEnt        : fltarr(num_galcat_trs)  ,$
               maxDens       : fltarr(num_galcat_trs)  ,$
               maxMachNum    : fltarr(num_galcat_trs)   }
               
    endif else begin
      ; restart
      if ~keyword_set(restart) then message,'Error: Restart file exists but restart not requested.'
      restore,resFilename,/verbose
      snapRange[0] = m
    endelse 

    for m=snapRange[0],snapRange[1],snapStep do begin
      sP.snap = m
  
      ; save restart?
      if m mod 20 eq 0 and m gt snapRange[0] and keyword_set(restart) then begin
        print,' --- Writing restart! ---'
        save,rtr,galcat_trids,num_galcat_trs,m,accMask,filename=resFilename
        print,' --- Done! ---'
        ;exit,status=33 ; requeue
      endif
      
      if maxTempsAllFlag eq 0 then begin ; normal
  
        ; load tracer ids and match to child ids from zMin
        tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
       
        idIndexMap = getIDIndexMap(tr_ids,minid=minid)
          
        trids_ind = idIndexMap[galcat_trids-minid]
        
        idIndexMap = !NULL
        tr_ids     = !NULL
        galcat_ind = !NULL
        
        ; update mask for any tracers passing 0.15rvir point (stop recording maxvals)
        h = loadSnapshotHeader(sP=sP)
        w = where(at lt h.time,countAccMask)
        if countAccMask gt 0 then accMask[w] = 1B
        
        ; update mask if windcounter>0
        ; NO
        countAccWind = 0
        ;if sP.gfmWinds then begin
        ;  windc = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_windcounter')
        ;  windc = windc[trids_ind]
        ;  w = where(windc gt 0,countAccWind)
        ;  if countAccWind gt 0 then accMask[w] = 1B
        ;endif
        
        ; maxtemp
        ; -------
        tr_maxval = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxtemp')
        
        ; tracerMC maxtemp field in Kelvin, tracerVEL still in unit system, convert to log
        tr_maxval = tr_maxval[trids_ind]   
        tr_maxval = mylog10(tr_maxval)

        ; replace existing values if current snapshot has higher temps
        ; note: if tracers are still inside a star particle, their maxtemp entry will be zero
        w1 = where(tr_maxval gt rtr.maxTemps and accMask eq 0B,countTemp)
        
        if countTemp gt 0 then begin
          rtr.maxTemps[w1] = tr_maxval[w1]
          
          ; sub-snapshot timing?
          if sP.trMCFields[7] ge 0 then begin
            tr_maxval_time = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxtemp_time')
            rtr.maxTempTime[w1] = tr_maxval_time[trids_ind[w1]]
            tr_maxval_time = !NULL
          endif else begin
            ; if not, just use constant time corresponding to this snapshot
            rtr.maxTempTime[w1] = snapNumToRedshift(sP=sP,/time)
          endelse
        endif
        
        ; maxent
        ; ------
        tr_maxval = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxent')
        
        ; convert entropy to log(cgs)
        tr_maxval = tr_maxval[trids_ind]  
        tr_maxval = convertTracerEntToCGS(tr_maxval,/log,sP=sP)
        
        w1 = where(tr_maxval gt rtr.maxEnt and accMask eq 0B,countEnt)
        if countEnt gt 0 then rtr.maxEnt[w1] = tr_maxval[w1]
        
        ; maxdens
        ; -------
        tr_maxval = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxdens')
        tr_maxval = tr_maxval[trids_ind]  
        
        w1 = where(tr_maxval gt rtr.maxDens and accMask eq 0B,countDens)
        if countDens gt 0 then rtr.maxDens[w1] = tr_maxval[w1]
        
        ; maxmachnum
        ; ----------
        tr_maxval = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxmachnum')
        
        w1 = where(tr_maxval gt rtr.maxMachNum and accMask eq 0B,countMach)
        if countMach gt 0 then rtr.maxMachNum[w1] = tr_maxval[w1]
        
        tr_maxval      = !NULL
        tr_maxval_time = !NULL
        
        ; output
        fracInMask = string( float(countAccMask)/num_galcat_trs * 100, format='(f5.1)')
        fracInWind = string( float(countAccWind)/num_galcat_trs * 100, format='(f5.1)')
        fracTemp   = string( float(countTemp)/num_galcat_trs * 100,    format='(f5.1)')
        fracEnt    = string( float(countEnt)/num_galcat_trs * 100,     format='(f5.1)')
        fracDens   = string( float(countDens)/num_galcat_trs * 100,    format='(f5.1)')
        fracMach   = string( float(countMach)/num_galcat_trs * 100,    format='(f5.1)')
        
        print,'['+string(m,format='(I3.3)')+'] z='+string(snapNumToRedshift(sP=sP),format='(f5.1)')+$
          ' fracInMask: '+fracInMask+' fracInWind: '+fracInWind+' fracTemp: '+fracTemp+$
          ' fracEnt: '+fracEnt+' fracDens: '+fracDens+' fracMach: '+fracMach
        
      endif else begin ; we have a maxTempsAll save, just load and move immediately to save in correct format
        
        ; load tracer ids and match to maxTempsAll save order (sorted ascending)
        tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
        tr_ids = tr_ids[sort(tr_ids)]

        idIndexMap = getIDIndexMap(tr_ids,minid=minid)
          
        trids_ind  = idIndexMap[galcat_trids-minid]
        
        idIndexMap = !NULL
        tr_ids     = !NULL
        galcat_ind = !NULL
        
        ; fill rtr
        restore,maxTempsAllSaveFilename,/verbose
        
        if n_tags(rtr) ne n_tags(rtr_all) then message,'Error: Structure mismatch.'
        
        for i=0,n_tags(rtr)-1 do rtr.(i) = (rtr_all.(i))[ trids_ind ]
        
        rtr_all = !NULL
        
        ; modify current snapshot
        sP.snap = maxSnap
        m = maxSnap
        
      endelse
      
      ; SAVE?
      if sP.snap eq maxSnap then begin
        ; (1) full tracer information (galaxy members) - set savefilename
        saveFilename = sP.derivPath + 'maxVals.' + sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(maxSnap)+'.sav'

        ; output all tracer values at this point
        save,rtr,filename=saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        
        ; NOTE: values condensed to gas parents (unused and removed for now, see r99)
      endif ;save
      
    endfor ;m  
    
  endif
  
  ; VELOCITY TRACERS case - will be similar to above since there could be multiple
  ; ---------------------
  if sP.trMCPerCell eq -1 then begin
  
    print,'Calculating new maxVals using ( TracerVEL ) res = '+str(sP.res)+$
      ' in range ['+str(minSnap)+'-'+str(maxSnap)+'].'

    if ~file_test(resFilename) then begin ; no restart   
    
      ; list of velocity tracers
      galcat_trids   = galcat.trVel_ids ;tr_ids[galcat_trids]
      galcat_num_trs = n_elements(galcat_trids)
      accMask        = bytarr(galcat_num_trs)
  
      ; store the main arrays for all tracers as structures so we can write them directly
      rtr  = { maxTemps      : fltarr(galcat_num_trs)  ,$
               maxTempTime   : fltarr(galcat_num_trs)  ,$
               maxEnt        : fltarr(galcat_num_trs)  ,$
               maxDens       : fltarr(galcat_num_trs)  ,$
               maxMachNum    : fltarr(galcat_num_trs)   }
               
    endif else begin
      ; restart
      if ~keyword_set(restart) then message,'Error: Restart file exists but restart not requested.'
      restore,resFilename,/verbose
      snapRange[0] = m
    endelse 
    
    for m=snapRange[0],snapRange[1],snapStep do begin
      sP.snap = m
  
      ; save restart?
      if m mod 20 eq 0 and m gt snapRange[0] and keyword_set(restart) then begin
        print,' --- Writing restart! ---'
        save,rtr,galcat_trids,galcat_num_trs,m,accMask,filename=resFilename
        print,' --- Done! ---'
        ;exit,status=33 ; requeue
      endif
  
      ; load tracer ids and match to child ids from zMin
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerVel',field='ids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs, need this if we want ids[ids_ind], 
      ; temp[ids_ind], etc to be in the same order as the group catalog id list    
      calcMatch,galcat_trids,tr_ids,galcat_ind,trids_ind,count=countMatch
      trids_ind = trids_ind[sort(galcat_ind)]
      
      tr_ids     = !NULL
      galcat_ind = !NULL
      
      ; update mask for any tracers passing 0.15rvir point (stop recording maxvals)
      h = loadSnapshotHeader(sP=sP)
      w = where(at lt h.time,count)
      if count gt 0 then accMask[w] = 1B
      
      ; maxtemp
      ; -------
      tr_maxval = loadSnapshotSubset(sP=sP,partType='tracerVel',field='tracer_maxtemp')
      tr_maxval_time = loadSnapshotSubset(sP=sP,partType='tracerVel',field='tracer_maxtemp_time')
      
      tr_maxval = tr_maxval[trids_ind]
      tr_maxval_time = tr_maxval_time[trids_ind]
      
      tr_maxval = codeTempToLogK(tr_maxval) ; tracerVel output still in unit system
      
      ; replace existing values if current snapshot has higher temps (galaxy members)
      w1 = where(tr_maxval gt rtr.maxTemps and accMask eq 0B,countTemp)
      if countTemp gt 0 then begin
        rtr.maxTemps[w1]    = tr_maxval[w1]
        rtr.maxTempTime[w1] = tr_maxval_time[w1] ; sub-snapshot timing
      endif
      
      ; entropy
      ; -------
      tr_maxval = loadSnapshotSubset(sP=sP,partType='tracerVel',field='tracer_maxent')
      
      ; convert entropy to log(cgs)
      tr_maxval = tr_maxval[trids_ind]  
      tr_maxval = convertTracerEntToCGS(tr_maxval,/log,sP=sP)
      
      w1 = where(tr_maxval gt rtr.maxEnt and accMask eq 0B,countEnt)
      if countEnt gt 0 then rtr.maxEnt[w1] = tr_maxval[w1]
        
      ; density
      ; -------
      tr_maxval = loadSnapshotSubset(sP=sP,partType='tracerVel',field='tracer_maxdens')
      tr_maxval = tr_maxval[trids_ind]
        
      w1 = where(tr_maxval gt rtr.maxDens and accMask eq 0B,countDens)
      if countDens gt 0 then rtr.maxDens[w1] = tr_maxval[w1]
        
      ; maxmachnum
      ; ----------
      tr_maxval = loadSnapshotSubset(sP=sP,partType='tracerVel',field='tracer_maxmachnum')
        
      w1 = where(tr_maxval gt rtr.maxMachNum and accMask eq 0B,countMach)
      if countMach gt 0 then rtr.maxMachNum[w1] = tr_maxval[w1]
        
      tr_maxval      = !NULL
      tr_maxval_time = !NULL
      
      ; output
      fracInMask = string( float(countAccMask)/galcat_num_trs * 100, format='(f5.1)')
      fracInWind = string( 0.0,                                       format='(f5.1)')
      fracTemp   = string( float(countTemp)/galcat_num_trs * 100,    format='(f5.1)')
      fracEnt    = string( float(countEnt)/galcat_num_trs * 100,     format='(f5.1)')
      fracDens   = string( float(countDens)/galcat_num_trs * 100,    format='(f5.1)')
      fracMach   = string( float(countMach)/galcat_num_trs * 100,    format='(f5.1)')
        
      print,'['+string(m,format='(I3.3)')+'] z='+string(snapNumToRedshift(sP=sP),format='(f5.1)')+$
        ' fracInMask: '+fracInMask+' fracInWind: '+fracInWind+' fracTemp: '+fracTemp+$
        ' fracEnt: '+fracEnt+' fracDens: '+fracDens+' fracMach: '+fracMach
      
      ; SAVE?
      if sP.snap eq maxSnap then begin
        ; (1) full tracer information (galaxy members) - set savefilename
        saveFilename = sP.derivPath + 'maxVals.' + sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(maxSnap)+'.sav'

        ; output all values at this point
        save,rtr,filename=saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        
        ; NOTE: values condensed to gas parents (unused and removed for now, see r99)
      endif ;save
    endfor ;m
  endif
  
  sP.snap = maxSnap + 1 ; restore sP.snap
  return, rtr

end

; maxValsAll(): record all max vals of all tracers (to save time, can start from snap>0 given an existing save)

pro maxValsAll, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; set minimum snapshot (maxmimum redshift) and maximum snapshot (minimum redshift)
  minSnap = 0
  maxSnap = sP.snap - 1
  
  ; search for pre-existing maxTempAll save
  results = file_search(sP.derivPath+'maxVals.All.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.*.sav')
  
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
  saveFilename = sP.derivPath + 'maxVals.All.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+str(maxSnap)+'.sav'

  if file_test(saveFilename) then message,'Error: Save file already exists.'

  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    message,'todo'
  endif
  
  ; MONTE CARLO TRACERS CASE - track all tracers (constant number) by ID
  ; ------------------------
  if sP.trMCPerCell gt 0 then begin
      
    h = loadSnapshotHeader(sP=sP)
    nTracers = h.nPartTot[partTypeNum('tracerMC')]
    
    ; store the main arrays for all tracers as structures so we can write them directly
    if minSnap eq 0 then begin
      accMask = bytarr(nTracers)
      
      rtr_all  = { maxTemps      : fltarr(nTracers)  ,$
                   maxTempTime   : fltarr(nTracers)  ,$
                   maxEnt        : fltarr(nTracers)  ,$
                   maxDens       : fltarr(nTracers)  ,$
                   maxMachNum    : fltarr(nTracers)   }
    endif else begin
      ; starting with a previous save, just load rtr_all
      loadFilename = sP.derivPath + 'maxVals.All.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+str(prevSaveSnap)+'.sav'
      restore,loadFilename,/verbose
    endelse
    
    for m=snapRange[0],snapRange[1],snapStep do begin
      sP.snap = m
      sP.redshift = snapNumToRedshift(sP=sP)
      
      ; load tracer ids and match to original list
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
      trids_ind = sort(tr_ids) ; just put them in order
      tr_ids = !NULL
      
      ; update mask for any tracers passing 0.15rvir point (stop recording maxvals)
      h = loadSnapshotHeader(sP=sP)
      w = where(at lt h.time,countAccMask)
      if countAccMask gt 0 then accMask[w] = 1B
        
      ; update mask if windcounter>0
      if sP.gfmWinds then begin
        windc = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_windcounter')
        windc = windc[trids_ind]
        w = where(windc gt 0,countAccWind)
        if countAccWind gt 0 then accMask[w] = 1B
      endif
      
      ; maxtemp
      ; -------
      tr_maxval = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxtemp')
        
      ; tracerMC maxtemp field in Kelvin, tracerVEL still in unit system, convert to log
      tr_maxval = tr_maxval[trids_ind]   
      tr_maxval = mylog10(tr_maxval)
        
      ; replace existing values if current snapshot has higher temps
      ; note: if tracers are inside a star particle, their maxtemp entry will be zero
      w1 = where(tr_maxval gt rtr_all.maxTemps and accMask eq 0B,countTemp)
        
      if countTemp gt 0 then begin
        rtr_all.maxTemps[w1] = tr_maxval[w1]
        
        ; sub-snapshot timing?
        if sP.trMCFields[7] ge 0 then begin
          tr_maxval_time = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxtemp_time')
          rtr_all.maxTempTime[w1] = tr_maxval_time[trids_ind]
          tr_maxval_time = !NULL
        endif else begin
          ; if not, just use constant time corresponding to this snapshot
          rtr_all.maxTempTime[w1] = snapNumToRedshift(sP=sP,/time)
        endelse
      endif
        
      ; maxent
      ; ------
      tr_maxval = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxent')
      tr_maxval = tr_maxval[trids_ind]
       
      ; convert entropy to log(cgs)
      tr_maxval = convertTracerEntToCGS(tr_maxval,/log,sP=sP)
        
      w1 = where(tr_maxval gt rtr_all.maxEnt and accMask eq 0B,countEnt)
      if countEnt gt 0 then rtr_all.maxEnt[w1] = tr_maxval[w1]
        
      ; maxdens
      ; -------
      tr_maxval = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxdens')
      tr_maxval = tr_maxval[trids_ind]
        
      w1 = where(tr_maxval gt rtr_all.maxDens and accMask eq 0B,countDens)
      if countDens gt 0 then rtr_all.maxDens[w1] = tr_maxval[w1]
        
      ; maxmachnum
      ; ----------
      tr_maxval = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxmachnum')
      tr_maxval = tr_maxval[trids_ind]
       
      w1 = where(tr_maxval gt rtr_all.maxMachNum and accMask eq 0B,countMach)
      if countMach gt 0 then rtr_all.maxMachNum[w1] = tr_maxval[w1]
        
      tr_maxval      = !NULL
      tr_maxval_time = !NULL
            
      trids_ind = !NULL
      
      ; output
      fracInMask = string( float(countAccMask)/nTracers * 100, format='(f5.1)')
      fracInWind = string( float(countAccWind)/nTracers * 100, format='(f5.1)')
      fracTemp   = string( float(countTemp)/nTracers * 100,    format='(f5.1)')
      fracEnt    = string( float(countEnt)/nTracers * 100,     format='(f5.1)')
      fracDens   = string( float(countDens)/nTracers * 100,    format='(f5.1)')
      fracMach   = string( float(countMach)/nTracers * 100,    format='(f5.1)')
        
      print,'['+string(m,format='(I3.3)')+'] z='+string(sP.redshift,format='(f5.1)')+$
        ' fracInMask: '+fracInMask+' fracInWind: '+fracInWind+' fracTemp: '+fracTemp+$
        ' fracEnt: '+fracEnt+' fracDens: '+fracDens+' fracMach: '+fracMach
      
      ; SAVE?
      if sP.snap eq maxSnap then begin
        save,rtr_all,accMask,filename=saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
      endif ; save
      
    endfor
    
  endif

  ; VELOCITY TRACERS case - track all tracers (constant number) by ID
  ; ---------------------
  if sP.trMCPerCell eq -1 then begin
    message,'todo'
  endif
  
  sP.snap = maxSnap + 1 ; restore sP.snap

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
