; cosmoAnalysis.pro
; gas accretion project - main
; dnelson jan.2012

; galaxyCat(): if snap not specified, create and save complete galaxy catalog from the group catalog by 
;              imposing additional cut in the (rho,temp) plane (same as that used by Torrey+ 2011)
;              if snap is specified, create only for one snapshot number or return previously saved
;              results for that snapshot
; Note: the gal/gmem catalogs have the same size (and indexing) as the subgroups

function galaxyCat, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  galcut_T   = 6.0
  galcut_rho = 0.25

  ; if snap specified, run only one snapshot (and/or just return previous results)
  if (sP.snap ne -1) then begin
    saveFilename1 = sP.derivPath + 'galcat.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) +'.sav'
    saveFilename2 = sP.derivPath + 'groupmemcat.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.sav'
    
    ; results exist, return
    if (file_test(saveFilename1) and file_test(saveFilename2)) then begin
      restore,saveFilename1
      restore,saveFilename2
      r = {galaxyOff:galaxyOff,galaxyLen:galaxyLen,galaxyIDs:galaxyIDs,$
           groupmemOff:groupmemOff,groupmemLen:groupmemLen,groupmemIDs:groupmemIDs}
      return,r
    endif
    
    ; need to compute, set restricted range of snapshots to process
    snapRange = [sP.snap,sP.snap]
  endif else begin
    ;snapRange = [210,219] ;running 512 gadget
    snapRange = sP.groupCatRange
  endelse
  
  for m=snapRange[0],snapRange[1],1 do begin
    sP.snap = m
    ; skip if previous results exist
    saveFilename1 = sP.derivPath + 'galcat.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.sav'
    saveFilename2 = sP.derivPath + 'groupmemcat.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.sav'
    
    if (file_test(saveFilename1) and file_test(saveFilename2)) then begin
      print,'Skipping: '+strmid(saveFilename1,strlen(sP.derivPath))
      continue
    endif
    
    ; load ids of particles in all subfind groups
    gc = loadGroupCat(sP=sP,/readIDs)
    gcPIDs = gcPIDList(gc=gc,select='all',partType='gas')
    
    ; load gas ids and match to catalog
    ids  = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    match,gcPIDs,ids,gc_ind,ids_ind,count=countMatch,/sort
    ids_ind = ids_ind[sort(gc_ind)] ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs
                                    ; need this if we want ids[ids_ind], temp[ids_ind], etc to be
                                    ; in the same order as the group catalog id list

    gcPIDs = !NULL
    gc_ind = !NULL
    ids = ids[ids_ind]

    ; load u,nelec and calculate temp of gas
    u = loadSnapshotSubset(sP=sP,partType='gas',field='u')
    u = u[ids_ind]
    nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
    nelec = nelec[ids_ind]
    
    temp = convertUtoTemp(u,nelec)
    
    u     = !NULL
    nelec = !NULL
    
    ; load rho of gas and make galaxy (rho,temp) plane cut
    dens = loadSnapshotSubset(sP=sP,partType='gas',field='density')
    dens = dens[ids_ind]
    
    ; scale Torrey+ (2011) galaxy cut to physical density
    scalefac = snapNumToRedshift(sP=sP,/time) ; time flag gives simulation time = scale factor
    a3inv = 1.0 / (scalefac*scalefac*scalefac)
    dens *= a3inv
    
    w = where(alog10(temp) - galcut_rho * alog10(dens) lt galcut_T,$
              countCut,comp=wComp,ncomp=countComp)

    if (countCut eq 0 or countComp eq 0) then begin
      print,'Warning: Empty galaxy cut or comp. Skipping: ' + strmid(saveFilename1,strlen(sP.derivPath))
      continue
    endif
    
    temp = !NULL
    dens = !NULL
    
    ; make subsets of ids matching galaxy cut and complement
    ids_groupmem = ids[wComp]
    ids = ids[w]    

    ; construct group member catalog
    if (not file_test(saveFilename2)) then begin
    
      groupmemLen = ulonarr(gc.nSubgroupsTot)
      groupmemOff = ulonarr(gc.nSubgroupsTot)
      groupmemIDs = ulonarr(n_elements(ids_groupmem))
      
      nextOff = 0UL
      
      ; match indices between gas ids and group member ids
      match,ids_groupmem,gc.IDs,ids_ind,gcIDs_ind,count=countID,/sort
      
      for gcID=0,gc.nSubgroupsTot-1 do begin
        ; select ids in group
        groupStart = gc.subGroupOffset[gcID]
        groupEnd   = groupStart + gc.subGroupLen[gcID]
        
        w = where(gcIDs_ind ge groupStart and gcIDs_ind lt groupEnd,count)
        
        ; save in similar 1D offset format
        if (count gt 0) then begin
          groupmemLen[gcID] = count
          groupmemOff[gcID] = nextOff
          nextOff = groupmemOff[gcID] + groupmemLen[gcID]
    
          groupmemIDs[groupmemOff[gcID]:nextOff-1] = ids_groupmem[ids_ind[w]]
        endif      
        
      end
  
      ; debug: make sure all gas particles were found in the group catalog
      match,ids_groupmem,groupmemIDs,ind1,ind2,count=countCheck
      if (countCheck ne n_elements(ids_groupmem)) then begin
        print,'Uhoh, check2.',countCheck,n_elements(ids_groupmem)
        nuniq = n_elements(uniq(ids_groupmem[sort(ids_groupmem)]))
        stop
      endif

      ; save group membership catalog
      save,groupmemLen,groupmemOff,groupmemIDs,filename=saveFilename2    
      print,'Saved: '+strmid(saveFilename2,strlen(sP.derivPath))+' ['+str(countComp)+'/'+str(countMatch)+']'    
    
    endif

    ; construct galaxy catalog
    if (not file_test(saveFilename1)) then begin
      galaxyLen = ulonarr(gc.nSubgroupsTot)
      galaxyOff = ulonarr(gc.nSubgroupsTot)
      galaxyIDs = ulonarr(n_elements(ids))
      
      nextOff = 0UL
      
      ; match indices between galaxy ids and group ids
      match,ids,gc.IDs,ids_ind,gcIDs_ind,count=countID,/sort
      
      for gcID=0,gc.nSubgroupsTot-1 do begin
        ; select ids in group
        groupStart = gc.subGroupOffset[gcID]
        groupEnd   = groupStart + gc.subGroupLen[gcID]
        
        w = where(gcIDs_ind ge groupStart and gcIDs_ind lt groupEnd,count)
        
        ; save in similar 1D offset format
        if (count gt 0) then begin
          galaxyLen[gcID] = count
          galaxyOff[gcID] = nextOff
          nextOff = galaxyOff[gcID] + galaxyLen[gcID]
    
          galaxyIDs[galaxyOff[gcID]:nextOff-1] = ids[ids_ind[w]]
        endif      
        
      end
      
      ; save galaxy catalog
      save,galaxyLen,galaxyOff,galaxyIDs,filename=saveFilename1
      print,'Saved: '+strmid(saveFilename1,strlen(sP.derivPath))+' ['+str(countCut)+'/'+str(countMatch)+']'
    endif

  endfor ;snapRange
  
end

; galaxyCatRadii(): find radial distance of all group member particles wrt the group they belong to
;                   as well as the rad to the primary group if this is a secondary group

function galaxyCatRadii, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; if snap specified, run only one snapshot (and/or just return previous results)
  if (sP.snap ne -1) then begin
    saveFilename = sP.derivPath + 'galradii.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.sav'
    
    ; results exist, return
    if (file_test(saveFilename)) then begin
      restore,saveFilename
      r = {gal_pri:gal_pri,gal_sec:gal_sec,gmem_pri:gmem_pri,gmem_sec:gmem_sec}
      return,r
    endif
    
    ; need to compute, set restricted range of snapshots to process
    snapRange = [sP.snap,sP.snap]
  endif else begin
    snapRange = sP.groupCatRange
  endelse
  
  ; loop from target redshift to beginning of group catalogs
  for m=snapRange[0],snapRange[1],1 do begin
  
    ; set save filename and check for existence
    sP.snap  = m
    saveFilename = sP.derivPath + 'galradii.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.sav'
    
    if (file_test(saveFilename)) then begin
      print,'Skipping: '+strmid(saveFilename,strlen(sP.derivPath))
      continue
    endif
    
    ; load galaxy and group membership catalogs
    galcat = galaxyCat(sP=sP)
    
    ; restrict gas particle positions to gal/gmem gas only
    ids  = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    
    ; IMPORTANT! rearrange ids_ind to be in the order of gc.xIDs, need this if we want ids[ids_ind], 
    ; temp[ids_ind], etc to be in the same order as the galaxy catalog id list    
    match,galcat.galaxyIDs,ids,galcat_ind,ids_gal_ind,count=countGal,/sort
    ids_gal_ind = ids_gal_ind[sort(galcat_ind)]
    
    match,galcat.groupmemIDs,ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
    ids_gmem_ind = ids_gmem_ind[sort(galcat_ind)]
    
    ids        = !NULL
    galcat_ind = !NULL
    
    ; calculate radial distances of gas elements to primary and secondary parents
    pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
    
    pos_gal  = pos[*,ids_gal_ind]
    pos_gmem = pos[*,ids_gmem_ind]
    
    pos = !NULL
    
    ; load subhalo catalog for mostBoundParticleID and for priParentIDs
    gc = loadGroupCat(sP=sP,/skipIDs)
    
    ; find group center positions with most bound particles for each group
    subgroupCen = subgroupPosByMostBoundID(sP=sP)

    ; replicate parent IDs (of SECONDARY/direct)
    gcIDs = galCatRepParentIDs(galcat=galcat)
    
    ; calulate radial vector of gas from group center (SEC)
    gal_sec  = periodicDists(subgroupCen[*,gcIDs.gal],pos_gal,sP=sP)
    gmem_sec = periodicDists(subgroupCen[*,gcIDs.gmem],pos_gmem,sP=sP)
    
    ; replicate parent IDs (of PRIMARY/parent)
    priParentIDs = gcIDList(gc=gc,select='pri')
    gc = !NULL
    
    gcIDs = galCatRepParentIDs(galcat=galcat,priParentIDs=priParentIDs)
    
    ; calulate radial vector of gas from group center (PRI)
    gal_pri  = periodicDists(subgroupCen[*,gcIDs.gal],pos_gal,sP=sP)
    gmem_pri = periodicDists(subgroupCen[*,gcIDs.gmem],pos_gmem,sP=sP)

    ; save radial distances (and group centers)
    save,gal_pri,gal_sec,gmem_pri,gmem_sec,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))

  endfor
  
end

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

function maxTemps, sP=sP, zStart=zStart, saveRedshifts=saveRedshifts

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
  saveFilename = sP.derivPath + 'maxtemp.'+str(sP.res)+'.'+str(minSnap)+'-'+str(maxSnap)+'.sav'
  
  if (file_test(saveFilename)) then begin
    restore, saveFilename
    return, r
  endif
  
  ; load galaxy/group member catalogs at zMin for gas ids to search for
  galcat = galaxyCat(sP=sP) ;sP.snap is still at zMin
  
  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    print,'Calculating new maxtemp using ( SPH Particles ) res = '+str(sP.res)+$
      ' in range ['+str(minSnap)+'-'+str(maxSnap)+'].'
      
    ; arrays store exactly one quantity per galaxy/groupmem gas element
    maxTemps_gal      = fltarr(n_elements(galcat.galaxyIDs))
    maxTempTime_gal   = fltarr(n_elements(galcat.galaxyIDs))
    maxTemps_gmem     = fltarr(n_elements(galcat.groupmemIDs))
    maxTempTime_gmem  = fltarr(n_elements(galcat.groupmemIDs))
    
    ; no dispersion estimate for the maximum temperatures available
    maxTempDisp_gal  = []
    maxTempDisp_gmem = []
    
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
      w1 = where(temp_gal gt maxTemps_gal and sfr_gal eq 0.0,count1)
      if (count1 gt 0) then begin
        maxTemps_gal[w1]    = temp_gal[w1]
        maxTempTime_gal[w1] = snapNumToRedshift(sP=sP,/time)
      endif
      
      ; replace existing values if current snapshot has higher temps (group members)
      w2 = where(temp_gmem gt maxTemps_gmem and sfr_gmem eq 0.0,count2)
      if (count2 gt 0) then begin
        maxTemps_gmem[w2]    = temp_gmem[w2]
        maxTempTime_gmem[w2] = snapNumToRedshift(sP=sP,/time)
      endif
      
      ; SAVE?
      if total(sP.snap eq saveSnaps) gt 0 then begin
        ; set savefilename
        w = where(sP.snap eq saveSnaps)
        saveFilename = sP.derivPath + 'maxtemp.'+str(sP.res)+'.'+str(minSnap)+'-'+str(saveSnaps[w[0]])+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: saveRedshifts but ['+saveFilename+'] exists, skipping write!'
          continue
        endif
      
        ; save for future lookups
        r = {maxTemps_gal:maxTemps_gal,maxTemps_gmem:maxTemps_gmem,$
             maxTempTime_gal:maxTempTime_gal,maxTempTime_gmem:maxTempTime_gmem,$
             maxTempDisp_gal:maxTempDisp_gal,maxTempDisp_gmem:maxTempDisp_gmem}
             
        save,r,filename=saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        r = !NULL
      endif ; save
    endfor ;m
  endif
  
  ; MONTE CARLO TRACERS CASE - for each original gas cell, determine a single statistic of its 
  ; population of tracers and an estimate for the dispersion in that statistic
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
    
    ; arrays store all tracers (compact into max and dispersion at the end)
    maxTemps_gal      = fltarr(n_elements(galcat_gal_trids))
    maxTempTime_gal   = fltarr(n_elements(galcat_gal_trids))
    maxTemps_gmem     = fltarr(n_elements(galcat_gmem_trids))
    maxTempTime_gmem  = fltarr(n_elements(galcat_gmem_trids))  

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
      
      ; replace existing values if current snapshot has higher temps (enforce off effective EOS)
      w1 = where(temp_gal gt maxTemps_gal,count1)
      if (count1 gt 0) then begin
        maxTemps_gal[w1]    = temp_gal[w1]
        maxTempTime_gal[w1] = temp_time_gal[w1] ; sub-snapshot timing
      endif
      
      ; replace existing values if current snapshot has higher temps (group members)
      w2 = where(temp_gmem gt maxTemps_gmem,count2)
      if (count2 gt 0) then begin
        maxTemps_gmem[w2]    = temp_gmem[w2]
        maxTempTime_gmem[w2] = temp_time_gmem[w2] ; sub-snapshot timing
      endif
      
      temp_gal  = !NULL
      temp_gmem = !NULL
      temp_time_gal  = !NULL
      temp_time_gmem = !NULL
      
      ; SAVE?
      if total(sP.snap eq saveSnaps) gt 0 then begin
        ; set savefilename
        w = where(sP.snap eq saveSnaps)
        saveFilename = sP.derivPath + 'maxtemp.'+str(sP.res)+'.'+str(minSnap)+'-'+str(saveSnaps[w[0]])+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: saveRedshifts but ['+saveFilename+'] exists, skipping write!'
          continue
        endif
        
        ; prepare output
        offset = 0L
        
        r = {maxTemps_gal     : fltarr(n_elements(galcat_gal_cc))  ,$
             maxTemps_gmem    : fltarr(n_elements(galcat_gmem_cc)) ,$
             maxTempTime_gal  : fltarr(n_elements(galcat_gal_cc))  ,$
             maxTempTime_gmem : fltarr(n_elements(galcat_gmem_cc)) ,$
             maxTempDisp_gal  : fltarr(n_elements(galcat_gal_cc))  ,$
             maxTempDisp_gmem : fltarr(n_elements(galcat_gmem_cc))  }
        
        for i=0,n_elements(galcat_gal_cc)-1 do begin
          if galcat_gal_cc[i] gt 0 then begin
            ; for each gas cell, collect its child tracers
            childInds = lindgen(galcat_gal_cc[i]) + offset
            
            locMaxTemps = maxTemps_gal[childInds]
            locMaxTempTime = maxTempTime_gal[childInds] ; this is really time (scale factor) for tracerMC
            
            ; save the maximum value reached by any child tracer and the stddev of the population maxima
            r.maxTemps_gal[i]    = max(locMaxTemps,indmax)
            r.maxTempTime_gal[i] = locMaxTempTime[indmax]
            r.maxTempDisp_gal[i] = stddev(locMaxTemps)
            
            offset += galcat_gal_cc[i]
          endif
        endfor
        
        offset = 0L
        
        for i=0,n_elements(galcat_gmem_cc)-1 do begin
          if galcat_gmem_cc[i] gt 0 then begin
            ; for each gas cell, collect its child tracers
            childInds = lindgen(galcat_gmem_cc[i]) + offset
            
            locMaxTemps = maxTemps_gmem[childInds]
            locMaxTempTime = maxTempTime_gmem[childInds]
            
            ; save the maximum value reached by any child tracer and the stddev of the population maxima
            r.maxTemps_gmem[i]    = max(locMaxTemps,indmax)
            r.maxTempTime_gmem[i] = locMaxTempTime[indmax]
            r.maxTempDisp_gmem[i] = stddev(locMaxTemps)
            
            offset += galcat_gmem_cc[i]
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
    stop ; TODO START - similar to tracerMC since there could be multiple  
  
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
      
      ; find gas parents of both trids_gal and trids_gmem
      tr_parinds = cosmoTracerVelParents(sP=sP,/getInds)
      
      ; use map to convert both tr_parids_gal and tr_parids_gmem to indices
      tr_parinds_gal  = tr_parinds[trids_gal_ind]
      tr_parinds_gmem = tr_parinds[trids_gmem_ind]
      tr_parinds = !NULL
      
      ; load u,nelec to calculate temperatures
      u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='ne')
      
      u_gal  = u[tr_parinds_gal]
      u_gmem = u[tr_parinds_gmem]
      u = !NULL
      
      nelec_gal  = nelec[tr_parinds_gal]
      nelec_gmem = nelec[tr_parinds_gmem]
      nelec = !NULL
      
      temp_gal  = convertUtoTemp(u_gal,nelec_gal)
      temp_gmem = convertUtoTemp(u_gal,nelec_gal)
      
      ; load gas SFR and select off the effective EOS (SFR must be zero)
      sfr = loadSnapshotSubset(sP=sP,partType='gas',field='sfr')
      
      sfr_gal  = sfr[tr_parinds_gal]
      sfr_gmem = sfr[tr_parinds_gmem]
      sfr = !NULL
      
      ; take log of temperatures
      w = where(temp_gal le 0,count_gal)
      if (count_gal ne 0) then temp_gal[w] = 1.0
      temp_gal = alog10(temp_gal)
      
      w = where(temp_gmem le 0,count_gmem)
      if (count_gmem ne 0) then temp_gmem[w] = 1.0
      temp_gmem = alog10(temp_gmem)
      
      ; replace existing values if current snapshot has higher temps (enforce off effective EOS)
      w1 = where(temp_gal gt maxTemps_gal and sfr_gal eq 0.0,count1)
      if (count1 gt 0) then begin
        maxTemps_gal[w1]    = temp_gal[w1]
        maxTempTime_gal[w1] = snapNumToRedshift(sP=sP,/time)
      endif
      
      ; replace existing values if current snapshot has higher temps (group members)
      w2 = where(temp_gmem gt maxTemps_gmem and sfr_gmem eq 0.0,count2)
      if (count2 gt 0) then begin
        maxTemps_gmem[w2]    = temp_gmem[w2]
        maxTempTime_gmem[w2] = snapNumToRedshift(sP=sP,/time)
      endif
     
      ; SAVE?
      if total(sP.snap eq saveSnaps) gt 0 then begin
        ; set savefilename
        w = where(sP.snap eq saveSnaps)
        saveFilename = sP.derivPath + 'maxtemp.'+str(sP.res)+'.'+str(minSnap)+'-'+str(saveSnaps[w[0]])+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: saveRedshifts but ['+saveFilename+'] exists, skipping write!'
          continue
        endif
        
        stop ;TODO FINISH similar to tracerMC since there could be multiple
        
        ; save for future lookups
        r = {maxTemps_gal:maxTemps_gal,maxTemps_gmem:maxTemps_gmem,$
             maxTempTime_gal:maxTempTime_gal,maxTempTime_gmem:maxTempTime_gmem,$
             maxTempDisp_gal:maxTempDisp_gal,maxTempDisp_gmem:maxTempDisp_gmem}
             
        save,r,filename=saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        r = !NULL
      endif ;save
    endfor ;m
  endif

  return,[]
end

; gcSubsetProp(): read galaxy catalog for a specific subgroup selection (pri,sec,all) and
;                 return properties for each gas element (may or may not depend on parent halo)
;                 
;                 note: target redshift is read from sP.snap and values are returned for that
;                       redshift unless oSnap is specified, in which case gasOrigins() is used and the
;                       gas elements are matched to those existing in the subhalo catalog at the
;                       earlier snapshot
;
; rVirNorm=1    : radial distances normalized by r_vir of either primary or secondary parent
; virTemp=1     : virial temperatures of parent halos
; parMass=1     : total mass (dm+baryonic) of parent halos (from catalog)
; curTemp=1     : current temperature of each element
; curDens=1     : current density of each gas element
; maxPastTemp=1 : maximum past previous temperature of each element

function gcSubsetProp, sP=sP, select=select, oSnap=oSnap, $
                       rVirNorm=rVirNorm, virTemp=virTemp, parMass=parMass, $
                       curTemp=curTemp, curDens=curDens, maxPastTemp=maxPastTemp, $
                       parNorm=parNorm ;this row for rVirNorm only

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; select primary,secondary,or all subhalos subject to minimum number of particles
  gcIDList = gcIDList(sP=sP,select=select)

  ; select galaxycat indices corresponding to this list of subgroup ids
  galcatInds = galcatINDList(sP=sP,gcIDList=gcIDList)

  if keyword_set(rVirNorm) then begin
    ; gasOrigins - return this quantity for the gas elements at the different snapshot oSnap
    if (keyword_set(oSnap)) then begin
    
      snapSwap = sP.snap
      sP.snap = oSnap
      gasOrig = gasOrigins(sP=sP)
      
      ; load parent r_vir and galaxy radii catalog at oSnap     
      r_vir = galCatParentProperties(sP=sP, /rVir)
      gcr   = galaxyCatRadii(sP=sP)
      
      sP.snap = snapSwap
      
      ; arrays
      rad_gal  = fltarr(n_elements(gcr.gal_pri)) - 1.0
      rad_gmem = fltarr(n_elements(gcr.gmem_pri)) - 1.0
      
      ; store radial distance of gas elements to parents for matching IDs
      if (parNorm eq 'pri') then begin
        rad_gal[gasOrig.indMatch.gc_ind_gal]  = gcr.gal_pri[gasOrig.indMatch.gcCur_ind_gal] / $
                                                r_vir.gal[gasOrig.indMatch.gcCur_ind_gal]
                                                
        if (n_elements(gasOrig.indMatch.gc_ind_gal2) gt 0) then $
          rad_gal[gasOrig.indMatch.gc_ind_gal2] = gcr.gmem_pri[gasOrig.indMatch.gcCur_ind_gal2] / $
                                                  r_vir.gmem[gasOrig.indMatch.gcCur_ind_gal2]
        
        rad_gmem[gasOrig.indMatch.gc_ind_gmem]  = gcr.gmem_pri[gasOrig.indMatch.gcCur_ind_gmem] / $
                                                r_vir.gmem[gasOrig.indMatch.gcCur_ind_gmem]
                                                
        if (n_elements(gasOrig.indMatch.gc_ind_gmem2) gt 0) then $
          rad_gmem[gasOrig.indMatch.gc_ind_gmem2] = gcr.gal_pri[gasOrig.indMatch.gcCur_ind_gmem2] / $
                                                  r_vir.gal[gasOrig.indMatch.gcCur_ind_gmem2]
      endif
      
      if (parNorm eq 'sec') then begin
        rad_gal[gasOrig.indMatch.gc_ind_gal]  = gcr.gal_sec[gasOrig.indMatch.gcCur_ind_gal] / $
                                                r_vir.gal[gasOrig.indMatch.gcCur_ind_gal]
        rad_gal[gasOrig.indMatch.gc_ind_gal2] = gcr.gmem_sec[gasOrig.indMatch.gcCur_ind_gal2] / $
                                                r_vir.gmem[gasOrig.indMatch.gcCur_ind_gal2]
        
        rad_gmem[gasOrig.indMatch.gc_ind_gmem]  = gcr.gmem_sec[gasOrig.indMatch.gcCur_ind_gmem] / $
                                                r_vir.gmem[gasOrig.indMatch.gcCur_ind_gmem]
        rad_gmem[gasOrig.indMatch.gc_ind_gmem2] = gcr.gal_sec[gasOrig.indMatch.gcCur_ind_gmem2] / $
                                                r_vir.gal[gasOrig.indMatch.gcCur_ind_gmem2]
      endif
    endif else begin
      ; load parent r_vir and galaxy radii catalog at sP.snap
      r_vir = galCatParentProperties(sP=sP, /rVir)
      gcr   = galaxyCatRadii(sP=sP)
      
      ; store radial distance normalized by parent r_vir
      if (parNorm eq 'pri') then begin
        rad_gal  = gcr.gal_pri / r_vir.gal
        rad_gmem = gcr.gmem_pri / r_vir.gmem
      endif
      
      if (parNorm eq 'sec') then begin
        rad_gal  = gcr.gal_sec / r_vir.gal
        rad_gmem = gcr.gmem_sec / r_vir.gmem
      endif
    endelse
    
    ; restrict to subset
    rad_gal  = rad_gal[galcatInds.gal]
    rad_gmem = rad_gmem[galcatInds.gmem]
    
    r = {gal:rad_gal,gmem:rad_gmem}
  endif
  
  if keyword_set(virTemp) then begin
    ; load parent t_vir
    t_vir = galCatParentProperties(sP=sP, /virTemp)
    
    ; create subsets for subhalo selection
    temp_gal  = t_vir.gal[galcatInds.gal]
    temp_gmem = t_vir.gmem[galcatInds.gmem]
    
    r = {gal:temp_gal,gmem:temp_gmem}
  endif
  
  if keyword_set(parMass) then begin
    ; load parent total mass
    mass = galCatParentProperties(sP=sP, /mass)
    
    ; create subsets for subhalo selection
    mass_gal  = mass.gal[galcatInds.gal]
    mass_gmem = mass.gmem[galcatInds.gmem]
    
    r = {gal:mass_gal,gmem:mass_gmem}
  endif
  
  if keyword_set(curTemp) then begin
    ; gasOrigins - return this quantity for the gas elements at the different snapshot oSnap
    if (keyword_set(oSnap)) then begin
      sP.snap = oSnap
      gasOrig = gasOrigins(sP=sP)
      
      temp_gal  = gasOrig.temp_gal[galcatInds.gal]
      temp_gmem = gasOrig.temp_gmem[galcatInds.gmem]
    endif else begin
      ; load gas ids and make ID->index map
      ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
      idsIndMap = getIDIndexMap(ids,minid=minid)
      ids = !NULL
      
      ; load galaxy catalog to change INDs to gas IDs
      galcat = galaxyCat(sP=sP)
      
      ; load gas u and restrict to subset of galaxy cat
      u      = loadSnapshotSubset(sP=sP,partType='gas',field='u')
      u_gal  = u[idsIndMap[galcat.galaxyIDs[galcatInds.gal]-minid]]
      u_gmem = u[idsIndMap[galcat.groupmemIDs[galcatInds.gmem]-minid]]
      u      = !NULL
      
      ; load gas nelec and restrict to subset of galaxy cat
      nelec      = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
      nelec_gal  = nelec[idsIndMap[galcat.galaxyIDs[galcatInds.gal]-minid]]
      nelec_gmem = nelec[idsIndMap[galcat.groupmemIDs[galcatInds.gmem]-minid]]
      nelec      = !NULL
      
      idsIndMap = !NULL
      galcat    = !NULL
      
      ; calculate temperature
      temp_gal  = convertUtoTemp(u_gal,nelec_gal)
      temp_gmem = convertUtoTemp(u_gmem,nelec_gmem)
        
      ; take log of temperatures
      w = where(temp_gal le 0,count)
      if (count ne 0) then temp_gal[w] = 1.0
      w = where(temp_gmem le 0,count)
      if (count ne 0) then temp_gmem[w] = 1.0
    endelse
  
    temp_gal  = alog10(temp_gal)
    temp_gmem = alog10(temp_gmem)
    
    r = {gal:temp_gal,gmem:temp_gmem}
  endif
  
  if keyword_set(curDens) then begin
    ; load gas ids and make ID->index map
    ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    idsIndMap = getIDIndexMap(ids,minid=minid)
    ids = !NULL
    
    ; load gas densities
    density = loadSnapshotSubset(sP=sP,partType='gas',field='dens')
    
    ; load galaxy catalog to change INDs to gas IDs
    galcat = galaxyCat(sP=sP)

    ; restrict densities to subset of galaxy cat
    dens_gal  = density[idsIndMap[galcat.galaxyIDs[galcatInds.gal]-minid]]
    dens_gmem = density[idsIndMap[galcat.groupmemIDs[galcatInds.gmem]-minid]]

    galcat    = !NULL
    density   = !NULL
    idsIndMap = !NULL
      
    r = {gal:dens_gal,gmem:dens_gmem}
  endif
  
  if keyword_set(maxPastTemp) then begin
    ; load maximum past temperature
    maxt = maxTemps(sP=sP)
    
    ; restrict temps to subset of galaxy cat
    temp_gal  = maxt.maxTemps_gal[galcatInds.gal]
    temp_gmem = maxt.maxTemps_gmem[galcatInds.gmem]
    
    maxt = !NULL
    
    r = {gal:temp_gal,gmem:temp_gmem}
  endif
  
  return,r
end
