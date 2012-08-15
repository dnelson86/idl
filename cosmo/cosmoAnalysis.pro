; cosmoAnalysis.pro
; gas accretion project - main
; dnelson aug.2012

; galaxyCat(): if snap not specified, create and save complete galaxy catalog from the group catalog by 
;              imposing additional cut in the (rho,temp) plane (same as that used by Torrey+ 2011)
;              if snap is specified, create only for one snapshot number or return previously saved
;              results for that snapshot
; Note: the gal/gmem catalogs have the same size (and indexing) as the subgroups

function galaxyCat, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; if snap specified, run only one snapshot (and/or just return previous results)
  if (sP.snap ne -1) then begin
    saveFilename1 = sP.derivPath + 'galcat.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) +'.sav'
    saveFilename2 = sP.derivPath + 'groupmemcat.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.sav'
    saveFilename3 = sP.derivPath + 'starcat.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.sav'
    
    ; results exist, return
    if (file_test(saveFilename1) and file_test(saveFilename2) and file_test(saveFilename3)) then begin
      restore,saveFilename1
      restore,saveFilename2
      restore,saveFilename3
      r = {galaxyOff:galaxyOff,galaxyLen:galaxyLen,galaxyIDs:galaxyIDs,$
           stellarOff:stellarOff,stellarLen:stellarLen,stellarIDs:stellarIDs,$
           groupmemOff:groupmemOff,groupmemLen:groupmemLen,groupmemIDs:groupmemIDs}
      return,r
    endif
    
    ; need to compute, set restricted range of snapshots to process
    snapRange = [sP.snap,sP.snap]
  endif else begin
    snapRange = sP.groupCatRange
  endelse
  
  for m=snapRange[0],snapRange[1],1 do begin
    sP.snap = m
    ; skip if previous results exist
    saveFilename1 = sP.derivPath + 'galcat.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.sav'
    saveFilename2 = sP.derivPath + 'groupmemcat.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.sav'
    saveFilename3 = sP.derivPath + 'starcat.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.sav'
    
    if (file_test(saveFilename1) and file_test(saveFilename2) and file_test(saveFilename3)) then begin
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

    if countMatch ne n_elements(gcPIDs) then message,'Error: Failed to locate all gas gcPIDs in gas ids.'

    gcPIDs = !NULL
    ids = ids[ids_ind]

    ; load star ids and match to catalog
    gcPIDs_stars = gcPIDList(gc=gc,select='all',partType='stars')
    ids_stars = loadSnapshotSubset(sP=sP,partType='stars',field='ids')
    match,gcPIDs_stars,ids_stars,gc_ind,ids_ind_stars,count=countMatch,/sort
    ids_ind_stars = ids_ind_stars[sort(gc_ind)] ; IMPORTANT! rearrange ids_ind_stars

    if countMatch ne n_elements(gcPIDs_stars) then message,'Error: Failed to locate all star gcPIDs in gas ids.'

    gcPIDs = !NULL
    gc_ind = !NULL
    ids_stars = ids_stars[ids_ind_stars]

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
    
    ; (rho,temp) cut
    wGal = where(alog10(temp) - sP.galcut_rho * alog10(dens) lt sP.galcut_T,$
                 countGal,comp=wGmem,ncomp=countGmem)

    ; calculate radial distances of gas elements to primary parents
    if sP.radcut_rvir ne 0.0 then begin
      print,'GALCAT: ENFORCING RADIAL CUT (AND STARS)!'
      
      pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
      pos = pos[*,ids_ind]
      
      ; find group center positions with most bound particles for each group
      subgroupCen = subgroupPosByMostBoundID(sP=sP)
  
      ; create PRIMARY parent -subgroup- ID list
      if gc.nSubgroupsTot ne total(gc.groupNsubs,/int) then message,'Error: Subgroup counts mismatch.'
      parIDs = lonarr(gc.nSubgroupsTot)
      offset = 0L
      
      for gID=0L,gc.nGroupsTot-1 do begin
        if gc.groupNSubs[gID] gt 0 then begin
          parIDs[offset:offset+gc.groupNSubs[gID]-1] = cmreplicate(gc.groupFirstSub[gID],gc.groupNSubs[gID])
          offset += gc.groupNSubs[gID]
        endif
      endfor
  
      ; replicate group parent IDs (of PRIMARY/parent) to each member particle
      ptNum = { gas : partTypeNum('gas'), stars : partTypeNum('stars') }
      
      sgParIDs_gas   = lonarr(total(gc.subgroupLenType[ptNum.gas,*],/int))
      sgParIDs_stars = lonarr(total(gc.subgroupLenType[ptNum.stars,*],/int))
      offset = { gas : 0L, stars: 0L }
      
      for sgID=0L,gc.nSubgroupsTot-1 do begin
        if gc.subgroupLenType[ptNum.gas,sgID] gt 0 then begin
          sgParIDs_gas[offset.gas:offset.gas+gc.subgroupLenType[ptNum.gas,sgID]-1] = $
            cmreplicate(parIDs[sgID],gc.subgroupLenType[ptNum.gas,sgID])
          offset.gas += gc.subgroupLenType[ptNum.gas,sgID]
        endif
        
        if gc.subgroupLenType[ptNum.stars,sgID] gt 0 then begin
          sgParIDs_stars[offset.stars:offset.stars+gc.subgroupLenType[ptNum.stars,sgID]-1] = $
            cmreplicate(parIDs[sgID],gc.subgroupLenType[ptNum.stars,sgID])
          offset.stars += gc.subgroupLenType[ptNum.stars,sgID]
        endif
      endfor
      
      if offset.gas ne n_elements(sgParIDs_gas) then message,'Error: Bad parent ID replication.'
      if offset.stars ne n_elements(sgParIDs_stars) then message,'Error: Bad parent ID star rep.'
      
      ; fof has rvir=0 (no SO values) if zero subgroups, marginal overdensity, or low total mass
      ; these are low mass halos which we aren't going to plot anyways
      w = where(gc.group_r_crit200 eq 0.0,count)
      if count gt 0 then gc.group_r_crit200[w] = 1e-8 ; remove with outer radial cut
      
      ; calculate radial vector of gas from group center (PRI)
      rad_pri  = periodicDists(subgroupCen[*,sgParIDs_gas],pos,sP=sP)
      par_rvir = gc.group_r_crit200[gc.subgroupGrNr[sgParIDs_gas]] ; normalize by fof parent rvir

      rad_pri /= par_rvir

      ; override (rho,temp) cut with (rho,temp,rad) cut
      wGal  = where(alog10(temp) - sP.galcut_rho*alog10(dens) lt sP.galcut_T and rad_pri lt sP.radcut_rvir,countGal)
      wGmem = where(alog10(temp) - sP.galcut_rho*alog10(dens) ge sP.galcut_T and $
                     rad_pri ge sP.radcut_rvir and rad_pri le sP.radcut_out,countGmem)
                     
      ; load stellar positions calculate radial vector of stars from group center (PRI)
      pos = loadSnapshotSubset(sP=sP,partType='stars',field='pos')
      pos = pos[*,ids_ind_stars]
      
      rad_pri  = periodicDists(subgroupCen[*,sgParIDs_stars],pos,sP=sP)
      par_rvir = gc.group_r_crit200[gc.subgroupGrNr[sgParIDs_stars]] ; normalize by fof parent rvir
      
      rad_pri /= par_rvir
      
      ; make stellar cut
      wStars = where(rad_pri lt sP.radcut_rvir,countStars)
      
    endif

    if countGal eq 0 or countGmem eq 0 then begin
      print,'Warning: Empty galaxy cut or comp. Skipping: ' + strmid(saveFilename1,strlen(sP.derivPath))
      continue
    endif
    
    temp = !NULL
    dens = !NULL
    
    ; make subsets of ids matching galaxy cut and complement
    ids_groupmem = ids[wGmem]
    ids          = ids[wGal]
    ids_stars    = ids_stars[wStars]

    ; construct group member catalog
    if (not file_test(saveFilename2)) then begin
    
      groupmemLen = lonarr(gc.nSubgroupsTot)
      groupmemOff = lonarr(gc.nSubgroupsTot)
      groupmemIDs = lon64arr(n_elements(ids_groupmem))
      
      nextOff = 0UL
      
      ; match indices between gas ids and group member ids
      match,ids_groupmem,gc.IDs,ids_ind,gcIDs_ind,count=countID,/sort
      
      for gcID=0L,gc.nSubgroupsTot-1 do begin
        ; select ids in group
        groupStart = gc.subGroupOffset[gcID]
        groupEnd   = groupStart + gc.subGroupLen[gcID]
        
        w = where(gcIDs_ind ge groupStart and gcIDs_ind lt groupEnd,count)
        
        ; save in similar 1D offset format
        groupmemLen[gcID] = count
        groupmemOff[gcID] = nextOff
        if (count gt 0) then begin
          nextOff = groupmemOff[gcID] + groupmemLen[gcID]
          groupmemIDs[groupmemOff[gcID]:nextOff-1] = ids_groupmem[ids_ind[w]]
        endif      
        
      end
  
      ; make sure all gas particles were found in the group catalog
      match,ids_groupmem,groupmemIDs,ind1,ind2,count=countCheck
      if countCheck ne n_elements(ids_groupmem) then message,'Error: Failed to locate all gmem.'

      ; save group membership catalog
      save,groupmemLen,groupmemOff,groupmemIDs,filename=saveFilename2    
      print,'Saved: '+strmid(saveFilename2,strlen(sP.derivPath))+' ['+str(countGmem)+'/'+str(countMatch)+']'    
    
      groupmemLen = !NULL
      groupmemOff = !NULL
      groupmemIDs = !NULL
    
    endif

    ; construct galaxy catalog
    if (not file_test(saveFilename1)) then begin
      galaxyLen = lonarr(gc.nSubgroupsTot)
      galaxyOff = lonarr(gc.nSubgroupsTot)
      galaxyIDs = lon64arr(n_elements(ids))
      
      nextOff = 0UL
      
      ; match indices between galaxy ids and group ids
      match,ids,gc.IDs,ids_ind,gcIDs_ind,count=countID,/sort
      
      for gcID=0L,gc.nSubgroupsTot-1 do begin
        ; select ids in group
        groupStart = gc.subGroupOffset[gcID]
        groupEnd   = groupStart + gc.subGroupLen[gcID]
        
        w = where(gcIDs_ind ge groupStart and gcIDs_ind lt groupEnd,count)
        
        ; save in similar 1D offset format
        galaxyLen[gcID] = count
        galaxyOff[gcID] = nextOff
        if (count gt 0) then begin
          nextOff = galaxyOff[gcID] + galaxyLen[gcID]
          galaxyIDs[galaxyOff[gcID]:nextOff-1] = ids[ids_ind[w]]
        endif      
        
      end
      
      ; make sure all gas particles were found in the group catalog
      match,ids,galaxyIDs,ind1,ind2,count=countCheck
      if countCheck ne n_elements(ids) then message,'Error: Failed to locate all gal.'

      ; save galaxy catalog
      save,galaxyLen,galaxyOff,galaxyIDs,filename=saveFilename1
      print,'Saved: '+strmid(saveFilename1,strlen(sP.derivPath))+' ['+str(countGal)+'/'+str(countMatch)+']'
    
      galaxyLen = !NULL
      galaxyOff = !NULL
      galaxyIDs = !NULL
      
    endif
    
    ; construct stellar catalog
    if (not file_test(saveFilename3)) then begin
      stellarLen = lonarr(gc.nSubgroupsTot)
      stellarOff = lonarr(gc.nSubgroupsTot)
      stellarIDs = lon64arr(n_elements(ids_stars))
      
      nextOff = 0UL
      
      ; match indices between star ids and group ids
      match,ids_stars,gc.IDs,ids_ind,gcIDs_ind,count=countID,/sort
      
      for gcID=0L,gc.nSubgroupsTot-1 do begin
        ; select ids in group
        groupStart = gc.subGroupOffset[gcID]
        groupEnd   = groupStart + gc.subGroupLen[gcID]
        
        w = where(gcIDs_ind ge groupStart and gcIDs_ind lt groupEnd,count)
        
        ; save in similar 1D offset format
        stellarLen[gcID] = count
        stellarOff[gcID] = nextOff
        if (count gt 0) then begin
          nextOff = stellarOff[gcID] + stellarLen[gcID]
          stellarIDs[stellarOff[gcID]:nextOff-1] = ids_stars[ids_ind[w]]
        endif      
        
      end
      
      ; make sure all stellar particles were found in the group catalog
      match,ids_stars,stellarIDs,ind1,ind2,count=countCheck
      if countCheck ne n_elements(ids_stars) then message,'Error: Failed to locate all stars.'

      ; save stellar catalog
      save,stellarLen,stellarOff,stellarIDs,filename=saveFilename3
      print,'Saved: '+strmid(saveFilename3,strlen(sP.derivPath))+' ['+str(countStars)+'/'+str(countMatch)+']'
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
      r = {gal_pri:gal_pri,gal_sec:gal_sec,gmem_pri:gmem_pri,gmem_sec:gmem_sec,$
           stars_pri:stars_pri,stars_sec:stars_sec}
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
    
    ; load stellar ids and match, calculate radial distances
    ids = loadSnapshotSubset(sP=sP,partType='stars',field='ids')
    
    match,galcat.stellarIDs,ids,galcat_ind,ids_stars_ind,count=countStars,/sort
    ids_stars_ind = ids_stars_ind[sort(galcat_ind)]
    
    ids        = !NULL
    galcat_ind = !NULL
    
    pos_stars = loadSnapshotSubset(sP=sP,partType='stars',field='pos')
    pos_stars = pos_stars[*,ids_stars_ind]
    
    ; load subhalo catalog for mostBoundParticleID and for priParentIDs
    gc = loadGroupCat(sP=sP,/skipIDs)
    
    ; find group center positions with most bound particles for each group
    subgroupCen = subgroupPosByMostBoundID(sP=sP)

    ; replicate parent IDs (of SECONDARY/direct)
    gcIDs = galCatRepParentIDs(galcat=galcat)
    
    ; calulate radial vector of gas from group center (SEC)
    gal_sec   = periodicDists(subgroupCen[*,gcIDs.gal],pos_gal,sP=sP)
    gmem_sec  = periodicDists(subgroupCen[*,gcIDs.gmem],pos_gmem,sP=sP)
    stars_sec = periodicDists(subgroupCen[*,gcIDs.stars],pos_stars,sP=sP)
    
    ; replicate parent IDs (of PRIMARY/parent)
    priParentIDs = gcIDList(gc=gc,select='pri')
    gc = !NULL
    
    gcIDs = galCatRepParentIDs(galcat=galcat,priParentIDs=priParentIDs)
    
    ; calulate radial vector of gas from group center (PRI)
    gal_pri   = periodicDists(subgroupCen[*,gcIDs.gal],pos_gal,sP=sP)
    gmem_pri  = periodicDists(subgroupCen[*,gcIDs.gmem],pos_gmem,sP=sP)
    stars_pri = periodicDists(subgroupCen[*,gcIDs.stars],pos_stars,sP=sP)

    ; save radial distances (and group centers)
    save,gal_pri,gal_sec,gmem_pri,gmem_sec,stars_pri,stars_sec,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))

  endfor
  
end

; gcSubsetProp(): read galaxy catalog for a specific subgroup selection (pri,sec,all) and
;                 return properties for each gas element/tracer (may or may not depend on parent halo)
;                 at the redshift specified by sP.snap (may or may not depend on previous time)
;
; rVirNorm=1    : radial distances normalized by r_vir of either primary or secondary parent
;  parNorm      : either 'pri' or 'sec' if rVirNorm requested
; virTemp=1     : current virial temperatures of parent halos
; parMass=1     : total mass (dm+baryonic) of parent halos (from catalog)
; curTemp=1     : current temperature of each element
; maxPastTemp=1 : maximum past previous temperature of each element
; maxTempTime=1 : time when maximum past previous temperature was reached (in redshift)
;  trPopMin,trPopMean,trPopMax : return the respective statistic for each gas cell for the tracers
; elemIDs=1     : ids of each element (either SPH particles or tracers)
; 
; curSingleVal=1  : current single quantity (e.g. mass, density) returned without manipulation
;  singleValField : field name in snapshot file for the above
;
; mergerTreeSubset    : return values only for the subset of halos tracked in the merger tree subset
; accretionTimeSubset : return values only for the subset of particles/tracers with recorded accretion times
;  accTime,accTvir : time of accretion (in redshift) or virial temp of parent halo at time of accretion
;  accMode : return values only for one accretionMode (all,smooth,bclumpy,sclumpy,smooth)

function gcSubsetProp, sP=sP, select=select, $
           rVirNorm=rVirNorm, virTemp=virTemp, parMass=parMass, $
           curTemp=curTemp, maxPastTemp=maxPastTemp, maxTempTime=maxTempTime, elemIDs=elemIDs, $
           trPopMin=trPopMin, trPopMax=trPopMax, trPopMean=trPopMean, $ ; for maxPastTemp/maxTempTime only
           curSingleVal=curSingleVal, singleValField=singleValField, $
           parNorm=parNorm, $ ; for rVirNorm,virTemp,parMass only
           mergerTreeSubset=mergerTreeSubset, accretionTimeSubset=accretionTimeSubset,$
           accTime=accTime,accTvir=accTvir,accMode=accMode ; for accretionTimeSubset only

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; check combinations of input options for validity
  if keyword_set(mergerTreeSubset) and select ne 'pri' then $
    print,'Warning: The merger tree subset actually contains only pri subgroups.'
  if (keyword_set(accTime) or keyword_set(accTvir)) and ~keyword_set(accretionTimeSubset) then $
    message,'Error: Can only return accretion time or Tvir at accretion time for atS.'
  if keyword_set(accretionTimeSubset) and ~keyword_set(mergerTreeSubset) then $
    message,'Error: Can only subselect the mtS with the atS.'
  if keyword_set(accMode) and ~keyword_set(accretionTimeSubset) then $
    message,'Error: Can only return accretion mode subsets of the accretionTime subset.'
  if keyword_set(elemIDs) and (keyword_set(trPopMin) or keyword_set(trPopMean) or keyword_set(trPopMax)) then $
    message,'Error: Cannot return pop stats of unique element IDs.'
    
  ; default behavior: return 1 value per tracer for tracer sims, 1 value per gas particle for sph
  allTR = 0 ; sph
  if sP.trMCPerCell ne 0 then allTR = 1
  
  ; select primary,secondary,or all subhalos subject to minimum number of particles
  gcIDList = gcIDList(sP=sP,select=select)

  ; subset subgroup id list by only those with good merger histories
  if keyword_set(mergerTreeSubset) then begin
    mt = mergerTreeSubset(sP=sP)
    
    ; use intersection of (pri,sec,all) list and tracked list
    match,gcIDList,mt.galcatIDList,ind1,ind2,count=count,/sort
    if count eq 0 then message,'Error: mtS and gcIDList intersection empty.'
    ;print,'gcSubsetProp intersect',n_elements(gcIDList),n_elements(mt.galcatIDList),count
    
    gcIDList = gcIDList[ind1]
    
    mt   = !NULL
    ind1 = !NULL
    ind2 = !NULL
  endif
  
  ; select galaxycat indices corresponding to the list of subgroup ids
  galcatInds = galcatINDList(sP=sP,gcIDList=gcIDList) ;identical to mt.galcatSub if mtS

  ; subset galcat member indlist by those with recorded accretion times (and associated properties)
  if keyword_set(accretionTimeSubset) then begin
    at = accretionTimes(sP=sP)
    
    ; this is used to access accretionTimes
    if n_elements(accMode) eq 0 then accMode = 'all'
    accTimeInds = accModeInds(at=at,sP=sP,accMode=accMode)

    ; sph case: modify galcatInds such that the accretionTimes subset is taken
    if ~allTR then galcatInds = { gal   : galcatInds.gal[accTimeInds.gal]    ,$
                                  gmem  : galcatInds.gmem[accTimeInds.gmem]  ,$
                                  stars : galcatInds.stars[accTimeInds.stars] }

    ; tracer case: handle only after we have child counts (after expansion or in allTR for maxTemps)
  endif
  
  ; ----- values -----
  
  if keyword_set(accTime) then begin
    ;convert scale factors -> redshift
    r = { gal   : 1/at.AccTime_gal[0,accTimeInds.gal]-1    ,$
          gmem  : 1/at.AccTime_gmem[0,accTimeInds.gmem]-1  ,$
          stars : 1/at.AccTime_stars[0,accTimeInds.stars]-1 } 
    return,r
  endif
  
  if keyword_set(accTvir) then begin
    ; take accretionTime subset and return
    r = { gal   : at.AccHaloTvir_gal[accTimeInds.gal]    ,$
          gmem  : at.AccHaloTvir_gmem[accTimeInds.gmem]  ,$
          stars : at.AccHaloTvir_stars[accTimeInds.stars] } 
    return,r
  endif

  if keyword_set(rVirNorm) then begin
    ; load parent r_vir and galaxy radii catalog at sP.snap
    r_vir = galCatParentProperties(sP=sP, /rVir)
    gcr   = galaxyCatRadii(sP=sP)

    ; store radial distance normalized by parent r_vir
    ; note: if SO values not calculated (no subgroup in fof group), rvir=0 and rad->Inf (not plotted)
    ; note: since using most bound particle for subgroup centers, one per subgroup will have r=0
    if (parNorm eq 'pri') then begin
      rad_gal   = gcr.gal_pri / r_vir.gal
      rad_gmem  = gcr.gmem_pri / r_vir.gmem
      rad_stars = gcr.stars_pri / r_vir.stars
    endif
    
    if (parNorm eq 'sec') then begin
      rad_gal   = gcr.gal_sec / r_vir.gal
      rad_gmem  = gcr.gmem_sec / r_vir.gmem
      rad_stars = gcr.stars_sec / r_vir.stars
    endif
      
    ; restrict to subgroup type / merger tree / accretion time subset
    r = { gal   : rad_gal[galcatInds.gal]    ,$
          gmem  : rad_gmem[galcatInds.gmem]  ,$
          stars : rad_stars[galcatInds.stars] }
  endif
  
  if keyword_set(virTemp) then begin
    ; load parent t_vir
    t_vir = galCatParentProperties(sP=sP, /virTemp, parNorm=parNorm)
    
    ; restrict to subgroup type / merger tree / accretion time subset
    r = { gal   : t_vir.gal[galcatInds.gal]    ,$
          gmem  : t_vir.gmem[galcatInds.gmem]  ,$
          stars : t_vir.stars[galcatInds.stars] }
  endif
  
  if keyword_set(parMass) then begin
    ; load parent total mass
    mass = galCatParentProperties(sP=sP, /mass, parNorm=parNorm)
    
    ; restrict to subgroup type / merger tree / accretion time subset
    r = { gal   : mass.gal[galcatInds.gal]    ,$
          gmem  : mass.gmem[galcatInds.gmem]  ,$
          stars : mass.stars[galcatInds.stars] }
  endif
  
  if keyword_set(curTemp) then begin
    ; load gas ids and make ID->index map
    ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    idsIndMap = getIDIndexMap(ids,minid=minid)
    ids = !NULL
    
    ; load galaxy catalog to change INDs to gas IDs
    galcat = galaxyCat(sP=sP)
    
    ; load gas u and restrict to subset of galaxy cat
    u       = loadSnapshotSubset(sP=sP,partType='gas',field='u')
    u_gal   = u[idsIndMap[galcat.galaxyIDs[galcatInds.gal]-minid]]
    u_gmem  = u[idsIndMap[galcat.groupmemIDs[galcatInds.gmem]-minid]]
    u       = !NULL
    
    ; load gas nelec and restrict to subset of galaxy cat
    nelec       = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
    nelec_gal   = nelec[idsIndMap[galcat.galaxyIDs[galcatInds.gal]-minid]]
    nelec_gmem  = nelec[idsIndMap[galcat.groupmemIDs[galcatInds.gmem]-minid]]
    nelec       = !NULL
    
    idsIndMap = !NULL
    
    ; calculate temperature
    temp_gal   = convertUtoTemp(u_gal,nelec_gal,/log)
    temp_gmem  = convertUtoTemp(u_gmem,nelec_gmem,/log)
    temp_stars = fltarr(n_elements(galcat.stellarIDs[galcatInds.stars])) ; zeros, stars have no current temperature
     
    galcat = !NULL
     
    r = {gal:temp_gal,gmem:temp_gmem,stars:temp_stars}
  endif
  
  if keyword_set(curSingleVal) then begin
    ; load gas ids and make ID->index map
    ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    idsIndMap = getIDIndexMap(ids,minid=minid)
    ids = !NULL
    
    ; load gas densities
    singleVal = loadSnapshotSubset(sP=sP,partType='gas',field=singleValField)
    
    ; load galaxy catalog to change INDs to gas IDs
    galcat = galaxyCat(sP=sP)

    ; restrict to subgroup type / merger tree / accretion time subset
    val_gal   = singleVal[idsIndMap[galcat.galaxyIDs[galcatInds.gal]-minid]]
    val_gmem  = singleVal[idsIndMap[galcat.groupmemIDs[galcatInds.gmem]-minid]]
    val_stars = fltarr(n_elements(galcat.stellarIDs[galcatInds.stars])) ; zeros, stars don't have the same values as gas
    
    galcat = !NULL & density = !NULL & idsIndMap = !NULL
    
    r = {gal:val_gal,gmem:val_gmem,stars:val_stars}
  endif
  
  if keyword_set(elemIDs) then begin
    if sP.trMCPerCell eq -1 then message,'Error: Not implemented.'
    
    ; load galaxy catalog for element IDs
    galcat = galaxyCat(sP=sP)
    
    ; if all tracers requested, find children of all the gas IDs in the groupcat
    if keyword_set(allTR) then begin
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
      ids_gal   = tr_ids[cosmoTracerChildren(sP=sP,gasIDs=galcat.galaxyIDs[galcatInds.gal],/getInds)]
      ids_gmem  = tr_ids[cosmoTracerChildren(sP=sP,gasIDs=galcat.groupmemIDs[galcatInds.gmem],/getInds)]
      ids_stars = tr_ids[cosmoTracerChildren(sP=sP,starIDs=galcat.stellarIDs[galcatInds.stars],/getInds)]
      tr_ids = !NULL
      
      r = {gal:ids_gal,gmem:ids_gmem,stars:ids_stars}
      
      ids_gal = !NULL & ids_gmem = !NULL & ids_stars = !NULL
      
      ; take accretionTime subset of mtS of all tracers and return
      if keyword_set(accretionTimeSubset) then begin
        r = { gal:r.gal[accTimeInds.gal], gmem:r.gmem[accTimeInds.gmem], stars:r.stars[acctimeInds.stars] }
      endif
      return,r
    endif
    
    ; otherwise, return ids for each gas particle
    r = { gal   : galcat.galaxyIDs[galcatInds.gal]    ,$
          gmem  : galcat.groupmemIDs[galcatInds.gmem] ,$
          stars : galcat.stellarIDs[galcatInds.stars]  }
  endif
  
  if keyword_set(maxPastTemp) or keyword_set(maxTempTime) then begin
    ; if all tracers requested, load directly and immediately return
    if keyword_set(allTR) then begin
      maxt_gal   = maxTemps(sP=sP,/loadAllTRGal)
      maxt_gmem  = maxTemps(sP=sP,/loadAllTRGmem)
      maxt_stars = maxTemps(sP=sP,/loadAllTRStars)
      
      ; make indices for mergerTreeSubset
      galcatInds = galcatINDList(sP=sP,gcIDList=gcIDList,$
                     child_counts={gal:maxt_gal.child_counts,gmem:maxt_gmem.child_counts,stars:maxt_stars.child_counts})
                       
      ; return temps (logK) or times (converted to redshift)
      if keyword_set(maxPastTemp) then $
        r = {gal:maxt_gal.maxTemps[galcatInds.gal],gmem:maxt_gmem.maxTemps[galcatInds.gmem],stars:maxt_stars.maxTemps[galcatInds.stars]}
      if keyword_set(maxTempTime) then $
        r = {gal:1/maxt_gal.maxTempTime[galcatInds.gal]-1,gmem:1/maxt_gmem.maxTempTime[galcatInds.gmem]-1,stars:1/maxt_stars.maxTempTime[galcatInds.stars]-1}
      
      maxt_gal = !NULL & maxt_gmem = !NULL & maxt_stars = !NULL
      
      ; take accretionTime subset of mtS of all tracers and return
      if keyword_set(accretionTimeSubset) then begin
        r = { gal:r.gal[accTimeInds.gal], gmem:r.gmem[accTimeInds.gmem], stars:r.stars[acctimeInds.stars] }
      endif
      return,r
    endif
    
    ; otherwise, load maximum past temperature per gas particle
    ; (or statistics for the child population of each gas cell)
    maxt = maxTemps(sP=sP,/loadByGas)

    ; restrict maxPastTemp to subset of galaxycat (tracers)
    if keyword_set(maxPastTemp) then begin
      if keyword_set(trPopMax) then begin
        val_gal   = maxt.maxTemps_gal[galcatInds.gal]
        val_gmem  = maxt.maxTemps_gmem[galcatInds.gmem]
        val_stars = maxt.maxTemps_stars[galcatInds.stars]
      endif
      if keyword_set(trPopMin) then begin
        val_gal   = maxt.maxTemps_min_gal[galcatInds.gal]
        val_gmem  = maxt.maxTemps_min_gmem[galcatInds.gmem]
        val_stars = maxt.maxTemps_min_stars[galcatInds.stars]
      endif
      if keyword_set(trPopMean) then begin
        val_gal   = maxt.maxTemps_mean_gal[galcatInds.gal]
        val_gmem  = maxt.maxTemps_mean_gmem[galcatInds.gmem]
        val_stars = maxt.maxTemps_mean_stars[galcatInds.stars]
      endif
    endif
    
    ; restrict maxTempTime to subset of galaxycat (and convert to redshift)
    if keyword_set(maxTempTime) then begin
      if keyword_set(trPopMax) then begin
        val_gal   = 1/maxt.maxTempTime_gal[galcatInds.gal]-1
        val_gmem  = 1/maxt.maxTempTime_gmem[galcatInds.gmem]-1
        val_stars = 1/maxt.maxTempTime_stars[galcatInds.stars]-1
      endif
      if keyword_set(trPopMin) then begin
        val_gal   = 1/maxt.maxTempTime_min_gal[galcatInds.gal]-1
        val_gmem  = 1/maxt.maxTempTime_min_gmem[galcatInds.gmem]-1
        val_stars = 1/maxt.maxTempTime_min_stars[galcatInds.stars]-1
      endif
      if keyword_set(trPopMean) then begin
        val_gal   = 1/maxt.maxTempTime_mean_gal[galcatInds.gal]-1
        val_gmem  = 1/maxt.maxTempTime_mean_gmem[galcatInds.gmem]-1
        val_stars = 1/maxt.maxTempTime_mean_stars[galcatInds.stars]-1
      endif
    endif
    
    ; restrict temps/times to subset of galaxy cat (sph)
    if sP.trMCPerCell eq 0 then begin
      if keyword_set(maxPastTemp) then val_gal   = maxt.maxTemps_gal[galcatInds.gal]
      if keyword_set(maxPastTemp) then val_gmem  = maxt.maxTemps_gmem[galcatInds.gmem]
      if keyword_set(maxPastTemp) then val_stars = maxt.maxTemps_stars[galcatInds.stars]
      
      if keyword_set(maxTempTime) then val_gal   = 1/maxt.maxTempTime_gal[galcatInds.gal]-1
      if keyword_set(maxTempTime) then val_gmem  = 1/maxt.maxTempTime_gmem[galcatInds.gmem]-1
      if keyword_set(maxTempTime) then val_stars = 1/maxt.maxTempTime_stars[galcatInds.stars]-1
    endif
    
    maxt = !NULL
    
    r = {gal:val_gal,gmem:val_gmem,stars:val_stars}
  endif
  
  ; tracer expansion if requested
  if keyword_set(allTR) then begin
    ; load child counts for both galaxy members and group members
    rtr = maxTemps(sP=sP,/loadAllTRGal)
    gal_child_counts = rtr.child_counts[galcatInds.gal]
    
    rtr = maxTemps(sP=sP,/loadAllTRGmem)
    gmem_child_counts = rtr.child_counts[galcatInds.gmem]
    
    rtr = maxTemps(sP=sP,/loadAllTRStars)
    stars_child_counts = rtr.child_counts[galcatInds.stars]
    
    rtr = !NULL  
    
    ; create new return arrays
    if size(r.gal,/tname) eq 'LONG64' then $ ; to support curSingleVal=ids
      rr = { gal   : lon64arr(total(gal_child_counts,/int))  ,$
             gmem  : lon64arr(total(gmem_child_counts,/int)) ,$
             stars : lon64arr(total(stars_child_counts,/int)) }
             
    if size(r.gal,/tname) eq 'LONG' then $
      rr = { gal   : lonarr(total(gal_child_counts,/int))  ,$
             gmem  : lonarr(total(gmem_child_counts,/int)) ,$
             stars : lonarr(total(stars_child_counts,/int)) }      
             
    if size(r.gal,/tname) eq 'FLOAT' then $
      rr = { gal   : fltarr(total(gal_child_counts,/int))  ,$
             gmem  : fltarr(total(gmem_child_counts,/int)) ,$
             stars : fltarr(total(stars_child_counts,/int)) }
    
    if n_elements(rr) eq 0 then message,'Error: Unknown type for tracer replication.'
    
    if n_elements(gal_child_counts) ne n_elements(r.gal) then message,'Error1'
    if n_elements(gmem_child_counts) ne n_elements(r.gmem) then message,'Error2'
    if n_elements(stars_child_counts) ne n_elements(r.stars) then message,'Error3'
    
    ; replicate byGas arrays into byTracer arrays
    offset = 0L
    
    for i=0L,n_elements(r.gal)-1 do begin
      if gal_child_counts[i] gt 0 then begin
        rr.gal[offset:offset+gal_child_counts[i]-1] = replicate(r.gal[i],gal_child_counts[i])
        offset += gal_child_counts[i]
      endif
    endfor
    
    ; replicate byGmem arrays into byTracer arrays
    offset = 0L
    
    for i=0L,n_elements(r.gmem)-1 do begin
      if gmem_child_counts[i] gt 0 then begin
        rr.gmem[offset:offset+gmem_child_counts[i]-1] = replicate(r.gmem[i],gmem_child_counts[i])
        offset += gmem_child_counts[i]
      endif
    endfor
    
    ; replicate byStars arrays into byTracer arrays
    offset = 0L
    
    for i=0L,n_elements(r.stars)-1 do begin
      if stars_child_counts[i] gt 0 then begin
        rr.stars[offset:offset+stars_child_counts[i]-1] = replicate(r.stars[i],stars_child_counts[i])
        offset += stars_child_counts[i]
      endif
    endfor
    
    ; take accretionTime subset and return tracer expanded array
    if keyword_set(accretionTimeSubset) then begin
      rr = { gal:rr.gal[accTimeInds.gal], gmem:rr.gmem[accTimeInds.gmem], stars:rr.stars[acctimeInds.stars] }
    endif
    return,rr
  endif ; allTR
  
  return,r
end


