; galaxyCat.pro
; gas accretion project - gas selections of interest (galaxy/halo catalogs)
; dnelson mar.2013

; galaxyCat(): if snap not specified, create and save complete galaxy catalog from the group catalog by 
;              imposing additional cut in the (rho,temp) plane (same as that used by Torrey+ 2011)
;              if snap is specified, create only for one snapshot number or return previously saved
;              results for that snapshot
; Note: the gal/gmem catalogs have the same size (and indexing) as the subgroups
; galaxyOnly=1 : calculate on the fly and return for galaxy only (don't save)

function galaxyCat, sP=sP, galaxyOnly=galaxyOnly

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
    calcMatch,gcPIDs,ids,gc_ind,ids_ind,count=countMatch
    ids_ind = ids_ind[calcSort(gc_ind)] ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs
                                    ; need this if we want ids[ids_ind], temp[ids_ind], etc to be
                                    ; in the same order as the group catalog id list

    if countMatch ne n_elements(gcPIDs) then message,'Error: Failed to locate all gas gcPIDs in gas ids.'

    gcPIDs = !NULL
    ids = ids[ids_ind]

    ; load star ids and match to catalog
    gcPIDs_stars = gcPIDList(gc=gc,select='all',partType='stars')
    ids_stars = loadSnapshotSubset(sP=sP,partType='stars',field='ids')
    calcMatch,gcPIDs_stars,ids_stars,gc_ind,ids_ind_stars,count=countMatch
    ids_ind_stars = ids_ind_stars[calcSort(gc_ind)] ; IMPORTANT! rearrange ids_ind_stars

    if countMatch ne n_elements(gcPIDs_stars) then message,'Error: Failed to locate all star gcPIDs in star ids.'
    
    ; if GFM_WINDS, remove any wind particles from these stars
    if sP.gfmWinds ne 0 then begin
        btime_stars = loadSnapshotSubset(sP=sP,partType='stars',field='gfm_sftime')
        
        w_nonwind = where(btime_stars[ids_ind_stars] ge 0.0,count_nonwind,ncomp=count_wind)
        
        if count_nonwind eq 0 then message,'Error'
        print,'Excluding ['+str(count_wind)+'] wind particles ('+$
          string(float(count_wind)/countMatch,format='(f5.2)')+'%).'
        
        ; take subset
        ids_ind_stars = ids_ind_stars[w_nonwind]
    endif
    
    gcPIDs = !NULL
    gcPIDs_stars = !NULL
    gc_ind = !NULL
    ids_stars = ids_stars[ids_ind_stars]

    ; load u,nelec and calculate temp of gas
    u = loadSnapshotSubset(sP=sP,partType='gas',field='u',inds=ids_ind)
    nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec',inds=ids_ind)
    
    temp = convertUtoTemp(u,nelec)
    
    u     = !NULL
    nelec = !NULL
    
    ; load rho of gas and make galaxy (rho,temp) plane cut
    dens = loadSnapshotSubset(sP=sP,partType='gas',field='density',inds=ids_ind)
    
    ; scale Torrey+ (2011) galaxy cut to physical density
    scalefac = snapNumToRedshift(sP=sP,/time) ; time flag gives simulation time = scale factor
    a3inv = 1.0 / (scalefac*scalefac*scalefac)
    dens *= a3inv
    
    ; (rho,temp) cut
    wGal = where(alog10(temp) - sP.galcut_rho * alog10(dens) lt sP.galcut_T,$
                 countGal,comp=wGmem,ncomp=countGmem)
                 
    ; calculate radial distances of gas elements to primary parents
    if sP.radcut_rvir ne 0.0 then begin
      pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos',inds=ids_ind)
      
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
      
      ; if GFM_WINDS, restrict sgParIDs_stars to non-wind particles
      if sP.gfmWinds ne 0 then begin
        sgParIDs_stars = sgParIDs_stars[w_nonwind]
        offset.stars -= count_wind
      endif
      
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
      wGal  = where(alog10(temp) - sP.galcut_rho*alog10(dens) lt sP.galcut_T and $
                     rad_pri lt sP.radcut_rvir,countGal)
      wGmem = where(alog10(temp) - sP.galcut_rho*alog10(dens) ge sP.galcut_T and $
                     rad_pri ge sP.radcut_rvir and rad_pri le sP.radcut_out,countGmem)
                     
      ; load stellar positions calculate radial vector of stars from group center (PRI)
      pos = loadSnapshotSubset(sP=sP,partType='stars',field='pos',inds=ids_ind_stars)
      
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

    ; construct galaxy catalog
    if (not file_test(saveFilename1)) then begin
      galaxyLen = lonarr(gc.nSubgroupsTot)
      galaxyOff = lonarr(gc.nSubgroupsTot)
      galaxyIDs = lon64arr(n_elements(ids))
      
      nextOff = 0UL
      
      ; match indices between galaxy ids and group ids
      calcMatch,ids,gc.IDs,ids_ind,gcIDs_ind,count=countID
      
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
      calcMatch,ids,galaxyIDs,ind1,ind2,count=countCheck
      if countCheck ne n_elements(ids) then message,'Error: Failed to locate all gal.'

      ; immediate return for galaxy only (accretionRates)?
      if keyword_set(galaxyOnly) then begin
        r = {galaxyLen:galaxyLen, galaxyOff:galaxyOff, galaxyIDs:galaxyIDs}
        return, r
      endif
      
      ; save galaxy catalog
      save,galaxyLen,galaxyOff,galaxyIDs,filename=saveFilename1
      print,'Saved: '+strmid(saveFilename1,strlen(sP.derivPath))+' ['+str(countGal)+'/'+str(countMatch)+']'
    
      galaxyLen = !NULL
      galaxyOff = !NULL
      galaxyIDs = !NULL
      
    endif
    
    ; construct group member catalog
    if (not file_test(saveFilename2)) then begin
    
      groupmemLen = lonarr(gc.nSubgroupsTot)
      groupmemOff = lonarr(gc.nSubgroupsTot)
      groupmemIDs = lon64arr(n_elements(ids_groupmem))
      
      nextOff = 0UL
      
      ; match indices between gas ids and group member ids
      calcMatch,ids_groupmem,gc.IDs,ids_ind,gcIDs_ind,count=countID
      
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
      calcMatch,ids_groupmem,groupmemIDs,ind1,ind2,count=countCheck
      if countCheck ne n_elements(ids_groupmem) then message,'Error: Failed to locate all gmem.'
      
      ; save group membership catalog
      save,groupmemLen,groupmemOff,groupmemIDs,filename=saveFilename2    
      print,'Saved: '+strmid(saveFilename2,strlen(sP.derivPath))+' ['+str(countGmem)+'/'+str(countMatch)+']'    
    
      groupmemLen = !NULL
      groupmemOff = !NULL
      groupmemIDs = !NULL
    
    endif
    
    ; construct stellar catalog
    if (not file_test(saveFilename3)) then begin
      stellarLen = lonarr(gc.nSubgroupsTot)
      stellarOff = lonarr(gc.nSubgroupsTot)
      stellarIDs = lon64arr(n_elements(ids_stars))
      
      nextOff = 0UL
      
      ; match indices between star ids and group ids
      calcMatch,ids_stars,gc.IDs,ids_ind,gcIDs_ind,count=countID
      
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
      calcMatch,ids_stars,stellarIDs,ind1,ind2,count=countCheck
      if countCheck ne n_elements(ids_stars) then message,'Error: Failed to locate all stars.'

      ; save stellar catalog
      save,stellarLen,stellarOff,stellarIDs,filename=saveFilename3
      print,'Saved: '+strmid(saveFilename3,strlen(sP.derivPath))+' ['+str(countStars)+'/'+str(countMatch)+']'
    endif

  endfor ;snapRange
  
  ; if just one snap requested and just calculated, return it now
  if snapRange[0] eq snapRange[1] then begin
    restore,saveFilename1
    restore,saveFilename2
    restore,saveFilename3
    r = {galaxyOff:galaxyOff,galaxyLen:galaxyLen,galaxyIDs:galaxyIDs,$
         stellarOff:stellarOff,stellarLen:stellarLen,stellarIDs:stellarIDs,$
         groupmemOff:groupmemOff,groupmemLen:groupmemLen,groupmemIDs:groupmemIDs}
    return,r
  endif
  
end

; galaxyCatRadii(): find radial distance and radial velocities of all group member particles wrt 
;  the group they belong to as well as the rad to the primary group if this is a secondary group

function galaxyCatRadii, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; if snap specified, run only one snapshot (and/or just return previous results)
  if (sP.snap ne -1) then begin
    saveFilename = sP.derivPath + 'galradii.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.sav'
    
    ; results exist, return
    if (file_test(saveFilename)) then begin
      restore,saveFilename
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
    calcMatch,galcat.galaxyIDs,ids,galcat_ind,ids_gal_ind,count=countGal
    ids_gal_ind = ids_gal_ind[calcSort(galcat_ind)]
    
    calcMatch,galcat.groupmemIDs,ids,galcat_ind,ids_gmem_ind,count=countGmem
    ids_gmem_ind = ids_gmem_ind[calcSort(galcat_ind)]
    
    ids        = !NULL
    galcat_ind = !NULL
    
    ; calculate radial distances of gas elements to primary and secondary parents
    pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
    
    pos_gal  = pos[*,ids_gal_ind]
    pos_gmem = pos[*,ids_gmem_ind]
    
    pos = !NULL
    
    ; load stellar ids and match, calculate radial distances
    ids = loadSnapshotSubset(sP=sP,partType='stars',field='ids')
    
    calcMatch,galcat.stellarIDs,ids,galcat_ind,ids_stars_ind,count=countStars
    ids_stars_ind = ids_stars_ind[calcSort(galcat_ind)]
    
    ids        = !NULL
    galcat_ind = !NULL
    
    pos_stars = loadSnapshotSubset(sP=sP,partType='stars',field='pos',ind=ids_stars_ind)
    
    ; load subhalo catalog for mostBoundParticleID and for priParentIDs
    gc = loadGroupCat(sP=sP,/skipIDs)
    
    ; find group center positions with most bound particles for each group
    subgroupCen = subgroupPosByMostBoundID(sP=sP)

    ; replicate parent IDs (of SECONDARY/direct)
    gcIDs = galCatRepParentIDs(galcat=galcat)
    
    ; allocate save structure
    r = { gal_pri   : fltarr(n_elements(gcIDs.gal))   ,$
          gal_sec   : fltarr(n_elements(gcIDs.gal))   ,$
          gmem_pri  : fltarr(n_elements(gcIDs.gmem))  ,$
          gmem_sec  : fltarr(n_elements(gcIDs.gmem))  ,$
          stars_pri : fltarr(n_elements(gcIDs.stars)) ,$
          stars_sec : fltarr(n_elements(gcIDs.stars)) ,$
          gal_vrad_pri   : fltarr(n_elements(gcIDs.gal))  ,$
          gmem_vrad_pri  : fltarr(n_elements(gcIDs.gmem)) ,$
          stars_vrad_pri : fltarr(n_elements(gcIDs.stars)) }
    
    ; calulate radial vector of gas from group center (SEC)
    r.gal_sec   = periodicDists(subgroupCen[*,gcIDs.gal],pos_gal,sP=sP)
    r.gmem_sec  = periodicDists(subgroupCen[*,gcIDs.gmem],pos_gmem,sP=sP)
    r.stars_sec = periodicDists(subgroupCen[*,gcIDs.stars],pos_stars,sP=sP)
    
    ; replicate parent IDs (of PRIMARY/parent)
    priParentIDs = gcIDList(gc=gc,select='pri')

    gcIDs = galCatRepParentIDs(galcat=galcat,priParentIDs=priParentIDs)
    priParentIDs = !NULL
    galcat = !NULL
    
    ; calulate radial vector of gas from group center (PRI)
    r.gal_pri   = periodicDists(subgroupCen[*,gcIDs.gal],pos_gal,sP=sP)
    r.gmem_pri  = periodicDists(subgroupCen[*,gcIDs.gmem],pos_gmem,sP=sP)
    r.stars_pri = periodicDists(subgroupCen[*,gcIDs.stars],pos_stars,sP=sP)
    
    ; replace coordinates by relative coordinates (radial vectors)
    for i=0,2 do begin
      pos_rel = reform(pos_gal[i,*] - subgroupCen[i,gcIDs.gal])
      correctPeriodicDistVecs, pos_rel, sP=sP
      pos_gal[i,*] = pos_rel
      
      pos_rel = reform(pos_gmem[i,*] - subgroupCen[i,gcIDs.gmem])
      correctPeriodicDistVecs, pos_rel, sP=sP
      pos_gmem[i,*] = pos_rel
      
      pos_rel = reform(pos_stars[i,*] - subgroupCen[i,gcIDs.stars])
      correctPeriodicDistVecs, pos_rel, sP=sP
      pos_stars[i,*] = pos_rel
    endfor
    
    ; load velocities
    vel = loadSnapshotSubset(sP=sP,partType='gas',field='vel')
    vel_gal   = vel[*,ids_gal_ind]
    vel_gmem  = vel[*,ids_gmem_ind]
    vel_stars = vel[*,ids_stars_ind]
    vel = !NULL
    
    r.gal_vrad_pri = ((vel_gal[0,*] - gc.subgroupVel[0,gcIDs.gal]) * pos_gal[0,*] + $
                      (vel_gal[1,*] - gc.subgroupVel[1,gcIDs.gal]) * pos_gal[1,*] + $
                      (vel_gal[2,*] - gc.subgroupVel[2,gcIDs.gal]) * pos_gal[2,*]) $
                      / r.gal_pri
                  
    r.gmem_vrad_pri = ((vel_gmem[0,*] - gc.subgroupVel[0,gcIDs.gmem]) * pos_gmem[0,*] + $
                       (vel_gmem[1,*] - gc.subgroupVel[1,gcIDs.gmem]) * pos_gmem[1,*] + $
                       (vel_gmem[2,*] - gc.subgroupVel[2,gcIDs.gmem]) * pos_gmem[2,*]) $
                       / r.gmem_pri      

    r.stars_vrad_pri = ((vel_stars[0,*] - gc.subgroupVel[0,gcIDs.stars]) * pos_stars[0,*] + $
                        (vel_stars[1,*] - gc.subgroupVel[1,gcIDs.stars]) * pos_stars[1,*] + $
                        (vel_stars[2,*] - gc.subgroupVel[2,gcIDs.stars]) * pos_stars[2,*]) $
                        / r.stars_pri
    
    ; save radial distances (and group centers)
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))

  endfor
  
  if snapRange[0] eq snapRange[1] then return,r
  
end

; galCatRepParentIDs(): for the galaxy catalog, replicate the list of ordered parent IDs such that
;                       the return array is the same size as the number of gas particles with
;                       each element the id of its parent subgroup/galaxy group/groupmember group
;
; priParentIDs=1 : should be the primary parent group ID for each group,
;                  e.g. priParentIDs=gcIDList(gc=gc,select='pri')
;                  if set, for secondary groups, return instead of the ID of the primary parent
;                  (for primary groups the return is unchanged)
;
; gcIDList : return only a replicated parent ID list of the specified subgroups in the groupcat
; child_counts: return a replicated list for tracers which have multiplicity inside each gas cell
;               or star particle as specified in child_counts

function galCatRepParentIDs, galcat=galcat, priParentIDs=priParentIDs, $
                             gcIDList=gcIDList, child_counts=child_counts

  compile_opt idl2, hidden, strictarr, strictarrsubs

    if not keyword_set(priParentIDs) then $
      priParentIDs = lindgen(n_elements(galcat.galaxyLen)) ; valid id list set to all
    
    if not keyword_set(gcIDList) then $
      gcIDList = lindgen(n_elements(galcat.galaxyLen)) ; id list to process set to all
    
    if keyword_set(child_counts) then $
      if n_elements(child_counts.gal) ne total(galcat.galaxyLen[gcIDList],/int) or $
         n_elements(child_counts.gmem) ne total(galcat.groupmemLen[gcIDList],/int) or $
         n_elements(child_counts.stars) ne total(galcat.stellarLen[gcIDList],/int) then $
         message,'Error: Child_counts gal/gmem/stars should have same size as galcat gcIDList subset.'
         
    if keyword_set(child_counts) then begin
      if total(child_counts.gal,/int) gt 2e9 then stop ; consider lon64/removing /int
      if total(child_counts.gmem,/int) gt 2e9 then stop
      if total(child_counts.stars,/int) gt 2e9 then stop
    endif else begin
      child_counts = { gal   : lonarr(total(galcat.galaxyLen[gcIDList],/int))+1   ,$
                       gmem  : lonarr(total(galcat.groupmemLen[gcIDList],/int))+1 ,$
                       stars : lonarr(total(galcat.stellarLen[gcIDList],/int))+1   }
    endelse
    
    r = { gal   : lonarr(total(child_counts.gal,/int))  ,$
          gmem  : lonarr(total(child_counts.gmem,/int)) ,$
          stars : lonarr(total(child_counts.stars,/int)) }
    
    offset   = { gal : 0L, gmem: 0L, stars: 0L }
    offset_c = { gal : 0L, gmem: 0L, stars: 0L }
    
    foreach gcID,gcIDList do begin
    
        ; parent ID (pri or sec)
        if (total(gcID eq priParentIDs) eq 0) then begin
          dists = gcID - priParentIDs
          dists = dists[where(dists gt 0)]
          w = where(dists eq min(dists),count)
          parID = priParentIDs[w[0]]
        endif else begin
          parID = gcID
        endelse
        
        ; galaxies
        if galcat.galaxyLen[gcID] gt 0 then begin
          tot_children = total(child_counts.gal[offset_c.gal:offset_c.gal+galcat.galaxyLen[gcID]-1],/int)
          gal_ind_end = offset.gal+tot_children-1

          if tot_children gt 0 then r.gal[offset.gal:gal_ind_end] = cmreplicate(parID,tot_children)
          offset.gal += tot_children
          offset_c.gal += galcat.galaxyLen[gcID]
        endif
        
        ; group members
        if galcat.groupmemLen[gcID] gt 0 then begin
          tot_children = total(child_counts.gmem[offset_c.gmem:offset_c.gmem+galcat.groupmemLen[gcID]-1],/int)
          gmem_ind_end = offset.gmem+tot_children-1

          if tot_children gt 0 then r.gmem[offset.gmem:gmem_ind_end] = cmreplicate(parID,tot_children)
          offset.gmem += tot_children
          offset_c.gmem += galcat.groupmemLen[gcID]
        endif
        
        ; stars
        if galcat.stellarLen[gcID] gt 0 then begin
          tot_children = total(child_counts.stars[offset_c.stars:offset_c.stars+galcat.stellarLen[gcID]-1],/int)
          stars_ind_end = offset.stars+tot_children-1

          if tot_children gt 0 then r.stars[offset.stars:stars_ind_end] = cmreplicate(parID,tot_children)
          offset.stars += tot_children
          offset_c.stars += galcat.stellarLen[gcID]
        endif
        
    endforeach

    return,r
end

; galCatParentProperties: calculate some property of the parent galaxy/group for every gas elements
;                         in the galaxy catalog at some snapshot
; virTemp=1 : virial temperature
; mass=1    : total mass (from catalog, dm+baryon)
; rVir=1    : virial radius (r_200 critical)
; parNorm   : 'pri' or 'sec' (if pri then return properties of primary parent even for gas elements
;             in secondary/"satelitte" subgroups) (if sec effectively ignored, this is default behavior)

function galCatParentProperties, sP=sP, virTemp=virTemp, mass=mass, rVir=rVir, vCirc=vCirc, $
                                 parNorm=parNorm

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function galaxyCat, snapNumToRedshift, codeMassToLogMsun

  ; load group catalog for masses
  gc = loadGroupCat(sP=sP,/skipIDs)

  ; load galaxy catalog
  galcat = galaxyCat(sP=sP)

  ; replicate parent IDs
  if keyword_set(parNorm) then begin
    if parNorm eq 'pri' then begin
      ; replicate parent IDs (of PRIMARY)
      priParentIDs = gcIDList(gc=gc,select='pri')
      gcInd = galCatRepParentIDs(galcat=galcat,priParentIDs=priParentIDs)
      priParentIDs = !NULL
    endif
    if parNorm eq 'sec' then begin
      ; do the usual (use most direct, e.g. secondary, parent)
      gcInd = galCatRepParentIDs(galcat=galcat)
    endif
    if parNorm ne 'pri' and parNorm ne 'sec' then message,'Error! Unrecognized parNorm'
  endif else begin
    ; parNorm not specified, do the usual (use most direct, e.g. secondary, parent)
    gcInd = galCatRepParentIDs(galcat=galcat)
  endelse
   
  ; arrays
  r = { gal   : fltarr(n_elements(galcat.galaxyIDs))    ,$
        gmem  : fltarr(n_elements(galcat.groupmemIDs))  ,$
        stars : fltarr(n_elements(galcat.stellarIDs))    }
  
  ; masses (log msun)
  if keyword_set(mass) then begin
    r.gal   = gc.subgroupMass[gcInd.gal]
    r.gmem  = gc.subgroupMass[gcInd.gmem]
    r.stars = gc.subgroupMass[gcInd.stars]
    
    r.gal   = codeMassToLogMsun(r.gal)
    r.gmem  = codeMassToLogMsun(r.gmem)
    r.stars = codeMassToLogMsun(r.stars)
  endif

  ; calculate virial temperatures (K)
  if keyword_set(virTemp) then begin
    r.gal   = gc.subgroupMass[gcInd.gal]
    r.gmem  = gc.subgroupMass[gcInd.gmem]
    r.stars = gc.subgroupMass[gcInd.stars]
    
    redshift = snapNumToRedshift(sP=sP)
  
    r.gal   = alog10(codeMassToVirTemp(r.gal,sP=sP))
    r.gmem  = alog10(codeMassToVirTemp(r.gmem,sP=sP))
    r.stars = alog10(codeMassToVirTemp(r.stars,sP=sP))
  endif
  
  if keyword_set(rVir) then begin
    r.gal   = gc.subgroupGrnr[gcInd.gal]
    r.gmem  = gc.subgroupGrnr[gcInd.gmem]
    r.stars = gc.subgroupGrnr[gcInd.stars]
    
    r.gal   = gc.group_r_crit200[r.gal]
    r.gmem  = gc.group_r_crit200[r.gmem]
    r.stars = gc.group_r_crit200[r.stars]
  endif
  
  if keyword_set(vCirc) then begin
    units = getUnits()
    r.gal   = gc.subgroupGrnr[gcInd.gal]
    r.gmem  = gc.subgroupGrnr[gcInd.gmem]
    r.stars = gc.subgroupGrnr[gcInd.stars]
    
    r.gal   = sqrt(units.G * gc.group_m_crit200[r.gal] / gc.group_r_crit200[r.gal])
    r.gmem  = sqrt(units.G * gc.group_m_crit200[r.gmem] / gc.group_r_crit200[r.gmem])
    r.stars = sqrt(units.G * gc.group_m_crit200[r.stars] / gc.group_r_crit200[r.stars])
  endif

  return,r
end

; galcatINDList(): return a list of indices into the galaxy/groupmem/subgroup catalog for a subset of the
;                  members defined by the subgroup ID list gcIDList
;                  
; child_counts: return a replicated list for tracers which have multiplicity inside each gas cell
;               as specified in child_counts

function galcatINDList, sP=sP, galcat=galcat, gcIDList=gcIDList, child_counts=child_counts

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function galaxyCat

  ; load galaxy cat if necessary
  if not keyword_set(galcat) then begin
    if not keyword_set(sP) then begin
      print,'Error: Must specify galcat or sP.' & stop
    endif
    galcat = galaxyCat(sP=sP)
  endif
  
  if max(galcat.galaxyLen+galcat.galaxyOff) gt 2e9 then stop ; change to lon64arr
  if max(galcat.groupmemLen+galcat.groupmemOff) gt 2e9 then stop
  if max(galcat.stellarLen+galcat.stellarOff) gt 2e9 then stop
  
  ; make mask for requested indices
  gcIDMask = bytarr(n_elements(galcat.galaxyLen))
  if keyword_set(gcIDList) then gcIDMask[gcIDList] = 1B  
  if ~keyword_set(gcIDList) then begin
    gcIDMask[*] = 1B
    gcIDList = lindgen(n_elements(galcat.galaxyLen))
  endif
  
  ; normal indices return
  r = {gal   : ulonarr(total(galcat.galaxyLen[gcIDList],/int))   ,$
       gmem  : ulonarr(total(galcat.groupmemLen[gcIDList],/int)) ,$
       stars : ulonarr(total(galcat.stellarLen[gcIDList],/int))   }  
  
  offset = { gal : 0L, gmem : 0L, stars: 0L }
  
  ; (1) make list for gas cells/particles
  foreach gcID, gcIDList do begin
    ; galaxy
    if (galcat.galaxyLen[gcID] gt 0) then begin
      galInds    = ulindgen(galcat.galaxyLen[gcID]) + galcat.galaxyOff[gcID]
      r.gal[offset.gal:offset.gal+galcat.galaxyLen[gcID]-1] = galInds
      offset.gal += galcat.galaxyLen[gcID]
    endif
    
    ; group member
    if (galcat.groupmemLen[gcID] gt 0) then begin
      gmemInds    = ulindgen(galcat.groupmemLen[gcID]) + galcat.groupmemOff[gcID]
      r.gmem[offset.gmem:offset.gmem+galcat.groupmemLen[gcID]-1] = gmemInds
      offset.gmem += galcat.groupmemLen[gcID]
    endif
    
    ; stars
    if (galcat.stellarLen[gcID] gt 0) then begin
      starsInds    = ulindgen(galcat.stellarLen[gcID]) + galcat.stellarOff[gcID]
      r.stars[offset.stars:offset.stars+galcat.stellarLen[gcID]-1] = starsInds
      offset.stars += galcat.stellarLen[gcID]
    endif
  endforeach
  
  ; (2) make list including child counts if requested
  if ~keyword_set(child_counts) then return,r
  
  if n_elements(child_counts.gal) ne total(galcat.galaxyLen,/int) or $
     n_elements(child_counts.gmem) ne total(galcat.groupmemLen,/int) or $
     n_elements(child_counts.stars) ne total(galcat.stellarLen,/int) then $
     message,'Error: Child_counts gal/gmem/stars should have same size as full galcat subset.'

  if total(child_counts.gal,/int) gt 2e9 then stop ; consider lon64/removing /int
  if total(child_counts.gmem,/int) gt 2e9 then stop
  if total(child_counts.stars,/int) gt 2e9 then stop

  rcc = { gal   : ulonarr(total(child_counts.gal[r.gal],/int))    ,$
          gmem  : ulonarr(total(child_counts.gmem[r.gmem],/int))  ,$
          stars : ulonarr(total(child_counts.stars[r.stars],/int)) }
       
  offset     = { gal : 0L, gmem : 0L, stars : 0L }
  offset_all = { gal : 0L, gmem : 0L, stars : 0L }
  
  for gcID=0UL,n_elements(galcat.galaxyLen)-1 do begin
    ; galaxy
    if galcat.galaxyLen[gcID] gt 0 then begin
      ; calculate total number of children in this subgroup
      tot_children_gal  = total(child_counts.gal[galcat.galaxyOff[gcID]:$
                                                 galcat.galaxyOff[gcID]+galcat.galaxyLen[gcID]-1],/int)
      ; add indices only for specified galaxy IDs
      if gcIDMask[gcID] eq 1B and tot_children_gal gt 0 then begin
        ; calculate place and store indices
        galInds = ulindgen(tot_children_gal) + offset_all.gal
        rcc.gal[offset.gal:offset.gal+tot_children_gal-1] = galInds
  
        offset.gal += tot_children_gal
      endif
      
      ; add child counts to indice array offset
      offset_all.gal  += tot_children_gal
    endif      
                  
    ; group member
    if galcat.groupmemLen[gcID] gt 0 then begin
      ; calculate total number of children in this subgroup             
      tot_children_gmem = total(child_counts.gmem[galcat.groupmemOff[gcID]:$
                                                  galcat.groupmemOff[gcID]+galcat.groupmemLen[gcID]-1],/int)
            
      ; add indices only for specified group member IDs
      if gcIDMask[gcID] eq 1B and tot_children_gmem gt 0 then begin
        ; calculate place and store indices
        gmemInds     = ulindgen(tot_children_gmem) + offset_all.gmem
        rcc.gmem[offset.gmem:offset.gmem+tot_children_gmem-1] = gmemInds
        
        offset.gmem += tot_children_gmem
      endif
      
      offset_all.gmem += tot_children_gmem
    endif
    
    ; stars
    if galcat.stellarLen[gcID] gt 0 then begin
      ; calculate total number of children in this subgroup             
      tot_children_stars = total(child_counts.stars[galcat.stellarOff[gcID]:$
                                                   galcat.stellarOff[gcID]+galcat.stellarLen[gcID]-1],/int)
            
      ; add indices only for specified stellar IDs
      if gcIDMask[gcID] eq 1B and tot_children_stars gt 0 then begin
        ; calculate place and store indices
        starInds     = ulindgen(tot_children_stars) + offset_all.stars
        rcc.stars[offset.stars:offset.stars+tot_children_stars-1] = starInds
        
        offset.stars  += tot_children_stars
      endif
      
      offset_all.stars += tot_children_stars
    endif
    
  endfor
  
  if offset.gal   ne n_elements(rcc.gal) or $
     offset.gmem  ne n_elements(rcc.gmem) or $
     offset.stars ne n_elements(rcc.stars) then message,'Error.'

  return,rcc
  
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
    calcMatch,gcIDList,mt.galcatIDList,ind1,ind2,count=count
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
      tr_ids    = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
      tr_parids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentids')
      
      ids_gal   = tr_ids[cosmoTracerChildren(sP=sP,gasIDs=galcat.galaxyIDs[galcatInds.gal],tr_parids=tr_parids,/getInds)]
      ids_gmem  = tr_ids[cosmoTracerChildren(sP=sP,gasIDs=galcat.groupmemIDs[galcatInds.gmem],tr_parids=tr_parids,/getInds)]
      ids_stars = tr_ids[cosmoTracerChildren(sP=sP,starIDs=galcat.stellarIDs[galcatInds.stars],tr_parids=tr_parids,/getInds)]
      
      tr_ids = !NULL
      tr_parids = !NULL
      
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


