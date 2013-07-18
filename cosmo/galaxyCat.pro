; galaxyCat.pro
; gas accretion project - gas selections of interest (galaxy/halo catalogs)
; dnelson jun.2013

; galaxyCat(): if snap not specified, create and save complete galaxy catalog from the group catalog by 
;              imposing additional cut in the (rho,temp) plane (same as that used by Torrey+ 2011)
;              if snap is specified, create only for one snapshot number or return previously saved
;              results for that snapshot
; Note: the galaxyCat catalogs have the same size (and indexing) as the subgroups
; skipSave=1 : calculate on the fly and return IDs (don't save)

function galaxyCat, sP=sP, skipSave=skipSave

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; if snap specified, run only one snapshot (and/or just return previous results)
  if (sP.snap ne -1) then begin
    saveFilename = sP.derivPath + 'galaxyCat.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) +'.sav'
    
    ; results exist, return
    if file_test(saveFilename) then begin
      restore,saveFilename
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
    saveFilename = sP.derivPath + 'galaxyCat.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.sav'
        
    ; load ids of particles in all subfind groups
    gc = loadGroupCat(sP=sP,/readIDs)
    gcPIDs = gcPIDList(gc=gc,select='all',partType='gas')

    ; load gas ids and match to catalog
    ; --------
    ids_gas = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    calcMatch,gcPIDs,ids_gas,gc_ind,ids_ind,count=countMatch
    ids_ind = ids_ind[calcSort(gc_ind)] ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs

    if countMatch ne n_elements(gcPIDs) then message,'Error: Failed to locate all gas gcPIDs in gas ids.'

    gcPIDs = !NULL
    ids_gas = ids_gas[ids_ind]
    
    ; load star ids and match to catalog
    ; ---------
    gcPIDs_stars = gcPIDList(gc=gc,select='all',partType='stars')
    ids_stars = loadSnapshotSubset(sP=sP,partType='stars',field='ids')
    calcMatch,gcPIDs_stars,ids_stars,gc_ind,ids_ind_stars,count=countMatch
    ids_ind_stars = ids_ind_stars[calcSort(gc_ind)] ; IMPORTANT! rearrange ids_ind_stars

    if countMatch ne n_elements(gcPIDs_stars) then message,'Error: Failed to locate all star gcPIDs in star ids.'
    
    ; if GFM_WINDS, remove any wind particles from these stars
    if sP.gfmWinds ne 0 then begin
        btime_stars = loadSnapshotSubset(sP=sP,partType='stars',field='gfm_sftime')
        
        w_nonwind = where(btime_stars[ids_ind_stars] ge 0.0,count_nonwind,ncomp=count_wind,comp=w_wind)
        
        if count_nonwind eq 0 then message,'Error'
        print,'Excluding ['+str(count_wind)+'] wind particles ('+$
          string(float(count_wind)/countMatch,format='(f5.2)')+'%) from galaxy.'
        
        ; take subset of stellar ids which are in the wind (for inter)
        ids_stars_wind = ids_stars[ids_ind_stars[w_wind]]
        
        ; modify stellar indices to keep considering for gal(stars) membership
        ids_ind_stars = ids_ind_stars[w_nonwind]
    endif
    
    gcPIDs = !NULL
    gcPIDs_stars = !NULL
    ids_stars = ids_stars[ids_ind_stars]
    
    ; load BHs ids and match to catalog
    ; --------
    if sP.gfmBHs ne 0 then begin
      gcPIDs_BHs = gcPIDList(gc=gc,select='all',partType='BHs')
      ids_BHs = loadSnapshotSubset(sP=sP,partType='BHs',field='ids')
      calcMatch,gcPIDs_BHs,ids_BHs,gc_ind,ids_ind_BHs,count=countMatch
      ids_ind_BHs = ids_ind_BHs[calcSort(gc_ind)] ; IMPORTANT! rearrange ids_ind_BHs

      if countMatch ne n_elements(gcPIDs_BHs) then message,'Error: Failed to locate all BH gcPIDs in star ids.'
      
      gcPIDs_BHs = !NULL
      ids_BHs = ids_BHs[ids_ind_BHs]
    endif
    
    ; load u,nelec and calculate temp of gas
    ; ------------
    u     = loadSnapshotSubset(sP=sP,partType='gas',field='u',inds=ids_ind)
    nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec',inds=ids_ind)
    
    temp  = convertUtoTemp(u,nelec,/log)
    u     = !NULL
    nelec = !NULL
    
    ; load rho of gas and make galaxy (rho,temp) plane cut
    dens = loadSnapshotSubset(sP=sP,partType='gas',field='density',inds=ids_ind)
    
    ; scale Torrey+ (2011) galaxy cut to physical density
    scalefac = snapNumToRedshift(sP=sP,/time) ; time flag gives simulation time = scale factor
    a3inv = 1.0 / (scalefac*scalefac*scalefac)
    dens = alog10( dens * a3inv )
    
    ; calculate radial distances of gas elements to primary parents
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
    ptNum = { gas : partTypeNum('gas'), stars : partTypeNum('stars'), BHs : partTypeNum('BHs') }
      
    sgParIDs_gas   = lonarr(total(gc.subgroupLenType[ptNum.gas,*],/int))
    sgParIDs_stars = lonarr(total(gc.subgroupLenType[ptNum.stars,*],/int))
    
    if sP.gfmBHs ne 0 then $
      sgParIDs_BHs = lonarr(total(gc.subgroupLenType[ptNum.BHs,*],/int))
    
    offset = { gas : 0L, stars: 0L, BHs: 0L }
      
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
      
      if gc.subgroupLenType[ptNum.BHs,sgID] gt 0 then begin
        sgParIDs_BHs[offset.BHs:offset.BHs+gc.subgroupLenType[ptNum.BHs,sgID]-1] = $
          cmreplicate(parIDs[sgID],gc.subgroupLenType[ptNum.BHs,sgID])
        offset.BHs += gc.subgroupLenType[ptNum.BHs,sgID]
      endif
    endfor
      
    ; if GFM_WINDS, restrict sgParIDs_stars to non-wind particles
    if sP.gfmWinds ne 0 then begin
      sgParIDs_stars = sgParIDs_stars[w_nonwind]
      offset.stars -= count_wind
    endif
      
    if offset.gas   ne n_elements(sgParIDs_gas)   then message,'Error: Bad parent ID replication.'
    if offset.stars ne n_elements(sgParIDs_stars) then message,'Error: Bad parent ID star rep.'
    if offset.BHs   ne n_elements(sgParIDs_BHs)   then message,'Error: Bad parent ID BH rep.'
      
    ; fof has rvir=0 (no SO values) if zero subgroups, marginal overdensity, or low total mass
    ; these are low mass halos which we aren't going to plot anyways
    w = where(gc.group_r_crit200 eq 0.0,count)
    if count gt 0 then gc.group_r_crit200[w] = 1e-8 ; remove with outer radial cut
      
    ; calculate radial vector of gas from group center (PRI)
    rad_pri  = periodicDists(subgroupCen[*,sgParIDs_gas],pos,sP=sP)
    par_rvir = gc.group_r_crit200[gc.subgroupGrNr[sgParIDs_gas]] ; normalize by fof parent rvir

    rad_pri /= par_rvir

    ; make a mask for "inter": all gas failing both the gal and gmem cuts, and all stars
    ; failing the star cut, plus any current wind particles
    interMask_gas   = intarr(n_elements(temp))
    interMask_stars = intarr(n_elements(sgParIDs_stars))
    
    ; (rho,temp,rad) cut
    wGal  = where(temp - sP.galcut_rho*dens lt sP.galcut_T and rad_pri lt sP.radcut_rvir,countGal)
    
    wGmem = where(temp - sP.galcut_rho*dens ge sP.galcut_T and $
                   rad_pri ge sP.radcut_rvir and rad_pri le sP.radcut_out,countGmem)
                   
    wInter = where( ((temp - sP.galcut_rho*dens gt sP.galcut_T and rad_pri lt sP.radcut_rvir) or $ ; fail gal
                     (temp - sP.galcut_rho*dens lt sP.galcut_T and rad_pri gt sP.radcut_rvir)) and $ ; fail gmem
                     rad_pri le sP.radcut_out, countInter) ; keep outer rad cut
    temp = !NULL
    dens = !NULL
             
    ; load stellar positions, calculate radial vector from group center (PRI)
    pos = loadSnapshotSubset(sP=sP,partType='stars',field='pos',inds=ids_ind_stars)
     
    rad_pri  = periodicDists(subgroupCen[*,sgParIDs_stars],pos,sP=sP)
    par_rvir = gc.group_r_crit200[gc.subgroupGrNr[sgParIDs_stars]] ; normalize by fof parent rvir
      
    rad_pri /= par_rvir
    
    ; make stellar cut
    wStars = where(rad_pri lt sP.radcut_rvir,countStars,$
                   comp=wStarsComp,ncomp=wStarsComp_num)

    if countGal eq 0 or countGmem eq 0 then begin
      print,'Warning: Empty galaxy cut or comp. Skipping: ' + strmid(saveFilename1,strlen(sP.derivPath))
      continue
    endif
        
    ; load BH positions, calculate radial vector from group center (PRI)
    countBHs = 0

    if sP.gfmBHs ne 0 then begin
      pos = loadSnapshotSubset(sP=sP,partType='BHs',field='pos',inds=ids_ind_bhs)
     
      rad_pri  = periodicDists(subgroupCen[*,sgParIDs_BHs],pos,sP=sP)
      par_rvir = gc.group_r_crit200[gc.subgroupGrNr[sgParIDs_BHs]] ; normalize by fof parent rvir
      
      rad_pri /= par_rvir
      
      wBHs = where(rad_pri lt sP.radcut_rvir,countBHs)
    endif

    ; "inter" is all non-galaxy, non-gmem gas in the subgroup, plus any current wind
    ids_inter = []
    
    if countInter gt 0 then ids_inter = [ ids_inter,ids_gas[wInter] ]
    
    if wStarsComp_num gt 0 then ids_inter = [ ids_inter,ids_stars[wStarsComp] ]
    
    if sP.gfmWinds ne 0 then $
      if count_wind gt 0 then ids_inter = [ ids_inter,ids_stars_wind ]
    
    countInter = n_elements(ids_inter)
    
    ; make subsets of ids matching galaxy cut and complement
    countTot = countGal + countGmem + countStars + countInter + countBHs

    types = {gal:1, gmem:2, inter:3, stars:4}
    if countBHs gt 0 then types = mod_struct( types, 'bhs', 5)
    
    ; construct galaxy catalog
    r = { len   : lonarr(gc.nSubgroupsTot) ,$
          off   : lonarr(gc.nSubgroupsTot) ,$
          ids   : lon64arr(countTot)       ,$
          type  : intarr(countTot)         ,$
          types : types                    ,$
          countGal: countGal, countGmem: countGmem, countStars: countStars ,$
          countInter: countInter, countBHs: countBHs, countTot:countTot,$
          nGroups: gc.nSubgroupsTot }
    
    ; insert IDs and types
    nextOff = 0L
    
    if countGal gt 0 then begin
      ids_gal = ids_gas[wGal]
      r.ids[nextOff : nextOff+countGal-1] = ids_gal
      r.type[nextOff : nextOff+countGal-1] = 1
      nextOff += countGal
    endif

    if countGmem gt 0 then begin
      ids_groupmem = ids_gas[wGmem]
      r.ids[nextOff : nextOff+countGmem-1] = ids_groupmem
      r.type[nextOff : nextOff+countGmem-1] = 2
      nextOff += countGmem
    endif
    
    if countInter gt 0 then begin
      r.ids[nextOff : nextOff+countInter-1] = ids_inter
      r.type[nextOff : nextOff+countInter-1] = 3
      nextOff += countInter
    endif
    
    if countStars gt 0 then begin
      ids_stars = ids_stars[wStars]
      r.ids[nextOff : nextOff+countStars-1] = ids_stars
      r.type[nextOff : nextOff+countStars-1] = 4
      nextOff += countStars
    endif
    
    if countBHs gt 0 then begin
      ids_bh = ids_BHs[wBHs]
      r.ids[nextOff : nextOff+countBHs-1] = ids_bh
      r.type[nextOff : nextOff+countBHs-1]   = 5
      nextOff += countBHs
    endif
    
    if nuniq(r.ids) ne n_elements(r.ids) then message,'Error: r.ids is not unique.'

    calcMatch,r.ids,gc.IDs,r_ind,gc_ind,count=countCheck
    gc_ind = gc_ind[sort(r_ind)]
    if countCheck ne n_elements(r.ids) then message,'Error: Count mismatch.'
    
    ; calculate which parent subgroup each belongs to, and rearrange in order of SG parent
    par_inds = value_locate(gc.subgroupOffset,gc_ind)
      
    if min(par_inds) lt 0 or max(par_inds) ge n_elements(gc.subgroupOffset) then message,'Error'
    
    print,'TODO: Change this to a stable sort (preserve order within each parent)'
    par_inds_sort = sort(par_inds)
    
    r.ids  = r.ids[par_inds_sort]
    r.type = r.type[par_inds_sort]
    
    ; calculate the galaxyCat length and offset for each subgroup
    r.len = histogram(par_inds,min=0,max=n_elements(gc.subgroupOffset)-1)
    r.off = [0,(total(r.len,/int,/cum))[0:-2]]

    ; check that len of all secondary groups is zero
    gcIDs_sec = gcIDList(gc=gc,select='sec')
    if total(r.len[gcIDs_sec]) gt 0 then message,'Error'
    
    ; calculate number of tracer children of each parent (trMC)
    if sP.trMCPerCell gt 0 then begin
      galcat_trids = cosmoTracerChildren(sP=sP, /getInds, gasIDs=r.ids, child_counts=galcat_cc)
      r = mod_struct( r, 'ccMC', galcat_cc )
    endif
    
    stop
    
    ; immediate return and skip save?
    if keyword_set(skipSave) then return, r
      
    ; save galaxy catalog
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))

  endfor ;snapRange
  
  ; if just one snap requested and just calculated, return it now
  if snapRange[0] eq snapRange[1] then return, r
  
end

; galaxyCatRadii(): find radial distance and radial velocities of all group member particles wrt 
;  the group they belong to as well as the rad to the primary group if this is a secondary group

function galaxyCatRadii, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; if snap specified, run only one snapshot (and/or just return previous results)
  if (sP.snap ne -1) then begin
    saveFilename = sP.derivPath + 'galaxyCatRadii.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.sav'
    
    ; results exist, return
    if file_test(saveFilename) then begin
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
  
    ; set save filename
    sP.snap  = m
    saveFilename = sP.derivPath + 'galaxyCatRadii.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.sav'
    
    ; load galaxy and group membership catalogs
    h = loadSnapshotHeader(sP=sP)
    galcat = galaxyCat(sP=sP)
    
    ; load subhalo catalog for mostBoundParticleID and for priParentIDs
    gc = loadGroupCat(sP=sP,/skipIDs)
    
    ; find group center positions with most bound particles for each group
    subgroupCen = subgroupPosByMostBoundID(sP=sP)

    ; replicate parent IDs (of SECONDARY/direct and PRIMARY/parent)
    gcIDs_sec = galCatRepParentIDs(galcat=galcat)
    
    priParentIDs = gcIDList(gc=gc,select='pri')
    gcIDs_pri = galCatRepParentIDs(galcat=galcat,priParentIDs=priParentIDs)
    priParentIDs = !NULL
      
    r = { rad_pri   : fltarr(galcat.countTot)  ,$
          rad_sec   : fltarr(galcat.countTot)  ,$
          vrad_pri  : fltarr(galcat.countTot)  ,$
          vrad_sec  : fltarr(galcat.countTot)   }
    
    ; restrict gas particle positions to gal/gmem gas only
    parTypes = ['gas','stars','BHs']
    
    foreach parType,parTypes do begin
    
      ; skip non-existent particle types
      if h.nPartTot[partTypeNum(parType)] eq 0 then continue

      ids_type = loadSnapshotSubset(sP=sP,partType=parType,field='ids')
    
      ; IMPORTANT! rearrange ids_ind to be in the order of gc.xIDs, need this if we want ids[ids_ind], 
      ; temp[ids_ind], etc to be in the same order as the galaxy catalog id list    
      calcMatch,galcat.ids,ids_type,galcat_ind,ids_type_ind,count=countType
      ids_type_ind = ids_type_ind[calcSort(galcat_ind)]
      
      if countType eq 0 then print,'Warning: Skipping parent type '+parType
    
      ids_type = !NULL
      galcat_ind = !NULL
    
      ; calculate radial distances to primary and secondary parents
      pos = loadSnapshotSubset(sP=sP,partType=parType,field='pos')
      pos = pos[*,ids_type_ind]
            
      ; calulate radial vector of gas from group center (PRI and SEC)
      r.rad_sec   = periodicDists(subgroupCen[*,gcIDs_sec],pos,sP=sP)
      r.rad_pri   = periodicDists(subgroupCen[*,gcIDs_pri],pos,sP=sP)
    
      ; replace coordinates by relative coordinates (radial vectors) to direct parent
      for i=0,2 do begin
        pos_rel = reform(pos[i,*] - subgroupCen[i,gcIDs_pri])
        correctPeriodicDistVecs, pos_rel, sP=sP
        pos[i,*] = pos_rel
      endfor
    
      ; load velocities
      vel = loadSnapshotSubset(sP=sP,partType=parType,field='vel')
      vel = vel[*,ids_type_ind]
    
      r.vrad_pri = ((vel[0,*] - gc.subgroupVel[0,gcIDs_pri]) * pos[0,*] + $
                    (vel[1,*] - gc.subgroupVel[1,gcIDs_pri]) * pos[1,*] + $
                    (vel[2,*] - gc.subgroupVel[2,gcIDs_pri]) * pos[2,*]) $
                    / r.rad_pri
    
    endforeach
    
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

    if not keyword_set(priParentIDs) then priParentIDs = lindgen(galcat.nGroups) ; valid id list set to all
    
    if not keyword_set(gcIDList) then gcIDList = lindgen(galcat.nGroups) ; id list to process set to all
    
    if keyword_set(child_counts) then begin
      if n_elements(child_counts) ne total(galcat.len[gcIDList],/int) then $
         message,'Error: Child_counts should have same size as galcat gcIDList subset.'
         
      if total(child_counts,/int) gt 2e9 then stop ; consider lon64/removing /int
    endif else begin
      child_counts = lonarr(total(galcat.len[gcIDList],/int))+1
    endelse
    
    r = lonarr(total(child_counts,/int))
    
    offset   = 0L
    offset_c = 0L
    
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
        
        ; replicate for children
        if galcat.len[gcID] gt 0 then begin
          tot_children = total(child_counts[offset_c:offset_c+galcat.len[gcID]-1],/int)
          ind_end = offset + tot_children - 1

          if tot_children gt 0 then r[offset:ind_end] = cmreplicate(parID,tot_children)
          offset += tot_children
          offset_c += galcat.len[gcID]
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
  r = fltarr(galcat.nGroups)
  
  ; masses (log msun)
  if keyword_set(mass) then begin
    r = gc.subgroupMass[gcInd]
    r = codeMassToLogMsun(r)
  endif

  ; calculate virial temperatures (K)
  if keyword_set(virTemp) then begin
    r = gc.subgroupMass[gcInd]
    r = alog10(codeMassToVirTemp(r,sP=sP))
  endif
  
  if keyword_set(rVir) then begin
    r = gc.subgroupGrnr[gcInd]
    r = gc.group_r_crit200[r]
  endif
  
  if keyword_set(vCirc) then begin
    units = getUnits()
    
    r = gc.subgroupGrnr[gcInd]
    r = sqrt(units.G * gc.group_m_crit200[r] / gc.group_r_crit200[r])
  endif

  return,r
end

; galcatINDList(): return a list of indices into the galaxyCat for a subset of the
;                  members defined by the subgroup ID list gcIDList
;                  
; trRep: return a replicated list for tracers which have some multiplicity inside each gas cell
; gcSep: separate out indices by galaxyCat.type

function galcatINDList, sP=sP, galcat=galcat, gcIDList=gcIDList, gcSep=gcSep, trRep=trRep

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function galaxyCat

  ; load galaxy cat if necessary
  if not keyword_set(galcat) then begin
    if not keyword_set(sP) then begin
      print,'Error: Must specify galcat or sP.' & stop
    endif
    galcat = galaxyCat(sP=sP)
  endif
  
  if keyword_set(trRep) and sP.trMCPerCell le 0 then message,'Error.'
  
  if max(galcat.len+galcat.off) gt 2e9 then stop ; change to lon64arr
  
  ; make mask for requested indices
  gcIDMask = bytarr(galcat.nGroups)
  if keyword_set(gcIDList) then gcIDMask[gcIDList] = 1B  
  if ~keyword_set(gcIDList) then begin
    gcIDMask[*] = 1B
    gcIDList = lindgen(galcat.nGroups)
  endif
  
  ; normal indices return
  r = ulonarr(total(galcat.len[gcIDList],/int))
  
  offset = 0L
  
  ; (1) make list for gas cells/particles
  foreach gcID, gcIDList do begin
  
    if (galcat.len[gcID] gt 0) then begin
      galInds = ulindgen(galcat.len[gcID]) + galcat.off[gcID]
      r[offset:offset+galcat.len[gcID]-1] = galInds
      offset += galcat.len[gcID]
    endif
    
  endforeach
  
  if ~keyword_set(trRep) then begin
    if ~keyword_set(gcSep) then return,r
    
    ; split into gal,gmem,stars,inter,bhs
    type = galcat.type[r]
    
    rr = {}
    totCount = 0L
    
    for i=0,n_tags(galcat.types)-1 do begin
      w_type = where(type eq galcat.types.(i),count)
      if count gt 0 then $
        rr = mod_struct(rr, (tag_names(galcat.types))[i], r[w_type]) $
      else $
        rr = mod_struct(rr, (tag_names(galcat.types))[i], -1)
      totCount += count
    endfor
    
    if totCount ne n_elements(r) then message,'Error in gcSep.'
    
    return,rr
  endif
  
  ; (2) make list including child counts if requested
  if sP.trMCPerCell gt 0 then child_counts = galcat.ccMC
  if sP.trMCPerCell lt 0 then message,'todo'

  if total(child_counts,/int) gt 2e9 then stop ; consider lon64/removing /int

  rcc = ulonarr(total(child_counts[r],/int))
       
  offset     = 0L
  offset_all = 0L
  
  for gcID=0L,galcat.nGroups-1 do begin

    if galcat.len[gcID] gt 0 then begin
      ; calculate total number of children in this subgroup
      tot_children = total(child_counts[galcat.off[gcID]:galcat.off[gcID]+galcat.len[gcID]-1],/int)
      
      ; add indices only for specified galaxy IDs
      if gcIDMask[gcID] eq 1B and tot_children gt 0 then begin
      
        ; calculate place and store indices
        gInds = ulindgen(tot_children) + offset_all
        rcc[ offset:offset+tot_children-1 ] = gInds
  
        offset += tot_children
      endif
      
      ; add child counts to indice array offset
      offset_all  += tot_children
    endif
    
  endfor
  
  if offset ne n_elements(rcc) then message,'Error.'
  
  if ~keyword_set(gcSep) then return,rcc
  
  ; split into gal,gmem,stars,inter,bhs
  if sP.trMCPerCell lt 0 then message,'Error: trMC hard coded.'
  xx = replicate_var(galcat.ccMC, subset_inds=r)
  type = galcat.type[ xx.parent_inds ]
  
  rr = {}
  totCount = 0L
    
  for i=0,n_tags(galcat.types)-1 do begin
    w_type = where(type eq galcat.types.(i),count)
    if count gt 0 then $
      rr = mod_struct(rr, (tag_names(galcat.types))[i], rcc[w_type]) $
    else $
      rr = mod_struct(rr, (tag_names(galcat.types))[i], -1)
    totCount += count
  endfor
    
  if totCount ne n_elements(rcc) then message,'Error in gcSep.'

  return,rr

end

; gcSubsetProp(): read galaxy catalog for a specific subgroup selection (pri,sec,all) and
;                 return properties for each gas element/tracer (may or may not depend on parent halo)
;                 at the redshift specified by sP.snap (may or may not depend on previous time)
;
; select=''     : 'pri', 'sec' or 'all', the type of subgroups to select
; gcIDList=[]   : alternatively, an explicit list of subgroups to select
;
; rVirNorm=1    : radial distances normalized by r_vir of either primary or secondary parent
;  parNorm      : either 'pri' or 'sec' if rVirNorm requested
; virTemp=1     : current virial temperatures of parent halos
; parMass=1     : total mass (dm+baryonic) of parent halos (from catalog)
; curTemp=1     : current temperature of each element
;
; maxPastTemp=1 : maximum past previous temperature of each element
; maxTempTime=1 : time when maximum past previous temperature was reached (in redshift)
; maxPastEnt=1  : maximum past entropy
; maxPastDens=1 : maximum past density
; maxPastMach=1 : maximum past mach number
;  trPopMin,trPopMean,trPopMax : return the respective statistic for each gas cell for the tracers
;
; elemIDs=1       : ids of each element (either SPH particles or tracers)
;
; tracksFluid=1   : temp,dens,ent history of each tracer (trMC only)
;
; curSingleVal=1  : current single quantity (e.g. mass, density) returned without manipulation (gas, or replicated)
;  singleValField : field name in snapshot file for the above
; curTracerVal=1    : as above, but individual value for each tracer (e.g. wind_counter)
; singleTracerField : ...
;
; mergerTreeSubset    : return values only for the subset of halos tracked in the merger tree subset
;  jun2013: this option removed, since the mtS is exactly all primary halos, and this function only works
;           with exactly that selection, in the case that we want the atS
; accretionTimeSubset : return values only for the subset of particles/tracers with recorded accretion times
;  accTime,accTvir : time of accretion (in redshift) or virial temp of parent halo at time of accretion
;  accMode : return values only for one accretionMode (all,smooth,bclumpy,sclumpy,smooth)

function gcSubsetProp, sP=sP, select=select, gcIDList=gcIDList, $
           rVirNorm=rVirNorm, virTemp=virTemp, parMass=parMass, $
           curTemp=curTemp, maxPastTemp=maxPastTemp, maxTempTime=maxTempTime, $
           maxPastEnt=maxPastEnt, maxPastDens=maxPastDens, maxPastMach=maxPastMach, $
           trPopMin=trPopMin, trPopMax=trPopMax, trPopMean=trPopMean, $ ; for maxPastXX/maxXXTime only
           elemIDs=elemIDs, tracksFluid=tracksFluid, $
           curSingleVal=curSingleVal, singleValField=singleValField, $
           curTracerVal=curTracerVal, singleTracerField=singleTracerField, $ ; trMC only
           parNorm=parNorm, $ ; for rVirNorm,virTemp,parMass only
           accretionTimeSubset=accretionTimeSubset, $
           accTime=accTime,accTvir=accTvir,accMode=accMode ; for accretionTimeSubset only

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; check combinations of input options for validity
  if keyword_set(accretionTimeSubset) then if select ne 'pri' then $
    message,'Error: The atS (mtS) contains only pri subgroups (select must be pri).'
  if (keyword_set(accTime) or keyword_set(accTvir)) and ~keyword_set(accretionTimeSubset) then $
    message,'Error: Can only return accretion time or Tvir at accretion time for atS.'
  if keyword_set(accMode) and ~keyword_set(accretionTimeSubset) then $
    message,'Error: Can only return accretion mode subsets of the accretionTime subset.'
  if keyword_set(elemIDs) and (keyword_set(trPopMin) or keyword_set(trPopMean) or keyword_set(trPopMax)) then $
    message,'Error: Cannot return pop stats of unique element IDs.'
  if (~keyword_set(select) and ~keyword_set(gcIDList)) or (keyword_set(select) and keyword_set(gcIDList)) then $
    message,'Error: Should specify either group type OR explicit list of groups.'
  if keyword_set(gcIDList) and keyword_set(accretionTimeSubset) then $
    message,'Error: since we take accModeInds must subset all of galcat, not just one/few halos.'
  if keyword_set(rVirNorm) and ~keyword_set(parNorm) then $
    message,'Error: Must specify parNorm (pri or sec) with rVirNorm.'
    
  ; default behavior: return 1 value per tracer for tracer sims, 1 value per gas particle for sph
  allTR = 0 ; sph
  if sP.trMCPerCell ne 0 then allTR = 1
  
  ; check input options vs. simulation type
  if keyword_set(curTracerVal) and sP.trMCPerCell le 0 then message,'Error: curTracerVal is trMC only.'
  if keyword_set(tracksFluid) and sP.trMCPerCell le 0 then message,'Error: tracksFluid is trMC only.'
  if keyword_set(elemIDs) and sP.trMCPerCell eq -1 then message,'Error: Not implemented.'
  
  ; ----- selection subset -----
   
  ; select primary,secondary,or all subhalos
  if ~keyword_set(gcIDList) then gcIDList = gcIDList(sP=sP,select=select)
  
  ; subset galcat member indlist by those with recorded accretion times (and associated properties)
  if keyword_set(accretionTimeSubset) then begin
    at = accretionTimes(sP=sP)
    
    ; select galaxyCat indices (separated for TR, not for SPH)
    galcatInds = galcatINDList(sP=sP,gcIDList=gcIDList,gcSep=allTR)
    
    ; this is used to access accretionTimes
    if n_elements(accMode) eq 0 then accMode = 'all'
    accTimeInds = accModeInds(at=at,sP=sP,accMode=accMode)

    ; sph case: modify galcatInds such that the accretionTimes subset is taken
    if ~allTR then begin
      for i=0,n_tags(accTimeInds)-1 do $
        galcatIndsNew = mod_struct( galcatIndsNew, (tag_names(accTimeInds))[i], $
          galcatInds[accTimeInds.(i)] )
          
      galcatInds = galcatIndsNew
    endif

    ; tracer case: handle only after we have child counts (after expansion or in allTR for maxTemps)
  endif else begin
    ; not taking atS (perhaps mtS/all pri, or just one or a few gcIDs)
    ; select galaxycat indices corresponding to the list of subgroup ids (no gcSep for SPH)
    galcatInds = galcatINDList(sP=sP,gcIDList=gcIDList,/gcSep)
  endelse
  
  ; load galaxy catalog to change INDs to gas IDs, or for element IDs
  galcat = galaxyCat(sP=sP)
    
  ; ----- values -----
  
  if keyword_set(accTime) then begin
    ;convert scale factors -> redshift
    for i=0,n_tags(galcat.types)-1 do $
      r = mod_struct( r, (tag_names(galcat.types))[i], $
        reform( 1/at.accTime[sP.atIndMode,accTimeInds.(i)]-1 ) )

    return,r
  endif
  
  if keyword_set(accTvir) then begin
    ; take accretionTime subset and return
    for i=0,n_tags(galcat.types)-1 do $
      r = mod_struct( r, (tag_names(galcat.types))[i], $
        at.accHaloTvir[accTimeInds.(i)] )

    return,r
  endif

  if keyword_set(rVirNorm) or keyword_set(virTemp) or keyword_set(parMass) then begin
    if keyword_set(rVirNorm) then begin
      ; load parent r_vir and galaxy radii catalog at sP.snap
      r_vir = galCatParentProperties(sP=sP, /rVir)
      gcr   = galaxyCatRadii(sP=sP)

      ; store radial distance normalized by parent r_vir
      ; note: if SO values not calculated (no subgroup in fof group), rvir=0 and rad->Inf (not plotted)
      ; note: since using most bound particle for subgroup centers, one per subgroup will have r=0
      if parNorm eq 'pri' then val = gcr.rad_pri / r_vir
      if parNorm eq 'sec' then val = gcr.rad_sec / r_vir
    endif

    if keyword_set(virTemp) then begin
      val = galCatParentProperties(sP=sP, /virTemp, parNorm=parNorm) 
    endif

    if keyword_set(parMass) then begin
      val = galCatParentProperties(sP=sP, /mass, parNorm=parNorm)
    endif

    ; restrict to subgroup type / merger tree / accretion time subset
    if ~allTR then begin
      for i=0,n_tags(galcat.types)-1 do $
        r = mod_struct( r, (tag_names(galcat.types))[i], val[galcatInds.(i)] )
      return,r
    endif

    ; allTR: replicate
    if sP.trMCPerCell gt 0 then child_counts = galcat.ccMC
    if sP.trMCPerCell lt 0 then message,'todo'

    ; create new return arrays
    r = fltarr( total(child_counts,/int) )

    ; take accretionTime subset and return tracer expanded array
    if keyword_set(accretionTimeSubset) then begin

      offset = 0L

      for i=0L,n_elements(val)-1 do begin
        if child_counts[i] gt 0 then begin
          r[offset:offset+child_counts[i]-1] = replicate(val[i],child_counts[i])
          offset += child_counts[i]
        endif
      endfor

      for i=0,n_tags(galcat.types)-1 do $
        rr = mod_struct(rr, (tag_names(galcat.types))[i], r[accTimeInds.(i)])
      return, rr
    endif

    ; not taking atS, so replicate into tracer expanded array by type, then take galcatInds
    for j=0,n_tags(galcat.types)-1 do begin
      r = mod_struct(r, (tag_names(galcat.types))[j], $
        fltarr(total(child_counts[galcatInds.(j)],/int)) )

      offset = 0L
      type_child_counts = child_counts[galcatInds.(j)]
      type_val = val[galcatInds.(i)]

      for i=0L,n_elements(galcatInds.(j))-1 do begin
        if type_child_counts[i] gt 0 then begin
          r.(j)[offset:offset+type_child_counts[i]-1] = replicate(type_val[i],type_child_counts[i])
          offset += type_child_counts[i]
        endif
      endfor

    endfor
    return, r

  endif
  
  if keyword_set(curTemp) then begin
    message,'todo'
    ; load gas ids and make ID->index map
    ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    idsIndMap = getIDIndexMap(ids,minid=minid)
    ids = !NULL
    
    ; load gas u and restrict to subset of galaxy cat
    u       = loadSnapshotSubset(sP=sP,partType='gas',field='u')
    u_gal   = u[idsIndMap[galcat.ids[galcatInds.gal]-minid]]
    u_gmem  = u[idsIndMap[galcat.ids[galcatInds.gmem]-minid]]
    u       = !NULL
    
    ; load gas nelec and restrict to subset of galaxy cat
    nelec       = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
    nelec_gal   = nelec[idsIndMap[galcat.ids[galcatInds.gal]-minid]]
    nelec_gmem  = nelec[idsIndMap[galcat.ids[galcatInds.gmem]-minid]]
    nelec       = !NULL
    
    idsIndMap = !NULL
    
    ; calculate temperature
    temp_gal   = convertUtoTemp(u_gal,nelec_gal,/log)
    temp_gmem  = convertUtoTemp(u_gmem,nelec_gmem,/log)
    temp_stars = fltarr(n_elements(galcat.stellarIDs[galcatInds.stars])) ; zeros, stars have no current temperature
     
    r = {gal:temp_gal,gmem:temp_gmem,stars:temp_stars,inter:-1,bhs:-1}
  endif
  
  if keyword_set(curSingleVal) then begin
    message,'todo'
    ; load gas ids and make ID->index map
    ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    idsIndMap = getIDIndexMap(ids,minid=minid)
    ids = !NULL
    
    ; load gas densities
    singleVal = loadSnapshotSubset(sP=sP,partType='gas',field=singleValField)
    
    ; restrict to subgroup type / merger tree / accretion time subset
    val_gal   = singleVal[idsIndMap[galcat.ids[galcatInds.gal]-minid]]
    val_gmem  = singleVal[idsIndMap[galcat.ids[galcatInds.gmem]-minid]]
    val_stars = fltarr(n_elements(galcatInds.stars)) ; zeros, stars don't have the same values as gas
    
    density = !NULL & idsIndMap = !NULL
    
    r = {gal:val_gal,gmem:val_gmem,stars:val_stars}
  endif
  
  if keyword_set(curTracerVal) then begin
    ; load some tracer field, single value per MC tracer
    tr_field  = loadSnapshotSubset(sP=sP,partType='tracerMC',field=singleTracerField)
    tr_parids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentids')
      
    field = tr_field[cosmoTracerChildren(sP=sP,gasIDs=galcat.ids,tr_parids=tr_parids,/getInds)]
      
    tr_field = !NULL
    tr_parids = !NULL
     
    message,'todo, note the cosmoTracerChildren splits gasIDs and starIDs, the above not going to work'
 
    r = {gal:field_gal,gmem:field_gmem,stars:field_stars}
      
      
    ; take accretionTime subset of mtS of all tracers and return
    if keyword_set(accretionTimeSubset) then begin
      r = { gal:r.gal[accTimeInds.gal], gmem:r.gmem[accTimeInds.gmem], stars:r.stars[acctimeInds.stars] }
    endif
    return,r
  endif
  
  if keyword_set(elemIDs) then begin
    message,'todo'
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
  
  if keyword_set(tracksFluid) then begin
    message,'todo'
    ; make indices for accretionTreeSubset (tracksFluid stored only for mtS/atS)  
    mtsInds = mergerTreeINDList(sP=sP,galcat=galcat,gcIDList=gcIDList)
    
    ; take accretionTime subset of mtS? if so modify indices now
    if keyword_set(accretionTimeSubset) then begin
      mtsInds = { gal   : mtsInds.gal[accTimeInds.gal]    ,$
                  gmem  : mtsInds.gmem[accTimeInds.gmem]  ,$
                  stars : mtsInds.stars[acctimeInds.stars] }
    endif
    
    ; all tracers requested, load directly and immediately return
    ; (could modify this to be like maxPastTemp when sph/trVel added to tracksFluid)
    tracks_gal   = tracksFluid(sP=sP,/loadAllTRGal)
    tracks_gmem  = tracksFluid(sP=sP,/loadAllTRGmem)
    tracks_stars = tracksFluid(sP=sP,/loadAllTRStars)

    ; return temps (logK), entropy (log CGS), or density (log code)
    r = {gal   : { temp : tracks_gal.temp[*,mtsInds.gal] ,$
                   ent  : tracks_gal.ent[*,mtsInds.gal]  ,$
                   dens : tracks_gal.dens[*,mtsInds.gal] ,$
                   rad  : tracks_gal.rad[*,mtsInds.gal] ,$
                   flag : tracks_gal.flag[*,mtsInds.gal]  }      ,$
         gmem  : { temp : tracks_gmem.temp[*,mtsInds.gmem] ,$
                   ent  : tracks_gmem.ent[*,mtsInds.gmem]  ,$
                   dens : tracks_gmem.dens[*,mtsInds.gmem] ,$
                   rad  : tracks_gmem.rad[*,mtsInds.gmem] ,$
                   flag : tracks_gmem.flag[*,mtsInds.gmem]  }    ,$
         stars : { temp : tracks_stars.temp[*,mtsInds.stars] ,$
                   ent  : tracks_stars.ent[*,mtsInds.stars]  ,$
                   dens : tracks_stars.dens[*,mtsInds.stars] ,$
                   rad  : tracks_stars.rad[*,mtsInds.stars] ,$
                   flag : tracks_stars.flag[*,mtsInds.stars]  }  ,$
         rr    : tracks_gal.rr                                       }
        
    tracks_gal = !NULL & tracks_gmem = !NULL & tracks_stars = !NULL
      
    return,r
  endif
  
  if keyword_set(maxPastTemp) or keyword_set(maxTempTime) or $
     keyword_set(maxPastEnt)  or keyword_set(maxPastMach)  or $
     keyword_set(maxPastDens) then begin
     
    ; if all tracers requested, load directly and immediately return
    if keyword_set(allTR) then begin
      maxval = maxVals(sP=sP)
          
      ; return temps (logK), entropy (log CGS), or density (log code)
      if keyword_set(maxPastTemp) then maxvalInd = ( where(tag_names(maxval) eq 'MAXTEMPS') )[0]
      if keyword_set(maxPastEnt)  then maxvalInd = ( where(tag_names(maxval) eq 'MAXENT') )[0]
      if keyword_set(maxPastDens) then maxvalInd = ( where(tag_names(maxval) eq 'MAXDENS') )[0]
      if keyword_set(maxPastMach) then maxvalInd = ( where(tag_names(maxval) eq 'MAXMACHNUM') )[0]
      if keyword_set(maxTempTime) then maxvalInd = ( where(tag_names(maxval) eq 'MAXTEMPTIME') )[0]

      r = maxval.(maxvalInd)
      if keyword_set(maxTempTime) then r = 1/r-1 ; scalefac to redshift

      maxval = !NULL
      
      ; take accretionTime subset of mtS of all tracers and return
      if keyword_set(accretionTimeSubset) then begin
        galcatInds = galcatINDList(sP=sP,gcIDList=gcIDList,/trRep)

        ; index r with galcatInds[accTimeInds.X]
        for i=0,n_tags(galcat.types)-1 do $ ; n_tags(accTimeInds) includes mask
          rr = mod_struct( rr, (tag_names(accTimeInds))[i], r[galcatInds[accTimeInds.(i)]] )
 
       return,rr
      endif
      
      ; not taking atS, get separated galcatInds and split
      galcatInds = galcatINDList(sP=sP,gcIDList=gcIDList,/gcSep,/trRep)

      for i=0,n_tags(accTimeInds)-1 do $
        rr = mod_struct( rr, (tag_names(galcatInds))[i], r[galcatInds.(i)] )

      return,rr
    endif
   
    ; otherwise, load maximum past temperature per gas particle
    ; (or statistics for the child population of each gas cell, this code removed for now)
    maxval = maxVals(sP=sP)

    ; return temps (logK), entropy (log CGS), or density (log code)
    if keyword_set(maxPastTemp) then maxvalInd = ( where(tag_names(maxval) eq 'MAXTEMPS') )[0]
    if keyword_set(maxPastEnt)  then maxvalInd = ( where(tag_names(maxval) eq 'MAXENT') )[0]
    if keyword_set(maxPastDens) then maxvalInd = ( where(tag_names(maxval) eq 'MAXDENS') )[0]
    if keyword_set(maxPastMach) then return,-1 ; no max mach number for SPH
    if keyword_set(maxTempTime) then maxvalInd = ( where(tag_names(maxval) eq 'MAXTEMPTIME') )[0]
    
    val = maxval.(maxvalInd)
    if keyword_set(maxTempTime) then val = 1/val-1 ; scalefac to redshift
    
    maxval = !NULL
    
    if sP.trMCPerCell ne 0 then message,'error: how did we get here, todo'
    
    ; restrict temps/times to subset of galaxy cat (sph)
    for i=0,n_tags(galcat.types)-1 do $
      r = mod_struct( r, (tag_names(galcat.types))[i], $
        val[galcatInds.(i)] )
        
    return,r
        
  endif
  
  ; ----- tracer expansion if requested -----
 
  message,'move return inside'
 
  if keyword_set(allTR) then begin
    ; load child counts for both galaxy members and group members
    if sP.trMCPerCell gt 0 then child_counts = galcat.ccMC
    if sP.trMCPerCell lt 0 then message,'todo'
    
    ; create new return arrays
    for i=0,n_tags(galcat.types)-1 do begin
    
      if size(r.gal,/tname) eq 'LONG64' then $ ; to support curSingleVal=ids
        rr = mod_struct(rr, (tag_names(galcat.types))[i], $
        lon64arr(total(child_counts[galcatInds.(i)],/int)) )
             
      if size(r.gal,/tname) eq 'LONG' then $
        rr = mod_struct(rr, (tag_names(galcat.types))[i], $
        lonarr(total(child_counts[galcatInds.(i)],/int)) )  
             
      if size(r.gal,/tname) eq 'FLOAT' then $
        rr = mod_struct(rr, (tag_names(galcat.types))[i], $
        fltarr(total(child_counts[galcatInds.(i)],/int)) )
        
    endfor
    
    if n_elements(rr) eq 0 then message,'Error: Unknown type for tracer replication.'
    
    ; replicate byGas arrays into byTracer arrays
    for j=0,n_tags(galcat.types)-1 do begin
      offset = 0L
      type_child_counts = child_counts[galcatInds.(j)]
      
      for i=0L,n_elements(r.(j))-1 do begin
        if type_child_counts[i] gt 0 then begin
          rr.(j)[offset:offset+type_child_counts[i]-1] = replicate(r.(j)[i],type_child_counts[i])
          offset += type_child_counts[i]
        endif
      endfor
      
    endfor
    
    stop
    
    ; take accretionTime subset and return tracer expanded array
    if keyword_set(accretionTimeSubset) then $
      for j=0,n_tags(galcat.types)-1 do $
        rr.(j) = rr.(j)[accTimeInds.(j)]

    return,rr
  endif ; allTR
  
  return,r
end


