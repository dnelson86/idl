; galaxyCat.pro
; gas accretion project - gas selections of interest (galaxy/halo catalogs)
; dnelson jul.2013

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
        
    ; load ids of particles in all primary subfind groups
    gc = loadGroupCat(sP=sP,/readIDs)
    gcPIDs = gcPIDList(gc=gc,select='pri',partType='gas')

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
    gcPIDs_stars = gcPIDList(gc=gc,select='pri',partType='stars')
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
      gcPIDs_BHs = gcPIDList(gc=gc,select='pri',partType='BHs')
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
  
    ; replicate group parent IDs (of PRIMARY/parent) to each member particle
    ptNum = { gas : partTypeNum('gas'), stars : partTypeNum('stars'), BHs : partTypeNum('BHs') }
    
    sgIDs_pri = gcIDList(gc=gc,select='pri')    
    
    sgParIDs_gas   = lonarr(total(gc.subgroupLenType[ptNum.gas,sgIDs_pri],/int))
    sgParIDs_stars = lonarr(total(gc.subgroupLenType[ptNum.stars,sgIDs_pri],/int))
    
    if sP.gfmBHs ne 0 then $
      sgParIDs_BHs = lonarr(total(gc.subgroupLenType[ptNum.BHs,sgIDs_pri],/int))
    
    offset = { gas : 0L, stars: 0L, BHs: 0L }
      
    foreach sgID,sgIDs_pri do begin
      if gc.subgroupLenType[ptNum.gas,sgID] gt 0 then begin
        sgParIDs_gas[offset.gas:offset.gas+gc.subgroupLenType[ptNum.gas,sgID]-1] = $
          cmreplicate(sgID,gc.subgroupLenType[ptNum.gas,sgID])
        offset.gas += gc.subgroupLenType[ptNum.gas,sgID]
      endif
       
      if gc.subgroupLenType[ptNum.stars,sgID] gt 0 then begin
        sgParIDs_stars[offset.stars:offset.stars+gc.subgroupLenType[ptNum.stars,sgID]-1] = $
          cmreplicate(sgID,gc.subgroupLenType[ptNum.stars,sgID])
        offset.stars += gc.subgroupLenType[ptNum.stars,sgID]
      endif
      
      if gc.subgroupLenType[ptNum.BHs,sgID] gt 0 then begin
        sgParIDs_BHs[offset.BHs:offset.BHs+gc.subgroupLenType[ptNum.BHs,sgID]-1] = $
          cmreplicate(sgID,gc.subgroupLenType[ptNum.BHs,sgID])
        offset.BHs += gc.subgroupLenType[ptNum.BHs,sgID]
      endif
    endforeach
      
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
    wStars = where(rad_pri lt sP.radcut_rvir,countStars,comp=wStarsComp,ncomp=wStarsComp_num)

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
          ;ids   : lon64arr(countTot)      ,$ ; add below
          type  : intarr(countTot)         ,$
          rad   : fltarr(countTot)         ,$
          vrad  : fltarr(countTot)         ,$
          types : types                    ,$
          trMC_cc   : -1                   ,$
          trMC_ids  : -1                   ,$
          trVel_cc  : -1                   ,$
          trVel_ids : -1                   ,$
          countGal: countGal, countGmem: countGmem, countStars: countStars ,$
          countInter: countInter, countBHs: countBHs, countTot:countTot,$
          nGroups: gc.nSubgroupsTot }
          
    if size(gc.IDs,/tname) eq 'LONG'    then r = create_struct(r, {IDs:lonarr(countTot)})
    if size(gc.IDs,/tname) eq 'LONG64'  then r = create_struct(r, {IDs:lon64arr(countTot)})
    if size(gc.IDs,/tname) eq 'ULONG64' then r = create_struct(r, {IDs:ulon64arr(countTot)})
    
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
    gcIDs_pri = gcIDList(gc=gc,select='pri')
    par_inds = gcIDs_pri[ value_locate(gc.subgroupOffset[gcIDs_pri],gc_ind) ]
      
    if min(par_inds) lt 0 or max(par_inds) ge n_elements(gc.subgroupOffset) then message,'Error'
    
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
      galcat_trids = cosmoTracerChildren(sP=sP, /getIDs, gasIDs=r.ids, child_counts=galcat_cc)
      r = mod_struct( r, 'trMC_cc', galcat_cc )
      r = mod_struct( r, 'trMC_ids', galcat_trids )
    endif
    
    if sP.trVelPerCell gt 0 then begin
      galcat_trids = cosmoTracerVelChildren(sP=sP,/getIDs,gasIDs=r.ids,child_counts=galcat_cc)
      r = mod_struct( r, 'trVel_cc', galcat_cc )
      r = mod_struct( r, 'trVel_ids', galcat_trids )
    endif
    
    ; calculate radius and radial velocity for each member
    h = loadSnapshotHeader(sP=sP)
    gcParIDs = galCatRepParentIDs(galcat=r)
    
    ; restrict gas particle positions to gal/gmem gas only
    parTypes = ['gas','stars','BHs']
    
    foreach parType,parTypes do begin
    
      ; skip non-existent particle types
      if h.nPartTot[partTypeNum(parType)] eq 0 then continue

	  ; sub-match to parents of this type
      ids_type = loadSnapshotSubset(sP=sP,partType=parType,field='ids')
       
      calcMatch,r.ids,ids_type,galcat_ind,ids_type_ind,count=countType
      if countType eq 0 then continue
      
      ids_type_ind = ids_type_ind[calcSort(galcat_ind)] ; rearrange indices to be in the order of r.ids
      galcat_ind   = galcat_ind[calcSort(galcat_ind)]
      
      gcParIDs_type = gcParIDs[galcat_ind]
      
      if countType eq 0 then print,'Warning: Skipping parent type '+parType
    
      ids_type = !NULL
    
      ; calculate radial distances to parent
      pos = loadSnapshotSubset(sP=sP,partType=parType,field='pos')
      pos = pos[*,ids_type_ind]
            
      ; calulate radial vector of gas from group center
      r.rad[galcat_ind] = periodicDists(subgroupCen[*,gcParIDs_type],pos,sP=sP)
      
      ; replace coordinates by relative coordinates (radial vectors) to direct parent
      for i=0,2 do begin
        pos_rel = reform(pos[i,*] - subgroupCen[i,gcParIDs_type])
        correctPeriodicDistVecs, pos_rel, sP=sP
        pos[i,*] = pos_rel
      endfor
    
      ; load velocities
      vel = loadSnapshotSubset(sP=sP,partType=parType,field='vel')
      vel = vel[*,ids_type_ind]

      r.vrad[galcat_ind] = ((vel[0,*] - gc.subgroupVel[0,gcParIDs_type]) * pos[0,*] + $
                            (vel[1,*] - gc.subgroupVel[1,gcParIDs_type]) * pos[1,*] + $
                            (vel[2,*] - gc.subgroupVel[2,gcParIDs_type]) * pos[2,*]) $
                            / r.rad
    
    endforeach
    
    ; immediate return and skip save?
    if keyword_set(skipSave) then return, r
      
    ; save galaxy catalog
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))

  endfor ;snapRange
  
  ; if just one snap requested and just calculated, return it now
  if snapRange[0] eq snapRange[1] then return, r
  
end

; galCatRepParentIDs(): for the galaxy catalog, replicate the list of ordered parent IDs such that
;                       the return array is the same size as the number of gas particles with
;                       each element the id of its parent subgroup/galaxy group/groupmember group
;
;
; gcIDList : return only a replicated parent ID list of the specified subgroups in the groupcat
; child_counts: return a replicated list for tracers which have multiplicity inside each gas cell
;               or star particle as specified in child_counts

function galCatRepParentIDs, galcat=galcat, gcIDList=gcIDList, child_counts=child_counts

  compile_opt idl2, hidden, strictarr, strictarrsubs

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
            
        ; replicate for children
        if galcat.len[gcID] gt 0 then begin
          tot_children = total(child_counts[offset_c:offset_c+galcat.len[gcID]-1],/int)
          ind_end = offset + tot_children - 1

          if tot_children gt 0 then r[offset:ind_end] = cmreplicate(gcID,tot_children)
          offset += tot_children
          offset_c += galcat.len[gcID]
        endif
        
    endforeach

    return,r
end

; galCatParentProperties: calculate some property of the parent galaxy/group for every gas elements
;                         in the galaxy catalog at some snapshot
;
; trRep=1 : return a replicated list for tracers which have some multiplicity inside each gas cell
;
; virTemp=1  : virial temperature
; mass=1     : total mass (from catalog, dm+baryon)
; rVir=1     : virial radius (r_200 critical)
; rVirNorm=1 : virial radius / r200
; vCirc=1    : circular velocity (v_200)

function galCatParentProperties, sP=sP, galcat=galcat, trRep=trRep, $
                                 virTemp=virTemp, mass=mass, rVir=rVir, rVirNorm=rVirNorm, vCirc=vCirc

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function galaxyCat, snapNumToRedshift, codeMassToLogMsun

  ; load group catalog for masses, galaxy catalog
  gc = loadGroupCat(sP=sP,/skipIDs)
  if ~keyword_set(galcat) then galcat = galaxyCat(sP=sP)

  ; replicate parent IDs
  gcInd = galCatRepParentIDs(galcat=galcat)
   
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
  
  if keyword_set(rVirNorm) then begin
    ; note: if SO values not calculated (no subgroup in fof group), rvir=0 and rad->Inf (not plotted)
    ; note: since using most bound particle for subgroup centers, one per subgroup will have r=0
    r = gc.subgroupGrnr[gcInd]
    r = gc.group_r_crit200[r]
	r = galcat.rad / r ; r/rvir
  endif
  
  if keyword_set(vCirc) then begin
    units = getUnits()
    
    r = gc.subgroupGrnr[gcInd]
    r = sqrt(units.G * gc.group_m_crit200[r] / gc.group_r_crit200[r])
  endif
  
  ; replicate for tracers?
  if keyword_set(trRep) then begin
    if sP.trMCPerCell gt 0 then r = r[ replicate_var(galcat.trMC_cc) ]
	if sP.trMCPerCell lt 0 then r = r[ replicate_var(galcat.trVel_cc) ]
	if sP.trMCPerCell eq 0 then message,'Error: No trRep on SPH run.'
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
  
  if keyword_set(trRep) and sP.trMCPerCell eq 0 then message,'Error.'
  
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
  if sP.trMCPerCell gt 0 then child_counts = galcat.trMC_cc
  if sP.trMCPerCell lt 0 then child_counts = galcat.trVel_cc

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
  if sP.trMCPerCell gt 0 then xx = replicate_var(galcat.trMC_cc, subset_inds=r)
  if sP.trMCPerCell gt 0 then xx = replicate_var(galcat.trVel_cc, subset_inds=r)
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

; gcSubsetProp(): read galaxy catalog for a specific subgroup selection (gcIDList) or the atS,
;                 for each gas element/tracer (may or may not depend on parent halo)
;                 at the redshift specified by sP.snap (may or may not depend on previous time)
; default behavior: return 1 value per tracer for tracer sims, 1 value per gas particle for sph
;
; gcIDList=[]   : an explicit list of (primary) subgroups to select, otherwise return for all
;
; rVirNorm=1    : radial distances normalized by r_vir of either primary or secondary parent
; virTemp=1     : current virial temperatures of parent halos
; parMass=1     : total mass (dm+baryonic) of parent halos (from catalog)
; elemIDs=1       : ids of each element (either SPH particles or tracers)
; tracksFluid=1   : temp,dens,ent history of each tracer (trMC only)
;
; maxPastTemp=1 : maximum past previous temperature of each element
; maxTempTime=1 : time when maximum past previous temperature was reached (in redshift)
; maxPastEnt=1  : maximum past entropy
; maxPastDens=1 : maximum past density
; maxPastMach=1 : maximum past mach number
;  trPopMin,trPopMean,trPopMax : return the respective statistic for each gas cell for the tracers
;
; curGasVal=1    : current single quantity (e.g. mass, density) returned without manipulation (gas, or replicated)
; curTracerVal=1 : as above, but individual value for each tracer (e.g. wind_counter)
;  curField : field name in snapshot file for the above
;
; accretionTimeSubset : return values only for the subset of particles/tracers with recorded accretion times
;  timeWindow : either in Myr, or 'all', the time over which to include accretion events
;  accTime,accTvir : time of accretion (in redshift) or virial temp of parent halo at time of accretion
;  accMode : return values only for one accretionMode (all,smooth,bclumpy,sclumpy,smooth)

function gcSubsetProp, sP=sP, gcIDList=gcIDList, $
           rVirNorm=rVirNorm, virTemp=virTemp, parMass=parMass, $
           curTemp=curTemp, maxPastTemp=maxPastTemp, maxTempTime=maxTempTime, $
           maxPastEnt=maxPastEnt, maxPastDens=maxPastDens, maxPastMach=maxPastMach, $
           elemIDs=elemIDs, tracksFluid=tracksFluid, $
           curGasVal=curGasVal, curTracerVal=curTracerVal, curField=curField, $
           accretionTimeSubset=accretionTimeSubset, timeWindow=TW, $
           accTime=accTime,accTvir=accTvir,accMode=accMode ; for accretionTimeSubset only

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; check combinations of input options for validity
  if (keyword_set(accTime) or keyword_set(accTvir)) and ~keyword_set(accretionTimeSubset) then $
    message,'Error: Can only return accretion time or Tvir at accretion time for atS.'
  if keyword_set(accMode) and ~keyword_set(accretionTimeSubset) then $
    message,'Error: Can only return accretion mode subsets of the accretionTime subset.'
  if keyword_set(TW) and ~keyword_set(accretionTimeSubset) then $
    message,'Error: Can only enforce timeWindow restriction on the atS.'
  if keyword_set(curField) and (~keyword_set(curGasVal) and ~keyword_set(curTracerVal)) then $
    message,'Error: Must specify either gas or tracer for a single current field.'
	
  ; check input options vs. simulation type
  if keyword_set(curTracerVal) and sP.trMCPerCell le 0 then message,'Error: curTracerVal is trMC only.'
  if keyword_set(tracksFluid) and sP.trMCPerCell le 0 then message,'Error: tracksFluid is trMC only.'
  if keyword_set(elemIDs) and sP.trMCPerCell eq -1 then message,'Error: Not implemented.'
  
  ; load galaxy catalog to change INDs to gas IDs, or for element IDs
  galcat = galaxyCat(sP=sP)
  
  ; ----- selection subset -----
  
  ; subset galcat member indlist by those with recorded accretion times (and associated properties)
  if keyword_set(accretionTimeSubset) then begin
    at = accretionTimes(sP=sP)
    
    ; unseparated, for either each sph particle or each tracer
    if n_elements(accMode) eq 0 then accMode = 'all'
    galcatInds = accModeInds(at=at,sP=sP,accMode=accMode)

    ; do the timewindow restriction immediately:
    if keyword_set(TW) then begin
	
      ; NOTE: this TW option only used for maxVals binning, is not consistent with TW approach for accRates
      
      ; current time
      h = loadSnapshotHeader(sP=sP)
      curtime = 1/h.time - 1 ; redshift
      curtime = redshiftToAgeFlat(curtime)*1e9 ; yr
  
      if str(TW) eq 'all' then begin
        timeWindow = curtime - redshiftToAgeFlat(6.0)*1e9 ; go back to z=6 (in yr)
      endif else begin
        timeWindow = TW * 1e6 ; convert input Myr to yr
      endelse

      ; GALAXY (1): accretion defined as (rho,temp) joining time or 0.15rvir crossing time (most recent)
      loc_atime = reform(at.accTimeRT[galcatInds.gal])
    
      r_crossing_time = reform(at.accTime[sP.radIndGalAcc,galcatInds.gal])
      w = where(r_crossing_time gt loc_atime,count)
      if count gt 0 then loc_atime[w] = r_crossing_time[w]
      
      loc_atime = 1/loc_atime - 1 ; redshift
      loc_atime = redshiftToAgeFlat(loc_atime)*1e9 ; yr
    
      ; make a count of those falling in the time window
      w_type = where(curtime - loc_atime le timeWindow,nloc)
      galcatInds = mod_struct( galcatInds, 'gal', galcatInds.gal[w_type] )
	
      ; GMEM (2)
      loc_atime = reform(at.accTime[sP.radIndHaloAcc,galcatInds.gmem])
      loc_atime = 1/loc_atime - 1 ; redshift
      loc_atime = redshiftToAgeFlat(loc_atime)*1e9 ; yr
    
      ; make a count of those falling in the time window
      w_type = where(curtime - loc_atime le timeWindow,nloc)
      galcatInds = mod_struct( galcatInds, 'gmem', galcatInds.gmem[w_type] )
	
      ; INTER (3) - look at outflow time within this TW
      loc_atime = reform(at.outTimeRT[galcatInds.inter])
    
      r_crossing_time = reform(at.outTime[sP.radIndGalAcc,galcatInds.inter])
      w = where(r_crossing_time gt loc_atime,count)
      if count gt 0 then loc_atime[w] = r_crossing_time[w]
      
      loc_atime = 1/loc_atime - 1 ; redshift
      loc_atime = redshiftToAgeFlat(loc_atime)*1e9 ; yr
    
      ; make a count of those falling in the time window
      w_type = where(curtime - loc_atime le timeWindow,nloc)
	galcatInds = mod_struct( galcatInds, 'inter', galcatInds.inter[w_type] )
    
      ; STARS (4)
      loc_atime = reform(at.accTimeRT[galcatInds.stars])
    
      r_crossing_time = reform(at.accTime[sP.radIndGalAcc,galcatInds.stars])
      w = where(r_crossing_time gt loc_atime,count)
      if count gt 0 then loc_atime[w] = r_crossing_time[w]
    
      ; convert from scale factor to age of the universe
      loc_atime = 1/loc_atime - 1 ; redshift
      loc_atime = redshiftToAgeFlat(loc_atime)*1e9 ; yr
    
      ; make a count of those falling in the time window
      w_type = where(curtime - loc_atime le timeWindow,nloc)
      galcatInds = mod_struct( galcatInds, 'stars', galcatInds.stars[w_type] )
    
      ; BHs (5)
      if sP.gfmBHs ne 0 then begin
        loc_atime = reform(at.accTimeRT[galcatInds.bhs])
    
        r_crossing_time = reform(at.accTime[sP.radIndGalAcc,galcatInds.bhs])
        w = where(r_crossing_time gt loc_atime,count)
        if count gt 0 then loc_atime[w] = r_crossing_time[w]
    
        loc_atime = 1/loc_atime - 1 ; redshift
        loc_atime = redshiftToAgeFlat(loc_atime)*1e9 ; yr
    
        ; make a count of those falling in the time window
        w_type = where(curtime - loc_atime le timeWindow,nloc)
	  galcatInds = mod_struct( galcatInds, 'bhs', galcatInds.bhs[w_type] )
      endif
	
    endif ; TW
	
    ; handle if not selecting all primary subgroups
    if keyword_set(gcIDList) then begin
      galcatInds_list = galcatINDList(sP=sP,gcIDList=gcIDList,/gcSep)
    
      ; intersect with galcatInds, and take subset
      for i=0,n_tags(galcatInds)-1 do begin
        inds_type = intersection( galcatInds.(i), galcatInds_list.(i) )
        galcatInds = mod_struct( galcatInds, (tag_names(galcatInds))[i], inds_type )
      endfor
      
    endif
      
  endif else begin
    ; not taking atS: take either gcIDList or all in galaxyCat
    if n_elements(gcIDList) eq 0 then gcIDList = lindgen(galcat.nGroups)
    galcatInds = galcatINDList(sP=sP,gcIDList=gcIDList,/gcSep)
  endelse
  
  ; check for blank entries?
  for i=0,n_tags(galcatInds)-1 do $
    if n_elements(galcatInds.(i)) eq 1 then $
      if galcatInds.(i) eq -1 then message,'handle this'
    
  ; ----- values -----
  
  if keyword_set(accTime) then begin
    ;convert scale factors -> redshift
    for i=0,n_tags(galcat.types)-1 do $
      r = mod_struct( r, (tag_names(galcat.types))[i], $
        reform( 1/at.accTime[sP.atIndMode,galcatInds.(i)]-1 ) )

    return,r
  endif
  
  if keyword_set(accTvir) then begin
    ; take accretionTime subset and return
    for i=0,n_tags(galcat.types)-1 do $
      r = mod_struct( r, (tag_names(galcat.types))[i], $
        at.accHaloTvir[galcatInds.(i)] )

    return,r
  endif

  if keyword_set(rVirNorm) or keyword_set(virTemp) or keyword_set(parMass) then begin
  
    trRep = sP.trMCPerCell ne 0
  
    if keyword_set(rVirNorm) then val = galCatParentProperties(sP=sP, galcat=galcat, trRep=trRep, /rVirNorm)
    if keyword_set(virTemp)  then val = galCatParentProperties(sP=sP, galcat=galcat, trRep=trRep, /virTemp) 
    if keyword_set(parMass)  then val = galCatParentProperties(sP=sP, galcat=galcat, trRep=trRep, /mass)

    ; separate and return
    for i=0,n_tags(galcat.types)-1 do $
      r = mod_struct( r, (tag_names(galcat.types))[i], val[galcatInds.(i)] )
	  
    return,r

  endif
  
  if keyword_set(curGasVal) then begin
    ; single current gas field: load gas ids and do match
    ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
	
    calcMatch,ids,galcat.ids,ids_ind,galcat_ind,count=countMatch
    ids_ind    = ids_ind[sort(galcat_ind)]
    galcat_ind = galcat_ind[sort(galcat_ind)]
    
    val = fltarr(n_elements(galcat.ids))
	
    if curField eq 'temp' then begin
      ; load gas u, nelec and calculate temperature
      u     = loadSnapshotSubset(sP=sP,partType='gas',field='u',inds=ids_ind)
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec',inds=ids_ind)
	val[galcat_ind] = convertUtoTemp(u,nelec,/log)
	  
	u = !NULL
	nelec = !NULL
    endif else begin
	val[galcat_ind] = loadSnapshotSubset(sP=sP,partType='gas',field=curField,inds=ids_ind)
    endelse
	
    ; replicate if curGasVal (one per cell) and using tracers (either type)
    if sP.trMCPerCell gt 0 then val = val[ replicate_var(galcat.trMC_cc) ]
    if sP.trMCPerCell lt 0 then val = val[ replicate_var(galcat.trVel_cc) ]
	
    ; separate and return
    for i=0,n_tags(galcat.types)-1 do $
      r = mod_struct( r, (tag_names(galcat.types))[i], val[galcatInds.(i)] )
	  
    return,r

  endif
  
  if keyword_set(curTracerVal) then begin
    ; load some tracer field, single value per MC tracer
    ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
	
    calcMatch,ids,galcat.trMC_ids,ids_ind,galcat_ind,count=countMatch
    ids_ind    = ids_ind[sort(galcat_ind)]
    galcat_ind = galcat_ind[sort(galcat_ind)]
    
    val = fltarr(n_elements(galcat.trMC_ids))
	
    load = loadSnapshotSubset(sP=sP,partType='tracerMC',field=curField)
    val[galcat_ind] = load[ids_ind]
	
    ; separate and return
    for i=0,n_tags(galcat.types)-1 do $
      r = mod_struct( r, (tag_names(galcat.types))[i], val[galcatInds.(i)] )
	  
    return,r
	
  endif
  
  if keyword_set(elemIDs) then begin
    ; which ids to subselect within?
    if sP.trMCPerCell eq 0 then ids = galcat.ids
    if sP.trMCPerCell gt 0 then ids = galcat.trMC_ids
    if sP.trMCPerCell lt 0 then ids = galcat.trVel_ids
    
    ; separate and return
    for i=0,n_tags(galcat.types)-1 do $
      r = mod_struct( r, (tag_names(galcat.types))[i], ids[galcatInds.(i)] )
	  
    return,r
    
  endif
  
  if keyword_set(tracksFluid) then begin
    ; all tracers requested, load directly and immediately return
    ; (could modify this to be like maxPastTemp when sph/trVel added to tracksFluid)
    tracks = tracksFluid(sP=sP)

    ; separate and return
    for i=0,n_tags(galcat.types)-1 do $
      r = mod_struct( r, (tag_names(galcat.types))[i], $
                      { temp : tracks.temp[*,galcatInds.(i)] ,$
                        ent  : tracks.ent [*,galcatInds.(i)] ,$
                        dens : tracks.dens[*,galcatInds.(i)] ,$
                        rad  : tracks.rad [*,galcatInds.(i)] ,$
                        flag : tracks.flag[*,galcatInds.(i)]  } )
                        
    r = mod_struct( r, 'rr', tracks.rr )
    
    return,r
  endif
  
  if keyword_set(maxPastTemp) or keyword_set(maxTempTime) or $
     keyword_set(maxPastEnt)  or keyword_set(maxPastMach)  or $
     keyword_set(maxPastDens) then begin
	 
    maxval = maxVals(sP=sP)
          
    ; return temps (logK), entropy (log CGS), or density (log code)
    if keyword_set(maxPastTemp) then maxvalInd = ( where(tag_names(maxval) eq 'MAXTEMPS') )[0]
    if keyword_set(maxPastEnt)  then maxvalInd = ( where(tag_names(maxval) eq 'MAXENT') )[0]
    if keyword_set(maxPastDens) then maxvalInd = ( where(tag_names(maxval) eq 'MAXDENS') )[0]
    if keyword_set(maxPastMach) then maxvalInd = ( where(tag_names(maxval) eq 'MAXMACHNUM') )[0]
    if keyword_set(maxTempTime) then maxvalInd = ( where(tag_names(maxval) eq 'MAXTEMPTIME') )[0]

    maxval = maxval.(maxvalInd)
    if keyword_set(maxTempTime) then maxval = 1/maxval-1 ; scalefac to redshift
      
    ; if SPH and machNum (doesn't exist), return zeros
    if sP.trMCPerCell eq 0 and keyword_set(maxPastMach) then begin
      for i=0,n_tags(galcat.types)-1 do $
        r = mod_struct( r, (tag_names(galcat.types))[i], fltarr(n_elements(galcatInds.(i))) )
	return, r
    endif
	
    ; otherwise, separate and return
    for i=0,n_tags(galcat.types)-1 do $
      r = mod_struct( r, (tag_names(galcat.types))[i], maxval[galcatInds.(i)] )
	
    return,r
      
  endif
  
end
