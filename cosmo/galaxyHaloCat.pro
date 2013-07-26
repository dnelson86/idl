; galaxyHaloCat.pro
; halo cooling project - gas selections of interest (halo using full fof)
; dnelson may.2013

; galaxyHaloCat(): make a catalog similar to gmem, but include all of the fof minus (rho,T) cut, 
;                  radial restrictions  and all non-primary subhalos
; Note: the halo catalog has the same size/indexing as the subgroups and the original galaxyCat(), 
;       but all non-primary have zero length

function galaxyHaloCat, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs

  radCutOuter = 2.0 ; fraction of rvir, instead of sP.radcut_out which is 1.5 for gas acc paper
  
  ; check existence of save
  saveFilename = sP.derivPath + 'haloCat.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) +'.sav'
    
  ; results exist, return
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif
    
  ; galaxyHaloCat requires a minimum rvir cut to be defined, or else it need be hardcoded here
  if sP.radcut_rvir eq 0.0 then message,'Error: sP has no minimum rVir radial cut.'
    
  ; load ids of particles in all subfind groups
  gc = loadGroupCat(sP=sP,/readIDs)
  
  ; make candidate gcPIDs list: exclude all gas in secondary subhalos
  ptGas = partTypeNum('gas')
  
  gcIDs_sec   = gcIDList(gc=gc,select='sec')
  gcPIDs_mask = bytarr(n_elements(gc.IDs))
  
  foreach gcID,gcIDs_sec do begin
    if gc.subgroupLenType[ptGas,gcID] gt 0 then begin
      gcPIDs_mask[gc.subGroupOffsetType[ptGas,gcID] : gc.subGroupOffsetType[ptGas,gcID] + $
        gc.subGroupLenType[ptGas,gcID] - 1] = 1B
    endif
  endforeach
  
  gcIDs_ind = where(gcPIDs_mask eq 0B,count)
  if count eq 0 then message,'Error: No halo gas remaining after removing secondary subhalos.'
  gcPIDs = gc.IDs[gcIDs_ind]
  
  gcPIDs_mask = !NULL
  gcIDs_sec   = !NULL
  
  ; load gas ids and match to catalog
  ids  = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
  calcMatch,gcPIDs,ids,gc_ind,ids_ind,count=countMatch
  
  ids_ind = ids_ind[calcSort(gc_ind)] ; rearrange ids_ind to be in the order of gcPIDs
  gc_ind  = gc_ind[calcSort(gc_ind)] ; rearrange gc_ind the same way
  
  gcIDs_ind = gcIDs_ind[gc_ind] ; take gcIDs_ind subset in original order
  ids = ids[ids_ind] ; take gas ids subset in same order
  
  if ~array_equal(gc.IDs[gcIDs_ind],ids) then message,'Error: id mismatch.'
  
  ; sanity check: number found should be total gas in fof - total gas in secondary subhalos
  totGas_sec = n_elements(gcPIDList(gc=gc,select='sec',partType='gas'))
  totGas_fof = total(gc.groupLenType[partTypeNum('gas'),*],/int)
  if countMatch ne totGas_fof-totGas_sec then message,'Error: Failed to locate all gas gcPIDs in gas ids.'

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
    
  ; calculate radial distances of gas elements to primary parents
  pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos',inds=ids_ind)
  vel = loadSnapshotSubset(sP=sP,partType='gas',field='vel',inds=ids_ind)
  
  ; find group center positions with most bound particles for each group
  subgroupCen = subgroupPosByMostBoundID(sP=sP)
  
  ; fof has rvir=0 (no SO values) if zero subgroups, marginal overdensity, or low total mass
  ; these are low mass halos which we aren't going to plot anyways
  w = where(gc.group_r_crit200 eq 0.0,count)
  if count gt 0 then gc.group_r_crit200[w] = 1e-8 ; remove with outer radial cut
  
  ; replicate group parent IDs (of PRIMARY/parent) to each member particle
  groupParIDs = value_locate(gc.groupOffset, gcIDs_ind)
  
  if min(groupParIDs) lt 0 or max(groupParIDs) ge gc.nGroupsTot then message,'Error: Bad parent ID replication.'
      
  ; make subgroup parent list too, but note that groupFirstSub can be wrong if e.g. group has 
  ; no subgroups, and can overrun number of subgroups for the last group
  sgParIDs = gc.groupFirstSub[groupParIDs] < (gc.nSubgroupsTot-1L)

  wSGOk = where(gc.subgroupGrNr[sgParIDs] eq groupParIDs,countSGOk,$
                comp=wSGNo,ncomp=countSGNo)     
      
  ; calculate radial distance of each gas cell (from SG center if possible, else from fof center)
  rad_pri = fltarr(n_elements(sgParIDs)) - 1.0
  
  if countSGOk gt 0 then $
    rad_pri[wSGOk] = periodicDists(subgroupCen[*,sgParIDs[wSGOk]],pos[*,wSGOk],sP=sP)
  if countSGNo gt 0 then $
    rad_pri[wSGNo] = periodicDists(gc.groupPos[*,groupParIDs[wSGNo]],pos[*,wSGNo],sP=sP)
  
  if min(rad_pri) lt 0.0 then message,'Error: Unassigned rad.'
  
  ; for vrad: replace coordinates by relative coordinates (radial vectors)
  for i=0,2 do begin
    pos_rel = fltarr(n_elements(sgParIDs))
    
    if countSGOk gt 0 then $
      pos_rel[wSGOk] = reform(pos[i,wSGOk] - subgroupCen[i,sgParIDs[wSGOk]])
    if countSGNo gt 0 then $
      pos_rel[wSGNo] = reform(pos[i,wSGNo] - gc.groupPos[i,groupParIDs[wSGNo]])
    
    correctPeriodicDistVecs, pos_rel, sP=sP
    pos[i,*] = pos_rel
  endfor
  
  pos_rel = !NULL
  
  ; calculate radial velocity
  vrad_pri = fltarr(n_elements(sgParIDs))
  
  if countSGOk gt 0 then $
    vrad_pri[wSGOk] = ((vel[0,wSGOk] - gc.subgroupVel[0,sgParIDs[wSGOk]]) * pos[0,wSGOk] + $
                       (vel[1,wSGOk] - gc.subgroupVel[1,sgParIDs[wSGOk]]) * pos[1,wSGOk] + $
                       (vel[2,wSGOk] - gc.subgroupVel[2,sgParIDs[wSGOk]]) * pos[2,wSGOk]) $
                       / rad_pri[wSGOk]
  if countSGNo gt 0 then $
    vrad_pri[wSGNo] = ((vel[0,wSGNo] - gc.groupVel[0,groupParIDs[wSGNo]]) * pos[0,wSGNo] + $
                       (vel[1,wSGNo] - gc.groupVel[1,groupParIDs[wSGNo]]) * pos[1,wSGNo] + $
                       (vel[2,wSGNo] - gc.groupVel[2,groupParIDs[wSGNo]]) * pos[2,wSGNo]) $
                       / rad_pri[wSGNo]

  ; make (rho,temp,rad) cut
  par_rvir = gc.group_r_crit200[groupParIDs] ; fof parent rvir normalization
  
  wHalo = where(alog10(temp) - sP.galcut_rho*alog10(dens) ge sP.galcut_T and $
                rad_pri/par_rvir ge sP.radcut_rvir and rad_pri/par_rvir le radCutOuter,countHalo)

  if countHalo eq 0 then message,'Warning: Empty halo.'
    
  par_rvir    = !NULL
  temp        = !NULL
  dens        = !NULL
  groupParIDs = !NULL
  gcIDs_ind   = !NULL
  
  ; make subsets matching halo cut
  sgParIDs    = sgParIDs[wHalo]
  rad_pri     = rad_pri[wHalo]
  vrad_pri    = vrad_pri[wHalo]
  
  ; construct galaxyHalo catalog
  fullhaloLen = histogram(sgParIDs,min=0L,max=gc.nSubgroupsTot-1)
  fullhaloOff = ( total([0,fullhaloLen],/int,/cum) )[0:-2] ; remove last
  fullhaloIDs = ids[wHalo]

  if total(fullhaloLen,/int) ne n_elements(fullhaloIDs) then message,'Error: Catalog fail.'
  if n_elements(fullhaloLen) ne n_elements(fullHaloOff) then message,'Error: Catalog fail2.'
  if n_elements(fullHaloLen) ne gc.nSubgroupsTot then message,'Error: Catalog fail3.'
  
  ; save group membership catalog
  r = { fullhaloLen:fullhaloLen, fullhaloOff:fullhaloOff, fullhaloIDs:fullhaloIDs, $
        fullhaloRad:rad_pri, fullhaloVRad:vrad_pri, radCutOuter:radCutOuter }
        
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))+' ['+str(countHalo)+'/'+str(countMatch)+']'    
  
  return, r
  
end

; galHaloCatRepParentIDs(): for the galaxyHalo catalog, replicate the list of ordered parent IDs such that
;                           the return array is the same size as the number of gas particles with
;                           each element the id of its parent subgroup/galaxy group/groupmember group

; gcIDList : return only a replicated parent ID list of the specified subgroups in the groupcat
; child_counts: return a replicated list for tracers which have multiplicity inside each gas cell
;               or star particle as specified in child_counts

function galHaloCatRepParentIDs, galHaloCat=galHaloCat, gcIDList=gcIDList, child_counts=child_counts

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if not keyword_set(gcIDList) then $
    gcIDList = lindgen(n_elements(galHaloCat.fullhaloLen)) ; id list to process set to all
    
  if keyword_set(child_counts) then $
    if n_elements(child_counts) ne total(galHaloCat.fullhaloLen[gcIDList],/int) then $
       message,'Error: Child_counts should have same size as galHaloCat gcIDList subset.'
         
  if keyword_set(child_counts) then begin
    if total(child_counts,/int) gt 2e9 then stop ; consider lon64/removing /int
  endif else begin
    child_counts = lonarr(total(galHaloCat.fullhaloLen[gcIDList],/int))+1
  endelse
    
  r = lonarr(total(child_counts,/int))
    
  offset   = 0L
  offset_c = 0L
    
  ; loop over all subhalos
  foreach gcID,gcIDList do begin
    
      if galHaloCat.fullhaloLen[gcID] gt 0 then begin
        tot_children = total(child_counts[offset_c:offset_c+galHaloCat.fullhaloLen[gcID]-1],/int)
        ind_end = offset+tot_children-1

        if tot_children gt 0 then r[offset:ind_end] = cmreplicate(gcID,tot_children)
        offset += tot_children
        offset_c += galHaloCat.fullhaloLen[gcID]
      endif
        
  endforeach

  return,r
end

; galHaloCatParentProperties: calculate some property of the parent for every gas elements
;                             in the galaxyHalo catalog at some snapshot
; virTemp=1 : virial temperature
; mass=1    : total mass (from catalog, dm+baryon)
; rVir=1    : virial radius (r_200 critical)

function galHaloCatParentProperties, sP=sP, virTemp=virTemp, mass=mass, rVir=rVir, vCirc=vCirc

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function galaxyHaloCat, snapNumToRedshift, codeMassToLogMsun

  ; load group catalog for masses
  gc = loadGroupCat(sP=sP,/skipIDs)

  ; load galaxy catalog
  galHaloCat = galaxyHaloCat(sP=sP)

  ; load parent indices for each member
  gcInd = galHaloCatRepParentIDs(galHaloCat=galHaloCat)
   
  ; arrays
  r = fltarr(n_elements(galHaloCat.fullhaloIDs))
  
  ; masses (log msun)
  if keyword_set(mass) then begin
    r = gc.subgroupMass[gcInd]
    r = codeMassToLogMsun(r)
  endif

  ; calculate virial temperatures (K)
  if keyword_set(virTemp) then begin
    redshift = snapNumToRedshift(sP=sP)
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

; galHaloCatINDList(): return a list of indices into the galaxyHalo catalog for a subset of the
;                      members defined by the subgroup ID list gcIDList
;                  
; child_counts: return a replicated list for tracers which have multiplicity inside each gas cell
;               as specified in child_counts

function galHaloCatINDList, sP=sP, galHaloCat=galHaloCat, gcIDList=gcIDList, child_counts=child_counts

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function galaxyCat

  ; load galaxy cat if necessary
  if not keyword_set(galHaloCat) then begin
    if not keyword_set(sP) then begin
      print,'Error: Must specify galcat or sP.' & stop
    endif
    galHaloCat = galaxyHaloCat(sP=sP)
  endif
  
  if max(galHaloCat.fullhaloLen+galHaloCat.fullhaloOff) gt 2e9 then stop ; change to lon64arr
  
  ; make mask for requested indices
  gcIDMask = bytarr(n_elements(galHaloCat.fullhaloLen))
  if keyword_set(gcIDList) then gcIDMask[gcIDList] = 1B  
  if ~keyword_set(gcIDList) then begin
    gcIDMask[*] = 1B
    gcIDList = lindgen(n_elements(galHaloCat.fullhaloLen))
  endif
  
  ; normal indices return
  r = ulonarr(total(galHaloCat.fullhaloLen[gcIDList],/int))
  
  offset = 0L
  
  ; (1) make list for gas cells/particles
  foreach gcID, gcIDList do begin
    ; galaxy
    if (galHaloCat.fullhaloLen[gcID] gt 0) then begin
      hInds = ulindgen(galHaloCat.fullhaloLen[gcID]) + galHaloCat.fullhaloOff[gcID]
      r[ offset:offset+galHaloCat.fullhaloLen[gcID]-1 ] = hInds
      offset += galHaloCat.fullhaloLen[gcID]
    endif
  endforeach
  
  ; (2) make list including child counts if requested
  if ~keyword_set(child_counts) then return,r
  
  if n_elements(child_counts) ne total(galHaloCat.fullhaloLen,/int) then $
     message,'Error: Child_counts should have same size as full galHaloCat subset.'

  if total(child_counts,/int) gt 2e9 then stop ; consider lon64/removing /int

  rcc = ulonarr(total(child_counts[r],/int))
       
  offset     = 0L
  offset_all = 0L
  
  for gcID=0UL,n_elements(galHaloCat.fullhaloLen)-1 do begin
    ; galaxy
    if galcat.galaxyLen[gcID] gt 0 then begin
      ; calculate total number of children in this subgroup
      tot_children = total(child_counts[galHaloCat.fullhaloOff[gcID]:$
                           galHaloCat.fullhaloOff[gcID]+galHaloCat.fullhaloLen[gcID]-1],/int)
                           
      ; add indices only for specified galaxy IDs
      if gcIDMask[gcID] eq 1B and tot_children gt 0 then begin
        ; calculate place and store indices
        hInds = ulindgen(tot_children) + offset_all
        rcc[ offset:offset+tot_children-1 ] = hInds
  
        offset += tot_children
      endif
      
      ; add child counts to indice array offset
      offset_all  += tot_children
    endif      
                  
    
  endfor
  
  if offset ne n_elements(rcc) then message,'Error.'

  return,rcc
  
end
