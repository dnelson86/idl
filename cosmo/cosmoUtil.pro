; cosmoUtil.pro
; cosmological simulations - utility functions
; dnelson apr.2012

; redshiftToSnapNum(): convert redshift to the nearest snapshot number

function redshiftToSnapNum, redshiftList, sP=sP, verbose=verbose

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if not keyword_set(verbose) then verbose = 0
  
  saveFileName = sP.derivPath + sP.savPrefix + '_snapnum.redshift.sav'

  if not (file_test(saveFileName)) then begin

    nSnaps = 400 ; maximum
    
    redshifts = fltarr(nSnaps) - 1
    times     = fltarr(nSnaps) - 1
    
    for m=0,nSnaps-1 do begin
      ; format filename
      ext = str(string(m,format='(I3.3)'))
      f = sP.simPath + 'snapdir_' + ext + '/snap_' + ext + '.0.hdf5'
      
      ; single file per group catalog
      if (not file_test(f)) then $
        f = sP.simPath + 'snap_' + ext + '.hdf5'
        
      ; single groupordered file per group catalog
      if (not file_test(f)) then $
        f = sP.simPath + 'snap-groupordered_' + ext + '.hdf5'
      
      ; if file doesn't exist yet, skip (e.g. every other file deleted)
      if (not file_test(f)) then continue
    
      ; load hdf5 header and save time+redshift
      fileID   = h5f_open(f)
      headerID = h5g_open(fileID,"Header")
      
      redshifts[m] = h5a_read(h5a_open_name(headerID,"Redshift"))
      times[m]     = h5a_read(h5a_open_name(headerID,"Time"))
      
      h5g_close, headerID
      h5f_close, fileID
    endfor
  
    ; save/restore
    save,nSnaps,redshifts,times,filename=saveFileName
  endif else begin
    restore,saveFileName
  endelse

  ; for -1 return !NULL
  if (total(redshiftList) eq -1) then return,!NULL
  
  ; return array
  snapNum = intarr(n_elements(redshiftList))
  
  foreach redshift,redshiftList,i do begin

    dists = abs(redshifts - redshift)
    w = where(dists eq min(dists),count)
    if count eq 2 then w = w[0]
    if count eq 0 then stop
    if min(dists) gt 0.1 then print,'Warning! Snapshot selected with redshift error = ',min(dists)
    snapNum[i] = w[0]
    
    if (verbose) then $
      print,'Found nearest snapshot to z = ' + str(redshift) + ' (num = ' + str(snapNum[i]) + ')'
  endforeach

  if (n_elements(snapNum) eq 1) then snapNum = snapNum[0]

  return,snapNum
end

; mergerTreeChild(): construct a Child array given a Parent array
; 
; childPrev : if specified, compose the two mappings such that the returned Child points not at the
;             immediate children of Parent but to a prior child of Parent (say, the first)

function mergerTreeChild, Parent, ChildPrev=ChildPrev

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  Child = lonarr(max(Parent)+1)-1
  for i=0UL,n_elements(Parent)-1L do if Parent[i] ne -1 then Child[Parent[i]] = i

  ; compose with previous Child mapping?
  if n_elements(ChildPrev) gt 0 then begin
    ChildNew = lonarr(n_elements(Child))-1
    
    ; for each immediate Child, replace Cur->Child pointer by Cur->Child->PrevChild pointer
    for i=0UL,n_elements(Child)-1 do begin
      ; discard pointers to children greater than the maximum in ChildPrev, which can happen when
      ; the Child pointed to is not part of the halo set being tracked (the composed mapping is
      ; undefined for this Child)
      if Child[i] lt n_elements(ChildPrev) then $
        if Child[i] ne -1 then ChildNew[i] = ChildPrev[Child[i]] ; only follow valid Child links
    endfor
    
    return, ChildNew
  endif

  return, Child
end

; gcPriChildInd(): get the subgroup index of the primary subgroup of a given group, or -1 if none

function gcPriChildInd, gc=gc, haloID=haloID

  if (n_elements(gc) eq 0 or n_elements(haloID) eq 0) then stop

  gcInd = min(where(gc.subgroupGrNr eq haloID,count))
  fsInd = gc.groupFirstSub[haloID]
  
  ; skip if group has no subgroup to obtain this value from
  if (count eq 0) then return,-1
  
  ; double-check (failing either should indicate a problem in the group catalogs, right?)
  if (gcInd ne fsInd) then return,-1 ;stop
  if (gc.subgroupGrNr[gcInd] ne haloID) then return,-1 ;stop

  return,gcInd

end

; gcIDList(): return a sub-list of subgroup IDs from a given group catalog gc
; 
; select: one of the following
;   'pri' : members of the first subgroup of each group only ("background"/"main subhalo"/"halo")
;   'sec' : members of the non-first subgroups of each group only ("satellites"/"subhalos")
;   'all' : all subgroups

function gcIDList, sP=sP, gc=gc, select=select

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function loadGroupCat

  if (not keyword_set(select)) then begin
    print,'Error: Must specify select.' & stop
  endif
    
  ; load galaxy cat if necessary
  if not keyword_set(gc) then begin
    if not keyword_set(sP) then begin
      print,'Error: Must specify gc or sP.' & stop
    endif
    gc = loadGroupCat(sP=sP,/skipIDs)
  endif

  ; require a minimum number of gas particles in (subfind) group to include
  ;if sP.minNumGasPart gt 0 then stop ; don't always have access to sP here

  prevGrNr   = -1
  valGCids   = []

  if (select eq 'pri') then begin
  
    ; "background"/"main subhalos"/"halos" only
    for i=0UL,n_elements(gc.subgroupLen)-1 do begin
      if (gc.subgroupGrnr[i] eq prevGrNr) then begin
        prevGrNr = gc.subgroupGrnr[i]
      endif else begin
      
        if (keyword_set(minNumPart)) then begin
          if (gc.subgroupLen[i] ge minNumPartVal) then $
            valGCids = [valGCids,i]
        endif else begin
          valGCids = [valGCids,i]
        endelse ;minNumPart
        
        prevGrNr = gc.subgroupGrnr[i]
      endelse
    endfor
    return,valGCids
  endif
  
  if (select eq 'sec') then begin
  
    ; "satellites"/"subhalos" only
    for i=0UL,n_elements(gc.subgroupLen)-1 do begin
      if (gc.subgroupGrnr[i] ne prevGrNr) then begin
        prevGrNr = gc.subgroupGrnr[i]
      endif else begin
      
        if (keyword_set(minNumPart)) then begin
          if (gc.subgroupLen[i] ge minNumPartVal) then $
            valGCids = [valGCids,i]
        endif else begin
          valGCids = [valGCids,i]
        endelse ;minNumPart
        
      endelse
    endfor
    return,valGCids
  endif
  
  if (select eq 'all') then begin
  
    ; both primary and secondary
    for i=0UL,n_elements(gc.subgroupLen)-1 do begin
    
      if (keyword_set(minNumPart)) then begin
        if (gc.subgroupLen[i] ge minNumPartVal) then $
          valGCids = [valGCids,i]
      endif else begin
        valGCids = [valGCids,i]
      endelse ;minNumPart
      
    endfor
    return,valGCids
  endif
  
  message,'Error: Unrecognized select in gcIDList?'

end

; gcPIDList(): return a list of member particle IDs from a given group catalog gc
;
; select: one of the following
;   'pri' : members of the first subgroup of each group only ("background"/"main subhalo"/"halo")
;   'sec' : members of the non-first subgroups of each group only ("satellites"/"subhalos")
;   'all' : all subgroups
; valGCids: return member ids from only these groups
;
; partType=N     : use subgroup type offset table to return only IDs of a specific particle type
; partType='all' : return all particle IDs in that subgroup regardless of type

function gcPIDList, gc=gc, select=select, valGCids=valGCids, partType=PT

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if (n_elements(gc) eq 0 or (n_elements(select) eq 0 and n_elements(valGCids) eq 0) or $
      n_elements(PT) eq 0) then message,'Error: gcPIDList: Bad inputs.'

  ; get list of appropriate group ids
  if n_elements(valGCids) eq 0 then valGCids = gcIDList(gc=gc,select=select)
  
  ; make list of particle ids in these groups
  start = 0L
  
  partType = strlowcase(string(PT)) ; so we don't change the input
  
  if strcmp(partType,'all') then begin
    subgroupPIDs = lonarr(total(gc.subGroupLen[valGCids],/int))
    
    foreach gcID, valGCids do begin
      ; select particle IDs in subgroup
      subgroupPIDs[start:start+gc.subGroupLen[gcID]-1] = $
        gc.IDs[gc.subGroupOffset[gcID] : gc.subGroupOffset[gcID] + gc.subGroupLen[gcID] - 1]
      start += gc.subGroupLen[gcID]
    endforeach
    
    if min(subgroupPIDs) lt 0 then stop ; check 32 bit long overflow
    return, subgroupPIDs
  endif
  
  ; convert partType to number if it is a name
  partType = partTypeNum(partType)

  ; check if this particle type is present in the subgroup selection
  if total(gc.subgroupLenType[partType,valGCids] gt 0) then begin

    if max(gc.IDs) gt 2e9 then stop ; change to lon64arr
    subgroupPIDs = lonarr(total(gc.subGroupLenType[partType,valGCids],/int))

    ; store particle IDs of this type from each subgroup
    foreach gcID, valGCids do begin
      ; may be none
      if gc.subgroupLenType[partType,gcID] gt 0 then begin
        subgroupPIDs[start:start+gc.subGroupLenType[partType,gcID]-1] = $
          gc.IDs[gc.subGroupOffsetType[partType,gcID] : gc.subGroupOffsetType[partType,gcID] + $
          gc.subGroupLenType[partType,gcID] - 1]
        start += gc.subGroupLenType[partType,gcID]
      endif
    endforeach
    
    if min(subgroupPIDs) lt 0 then stop ; check 32 bit long overflow
    return, subgroupPIDs    
  endif
  
  print,'Warning! Empty gcPIDList return.'
  return,[]

end

; massTargetToHaloID(): return the halo ID (subgroupInd) nearest to the target mass in log(Msun)

function massTargetToHaloID, hMassTargets, sP=sP, verbose=verbose

  if ~keyword_set(sP) then message,'Error'
  
  ; load group catalog and calculate log(M) masses of all subgroups of the requested type
  gc = loadGroupCat(sP=sP,/skipIDs)
  priSGIDs = gcIDList(gc=gc,select='pri')
  hMasses = codeMassToLogMsun(gc.subgroupMass[priSGIDs])

  ; locate nearest masses to requested masses in a consistent way
  hInds = value_locate(hMasses,hMassTargets)
  w = where(hInds eq -1,count)
  if count gt 0 then hInds[w] = 0 ; largest mass halo available
 
  ; convert indices in hMasses to sgIDs
  subgroupIDs = priSGIDs[hInds] 

  if keyword_set(verbose) then $
      print,'selected halo ind ['+str(hInd)+'] sgID ['+str(subgroupID)+'] mass = '+$
        string(hMasses[hInd],format='(f5.2)')

  return,subgroupIDs
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
;               as specified in child_counts

function galCatRepParentIDs, galcat=galcat, priParentIDs=priParentIDs, $
                             gcIDList=gcIDList, child_counts=child_counts

  compile_opt idl2, hidden, strictarr, strictarrsubs

    if not keyword_set(priParentIDs) then $
      priParentIDs = lindgen(n_elements(galcat.galaxyLen)) ; valid id list set to all
    
    if not keyword_set(gcIDList) then $
      gcIDList = lindgen(n_elements(galcat.galaxyLen)) ; id list to process set to all
    
    if keyword_set(child_counts) then $
      if n_elements(child_counts.gal) ne total(galcat.galaxyLen[gcIDList],/int) or $
         n_elements(child_counts.gmem) ne total(galcat.groupmemLen[gcIDList],/int) then $
         message,'Error: Child_counts gal/gmem should have same size as galcat gcIDList subset.'
         
    if keyword_set(child_counts) then begin
      if total(child_counts.gal,/int) gt 2e9 then stop ; consider lon64/removing /int
      if total(child_counts.gmem,/int) gt 2e9 then stop
    endif else begin
      child_counts = { gal  : lonarr(total(galcat.galaxyLen[gcIDList],/int))+1   ,$
                       gmem : lonarr(total(galcat.groupmemLen[gcIDList],/int))+1  }
    endelse
    
    ; old version without child_counts
    ;r = { gal  : lonarr(total(galcat.galaxyLen[gcIDList],/int))  ,$
    ;      gmem : lonarr(total(galcat.groupmemLen[gcIDList],/int)) }
    
    r = { gal  : lonarr(total(child_counts.gal,/int))  ,$
          gmem : lonarr(total(child_counts.gmem,/int)) }
    
    offsetGal  = 0L
    offsetGmem = 0L
    
    offsetGal_c  = 0L
    offsetGmem_c = 0L
    
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
          tot_children = total(child_counts.gal[offsetGal_c:offsetGal_c+galcat.galaxyLen[gcID]-1],/int)
          gal_ind_end = offsetGal+tot_children-1

          if tot_children gt 0 then r.gal[offsetGal:gal_ind_end] = cmreplicate(parID,tot_children)
          offsetGal += tot_children
          offsetGal_c += galcat.galaxyLen[gcID]
          
          ;r.gal[offsetGal:offsetGal+galcat.galaxyLen[gcID]-1] = $
          ;  cmreplicate(parID,galcat.galaxyLen[gcID])
          ;offsetGal += galcat.galaxyLen[gcID]
        endif
        
        ; group members
        if galcat.groupmemLen[gcID] gt 0 then begin
          tot_children = total(child_counts.gmem[offsetGmem_c:offsetGmem_c+galcat.groupmemLen[gcID]-1],/int)
          gmem_ind_end = offsetGmem+tot_children-1

          if tot_children gt 0 then r.gmem[offsetGmem:gmem_ind_end] = cmreplicate(parID,tot_children)
          offsetGmem += tot_children
          offsetGmem_c += galcat.groupmemLen[gcID]
          ;r.gmem[offsetGamem:offsetGmem+galcat.groupmemLen[gcID]-1] = $
          ;  cmreplicate(parID,galcat.groupmemLen[gcID])
          ;offsetGmem += galcat.groupmemLen[gcID]
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

function galCatParentProperties, sP=sP, virTemp=virTemp, mass=mass, rVir=rVir, parNorm=parNorm

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
  gal  = fltarr(n_elements(galcat.galaxyIDs))
  gmem = fltarr(n_elements(galcat.groupmemIDs))
  
  ; masses (log msun)
  if keyword_set(mass) then begin
    gal  = gc.subgroupMass[gcInd.gal]
    gmem = gc.subgroupMass[gcInd.gmem]
    
    gal  = codeMassToLogMsun(gal)
    gmem = codeMassToLogMsun(gmem)
  endif

  ; calculate virial temperatures (K)
  if keyword_set(virTemp) then begin
    gal  = gc.subgroupMass[gcInd.gal]
    gmem = gc.subgroupMass[gcInd.gmem]
    
    redshift = snapNumToRedshift(sP=sP)
  
    gal  = alog10(codeMassToVirTemp(gal,sP=sP))
    gmem = alog10(codeMassToVirTemp(gmem,sP=sP))
  endif
  
  if keyword_set(rVir) then begin
    gal  = gc.subgroupGrnr[gcInd.gal]
    gmem = gc.subgroupGrnr[gcInd.gmem]
    
    gal  = gc.group_r_crit200[gal]
    gmem = gc.group_r_crit200[gmem]
  endif

  r = {gal:gal,gmem:gmem}
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
  
  ; make mask for requested indices
  gcIDMask = bytarr(n_elements(galcat.galaxyLen))
  if keyword_set(gcIDList) then gcIDMask[gcIDList] = 1B  
  if ~keyword_set(gcIDList) then gcIDMask[*] = 1B
  
  ; normal indices return
  r = {gal  : ulonarr(total(galcat.galaxyLen[gcIDList],/int)) ,$
       gmem : ulonarr(total(galcat.groupmemLen[gcIDList],/int)) }  
  
  offsetGal  = 0L
  offsetGmem = 0L  
  
  ; (1) make list for gas cells/particles
  foreach gcID, gcIDList do begin
    ; galaxy
    if (galcat.galaxyLen[gcID] gt 0) then begin
      galInds    = ulindgen(galcat.galaxyLen[gcID]) + galcat.galaxyOff[gcID]
      r.gal[offsetGal:offsetGal+galcat.galaxyLen[gcID]-1] = galInds
      offsetGal += galcat.galaxyLen[gcID]
    endif
    
    ; group member
    if (galcat.groupmemLen[gcID] gt 0) then begin
      gmemInds    = ulindgen(galcat.groupmemLen[gcID]) + galcat.groupmemOff[gcID]
      r.gmem[offsetGmem:offsetGmem+galcat.groupmemLen[gcID]-1] = gmemInds
      offsetGmem += galcat.groupmemLen[gcID]
    endif
  endforeach
  
  ; (2) make list including child counts if requested
  if ~keyword_set(child_counts) then return,r
  
  if n_elements(child_counts.gal) ne total(galcat.galaxyLen,/int) or $
     n_elements(child_counts.gmem) ne total(galcat.groupmemLen,/int) then $
     message,'Error: Child_counts gal/gmem should have same size as full galcat subset.'

  if total(child_counts.gal,/int) gt 2e9 then stop ; consider lon64/removing /int
  if total(child_counts.gmem,/int) gt 2e9 then stop

  rcc = { gal  : ulonarr(total(child_counts.gal[r.gal],/int))  ,$
          gmem : ulonarr(total(child_counts.gmem[r.gmem],/int)) }
       
  offsetGal  = 0L
  offsetGmem = 0L
  
  offsetGal_all  = 0UL
  offsetGmem_all = 0UL
  
  for gcID=0UL,n_elements(galcat.galaxyLen)-1 do begin
    ; galaxy
    if galcat.galaxyLen[gcID] gt 0 then begin
      ; calculate total number of children in this subgroup
      tot_children_gal  = total(child_counts.gal[galcat.galaxyOff[gcID]:$
                                                 galcat.galaxyOff[gcID]+galcat.galaxyLen[gcID]-1],/int)
      ; add indices only for specified galaxy IDs
      if gcIDMask[gcID] then begin
        ; calculate place and store indices
        galInds = ulindgen(tot_children_gal) + offsetGal_all
        rcc.gal[offsetGal:offsetGal+tot_children_gal-1] = galInds
  
        offsetGal += tot_children_gal
      endif
      
      ; add child counts to indice array offset
      offsetGal_all  += tot_children_gal
    endif      
                  
    ; group member
    if galcat.groupmemLen[gcID] gt 0 then begin
      ; calculate total number of children in this subgroup             
      tot_children_gmem = total(child_counts.gmem[galcat.groupmemOff[gcID]:$
                                                  galcat.groupmemOff[gcID]+galcat.groupmemLen[gcID]-1],/int)
            
      ; add indices only for specified group member IDs
      if gcIDMask[gcID] then begin
        ; calculate place and store indices
        gmemInds     = ulindgen(tot_children_gmem) + offsetGmem_all
        rcc.gmem[offsetGmem:offsetGmem+tot_children_gmem-1] = gmemInds
        
        offsetGmem += tot_children_gmem
      endif
      
      offsetGmem_all += tot_children_gmem
    endif
    
    ;if gcID lt 100 then $
    ;print,gcID,gcIDMask[gcID],offsetGal,offsetGal_all,offsetGmem,offsetGmem_all
    
  endfor

    ; this code isn't useful here but is the variable-count replicate (better verify though):
    ;galInds = ulonarr(tot_children)
    ;cc = total(child_counts.gal[offsetGal_c:offsetGal_c+galcat.galaxyLen[gcID]-1],/int,/cum)
    ;if n_elements(cc) gt 1 then begin
    ;  if cc[n_elements(cc)-1] eq cc[n_elements(cc)-2] then cc=cc[0:n_elements(cc)-2]
    ;  galInds[cc[0:n_elements(cc)-2]] = 1
    ;endif
    ;galInds = total(galInds,/cum,/int)
    ;galInds += galcat.galaxyOff[gcID]  
  
  return,rcc
  
end

; subgroupPosByMostBoundID(): compute a "best" center position in space for all subfind groups by 
;                             using the position of the most bound particles, whose IDs are stored in 
;                             the group catalog but without knowing the particle types we have to 
;                             load all gas+dm+stars particle positions

function subgroupPosByMostBoundID, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function loadSnapshotSubset, loadGroupCat
  
  ; save/restore
  saveFilename = sP.derivPath+'sgCen.mbID/sgCen.mbID.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.sav'
                 
  if file_test(saveFilename) then begin
    restore,saveFilename
  endif else begin 
  
    gc = loadGroupCat(sP=sP,/skipIDs)
  
    groupCen = fltarr(3,gc.nSubgroupsTot)
    
    ; load gas ids and pos, find matches
    ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    match,gc.subgroupIdMostBound,ids,gc_ind,ids_ind,count=count1,/sort
    
    ids_ind = ids_ind[sort(gc_ind)]
    gc_ind  = gc_ind[sort(gc_ind)]
    
    pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
    
    groupCen[*,gc_ind] = pos[*,ids_ind]
    
    ; load stars ids and pos, find matches
    ids = loadSnapshotSubset(sP=sP,partType='stars',field='ids')
    
    ; skip stars for non-radiative runs or early times
    if (n_elements(ids) ne 1) then begin
      match,gc.subgroupIdMostBound,ids,gc_ind,ids_ind,count=count2,/sort
      
      ids_ind = ids_ind[sort(gc_ind)]
      gc_ind  = gc_ind[sort(gc_ind)]
      
      pos = loadSnapshotSubset(sP=sP,partType='stars',field='pos')
      
      groupCen[*,gc_ind] = pos[*,ids_ind]
    endif else begin
      count2 = 0
    endelse
    
    ; load dm ids and pos, find matches
    ids = loadSnapshotSubset(sP=sP,partType='dm',field='ids')
    match,gc.subgroupIdMostBound,ids,gc_ind,ids_ind,count=count3,/sort
    
    ids_ind = ids_ind[sort(gc_ind)]
    gc_ind  = gc_ind[sort(gc_ind)]
    
    pos = loadSnapshotSubset(sP=sP,partType='dm',field='pos')
    
    groupCen[*,gc_ind] = pos[*,ids_ind]
  
    if ((count1 + count2 + count3) ne gc.nSubgroupsTot) then begin
      print,'ERROR: Failed to find all most bound IDs.' & stop
    endif
  
    ; save
    save,groupCen,filename=saveFilename
  endelse
  
  return, groupCen
end

; correctPeriodicDistVecs(): enforce periodic B.C. for distance vecotrs (effectively component by 
;                            component), input vecs in format fltarr[3,n]

pro correctPeriodicDistVecs, vecs, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  w = where(vecs gt sP.boxSize*0.5,count)
  if (count ne 0) then $
    vecs[w] = vecs[w] - sP.boxSize
    
  w = where(vecs lt -sP.boxSize*0.5,count)
  if (count ne 0) then $
    vecs[w] = sP.boxSize + vecs[w]

end

; periodicDists(): calculate distances correctly taking into account periodic B.C.
; 
; if pt is one point: distance from pt to all vecs
; if pt is several points: distance from each pt to each vec (must have same number of points)

function periodicDists, pt, vecs, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  nDimsPt = (size(pt))[0]

  if ( ((size(vecs))[0] ne 1 and (size(vecs))[0] ne 2) or $
       (size(vecs))[1] ne 3) then message,'Error: vecs not in expected shape'
  if (nDimsPt ne 1 and nDimsPt ne 2) then message,'Error: something is wrong'
  if not keyword_set(sP) then message,'Error: need sP for boxSize'

  ; distances from one point to a vector of other points [3,n]
  if (nDimsPt eq 1) then begin
    xDist = vecs[0,*]-pt[0]
    yDist = vecs[1,*]-pt[1]
    zDist = vecs[2,*]-pt[2]
  endif
  
  ; distances from a vector of points [3,n] to another vector of other points [3,n]
  if (nDimsPt eq 2) then begin
    xDist = vecs[0,*]-pt[0,*]
    yDist = vecs[1,*]-pt[1,*]
    zDist = vecs[2,*]-pt[2,*]
  endif
  
  correctPeriodicDistVecs, xDist, sP=sP
  correctPeriodicDistVecs, yDist, sP=sP
  correctPeriodicDistVecs, zDist, sP=sP
  
  dists = reform( sqrt( xDist*xDist + yDist*yDist + zDist*zDist ) )
  
  return, dists

end

; periodicPairwiseDists(): calculate pairwise distances between all 3D points, correctly
;                          taking into account periodic B.C.

function periodicPairwiseDists, pts, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; construct indices for vector computation of pairwise distances
  dims = size(pts,/dimensions)
  
  num = dims[1]*(dims[1]-1)/2
  
  ii = 0L
  index0 = lindgen(dims[1]-1) + 1
  index1 = lonarr(num,/nozero)
  index2 = lonarr(num,/nozero)
  
  for i=0,dims[1]-2 do begin
    n1 = dims[1] - (i+1)
    index1[ii:ii+n1-1] = i
    index2[ii] = index0[0:n1-1] + i
    ii += n1
  endfor
  
  ; component wise difference
  xDist = pts[0,index1] - pts[0,index2]
  yDist = pts[1,index1] - pts[1,index2]
  zDist = pts[2,index1] - pts[2,index2]
  
  ; account for periodic distance function
  correctPeriodicDistVecs, xDist, sP=sP
  correctPeriodicDistVecs, yDist, sP=sP
  correctPeriodicDistVecs, zDist, sP=sP
  
  dists = reform( sqrt( xDist*xDist + yDist*yDist + zDist*zDist ) )

  return,dists

end

; groupCenterPosByIterativeCM(): compute a "better" center position in space for all FoF groups by
;                                iteratively searching for the center of mass

function groupCenterPosByIterativeCM, sP=sP, gc=gc, haloIDs=haloIDs

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function loadSnapshotSubset, loadSnapshotHeader

  if (n_elements(haloIDs) eq 0) then stop

  ; config
  radStartStop = [1.0,0.2] ; fractions of r_vir
  radStep      = 0.01 ; decrease search radius by this fraction each step

  nSteps = fix((radStartStop[0]-radStartStop[1]) / radStep + 1)
  radSteps = linspace(radStartStop[0],radStartStop[1],nSteps)
  
  iterDM  = fltarr(3,n_elements(haloIDs))
  iterGAS = fltarr(3,n_elements(haloIDs))

  ; load
  h  = loadSnapshotHeader(sP=sP)
  
  pos_gas   = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  pos_dm    = loadSnapshotSubset(sP=sP,partType='dm',field='pos')
  
  gas_mass   = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
  dm_mass    = h.massTable[1]  
  
  ; locate halo
  foreach haloID,haloIDs,j do begin
    haloPos = gc.groupPos[*,haloID]

    ; NaNs in gc?
    if (not finite(haloPos[0])) then begin
      iterDM[*,j]  = [-1,0,0]
      iterGAS[*,j] = [-1,0,0]
      continue
    endif

    ;haloRad = 100.0
    ;if tag_exist(gc,'group_r_crit200') then stop ; catch
    haloRad = gc.group_r_crit200[haloID]

    ; calculate radii and make radial cut
    rad = periodicDists(haloPos,pos_gas,sP=sP)
    
    gas_ind = where(rad le haloRad*radStartStop[0],gas_count)
    pos_gas_sub = pos_gas[*,gas_ind]

    rad = periodicDists(haloPos,pos_dm,sP=sP)

    dm_ind = where(rad le haloRad*radStartStop[0],dm_count)
    pos_dm_sub = pos_dm[*,dm_ind]

    rad = !NULL
    
    ;print,'('+str(haloID)+') Found ['+str(gas_count)+'] gas and ['+str(dm_count)+'] dm  stars inside cut.'
                     
    ; subselect gas masses
    gas_mass_sub = gas_mass[gas_ind]
  
    cm_gas = fltarr(3,nSteps) ;xyz
    cm_dm  = fltarr(3,nSteps) ;xyz
    
    cm_diffs = fltarr(2,nSteps) ;gas,dm
    cm_drift = fltarr(2,nSteps) ;gas,dm
    
    ; compute a center of mass (gas) and compare to haloPos
    foreach radCut,radSteps,k do begin

      ; calculate radii based on last CM
      rad_gas = periodicDists(cm_gas[*,k-1>0],pos_gas_sub,sP=sP)
      rad_dm  = periodicDists(cm_dm[*,k-1>0],pos_dm_sub,sP=sP)

      ; make selection inside this radial cut
      w_gas = where(rad_gas le radCut*haloRad,count_gas)
      w_dm  = where(rad_dm  le radCut*haloRad,count_dm)
      
      ;if (count_dm lt 100) then continue
      
      ; re-calculate center of mass and distance from FoF center
      for i=0,2 do cm_gas[i,k] = total(gas_mass_sub[w_gas] * pos_gas_sub[i,w_gas]) / $
                                 total(gas_mass_sub[w_gas])
      for i=0,2 do cm_dm[i,k]  = mean(pos_dm_sub[i,w_dm])
      
      cm_diffs[*,k] = [sqrt(total((haloPos-cm_gas[*,k])^2.0)),sqrt(total((haloPos-cm_dm[*,k])^2.0))]
      cm_drift[0,k] = sqrt(total((cm_gas[*,k-1>0]-cm_gas[*,k])^2.0))
      cm_drift[1,k] = sqrt(total((cm_dm[*,k-1>0]-cm_dm[*,k])^2.0))
      
      ;debug:
      ;print,k,radCut,count_gas,count_dm
      ;print,'haloPos: ',haloPos
      ;print,'gasCM:   ',cm_gas[*,k]
      ;print,'dmCM:    ',cm_dm[*,k]
      ;print,'cmDiffs:  ',cm_diffs[*,k]/haloRad
    endforeach ;radSteps

    ; save final
    iterDM[*,j]  = cm_dm[*,nSteps-1]
    iterGAS[*,j] = cm_gas[*,nSteps-1]
  endforeach ;haloIDs
  
  ; debug
  ;start_PS,'itercm_'+str(haloID)+'.eps'
  ;  plotsym,0,/fill
  ;  fsc_plot,[0],[0],/nodata,xrange=radStartStop,yrange=[0.0002,4.0],/xs,/ys,/ylog,$
  ;    xtitle="radial search cut [r/"+textoidl("r_{vir}")+"]",$
  ;    ytitle="difference [r/"+textoidl("r_{vir}")+"]",$
  ;    title="iterative CoM convergence (haloID="+str(haloID)+")"
  ;  
  ;  fsc_plot,radSteps,cm_diffs[0,*]/haloRad,line=0,color=getColor(1),/overplot
  ;  fsc_plot,radSteps,cm_diffs[1,*]/haloRad,line=0,color=getColor(2),/overplot
  ;  fsc_plot,radSteps[1:*],cm_drift[0,1:*]/haloRad,line=0,color=getColor(3),/overplot
  ;  fsc_plot,radSteps[1:*],cm_drift[1,1:*]/haloRad,line=0,color=getColor(7),/overplot
  ;  
  ;  legend,['gas FoF delta','dm FoF delta','gas drift','dm drift'],textcolors=getColor([1,2,3,7],/name),$
  ;    box=0,/left,/top
  ;end_PS
  
  r = {iterDM:iterDM,iterGAS:iterGAS}
  return,r
end

; findMatchedHalos(): cross-match halos between two group catalogs
; (different runs or resolutions of the same ICs where the same halos will be in nearly the same spots)

function findMatchedHalos, sP1=sP1, sP2=sP2

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function loadGroupCat

  ; config
  distTol = 40.0 ; kpc
  massTol = 0.4  ; 40%
  
  ; save/restore
  saveFilename = sP1.derivPath + 'matchCat_'+sP2.run+'_'+str(sP2.res)+'-'+str(sP1.snap)+'_m'+$
                 str(fix(massTol*10))+'_d'+str(fix(distTol))+'.sav'
    
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif  
  
  ; load group catalogs
  gc1 = loadGroupCat(sP=sP1,/skipIDs)
  gc2 = loadGroupCat(sP=sP2,/skipIDs)

  ; matched Ind in groupcat2 for each FoF halo in groupcat1
  matchedInds = lonarr(gc1.nGroupsTot) - 1
  massDiffs   = fltarr(gc1.nGroupsTot)
  posDiffs    = fltarr(gc1.nGroupsTot)
  
  targetPos   = fltarr(3,gc2.nGroupsTot)
 
  ; load most bound particle ID positions for subfind groups
  idMBCM1 = subgroupPosByMostBoundID(sP=sP1)
  idMBCM2 = subgroupPosByMostBoundID(sP=sP2)
  
  ; fill target positions with mbID positions for each first subgroup of each FoF group in groupcat2
  for i=0,gc2.nGroupsTot-1 do begin
    gcInd = gcPriChildInd(gc=gc2,haloID=i)
    targetPos[*,i] = idMBCM2[*,gcInd]
  endfor

  ; loop over all FoF halos in groupcat1
  for i=0,gc1.nGroupsTot-1 do begin
    ; find primary subgroup child and so position, skip if none
    gcInd = gcPriChildInd(gc=gc1,haloID=i)
    if (gcInd eq -1) then continue
    
    ; compute distances from subfind CM center and cut at maximum
    hPos  = idMBCM1[*,gcInd]
    hMass = gc1.groupMass[i]
    
    dists = periodicDists(hPos,targetPos,sP=sP1)
    minDist = min(dists)
    
    if (minDist le distTol) then begin
      ; select closest halo
      w = where(dists eq minDist,count)
      if (count eq 1) then begin 
        ; check mass agreement
        targInd = w[0]
        massDiff = abs(hMass-gc2.groupMass[targInd]) / hMass
        
        massDiffs[i]   = massDiff
        posDiffs[i]    = minDist        
        
        if (massDiff lt massTol) then begin
          ; select halo
          ;print,'['+str(i)+'] Matched ID ['+str(targInd)+'] at distance = '+string(minDist,format='(f5.2)')+$
          ;      ' Mass Diff = '+string(massDiff,format='(f5.2)')+' '+str(massTol)
          matchedInds[i] = targInd
        endif else begin;mass
          ;print,massDiff,massTol
        endelse
      endif ;count
    endif ;dist
  endfor
  
  wMatch = where(matchedInds ne -1,nMatched)

  r = {matchedInds:matchedInds,massDiffs:massDiffs,posDiffs:posDiffs,wMatch:wMatch,nMatched:nMatched}
  
  ;save
  ; save
  save,r,filename=saveFilename
  
  return,r
end

; cosmoCompareHaloCenters(): compare relative differences between different calcuations of the centers
;                            of halos, and verify group catalog consistency (group-subgroup mappings)

pro cosmoCompareHaloCenters

  compile_opt idl2, hidden, strictarr, strictarrsubs

  res = 128
  run = 'dev.tracer.nograd'
  redshift = 3.0
  
  sP = simParams(res=res,run=run,redshift=redshift)
  
  gc    = loadGroupCat(sP=sP,/skipIDs)
  gcCen = subgroupPosByMostBoundID(sP=sP)
  
  binSize = 1.0
  min = 0.0
  max = 50.0
  
  num = 3500
  
  ; group catalog consistency checks
  print,'nGroups nSubgroups',gc.nGroupsTot,gc.nSubgroupsTot
  ; should pass this
  for i=0,gc.nGroupsTot-1 do begin
    gcInd = gcPriChildInd(gc=gc,haloID=i)
    if (gcInd eq -1) then continue
    grNr = gc.subgroupGrNr[gcInd]
    if (grNr ne i) then stop
  endfor
  
  ; centers
  cen_fof  = gc.groupPos[*,0:num]
  
  cen_sfcm = fltarr(3,num)
  cen_mb   = fltarr(3,num)
  
  w = []
  
  for i=0,num-1 do begin
    gcInd = gcPriChildInd(gc=gc,haloID=i)
    if (gcInd ne -1) then begin
      w = [w,i]
      cen_sfcm[*,i] = gc.subgroupCM[*,gcInd]
      cen_mb[*,i] = gcCen[*,gcInd]
    endif
  endfor
  
  ; distances
  dist_fof_sfcm = periodicDists(cen_fof,cen_sfcm,sP=sP)
  dist_fof_mb   = periodicDists(cen_fof,cen_mb,sP=sP)
  dist_sfcm_mb  = periodicDists(cen_sfcm,cen_mb,sP=sP)

  ; plot histos
  hist_fof_sfcm = histogram(dist_fof_sfcm[w],binsize=binSize,min=min,max=max,loc=loc1)
  hist_fof_mb   = histogram(dist_fof_mb[w],binsize=binSize,min=min,max=max,loc=loc2)
  hist_sfcm_mb  = histogram(dist_sfcm_mb[w],binsize=binSize,min=min,max=max,loc=loc3)
  
  start_PS,'hist_fof_sfcm.eps'
    fsc_plot,loc1,hist_fof_sfcm,xtitle="dist [kpc]",ytitle="N"
  end_PS
  
  start_PS,'hist_fof_mb.eps'
    fsc_plot,loc2,hist_fof_mb,xtitle="dist [kpc]",ytitle="N"
  end_PS
  
  start_PS,'hist_sfcm_mb.eps'
    fsc_plot,loc3,hist_sfcm_mb,xtitle="dist [kpc]",ytitle="N"
  end_PS
  
  stop

end

; snapNumToRedshift(): convert snapshot number to redshift or time (scale factor)

function snapNumToRedshift, time=time, all=all, sP=sP, snap=snap

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if not keyword_set(sP) then message,'Error: Need sP to convert snapshot number to redshift!'

  saveFileName = sP.derivPath + sP.savPrefix + '_snapnum.redshift.sav'

  if not file_test(saveFileName) then stop
  if not keyword_set(snap) then snap = sP.snap
  if snap eq -1 then stop
  
  ; restore
  restore,saveFilename

  if (not keyword_set(time)) then begin
    if (keyword_set(all)) then return,redshifts
    
    if (snap ge 0 and snap lt n_elements(redshifts)) then $
      return,redshifts[snap]
  endif else begin
    if (keyword_set(all)) then return,times
      
    if (snap ge 0 and snap lt n_elements(redshifts)) then $
      return,times[snap]
  endelse

end

; snapNumToAge(): convert snapshot number to approximate age of the universe at that point

function snapNumToAge, snap=snap, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if not keyword_set(snap) then snap = sP.snap
  z = snapNumToRedshift(snap=snap,sP=sP)
  return, redshiftToAge(z)
end
  
; rhoTHisto(): make mass-weighted density-temperature 2d histogram

function rhoTHisto, dens_in, temp_in, mass=mass, nbins=nbins, plot=plot

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; config
  if not keyword_set(nbins) then nbins = 20
  dens = alog10(rhoRatioToCrit(dens_in))
  temp = alog10(temp_in)
  rMinMax = [-2.0,8.0]
  tMinMax = [3.0,7.0]
  
  ; calculate bin sizes
  binSizeRho  = (rMinMax[1]-rMinMax[0]) / nbins
  binSizeTemp = (tMinMax[1]-tMinMax[0]) / nbins
  
  ; if not mass weighting assign equal weights
  if not keyword_set(mass) then mass = replicate(1.0,n_elements(dens))
  
  ; use hist_nd with centered bins and boundary fixes
  h2rt = hist_nd_weight(transpose([[dens],[temp]]),weight=mass,[binSizeRho,binSizeTemp],$
                        min=[rMinMax[0]-binSizeRho*0.5,tMinMax[0]-binSizeTemp*0.5],$
                        max=[rMinMax[1]+binSizeRho*0.49,tMinMax[1]+binSizeTemp*0.49])
  
  ; plot
  if keyword_set(plot) then begin
    ; color table
    loadct, 2, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
      
    tvim,h2rt,pcharsize=!p.charsize-1.0,scale=1,clip=[10,100],$;,/c_map
         xtitle="log ("+textoidl("\rho / \rho_{crit}")+")",ytitle="log (T [K])",$
         stitle="Total Mass ("+textoidl("M_{sun}")+")",barwidth=0.5,lcharsize=!p.charsize-1.5,$
         xrange=[-2.0,8.0],yrange=[3.0,7.0],$;xrange=rMinMax,yrange=tMinMax,$
         /rct;,nodata=0,rgb_nodata=[1.0,1.0,1.0] ;display zeros as white not black
  endif
  
  return,h2rt

end

; redshift_axis(): draw redshift axis

pro redshift_axis, xRange, yRange, ylog=ylog, sP=sP, dotted=dotted, zTicknames=zTicknames

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  logFac = 1.0
  if keyword_set(ylog) then logFac = 10.0
  
  yTickLen = (yRange[1]-yRange[0]) / 40.0 * logFac
  yTextOff = (yRange[1]-yRange[0]) / 50.0 * logFac

  if (not keyword_set(zTicknames)) then $
    zTicknames = ['30','6','4','3','2','1','0.5','0.25','0']
  nZ = n_elements(zTicknames)
  
  zXPos = fltarr(nZ)
  
  ; plot "maximum" redshift label
  if (xRange[0] le 0.0) then $
    fsc_text,0.0,yRange[1]+yTextOff,zTicknames[0],alignment=0.5
  
  ; skip z=30 (highest) at t=0
  for i=1,nZ-1 do begin
    if (sP.snap ne -1) then begin ;x-axis in snapshot number
      zXPos[i] = redshiftToSnapNum(float(zTicknames[i]),sP=sP)
    endif else begin ;x-axis in time [gyr]
      zXPos[i] = redshiftToAge(float(zTicknames[i]))
    endelse
    
    ; plot tick mark and label if inside plotrange
    if (zXPos[i] ge xRange[0] and zXPos[i] le xRange[1]) then begin
      fsc_plot,[zXPos[i],zXPos[i]],[yRange[1],yRange[1]-yTickLen],/overplot
      fsc_text,zXPos[i],yRange[1]+yTextOff,zTicknames[i],alignment=0.5
    endif
    
    ; plot vertical dotted line at each redshift mark if requested
    if keyword_set(dotted) then $
      fsc_plot,[zXPos[i],zXPos[i]],yRange,line=dotted,/overplot,thick=!p.thick-0.5
  endfor
  
  fsc_plot,xRange,[yRange[1],yRange[1]],/overplot
  fsc_text,0.5,0.94,"Redshift",/normal

end

; universeage_axis(): draw age of universe (elapsed) axis

pro universeage_axis, xRange, yRange, ylog=ylog, dotted=dotted

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  logFac = 1.0
  if keyword_set(ylog) then logFac = 10.0
  
  yTickLen = (yRange[1]-yRange[0]) / 40.0 * logFac
  yTextOff = (yRange[1]-yRange[0]) / 50.0 * logFac

  if (not keyword_set(zTicknames)) then begin
    zTicknames = ['0.1','1','1.5','2','3','4','5','6','8','13.8'] ;10=0.33
    zRedshifts = [30,5.7,4.21,3.25,2.23,1.65,1.265,0.98,0.595,0.0]
  endif
  nZ = n_elements(zTicknames)
  
  zXPos = fltarr(nZ)
  
  ; plot "maximum" redshift label
  if (xRange[0] le 0.0) then $
    fsc_text,0.0,yRange[1]+yTextOff,zTicknames[0],alignment=0.5
  
  ; skip t=0 (highest) at z=inf (?)
  for i=1,nZ-1 do begin

    zXPos[i] = zRedshifts[i]
    
    ; plot tick mark and label if inside plotrange
    if (zXPos[i] le xRange[0] and zXPos[i] ge xRange[1]) then begin
      fsc_plot,[zXPos[i],zXPos[i]],[yRange[1],yRange[1]-yTickLen],/overplot
      fsc_text,zXPos[i],yRange[1]+yTextOff,zTicknames[i],alignment=0.5
    endif
    
    ; plot vertical dotted line at each redshift mark if requested
    if keyword_set(dotted) then $
      fsc_plot,[zXPos[i],zXPos[i]],yRange,line=dotted,/overplot,thick=!p.thick-0.5
  endfor
  
  fsc_plot,xRange,[yRange[1],yRange[1]],/overplot
  fsc_text,(xRange[0]-xRange[1])/2.0,yRange[1]+yTextOff*5.0,"Time [Gyr]",alignment=0.5

end
