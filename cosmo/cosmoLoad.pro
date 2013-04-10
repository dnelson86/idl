; cosmoLoad.pro
; cosmological simulations - loading procedures (snapshots, fof/subhalo group cataloges)
; dnelson mar.2013

; getTypeSortedIDList(): within the group catalog ID list rearrange the IDs for each FOF group to be 
;                        ordered first by type (not by SubNr since Subfind was run) such that the
;                        groupOffsetType indexing can be used

function getTypeSortedIDList, sP=sP, gc=gc

  ; load header
  h  = loadsnapshotHeader(sP=sP)
  
  if (h.nPartTot[5] ne 0) then stop ; not implemented
  
  ; save/restore
  saveFilename = sP.derivPath + 'gcSortedIDs.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.sav'
                 
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,sortedIDList
  endif
    
  ; for each type, load IDs, match to group cat IDs
  count_gas = 0
  if h.nPartTot[partTypeNum('gas')] gt 0 then begin
    gas_ids  = loadsnapshotSubset(sP=sP,partType='gas',field='ids')
    calcMatch,gas_ids,gc.IDs,gas_ind,gc_ind_gas,count=count_gas
    gas_ids  = !NULL
    
    ; reorder indices into ID arrays to match order found in gc.IDs
    gc_ind_gas   = gc_ind_gas[sort(gc_ind_gas)]
    
    ; verify counts
    if (count_gas ne total(gc.groupLenType[partTypeNum('gas'),*])) then stop
  endif
  
  count_dm = 0
  if h.nPartTot[partTypeNum('dm')] gt 0 then begin
    dm_ids = loadsnapshotSubset(sP=sP,partType='dm',field='ids')
    calcMatch,dm_ids,gc.IDs,dm_ind,gc_ind_dm,count=count_dm
    dm_ids = !NULL
    gc_ind_dm = gc_ind_dm[sort(gc_ind_dm)]
    if (count_dm ne total(gc.groupLenType[partTypeNum('dm'),*])) then stop
  endif
  
  count_trvel = 0
  if h.nPartTot[partTypeNum('tracerVel')] gt 0 then begin
    trvel_ids = loadsnapshotSubset(sP=sP,partType='tracerVel',field='ids')
    calcMatch,trvel_ids,gc.IDs,trvel_ind,gc_ind_trvel,count=count_trvel
    trvel_ids = !NULL
    gc_ind_trvel = gc_ind_trvel[sort(gc_ind_trvel)]
    if (count_trvel ne total(gc.groupLenType[partTypeNum('tracerVel'),*])) then stop
  endif
  
  count_star = 0
  if h.nPartTot[partTypeNum('stars')] gt 0 then begin
    star_ids = loadsnapshotSubset(sP=sP,partType='star',field='ids')
    calcMatch,star_ids,gc.IDs,star_ind,gc_ind_star,count=count_star
    star_ids = !NULL
    gc_ind_star  = gc_ind_star[sort(gc_ind_star)]
    if (count_star ne total(gc.groupLenType[partTypeNum('stars'),*])) then stop
  endif
  
  if (count_gas + count_dm + count_trvel + count_star ne total(gc.groupLen)) then stop
  
  start_gas   = 0ULL
  start_dm    = 0ULL
  start_trvel = 0ULL
  start_star  = 0ULL
  offset      = 0ULL
  
  ; DEBUG:
  ;mask_gas   = intarr(max(gc.IDs[gc_ind_gas])+1)
  ;mask_dm    = intarr(max(gc.IDs[gc_ind_dm])+1)
  ;mask_trvel = intarr(max(gc.IDs[gc_ind_trvel])+1)
  ;mask_star  = intarr(max(gc.IDs[gc_ind_star])+1)
  
  sortedIDList = ulon64arr(gc.nIDsTot)
  
  ; loop over each fof group
  for i=0L,gc.nGroupsTot-1 do begin
    ; gas
    if (gc.groupLenType[partTypeNum('gas'),i] gt 0) then begin
      halo_gas_ids  = gc.IDs[gc_ind_gas[start_gas:start_gas+gc.groupLenType[partTypeNum('gas'),i]-1]]
      if n_elements(halo_gas_ids)  ne gc.groupLenType[partTypeNum('gas'),i] then stop
      
      ; DEBUG: fill mask as check for duplicates
      ;mask_gas[halo_gas_ids] += 1
      
      ; fill in sorted ID list
      sortedIDList[offset:offset+n_elements(halo_gas_ids)-1] = halo_gas_ids
      offset += n_elements(halo_gas_ids)    
      
      ; increment
      start_gas  += n_elements(halo_gas_ids)
    endif
    
    ; dm
    if (gc.groupLenType[partTypeNum('dm'),i] gt 0) then begin
      halo_dm_ids   = gc.IDs[gc_ind_dm[start_dm:start_dm+gc.groupLenType[partTypeNum('dm'),i]-1]]
      if n_elements(halo_dm_ids)   ne gc.groupLenType[partTypeNum('dm'),i] then stop
      ;mask_dm[halo_dm_ids] += 1
      sortedIDList[offset:offset+n_elements(halo_dm_ids)-1] = halo_dm_ids
      offset += n_elements(halo_dm_ids)
      start_dm   += n_elements(halo_dm_ids)
    endif
    
    ; velocity tracers
    if (gc.groupLenType[partTypeNum('tracerVel'),i] gt 0) then begin
      halo_trvel_ids = gc.IDs[gc_ind_trvel[start_trvel:start_trvel+gc.groupLenType[partTypeNum('tracerVel'),i]-1]]
      if n_elements(halo_trvel_ids) ne gc.groupLenType[partTypeNum('tracerVel'),i] then stop
      ;mask_trvel[halo_trvel_ids] += 1
      sortedIDList[offset:offset+n_elements(halo_trvel_ids)-1] = halo_trvel_ids
      offset += n_elements(halo_trvel_ids)
      start_trvel += n_elements(halo_trvel_ids)
    endif
    
    ; stars
    if (gc.groupLenType[partTypeNum('stars'),i] gt 0) then begin
      halo_star_ids = gc.IDs[gc_ind_star[start_star:start_star+gc.groupLenType[partTypeNum('stars'),i]-1]]
      if n_elements(halo_star_ids) ne gc.groupLenType[partTypeNum('stars'),i] then stop
      ;mask_star[halo_star_ids] += 1
      sortedIDList[offset:offset+n_elements(halo_star_ids)-1] = halo_star_ids
      offset += n_elements(halo_star_ids)
      start_star += n_elements(halo_star_ids)
    endif
  endfor
  
  ; DEBUG: check masks
  ;w_gas   = where(mask_gas gt 1,count_gas)
  ;w_dm    = where(mask_dm gt 1,count_dm)
  ;w_trvel = where(mask_trvel gt 1,count_trvel)
  ;w_star  = where(mask_star gt 1,count_star)
  
  ;if (count_gas gt 0 or count_dm gt 0 or count_trvel gt 0 or count_star gt 0) then stop
  
  ; DEBUG: match old and sorted ID lists for consistency
  calcMatch,gc.IDs,sortedIDList,ind1,ind2,count=count
  if (count ne n_elements(sortedIDList)) then stop
  
  if min(sortedIDList) lt 0 then stop ; check for 32 bit long overflow
  
  ; save
  save,sortedIDList,filename=saveFilename

  return, sortedIDList

end

; getGroupCatFilename(): take input path and snapshot number and find fof/subfind group catalog filename

function getGroupCatFilename, fileBase, snapNum=m

  flag = 0
  
  ; m='none' if directly specifying file name
  if (str(m) eq 'none') then begin
    ext = ''
    f = fileBase
    
    if file_test(f) then return,f
  endif
  
  ; check for '/' on end of fileBase
  lastChar = strmid(fileBase,strlen(fileBase)-1,1)
  if (lastChar ne '/') then fileBase += '/'
    
  ; format snap number
  if (m le 999) then $
    ext = string(m,format='(I3.3)')
  if (m gt 999) then $
    ext = string(m,format='(I4.4)')
      
  ; check for single (non-split) in root directory
  if file_test(fileBase+'/fof_subhalo_tab_'+ext+'.hdf5') then $
    return, [fileBase + 'fof_subhalo_tab_' + ext + '.hdf5']
    
  ; check for both fof+subfind and just fof
  fBases = [ fileBase + 'groups_' + ext + '/fof_subhalo_tab_' + ext ,$
             fileBase + 'groups_' + ext + '/fof_tab_' + ext ]
             
  foreach f,fBases do begin

    ; check for single inside a directory?
    if file_test(f+'.hdf5') then $
      return, [f+'.hdf5']
    
    ; check for split
    if file_test(f+'.0.hdf5') then begin
      ff = file_search(f+'.*hdf5')
      return, f + '.' + str(indgen(n_elements(ff))) + '.hdf5'
    endif
  
  endforeach
    
  print,'Error: Failed to find group catalog file(s).'

end

; loadGroupCat(): load new HDF5 fof/subfind group catalog for a given snapshot
;                 construct offset tables and return structure
;                     
; readIDs=1 : by default, skip IDs since we operate under the group ordered snapshot assumption, but
;             if this flag is set then read IDs and include them (if they exist)
;             also generate (GrNr,Type) sorted id list
; skipIDs=1 : acknowledge we are working with a STOREIDS type .hdf5 group cat and don't warn
; getSortedIDs=1 : create a second ID list sorted by GrNr->Type for use with GroupOffsetType

function loadGroupCat, sP=sP, readIDs=readIDs, skipIDs=skipIDsFlag, getSortedIDs=getSortedIDs, verbose=verbose

  forward_function loadSubhaloGroups
  if not keyword_set(verbose) then verbose = 0
  !except = 0 ;suppress floating point underflow/overflow errors

  fileList = getGroupCatFileName(sP.simPath,snapNum=sP.snap)

  nFiles = n_elements(fileList)

  ; if new group catalog not found, look for old format
  if nFiles eq 0 then return,loadSubhaloGroups(sP=sP,skipIDs=skipIDs,verbose=verbose)

  ; load number of split files from header of first part
  hdf5s    = h5_parse(fileList[0]) ;structure only
  NumFiles = hdf5s.Header.NumFiles._DATA
  
  ; what additional fields do we have?
  SubfindExistsFlag    = tag_exist(hdf5s,'SubhaloLen')
  NewSubfindExistsFlag = tag_exist(hdf5s,'SubhaloSFRinRad')
  GFMExistsFlag        = tag_exist(hdf5s,'GroupGasMetalFractions')
  if GFMExistsFlag then gfmNumElements = hdf5s.Group.GroupGasMetalFractions._DIMENSIONS[0]
  if GFMExistsFlag then gfmNumPhotometrics = hdf5s.Subhalo.SubhaloStellarPhotometrics._DIMENSIONS[0]
  
  if (NumFiles ne nFiles) then $
    message,'ERROR: NumFiles ['+str(NumFiles)+'] differs from number of files found ['+str(nFiles)+'].'
  
  if (verbose) then $
    print,'Loading group catalog from snapshot ('+str(sP.snap)+') in [' + str(NumFiles) + '] files.'  
  
  ; counters
  nGroupsTot    = 0L
  nIDsTot       = 0L
  nSubgroupsTot = 0L
  
  skip    = 0L
  skipSub = 0L
  skipIDs = 0L
  
  ; load across all file parts
  for i=0,nFiles-1 do begin
    ; load header
    fileID = h5f_open(fileList[i])
    
    s = h5_parse(fileID,"Header")
    
    h = {                                                          $
          nGroups             : s.nGroups_ThisFile._DATA          ,$
          nGroupsTot          : s.nGroups_Total._DATA             ,$
          nIDs                : s.nIDs_ThisFile._DATA             ,$
          nIDsTot             : s.nIDs_Total._DATA                ,$
          nSubgroups          : s.nSubgroups_ThisFile._DATA       ,$
          nSubgroupsTot       : s.nSubgroups_Total._DATA          ,$
          numFiles            : s.NumFiles._DATA                  ,$
          flagDoublePrecision : s.flagDoublePrecision._DATA       $
        }
         
    ; add counters
    nGroupsTot    += h.nGroups
    nIDsTot       += h.nIDs
    nSubgroupsTot += h.nSubgroups
          
    ; allocate storage if this is the first iteration
    if (i eq 0) then begin
      sf  = {                                      $
        GroupLen        : ulonarr(h.nGroupsTot)   ,$
        GroupLenType    : ulonarr(6,h.nGroupsTot) ,$
        GroupMass       : fltarr(h.nGroupsTot)    ,$
        GroupMassType   : fltarr(6,h.nGroupsTot)  ,$   
        GroupPos        : fltarr(3,h.nGroupsTot)  ,$
        GroupVel        : fltarr(3,h.nGroupsTot)  ,$
        GroupSFR        : fltarr(h.nGroupsTot)    ,$
                                                   $
        GroupOffset     : ulonarr(h.nGroupsTot)   ,$
        GroupOffsetType : ulonarr(6,h.nGroupsTot) ,$
                                                   $
        nGroupsTot          : h.nGroupsTot                ,$
        nSubgroupsTot       : h.nSubgroupsTot             ,$
        nIDsTot             : h.nIDsTot                    }
      
      ; SUBFIND products
      if (SubfindExistsFlag eq 1) then begin
        sfsub = {                                  $
        Group_M_Mean200 : fltarr(h.nGroupsTot)    ,$
        Group_R_Mean200 : fltarr(h.nGroupsTot)    ,$
        Group_M_Crit200 : fltarr(h.nGroupsTot)    ,$
        Group_R_Crit200 : fltarr(h.nGroupsTot)    ,$
        Group_M_TH200   : fltarr(h.nGroupsTot)    ,$
        Group_R_TH200   : fltarr(h.nGroupsTot)    ,$
        GroupNsubs      : ulonarr(h.nGroupsTot)   ,$
        GroupFirstSub   : ulonarr(h.nGroupsTot)   ,$
                                                   $
        SubgroupLen         : ulonarr(h.nSubgroupsTot)    ,$
        SubgroupLenType     : ulonarr(6,h.nSubgroupsTot)  ,$
        SubgroupOffset      : ulonarr(h.nSubgroupsTot)    ,$
        SubgroupOffsetType  : ulonarr(6,h.nSubgroupsTot)  ,$
        
        SubgroupMass        : fltarr(h.nSubgroupsTot)     ,$
        SubgroupMassType    : fltarr(6,h.nSubgroupsTot)   ,$  
        SubgroupPos         : fltarr(3,h.nSubgroupsTot)   ,$
        SubgroupVel         : fltarr(3,h.nSubgroupsTot)   ,$
        SubgroupCM          : fltarr(3,h.nSubgroupsTot)   ,$
        SubgroupSpin        : fltarr(3,h.nSubgroupsTot)   ,$
        SubgroupVelDisp     : fltarr(h.nSubgroupsTot)     ,$
        SubgroupVmax        : fltarr(h.nSubgroupsTot)     ,$
        SubgroupVmaxRad     : fltarr(h.nSubgroupsTot)     ,$
        SubgroupHalfMassRad : fltarr(h.nSubgroupsTot)     ,$
        SubgroupIDMostBound : lon64arr(h.nSubgroupsTot)   ,$
        SubgroupGrNr        : lonarr(h.nSubgroupsTot)     ,$
        SubgroupParent      : ulonarr(h.nSubgroupsTot)     }
        
        sf = create_struct(sf,sfsub) ;concat
      endif
      
      ; SUBFIND newer products
      if NewSubfindExistsFlag then begin
        sfsub = {                                                     $
        SubgroupHalfMassRadType       : fltarr(6,h.nSubgroupsTot)    ,$
        SubgroupMassInRad             : fltarr(h.nSubgroupsTot)      ,$
        SubgroupMassInRadType         : fltarr(6,h.nSubgroupsTot)    ,$
        SubgroupSFR                   : fltarr(h.nSubgroupsTot)      ,$
        SubgroupSFRInRad              : fltarr(h.nSubgroupsTot)       }
        
        sf = create_struct(sf,sfsub)
      endif
      
      ; GFM products
      if GFMExistsFlag then begin
        sfsub = {                                                         $
        GroupGasMetalFractions  : fltarr(gfmNumElements,h.nGroupsTot)    ,$
        GroupGasMetallicity     : fltarr(h.nGroupsTot)                   ,$
        GroupStarMetalFractions : fltarr(gfmNumElements,h.nGroupsTot)    ,$
        GroupStarMetallicity    : fltarr(h.nGroupsTot)                   ,$
                                                                          $
        SubgroupGasMetalFractions     : fltarr(gfmNumElements,h.nSubgroupsTot)    ,$
        SubgroupGasMetalFractionsSFR  : fltarr(gfmNumElements,h.nSubgroupsTot)    ,$
        SubgroupGasMetallicity        : fltarr(h.nSubgroupsTot)                   ,$
        SubgroupGasMetallicitySFR     : fltarr(h.nSubgroupsTot)                   ,$
        SubgroupStarMetalFractions    : fltarr(gfmNumElements,h.nSubgroupsTot)       ,$
        SubgroupStarMetallicity       : fltarr(h.nSubgroupsTot)                      ,$
        SubgroupStellarPhotometrics   : fltarr(gfmNumPhotometrics,h.nSubgroupsTot)    }
        
        sf = create_struct(sf,sfsub)
      endif
      
      ; ID load requested?
      if keyword_set(readIDs) then begin
        if (h.nIDsTot eq 0) then begin
          print,'Warning: readIDs requested but no IDs in group catalog!' & stop
        endif
        
        sfids = { IDs:lon64arr(h.nIDsTot) }
        sf = create_struct(sf,sfids) ;concat
      endif
      
    endif ;i=0
    
    ; fill sf with group data from this part
    sf.GroupLen        [skip:(skip+h.nGroups-1)]   = h5d_read(h5d_open(fileID,"Group/GroupLen"))
    sf.GroupLenType    [*,skip:(skip+h.nGroups-1)] = h5d_read(h5d_open(fileID,"Group/GroupLenType"))
    sf.GroupMass       [skip:(skip+h.nGroups-1)]   = h5d_read(h5d_open(fileID,"Group/GroupMass"))
    sf.GroupMassType   [*,skip:(skip+h.nGroups-1)] = h5d_read(h5d_open(fileID,"Group/GroupMassType"))
    sf.GroupPos        [*,skip:(skip+h.nGroups-1)] = h5d_read(h5d_open(fileID,"Group/GroupPos"))
    sf.GroupVel        [*,skip:(skip+h.nGroups-1)] = h5d_read(h5d_open(fileID,"Group/GroupVel"))
    if tag_exist(hdf5s,'GroupSFR') then $
      sf.GroupSFR        [skip:(skip+h.nGroups-1)]   = h5d_read(h5d_open(fileID,"Group/GroupSFR"))
    
    if SubfindExistsFlag then begin
      ; these group properties only exist if subfind was run
      sf.Group_M_Mean200 [skip:(skip+h.nGroups-1)]   = h5d_read(h5d_open(fileID,"Group/Group_M_Mean200"))
      sf.Group_R_Mean200 [skip:(skip+h.nGroups-1)]   = h5d_read(h5d_open(fileID,"Group/Group_R_Mean200"))
      sf.Group_M_Crit200 [skip:(skip+h.nGroups-1)]   = h5d_read(h5d_open(fileID,"Group/Group_M_Crit200"))
      sf.Group_R_Crit200 [skip:(skip+h.nGroups-1)]   = h5d_read(h5d_open(fileID,"Group/Group_R_Crit200"))
      sf.Group_M_TH200   [skip:(skip+h.nGroups-1)]   = h5d_read(h5d_open(fileID,"Group/Group_M_TopHat200"))
      sf.Group_R_TH200   [skip:(skip+h.nGroups-1)]   = h5d_read(h5d_open(fileID,"Group/Group_R_TopHat200"))
      
      sf.GroupNsubs      [skip:(skip+h.nGroups-1)]   = h5d_read(h5d_open(fileID,"Group/GroupNsubs"))
      sf.GroupFirstsub   [skip:(skip+h.nGroups-1)]   = h5d_read(h5d_open(fileID,"Group/GroupFirstSub"))
    endif
    
    if GFMExistsFlag then begin
      sf.GroupGasMetalFractions  [*,skip:(skip+h.nGroups-1)] = h5d_read(h5d_open(fileID,"Group/GroupGasMetalFractions"))
      sf.GroupGasMetallicity     [skip:(skip+h.nGroups-1)] = h5d_read(h5d_open(fileID,"Group/GroupGasMetallicity"))
      sf.GroupStarMetalFractions [*,skip:(skip+h.nGroups-1)] = h5d_read(h5d_open(fileID,"Group/GroupStarMetalFractions"))
      sf.GroupStarMetallicity    [skip:(skip+h.nGroups-1)] = h5d_read(h5d_open(fileID,"Group/GroupStarMetallicity"))    
    endif
    
    skip += h.nGroups
    
    ; fill sf with subhalo data from this part
    if SubfindExistsFlag then begin
      sf.SubgroupLen     [skipSub:(skipSub+h.nSubgroups-1)]   = h5d_read(h5d_open(fileID,"Subhalo/SubhaloLen"))
      sf.SubgroupLenType [*,skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloLenType"))
      sf.SubgroupMass    [skipSub:(skipSub+h.nSubgroups-1)]   = h5d_read(h5d_open(fileID,"Subhalo/SubhaloMass"))
      sf.SubgroupMassType[*,skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloMassType"))
      sf.SubgroupPos     [*,skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloPos"))
      sf.SubgroupVel     [*,skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloVel"))
      sf.SubgroupCM      [*,skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloCM"))
      sf.SubgroupSpin    [*,skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloSpin"))
      
      sf.SubgroupVelDisp [skipSub:(skipSub+h.nSubgroups-1)]   = h5d_read(h5d_open(fileID,"Subhalo/SubhaloVelDisp"))
      sf.SubgroupVmax    [skipSub:(skipSub+h.nSubgroups-1)]   = h5d_read(h5d_open(fileID,"Subhalo/SubhaloVmax"))
      sf.SubgroupVmaxRad [skipSub:(skipSub+h.nSubgroups-1)]   = h5d_read(h5d_open(fileID,"Subhalo/SubhaloVmaxRad"))
      sf.SubgroupHalfMassRad[skipSub:(skipSub+h.nSubgroups-1)]   = h5d_read(h5d_open(fileID,"Subhalo/SubhaloHalfmassRad"))
      sf.SubgroupIDMostBound[skipSub:(skipSub+h.nSubgroups-1)]   = h5d_read(h5d_open(fileID,"Subhalo/SubhaloIDMostbound"))
      
      sf.SubgroupGrnr    [skipSub:(skipSub+h.nSubgroups-1)]   = h5d_read(h5d_open(fileID,"Subhalo/SubhaloGrNr"))
      sf.SubgroupParent  [skipSub:(skipSub+h.nSubgroups-1)]   = h5d_read(h5d_open(fileID,"Subhalo/SubhaloParent"))
      
      if NewSubfindExistsFlag then begin
        sf.SubgroupHalfMassRadType[*,skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloHalfmassRadType"))
        sf.SubgroupMassInRad      [skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloHalfmassRad"))
        sf.SubgroupMassInRadType  [*,skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloMassInRadType"))
        sf.SubgroupSFR            [skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloSFR"))
        sf.SubgroupSFRInRad       [skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloSFRinRad"))
      endif
    
      if GFMExistsFlag then begin
        sf.SubgroupGasMetalFractions   [*,skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloGasMetalFractions"))
        sf.SubgroupGasMetalFractionsSFR[*,skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloGasMetalFractionsSfr"))
        sf.SubgroupGasMetallicity      [skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloGasMetallicity"))
        sf.SubgroupGasMetallicitySFR   [skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloGasMetallicitySfr"))
        sf.SubgroupStarMetalFractions  [*,skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloStarMetalFractions"))
        sf.SubgroupStarMetallicity     [skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloStarMetallicity"))
        sf.SubgroupStellarPhotometrics [*,skipSub:(skipSub+h.nSubgroups-1)] = h5d_read(h5d_open(fileID,"Subhalo/SubhaloStellarPhotometrics"))
      endif
      
      skipSub += h.nSubgroups
    endif
    
    ; fill sf with IDs from this part (if requested)
    if keyword_set(readIDs) then begin
      sf.IDs[skipIDs:(skipIDs+h.nIDs-1)] = h5d_read(h5d_open(fileID,"IDs/ID"))
      skipIDs += h.nIDs
    endif
    
    ; close file
    h5f_close, fileID
  
  endfor
  
  ; create group offset table
  ; when subfind is run, sort to create ID list is: (1) GrNr, (2) SubNr, (3) Type, (4) BindingEnergy
  for GrNr=1L, h.nGroupsTot-1 do begin
    sf.GroupOffset[GrNr] = sf.GroupOffset[GrNr-1] + sf.GroupLen[GrNr-1]
  endfor
  
  ; NOTE: given SubNr sorted before Type, GroupOffsetType can ONLY be used on IDsSorted (NOT IDs!)
  for GrNr=0L, h.nGroupsTot-1 do begin
    ; create group offset type table
    typeCumSum = [0L,total(sf.groupLenType[0:4,GrNr],/cum,/int)]
    sf.GroupOffsetType[*,GrNr] = sf.GroupOffset[GrNr] + typeCumSum
    
    ; create subgroup offset tables (for groups with at least one subgroup)
    ; NOTE: "fuzz" group IDs for each group are sorted last within each group
    if (SubfindExistsFlag eq 1) then begin
      if (sf.groupNSubs[GrNr] gt 0) then begin
        ; first subgroup
        SubNr = sf.groupFirstSub[GrNr]
        sf.subgroupOffset[SubNr] = sf.groupOffset[GrNr]
        
        ; subsequent subgroups
        if (sf.groupNSubs[GrNr] gt 1) then begin
          SubOffsets = total(sf.subgroupLen[SubNr:SubNr+sf.groupNSubs[GrNr]-2],/cum,/int)
          sf.subgroupOffset[SubNr+1:SubNr+sf.groupNSubs[GrNr]-1] = SubOffsets + sf.groupOffset[GrNr]
        endif
        
        ; construct subgroup type offset table in a loop
        for SubNr=sf.groupFirstSub[GrNr],sf.groupFirstSub[GrNr]+sf.groupNSubs[GrNr]-1 do begin
          subTypeCumSum = [0L,total(sf.subgroupLenType[0:4,SubNr],/cum,/int)]
          sf.subgroupOffsetType[*,SubNr] = sf.subgroupOffset[SubNr] + subTypeCumSum
        endfor
        
      endif
    endif ;SubfindExistsFlag
    
  endfor
  
  ; if ID read requested, create typeSortedIDList (and save), add to return structure
  if keyword_set(readIDs) and keyword_set(getSortedIDs) then begin
    sfsorted = { IDsSorted:getTypeSortedIDList(sP=sP,gc=sf) }
    sf = create_struct(sf,sfsorted) ;concat
  endif
  
  ; verify accumulated totals with last header totals
  if ((nGroupsTot ne h.nGroupsTot) or (nSubgroupsTot ne h.nSubgroupsTot) or $ 
      (nIDsTot ne h.nIDsTot and keyword_set(readIDs))) then begin
    print,'ERROR: Totals do not add up.'
    stop
  endif
  
  ; check for 32 bit long overflow
  if keyword_set(readIDs) then if min(sf.IDs) lt 0 then message,'Error: ID overflow.'
  if (SubfindExistsFlag eq 1) then if min(sf.subgroupIDMostBound) lt 0 then message,'Error: ID overflow2.'

  ; if ID read was not requested but IDs exist, stop for now (possibly under the group ordered assumption)
  if (nIDsTot gt 0 and n_elements(readIDs) eq 0 and n_elements(skipIDsFlag) eq 0) then begin
    print,'Warning: readIDs not requested, but IDs present in group catalog!'
  endif
  
  !except = 1
  
  return,sf
end

; loadSubhaloGroups(): load (OLD, not HDF5) complete subfind group catalog for a given snapshot
; skipIDs=1 : don't load actual group member particle IDs

function loadSubhaloGroups, sP=sP, verbose=verbose, skipIDs=skipIDs

  print,'WARNING: OLD GROUP CATALOG'
  if not keyword_set(verbose) then verbose = 0

  ; set filename
  ext = string(sP.snap,format='(I3.3)')
  fIDs = sP.simPath + 'groups_' + ext + '/subhalo_ids_' + ext
  fTab = sP.simPath + 'groups_' + ext + '/subhalo_tab_' + ext
  
  ; check existance and multiple outputs
  if not file_test(fIDs) then begin
    if (file_test(fIDs+'.0')) then begin
      ; split into multiples, get count
      nSplit_IDs = n_elements(file_search(fIDs+".*"))
      nSplit_tab = n_elements(file_search(fTab+".*"))
    endif else begin
      message, 'ERROR: group_ids file ' + sP.simPath + str(sP.snap) + ' does not exist!'
    endelse
  endif
  
  if (nSplit_IDs ne nSplit_tab) then message,'ERROR: different number of ids and tab files'
  
  if (verbose) then $
    print,'Loading subhalo groups from snapshot ('+str(sP.snap)+') in [' + str(nSplit_IDs) + '] files.'
  
  ; counters
  nGroupsTot    = 0L
  nIDsTot       = 0L
  nSubgroupsTot = 0L
  
  skip    = 0L
  skipSub = 0L
  
  ; headers
  h  = { headerTab,             $
               nGroups:        0L,  $
               nGroupsTot:     0L,  $
               nIDs:           0L,  $
               nIDsTot:        0LL, $
               nTask:          0L,  $
               nSubgroups:     0L,  $
               nSubgroupsTot:  0L   $
             }
             
  ; load 0 for header
  fName = fTab + '.' + str(0)
  openr,lun,fName,/GET_LUN
  readu,lun,h
  close,lun
  free_lun,lun
  
  ;if (h.nTask ne nSplit_tab) then begin
    ;print,'WARNING: h.nTask='+str(h.nTask)+' (m='+str(m)+$
    ;      ') differs from number of TAB split files ('+str(nSplit_tab)+'.'
    ;return,0
  ;endif
  
  for i=0,h.nTask-1 do begin
  
    fName = fTab + '.' + str(i)
    
    ; open and read header
    openr,lun,fName,/GET_LUN
    
    readu,lun,h
      
    ; skip loading actual particle IDs
    if keyword_set(skipIDs) then h.nIDsTot = 1
  
    ; add counters and error check
    nGroupsTot    += h.nGroups
    nIDsTot       += h.nIDs
    nSubgroupsTot += h.nSubgroups
    
    ; allocate storage if this is the first iteration
    if (i eq 0) then begin
      sf  = {                                      $
        GroupLen        : ulonarr(h.nGroupsTot)   ,$
        GroupOffset     : ulonarr(h.nGroupsTot)   ,$
        GroupMass       : fltarr(h.nGroupsTot)    ,$
        GroupPos        : fltarr(3,h.nGroupsTot)  ,$
        Group_M_Mean200 : fltarr(h.nGroupsTot)    ,$
        Group_R_Mean200 : fltarr(h.nGroupsTot)    ,$
        Group_M_Crit200 : fltarr(h.nGroupsTot)    ,$
        Group_R_Crit200 : fltarr(h.nGroupsTot)    ,$
        Group_M_TH200   : fltarr(h.nGroupsTot)    ,$
        Group_R_TH200   : fltarr(h.nGroupsTot)    ,$
        Group_VD_M200   : fltarr(h.nGroupsTot)    ,$
        Group_VD_C200   : fltarr(h.nGroupsTot)    ,$
        Group_VD_TH200  : fltarr(h.nGroupsTot)    ,$
        GroupContCount  : ulonarr(h.nGroupsTot)   ,$
        GroupContMass   : fltarr(h.nGroupsTot)    ,$
        GroupNsubs      : ulonarr(h.nGroupsTot)   ,$
        GroupFirstSub   : ulonarr(h.nGroupsTot)   ,$
                                                         $
        SubgroupLen         : ulonarr(h.nSubgroupsTot)  ,$
        SubgroupOffset      : ulonarr(h.nSubgroupsTot)  ,$
        SubgroupParent      : ulonarr(h.nSubgroupsTot)  ,$
        SubgroupMass        : fltarr(h.nSubgroupsTot)   ,$
        SubgroupPos         : fltarr(3,h.nSubgroupsTot) ,$
        SubgroupVel         : fltarr(3,h.nSubgroupsTot) ,$
        SubgroupCM          : fltarr(3,h.nSubgroupsTot) ,$
        SubgroupSpin        : fltarr(3,h.nSubgroupsTot) ,$
        SubgroupVelDisp     : fltarr(h.nSubgroupsTot)   ,$
        SubgroupVmax        : fltarr(h.nSubgroupsTot)   ,$
        SubgroupVmaxRad     : fltarr(h.nSubgroupsTot)   ,$
        SubgroupHalfMassRad : fltarr(h.nSubgroupsTot)   ,$
        SubgroupIDMostBound : ulonarr(h.nSubgroupsTot)  ,$
        SubgroupGrNr        : lonarr(h.nSubgroupsTot)   ,$
                                                         $
        SubgroupIDs         : lon64arr(h.nIDsTot)       ,$
        nGroupsTot          : h.nGroupsTot                ,$
        nSubgroupsTot       : h.nSubgroupsTot             ,$
        nIDsTot             : h.nIDsTot                    $
      }
    endif
    
    ; allocate temporary storage for groups in this part
    if (h.nGroups gt 0) then begin
      part = {Len        : ulonarr(h.nGroups)   ,$
              Offset     : ulonarr(h.nGroups)   ,$
              Mass       : fltarr(h.nGroups)    ,$
              Pos        : fltarr(3,h.nGroups)  ,$
              M_Mean200  : fltarr(h.nGroups)    ,$
              R_Mean200  : fltarr(h.nGroups)    ,$
              M_Crit200  : fltarr(h.nGroups)    ,$
              R_Crit200  : fltarr(h.nGroups)    ,$
              M_TH200    : fltarr(h.nGroups)    ,$
              R_TH200    : fltarr(h.nGroups)    ,$
              ;VD_M200    : fltarr(h.nGroups)    ,$ ;FLAG_Group_VelDisp
              ;VD_C200    : fltarr(h.nGroups)    ,$ ;FLAG_Group_VelDisp
              ;VD_TH200   : fltarr(h.nGroups)    ,$ ;FLAG_Group_VelDisp
              ContCount  : ulonarr(h.nGroups)   ,$
              ContMass   : fltarr(h.nGroups)    ,$
              Nsubs      : ulonarr(h.nGroups)   ,$
              FirstSub   : ulonarr(h.nGroups)    $
              }
    endif
    
    ; load group data from this part
    readu,lun,part
    
    ; fill sf with group data from this part
    sf.GroupLen        [skip:(skip+h.nGroups-1)] = part.Len
    sf.GroupOffset     [skip:(skip+h.nGroups-1)] = part.Offset
    sf.GroupMass       [skip:(skip+h.nGroups-1)] = part.Mass
    sf.GroupPos        [*,skip:(skip+h.nGroups-1)] = part.Pos
    
    sf.Group_M_Mean200 [skip:(skip+h.nGroups-1)] = part.M_Mean200
    sf.Group_R_Mean200 [skip:(skip+h.nGroups-1)] = part.R_Mean200
    sf.Group_M_Crit200 [skip:(skip+h.nGroups-1)] = part.M_Crit200
    sf.Group_R_Crit200 [skip:(skip+h.nGroups-1)] = part.R_Crit200
    sf.Group_M_TH200   [skip:(skip+h.nGroups-1)] = part.M_TH200
    sf.Group_R_TH200   [skip:(skip+h.nGroups-1)] = part.R_TH200
        
    sf.GroupContCount  [skip:(skip+h.nGroups-1)] = part.ContCount
    sf.GroupContMass   [skip:(skip+h.nGroups-1)] = part.ContMass
    sf.GroupNsubs      [skip:(skip+h.nGroups-1)] = part.Nsubs
    sf.GroupFirstsub   [skip:(skip+h.nGroups-1)] = part.FirstSub
        
    skip += h.nGroups
    
    ; allocate temporary storage for subgroups in this part
    if (h.nSubgroups gt 0) then begin
      part = {Len         : ulonarr(h.nSubgroups)   ,$
              Offset      : ulonarr(h.nSubgroups)   ,$
              Parent      : ulonarr(h.nSubgroups)   ,$
              Mass        : fltarr(h.nSubgroups)    ,$
              
              Pos         : fltarr(3,h.nSubgroups)  ,$
              Vel         : fltarr(3,h.nSubgroups)  ,$
              CM          : fltarr(3,h.nSubgroups)  ,$
              Spin        : fltarr(3,h.nSubgroups)    ,$
              
              VelDisp     : fltarr(h.nSubgroups)      ,$
              Vmax        : fltarr(h.nSubgroups)      ,$
              VmaxRad     : fltarr(h.nSubgroups)      ,$
              
              HalfMassRad : fltarr(h.nSubgroups)      ,$
              IDMostBound : ulon64arr(h.nSubgroups)   ,$
              GrNr        : ulonarr(h.nSubgroups)     ,$
              masstab     : fltarr(6,h.nSubgroups)     $
              }
    endif
    
    ; load subgroup data from this part
    readu,lun,part
    
    ; fill sf with Subgroup data from this part
    sf.SubgroupLen         [skipSub:(skipSub+h.nSubgroups-1)]    = part.Len
    sf.SubgroupOffset      [skipSub:(skipSub+h.nSubgroups-1)]    = part.Offset
    sf.SubgroupParent      [skipSub:(skipSub+h.nSubgroups-1)]    = part.Parent
    sf.SubgroupMass        [skipSub:(skipSub+h.nSubgroups-1)]    = part.Mass
    
    sf.SubgroupPos         [*,skipSub:(skipSub+h.nSubgroups-1)]  = part.Pos
    sf.SubgroupVel         [*,skipSub:(skipSub+h.nSubgroups-1)]  = part.Vel
    sf.SubgroupCM          [*,skipSub:(skipSub+h.nSubgroups-1)]  = part.CM
    sf.SubgroupSpin        [*,skipSub:(skipSub+h.nSubgroups-1)]  = part.Spin
    
    sf.SubgroupVelDisp     [skipSub:(skipSub+h.nSubgroups-1)]    = part.VelDisp
    sf.SubgroupVmax        [skipSub:(skipSub+h.nSubgroups-1)]    = part.Vmax
    sf.SubgroupVmaxRad     [skipSub:(skipSub+h.nSubgroups-1)]    = part.VmaxRad
        
    sf.SubgroupHalfMassRad [skipSub:(skipSub+h.nSubgroups-1)]    = part.HalfMassRad
    sf.SubgroupIDMostBound [skipSub:(skipSub+h.nSubgroups-1)]    = part.IDMostBound
    sf.SubgroupGrNr        [skipSub:(skipSub+h.nSubgroups-1)]    = part.GrNr
    ;sf.SubgroupMasstab    [*,skipSub:(skipSub+h.nSubgroups-1)]  = part.masstab
    
    skipSub += h.nSubgroups
    
    ; close file
    close,lun
    free_lun, lun
  
  endfor
  
  ; verify accumulated totals with last header totals
  if (nGroupsTot ne h.nGroupsTot or $
      (nIDsTot ne h.nIDsTot and not keyword_set(skipIDs)) or $
      nSubgroupsTot ne h.nSubgroupsTot) then begin
    print,'ERROR: Totals do not add up.'
    return,0
  endif
  
  if not keyword_set(skipIDs) then begin
    ; IDs header
    hIDs  = { headerIDs,              $
                 nGroups:        0L,  $
                 nGroupsTot:     0L,  $
                 nIDs:           0L,  $
                 nIDsTot:        0LL, $
                 nTask:          0L,  $
                 nPrevIDs:       0L   $
               }
    
    ; IDs counters
    nGroupsTot_IDs = 0L
    nIDsTot_IDs    = 0L
    
    count = 0
    found = 0
    
    ; load 0 for header
    fName = fIDs + '.' + str(0)
    openr,lun,fName,/GET_LUN
    readu,lun,hIDs
    
    if (hIDs.nTask ne nSplit_IDs) then begin
      print,'WARNING: h.nTask='+str(h.nTask)+' (sP.snap='+str(sP.snap)+$
            ') differs from number of IDS split files ('+str(nSplit_IDs)+'.'
      ;return,0
    endif
    
    close,lun
    free_lun,lun
    
    ; load IDs
    for i=0,hIDs.nTask-1 do begin
    
      fName = fIDs + '.' + str(i)
      
      ; open and read header
      openr,lun,fName,/GET_LUN
      
      readu,lun,hIDs
    
      ; add counters and error check
      nGroupsTot_IDs    += hIDs.nGroups
      nIDsTot_IDs       += hIDs.nIDs
      
      ; allocate storage if this is the first iteration
      if (i eq 0) then begin
        subgroupLen = hIDs.nIDsTot ; only for i==0
        ;subgroupIDs = lon64arr(hIDs.nIDsTot)
      endif
      
      ; determine number of IDs to read for this subgroup
      if (count lt hIDs.nPrevIDs + hIDs.nIDs) then begin
        nSkip       = count - hIDs.nPrevIDs
        nRemaining  = hIDs.nPrevIDs + hIDs.nIDs - count
        
        if (subgroupLen gt nRemaining) then begin
          nToRead = nRemaining
        endif else begin
          nToRead = subgroupLen
        endelse
      endif
      
      ; read IDs for this subgroup
      if (nToRead gt 0) then begin
        if (nSkip gt 0) then begin
          dummy = lon64arr(nSkip)
          readu,lun,dummy
          print,dummy
          return,0 ;die
        endif
        
        ; read IDs for this file part
        partIDs = ulon64arr(hIDs.nIDs)
        readu,lun,partIDs
        
        ; fill sf.SubgroupIDs with this part
        sf.SubgroupIDs[hIDs.nPrevIDs : (hIDs.nPrevIDs + hIDs.nIDs - 1)] = partIDs
        
        ; update counters
        count       += nToRead
        subgroupLen -= nToRead
      endif
      
      ; close file
      close,lun
      free_lun, lun    
      
    endfor
    
    ; verify accumulated totals with last header totals
    if (nGroupsTot_IDs ne hIDs.nGroupsTot or nIDsTot_IDs ne hIDs.nIDsTot) then begin
      print,'ERROR: IDs totals do not add up.'
      return,0
    endif
  
  endif ;skipIDs
  
  if (verbose) then $
    print,'Load complete. (nGroupsTot = ' + str(nGroupsTot) + ' nSubgroupsTot = ' + str(nSubgroupsTot) + $
          ' nIDsTot = ' + str(nIDsTot) + ')'
        
  return,sf
end

; loadHsmlDir(): load (OLD, not HDF5) "hsmldir" with Hsml,Density,VelDisp of -all- particles
;                may be ID ordered or "IC ID" (effectively snapshot index) ordered

function loadHsmlDir, sP=sP, verbose=verbose, partType=partType, $
                      readHsml=readHsml, readDens=readDens, readVelDisp=readVelDisp

  if not keyword_set(verbose) then verbose = 0
  if ~keyword_set(readHsml) and ~keyword_set(readDens) and ~keyword_set(readVelDisp) then $
    message,'Error: Need to specify at least one field to read.'

  ; set filename
  ext = string(sP.snap,format='(I3.3)')
  fName = sP.simPath + 'hsmldir_' + ext + '/hsml_' + ext
  
  ; check existance and multiple outputs
  if not file_test(fName) then begin
    if (file_test(fName+'.0')) then begin
      ; split into multiples, get count
      nSplit = n_elements(file_search(fName+".*"))
    endif else begin
      message, 'ERROR: hsmldir file ' + sP.simPath + str(sP.snap) + ' does not exist!'
    endelse
  endif
  
  if (verbose) then $
    print,'Loading hsmldir from snapshot ('+str(sP.snap)+') in [' + str(nSplit) + '] files.'
  
  ; counters
  nHSMLTot = 0L
  skip     = 0L
  
  ; headers
  h  = { nHSML:          0L,  $
         SendOffset:     0L,  $
         nTotal:        0LL,  $
         nTask:          0L   }
             
  ; load 0 for header
  openr,lun,fName+'.0',/GET_LUN
  readu,lun,h
  close,lun
  free_lun,lun
  
  if (h.nTask ne nSplit) then $
     print,'WARNING: h.nTask='+str(h.nTask)+' differs from number of split files ('+str(nSplit)
  
  for i=0,h.nTask-1 do begin
  
    fNameCur = fName + '.' + str(i)
    
    ; open and read header
    openr,lun,fNameCur,/GET_LUN
    readu,lun,h
      
    ; skip loading actual particle IDs
    if keyword_set(skipIDs) then h.nIDsTot = 1
  
    ; add counters and error check
    nHSMLTot    += h.nHSML
    
    ; allocate storage if this is the first iteration
    if (i eq 0) then begin
      if keyword_set(readHsml) and keyword_set(readDens) and keyword_set(readVelDisp) then $
        hd = { hsml : fltarr(h.nTotal), density : fltarr(h.nTotal), veldisp : fltarr(h.nTotal) }
      if keyword_set(readHsml) and keyword_set(readDens) and ~keyword_set(readVelDisp) then $
        hd = { hsml : fltarr(h.nTotal), density : fltarr(h.nTotal) }
      if keyword_set(readHsml) and ~keyword_set(readDens) and keyword_set(readVelDisp) then $
        hd = { hsml : fltarr(h.nTotal), veldisp : fltarr(h.nTotal) }
      if ~keyword_set(readHsml) and keyword_set(readDens) and keyword_set(readVelDisp) then $
        hd = { density : fltarr(h.nTotal), veldisp : fltarr(h.nTotal) }
      if ~keyword_set(readHsml) and ~keyword_set(readDens) and keyword_set(readVelDisp) then $
        hd = { veldisp : fltarr(h.nTotal) }
      if keyword_set(readHsml) and ~keyword_set(readDens) and ~keyword_set(readVelDisp) then $
        hd = { hsml : fltarr(h.nTotal) }
      if ~keyword_set(readHsml) and keyword_set(readDens) and ~keyword_set(readVelDisp) then $
        hd = { density : fltarr(h.nTotal) }
    endif
    
    ; allocate temporary storage for groups in this part
    if (h.nHSML gt 0) then $
      part = { hsml : fltarr(h.nHSML), density : fltarr(h.nHSML), veldisp : fltarr(h.nHSML) }
    
    ; load group data from this part
    readu,lun,part
    
    ; fill hd with group data from this part
    if keyword_set(readHsml) then hd.hsml[skip:(skip+h.nHSML-1)]       = part.hsml
    if keyword_set(readDens) then hd.density[skip:(skip+h.nHSML-1)]    = part.density
    if keyword_set(readVelDisp) then hd.veldisp[skip:(skip+h.nHSML-1)] = part.veldisp
    
    skip += h.nHSML
    
    ; close file
    close,lun
    free_lun, lun
  
  endfor
  
  ; verify accumulated totals with last header totals
  if nHSMLTot ne h.nTotal then message,'ERROR: Totals do not add up.'
  
  ; restrict to a particle type if requested
  if keyword_set(partType) then begin
    if partType eq 'dm' then begin
      ; ID ordered assumption
      ids_dm    = loadSnapshotSubset(sP=sP,partType='dm',field='ids')
      ids_gas   = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
      ids_stars = loadSnapshotSubset(sP=sP,partType='stars',field='ids')
      
      ids = [ids_dm,ids_gas,ids_stars] ; concat
      ids_gas   = !NULL
      ids_stars = !NULL
      
      idsMask = bytarr(n_elements(ids)) ; mask to keep track of DM
      idsMask[0:n_elements(ids_dm)-1] = 1B ; DM first
      ids_sort = sort(ids) ; ID sort ascending
      idsMask = idsMask[ids_sort] ; sort mask
      ids = !NULL
      ids_sort = !NULL
      
      dm_hsml_inds = where(idsMask eq 1B,count) ; dm locations in hsmldir (ascending id order)
      dm_hsml_inds = dm_hsml_inds[ids_dm-min(ids_dm)] ; put back in snapshot order
      
      if count ne n_elements(ids_dm) then message,'Error: Bad location search.'
      idsMask = !NULL
      ids_dm  = !NULL
      
      ; if just one field requested, return it as an array
      if ~keyword_set(readHsml) and ~keyword_set(readDens) and keyword_set(readVelDisp) then $
        return, hd.veldisp[dm_hsml_inds]
      if keyword_set(readHsml) and ~keyword_set(readDens) and ~keyword_set(readVelDisp) then $
        return, hd.hsml[dm_hsml_inds]
      if ~keyword_set(readHsml) and keyword_set(readDens) and ~keyword_set(readVelDisp) then $
        return, hd.density[dm_hsml_inds]
      
      ; otherwise modify hd and return
      if keyword_set(readHsml) and keyword_set(readDens) and keyword_set(readVelDisp) then $
        hd2 = { hsml : hd.hsml[dm_hsml_inds], density : hd.density[dm_hsml_inds], veldisp : hd.veldisp[dm_hsml_inds] }
      if keyword_set(readHsml) and keyword_set(readDens) and ~keyword_set(readVelDisp) then $
        hd2 = {  hsml : hd.hsml[dm_hsml_inds], density : hd.density[dm_hsml_inds] }
      if keyword_set(readHsml) and ~keyword_set(readDens) and keyword_set(readVelDisp) then $
        hd2 = {  hsml : hd.hsml[dm_hsml_inds], veldisp : hd.veldisp[dm_hsml_inds] }
      if ~keyword_set(readHsml) and keyword_set(readDens) and keyword_set(readVelDisp) then $
        hd2 = { density : hd.density[dm_hsml_inds], veldisp : hd.veldisp[dm_hsml_inds] }
        
      return,hd2
    endif
    message,'Error: Other part type returns not generalized yet.'
  endif
  
  if (verbose) then print,'Load complete for all particle types. (nHSMLTot = ' + str(nHSMLTot) + ')'
        
  return,hd
end

; getSnapFilelist(): take input path and snapshot number and find snapshot filename

function getSnapFilelist, fileBase, snapNum=m, subBox=subBox

  sbstr = ''
  if keyword_set(subBox) then sbstr = 'subbox_'

  ; format snapNum and initial guess
  if (str(m) eq 'none') then begin
    ext = ''
    f = fileBase
    if file_test(f) then return,[f]
  endif else begin
    ; check for '/' on end of fileBase
    lastChar = strmid(fileBase,strlen(fileBase)-1,1)
    if (lastChar ne '/') then fileBase += '/'
    
    ; format snapshot number
    if (m le 999) then $
      ext = string(m,format='(I3.3)')
    if (m gt 999) then $
      ext = string(m,format='(I4.4)')
      
    f = fileBase + 'snapdir_' + ext + '/snap_' + sbstr + ext
  endelse

  ; check for single (non-split)
  if file_test(fileBase+'snap_'+sbstr+ext+'.hdf5') then $
    return, [fileBase+'snap_'+sbstr+ext+'.hdf5']

  ; check for single groupordered
  if file_test(fileBase+'snap-groupordered_'+ext+'.hdf5') then $
    return, [fileBase+'snap-groupordered_'+ext+'.hdf5']
  
  ; check for multiple groupordered
  if file_test(fileBase+'snapdir_'+ext+'/snap-groupordered_'+ext+'.0.hdf5') then begin
    ff = file_search(fileBase + 'snapdir_' + ext + '/snap-groupordered_' + ext+".*hdf5")
    return, fileBase + 'snapdir_' + ext + '/snap-groupordered_' + ext + '.' + $
            str(indgen(n_elements(ff))) + '.hdf5'
  endif
  
  ; check for exact name or multiple outputs
  if file_test(f+'.hdf5') then return, [f+'.hdf5']
  
  ; check for multiple outputs (general cosmo case)
  if file_test(f+'.0.hdf5') then begin
    ff = file_search(f+".*hdf5")
    return, f + '.' + str(indgen(n_elements(ff))) + '.hdf5'
  endif

  message,'Error: Failed to find snapshot.'
end

; loadSnapshotHeader(): load header

function loadSnapshotHeader, sP=sP, verbose=verbose, subBox=subBox, fileName=fileName

  if not keyword_set(verbose) then verbose = 0

  ; get matching filename
  if n_elements(fileName) eq 0 then $
    fileList = getSnapFilelist(sP.simPath,snapNum=sP.snap,subBox=subBox)
  if n_elements(fileName) gt 0 then $
    fileList = [fileName]

  ; read header from first part
  fileID   = h5f_open(fileList[0])
  s = h5_parse(fileID,"Header")
  
  ; fill header struct
  h = {                                                          $
        nPartThisFile       : s.numPart_ThisFile._DATA          ,$
        nPartTot            : s.numPart_Total._DATA             ,$
        nPartTotHighword    : s.numPart_Total_Highword._DATA    ,$
        massTable           : s.massTable._DATA                 ,$
        time                : s.time._DATA                      ,$
        redshift            : s.redshift._DATA                  ,$
        boxSize             : s.boxSize._DATA                   ,$
        numFilesPerSnapshot : s.numFilesPerSnapshot._DATA       ,$
        Omega0              : s.Omega0._DATA                    ,$
        OmegaLambda         : s.OmegaLambda._DATA               ,$
        hubbleParam         : s.hubbleParam._DATA               ,$
        flagSFR             : s.flag_SFR._DATA                  ,$
        flagCooling         : s.flag_Cooling._DATA              ,$
        flagStellarAge      : s.flag_StellarAge._DATA           ,$
        flagMetals          : s.flag_Metals._DATA               ,$
        flagFeedback        : s.flag_Feedback._DATA             ,$
        flagDoublePrecision : s.flag_DoublePrecision._DATA       $
        ;flagICInfo          : s.flag_IC_Info._DATA              $
        ;compVectorLen      : s.composition_vector_length        $
      }
      
  h5f_close, fileID
  
  return, h
end

; readStrategyIO(): determine the blocking read strategy given the requested read pattern

function readStrategyIO, inds=inds, indRange=indRange, $
                         nPartTot=nPartTot, partType=partType, verbose=verbose, $
                         sort_inds=sort_inds, sorted_inds=sorted_inds ; output for MultiBlock

  if keyword_set(inds) and keyword_set(indRange) then message,'Error: Not both.'
  
  ; no strategy: whole data field for all particles
  if ~keyword_set(inds) and ~keyword_set(indRange) then begin
    
    rs = { nBlocks : 1, blockStart: [0LL], blockEnd: [nPartTot[partType]-1LL] }
    
    return, rs
  endif
  
  ; indRange: whole block from min to max
  if keyword_set(indRange) then begin
    if indRange[0] ge indRange[1] or indRange[0] lt 0 or $
       indRange[1] ge nPartTot[partType] then message,'Error: Bad indRange.'
       
    rs = { nBlocks     : 1                      ,$
           nTotRead    : 0LL                    ,$
           blockStart  : [long(indRange[0])]    ,$ ; snapshot index start
           blockEnd    : [long(indRange[1])]    ,$ ; snapshot index end
           indexMin    : [0LL]                  ,$ ; unused
           indexMax    : [0LL]                   } ; unused
           
    rs.nTotRead = rs.blockEnd[0] - rs.blockStart[0] + 1
    return, rs
  endif
  
  ; for simple or multiblock we need to sort the indices to load in contiguous chunks
  sort_inds = calcSort(inds)
  sorted_inds = inds[sort_inds]
         
  ; simple: whole block from min to max, then take subset
  rsSimple = { nBlocks     : 1                      ,$
               nTotRead    : 0LL                    ,$
               blockStart  : [min(inds)]            ,$ ; snapshot index start
               blockEnd    : [max(inds)]            ,$ ; snapshot index end
               indexMin    : [0LL]                  ,$ ; start of subset of inds handled
               indexMax    : [n_elements(inds)-1LL]  } ; end of subset of inds handled
                     
  rsSimple.nTotRead = rsSimple.blockEnd[0] - rsSimple.blockStart[0] + 1
  efficiencyFracSimple = n_elements(inds) / float(rsSimple.nTotRead)
    
  if verbose then $
    print,'I/O Strategy: [SimpleRead] reading ['+str(rsSimple.nTotRead)+$
      ' / '+str(nPartTot[partType])+'] fracOfSnap: '+$
      string(rsSimple.nTotRead*100.0/nPartTot[partType],format='(f5.2)')+'% efficiencyFrac: '+$
      string(efficiencyFracSimple*100,format='(f5.2)')+'%'    
    
  ; if efficiency of the simple read is 40% or greater just go for that
  if efficiencyFracSimple ge 0.4 then $
    return, rsSimple
 
  ; lowN: for an extremely small number of indices (<1000 or so) just read those exact elements
  ; TODO
      
  ; multiBlock: analyze sparsity of indices, decompose read into a reasonable number of blocks
  ;             such that we don't read too much unwanted data, but also do enough sequential
  ;             I/O and not too much random I/O
  binSize = round(sqrt(max(inds))/10000.0) * 10000 > 10000 < 10000000 ; adjust to range
      
  hist = histogram(inds,binsize=binsize,min=0,loc=loc)
  occBlocks = where(hist ne 0,count)
      
  nBlocks = 0
  blockStart = []
  blockLen   = []
  prevBlock = ulong64(0)
  count = 0
      
  ; search for contiguous id blocks
  for i=0,n_elements(occBlocks)-1 do begin
    count += 1
    
    ; jump detected (zero block length ge binsize)
    if occBlocks[i] ne (prevBlock+1) then begin
      nBlocks += 1
      blockStart = [blockStart,loc[occBlocks[i]]]
      blockLen   = [blockLen,count]
      if i ne n_elements(occBlocks)-1 then count = 0
    endif
    
    prevBlock = occBlocks[i]
  endfor
      
  ; remove first (zero) block length and add final
  count += 1
  if nBlocks gt 1 then blockLen = [blockLen[1:*],count]
  if nBlocks eq 1 then blockLen = [count]
  if nBlocks le 0 then message,'Badness. Revert to simple?'
  
  ; determine rest of block structure
  rs = { nBlocks     : nBlocks                                ,$
         nTotRead    : round(total(blockLen * binSize))       ,$
         blockStart  : blockStart                             ,$ ; snapshot index start
         blockEnd    : blockStart + blockLen * binSize - 1    ,$ ; snapshot index end
         indexMin    : lon64arr(nBlocks)                      ,$ ; start of subset of inds handled
         indexMax    : lon64arr(nBlocks)                       } ; end of subset of inds handled
  
  for i=0,nBlocks-1 do begin
    w = where(sorted_inds ge rs.blockStart[i] and sorted_inds le rs.blockEnd[i],count)
    if count eq 0 then message,'Error: MultiBlock failure.'
    rs.indexMin[i] = w[0]
    rs.indexMax[i] = w[count-1]
  endfor
  
  efficiencyFracMulti = n_elements(inds) / float(rs.nTotRead)
  fracOfSnap = rs.nTotRead*100.0/nPartTot[partType]
      
  if fracOfSnap gt 90 then return, rsSimple ; just a simple read if we already do most of the data
      
  if verbose then $
    print,'I/O Strategy: [MultiBlock] reading ['+str(rs.nTotRead)+$
      ' / '+str(nPartTot[partType])+'] fracOfSnap: '+$
      string(fracOfSnap,format='(f5.2)')+'% efficiencyFrac: '+$
      string(efficiencyFracMulti*100,format='(f5.2)')+'% (nBlocks = '+str(rs.nBlocks)+')'    

  return, rs
  
end

; snapshotFieldExists(): check if a specified field exists in the snapshot
function snapshotFieldExists, sP=sP, fileName=fileName, field=field

  if n_elements(fileName) eq 0 then $
    fileList = getSnapFilelist(sP.simPath,snapNum=sP.snap,subBox=subBox)
  if n_elements(fileName) gt 0 then $
    fileList = [fileName]

  ; input config: set partType number if input in string
  ;partType = strlowcase(string(PT)) ; so we don't change the input
  ;partType = partTypeNum(partType)
  
  ; parse for this field name (over all particle types)
  hdf5s = h5_parse(fileList[0]) ;structure only
  existsFlag = tag_exist(hdf5s,field)

  return, existsFlag
end

; loadSnapshotSubset(): for a given snapshot load only one field for one particle type
;                       partType = [0,1,2,4] or ('gas','dm','tracer','stars') (case insensitive)
;                       field    = ['ParticleIDs','coordinates','xyz',...] (case insensitive)
;                       inds     (optional) : known indices requested, optimize the load
;                       indRange (optional) : same, but specify only min and max indices
; - specify either sP (with .snap) or a direct fileName

function loadSnapshotSubset, sP=sP, fileName=fileName, partType=PT, field=field, verbose=verbose, $
                             inds=inds, indRange=indRange, doublePrec=doublePrec, subBox=subBox

  if not keyword_set(verbose) then verbose = 0

  if n_elements(fileName) eq 0 then $
    fileList = getSnapFilelist(sP.simPath,snapNum=sP.snap,subBox=subBox)
  if n_elements(fileName) gt 0 then $
    fileList = [fileName]

  nFiles = n_elements(fileList)

  ; input config: set partType number if input in string
  partType = strlowcase(string(PT)) ; so we don't change the input
  partType = partTypeNum(partType)
  
  ; load particle array sizes from header of first part
  fileID   = h5f_open(fileList[0])
  headerID = h5g_open(fileID,"Header")
  nPartTot = h5a_read(h5a_open_name(headerID,"NumPart_Total"))
  nSplits  = h5a_read(h5a_open_name(headerID,"NumFilesPerSnapshot"))
  flagDbl  = h5a_read(h5a_open_name(headerID,"Flag_DoublePrecision"))
  h5g_close, headerID
  h5f_close, fileID
  
  ; parse for new GFM related fields
  hdf5s = h5_parse(fileList[0]) ;structure only
  GFMExistsFlag  = tag_exist(hdf5s,'GFM_Metallicity')
  if GFMExistsFlag then begin
    gfmNumElements = hdf5s.PartType0.GFM_Metals._DIMENSIONS[0]
    if PT eq 'star' or PT eq 'stars' and hdf5s.header.numPart_Total._DATA[partType] gt 0 then $
      gfmNumPhotometrics = hdf5s.PartType4.GFM_StellarPhotometrics._DIMENSIONS[0]
  endif
  
  if (nSplits ne nFiles) then $
    message,'ERROR: NumFilesPerSnapshot ['+str(nSplits)+'] differs from number of split files found ['+$
           str(nFiles)+'].'
  
  ; double precision
  if (flagDbl eq 1 and not keyword_set(doublePrec)) then $
    print,'Warning: Snapshot is double precision but only singlePrec load requested.'
  if (flagDbl eq 0 and keyword_set(doublePrec)) then $
    print,'Warning: Snapshot is single precision but doublePrec load requested.'
  
  ; early exit: no particles of requested type
  if (nPartTot[partType] eq 0) then begin
    print,'Warning: No particles of requested type present.'
    return,[]
  endif
  
  ; input config: set fieldName and return array
  count = 0L
  field = strlowcase(field)
  rDims = 1 ;override if needed
  fieldName = ''

  ; common fields (all besides tracersMC)
  ; -------------------------------------
  
  if (field eq 'coordinates' or field eq 'xyz' or field eq 'positions' or field eq 'pos' or $
      field eq 'x' or field eq 'y' or field eq 'z') then begin
    fieldName = 'Coordinates'
    rDims = 3
    rType = 'float'
  endif
  if (field eq 'particleids' or field eq 'ids') then begin
    rType = 'long'
    fieldName = 'ParticleIDs'
    if partType eq partTypeNum('tracerMC') then message,'Error: tracersMC use tracerids not ids'
  endif
  if (field eq 'potential' or field eq 'phi') then begin
    rType = 'float'
    fieldName = 'Potential'
  endif
  if (field eq 'velocities' or field eq 'vel' or $
      field eq 'velx' or field eq 'vely' or field eq 'velz') then begin
    fieldName = 'Velocities'
    rDims = 3
    rType = 'float'
  endif
  
  ; gas only
  ; --------
  if (field eq 'center_of_mass' or field eq 'centerofmass' or field eq 'com' or field eq 'cm' or $
      field eq 'cmx' or field eq 'cmy' or field eq 'cmz') then begin
    fieldName = 'Center-of-Mass'
    rDims = 3
    rType = 'float'
    if (partType ne 0) then message,'Error: CoM is gas only!'
  endif
  if (field eq 'coolingrate' or field eq 'coolrate') then begin
    rType = 'float'
    fieldName = 'CoolingRate'
    if (partType ne 0) then message,'Error: CoolingRate is gas only!'
    if GFMExistsFlag then message,'Error: Old field. Use GFM_CoolingRate.'
  endif
  if (field eq 'density' or field eq 'rho' or field eq 'dens') then begin
    rType = 'float'
    fieldName = 'Density'
    if (partType ne 0) then message,'Error: Density is gas only!'
  endif
  if (field eq 'electronabundance' or field eq 'ne' or field eq 'nelec') then begin
    rType = 'float'
    fieldName = 'ElectronAbundance'
    if (partType ne 0) then message,'Error: NE is gas only!'
  endif
  if (field eq 'internalenergy' or field eq 'u') then begin
    rType = 'float'
    fieldName = 'InternalEnergy'
    if (partType ne 0) then message,'Error: U is gas only!'
  endif
  if (field eq 'machnumber' or field eq 'machnum') then begin
    rType = 'float'
    fieldName = 'MachNumber'
    if (partType ne 0) then message,'Error: MachNum is gas only!'
  endif
  if (field eq 'maxfaceangle') then begin
    rType = 'float'
    fieldName = 'MaxFaceAngle'
    if (partType ne 0) then message,'Error: MaxFaceAngle is gas only!'
  endif
  if (field eq 'neutralhydrogenabundance' or field eq 'nh') then begin
    rType = 'float'
    fieldName = 'NeutralHydrogenAbundance'
    if (partType ne 0) then message,'Error: NH is gas only!'
  endif
  if (field eq 'number_of_faces_of_cell' or field eq 'num_cell_faces' or field eq 'numcellfaces' or $
      field eq 'numfaces') then begin
    rType = 'float'
    fieldName = 'Number of faces of cell'
    if (partType ne 0) then message,'Error: NumFaces is gas only!'
  endif
  if (field eq 'pressure' or field eq 'pres') then begin
    rType = 'float'
    if (partType ne 0) then message,'Error: Pressure is gas only!'
    fieldName = 'Pressure'
  endif
  if (field eq 'smoothinglength' or field eq 'hsml') then begin
    rType = 'float'
    fieldName = 'SmoothingLength'
    if (partType ne 0) then message,'Error: HSML is gas only!'
  endif
  if (field eq 'starformationrate' or field eq 'sfr') then begin
    rType = 'float'
    fieldName = 'StarFormationRate'
    if (partType ne 0) then message,'Error: SFR is gas only!'
  endif
  if (field eq 'surface_area' or field eq 'surfarea') then begin
    rType = 'float'
    fieldName = 'Surface Area'
    if (partType ne 0) then message,'Error: SurfArea is gas only!'
  endif
  if (field eq 'volume' or field eq 'vol') then begin
    rType = 'float'
    fieldName = 'Volume'
    if (partType ne 0) then message,'Error: Vol is gas only!'
  endif
  
  ; newer
  if (field eq 'gfm_coolingrate' or field eq 'gfm_coolrate') then begin
    rType = 'float'
    fieldName = 'GFM_CoolingRate'
    if (partType ne 0) then message,'Error: GFM_CoolingRate is gas only!'
    if ~GFMExistsFlag then message,'Error: Snapshot does not have GFM fields!'
  endif
  if (field eq 'gfm_winddmveldisp') then begin
    rType = 'float'
    fieldName = 'GFM_WindDMVelDisp'
    if (partType ne 0) then message,'Error: GFM_WindDMVelDisp is gas only!'
    if ~GFMExistsFlag then message,'Error: Snapshot does not have GFM fields!'
  endif
  if (field eq 'gfm_agnradiation' or field eq 'gfm_agnrad') then begin
    rType = 'float'
    fieldName = 'GFM_AGNRadiation'
    if (partType ne 0 and partType ne 5) then message,'Error: HostHaloMass is gas/BHs only!'
  endif
  
  ; gas/stars only
  ; --------------
  if (field eq 'metallicity' or field eq 'metal') then begin ; OLD
    rType = 'float'
    fieldName = 'Metallicity'
    if (partType ne 0 and partType ne 4) then message,'Error: Z is gas/stars only!'
    if GFMExistsFlag then fieldName = 'GFM_Metallicity' ; use GFM field
  endif
  if (field eq 'gfm_metallicity' or field eq 'gfm_z') then begin
    rType = 'float'
    fieldName = 'GFM_Metallicity'
    if (partType ne 0 and partType ne 4) then message,'Error: GFM_Metallicity is gas/stars only!'
    if ~GFMExistsFlag then message,'Error: Snapshot does not have GFM fields!'
  endif
  if (field eq 'gfm_metals') then begin
    rType = 'float'
    rDims = gfmNumElements
    fieldName = 'GFM_Metals'
    if (partType ne 0 and partType ne 4) then message,'Error: GFM_Metals is gas/stars only!'
    if ~GFMExistsFlag then message,'Error: Snapshot does not have GFM fields!'
  endif
  
  ; gas/stars/bh only
  ; --------------
  if (field eq 'masses' or field eq 'mass') then begin
    rType = 'float'
    fieldName = 'Masses'
    if (partType ne 0 and partType ne 4 and partType ne 5) then message,'Error: Mass is gas/stars/BH only!'
  endif
  if (field eq 'numtr' or field eq 'numtracers') then begin
    rType = 'int'
    fieldName = 'NumTracers'
    if (partType ne 0 and partType ne 4 and partType ne 5) then message,'Error: NumTracers is gas/stars/BH only!'
  endif
  
  ; stars only
  ; ----------
  if (field eq 'gfm_initialmass' or field eq 'gfm_inimass') then begin
    rType = 'float'
    fieldName = 'GFM_InitialMass'
    if (partType ne 4) then message,'Error: GFM_InitialMass is stars only!'
    if ~GFMExistsFlag then message,'Error: Snapshot does not have GFM fields!'
  endif
  if (field eq 'gfm_stellarformationtime' or field eq 'gfm_sftime') then begin
    rType = 'float'
    fieldName = 'GFM_StellarFormationTime'
    if (partType ne 4) then message,'Error: GFM_StellarFormationTime is stars only!'
    if ~GFMExistsFlag then message,'Error: Snapshot does not have GFM fields!'
  endif
  if (field eq 'gfm_stellarphotometrics' or field eq 'gfm_sphot') then begin
    rType = 'float'
    rDims = gfmNumPhotometrics
    fieldName = 'GFM_StellarPhotometrics'
    if (partType ne 4) then message,'Error: GFM_StellarPhotometrics is stars only!'
    if ~GFMExistsFlag then message,'Error: Snapshot does not have GFM fields!'
  endif
  
  ; bh only
  ; -------
  if (field eq 'bh_cumegyinjection_qm' or field eq 'cumegyinjection_qm') then begin
    rType = 'float'
    fieldName = 'BH_CumEgyInjection_QM'
  endif
  if (field eq 'bh_cummassgrowth_qm' or field eq 'cummassgrowth_qm') then begin
    rType = 'float'
    fieldName = 'BH_CumMassGrowth_QM'
  endif
  if (field eq 'bh_coolingluminosity' or field eq 'coolingluminosity') then begin
    rType = 'float'
    fieldName = 'BH_CoolingLuminosity'
  endif
  if (field eq 'bh_density' or field eq 'bh_dens' or field eq 'bh_rho') then begin
    rType = 'float'
    fieldName = 'BH_Density'
  endif
  if (field eq 'bh_hsml' or field eq 'bh_smoothinglength') then begin
    rType = 'float'
    fieldName = 'BH_HSML'
  endif
  if (field eq 'bh_mass' or field eq 'bh_masses') then begin
    rType = 'float'
    fieldName = 'BH_Mass'
  endif
  if (field eq 'bh_mass_bubbles') then begin
    rType = 'float'
    fieldName = 'BH_Mass_Bubbles'
  endif
  if (field eq 'bh_mass_ini') then begin
    rType = 'float'
    fieldName = 'BH_Mass_Ini'
  endif
  if (field eq 'bh_mdot') then begin
    rType = 'float'
    fieldName = 'BH_MDot'
  endif
  if (field eq 'bh_mdotradio' or field eq 'bh_mdot_radio') then begin
    rType = 'float'
    fieldName = 'BH_MDotRadio'
  endif
  if (field eq 'bh_pressure' or field eq 'bh_pres') then begin
    rType = 'float'
    fieldName = 'BH_Pressure'
  endif
  if (field eq 'bh_prog' or field eq 'bh_progs' or field eq 'bh_progenitors') then begin
    rType = 'float'
    fieldName = 'BH_Progs'
  endif
  if (field eq 'bh_u' or field eq 'bh_energy') then begin
    rType = 'float'
    fieldName = 'BH_U'
  endif
  if (field eq 'bh_hosthalomass' or field eq 'hosthalomass') then begin
    rType = 'float'
    fieldName = 'HostHaloMass'
    if partType ne 5 then message,'Error: Hosthalo mass is BH only!'
  endif
  
  if (strmid(field,0,3) eq 'bh_' and partType ne 5) then message,'Error: BH fields are BH only!'
  
  ; subfind SAVE_HSML_IN_SNAPSHOT
  ; -------
  if (field eq 'subfind_dens' or field eq 'subfind_density') then begin
    rType = 'float'
    fieldName = 'SubfindDensity'
  endif
  if (field eq 'subfind_hsml' or field eq 'subfind_smoothinglength') then begin
    rType = 'float'
    fieldName = 'SubfindHSML'
  endif
  if (field eq 'subfind_veldisp' or field eq 'subfind_dmveldisp') then begin
    rType = 'float'
    fieldName = 'SubfindVelDisp'
  endif
  
  ; tracers (Monte Carlo)
  ; ---------------------
  if (field eq 'parentid' or field eq 'parentids') then begin
    rType = 'long'
    fieldName = 'ParentID'
    if (partType ne 3) then message,'Error: ParentID is tracerMC only!'
  endif
  if (field eq 'tracerid' or field eq 'tracerids') then begin
    rType = 'long'
    fieldName = 'TracerID'
    if (partType ne 3) then message,'Error: TracerID is tracerMC only!'
  endif
  
  ; tracers (common)
  ; ----------------
  if (field eq 'properties' or field eq 'quants' or field eq 'quantities' or $
      field eq 'tracer_maxtemp' or field eq 'tracer_maxtemp_time' or field eq 'tracer_maxtemp_dens' or $
      field eq 'tracer_maxdens' or field eq 'tracer_maxdens_time' or field eq 'tracer_maxmachnum' or $
      field eq 'tracer_maxent' or field eq 'tracer_maxent_time' or $
      field eq 'tracer_laststartime' or field eq 'tracer_exchdist' or field eq 'tracer_exchdisterr') then begin
    fieldName = 'FluidQuantities'
    rDims = 5 ; must match to setup in Arepo run if loading all at once (fine to take a slice of one quantity)
    rType = 'float'
    if (partType ne 3 and partType ne 2) then message,'Error: Fluid quantities are tracerMC/Vel only!'
  endif
  
  if (field eq 'tracer_windcounter' or field eq 'tracer_exchcounter') then begin
    fieldName = 'FluidQuantities'
    rDims = 5
    rType = 'int'
    if (partType ne 3 and partType ne 2) then message,'Error: Fluid quantities are tracerMC/Vel only!'
  endif
  
  if field eq 'tracer_maxent' then print,'WARNING: a3inv factor!' ; temp warning
  
  if (fieldName eq '') then begin
    print,'ERROR: Requested field -- ' + strlowcase(field) + ' -- not recognized!'
    stop
  endif
  
  ; I/0 planning
  ; ------------
  
  ; multidim slice (hyperslab selection) request
  multiDimSliceFlag = 0
  if (field eq 'x' or field eq 'y' or field eq 'z' or $
      field eq 'velx' or field eq 'vely' or field eq 'velz' or $
      field eq 'cmx' or field eq 'cmy' or field eq 'cmz' or $
      field eq 'tracer_maxtemp' or field eq 'tracer_maxtemp_time' or field eq 'tracer_maxtemp_dens' or $
      field eq 'tracer_maxdens' or field eq 'tracer_maxdens_time' or field eq 'tracer_maxmachnum' or $
      field eq 'tracer_maxent' or field eq 'tracer_maxent_time' or field eq 'tracer_laststartime' or $
      field eq 'tracer_windcounter' or field eq 'tracer_exchcounter' or field eq 'tracer_exchdist' or $
      field eq 'tracer_exchdisterr') then begin
    multiDimSliceFlag = 1
    
    ; override rDims to one
    rDims = 1
  endif

  ; decide I/O strategy if we have known indices to read (simple, lowN, complex)
  readStrategy = readStrategyIO(inds=inds, indRange=indRange, nPartTot=nPartTot, $
                                partType=partType, verbose=verbose, $
                                sort_inds=sort_inds, sorted_inds=sorted_inds)
    
  ; decide size of return array
  if keyword_set(inds) or keyword_set(indRange) then $
    readSize = readStrategy.nTotRead $
  else $
    readSize = nPartTot[partType]
  
  ; use rType to make return array
  if rType eq 'float' then begin
    ; double precision requested?
    if keyword_set(doublePrec) then $
      r = dblarr(rDims,readSize) $
    else $
      r = fltarr(rDims,readSize)
  endif
  
  if rType eq 'int'   then r = intarr(rDims,readSize)
  
  if rType eq 'long'  then begin
    if fieldName eq 'ParticleIDs' or fieldName eq 'ParentID' or fieldName eq 'TracerID' then begin
      ; IDs are 64bit?
      longIDsBits = h5_parse(fileList[0]) ; just parse full structure
      longIDsBits = longIDsBits.partType0.particleIDs._precision
      if longIDsBits eq 32 then r = lonarr(rDims,readSize)
      if longIDsBits eq 64 then r = lon64arr(rDims,readSize)
      if longIDsBits ne 32 and longIDsBits ne 64 then message,'Error: Unexpected IDs precision.'
    endif else begin
      ; non-ID long field
      r = lonarr(rDims,readSize)
    endelse
  endif
  
  r = reform(r)
  
  if verbose then $
    print,'Loading "' + str(fieldName) + '" for partType=' + str(partType) + ' from snapshot (' + $
          str(sP.snap) + ') in [' + str(nFiles) + '] files. (numPartTot=' + $
          str(nPartTot[partType]) + ' nReadSize=' + str(readSize) + ')' 
   
  ; load requested field from particle type across all file parts
  pOffset = 0UL
  
  for i=0,nFiles-1 do begin
      fileID   = h5f_open(fileList[i])
      
      ; get number of requested particle type in this file part
      headerID = h5g_open(fileID,"Header")
      nPart = h5a_read(h5a_open_name(headerID,"NumPart_ThisFile"))
      
      if (nPart[partType] eq 0) then continue
      
      ; decide: reading anything from this file?
      readFlag = 0
      
      for j=0,readStrategy.nBlocks-1 do begin
        ; if block starts in this file, or file starts in block (i.e. block ends in file), we have size>0 to read
        if readStrategy.blockStart[j] ge pOffset and $
           readStrategy.blockStart[j] le pOffset+nPart[partType]-1 then readFlag = 1
        if pOffset ge readStrategy.blockStart[j] and $
           pOffset le readStrategy.blockEnd[j] then readFlag = 1
      endfor
      
      ; debug
      if verbose then begin
        if readFlag eq 0 then print,' ['+str(i)+'] readFlag = 0, skipping'
        if readFlag eq 1 then print,' ['+str(i)+'] readFlag = 1, reading'
      endif
      
      pOffset += nPart[partType]
      
      if readFlag eq 0 then continue
      
      ; get field _data
      groupName = 'PartType'+str(partType)
      
      groupID = h5g_open(fileID,groupName)
      dataSetID = h5d_open(groupID,fieldName)
      
      ; load dataset size
      dataSpaceID = h5d_get_space(dataSetID)
      dataSpaceDims = h5s_get_simple_extent_dims(dataSpaceID)
      
      ; for defining HDF5 subsets (start and length of requested data in each dimension)
      nDims = n_elements(dataSpaceDims)
      start  = ulon64arr(nDims)
      length = ulon64arr(nDims)
      
      ; if multiDimSlice requested, select hyperslab
      if multiDimSliceFlag eq 1 then begin
        ; these are all 2D arrays, determine which column (IDL) to read (row in hdf5/C)
        case field of
          'x'   : fN = 0
          'velx': fN = 0
          'cmx' : fN = 0
          'y'   : fN = 1
          'vely': fN = 1
          'cmy' : fN = 1
          'z'   : fN = 2
          'velz': fN = 2
          'cmz' : fN = 2
          
          'tracer_maxtemp'         : fN = 0
          'tracer_maxtemp_time'    : fN = 1
          'tracer_maxtemp_dens'    : fN = 2
          'tracer_maxdens'         : fN = 3
          'tracer_maxdens_time'    : fN = 4
          'tracer_maxmachnum'      : fN = 5
          'tracer_maxent'          : fN = 6
          'tracer_maxent_time'     : fN = 7
          'tracer_laststartime'    : fN = 8
          'tracer_windcounter'     : fN = 9
          'tracer_exchcounter'     : fN = 10
          'tracer_exchdist'        : fN = 11
          'tracer_exchdisterr'     : fN = 12
        endcase
        
        ; multiDimSlice not yet unified with inds/indRange
        if keyword_set(inds) or keyword_set(indRange) then message,'Error: no multidim yet with inds.'
        
        ; start at this column with length equal to the dataset size
        start[0]  = fN
        length[0] = 1
        length[1] = dataSpaceDims[1]
        
        ; select hyperslab and create a memory space to hold the result (otherwise it is full size+sparse)
        h5s_select_hyperslab, dataSpaceID, start, length, /reset
        memSpaceID = h5s_create_simple(length)
        
        ; read the data in the selected hyperslab
        groupData = reform(h5d_read(dataSetID, file_space=dataSpaceID, memory_space=memSpaceID))
        
        totReadSizeLocal = nPart[partType]
      endif else begin
        ; if more than one dimension, want all of the first (e.g. all of xyz in pos)
        if nDims gt 1 then begin
          start[0] = 0
          length[0] = dataSpaceDims[0]
        endif
        
        ; for the last dimension (second if multidim), add simple hyperslabs based on our read strategy
        totReadSizeLocal = 0L
		
	  firstIndexLocal = pOffset - nPart[partType]
	  lastIndexLocal = pOffset - 1
        
        for j=0,readStrategy.nBlocks-1 do begin
          ; skip this block if it has no intersection with this file
          if (readStrategy.blockStart[j] lt firstIndexLocal and $ ; block does not start in this file
              readStrategy.blockEnd[j] lt firstIndexLocal) or $
             (readStrategy.blockStart[j] gt lastIndexLocal and $ ; file does not start in this block
              readStrategy.blockEnd[j] gt lastIndexLocal) then continue
           
          ; determine offset within this file of blockStart and blockEnd
          blockLen     = readStrategy.blockEnd[j] - readStrategy.blockStart[j] + 1
          bOffsetStart = readStrategy.blockStart[j] - firstIndexLocal
          bOffsetEnd   = bOffsetStart + blockLen - 1
          
          ; truncate start and end to the local extent of this file
          if bOffsetStart lt 0 then bOffsetStart = 0
          if bOffsetEnd ge nPart[partType] then bOffsetEnd = nPart[partType] - 1
          
          readSizeBlock = (bOffsetEnd - bOffsetStart + 1)
          
          ; sanity check
          if bOffsetStart ge nPart[partType] then message,'Error: Bad starting local block offset.'
          if bOffsetEnd lt 0 then message,'Error: Bad ending local block offset.'
          
          ; set slab dimensions and add to dataSpaceID
          start[nDims-1] = bOffsetStart
          length[nDims-1] = readSizeBlock
          
          if verbose then print,'  Block ['+str(j)+'] local start = ' + $
            str(bOffsetStart) + ' length = ' + str(readSizeBlock)
          
          h5s_select_hyperslab, dataSpaceID, start, length, reset=(totReadSizeLocal eq 0) ; clear on first
		  
          totReadSizeLocal += readSizeBlock
        endfor
        
        if totReadSizeLocal eq 0 then message,'Error: Zero read size for this file, but not skipped.'
        
        if h5s_get_select_npoints(dataSpaceID) ne totReadSizeLocal*rDims then $
	    message,'Error: HDF5 dataspace plan mismatch with blocking strategy.'
		  
        ; create a memory space to hold the total local result
        length[nDims-1] = totReadSizeLocal ; override last block size with total
        memSpaceID = h5s_create_simple(length)
        
        ; normal read of all data in the dataSet
        if ~h5s_select_valid(dataSpaceID) then message,'Error: Invalid hdf5 selection.'

        groupData = reform(h5d_read(dataSetID, file_space=dataSpaceID, memory_space=memSpaceID))
      endelse

      ; close file
	h5s_close, dataSpaceID
      h5d_close, dataSetID
      h5g_close, groupID
      h5g_close, headerID
      h5f_close, fileID
      
      if ((size(groupData))[0] ne rDims and (size(groupData))[1] ne rDims) then $
        message,'ERROR: Return dimensionality of requested field not expected.'
  
      ; fill return array
      if rDims eq 1 then r[count : (count + totReadSizeLocal - 1)] = groupData
      if rDims gt 1 then r[*,count : (count + totReadSizeLocal - 1)] = groupData
      
      count += totReadSizeLocal
  endfor
  
  if pOffset ne nPartTot[partType] then message,'Error: Read failure.'
  ;if count ne readSize then message,'Error: Read failure.' ; not true if block is truncated in file
  
  ; if we have known indices and an I/O strategy, take the inds subset
  groupData = !NULL
  
  if keyword_set(inds) then begin
    ; make a return array with the size of inds, keep type
    if rDims eq 1 then rr = r[0:n_elements(inds)-1]
    if rDims gt 1 then rr = r[*,0:n_elements(inds)-1]
    
    ; stamp in block by block
    blockOffset = 0LL
    
    for j=0,readStrategy.nBlocks-1 do begin
      
      ; where in the output rr to place the reads from this block (original inds ordering)
      block_output_inds = sort_inds[readStrategy.indexMin[j] : readStrategy.indexMax[j]]
      
      ; where in the full read r to obtain the reads of these elements
      block_read_inds  = sorted_inds[readStrategy.indexMin[j] : readStrategy.indexMax[j]] $
                         - readStrategy.blockStart[j] + blockOffset
                            
      if rDims eq 1 then rr[block_output_inds] = r[ block_read_inds ]
      if rDims gt 1 then rr[*,block_output_inds] = r[ *, block_read_inds]
        
      blockLen = readStrategy.blockEnd[j] - readStrategy.blockStart[j] + 1
      blockOffset += blockLen
    endfor
    
    ; DEBUG (painful, recursively call the snapshotLoad for all particles then directly compare the subset)
    ;aaa = loadSnapshotSubset(sP=sP,partType=PT,field=field)
    ;if rDims eq 1 then aaa = aaa[inds]
    ;if rDims gt 1 then aaa = aaa[*,inds]
    ;if ~array_equal(aaa,rr) then message,'Error: Blocking strategy fail.'
    ; END DEBUG
    
    return,rr
  endif
  
  return,r
end

; createSnapshotCutout(): create a new HDF5 snapshot-format file containing a spatial subregion

pro createSnapshotCutout, sP=sP, fOut=fOut, cenPos=cenPos, boxSize=boxSize, $
                          includeGas=includeGas, includeStars=includeStars, includeDM=includeDM, $
                          excludeSize=excludeSize, convertUtoTemp=convertUtoTemp, verbose=verbose

  ; config
  if (not keyword_set(fOut)   or not keyword_set(sP) or $
      not keyword_set(cenPos) or not keyword_set(boxSize)) then begin
    print,'Error: Missing required inputs.'
    return
  endif
  if (not keyword_set(includeGas) and not keyword_set(includeStars) and not keyword_set(includeDM)) then begin
    print,'Error: No particle type requested.'
    return
  endif
  
  if not keyword_set(includeGas)   then countGas = 0
  if not keyword_set(includeStars) then countStars = 0
  if not keyword_set(includeDM)    then countDM = 0
  if not keyword_set(verbose)      then verbose = 0

  ; load original (structure only)
  ext = string(sP.snap,format='(I3.3)')
  f = sP.simPath + 'snapdir_' + ext + '/snap_' + ext + '.0.hdf5'
  s = h5_parse(f)

  ; modify base and header
  s._NAME    = fOut
  s._FILE    = fOut
  s._COMMENT = "dnelson snapshot cutout"
  s.header.NumFilesPerSnapshot._DATA = 1
  
  ; load gas positions and make selection
  ; ---
  if keyword_set(includeGas) then begin
    gasfield = loadSnapshotSubset(sP=sP,partType='gas',field='pos',verbose=verbose)
        
    wGas = where( abs(gasfield[0,*]-cenPos[0]) le boxSize[0]/2.0 and $
                  abs(gasfield[1,*]-cenPos[1]) le boxSize[1]/2.0 and $
                  abs(gasfield[2,*]-cenPos[2]) le boxSize[2]/2.0, countGas)
                  
    print,'Requested cutout contains ['+str(countGas)+'] gas particles.'

    gasfield = gasfield[*,wGas]
    
    if keyword_set(excludeSize) then begin
      dists = periodicDists(cenPos,gasfield,sP=sP)
      wKeep = where( dists ge excludeSize, countGas )
      wGas = wGas[wKeep]
      print,'After central exclusion have ['+str(countGas)+'] remaining.'
    endif
    
    ; gas - coordinates
    s.PARTTYPE0.COORDINATES._DIMENSIONS    = [3,countGas]
    s.PARTTYPE0.COORDINATES._NELEMENTS     = countGas*3
    s1 = mod_struct(s.PARTTYPE0.COORDINATES,'_DATA',gasfield) ;change _DATA size
    s2 = mod_struct(s.PARTTYPE0,'COORDINATES',s1) ;update PARTTYPE0 with child  
  
    ; gas - CM
    ;if tag_exist(s.parttype0,'center_of_mass') then begin
    ;  gasfield = loadSnapshotSubset(sP=sP,partType='gas',field='cm',verbose=verbose)
    ;  gasfield = gasfield[*,wGas]
    ; 
    ;  s.PARTTYPE0.CENTER_OF_MASS._DIMENSIONS    = [3,countGas]
    ;  s.PARTTYPE0.CENTER_OF_MASS._NELEMENTS     = countGas*3
    ;  s1 = mod_struct(s.PARTTYPE0.CENTER_OF_MASS,'_DATA',gasfield)
    ;  s2 = mod_struct(s2,'CENTER_OF_MASS',s1)
    ;endif
    
    ; gas - coolingrate
    ;if tag_exist(s.parttype0,'coolingrate') then begin
    ;  gasfield = loadSnapshotSubset(sP=sP,partType='gas',field='coolingrate',verbose=verbose)
    ;  gasfield = gasfield[wGas]
    ;  
    ;  s.PARTTYPE0.COOLINGRATE._DIMENSIONS    = [countGas]
    ;  s.PARTTYPE0.COOLINGRATE._NELEMENTS     = countGas
    ;  s1 = mod_struct(s.PARTTYPE0.COOLINGRATE,'_DATA',gasfield)
    ;  s2 = mod_struct(s2,'COOLINGRATE',s1)
    ;endif
    
    ; gas - density
    if tag_exist(s.parttype0,'density') then begin
      gasfield = loadSnapshotSubset(sP=sP,partType='gas',field='density',verbose=verbose)
      gasfield = gasfield[wGas]
      
      s.PARTTYPE0.DENSITY._DIMENSIONS    = [countGas]
      s.PARTTYPE0.DENSITY._NELEMENTS     = countGas
      s1 = mod_struct(s.PARTTYPE0.DENSITY,'_DATA',gasfield)
      s2 = mod_struct(s2,'DENSITY',s1)
    endif
    
    ; gas - electronabundance
    ;if tag_exist(s.parttype0,'electronabundance') then begin
    ;  gasfield = loadSnapshotSubset(sP=sP,partType='gas',field='nelec',verbose=verbose)
    ;  gasfield = gasfield[wGas]
    ;  
    ;  s.PARTTYPE0.ELECTRONABUNDANCE._DIMENSIONS    = [countGas]
    ;  s.PARTTYPE0.ELECTRONABUNDANCE._NELEMENTS     = countGas
    ;  s1 = mod_struct(s.PARTTYPE0.ELECTRONABUNDANCE,'_DATA',gasfield)
    ;  s2 = mod_struct(s2,'ELECTRONABUNDANCE',s1)
    ;endif
    
    ; gas - internal energy
    if tag_exist(s.parttype0,'internalenergy') then begin
      if keyword_set(convertUtoTemp) and tag_exist(s.parttype0,'electronabundance') then begin
        utherm = loadSnapshotSubset(sP=sP,partType='gas',field='u',verbose=verbose)
        utherm = utherm[wGas]
        gasfield = convertUtoTemp(utherm,gasfield) ; Kelvin
      endif else begin
        gasfield = loadSnapshotSubset(sP=sP,partType='gas',field='u',verbose=verbose)
        gasfield = gasfield[wGas]
      endelse
      
      s.PARTTYPE0.INTERNALENERGY._DIMENSIONS    = [countGas]
      s.PARTTYPE0.INTERNALENERGY._NELEMENTS     = countGas
      s1 = mod_struct(s.PARTTYPE0.INTERNALENERGY,'_DATA',gasfield)
      s2 = mod_struct(s2,'INTERNALENERGY',s1)
    endif
    
    ; gas - machnumber
    ;if tag_exist(s.parttype0,'machnumber') then begin
    ;  gasfield = loadSnapshotSubset(sP=sP,partType='gas',field='machnum',verbose=verbose)
    ;  gasfield = gasfield[wGas]
    ;  
    ;  s.PARTTYPE0.MACHNUMBER._DIMENSIONS    = [countGas]
    ;  s.PARTTYPE0.MACHNUMBER._NELEMENTS     = countGas
    ;  s1 = mod_struct(s.PARTTYPE0.MACHNUMBER,'_DATA',gasfield)
    ;  s2 = mod_struct(s2,'MACHNUMBER',s1)
    ;endif
    
    ; gas - massses
    if tag_exist(s.parttype0,'masses') then begin
      gasfield = loadSnapshotSubset(sP=sP,partType='gas',field='mass',verbose=verbose)
      gasfield = gasfield[wGas]
      
      s.PARTTYPE0.MASSES._DIMENSIONS    = [countGas]
      s.PARTTYPE0.MASSES._NELEMENTS     = countGas
      s1 = mod_struct(s.PARTTYPE0.MASSES,'_DATA',gasfield)
      s2 = mod_struct(s2,'MASSES',s1)
    endif
    
    ; gas - metallicity
    if tag_exist(s.parttype0,'metallicity') then begin
      gasfield = loadSnapshotSubset(sP=sP,partType='gas',field='metallicity',verbose=verbose)
      gasfield = gasfield[wGas]
      
      s.PARTTYPE0.METALLICITY._DIMENSIONS    = [countGas]
      s.PARTTYPE0.METALLICITY._NELEMENTS     = countGas
      s1 = mod_struct(s.PARTTYPE0.METALLICITY,'_DATA',gasfield)
      s2 = mod_struct(s2,'METALLICITY',s1)
    endif
    
    ; gas - nh
    ;if tag_exist(s.parttype0,'neutralhydrogenabundance') then begin
    ;  gasfield = loadSnapshotSubset(sP=sP,partType='gas',field='nh',verbose=verbose)
    ;  gasfield = gasfield[wGas]
    ;  
    ;  s.PARTTYPE0.NEUTRALHYDROGENABUNDANCE._DIMENSIONS    = [countGas]
    ;  s.PARTTYPE0.NEUTRALHYDROGENABUNDANCE._NELEMENTS     = countGas
    ;  s1 = mod_struct(s.PARTTYPE0.NEUTRALHYDROGENABUNDANCE,'_DATA',gasfield)
    ;  s2 = mod_struct(s2,'NEUTRALHYDROGENABUNDANCE',s1)
    ;endif
    
    ; gas - num tracers per cell
    ;if tag_exist(s.parttype0,'numtracers') then begin
    ;  gasfield = loadSnapshotSubset(sP=sP,partType='gas',field='numtracers',verbose=verbose)
    ;  gasfield = gasfield[wGas]
    ;  
    ;  s.PARTTYPE0.NUMTRACERS._DIMENSIONS    = [countGas]
    ;  s.PARTTYPE0.NUMTRACERS._NELEMENTS     = countGas
    ;  s1 = mod_struct(s.PARTTYPE0.NUMTRACERS,'_DATA',gasfield)
    ;  s2 = mod_struct(s2,'NUMTRACERS',s1)
    ;endif
    
    ; gas - particleIDs
    if tag_exist(s.parttype0,'particleids') then begin
      gasfield = loadSnapshotSubset(sP=sP,partType='gas',field='ids',verbose=verbose)
      gasfield = gasfield[wGas]
      
      s.PARTTYPE0.PARTICLEIDS._DIMENSIONS    = [countGas]
      s.PARTTYPE0.PARTICLEIDS._NELEMENTS     = countGas
      s1 = mod_struct(s.PARTTYPE0.PARTICLEIDS,'_DATA',gasfield)
      s2 = mod_struct(s2,'PARTICLEIDS',s1)
    endif
    
    ; gas - potential
    ;if tag_exist(s.parttype0,'potential') then begin
    ;  gasfield = loadSnapshotSubset(sP=sP,partType='gas',field='potential',verbose=verbose)
    ;  gasfield = gasfield[wGas]
    ;  
    ;  s.PARTTYPE0.POTENTIAL._DIMENSIONS    = [countGas]
    ;  s.PARTTYPE0.POTENTIAL._NELEMENTS     = countGas
    ;  s1 = mod_struct(s.PARTTYPE0.POTENTIAL,'_DATA',gasfield)
    ;  s2 = mod_struct(s2,'POTENTIAL',s1)
    ;endif
    
    ; gas - hsml
    if tag_exist(s.parttype0,'smoothinglength') then begin
      gasfield = loadSnapshotSubset(sP=sP,partType='gas',field='hsml',verbose=verbose)
      gasfield = gasfield[wGas]
      
      s.PARTTYPE0.SMOOTHINGLENGTH._DIMENSIONS    = [countGas]
      s.PARTTYPE0.SMOOTHINGLENGTH._NELEMENTS     = countGas
      s1 = mod_struct(s.PARTTYPE0.SMOOTHINGLENGTH,'_DATA',gasfield)
      s2 = mod_struct(s2,'SMOOTHINGLENGTH',s1)
    endif
    
    ; gas - sfr
    ;if tag_exist(s.parttype0,'starformationrate') then begin
    ;  gasfield = loadSnapshotSubset(sP=sP,partType='gas',field='sfr',verbose=verbose)
    ;  gasfield = gasfield[wGas]
    ;  
    ;  s.PARTTYPE0.STARFORMATIONRATE._DIMENSIONS    = [countGas]
    ;  s.PARTTYPE0.STARFORMATIONRATE._NELEMENTS     = countGas
    ;  s1 = mod_struct(s.PARTTYPE0.STARFORMATIONRATE,'_DATA',gasfield)
    ;  s2 = mod_struct(s2,'STARFORMATIONRATE',s1)
    ;endif
    
    ; gas - velocities
    if tag_exist(s.parttype0,'velocities') then begin
      gasfield = loadSnapshotSubset(sP=sP,partType='gas',field='vel',verbose=verbose)
      gasfield = gasfield[*,wGas]
      
      s.PARTTYPE0.VELOCITIES._DIMENSIONS    = [3,countGas]
      s.PARTTYPE0.VELOCITIES._NELEMENTS     = countGas*3
      s1 = mod_struct(s.PARTTYPE0.VELOCITIES,'_DATA',gasfield)
      s2 = mod_struct(s2,'VELOCITIES',s1)
    endif
    
    ; gas - volume
    if tag_exist(s.parttype0,'volume') then begin
      gasfield = loadSnapshotSubset(sP=sP,partType='gas',field='volume',verbose=verbose)
      gasfield = gasfield[wGas]
      
      s.PARTTYPE0.VOLUME._DIMENSIONS    = [countGas]
      s.PARTTYPE0.VOLUME._NELEMENTS     = countGas
      s1 = mod_struct(s.PARTTYPE0.VOLUME,'_DATA',gasfield)
      s2 = mod_struct(s2,'VOLUME',s1)
    endif
    
    ;import new PARTTYPE0 structure
    s = mod_struct(s,'PARTTYPE0',s2)
    
    gasfield = !NULL
    s1 = !NULL
    s2 = !NULL
    
  endif ;gas
  
  ; load DM positions and make selection
  ; ---
  if keyword_set(includeDM) then begin
    dmfield = loadSnapshotSubset(sP=sP,partType='dm',field='pos',verbose=verbose)
  
    wDM = where( abs(dmfield[0,*]-cenPos[0]) le boxSize[0]/2.0 and $
                 abs(dmfield[1,*]-cenPos[1]) le boxSize[1]/2.0 and $
                 abs(dmfield[2,*]-cenPos[2]) le boxSize[2]/2.0, countDM)
                    
    print,'Requested cutout contains ['+str(countDM)+'] dark matter particles.'              
                  
    dmfield = dmfield[*,wDM]
    
    ; modify data - coordinates
    s.PARTTYPE1.COORDINATES._DIMENSIONS    = [3,countDM]
    s.PARTTYPE1.COORDINATES._NELEMENTS     = countDM*3
    s1 = mod_struct(s.PARTTYPE1.COORDINATES,'_DATA',dmfield)
    s2 = mod_struct(s.PARTTYPE1,'COORDINATES',s1)
  
    ; dm - particleIDs
    dmfield = loadSnapshotSubset(sP=sP,partType='dm',field='ids',verbose=verbose)
    dmfield = dmfield[wDM]
    
    s.PARTTYPE1.PARTICLEIDS._DIMENSIONS    = [countDM]
    s.PARTTYPE1.PARTICLEIDS._NELEMENTS     = countDM
    s1 = mod_struct(s.PARTTYPE1.PARTICLEIDS,'_DATA',dmfield)
    s2 = mod_struct(s2,'PARTICLEIDS',s1)
    
    ; dm - potential
    dmfield = loadSnapshotSubset(sP=sP,partType='dm',field='potential',verbose=verbose)
    dmfield = dmfield[wDM]
    
    s.PARTTYPE1.POTENTIAL._DIMENSIONS    = [countDM]
    s.PARTTYPE1.POTENTIAL._NELEMENTS     = countDM
    s1 = mod_struct(s.PARTTYPE1.POTENTIAL,'_DATA',dmfield)
    s2 = mod_struct(s2,'POTENTIAL',s1)
  
    ; dm - velocities
    dmfield = loadSnapshotSubset(sP=sP,partType='dm',field='vel',verbose=verbose)
    dmfield = dmfield[*,wDM]
    
    s.PARTTYPE1.VELOCITIES._DIMENSIONS    = [3,countDM]
    s.PARTTYPE1.VELOCITIES._NELEMENTS     = countDM*3
    s1 = mod_struct(s.PARTTYPE1.VELOCITIES,'_DATA',dmfield)
    s2 = mod_struct(s2,'VELOCITIES',s1)
  
    ;import new PARTTYPE1 structure
    s = mod_struct(s,'PARTTYPE1',s2)
    
    dmfield = !NULL
    s1 = !NULL
    s2 = !NULL
  
  endif ;dm
  
  ; load star positions and make selection
  ; ---
  if keyword_set(includeStars) then begin
    starfield = loadSnapshotSubset(sP=sP,partType='stars',field='pos',verbose=verbose)
  
    wStars = where( abs(starfield[0,*]-cenPos[0]) le boxSize[0]/2.0 and $
                    abs(starfield[1,*]-cenPos[1]) le boxSize[1]/2.0 and $
                    abs(starfield[2,*]-cenPos[2]) le boxSize[2]/2.0, countStars)
                    
    print,'Requested cutout contains ['+str(countStars)+'] star particles.'              
                  
    starfield = starfield[*,wStars]
    
    ; modify data - coordinates
    s.PARTTYPE4.COORDINATES._DIMENSIONS    = [3,countStars]
    s.PARTTYPE4.COORDINATES._NELEMENTS     = countStars*3
    s1 = mod_struct(s.PARTTYPE4.COORDINATES,'_DATA',starfield)
    s2 = mod_struct(s.PARTTYPE4,'COORDINATES',s1)
  
    ; stars - masses
    starfield = loadSnapshotSubset(sP=sP,partType='stars',field='mass',verbose=verbose)
    starfield = starfield[wStars]
    
    s.PARTTYPE4.MASSES._DIMENSIONS    = [countStars]
    s.PARTTYPE4.MASSES._NELEMENTS     = countStars
    s1 = mod_struct(s.PARTTYPE4.MASSES,'_DATA',starfield)
    s2 = mod_struct(s2,'MASSES',s1)
   
    ; stars - metallicity
    starfield = loadSnapshotSubset(sP=sP,partType='stars',field='metallicity',verbose=verbose)
    starfield = starfield[wStars]
    
    s.PARTTYPE4.METALLICITY._DIMENSIONS    = [countStars]
    s.PARTTYPE4.METALLICITY._NELEMENTS     = countStars
    s1 = mod_struct(s.PARTTYPE4.METALLICITY,'_DATA',starfield)
    s2 = mod_struct(s2,'METALLICITY',s1)
  
    ; stars - particleIDs
    starfield = loadSnapshotSubset(sP=sP,partType='stars',field='ids',verbose=verbose)
    starfield = starfield[wStars]
    
    s.PARTTYPE4.PARTICLEIDS._DIMENSIONS    = [countStars]
    s.PARTTYPE4.PARTICLEIDS._NELEMENTS     = countStars
    s1 = mod_struct(s.PARTTYPE4.PARTICLEIDS,'_DATA',starfield)
    s2 = mod_struct(s2,'PARTICLEIDS',s1)
    
    ; stars - potential
    starfield = loadSnapshotSubset(sP=sP,partType='stars',field='potential',verbose=verbose)
    starfield = starfield[wStars]
    
    s.PARTTYPE4.POTENTIAL._DIMENSIONS    = [countStars]
    s.PARTTYPE4.POTENTIAL._NELEMENTS     = countStars
    s1 = mod_struct(s.PARTTYPE4.POTENTIAL,'_DATA',starfield)
    s2 = mod_struct(s2,'POTENTIAL',s1)
  
    ; stars - particleIDs
    starfield = loadSnapshotSubset(sP=sP,partType='stars',field='sftime',verbose=verbose)
    starfield = starfield[wStars]
    
    s.PARTTYPE4.STELLARFORMATIONTIME._DIMENSIONS    = [countStars]
    s.PARTTYPE4.STELLARFORMATIONTIME._NELEMENTS     = countStars
    s1 = mod_struct(s.PARTTYPE4.STELLARFORMATIONTIME,'_DATA',starfield)
    s2 = mod_struct(s2,'STELLARFORMATIONTIME',s1)
  
    ; stars - velocities
    starfield = loadSnapshotSubset(sP=sP,partType='stars',field='vel',verbose=verbose)
    starfield = starfield[*,wStars]
    
    s.PARTTYPE4.VELOCITIES._DIMENSIONS    = [3,countStars]
    s.PARTTYPE4.VELOCITIES._NELEMENTS     = countStars*3
    s1 = mod_struct(s.PARTTYPE4.VELOCITIES,'_DATA',starfield)
    s2 = mod_struct(s2,'VELOCITIES',s1)
  
    ;import new PARTTYPE4 structure
    s = mod_struct(s,'PARTTYPE4',s2)
    
    starfield = !NULL
    s1 = !NULL
    s2 = !NULL  
  
  endif ;stars
  
  ; modify HEADER
  s.HEADER._FILE                      = fOut
  s.HEADER.NUMPART_THISFILE._DATA     = [countGas,countDM,0,0,countStars,0]
  s.HEADER.NUMPART_TOTAL._DATA        = [countGas,countDM,0,0,countStars,0]  
  
  ;delete unused structures
  if ~keyword_set(includeGas) then $
    s = mod_struct(s,'PARTTYPE0',/delete)
  if ~keyword_set(includeDM) then $
    s = mod_struct(s,'PARTTYPE1',/delete)
  if ~keyword_set(includeStars) then $
    s = mod_struct(s,'PARTTYPE4',/delete)
    
  if tag_exist(s,'parttype2') then $ ; tracerVEL
    s = mod_struct(s,'PARTTYPE2',/delete)
  if tag_exist(s,'parttype3') then $ ; tracerMC
    s = mod_struct(s,'PARTTYPE3',/delete)
  
  ; output
  h5_create, fOut, s

end

; checkIDRanges()

pro checkIDRanges, res=res

  ; config
  gadgetPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'

  ; loop over snapshots and load ids
  mList = indgen(314)
  foreach m,mList do begin
  
    gas_ids  = loadSnapshotSubset(gadgetPath,snapNum=m,partType='gas',field='ids')
    dm_ids   = loadSnapshotSubset(gadgetPath,snapNum=m,partType='dm',field='ids')
    star_ids = loadSnapshotSubset(gadgetPath,snapNum=m,partType='stars',field='ids')    
  
    print,"["+str(m)+"] z="+str(snapNumToRedshift(snapNum=m))+" nStars = "+str(n_elements(star_ids)) + $
          " ("+str(min(star_ids))+"-"+str(max(star_ids)) + ") nGas = "+str(n_elements(gas_ids)) + $
          " ("+str(min(gas_ids))+"-"+str(max(gas_ids)) + ") nDM = "+str(n_elements(dm_ids)) + $
          " ("+str(min(dm_ids))+"-"+str(max(dm_ids)) + ") nStarsAndGas = " + $
          str(n_elements(star_ids) + n_elements(gas_ids)) 
   
  endforeach

end
