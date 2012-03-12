; cosmoLoad.pro
; cosmological simulations - loading procedures (snapshots, fof/subhalo group cataloges)
; dnelson jan.2012

; getTypeSortedIDList(): within the group catalog ID list rearrange the IDs for each FOF group to be 
;                        ordered first by type (not by SubNr since Subfind was run) such that the
;                        groupOffsetType indexing can be used

function getTypeSortedIDList, sP=sP, gc=gc

  ; load header
  h  = loadsnapshotHeader(sP=sP)
  
  if (h.nPartTot[2] ne 0 or h.nPartTot[5] ne 0) then stop ; not implemented
  
  ; save/restore
  saveFilename = sP.derivPath + sP.savPrefix + str(sP.res) + '.gcSortedIDs.snap=' + str(sP.snap) + '.sav'
                 
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,sortedIDList
  endif
    
  ; for each type, load IDs, match to group cat IDs
  gas_ids  = loadsnapshotSubset(sP=sP,partType='gas',field='ids')
  match,gas_ids,gc.IDs,gas_ind,gc_ind_gas,count=count_gas
  gas_ids  = !NULL
  
  dm_ids   = loadsnapshotSubset(sP=sP,partType='dm',field='ids')
  match,dm_ids,gc.IDs,dm_ind,gc_ind_dm,count=count_dm
  dm_ids   = !NULL
  
  star_ids = loadsnapshotSubset(sP=sP,partType='star',field='ids')
  match,star_ids,gc.IDs,star_ind,gc_ind_star,count=count_star
  star_ids = !NULL
  
  ; reorder indices into ID arrays to match order found in gc.IDs
  gc_ind_gas  = gc_ind_gas[sort(gc_ind_gas)]
  gc_ind_dm   = gc_ind_dm[sort(gc_ind_dm)]
  gc_ind_star = gc_ind_star[sort(gc_ind_star)]
  
  ; verify counts
  if (count_gas + count_dm + count_star ne total(gc.groupLen)) then stop
  if (count_gas ne total(gc.groupLenType[0,*])) then stop
  if (count_dm ne total(gc.groupLenType[1,*])) then stop
  if (count_star ne total(gc.groupLenType[4,*])) then stop
  
  start_gas  = 0L
  start_dm   = 0L
  start_star = 0L
  offset     = 0L
  
  ; DEBUG:
  mask_gas  = intarr(max(gc.IDs[gc_ind_gas])+1)
  mask_dm   = intarr(max(gc.IDs[gc_ind_dm])+1)
  mask_star = intarr(max(gc.IDs[gc_ind_star])+1)
  
  sortedIDList = lon64arr(gc.nIDsTot)
  
  ; loop over each fof group
  for i=0L,gc.nGroupsTot-1 do begin
    ; gas
    if (gc.groupLenType[0,i] gt 0) then begin
      halo_gas_ids  = gc.IDs[gc_ind_gas[start_gas:start_gas+gc.groupLenType[0,i]-1]]
      if n_elements(halo_gas_ids)  ne gc.groupLenType[0,i] then stop
      
      ; DEBUG: fill mask as check for duplicates
      mask_gas[halo_gas_ids] += 1
      
      ; fill in sorted ID list
      sortedIDList[offset:offset+n_elements(halo_gas_ids)-1] = halo_gas_ids
      offset += n_elements(halo_gas_ids)    
      
      ; increment
      start_gas  += n_elements(halo_gas_ids)
    endif
    
    ; dm
    if (gc.groupLenType[1,i] gt 0) then begin
      halo_dm_ids   = gc.IDs[gc_ind_dm[start_dm:start_dm+gc.groupLenType[1,i]-1]]
      if n_elements(halo_dm_ids)   ne gc.groupLenType[1,i] then stop
      mask_dm[halo_dm_ids] += 1
      sortedIDList[offset:offset+n_elements(halo_dm_ids)-1] = halo_dm_ids
      offset += n_elements(halo_dm_ids)
      start_dm   += n_elements(halo_dm_ids)
    endif
    
    ; stars
    if (gc.groupLenType[4,i] gt 0) then begin
      halo_star_ids = gc.IDs[gc_ind_star[start_star:start_star+gc.groupLenType[4,i]-1]]
      if n_elements(halo_star_ids) ne gc.groupLenType[4,i] then stop
      mask_star[halo_star_ids] += 1
      sortedIDList[offset:offset+n_elements(halo_star_ids)-1] = halo_star_ids
      offset += n_elements(halo_star_ids)
      start_star += n_elements(halo_star_ids)
    endif
  endfor
  
  ; DEBUG: check masks
  w_gas  = where(mask_gas gt 1,count_gas)
  w_dm   = where(mask_dm gt 1,count_dm)
  w_star = where(mask_star gt 1,count_star)
  
  if (count_gas gt 0 or count_dm gt 0 or count_star gt 0) then stop
  
  ; DEBUG: match old and sorted ID lists for consistency
  match,gc.IDs,sortedIDList,ind1,ind2,count=count
  if (count ne n_elements(sortedIDList)) then stop
  
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
  endif else begin
    ; check for '/' on end of fileBase
    lastChar = strmid(fileBase,strlen(fileBase)-1,1)
    if (lastChar ne '/') then fileBase += '/'
    
    ; format snap number
    if (m le 999) then $
      ext = string(m,format='(I3.3)')
    if (m gt 999) then $
      ext = string(m,format='(I4.4)')
      
    f = fileBase + 'groups_' + ext + '/fof_subhalo_tab_' + ext
  endelse

  ; check for single (non-split)
  if file_test(fileBase+'/fof_subhalo_tab_'+ext+'.hdf5') then begin
    f = fileBase + 'fof_subhalo_tab_' + ext + '.hdf5'
  endif else begin
  
    ; check existance and multiple outputs
    if not file_test(f+'.hdf5') then begin
      if not file_test(f+'.0.hdf5') then begin
        ;print, 'ERROR: group catalog [' + str(m) + '] at ' + fileBase + ' does not exist!'
        flag = 1
      endif
    endif
    
    f = f + '.*hdf5'
    
  endelse
  
  ; look for FoF only results (no subhalos)
  if (flag eq 1) then begin
    f = fileBase + 'groups_' + ext + '/fof_tab_' + ext
  
    ; check for single (non-split)
    if file_test(fileBase+'/fof_tab_'+ext+'.hdf5') then begin
      f = fileBase + 'fof_tab_' + ext + '.hdf5'
    endif else begin
    
      ; check existance and multiple outputs
      if not file_test(f+'.hdf5') then begin
        if not file_test(f+'.0.hdf5') then begin
          print, 'ERROR: group catalog [' + str(m) + '] at ' + fileBase + ' does not exist!'
          stop
        endif
      endif
      
    endelse
  endif ;flag
  
  return, file_search(f)

end

; loadGroupCat(): load new HDF5 fof/subfind group catalog for a given snapshot
;                 construct offset tables and return structure
;                     
; readIDs=1 : by default, skip IDs since we operate under the group ordered snapshot assumption, but
;             if this flag is set then read IDs and include them (if they exist)
;             also generate (GrNr,Type) sorted id list

function loadGroupCat, sP=sP, readIDs=readIDs, verbose=verbose

  if not keyword_set(verbose) then verbose = 0
  !except = 0 ;suppress floating point underflow/overflow errors

  fileList = getGroupCatFileName(sP.simPath,snapNum=sP.snap)

  nFiles = n_elements(fileList)
  
  ; load number of split files from header of first part
  hdf5s    = h5_parse(fileList[0]) ;structure only
  NumFiles = hdf5s.Header.NumFiles._DATA
  SubfindExistsFlag = tag_exist(hdf5s,'SubhaloLen')
  
  if (NumFiles ne nFiles) then begin
    print,'ERROR: NumFiles ['+str(NumFiles)+'] differs from number of files found ['+str(nFiles)+'].'
    stop
  endif
  
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
        nIDsTot             : h.nIDsTot                    $
      }
      
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
        SubgroupIDMostBound : ulonarr(h.nSubgroupsTot)    ,$
        ;SubgroupIDMostBound : ulon64arr(h.nSubgroupsTot)  ,$ ; LONGIDS maybe? >4B
        SubgroupGrNr        : lonarr(h.nSubgroupsTot)     ,$
        SubgroupParent      : ulonarr(h.nSubgroupsTot)     $
      }
        sf = create_struct(sf,sfsub) ;concat
      endif
      
      ; ID load requested?
      if keyword_set(readIDs) then begin
        if (h.nIDsTot eq 0) then begin
          print,'Warning: readIDs requested but no IDs in group catalog!'
          stop
        endif
        
        sfids = { IDs:lonarr(h.nIDsTot) }
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
    
    if (SubfindExistsFlag eq 1) then begin
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
    
    skip += h.nGroups
    
    ; fill sf with subhalo data from this part
    if (SubfindExistsFlag eq 1) then begin
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
  
  ; create offset tables (sort to create ID list is: (1) GrNr, (2) SubNr, (3) Type, (4) BindingEnergy)
  for i=1L, h.nGroupsTot-1 do begin
    sf.GroupOffset[i] = sf.GroupOffset[i-1] + sf.GroupLen[i-1]
  endfor
  
  ; given SubNr sorted before Type, GroupOffsetType can ONLY be used on getTypeSortedIDList
  for i=0L, h.nGroupsTot-1 do begin
    typeCumSum = [0,total(sf.groupLenType[0:4,i],/cum,/pres)]
    sf.GroupOffsetType[*,i] = sf.GroupOffset[i] + typeCumSum
  endfor
  
  ; if ID read requested, create typeSortedIDList (and save), add to return structure
  if keyword_set(readIDs) then begin
    sfsorted = { IDsSorted:getTypeSortedIDList(sP=sP,gc=sf) }
    sf = create_struct(sf,sfsorted) ;concat
  endif
  
  ; TODO: SubgroupOffset and SubgroupOffsetType
  
  ; verify accumulated totals with last header totals
  if ((nGroupsTot ne h.nGroupsTot) or (nSubgroupsTot ne h.nSubgroupsTot) or $ 
      (nIDsTot ne h.nIDsTot and keyword_set(readIDs))) then begin
    print,'ERROR: Totals do not add up.'
    stop
  endif
  
  ; if ID read was not requested but IDs exist, stop for now (possibly under the group ordered assumption)
  if (nIDsTot gt 0 and not keyword_set(readIDs)) then begin
    print,'Warning: readIDs not requested, but IDs present in group catalog!'
    ;stop
  endif
  
  !except = 1
  
  return,sf
end

; loadSubhaloGroups(): load (OLD, not HDF5) complete subfind group catalog for a given snapshot
;
; skipIDs=1 : don't load actual group member particle IDs

function loadSubhaloGroups, sP=sP, verbose=verbose, skipIDs=skipIDs

  if not keyword_set(verbose) then verbose = 0

  ; set filename
  ext = string(m,format='(I3.3)')
  fIDs = sP.simPath + 'groups_' + ext + '/subhalo_ids_' + ext
  fTab = sP.simPath + 'groups_' + ext + '/subhalo_tab_' + ext
  
  ; check existance and multiple outputs
  if not file_test(fIDs) then begin
    if (file_test(fIDs+'.0')) then begin
      ; split into multiples, get count
      nSplit_IDs = n_elements(file_search(fIDs+".*"))
      nSplit_tab = n_elements(file_search(fTab+".*"))
    endif else begin
      print, 'ERROR: group_ids file ' + sP.simPath + str(sP.snap) + ' does not exist!'
      return,0
    endelse
  endif
  
  if (nSplit_IDs ne nSplit_tab) then begin
    print, 'ERROR: different number of ids and tab files'
    return,0
  endif
  
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

; getSnapFilelist(): take input path and snapshot number and find snapshot filename

function getSnapFilelist, fileBase, snapNum=m, groupOrdered=groupOrdered, subBox=subBox

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
  if not keyword_set(groupOrdered) then $
  if file_test(fileBase+'snap_'+sbstr+ext+'.hdf5') then $
    return, file_search(fileBase + 'snap_' + sbstr + ext + ".*hdf5")

  ; check for single groupordered
  if file_test(fileBase+'snap-groupordered_'+ext+'.hdf5') then $
    return, file_search(fileBase + 'snap-groupordered_' + ext+".*hdf5")
  
  ; check for multiple groupordered
  if file_test(fileBase+'snapdir_'+ext+'/snap-groupordered_'+ext+'.0.hdf5') then $
    return, file_search(fileBase + 'snapdir_' + ext + '/snap-groupordered_' + ext+".*hdf5")
  
  ; check for exact name or multiple outputs
  if (file_test(f+'.hdf5') or file_test(f+'.0.hdf5')) then $
    return, file_search(f+".*hdf5")

  print,'Error: Failed to find snapshot.'
  stop
end

; loadSnapshotHeader(): load header

function loadSnapshotHeader, sP=sP, verbose=verbose, subBox=subBox

  if not keyword_set(verbose) then verbose = 0

  ; get matching filename (return -1 if not found)
  fileList = getSnapFilelist(sP.simPath,snapNum=sP.snap,subBox=subBox)

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

; loadSnapshotSubset(): for a given snapshot load only one field for one particle type
;                       partType = [0,1,2,4] or ('gas','dm','tracer','stars') (case insensitive)
;                       field    = ['ParticleIDs','coordinates','xyz',...] (case insensitive)

function loadSnapshotSubset, sP=sP, partType=PT, field=field, $
                             verbose=verbose, $
                             doublePrec=doublePrec, groupOrdered=groupOrdered, subBox=subBox

  if not keyword_set(verbose) then verbose = 0

  fileList = getSnapFilelist(sP.simPath,snapNum=sP.snap,groupOrdered=groupOrdered,subBox=subBox)

  nFiles = n_elements(fileList)
  
  ; input config: set partType number if input in string
  partType = strlowcase(string(PT)) ; so we don't change the input
  if (strcmp(partType,'gas')       or strcmp(partType,'hydro'))      then partType = 0
  if (strcmp(partType,'dm')        or strcmp(partType,'darkmatter')) then partType = 1
  if (strcmp(partType,'tracervel') or strcmp(partType,'tracersvel')) then partType = 2
  if (strcmp(partType,'tracermc')  or strcmp(partType,'tracersmc'))  then partType = 3
  if (strcmp(partType,'stars')     or strcmp(partType,'star'))       then partType = 4
  
  ; error checking
  if (strcmp(partType,'tracer') or strcmp(partType,'tracers')) then begin
    print,'ERROR: Please specify which type of tracers!'
    stop
  endif
  if (not isnumeric(partType)) then begin
    print,'ERROR: Bad partType = ' + partType
    stop
  endif
  if (partType lt 0 or partType gt 4) then begin
    print,'ERROR: partType = ' + str(partType) + ' out of bounds!'
    stop
  endif
  
  ; load particle array sizes from header of first part
  fileID   = h5f_open(fileList[0])
  headerID = h5g_open(fileID,"Header")
  nPartTot = h5a_read(h5a_open_name(headerID,"NumPart_Total"))
  nSplits  = h5a_read(h5a_open_name(headerID,"NumFilesPerSnapshot"))
  flagDbl  = h5a_read(h5a_open_name(headerID,"Flag_DoublePrecision"))
  h5g_close, headerID
  h5f_close, fileID
  
  if (nSplits ne nFiles) then begin
    print,'ERROR: NumFilesPerSnapshot ['+str(nSplits)+'] differs from number of split files found ['+$
           str(nFiles)+'].'
    stop
  endif
  
  ; double precision
  if (flagDbl eq 1 and not keyword_set(doublePrec)) then $
    print,'Warning: Snapshot is double precision but only singlePrec load requested.'
  if (flagDbl eq 0 and keyword_set(doublePrec)) then $
    print,'Warning: Snapshot is single precision but doublePrec load requested.'
  
  ; early exit: no particles of requested type
  if (nPartTot[partType] eq 0) then return, 0
  
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
    r = fltarr(rDims,nPartTot[partType])
  endif
  if (field eq 'particleids' or field eq 'ids') then begin
    r = lonarr(nPartTot[partType])
    fieldName = 'ParticleIDs'
  endif
  if (field eq 'potential' or field eq 'phi') then begin
    r = fltarr(nPartTot[partType])
    fieldName = 'Potential'
  endif
  if (field eq 'velocities' or field eq 'vel' or $
      field eq 'velx' or field eq 'vely' or field eq 'velz') then begin
    fieldName = 'Velocities'
    rDims = 3
    r = dblarr(rDims,nPartTot[partType])
  endif
  
  ; gas only
  ; --------
  if (field eq 'center_of_mass' or field eq 'centerofmass' or field eq 'com' or field eq 'cm' or $
      field eq 'cmx' or field eq 'cmy' or field eq 'cmz') then begin
    fieldName = 'Center-of-Mass'
    rDims = 3
    r = fltarr(rDims,nPartTot[partType])
    if (partType ne 0) then begin & print,'Error: CoM is gas only!' & return,0 & endif
  endif
  if (field eq 'coolingrate' or field eq 'coolrate') then begin
    r = fltarr(nPartTot[partType])
    fieldName = 'CoolingRate'
    if (partType ne 0) then begin & print,'Error: CoolingRate is gas only!' & return,0 & endif
  endif
  if (field eq 'density' or field eq 'rho' or field eq 'dens') then begin
    r = fltarr(nPartTot[partType])
    fieldName = 'Density'
    if (partType ne 0) then begin & print,'Error: Density is gas only!' & return,0 & endif
  endif
  if (field eq 'electronabundance' or field eq 'ne' or field eq 'nelec') then begin
    r = fltarr(nPartTot[partType])
    fieldName = 'ElectronAbundance'
    if (partType ne 0) then begin & print,'Error: NE is gas only!' & return,0 & endif
  endif
  if (field eq 'internalenergy' or field eq 'u') then begin
    r = fltarr(nPartTot[partType])
    fieldName = 'InternalEnergy'
    if (partType ne 0) then begin & print,'Error: U is gas only!' & return,0 & endif
  endif
  if (field eq 'machnumber' or field eq 'machnum') then begin
    r = fltarr(nPartTot[partType])
    fieldName = 'MachNumber'
    if (partType ne 0) then begin & print,'Error: MachNum is gas only!' & return,0 & endif
  endif
  if (field eq 'maxfaceangle') then begin
    r = fltarr(nPartTot[partType])
    fieldName = 'MaxFaceAngle'
    if (partType ne 0) then begin & print,'Error: MaxFaceAngle is gas only!' & return,0 & endif
  endif
  if (field eq 'neutralhydrogenabundance' or field eq 'nh') then begin
    r = fltarr(nPartTot[partType])
    fieldName = 'NeutralHydrogenAbundance'
    if (partType ne 0) then begin & print,'Error: NH is gas only!' & return,0 & endif
  endif
  if (field eq 'number_of_faces_of_cell' or field eq 'num_cell_faces' or field eq 'numcellfaces' or $
      field eq 'numfaces') then begin
    r = fltarr(nPartTot[partType])
    fieldName = 'Number of faces of cell'
    if (partType ne 0) then begin & print,'Error: NumFaces is gas only!' & return,0 & endif
  endif
  if (field eq 'pressure' or field eq 'pres') then begin
    r = fltarr(nPartTot[partType])
    if (partType ne 0) then begin & print,'Error: Pressure is gas only!' & return,0 & endif
    fieldName = 'Pressure'
  endif
  if (field eq 'smoothinglength' or field eq 'hsml') then begin
    r = fltarr(nPartTot[partType])
    fieldName = 'SmoothingLength'
    if (partType ne 0) then begin & print,'Error: HSML is gas only!' & return,0 & endif
  endif
  if (field eq 'starformationrate' or field eq 'sfr') then begin
    r = fltarr(nPartTot[partType])
    fieldName = 'StarFormationRate'
    if (partType ne 0) then begin & print,'Error: SFR is gas only!' & return,0 & endif
  endif
  if (field eq 'surface_area' or field eq 'surfarea') then begin
    r = fltarr(nPartTot[partType])
    fieldName = 'Surface Area'
    if (partType ne 0) then begin & print,'Error: SurfArea is gas only!' & return,0 & endif
  endif
  if (field eq 'volume' or field eq 'vol') then begin
    r = fltarr(nPartTot[partType])
    fieldName = 'Volume'
    if (partType ne 0) then begin & print,'Error: Vol is gas only!' & return,0 & endif
  endif
  
  ; gas/stars only
  ; --------------
  if (field eq 'masses' or field eq 'mass') then begin
    r = fltarr(nPartTot[partType])
    fieldName = 'Masses'
    if (partType ne 0 and partType ne 4) then begin & print,'Error: Mass is gas/stars only!' & return,0 & endif
  endif
  if (field eq 'metallicity' or field eq 'metal') then begin
    r = fltarr(nPartTot[partType])
    fieldName = 'Metallicity'
    if (partType ne 0 and partType ne 4) then begin & print,'Error: Z is gas/stars only!' & return,0 & endif
  endif
  if (field eq 'numtr' or field eq 'numtracers') then begin
    r = lonarr(nPartTot[partType])
    fieldName = 'NumTracers'
    if (partType ne 0 and partType ne 4) then begin & print,'Error: NumTracers is gas/stars only!' & return,0 & endif
  endif
  
  ; stars only (TODO: GFM)
  ; ----------
  if (field eq 'stellarformationtime' or field eq 'sftime') then begin
    r = fltarr(nPartTot[partType])
    fieldName = 'StellarFormationTime'
    if (partType ne 4) then begin & print,'Error: SFTime is stars only!' & return,0 & endif
  endif
  
  ; tracers (Monte Carlo)
  ; ---------------------
  if (field eq 'parentid' or field eq 'parentids') then begin
    r = lonarr(nPartTot[partType])
    fieldName = 'ParentID'
    if (partType ne 3) then begin & print,'Error: ParentID is tracerMC only!' & return,0 & endif
  endif
  if (field eq 'properties' or field eq 'quants' or field eq 'quantities' or $
      field eq 'tracer_maxtemp' or field eq 'tracer_maxtemp_time' or field eq 'tracer_maxdens' or $
      field eq 'tracer_maxmachnum' or field eq 'tracer_maxentropy') then begin
    fieldName = 'FluidQuantities'
    rDims = 5 ; WARNING: must match to setup in Arepo run
    r = fltarr(rDims,nPartTot[partType])
    if (partType ne 3) then begin & print,'Error: Fluid quantities are tracerMC/Vel only!' & return,0 & endif
  endif
  if (field eq 'tracerid' or field eq 'tracerids') then begin
    r = lonarr(nPartTot[partType])
    fieldName = 'TracerID'
    if (partType ne 3) then begin & print,'Error: TracerID is tracerMC only!' & return,0 & endif
  endif
  
  if (fieldName eq '') then begin
    print,'ERROR: Requested field -- ' + strlowcase(field) + ' -- not recognized!'
    return,0
  endif
  
  if (verbose) then $
    print,'Loading "' + str(fieldName) + '" for partType=' + str(partType) + ' from snapshot (' + $
          str(m) + ') in [' + str(nFiles) + '] files. (nGas=' + $
          str(nPartTot[0]) + ' nDM=' + str(nPartTot[1]) + ' nStars=' + str(nPartTot[4]) + ')' 
   
  ; double precision requested?
  if (size(r,/tname) eq 'FLOAT' and keyword_set(doublePrec)) then $
    r = double(r)
   
  ; load requested field from particle type across all file parts
  for i=0,nFiles-1 do begin
      fileID   = h5f_open(fileList[i])
      
      ; get number of requested particle type in this file part
      headerID = h5g_open(fileID,"Header")
      nPart = h5a_read(h5a_open_name(headerID,"NumPart_ThisFile"))
      
      ; get field _data
      if (nPart[partType] eq 0) then continue
      
      groupName = 'PartType'+str(partType)
      
      groupID = h5g_open(fileID,groupName)
      groupData = h5d_read(h5d_open(groupID,fieldName))

      ; close file
      h5g_close, groupID
      h5g_close, headerID
      h5f_close, fileID
      
      if ((size(groupData))[0] ne rDims and (size(groupData))[1] ne rDims) then begin
        print,'ERROR: Return dimensionality of requested field not expected.'
        return,0
      endif
  
      ; fill return array
      if (rDims eq 1) then $
        r[count : (count + nPart[partType] - 1)] = groupData
      if (rDims gt 1) then $
        r[*,count : (count + nPart[partType] - 1)] = groupData
      
      count += nPart[partType]
  
  endfor
  
  ; if requesting a slice of a multidim array, do it now
  if (field eq 'x' or field eq 'y' or field eq 'z' or $
      field eq 'velx' or field eq 'vely' or field eq 'velz' or $
      field eq 'cmx' or field eq 'cmy' or field eq 'cmz' or $
      field eq 'tracer_maxtemp' or field eq 'tracer_maxtemp_time' or field eq 'tracer_maxdens' or $
      field eq 'tracer_maxmachnum' or field eq 'tracer_maxentropy') then begin
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
      
      'tracer_maxtemp'      : fN = 0
      'tracer_maxtemp_time' : fN = 1
      'tracer_maxdens'      : fN = 2
      'tracer_maxmachnum'   : fN = 3
      'tracer_maxentropy'   : fN = 4
    endcase
    r = reform(r[fN,*])
  endif
  
  return,r
end

; createSnapshotCutout(): create a new HDF5 snapshot-format file containing a spatial subregion

pro createSnapshotCutout, snapPath=snapPath, snapNum=snapNum, fOut=fOut, cenPos=cenPos, boxSize=boxSize, $
                          includeGas=includeGas, includeStars=includeStars, includeDM=includeDM, $
                          verbose=verbose

  ; config
  if (not keyword_set(fOut)   or not keyword_set(snapPath) or not keyword_set(snapNum) or $
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
  ext = string(snapNum,format='(I3.3)')
  f = snapPath + 'snapdir_' + ext + '/snap_' + ext + '.0.hdf5'
  s = h5_parse(f)
  
  ; modify base and header
  s._NAME    = fOut
  s._FILE    = fOut
  s._COMMENT = "dnelson snapshot cutout"
  s.header.NumFilesPerSnapshot._DATA = 1
  
  ; load gas positions and make selection
  ; ---
  if keyword_set(includeGas) then begin
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='pos',verbose=verbose)
    
    wGas = where( abs(gasfield[0,*]-cenPos[0]) le boxSize[0]/2.0 and $
                  abs(gasfield[1,*]-cenPos[1]) le boxSize[1]/2.0 and $
                  abs(gasfield[2,*]-cenPos[2]) le boxSize[2]/2.0, countGas)
                  
    print,'Requested cutout contains ['+str(countGas)+'] gas particles.'              
                  
    gasfield = gasfield[*,wGas]
    
    ; modify data
    s.PARTTYPE0.COORDINATES._DIMENSIONS    = [3,countGas]
    s.PARTTYPE0.COORDINATES._NELEMENTS     = countGas*3
    s1 = mod_struct(s.PARTTYPE0.COORDINATES,'_DATA',gasfield) ;change _DATA size
    s2 = mod_struct(s.PARTTYPE0,'COORDINATES',s1) ;update PARTTYPE0 with child  
  
    ; gas - CM
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='cm',verbose=verbose)
    gasfield = gasfield[*,wGas]
   
    s.PARTTYPE0.CENTER_OF_MASS._DIMENSIONS    = [3,countGas]
    s.PARTTYPE0.CENTER_OF_MASS._NELEMENTS     = countGas*3
    s1 = mod_struct(s.PARTTYPE0.CENTER_OF_MASS,'_DATA',gasfield)
    s2 = mod_struct(s2,'CENTER_OF_MASS',s1)
    
    ; gas - coolingrate
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='coolingrate',verbose=verbose)
    gasfield = gasfield[wGas]
    
    s.PARTTYPE0.COOLINGRATE._DIMENSIONS    = [countGas]
    s.PARTTYPE0.COOLINGRATE._NELEMENTS     = countGas
    s1 = mod_struct(s.PARTTYPE0.COOLINGRATE,'_DATA',gasfield)
    s2 = mod_struct(s2,'COOLINGRATE',s1)
    
    ; gas - density
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='density',verbose=verbose)
    gasfield = gasfield[wGas]
    
    s.PARTTYPE0.DENSITY._DIMENSIONS    = [countGas]
    s.PARTTYPE0.DENSITY._NELEMENTS     = countGas
    s1 = mod_struct(s.PARTTYPE0.DENSITY,'_DATA',gasfield)
    s2 = mod_struct(s2,'DENSITY',s1)
    
    ; gas - electronabundance
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='nelec',verbose=verbose)
    gasfield = gasfield[wGas]
    
    s.PARTTYPE0.ELECTRONABUNDANCE._DIMENSIONS    = [countGas]
    s.PARTTYPE0.ELECTRONABUNDANCE._NELEMENTS     = countGas
    s1 = mod_struct(s.PARTTYPE0.ELECTRONABUNDANCE,'_DATA',gasfield)
    s2 = mod_struct(s2,'ELECTRONABUNDANCE',s1)
    
    ; gas - internal energy
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='u',verbose=verbose)
    gasfield = gasfield[wGas]
    
    s.PARTTYPE0.INTERNALENERGY._DIMENSIONS    = [countGas]
    s.PARTTYPE0.INTERNALENERGY._NELEMENTS     = countGas
    s1 = mod_struct(s.PARTTYPE0.INTERNALENERGY,'_DATA',gasfield)
    s2 = mod_struct(s2,'INTERNALENERGY',s1)
    
    ; gas - machnumber
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='machnum',verbose=verbose)
    gasfield = gasfield[wGas]
    
    s.PARTTYPE0.MACHNUMBER._DIMENSIONS    = [countGas]
    s.PARTTYPE0.MACHNUMBER._NELEMENTS     = countGas
    s1 = mod_struct(s.PARTTYPE0.MACHNUMBER,'_DATA',gasfield)
    s2 = mod_struct(s2,'MACHNUMBER',s1)
    
    ; gas - massses
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='mass',verbose=verbose)
    gasfield = gasfield[wGas]
    
    s.PARTTYPE0.MASSES._DIMENSIONS    = [countGas]
    s.PARTTYPE0.MASSES._NELEMENTS     = countGas
    s1 = mod_struct(s.PARTTYPE0.MASSES,'_DATA',gasfield)
    s2 = mod_struct(s2,'MASSES',s1)
    
    ; gas - maxfaceangle
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='maxfaceangle',verbose=verbose)
    gasfield = gasfield[wGas]
    
    s.PARTTYPE0.MAXFACEANGLE._DIMENSIONS    = [countGas]
    s.PARTTYPE0.MAXFACEANGLE._NELEMENTS     = countGas
    s1 = mod_struct(s.PARTTYPE0.MAXFACEANGLE,'_DATA',gasfield)
    s2 = mod_struct(s2,'MAXFACEANGLE',s1)
    
    ; gas - metallicity
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='metallicity',verbose=verbose)
    gasfield = gasfield[wGas]
    
    s.PARTTYPE0.METALLICITY._DIMENSIONS    = [countGas]
    s.PARTTYPE0.METALLICITY._NELEMENTS     = countGas
    s1 = mod_struct(s.PARTTYPE0.METALLICITY,'_DATA',gasfield)
    s2 = mod_struct(s2,'METALLICITY',s1)
  
    ; gas - nh
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='nh',verbose=verbose)
    gasfield = gasfield[wGas]
    
    s.PARTTYPE0.NEUTRALHYDROGENABUNDANCE._DIMENSIONS    = [countGas]
    s.PARTTYPE0.NEUTRALHYDROGENABUNDANCE._NELEMENTS     = countGas
    s1 = mod_struct(s.PARTTYPE0.NEUTRALHYDROGENABUNDANCE,'_DATA',gasfield)
    s2 = mod_struct(s2,'NEUTRALHYDROGENABUNDANCE',s1)
    
    ; gas - num faces per cell
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='numfaces',verbose=verbose)
    gasfield = gasfield[wGas]
    
    s.PARTTYPE0.NUMBER_OF_FACES_OF_CELL._DIMENSIONS    = [countGas]
    s.PARTTYPE0.NUMBER_OF_FACES_OF_CELL._NELEMENTS     = countGas
    s1 = mod_struct(s.PARTTYPE0.NUMBER_OF_FACES_OF_CELL,'_DATA',gasfield)
    s2 = mod_struct(s2,'NUMBER_OF_FACES_OF_CELL',s1)
    
    ; gas - particleIDs
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='ids',verbose=verbose)
    gasfield = gasfield[wGas]
    
    s.PARTTYPE0.PARTICLEIDS._DIMENSIONS    = [countGas]
    s.PARTTYPE0.PARTICLEIDS._NELEMENTS     = countGas
    s1 = mod_struct(s.PARTTYPE0.PARTICLEIDS,'_DATA',gasfield)
    s2 = mod_struct(s2,'PARTICLEIDS',s1)
    
    ; gas - potential
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='potential',verbose=verbose)
    gasfield = gasfield[wGas]
    
    s.PARTTYPE0.POTENTIAL._DIMENSIONS    = [countGas]
    s.PARTTYPE0.POTENTIAL._NELEMENTS     = countGas
    s1 = mod_struct(s.PARTTYPE0.POTENTIAL,'_DATA',gasfield)
    s2 = mod_struct(s2,'POTENTIAL',s1)
    
    ; gas - hsml
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='hsml',verbose=verbose)
    gasfield = gasfield[wGas]
    
    s.PARTTYPE0.SMOOTHINGLENGTH._DIMENSIONS    = [countGas]
    s.PARTTYPE0.SMOOTHINGLENGTH._NELEMENTS     = countGas
    s1 = mod_struct(s.PARTTYPE0.SMOOTHINGLENGTH,'_DATA',gasfield)
    s2 = mod_struct(s2,'SMOOTHINGLENGTH',s1)
    
    ; gas - sfr
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='sfr',verbose=verbose)
    gasfield = gasfield[wGas]
    
    s.PARTTYPE0.STARFORMATIONRATE._DIMENSIONS    = [countGas]
    s.PARTTYPE0.STARFORMATIONRATE._NELEMENTS     = countGas
    s1 = mod_struct(s.PARTTYPE0.STARFORMATIONRATE,'_DATA',gasfield)
    s2 = mod_struct(s2,'STARFORMATIONRATE',s1)
    
    ; gas - surface area
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='surfarea',verbose=verbose)
    gasfield = gasfield[wGas]
    
    s.PARTTYPE0.SURFACE_AREA._DIMENSIONS    = [countGas]
    s.PARTTYPE0.SURFACE_AREA._NELEMENTS     = countGas
    s1 = mod_struct(s.PARTTYPE0.SURFACE_AREA,'_DATA',gasfield)
    s2 = mod_struct(s2,'SURFACE_AREA',s1)
    
    ; gas - velocities
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='vel',verbose=verbose)
    gasfield = gasfield[*,wGas]
    
    s.PARTTYPE0.VELOCITIES._DIMENSIONS    = [3,countGas]
    s.PARTTYPE0.VELOCITIES._NELEMENTS     = countGas*3
    s1 = mod_struct(s.PARTTYPE0.VELOCITIES,'_DATA',gasfield)
    s2 = mod_struct(s2,'VELOCITIES',s1)
    
    ; gas - volume
    gasfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='gas',field='volume',verbose=verbose)
    gasfield = gasfield[wGas]
    
    s.PARTTYPE0.VOLUME._DIMENSIONS    = [countGas]
    s.PARTTYPE0.VOLUME._NELEMENTS     = countGas
    s1 = mod_struct(s.PARTTYPE0.VOLUME,'_DATA',gasfield)
    s2 = mod_struct(s2,'VOLUME',s1)
    
    ;import new PARTTYPE0 structure
    s = mod_struct(s,'PARTTYPE0',s2)
    
    gasfield = !NULL
    s1 = !NULL
    s2 = !NULL
    
  endif ;gas
  
  ; load DM positions and make selection
  ; ---
  if keyword_set(includeDM) then begin
    dmfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='dm',field='pos',verbose=verbose)
  
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
    dmfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='dm',field='ids',verbose=verbose)
    dmfield = dmfield[wDM]
    
    s.PARTTYPE1.PARTICLEIDS._DIMENSIONS    = [countDM]
    s.PARTTYPE1.PARTICLEIDS._NELEMENTS     = countDM
    s1 = mod_struct(s.PARTTYPE1.PARTICLEIDS,'_DATA',dmfield)
    s2 = mod_struct(s2,'PARTICLEIDS',s1)
    
    ; dm - potential
    dmfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='dm',field='potential',verbose=verbose)
    dmfield = dmfield[wDM]
    
    s.PARTTYPE1.POTENTIAL._DIMENSIONS    = [countDM]
    s.PARTTYPE1.POTENTIAL._NELEMENTS     = countDM
    s1 = mod_struct(s.PARTTYPE1.POTENTIAL,'_DATA',dmfield)
    s2 = mod_struct(s2,'POTENTIAL',s1)
  
    ; dm - velocities
    dmfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='dm',field='vel',verbose=verbose)
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
    starfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='stars',field='pos',verbose=verbose)
  
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
    starfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='stars',field='mass',verbose=verbose)
    starfield = starfield[wStars]
    
    s.PARTTYPE4.MASSES._DIMENSIONS    = [countStars]
    s.PARTTYPE4.MASSES._NELEMENTS     = countStars
    s1 = mod_struct(s.PARTTYPE4.MASSES,'_DATA',starfield)
    s2 = mod_struct(s2,'MASSES',s1)
   
    ; stars - metallicity
    starfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='stars',field='metallicity',verbose=verbose)
    starfield = starfield[wStars]
    
    s.PARTTYPE4.METALLICITY._DIMENSIONS    = [countStars]
    s.PARTTYPE4.METALLICITY._NELEMENTS     = countStars
    s1 = mod_struct(s.PARTTYPE4.METALLICITY,'_DATA',starfield)
    s2 = mod_struct(s2,'METALLICITY',s1)
  
    ; stars - particleIDs
    starfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='stars',field='ids',verbose=verbose)
    starfield = starfield[wStars]
    
    s.PARTTYPE4.PARTICLEIDS._DIMENSIONS    = [countStars]
    s.PARTTYPE4.PARTICLEIDS._NELEMENTS     = countStars
    s1 = mod_struct(s.PARTTYPE4.PARTICLEIDS,'_DATA',starfield)
    s2 = mod_struct(s2,'PARTICLEIDS',s1)
    
    ; stars - potential
    starfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='stars',field='potential',verbose=verbose)
    starfield = starfield[wStars]
    
    s.PARTTYPE4.POTENTIAL._DIMENSIONS    = [countStars]
    s.PARTTYPE4.POTENTIAL._NELEMENTS     = countStars
    s1 = mod_struct(s.PARTTYPE4.POTENTIAL,'_DATA',starfield)
    s2 = mod_struct(s2,'POTENTIAL',s1)
  
    ; stars - particleIDs
    starfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='stars',field='sftime',verbose=verbose)
    starfield = starfield[wStars]
    
    s.PARTTYPE4.STELLARFORMATIONTIME._DIMENSIONS    = [countStars]
    s.PARTTYPE4.STELLARFORMATIONTIME._NELEMENTS     = countStars
    s1 = mod_struct(s.PARTTYPE4.STELLARFORMATIONTIME,'_DATA',starfield)
    s2 = mod_struct(s2,'STELLARFORMATIONTIME',s1)
  
    ; stars - velocities
    starfield = loadSnapshotSubset(snapPath,snapNum=snapNum,partType='stars',field='vel',verbose=verbose)
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
  if not keyword_set(includeGas) then $
    s = mod_struct(s,'PARTTYPE0',/delete)
  if not keyword_set(includeDM) then $
    s = mod_struct(s,'PARTTYPE1',/delete)
  if not keyword_set(includeStars) then $
    s = mod_struct(s,'PARTTYPE4',/delete)
  
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