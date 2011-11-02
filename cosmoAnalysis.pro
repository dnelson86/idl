; cosmoAnalysis.pro
; gas accretion project - main
; dnelson nov.2011

@helper
@cosmoUtil
@cosmoLoad

; galaxyCat(): if snap not specified, create and save complete galaxy catalog from the group catalog by 
;              imposing additional cut in the (rho,temp) plane (same as that used by Torrey+ 2011)
;              if snap is specified, create only for one snapshot number or return previously saved
;              results for that snapshot

function galaxyCat, res=res, run=run, snap=snap

  if not keyword_set(res) or not keyword_set(run) then begin
     print,'Error: galaxyCat() arguments not specified.'
     return,0
  endif

  ; config
  workingPath = '/n/home07/dnelson/coldflows/galaxycat/'
    
  galcut_T   = 6.0
  galcut_rho = 0.25
    
  snapRange = [50,314]    
    
  if (run eq 'gadget') then begin
    simPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'
    savPrefix = 'G'
  endif
  
  ; if snap specified, run only one snapshot (and/or just return previous results)
  if keyword_set(snap) then begin
    saveFilename = workingPath + 'galcat.' + savPrefix + str(res) + '.' + str(snap) + '.sav'
    
    ; results exist, return
    if (file_test(saveFilename)) then begin
      restore,saveFilename
      r = {galaxyOff:galaxyOff,galaxyLen:galaxyLen,galaxyIDs:galaxyIDs}
      return,r
    endif
    
    ; need to compute, set restricted range of snapshots to process
    snapRange = [snap,snap]
  endif
  
  for m=snapRange[0],snapRange[1],1 do begin
    ; skip if previous results exist
    saveFilename = workingPath + 'galcat.' + savPrefix + str(res) + '.' + str(m) + '.sav'
    
    if (file_test(saveFilename)) then begin
      print,'Skipping: '+strmid(saveFilename,strlen(workingPath))
      continue
    endif
    
    ; load ids of particles in all subfind groups
    sg = loadSubhaloGroups(simPath,m)
    sgPIDs = sgPIDList(sg=sg,/all)
    
    ; load gas ids and match to catalog
    ids  = loadSnapshotSubset(simPath,snapNum=m,partType='gas',field='ids')
    match,sgPIDs,ids,sg_ind,ids_ind,count=countMatch,/sort
   
    ids = ids[ids_ind]
    sgPIDs = !NULL

    ; load u,nelec and calculate temp of gas
    u = loadSnapshotSubset(simPath,snapNum=m,partType='gas',field='u')
    u = u[ids_ind]
    nelec = loadSnapshotSubset(simPath,snapNum=m,partType='gas',field='nelec')
    nelec = nelec[ids_ind]
    
    temp = convertUtoTemp(u,nelec)
    
    u     = !NULL
    nelec = !NULL
    
    ; load rho of gas and make galaxy (rho,temp) plane cut
    dens = loadSnapshotSubset(simPath,snapNum=m,partType='gas',field='density')
    dens = dens[ids_ind]
    
    w = where(alog10(temp) - galcut_rho * alog10(dens) lt galcut_T,countCut)

    if (countCut eq 0) then begin
      print,'Warning: Empty galaxy cut. Skipping: ' + strmid(saveFilename,strlen(workingPath))
      continue
    endif
    
    temp = !NULL
    dens = !NULL
    
    ; construct galaxy catalog
    ids = ids[w]
    
    galaxyLen = ulonarr(sg.nSubgroupsTot)
    galaxyOff = ulonarr(sg.nSubgroupsTot)
    galaxyIDs = ulonarr(n_elements(ids))
    
    nextOff = 0UL
    
    ; match indices between galaxy ids and group ids
    match,ids,sg.subgroupIDs,ids_ind,sgIDs_ind,count=countID,/sort
    
    for sgID=0,sg.nSubgroupsTot-1 do begin
      ; select ids in group
      groupStart = sg.subGroupOffset[sgID]
      groupEnd   = groupStart + sg.subGroupLen[sgID]
      
      w = where(sgIDs_ind ge groupStart and sgIDs_ind lt groupEnd,count)
      
      ; save in similar 1D offset format
      if (count gt 0) then begin
        galaxyLen[sgID] = count
        galaxyOff[sgID] = nextOff
        nextOff = galaxyOff[sgID] + galaxyLen[sgID]
  
        galaxyIDs[galaxyOff[sgID]:nextOff-1] = ids[ids_ind[w]]
      endif      
      
    end
    
    ; debug: make sure all galaxy gas particles were found in the group catalog
    ;match,ids,galaxyIDs,ind1,ind2,count=countCheck
    ;if (countCheck ne n_elements(ids)) then begin
    ;  print,'Uhoh, check.',countCheck,n_elements(ids)
    ;  nuniq = n_elements(uniq(ids[sort(ids)]))
    ;  stop
    ;endif
    
    ; save galaxy catalog
    save,galaxyLen,galaxyOff,galaxyIDs,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(workingPath))+' ['+str(countCut)+' / '+str(countMatch)+']'

  endfor ;m
  
end