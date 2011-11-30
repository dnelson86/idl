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
;
; note: changed to scale Torrey+ galaxy cut to physical density by default

function galaxyCat, res=res, run=run, snap=snap

  sP = simParams(res=res,run=run)

  ; config
  galcut_T   = 6.0
  galcut_rho = 0.25

  ; if snap specified, run only one snapshot (and/or just return previous results)
  if (keyword_set(snap)) then begin
    saveFilename1 = sP.galCatPath + 'galcat.' + sP.savPrefix + str(res) + '.' + str(snap) +'.sav'
    saveFilename2 = sP.galCatPath + 'groupmemcat.' + sP.savPrefix + str(res) + '.' + str(snap) + '.sav'
    
    ; results exist, return
    if (file_test(saveFilename1) and file_test(saveFilename2)) then begin
      restore,saveFilename1
      restore,saveFilename2
      r = {galaxyOff:galaxyOff,galaxyLen:galaxyLen,galaxyIDs:galaxyIDs,$
           groupmemOff:groupmemOff,groupmemLen:groupmemLen,groupmemIDs:groupmemIDs}
      return,r
    endif
    
    ; need to compute, set restricted range of snapshots to process
    snapRange = [snap,snap]
  endif else begin
    ;snapRange = [210,219] ;running 512 gadget
    snapRange = sP.groupCatRange
  endelse
  
  for m=snapRange[0],snapRange[1],1 do begin
    ; skip if previous results exist
    saveFilename1 = sP.galCatPath + 'galcat.' + sP.savPrefix + str(res) + '.' + str(m) + '.sav'
    saveFilename2 = sP.galCatPath + 'groupmemcat.' + sP.savPrefix + str(res) + '.' + str(m) + '.sav'
    
    if (file_test(saveFilename1) and file_test(saveFilename2)) then begin
      print,'Skipping: '+strmid(saveFilename1,strlen(sP.galCatPath))
      continue
    endif
    
    ; load ids of particles in all subfind groups
    sg = loadSubhaloGroups(sP.simPath,m)
    sgPIDs = sgPIDList(sg=sg,/all)
    
    ; load gas ids and match to catalog
    ids  = loadSnapshotSubset(sP.simPath,snapNum=m,partType='gas',field='ids')
    match,sgPIDs,ids,sg_ind,ids_ind,count=countMatch,/sort
    ids_ind = ids_ind[sort(sg_ind)] ; IMPORTANT! rearrange ids_ind to be in the order of sgPIDs
                                    ; need this if we want ids[ids_ind], temp[ids_ind], etc to be
                                    ; in the same order as the group catalog id list

    sgPIDs = !NULL
    sg_ind = !NULL
    ids = ids[ids_ind]

    ; load u,nelec and calculate temp of gas
    u = loadSnapshotSubset(sP.simPath,snapNum=m,partType='gas',field='u')
    u = u[ids_ind]
    nelec = loadSnapshotSubset(sP.simPath,snapNum=m,partType='gas',field='nelec')
    nelec = nelec[ids_ind]
    
    temp = convertUtoTemp(u,nelec)
    
    u     = !NULL
    nelec = !NULL
    
    ; load rho of gas and make galaxy (rho,temp) plane cut
    dens = loadSnapshotSubset(sP.simPath,snapNum=m,partType='gas',field='density')
    dens = dens[ids_ind]
    
    ; scale Torrey+ (2011) galaxy cut to physical density
    scalefac = snapNumToRedshift(snap=m,/time) ; time flag gives simulation time = scale factor
    a3inv = 1.0 / (scalefac*scalefac*scalefac)
    dens *= a3inv
    
    w = where(alog10(temp) - galcut_rho * alog10(dens) lt galcut_T,$
              countCut,comp=wComp,ncomp=countComp)

    if (countCut eq 0 or countComp eq 0) then begin
      print,'Warning: Empty galaxy cut or comp. Skipping: ' + strmid(saveFilename,strlen(sP.galCatPath))
      continue
    endif
    
    temp = !NULL
    dens = !NULL
    
    ; make subsets of ids matching galaxy cut and complement
    ids_groupmem = ids[wComp]
    ids = ids[w]    

    ; construct group member catalog
    if (not file_test(saveFilename2)) then begin
    
      groupmemLen = ulonarr(sg.nSubgroupsTot)
      groupmemOff = ulonarr(sg.nSubgroupsTot)
      groupmemIDs = ulonarr(n_elements(ids_groupmem))
      
      nextOff = 0UL
      
      ; match indices between gas ids and group member ids
      match,ids_groupmem,sg.subgroupIDs,ids_ind,sgIDs_ind,count=countID,/sort
      
      for sgID=0,sg.nSubgroupsTot-1 do begin
        ; select ids in group
        groupStart = sg.subGroupOffset[sgID]
        groupEnd   = groupStart + sg.subGroupLen[sgID]
        
        w = where(sgIDs_ind ge groupStart and sgIDs_ind lt groupEnd,count)
        
        ; save in similar 1D offset format
        if (count gt 0) then begin
          groupmemLen[sgID] = count
          groupmemOff[sgID] = nextOff
          nextOff = groupmemOff[sgID] + groupmemLen[sgID]
    
          groupmemIDs[groupmemOff[sgID]:nextOff-1] = ids_groupmem[ids_ind[w]]
        endif      
        
      end
  
      ; debug: make sure all gas particles were found in the group catalog
      ;match,ids_groupmem,groupmemIDs,ind1,ind2,count=countCheck
      ;if (countCheck ne n_elements(ids_groupmem)) then begin
      ;  print,'Uhoh, check2.',countCheck,n_elements(ids_groupmem)
      ;  nuniq = n_elements(uniq(ids_groupmem[sort(ids_groupmem)]))
      ;  stop
      ;endif

      ; save group membership catalog
      save,groupmemLen,groupmemOff,groupmemIDs,filename=saveFilename2    
      print,'Saved: '+strmid(saveFilename2,strlen(sP.galCatPath))+' ['+str(countComp)+'/'+str(countMatch)+']'    
    
    endif
    
    ; construct galaxy catalog
    if (not file_test(saveFilename1)) then begin
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
      
      ; save galaxy catalog
      save,galaxyLen,galaxyOff,galaxyIDs,filename=saveFilename1
      print,'Saved: '+strmid(saveFilename1,strlen(sP.galCatPath))+' ['+str(countCut)+'/'+str(countMatch)+']'
    endif

  endfor ;m
  
end

; gasOrigins(): loop through all snapshots from the beginning of the group catalogs (z=6) and
;               categorize gas into different modes/origins:
;                   
;               1. smooth (in the DM sense) - no membership in any subfind group at any previous time
;               2-4. TODO

function gasOrigins, res=res, run=run

  sP = simParams(res=res,run=run)

  targetRedshifts = [3.0,2.0,1.0,0.0]
  targetSnapshots = redshiftToSnapNum(targetRedshifts)

  ; masks
  h = loadSnapshotHeader(sP.simPath,snapNum=0)
  
  nTot   = h.nPartTot[0] + h.nPartTot[1] + h.nPartTot[3] ; gas + dm + stars
  padFac = 1.1 ; allow 10% more ids to exist due to refinement (Arepo)
  nMask  = round(nTot * padFac)
  
  cat1Mask = bytarr(nMask) ;smooth
  ;cat2Mask = bytarr(nMask) ;quasi-smooth
  ;cat3Mask = bytarr(nMask) ;merger/clumpy
  ;cat4Mask = bytarr(nMask) ;recycled/stripped
  
  for m=sP.groupCatRange[0],sP.groupCatRange[1],1 do begin
  
    ; load catalog
    gc = galaxyCat(res=res,run=run,snap=m)
    
    ; if target redshift, save results
    w = where(targetSnapshots eq m,count)
    if (count ne 0) then begin
        ; find all particle ids that have not been added to masks
        pIDs_Cat1 = where(cat1Mask eq 0B,count1)
        
        ; match with ids in catalog at this snapshot
        match,pIDs_Cat1,gc.galaxyIDs,pids_ind_gal,gc_ind_gal,count=count_gal,/sort
        match,pIDs_Cat1,gc.groupmemIDs,pids_ind_gmem,gc_ind_gmem,count=count_gmem,/sort
        
        pIDs_Cat1_gal  = pIDs_Cat1[pids_ind_gal]
        pIDs_Cat1_gmem = pIDs_Cat1[pids_ind_gmem]
        
        ; set saveFilename and save
        saveFilename = sP.derivPath + 'gas.smooth.'+str(res)+'.'+str(sP.snapRange[0])+'-'+str(m)+'.sav'

        save,pIDs_Cat1_gal,pIDs_Cat1_gmem,filename=saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    endif
    
    ; check overflow on masks
    if (max(gc.galaxyIDs) ge nMask or max(gc.groupmemIDs) ge nMask) then begin
      print,'ERROR'
      stop
    endif
    
    ; add catalog members to masks
    cat1Mask[gc.galaxyIDs]   = 1B
    cat1Mask[gc.groupmemIDs] = 1B
  
  endfor

end

; galaxyCatRadii(): find radial distance of all group member particles wrt the group they belong to
;                   as well as the rad to the primary group if this is a secondary group

function galaxyCatRadii, res=res, run=run, snap=snap

  sP = simParams(res=res,run=run)

  ; if snap specified, run only one snapshot (and/or just return previous results)
  if (keyword_set(snap)) then begin
    saveFilename = sP.galCatPath + 'galradii.' + sP.savPrefix + str(res) + '.' + str(snap) + '.sav'
    
    ; results exist, return
    if (file_test(saveFilename)) then begin
      restore,saveFilename
      r = {gal_pri:gal_pri,gal_sec:gal_sec,gmem_pri:gmem_pri,gmem_sec:gmem_sec}
      return,r
    endif
    
    ; need to compute, set restricted range of snapshots to process
    snapRange = [snap,snap]
  endif else begin
    snapRange = sP.groupCatRange
  endelse
  
  ; loop from target redshift to beginning of group catalogs
  for m=snapRange[0],snapRange[1],1 do begin
    saveFilename = sP.galCatPath + 'galradii.' + sP.savPrefix + str(res) + '.' + str(m) + '.sav'
    
    if (file_test(saveFilename)) then begin
      print,'Skipping: '+strmid(saveFilename,strlen(sP.galCatPath))
      continue
    endif
    
    ; load ids of particles in all subfind groups
    sg = loadSubhaloGroups(sP.simPath,m)
    sgPIDs = sgPIDList(sg=sg,/all)
    sgIDListPri = sgIDList(sg=sg,/pri)
    
    sg = !NULL
  
    ; load galaxy and group membership catalogs
    gc = galaxyCat(res=res,run=run,snap=m)
    
    ; load group catalog for positions and gas particle positions
    sg  = loadSubhaloGroups(sP.simPath,m,/skipIDs)
    
    ; restrict gas particle positions to ids_ind
    ids  = loadSnapshotSubset(sP.simPath,snapNum=m,partType='gas',field='ids')
    
    ; IMPORTANT! rearrange ids_ind to be in the order of sgPIDs, need this if we want ids[ids_ind], 
    ; temp[ids_ind], etc to be in the same order as the group catalog id list    
    match,gc.galaxyIDs,ids,gc_ind,ids_gal_ind,count=countGal,/sort
    ids_gal_ind = ids_gal_ind[sort(gc_ind)]
    
    match,gc.groupmemIDs,ids,gc_ind,ids_gmem_ind,count=countGmem,/sort
    ids_gmem_ind = ids_gmem_ind[sort(gc_ind)]
    
    ids    = !NULL
    gc_ind = !NULL
    
    pos = loadSnapshotSubset(sP.simPath,snapNum=m,partType='gas',field='pos')
    
    pos_gal  = pos[*,ids_gal_ind]
    pos_gmem = pos[*,ids_gmem_ind]
    
    pos = !NULL
    
    ; find positions of most bound particles for each group
    sP.snap  = m
    groupCen = groupCenterPosByMostBoundID(sP=sP,sg=sg)
    
    ; arrays to hold r for all remaining gas particles
    rvec_gal  = fltarr(3,countGal)
    rvec_gmem = fltarr(3,countGmem)
    
    ; replicate parent IDs (of SECONDARY/direct)
    sgInd = galCatRepParentIDs(gc=gc)
    
    ; calulate radial vector of gas from group center, correct for periodic B.C.
    rvec_gal  = groupCen[*,sgInd.gal] - pos_gal
    rvec_gmem = groupCen[*,sgInd.gmem] - pos_gmem

    correctPeriodicDistVecs, rvec_gal
    correctPeriodicDistVecs, rvec_gmem

    ; calculate radial distance
    gal_sec = reform( sqrt(rvec_gal[0,*]*rvec_gal[0,*] + $
                           rvec_gal[1,*]*rvec_gal[1,*] + $
                           rvec_gal[2,*]*rvec_gal[2,*]) )
    
    rvec_gal = !NULL
    
    gmem_sec = reform( sqrt(rvec_gmem[0,*]*rvec_gmem[0,*] + $
                            rvec_gmem[1,*]*rvec_gmem[1,*] + $
                            rvec_gmem[2,*]*rvec_gmem[2,*]) )
    
    rvec_gmem = !NULL
    
    ; replicate parent IDs (of PRIMARY/parent)
    sgInd = galCatRepParentIDs(gc=gc,sgIDListPri=sgIDListPri)
    
    ; calulate radial vector of gas from group center, correct for periodic B.C.
    rvec_gal  = groupCen[*,sgInd.gal] - pos_gal
    rvec_gmem = groupCen[*,sgInd.gmem] - pos_gmem

    correctPeriodicDistVecs, rvec_gal
    correctPeriodicDistVecs, rvec_gmem

    ; calculate radial distance
    gal_pri = reform( sqrt(rvec_gal[0,*]*rvec_gal[0,*] + $
                           rvec_gal[1,*]*rvec_gal[1,*] + $
                           rvec_gal[2,*]*rvec_gal[2,*]) )
    
    rvec_gal = !NULL
    
    gmem_pri = reform( sqrt(rvec_gmem[0,*]*rvec_gmem[0,*] + $
                            rvec_gmem[1,*]*rvec_gmem[1,*] + $
                            rvec_gmem[2,*]*rvec_gmem[2,*]) )
    
    rvec_gmem = !NULL
    
    ; save galaxy catalog
    save,gal_pri,gal_sec,gmem_pri,gmem_sec,groupCen,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.galCatPath))

  endfor
  
end

; maxTemps(): find maximum temperature for gas particles in galaxy/group member catalogs at redshift
;             through the redshift range (redshift,zStart] where zStart is typically the start of 
;             the simulation
;
; note: currently temps are only saved for gas in groups at the end of the interval, which will
;       prevent us from making the maxtemp vs. redshift evaluated at z=0 plot
;
; inclEffEOS=1 : include temperature of gas when on the effective EOS (star forming)

function maxTemps, res=res, run=run, redshift=redshift, zStart=zStart, inclEffEOS=inclEffEOS

  sP = simParams(res=res,run=run)

  ; restrict temp array to redshift range
  if not keyword_set(zStart) then zStart = 30.0
  
  maxSnap = redShiftToSnapnum(redshift) - 1 ; do not include current time in Tmax search
  minSnap = redShiftToSnapnum(zStart)

  ; set saveFilename and check for existence
  saveFilename = sP.derivPath + 'maxtemp.'+str(res)+'.'+str(minSnap)+'-'+str(maxSnap)+'.sav'
  
  if (keyword_set(inclEffEOS)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.inclEffEOS.sav'
      
  if (file_test(saveFilename)) then begin
    restore, saveFilename
    r = {maxTemps_gal:maxTemps_gal,maxTemps_gmem:maxTemps_gmem,$
         maxTempSnap_gal:maxTempSnap_gal,maxTempSnap_gmem:maxTempSnap_gmem}
    return, r
  endif
  
  print,'Calculating new maxtemp for res = '+str(res)+' in range ['+str(minSnap)+'-'+str(maxSnap)+'].'
  if keyword_set(inclEffEOS) then print,' Including gas temperature on effective EOS!'
  
  ; load galaxy/group member catalogs at zMin for gas ids to search for
  gc = galaxyCat(res=res,run=run,snap=maxSnap+1)

  maxTemps_gal      = fltarr(n_elements(gc.galaxyIDs))
  maxTempSnap_gal   = intarr(n_elements(gc.galaxyIDs))
  maxTemps_gmem     = fltarr(n_elements(gc.groupmemIDs))
  maxTempSnap_gmem  = intarr(n_elements(gc.groupmemIDs))
  
  for m=minSnap,maxSnap,1 do begin
    ; load thermal state
    thFilename = sP.thistPath + 'thermhist.gas.'+str(res)+'_'+str(m)+'.sav'
    restore,thFilename
    
    ; restriction from effective EOS
    if not keyword_set(inclEffEOS) then begin
      ; load gas ids and match to catalog
      ids  = loadSnapshotSubset(sP.simPath,snapNum=m,partType='gas',field='ids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of sgPIDs, need this if we want ids[ids_ind], 
      ; temp[ids_ind], etc to be in the same order as the group catalog id list    
      match,gc.galaxyIDs,ids,gc_ind,ids_gal_ind,count=countGal,/sort
      ids_gal_ind = ids_gal_ind[sort(gc_ind)]
      
      match,gc.groupmemIDs,ids,gc_ind,ids_gmem_ind,count=countGmem,/sort
      ids_gmem_ind = ids_gmem_ind[sort(gc_ind)]
      
      ids    = !NULL
      gc_ind = !NULL
      
      ; load gas SFR and select for any non-zero
      sfr  = loadSnapshotSubset(sP.simPath,snapNum=m,partType='gas',field='sfr')
      
      sfr_gal  = sfr[ids_gal_ind]
      sfr_gmem = sfr[ids_gmem_ind]
      
      sfr = !NULL
      
    endif else begin
      ; if not restricting on effective EOS, create true arrays for where selects
      sfr_gal  = bytarr(n_elements(gc.galaxyIDs))
      sfr_gmem = bytarr(n_elements(gc.groupmemIDs))
    endelse
    
    ; take log of temperatures
    w = where(temp le 0,count)
    if (count ne 0) then temp[w] = 1.0
    
    temp = alog10(temp)
    
    ; calc max temp and time for galaxy members
    w1 = where(temp[gc.galaxyIDs] gt maxTemps_gal and sfr_gal eq 0.0,count1)
    if (count1 gt 0) then begin
      maxTemps_gal[w1]    = temp[gc.galaxyIDs[w1]]
      maxTempSnap_gal[w1] = m
    endif
    
    ; calc max temp and time for group members
    w2 = where(temp[gc.groupmemIDs] gt maxTemps_gmem and sfr_gmem eq 0.0,count2)
    if (count2 gt 0) then begin
      maxTemps_gmem[w2]    = temp[gc.groupmemIDs[w2]]
      maxTempSnap_gmem[w2] = m
    endif
    
  endfor ;m

  ; save for future lookups
  save,maxTemps_gal,maxTemps_gmem,maxTempSnap_gal,maxTempSnap_gmem,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))

  r = {maxTemps_gal:maxTemps_gal,maxTemps_gmem:maxTemps_gmem,$
       maxTempSnap_gal:maxTempSnap_gal,maxTempSnap_gmem:maxTempSnap_gmem}
  return,r

end

@cosmoPlot
