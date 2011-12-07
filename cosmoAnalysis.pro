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
    sgPIDs = sgPIDList(sg=sg,select='all')
    
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
      print,'Warning: Empty galaxy cut or comp. Skipping: ' + strmid(saveFilename1,strlen(sP.galCatPath))
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
  
    ; set save filename and check for existence
    sP.snap  = m
    saveFilename = sP.galCatPath + 'galradii.' + sP.savPrefix + str(res) + '.' + str(m) + '.sav'
    
    if (file_test(saveFilename)) then begin
      print,'Skipping: '+strmid(saveFilename,strlen(sP.galCatPath))
      continue
    endif
    
    ; load galaxy and group membership catalogs
    gc = galaxyCat(res=res,run=run,snap=m)
    
    ; restrict gas particle positions to ids_ind
    ids  = loadSnapshotSubset(sP.simPath,snapNum=m,partType='gas',field='ids')
    
    ; IMPORTANT! rearrange ids_ind to be in the order of gc.xIDs, need this if we want ids[ids_ind], 
    ; temp[ids_ind], etc to be in the same order as the galaxy catalog id list    
    match,gc.galaxyIDs,ids,gc_ind,ids_gal_ind,count=countGal,/sort
    ids_gal_ind = ids_gal_ind[sort(gc_ind)]
    
    match,gc.groupmemIDs,ids,gc_ind,ids_gmem_ind,count=countGmem,/sort
    ids_gmem_ind = ids_gmem_ind[sort(gc_ind)]
    
    ids    = !NULL
    gc_ind = !NULL
    
    ; calculate radial distances of gas elements to primary and secondary parents
    pos = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='pos')
    
    pos_gal  = pos[*,ids_gal_ind]
    pos_gmem = pos[*,ids_gmem_ind]
    
    pos = !NULL
    
    ; load subhalo catalog for mostBoundParticleID and for sgIDListPri
    sg  = loadSubhaloGroups(sP.simPath,sP.snap)
    
    ; find group center positions with most bound particles for each group
    groupCen = groupCenterPosByMostBoundID(sP=sP,sg=sg)
    
    ; arrays to hold r for all remaining gas particles
    rvec_gal  = fltarr(3,countGal)
    rvec_gmem = fltarr(3,countGmem)
    
    ; replicate parent IDs (of SECONDARY/direct)
    sgInd = galCatRepParentIDs(gc=gc)
    
    ; if requested, match this snapshot galCat against the IDs of another galCat
    
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
    sgIDListPri = sgIDList(sg=sg,select='pri')
    sg = !NULL
    
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
    
    ; save radial distances (and group centers)
    save,gal_pri,gal_sec,gmem_pri,gmem_sec,groupCen,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.galCatPath))

  endfor
  
end

; gasOrigins(): from a target redshift load the galaxy catalog and consider the evolution of all the gas 
;               elements snapshot by snapshot backwards in time. at each step save:
;                 1. the radial distance of each from its primary/secondary parent
;                 2. temperature and entropy
;                 3. implicitly, whether the gas is in a primary or secondary subhalo (corresponding to 
;                    the two radial distances being equal or different)              

function gasOrigins, res=res, run=run, snap=snap

  sP = simParams(res=res,run=run)

  ; config
  redshift = 2.0
  
  targetSnap = redshiftToSnapNum(redshift)
  
  ; if snap specified, run only one snapshot (and/or just return previous results)
  if (keyword_set(snap)) then begin
    saveFilename = sP.derivPath + 'gas.origins.'+str(res)+'.'+str(targetSnap)+'-'+str(snap)+'.sav'
    
    ; results exist, return
    if (file_test(saveFilename)) then begin
      restore,saveFilename
      r = {temp_gal:temp_gal,temp_gmem:temp_gmem,entropy_gal:entropy_gal,entropy_gmem:entropy_gmem,$
           gal_pri:gal_pri,gal_sec:gal_sec,gmem_pri:gmem_pri,gmem_sec:gmem_sec}
      return,r
    endif
    
    ; need to compute, set restricted range of snapshots to process
    snapRange = [snap,snap]
  endif else begin
    ; default config
    numSnapsBack = 5
    
    snapRange = [targetSnap,targetSnap-numSnapsBack]
  endelse  
  
  ; load galaxy catalog at target redshift
  gc = galaxyCat(res=res,run=run,snap=targetSnap)
  
  print,'Loaded  ['+str(n_elements(gc.galaxyIDs))+'] ['+str(n_elements(gc.groupmemIDs))+'] from galCat.'
  
  for m=snapRange[0],snapRange[1],-1 do begin
  
    ; set save filename and check for existence
    sP.snap = m
    saveFilename = sP.derivPath + 'gas.origins.'+str(res)+'.'+str(targetSnap)+'-'+str(m)+'.sav'
    
    if (file_test(saveFilename)) then begin
      print,'Skipping: '+strmid(saveFilename,strlen(sP.derivPath))
      continue
    endif  
  
    ; load gas IDs and match
    ids  = loadSnapshotSubset(sP.simPath,snapNum=m,partType='gas',field='ids')
    
    ; IMPORTANT! rearrange ids_ind to be in the order of gc.xIDs, need this if we want ids[ids_ind], 
    ; temp[ids_ind], etc to be in the same order as the galaxy catalog id list     
    match,gc.galaxyIDs,ids,gc_ind,ids_gal_ind,count=countGal,/sort
    ids_gal_ind = ids_gal_ind[sort(gc_ind)]
    
    match,gc.groupmemIDs,ids,gc_ind,ids_gmem_ind,count=countGmem,/sort
    ids_gmem_ind = ids_gmem_ind[sort(gc_ind)]
    
    ids    = !NULL
    gc_ind = !NULL
    
    ; load u,nelec and calculate temp of gas
    u     = loadSnapshotSubset(sP.simPath,snapNum=m,partType='gas',field='u')
    nelec = loadSnapshotSubset(sP.simPath,snapNum=m,partType='gas',field='nelec')

    temp_gal  = convertUtoTemp(u[ids_gal_ind], nelec[ids_gal_ind])
    temp_gmem = convertUtoTemp(u[ids_gmem_ind],nelec[ids_gmem_ind])
    
    nelec = !NULL
    
    ; load gas density to calculate entropy
    dens = loadSnapshotSubset(sP.simPath,snapNum=m,partType='gas',field='density')
    
    entropy_gal  = calcEntropy(u[ids_gal_ind], dens[ids_gal_ind])
    entropy_gmem = calcEntropy(u[ids_gmem_ind],dens[ids_gmem_ind])
    
    u    = !NULL
    dens = !NULL
    
    ; load the galaxy catalog at this redshift and match IDs from the galCat at the target redshift
    gcCur = galaxyCat(res=res,run=run,snap=m)
    
    match,gc.galaxyIDs,gcCur.galaxyIDs,gc_ind_gal,gcCur_ind_gal,count=countGal,/sort
    match,gc.galaxyIDs,gcCur.groupmemIDs,gc_ind_gal2,gcCur_ind_gal2,count=countGal2,/sort
    match,gc.groupmemIDs,gcCur.groupmemIDs,gc_ind_gmem,gcCur_ind_gmem,count=countGmem,/sort
    match,gc.groupmemIDs,gcCur.galaxyIDs,gc_ind_gmem2,gcCur_ind_gmem2,count=countGmem2,/sort
    
    gcCur = !NULL
    
    print,'['+str(targetSnap-m)+'] Matched ['+str(countGal)+' + '+str(countGal2)+'] ['+$
          str(countGmem)+' + '+str(countGmem2)+'].'
    
    ; allocate space for radial distances
    gal_pri = fltarr(n_elements(gc.galaxyIDs)) - 1.0
    gal_sec = fltarr(n_elements(gc.galaxyIDs)) - 1.0
    
    gmem_pri = fltarr(n_elements(gc.groupmemIDs)) - 1.0
    gmem_sec = fltarr(n_elements(gc.groupmemIDs)) - 1.0
    
    ; load the galaxy radii catalog at this redshift
    galRad = galaxyCatRadii(res=res,run=run,snap=m)

    ; store radial distances of gas elements to primary and secondary parents for matching IDs
    gal_pri[gc_ind_gal]  = galRad.gal_pri[gcCur_ind_gal]
    gal_sec[gc_ind_gal]  = galRad.gal_sec[gcCur_ind_gal]
    
    if (countGal2 gt 0) then begin
      gal_pri[gc_ind_gal2] = galRad.gmem_pri[gcCur_ind_gal2]
      gal_sec[gc_ind_gal2] = galRad.gmem_sec[gcCur_ind_gal2]
    endif
      
    gmem_pri[gc_ind_gmem]  = galRad.gmem_pri[gcCur_ind_gmem]
    gmem_sec[gc_ind_gmem]  = galRad.gmem_sec[gcCur_ind_gmem]
    
    if (countGmem2 gt 0) then begin
      gmem_pri[gc_ind_gmem2] = galRad.gal_pri[gcCur_ind_gmem2]
      gmem_sec[gc_ind_gmem2] = galRad.gal_sec[gcCur_ind_gmem2]
    endif

    ; save
    save,temp_gal,temp_gmem,entropy_gal,entropy_gmem,gal_pri,gal_sec,gmem_pri,gmem_sec,filename=saveFilename
    print,'    Saved: '+strmid(saveFilename,strlen(sP.derivPath))

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



; gcSubsetProp(): read galaxy catalog for a specific subgroup selection (pri,sec,all) and
;                 return properties for each gas element (may or may not depend on parent halo)
;
; rVirNorm=1    : radial distances normalized by r_vir of either primary or secondary parent
; virTemp=1     : virial temperatures of parent halos
; parMass=1     : total mass (dm+baryonic) of parent halos (from catalog)
; curTemp=1     : current temperature of each element
; curDens=1     : current density of each gas element
; maxPastTemp=1 : maximum past previos temperature of each element

function gcSubsetProp, sP=sP, select=select, $
                       rVirNorm=rVirNorm, virTemp=virTemp, parMass=parMass, $
                       curTemp=curTemp, curDens=curDens, maxPastTemp=maxPastTemp, $
                       parNorm=parNorm, minNumPart=minNumPart ;this row for rVirNorm only

  ; select primary,secondary,or all subhalos subject to minimum number of particles
  sgIDList = sgIDList(sP=sP,select=select,minNumPart=minNumPart)

  ; select galaxycat indices corresponding to this list of subgroup ids
  gcInds = gcINDList(sP=sP,sgIDList=sgIDList)

  if keyword_set(rVirNorm) then begin
    ; load parent r_vir and galaxy radii catalog
    r_vir = galCatParentProperties(sP=sP, /rVir)
    gcr   = galaxyCatRadii(res=sP.res,run=sP.run,snap=sP.snap)
    
    ; create subsets for subhalo selection and normalized by parent r_vir
    if (parNorm eq 'pri') then begin
      rad_gal  = gcr.gal_pri[gcInds.gal] / r_vir.gal[gcInds.gal]
      rad_gmem = gcr.gmem_pri[gcInds.gmem] / r_vir.gmem[gcInds.gmem]
    endif
    
    if (parNorm eq 'sec') then begin
      rad_gal  = gcr.gal_sec[gcInds.gal] / r_vir.gal[gcInds.gal]
      rad_gmem = gcr.gmem_sec[gcInds.gmem] / r_vir.gmem[gcInds.gmem]
    endif
    
    r = {gal:rad_gal,gmem:rad_gmem}
  endif
  
  if keyword_set(virTemp) then begin
    ; load parent t_vir
    t_vir = galCatParentProperties(sP=sP, /virTemp)
    
    ; create subsets for subhalo selection
    temp_gal  = t_vir.gal[gcInds.gal]
    temp_gmem = t_vir.gmem[gcInds.gmem]
    
    r = {gal:temp_gal,gmem:temp_gmem}
  endif
  
  if keyword_set(parMass) then begin
    ; load parent total mass
    mass = galCatParentProperties(sP=sP, /mass)
    
    ; create subsets for subhalo selection
    mass_gal  = mass.gal[gcInds.gal]
    mass_gmem = mass.gmem[gcInds.gmem]
    
    r = {gal:mass_gal,gmem:mass_gmem}
  endif
  
  if keyword_set(curTemp) then begin
    ; load intermediate save for thermal state
    thFilename = sP.thistPath + 'thermhist.gas.'+str(sP.res)+'_'+str(sP.snap)+'.sav'
    restore,thFilename
    
    ; load galaxy catalog to change INDs to gas IDs
    gc = galaxyCat(res=sP.res,run=sP.run,snap=sP.snap)
    
    ; restrict temps to subset of galaxy cat
    temp_gal  = temp[gc.galaxyIDs[gcInds.gal]]
    temp_gmem = temp[gc.groupmemIDs[gcInds.gmem]]
    
    gc      = !NULL
    density = !NULL
    temp    = !NULL
    gcInds  = !NULL
      
    ; take log of temperatures
    w = where(temp_gal le 0,count)
    if (count ne 0) then temp_gal[w] = 1.0
    w = where(temp_gmem le 0,count)
    if (count ne 0) then temp_gmem[w] = 1.0
    
    temp_gal  = alog10(temp_gal)
    temp_gmem = alog10(temp_gmem)
    
    r = {gal:temp_gal,gmem:temp_gmem}
  endif
  
  if keyword_set(curDens) then begin
    ; load intermediate save for gas densities
    thFilename = sP.thistPath + 'thermhist.gas.'+str(sP.res)+'_'+str(sP.snap)+'.sav'
    restore,thFilename
    
    ; load galaxy catalog to change INDs to gas IDs
    gc = galaxyCat(res=sP.res,run=sP.run,snap=sP.snap)

    ; restrict temps to subset of galaxy cat
    dens_gal  = density[gc.galaxyIDs[gcInds.gal]]
    dens_gmem = density[gc.groupmemIDs[gcInds.gmem]]

    gc      = !NULL
    density = !NULL
    temp    = !NULL
    gcInds  = !NULL
      
    r = {gal:dens_gal,gmem:dens_gmem}
  endif
  
  if keyword_set(maxPastTemp) then begin
    redshift = snapNumToRedshift(snap=sP.snap)
    
    ; load maximum past temperature
    maxt = maxTemps(res=sP.res,run=sP.run,redshift=redshift)
    
    ; restrict temps to subset of galaxy cat
    temp_gal  = maxt.maxTemps_gal[gcInds.gal]
    temp_gmem = maxt.maxTemps_gmem[gcInds.gmem]
    
    maxt = !NULL
    
    r = {gal:temp_gal,gmem:temp_gmem}
  endif
  
  return,r
end

@cosmoPlot
