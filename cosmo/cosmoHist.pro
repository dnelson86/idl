; cosmoHist.pro
; gas accretion project - past history of gas (all routines that step through multiple snapshots)
; dnelson apr.2012

; -----------------------------------------------------------------------------------------------------
; accretionMode(): for eaching gas particle/tracer with a recorded accretion time, starting at some 
;                  redshift, track backwards in time with respect to the tracked parent halos (using 
;                  mergerTree) and classify the mode of accretion by considering its group membership
;                  at the outer bracketing snapshot just prior to accretion (r>~r_vir) as:
;                  
;  1. smooth  = in the (primary) parent of the original halo or not in any subgroup
;  2. sclumpy = in any subgroup other than the (primary) parent halo that is a secondary subgroup ("small")
;  3. bclumpy = in any subgroup other than the (primary) parent halo that is a primary subgroup ("big")
; -----------------------------------------------------------------------------------------------------

function accretionMode, sP=sP

  forward_function cosmoTracerChildren, cosmoTracerVelParents, accretionTimes
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; first, walk back through the merger tree and find primary subhalos with good parent histories
  mt = mergerTreeSubset(sP=sP,/verbose)
  at = accretionTimes(sP=sP)

  ; set saveFilename and check for existence
  saveTag = ''
  if sP.trMCPerCell eq -1 then saveTag = '.trVel'
  if sP.trMCPerCell gt 0  then saveTag = '.trMC'
  if sP.trMCPerCell eq 0  then saveTag = '.SPH'
  
  saveFilename = sP.derivPath + 'accMode'+saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(mt.maxSnap)+'-'+str(mt.minSnap)+'.sav'  
  
  if file_test(saveFilename) then begin
    restore, saveFilename
    return, r
  endif
  
  ; select those particles/tracers with recorded accretion times
  gal_w_at  = where(at.AccTime_gal ne -1,count_gal)
  gmem_w_at = where(at.AccTime_gmem ne -1,count_gmem)
  
  ; load galaxy/group member catalogs at zMin for gas ids to search for
  galcat = galaxyCat(sP=sP)  
  
  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    print,'Calculating new accretion mode using ( SPH Particles ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
      
    galcatSub = { gal  : mt.galcatSub.gal[gal_w_at]    ,$
                  gmem : mt.galcatSub.gmem[gmem_w_at]   }
      
    ; convert the accretion time (scale factor) into the outer bracketing snapshot number for each member
    bracketSnap = { gal  : intarr(n_elements(galcatSub.gal))   ,$
                    gmem : intarr(n_elements(galcatSub.gmem))   }
                    
    snapTimes = snapNumToRedshift(sP=sP,/all,/time)
    snapTimes = snapTimes[where(snapTimes ne -1)]

    bracketSnap.gal  = value_locate(snapTimes,at.accTime_gal[gal_w_at])
    bracketSnap.gmem = value_locate(snapTimes,at.accTime_gmem[gmem_w_at])
    
    w = where(bracketSnap.gal eq -1 or bracketSnap.gal eq n_elements(snapTimes)-1,count)
    if count then message,'Error: Bad gal bracketing.'
    w = where(bracketSnap.gmem eq -1 or bracketSnap.gmem eq n_elements(snapTimes)-1,count)
    if count then message,'Error: Bad gmem bracketing.'
    
    ; store the main arrays as a structure so we can write them directly
    r = {accMode_gal       : intarr(n_elements(galcatSub.gal))-1  ,$
         accMode_gmem      : intarr(n_elements(galcatSub.gmem))-1  }
    
    for m=mt.maxSnap,mt.minSnap,-1 do begin
      sP.snap = m
      ; load local group catalog
      h = loadSnapshotHeader(sP=sP)
      gc = loadGroupCat(sP=sP,/readIDs)
      
      ; create a list of IDs of all subgroups (exclude fuzz)
     ; gcGasIDs = gcPIDList(gc=gc,select='all',partType='gas')

      ; select those gas particles whose accretion modes are determined at this snapshot
      w_gal  = where(bracketSnap.gal eq sP.snap,count_gal)
      w_gmem = where(bracketSnap.gmem eq sP.snap,count_gmem)
      
      print,'['+str(m)+'] bracketed gal: '+str(count_gal)+' gmem: '+str(count_gmem)
      
      stop
      
      if count_gal gt 0 then begin
        ; global match against subgroup ID list (this includes fuzz and unbound gas)
        match,galcat.galaxyIDs[galcatSub.gal[w_gal]],gc.IDs,galcat_ind,gc_gal_ind,count=countGal,/sort
        
        ; those that don't match, assign to category 1. smooth
        all = bytarr(count_gal)
        if countGal lt count_gal then all[galcat_ind] = 1B
        w = where(all eq 0B, ncomp)
        
        r.accMode_gal[w_gal[w]] = 1
        print,'['+str(m)+'] immediate smooth: '+str(ncomp)
        
        ; those that match, calculate subgroup ID they belong to
        if countGal gt 0 then begin
          ; calculate parent subgroup ID each gas particle belongs to
          parIDs = value_locate(gc.subgroupOffsetType[partNum('gas'),*],gc_gal_ind)
          
          ; exclude those that aren't actually within the offsetLenType range
          diff = gc_gal_ind - gc.subgroupOffsetType[partNum('gas'),parIDs]
          w = where(diff lt gc.subgroupLenType[partNum('gas'),parIDs],count,comp=wc,ncomp=ncomp)
          
          print,'['+str(m)+'] sg bounds in: '+str(count)+' out: '+str(ncomp)
          
          ; those failing the subgroup bound are smooth
          
          ; loop over each matched particle
          for i=0,countGal-1 do begin
            ; calculate subgroup ID it belongs to
            loc_pos = gc_gal_ind[i]
            ;TODO
            ; if their subgroup ID is the traced parent, assign to category 1. smooth
  
            ; if their subgroup ID is not the parent but is a primary sg, assign to category 3. bclumpy
  
            ; else (subgroup ID is not a primary sg), assign to category 2. sclumpy
        
          endfor
        endif ; countGal
      endif ; count_gal
        
      if count_gmem gt 0 then begin
        ; global match against subgroup ID list
        match,galcat.groupmemIDs[galcatSub.gmem[w_gmem]],gcGasIDs,galcat_ind,gc_gmem_ind,count=countGmem,/sort

        ; those that don't match, assign to category 1. smooth
        all = bytarr(count_gmem)
        if countGal lt count_gmem then all[galcat_ind] = 1B
        w = where(all eq 0B, ncomp)
        
        r.accMode_gmem[w_gmem[w]] = 1
        
        ; those that match, calculate subgroup ID they belong to
        if countGal gt 0 then begin
          ; loop over each matched particle
          for i=0,countGal-1 do begin
            ; calculate subgroup ID it belongs to
            ; TODO
            ; if their subgroup ID is the traced parent, assign to category 1. smooth
  
            ; if their subgroup ID is not the parent but is a primary sg, assign to category 3. bclumpy
  
            ; else (subgroup ID is not a primary sg), assign to category 2. sclumpy
        
          endfor
        endif ; countGmem
      endif ; count_gmem
      
      ; load mergerTree and move to Parent
      
      ; sanity check no IDs are -1 (we should only be processing the mergerTreeSubset)
      
      print,' ['+string(m,format='(i3)')+'] accreted mode counts ',count_smooth,count_sclumpy,count_bclumpy
      
      ; free some memory for next load
      w_gal  = !NULL
      w_gmem = !NULL
    endfor
    
    ; save
    ;save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
  endif  
  
end

; -----------------------------------------------------------------------------------------------------
; accretionTimes(): for each gas particle/tracer, starting at some redshift, track backwards in time
;                   with respect to the tracked parent halos (using mergerTree) and determine the
;                   time when the particle radius = the virial radius (and record the virial temp of
;                   the parent halo at that time)
; -----------------------------------------------------------------------------------------------------

function accretionTimes, sP=sP

  forward_function cosmoTracerChildren, cosmoTracerVelParents
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; first, walk back through the merger tree and find primary subhalos with good parent histories
  mt = mergerTreeSubset(sP=sP,/verbose)

  ; set saveFilename and check for existence
  saveTag = ''
  if sP.trMCPerCell eq -1 then saveTag = '.trVel'
  if sP.trMCPerCell gt 0  then saveTag = '.trMC'
  if sP.trMCPerCell eq 0  then saveTag = '.SPH'
  
  saveFilename = sP.derivPath + 'accTimes'+saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(mt.maxSnap)+'-'+str(mt.minSnap)+'.sav'  
  
  if file_test(saveFilename) then begin
    restore, saveFilename
    return, r
  endif
  
  ; load galaxy/group member catalogs at zMin for gas ids to search for
  galcat = galaxyCat(sP=sP)
  
  accCount = { gal : 0UL, gmem : 0UL } ; count of how many particles/tracers we tracked through r_vir
  prevTime = 0 ; scale factor at previous snapshot

  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    print,'Calculating new accretion time using ( SPH Particles ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
      
    ; store the r/rvir of each at the previous snapshot for interpolation, and as a monotonic check
    prevRad = { gal  : fltarr(n_elements(mt.galcatSub.gal)) ,$
                gmem : fltarr(n_elements(mt.galcatSub.gmem)) }
    
    accMask = { gal  : bytarr(n_elements(mt.galcatSub.gal)), $
                gmem : bytarr(n_elements(mt.galcatSub.gmem)) }
    
    ; store the main arrays as a structure so we can write them directly
    r = {accTime_gal       : fltarr(n_elements(mt.galcatSub.gal))-1  ,$
         accTime_gmem      : fltarr(n_elements(mt.galcatSub.gmem))-1 ,$
         accHaloTvir_gal   : fltarr(n_elements(mt.galcatSub.gal))    ,$
         accHaloTvir_gmem  : fltarr(n_elements(mt.galcatSub.gmem))    }
    
    ; debugging r(t)
    ;radtemp = { gal  : fltarr(mt.maxSnap-mt.minSnap+1,n_elements(mt.galcatSub.gal)) ,$
    ;            gmem : fltarr(mt.maxSnap-mt.minSnap+1,n_elements(mt.galcatSub.gmem)) }
    
    for m=mt.maxSnap,mt.minSnap,-1 do begin
      sP.snap = m
      ; load gas ids and match to catalog
      h = loadSnapshotHeader(sP=sP)
      ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs   
      match,galcat.galaxyIDs[mt.galcatSub.gal],ids,galcat_ind,ids_gal_ind,count=countGal,/sort
      ids_gal_ind = ids_gal_ind[sort(galcat_ind)]
      
      match,galcat.groupmemIDs[mt.galcatSub.gmem],ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
      ids_gmem_ind = ids_gmem_ind[sort(galcat_ind)]
      
      ids        = !NULL
      galcat_ind = !NULL
      
      ; load pos to calculate radii
      pos   = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
      
      pos_gal  = pos[*,ids_gal_ind]
      pos_gmem = pos[*,ids_gmem_ind]
      
      pos = !NULL

      ; calculate current distance of gas particle from smoothed halo center position for galaxy members
      gal_pri  = periodicDists(reform(mt.hPos[mt.maxSnap-m,*,mt.gcIndOrig.gal]),pos_gal,sP=sP)
      gal_pri /= mt.hVirRad[mt.maxSnap-m,mt.gcIndOrig.gal]
      
      ; for group members
      gmem_pri = periodicDists(reform(mt.hPos[mt.maxSnap-m,*,mt.gcIndOrig.gmem]),pos_gmem,sP=sP)
      gmem_pri /= mt.hVirRad[mt.maxSnap-m,mt.gcIndOrig.gmem]
      
      pos_gal  = !NULL
      pos_gmem = !NULL
      
      ; for particles who are past r_vir, sanity check that they are not now within
      ;gal_err  = where(gal_pri lt 1.0 and prevRad.gal ge 1.0,count_gal_err)
      ;gmem_err = where(gmem_pri lt 1.0 and prevRad.gmem ge 1.0,count_gmem_err)
      ;print,' warning counts ',count_gal_err,count_gmem_err
      
      ; for particles who are still within r_vir, check if they have passed beyond
      gal_w  = where(gal_pri ge 1.0 and prevRad.gal lt 1.0 and accMask.gal eq 0B,count_gal)
      gmem_w = where(gmem_pri ge 1.0 and prevRad.gmem lt 1.0 and accMask.gmem eq 0B,count_gmem)
      
      print,' ['+string(m,format='(i3)')+'] accreted now counts '+string(count_gal,format='(i5)')+' ('+$
        string(float(count_gal)/n_elements(gal_pri)*100,format='(f4.1)')+'%) '+$
        string(count_gmem,format='(i5)')+' ('+$
        string(float(count_gmem)/n_elements(gmem_pri)*100,format='(f4.1)')+'%)'
      
      ; interpolate these (time,radii) to find time crossing the virial radius
      times = [prevTime,h.time]
      
      for i=0,count_gal-1 do begin
        radii = [ prevRad.gal[gal_w[i]],gal_pri[gal_w[i]] ]
        time = interpol(times,radii,1.0) ; lerp time to r/rvir=1
        tvir = [ mt.hVirTemp[mt.maxSnap-m-1,mt.gcIndOrig.gal[gal_w[i]]], $
                 mt.hVirTemp[mt.maxSnap-m,mt.gcIndOrig.gal[gal_w[i]]] ]
        tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
        r.accTime_gal[gal_w[i]] = time
        r.accHaloTvir_gal[gal_w[i]] = tvir
      endfor
      
      for i=0,count_gmem-1 do begin
        radii = [ prevRad.gmem[gmem_w[i]],gmem_pri[gmem_w[i]] ]
        time = interpol(times,radii,1.0) ; lerp time to r/rvir=1
        tvir = [ mt.hVirTemp[mt.maxSnap-m-1,mt.gcIndOrig.gmem[gmem_w[i]]], $
                 mt.hVirTemp[mt.maxSnap-m,mt.gcIndOrig.gmem[gmem_w[i]]] ]
        tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
        r.accTime_gmem[gmem_w[i]] = time
        r.accHaloTvir_gmem[gmem_w[i]] = tvir
      endfor
      
      ; if we are on the first snapshot, override accretion times with -1 to indicate always outside rvir
      if m eq mt.maxSnap then r.accTime_gal[gal_w] = -1
      if m eq mt.maxSnap then r.accTime_gmem[gmem_w] = -1
      
      ; update counters for the number of particles we have found the accretion times of
      accMask.gal[gal_w]   = 1B
      accMask.gmem[gmem_w] = 1B
      accCount.gal  += count_gal
      accCount.gmem += count_gmem
      
      ; store current radius of particles
      prevRad.gal  = gal_pri
      prevRad.gmem = gmem_pri
     
      ;radtemp.gal[mt.maxSnap-m,*] = gal_pri ;debug
      ;radtemp.gmem[mt.maxSnap-m,*] = gmem_pri ;debug
      
      prevTime = h.time
      
      ; free some memory for next load
      gal_w    = !NULL
      gmem_w   = !NULL
      gal_pri  = !NULL
      gmem_pri = !NULL
    endfor
    
    print,'found accretion times for ['+str(accCount.gal)+' of '+str(n_elements(mt.galcatSub.gal))+$
      '] gal, ['+str(accCount.gmem)+' of '+str(n_elements(mt.galcatSub.gmem))+'] gmem'
    
    ;save,mt.maxSnap,mt.minSnap,radtemp,filename=sP.plotPath+'temprad.sav' ;debug
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
  endif
  
  ; MONTE CARLO TRACERS CASE - for all gas cells, track back all child tracers
  ; ------------------------
  if sP.trMCPerCell gt 0 then begin

    print,'Calculating new accretion time using ( TracerMC ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
      
    ; load gas ids
    gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')

    ; match galcat IDs to gas_ids
    match,galcat.galaxyIDs[mt.galcatSub.gal],gas_ids,galcat_ind,ids_gal_ind,count=countGal,/sort
    ids_gal = gas_ids[ids_gal_ind[sort(galcat_ind)]]
    
    match,galcat.groupmemIDs[mt.galcatSub.gmem],gas_ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
    ids_gmem = gas_ids[ids_gmem_ind[sort(galcat_ind)]]
    
    gas_ids = !NULL
    
    ; locate tracer children (indices) of gas id subsets
    galcat_gal_trids  = cosmoTracerChildren(sP=sP, /getInds, gasIDs=ids_gal, child_counts=galcat_gal_cc)
    galcat_gmem_trids = cosmoTracerChildren(sP=sP, /getInds, gasIDs=ids_gmem, child_counts=galcat_gmem_cc)
    
    ; convert tracer children indices to tracer IDs at this zMin
    tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
    
    galcat_gal_trids  = tr_ids[galcat_gal_trids]
    galcat_gmem_trids = tr_ids[galcat_gmem_trids]

    tr_ids   = !NULL
    ids_gal  = !NULL
    ids_gmem = !NULL
    
    ; create a gcIndOrig for the tracers
    gcIndOrigTr = galCatRepParentIDs(galcat=galcat,gcIDList=mt.galcatIDList,$
                                     child_counts={gal:galcat_gal_cc,gmem:galcat_gmem_cc}) 
                  
    galcat = !NULL ; not used past this point
    
    ; want to use these parent IDs to access hVirRad,etc so compact the same way (ascending ID->index)
    placeMap = getIDIndexMap(mt.galcatIDList,minid=minid)
    gcIndOrigTr.gal = placeMap[gcIndOrigTr.gal-minid]
    gcIndOrigTr.gmem = placeMap[gcIndOrigTr.gmem-minid]
    placeMap = !NULL
        
    ; store the r/rvir of each at the previous snapshot for interpolation, and as a monotonic check
    prevRad = { gal  : fltarr(n_elements(galcat_gal_trids)) ,$
                gmem : fltarr(n_elements(galcat_gmem_trids)) }
    
    accMask = { gal  : bytarr(n_elements(galcat_gal_trids)), $
                gmem : bytarr(n_elements(galcat_gmem_trids)) }
    
    ; store the main arrays as a structure so we can write them directly for all tracers
    r = {accTime_gal       : fltarr(n_elements(galcat_gal_trids))-1  ,$
         accTime_gmem      : fltarr(n_elements(galcat_gmem_trids))-1 ,$
         accHaloTvir_gal   : fltarr(n_elements(galcat_gal_trids))    ,$
         accHaloTvir_gmem  : fltarr(n_elements(galcat_gmem_trids))   ,$
         child_counts_gal  : galcat_gal_cc                           ,$
         child_counts_gmem : galcat_gmem_cc                           }

    for m=mt.maxSnap,mt.minSnap,-1 do begin
      sP.snap = m
      ; load tracer ids and match to child ids from zMin
      h = loadSnapshotHeader(sP=sP)
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs, need this if we want ids[ids_ind], 
      ; temp[ids_ind], etc to be in the same order as the group catalog id list    
      match,galcat_gal_trids,tr_ids,galcat_ind,trids_gal_ind,count=countGal,/sort
      trids_gal_ind = trids_gal_ind[sort(galcat_ind)]
      
      match,galcat_gmem_trids,tr_ids,galcat_ind,trids_gmem_ind,count=countGmem,/sort
      trids_gmem_ind = trids_gmem_ind[sort(galcat_ind)]
      
      tr_ids     = !NULL
      galcat_ind = !NULL
      
      ; load tracer parents to match to gas
      tr_parids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentid')
      tr_parids_gal  = tr_parids[trids_gal_ind]
      tr_parids_gmem = tr_parids[trids_gmem_ind]
      tr_parids = !NULL
      
      ; load gas IDs and convert tracer parent IDs -> indices
      gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
      gasIDMap = getIDIndexMap(gas_ids,minid=minid)
      gas_ids = !NULL
      
      tr_parids_gal  = gasIDMap[tr_parids_gal-minid]  ; convert ID->index
      tr_parids_gmem = gasIDMap[tr_parids_gmem-minid] ; convert ID->index
      gasIDMap = !NULL
      
      ; load gas positions and convert to tracer positions
      gas_pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
      tr_pos_gal  = gas_pos[*,tr_parids_gal]
      tr_pos_gmem = gas_pos[*,tr_parids_gmem]
      
      gas_pos = !NULL
      tr_parids_gal  = !NULL
      tr_parids_gmem = !NULL

      ; calculate current distance of gas particle from smoothed halo center position for galaxy members
      gal_pri  = periodicDists(reform(mt.hPos[mt.maxSnap-m,*,gcIndOrigTr.gal]),tr_pos_gal,sP=sP)
      gal_pri /= mt.hVirRad[mt.maxSnap-m,gcIndOrigTr.gal]
      
      ; for group members
      gmem_pri = periodicDists(reform(mt.hPos[mt.maxSnap-m,*,gcIndOrigTr.gmem]),tr_pos_gmem,sP=sP)
      gmem_pri /= mt.hVirRad[mt.maxSnap-m,gcIndOrigTr.gmem]
      
      tr_pos_gal  = !NULL
      tr_pos_gmem = !NULL
      
      ; for particles who are past r_vir, sanity check that they are not now within
      gal_err  = where(gal_pri lt 1.0 and prevRad.gal ge 1.0,count_gal_err)
      gmem_err = where(gmem_pri lt 1.0 and prevRad.gmem ge 1.0,count_gmem_err)
      print,' warning counts ',count_gal_err,count_gmem_err
      
      ; for particles who are still within r_vir, check if they have passed beyond
      gal_w  = where(gal_pri ge 1.0 and prevRad.gal lt 1.0 and accMask.gal eq 0B,count_gal)
      gmem_w = where(gmem_pri ge 1.0 and prevRad.gmem lt 1.0 and accMask.gmem eq 0B,count_gmem)
      
      print,' ['+string(m,format='(i3)')+'] accreted now counts '+string(count_gal,format='(i5)')+' ('+$
        string(float(count_gal)/n_elements(gal_pri)*100,format='(f4.1)')+'%) '+$
        string(count_gmem,format='(i5)')+' ('+$
        string(float(count_gmem)/n_elements(gmem_pri)*100,format='(f4.1)')+'%)'
      
      ; interpolate these (time,radii) to find time crossing the virial radius
      times = [prevTime,h.time]
      
      for i=0,count_gal-1 do begin
        radii = [ prevRad.gal[gal_w[i]],gal_pri[gal_w[i]] ]
        time = interpol(times,radii,1.0) ; lerp time to r/rvir=1
        tvir = [ mt.hVirTemp[mt.maxSnap-m-1,gcIndOrigTr.gal[gal_w[i]]], $
                 mt.hVirTemp[mt.maxSnap-m,gcIndOrigTr.gal[gal_w[i]]] ]
        tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
        r.accTime_gal[gal_w[i]] = time
        r.accHaloTvir_gal[gal_w[i]] = tvir
      endfor
      
      for i=0,count_gmem-1 do begin
        radii = [ prevRad.gmem[gmem_w[i]],gmem_pri[gmem_w[i]] ]
        time = interpol(times,radii,1.0) ; lerp time to r/rvir=1
        tvir = [ mt.hVirTemp[mt.maxSnap-m-1,gcIndOrigTr.gmem[gmem_w[i]]], $
                 mt.hVirTemp[mt.maxSnap-m,gcIndOrigTr.gmem[gmem_w[i]]] ]
        tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
        r.accTime_gmem[gmem_w[i]] = time
        r.accHaloTvir_gmem[gmem_w[i]] = tvir
      endfor
      
      ; if we are on the first snapshot, override accretion times with -1 to indicate always outside rvir
      if m eq mt.maxSnap then r.accTime_gal[gal_w] = -1
      if m eq mt.maxSnap then r.accTime_gmem[gmem_w] = -1
      
      ; update counters for the number of particles we have found the accretion times of
      accMask.gal[gal_w]   = 1B
      accMask.gmem[gmem_w] = 1B
      accCount.gal  += count_gal
      accCount.gmem += count_gmem
      
      ; store current radius of particles
      prevRad.gal  = gal_pri
      prevRad.gmem = gmem_pri
     
      ;radtemp.gal[mt.maxSnap-m,*] = gal_pri ;debug
      ;radtemp.gmem[mt.maxSnap-m,*] = gmem_pri ;debug
      
      prevTime = h.time
      
      ; free some memory for next load
      gal_w    = !NULL
      gmem_w   = !NULL
      gal_pri  = !NULL
      gmem_pri = !NULL
    endfor
    
    print,'found accretion times for ['+str(accCount.gal)+' of '+str(n_elements(galcat_gal_trids))+$
      '] galtr, ['+str(accCount.gmem)+' of '+str(n_elements(galcat_gmem_trids))+'] gmemtr'
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
  endif
  
  
  ; VELOCITY TRACERS case - will be similar to above since there could be multiple
  ; ---------------------
  if sP.trMCPerCell eq -1 then begin
  
  endif
  
end

; -----------------------------------------------------------------------------------------------------
; maxTemps(): find maximum temperature for gas particles in galaxy/group member catalogs at redshift
;             through the redshift range (redshift,zStart] where zStart is typically the start of 
;             the simulation
;
; NOTE: currently temps are only saved for gas in groups at the end of the interval (not all gas)
; saveRedshifts: if not set, then the history trace runs down to sP.snap and is saved
;                if set, save at each redshift as it is reached, running to the lowest
; -----------------------------------------------------------------------------------------------------

function maxTemps, sP=sP, zStart=zStart, saveRedshifts=saveRedshifts, $
                   loadByGas=loadByGas, loadAllTrGal=loadAllTrGal, loadAllTrGmem=loadAllTrGmem ; load opt

  forward_function cosmoTracerChildren, cosmoTracerVelParents
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; set minimum snapshot (maxmimum redshift)
  if not keyword_set(zStart) then zStart = 30.0
  minSnap = redshiftToSnapnum(zStart,sP=sP)

  ; set maximum snapshot (minimum redshift) (-1 do not include final snapshot in Tmax search)
  if not keyword_set(saveRedshifts) then $
    saveSnaps = [sP.snap - 1]
  if keyword_set(saveRedshifts) then $
    saveSnaps = redshiftToSnapnum(saveRedshifts,sP=sP) - 1

  if keyword_set(saveRedshifts) then stop ; DO NOT USE (in order to know which tracers to follow the
                                          ; properties of we have to first know the galaxy catalog of
                                          ; interest, but for multiple save times the only realistic
                                          ; approach would be to process all tracers and at each save
                                          ; time extract the subset in that galaxy catalog (TODO)
  maxSnap = max(saveSnaps)
  
  ; set saveFilename and check for existence
  if keyword_set(loadByGas) or keyword_set(loadAllTrGal) or keyword_set(loadAllTrGmem) then begin
    if ((keyword_set(loadAllTrGas) or keyword_set(loadAllTrGmem)) and sP.trMCPerCell eq 0) then $
      message,'Cannot load all tracers for SPH type.'
    
    saveTag = ''
    if sP.trMCPerCell eq -1 then saveTag = '.trVel'
    if sP.trMCPerCell gt 0 then  saveTag = '.trMC'
    if sP.trMCPerCell eq 0 then  saveTag = '.SPH'
    
    ; load the condensed results (one per gas particle in galcat)
    if keyword_set(loadByGas) then begin
      saveFilename = sP.derivPath+'maxtemp'+saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                     str(minSnap)+'-'+str(maxSnap)+'.sav'
      
      if not file_test(saveFilename) then message,'Error: Specified maxTemps not found!'
      restore, saveFilename
      return, r
    endif
    
    ; load the results for all tracers (galaxy)
    if keyword_set(loadAllTrGal) then begin
      saveFilename = sP.derivPath + 'maxtemp'+saveTag+'.gal.'+sP.savPrefix+str(sP.res)+'.'+$
                     str(minSnap)+'-'+str(maxSnap)+'.sav'
                     
      if not file_test(saveFilename) then message,'Error: Specified maxTemps gal not found!'
      restore, saveFilename
      return, rtr_gal
    endif
    
    ; load the results for all tracers (groupmem)
    if keyword_set(loadAllTrGmem) then begin
        saveFilename = sP.derivPath + 'maxtemp'+saveTag+'.gmem.'+sP.savPrefix+str(sP.res)+'.'+$
                        str(minSnap)+'-'+str(maxSnap)+'.sav'
                     
      if not file_test(saveFilename) then message,'Error: Specified maxTemps gmem not found!'
      restore, saveFilename
      return, rtr_gmem
    endif
  endif
  
  ; load galaxy/group member catalogs at zMin for gas ids to search for
  galcat = galaxyCat(sP=sP) ;sP.snap is still at zMin
  
  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    print,'Calculating new maxtemp using ( SPH Particles ) res = '+str(sP.res)+$
      ' in range ['+str(minSnap)+'-'+str(maxSnap)+'].'
      
    ; store the main arrays for all tracers as structures so we can write them directly
      r = {maxTemps_gal       : fltarr(n_elements(galcat.galaxyIDs))   ,$
           maxTemps_gmem      : fltarr(n_elements(galcat.groupmemIDs)) ,$
           maxTempTime_gal    : fltarr(n_elements(galcat.galaxyIDs))   ,$
           maxTempTime_gmem   : fltarr(n_elements(galcat.groupmemIDs)) ,$
           maxTempDisp_gal    : 0 ,$
           maxTempDisp_gmem   : 0 ,$ ; no dispersions for sph
           maxTemps_min_gal   : 0 ,$
           maxTemps_min_gmem  : 0 ,$ ; no minimum values for sph
           maxTemps_mean_gal  : 0 ,$
           maxTemps_mean_gmem : 0}   ; no mean values for sph
    
    for m=minSnap,maxSnap,1 do begin
      sP.snap = m
      print,m
      
      ; load gas ids and match to catalog
      ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs, need this if we want ids[ids_ind], 
      ; temp[ids_ind], etc to be in the same order as the group catalog id list    
      match,galcat.galaxyIDs,ids,galcat_ind,ids_gal_ind,count=countGal,/sort
      ids_gal_ind = ids_gal_ind[sort(galcat_ind)]
      
      match,galcat.groupmemIDs,ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
      ids_gmem_ind = ids_gmem_ind[sort(galcat_ind)]
      
      ids        = !NULL
      galcat_ind = !NULL
      
      ; load u,nelec to calculate temperatures
      u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='ne')
      
      u_gal  = u[ids_gal_ind]
      u_gmem = u[ids_gmem_ind]
      u = !NULL
      
      nelec_gal  = nelec[ids_gal_ind]
      nelec_gmem = nelec[ids_gmem_ind]
      nelec = !NULL
      
      temp_gal  = convertUtoTemp(u_gal,nelec_gal,/log)
      temp_gmem = convertUtoTemp(u_gal,nelec_gal,/log)
      
      ; load gas SFR and select off the effective EOS (SFR must be zero)
      sfr = loadSnapshotSubset(sP=sP,partType='gas',field='sfr')
      
      sfr_gal  = sfr[ids_gal_ind]
      sfr_gmem = sfr[ids_gmem_ind]
      
      sfr = !NULL
      
      ; replace existing values if current snapshot has higher temps (enforce off effective EOS)
      w1 = where(temp_gal gt r.maxTemps_gal and sfr_gal eq 0.0,count1)
      if (count1 gt 0) then begin
        r.maxTemps_gal[w1]    = temp_gal[w1]
        r.maxTempTime_gal[w1] = snapNumToRedshift(sP=sP,/time)
      endif
      
      ; replace existing values if current snapshot has higher temps (group members)
      w2 = where(temp_gmem gt r.maxTemps_gmem and sfr_gmem eq 0.0,count2)
      if (count2 gt 0) then begin
        r.maxTemps_gmem[w2]    = temp_gmem[w2]
        r.maxTempTime_gmem[w2] = snapNumToRedshift(sP=sP,/time)
      endif
      
      ; SAVE?
      if total(sP.snap eq saveSnaps) gt 0 then begin
        ; set savefilename
        w = where(sP.snap eq saveSnaps)
        saveFilename = sP.derivPath + 'maxtemp.SPH.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(saveSnaps[w[0]])+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
          continue
        endif
      
        save,r,filename=saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
      endif ; save
    endfor ;m
  endif
  
  ; MONTE CARLO TRACERS CASE - for each original gas cell, determine some statistics of its
  ; population of tracers and an estimate for the dispersion in those statistics
  ; ------------------------
  if sP.trMCPerCell gt 0 then begin
  
    print,'Calculating new maxtemp using ( TracerMC ) res = '+str(sP.res)+$
      ' in range ['+str(minSnap)+'-'+str(maxSnap)+'].'
      
    ; load gas ids
    gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')

    ; match galcat IDs to gas_ids
    match,galcat.galaxyIDs,gas_ids,galcat_ind,ids_gal_ind,count=countGal,/sort
    ids_gal = gas_ids[ids_gal_ind[sort(galcat_ind)]]
    
    match,galcat.groupmemIDs,gas_ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
    ids_gmem = gas_ids[ids_gmem_ind[sort(galcat_ind)]]
    
    gas_ids = !NULL
    
    ; locate tracer children (indices) of gas id subsets
    galcat_gal_trids  = cosmoTracerChildren(sP=sP, /getInds, gasIDs=ids_gal, child_counts=galcat_gal_cc)
    galcat_gmem_trids = cosmoTracerChildren(sP=sP, /getInds, gasIDs=ids_gmem, child_counts=galcat_gmem_cc)
    
    ; convert tracer children indices to tracer IDs at this zMin
    tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
    
    galcat_gal_trids  = tr_ids[galcat_gal_trids]
    galcat_gmem_trids = tr_ids[galcat_gmem_trids]

    tr_ids   = !NULL
    ids_gal  = !NULL
    ids_gmem = !NULL
    galcat   = !NULL ; not used past this point
    
    ; store the main arrays for all tracers as structures so we can write them directly
    rtr_gal  = { maxTemps      : fltarr(n_elements(galcat_gal_trids))  ,$
                 maxTempTime   : fltarr(n_elements(galcat_gal_trids))  ,$
                 child_counts  : galcat_gal_cc }
    galcat_gal_cc  = !NULL
           
    rtr_gmem = { maxTempTime  : fltarr(n_elements(galcat_gmem_trids)) ,$
                 maxTemps     : fltarr(n_elements(galcat_gmem_trids)) ,$
                 child_counts : galcat_gmem_cc }
    galcat_gmem_cc = !NULL

    for m=minSnap,maxSnap,1 do begin
      sP.snap = m
      print,m  
  
      ; load tracer ids and match to child ids from zMin
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs, need this if we want ids[ids_ind], 
      ; temp[ids_ind], etc to be in the same order as the group catalog id list    
      match,galcat_gal_trids,tr_ids,galcat_ind,trids_gal_ind,count=countGal,/sort
      trids_gal_ind = trids_gal_ind[sort(galcat_ind)]
      
      match,galcat_gmem_trids,tr_ids,galcat_ind,trids_gmem_ind,count=countGmem,/sort
      trids_gmem_ind = trids_gmem_ind[sort(galcat_ind)]
      
      tr_ids     = !NULL
      galcat_ind = !NULL
      
      ; load tracer maximum temperature and sub-snapshot time at this snapshot
      tr_maxtemp = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxtemp')
      
      temp_gal  = codeTempToLogK(tr_maxtemp[trids_gal_ind]) ; tracer output still in unit system
      temp_gmem = codeTempToLogK(tr_maxtemp[trids_gmem_ind])
      tr_maxtemp = !NULL
      
      tr_maxtemp_time = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxtemp_time')
      
      temp_time_gal  = tr_maxtemp_time[trids_gal_ind]
      temp_time_gmem = tr_maxtemp_time[trids_gmem_ind]
      tr_maxtemp_time = !NULL
      
      ; replace existing values if current snapshot has higher temps (galaxy members)
      w1 = where(temp_gal gt rtr_gal.maxTemps,count1)
      if (count1 gt 0) then begin
        rtr_gal.maxTemps[w1]    = temp_gal[w1]
        rtr_gal.maxTempTime[w1] = temp_time_gal[w1] ; sub-snapshot timing
      endif
      
      ; replace existing values if current snapshot has higher temps (group members)
      w2 = where(temp_gmem gt rtr_gmem.maxTemps,count2)
      if (count2 gt 0) then begin
        rtr_gmem.maxTemps[w2]    = temp_gmem[w2]
        rtr_gmem.maxTempTime[w2] = temp_time_gmem[w2] ; sub-snapshot timing
      endif
      
      temp_gal  = !NULL
      temp_gmem = !NULL
      temp_time_gal  = !NULL
      temp_time_gmem = !NULL
      
      ; SAVE?
      if total(sP.snap eq saveSnaps) gt 0 then begin
        ; (1) full tracer information (galaxy members) - set savefilename
        w = where(sP.snap eq saveSnaps)
        saveFilename = sP.derivPath + 'maxtemp.trMC.gal.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(saveSnaps[w[0]])+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
        endif else begin
          ; output all gal/gmem tracer values at this point and the child counts so they can be 
          ; matched back up to their parents later
          save,rtr_gal,filename=saveFilename
          print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        endelse
        
        ; (2) full tracer information (group members) - set savefilename
        w = where(sP.snap eq saveSnaps)
        saveFilename = sP.derivPath + 'maxtemp.trMC.gmem.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(saveSnaps[w[0]])+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
        endif else begin
          ; output all gal/gmem tracer values at this point and the child counts so they can be 
          ; matched back up to their parents later
          save,rtr_gmem,filename=saveFilename
          print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        endelse
        
        ; (3) values condensed to gas parents - set savefilename
        w = where(sP.snap eq saveSnaps)
        saveFilename = sP.derivPath + 'maxtemp.trMC.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(saveSnaps[w[0]])+'.sav'
                       
        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
          continue
        endif
        
        ; prepare output
        offset = 0L
        
        r = {maxTemps_gal       : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; max(maxTemps)
             maxTempTime_gal    : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; time of max(maxTemps)
             maxTempDisp_gal    : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; stddev(maxTemps)
             maxTemps_gmem      : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTempTime_gmem   : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTempDisp_gmem   : fltarr(n_elements(rtr_gmem.child_counts)) ,$
                                                                      $
             maxTemps_min_gal   : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; min(maxTemps)
             maxTemps_mean_gal  : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; mean(maxTemps)
             maxTemps_min_gmem  : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTemps_mean_gmem : fltarr(n_elements(rtr_gmem.child_counts)) ,$
                                                                      $
             maxTempTime_min_gal   : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; time of min(maxTemps)
             maxTempTime_mean_gal  : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; mean(times maxTemps)
             maxTempTime_min_gmem  : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTempTime_mean_gmem : fltarr(n_elements(rtr_gmem.child_counts))  }
        
        for i=0,n_elements(rtr_gal.child_counts)-1 do begin
          if rtr_gal.child_counts[i] gt 0 then begin
            ; for each gas cell, collect its child tracers
            childInds = lindgen(rtr_gal.child_counts[i]) + offset
            
            locMaxTemps    = rtr_gal.maxTemps[childInds]
            locMaxTempTime = rtr_gal.maxTempTime[childInds] ; this is really time (scale factor) for tracers
            
            ; save the maximum value reached by any child tracer and the stddev of the population maxima
            r.maxTemps_gal[i]    = max(locMaxTemps,indmax)
            r.maxTempTime_gal[i] = locMaxTempTime[indmax]
            r.maxTempDisp_gal[i] = stddev(locMaxTemps)
            
            ; save minimum and means (and associated times)
            r.maxTemps_min_gal[i]     = min(locMaxTemps,indmin)
            r.maxTempTime_min_gal[i]  = locMaxTempTime[indmin]
            r.maxTemps_mean_gal[i]    = mean(locMaxTemps)
            r.maxTempTime_mean_gal[i] = mean(locMaxTempTime)
            
            offset += rtr_gal.child_counts[i]
          endif
        endfor
        
        offset = 0L
        
        for i=0,n_elements(rtr_gmem.child_counts)-1 do begin
          if rtr_gmem.child_counts[i] gt 0 then begin
            ; for each gas cell, collect its child tracers
            childInds = lindgen(rtr_gmem.child_counts[i]) + offset
            
            locMaxTemps    = rtr_gmem.maxTemps[childInds]
            locMaxTempTime = rtr_gmem.maxTempTime[childInds]
            
            ; save the maximum value reached by any child tracer and the stddev of the population maxima
            r.maxTemps_gmem[i]    = max(locMaxTemps,indmax)
            r.maxTempTime_gmem[i] = locMaxTempTime[indmax]
            r.maxTempDisp_gmem[i] = stddev(locMaxTemps)
            
            ; save minimum and means (and associated times)
            r.maxTemps_min_gmem[i]     = min(locMaxTemps,indmin)
            r.maxTempTime_min_gmem[i]  = locMaxTempTime[indmin]
            r.maxTemps_mean_gmem[i]    = mean(locMaxTemps)
            r.maxTempTime_mean_gmem[i] = mean(locMaxTempTime)
            
            offset += rtr_gmem.child_counts[i]
          endif
        endfor
        
        ; save for future lookups        
        save,r,filename=saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        r = !NULL
      endif ;save
    endfor ;m  
  endif
  
  ; VELOCITY TRACERS case - will be similar to above since there could be multiple
  ; ---------------------
  if sP.trMCPerCell eq -1 then begin
  
    print,'Calculating new maxtemp using ( TracerVEL ) res = '+str(sP.res)+$
      ' in range ['+str(minSnap)+'-'+str(maxSnap)+'].'

    ; load gas ids
    gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')

    ; match galcat IDs to gas_ids
    match,galcat.galaxyIDs,gas_ids,galcat_ind,ids_gal_ind,count=countGal,/sort
    inds_gal = ids_gal_ind[sort(galcat_ind)]
    
    match,galcat.groupmemIDs,gas_ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
    inds_gmem = ids_gmem_ind[sort(galcat_ind)]

    gas_ids = !NULL
    
    ; locate tracer children (indices) of gas id subsets
    galcat_gal_trids  = cosmoTracerVelChildren(sP=sP,/getInds,gasInds=inds_gal,child_counts=galcat_gal_cc)
    galcat_gmem_trids = cosmoTracerVelChildren(sP=sP,/getInds,gasInds=inds_gmem,child_counts=galcat_gmem_cc)
    
    ; convert tracer children indices to tracer IDs at this zMin
    tr_ids = loadSnapshotSubset(sP=sP,partType='tracerVel',field='ids')
    
    galcat_gal_trids  = tr_ids[galcat_gal_trids]
    galcat_gmem_trids = tr_ids[galcat_gmem_trids]

    tr_ids   = !NULL
    inds_gal  = !NULL
    inds_gmem = !NULL
    galcat   = !NULL ; not used past this point

    ; store the main arrays for all tracers as structures so we can write them directly
    rtr_gal  = { maxTemps      : fltarr(n_elements(galcat_gal_trids))  ,$
                 maxTempTime   : fltarr(n_elements(galcat_gal_trids))  ,$
                 child_counts  : galcat_gal_cc }
    galcat_gal_cc  = !NULL
           
    rtr_gmem = { maxTempTime  : fltarr(n_elements(galcat_gmem_trids)) ,$
                 maxTemps     : fltarr(n_elements(galcat_gmem_trids)) ,$
                 child_counts : galcat_gmem_cc }
    galcat_gmem_cc = !NULL
  
    for m=minSnap,maxSnap,1 do begin
      sP.snap = m
      print,m  
  
      ; load tracer ids and match to child ids from zMin
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerVel',field='ids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs, need this if we want ids[ids_ind], 
      ; temp[ids_ind], etc to be in the same order as the group catalog id list    
      match,galcat_gal_trids,tr_ids,galcat_ind,trids_gal_ind,count=countGal,/sort
      trids_gal_ind = trids_gal_ind[sort(galcat_ind)]
      
      match,galcat_gmem_trids,tr_ids,galcat_ind,trids_gmem_ind,count=countGmem,/sort
      trids_gmem_ind = trids_gmem_ind[sort(galcat_ind)]
      
      tr_ids     = !NULL
      galcat_ind = !NULL
      
      ; load tracer maximum temperature and sub-snapshot time at this snapshot
      tr_maxtemp = loadSnapshotSubset(sP=sP,partType='tracerVel',field='tracer_maxtemp')
      
      temp_gal  = codeTempToLogK(tr_maxtemp[trids_gal_ind]) ; tracer output still in unit system
      temp_gmem = codeTempToLogK(tr_maxtemp[trids_gmem_ind])
      tr_maxtemp = !NULL
      
      tr_maxtemp_time = loadSnapshotSubset(sP=sP,partType='tracerVel',field='tracer_maxtemp_time')
      
      temp_time_gal  = tr_maxtemp_time[trids_gal_ind]
      temp_time_gmem = tr_maxtemp_time[trids_gmem_ind]
      tr_maxtemp_time = !NULL
      
      ; replace existing values if current snapshot has higher temps (galaxy members)
      w1 = where(temp_gal gt rtr_gal.maxTemps,count1)
      if (count1 gt 0) then begin
        rtr_gal.maxTemps[w1]    = temp_gal[w1]
        rtr_gal.maxTempTime[w1] = temp_time_gal[w1] ; sub-snapshot timing
      endif
      
      ; replace existing values if current snapshot has higher temps (group members)
      w2 = where(temp_gmem gt rtr_gmem.maxTemps,count2)
      if (count2 gt 0) then begin
        rtr_gmem.maxTemps[w2]    = temp_gmem[w2]
        rtr_gmem.maxTempTime[w2] = temp_time_gmem[w2] ; sub-snapshot timing
      endif
      
      temp_gal  = !NULL
      temp_gmem = !NULL
      temp_time_gal  = !NULL
      temp_time_gmem = !NULL
     
      ; SAVE?
      if total(sP.snap eq saveSnaps) gt 0 then begin
        ; (1) full tracer information (galaxy members) - set savefilename
        w = where(sP.snap eq saveSnaps)
        saveFilename = sP.derivPath + 'maxtemp.trVel.gal.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(saveSnaps[w[0]])+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
        endif else begin
          ; output all gal/gmem tracer values at this point and the child counts so they can be 
          ; matched back up to their parents later
          save,rtr_gal,filename=saveFilename
          print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        endelse
        
        ; (2) full tracer information (group members) - set savefilename
        w = where(sP.snap eq saveSnaps)
        saveFilename = sP.derivPath + 'maxtemp.trVel.gmem.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(saveSnaps[w[0]])+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
        endif else begin
          ; output all gal/gmem tracer values at this point and the child counts so they can be 
          ; matched back up to their parents later
          save,rtr_gmem,filename=saveFilename
          print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        endelse
        
        ; (3) values condensed to gas parents - set savefilename
        w = where(sP.snap eq saveSnaps)
        saveFilename = sP.derivPath + 'maxtemp.trVel.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(saveSnaps[w[0]])+'.sav'

        if file_test(saveFilename) then begin
          print,'WARNING: ['+saveFilename+'] exists, skipping write!'
          continue
        endif
        
        ; prepare output
        offset = 0L
        
        r = {maxTemps_gal       : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; max(maxTemps)
             maxTempTime_gal    : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; time of max(maxTemps)
             maxTempDisp_gal    : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; stddev(maxTemps)
             maxTemps_gmem      : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTempTime_gmem   : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTempDisp_gmem   : fltarr(n_elements(rtr_gmem.child_counts)) ,$
                                                                      $
             maxTemps_min_gal   : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; min(maxTemps)
             maxTemps_mean_gal  : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; mean(maxTemps)
             maxTemps_min_gmem  : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTemps_mean_gmem : fltarr(n_elements(rtr_gmem.child_counts)) ,$
                                                                      $
             maxTempTime_min_gal   : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; time of min(maxTemps)
             maxTempTime_mean_gal  : fltarr(n_elements(rtr_gal.child_counts))  ,$ ; mean(times maxTemps)
             maxTempTime_min_gmem  : fltarr(n_elements(rtr_gmem.child_counts)) ,$
             maxTempTime_mean_gmem : fltarr(n_elements(rtr_gmem.child_counts))  }
        
        for i=0,n_elements(rtr_gal.child_counts)-1 do begin
          if rtr_gal.child_counts[i] gt 0 then begin
            ; for each gas cell, collect its child tracers
            childInds = lindgen(rtr_gal.child_counts[i]) + offset
            
            locMaxTemps    = rtr_gal.maxTemps[childInds]
            locMaxTempTime = rtr_gal.maxTempTime[childInds] ; this is really time (scale factor) for tracers
            
            ; save the maximum value reached by any child tracer and the stddev of the population maxima
            r.maxTemps_gal[i]    = max(locMaxTemps,indmax)
            r.maxTempTime_gal[i] = locMaxTempTime[indmax]
            r.maxTempDisp_gal[i] = stddev(locMaxTemps)
            
            ; save minimum and means (and associated times)
            r.maxTemps_min_gal[i]     = min(locMaxTemps,indmin)
            r.maxTempTime_min_gal[i]  = locMaxTempTime[indmin]
            r.maxTemps_mean_gal[i]    = mean(locMaxTemps)
            r.maxTempTime_mean_gal[i] = mean(locMaxTempTime)
            
            offset += rtr_gal.child_counts[i]
          endif
        endfor
        
        offset = 0L
        
        for i=0,n_elements(rtr_gmem.child_counts)-1 do begin
          if rtr_gmem.child_counts[i] gt 0 then begin
            ; for each gas cell, collect its child tracers
            childInds = lindgen(rtr_gmem.child_counts[i]) + offset
            
            locMaxTemps    = rtr_gmem.maxTemps[childInds]
            locMaxTempTime = rtr_gmem.maxTempTime[childInds]
            
            ; save the maximum value reached by any child tracer and the stddev of the population maxima
            r.maxTemps_gmem[i]    = max(locMaxTemps,indmax)
            r.maxTempTime_gmem[i] = locMaxTempTime[indmax]
            r.maxTempDisp_gmem[i] = stddev(locMaxTemps)
            
            ; save minimum and means (and associated times)
            r.maxTemps_min_gmem[i]     = min(locMaxTemps,indmin)
            r.maxTempTime_min_gmem[i]  = locMaxTempTime[indmin]
            r.maxTemps_mean_gmem[i]    = mean(locMaxTemps)
            r.maxTempTime_mean_gmem[i] = mean(locMaxTempTime)
            
            offset += rtr_gmem.child_counts[i]
          endif
        endfor
        
        ; save for future lookups        
        save,r,filename=saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        r = !NULL
      endif ;save
    endfor ;m
  endif

end
