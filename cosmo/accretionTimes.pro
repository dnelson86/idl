; accretionTimes.pro
; gas accretion project - past radial history of gas elements (virial radius crossing)
; dnelson may.2012

; -----------------------------------------------------------------------------------------------------
; accretionTimes(): for each gas particle/tracer, starting at some redshift, track backwards in time
;                   with respect to the tracked parent halos (using mergerTree) and determine the
;                   time when the particle radius = the virial radius (and record the virial temp of
;                   the parent halo at that time)
; -----------------------------------------------------------------------------------------------------

function accretionTimes, sP=sP, restart=restart

  forward_function cosmoTracerChildren, cosmoTracerVelParents
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; config
  rVirFacs = [1.0,0.75,0.5,0.25,0.1,0.05,0.01]
  nVirFacs = n_elements(rVirFacs)

  ; first, walk back through the merger tree and find primary subhalos with good parent histories
  mt = mergerTreeSubset(sP=sP,/verbose)
  snapRange = [mt.maxSnap,mt.minSnap]
  
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
  
  resFilename = sP.derivPath + 'accTimes'+saveTag+'.restart.'+sP.savPrefix+str(sP.res)+'.'+$
                str(mt.maxSnap)+'-'+str(mt.minSnap)+'.sav'
  
  ; load galaxy/group member catalogs at zMin for gas ids to search for
  galcat = galaxyCat(sP=sP)
  
  accCount = { gal : ulonarr(nVirFacs), gmem : ulonarr(nVirFacs) } ; count of how many particles/tracers we tracked through r_vir
  prevTime = 0 ; scale factor at previous snapshot

  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    print,'Calculating new accretion time using ( SPH Particles ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
      
    if ~file_test(resFilename) then begin
      ; no restart
      ; store the r/rvir of each at the previous snapshot for interpolation, and as a monotonic check
      prevRad = { gal  : fltarr(n_elements(mt.galcatSub.gal)) ,$
                  gmem : fltarr(n_elements(mt.galcatSub.gmem)) }
      
      accMask = { gal  : bytarr(nVirFacs,n_elements(mt.galcatSub.gal)), $
                  gmem : bytarr(nVirFacs,n_elements(mt.galcatSub.gmem)) }
      
      ; store the main arrays as a structure so we can write them directly
      r = {accTime_gal       : fltarr(nVirFacs,n_elements(mt.galcatSub.gal))-1  ,$
           accTime_gmem      : fltarr(nVirFacs,n_elements(mt.galcatSub.gmem))-1 ,$
           accHaloTvir_gal   : fltarr(n_elements(mt.galcatSub.gal))             ,$
           accHaloTvir_gmem  : fltarr(n_elements(mt.galcatSub.gmem))            ,$
           rVirFacs          : rVirFacs                                          }
    endif else begin
      ; restart
      if ~keyword_set(restart) then message,'Error: Restart file exists but restart not requested.'
      restore,resFilename,/verbose
      snapRange[0] = m
    endelse
    
    for m=snapRange[0],snapRange[1],-1 do begin
      sP.snap = m
      
      ; save restart?
      if m mod 10 eq 0 and m lt snapRange[0] and keyword_set(restart) then begin
        print,' --- Writing restart! ---'
        save,prevRad,accMask,r,prevTime,m,filename=resFilename
        print,' --- Done! ---'
      endif
      
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
      
      if countGal ne n_elements(mt.galcatSub.gal) or countGmem ne n_elements(mt.galcatSub.gmem) then $
        message,'Error: Failed to locate all of galcat in gas_ids (overflow64?).'      
      
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
      
      ; loop over each target radius
      foreach rVirFac,rVirFacs,k do begin
        ; for particles who are still within r_vir, check if they have passed beyond
        gal_w  = where(gal_pri ge rVirFac and prevRad.gal lt rVirFac and accMask.gal[k,*] eq 0B,count_gal)
        gmem_w = where(gmem_pri ge rVirFac and prevRad.gmem lt rVirFac and accMask.gmem[k,*] eq 0B,count_gmem)
        
        print,' ['+string(m,format='(i3)')+'] ['+str(k)+'] r='+string(rVirFac,format='(f4.2)')+$
          ' accreted now counts '+string(count_gal,format='(i5)')+' ('+$
          string(float(count_gal)/n_elements(gal_pri)*100,format='(f4.1)')+'%) '+$
          string(count_gmem,format='(i5)')+' ('+$
          string(float(count_gmem)/n_elements(gmem_pri)*100,format='(f4.1)')+'%)'
        
        ; interpolate these (time,radii) to find time crossing the virial radius
        times = [prevTime,h.time]
        
        for i=0,count_gal-1 do begin
          radii = [ prevRad.gal[gal_w[i]],gal_pri[gal_w[i]] ]
          time = interpol(times,radii,rVirFac) ; lerp time to r/rvir=rVirFac
          r.accTime_gal[k,gal_w[i]] = time
          ; record the tvir of the halo at the rvir crossing time
          if k eq 0 then begin
            tvir = [ mt.hVirTemp[mt.maxSnap-m-1,mt.gcIndOrig.gal[gal_w[i]]], $
                     mt.hVirTemp[mt.maxSnap-m,mt.gcIndOrig.gal[gal_w[i]]] ]
            tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
            r.accHaloTvir_gal[gal_w[i]] = tvir
          endif
        endfor
        
        for i=0,count_gmem-1 do begin
          radii = [ prevRad.gmem[gmem_w[i]],gmem_pri[gmem_w[i]] ]
          time = interpol(times,radii,rVirFac) ; lerp time to r/rvir=rVirFac
          r.accTime_gmem[k,gmem_w[i]] = time
          if k eq 0 then begin
            tvir = [ mt.hVirTemp[mt.maxSnap-m-1,mt.gcIndOrig.gmem[gmem_w[i]]], $
                     mt.hVirTemp[mt.maxSnap-m,mt.gcIndOrig.gmem[gmem_w[i]]] ]
            tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
            r.accHaloTvir_gmem[gmem_w[i]] = tvir
          endif
        endfor
        
        ; if we are on the first snapshot, override accretion times with -1 to indicate always outside rad
        if m eq mt.maxSnap then begin
          r.accTime_gal[k,gal_w] = -1
          r.accTime_gmem[k,gmem_w] = -1
        endif else begin
          ; otherwise, update counters for the number of particles we have found the accretion times of
          accCount.gal[k]  += count_gal
          accCount.gmem[k] += count_gmem
        endelse
        
        ; update mask for particles we no longer search for
        accMask.gal[k,gal_w]   = 1B
        accMask.gmem[k,gmem_w] = 1B

      endforeach
      
      ; store current radius of particles
      prevRad.gal  = gal_pri
      prevRad.gmem = gmem_pri
      
      prevTime = h.time
      
      ; free some memory for next load
      gal_w    = !NULL
      gmem_w   = !NULL
      gal_pri  = !NULL
      gmem_pri = !NULL
    endfor
    
    foreach rVirFac,rVirFacs,k do $
    print,'['+str(k)+'] r='+string(rVirFacs[k],format='(f4.2)')+' found accretion times for ['+$
      str(accCount.gal[k])+' of '+str(n_elements(mt.galcatSub.gal))+$
      '] gal, ['+str(accCount.gmem[k])+' of '+str(n_elements(mt.galcatSub.gmem))+'] gmem'
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
  endif
  
  ; MONTE CARLO TRACERS CASE - for all gas cells, track back all child tracers
  ; ------------------------
  if sP.trMCPerCell gt 0 then begin
                  
    if ~file_test(resFilename) then begin
      ; no restart
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
      
      if countGal ne n_elements(mt.galcatSub.gal) or countGmem ne n_elements(mt.galcatSub.gmem) then $
        message,'Error: Failed to locate all of galcat in gas_ids (overflow64?).'    
      
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
      
      accMask = { gal  : bytarr(nVirFacs,n_elements(galcat_gal_trids)), $
                  gmem : bytarr(nVirFacs,n_elements(galcat_gmem_trids)) }
      
      ; store the main arrays as a structure so we can write them directly for all tracers
      r = {accTime_gal       : fltarr(nVirFacs,n_elements(galcat_gal_trids))-1  ,$
           accTime_gmem      : fltarr(nVirFacs,n_elements(galcat_gmem_trids))-1 ,$
           accHaloTvir_gal   : fltarr(n_elements(galcat_gal_trids))             ,$
           accHaloTvir_gmem  : fltarr(n_elements(galcat_gmem_trids))            ,$
           child_counts_gal  : galcat_gal_cc                                    ,$
           child_counts_gmem : galcat_gmem_cc                                   ,$
           rVirFacs          : rVirFacs                                          }
    endif else begin
      ; restart
      if ~keyword_set(restart) then message,'Error: Restart file exists but restart not requested.'
      restore,resFilename,/verbose
      snapRange[0] = m
    endelse

    for m=snapRange[0],snapRange[1],-1 do begin
      sP.snap = m
      
      ; save restart?
      if m mod 10 eq 0 and m lt snapRange[0] and keyword_set(restart) then begin
        print,' --- Writing restart! ---'
        save,prevRad,accMask,r,gcIndOrigTr,galcat_gal_trids,galcat_gmem_trids,prevTime,m,filename=resFilename
        print,' --- Done! ---'
      endif
      
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
      
      ; loop over each critical radius
      foreach rVirFac,rVirFacs,k do begin
        ; for particles who are still within r_vir, check if they have passed beyond
        gal_w  = where(gal_pri ge rVirFac and prevRad.gal lt rVirFac and accMask.gal[k,*] eq 0B,count_gal)
        gmem_w = where(gmem_pri ge rVirFac and prevRad.gmem lt rVirFac and accMask.gmem[k,*] eq 0B,count_gmem)
        
        print,' ['+string(m,format='(i3)')+'] ['+str(k)+'] r='+string(rVirFac,format='(f4.2)')+$
          ' accreted now counts '+string(count_gal,format='(i5)')+' ('+$
          string(float(count_gal)/n_elements(gal_pri)*100,format='(f4.1)')+'%) '+$
          string(count_gmem,format='(i5)')+' ('+$
          string(float(count_gmem)/n_elements(gmem_pri)*100,format='(f4.1)')+'%)'
        
        ; interpolate these (time,radii) to find time crossing the virial radius
        times = [prevTime,h.time]
        
        for i=0,count_gal-1 do begin
          radii = [ prevRad.gal[gal_w[i]],gal_pri[gal_w[i]] ]
          time = interpol(times,radii,rVirFac) ; lerp time to r/rvir=rVirFac
          r.accTime_gal[k,gal_w[i]] = time
          if k eq 0 then begin
            tvir = [ mt.hVirTemp[mt.maxSnap-m-1,gcIndOrigTr.gal[gal_w[i]]], $
                     mt.hVirTemp[mt.maxSnap-m,gcIndOrigTr.gal[gal_w[i]]] ]
            tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
            r.accHaloTvir_gal[gal_w[i]] = tvir
          endif
        endfor
        
        for i=0,count_gmem-1 do begin
          radii = [ prevRad.gmem[gmem_w[i]],gmem_pri[gmem_w[i]] ]
          time = interpol(times,radii,rVirFac) ; lerp time to r/rvir=rVirFac
          r.accTime_gmem[k,gmem_w[i]] = time
          if k eq 0 then begin
            tvir = [ mt.hVirTemp[mt.maxSnap-m-1,gcIndOrigTr.gmem[gmem_w[i]]], $
                     mt.hVirTemp[mt.maxSnap-m,gcIndOrigTr.gmem[gmem_w[i]]] ]
            tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
            r.accHaloTvir_gmem[gmem_w[i]] = tvir
          endif
        endfor
        
        ; if we are on the first snapshot, override accretion times with -1 to indicate always outside rad
        if m eq mt.maxSnap then begin
          r.accTime_gal[k,gal_w] = -1
          r.accTime_gmem[k,gmem_w] = -1
        endif else begin
          ; otherwise, update counters for the number of particles we have found the accretion times of
          accCount.gal[k]  += count_gal
          accCount.gmem[k] += count_gmem
        endelse
        
        ; update mask for particles we no longer search for
        accMask.gal[k,gal_w]   = 1B
        accMask.gmem[k,gmem_w] = 1B
      endforeach
      
      ; store current radius of particles
      prevRad.gal  = gal_pri
      prevRad.gmem = gmem_pri
     
      prevTime = h.time
      
      ; free some memory for next load
      gal_w    = !NULL
      gmem_w   = !NULL
      gal_pri  = !NULL
      gmem_pri = !NULL
    endfor
    
    foreach rVirFac,rVirFacs,k do $
    print,'['+str(k)+'] r='+string(rVirFacs[k],format='(f4.2)')+' found accretion times for ['+$
      str(accCount.gal[k])+' of '+str(n_elements(galcat_gal_trids))+$
      '] gal, ['+str(accCount.gmem[k])+' of '+str(n_elements(galcat_gmem_trids))+'] gmem'
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
  endif
  
  
  ; VELOCITY TRACERS case - will be similar to above since there could be multiple
  ; ---------------------
  if sP.trMCPerCell eq -1 then begin
      print,'Calculating new accretion time using ( TracerVEL ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
  
    ; load gas ids
    gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')

    ; match galcat IDs to gas_ids
    match,galcat.galaxyIDs[mt.galcatSub.gal],gas_ids,galcat_ind,ids_gal_ind,count=countGal,/sort
    inds_gal = ids_gal_ind[sort(galcat_ind)]
    
    match,galcat.groupmemIDs[mt.galcatSub.gmem],gas_ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
    inds_gmem = ids_gmem_ind[sort(galcat_ind)]

    if countGal ne n_elements(mt.galcatSub.gal) or countGmem ne n_elements(mt.galcatSub.gmem) then $
      message,'Error: Failed to locate all of galcat in gas_ids (overflow64?).'

    gas_ids = !NULL
    
    ; locate tracer children (indices) of gas id subsets
    galcat_gal_trids  = cosmoTracerVelChildren(sP=sP,/getInds,gasInds=inds_gal,child_counts=galcat_gal_cc)
    galcat_gmem_trids = cosmoTracerVelChildren(sP=sP,/getInds,gasInds=inds_gmem,child_counts=galcat_gmem_cc)
    
    ; convert tracer children indices to tracer IDs at this zMin
    tr_ids = loadSnapshotSubset(sP=sP,partType='tracerVel',field='ids')
    
    galcat_gal_trids  = tr_ids[galcat_gal_trids]
    galcat_gmem_trids = tr_ids[galcat_gmem_trids]

    tr_ids    = !NULL
    inds_gal  = !NULL
    inds_gmem = !NULL

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
    
    accMask = { gal  : bytarr(nVirFacs,n_elements(galcat_gal_trids)), $
                gmem : bytarr(nVirFacs,n_elements(galcat_gmem_trids)) }
                
    ; store the main arrays as a structure so we can write them directly for all tracers
    r = {accTime_gal       : fltarr(nVirFacs,n_elements(galcat_gal_trids))-1  ,$
         accTime_gmem      : fltarr(nVirFacs,n_elements(galcat_gmem_trids))-1 ,$
         accHaloTvir_gal   : fltarr(n_elements(galcat_gal_trids))             ,$
         accHaloTvir_gmem  : fltarr(n_elements(galcat_gmem_trids))            ,$
         child_counts_gal  : galcat_gal_cc                                    ,$
         child_counts_gmem : galcat_gmem_cc                                   ,$
         rVirFacs          : rVirFacs                                          }
            
    for m=snapRange[0],snapRange[1],-1 do begin
      sP.snap = m
      print,m
      ; load tracer ids and match to child ids from zMin
      h = loadSnapshotHeader(sP=sP)
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerVel',field='ids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs, need this if we want ids[ids_ind], 
      ; temp[ids_ind], etc to be in the same order as the group catalog id list    
      match,galcat_gal_trids,tr_ids,galcat_ind,trids_gal_ind,count=countGal,/sort
      trids_gal_ind = trids_gal_ind[sort(galcat_ind)]
      
      match,galcat_gmem_trids,tr_ids,galcat_ind,trids_gmem_ind,count=countGmem,/sort
      trids_gmem_ind = trids_gmem_ind[sort(galcat_ind)]
      
      tr_ids     = !NULL
      galcat_ind = !NULL
      
      ; load pos to calculate tracer positions relative to halo centers
      pos = loadSnapshotSubset(sP=sP,partType='tracerVel',field='pos')
      
      tr_pos_gal  = pos[*,trids_gal_ind]
      tr_pos_gmem = pos[*,trids_gmem_ind]

      pos = !NULL
      
      ; calculate current distance of gas particle from smoothed halo center position for galaxy members
      gal_pri  = periodicDists(reform(mt.hPos[mt.maxSnap-m,*,gcIndOrigTr.gal]),tr_pos_gal,sP=sP)
      gal_pri /= mt.hVirRad[mt.maxSnap-m,gcIndOrigTr.gal]
      
      ; for group members
      gmem_pri = periodicDists(reform(mt.hPos[mt.maxSnap-m,*,gcIndOrigTr.gmem]),tr_pos_gmem,sP=sP)
      gmem_pri /= mt.hVirRad[mt.maxSnap-m,gcIndOrigTr.gmem]
      
      tr_pos_gal  = !NULL
      tr_pos_gmem = !NULL
      
      ; loop over each critical radius
      foreach rVirFac,rVirFacs,k do begin
        ; for particles who are still within r_vir, check if they have passed beyond
        gal_w  = where(gal_pri ge rVirFac and prevRad.gal lt rVirFac and accMask.gal[k,*] eq 0B,count_gal)
        gmem_w = where(gmem_pri ge rVirFac and prevRad.gmem lt rVirFac and accMask.gmem[k,*] eq 0B,count_gmem)
        
        print,' ['+string(m,format='(i3)')+'] ['+str(k)+'] r='+string(rVirFac,format='(f4.2)')+$
          ' accreted now counts '+string(count_gal,format='(i5)')+' ('+$
          string(float(count_gal)/n_elements(gal_pri)*100,format='(f4.1)')+'%) '+$
          string(count_gmem,format='(i5)')+' ('+$
          string(float(count_gmem)/n_elements(gmem_pri)*100,format='(f4.1)')+'%)'
        
        ; interpolate these (time,radii) to find time crossing the virial radius
        times = [prevTime,h.time]
        
        for i=0,count_gal-1 do begin
          radii = [ prevRad.gal[gal_w[i]],gal_pri[gal_w[i]] ]
          time = interpol(times,radii,rVirFac) ; lerp time to r/rvir=rVirFac
          r.accTime_gal[k,gal_w[i]] = time
          if k eq 0 then begin
            tvir = [ mt.hVirTemp[mt.maxSnap-m-1,gcIndOrigTr.gal[gal_w[i]]], $
                     mt.hVirTemp[mt.maxSnap-m,gcIndOrigTr.gal[gal_w[i]]] ]
            tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
            r.accHaloTvir_gal[gal_w[i]] = tvir
          endif
        endfor
        
        for i=0,count_gmem-1 do begin
          radii = [ prevRad.gmem[gmem_w[i]],gmem_pri[gmem_w[i]] ]
          time = interpol(times,radii,rVirFac) ; lerp time to r/rvir=rVirFac
          r.accTime_gmem[k,gmem_w[i]] = time
          if k eq 0 then begin
            tvir = [ mt.hVirTemp[mt.maxSnap-m-1,gcIndOrigTr.gmem[gmem_w[i]]], $
                     mt.hVirTemp[mt.maxSnap-m,gcIndOrigTr.gmem[gmem_w[i]]] ]
            tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
            r.accHaloTvir_gmem[gmem_w[i]] = tvir
          endif
        endfor
        
        ; if we are on the first snapshot, override accretion times with -1 to indicate always outside rad
        if m eq mt.maxSnap then begin
          r.accTime_gal[k,gal_w] = -1
          r.accTime_gmem[k,gmem_w] = -1
        endif else begin
          ; otherwise, update counters for the number of particles we have found the accretion times of
          accCount.gal[k]  += count_gal
          accCount.gmem[k] += count_gmem
        endelse
        
        ; update mask for particles we no longer search for
        accMask.gal[k,gal_w]   = 1B
        accMask.gmem[k,gmem_w] = 1B
      endforeach
      
      ; store current radius of particles
      prevRad.gal  = gal_pri
      prevRad.gmem = gmem_pri
     
      prevTime = h.time
      
      ; free some memory for next load
      gal_w    = !NULL
      gmem_w   = !NULL
      gal_pri  = !NULL
      gmem_pri = !NULL
    endfor
              
    foreach rVirFac,rVirFacs,k do $
    print,'['+str(k)+'] r='+string(rVirFacs[k],format='(f4.2)')+' found accretion times for ['+$
      str(accCount.gal[k])+' of '+str(n_elements(galcat_gal_trids))+$
      '] gal, ['+str(accCount.gmem[k])+' of '+str(n_elements(galcat_gmem_trids))+'] gmem'
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
                
  endif
  
end
