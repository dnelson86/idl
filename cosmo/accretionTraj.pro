; accretionTraj.pro
; particle trajectories with time for visualization
; dnelson apr.2012

; -----------------------------------------------------------------------------------------------------
; accretionTraj(): for each gas particle/tracer, starting at some redshift, track backwards in time
;                  with respect to the tracked parent halos (using mergerTree) and save the relative
;                  position and temperature (for visualization)
; -----------------------------------------------------------------------------------------------------

function accretionTraj, sP=sP, getVel=getVel

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
  
  saveFilenameVel = sP.derivPath + 'accTrajVel'+saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                    str(mt.maxSnap)+'-'+str(mt.minSnap)+'.sav'   
  
  if keyword_set(getVel) then begin            
    if file_test(saveFilenameVel) then begin
      restore, saveFilenameVel
      return, rvel
    endif
    message,'Error: Velocity save not found.'
  endif  
  
  saveFilename = sP.derivPath + 'accTraj'+saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(mt.maxSnap)+'-'+str(mt.minSnap)+'.sav'          
  
  if file_test(saveFilename) and ~keyword_set(dovel) then begin
    restore, saveFilename
    return, r
  endif

  ; load galaxy/group member catalogs at zMin for gas ids to search for
  galcat = galaxyCat(sP=sP)
  
  nSnapsTot = mt.maxSnap - mt.minSnap + 1 ; how many values to store in time dimension

  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    print,'Calculating new accretion trajectories using ( SPH Particles ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'

    ; store the main arrays as a structure so we can write them directly
    r = {relPos_gal    : fltarr(nSnapsTot,3,n_elements(mt.galcatSub.gal))    ,$
         relPos_gmem   : fltarr(nSnapsTot,3,n_elements(mt.galcatSub.gmem))   ,$
         curTemp_gal   : fltarr(nSnapsTot,n_elements(mt.galcatSub.gal))      ,$
         curTemp_gmem  : fltarr(nSnapsTot,n_elements(mt.galcatSub.gmem))      }
         
    ; save the velocities for hermite (known derivative) interpolation, separately for smaller files
    rvel = {vel_gal    : fltarr(nSnapsTot,3,n_elements(mt.galcatSub.gal))    ,$
            vel_gmem   : fltarr(nSnapsTot,3,n_elements(mt.galcatSub.gmem))   }
    
    for m=mt.maxSnap,mt.minSnap,-1 do begin
      sP.snap = m
      print,m
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
      
      ; load velocities and store
      vel   = loadSnapshotSubset(sP=sP,partType='gas',field='vel')
      
      rvel.vel_gal[mt.maxSnap-m,*,*]  = vel[*,ids_gal_ind]
      rvel.vel_gmem[mt.maxSnap-m,*,*] = vel[*,ids_gmem_ind]
      
      vel = !NULL
      
      ; load pos to calculate positions relative to halo centers
      pos   = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
      
      pos_gal  = pos[*,ids_gal_ind]
      pos_gmem = pos[*,ids_gmem_ind]
      
      pos = !NULL

      ; calculate current distance of gas particle from smoothed halo center position for galaxy members
      pos_gal[0,*] = reform(mt.hPos[mt.maxSnap-m,0,mt.gcIndOrig.gal]) - pos_gal[0,*]
      pos_gal[1,*] = reform(mt.hPos[mt.maxSnap-m,1,mt.gcIndOrig.gal]) - pos_gal[1,*]
      pos_gal[2,*] = reform(mt.hPos[mt.maxSnap-m,2,mt.gcIndOrig.gal]) - pos_gal[2,*]

      correctPeriodicDistVecs, pos_gal, sP=sP ; account for periodic distance function
      
      r.relPos_gal[mt.maxSnap-m,0,*] = pos_gal[0,*]
      r.relPos_gal[mt.maxSnap-m,1,*] = pos_gal[1,*]
      r.relPos_gal[mt.maxSnap-m,2,*] = pos_gal[2,*]
      pos_gal = !NULL
      
      ; for group members
      pos_gmem[0,*] = reform(mt.hPos[mt.maxSnap-m,0,mt.gcIndOrig.gmem]) - pos_gmem[0,*]
      pos_gmem[1,*] = reform(mt.hPos[mt.maxSnap-m,1,mt.gcIndOrig.gmem]) - pos_gmem[1,*]
      pos_gmem[2,*] = reform(mt.hPos[mt.maxSnap-m,2,mt.gcIndOrig.gmem]) - pos_gmem[2,*]

      correctPeriodicDistVecs, pos_gmem, sP=sP ; account for periodic distance function
      
      r.relPos_gmem[mt.maxSnap-m,0,*] = pos_gmem[0,*]
      r.relPos_gmem[mt.maxSnap-m,1,*] = pos_gmem[1,*]
      r.relPos_gmem[mt.maxSnap-m,2,*] = pos_gmem[2,*]
      pos_gmem = !NULL
      
      ; load gas u,nelec and calculate temperatures
      u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='ne')
      
      u_gal  = u[ids_gal_ind]
      u_gmem = u[ids_gmem_ind]
      u = !NULL
      
      nelec_gal  = nelec[ids_gal_ind]
      nelec_gmem = nelec[ids_gmem_ind]
      nelec = !NULL
      
      ; load gas SFR and select off the effective EOS (SFR must be zero)
      sfr = loadSnapshotSubset(sP=sP,partType='gas',field='sfr')
      
      sfr_gal  = sfr[ids_gal_ind]
      sfr_gmem = sfr[ids_gmem_ind]
      sfr = !NULL
      
      ; replace existing values if current snapshot has higher temps (enforce off effective EOS)
      w = where(sfr_gal eq 0.0,count)
      if count gt 0 then $
        r.curTemp_gal[mt.maxSnap-m,w]  = convertUtoTemp(u_gal[w],nelec_gal[w],/log)
        
      w = where(sfr_gmem eq 0.0,count)
      if count gt 0 then $
        r.curTemp_gmem[mt.maxSnap-m,w] = convertUtoTemp(u_gmem[w],nelec_gmem[w],/log)

      u_gal      = !NULL
      u_gmem     = !NULL
      nelec_gal  = !NULL
      nelec_gmem = !NULL
      sfr_gal    = !NULL
      sfr_gmem   = !NULL
    endfor
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
    save,rvel,filename=saveFilenameVel
    print,'Saved: '+strmid(saveFilenameVel,strlen(sp.derivPath))

  endif
  
  ; MONTE CARLO TRACERS CASE - for all gas cells, track back all child tracers
  ; ------------------------
  if sP.trMCPerCell gt 0 then begin
    message,'todo'
  endif
  
  
  ; VELOCITY TRACERS case - will be similar to above since there could be multiple
  ; ---------------------
  if sP.trMCPerCell eq -1 then begin
      print,'Calculating new accretion trajectories using ( TracerVEL ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
  
    ; load gas ids
    gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')

    ; match galcat IDs to gas_ids
    match,galcat.galaxyIDs[mt.galcatSub.gal],gas_ids,galcat_ind,ids_gal_ind,count=countGal,/sort
    inds_gal = ids_gal_ind[sort(galcat_ind)]
    
    match,galcat.groupmemIDs[mt.galcatSub.gmem],gas_ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
    inds_gmem = ids_gmem_ind[sort(galcat_ind)]

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
    
    ; store the main arrays as a structure so we can write them directly for all tracers
    r = {relPos_gal    : fltarr(nSnapsTot,3,n_elements(galcat_gal_trids))    ,$
         relPos_gmem   : fltarr(nSnapsTot,3,n_elements(galcat_gmem_trids))   ,$
         curTemp_gal   : fltarr(nSnapsTot,n_elements(galcat_gal_trids))      ,$
         curTemp_gmem  : fltarr(nSnapsTot,n_elements(galcat_gmem_trids))      }

    ; save the velocities for hermite (known derivative) interpolation, separately for smaller files
    rvel = {vel_gal    : fltarr(nSnapsTot,3,n_elements(galcat_gal_trids))    ,$
            vel_gmem   : fltarr(nSnapsTot,3,n_elements(galcat_gmem_trids))   }

    for m=mt.maxSnap,mt.minSnap,-1 do begin
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
      
      ; load velocities and store
      vel = loadSnapshotSubset(sP=sP,partType='tracerVel',field='vel')
      
      rvel.vel_gal[mt.maxSnap-m,*,*]  = vel[*,trids_gal_ind]
      rvel.vel_gmem[mt.maxSnap-m,*,*] = vel[*,trids_gmem_ind]
      
      vel = !NULL

      ; load tracer maxtemps for this timestep (best proxy for current temp without doing NN parent finds)
      tr_maxtemp = loadSnapshotSubset(sP=sP,partType='tracerVel',field='tracer_maxtemp')
      
      r.curTemp_gal[mt.maxSnap-m,*]  = codeTempToLogK(tr_maxtemp[trids_gal_ind])
      r.curTemp_gmem[mt.maxSnap-m,*] = codeTempToLogK(tr_maxtemp[trids_gmem_ind])

      tr_maxtemp = !NULL
      
      ; load pos to calculate tracer positions relative to halo centers
      pos = loadSnapshotSubset(sP=sP,partType='tracerVel',field='pos')
      
      tr_pos_gal  = pos[*,trids_gal_ind]
      tr_pos_gmem = pos[*,trids_gmem_ind]
      
      tr_parids_gal  = !NULL
      tr_parids_gmem = !NULL
      pos = !NULL

      ; calculate current distance of all tracers from smoothed halo center position for galaxy members
      tr_pos_gal[0,*] = reform(mt.hPos[mt.maxSnap-m,0,gcIndOrigTr.gal]) - tr_pos_gal[0,*]
      tr_pos_gal[1,*] = reform(mt.hPos[mt.maxSnap-m,1,gcIndOrigTr.gal]) - tr_pos_gal[1,*]
      tr_pos_gal[2,*] = reform(mt.hPos[mt.maxSnap-m,2,gcIndOrigTr.gal]) - tr_pos_gal[2,*]
      
      correctPeriodicDistVecs, tr_pos_gal, sP=sP ; account for periodic distance function
      
      r.relPos_gal[mt.maxSnap-m,0,*] = tr_pos_gal[0,*]
      r.relPos_gal[mt.maxSnap-m,1,*] = tr_pos_gal[1,*]
      r.relPos_gal[mt.maxSnap-m,2,*] = tr_pos_gal[2,*]
      tr_pos_gal = !NULL
      
      ; for group members
      tr_pos_gmem[0,*]  = reform(mt.hPos[mt.maxSnap-m,0,gcIndOrigTr.gmem]) - tr_pos_gmem[0,*]
      tr_pos_gmem[1,*]  = reform(mt.hPos[mt.maxSnap-m,1,gcIndOrigTr.gmem]) - tr_pos_gmem[1,*]
      tr_pos_gmem[2,*]  = reform(mt.hPos[mt.maxSnap-m,2,gcIndOrigTr.gmem]) - tr_pos_gmem[2,*]

      correctPeriodicDistVecs, tr_pos_gmem, sP=sP ; account for periodic distance function
      
      r.relPos_gmem[mt.maxSnap-m,0,*] = tr_pos_gmem[0,*]
      r.relPos_gmem[mt.maxSnap-m,1,*] = tr_pos_gmem[1,*]
      r.relPos_gmem[mt.maxSnap-m,2,*] = tr_pos_gmem[2,*]
      tr_pos_gmem = !NULL
    endfor
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
    save,rvel,filename=saveFilenameVel
    print,'Saved: '+strmid(saveFilenameVel,strlen(sp.derivPath))  
    
  endif

end

function smoothHaloPos, mt=mt, hInd=hInd, sP=sP

  ; apply smoothing to halo position over time
  smoothKer = 5     ; number of snapshots for boxcar window
  wrapTol   = 100.0 ; kpc 
  polyOrder = 3     ; order polynomial (2,3)
  muxCoeff  = [0.8,0.2] ; [smooth,polyfit] linear combination coefficients, should add to 1
  
  hPos = mt.hPos[*,*,hInd]
  
  ; uniform errors for cubic fit weighted towards the ends
  merrs = findgen(n_elements(mt.times))+0.1
  merrs[0:4] = [0.09,0.0925,0.095,0.0975,0.1]
  merrs[n_elements(merrs)-5:n_elements(merrs)-1] = [0.1,0.0975,0.095,0.0925,0.09]
  
  ; debug: plot fit
  start_PS,'test_smoothHaloPos.eps',xs=10,ys=8
    !p.multi = [0,1,3]
    for i=0,2 do begin
      cgPlot,mt.times,hPos[*,i],psym=4,xrange=minmax(mt.times),yrange=minmax(hPos[*,i]),/xs,/ys
      for j=2,3 do begin
        res = svdfit(mt.times,hPos[*,i],j,measure_errors=merrs,yfit=yfit)
        hPosXyz = muxCoeff[0]*smooth(hPos[*,i],smoothKer) + muxCoeff[1]*yfit
        cgPlot,mt.times,yfit,color=getColor(j),line=0,/overplot
        cgPlot,mt.times,hPosXyz,color=getColor(j),line=1,/overplot
      endfor
      legend,['m='+['2','3']],textcolor=getColor([2,3],/name),box=0,/top,/left
    endfor
  end_PS  
  
  ; extract single coordinate separately
  for j=0,2 do begin
    hPosXyz = hPos[*,j]
    
    ; detect periodic wrap in x,y,z
    if min(hPosXyz) lt wrapTol and max(hPosXyz) gt sP.boxSize-wrapTol then begin
      ; middle should be empty if both edges are occupied, otherwise this halo moved too much
      w = where(hPosXyz gt sP.boxSize*0.4 and hPosXyz lt sP.boxSize*0.6,count)
      if count ne 0 then message,'Error: Something odd in periodic wrap.'
      
      ; move lower points up by one boxSize, smooth, and move back
      w = where(hPosXyz lt sP.boxSize*0.5,count)

      hPosXyz[w] += sP.boxSize
      ; combine tophat smooth and n=3 cubic polyfit
      result = svdfit(mt.times,hPosXyz,polyOrder,measure_errors=merrs,yfit=yfit)
      hPosXyz = muxCoeff[0]*smooth(hPosXyz,smoothKer) + muxCoeff[1]*yfit
      hPosXyz[w] -= sP.boxSize
    endif else begin
      ; no wrap, just smooth
      result = svdfit(mt.times,hPosXyz,polyOrder,measure_errors=merrs,yfit=yfit)
      hPosXyz = muxCoeff[0]*smooth(hPosXyz,smoothKer) + muxCoeff[1]*yfit
    endelse
    
    ; save over original time sequence for this coordinate
    if max(abs(hPos[*,j]-hPosXyz)) gt wrapTol then $
      message,'Error: Big offset requested ('+str(j)+' max '+string(max(abs(hPos[*,j]-hPosXyz)))+').'
    hPos[*,j] = hPosXyz

  endfor ; j
  
  return, hPos

end

; -----------------------------------------------------------------------------------------------------
; accretionTrajSingle(): for each gas particle/tracer, starting at some redshift, track backwards in time
;                        with respect to the tracked parent halo (using mergerTree) and save the relative
;                        position and temperature (for visualization) of gas and positions of DM
; -----------------------------------------------------------------------------------------------------

function accretionTrajSingle, sP=sP, hInd=hInd

  forward_function cosmoTracerChildren, cosmoTracerVelParents
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; first, walk back through the merger tree and find primary subhalos with good parent histories
  mt = mergerTreeSubset(sP=sP,/verbose)

  ; set saveFilename and check for existence
  saveTag = '.h'+str(hInd)
  
  saveFilename = sP.derivPath + 'accTraj'+saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(mt.maxSnap)+'-'+str(mt.minSnap)+'.sav'          
  
  if file_test(saveFilename) and ~keyword_set(dovel) then begin
    restore, saveFilename
    return, r
  endif
  
  ; load group catalog at zMin to find particle ids to track
  gc = loadGroupCat(sP=sP,/readIDs)
  
  print,'halo ['+str(hInd)+'] mass: ',codeMassToLogMsun(gc.subgroupMass[mt.galcatIDList[hInd]])
  
  nSnapsTot = mt.maxSnap - mt.minSnap + 1 ; how many values to store in time dimension

  ; get smoothed halo position with time
  hPos = smoothHaloPos(mt=mt,hInd=hInd,sP=sP)

  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    print,'Calculating single accretion trajectories using ( SPH Particles ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'

    ; get particle ids
    gas_ids = gcPIDList(gc=gc,valGCids=[mt.galcatIDList[hInd]],partType='gas')
    dm_ids  = gcPIDList(gc=gc,valGCids=[mt.galcatIDList[hInd]],partType='dm')

    ; store the main arrays as a structure so we can write them directly
    r = {relPos_gas  : fltarr(nSnapsTot,3,n_elements(gas_ids)) ,$
         curTemp_gas : fltarr(nSnapsTot,n_elements(gas_ids))   ,$
         relPos_dm   : fltarr(nSnapsTot,3,n_elements(dm_ids))   }

    for m=mt.maxSnap,mt.minSnap,-1 do begin
      sP.snap = m
      print,m
      ; load gas ids and match to catalog
      h = loadSnapshotHeader(sP=sP)
      ids_local = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of original gas_ids 
      match,gas_ids,ids_local,orig_ids_ind,ids_gas_ind,count=countGas,/sort
      ids_gas_ind = ids_gas_ind[sort(orig_ids_ind)]
      
      ids_local    = !NULL
      orig_ids_ind = !NULL
      
      ; load pos to calculate positions relative to halo centers
      pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
      pos = pos[*,ids_gas_ind]

      ; calculate current distance of gas particle from smoothed halo center position for galaxy members
      pos[0,*] = hPos[mt.maxSnap-m,0] - pos[0,*]
      pos[1,*] = hPos[mt.maxSnap-m,1] - pos[1,*]
      pos[2,*] = hPos[mt.maxSnap-m,2] - pos[2,*]

      correctPeriodicDistVecs, pos, sP=sP ; account for periodic distance function
      
      r.relPos_gas[mt.maxSnap-m,0,*] = pos[0,*]
      r.relPos_gas[mt.maxSnap-m,1,*] = pos[1,*]
      r.relPos_gas[mt.maxSnap-m,2,*] = pos[2,*]
      pos = !NULL
      
      ; load gas u,nelec and calculate temperatures
      u = loadSnapshotSubset(sP=sP,partType='gas',field='u')
      u = u[ids_gas_ind]
      
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='ne')
      nelec = nelec[ids_gas_ind]
      
      ; load gas SFR and select off the effective EOS (SFR must be zero)
      sfr = loadSnapshotSubset(sP=sP,partType='gas',field='sfr')
      sfr  = sfr[ids_gas_ind]

      w = where(sfr eq 0.0,count)
      if count gt 0 then r.curTemp_gas[mt.maxSnap-m,w] = convertUtoTemp(u[w],nelec[w],/log)

      u     = !NULL
      nelec = !NULL
      sfr   = !NULL
      
      ; dark matter positions
      ids_local = loadSnapshotSubset(sP=sP,partType='dm',field='ids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of original gas_ids 
      match,dm_ids,ids_local,orig_ids_ind,ids_dm_ind,count=countDM,/sort
      ids_dm_ind = ids_dm_ind[sort(orig_ids_ind)]
      
      ids_local    = !NULL
      orig_ids_ind = !NULL
      
      ; load pos to calculate positions relative to halo centers
      pos = loadSnapshotSubset(sP=sP,partType='dm',field='pos')
      pos = pos[*,ids_dm_ind]

      ; calculate current distance of gas particle from smoothed halo center position for galaxy members
      pos[0,*] = hPos[mt.maxSnap-m,0] - pos[0,*]
      pos[1,*] = hPos[mt.maxSnap-m,1] - pos[1,*]
      pos[2,*] = hPos[mt.maxSnap-m,2] - pos[2,*]

      correctPeriodicDistVecs, pos, sP=sP ; account for periodic distance function
      
      r.relPos_dm[mt.maxSnap-m,0,*] = pos[0,*]
      r.relPos_dm[mt.maxSnap-m,1,*] = pos[1,*]
      r.relPos_dm[mt.maxSnap-m,2,*] = pos[2,*]
      pos = !NULL

    endfor
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))

  endif
  
end
