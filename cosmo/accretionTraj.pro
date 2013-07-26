; accretionTraj.pro
; particle trajectories with time for visualization
; dnelson jul.2013

; smoothHaloPos(): smooth the center position or bulk velocity over time of a halo

function smoothHaloPos, mt=mt, hInd=hInd, hAll=hAll, sP=sP, vel=vel

  ; apply smoothing to halo position over time
  smoothKer = 5     ; number of snapshots for boxcar window
  wrapTol   = 200.0 ; kpc 
  polyOrder = 14    ; order polynomial (2,3)
  muxCoeff  = [0.0,1.0] ; [smooth,polyfit] linear combination coefficients, should add to 1
  
  ; do one halo or all halos?
  if ~keyword_set(hAll) then hInds = [hInd]
  if keyword_set(hAll) then hInds = lindgen(n_elements(mt.hMinSnap))
  
  ; if smoothing velocity, do directly and return
  if keyword_set(vel) then return,smooth(mt.hVel[*,*,hInds],[smoothKer,1,1])  
  
  ; form subset of positions
  hPos = mt.hPos[*,*,hInds]  
  
  hPos = reform(hPos,[n_elements(mt.times),3,n_elements(hInds)]) ; [nSnaps,3,nHalos]  

  ; fit each halo
  foreach hInd,hInds,k do begin
    ; uniform errors for cubic fit weighted towards the ends
    merrs = findgen(n_elements(mt.times))+0.1
    ;merrs[0:4] = [0.09,0.0925,0.095,0.0975,0.1]
    ;merrs[-5:-1] = [0.1,0.0975,0.095,0.0925,0.09]
    
    ; debug: plot fit
    ;start_PS,'test_smoothHaloPos.eps',xs=10,ys=8
    ;  !p.multi = [0,1,3]
        
    ; extract single coordinate separately
    for j=0,2 do begin
      hPosXyz = hPos[*,j,k]
     
      ; detect periodic wrap in x,y,z
      if min(hPosXyz) lt wrapTol and max(hPosXyz) gt sP.boxSize-wrapTol then begin
        ; middle should be empty if both edges are occupied, otherwise this halo moved too much
        w = where(hPosXyz gt sP.boxSize*0.4 and hPosXyz lt sP.boxSize*0.6,count)
        if count ne 0 then message,'Error: Something odd in periodic wrap.'
        
        ; move lower points up by one boxSize, smooth, and move back
        w = where(hPosXyz lt sP.boxSize*0.5,count)
  
        hPosXyz[w] += sP.boxSize
  
        ; combine tophat smooth and n=3 cubic polyfit
        result = svdfit(mt.times,hPosXyz,polyOrder,yfit=yfit,measure_errors=merrs,status=status)
        ;  cgPlot,mt.times,hPosXyz,psym=4,xrange=minmax(mt.times),yrange=minmax(hPosXyz)+[-100,100],/xs,/ys
        hPosXyz = muxCoeff[0]*smooth(hPosXyz,smoothKer) + muxCoeff[1]*yfit
        ;  cgPlot,mt.times,yfit,color=cgColor('red'),line=0,/overplot
        ;  cgPlot,mt.times,hPosXyz,color=cgColor('blue'),line=1,/overplot
        hPosXyz[w] -= sP.boxSize
      endif else begin
        ; no wrap, just smooth
        result = svdfit(mt.times,hPosXyz,polyOrder,yfit=yfit,measure_errors=merrs,status=status)
        ;  cgPlot,mt.times,hPosXyz,psym=4,xrange=minmax(mt.times),yrange=minmax(hPosXyz)+[-100,100],/xs,/ys
        hPosXyz = muxCoeff[0]*smooth(hPosXyz,smoothKer) + muxCoeff[1]*yfit
        ;  cgPlot,mt.times,yfit,color=cgColor('red'),line=0,/overplot
        ;  cgPlot,mt.times,hPosXyz,color=cgColor('blue'),line=1,/overplot
      endelse
      
      ; save over original time sequence for this coordinate
      if max(abs(hPos[*,j,k]-hPosXyz)) gt wrapTol then $
        print,'Error: Big offset requested ('+str(j)+' max '+string(max(abs(hPos[*,j]-hPosXyz)))+').'
        
      hPos[*,j,k] = hPosXyz
  
    endfor ; j
    
    ;end_PS
  endforeach ;hInds
  
  if ~keyword_set(hAll) then hPos = reform(hPos,[n_elements(mt.times),3]) ; [nSnaps,3] if only one halo
  return, hPos
end

; -----------------------------------------------------------------------------------------------------
; accretionTraj(): for each gas particle/tracer, starting at some redshift, track backwards in time
;                  with respect to the tracked parent halos (using mergerTree) and save the relative
;                  position and temperature (for visualization)
; -----------------------------------------------------------------------------------------------------

function accretionTraj, sP=sP

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
    
  saveFilename = sP.derivPath + 'accTraj'+saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(mt.maxSnap)+'-'+str(mt.minSnap)+'.sav'          
  
  if file_test(saveFilename) then begin
    restore, saveFilename
    return, r
  endif

  ; load galaxy/group member catalogs at zMin for gas ids to search for
  galcat = galaxyCat(sP=sP)
  
  nSnapsTot = mt.maxSnap - mt.minSnap + 1 ; how many values to store in time dimension

  ; get smoothed halo position and bulk velocity with time
  hPos = smoothHaloPos(mt=mt,/hAll,sP=sP)
  hVel = smoothHaloPos(mt=mt,/hAll,sP=sP,/vel)

  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    print,'Calculating new accretion trajectories using ( SPH Particles ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'

    ; store the main arrays as a structure so we can write them directly
    r = {relPos   : fltarr(nSnapsTot,3,n_elements(galcat.ids)) ,$
         relVel   : fltarr(nSnapsTot,3,n_elements(galcat.ids)) ,$
         curTemp  : fltarr(nSnapsTot,n_elements(galcat.ids))    }
         
    for m=mt.maxSnap,mt.minSnap,-1 do begin
      sP.snap = m
      print,m
      
      message,'TODO: add star,bh parent types, skip temp for them'
      
      ; load gas ids and match to catalog
      h = loadSnapshotHeader(sP=sP)
      ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs   
      calcMatch,galcat.ids,ids,galcat_ind,ids_ind,count=countMatch
      ids_ind = ids_ind[sort(galcat_ind)]
      
      ids        = !NULL
      galcat_ind = !NULL
      
      ; load velocities and store relative to smoothed halo CM velocity
      vel   = loadSnapshotSubset(sP=sP,partType='gas',field='vel')
      
      r.relVel[mt.maxSnap-m,0,*] = vel[0,ids_ind] - hVel[mt.maxSnap-m,0,mt.gcIndOrig] 
      r.relVel[mt.maxSnap-m,1,*] = vel[1,ids_ind] - hVel[mt.maxSnap-m,1,mt.gcIndOrig] 
      r.relVel[mt.maxSnap-m,2,*] = vel[2,ids_ind] - hVel[mt.maxSnap-m,2,mt.gcIndOrig] 
      
      vel = !NULL
      
      ; load pos to calculate positions relative to halo centers
      pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
      pos = pos[*,ids_ind]
      
      ; calculate current distance of gas particle from smoothed halo center position for galaxy members
      pos[0,*] = hPos[mt.maxSnap-m,0,mt.gcIndOrig] - pos[0,*]
      pos[1,*] = hPos[mt.maxSnap-m,1,mt.gcIndOrig] - pos[1,*]
      pos[2,*] = hPos[mt.maxSnap-m,2,mt.gcIndOrig] - pos[2,*]

      correctPeriodicDistVecs, pos, sP=sP ; account for periodic distance function
      
      r.relPos[mt.maxSnap-m,0,*] = pos[0,*]
      r.relPos[mt.maxSnap-m,1,*] = pos[1,*]
      r.relPos[mt.maxSnap-m,2,*] = pos[2,*]
      pos = !NULL
      
      ; load gas u,nelec and calculate temperatures
      u = loadSnapshotSubset(sP=sP,partType='gas',field='u')
      u = u[ids_ind]
      
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='ne')
      nelec = nelec[ids_ind]
      
      ; load gas SFR and select off the effective EOS (SFR must be zero)
      sfr = loadSnapshotSubset(sP=sP,partType='gas',field='sfr')
      sfr = sfr[ids_ind]
      
      ; replace existing values if current snapshot has higher temps (enforce off effective EOS)
      w = where(sfr eq 0.0,count)
      if count gt 0 then $
        r.curTemp[mt.maxSnap-m,w]  = convertUtoTemp(u[w],nelec[w],/log)
        
      u     = !NULL
      nelec = !NULL
      sfr   = !NULL
    endfor
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
  endif
  
  ; MONTE CARLO TRACERS CASE - for all gas cells, track back all child tracers
  ; ------------------------
  if sP.trMCPerCell gt 0 then begin
    print,'Calculating new accretion trajectories using ( TracerMC ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
      
    galcat_trids  = mt.trMC_ids

    ; store the main arrays as a structure so we can write them directly for all tracers
    r = {relPos  : fltarr(nSnapsTot,3,n_elements(galcat.trMC_ids))   ,$
         relVel  : fltarr(nSnapsTot,3,n_elements(galcat.trMC_ids))   ,$
         curTemp : fltarr(nSnapsTot,n_elements(galcat.trMC_ids))      }
 
    for m=mt.maxSnap,mt.minSnap,-1 do begin
      sP.snap = m
      print,m
      ; load tracer ids and match to child ids from zMin
      h = loadSnapshotHeader(sP=sP)
      
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
      
      idIndexMap = getIDIndexMap(tr_ids,minid=minid)
      trids_ind  = idIndexMap[ galcat.trMC_ids - minid ]
        
      idIndexMap = !NULL
      tr_ids     = !NULL
      
      message,'TODO: handle other parent types (stars,bh) or restrict input tracer search to gas'
      
      ; load tracer parents to match to gas
      tr_parids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentid')
      tr_parids = tr_parids[trids_ind]
      
      ; load gas IDs and convert tracer parent IDs -> indices
      gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
      gasIDMap = getIDIndexMap(gas_ids,minid=minid)
      gas_ids = !NULL
      
      tr_parids  = gasIDMap[tr_parids-minid]  ; convert ID->index
      gasIDMap = !NULL
      
      ; load velocities and store relative to smoothed halo CM velocity
      vel   = loadSnapshotSubset(sP=sP,partType='gas',field='vel')
      
      r.relVel[mt.maxSnap-m,0,*] = vel[0,tr_parids] - hVel[mt.maxSnap-m,0,gcIndOrigTrMC] 
      r.relVel[mt.maxSnap-m,1,*] = vel[1,tr_parids] - hVel[mt.maxSnap-m,1,gcIndOrigTrMC] 
      r.relVel[mt.maxSnap-m,2,*] = vel[2,tr_parids] - hVel[mt.maxSnap-m,2,gcIndOrigTrMC] 
      
      vel = !NULL

      ; load gas u,nelec and calculate temperatures
      u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='ne')
      temp  = convertUtoTemp(u[tr_parids],nelec[tr_parids],/log)
      u     = !NULL
      nelec = !NULL
      
      ; load gas SFR and select off the effective EOS (SFR must be zero)
      sfr = loadSnapshotSubset(sP=sP,partType='gas',field='sfr')
      sfr = sfr[tr_parids]
      
      ; replace existing values if current snapshot has higher temps (enforce off effective EOS)
      w = where(sfr eq 0.0,count)
      if count gt 0 then r.curTemp[mt.maxSnap-m,w]  = temp[w]
      temp = !NULL
      sfr  = !NULL
      
      ; load pos to calculate tracer positions relative to halo centers
      pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
      
      tr_pos  = pos[*,tr_parids]
      
      tr_parids  = !NULL
      pos = !NULL

      ; calculate current distance of all tracers from smoothed halo center position for galaxy members
      tr_pos[0,*] = hPos[mt.maxSnap-m,0,gcIndOrigTrMC] - tr_pos[0,*]
      tr_pos[1,*] = hPos[mt.maxSnap-m,1,gcIndOrigTrMC] - tr_pos[1,*]
      tr_pos[2,*] = hPos[mt.maxSnap-m,2,gcIndOrigTrMC] - tr_pos[2,*]
      
      correctPeriodicDistVecs, tr_pos, sP=sP ; account for periodic distance function
      
      r.relPos[mt.maxSnap-m,0,*] = tr_pos[0,*]
      r.relPos[mt.maxSnap-m,1,*] = tr_pos[1,*]
      r.relPos[mt.maxSnap-m,2,*] = tr_pos[2,*]
      tr_pos = !NULL
    endfor
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
  endif
  
  ; VELOCITY TRACERS case - will be similar to above since there could be multiple
  ; ---------------------
  if sP.trMCPerCell eq -1 then begin
      print,'Calculating new accretion trajectories using ( TracerVEL ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
  
    ; store the main arrays as a structure so we can write them directly for all tracers
    r = {relPos  : fltarr(nSnapsTot,3,n_elements(galcat.trVel_ids))    ,$
         relVel  : fltarr(nSnapsTot,3,n_elements(galcat.trVel_ids))   ,$
         curTemp : fltarr(nSnapsTot,n_elements(galcat.trVel_ids))       }

    for m=mt.maxSnap,mt.minSnap,-1 do begin
      sP.snap = m
      print,m
      ; load tracer ids and match to child ids from zMin
      h = loadSnapshotHeader(sP=sP)
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerVel',field='ids')
      
      idIndexMap = getIDIndexMap(tr_ids,minid=minid)
      trids_ind  = idIndexMap[ galcat.trVel_ids - minid ]
        
      idIndexMap = !NULL
      tr_ids     = !NULL
      
      message,'TODO: handle other parent types (stars,bh) or restrict input tracer search to gas'
      
      ; load velocities and store relative to smoothed halo CM velocity
      vel   = loadSnapshotSubset(sP=sP,partType='tracerVel',field='vel')
      
      r.relVel[mt.maxSnap-m,0,*] = vel[0,trids_ind] - hVel[mt.maxSnap-m,0,gcIndOrigTrVel] 
      r.relVel[mt.maxSnap-m,1,*] = vel[1,trids_ind] - hVel[mt.maxSnap-m,1,gcIndOrigTrVel] 
      r.relVel[mt.maxSnap-m,2,*] = vel[2,trids_ind] - hVel[mt.maxSnap-m,2,gcIndOrigTrVel] 
      
      vel = !NULL

      ; load tracer maxtemps for this timestep (best proxy for current temp without doing NN parent finds)
      tr_maxtemp = loadSnapshotSubset(sP=sP,partType='tracerVel',field='tracer_maxtemp')
      
      r.curTemp[mt.maxSnap-m,*]  = codeTempToLogK(tr_maxtemp[trids_ind])

      tr_maxtemp = !NULL
      
      ; load pos to calculate tracer positions relative to halo centers
      pos = loadSnapshotSubset(sP=sP,partType='tracerVel',field='pos')
      
      tr_pos  = pos[*,trids_ind]
      
      tr_parids  = !NULL
      pos = !NULL

      ; calculate current distance of all tracers from smoothed halo center position for galaxy members
      tr_pos[0,*] = hPos[mt.maxSnap-m,0,gcIndOrigTrVel] - tr_pos[0,*]
      tr_pos[1,*] = hPos[mt.maxSnap-m,1,gcIndOrigTrVel] - tr_pos[1,*]
      tr_pos[2,*] = hPos[mt.maxSnap-m,2,gcIndOrigTrVel] - tr_pos[2,*]
      
      correctPeriodicDistVecs, tr_pos, sP=sP ; account for periodic distance function
      
      r.relPos[mt.maxSnap-m,0,*] = tr_pos[0,*]
      r.relPos[mt.maxSnap-m,1,*] = tr_pos[1,*]
      r.relPos[mt.maxSnap-m,2,*] = tr_pos[2,*]
      tr_pos = !NULL
      
    endfor
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
  endif

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
  
  print,'halo ['+str(hInd)+'] mass: ',codeMassToLogMsun(gc.subgroupMass[hInd])
  
  nSnapsTot = mt.maxSnap - mt.minSnap + 1 ; how many values to store in time dimension

  ; get smoothed halo position with time
  hPos = smoothHaloPos(mt=mt,hInd=hInd,sP=sP)

  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    print,'Calculating single accretion trajectories using ( SPH Particles ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'

    ; get particle ids
    gas_ids = gcPIDList(gc=gc,valGCids=[hInd],partType='gas')
    dm_ids  = gcPIDList(gc=gc,valGCids=[hInd],partType='dm')

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
      calcMatch,gas_ids,ids_local,orig_ids_ind,ids_gas_ind,count=countGas
      ids_gas_ind = ids_gas_ind[calcSort(orig_ids_ind)]
      
      ids_local    = !NULL
      orig_ids_ind = !NULL
      
      ; load pos to calculate positions relative to halo centers
      pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos',inds=ids_gas_ind)

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
      u     = loadSnapshotSubset(sP=sP,partType='gas',field='u',inds=ids_gas_ind)
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='ne',inds=ids_gas_ind)
      
      ; load gas SFR and select off the effective EOS (SFR must be zero)
      sfr = loadSnapshotSubset(sP=sP,partType='gas',field='sfr',inds=ids_gas_ind)

      w = where(sfr eq 0.0,count)
      if count gt 0 then r.curTemp_gas[mt.maxSnap-m,w] = convertUtoTemp(u[w],nelec[w],/log)

      u     = !NULL
      nelec = !NULL
      sfr   = !NULL
      
      ; dark matter positions
      ids_local = loadSnapshotSubset(sP=sP,partType='dm',field='ids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of original gas_ids 
      calcMatch,dm_ids,ids_local,orig_ids_ind,ids_dm_ind,count=countDM
      ids_dm_ind = ids_dm_ind[calcSort(orig_ids_ind)]
      
      ids_local    = !NULL
      orig_ids_ind = !NULL
      
      ; load pos to calculate positions relative to halo centers
      pos = loadSnapshotSubset(sP=sP,partType='dm',field='pos',inds=ids_dm_ind)

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
