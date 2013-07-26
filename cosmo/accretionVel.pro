; accretionVel.pro
; gas accretion project - past velocity history of gas elements
; dnelson jul.2013

; -----------------------------------------------------------------------------------------------------
; accretionVel(): for each gas particle/tracer, consider it over the time range between it accretes on
;                 to the galaxy (or where it starts in the halo for gmem) back to its 1.0rvir crossing 
;                 time and calculate the evolution of the ratio of the radial to tangential velocity 
;                 component over that period
; -----------------------------------------------------------------------------------------------------

function accretionVel, sP=sP, restart=restart

  forward_function cosmoTracerChildren, cosmoTracerVelParents, accretionTimes
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  halo_c = 5.0 ; good enough (constant) halo concentration for this mass/redshift range

  ; first, walk back through the merger tree and find primary subhalos with good parent histories
  mt = mergerTreeSubset(sP=sP)
  at = accretionTimes(sP=sP)

  ; set saveFilename and check for existence  
  saveFilename = sP.derivPath + 'accVel.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(mt.maxSnap)+'-'+str(mt.minSnap)+'.sav'  
  
  if file_test(saveFilename) then begin
    restore, saveFilename
    return, r
  endif
  
  resFilename = sP.derivPath + 'accVel.'+sP.saveTag+'.restart.'+sP.savPrefix+str(sP.res)+'.'+$
                str(mt.maxSnap)+'-'+str(mt.minSnap)+'.sav'
  
  ; select those particles/tracers with recorded halo accretion times
  w_at = where(at.AccTime[sP.radIndHaloAcc,*] ne -1,count_w_at)
  
  ; load galaxy/group member catalogs at zMin for gas ids to search for
  galcat = galaxyCat(sP=sP)
  
  count_bracketed = 0UL
  
  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    message,'todo'
  endif  
  
  ; MONTE CARLO TRACERS CASE - for all gas cells, track back all child tracers
  ; ------------------------
  if sP.trMCPerCell gt 0 then begin
  
    if ~file_test(resFilename) then begin
      ; no restart
      print,'Calculating new accretion velocities using ( TracerMC ) res = '+str(sP.res)+$
        ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'

      message,'TODO: check that this works for the global galcat.trMC_ids (before was only for gal,stars)'
        
      ; parents
      origParIDs = gcIndOrigTrMC[w_at]
                     
      if min(origParIDs.gal) lt 0 or min(origParIDs.stars) lt 0 then message,'Error: Bad starting parents.'
      
      ; only find accretion mode for those tracers with accretion times
      galcatSub = galcat.trMC_ids[w_at]
        
      ; convert the halo accretion time (scale factor) into the outer bracketing snapshot number  
      snapTimes = snapNumToRedshift(sP=sP,/all,/time)
      snapTimes = snapTimes[where(snapTimes ne -1)]
  
      bracketSnap = value_locate(snapTimes,at.accTime[sP.radIndHaloAcc,w_at])
      
      ; avoid first snapshot
      w = where(bracketSnap eq mt.maxSnap,count)
      if count gt 0 then bracketSnap[w] = mt.maxSnap - 1
      
      w = where(bracketSnap eq -1 or bracketSnap eq n_elements(snapTimes)-1,count)
      if count then message,'Error: Bad bracketing.'
      
      ; convert the galaxy accretion time (or current time for gmem) into the inner bracketing snap
      innerSnap = value_locate(snapTimes,at.accTime[sP.radIndGalAcc,w_at])
      
      w = where(innerSnap eq -1 or innerSnap eq n_elements(snapTimes)-1,count)
      if count then message,'Error: Bad bracketing.'

      ; we store the full alpha distribution using an offset list, calculate total length now
      alphaDistSize   = fix(innerSnap - bracketSnap + 1)
      alphaDistOffset = [0,(total(alphaDistSize,/cum,/int))[0:-2]]
      
      if min(alphaDistSize) lt 1 then message,'Error: Bad alpha dist size.'      
      if min(alphaDistOffset) lt 0 then message,'Error: Bad alpha dist offset.'
      
      ; store the main arrays as a structure so we can write them directly
      r = {accAlphaDist      : fltarr(total(alphaDistSize,/int))-1       ,$
           alphaDistSize     : alphaDistSize                             ,$
           alphaDistOffset   : alphaDistOffset                            }
           
      alphaDistSize = !NULL
      alphaDistOffset = !NULL
      
      ; keep track of the number of times recorded for each tracer
      numFound = intarr(n_elements(galcat.trMC_ids))
                
      snapRange = [mt.maxSnap-1,mt.minSnap]
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
        save,origParIDs,snapRange,r,numFound,bracketSnap,innerSnap,m,filename=resFilename
        print,' --- Done! ---'
      endif
      
      ; select those tracers whose accretionvels are calculated at this snapshot
      w_now = where(bracketSnap le sP.snap and innerSnap ge sP.snap,count_now)
      
      print,m,count_now,count_stars,$
        string(total(numFound)/total(r.alphaDistSize)*100,format='(f5.1)')+'%'
      if count_now eq 0 then message,'ERROR' ; not handled correctly below
      
      ; load tracer ids and match to child ids from zMin
      h = loadSnapshotHeader(sP=sP)
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')

      ; match
      calcMatch,galcat.trMC_ids[w_now],tr_ids,galcat_ind,trids_ind,count=countMatch
      trids_ind = trids_ind[calcSort(galcat_ind)]
      galcat_ind = galcat_ind[calcSort(galcat_ind)]

      ; load tracer parents to match to gas
      tr_parids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentid')
      tr_parids = tr_parids[trids_ind]
      
      trids_ind = !NULL
      
      message,'TODO: handle different parent types (we need only pos,vel, so can get for all)'
      
      ; load gas IDs and convert tracer parent IDs -> indices
      gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
      gasIDMap = getIDIndexMap(gas_ids,minid=minid)

      ; match tracer parents to gas ids (keep only those with gas parents at this snapshot)
      calcMatch,gas_ids,tr_parids,ind1,ind2,count=countNow ; override countStars
      ind2 = ind2[sort(ind2)]
      
      if countNow eq 0 then message,'Error: Not handled correctly below.'
      
      tr_parids  = tr_parids[ind2]
      galcat_ind = galcat_ind[ind2]
      w_now = w_now[ind2]
      
      ind1 = !NULL
      ind2 = !NULL
      gas_ids = !NULL
      
      tr_parids = gasIDMap[tr_parids-minid]  ; convert ID->index
      gasIDMap = !NULL
      
      ; load gas positions and convert to tracer positions
      gas_pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
      tr_pos  = gas_pos[*,tr_parids]
      gas_pos = !NULL

      ; calculate radial unit vector of each tracer from smoothed halo center position
      for i=0,2 do $
        tr_pos[i,*] = mt.hPos[mt.maxSnap-m,i,origParIDs[w_now]] - tr_pos[i,*]
      correctPeriodicDistVecs, tr_pos, sP=sP
            
      ; calculate r/rvir for each
      tr_rad = reform(sqrt(tr_pos[0,*]^2.0 + tr_pos[1,*]^2.0 + tr_pos[2,*]^2.0))
      tr_rad = reform(tr_rad / mt.hVirRad[mt.maxSnap-m,origParIDs[w_now]])
      
      ; load gas velocities and convert to tracer velocities
      gas_vel = loadSnapshotSubset(sP=sP,partType='gas',field='vel')
      tr_vel  = gas_vel[*,tr_parids]
      gas_vel = !NULL
      
      tr_parids = !NULL
      
      ; make tracer velocities related to bulk velocity of primary subhalo
      for i=0,2 do $
        tr_vel[i,*] = mt.hVel[mt.maxSnap-m,i,origParIDs[w_now]] - tr_vel[i,*]
      
      ; normalize relative velocities by v_c(r)
      vhalo = sqrt(units.G * mt.hMass[mt.maxSnap-m,origParIDs[w_now]] / $
                   mt.hVirRad[mt.maxSnap-m,origParIDs[w_now]]) ; v200 (km/s)
      vhalo = sqrt( (alog(1+halo_c*tr_rad)-halo_c*tr_rad/(1+halo_c*tr_rad)) / $
                        (alog(1+halo_c) - halo_c / (1+halo_c)) / tr_rad ) * vhalo ; v_c(r) (km/s)
                        
      for i=0,2 do tr_vel[i,*] /= vhalo

      vhalo = !NULL
            
      ; convert radial distances back into kpc (rnorm)
      tr_rad *= mt.hVirRad[mt.maxSnap-m,origParIDs[w_now]]
      
      ; decompose velocity vectors into radial and tangential projections (gal)
      ; vradvec = (vvec dot rvec) / rnorm * rvec, vtanvec = vvec - vradvec
      const_now = (tr_vel[0,*] * tr_pos[0,*] + $
               tr_vel[1,*] * tr_pos[1,*] + $
               tr_vel[2,*] * tr_pos[2,*]) / (tr_rad*tr_rad)
                  
      vrad = fltarr(3,countNow)
      for i=0,2 do vrad[i,*] = const_now * tr_pos[i,*]
      const_now = !NULL
      
      vtan = fltarr(3,countNow)
      for i=0,2 do vtan[i,*] = tr_vel[i,*] - vrad[i,*]
      
      ; override projection vectors with norms
      vrad = sqrt(vrad[0,*]^2.0 + vrad[1,*]^2.0 + vrad[2,*]^2.0)
      vtan = sqrt(vtan[0,*]^2.0 + vtan[1,*]^2.0 + vtan[2,*]^2.0)
      
      ; save alpha=radial/tangential ratio for each tracer
      alpha = reform(vrad / vtan)
      
      vrad = !NULL
      vtan = !NULL
      
      ; get offsets into save table and increment found counters
      curOffsets = r.alphaDistOffset[w_now] + numFound[w_now]
      numFound[w_now] += 1
      
      r.accAlphaDist[curOffsets] = alpha
      curOffsets = !NULL
      alpha = !NULL
            
    endfor ;m
    
    ; verify we found an accretion mode for every gas particle
    if ~array_equal(numFound,r.alphaDistSize) then message,'Error: Not all found.'
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
  endif
  
  ; VELOCITY TRACERS case - will be similar to above since there could be multiple
  ; ---------------------
  if sP.trMCPerCell eq -1 then begin
    message,'todo'
  endif  
  
end

; velSingleHalo(): explore above for a single halo

pro velSingleHalo

  sP = simParams(res=128,run='tracer',redshift=2.0)
  
  haloID = 304 ;z2.304 z2.301 z2.130 z2.64
  gcID = getMatchedIDs(simParams=sP,haloID=haloID) 

  accMode = 'smooth'
  
  saveFilename = sP.derivPath + 'accVel.h'+str(haloID)+'.am-'+accMode+'.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.sav'
  
  if file_test(saveFilename) then begin
    restore, saveFilename, /verbose
  endif else begin
    
    ; load
    mt = mergerTreeSubset(sP=sP)
    at = accretionTimes(sP=sP)
    galcat = galaxyCat(sP=sP)
    
    wAm = accModeInds(at=at,accMode=accMode,sP=sP)
    
    origParIDs = mt.gcIndOrigTrMC[wAm.gal]
    
    av = accretionVel(sP=sP)
    
    message,'Move av into gcSubsetProp, unify with maxTemp'
    
    ; select this halo
    w_halo = where(origParIDs eq gcID,count)
    
    ; cut out our subset
    accAlphaDist = fltarr(total(av.alphaDistSize[w_halo],/int))
                     
    ; make sure subsets are contiguous
    if max(w_halo)-min(w_halo) ne count-1 then message,'uhoh'
    
    accAlphaDist = av.accAlphaDist[ av.alphaDistOffset[w_halo[0]] : $
      av.alphaDistOffset[w_halo[-1]] + av.alphaDistSize[w_halo[-1]] - 1]
    
    alphaDistSize   = av.alphaDistSize[w_halo]
    alphaDistOffset = av.alphaDistOffset[w_halo] - av.alphaDistOffset[w_halo[0]]
      
    ; load maxtemps
    maxTemp = gcSubsetProp(sP=sP,gcIDList=[gcID],/maxPastTemp,/accretionTimeSubset,accMode=accMode)
    
    ; calculate elapsed time 1.0-0.15 rvir in Myr
    accTime = redshiftToAgeFlat(1/at.AccTime[sP.radIndGalAcc,wAm.gal[w_gal]]-1) - $
              redshiftToAgeFlat(1/at.AccTime[sP.radIndHaloAcc,wAm.gal[w_gal]]-1)
    save,accAlphaDist,alphaDistSize,alphaDistOffset,maxTemps,accTime,filename=saveFilename
    
  endelse
  
  ; plots
  binsize = 0.1
  
  start_PS, sP.plotPath + 'vel.h'+str(haloID)+'.'+sP.plotPrefix+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=[0,5],yrange=[0.0,1.1],xs=1,ys=1,$
      ytitle="N",xtitle="alpha"
    
    
    hist = histogram(accAlphaDist,binsize=binsize,loc=loc)
    cgPlot,loc+0.5*binsize,hist/float(max(hist)),line=0,/overplot
    
  end_PS
  
  ; individual alpha distributions (random)
  binsize = 0.2
  
  start_PS, sP.plotPath + 'vel2.h'+str(haloID)+'.'+sP.plotPrefix+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    !p.multi = [0,3,3]
    
    inds = fix(randomu(seed,!p.multi[1]*!p.multi[2])*n_elements(alphaDistSize))
    print,inds
    
    densRes = 100
    
    foreach ind,inds do begin
      ; plot the alpha distribution of this tracer's history
      tar_inds = lindgen(alphaDistSize[ind]) + alphaDistOffset[ind]
      
      hist = histogram(accAlphaDist[tar_inds],binsize=binsize,loc=loc)
      cgPlot,loc+0.5*binsize,hist,line=0,charsize=!p.charsize-0.5,xrange=[0,max(loc)],/xs
      
      ; use an adaptive kernel density estimator
      ;densPts = findgen(densRes)/(densRes-1) * $
      ;  (max(accAlphaDist[tar_inds])-min(accAlphaDist[tar_inds])) + $
      ;  min(accAlphaDist[tar_inds])
      ;dens = akde(accAlphaDist[tar_inds],densPts)
      ;cgPlot,densPts,dens,line=1,color=cgColor('red'),/overplot
      
      ; weighted centroid
      wtCen = total(loc*hist)/total(hist)
      wtCen2 = mean(accAlphaDist[tar_inds])
      cgPlot,[wtCen,wtCen],[0,max(hist)],line=1,color=cgColor('red'),/overplot
      cgPlot,[wtCen2,wtCen2],[0,max(hist)],line=1,color=cgColor('blue'),/overplot
      cgPlot,[1.0,1.0],[0,max(hist)],line=0,color=cgColor('green'),/overplot
    endforeach
    
  end_PS
  
  ; calculate weighted centroid distribution
  wtCens = fltarr(n_elements(alphaDistSize))
  wtMin  = fltarr(n_elements(alphaDistSize))
  
  for i=0,n_elements(alphaDistSize)-1 do begin
    tar_inds = lindgen(alphaDistSize[i]) + alphaDistOffset[i]
    wtCens[i] = mean(accAlphaDist[tar_inds])
    wtMin[i]  = min(accAlphaDist[tar_inds])
  endfor
  
  ; tmax vs. alpha correlation?
  start_PS, sP.plotPath + 'vel.alpha.mean.tmax.h'+str(haloID)+'.'+sP.plotPrefix+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    cgPlot,maxTemps.gal,wtCens,psym=3,yrange=[0,5],xtitle="log Tmax",ytitle="alpha mean"
  end_PS
  
  start_PS, sP.plotPath + 'vel.alpha.min.tmax.h'+str(haloID)+'.'+sP.plotPrefix+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    cgPlot,maxTemps.gal,wtMin,psym=3,yrange=[0,2],xtitle="log Tmax",ytitle="alpha min"
  end_PS
  
  ; alpha mean/min distributions
  start_PS, sP.plotPath + 'vel.alpha.mean.h'+str(haloID)+'.'+sP.plotPrefix+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    hist=histogram(wtCens,binsize=binsize,loc=loc)
    cgPlot,loc+0.5*binsize,wtCens,line=0,xrange=[0,10]
  end_PS
  
  start_PS, sP.plotPath + 'vel.alpha.min.h'+str(haloID)+'.'+sP.plotPrefix+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    hist=histogram(wtMin,binsize=binsize,loc=loc)
    cgPlot,loc+0.5*binsize,wtCens,line=0,xrange=[0,5]
  end_PS
  
  stop
  
end
