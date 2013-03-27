; accretionVel.pro
; gas accretion project - past velocity history of gas elements
; dnelson mar.2013

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
  gal_w_at   = where(at.AccTime_gal[sP.radIndHaloAcc,*] ne -1,count_gal)
  gmem_w_at  = where(at.AccTime_gmem[sP.radIndHaloAcc,*] ne -1,count_gmem)
  stars_w_at = where(at.AccTime_stars[sP.radIndHaloAcc,*] ne -1,count_stars)
  
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

      ; make gcIndOrigTr and ids of sorted tracer children
      gcIndOrigTr = mergerTreeRepParentIDs(mt=mt,galcat=galcat,sP=sP,$
                      trids_gal=galcat_gal_trids, trids_gmem=galcat_gmem_trids, trids_stars=galcat_stars_trids)

      galcat = !NULL
      
      origParIDs = { gal   : gcIndOrigTr.gal[gal_w_at]    ,$
                     stars : gcIndOrigTr.stars[stars_w_at] }
      
      gcIndOrigTr = !NULL
                     
      if min(origParIDs.gal) lt 0 or min(origParIDs.stars) lt 0 then message,'Error: Bad starting parents.'

      ; do compactMtS on origParIDs (which is just gcIndOrigTr restricted to the AtS)
      placeMap = getIDIndexMap(mt.galcatIDList,minid=minid)
      origParIDs.gal   = placeMap[origParIDs.gal-minid]
      origParIDs.stars = placeMap[origParIDs.stars-minid]
      placeMap = !NULL
      
      if min(origParIDs.gal) lt 0 or min(origParIDs.stars) lt 0 then message,'Error: Bad starting parents.'      
      
      ; only find accretion mode for those tracers with accretion times
      galcatSub = { gal   : galcat_gal_trids[gal_w_at]    ,$
                    stars : galcat_stars_trids[stars_w_at] }

      galcat_gal_trids   = !NULL
      galcat_stars_trids = !NULL
        
      ; convert the halo accretion time (scale factor) into the outer bracketing snapshot number
      bracketSnap = { gal   : intarr(n_elements(galcatSub.gal))   ,$
                      stars : intarr(n_elements(galcatSub.stars))  }
                      
      snapTimes = snapNumToRedshift(sP=sP,/all,/time)
      snapTimes = snapTimes[where(snapTimes ne -1)]
  
      bracketSnap.gal   = value_locate(snapTimes,at.accTime_gal[sP.radIndHaloAcc,gal_w_at])
      bracketSnap.stars = value_locate(snapTimes,at.accTime_stars[sP.radIndHaloAcc,stars_w_at])
      
      ; avoid first snapshot
      w = where(bracketSnap.gal eq mt.maxSnap,count)
      if count gt 0 then bracketSnap.gal[w] = mt.maxSnap - 1
      w = where(bracketSnap.stars eq mt.maxSnap,count)
      if count gt 0 then bracketSnap.stars[w] = mt.maxSnap - 1
      
      w = where(bracketSnap.gal eq -1 or bracketSnap.gal eq n_elements(snapTimes)-1,count)
      if count then message,'Error: Bad gal bracketing.'
      w = where(bracketSnap.stars eq -1 or bracketSnap.stars eq n_elements(snapTimes)-1,count)
      if count then message,'Error: Bad stars bracketing.'
      
      ; convert the galaxy accretion time (or current time for gmem) into the inner bracketing snap
      innerSnap = { gal   : intarr(n_elements(galcatSub.gal))   ,$
                    stars : intarr(n_elements(galcatSub.stars))  }
  
      innerSnap.gal   = value_locate(snapTimes,at.accTime_gal[sP.radIndGalAcc,gal_w_at])
      innerSnap.stars = value_locate(snapTimes,at.accTime_stars[sP.radIndGalAcc,stars_w_at])
      
      w = where(innerSnap.gal eq -1 or innerSnap.gal eq n_elements(snapTimes)-1,count)
      if count then message,'Error: Bad gal bracketing.'
      w = where(innerSnap.stars eq -1 or innerSnap.stars eq n_elements(snapTimes)-1,count)
      if count then message,'Error: Bad stars bracketing.'
      
      ; we store the full alpha distribution using an offset list, calculate total length now
      alphaDistSize = { gal   : fix(innerSnap.gal - bracketSnap.gal + 1)    ,$
                        stars : fix(innerSnap.stars - bracketSnap.stars + 1) }
      
      if min(alphaDistSize.gal) lt 1 or min(alphaDistSize.stars) lt 1 then $
        message,'Error: Bad alpha dist size.'
        
      alphaDistOffset = { gal   : [0,(total(alphaDistSize.gal,/cum,/int))[0:-2]]    ,$
                          stars : [0,(total(alphaDistSize.stars,/cum,/int))[0:-2]]   }
                          
      if min(alphaDistOffset.gal) lt 0 or min(alphaDistOffset.stars) lt 0 then $
        message,'Error: Bad alpha dist offset.'
      
      ; store the main arrays as a structure so we can write them directly
      r = {accAlphaDist_gal      : fltarr(total(alphaDistSize.gal,/int))-1       ,$
           accAlphaDist_stars    : fltarr(total(alphaDistSize.stars,/int))-1     ,$
           alphaDistSize         : alphaDistSize                                 ,$
           alphaDistOffset       : alphaDistOffset                                }
           
      alphaDistSize = !NULL
      alphaDistOffset = !NULL
      
      ; keep track of the number of times recorded for each tracer
      numFound = { gal   : intarr(n_elements(galcatSub.gal))   ,$
                   stars : intarr(n_elements(galcatSub.stars))  }
                
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
        save,origParIDs,snapRange,r,numFound,galcatSub,bracketSnap,innerSnap,m,filename=resFilename
        print,' --- Done! ---'
      endif
      
      ; select those tracers whose accretionvels are calculated at this snapshot
      w_gal   = where(bracketSnap.gal le sP.snap and innerSnap.gal ge sP.snap,count_gal)
      w_stars = where(bracketSnap.stars le sP.snap and innerSnap.stars ge sP.snap,count_stars)
      
      print,m,count_gal,count_stars,$
        string(total(numFound.gal)/total(r.alphaDistSize.gal)*100,format='(f5.1)')+'%',$
        string(total(numFound.stars)/total(r.alphaDistSize.stars)*100,format='(f5.1)')+'%'
      if count_gal eq 0 or count_stars eq 0 then message,'ERROR' ; not handled correctly below
      
      ; load tracer ids and match to child ids from zMin
      h = loadSnapshotHeader(sP=sP)
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')

      ; do a single match (re-sort results into the order of single_match_ids)
      single_match_ids = [galcatSub.gal[w_gal],galcatSub.stars[w_stars]]
      
      calcMatch,single_match_ids,tr_ids,galcat_ind,trids_single_ind,count=countSingle
      trids_single_ind = trids_single_ind[calcSort(galcat_ind)]
      galcat_ind = galcat_ind[calcSort(galcat_ind)]

      ; explode results back into gal,gmem,stars
      countGal   = total(galcat_ind lt count_gal,/int)
      countStars = total(galcat_ind ge count_gal,/int)
      
      if countGal ne count_gal or countStars ne count_stars then message,'Count mismatch.'
      
      trids_gal_ind    = trids_single_ind[0:countGal-1]
      trids_stars_ind  = trids_single_ind[countGal:countGal+countStars-1]
      trids_single_ind = !NULL
      
      ;galcat_gal_ind   = galcat_ind[0:countGal-1] ; unused
      ;galcat_gal_ind   -= min(galcat_gal_ind)
      galcat_stars_ind = galcat_ind[countGal:countGal+countStars-1]
      galcat_stars_ind -= min(galcat_stars_ind)
      galcat_ind = !NULL
      
      ; load tracer parents to match to gas
      tr_parids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentid')
      tr_parids_gal   = tr_parids[trids_gal_ind]
      tr_parids_stars = tr_parids[trids_stars_ind]
      
      tr_parids = !NULL
      trids_gal_ind = !NULL
      trids_stars_ind = !NULL
      
      ; load gas IDs and convert tracer parent IDs -> indices
      gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
      gasIDMap = getIDIndexMap(gas_ids,minid=minid)

      ; match star tracer parents to gas ids (keep only those with gas parents at this snapshot)
      calcMatch,gas_ids,tr_parids_stars,ind1,ind2,count=countStars ; override countStars
      ind2 = ind2[sort(ind2)]
      
      if countStars gt 0 then begin
        tr_parids_stars  = tr_parids_stars[ind2]
        galcat_stars_ind = galcat_stars_ind[ind2]
        w_stars = w_stars[ind2]
      endif
      
      ind1 = !NULL
      ind2 = !NULL
      gas_ids = !NULL
      
      tr_parids_gal  = gasIDMap[tr_parids_gal-minid]  ; convert ID->index
      if countStars gt 0 then tr_parids_stars = gasIDMap[tr_parids_stars-minid]
      gasIDMap = !NULL
      
      ; load gas positions and convert to tracer positions
      gas_pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
      tr_pos_gal  = gas_pos[*,tr_parids_gal]
      if countStars gt 0 then tr_pos_stars = gas_pos[*,tr_parids_stars]
      
      gas_pos = !NULL

      ; calculate radial unit vector of each tracer from smoothed halo center position
      for i=0,2 do $
        tr_pos_gal[i,*] = mt.hPos[mt.maxSnap-m,i,origParIDs.gal[w_gal]] - tr_pos_gal[i,*]
      correctPeriodicDistVecs, tr_pos_gal, sP=sP
      
      if countStars gt 0 then begin
        for i=0,2 do $
          tr_pos_stars[i,*] = mt.hPos[mt.maxSnap-m,i,origParIDs.stars[w_stars]] - tr_pos_stars[i,*]
        correctPeriodicDistVecs, tr_pos_stars, sP=sP
      endif
      
      ; calculate r/rvir for each
      tr_rad_gal = reform(sqrt(tr_pos_gal[0,*]^2.0 + tr_pos_gal[1,*]^2.0 + tr_pos_gal[2,*]^2.0))
      tr_rad_gal = reform(tr_rad_gal / mt.hVirRad[mt.maxSnap-m,origParIDs.gal[w_gal]])
      
      if countStars gt 0 then begin
        tr_rad_stars = reform(sqrt(tr_pos_stars[0,*]^2.0 + tr_pos_stars[1,*]^2.0 + tr_pos_stars[2,*]^2.0))
        tr_rad_stars = reform(tr_rad_stars / mt.hVirRad[mt.maxSnap-m,origParIDs.stars[w_stars]])
      endif

      ; load gas velocities and convert to tracer velocities
      gas_vel = loadSnapshotSubset(sP=sP,partType='gas',field='vel')
      tr_vel_gal  = gas_vel[*,tr_parids_gal]
      if countStars gt 0 then tr_vel_stars = gas_vel[*,tr_parids_stars]
      
      gas_vel = !NULL
      tr_parids_gal  = !NULL
      tr_parids_stars = !NULL
      
      ; make tracer velocities related to bulk velocity of primary subhalo
      for i=0,2 do $
        tr_vel_gal[i,*] = mt.hVel[mt.maxSnap-m,i,origParIDs.gal[w_gal]] - tr_vel_gal[i,*]
      
      if countStars gt 0 then $
        for i=0,2 do $
          tr_vel_stars[i,*] = mt.hVel[mt.maxSnap-m,i,origParIDs.stars[w_stars]] - tr_vel_stars[i,*]     
      
      ; normalize relative velocities by v_c(r)
      vhalo_gal = sqrt(units.G * mt.hMass[mt.maxSnap-m,origParIDs.gal[w_gal]] / $
                      mt.hVirRad[mt.maxSnap-m,origParIDs.gal[w_gal]]) ; v200 (km/s)
      vhalo_gal = sqrt( (alog(1+halo_c*tr_rad_gal)-halo_c*tr_rad_gal/(1+halo_c*tr_rad_gal)) / $
                        (alog(1+halo_c) - halo_c / (1+halo_c)) / tr_rad_gal ) * vhalo_gal ; v_c(r) (km/s)
                        
      for i=0,2 do tr_vel_gal[i,*] /= vhalo_gal

      vhalo_gal = !NULL
      
      if countStars gt 0 then begin
        vhalo_stars = sqrt(units.G * mt.hMass[mt.maxSnap-m,origParIDs.stars[w_stars]] / $
                        mt.hVirRad[mt.maxSnap-m,origParIDs.stars[w_stars]]) ; v200 (km/s)
        vhalo_stars = sqrt( (alog(1+halo_c*tr_rad_stars)-halo_c*tr_rad_stars/(1+halo_c*tr_rad_stars)) / $
                          (alog(1+halo_c) - halo_c / (1+halo_c)) / tr_rad_stars ) * vhalo_stars ; v_c(r) (km/s)
                          
        for i=0,2 do tr_vel_stars[i,*] /= vhalo_stars
        vhalo_stars = !NULL
      endif
      
      ; convert radial distances back into kpc (rnorm)
      tr_rad_gal   *= mt.hVirRad[mt.maxSnap-m,origParIDs.gal[w_gal]]
      tr_rad_stars *= mt.hVirRad[mt.maxSnap-m,origParIDs.stars[w_stars]]
      
      ; decompose velocity vectors into radial and tangential projections (gal)
      ; vradvec = (vvec dot rvec) / rnorm * rvec, vtanvec = vvec - vradvec
      const_gal = (tr_vel_gal[0,*] * tr_pos_gal[0,*] + $
                   tr_vel_gal[1,*] * tr_pos_gal[1,*] + $
                   tr_vel_gal[2,*] * tr_pos_gal[2,*]) / (tr_rad_gal*tr_rad_gal)
                  
      vrad_gal = fltarr(3,countGal)
      for i=0,2 do vrad_gal[i,*] = const_gal * tr_pos_gal[i,*]
      const_gal = !NULL
      
      vtan_gal = fltarr(3,countGal)
      for i=0,2 do vtan_gal[i,*] = tr_vel_gal[i,*] - vrad_gal[i,*]
      
      ; override projection vectors with norms
      vrad_gal = sqrt(vrad_gal[0,*]^2.0 + vrad_gal[1,*]^2.0 + vrad_gal[2,*]^2.0)
      vtan_gal = sqrt(vtan_gal[0,*]^2.0 + vtan_gal[1,*]^2.0 + vtan_gal[2,*]^2.0)
      
      ; save alpha=radial/tangential ratio for each tracer
      alpha_gal = reform(vrad_gal / vtan_gal)
      
      vrad_gal = !NULL
      vtan_gal = !NULL
      
      ; get offsets into save table and increment found counters
      curOffsets = r.alphaDistOffset.gal[w_gal] + numFound.gal[w_gal]
      numFound.gal[w_gal] += 1
      
      r.accAlphaDist_gal[curOffsets] = alpha_gal
      curOffsets = !NULL
      alpha_gal = !NULL
      
      ; now do stars
      const_stars = (tr_vel_stars[0,*] * tr_pos_stars[0,*] + $
                     tr_vel_stars[1,*] * tr_pos_stars[1,*] + $
                     tr_vel_stars[2,*] * tr_pos_stars[2,*]) / (tr_rad_stars*tr_rad_stars)
                  
      vrad_stars = fltarr(3,countStars)
      for i=0,2 do vrad_stars[i,*] = const_stars * tr_pos_stars[i,*]
      const_stars = !NULL
      
      vtan_stars = fltarr(3,countStars)
      for i=0,2 do vtan_stars[i,*] = tr_vel_stars[i,*] - vrad_stars[i,*]
      
      ; override projection vectors with norms
      vrad_stars = sqrt(vrad_stars[0,*]^2.0 + vrad_stars[1,*]^2.0 + vrad_stars[2,*]^2.0)
      vtan_stars = sqrt(vtan_stars[0,*]^2.0 + vtan_stars[1,*]^2.0 + vtan_stars[2,*]^2.0)
      
      ; save alpha=radial/tangential ratio for each tracer
      alpha_stars = reform(vrad_stars / vtan_stars)
      
      vrad_stars = !NULL
      vtan_stars = !NULL
      
      ; get offsets into save table and increment found counters
      curOffsets = r.alphaDistOffset.stars[w_stars] + numFound.stars[w_stars]
      numFound.stars[w_stars] += 1
      
      r.accAlphaDist_stars[curOffsets] = alpha_stars
      curOffsets = !NULL
      alpha_stars = !NULL
      
    endfor ;m
    
    ; verify we found an accretion mode for every gas particle
    if ~array_equal(numFound.gal,r.alphaDistSize.gal) then message,'Error: Not all gal found.'
    ;if ~array_equal(numFound.stars,r.alphaDistSize.stars) then message,'Error: Not all stars found (gassrconly).'
    
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
  sPg = simParams(res=128,run='gadget',redshift=2.0)
  
  haloID = 304 ;z2.304 z2.301 z2.130 z2.64
  gcID = getMatchedIDs(sPa=sP,sPg=sPg,haloID=haloID) 
  
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
    
    ; make parent list for groupcat
    gcIndOrigTr = mergerTreeRepParentIDs(mt=mt,galcat=galcat,sP=sP,$
      trids_gal=galcat_gal_trids, trids_gmem=galcat_gmem_trids, trids_stars=galcat_stars_trids)
  
    galcat = !NULL
    
    origParIDs = { gal   : gcIndOrigTr.gal[wAm.gal]    ,$
                   stars : gcIndOrigTr.stars[wAm.stars] }

    gcIndOrigTr = !NULL
  
    ; only have information for those tracers with halo accretion times
    galcatSub = { gal   : galcat_gal_trids[wAm.gal]    ,$
                  stars : galcat_stars_trids[wAm.stars] }
                  
    galcat_gal_trids = !NULL
    galcat_gmem_trids = !NULL
    galcat_star_trids = !NULL
    
    av = accretionVel(sP=sP)
    
    ; select this halo
    w_gal   = where(origParIDs.gal eq gcID.a,count_gal)
    w_stars = where(origParIDs.stars eq gcID.a,count_stars)
    
    ; cut out our subset
    accAlphaDist = { gal   : fltarr(total(av.alphaDistSize.gal[w_gal],/int))    ,$
                     stars : fltarr(total(av.alphaDistSize.stars[w_stars],/int)) }
                     
    ; make sure subsets are contiguous
    if max(w_gal)-min(w_gal) ne count_gal-1 or max(w_stars)-min(w_stars) ne count_stars-1 then message,'uhoh'
    
    accAlphaDist.gal = av.accAlphaDist_gal[ av.alphaDistOffset.gal[w_gal[0]] : $
      av.alphaDistOffset.gal[w_gal[-1]] + av.alphaDistSize.gal[w_gal[-1]] - 1]
    
    accAlphaDist.stars = av.accAlphaDist_stars[ av.alphaDistOffset.stars[w_stars[0]] : $
      av.alphaDistOffset.stars[w_stars[-1]] + av.alphaDistSize.stars[w_stars[-1]] - 1]
      
    alphaDistSize   = { gal   : av.alphaDistSize.gal[w_gal], stars : av.alphaDistSize.stars[w_stars] }
    alphaDistOffset = { gal   : av.alphaDistOffset.gal[w_gal] - av.alphaDistOffset.gal[w_gal[0]]        ,$
                        stars : av.alphaDistOffset.stars[w_stars] - av.alphaDistOffset.stars[w_stars[0]] }
      
    ; load maxtemps
    maxTemp = gcSubsetProp(sP=sP,select='pri',/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
    
    maxTemps = { gal : maxTemp.gal[w_gal], stars : maxTemp.stars[w_stars] }

    ; calculate elapsed time 1.0-0.15 rvir in Myr
    accTime = { gal   : redshiftToAgeFlat(1/at.AccTime_gal[sP.radIndGalAcc,wAm.gal[w_gal]]-1) - $
                        redshiftToAgeFlat(1/at.AccTime_gal[sP.radIndHaloAcc,wAm.gal[w_gal]]-1)   ,$
                stars : redshiftToAgeFlat(1/at.AccTime_stars[sP.radIndGalAcc,wAm.stars[w_stars]]-1) - $
                        redshiftToAgeFlat(1/at.AccTime_stars[sP.radIndHaloAcc,wAm.stars[w_stars]]-1)    }
    save,accAlphaDist,alphaDistSize,alphaDistOffset,maxTemps,accTime,filename=saveFilename
    
  endelse
  
  ; plots
  binsize = 0.1
  
  start_PS, sP.plotPath + 'vel.h'+str(haloID)+'.'+sP.plotPrefix+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=[0,5],yrange=[0.0,1.1],xs=1,ys=1,$
      ytitle="N",xtitle="alpha"
    
    
    hist = histogram(accAlphaDist.gal,binsize=binsize,loc=loc)
    cgPlot,loc+0.5*binsize,hist/float(max(hist)),line=0,/overplot
    
  end_PS
  
  ; individual alpha distributions (random)
  binsize = 0.2
  
  start_PS, sP.plotPath + 'vel2.h'+str(haloID)+'.'+sP.plotPrefix+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    !p.multi = [0,3,3]
    
    inds = fix(randomu(seed,!p.multi[1]*!p.multi[2])*n_elements(alphaDistSize.gal))
    print,inds
    
    densRes = 100
    
    foreach ind,inds do begin
      ; plot the alpha distribution of this tracer's history
      tar_inds = lindgen(alphaDistSize.gal[ind]) + alphaDistOffset.gal[ind]
      
      hist = histogram(accAlphaDist.gal[tar_inds],binsize=binsize,loc=loc)
      cgPlot,loc+0.5*binsize,hist,line=0,charsize=!p.charsize-0.5,xrange=[0,max(loc)],/xs
      
      ; use an adaptive kernel density estimator
      ;densPts = findgen(densRes)/(densRes-1) * $
      ;  (max(accAlphaDist.gal[tar_inds])-min(accAlphaDist.gal[tar_inds])) + $
      ;  min(accAlphaDist.gal[tar_inds])
      ;dens = akde(accAlphaDist.gal[tar_inds],densPts)
      ;cgPlot,densPts,dens,line=1,color=cgColor('red'),/overplot
      
      ; weighted centroid
      wtCen = total(loc*hist)/total(hist)
      wtCen2 = mean(accAlphaDist.gal[tar_inds])
      cgPlot,[wtCen,wtCen],[0,max(hist)],line=1,color=cgColor('red'),/overplot
      cgPlot,[wtCen2,wtCen2],[0,max(hist)],line=1,color=cgColor('blue'),/overplot
      cgPlot,[1.0,1.0],[0,max(hist)],line=0,color=cgColor('green'),/overplot
    endforeach
    
  end_PS
  
  ; calculate weighted centroid distribution
  wtCens = fltarr(n_elements(alphaDistSize.gal))
  wtMin  = fltarr(n_elements(alphaDistSize.gal))
  
  for i=0,n_elements(alphaDistSize.gal)-1 do begin
    tar_inds = lindgen(alphaDistSize.gal[i]) + alphaDistOffset.gal[i]
    wtCens[i] = mean(accAlphaDist.gal[tar_inds])
    wtMin[i]  = min(accAlphaDist.gal[tar_inds])
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
