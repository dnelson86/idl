; accretionTimes.pro
; gas accretion project - past radial history of gas elements (virial radius crossing)
; dnelson aug.2012

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

  ; first, walk back through the merger tree and find primary subhalos with good parent histories
  nVirFacs = n_elements(sP.rVirFacs)
  mt = mergerTreeSubset(sP=sP,/verbose)
  snapRange = [mt.maxSnap,mt.minSnap]
  
  ; set saveFilename and check for existence
  saveFilename = sP.derivPath + 'accTimesAda.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(mt.maxSnap)+'-'+str(mt.minSnap)+'.sav'  
  
  if file_test(saveFilename) then begin
    restore, saveFilename
    return, r
  endif
  
  resFilename = sP.derivPath + 'accTimesAda.'+sP.saveTag+'.restart.'+sP.savPrefix+str(sP.res)+'.'+$
                str(mt.maxSnap)+'-'+str(mt.minSnap)+'.sav'
  
  ; load galaxy/group member catalogs at zMin for gas ids to search for
  galcat = galaxyCat(sP=sP)
  
  ; count of how many particles/tracers we tracked through r_vir
  accCount = { gal   : ulonarr(nVirFacs), gmem    : ulonarr(nVirFacs), stars : ulonarr(nVirFacs), $
               galRT : ulonarr(nVirFacs), starsRT : ulonarr(nVirFacs) }
  prevTime = 0 ; scale factor at previous snapshot

  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    print,'Calculating new accretion time using ( SPH Particles ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
      
    if ~file_test(resFilename) then begin ; no restart
      ; replicate hMinSnap for each child gas element
      gasMinSnap = { gal   : mt.hMinSnap[mt.gcIndOrig.gal]   ,$
                     gmem  : mt.hMinSnap[mt.gcIndOrig.gmem]  ,$
                     stars : mt.hMinSnap[mt.gcIndOrig.stars]  }  
  
      ; store the r/rvir of each at the previous snapshot for interpolation, and as a monotonic check
      prevRad = { gal   : fltarr(n_elements(mt.galcatSub.gal))  ,$
                  gmem  : fltarr(n_elements(mt.galcatSub.gmem)) ,$
                  stars : fltarr(n_elements(mt.galcatSub.stars)) }
      
      accMask = { gal   : bytarr(nVirFacs,n_elements(mt.galcatSub.gal))  ,$
                  gmem  : bytarr(nVirFacs,n_elements(mt.galcatSub.gmem)) ,$
                  stars : bytarr(nVirFacs,n_elements(mt.galcatSub.stars)) }
      
      ; store the main arrays as a structure so we can write them directly
      r = {accTime_gal       : fltarr(nVirFacs,n_elements(mt.galcatSub.gal))-1   ,$
           accTime_gmem      : fltarr(nVirFacs,n_elements(mt.galcatSub.gmem))-1  ,$
           accTime_stars     : fltarr(nVirFacs,n_elements(mt.galcatSub.stars))-1 ,$
           accTimeRT_gal     : fltarr(n_elements(mt.galcatSub.gal))-1            ,$
           accTimeRT_stars   : fltarr(n_elements(mt.galcatSub.stars))-1          ,$
           accHaloTvir_gal   : fltarr(n_elements(mt.galcatSub.gal))              ,$
           accHaloTvir_gmem  : fltarr(n_elements(mt.galcatSub.gmem))             ,$
           accHaloTvir_stars : fltarr(n_elements(mt.galcatSub.stars))            ,$
           rVirFacs          : sP.rVirFacs                                        }
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
        save,prevRad,gasMinSnap,accMask,r,prevTime,accCount,m,filename=resFilename
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
      
      match,galcat.stellarIDs[mt.galcatSub.stars],ids,galcat_stars_ind,ids_stars_ind,count=countStars,/sort
      ids_stars_ind    = ids_stars_ind[sort(galcat_stars_ind)]
      galcat_stars_ind = galcat_stars_ind[sort(galcat_stars_ind)] ; use to incompletely fill r.x_stars

      ids        = !NULL
      galcat_ind = !NULL
      
      if countGal ne n_elements(mt.galcatSub.gal) or countGmem ne n_elements(mt.galcatSub.gmem) then $
        message,'Error: Failed to locate all of gal/gmem in gas_ids (overflow64?).'      
      
      ; load density,temp to apply galaxy cut
      u = loadSnapshotSubset(sP=sP,partType='gas',field='u')
      u_gal   = u[ids_gal_ind]
      if countStars gt 0 then u_stars = u[ids_stars_ind]
      u = !NULL
      
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
      nelec_gal   = nelec[ids_gal_ind]
      if countStars gt 0 then nelec_stars = nelec[ids_stars_ind]
      nelec = !NULL
      
      temp_gal = convertUtoTemp(u_gal,nelec_gal)
      u_gal = !NULL
      nelec_gal = !NULL
      
      if countStars gt 0 then begin
        temp_stars = convertUtoTemp(u_stars,nelec_stars)
        u_stars = !NULL
        nelec_stars = !NULL
      endif
      
      ; load rho of gas and make galaxy (rho,temp) plane cut
      dens = loadSnapshotSubset(sP=sP,partType='gas',field='density')
      dens_gal = dens[ids_gal_ind]
      if countStars gt 0 then dens_stars = dens[ids_stars_ind]
      dens = !NULL
      
      ; scale Torrey+ (2011) galaxy cut to physical density
      scalefac = snapNumToRedshift(sP=sP,/time) ; time flag gives simulation time = scale factor
      a3inv = 1.0 / (scalefac*scalefac*scalefac)
      dens_gal *= a3inv
      if countStars gt 0 then dens_stars *= a3inv
      
      ; mark any galaxy gas failing cut as accreted at this snapshot (if not previously marked)
      wGalCut = where(alog10(temp_gal) - sP.galcut_rho * alog10(dens_gal) ge sP.galcut_T and $
                      r.accTimeRT_gal eq -1,countGalCut)
               
      if countGalCut gt 0 then r.accTimeRT_gal[wGalCut] = h.time
      accCount.galRT += countGalCut

      wGalCut  = !NULL
      dens_gal = !NULL
      temp_gal = !NULL
      
      ; mark any gas hosting stellar tracers as accreted as above
      if countStars gt 0 then begin
        wStarsCut = where(alog10(temp_stars) - sP.galcut_rho * alog10(dens_stars) ge sP.galcut_T and $
                        r.accTimeRT_stars eq -1,countStarsCut)
                 
        if countStarsCut gt 0 then r.accTimeRT_stars[galcat_stars_ind[wStarsCut]] = h.time
        accCount.starsRT += countStarsCut
        
        wStarsCut  = !NULL
        dens_stars = !NULL
        temp_stars = !NULL
      endif else begin
        countStarsCut = 0
      endelse
      
      print,' ['+string(m,format='(i3)')+'] [-] (rhot)'+$
        ' accreted now counts '+string(countGalCut,format='(i7)')+' ('+$
        string(float(countGalCut)/n_elements(mt.gcIndOrig.gal)*100,format='(f4.1)')+'%)   stars '+$
        string(countStarsCut,format='(i7)')+' ('+$
        string(float(countStarsCut)/n_elements(mt.gcIndOrig.stars)*100,format='(f4.1)')+'%) || cum '+$
        string(accCount.galRT,format='(i7)')+' ('+$
        string(float(accCount.galRT)/n_elements(mt.gcIndOrig.gal)*100,format='(f4.1)')+'%)   stars '+$
        string(accCount.starsRT,format='(i7)')+' ('+$
        string(float(accCount.starsRT)/n_elements(mt.gcIndOrig.stars)*100,format='(f4.1)')+'%)'
      
      ; load pos to calculate radii
      pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
      
      pos_gal  = pos[*,ids_gal_ind]
      pos_gmem = pos[*,ids_gmem_ind]
      if countStars gt 0 then pos_stars = pos[*,ids_stars_ind]
      
      pos = !NULL

      ; calculate current distance of gas particle from smoothed halo center position for galaxy members
      gal_pri  = periodicDists(reform(mt.hPos[mt.maxSnap-m,*,mt.gcIndOrig.gal]),pos_gal,sP=sP)
      gal_pri /= mt.hVirRad[mt.maxSnap-m,mt.gcIndOrig.gal]
      pos_gal = !NULL
      
      ; for group members
      gmem_pri = periodicDists(reform(mt.hPos[mt.maxSnap-m,*,mt.gcIndOrig.gmem]),pos_gmem,sP=sP)
      gmem_pri /= mt.hVirRad[mt.maxSnap-m,mt.gcIndOrig.gmem]
      pos_gmem = !NULL
      
      ; stars
      if countStars gt 0 then begin
        stars_pri = periodicDists($
          reform(mt.hPos[mt.maxSnap-m,*,mt.gcIndOrig.stars[galcat_stars_ind]]),pos_stars,sP=sP)
        stars_pri /= mt.hVirRad[mt.maxSnap-m,mt.gcIndOrig.stars[galcat_stars_ind]]
        pos_stars = !NULL
      endif
      
      ; loop over each target radius
      foreach rVirFac,sP.rVirFacs,k do begin
        ; for particles who are still within r_vir, check if they have passed beyond
        gal_w   = where(gal_pri ge rVirFac and prevRad.gal lt rVirFac and accMask.gal[k,*] eq 0B,count_gal)
        gmem_w  = where(gmem_pri ge rVirFac and prevRad.gmem lt rVirFac and accMask.gmem[k,*] eq 0B,count_gmem)

        if countStars gt 0 then begin
          stars_w = where(stars_pri ge rVirFac and prevRad.stars[galcat_stars_ind] lt rVirFac and $
                          accMask.stars[k,galcat_stars_ind] eq 0B,count_stars)
        endif else begin
          count_stars = 0 & stars_pri = 0
        endelse
        
        print,' ['+string(m,format='(i3)')+'] ['+str(k)+'] r='+string(rVirFac,format='(f4.2)')+$
          ' accreted now counts '+string(count_gal,format='(i7)')+' ('+$
          string(float(count_gal)/n_elements(mt.gcIndOrig.gal)*100,format='(f4.1)')+'%) '+$
          string(count_gmem,format='(i7)')+' ('+$
          string(float(count_gmem)/n_elements(mt.gcIndOrig.gmem)*100,format='(f4.1)')+'%)'+$
          string(count_stars,format='(i7)')+' ('+$
          string(float(count_stars)/n_elements(mt.gcIndOrig.stars)*100,format='(f4.1)')+'%)'+' || cum '+$
          string(accCount.gal[k],format='(i7)')+' ('+$
          string(float(accCount.gal[k])/n_elements(mt.gcIndOrig.gal)*100,format='(f4.1)')+'%) '+$
          string(accCount.gmem[k],format='(i7)')+' ('+$
          string(float(accCount.gmem[k])/n_elements(mt.gcIndOrig.gmem)*100,format='(f4.1)')+'%)'+$
          string(accCount.stars[k],format='(i7)')+' ('+$
          string(float(accCount.stars[k])/n_elements(mt.gcIndOrig.stars)*100,format='(f4.1)')+'%)'
        
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
                     
            ; fix: interpolating to the m-1 snapshot could be untracked, in which case mt.hVirTemp=0
            if tvir[0] eq 0.0 or abs(tvir[0]-tvir[1]) gt 0.5 then tvir[0]=tvir[1]
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
                     
            ; fix: interpolating to the m-1 snapshot could be untracked, in which case mt.hVirTemp=0
            if tvir[0] eq 0.0 or abs(tvir[0]-tvir[1]) gt 0.5 then tvir[0]=tvir[1]
            tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
            r.accHaloTvir_gmem[gmem_w[i]] = tvir            
          endif
        endfor
        
        if countStars gt 0 then begin
          for i=0,count_stars-1 do begin
            ; stars_w indexes stars_pri (e.g. galcat_stars_ind) which indexes r.x_stars
            curInd = galcat_stars_ind[stars_w[i]]
            
            radii = [ prevRad.stars[curInd],stars_pri[stars_w[i]] ]
            time = interpol(times,radii,rVirFac) ; lerp time to r/rvir=rVirFac
            r.accTime_stars[k,curInd] = time
            if k eq 0 then begin
              tvir = [ mt.hVirTemp[mt.maxSnap-m-1,mt.gcIndOrig.stars[curInd]], $
                       mt.hVirTemp[mt.maxSnap-m,mt.gcIndOrig.stars[curInd]] ]
                       
              ; fix: interpolating to the m-1 snapshot could be untracked, in which case mt.hVirTemp=0
              if tvir[0] eq 0.0 or abs(tvir[0]-tvir[1]) gt 0.5 then tvir[0]=tvir[1]
              tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
              r.accHaloTvir_stars[curInd] = tvir
            endif
          endfor
        endif
        
        ; if we are on the first snapshot, override accretion times with -1 to indicate always outside rad
        if m eq mt.maxSnap then begin
          r.accTime_gal[k,gal_w] = -1
          r.accTime_gmem[k,gmem_w] = -1
          ; stars: no stars will be found in gas at first snapshot
        endif else begin
          ; otherwise, update counters for the number of particles we have found the accretion times of
          accCount.gal[k]  += count_gal
          accCount.gmem[k] += count_gmem
          if countStars gt 0 then accCount.stars[k] += count_stars
        endelse
        
        ; update mask for particles we no longer search for
        accMask.gal[k,gal_w]   = 1B
        accMask.gmem[k,gmem_w] = 1B
        if countStars gt 0 then accMask.stars[k,galcat_stars_ind[stars_w]] = 1B

      endforeach
      
      ; adaptive: update mask, mark all children of halos whose tracking ends at this snapshot
      w = where(gasMinSnap.gal eq sP.snap,count)
      if count gt 0 then accMask.gal[*,w] = 1B
      w = where(gasMinSnap.gmem eq sP.snap,count)
      if count gt 0 then accMask.gmem[*,w] = 1B
      w = where(gasMinSnap.stars eq sP.snap,count)
      if count gt 0 then accMask.stars[*,w] = 1B
      
      ; store current radius of particles
      prevRad.gal  = gal_pri
      prevRad.gmem = gmem_pri
      if countStars gt 0 then prevRad.stars[galcat_stars_ind] = stars_pri
      
      prevTime = h.time
      
      ; free some memory for next load
      gal_w    = !NULL
      gmem_w   = !NULL
      gal_pri  = !NULL
      gmem_pri = !NULL
      if countStars gt 0 then begin
        stars_w   = !NULL
        stars_pri = !NULL
      endif
    endfor
    
    ; print final results
    print,'[-] (rhot) found accretion times for ['+$
      str(accCount.galRT)+' of '+str(n_elements(mt.galcatSub.gal))+$
      '] gal, ['+str(accCount.starsRT)+' of '+str(n_elements(mt.galcatSub.stars))+'] stars'
    foreach rVirFac,sP.rVirFacs,k do $
    print,'['+str(k)+'] r='+string(sP.rVirFacs[k],format='(f4.2)')+' found accretion times for ['+$
      str(accCount.gal[k])+' of '+str(n_elements(mt.galcatSub.gal))+$
      '] gal, ['+str(accCount.gmem[k])+' of '+str(n_elements(mt.galcatSub.gmem))+'] gmem'+$
      ', ['+str(accCount.stars[k])+' of '+str(n_elements(mt.galcatSub.stars))+'] stars'
    
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
      
      ; locate tracer children at starting snapshot
      gcIndOrigTr = mergerTreeRepParentIDs(mt=mt,galcat=galcat,sP=sP,/compactMtS,$
          trids_gal=galcat_gal_trids,trids_gmem=galcat_gmem_trids,trids_stars=galcat_stars_trids,$
          gc_gal_cc=galcat_gal_cc,gc_gmem_cc=galcat_gmem_cc,gc_stars_cc=galcat_stars_cc)
          
      galcat = !NULL ; not used past this point    
          
      ; replicate hMinSnap for each child gas/star tracer
      trMinSnap = { gal   : mt.hMinSnap[gcIndOrigTr.gal]  ,$
                    gmem  : mt.hMinSnap[gcIndOrigTr.gmem] ,$
                    stars : mt.hMinSnap[gcIndOrigTr.stars] }  
          
      ; store the r/rvir of each at the previous snapshot for interpolation, and as a monotonic check
      prevRad = { gal   : fltarr(n_elements(galcat_gal_trids))  ,$
                  gmem  : fltarr(n_elements(galcat_gmem_trids)) ,$
                  stars : fltarr(n_elements(galcat_stars_trids))  }
      
      accMask = { gal   : bytarr(nVirFacs,n_elements(galcat_gal_trids))  ,$
                  gmem  : bytarr(nVirFacs,n_elements(galcat_gmem_trids)) ,$
                  stars : bytarr(nVirFacs,n_elements(galcat_stars_trids)) }
      
      ; store the main arrays as a structure so we can write them directly for all tracers
      r = {accTime_gal        : fltarr(nVirFacs,n_elements(galcat_gal_trids))-1   ,$
           accTime_gmem       : fltarr(nVirFacs,n_elements(galcat_gmem_trids))-1  ,$
           accTime_stars      : fltarr(nVirFacs,n_elements(galcat_stars_trids))-1 ,$
           accTimeRT_gal      : fltarr(n_elements(galcat_gal_trids))-1            ,$
           accTimeRT_stars    : fltarr(n_elements(galcat_stars_trids))-1          ,$
           accHaloTvir_gal    : fltarr(n_elements(galcat_gal_trids))              ,$
           accHaloTvir_gmem   : fltarr(n_elements(galcat_gmem_trids))             ,$
           accHaloTvir_stars  : fltarr(n_elements(galcat_stars_trids))            ,$
           child_counts_gal   : galcat_gal_cc                                     ,$
           child_counts_gmem  : galcat_gmem_cc                                    ,$
           child_counts_stars : galcat_stars_cc                                   ,$
           rVirFacs           : sP.rVirFacs                                        }
           
      galcat_gal_cc   = !NULL
      galcat_gmem_cc  = !NULL
      galcat_stars_cc = !NULL
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
        save,prevRad,trMinSnap,accMask,accCount,r,gcIndOrigTr,$
             galcat_gal_trids,galcat_gmem_trids,galcat_stars_ind,prevTime,m,filename=resFilename
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
      
      match,galcat_stars_trids,tr_ids,galcat_stars_ind,trids_stars_ind,count=countStars,/sort
      trids_stars_ind  = trids_stars_ind[sort(galcat_stars_ind)]
      galcat_stars_ind = galcat_stars_ind[sort(galcat_stars_ind)]
      
      tr_ids     = !NULL
      galcat_ind = !NULL
      
      if countGal ne n_elements(galcat_gal_trids) or countGmem ne n_elements(galcat_gmem_trids) then message,'Error: Counts'
      
      ; load tracer parents to match to gas
      tr_parids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentid')
      tr_parids_gal   = tr_parids[trids_gal_ind]
      tr_parids_gmem  = tr_parids[trids_gmem_ind]
      tr_parids_stars = tr_parids[trids_stars_ind]
      tr_parids = !NULL
      
      ; load gas IDs and convert tracer parent IDs -> indices
      gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
      gasIDMap = getIDIndexMap(gas_ids,minid=minid)
      
      ; match star tracer parents to gas ids (keep only those with gas parents at this snapshot)
      match,gas_ids,tr_parids_stars,ind1,ind2,count=countStars ; override countStars
      ind2 = ind2[sort(ind2)]
      
      if countStars gt 0 then begin
        tr_parids_stars  = tr_parids_stars[ind2]
        galcat_stars_ind = galcat_stars_ind[ind2]
      endif
      
      ind1 = !NULL
      ind2 = !NULL
      gas_ids = !NULL
      
      tr_parids_gal  = gasIDMap[tr_parids_gal-minid]  ; convert ID->index
      tr_parids_gmem = gasIDMap[tr_parids_gmem-minid]
      if countStars gt 0 then tr_parids_stars = gasIDMap[tr_parids_stars-minid]
      gasIDMap = !NULL
      
      ; load density,temp to apply galaxy cut
      u = loadSnapshotSubset(sP=sP,partType='gas',field='u')
      u_gal = u[tr_parids_gal]
      if countStars gt 0 then u_stars = u[tr_parids_stars]
      u = !NULL
      
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
      nelec_gal = nelec[tr_parids_gal]
      if countStars gt 0 then nelec_stars = nelec[tr_parids_stars]
      nelec = !NULL
      
      temp_gal = convertUtoTemp(u_gal,nelec_gal)
      u_gal = !NULL
      nelec_gal = !NULL
      
      if countStars gt 0 then begin
        temp_stars = convertUtoTemp(u_stars,nelec_stars)
        u_stars = !NULL
        nelec_stars = !NULL
      endif
      
      ; load rho of gas and make galaxy (rho,temp) plane cut
      dens = loadSnapshotSubset(sP=sP,partType='gas',field='density')
      dens_gal = dens[tr_parids_gal]
      if countStars gt 0 then dens_stars = dens[tr_parids_stars]
      dens = !NULL
      
      ; scale Torrey+ (2011) galaxy cut to physical density
      scalefac = snapNumToRedshift(sP=sP,/time) ; time flag gives simulation time = scale factor
      a3inv = 1.0 / (scalefac*scalefac*scalefac)
      dens_gal *= a3inv
      if countStars gt 0 then dens_stars *= a3inv
      
      ; mark any galaxy gas failing cut as accreted at this snapshot (if not previously marked)
      wGalCut = where(alog10(temp_gal) - sP.galcut_rho * alog10(dens_gal) ge sP.galcut_T and $
                      r.accTimeRT_gal eq -1,countGalCut)
               
      if countGalCut gt 0 then r.accTimeRT_gal[wGalCut] = h.time
      accCount.galRT += countGalCut

      wGalCut  = !NULL
      dens_gal = !NULL
      temp_gal = !NULL
      
      ; mark any gas hosting stellar tracers as accreted as above
      if countStars gt 0 then begin
        wStarsCut = where(alog10(temp_stars) - sP.galcut_rho * alog10(dens_stars) ge sP.galcut_T and $
                        r.accTimeRT_stars eq -1,countStarsCut)
                 
        if countStarsCut gt 0 then r.accTimeRT_stars[galcat_stars_ind[wStarsCut]] = h.time
        accCount.starsRT += countStarsCut
        
        wStarsCut  = !NULL
        dens_stars = !NULL
        temp_stars = !NULL
      endif else begin
        countStarsCut = 0
      endelse
      
      print,' ['+string(m,format='(i3)')+'] [-] (rhot)'+$
        ' accreted now counts '+string(countGalCut,format='(i7)')+' ('+$
        string(float(countGalCut)/n_elements(galcat_gal_trids)*100,format='(f4.1)')+'%)   stars '+$
        string(countStarsCut,format='(i7)')+' ('+$
        string(float(countStarsCut)/n_elements(galcat_stars_trids)*100,format='(f4.1)')+'%) || cum '+$
        string(accCount.galRT,format='(i7)')+' ('+$
        string(float(accCount.galRT)/n_elements(galcat_gal_trids)*100,format='(f4.1)')+'%)   stars '+$
        string(accCount.starsRT,format='(i7)')+' ('+$
        string(float(accCount.starsRT)/n_elements(galcat_stars_trids)*100,format='(f4.1)')+'%)'
      
      ; load gas positions and convert to tracer positions
      gas_pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
      tr_pos_gal  = gas_pos[*,tr_parids_gal]
      tr_pos_gmem = gas_pos[*,tr_parids_gmem]
      if countStars gt 0 then tr_pos_stars = gas_pos[*,tr_parids_stars]
      
      gas_pos = !NULL
      tr_parids_gal  = !NULL
      tr_parids_gmem = !NULL
      tr_parids_stars = !NULL

      ; calculate current distance of gas particle from smoothed halo center position for galaxy members
      gal_pri  = periodicDists(reform(mt.hPos[mt.maxSnap-m,*,gcIndOrigTr.gal]),tr_pos_gal,sP=sP)
      gal_pri /= mt.hVirRad[mt.maxSnap-m,gcIndOrigTr.gal]
      tr_pos_gal  = !NULL
      
      ; for group members
      gmem_pri = periodicDists(reform(mt.hPos[mt.maxSnap-m,*,gcIndOrigTr.gmem]),tr_pos_gmem,sP=sP)
      gmem_pri /= mt.hVirRad[mt.maxSnap-m,gcIndOrigTr.gmem]
      tr_pos_gmem = !NULL
      
      ; stars
      if countStars gt 0 then begin
        stars_pri = periodicDists($
          reform(mt.hPos[mt.maxSnap-m,*,gcIndOrigTr.stars[galcat_stars_ind]]),tr_pos_stars,sP=sP)
        stars_pri /= mt.hVirRad[mt.maxSnap-m,gcIndOrigTr.stars[galcat_stars_ind]]
        pos_stars = !NULL
      endif
      
      ; loop over each critical radius
      foreach rVirFac,sP.rVirFacs,k do begin
        ; for particles who are still within r_vir, check if they have passed beyond
        gal_w  = where(gal_pri ge rVirFac and prevRad.gal lt rVirFac and accMask.gal[k,*] eq 0B,count_gal)
        gmem_w = where(gmem_pri ge rVirFac and prevRad.gmem lt rVirFac and accMask.gmem[k,*] eq 0B,count_gmem)
        
        if countStars gt 0 then begin
          stars_w = where(stars_pri ge rVirFac and prevRad.stars[galcat_stars_ind] lt rVirFac and $
                          accMask.stars[k,galcat_stars_ind] eq 0B,count_stars)
        endif else begin
          count_stars = 0 & stars_pri = 0
        endelse
        
        print,' ['+string(m,format='(i3)')+'] ['+str(k)+'] r='+string(rVirFac,format='(f4.2)')+$
          ' accreted now counts '+string(count_gal,format='(i7)')+' ('+$
          string(float(count_gal)/n_elements(galcat_gal_trids)*100,format='(f4.1)')+'%) '+$
          string(count_gmem,format='(i7)')+' ('+$
          string(float(count_gmem)/n_elements(galcat_gmem_trids)*100,format='(f4.1)')+'%)'+$
          string(count_stars,format='(i7)')+' ('+$
          string(float(count_stars)/n_elements(galcat_stars_trids)*100,format='(f4.1)')+'%)'+' || cum '+$
          string(accCount.gal[k],format='(i7)')+' ('+$
          string(float(accCount.gal[k])/n_elements(galcat_gal_trids)*100,format='(f4.1)')+'%) '+$
          string(accCount.gmem[k],format='(i7)')+' ('+$
          string(float(accCount.gmem[k])/n_elements(galcat_gmem_trids)*100,format='(f4.1)')+'%)'+$
          string(accCount.stars[k],format='(i7)')+' ('+$
          string(float(accCount.stars[k])/n_elements(galcat_stars_trids)*100,format='(f4.1)')+'%)'
        
        ; interpolate these (time,radii) to find time crossing the virial radius
        times = [prevTime,h.time]
        
        for i=0,count_gal-1 do begin
          radii = [ prevRad.gal[gal_w[i]],gal_pri[gal_w[i]] ]
          time = interpol(times,radii,rVirFac) ; lerp time to r/rvir=rVirFac
          r.accTime_gal[k,gal_w[i]] = time
          
          if k eq 0 then begin
            tvir = [ mt.hVirTemp[mt.maxSnap-m-1,gcIndOrigTr.gal[gal_w[i]]], $
                     mt.hVirTemp[mt.maxSnap-m,gcIndOrigTr.gal[gal_w[i]]] ]
            if tvir[0] eq 0.0 or abs(tvir[0]-tvir[1]) gt 0.5 then tvir[0]=tvir[1]
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
            if tvir[0] eq 0.0 or abs(tvir[0]-tvir[1]) gt 0.5 then tvir[0]=tvir[1]
            tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
            r.accHaloTvir_gmem[gmem_w[i]] = tvir
          endif
        endfor
        
        if countStars gt 0 then begin
          for i=0,count_stars-1 do begin
            ; stars_w indexes stars_pri (e.g. galcat_stars_ind) which indexes r.x_stars
            curInd = galcat_stars_ind[stars_w[i]]
            
            radii = [ prevRad.stars[curInd],stars_pri[stars_w[i]] ]
            time = interpol(times,radii,rVirFac) ; lerp time to r/rvir=rVirFac
            r.accTime_stars[k,curInd] = time
            
            if k eq 0 then begin
              tvir = [ mt.hVirTemp[mt.maxSnap-m-1,gcIndOrigTr.stars[curInd]], $
                       mt.hVirTemp[mt.maxSnap-m,gcIndOrigTr.stars[curInd]] ]
              if tvir[0] eq 0.0 or abs(tvir[0]-tvir[1]) gt 0.5 then tvir[0]=tvir[1]
              tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
              r.accHaloTvir_stars[curInd] = tvir
            endif
          endfor
        endif
        
        ; if we are on the first snapshot, override accretion times with -1 to indicate always outside rad
        if m eq mt.maxSnap then begin
          r.accTime_gal[k,gal_w] = -1
          r.accTime_gmem[k,gmem_w] = -1
          ; stars: no stars will be found in gas at first snapshot
        endif else begin
          ; otherwise, update counters for the number of particles we have found the accretion times of
          accCount.gal[k]  += count_gal
          accCount.gmem[k] += count_gmem
          if countStars gt 0 then accCount.stars[k] += count_stars
        endelse
        
        ; update mask for particles we no longer search for
        accMask.gal[k,gal_w]   = 1B
        accMask.gmem[k,gmem_w] = 1B
        if countStars gt 0 then accMask.stars[k,galcat_stars_ind[stars_w]] = 1B
      endforeach
      
      ; adaptive: update mask, mark all children of halos whose tracking ends at this snapshot
      w = where(trMinSnap.gal eq sP.snap,count)
      if count gt 0 then accMask.gal[*,w] = 1B
      w = where(trMinSnap.gmem eq sP.snap,count)
      if count gt 0 then accMask.gmem[*,w] = 1B
      w = where(trMinSnap.stars eq sP.snap,count)
      if count gt 0 then accMask.stars[*,w] = 1B
      
      ; store current radius of particles
      prevRad.gal  = gal_pri
      prevRad.gmem = gmem_pri
      if countStars gt 0 then prevRad.stars[galcat_stars_ind] = stars_pri
      
      prevTime = h.time
      
      ; free some memory for next load
      gal_w    = !NULL
      gmem_w   = !NULL
      gal_pri  = !NULL
      gmem_pri = !NULL
      if countStars gt 0 then begin
        stars_w   = !NULL
        stars_pri = !NULL
      endif
    endfor
    
    ; print final results
    print,'[-] (rhot) found accretion times for ['+$
      str(accCount.galRT)+' of '+str(n_elements(galcat_gal_trids))+$
      '] gal, ['+str(accCount.starsRT)+' of '+str(n_elements(galcat_stars_trids))+'] stars'
    foreach rVirFac,sP.rVirFacs,k do $
    print,'['+str(k)+'] r='+string(sP.rVirFacs[k],format='(f4.2)')+' found accretion times for ['+$
      str(accCount.gal[k])+' of '+str(n_elements(galcat_gal_trids))+$
      '] gal, ['+str(accCount.gmem[k])+' of '+str(n_elements(galcat_gmem_trids))+'] gmem'+$
      ', ['+str(accCount.stars[k])+' of '+str(n_elements(galcat_stars_trids))+'] stars'
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
  endif
  
  
  ; VELOCITY TRACERS case - will be similar to above since there could be multiple
  ; ---------------------
  if sP.trMCPerCell eq -1 then begin
  
    if ~file_test(resFilename) then begin
      ; no restart
      print,'Calculating new accretion time using ( TracerVEL ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
  
      print,'NOTE: no stars (or accTimeRT) for tracerVEL'  
      
      ; locate tracer children at starting snapshot
      gcIndOrigTr = mergerTreeRepParentIDs(mt=mt,galcat=galcat,sP=sP,/compactMtS,$
                      trids_gal=galcat_gal_trids,trids_gmem=galcat_gmem_trids,trids_stars=galcat_stars_trids)
          
      galcat = !NULL ; not used past this point
      
      ; replicate hMinSnap for each child gas element
      trMinSnap = { gal : mt.hMinSnap[gcIndOrigTr.gal] , gmem : mt.hMinSnap[gcIndOrigTr.gmem] }  
      
      ; store the r/rvir of each at the previous snapshot for interpolation, and as a monotonic check
      prevRad = { gal  : fltarr(n_elements(galcat_gal_trids)) ,$
                  gmem : fltarr(n_elements(galcat_gmem_trids)) }
      
      accMask = { gal  : bytarr(nVirFacs,n_elements(galcat_gal_trids)), $
                  gmem : bytarr(nVirFacs,n_elements(galcat_gmem_trids)) }
  
      ; store the main arrays as a structure so we can write them directly for all tracers
      r = {accTime_gal       : fltarr(nVirFacs,n_elements(galcat_gal_trids))-1  ,$
           accTimeRT_gal     : fltarr(n_elements(galcat_gal_trids))-1           ,$
           accTime_gmem      : fltarr(nVirFacs,n_elements(galcat_gmem_trids))-1 ,$
           accHaloTvir_gal   : fltarr(n_elements(galcat_gal_trids))             ,$
           accHaloTvir_gmem  : fltarr(n_elements(galcat_gmem_trids))            ,$
           child_counts_gal  : galcat_gal_cc                                    ,$
           child_counts_gmem : galcat_gmem_cc                                   ,$
           rVirFacs          : sP.rVirFacs                                          }
 
    endif else begin
      ; restart
      if ~keyword_set(restart) then message,'Error: Restart file exists but restart not requested.'
      restore,resFilename,/verbose
      snapRange[0] = m
    endelse

    for m=snapRange[0],snapRange[1],-1 do begin
      sP.snap = m
      print,m
      
      ; save restart?
      if m mod 10 eq 0 and m lt snapRange[0] and keyword_set(restart) then begin
        print,' --- Writing restart! ---'
        save,prevRad,trMinSnap,accMask,accCount,r,gcIndOrigTr,$
             galcat_gal_trids,galcat_gmem_trids,prevTime,m,filename=resFilename
        print,' --- Done! ---'
      endif
      
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
      
      ; no (rho,temp) cut (would need to load tracer parents which isn't currently done)
      
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
      foreach rVirFac,sP.rVirFacs,k do begin
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
            if tvir[0] eq 0.0 or abs(tvir[0]-tvir[1]) gt 0.5 then tvir[0]=tvir[1]
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
            if tvir[0] eq 0.0 or abs(tvir[0]-tvir[1]) gt 0.5 then tvir[0]=tvir[1]
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
      
      ; adaptive: update mask, mark all children of halos whose tracking ends at this snapshot
      w = where(trMinSnap.gal eq sP.snap,count)
      if count gt 0 then accMask.gal[*,w] = 1B
      w = where(trMinSnap.gmem eq sP.snap,count)
      if count gt 0 then accMask.gmem[*,w] = 1B
      
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
              
    foreach rVirFac,sP.rVirFacs,k do $
    print,'['+str(k)+'] r='+string(sP.rVirFacs[k],format='(f4.2)')+' found accretion times for ['+$
      str(accCount.gal[k])+' of '+str(n_elements(galcat_gal_trids))+$
      '] gal, ['+str(accCount.gmem[k])+' of '+str(n_elements(galcat_gmem_trids))+'] gmem'
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
                
  endif
  
end
