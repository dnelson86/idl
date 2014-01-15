; accretionTimes.pro
; gas accretion project - past radial history of gas elements (virial radius crossing)
; dnelson dec.2013

; -----------------------------------------------------------------------------------------------------
; accretionTimes(): for each gas particle/tracer, starting at some redshift, track backwards in time
;                   with respect to the tracked parent halos (using mergerTree) and determine the
;                   time when the particle radius = some fraction of the virial radius (and record the 
;                   virial temp of the parent halo at the rvir crossing time). 
; -----------------------------------------------------------------------------------------------------

function accretionTimes, sP=sP, restart=restart

  forward_function cosmoTracerChildren, cosmoTracerVelParents
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; first, walk back through the merger tree and find primary subhalos with good parent histories
  nVirFacs = n_elements(sP.rVirFacs) + 2 ; two extra for first 0.15/1.0 rvir crossings
  mt = mergerTreeSubset(sP=sP,/verbose)
  
  snapStep = -1 ; process every snapshot
  ;if sP.snapRange[1] ge 300 then snapStep = -2 ; process every other snap for gadget/tracer runs
  
  snapRange = [mt.maxSnap,mt.minSnap]
  
  ; set saveFilename and check for existence
  saveFilename = sP.derivPath + 'accTimes.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(mt.maxSnap)+'-'+str(mt.minSnap)+'.sav'  
  
  if file_test(saveFilename) then begin
    restore, saveFilename
    return, r
  endif
  
  resFilename = sP.derivPath + 'accTimes.'+sP.saveTag+'.restart.'+sP.savPrefix+str(sP.res)+'.'+$
                str(mt.maxSnap)+'-'+str(mt.minSnap)+'.sav'  
  
  ; load galaxy/group member catalogs at zMin for gas ids to search for
  origSnap = sP.snap
  galcat = galaxyCat(sP=sP)
  
  ; count of how many particles/tracers we tracked through r_vir
  accCount   = ulonarr(nVirFacs)
  outCount   = ulonarr(nVirFacs)
  accCountRT = 0L
  outCountRT = 0L
  
  prevTime = 0 ; scale factor at previous snapshot

  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    print,'Calculating new accretion time using ( SPH Particles ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
      
    if ~file_test(resFilename) then begin ; no restart
      ; replicate hMinSnap for each child gas element
      gasMinSnap = mt.hMinSnap[mt.gcIndOrig]
  
      ; store the r/rvir of each at the previous snapshot for interpolation, and as a monotonic check
      prevRad = fltarr(galcat.countTot)
      accMask = bytarr(galcat.countTot)

      ; store the main arrays as a structure so we can write them directly
      r = {accTime       : fltarr(nVirFacs,galcat.countTot)-1   ,$
           accTimeRT     : fltarr(galcat.countTot)-1            ,$
           outTime       : fltarr(nVirFacs,galcat.countTot)-1   ,$
           outTimeRT     : fltarr(galcat.countTot)-1            ,$
           accHaloTvir   : fltarr(galcat.countTot)              ,$
           rVirFacs      : sP.rVirFacs                           }

    endif else begin
      ; restart
      if ~keyword_set(restart) then message,'Error: Restart file exists but restart not requested.'
      restore,resFilename,/verbose
      snapRange[0] = m
    endelse
    
    for m=snapRange[0],snapRange[1],snapStep do begin
      sP.snap = m
      
      ; save restart?
      if m mod 20 eq 0 and m lt snapRange[0] and keyword_set(restart) then begin
        print,' --- Writing restart! ---'
        save,prevRad,gasMinSnap,accMask,r,prevTime,accCount,outCount,$
          accCountRT,outCountRT,m,filename=resFilename
        print,' --- Done! ---'
        ;exit,status=33 ; requeue
      endif
      
      ; load gas ids and match to catalog
      h = loadSnapshotHeader(sP=sP)
            
      ; --- for each each possible parent particle type, match child tracers and save times ---
      parPartTypes = ['gas','stars']
      
      foreach partType,parPartTypes do begin
      
        par_ids = loadSnapshotSubset(sP=sP,partType=partType,field='ids')
            
        ; note: tr_parids are NOT UNIQUE, use a value_locate approach (not match)
        ; they are in fact unique for each partType for this SPH case, but keep same approach
        sort_inds = calcSort(par_ids)
        par_ids_sorted = par_ids[sort_inds]
        
        ; locate
        par_ind = value_locate(par_ids_sorted,galcat.ids) ; indices to par_ids_sorted
        par_ind = sort_inds[par_ind>0] ; indices to par_ids (>0 removes -1 entries, which are removed next line)
        galcat_ind_inPar = where(par_ids[par_ind] eq galcat.ids,count_inPar) ; verify we actually matched the ID
        par_ind = par_ind[galcat_ind_inPar] ; indices of matched parents
        
        sort_inds = !NULL
        par_ids_sorted = !NULL
              
        ; for tracers with gas parents: apply galaxy cut in (rho,T) plane
        if partType eq 'gas' then begin
          ; load density,temp
          u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
          nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
          temp  = convertUtoTemp(u,nelec)
          u     = !NULL
          nelec = !NULL
          
          if count_inPar gt 0 then temp = temp[par_ind]
      
          ; scale Torrey+ (2012) galaxy cut to physical density
          scalefac = snapNumToRedshift(sP=sP,/time) ; time flag gives simulation time = scale factor
          a3inv = 1.0 / (scalefac*scalefac*scalefac)
      
          dens = loadSnapshotSubset(sP=sP,partType='gas',field='density')
          if count_inPar gt 0 then dens = dens[par_ind] * a3inv
      
          ; mark any galaxy gas failing cut as accreted at this snapshot, if not previously
          ; marked due to RT cut (note: for gmem accTimeRT is flagged at mt.maxSnap by definition)
          if count_inPar gt 0 then begin
            wCut = where(alog10(temp) - sP.galcut_rho * alog10(dens) ge sP.galcut_T and $
                        r.accTimeRT[galcat_ind_inPar] eq -1,countRTCut)
      
            if countRTCut gt 0 then r.accTimeRT[galcat_ind_inPar[wCut]] = h.time
            accCountRT += countRTCut
            
            ; similarly, look for outflow causing the tracer to satisfy the RT cut
            wCut = where(alog10(temp) - sP.galcut_rho * alog10(dens) lt sP.galcut_T and $
                        r.outTimeRT[galcat_ind_inPar] eq -1,countRTCutOut)
               
            if countRTCutOut gt 0 then r.outTimeRT[galcat_ind_inPar[wCut]] = h.time
            outCountRT += countRTCutOut
            
            wCut  = !NULL
            dens = !NULL
            temp = !NULL
          endif
          
          print,' ['+string(m,format='(i3)')+'] [-] (rhot)'+$
            ' '+strpad(partType,5)+' accNowCounts '+string(countRTCut,format='(i8)')+' ('+$
            string(float(countRTCut)/n_elements(galcat.ids)*100,format='(f4.1)')+'%)'+$
            ' outNowCounts '+string(countRTCutOut,format='(i8)')+' ('+$
            string(float(countRTCutOut)/n_elements(galcat.ids)*100,format='(f4.1)')+'%) || cumIn '+$
            string(accCountRT,format='(i8)')+' ('+$
            string(float(accCountRT)/n_elements(galcat.ids)*100,format='(f4.1)')+'%)'+$
            ' cumOut '+string(outCountRT,format='(i8)')+' ('+$
            string(float(outCountRT)/n_elements(galcat.ids)*100,format='(f4.1)')+'%)'
          
        endif ; partType eq 'gas'
            
        ; load parent positions and convert to tracer positions
        par_pos = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
        if count_inPar gt 0 then child_pos = par_pos[*,par_ind]
      
        par_pos = !NULL
        par_ind = !NULL

        ; calculate current distance of parent from smoothed halo center position for galaxy members
        if count_inPar gt 0 then begin
          rad_pri = periodicDists( reform(mt.hPos[mt.maxSnap-m,*,mt.gcIndOrig[galcat_ind_inPar]]),child_pos,sP=sP )
            
          rad_pri /= mt.hVirRad[mt.maxSnap-m,mt.gcIndOrig[galcat_ind_inPar]]
          child_pos  = !NULL
        endif
        
        ; loop over each target radius
        foreach rVirFac,[sP.rVirFacs,sP.radcut_rvir,1.0],k do begin
          count_rad = 0
          count_out = 0
        
          if k eq n_elements(sP.rVirFacs) then begin
            ; for the last two iterations, take 0.15/1.0rvir and do not use accMask
            ; thereby recording the earliest/highest redshift crossing
            if count_inPar gt 0 then begin
              rad_w = where(rad_pri ge rVirFac and prevRad[galcat_ind_inPar] lt rVirFac and $
                            accMask[galcat_ind_inPar] eq 0B,count_rad)
              out_w = where(rad_pri le rVirFac and prevRad[galcat_ind_inPar] gt rVirFac and $
                            accMask[galcat_ind_inPar] eq 0B,count_out)
            endif
          endif else begin
            ; take rVirFac and use accMask (skip if past the end of halo tracking)
            if count_inPar gt 0 then begin
              rad_w = where(rad_pri ge rVirFac and prevRad[galcat_ind_inPar] lt rVirFac and $
                            accMask[galcat_ind_inPar] eq 0B and r.accTime[k,galcat_ind_inPar] eq -1,count_rad)
              out_w = where(rad_pri le rVirFac and prevRad[galcat_ind_inPar] gt rVirFac and $
                            accMask[galcat_ind_inPar] eq 0B and r.outTime[k,galcat_ind_inPar] eq -1,count_out)          
            endif
          endelse
        
          print,' ['+string(m,format='(i3)')+'] ['+str(k)+'] r='+string(rVirFac,format='(f4.2)')+$
            ' '+strpad(partType,5)+' accNowCounts '+string(count_rad,format='(i7)')+' ('+$
            string(float(count_rad)/n_elements(galcat.ids)*100,format='(f4.1)')+'%)'+$
            ' outNowCounts '+string(count_out,format='(i7)')+' ('+$
            string(float(count_out)/n_elements(galcat.ids)*100,format='(f4.1)')+'%) || cumIn '+$
            string(accCount[k],format='(i7)')+' ('+$
            string(float(accCount[k])/n_elements(galcat.ids)*100,format='(f4.1)')+'%)'+$
            ' cumOut '+string(outCount[k],format='(i7)')+' ('+$
            string(float(outCount[k])/n_elements(galcat.ids)*100,format='(f4.1)')+'%)'
        
          ; interpolate these (time,radii) to find time crossing the radius
          times = [prevTime,h.time]
        
          if count_inPar gt 0 then begin
            ; accretion/inflow
            for i=0,count_rad-1 do begin
              curInd = galcat_ind_inPar[rad_w[i]]
          
              radii = [ prevRad[curInd],rad_pri[rad_w[i]] ]
              time = interpol(times,radii,rVirFac) ; lerp time to r/rvir=rVirFac
              r.accTime[k,curInd] = time
          
              ; record Tvir at the first rvir crossing
              if k eq n_elements(sP.rVirFacs) then begin
                tvir = [ mt.hVirTemp[mt.maxSnap-m-1,mt.gcIndOrig[curInd]], $
                         mt.hVirTemp[mt.maxSnap-m,mt.gcIndOrig[curInd]] ]
                if tvir[0] eq 0.0 or abs(tvir[0]-tvir[1]) gt 0.5 then tvir[0]=tvir[1]
                tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
                r.accHaloTvir[curInd] = tvir
              endif
            endfor
            
            ; outflow
            for i=0,count_out-1 do begin
              curInd = galcat_ind_inPar[out_w[i]]
          
              radii = [ prevRad[curInd],rad_pri[out_w[i]] ]
              time = interpol(times,radii,rVirFac) ; lerp time to r/rvir=rVirFac
              r.outTime[k,curInd] = time
          
            endfor
          endif
        
          ; if we are on the first snapshot, override accretion times with -1 to indicate always outside rad
          if m eq mt.maxSnap then begin
            if count_inPar gt 0 then begin
              r.accTime[k,galcat_ind_inPar[rad_w]] = -1
              r.outTime[k,galcat_ind_inPar[out_w]] = -1
            endif
          endif else begin
            ; otherwise, update counters for the number of particles we have found the accretion times of
            accCount[k] += count_rad
            outCount[k] += count_out
          endelse
                  
      endforeach ; rVirFac
      
      ; adaptive: update mask, mark all children of halos whose tracking ends at this snapshot
      w = where(gasMinSnap eq sP.snap,count)
      
      if count gt 0 then accMask[w] = 1B
      
      prevTime = h.time
      
      ; store current radius of particles
      if count_inPar gt 0 then prevRad[galcat_ind_inPar] = rad_pri
      
      ; free some memory for next load
      rad_w   = !NULL
      out_w   = !NULL
      rad_pri = !NULL
      endforeach ; partType
      
      ; DEBUG: check for odd results
      w = where(r.accTime[0,*] lt r.accTime[-1,*],count)
      if count gt 0 then message,'How can the latest 1rvir crossing be earlier than the earliest?'
      
    endfor ;m
    
    ; print final results
    print,'[-] (rhot) found accretion/outflow times for ['+$
      str(accCountRT)+' / '+str(outCountRT)+' of '+str(n_elements(galcat.ids))+']'
    foreach rVirFac,[sP.rVirFacs,1.0],k do $
      print,'['+str(k)+'] r='+string(rVirFac,format='(f4.2)')+' found accretion times for ['+$
        str(accCount[k])+' or '+string(float(accCount[k])/n_elements(galcat.ids)*100,format='(f4.1)')+'%], '+$
        ' outflow times for ['+$
        str(outCount[k])+' or '+string(float(outCount[k])/n_elements(galcat.ids)*100,format='(f4.1)')+'%]'
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
  endif
  
  ; MONTE CARLO TRACERS CASE - for all galaxy catalog members, track back all child tracers
  ; NOTE: both gas and star parents are considered at each snapshot, unlike with TRVEL or SPH
  ; ------------------------
  if sP.trMCPerCell gt 0 then begin
                  
    if ~file_test(resFilename) then begin
      ; no restart
      print,'Calculating new accretion time using ( TracerMC ) res = '+str(sP.res)+$
        ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
      
      ; replicate hMinSnap for each child gas/star tracer
      trMinSnap = mt.hMinSnap[mt.gcIndOrigTrMC]
          
      ; store the r/rvir of each at the previous snapshot for interpolation, and as a monotonic check
      num_tracers = n_elements(galcat.trMC_ids)
      
      prevRad = fltarr(num_tracers)
      accMask = bytarr(num_tracers)
      
      ; store the main arrays as a structure so we can write them directly for all tracers
      r = {accTime        : fltarr(nVirFacs,num_tracers)-1   ,$
           accTimeRT      : fltarr(num_tracers)-1            ,$
           outTime        : fltarr(nVirFacs,num_tracers)-1   ,$
           outTimeRT      : fltarr(num_tracers)-1            ,$
           accHaloTvir    : fltarr(num_tracers)              ,$
           rVirFacs       : sP.rVirFacs                       }
           
    endif else begin
      ; restart
      if ~keyword_set(restart) then message,'Error: Restart file exists but restart not requested.'
      restore,resFilename,/verbose
      snapRange[0] = m
    endelse

    for m=snapRange[0],snapRange[1],snapStep do begin
      sP.snap = m
      
      ; save restart?
      if m mod 20 eq 0 and m lt snapRange[0] and keyword_set(restart) then begin
        print,' --- Writing restart! ---'
        save,prevRad,accMask,trMinSnap,accCount,accCountRT,r,prevTime,m,filename=resFilename
        print,' --- Done! ---'
        ;exit,status=33 ; requeue
      endif
      
      ; load tracer ids and match to child ids from zMin
      h = loadSnapshotHeader(sP=sP)
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
      
      idIndexMap = getIDIndexMap(tr_ids,minid=minid)
      trids_ind  = idIndexMap[galcat.trMC_ids - minid]
      idIndexMap = !NULL
      tr_ids     = !NULL
      
      ; load tracer parents IDs
      tr_parids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentid')
      tr_parids = tr_parids[trids_ind]
      trids_ind = !NULL
      
      ; --- for each each possible parent particle type, match child tracers and save times ---
      parPartTypes = ['gas','stars','BHs']
      
      foreach partType,parPartTypes do begin
        if partType eq 'BHs' and sP.gfmBHs eq 0 then continue
        
        ; load parent IDs and convert tracer parent IDs -> indices (for tracers with this partType parent now)
        par_ids = loadSnapshotSubset(sP=sP,partType=partType,field='ids')
      
        ; note: tr_parids are NOT UNIQUE, use a value_locate approach (not match)
        sort_inds = calcSort(par_ids)
        par_ids_sorted = par_ids[sort_inds]
        
        ; locate
        par_ind = value_locate(par_ids_sorted,tr_parids) ; indices to par_ids_sorted
        par_ind = sort_inds[par_ind>0] ; indices to par_ids (>0 removes -1 entries, which are removed next line)
        galcat_ind_inPar = where(par_ids[par_ind] eq tr_parids,count_inPar) ; verify we actually matched the ID
        par_ind = par_ind[galcat_ind_inPar] ; indices of matched parents
        
        sort_inds = !NULL
        par_ids_sorted = !NULL
        
        ; for tracers with gas parents: apply galaxy cut in (rho,T) plane
        if partType eq 'gas' then begin
          ; load density,temp
          u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
          nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
          temp  = convertUtoTemp(u,nelec)
          u     = !NULL
          nelec = !NULL
          
          if count_inPar gt 0 then temp = temp[par_ind]
      
          ; scale Torrey+ (2012) galaxy cut to physical density
          scalefac = snapNumToRedshift(sP=sP,/time) ; time flag gives simulation time = scale factor
          a3inv = 1.0 / (scalefac*scalefac*scalefac)
      
          dens = loadSnapshotSubset(sP=sP,partType='gas',field='density')
          if count_inPar gt 0 then dens = dens[par_ind] * a3inv
      
          ; mark any galaxy gas failing cut as accreted at this snapshot, if not previously
          ; marked due to RT cut (note: for gmem accTimeRT is flagged at mt.maxSnap by definition)
          if count_inPar gt 0 then begin
            wCut = where(alog10(temp) - sP.galcut_rho * alog10(dens) ge sP.galcut_T and $
                        r.accTimeRT[galcat_ind_inPar] eq -1,countRTCut)
               
            if countRTCut gt 0 then r.accTimeRT[galcat_ind_inPar[wCut]] = h.time
            accCountRT += countRTCut
            
            ; similarly, look for outflow causing the tracer to satisfy the RT cut (used for inter)
            wCut = where(alog10(temp) - sP.galcut_rho * alog10(dens) lt sP.galcut_T and $
                        r.outTimeRT[galcat_ind_inPar] eq -1,countRTCutOut)
               
            if countRTCutOut gt 0 then r.outTimeRT[galcat_ind_inPar[wCut]] = h.time
            outCountRT += countRTCutOut
            
            wCut  = !NULL
            dens = !NULL
            temp = !NULL
          endif
          
          print,' ['+string(m,format='(i3)')+'] [-] (rhot)'+$
            ' '+strpad(partType,5)+' accNowCounts '+string(countRTCut,format='(i9)')+' ('+$
            string(float(countRTCut)/n_elements(galcat.trMC_ids)*100,format='(f4.1)')+'%)'+$
            ' outNowCounts '+string(countRTCutOut,format='(i9)')+' ('+$
            string(float(countRTCutOut)/n_elements(galcat.trMC_ids)*100,format='(f4.1)')+'%) || cumIn '+$
            string(accCountRT,format='(i9)')+' ('+$
            string(float(accCountRT)/n_elements(galcat.trMC_ids)*100,format='(f4.1)')+'%)'+$
            ' cumOut '+string(outCountRT,format='(i9)')+' ('+$
            string(float(outCountRT)/n_elements(galcat.trMC_ids)*100,format='(f4.1)')+'%)'
          
        endif ; partType eq 'gas'
      
        ; load parent positions and convert to tracer positions
        par_pos = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
        if count_inPar gt 0 then tr_pos = par_pos[*,par_ind]
      
        par_pos = !NULL
        par_ind = !NULL

        ; calculate current distance of parent from smoothed halo center position for galaxy members
        if count_inPar gt 0 then begin
          rad_pri  = periodicDists( $
            reform(mt.hPos[mt.maxSnap-m,*,mt.gcIndOrigTrMC[galcat_ind_inPar]]),tr_pos,sP=sP)
            
          rad_pri /= mt.hVirRad[mt.maxSnap-m,mt.gcIndOrigTrMC[galcat_ind_inPar]]
          tr_pos  = !NULL
        endif
      
      ; loop over each critical radius
      foreach rVirFac,[sP.rVirFacs,sP.radcut_rvir,1.0],k do begin
        ; for particles who are still within r_vir, check if they have passed beyond
        count_rad = 0
        count_out = 0
        
        if k eq n_elements(sP.rVirFacs) then begin
          ; for the last two iterations, take 0.15/1.0rvir and do not use accMask
          ; thereby recording the earliest/highest redshift crossing
          if count_inPar gt 0 then begin
            rad_w = where(rad_pri ge rVirFac and prevRad[galcat_ind_inPar] lt rVirFac and $
                          accMask[galcat_ind_inPar] eq 0B,count_rad)
            out_w = where(rad_pri le rVirFac and prevRad[galcat_ind_inPar] gt rVirFac and $
                          accMask[galcat_ind_inPar] eq 0B,count_out)          
          endif
        endif else begin
          ; take rVirFac and use accMask (skip if past the end of halo tracking)
          if count_inPar gt 0 then begin
            rad_w = where(rad_pri ge rVirFac and prevRad[galcat_ind_inPar] lt rVirFac and $
                          accMask[galcat_ind_inPar] eq 0B and r.accTime[k,galcat_ind_inPar] eq -1,count_rad)
            out_w = where(rad_pri le rVirFac and prevRad[galcat_ind_inPar] gt rVirFac and $
                          accMask[galcat_ind_inPar] eq 0B and r.outTime[k,galcat_ind_inPar] eq -1,count_out)          
          endif
        endelse
                       
          print,' ['+string(m,format='(i3)')+'] ['+str(k)+'] r='+string(rVirFac,format='(f4.2)')+$
            ' '+strpad(partType,5)+' accNowCounts '+string(count_rad,format='(i7)')+' ('+$
            string(float(count_rad)/n_elements(galcat.trMC_ids)*100,format='(f4.1)')+'%)'+$
            ' outNowCounts '+string(count_out,format='(i7)')+' ('+$
            string(float(count_out)/n_elements(galcat.trMC_ids)*100,format='(f4.1)')+'%) || cumIn '+$
            string(accCount[k],format='(i7)')+' ('+$
            string(float(accCount[k])/n_elements(galcat.trMC_ids)*100,format='(f4.1)')+'%)'+$
            ' cumOut '+string(outCount[k],format='(i7)')+' ('+$
            string(float(outCount[k])/n_elements(galcat.trMC_ids)*100,format='(f4.1)')+'%)'
        
          ; interpolate these (time,radii) to find time crossing the radius
          times = [prevTime,h.time]
        
          if count_inPar gt 0 then begin
            ; accretion/inflow
            for i=0,count_rad-1 do begin
              curInd = galcat_ind_inPar[rad_w[i]]
          
              radii = [ prevRad[curInd],rad_pri[rad_w[i]] ]
              time = interpol(times,radii,rVirFac) ; lerp time to r/rvir=rVirFac
              r.accTime[k,curInd] = time
          
              ; record Tvir at the first rvir crossing
              if k eq n_elements(sP.rVirFacs) then begin
                tvir = [ mt.hVirTemp[mt.maxSnap-m-1,mt.gcIndOrigTrMC[curInd]], $
                         mt.hVirTemp[mt.maxSnap-m,mt.gcIndOrigTrMC[curInd]] ]
                if tvir[0] eq 0.0 or abs(tvir[0]-tvir[1]) gt 0.5 then tvir[0]=tvir[1]
                tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
                r.accHaloTvir[curInd] = tvir
              endif
            endfor
            
            ; outflow
            for i=0,count_out-1 do begin
              curInd = galcat_ind_inPar[out_w[i]]
          
              radii = [ prevRad[curInd],rad_pri[out_w[i]] ]
              time = interpol(times,radii,rVirFac) ; lerp time to r/rvir=rVirFac
              r.outTime[k,curInd] = time
          
            endfor
          endif
        
          ; if we are on the first snapshot, override accretion times with -1 to indicate always outside rad
          if m eq mt.maxSnap then begin
            if count_inPar gt 0 then begin
              r.accTime[k,galcat_ind_inPar[rad_w]] = -1
              r.outTime[k,galcat_ind_inPar[out_w]] = -1
            endif
          endif else begin
            ; otherwise, update counters for the number of particles we have found the accretion times of
            accCount[k] += count_rad
            outCount[k] += count_out
          endelse
        
        endforeach
      
        ; adaptive: update mask, mark all children of halos whose tracking ends at this snapshot
        w = where(trMinSnap eq sP.snap,count)
        if count gt 0 then accMask[w] = 1B
      
        ; store current radius of particles
        if count_inPar gt 0 then prevRad[galcat_ind_inPar] = rad_pri
      
        prevTime = h.time
      
        ; free some memory for next load
        rad_w    = !NULL
        out_w    = !NULL
        rad_pri  = !NULL
        
      endforeach ; partType
    endfor ; snapshot
    
    ; print final results
    print,'[-] (rhot) found accretion/outflow times for ['+$
      str(accCountRT)+' / '+str(outCountRT)+' of '+str(n_elements(galcat.trMC_ids))+']'
    foreach rVirFac,[sP.rVirFacs,1.0],k do $
      print,'['+str(k)+'] r='+string(rVirFac,format='(f4.2)')+' found accretion times for ['+$
        str(accCount[k])+' or '+string(float(accCount[k])/n_elements(galcat.trMC_ids)*100,format='(f4.1)')+'%], '+$
        ' outflow times for ['+$
        str(outCount[k])+' or '+string(float(outCount[k])/n_elements(galcat.trMC_ids)*100,format='(f4.1)')+'%]'
    
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
  
      ; replicate hMinSnap for each child gas element
      trMinSnap = mt.hMinSnap[mt.gcIndOrigTr]
      
      ; store the r/rvir of each at the previous snapshot for interpolation, and as a monotonic check
      num_tracers = n_elements(galcat.trVel_ids)
      
      prevRad = fltarr(num_tracers)
      accMask = bytarr(num_tracers)
  
      ; store the main arrays as a structure so we can write them directly for all tracers
      r = {accTime       : fltarr(nVirFacs,num_tracers)-1  ,$
           accTimeRT     : fltarr(num_tracers)-1           ,$
           outTime       : fltarr(nVirFacs,num_tracers)-1  ,$
           outTimeRT     : fltarr(num_tracers)-1           ,$
           accHaloTvir   : fltarr(num_tracers)             ,$
           rVirFacs      : sP.rVirFacs                      }
 
    endif else begin
      ; restart
      if ~keyword_set(restart) then message,'Error: Restart file exists but restart not requested.'
      restore,resFilename,/verbose
      snapRange[0] = m
    endelse

    for m=snapRange[0],snapRange[1],snapStep do begin
      sP.snap = m
      print,m
      
      ; save restart?
      if m mod 10 eq 0 and m lt snapRange[0] and keyword_set(restart) then begin
        print,' --- Writing restart! ---'
        save,prevRad,trMinSnap,accMask,accCount,accCountRT,outCount,outCountRT,r,$
             prevTime,m,filename=resFilename
        print,' --- Done! ---'
        exit,status=33 ; requeue
      endif
      
      ; load tracer ids and match to child ids from zMin
      h = loadSnapshotHeader(sP=sP)
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerVel',field='ids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs, need this if we want ids[ids_ind], 
      ; temp[ids_ind], etc to be in the same order as the group catalog id list    
      calcMatch,galcat.trVel_ids,tr_ids,galcat_ind,trids_ind,count=countMatch
      trids_ind = trids_ind[sort(galcat_ind)]
      
      tr_ids     = !NULL
      galcat_ind = !NULL
      
      ; no (rho,temp) cut (would need to load tracer parents which isn't currently done)
      
      ; load pos to calculate tracer positions relative to halo centers
      pos = loadSnapshotSubset(sP=sP,partType='tracerVel',field='pos')
      
      tr_pos = pos[*,trids_ind]
      pos = !NULL
      
      ; calculate current distance of gas particle from smoothed halo center position for galaxy members
      rad_pri  = periodicDists(reform(mt.hPos[mt.maxSnap-m,*,mt.gcIndOrigTrVel]),tr_pos,sP=sP)
      rad_pri /= mt.hVirRad[mt.maxSnap-m,mt.gcIndOrigTrVel]
      
      tr_pos  = !NULL
      
      foreach rVirFac,[sP.rVirFacs,sP.radcut_rvir,1.0],k do begin
        ; for particles who are still within each radius, check if they have passed beyond
        count_rad = 0
        count_out = 0
        
        if k eq n_elements(sP.rVirFacs) then begin
          ; for the last two iterations, take 0.15/1.0rvir and do not use accMask
          ; thereby recording the earliest/highest redshift crossing
          if count_inPar gt 0 then begin
            rad_w = where(rad_pri ge rVirFac and prevRad lt rVirFac and accMask eq 0B,count_rad)
            out_w = where(rad_pri le rVirFac and prevRad gt rVirFac and accMask eq 0B,count_out)          
          endif
        endif else begin
          ; take rVirFac and use accMask (skip if past the end of halo tracking)
          if count_inPar gt 0 then begin
            rad_w = where(rad_pri ge rVirFac and prevRad lt rVirFac and accMask eq 0B and $
                          r.accTime[k,*] eq -1,count_rad)
            out_w = where(rad_pri le rVirFac and prevRad gt rVirFac and accMask eq 0B and $
                          r.outTime[k,*] eq -1,count_out)          
          endif
        endelse
        
        print,' ['+string(m,format='(i3)')+'] ['+str(k)+'] r='+string(rVirFac,format='(f4.2)')+$
          ' accreted now counts '+string(count_gal,format='(i7)')+' ('+$
          string(float(count_gal)/n_elements(gal_pri)*100,format='(f4.1)')+'%) '+$
          string(count_gmem,format='(i7)')+' ('+$
          string(float(count_gmem)/n_elements(gmem_pri)*100,format='(f4.1)')+'%)'
        
        ; interpolate these (time,radii) to find time crossing the virial radius
        times = [prevTime,h.time]
        
        for i=0,count_rad-1 do begin
          radii = [ prevRad[rad_w[i]],rad_pri[rad_w[i]] ]
          time = interpol(times,radii,rVirFac) ; lerp time to r/rvir=rVirFac
          r.accTime[k,rad_w[i]] = time
          
          if k eq 0 then begin
            tvir = [ mt.hVirTemp[mt.maxSnap-m-1,mt.gcIndOrigTrVel[rad_w[i]]], $
                     mt.hVirTemp[mt.maxSnap-m,mt.gcIndOrigTrVel[rad_w[i]]] ]
            if tvir[0] eq 0.0 or abs(tvir[0]-tvir[1]) gt 0.5 then tvir[0]=tvir[1]
            tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
            r.accHaloTvir[rad_w[i]] = tvir
          endif
        endfor
        
        ; if we are on the first snapshot, override accretion times with -1 to indicate always outside rad
        if m eq mt.maxSnap then begin
          r.accTime[k,rad_w] = -1
          r.outTime[k,rad_w] = -1
        endif else begin
          ; otherwise, update counters for the number of particles we have found the accretion times of
          accCount[k] += count_rad
          outCount[k] += count_out
        endelse
        
      endforeach
      
      ; adaptive: update mask, mark all children of halos whose tracking ends at this snapshot
      w = where(trMinSnap eq sP.snap,count)
      if count gt 0 then accMask[w] = 1B
      
      ; store current radius of particles
      prevRad  = rad_pri
      prevTime = h.time
      
      ; free some memory for next load
      rad_w    = !NULL
      rad_pri  = !NULL
    endfor
              
    foreach rVirFac,[sP.rVirFacs,1.0],k do $
      print,'['+str(k)+'] r='+string(sP.rVirFacs[k],format='(f4.2)')+' found accretion times for ['+$
        str(accCount[k])+' of '+str(n_elements(galcat.trVel_ids))+'], outflow times for ['+$
        str(outCount[k])+' of '+str(n_elements(galcat.trVel_ids))+']'
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
                
  endif
  
  sP.snap = origSnap ; restore sP.snap
  return,r
  
end

; shyPlot(): at a target redshift, make a selection of galaxy star tracers which were in gas cells
;            in the previous timestep (newly formed stars). normalize the time elapsed since the 
;            first 1.0rvir crossing of each such tracer by the age of the universe at the time of
;            that crossing. plot the resulting normalized 'accTime' distribution, comparing
;            TRACER to FEEDBACK at z=0,1,2,3

function shyPlotBin, sP=sP, aMode=aMode, allTypes=allTypes

  if keyword_set(allTypes) then allTypesStr = '.allTypes' else allTypesStr = ''

  ; set saveFilename and check for existence
  saveFilename = sP.derivPath + 'shyPlot.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(sP.snap)+'.mode'+str(aMode)+allTypesStr+'.sav'
  
  if file_test(saveFilename) then begin
    restore, saveFilename
    return, r
  endif

  nSnapsBack = 1
  
  ; load at target redshift
  galcat = galaxyCat(sP=sP)
  
  ; make tracer selection (children of stars in galaxies only)
  if keyword_set(allTypes) then begin
    print,'AllTypes = 1'
    
    ; load all star ids
    star_ids = loadSnapshotSubset(sP=sP, partType='star', field='ids')
    
    ; replicate parent id list for all galcat tracers
    tr_parids = galcat.IDs[ replicate_var(galcat.trMC_cc) ]
    
    ; decide overlap, use those tracers (tr_parids not unique, use value_locate approach)
    sort_inds = calcSort(star_ids)
    star_ids_sorted = star_ids[sort_inds]
        
    par_ind = value_locate(star_ids_sorted,tr_parids) ; indices to star_ids_sorted
    par_ind = sort_inds[par_ind>0] ; indices to star_ids (>0 removes -1 entries, which are removed next line)
    tr_parids_ind = where(star_ids[par_ind] eq tr_parids,count_inPar) ; verify we actually matched the ID

    ; save master list of ids/inds to search for
    tr_inds_orig = par_ind[tr_parids_ind]
    tr_ids_orig  = galcat.trMC_ids[ tr_inds_orig ]
  endif else begin
    print,'AllTypes = 0'
    type = ( galcat.type[ replicate_var(galcat.trMC_cc) ] )
  
    tr_inds_orig = where(type eq galcat.types.stars)
    tr_ids_orig = galcat.trMC_ids[ tr_inds_orig ]
  endelse
  
  ; move back requested number of snapshots in a row, load tracer parent ids for this subset
  tr_ids_final  = []
  tr_inds_final = []

  for i=0,nSnapsBack-1 do begin
    sP.snap -= 1
    tr_ids = loadSnapshotSubset(sP=sP, partType='tracerMC', field='tracerids')
  
    idIndexMap = getIDIndexMap(tr_ids, minid=minid)
    tr_ids = !NULL
    tr_inds = idIndexMap[ tr_ids_orig - minid ]
  
    tr_parids = loadSnapshotSubset(sP=sP, partType='tracerMC', field='parentids')
    tr_parids = tr_parids[ tr_inds ]
  
    ; load gas ids and crossmatch
    gas_ids = loadSnapshotSubset(sP=sP, partType='gas', field='ids')
  
    sort_inds = calcSort(gas_ids)
    gas_ids_sorted = gas_ids[sort_inds]
    gas_ind = value_locate(gas_ids_sorted,tr_parids) ; indices to gas_ids_sorted
    gas_ind = sort_inds[gas_ind>0] ; indices to gas_ids (>0 removes -1 entries, which are removed next line)
    w = where(gas_ids[gas_ind] eq tr_parids,count_inPar) ; verify we actually matched the ID

    tr_ids_final  = [ tr_ids_final, tr_ids_orig[w] ]
    tr_inds_final = [ tr_inds_final, tr_inds_orig[w] ] ; index galcat.trMC_ids
    print,count_inPar
  endfor
  
  ; now just insure uniqueness
  if nSnapsBack gt 1 then begin
    tr_ids_final = tr_ids_final[sort(tr_ids_final)]
    tr_inds_final = tr_inds_final[sort(tr_inds_final)]
    
    w = where( tr_ids_final ne shift(tr_ids_final,1), count)
    if count eq 0 then message,'Error, expect some duplicates.'
    tr_ids_final = tr_ids_final[w]
    
    w = where( tr_inds_final ne shift(tr_inds_final,1), count)
    if count eq 0 then message,'Error, expect some duplicates.'
    tr_inds_final = tr_inds_final[w]
  endif
  
  print,' Found: ['+str(n_elements(tr_ids_final))+'] of ['+str(n_elements(tr_ids_orig))+'] tracers.'

  gas_ids         = !NULL
  star_ids        = !NULL
  gas_ind         = !NULL
  gas_ids_sorted  = !NULL
  star_ids_sorted = !NULL
  tr_parids_ind   = !NULL
  par_ind         = !NULL
  tr_parids       = !NULL
  tr_inds         = !NULL
  sort_inds       = !NULL
  w               = !NULL
  
  ; return to original snapshot
  sP.snap += nSnapsBack
  
  ; aMode=1: load accTimes of first 1.0rvir crossing of main progenitor branch
  if aMode eq 1 then begin
    at = accretionTimes(sP=sP)
    at = reform( at.accTime[-1,tr_inds_final] )
  
    ; restrict to accretion mode (all are +_rec) and valid accTimes
    ;am = accretionMode(sP=sP)
    ;am = reform( am[tr_inds_final] )
    ;if accMode eq 'smooth'   then $
    ;  w_at = where(at ne -1 and (am eq 1 or am eq 11),count)
    ;if accMode eq 'clumpy'   then $
    ;  w_at = where(at ne -1 and (am eq 2 or am eq 3 or am eq 12 or am eq 13),count)
    ;if accMode eq 'stripped' then $
    ;  w_at = where(at ne -1 and (am eq 4 or am eq 14),count)
    ;if accMode eq 'all' then $
    ;  w_at = where(at ne -1,count)
    ;am = !NULL
  
    ; filter out -1 values
    at = at[where(at gt 0.0)]
    
    ; convert scalefac to redshift
    at = 1.0/at-1 ; replace at by at[w_at] to use accMode selection
  endif
  
  ; aMode=2: group membership, use mergerTreeSubset() of the subgroups, walk back through 
  ; the snapshots, update Parent (for each tracer) at each step, calculate current Parent (for each
  ; tracer) at each step, or -2 for not in group catalog, save first mismatch time as 
  ; the accretion time entering the subgroup
  if aMode eq 2 then begin
    mt = mergerTreeSubset(sP=sP)
    
    ; at each snapshot holds the tracked main progenitor parents, for each tracer
    curParInds = mt.gcIndOrigTrMC[ tr_inds_final ]
    
    ; the calculated accretion time
    at = fltarr( n_elements(curParInds) ) - 1.0
    
    origSnap   = sP.snap
    targetSnap = sP.groupCatRange[0]
    
    for m=origSnap,targetSnap+1,-1 do begin
      sP.snap = m
          
      ; at each snapshot holds the actual subhalo parent, for each tracer
      actualParInds = lonarr( n_elements(tr_inds_final) )
          
      ; load group catalog
      gcCur = loadGroupCat(sP=sP,/readIDs)
      
      ; calculate actualParInds
        ; load all tracer IDs and re-locate starting IDs
        tr_ids = loadSnapshotSubset(sP=sP, partType='tracerMC', field='tracerids')
        
        idIndexMap = getIDIndexMap(tr_ids, minid=minid)
        tr_ids = !NULL
        tr_inds = idIndexMap[ tr_ids_final - minid ] ; overwrite previous
        
        ; load parent IDs for this subset (don't care about type)
        tr_parids = loadSnapshotSubset(sP=sP, partType='tracerMC', field='parentids', inds=tr_inds)
        
        ; note: tr_parids are NOT UNIQUE, use a value_locate approach (not match)
          sort_inds = calcSort(gcCur.IDs)
          GC_ids_sorted = gcCur.IDs[sort_inds]
        
          par_ind = value_locate(GC_ids_sorted,tr_parids) ; indices to GC_ids_sorted
          par_ind = sort_inds[par_ind>0] ; indices to gc.IDs (>0 removes -1 entries, which are removed next line)
          tr_parids_ind = where(gcCur.IDs[par_ind] eq tr_parids,count_inPar) ; verify we actually matched the ID
          gcCur_ids_ind = par_ind[tr_parids_ind] ; indices of matched parents
          countMatch = n_elements(tr_parids_ind)
          
          ; for those tracers with non-matching parent IDs (not in group catalog), set to -2
          matchMask = bytarr( n_elements(tr_parids) )
          matchMask[tr_parids_ind] = 1B
          
          w_nonMatch = where(matchMask eq 0B,count)
          
          if count gt 0 then begin
            actualParInds[w_nonMatch] = -2
          endif
          
          ; for those tracers with matching parent IDs
          if countMatch gt 0 then begin
            ; do a value locate of the index (of gcCur.IDs) with group.offset
            ; use whole FOF!
            parCandidates = value_locate( gcCur.groupOffset, gcCur_ids_ind )

            ; calculate delta(index,offset) to check if < group.len, to verify actually in group
            parDelta = gcCur_ids_ind - gcCur.groupOffset[ parCandidates ]
            
            w = where(parDelta ge 0 and parDelta lt gcCur.groupLen[ parCandidates ],$
              countCand,comp=wc,ncomp=countOutsideCand)
            
            if countCand gt 0 then begin
              subgroupParInds = gcCur.groupFirstSub[ parCandidates[w] ]
              actualParInds[tr_parids_ind[w]] = subgroupParInds
            endif
            
            ; those outside their candidates (should not happen since considering full fof groups)
            if countOutsideCand gt 0 then message,'Error, reconsider'
          endif
      
      ; for now just mismatching parInds, set accretionTime
      w = where(curParInds ne actualParInds and curParInds ne -1,count)

      if max(at[w]) gt 0.0 then message,'Error: About to set at for tracers with already set at.'
      if count gt 0 and sP.snap eq origSnap then message,'Error: Bad parents at sP.snap'

      at[w] = snapNumToRedshift(sP=sP)
      curParInds[w] = -1 ; no longer search/set on these tracers
           
      frac = float(count)*100/n_elements(tr_ids_final)
      print,'['+str(m)+'] frac now set = '+string(frac,format='(f4.1)')+'% ('+$
        str(count)+' of '+str(n_elements(tr_ids_final))+')'

      ; load mergerTree and move to Parent
      Parent = mergerTree(sP=sP)
        
      ; change to parent IDs for the next snapshot
      w = where(curParInds ne -1,count)
      if count eq 0 then message,'error'
        
      curParInds[w] = Parent[curParInds[w]]
        
      frac = float(count)*100/n_elements(tr_ids_final)
      print,'['+str(m)+'] frac remaining = '+string(frac,format='(f4.1)')+'% ('+$
        str(count)+' of '+str(n_elements(tr_ids_final))+')'
          
    endfor ; snapshot
    
    ; restrict times to those we found (already in redshift)
    at = at[where(at gt 0.0)]
  endif
  
  if aMode ne 1 and aMode ne 2 then message,'Error'
  
  ; age of universe (in Gyr) at each accTime
  r = { age_at  : redshiftToAgeFlat(at) ,$
        age_cur : redshiftToAgeFlat(sP.redshift) }
        
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        
  
  return, r
  
end

; shyPlot():
;  plot (1): delta time between "halo accretion time" and "star forming time" (Gyr)
;  plot (2): delta time as above, normalized by age of universe at "halo accretion time"
; aMode=1: "halo accretion time" = first 1.0 rvir crossing from accretionTimes()
; aMode=2: "halo accretion time" = group membership
; allTypes=1: all tracers in galaxyCat, otherwise just 'star' (e.g. inside 0.15rvir cut) tracers

pro shyPlot

  ; config
  redshifts = [3.0,2.0,1.0,0.0]
  runs      = ['feedback','tracer']
  res       = 512
  aMode     = 1
  allTypes  = 1
  
  ; plot config
  yrangeRaw  = [0.0,0.2] ;[0.0,0.1]
  yrangeNorm = [0.0,0.3]
  nBins     = 20
  colors    = ['red','blue'] ; one per run
  
  start_PS,'shyPlot_aMode'+str(aMode)+'_allTypes'+str(allTypes)+'_'+str(res)+'b.eps', xs=3*3.5, ys=4*3.5
  
    pos = plot_pos(rows=4,cols=2,/gap)
    
    foreach redshift,redshifts,i do begin
    
    ;xrange = [-0.1, redshiftToAgeFlat(redshift)*1.0]
    xrange = [-0.1, 2.0]
    
    ; plot (1)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeNorm,title="",$
      xtitle=textoidl("( t_{SF} - t_{acc} ) / t_{H}(t_{acc})"),$
      ytitle="Fraction",/xs,/ys,pos=pos[i*2],/noerase
	
	; plot histogram for each run/redshift combination
      foreach run,runs,j do begin
        sP_cur = simParams(res=res,run=run,redshift=redshift)
        bV     = shyPlotBin(sP=sP_cur, aMode=aMode, allTypes=allTypes)
        
        val = (bV.age_cur - bV.age_at) / bV.age_at

        if min(val) lt xrange[0] or max(val) gt xrange[1] then begin
          print,'Warning: normalized val min = ',min(val),' max = ',max(val)
        endif
        
        
        binSize = (xrange[1] - xrange[0]) / float(nBins)
	  h = histogram(val, bin=binSize, loc=loc, min=xrange[0], max=xrange[1])
	  h = float(h)/total(h)
	  
	  cgPlot,loc+binSize*0.5,h,color=cgColor(colors[j]),/overplot
	endforeach
	
	; legend
	legend,['z = '+string(redshift,format='(f3.1)'),runs],textcolors=['black',colors],/top,/right
	
    ; plot (2)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeRaw,title="",$
      xtitle=textoidl("( t_{SF} - t_{acc} ) [Gyr]"),$
      ytitle="Fraction",/xs,/ys,pos=pos[i*2+1],/noerase
	
	; plot histogram for each run/redshift combination
      foreach run,runs,j do begin
        sP_cur = simParams(res=res,run=run,redshift=redshift)
        bV     = shyPlotBin(sP=sP_cur, aMode=aMode, allTypes=allTypes)
        
        val = (bV.age_cur - bV.age_at)
        
        if min(val) lt xrange[0] or max(val) gt xrange[1] then begin
          print,'Warning: raw val min = ',min(val),' max = ',max(val)
        endif
        
        binSize = (xrange[1] - xrange[0]) / float(nBins)
	  h = histogram(val, bin=binSize, loc=loc, min=xrange[0], max=xrange[1])
	  h = float(h)/total(h)
	  
	  cgPlot,loc+binSize*0.5,h,color=cgColor(colors[j]),/overplot
	endforeach
	
	; legend
	legend,['z = '+string(redshift,format='(f3.1)'),runs],textcolors=['black',colors],/top,/right
	
    endforeach ;redshifts
      
  end_PS
  
  stop
  
end
