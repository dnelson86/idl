; accretionTimes.pro
; gas accretion project - past radial history of gas elements (virial radius crossing)
; dnelson jul.2013

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
  nVirFacs = n_elements(sP.rVirFacs) + 1 ; one extra for first rvir crossing
  mt = mergerTreeSubset(sP=sP,/verbose)
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
  
  ; OLD START
  restore,'/n/home07/dnelson/data3/sims.gadget/128_20Mpc/data.files.old/galcat.G128.189.sav'
  restore,'/n/home07/dnelson/data3/sims.gadget/128_20Mpc/data.files.old/groupmemcat.G128.189.sav'
  restore,'/n/home07/dnelson/data3/sims.gadget/128_20Mpc/data.files.old/starcat.G128.189.sav'
  galcat2 = { galaxyIDs:galaxyids, groupmemIDs:groupmemids, stellarIDs:stellarIDs ,$
              galaxyOff:galaxyoff, groupmemOff:groupmemoff, stellarOff:stellaroff ,$
              galaxyLen:galaxylen, groupmemLen:groupmemlen, stellarLen:stellarlen }
  restore,'/n/home07/dnelson/data3/sims.gadget/128_20Mpc/data.files.old/mTreeAdaSub.G128.189-51.sK3.sav'
  mt2 = r
  
  accCount2 = { gal   : ulonarr(nVirFacs), gmem    : ulonarr(nVirFacs), stars : ulonarr(nVirFacs), $
               galRT : 0L, starsRT : 0L }  
  ; OLD END
  
  ; NEWGAL START
  mts_type = galcat.type[mt.galcatSub]
  mts_galcatSub_gal  = where(mts_type eq 1)
  mts_galcat_gal     = mt.galcatSub[mts_galcatSub_gal]
  
  new_gal_ids = galcat.ids[mts_galcat_gal]
  old_gal_ids = galcat2.galaxyids[mt2.galcatsub.gal]
  
  calcMatch,new_gal_ids,old_gal_ids,ind_new,ind_old,count=countMatch
  if countMatch ne n_elements(new_gal_ids) then message,'Fail'
  
  mts_galcat_gal = mts_galcat_gal[ind_new[sort(ind_old)]] ; rearrange new mts_X_gal to match old order
  mts_galcatSub_gal = mts_galcatSub_gal[ind_new[sort(ind_old)]]
  
  if ~array_equal(galcat.ids[mts_galcat_gal],old_gal_ids) then message,'Fail' ; gal ids match
  ; NEWGAL END
  
  ; count of how many particles/tracers we tracked through r_vir
  accCount   = ulonarr(nVirFacs)
  outCount   = ulonarr(nVirFacs)
  accCountRT = 0L
  outCountRT = 0L
  
  prevTime = 0 ; scale factor at previous snapshot
  
  reportMemory

  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    print,'Calculating new accretion time using ( SPH Particles ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
      
    if ~file_test(resFilename) then begin ; no restart
      ; replicate hMinSnap for each child gas element
      gasMinSnap = mt.hMinSnap[mt.gcIndOrig]
  
      ; store the r/rvir of each at the previous snapshot for interpolation, and as a monotonic check
      prevRad = fltarr(n_elements(mt.galcatSub))
      
      accMaskIn = bytarr(nVirFacs,n_elements(mt.galcatSub))
      accMaskOut = bytarr(nVirFacs,n_elements(mt.galcatSub))      

      ; store the main arrays as a structure so we can write them directly
      r = {accTime       : fltarr(nVirFacs,n_elements(mt.galcatSub))-1   ,$
           accTimeRT     : fltarr(n_elements(mt.galcatSub))-1            ,$
           outTime       : fltarr(nVirFacs,n_elements(mt.galcatSub))-1   ,$
           outTimeRT     : fltarr(n_elements(mt.galcatSub))-1            ,$
           accHaloTvir   : fltarr(n_elements(mt.galcatSub))              ,$
           rVirFacs      : sP.rVirFacs                                        }
           
      ; OLD START
      ; replicate hMinSnap for each child gas element
      gasMinSnap2 = { gal   : mt2.hMinSnap[mt2.gcIndOrig.gal]   ,$
                      gmem  : mt2.hMinSnap[mt2.gcIndOrig.gmem]  ,$
                      stars : mt2.hMinSnap[mt2.gcIndOrig.stars]  }  
  
      ; store the r/rvir of each at the previous snapshot for interpolation, and as a monotonic check
      prevRad2 = { gal   : fltarr(n_elements(mt2.galcatSub.gal))  ,$
                   gmem  : fltarr(n_elements(mt2.galcatSub.gmem)) ,$
                   stars : fltarr(n_elements(mt2.galcatSub.stars)) }
      
      accMask2 = { gal   : bytarr(nVirFacs,n_elements(mt2.galcatSub.gal))  ,$
                   gmem  : bytarr(nVirFacs,n_elements(mt2.galcatSub.gmem)) ,$
                   stars : bytarr(nVirFacs,n_elements(mt2.galcatSub.stars)) }
      
      ; store the main arrays as a structure so we can write them directly
      r2 = {accTime_gal       : fltarr(nVirFacs,n_elements(mt2.galcatSub.gal))-1   ,$
            accTime_gmem      : fltarr(nVirFacs,n_elements(mt2.galcatSub.gmem))-1  ,$
            accTime_stars     : fltarr(nVirFacs,n_elements(mt2.galcatSub.stars))-1 ,$
            accTimeRT_gal     : fltarr(n_elements(mt2.galcatSub.gal))-1            ,$
            accTimeRT_stars   : fltarr(n_elements(mt2.galcatSub.stars))-1          ,$
            accHaloTvir_gal   : fltarr(n_elements(mt2.galcatSub.gal))              ,$
            accHaloTvir_gmem  : fltarr(n_elements(mt2.galcatSub.gmem))             ,$
            accHaloTvir_stars : fltarr(n_elements(mt2.galcatSub.stars))            ,$
            rVirFacs          : sP.rVirFacs                                        }
      ; OLD END
      
      ; NEWGAL START
      if ~array_equal(gasMinSnap[mts_galcatSub_gal],gasMinSnap2.gal) then message,'Fail'
      ; NEWGAL END

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
        save,prevRad,gasMinSnap,accMaskIn,accMaskOut,r,prevTime,accCount,outCount,$
          accCountRT,outCountRT,m,filename=resFilename
        print,' --- Done! ---'
      endif
      
      ; load gas ids and match to catalog
      h = loadSnapshotHeader(sP=sP)
      
      ; OLD START
      ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs   
      calcMatch,galcat2.galaxyIDs[mt2.galcatSub.gal],ids,galcat_ind_gal,ids_gal_ind,count=countGal
      ids_gal_ind = ids_gal_ind[sort(galcat_ind_gal)]
      
      calcMatch,galcat2.groupmemIDs[mt2.galcatSub.gmem],ids,galcat_ind_gmem,ids_gmem_ind,count=countGmem
      ids_gmem_ind = ids_gmem_ind[sort(galcat_ind_gmem)]
      
      calcMatch,galcat2.stellarIDs[mt2.galcatSub.stars],ids,galcat_stars_ind,ids_stars_ind,count=countStars
      ids_stars_ind    = ids_stars_ind[sort(galcat_stars_ind)]
      galcat_stars_ind = galcat_stars_ind[sort(galcat_stars_ind)] ; use to incompletely fill r.x_stars
      
      if countGal ne n_elements(mt2.galcatSub.gal) or countGmem ne n_elements(mt2.galcatSub.gmem) then $
        message,'Error: Failed to locate all of gal/gmem in gas_ids (overflow64?).'  
        
      u = loadSnapshotSubset(sP=sP,partType='gas',field='u')
      u_gal   = u[ids_gal_ind]
      
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
      nelec_gal   = nelec[ids_gal_ind]
      
      temp_gal = convertUtoTemp(u_gal,nelec_gal)
        
      dens = loadSnapshotSubset(sP=sP,partType='gas',field='density')
      dens_gal = dens[ids_gal_ind]
      dens = !NULL
      
      ; scale Torrey+ (2011) galaxy cut to physical density
      scalefac = snapNumToRedshift(sP=sP,/time) ; time flag gives simulation time = scale factor
      a3inv = 1.0 / (scalefac*scalefac*scalefac)
      dens_gal *= a3inv
      
      ; mark any galaxy gas failing cut as accreted at this snapshot (if not previously marked)
      wGalCut = where(alog10(temp_gal) - sP.galcut_rho * alog10(dens_gal) ge sP.galcut_T and $
                      r2.accTimeRT_gal eq -1,countGalCut)
               
      if countGalCut gt 0 then r2.accTimeRT_gal[wGalCut] = h.time
      accCount2.galRT += countGalCut
      
      print,' ['+string(m,format='(i3)')+'] [-] (rhot)'+$
        ' '+strpad('gas',5)+' accNowCounts '+string(countGalCut,format='(i7)')+' ('+$
        string(float(countGalCut)/n_elements(mt2.galcatSub.gal)*100,format='(f4.1)')+'%)'+$
        ' outNowCounts '+string(0,format='(i7)')+' ('+$
        string(float(0)/n_elements(mt.galcatSub)*100,format='(f4.1)')+'%) || cumIn '+$
        string(accCount2.galRT,format='(i7)')+' ('+$
        string(float(accCount2.galRT)/n_elements(mt2.galcatSub)*100,format='(f4.1)')+'%)'+$
        ' cumOut '+string(0,format='(i7)')+' ('+$
        string(float(0)/n_elements(mt.galcatSub)*100,format='(f4.1)')+'%)'
        
      ; load pos to calculate radii
      pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
      
      pos_gal  = pos[*,ids_gal_ind]
      pos_gmem = pos[*,ids_gmem_ind]
      
      pos = !NULL

      ; calculate current distance of gas particle from smoothed halo center position for galaxy members
      gal_pri  = periodicDists(reform(mt2.hPos[mt2.maxSnap-m,*,mt2.gcIndOrig.gal]),pos_gal,sP=sP)
      gal_pri /= mt2.hVirRad[mt2.maxSnap-m,mt2.gcIndOrig.gal]
      pos_gal = !NULL
      
      ; for group members
      gmem_pri = periodicDists(reform(mt2.hPos[mt2.maxSnap-m,*,mt2.gcIndOrig.gmem]),pos_gmem,sP=sP)
      gmem_pri /= mt2.hVirRad[mt2.maxSnap-m,mt2.gcIndOrig.gmem]
      pos_gmem = !NULL
      
      ; loop over each target radius
      foreach rVirFac,[sP.rVirFacs,1.0],k do begin
        ; for particles who are still within r_vir, check if they have passed beyond


          if k eq n_elements(sP.rVirFacs) then begin
            ; for the last iteration, take 1.0rvir and do not use accMask
            ; thereby recording the earliest/highest redshift crossing
            gal_w   = where(gal_pri ge rVirFac and prevRad2.gal lt rVirFac,count_gal)
            gmem_w  = where(gmem_pri ge rVirFac and prevRad2.gmem lt rVirFac,count_gmem)  
          endif else begin
            ; take rVirFac and use accMask (skip if past the end of halo tracking)
            gal_w   = where(gal_pri ge rVirFac and prevRad2.gal lt rVirFac and accMask2.gal[k,*] eq 0B,count_gal)
            gmem_w  = where(gmem_pri ge rVirFac and prevRad2.gmem lt rVirFac and accMask2.gmem[k,*] eq 0B,count_gmem)  
          endelse
        
        count_stars = 0 & stars_pri = 0
        
        print,' ['+string(m,format='(i3)')+'] ['+str(k)+'] r='+string(rVirFac,format='(f4.2)')+$
          ' accreted now counts '+string(count_gal,format='(i7)')+' ('+$
          string(float(count_gal)/n_elements(mt2.galcatSub.gal)*100,format='(f4.1)')+'%) '+$
          string(count_gmem,format='(i7)')+' ('+$
          string(float(count_gmem)/n_elements(mt2.galcatSub.gmem)*100,format='(f4.1)')+'%)'+$
          string(count_stars,format='(i7)')+' ('+$
          string(float(count_stars)/n_elements(mt2.galcatSub.stars)*100,format='(f4.1)')+'%)'+' || cum '+$
          string(accCount2.gal[k],format='(i7)')+' ('+$
          string(float(accCount2.gal[k])/n_elements(mt2.galcatSub.gal)*100,format='(f4.1)')+'%) '+$
          string(accCount2.gmem[k],format='(i7)')+' ('+$
          string(float(accCount2.gmem[k])/n_elements(mt2.galcatSub.gmem)*100,format='(f4.1)')+'%)'+$
          string(accCount2.stars[k],format='(i7)')+' ('+$
          string(float(accCount2.stars[k])/n_elements(mt2.galcatSub.stars)*100,format='(f4.1)')+'%)'
        
        ; interpolate these (time,radii) to find time crossing the virial radius
        times = [prevTime,h.time]
        
        for i=0,count_gal-1 do begin
          radii = [ prevRad2.gal[gal_w[i]],gal_pri[gal_w[i]] ]
          time = interpol(times,radii,rVirFac) ; lerp time to r/rvir=rVirFac
          r2.accTime_gal[k,gal_w[i]] = time
          ; record the tvir of the halo at the rvir crossing time
          if k eq 0 then begin
            tvir = [ mt2.hVirTemp[mt2.maxSnap-m-1,mt2.gcIndOrig.gal[gal_w[i]]], $
                     mt2.hVirTemp[mt2.maxSnap-m,mt2.gcIndOrig.gal[gal_w[i]]] ]
                     
            ; fix: interpolating to the m-1 snapshot could be untracked, in which case mt.hVirTemp=0
            if tvir[0] eq 0.0 or abs(tvir[0]-tvir[1]) gt 0.5 then tvir[0]=tvir[1]
            tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
            r2.accHaloTvir_gal[gal_w[i]] = tvir
          endif
        endfor
        
        for i=0,count_gmem-1 do begin
          radii = [ prevRad2.gmem[gmem_w[i]],gmem_pri[gmem_w[i]] ]
          time = interpol(times,radii,rVirFac) ; lerp time to r/rvir=rVirFac
          r2.accTime_gmem[k,gmem_w[i]] = time
          if k eq 0 then begin
            tvir = [ mt2.hVirTemp[mt2.maxSnap-m-1,mt2.gcIndOrig.gmem[gmem_w[i]]], $
                     mt2.hVirTemp[mt2.maxSnap-m,mt2.gcIndOrig.gmem[gmem_w[i]]] ]
                     
            if tvir[0] eq 0.0 or abs(tvir[0]-tvir[1]) gt 0.5 then tvir[0]=tvir[1]
            tvir = interpol(tvir,times,time) ; lerp tvir to time=tcross
            r2.accHaloTvir_gmem[gmem_w[i]] = tvir            
          endif
        endfor
        
        ; if we are on the first snapshot, override accretion times with -1 to indicate always outside rad
        if m eq mt.maxSnap then begin
          r2.accTime_gal[k,gal_w] = -1
          r2.accTime_gmem[k,gmem_w] = -1
        endif else begin
          ; otherwise, update counters for the number of particles we have found the accretion times of
          accCount2.gal[k]  += count_gal
          accCount2.gmem[k] += count_gmem
        endelse
        
        ; update mask for particles we no longer search for
        accMask2.gal[k,gal_w]   = 1B
        accMask2.gmem[k,gmem_w] = 1B

      endforeach
      
      ; adaptive: update mask, mark all children of halos whose tracking ends at this snapshot
      w = where(gasMinSnap2.gal eq sP.snap,count)
      if count gt 0 then accMask2.gal[*,w] = 1B
      w = where(gasMinSnap2.gmem eq sP.snap,count)
      if count gt 0 then accMask2.gmem[*,w] = 1B
      
      ; store current radius of particles
      prevRad2.gal  = gal_pri
      prevRad2.gmem = gmem_pri
      
      print,''
      ; OLD END
      
      ; --- for each each possible parent particle type, match child tracers and save times ---
      ;parPartTypes = ['gas','stars']
      parPartTypes = ['gas']
      
      foreach partType,parPartTypes do begin
      
        par_ids = loadSnapshotSubset(sP=sP,partType=partType,field='ids')
            
        ; note: tr_parids are NOT UNIQUE, use a value_locate approach (not match)
        ; they are in fact unique for each partType for this SPH approach, but keep same approach
        sort_inds = calcSort(par_ids)
        par_ids_sorted = par_ids[sort_inds]
        
        ; locate
        par_ind = value_locate(par_ids_sorted,galcat.ids[mt.galcatSub]) ; indices to par_ids_sorted
        par_ind = sort_inds[par_ind>0] ; indices to par_ids (>0 removes -1 entries, which are removed next line)
        w = where(par_ids[par_ind] eq galcat.ids[mt.galcatSub],count_inPar) ; verify we actually matched the ID
        
        if count_inPar gt 0 then begin
          tr_parids_inPar = par_ind[w]
          galcat_ind_inPar = w ; used to be galcat_par_ind[w] but this is unnecessary
          par_ind = !NULL
        endif
        
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
          
          if count_inPar gt 0 then temp = temp[tr_parids_inPar]
      
          ; scale Torrey+ (2012) galaxy cut to physical density
          scalefac = snapNumToRedshift(sP=sP,/time) ; time flag gives simulation time = scale factor
          a3inv = 1.0 / (scalefac*scalefac*scalefac)
      
          dens = loadSnapshotSubset(sP=sP,partType='gas',field='density')
          if count_inPar gt 0 then dens = dens[tr_parids_inPar] * a3inv
      
          ; mark any galaxy gas failing cut as accreted at this snapshot, if not previously
          ; marked due to RT cut or due to accreting across 0.15rvir
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
            ' '+strpad(partType,5)+' accNowCounts '+string(countRTCut,format='(i7)')+' ('+$
            string(float(countRTCut)/n_elements(mt.galcatSub)*100,format='(f4.1)')+'%)'+$
            ' outNowCounts '+string(countRTCutOut,format='(i7)')+' ('+$
            string(float(countRTCutOut)/n_elements(mt.galcatSub)*100,format='(f4.1)')+'%) || cumIn '+$
            string(accCountRT,format='(i7)')+' ('+$
            string(float(accCountRT)/n_elements(mt.galcatSub)*100,format='(f4.1)')+'%)'+$
            ' cumOut '+string(outCountRT,format='(i7)')+' ('+$
            string(float(outCountRT)/n_elements(mt.galcatSub)*100,format='(f4.1)')+'%)'
          
        endif ; partType eq 'gas'
            
        ; load parent positions and convert to tracer positions
        par_pos = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
        if count_inPar gt 0 then child_pos = par_pos[*,tr_parids_inPar]
      
        ;par_pos = !NULL
        ;tr_parids_inPar  = !NULL

        ; calculate current distance of parent from smoothed halo center position for galaxy members
        if count_inPar gt 0 then begin
          rad_pri = periodicDists( reform(mt.hPos[mt.maxSnap-m,*,mt.gcIndOrig[galcat_ind_inPar]]),child_pos,sP=sP )
            
          rad_pri /= mt.hVirRad[mt.maxSnap-m,mt.gcIndOrig[galcat_ind_inPar]]
         ; child_pos  = !NULL
        endif
        
        ; loop over each target radius
        foreach rVirFac,[sP.rVirFacs,1.0],k do begin
          count_rad = 0
          count_out = 0
        
          if k eq n_elements(sP.rVirFacs) then begin
            ; for the last iteration, take 1.0rvir and do not use accMask
            ; thereby recording the earliest/highest redshift crossing
            if count_inPar gt 0 then begin
              rad_w = where(rad_pri ge rVirFac and prevRad[galcat_ind_inPar] lt rVirFac,count_rad)
              out_w = where(rad_pri le rVirFac and prevRad[galcat_ind_inPar] gt rVirFac,count_out)          
            endif
          endif else begin
            ; take rVirFac and use accMask (skip if past the end of halo tracking)
            if count_inPar gt 0 then begin
              rad_w = where(rad_pri ge rVirFac and prevRad[galcat_ind_inPar] lt rVirFac and $
                            accMaskIn[k,galcat_ind_inPar] eq 0B,count_rad)
              out_w = where(rad_pri le rVirFac and prevRad[galcat_ind_inPar] gt rVirFac and $
                            accMaskOut[k,galcat_ind_inPar] eq 0B,count_out)          
            endif
          endelse
        
          print,' ['+string(m,format='(i3)')+'] ['+str(k)+'] r='+string(rVirFac,format='(f4.2)')+$
            ' '+strpad(partType,5)+' accNowCounts '+string(count_rad,format='(i7)')+' ('+$
            string(float(count_rad)/n_elements(mt.galcatSub)*100,format='(f4.1)')+'%)'+$
            ' outNowCounts '+string(count_out,format='(i7)')+' ('+$
            string(float(count_out)/n_elements(mt.galcatSub)*100,format='(f4.1)')+'%) || cumIn '+$
            string(accCount[k],format='(i7)')+' ('+$
            string(float(accCount[k])/n_elements(mt.galcatSub)*100,format='(f4.1)')+'%)'+$
            ' cumOut '+string(outCount[k],format='(i7)')+' ('+$
            string(float(outCount[k])/n_elements(mt.galcatSub)*100,format='(f4.1)')+'%)'
        
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
        
          ; update mask for particles we no longer search for
          if count_inPar gt 0 then begin
            if count_rad gt 0 then accMaskIn[k,galcat_ind_inPar[rad_w]] = 1B
            if count_out gt 0 then accMaskOut[k,galcat_ind_inPar[out_w]] = 1B
          endif
                
      endforeach ; rVirFac
      
      ; adaptive: update mask, mark all children of halos whose tracking ends at this snapshot
      w = where(gasMinSnap eq sP.snap,count)
      
      if count gt 0 then begin
        accMaskIn[*,w] = 1B
        accMaskOut[*,w] = 1B
      endif
      
      prevTime = h.time
      
      ; store current radius of particles
      if count_inPar gt 0 then prevRad[galcat_ind_inPar] = rad_pri
      
      ; free some memory for next load
      rad_w   = !NULL
      out_w   = !NULL
      rad_pri = !NULL
      endforeach ; partType
      
      ; NEWGAL START
      
      ; check prevRad
      prevRad_new_gal = prevRad[mts_galcatSub_gal]
      prevRad_old_gal = prevRad2.gal
      
      w = where(prevRad_new_gal ne prevRad_old_gal,count)
      if count gt 0 then begin
        foreach ind,w do if finite(prevRad_new_gal[ind]) or finite(prevRad_old_gal[ind]) then $
          message,'Fail prevRad'
      endif
      
      ; check accMask
      accMask_new_gal = accMaskIn[*,mts_galcatSub_gal]
      accMask_old_gal = accMask2.gal
      
      if ~array_equal(accMask_new_gal,accMask_old_gal) then begin
        ; which dimension fails?
        for i=0,n_elements(accMask_new_gal[*,0])-1 do begin
          if ~array_equal(accMask_new_gal[i,*],accMask_old_gal[i,*]) then print,'Fail ',i
        endfor
        message,'Fail accMask'
      endif
      
      ; check accTime (rad)
      accTime_new_gal = r.accTime[*,mts_galcatSub_gal]
      accTime_old_gal = r2.accTime_gal
      
      if ~array_equal(accTime_new_gal,accTime_old_gal) then message,'Fail accTime'
      
      ; check accTimeRT (dens/temp)
      accTimeRT_new_gal = r.accTimeRT[mts_galcatSub_gal]
      accTimeRT_old_gal = r2.accTimeRT_gal
      
      if ~array_equal(accTimeRT_new_gal,accTimeRT_old_gal) then message,'Fail accTimeRT'
      
      ; check accHaloTvir
      accHaloTvir_new_gal = r.accHaloTvir[mts_galcatSub_gal]
      accHaloTvir_old_gal = r2.accHaloTvir_gal
      
      if ~array_equal(accHaloTvir_new_gal,accHaloTvir_old_gal) then message,'Fail accHaloTvir'
      
      stop
      ; NEWGAL END
      
    endfor ;m
    
    ; print final results
    print,'[-] (rhot) found accretion/outflow times for ['+$
      str(accCountRT)+' / '+str(outCountRT)+' of '+str(n_elements(mt.galcatSub))+']'
    foreach rVirFac,[sP.rVirFacs,1.0],k do $
      print,'['+str(k)+'] r='+string(rVirFac,format='(f4.2)')+' found accretion times for ['+$
        str(accCount[k])+' or '+string(float(accCount[k])/n_elements(mt.galcatSub)*100,format='(f4.1)')+'%], '+$
        ' outflow times for ['+$
        str(outCount[k])+' or '+string(float(outCount[k])/n_elements(mt.galcatSub)*100,format='(f4.1)')+'%]'
    
    ; save
    ;save,r,filename=saveFilename
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
      
      ; locate tracer children at starting snapshot
      gcIndOrigTr = mergerTreeRepParentIDs(mt=mt,galcat=galcat,sP=sP,/compactMtS,$
        trids=galcat_trids,gc_cc=galcat_cc)
          
      galcat = !NULL ; not used past this point    
          
      ; replicate hMinSnap for each child gas/star tracer
      trMinSnap = mt.hMinSnap[gcIndOrigTr]
          
      ; store the r/rvir of each at the previous snapshot for interpolation, and as a monotonic check
      prevRad = fltarr(n_elements(galcat_trids))
      
      accMaskIn  = bytarr(nVirFacs,n_elements(galcat_trids))
      accMaskOut = bytarr(nVirFacs,n_elements(galcat_trids))
      
      ; store the main arrays as a structure so we can write them directly for all tracers
      r = {accTime        : fltarr(nVirFacs,n_elements(galcat_trids))-1   ,$
           accTimeRT      : fltarr(n_elements(galcat_trids))-1            ,$
           outTime        : fltarr(nVirFacs,n_elements(galcat_trids))-1   ,$
           outTimeRT      : fltarr(n_elements(galcat_trids))-1            ,$
           accHaloTvir    : fltarr(n_elements(galcat_trids))              ,$
           child_counts   : galcat_cc                                     ,$
           rVirFacs       : sP.rVirFacs                                    }
           
      galcat_cc   = !NULL
    endif else begin
      ; restart
      if ~keyword_set(restart) then message,'Error: Restart file exists but restart not requested.'
      restore,resFilename,/verbose
      snapRange[0] = m
    endelse

    for m=snapRange[0],snapRange[1],-1 do begin
      sP.snap = m
      
      reportMemory
      
      ; save restart?
      if m mod 10 eq 0 and m lt snapRange[0] and keyword_set(restart) then begin
        print,' --- Writing restart! ---'
        save,prevRad,accMaskIn,accMaskOut,trMinSnap,accCount,accCountRT,r,gcIndOrigTr,$
             galcat_trids,prevTime,m,filename=resFilename
        print,' --- Done! ---'
      endif
      
      ; load tracer ids and match to child ids from zMin
      h = loadSnapshotHeader(sP=sP)
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
      
      idIndexMap = getIDIndexMap(tr_ids,minid=minid)
      trids_ind  = idIndexMap[galcat_trids-minid]
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
        w = where(par_ids[par_ind] eq tr_parids,count_inPar) ; verify we actually matched the ID
        
        if count_inPar gt 0 then begin
          tr_parids_inPar = par_ind[w]
          galcat_ind_inPar = w ; used to be galcat_par_ind[w] but this is unnecessary
          par_ind = !NULL
        endif
        
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
          
          if count_inPar gt 0 then temp = temp[tr_parids_inPar]
      
          ; scale Torrey+ (2012) galaxy cut to physical density
          scalefac = snapNumToRedshift(sP=sP,/time) ; time flag gives simulation time = scale factor
          a3inv = 1.0 / (scalefac*scalefac*scalefac)
      
          dens = loadSnapshotSubset(sP=sP,partType='gas',field='density')
          if count_inPar gt 0 then dens = dens[tr_parids_inPar] * a3inv
      
          ; mark any galaxy gas failing cut as accreted at this snapshot, if not previously
          ; marked due to RT cut or due to accreting across 0.15rvir
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
            ' '+strpad(partType,5)+' accNowCounts '+string(countRTCut,format='(i7)')+' ('+$
            string(float(countRTCut)/n_elements(galcat_trids)*100,format='(f4.1)')+'%)'+$
            ' outNowCounts '+string(countRTCutOut,format='(i7)')+' ('+$
            string(float(countRTCutOut)/n_elements(galcat_trids)*100,format='(f4.1)')+'%) || cumIn '+$
            string(accCountRT,format='(i7)')+' ('+$
            string(float(accCountRT)/n_elements(galcat_trids)*100,format='(f4.1)')+'%)'+$
            ' cumOut '+string(outCountRT,format='(i7)')+' ('+$
            string(float(outCountRT)/n_elements(galcat_trids)*100,format='(f4.1)')+'%)'
          
        endif ; partType eq 'gas'
      
        ; load parent positions and convert to tracer positions
        par_pos = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
        if count_inPar gt 0 then tr_pos = par_pos[*,tr_parids_inPar]
      
        par_pos = !NULL
        tr_parids_inPar  = !NULL

        ; calculate current distance of parent from smoothed halo center position for galaxy members
        if count_inPar gt 0 then begin
          rad_pri  = periodicDists( $
            reform(mt.hPos[mt.maxSnap-m,*,gcIndOrigTr[galcat_ind_inPar]]),tr_pos,sP=sP)
            
          rad_pri /= mt.hVirRad[mt.maxSnap-m,gcIndOrigTr[galcat_ind_inPar]]
          tr_pos  = !NULL
        endif
      
      ; loop over each critical radius
      foreach rVirFac,[sP.rVirFacs,1.0],k do begin
        ; for particles who are still within r_vir, check if they have passed beyond
        count_rad = 0
        count_out = 0
        
        if k eq n_elements(sP.rVirFacs) then begin
          ; for the last iteration, take 1.0rvir and do not use accMask
          ; thereby recording the earliest/highest redshift crossing
          if count_inPar gt 0 then begin
            rad_w = where(rad_pri ge rVirFac and prevRad[galcat_ind_inPar] lt rVirFac,count_rad)
            out_w = where(rad_pri le rVirFac and prevRad[galcat_ind_inPar] gt rVirFac,count_out)          
          endif
        endif else begin
          ; take rVirFac and use accMask (skip if past the end of halo tracking)
          if count_inPar gt 0 then begin
            rad_w = where(rad_pri ge rVirFac and prevRad[galcat_ind_inPar] lt rVirFac and $
                          accMaskIn[k,galcat_ind_inPar] eq 0B,count_rad)
            out_w = where(rad_pri le rVirFac and prevRad[galcat_ind_inPar] gt rVirFac and $
                          accMaskOut[k,galcat_ind_inPar] eq 0B,count_out)          
          endif
        endelse
                       
          print,' ['+string(m,format='(i3)')+'] ['+str(k)+'] r='+string(rVirFac,format='(f4.2)')+$
            ' '+strpad(partType,5)+' accNowCounts '+string(count_rad,format='(i7)')+' ('+$
            string(float(count_rad)/n_elements(galcat_trids)*100,format='(f4.1)')+'%)'+$
            ' outNowCounts '+string(count_out,format='(i7)')+' ('+$
            string(float(count_out)/n_elements(galcat_trids)*100,format='(f4.1)')+'%) || cumIn '+$
            string(accCount[k],format='(i7)')+' ('+$
            string(float(accCount[k])/n_elements(galcat_trids)*100,format='(f4.1)')+'%)'+$
            ' cumOut '+string(outCount[k],format='(i7)')+' ('+$
            string(float(outCount[k])/n_elements(galcat_trids)*100,format='(f4.1)')+'%)'
        
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
                tvir = [ mt.hVirTemp[mt.maxSnap-m-1,gcIndOrigTr[curInd]], $
                         mt.hVirTemp[mt.maxSnap-m,gcIndOrigTr[curInd]] ]
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
        
          ; update mask for particles we no longer search for
          if count_inPar gt 0 then begin
            if count_rad gt 0 then accMaskIn[k,galcat_ind_inPar[rad_w]] = 1B
            if count_out gt 0 then accMaskOut[k,galcat_ind_inPar[out_w]] = 1B
          endif
        
        endforeach
      
        ; adaptive: update mask, mark all children of halos whose tracking ends at this snapshot
        w = where(trMinSnap eq sP.snap,count)
        if count gt 0 then begin
          accMaskIn[*,w] = 1B
          accMaskOut[*,w] = 1B
        endif
      
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
      str(accCountRT)+' / '+str(outCountRT)+' of '+str(n_elements(galcat_trids))+']'
    foreach rVirFac,[sP.rVirFacs,1.0],k do $
      print,'['+str(k)+'] r='+string(rVirFac,format='(f4.2)')+' found accretion times for ['+$
        str(accCount[k])+' or '+string(float(accCount[k])/n_elements(galcat_trids)*100,format='(f4.1)')+'%], '+$
        ' outflow times for ['+$
        str(outCount[k])+' or '+string(float(outCount[k])/n_elements(galcat_trids)*100,format='(f4.1)')+'%]'
    
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
  
      ; locate tracer children at starting snapshot
      gcIndOrigTr = mergerTreeRepParentIDs(mt=mt,galcat=galcat,sP=sP,/compactMtS,trids=galcat_trids)
          
      galcat = !NULL ; not used past this point
      
      ; replicate hMinSnap for each child gas element
      trMinSnap = mt.hMinSnap[gcIndOrigTr]
      
      ; store the r/rvir of each at the previous snapshot for interpolation, and as a monotonic check
      prevRad = fltarr(n_elements(galcat_trids))
      
      accMaskIn = bytarr(nVirFacs,n_elements(galcat_trids))
      accMaskOut = bytarr(nVirFacs,n_elements(galcat_trids))
  
      ; store the main arrays as a structure so we can write them directly for all tracers
      r = {accTime       : fltarr(nVirFacs,n_elements(galcat_trids))-1  ,$
           accTimeRT     : fltarr(n_elements(galcat_trids))-1           ,$
           outTime       : fltarr(nVirFacs,n_elements(galcat_trids))-1  ,$
           outTimeRT     : fltarr(n_elements(galcat_trids))-1           ,$
           accHaloTvir   : fltarr(n_elements(galcat_trids))             ,$
           child_counts  : galcat_cc                                    ,$
           rVirFacs      : sP.rVirFacs                               }
 
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
        save,prevRad,trMinSnap,accMaskIn,accMaskOut,accCount,accCountRT,outCount,outCountRT,r,$
             gcIndOrigTr,galcat_trids,prevTime,m,filename=resFilename
        print,' --- Done! ---'
      endif
      
      ; load tracer ids and match to child ids from zMin
      h = loadSnapshotHeader(sP=sP)
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerVel',field='ids')
      
      ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs, need this if we want ids[ids_ind], 
      ; temp[ids_ind], etc to be in the same order as the group catalog id list    
      calcMatch,galcat_trids,tr_ids,galcat_ind,trids_ind,count=countMatch
      trids_ind = trids_ind[sort(galcat_ind)]
      
      tr_ids     = !NULL
      galcat_ind = !NULL
      
      ; no (rho,temp) cut (would need to load tracer parents which isn't currently done)
      
      ; load pos to calculate tracer positions relative to halo centers
      pos = loadSnapshotSubset(sP=sP,partType='tracerVel',field='pos')
      
      tr_pos = pos[*,trids_ind]
      pos = !NULL
      
      ; calculate current distance of gas particle from smoothed halo center position for galaxy members
      rad_pri  = periodicDists(reform(mt.hPos[mt.maxSnap-m,*,gcIndOrigTr]),tr_pos,sP=sP)
      rad_pri /= mt.hVirRad[mt.maxSnap-m,gcIndOrigTr]
      
      tr_pos  = !NULL
      
      foreach rVirFac,[sP.rVirFacs,1.0],k do begin
        ; for particles who are still within each radius, check if they have passed beyond
        count_rad = 0
        count_out = 0
        
        if k eq n_elements(sP.rVirFacs) then begin
          ; for the last iteration, take 1.0rvir and do not use accMask
          ; thereby recording the earliest/highest redshift crossing
          if count_inPar gt 0 then begin
            rad_w = where(rad_pri ge rVirFac and prevRad lt rVirFac,count_rad)
            out_w = where(rad_pri le rVirFac and prevRad gt rVirFac,count_out)          
          endif
        endif else begin
          ; take rVirFac and use accMask (skip if past the end of halo tracking)
          if count_inPar gt 0 then begin
            rad_w = where(rad_pri ge rVirFac and prevRad lt rVirFac and accMaskIn[k,*] eq 0B,count_rad)
            out_w = where(rad_pri le rVirFac and prevRad gt rVirFac and accMaskOut[k,*] eq 0B,count_out)          
          endif
        endelse
        
        print,' ['+string(m,format='(i3)')+'] ['+str(k)+'] r='+string(rVirFac,format='(f4.2)')+$
          ' accreted now counts '+string(count_gal,format='(i5)')+' ('+$
          string(float(count_gal)/n_elements(gal_pri)*100,format='(f4.1)')+'%) '+$
          string(count_gmem,format='(i5)')+' ('+$
          string(float(count_gmem)/n_elements(gmem_pri)*100,format='(f4.1)')+'%)'
        
        ; interpolate these (time,radii) to find time crossing the virial radius
        times = [prevTime,h.time]
        
        for i=0,count_rad-1 do begin
          radii = [ prevRad[rad_w[i]],rad_pri[rad_w[i]] ]
          time = interpol(times,radii,rVirFac) ; lerp time to r/rvir=rVirFac
          r.accTime[k,rad_w[i]] = time
          
          if k eq 0 then begin
            tvir = [ mt.hVirTemp[mt.maxSnap-m-1,gcIndOrigTr[rad_w[i]]], $
                     mt.hVirTemp[mt.maxSnap-m,gcIndOrigTr[rad_w[i]]] ]
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
      if count gt 0 then begin
        accMask[*,w] = 1B
        outMask[*,w] = 1B
      endif
      
      ; store current radius of particles
      prevRad  = rad_pri
      prevTime = h.time
      
      ; free some memory for next load
      rad_w    = !NULL
      rad_pri  = !NULL
    endfor
              
    foreach rVirFac,[sP.rVirFacs,1.0],k do $
      print,'['+str(k)+'] r='+string(sP.rVirFacs[k],format='(f4.2)')+' found accretion times for ['+$
        str(accCount[k])+' of '+str(n_elements(galcat_trids))+'], outflow times for ['+$
        str(outCount[k])+' of '+str(n_elements(galcat_trids))+']'
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
                
  endif
  
  sP.snap = origSnap ; restore sP.snap
  return,r
  
end
