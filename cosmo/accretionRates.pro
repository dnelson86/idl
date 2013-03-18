; accretionRates.pro
; halo cooling project - future accretion history of hot halo gas onto the central galaxy
; dnelson march.2013

; -----------------------------------------------------------------------------------------------------
; accretionRates(): description
; -----------------------------------------------------------------------------------------------------

function accretionRates, sP=sP, restart=restart

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  dynFac = 2.0 ; how many dynamical times to track forward for

  ; load galaxy/group member catalogs at zMin for gas ids to search for
  gc = loadGroupCat(sP=sP,/skipIDs)
  galcat = galaxyCat(sP=sP)
  
  ; target halo list
  galcat_ind = galcatINDList(galcat=galcat) ;all (pri+sec)
  galcat_gmem_ids = galcat.groupmemIDs[galcat_ind.gmem]
  galcat_ind = !NULL
  
  ; load dynamical times of this halo gas
  ts = loadFitTimescales(sP=sP, gcIDList=lindgen(n_elements(galcat.galaxyLen))) ; gcIDList=all
  maxDynTime = float(max(ts.dynTime,/nan))
  galcat_gmem_tdyn = ts.dynTime
  ts = !NULL
  
  ; determine snapshot range
  curAge     = redshiftToAgeFlat(sP.redshift)
  targetAge  = curAge + dynFac * maxDynTime
  
  redshifts = snapNumToRedshift(sP=sP,/all)
  ages      = redshiftToAgeFlat(redshifts)

  minSnap = min(where(ages lt targetAge eq 0,count))
  if count eq 0 or redshifts[minSnap] lt 0.0 then message,'Error: Bad bracket.'
  
  snapRange = [sP.snap+1,minSnap]
  
  ; set saveFilename and check for existence
  saveFilename = sP.derivPath + 'accRates.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(snapRange[0])+'-'+str(snapRange[1])+'.sav'  
  
  if file_test(saveFilename) then begin
    restore, saveFilename
    return, r
  endif

  resFilename = sP.derivPath + 'accRates.'+sP.saveTag+'.restart.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(snapRange[0])+'-'+str(snapRange[1])+'.sav'
  
  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    message,'todo'
  endif
  
  ; MONTE CARLO TRACERS CASE - for all tracers in gmem parents, track forward until inside galaxy
  ; NOTE: both gas and star parents are considered for the galaxy, unlike with TRVEL or SPH
  ; ------------------------
  if sP.trMCPerCell gt 0 then begin

    if ~file_test(resFilename) then begin
      ; no restart

      print,'Calculating new accretion rates using ( TracerMC ) res = '+str(sP.res)+$
            ' in range ['+str(snapRange[0])+'-'+str(snapRange[1])+' down to z = '+$
            string(redshifts[minSnap],format='(f4.2)')+'].'
        
      ; load all tracer children of gas_ids_gmem parents
      galcat_gmem_trids = cosmoTracerChildren(sP=sP, /getInds, gasIDs=galcat_gmem_ids, child_counts=galcat_gmem_cc)
    
      ; convert tracer children indices to tracer IDs
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
      galcat_gmem_trids = tr_ids[galcat_gmem_trids]
      tr_ids = !NULL
    
      ; replicate tdyn at target redshift for each tracer
      galcat_gmem_tr_tdyn = galcat_gmem_tdyn[replicate_var(galcat_gmem_cc)]
      galcat_gmem_tdyn = !NULL
    
      ; list of indices into galcat_gmem_trids that we still search for (start as complete index list)
      galcat_gmem_tr_src = lindgen(n_elements(galcat_gmem_trids))
      galcat_gmem_tr_mask = intarr(n_elements(galcat_gmem_trids)) ; debugging only
    
      ; store the main arrays as a structure so we can write them directly
      r = {accSnap_gmem      : intarr(n_elements(galcat_gmem_trids))-1  ,$
           galcat_gmem_cc    : galcat_gmem_cc                           ,$
           dynFac            : dynFac                                    }    
      
    endif else begin
      ; restart
      if ~keyword_set(restart) then message,'Error: Restart file exists but restart not requested.'
      restore,resFilename,/verbose
      snapRange[0] = m
    endelse
   
    for m=snapRange[0],snapRange[1],1 do begin
      sP.snap = m

      ; save restart?
      if m mod 5 eq 0 and m gt snapRange[0] and keyword_set(restart) then begin
        print,' --- Writing restart! ---'
        save,r,m,galcat_gmem_trids,galcat_gmem_tr_tdyn,galcat_gmem_tr_src,galcat_gmem_tr_mask,$
             filename=resFilename
        print,' --- Done! ---'
      endif
      
      ; load tracer ids and match to child ids from zMin
      h = loadSnapshotHeader(sP=sP)
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
        
      match,galcat_gmem_trids[galcat_gmem_tr_src],tr_ids,galcat_ind,trids_gmem_ind,count=countGmem,/sort
      trids_gmem_ind  = trids_gmem_ind[sort(galcat_ind)] ; rearrange trids_gmem_ind to be ordered as galcat_gmem_trids
      
      tr_ids     = !NULL
      galcat_ind = !NULL
      
      if countGmem ne n_elements(galcat_gmem_tr_src) then message,'Error: Tracer id match counts'
      
      ; load tracer parents IDs, take those only of the tracers we are searching for
      tr_parids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentid')
      tr_parids_gmem  = tr_parids[trids_gmem_ind]
      
      tr_parids       = !NULL
      trids_gmem_ind  = !NULL
      
      ; load galaxy catalog at this redshift
      loc_galcat = galaxyCat(sP=sP, /galaxyOnly)
      
      ; --- for each each possible parent particle type, match child tracers and save times ---
      ; note: we never actually find any in stars in sims.tracers, I assume since the snapshot spacing is fine enough
      parPartTypes = ['gas','stars']
      if sP.gfmWinds then parPartTypes = ['gas','stars','BHs']
      
      accFound = lonarr(n_elements(parPartTypes))
      
      foreach partType,parPartTypes,k do begin
        ; load parent IDs and convert tracer parent IDs -> indices (for those tracers with this partType parent now)
        par_ids = loadSnapshotSubset(sP=sP,partType=partType,field='ids')
        if n_elements(par_ids) eq 0 then continue ; no particles of this type in snapshot
        
        ; note: tr_parids_gal,gmem,stars are NOT UNIQUE, use a value_locate approach (not match)
        sort_inds = sort(par_ids)
        
        ; gmem
        gmem_ind = value_locate(par_ids[sort_inds],tr_parids_gmem)
        gmem_ind = sort_inds[gmem_ind>0]
        w = where(par_ids[gmem_ind] eq tr_parids_gmem,countGmem_inPar)
        
        if countGmem_inPar eq 0 then continue ; no tracers have parents of this type
        
        tr_parids_gmem_inPar = gmem_ind[w] ; particle indices for this snapshot
        galcat_gmem_ind_inPar = w ; indices into galcat_gmem_tr_src array
        
        ; note: now we have par_ids[tr_parids_gmem_inPar] = tr_parids_gmem[galcat_gmem_ind_inPar]
        
        ; determine which of these matched parents is in a galaxy at this snapshot
        gmem_parids = par_ids[tr_parids_gmem_inPar]
        
        ; note: gmem_parids are NOT UNIQUE, use value_locate approach (not match)
        sort_inds = sort(loc_galcat.galaxyIDs)
        
        gmem_ind = value_locate(loc_galcat.galaxyIDs[sort_inds],gmem_parids)
        gmem_ind = sort_inds[gmem_ind>0]
        
        w = where(loc_galcat.galaxyIDs[gmem_ind] eq gmem_parids,countMatch)
        
        ; note: now we have loc_galcat.galaxyIDs[gmem_ind[w]] = gmem_parids[w]
        if countMatch eq 0 then continue ; no parents are in galaxies
                
        ; mark accretion time for those located tracers
        insert_inds = galcat_gmem_tr_src[galcat_gmem_ind_inPar[w]]
        
        r.accSnap_gmem[insert_inds] = sP.snap
        galcat_gmem_tr_mask[insert_inds] += 1
        
        accFound[k] = countMatch
      endforeach ; parPartTypes
      
      ; update search list, restrict each tracer to dynFac*tdyn maximum search time
      galcat_gmem_tr_src = where( (r.accSnap_gmem eq -1) and $
                                  (galcat_gmem_tr_tdyn ge (ages[m]-curAge)/dynFac) ,$
                                  count)
      
      ; sanity checks
      countMask    = total(galcat_gmem_tr_mask) ; number found
      countSrcList = n_elements(galcat_gmem_tr_src) ; number remaining (not yet found)
      countOrig    = n_elements(r.accSnap_gmem)
    
      print,'['+string(m,format='(i3)')+'] Found: '+$
        string(float(countMask*100)/n_elements(r.accSnap_gmem),format='(f4.1)')+$
        '% Remain: '+string(countSrcList,format='(i9)')+' ThisSnapFound: ['+string(accFound[0],format='(i8)')+' gas, '+$
        string(accFound[1],format='(i8)')+' stars]' ; '+string(accFound[2])+' BHs]'
    
      if countOrig-countMask lt countSrcList then message,'Error'
      if max(galcat_gmem_tr_mask) gt 1 then message,'Error: Mask exceeds one.'
      if max(galcat_gmem_tr_mask[galcat_gmem_tr_src]) gt 0 then message,'Error: Active search list has nonzero mask.'
    
    endfor ;m
    
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))    
                                     
  endif
  
  ; VELOCITY TRACERS case
  ; ---------------------
  if sP.trMCPerCell eq -1 then begin
    message,'todo'
  endif
  
end

; binAccretionRates(): turn the "snapshot when accreted onto galaxy" data into actual accretion rates on 
;                      a halo by halo basis, and do a median curve

function binAccretionRates, sP=sP, sgSelect=sgSelect
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  nWindows = 10 ; go back i=1,2,...,nWindows and calculate accretion rates using each
  
  ; check if save exists
  saveFilename = sP.derivPath + 'binnedVals/binAR.' + sP.saveTag + '.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.nw' + str(nWindows) + '.' + sgSelect + '.sav'
  
  ; results exist, return
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif  
  
  ; load
  ar     = accretionRates(sP=sP)
  gc     = loadGroupCat(sP=sP,/skipIDs)
  galcat = galaxyCat(sP=sP)
  
  ; target halo list (and make mask)
  gcIDList = gcIDList(gc=gc,select=sgSelect)
  
  gcIDMask = bytarr(n_elements(galcat.galaxyLen))
  gcIDMask[gcIDList] = 1B
  
  ; halo count and halo masses for median binning
  nHalos   = n_elements(gcIDList)
  gcMasses = codeMassToLogMsun(gc.subgroupMass[gcIDList])
  
  ; time in Gyr back to subsequent snapshots
  curAge     = redshiftToAgeFlat(sP.redshift)
  redshifts  = snapNumToRedshift(sP=sP,/all)
  ages       = redshiftToAgeFlat(redshifts)
  windowsGyr = ages[sP.snap + indgen(nWindows) + 1] - curAge ; Gyr
  
  if sP.trMCPerCell le 0 then massPerPart = sP.targetGasMass ; SPH or vel tracer
  if sP.trMCPerCell gt 0 then massPerPart = sP.trMassConst ; MC tracer  

  ; bin fractions into halo mass bins and make median lines
  logMassBins   = [8.5,linspace(9.0,10.0,21),linspace(10.1,11.0,10),11.1,11.25,11.5,11.75,11.9,13.1]
  logMassNBins  = n_elements(logMassBins)-1
  logMassBinCen = 0.5 * (logMassBins + shift(logMassBins,-1))
  logMassBinCen = logMassBinCen[0:-2]
  
  ; save arrays (per halo and median)
  r = { gmem_accRate_halo    : fltarr(nWindows,n_elements(gcIDList)) + !values.f_nan  ,$
        gmem_accRate_median  : fltarr(nWindows,logMassNBins) + !values.f_nan          ,$
        gcMasses             : gcMasses                                               ,$
        logMassBins          : logMassBins                                            ,$
        logMassNBins         : logMassNBins                                           ,$
        logMassBinCen        : logMassBinCen                                          ,$
        nWindows             : nWindows                                               ,$
        windowsGyr           : windowsGyr                                              }
  
  offset_tr   = 0L ; we don't have an offset table for the tracer inds, so walk all gcIDLIst in order
  offset_gmem = 0L ; as is done in galcatINDList(), same to get primary indices and skip secondary gmem
  
  ; loop over all halos
  for i=0UL,n_elements(galcat.galaxyLen)-1 do begin
    if i mod fix(nHalos/20) eq 0 then print,'i = '+string(i,format='(i6)')+' '+str(float(i)*100/nHalos)+'%'
    
    ; indices for this halo
    gcInd = gcIDList[i]
    
    if galcat.groupmemLen[gcInd] eq 0 then continue ; no gmem gas in this halo
       
    inds = lindgen(galcat.groupmemLen[gcInd]) + offset_gmem
    
    offset_gmem += galcat.groupmemLen[gcInd]
    
    ; replicate child_counts for tracer indices
    tot_children_tr = total(ar.galcat_gmem_cc[inds],/int)
    
    offset_tr += tot_children_tr
    
    if tot_children_tr eq 0 then continue ; although have parents, they have no children
    if gcIDMask[gcID] eq 0B then continue ; don't want to process this halo
    
    inds_tr = lindgen(tot_children_tr) + offset_tr
    
    ; get child tracers in this halo
    loc_accSnapOffset = ar.accSnap_gmem[inds_tr]
    w = where(loc_accSnapOffset ne -1,count)
    if count eq 0 then continue ; no accretion times found for any children
    
    loc_accSnapOffset = loc_accSnapOffset[w] - sP.snap
    
    if min(loc_accSnapOffset) lt 1 then message,'Strange'
    
    ; count for each number of snaps back
    for j=0,nWindows-1 do begin
      w = where(loc_AccSnapOffset le j+1,count)
      if count gt 0 then begin
        r.gmem_accRate_halo[j,i] = count / windowsGyr[j]
      endif
    endfor
    
  endfor
  
  if offset_gmem ne n_elements(ar.galcat_gmem_cc) then message,'Error: Bad gmem indexing.'
  if offset_tr ne n_elements(ar.accSnap_gmem) then message,'Error: Bad tracer indexing.'
  
  ; convert (numPart / Gyr) to (msun (h^-1) / yr)
  r.gmem_accRate_halo *= (massPerPart * float(units.UnitMass_in_Msun)) / 1e9
  
  ; bin medians vs. halo mass
  for i=0,logMassNbins-1 do begin
  
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1],count)
    
    if count gt 0 then begin
      for j=0,nWindows-1 do begin
        r.gmem_accRate_median[j,i] = median(r.gmem_accRate_halo[j,w])
      endfor
    endif
    
  endfor
  
  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  return, r
end

; plotAccretionRates(): vs halo mass

pro plotAccretionRates
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  sP = simParams(res=128,run='tracer',redshift=2.0)
  
  if sP.run eq 'tracer'   then snapsBackMinusOne = 8 ; back 9 steps, ~300Myr on 313 snaps
  if sP.run eq 'feedback' then snapsBackMinusOne = 1 ; back 2 snaps, ~300Myr on 139 snaps
  if n_elements(snapsBackMinusOne) eq 0 then message,'Error'
  
  ; load
  bAr = binAccretionRates(sP=sP,sgSelect='pri')
  
  ; plot
  sK     = 5 ; smoothing kernel
  xrange = [9.0,12.0]
  yrange = [0.01,50.0]

  start_PS, sP.plotPath + 'accrate_gmem.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/ylog,/xs,/ys,yminor=0,$
      xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
      ytitle=textoidl("Gas Accretion Rate [M_{sun } h^{-1} yr^{-1}]")
         
    ; set a contour palette
    loadColorTable,'brewerc-blues'
    tvlct, rr, gg, bb, /get
    palette = [[rr],[gg],[bb]]
    
    ; contour
    cgd = calcGridData(xx=bAr.gcMasses,yy=bAr.gmem_accRate_halo[snapsBackMinusOne,*],$
      xMinMax=xrange,yMinMax=yrange,nPixels=[50,50],/logY)
    
    cgContour, smooth(cgd.dens_out,[3,3]), cgd.xPts, 10.0^cgd.yPts, $
      /overplot, /fill, palette=palette, levels=[0.25,0.5,1,2,5,10,20],c_colors=(indgen(7)*20+50)
    ;cgPlot,bAr.gcMasses,bAr.gmem_accRate_halo[1,*],psym=4,color=cgColor('orange'),/overplot
    
    ; median lines
    loadct,0,/silent
    
    cgPlot,bAr.logMassBinCen,smooth(bAr.gmem_accRate_median[snapsBackMinusOne,*],sK,/nan),line=j,/overplot
    
    ;for j=0,bAr.nWindows-1 do $
    ;  cgPlot,bAr.logMassBinCen,bAr.gmem_accRate_median[j,*],line=j,/overplot
      
    ; legend, redo plot borders
    ;legend,string(indgen(bAr.nWindows)+1,format='(i2)'),linestyle=indgen(bAr.nWindows),box=0,/top,/left
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/ylog,yminor=0,/xs,/ys,/noerase
    
  end_PS
  stop
end


