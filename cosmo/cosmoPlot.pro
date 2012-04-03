; cosmoPlot.pro
; gas accretion project - plots
; dnelson apr.2012

; plotColdFracVsHaloMass(): plot the "cold mode fraction" vs halo mass in a few different ways

pro plotColdFracVsHaloMass, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect = 'pri'
  TcutVals = [5.3,5.4,5.5,5.6,5.7] ; for constant threshold
  nCuts = n_elements(TcutVals)
  
  minNum = 32
  xrange = [9.0,12.5]
  yrange = [0.0,1.1]
  
  logMassBinSize = 0.1
  
  ; make a uniform gas selection at the start
  at = accretionTimes(sP=sP)
  mt = mergerTreeSubset(sP=sP)
  
  gal_w  = where(at.AccTime_gal ne -1,count_gal)
  gmem_w = where(at.AccTime_gmem ne -1,count_gmem)
    
  ; reverse histogram parent IDs of all particles/tracers in this selection
  if sP.trMCPerCell eq 0 then begin
    hist_gal  = histogram(mt.gcIndOrig.gal[gal_w],min=0,loc=loc_gal,rev=rev_gal)
    hist_gmem = histogram(mt.gcIndOrig.gmem[gmem_w],min=0,loc=loc_gmem,rev=rev_gmem)
  endif else begin
    ; create a gcIndOrig for the tracers
    galcat = galaxyCat(sP=sP)
    
    gcIndOrigTr = galCatRepParentIDs(galcat=galcat,gcIDList=mt.galcatIDList,$
                    child_counts={gal:at.child_counts_gal,gmem:at.child_counts_gmem}) 

    galcat = !NULL     
    
    ; want to use these parent IDs to access mergerTreeSubset so compact the same way (ascending ID->index)
    placeMap = getIDIndexMap(mt.galcatIDList,minid=minid)
    gcIndOrigTr.gal = placeMap[gcIndOrigTr.gal-minid]
    gcIndOrigTr.gmem = placeMap[gcIndOrigTr.gmem-minid]
    placeMap = !NULL
    
    hist_gal  = histogram(gcIndOrigTr.gal[gal_w],min=0,loc=loc_gal,rev=rev_gal)
    hist_gmem = histogram(gcIndOrigTr.gmem[gmem_w],min=0,loc=loc_gmem,rev=rev_gmem)

    gcIndOrigTr = !NULL
  endelse

  ; load max temps, current tvir, tvir at accretion, and parent masses
  accTvir = gcSubsetProp(sP=sP,select=sgSelect,/accTvir,/mergerTreeSubset,/accretionTimeSubset)
  curTvir = gcSubsetProp(sP=sP,select=sgSelect,/virTemp,/mergerTreeSubset,/accretionTimeSubset)
  maxTemp = gcSubsetProp(sP=sP,select=sgSelect,/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset)
  parentMass  = gcSubsetProp(sP=sP,select=sgSelect,/parMass,/mergerTreeSubset,/accretionTimeSubset)

  ; load current temps and current SFR
  curTemp = gcSubsetProp(sP=sP,select=sgSelect,/curTemp,/mergerTreeSubset,/accretionTimeSubset)
  curSFR  = gcSubsetProp(sP=sP,select=sgSelect,/curSingleVal,singleValField='sfr',$
                         /mergerTreeSubset,/accretionTimeSubset)

  ; structures to store results (Tmax)
  coldFrac = { gal_const    : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
               gmem_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
               gal_tvircur  : fltarr(n_elements(mt.galcatIDList))         ,$
               gmem_tvircur : fltarr(n_elements(mt.galcatIDList))         ,$
               gal_tviracc  : fltarr(n_elements(mt.galcatIDList))         ,$
               gmem_tviracc : fltarr(n_elements(mt.galcatIDList))         ,$
               both_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
               both_tvircur : fltarr(n_elements(mt.galcatIDList))         ,$
               both_tviracc : fltarr(n_elements(mt.galcatIDList))         ,$
               gal_num      : lonarr(n_elements(mt.galcatIDList))         ,$
               gmem_num     : lonarr(n_elements(mt.galcatIDList))          }
  
  ; loop over all tracked subgroups (galaxy, Tmax)
  for i=0,n_elements(hist_gal)-1 do begin
    if hist_gal[i] gt 0 then begin
      ; list of indices of galaxy gas particles in this subgroup
      loc_inds_gal = rev_gal[rev_gal[i]:rev_gal[i+1]-1]
      loc_maxt_gal = maxTemp.gal[loc_inds_gal]
      nloc = n_elements(loc_maxt_gal)
      
      coldFrac.gal_num[i] = nloc
      
      ; count fraction Tmax below each constant temperature threshold
      for j=0,nCuts-1 do begin
        w = where(loc_maxt_gal le TcutVals[j],count_below)
        coldFrac.gal_const[j,i] = float(count_below) / nloc
      endfor
      
      ; count fraction Tmax below Tvir at current time
      w = where(loc_maxt_gal le curTvir.gal[loc_inds_gal],count_below)
      coldFrac.gal_tvircur[i] = float(count_below) / nloc
      
      ; count fraction Tmax below Tvir at accretion time
      w = where(loc_maxt_gal le accTvir.gal[loc_inds_gal],count_below)
      coldFrac.gal_tviracc[i] = float(count_below) / nloc
    endif
  endfor
  
  ; loop over all tracked subgroups (groupmem, Tmax)
  for i=0,n_elements(hist_gmem)-1 do begin
    if hist_gmem[i] gt 0 then begin
      ; list of indices of groupmem gas particles in this subgroup
      loc_inds_gmem = rev_gmem[rev_gmem[i]:rev_gmem[i+1]-1]
      loc_maxt_gmem = maxTemp.gmem[loc_inds_gmem]
      nloc = n_elements(loc_maxt_gmem)
      
      coldFrac.gmem_num[i] = nloc
      
      ; count fraction Tmax below each constant temperature threshold
      for j=0,nCuts-1 do begin
        w = where(loc_maxt_gmem le TcutVals[j],count_below)
        coldFrac.gmem_const[j,i] = float(count_below) / nloc
      endfor
      
      ; count fraction Tmax below Tvir at current time
      w = where(loc_maxt_gmem le curTvir.gmem[loc_inds_gmem],count_below)
      coldFrac.gmem_tvircur[i] = float(count_below) / nloc
      
      ; count fraction Tmax below Tvir at accretion time
      w = where(loc_maxt_gmem le accTvir.gmem[loc_inds_gmem],count_below)
      coldFrac.gmem_tviracc[i] = float(count_below) / nloc
    endif
  endfor
  
  ; create composite gal+gmem (Tmax)
  coldFrac.both_tvircur = (coldFrac.gal_tvircur * coldFrac.gal_num + $
                           coldFrac.gmem_tvircur * coldFrac.gmem_num) / $
                          (coldFrac.gal_num + coldFrac.gmem_num)
  coldFrac.both_tviracc = (coldFrac.gal_tviracc * coldFrac.gal_num + $
                           coldFrac.gmem_tviracc * coldFrac.gmem_num) / $
                          (coldFrac.gal_num + coldFrac.gmem_num)
  for j=0,nCuts-1 do $
    coldFrac.both_const[j,*] = (coldFrac.gal_const[j,*] * coldFrac.gal_num + $
                                coldFrac.gmem_const[j,*] * coldFrac.gmem_num) / $
                               (coldFrac.gal_num + coldFrac.gmem_num)  
  
  ; structures to store results (Tcur)
  coldFrac_cur = { gal_const    : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                   gmem_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                   gal_tvircur  : fltarr(n_elements(mt.galcatIDList))         ,$
                   gmem_tvircur : fltarr(n_elements(mt.galcatIDList))         ,$
                   gal_tviracc  : fltarr(n_elements(mt.galcatIDList))         ,$
                   gmem_tviracc : fltarr(n_elements(mt.galcatIDList))         ,$
                   both_const   : fltarr(nCuts,n_elements(mt.galcatIDList))   ,$
                   both_tvircur : fltarr(n_elements(mt.galcatIDList))         ,$
                   both_tviracc : fltarr(n_elements(mt.galcatIDList))         ,$
                   gal_num      : lonarr(n_elements(mt.galcatIDList))         ,$
                   gmem_num     : lonarr(n_elements(mt.galcatIDList))          }
  
  ; loop over all tracked subgroups (galaxy, Tcur)
  for i=0,n_elements(hist_gal)-1 do begin
    if hist_gal[i] gt 0 then begin
      ; list of indices of galaxy gas particles in this subgroup
      loc_inds_gal   = rev_gal[rev_gal[i]:rev_gal[i+1]-1]
      loc_curt_gal   = curTemp.gal[loc_inds_gal]
      loc_cursfr_gal = curSFR.gal[loc_inds_gal]
      
      ; select only non-eff eos gas
      w_sfr = where(loc_cursfr_gal ne 0.0,count_nosfr)
      if count_nosfr gt 0 then begin
        loc_curt_gal = loc_curt_gal[w_sfr] ; select
        loc_inds_gal = loc_inds_gal[w_sfr]
        
        coldFrac_cur.gal_num[i] = count_nosfr
      
        ; count fraction Tcur below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_curt_gal le TcutVals[j],count_below)
          coldFrac_cur.gal_const[j,i] = float(count_below) / count_nosfr
        endfor
        
        ; count fraction Tmax below Tvir at current time
        w = where(loc_curt_gal le curTvir.gal[loc_inds_gal],count_below)
        coldFrac_cur.gal_tvircur[i] = float(count_below) / count_nosfr
        
        ; count fraction Tmax below Tvir at accretion time
        w = where(loc_curt_gal le accTvir.gal[loc_inds_gal],count_below)
        coldFrac_cur.gal_tviracc[i] = float(count_below) / count_nosfr
      endif ; cursfr!=0
    endif ; hist_gal>0
  endfor  
  
  ; loop over all tracked subgroups (groupmem, Tcur)
  for i=0,n_elements(hist_gmem)-1 do begin
    if hist_gmem[i] gt 0 then begin
      ; list of indices of groupmem gas particles in this subgroup
      loc_inds_gmem   = rev_gmem[rev_gmem[i]:rev_gmem[i+1]-1]
      loc_curt_gmem   = curTemp.gmem[loc_inds_gmem]
      loc_cursfr_gmem = curSFR.gmem[loc_inds_gmem]
      
      ; select only non-eff eos gas
      w_sfr = where(loc_cursfr_gmem ne 0.0,count_nosfr)
      if count_nosfr gt 0 then begin
        loc_curt_gmem = loc_curt_gmem[w_sfr] ; select
        loc_inds_gmem = loc_inds_gmem[w_sfr]
        
        coldFrac_cur.gmem_num[i] = count_nosfr
      
        ; count fraction Tcur below each constant temperature threshold
        for j=0,nCuts-1 do begin
          w = where(loc_curt_gmem le TcutVals[j],count_below)
          coldFrac_cur.gmem_const[j,i] = float(count_below) / count_nosfr
        endfor
        
        ; count fraction Tmax below Tvir at current time
        w = where(loc_curt_gmem le curTvir.gmem[loc_inds_gmem],count_below)
        coldFrac_cur.gmem_tvircur[i] = float(count_below) / count_nosfr
        
        ; count fraction Tmax below Tvir at accretion time
        w = where(loc_curt_gmem le accTvir.gmem[loc_inds_gmem],count_below)
        coldFrac_cur.gmem_tviracc[i] = float(count_below) / count_nosfr
      endif ; cursfr!=0
    endif ; hist_gal>0
  endfor
                        
  ; create composite gal+gmem (Tcur)
  coldFrac_cur.both_tvircur = (coldFrac_cur.gal_tvircur * coldFrac_cur.gal_num + $
                               coldFrac_cur.gmem_tvircur * coldFrac_cur.gmem_num) / $
                              (coldFrac_cur.gal_num + coldFrac_cur.gmem_num)
  coldFrac_cur.both_tviracc = (coldFrac_cur.gal_tviracc * coldFrac_cur.gal_num + $
                               coldFrac_cur.gmem_tviracc * coldFrac_cur.gmem_num) / $
                              (coldFrac_cur.gal_num + coldFrac_cur.gmem_num)
  for j=0,nCuts-1 do $
    coldFrac_cur.both_const[j,*] = (coldFrac_cur.gal_const[j,*] * coldFrac_cur.gal_num + $
                                    coldFrac_cur.gmem_const[j,*] * coldFrac_cur.gmem_num) / $
                                   (coldFrac_cur.gal_num + coldFrac_cur.gmem_num)                       
                                 
  ; load group cat for subgroup masses
  gc = loadGroupCat(sP=sP,/skipIDs)
  gcMasses = codeMassToLogMsun(gc.subgroupMass[mt.galcatIDList])
  gc = !NULL
  
  ; bin fractions into halo mass bins and make median lines
  logMassNbins  = (xrange[1]-xrange[0]) / logMassBinSize
  logMassBins   = linspace(xrange[0],xrange[1],logMassNbins+1) ; edges
  logMassBinCen = linspace(xrange[0],xrange[1],logMassNbins+1) + logMassBinSize/2.0
  
  medianVals = { const_gal    : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 const_gmem   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 const_both   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                 tVircur_gal  : fltarr(logMassNbins)       + !values.f_nan ,$
                 tVircur_gmem : fltarr(logMassNbins)       + !values.f_nan ,$
                 tVircur_both : fltarr(logMassNbins)       + !values.f_nan ,$
                 tViracc_gal  : fltarr(logMassNbins)       + !values.f_nan ,$
                 tViracc_gmem : fltarr(logMassNbins)       + !values.f_nan ,$
                 tViracc_both : fltarr(logMassNbins)       + !values.f_nan  }
                 
  medianVals_cur = { const_gal    : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                     const_gmem   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                     const_both   : fltarr(nCuts,logMassNbins) + !values.f_nan ,$
                     tVircur_gal  : fltarr(logMassNbins)       + !values.f_nan ,$
                     tVircur_gmem : fltarr(logMassNbins)       + !values.f_nan ,$
                     tVircur_both : fltarr(logMassNbins)       + !values.f_nan ,$
                     tViracc_gal  : fltarr(logMassNbins)       + !values.f_nan ,$
                     tViracc_gmem : fltarr(logMassNbins)       + !values.f_nan ,$
                     tViracc_both : fltarr(logMassNbins)       + !values.f_nan  }
  
  ; calculate median in bins (Tmax) and enforce minimum particle numbers
  for i=0,logMassNbins-2 do begin
    ; gal (Tmax)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac.gal_num ge minNum,count)
    if count gt 0 then begin
      medianVals.tVircur_gal[i] = median(coldFrac.gal_tvircur[w])
      medianVals.tViracc_gal[i] = median(coldFrac.gal_tviracc[w])
      for j=0,nCuts-1 do medianVals.const_gal[j,i] = median(coldFrac.gal_const[j,w])
    endif
    
    ; gmem (Tmax)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac.gmem_num ge minNum,count)
    if count gt 0 then begin
      medianVals.tVircur_gmem[i] = median(coldFrac.gmem_tvircur[w])
      medianVals.tViracc_gmem[i] = median(coldFrac.gmem_tviracc[w])
      for j=0,nCuts-1 do medianVals.const_gmem[j,i] = median(coldFrac.gmem_const[j,w])
    endif
    
    ; both (Tmax)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac.gal_num+coldFrac.gmem_num ge minNum,count)
    if count gt 0 then begin
      medianVals.tVircur_both[i] = median(coldFrac.both_tvircur[w])
      medianVals.tViracc_both[i] = median(coldFrac.both_tviracc[w])
      for j=0,nCuts-1 do medianVals.const_both[j,i] = median(coldFrac.both_const[j,w])
    endif
    
    ; gal (Tcur)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac_cur.gal_num ge minNum,count)
    if count gt 0 then begin
      medianVals_cur.tVircur_gal[i] = median(coldFrac_cur.gal_tvircur[w])
      medianVals_cur.tViracc_gal[i] = median(coldFrac_cur.gal_tviracc[w])
      for j=0,nCuts-1 do medianVals_cur.const_gal[j,i] = median(coldFrac_cur.gal_const[j,w])
    endif
    
    ; gmem (Tcur)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac_cur.gmem_num ge minNum,count)
    if count gt 0 then begin
      medianVals_cur.tVircur_gmem[i] = median(coldFrac_cur.gmem_tvircur[w])
      medianVals_cur.tViracc_gmem[i] = median(coldFrac_cur.gmem_tviracc[w])
      for j=0,nCuts-1 do medianVals_cur.const_gmem[j,i] = median(coldFrac_cur.gmem_const[j,w])
    endif
    
    ; both (Tcur)
    w = where(gcMasses gt logMassBins[i] and gcMasses le logMassBins[i+1] and $
              coldFrac_cur.gal_num+coldFrac_cur.gmem_num ge minNum,count)
    if count gt 0 then begin
      medianVals_cur.tVircur_both[i] = median(coldFrac_cur.both_tvircur[w])
      medianVals_cur.tViracc_both[i] = median(coldFrac_cur.both_tviracc[w])
      for j=0,nCuts-1 do medianVals_cur.const_both[j,i] = median(coldFrac_cur.both_const[j,w])
    endif
  endfor
  
  ; plot (1) - all data points and median lines (Tmax)
  start_PS, sP.plotPath + 'coldFrac.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("log(M_{halo})")+"",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+")"
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; plot raw points
    psym = 4
    symsize = 0.5
    
    w = where(coldFrac.gal_num+coldFrac.gmem_num ge minNum,count)
    print,count,n_elements(coldFrac.gal_num)
    
    cgPlot,gcMasses,coldFrac.both_tVircur[w],psym=psym,symsize=symsize,thick=1.0,$
      color=getColor(1),/overplot
    cgPlot,gcMasses,coldFrac.both_tViracc[w],psym=psym,symsize=symsize,thick=1.0,$
      color=getColor(2),/overplot
    
    for j=0,nCuts-1 do $
      cgPlot,gcMasses,coldFrac.both_const[j,w],psym=psym,symsize=symsize,thick=1.0,$
        color=getColor(3+j),/overplot
        
    ; plot median lines
    cgPlot,logMassBinCen,medianVals.tVircur_both,color=getColor(1),/overplot
    cgPlot,logMassBinCen,medianVals.tViracc_both,color=getColor(2),/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen,medianVals.const_both[j,*],color=getColor(3+j),/overplot
    
    ; legend
    strings = [textoidl("T_{vir,cur}"),textoidl("T_{vir,acc}"),$
               textoidl("T_{cut} = ")+string(TcutVals,format='(f4.1)')]
    colors  = getColor([[1,2],indgen(nCuts)+3],/name)
    legend,strings,textcolor=colors,box=0,/bottom,/left
    
  end_PS
  
  ; plot (2) - just median lines (Tmax)
  start_PS, sP.plotPath + 'coldFrac2.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("log(M_{halo})")+"",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+")"
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; plot raw points
    psym = 4
    symsize = 0.5
    
    ; gal
    cgPlot,logMassBinCen,medianVals.tVircur_gal,color=getColor(1),line=1,/overplot
    cgPlot,logMassBinCen,medianVals.tViracc_gal,color=getColor(2),line=1,/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen,medianVals.const_gal[j,*],color=getColor(3+j),line=1,/overplot
      
    ; gmem
    cgPlot,logMassBinCen,medianVals.tVircur_gmem,color=getColor(1),line=2,/overplot
    cgPlot,logMassBinCen,medianVals.tViracc_gmem,color=getColor(2),line=2,/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen,medianVals.const_gmem[j,*],color=getColor(3+j),line=2,/overplot 
    
    ; both
    cgPlot,logMassBinCen,medianVals.tVircur_both,color=getColor(1),/overplot
    cgPlot,logMassBinCen,medianVals.tViracc_both,color=getColor(2),/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen,medianVals.const_both[j,*],color=getColor(3+j),/overplot
    
    ; legend
    strings = [textoidl("T_{vir,cur}"),textoidl("T_{vir,acc}"),$
               textoidl("T_{cut} = ")+string(TcutVals,format='(f4.1)')]
    colors  = getColor([[1,2],indgen(nCuts)+3],/name)
    legend,strings,textcolor=colors,box=0,/bottom,/left
    legend,['gal','gmem','both'],linestyle=[1,2,0],box=0,linesize=0.25,/top,/left
    
  end_PS
  
  ; plot (3) - all data points and median lines (Tcur)
  start_PS, sP.plotPath + 'coldFrac_tcur.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("log(M_{halo})")+"",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+")"
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; plot raw points
    psym = 4
    symsize = 0.5
    
    w = where(coldFrac_cur.gal_num+coldFrac_cur.gmem_num ge minNum,count)
    print,count,n_elements(coldFrac_cur.gal_num)
    
    cgPlot,gcMasses,coldFrac_cur.both_tVircur[w],psym=psym,symsize=symsize,thick=1.0,$
      color=getColor(1),/overplot
    cgPlot,gcMasses,coldFrac_cur.both_tViracc[w],psym=psym,symsize=symsize,thick=1.0,$
      color=getColor(2),/overplot
    
    for j=0,nCuts-1 do $
      cgPlot,gcMasses,coldFrac_cur.both_const[j,w],psym=psym,symsize=symsize,thick=1.0,$
        color=getColor(3+j),/overplot
        
    ; plot median lines
    cgPlot,logMassBinCen,medianVals_cur.tVircur_both,color=getColor(1),/overplot
    cgPlot,logMassBinCen,medianVals_cur.tViracc_both,color=getColor(2),/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen,medianVals_cur.const_both[j,*],color=getColor(3+j),/overplot
    
    ; legend
    strings = [textoidl("T_{vir,cur}"),textoidl("T_{vir,acc}"),$
               textoidl("T_{cut} = ")+string(TcutVals,format='(f4.1)')]
    colors  = getColor([[1,2],indgen(nCuts)+3],/name)
    legend,strings,textcolor=colors,box=0,/bottom,/left
    
  end_PS
  
  ; plot (4) - just median lines (Tcur)
  start_PS, sP.plotPath + 'coldFrac_tcur2.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("log(M_{halo})")+"",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+")"
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; plot raw points
    psym = 4
    symsize = 0.5
    
    ; gal
    cgPlot,logMassBinCen,medianVals_cur.tVircur_gal,color=getColor(1),line=1,/overplot
    cgPlot,logMassBinCen,medianVals_cur.tViracc_gal,color=getColor(2),line=1,/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen,medianVals_cur.const_gal[j,*],color=getColor(3+j),line=1,/overplot
      
    ; gmem
    cgPlot,logMassBinCen,medianVals_cur.tVircur_gmem,color=getColor(1),line=2,/overplot
    cgPlot,logMassBinCen,medianVals_cur.tViracc_gmem,color=getColor(2),line=2,/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen,medianVals_cur.const_gmem[j,*],color=getColor(3+j),line=2,/overplot 
    
    ; both
    cgPlot,logMassBinCen,medianVals_cur.tVircur_both,color=getColor(1),/overplot
    cgPlot,logMassBinCen,medianVals_cur.tViracc_both,color=getColor(2),/overplot
    for j=0,nCuts-1 do $
      cgPlot,logMassBinCen,medianVals_cur.const_both[j,*],color=getColor(3+j),/overplot
    
    ; legend
    strings = [textoidl("T_{vir,cur}"),textoidl("T_{vir,acc}"),$
               textoidl("T_{cut} = ")+string(TcutVals,format='(f4.1)')]
    colors  = getColor([[1,2],indgen(nCuts)+3],/name)
    legend,strings,textcolor=colors,box=0,/bottom,/left
    legend,['gal','gmem','both'],linestyle=[1,2,0],box=0,linesize=0.25,/top,/left
    
  end_PS
end

; plotTmaxVsTvirAccComp(); plot the previous max temp vs. the virial temperature of the parent halos at the
;                          time of accretion for arepo vs. gadget

pro plotTmaxVsTvirAccComp

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  res = 256
  redshift = 3.0

  sP1 = simParams(res=res,run='gadget',redshift=redshift)
  sP2 = simParams(res=res,run='tracer',redshift=redshift)
  
  binSizeLog = 0.1 / (res/128)
  
  sgSelect = 'pri'
  ;massBins = [0.0,100.0] ; no massbins
  massBins = [9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5] ; log(M)   
  
  ; load sP1 (gadget)
  accTvir_gadget = gcSubsetProp(sP=sP1,select=sgSelect,/accTvir,/mergerTreeSubset,/accretionTimeSubset)
  curTvir_gadget = gcSubsetProp(sP=sP1,select=sgSelect,/virTemp,/mergerTreeSubset,/accretionTimeSubset)
  maxTemp_gadget = gcSubsetProp(sP=sP1,select=sgSelect,/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset)  

  ; load parent halo masses so we can make halo massbins
  parentMass_ga = gcSubsetProp(sP=sP1,select=sgSelect,/parMass,/mergerTreeSubset,/accretionTimeSubset)

  ; load sP2 (tracer)
  accTvir_tracer = gcSubsetProp(sP=sP2,select=sgSelect,/accTvir,/mergerTreeSubset,/accretionTimeSubset)
  curTvir_tracer = gcSubsetProp(sP=sP2,select=sgSelect,/virTemp,/mergerTreeSubset,/accretionTimeSubset)
  maxTemp_tracer = gcSubsetProp(sP=sP2,select=sgSelect,/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset)

  parentMass_tr = gcSubsetProp(sP=sP2,select=sgSelect,/parMass,/mergerTreeSubset,/accretionTimeSubset)
    
  ; plot (1) - compare tmax to tvir at time of accretion
  start_PS, sP1.plotPath + 'tmax_tviracc_comp.'+str(sP1.res)+'_'+str(sP1.snap)+'.eps'
    !p.thick += 1
    xrange = [-2.5,1.0]
    yrange = [1e-3,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="fraction",xtitle=textoidl("log( T_{max} / T_{vir,acc} )")+"",$
      title=str(res)+textoidl("^3")+" Gadget vs. ArepoMC (z="+string(sP1.redshift,format='(f3.1)')+")"
    cgPlot,[0,0],yrange,line=0,color=fsc_color('light gray'),/overplot
    
  strings = []
  colors = []
  
  for j=0,n_elements(massBins)-2 do begin
    
    ; select members of this parent mass bins and r>0<inf
    wGadget_gal  = where(parentMass_ga.gal gt massBins[j] and parentMass_ga.gal le massBins[j+1],count1)
    wGadget_gmem = where(parentMass_ga.gmem gt massBins[j] and parentMass_ga.gmem le massBins[j+1],count2)
    wTracer_gal  = where(parentMass_tr.gal gt massBins[j] and parentMass_tr.gal le massBins[j+1],count3)
    wTracer_gmem = where(parentMass_tr.gmem gt massBins[j] and parentMass_tr.gmem le massBins[j+1],count4)
    
    print,j,count1,count2,count3,count4
    
    if ~count1 or ~count2 or ~count3 or ~count4 then continue ; no halos in this mass bin  
  
    massBinStr = string(massBins[j],format='(f4.1)') + '-' + string(massBins[j+1],format='(f4.1)')
    strings = [strings,massBinStr]
    colors = [colors,getColor(j,/name)]
    
    ; histogram gadget (gal+gmem) differences
    vals = [10.0^maxTemp_gadget.gal[wGadget_gal]/10.0^accTvir_gadget.gal[wGadget_gal],$
            10.0^maxTemp_gadget.gmem[wGadget_gmem]/10.0^accTvir_gadget.gmem[wGadget_gmem]]
    hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=0,color=getColor(j),/overplot

    ; histogram tracer (gal+gmem) differences
    vals = [10.0^maxTemp_tracer.gal[wTracer_gal]/10.0^accTvir_tracer.gal[wTracer_gal],$
            10.0^maxTemp_tracer.gmem[wTracer_gmem]/10.0^accTvir_tracer.gmem[wTracer_gmem]]
    hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=2,color=getColor(j),/overplot
  
  endfor
  
  ; legend
  legend,['gadget','arepoMC'],linestyle=[0,2],linesize=0.25,box=0,/right,/top
  if n_elements(massBins) gt 2 then $
    legend,strings,textcolors=colors,box=0,/left,/top
  
  end_PS
  
  ; plot (2) - same plot using current Tvir instead of Tvir at accretion time
  start_PS, sP1.plotPath + 'tmax_tvircur_comp.'+str(sP1.res)+'_'+str(sP1.snap)+'.eps'
    !p.thick += 1
    xrange = [-2.5,1.0]
    yrange = [1e-3,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="fraction",xtitle=textoidl("log ( T_{max} / T_{vir,cur} )")+"",$
      title=str(res)+textoidl("^3")+" Gadget vs. ArepoMC (z="+string(sP1.redshift,format='(f3.1)')+")"
    cgPlot,[1,1],yrange,line=0,color=fsc_color('light gray'),/overplot
    
  strings = []
  colors = []
  
  for j=0,n_elements(massBins)-2 do begin
    
    ; select members of this parent mass bins and r>0<inf
    wGadget_gal  = where(parentMass_ga.gal gt massBins[j] and parentMass_ga.gal le massBins[j+1],count1)
    wGadget_gmem = where(parentMass_ga.gmem gt massBins[j] and parentMass_ga.gmem le massBins[j+1],count2)
    wTracer_gal  = where(parentMass_tr.gal gt massBins[j] and parentMass_tr.gal le massBins[j+1],count3)
    wTracer_gmem = where(parentMass_tr.gmem gt massBins[j] and parentMass_tr.gmem le massBins[j+1],count4)
    
    print,j,count1,count2,count3,count4
    
    if ~count1 or ~count2 or ~count3 or ~count4 then continue ; no halos in this mass bin  
  
    massBinStr = string(massBins[j],format='(f4.1)') + '-' + string(massBins[j+1],format='(f4.1)')
    strings = [strings,massBinStr]
    colors = [colors,getColor(j,/name)]
    
    ; histogram gadget (gal+gmem) differences
    vals = [10.0^maxTemp_gadget.gal[wGadget_gal]/10.0^curTvir_gadget.gal[wGadget_gal],$
            10.0^maxTemp_gadget.gmem[wGadget_gmem]/10.0^curTvir_gadget.gmem[wGadget_gmem]]
    hist = histogram(alog10(vals),binsize=binSizeLog,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=0,color=getColor(j),/overplot

    ; histogram tracer (gal+gmem) differences
    vals = [10.0^maxTemp_tracer.gal[wTracer_gal]/10.0^curTvir_tracer.gal[wTracer_gal],$
            10.0^maxTemp_tracer.gmem[wTracer_gmem]/10.0^curTvir_tracer.gmem[wTracer_gmem]]
    hist = histogram(alog10(vals),binsize=binSizeLog,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=2,color=getColor(j),/overplot
  
  endfor
  
  ; legend
  legend,['gadget','arepoMC'],linestyle=[0,2],linesize=0.25,box=0,/right,/top
  if n_elements(massBins) gt 2 then $
    legend,strings,textcolors=colors,box=0,/left,/top
  
  end_PS
  
  ; plot (3) - unscaled tmax
  start_PS, sP1.plotPath + 'tmax_nonorm.'+str(sP1.res)+'_'+str(sP1.snap)+'.eps'
    !p.thick += 1
    xrange = [4.0,7.0]
    yrange = [1e-3,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="fraction",xtitle=textoidl("log (T_{max})")+"",$
      title=str(res)+textoidl("^3")+" Gadget vs. ArepoMC (z="+string(sP1.redshift,format='(f3.1)')+")"
    cgPlot,[1,1],yrange,line=0,color=fsc_color('light gray'),/overplot
    
  strings = []
  colors = []
  
  for j=0,n_elements(massBins)-2 do begin
    
    ; select members of this parent mass bins and r>0<inf
    wGadget_gal  = where(parentMass_ga.gal gt massBins[j] and parentMass_ga.gal le massBins[j+1],count1)
    wGadget_gmem = where(parentMass_ga.gmem gt massBins[j] and parentMass_ga.gmem le massBins[j+1],count2)
    wTracer_gal  = where(parentMass_tr.gal gt massBins[j] and parentMass_tr.gal le massBins[j+1],count3)
    wTracer_gmem = where(parentMass_tr.gmem gt massBins[j] and parentMass_tr.gmem le massBins[j+1],count4)
    
    print,j,count1,count2,count3,count4
    
    if ~count1 or ~count2 or ~count3 or ~count4 then continue ; no halos in this mass bin  
  
    massBinStr = string(massBins[j],format='(f4.1)') + '-' + string(massBins[j+1],format='(f4.1)')
    strings = [strings,massBinStr]
    colors = [colors,getColor(j,/name)]
    
    ; histogram gadget (gal+gmem) differences
    vals = [maxTemp_gadget.gal[wGadget_gal],maxTemp_gadget.gmem[wGadget_gmem]]
    hist = histogram(vals,binsize=binSizeLog,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=0,color=getColor(j),/overplot

    ; histogram tracer (gal+gmem) differences
    vals = [maxTemp_tracer.gal[wTracer_gal],maxTemp_tracer.gmem[wTracer_gmem]]
    hist = histogram(vals,binsize=binSizeLog,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=2,color=getColor(j),/overplot
  
  endfor
  
  ; legend
  legend,['gadget','arepoMC'],linestyle=[0,2],linesize=0.25,box=0,/right,/top
  if n_elements(massBins) gt 2 then $
    legend,strings,textcolors=colors,box=0,/left,/top
  
  end_PS
  
end

; plotTmaxVsTvirAccCur(): evaluate how the ratio of Tmax/Tvir changes when using either the current
;                         Tvir of the halo or the Tvir of the halo at the time of accretion

pro plotTmaxVsTvirAccCur, sP=sP
  
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  binSizeLog = 0.1 / (sP.res/128)
  
  sgSelect = 'pri'
  ;massBins = [0.0,100.0] ; no massbins
  massBins = [9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5] ; log(M)   
  
  ; load
  accTvir = gcSubsetProp(sP=sP,select=sgSelect,/accTvir,/mergerTreeSubset,/accretionTimeSubset)
  curTvir = gcSubsetProp(sP=sP,select=sgSelect,/virTemp,/mergerTreeSubset,/accretionTimeSubset)
  maxTemp = gcSubsetProp(sP=sP,select=sgSelect,/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset)  

  ; load parent halo masses so we can make halo massbins
  parentMass = gcSubsetProp(sP=sP,select=sgSelect,/parMass,/mergerTreeSubset,/accretionTimeSubset)  
  
  ; plot (1)
  start_PS, sP.plotPath + 'tmax_tvircur_tviracc.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    !p.thick += 1
    xrange = [-2.5,1.0]
    yrange = [1e-3,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="fraction",xtitle=textoidl("log ( T_{max} / T_{vir} )")+"",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+")"
    cgPlot,[0,0],yrange,line=0,color=fsc_color('light gray'),/overplot
    
  strings = []
  colors = []
  
  for j=0,n_elements(massBins)-2 do begin
    
    ; select members of this parent mass bins and r>0<inf
    w_gal  = where(parentMass.gal gt massBins[j] and parentMass.gal le massBins[j+1],count1)
    w_gmem = where(parentMass.gmem gt massBins[j] and parentMass.gmem le massBins[j+1],count2)
    
    print,j,count1,count2
    if ~count1 or ~count2 then continue ; no halos in this mass bin  
  
    massBinStr = string(massBins[j],format='(f4.1)') + '-' + string(massBins[j+1],format='(f4.1)')
    strings = [strings,massBinStr]
    colors = [colors,getColor(j,/name)]
    
    ; histogram gadget (gal+gmem) differences for current Tvir
    vals = [10.0^maxTemp.gal[w_gal]/10.0^accTvir.gal[w_gal],$
            10.0^maxTemp.gmem[w_gmem]/10.0^accTvir.gmem[w_gmem]]
    hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=0,color=getColor(j),/overplot

    ; histogram gadget (gal+gmem) differences for accretion Tvir
    vals = [10.0^maxTemp.gal[w_gal]/10.0^curTvir.gal[w_gal],$
            10.0^maxTemp.gmem[w_gmem]/10.0^curTvir.gmem[w_gmem]]
    hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=2,color=getColor(j),/overplot
  
  endfor
  
  ; legend
  legend,[textoidl("T_{vir,acc}"),textoidl("T_{vir,cur}")],linestyle=[0,2],linesize=0.25,box=0,/left,/top
  if n_elements(massBins) gt 2 then $
    legend,strings,textcolors=colors,box=0,/right,/top
  
  end_PS
  
  ; plot (2) - just tviracc but separate out gal and gmem
  start_PS, sP.plotPath + 'tmax_tviracc_gal_gmem.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    !p.thick += 1
    xrange = [-2.5,1.0]
    yrange = [1e-3,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="fraction",xtitle=textoidl("log ( T_{max} / T_{vir} )")+"",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+")"
    cgPlot,[0,0],yrange,line=0,color=fsc_color('light gray'),/overplot
    
  strings = []
  colors = []
  
  for j=0,n_elements(massBins)-2 do begin
    
    ; select members of this parent mass bins and r>0<inf
    w_gal  = where(parentMass.gal gt massBins[j] and parentMass.gal le massBins[j+1],count1)
    w_gmem = where(parentMass.gmem gt massBins[j] and parentMass.gmem le massBins[j+1],count2)
    
    print,j,count1,count2
    if ~count1 or ~count2 then continue ; no halos in this mass bin  
  
    massBinStr = string(massBins[j],format='(f4.1)') + '-' + string(massBins[j+1],format='(f4.1)')
    strings = [strings,massBinStr]
    colors = [colors,getColor(j,/name)]
    
    ; histogram gadget (both) differences for current Tvir
    vals = [10.0^maxTemp.gal[w_gal]/10.0^accTvir.gal[w_gal],$
            10.0^maxTemp.gmem[w_gmem]/10.0^accTvir.gmem[w_gmem]]
    hist_both = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
    cgPlot,loc,float(hist_both)/total(hist_both),line=0,color=getColor(j),/overplot
    
    ; histogram gadget (gal) differences for current Tvir
    vals = [10.0^maxTemp.gal[w_gal]/10.0^accTvir.gal[w_gal]]
    hist_gal = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
    cgPlot,loc,float(hist_gal)/total(hist_both),line=1,color=getColor(j),/overplot
  
    ; histogram gadget (gmem) differences for current Tvir
    vals = [10.0^maxTemp.gmem[w_gmem]/10.0^accTvir.gmem[w_gmem]]
    hist_gmem = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
    cgPlot,loc,float(hist_gmem)/total(hist_both),line=2,color=getColor(j),/overplot
  
  endfor
  
  ; legend
  legend,['gal','gmem','both'],linestyle=[1,2,0],linesize=0.25,box=0,/left,/top
  if n_elements(massBins) gt 2 then $
    legend,strings,textcolors=colors,box=0,/right,/top
  
  end_PS
  stop
end

; plotAccTimeVsTmaxTime(): plot the offset between the time of accretion and the time of max prev. temp

pro plotAccTimeVsTmaxTime, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  binsizeGyr = 0.15 / (sP.res/128)

  ; config
  sgSelect = 'pri'
  ;massBins = [0.0,100.0] ; no massbins
  massBins = [9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5] ; log(M)  
  
  accTime = gcSubsetProp(sP=sP,select=sgSelect,/accTime,/mergerTreeSubset,/accretionTimeSubset)
  maxTempTime = gcSubsetProp(sP=sP,select=sgSelect,/maxTempTime,/mergerTreeSubset,/accretionTimeSubset)
  
  ; convert redshifts to age of universe
  accTime.gal      = redshiftToAgeFlat(accTime.gal)
  accTime.gmem     = redshiftToAgeFlat(accTime.gmem)
  maxTempTime.gal  = redshiftToAgeFlat(maxTempTime.gal)
  maxTempTime.gmem = redshiftToAgeFlat(maxTempTime.gmem)
    
  ; load parent halo masses so we can make halo massbins
  parentMass = gcSubsetProp(sP=sP,select=sgSelect,/parMass,/mergerTreeSubset,/accretionTimeSubset)

  ; plot
  start_PS, sP.plotPath + 'acctime_tmaxtime_' + sP.run + '.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    xrange = [-1,1.5]
    yrange = [1e-3,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="fraction",xtitle=textoidl("t_{acc} - t_{T_{max}}")+" [Gyr]",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+$
      ") gal=solid gmem=dot"
    cgPlot,[0,0],yrange,line=0,color=fsc_color('light gray'),/overplot
    
  strings = []
  colors = []
  
  for j=0,n_elements(massBins)-2 do begin
    
    ; select members of this parent mass bins and r>0<inf
    wGal  = where(parentMass.gal gt massBins[j] and parentMass.gal le massBins[j+1],count1)
    wGmem = where(parentMass.gmem gt massBins[j] and parentMass.gmem le massBins[j+1],count2)
    
    print,j,count1,count2
    
    if count1 eq 0 or count2 eq 0 then continue ; no halos in this mass bin  
  
    massBinStr = string(massBins[j],format='(f4.1)') + '-' + string(massBins[j+1],format='(f4.1)')
    strings = [strings,massBinStr]
    colors = [colors,j]
  
    ; histogram gal differences
    hist = histogram(accTime.gal[wGal]-maxTempTime.gal[wGal],binsize=binsizeGyr,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=0,color=getColor(j),/overplot

    ; histogram gmem differences
    hist = histogram(accTime.gmem[wGmem]-maxTempTime.gmem[wGmem],binsize=binsizeGyr,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=1,color=getColor(j),/overplot
  
  endfor
  
  ; legend
  legend,strings,textcolor=getColor(colors,/name),box=0,/right,/top
  
  end_PS
end

; plotVerticalSlices(): helper called by the subsequent 4 plot routines
  
pro plotVerticalSlices,rad_gal,rad_gmem,temp_gal,temp_gmem,plotName,binSizeTemp,temprange,xtitle

  compile_opt idl2, hidden, strictarr, strictarrsubs
  !except = 0 ;suppress floating point underflow/overflow errors
  
  ; config
  radBinEdges = [0.0,0.1,0.2,0.5,0.8,1.0,1.5] ;quasi-log spacing
  
  start_PS, plotName, xs=8, ys=6
  
    !p.multi = [0,3,2]
    
    !x.margin -= [4.0,1.5]
    !y.margin -= [0.5,0.5]
  
    for i=0,n_elements(radBinEdges)-2 do begin
      radBin = [radBinEdges[i],radBinEdges[i+1]]
      fsc_plot,[0],[0],/nodata,xrange=temprange,yrange=[0,1],/xs,/ys,ytitle="",xtitle=xtitle               
    
      ; select in radial bin
      w = where(rad_gal ge radBin[0] and rad_gal lt radBin[1],count_gal)
      if count_gal gt 0 then $  
        h1_gal = histogram(temp_gal[w],binsize=binSizeTemp,min=temprange[0],max=temprange[1],locations=xpts)
                          
      w = where(rad_gmem ge radBin[0] and rad_gmem lt radBin[1],count_gmem)
      if count_gmem gt 0 then $                 
        h1_gmem = histogram(temp_gmem[w],binsize=binSizeTemp,min=temprange[0],max=temprange[1],locations=xpts)

      ; skip when we have no points in this bin
      if ~count_gal and ~count_gmem then continue
      if n_elements(h1_gal) eq 0 then h1_gal = xpts*0
      if n_elements(h1_gmem) eq 0 then h1_gmem = xpts*0

      ; move xpts to bin centers and normalize histograms
      xpts += binSizeTemp/2.0
      
      h1_both = h1_gal + h1_gmem
      normfac = 1.0 / max(h1_both)
      
      h1_gal  *= normfac
      h1_gmem *= normfac
      h1_both *= normfac
      
      ; plot
      fsc_plot,xpts,h1_gal,color=getColor(16),thick=!p.thick+2,/overplot
      fsc_plot,xpts,h1_gmem,color=getColor(17),thick=!p.thick+2,/overplot
      fsc_plot,xpts,h1_both,color=getColor(0),line=1,thick=!p.thick+2,/overplot
      
      binStr = string(radBin[0],format='(f3.1)') + "-" + string(radBin[1],format='(f3.1)')
      fsc_text,(temprange[1]-temprange[0])*0.82,0.9,binStr,alignment=0.5,charsize=!p.charsize-0.5
    endfor
  
  !p.multi = 0
  
  end_PS
  !except = 1
end

; plot2DRadHistos(): helper function used by the subsequent 4 plot routines
  
pro plot2DRadHistos,plotBase,sP,h2rt_gal,h2rt_gmem,xrange,yrange,ytitle,$
                    virTempRange=virTempRange,massBinStr=massBinStr

  compile_opt idl2, hidden, strictarr, strictarrsubs
  !except = 0 ;suppress floating point underflow/overflow errors
  
  if sP.trMCPerCell eq 0 then  titleName = 'gadget'
  if sP.trMCPerCell eq -1 then titleName = 'tracerVel'
  if sP.trMCPerCell gt 0 then  titleName = 'tracerMC'
  
  h2rt_both = h2rt_gal + h2rt_gmem
  
  exp = 0.5 ; gamma exponent for non-linear color scaling
  ndivs = 5 ; number of divisions on colorbar  
  
  redshift = snapNumToRedshift(sP=sP)

  start_PS, sP.plotPath + plotBase + '_rad_both.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    loadct, 33, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
    ;cgLoadCT,13,/brewer,bottom=1
    ;loadColorTable, 'helix'
    
    tvim,h2rt_both^exp,pcharsize=!p.charsize,scale=0,clip=-1,$;,/c_map,range=[5e10,1e8,5e9];,/rct
         xtitle="r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}"),ytitle=ytitle,$
         title=str(sP.res)+textoidl("^3")+" "+titleName+" (z="+string(redshift,format='(f3.1)')+") gal+gmem",$
         barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=xrange,yrange=yrange,$
         /xlog,xticks=6,xtickv=[0.01,0.05,0.1,0.2,0.5,1.0,3.0],$
         xtickname=['0.01','0.05','0.1','0.2','0.5','1','3'],xmargin=2.0

     barvals = findgen(ndivs+1)/ndivs*(max(h2rt_both^exp)-min(h2rt_both^exp)) + min(h2rt_both^exp)
     ticknames = textoidl('10^{'+str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)'))+'}')
     ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
     colorbar,bottom=1,range=minmax(h2rt_both),position=[0.83,0.155,0.87,0.925],$
       /vertical,/right,title=textoidl("M_{sun,tot}"),divisions=ndivs,ticknames=ticknames,ncolors=255
            
      if keyword_set(virTempRange) then begin
        fsc_text,xrange[1]*0.45,(yrange[1]-yrange[0])*0.92+yrange[0],massBinStr,alignment=0.5,color=fsc_color('white')
        fsc_text,xrange[0]*1.6,virTempRange[1]*1.02,"T"+textoidl("_{vir}"),alignment=0.5,color=fsc_color('yellow')
        fsc_plot,xrange,[virTempRange[0],virTempRange[0]],line=0,color=fsc_color('yellow'),/overplot
        fsc_plot,xrange,[virTempRange[1],virTempRange[1]],line=0,color=fsc_color('yellow'),/overplot
      endif
             
  end_PS
  
  start_PS, sP.plotPath + plotBase + '_rad_gal.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    loadct, 33, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
    ;cgLoadCT,13,/brewer,bottom=1
    ;loadColorTable, 'helix'
    
    tvim,h2rt_gal^exp,pcharsize=!p.charsize,scale=0,clip=-1,$;,/c_map
         xtitle="r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}"),ytitle=ytitle,$
         title=str(sP.res)+textoidl("^3")+" "+titleName+" (z="+string(redshift,format='(f3.1)')+") gal",$
         barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=xrange,yrange=yrange,$
         /xlog,xticks=6,xtickv=[0.01,0.05,0.1,0.2,0.5,1.0,3.0],$
         xtickname=['0.01','0.05','0.1','0.2','0.5','1','3'],xmargin=2.0

     barvals = findgen(ndivs+1)/ndivs*(max(h2rt_gal^exp)-min(h2rt_gal^exp)) + min(h2rt_gal^exp)
     ticknames = textoidl('10^{'+str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)'))+'}')
     ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
     colorbar,bottom=1,range=minmax(h2rt_gal),position=[0.83,0.155,0.87,0.925],$
       /vertical,/right,title=textoidl("M_{sun,tot}"),divisions=ndivs,ticknames=ticknames,ncolors=255
         
      if keyword_set(virTempRange) then begin
        fsc_text,xrange[1]*0.45,(yrange[1]-yrange[0])*0.92+yrange[0],massBinStr,alignment=0.5,color=fsc_color('white')
        fsc_text,xrange[0]*1.6,virTempRange[1]*1.02,"T"+textoidl("_{vir}"),alignment=0.5,color=fsc_color('yellow')
        fsc_plot,xrange,[virTempRange[0],virTempRange[0]],line=0,color=fsc_color('orange'),/overplot
        fsc_plot,xrange,[virTempRange[1],virTempRange[1]],line=0,color=fsc_color('orange'),/overplot
      endif
      
  end_PS
  
  start_PS, sP.plotPath + plotBase + '_rad_gmem.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    loadct, 33, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
    ;cgLoadCT,13,/brewer,bottom=1
    ;loadColorTable, 'helix'
    
    tvim,h2rt_gmem^exp,pcharsize=!p.charsize,scale=0,clip=-1,$;,/c_map
         xtitle="r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}"),ytitle=ytitle,$
         title=str(sP.res)+textoidl("^3")+" "+titleName+" (z="+string(redshift,format='(f3.1)')+") gmem",$
         barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=xrange,yrange=yrange,$
         /xlog,xticks=6,xtickv=[0.01,0.05,0.1,0.2,0.5,1.0,3.0],$
         xtickname=['0.01','0.05','0.1','0.2','0.5','1','3'],xmargin=2.0

     barvals = findgen(ndivs+1)/ndivs*(max(h2rt_gmem^exp)-min(h2rt_gmem^exp)) + min(h2rt_gmem^exp)
     ticknames = textoidl('10^{'+str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)'))+'}')
     ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
     colorbar,bottom=1,range=minmax(h2rt_gmem),position=[0.83,0.155,0.87,0.925],$
       /vertical,/right,title=textoidl("M_{sun,tot}"),divisions=ndivs,ticknames=ticknames,ncolors=255
         
      if keyword_set(virTempRange) then begin
        fsc_text,xrange[1]*0.45,(yrange[1]-yrange[0])*0.92+yrange[0],massBinStr,alignment=0.5,color=fsc_color('white')
        fsc_text,xrange[0]*1.6,virTempRange[1]*1.02,"T"+textoidl("_{vir}"),alignment=0.5,color=fsc_color('yellow')
        fsc_plot,xrange,[virTempRange[0],virTempRange[0]],line=0,color=fsc_color('orange'),/overplot
        fsc_plot,xrange,[virTempRange[1],virTempRange[1]],line=0,color=fsc_color('orange'),/overplot
      endif
      
  end_PS
  !except = 1
end

; plotTempRad2DHisto(): plot 2d histogram of gas temperature (log Kelvin) as a function of r_gas/r_vir
;                       as well as vertical slices

pro plotTempRad2DHisto, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; config
  ;massBins = [0.0,100.0] ; no massbins
  massBins = [9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5] ; log(M)
  
  ; select one of:
  curTemp     = 0 ; current temperature
  maxPastTemp = 1 ; maximum previous temperature
  maxTempTime = 0 ; time of maximum previous temperature
  accTime     = 0 ; redshift of accretion across virial radius of parent halo
  accTvir     = 0 ; virial temperature of parent halo at time of accretion
  
   ; if maxPastTemp or maxTempTime, choose one of: population min,mean,max to one statistic per
   ; gas cell instead of all the tracers individually
   trPopMax   = 0
   trPopMin   = 0 
   trPopMean  = 0
   
  tVirNorm    = 0  ; normalize temperature by virial temp of parent halo at the starting time
  tVirAccNorm = 1  ; normalize temperature by virial temp of parent halo at the time of accretion
   
  parNorm  = 'pri' ; pri,sec (normalize r_vir and temp (if doing tVirNorm) by primary or secondary parent)
  sgSelect = 'pri' ; pri,sec,all (subhalo selection config, only plot gas in this category of subhalos)
  
  ; consider only the subset with recorded accretion times for certain types of plots
  mergerTreeSubset    = 0
  accretionTimeSubset = 0
  
  if accTime or accTvir or tVirAccNorm then begin
    mergerTreeSubset    = 1
    accretionTimeSubset = 1
  endif
  
  ; sanity check config parameters
  if accTime or accTvir then if sP.trMCPerCell eq -1 then message,'vel tracers not implemented yet'
  if accTime or accTvir then if tVirNorm then message,'no'
  if accTime or accTvir or tVirAccNorm then if parNorm ne 'pri' or sgSelect ne 'pri' then message,'better check you mean this'
  if trPopMax or trPopMin or trPopMean then if ~maxPastTemp and ~maxTempTime then message,'Error'
  if total(curTemp+maxPastTemp+maxTempTime+accTime+accTvir) ne 1 then message,'Just one'
  if total(tVirNorm+tVirAccNorm) gt 1 then message,'Only one allowed'
  if tVirNorm or tVirAccNorm then if ~curTemp and ~maxPastTemp then message,'Normalize what?'

  ; load group catalog just for counts of objects in each mass bin
  gc = loadGroupCat(sP=sP,/skipIDs)
  subgroupMasses = codeMassToLogMsun(gc.subgroupMass)
  gc = !NULL
  
  ; get normalized r/r_vir
  gcRad = gcSubsetProp(sP=sP,select=sgSelect,/rVirNorm,parNorm=parNorm,$
                       mergerTreeSubset=mergerTreeSubset,accretionTimeSubset=accretionTimeSubset)
  stop
  ; get current or maxPast temperature, or accretion time / tvir at accretion
  gcTemp = gcSubsetProp(sP=sP,select=sgSelect,curTemp=curTemp,maxPastTemp=maxPastTemp,$
                        maxTempTime=maxTempTime,accTime=accTime,accTvir=accTvir,$
                        mergerTreeSubset=mergerTreeSubset,accretionTimeSubset=accretionTimeSubset,$
                        trPopMin=trPopMin,trPopMax=trPopMax,trPopMean=trPopMean)

  ; normalize by halo Tvir at current time
  if tVirNorm then begin
    ; calculate temperatures of parents and take (log) ratio
    gcVirTemp = gcSubsetProp(sP=sP,select=sgSelect,/virTemp,parNorm=parNorm)
    
    gcTemp.gal  = alog10(10.0^gcTemp.gal / 10.0^gcVirTemp.gal)
    gcTemp.gmem = alog10(10.0^gcTemp.gmem / 10.0^gcVirTemp.gmem)
  endif
  
  ; normalize by halo Tvir at time of accretion
  if tVirAccNorm then begin
    ; calculate temperatures of parents at time of accretion and take (log) ratio
    gcAccTvir = gcSubsetProp(sP=sP,select=sgSelect,/accTvir,$
                             mergerTreeSubset=mergerTreeSubset,accretionTimeSubset=accretionTimeSubset)
    
    gcTemp.gal  = alog10(10.0^gcTemp.gal / 10.0^gcAccTvir.gal)
    gcTemp.gmem = alog10(10.0^gcTemp.gmem / 10.0^gcAccTvir.gmem)
  endif
  
  ; load gas masses if necessary (sph or byGas arepo)
  if sP.trMCPerCell eq 0 then $
    gcMass = gcSubsetProp(sP=sP,select=sgSelect,/curSingleVal,singleValField='mass',$
                          mergerTreeSubset=mergerTreeSubset,accretionTimeSubset=accretionTimeSubset)
  
  ; calculate masses of parents (for mass binning halos only)
  parentMass = gcSubsetProp(sP=sP,select=sgSelect,/parMass,$
                            mergerTreeSubset=mergerTreeSubset,accretionTimeSubset=accretionTimeSubset)

  for j=0,n_elements(massBins)-2 do begin
  
    ; plot config
    xrange = alog10([0.01,3.0])
    yrange = [4.0,7.0]
  
    binSizeRad  = 0.04 / (sP.res/128) ;0.04
    binSizeTemp = 0.05 / (sP.res/128) ;0.04
    
    ; preserve number of bins in log(rad) histogram and if changing yrange
    nBinsRad_linear = ceil((10.0^xrange[1]-10.0^xrange[0])/binSizeRad)+1
    nBinsTemp = ((yrange[1]-yrange[0])/binSizeTemp)+1
    binSizeRad_log  = (xrange[1]-xrange[0])/(nBinsRad_linear-1)
    
    if tVirNorm or tVirAccNorm then begin
      yrange = [-2.5,1.0]
      binSizeTemp = (yrange[1]-yrange[0])/(nBinsTemp-1)
    endif
    
    if accTime then begin
      yrange = [sP.redshift,4.0]
      binSizeTemp = (yrange[1]-yrange[0])/(nBinsTemp-1)
    endif
    
    ; select members of this parent mass bins and r>0<inf
    wGal  = where(parentMass.gal gt massBins[j] and parentMass.gal le massBins[j+1] and $
                  gcRad.gal gt 0.0 and finite(gcRad.gal),count1)
    wGmem = where(parentMass.gmem gt massBins[j] and parentMass.gmem le massBins[j+1] and $
                  gcRad.gmem gt 0.0 and finite(gcRad.gmem),count2)
    
    ; count for output
    wGCMassBin = where(subgroupMasses gt massBins[j] and subgroupMasses le massBins[j+1],count_gc)
    
    count_mt = 0
    if n_elements(mt) gt 0 then $
      wMTMassBin = where(hMass gt massBins[j] and hMass le massBins[j+1],count_mt)
    
    print,j,count1,count2,count_gc,count_mt
    
    ; select all particles/tracers in this mass bin
    temp_gal  = gcTemp.gal[wGal]
    temp_gmem = gcTemp.gmem[wGmem]
    
    rad_gal   = alog10( gcRad.gal[wGal] )
    rad_gmem  = alog10( gcRad.gmem[wGmem] )

    if count1 eq 0 or count2 eq 0 then continue ; no halos in this mass bin
    
    ; do mass weighting
    if sP.trMCPerCell eq 0 then begin
      ; plotting by gas cell, load gas mass subsets
      mass_gal  = gcMass.gal[wGal] * units.UnitMass_in_Msun
      mass_gmem = gcMass.gmem[wGmem] * units.UnitMass_in_Msun
      
      h2rt_gal = hist_nd_weight(transpose([[rad_gal],[temp_gal]]),weight=mass_gal,$
        [binSizeRad_log,binSizeTemp],$
        min=[xrange[0]-binSizeRad_log*0.5,yrange[0]-binSizeTemp*0.5],$
        max=[xrange[1]+binSizeRad_log*0.49,yrange[1]+binSizeTemp*0.49])
      h2rt_gmem = hist_nd_weight(transpose([[rad_gmem],[temp_gmem]]),weight=mass_gmem,$
        [binSizeRad_log,binSizeTemp],$
        min=[xrange[0]-binSizeRad_log*0.5,yrange[0]-binSizeTemp*0.5],$
        max=[xrange[1]+binSizeRad_log*0.49,yrange[1]+binSizeTemp*0.49])
    endif else begin
      ; create unweighted 2d histo (separate for galaxy/group member and composite)
      h2rt_gal  = hist_nd(transpose([[rad_gal],[temp_gal]]),[binSizeRad_log,binSizeTemp],$
        min=[xrange[0]-binSizeRad_log*0.5,yrange[0]-binSizeTemp*0.5],$
        max=[xrange[1]+binSizeRad_log*0.49,yrange[1]+binSizeTemp*0.49])
      h2rt_gmem = hist_nd(transpose([[rad_gmem],[temp_gmem]]),[binSizeRad_log,binSizeTemp],$
        min=[xrange[0]-binSizeRad_log*0.5,yrange[0]-binSizeTemp*0.5],$
        max=[xrange[1]+binSizeRad_log*0.49,yrange[1]+binSizeTemp*0.49])
                
      ; plotting all tracers, multiply number histogram by constant tracer mass
      if sP.trMCPerCell gt 0 then begin
        h2rt_gal  *= sP.trMassConst * units.UnitMass_in_Msun
        h2rt_gmem *= sP.trMassConst * units.UnitMass_in_Msun
      endif
      if sP.trMCPerCell eq -1 then begin
        h2rt_gal  *= sP.targetGasMass * units.UnitMass_in_Msun
        h2rt_gmem *= sP.targetGasMass * units.UnitMass_in_Msun
      endif
    endelse

    ; 2d histo plot config
    if curTemp then begin
      plotBase = "tcur_"+sgSelect+'_'+parNorm
      ytitle   = "log ( T"+textoidl("_{cur}")+" )"
    endif
    if maxPastTemp then begin
      plotBase = "tmax_"+sgSelect+'_'+parNorm
      ytitle   = "log ( T"+textoidl("_{max}")+" )"
    endif
    
    if maxTempTime then begin
      plotBase = "tmaxtime_"+sgSelect+'_'+parNorm
      ytitle   = "Maximum Temperature Redshift"
    endif
    if accTime then begin
      plotBase = "acctime_"+sgSelect+'_'+parNorm
      ytitle   = "Accretion Redshift"
    endif
    if accTvir then begin
      plotBase = "acctvir_"+sgSelect+'_'+parNorm
      ytitle   = "Parent "+textoidl("T_{vir}")+" at "+textoidl("t_{acc}")
    endif
    
    if tVirNorm then begin
      plotBase = strmid(plotBase,0,4)+"_tvirNorm_"+sgSelect+'_'+parNorm
      ytitle   = strmid(ytitle,0,strlen(ytitle)-2) + " / T"+textoidl("_{vir,cur}"+" )")
    endif
    if tVirAccNorm then begin
      plotBase = strmid(plotBase,0,4)+"_tvirAccNorm_"+sgSelect+'_'+parNorm
      ytitle   = strmid(ytitle,0,strlen(ytitle)-2) + " / T"+textoidl("_{vir,acc}"+" )")
    endif
    
    ; extra config for mass bins
    massBinStr   = !NULL
    virTempRange = !NULL
    
    if n_elements(massBins) gt 2 then begin
      plotBase = plotBase+"_mbin="+str(j)
      
      massBinStr = string(massBins[j],format='(f4.1)') + '-' + string(massBins[j+1],format='(f4.1)')
      
      massRangeCode = 10.0^[massBins[j],massBins[j+1]] / 1e10
      virTempRange = alog10( codeMassToVirTemp(massRangeCode,sP=sP) )
    endif
    
    ; plot 2d histo
    plot2DRadHistos,plotBase,sP,h2rt_gal,h2rt_gmem,10.0^xrange,yrange,ytitle,$
                    virTempRange=virTempRange,massBinStr=massBinStr    
    
    ; vertical slices plot config
    plotName = sP.plotPath + plotBase + '_slices.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    if n_elements(massBins) gt 2 then $
      plotName = strmid(plotName,0,strlen(plotName)-4) + '.mbin='+str(j)+'.eps'
    
    xtitle = ytitle
    
    ; plot vertical slices
    plotVerticalSlices,rad_gal,rad_gmem,temp_gal,temp_gmem,plotName,binSizeTemp*2.0,yrange,xtitle
    
  endfor ;j
  
end
