; plotVsRedshift.pro
; feedback - plots skipping tconst/tvircur/tviracc definitions in favor of redshift panels
; dnelson jun.2013

; plotRatesFracsInRedshift():

pro plotRatesFracsInRedshift
  message,'unchecked with new format'
  ; config
  runs       = ['feedback','tracer','gadget']
  redshifts  = [3.0,2.0,1.0,0.0]
  res        = 256
  accMode    = 'smooth' ; accretion mode: all, smooth, clumpy, stripped, recycled
  timeWindow = 1000.0   ; consider accretion over this past time range (Myr)
                        ; 250.0 500.0 1000.0 "all" "tVir_tIGM" or "tVir_tIGM_bin"
  tVirInd    = 0        ; use Tmax/Tviracc=1 to separate hot vs. cold
  massBinInd = 4        ; plot 11.0<logM<11.5 halos for Tmax histos
  
  accModes   = ['all','smooth','clumpy','stripped','recycled'] ; for fractional plot
  
  ; plot config
  lines   = [0,1]       ; gal/both,gmem
  linesAM = [0,1,2,3,4] ; for each of accModes
  sK      = 3           ; smoothing kernel size  
  cInd    = 1           ; color index

  xrange_halo = [10.0,12.0]
  yrange_frac = [0.0,1.0]
  yrange_rate = [-5,-1.0]
  
  xrange_tmax = [-2.0,1.5]
  yrange_hist = [6e-4,2e-1]
  
  ; load
  foreach run,runs,i do begin
    sP_z = {} & mbv_z = {} & bth_z = {} & mode_z = {}
    
    ; make for all the redshifts
    foreach redshift,redshifts,j do begin
      sP_z  = mod_struct(sP_z, 'redshift'+str(j), simParams(res=res,run=run,redshift=redshift))
      mbv_z = mod_struct(mbv_z, 'redshift'+str(j), $
        haloMassBinValues(sP=sP_z.(j),accMode=accMode,timeWindow=timeWindow))
      bth_z = mod_struct(bth_z, 'redshift'+str(j), $
        binValMaxHistos(sP=sP_z.(j),accMode=accMode,timeWindow=timeWindow))
        
      ; for each redshift, make for all the accretion modes
      mbv_mode = {}
      
      foreach amCur,accModes,k do begin
        if amCur eq 'recycled' and sP_z.(j).gfmWinds eq 0 then continue ; skip recycled for nonWinds  
        mbv_mode = mod_struct(mbv_mode, 'mode'+str(k), $
          haloMassBinValues(sP=sP_z.(j),accMode=amCur,timeWindow=timeWindow))
      endforeach
      
      mode_z = mod_struct(mode_z, 'redshift'+str(j), mbv_mode)
    endforeach
    
    ; put this mode collection into mbv, once per run, and sP collection also
    amv = mod_struct(amv, 'amv'+str(i), mode_z)
    mbv = mod_struct(mbv, 'mbv'+str(i), mbv_z)
    bth = mod_struct(bth, 'bth'+str(i), bth_z)
    sP  = mod_struct(sP, 'sP'+str(i), sP_z)
  endforeach

  ; strings
  if isnumeric(timeWindow) then twStr = str(fix(timeWindow))
  if ~isnumeric(timeWindow) then twStr = "-"+timeWindow
  
  plotStr   = ''
  simNames  = []
  simColors = []
  
  foreach run,runs,i do begin
    plotStr   = plotStr + sP.(i).(0).plotPrefix + '.'
    simNames  = [simNames, sP.(i).(0).simName]
    simColors = [simColors, sP.(i).(0).colors[cInd]]
  endforeach

  plotStr += str(res) + '_' + str(sP.(0).(0).snap) + '_tw' + twStr + '_am-' + accMode
  
  pos = plot_pos(total=n_elements(redshifts),/gap) ; plot positioning (3x2, 2x2, or 1x2 with gaps)
  
  ; cold fraction (allgal,gmem)
  start_PS, sP.(0).(0).plotPath + 'coldFracRedshift.allgal-gmem.' + plotStr + '.eps', /big
    
    for zind=0,3 do begin
      cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_frac,/xs,/ys,$
        ytitle="Cold Fraction",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
        /noerase,pos=pos[zind]
    
      for i=0,n_tags(sP)-1 do begin
        xvals = mbv.(i).(zind).logMassBinCen
        yvals = smooth(mbv.(i).(zind).fracMedian.allGal.tVirAcc[tVirInd,*],sK,/nan)
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[0],/overplot ; allgal
        
        yvals = smooth(mbv.(i).(zind).fracMedian.gmem.tVirAcc[tVirInd,*],sK,/nan)
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[1],/overplot ; gmem
      endfor
      
      ; redshift legend
      legend,textoidl('z_{ }=_{ }'+string(redshifts[zind],format='(f3.1)')),/top,/right,$
        charsize=!p.charsize-0.2
      
      ; tvir axis on top
      ;tvirrange = alog10(codeMassToVirTemp( logMsunToCodeMass(xrange),redshift=redshifts[zind] ))
      ;cgAxis,/xaxis,xrange=tvirrange,/xs,color=cgColor('green')
      ;cgText,mean([x0,x1]),y2+0.06,textoidl("T_{vir} [_{ }log K_{ }] (z="+str(fix(redshifts[zind]))+")"),$
      ;  alignment=0.5,/normal
      
      ; panel-specific legends
      if zind eq 2 then $
        legend,['galaxy','halo'],linestyle=lines,$
        color=cgColor('dark gray'),textcolors=['dark gray','dark gray'],$
        linesize=0.4,box=0,/top,/left,charsize=!p.charsize-0.2
        
      if zind eq 3 then $
        legend,simNames,textcolors=simColors,box=0,/top,/left,charsize=!p.charsize-0.2
    endfor

  end_PS  
  
  ; accretion rates (allgal)
  start_PS, sP.(0).(0).plotPath + 'accRateRedshift.allgal.' + plotStr + '.eps', /big
    
    for zind=0,3 do begin
      cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_rate,/xs,/ys,$
        ytitle=textoidl("dM_{gas}/dt  M_{halo}^{-1} [_{ }Gyr^{-1 }]"),$
        xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
        /noerase,pos=pos[zind]
    
      for i=0,n_tags(sP)-1 do begin
        xvals = mbv.(i).(zind).logMassBinCen
        yvals = smooth(alog10(mbv.(i).(zind).hotMedian.allGal.tVirAcc[tVirInd,*] / $
                       10.0^mbv.(i).(zind).logMassBinCen * 1e9),sK,/nan) ; yr->Gyr
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[0],/overplot ; allgal, hot
        
        yvals = smooth(alog10(mbv.(i).(zind).coldMedian.allGal.tVirAcc[tVirInd,*] / $
                       10.0^mbv.(i).(zind).logMassBinCen * 1e9),sK,/nan) ; yr->Gyr
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[1],/overplot ; allgal, cold
      endfor
      
      ; redshift legend
      legend,textoidl('z_{ }=_{ }'+string(redshifts[zind],format='(f3.1)')),/bottom,/left,$
        charsize=!p.charsize-0.2
      
      ; panel-specific legends
      if zind eq 2 then $
        legend,['hot mode (galaxy)','cold mode (galaxy)'],linestyle=lines,color=cgColor('dark gray'),$
        textcolors=replicate('dark gray',2),box=0,/top,/left,charsize=!p.charsize-0.2
        
      if zind eq 3 then $
        legend,simNames,textcolors=simColors,box=0,/top,/left,charsize=!p.charsize-0.2
    endfor

  end_PS
  
  ; accretion rates (gmem)
  start_PS, sP.(0).(0).plotPath + 'accRateRedshift.gmem.' + plotStr + '.eps', /big
    
    for zind=0,3 do begin
      cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_rate,/xs,/ys,$
        ytitle=textoidl("dM_{gas}/dt  M_{halo}^{-1} [_{ }Gyr^{-1 }]"),$
        xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
        /noerase,pos=pos[zind]
    
      for i=0,n_tags(sP)-1 do begin
        xvals = mbv.(i).(zind).logMassBinCen
        yvals = smooth(alog10(mbv.(i).(zind).hotMedian.gmem.tVirAcc[tVirInd,*] / $
                       10.0^mbv.(i).(zind).logMassBinCen * 1e9),sK,/nan) ; yr->Gyr
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[0],/overplot ; gmem, hot
        
        yvals = smooth(alog10(mbv.(i).(zind).coldMedian.gmem.tVirAcc[tVirInd,*] / $
                       10.0^mbv.(i).(zind).logMassBinCen * 1e9),sK,/nan) ; yr->Gyr
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[1],/overplot ; gmem, cold
      endfor
      
      ; redshift legend
      legend,textoidl('z_{ }=_{ }'+string(redshifts[zind],format='(f3.1)')),/bottom,/left,$
        charsize=!p.charsize-0.2
      
      ; panel-specific legends
      if zind eq 2 then $
        legend,['hot mode (halo)','cold mode (halo)'],linestyle=lines,color=cgColor('dark gray'),$
        textcolors=replicate('dark gray',2),box=0,/top,/left,charsize=!p.charsize-0.2
        
      if zind eq 3 then $
        legend,simNames,textcolors=simColors,box=0,/top,/left,charsize=!p.charsize-0.2
    endfor

  end_PS
  
  ; maximum temp histos (allgal,gmem) (one halo mass bin)
  start_PS, sP.(0).(0).plotPath + 'maxTempHistosRedshift.massBin' + $
            str(massBinInd) + '.' + plotStr + '.eps', /big
    
    for zind=0,3 do begin
      cgPlot,[0],[0],/nodata,xrange=xrange_tmax,yrange=yrange_hist,/xs,/ys,/ylog,yminor=0,$
        ytitle=textoidl("Gas Mass Fraction"),$
        xtitle=textoidl("log ( T_{max} / T_{vir,acc} )"),$
        /noerase,pos=pos[zind],xticks=4,xtickv=[-2,-1,0,1,1.5]
    
      for i=0,n_tags(sP)-1 do begin
      
        ; gal
        xvals = bth.(i).(zind).binLocRatio
        yvals = float( bth.(i).(zind).allGal.TmaxTviracc[massBinInd,*] ) / $
                total( bth.(i).(zind).allGal.TmaxTviracc[j,*] )
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[0],/overplot ; allgal
        
        ; gmem
        yvals = float( bth.(i).(zind).gmem.TmaxTviracc[massBinInd,*] ) / $
                total( bth.(i).(zind).gmem.TmaxTviracc[j,*] )
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[1],/overplot ; gmem
      endfor
      
      ; redshift legend
      legend,textoidl('z_{ }=_{ }'+string(redshifts[zind],format='(f3.1)')),/bottom,/left,$
        charsize=!p.charsize-0.2
      
      cgPlot,[0,0],[9e-4,0.13],line=2,color=cgColor('black'),/overplot
      
      ; panel-specific legends
      if zind eq 2 then $
        legend,['galaxy','halo'],linestyle=lines,color=cgColor('dark gray'),$
        textcolors=replicate('dark gray',2),box=0,/top,/left,charsize=!p.charsize-0.2
        
      if zind eq 3 then $
        legend,simNames,textcolors=simColors,box=0,/top,/left,charsize=!p.charsize-0.2
    endfor

  end_PS
  
  ; fraction of accretion by mode (allgal)
  pos_local = plot_pos(total=2*n_elements(redshifts),/gap)
  if accModes[0] ne 'all' then message,'Error: Not going to work.'
  
  start_PS, sP.(0).(0).plotPath + 'accRateFracsRedshift.allgal.' + plotStr + '.eps', xs=9, ys=12
    
    for zind=0,3 do begin
      cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_frac,/xs,/ys,$
        ytitle=textoidl("Fraction"),$
        xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
        /noerase,pos=pos_local[2*zind]
    
      ; hot (left column)
      for i=0,n_tags(sP)-1 do begin
        for j=1,n_tags(amv.(i).(zind))-1 do begin ; take ratio to first
        
          xvals = amv.(i).(zind).(0).logMassBinCen
          yvals = smooth(amv.(i).(zind).(j).hotMedian.allGal.tVirAcc[tVirInd,*] / $
                         amv.(i).(zind).(0).hotMedian.allGal.tVirAcc[tVirInd,*],sK,/nan)
          cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=linesAM[j],/overplot
          
        endfor
      endfor
      
      ; redshift legend
      legend,textoidl('z_{ }=_{ }'+string(redshifts[zind],format='(f3.1)')),/top,/right,$
        charsize=!p.charsize-0.2
      
      cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_frac,/xs,/ys,$
        ytitle=textoidl("Fraction"),$
        xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
        /noerase,pos=pos_local[2*zind+1]
      
      ; cold (right column)
      for i=0,n_tags(sP)-1 do begin
        for j=1,n_tags(amv.(i).(zind))-1 do begin ; take ratio to first
        
          xvals = amv.(i).(zind).(0).logMassBinCen
          yvals = smooth(amv.(i).(zind).(j).coldMedian.allGal.tVirAcc[tVirInd,*] / $
                         amv.(i).(zind).(0).coldMedian.allGal.tVirAcc[tVirInd,*],sK,/nan)
          cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=linesAM[j],/overplot
          
        endfor
      endfor
      
      ; redshift legend
      legend,textoidl('z_{ }=_{ }'+string(redshifts[zind],format='(f3.1)')),/top,/right,$
        charsize=!p.charsize-0.2
      
      ; panel-specific legends
      if zind eq 2 then $
        legend,accModes[1:*],linestyle=linesAM[1:*],box=0,linesize=0.4,/top,/left,charsize=!p.charsize-0.2
      if zind eq 3 then $
        legend,simNames,textcolors=simColors,box=0,/top,/left,charsize=!p.charsize-0.2
    endfor
    
    ; hot/cold labels
    cgText,mean( (pos_local[0])[0:2:2] ), (pos_local[0])[3]+0.02, "Hot Mode",alignment=0.5,/normal
    cgText,mean( (pos_local[1])[0:2:2] ), (pos_local[1])[3]+0.02, "Cold Mode",alignment=0.5,/normal

  end_PS
  stop
end

; binAccRateVsRedshift():

function binAccRateVsRedshift, sP=sP, accMode=accMode, binWindow=BW
  message,'needs update to new format'
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  if ~sP.gfmWinds and accMode eq 'recycled' then message,'Error: Request recycled on non-winds run.'
  if ~keyword_set(BW) then message,'time window required (in Myr)'
  
  if sP.gfmWinds eq 1 and (accMode eq 'smooth' or accMode eq 'stripped' or accMode eq 'clumpy') then begin
    print,'Switching [' + accMode + '] to [' + accMode + '_rec] for GFM run.'
    accMode = accMode + '_rec'
  endif
  
  ; current time, earliest time, number of time steps
  h = loadSnapshotHeader(sP=sP)
  curtime = 1/h.time - 1 ; redshift
  curtime = redshiftToAgeFlat(curtime)*1e9 ; yr
  
  snapTimes = redshiftToAgeFlat(snapNumToRedshift(/all,sP=sP))*1e9 ; yr
  
  startTime = snapNumToAgeFlat(snap=sP.groupCatRange[0],sP=sP)*1e9 ;yr
  totalTime = curtime - startTime
  
  timeWindow = BW * 1e6 ; convert input Myr to yr
  
  nTimeBins = ceil(totalTime / timeWindow)
  
  ; config
  maxTempTvirFac = 1.0 ; differentiate hot and cold
  
  if sP.trMCPerCell le 0 then massPerPart = sP.targetGasMass ; SPH or vel tracer
  if sP.trMCPerCell gt 0 then massPerPart = sP.trMassConst ; MC tracer
  
  ; check if save exists
  saveFilename = sP.derivPath + 'binnedVals/binRateZ.' + sP.saveTag + '.' + sP.savPrefix + str(sP.res) + $
    '.' + str(sP.snap) + '.' + accMode + '_bw' + str(BW) + '.sav'
  
  ; results exist, return
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif   
  
  ; make a uniform gas selection at the start
  at = accretionTimes(sP=sP)
  mt = mergerTreeSubset(sP=sP)
  
  wAm = accModeInds(at=at,accMode=accMode,sP=sP,/mask)
    
  ; reverse histogram parent IDs of all particles/tracers in this selection
  gcIndOrig = mergerTreeRepParentIDs(mt=mt,sP=sP,/compactMtS)
  
  maxHist = max([gcIndOrig.gal[wAm.gal],gcIndOrig.gmem[wAm.gmem],gcIndOrig.stars[wAm.stars]])
  hist_gal   = histogram(gcIndOrig.gal[wAm.gal],min=0,max=maxHist,loc=loc_gal,rev=rev_gal)
  hist_gmem  = histogram(gcIndOrig.gmem[wAm.gmem],min=0,max=maxHist,loc=loc_gmem,rev=rev_gmem)
  hist_stars = histogram(gcIndOrig.stars[wAm.stars],min=0,max=maxHist,loc=loc_stars,rev=rev_stars)
  gcIndOrig = !NULL
  
  ; load max temps, current tvir, tvir at accretion
  accTvir = gcSubsetProp(sP=sP,select='pri',/accTvir,/accretionTimeSubset,accMode=accMode)
  maxTemp = gcSubsetProp(sP=sP,select='pri',/maxPastTemp,/accretionTimeSubset,accMode=accMode)

  ; load group cat for subgroup masses
  gc = loadGroupCat(sP=sP,/skipIDs)
  gcMasses = codeMassToLogMsun(gc.subgroupMass[mt.galcatIDList])
  gc = !NULL
  
  ; structure to store results
  accRate = { gal_cold   : fltarr(nTimeBins,n_elements(mt.galcatIDList))+!values.f_nan  ,$
              gmem_cold  : fltarr(nTimeBins,n_elements(mt.galcatIDList))+!values.f_nan  ,$
              stars_cold : fltarr(nTimeBins,n_elements(mt.galcatIDList))+!values.f_nan  ,$
              both_cold  : fltarr(nTimeBins,n_elements(mt.galcatIDList))+!values.f_nan  ,$
              gal_hot    : fltarr(nTimeBins,n_elements(mt.galcatIDList))+!values.f_nan  ,$
              gmem_hot   : fltarr(nTimeBins,n_elements(mt.galcatIDList))+!values.f_nan  ,$
              stars_hot  : fltarr(nTimeBins,n_elements(mt.galcatIDList))+!values.f_nan  ,$
              both_hot   : fltarr(nTimeBins,n_elements(mt.galcatIDList))+!values.f_nan   }
                  
  ; loop over all tracked subgroups
  for i=0L,maxHist do begin
   
    if hist_gal[i] gt 0 then begin
      ; list of indices of galaxy gas particles in this subgroup
      loc_inds_gal = rev_gal[rev_gal[i]:rev_gal[i+1]-1]
      
      ; galaxy accretion defined as (rho,temp) joining time or 0.15rvir crossing time
      loc_atime_gal = reform(at.accTimeRT_gal[wAm.gal[loc_inds_gal]])
      
      r_crossing_time = reform(at.accTime_gal[sP.radIndGalAcc,wAm.gal[loc_inds_gal]])
      w = where(r_crossing_time gt loc_atime_gal,count)
      if count gt 0 then loc_atime_gal[w] = r_crossing_time[w]
      
      loc_atime_gal = 1/loc_atime_gal - 1 ; redshift
      loc_atime_gal = redshiftToAgeFlat(loc_atime_gal)*1e9 ; yr
      
      ; split into hot and cold
      loc_maxt_gal    = maxTemp.gal[loc_inds_gal]
      loc_accTvir_gal = accTvir.gal[loc_inds_gal]
      
      w_cold = where(10.0^loc_maxt_gal / 10.0^loc_accTvir_gal le maxTempTvirFac,$
                     count_cold,ncomp=count_hot,comp=w_hot)
      
      ; histogram in time
      if count_cold gt 0 then begin
        hh_cold = histogram(loc_atime_gal[w_cold],min=startTime,max=curtime,binsize=timeWindow,loc=locc)
        accRate.gal_cold[*,i] = hh_cold
      endif
      
      if count_hot gt 0 then begin
        hh_hot = histogram(loc_atime_gal[w_hot],min=startTime,max=curtime,binsize=timeWindow,loc=loch)
        accRate.gal_hot[*,i]  = hh_hot
      endif
      
    endif ; hist_gal[i]

    if hist_gmem[i] gt 0 then begin
      ; list of indices of group member gas particles in this subgroup
      loc_inds_gmem = rev_gmem[rev_gmem[i]:rev_gmem[i+1]-1]
      
      ; corresponding accretion times for these particles
      loc_atime_gmem = reform(at.accTime_gmem[sP.radIndHaloAcc,wAm.gmem[loc_inds_gmem]])
      loc_atime_gmem = 1/loc_atime_gmem - 1 ; redshift
      loc_atime_gmem = redshiftToAgeFlat(loc_atime_gmem)*1e9 ; yr
      
      ; split into hot and cold
      loc_maxt_gmem    = maxTemp.gmem[loc_inds_gmem]
      loc_accTvir_gmem = accTvir.gmem[loc_inds_gmem]
      
      w_cold = where(10.0^loc_maxt_gmem / 10.0^loc_accTvir_gmem le maxTempTvirFac,$
                     count_cold,ncomp=count_hot,comp=w_hot)
      
      ; histogram in time
      if count_cold gt 0 then begin
        hh_cold = histogram(loc_atime_gmem[w_cold],min=startTime,max=curtime,binsize=timeWindow,loc=locc)
        accRate.gmem_cold[*,i] = hh_cold
      endif
      
      if count_hot gt 0 then begin
        hh_hot = histogram(loc_atime_gmem[w_hot],min=startTime,max=curtime,binsize=timeWindow,loc=loch)
        accRate.gmem_hot[*,i]  = hh_hot
      endif
      
    endif ; hist_gmem[i]

    if hist_stars[i] gt 0 then begin
      ; list of indices of star gas particles in this subgroup
      loc_inds_stars = rev_stars[rev_stars[i]:rev_stars[i+1]-1]
      
      ; corresponding accretion times for these particles:
      ; stellar accretion defined as (rho,temp) joining time or 0.15rvir crossing time (most recent)
      loc_atime_stars = reform(at.accTimeRT_stars[wAm.stars[loc_inds_stars]])
      
      r_crossing_time = reform(at.accTime_stars[sP.radIndGalAcc,wAm.stars[loc_inds_stars]])
      w = where(r_crossing_time gt loc_atime_stars,count)
      if count gt 0 then loc_atime_stars[w] = r_crossing_time[w]
      
      ; convert from scale factor to age of the universe
      loc_atime_stars = 1/loc_atime_stars - 1 ; redshift
      loc_atime_stars = redshiftToAgeFlat(loc_atime_stars)*1e9 ; yr
      
      ; split into hot and cold
      loc_maxt_stars    = maxTemp.stars[loc_inds_stars]
      loc_accTvir_stars = accTvir.stars[loc_inds_stars]
      
      w_cold = where(10.0^loc_maxt_stars / 10.0^loc_accTvir_stars le maxTempTvirFac,$
                     count_cold,ncomp=count_hot,comp=w_hot)
      
      ; histogram in time
      if count_cold gt 0 then begin
        hh_cold = histogram(loc_atime_stars[w_cold],min=startTime,max=curtime,binsize=timeWindow,loc=locc)
        accRate.stars_cold[*,i] = hh_cold
      endif
      
      if count_hot gt 0 then begin
        hh_hot = histogram(loc_atime_stars[w_hot],min=startTime,max=curtime,binsize=timeWindow,loc=loch)
        accRate.stars_hot[*,i]  = hh_hot
      endif
      
    endif ; hist_stars[i]
    
    ; convert accretion total(counts) to msun/year
    accRate.gal_hot[*,i]     *= massPerPart * units.UnitMass_in_Msun / timeWindow
    accRate.gmem_hot[*,i]    *= massPerPart * units.UnitMass_in_Msun / timeWindow
    accRate.stars_hot[*,i]   *= massPerPart * units.UnitMass_in_Msun / timeWindow
    
    accRate.gal_cold[*,i]     *= massPerPart * units.UnitMass_in_Msun / timeWindow
    accRate.gmem_cold[*,i]    *= massPerPart * units.UnitMass_in_Msun / timeWindow
    accRate.stars_cold[*,i]   *= massPerPart * units.UnitMass_in_Msun / timeWindow
    
    accRate.both_hot[*,i]  = accRate.gal_hot[*,i] + accRate.stars_hot[*,i]
    accRate.both_cold[*,i] = accRate.gal_cold[*,i] + accRate.stars_cold[*,i]
    
  endfor ; i

  ; for a number of halo mass bins, calculate median acc rates as a function of redshift
  timeBinCen = ( locc + timeWindow*0.5 ) / 1e9 ; centers of time bins for plotting, tage(Gyr)
  logMassBins = list([8.0,13.0],[9.0,9.5],[9.5,10.0],[10.0,10.5],[10.5,11.0],[11.0,11.5],[11.5,12.5])
  
  medianRate = { gal_cold   : fltarr(nTimeBins,n_elements(logMassBins))  ,$
                 gmem_cold  : fltarr(nTimeBins,n_elements(logMassBins))  ,$
                 stars_cold : fltarr(nTimeBins,n_elements(logMassBins))  ,$
                 both_cold  : fltarr(nTimeBins,n_elements(logMassBins))  ,$
                 gal_hot    : fltarr(nTimeBins,n_elements(logMassBins))  ,$
                 gmem_hot   : fltarr(nTimeBins,n_elements(logMassBins))  ,$
                 stars_hot  : fltarr(nTimeBins,n_elements(logMassBins))  ,$
                 both_hot   : fltarr(nTimeBins,n_elements(logMassBins))   }
  
  foreach logMassBin,logMassBins,k do begin
    ; select halos in this mass bin
    wHalo = where(gcMasses ge logMassBin[0] and gcMasses lt logMassBin[1],countHalo)
    if countHalo eq 0 then continue
  
    ; loop over each time bin
    for i=0,nTimeBins-1 do begin
      ; median accretion rates of all halos in this mass bin at this time
      medianRate.gal_cold[i,k]   = median(accRate.gal_cold[i,wHalo])
      medianRate.gmem_cold[i,k]  = median(accRate.gmem_cold[i,wHalo])
      medianRate.stars_cold[i,k] = median(accRate.stars_cold[i,wHalo])
      medianRate.both_cold[i,k]  = median(accRate.both_cold[i,wHalo])
      
      medianRate.gal_hot[i,k]   = median(accRate.gal_hot[i,wHalo])
      medianRate.gmem_hot[i,k]  = median(accRate.gmem_hot[i,wHalo])
      medianRate.stars_hot[i,k] = median(accRate.stars_hot[i,wHalo])
      medianRate.both_hot[i,k]  = median(accRate.both_hot[i,wHalo])
    endfor
  
  endforeach
  
  r = {accRateByHalo:accRate,accRateByZ:medianRate,timeRange:[startTime,curtime],$
       logMassBins:logMassBins,timeBinCen:timeBinCen,gcMasses:gcMasses}

  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  return, r  

end

; plotAccRateVsRedshift():

pro plotAccRateVsRedshift

  ; config
  sP = simParams(res=256,run='tracer',redshift=2.0)
  accMode = 'all'
  binWindow = 200.0 ; Myr
  
  massBinInd = 0 ; 0=all
  
  xrange_time = [0.0,4.0]
  yrange_rate = [0.001,1.0]
  
  colors = { cold:'slate blue', hot:'red' }
  
  ; get binned values
  bar = binAccRateVsRedshift(sP=sP,accMode=accMode,binWindow=binWindow)
  
  ; plot (1)
  plotStr = str(sP.res) + '_' + str(sP.snap) + '_bw' + str(fix(binWindow)) + '_am-' + accMode
  
  start_PS, sP.plotPath + 'accRateVsRedshift.' + plotStr + '.eps', /big
  
    cgPlot,[0],[0],/nodata,xrange=xrange_time,yrange=yrange_rate,xs=9,/ys,/ylog,yminor=0,$
      xtitle=textoidl("t_{age} [Gyr]"),ytitle="Accretion Rate"
      
    cgPlot,bar.timeBinCen,bar.accRateByZ.gal_cold[*,massBinInd],color=colors.cold,line=0,/overplot
    cgPlot,bar.timeBinCen,bar.accRateByZ.gal_hot[*,massBinInd],color=colors.hot,line=0,/overplot
    
    cgPlot,bar.timeBinCen,bar.accRateByZ.gmem_cold[*,massBinInd],color=colors.cold,line=1,/overplot
    cgPlot,bar.timeBinCen,bar.accRateByZ.gmem_hot[*,massBinInd],color=colors.hot,line=1,/overplot
      
    redshift_axis, xrange_time, yrange_rate, /ylog, sP=sP
  end_PS
  
  stop
  

end