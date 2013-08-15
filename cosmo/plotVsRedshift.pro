; plotVsRedshift.pro
; feedback - plots skipping tconst/tvircur/tviracc definitions in favor of redshift panels
; dnelson aug.2013

; plotRatesFracsInRedshift():

pro plotRatesFracsInRedshift

  ; config
  runs       = ['feedback','gadget']
  redshifts  = [3.0,2.0,1.0]
  res        = 512
  accMode    = 'smooth' ; accretion mode: all, smooth, clumpy, stripped, recycled
  timeWindow = 500.0    ; consider accretion over this past time range (Myr)
                        ; 250.0 500.0 1000.0 "all" "tVir_tIGM" or "tVir_tIGM_bin"
  tVirInd    = 0        ; use Tmax/Tviracc=1 to separate hot vs. cold
  massBinInd = 4        ; plot 11.0<logM<11.5 halos for Tmax histos
  
  accModes   = ['all','smooth','clumpy','stripped','recycled'] ; for fractional plot
  
  ; plot config
  lines   = [0,1]       ; gal/both,gmem
  linesAM = [0,1,2,3,4] ; for each of accModes
  sK      = 3           ; smoothing kernel size  
  cInd    = 1           ; color index

  xrange_halo = [9.0,12.0]
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
  
  ; cold fraction (galaxy,halo)
  start_PS, sP.(0).(0).plotPath + 'coldFracRedshift.galaxy-halo.' + plotStr + '.eps', /big
    
    for zind=0,3 do begin
      cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_frac,/xs,/ys,$
        ytitle="Cold Fraction",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
        /noerase,pos=pos[zind]
    
      for i=0,n_tags(sP)-1 do begin
        xvals = mbv.(i).(zind).logMassBinCen
        yvals = smooth(mbv.(i).(zind).galaxyMedian.coldFrac.total.tVirAcc[tVirInd,*],sK,/nan)
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[0],/overplot ; allgal
        
        yvals = smooth(mbv.(i).(zind).haloMedian.coldFrac.total.tVirAcc[tVirInd,*],sK,/nan)
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
  
  ; net accretion rates (galaxy)
  start_PS, sP.(0).(0).plotPath + 'netRateRedshift.galaxy.' + plotStr + '.eps', /big
    
    for zind=0,3 do begin
      cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_rate,/xs,/ys,$
        ytitle=textoidl("dM_{gas}/dt  M_{halo}^{-1} [_{ }Gyr^{-1 }]"),$
        xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
        /noerase,pos=pos[zind]
    
      for i=0,n_tags(sP)-1 do begin
        xvals = mbv.(i).(zind).logMassBinCen
        yvals = smooth(alog10(mbv.(i).(zind).galaxyMedian.netRate.total.hot.tVirAcc[tVirInd,*] / $
                       10.0^mbv.(i).(zind).logMassBinCen * 1e9),sK,/nan) ; yr->Gyr
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[0],/overplot ; galaxy, hot
        
        yvals = smooth(alog10(mbv.(i).(zind).galaxyMedian.netRate.total.cold.tVirAcc[tVirInd,*] / $
                       10.0^mbv.(i).(zind).logMassBinCen * 1e9),sK,/nan) ; yr->Gyr
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[1],/overplot ; galaxy, cold
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
  
  ; net accretion rates (halo)
  start_PS, sP.(0).(0).plotPath + 'netRateRedshift.halo.' + plotStr + '.eps', /big
    
    for zind=0,3 do begin
      cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_rate,/xs,/ys,$
        ytitle=textoidl("dM_{gas}/dt  M_{halo}^{-1} [_{ }Gyr^{-1 }]"),$
        xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
        /noerase,pos=pos[zind]
    
      for i=0,n_tags(sP)-1 do begin
        xvals = mbv.(i).(zind).logMassBinCen
        yvals = smooth(alog10(mbv.(i).(zind).haloMedian.netRate.total.hot.tVirAcc[tVirInd,*] / $
                       10.0^mbv.(i).(zind).logMassBinCen * 1e9),sK,/nan) ; yr->Gyr
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[0],/overplot ; halo, hot
        
        yvals = smooth(alog10(mbv.(i).(zind).haloMedian.netRate.total.cold.tVirAcc[tVirInd,*] / $
                       10.0^mbv.(i).(zind).logMassBinCen * 1e9),sK,/nan) ; yr->Gyr
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[1],/overplot ; halo, cold
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
        xvals = bth.(i).(zind).params.binLoc.TmaxTviracc
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
  
  ; fraction of accretion by mode (galaxy)
  pos_local = plot_pos(total=2*n_elements(redshifts),/gap)
  if accModes[0] ne 'all' then message,'Error: Not going to work.'
  
  start_PS, sP.(0).(0).plotPath + 'accRateFracsRedshift.galaxy.' + plotStr + '.eps', xs=9, ys=12
    
    for zind=0,3 do begin
      cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_frac,/xs,/ys,$
        ytitle=textoidl("Fraction"),$
        xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
        /noerase,pos=pos_local[2*zind]
    
      ; hot (left column)
      for i=0,n_tags(sP)-1 do begin
        for j=1,n_tags(amv.(i).(zind))-1 do begin ; take ratio to first
        
          xvals = amv.(i).(zind).(0).logMassBinCen
          yvals = smooth(amv.(i).(zind).(j).galaxyMedian.netRate.total.hot.tVirAcc[tVirInd,*] / $
                         amv.(i).(zind).(0).galaxyMedian.netRate.total.hot.tVirAcc[tVirInd,*],sK,/nan)
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
          yvals = smooth(amv.(i).(zind).(j).galaxyMedian.netRate.total.cold.tVirAcc[tVirInd,*] / $
                         amv.(i).(zind).(0).galaxyMedian.netRate.total.cold.tVirAcc[tVirInd,*],sK,/nan)
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

; binAccRateVsRedshift(): NOTE: lots of duplication with binHaloMassVals()

function binAccRateVsRedshift, sP=sP, accMode=accMode, binWindow=BW

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
  
  wAm = accModeInds(at=at,accMode=accMode,sP=sP)
    
  ; reverse histogram parent IDs of all particles/tracers in this selection
  if sP.trMCPerCell gt 0 then gcIndOrig = mt.gcIndOrigTrMC
  if sP.trMCPerCell eq 0 then gcIndOrig = mt.gcIndOrig
  if sP.trMCPerCell lt 0 then gcIndOrig = mt.gcIndOrigTrVel
  
  maxHist = max(gcIndOrig)

  for i=0,n_tags(wAm)-1 do begin
    hist_type = histogram(gcIndOrig[wAm.(i)],min=0,max=maxHist,loc=loc,rev=rev_type)
    hist = mod_struct( hist, (tag_names(wAm))[i], hist_type )
    rev  = mod_struct( rev,  (tag_names(wAm))[i], rev_type )
  endfor
  
  gcIndOrig = !NULL
  
  ; load max temps, current tvir, tvir at accretion
  accTvir = gcSubsetProp(sP=sP,/accTvir,/accretionTimeSubset,accMode=accMode)
  maxTemp = gcSubsetProp(sP=sP,/maxPastTemp,/accretionTimeSubset,accMode=accMode)

  ; load group cat for subgroup masses, zero secondary so they don't get median binned
  galcat = galaxyCat(sP=sP)
  
  gc = loadGroupCat(sP=sP,/skipIDs)
  gcMasses = codeMassToLogMsun(gc.subgroupMass)
  gcIDs_sec = gcIDList(gc=gc,select='sec')
  gcMasses[gcIDs_sec] = -1
  gc = !NULL
  
  ; which galaxyCat types contribute to "total"?
  typeLabels = [ tag_names(galcat.types), 'total' ]
  
  totalInds = where( typeLabels ne 'total' ) ; where to take total from
  totalInd  = ( where( typeLabels eq 'total' ) )[0] ; where to store total  
  
  ; structure to store results
  tempHC = { cold : fltarr(nTimeBins,galcat.nGroups) ,$
             hot  : fltarr(nTimeBins,galcat.nGroups)  }
               
  foreach typeLabel,typeLabels do $  
    templateTypes = mod_struct(templateTypes, typeLabel, tempHC)
               
  ; for (1) galaxy and (2) halo, duplicate template for: accretion (inflow), outflow, and net
  galaxy = { accRate : templateTypes, outRate : templateTypes, netRate : templateTypes }
  halo   = { accRate : templateTypes, outRate : templateTypes, netRate : templateTypes }
                    
  ; loop over all tracked subgroups
  for i=0L,maxHist do begin
  
    if i mod (maxHist/10) eq 0 then print,string(float(i)/maxHist*100,format='(f5.2)')+'%'
    
    ; loop over types
    for k=0,n_tags(wAm)-1 do begin
      if (tag_names(wAm))[k] ne (tag_names(rev))[k] then message,'Check'
      
      ; enforce local timewindow
      loc_rev  = rev.(k)
	  
	if loc_rev[i] eq loc_rev[i+1] then continue ; no particles of type in this group, or whole group empty (secondary)
	  
      loc_inds = loc_rev[ loc_rev[i] : loc_rev[i+1]-1 ]
      
      ; accretion time defined as:
      ;  switching side of (rho,temp) cut or 0.15rvir crossing (headed in)
      ; outflow time defined as:
      ;  switching side of (rho,temp) cut or 0.15rvir crossing (headed out)
      
      ; --- GALAXY ACCRETION ---

      ; decide between using outTimeRT and accTimeRT (satisfying/failing the rho,temp cut, respectively)
      ; based on status at sP.snap (e.g. gmem by definition fails rho,temp cut)
      curTimeRT = max(at.accTimeRT)
      
      loc_atime = reform( at.accTimeRT[ (wAm.(k))[loc_inds] ] )
      
      ; those that fail at sP.snap (e.g. gmem), use satisfying the cut as the accretion time
      w = where( at.accTimeRT[ (wAm.(k))[loc_inds] ] eq curTimeRT, count)
      if count gt 0 then loc_atime[w] = at.outTimeRT[ (wAm.(k))[loc_inds[w]] ]
      
      ; for those that neither satisfy nor fail the cut at sP.snap (all stars, stars in inter) we
      ; cannot use the (rho,temp) cut, since don't know where they will fall when switching back
      ; to a gas parent, so just use the radial cut
      w_noCut = where( at.outTimeRT[ (wAm.(k))[loc_inds] ] ne max(at.outTimeRT) and $
                       at.accTimeRT[ (wAm.(k))[loc_inds] ] ne max(at.accTimeRT), count_noCut)
      if count_noCut gt 0 then loc_atime[w_noCut] = -1
		
      ; modify based on 0.15rvir crossing time (the latest / lowest redshift of the two criteria)
      r_crossing_time = reform( at.accTime[ sP.radIndGalAcc, (wAm.(k))[ loc_inds ] ] )
      w = where(r_crossing_time gt loc_atime, count)
      if count gt 0 then loc_atime[w] = r_crossing_time[w]
	  
      loc_atime = 1/loc_atime - 1 ; redshift
      loc_atime = redshiftToAgeFlat(loc_atime) * 1e9 ; yr
      
      ; split into hot and cold
      loc_maxTemp = maxTemp.(k)[loc_inds]
      loc_accTvir = accTvir.(k)[loc_inds]
      
      w_cold = where(10.0^loc_maxTemp / 10.0^loc_accTvir le maxTempTvirFac,$
                     count_cold,ncomp=count_hot,comp=w_hot)
      
      ; histogram in time
      if count_cold gt 0 then begin
        hh_cold = histogram(loc_atime[w_cold],min=startTime,max=curtime,binsize=timeWindow,loc=locc)
        galaxy.accRate.(k).cold[*,i] = hh_cold
      endif
      
      if count_hot gt 0 then begin
        hh_hot = histogram(loc_atime[w_hot],min=startTime,max=curtime,binsize=timeWindow,loc=loch)
        galaxy.accRate.(k).hot[*,i]  = hh_hot
      endif
      
      ; --- GALAXY OUTFLOW ---
      curTimeRT = max(at.outTimeRT)
      
      loc_atime = reform( at.outTimeRT[ (wAm.(k))[loc_inds] ] )
      
      ; those that satisfy at sP.snap (e.g. gal gas), use failing the cut as the outflow time
      w = where( at.outTimeRT[ (wAm.(k))[loc_inds] ] eq curTimeRT, count)
      if count gt 0 then loc_atime[w] = at.accTimeRT[ (wAm.(k))[loc_inds[w]] ]
      
      if count_noCut gt 0 then loc_atime[w_noCut] = -1 ; see above
      
      ; modify by 0.15rvir crossing time
      r_crossing_time = reform( at.outTime[ sP.radIndGalAcc, (wAm.(k))[ loc_inds ] ] )
      w = where(r_crossing_time gt loc_atime, count)
      if count gt 0 then loc_atime[w] = r_crossing_time[w]
      
      loc_atime = 1/loc_atime - 1 ; redshift
      loc_atime = redshiftToAgeFlat(loc_atime) * 1e9 ; yr
      
      ; histogram in time (using hot/cold split from above)
      if count_cold gt 0 then begin
        hh_cold = histogram(loc_atime[w_cold],min=startTime,max=curtime,binsize=timeWindow,loc=locc)
        galaxy.outRate.(k).cold[*,i] = hh_cold
      endif
      
      if count_hot gt 0 then begin
        hh_hot = histogram(loc_atime[w_hot],min=startTime,max=curtime,binsize=timeWindow,loc=loch)
        galaxy.outRate.(k).hot[*,i]  = hh_hot
      endif
            
      ; --- HALO ACCRETION ---
      ; accretion time defined as 1.0rvir crossing (headed in)
      ; outflow time defined as 1.0rvir crossing (headed out)
      loc_atime = reform( at.accTime[ sP.radIndHaloAcc, (wAm.(k))[ loc_inds ] ] )      
      loc_atime = 1/loc_atime - 1 ; redshift
      loc_atime = redshiftToAgeFlat(loc_atime) * 1e9 ; yr
      
      ; histogram in time (using hot/cold split from above)
      if count_cold gt 0 then begin
        hh_cold = histogram(loc_atime[w_cold],min=startTime,max=curtime,binsize=timeWindow,loc=locc)
        halo.accRate.(k).cold[*,i] = hh_cold
      endif
      
      if count_hot gt 0 then begin
        hh_hot = histogram(loc_atime[w_hot],min=startTime,max=curtime,binsize=timeWindow,loc=loch)
        halo.accRate.(k).hot[*,i]  = hh_hot
      endif
            
      ; --- HALO OUTFLOW ---
      loc_atime = reform( at.outTime[ sP.radIndHaloAcc, (wAm.(k))[ loc_inds ] ] )      
      loc_atime = 1/loc_atime - 1 ; redshift
      loc_atime = redshiftToAgeFlat(loc_atime) * 1e9 ; yr
      
      ; histogram in time (using hot/cold split from above)
      if count_cold gt 0 then begin
        hh_cold = histogram(loc_atime[w_cold],min=startTime,max=curtime,binsize=timeWindow,loc=locc)
        halo.outRate.(k).cold[*,i] = hh_cold
      endif
      
      if count_hot gt 0 then begin
        hh_hot = histogram(loc_atime[w_hot],min=startTime,max=curtime,binsize=timeWindow,loc=loch)
        halo.outRate.(k).hot[*,i]  = hh_hot
      endif
      	  
	; right now, accRate/outRate hold just particle/tracer counts, convert to mass rates
	factor = massPerPart * units.UnitMass_in_Msun / timeWindow
		
	for hotCold=0,1 do begin
	  galaxy.accRate.(k).(hotCold)[*,i] *= factor
	  galaxy.outRate.(k).(hotCold)[*,i] *= factor
	  halo.accRate.(k).(hotCold)[*,i] *= factor
	  halo.outRate.(k).(hotCold)[*,i] *= factor
		
	  ; compute net=(inflow-outflow) rates by type
	  galaxy.netRate.(k).(hotCold)[*,i] = $
	    galaxy.accRate.(k).(hotCold)[*,i] - galaxy.outRate.(k).(hotCold)[*,i]
	  halo.netRate.(k).(hotCold)[*,i] = $
	    halo.accRate.(k).(hotCold)[*,i] - halo.outRate.(k).(hotCold)[*,i]
		
	  ; calculate total rates (not by type)
	  galaxy.accRate.total.(hotCold)[*,i] += galaxy.accRate.(k).(hotCold)[*,i]
	  galaxy.outRate.total.(hotCold)[*,i] += galaxy.outRate.(k).(hotCold)[*,i]
	  galaxy.netRate.total.(hotCold)[*,i] += galaxy.netRate.(k).(hotCold)[*,i]
	  halo.accRate.total.(hotCold)[*,i] += halo.accRate.(k).(hotCold)[*,i]
	  halo.outRate.total.(hotCold)[*,i] += halo.outRate.(k).(hotCold)[*,i]
	  halo.netRate.total.(hotCold)[*,i] += halo.netRate.(k).(hotCold)[*,i]
	endfor

    endfor ; n_tags(wAm)
    
  endfor ; i

  ; for a number of halo mass bins, calculate median acc rates as a function of redshift
  timeBinCen = ( locc + timeWindow*0.5 ) / 1e9 ; centers of time bins for plotting, tage(Gyr)
  logMassBins = list([8.0,13.0],[9.0,9.5],[9.5,10.0],[10.0,10.5],[10.5,11.0],[11.0,11.5],[11.5,12.5])
  logMassNBins  = n_elements(logMassBins)
  
  ; structures to store the binned values
  tempHC = { hot  : fltarr(nTimeBins,logMassNbins) ,$
             cold : fltarr(nTimeBins,logMassNbins)  }
  
  foreach typeLabel,typeLabels do $  
    templateTypes = mod_struct(templateTypes, typeLabel, tempHC)
  
  ; for (1) galaxy and (2) halo, duplicate template for: accretion (inflow), outflow, and net
  galaxyMedian = { accRate : templateTypes, outRate : templateTypes, netRate : templateTypes }
  haloMedian   = { accRate : templateTypes, outRate : templateTypes, netRate : templateTypes }
    
  foreach logMassBin,logMassBins,j do begin
    ; select halos in this mass bin (primary only)
    w = where(gcMasses ge logMassBin[0] and gcMasses lt logMassBin[1],countHalo)
    if countHalo eq 0 then continue
    
    print,j,countHalo
  
    ; loop over each galaxyCat type (gal,gmem,stars,inter,bhs,total)
    for k=0,n_tags(galaxy.accRate)-1 do begin
      ; loop over each time bin
      for i=0,nTimeBins-1 do begin
	  ; loop  over both hot and cold modes
	  for hotCold=0,1 do begin
          ; median accretion, outflow and net rates
          galaxyMedian.accRate.(k).(hotCold)[i,j] = mean(galaxy.accRate.(k).(hotCold)[i,w])
	    galaxyMedian.outRate.(k).(hotCold)[i,j] = mean(galaxy.outRate.(k).(hotCold)[i,w])
	    galaxyMedian.netRate.(k).(hotCold)[i,j] = mean(galaxy.netRate.(k).(hotCold)[i,w])
		  
          haloMedian.accRate.(k).(hotCold)[i,j] = mean(halo.accRate.(k).(hotCold)[i,w])
	    haloMedian.outRate.(k).(hotCold)[i,j] = mean(halo.outRate.(k).(hotCold)[i,w])
	    haloMedian.netRate.(k).(hotCold)[i,j] = mean(halo.netRate.(k).(hotCold)[i,w])
        endfor
      endfor
    endfor
  
  endforeach
  
  r = {galaxy:galaxy, halo:halo, galaxyMedian:galaxyMedian, haloMedian:haloMedian, $
       timeRange:[startTime,curtime], logMassBins:logMassBins, timeBinCen:timeBinCen, gcMasses:gcMasses}

  ; save
  ;save,r,filename=saveFilename
  print,'SKIP Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  return, r  

end

; plotAccRateVsRedshift():

pro plotAccRateVsRedshift

  ; config
  sP = simParams(res=128,run='gadget',redshift=2.0)
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
      
    cgPlot,bar.timeBinCen,bar.galaxyMedian.netRate.total.cold[*,massBinInd],color=colors.cold,line=0,/overplot
    cgPlot,bar.timeBinCen,bar.galaxyMedian.netRate.total.hot[*,massBinInd],color=colors.hot,line=0,/overplot
    
    cgPlot,bar.timeBinCen,bar.haloMedian.netRate.total.cold[*,massBinInd],color=colors.cold,line=1,/overplot
    cgPlot,bar.timeBinCen,bar.haloMedian.netRate.total.hot[*,massBinInd],color=colors.hot,line=1,/overplot
      
    redshift_axis, xrange_time, yrange_rate, /ylog, sP=sP
  end_PS
  
  stop
  

end
