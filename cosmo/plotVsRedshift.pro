; plotVsRedshift.pro
; feedback - plots skipping tconst/tvircur/tviracc definitions in favor of redshift panels
; dnelson jun.2013

; plotRatesFracsInRedshift():

pro plotRatesFracsInRedshift

  ; config
  runs       = ['feedback','tracer','gadget']
  redshifts  = [3.0,2.0,1.0,0.0]
  res        = 256
  
  sgSelect   = 'pri'    ; subhalo type: primary, secondary, all
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
        haloMassBinValues(sP=sP_z.(j),sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow))
      bth_z = mod_struct(bth_z, 'redshift'+str(j), $
        binTmaxHistos(sP=sP_z.(j),sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow))
        
      ; for each redshift, make for all the accretion modes
      mbv_mode = {}
      
      foreach amCur,accModes,k do begin
        if amCur eq 'recycled' and sP_z.(j).gfmWinds eq 0 then continue ; skip recycled for nonWinds  
        mbv_mode = mod_struct(mbv_mode, 'mode'+str(k), $
          haloMassBinValues(sP=sP_z.(j),sgSelect=sgSelect,accMode=amCur,timeWindow=timeWindow))
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
  
  ; cold fraction (both,gmem)
  start_PS, sP.(0).(0).plotPath + 'coldFracRedshift.both-gmem.' + plotStr + '.eps', /big
    
    for zind=0,3 do begin
      cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_frac,/xs,/ys,$
        ytitle="Cold Fraction",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
        /noerase,pos=pos[zind]
    
      for i=0,n_tags(sP)-1 do begin
        xvals = mbv.(i).(zind).logMassBinCen
        yvals = smooth(mbv.(i).(zind).fracMedian.both_tVirAcc[tVirInd,*],sK,/nan)
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[0],/overplot ; both
        
        yvals = smooth(mbv.(i).(zind).fracMedian.gmem_tVirAcc[tVirInd,*],sK,/nan)
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
  
  ; accretion rates (both)
  start_PS, sP.(0).(0).plotPath + 'accRateRedshift.both.' + plotStr + '.eps', /big
    
    for zind=0,3 do begin
      cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_rate,/xs,/ys,$
        ytitle=textoidl("dM_{gas}/dt  M_{halo}^{-1} [_{ }Gyr^{-1 }]"),$
        xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
        /noerase,pos=pos[zind]
    
      for i=0,n_tags(sP)-1 do begin
        xvals = mbv.(i).(zind).logMassBinCen
        yvals = smooth(alog10(mbv.(i).(zind).hotMedian.both_tVirAcc[tVirInd,*] / $
                       10.0^mbv.(i).(zind).logMassBinCen * 1e9),sK,/nan) ; yr->Gyr
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[0],/overplot ; both, hot
        
        yvals = smooth(alog10(mbv.(i).(zind).coldMedian.both_tVirAcc[tVirInd,*] / $
                       10.0^mbv.(i).(zind).logMassBinCen * 1e9),sK,/nan) ; yr->Gyr
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[1],/overplot ; both, cold
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
        yvals = smooth(alog10(mbv.(i).(zind).hotMedian.gmem_tVirAcc[tVirInd,*] / $
                       10.0^mbv.(i).(zind).logMassBinCen * 1e9),sK,/nan) ; yr->Gyr
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[0],/overplot ; gmem, hot
        
        yvals = smooth(alog10(mbv.(i).(zind).coldMedian.gmem_tVirAcc[tVirInd,*] / $
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
  
  ; maximum temp histos (both,gmem) (one halo mass bin)
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
        yvals = float( bth.(i).(zind).hBothTmaxTviracc[massBinInd,*] ) / $
                total( bth.(i).(zind).hBothTmaxTviracc[j,*] )
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[0],/overplot ; both / GAL ???
        
        ; gmem
        yvals = float( bth.(i).(zind).hGmemTmaxTviracc[massBinInd,*] ) / $
                total( bth.(i).(zind).hGmemTmaxTviracc[j,*] )
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
  
  ; fraction of accretion by mode (both)
  pos_local = plot_pos(total=2*n_elements(redshifts),/gap)
  if accModes[0] ne 'all' then message,'Error: Not going to work.'
  
  start_PS, sP.(0).(0).plotPath + 'accRateFracsRedshift.both.' + plotStr + '.eps', xs=9, ys=12
    
    for zind=0,3 do begin
      cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_frac,/xs,/ys,$
        ytitle=textoidl("Fraction"),$
        xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
        /noerase,pos=pos_local[2*zind]
    
      ; hot (left column)
      for i=0,n_tags(sP)-1 do begin
        for j=1,n_tags(amv.(i).(zind))-1 do begin ; take ratio to first
        
          xvals = amv.(i).(zind).(0).logMassBinCen
          yvals = smooth(amv.(i).(zind).(j).hotMedian.both_tVirAcc[tVirInd,*] / $
                         amv.(i).(zind).(0).hotMedian.both_tVirAcc[tVirInd,*],sK,/nan)
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
          yvals = smooth(amv.(i).(zind).(j).coldMedian.both_tVirAcc[tVirInd,*] / $
                         amv.(i).(zind).(0).coldMedian.both_tVirAcc[tVirInd,*],sK,/nan)
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