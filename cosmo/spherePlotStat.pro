; spherePlotStat.pro
; feedback project - plotting of (statistical/quantitative) quantities derived from healpix spheres
; dnelson apr.2014

; processShellInfallFrac(): take statistics over all individual halos and bin by halo mass

function processShellInfallFrac, shellFrac=shellFrac, logMassBins=logMassBins
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  logMassNBins = n_elements(logMassBins)-1
  fieldNames   = tag_names(shellFrac)
  numFields    = n_tags(shellFrac)
    
  ; loop over each quantity that is measured as a function of radius
  for i=0,numFields-1 do begin
    ; skip meta-data
    if ( size(shellFrac.(i)) )[0] ne 2 then continue
      
    ; add a median line for this field to the original structure
    shellFrac = mod_struct( shellFrac, fieldNames[i]+'_median', $
                            fltarr(shellFrac.nRadFacs,logMassNBins) + !values.f_nan )
      
    for j=0,logMassNbins-1 do begin
      w = where(shellFrac.priMasses gt logMassBins[j] and shellFrac.priMasses le logMassBins[j+1],count)

      if count gt 0 then begin
        ; loop over radial facs and compute it
        for k=0,shellFrac.nRadFacs-1 do begin
          shellFrac.( n_tags(shellFrac)-1 )[k,j] = median( (shellFrac.(i))[k,w] )
        endfor
          
      endif ;count>0
    endfor ;j
  endfor ;i
    
  ; add quasi-static calculation for RelTh and AbsTh
  shellFrac = mod_struct( shellFrac, 'fracStatic_RelTh_median', $
                          fltarr(shellFrac.nRadFacs,logMassNBins) + !values.f_nan )
  shellFrac = mod_struct( shellFrac, 'fracStatic_AbsTh_median', $
                          fltarr(shellFrac.nRadFacs,logMassNBins) + !values.f_nan )
                            
  for j=0,logMassNbins-1 do begin
    w = where(shellFrac.priMasses gt logMassBins[j] and shellFrac.priMasses le logMassBins[j+1],count)
      
    if count gt 0 then begin
      ; loop over radial facs and compute it
      for k=0,shellFrac.nRadFacs-1 do begin
        shellFrac.fracStatic_RelTh_median[k,j] = $
          median( 1.0 - shellFrac.fracInfall_RelTh[k,w] - shellFrac.fracOutflow_RelTh[k,w] )
        shellFrac.fracStatic_AbsTh_median[k,j] = $
          median( 1.0 - shellFrac.fracInfall_AbsTh[k,w] - shellFrac.fracOutflow_AbsTh[k,w] )
      endfor
          
    endif ;count>0
  endfor ;j 

  return, shellFrac
end  

; plotShellInOutFracRateComp(): compare inflow/outflow fractions/rates 
;   vs halo mass and radius at fixed redshift

pro plotShellInOutFracRateComp
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  res      = 512
  redshift = 2.0
  runs     = ['gadget','tracer','feedback']
  
  ; plot config
  lines      = [0,2,1]    ; line style for NoTh/RelTh/AbsTh or In/Out/QS
  cInd       = 1          ; color index
  sK         = 3          ; smoothing kernel
  xrange     = [9.0,12.0] ; log halo mass
  xrangeR    = [0.25,1.5] ; r/rvir
  yrange     = [0.0,1.0]  ; fraction
  yrangeM    = [-1.5,2.5]  ; log mdot
  radInd     = 0          ; when only plotting one radius, use this one
  radInds    = [0,1,3,4]  ; when plotting 2x2 panels split by rad
  massInds   = indgen(26) ; when plotting vs radius, which mass bins to include (all)
  massLabels = [0,5,10,15,20,22,25] ; which to label
  
  plotStr   = ""
  simNames  = []
  simColors = []
  
  ; bin fractions into halo mass bins and make median lines
  logMassBins   = [linspace(9.0,10.9,20),11.0,11.1,11.3,11.5,11.7,11.9,12.21]
  logMassBinCen = 0.5 * (logMassBins + shift(logMassBins,-1))
  logMassBinCen = logMassBinCen[0:-2]  
  
  foreach run,runs do begin
    sP = simParams(res=res,run=run,redshift=redshift)
    
    ; save info for plots
    plotStr  += sP.plotPrefix + "_"
    simNames  = [simNames, sP.simName]
    simColors = [simColors, sP.colors[cInd]]
    
    ; load, process and store
    shellFrac  = calcShellInfallFrac( sP=sP )  
    shellFrac  = processShellInfallFrac( shellFrac=shellFrac, logMassBins=logMassBins )
    shellFracs = mod_struct( shellFracs, run, shellFrac )
  endforeach
    
  ; plot init
  set_plot,'ps'
  plotStr   += str(res) + '.z' + string(redshift,format='(f3.1)')
  massColors = reverse( sampleColorTable('red-temp', n_elements(massInds), bounds=[0.15,0.8]) )
  
  ; plot (1) DEBUG vs mass (one rad, one threshold, all runs, all points)
  start_PS, sP.plotPath+'fracInflowPointsVsMass.NoTh.radInd' + str(radInd) + '.' + plotStr + '.eps'
  
    ; plot
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Angular Covering Fraction of Inflow",$
      xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
      title=textoidl("r/r_{vir} = " + string(shellFracs.(0).radFacs[radInd],format='(f4.2)') + $
                     "  z = " + string(redshift,format='(f3.1)') + " (NoTh)")
      
    cgPlot,xrange,[0.5,0.5],line=0,color=cgColor('light gray'),/overplot
    
    for i=0,n_tags(shellFracs)-1 do begin
      cgPlot,shellFracs.(i).priMasses,shellFracs.(i).fracInfall_noTh[radInd,*],$
             psym=4,color=simColors[i],/overplot
      cgPlot,logMassBinCen,shellFracs.(i).fracInfall_noTh_median[radInd,*],$
             line=0,color=simColors[i],/overplot
    endfor
    
    ; legend
    legend,simNames,textcolors=simColors,/bottom,/right
    
  end_PS
  
  ; plot (1B) - direct test
  if 0 then begin
  pp = plot([0],/nodata,/buffer,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,$
      font_size = 20,$
      ytitle="Angular Covering Fraction of Inflow",$
      xtitle="$M_{halo} [_{ }log h^{-1} M_{sun }]$",$
      title=textoidl("r/r_{vir} = " + string(shellFracs.(0).radFacs[radInd],format='(f4.2)') + $
                     "  z = " + string(redshift,format='(f3.1)') + " (NoTh)"))
  
  pp.xthick = 1.6
  pp.ythick = 1.6
  pp.xminor = 2
  pp.yminor = 2
  
  pp = plot(xrange,[0.5,0.5],linestyle='-',color='light gray',overplot=pp)
  
  for i=0,n_tags(shellFracs)-1 do begin
    pp = plot(shellFracs.(i).priMasses,shellFracs.(i).fracInfall_noTh[radInd,*],$
              symbol='D',linestyle='none',color=simColors[i],/overplot)
    ;cgPlot,logMassBinCen,shellFracs.(i).fracInfall_noTh_median[radInd,*],$
    ;       line=0,color=simColors[i],/overplot
  endfor
  
  q = polygon( [10,10.5,10.4,9.8,10], [0.2,0.22,0.33,0.31,0.2], /data, $
               /fill_background, fill_color='blue', linestyle='none', fill_transparency=50 )
  q = polygon( [10,10.5,10.4,9.8,10]+0.2, [0.2,0.22,0.33,0.31,0.2]+0.03, /data, $
               /fill_background, fill_color='green', linestyle='none', fill_transparency=50 )   
  legend,simNames,textcolors=simColors,/bottom,/right 
  pp.save, sP.plotPath+'test.eps'
  endif ;0
  
  ; plot (2) - vs mass (one rad, all thresholds, all runs, median lines)
  start_PS, sP.plotPath + 'fracInflowVsMass.vsTh.radInd' + str(radInd) + '.' + plotStr + '.eps'
    ; plot    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Angular Covering Fraction of Inflow",$
      xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
      title=textoidl("r/r_{vir} = " + string(shellFracs.(0).radFacs[radInd],format='(f4.2)') + $
                     "  z = " + string(redshift,format='(f3.1)') + "")
      
    cgPlot,xrange,[0.5,0.5],line=0,color=cgColor('light gray'),/overplot
    
    for i=0,n_tags(shellFracs)-1 do begin
      cgPlot,logMassBinCen,smooth(shellFracs.(i).fracInfall_NoTh_median[radInd,*],sK),$
             line=lines[0],color=simColors[i],/overplot
      cgPlot,logMassBinCen,smooth(shellFracs.(i).fracInfall_RelTh_median[radInd,*],sK),$
             line=lines[1],color=simColors[i],/overplot
      cgPlot,logMassBinCen,smooth(shellFracs.(i).fracInfall_AbsTh_median[radInd,*],sK),$
             line=lines[2],color=simColors[i],/overplot
    endfor
    
    ; legends
    legend,simNames,textcolors=simColors,/bottom,/right
    legend,['NoTh','RelTh','AbsTh'],linestyle=lines,/bottom,/left
  end_PS
  
  ; plot (3) - vs mass 2x2 (multiple rad, all thresholds, all runs, median lines)
  start_PS, sP.plotPath + 'fracInflowVsMass.vsTh.rad2x2.' + plotStr + '.eps', /extrabig
    pos = plot_pos(col=2, row=2, /gap)
    
    for k=0,3 do begin
      ; plot    
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,pos=pos[k],noerase=(k gt 0),$
        ytitle=textoidl("f_\alpha ( Inflow )"),$
        xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
        title=textoidl("r/r_{vir} = " + string(shellFracs.(0).radFacs[radInds[k]],format='(f4.2)') + $
                       "  z = " + string(redshift,format='(f3.1)') + "")
      
      cgPlot,xrange,[0.5,0.5],line=0,color=cgColor('light gray'),/overplot
    
      for i=0,n_tags(shellFracs)-1 do begin
        cgPlot,logMassBinCen,smooth(shellFracs.(i).fracInfall_NoTh_median[radInds[k],*],sK),$
               line=lines[0],color=simColors[i],/overplot
        cgPlot,logMassBinCen,smooth(shellFracs.(i).fracInfall_RelTh_median[radInds[k],*],sK),$
               line=lines[1],color=simColors[i],/overplot
        cgPlot,logMassBinCen,smooth(shellFracs.(i).fracInfall_AbsTh_median[radInds[k],*],sK),$
               line=lines[2],color=simColors[i],/overplot
      endfor
    
      ; legends
      legend,simNames,textcolors=simColors,/bottom,/right
      legend,['NoTh','RelTh','AbsTh'],linestyle=lines,/bottom,/left
    endfor ;k
  end_PS
  
  ; plot (4) - vs mass 2x2 (multiple rad, RelTh, all runs, median lines, inflow+outflow+static)
  start_PS, sP.plotPath + 'fracInOutStatic.VsMass.RelTh.rad2x2.' + plotStr + '.eps', /extrabig
    pos = plot_pos(col=2, row=2, /gap)
    
    for k=0,3 do begin
      ; plot    
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,pos=pos[k],noerase=(k gt 0),$
        ytitle=textoidl("f_\alpha"),$
        xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
        title=textoidl("r/r_{vir} = " + string(shellFracs.(0).radFacs[radInds[k]],format='(f4.2)') + $
                       "  z = " + string(redshift,format='(f3.1)') + " (RelTh)")
      
      cgPlot,xrange,[0.5,0.5],line=0,color=cgColor('light gray'),/overplot
    
      for i=0,n_tags(shellFracs)-1 do begin
        cgPlot,logMassBinCen,smooth(shellFracs.(i).fracInfall_RelTh_median[radInds[k],*],sK),$
               line=lines[0],color=simColors[i],/overplot
        cgPlot,logMassBinCen,smooth(shellFracs.(i).fracOutflow_RelTh_median[radInds[k],*],sK),$
               line=lines[1],color=simColors[i],/overplot 
        ;cgPlot,logMassBinCen,smooth(shellFracs.(i).fracStatic_RelTh_median[radInds[k],*],sK),$
        ;       line=lines[2],color=simColors[i],/overplot
      endfor
    
      ; legends
      legend,simNames,textcolors=simColors,/bottom,/right
      legend,['Inflow','Outflow'],linestyle=lines[0:1],/bottom,/left
      ;legend,['Inflow','Outflow','Static'],linestyle=lines,/bottom,/left
    endfor ;k
  end_PS
  
  ; plot (5) - vs rad (many mass bins, RelTh, all runs, median lines)
  start_PS, sP.plotPath + 'fracInflowVsRad.RelTh.' + plotStr + '.eps', xs=4.5*n_tags(shellFracs), y=5.0
    pos = plot_pos(col=n_tags(shellFracs), row=1, /gap)
    
    for i=0,n_tags(shellFracs)-1 do begin
      ; plot    
      cgPlot,[0],[0],/nodata,xrange=xrangeR,yrange=yrange,/xs,/ys,pos=pos[i],noerase=(i gt 0),$
        ytitle="Angular Covering Fraction of Inflow",$
        xtitle=textoidl("r/r_{vir}"),$
        title=textoidl("z = " + string(redshift,format='(f3.1)') + ""),$
      xtickv=shellFracs.(0).radFacs,xticks=n_elements(shellFracs.(0).radFacs)-1,xminor=0
      
      cgPlot,xrange,[0.5,0.5],line=0,color=cgColor('light gray'),/overplot
    
      massStrs = []
      for j=0,n_elements(massInds)-1 do begin
        cgPlot,shellFracs.(i).radFacs,smooth(shellFracs.(i).fracInfall_RelTh_median[*,massinds[j]],1),$
               line=lines[0],color=massColors[j],/overplot
        ;cgPlot,shellFracs.(i).radFacs,smooth(shellFracs.(i).fracOutflow_RelTh_median[*,massinds[j]],1),$
        ;       line=lines[2],color=massColors[j],/overplot
               
        minMass = string(logMassBins[massInds[j]],format='(f4.1)')
        maxMass = string(logMassBins[massInds[j]+1],format='(f4.1)')
        massStrs = [massStrs, textoidl(minMass + " < M_{halo} < " + maxMass)]
        if logMassBins[massInds[j]] lt 10.0 then massStrs[-1] = " " + massStrs[-1] + " "
      endfor
    
      ; legends
      legend,simNames[i],textcolors=simColors[i],/bottom,/right
      legend,massStrs[massLabels],textcolors=massColors[massLabels],/bottom,/left
    endfor ;k
  end_PS
  
  ; plot (6) - RATE vs mass (one rad, in+out, all runs, median lines)
  Myr_to_yr = 1e6
  
  start_PS, sP.plotPath + 'rateInOutVsMass.radInd' + str(radInd) + '.' + plotStr + '.eps'
    ; plot    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeM,/xs,/ys,$
      ytitle=textoidl("log ( dM_{rad }/dt ) [_{ }h^{-1} M_{sun} yr^{-1 }]"),$
      xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
      title=textoidl("r/r_{vir} = " + string(shellFracs.(0).radFacs[radInd],format='(f4.2)') + $
                     "  z = " + string(redshift,format='(f3.1)') + "")
    
    cgPlot,xrange,[2.0,2.0],line=0,color=cgColor('light gray'),/overplot
    
    for i=0,n_tags(shellFracs)-1 do begin
      cgPlot,logMassBinCen,alog10(-smooth(shellFracs.(i).massFluxIn_median[radInd,*]/Myr_to_yr,sK)),$
             line=lines[0],color=simColors[i],/overplot
      cgPlot,logMassBinCen,alog10(smooth(shellFracs.(i).massFluxOut_median[radInd,*]/Myr_to_yr,sK)),$
             line=lines[2],color=simColors[i],/overplot
    endfor
    
    ; legends
    legend,simNames,textcolors=simColors,/bottom,/right
    legend,['Inflow','Outflow'],linestyle=lines[[0,2]],/top,/left
  end_PS
  
  ; plot (7) - RATE vs mass 2x2 (multiple rad, in+out, all runs, median lines)
  start_PS, sP.plotPath + 'rateInOutVsMass.rad2x2.' + plotStr + '.eps', /extrabig
    pos = plot_pos(col=2, row=2, /gap)
    
    for k=0,3 do begin
      ; plot    
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeM-0.5*(k ge 2),/xs,/ys,pos=pos[k],noerase=(k gt 0),$
        ytitle=textoidl("log ( dM_{rad }/dt ) [_{ }h^{-1} M_{sun} Myr^{-1 }]"),$
        xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
        title=textoidl("r/r_{vir} = " + string(shellFracs.(0).radFacs[radInds[k]],format='(f4.2)') + $
                       "  z = " + string(redshift,format='(f3.1)') + "")
    
      cgPlot,xrange,[2.0,2.0],line=0,color=cgColor('light gray'),/overplot
    
      for i=0,n_tags(shellFracs)-1 do begin
        cgPlot,logMassBinCen,alog10(-smooth(shellFracs.(i).massFluxIn_median[radInds[k],*]/Myr_to_yr,sK)),$
               line=lines[0],color=simColors[i],/overplot
        cgPlot,logMassBinCen,alog10(smooth(shellFracs.(i).massFluxOut_median[radInds[k],*]/Myr_to_yr,sK)),$
               line=lines[2],color=simColors[i],/overplot
      endfor
    
      ; legends
      legend,simNames,textcolors=simColors,/bottom,/right
      legend,['Inflow','Outflow'],linestyle=lines[[0,2]],/top,/left
    endfor ;k
  end_PS
  
  ; plot (8) - RATE vs rad (many mass bins, in+out, all runs, median lines)
  start_PS, sP.plotPath + 'rateInOutVsRad.' + plotStr + '.eps', xs=5.0*n_tags(shellFracs), y=2*5.0
    pos = plot_pos(col=n_tags(shellFracs), row=2, /gap)
    
    ; INFLOW
    for i=0,n_tags(shellFracs)-1 do begin
      ; plot    
      cgPlot,[0],[0],/nodata,xrange=xrangeR,yrange=yrangeM-0.5,/xs,/ys,pos=pos[i],noerase=(i gt 0),$
        ytitle=textoidl("log ( dM_{rad,IN }/dt ) [_{ }h^{-1} M_{sun} Myr^{-1 }]"),$
        xtitle=textoidl("r/r_{vir}"),$
        title=textoidl("z = " + string(redshift,format='(f3.1)') + ""),$
      xtickv=shellFracs.(0).radFacs,xticks=n_elements(shellFracs.(0).radFacs)-1,xminor=0
      
      cgPlot,xrange,[2.0,2.0],line=0,color=cgColor('light gray'),/overplot
    
      massStrs = []
      for j=0,n_elements(massInds)-1 do begin
        cgPlot,shellFracs.(i).radFacs,$
               alog10(-smooth(shellFracs.(i).massFluxIn_median[*,massinds[j]]/Myr_to_yr,sK)),$
               line=lines[0],color=massColors[j],/overplot
               
        minMass = string(logMassBins[massInds[j]],format='(f4.1)')
        maxMass = string(logMassBins[massInds[j]+1],format='(f4.1)')
        massStrs = [massStrs, textoidl(minMass + " < M_{halo} < " + maxMass)]
        if logMassBins[massInds[j]] lt 10.0 then massStrs[-1] = " " + massStrs[-1] + " "
      endfor
    
      ; legends
      legend,simNames[i],textcolors=simColors[i],/bottom,/left
      legend,massStrs[massLabels],textcolors=massColors[massLabels],/bottom,/right
    endfor ;k
    
    ; OUTFLOW
    for i=0,n_tags(shellFracs)-1 do begin
      ; plot    
      cgPlot,[0],[0],/nodata,xrange=xrangeR,yrange=yrangeM-0.5,/xs,/ys,pos=pos[3+i],/noerase,$
        ytitle=textoidl("log ( dM_{rad,OUT }/dt ) [_{ }h^{-1} M_{sun} Myr^{-1 }]"),$
        xtitle=textoidl("r/r_{vir}"),$
        title=textoidl("z = " + string(redshift,format='(f3.1)') + ""),$
      xtickv=shellFracs.(0).radFacs,xticks=n_elements(shellFracs.(0).radFacs)-1,xminor=0
      
      cgPlot,xrange,[2.0,2.0],line=0,color=cgColor('light gray'),/overplot
    
      massStrs = []
      for j=0,n_elements(massInds)-1 do begin
        cgPlot,shellFracs.(i).radFacs,$
               alog10(smooth(shellFracs.(i).massFluxOut_median[*,massinds[j]]/Myr_to_yr,sK)),$
               line=lines[2],color=massColors[j],/overplot
               
        minMass = string(logMassBins[massInds[j]],format='(f4.1)')
        maxMass = string(logMassBins[massInds[j]+1],format='(f4.1)')
        massStrs = [massStrs, textoidl(minMass + " < M_{halo} < " + maxMass)]
        if logMassBins[massInds[j]] lt 10.0 then massStrs[-1] = " " + massStrs[-1] + " "
      endfor
    
      ; legends
      legend,simNames[i],textcolors=simColors[i],/bottom,/left
      legend,['Inflow','Outflow'],linestyle=lines[[0,2]],/top,/right
      ;legend,massStrs[massLabels],textcolors=massColors[massLabels],/bottom,/right
    endfor ;k
  end_PS
  
  stop

end

; plotShellInOutFracRateVsRedshift(): compare inflow/outflow fractions/rates 
;   vs redshift at a few bins of halo mass and radius

pro plotShellInOutFracRateVsRedshift
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  res  = 512
  runs = ['gadget','tracer','feedback']
  
  ; plot config
  xrange     = [0.0,6.0]  ; redshift
  yrange     = [0.0,1.0]  ; fraction
  yrangeM    = [-2.5,2.5] ; log mdot
  radInd     = 0          ; when only plotting one radius, use this one
  radInds    = [0,1,3,4]  ; when plotting 2x2 panels split by rad
  massInd    = 4          ; when only plotting one massbin, use this one
  sK         = 3          ; smoothing kernel
  cInd       = 1          ; color index
  lines      = [0,2,1]
  massStrSkip = 2         ; label every other mass bin
  
  plotStr   = ""
  simNames  = []
  simColors = []
  
  foreach run,runs do begin
    sP = simParams(res=res,run=run,redshift=redshift)
    
    ; save info for plots
    plotStr  += sP.plotPrefix + "_"
    simNames  = [simNames, sP.simName]
    simColors = [simColors, sP.colors[cInd]]
  endforeach
  
  ; check for existing binned file (just to save time, load is ~5-10 minutes)
  runStr = ''
  foreach run,runs do runStr += '_' + str(run)
  sP = simParams(res=res,run=runs[0],redshift=0.0)
  saveFilename = sP.derivPath + 'hShells/shellFracRedshift.' + str(res) + runStr + '.snapFac2.sav'
    
  if file_test(saveFilename) then begin
    restore,saveFilename
  endif else begin
    ; compute new and save
    foreach run,runs do begin
      pp = calcShellValuesAcrossRedshifts(run=run,/getParams)
      sP = simParams(res=res,run=run)
    
      runShells = {}
      snaps     = lonarr( n_elements(pp.targetSnaps) )
      redshifts = fltarr( n_elements(pp.targetSnaps) )

      foreach snap, pp.targetSnaps, k do begin
        sP.snap = snap
      
        ; load
        shellFrac = calcShellInfallFrac(sP=sP,massBinCens=pp.massBinCens,/blindSearch)
    
        if n_elements(shellFrac) eq 0 then message,'Error: Will have index alignment issues, do all.'
        print,run,snap
      
        ; bin into pre-chosen pp.massBinCens
        shellFrac = processShellInfallFrac( shellFrac=shellFrac, logMassBins=[pp.massBinCens,13.5] )
        runShells = mod_struct( runShells, 'shellFrac_'+str(snap), shellFrac )
        
        snaps[k]     = snap
        redshifts[k] = snapNumToRedshift(snap=snap,sP=sP)
      endforeach
    
      ; store
      shellFracs   = mod_struct( shellFracs, run, runShells )
      runSnaps     = mod_struct( runSnaps, run, snaps )
      runRedshifts = mod_struct( runRedshifts, run, redshifts )
      pParams      = mod_struct( pParams, run, pp )
    endforeach
    
    ; save
    save,shellFracs,runSnaps,runRedshifts,pParams,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endelse
  
  ; plot init
  set_plot,'ps'
  massColors = reverse( sampleColorTable('red-temp', n_elements(pParams.(0).massBinCens), bounds=[0.15,0.8]) )
  nMassBins = n_elements( pParams.(0).massBinCens )
    
  ; plot (1) - vs redshift (one rad, one massBin, all thresholds, all runs, median lines)
  start_PS, sP.plotPath + 'fracInflow.VsRedshift.vsTh.radInd'+str(radInd)+'.massInd'+str(massInd)+ '.' + plotStr + '.eps'
    ; plot    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Angular Covering Fraction of Inflow",$
      xtitle=textoidl("redshift"),$
      title=textoidl("r/r_{vir} = " + string(shellFracs.(0).(0).radFacs[radInd],format='(f4.2)') + $
                     "  M = " + string(pParams.(0).massBinCens[massInd],format='(f4.1)') + "")
      
    cgPlot,xrange,[0.5,0.5],line=0,color=cgColor('light gray'),/overplot
    
    for i=0,n_tags(shellFracs)-1 do begin
      yy_noTh  = fltarr( n_elements(runRedshifts.(i)), nMassBins )
      yy_relTh = fltarr( n_elements(runRedshifts.(i)), nMassBins )
      yy_absTh = fltarr( n_elements(runRedshifts.(i)), nMassBins )
      
      for j=0,n_elements(runRedshifts.(i))-1 do begin
        yy_noTh[j,*]  = shellFracs.(i).(j).fracInfall_NoTh_median[radInd,*]
        yy_relTh[j,*] = shellFracs.(i).(j).fracInfall_RelTh_median[radInd,*]
        yy_absTh[j,*] = shellFracs.(i).(j).fracInfall_AbsTh_median[radInd,*]
      endfor
      
      cgPlot,runRedshifts.(i),smooth(yy_noTh[*,massInd],sK,/nan),line=lines[0],color=simColors[i],/overplot
      cgPlot,runRedshifts.(i),smooth(yy_relTh[*,massInd],sK,/nan),line=lines[1],color=simColors[i],/overplot
      cgPlot,runRedshifts.(i),smooth(yy_absTh[*,massInd],sK,/nan),line=lines[2],color=simColors[i],/overplot
    endfor
    
    ; legends
    legend,simNames,textcolors=simColors,/bottom,/right
    legend,['NoTh','RelTh','AbsTh'],linestyle=lines,/bottom,/left
  end_PS
  
  ; plot (2) - vs redshift (2x2 rad, one massBin,all thresholds, all runs, median lines)
  start_PS, sP.plotPath + 'fracInflow.VsRedshift.vsTh.rad2x2.massInd'+str(massInd)+'.'+plotStr + '.eps', /extrabig
    pos = plot_pos(col=2, row=2, /gap)
    
    for k=0,3 do begin
      ; plot    
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,pos=pos[k],noerase=(k gt 0),$
        ytitle=textoidl("f_\alpha ( Inflow )"),xtitle=textoidl("redshift"),$
        title=textoidl("r/r_{vir} = " + string(shellFracs.(0).(0).radFacs[radInds[k]],format='(f4.2)')) + $
                     "  M = " + string(pParams.(0).massBinCens[massInd],format='(f4.1)') + ""
      
      cgPlot,xrange,[0.5,0.5],line=0,color=cgColor('light gray'),/overplot
    
      for i=0,n_tags(shellFracs)-1 do begin
        yy_noTh  = fltarr( n_elements(runRedshifts.(i)), nMassBins )
        yy_relTh = fltarr( n_elements(runRedshifts.(i)), nMassBins )
        yy_absTh = fltarr( n_elements(runRedshifts.(i)), nMassBins )
      
        for j=0,n_elements(runRedshifts.(i))-1 do begin
          yy_noTh[j,*]  = shellFracs.(i).(j).fracInfall_NoTh_median[radInds[k],*]
          yy_relTh[j,*] = shellFracs.(i).(j).fracInfall_RelTh_median[radInds[k],*]
          yy_absTh[j,*] = shellFracs.(i).(j).fracInfall_AbsTh_median[radInds[k],*]
        endfor
        
        cgPlot,runRedshifts.(i),smooth(yy_noTh[*,massInd],sK,/nan),line=lines[0],color=simColors[i],/overplot
        cgPlot,runRedshifts.(i),smooth(yy_relTh[*,massInd],sK,/nan),line=lines[1],color=simColors[i],/overplot
        cgPlot,runRedshifts.(i),smooth(yy_absTh[*,massInd],sK,/nan),line=lines[2],color=simColors[i],/overplot
      endfor
    
      ; legends
      legend,simNames,textcolors=simColors,/bottom,/right
      legend,['NoTh','RelTh','AbsTh'],linestyle=lines,/bottom,/left
    endfor ;k
  end_PS
  
  ; plot (3) - vs redshift (2x2 rad, one massBin, RelTh, all runs, median lines, inflow+outflow+static)
  start_PS, sP.plotPath + 'fracInOutStatic.VsRedshift.RelTh.rad2x2.massInd'+str(massInd)+'.'+plotStr + '.eps', /extrabig
    pos = plot_pos(col=2, row=2, /gap)
    
    for k=0,3 do begin
      ; plot    
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,pos=pos[k],noerase=(k gt 0),$
        ytitle=textoidl("f_\alpha"),xtitle=textoidl("redshift"),$
        title=textoidl("r/r_{vir} = " + string(shellFracs.(0).(0).radFacs[radInds[k]],format='(f4.2)') + $
                       " (RelTh)") + $
                     "  M = " + string(pParams.(0).massBinCens[massInd],format='(f4.1)') + ""
      
      cgPlot,xrange,[0.5,0.5],line=0,color=cgColor('light gray'),/overplot
    
      for i=0,n_tags(shellFracs)-1 do begin
        yy_in   = fltarr( n_elements(runRedshifts.(i)), nMassBins )
        yy_out  = fltarr( n_elements(runRedshifts.(i)), nMassBins )
        yy_stat = fltarr( n_elements(runRedshifts.(i)), nMassBins )
      
        for j=0,n_elements(runRedshifts.(i))-1 do begin
          yy_in[j,*]   = shellFracs.(i).(j).fracInfall_RelTh_median[radInds[k],*]
          yy_out[j,*]  = shellFracs.(i).(j).fracOutflow_RelTh_median[radInds[k],*]
          yy_stat[j,*] = shellFracs.(i).(j).fracStatic_RelTh_median[radInds[k],*]
        endfor
        
        cgPlot,runRedshifts.(i),smooth(yy_in[*,massInd],sK,/nan),line=lines[0],color=simColors[i],/overplot
        cgPlot,runRedshifts.(i),smooth(yy_out[*,massInd],sK,/nan),line=lines[1],color=simColors[i],/overplot
        ;cgPlot,runRedshifts.(i),smooth(yy_stat[*,massInd],sK,/nan),line=lines[2],color=simColors[i],/overplot
      endfor
    
      ; legends
      legend,simNames,textcolors=simColors,/bottom,/right
      legend,['Inflow','Outflow'],linestyle=lines[0:1],/bottom,/left
      ;legend,['Inflow','Outflow','Static'],linestyle=lines,/bottom,/left
    endfor ;k
  end_PS
  
  ; plot (4) - vs redshift (one rad, many mass bins, RelTh, all runs, median lines)
  foreach indSet, list( [0,1], [2,3] ) do begin
  start_PS, sP.plotPath + 'fracInflow.VsRedshift.RelTh.radInds-'+str(indSet[0])+str(indSet[1])+$
    '.allMB.' + plotStr + '.eps', xs=4.5*3, y=2*5.0
    
    pos = plot_pos(col=n_tags(shellFracs), row=2, /gap)
    
    foreach radInd,radInds[indSet],m do begin
    
    for i=0,n_tags(shellFracs)-1 do begin
      ; bin in redshift
      yy_in   = fltarr( n_elements(runRedshifts.(i)), nMassBins )
      yy_out  = fltarr( n_elements(runRedshifts.(i)), nMassBins )
      yy_stat = fltarr( n_elements(runRedshifts.(i)), nMassBins )
      
      for j=0,n_elements(runRedshifts.(i))-1 do begin
        yy_in[j,*]   = shellFracs.(i).(j).fracInfall_RelTh_median[radInd,*]
        yy_out[j,*]  = shellFracs.(i).(j).fracOutflow_RelTh_median[radInd,*]
        yy_stat[j,*] = shellFracs.(i).(j).fracStatic_RelTh_median[radInd,*]
      endfor
        
      ; plot    
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
        pos=pos[n_tags(shellFracs)*m+i],noerase=(n_tags(shellFracs)*m+i gt 0),$
        ytitle=textoidl("f_\alpha")+" (Inflow)",xtitle=textoidl("redshift"),$
        title=textoidl("r/r_{vir} = " + string(shellFracs.(0).(0).radFacs[radInd],format='(f4.2)') + "")
      
      cgPlot,xrange,[0.5,0.5],line=0,color=cgColor('light gray'),/overplot
    
      massStrs = []
      for j=0,nMassBins-1 do begin
        cgPlot,runRedshifts.(i),smooth(yy_in[*,j],sK,/nan),line=lines[0],color=massColors[j],/overplot
               
        massStrs = [massStrs, textoidl("M_{halo} = " + string(pParams.(i).massBinCens[j],format='(f4.1)'))]
        if pParams.(i).massBinCens[j] lt 10.0 then massStrs[-1] = massStrs[-1] + " "
      endfor
    
      ; legends
      legend,simNames[i],textcolors=simColors[i],/bottom,/left
      legend,massStrs,textcolors=massColors,/bottom,/right
    endfor ;shellFracs
    
    endforeach ;radInds
  end_PS
  
  endforeach ;indSet
  
  ; plot (5) - RATE vs redshift (2x2 radius, one massBin, in+out, all runs, median lines)
  Myr_to_yr = 1e6
  
  for massInd=0,nMassBins-1 do begin
  start_PS, sP.plotPath + 'rateInOut.VsRedshift.rad2x2.massInd'+str(massInd)+'.'+plotStr+'.eps', /extrabig
    pos = plot_pos(col=2, row=2, /gap)
    
    for k=0,3 do begin
    ; plot    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeM,/xs,/ys,pos=pos[k],noerase=(k gt 0),$
      ytitle=textoidl("log ( dM_{rad }/dt ) [_{ }h^{-1} M_{sun} yr^{-1 }]"),$
      xtitle=textoidl("redshift"),$
      title=textoidl("r/r_{vir} = " + string(shellFracs.(0).(0).radFacs[radInds[k]],format='(f4.2)') + $
                     "  M = " + string(pParams.(0).massBinCens[massInd],format='(f4.1)'))
    
    cgPlot,xrange,[0.0,0.0],line=0,color=cgColor('light gray'),/overplot
    
    for i=0,n_tags(shellFracs)-1 do begin
      ; bin in redshift
      yy_in   = fltarr( n_elements(runRedshifts.(i)), nMassBins )
      yy_out  = fltarr( n_elements(runRedshifts.(i)), nMassBins )
      
      for j=0,n_elements(runRedshifts.(i))-1 do begin
        yy_in[j,*]   = shellFracs.(i).(j).massFluxIn_median[radInds[k],*] / Myr_to_yr
        yy_out[j,*]  = shellFracs.(i).(j).massFluxOut_median[radInds[k],*] / Myr_to_yr
      endfor
      
      cgPlot,runRedshifts.(i),alog10(-smooth(yy_in[*,massInd],sK,/nan)),$
             line=lines[0],color=simColors[i],/overplot
      cgPlot,runRedshifts.(i),alog10(smooth(yy_out[*,massInd],sK,/nan)),$
             line=lines[2],color=simColors[i],/overplot
    endfor
    
    ; legends
    legend,simNames,textcolors=simColors,/bottom,/right
    legend,['Inflow','Outflow'],linestyle=lines[[0,2]],/top,/left
    endfor ;k
  end_PS
  endfor ;massInd
  
  ; plot (7) - RATE vs redshift (one rad, many mass bins, in+out, all runs, median lines)
  foreach indSet, list( [0,1], [2,3] ) do begin
  start_PS, sP.plotPath + 'rateInOut.VsRedshift.radInds-'+str(indSet[0])+str(indSet[1])+'.allMB.'+plotStr + '.eps', xs=5.0*n_tags(shellFracs), y=2*5.0
    pos = plot_pos(col=n_tags(shellFracs), row=2, /gap)
    
    foreach radInd,radInds[indSet],m do begin
    
    for i=0,n_tags(shellFracs)-1 do begin
      ; bin in redshift
      yy_in   = fltarr( n_elements(runRedshifts.(i)), nMassBins )
      yy_out  = fltarr( n_elements(runRedshifts.(i)), nMassBins )
      
      for j=0,n_elements(runRedshifts.(i))-1 do begin
        yy_in[j,*]   = shellFracs.(i).(j).massFluxIn_median[radInd,*] / Myr_to_yr
        yy_out[j,*]  = shellFracs.(i).(j).massFluxOut_median[radInd,*] / Myr_to_yr
      endfor
      
      ; plot    
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeM+[-1.5,1.0],/xs,/ys,$
        pos=pos[n_tags(shellFracs)*m+i],noerase=(n_tags(shellFracs)*m+i gt 0),$
        ytitle=textoidl("log ( dM_{rad,IN/OUT }/dt ) [_{ }h^{-1} M_{sun} Myr^{-1 }]"),$
        xtitle=textoidl("redshift"),$
        title=textoidl("r/r_{vir} = " + string(shellFracs.(0).(0).radFacs[radInd],format='(f4.2)'))
      
      cgPlot,xrange,[0.0,0.0],line=0,color=cgColor('light gray'),/overplot
    
      massStrs = []
      for j=0,nMassBins-1 do begin
        cgPlot,runRedshifts.(i),alog10(-smooth(yy_in[*,j],sK,/nan)),line=lines[0],color=massColors[j],/overplot
        cgPlot,runRedshifts.(i),alog10(smooth(yy_out[*,j],sK,/nan)),line=lines[2],color=massColors[j],/overplot
                  
        massStrs = [massStrs, textoidl("M_{halo} = " + string(pParams.(i).massBinCens[j],format='(f4.1)'))]
        if pParams.(i).massBinCens[j] lt 10.0 then massStrs[-1] = massStrs[-1] + " "
      endfor
    
      ; legends
      legend,simNames[i],textcolors=simColors[i],/top,/left
      legend,massStrs[0:-1:massStrSkip],textcolors=massColors[0:-1:massStrSkip],/bottom,/right
      legend,['Inflow','Outflow'],linestyle=lines[[0,2]],/top,/right
    endfor ;i
    
    endforeach ;radInds
    
  end_PS
  endforeach ;indSet
  
  stop

end
