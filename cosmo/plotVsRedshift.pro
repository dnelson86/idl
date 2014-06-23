; plotVsRedshift.pro
; feedback - plots skipping tconst/tvircur/tviracc definitions and timewindow variations
;   in favor of redshift evolution and redshift panels
; dnelson jun.2014

; plotsVsRedshift(): plot fractional rate contributions by mode, as a function of redshift
;  and also at 4 discrete redshifts, as a function of halo mass

pro plotsVsRedshift
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  runs       = ['feedback','tracer','gadget']
  redshifts  = [5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.75, 0.5, 0.25, 0.0]
  res        = 512
  timeWindow = 500.0    ; consider accretion over this past time range (Myr)
                        ; 250.0 500.0 1000.0 "all" "tVir_tIGM" or "tVir_tIGM_bin"
  tVirInd    = 0        ; use Tmax/Tviracc=1 to separate hot vs. cold
  accModes   = ['all','smooth','clumpy','stripped','recycled']
  
  ; plot config
  massInds = [22]        ; for plot vs redshift
  zInds    = [4,6,8,12]  ; for 2x2 panels vs halo mass
  
  xrange_z     = [-0.15,5.15] ; redshift
  xrange_halo  = [9.0,12.0]   ; log halo mass
  yrange_frac  = [0.01,2.0]   ; log fraction
  yrange_frac2 = [0.0,1.0]    ; non-log fraction
  yrange_spec  = [-4.0,-1.0]  ; log specific net rate
  yrange_rate  = [1e-3,1e1]   ; net rate
  
  lines   = [4,0,1,2,3]  ; for each accmode ratio to all
  linesHC = [0,1]        ; hot,cold or gal,halo
  linesIO = [1,2,0]      ; in,out,net
  sK      = 3            ; smoothing kernel size
  psyms   = [-9,-16]     ; symbols for massbins (open/closed circles)
  cts     = ['brewerC-blues','brewerC-greens']
  cInd    = 1     ; color index into sP.colors
  symsize = 0.8   ; symbol size
  addGA   = 0     ; add gadget z=2 line to the one coldfrac plot?
  
  py_min  = 0.25  ; gray band min in fraction
  py_max  = 1.0   ; gray band max in fraction
  py_min2 = 0.25  ; 2nd gray band min in fraction
  py_max2 = 0.5   ; 2nd gray band max in fraction
  px_min  = 11.25 ; gray band min in halo mass (must correspond to massInds)
  px_max  = 11.5  ; gray band max in halo mass (must correspond to massInds)
  px_cl   = 'light gray'  ; mass bin
  py_cl   = [240,240,240] ; fraction, lighter
  py_cl2  = [230,230,230] ; fraction, medium
  
  ; load
  if accModes[0] ne 'all' then message,'Not going to work.'
  if runs[0] ne 'feedback' then message,'Maybe a problem.'
  if runs[-1] ne 'gadget' then message,'Problem, since we skip plotting the last.'
  
  tempFilename = '/n/home07/dnelson/plots/temp_plotsVsRedshift_'+str(n_elements(redshifts))+'.sav'
  if ~file_test(tempFilename) then begin
  
    foreach run,runs,i do begin
      sP_loc = simParams(res=res,run=run)
      
      sP_z  = {}
      mbv_z = {}
      nModes = 0
    
      ; loop over each redshift
      foreach redshift,redshifts,j do begin
      
        sP_mode  = {}
        mbv_mode = {}
        
        if run eq 'gadget' and redshift ne 2.0 then continue ; only want z=2 for gadget
      
        ; loop over each accretion mode
        foreach accMode,accModes,k do begin
          if accMode eq 'recycled' and sP_loc.gfmWinds eq 0 then continue ; skip recycled for nonWinds 

          sP_loc = simParams(res=res,run=run,redshift=redshift)
      
          mbv_loc = haloMassBinValues(sP=sP_loc,accMode=accMode,timeWindow=timeWindow,/search)
        
          if n_tags(mbv_loc) eq 0 then begin
            print,run,' ',redshift,' ',accMode,' (SKIP)'
            continue
          endif
          
          ; empty structure to hold ratios
          numMassBins = n_elements(mbv_loc.logMassBinCen)
          
          mode_ratio   = { gal_hot   : fltarr(numMassBins) + !values.f_nan ,$
                           gal_cold  : fltarr(numMassBins) + !values.f_nan ,$
                           halo_hot  : fltarr(numMassBins) + !values.f_nan ,$
                           halo_cold : fltarr(numMassBins) + !values.f_nan  } 
          specific_net = { gal_hot   : fltarr(numMassBins) + !values.f_nan ,$
                           gal_cold  : fltarr(numMassBins) + !values.f_nan ,$
                           halo_hot  : fltarr(numMassBins) + !values.f_nan ,$
                           halo_cold : fltarr(numMassBins) + !values.f_nan  }
          coldfrac_net = { gal  : fltarr(numMassBins) + !values.f_nan ,$
                           halo : fltarr(numMassBins) + !values.f_nan  }
                         
          ; extract fields we care about
          mbv_locAll = mbv_loc
          mbv_loc = {galaxyMedian  : mbv_loc.galaxyMedian  ,$
                     haloMedian    : mbv_loc.haloMedian    ,$
                     logMassBinCen : mbv_loc.logMassBinCen ,$
                     mode_ratio    : mode_ratio            ,$
                     specific_net  : specific_net          ,$
                     coldfrac_net  : coldfrac_net           }
                   
          print,run,' ',redshift,' ',accMode
          sP_mode  = mod_struct(sP_mode, 'mode'+str(k)+'_'+accMode, sP_loc)
          mbv_mode = mod_struct(mbv_mode, 'mode'+str(k)+'_'+accMode, mbv_loc)
          mbv_modeAll = mod_struct(mbv_modeAll, 'mode'+str(k)+'_'+accMode, mbv_locAll) ; not saved
        endforeach ; accModes,k
      
        if nModes eq 0 then nModes = n_tags(mbv_mode) ; set
        
        if n_tags(sP_mode) eq 0 or n_tags(sP_mode) lt nModes then begin
          print,'SKIP REDSHIFT: ',redshift,' ',run,' (missing mode)'
          continue
        endif
        
        ; load halo masses to re-bin median lines of some quantities
        gc = loadGroupCat(sP=sP_loc,/skipIDs)
        gcMasses = codeMassToLogMsun(gc.subgroupMass)
        
        maxHist = n_elements(mbv_locAll.galaxy.netRate.total.hot.num)
        gcMasses = gcMasses[0:maxHist-1]
          
        for k=0,nModes-1 do begin
          
          ratio_gal_hot   = reform(mbv_modeAll.(k).galaxy.netRate.total.hot.tVirAcc[tVirInd,*] / $
                                   mbv_modeAll.(0).galaxy.netRate.total.hot.tVirAcc[tVirInd,*])
          ratio_gal_cold  = reform(mbv_modeAll.(k).galaxy.netRate.total.cold.tVirAcc[tVirInd,*] / $
                                   mbv_modeAll.(0).galaxy.netRate.total.cold.tVirAcc[tVirInd,*])
          ratio_halo_hot  = reform(mbv_modeAll.(k).halo.netRate.total.hot.tVirAcc[tVirInd,*] / $
                                   mbv_modeAll.(0).halo.netRate.total.hot.tVirAcc[tVirInd,*])
          ratio_halo_cold = reform(mbv_modeAll.(k).halo.netRate.total.cold.tVirAcc[tVirInd,*] / $
                                   mbv_modeAll.(0).halo.netRate.total.cold.tVirAcc[tVirInd,*])

          gal_hot   = reform(mbv_modeAll.(k).galaxy.netRate.total.hot.tVirAcc[tVirInd,*]) / $
                      10.0^gcMasses * 1e9
          gal_cold  = reform(mbv_modeAll.(k).galaxy.netRate.total.cold.tVirAcc[tVirInd,*]) / $
                      10.0^gcMasses * 1e9
          halo_hot  = reform(mbv_modeAll.(k).halo.netRate.total.hot.tVirAcc[tVirInd,*]) / $
                      10.0^gcMasses * 1e9
          halo_cold = reform(mbv_modeAll.(k).halo.netRate.total.cold.tVirAcc[tVirInd,*]) / $
                      10.0^gcMasses * 1e9
                                   
          for m=0,n_elements(mbv_loc.logMassBinCen)-1 do begin
            ; restrict medians to primary only (num>0)
            w = where(gcMasses gt mbv_modeAll.(k).logMassBins[m] and $
                      gcMasses le mbv_modeAll.(k).logMassBins[m+1] and $
                      (mbv_modeAll.(0).galaxy.coldFrac.total.num+$
                       mbv_modeAll.(0).halo.coldFrac.total.num) gt 0,count)
                      
            if count eq 0 then continue
            if max(w) ge maxHist then message,'Error'
            
            ; (1) mode_ratio: bin ratios of modes vs 'all' mode
            mbv_mode.(k).mode_ratio.gal_hot[m]   = median( ratio_gal_hot[w] )
            mbv_mode.(k).mode_ratio.gal_cold[m]  = median( ratio_gal_cold[w] )
            mbv_mode.(k).mode_ratio.halo_hot[m]  = median( ratio_halo_hot[w] )
            mbv_mode.(k).mode_ratio.halo_cold[m] = median( ratio_halo_cold[w] )
            
            ; (2) specific_net and coldfrac_net: specific net rates/coldfracs, galaxy by galaxy (yr -> Gyr)
            mbv_mode.(k).specific_net.gal_hot[m]   = median( gal_hot[w] )
            mbv_mode.(k).specific_net.gal_cold[m]  = median( gal_cold[w] )
            mbv_mode.(k).specific_net.halo_hot[m]  = median( halo_hot[w] )
            mbv_mode.(k).specific_net.halo_cold[m] = median( halo_cold[w] )
            
            mbv_mode.(k).coldfrac_net.gal[m]  = median( gal_cold[w] / (gal_hot[w]+gal_cold[w]) )
            mbv_mode.(k).coldfrac_net.halo[m] = median( halo_cold[w] / (halo_hot[w]+halo_cold[w]) )
          endfor ;m
          
        endfor ; k
        
        sP_z  = mod_struct(sP_z, 'redshift'+str(j)+'_'+str(round(redshift*100)), sP_mode)
        mbv_z = mod_struct(mbv_z, 'redshift'+str(j)+'_'+str(round(redshift*100)), mbv_mode)
      endforeach ; redshifts,j
      
      ; put this mode collection into mbv, once per run, and sP collection also
      mbv = mod_struct(mbv, 'mbv'+str(i)+'_'+run, mbv_z)
      sP  = mod_struct(sP, 'sP'+str(i)+'_'+run, sP_z)
    endforeach ;runs,i
  
    save,mbv,sP,filename=tempFilename
    print,'Saved TEMP file: ['+tempFilename+']'
  endif else begin
    restore,tempFilename,/verbose
  endelse
  
  ; strings
  if isnumeric(timeWindow) then twStr = str(fix(timeWindow))
  if ~isnumeric(timeWindow) then twStr = "-"+timeWindow
  
  plotStr   = ''
  simNames  = []
  simColors = []
  simNames2  = []
  simColors2 = []
  
  for i=0,n_elements(runs)-2 do begin
    plotStr   = plotStr + sP.(i).(0).(0).plotPrefix + '.'
    simNames  = [simNames, sP.(i).(0).(0).simName]
    simColors = [simColors, sP.(i).(0).(0).colors[cInd]]
    simNames2  = [simNames2, sP.(i).(0).(0).simName, sP.(i).(0).(0).simName]
    simColors2 = [simColors2, sampleColorTable(cts[i],2,bounds=[0.4,0.9])]
  endfor

  plotStr  += str(res) + '_' + str(sP.(0).(0).(0).snap) + '_tw' + twStr             
  plotStrMB = plotStr + '_mB-' + strjoin(str(massInds),",")
  print,'Plotting mass bin centered on: '+str(mbv.(0).(0).(0).logMassBinCen[massInds])
  print,'Plotting redshifts: '+strjoin(redshifts[zInds]," ")
             
  ; plots vs redshift (all accModes together, @ one halo mass bin)
  ; ---------------------------------------
  pos_single = [0.15,0.11,0.94,0.87] ; single plot with second x-axis above
  
  ; plot (1) - fractional rates, gal hot/cold (all modes at once)
  start_PS, sP.(0).(0).(0).plotPath + 'rateFracsVsRedshift.galaxy.' + plotStrMB + '-allModes.eps'
  
    yrange2 = yrange_frac * [1.0,0.7]

    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=5,ys=5,/ylog,/noerase,pos=pos_single
      
    cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
      [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
    cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
      [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
        
    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=9,/ys,/ylog,yminor=0,$
      ytitle=textoidl("(dM/dt)_{mode} / (dM/dt)_{total}"),xtitle="Redshift",/noerase,pos=pos_single
    universeage_axis, xrange_z, yrange2, /ylog
        
    ; loop over runs
    for i=0,n_tags(sP)-2 do begin
      ; loop over mass bins (just one now)
      massStrs = []
      foreach massInd,massInds,k do begin
        
        massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
          
        ; loop over each accMode
        for j=1,n_tags(sP.(i).(0))-1 do begin
        
          xvals           = []
          yvals_gal_hot   = []
          yvals_gal_cold  = []
          yvals_halo_hot  = []
          yvals_halo_cold = []
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals = [xvals, sP.(i).(m).(j).redshift]
            yvals_gal_hot   = [yvals_gal_hot,   mbv.(i).(m).(j).mode_ratio.gal_hot[massInd]]
            yvals_gal_cold  = [yvals_gal_cold,  mbv.(i).(m).(j).mode_ratio.gal_cold[massInd]]
            yvals_halo_hot  = [yvals_halo_hot,  mbv.(i).(m).(j).mode_ratio.halo_hot[massInd]]
            yvals_halo_cold = [yvals_halo_cold, mbv.(i).(m).(j).mode_ratio.halo_cold[massInd]]
          endfor ;m
        
          ; smooth and set zeros to a plottable value
          yvals_gal_hot   = smooth(yvals_gal_hot,sK,/nan) < 1.0
          yvals_gal_cold  = smooth(yvals_gal_cold,sK,/nan) < 1.0
          yvals_halo_hot  = smooth(yvals_halo_hot,sK,/nan) < 1.0
          yvals_halo_cold = smooth(yvals_halo_cold,sK,/nan) < 1.0
      
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal_hot,color=simColors2[i*2+0],$
            line=lines[j],psym=-9,symsize=0.8,/overplot ;psym=psyms[k],
          cgPlot,xvals,yvals_gal_cold,color=simColors2[i*2+1],$
            line=lines[j],psym=-16,symsize=0.8,/overplot ;psym=psyms[k],
          
        endfor ; accModes 
      endforeach ; massBins
    endfor ; runs
      
    ;legend,textoidl("M_{halo} = "+string(massStrs,format='(f4.1)')),/top,/right
    legend,accModes[1:*],linestyle=lines[1:*],pos=[1.2,yrange2[0]*3.4]
      
    tvirStrs = textoidl(["T_{max} > T_{vir,acc}","T_{max} < T_{vir,acc}",$
                         "T_{max} > T_{vir,acc}","T_{max} < T_{vir,acc}"])
    legend,tvirStrs+" ("+simNames2+")",textcolor=simColors2,/bottom,/right,spacing=1.6
  
  end_PS
    
  ; plot (2) - fractional rates, gal hot/cold (2x2 panels, one per mode)
  start_PS, sP.(0).(0).(0).plotPath + 'rateFracsVsRedshiftPanels.galaxy.' + plotStrMB + '-allModes.eps', /extrabig
  
    pos = plot_pos(rows=2,cols=2,/gap)
    
    ; loop over each accMode
    for j=1,n_tags(sP.(0).(0))-1 do begin
    
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_frac,xs=5,ys=5,/ylog,/noerase,$
        pos=pos[j-1]-[0,0.02,0,0.04]
        
      cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
        [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
      cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
        [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
        
      barStr = str(round(py_max2*100))+'% - '+str(round(py_max*100))+'%'
      if j eq 3 then $
        cgText, xrange_z[1]-1.0, py_min+0.39, barStr, alignment=0.5, color=cgColor('gray')
        
      barStr = str(round(py_min*100))+'% - '+str(round(py_max2*100))+'% '
      if j eq 3 then $
        cgText, xrange_z[1]-1.0, py_min+0.06, barStr, alignment=0.5, color=cgColor('gray')
        
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_frac,xs=9,/ys,/ylog,yminor=0,$
        ytitle=textoidl("(dM/dt)_{"+accModes[j]+"} / (dM/dt)_{total}"),xtitle="Redshift",/noerase,$
        pos=pos[j-1]-[0,0.02,0,0.04]
      universeage_axis, xrange_z, yrange_frac, /ylog
    
      ; loop over runs
      for i=0,n_tags(sP)-2 do begin
        ; accMode does not exist for this run? then skip
        if j ge n_tags(sP.(i).(0)) then continue
        
        ; loop over mass bins (just one now)
        massStrs = []
        foreach massInd,massInds,k do begin
        
          massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
          
          xvals           = []
          yvals_gal_hot   = []
          yvals_gal_cold  = []
          yvals_halo_hot  = []
          yvals_halo_cold = []
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals = [xvals, sP.(i).(m).(j).redshift]
            yvals_gal_hot   = [yvals_gal_hot,   mbv.(i).(m).(j).mode_ratio.gal_hot[massInd]]
            yvals_gal_cold  = [yvals_gal_cold,  mbv.(i).(m).(j).mode_ratio.gal_cold[massInd]]
            yvals_halo_hot  = [yvals_halo_hot,  mbv.(i).(m).(j).mode_ratio.halo_hot[massInd]]
            yvals_halo_cold = [yvals_halo_cold, mbv.(i).(m).(j).mode_ratio.halo_cold[massInd]]
          endfor ;m
        
          ; smooth and set zeros to a plottable value
          yvals_gal_hot   = smooth(yvals_gal_hot,sK) < 1.0
          yvals_gal_cold  = smooth(yvals_gal_cold,sK) < 1.0
          yvals_halo_hot  = smooth(yvals_halo_hot,sK) < 1.0
          yvals_halo_cold = smooth(yvals_halo_cold,sK) < 1.0
      
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal_hot,color=simColors[i],psym=psyms[k],line=linesHC[0],/overplot
          cgPlot,xvals,yvals_gal_cold,color=simColors[i],psym=psyms[k],line=linesHC[1],/overplot
          
        endforeach ; massBins
      endfor ; runs
    
      if j eq 4 then $
      legend,textoidl(["T_{max} > T_{vir,acc}","T_{max} < T_{vir,acc}"]),$
        linestyle=linesHC,color=cgColor('black'),textcolor=cgColor('black'),/bottom,/right
      if j eq 2 then $
      legend,simNames,textcolor=simColors,/bottom,/left
    
    endfor ; accModes
  
  end_PS
  
  ; plot (2b) - coldfrac_net vs redshift, galaxy only, all accModes on plot
  start_PS, sP.(0).(0).(0).plotPath + 'coldFracNetVsRedshift.galaxy.' + plotStrMB + '-allModes.eps'
  
    yrange2 = yrange_frac * [1.0,0.7]

    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=5,ys=5,/ylog,/noerase,pos=pos_single
      
    cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
      [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
    cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
      [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
        
    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=9,/ys,/ylog,yminor=0,$
      ytitle=textoidl("(dM/dt)_{T_{max}<T_{vir,acc}} / (dM/dt)_{total,mode }"),$
      xtitle="Redshift",/noerase,pos=pos_single
    universeage_axis, xrange_z, yrange2, /ylog
        
    ; loop over runs
    for i=0,n_tags(sP)-2 do begin
      ; loop over mass bins (just one now)
      massStrs = []
      foreach massInd,massInds,k do begin
        
        massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
          
        ; loop over each accMode
        for j=0,n_tags(sP.(i).(0))-1 do begin
        
          xvals      = []
          yvals_gal  = []
          yvals_halo = []
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals = [xvals, sP.(i).(m).(j).redshift]
            
            yvals_gal  = [yvals_gal,  mbv.(i).(m).(j).coldfrac_net.gal[massInd]]
            yvals_halo = [yvals_halo, mbv.(i).(m).(j).coldfrac_net.halo[massInd]]
          endfor ;m
        
          ; smooth and set zeros to a plottable value
          yvals_gal  = smooth(yvals_gal,sK,/nan) < 1.0
          yvals_halo = smooth(yvals_halo,sK,/nan) < 1.0
      
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal,color=simColors[i],psym=psyms[1],line=lines[j],symsize=symsize,/overplot
          
        endfor ; accModes 
      endforeach ; massBins
    endfor ; runs
      
    ;legend,textoidl("M_{halo} = "+string(massStrs,format='(f4.1)')),/top,/right
    legend,accModes,linestyle=lines,linesize=0.6,/bottom,/right
    legend,simNames,textcolor=simColors,/top,/left
  
  end_PS
  
  ; plots vs redshift (one plot per accMode, @ one halo mass bin)
  ; --------------------------------------------------------
  temp = {}
  
  for j=0,n_tags(sP.(1).(0))-1 do begin
    plotStrAM = plotStrMB + '_am-'+accModes[j] 
    
    ; plot (3) - net rates vs redshift, galaxy hot/cold
    start_PS, sP.(0).(0).(0).plotPath + 'netRateVsRedshift.galaxy.' + plotStrAM + '.eps'
    
      yrange2 = yrange_rate
      if accModes[j] eq 'all' then yrange2 *= 4
      
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=9,/ys,/ylog,yminor=0,$
        ytitle=textoidl("dM_{gas,"+accModes[j]+" }/dt [_{ }M_{sun }/_{ }yr_{ }]"),$
        xtitle="Redshift",/noerase,pos=pos_single
        
      universeage_axis, xrange_z, yrange2, /ylog
      
      ; loop over runs      
      for i=0,n_tags(sP)-2 do begin
        ; loop over mass bins (just one now)
        massStrs = []
        foreach massInd,massInds,k do begin
        
          massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
        
          xvals           = []
          yvals_gal_hot   = []
          yvals_gal_cold  = []
          yvals_halo_hot  = []
          yvals_halo_cold = []
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals = [xvals, sP.(i).(m).(j).redshift]
            yvals_gal_hot   = [yvals_gal_hot, $
              mbv.(i).(m).(j).galaxyMedian.netRate.total.hot.tVirAcc[tVirInd,massInd]]
            yvals_gal_cold  = [yvals_gal_cold, $
              mbv.(i).(m).(j).galaxyMedian.netRate.total.cold.tVirAcc[tVirInd,massInd]]
            yvals_halo_hot  = [yvals_halo_hot, $
              mbv.(i).(m).(j).haloMedian.netRate.total.hot.tVirAcc[tVirInd,massInd]]
            yvals_halo_cold = [yvals_halo_cold, $
              mbv.(i).(m).(j).haloMedian.netRate.total.cold.tVirAcc[tVirInd,massInd]]
          endfor ;m
          
          ;temp = mod_struct( temp, 'total_gal_'+sP.(i).(0).(0).run+'_'+accModes[j], (yvals_gal_hot+yvals_gal_cold)/units.hubbleParam)
        
          ; smooth and set zeros to a plottable value
          yvals_gal_hot   = smooth(yvals_gal_hot/units.hubbleParam > yrange2[0]/10,sK,/nan)
          yvals_gal_cold  = smooth(yvals_gal_cold/units.hubbleParam > yrange2[0]/10,sK,/nan)
          yvals_halo_hot  = smooth(yvals_halo_hot/units.hubbleParam > yrange2[0]/10,sK,/nan)
          yvals_halo_cold = smooth(yvals_halo_cold/units.hubbleParam > yrange2[0]/10,sK,/nan)
        
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal_hot,color=simColors[i],psym=psyms[k],line=linesHC[0],/overplot
          cgPlot,xvals,yvals_gal_cold,color=simColors[i],psym=psyms[k],line=linesHC[1],/overplot
          
          ;cgPlot,xvals,yvals_halo_hot,color=simColors[i],psym=psyms[k],line=linesHC[0],/overplot
          ;cgPlot,xvals,yvals_halo_cold,color=simColors[i],psym=psyms[k],line=linesHC[1],/overplot
        endforeach ; massBins
        
      endfor ;i
      
      ;legend,textoidl("M_{halo} = "+string(massStrs,format='(f4.1)')),/top,/right
      legend,textoidl(["T_{max} > T_{vir,acc}","T_{max} < T_{vir,acc}"]),$
        linestyle=linesHC,color=cgColor('black'),textcolor=cgColor('black'),/bottom,/right
      legend,simNames,textcolor=simColors,/top,/left
    
    end_PS
    
    ; plot (4) - net/in/out rates vs redshift, galaxy hot/cold
    start_PS, sP.(0).(0).(0).plotPath + 'accOutNetRatesVsRedshift.galaxy.' + plotStrAM + '.eps'
    
      yrange2 = yrange_rate*[10,3]
      if accModes[j] eq 'all' then yrange2 *= 4
      
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=9,/ys,/ylog,yminor=0,$
        ytitle=textoidl("dM_{gas,"+accModes[j]+" }/dt [_{ }M_{sun }/_{ }yr_{ }]"),$
        xtitle="Redshift",/noerase,pos=pos_single
        
      universeage_axis, xrange_z, yrange2, /ylog
        
      ; loop over runs      
      for i=0,n_tags(sP)-2 do begin
        ; loop over mass bins (just one now)
        massStrs = []
        foreach massInd,massInds,k do begin
        
          massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
        
          xvals          = []
          yvals_gal_acc  = []
          yvals_gal_out  = []
          yvals_gal_net  = []
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals = [xvals, sP.(i).(m).(j).redshift]
            yvals_gal_acc = [yvals_gal_acc, $
              mbv.(i).(m).(j).galaxyMedian.accRate.total.hot.tVirAcc[tVirInd,massInd] + $
              mbv.(i).(m).(j).galaxyMedian.accRate.total.cold.tVirAcc[tVirInd,massInd]]
            yvals_gal_out = [yvals_gal_out, $
              mbv.(i).(m).(j).galaxyMedian.outRate.total.hot.tVirAcc[tVirInd,massInd] + $
              mbv.(i).(m).(j).galaxyMedian.outRate.total.cold.tVirAcc[tVirInd,massInd]]
            yvals_gal_net = [yvals_gal_net, $
              mbv.(i).(m).(j).galaxyMedian.netRate.total.hot.tVirAcc[tVirInd,massInd] + $
              mbv.(i).(m).(j).galaxyMedian.netRate.total.cold.tVirAcc[tVirInd,massInd]]
          endfor ;m
        
          ; smooth and set zeros to a plottable value
          yvals_gal_acc = smooth(yvals_gal_acc/units.hubbleParam,sK,/nan)
          yvals_gal_out = smooth(yvals_gal_out/units.hubbleParam,sK,/nan)
          yvals_gal_net = smooth(yvals_gal_net/units.hubbleParam,sK,/nan)
          
          temp = mod_struct( temp, 'ratio_gal_'+sP.(i).(0).(0).run+'_'+accModes[j], (yvals_gal_out/yvals_gal_acc) )
        
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal_acc,color=simColors[i],psym=psyms[k],line=linesIO[0],/overplot
          cgPlot,xvals,yvals_gal_out,color=simColors[i],psym=psyms[k],line=linesIO[1],/overplot
          cgPlot,xvals,yvals_gal_net,color=simColors[i],psym=psyms[k],line=linesIO[2],/overplot
          
        endforeach ; massBins
        
      endfor ;i
      
      ;legend,textoidl("M_{halo} = "+string(massStrs,format='(f4.1)')),/top,/right
      legend,textoidl(["acc","out","net"]),$
        linestyle=linesIO,color=cgColor('black'),textcolor=cgColor('black'),/bottom,/right
      legend,simNames,textcolor=simColors,/top,/left
    
    end_PS
    
    ; plot (5) - net cold fractions vs redshift, galaxy hot/cold
    start_PS, sP.(0).(0).(0).plotPath + 'coldFracNetVsRedshift.galaxy-halo.' + plotStrAM + '.eps'
    
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_frac,xs=5,ys=5,/ylog,/noerase,$
        pos=pos_single
        
      cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
        [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
      cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
        [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
    
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_frac,xs=9,/ys,/ylog,yminor=0,$
        ytitle=textoidl("(dM/dt)_{T_{max}<T_{vir,acc}} / (dM/dt)_{total,"+accModes[j]+" }"),$
        xtitle="Redshift",/noerase,pos=pos_single
        
      universeage_axis, xrange_z, yrange_frac, /ylog
        
      ; loop over runs      
      for i=0,n_tags(sP)-2 do begin
        ; loop over mass bins (just one now)
        massStrs = []
        foreach massInd,massInds,k do begin
        
          massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
        
          xvals      = []
          yvals_gal  = []
          yvals_halo = []
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals = [xvals, sP.(i).(m).(j).redshift]
            
            yvals_gal  = [yvals_gal,  mbv.(i).(m).(j).coldfrac_net.gal[massInd]]
            yvals_halo = [yvals_halo, mbv.(i).(m).(j).coldfrac_net.halo[massInd]]
          endfor ;m
        
          ; smooth and set zeros to a plottable value
          yvals_gal  = smooth(yvals_gal,sK,/nan)
          yvals_halo = smooth(yvals_halo,sK,/nan)
        
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal,color=simColors[i],psym=psyms[k],line=linesHC[0],/overplot
          cgPlot,xvals,yvals_halo,color=simColors[i],psym=psyms[k],line=linesHC[1],/overplot
          
        endforeach ; massBins
        
      endfor ;i
      
      ;legend,textoidl("M_{halo} = "+string(massStrs,format='(f4.1)')),/top,/right
      legend,['galaxy','halo'],linestyle=linesHC,$
        color=cgColor('black'),textcolors=['black','black'],linesize=0.4,box=0,/bottom,/right
      legend,simNames,textcolor=simColors,/top,/right
  
    end_PS

  endfor ;accModes,j
     
  ; plots vs halo mass (all accModes together, one plot per redshift)
  ; ------------------------------------------------------------
  foreach zInd,zInds,k do begin
  
    ; plot (6) - fraction of accretion by mode, galaxy hot/cold
    start_PS, sP.(0).(0).(0).plotPath + $
      'rateFracs4Redshifts.zInd-'+str(zInd)+'.galaxy.'+plotStr+'-allModes.eps', /extrabig
      
      pos = plot_pos(rows=2,cols=2,/gap)
      
      ; loop over accretion modes
      for m=1,n_elements(accModes)-1 do begin
        
        cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_frac,xs=5,ys=5,/ylog,/noerase,$
          pos=pos[m-1]-[0,0.02,0,0.04]
          
        cgPolygon,[xrange_halo[0],xrange_halo[1],xrange_halo[1],xrange_halo[0],xrange_halo[0]],$
          [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
        cgPolygon,[xrange_halo[0],xrange_halo[1],xrange_halo[1],xrange_halo[0],xrange_halo[0]],$
          [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
        
        cgPolygon,[px_min,px_max,px_max,px_min,px_min],$
          [yrange_frac[0],yrange_frac[0],yrange_frac[1],yrange_frac[1],yrange_frac[0]],$
          color=cgColor(px_cl)
        
        barStr = str(round(py_max2*100))+'% - '+str(round(py_max*100))+'%'
        if m eq 4 then $
          cgText, xrange_halo[0]+0.6, py_min+0.39, barStr, alignment=0.5, color=cgColor('gray')
          
        barStr = str(round(py_min*100))+'% - '+str(round(py_max2*100))+'% '
        if m eq 4 then $
          cgText, xrange_halo[0]+0.6, py_min+0.06, barStr, alignment=0.5, color=cgColor('gray')
        
        cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_frac,/xs,/ys,/ylog,yminor=0,$
          ytitle=textoidl("(dM/dt)_{"+accModes[m]+"} / (dM/dt)_{total}"),$
          xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),$
          /noerase,pos=pos[m-1]-[0,0.02,0,0.04]
      
        for i=0,n_tags(sP)-2 do begin
          if m ge n_tags(mbv.(i).(zInd)) then continue ; mode does not exist for this run
          
          ; gal hot
          xvals = alog10( 10.0^mbv.(i).(zind).(0).logMassBinCen / units.hubbleParam )
          yvals = smooth( mbv.(i).(zInd).(m).mode_ratio.gal_hot, sK, /nan )
          cgPlot,xvals,yvals,color=sP.(i).(zind).(0).colors[cInd],line=linesHC[0],/overplot
          
          ; gal cold
          xvals = alog10( 10.0^mbv.(i).(zind).(0).logMassBinCen / units.hubbleParam )
          yvals = smooth( mbv.(i).(zInd).(m).mode_ratio.gal_cold, sK, /nan )
          cgPlot,xvals,yvals,color=sP.(i).(zind).(0).colors[cInd],line=linesHC[1],/overplot
            
        endfor
        
        ; panel-specific legends
        if m eq 0 then $
          legend,textoidl('z_{ }=_{ }'+string(redshifts[zind],format='(f3.1)')),/top,/right
        if m eq 2 then $
          legend,textoidl(["T_{max} > T_{vir,acc}","T_{max} < T_{vir,acc}"]),$
            linestyle=linesHC,color=cgColor('black'),textcolor=cgColor('black'),/bottom,/right
        if m eq 3 then $
          legend,simNames,textcolors=simColors,box=0,/top,/left
          
      endfor ;m

    end_PS
  endforeach
       
  ; plots vs halo mass (all 2x2 redshifts together, one plot per accMode)
  ; ------------------------------------------------------------
  for j=0,n_tags(sP.(1).(0))-1 do begin
    pos = plot_pos(total=n_elements(zInds),/gap) ; plot positioning (3x2, 2x2, or 1x2 with gaps)
    
    plotStrAM = plotStr + '_am-'+accModes[j]
  
    ; plot (7) - net cold fractions, galaxy/halo
    start_PS, sP.(0).(0).(0).plotPath + 'coldFracNet4Redshifts.galaxy-halo.' + plotStrAM + '.eps', /extrabig
      
      foreach zInd,zInds,k do begin
        cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_frac2,xs=5,ys=5,$
          /noerase,pos=pos[k] ; suppress all axes, only establish coordinate system
      
        ;cgPolygon,[xrange_halo[0],xrange_halo[1],xrange_halo[1],xrange_halo[0],xrange_halo[0]],$
        ;  [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
        ;cgPolygon,[xrange_halo[0],xrange_halo[1],xrange_halo[1],xrange_halo[0],xrange_halo[0]],$
        ;  [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
          
        cgPolygon,[px_min,px_max,px_max,px_min,px_min],$
          [yrange_frac2[0],yrange_frac2[0],yrange_frac2[1],yrange_frac2[1],yrange_frac2[0]],$
          color=cgColor(px_cl),/fill
          
        cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_frac2,/xs,/ys,$
          ytitle=textoidl("(dM/dt)_{T_{max}<T_{vir,acc}} / (dM/dt)_{total,"+accModes[j]+" }"),$
          xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),/noerase,pos=pos[k]
      
        for i=0,n_tags(sP)-2 do begin
          xvals = alog10( 10.0^mbv.(i).(zInd).(j).logMassBinCen / units.hubbleParam )
          
          yvals = smooth(mbv.(i).(zInd).(j).coldfrac_net.gal,sK,/nan) < 1.0
          cgPlot,xvals,yvals,color=simColors[i],line=linesHC[0],/overplot ; gal
          
          yvals = smooth(mbv.(i).(zInd).(j).coldfrac_net.halo,sK,/nan) < 1.0
          cgPlot,xvals,yvals,color=simColors[i],line=linesHC[1],/overplot ; halo
        endfor
        
        ; add gadget line at z=2
        if k eq 1 and addGA eq 1 then begin
          gaColor = sP.sP2_gadget.(0).(0).colors[cInd]
          
          xvals = alog10( 10.0^mbv.mbv2_Gadget.redshift6_200.(j).logMassBinCen / units.hubbleParam )
          
          yvals = smooth(mbv.mbv2_Gadget.redshift6_200.(j).coldfrac_net.gal,sK,/nan) < 1.0
          cgPlot,xvals,yvals,color=gaColor,line=linesHC[0],/overplot ; gal
          
          yvals = smooth(mbv.mbv2_Gadget.redshift6_200.(j).coldfrac_net.halo,sK,/nan) < 1.0
          cgPlot,xvals,yvals,color=gaColor,line=linesHC[1],/overplot ; halo
          
          legend,['GADGET'],textcolors=[gaColor],box=0,/top,/left
        endif
        
        ; redshift legend
        legend,textoidl('z_{ }=_{ }'+string(redshifts[zInd],format='(f3.1)')),/top,/right
          
        ; panel-specific legends
        if k eq 2 then $
          legend,['galaxy','halo'],linestyle=linesHC,$
          color=cgColor('black'),textcolors=['black','black'],linesize=0.4,box=0,/top,/left
          
        if k eq 3 then $
          legend,simNames,textcolors=simColors,box=0,/top,/left
      endforeach

    end_PS
    
    ; plot (8) - specific net accretion rates, galaxy hot/cold
    start_PS, sP.(0).(0).(0).plotPath + 'netSpecificRate4Redshifts.galaxy.' + plotStrAM + '.eps', /extrabig
      
      foreach zInd,zInds,k do begin
        yrange2 = yrange_spec - 0.5*k
        
        cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange2,xs=5,ys=5,$
          /noerase,pos=pos[k] ; suppress all axes, only establish coordinate system
      
        cgPolygon,[px_min,px_max,px_max,px_min,px_min],$
          [yrange2[0],yrange2[0],yrange2[1],yrange2[1],yrange2[0]],$
          color=cgColor(px_cl),/fill
          
        cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange2,/xs,/ys,$
          ytitle=textoidl("dM_{gas,"+accModes[j]+"}/dt  M_{halo}^{-1} [_{ }log Gyr^{-1 }]"),$
          xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),/noerase,pos=pos[k]
      
        for i=0,n_tags(sP)-2 do begin
          xvals = alog10( 10.0^mbv.(i).(zInd).(j).logMassBinCen / units.hubbleParam )
          
          yvals = mbv.(i).(zInd).(j).specific_net.gal_hot
          yvals = alog10( smooth(yvals,sK,/nan) )
          cgPlot,xvals,yvals,color=simColors[i],line=linesHC[0],/overplot ; galaxy, hot
          
          yvals = mbv.(i).(zInd).(j).specific_net.gal_cold
          yvals = alog10( smooth(yvals,sK,/nan) )
          cgPlot,xvals,yvals,color=simColors[i],line=linesHC[1],/overplot ; galaxy, cold
        endfor
        
        ; redshift legend
        legend,textoidl('z_{ }=_{ }'+string(redshifts[zind],format='(f3.1)')),/top,/right
        
        ; panel-specific legends
        if k eq 0 then $
          legend,textoidl(['T_{max} > T_{vir,acc}','T_{max} < T_{vir,acc}']),$
          linestyle=linesHC,color=cgColor('black'),textcolors=replicate('black',2),box=0,/top,/left
          
        if k eq 3 then $
          legend,simNames,textcolors=simColors,box=0,/top,/left
      endforeach

    end_PS
    
  endfor ;accModes,j

  stop       
end

; plotMaxValsVsRedshift(): TODO

pro plotMaxValsVsRedshift
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  message,'revamp'
  ; config
  runs       = ['feedback','tracer'] ;,'gadget']
  redshifts  = [3.0,2.0,1.0,0.0]
  res        = 512
  accMode    = 'smooth' ; accretion mode: all, smooth, clumpy, stripped, recycled
  timeWindow = 500.0    ; consider accretion over this past time range (Myr)
                        ; 250.0 500.0 1000.0 "all" "tVir_tIGM" or "tVir_tIGM_bin"
  tVirInd    = 0        ; use Tmax/Tviracc=1 to separate hot vs. cold
  
  ; max histograms config
  massBinInd = 4        ; plot 11.0<logM<11.5 halos for ValMax histos
  tagNames   = ['Tmax','TmaxTviracc','EntMax','DensMax','MachMax'] ; plot each quantity
  pParams = { TmaxTviracc : {xrange:[-2.2,1.4], xtickv:[-2,-1,0,1],   xlabel : "log ( T_{max} / T_{vir,acc} )"} ,$
              Tmax        : {xrange:[3.5,7.5],  xtickv:[4,5,6,7],     xlabel : "log ( T_{max} )"}               ,$
              EntMax      : {xrange:[4.5,10.0], xtickv:[5,6,7,8,9],   xlabel : "log ( S_{max} ) [K cm^{2 }]"}   ,$
              DensMax     : {xrange:[-10,0],    xtickv:[-8,-6,-4,-2], xlabel : "log ( \rho_{max} )"}            ,$
              MachMax     : {xrange:[0,100],    xtickv:[10,30,50,80], xlabel : "M_{max}"} }
  
  ; plot config
  lines   = [0,1]       ; gal/both,gmem
  sK      = 3           ; smoothing kernel size  
  cInd    = 1           ; color index

  xrange_halo = [9.0,12.0]
  yrange_hist = [6e-4,2e-1]
  
  ; load
  tempFilename = '/n/home07/dnelson/temp_plotRatesMaxValsAt4Redshifts.sav'
  if ~file_test(tempFilename) then begin
  
  foreach run,runs,i do begin
    sP_z = {} & mbv_z = {} & bth_z = {}
    
    ; make for all the redshifts
    foreach redshift,redshifts,j do begin
      sP_loc = simParams(res=res,run=run,redshift=redshift)
      mbv_loc = haloMassBinValues(sP=sP_loc,accMode=accMode,timeWindow=timeWindow)

      ; load group cat for subgroup masses
      gc = loadGroupCat(sP=sP_loc,/skipIDs)
      gcMasses = codeMassToLogMsun(gc.subgroupMass)
      
      ; calculate specific net rates, galaxy by galaxy (yr -> Gyr)
      maxHist = n_elements(mbv_loc.galaxy.netRate.total.hot.num)
      gcMasses = gcMasses[0:maxHist-1]
        
      gal_hot   = reform(mbv_loc.galaxy.netRate.total.hot.tVirAcc[tVirInd,*]) / 10.0^gcMasses * 1e9
      gal_cold  = reform(mbv_loc.galaxy.netRate.total.cold.tVirAcc[tVirInd,*]) / 10.0^gcMasses * 1e9
      halo_hot  = reform(mbv_loc.halo.netRate.total.hot.tVirAcc[tVirInd,*]) / 10.0^gcMasses * 1e9
      halo_cold = reform(mbv_loc.halo.netRate.total.cold.tVirAcc[tVirInd,*]) / 10.0^gcMasses * 1e9

      specific_net = { gal_hot   : fltarr(n_elements(mbv_loc.logMassBinCen)) + !values.f_nan ,$
                       gal_cold  : fltarr(n_elements(mbv_loc.logMassBinCen)) + !values.f_nan ,$
                       halo_hot  : fltarr(n_elements(mbv_loc.logMassBinCen)) + !values.f_nan ,$
                       halo_cold : fltarr(n_elements(mbv_loc.logMassBinCen)) + !values.f_nan  }
      coldfrac_net = { gal  : fltarr(n_elements(mbv_loc.logMassBinCen)) + !values.f_nan ,$
                       halo : fltarr(n_elements(mbv_loc.logMassBinCen)) + !values.f_nan  }
      
      for k=0,n_elements(mbv_loc.logMassBinCen)-1 do begin
        ; restrict medians to primary only (num>0)
        w = where(gcMasses gt mbv_loc.logMassBins[k] and gcMasses le mbv_loc.logMassBins[k+1] and $
                  (mbv_loc.galaxy.coldFrac.total.num+mbv_loc.halo.coldFrac.total.num) gt 0,count)
                  
        if count eq 0 then continue
        if max(w) ge maxHist then message,'Error'
        
        specific_net.gal_hot[k]   = median( gal_hot[w] )
        specific_net.gal_cold[k]  = median( gal_cold[w] )
        specific_net.halo_hot[k]  = median( halo_hot[w] )
        specific_net.halo_cold[k] = median( halo_cold[w] )
        
        coldfrac_net.gal[k]  = median( gal_cold[w] / (gal_hot[w]+gal_cold[w]) )
        coldfrac_net.halo[k] = median( halo_cold[w] / (halo_hot[w]+halo_cold[w]) )
      endfor ;k
      
      mbv_loc = { galaxyMedian  : mbv_loc.galaxyMedian  ,$
                  haloMedian    : mbv_loc.haloMedian    ,$
                  logMassBinCen : mbv_loc.logMassBinCen ,$
                  specific_net  : specific_net          ,$
                  coldfrac_net  : coldfrac_net           }
      
      ; add to keeper array for this redshift
      sP_z  = mod_struct(sP_z, 'redshift'+str(j), sP_loc)
      mbv_z = mod_struct(mbv_z, 'redshift'+str(j), mbv_loc)
      ;bth_z = mod_struct(bth_z, 'redshift'+str(j), $
      ;  binValMaxHistos(sP=sP_z.(j),accMode=accMode,timeWindow=timeWindow))
    endforeach
    
    ; put this mode collection into mbv, once per run, and sP collection also
    mbv = mod_struct(mbv, 'mbv'+str(i), mbv_z)
    ;bth = mod_struct(bth, 'bth'+str(i), bth_z)
    sP  = mod_struct(sP, 'sP'+str(i), sP_z)
  endforeach
  
    save,mbv,sP,filename=tempFilename
    print,'Saved TEMP file: ['+tempFilename+']'
  endif else begin
    restore,tempFilename,/verbose
  endelse
  
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
  
  ; maximum val (temp,temp/tvir,ent,dens,mach) histos (allgal,gmem) (one halo mass bin)
  foreach tagName,tagNames do begin
    bVind   = ( where(tag_names(bth.(0).(0).allGal) eq strupcase(tagName)) )[0]
    pPind   = ( where(tag_names(pParams) eq strupcase(tagName)) )[0]
    
    xrange = pParams.(pPind).xrange
    xtickv = pParams.(pPind).xtickv
    
    start_PS, sP.(0).(0).plotPath + 'maxHistosRedshift_' + tagName + '.massBin' + $
              str(massBinInd) + '.' + plotStr + '.eps', /big
    
    for zind=0,3 do begin
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_hist,/xs,/ys,/ylog,yminor=0,$
        ytitle=textoidl("Gas Mass Fraction"),$
        xtitle=textoidl( pParams.(pPind).xlabel ),$
        /noerase,pos=pos[zind],xticks=n_elements(xtickv)-1,xtickv=xtickv
    
      for i=0,n_tags(sP)-2 do begin
      
        ; gal
        xvals = bth.(i).(zind).params.binLoc.(bVind)
        yvals = float( bth.(i).(zind).allGal.(bVind)[massBinInd,*] ) / $
                total( bth.(i).(zind).allGal.(bVind)[j,*] )
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[0],/overplot ; allgal
        
        ; gmem
        yvals = float( bth.(i).(zind).gmem.(bVind)[massBinInd,*] ) / $
                total( bth.(i).(zind).gmem.(bVind)[j,*] )
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
  endforeach ;tagName
   
  stop
end

; plotEvolIGMTemp(): plot evolving temperature of the IGM gas (outside all FOFs) vs redshift

pro plotEvolIGMTemp

  ; config
  sP_FB   = simParams(res=512,run='feedback')
  sP_noFB = simParams(res=512,run='tracer')
  
  igmEvol_FB   = evolIGMTemp(sP=sP_FB)
  igmEvol_noFB = evolIGMTemp(sP=sP_noFB)
  
  lines    = [0,2,1,3]    ; mean,median,min,max
  cInd     = 1            ; color index
  xrange_z = [-0.15,6.15] ; redshift
  yrange_T = [1.5,9.8]    ; temp (log K)
  yrange_P = [1e-7,1e-1]  ; pdf
  zDists   = [5.0,3.0,2.0,1.0,0.5,0.0] ; 6 redshifts to plot full distributions
  
  ; plot (1) - mean/median evolution of IGM gas temperature with redshift
  pos_single = [0.14,0.11,0.94,0.87] ; single plot with second x-axis above
  
  plotStr = sP_FB.savPrefix + str(sP_FB.res) + "_" + sP_noFB.savPrefix + str(sP_noFB.res)
  
  start_PS, sP_FB.plotPath + 'igmTemp.comp.'+plotStr+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_T,xs=9,/ys,$
      ytitle="IGM Gas Temperature [ log K ]",xtitle="Redshift",pos=pos_single
    
    universeage_axis, xrange_z, yrange_T
    
    ; with feedback
    cgPlot,igmEvol_FB.redshifts,igmEvol_FB.meanTemp,$
      color=sP_FB.colors[cInd],line=lines[0],/overplot
    cgPlot,igmEvol_FB.redshifts,igmEvol_FB.medianTemp,$
      color=sP_FB.colors[cInd],line=lines[1],/overplot
    cgPlot,igmEvol_FB.redshifts,igmEvol_FB.minTemp,$
      color=sP_FB.colors[cInd],line=lines[2],/overplot
    cgPlot,igmEvol_FB.redshifts,igmEvol_FB.maxTemp,$
      color=sP_FB.colors[cInd],line=lines[3],/overplot

    ; without feedback
    cgPlot,igmEvol_noFB.redshifts,igmEvol_noFB.meanTemp,$
      color=sP_noFB.colors[cInd],line=lines[0],/overplot
    cgPlot,igmEvol_noFB.redshifts,igmEvol_noFB.medianTemp,$
      color=sP_noFB.colors[cInd],line=lines[1],/overplot
    cgPlot,igmEvol_noFB.redshifts,igmEvol_noFB.minTemp,$
      color=sP_noFB.colors[cInd],line=lines[2],/overplot
    cgPlot,igmEvol_noFB.redshifts,igmEvol_noFB.maxTemp,$
      color=sP_noFB.colors[cInd],line=lines[3],/overplot
    
    ; legends
    legend,[sP_FB.simName,sP_noFB.simName],textcolors=[sP_FB.colors[1],sP_noFB.colors[1]],$
      pos=[xrange_z[0]+2.6,yrange_T[1]-0.2]
      
    legend,['mean','median','min','max'],linestyle=lines,/top,/right
    
  end_PS
  
  ; plot (2) - distributions
  start_PS, sP_FB.plotPath + 'igmTemp.compdist.'+plotStr+'.eps', xs=10*1.4, ys=6*1.4
  
    pos = plot_pos(total=n_elements(zDists),/gap)
    
    for i=0,n_elements(zDists)-1 do begin
     
      ; locate index
      curInd_FB = closest( igmEvol_FB.redshifts, zDists[i] )
      curInd_noFB = closest( igmEvol_FB.redshifts, zDists[i] )
    
      ; plot
      cgPlot,[0],[0],/nodata,xrange=yrange_T,yrange=yrange_P,xs=1,ys=1,/ylog,yminor=0,$
        ytitle="PDF",xtitle="IGM Gas Temperature [ log K ]",pos=pos[i],/noerase
  
      yy = igmEvol_FB.tempDist[curInd_FB,*]
      cgPlot,igmEvol_FB.tempBinCen,yy/total(yy),color=sP_FB.colors[cInd],/overplot
      
      yy = igmEvol_noFB.tempDist[curInd_noFB,*]
      cgPlot,igmEvol_noFB.tempBinCen,yy/total(yy),color=sP_noFB.colors[cInd],/overplot
      
      legend,["z="+string(zDists[i],format='(f3.1)')],/top,/left
  
      if i eq 2 then $
        legend,[sP_FB.simName,sP_noFB.simName],textcolors=[sP_FB.colors[1],sP_noFB.colors[1]],/top,/right
    endfor
  
  end_PS
  
  stop

end
