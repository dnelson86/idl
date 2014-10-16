; plotVsHaloMass.pro
; gas accretion project - plots as a function of halo mass
; dnelson jan.2014

; plotPreBin():

pro plotPreBin, sP=sP, redshifts=redshifts

  compile_opt idl2, hidden, strictarr, strictarrsubs
  if n_elements(redshifts) eq 0 then message,'Error'
  
  ; config
  ;redshift   = sP.redshift ;2.0
  ;redshifts  = [2.0,3.0]
  res        = sP.res ;256
  run        = sP.run

  runs        = [sP.run] ;['feedback','gadget'] ;tracer
  resolutions = [sP.res]
  timeWindows = list(250.0) ;,'all') ;list('all','tVir_tIGM','tVir_tIGM_bin') ; Myr
  accModes    = list('smooth') ;list('all','smooth','clumpy','stripped','recycled')
  
  foreach redshift,redshifts do begin
    foreach timeWindow,timeWindows do begin
      foreach run,runs do begin
        foreach res,resolutions do begin
          foreach accMode,accModes do begin
            if sP.gfmWinds eq 0 and accMode eq 'recycled' then continue ; skip
            
            print,run,res,redshift,timeWindow,accMode
            sP = simParams(res=res,run=run,redshift=redshift,hind=sP.hInd)
        
            binv = haloMassBinValues(sP=sP,timeWindow=timeWindow,accMode=accMode)
            th   = binValMaxHistos(sP=sP,timeWindow=timeWindow,accMode=accMode)
          endforeach
        endforeach
      endforeach
    endforeach
  endforeach
  
  print,'Done.'
  
end

; plotByMethod(): plot the "cold mode fraction" and accretion rates vs halo mass in a few different ways
;                 with accretionTimes and accretionModes

pro plotByMethod

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  runs       = ['gadget','tracer','feedback'] ;['feedback','feedback_noZ','feedback_noFB']
  accMode    = 'smooth' ; accretion mode: all, smooth, bclumpy, sclumpy
  timeWindow = 500.0 ; consider accretion over this past time range (Myr)
                      ; 250.0 500.0 1000.0 "all" "tVir_tIGM" or "tVir_tIGM_bin"
  res        = 512
  redshift   = 2.0
  
  lines   = [1,0,2] ; tvircur,tviracc,const5.5
  sK      = 3 ; smoothing kernel size
  cInd    = 1 ; color index
  zoomSym = ['open circle','open square','open diamond'] ;cutVals/virVals

  xrange = [10.0,12.0]
  yrange = [0.0,1.025]
  yrange_rate = [0.001,20.0]
  
  ; add any zoom run single points?
  ;sPz = mod_struct( sPz, 'sPz0', simParams(run='zoom_20Mpc',res=9,hInd=0,redshift=redshift) )
  ;sPz = mod_struct( sPz, 'sPz1', simParams(run='zoom_20Mpc',res=10,hInd=0,redshift=redshift) )
  ;sPz = mod_struct( sPz, 'sPz2', simParams(run='zoom_20Mpc',res=11,hInd=0,redshift=redshift) )

  ; load
  foreach run,runs,i do begin
    sP  = mod_struct(sP, 'sP'+str(i), simParams(res=res,run=run,redshift=redshift))
    mbv = mod_struct(mbv, 'mbv'+str(i), haloMassBinValues(sP=sP.(i),accMode=accMode,timeWindow=timeWindow))
  endforeach
  
  for i=0,n_tags(sPz)-1 do $
    mbvZ = mod_struct( mbvZ, 'mbvZ'+str(i), haloMassBinValues(sP=sPz.(i),accMode=accMode,timeWindow=timeWindow))
  
  ; strings
  if isnumeric(timeWindow) then twStr = str(fix(timeWindow))
  if ~isnumeric(timeWindow) then twStr = "-"+timeWindow
  
  plotStr   = ''
  simNames  = []
  simColors = []
  
  foreach run,runs,i do begin
    plotStr   = plotStr + sP.(i).plotPrefix + '.'
    simNames  = [simNames, sP.(i).simName]
    simColors = [simColors, sP.(i).colors[cInd]]
  endforeach
  for i=0,n_tags(sPz)-1 do begin
    plotStr   = plotStr + sPz.(i).plotPrefix + '.'
    simNames  = [simNames, sPz.(i).simName]
    simColors = [simColors, sPz.(i).colors[cInd]]
  endfor

  plotStr += str(res) + '_' + str(sP.(0).snap) + '_tw' + twStr + $
             '_am-' + accMode + '_model-' + str(sP.(0).accRateModel)
  
  ; --- cold fractions ---
  
  ; plot (1) - 2x2 of (galaxy,halo)x(const,tvir)
  start_PS, sP.(0).plotPath + 'coldFracByMethod.galaxy-halo.' + plotStr + '.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.10 & y1 = 0.50 & y2 = 0.90
    
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: allgal tvir
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4,pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).galaxyMedian.coldFrac.total.tVirAcc[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
        
    for i=0,n_tags(sPz)-1 do $
      for j=0,n_elements(mbvZ.(i).TvirVals)-1 do $
        cgPlot,mbvZ.(i).logMassBinCen,mbvZ.(i).galaxy.coldFrac.total.tVirAcc[j],$
        color=sPz.(i).colors[cInd],psym=zoomSym[j],/overplot
    
    ; top left axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=redshift))
    cgAxis,/xaxis,xrange=tvirrange,/xs
    
    ; symbol legend for zooms
    if n_tags(sPz) gt 0 then legend,string(mbvZ.(0).TvirVals,format='(f4.1)'),psym=zoomSym,/top,/left
    
    ; ur: allgal const
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TcutVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).galaxyMedian.coldFrac.total.tConst[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
        
    for i=0,n_tags(sPz)-1 do $
      for j=0,n_elements(mbvZ.(i).TvirVals)-1 do $
        cgPlot,mbvZ.(i).logMassBinCen,mbvZ.(i).galaxy.coldFrac.total.tConst[j],$
        color=sPz.(i).colors[cInd],psym=zoomSym[j],/overplot
      
    ; top right axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=redshift))
    tvirrange[1] = round(10*tvirrange[1])/10.0
    cgAxis,/xaxis,xrange=tvirrange,/xs
    
    ; symbol legend for zooms
    if n_tags(sPz) gt 0 then legend,string(mbvZ.(0).TcutVals,format='(f4.1)'),psym=zoomSym,/top,/right

    ; ll: gmem tvir
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).haloMedian.coldFrac.total.tViracc[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
      
    for i=0,n_tags(sPz)-1 do $
      for j=0,n_elements(mbvZ.(i).TvirVals)-1 do $
        cgPlot,mbvZ.(i).logMassBinCen,mbvZ.(i).halo.coldFrac.total.tVirAcc[j],$
        color=sPz.(i).colors[cInd],psym=zoomSym[j],/overplot
      
    ; legend
    strings = textoidl("T_{max} / T_{vir,acc} < ")+string(mbv.(0).TvirVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(mbv.(0).TvirVals)),box=0,linesize=0.25,/top,/left

    ; lr: gmem const
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).haloMedian.coldFrac.total.tConst[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    for i=0,n_tags(sPz)-1 do $
      for j=0,n_elements(mbvZ.(i).TvirVals)-1 do $
        cgPlot,mbvZ.(i).logMassBinCen,mbvZ.(i).halo.coldFrac.total.tConst[j],$
        color=sPz.(i).colors[cInd],psym=zoomSym[j],/overplot
    
    ; legend
    strings = textoidl("T_{max} < T_{c} = 10^{"+string(mbv.(0).TcutVals,format='(f3.1)')+"}")
    legend,strings,linestyle=indgen(n_elements(mbv.(0).TcutVals)),box=0,linesize=0.25,/top,/right

    legend,simNames,textcolors=simColors,box=0,/bottom,/left

    ; labels
    cgText,0.05,y1,"Cold Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.01,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    cgText,x1,y2+0.06,textoidl("T_{vir} [_{ }log K_{ }] (z="+str(fix(redshift))+")"),alignment=0.5,/normal
    
    cgText,x2+0.02,mean([y0,y1]),"Halo",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.02,mean([y1,y2]),"Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x1,x2])+0.1,y2+0.06,"Constant "+textoidl("T_c"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1])-0.1,y2+0.06,"Variable "+textoidl("T_{max} / T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; --- accretion rates ---
  if accMode eq 'all' then yrange_rate *= [2,5]
  
  labels  = ['accRate','outRate','netRate'] ; each a separate plot
  ytitles = ['Accretion','Outflow','Net Inflow']
  
  foreach label,labels,k do begin
  
    rateInd = where( tag_names(mbv.(0).galaxyMedian) eq strupcase(label) )
  
  ; plot (1) - galaxy/halo, variable*const
  if 0 then begin
  start_PS, sP.(0).plotPath + label + 'ByMethod.galaxy-halo.const.' + plotStr + '.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot galaxy
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_rate,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TcutVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).galaxyMedian.(rateInd).total.hot.tConst[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    for i=0,n_tags(sPz)-1 do $
      for j=0,n_elements(mbvZ.(i).TvirVals)-1 do $
        cgPlot,mbvZ.(i).logMassBinCen,mbvZ.(i).galaxy.(rateInd).total.hot.tConst[j],$
        color=sPz.(i).colors[cInd],psym=zoomSym[j],/overplot
    
    ; symbol legend for zooms
    if n_tags(sPz) gt 0 then legend,string(mbvZ.(0).TcutVals,format='(f4.1)'),psym=zoomSym,/top,/left
    
    ; ur: cold galaxy
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_rate,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TcutVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).galaxyMedian.(rateInd).total.cold.tConst[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    for i=0,n_tags(sPz)-1 do $
      for j=0,n_elements(mbvZ.(i).TvirVals)-1 do $
        cgPlot,mbvZ.(i).logMassBinCen,mbvZ.(i).galaxy.(rateInd).total.cold.tConst[j],$
        color=sPz.(i).colors[cInd],psym=zoomSym[j],/overplot
    
    ; legend
    legend,simNames,textcolors=simColors,box=0,/top,/left
    
    ; ll: hot halo
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_rate,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TcutVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).haloMedian.(rateInd).total.hot.tConst[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot

    for i=0,n_tags(sPz)-1 do $
      for j=0,n_elements(mbvZ.(i).TvirVals)-1 do $
        cgPlot,mbvZ.(i).logMassBinCen,mbvZ.(i).halo.(rateInd).total.hot.tConst[j],$
        color=sPz.(i).colors[cInd],psym=zoomSym[j],/overplot
        
    ; lr: cold halo
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_rate,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TcutVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).haloMedian.(rateInd).total.cold.tConst[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    for i=0,n_tags(sPz)-1 do $
      for j=0,n_elements(mbvZ.(i).TvirVals)-1 do $
        cgPlot,mbvZ.(i).logMassBinCen,mbvZ.(i).halo.(rateInd).total.cold.tConst[j],$
        color=sPz.(i).colors[cInd],psym=zoomSym[j],/overplot
    
    ; legend
    strings = textoidl("T_{max} < T_{c} = 10^{"+string(mbv.(0).TcutVals,format='(f3.1)')+"}")
    legend,strings,linestyle=indgen(n_elements(mbv.(0).TcutVals)),linesize=0.25,/top,/left

    ; labels
    cgText,0.05,y1,"Gas "+ytitles[k]+" Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Halo",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,"Hot",alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,"Cold",alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  endif ;0
  
  ; plot (2) - galaxy/halo, variable*tvir
  start_PS, sP.(0).plotPath + label + 'ByMethod.galaxy-halo.tvir.'+plotStr+'.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot galaxy
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_rate,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).galaxyMedian.(rateInd).total.hot.tVirAcc[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    for i=0,n_tags(sPz)-1 do $
      for j=0,n_elements(mbvZ.(i).TvirVals)-1 do $
        cgPlot,mbvZ.(i).logMassBinCen,mbvZ.(i).galaxy.(rateInd).total.hot.tVirAcc[j],$
        color=sPz.(i).colors[cInd],psym=zoomSym[j],/overplot
    
    ; symbol legend for zooms
    if n_tags(sPz) gt 0 then legend,string(mbvZ.(0).TvirVals,format='(f4.1)'),psym=zoomSym,/top,/left
    
    ; ur: cold galaxy
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_rate,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).galaxyMedian.(rateInd).total.cold.tVirAcc[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    for i=0,n_tags(sPz)-1 do $
      for j=0,n_elements(mbvZ.(i).TvirVals)-1 do $
        cgPlot,mbvZ.(i).logMassBinCen,mbvZ.(i).galaxy.(rateInd).total.cold.tVirAcc[j],$
        color=sPz.(i).colors[cInd],psym=zoomSym[j],/overplot
     
    ; ll: hot halo
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_rate,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).haloMedian.(rateInd).total.hot.tVirAcc[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    for i=0,n_tags(sPz)-1 do $
      for j=0,n_elements(mbvZ.(i).TvirVals)-1 do $
        cgPlot,mbvZ.(i).logMassBinCen,mbvZ.(i).halo.(rateInd).total.hot.tVirAcc[j],$
        color=sPz.(i).colors[cInd],psym=zoomSym[j],/overplot
    
    ; legend
    strings = textoidl("T_{max} / T_{vir,acc} < ")+string(mbv.(0).TvirVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(mbv.(0).TvirVals)),linesize=0.25,/bottom,/right
    
    ; lr: cold halo
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_rate,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).haloMedian.(rateInd).total.cold.tVirAcc[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot

    for i=0,n_tags(sPz)-1 do $
      for j=0,n_elements(mbvZ.(i).TvirVals)-1 do $
        cgPlot,mbvZ.(i).logMassBinCen,mbvZ.(i).halo.(rateInd).total.cold.tVirAcc[j],$
        color=sPz.(i).colors[cInd],psym=zoomSym[j],/overplot
        
    ; legend
    legend,simNames,textcolors=simColors,box=0,/bottom,/right
       
    ; labels
    cgText,0.05,y1,"Gas "+ytitles[k]+" Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Halo",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,"Hot",alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,"Cold",alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  endforeach ; labels
  
end

; plotByMode(): plot cold fractions and accretion rates separated out into the four modes

pro plotByMode

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  redshift   = 0.0
  res        = 512
  timeWindow = 500.0
  runs       = ['gadget','tracer','feedback'] ;['feedback','feedback_noZ','feedback_noFB']
  accModes   = ['all','smooth','clumpy','stripped','recycled']
  
  ; plot config
  sK     = 3 ; smoothing kernel size
  cInd   = 1 ; color index
  virInd = 0 ; 1.0 Tvir
  conInd = 1 ; 5.5 const
  lines  = [0,1,2,3,4] ; for each accMode
  
  xrange = [10.0,12.0]
  yrange = [0.0,1.025]
  
  if accModes[0] ne 'all' then message,'Error: Fraction plot not going to work.'

  ; load
  foreach run,runs,i do begin
    sP  = mod_struct(sP, 'sP'+str(i), simParams(res=res,run=run,redshift=redshift))
    mbv_mode = {}
    
    ; make for all the accretion modes
    foreach accMode,accModes,j do begin
      if accMode eq 'recycled' and sP.(i).gfmWinds eq 0 then continue ; skip recycled for nonWinds     
      mbv_mode = mod_struct(mbv_mode, 'mode'+str(j), $
        haloMassBinValues(sP=sP.(i),accMode=accMode,timeWindow=timeWindow))
    endforeach
    
    ; put this mode collection into mbv, once per run
    mbv = mod_struct(mbv, 'mbv'+str(i), mbv_mode)
  endforeach
  
  ; strings
  if isnumeric(timeWindow) then twStr = str(fix(timeWindow))
  if ~isnumeric(timeWindow) then twStr = "-"+timeWindow
  
  plotStr   = ''
  simNames  = []
  simColors = []
  twNames   = []
  
  foreach run,runs,i do begin
    plotStr   = plotStr + sP.(i).plotPrefix + '.'
    simNames  = [simNames, sP.(i).simName]
    simColors = [simColors, sP.(i).colors[cInd]]
  endforeach

  plotStr += str(res) + '_' + str(sP.(0).snap) + '_tw' + twStr
  
  virIndName = string(mbv.(0).(0).TvirVals[virInd],format='(f3.1)')
  conIndName = string(mbv.(0).(0).TcutVals[conInd],format='(f3.1)')  
  
  ; --- cold fractions ---

  ; plot (1) - 2x2 of (galaxy,halo)x(const,tvir)
  start_PS, sP.(0).plotPath + 'coldFracByMode.galaxy-halo.'+plotStr+'.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.10 & y1 = 0.50 & y2 = 0.90
    
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: galaxy tvir
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4,pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).galaxyMedian.coldFrac.total.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; top left axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=redshift))
    cgAxis,/xaxis,xrange=tvirrange,/xs
    
    ; ur: galaxy const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).galaxyMedian.coldFrac.total.tConst[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; top right axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=redshift))
    tvirrange[1] = round(10*tvirrange[1])/10.0
    cgAxis,/xaxis,xrange=tvirrange,/xs

    ; ll: halo tvir
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).haloMedian.coldFrac.total.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; legend
    legend,accModes,linestyle=lines,box=0,linesize=0.4,/bottom,/left
    
    ; lr: halo const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).haloMedian.coldFrac.total.tConst[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; legend
    legend,simNames,textcolors=simColors,box=0,/top,/right
      
    ; labels
    cgText,0.05,y1,"Cold Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.01,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    cgText,x1,y2+0.06,textoidl("T_{vir} [_{ }log K_{ }] (z="+str(fix(redshift))+")"),alignment=0.5,/normal
    
    cgText,x2+0.02,mean([y0,y1]),"Halo",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.02,mean([y1,y2]),"Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x1,x2])+0.05,y2+0.06,textoidl("T_{max} < T_c = " + conIndName),$
      alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1])-0.05,y2+0.06,textoidl("T_{max} < " + virIndName + " T_{vir,acc}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; --- accretion rates ---

  yrange = [0.01,20.0]

  ; plot (1) - net accretion rate galaxy/halo atmosphere variable*tvir
  start_PS, sP.(0).plotPath + 'netRateByMode.galaxy-halo.tvir.'+plotStr+'.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot galaxy
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).galaxyMedian.netRate.total.hot.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; legend
    legend,simNames,textcolors=simColors,box=0,/top,/left
    
    ; ur: cold galaxy
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).galaxyMedian.netRate.total.cold.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; ll: hot halo
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).haloMedian.netRate.total.hot.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; lr: cold halo
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).haloMedian.netRate.total.cold.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot

    ; legend
    legend,accModes,linestyle=lines,box=0,linesize=0.4,/top,/left
      
    ; labels
    cgText,0.05,y1,"Net Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Halo",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,textoidl("T_{max} > " + virIndName + " T_{vir,acc}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,textoidl("T_{max} < " + virIndName + " T_{vir,acc}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  yrange = [0.0,1.0]
  
  ; plot (2) - fraction, net accretion rate galaxy/halo atmosphere variable*tviracc
  start_PS, sP.(0).plotPath + 'netRateFracByMode.allgal-gmem.tvir.'+plotStr+'.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot galaxy
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=1,n_tags(mbv.(i))-1 do $ ; take ratio to first
        cgPlot,mbv.(i).(j).logMassBinCen,$
        smooth(mbv.(i).(j).galaxyMedian.netRate.total.hot.tVirAcc[virInd,*]/mbv.(i).(0).galaxyMedian.netRate.total.hot.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; legend
    legend,simNames,textcolors=simColors,box=0,/top,/left
    
    ; ur: cold galaxy (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=1,n_tags(mbv.(i))-1 do $ ; take ratio to first
        cgPlot,mbv.(i).(j).logMassBinCen,$
        smooth(mbv.(i).(j).galaxyMedian.netRate.total.cold.tVirAcc[virInd,*]/mbv.(i).(0).galaxyMedian.netRate.total.cold.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; ll: hot halo (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=1,n_tags(mbv.(i))-1 do $ ; take ratio to first
        cgPlot,mbv.(i).(j).logMassBinCen,$
        smooth(mbv.(i).(j).haloMedian.netRate.total.hot.tVirAcc[virInd,*]/mbv.(i).(0).haloMedian.netRate.total.hot.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; lr: cold halo (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=1,n_tags(mbv.(i))-1 do $ ; take ratio to first
        cgPlot,mbv.(i).(j).logMassBinCen,$
        smooth(mbv.(i).(j).haloMedian.netRate.total.cold.tVirAcc[virInd,*]/mbv.(i).(0).haloMedian.netRate.total.cold.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; legend
    legend,accModes[1:*],linestyle=lines[1:*],box=0,linesize=0.4,position=[11.2,0.95]
      
    ; labels
    cgText,0.05,y1,"Fraction of Net Accretion Rate",alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Halo",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,textoidl("T_{max} > " + virIndName + " T_{vir,acc}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,textoidl("T_{max} < " + virIndName + " T_{vir,acc}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
end

; plotByTimeWindow(): plot cold fractions and accretion rates for multiple time windows (one resolution only)

pro plotByTimeWindow

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  accMode  = 'all'
  redshift = 2.0
  res      = 512
  runs     = ['tracer','feedback'] ;['feedback','feedback_noZ','feedback_noFB']
  tws      = list(250.0,500.0,1000.0,'all')
  
  ; plot config
  sK     = 3 ; smoothing kernel size
  cInd   = 1 ; color index
  virInd = 0 ; 1.0 Tvir
  conInd = 1 ; 5.5 const
  lines  = [1,2,3,0] ; for tws
  
  xrange = [10.0,12.0]
  yrange = [0.0,1.025]
  
  ; load
  foreach run,runs,i do begin
    sP  = mod_struct(sP, 'sP'+str(i), simParams(res=res,run=run,redshift=redshift))
    mbv_tw = {}
    
    ; make for all the timeWindows
    foreach timeWindow,tws,j do $
      mbv_tw = mod_struct(mbv_tw, 'tw'+str(j), $
        haloMassBinValues(sP=sP.(i),accMode=accMode,timeWindow=timeWindow))
    
    ; put this timeWindow collection into mbv, once per run
    mbv = mod_struct(mbv, 'mbv'+str(i), mbv_tw)
  endforeach
  
  ; strings
  plotStr   = ''
  simNames  = []
  simColors = []
  twNames   = []
  
  foreach run,runs,i do begin
    plotStr   = plotStr + sP.(i).plotPrefix + '.'
    simNames  = [simNames, sP.(i).simName]
    simColors = [simColors, sP.(i).colors[cInd]]
  endforeach

  plotStr += str(res) + '_' + str(sP.(0).snap) + '_am-' + accMode

  foreach timeWindow,tws do begin
    if isnumeric(timeWindow) then begin
      twNames = [twNames,str(fix(timeWindow))+' Myr']
    endif else begin
      if timeWindow eq 'all' then twNames = [twNames,'Integrated']
    endelse
  endforeach
  
  virIndName = string(mbv.(0).(0).TvirVals[virInd],format='(f3.1)')
  conIndName = string(mbv.(0).(0).TcutVals[conInd],format='(f3.1)')
  
  ; --- cold fractions ---

  ; plot (1) - 2x2 of (galaxy,halo)x(const,tvir)
  start_PS, sP.(0).plotPath + 'coldFracByTW.galaxy-halo.'+plotStr+'.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.10 & y1 = 0.50 & y2 = 0.90
    
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: galaxy tvir
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4,pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).galaxyMedian.coldFrac.total.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot

    ; top left axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=redshift))
    cgAxis,/xaxis,xrange=tvirrange,/xs
    
    ; ur: galaxy const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).galaxyMedian.coldFrac.total.const[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot

    ; top right axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=redshift))
    tvirrange[1] = round(10*tvirrange[1])/10.0
    cgAxis,/xaxis,xrange=tvirrange,/xs

    ; legend
    legend,twNames,linestyle=lines,box=0,linesize=0.4,position=[10.1,0.4]

    ; ll: halo tvir
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).haloMedian.coldFrac.total.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot

    ; lr: halo const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
      
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).haloMedian.coldFrac.total.tConst[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot

    legend,simNames,textcolors=simColors,box=0,/bottom,/left

    ; labels
    cgText,0.05,y1,"Cold Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.01,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    cgText,x1,y2+0.06,textoidl("T_{vir} [_{ }log K_{ }] (z="+str(fix(redshift))+")"),alignment=0.5,/normal
    
    cgText,x2+0.02,mean([y0,y1]),"Halo",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.02,mean([y1,y2]),"Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x1,x2])+0.05,y2+0.06,textoidl("T_{max} < T_c = " + conIndName),$
      alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1])-0.05,y2+0.06,textoidl("T_{max} < " + virIndName + " T_{vir,acc}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; --- accretion rates ---
  
  yrange = [0.01,20.0]

  ; plot (2) - net accretion rate galaxy/halo variable*tvir
  start_PS, sP.(0).plotPath + 'netRateByTW.galaxy-halo.'+plotStr+'.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot galaxy (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).galaxyMedian.netRate.total.hot.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; legend
    legend,simNames,textcolors=simColors,box=0,/top,/left
    
    ; ur: cold galaxy (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).galaxyMedian.netRate.total.cold.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; legend
    if accMode eq 'smooth' then $
      legend,twNames,linestyle=lines,box=0,linesize=0.4,/top,/left
    if accMode eq 'all' then $
      legend,twNames,linestyle=lines,box=0,linesize=0.4,/bottom,/right
    
    ; ll: hot halo (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).haloMedian.netRate.total.hot.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; lr: cold halo (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).haloMedian.netRate.total.cold.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
      
    ; labels
    cgText,0.05,y1,"Net Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Halo",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,textoidl("T_{max} > "+virIndName+" T_{vir,acc}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,textoidl("T_{max} < "+virIndName+" T_{vir,acc}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
end

; plotByRes(): plot cold fractions and accretion rates for multiple resolutions

pro plotByRes

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  accMode     = 'smooth'
  redshift    = 0.0
  timeWindow  = 500.0 ; Myr
  resolutions = [512,256,128]
  runs        = ['gadget','tracer','feedback']

  ; plot config
  lines  = [0,2,3] ; 512,256,128
  sK     = 3 ; smoothing kernel size
  cInd   = 1 ; color index
  virInd = 0 ; 1.0 Tvir
  conInd = 1 ; 5.5 const
  nElem  = 5000 ; draw this many gas elements as our "resolved" lines
  
  xrange = [10.0,12.0]
  yrange = [0.0,1.025]
  
  ; load
  foreach run,runs,i do begin
    sP  = mod_struct(sP, 'sP'+str(i), simParams(res=resolutions[0],run=run,redshift=redshift))
    mbv_res = {}
    
    ; make for all the resolutions
    foreach res,resolutions,j do begin
      sP_res = simParams(res=res,run=run,redshift=redshift)
      mbv_res = mod_struct(mbv_res, 'res'+str(j), $
        haloMassBinValues(sP=sP_res,accMode=accMode,timeWindow=timeWindow))
    endforeach
    
    ; put this res collection into mbv, once per run
    mbv = mod_struct(mbv, 'mbv'+str(i), mbv_res)
  endforeach
  
  if n_elements(resolutions) gt n_elements(sP.(0).colors) then message,'Going to be a problem.'
  
  ; strings
  if isnumeric(timeWindow) then twStr = str(fix(timeWindow))
  if ~isnumeric(timeWindow) then twStr = "-"+timeWindow
  
  plotStr   = ''
  simNames  = []
  simColors = []
  twNames   = []
  
  foreach run,runs,i do begin
    plotStr   = plotStr + sP.(i).plotPrefix + '.'
    simNames  = [simNames, sP.(i).simName]
    simColors = [simColors, sP.(i).colors[cInd]]
  endforeach

  plotStr += str(sP.(0).snap) + '_tw' + twStr + '_am-' + accMode
  
  virIndName = string(mbv.(0).(0).TvirVals[virInd],format='(f3.1)')
  conIndName = string(mbv.(0).(0).TcutVals[conInd],format='(f3.1)')
  
  ; calculate "resolved" halo masses
  yrangeTicks = [2.0,4.0]
  yrangeTicks2 = [0.3,0.7]
  elemMass512 = sP.(0).targetGasMass
  elemMasses = codeMassToLogMsun([elemMass512,elemMass512*8,elemMass512*8^2] * nElem / units.omega_b)
  
  ; --- cold fractions ---
  
  ; plot (1) - 3x2 panel split by definition type, galaxy/halo (resolution lines)
  start_PS, sP.(0).plotPath + 'coldFracByRes.3x2.'+plotStr+'.eps',/big
    
    x0 = 0.13 & x1 = 0.40 & x2 = 0.67 & x3 = 0.94
    y0 = 0.13 & y1 = 0.53 & y2 = 0.93
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; uc
                [x2,y1,x3,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1] ,$ ; lc
                [x2,y0,x3,y1] )  ; lr 

    ; ul - galaxy const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4,pos=pos[0]

    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).galaxyMedian.coldFrac.total.tConst[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
     
    ; uc - galaxy tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),yticks=4,pos=pos[1]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).galaxyMedian.coldFrac.total.tVirCur[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
  
    ; ur - galaxy tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[2]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).galaxyMedian.coldFrac.total.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; ll - halo const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/noerase,$
      ytitle="",xtitle="",ytickv=[0.0,0.2,0.4,0.6,0.8,1.0],yticks=5,pos=pos[3]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).haloMedian.coldFrac.total.tConst[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    legend,simNames,textcolors=simColors,box=0,/bottom,/left
    
    ; lc - halo tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[4]
      
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).haloMedian.coldFrac.total.tVirCur[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
  
    ; lr - halo tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[5]
      
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).haloMedian.coldFrac.total.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    legend,textoidl(str(resolutions)+'^3'),linestyle=lines,box=0,linesize=0.25,/bottom,/left
    
    ; labels    
    cgText,0.05,y1,"Cold Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,mean([x0,x3]),0.02,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x3+0.02,mean([y0,y1]),"Halo",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x3+0.02,mean([y1,y2]),"Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal

    cgText,mean([x2,x3]),y2+0.03,textoidl("T_{max} < " + virIndName + " T_{vir,acc}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.03,textoidl("T_{max} < " + virIndName + " T_{vir,cur}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1]),y2+0.03,textoidl("T_{max} < T_{c} = 10^{" + conIndName + "} K"),$
      alignment=0.5,color=cgColor('dark gray'),/normal

   end_PS
  
  ; --- accretion rates ---
  
  yrange = [0.01,20.0]
  
  ; plot (2) - 3x2 panel split by definition type, galaxy hot/cold (resolution lines)
  start_PS, sP.(0).plotPath + 'netRateByRes.galaxy.3x2.'+plotStr+'.eps',/big
    
    x0 = 0.13 & x1 = 0.40 & x2 = 0.67 & x3 = 0.94
    y0 = 0.13 & y1 = 0.53 & y2 = 0.93
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; uc
                [x2,y1,x3,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1] ,$ ; lc
                [x2,y0,x3,y1] )  ; lr 

    ; ul - galaxy cold const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]

    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).galaxyMedian.netRate.total.cold.tConst[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot

    legend,simNames,textcolors=simColors,box=0,/top,/left

    ; uc - galaxy cold tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),yticks=4,pos=pos[1]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).galaxyMedian.netRate.total.cold.tVirCur[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
     
    ; ur - galaxy cold tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[2]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).galaxyMedian.netRate.total.cold.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    if accMode eq 'smooth' then $
      legend,textoidl(str(resolutions)+'^3'),linestyle=lines,box=0,linesize=0.34,position=[10.0,14.0]
    if accMode eq 'all' then $
      legend,textoidl(str(resolutions)+'^3'),linestyle=lines,box=0,linesize=0.34,/bottom,/right

    ; ll - galaxy hot const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",pos=pos[3]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).galaxyMedian.netRate.total.hot.tConst[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot

    ; draw some resolution Nx element lines
    cgPlot,[elemMasses[0],elemMasses[0]],yrangeTicks,line=lines[0],color=cgColor('black'),thick=!p.thick+1,/overplot ; 512
    cgPlot,[elemMasses[1],elemMasses[1]],yrangeTicks2,line=lines[1],color=cgColor('black'),thick=!p.thick+1,/overplot ; 256
 
    ; lc - galaxy hot tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[4]
      
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).galaxyMedian.netRate.total.hot.tVirCur[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
 
    cgPlot,[elemMasses[0],elemMasses[0]],yrangeTicks,line=lines[0],color=cgColor('black'),thick=!p.thick+1,/overplot ; 512
    cgPlot,[elemMasses[1],elemMasses[1]],yrangeTicks2,line=lines[1],color=cgColor('black'),thick=!p.thick+1,/overplot ; 256
 
    ; lr - galaxy hot tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[5]
      
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).galaxyMedian.netRate.total.hot.tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    cgPlot,[elemMasses[0],elemMasses[0]],yrangeTicks,line=lines[0],color=cgColor('black'),thick=!p.thick+1,/overplot ; 512
    cgPlot,[elemMasses[1],elemMasses[1]],yrangeTicks2,line=lines[1],color=cgColor('black'),thick=!p.thick+1,/overplot ; 256
 
    ; labels    
    cgText,0.05,y1,"Net Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,mean([x0,x3]),0.02,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x3+0.02,mean([y0,y1]),"Galaxy (Hot)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x3+0.02,mean([y1,y2]),"Galaxy (Cold)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal

    cgText,mean([x2,x3]),y2+0.03,textoidl("T_{max} < " + virIndName + " T_{vir,acc}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.03,textoidl("T_{max} < " + virIndName + " T_{vir,cur}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1]),y2+0.03,textoidl("T_{max} < T_{c} = 10^{" + conIndName + "} K"),$
      alignment=0.5,color=cgColor('dark gray'),/normal

   end_PS
  
end
