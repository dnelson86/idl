; plotVsHaloMass.pro
; gas accretion project - plots as a function of halo mass
; dnelson apr.2013

; plotPreBin():

pro plotPreBin

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect   = 'pri'
  redshift   = 2.0 ;redshifts  = [2.0,3.0]
  res        = 128
  
  timeWindows = list(250.0,500.0,1000.0,'all') ;list('all','tVir_tIGM','tVir_tIGM_bin') ; Myr
  accModes = list('all','smooth','clumpy','stripped') ;,'recycled'
  runs = list('gadget') ;list('tracer','gadget')
  
  ;foreach redshift,redshifts do begin
    foreach timeWindow,timeWindows do begin
      foreach run,runs do begin
        foreach accMode,accModes do begin
          print,run,res,timeWindow,accMode
          sP = simParams(res=res,run=run,redshift=redshift)
        
          binv = haloMassBinValues(sP=sP,sgSelect=sgSelect,timeWindow=timeWindow,accMode=accMode)
          th   = binTmaxHistos(sP=sP,sgSelect=sgSelect,timeWindow=timeWindow,accMode=accMode)
          h2   = binTmaxHisto2D(sP=sP,sgSelect=sgSelect,timeWindow=timeWindow,accMode=accMode)
        endforeach
      endforeach
    endforeach
  ;endforeach
  
  print,'Done.'
  
end

; markPlots():

pro markPlots

  sgSelect   = 'pri'
  accMode    = 'smooth'
  timeWindow = 1000.0 ; Myr
  tVirInd    = 0 ; plot Tmax/Tvir=1
  res        = 512
  redshift   = 2.0

  sK      = 3 ; smoothing kernel size  
  cInd    = 1 ; color index
  
  xrange = [10.0,12.0]
  yrange = [0.0,1.0] 
  
  ; load
  sPg = simParams(res=res,run='gadget',redshift=redshift)
  sPa = simParams(res=res,run='tracer',redshift=redshift) ; f=-1 use velocity tracers

  GA = haloMassBinValues(sP=sPg,sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow)
  AR = haloMassBinValues(sP=sPa,sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow)

  start_PS, 'plot_coldfrac.eps'
  
    !p.thick += 2
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,$
      ytitle="Cold Fraction",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
      ytickv=[0.0,0.2,0.4,0.6,0.8,1.0],yticks=5
    
    cgPlot,GA.logMassBinCen,smooth(GA.fracMedian.both_tViracc[tVirInd,*],sK,/nan),color='blue',line=0,/overplot
    cgPlot,AR.logMassBinCen,smooth(AR.fracMedian.both_tViracc[tVirInd,*],sK,/nan),color='red',line=0,/overplot  
  
    legend,['gadget','arepo'],textcolors=['blue','red'],box=0,/top,/left
  
  end_PS
  
  yrange = [0.01,20.0]
  
  start_PS, 'plot_accrates.eps', xs=10.0, ys=5.0
  
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y0,x1,y2] ,$ ;left
                [x1,y0,x2,y2]  ) ;right
                
    !p.thick += 2
   
    ; left: hot both (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),xtitle="",xticks=4,xtickv=[10.0,10.5,11.0,11.5,12.0],pos=pos[0]
    
    cgPlot,GA.logMassBinCen,smooth(GA.hotMedian.both_tviracc[tVirInd,*],sK,/nan),color='blue',line=0,/overplot
    cgPlot,AR.logMassBinCen,smooth(AR.hotMedian.both_tviracc[tVirInd,*],sK,/nan),color='red',line=0,/overplot
    
    legend,['gadget','arepo'],textcolors=['blue','red'],box=0,/top,/left
    
    ; right: cold both (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xticks=3,xtickv=[10.5,11.0,11.5,12.0],ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    cgPlot,GA.logMassBinCen,smooth(GA.coldMedian.both_tviracc[tVirInd,*],sK,/nan),color='blue',line=0,/overplot
    cgPlot,AR.logMassBinCen,smooth(AR.coldMedian.both_tviracc[tVirInd,*],sK,/nan),color='red',line=0,/overplot
    
    ; labels
    cgText,x1,0.03,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    cgText,mean([x0,x1]),y2+0.015,"Hot",alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,"Cold",alignment=0.5,color=cgColor('dark gray'),/normal
    
  end_PS
  
  stop
  
end

; plotByMethod(): plot the "cold mode fraction" and accretion rates vs halo mass in a few different ways
;                 using the mergerTreeSubset with accretionTimes and accretionModes

pro plotByMethod

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  runs       = ['feedback','feedback_noZ']
  sgSelect   = 'pri'
  accMode    = 'smooth' ; accretion mode: all, smooth, bclumpy, sclumpy
  timeWindow = 1000.0 ; consider accretion over this past time range (Myr)
                      ; 250.0 500.0 1000.0 "all" "tVir_tIGM" or "tVir_tIGM_bin"
  tVirInd    = 0 ; plot Tmax/Tvir=1
  res        = 256
  redshift   = 2.0
  
  lines   = [1,0,2] ; tvircur,tviracc,const5.5
  thicks  = [4,6,8] ; 128,256,512
  sK      = 3 ; smoothing kernel size  
  cInd    = 1 ; color index

  xrange = [10.0,12.0]
  yrange = [0.0,1.025]
  
  ; load
  foreach run,runs,i do begin
    sP  = mod_struct(sP, 'sP'+str(i), simParams(res=res,run=run,redshift=redshift))
    mbv = mod_struct(mbv, 'mbv'+str(i), $
      haloMassBinValues(sP=sP.(i),sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow))
  endforeach

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

  plotStr += str(res) + '_' + str(sP.(0).snap) + '_tw' + twStr + '_am-' + accMode
  
  ; --- cold fractions ---
  
  ; plot (1) - 2x2 of (both,gmem)x(const,tvir)
  start_PS, sP.(0).plotPath + 'coldFracByMethod.both-gmem.' + plotStr + '.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.10 & y1 = 0.50 & y2 = 0.90
    
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: both tvir
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4,pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).fracMedian.both_tViracc[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    ; top left axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=redshift))
    cgAxis,/xaxis,xrange=tvirrange,/xs
    
    ; ur: both const
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TcutVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).fracMedian.both_const[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
      
    ; top right axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=redshift))
    tvirrange[1] = round(10*tvirrange[1])/10.0
    cgAxis,/xaxis,xrange=tvirrange,/xs

    ; ll: gmem tvir
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).fracMedian.gmem_tViracc[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
      
    ; legend
    strings = textoidl("T_{max} / T_{vir,acc} < ")+string(mbv.(0).TvirVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(mbv.(0).TvirVals)),box=0,linesize=0.25,$
      position=[10.8,0.97],charsize=!p.charsize-0.27

    ; lr: gmem const
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).fracMedian.gmem_const[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    ; legend
    strings = textoidl("T_{max} < T_{c} = 10^{"+string(mbv.(0).TcutVals,format='(f3.1)')+"}")
    legend,strings,linestyle=indgen(n_elements(mbv.(0).TcutVals)),box=0,linesize=0.25,$
      position=[10.92,0.97],charsize=!p.charsize-0.27

    legend,simNames,textcolors=simColors,box=0,/bottom,/left,charsize=!p.charsize-0.2

    ; labels
    cgText,0.05,y1,"Cold Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.01,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    cgText,x1,y2+0.06,textoidl("T_{vir} [_{ }log K_{ }] (z="+str(fix(redshift))+")"),alignment=0.5,/normal
    
    cgText,x2+0.02,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.02,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x1,x2])+0.1,y2+0.06,"Variable Constant "+textoidl("T_c"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1])-0.1,y2+0.06,"Variable "+textoidl("T_{max} / T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; --- accretion rates ---
  
  yrange = [0.01,20.0]
  yrange_stars = 0.1*yrange
  
  ; plot (1) - both/gmem, variable*const
  start_PS, sP.(0).plotPath + 'accRateByMethod.both-gmem.const.' + plotStr + '.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot both
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TcutVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).hotMedian.both_const[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    ; legend
    legend,simNames,textcolors=simColors,box=0,/top,/left,charsize=!p.charsize-0.2
    
    ; ur: cold both
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TcutVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).coldMedian.both_const[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    ; ll: hot gmem
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TcutVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).hotMedian.gmem_const[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot

    ; lr: cold gmem
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TcutVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).coldMedian.gmem_const[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    ; legend
    strings = textoidl("T_{max} < T_{c} = 10^{"+string(mbv.(0).TcutVals,format='(f3.1)')+"}")
    legend,strings,linestyle=indgen(n_elements(mbv.(0).TcutVals)),$
      box=0,linesize=0.25,position=[10.05,13.0],charsize=!p.charsize-0.27

    ; labels
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,"Hot",alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,"Cold",alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; plot (2) - both/gmem, variable*tvir
  start_PS, sP.(0).plotPath + 'accRateByMethod.both-gmem.tvir.'+plotStr+'.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot both
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).hotMedian.both_tVirAcc[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    ; legend
    legend,simNames,textcolors=simColors,box=0,/top,/left,charsize=!p.charsize-0.2
    
    ; ur: cold both
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).coldMedian.both_tVirAcc[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    ; ll: hot gmem
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).hotMedian.gmem_tVirAcc[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    ; lr: cold gmem
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).coldMedian.gmem_tVirAcc[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot

    ; legend
    strings = textoidl("T_{max} / T_{vir,acc} < ")+string(mbv.(0).TvirVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(mbv.(0).TvirVals)),$
      box=0,linesize=0.25,position=[10.05,13.0],charsize=!p.charsize-0.27

    ; labels
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,"Hot",alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,"Cold",alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; --- total mass ---  
  yrange_gal = [2e7,2e10]
  yrange_halo = [2e7,2e10]
  
  if accMode eq 'all' then yrange_gal *= 3
  if accMode eq 'all' then yrange_halo *= 3
  
  if 0 then begin
  ; plot (1) - both/gmem, variable*const
  start_PS, sP.(0).plotPath + 'totalMassByMethod.both-gmem.const.'+plotStr+'.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot both
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_gal,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TcutVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).hotTotal.both_const[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    ; legend
    legend,simNames,textcolors=simColors,box=0,/top,/left,charsize=!p.charsize-0.2
    
    ; ur: cold both
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_gal,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TcutVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).coldTotal.both_const[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    ; ll: hot gmem
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_halo,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TcutVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).hotTotal.gmem_const[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot

    ; lr: cold gmem
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_halo,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TcutVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).coldTotal.gmem_const[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
 
    ; legend
    strings = textoidl("T_{max} < T_{c} = 10^{"+string(mbv.(0).TcutVals,format='(f3.1)')+"}")
    legend,strings,linestyle=indgen(n_elements(mbv.(0).TcutVals)),$
      box=0,linesize=0.25,position=[10.05,0.7*yrange_halo[1]],charsize=!p.charsize-0.27

    ; labels
    cgText,0.05,y1,"Total Accreted Mass "+textoidl("[_{ }h^{-1} M_{sun }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,"Hot",alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,"Cold",alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; plot (2) - both/gmem, variable*tvir
  start_PS, sP.(0).plotPath + 'totalMassByMethod.both-gmem.tvir.'+plotStr+'.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot both
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_gal,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).hotTotal.both_tVirAcc[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    ; legend
    legend,simNames,textcolors=simColors,box=0,/top,/left,charsize=!p.charsize-0.2
    
    ; ur: cold both
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_gal,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).coldTotal.both_tVirAcc[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
    
    ; ll: hot gmem
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_halo,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).hotTotal.gmem_tVirAcc[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot
      
    ; lr: cold gmem
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_halo,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_elements(mbv.(i).TvirVals)-1 do $
        cgPlot,mbv.(i).logMassBinCen,smooth(mbv.(i).coldTotal.gmem_tVirAcc[j,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=j,/overplot

    ; legend
    strings = textoidl("T_{max} / T_{vir,acc} < ")+string(mbv.(0).TvirVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(mbv.(0).TvirVals)),$
      box=0,linesize=0.25,position=[10.05,0.7*yrange_halo[1]],charsize=!p.charsize-0.27

    ; labels
    cgText,0.05,y1,"Total Accreted Mass "+textoidl("[_{ }h^{-1} M_{sun }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,"Hot",alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,"Cold",alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  endif ;0
  
  stop
end

; plotByMode(): plot cold fractions and accretion rates separated out into the four modes and also
;               a view with decomposed galaxy=gal+stars (one resolution only)

pro plotByMode

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect   = 'pri'
  redshift   = 2.0
  res        = 256
  timeWindow = 1000.0
  runs       = ['feedback','feedback_noZ']
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
    
    ; make for all the accretion modes
    foreach accMode,accModes,j do begin
      if accMode eq 'recycled' and sP.(i).gfmWinds eq 0 then continue ; skip recycled for nonWinds
      mbv_mode = mod_struct(mbv_mode, 'res'+str(j), $
        haloMassBinValues(sP=sP.(i),sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow))
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

  ; plot (1) - 2x2 of (both,gmem)x(const,tvir)
  start_PS, sP.(0).plotPath + 'coldFracByMode.both-gmem.'+plotStr+'.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.10 & y1 = 0.50 & y2 = 0.90
    
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: gal tvir
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4,pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).fracMedian.both_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; top left axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=redshift))
    cgAxis,/xaxis,xrange=tvirrange,/xs
    
    ; ur: gal const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).fracMedian.both_const[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; top right axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=redshift))
    tvirrange[1] = round(10*tvirrange[1])/10.0
    cgAxis,/xaxis,xrange=tvirrange,/xs

    ; ll: gmem tvir
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).fracMedian.gmem_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; legend
    legend,accModes,linestyle=lines,box=0,linesize=0.4,/bottom,/left,charsize=!p.charsize-0.2
    
    ; lr: gmem const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).fracMedian.gmem_const[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; legend
    legend,simNames,textcolors=simColors,box=0,/top,/right,charsize=!p.charsize-0.2
      
    ; labels
    cgText,0.05,y1,"Cold Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.01,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    cgText,x1,y2+0.06,textoidl("T_{vir} [_{ }log K_{ }] (z="+str(fix(redshift))+")"),alignment=0.5,/normal
    
    cgText,x2+0.02,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.02,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x1,x2])+0.05,y2+0.06,textoidl("T_{max} < T_c = " + conIndName),$
      alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1])-0.05,y2+0.06,textoidl("T_{max} < " + virIndName + " T_{vir,acc}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; --- accretion rates ---

  yrange = [0.01,20.0]
  yrange_stars = 0.1*yrange

  ; plot (1) - accretion rate galaxy/halo atmosphere variable*tvir
  start_PS, sP.(0).plotPath + 'accRateByMode.both-gmem.tvir.'+plotStr+'.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot both
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).hotMedian.both_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; legend
    legend,simNames,textcolors=simColors,box=0,/top,/left,charsize=!p.charsize-0.2
    
    ; ur: cold both
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).coldMedian.both_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; ll: hot gmem
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).hotMedian.gmem_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; lr: cold gmem
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).coldMedian.gmem_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot

    
    ; legend
    legend,accModes,linestyle=lines,box=0,linesize=0.4,/top,/left,charsize=!p.charsize-0.2
      
    ; labels
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,textoidl("T_{max} > " + virIndName + " T_{vir,acc}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,textoidl("T_{max} < " + virIndName + " T_{vir,acc}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  yrange = [0.0,1.0]
  
  ; plot (2) - fraction, accretion rate galaxy/halo atmosphere variable*tviracc
  start_PS, sP.(0).plotPath + 'accRateFracByMode.both-gmem.tvir.'+plotStr+'.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot gal
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=1,n_tags(mbv.(i))-1 do $ ; take ratio to first
        cgPlot,mbv.(i).(j).logMassBinCen,$
        smooth(mbv.(i).(j).hotMedian.both_tVirAcc[virInd,*]/mbv.(i).(0).hotMedian.both_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; legend
    legend,simNames,textcolors=simColors,box=0,/top,/left,charsize=!p.charsize-0.2
    
    ; ur: cold gal (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=1,n_tags(mbv.(i))-1 do $ ; take ratio to first
        cgPlot,mbv.(i).(j).logMassBinCen,$
        smooth(mbv.(i).(j).coldMedian.both_tVirAcc[virInd,*]/mbv.(i).(0).coldMedian.both_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; ll: hot gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=1,n_tags(mbv.(i))-1 do $ ; take ratio to first
        cgPlot,mbv.(i).(j).logMassBinCen,$
        smooth(mbv.(i).(j).hotMedian.gmem_tVirAcc[virInd,*]/mbv.(i).(0).hotMedian.gmem_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; lr: cold gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=1,n_tags(mbv.(i))-1 do $ ; take ratio to first
        cgPlot,mbv.(i).(j).logMassBinCen,$
        smooth(mbv.(i).(j).coldMedian.gmem_tVirAcc[virInd,*]/mbv.(i).(0).coldMedian.gmem_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; legend
    legend,accModes[1:*],linestyle=lines[1:*],box=0,linesize=0.4,position=[11.2,0.95],charsize=!p.charsize-0.2
      
    ; labels
    cgText,0.05,y1,"Fraction of Gas Accretion Rate",alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
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
  sgSelect = 'pri'
  accMode  = 'smooth'
  redshift = 2.0
  res      = 512
  runs     = ['gadget','tracer','feedback']
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
    
    ; make for all the timeWindows
    foreach timeWindow,tws,j do $
      mbv_tw = mod_struct(mbv_tw, 'tw'+str(j), $
        haloMassBinValues(sP=sP.(i),sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow))
    
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

  ; plot (1) - 2x2 of (both,gmem)x(const,tvir)
  start_PS, sP.(0).plotPath + 'coldFracByTW.both-gmem.'+plotStr+'.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.10 & y1 = 0.50 & y2 = 0.90
    
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: gal tvir
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4,pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).fracMedian.both_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot

    ; top left axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=redshift))
    cgAxis,/xaxis,xrange=tvirrange,/xs
    
    ; ur: gal const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).fracMedian.both_const[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot

    ; top right axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=redshift))
    tvirrange[1] = round(10*tvirrange[1])/10.0
    cgAxis,/xaxis,xrange=tvirrange,/xs

    ; legend
    legend,twNames,linestyle=lines,box=0,linesize=0.4,position=[10.1,0.4],charsize=!p.charsize-0.27

    ; ll: gmem tvir
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).fracMedian.gmem_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot

    ; lr: gmem const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
      
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).fracMedian.gmem_const[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot

    legend,simNames,textcolors=simColors,box=0,/bottom,/left,charsize=!p.charsize-0.2

    ; labels
    cgText,0.05,y1,"Cold Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.01,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    cgText,x1,y2+0.06,textoidl("T_{vir} [_{ }log K_{ }] (z="+str(fix(redshift))+")"),alignment=0.5,/normal
    
    cgText,x2+0.02,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.02,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x1,x2])+0.05,y2+0.06,textoidl("T_{max} < T_c = " + conIndName),$
      alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1])-0.05,y2+0.06,textoidl("T_{max} < " + virIndName + " T_{vir,acc}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; --- accretion rates ---
  
  yrange = [0.01,20.0]
  yrange_stars = 0.1*yrange

  ; plot (2) - accretion rate both/gmem variable*tvir
  start_PS, sP.(0).plotPath + 'accRateByTW.both-gmem.'+plotStr+'.eps', /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot both (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).hotMedian.both_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; legend
    legend,simNames,textcolors=simColors,box=0,/top,/left,charsize=!p.charsize-0.2
    
    ; ur: cold both (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).coldMedian.both_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; legend
    if accMode eq 'smooth' then $
      legend,twNames,linestyle=lines,box=0,linesize=0.4,/top,/left,charsize=!p.charsize-0.27
    if accMode eq 'all' then $
      legend,twNames,linestyle=lines,box=0,linesize=0.4,/bottom,/right,charsize=!p.charsize-0.27
    
    ; ll: hot gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).hotMedian.gmem_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; lr: cold gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).coldMedian.gmem_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
      
    ; labels
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
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
  sgSelect    = 'pri'
  accMode     = 'smooth'
  redshift    = 3.0
  timeWindow  = 1000.0 ; Myr
  resolutions = [512,256,128]
  runs        = ['gadget','tracer','feedback']

  ; plot config
  lines  = [0,2,3] ; tvircur,tviracc,const5.5 or 128,512,256
  thicks = [4,6,8] ; 128,256,512  
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
    
    ; make for all the resolutions
    foreach res,resolutions,j do begin
      sP_res = simParams(res=res,run=run,redshift=redshift)
      mbv_res = mod_struct(mbv_res, 'res'+str(j), $
        haloMassBinValues(sP=sP_res,sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow))
    endforeach
    
    ; put this res collection into mbv, once per run
    mbv = mod_struct(mbv, 'mbv'+str(i), mbv_res)
  endforeach
  
  if n_elements(resolutions) ne n_elements(sP.(0).colors) then message,'Going to be a problem.'
  
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
  
  ; plot (1) - 3x2 panel split by definition type, both/gmem (resolution lines)
  start_PS, sP.(0).plotPath + 'coldFracByRes.3x2.'+plotStr+'.eps',/big
    
    x0 = 0.13 & x1 = 0.40 & x2 = 0.67 & x3 = 0.94
    y0 = 0.13 & y1 = 0.53 & y2 = 0.93
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; uc
                [x2,y1,x3,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1] ,$ ; lc
                [x2,y0,x3,y1] )  ; lr 

    ; ul - gal const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4,pos=pos[0]

    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).fracMedian.both_const[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
     
    ; uc - gal tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),yticks=4,pos=pos[1]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).fracMedian.both_tVirCur[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
  
    ; ur - gal tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[2]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).fracMedian.both_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; ll - gmem const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/noerase,$
      ytitle="",xtitle="",ytickv=[0.0,0.2,0.4,0.6,0.8,1.0],yticks=5,pos=pos[3]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).fracMedian.gmem_const[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    legend,simNames,textcolors=simColors,box=0,/bottom,/left,charsize=!p.charsize-0.2
    
    ; lc - gmem tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[4]
      
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).fracMedian.gmem_tVirCur[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
  
    ; lr - gmem tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[5]
      
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).fracMedian.gmem_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    legend,textoidl(str(resolutions)+'^3'),linestyle=lines,box=0,linesize=0.25,$
      /bottom,/left,charsize=!p.charsize-0.27
    
    ; labels    
    cgText,0.05,y1,"Cold Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,mean([x0,x3]),0.02,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x3+0.02,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x3+0.02,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal

    cgText,mean([x2,x3]),y2+0.03,textoidl("T_{max} < " + virIndName + " T_{vir,acc}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.03,textoidl("T_{max} < " + virIndName + " T_{vir,cur}"),$
      alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1]),y2+0.03,textoidl("T_{max} < T_{c} = 10^{" + conIndName + "} K"),$
      alignment=0.5,color=cgColor('dark gray'),/normal

   end_PS
  
  ; --- accretion rates ---
  
  yrange = [0.01,20.0]
  yrange_stars = 0.1*yrange
  
  ; plot (2) - 3x2 panel split by definition type, both hot/cold (resolution lines)
  start_PS, sP.(0).plotPath + 'accRateByRes.both.3x2.'+plotStr+'.eps',/big
    
    x0 = 0.13 & x1 = 0.40 & x2 = 0.67 & x3 = 0.94
    y0 = 0.13 & y1 = 0.53 & y2 = 0.93
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; uc
                [x2,y1,x3,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1] ,$ ; lc
                [x2,y0,x3,y1] )  ; lr 

    ; ul - both cold const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]

    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).coldMedian.both_const[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot

    legend,simNames,textcolors=simColors,box=0,/top,/left,charsize=!p.charsize-0.2

    ; uc - both cold tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),yticks=4,pos=pos[1]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).coldMedian.both_tVirCur[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
     
    ; ur - both cold tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[2]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).coldMedian.both_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    if accMode eq 'smooth' then $
      legend,textoidl(str(resolutions)+'^3'),linestyle=lines,box=0,linesize=0.34,position=[10.0,14.0],charsize=!p.charsize-0.27
    if accMode eq 'all' then $
      legend,textoidl(str(resolutions)+'^3'),linestyle=lines,box=0,linesize=0.34,/bottom,/right,charsize=!p.charsize-0.27

    ; ll - both hot const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",pos=pos[3]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).hotMedian.both_const[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot

    ; draw some resolution Nx element lines
    cgPlot,[elemMasses[0],elemMasses[0]],yrangeTicks,line=lines[0],color=cgColor('black'),thick=!p.thick+1,/overplot ; 512
    cgPlot,[elemMasses[1],elemMasses[1]],yrangeTicks2,line=lines[1],color=cgColor('black'),thick=!p.thick+1,/overplot ; 256
 
    ; lc - both hot tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[4]
      
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).hotMedian.both_tVirCur[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
 
    cgPlot,[elemMasses[0],elemMasses[0]],yrangeTicks,line=lines[0],color=cgColor('black'),thick=!p.thick+1,/overplot ; 512
    cgPlot,[elemMasses[1],elemMasses[1]],yrangeTicks2,line=lines[1],color=cgColor('black'),thick=!p.thick+1,/overplot ; 256
 
    ; lr - both hot tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[5]
      
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).hotMedian.both_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    cgPlot,[elemMasses[0],elemMasses[0]],yrangeTicks,line=lines[0],color=cgColor('black'),thick=!p.thick+1,/overplot ; 512
    cgPlot,[elemMasses[1],elemMasses[1]],yrangeTicks2,line=lines[1],color=cgColor('black'),thick=!p.thick+1,/overplot ; 256
 
    ; labels    
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
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

  ; plot (3) - both hot/cold/total (resolution lines)
  start_PS, sP.(0).plotPath + 'accRateByRes.both.3x3.'+plotStr+'.eps',xs=9,ys=9
    
    x0 = 0.13 & x1 = 0.40 & x2 = 0.67 & x3 = 0.94
    y0 = 0.13 & y1 = 0.40 & y2 = 0.67 & y3 = 0.94
    pos = list( [x0,y2,x1,y3] ,$ ; ul
                [x1,y2,x2,y3] ,$ ; uc
                [x2,y2,x3,y3] ,$ ; ur
                [x0,y1,x1,y2] ,$ ; ml
                [x1,y1,x2,y2] ,$ ; mc
                [x2,y1,x3,y2] ,$ ; mr
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1] ,$ ; lc
                [x2,y0,x3,y1] )  ; lr 

    ; ul - both total const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]

    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).totalHCMedian.both,sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot

    ; uc - both total tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),yticks=4,pos=pos[1]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).totalHCMedian.both,sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; ur - both total tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[2]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).totalHCMedian.both,sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    ; ml - both cold const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[3]

    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).coldMedian.both_const[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot

    ; mc - both cold tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),yticks=4,pos=pos[4]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).coldMedian.both_tVirCur[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
 
    ; mr - both cold tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[5]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).coldMedian.both_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    if accMode eq 'smooth' then $
      legend,textoidl(str(resolutions)+'^3'),linestyle=lines,box=0,linesize=0.25,/top,/left,charsize=!p.charsize-0.27
    if accMode eq 'all' then $
      legend,textoidl(str(resolutions)+'^3'),linestyle=lines,box=0,linesize=0.25,/bottom,/right,charsize=!p.charsize-0.27

    ; ll - both hot const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",pos=pos[6]
    
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).hotMedian.both_const[conInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot

    legend,simNames,textcolors=simColors,box=0,/top,/left,charsize=!p.charsize-0.2
    
    cgPlot,[elemMasses[0],elemMasses[0]],yrangeTicks,line=lines[0],color=cgColor('black'),thick=!p.thick+1,/overplot ; 512
    cgPlot,[elemMasses[1],elemMasses[1]],yrangeTicks2,line=lines[1],color=cgColor('black'),thick=!p.thick+1,/overplot ; 256
    
    ; lc - both hot tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[7]
      
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).hotMedian.both_tVirCur[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
 
    cgPlot,[elemMasses[0],elemMasses[0]],yrangeTicks,line=lines[0],color=cgColor('black'),thick=!p.thick+1,/overplot ; 512
    cgPlot,[elemMasses[1],elemMasses[1]],yrangeTicks2,line=lines[1],color=cgColor('black'),thick=!p.thick+1,/overplot ; 256
 
    ; lr - both hot tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[8]
      
    for i=0,n_tags(sP)-1 do $
      for j=0,n_tags(mbv.(i))-1 do $
        cgPlot,mbv.(i).(j).logMassBinCen,smooth(mbv.(i).(j).hotMedian.both_tVirAcc[virInd,*],sK,/nan),$
        color=sP.(i).colors[cInd],line=lines[j],/overplot
    
    cgPlot,[elemMasses[0],elemMasses[0]],yrangeTicks,line=lines[0],color=cgColor('black'),thick=!p.thick+1,/overplot ; 512
    cgPlot,[elemMasses[1],elemMasses[1]],yrangeTicks2,line=lines[1],color=cgColor('black'),thick=!p.thick+1,/overplot ; 256
    
    ; labels    
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
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

; plotHotColdMassRatio(): plot the mass ratio of hot to cold material in the "halo atmosphere"

pro plotHotColdMassRatio

  compile_opt idl2, hidden, strictarr, strictarrsubs

  foreach res,[128,256,512] do begin
  
    ; config
    redshift = 3.0
    sPa = simParams(res=res,run='feedback',redshift=redshift)
    sPg = simParams(res=res,run='gadget',redshift=redshift)

    haloID = 304 ;z2.304 z2.301 z2.130 z2.64
    gcID = getMatchedIDs(sPa=sPa,sPg=sPg,haloID=haloID)   
    
    ; load arepo
    gcA = loadGroupCat(sP=sPa,/readIDs)
    galcatA = galaxyCat(sP=sPa)
    parIDsA = galCatRepParentIDs(galcat=galcatA)
    
    u     = loadSnapshotSubset(sP=sPa,partType='gas',field='u')
    nelec = loadSnapshotSubset(sP=sPa,partType='gas',field='ne')
    temp  = convertUtoTemp(u,nelec)
    u     = !NULL
    nelec = !NULL
    
    mass = loadSnapshotSubset(sP=sPa,partType='gas',field='mass')
    ids  = loadSnapshotSubset(sP=sPa,partType='gas',field='ids')
    
    ; find gas elements for single halo in galcat
    w = where(parIDsA.gmem eq gcID.a,count)
    gmem_ids = galcatA.groupMemIDs[w]
    
    calcMatch,ids,gmem_ids,ind_ids,ind_gmem_ids,count=countMatch
    if countMatch ne count then message,'Error: Failed to find all ids.'
    
    gmem_temp = temp[ind_ids]
    
    ; find gas elements for single halo using primary subgroup and min radius
    sh_ids = gcPIDList(gc=gcA,select='pri',valGCids=[gcID.a],partType='gas')
    
      ; remove galaxy (<0.15 rvir)
      w = where(parIDsA.gal eq gcID.a,count)
      sh_ids_gal = galcatA.galaxyIDs[w]
      
      sh_ids = removeIntersectionFromB(sh_ids_gal,sh_ids)
    
    calcMatch,ids,sh_ids,ind_ids_sh,ind_sh_ids,count=countMatchSh
    if countMatchSh ne n_elements(sh_ids) then message,'Error: Failed to find all sh ids.'
    
    sh_temp = temp[ind_ids_sh]
    
    ; calculate fraction for single halo
    tvir = codeMassToVirTemp(gcA.subgroupMass[gcID.a],sP=sPa)
    
    foreach tVirRatio,sPa.tVirVals do begin
    
      w_cold = where(gmem_temp/tvir le tVirRatio,count_cold)
      w_hot  = where(gmem_temp/tvir gt tVirRatio,count_hot)
      
      if (count_cold+count_hot) ne countMatch then message,'Error: Failed to classify all.'
      
      frac = float(count_cold)/(count_hot+count_cold)
      frac_mass = total(mass[ind_ids[w_cold]]) / (total(mass[ind_ids[w_hot]])+total(mass[ind_ids[w_cold]]))
      print,res,tVirRatio,frac,frac_mass
      
      w_cold = where(sh_temp/tvir le tVirRatio,count_cold)
      w_hot  = where(sh_temp/tvir gt tVirRatio,count_hot)
      
      if (count_cold+count_hot) ne countMatchSh then message,'Error: Failed to classify sh all.'
      
      frac = float(count_cold)/(count_hot+count_cold)
      frac_mass = total(mass[ind_ids_sh[w_cold]]) / (total(mass[ind_ids_sh[w_hot]])+total(mass[ind_ids_sh[w_cold]]))
      print,res,tVirRatio,frac,frac_mass,'sh'
    
    endforeach
    
  endforeach

end
