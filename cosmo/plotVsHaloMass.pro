; plotVsHaloMass.pro
; gas accretion project - plots as a function of halo mass
; dnelson sep.2012

; plotPreBin():

pro plotPreBin

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect   = 'pri'
  redshift   = 2.0
  res        = 256
  
  timeWindows = list(250.0,500.0,1000.0,'all') ; Myr
  
  foreach timeWindow,timeWindows do begin
  
    ; load
    sPg = simParams(res=res,run='gadget',redshift=redshift)
    sPa = simParams(res=res,run='tracer',redshift=redshift)
  
    GA_all      = haloMassBinValues(sP=sPg,sgSelect=sgSelect,timeWindow=timeWindow,accMode='all')
    AR_all      = haloMassBinValues(sP=sPa,sgSelect=sgSelect,timeWindow=timeWindow,accMode='all')
    GA_smooth   = haloMassBinValues(sP=sPg,sgSelect=sgSelect,timeWindow=timeWindow,accMode='smooth')
    AR_smooth   = haloMassBinValues(sP=sPa,sgSelect=sgSelect,timeWindow=timeWindow,accMode='smooth')
    GA_clumpy   = haloMassBinValues(sP=sPg,sgSelect=sgSelect,timeWindow=timeWindow,accMode='clumpy')
    AR_clumpy   = haloMassBinValues(sP=sPa,sgSelect=sgSelect,timeWindow=timeWindow,accMode='clumpy') 
    GA_stripped = haloMassBinValues(sP=sPg,sgSelect=sgSelect,timeWindow=timeWindow,accMode='stripped')
    AR_stripped = haloMassBinValues(sP=sPa,sgSelect=sgSelect,timeWindow=timeWindow,accMode='stripped')
    
  endforeach
  
  print,'Done.'
  
end

; plotByMethod(): plot the "cold mode fraction" and accretion rates vs halo mass in a few different ways
;                 using the mergerTreeSubset with accretionTimes and accretionModes

pro plotByMethod

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect   = 'pri'
  accMode    = 'smooth' ; accretion mode: all, smooth, bclumpy, sclumpy
  timeWindow = 1000.0 ; consider accretion over this past time range (Myr)
  tVirInd    = 0 ; plot Tmax/Tvir=1
  res        = 512
  redshift   = 2.0
  
  lines   = [1,0,2] ; tvircur,tviracc,const5.5
  thicks  = [4,6,8] ; 128,256,512
  sK      = 3 ; smoothing kernel size  
  cInd    = 1 ; color index

  xrange = [10.0,12.0]
  yrange = [0.0,1.025]  
  
  ; load
  sPg = simParams(res=res,run='gadget',redshift=redshift)
  sPa = simParams(res=res,run='tracer',redshift=redshift) ; f=-1 use velocity tracers

  GA = haloMassBinValues(sP=sPg,sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow)
  AR = haloMassBinValues(sP=sPa,sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow)
  
  ; --- cold fractions ---
  
  ; plot (1) - 2x2 of (both,gmem)x(const,tvir)
  pName = sPg.plotPath + 'coldFracByMethod.both-gmem.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps'
  start_PS, pName, /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.10 & y1 = 0.50 & y2 = 0.90
    
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: both tvir (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4,pos=pos[0]
    
    for j=0,n_elements(GA.TvirVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.fracMedian.both_tViracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(AR.TvirVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.fracMedian.both_tViracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; top left axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    cgAxis,/xaxis,xrange=tvirrange,/xs
    
    ; ur: both const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    ; gadget and arepo gal
    for j=0,n_elements(GA.TcutVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.fracMedian.both_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(GA.TcutVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.fracMedian.both_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
      
    ; top right axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    tvirrange[1] = round(10*tvirrange[1])/10.0
    cgAxis,/xaxis,xrange=tvirrange,/xs

    ; ll: gmem tvir (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for j=0,n_elements(GA.TvirVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.fracMedian.gmem_tViracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(GA.TvirVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.fracMedian.gmem_tViracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
      
    ; legend
    strings = textoidl("T_{max} / T_{vir,acc} < ")+string(GA.TvirVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(GA.TvirVals)),box=0,linesize=0.25,$
      position=[10.8,0.97],charsize=!p.charsize-0.27

    ; lr: gmem const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for j=0,n_elements(GA.TcutVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.fracMedian.gmem_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(GA.TcutVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.fracMedian.gmem_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot

    ; legend
    strings = textoidl("T_{max} < T_{c} = ")+string(GA.TcutVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(GA.TcutVals)),box=0,linesize=0.25,$
      position=[10.95,0.97],charsize=!p.charsize-0.27

    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/bottom,/left

    ; labels
    cgText,0.05,y1,"Cold Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.01,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    cgText,x1,y2+0.06,textoidl("T_{vir} [_{ }log K_{ }] (z=2)"),alignment=0.5,/normal
    
    cgText,x2+0.02,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.02,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x1,x2])+0.1,y2+0.06,"Variable Constant "+textoidl("T_c"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1])-0.1,y2+0.06,"Variable "+textoidl("T_{max} / T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; plot (1) - 2x2 of (gal,stars)x(const,tvir)
  pName = sPg.plotPath + 'coldFracByMethod.gal-stars.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps'
  start_PS, pName, /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.10 & y1 = 0.50 & y2 = 0.90
    
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: gal tvir (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4,pos=pos[0]
    
    for j=0,n_elements(GA.TvirVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.fracMedian.gal_tViracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(AR.TvirVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.fracMedian.gal_tViracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; top left axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    cgAxis,/xaxis,xrange=tvirrange,/xs
    
    ; ur: gal const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    ; gadget and arepo gal
    for j=0,n_elements(GA.TcutVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.fracMedian.gal_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(GA.TcutVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.fracMedian.gal_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
      
    ; top right axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    tvirrange[1] = round(10*tvirrange[1])/10.0
    cgAxis,/xaxis,xrange=tvirrange,/xs

    ; ll: stars tvir (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for j=0,n_elements(GA.TvirVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.fracMedian.stars_tViracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(GA.TvirVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.fracMedian.stars_tViracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
      
    ; legend
    strings = textoidl("T_{max} / T_{vir} < ")+string(GA.TvirVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(GA.TvirVals)),box=0,linesize=0.4,$
      position=[10.95,0.95],charsize=!p.charsize-0.2

    ; lr: stars const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for j=0,n_elements(GA.TcutVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.fracMedian.stars_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(GA.TcutVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.fracMedian.stars_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot

    ; legend
    strings = textoidl("T_{max} < T_{c} = ")+string(GA.TcutVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(GA.TcutVals)),box=0,linesize=0.4,$
      position=[10.95,0.95],charsize=!p.charsize-0.2

    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/bottom,/left

    ; labels
    cgText,0.05,y1,"Cold Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.01,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    cgText,x1,y2+0.06,textoidl("T_{vir} [_{ }log K_{ }] (z=2)"),alignment=0.5,/normal
    
    cgText,x2+0.02,mean([y0,y1]),"Galaxy (Stars Only)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.02,mean([y1,y2]),"Galaxy (Gas Only)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x1,x2])+0.05,y2+0.06,"Variable Constant "+textoidl("T_c"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1])-0.05,y2+0.06,"Variable "+textoidl("T_{max} / T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS  
  
  ; --- accretion rates ---
  
  yrange = [0.01,20.0]
  yrange_stars = 0.1*yrange
  
  ; plot (1) - both/gmem, variable*const
  pName = sPg.plotPath + 'accRateByMethod.both-gmem.const.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps'
  start_PS, pName, /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot both (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for j=0,n_elements(GA.TcutVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.hotMedian.both_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(AR.TcutVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.hotMedian.both_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; ur: cold both (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for j=0,n_elements(GA.TcutVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.coldMedian.both_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(AR.TcutVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.coldMedian.both_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; ll: hot gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for j=0,n_elements(GA.TcutVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.hotMedian.gmem_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(AR.TcutVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.hotMedian.gmem_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot

    ; lr: cold gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for j=0,n_elements(GA.TcutVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.coldMedian.gmem_const[j,*],sK,/nan),color=sPg.colorsG[1],line=j,/overplot
    for j=0,n_elements(AR.TcutVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.coldMedian.gmem_const[j,*],sK,/nan),color=sPa.colorsA[1],line=j,/overplot

    ; legend
    strings = textoidl("T_{max} / T_{c} = ")+string(GA.TcutVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(GA.TcutVals)),box=0,linesize=0.25,position=[10.05,13.0]

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
  pName = sPg.plotPath + 'accRateByMethod.both-gmem.tvir.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps'
  start_PS, pName, /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot both (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for j=0,n_elements(GA.TvirVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.hotMedian.both_tviracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(AR.TvirVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.hotMedian.both_tviracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; ur: cold both (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for j=0,n_elements(GA.TvirVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.coldMedian.both_tviracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(AR.TvirVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.coldMedian.both_tviracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; ll: hot gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for j=0,n_elements(GA.TvirVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.hotMedian.gmem_tviracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(AR.TvirVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.hotMedian.gmem_tviracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    ; lr: cold gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for j=0,n_elements(GA.TvirVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.coldMedian.gmem_tviracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(AR.TvirVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.coldMedian.gmem_tviracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot

    ; legend
    strings = textoidl("T_{max} / T_{vir,acc} < ")+string(GA.TvirVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(GA.TvirVals)),box=0,linesize=0.25,position=[10.05,13.0],charsize=!p.charsize-0.1

    ; labels
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,"Hot",alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,"Cold",alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; plot (3) - gal/stars, variable*const
  pName = sPg.plotPath + 'accRateByMethod.gal-stars.const.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps'
  start_PS, pName, /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot gal (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for j=0,n_elements(GA.TcutVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.hotMedian.gal_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(AR.TcutVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.hotMedian.gal_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; ur: cold gal (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for j=0,n_elements(GA.TcutVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.coldMedian.gal_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(AR.TcutVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.coldMedian.gal_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; ll: hot stars (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_stars,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for j=0,n_elements(GA.TcutVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.hotMedian.stars_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(AR.TcutVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.hotMedian.stars_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot

    ; lr: cold stars (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_stars,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for j=0,n_elements(GA.TcutVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.coldMedian.stars_const[j,*],sK,/nan),color=sPg.colorsG[1],line=j,/overplot
    for j=0,n_elements(AR.TcutVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.coldMedian.stars_const[j,*],sK,/nan),color=sPa.colorsA[1],line=j,/overplot

    ; legend
    strings = textoidl("T_{max} / T_{c} = ")+string(GA.TcutVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(GA.TcutVals)),box=0,linesize=0.25,position=[10.05,13.0]

    ; labels
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Galaxy (Stars Only)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Galaxy (Gas Only)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,"Hot",alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,"Cold",alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; plot (4) - gal/stars, variable*tvir
  pName = sPg.plotPath + 'accRateByMethod.gal-stars.tvir.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps'
  start_PS, pName, /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot gal (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for j=0,n_elements(GA.TvirVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.hotMedian.gal_tviracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(AR.TvirVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.hotMedian.gal_tviracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; ur: cold gal (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for j=0,n_elements(GA.TvirVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.coldMedian.gal_tviracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(AR.TvirVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.coldMedian.gal_tviracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; ll: hot stars (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_stars,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for j=0,n_elements(GA.TvirVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.hotMedian.stars_tviracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(AR.TvirVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.hotMedian.stars_tviracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot

    ; lr: cold gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_stars,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for j=0,n_elements(GA.TvirVals)-1 do $
      cgPlot,GA.logMassBinCen,smooth(GA.coldMedian.stars_tviracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(AR.TvirVals)-1 do $
      cgPlot,AR.logMassBinCen,smooth(AR.coldMedian.stars_tviracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot

    ; legend
    strings = textoidl("T_{max} / T_{vir,acc} < ")+string(GA.TvirVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(GA.TvirVals)),box=0,linesize=0.25,position=[10.05,13.0],charsize=!p.charsize-0.1

    ; labels
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Galaxy (Stars Only)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Galaxy (Gas Only)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,"Hot",alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,"Cold",alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS  
  
  ; plot (5) - median lines (gadget vs arepo) both/gmem separated with specific inset
  if 0 then begin
  pName = sPg.plotPath + 'accRateByMethod.both-gmem.sp.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps'
  start_PS, pName, xs = 7.5, ys = 10.0
    
    x0 = 0.15 & x1 = 0.95
    y0 = 0.1 & y1 = 0.525 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; uc
                [x0,y0,x1,y1]  ) ; lc
                
    possub = list( [0.615,0.57,0.895,0.74] ,$ ; upper sub
                   [0.615,0.14,0.895,0.31]  ) ; lower sub
   
    ; both
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    ; gadget (both)
    cgPlot,GA.logMassBinCen,smooth(GA.coldMedian.both_tViracc[tVirInd,*],sK,/nan),color=sPg.colorsG[cInd],line=lines[0],/overplot
    cgPlot,GA.logMassBinCen,smooth(GA.hotMedian.both_tViracc[tVirInd,*],sK,/nan),color=sPg.colorsG[cInd],line=lines[2],/overplot
    cgPlot,GA.logMassBinCen,smooth(GA.totalMedian.both,sK,/nan),color=sPg.colorsG[cInd],line=lines[1],/overplot
   
    ; arepo (both)
    cgPlot,AR.logMassBinCen,smooth(AR.coldMedian.both_tViracc[tVirInd,*],sK,/nan),color=sPa.colorsA[cInd],line=lines[0],/overplot
    cgPlot,AR.logMassBinCen,smooth(AR.hotMedian.both_tViracc[tVirInd,*],sK,/nan),color=sPa.colorsA[cInd],line=lines[2],/overplot
    cgPlot,AR.logMassBinCen,smooth(AR.totalMedian.both,sK,/nan),color=sPa.colorsA[cInd],line=lines[1],/overplot
    
    ; legend
    legend,[textoidl("T_{max} < 1.0 T_{vir,acc}"),textoidl("T_{max} > 1.0 T_{vir,acc}"),'total'],$
      linestyle=[lines[0],lines[2],lines[1]],box=0,linesize=0.25,/top,/left,spacing=!p.charsize+0.5
    
    ; SUBPLOT: specific (both)
    polyfill,[possub[0,0]-0.03,possub[0,0]-0.03,possub[0,2]+0.002,possub[0,2]+0.002,possub[0,0]-0.03],$
             [possub[0,1],possub[0,3],possub[0,3],possub[0,1],possub[0,1]],$
             /normal,color=cgColor('white')
             
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=[-13.5,-10.0],/xs,/ys,ytitle="",xtitle="",$
      pos=possub[0],/noerase,charsize=!p.charsize-0.4,ytickv=[-10,-11,-12,-13],yticks=3

    ; gadget (both)
    cgPlot,GA.logMassBinCen,smooth(alog10(GA.coldMedian.both_tViracc[tVirInd,*]/10.0^GA.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[0],/overplot
    cgPlot,GA.logMassBinCen,smooth(alog10(GA.hotMedian.both_tViracc[tVirInd,*]/10.0^GA.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[2],/overplot
    cgPlot,GA.logMassBinCen,smooth(alog10(GA.totalMedian.both/10.0^GA.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[1],/overplot
   
    ; arepo (both)
    cgPlot,AR.logMassBinCen,smooth(alog10(AR.coldMedian.both_tViracc[tVirInd,*]/10.0^GA.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[0],/overplot
    cgPlot,AR.logMassBinCen,smooth(alog10(AR.hotMedian.both_tViracc[tVirInd,*]/10.0^GA.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[2],/overplot
    cgPlot,AR.logMassBinCen,smooth(alog10(AR.totalMedian.both/10.0^GA.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[1],/overplot
    
    ; gmem
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),pos=pos[1],/noerase
    
    ; gadget (gmem)
    cgPlot,GA.logMassBinCen,smooth(GA.coldMedian.gmem_tViracc[tVirInd,*],sK,/nan),color=sPg.colorsG[cInd],line=lines[0],/overplot
    cgPlot,GA.logMassBinCen,smooth(GA.hotMedian.gmem_tViracc[tVirInd,*],sK,/nan),color=sPg.colorsG[cInd],line=lines[2],/overplot
    cgPlot,GA.logMassBinCen,smooth(GA.totalMedian.gmem,sK,/nan),color=sPg.colorsG[cInd],line=lines[1],/overplot
    
    ; arepo (gmem)
    cgPlot,AR.logMassBinCen,smooth(AR.coldMedian.gmem_tViracc[tVirInd,*],sK,/nan),color=sPa.colorsA[cInd],line=lines[0],/overplot
    cgPlot,AR.logMassBinCen,smooth(AR.hotMedian.gmem_tViracc[tVirInd,*],sK,/nan),color=sPa.colorsA[cInd],line=lines[2],/overplot
    cgPlot,AR.logMassBinCen,smooth(AR.totalMedian.gmem,sK,/nan),color=sPa.colorsA[cInd],line=lines[1],/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],$
      box=0,/top,/left,spacing=!p.charsize+0.5
    
    ; SUBPLOT: specific (gmem)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=[-12.0,-10.5],/xs,/ys,ytitle="",xtitle="",$
      pos=possub[1],/noerase,charsize=!p.charsize-0.4,ytickv=[-10,-11,-12],yticks=2
    
    ; gadget (gmem)
    cgPlot,GA.logMassBinCen,smooth(alog10(GA.coldMedian.gmem_tViracc[tVirInd,*]/10.0^GA.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[0],/overplot
    cgPlot,GA.logMassBinCen,smooth(alog10(GA.hotMedian.gmem_tViracc[tVirInd,*]/10.0^GA.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[2],/overplot
    cgPlot,GA.logMassBinCen,smooth(alog10(GA.totalMedian.gmem/10.0^GA.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[1],/overplot
    
    ; arepo (gmem)
    cgPlot,AR.logMassBinCen,smooth(alog10(AR.coldMedian.gmem_tViracc[tVirInd,*]/10.0^GA.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[0],/overplot
    cgPlot,AR.logMassBinCen,smooth(alog10(AR.hotMedian.gmem_tViracc[tVirInd,*]/10.0^GA.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[2],/overplot
    cgPlot,AR.logMassBinCen,smooth(alog10(AR.totalMedian.gmem/10.0^GA.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[1],/overplot
    
    ; labels
    cgText,0.05,0.5,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1+0.02,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x1+0.02,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    
  end_PS
  endif ;0  
  
end

; plotByMode(): plot cold fractions and accretion rates separated out into the four modes and also
;               a view with decomposed galaxy=gal+stars (one resolution only)

pro plotByMode

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect   = 'pri'
  redshift   = 2.0
  res        = 512
  timeWindow = 1000.0

  sK     = 3 ; smoothing kernel size
  cInd   = 1 ; color index
  virInd = 0 ; 1.0 Tvir
  conInd = 1 ; 5.5 const
  
  xrange = [10.0,12.0]
  yrange = [0.0,1.025]

  ; load
  sPg = simParams(res=res,run='gadget',redshift=redshift)
  sPa = simParams(res=res,run='tracer',redshift=redshift)

  GA_all      = haloMassBinValues(sP=sPg,sgSelect=sgSelect,timeWindow=timeWindow,accMode='all')
  AR_all      = haloMassBinValues(sP=sPa,sgSelect=sgSelect,timeWindow=timeWindow,accMode='all')
  GA_smooth   = haloMassBinValues(sP=sPg,sgSelect=sgSelect,timeWindow=timeWindow,accMode='smooth')
  AR_smooth   = haloMassBinValues(sP=sPa,sgSelect=sgSelect,timeWindow=timeWindow,accMode='smooth')
  GA_stripped = haloMassBinValues(sP=sPg,sgSelect=sgSelect,timeWindow=timeWindow,accMode='stripped')
  AR_stripped = haloMassBinValues(sP=sPa,sgSelect=sgSelect,timeWindow=timeWindow,accMode='stripped')
  GA_clumpy   = haloMassBinValues(sP=sPg,sgSelect=sgSelect,timeWindow=timeWindow,accMode='clumpy')
  AR_clumpy   = haloMassBinValues(sP=sPa,sgSelect=sgSelect,timeWindow=timeWindow,accMode='clumpy')

  ; --- cold fractions ---

  ; plot (1) - 2x2 of (both,gmem)x(const,tvir)
  pName = sPg.plotPath + 'coldFracByMode.both-gmem.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(fix(timeWindow))+'.eps'
  start_PS, pName, /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.10 & y1 = 0.50 & y2 = 0.90
    
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: gal tvir (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4,pos=pos[0]
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_all.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_all.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_smooth.logMassBinCen,smooth(GA_smooth.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_smooth.logMassBinCen,smooth(AR_smooth.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,GA_clumpy.logMassBinCen,smooth(GA_clumpy.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_clumpy.logMassBinCen,smooth(AR_clumpy.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_stripped.logMassBinCen,smooth(GA_stripped.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_stripped.logMassBinCen,smooth(AR_stripped.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    ; top left axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    cgAxis,/xaxis,xrange=tvirrange,/xs
    
    ; ur: gal const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_all.fracMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_all.fracMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_smooth.logMassBinCen,smooth(GA_smooth.fracMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_smooth.logMassBinCen,smooth(AR_smooth.fracMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,GA_clumpy.logMassBinCen,smooth(GA_clumpy.fracMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_clumpy.logMassBinCen,smooth(AR_clumpy.fracMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_stripped.logMassBinCen,smooth(GA_stripped.fracMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_stripped.logMassBinCen,smooth(AR_stripped.fracMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    ; top right axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    tvirrange[1] = round(10*tvirrange[1])/10.0
    cgAxis,/xaxis,xrange=tvirrange,/xs

    ; ll: gmem tvir (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_all.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_all.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_smooth.logMassBinCen,smooth(GA_smooth.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_smooth.logMassBinCen,smooth(AR_smooth.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,GA_clumpy.logMassBinCen,smooth(GA_clumpy.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_clumpy.logMassBinCen,smooth(AR_clumpy.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_stripped.logMassBinCen,smooth(GA_stripped.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_stripped.logMassBinCen,smooth(AR_stripped.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    ; legend
    legend,['all','smooth','clumpy','stripped'],linestyle=[0,1,2,3],box=0,linesize=0.4,$
      position=[11.2,0.95],charsize=!p.charsize-0.2

    ; lr: gmem const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_all.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_all.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_smooth.logMassBinCen,smooth(GA_smooth.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_smooth.logMassBinCen,smooth(AR_smooth.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,GA_clumpy.logMassBinCen,smooth(GA_clumpy.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_clumpy.logMassBinCen,smooth(AR_clumpy.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_stripped.logMassBinCen,smooth(GA_stripped.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_stripped.logMassBinCen,smooth(AR_stripped.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/bottom,/left

    ; labels
    cgText,0.05,y1,"Cold Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.01,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    cgText,x1,y2+0.06,textoidl("T_{vir} [_{ }log K_{ }] (z=2)"),alignment=0.5,/normal
    
    cgText,x2+0.02,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.02,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x1,x2])+0.05,y2+0.06,textoidl("T_{max} < T_c = 5.5"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1])-0.05,y2+0.06,textoidl("T_{max} < 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; plot (2) - 2x2 of (gal,stars)x(const,tvir)
  pName = sPg.plotPath + 'coldFracByMode.gal-stars.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(fix(timeWindow))+'.eps'
  start_PS, pName, /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.10 & y1 = 0.50 & y2 = 0.90
    
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: gal gas only tvir (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4,pos=pos[0]
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_all.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_all.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_smooth.logMassBinCen,smooth(GA_smooth.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_smooth.logMassBinCen,smooth(AR_smooth.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,GA_clumpy.logMassBinCen,smooth(GA_clumpy.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_clumpy.logMassBinCen,smooth(AR_clumpy.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_stripped.logMassBinCen,smooth(GA_stripped.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_stripped.logMassBinCen,smooth(AR_stripped.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    ; top left axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    cgAxis,/xaxis,xrange=tvirrange,/xs
    
    ; ur: gal gas only const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_all.fracMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_all.fracMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_smooth.logMassBinCen,smooth(GA_smooth.fracMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_smooth.logMassBinCen,smooth(AR_smooth.fracMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,GA_clumpy.logMassBinCen,smooth(GA_clumpy.fracMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_clumpy.logMassBinCen,smooth(AR_clumpy.fracMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_stripped.logMassBinCen,smooth(GA_stripped.fracMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_stripped.logMassBinCen,smooth(AR_stripped.fracMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    ; top right axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    tvirrange[1] = round(10*tvirrange[1])/10.0
    cgAxis,/xaxis,xrange=tvirrange,/xs

    ; ll: stars only tvir (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_all.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_all.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_smooth.logMassBinCen,smooth(GA_smooth.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_smooth.logMassBinCen,smooth(AR_smooth.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,GA_clumpy.logMassBinCen,smooth(GA_clumpy.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_clumpy.logMassBinCen,smooth(AR_clumpy.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_stripped.logMassBinCen,smooth(GA_stripped.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_stripped.logMassBinCen,smooth(AR_stripped.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    ; legend
    legend,['all','smooth','clumpy','stripped'],linestyle=[0,1,2,3],box=0,linesize=0.4,$
      position=[11.2,0.95],charsize=!p.charsize-0.2

    ; lr: stars only const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_all.fracMedian.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_all.fracMedian.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_smooth.logMassBinCen,smooth(GA_smooth.fracMedian.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_smooth.logMassBinCen,smooth(AR_smooth.fracMedian.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,GA_clumpy.logMassBinCen,smooth(GA_clumpy.fracMedian.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_clumpy.logMassBinCen,smooth(AR_clumpy.fracMedian.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_stripped.logMassBinCen,smooth(GA_stripped.fracMedian.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_stripped.logMassBinCen,smooth(AR_stripped.fracMedian.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/bottom,/left

    ; labels
    cgText,0.05,y1,"Cold Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.01,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    cgText,x1,y2+0.06,textoidl("T_{vir} [_{ }log K_{ }] (z=2)"),alignment=0.5,/normal
    
    cgText,x2+0.02,mean([y0,y1]),"Galaxy (Stars Only)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.02,mean([y1,y2]),"Galaxy (Gas Only)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x1,x2])+0.05,y2+0.06,textoidl("T_{max} < T_c = 5.5"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1])-0.05,y2+0.06,textoidl("T_{max} < 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS

  ; --- accretion rates ---

  yrange = [0.01,20.0]
  yrange_stars = 0.1*yrange

  ; plot (3) - accretion rate decomposed galaxy=gas+stars variable*tviracc
  pName = sPg.plotPath + 'accRateByMode.gal-stars.tvir.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(fix(timeWindow))+'.eps'
  start_PS, pName, /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot gal gas (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_all.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_all.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_smooth.logMassBinCen,smooth(GA_smooth.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_smooth.logMassBinCen,smooth(AR_smooth.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,GA_clumpy.logMassBinCen,smooth(GA_clumpy.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_clumpy.logMassBinCen,smooth(AR_clumpy.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_stripped.logMassBinCen,smooth(GA_stripped.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_stripped.logMassBinCen,smooth(AR_stripped.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; ur: cold gal gas (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_all.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_all.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_smooth.logMassBinCen,smooth(GA_smooth.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_smooth.logMassBinCen,smooth(AR_smooth.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,GA_clumpy.logMassBinCen,smooth(GA_clumpy.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_clumpy.logMassBinCen,smooth(AR_clumpy.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_stripped.logMassBinCen,smooth(GA_stripped.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_stripped.logMassBinCen,smooth(AR_stripped.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; ll: hot stars (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_stars,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_all.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_all.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_smooth.logMassBinCen,smooth(GA_smooth.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_smooth.logMassBinCen,smooth(AR_smooth.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,GA_clumpy.logMassBinCen,smooth(GA_clumpy.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_clumpy.logMassBinCen,smooth(AR_clumpy.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_stripped.logMassBinCen,smooth(GA_stripped.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_stripped.logMassBinCen,smooth(AR_stripped.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; lr: cold stars (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_stars,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_all.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_all.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_smooth.logMassBinCen,smooth(GA_smooth.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_smooth.logMassBinCen,smooth(AR_smooth.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,GA_clumpy.logMassBinCen,smooth(GA_clumpy.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_clumpy.logMassBinCen,smooth(AR_clumpy.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_stripped.logMassBinCen,smooth(GA_stripped.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_stripped.logMassBinCen,smooth(AR_stripped.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; legend
    legend,['all','smooth','clumpy','stripped'],linestyle=[0,1,2,3],box=0,linesize=0.4,$
      /top,/left,charsize=!p.charsize-0.2
      
    ; labels
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Galaxy (Stars Only)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Galaxy (Gas Only)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,textoidl("T_{max} > 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,textoidl("T_{max} < 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; plot (4) - accretion rate galaxy/halo atmosphere variable*tvir
  pName = sPg.plotPath + 'accRateByMode.both-gmem.tvir.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(fix(timeWindow))+'.eps'
  start_PS, pName, /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot both (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_all.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_all.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_smooth.logMassBinCen,smooth(GA_smooth.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_smooth.logMassBinCen,smooth(AR_smooth.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,GA_clumpy.logMassBinCen,smooth(GA_clumpy.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_clumpy.logMassBinCen,smooth(AR_clumpy.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_stripped.logMassBinCen,smooth(GA_stripped.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_stripped.logMassBinCen,smooth(AR_stripped.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; ur: cold both (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_all.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_all.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_smooth.logMassBinCen,smooth(GA_smooth.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_smooth.logMassBinCen,smooth(AR_smooth.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,GA_clumpy.logMassBinCen,smooth(GA_clumpy.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_clumpy.logMassBinCen,smooth(AR_clumpy.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_stripped.logMassBinCen,smooth(GA_stripped.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_stripped.logMassBinCen,smooth(AR_stripped.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; ll: hot gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_all.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_all.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_smooth.logMassBinCen,smooth(GA_smooth.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_smooth.logMassBinCen,smooth(AR_smooth.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,GA_clumpy.logMassBinCen,smooth(GA_clumpy.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_clumpy.logMassBinCen,smooth(AR_clumpy.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_stripped.logMassBinCen,smooth(GA_stripped.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_stripped.logMassBinCen,smooth(AR_stripped.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; lr: cold gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_all.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_all.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_smooth.logMassBinCen,smooth(GA_smooth.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_smooth.logMassBinCen,smooth(AR_smooth.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,GA_clumpy.logMassBinCen,smooth(GA_clumpy.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_clumpy.logMassBinCen,smooth(AR_clumpy.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_stripped.logMassBinCen,smooth(GA_stripped.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_stripped.logMassBinCen,smooth(AR_stripped.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; legend
    legend,['all','smooth','clumpy','stripped'],linestyle=[0,1,2,3],box=0,linesize=0.4,$
      /top,/left,charsize=!p.charsize-0.2
      
    ; labels
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,textoidl("T_{max} > 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,textoidl("T_{max} < 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  yrange = [0.0,1.0]

  ; plot (5) - fraction, accretion rate decomposed galaxy=gas+stars variable*tvir
  pName = sPg.plotPath + 'accRateFracByMode.gal-stars.tvir.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(fix(timeWindow))+'.eps'
  start_PS, pName, /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot gal gas (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_smooth.hotMedian.gal_tviracc[virInd,*]/GA_all.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_smooth.hotMedian.gal_tviracc[virInd,*]/AR_all.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot
    cgPlot,GA_all.logMassBinCen,smooth(GA_clumpy.hotMedian.gal_tviracc[virInd,*]/GA_all.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_clumpy.hotMedian.gal_tviracc[virInd,*]/AR_all.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_all.logMassBinCen,smooth(GA_stripped.hotMedian.gal_tviracc[virInd,*]/GA_all.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_stripped.hotMedian.gal_tviracc[virInd,*]/AR_all.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; ur: cold gal gas (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_smooth.coldMedian.gal_tviracc[virInd,*]/GA_all.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_smooth.coldMedian.gal_tviracc[virInd,*]/AR_all.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot
    cgPlot,GA_all.logMassBinCen,smooth(GA_clumpy.coldMedian.gal_tviracc[virInd,*]/GA_all.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_clumpy.coldMedian.gal_tviracc[virInd,*]/AR_all.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_all.logMassBinCen,smooth(GA_stripped.coldMedian.gal_tviracc[virInd,*]/GA_all.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_stripped.coldMedian.gal_tviracc[virInd,*]/AR_all.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; ll: hot stars (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="",xtitle="",pos=pos[2],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_smooth.hotMedian.stars_tviracc[virInd,*]/GA_all.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_smooth.hotMedian.stars_tviracc[virInd,*]/AR_all.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot
    cgPlot,GA_all.logMassBinCen,smooth(GA_clumpy.hotMedian.stars_tviracc[virInd,*]/GA_all.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_clumpy.hotMedian.stars_tviracc[virInd,*]/AR_all.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_all.logMassBinCen,smooth(GA_stripped.hotMedian.stars_tviracc[virInd,*]/GA_all.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_stripped.hotMedian.stars_tviracc[virInd,*]/AR_all.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; lr: cold stars (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_smooth.coldMedian.stars_tviracc[virInd,*]/GA_all.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_smooth.coldMedian.stars_tviracc[virInd,*]/AR_all.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot
    cgPlot,GA_all.logMassBinCen,smooth(GA_clumpy.coldMedian.stars_tviracc[virInd,*]/GA_all.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_clumpy.coldMedian.stars_tviracc[virInd,*]/AR_all.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_all.logMassBinCen,smooth(GA_stripped.coldMedian.stars_tviracc[virInd,*]/GA_all.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_stripped.coldMedian.stars_tviracc[virInd,*]/AR_all.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; legend
    legend,['smooth','clumpy','stripped'],linestyle=[1,2,3],box=0,linesize=0.4,position=[11.2,0.95],charsize=!p.charsize-0.2
      
    ; labels
    cgText,0.05,y1,"Fraction of Gas Accretion Rate",alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Galaxy (Stars Only)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Galaxy (Gas Only)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,textoidl("T_{max} > 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,textoidl("T_{max} < 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; plot (6) - fraction, accretion rate galaxy/halo atmosphere variable*tviracc
  pName = sPg.plotPath + 'accRateFracByMode.both-gmem.tvir.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(fix(timeWindow))+'.eps'
  start_PS, pName, /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot gal (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_smooth.hotMedian.both_tviracc[virInd,*]/GA_all.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_smooth.hotMedian.both_tviracc[virInd,*]/AR_all.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot
    cgPlot,GA_all.logMassBinCen,smooth(GA_clumpy.hotMedian.both_tviracc[virInd,*]/GA_all.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_clumpy.hotMedian.both_tviracc[virInd,*]/AR_all.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_all.logMassBinCen,smooth(GA_stripped.hotMedian.both_tviracc[virInd,*]/GA_all.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_stripped.hotMedian.both_tviracc[virInd,*]/AR_all.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; ur: cold gal (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_smooth.coldMedian.both_tviracc[virInd,*]/GA_all.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_smooth.coldMedian.both_tviracc[virInd,*]/AR_all.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot
    cgPlot,GA_all.logMassBinCen,smooth(GA_clumpy.coldMedian.both_tviracc[virInd,*]/GA_all.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_clumpy.coldMedian.both_tviracc[virInd,*]/AR_all.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_all.logMassBinCen,smooth(GA_stripped.coldMedian.both_tviracc[virInd,*]/GA_all.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_stripped.coldMedian.both_tviracc[virInd,*]/AR_all.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; ll: hot gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="",xtitle="",pos=pos[2],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_smooth.hotMedian.gmem_tviracc[virInd,*]/GA_all.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_smooth.hotMedian.gmem_tviracc[virInd,*]/GA_all.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot
    cgPlot,GA_all.logMassBinCen,smooth(GA_clumpy.hotMedian.gmem_tviracc[virInd,*]/GA_all.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_clumpy.hotMedian.gmem_tviracc[virInd,*]/GA_all.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_all.logMassBinCen,smooth(GA_stripped.hotMedian.gmem_tviracc[virInd,*]/GA_all.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_stripped.hotMedian.gmem_tviracc[virInd,*]/GA_all.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; lr: cold gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    cgPlot,GA_all.logMassBinCen,smooth(GA_smooth.coldMedian.gmem_tviracc[virInd,*]/GA_all.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_smooth.coldMedian.gmem_tviracc[virInd,*]/AR_all.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot
    cgPlot,GA_all.logMassBinCen,smooth(GA_clumpy.coldMedian.gmem_tviracc[virInd,*]/GA_all.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_clumpy.coldMedian.gmem_tviracc[virInd,*]/AR_all.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_all.logMassBinCen,smooth(GA_stripped.coldMedian.gmem_tviracc[virInd,*]/GA_all.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_all.logMassBinCen,smooth(AR_stripped.coldMedian.gmem_tviracc[virInd,*]/AR_all.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; legend
    legend,['smooth','clumpy','stripped'],linestyle=[1,2,3],box=0,linesize=0.4,position=[11.2,0.95],charsize=!p.charsize-0.2
      
    ; labels
    cgText,0.05,y1,"Fraction of Gas Accretion Rate",alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,textoidl("T_{max} > 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,textoidl("T_{max} < 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
end

; plotByTimeWindow(): plot cold fractions and accretion rates for multiple time windows (one resolution only)

pro plotByTimeWindow

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect = 'pri'
  accMode  = 'smooth'
  redshift = 2.0
  res      = 256

  ; load
  sPg = simParams(res=res,run='gadget',redshift=redshift)
  sPa = simParams(res=res,run='tracer',redshift=redshift)

  GA_tw0      = haloMassBinValues(sP=sPg,sgSelect=sgSelect,timeWindow='all',accMode=accMode)
  AR_tw0      = haloMassBinValues(sP=sPa,sgSelect=sgSelect,timeWindow='all',accMode=accMode)
  GA_tw250    = haloMassBinValues(sP=sPg,sgSelect=sgSelect,timeWindow=250.0,accMode=accMode)
  AR_tw250    = haloMassBinValues(sP=sPa,sgSelect=sgSelect,timeWindow=250.0,accMode=accMode)
  ;GA_tw500    = haloMassBinValues(sP=sPg,sgSelect=sgSelect,timeWindow=500.0,accMode=accMode)
  ;AR_tw500    = haloMassBinValues(sP=sPa,sgSelect=sgSelect,timeWindow=500.0,accMode=accMode)
  GA_tw1000   = haloMassBinValues(sP=sPg,sgSelect=sgSelect,timeWindow=1000.0,accMode=accMode)
  AR_tw1000   = haloMassBinValues(sP=sPa,sgSelect=sgSelect,timeWindow=1000.0,accMode=accMode)
  
  sK     = 3 ; smoothing kernel size
  cInd   = 1 ; color index
  virInd = 0 ; 1.0 Tvir
  conInd = 1 ; 5.5 const
  
  xrange = [10.0,12.0]
  yrange = [0.0,1.025]

  ; --- cold fractions ---

  ; plot (1) - 2x2 of (both,gmem)x(const,tvir)
  pName = sPg.plotPath + 'coldFracByTW.both-gmem.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_am-'+accMode+'.eps'
  start_PS, pName, /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.10 & y1 = 0.50 & y2 = 0.90
    
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: gal tvir (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4,pos=pos[0]
    
    cgPlot,GA_tw0.logMassBinCen,smooth(GA_tw0.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_tw0.logMassBinCen,smooth(AR_tw0.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_tw250.logMassBinCen,smooth(GA_tw250.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_tw250.logMassBinCen,smooth(AR_tw250.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot
    ;cgPlot,GA_tw500.logMassBinCen,smooth(GA_tw500.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    ;cgPlot,AR_tw500.logMassBinCen,smooth(AR_tw500.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_tw1000.logMassBinCen,smooth(GA_tw1000.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_tw1000.logMassBinCen,smooth(AR_tw1000.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    ; top left axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    cgAxis,/xaxis,xrange=tvirrange,/xs
    
    ; ur: gal const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    cgPlot,GA_tw0.logMassBinCen,smooth(GA_tw0.fracMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_tw0.logMassBinCen,smooth(AR_tw0.fracMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_tw250.logMassBinCen,smooth(GA_tw250.fracMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_tw250.logMassBinCen,smooth(AR_tw250.fracMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot
    ;cgPlot,GA_tw500.logMassBinCen,smooth(GA_tw500.fracMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    ;cgPlot,AR_tw500.logMassBinCen,smooth(AR_tw500.fracMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_tw1000.logMassBinCen,smooth(GA_tw1000.fracMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_tw1000.logMassBinCen,smooth(AR_tw1000.fracMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    ; top right axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    tvirrange[1] = round(10*tvirrange[1])/10.0
    cgAxis,/xaxis,xrange=tvirrange,/xs

    ; legend
    legend,['250 Myr','500 Myr','1000 Myr','integrated'],linestyle=[1,2,3,0],box=0,linesize=0.4,$
      position=[10.1,0.4],charsize=!p.charsize-0.2

    ; ll: gmem tvir (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    cgPlot,GA_tw0.logMassBinCen,smooth(GA_tw0.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_tw0.logMassBinCen,smooth(AR_tw0.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_tw250.logMassBinCen,smooth(GA_tw250.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_tw250.logMassBinCen,smooth(AR_tw250.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot
    ;cgPlot,GA_tw500.logMassBinCen,smooth(GA_tw500.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    ;cgPlot,AR_tw500.logMassBinCen,smooth(AR_tw500.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_tw1000.logMassBinCen,smooth(GA_tw1000.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_tw1000.logMassBinCen,smooth(AR_tw1000.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    ; lr: gmem const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    cgPlot,GA_tw0.logMassBinCen,smooth(GA_tw0.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_tw0.logMassBinCen,smooth(AR_tw0.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_tw250.logMassBinCen,smooth(GA_tw250.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_tw250.logMassBinCen,smooth(AR_tw250.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot
    ;cgPlot,GA_tw500.logMassBinCen,smooth(GA_tw500.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    ;cgPlot,AR_tw500.logMassBinCen,smooth(AR_tw500.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_tw1000.logMassBinCen,smooth(GA_tw1000.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_tw1000.logMassBinCen,smooth(AR_tw1000.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/bottom,/left

    ; labels
    cgText,0.05,y1,"Cold Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.01,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    cgText,x1,y2+0.06,textoidl("T_{vir} [_{ }log K_{ }] (z=2)"),alignment=0.5,/normal
    
    cgText,x2+0.02,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.02,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x1,x2])+0.05,y2+0.06,textoidl("T_{max} < T_c = 5.5"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1])-0.05,y2+0.06,textoidl("T_{max} < 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; --- accretion rates ---
  
  yrange = [0.01,20.0]
  yrange_stars = 0.1*yrange

  ; plot (4) - accretion rate both/gmem variable*tvir
  pName = sPg.plotPath + 'accRateByTW.both-gmem.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_am-'+accMode+'.eps'
  start_PS, pName, /big
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot both (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    cgPlot,GA_tw0.logMassBinCen,smooth(GA_tw0.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_tw0.logMassBinCen,smooth(AR_tw0.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_tw250.logMassBinCen,smooth(GA_tw250.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_tw250.logMassBinCen,smooth(AR_tw250.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot
    ;cgPlot,GA_tw500.logMassBinCen,smooth(GA_tw500.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    ;cgPlot,AR_tw500.logMassBinCen,smooth(AR_tw500.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_tw1000.logMassBinCen,smooth(GA_tw1000.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_tw1000.logMassBinCen,smooth(AR_tw1000.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; ur: cold both (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    cgPlot,GA_tw0.logMassBinCen,smooth(GA_tw0.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_tw0.logMassBinCen,smooth(AR_tw0.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_tw250.logMassBinCen,smooth(GA_tw250.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_tw250.logMassBinCen,smooth(AR_tw250.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot
    ;cgPlot,GA_tw500.logMassBinCen,smooth(GA_tw500.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    ;cgPlot,AR_tw500.logMassBinCen,smooth(AR_tw500.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_tw1000.logMassBinCen,smooth(GA_tw1000.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_tw1000.logMassBinCen,smooth(AR_tw1000.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; legend
    if accMode eq 'smooth' then $
      legend,['250 Myr','500 Myr','1000 Myr','integrated'],linestyle=[1,2,3,0],box=0,linesize=0.4,$
        /top,/left,charsize=!p.charsize-0.2
    if accMode eq 'all' then $
      legend,['250 Myr','500 Myr','1000 Myr','integrated'],linestyle=[1,2,3,0],box=0,linesize=0.4,$
        /bottom,/right,charsize=!p.charsize-0.2
    
    ; ll: hot gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    cgPlot,GA_tw0.logMassBinCen,smooth(GA_tw0.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_tw0.logMassBinCen,smooth(AR_tw0.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_tw250.logMassBinCen,smooth(GA_tw250.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_tw250.logMassBinCen,smooth(AR_tw250.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot
    ;cgPlot,GA_tw500.logMassBinCen,smooth(GA_tw500.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    ;cgPlot,AR_tw500.logMassBinCen,smooth(AR_tw500.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_tw1000.logMassBinCen,smooth(GA_tw1000.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_tw1000.logMassBinCen,smooth(AR_tw1000.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; lr: cold gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    cgPlot,GA_tw0.logMassBinCen,smooth(GA_tw0.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,AR_tw0.logMassBinCen,smooth(AR_tw0.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,GA_tw250.logMassBinCen,smooth(GA_tw250.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,AR_tw250.logMassBinCen,smooth(AR_tw250.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot
    ;cgPlot,GA_tw500.logMassBinCen,smooth(GA_tw500.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    ;cgPlot,AR_tw500.logMassBinCen,smooth(AR_tw500.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,GA_tw1000.logMassBinCen,smooth(GA_tw1000.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,AR_tw1000.logMassBinCen,smooth(AR_tw1000.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
      
    ; labels
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,textoidl("T_{max} > 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,textoidl("T_{max} < 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
end

; plotByRes(): plot cold fractions and accretion rates for multiple resolutions

pro plotByRes

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect    = 'pri'
  accMode     = 'smooth'
  redshift    = 2.0
  timeWindow  = 1000.0 ; Myr

  ; load
  sPg = simParams(res=512,run='gadget',redshift=redshift)
  sPa = simParams(res=512,run='tracer',redshift=redshift)
  GA_512 = haloMassBinValues(sP=sPg,sgSelect=sgSelect,timeWindow=timeWindow,accMode=accMode)
  AR_512 = haloMassBinValues(sP=sPa,sgSelect=sgSelect,timeWindow=timeWindow,accMode=accMode)
  
  sPg = simParams(res=256,run='gadget',redshift=redshift)
  sPa = simParams(res=256,run='tracer',redshift=redshift)
  GA_256 = haloMassBinValues(sP=sPg,sgSelect=sgSelect,timeWindow=timeWindow,accMode=accMode)
  AR_256 = haloMassBinValues(sP=sPa,sgSelect=sgSelect,timeWindow=timeWindow,accMode=accMode)
  
  sPg = simParams(res=128,run='gadget',redshift=redshift)
  sPa = simParams(res=128,run='tracer',redshift=redshift)
  GA_128 = haloMassBinValues(sP=sPg,sgSelect=sgSelect,timeWindow=timeWindow,accMode=accMode)
  AR_128 = haloMassBinValues(sP=sPa,sgSelect=sgSelect,timeWindow=timeWindow,accMode=accMode)

  lines  = [1,0,2] ; tvircur,tviracc,const5.5 or 128,512,256
  thicks = [4,6,8] ; 128,256,512  
  sK     = 3 ; smoothing kernel size
  cInd   = 1 ; color index
  virInd = 0 ; 1.0 Tvir
  conInd = 1 ; 5.5 const
  
  xrange = [10.0,12.0]
  yrange = [0.0,1.025]
  
  ; --- cold fractions ---
  
  ; plot (0)
  start_PS, sPg.plotPath + 'coldFrac0.gal.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps'
            
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),$
      title="Gas Accretion on to Galaxies at z=2"
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
        
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot  
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[2],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[2],/overplot 
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[0],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[0],/overplot 
    
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/bottom,/left
            
  end_PS
  
  ; plot (1) - gal median lines (Tmax) 1*tvircur,1*tviracc,tcut=5.5
  start_PS, sPg.plotPath + 'coldFracByRes.gal.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]")
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; gadget gal 128
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[2],thick=thicks[0],/overplot
    
    ; arepo gal 128
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[2],thick=thicks[0],/overplot
        
    ; gadget gal 256
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],thick=thicks[1],/overplot
    
    ; arepo gal 256
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],thick=thicks[1],/overplot
    
    ; gadget gal 512
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[2],thick=thicks[2],/overplot
    
    ; arepo gal 512
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[2],thick=thicks[2],/overplot
    
    ; style legend
    strings = [textoidl("T_{max} < 1.0 T_{vir,cur}"),$
               textoidl("T_{max} < 1.0 T_{vir,acc}"),$
               textoidl("T_{max} < T_{c} = 10^{5.5} K")]
    legend,strings,linestyle=[1,0,2],box=0,linesize=0.25,spacing=!p.charsize+0.5,position=[10.05,0.3],charsize=!p.charsize-0.1
    
    ; color/thick legend
    strings = [textoidl(['128^3','256^3','512^3']+' GADGET'),textoidl(['128^3','256^3','512^3']+' AREPO')]
    legend,strings,textcolors=[sPg.colorsG,sPa.colorsA],box=0,position=[11.45,0.82],charsize=!p.charsize-0.3
    
  end_PS
  
  ; plot (1b) - gmem
  start_PS, sPg.plotPath + 'coldFracByRes.gmem.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]")
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; gadget gmem 128
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[2],thick=thicks[0],/overplot
    
    ; arepo gmem 128
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[2],thick=thicks[0],/overplot
        
    ; gadget gmem 256
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],thick=thicks[1],/overplot
    
    ; arepo gmem 256
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],thick=thicks[1],/overplot
    
    ; gadget gmem 512
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[2],thick=thicks[2],/overplot
    
    ; arepo gmem 512
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[2],thick=thicks[2],/overplot
    
    ; style legend
    strings = [textoidl("T_{max} < 1.0 T_{vir,cur}"),$
               textoidl("T_{max} < 1.0 T_{vir,acc}"),$
               textoidl("T_{max} < T_{c} = 10^{5.5} K")]
    legend,strings,linestyle=[1,0,2],box=0,linesize=0.25,spacing=!p.charsize+0.5,position=[10.05,0.3],charsize=!p.charsize-0.2
    
    ; color/thick legend
    strings = [textoidl(['128^3','256^3','512^3']+' GADGET'),textoidl(['128^3','256^3','512^3']+' AREPO')]
    legend,strings,textcolors=[sPg.colorsG,sPa.colorsA],box=0,position=[11.45,0.95],charsize=!p.charsize-0.3
    
  end_PS
  
  ; plot (1c) - stars median lines (Tmax) 1*tvircur,1*tviracc,tcut=5.5
  start_PS, sPg.plotPath + 'coldFracByRes.stars.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]")
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; gadget gal 128
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.stars_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[2],thick=thicks[0],/overplot
    
    ; arepo gal 128
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.stars_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[2],thick=thicks[0],/overplot
        
    ; gadget gal 256
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.stars_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],thick=thicks[1],/overplot
    
    ; arepo gal 256
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.stars_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],thick=thicks[1],/overplot
    
    ; gadget gal 512
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.stars_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[2],thick=thicks[2],/overplot
    
    ; arepo gal 512
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.stars_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[2],thick=thicks[2],/overplot
    
    ; style legend
    strings = [textoidl("T_{max} < 1.0 T_{vir,cur}"),$
               textoidl("T_{max} < 1.0 T_{vir,acc}"),$
               textoidl("T_{max} < T_{c} = 10^{5.5} K")]
    legend,strings,linestyle=[1,0,2],box=0,linesize=0.25,spacing=!p.charsize+0.5,position=[10.05,0.3],charsize=!p.charsize-0.1
    
    ; color/thick legend
    strings = [textoidl(['128^3','256^3','512^3']+' GADGET'),textoidl(['128^3','256^3','512^3']+' AREPO')]
    legend,strings,textcolors=[sPg.colorsG,sPa.colorsA],box=0,position=[11.45,0.82],charsize=!p.charsize-0.3
    
  end_PS
  
  ; plot (1d) - both=gal+stars median lines (Tmax) 1*tvircur,1*tviracc,tcut=5.5
  start_PS, sPg.plotPath + 'coldFracByRes.both.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]")
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; gadget gal 128
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[2],thick=thicks[0],/overplot
    
    ; arepo gal 128
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[2],thick=thicks[0],/overplot
        
    ; gadget gal 256
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],thick=thicks[1],/overplot
    
    ; arepo gal 256
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],thick=thicks[1],/overplot
    
    ; gadget gal 512
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[2],thick=thicks[2],/overplot
    
    ; arepo gal 512
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[2],thick=thicks[2],/overplot
    
    ; style legend
    strings = [textoidl("T_{max} < 1.0 T_{vir,cur}"),$
               textoidl("T_{max} < 1.0 T_{vir,acc}"),$
               textoidl("T_{max} < T_{c} = 10^{5.5} K")]
    legend,strings,linestyle=[1,0,2],box=0,linesize=0.25,spacing=!p.charsize+0.5,position=[10.05,0.3],charsize=!p.charsize-0.1
    
    ; color/thick legend
    strings = [textoidl(['128^3','256^3','512^3']+' GADGET'),textoidl(['128^3','256^3','512^3']+' AREPO')]
    legend,strings,textcolors=[sPg.colorsG,sPa.colorsA],box=0,position=[11.45,0.82],charsize=!p.charsize-0.3
    
  end_PS
  
  ; plot (1e) - combined 1a+1b in a 3x2 panel split by definition type, gal/gmem (resolution lines)
  start_PS, sPg.plotPath + 'coldFracByRes.3x2.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps',/big
    
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

    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot

    ; uc - gal tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),yticks=4,pos=pos[1]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot
  
    ; ur - gal tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[2]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot  
    
    ; ll - gmem const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/noerase,$
      ytitle="",xtitle="",ytickv=[0.0,0.2,0.4,0.6,0.8,1.0],yticks=5,pos=pos[3]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot
    
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/bottom,/left
    
    ; lc - gmem tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[4]
      
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot
  
    ; lr - gmem tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[5]
      
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot
    
    legend,textoidl(['128^3','256^3','512^3']),linestyle=[lines[0],lines[2],lines[1]],box=0,linesize=0.25,$
      /bottom,/left
    
    ; labels    
    cgText,0.05,y1,"Cold Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,mean([x0,x3]),0.02,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x3+0.02,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x3+0.02,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal

    cgText,mean([x2,x3]),y2+0.03,textoidl("T_{max} < 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.03,textoidl("T_{max} < 1.0 T_{vir,cur}"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1]),y2+0.03,textoidl("T_{max} < T_{c} = 10^{5.5} K"),alignment=0.5,color=cgColor('dark gray'),/normal

   end_PS
  
  ; --- accretion rates ---
  
  yrange = [0.01,20.0]
  yrange_stars = 0.1*yrange
  
  ; plot (2e) - combined 1a+1b in a 3x2 panel split by definition type, both hot/cold (resolution lines)
  start_PS, sPg.plotPath + 'accRateByRes.both.3x2.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps',/big
    
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

    cgPlot,GA_128.logMassBinCen,smooth(GA_128.coldMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.coldMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.coldMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.coldMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.coldMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.coldMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot

    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left

    ; uc - both cold tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),yticks=4,pos=pos[1]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.coldMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.coldMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.coldMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.coldMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.coldMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.coldMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot
 
    ; ur - both cold tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[2]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.coldMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.coldMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.coldMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.coldMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.coldMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.coldMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot  
    
    if accMode eq 'smooth' then $
      legend,textoidl(['128^3','256^3','512^3']),linestyle=[lines[0],lines[2],lines[1]],box=0,linesize=0.34,position=[10.0,14.0]
    if accMode eq 'all' then $
      legend,textoidl(['128^3','256^3','512^3']),linestyle=[lines[0],lines[2],lines[1]],box=0,linesize=0.34,/bottom,/right

    ; ll - both hot const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",pos=pos[3]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.hotMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.hotMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.hotMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.hotMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.hotMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.hotMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot

    ; draw some resolution Nx element lines
    omegaBaryon = 0.044
    nElem = 5000
    yrangeTicks = [2.0,4.0]
    yrangeTicks2 = [1.5,4.0]
    elemMass128 = sPg.targetGasMass
    
    elemMasses = codeMassToLogMsun([elemMass128,elemMass128/8.0,elemMass128/(8.0)^2] * nElem / omegaBaryon)
    
    cgPlot,[elemMasses[2],elemMasses[2]],yrangeTicks,line=lines[1],color=cgColor('black'),thick=!p.thick+1,/overplot ; 512
    cgPlot,[elemMasses[1],elemMasses[1]],yrangeTicks2,line=lines[2],color=cgColor('black'),thick=!p.thick+1,/overplot ; 256
 
    ; lc - both hot tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[4]
      
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.hotMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.hotMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.hotMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.hotMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.hotMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.hotMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot
 
    cgPlot,[elemMasses[2],elemMasses[2]],yrangeTicks,line=lines[1],color=cgColor('black'),thick=!p.thick+1,/overplot ; 512
    cgPlot,[elemMasses[1],elemMasses[1]],yrangeTicks2,line=lines[2],color=cgColor('black'),thick=!p.thick+1,/overplot ; 256
 
    ; lr - both hot tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[5]
      
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.hotMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.hotMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.hotMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.hotMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.hotMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.hotMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot  
    
    cgPlot,[elemMasses[2],elemMasses[2]],yrangeTicks,line=lines[1],color=cgColor('black'),thick=!p.thick+1,/overplot ; 512
    cgPlot,[elemMasses[1],elemMasses[1]],yrangeTicks2,line=lines[2],color=cgColor('black'),thick=!p.thick+1,/overplot ; 256
 
    ; labels    
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,mean([x0,x3]),0.02,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x3+0.02,mean([y0,y1]),"Galaxy (Hot)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x3+0.02,mean([y1,y2]),"Galaxy (Cold)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal

    cgText,mean([x2,x3]),y2+0.03,textoidl("T_{max} < 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.03,textoidl("T_{max} < 1.0 T_{vir,cur}"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1]),y2+0.03,textoidl("T_{max} < T_{c} = 10^{5.5} K"),alignment=0.5,color=cgColor('dark gray'),/normal

   end_PS

  ; plot (2ee) - both hot/cold/total (resolution lines)
  start_PS, sPg.plotPath + 'accRateByRes.both.3x3.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps',xs=9,ys=9
    
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

    cgPlot,GA_128.logMassBinCen,smooth(GA_128.totalHCMedian.both,sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.totalHCMedian.both,sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.totalHCMedian.both,sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.totalHCMedian.both,sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.totalHCMedian.both,sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.totalHCMedian.both,sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot

    ; uc - both total tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),yticks=4,pos=pos[1]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.totalHCMedian.both,sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.totalHCMedian.both,sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.totalHCMedian.both,sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.totalHCMedian.both,sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.totalHCMedian.both,sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.totalHCMedian.both,sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot
 
    ; ur - both total tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[2]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.totalHCMedian.both,sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.totalHCMedian.both,sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.totalHCMedian.both,sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.totalHCMedian.both,sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.totalHCMedian.both,sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.totalHCMedian.both,sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot  

    ; ml - both cold const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[3]

    cgPlot,GA_128.logMassBinCen,smooth(GA_128.coldMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.coldMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.coldMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.coldMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.coldMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.coldMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot

    ; mc - both cold tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),yticks=4,pos=pos[4]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.coldMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.coldMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.coldMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.coldMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.coldMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.coldMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot
 
    ; mr - both cold tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[5]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.coldMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.coldMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.coldMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.coldMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.coldMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.coldMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot  
    
    if accMode eq 'smooth' then $
      legend,textoidl(['128^3','256^3','512^3']),linestyle=[lines[0],lines[2],lines[1]],box=0,linesize=0.25,/top,/left
    if accMode eq 'all' then $
      legend,textoidl(['128^3','256^3','512^3']),linestyle=[lines[0],lines[2],lines[1]],box=0,linesize=0.25,/bottom,/right

    ; ll - both hot const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",pos=pos[6]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.hotMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.hotMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.hotMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.hotMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.hotMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.hotMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot

    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; lc - both hot tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[7]
      
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.hotMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.hotMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.hotMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.hotMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.hotMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.hotMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot
 
    ; lr - both hot tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[8]
      
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.hotMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.hotMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.hotMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.hotMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.hotMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.hotMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot  
    
    ; labels    
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,mean([x0,x3]),0.02,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x3+0.02,mean([y0,y1]),"Galaxy (Hot)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x3+0.02,mean([y1,y2]),"Galaxy (Cold)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal

    cgText,mean([x2,x3]),y2+0.03,textoidl("T_{max} < 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.03,textoidl("T_{max} < 1.0 T_{vir,cur}"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1]),y2+0.03,textoidl("T_{max} < T_{c} = 10^{5.5} K"),alignment=0.5,color=cgColor('dark gray'),/normal

   end_PS   
   
  ; plot (2f) - gal hot/cold (resolution lines)
  start_PS, sPg.plotPath + 'accRateByRes.gal.3x2.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps',/big
    
    x0 = 0.13 & x1 = 0.40 & x2 = 0.67 & x3 = 0.94
    y0 = 0.13 & y1 = 0.53 & y2 = 0.93
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; uc
                [x2,y1,x3,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1] ,$ ; lc
                [x2,y0,x3,y1] )  ; lr 

    ; ul - gal cold const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]

    cgPlot,GA_128.logMassBinCen,smooth(GA_128.coldMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.coldMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.coldMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.coldMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.coldMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.coldMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot

    ; uc - gal cold tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),yticks=4,pos=pos[1]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.coldMedian.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.coldMedian.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.coldMedian.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.coldMedian.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.coldMedian.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.coldMedian.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot
 
    ; ur - gal cold tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[2]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.coldMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.coldMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.coldMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.coldMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.coldMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.coldMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot  
    
    if accMode eq 'smooth' then $
      legend,textoidl(['128^3','256^3','512^3']),linestyle=[lines[0],lines[2],lines[1]],box=0,linesize=0.25,/top,/left
    if accMode eq 'all' then $
      legend,textoidl(['128^3','256^3','512^3']),linestyle=[lines[0],lines[2],lines[1]],box=0,linesize=0.25,/bottom,/right

    ; ll - gal hot const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",pos=pos[3]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.hotMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.hotMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.hotMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.hotMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.hotMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.hotMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot

    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; lc - gal hot tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[4]
      
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.hotMedian.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.hotMedian.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.hotMedian.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.hotMedian.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.hotMedian.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.hotMedian.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot
 
    ; lr - gal hot tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[5]
      
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.hotMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.hotMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.hotMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.hotMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.hotMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.hotMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot  
    

    ; labels    
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,mean([x0,x3]),0.02,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x3+0.02,mean([y0,y1]),"Gal Gas Only (Hot)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x3+0.02,mean([y1,y2]),"Gal Gas Only (Cold)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal

    cgText,mean([x2,x3]),y2+0.03,textoidl("T_{max} < 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.03,textoidl("T_{max} < 1.0 T_{vir,cur}"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1]),y2+0.03,textoidl("T_{max} < T_{c} = 10^{5.5} K"),alignment=0.5,color=cgColor('dark gray'),/normal

   end_PS
   
   ; plot (2g) - stars hot/cold (resolution lines)
  start_PS, sPg.plotPath + 'accRateByRes.stars.3x2.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps',/big
    
    x0 = 0.13 & x1 = 0.40 & x2 = 0.67 & x3 = 0.94
    y0 = 0.13 & y1 = 0.53 & y2 = 0.93
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; uc
                [x2,y1,x3,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1] ,$ ; lc
                [x2,y0,x3,y1] )  ; lr 

    ; ul - gal cold const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]

    cgPlot,GA_128.logMassBinCen,smooth(GA_128.coldMedian.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.coldMedian.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.coldMedian.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.coldMedian.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.coldMedian.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.coldMedian.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot

    ; uc - gal cold tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),yticks=4,pos=pos[1]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.coldMedian.stars_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.coldMedian.stars_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.coldMedian.stars_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.coldMedian.stars_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.coldMedian.stars_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.coldMedian.stars_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot
 
    ; ur - gal cold tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[2]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.coldMedian.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.coldMedian.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.coldMedian.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.coldMedian.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.coldMedian.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.coldMedian.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot  
    
    if accMode eq 'smooth' then $
      legend,textoidl(['128^3','256^3','512^3']),linestyle=[lines[0],lines[2],lines[1]],box=0,linesize=0.25,/top,/left
    if accMode eq 'all' then $
      legend,textoidl(['128^3','256^3','512^3']),linestyle=[lines[0],lines[2],lines[1]],box=0,linesize=0.25,/bottom,/right

    ; ll - gal hot const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",pos=pos[3]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.hotMedian.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.hotMedian.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.hotMedian.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.hotMedian.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.hotMedian.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.hotMedian.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot

    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; lc - gal hot tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[4]
      
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.hotMedian.stars_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.hotMedian.stars_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.hotMedian.stars_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.hotMedian.stars_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.hotMedian.stars_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.hotMedian.stars_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot
 
    ; lr - gal hot tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[5]
      
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.hotMedian.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.hotMedian.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.hotMedian.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.hotMedian.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.hotMedian.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.hotMedian.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot  
    
    ; labels    
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,mean([x0,x3]),0.02,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x3+0.02,mean([y0,y1]),"Gal Stars Only (Hot)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x3+0.02,mean([y1,y2]),"Gal Stars Only (Cold)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal

    cgText,mean([x2,x3]),y2+0.03,textoidl("T_{max} < 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.03,textoidl("T_{max} < 1.0 T_{vir,cur}"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1]),y2+0.03,textoidl("T_{max} < T_{c} = 10^{5.5} K"),alignment=0.5,color=cgColor('dark gray'),/normal

   end_PS  
   
   ; plot (2h) - gmem hot/cold (resolution lines)
  start_PS, sPg.plotPath + 'accRateByRes.gmem.3x2.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.snap)+'_tw'+str(fix(timeWindow))+'_am-'+accMode+'.eps',/big
    
    x0 = 0.13 & x1 = 0.40 & x2 = 0.67 & x3 = 0.94
    y0 = 0.13 & y1 = 0.53 & y2 = 0.93
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; uc
                [x2,y1,x3,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1] ,$ ; lc
                [x2,y0,x3,y1] )  ; lr 

    ; ul - gal cold const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]

    cgPlot,GA_128.logMassBinCen,smooth(GA_128.coldMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.coldMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.coldMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.coldMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.coldMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.coldMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot

    ; uc - gal cold tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),yticks=4,pos=pos[1]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.coldMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.coldMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.coldMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.coldMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.coldMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.coldMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot
 
    ; ur - gal cold tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[2]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.coldMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.coldMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.coldMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.coldMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.coldMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.coldMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot  
    
    ; ll - gal hot const
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",pos=pos[3]
    
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.hotMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.hotMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.hotMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.hotMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.hotMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.hotMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot

    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; lc - gal hot tvircur
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[4]
      
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.hotMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.hotMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.hotMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.hotMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.hotMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.hotMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot
 
    ; lr - gal hot tviracc
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=1,ys=1,/ylog,yminor=0,/noerase,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),xticks=3,xtickv=[10.5,11.0,11.5,12.0],pos=pos[5]
      
    cgPlot,GA_128.logMassBinCen,smooth(GA_128.hotMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],/overplot
    cgPlot,AR_128.logMassBinCen,smooth(AR_128.hotMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],/overplot
    cgPlot,GA_256.logMassBinCen,smooth(GA_256.hotMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],/overplot
    cgPlot,AR_256.logMassBinCen,smooth(AR_256.hotMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],/overplot
    cgPlot,GA_512.logMassBinCen,smooth(GA_512.hotMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],/overplot
    cgPlot,AR_512.logMassBinCen,smooth(AR_512.hotMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],/overplot  
    
    legend,textoidl(['128^3','256^3','512^3']),linestyle=[lines[0],lines[2],lines[1]],box=0,linesize=0.25,/top,/left
    
    ; labels    
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,mean([x0,x3]),0.02,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x3+0.02,mean([y0,y1]),"Halo (Hot)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x3+0.02,mean([y1,y2]),"Halo (Cold)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal

    cgText,mean([x2,x3]),y2+0.03,textoidl("T_{max} < 1.0 T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.03,textoidl("T_{max} < 1.0 T_{vir,cur}"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1]),y2+0.03,textoidl("T_{max} < T_{c} = 10^{5.5} K"),alignment=0.5,color=cgColor('dark gray'),/normal

   end_PS  
  
end
