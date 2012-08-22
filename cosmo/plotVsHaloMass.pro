; plotVsHaloMass.pro
; gas accretion project - plots as a function of halo mass
; dnelson aug.2012

; plotAccRateVsHaloMass(): plot the "cold mode fraction" vs halo mass in a few different ways
;                           using the mergerTreeSubset with accretionTimes and accretionModes

pro plotAccRateVsHaloMass

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect   = 'pri'
  accMode    = 'all' ; accretion mode: all, smooth, bclumpy, sclumpy
  timeWindow = 250.0 ; consider accretion over this past time range (Myr)
  tVirInd    = 0 ; plot Tmax/Tvir=1
  
  xrange = [10.0,12.0]
  yrange = [0.01,20.0]
  yrange_stars = 0.1*yrange
  
  lines   = [1,0,2] ; tvircur,tviracc,const5.5
  thicks  = [4,6,8] ; 128,256,512
  sK      = 3 ; smoothing kernel size  
  cInd    = 1 ; color index
  
  ; load
  sPg = simParams(res=256,run='gadget',redshift=2.0)
  sPa = simParams(res=256,run='tracer',redshift=2.0) ; f=-1 use velocity tracers

  arG = haloMassBinValues(sP=sPg,sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow)
  arA = haloMassBinValues(sP=sPa,sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow)
  
  ; plot (1) - both/gmem, variable*const
  pName = sPg.plotPath + 'accRate1.'+accMode+'.comp.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(timeWindow)+'.eps'
  start_PS, pName, /extrabig
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot both (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for j=0,n_elements(arG.TcutVals)-1 do $
      cgPlot,arG.logMassBinCen,smooth(arG.hotMedian.both_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(arA.TcutVals)-1 do $
      cgPlot,arA.logMassBinCen,smooth(arA.hotMedian.both_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; ur: cold both (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for j=0,n_elements(arG.TcutVals)-1 do $
      cgPlot,arG.logMassBinCen,smooth(arG.coldMedian.both_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(arA.TcutVals)-1 do $
      cgPlot,arA.logMassBinCen,smooth(arA.coldMedian.both_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; ll: hot gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for j=0,n_elements(arG.TcutVals)-1 do $
      cgPlot,arG.logMassBinCen,smooth(arG.hotMedian.gmem_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(arA.TcutVals)-1 do $
      cgPlot,arA.logMassBinCen,smooth(arA.hotMedian.gmem_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot

    ; lr: cold gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for j=0,n_elements(arG.TcutVals)-1 do $
      cgPlot,arG.logMassBinCen,smooth(arG.coldMedian.gmem_const[j,*],sK,/nan),color=sPg.colorsG[1],line=j,/overplot
    for j=0,n_elements(arA.TcutVals)-1 do $
      cgPlot,arA.logMassBinCen,smooth(arA.coldMedian.gmem_const[j,*],sK,/nan),color=sPa.colorsA[1],line=j,/overplot

    ; legend
    strings = textoidl("T_{max} / T_{c} = ")+string(arG.TcutVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(arG.TcutVals)),box=0,linesize=0.25,position=[10.05,13.0]

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
  pName = sPg.plotPath + 'accRate2.'+accMode+'.comp.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(timeWindow)+'.eps'
  start_PS, pName, /extrabig
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot both (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for j=0,n_elements(arG.TvirVals)-1 do $
      cgPlot,arG.logMassBinCen,smooth(arG.hotMedian.both_tviracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(arA.TvirVals)-1 do $
      cgPlot,arA.logMassBinCen,smooth(arA.hotMedian.both_tviracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; ur: cold both (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for j=0,n_elements(arG.TvirVals)-1 do $
      cgPlot,arG.logMassBinCen,smooth(arG.coldMedian.both_tviracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(arA.TvirVals)-1 do $
      cgPlot,arA.logMassBinCen,smooth(arA.coldMedian.both_tviracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; ll: hot gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for j=0,n_elements(arG.TvirVals)-1 do $
      cgPlot,arG.logMassBinCen,smooth(arG.hotMedian.gmem_tviracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(arA.TvirVals)-1 do $
      cgPlot,arA.logMassBinCen,smooth(arA.hotMedian.gmem_tviracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    ; lr: cold gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for j=0,n_elements(arG.TvirVals)-1 do $
      cgPlot,arG.logMassBinCen,smooth(arG.coldMedian.gmem_tviracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(arA.TvirVals)-1 do $
      cgPlot,arA.logMassBinCen,smooth(arA.coldMedian.gmem_tviracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot

    ; legend
    strings = textoidl("T_{max} / T_{vir} < ")+string(arG.TvirVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(arG.TvirVals)),box=0,linesize=0.25,position=[10.05,13.0],charsize=!p.charsize-0.1

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
  pName = sPg.plotPath + 'accRate3.'+accMode+'.comp.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(timeWindow)+'.eps'
  start_PS, pName, /extrabig
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot gal (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for j=0,n_elements(arG.TcutVals)-1 do $
      cgPlot,arG.logMassBinCen,smooth(arG.hotMedian.gal_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(arA.TcutVals)-1 do $
      cgPlot,arA.logMassBinCen,smooth(arA.hotMedian.gal_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; ur: cold gal (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for j=0,n_elements(arG.TcutVals)-1 do $
      cgPlot,arG.logMassBinCen,smooth(arG.coldMedian.gal_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(arA.TcutVals)-1 do $
      cgPlot,arA.logMassBinCen,smooth(arA.coldMedian.gal_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; ll: hot stars (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_stars,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for j=0,n_elements(arG.TcutVals)-1 do $
      cgPlot,arG.logMassBinCen,smooth(arG.hotMedian.stars_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(arA.TcutVals)-1 do $
      cgPlot,arA.logMassBinCen,smooth(arA.hotMedian.stars_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot

    ; lr: cold stars (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_stars,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for j=0,n_elements(arG.TcutVals)-1 do $
      cgPlot,arG.logMassBinCen,smooth(arG.coldMedian.stars_const[j,*],sK,/nan),color=sPg.colorsG[1],line=j,/overplot
    for j=0,n_elements(arA.TcutVals)-1 do $
      cgPlot,arA.logMassBinCen,smooth(arA.coldMedian.stars_const[j,*],sK,/nan),color=sPa.colorsA[1],line=j,/overplot

    ; legend
    strings = textoidl("T_{max} / T_{c} = ")+string(arG.TcutVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(arG.TcutVals)),box=0,linesize=0.25,position=[10.05,13.0]

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
  pName = sPg.plotPath + 'accRate4.'+accMode+'.comp.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(timeWindow)+'.eps'
  start_PS, pName, /extrabig
    
    x0 = 0.13 & x1 = 0.53 & x2 = 0.93
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1]  ) ; lr
   
    ; ul: hot gal (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    for j=0,n_elements(arG.TvirVals)-1 do $
      cgPlot,arG.logMassBinCen,smooth(arG.hotMedian.gal_tviracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(arA.TvirVals)-1 do $
      cgPlot,arA.logMassBinCen,smooth(arA.hotMedian.gal_tviracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; ur: cold gal (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    for j=0,n_elements(arG.TvirVals)-1 do $
      cgPlot,arG.logMassBinCen,smooth(arG.coldMedian.gal_tviracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(arA.TvirVals)-1 do $
      cgPlot,arA.logMassBinCen,smooth(arA.coldMedian.gal_tviracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; ll: hot stars (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_stars,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for j=0,n_elements(arG.TvirVals)-1 do $
      cgPlot,arG.logMassBinCen,smooth(arG.hotMedian.stars_tviracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(arA.TvirVals)-1 do $
      cgPlot,arA.logMassBinCen,smooth(arA.hotMedian.stars_tviracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot

    ; lr: cold gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_stars,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for j=0,n_elements(arG.TvirVals)-1 do $
      cgPlot,arG.logMassBinCen,smooth(arG.coldMedian.stars_tviracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(arA.TvirVals)-1 do $
      cgPlot,arA.logMassBinCen,smooth(arA.coldMedian.stars_tviracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot

    ; legend
    strings = textoidl("T_{max} / T_{vir} < ")+string(arG.TvirVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(arG.TvirVals)),box=0,linesize=0.25,position=[10.05,13.0],charsize=!p.charsize-0.1

    ; labels
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Galaxy (Stars Only)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Galaxy (Gas Only)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,"Hot",alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,"Cold",alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS  
  
  ; plot (5) - median lines (gadget vs arepo) gal/gmem separated with specific inset
  if 0 then begin
  pName = sPg.plotPath + 'accRateSp2.'+accMode+'.comp.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(timeWindow)+'.eps'
  start_PS, pName, xs = 7.5, ys = 10.0
    
    x0 = 0.15 & x1 = 0.95
    y0 = 0.1 & y1 = 0.525 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; uc
                [x0,y0,x1,y1]  ) ; lc
                
    possub = list( [0.615,0.57,0.895,0.74] ,$ ; upper sub
                   [0.615,0.14,0.895,0.31]  ) ; lower sub
   
    ; gal
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    ; gadget (gal)
    cgPlot,arG.logMassBinCen,smooth(arG.coldMedian.gal_tViracc[tVirInd,*],sK,/nan),color=sPg.colorsG[cInd],line=lines[0],/overplot
    cgPlot,arG.logMassBinCen,smooth(arG.hotMedian.gal_tViracc[tVirInd,*],sK,/nan),color=sPg.colorsG[cInd],line=lines[2],/overplot
    cgPlot,arG.logMassBinCen,smooth(arG.totalMedian.gal,sK,/nan),color=sPg.colorsG[cInd],line=lines[1],/overplot
   
    ; arepo (gal)
    cgPlot,arA.logMassBinCen,smooth(arA.coldMedian.gal_tViracc[tVirInd,*],sK,/nan),color=sPa.colorsA[cInd],line=lines[0],/overplot
    cgPlot,arA.logMassBinCen,smooth(arA.hotMedian.gal_tViracc[tVirInd,*],sK,/nan),color=sPa.colorsA[cInd],line=lines[2],/overplot
    cgPlot,arA.logMassBinCen,smooth(arA.totalMedian.gal,sK,/nan),color=sPa.colorsA[cInd],line=lines[1],/overplot
    
    ; legend
    legend,[textoidl("T_{max} < 1.0 T_{vir,acc}"),textoidl("T_{max} > 1.0 T_{vir,acc}"),'total'],$
      linestyle=[lines[0],lines[2],lines[1]],box=0,linesize=0.25,/top,/left,spacing=!p.charsize+0.5
    
    ; SUBPLOT: specific (gal)
    polyfill,[possub[0,0]-0.03,possub[0,0]-0.03,possub[0,2]+0.002,possub[0,2]+0.002,possub[0,0]-0.03],$
             [possub[0,1],possub[0,3],possub[0,3],possub[0,1],possub[0,1]],$
             /normal,color=cgColor('white')
             
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=[-13.5,-10.0],/xs,/ys,ytitle="",xtitle="",$
      pos=possub[0],/noerase,charsize=!p.charsize-0.4,ytickv=[-10,-11,-12,-13],yticks=3

    ; gadget (gal)
    cgPlot,arG.logMassBinCen,smooth(alog10(arG.coldMedian.gal_tViracc[tVirInd,*]/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[0],/overplot
    cgPlot,arG.logMassBinCen,smooth(alog10(arG.hotMedian.gal_tViracc[tVirInd,*]/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[2],/overplot
    cgPlot,arG.logMassBinCen,smooth(alog10(arG.totalMedian.gal/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[1],/overplot
   
    ; arepo (gal)
    cgPlot,arA.logMassBinCen,smooth(alog10(arA.coldMedian.gal_tViracc[tVirInd,*]/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[0],/overplot
    cgPlot,arA.logMassBinCen,smooth(alog10(arA.hotMedian.gal_tViracc[tVirInd,*]/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[2],/overplot
    cgPlot,arA.logMassBinCen,smooth(alog10(arA.totalMedian.gal/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[1],/overplot
    
    ; gmem
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),pos=pos[1],/noerase
    
    ; gadget (gmem)
    cgPlot,arG.logMassBinCen,smooth(arG.coldMedian.gmem_tViracc[tVirInd,*],sK,/nan),color=sPg.colorsG[cInd],line=lines[0],/overplot
    cgPlot,arG.logMassBinCen,smooth(arG.hotMedian.gmem_tViracc[tVirInd,*],sK,/nan),color=sPg.colorsG[cInd],line=lines[2],/overplot
    cgPlot,arG.logMassBinCen,smooth(arG.totalMedian.gmem,sK,/nan),color=sPg.colorsG[cInd],line=lines[1],/overplot
    
    ; arepo (gmem)
    cgPlot,arA.logMassBinCen,smooth(arA.coldMedian.gmem_tViracc[tVirInd,*],sK,/nan),color=sPa.colorsA[cInd],line=lines[0],/overplot
    cgPlot,arA.logMassBinCen,smooth(arA.hotMedian.gmem_tViracc[tVirInd,*],sK,/nan),color=sPa.colorsA[cInd],line=lines[2],/overplot
    cgPlot,arA.logMassBinCen,smooth(arA.totalMedian.gmem,sK,/nan),color=sPa.colorsA[cInd],line=lines[1],/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],$
      box=0,/top,/left,spacing=!p.charsize+0.5
    
    ; SUBPLOT: specific (gmem)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=[-12.0,-10.5],/xs,/ys,ytitle="",xtitle="",$
      pos=possub[1],/noerase,charsize=!p.charsize-0.4,ytickv=[-10,-11,-12],yticks=2
    
    ; gadget (gmem)
    cgPlot,arG.logMassBinCen,smooth(alog10(arG.coldMedian.gmem_tViracc[tVirInd,*]/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[0],/overplot
    cgPlot,arG.logMassBinCen,smooth(alog10(arG.hotMedian.gmem_tViracc[tVirInd,*]/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[2],/overplot
    cgPlot,arG.logMassBinCen,smooth(alog10(arG.totalMedian.gmem/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[1],/overplot
    
    ; arepo (gmem)
    cgPlot,arA.logMassBinCen,smooth(alog10(arA.coldMedian.gmem_tViracc[tVirInd,*]/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[0],/overplot
    cgPlot,arA.logMassBinCen,smooth(alog10(arA.hotMedian.gmem_tViracc[tVirInd,*]/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[2],/overplot
    cgPlot,arA.logMassBinCen,smooth(alog10(arA.totalMedian.gmem/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[1],/overplot
    
    ; labels
    cgText,0.05,0.5,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1+0.02,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x1+0.02,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    
  end_PS
  endif ;0  
  
end

; plotColdFracVsHaloMass(): plot the "cold mode fraction" vs halo mass in a few different ways
;                           using the mergerTreeSubset with accretionTimes and accretionModes

pro plotColdFracVsHaloMass

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect   = 'pri'
  accMode    = 'all'
  redshift   = 2.0
  timeWindow = 250.0 ; Myr
 
  lines   = [1,0,2] ; tvircur,tviracc,const5.5
  thicks  = [4,6,8] ; 128,256,512
  sK      = 1 ; smoothing kernel size
  cInd    = 1 ; color index

  virInd = 0 ;1.0 Tvir
  conInd = 1 ;5.5 const
  
  xrange = [10.0,12.0]
  yrange = [0.0,1.025]
 
  ; 512
  print,'todo change to 512'
  sPg = simParams(res=256,run='gadget',redshift=redshift)
  sPa = simParams(res=256,run='tracer',redshift=redshift) ; f=-1 use velocity tracers
  
  cfG = haloMassBinValues(sP=sPg,sgSelect=sgSelect,timeWindow=timeWindow,accMode=accMode)
  cfA = haloMassBinValues(sP=sPa,sgSelect=sgSelect,timeWindow=timeWindow,accMode=accMode)

  ; load 256
  cfG_256 = haloMassBinValues(sP=simParams(res=256,run='gadget',redshift=redshift),sgSelect=sgSelect,timeWindow=timeWindow,accMode=accMode)
  cfA_256 = haloMassBinValues(sP=simParams(res=256,run='tracer',redshift=redshift),sgSelect=sgSelect,timeWindow=timeWindow,accMode=accMode)
  
  ; load 128
  cfG_128 = haloMassBinValues(sP=simParams(res=128,run='gadget',redshift=redshift),sgSelect=sgSelect,timeWindow=timeWindow,accMode=accMode)
  cfA_128 = haloMassBinValues(sP=simParams(res=128,run='tracer',redshift=redshift),sgSelect=sgSelect,timeWindow=timeWindow,accMode=accMode)

  ; plot (1) - gal median lines (Tmax) 1*tvircur,1*tviracc,tcut=5.5
  start_PS, sPg.plotPath + 'coldFrac1.gal.'+accMode+'.comp.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.res)+'_'+str(sPg.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]")
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; gadget gal 128
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.fracMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[2],thick=thicks[0],/overplot
    
    ; arepo gal 128
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.fracMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[2],thick=thicks[0],/overplot
        
    ; gadget gal 256
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.fracMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],thick=thicks[1],/overplot
    
    ; arepo gal 256
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.fracMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],thick=thicks[1],/overplot
    
    ; gadget gal 512
    cgPlot,cfG.logMassBinCen,smooth(cfG.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,cfG.logMassBinCen,smooth(cfG.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,cfG.logMassBinCen,smooth(cfG.fracMedian.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[2],thick=thicks[2],/overplot
    
    ; arepo gal 512
    cgPlot,cfA.logMassBinCen,smooth(cfA.fracMedian.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,cfA.logMassBinCen,smooth(cfA.fracMedian.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,cfA.logMassBinCen,smooth(cfA.fracMedian.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[2],thick=thicks[2],/overplot
    
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
  start_PS, sPg.plotPath + 'coldFrac1.gmem.'+accMode+'.comp.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.res)+'_'+str(sPg.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]")
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; gadget gal 128
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[2],thick=thicks[0],/overplot
    
    ; arepo gal 128
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[2],thick=thicks[0],/overplot
        
    ; gadget gal 256
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],thick=thicks[1],/overplot
    
    ; arepo gal 256
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],thick=thicks[1],/overplot
    
    ; gadget gal 512
    cgPlot,cfG.logMassBinCen,smooth(cfG.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,cfG.logMassBinCen,smooth(cfG.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,cfG.logMassBinCen,smooth(cfG.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[2],thick=thicks[2],/overplot
    
    ; arepo gal 512
    cgPlot,cfA.logMassBinCen,smooth(cfA.fracMedian.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,cfA.logMassBinCen,smooth(cfA.fracMedian.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,cfA.logMassBinCen,smooth(cfA.fracMedian.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[2],thick=thicks[2],/overplot
    
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
  start_PS, sPg.plotPath + 'coldFrac1.stars.'+accMode+'.comp.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.res)+'_'+str(sPg.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]")
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; gadget gal 128
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.fracMedian.stars_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.fracMedian.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[2],thick=thicks[0],/overplot
    
    ; arepo gal 128
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.fracMedian.stars_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.fracMedian.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[2],thick=thicks[0],/overplot
        
    ; gadget gal 256
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.fracMedian.stars_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.fracMedian.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],thick=thicks[1],/overplot
    
    ; arepo gal 256
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.fracMedian.stars_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.fracMedian.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],thick=thicks[1],/overplot
    
    ; gadget gal 512
    cgPlot,cfG.logMassBinCen,smooth(cfG.fracMedian.stars_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,cfG.logMassBinCen,smooth(cfG.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,cfG.logMassBinCen,smooth(cfG.fracMedian.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[2],thick=thicks[2],/overplot
    
    ; arepo gal 512
    cgPlot,cfA.logMassBinCen,smooth(cfA.fracMedian.stars_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,cfA.logMassBinCen,smooth(cfA.fracMedian.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,cfA.logMassBinCen,smooth(cfA.fracMedian.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[2],thick=thicks[2],/overplot
    
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
  start_PS, sPg.plotPath + 'coldFrac1.both.'+accMode+'.comp.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.res)+'_'+str(sPg.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]")
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; gadget gal 128
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.fracMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.fracMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[2],thick=thicks[0],/overplot
    
    ; arepo gal 128
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.fracMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.fracMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[2],thick=thicks[0],/overplot
        
    ; gadget gal 256
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.fracMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.fracMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],thick=thicks[1],/overplot
    
    ; arepo gal 256
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.fracMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.fracMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],thick=thicks[1],/overplot
    
    ; gadget gal 512
    cgPlot,cfG.logMassBinCen,smooth(cfG.fracMedian.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,cfG.logMassBinCen,smooth(cfG.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,cfG.logMassBinCen,smooth(cfG.fracMedian.both_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[2],thick=thicks[2],/overplot
    
    ; arepo gal 512
    cgPlot,cfA.logMassBinCen,smooth(cfA.fracMedian.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,cfA.logMassBinCen,smooth(cfA.fracMedian.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,cfA.logMassBinCen,smooth(cfA.fracMedian.both_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[2],thick=thicks[2],/overplot
    
    ; style legend
    strings = [textoidl("T_{max} < 1.0 T_{vir,cur}"),$
               textoidl("T_{max} < 1.0 T_{vir,acc}"),$
               textoidl("T_{max} < T_{c} = 10^{5.5} K")]
    legend,strings,linestyle=[1,0,2],box=0,linesize=0.25,spacing=!p.charsize+0.5,position=[10.05,0.3],charsize=!p.charsize-0.1
    
    ; color/thick legend
    strings = [textoidl(['128^3','256^3','512^3']+' GADGET'),textoidl(['128^3','256^3','512^3']+' AREPO')]
    legend,strings,textcolors=[sPg.colorsG,sPa.colorsA],box=0,position=[11.45,0.82],charsize=!p.charsize-0.3
    
  end_PS
  
  ; plot (2) - 2x2 of (gal,gmem)x(const,tvir)
  pName = sPg.plotPath + 'coldFrac4.'+accMode+'.comp.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'.eps'
  start_PS, pName, /extrabig
    
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
    
    for j=0,n_elements(cfG.TvirVals)-1 do $
      cgPlot,cfG.logMassBinCen,smooth(cfG.fracMedian.gal_tViracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(cfA.TvirVals)-1 do $
      cgPlot,cfA.logMassBinCen,smooth(cfA.fracMedian.gal_tViracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; top left axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    cgAxis,/xaxis,xrange=tvirrange,/xs
    
    ; ur: gal const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    ; gadget and arepo gal
    for j=0,n_elements(cfG.TcutVals)-1 do $
      cgPlot,cfG.logMassBinCen,smooth(cfG.fracMedian.gal_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(cfG.TcutVals)-1 do $
      cgPlot,cfA.logMassBinCen,smooth(cfA.fracMedian.gal_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
      
    ; top right axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    tvirrange[1] = round(10*tvirrange[1])/10.0
    cgAxis,/xaxis,xrange=tvirrange,/xs

    ; ll: gmem tvir (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for j=0,n_elements(cfG.TvirVals)-1 do $
      cgPlot,cfG.logMassBinCen,smooth(cfG.fracMedian.gmem_tViracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(cfG.TvirVals)-1 do $
      cgPlot,cfA.logMassBinCen,smooth(cfA.fracMedian.gmem_tViracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
      
    ; legend
    strings = textoidl("T_{max} / T_{vir} < ")+string(cfG.TvirVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(cfG.TvirVals)),box=0,linesize=0.4,$
      position=[10.95,0.95],charsize=!p.charsize-0.2

    ; lr: gmem const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for j=0,n_elements(cfG.TcutVals)-1 do $
      cgPlot,cfG.logMassBinCen,smooth(cfG.fracMedian.gmem_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(cfG.TcutVals)-1 do $
      cgPlot,cfA.logMassBinCen,smooth(cfA.fracMedian.gmem_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot

    ; legend
    strings = textoidl("T_{max} < T_{c} = ")+string(cfG.TcutVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(cfG.TcutVals)),box=0,linesize=0.4,$
      position=[10.95,0.95],charsize=!p.charsize-0.2

    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/bottom,/left

    ; labels
    cgText,0.05,y1,"Cold Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.01,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    cgText,x1,y2+0.06,textoidl("T_{vir} [_{ }log K_{ }] (z=2)"),alignment=0.5,/normal
    
    cgText,x2+0.02,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.02,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x1,x2])+0.05,y2+0.06,"Variable Constant "+textoidl("T_c"),alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x0,x1])-0.05,y2+0.06,"Variable "+textoidl("T_{max} / T_{vir,acc}"),alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
end

; plotByMode(): plot cold fractions and accretion rates separated out into the four modes and also
;               a view with decomposed galaxy=gal+stars (one resolution only)

pro plotByMode, res=res

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect = 'pri'
  redshift = 2.0
  if ~keyword_set(res) then message,'need res'
  
  timeWindow = 250.0
  ;timeWindows = list(250.0,1000.0,'all') ; Myr
  ;foreach timeWindow,timeWindows do begin
  
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

  ;endforeach
  ;stop

  sK     = 3 ; smoothing kernel size
  cInd   = 1 ; color index
  virInd = 0 ; 1.0 Tvir
  conInd = 1 ; 5.5 const
  
  xrange = [10.0,12.0]
  yrange = [0.0,1.025]

  ; plot (1) - 2x2 of (gal,gmem)x(const,tvir)
  pName = sPg.plotPath + 'coldFracByMode.comp1.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'.eps'
  start_PS, pName, /extrabig
    
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
    cgText,mean([x0,x1])-0.05,y2+0.06,textoidl("T_{max} / T_{vir,acc} < 1.0"),alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; plot (2) - 2x2 of (gal,gmem)x(const,tvir)
  pName = sPg.plotPath + 'coldFracByMode.comp2.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'.eps'
  start_PS, pName, /extrabig
    
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
    cgText,mean([x0,x1])-0.05,y2+0.06,textoidl("T_{max} / T_{vir,acc} < 1.0"),alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS

  yrange = [0.01,20.0]
  yrange_stars = 0.1*yrange

  ; plot (3) - accretion rate decomposed galaxy=gas+stars variable*tvir
  pName = sPg.plotPath + 'accRateByMode.comp1.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(timeWindow)+'.eps'
  start_PS, pName, /extrabig
    
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
      position=[11.2,0.95],charsize=!p.charsize-0.2
      
    ; labels
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Galaxy (Stars Only)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Galaxy (Gas Only)",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,"Hot",alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,"Cold",alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; plot (4) - accretion rate galaxy/halo atmosphere variable*tvir
  pName = sPg.plotPath + 'accRateByMode.comp2.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(timeWindow)+'.eps'
  start_PS, pName, /extrabig
    
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
      position=[11.2,0.95],charsize=!p.charsize-0.2
      
    ; labels
    cgText,0.05,y1,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1,0.05,textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),alignment=0.5,/normal
    
    cgText,x2+0.015,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x2+0.015,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,mean([x0,x1]),y2+0.015,"Hot",alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,"Cold",alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  yrange = [0.0,1.0]

  ; plot (5) - fraction, accretion rate decomposed galaxy=gas+stars variable*tvir
  pName = sPg.plotPath + 'accRateFracByMode1.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(timeWindow)+'.eps'
  start_PS, pName, /extrabig
    
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
    cgText,mean([x0,x1]),y2+0.015,"Hot",alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,"Cold",alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
  ; plot (6) - fraction, accretion rate galaxy/halo atmosphere variable*tvir
  pName = sPg.plotPath + 'accRateFracByMode2.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(timeWindow)+'.eps'
  start_PS, pName, /extrabig
    
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
    cgText,mean([x0,x1]),y2+0.015,"Hot",alignment=0.5,color=cgColor('dark gray'),/normal
    cgText,mean([x1,x2]),y2+0.015,"Cold",alignment=0.5,color=cgColor('dark gray'),/normal
  end_PS
  
end
