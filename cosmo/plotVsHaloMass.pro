; plotVsHaloMass.pro
; gas accretion project - plots as a function of halo mass
; dnelson aug.2012

; plotAccRateVsHaloMass(): plot the "cold mode fraction" vs halo mass in a few different ways
;                           using the mergerTreeSubset with accretionTimes and accretionModes

pro plotAccRateVsHaloMass

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect = 'pri'
  accMode  = 'all' ; accretion mode: all, smooth, bclumpy, sclumpy
  
  timeWindow = 1000.0 ; consider accretion over this past time range (Myr)
  tVirInd = 0 ; plot Tmax/Tvir=1
  
  print,'gal not changed to both yet for accrate'
  
  sPg = simParams(res=256,run='gadget',redshift=2.0)
  sPa = simParams(res=256,run='tracer',redshift=2.0) ; f=-1 use velocity tracers

  arG = haloMassBinAccRate(sP=sPg,sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow)
  arA = haloMassBinAccRate(sP=sPa,sgSelect=sgSelect,accMode=accMode,timeWindow=timeWindow)
  
  ; minimum halo mass to plot
  w = where(arG.logMassBinCen gt 8.0)
  
  xrange = [10.0,12.0]
  yrange = [0.01,20.0]
  
  lines   = [1,0,2] ; tvircur,tviracc,const5.5
  thicks  = [4,6,8] ; 128,256,512
  
  sK = 3 ; smoothing kernel size  
  cInd = 1
  
  ; plot (1c) - median lines (gadget vs arepo) gal/gmem separated with specific inset
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
    cgPlot,arG.logMassBinCen[w],smooth(arG.coldMedian.gal_tViracc[tVirInd,w],sK,/nan),color=sPg.colorsG[cInd],line=lines[0],/overplot
    cgPlot,arG.logMassBinCen[w],smooth(arG.hotMedian.gal_tViracc[tVirInd,w],sK,/nan),color=sPg.colorsG[cInd],line=lines[2],/overplot
    cgPlot,arG.logMassBinCen[w],smooth(arG.totalMedian.gal[w],sK,/nan),color=sPg.colorsG[cInd],line=lines[1],/overplot
   
    ; arepo (gal)
    cgPlot,arA.logMassBinCen[w],smooth(arA.coldMedian.gal_tViracc[tVirInd,w],sK,/nan),color=sPa.colorsA[cInd],line=lines[0],/overplot
    cgPlot,arA.logMassBinCen[w],smooth(arA.hotMedian.gal_tViracc[tVirInd,w],sK,/nan),color=sPa.colorsA[cInd],line=lines[2],/overplot
    cgPlot,arA.logMassBinCen[w],smooth(arA.totalMedian.gal[w],sK,/nan),color=sPa.colorsA[cInd],line=lines[1],/overplot
    
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
    cgPlot,arG.logMassBinCen[w],smooth(alog10(arG.coldMedian.gal_tViracc[tVirInd,w]/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[0],/overplot
    cgPlot,arG.logMassBinCen[w],smooth(alog10(arG.hotMedian.gal_tViracc[tVirInd,w]/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[2],/overplot
    cgPlot,arG.logMassBinCen[w],smooth(alog10(arG.totalMedian.gal[w]/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[1],/overplot
   
    ; arepo (gal)
    cgPlot,arA.logMassBinCen[w],smooth(alog10(arA.coldMedian.gal_tViracc[tVirInd,w]/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[0],/overplot
    cgPlot,arA.logMassBinCen[w],smooth(alog10(arA.hotMedian.gal_tViracc[tVirInd,w]/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[2],/overplot
    cgPlot,arA.logMassBinCen[w],smooth(alog10(arA.totalMedian.gal[w]/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[1],/overplot
    
    ; gmem
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),pos=pos[1],/noerase
    
    ; gadget (gmem)
    cgPlot,arG.logMassBinCen[w],smooth(arG.coldMedian.gmem_tViracc[tVirInd,w],sK,/nan),color=sPg.colorsG[cInd],line=lines[0],/overplot
    cgPlot,arG.logMassBinCen[w],smooth(arG.hotMedian.gmem_tViracc[tVirInd,w],sK,/nan),color=sPg.colorsG[cInd],line=lines[2],/overplot
    cgPlot,arG.logMassBinCen[w],smooth(arG.totalMedian.gmem[w],sK,/nan),color=sPg.colorsG[cInd],line=lines[1],/overplot
    
    ; arepo (gmem)
    cgPlot,arA.logMassBinCen[w],smooth(arA.coldMedian.gmem_tViracc[tVirInd,w],sK,/nan),color=sPa.colorsA[cInd],line=lines[0],/overplot
    cgPlot,arA.logMassBinCen[w],smooth(arA.hotMedian.gmem_tViracc[tVirInd,w],sK,/nan),color=sPa.colorsA[cInd],line=lines[2],/overplot
    cgPlot,arA.logMassBinCen[w],smooth(arA.totalMedian.gmem[w],sK,/nan),color=sPa.colorsA[cInd],line=lines[1],/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],$
      box=0,/top,/left,spacing=!p.charsize+0.5
    
    ; SUBPLOT: specific (gmem)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=[-12.0,-10.5],/xs,/ys,ytitle="",xtitle="",$
      pos=possub[1],/noerase,charsize=!p.charsize-0.4,ytickv=[-10,-11,-12],yticks=2
    
    ; gadget (gmem)
    cgPlot,arG.logMassBinCen[w],smooth(alog10(arG.coldMedian.gmem_tViracc[tVirInd,w]/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[0],/overplot
    cgPlot,arG.logMassBinCen[w],smooth(alog10(arG.hotMedian.gmem_tViracc[tVirInd,w]/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[2],/overplot
    cgPlot,arG.logMassBinCen[w],smooth(alog10(arG.totalMedian.gmem[w]/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[1],/overplot
    
    ; arepo (gmem)
    cgPlot,arA.logMassBinCen[w],smooth(alog10(arA.coldMedian.gmem_tViracc[tVirInd,w]/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[0],/overplot
    cgPlot,arA.logMassBinCen[w],smooth(alog10(arA.hotMedian.gmem_tViracc[tVirInd,w]/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[2],/overplot
    cgPlot,arA.logMassBinCen[w],smooth(alog10(arA.totalMedian.gmem[w]/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[1],/overplot
    
    ; labels
    cgText,0.05,0.5,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1+0.02,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x1+0.02,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    
  end_PS
  
  ; plot (2) - variable*const
  pName = sPg.plotPath + 'accRate2.'+accMode+'.comp.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
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
  
  ; plot (2) - variable*tvir
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
  
  ; plot (1) - median lines (gadget vs arepo) gal/gmem separated
  if 0 then begin
  pName = sPg.plotPath + 'accRate.'+accMode+'.comp.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(timeWindow)+'.eps'
  start_PS, pName, xs = 7.5, ys = 10.0
    
    x0 = 0.15 & x1 = 0.95
    y0 = 0.1 & y1 = 0.525 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; uc
                [x0,y0,x1,y1]  ) ; lc
   
    ; gal
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    ; gadget (gal)
    cgPlot,arG.logMassBinCen[w],smooth(arG.coldMedian.gal_tViracc[tVirInd,w],sK,/nan),color=sPg.colorsG[cInd],line=lines[0],/overplot
    cgPlot,arG.logMassBinCen[w],smooth(arG.hotMedian.gal_tViracc[tVirInd,w],sK,/nan),color=sPg.colorsG[cInd],line=lines[2],/overplot
    cgPlot,arG.logMassBinCen[w],smooth(arG.totalMedian.gal[w],sK,/nan),color=sPg.colorsG[cInd],line=lines[1],/overplot
   
    ; arepo (gal)
    cgPlot,arA.logMassBinCen[w],smooth(arA.coldMedian.gal_tViracc[tVirInd,w],sK,/nan),color=sPa.colorsA[cInd],line=lines[0],/overplot
    cgPlot,arA.logMassBinCen[w],smooth(arA.hotMedian.gal_tViracc[tVirInd,w],sK,/nan),color=sPa.colorsA[cInd],line=lines[2],/overplot
    cgPlot,arA.logMassBinCen[w],smooth(arA.totalMedian.gal[w],sK,/nan),color=sPa.colorsA[cInd],line=lines[1],/overplot
    
    ; legend
    legend,[textoidl("T_{max} < 1.0 T_{vir,acc}"),textoidl("T_{max} > 1.0 T_{vir,acc}"),'total'],$
      linestyle=[lines[0],lines[2],lines[1]],box=0,linesize=0.25,/top,/left,spacing=!p.charsize+0.5
    
    ; gmem
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),pos=pos[1],/noerase
    
    ; gadget (gmem)
    cgPlot,arG.logMassBinCen[w],smooth(arG.coldMedian.gmem_tViracc[tVirInd,w],sK,/nan),color=sPg.colorsG[cInd],line=lines[0],/overplot
    cgPlot,arG.logMassBinCen[w],smooth(arG.hotMedian.gmem_tViracc[tVirInd,w],sK,/nan),color=sPg.colorsG[cInd],line=lines[2],/overplot
    cgPlot,arG.logMassBinCen[w],smooth(arG.totalMedian.gmem[w],sK,/nan),color=sPg.colorsG[cInd],line=lines[1],/overplot
    
    ; arepo (gmem)
    cgPlot,arA.logMassBinCen[w],smooth(arA.coldMedian.gmem_tViracc[tVirInd,w],sK,/nan),color=sPa.colorsA[cInd],line=lines[0],/overplot
    cgPlot,arA.logMassBinCen[w],smooth(arA.hotMedian.gmem_tViracc[tVirInd,w],sK,/nan),color=sPa.colorsA[cInd],line=lines[2],/overplot
    cgPlot,arA.logMassBinCen[w],smooth(arA.totalMedian.gmem[w],sK,/nan),color=sPa.colorsA[cInd],line=lines[1],/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPa.colorsG[cInd],sPg.colorsA[cInd]],box=0,/top,/left,spacing=!p.charsize+0.5
    
    ; labels
    cgText,0.05,0.5,"Gas Accretion Rate "+textoidl("[_{ }h^{-1} M_{sun } yr^{-1 }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,x1+0.02,mean([y0,y1]),"Halo Atmosphere",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    cgText,x1+0.02,mean([y1,y2]),"Central Galaxy",alignment=0.5,color=cgColor('dark gray'),orientation=-90.0,/normal
    
  end_PS

  ; plot (1b) - specific median lines (gadget vs arepo) gal/gmem separated
  pName = sPg.plotPath + 'accRateSp.'+accMode+'.comp.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
          str(sPg.res)+'_'+str(sPg.snap)+'_tw'+str(timeWindow)+'.eps'
  start_PS, pName, xs = 7.5, ys = 10.0
    
    x0 = 0.15 & x1 = 0.95
    y0 = 0.1 & y1 = 0.525 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; uc
                [x0,y0,x1,y1]  ) ; lc
   
    ; gal
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=[-13.5,-10.0],/xs,/ys,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    ; gadget (gal)
    cgPlot,arG.logMassBinCen[w],smooth(alog10(arG.coldMedian.gal_tViracc[tVirInd,w]/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[0],/overplot
    cgPlot,arG.logMassBinCen[w],smooth(alog10(arG.hotMedian.gal_tViracc[tVirInd,w]/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[2],/overplot
    cgPlot,arG.logMassBinCen[w],smooth(alog10(arG.totalMedian.gal[w]/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[1],/overplot
   
    ; arepo (gal)
    cgPlot,arA.logMassBinCen[w],smooth(alog10(arA.coldMedian.gal_tViracc[tVirInd,w]/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[0],/overplot
    cgPlot,arA.logMassBinCen[w],smooth(alog10(arA.hotMedian.gal_tViracc[tVirInd,w]/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[2],/overplot
    cgPlot,arA.logMassBinCen[w],smooth(alog10(arA.totalMedian.gal[w]/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[1],/overplot
    
    ; legend
    legend,[textoidl("T_{max} < 1.0 T_{vir,acc}"),textoidl("T_{max} > 1.0 T_{vir,acc}"),'total'],$
      linestyle=[lines[0],lines[2],lines[1]],box=0,linesize=0.25,/bottom,/left,spacing=!p.charsize+0.5
    
    ; gmem
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=[-12.0,-10.5],/xs,/ys,$
      ytitle="",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),pos=pos[1],/noerase
    
    ; gadget (gmem)
    cgPlot,arG.logMassBinCen[w],smooth(alog10(arG.coldMedian.gmem_tViracc[tVirInd,w]/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[0],/overplot
    cgPlot,arG.logMassBinCen[w],smooth(alog10(arG.hotMedian.gmem_tViracc[tVirInd,w]/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[2],/overplot
    cgPlot,arG.logMassBinCen[w],smooth(alog10(arG.totalMedian.gmem[w]/10.0^arG.logMassBinCen),sK,/nan),color=sPg.colorsG[cInd],line=lines[1],/overplot
    
    ; arepo (gmem)
    cgPlot,arA.logMassBinCen[w],smooth(alog10(arA.coldMedian.gmem_tViracc[tVirInd,w]/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[0],/overplot
    cgPlot,arA.logMassBinCen[w],smooth(alog10(arA.hotMedian.gmem_tViracc[tVirInd,w]/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[2],/overplot
    cgPlot,arA.logMassBinCen[w],smooth(alog10(arA.totalMedian.gmem[w]/10.0^arG.logMassBinCen),sK,/nan),color=sPa.colorsA[cInd],line=lines[1],/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left,spacing=!p.charsize+0.5
    
    ; labels
    cgText,0.05,0.5,"Specific Gas Accretion Rate "+textoidl("[_{ }h^{-1} yr^{-1 }]"),$
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
  sgSelect = 'pri'
  accMode = 'all'
  redshift = 2.0
 
  sK = 1 ; smoothing kernel size
  cInd = 1
  
  lines   = [1,0,2] ; tvircur,tviracc,const5.5
  thicks  = [4,6,8] ; 128,256,512
  
  virInd = 0 ;1.0 Tvir
  conInd = 1 ;5.5 const
  
  xrange = [10.0,12.0]
  yrange = [0.0,1.025]
 
  ; 512
  print,'todo change to 512'
  sPg = simParams(res=256,run='gadget',redshift=redshift)
  sPa = simParams(res=128,run='tracer',redshift=redshift) ; f=-1 use velocity tracers
  
  ; old
  ;sPg.derivPath = '/n/home07/dnelson/data3/sims.gadget/256_20Mpc/data.files.bak/'
  ;sPa.derivPath = '/n/home07/dnelson/data3/sims.tracers/256_20Mpc/data.files.bak/'

  cfG = haloMassBinColdFracs(sP=sPg,sgSelect=sgSelect,accMode=accMode)
  cfA = haloMassBinColdFracs(sP=sPa,sgSelect=sgSelect,accMode=accMode)

  ; load 256
  cfG_256 = haloMassBinColdFracs(sP=simParams(res=256,run='gadget',redshift=redshift),sgSelect=sgSelect,accMode=accMode)
  cfA_256 = haloMassBinColdFracs(sP=simParams(res=128,run='tracer',redshift=redshift),sgSelect=sgSelect,accMode=accMode)
  
  ; load 128
  cfG_128 = haloMassBinColdFracs(sP=simParams(res=128,run='gadget',redshift=redshift),sgSelect=sgSelect,accMode=accMode)
  cfA_128 = haloMassBinColdFracs(sP=simParams(res=128,run='tracer',redshift=redshift),sgSelect=sgSelect,accMode=accMode)

  ; plot (1) - gal median lines (Tmax) 1*tvircur,1*tviracc,tcut=5.5
  start_PS, sPg.plotPath + 'coldFrac1.gal.'+accMode+'.comp.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.res)+'_'+str(sPg.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]")
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; gadget gal 128
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.medianVals.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.medianVals.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.medianVals.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[2],thick=thicks[0],/overplot
    
    ; arepo gal 128
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.medianVals.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.medianVals.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.medianVals.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[2],thick=thicks[0],/overplot
        
    ; gadget gal 256
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.medianVals.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.medianVals.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.medianVals.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],thick=thicks[1],/overplot
    
    ; arepo gal 256
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.medianVals.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.medianVals.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.medianVals.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],thick=thicks[1],/overplot
    
    ; gadget gal 512
    cgPlot,cfG.logMassBinCen,smooth(cfG.medianVals.gal_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,cfG.logMassBinCen,smooth(cfG.medianVals.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,cfG.logMassBinCen,smooth(cfG.medianVals.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[2],thick=thicks[2],/overplot
    
    ; arepo gal 512
    cgPlot,cfA.logMassBinCen,smooth(cfA.medianVals.gal_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,cfA.logMassBinCen,smooth(cfA.medianVals.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,cfA.logMassBinCen,smooth(cfA.medianVals.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[2],thick=thicks[2],/overplot
    
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
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.medianVals.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.medianVals.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.medianVals.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[2],thick=thicks[0],/overplot
    
    ; arepo gal 128
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.medianVals.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.medianVals.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.medianVals.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[2],thick=thicks[0],/overplot
        
    ; gadget gal 256
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.medianVals.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.medianVals.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.medianVals.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],thick=thicks[1],/overplot
    
    ; arepo gal 256
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.medianVals.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.medianVals.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.medianVals.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],thick=thicks[1],/overplot
    
    ; gadget gal 512
    cgPlot,cfG.logMassBinCen,smooth(cfG.medianVals.gmem_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,cfG.logMassBinCen,smooth(cfG.medianVals.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,cfG.logMassBinCen,smooth(cfG.medianVals.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[2],thick=thicks[2],/overplot
    
    ; arepo gal 512
    cgPlot,cfA.logMassBinCen,smooth(cfA.medianVals.gmem_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,cfA.logMassBinCen,smooth(cfA.medianVals.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,cfA.logMassBinCen,smooth(cfA.medianVals.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[2],thick=thicks[2],/overplot
    
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
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.medianVals.stars_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.medianVals.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.medianVals.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[2],thick=thicks[0],/overplot
    
    ; arepo gal 128
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.medianVals.stars_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.medianVals.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.medianVals.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[2],thick=thicks[0],/overplot
        
    ; gadget gal 256
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.medianVals.stars_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.medianVals.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.medianVals.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],thick=thicks[1],/overplot
    
    ; arepo gal 256
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.medianVals.stars_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.medianVals.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.medianVals.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],thick=thicks[1],/overplot
    
    ; gadget gal 512
    cgPlot,cfG.logMassBinCen,smooth(cfG.medianVals.stars_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,cfG.logMassBinCen,smooth(cfG.medianVals.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,cfG.logMassBinCen,smooth(cfG.medianVals.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[2],thick=thicks[2],/overplot
    
    ; arepo gal 512
    cgPlot,cfA.logMassBinCen,smooth(cfA.medianVals.stars_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,cfA.logMassBinCen,smooth(cfA.medianVals.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,cfA.logMassBinCen,smooth(cfA.medianVals.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[2],thick=thicks[2],/overplot
    
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
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.medianVals.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.medianVals.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,cfG_128.logMassBinCen,smooth(cfG_128.medianVals.both_const[conInd,*],sK,/nan),color=sPg.colorsG[0],line=lines[2],thick=thicks[0],/overplot
    
    ; arepo gal 128
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.medianVals.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[0],thick=thicks[0],/overplot
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.medianVals.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[1],thick=thicks[0],/overplot
    cgPlot,cfA_128.logMassBinCen,smooth(cfA_128.medianVals.both_const[conInd,*],sK,/nan),color=sPa.colorsA[0],line=lines[2],thick=thicks[0],/overplot
        
    ; gadget gal 256
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.medianVals.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.medianVals.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,cfG_256.logMassBinCen,smooth(cfG_256.medianVals.both_const[conInd,*],sK,/nan),color=sPg.colorsG[1],line=lines[2],thick=thicks[1],/overplot
    
    ; arepo gal 256
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.medianVals.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[0],thick=thicks[1],/overplot
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.medianVals.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[1],thick=thicks[1],/overplot
    cgPlot,cfA_256.logMassBinCen,smooth(cfA_256.medianVals.both_const[conInd,*],sK,/nan),color=sPa.colorsA[1],line=lines[2],thick=thicks[1],/overplot
    
    ; gadget gal 512
    cgPlot,cfG.logMassBinCen,smooth(cfG.medianVals.both_tVircur[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,cfG.logMassBinCen,smooth(cfG.medianVals.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,cfG.logMassBinCen,smooth(cfG.medianVals.both_const[conInd,*],sK,/nan),color=sPg.colorsG[2],line=lines[2],thick=thicks[2],/overplot
    
    ; arepo gal 512
    cgPlot,cfA.logMassBinCen,smooth(cfA.medianVals.both_tVircur[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[0],thick=thicks[2],/overplot
    cgPlot,cfA.logMassBinCen,smooth(cfA.medianVals.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[1],thick=thicks[2],/overplot
    cgPlot,cfA.logMassBinCen,smooth(cfA.medianVals.both_const[conInd,*],sK,/nan),color=sPa.colorsA[2],line=lines[2],thick=thicks[2],/overplot
    
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
      cgPlot,cfG.logMassBinCen,smooth(cfG.medianVals.gal_tViracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(cfA.TvirVals)-1 do $
      cgPlot,cfA.logMassBinCen,smooth(cfA.medianVals.gal_tViracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
    
    ; top left axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    cgAxis,/xaxis,xrange=tvirrange,/xs
    
    ; ur: gal const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    ; gadget and arepo gal
    for j=0,n_elements(cfG.TcutVals)-1 do $
      cgPlot,cfG.logMassBinCen,smooth(cfG.medianVals.gal_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(cfG.TcutVals)-1 do $
      cgPlot,cfA.logMassBinCen,smooth(cfA.medianVals.gal_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
      
    ; top right axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    tvirrange[1] = round(10*tvirrange[1])/10.0
    cgAxis,/xaxis,xrange=tvirrange,/xs

    ; ll: gmem tvir (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    for j=0,n_elements(cfG.TvirVals)-1 do $
      cgPlot,cfG.logMassBinCen,smooth(cfG.medianVals.gmem_tViracc[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(cfG.TvirVals)-1 do $
      cgPlot,cfA.logMassBinCen,smooth(cfA.medianVals.gmem_tViracc[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot
      
    ; legend
    strings = textoidl("T_{max} / T_{vir} < ")+string(cfG.TvirVals,format='(f4.1)')
    legend,strings,linestyle=indgen(n_elements(cfG.TvirVals)),box=0,linesize=0.4,$
      position=[10.95,0.95],charsize=!p.charsize-0.2

    ; lr: gmem const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    for j=0,n_elements(cfG.TcutVals)-1 do $
      cgPlot,cfG.logMassBinCen,smooth(cfG.medianVals.gmem_const[j,*],sK,/nan),color=sPg.colorsG[cInd],line=j,/overplot
    for j=0,n_elements(cfG.TcutVals)-1 do $
      cgPlot,cfA.logMassBinCen,smooth(cfA.medianVals.gmem_const[j,*],sK,/nan),color=sPa.colorsA[cInd],line=j,/overplot

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

; ---------------------------------------------------------------------------------------------

; plotByMode(): plot cold fractions and accretion rates separated out into the four modes and also
;               a view with decomposed galaxy=gal+stars (one resolution only)

pro plotByMode

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect = 'pri'
  redshift = 2.0
  res = 256
  
  timeWindow = 1000.0 ; Myr
 
  ; load
  sPg = simParams(res=res,run='gadget',redshift=redshift)
  sPa = simParams(res=res,run='tracer',redshift=redshift)
  
  cfG_all      = haloMassBinColdFracs(sP=sPg,sgSelect=sgSelect,accMode='all')
  cfA_all      = haloMassBinColdFracs(sP=sPa,sgSelect=sgSelect,accMode='all')
  cfG_smooth   = haloMassBinColdFracs(sP=sPg,sgSelect=sgSelect,accMode='smooth')
  cfA_smooth   = haloMassBinColdFracs(sP=sPa,sgSelect=sgSelect,accMode='smooth')
  cfG_clumpy   = haloMassBinColdFracs(sP=sPg,sgSelect=sgSelect,accMode='clumpy')
  cfA_clumpy   = haloMassBinColdFracs(sP=sPa,sgSelect=sgSelect,accMode='clumpy')
  cfG_stripped = haloMassBinColdFracs(sP=sPg,sgSelect=sgSelect,accMode='stripped')
  cfA_stripped = haloMassBinColdFracs(sP=sPa,sgSelect=sgSelect,accMode='stripped')

  arG_all      = haloMassBinAccRate(sP=sPg,sgSelect=sgSelect,timeWindow=timeWindow,accMode='all')
  arA_all      = haloMassBinAccRate(sP=sPa,sgSelect=sgSelect,timeWindow=timeWindow,accMode='all')
  arG_smooth   = haloMassBinAccRate(sP=sPg,sgSelect=sgSelect,timeWindow=timeWindow,accMode='smooth')
  arA_smooth   = haloMassBinAccRate(sP=sPa,sgSelect=sgSelect,timeWindow=timeWindow,accMode='smooth')
  arG_clumpy   = haloMassBinAccRate(sP=sPg,sgSelect=sgSelect,timeWindow=timeWindow,accMode='clumpy')
  arA_clumpy   = haloMassBinAccRate(sP=sPa,sgSelect=sgSelect,timeWindow=timeWindow,accMode='clumpy') 
  arG_stripped = haloMassBinAccRate(sP=sPg,sgSelect=sgSelect,timeWindow=timeWindow,accMode='stripped')
  arA_stripped = haloMassBinAccRate(sP=sPa,sgSelect=sgSelect,timeWindow=timeWindow,accMode='stripped')
  stop
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
    
    cgPlot,cfG_all.logMassBinCen,smooth(cfG_all.medianVals.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,cfA_all.logMassBinCen,smooth(cfA_all.medianVals.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,cfG_smooth.logMassBinCen,smooth(cfG_smooth.medianVals.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,cfA_smooth.logMassBinCen,smooth(cfA_smooth.medianVals.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,cfG_clumpy.logMassBinCen,smooth(cfG_clumpy.medianVals.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,cfA_clumpy.logMassBinCen,smooth(cfA_clumpy.medianVals.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,cfG_stripped.logMassBinCen,smooth(cfG_stripped.medianVals.both_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,cfA_stripped.logMassBinCen,smooth(cfA_stripped.medianVals.both_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    ; top left axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    cgAxis,/xaxis,xrange=tvirrange,/xs
    
    ; ur: gal const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    cgPlot,cfG_all.logMassBinCen,smooth(cfG_all.medianVals.both_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,cfA_all.logMassBinCen,smooth(cfA_all.medianVals.both_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,cfG_smooth.logMassBinCen,smooth(cfG_smooth.medianVals.both_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,cfA_smooth.logMassBinCen,smooth(cfA_smooth.medianVals.both_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,cfG_clumpy.logMassBinCen,smooth(cfG_clumpy.medianVals.both_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,cfA_clumpy.logMassBinCen,smooth(cfA_clumpy.medianVals.both_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,cfG_stripped.logMassBinCen,smooth(cfG_stripped.medianVals.both_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,cfA_stripped.logMassBinCen,smooth(cfA_stripped.medianVals.both_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    ; top right axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    tvirrange[1] = round(10*tvirrange[1])/10.0
    cgAxis,/xaxis,xrange=tvirrange,/xs

    ; ll: gmem tvir (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    cgPlot,cfG_all.logMassBinCen,smooth(cfG_all.medianVals.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,cfA_all.logMassBinCen,smooth(cfA_all.medianVals.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,cfG_smooth.logMassBinCen,smooth(cfG_smooth.medianVals.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,cfA_smooth.logMassBinCen,smooth(cfA_smooth.medianVals.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,cfG_clumpy.logMassBinCen,smooth(cfG_clumpy.medianVals.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,cfA_clumpy.logMassBinCen,smooth(cfA_clumpy.medianVals.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,cfG_stripped.logMassBinCen,smooth(cfG_stripped.medianVals.gmem_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,cfA_stripped.logMassBinCen,smooth(cfA_stripped.medianVals.gmem_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    ; legend
    legend,['all','smooth','clumpy','stripped'],linestyle=[0,1,2,3],box=0,linesize=0.4,$
      position=[11.2,0.95],charsize=!p.charsize-0.2

    ; lr: gmem const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    cgPlot,cfG_all.logMassBinCen,smooth(cfG_all.medianVals.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,cfA_all.logMassBinCen,smooth(cfA_all.medianVals.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,cfG_smooth.logMassBinCen,smooth(cfG_smooth.medianVals.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,cfA_smooth.logMassBinCen,smooth(cfA_smooth.medianVals.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,cfG_clumpy.logMassBinCen,smooth(cfG_clumpy.medianVals.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,cfA_clumpy.logMassBinCen,smooth(cfA_clumpy.medianVals.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,cfG_stripped.logMassBinCen,smooth(cfG_stripped.medianVals.gmem_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,cfA_stripped.logMassBinCen,smooth(cfA_stripped.medianVals.gmem_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

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
    
    cgPlot,cfG_all.logMassBinCen,smooth(cfG_all.medianVals.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,cfA_all.logMassBinCen,smooth(cfA_all.medianVals.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,cfG_smooth.logMassBinCen,smooth(cfG_smooth.medianVals.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,cfA_smooth.logMassBinCen,smooth(cfA_smooth.medianVals.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,cfG_clumpy.logMassBinCen,smooth(cfG_clumpy.medianVals.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,cfA_clumpy.logMassBinCen,smooth(cfA_clumpy.medianVals.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,cfG_stripped.logMassBinCen,smooth(cfG_stripped.medianVals.gal_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,cfA_stripped.logMassBinCen,smooth(cfA_stripped.medianVals.gal_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    ; top left axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    cgAxis,/xaxis,xrange=tvirrange,/xs
    
    ; ur: gal gas only const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    cgPlot,cfG_all.logMassBinCen,smooth(cfG_all.medianVals.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,cfA_all.logMassBinCen,smooth(cfA_all.medianVals.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,cfG_smooth.logMassBinCen,smooth(cfG_smooth.medianVals.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,cfA_smooth.logMassBinCen,smooth(cfA_smooth.medianVals.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,cfG_clumpy.logMassBinCen,smooth(cfG_clumpy.medianVals.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,cfA_clumpy.logMassBinCen,smooth(cfA_clumpy.medianVals.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,cfG_stripped.logMassBinCen,smooth(cfG_stripped.medianVals.gal_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,cfA_stripped.logMassBinCen,smooth(cfA_stripped.medianVals.gal_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    ; top right axis
    tvirrange = alog10(codeMassToVirTemp((10.0^xrange)/1e10,redshift=sPg.redshift))
    tvirrange[1] = round(10*tvirrange[1])/10.0
    cgAxis,/xaxis,xrange=tvirrange,/xs

    ; ll: stars only tvir (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    cgPlot,cfG_all.logMassBinCen,smooth(cfG_all.medianVals.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,cfA_all.logMassBinCen,smooth(cfA_all.medianVals.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,cfG_smooth.logMassBinCen,smooth(cfG_smooth.medianVals.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,cfA_smooth.logMassBinCen,smooth(cfA_smooth.medianVals.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,cfG_clumpy.logMassBinCen,smooth(cfG_clumpy.medianVals.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,cfA_clumpy.logMassBinCen,smooth(cfA_clumpy.medianVals.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,cfG_stripped.logMassBinCen,smooth(cfG_stripped.medianVals.stars_tViracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,cfA_stripped.logMassBinCen,smooth(cfA_stripped.medianVals.stars_tViracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

    ; legend
    legend,['all','smooth','clumpy','stripped'],linestyle=[0,1,2,3],box=0,linesize=0.4,$
      position=[11.2,0.95],charsize=!p.charsize-0.2

    ; lr: stars only const (gadget/arepo)
    ; ---
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    cgPlot,cfG_all.logMassBinCen,smooth(cfG_all.medianVals.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,cfA_all.logMassBinCen,smooth(cfA_all.medianVals.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,cfG_smooth.logMassBinCen,smooth(cfG_smooth.medianVals.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,cfA_smooth.logMassBinCen,smooth(cfA_smooth.medianVals.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,cfG_clumpy.logMassBinCen,smooth(cfG_clumpy.medianVals.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,cfA_clumpy.logMassBinCen,smooth(cfA_clumpy.medianVals.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,cfG_stripped.logMassBinCen,smooth(cfG_stripped.medianVals.stars_const[conInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,cfA_stripped.logMassBinCen,smooth(cfA_stripped.medianVals.stars_const[conInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot

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
    
    cgPlot,arG_all.logMassBinCen,smooth(arG_all.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,arA_all.logMassBinCen,smooth(arA_all.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,arG_smooth.logMassBinCen,smooth(arG_smooth.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,arA_smooth.logMassBinCen,smooth(arA_smooth.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,arG_clumpy.logMassBinCen,smooth(arG_clumpy.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,arA_clumpy.logMassBinCen,smooth(arA_clumpy.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,arG_stripped.logMassBinCen,smooth(arG_stripped.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,arA_stripped.logMassBinCen,smooth(arA_stripped.hotMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; ur: cold gal gas (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    cgPlot,arG_all.logMassBinCen,smooth(arG_all.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,arA_all.logMassBinCen,smooth(arA_all.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,arG_smooth.logMassBinCen,smooth(arG_smooth.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,arA_smooth.logMassBinCen,smooth(arA_smooth.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,arG_clumpy.logMassBinCen,smooth(arG_clumpy.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,arA_clumpy.logMassBinCen,smooth(arA_clumpy.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,arG_stripped.logMassBinCen,smooth(arG_stripped.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,arA_stripped.logMassBinCen,smooth(arA_stripped.coldMedian.gal_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; ll: hot stars (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    cgPlot,arG_all.logMassBinCen,smooth(arG_all.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,arA_all.logMassBinCen,smooth(arA_all.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,arG_smooth.logMassBinCen,smooth(arG_smooth.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,arA_smooth.logMassBinCen,smooth(arA_smooth.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,arG_clumpy.logMassBinCen,smooth(arG_clumpy.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,arA_clumpy.logMassBinCen,smooth(arA_clumpy.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,arG_stripped.logMassBinCen,smooth(arG_stripped.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,arA_stripped.logMassBinCen,smooth(arA_stripped.hotMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; lr: cold stars (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    cgPlot,arG_all.logMassBinCen,smooth(arG_all.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,arA_all.logMassBinCen,smooth(arA_all.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,arG_smooth.logMassBinCen,smooth(arG_smooth.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,arA_smooth.logMassBinCen,smooth(arA_smooth.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,arG_clumpy.logMassBinCen,smooth(arG_clumpy.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,arA_clumpy.logMassBinCen,smooth(arA_clumpy.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,arG_stripped.logMassBinCen,smooth(arG_stripped.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,arA_stripped.logMassBinCen,smooth(arA_stripped.coldMedian.stars_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
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
   
    ; ul: hot gal (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    cgPlot,arG_all.logMassBinCen,smooth(arG_all.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,arA_all.logMassBinCen,smooth(arA_all.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,arG_smooth.logMassBinCen,smooth(arG_smooth.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,arA_smooth.logMassBinCen,smooth(arA_smooth.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,arG_clumpy.logMassBinCen,smooth(arG_clumpy.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,arA_clumpy.logMassBinCen,smooth(arA_clumpy.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,arG_stripped.logMassBinCen,smooth(arG_stripped.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,arA_stripped.logMassBinCen,smooth(arA_stripped.hotMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/top,/left
    
    ; ur: cold gal (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),pos=pos[1],/noerase
    
    cgPlot,arG_all.logMassBinCen,smooth(arG_all.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,arA_all.logMassBinCen,smooth(arA_all.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,arG_smooth.logMassBinCen,smooth(arG_smooth.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,arA_smooth.logMassBinCen,smooth(arA_smooth.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,arG_clumpy.logMassBinCen,smooth(arG_clumpy.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,arA_clumpy.logMassBinCen,smooth(arA_clumpy.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,arG_stripped.logMassBinCen,smooth(arG_stripped.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,arA_stripped.logMassBinCen,smooth(arA_stripped.coldMedian.both_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; ll: hot gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",pos=pos[2],/noerase
    
    cgPlot,arG_all.logMassBinCen,smooth(arG_all.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,arA_all.logMassBinCen,smooth(arA_all.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,arG_smooth.logMassBinCen,smooth(arG_smooth.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,arA_smooth.logMassBinCen,smooth(arA_smooth.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,arG_clumpy.logMassBinCen,smooth(arG_clumpy.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,arA_clumpy.logMassBinCen,smooth(arA_clumpy.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,arG_stripped.logMassBinCen,smooth(arG_stripped.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,arA_stripped.logMassBinCen,smooth(arA_stripped.hotMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
    ; lr: cold gmem (gadget/arepo)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle="",xtitle="",ytickname=replicate(' ',10),pos=pos[3],xticks=3,xtickv=[10.5,11.0,11.5,12.0],/noerase
    
    cgPlot,arG_all.logMassBinCen,smooth(arG_all.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=0,/overplot
    cgPlot,arA_all.logMassBinCen,smooth(arA_all.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=0,/overplot
    cgPlot,arG_smooth.logMassBinCen,smooth(arG_smooth.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=1,/overplot
    cgPlot,arA_smooth.logMassBinCen,smooth(arA_smooth.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=1,/overplot

    cgPlot,arG_clumpy.logMassBinCen,smooth(arG_clumpy.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=2,/overplot
    cgPlot,arA_clumpy.logMassBinCen,smooth(arA_clumpy.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=2,/overplot
    cgPlot,arG_stripped.logMassBinCen,smooth(arG_stripped.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPg.colorsG[cInd],line=3,/overplot
    cgPlot,arA_stripped.logMassBinCen,smooth(arA_stripped.coldMedian.gmem_tviracc[virInd,*],sK,/nan),color=sPa.colorsA[cInd],line=3,/overplot
    
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

end

; plotCurAngMomRad

pro plotCurAngMomRad

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  res = 128
  sP = simParams(res=res,run='tracer',redshift=2.0)

  ; load group catalog and make halo selection
  hMassRange = [11.75,12.0]
  
  gc = loadGroupCat(sP=sP,/readIDs)
  gcIDs = gcIDList(gc=gc,select='pri')
  gcMasses = codeMassToLogMsun(gc.subgroupMass[gcIDs])
  
  w = where(gcMasses ge hMassRange[0] and gcMasses lt hMassRange[1],count)
  
  gcIDsKeep = gcIDs[w]
  
  ; load subgroup centers
  sgcen = subgroupPosByMostBoundID(sP=sP)
  
  ; load gas IDs and match
  gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
  gas_pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  gas_vel = loadSnapshotSubset(sP=sP,partType='gas',field='vel')
  
  ; keeper arrays
  if sP.trMCPerCell eq 0 then begin
    ; one per gas particle
    jnorms    = fltarr(total(gc.subgroupLenType[partTypeNum('gas'),gcIDsKeep],/int))
    radii     = fltarr(total(gc.subgroupLenType[partTypeNum('gas'),gcIDsKeep],/int))
    parmasses = fltarr(total(gc.subgroupLenType[partTypeNum('gas'),gcIDsKeep],/int))
  endif else begin
    ; find total number of MC tracers in these gas cells, one per tracerMC
    gcPIDsGas = gcPIDList(gc=gc,valGCids=gcIDsKeep,partType='gas',select='all')
    tr_inds = cosmoTracerChildren(sP=sP, /getInds, gasIDs=gcPIDsGas, child_counts=gas_cc)
    
    jnorms    = fltarr(n_elements(tr_inds))
    radii     = fltarr(n_elements(tr_inds))
    parmasses = fltarr(n_elements(tr_inds))
  endelse
  count = 0L
  
  ; loop over each halo
  foreach gcID,gcIDsKeep,k do begin
    print,k,n_elements(gcIDsKeep)
  
    ; get local ID list and match
    gcPIDsGas = gcPIDList(gc=gc,valGCids=[gcID],partType='gas',select='all') ; effectively pri only
    match,gcPIDsGas,gas_ids,gc_ind,ids_ind,count=countMatch,/sort
    ids_ind = ids_ind[sort(gc_ind)] ; reorder ids_ind to be in the order of gcPIDsGas
    if n_elements(gcPIDsGas) ne countMatch then message,'Error'

    ; if doing tracers, modify countMatch to number of child tracers in gcPIDsGas
    ; and set ids_ind to replicated indices according to child counts
    if sP.trMCPerCell ne 0 then begin
      tr_inds = cosmoTracerChildren(sP=sP, /getInds, gasIDs=gcPIDsGas, child_counts=gas_cc)
      gas_ids_ind = ids_ind ; preserve
      ids_ind = lonarr(total(gas_cc,/int))
      cc = 0L
      for i=0,n_elements(gas_cc)-1 do begin
        if gas_cc[i] gt 0 then ids_ind[cc:cc+gas_cc[i]-1] = gas_ids_ind[i]
        cc += gas_cc[i]
      endfor
      countMatch = total(gas_cc,/int)
    endif
    
    ; make positions and velocities relative to halo
    loc_pos = fltarr(3,countMatch)
    loc_vel = fltarr(3,countMatch)
    
    loc_pos[0,*] = sgcen[0,gcID] - gas_pos[0,ids_ind]
    loc_pos[1,*] = sgcen[1,gcID] - gas_pos[1,ids_ind]
    loc_pos[2,*] = sgcen[2,gcID] - gas_pos[2,ids_ind]
    
    loc_vel[0,*] = gas_vel[0,ids_ind] - gc.subgroupVel[0,gcID]
    loc_vel[1,*] = gas_vel[1,ids_ind] - gc.subgroupVel[1,gcID]
    loc_vel[2,*] = gas_vel[2,ids_ind] - gc.subgroupVel[2,gcID]
    
    ; calculate angular momenta
    jvec = fltarr(3,countMatch)
    jvec[0,*] = loc_pos[1,*] * loc_vel[2,*] - loc_pos[2,*] * loc_vel[1,*]
    jvec[1,*] = loc_pos[2,*] * loc_vel[0,*] - loc_pos[0,*] * loc_vel[2,*]
    jvec[2,*] = loc_pos[0,*] * loc_vel[1,*] - loc_pos[1,*] * loc_vel[0,*]
    
    halo_rvir = (gc.group_r_crit200[gc.subgroupGrNr[gcID]])[0]
    
    radii_loc = reform(sqrt(loc_pos[0,*]*loc_pos[0,*] + $
                            loc_pos[1,*]*loc_pos[1,*] + $
                            loc_pos[2,*]*loc_pos[2,*])) / halo_rvir
                            
    jnorms_loc = reform(sqrt(jvec[0,*]*jvec[0,*] + jvec[1,*]*jvec[1,*] + jvec[2,*]*jvec[2,*]))
    
    ; store
    jnorms[count:count+countMatch-1] = jnorms_loc
    radii[count:count+countMatch-1] = radii_loc
    parmasses[count:count+countMatch-1] = codeMassToLogMsun(gc.subgroupMass[gcID])
    
    count += countMatch
    
  endforeach
  
  ; plot
  start_PS, sP.plotPath + 'angmom.test.'+sP.plotPrefix+'.eps'
    xrange = [0,1.5]
    yrange = [0.001,1.0]
    binsize = 0.05
      
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      xtitle="radius [kpc]",ytitle="frac"  
    
    hist = histogram(radii,binsize=binsize,loc=loc)
    cgPlot,loc+binsize*0.5,hist/float(total(hist)),line=0,/overplot
  end_PS
  
  start_PS, sP.plotPath + 'angmom.test2.'+sP.plotPrefix+'.eps'
    xrange = [1,6]
    yrange = [0.001,1.0]
    binsize = 0.1
      
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      xtitle="log j [kpc km/s]",ytitle="frac"  
    
    hist = histogram(alog10(jnorms),binsize=binsize,loc=loc)
    cgPlot,loc+binsize*0.5,hist/float(total(hist)),line=0,/overplot
  end_PS
  
  start_PS, sP.plotPath + 'angmom.test3.'+sP.plotPrefix+'.eps'
    xrange = [0,1.5]
    yrange = [1,6]
    binsize = 0.05
      
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      xtitle="radius [kpc]",ytitle="log j [kpc km/s]"  
    
    nbins = fix((xrange[1]-xrange[0])/binsize)
    
    vals = fltarr(nbins)
    bincens = fltarr(nbins)
    
    for i=0,nbins-1 do begin
      binstart = 0 + (i+0)*binsize
      binend   = 0 + (i+1)*binsize
      bincens[i] = mean([binstart,binend])
      
      w = where(radii ge binstart and radii lt binend,count)
      if count gt 0 then vals[i] = median(alog10(jnorms[w]))
    endfor
    
    cgPlot,bincens,vals,psym=-4,/overplot
  end_PS
end

; plotAngMomRadVsHaloMass(): plot the median angular momentum of hot/cold modes at different radial
;                            times as a function of halo mass

pro plotAngMomRadVsHaloMass

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  sgSelect = 'pri'
  accMode = 'all'
  
  res = 256
  sPg = simParams(res=res,run='gadget',redshift=2.0)
  sPa = simParams(res=res,run='tracer',redshift=2.0)
  
  binnedGadget = haloMassBinAngMom(sP=sPg,sgSelect=sgSelect,accMode=accMode)
  binnedArepo  = haloMassBinAngMom(sP=sPa,sgSelect=sgSelect,accMode=accMode)
  
  ; load gc and calculate jcirc normalizations
  gc = loadGroupCat(sP=sPg,/skipIDs)
  gcMasses = codeMassToLogMsun(gc.subgroupMass)
  
  rvirs  = fltarr(n_elements(binnedGadget.logMassBinCen))
  vcircs = fltarr(n_elements(binnedGadget.logMassBinCen))
  jcircs = fltarr(n_elements(binnedGadget.logMassBinCen))
  
  mWidth = binnedGadget.logMassBinCen[1] - binnedGadget.logMassBinCen[0]
  for i=0,n_elements(binnedGadget.logMassBinCen)-1 do begin
    ; locate halos in this mass bin and get a mean rvir
    w = where(gcMasses ge binnedGadget.logMassBinCen[i]-mWidth*0.5 and $
              gcMasses ge binnedGadget.logMassBinCen[i]+mWidth*0.5,count)
    if count gt 0 then begin
      rvirs[i] = mean(gc.group_r_crit200[gc.subgroupGrNr[w]])
      vcircs[i] = sqrt(units.G * (10.0^binnedGadget.logMassBinCen[i])/1e10 / rvirs[i]) ;v_circ = sqrt(GM/R)
      jcircs[i] = rvirs[i] * vcircs[i]
    endif
  endfor
  
  if jcircs[-1] eq 0 then jcircs[-1] = jcircs[-2] ; most massive bin empty (too small width)
  
  ; minimum halo mass to plot
  w = where(binnedGadget.logMassBinCen gt 8.0)
  sK = 3 ; smoothing kernel size

  ; plot (1) - lines for all rVirFacs vs halo mass
  start_PS, sPg.plotPath + 'angmom.vshalo.comp.'+accMode+'.'+str(sPg.res)+'_'+str(sPg.snap)+'.eps', /big
    
    xrange = [9.75,12.5]
    yrange = [5e1,1e5]
    
    xtickv = [10,11,12]
    
    x0 = 0.15 & x1 = 0.42 & x2 = 0.69 & x3 = 0.96
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; uc
                [x2,y1,x3,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1] ,$ ; lc
                [x2,y0,x3,y1] )  ; lr
    
    for j=1,n_elements(binnedGadget.rVirFacs)-1 do begin
      
      if j eq 1 or j eq 4 then ytickname = '' else ytickname = replicate(' ',10)
      if j ge 4 then xtickname = '' else xtickname = replicate(' ',10)
      if j gt 1 then noerase = 1 else noerase = 0
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,pos=pos[j-1],$
        ytitle="",xtitle="",xticks=n_elements(xtickv)-1,xtickv=xtickv,$
        xtickname=xtickname,ytickname=ytickname,noerase=noerase       
     
      ; gadget hot both
      cgPlot,binnedGadget.logMassBinCen[w],smooth(binnedGadget.hotMode.median_both[j,w],sK,/nan),color=getColor(1),line=1,/overplot
      cgPlot,binnedGadget.logMassBinCen[w],smooth(binnedGadget.coldMode.median_both[j,w],sK,/nan),color=getColor(1),line=2,/overplot

      ; arepo hot both
      cgPlot,binnedArepo.logMassBinCen[w],smooth(binnedArepo.hotMode.median_both[j,w],sK,/nan),color=getColor(3),line=1,/overplot
      cgPlot,binnedArepo.logMassBinCen[w],smooth(binnedArepo.coldMode.median_both[j,w],sK,/nan),color=getColor(3),line=2,/overplot
      
      ; legend
      massBinStr = textoidl('r/r_{vir}')+' = ' + string(binnedGadget.rVirFacs[j],format='(f4.2)')
      cgText,mean(xrange),yrange[1]*0.45,massBinStr,charsize=!p.charsize-0.2,alignment=0.5,$
        color=cgColor('forest green')
        
      if j eq 1 then legend,['cold','hot'],linestyle=[1,2],box=0,linesize=0.25,/bottom,/left
    
    endfor
    
    ; axis labels
    cgText,0.05,0.5,textoidl("< j_{gas} > [kpc km/s]"),alignment=0.5,orientation=90.0,/normal
    cgText,0.5,0.05,textoidl("M_{halo} [_{ }log M_{sun }]"),alignment=0.5,/normal
    
    ; annotations
    legend,['gadget','arepo'],textcolor=getColor([1,3],/name),box=0,position=[9.8,3e4]
    
  end_PS
  
  ; plot (2) - ratio of arepo/gadget for all radii
  start_PS, sPg.plotPath + 'angmom.vshalo.comp2.'+accMode+'.'+str(sPg.res)+'_'+str(sPg.snap)+'.eps'
    
    xrange = [9.75,12.5]
    yrange = [0.5,3.0]
    xtickv = [10.0,11.0,12.0]
    
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
        ytitle=textoidl("< j_{arepo} > / < j_{gadget} >"),$
        xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),$
        xticks=n_elements(xtickv)-1,xtickv=xtickv 
    
      cgPlot,xrange,[1.0,1.0],line=0,color=cgColor('light gray'),/overplot
    
    for j=0,n_elements(binnedGadget.rVirFacs)-1 do begin
      
      w = where(binnedGadget.logMassBinCen gt 8.0)
      
      ; arepo/gadget hot both
      vals = binnedArepo.hotMode.median_both[j,w] / binnedGadget.hotMode.median_both[j,w]
      cgPlot,binnedGadget.logMassBinCen[w],smooth(vals,sK,/nan),color=getColor(j),line=1,/overplot
      
      ; arepo/gadget cold both
      vals = binnedArepo.coldMode.median_both[j,w] / binnedGadget.coldMode.median_both[j,w]
      cgPlot,binnedGadget.logMassBinCen[w],smooth(vals,sK,/nan),color=getColor(j),line=2,/overplot
    
    endfor
    
    ; legend
    strings = string(binnedGadget.rVirFacs,format='(f4.2)')
    colors = getColor(indgen(n_elements(binnedGadget.rVirFacs)),/name)
    legend,strings,textcolor=colors,box=0,/top,/right
    legend,['cold','hot'],linestyle=[1,2],box=0,linesize=0.25,/top,/left
    
  end_PS
  
  ; plot (3) - arepo and gadget for each radii on a single plot (since they are offset)
  start_PS, sPg.plotPath + 'angmom.vshalo.comp3.'+accMode+'.'+str(sPg.res)+'_'+str(sPg.snap)+'.eps',$
    xs=7.5,ys=9.0
    
    xrange = [9.75,12.5]
    yrange = [0.01,2]
    xtickv = [10.0,10.5,11.0,11.5,12.0,12.5]
    
    yposl = [1.4,5e4,0.4,2e4,0.13,8e3,0.025]
    
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
        ytitle=textoidl("< j_{gas} > / j_{circ}(r_{vir})"),$
        xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),$
        xticks=n_elements(xtickv)-1,xtickv=xtickv
    
    sK=1
    
    for j=0,n_elements(binnedGadget.rVirFacs)-1,2 do begin
      
      w = where(binnedGadget.logMassBinCen gt 8.0)
      
      ; gadget hot both
      cgPlot,binnedGadget.logMassBinCen[w],smooth(binnedGadget.hotMode.median_both[j,w]/jcircs[w],sK,/nan),color=getColor(1),line=1,/overplot
      cgPlot,binnedGadget.logMassBinCen[w],smooth(binnedGadget.coldMode.median_both[j,w]/jcircs[w],sK,/nan),color=getColor(1),line=2,/overplot

      ; arepo hot both
      cgPlot,binnedArepo.logMassBinCen[w],smooth(binnedArepo.hotMode.median_both[j,w]/jcircs[w],sK,/nan),color=getColor(3),line=1,/overplot
      cgPlot,binnedArepo.logMassBinCen[w],smooth(binnedArepo.coldMode.median_both[j,w]/jcircs[w],sK,/nan),color=getColor(3),line=2,/overplot
      
      ; legend
      massBinStr = textoidl("r/r_{vir} = ")+string(binnedGadget.rVirFacs[j],format='(f4.2)')
      cgText,12.2,yposl[j],massBinStr,charsize=!p.charsize-0.2,alignment=0.5,$
        color=cgColor('forest green')
    
    endfor
    
    ; annotations
    legend,['gadget','arepo'],textcolor=getColor([1,3],/name),box=0,/bottom,/left
    legend,['cold','hot'],linestyle=[1,2],box=0,linesize=0.25,/bottom,/right
    
  end_PS
  stop
end

; plotDeltaAccTimeVsHaloMass(): plot the mean time for accreting gas to reach various radii from the 
;                               virial radius (normalized by t_circ)

pro plotDeltaAccTimeVsHaloMass

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect = 'pri'
  accMode = 'all'
  
  res = 256
  sPg = simParams(res=res,run='gadget',redshift=2.0)
  sPa = simParams(res=res,run='tracer',redshift=2.0)
  
  binnedGadget = haloMassBinDeltaAccTime(sP=sPg,sgSelect=sgSelect,accMode=accMode)
  binnedArepo  = haloMassBinDeltaAccTime(sP=sPa,sgSelect=sgSelect,accMode=accMode)

  ; minimum halo mass to plot
  w = where(binnedGadget.logMassBinCen gt 8.0)
  sK = 3 ; smoothing kernel size

  ; plot (1) - lines for all rVirFacs vs halo mass
  start_PS, sPg.plotPath + 'accdt.vshalo.comp.'+accMode+'.'+str(sPg.res)+'_'+str(sPg.snap)+'.eps', /big
    
    xrange = [9.5,12.5]
    yrange = [0.0,0.22]
    
    xtickv = [10,11,12]
    
    x0 = 0.15 & x1 = 0.42 & x2 = 0.69 & x3 = 0.96
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; uc
                [x2,y1,x3,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1] ,$ ; lc
                [x2,y0,x3,y1] )  ; lr
    
    for j=1,n_elements(binnedGadget.rVirFacs)-1 do begin
      
      if j eq 1 or j eq 4 then ytickname = '' else ytickname = replicate(' ',10)
      if j ge 4 then xtickname = '' else xtickname = replicate(' ',10)
      if j gt 1 then noerase = 1 else noerase = 0
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,pos=pos[j-1],$
        ytitle="",xtitle="",xticks=n_elements(xtickv)-1,xtickv=xtickv,$
        xtickname=xtickname,ytickname=ytickname,noerase=noerase       
      
      cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
      
      ; gadget hot both
      cgPlot,binnedGadget.logMassBinCen[w],smooth(binnedGadget.hotMode.median_both[j,w],sK,/nan),color=getColor(1),line=1,/overplot
      cgPlot,binnedGadget.logMassBinCen[w],smooth(binnedGadget.coldMode.median_both[j,w],sK,/nan),color=getColor(1),line=2,/overplot

      ; arepo hot both
      cgPlot,binnedArepo.logMassBinCen[w],smooth(binnedArepo.hotMode.median_both[j,w],sK,/nan),color=getColor(3),line=1,/overplot
      cgPlot,binnedArepo.logMassBinCen[w],smooth(binnedArepo.coldMode.median_both[j,w],sK,/nan),color=getColor(3),line=2,/overplot
      
      ; legend
      massBinStr = textoidl('r/r_{vir}')+' = ' + string(binnedGadget.rVirFacs[j],format='(f4.2)')
      cgText,mean(xrange),yrange[1]*0.89,massBinStr,charsize=!p.charsize-0.2,alignment=0.5,$
        color=cgColor('forest green')
      if j eq 1 then legend,['gadget','arepo'],textcolor=getColor([1,3],/name),box=0,position=[9.6,0.13]
      
    endfor
    
    ; axis labels
    cgText,0.05,0.5,textoidl("\Delta t_{acc} / \tau_{circ}"),alignment=0.5,orientation=90.0,/normal
    cgText,0.5,0.05,textoidl("log ( M_{halo} ) [_{ }M_{sun }]"),alignment=0.5,/normal
    
    ; annotations
    
    legend,['cold','hot'],linestyle=[1,2],box=0,linesize=0.25,/bottom,/left
    
  end_PS
  
  ; plot (2) - ratio of arepo/gadget for all radii
  start_PS, sPg.plotPath + 'accdt.vshalo.comp2.'+accMode+'.'+str(sPg.res)+'_'+str(sPg.snap)+'.eps'
    
    xrange = [10.0,12.5]
    yrange = [0.5,2.5]
    xtickv = [10.0,10.5,11.0,11.5,12.0,12.5]
    
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
        ytitle=textoidl("< \Delta t_{acc,arepo} > / < \Delta t_{acc,gadget} >"),$
        xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),$
        xticks=n_elements(xtickv)-1,xtickv=xtickv 
    
      cgPlot,xrange,[1.0,1.0],line=0,color=cgColor('light gray'),/overplot
    
    for j=1,n_elements(binnedGadget.rVirFacs)-1 do begin
      
      w = where(binnedGadget.logMassBinCen gt 8.0)
      
      ; arepo/gadget hot both
      vals = binnedArepo.hotMode.median_both[j,w] / binnedGadget.hotMode.median_both[j,w]
      cgPlot,binnedGadget.logMassBinCen[w],smooth(vals,sK,/nan),color=getColor(j),line=1,/overplot
      
      ; arepo/gadget cold both
      vals = binnedArepo.coldMode.median_both[j,w] / binnedGadget.coldMode.median_both[j,w]
      cgPlot,binnedGadget.logMassBinCen[w],smooth(vals,sK,/nan),color=getColor(j),line=2,/overplot
    
    endfor
    
    ; legend
    strings = string(binnedGadget.rVirFacs[1:*],format='(f4.2)')
    colors = getColor(indgen(n_elements(binnedGadget.rVirFacs)-1)+1,/name)
    legend,strings,textcolor=colors,box=0,/top,/right
    legend,['cold','hot'],linestyle=[1,2],box=0,linesize=0.25,/top,/left
    
  end_PS
  stop
end

; plotModeMassesVsHaloMass(): plot the "cold mass" and "hot mass" (not fraction) vs halo mass

pro plotModeMassesVsHaloMass

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect = 'pri'
  accMode  = 'smooth' ; accretion mode: all, smooth, bclumpy, sclumpy
  
  res = 256
  sPg = simParams(res=res,run='gadget',redshift=2.0)
  sPa = simParams(res=res,run='tracer',redshift=2.0)
  
  mmG = haloMassBinModeMasses(sP=sPg,sgSelect=sgSelect,accMode=accMode)
  mmA = haloMassBinModeMasses(sP=sPa,sgSelect=sgSelect,accMode=accMode)
  
  ; plot - cold,hot,total masses
  xrange = [10.5,mmG.xrange[1]]
  w = where(mmG.logMassBinCen gt 10.5)
    
  ; plot (1) - tviracc only, cold+hot+total
  start_PS, sPg.plotPath + 'massBudget.'+accMode+'.comp.'+str(sPg.res)+'_'+str(sPg.snap)+'.eps',$
    xs = 7.5, ys = 10
    
    x0 = 0.15 & x1 = 0.96
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; uc
                [x0,y0,x1,y1]  ) ; lc
   
    ; gal
    yrange = [min(mmG.hotMass.tviracc_gal,/nan)*0.97,max(mmA.totalMass.tviracc_gal,/nan)*1.03]
    yrange = [8.0,12.0]
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle="",xtickname=replicate(' ',10),pos=pos[0]
    
    ;xpts = findgen(50)/50.0*2.0 + 10.5
    ;cgPlot,xpts,mmG.tt_gal,line=3,/overplot
    
    cgPlot,mmG.logMassBinCen[w],mmG.coldMass.tViracc_gal[w],color=getColor(1),line=1,/overplot
    cgPlot,mmG.logMassBinCen[w],mmG.hotMass.tViracc_gal[w],color=getColor(1),line=2,/overplot
    cgPlot,mmG.logMassBinCen[w],mmG.totalMass.tViracc_gal[w],color=getColor(1),line=0,/overplot
    
    cgPlot,mmA.logMassBinCen[w],mmA.coldMass.tViracc_gal[w],color=getColor(3),line=1,/overplot
    cgPlot,mmA.logMassBinCen[w],mmA.hotMass.tViracc_gal[w],color=getColor(3),line=2,/overplot
    cgPlot,mmA.logMassBinCen[w],mmA.totalMass.tViracc_gal[w],color=getColor(3),line=0,/overplot
    
    ; legend
    legend,['cold','hot','total'],linestyle=[1,2,0],box=0,linesize=0.25,/bottom,/left
    
    ; gmem
    yrange = [min(mmA.hotMass.tviracc_gmem,/nan)*0.97,max(mmG.totalMass.tviracc_gmem,/nan)*1.03]
    yrange = [10.0,11.75]
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="",xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),pos=pos[1],/noerase
    
    ;xpts = findgen(50)/50.0*2.0 + 10.5
    ;cgPlot,xpts,mmG.tt_gmem,line=3,/overplot
    
    cgPlot,mmG.logMassBinCen[w],mmG.coldMass.tViracc_gmem[w],color=getColor(1),line=1,/overplot
    cgPlot,mmG.logMassBinCen[w],mmG.hotMass.tViracc_gmem[w],color=getColor(1),line=2,/overplot
    cgPlot,mmG.logMassBinCen[w],mmG.totalMass.tViracc_gmem[w],color=getColor(1),line=0,/overplot
    
    cgPlot,mmA.logMassBinCen[w],mmA.coldMass.tViracc_gmem[w],color=getColor(3),line=1,/overplot
    cgPlot,mmA.logMassBinCen[w],mmA.hotMass.tViracc_gmem[w],color=getColor(3),line=2,/overplot
    cgPlot,mmA.logMassBinCen[w],mmA.totalMass.tViracc_gmem[w],color=getColor(3),line=0,/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=getColor([1,3],/name),box=0,/top,/left
    
    ; labels
    cgText,0.05,0.5,"Total Accreted Gas Mass "+textoidl("[_{ }log M_{sun }]"),$
      alignment=0.5,orientation=90.0,/normal
    cgText,mean([x0,x1]),0.52,"Halo Atmosphere",alignment=0.5,color=cgColor('forest green'),/normal
    cgText,mean([x0,x1]),0.92,"Central Galaxy",alignment=0.5,color=cgColor('forest green'),/normal
    
  end_PS
  
end

; plotColdFracVsHaloMassAll(): plot the "cold mode fraction" vs halo mass in a few different ways
;                              without making the mtS/atS cuts

pro plotColdFracVsHaloMassAll, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect = 'pri'
  TcutVals = [5.3,5.4,5.5,5.6,5.7] ; for constant threshold
  nCuts = n_elements(TcutVals)
  
  minNum = 32
  xrange = [9.0,12.5]
  yrange = [0.0,1.1]
  
  logMassBinSize = 0.1
    
  ; load galaxy catalog
  galcat = galaxyCat(sP=sP)
  gcIDList = gcIDList(sP=sP,select=sgSelect)    
    
  ; for all the child gas particles/tracers in the halo selection, replicate parent IDs
  if sP.trMCPerCell eq 0 then begin
    gcIndOrig = galcatRepParentIDs(galcat=galcat,gcIDList=gcIDList)
  endif else begin
    gcIndOrig = galCatRepParentIDs(galcat=galcat,gcIDList=gcIDList,$
                  child_counts={gal:at.child_counts_gal,gmem:at.child_counts_gmem}) 
  endelse
  
  ; compact parents (ascending ID->index)
  placeMap = getIDIndexMap(gcIDList,minid=minid)
  gcIndOrig.gal = placeMap[gcIndOrig.gal-minid]
  gcIndOrig.gmem = placeMap[gcIndOrig.gmem-minid]
  placeMap = !NULL  
  
  ; reverse histogram parent IDs of all particles/tracers in this selection
  hist_gal  = histogram(gcIndOrig.gal,min=0,loc=loc_gal,rev=rev_gal)
  hist_gmem = histogram(gcIndOrig.gmem,min=0,loc=loc_gmem,rev=rev_gmem)

  gcIndOrig = !NULL
  galcat = !NULL
  
  ; load max temps, current tvir, tvir at accretion
  curTvir = gcSubsetProp(sP=sP,select=sgSelect,/virTemp)
  maxTemp = gcSubsetProp(sP=sP,select=sgSelect,/maxPastTemp)

  ; load group cat for subgroup masses
  gc = loadGroupCat(sP=sP,/skipIDs)
  gcMasses = codeMassToLogMsun(gc.subgroupMass[gcIDList])
  gc = !NULL

  ; structures to store results (Tmax)
  coldFrac = { gal_const    : fltarr(nCuts,n_elements(gcIDList))   ,$
               gmem_const   : fltarr(nCuts,n_elements(gcIDList))   ,$
               gal_tvircur  : fltarr(n_elements(gcIDList))         ,$
               gmem_tvircur : fltarr(n_elements(gcIDList))         ,$
               gal_tviracc  : fltarr(n_elements(gcIDList))         ,$
               gmem_tviracc : fltarr(n_elements(gcIDList))         ,$
               both_const   : fltarr(nCuts,n_elements(gcIDList))   ,$
               both_tvircur : fltarr(n_elements(gcIDList))         ,$
               both_tviracc : fltarr(n_elements(gcIDList))         ,$
               gal_num      : lonarr(n_elements(gcIDList))         ,$
               gmem_num     : lonarr(n_elements(gcIDList))          }
  
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
  endfor
  
  ; plot (1) - all data points and median lines (Tmax)
  start_PS, sP.plotPath + 'coldFrac_noMTs.'+sP.plotPrefix+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("log(M_{halo})")+"",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+") noMTs"
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
  start_PS, sP.plotPath + 'coldFrac_noMTs2.'+sP.plotPrefix+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("log(M_{halo})")+"",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+") noMTs"
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
end
