; timeScalesPlot.pro
; plotting: cooling times of halo gas vs. dynamical/hubble timescales
; dnelson feb.2013

; plotTSFracsVsHaloMass(): plot gas mass fractions for timescale ratios vs halo mass

pro plotTSFracsVsHaloMass

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  sgSelect     = 'pri'
  simHotMasses = 1 ; 0=analytical, 1=derive hot halo masses from sim
  
  sP = simParams(res=512,run='tracer',redshift=2.0)
  
  xrange = [9.0,12.0]
  yrange = [0.0,1.03]
  lines  = [1,0,2]
  cInd   = 1 ; colorindex
  
  ; 2d config
  binSizeMass = 0.07 ; 10-12 0.055
  binSizeFrac = 0.04
  
  hsp = [0.004,0.0035] ; mass (0.003 10-12), frac
  nc  = 110 ; number of colors (of 255) to use for background 2d histo

  ; load timescales
  ts = timescaleFracsVsHaloMass(sP=sP,sgSelect=sgSelect)
  
  ; calculate mass fraction with tcool<tdyn for SIS
  mmf = modelMassFracs(sP=sP, ts=ts, simHotMasses=simHotMasses)
  
  simHotMassesTag = 'mAna.'
  if simHotMasses then simHotMassesTag = 'mSim.'

  ; plot (1) - tcool<>tdyn
  start_PS, sP.plotPath + 'tsfrac_tdyn_vsmass.' + simHotMassesTag + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),ytitle="Gas Mass Fraction "+textoidl("(t_{cool} < t_{dyn})")
      
    ; sim 2d histo
    yy = ts.tsFracs.gmem_tcool_tdyn[1,*]
    f2d = binHisto2D(xx=ts.gcMasses, yy=yy, xmm=xrange, ymm=yrange, xbs=binSizeMass, ybs=binSizeFrac)      
    oplot2DHistoSq, f2d, hsp=hsp, nc=nc, xrange=xrange, yrange=yrange, /colNorm, /green
             
    ; sim lines
    j=1
    cgPlot,ts.gcMasses,ts.tsFracs.gmem_tcool_tdyn[j,*],psym=4,color=sP.colorsA[cInd],/overplot
      
    for j=0,n_elements(ts.tsRatioVals)-1 do $
      cgPlot,ts.logMassBinCen,ts.tsMedian.gmem_tcool_tdyn[j,*],color=cgColor('black'),line=lines[j],/overplot
      
    ; models
    mMasses = codeMassToLogMsun(mmf.modelMasses)
    sK = 1
    for j=0,n_elements(ts.tsRatioVals)-1 do begin
      cgPlot,mMasses,smooth(mmf.SIS2_tDyn[j,*],sK),line=lines[j],color=cgColor('orange'),/overplot
      cgPlot,mMasses,smooth(mmf.NFW_iso_tDyn[j,*],sK),line=lines[j],color=cgColor('blue'),/overplot
      cgPlot,mMasses,smooth(mmf.NFW_poly_tDyn[j,*],sK),line=lines[j],color=cgColor('red'),/overplot
    endfor
      
    ; legend
    strings = textoidl("t_{cool} / t_{dyn} < "+string(ts.tsRatioVals,format='(f3.1)'))
    legend,strings,linestyle=lines,box=0,linesize=0.2,/top,/right,$
      textcolors=['black'],color=cgColor('black')
    
    strings = ['SIS','NFW iso','NFW poly','median']
    legend,strings,textcolors=['orange','blue','red','black'],box=0,/top,/right,pos=[11.96,0.82]
    
    ; redo plot borders
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xtitle="",ytitle="",/noerase
    
  end_PS
  
  ; plot (2) - tcool<>tHubble
  start_PS, sP.plotPath + 'tsfrac_tH_vsmass.' + simHotMassesTag + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),ytitle="Gas Mass Fraction "+textoidl("(t_{cool} < t_{H})")
      
    ; sim 2d histo
    yy = ts.tsFracs.gmem_tcool_tage[1,*]
    f2d = binHisto2D(xx=ts.gcMasses, yy=yy, xmm=xrange, ymm=yrange, xbs=binSizeMass, ybs=binSizeFrac)      
    oplot2DHistoSq, f2d, hsp=hsp, nc=nc, xrange=xrange, yrange=yrange, /colNorm, /green
    
    ; sim lines
    j=1
    cgPlot,ts.gcMasses,ts.tsFracs.gmem_tcool_tage[j,*],psym=4,color=sP.colorsA[cInd],/overplot
      
    for j=0,n_elements(ts.tsRatioVals)-1 do $
      cgPlot,ts.logMassBinCen,ts.tsMedian.gmem_tcool_tage[j,*],color=cgColor('black'),line=lines[j],/overplot
      
    ; models
    mMasses = codeMassToLogMsun(mmf.modelMasses)
    sK = 1
    for j=0,n_elements(ts.tsRatioVals)-1 do begin
      cgPlot,mMasses,smooth(mmf.SIS2_tH[j,*],sK),line=lines[j],color=cgColor('orange'),/overplot
      cgPlot,mMasses,smooth(mmf.NFW_iso_tH[j,*],sK),line=lines[j],color=cgColor('blue'),/overplot
      cgPlot,mMasses,smooth(mmf.NFW_poly_tH[j,*],sK),line=lines[j],color=cgColor('red'),/overplot
    endfor
      
    ; legend
    strings = textoidl("t_{cool} / t_{H} < "+string(ts.tsRatioVals,format='(f3.1)'))
    legend,strings,linestyle=lines,box=0,linesize=0.2,/bottom,/left
      
    strings = ['SIS','NFW iso','NFW poly','median']
    legend,strings,textcolors=['orange','blue','red','black'],box=0,/bottom,/left,pos=[9.06,0.22]
      
    ; redo plot borders
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xtitle="",ytitle="",/noerase
  end_PS
  
  ; plot (3) - do hot halo masses (measured) agree with those used for SIS/NFW models?
  yrange = [0.01,2.0]
  
  start_PS, sP.plotPath + 'tsfrac_hotmass.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/ylog,/xs,/ys,yminor=0,$
      xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),ytitle=textoidl("M_{hot,sim} / (f_b M_{DM})")
      
    cgPlot,xrange,[1.0,1.0],line=2,color=cgColor('light gray'),/overplot
      
    ; load gcMasses for DM only
    gc = loadGroupCat(sP=sP,/skipIDs)
    gcIDList = gcIDList(gc=gc,select=sgSelect)
    gcMasses_DM = codeMassToLogMsun(gc.subgroupMassType[partTypeNum('dm'),gcIDList])
      
    ; 2d histo
    yy1 = ts.gmem_hotmasses/(units.f_b*logMsunToCodeMass(ts.gcMasses))
    yy2 = ts.gmem_hotmasses/(units.f_b*logMsunToCodeMass(gcMasses_DM))
    
    ; set a contour palette
    loadColorTable,'brewerc-blues'
    tvlct, rr, gg, bb, /get
    palette = [[rr],[gg],[bb]]
    
    ; contour
    cgd = calcGridData(xx=ts.gcMasses,yy=yy2,xMinMax=xrange,yMinMax=yrange,$
      nPixels=[50,50],/logY)
    
    cgContour, smooth(cgd.dens_out,[3,3]), cgd.xPts, 10.0^cgd.yPts, $
      /overplot, /fill, palette=palette, levels=[0.25,0.5,1,2,5,10,20],c_colors=(indgen(7)*20+50)
    
    ; median lines
    mm1 = fitRadProfile(radii=ts.gcMasses, vals=yy1, range=xrange, radBins=20)
    mm2 = fitRadProfile(radii=ts.gcMasses, vals=yy2, range=xrange, radBins=20)
    cgPlot,mm2.binCen,mm2.radMedian,line=0,/overplot
    cgPlot,mm1.binCen,mm1.radMedian,line=2,/overplot
    
    ; how well is quadratic fit to hotHaloMasses doing?
    ;cgPlot,codeMassToLogMsun(mmf.modelMasses),mmf.hotHaloMasses/(units.f_b*mmf.modelMasses),/overplot
    
    ; redo plot borders
    cgPlot,[0.1],[0.1],/nodata,xrange=xrange,yrange=yrange,/ylog,yminor=0,/xs,/ys,/noerase
    
  end_PS
  
  stop

end

; timescaleRadStack(): 2D bin cooling time and overplot other timescales for 4 stacked halo mass bins
pro timescaleRadStack

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  sP = simParams(res=512,run='tracer',redshift=2.0)
  gc = loadGroupCat(sP=sP,/skipIDs)
  h  = loadSnapshotHeader(sP=sP)
  
  massRanges = list([10.5,10.7], [10.9,11.1], [11.3,11.5], [11.7,11.9])
  
  ; 2d histo config
  binSizeRad  = 0.036
  binSizeTime = 0.12
  
  mmRad  = [0.0,1.5] ; r/rvir
  mmTime = [-3.0,2.0]  ; log(Gyr)
  
  hsp = [0.0036,0.018] ; r/rvir, Gyr
  nc  = 160 ; number of colors (of 255) to use for background 2d histo
  
  ; plot config - 2x2 layout
  radRange = [0.05,1.6]  ; r/rvir
  tsRange  = [0.004,40.0] ; Gyr
  
  x0 = 0.14 & x1 = 0.50 & x2 = 0.85
  y0 = 0.14 & y1 = 0.54 & y2 = 0.94
    
  pos = list( [x0,y1,x1,y2] ,$ ; upper left
              [x1,y1,x2,y2] ,$ ; upper right
              [x0,y0,x1,y1] ,$ ; lower left
              [x1,y0,x2,y1]  ) ; lower right
  
  ; start plot
  start_PS, sP.plotPath + 'timescales_2x2.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
      
  foreach massRange,massRanges,k do begin
    ; make list of halos
    gcIDList = gcIDList(gc=gc,select='pri',massRange=massRange)
                    
    print,k,n_elements(gcIDList)
    
    ; start plot
    if k eq 0 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=radRange,yrange=tsRange,/ylog,/xs,/ys,yminor=0,pos=pos[k],$
      xtickname=replicate(' ',10),/noerase,xtickv=[0.15,0.5,1.0,1.5],xticks=3
    if k eq 1 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=radRange,yrange=tsRange,/ylog,/xs,/ys,yminor=0,pos=pos[k],$
      xtickname=replicate(' ',10),ytickname=replicate(' ',10),/noerase,xtickv=[0.15,0.5,1.0,1.5],xticks=3
    if k eq 2 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=radRange,yrange=tsRange,/ylog,/xs,/ys,yminor=0,pos=pos[k],$
      /noerase,xtickv=[0.15,0.5,1.0,1.5],xticks=3,xtickname=['0.15','0.5','1.0','1.5']
    if k eq 3 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=radRange,yrange=tsRange,/ylog,/xs,/ys,yminor=0,pos=pos[k],$
      ytickname=replicate(' ',10),/noerase,xtickv=[0.15,0.5,1.0,1.5],xticks=3,xtickname=['0.15','0.5','1.0','1.5']
    
    ; load gas timescales and best model fits
    ts = loadFitTimescales(sP=sP,gcIDList=gcIDList)
    
    ; 2d bin cooling times and plot
    ct2d = binHisto2D(xx=ts.gasRadii/ts.gasRvir, yy=alog10(ts.coolTime), wt=ts.masses, $
                      xmm=mmRad, ymm=mmTime, xbs=binSizeRad, ybs=binSizeTime)
    
    oplot2DHistoSq, ct2d, hsp=hsp, nc=nc, xrange=radRange, yrange=tsRange, /nonZero, /logY, /blue
    
    ; radial fits: cooling and dynamical
    cgPlot,ts.radCt.binCen,ts.radCt.radMedian,line=0,color=cgColor('black'),/overplot
    cgPlot,ts.radDt.binCen,ts.radDt.radMedian,line=2,color=cgColor('black'),/overplot
                   
    ; SIS/NFW models: cooling and dynamical
    cgPlot,ts.x,ts.sis_gas.dynTime,line=2,color=cgColor('orange'),/overplot
    cgPlot,ts.x,ts.sis_gas.coolTime,line=0,color=cgColor('orange'),/overplot
    
    cgPlot,ts.x,ts.nfw_gas.dynTime,line=2,color=cgColor('red'),/overplot
    cgPlot,ts.x,ts.nfw_gas.coolTime,line=0,color=cgColor('red'),/overplot
           
    cgText,radRange[1]*0.68,tsRange[0]*2,"M = "+string(mean(massRange),format='(f4.1)'),/data
           
    ; redo clean plot lines
    loadct,0,/silent
    
    if k eq 0 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=radRange,yrange=tsRange,/ylog,/xs,/ys,yminor=0,pos=pos[k],$
      xtickname=replicate(' ',10),/noerase,xtickv=[0.15,0.5,1.0,1.5],xticks=3
    if k eq 1 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=radRange,yrange=tsRange,/ylog,/xs,/ys,yminor=0,pos=pos[k],$
      xtickname=replicate(' ',10),ytickname=replicate(' ',10),/noerase,xtickv=[0.15,0.5,1.0,1.5],xticks=3
    if k eq 2 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=radRange,yrange=tsRange,/ylog,/xs,/ys,yminor=0,pos=pos[k],$
      /noerase,xtickv=[0.15,0.5,1.0,1.5],xticks=3,xtickname=['0.15','0.5','1.0','1.5']
    if k eq 3 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=radRange,yrange=tsRange,/ylog,/xs,/ys,yminor=0,pos=pos[k],$
      ytickname=replicate(' ',10),/noerase,xtickv=[0.15,0.5,1.0,1.5],xticks=3,xtickname=['0.15','0.5','1.0','1.5']
           
    if k eq 0 then $
    legend,textoidl(['\tau_{dyn}','\tau_{cool}']),linestyle=[2,0],box=0,pos=[radRange[0],tsRange[1]*0.7],linesize=0.2
           
  endforeach
  
  ; colorbar
  cbTickNames = linspace(0.0,max(ct2d.h2),7)
  cbTickNames = ['0',string(alog10(cbTickNames[1:*]),format='(f4.1)')]
  
  loadColorTable,'brewerc-blues'
  cgColorbar, range=[0,max(ct2d.h2)], /vertical, ncolors=nc, pos=[x2+0.01,y0,x2+0.04,y2], $
    /right,ticknames=cbTickNames
  
  ; labels
  cgText,mean([x0,x2]),0.03,textoidl("r / r_{vir}"),/normal,alignment=0.5
  cgText,0.04,mean([y0,y2]),textoidl("Timescale [Gyr]"),/normal,alignment=0.5,orientation=90.0
  cgText,0.96,mean([y0,y2]),textoidl("Total Gas Mass [_{ }log h^{-1} M_{sun }]"),/normal,alignment=0.5,orientation=-90.0
  
  end_PS
  
end

; timescaleVRadStack(): 2D bin cooling/dynamical time ratio vs radial vel for 4 stacked halo mass bins
pro timescaleVRadStack

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  sP = simParams(res=512,run='tracer',redshift=2.0)
  gc = loadGroupCat(sP=sP,/skipIDs)
  h  = loadSnapshotHeader(sP=sP)
  
  massRanges = list([10.5,10.7], [10.9,11.1], [11.3,11.5], [11.7,11.9])
  
  ; plot config - 2x2 layout
  vradRange  = [-5.0,2.5]  ; vrad/v200
  ratioRange = [0.1,80.0] ; tcool/tdyn (log)
  
  x0 = 0.14 & x1 = 0.55 & x2 = 0.95
  y0 = 0.14 & y1 = 0.54 & y2 = 0.94
    
  pos = list( [x0,y1,x1,y2] ,$ ; upper left
              [x1,y1,x2,y2] ,$ ; upper right
              [x0,y0,x1,y1] ,$ ; lower left
              [x1,y0,x2,y1]  ) ; lower right
  
  ; start plot
  start_PS, sP.plotPath + 'timescales_vrad_2x2.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
      
  foreach massRange,massRanges,k do begin
    ; make list of halos
    gcIDList = gcIDList(gc=gc,select='pri',massRange=massRange)
                    
    print,k,n_elements(gcIDList)
    
    ; start plot
    if k eq 0 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=vradRange,yrange=ratioRange,/ylog,/xs,/ys,yminor=0,pos=pos[k],xtickname=replicate(' ',10),/noerase,xtickv=[-4,-3,-2,-1,0,1,2],xticks=6
    if k eq 1 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=vradRange,yrange=ratioRange,/ylog,/xs,/ys,yminor=0,pos=pos[k],xtickname=replicate(' ',10),ytickname=replicate(' ',10),/noerase,xtickv=[-4,-3,-2,-1,0,1,2],xticks=6
    if k eq 2 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=vradRange,yrange=ratioRange,/ylog,/xs,/ys,yminor=0,pos=pos[k],/noerase,xtickv=[-4,-3,-2,-1,0,1,2],xticks=6
    if k eq 3 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=vradRange,yrange=ratioRange,/ylog,/xs,/ys,yminor=0,pos=pos[k],ytickname=replicate(' ',10),/noerase,xtickv=[-4,-3,-2,-1,0,1,2],xticks=6
    
    ; load gas timescales and best model fits
    ts = loadFitTimescales(sP=sP,gcIDList=gcIDList)
    
    ; radial fit
    cgPlot,ts.radVrad.binCen,ts.radVrad.radMedian,line=2,color=cgColor('black'),/overplot
                
    ; set a contour palette
    loadColorTable,'brewerc-blues'
    tvlct, rr, gg, bb, /get
    palette = [[rr],[gg],[bb]]
    
    ; contour
    cgd = calcGridData(xx=ts.gasVrad/ts.gasVcirc,yy=ts.coolTime/ts.dynTime,$
                       xMinMax=vradRange*1.1,yMinMax=ratioRange*[0.8,1.1],nPixels=[60,60],/logY)
    
    hh = smooth(cgd.dens_out,[5,5])
    cgContour, hh/max(hh), cgd.xPts, 10.0^cgd.yPts, $
      /overplot, /fill, palette=palette, levels=[0.05,0.1,0.2,0.4,0.6,0.8,0.9],c_colors=(indgen(7)*20+50)
                
    ; mass label
    cgPlot,[0.0,0.0],[0.2,50.0],line=2,color=cgColor('black'),/overplot
    cgPlot,[-4.2,1.8],[1.0,1.0],line=2,color=cgColor('black'),/overplot
    
    cgText,vradRange[1]*0.5,ratioRange[0]*2,$
      "M = "+string(mean(massRange),format='(f4.1)'),alignment=0.5,/data
           
    ; redo clean plot lines
    loadct,0,/silent
    
    if k eq 0 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=vradRange,yrange=ratioRange,/ylog,/xs,/ys,yminor=0,pos=pos[k],xtickname=replicate(' ',10),/noerase,xtickv=[-4,-3,-2,-1,0,1,2],xticks=6
    if k eq 1 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=vradRange,yrange=ratioRange,/ylog,/xs,/ys,yminor=0,pos=pos[k],xtickname=replicate(' ',10),ytickname=replicate(' ',10),/noerase,xtickv=[-4,-3,-2,-1,0,1,2],xticks=6
    if k eq 2 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=vradRange,yrange=ratioRange,/ylog,/xs,/ys,yminor=0,pos=pos[k],/noerase,xtickv=[-4,-3,-2,-1,0,1,2],xticks=6
    if k eq 3 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=vradRange,yrange=ratioRange,/ylog,/xs,/ys,yminor=0,pos=pos[k],ytickname=replicate(' ',10),/noerase,xtickv=[-4,-3,-2,-1,0,1,2],xticks=6
    
  endforeach
  
  ; labels
  cgText,mean([x0,x2]),0.03,textoidl("v_{rad} / v_{circ}"),/normal,alignment=0.5
  cgText,0.04,mean([y0,y2]),textoidl("\tau_{cool} / \tau_{dyn}"),/normal,alignment=0.5,orientation=90.0
  
  end_PS

end

; timescaleAccTimeStack(): 2d bin cooling and dynamical times vs. acc times in mass bins

pro timescaleAccTimeStack
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  sP = simParams(res=512,run='tracer',redshift=2.0)
  gc = loadGroupCat(sP=sP,/skipIDs)
  h  = loadSnapshotHeader(sP=sP)
  massRanges = list([10.5,10.7], [10.9,11.1], [11.3,11.5], [11.7,11.9])
    
  yaxisWhat = 'tcool' ; tdyn,tcool,ratio,vrad
    
  ; median line
  fitRange  = [0.05,1.75]
  fitNBins  = 10
  fitSK     = 3
  
  ; calcGridData
  cgdRes = 40 ; NxN 2d sph mapping
  cgdSK  = 1  ; NxN post-smoothing
  
  ; plot config
  pConfig = hash()
  pConfig['tdyn']  = { yrange : [0.03,2.0] , ylog : 1, yaxisLabel : "\tau_{dyn} [Gyr]" }
  pConfig['tcool'] = { yrange : [0.03,2.0] , ylog : 1, yaxisLabel : "\tau_{cool} [Gyr]" }
  pConfig['ratio'] = { yrange : [0.1,80.0] , ylog : 1, yaxisLabel : "\tau_{cool} / \tau_{dyn}" }
  pConfig['vrad']  = { yrange : [-4.0,2.5] , ylog : 0, yaxisLabel : "v_{rad} / v_{circ}" }
  
  accTimeRange  = [0.03,2.0]
  
  ; 2x2 layout
  x0 = 0.14 & x1 = 0.55 & x2 = 0.95
  y0 = 0.14 & y1 = 0.54 & y2 = 0.94
    
  pos = list( [x0,y1,x1,y2] ,$ ; upper left
              [x1,y1,x2,y2] ,$ ; upper right
              [x0,y0,x1,y1] ,$ ; lower left
              [x1,y0,x2,y1]  ) ; lower right
  
  ; start plot
  start_PS, sP.plotPath+'tacc_'+yaxisWhat+'_2x2.'+str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
      
  foreach massRange,massRanges,k do begin
    ; make list of halos
    gcIDList = gcIDList(gc=gc,select='pri',massRange=massRange)

    ; start plot
    yrange = pConfig[yaxisWhat].yrange
    ylog   = pConfig[yaxisWhat].ylog
    yminor = !y.minor - !y.minor*(pConfig[yaxisWhat].ylog ne 0) ; 0 if log, !y.minor (2) otherwise
    
    if k eq 0 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=accTimeRange,yrange=yrange,/xs,/ys,$
      ylog=ylog,yminor=yminor,/xlog,xminor=0,pos=pos[k],xtickname=replicate(' ',10),/noerase
    if k eq 1 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=accTimeRange,yrange=yrange,/xs,/ys,$
      ylog=ylog,yminor=yminor,/xlog,xminor=0,pos=pos[k],xtickname=replicate(' ',10),$
      ytickname=replicate(' ',10),/noerase
    if k eq 2 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=accTimeRange,yrange=yrange,/xs,/ys,$
      ylog=ylog,yminor=yminor,/xlog,xminor=0,pos=pos[k],/noerase
    if k eq 3 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=accTimeRange,yrange=yrange,/xs,/ys,$
      ylog=ylog,yminor=yminor,/xlog,xminor=0,pos=pos[k],ytickname=replicate(' ',10),/noerase
    
    ; load gas timescales and best model fits
    ts = loadFitTimescales(sP=sP,gcIDList=gcIDList,/accTimesRepTR)
         
    print,'['+str(k)+'] nHalos = '+string(n_elements(gcIDList),format='(i4)')+' with avg nGas = '+$
      string(ts.nGas/n_elements(gcIDList),format='(i7)')
         
    ; set a contour palette
    loadColorTable,'brewerc-blues'
    tvlct, rr, gg, bb, /get
    palette = [[rr],[gg],[bb]]
    
    ; select yaxis data
    if yaxisWhat eq 'tdyn'  then yy = ts.dynTime
    if yaxisWhat eq 'tcool' then yy = ts.coolTime
    if yaxisWhat eq 'ratio' then yy = ts.coolTime/ts.dynTime
    if yaxisWhat eq 'vrad'  then yy = ts.gasVrad/ts.gasVcirc
    
    ; contour
    if  ylog then yMinMax = yrange * [0.8,1.1]
    if ~ylog then yMinMax = yrange * [1.1,1.1]
    
    cgd = calcGridData(xx=ts.accTime,yy=yy,xMinMax=accTimeRange*[0.8,1.1],yMinMax=yMinMax,$
                       nPixels=[cgdRes,cgdRes],logY=ylog,/logX)
    
    hh = smooth(cgd.dens_out,[cgdSK,cgdSK])
    hh = hh/max(hh) > 0.049
    
    yPts = cgd.yPts
    if ylog then yPts = 10.0^cgd.yPts
    cgContour, hh, 10.0^cgd.xPts, yPts, $
      /overplot, /fill, palette=palette, levels=[0.05,0.1,0.2,0.4,0.6,0.8,0.9],c_colors=(indgen(7)*20+50)
                
    ; mass label
    cgText,pos[k,2]-0.07,pos[k,1]+0.03,"M = "+string(mean(massRange),format='(f4.1)'),alignment=0.5,/normal
    
    ; median fit
    radFit = fitRadProfile(radii=alog10(ts.accTime),vals=yy,$
      range=alog10(fitRange),radBins=fitNBins)
    cgPlot,10.0^radFit.binCen,smooth(radFit.radMedian,fitSK),psym=-16,/overplot ; filled circle
    
    if  ylog then lineY = [0.04,1.76] ; 1-to-1 line
    if ~ylog then lineY = [0.0,0.0] ; zero line
    
    cgPlot,[0.04,1.76],lineY,line=2,color=cgColor('black'),/overplot ; 1-to-1 line
           
    ; redo clean plot lines
    loadct,0,/silent
    
    if k eq 0 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=accTimeRange,yrange=yrange,/xs,/ys,$
      ylog=ylog,yminor=yminor,/xlog,xminor=0,pos=pos[k],xtickname=replicate(' ',10),/noerase
    if k eq 1 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=accTimeRange,yrange=yrange,/xs,/ys,$
      ylog=ylog,yminor=yminor,/xlog,xminor=0,pos=pos[k],xtickname=replicate(' ',10),$
      ytickname=replicate(' ',10),/noerase
    if k eq 2 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=accTimeRange,yrange=yrange,/xs,/ys,$
      ylog=ylog,yminor=yminor,/xlog,xminor=0,pos=pos[k],/noerase
    if k eq 3 then $
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=accTimeRange,yrange=yrange,/xs,/ys,$
      ylog=ylog,yminor=yminor,/xlog,xminor=0,pos=pos[k],ytickname=replicate(' ',10),/noerase
    
  endforeach
  
  ; labels
  cgText,mean([x0,x2]),0.03,textoidl("\tau_{acc} [Gyr]"),/normal,alignment=0.5
  cgText,0.04,mean([y0,y2]),textoidl(pConfig[yaxisWhat].yaxisLabel),/normal,alignment=0.5,orientation=90.0
  
  end_PS
  
  stop
end

@akde ; required below

; timescaleVRadSlice(): slice timescaleVradStack in tsRatio (fixed Mhalo) and in Mhalo (at fixed tsRatio)
pro timescaleVRadSlice

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  sP = simParams(res=512,run='tracer',redshift=2.0)
  gc = loadGroupCat(sP=sP,/skipIDs)
  h  = loadSnapshotHeader(sP=sP)
  
  massRanges = list([10.5,10.7], [10.9,11.1], [11.3,11.5], [11.7,11.9])
  tsRanges   = list([0.1,0.5], [0.6, 0.9], [2.0, 5.0], [6.0, 10.0], [10.0, 30.0])
  
  colors = ['red','blue','orange','forest green','purple']
  
  massRangeTS  = [0.1,0.9] ; timescale ratio to take halo mass slices at
  tsRangeMhalo = [11.3,11.5] ; log msun of halo mass to take tsRange slices at
  
  ; plot config - 2x2 layout
  vradRange  = [-5.0,2.5]  ; vrad/v200
  ratioRange = [0.1,80.0] ; tcool/tdyn (log)
  
  fracCutoff = 0.01 ; don't plot below this
  
  ; start plot
  start_PS,$
    sP.plotPath+'timescales_vrad_slice.'+str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)+'.eps',$
    xs=8,ys=8
  
  pos = list([0.15,0.56,0.95,0.95], $
             [0.15,0.15,0.95,0.54])
  
  ; timescale ratio slices
  cgPlot,[0],[0],/nodata,xtitle="",ytitle="",xrange=vradRange,yrange=[0.0,0.5],$
    /xs,/ys,position=pos[0],xtickname=replicate(' ',10)
  
  gcIDList = where(codeMassToLogMsun(gc.subgroupMass) ge tsRangeMhalo[0] and $
                   codeMassToLogMsun(gc.subgroupMass) lt tsRangeMhalo[1],count)
                   
  strings = []
  
  ; load gas timescales
  ts = loadFitTimescales(sP=sP,gcIDList=gcIDList)  
  
  tsRatio = ts.coolTime / ts.dynTime
  vRadNorm = ts.gasVrad / ts.gasVcirc
  
  foreach tsRange,tsRanges,k do begin
    ; restrict to timescale slice
    w = where(tsRatio ge tsRange[0] and tsRatio lt tsRange[1],count)
    print,k,tsRange[0],tsRange[1],count
    
    ; histogram and plot
    ;binsize = 10.0 / sqrt(count) / (sP.res/128)
    ;h = histogram(vRadNorm[w],binsize=binsize,loc=loc)
    ;cgPlot,loc+binsize*0.5,float(h)/max(h)/2,color=colors[k],line=2,/overplot
    
    ; KDE method
    xpts = linspace(vradRange[0],vradRange[1],200)
    dd = kde(vRadNorm[w],xpts)
    
    w = where(dd ge fracCutoff)
    cgPlot,xpts[w],dd[w],color=colors[k],line=0,/overplot
    
    ; AKDE
    ;dd = akde(vRadNorm[w],xpts)
    
    ;w = where(dd ge fracCutoff)
    ;cgPlot,xpts[w],dd[w],color=colors[k],line=0,/overplot
    
    ; save for legend
    strings = [strings,string(tsRange[0],format='(f4.1)') + ' ' + $
               textoidl('< \tau_{c} / \tau_{d} <') + ' ' + string(tsRange[1],format='(f4.1)')]
  endforeach
  
  legend,strings,textcolors=colors,box=0,/top,/right
  cgText,-3.6,0.43,string(tsRangeMhalo[0],format='(f4.1)') + ' ' + $
               textoidl('< M_{halo} <') + ' ' + string(tsRangeMhalo[1],format='(f4.1)'),/data,alignment=0.5
  
  ; halo mass slices
  cgPlot,[0],[0],/nodata,xtitle=textoidl('v_{rad} / v_{circ}'),ytitle="",$
    xrange=vradRange,yrange=[0.0,0.5],/xs,/ys,position=pos[1],/noerase
  
  strings = []
  
  foreach massRange,massRanges,k do begin
    ; make list of halos
    gcIDList = where(codeMassToLogMsun(gc.subgroupMass) ge massRange[0] and $
                     codeMassToLogMsun(gc.subgroupMass) lt massRange[1],count)
                    
    print,k,count
    
    ; load gas timescales and best model fits
    ts = loadFitTimescales(sP=sP,gcIDList=gcIDList)

    tsRatio = ts.coolTime / ts.dynTime
    vRadNorm = ts.gasVrad / ts.gasVcirc
    
    ; restrict to timescale slice
    w = where(tsRatio ge massRangeTS[0] and tsRatio lt massRangeTS[1],count)
    
    ; KDE method
    xpts = linspace(vradRange[0],vradRange[1],200)
    dd = kde(vRadNorm[w],xpts)
    
    w = where(dd ge fracCutoff)
    cgPlot,xpts[w],dd[w],color=colors[k],line=0,/overplot
    
    ; save for legend
    strings = [strings,string(massRange[0],format='(f4.1)') + ' ' + $
               textoidl('< M_{halo} <') + ' ' + string(massRange[1],format='(f4.1)')]
    
  endforeach
  
  legend,strings,textcolors=colors,box=0,/top,/right
  cgText,-3.7,0.43,string(massRangeTS[0],format='(f4.1)') + ' ' + $
               textoidl('< \tau_c / \tau_d <') + ' ' + string(massRangeTS[1],format='(f4.1)'),/data,alignment=0.5

  cgText,0.05,mean([0.15,0.95]),"Gas Fraction",/normal,alignment=0.5,orientation=90.0
  
  end_PS

end

; compareTimescalesHalo(): explore for a single halo or a mass range stacked

pro compareTimescalesHalo

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  sP = simParams(res=256,run='tracer',redshift=2.0)
  gc = loadGroupCat(sP=sP,/skipIDs)

  ; do for tracers (and include tacc comparison plots) or instead do for gas?
  accTimesRepTR = 1
  
  ; single halo by haloID
  haloID = 304 ;z2.304 z2.301 z2.130 z2.64
  gcID = getMatchedIDs(sPa=sP,sPg=sP,haloID=haloID)
  gcIDList = [gcID.a]
  hTag = 'h'+str(haloID)+'_m='+string(codeMassToLogMsun(gc.subgroupMass[gcID.a]),format='(f4.1)')+$
         '_tr'+str(accTimesRepTR)

  ; load gas timescales and best model fits
  ts = loadFitTimescales(sP=sP,gcIDList=gcIDList,accTimesRepTR=accTimesRepTR)
  
  ; plot config
  coolingRange = [0.03,20]  ; Gyr
  dynRange     = [0.03,2.0] ; Gyr
  accRange     = dynRange
  tempRange    = [4.25,6.5] ; log K
  densRange    = [0.1,2000] ; log ratioToCrit
  radRange     = [0.05,1.6] ; r/rvir
  vradRange    = [-5.0,3.0] ; vrad/vcirc
  ratioRange   = [0.01,80.0] ; tcool/tdyn (log)
  
  binsize = 0.1 / (sP.res/128)
  
  psym = 4
  
  ; 2d config
  do2DBin   = sP.res gt 128 ; 0=plot individual gas elements, 1=do 2d binning
  binSizeXY = [0.06,0.2,0.16,0.24] ; rad,timescale,temp,rhoRatio
  hsp       = [0.005,0.02] ; rad, gyr
  nc        = 120 ; number of colors (of 255) to use for background 2d histo
  
  ; plot (1) - cooling time vs current temperature scatter
  start_PS, $
    sP.plotPath + 'timescales_vs_tcool.'+hTag+'.'+str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)+'.eps', ys=8, xs=6
    cgPlot,ts.coolTime,ts.curTemp,psym=psym,xtitle="",ytitle=textoidl('T_{gas} [log K]'),$
      xrange=coolingRange,yrange=tempRange,/xlog,/xs,/ys,xminor=0,position=(sP.pos_3x1)[0],xtickname=replicate(' ',10)

    cgPlot,ts.coolTime,alog10(rhoRatioToCrit(ts.curDens,redshift=sP.redshift)),psym=psym,xtitle="",$
      ytitle=textoidl('log( \rho_{gas} / \rho_{crit,z} )'),xtickname=replicate(' ',10),$
      xrange=coolingRange,yrange=alog10(densRange),/xlog,/xs,/ys,xminor=0,/noerase,position=(sP.pos_3x1)[1]
      
    cgPlot,[0],[0],/nodata,xtitle=textoidl('\tau_{cool} [Gyr]'),ytitle=textoidl('\tau_{dyn} [Gyr]'),$
      xrange=coolingRange,yrange=dynRange,/ys,/ylog,yminor=0,/xs,/xlog,xminor=0,position=(sP.pos_3x1)[2],/noerase
    cgPlot,[0.04,1.8],[0.04,1.8],line=0,color=cgColor('light gray'),/overplot
    cgPlot,ts.coolTime,ts.dynTime,psym=psym,/overplot
  end_PS
  
  ; plot (2) - gas density/temperature vs radius
  start_PS,$
    sP.plotPath+'timescales_vs_rad.'+hTag+'.'+str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)+'.eps',xs=6,ys=10
    ; timescales
    ; ----------
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="Timescale [Gyr]",$
      xrange=radRange,yrange=coolingRange+[0.0,10.0],/ylog,/xs,/ys,yminor=0,$
      xtickv=[0.15,0.5,1.0,1.5],xticks=3,xtickname=replicate(' ',10),position=(sP.pos_3x1)[0]
      
    ; individual gas elements
    if ~keyword_set(do2DBin) then begin
      cgPlot,ts.gasRadii/ts.gasRvir,ts.coolTime,psym=psym,color=cgColor('light gray'),/overplot
    endif else begin
      f2d = binHisto2D(xx=ts.gasRadii/ts.gasRvir, yy=alog10(ts.coolTime), $
                       xmm=radRange, ymm=alog10(coolingRange), xbs=binSizeXY[0], ybs=binSizeXY[1])        
      oplot2DHistoSq, f2d, hsp=hsp, nc=nc, xrange=radRange, yrange=coolingRange, /colNorm, /logY, /gray
    endelse
    
    ; dynTime(r) and hubbleTime
    ;cgPlot,ts.gasRadii/ts.gasRvir,ts.dynTime,psym=psym,line=2,color=cgColor('dark gray'),/overplot
    cgPlot,[0.2,1.45],[ts.hubbleTime,ts.hubbleTime],line=1,color=cgColor('forest green'),/overplot
     
    ; radial fits
    cgPlot,ts.radCt.binCen,ts.radCt.radMedian,line=0,psym=-4,color=cgColor('black'),/overplot
    cgPlot,ts.radDt.binCen,ts.radDt.radMedian,line=2,psym=-4,color=cgColor('black'),/overplot
     
    ; SIS/NFW models
    cgPlot,ts.x,ts.sis_gas.coolTime,line=0,color=cgColor('orange'),/overplot
    cgPlot,ts.x,ts.sis_gas.dynTime,line=2,color=cgColor('orange'),/overplot
    
    cgPlot,ts.x,ts.nfw_gas.dynTime,line=2,color=cgColor('red'),/overplot
    cgPlot,ts.x,ts.nfw_gas.coolTime,line=0,color=cgColor('red'),/overplot
     
    ; mark various r_cool
    yr_rcool = [0.04,0.1]
    
    cgPlot,[ts.sis_gas.r_cool_hubble/ts.sis_dm.r200,ts.sis_gas.r_cool_hubble/ts.sis_dm.r200],yr_rcool,$
      line=0,color=cgColor('forest green'),/overplot
    cgPlot,[ts.sis_gas.r_cool/ts.sis_dm.r200,ts.sis_gas.r_cool/ts.sis_dm.r200],yr_rcool,$
      line=0,color=cgColor('orange'),/overplot
    cgPlot,[ts.nfw_gas.r_cool/ts.nfw_dm.r200,ts.nfw_gas.r_cool/ts.nfw_dm.r200],yr_rcool,$
      line=0,color=cgColor('red'), /overplot
     
    ; legend
    legend,textoidl(['\tau_{dyn}','\tau_{cool}']),linestyle=[2,0],box=0,/bottom,/right,linesize=0.2
    cgText,0.2,6.0,textoidl(['t_{age}']),color=cgColor('forest green'),/data,alignment=0.5
    
    ; temp
    ; ----
    cgPlot,[0],[0],/nodata,xtitle="",ytitle=textoidl('T_{gas} [log K]'),xrange=radRange,yrange=tempRange+[0.1,0.0],$
      /xs,/ys,position=(sP.pos_3x1)[1],xtickv=[0.15,0.5,1.0,1.5],xticks=3,xtickname=replicate(' ',10),/noerase
      
    ; individual gas elements and median binned
    if ~keyword_set(do2DBin) then begin
      cgPlot,ts.gasRadii/ts.gasRvir,ts.curTemp,psym=psym,color=cgColor('light gray'),/overplot
    endif else begin
      f2d = binHisto2D(xx=ts.gasRadii/ts.gasRvir, yy=ts.curTemp, $
                       xmm=radRange, ymm=tempRange, xbs=binSizeXY[0], ybs=binSizeXY[2])        
      oplot2DHistoSq, f2d, hsp=hsp, nc=nc, xrange=radRange, yrange=coolingRange, /colNorm, /gray
    endelse
    
    cgPlot,ts.radTemp.binCen,ts.radTemp.radMedian,psym=-4,color=cgColor('black'),/overplot
    
    ; SIS model and SIS fit
    cgPlot,ts.x,alog10(ts.sis_gas.temp_gas),line=0,color=cgColor('orange'),/overplot
    cgPlot,ts.x,ts.sis_fit.f_temp,line=2,color=cgColor('orange'),/overplot
    
    ; NFW iso
    cgPlot,ts.x,alog10(replicate(ts.nfw_gas.T_0,n_elements(ts.x))),line=0,color=cgColor('blue'),/overplot
    
    ; NFW model and NFW fit (poly)
    cgPlot,ts.x,alog10(ts.nfw_gas.temp_gas),line=0,color=cgColor('red'),/overplot    
    cgPlot,ts.x,ts.nfw_fit.f_temp,line=2,color=cgColor('red'),/overplot
    
    ; dens
    ; ----
    cgPlot,[0],[0],/nodata,xtitle=textoidl("r / r_{vir}"),ytitle=textoidl('log( \rho_{gas} / \rho_{crit,z} )'),$
      xrange=radRange,yrange=alog10(densRange),/xs,/ys,/noerase,position=(sP.pos_3x1)[2],$
      xtickv=[0.15,0.5,1.0,1.5],xticks=3,xtickname=['0.15','0.5','1.0','1.5']
      
    ; individual gas elements and median binned
    if ~keyword_set(do2DBin) then begin
      cgPlot,ts.gasRadii/ts.gasRvir,alog10(rhoRatioToCrit(ts.curDens,sP=sP)),$
        psym=psym,color=cgColor('light gray'),/overplot
    endif else begin
      f2d = binHisto2D(xx=ts.gasRadii/ts.gasRvir, yy=alog10(rhoRatioToCrit(ts.curDens,sP=sP)), $
                       xmm=radRange, ymm=alog10(densRange), xbs=binSizeXY[0], ybs=binSizeXY[3])        
      oplot2DHistoSq, f2d, hsp=hsp, nc=nc, xrange=radRange, yrange=alog10(densRange), /colNorm, /gray
    endelse

    cgPlot,ts.radDens.binCen,alog10(rhoRatioToCrit(ts.radDens.radMedian,sP=sP)),psym=-4,color=cgColor('black'),/overplot
    
    ; SIS model and SIS fit
    cgPlot,ts.x,alog10(rhoRatioToCrit(ts.sis_gas.rho_gas,sP=sP)),line=0,color=cgColor('orange'),/overplot
    cgPlot,ts.x,alog10(rhoRatioToCrit(ts.sis_fit.f_dens,sP=sP)),line=2,color=cgColor('orange'),/overplot
    
    ; NFW iso dens
    cgPlot,ts.x,alog10(rhoRatioToCrit(ts.nfw_gas.rho_gas_iso,sP=sP)),line=0,color=cgColor('blue'),/overplot
    
    ; NFW model and NFW fit (poly)
    cgPlot,ts.x,alog10(rhoRatioToCrit(ts.nfw_gas.rho_gas,sP=sP)),line=0,color=cgColor('red'),/overplot
    cgPlot,ts.x,alog10(rhoRatioToCrit(ts.nfw_fit.f_dens,sP=sP)),line=2,color=cgColor('red'),/overplot
    
    legend,['median','SIS','NFW poly','NFW iso'],textcolors=['black','orange','red','blue'],$
      box=0,/bottom,/left
    legend,['model','fit'],linestyle=[0,2],linesize=0.2,box=0,/top,/right
    
  end_PS
  
  ; accretion time plots
  ; --------------------
  if accTimesRepTR then begin
    binSizeXY = [0.08,0.16]
    do2DBin = 0
    
    ; plot (1) - cooling, dynamical, and ratio vs. tacc
    start_PS, $
    sP.plotPath + 'timescales_vs_tacc.'+hTag+'.'+str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)+'.eps', ys=8,xs=6
    
    if sP.res eq 128 then pThick = !p.thick else pThick = 1.0

    ; top: cooling
    cgPlot,[0],[0],/nodata,xtitle="",ytitle=textoidl('\tau_{cool} [Gyr]'),$
      xrange=accRange,yrange=coolingRange,/xlog,/xs,/ys,xminor=0,/ylog,yminor=0,$
      position=(sP.pos_3x1)[0],xtickname=replicate(' ',10)
    
    if ~keyword_set(do2DBin) then begin
      cgPlot,ts.accTime,ts.coolTime,psym=psym,/overplot,thick=pThick
    endif else begin
      f2d = binHisto2D(xx=alog10(ts.accTime), yy=alog10(ts.coolTime), $
                       xmm=alog10(accRange), ymm=alog10(coolingRange), $
                       xbs=binSizeXY[0], ybs=binSizeXY[1])        
      oplot2DHistoSq, f2d, hsp=hsp, nc=nc, xrange=accRange, yrange=coolingRange, /logX, /logY, /gray
    endelse
    
    cgPlot,[0.04,1.8],[0.04,1.8],line=0,color=cgColor('orange'),/overplot
    
    ; middle: dynamical
    cgPlot,[0],[0],/nodata,xtitle="",ytitle=textoidl('\tau_{dyn} [Gyr]'),$
      xrange=accRange,yrange=dynRange,/xlog,/xs,/ys,xminor=0,/ylog,yminor=0,$
      /noerase,position=(sP.pos_3x1)[1],xtickname=replicate(' ',10)

    ; set a contour palette
    loadColorTable,'brewerc-blues'
    tvlct, rr, gg, bb, /get
    palette = [[rr],[gg],[bb]]
    
    ; contour
    cgd = calcGridData(xx=ts.accTime,yy=ts.dynTime,$
                       xMinMax=accRange*[0.8,1.1],yMinMax=dynRange*[0.8,1.1],$
                       nPixels=[40,40],/logY,/logX)
    
    hh = smooth(cgd.dens_out,[3,3])
    cgContour, hh/max(hh) > 0, 10.0^cgd.xPts, 10.0^cgd.yPts, $
      /overplot, /fill, palette=palette, levels=[0.05,0.1,0.2,0.4,0.6,0.8,0.9],c_colors=(indgen(7)*20+50)
      
    ; individual elements
    if ~keyword_set(do2DBin) then begin
      cgPlot,ts.accTime,ts.dynTime,psym=psym,/overplot,thick=pThick
    endif else begin
      f2d = binHisto2D(xx=alog10(ts.accTime), yy=alog10(ts.dynTime), $
                       xmm=alog10(accRange), ymm=alog10(dynRange), $
                       xbs=binSizeXY[0], ybs=binSizeXY[1])        
      oplot2DHistoSq, f2d, hsp=hsp, nc=nc, xrange=accRange, yrange=dynRange, /logX, /logY, /gray
    endelse
    
    cgPlot,[0.04,1.8],[0.04,1.8],line=0,color=cgColor('orange'),/overplot
    
    ; 2tdyn envelope (no accTime greater than this can be found)
    cgPlot,[0.08,1.8],0.5*[0.08,1.8],line=2,color=cgColor('orange'),/overplot
      
    ; bottom: ratio
    cgPlot,[0],[0],/nodata,xtitle=textoidl('\tau_{acc} [Gyr]'),ytitle=textoidl('\tau_{cool}/\tau_{dyn}'),$
      xrange=accRange,yrange=ratioRange,/ys,/xs,/xlog,xminor=0,/ylog,yminor=0,$
      position=(sP.pos_3x1)[2],/noerase
      
    ; individual elements
    if ~keyword_set(do2DBin) then begin
      cgPlot,ts.accTime,ts.coolTime/ts.dynTime,psym=psym,/overplot,thick=pThick
    endif else begin
      f2d = binHisto2D(xx=alog10(ts.accTime), yy=ts.coolTime/ts.dynTime, $
                       xmm=alog10(accRange), ymm=ratioRange, xbs=binSizeXY[0], ybs=binSizeXY[1])        
      oplot2DHistoSq, f2d, hsp=hsp, nc=nc, xrange=accRange, yrange=ratioRange, /logX, /gray
    endelse
    
    cgPlot,[0.04,1.8],[1.0,1.0],line=0,color=cgColor('orange'),/overplot
    
    end_PS
    
    ; plot (2) - radius, vrad vs. tacc
    start_PS, $
    sP.plotPath + 'timescales_vs_tacc2.'+hTag+'.'+str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)+'.eps', ys=8,xs=6
    
    if sP.res eq 128 then pThick = !p.thick else pThick = 1.0

    ; top: radius
    cgPlot,[0],[0],/nodata,xtitle="",ytitle=textoidl('r / r_{vir}'),$
      xrange=accRange,yrange=radRange,/xlog,/xs,/ys,xminor=0,$
      position=(sP.pos_3x1)[0],xtickname=replicate(' ',10),$
      ytickv=[0.15,0.5,1.0,1.5],yticks=3,ytickname=['0.15','0.5','1.0','1.5']
    
    if ~keyword_set(do2DBin) then begin
      cgPlot,ts.accTime,ts.gasRadii/ts.gasRvir,psym=psym,/overplot,thick=pThick
    endif else begin
      f2d = binHisto2D(xx=alog10(ts.accTime), yy=ts.gasRadii/ts.gasRvir, $
                       xmm=alog10(accRange), ymm=radRange, xbs=binSizeXY[0], ybs=binSizeXY[1])        
      oplot2DHistoSq, f2d, hsp=hsp, nc=nc, xrange=accRange, yrange=radRange, /logX, /gray
    endelse
    
    ; middle: vrad
    cgPlot,[0],[0],/nodata,xtitle="",ytitle=textoidl('v_{rad} / v_{circ}'),$
      xrange=accRange,yrange=vradRange,/xlog,/xs,/ys,xminor=0,$
      /noerase,position=(sP.pos_3x1)[1],xtickname=replicate(' ',10)

    ; individual elements
    if ~keyword_set(do2DBin) then begin
      cgPlot,ts.accTime,ts.gasVrad/ts.gasVcirc,psym=psym,/overplot,thick=pThick
    endif else begin
      f2d = binHisto2D(xx=alog10(ts.accTime), yy=ts.gasVrad/ts.gasVcirc, $
                       xmm=alog10(accRange), ymm=vradRange, $
                       xbs=binSizeXY[0], ybs=binSizeXY[1])        
      oplot2DHistoSq, f2d, hsp=hsp, nc=nc, xrange=accRange, yrange=vradRange, /logX, /gray
    endelse
    
    cgPlot,[0.04,1.8],[0.0,0.0],line=0,color=cgColor('orange'),/overplot
      
    ; bottom: r / vrad
    rvradRatioRange = [-20.0,5.0]
    cgPlot,[0],[0],/nodata,xtitle=textoidl('\tau_{acc} [Gyr]'),$
      ytitle=textoidl('(v_{rad}/v_{circ}) / (r/r_{vir})'),$
      xrange=accRange,yrange=rvradRatioRange,/ys,/xs,/xlog,xminor=0,position=(sP.pos_3x1)[2],/noerase
      
    yy = (ts.gasVrad/ts.gasVcirc) / (ts.gasRadii/ts.gasRvir)
    
    ; individual elements
    if ~keyword_set(do2DBin) then begin
      cgPlot,ts.accTime,yy,psym=psym,/overplot,thick=pThick
    endif else begin
      f2d = binHisto2D(xx=alog10(ts.accTime), yy=yy, $
                       xmm=alog10(accRange), ymm=vradRange, $
                       xbs=binSizeXY[0], ybs=binSizeXY[1])        
      oplot2DHistoSq, f2d, hsp=hsp, nc=nc, xrange=accRange, yrange=rvradRatioRange, /logX, /gray
    endelse
    
    cgPlot,[0.04,1.8],[1.0,1.0],line=0,color=cgColor('orange'),/overplot
    
    end_PS
  
  endif
  
  stop
  
end

; compCoolTimes(): KWH primordial network vs. SD93 tables

pro compCoolTimes
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; simulation (KWH)
  sP = simParams(res=128,run='tracer',redshift=2.0)
  ct = coolingTime(sP=sP)

  ; load, convert -> physical
  rho = codeDensToPhys(ct.dens, sP=sP, /cgs) ; cgs
  
  ; load galaxy/group member catalogs for gas ids to search for
  h = loadSnapshotHeader(sP=sP)
  galcat = galaxyCat(sP=sP)
  
  ; load gas ids and match to catalog
  ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')

  calcMatch,galcat.groupmemIDs,ids,galcat_ind,ids_gmem_ind,count=countGmem
  ids_gmem_ind = ids_gmem_ind[calcSort(galcat_ind)]
  
  ids = !NULL
  galcat_ind = !NULL

  u = loadSnapshotSubset(sP=sP,partType='gas',field='u')
  u = u[ids_gmem_ind]
  u *= units.UnitPressure_in_cgs / units.UnitDensity_in_cgs
  
  nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
  nelec = nelec[ids_gmem_ind]
  ids_gmem_ind = !NULL
  
  ; SD93 interp
  tables = interpLambdaSD93()
  lambdaNetNorm = 10.0^(interpLambdaSD93(Z=0.0,logT=ct.temp,tables=tables,/norm))
  
  ; ratefact (used in code for KWH CoolingRate())
  nHcgs = units.hydrogen_massfrac * rho / units.mass_proton; /* hydrogen number dens in cgs units */
  ratefact = nHcgs * nHcgs / rho;
  
  ; SD93 ratefact for lambda_N
  n_e = nelec * nHcgs ; nelec is relative to hydrogen number density
  n_t = rho / units.mass_proton
  ratefact2 = n_t * n_e / rho
  
  sd_coolTime = u / (rateFact2 * lambdaNetNorm)
  sd_coolTime *= units.HubbleParam / units.UnitTime_in_s ; Gyr
  
  ; plot histos
  xrange = [0.01,100]
  bin = 0.1
  
  start_PS,'comp_ct_hist.eps'
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=[0.01,1.1],/xs,/ys,/ylog,/xlog,xtitle="Cooling Time [Gyr]",ytitle="N / max"
    
    h = histogram(alog10(ct.coolTime),binsize=bin,min=alog10(xrange[0]),loc=loc)
    cgPlot,10.0^(loc+bin*0.5),float(h)/max(h),psym=-cgSymCat('filled circle'),color=cgColor('red'),/overplot
    
    h = histogram(alog10(sd_coolTime),binsize=bin,min=alog10(xrange[0]),loc=loc)
    cgPlot,10.0^(loc+bin*0.5),float(h)/max(h),psym=-8,color=cgColor('blue'),/overplot
    
    legend,['KWH','SD'],textcolors=['red','blue'],box=0,/top,/right
  end_PS
  
  start_PS,'comp_ct_vs.eps', xs=7, ys=5
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=xrange,/xs,/ys,xtitle="SD CT [Gyr]",ytitle="KWH CT [Gyr]", /xlog, /ylog
    cgPlot,xrange,xrange,line=0,color=cgColor('orange'),/overplot
    cgPlot,sd_coolTime,ct.coolTime,psym=3,/overplot
  end_PS
  
  start_PS,'comp_ct_vs2.eps'
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=[0.6,1.6],/xs,/ys,xtitle="SD CT [Gyr]",ytitle="KWH CT / SD CT", /xlog
    cgPlot,xrange,[1.0,1.0],line=0,color=cgColor('orange'),/overplot
    cgPlot,sd_coolTime,ct.coolTime/sd_coolTime,psym=3,/overplot
  end_PS
  
  ; compare lambdaNet
  logT = linspace(4.0,8.0,200)
  ddd = fltarr(200) + 1.0 ; dens
  nnn = fltarr(200) + 1.0 ; nelec
  
  lambdaNet_KWH = CalcCoolTime(logT,ddd,nnn,scalefac=0.333,flag=3)
  lambdaNet_SD  = interpLambdaSD93(Z=0.0,logT=logT,tables=tables,/norm)
  lambdaNet_SD  = alog10(abs(nnn * 10.0^lambdaNet_SD)) ; do n_e normalization
  
  start_PS,'comp_lambdaNet.eps'
    cgPlot,[0],[0],/nodata,xrange=[4.0,8.0],yrange=[-24.0,-21.0],/xs,/ys,xtitle="log T",ytitle="Lambda / n^2"
    
    cgPlot,logT,alog10(abs(lambdaNet_KWH)),color=cgColor('red'),/overplot
    cgPlot,logT,lambdaNet_SD,color=cgColor('blue'),/overplot
    
    legend,['KWH','SD'],textcolors=['red','blue'],box=0,/top,/right
  end_PS
  
  stop
  
end

; reesOstrikerFig1(): recreate Fig.1 of Rees & Ostriker (1977) to verify our cooling time calculations

pro reesOstrikerFig1

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  tempMinMax = [2.0,9.0] ; log(K)
  densMinMaxPlot = [-7,3] ; log(cm^-3)
  densMinMaxCalc = [-5,3] ; log(cm^-3)
  nTemps = 200
  nDens  = 200
  
  redshifts = [0.0,1.0,2.0]
  colors = ['red','blue','green']
  
  ; start plot
  start_PS,'fig1.eps'
    cgPlot,[0],[0],/nodata,xrange=10.0^densMinMaxPlot,yrange=10.0^tempMinMax,xtitle=textoidl("n [cm^{-3}]"),$
      ytitle="Temp [K]",/xlog,/ylog,/xs,/ys,xminor=1,yminor=1
  
  ; generate arrays
  temps = 10.0^(findgen(nTemps)/(nTemps-1) * (tempMinMax[1]-tempMinMax[0]) + tempMinMax[0]) ; K
  dens  = 10.0^(findgen(nDens)/(nDens-1) * (densMinMaxCalc[1]-densMinMaxCalc[0]) + densMinMaxCalc[0]) ; 1/cm^3
  dens_gcm3 = dens * units.mass_proton * units.HubbleParam * units.HubbleParam ; g/cm^3
  dens_code = float(dens_gcm3 / units.UnitDensity_in_cgs) ; code units (10^10 msun/kpc^3)
  nelec = fltarr(nDens) + 1.0

  ; constant jeans mass lines (M_jeans = 1e8 * (temps/1e4)^(1.5) * dens^(-0.5) ; Msun)
  M_jeans_targets = [1e12,1e11,1e9]
  temp_jeans0 = ( (M_jeans_targets[0]/1e8) / dens^(-0.5) ) ^ (2.0/3.0) * 1e4
  temp_jeans1 = ( (M_jeans_targets[1]/1e8) / dens^(-0.5) ) ^ (2.0/3.0) * 1e4
  temp_jeans2 = ( (M_jeans_targets[2]/1e8) / dens^(-0.5) ) ^ (2.0/3.0) * 1e4
  
  ; plot: constant jeans mass lines
  w = where(dens ge 10.0^(-5.5))
  cgPlot,dens[w],temp_jeans0[w],line=2,/overplot
  cgPlot,dens[w],temp_jeans1[w],line=2,/overplot
  cgPlot,dens[w],temp_jeans2[w],line=2,/overplot
  cgText,10^(-5.9),10^(5.2),textoidl('10^{12} M_{sun}'),/data,alignment=0.5
  cgText,10^(-5.5),10^(3.7),textoidl('10^{11} M_{sun}'),/data,alignment=0.5
  cgText,10^(-4.0),10^(2.5),textoidl('M_J=10^{9} M_{sun}'),/data,alignment=0.5
  
  ; ABC regions
  cgText,10^(-4.0),1e8,"A",/data,alignment=0.5
  cgText,10^(-0.6),1e8,"B",/data,alignment=0.5
  cgText,10^(1.0),10^(5.5),"C",/data,alignment=0.5
  
  foreach redshift,redshifts,m do begin
  
    ; age of universe
    scalefac = 1.0/(1+redshift)
    print,scalefac
  
    tage = redshiftToAgeFlat(redshift)

    ; tcool = tage and tcool = tdyn lines
    dens_tcool_tage = fltarr(nTemps)
    dens_tcool_tdyn = fltarr(nTemps)
  
    ; estimate the dynamical (free-fall) timescale for these densities
    tdyn_eval = sqrt( 3*!pi / (32 * float(units.G) * dens_code/scalefac^3.0) ) ; code units (Gyr)
  
    for i=0,nTemps-1 do begin
      ; for each temperature, calculate the tcool across all densities
      temps_loc = replicate(alog10(temps[i]),nDens)
    
      tcool_eval = calcCoolTime(temps_loc,dens_gcm3,nelec,scalefac=scalefac,flag=1)
      tcool_eval *= float(units.HubbleParam / units.UnitTime_in_s) ; convert cgs -> code units (Gyr)
  
      ; interpolate the evaluated tcool to tage and save 
      dens_tcool_tage[i] = interpol(dens,tcool_eval-tage,0.0)

      ; interpolate the evaluate tcool to tdyn and save
      dens_tcool_tdyn[i] = interpol(dens,tcool_eval-tdyn_eval,0.0)
      
      ;if i eq 150 then stop
    endfor
      
    ; plot tcool=tage/tdyn lines
    w = where(temps ge 10.0^(5.1))
    cgPlot,dens_tcool_tage[w],temps[w],line=1,/overplot,color=cgColor(colors[m])
    cgPlot,dens_tcool_tdyn[w],temps[w],line=0,/overplot,color=cgColor(colors[m])
  
  endforeach ;redshifts

  legend,'z='+string(redshifts,format='(f3.1)'),textcolors=colors,box=0,/top,/left
  
  ; end plot
  end_PS
  
  ; for a fixed temperature, plot the cooling time vs density at different redshifts
  redshifts = [0.0,0.5,1.0,1.5,2.0]
  colors = ['red','blue','green','magenta','orange']
  log_temp = [7.0,8.0]
  
  start_PS,'fig1_b.eps'
    cgPlot,[0],[0],/nodata,xrange=10.0^densMinMaxPlot,yrange=[0.01,1e4],xtitle=textoidl("n [cm^{-3}]"),$
      ytitle="Cooling/Dynamical Time [Gyr]",/xlog,/ylog,/xs,/ys,xminor=1,yminor=1
      
  foreach redshift,redshifts,m do begin
    scalefac = 1.0/(1+redshift)
    
    ; temp 1
    temps_loc  = replicate(log_temp[0],nDens)
    tcool_eval = calcCoolTime(temps_loc,dens_gcm3,nelec,scalefac=scalefac,flag=1)
    tcool_eval *= float(units.HubbleParam / units.UnitTime_in_s)
     
    cgPlot,dens,tcool_eval,color=cgColor(colors[m]),/overplot
    
    ; temp 2
    temps_loc  = replicate(log_temp[1],nDens)
    tcool_eval = calcCoolTime(temps_loc,dens_gcm3,nelec,scalefac=scalefac,flag=1)
    tcool_eval *= float(units.HubbleParam / units.UnitTime_in_s)
     
    cgPlot,dens,tcool_eval,color=cgColor(colors[m]),line=2,/overplot
    
    ; dynamical
    tdyn_eval = sqrt( 3*!pi / (32 * float(units.G) * dens_code/scalefac^3.0) ) * float(units.HubbleParam)
    cgPlot,dens,tdyn_eval,color=cgColor(colors[m]),line=1,/overplot
  endforeach
  
  legend,'z='+string(redshifts,format='(f3.1)'),textcolors=colors,box=0,/top,/right
  
  end_PS
  
  stop

end
