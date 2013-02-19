; timeScalesPlot.pro
; plotting: cooling times of halo gas vs. dynamical/hubble timescales
; dnelson feb.2013

; oplot2DHistoSq(): plot a 2d histogram in the separated squares style
    
pro oplot2DHistoSq, ct2d, hsp=hsp, nc=nc, xRange=xRange, yRange=yRange, logMass=logMass, logY=logY

  if ~keyword_set(hsp) or ~keyword_set(nc) then message,'error'

  ; store current table and load green-white for 2d histo
  tvlct, rr, gg, bb, /get
  
  loadColorTable,'green-white linear', /reverse
    
  ; process data
  w = where(ct2d.h2 eq 0.0,count,comp=wc)
    
  if keyword_set(logMass) then begin
    ct2d.h2 = ct2d.h2*1e10
    if count gt 0 then ct2d.h2[w] = min(ct2d.h2,/nan)
  endif
    
  fieldMinMax = [min(ct2d.h2[wc],/nan),max(ct2d.h2,/nan)]
    
  for i=0,ct2d.nXBins-1 do begin
    x = [ct2d.binCenX[i]-ct2d.binSizeX/2+hsp[0], ct2d.binCenX[i]+ct2d.binSizeX/2-hsp[0], $ ; ll, lr
         ct2d.binCenX[i]+ct2d.binSizeX/2-hsp[0], ct2d.binCenX[i]-ct2d.binSizeX/2+hsp[0]]   ; ur, ul
           
    if min(x) lt xRange[0] or max(x) gt xRange[1] then continue
             
    for j=0,ct2d.nYBins-1 do begin
      y = [ct2d.binCenY[j]-ct2d.binSizeY/2+hsp[1], ct2d.binCenY[j]-ct2d.binSizeY/2+hsp[1], $ ; ll, lr
           ct2d.binCenY[j]+ct2d.binSizeY/2-hsp[1], ct2d.binCenY[j]+ct2d.binSizeY/2-hsp[1]]   ; ur, ul
             
      if keyword_set(logY) then y = 10.0^y
             
      if min(y) lt yRange[0] or max(y) gt yRange[1] then continue
             
      ; determine color and make polygon
      colorind = (ct2d.h2[i,j]-fieldMinMax[0])*nc / (fieldMinMax[1]-fieldMinMax[0]) ;0-nc
      colorind = fix(colorind + 0.0) > 0 < 255 ;0-nc
        
      cgPolygon,x,y,color=colorind,/fill
    endfor
  endfor
  
  ; restore original CT
  tvlct, rr, gg, bb
  
end

; plotTSFracsVsHaloMass(): plot gas mass fractions for timescale ratios vs halo mass

pro plotTSFracsVsHaloMass, res=res

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  if ~keyword_set(res) then message,'error: set res'
  
  ; config
  redshift = 2.0
  sgSelect = 'pri'
  
  sP = simParams(res=res,run='tracer',redshift=redshift)
  
  xrange = [10.0,12.0]
  yrange = [0.0,1.03]
  lines  = [1,0,2]
  cInd   = 1 ; colorindex
  
  ; 2d config
  binSizeMass = 0.07 / (sP.res/128)
  binSizeFrac = 0.05 / (sP.res/128)
  
  hsp = [0.003,0.0035] ; mass, frac
  nc  = 150 ; number of colors (of 255) to use for background 2d histo

  ; load timescales
  ts = timescaleFracsVsHaloMass(sP=sP,sgSelect=sgSelect)
  
  ; calculate mass fraction with tcool<tdyn for SIS
  mmf = modelMassFracs(sP=sP, ts=ts)

  ; plot (1) - tcool<>tdyn
  start_PS, sP.plotPath + 'tsfrac_tdyn_vsmass.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),ytitle="Gas Mass Fraction "+textoidl("(t_{cool} < t_{dyn})")
      
    ; sim 2d histo
    yy = ts.tsFracs.gmem_tcool_tdyn[1,*]
    
    f2d = binHisto2D(xx=ts.gcMasses, yy=yy, xmm=xrange, ymm=yrange, xbs=binSizeMass, ybs=binSizeFrac)
             
    oplot2DHistoSq, f2d, hsp=hsp, nc=nc, xrange=xrange, yrange=yrange
             
    ; sim lines
    j=1
    cgPlot,ts.gcMasses,ts.tsFracs.gmem_tcool_tdyn[j,*],psym=4,color=sP.colorsA[cInd],/overplot
      
    for j=0,n_elements(ts.tsRatioVals)-1 do $
      cgPlot,ts.logMassBinCen,ts.tsMedian.gmem_tcool_tdyn[j,*],color=cgColor('black'),line=lines[j],/overplot
      
    ; models
    mMasses = codeMassToLogMsun(mmf.modelMasses)
    sK = 5
    for j=0,n_elements(ts.tsRatioVals)-1 do begin
      cgPlot,mMasses,smooth(mmf.modelMassFrac_SIS2[j,*],sK),line=lines[j],color=cgColor('orange'),/overplot
      cgPlot,mMasses,smooth(mmf.modelMassFrac_NFW_iso[j,*],sK),line=lines[j],color=cgColor('blue'),/overplot
      cgPlot,mMasses,smooth(mmf.modelMassFrac_NFW_poly[j,*],sK),line=lines[j],color=cgColor('red'),/overplot
    endfor
      
    ; legend
    strings = textoidl("t_{cool} / t_{dyn} < "+string(ts.tsRatioVals,format='(f3.1)'))
    legend,strings,linestyle=lines,box=0,linesize=0.25,/top,/left,charsize=!p.charsize-0.1,$
      textcolors=['black'],color=cgColor('black')
    
    strings = ['sim','SIS','NFW iso','NFW poly']
    legend,strings,textcolors=['black','orange','blue','red'],box=0,/top,/right,charsize=!p.charsize-0.1
    
    ; redo plot borders
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xtitle="",ytitle="",/noerase
    
  end_PS
  
  ; plot (2) - tcool<>tage
  start_PS, sP.plotPath + 'tsfrac_tage_vsmass.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),ytitle="Gas Mass Fraction"
      
    j=1
    cgPlot,ts.gcMasses,ts.tsFracs.gmem_tcool_tage[j,*],psym=4,color=sP.colorsA[cInd],/overplot
      
    for j=0,n_elements(ts.tsRatioVals)-1 do $
      cgPlot,ts.logMassBinCen,ts.tsMedian.gmem_tcool_tage[j,*],color=cgColor('black'),line=j,/overplot
      
    ; legend
    strings = textoidl("t_{cool} / t_{age} < "+string(ts.tsRatioVals,format='(f3.1)'))
    legend,strings,linestyle=indgen(n_elements(ts.tsRatioVals)),box=0,linesize=0.25,$
      position=[11.25,0.97],charsize=!p.charsize-0.1
  end_PS
  
  ; plot (3) - do hot halo masses (measured) agree with omega_b * M_DM
  start_PS, sP.plotPath + 'tsfrac_hotmass.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=[0.001,3.0],/ylog,/xs,/ys,$
      xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),ytitle=textoidl("M_{hot,sim} / (f_b M_{DM})")
      
    cgPlot,xrange,[1.0,1.0],line=2,/overplot
      
    cgPlot,ts.gcMasses,ts.gmem_hotmasses/(units.f_b*ts.gcMasses),psym=4,color=cgColor('blue'),/overplot
    
    ; recalculated gcMasses for DM only
    gc = loadGroupCat(sP=sP,/skipIDs)
    gcIDList = gcIDList(gc=gc,select=sgSelect)
    gcMasses_DM = codeMassToLogMsun(gc.subgroupMassType[partTypeNum('dm'),gcIDList])
    
    cgPlot,ts.gcMasses,ts.gmem_hotmasses/(units.f_b*gcMasses_DM),psym=4,color=cgColor('red'),/overplot
  end_PS
  
  stop

end

; timescaleRadStack(): 2D bin cooling time and overplot other timescales for 4 stacked halo mass bins
pro timescaleRadStack

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  sP = simParams(res=128,run='tracer',redshift=2.0)
  gc = loadGroupCat(sP=sP,/skipIDs)
  h  = loadSnapshotHeader(sP=sP)
  
  massRanges = list([10.5,10.7], [10.9,11.1], [11.3,11.5], [11.7,11.9])
  
  ; 2d histo config
  binSizeRad  = 0.045 / (sP.res/128) > 0.0225
  binSizeTime = 0.15 / (sP.res/128) > 0.075
  
  mmRad  = [0.0,1.5] ; r/rvir
  mmTime = [-3.0,2.0]  ; log(Gyr)
  
  hsp = [0.004,0.025] ; r/rvir, Gyr
  nc  = 150 ; number of colors (of 255) to use for background 2d histo
  
  ; plot config - 2x2 layout
  radRange = [0.05,1.6]  ; r/rvir
  tsRange  = [0.004,40.0] ; Gyr
  
  x0 = 0.12 & x1 = 0.50 & x2 = 0.88
  y0 = 0.14 & y1 = 0.54 & y2 = 0.94
    
  pos = list( [x0,y1,x1,y2] ,$ ; upper left
              [x1,y1,x2,y2] ,$ ; upper right
              [x0,y0,x1,y1] ,$ ; lower left
              [x1,y0,x2,y1]  ) ; lower right
  
  ; start plot
  start_PS, sP.plotPath + 'timescales_2x2.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
      
  foreach massRange,massRanges,k do begin
    ; make list of halos
    gcIDList = where(codeMassToLogMsun(gc.subgroupMass) ge massRange[0] and $
                    codeMassToLogMsun(gc.subgroupMass) lt massRange[1],count)
                    
    print,k,count
    
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
    
    oplot2DHistoSq, ct2d, hsp=hsp, nc=nc, xrange=radRange, yrange=tsRange, /logMass, /logY
    
    ; radial fits: cooling and dynamical
    cgPlot,ts.radCt.binCen,ts.radCt.radMedian,line=0,color=cgColor('black'),/overplot
    cgPlot,ts.radDt.binCen,ts.radDt.radMedian,line=2,color=cgColor('black'),/overplot
                   
    ; SIS/NFW models: cooling and dynamical
    cgPlot,[0.15,1.5],[ts.sis_gas.dynTime,ts.sis_gas.dynTime],line=2,color=cgColor('slate blue'),/overplot
    cgPlot,ts.x,ts.sis_gas.coolTime,line=0,color=cgColor('slate blue'),/overplot
    
    cgPlot,ts.x,ts.nfw_gas.dynTime,line=2,color=cgColor('orange'),/overplot
    cgPlot,ts.x,ts.nfw_gas.coolTime,line=0,color=cgColor('orange'),/overplot
           
    cgText,radRange[1]*0.7,tsRange[0]*2,"M = "+string(mean(massRange),format='(f4.1)'),/data
           
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
    legend,textoidl(['\tau_{dyn}','\tau_{cool}']),linestyle=[2,0],box=0,pos=[radRange[0],tsRange[1]*0.7],linesize=0.25
           
  endforeach
  
  ; colorbar
  cbTickNames = linspace(0.0,fieldMinMax[1],7)
  cbTickNames = ['0',string(alog10(cbTickNames[1:*]),format='(f4.1)')]
  
  loadColorTable,'green-white linear', /reverse
  cgColorbar, range=fieldMinMax, /vertical, ncolors=nc, pos=[x2+0.01,y0,x2+0.04,y2], $
    ticknames=cbTickNames, /right
  
  ; labels
  cgText,mean([x0,x2]),0.03,textoidl("r / r_{vir}"),/normal,alignment=0.5
  cgText,0.03,mean([y0,y2]),textoidl("Timescale [Gyr]"),/normal,alignment=0.5,orientation=90.0
  cgText,0.97,mean([y0,y2]),textoidl("Total Gas Mass [_{ } log M_{sun }]"),/normal,alignment=0.5,orientation=-90.0
  
  end_PS
  
  stop

end

; compareTimescalesHalo(): explore for a single halo or a mass range stacked

pro compareTimescalesHalo, redshift=redshift

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  if ~keyword_set(redshift) then message,'set redshift'
  sP = simParams(res=128,run='tracer',redshift=redshift)
  gc = loadGroupCat(sP=sP,/skipIDs)

  ; single halo by haloID
  haloID = 130 ;z2.304 z2.301 z2.130 z2.64
  gcID = getMatchedIDs(sPa=sP,sPg=sP,haloID=haloID)
  gcIDList = [gcID.a]
  hTag = 'h'+str(haloID)+'_m='+string(codeMassToLogMsun(gc.subgroupMass[gcID.a]),format='(f4.1)')

  ; load gas timescales and best model fits
  ts = loadFitTimescales(sP=sP,gcIDList=gcIDList)
  
  ; plot config
  coolingRange = [0.05,20]
  tempRange = [4.0,7.0]
  densRange = 10.0^[-9.0,-5.5]
  radRange = [0.1,1.75]
  binsize = 0.1 / (sP.res/128)
  
  psym = 4
  if n_elements(ts.coolTime) gt 10000 then psym = 3
  
  ; plot (1) - histograms
  start_PS, sP.plotPath + 'timescales_histo.' + hTag + '.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps'
    
    w = where(ts.coolTime eq 0.0,count)
    logCt = alog10(ts.coolTime)
    logCt[w] = 0.0
    hist = histogram(logCt,binsize=binsize,loc=loc)
    
    cgPlot,10.0^loc,hist/float(total(hist)),color=cgColor('red'),/xs,$
      xtitle="Timescale [Gyr]",ytitle="Halo Gas Mass Fraction",xrange=coolingRange,yrange=[0.0,0.7]/(sP.res/128),/xlog
    
    hist = histogram(alog10(ts.dynTime[where(finite(ts.dynTime))]),binsize=binsize,loc=loc)
    cgPlot,10.0^loc,hist/float(total(hist)),color=cgColor('blue'),/overplot
    
    cgPlot,[ts.hubbleTime,ts.hubbleTime],[0,1],color=cgColor('green'),/overplot
    
    legend,textoidl(['t_{cool}','t_{dyn}','t_{age}']),textcolors=['red','blue','green'],/top,/right,box=0
    
  end_PS
  
  ; plot (2) - cooling/dynamical scatter
  start_PS, sP.plotPath + 'timescales_scat.' + hTag + '.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps'
    cgPlot,coolTime,dynTime,psym=psym,xtitle="Cooling Time [Gyr]",ytitle="Dynamical Time [Gyr]",xrange=coolingRange
  end_PS

  ; plot (3) - cooling time vs current temperature scatter
  start_PS, sP.plotPath + 'timescales_cooltemp_scat.'+hTag+'.'+str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)+'.eps'
    cgPlot,ts.coolTime,ts.curTemp,psym=psym,xtitle="Cooling Time [Gyr]",ytitle="Gas Temperature [log K]",$
      xrange=coolingRange,yrange=tempRange
  end_PS
  
  ; plot (4) - cooling time vs current density scatter
  start_PS, sP.plotPath + 'timescales_cooldens_scat.' + hTag + '.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps'
    cgPlot,ts.coolTime,ts.curDens,psym=psym,xtitle="Cooling Time [Gyr]",ytitle="Density [Code]",$
      xrange=coolingRange,yrange=densRange,/ylog
  end_PS
  
  ; plot (5) - cooling/dynamical vs radius (temp/dens colormap)
  start_PS, sP.plotPath + 'timescales_vsrad1.' + hTag + '.' + str(sP.res) + '.' + sP.plotPrefix+'.'+str(sP.snap)+'.eps'
    cgPlot,[0],[0],/nodata,xtitle=textoidl("r / r_{vir}"),ytitle="Timescale [Gyr]",$
      xrange=radRange,yrange=coolingRange,/ylog,/xs,/ys,yminor=0,$
      title=hTag + ' (z='+string(sP.redshift,format='(f3.1)')+')',position=[0.15,0.13,0.85,0.9]
    
    ; individual gas elements: cooling time color mapped by temperature
    loadColorTable,'blue-red';,/reverse
    TVLCT, rr, gg, bb, /GET
    
    ; temp colormap:
    fieldMinMax = [5.0,6.5]
    colorinds = (ts.curTemp-fieldMinMax[0])*235.0 / (fieldMinMax[1]-fieldMinMax[0]) ;0-235
    
    ; dens colormap:
    ;fieldMinMax = [-8.0,-6.0]
    ;colorinds = (alog10(curDens)-fieldMinMax[0])*235.0 / (fieldMinMax[1]-fieldMinMax[0]) ;0-235
    
    ; individual gas elements
    colorinds = fix(colorinds + 20.0) > 0 < 255 ;20-255
    colorinds = getColor24([[rr[colorinds]], [gg[colorinds]], [bb[colorinds]]])
    
    for i=0L,n_elements(ts.coolTime)-1 do $
      oplot,[ts.gasRadii[i]/ts.gasRvir[i]],[ts.coolTime[i]],psym=psym,color=colorinds[i]
    
    ;cgPlot,gasRadii/gasRvir,coolTime,psym=psym,color=cgColor('red'),/overplot
    
    ; dynTime(r) and hubbleTime
    cgPlot,ts.gasRadii/ts.gasRvir,ts.dynTime,psym=psym,color=cgColor('blue'),/overplot
    cgPlot,[0.05,1.5],[ts.hubbleTime,ts.hubbleTime],color=cgColor('green'),/overplot
    
    ; mean halo
    cgPlot,[0.05,1.0],[ts.coolTime_halo,ts.coolTime_halo],color=cgColor('magenta'),line=0,/overplot
    cgPlot,[0.05,1.0],[ts.dynTime_halo,ts.dynTime_halo],color=cgColor('orange'),line=0,/overplot
     
    ; radial fits
    cgPlot,ts.radCt.binCen,ts.radCt.radMedian,line=0,color=cgColor('magenta'),/overplot
    cgPlot,ts.radDt.binCen,ts.radDt.radMedian,line=0,color=cgColor('orange'),/overplot
     
    ; SIS/NFW models
    cgPlot,ts.x,ts.sis_gas.coolTime,line=0,color=cgColor('brown'),/overplot
    cgPlot,[0.15,1.5],[ts.sis_gas.dynTime,ts.sis_gas.dynTime],line=0,color=cgColor('pink'),/overplot
    
    cgPlot,ts.x,ts.nfw_gas.dynTime,line=0,color=cgColor('dark gray'),/overplot
    cgPlot,ts.x,ts.nfw_gas.coolTime,line=0,color=cgColor('forest green'),/overplot
     
    ; mark various r_cool
    cgPlot,[ts.sis_gas.r_cool_h/ts.sis_dm.r200,ts.sis_gas.r_cool_h/ts.sis_dm.r200],$
      [0.06,0.1],line=2,color=cgColor('brown'),/overplot
    cgPlot,[ts.sis_gas.r_cool/ts.sis_dm.r200,ts.sis_gas.r_cool/ts.sis_dm.r200],$
      [0.06,0.1],line=0,color=cgColor('brown'),/overplot
     
    cgPlot,[ts.nfw_gas.r_cool/ts.nfw_dm.r200,ts.nfw_gas.r_cool/ts.nfw_dm.r200],$
      [0.06,0.1],line=0,color=cgColor('dark gray'), /overplot
     
    ; legend
    legend,textoidl(['t_{cool}','t_{dyn}','t_{age}','t_{cool,halo}','t_{dyn,halo}',$
                     't_{cool,sis}','t_{dyn,sis}','t_{cool,nfw}','t_{dyn,nfw}']),$
      textcolors=['red','blue','green','magenta','orange','brown','pink','forest green','dark gray'],$
      /bottom,/right,box=0
      
    ; colorbar
    cgColorbar,/right,/vertical,position=[0.88,0.13,0.93,0.9],range=fieldMinMax,title="log Temp"
  end_PS
  
  ; plot (6) - gas density/temperature vs radius
  start_PS,sP.plotPath+'timescales_denstemp_rad.'+hTag+'.'+str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)+'.eps',xs=6,ys=8
    ; temp
    cgPlot,[0],[0],/nodata,xtitle="",ytitle="log(T) [K]",$
      xrange=radRange,yrange=tempRange+[0.1,0.0],/xs,/ys,/xlog,xminor=0,$
      position=[0.15,0.55,0.95,0.95],xtickname=replicate(' ',10)
      
    ; individual gas elements and median binned
    cgPlot,ts.gasRadii/ts.gasRvir,ts.curTemp,psym=psym,color=cgColor('blue'),/overplot
    cgPlot,ts.radTemp.binCen,ts.radTemp.radMedian,psym=-4,color=cgColor('black'),/overplot
    
    ; SIS model and SIS fit
    cgPlot,ts.x,alog10(ts.sis_gas.temp_gas),line=0,color=cgColor('forest green'),/overplot
    cgPlot,ts.x,ts.sis_fit.f_temp,line=2,color=cgColor('green'),/overplot
    
    ; NFW iso
    cgPlot,ts.x,alog10(replicate(ts.nfw_gas.T_0,n_elements(ts.x))),line=0,color=cgColor('red'),/overplot
    
    ; NFW model and NFW fit (poly)
    cgPlot,ts.x,alog10(ts.nfw_gas.temp_gas),line=0,color=cgColor('orange'),/overplot    
    cgPlot,ts.x,ts.nfw_fit.f_temp,line=2,color=cgColor('orange'),/overplot
    
    ; dens
    cgPlot,[0],[0],/nodata,xtitle=textoidl("r / r_{vir}"),ytitle="Density [Code]",$
      xrange=radRange,yrange=densRange,/ylog,/xs,/ys,/xlog,yminor=0,xminor=0,$
      /noerase,position = [0.15,0.15,0.95,0.55]
      
    ; individual gas elements and median binned
    cgPlot,ts.gasRadii/ts.gasRvir,ts.curDens,psym=psym,color=cgColor('blue'),/overplot
    cgPlot,ts.radDens.binCen,ts.radDens.radMedian,psym=-4,color=cgColor('black'),/overplot
    
    ; SIS model and SIS fit
    cgPlot,ts.x,ts.sis_gas.rho_gas,line=0,color=cgColor('forest green'),/overplot
    cgPlot,ts.x,ts.sis_fit.f_dens,line=2,color=cgColor('forest green'),/overplot
    
    ; NFW iso dens
    cgPlot,ts.x,ts.nfw_gas.rho_gas_iso,line=0,color=cgColor('red'),/overplot
    
    ; NFW model and NFW fit (poly)
    cgPlot,ts.x,ts.nfw_gas.rho_gas,line=0,color=cgColor('orange'),/overplot
    cgPlot,ts.x,ts.nfw_fit.f_dens,line=2,color=cgColor('orange'),/overplot
    
    legend,['median','SIS','NFW iso','NFW poly'],textcolors=['black','forest green','red','orange'],$
      box=0,/bottom,/left
    legend,['model','fit'],linestyle=[0,2],linesize=0.25,box=0,/top,/right
    
  end_PS
  
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

  match,galcat.groupmemIDs,ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
  ids_gmem_ind = ids_gmem_ind[sort(galcat_ind)]
  
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
