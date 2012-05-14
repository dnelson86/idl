; plotVsHaloMass.pro
; gas accretion project - plots as a function of halo mass
; dnelson may.2012

; plotAngMomRadVsHaloMass(): plot the median angular momentum of hot/cold modes at different radial
;                            times as a function of halo mass

pro plotAngMomRadVsHaloMass

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect = 'pri'
  accMode = 'all'
  
  res = 128
  sPg = simParams(res=res,run='gadget',redshift=2.0)
  sPa = simParams(res=res,run='tracer',redshift=2.0)
  
  binnedGadget = haloMassBinAngMom(sP=sPg,sgSelect=sgSelect,accMode=accMode)
  binnedArepo  = haloMassBinAngMom(sP=sPa,sgSelect=sgSelect,accMode=accMode)

  ; plot (1) - lines for all rVirFacs vs halo mass
  start_PS, sPg.plotPath + 'angmom.vshalo.comp.'+accMode+'.'+str(sPg.res)+'_'+str(sPg.snap)+'.eps', /big
    
    xrange = [9.5,12.5]
    yrange = [0.3,1.0]
    
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
      
      w = where(binnedGadget.logMassBinCen gt 8.0)
      
      ; gadget hot both
      cgPlot,binnedGadget.logMassBinCen[w],binnedGadget.hotMode.median_both[j,w],color=getColor(1),line=1,/overplot
      cgPlot,binnedGadget.logMassBinCen[w],binnedGadget.coldMode.median_both[j,w],color=getColor(1),line=2,/overplot

      ; arepo hot both
      cgPlot,binnedArepo.logMassBinCen[w],binnedArepo.hotMode.median_both[j,w],color=getColor(3),line=1,/overplot
      cgPlot,binnedArepo.logMassBinCen[w],binnedArepo.coldMode.median_both[j,w],color=getColor(3),line=2,/overplot
      
      ; legend
      massBinStr = textoidl('r/r_{vir}')+' = ' + string(binnedGadget.rVirFacs[j],format='(f4.2)')
      cgText,mean(xrange),yrange[1]*0.9,massBinStr,charsize=!p.charsize-0.2,alignment=0.5,$
        color=cgColor('forest green')
    
    endfor
    
    ; axis labels
    cgText,0.05,0.5,textoidl("< j / j_{circ} >"),alignment=0.5,orientation=90.0,/normal
    cgText,0.5,0.05,textoidl("log ( M_{halo} ) [_{ }M_{sun }]"),alignment=0.5,/normal
    
    ; annotations
    cgText,0.95,0.06,"gadget",alignment=0.5,charsize=!p.charsize-0.2,color=getColor(1),/normal
    cgText,0.95,0.02,"arepo",alignment=0.5,charsize=!p.charsize-0.2,color=getColor(3),/normal
    legend,['cold','hot'],linestyle=[1,2],box=0,linesize=0.25,/bottom,/left
    
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
      
      w = where(binnedGadget.logMassBinCen gt 8.0)
      
      ; gadget hot both
      cgPlot,binnedGadget.logMassBinCen[w],binnedGadget.hotMode.median_both[j,w],color=getColor(1),line=1,/overplot
      cgPlot,binnedGadget.logMassBinCen[w],binnedGadget.coldMode.median_both[j,w],color=getColor(1),line=2,/overplot

      ; arepo hot both
      cgPlot,binnedArepo.logMassBinCen[w],binnedArepo.hotMode.median_both[j,w],color=getColor(3),line=1,/overplot
      cgPlot,binnedArepo.logMassBinCen[w],binnedArepo.coldMode.median_both[j,w],color=getColor(3),line=2,/overplot
      
      ; legend
      massBinStr = textoidl('r/r_{vir}')+' = ' + string(binnedGadget.rVirFacs[j],format='(f4.2)')
      cgText,mean(xrange),yrange[1]*0.9,massBinStr,charsize=!p.charsize-0.2,alignment=0.5,$
        color=cgColor('forest green')
    
    endfor
    
    ; axis labels
    cgText,0.05,0.5,textoidl("\Delta t_{acc} / \tau_{circ}"),alignment=0.5,orientation=90.0,/normal
    cgText,0.5,0.05,textoidl("log ( M_{halo} ) [_{ }M_{sun }]"),alignment=0.5,/normal
    
    ; annotations
    cgText,0.95,0.06,"gadget",alignment=0.5,charsize=!p.charsize-0.2,color=getColor(1),/normal
    cgText,0.95,0.02,"arepo",alignment=0.5,charsize=!p.charsize-0.2,color=getColor(3),/normal
    legend,['cold','hot'],linestyle=[1,2],box=0,linesize=0.25,/bottom,/left
    
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

; plotColdFracVsHaloMass(): plot the "cold mode fraction" vs halo mass in a few different ways
;                           using the mergerTreeSubset with accretionTimes and accretionModes

pro plotColdFracVsHaloMass

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sgSelect = 'pri'
  accMode = 'all' ; accretion mode: all, smooth, bclumpy, sclumpy
  
  res = 256
  sPg = simParams(res=res,run='gadget',redshift=2.0)
  sPa = simParams(res=res,run='tracer',redshift=2.0)
  
  cfG = haloMassBinColdFracs(sP=sPg,sgSelect=sgSelect,accMode=accMode)
  cfA = haloMassBinColdFracs(sP=sPa,sgSelect=sgSelect,accMode=accMode)
  
  ; minimum halo mass to plot
  w = where(cfG.logMassBinCen gt 10.0)  
  
  ; plot (1) - median lines (Tmax) 1*tvircur,1*tviracc,tcut=5.5
  start_PS, sPg.plotPath + 'coldFrac1.'+accMode+'.comp.'+str(sPg.res)+'_'+str(sPg.snap)+'.eps'
    
    xrange = [10.0,12.5]
    yrange = [0.0,1.1]
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("log ( M_{halo} ) [_{ }M_{sun }]")
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; gadget gal
    cgPlot,cfG.logMassBinCen[w],cfG.medianVals.tVircur_gal[1,w],color=getColor(1),line=1,/overplot
    cgPlot,cfG.logMassBinCen[w],cfG.medianVals.tViracc_gal[1,w],color=getColor(1),line=0,/overplot
    cgPlot,cfG.logMassBinCen[w],cfG.medianVals.const_gal[2,w],color=getColor(1),line=2,/overplot
    
    ; arepo gal
    cgPlot,cfA.logMassBinCen[w],cfA.medianVals.tVircur_gal[1,w],color=getColor(3),line=1,/overplot
    cgPlot,cfA.logMassBinCen[w],cfA.medianVals.tViracc_gal[1,w],color=getColor(3),line=0,/overplot
    cgPlot,cfA.logMassBinCen[w],cfA.medianVals.const_gal[2,w],color=getColor(3),line=2,/overplot
      
    ; gadget gmem
    ;cgPlot,cfG.logMassBinCen[w],cfG.medianVals.tVircur_gmem[1,w],color=getColor(1),line=1,/overplot
    ;cgPlot,cfG.logMassBinCen[w],cfG.medianVals.tViracc_gmem[1,w],color=getColor(1),line=0,/overplot
    ;cgPlot,cfG.logMassBinCen[w],cfG.medianVals.const_gmem[2,w],color=getColor(1),line=2,/overplot 
    
    ; arepo gmem
    ;cgPlot,cfA.logMassBinCen[w],cfA.medianVals.tVircur_gmem[1,w],color=getColor(3),line=1,/overplot
    ;cgPlot,cfA.logMassBinCen[w],cfA.medianVals.tViracc_gmem[1,w],color=getColor(3),line=0,/overplot
    ;cgPlot,cfA.logMassBinCen[w],cfA.medianVals.const_gmem[2,w],color=getColor(31),line=2,/overplot
    
    ; legend
    strings = [textoidl("T_{max} < T_{vir,cur}"),$
               textoidl("T_{max} < T_{vir,acc}"),$
               textoidl("T_{max} < T_{c} = 10^{5.5} K")]
    legend,strings,linestyle=[1,0,2],box=0,linesize=0.25,/bottom,/left
    cgText,0.05,0.06,"gadget",alignment=0.5,charsize=!p.charsize-0.2,color=getColor(1),/normal
    cgText,0.05,0.02,"arepo",alignment=0.5,charsize=!p.charsize-0.2,color=getColor(3),/normal
    
  end_PS  
  
  ; plot (2) - median lines (Tmax) variable*tviracc
  start_PS, sPg.plotPath + 'coldFrac2.'+accMode+'.comp.'+str(sPg.res)+'_'+str(sPg.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("log ( M_{halo} ) [_{ }M_{sun }]")
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; gadget and arepo gal
    for j=0,cfG.nVirs-1 do $
      cgPlot,cfG.logMassBinCen[w],cfG.medianVals.tViracc_gal[j,w],color=getColor(1),line=j,/overplot
    for j=0,cfA.nVirs-1 do $
      cgPlot,cfA.logMassBinCen[w],cfA.medianVals.tViracc_gal[j,w],color=getColor(3),line=j,/overplot
      
    ; gadget and arepo gmem
    ;for j=0,cfG.nVirs-1 do $
    ;  cgPlot,cfG.logMassBinCen[w],cfG.medianVals.tViracc_gmem[j,w],color=getColor(4),line=j,/overplot
    ;for j=0,cfG.nVirs-1 do $
    ;  cgPlot,cfA.logMassBinCen[w],cfA.medianVals.tViracc_gmem[j,w],color=getColor(5),line=j,/overplot
    
    ; legend
    strings = textoidl("T_{max} / T_{vir,acc} < ")+string(cfG.TvirVals,format='(f4.1)')
    legend,strings,linestyle=indgen(cfG.nVirs),box=0,linesize=0.25,position=[11.35,0.87]
    cgText,0.05,0.06,"gadget",alignment=0.5,charsize=!p.charsize-0.2,color=getColor(1),/normal
    cgText,0.05,0.02,"arepo",alignment=0.5,charsize=!p.charsize-0.2,color=getColor(3),/normal
    
  end_PS
  
  ; plot (3) - median lines (Tmax) variable*tconst
  start_PS, sPg.plotPath + 'coldFrac3.'+accMode+'.comp.'+str(sPg.res)+'_'+str(sPg.snap)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("log ( M_{halo} ) [_{ }M_{sun }]")
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; gadget and arepo gal
    for j=0,cfG.nCuts-1 do $
      cgPlot,cfG.logMassBinCen[w],cfG.medianVals.const_gal[j,w],color=getColor(1),line=j,/overplot
    for j=0,cfG.nCuts-1 do $
      cgPlot,cfA.logMassBinCen[w],cfA.medianVals.const_gal[j,w],color=getColor(3),line=j,/overplot
      
    ; gadget and arepo gmem
    ;for j=0,cfG.nCuts-1 do $
    ;  cgPlot,cfG.logMassBinCen[w],cfG.medianVals.const_gmem[j,w],color=getColor(4),line=j,/overplot
    ;for j=0,cfG.nCuts-1 do $
    ;  cgPlot,cfA.logMassBinCen[w],cfA.medianVals.const_gmem[j,w],color=getColor(5),line=j,/overplot
    
    ; legend
    strings = textoidl("T_{max} / T_{c} = ")+string(cfG.TcutVals,format='(f4.1)')
    legend,strings,linestyle=indgen(cfG.nCuts),box=0,linesize=0.25,/bottom,/left
    cgText,0.05,0.06,"gadget",alignment=0.5,charsize=!p.charsize-0.2,color=getColor(1),/normal
    cgText,0.05,0.02,"arepo",alignment=0.5,charsize=!p.charsize-0.2,color=getColor(3),/normal
    
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
  start_PS, sP.plotPath + 'coldFrac_noMTs.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
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
  start_PS, sP.plotPath + 'coldFrac_noMTs2.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
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