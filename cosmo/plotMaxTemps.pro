; plotMaxTemps.pro
; gas accretion project - plots related to maximum past temperature of the gas
; dnelson may.2012

; plotTempVsRad

pro plotTempVsRad

  sP = simParams(res=128,run='gadget',redshift=2.0)

  ; get gas/tracer positions with time
  at  = accretionTraj(sP=sP)
  mt  = mergerTreeSubset(sP=sP)

  hInd = 50
  haloIDs = [mt.galcatIDList[hInd]]

  ; get indices for the subset of the mergerTreeSubset corresponding to the selected halo(s)
  mergerTreeSub = mergerTreeINDList(sP=sP,gcIDList=haloIDs)

  ; count
  nPart = { gal  : n_elements(at.relPos_gal[0,0,mergerTreeSub.gal])   ,$
            gmem : n_elements(at.relPos_gmem[0,0,mergerTreeSub.gmem])  }

  rad = sqrt(at.relPos_gal[*,0,mergerTreeSub.gal] * at.relPos_gal[*,0,mergerTreeSub.gal] + $
             at.relPos_gal[*,1,mergerTreeSub.gal] * at.relPos_gal[*,1,mergerTreeSub.gal] + $
             at.relPos_gal[*,2,mergerTreeSub.gal] * at.relPos_gal[*,2,mergerTreeSub.gal])
  rad = reform(rad)
  
  temp = reform(at.curTemp_gal[*,mergerTreeSub.gal])
  
  ; normalize by rvir(t)
  for i=0,n_elements(mt.times)-1 do begin
    rad[i,*] /= (mt.hVirRad[i,hInd])[0]
  endfor

  ; plot (0)
  start_PS, sP.plotPath + 'tempvsrad.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    !p.thick -= 2
    xrange = [0,2]
    yrange = [3,7]
    
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="temp",xtitle="rad / rvir"      
      
      cgplot,reverse(mt.hvirrad[*,hInd]),reverse(mt.hvirtemp[*,hInd]),line=1,/overplot
      for j=0,50 do begin
        w = where(temp[*,j] ne 0.0,count,comp=wc,ncomp=ncomp)
        if count gt 0 then cgPlot,rad[w,j],temp[w,j],line=0,color=getColor(j),/overplot
        ;if ncomp gt 0 then cgPlot,rad[wc,j],temp[wc,j],line=1,color=getColor(j),/overplot
      endfor
    
  end_PS
  
  ; plot (0)
  start_PS, sP.plotPath + 'radvstime.'+str(sP.res)+'_'+str(sP.snap)+'.eps'

    xrange = [0,80]
    yrange = [0,2]
    
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="rad / rvir",xtitle="snaps back"      
      
      ;cgplot,reverse(mt.hvirrad[*,hInd]),line=1,/overplot
      for j=0,50 do $
        cgPlot,rad[*,j],line=0,color=getColor(j),/overplot
    
  end_PS
  
  stop
end

; plotTmaxVsTvirAccComp(); plot the previous max temp vs. the virial temperature of the parent halos at the
;                          time of accretion for arepo vs. gadget

pro plotTmaxVsTvirAccComp

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  res = 256
  redshift = 2.0

  sP1 = simParams(res=res,run='gadget',redshift=redshift)
  sP2 = simParams(res=res,run='tracer',redshift=redshift)
  
  binSizeLog = 0.15 / (res/128)
  
  sgSelect = 'pri'
  accMode  = 'smooth'
  massBins = [9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5] ; log(M)   
  
  ; load sP1 (gadget)
  accTvir_gadget = gcSubsetProp(sP=sP1,select=sgSelect,/accTvir,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  curTvir_gadget = gcSubsetProp(sP=sP1,select=sgSelect,/virTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  maxTemp_gadget = gcSubsetProp(sP=sP1,select=sgSelect,/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)  
 
  ; load parent halo masses so we can make halo massbins
  parentMass_ga = gcSubsetProp(sP=sP1,select=sgSelect,/parMass,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)

  ; load sP2 (tracer)
  accTvir_tracer = gcSubsetProp(sP=sP2,select=sgSelect,/accTvir,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  curTvir_tracer = gcSubsetProp(sP=sP2,select=sgSelect,/virTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  maxTemp_tracer = gcSubsetProp(sP=sP2,select=sgSelect,/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)

  parentMass_tr = gcSubsetProp(sP=sP2,select=sgSelect,/parMass,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
    
  ; plot (0) - 3x2 mass bins separated out and each panel with gadget+arepo, gal vs. gmem
  start_PS, sP1.plotPath + 'tmax_tviracc_3x2.'+accMode+'.'+str(sP1.res)+'_'+str(sP1.snap)+'.eps', /big
    !p.thick += 1
    xrange = [-2.2,1.2]
    yrange = [6e-4,1.0]
    
    xtickv = [-2.0,-1.0,0.0,1.0]
    
    x0 = 0.15 & x1 = 0.42 & x2 = 0.69 & x3 = 0.96
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; uc
                [x2,y1,x3,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1] ,$ ; lc
                [x2,y0,x3,y1] )  ; lr
    
    for j=1,n_elements(massBins)-2 do begin
      
      if j eq 1 or j eq 4 then ytickname = '' else ytickname = replicate(' ',10)
      if j ge 4 then xtickname = '' else xtickname = replicate(' ',10)
      if j gt 1 then noerase = 1 else noerase = 0
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,pos=pos[j-1],$
        ytitle="",xtitle="",xticks=3,xtickv=xtickv,xtickname=xtickname,ytickname=ytickname,noerase=noerase       
      
      cgPlot,[0,0],[8e-4,0.3],line=0,color=fsc_color('light gray'),thick=!p.thick-2.0,/overplot
      
      ; select members of this parent mass bins and r>0<inf
      wGadget_gal  = where(parentMass_ga.gal gt massBins[j] and parentMass_ga.gal le massBins[j+1],count1)
      wGadget_gmem = where(parentMass_ga.gmem gt massBins[j] and parentMass_ga.gmem le massBins[j+1],count2)
      wTracer_gal  = where(parentMass_tr.gal gt massBins[j] and parentMass_tr.gal le massBins[j+1],count3)
      wTracer_gmem = where(parentMass_tr.gmem gt massBins[j] and parentMass_tr.gmem le massBins[j+1],count4)
      
      print,j,count1,count2,count3,count4
      
      if ~count1 or ~count2 or ~count3 or ~count4 then continue ; no halos in this mass bin  
      
      ; histogram gadget (gal) differences
      vals = [10.0^maxTemp_gadget.gal[wGadget_gal]/10.0^accTvir_gadget.gal[wGadget_gal]]
      hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=0,color=getColor(1),/overplot
  
      ; histogram tracer (gal) differences
      vals = [10.0^maxTemp_tracer.gal[wTracer_gal]/10.0^accTvir_tracer.gal[wTracer_gal]]
      hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=0,color=getColor(3),/overplot
    
      ; histogram gadget (gmem) differences
      vals = [10.0^maxTemp_gadget.gmem[wGadget_gmem]/10.0^accTvir_gadget.gmem[wGadget_gmem]]
      hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=1,color=getColor(1),/overplot
  
      ; histogram tracer (gmem) differences
      vals = [10.0^maxTemp_tracer.gmem[wTracer_gmem]/10.0^accTvir_tracer.gmem[wTracer_gmem]]
      hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=1,color=getColor(3),/overplot
    
      ; legend
      massBinStr = string(massBins[j],format='(f4.1)') + ' < log(M) < ' + $
                   string(massBins[j+1],format='(f4.1)')

      cgText,mean(xrange),yrange[1]*0.4,massBinStr,charsize=!p.charsize-0.2,alignment=0.5
    
    endfor
    
    ; axis labels
    cgText,0.05,0.5,"Gas Mass Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,0.5,0.05,textoidl("log ( T_{max} / T_{vir,acc} )"),alignment=0.5,/normal
    
    ; annotations
    cgText,0.95,0.06,"gadget",alignment=0.5,charsize=!p.charsize-0.2,color=getColor(1),/normal
    cgText,0.95,0.02,"arepo",alignment=0.5,charsize=!p.charsize-0.2,color=getColor(3),/normal
    
    legend,['central galaxy','halo atmosphere'],linestyle=[0,1],$
      color=cgColor('forest green'),textcolors=['forest green','forest green'],$
      linesize=0.25,box=0,position=[0.0,0.09],/normal,charsize=!p.charsize-0.2
      
  end_PS
stop
  ; plot (1) - compare tmax to tvir at time of accretion (both gal+gmem)
  start_PS, sP1.plotPath + 'tmax_tviracc_comp_both.'+str(sP1.res)+'_'+str(sP1.snap)+'.eps'
    !p.thick += 1
    xrange = [-2.5,1.2]
    yrange = [1e-3,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="Mass Fraction",xtitle=textoidl("log( T_{max} / T_{vir,acc} )")
      ;title=str(res)+textoidl("^3")+" Gadget vs. ArepoMC (z="+string(sP1.redshift,format='(f3.1)')+")"
    cgPlot,[0,0],yrange,line=0,color=fsc_color('light gray'),/overplot
    
    strings = []
    colors = []
    
    for j=0,n_elements(massBins)-2 do begin
      
      ; select members of this parent mass bins and r>0<inf
      wGadget_gal  = where(parentMass_ga.gal gt massBins[j] and parentMass_ga.gal le massBins[j+1],count1)
      wGadget_gmem = where(parentMass_ga.gmem gt massBins[j] and parentMass_ga.gmem le massBins[j+1],count2)
      wTracer_gal  = where(parentMass_tr.gal gt massBins[j] and parentMass_tr.gal le massBins[j+1],count3)
      wTracer_gmem = where(parentMass_tr.gmem gt massBins[j] and parentMass_tr.gmem le massBins[j+1],count4)
      
      print,j,count1,count2,count3,count4
      
      if ~count1 or ~count2 or ~count3 or ~count4 then continue ; no halos in this mass bin  
    
      massBinStr = string(massBins[j],format='(f4.1)') + '-' + string(massBins[j+1],format='(f4.1)')
      strings = [strings,massBinStr]
      colors = [colors,getColor(j,/name)]
      
      ; histogram gadget (gal+gmem) differences
      vals = [10.0^maxTemp_gadget.gal[wGadget_gal]/10.0^accTvir_gadget.gal[wGadget_gal],$
              10.0^maxTemp_gadget.gmem[wGadget_gmem]/10.0^accTvir_gadget.gmem[wGadget_gmem]]
      hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=0,color=getColor(j),/overplot
  
      ; histogram tracer (gal+gmem) differences
      vals = [10.0^maxTemp_tracer.gal[wTracer_gal]/10.0^accTvir_tracer.gal[wTracer_gal],$
              10.0^maxTemp_tracer.gmem[wTracer_gmem]/10.0^accTvir_tracer.gmem[wTracer_gmem]]
      hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=2,color=getColor(j),/overplot
    
    endfor
    
    ; legend
    legend,['gadget','arepo'],linestyle=[0,2],linesize=0.25,box=0,/right,/top
    if n_elements(massBins) gt 2 then $
      legend,strings,textcolors=colors,box=0,/left,/top
  
  end_PS
  
  ; plot (2) - compare tmax to tvir at time of accretion (gal only)
  start_PS, sP1.plotPath + 'tmax_tviracc_comp_gal.'+str(sP1.res)+'_'+str(sP1.snap)+'.eps'
    !p.thick += 1
    xrange = [-2.5,1.2]
    yrange = [1e-3,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="Mass Fraction",xtitle=textoidl("log( T_{max} / T_{vir,acc} )")
      ;title=str(res)+textoidl("^3")+" Gadget vs. ArepoMC (z="+string(sP1.redshift,format='(f3.1)')+")"
    cgPlot,[0,0],yrange,line=0,color=fsc_color('light gray'),/overplot
    
    strings = []
    colors = []
    
    for j=0,n_elements(massBins)-2 do begin
      
      ; select members of this parent mass bins and r>0<inf
      wGadget_gal  = where(parentMass_ga.gal gt massBins[j] and parentMass_ga.gal le massBins[j+1],count1)
      wTracer_gal  = where(parentMass_tr.gal gt massBins[j] and parentMass_tr.gal le massBins[j+1],count3)
      
      print,j,count1,count3
      
      if ~count1 or ~count3 then continue ; no halos in this mass bin  
    
      massBinStr = string(massBins[j],format='(f4.1)') + '-' + string(massBins[j+1],format='(f4.1)')
      strings = [strings,massBinStr]
      colors = [colors,getColor(j,/name)]
      
      ; histogram gadget (gal) differences
      vals = [10.0^maxTemp_gadget.gal[wGadget_gal]/10.0^accTvir_gadget.gal[wGadget_gal]]
      hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=0,color=getColor(j),/overplot
  
      ; histogram tracer (gal+gmem) differences
      vals = [10.0^maxTemp_tracer.gal[wTracer_gal]/10.0^accTvir_tracer.gal[wTracer_gal]]
      hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=2,color=getColor(j),/overplot
    
    endfor
    
    ; legend
    legend,['gadget','arepo'],linestyle=[0,2],linesize=0.25,box=0,/right,/top
    if n_elements(massBins) gt 2 then $
      legend,strings,textcolors=colors,box=0,/left,/top
  
  end_PS
  
  ; plot (3) - compare tmax to tvir at time of accretion (gmem)
  start_PS, sP1.plotPath + 'tmax_tviracc_comp_gmem.'+str(sP1.res)+'_'+str(sP1.snap)+'.eps'
    !p.thick += 1
    xrange = [-2.5,1.2]
    yrange = [1e-3,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="Mass Fraction",xtitle=textoidl("log( T_{max} / T_{vir,acc} )")
      ;title=str(res)+textoidl("^3")+" Gadget vs. ArepoMC (z="+string(sP1.redshift,format='(f3.1)')+")"
    cgPlot,[0,0],yrange,line=0,color=fsc_color('light gray'),/overplot
    
    strings = []
    colors = []
    
    for j=0,n_elements(massBins)-2 do begin
      
      ; select members of this parent mass bins and r>0<inf
      wGadget_gmem = where(parentMass_ga.gmem gt massBins[j] and parentMass_ga.gmem le massBins[j+1],count2)
      wTracer_gmem = where(parentMass_tr.gmem gt massBins[j] and parentMass_tr.gmem le massBins[j+1],count4)
      
      print,j,count2,count4
      
      if ~count2 or ~count4 then continue ; no halos in this mass bin  
    
      massBinStr = string(massBins[j],format='(f4.1)') + '-' + string(massBins[j+1],format='(f4.1)')
      strings = [strings,massBinStr]
      colors = [colors,getColor(j,/name)]
      
      ; histogram gadget (gal+gmem) differences
      vals = [10.0^maxTemp_gadget.gmem[wGadget_gmem]/10.0^accTvir_gadget.gmem[wGadget_gmem]]
      hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=0,color=getColor(j),/overplot
  
      ; histogram tracer (gal+gmem) differences
      vals = [10.0^maxTemp_tracer.gmem[wTracer_gmem]/10.0^accTvir_tracer.gmem[wTracer_gmem]]
      hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=2,color=getColor(j),/overplot
    
    endfor
    
    ; legend
    legend,['gadget','arepo'],linestyle=[0,2],linesize=0.25,box=0,/right,/top
    if n_elements(massBins) gt 2 then $
      legend,strings,textcolors=colors,box=0,/left,/top
  
  end_PS
  
  ; plot (4) - same plot using current Tvir instead of Tvir at accretion time
  start_PS, sP1.plotPath + 'tmax_tvircur_comp.'+str(sP1.res)+'_'+str(sP1.snap)+'.eps'
    !p.thick += 1
    xrange = [-2.5,1.0]
    yrange = [1e-3,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="Mass Fraction",xtitle=textoidl("log ( T_{max} / T_{vir,cur} )")+"",$
      title=str(res)+textoidl("^3")+" Gadget vs. ArepoMC (z="+string(sP1.redshift,format='(f3.1)')+")"
    cgPlot,[1,1],yrange,line=0,color=fsc_color('light gray'),/overplot
    
    strings = []
    colors = []
    
    for j=0,n_elements(massBins)-2 do begin
      
      ; select members of this parent mass bins and r>0<inf
      wGadget_gal  = where(parentMass_ga.gal gt massBins[j] and parentMass_ga.gal le massBins[j+1],count1)
      wGadget_gmem = where(parentMass_ga.gmem gt massBins[j] and parentMass_ga.gmem le massBins[j+1],count2)
      wTracer_gal  = where(parentMass_tr.gal gt massBins[j] and parentMass_tr.gal le massBins[j+1],count3)
      wTracer_gmem = where(parentMass_tr.gmem gt massBins[j] and parentMass_tr.gmem le massBins[j+1],count4)
      
      print,j,count1,count2,count3,count4
      
      if ~count1 or ~count2 or ~count3 or ~count4 then continue ; no halos in this mass bin  
    
      massBinStr = string(massBins[j],format='(f4.1)') + '-' + string(massBins[j+1],format='(f4.1)')
      strings = [strings,massBinStr]
      colors = [colors,getColor(j,/name)]
      
      ; histogram gadget (gal+gmem) differences
      vals = [10.0^maxTemp_gadget.gal[wGadget_gal]/10.0^curTvir_gadget.gal[wGadget_gal],$
              10.0^maxTemp_gadget.gmem[wGadget_gmem]/10.0^curTvir_gadget.gmem[wGadget_gmem]]
      hist = histogram(alog10(vals),binsize=binSizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=0,color=getColor(j),/overplot
  
      ; histogram tracer (gal+gmem) differences
      vals = [10.0^maxTemp_tracer.gal[wTracer_gal]/10.0^curTvir_tracer.gal[wTracer_gal],$
              10.0^maxTemp_tracer.gmem[wTracer_gmem]/10.0^curTvir_tracer.gmem[wTracer_gmem]]
      hist = histogram(alog10(vals),binsize=binSizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=2,color=getColor(j),/overplot
    
    endfor
    
    ; legend
    legend,['gadget','arepoMC'],linestyle=[0,2],linesize=0.25,box=0,/right,/top
    if n_elements(massBins) gt 2 then $
      legend,strings,textcolors=colors,box=0,/left,/top
  
  end_PS
  
  ; plot (5) - unscaled tmax
  start_PS, sP1.plotPath + 'tmax_nonorm.'+str(sP1.res)+'_'+str(sP1.snap)+'.eps'
    !p.thick += 1
    xrange = [4.0,7.0]
    yrange = [1e-3,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="Mass Fraction",xtitle=textoidl("log (T_{max})")+"",$
      title=str(res)+textoidl("^3")+" Gadget vs. ArepoMC (z="+string(sP1.redshift,format='(f3.1)')+")"
    cgPlot,[1,1],yrange,line=0,color=fsc_color('light gray'),/overplot
    
    strings = []
    colors = []
    
    for j=0,n_elements(massBins)-2 do begin
      
      ; select members of this parent mass bins and r>0<inf
      wGadget_gal  = where(parentMass_ga.gal gt massBins[j] and parentMass_ga.gal le massBins[j+1],count1)
      wGadget_gmem = where(parentMass_ga.gmem gt massBins[j] and parentMass_ga.gmem le massBins[j+1],count2)
      wTracer_gal  = where(parentMass_tr.gal gt massBins[j] and parentMass_tr.gal le massBins[j+1],count3)
      wTracer_gmem = where(parentMass_tr.gmem gt massBins[j] and parentMass_tr.gmem le massBins[j+1],count4)
      
      print,j,count1,count2,count3,count4
      
      if ~count1 or ~count2 or ~count3 or ~count4 then continue ; no halos in this mass bin  
    
      massBinStr = string(massBins[j],format='(f4.1)') + '-' + string(massBins[j+1],format='(f4.1)')
      strings = [strings,massBinStr]
      colors = [colors,getColor(j,/name)]
      
      ; histogram gadget (gal+gmem) differences
      vals = [maxTemp_gadget.gal[wGadget_gal],maxTemp_gadget.gmem[wGadget_gmem]]
      hist = histogram(vals,binsize=binSizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=0,color=getColor(j),/overplot
  
      ; histogram tracer (gal+gmem) differences
      vals = [maxTemp_tracer.gal[wTracer_gal],maxTemp_tracer.gmem[wTracer_gmem]]
      hist = histogram(vals,binsize=binSizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=2,color=getColor(j),/overplot
    
    endfor
    
    ; legend
    legend,['gadget','arepoMC'],linestyle=[0,2],linesize=0.25,box=0,/right,/top
    if n_elements(massBins) gt 2 then $
      legend,strings,textcolors=colors,box=0,/left,/top
  
  end_PS
  
end

; plotTmaxVsTvirAccCur(): evaluate how the ratio of Tmax/Tvir changes when using either the current
;                         Tvir of the halo or the Tvir of the halo at the time of accretion

pro plotTmaxVsTvirAccCur, sP=sP
  
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  binSizeLog = 0.1 / (sP.res/128)
  
  sgSelect = 'pri'
  ;massBins = [0.0,100.0] ; no massbins
  massBins = [9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5] ; log(M)   
  
  ; load
  accTvir = gcSubsetProp(sP=sP,select=sgSelect,/accTvir,/mergerTreeSubset,/accretionTimeSubset)
  curTvir = gcSubsetProp(sP=sP,select=sgSelect,/virTemp,/mergerTreeSubset,/accretionTimeSubset)
  maxTemp = gcSubsetProp(sP=sP,select=sgSelect,/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset)  

  ; load parent halo masses so we can make halo massbins
  parentMass = gcSubsetProp(sP=sP,select=sgSelect,/parMass,/mergerTreeSubset,/accretionTimeSubset)  
  
  ; plot (1)
  start_PS, sP.plotPath + 'tmax_tvircur_tviracc.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    !p.thick += 1
    xrange = [-2.5,1.0]
    yrange = [1e-3,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="fraction",xtitle=textoidl("log ( T_{max} / T_{vir} )")+"",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+")"
    cgPlot,[0,0],yrange,line=0,color=fsc_color('light gray'),/overplot
    
  strings = []
  colors = []
  
  for j=0,n_elements(massBins)-2 do begin
    
    ; select members of this parent mass bins and r>0<inf
    w_gal  = where(parentMass.gal gt massBins[j] and parentMass.gal le massBins[j+1],count1)
    w_gmem = where(parentMass.gmem gt massBins[j] and parentMass.gmem le massBins[j+1],count2)
    
    print,j,count1,count2
    if ~count1 or ~count2 then continue ; no halos in this mass bin  
  
    massBinStr = string(massBins[j],format='(f4.1)') + '-' + string(massBins[j+1],format='(f4.1)')
    strings = [strings,massBinStr]
    colors = [colors,getColor(j,/name)]
    
    ; histogram gadget (gal+gmem) differences for current Tvir
    vals = [10.0^maxTemp.gal[w_gal]/10.0^accTvir.gal[w_gal],$
            10.0^maxTemp.gmem[w_gmem]/10.0^accTvir.gmem[w_gmem]]
    hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=0,color=getColor(j),/overplot

    ; histogram gadget (gal+gmem) differences for accretion Tvir
    vals = [10.0^maxTemp.gal[w_gal]/10.0^curTvir.gal[w_gal],$
            10.0^maxTemp.gmem[w_gmem]/10.0^curTvir.gmem[w_gmem]]
    hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=2,color=getColor(j),/overplot
  
  endfor
  
  ; legend
  legend,[textoidl("T_{vir,acc}"),textoidl("T_{vir,cur}")],linestyle=[0,2],linesize=0.25,box=0,/left,/top
  if n_elements(massBins) gt 2 then $
    legend,strings,textcolors=colors,box=0,/right,/top
  
  end_PS
  
  ; plot (2) - just tviracc but separate out gal and gmem
  start_PS, sP.plotPath + 'tmax_tviracc_gal_gmem.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    !p.thick += 1
    xrange = [-2.5,1.0]
    yrange = [1e-3,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="fraction",xtitle=textoidl("log ( T_{max} / T_{vir} )")+"",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+")"
    cgPlot,[0,0],yrange,line=0,color=fsc_color('light gray'),/overplot
    
  strings = []
  colors = []
  
  for j=0,n_elements(massBins)-2 do begin
    
    ; select members of this parent mass bins and r>0<inf
    w_gal  = where(parentMass.gal gt massBins[j] and parentMass.gal le massBins[j+1],count1)
    w_gmem = where(parentMass.gmem gt massBins[j] and parentMass.gmem le massBins[j+1],count2)
    
    print,j,count1,count2
    if ~count1 or ~count2 then continue ; no halos in this mass bin  
  
    massBinStr = string(massBins[j],format='(f4.1)') + '-' + string(massBins[j+1],format='(f4.1)')
    strings = [strings,massBinStr]
    colors = [colors,getColor(j,/name)]
    
    ; histogram gadget (both) differences for current Tvir
    vals = [10.0^maxTemp.gal[w_gal]/10.0^accTvir.gal[w_gal],$
            10.0^maxTemp.gmem[w_gmem]/10.0^accTvir.gmem[w_gmem]]
    hist_both = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
    cgPlot,loc,float(hist_both)/total(hist_both),line=0,color=getColor(j),/overplot
    
    ; histogram gadget (gal) differences for current Tvir
    vals = [10.0^maxTemp.gal[w_gal]/10.0^accTvir.gal[w_gal]]
    hist_gal = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
    cgPlot,loc,float(hist_gal)/total(hist_both),line=1,color=getColor(j),/overplot
  
    ; histogram gadget (gmem) differences for current Tvir
    vals = [10.0^maxTemp.gmem[w_gmem]/10.0^accTvir.gmem[w_gmem]]
    hist_gmem = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
    cgPlot,loc,float(hist_gmem)/total(hist_both),line=2,color=getColor(j),/overplot
  
  endfor
  
  ; legend
  legend,['gal','gmem','both'],linestyle=[1,2,0],linesize=0.25,box=0,/left,/top
  if n_elements(massBins) gt 2 then $
    legend,strings,textcolors=colors,box=0,/right,/top
  
  end_PS
  stop
end

; plotAccTimeVsTmaxTime(): plot the offset between the time of accretion and the time of max prev. temp

pro plotAccTimeVsTmaxTime, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  binsizeGyr = 0.15 / (sP.res/128)

  ; config
  sgSelect = 'pri'
  ;massBins = [0.0,100.0] ; no massbins
  massBins = [9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5] ; log(M)  
  
  accTime = gcSubsetProp(sP=sP,select=sgSelect,/accTime,/mergerTreeSubset,/accretionTimeSubset)
  maxTempTime = gcSubsetProp(sP=sP,select=sgSelect,/maxTempTime,/mergerTreeSubset,/accretionTimeSubset)
  
  ; convert redshifts to age of universe
  accTime.gal      = redshiftToAgeFlat(accTime.gal[0,*])
  accTime.gmem     = redshiftToAgeFlat(accTime.gmem[0,*])
  maxTempTime.gal  = redshiftToAgeFlat(maxTempTime.gal)
  maxTempTime.gmem = redshiftToAgeFlat(maxTempTime.gmem)
    
  ; load parent halo masses so we can make halo massbins
  parentMass = gcSubsetProp(sP=sP,select=sgSelect,/parMass,/mergerTreeSubset,/accretionTimeSubset)

  ; plot
  start_PS, sP.plotPath + 'acctime_tmaxtime_' + sP.run + '.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    xrange = [-1,1.5]
    yrange = [1e-3,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="fraction",xtitle=textoidl("t_{acc} - t_{T_{max}}")+" [Gyr]",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+$
      ") gal=solid gmem=dot"
    cgPlot,[0,0],yrange,line=0,color=fsc_color('light gray'),/overplot
    
  strings = []
  colors = []
  
  for j=0,n_elements(massBins)-2 do begin
    
    ; select members of this parent mass bins and r>0<inf
    wGal  = where(parentMass.gal gt massBins[j] and parentMass.gal le massBins[j+1],count1)
    wGmem = where(parentMass.gmem gt massBins[j] and parentMass.gmem le massBins[j+1],count2)
    
    print,j,count1,count2
    
    if count1 eq 0 or count2 eq 0 then continue ; no halos in this mass bin  
  
    massBinStr = string(massBins[j],format='(f4.1)') + '-' + string(massBins[j+1],format='(f4.1)')
    strings = [strings,massBinStr]
    colors = [colors,j]
  
    ; histogram gal differences
    hist = histogram(accTime.gal[wGal]-maxTempTime.gal[wGal],binsize=binsizeGyr,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=0,color=getColor(j),/overplot

    ; histogram gmem differences
    hist = histogram(accTime.gmem[wGmem]-maxTempTime.gmem[wGmem],binsize=binsizeGyr,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=1,color=getColor(j),/overplot
  
  endfor
  
  ; legend
  legend,strings,textcolor=getColor(colors,/name),box=0,/right,/top
  
  end_PS
end

; plotVerticalSlices(): helper called by the subsequent 4 plot routines
  
pro plotVerticalSlices,rad_gal,rad_gmem,temp_gal,temp_gmem,plotName,binSizeTemp,temprange,xtitle

  compile_opt idl2, hidden, strictarr, strictarrsubs
  !except = 0 ;suppress floating point underflow/overflow errors
  
  ; config
  radBinEdges = [0.0,0.1,0.2,0.5,0.8,1.0,1.5] ;quasi-log spacing
  
  start_PS, plotName, xs=8, ys=6
  
    !p.multi = [0,3,2]
    
    !x.margin -= [4.0,1.5]
    !y.margin -= [0.5,0.5]
  
    for i=0,n_elements(radBinEdges)-2 do begin
      radBin = [radBinEdges[i],radBinEdges[i+1]]
      fsc_plot,[0],[0],/nodata,xrange=temprange,yrange=[0,1],/xs,/ys,ytitle="",xtitle=xtitle               
    
      ; select in radial bin
      w = where(rad_gal ge radBin[0] and rad_gal lt radBin[1],count_gal)
      if count_gal gt 0 then $  
        h1_gal = histogram(temp_gal[w],binsize=binSizeTemp,min=temprange[0],max=temprange[1],locations=xpts)
                          
      w = where(rad_gmem ge radBin[0] and rad_gmem lt radBin[1],count_gmem)
      if count_gmem gt 0 then $                 
        h1_gmem = histogram(temp_gmem[w],binsize=binSizeTemp,min=temprange[0],max=temprange[1],locations=xpts)

      ; skip when we have no points in this bin
      if ~count_gal and ~count_gmem then continue
      if n_elements(h1_gal) eq 0 then h1_gal = xpts*0
      if n_elements(h1_gmem) eq 0 then h1_gmem = xpts*0

      ; move xpts to bin centers and normalize histograms
      xpts += binSizeTemp/2.0
      
      h1_both = h1_gal + h1_gmem
      normfac = 1.0 / max(h1_both)
      
      h1_gal  *= normfac
      h1_gmem *= normfac
      h1_both *= normfac
      
      ; plot
      fsc_plot,xpts,h1_gal,color=getColor(16),thick=!p.thick+2,/overplot
      fsc_plot,xpts,h1_gmem,color=getColor(17),thick=!p.thick+2,/overplot
      fsc_plot,xpts,h1_both,color=getColor(0),line=1,thick=!p.thick+2,/overplot
      
      binStr = string(radBin[0],format='(f3.1)') + "-" + string(radBin[1],format='(f3.1)')
      fsc_text,(temprange[1]-temprange[0])*0.82,0.9,binStr,alignment=0.5,charsize=!p.charsize-0.5
    endfor
  
  !p.multi = 0
  
  end_PS
  !except = 1
end

; plot2DRadHistos(): helper function used by the subsequent 4 plot routines
  
pro plot2DRadHistos,plotBase,sP,h2rt_gal,h2rt_gmem,xrange,yrange,ytitle,$
                    virTempRange=virTempRange,massBinStr=massBinStr

  compile_opt idl2, hidden, strictarr, strictarrsubs
  !except = 0 ;suppress floating point underflow/overflow errors
  
  if sP.trMCPerCell eq 0 then  titleName = 'gadget'
  if sP.trMCPerCell eq -1 then titleName = 'tracerVel'
  if sP.trMCPerCell gt 0 then  titleName = 'tracerMC'
  
  h2rt_both = h2rt_gal + h2rt_gmem
  
  exp = 0.5 ; gamma exponent for non-linear color scaling
  ndivs = 5 ; number of divisions on colorbar  
  
  redshift = snapNumToRedshift(sP=sP)

  start_PS, sP.plotPath + plotBase + '_rad_both.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    ;loadct, 33, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
    ;cgLoadCT,13,/brewer,bottom=1
    loadColorTable, 'helix', /reverse
    
    tvim,h2rt_both^exp,pcharsize=!p.charsize,scale=0,clip=-1,$;,/c_map,range=[5e10,1e8,5e9];,/rct
         xtitle="r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}"),ytitle=ytitle,$
         barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=xrange,yrange=yrange,$
         /xlog,xticks=5,xtickv=[0.01,0.05,0.1,0.2,0.5,1.0],$
         xtickname=['0.01','0.05','0.1','0.2','0.5','1'],xmargin=2.0
         ;title=str(sP.res)+textoidl("^3")+" "+titleName+" (z="+string(redshift,format='(f3.1)')+") gal+gmem",$

     barvals = findgen(ndivs+1)/ndivs*(max(h2rt_both^exp)-min(h2rt_both^exp)) + min(h2rt_both^exp)
     ticknames = textoidl('10^{'+str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)'))+'}')
     ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
     colorbar,bottom=1,range=minmax(h2rt_both),position=[0.83,0.155,0.87,0.925],$
       /vertical,/right,title=textoidl("log ( M_{tot} )  [_{ }M_{sun }]"),divisions=ndivs,ticknames=ticknames,ncolors=255
            
      if keyword_set(virTempRange) then begin
        fsc_text,xrange[1]*0.45,(yrange[1]-yrange[0])*0.92+yrange[0],massBinStr,alignment=0.5,color=fsc_color('black')
        fsc_text,xrange[0]*1.6,virTempRange[1]*1.02,"T"+textoidl("_{vir}"),alignment=0.5,color=fsc_color('yellow')
        fsc_plot,xrange,[virTempRange[0],virTempRange[0]],line=0,color=fsc_color('yellow'),/overplot
        fsc_plot,xrange,[virTempRange[1],virTempRange[1]],line=0,color=fsc_color('yellow'),/overplot
      endif
             
  end_PS
  
  start_PS, sP.plotPath + plotBase + '_rad_gal.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    ;loadct, 33, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
    ;cgLoadCT,13,/brewer,bottom=1
    loadColorTable, 'helix', /reverse
    
    tvim,h2rt_gal^exp,pcharsize=!p.charsize,scale=0,clip=-1,$;,/c_map
         xtitle="r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}"),ytitle=ytitle,$
         barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=xrange,yrange=yrange,$
         /xlog,xticks=5,xtickv=[0.01,0.05,0.1,0.2,0.5,1.0],$
         xtickname=['0.01','0.05','0.1','0.2','0.5','1'],xmargin=2.0
         ;title=str(sP.res)+textoidl("^3")+" "+titleName+" (z="+string(redshift,format='(f3.1)')+") gal",$

     barvals = findgen(ndivs+1)/ndivs*(max(h2rt_gal^exp)-min(h2rt_gal^exp)) + min(h2rt_gal^exp)
     ticknames = textoidl('10^{'+str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)'))+'}')
     ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
     colorbar,bottom=1,range=minmax(h2rt_gal),position=[0.83,0.155,0.87,0.925],$
       /vertical,/right,title=textoidl("log ( M_{tot} )  [_{ }M_{sun }]"),divisions=ndivs,ticknames=ticknames,ncolors=255
         
      if keyword_set(virTempRange) then begin
        fsc_text,xrange[1]*0.45,(yrange[1]-yrange[0])*0.92+yrange[0],massBinStr,alignment=0.5,color=fsc_color('black')
        fsc_text,xrange[0]*1.6,virTempRange[1]*1.02,"T"+textoidl("_{vir}"),alignment=0.5,color=fsc_color('yellow')
        fsc_plot,xrange,[virTempRange[0],virTempRange[0]],line=0,color=fsc_color('orange'),/overplot
        fsc_plot,xrange,[virTempRange[1],virTempRange[1]],line=0,color=fsc_color('orange'),/overplot
      endif
      
  end_PS
  
  start_PS, sP.plotPath + plotBase + '_rad_gmem.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    ;loadct, 33, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
    ;cgLoadCT,13,/brewer,bottom=1
    loadColorTable, 'helix', /reverse
    
    tvim,h2rt_gmem^exp,pcharsize=!p.charsize,scale=0,clip=-1,$;,/c_map
         xtitle="r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}"),ytitle=ytitle,$
         barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=xrange,yrange=yrange,$
         /xlog,xticks=5,xtickv=[0.01,0.05,0.1,0.2,0.5,1.0],$
         xtickname=['0.01','0.05','0.1','0.2','0.5','1'],xmargin=2.0
         ;title=str(sP.res)+textoidl("^3")+" "+titleName+" (z="+string(redshift,format='(f3.1)')+") gmem",$

     barvals = findgen(ndivs+1)/ndivs*(max(h2rt_gmem^exp)-min(h2rt_gmem^exp)) + min(h2rt_gmem^exp)
     ticknames = textoidl('10^{'+str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)'))+'}')
     ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
     colorbar,bottom=1,range=minmax(h2rt_gmem),position=[0.83,0.155,0.87,0.925],$
       /vertical,/right,title=textoidl("log ( M_{tot} )  [_{ }M_{sun }]"),divisions=ndivs,ticknames=ticknames,ncolors=255
         
      if keyword_set(virTempRange) then begin
        fsc_text,xrange[1]*0.45,(yrange[1]-yrange[0])*0.92+yrange[0],massBinStr,alignment=0.5,color=fsc_color('black')
        fsc_text,xrange[0]*1.6,virTempRange[1]*1.02,"T"+textoidl("_{vir}"),alignment=0.5,color=fsc_color('yellow')
        fsc_plot,xrange,[virTempRange[0],virTempRange[0]],line=0,color=fsc_color('orange'),/overplot
        fsc_plot,xrange,[virTempRange[1],virTempRange[1]],line=0,color=fsc_color('orange'),/overplot
      endif
      
  end_PS
  !except = 1
end

; plotTempRad2DHisto(): plot 2d histogram of gas temperature (log Kelvin) as a function of r_gas/r_vir
;                       as well as vertical slices

pro plotTempRad2DHisto, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; config
  ;massBins = [0.0,100.0] ; no massbins
  massBins = [9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5] ; log(M)
  
  ; select one of:
  curTemp     = 0 ; current temperature
  maxPastTemp = 1 ; maximum previous temperature
  maxTempTime = 0 ; time of maximum previous temperature
  accTime     = 0 ; redshift of accretion across virial radius of parent halo
  accTvir     = 0 ; virial temperature of parent halo at time of accretion
  
   ; if maxPastTemp or maxTempTime, choose one of: population min,mean,max to one statistic per
   ; gas cell instead of all the tracers individually
   trPopMax   = 0
   trPopMin   = 0 
   trPopMean  = 0
   
  tVirNorm    = 0  ; normalize temperature by virial temp of parent halo at the starting time
  tVirAccNorm = 1  ; normalize temperature by virial temp of parent halo at the time of accretion
   
  parNorm  = 'pri' ; pri,sec (normalize r_vir and temp (if doing tVirNorm) by primary or secondary parent)
  sgSelect = 'pri' ; pri,sec,all (subhalo selection config, only plot gas in this category of subhalos)
  
  ; consider only the subset with recorded accretion times for certain types of plots
  mergerTreeSubset    = 0
  accretionTimeSubset = 0
  
  if accTime or accTvir or tVirAccNorm then begin
    mergerTreeSubset    = 1
    accretionTimeSubset = 1
  endif
  
  ; sanity check config parameters
  if accTime or accTvir then if sP.trMCPerCell eq -1 then message,'vel tracers not implemented yet'
  if accTime or accTvir then if tVirNorm then message,'no'
  if accTime or accTvir or tVirAccNorm then if parNorm ne 'pri' or sgSelect ne 'pri' then message,'better check you mean this'
  if trPopMax or trPopMin or trPopMean then if ~maxPastTemp and ~maxTempTime then message,'Error'
  if total(curTemp+maxPastTemp+maxTempTime+accTime+accTvir) ne 1 then message,'Just one'
  if total(tVirNorm+tVirAccNorm) gt 1 then message,'Only one allowed'
  if tVirNorm or tVirAccNorm then if ~curTemp and ~maxPastTemp then message,'Normalize what?'

  ; load group catalog just for counts of objects in each mass bin
  gc = loadGroupCat(sP=sP,/skipIDs)
  subgroupMasses = codeMassToLogMsun(gc.subgroupMass)
  gc = !NULL
  
  ; get normalized r/r_vir
  gcRad = gcSubsetProp(sP=sP,select=sgSelect,/rVirNorm,parNorm=parNorm,$
                       mergerTreeSubset=mergerTreeSubset,accretionTimeSubset=accretionTimeSubset)

  ; get current or maxPast temperature, or accretion time / tvir at accretion
  gcTemp = gcSubsetProp(sP=sP,select=sgSelect,curTemp=curTemp,maxPastTemp=maxPastTemp,$
                        maxTempTime=maxTempTime,accTime=accTime,accTvir=accTvir,$
                        mergerTreeSubset=mergerTreeSubset,accretionTimeSubset=accretionTimeSubset,$
                        trPopMin=trPopMin,trPopMax=trPopMax,trPopMean=trPopMean)

  ; normalize by halo Tvir at current time
  if tVirNorm then begin
    ; calculate temperatures of parents and take (log) ratio
    gcVirTemp = gcSubsetProp(sP=sP,select=sgSelect,/virTemp,parNorm=parNorm)
    
    gcTemp.gal  = alog10(10.0^gcTemp.gal / 10.0^gcVirTemp.gal)
    gcTemp.gmem = alog10(10.0^gcTemp.gmem / 10.0^gcVirTemp.gmem)
  endif
  
  ; normalize by halo Tvir at time of accretion
  if tVirAccNorm then begin
    ; calculate temperatures of parents at time of accretion and take (log) ratio
    gcAccTvir = gcSubsetProp(sP=sP,select=sgSelect,/accTvir,$
                             mergerTreeSubset=mergerTreeSubset,accretionTimeSubset=accretionTimeSubset)
    
    gcTemp.gal  = alog10(10.0^gcTemp.gal / 10.0^gcAccTvir.gal)
    gcTemp.gmem = alog10(10.0^gcTemp.gmem / 10.0^gcAccTvir.gmem)
  endif
  
  ; load gas masses if necessary (sph or byGas arepo)
  if sP.trMCPerCell eq 0 then $
    gcMass = gcSubsetProp(sP=sP,select=sgSelect,/curSingleVal,singleValField='mass',$
                          mergerTreeSubset=mergerTreeSubset,accretionTimeSubset=accretionTimeSubset)
  
  ; calculate masses of parents (for mass binning halos only)
  parentMass = gcSubsetProp(sP=sP,select=sgSelect,/parMass,$
                            mergerTreeSubset=mergerTreeSubset,accretionTimeSubset=accretionTimeSubset)

  for j=0,n_elements(massBins)-2 do begin
  
    ; plot config
    xrange = alog10([0.01,1.0])
    yrange = [4.0,7.0]
  
    binSizeRad  = 0.014 / (sP.res/128) ;0.04
    binSizeTemp = 0.05 / (sP.res/128) ;0.04
    
    ; preserve number of bins in log(rad) histogram and if changing yrange
    nBinsRad_linear = ceil((10.0^xrange[1]-10.0^xrange[0])/binSizeRad)+1
    nBinsTemp = ((yrange[1]-yrange[0])/binSizeTemp)+1
    binSizeRad_log  = (xrange[1]-xrange[0])/(nBinsRad_linear-1)
    
    if tVirNorm or tVirAccNorm then begin
      yrange = [-2.0,1.0]
      binSizeTemp = (yrange[1]-yrange[0])/(nBinsTemp-1)
    endif
    
    if accTime then begin
      yrange = [sP.redshift,4.0]
      binSizeTemp = (yrange[1]-yrange[0])/(nBinsTemp-1)
    endif
    
    ; select members of this parent mass bins and r>0<inf
    wGal  = where(parentMass.gal gt massBins[j] and parentMass.gal le massBins[j+1] and $
                  gcRad.gal gt 0.0 and finite(gcRad.gal),count1)
    wGmem = where(parentMass.gmem gt massBins[j] and parentMass.gmem le massBins[j+1] and $
                  gcRad.gmem gt 0.0 and finite(gcRad.gmem),count2)
    
    ; count for output
    wGCMassBin = where(subgroupMasses gt massBins[j] and subgroupMasses le massBins[j+1],count_gc)
    
    count_mt = 0
    if n_elements(mt) gt 0 then $
      wMTMassBin = where(hMass gt massBins[j] and hMass le massBins[j+1],count_mt)
    
    print,j,count1,count2,count_gc,count_mt
    
    ; select all particles/tracers in this mass bin
    temp_gal  = gcTemp.gal[wGal]
    temp_gmem = gcTemp.gmem[wGmem]
    
    rad_gal   = alog10( gcRad.gal[wGal] )
    rad_gmem  = alog10( gcRad.gmem[wGmem] )

    if count1 eq 0 or count2 eq 0 then continue ; no halos in this mass bin
    
    ; do mass weighting
    if sP.trMCPerCell eq 0 then begin
      ; plotting by gas cell, load gas mass subsets
      mass_gal  = gcMass.gal[wGal] * units.UnitMass_in_Msun
      mass_gmem = gcMass.gmem[wGmem] * units.UnitMass_in_Msun
      
      h2rt_gal = hist_nd_weight(transpose([[rad_gal],[temp_gal]]),weight=mass_gal,$
        [binSizeRad_log,binSizeTemp],$
        min=[xrange[0]-binSizeRad_log*0.5,yrange[0]-binSizeTemp*0.5],$
        max=[xrange[1]+binSizeRad_log*0.49,yrange[1]+binSizeTemp*0.49])
      h2rt_gmem = hist_nd_weight(transpose([[rad_gmem],[temp_gmem]]),weight=mass_gmem,$
        [binSizeRad_log,binSizeTemp],$
        min=[xrange[0]-binSizeRad_log*0.5,yrange[0]-binSizeTemp*0.5],$
        max=[xrange[1]+binSizeRad_log*0.49,yrange[1]+binSizeTemp*0.49])
    endif else begin
      ; create unweighted 2d histo (separate for galaxy/group member and composite)
      h2rt_gal  = hist_nd(transpose([[rad_gal],[temp_gal]]),[binSizeRad_log,binSizeTemp],$
        min=[xrange[0]-binSizeRad_log*0.5,yrange[0]-binSizeTemp*0.5],$
        max=[xrange[1]+binSizeRad_log*0.49,yrange[1]+binSizeTemp*0.49])
      h2rt_gmem = hist_nd(transpose([[rad_gmem],[temp_gmem]]),[binSizeRad_log,binSizeTemp],$
        min=[xrange[0]-binSizeRad_log*0.5,yrange[0]-binSizeTemp*0.5],$
        max=[xrange[1]+binSizeRad_log*0.49,yrange[1]+binSizeTemp*0.49])
                
      ; plotting all tracers, multiply number histogram by constant tracer mass
      if sP.trMCPerCell gt 0 then begin
        h2rt_gal  *= sP.trMassConst * units.UnitMass_in_Msun
        h2rt_gmem *= sP.trMassConst * units.UnitMass_in_Msun
      endif
      if sP.trMCPerCell eq -1 then begin
        h2rt_gal  *= sP.targetGasMass * units.UnitMass_in_Msun
        h2rt_gmem *= sP.targetGasMass * units.UnitMass_in_Msun
      endif
    endelse

    ; 2d histo plot config
    if curTemp then begin
      plotBase = "tcur_"+sgSelect+'_'+parNorm
      ytitle   = "log ( T"+textoidl("_{cur}")+" )"
    endif
    if maxPastTemp then begin
      plotBase = "tmax_"+sgSelect+'_'+parNorm
      ytitle   = "log ( T"+textoidl("_{max}")+" )"
    endif
    
    if maxTempTime then begin
      plotBase = "tmaxtime_"+sgSelect+'_'+parNorm
      ytitle   = "Maximum Temperature Redshift"
    endif
    if accTime then begin
      plotBase = "acctime_"+sgSelect+'_'+parNorm
      ytitle   = "Accretion Redshift"
    endif
    if accTvir then begin
      plotBase = "acctvir_"+sgSelect+'_'+parNorm
      ytitle   = "Parent "+textoidl("T_{vir}")+" at "+textoidl("t_{acc}")
    endif
    
    if tVirNorm then begin
      plotBase = strmid(plotBase,0,4)+"_tvirNorm_"+sgSelect+'_'+parNorm
      ytitle   = strmid(ytitle,0,strlen(ytitle)-2) + " / T"+textoidl("_{vir,cur}"+" )")
    endif
    if tVirAccNorm then begin
      plotBase = strmid(plotBase,0,4)+"_tvirAccNorm_"+sgSelect+'_'+parNorm
      ytitle   = strmid(ytitle,0,strlen(ytitle)-2) + " / T"+textoidl("_{vir,acc}"+" )")
    endif
    
    ; extra config for mass bins
    massBinStr   = !NULL
    virTempRange = !NULL
    
    if n_elements(massBins) gt 2 then begin
      plotBase = plotBase+"_mbin="+str(j)
      
      massBinStr = string(massBins[j],format='(f4.1)') + '-' + string(massBins[j+1],format='(f4.1)')
      
      massRangeCode = 10.0^[massBins[j],massBins[j+1]] / 1e10
      virTempRange = alog10( codeMassToVirTemp(massRangeCode,sP=sP) )
    endif
    
    ; plot 2d histo
    plot2DRadHistos,plotBase,sP,h2rt_gal,h2rt_gmem,10.0^xrange,yrange,ytitle,$
                    virTempRange=virTempRange,massBinStr=massBinStr    
    
    ; vertical slices plot config
    plotName = sP.plotPath + plotBase + '_slices.'+sP.run+'.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    if n_elements(massBins) gt 2 then $
      plotName = strmid(plotName,0,strlen(plotName)-4) + '.mbin='+str(j)+'.eps'
    
    xtitle = ytitle
    
    ; plot vertical slices
    plotVerticalSlices,rad_gal,rad_gmem,temp_gal,temp_gmem,plotName,binSizeTemp*2.0,yrange,xtitle
    
  endfor ;j
  
end
