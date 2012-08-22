; plotMaxTemps.pro
; gas accretion project - plots related to maximum past temperature of the gas
; dnelson aug.2012

; plotTmaxHistos(); plot (1) the previous max temp normalized by tviracc for arepo vs. gadget, gal vs. 
;                   gmem, (2) same but unnormalized by tviracc, (3) global not binned by halo mass but
;                   normalized by each parent tviracc, (4) same global without normalization

pro plotTmaxHistos

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  sPg = simParams(res=256,run='gadget',redshift=2.0)
  sPa = simParams(res=256,run='tracer',redshift=2.0) ; f=-1 use velocity tracers

  binSizeLog = 0.15 / (sPg.res/128)
  
  sgSelect  = 'pri'
  accModes  = ['all','smooth','clumpy','stripped']
  massBins  = [9.0,9.5,10.0,10.5,11.0,11.5,12.0] ; log(M)    
  
  lines   = [1,0,2] ; tvircur,tviracc,const5.5
  thicks  = [4,6,8] ; 128,256,512
  cInd    = 1 ; color index
  
  foreach accMode,accModes do begin
   
    print,accMode
  
  ; load sPg (gadget)
  accTvir_gadget = gcSubsetProp(sP=sPg,select=sgSelect,/accTvir,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  curTvir_gadget = gcSubsetProp(sP=sPg,select=sgSelect,/virTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  maxTemp_gadget = gcSubsetProp(sP=sPg,select=sgSelect,/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)  
 
  ; load parent halo masses so we can make halo massbins
  parentMass_ga = gcSubsetProp(sP=sPg,select=sgSelect,/parMass,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)

  ; load sPa (tracer)
  accTvir_tracer = gcSubsetProp(sP=sPa,select=sgSelect,/accTvir,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  curTvir_tracer = gcSubsetProp(sP=sPa,select=sgSelect,/virTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
  maxTemp_tracer = gcSubsetProp(sP=sPa,select=sgSelect,/maxPastTemp,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)

  parentMass_tr = gcSubsetProp(sP=sPa,select=sgSelect,/parMass,/mergerTreeSubset,/accretionTimeSubset,accMode=accMode)
    
  ; plot (1) - 3x2 mass bins separated out and each panel with gadget+arepo, gal vs. gmem
  start_PS, sPg.plotPath + 'tmax_3x2_tviracc.'+accMode+'.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.res)+'_'+str(sPg.snap)+'.eps', /big
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
    
    for j=0,n_elements(massBins)-2 do begin
      
      if j eq 0 or j eq 3 then ytickname = '' else ytickname = replicate(' ',10)
      if j ge 3 then xtickname = '' else xtickname = replicate(' ',10)
      if j gt 0 then noerase = 1 else noerase = 0
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,pos=pos[j],$
        ytitle="",xtitle="",xticks=3,xtickv=xtickv,xtickname=xtickname,ytickname=ytickname,noerase=noerase       
      
      cgPlot,[0,0],[8e-4,0.25],line=2,color=fsc_color('black'),thick=!p.thick-0.0,/overplot
      
      ; select members of this parent mass bins and r>0<inf
      wGadget_gal  = where(parentMass_ga.gal gt massBins[j] and parentMass_ga.gal le massBins[j+1],count1)
      wGadget_gmem = where(parentMass_ga.gmem gt massBins[j] and parentMass_ga.gmem le massBins[j+1],count2)
      wTracer_gal  = where(parentMass_tr.gal gt massBins[j] and parentMass_tr.gal le massBins[j+1],count3)
      wTracer_gmem = where(parentMass_tr.gmem gt massBins[j] and parentMass_tr.gmem le massBins[j+1],count4)
      
      ;print,j,count1,count2,count3,count4
      
      if ~count1 or ~count2 or ~count3 or ~count4 then continue ; no halos in this mass bin  
      
      ; histogram gadget (gal) differences
      vals = [10.0^maxTemp_gadget.gal[wGadget_gal]/10.0^accTvir_gadget.gal[wGadget_gal]]
      hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=lines[1],color=sPg.colorsG[cInd],/overplot
  
      ; histogram tracer (gal) differences
      vals = [10.0^maxTemp_tracer.gal[wTracer_gal]/10.0^accTvir_tracer.gal[wTracer_gal]]
      hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=lines[1],color=sPa.colorsA[cInd],/overplot
    
      ; histogram gadget (gmem) differences
      vals = [10.0^maxTemp_gadget.gmem[wGadget_gmem]/10.0^accTvir_gadget.gmem[wGadget_gmem]]
      hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=lines[0],color=sPg.colorsG[cInd],/overplot
  
      ; histogram tracer (gmem) differences
      vals = [10.0^maxTemp_tracer.gmem[wTracer_gmem]/10.0^accTvir_tracer.gmem[wTracer_gmem]]
      hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=lines[0],color=sPa.colorsA[cInd],/overplot
    
      ; legends
      massBinStr = string(massBins[j],format='(f4.1)') + ' < log(M) < ' + $
                   string(massBins[j+1],format='(f4.1)')

      cgText,mean(xrange),yrange[1]*0.4,massBinStr,charsize=!p.charsize-0.0,alignment=0.5
          
      if j eq 0 then $
        legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,$
          position=[-2.0,0.1],charsize=!p.charsize-0.1,spacing=!p.charsize+0.5
    
    endfor
    
    legend,['central galaxy','halo atmosphere'],linestyle=[0,1],$
      color=cgColor('dark gray'),textcolors=['dark gray','dark gray'],$
      linesize=0.25,box=0,position=[0.08,0.09],/normal,charsize=!p.charsize-0.2
    
    ; axis labels
    cgText,0.05,mean([y0,y2]),"Gas Mass Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,mean([x0,x3]),0.05,textoidl("log ( T_{max} / T_{vir,acc} )"),alignment=0.5,/normal
    
  end_PS
  
  ; plot (2) - 3x2 mass bins separated out and each panel with gadget+arepo, gal vs. gmem
  start_PS, sPg.plotPath + 'tmax_3x2.'+accMode+'.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.res)+'_'+str(sPg.snap)+'.eps', /big
    !p.thick += 1
    xrange = [3.8,7.2]
    yrange = [6e-4,1.0]
    
    xtickv = [4.0,5.0,6.0,7.0]
    
    for j=0,n_elements(massBins)-2 do begin
      
      if j eq 0 or j eq 3 then ytickname = '' else ytickname = replicate(' ',10)
      if j ge 3 then xtickname = '' else xtickname = replicate(' ',10)
      if j gt 0 then noerase = 1 else noerase = 0
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,pos=pos[j],$
        ytitle="",xtitle="",xticks=3,xtickv=xtickv,xtickname=xtickname,ytickname=ytickname,noerase=noerase       
      
      cgPlot,[0,0],[8e-4,0.25],line=2,color=fsc_color('black'),thick=!p.thick-0.0,/overplot
      
      ; select members of this parent mass bins and r>0<inf
      wGadget_gal  = where(parentMass_ga.gal gt massBins[j] and parentMass_ga.gal le massBins[j+1],count1)
      wGadget_gmem = where(parentMass_ga.gmem gt massBins[j] and parentMass_ga.gmem le massBins[j+1],count2)
      wTracer_gal  = where(parentMass_tr.gal gt massBins[j] and parentMass_tr.gal le massBins[j+1],count3)
      wTracer_gmem = where(parentMass_tr.gmem gt massBins[j] and parentMass_tr.gmem le massBins[j+1],count4)
      
      ;print,j,count1,count2,count3,count4
      
      if ~count1 or ~count2 or ~count3 or ~count4 then continue ; no halos in this mass bin  
      
      ; histogram gadget (gal) differences
      hist = histogram(maxTemp_gadget.gal[wGadget_gal],binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=lines[1],color=sPg.colorsG[cInd],/overplot
  
      ; histogram tracer (gal) differences
      hist = histogram(maxTemp_tracer.gal[wTracer_gal],binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=lines[1],color=sPa.colorsA[cInd],/overplot
    
      ; histogram gadget (gmem) differences
      hist = histogram(maxTemp_gadget.gmem[wGadget_gmem],binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=lines[0],color=sPg.colorsG[cInd],/overplot
  
      ; histogram tracer (gmem) differences
      hist = histogram(maxTemp_tracer.gmem[wTracer_gmem],binsize=binsizeLog,loc=loc)
      cgPlot,loc+binsizeLog*0.5,float(hist)/total(hist),line=lines[0],color=sPa.colorsA[cInd],/overplot
    
      ; legends
      massBinStr = string(massBins[j],format='(f4.1)') + ' < log(M) < ' + $
                   string(massBins[j+1],format='(f4.1)')

      cgText,mean(xrange),yrange[1]*0.4,massBinStr,charsize=!p.charsize-0.0,alignment=0.5
          
      if j eq 0 then $
        legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,$
          position=[-2.0,0.1],charsize=!p.charsize-0.1,spacing=!p.charsize+0.5
    
    endfor
    
    legend,['central galaxy','halo atmosphere'],linestyle=[0,1],$
      color=cgColor('dark gray'),textcolors=['dark gray','dark gray'],$
      linesize=0.25,box=0,position=[0.08,0.09],/normal,charsize=!p.charsize-0.2
    
    ; axis labels
    cgText,0.05,mean([y0,y2]),"Gas Mass Fraction",alignment=0.5,orientation=90.0,/normal
    cgText,mean([x0,x3]),0.05,textoidl("log ( T_{max} )"),alignment=0.5,/normal
    
  end_PS  
  
  ; plot (3) - global histos of tmax/tviracc
  start_PS, sPg.plotPath + 'tmax_global_tviracc.'+accMode+'.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.res)+'_'+str(sPg.snap)+'.eps', /big
    ;!p.thick += 1
    
    xrange = [-2.0,1.5]
    yrange = [1e-3,1e3]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle=textoidl("M_{tot}"),xtitle=textoidl("log ( T_{max} / T_{vir,acc} )")
    
    ; histogram gadget (gal) differences
    vals = [10.0^maxTemp_gadget.gal/10.0^accTvir_gadget.gal]
    hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
    cgPlot,loc+binsizeLog*0.5,float(hist)*sPg.targetGasMass,line=lines[1],color=sPg.colorsG[cInd],/overplot

    ; histogram tracer (gal) differences
    vals = [10.0^maxTemp_tracer.gal/10.0^accTvir_tracer.gal]
    hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
    cgPlot,loc+binsizeLog*0.5,float(hist)*sPa.trMassConst,line=lines[1],color=sPa.colorsA[cInd],/overplot
  
    ; histogram gadget (gmem) differences
    vals = [10.0^maxTemp_gadget.gmem/10.0^accTvir_gadget.gmem]
    hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
    cgPlot,loc+binsizeLog*0.5,float(hist)*sPg.targetGasMass,line=lines[0],color=sPg.colorsG[cInd],/overplot

    ; histogram tracer (gmem) differences
    vals = [10.0^maxTemp_tracer.gmem/10.0^accTvir_tracer.gmem]
    hist = histogram(alog10(vals),binsize=binsizeLog,loc=loc)
    cgPlot,loc+binsizeLog*0.5,float(hist)*sPa.trMassConst,line=lines[0],color=sPa.colorsA[cInd],/overplot
  
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,$
      /top,/right,spacing=!p.charsize+0.5
    
    legend,['central galaxy','halo atmosphere'],linestyle=[0,1],$
      color=cgColor('dark gray'),textcolors=['dark gray','dark gray'],$
      linesize=0.25,box=0,/top,/left
    
  end_PS
  
  ; plot (4) - global histos of unnormalized tmax
  start_PS, sPg.plotPath + 'tmax_global.'+accMode+'.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
            str(sPg.res)+'_'+str(sPg.snap)+'.eps', /big
    ;!p.thick += 1
    
    xrange = [3.5,7.5]
    yrange = [1e-3,1e3]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,$
      ytitle=textoidl("M_{tot}"),xtitle=textoidl("log ( T_{max} )")
    
    ; histogram gadget (gal) differences
    hist = histogram(maxTemp_gadget.gal,binsize=binsizeLog,loc=loc)
    cgPlot,loc+binsizeLog*0.5,float(hist)*sPg.targetGasMass,line=lines[1],color=sPg.colorsG[cInd],/overplot

    ; histogram tracer (gal) differences
    hist = histogram(maxTemp_tracer.gal,binsize=binsizeLog,loc=loc)
    cgPlot,loc+binsizeLog*0.5,float(hist)*sPa.trMassConst,line=lines[1],color=sPa.colorsA[cInd],/overplot
  
    ; histogram gadget (gmem) differences
    hist = histogram(maxTemp_gadget.gmem,binsize=binsizeLog,loc=loc)
    cgPlot,loc+binsizeLog*0.5,float(hist)*sPg.targetGasMass,line=lines[0],color=sPg.colorsG[cInd],/overplot

    ; histogram tracer (gmem) differences
    hist = histogram(maxTemp_tracer.gmem,binsize=binsizeLog,loc=loc)
    cgPlot,loc+binsizeLog*0.5,float(hist)*sPa.trMassConst,line=lines[0],color=sPa.colorsA[cInd],/overplot
  
    legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,$
      /top,/right,spacing=!p.charsize+0.5
    
    legend,['central galaxy','halo atmosphere'],linestyle=[0,1],$
      color=cgColor('dark gray'),textcolors=['dark gray','dark gray'],$
      linesize=0.25,box=0,/top,/left
    
  end_PS
  
  ; plot (5) - same plot using current Tvir instead of Tvir at accretion time
  if 0 then begin
  start_PS, sPg.plotPath + 'tmax_tvircur_comp.'+str(sPg.res)+'_'+str(sPg.snap)+'.eps'
    !p.thick += 1
    xrange = [-2.5,1.0]
    yrange = [1e-3,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="Mass Fraction",xtitle=textoidl("log ( T_{max} / T_{vir,cur} )")+"",$
      title=str(res)+textoidl("^3")+" Gadget vs. ArepoMC (z="+string(sPg.redshift,format='(f3.1)')+")"
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
  endif ;0
  
  endforeach ; accModes
  
end
