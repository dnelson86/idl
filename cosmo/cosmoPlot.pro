; cosmoPlot.pro
; gas accretion project - plots
; dnelson mar.2012

; plotTmaxVsTvirAcc(); plot the previous max temp vs. the virial temperature of the parent halos at the
;                      time of accretion for arepo vs. gadget

pro plotTmaxVsTvirAcc

  ; config
  res = 256
  redshift = 2.0

  sP1 = simParams(res=res,run='gadget',redshift=redshift)
  sP2 = simParams(res=res,run='tracer',redshift=redshift)
  
  binSizeRatio   = 0.02 / (res/128)
  binSizeLogDiff = 0.1 / (res/128)
  
  ;massBins = [0.0,100.0] ; no massbins
  massBins = [9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5] ; log(M)   
  
  ; load sP1 (gadget)
  ; --------
  mt = mergerTreeSubset(sP=sP1)
  at = accretionTimes(sP=sP1)
  
  gal_w  = where(at.AccTime_gal ne -1,count_gal)
  gmem_w = where(at.AccTime_gmem ne -1,count_gmem)
    
  ; create a galcatSub which indexes the galcat for the subset of tracked halos
  galcatSub = galcatINDList(sP=sP1,gcIDList=mt.galcatIDList) ;identical to mt.galcatSub
  
  maxt = maxTemps(sP=sP1,/loadByGas)
  maxtemp_gadget = { gal  : maxt.maxTemps_gal[galcatSub.gal[gal_w]]     ,$
                     gmem : maxt.maxTemps_gmem[galcatSub.gmem[gmem_w]]   }
  
  acctvir_gadget = { gal  : at.accHaloTvir_gal[gal_w]     ,$
                     gmem : at.accHaloTvir_gmem[gmem_w]    }
  
  ; load parent halo masses so we can make halo massbins
  parentMass = gcSubsetProp(sP=sP1,select='pri',/parMass,allTR=0)
    
  parentMass_ga = { gal  : parentMass.gal[galcatSub.gal[gal_w]]    ,$
                    gmem : parentMass.gmem[galcatSub.gmem[gmem_w]]  }   
  
  maxt = !NULL
  at   = !NULL
  mt   = !NULL
  galcatSub  = !NULL
  parentMass = !NULL

  ; load sP2 (tracer)
  ; --------
  mt = mergerTreeSubset(sP=sP2)
  at = accretionTimes(sP=sP2)
  
  gal_w  = where(at.AccTime_gal ne -1,count_gal)
  gmem_w = where(at.AccTime_gmem ne -1,count_gmem)
  
  ; load child counts for all members of galcat
  maxt_gal = maxTemps(sP=sP2,/loadAllTRGal)
  child_counts_gal = maxt_gal.child_counts
  maxt_gmem = maxTemps(sP=sP2,/loadAllTRGmem)
  child_counts_gmem = maxt_gmem.child_counts
  
  ; get indices into subset of galcat found within accretionTimes()
  galcatSub = galCatINDList(sP=sP2,gcIDList=mt.galcatIDList,$
                             child_counts={gal:child_counts_gal,gmem:child_counts_gmem}) 
  
  maxtemp_tracer = { gal  : maxt_gal.maxTemps[galcatSub.gal[gal_w]]     ,$
                     gmem : maxt_gmem.maxTemps[galcatSub.gmem[gmem_w]]   }
                     
  acctvir_tracer = { gal  : at.accHaloTvir_gal[gal_w]     ,$
                     gmem : at.accHaloTvir_gmem[gmem_w]    }

  ; load parent halo masses so we can make halo massbins
  parentMass = gcSubsetProp(sP=sP2,select='all',/parMass,allTR=1)
    
  parentMass_tr = { gal  : parentMass.gal[galcatSub.gal[gal_w]]    ,$
                    gmem : parentMass.gmem[galcatSub.gmem[gmem_w]]  }   

  maxt = !NULL
  at   = !NULL
  mt   = !NULL
  galcatSub  = !NULL
  parentMass = !NULL

  ; plot (1)
  start_PS, sP1.plotPath + 'tmax_tviracc_comp.'+str(sP1.res)+'_'+str(sP1.snap)+'.eps'
    !p.thick += 1
    xrange = [-2.5,1.0]
    yrange = [5e-4,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="fraction",xtitle=textoidl("log(T_{max}) - log(T_{vir,acc})")+"",$
      title=str(res)+textoidl("^3")+" Gadget vs. ArepoMC (z="+string(sP1.redshift,format='(f3.1)')+")"
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
    vals = [maxtemp_gadget.gal[wGadget_gal]-acctvir_gadget.gal[wGadget_gal],$
            maxtemp_gadget.gmem[wGadget_gmem]-acctvir_gadget.gmem[wGadget_gmem]]
    hist = histogram(vals,binsize=binsizeLogDiff,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=0,color=getColor(j),/overplot

    ; histogram tracer (gal+gmem) differences
    vals = [maxtemp_tracer.gal[wTracer_gal]-acctvir_tracer.gal[wTracer_gal],$
            maxtemp_tracer.gmem[wTracer_gmem]-acctvir_tracer.gmem[wTracer_gmem]]
    hist = histogram(vals,binsize=binsizeLogDiff,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=2,color=getColor(j),/overplot
  
  endfor
  
  ; legend
  legend,['gadget','arepoMC'],linestyle=[0,2],linesize=0.25,box=0,/right,/top
  if n_elements(massBins) gt 2 then $
    legend,strings,textcolors=colors,box=0,/left,/top
  
  end_PS
  
  ; plot (2)
  start_PS, sP1.plotPath + 'tmax_tviracc_comp2.'+str(sP1.res)+'_'+str(sP1.snap)+'.eps'
    !p.thick += 1
    xrange = [0.6,1.2]
    yrange = [5e-4,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="fraction",xtitle=textoidl("log (T_{max}) / log(T_{vir,acc})")+"",$
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
    vals = [maxtemp_gadget.gal[wGadget_gal]/acctvir_gadget.gal[wGadget_gal],$
            maxtemp_gadget.gmem[wGadget_gmem]/acctvir_gadget.gmem[wGadget_gmem]]
    hist = histogram(vals,binsize=binsizeRatio,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=0,color=getColor(j),/overplot

    ; histogram tracer (gal+gmem) differences
    vals = [maxtemp_tracer.gal[wTracer_gal]/acctvir_tracer.gal[wTracer_gal],$
            maxtemp_tracer.gmem[wTracer_gmem]/acctvir_tracer.gmem[wTracer_gmem]]
    hist = histogram(vals,binsize=binsizeRatio,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=2,color=getColor(j),/overplot
  
  endfor
  
  ; legend
  legend,['gadget','arepoMC'],linestyle=[0,2],linesize=0.25,box=0,/right,/top
  if n_elements(massBins) gt 2 then $
    legend,strings,textcolors=colors,box=0,/left,/top
  
  end_PS
end

; plotAccTimeVsTmaxTime(): plot the offset between the time of accretion and the time of max prev. temp

pro plotAccTimeVsTmaxTime, sP=sP

  binsizeGyr = 0.15 / (sP.res/128)

  ; config
  sgSelect = 'all'
  allTR = 0
  if sP.trMCPerCell ne 0 then allTR = 1
  
  ;massBins = [0.0,100.0] ; no massbins
  massBins = [9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5] ; log(M)  
  
  ; select the subset of particles/tracers in the mergerTreeSubset
  mt = mergerTreeSubset(sP=sP)
  at = accretionTimes(sP=sP)
  
  gal_w  = where(at.AccTime_gal ne -1,count_gal)
  gmem_w = where(at.AccTime_gmem ne -1,count_gmem)
    
  ; create a galcatSub which indexes the galcat for the subset of tracked halos
  if sP.trMCPerCell eq 0 then begin
    galcatSub = galcatINDList(sP=sP,gcIDList=mt.galcatIDList) ;identical to mt.galcatSub
    
    maxt = maxTemps(sP=sP,/loadByGas)
    maxt_time_gal  = maxt.maxTempTime_gal[galcatSub.gal]
    maxt_time_gmem = maxt.maxTempTime_gmem[galcatSub.gmem]
  endif else begin
    ; load child counts for all members of galcat
    maxt_gal = maxTemps(sP=sP,/loadAllTRGal)
    child_counts_gal = maxt_gal.child_counts
    maxt_gmem = maxTemps(sP=sP,/loadAllTRGmem)
    child_counts_gmem = maxt_gmem.child_counts
    
    ; get indices into subset of galcat found within accretionTimes()
    galcatSub = galCatINDList(sP=sP,gcIDList=mt.galcatIDList,$
                               child_counts={gal:child_counts_gal,gmem:child_counts_gmem}) 
    
    maxt_time_gal  = maxt_gal.maxTempTime[galcatSub.gal]
    maxt_time_gmem = maxt_gmem.maxTempTime[galcatSub.gmem]
  endelse
  
  ; subset of accretion times and tmax times corresponding to those found within accretionTimes()
  ; convert to redshift
  acc_time_gal   = 1/at.accTime_gal[gal_w]-1
  acc_time_gmem  = 1/at.accTime_gmem[gmem_w]-1
  maxt_time_gal  = 1/maxt_time_gal[gal_w]-1
  maxt_time_gmem = 1/maxt_time_gmem[gmem_w]-1
  
  ; convert redshifts to age of universe
  acc_time_gal   = redshiftToAgeFlat(acc_time_gal)
  acc_time_gmem  = redshiftToAgeFlat(acc_time_gmem)
  maxt_time_gal  = redshiftToAgeFlat(maxt_time_gal)
  maxt_time_gmem = redshiftToAgeFlat(maxt_time_gmem)
    
  ; load parent halo masses so we can make halo massbins
  parentMass = gcSubsetProp(sP=sP,select=sgSelect,/parMass,allTR=allTR)
    
  parentMass = { gal  : parentMass.gal[galcatSub.gal[gal_w]]    ,$
                 gmem : parentMass.gmem[galcatSub.gmem[gmem_w]]  }    
    
  ; plot
  start_PS, sP.plotPath + 'acctime_tmaxtime_' + sP.run + '.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    xrange = [-1,1.5]
    yrange = [5e-4,1.0]
    
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
    hist = histogram(acc_time_gal[wGal]-maxt_time_gal[wGal],binsize=binsizeGyr,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=0,color=getColor(j),/overplot

    ; histogram gmem differences
    hist = histogram(acc_time_gmem[wGmem]-maxt_time_gmem[wGmem],binsize=binsizeGyr,loc=loc)
    cgPlot,loc,float(hist)/total(hist),line=1,color=getColor(j),/overplot
  
  endfor
  
  ; legend
  legend,strings,textcolor=getColor(colors,/name),box=0,/right,/top
  
  end_PS
stop
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
  
  plotName = sP.plotPath + plotBase + '_rad_both_'+str(sP.res)+'_'+str(sP.snap)+'.eps'
  
  if sP.trMCPerCell eq 0 then  titleName = 'gadget'
  if sP.trMCPerCell eq -1 then titleName = 'tracerVel'
  if sP.trMCPerCell gt 0 then  titleName = 'tracerMC'
  
  h2rt_both = h2rt_gal + h2rt_gmem
  
  exp = 0.5 ; gamma exponent for non-linear color scaling
  ndivs = 5 ; number of divisions on colorbar  
  
  redshift = snapNumToRedshift(sP=sP)

  start_PS, plotName
    
    loadct, 33, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
    ;cgLoadCT,13,/brewer,bottom=1
    ;loadColorTable, 'helix'
    
    tvim,h2rt_both^exp,pcharsize=!p.charsize,scale=0,clip=-1,$;,/c_map,range=[5e10,1e8,5e9];,/rct
         xtitle="r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}"),ytitle=ytitle,$
         title=str(sP.res)+textoidl("^3")+" "+titleName+" (z="+string(redshift,format='(f3.1)')+") gal+gmem",$
         barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=xrange,yrange=yrange,$
         /xlog,xticks=6,xtickv=[0.01,0.05,0.1,0.2,0.5,1.0,3.0],$
         xtickname=['0.01','0.05','0.1','0.2','0.5','1','3'],xmargin=2.0

     barvals = findgen(ndivs+1)/ndivs*(max(h2rt_both^exp)-min(h2rt_both^exp)) + min(h2rt_both^exp)
     ticknames = textoidl('10^{'+str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)'))+'}')
     ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
     colorbar,bottom=1,range=minmax(h2rt_both),position=[0.83,0.155,0.87,0.925],$
       /vertical,/right,title=textoidl("M_{sun,tot}"),divisions=ndivs,ticknames=ticknames,ncolors=255
            
      if keyword_set(virTempRange) then begin
        fsc_text,xrange[1]*0.45,(yrange[1]-yrange[0])*0.92+yrange[0],massBinStr,alignment=0.5,color=fsc_color('white')
        fsc_text,xrange[0]*1.6,virTempRange[1]*1.02,"T"+textoidl("_{vir}"),alignment=0.5,color=fsc_color('yellow')
        fsc_plot,xrange,[virTempRange[0],virTempRange[0]],line=0,color=fsc_color('yellow'),/overplot
        fsc_plot,xrange,[virTempRange[1],virTempRange[1]],line=0,color=fsc_color('yellow'),/overplot
      endif
             
  end_PS
  
  plotName = sP.plotPath + plotBase + '_rad_gal_'+str(sP.res)+'_'+str(sP.snap)+'.eps'
  
  start_PS, plotName
    
    loadct, 33, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
    ;cgLoadCT,13,/brewer,bottom=1
    ;loadColorTable, 'helix'
    
    tvim,h2rt_gal^exp,pcharsize=!p.charsize,scale=0,clip=-1,$;,/c_map
         xtitle="r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}"),ytitle=ytitle,$
         title=str(sP.res)+textoidl("^3")+" "+titleName+" (z="+string(redshift,format='(f3.1)')+") gal",$
         barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=xrange,yrange=yrange,$
         /xlog,xticks=6,xtickv=[0.01,0.05,0.1,0.2,0.5,1.0,3.0],$
         xtickname=['0.01','0.05','0.1','0.2','0.5','1','3'],xmargin=2.0

     barvals = findgen(ndivs+1)/ndivs*(max(h2rt_gal^exp)-min(h2rt_gal^exp)) + min(h2rt_gal^exp)
     ticknames = textoidl('10^{'+str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)'))+'}')
     ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
     colorbar,bottom=1,range=minmax(h2rt_gal),position=[0.83,0.155,0.87,0.925],$
       /vertical,/right,title=textoidl("M_{sun,tot}"),divisions=ndivs,ticknames=ticknames,ncolors=255
         
      if keyword_set(virTempRange) then begin
        fsc_text,xrange[1]*0.45,(yrange[1]-yrange[0])*0.92+yrange[0],massBinStr,alignment=0.5,color=fsc_color('white')
        fsc_text,xrange[0]*1.6,virTempRange[1]*1.02,"T"+textoidl("_{vir}"),alignment=0.5,color=fsc_color('yellow')
        fsc_plot,xrange,[virTempRange[0],virTempRange[0]],line=0,color=fsc_color('orange'),/overplot
        fsc_plot,xrange,[virTempRange[1],virTempRange[1]],line=0,color=fsc_color('orange'),/overplot
      endif
      
  end_PS
  
  plotName = sP.plotPath + plotBase + '_rad_gmem_'+str(sP.res)+'_'+str(sP.snap)+'.eps'
  
  start_PS, plotName
    
    loadct, 33, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
    ;cgLoadCT,13,/brewer,bottom=1
    ;loadColorTable, 'helix'
    
    tvim,h2rt_gmem^exp,pcharsize=!p.charsize,scale=0,clip=-1,$;,/c_map
         xtitle="r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}"),ytitle=ytitle,$
         title=str(sP.res)+textoidl("^3")+" "+titleName+" (z="+string(redshift,format='(f3.1)')+") gmem",$
         barwidth=0.75,lcharsize=!p.charsize-0.2,xrange=xrange,yrange=yrange,$
         /xlog,xticks=6,xtickv=[0.01,0.05,0.1,0.2,0.5,1.0,3.0],$
         xtickname=['0.01','0.05','0.1','0.2','0.5','1','3'],xmargin=2.0

     barvals = findgen(ndivs+1)/ndivs*(max(h2rt_gmem^exp)-min(h2rt_gmem^exp)) + min(h2rt_gmem^exp)
     ticknames = textoidl('10^{'+str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)'))+'}')
     ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
     colorbar,bottom=1,range=minmax(h2rt_gmem),position=[0.83,0.155,0.87,0.925],$
       /vertical,/right,title=textoidl("M_{sun,tot}"),divisions=ndivs,ticknames=ticknames,ncolors=255
         
      if keyword_set(virTempRange) then begin
        fsc_text,xrange[1]*0.45,(yrange[1]-yrange[0])*0.92+yrange[0],massBinStr,alignment=0.5,color=fsc_color('white')
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
  maxPastTemp = 0 ; maximum previous temperature
  accTime     = 1 ; redshift of accretion across virial radius of parent halo
  accTvir     = 0 ; virial temperature of parent halo at time of accretion
  
   ; if maxPastTemp=1, choose one of: population min,mean,max to return
   ; (choice of one requred for allTR=0, and ignored for allTR=1)
   trPopMax   = 0
   trPopMin   = 0 
   trPopMean  = 0
   
  tVirNorm    = 0  ; normalize temperature by virial temp of parent halo at the starting time
  tVirAccNorm = 0  ; normalize temperature by virial temp of parent halo at the time of accretion
   
  allTR    = 1     ; plot results for all tracers individually, not collapsed onto their parents
  parNorm  = 'pri' ; pri,sec (normalize r_vir and temp (if doing tVirNorm) by primary or secondary parent)
  sgSelect = 'all' ; pri,sec,all (subhalo selection config, only plot gas in this category of subhalos)
  
  ; sanity check config parameters
  if accTime or accTvir then if sP.trMCPerCell eq -1 then message,'vel tracers not implemented yet'
  if accTime or accTvir then if tVirNorm then message,'no'
  if accTime or accTvir or tVirAccNorm then if parNorm ne 'pri' or sgSelect ne 'all' then message,'better check you mean this'
  if trPopMax or trPopMin or trPopMean then if ~maxPastTemp then message,'Error'
  if total(curTemp+maxPastTemp+accTime+accTvir) ne 1 then message,'Just one'
  if allTR and sP.trMCPerCell eq 0 then message,'Error: AllTR set but no tracers in this run'
  if total(tVirNorm+tVirAccNorm) gt 1 then message,'Only one allowed'
  if tVirNorm or tVirAccNorm then if ~curTemp and ~maxPastTemp then message,'Normalize what?'
  if (tVirAccNorm or accTime or accTvir) and sP.trMCPerCell ne 0 and ~allTR then message,'tVirAccNorm with tracers only with allTR'
  
  ; load group catalog just for counts of objects in each mass bin
  gc = loadGroupCat(sP=sP,/skipIDs)
  subgroupMasses = codeMassToLogMsun(gc.subgroupMass)
  gc = !NULL
  
  ; get normalized r/r_vir
  gcRad = gcSubsetProp(sP=sP,select=sgSelect,/rVirNorm,allTR=allTR,parNorm=parNorm)
  
  ; get current or maxPast temperature
  if curTemp or maxPastTemp then $
    gcTemp = gcSubsetProp(sP=sP,select=sgSelect,curTemp=curTemp,maxPastTemp=maxPastTemp,allTR=allTR,$
                          trPopMin=trPopMin,trPopMax=trPopMax,trPopMean=trPopMean)

  if tVirNorm then begin
    ; calculate temperatures of parents and take (non-log) ratio
    gcVirTemp = gcSubsetProp(sP=sP,select=sgSelect,/virTemp,allTR=allTR,parNorm=parNorm)
    
    gcTemp.gal  = 10.0^gcTemp.gal / gcVirTemp.gal
    gcTemp.gmem = 10.0^gcTemp.gmem / gcVirTemp.gmem
  endif
  
  if ~allTR then $
    gcMass = gcSubsetProp(sP=sP,select=sgSelect,/curSingleVal,singleValField='mass',allTR=allTR)
  
  ; calculate masses of parents (for mass binning halos only)
  parentMass = gcSubsetProp(sP=sP,select=sgSelect,/parMass,allTR=allTR)

  ; get accretion time/tvir
  if accTime or accTvir or tVirAccNorm then begin
    mt = mergerTreeSubset(sP=sP)
    at = accretionTimes(sP=sP)
    
    gal_w  = where(at.AccTime_gal ne -1,count_gal)
    gmem_w = where(at.AccTime_gmem ne -1,count_gmem)

    ; accTime/accTvir: place the requested quantity into "gcTemp" for consistency
    if accTime then gcTemp = {gal :1/at.AccTime_gal[gal_w]-1    ,$
                              gmem:1/at.AccTime_gmem[gmem_w]-1   } ;convert scale factors -> redshift
    if accTvir then gcTemp = {gal :at.AccHaloTvir_gal[gal_w]    ,$
                              gmem:at.AccHaloTvir_gmem[gmem_w]   }
    
    ; create a galcatSub which indexes the galcat for the subset of tracked halos
    if sP.trMCPerCell eq 0 then begin
      galcatSub = galcatINDList(sP=sP,gcIDList=mt.galcatIDList) ;identical to mt.galcatSub
    endif else begin
      ; load child counts for all members of galcat
      maxt = maxTemps(sP=sP,/loadAllTRGal)
      child_counts_gal = maxt.child_counts
      maxt = maxTemps(sP=sP,/loadAllTRGmem)
      child_counts_gmem = maxt.child_counts
      maxt = !NULL
      
      ; get indices into subset of galcat found within accretionTimes()
      galcatSub = galCatINDList(sP=sP,gcIDList=mt.galcatIDList,$
                                 child_counts={gal:child_counts_gal,gmem:child_counts_gmem}) 
    endelse

    hMass = mt.hMass[0,*] ;temporary, using this or parentMass?
    ;mt     = !NULL
    
    ; tVirAccNorm: reform gcTemp as the merger tree subset and take ratio
    if tVirAccNorm then $
      gcTemp   = { gal  : gcTemp.gal[galcatSub.gal[gal_w]] / at.AccHaloTvir_gal[gal_w]      ,$
                   gmem : gcTemp.gmem[galcatSub.gmem[gmem_w]] / at.AccHaloTvir_gmem[gmem_w]  }

    ; all: reform gcRad,gcMass and parentMass as the merger tree subsets
    gcRad      = { gal  : gcRad.gal[galcatSub.gal[gal_w]]         ,$
                   gmem : gcRad.gmem[galcatSub.gmem[gmem_w]]       }
    parentMass = { gal  : parentMass.gal[galcatSub.gal[gal_w]]    ,$
                   gmem : parentMass.gmem[galcatSub.gmem[gmem_w]]  } ;TODO differs from mt.hMass(t=0) ??
    if ~allTR then $
      gcMass   = { gal  : gcMass.gal[galcatSub.gal[gal_w]]        ,$
                   gmem : gcMass.gmem[galcatSub.gmem[gmem_w]]      }
  endif

  for j=0,n_elements(massBins)-2 do begin
  
    ; plot config
    xrange = alog10([0.01,3.0])
    yrange = [4.0,7.0]
  
    binSizeRad  = 0.04 / (sP.res/128) ;0.04
    binSizeTemp = 0.05 / (sP.res/128) ;0.04
    
    ; preserve number of bins in log(rad) histogram and if changing yrange
    nBinsRad_linear = ceil((10.0^xrange[1]-10.0^xrange[0])/binSizeRad)+1
    nBinsTemp = ((yrange[1]-yrange[0])/binSizeTemp)+1
    binSizeRad_log  = (xrange[1]-xrange[0])/(nBinsRad_linear-1)
    
    if tVirNorm or tVirAccNorm then begin
      yrange = [0.0,1.4]
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
    if ~allTR then begin
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
      ytitle   = strmid(ytitle,5,strlen(ytitle)-6) + " / T"+textoidl("_{vir,cur}")
    endif
    if tVirAccNorm then begin
      plotBase = strmid(plotBase,0,4)+"_tvirAccNorm_"+sgSelect+'_'+parNorm
      ytitle   = strmid(ytitle,5,strlen(ytitle)-6) + " / T"+textoidl("_{vir,acc}")
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
    plotName = sP.plotPath + plotBase + '_slices.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    
    if n_elements(massBins) gt 2 then $
      plotName = strmid(plotName,0,strlen(plotName)-4) + '.mbin='+str(j)+'.eps'
    
    ;xrange = yrange
    xtitle = ytitle
    ;yrange = [0.0,1.05]
    
    ; plot vertical slices
    plotVerticalSlices,rad_gal,rad_gmem,temp_gal,temp_gmem,plotName,binSizeTemp*2.0,yrange,xtitle
    
  endfor ;j
  
end
