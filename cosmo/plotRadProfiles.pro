; plotRadProfiles.pro
; gas accretion project - plots of gas quantities vs radius (stacked in mass bins)
; dnelson sep.2012

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
       /vertical,/right,title=textoidl("M_{gas} [_{ }log M_{sun }]"),divisions=ndivs,ticknames=ticknames,ncolors=255
            
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
       /vertical,/right,title=textoidl("M_{gas} [_{ }log M_{sun }]"),divisions=ndivs,ticknames=ticknames,ncolors=255
         
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
       /vertical,/right,title=textoidl("M_{gas} [_{ }log M_{sun }]"),divisions=ndivs,ticknames=ticknames,ncolors=255
         
      if keyword_set(virTempRange) then begin
        fsc_text,xrange[1]*0.45,(yrange[1]-yrange[0])*0.92+yrange[0],massBinStr,alignment=0.5,color=fsc_color('black')
        fsc_text,xrange[0]*1.6,virTempRange[1]*1.02,"T"+textoidl("_{vir}"),alignment=0.5,color=fsc_color('yellow')
        fsc_plot,xrange,[virTempRange[0],virTempRange[0]],line=0,color=fsc_color('orange'),/overplot
        fsc_plot,xrange,[virTempRange[1],virTempRange[1]],line=0,color=fsc_color('orange'),/overplot
      endif
      
  end_PS
  !except = 1
end

; plot2DHisto(): plot 2d histogram of gas temperature (or other quantity) as a function of r_gas/r_vir
;                as well as vertical slices

pro plot2DHisto, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; config
  ;massBins = [0.0,1000.0] ; no massbins
  massBins = [9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5] ; log(M)
  
  ; select one of:
  curTemp      = 0 ; current temperature
  maxPastTemp  = 0 ; maximum previous temperature
  maxTempTime  = 0 ; time of maximum previous temperature
  accTime      = 0 ; redshift of accretion across virial radius of parent halo
  accTvir      = 0 ; virial temperature of parent halo at time of accretion
  
  curSingleVal = 1 ; current single gas value
    singleValField = 'coolingrate'
  
   ; if maxPastTemp or maxTempTime, choose one of: population min,mean,max to one statistic per
   ; gas cell instead of all the tracers individually
   trPopMax   = 0
   trPopMin   = 0 
   trPopMean  = 0
   
  tVirNorm    = 0  ; normalize temperature by virial temp of parent halo at the starting time
  tVirAccNorm = 0  ; normalize temperature by virial temp of parent halo at the time of accretion
   
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
  if total(curTemp+maxPastTemp+maxTempTime+accTime+accTvir+curSingleVal) ne 1 then message,'Just one'
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
                        curSingleVal=curSingleVal,singleValField=singleValField,$
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
  
  ; log coolingrate
  if curSingleVal then begin
    w = where(gcTemp.gal ne 0)
    gcTemp.gal[w] = alog10( gcTemp.gal[w] )
    w = where(gcTemp.gmem ne 0)
    gcTemp.gmem[w] = alog10( gcTemp.gmem[w] )
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
  
    binSizeRad  = 0.0035 ;0.014 / (sP.res/128) ;0.04
    binSizeTemp = 0.0125 ;0.05 / (sP.res/128) ;0.04
    
    ; preserve number of bins in log(rad) histogram and if changing yrange
    nBinsRad_linear = ceil((10.0^xrange[1]-10.0^xrange[0])/binSizeRad)+1
    nBinsTemp = ((yrange[1]-yrange[0])/binSizeTemp)+1
    binSizeRad_log  = (xrange[1]-xrange[0])/(nBinsRad_linear-1)
    
    if tVirNorm or tVirAccNorm then begin
      yrange = [-1.5,1.0]
      binSizeTemp = (yrange[1]-yrange[0])/(nBinsTemp-1)
    endif
    
    if accTime then begin
      yrange = [sP.redshift,4.0]
      binSizeTemp = (yrange[1]-yrange[0])/(nBinsTemp-1)
    endif
    
    if curSingleVal then begin
      yrange = [3.5,8.5]
      binSizeTemp = 0.022
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
    
    ; do mass weighting and 2D histogram
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
    
    if curSingleVal then begin
      plotBase = "coolrate_"+sgSelect+"_"+parNorm
      ytitle   = "Cooling Rate"
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

; binRadProfiles():

function binRadProfiles, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  ;massBins = [0.0,1000.0] ; no massbins
  massBins = [9.0,9.5,10.0,10.5,11.0,11.5,12.0] ; log(M)

  parNorm  = 'pri' ; pri,sec (normalize r_vir and temp (if doing tVirNorm) by primary or secondary parent)
  sgSelect = 'pri' ; pri,sec,all (subhalo selection config, only plot gas in this category of subhalos)
  
  ; consider only the subset with recorded accretion times for certain types of plots
  mergerTreeSubset    = 0
  accretionTimeSubset = 0
  
  xrange = [0.01,1.0]
  binSizeRad  = 0.03 ;0.014 / (sP.res/128) ;0.04
  nRadBins = (xrange[1]-xrange[0]) / binSizeRad
  radBinCen = findgen(nRadBins) / (nRadBins) * (xrange[1]-xrange[0]) + xrange[0] + binSizeRad*0.5
  
  ; check if save exists
  saveFilename = sP.derivPath + 'binRad.' + sP.saveTag + '.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.mb' + str(n_elements(massBins)) + '.pn-' + parNorm + '.' + sgSelect + '.' + $
    '.mts-' + str(mergerTreeSubset) + '.ats-' + str(accretionTimeSubset) + '.sav'
  
  ; results exist, return
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif 
  
  ; load group catalog just for counts of objects in each mass bin
  gc = loadGroupCat(sP=sP,/skipIDs)
  subgroupMasses = codeMassToLogMsun(gc.subgroupMass)
  gc = !NULL
  
  ; get normalized r/r_vir
  gcRad = gcSubsetProp(sP=sP,select=sgSelect,/rVirNorm,parNorm=parNorm,$
                       mergerTreeSubset=mergerTreeSubset,accretionTimeSubset=accretionTimeSubset)

  ; get gas properties
  gcCoolRate = gcSubsetProp(sP=sP,select=sgSelect,/curSingleVal,singleValField='coolingrate',$
                         mergerTreeSubset=mergerTreeSubset,accretionTimeSubset=accretionTimeSubset)
  gcTemp     = gcSubsetProp(sP=sP,select=sgSelect,/curTemp,$
                         mergerTreeSubset=mergerTreeSubset,accretionTimeSubset=accretionTimeSubset)               
  gcVirTemp  = gcSubsetProp(sP=sP,select=sgSelect,/virTemp,parNorm=parNorm,$
                         mergerTreeSubset=mergerTreeSubset,accretionTimeSubset=accretionTimeSubset)
    
  ; log (cooling rate)
  w = where(gcCoolRate.gal ne 0)
  gcCoolRate.gal[w] = alog10( gcCoolRate.gal[w] )
  w = where(gcCoolRate.gmem ne 0)
  gcCoolRate.gmem[w] = alog10( gcCoolRate.gmem[w] )

  ; temp/tvir
  gcTemp.gal  = 10.0^gcTemp.gal / 10.0^gcVirTemp.gal
  gcTemp.gmem = 10.0^gcTemp.gmem / 10.0^gcVirTemp.gmem
  
  ; save structure
  r = { gal  : { radCoolRateMB  : fltarr(n_elements(massBins)-1,nRadBins) + !values.f_nan  ,$
                 radCoolRateAll : fltarr(nRadBins) + !values.f_nan                         ,$
                 radTempMB      : fltarr(n_elements(massBins)-1,nRadBins) + !values.f_nan  ,$
                 radTempAll     : fltarr(nRadBins) + !values.f_nan                          } ,$
        gmem : { radCoolRateMB  : fltarr(n_elements(massBins)-1,nRadBins) + !values.f_nan  ,$
                 radCoolRateAll : fltarr(nRadBins) + !values.f_nan                         ,$
                 radTempMB      : fltarr(n_elements(massBins)-1,nRadBins) + !values.f_nan  ,$
                 radTempAll     : fltarr(nRadBins) + !values.f_nan                          } ,$     
        massBins : massBins, xrange : xrange, binSizeRad : binSizeRad, nRadBins : nRadBins ,$
        sP : sP, radBinCen : radBinCen }
  
  ; calculate masses of parents (for mass binning halos only)
  parentMass = gcSubsetProp(sP=sP,select=sgSelect,/parMass,$
                            mergerTreeSubset=mergerTreeSubset,accretionTimeSubset=accretionTimeSubset)

  for j=0,n_elements(massBins)-2 do begin

    ; select members of this parent mass bins and r>0<inf
    wGal  = where(parentMass.gal gt massBins[j] and parentMass.gal le massBins[j+1] and $
                  gcRad.gal gt 0.0 and finite(gcRad.gal),count1)
    wGmem = where(parentMass.gmem gt massBins[j] and parentMass.gmem le massBins[j+1] and $
                  gcRad.gmem gt 0.0 and finite(gcRad.gmem),count2)
    
    ;print,j,count1,count2
    
    ; select all particles/tracers in this mass bin
    coolrate_gal  = gcCoolRate.gal[wGal]
    coolrate_gmem = gcCoolRate.gmem[wGmem]
    temp_gal   = gcTemp.gal[wGal]
    temp_gmem  = gcTemp.gmem[wGmem]
    rad_gal    = alog10( gcRad.gal[wGal] )
    rad_gmem   = alog10( gcRad.gmem[wGmem] )

    if count1 eq 0 or count2 eq 0 then continue ; no halos in this mass bin
    
    ; do median profiles
    for i=0,nRadBins-1 do begin
      ; bin bounds and selection
      binMin = alog10(xrange[0] + i*binSizeRad)
      binMax = alog10(xrange[0] + (i+1)*binSizeRad)
      
    ; gal
    w = where(rad_gal ge binMin and rad_gal lt binMax,count)
    if count gt 0 then begin
      r.gal.radCoolRateMB[j,i] = median(coolrate_gal[w])
      r.gal.radTempMB[j,i]     = median(temp_gal[w])
    endif
    
    ; gmem
    w = where(rad_gmem ge binMin and rad_gmem lt binMax,count)
    if count gt 0 then begin
      r.gmem.radCoolRateMB[j,i] = median(coolrate_gmem[w])
      r.gmem.radTempMB[j,i]     = median(temp_gmem[w])
    endif
        
      ;print,10^binMin,10^binMax,count1,count2
    endfor

  endfor
  
  ; do non-mass bin
  wGal  = where(gcRad.gal gt 0.0 and finite(gcRad.gal),count1)
  wGmem = where(gcRad.gmem gt 0.0 and finite(gcRad.gmem),count2)
  coolrate_gal  = gcCoolRate.gal[wGal]
  coolrate_gmem = gcCoolRate.gmem[wGmem]
  temp_gal   = gcTemp.gal[wGal]
  temp_gmem  = gcTemp.gmem[wGmem]
  rad_gal    = alog10( gcRad.gal[wGal] )
  rad_gmem   = alog10( gcRad.gmem[wGmem] )
    
  for i=0,nRadBins-1 do begin
    ; bin bounds and selection
    binMin = alog10(xrange[0] + i*binSizeRad)
    binMax = alog10(xrange[0] + (i+1)*binSizeRad)
    
    ; gal
    w = where(rad_gal ge binMin and rad_gal lt binMax,count)
    if count gt 0 then begin
      r.gal.radCoolRateAll[i] = median(coolrate_gal[w])
      r.gal.radTempAll[i]     = median(temp_gal[w])
    endif
    
    ; gmem
    w = where(rad_gmem ge binMin and rad_gmem lt binMax,count)
    if count gt 0 then begin
      r.gmem.radCoolRateAll[i] = median(coolrate_gmem[w])
      r.gmem.radTempAll[i]     = median(temp_gmem[w])
    endif
  endfor
    
  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  return, r
  
end

; plotRadProfile():

pro plotRadProfile
  
  ; config
  redshift = 2.0

  lines  = [1,0,2] ; 128,512,256
  cInd   = 1 ; color index

  sPg = simParams(res=128,run='gadget',redshift=redshift)
  sPa = simParams(res=128,run='arepo',redshift=redshift)

  ga_128 = binRadProfiles(sP=sPg)
  ar_128 = binRadProfiles(sP=sPa)
  
  ga_256 = binRadProfiles(sP=simParams(res=256,run='gadget',redshift=redshift))
  ar_256 = binRadProfiles(sP=simParams(res=256,run='arepo',redshift=redshift))
  
  ga_512 = binRadProfiles(sP=simParams(res=512,run='gadget',redshift=redshift))
  ar_512 = binRadProfiles(sP=simParams(res=512,run='arepo',redshift=redshift))

  for j=0,n_elements(ga_128.massBins)-2 do begin

    yrange = [3.0,6.0]

    ; plot (0)
    start_PS, sPg.plotPath + 'coolrate.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
              str(sPg.snap)+'_mbin-'+str(j)+'.eps'
  
      ; extra config for mass bins
      massBinStr   = !NULL
      virTempRange = !NULL
        
      massBinStr = string(ga_128.massBins[j],format='(f4.1)') + '-' + $
                   string(ga_128.massBins[j+1],format='(f4.1)')
      
      massRangeCode = 10.0^[ga_128.massBins[j],ga_128.massBins[j+1]] / 1e10
      virTempRange = alog10( codeMassToVirTemp(massRangeCode,sP=sPg) )
  
      ; plot 1d profile
      cgPlot,[0],[0],/nodata,xrange=[0.1,1.0],yrange=yrange,title=massBinStr,$
        xtitle=textoidl('r / r_{vir}'),ytitle="Log ( Cooling Rate )",/xlog,/xs,/ys,xminor=0
        
      cgPlot,[0.15,0.15],yrange,line=1,color=cgColor('light gray'),/overplot
      
      cgPlot,ga_128.radBinCen,ga_128.gmem.radCoolRateMB[j,*],color=sPg.colorsG[cInd],line=lines[0],/overplot
      cgPlot,ga_256.radBinCen,ga_256.gmem.radCoolRateMB[j,*],color=sPg.colorsG[cInd],line=lines[2],/overplot
      cgPlot,ga_512.radBinCen,ga_512.gmem.radCoolRateMB[j,*],color=sPg.colorsG[cInd],line=lines[1],/overplot
        
      cgPlot,ar_128.radBinCen,ar_128.gmem.radCoolRateMB[j,*],color=sPa.colorsA[cInd],line=lines[0],/overplot
      cgPlot,ar_256.radBinCen,ar_256.gmem.radCoolRateMB[j,*],color=sPa.colorsA[cInd],line=lines[2],/overplot
      cgPlot,ar_512.radBinCen,ar_512.gmem.radCoolRateMB[j,*],color=sPa.colorsA[cInd],line=lines[1],/overplot
        
      ; labels
      legend,textoidl(['128^3','256^3','512^3']),linestyle=[lines[0],lines[2],lines[1]],$
        box=0,/top,/right,linesize=0.25
      legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/bottom,/left
  
    end_PS
    
    yrange = [0.1,2.0]
    
    ; plot (1)
    start_PS, sPg.plotPath + 'radtemp.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+$
              str(sPg.snap)+'_mbin-'+str(j)+'.eps'
  
      ; extra config for mass bins
      massBinStr   = !NULL
      virTempRange = !NULL
        
      massBinStr = string(ga_128.massBins[j],format='(f4.1)') + '-' + $
                   string(ga_128.massBins[j+1],format='(f4.1)')
      
      massRangeCode = 10.0^[ga_128.massBins[j],ga_128.massBins[j+1]] / 1e10
      virTempRange = alog10( codeMassToVirTemp(massRangeCode,sP=sPg) )
  
      ; plot 1d profile
      cgPlot,[0],[0],/nodata,xrange=[0.1,1.0],yrange=yrange,title=massBinStr,$
        xtitle=textoidl('r / r_{vir}'),ytitle=textoidl("T / T_{vir}"),/xlog,/xs,/ys,xminor=0
        
      cgPlot,[0.15,0.15],yrange,line=1,color=cgColor('light gray'),/overplot
      
      cgPlot,ga_128.radBinCen,ga_128.gmem.radTempMB[j,*],color=sPg.colorsG[cInd],line=lines[0],/overplot
      cgPlot,ga_256.radBinCen,ga_256.gmem.radTempMB[j,*],color=sPg.colorsG[cInd],line=lines[2],/overplot
      cgPlot,ga_512.radBinCen,ga_512.gmem.radTempMB[j,*],color=sPg.colorsG[cInd],line=lines[1],/overplot
        
      cgPlot,ar_128.radBinCen,ar_128.gmem.radTempMB[j,*],color=sPa.colorsA[cInd],line=lines[0],/overplot
      cgPlot,ar_256.radBinCen,ar_256.gmem.radTempMB[j,*],color=sPa.colorsA[cInd],line=lines[2],/overplot
      cgPlot,ar_512.radBinCen,ar_512.gmem.radTempMB[j,*],color=sPa.colorsA[cInd],line=lines[1],/overplot
        
      ; labels
      legend,textoidl(['128^3','256^3','512^3']),linestyle=[lines[0],lines[2],lines[1]],$
        box=0,/top,/right,linesize=0.25
      legend,['gadget','arepo'],textcolors=[sPg.colorsG[cInd],sPa.colorsA[cInd]],box=0,/bottom,/left
  
    end_PS
  
  endfor
  
  stop

end

