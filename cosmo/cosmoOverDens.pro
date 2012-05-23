; cosmoOverDens.pro
; cosmological boxes - selection and analysis based on DM field overdensity criteria
; dnelson apr.2012

; findDMOverDensAtGas(): foreach gas particle/cell estimate the local dark matter overdensity

function findDMOverDensAtGas, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; check for save file
  saveFilename = sP.derivPath + 'gasDMdens_' + sP.savPrefix + str(sP.res) + '_' + str(sP.snap) + '.sav'
  
  if file_test(saveFilename) then begin
    restore, saveFilename
  endif else begin
    ; load
    gas_pos = loadSnapshotSubset(sP=sP, field='pos', partType='gas')
    dm_pos  = loadSnapshotSubset(sP=sP, field='pos', partType='dm')
    
    ; EXTERNAL C - calculate nearest dm neighbor to each gas position
    nearest_dm_ind = calcNN(dm_pos,gas_pos,boxSize=sP.boxSize,ndims=3)

    ; DEBUG: for a random subset, verify nearest dm neighbors to each gas position
    nVerify = 100
    indVerify = floor(randomu(seed,nVerify) * n_elements(gas_pos[0,*]))

    for i=0,nVerify-1 do begin
      ; find minimum distance
      dists = periodicDists(gas_pos[*,indVerify[i]],dm_pos,sP=sP)
      w = where(dists eq min(dists),count)
      
      ; fail?
      if (count ne 1) then message,'ERROR: More than one mindist in NN debug search!'
      if (w[0] ne nearest_dm_ind[indVerify[i]]) then message,'ERROR: Nearest neighbor verification failed'
    endfor

    gas_pos = !NULL
    dm_pos  = !NULL
    
    ; load saved DM densities and associate with gas
    dens_dm = loadHsmlDir(sP=sP,partType='dm',/readDens)
    gas_dmDens = dens_dm[nearest_dm_ind]
    
    ; convert density to log(overdensity)
    h = loadSnapshotHeader(sP=sP)
    meanDensBox = float(h.nPartTot[partTypeNum('dm')] * h.masstable[partTypeNum('dm')] / (sP.boxSize)^3.0)
    gas_dmDens = alog10( gas_dmDens / meanDensBox )    
    
    ; save
    save,gas_dmDens,filename=saveFilename
    print,'Saved: ',saveFilename
  endelse
  
  return,gas_dmDens

end

; plotDMComps(): compare distances and densities between identical ID dm particles in gadget/arepo boxes
;                also compare global density distributions in the whole box

pro plotDMComps, redshift=redshift

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if n_elements(redshift) eq 0 then message,'Specify redshift.'

  ; config
  binsize = 1.0
  binsizeLog = 0.1
  res = 512
  
  sPa = simParams(res=res,run='arepo',redshift=redshift)
  sPg = simParams(res=res,run='gadget',redshift=redshift)
  
  ; load arepo
  ha = loadSnapshotHeader(sP=sPa)
  pos_ar  = loadSnapshotSubset(sP=sPa,partType='dm',field='pos')
  dens_ar = loadHsmlDir(sP=sPa,partType='dm',/readDens)
  ids_ar  = loadSnapshotSubset(sP=sPa,partType='dm',field='ids')
  
  ; load gadget
  hg = loadSnapshotHeader(sP=sPg)
  pos_ga  = loadSnapshotSubset(sP=sPg,partType='dm',field='pos')
  dens_ga = loadHsmlDir(sP=sPg,partType='dm',/readDens)
  ids_ga  = loadSnapshotSubset(sP=sPg,partType='dm',field='ids')
  
  ; verify all IDs matched
  match,ids_ar,ids_ga,ind1,ind2,count=count
  if count ne n_elements(ids_ar) or count ne n_elements(ids_ga) then message,'Error failed match.'
  
  ; sort positions by ID
  sort_ar = sort(ids_ar)
  pos_ar[0,*] = pos_ar[0,sort_ar]
  pos_ar[1,*] = pos_ar[1,sort_ar]
  pos_ar[2,*] = pos_ar[2,sort_ar]
  dens_ar     = dens_ar[sort_ar]
  sort_ar = !NULL
  ids_ar  = !NULL
  
  sort_ga = sort(ids_ga)
  pos_ga[0,*] = pos_ga[0,sort_ga]
  pos_ga[1,*] = pos_ga[1,sort_ga]
  pos_ga[2,*] = pos_ga[2,sort_ga]
  dens_ga     = dens_ga[sort_ga]
  sort_ga = !NULL
  ids_ga  = !NULL
  
  ; histogram particle by particle DM density field ratio
  dm_ratio = alog10(dens_ar / dens_ga)
  hist_ratio = histogram(dm_ratio,binsize=binsizeLog,loc=loc_ratio)
  hist_ratio /= float(total(hist_ratio))
  
  ; convert densities to log(overdensities)
  meanDensBox = float(ha.nPartTot[partTypeNum('dm')] * ha.masstable[partTypeNum('dm')] / (sPa.boxSize)^3.0)
  dens_ar = alog10( dens_ar / meanDensBox )
  dens_ga = alog10( dens_ga / meanDensBox )  
  
  ; calculate distances
  dists = periodicDists(pos_ar,pos_ga,sP=sPa)
  
  ; histogram distances distribution
  hist = histogram(dists,binsize=binSize,loc=loc)
  
  ; plot (1)
  start_PS, sPa.plotPath + 'dmpos.comp.'+str(sPa.res)+'_'+str(sPa.snap)+'.eps'
    xrange = [1.0,1000.0]
    yrange = [10,max(hist)*1.5]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/xlog,/ylog,$
      ytitle="Count",xtitle=textoidl("DM Distance Offset [ckpc]"),$
      title=str(sPa.res)+textoidl("^3")+" Gadget vs. Arepo (z="+string(sPa.redshift,format='(f3.1)')+")"

    ; mark maximum
    w = where(hist eq max(hist))
    cgPlot,[loc[w],loc[w]]+binsize*0.5,yrange,line=0,color=cgColor('light gray'),/overplot

    ; histograms    
    cgPlot,loc+binsize*0.5,hist,line=0,color=getColor(0),/overplot
  
  end_PS
  
  ; histogram density distributions and normalize by total
  hist_ar = histogram(dens_ar,binsize=binSizeLog,loc=loc_ar)
  hist_ga = histogram(dens_ga,binsize=binSizeLog,loc=loc_ga)
  
  hist_ar /= float(total(hist_ar))
  hist_ga /= float(total(hist_ga))
  
  ; plot (2)
  start_PS, sPa.plotPath + 'dmoverdens.comp.'+str(sPa.res)+'_'+str(sPa.snap)+'.eps'
    xrange = [-2.0,9.0]
    yrange = [1e-5,1e-1]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="Mass Fraction",xtitle=textoidl("log ( \rho_{DM} / <\rho_{DM}> )")+"",$
      title=str(sPa.res)+textoidl("^3")+" Gadget vs. Arepo (z="+string(sPa.redshift,format='(f3.1)')+")"
    cgPlot,[0,0],yrange,line=0,color=fsc_color('light gray'),/overplot

    ; histograms    
    cgPlot,loc_ar+binsizeLog*0.5,hist_ar,line=0,color=getColor(3),/overplot
    cgPlot,loc_ga+binsizeLog*0.5,hist_ga,line=0,color=getColor(1),/overplot
  
    ; legend
    legend,['gadget','arepo'],textcolors=getColor([1,3],/name),box=0,/right,/top
  
  end_PS

  ; plot (3)
  start_PS, sPa.plotPath + 'dmoverdens.pbyp.comp.'+str(sPa.res)+'_'+str(sPa.snap)+'.eps'
    xrange = [-5.0,5.0]
    yrange = [1e-5,5e-1]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="Mass Fraction",xtitle=textoidl("log ( \rho_{DM,ar} / \rho_{DM,ga} )"),$
      title=str(sPa.res)+textoidl("^3")+" Gadget vs. Arepo (z="+string(sPa.redshift,format='(f3.1)')+")"
    cgPlot,[0,0],yrange,line=0,color=fsc_color('light gray'),/overplot

    ; histograms    
    cgPlot,loc_ratio+binsizeLog*0.5,hist_ratio,line=0,color=getColor(0),/overplot
  
  end_PS
end

; plotGasPropsVsDMOverdens(): plot statistics of gas properties as a function of DM overdensity

pro plotGasPropsVsDMOverdens

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  redshift = 2.0
  res = 512
  xrange = [-1.0,5.0] ; log (dm overdens)
  
  dmBinsizeLog = 0.2 / (res/128) 
  sPa = simParams(res=res,run='arepo',redshift=redshift)
  sPg = simParams(res=res,run='gadget',redshift=redshift)  
  
  ; select a spatial subset instead of processing the whole box?
  sizeFac = 3.5  ; times rvir
  hID_ga  = 4422 ;(z2.198) ;4731 (z2.17) ;2289 (z2.301)
  hID_ar  = 3844 ;(z2.198) ;4296 (z2.17) ;2034 (z2.301)
  
  if hID_ga ne 0 then subsetTag = "hg"+str(hID_ga)+"_ha"+str(hID_ar)+"_" else subsetTag = "global_"
  if hID_ga ne 0 then dmBinsizeLog *= 4.0 ; increase binsize if subsetting
  
  ; dm overdensity bins
  dmNbins  = floor((xrange[1]-xrange[0]) / dmBinsizeLog)
  dmBins   = linspace(xrange[0],xrange[1],dmNbins+1) ; edges
  dmBinCen = linspace(xrange[0],xrange[1],dmNbins+1) + dmBinsizeLog/2.0
  dmBinCen = dmBinCen[0:-2] ; remove last
  
  saveFilename = sPa.derivPath + 'odhists_'+subsetTag+str(dmNbins)+'.sav'
  
  if file_test(saveFilename) then begin
    restore,saveFilename
  endif else begin
  
    ; select spatial subset if requested
    if hID_ga ne 0 then begin
      print,'subselection...'
      ; load subgroup positions and group catalog for halo position and rvir
      sgcen = subgroupPosByMostBoundID(sP=sPg) 
      gc = loadGroupCat(sP=sPg,/skipIDs)
      hPos_ga  = sgcen[*,hID_ga]
      hRVir_ga = gc.group_r_crit200[gc.subgroupGrNr[hID_ga]] ;ckpc
      
      sgcen = subgroupPosByMostBoundID(sP=sPa)
      gc = loadGroupCat(sP=sPa,/skipIDs)
      hPos_ar  = sgcen[*,hID_ar]
      hRVir_ar = gc.group_r_crit200[gc.subgroupGrNr[hID_ar]] ;ckpc
      
      ; load gas positions and make selection
      pos = loadSnapshotSubset(sP=sPg, field='pos', partType='gas')
      dists = periodicDists(hPos_ga,pos,sP=sPg)
      wGas_ga = where(dists lt sizeFac*mean([hRVir_ga,hRVir_ar]),count)
      
      pos = loadSnapshotSubset(sP=sPa, field='pos', partType='gas')
      dists = periodicDists(hPos_ar,pos,sP=sPa)
      wGas_ar = where(dists lt sizeFac*mean([hRVir_ga,hRVir_ar]),count)
      
      pos   = !NULL
      dists = !NULL
    endif
  
    ; histograms
    histDens_ga = {mean:fltarr(dmNBins),median:fltarr(dmNBins),stddev:fltarr(dmNBins),$
                   min:fltarr(dmNBins),max:fltarr(dmNBins),count:fltarr(dmNBins)}
    histDens_ar = {mean:fltarr(dmNBins),median:fltarr(dmNBins),stddev:fltarr(dmNBins),$
                   min:fltarr(dmNBins),max:fltarr(dmNBins),count:fltarr(dmNBins)}
    histEnt_ga  = {mean:fltarr(dmNBins),median:fltarr(dmNBins),stddev:fltarr(dmNBins),$
                   min:fltarr(dmNBins),max:fltarr(dmNBins),count:fltarr(dmNBins)}
    histEnt_ar  = {mean:fltarr(dmNBins),median:fltarr(dmNBins),stddev:fltarr(dmNBins),$
                   min:fltarr(dmNBins),max:fltarr(dmNBins),count:fltarr(dmNBins)}
    histTemp_ga = {mean:fltarr(dmNBins),median:fltarr(dmNBins),stddev:fltarr(dmNBins),$
                   min:fltarr(dmNBins),max:fltarr(dmNBins),count:fltarr(dmNBins)}
    histTemp_ar = {mean:fltarr(dmNBins),median:fltarr(dmNBins),stddev:fltarr(dmNBins),$
                   min:fltarr(dmNBins),max:fltarr(dmNBins),count:fltarr(dmNBins)}
    histVel_ga  = {mean:fltarr(dmNBins),median:fltarr(dmNBins),stddev:fltarr(dmNBins),$
                   min:fltarr(dmNBins),max:fltarr(dmNBins),count:fltarr(dmNBins)}
    histVel_ar  = {mean:fltarr(dmNBins),median:fltarr(dmNBins),stddev:fltarr(dmNBins),$
                   min:fltarr(dmNBins),max:fltarr(dmNBins),count:fltarr(dmNBins)}     
    
    ; load dark matter overdensities to make gas selection
    print,'gadget...'
    gasDMDens_ga = findDMOverDensAtGas(sP=sPg)
    if hID_ga ne 0 then gasDMDens_ga = gasDMDens_ga[wGas_ga]
   
    ; density: normalize to critical baryon density at this redshift (gadget)
    dens = loadSnapshotSubset(sP=sPg, field='dens', partType='gas')
    if hID_ga ne 0 then dens = dens[wGas_ga]
    gasProp_ga = alog10(rhoRatioToCrit(dens,redshift=sPg.redshift))
    
    for i=0,dmNbins-1 do begin
      w = where(gasDMDens_ga gt dmBins[i] and gasDMDens_ga le dmBins[i+1],count)
      if count gt 0 then begin
        histDens_ga.mean[i]   = mean(gasProp_ga[w])
        histDens_ga.median[i] = median(gasProp_ga[w])
        histDens_ga.stddev[i] = stddev(gasProp_ga[w])
        histDens_ga.min[i]    = min(gasProp_ga[w])
        histDens_ga.max[i]    = max(gasProp_ga[w])
        histDens_ga.count[i]  = count
      endif
    endfor
    
    ; entropy
    u = loadSnapshotSubset(sP=sPg,partType='gas',field='u')
    if hID_ga ne 0 then u = u[wGas_ga]
    gasProp_ga = calcEntropyCGS(u,dens,/log,sP=sPg)
    dens = !NULL
    
    for i=0,dmNbins-1 do begin
      w = where(gasDMDens_ga gt dmBins[i] and gasDMDens_ga le dmBins[i+1],count)
      if count gt 0 then begin
        histEnt_ga.mean[i]   = mean(gasProp_ga[w])
        histEnt_ga.median[i] = median(gasProp_ga[w])
        histEnt_ga.stddev[i] = stddev(gasProp_ga[w])
        histEnt_ga.min[i]    = min(gasProp_ga[w])
        histEnt_ga.max[i]    = max(gasProp_ga[w])
        histEnt_ga.count[i]  = count
      endif
    endfor
    
    ; temperature
    nelec = loadSnapshotSubset(sP=sPg,partType='gas',field='nelec')
    if hID_ga ne 0 then nelec = nelec[wGas_ga]
    gasProp_ga = convertUtoTemp(u,nelec,/log)
    u     = !NULL
    nelec = !NULL
    
    for i=0,dmNbins-1 do begin
      w = where(gasDMDens_ga gt dmBins[i] and gasDMDens_ga le dmBins[i+1],count)
      if count gt 0 then begin
        histTemp_ga.mean[i]   = mean(gasProp_ga[w])
        histTemp_ga.median[i] = median(gasProp_ga[w])
        histTemp_ga.stddev[i] = stddev(gasProp_ga[w])
        histTemp_ga.min[i]    = min(gasProp_ga[w])
        histTemp_ga.max[i]    = max(gasProp_ga[w])
        histTemp_ga.count[i]  = count
      endif
    endfor
    
    ; velocity (magnitude)
    gasProp_ga = loadSnapshotSubset(sP=sPg,partType='gas',field='vel')
    if hID_ga ne 0 then gasProp_ga = gasProp_ga[*,wGas_ga]
    gasProp_ga = sqrt(gasProp_ga[0,*]*gasProp_ga[0,*] + $
                      gasProp_ga[1,*]*gasProp_ga[1,*] + $
                      gasProp_ga[2,*]*gasProp_ga[2,*])
    
    for i=0,dmNbins-1 do begin
      w = where(gasDMDens_ga gt dmBins[i] and gasDMDens_ga le dmBins[i+1],count)
      if count gt 0 then begin
        histVel_ga.mean[i]   = mean(gasProp_ga[w])
        histVel_ga.median[i] = median(gasProp_ga[w])
        histVel_ga.stddev[i] = stddev(gasProp_ga[w])
        histVel_ga.min[i]    = min(gasProp_ga[w])
        histVel_ga.max[i]    = max(gasProp_ga[w])
        histVel_ga.count[i]  = count
      endif
    endfor
      
    ; repeat it all for arepo
    print,'arepo...'
    gasProp_ga   = !NULL
    gasDMDens_ga = !NULL
    
    gasDMDens_ar = findDMOverDensAtGas(sP=sPa)      
    if hID_ar ne 0 then gasDMDens_ar = gasDMDens_ar[wGas_ar]
    
    ; density: normalize to critical baryon density at this redshift (gadget)
    dens = loadSnapshotSubset(sP=sPa, field='dens', partType='gas')
    if hID_ar ne 0 then dens = dens[wGas_ar]
    gasProp_ar = alog10(rhoRatioToCrit(dens,redshift=sPa.redshift))
    
    for i=0,dmNbins-1 do begin
      w = where(gasDMDens_ar gt dmBins[i] and gasDMDens_ar le dmBins[i+1],count)
      if count gt 0 then begin
        histDens_ar.mean[i]   = mean(gasProp_ar[w])
        histDens_ar.median[i] = median(gasProp_ar[w])
        histDens_ar.stddev[i] = stddev(gasProp_ar[w])
        histDens_ar.min[i]    = min(gasProp_ar[w])
        histDens_ar.max[i]    = max(gasProp_ar[w])
        histDens_ar.count[i]  = count
      endif
    endfor
    
    ; entropy
    u = loadSnapshotSubset(sP=sPa,partType='gas',field='u')
    if hID_ar ne 0 then u = u[wGas_ar]
    gasProp_ar = calcEntropyCGS(u,dens,/log,sP=sPa)
    dens = !NULL
      
    for i=0,dmNbins-1 do begin
      w = where(gasDMDens_ar gt dmBins[i] and gasDMDens_ar le dmBins[i+1],count)
      if count gt 0 then begin
        histEnt_ar.mean[i]   = mean(gasProp_ar[w])
        histEnt_ar.median[i] = median(gasProp_ar[w])
        histEnt_ar.stddev[i] = stddev(gasProp_ar[w])
        histEnt_ar.min[i]    = min(gasProp_ar[w])
        histEnt_ar.max[i]    = max(gasProp_ar[w])
        histEnt_ar.count[i]  = count
      endif
    endfor
    
    ; temperature
    nelec = loadSnapshotSubset(sP=sPa,partType='gas',field='nelec')
    if hID_ar ne 0 then nelec = nelec[wGas_ar]
    gasProp_ar = convertUtoTemp(u,nelec,/log)
    u     = !NULL
    nelec = !NULL
    
    for i=0,dmNbins-1 do begin
      w = where(gasDMDens_ar gt dmBins[i] and gasDMDens_ar le dmBins[i+1],count)
      if count gt 0 then begin
        histTemp_ar.mean[i]   = mean(gasProp_ar[w])
        histTemp_ar.median[i] = median(gasProp_ar[w])
        histTemp_ar.stddev[i] = stddev(gasProp_ar[w])
        histTemp_ar.min[i]    = min(gasProp_ar[w])
        histTemp_ar.max[i]    = max(gasProp_ar[w])
        histTemp_ar.count[i]  = count
      endif
    endfor
    
    ; velocity
    gasProp_ar = loadSnapshotSubset(sP=sPa,partType='gas',field='vel')
    if hID_ar ne 0 then gasProp_ar = gasProp_ar[*,wGas_ar]
    gasProp_ar = sqrt(gasProp_ar[0,*]*gasProp_ar[0,*] + $
                      gasProp_ar[1,*]*gasProp_ar[1,*] + $
                      gasProp_ar[2,*]*gasProp_ar[2,*])
    
    for i=0,dmNbins-1 do begin
      w = where(gasDMDens_ar gt dmBins[i] and gasDMDens_ar le dmBins[i+1],count)
      if count gt 0 then begin
        histVel_ar.mean[i]   = mean(gasProp_ar[w])
        histVel_ar.median[i] = median(gasProp_ar[w])
        histVel_ar.stddev[i] = stddev(gasProp_ar[w])
        histVel_ar.min[i]    = min(gasProp_ar[w])
        histVel_ar.max[i]    = max(gasProp_ar[w])
        histVel_ar.count[i]  = count
      endif
    endfor
      
    gasProp_ar   = !NULL
    gasDMDens_ar = !NULL
  
    ; save the histograms
    config = {xrange:xrange,dmBinsizeLog:dmBinsizeLog,hID_ga:hID_ga,hID_ar:hID_ar,sPa:sPa,sPg:sPg}
    save,histDens_ga,histDens_ar,histEnt_ga,histEnt_ar,histTemp_ga,histTemp_ar,$
         histVel_ga,histVel_ar,config,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sPa.derivPath))
    
  endelse
  
  ; plot
  start_PS, sPa.plotPath + 'gasprops.'+subsetTag+str(sPa.res)+'_'+str(sPa.snap)+'.eps', /big
            
    xrange = minmax(dmBins)
    sigFac = 2
    
    !x.margin = [8.0,3.0]
    !y.margin = [4.0,1.0]
    
    x0 = 0.1 & x1 = 0.5 & x2 = 0.9
    y0 = 0.15 & y1 = 0.55 & y2 = 0.95
    pos = list( [x0,y1,x1,y2] ,$ ; ul
                [x1,y1,x2,y2] ,$ ; ur
                [x0,y0,x1,y1] ,$ ; ll
                [x1,y0,x2,y1] )  ; lr
    
    ; temperature
    yrange = [3.0,6.5]
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,pos=pos[0],$
      ytitle=textoidl("log ( T_{gas} ) [K]"),xtitle="",$
      title="",xtickname=replicate(' ',10)
      
    ; mark dm overdensity range
    cgPlot,[0.0,0.0],yrange,line=0,color=fsc_color('light gray'),/overplot
    cgPlot,[1.0,1.0],yrange,line=0,color=fsc_color('light gray'),/overplot

    cgPlot,dmbinCen,histTemp_ga.mean,line=1,color=getColor(1),/overplot
    cgPlot,dmbinCen,histTemp_ga.median,line=2,color=getColor(1),/overplot
    cgPlot,dmbinCen,histTemp_ga.mean+sigFac*histTemp_ga.stddev,line=0,color=getColor(1),/overplot
    cgPlot,dmbinCen,histTemp_ga.mean-sigFac*histTemp_ga.stddev,line=0,color=getColor(1),/overplot
    
    cgPlot,dmbinCen,histTemp_ar.mean,line=1,color=getColor(3),/overplot
    cgPlot,dmbinCen,histTemp_ar.median,line=2,color=getColor(3),/overplot
    cgPlot,dmbinCen,histTemp_ar.mean+sigFac*histTemp_ar.stddev,line=0,color=getColor(3),/overplot
    cgPlot,dmbinCen,histTemp_ar.mean-sigFac*histTemp_ar.stddev,line=0,color=getColor(3),/overplot
    
    ; legend
    legend,['gadget','arepo'],textcolors=getColor([1,3],/name),box=0,/top,/left
    
    ; density
    yrange = [-2.0,4.0]
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,ys=1,pos=pos[1],/noerase,$
      ytitle="",xtitle="",title="",xtickname=replicate(' ',10),ytickname=replicate(' ',10)
    cgAxis,/yaxis,yrange=yrange,ys=1
    cgText,x2+0.05,0.5*(y1+y2),textoidl("log ( \rho / \rho_{crit,b,z} )"),alignment=0.5,orientation=270.0,/normal
      
    ; mark dm overdensity range
    cgPlot,[0.0,0.0],yrange,line=0,color=fsc_color('light gray'),/overplot
    cgPlot,[1.0,1.0],yrange,line=0,color=fsc_color('light gray'),/overplot

    cgPlot,dmbinCen,histDens_ga.mean,line=1,color=getColor(1),/overplot
    cgPlot,dmbinCen,histDens_ga.median,line=2,color=getColor(1),/overplot
    cgPlot,dmbinCen,histDens_ga.mean+sigFac*histDens_ga.stddev,line=0,color=getColor(1),/overplot
    cgPlot,dmbinCen,histDens_ga.mean-sigFac*histDens_ga.stddev,line=0,color=getColor(1),/overplot
    
    cgPlot,dmbinCen,histDens_ar.mean,line=1,color=getColor(3),/overplot
    cgPlot,dmbinCen,histDens_ar.median,line=2,color=getColor(3),/overplot
    cgPlot,dmbinCen,histDens_ar.mean+sigFac*histDens_ar.stddev,line=0,color=getColor(3),/overplot
    cgPlot,dmbinCen,histDens_ar.mean-sigFac*histDens_ar.stddev,line=0,color=getColor(3),/overplot
    
    ; entropy
    yrange = [4.0,9.5]
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,pos=pos[2],/noerase,$
      ytitle=textoidl("log ( Entropy ) [K cm^2]"),xtitle="",$
      title="",xticks=6,xtickv=[-1,0,1,2,3,4]
      
    ; mark dm overdensity range
    cgPlot,[0.0,0.0],yrange,line=0,color=fsc_color('light gray'),/overplot
    cgPlot,[1.0,1.0],yrange,line=0,color=fsc_color('light gray'),/overplot

    cgPlot,dmbinCen,histEnt_ga.mean,line=1,color=getColor(1),/overplot
    cgPlot,dmbinCen,histEnt_ga.median,line=2,color=getColor(1),/overplot
    cgPlot,dmbinCen,histEnt_ga.mean+sigFac*histEnt_ga.stddev,line=0,color=getColor(1),/overplot
    cgPlot,dmbinCen,histEnt_ga.mean-sigFac*histEnt_ga.stddev,line=0,color=getColor(1),/overplot
    
    cgPlot,dmbinCen,histEnt_ar.mean,line=1,color=getColor(3),/overplot
    cgPlot,dmbinCen,histEnt_ar.median,line=2,color=getColor(3),/overplot
    cgPlot,dmbinCen,histEnt_ar.mean+sigFac*histEnt_ar.stddev,line=0,color=getColor(3),/overplot
    cgPlot,dmbinCen,histEnt_ar.mean-sigFac*histEnt_ar.stddev,line=0,color=getColor(3),/overplot
  
    ; velocity
    yrange = [50,450]
    sigFac = 1
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,ys=1,pos=pos[3],/noerase,$
      ytitle="",xtitle="",title="",ytickname=replicate(' ',10),$
      xticks=7,xtickv=[-1,0,1,2,3,4,5]
    cgAxis,/yaxis,yrange=yrange,ys=1
    cgText,x2+0.07,0.5*(y0+y1),textoidl("Velocity [km/s]"),alignment=0.5,orientation=270.0,/normal
      
    ; mark dm overdensity range
    cgPlot,[0.0,0.0],yrange,line=0,color=fsc_color('light gray'),/overplot
    cgPlot,[1.0,1.0],yrange,line=0,color=fsc_color('light gray'),/overplot

    cgPlot,dmbinCen,histVel_ga.mean,line=1,color=getColor(1),/overplot
    cgPlot,dmbinCen,histVel_ga.median,line=2,color=getColor(1),/overplot
    cgPlot,dmbinCen,histVel_ga.mean+sigFac*histVel_ga.stddev,line=0,color=getColor(1),/overplot
    cgPlot,dmbinCen,histVel_ga.mean-sigFac*histVel_ga.stddev,line=0,color=getColor(1),/overplot
    
    cgPlot,dmbinCen,histVel_ar.mean,line=1,color=getColor(3),/overplot
    cgPlot,dmbinCen,histVel_ar.median,line=2,color=getColor(3),/overplot
    cgPlot,dmbinCen,histVel_ar.mean+sigFac*histVel_ar.stddev,line=0,color=getColor(3),/overplot
    cgPlot,dmbinCen,histVel_ar.mean-sigFac*histVel_ar.stddev,line=0,color=getColor(3),/overplot
  
    ; x-axis label
    cgText,x1,y0-0.1,textoidl("log ( \rho_{DM} / <\rho_{DM}> )"),alignment=0.5,/normal
  
  end_PS
  
  ; plot (2) - plot ratio of arepo/gadget for each quantity over the filament range [0,1]
  start_PS, sPa.plotPath + 'gasprops.'+subsetTag+'ratio.'+str(sPa.res)+'_'+str(sPa.snap)+'.eps'
            
    xrange = [0,1]
    yrange = [0.25,2.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Arepo/Gadget Ratio (Mean/Median)",xtitle=textoidl("log ( \rho_{DM} / <\rho_{DM}> )")
      
    cgPlot,xrange,[1.0,1.0],line=0,color=cgColor('light gray'),/overplot
    cgPlot,dmBinCen,10.0^histTemp_ar.mean/10.0^histTemp_ga.mean,line=1,color=cgColor('brown'),/overplot
    cgPlot,dmBinCen,10.0^histDens_ar.mean/10.0^histDens_ga.mean,line=1,color=cgColor('forest green'),/overplot
    cgPlot,dmBinCen,10.0^histEnt_ar.mean/10.0^histEnt_ga.mean,line=1,color=cgColor('slate blue'),/overplot
    cgPlot,dmBinCen,histVel_ar.mean/histVel_ga.mean,line=1,color=cgColor('orange'),/overplot
    
    cgPlot,dmBinCen,10.0^histTemp_ar.median/10.0^histTemp_ga.median,line=2,color=cgColor('brown'),/overplot
    cgPlot,dmBinCen,10.0^histDens_ar.median/10.0^histDens_ga.median,line=2,color=cgColor('forest green'),/overplot
    cgPlot,dmBinCen,10.0^histEnt_ar.median/10.0^histEnt_ga.median,line=2,color=cgColor('slate blue'),/overplot
    cgPlot,dmBinCen,histVel_ar.median/histVel_ga.median,line=2,color=cgColor('orange'),/overplot
    
    ; legend
    legend,['temperature','density','entropy','velocity'],$
      textcolors=['brown','forest green','slate blue','orange'],box=0,/top,/left
      
  end_PS
stop

end

; plotGasSinglePropODRangeComp(): compare histograms of some gas property (dens,temp,...) selected within
;                                 a range of associated DM overdensity

pro plotGasSinglePropODRangeComp, redshift=redshift

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if n_elements(redshift) eq 0 then message,'Specify redshift.'

  ; config
  dmOverDensRange = [1.0,2.0] ; log (rho/mean rho)
  fieldName = 'temperature'
  res = 512
  
  sPa = simParams(res=res,run='arepo',redshift=redshift)
  sPg = simParams(res=res,run='gadget',redshift=redshift)
  
  ; make gas selection
  gasDMDens_ar = findDMOverDensAtGas(sP=sPa)
  wGas_ar = where(gasDMDens_ar ge dmOverDensRange[0] and gasDMDens_ar lt dmOverDensRange[1],count_ar)
  if count_ar eq 0 then message,'Error: No arepo gas found'
  gasDMDens_ar = !NULL
  
  gasDMDens_ga = findDMOverDensAtGas(sP=sPg)
  wGas_ga = where(gasDMDens_ga ge dmOverDensRange[0] and gasDMDens_ga lt dmOverDensRange[1],count_ga)
  if count_ga eq 0 then message,'Error: No gadget gas found'
  gasDMDens_ga = !NULL
  
  print,count_ar,count_ga
  
  ; load gas property
  if fieldName eq 'density' then begin
    gasProp_ar = loadSnapshotSubset(sP=sPa, field=fieldName, partType='gas')
    gasProp_ga = loadSnapshotSubset(sP=sPg, field=fieldName, partType='gas')
    
    ; normalize to critical baryon density at this redshift
    gasProp_ar = alog10(rhoRatioToCrit(gasProp_ar[wGas_ar],redshift=sPa.redshift))
    gasProp_ga = alog10(rhoRatioToCrit(gasProp_ga[wGas_ga],redshift=sPa.redshift))
    
    fieldBinSize = 0.1 ; for histogram
    fieldAxis    = "log ( \rho / \rho_{crit,b,z} )"
  endif
  
  if fieldName eq 'entropy' then begin
    ; load u,dens
    u     = loadSnapshotSubset(sP=sPa,partType='gas',field='u')
    dens  = loadSnapshotSubset(sP=sPa,partType='gas',field='dens')
    gasProp_ar = calcEntropyCGS(u[wGas_ar],dens[wGas_ar],/log,sP=sPa)
    u     = !NULL
    dens  = !NULL
    
    u     = loadSnapshotSubset(sP=sPg,partType='gas',field='u')
    dens  = loadSnapshotSubset(sP=sPg,partType='gas',field='dens')
    gasProp_ga = calcEntropyCGS(u[wGas_ga],dens[wGas_ga],/log,sP=sPg)
    u     = !NULL
    dens  = !NULL
    
    fieldBinSize = 0.1 ; for histogram
    fieldAxis    = "log ( entropy ) [K cm^2]"
  endif
  
  if fieldName eq 'temperature' then begin
    ; load u,nelec
    u     = loadSnapshotSubset(sP=sPa,partType='gas',field='u')
    nelec = loadSnapshotSubset(sP=sPa,partType='gas',field='nelec')
    gasProp_ar = convertUtoTemp(u[wGas_ar],nelec[wGas_ar],/log)
    u     = !NULL
    nelec = !NULL
    
    u     = loadSnapshotSubset(sP=sPg,partType='gas',field='u')
    nelec = loadSnapshotSubset(sP=sPg,partType='gas',field='nelec')
    gasProp_ga = convertUtoTemp(u[wGas_ga],nelec[wGas_ga],/log)
    u     = !NULL
    nelec = !NULL
    
    fieldBinSize = 0.1 ; for histogram
    fieldAxis    = "log ( T_{gas} ) [K]"
  endif
  
  ; histogram gas property distributions and normalize by total
  hist_ar = histogram(gasProp_ar,binsize=fieldBinSize,loc=loc_ar)
  hist_ga = histogram(gasProp_ga,binsize=fieldBinSize,loc=loc_ga)
  
  hist_ar /= float(total(hist_ar))
  hist_ga /= float(total(hist_ga))
  
  ; plot
  start_PS, sPa.plotPath + 'gasprop.'+fieldName+'.'+string(dmOverDensRange[0],format='(f4.1)')+'-'+$
            string(dmOverDensRange[1],format='(f4.1)')+'.'+str(sPa.res)+'_'+str(sPa.snap)+'.eps'
            
    xrange = [floor(min([loc_ar,loc_ga])),ceil(max([loc_ar,loc_ga]))]
    yrange = [1e-5,5e-1]
    
    title = string(dmOverDensRange[0],format='(f4.1)')+" < "+textoidl("log ( \rho_{DM} / <\rho_{DM}> )") + $
            " < "+string(dmOverDensRange[1],format='(f4.1)')+" ("+str(sPa.res)+textoidl("^3")+" z="+string(sPa.redshift,format='(f3.1)')+")"
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="Fraction",xtitle=textoidl(fieldAxis),title=title
      
    ; mark dm overdensity range
    ;cgPlot,[dmOverDensRange[0],dmOverDensRange[0]],yrange,line=0,color=fsc_color('light gray'),/overplot
    ;cgPlot,[dmOverDensRange[1],dmOverDensRange[1]],yrange,line=0,color=fsc_color('light gray'),/overplot

    ; histograms    
    cgPlot,loc_ar+fieldBinSize*0.5,hist_ar,line=0,color=getColor(3),/overplot
    cgPlot,loc_ga+fieldBinSize*0.5,hist_ga,line=0,color=getColor(1),/overplot
  
    ; legend
    legend,['gadget','arepo'],textcolors=getColor([1,3],/name),box=0,/right,/top
  
  end_PS
end