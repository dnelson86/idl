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

; plotGasSinglePropODRangeComp(): compare histograms of some gas property (dens,temp,...) selected within
;                                 a range of associated DM overdensity

pro plotGasSinglePropODRangeComp

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  dmOverDensRange = [0.0,1.0] ; log (rho/mean rho)
  fieldName = 'density'
  res = 512
  
  sPa = simParams(res=res,run='arepo',redshift=0.0)
  sPg = simParams(res=res,run='gadget',redshift=0.0)
  
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
    fieldAxis    = 0
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
    fieldAxis    = textoidl("log ( T_{gas} ) [K]")
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
      ytitle="Fraction",xtitle=textoidl("log ( \rho / \rho_{crit,b,z} )"),title=title
      
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