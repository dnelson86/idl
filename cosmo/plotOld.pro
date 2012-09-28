; plotOld.pro
; gas accretion project - unused/exploratory plots
; dnelson aug.2012

; checkBlobsHotCold

pro checkBlobsHotCold

  compile_opt idl2, hidden, strictarr, strictarrsubs

  sP = simParams(res=256,run='gadget',redshift=2.0)

  ; config
  sizeFac = 2.0 ;3.5  ; times rvir
  cutFac  = 1.0       ; times boxSize
  exclude_flag = 1
  
  gcMassRange = [11.5,11.75]
  
  ; target list
  gc    = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sP)
  
  gcIDsPri = gcIDList(gc=gc,select='pri')
  gcMassesPri = codeMassToLogMsun(gc.subgroupMass[gcIDsPri])
  
  w = where(gcMassesPri gt gcMassRange[0] and gcMassesPri le gcMassRange[1])
  gcIDs = gcIDsPri[w]
  print,gcIDs

  ; load gas positions and densities
  pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  dens = loadSnapshotSubset(sP=sP,partType='gas',field='dens')

  gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')

  ; exclude all secondary subgroups?
  if exclude_flag eq 1 then begin
    print,'excluding secondaries...'
    gc2 = loadGroupCat(sP=sP,/readIDs)
    secPIDList = gcPIDList(gc=gc2,select='sec',partType='all')
    
    ; TEMP load galaxycat and mtS, goal: get gas ID list in cold and hot modes
    galcat = galaxyCat(sP=sP)
    mt = mergerTreeSubset(sP=sP)
    mtSIDs = { gal : galcat.galaxyIDs[mt.galcatSub.gal], gmem : galcat.groupmemIDs[mt.galcatSub.gmem] }
    
    ;accTvir = gcSubsetProp(sP=sP,select='all',/accTvir,/mergerTreeSubset,/accretionTimeSubset)
    maxTemp = gcSubsetProp(sP=sP,select='all',/maxPastTemp)

    hotIDs = [galcat.galaxyIDs[where(maxTemp.gal gt 5.9,count1)], $
              galcat.groupmemIDs[where(maxTemp.gmem gt 5.9,count2)]]
    coldIDs = [galcat.galaxyIDs[where(maxTemp.gal le 5.9,count3)], $
              galcat.groupmemIDs[where(maxTemp.gmem le 5.9,count4)]]
    
    print,count1+count2,count3+count4,n_elements(gas_ids)
    
    ; make sure no collision
    match,hotIDs,coldIDs,ind1,ind2,count=countMatch
    if countMatch gt 0 then message,'Error: Collision.'
    
    ; make complement mask
    match,gas_ids,secPIDList,ind1,ind2,count=countMatch
    mask = bytarr(n_elements(gas_ids))
    mask[ind1] = 1B
    wInd = where(mask eq 0B,count)
    
    if count eq 0 then message,'Error'
    print,'Removing ['+str(countMatch)+'] of ['+str(n_elements(gas_ids))+'] have left: '+$
      string(float(count)*100/n_elements(gas_ids),format='(f4.1)')+'%'
    
    ; restrict
    dens = dens[wInd]
    pos = pos[*,wInd]
    gas_ids = gas_ids[wInd]
  endif 
  
  ; counters
  tot_hot = 0L
  tot_cold = 0L
  tot_count = 0L
  
  print,'calculating...'
  ; loop over all requested halos and image
  foreach gcID, gcIDs do begin
  
    ; get subhalo position and size of imaging box
    boxCen     = sgcen[*,gcID]
    boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcID]] / 10.0) * 10.0
    boxSizeImg = [boxSize,boxSize,boxSize] ; cube
  
    ; make conservative cutout greater than boxsize accounting for periodic (do cube not sphere)
    xDist = pos[0,*] - boxCen[0]
    yDist = pos[1,*] - boxCen[1]
    zDist = pos[2,*] - boxCen[2]
    
    correctPeriodicDistVecs, xDist, sP=sP
    correctPeriodicDistVecs, yDist, sP=sP
    correctPeriodicDistVecs, zDist, sP=sP
    
    rvir = gc.group_r_crit200[gc.subgroupGrNr[gcID]]
    rad = reform(sqrt(xDist*xDist + yDist*yDist + zDist*zDist)) / rvir[0]
  
    ; local (cube) cutout
    ;wCut = where(abs(xDist) le 0.5*cutFac*boxSize and abs(yDist) le 0.5*cutFac*boxSize and $
    ;             abs(zDist) le 0.5*cutFac*boxSize,nCutout)
                 
    ; local (cube) cutout - only blobs
    wCut = where(abs(xDist) le 0.5*cutFac*boxSize and abs(yDist) le 0.5*cutFac*boxSize and $
                 abs(zDist) le 0.5*cutFac*boxSize and alog10(dens*1e10) gt 5.0 and rad gt 0.125,nCutout)

    ; local (cube) cutout - only non-blobs
    ;wCut = where(abs(xDist) le 0.5*cutFac*boxSize and abs(yDist) le 0.5*cutFac*boxSize and $
    ;             abs(zDist) le 0.5*cutFac*boxSize and (alog10(dens*1e10) le 5.0 or rad le 0.125),nCutout)
    
    ; check hot or cold? get IDs inside selection and match to hot/cold IDs
    loc_ids = gas_ids[wCut]
    
    match,loc_ids,hotIDs,ind1,ind2,count=count_hot
    match,loc_ids,coldIDs,ind1,ind2,count=count_cold
    
    print,'['+string(gcID,format='(i5)')+'] hot: '+str(count_hot)+'  cold: '+str(count_cold)+'  tot: '+str(n_elements(loc_ids))
    
    tot_cold  += count_cold
    tot_hot   += count_hot
    tot_count += n_elements(loc_ids) 
    
    ; TEMP
    ;loc_dens = alog10(dens[wCut]*1e10)
    ;start_PS,'test_rhor_'+str(gcID)+'.eps'
    ;  cgplot,rad,loc_dens,psym=3,$
    ;    xtitle="radius [kpc]",ytitle="log (dens) [msun/kpc3]"
    ;  w = where(loc_dens gt 5.0 and rad gt 0.1)
    ;  cgplot,rad[w],loc_dens[w],psym=3,color=cgColor('red'),/overplot
    ;end_PS
    ; END TEMP
  
  endforeach ;gcIDs
  
  hot_frac  = float(tot_hot) / tot_count
  cold_frac = float(tot_cold) / tot_count
  print,'hot frac: '+str(hot_frac)+'  cold frac: '+str(cold_frac)

end

; plotTempVsRad

pro plotTempVsRad

  sP = simParams(res=128,run='gadget',redshift=2.0)

  ; get gas/tracer positions with time
  at  = accretionTraj(sP=sP)
  mt  = mergerTreeSubset(sP=sP)

  ; some radius/temp selection
  rad = sqrt(at.relPos_gal[*,0,*] * at.relPos_gal[*,0,*] + $
             at.relPos_gal[*,1,*] * at.relPos_gal[*,1,*] + $
             at.relPos_gal[*,2,*] * at.relPos_gal[*,2,*])
  rad = reform(rad)
  temp = at.curTemp_gal
  
  ; normalize by rvir(t)
  for i=0,n_elements(mt.times)-1 do begin
    rad[i,*] /= mt.hVirRad[i,mt.gcIndOrig.gal]
  endfor

  ; plot (0)
  start_PS, sP.plotPath + 'tempvsrad.'+str(sP.res)+'_'+str(sP.snap)+'.eps'
    !p.thick -= 2
    xrange = [0,2]
    yrange = [3,7]
    
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="temp",xtitle="rad / rvir"      
      
      w = where(temp ne 0.0,count,comp=wc,ncomp=ncomp)
      if count gt 0 then cgPlot,rad[w],temp[w],psym=3,/overplot
    
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

; -------------------------------------------------------------------------------------------------
; vsHaloMass:
; -------------------------------------------------------------------------------------------------

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
