; filamentSearchPlot.pro
; plotting: locating filamentary structures/streams in halos (inc. using healpix spheres)
; dnelson apr.2012

; plotFilamentProfile():

pro plotFilamentProfile

  sP = simParams(res=512,run='arepo',redshift=2.0)
  subgroupID = 2132 ;z2.304 g2342 a2132
  
  fp = makeFilamentProfile(sP=sP,subgroupID=subgroupID)

  ; make filament number mask
  filMask = intarr(n_elements(fp.fil_dist))
  offset = 0UL
  for i=0,fp.hfs.nFilaments-1 do begin
    filMask[offset:offset+fp.fil_num[i]-1] = i
    offset += fp.fil_num[i]
  endfor

  ; plot (0) - scatter
  start_PS, sP.plotPath + 'fil.scat.'+sP.run+'.eps'

    xrange = [0.0,2.0]
    yrange = [0,100]
      
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      xtitle="Radius",ytitle="Filament Distance"  
    
    cgPlot,fp.fil_rad,fp.fil_dist,psym=3,/overplot
  end_PS  

  ; plot (1) - cross sectional profiles near rvir
  radRange = [0.5,1.5]
  w = where(fp.fil_rad ge radRange[0] and fp.fil_rad lt radRange[1] and filMask eq 0,count)
  
  xrange = [0,100]
  xvals = fp.fil_dist[w]
  
  start_PS, sP.plotPath + 'fil.csec.dens.'+sP.run+'.eps'

    yrange = [1.0,5]
    
    yvals = alog10(fp.fil_dens[w]*1e10)
      
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      xtitle="Filament Distance [ckpc]",ytitle="Gas Density"  
    
    cgPlot,xvals,yvals,psym=3,/overplot
  end_PS  
  
  start_PS, sP.plotPath + 'fil.csec.temp.'+sP.run+'.eps'

    yrange = [4.0,8.0]
    
    yvals = fp.fil_temp[w]
      
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      xtitle="Filament Distance [ckpc]",ytitle="Gas Temperature"  
    
    cgPlot,xvals,yvals,psym=3,/overplot
  end_PS
  
  stop
  
end

; stackFilamentCrossSec():

function stackFilamentCrossSec

  start_PS, 'fil.stack.dens.eps'
 
    xrange = [0,150]
    ;yrange = [1.0,10.0]
    yrange = [2,6]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      xtitle="Filament Distance [ckpc]",ytitle="Gas Density"  
    
    ; arepo
    restore,'filstack.ar.sav',/verbose
    cgPlot,r.rBinCen,alog10(r.binnedVals.dens*1e10),line=0,/overplot,color=cgColor('red')
    ;cgPlot,r.rBinCen,r.binnedRel.dens,line=0,/overplot,color=cgColor('red')
    
    ; gadget
    restore,'filstack.ga.sav',/verbose
    cgPlot,r.rBinCen,alog10(r.binnedVals.dens*1e10),line=0,/overplot,color=cgColor('blue')
    ;cgPlot,r.rBinCen,r.binnedRel.dens,line=0,/overplot,color=cgColor('blue')
  end_PS
  
  start_PS, 'fil.stack.temp.eps'
 
    xrange = [0,150]
    ;yrange = [0.0,1.5]
    yrange = [4,7]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      xtitle="Filament Distance [ckpc]",ytitle="Gas Temp"  
    
    ; arepo
    restore,'filstack.ar.sav',/verbose
    cgPlot,r.rBinCen,r.binnedVals.temp,line=0,/overplot,color=cgColor('red')
    ;cgPlot,r.rBinCen,r.binnedRel.temp,line=0,/overplot,color=cgColor('red')
    
    ; gadget
    restore,'filstack.ga.sav',/verbose
    cgPlot,r.rBinCen,r.binnedVals.temp,line=0,/overplot,color=cgColor('blue')
    ;cgPlot,r.rBinCen,r.binnedRel.temp,line=0,/overplot,color=cgColor('red')
  end_PS


stop
  sP = simParams(res=512,run='arepo',redshift=2.0)
  radIndOffset = 0
  
  ; find shells that have been made and get subgroupIDs
  fileList = file_search(sP.derivPath+'hShells/hShells.gas.density.*.13.sav')
  
  subgroupIDs = lonarr(n_elements(fileList))
  foreach fileName,fileList,k do begin
    strArr = strsplit(fileList[k],".",count=count,/extract)
    if count ne 12 then message,'error'
    subgroupIDs[k] = long(strmid(strArr[-3],1))
  endforeach
  
  ; arrays for each property
  rad  = []
  dens = []
  temp = []
  
  dens_rel = []
  temp_rel = [] 
  
  ; loop over each halo
  foreach subgroupID,subgroupIDs,k do begin
    print,string(k,format='(i3)'),' ',string(n_elements(subgroupIDs),format='(i3)')
    ; load
    fcs = haloFilamentCrossSec(sP=sP, subgroupID=subgroupID, radIndOffset=radIndOffset)
    if n_elements(fcs) eq 0 then continue ; no filaments
    
    hsv_dens = haloShellValue(sP=sP,partType='gas',valName='density',subgroupID=subgroupID,/cutSubS)
    hsv_temp = haloShellValue(sP=sP,partType='gas',valName='temp',subgroupID=subgroupID,/cutSubS)
    
    ; save
    rad = [rad,fcs.ckpcDists]
    dens = [dens,hsv_dens.value[fcs.filPxInds,fcs.radInd]]
    temp = [temp,hsv_temp.value[fcs.filPxInds,fcs.radInd]]
    
    ; save relative
    densMed = median(hsv_dens.value[*,fcs.radInd])
    tempMed = median(10.0^hsv_temp.value[*,fcs.radInd])
    dens_rel = [dens_rel, hsv_dens.value[fcs.filPxInds,fcs.radInd]/densMed]
    temp_rel = [temp_rel, 10.0^hsv_temp.value[fcs.filPxInds,fcs.radInd]/tempMed]
  endforeach
  
  ; manual radial bins
  rBins = [0.0,5.0,10.0,15.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,120.0,140.0,160.0]
  rNBins = n_elements(rBins)-1
  rBinCen = 0.5 * (rBins + shift(rBins,-1))
  rBinCen = rBinCen[0:-2] 
  
  binnedVals = { dens : fltarr(rNBins)    ,$
                 temp : fltarr(rNBins)     }
  binnedRel = { dens : fltarr(rNBins)    ,$
                 temp : fltarr(rNBins)     }  
  
  for i=0,rNBins-1 do begin
    w = where(rad ge rBins[i] and rad lt rBins[i+1],count)
    if count gt 0 then begin
      binnedVals.dens[i] = median(dens[w])
      binnedVals.temp[i] = median(temp[w])
      binnedRel.dens[i]  = median(dens_rel[w])
      binnedRel.temp[i]  = median(temp_rel[w])
    endif
  endfor

  r = { rad:rad, dens:dens, temp:temp, rNBins:rNBins, rBinCen:rBinCen, $
        binnedVals:binnedVals, binnedRel:binnedRel }
  ;save,r,filename='filstack.ar.sav'
  return, r

end

; binFilamentCrossSec():

function binFilamentCrossSec, sP=sP, fcs=fcs, subgroupID=subgroupID

  hsv_dens = haloShellValue(sP=sP,partType='gas',valName='density',subgroupID=subgroupID,/cutSubS)
  hsv_temp = haloShellValue(sP=sP,partType='gas',valName='temp',subgroupID=subgroupID,/cutSubS)

  ; make filament number mask
  filMask = intarr(n_elements(fcs.ckpcDists))
  offset = 0UL
  for i=0,fcs.nFilaments-1 do begin
    filMask[offset:offset+fcs.pxNums[i]-1] = i
    offset += fcs.pxNums[i]
  endfor

  ; manual radial bins
  rBins = [0.0,5.0,10.0,15.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,120.0,140.0,160.0]
  rNBins = n_elements(rBins)-1
  rBinCen = 0.5 * (rBins + shift(rBins,-1))
  rBinCen = rBinCen[0:-2]
  
  ; arrays
  binnedVals = { dens : fltarr(rNBins)    ,$
                 temp : fltarr(rNBins)     }
  binnedRel = { dens : fltarr(rNBins)    ,$
                 temp : fltarr(rNBins)     }     

  for i=0,rNBins-1 do begin
    w = where(fcs.ckpcDists ge rBins[i] and fcs.ckpcDists lt rBins[i+1],count)
    if count gt 0 then begin
      binnedVals.dens[i] = median(hsv_dens.value[fcs.filPxInds[w],fcs.radInd])
      binnedVals.temp[i] = median(hsv_temp.value[fcs.filPxInds[w],fcs.radInd])
      binnedRel.dens[i]  = median(hsv_dens.value[fcs.filPxInds[w],fcs.radInd] / $
                                  median(hsv_dens.value[*,fcs.radInd]))
      binnedRel.temp[i]  = median(10.0^hsv_temp.value[fcs.filPxInds[w],fcs.radInd] / $
                                  median(10.0^hsv_temp.value[*,fcs.radInd]))
    endif
  endfor
  
  r = { hsv_dens:hsv_dens, hsv_temp:hsv_temp, filMask:filMask, rNBins:rNBins, rBinCen:rBinCen, $
        binnedVals:binnedVals, binnedRel:binnedRel }
  return, r

end

; plotFilamentCrossSec():

pro plotFilamentCrossSec

  sPa = simParams(res=512,run='arepo',redshift=2.0)
  sPg = simParams(res=512,run='gadget',redshift=2.0)
  subgroupIDs = [2132,2342] ;z2.304 g2342 a2132
  
  filIDs = list([0,1],[1,0],[2,3],[3,4]) ;[a,g] pairsfor z2.304
  filInd = 2
  
  radIndOffset = +1 ; 0=rvir, -1=0.9, -2=0.75
  
  ; load
  fcsA = haloFilamentCrossSec(sP=sPa, subgroupID=subgroupIDs[0], radIndOffset=radIndOffset)
  fcsG = haloFilamentCrossSec(sP=sPg, subgroupID=subgroupIDs[1], radIndOffset=radIndOffset)

  bvA = binFilamentCrossSec(sP=sPa, fcs=fcsA, subgroupID=subgroupIDs[0])
  bvG = binFilamentCrossSec(sP=sPg, fcs=fcsG, subgroupID=subgroupIDs[1])

  ; plot (0) - histogram
  start_PS, sPg.plotPath + 'fil2.hist.comp.eps'

    xrange = [0.0,150.0]
    yrange = [1,800]
      
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,xtitle="Radius",ytitle="N"  
    
    binsize = 2.0
    hist = histogram(fcsA.ckpcDists,binsize=binsize,loc=loc)
    cgPlot,loc+binsize*0.5,hist,line=0,/overplot,color=cgColor('red')
    
    hist = histogram(fcsG.ckpcDists,binsize=binsize,loc=loc)
    cgPlot,loc+binsize*0.5,hist,line=0,/overplot,color=cgColor('blue')
    
  end_PS  

  ; plot (1) - cross sectional profiles near rvir
  xrange = [0,150]
  
  start_PS, sPg.plotPath + 'fil2.dens.comp.fil_'+str(filInd)+'.eps'

    yrange = [1.0,4.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      xtitle="Filament Distance [ckpc]",ytitle="Gas Density"  
    
    ; arepo
    w = where(bvA.filMask eq filIDs[filInd,0],count)
    xvals = fcsA.ckpcDists[w]
    yvals = alog10(bvA.hsv_dens.value[fcsA.filPxInds[w],fcsA.radInd]*1e10)
    cgPlot,xvals,yvals,psym=3,/overplot,color=cgColor('red')
    cgPlot,bvA.rBinCen,alog10(bvA.binnedVals.dens*1e10),line=0,/overplot,color=cgColor('red')
    
    ; gadget
    w = where(bvG.filMask eq filIDs[filInd,1],count)
    xvals = fcsG.ckpcDists[w]
    yvals = alog10(bvG.hsv_dens.value[fcsG.filPxInds[w],fcsG.radInd]*1e10)
    cgPlot,xvals,yvals,psym=3,/overplot,color=cgColor('blue')
    cgPlot,bvG.rBinCen,alog10(bvG.binnedVals.dens*1e10),line=0,/overplot,color=cgColor('blue')
  end_PS
  
  start_PS, sPg.plotPath + 'fil2.densrel.comp.fil_'+str(filInd)+'.eps'

    yrange = [0.1,10.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      xtitle="Filament Distance [ckpc]",ytitle="Gas Density / Median"  
    
    ; arepo
    w = where(bvA.filMask eq filIDs[filInd,0],count)
    meanVal = median(bvA.hsv_dens.value[*,fcsA.radInd])
    xvals = fcsA.ckpcDists[w]
    yvals = bvA.hsv_dens.value[fcsA.filPxInds[w],fcsA.radInd] / meanVal
    cgPlot,xvals,yvals,psym=3,/overplot,color=cgColor('red')
    cgPlot,bvA.rBinCen,bvA.binnedRel.dens,line=0,/overplot,color=cgColor('red')
    
    ; gadget
    w = where(bvG.filMask eq filIDs[filInd,1],count)
    meanVal = median(bvG.hsv_dens.value[*,fcsG.radInd])
    xvals = fcsG.ckpcDists[w]
    yvals = bvG.hsv_dens.value[fcsG.filPxInds[w],fcsG.radInd] / meanVal
    cgPlot,xvals,yvals,psym=3,/overplot,color=cgColor('blue')
    cgPlot,bvG.rBinCen,bvG.binnedRel.dens,line=0,/overplot,color=cgColor('blue')
  end_PS  

  start_PS, sPg.plotPath + 'fil2.temp.comp.fil_'+str(filInd)+'.eps'

    yrange = [3,8]
      
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      xtitle="Filament Distance [ckpc]",ytitle="Gas Temperature"  
    
    ; arepo
    w = where(bvA.filMask eq filIDs[filInd,0],count)
    xvals = fcsA.ckpcDists[w]
    yvals = bvA.hsv_temp.value[fcsA.filPxInds[w],fcsA.radInd]
    cgPlot,xvals,yvals,psym=3,/overplot,color=cgColor('red')
    cgPlot,bvA.rBinCen,bvA.binnedVals.temp,line=0,/overplot,color=cgColor('red')
    
    ; gadget
    w = where(bvG.filMask eq filIDs[filInd,1],count)
    xvals = fcsG.ckpcDists[w]
    yvals = bvG.hsv_temp.value[fcsG.filPxInds[w],fcsG.radInd]
    cgPlot,xvals,yvals,psym=3,/overplot,color=cgColor('blue')
    cgPlot,bvG.rBinCen,bvG.binnedVals.temp,line=0,/overplot,color=cgColor('blue')
  end_PS  
  
  start_PS, sPg.plotPath + 'fil2.temprel.comp.fil_'+str(filInd)+'.eps'

    yrange = [0.05,2]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      xtitle="Filament Distance [ckpc]",ytitle="Gas Temp / Median"  
    
    ; arepo
    w = where(bvA.filMask eq filIDs[filInd,0],count)
    meanVal = median(10.0^bvA.hsv_temp.value[*,fcsA.radInd])
    xvals = fcsA.ckpcDists[w]
    yvals = 10.0^bvA.hsv_temp.value[fcsA.filPxInds[w],fcsA.radInd] / meanVal
    cgPlot,xvals,yvals,psym=3,/overplot,color=cgColor('red')
    cgPlot,bvA.rBinCen,bvA.binnedRel.temp,line=0,/overplot,color=cgColor('red')
    
    ; gadget
    w = where(bvG.filMask eq filIDs[filInd,1],count)
    meanVal = median(10.0^bvG.hsv_temp.value[*,fcsG.radInd])
    xvals = fcsG.ckpcDists[w]
    yvals = 10.0^bvG.hsv_temp.value[fcsG.filPxInds[w],fcsG.radInd] / meanVal
    cgPlot,xvals,yvals,psym=3,/overplot,color=cgColor('blue')
    cgPlot,bvG.rBinCen,bvG.binnedRel.temp,line=0,/overplot,color=cgColor('blue')
  end_PS 
  
  stop
  
end