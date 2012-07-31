; plotSphere.pro
; gas accretion project - visualization/plotting of quantities onto healpix spheres
; dnelson may.2012

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

; plotHaloShellDensComp(): compare shell density for four different mass halos at one redshift

pro plotHaloShellDensComp

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  redshift = 2
  partType = 'gas'
  radInds  = [6,11,15]  ; pre-saved radFacs
  minmax   = [-0.6,2.0] ; log (rho/mean rho)
  rot_ang  = [0,0]      ; [60,-45] ;[lat,long] center in deg (left,up)
  cutSubS  = 1          ; cut satellite substructures out from halo
  
  sP = simParams(res=512,run='gadget',redshift=float(redshift))

  ; get IDs of mass targets
  hMassTargets = [12.5,12.0,11.5,11.0]
  subgroupIDs  = massTargetToHaloID(hMassTargets,sP=sP)

  ; plot
  foreach radInd,radInds do begin
    print,radInd
    
    if cutSubS then csTag = '.cutSubS' else csTag = ''
    start_PS, sP.plotPath+'shell_rcomp_z'+str(redshift)+'_r'+str(radInd)+'_'+partType+csTag+'.eps', xs=6*1.5, ys=6
      
      pos = ['ul_nb','ur_nb','ll_nb','lr_nb']
      xtpos = [0.06,0.56,0.06,0.56]
      ytpos = [0.55,0.55,0.13,0.13]
      
      bartitle = "log ( "+textoidl("\rho / <\rho>")+" )"
        
      for i=0,3 do begin
        ; interpolate onto the shell
        hsd = haloShellValue(sP=sP,partType=partType,valName='density',subgroupID=subgroupIDs[i],cutSubS=cutSubS)
  
        ; convert densities into ratios to the mean
        healpix_data = alog10(10.0^hsd.val_dens[*,radInd] / mean(10.0^hsd.val_dens[*,radInd]))
        healpix_data = reform(healpix_data)
        
        title = sP.run+" "+str(sP.res)+textoidl("^3")+"  z = "+string(sP.redshift,format='(f3.1)')+" "+$
                textoidl("\rho_{"+hsd.partType+"} (r / r_{vir} = "+string(hsd.radFacs[radInd],format='(f4.2)'))+")"
  
        if i eq 0 then $
          plotMollweideProj,healpix_data,rot_ang=rot_ang,title=title,bartitle=bartitle,$
            minmax=minmax,/bigbar,pos=pos[i]
        if i gt 0 then $
          plotMollweideProj,healpix_data,rot_ang=rot_ang,title="",bartitle="",$
            minmax=minmax,pos=pos[i],/noerase
  
        cgText,xtpos[i],ytpos[i],"M = "+string(hMassTargets[i],format='(f4.1)'),$
          /normal,alignment=0.5,color=cgColor('dark gray'),charsize=1.0
      endfor
    end_PS, pngResize=60, /deletePS
  
  endforeach ;radInds
  
end

; plotHaloShellValueComp(): compare four different particle fields for one halo at one redshift

pro plotHaloShellValueComp

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  redshift = 2
  
  sP = simParams(res=512,run='gadget',redshift=float(redshift))
  
  subgroupIDs = [981] ;z2.304 g2342 a2132
                      ;z2.301 g2289 a2034
                      ;z2.314 g981 a927
  
  radInds  = [7]        ; pre-saved radFacs (3=0.25, 4=0.5, 7=rvir)
  rot_ang  = [0,45]      ; [lat,long] center in deg (left,up)
  cutSubS  = 1          ; cut satellite substructures out from halo
  
  ; value config
  partTypes = ['gas','gas','gas','gas',$
               'gas','gas','gas','gas'] ;'dm','dm'
  valNames  = ['temp','density','pressure','radmassflux',$
               'radvel','entropy','metallicity','angmom'] ;dm 'density','veldisp'
  ctNames   = ['helix','helix','helix','brewer-redblue',$
               'brewer-redpurple','helix','helix','helix']
  
  bartitles = ["T_{gas} [_{ }log K_{ }]",$
               "log ( \rho / <\rho> )",$
               "log ( P / k_B ) [_{ }K cm^{-3 }]",$
               "Radial Mass Flux [_{ }M_{sun} kpc^{-2} Myr^{-1 }]",$
               "v_{rad} [_{ }km s^{-1 }]",$
               "log ( Entropy ) [_{ }K cm^{2 }]",$
               ;"log ( \rho_{DM} / <\rho_{DM}> )",$
               ;"\sigma_{vel,DM} [_{ } km/s_{ }]"]
               "Metallicity / 0.0127",$
               "log ( Angular Momentum ) [kpc km/s]"]
               
  ; rvir
  ranges = list([4.3,6.5],[-1.0,1.5],[2.0,3.0],[-1.0,1.0],$
                [-400,400],[6.0,9.0],[0.0,1.0],[3.8,5.0]) ;dm[-0.5,1.0],[0.0,200.0]
        
  ; 0.5/0.25 rvir
  ;ranges = list([4.3,7.0],[-1.0,1.5],[3.0,4.5],[-1.5,1.5],$
  ;              [-400,400],[6.0,9.0],[0.0,2.0],[3.5,4.5]) ;dm[-0.5,1.0],[100.0,300.0]
               
  ratioToMean = [0,1,0,0,0,0,0,0] ; plot value/mean(value) ratio
  plotLog     = [0,1,1,0,0,1,0,1] ; plot log(value)
  symMinMax   = [0,0,0,1,1,0,0,0] ; symmetric range about zero
  mmRound     = [0.1,0.1,0.1,0.1,10,0.1,0,0.1] ; round range to a more equal number (slightly clip)
  
  pos = ['ul_nt','ur_nt','ll_nt','lr_nt']
  
  ; plot
  foreach radInd,radInds do begin
    print,radInd
    
    if cutSubS then csTag = '.cutSubS' else csTag = ''
    start_PS, sP.plotPath+'shell_valcomp_'+sP.savPrefix+str(sP.res)+'_z'+$
      str(redshift)+'_h'+str(subgroupIDs[0])+'_r'+str(radInd)+csTag+'.eps', xs=6*1.5, ys=6

      for i=0,3 do begin
        ; interpolate onto the shell
        hsv = haloShellValue(sP=sP,partType=partTypes[i],valName=valNames[i],$
                             subgroupID=subgroupIDs[0],cutSubS=cutSubS)

        ; convert values into ratios to the mean
        if ratioToMean[i] then healpix_data = reform(hsv.value[*,radInd] / mean(hsv.value[*,radInd])) $
        else healpix_data = reform(hsv.value[*,radInd])
        if plotLog[i] then healpix_data = alog10(healpix_data)
        
        ; calculate appropriate minmax and clip
        ;if symMinMax[i] then begin
        ;  minVal = -1.0 * max(abs(minmax(healpix_data)))
        ;  maxVal = max(abs(minmax(healpix_data)))
        ;endif else begin
        ;  minVal = min(healpix_data)
        ;  maxVal = max(healpix_data)
        ;endelse
       ; 
       ; if mmRound[i] ne 0 then begin
       ;   minVal = round(minVal/mmRound[i])*mmRound[i]
       ;   maxVal = round(maxVal/mmRound[i])*mmRound[i]
       ; endif
  
        minMaxVal = ranges[i]

        w = where(healpix_data gt minMaxVal[1]*0.99,count)
        if count gt 0 then healpix_data[w] = minMaxVal[1] * 0.99
        w = where(healpix_data lt minMaxVal[0]*0.99,count)
        if count gt 0 then healpix_data[w] = minMaxVal[0] * 0.99
        
        print,minMaxVal
  
        if i eq 0 then $
          plotMollweideProj,healpix_data,rot_ang=rot_ang,title="",bartitle=bartitles[i],pos=pos[i],$
            ctName=ctNames[i],minmax=ranges[i]
        if i gt 0 then $
          plotMollweideProj,healpix_data,rot_ang=rot_ang,title="",bartitle=bartitles[i],pos=pos[i],$
            /noerase,ctName=ctNames[i],minmax=ranges[i]
      endfor
    end_PS, pngResize=60;, /deletePS
  
    if cutSubS then csTag = '.cutSubS' else csTag = ''
    start_PS, sP.plotPath+'shell_valcomp2_'+sP.savPrefix+str(sP.res)+'_z'+$
      str(redshift)+'_h'+str(subgroupIDs[0])+'_r'+str(radInd)+csTag+'.eps', xs=6*1.5, ys=6
      
      for i=4,7 do begin
        ; interpolate onto the shell
        hsv = haloShellValue(sP=sP,partType=partTypes[i],valName=valNames[i],$
                             subgroupID=subgroupIDs[0],cutSubS=cutSubS)

        ; convert values into ratios to the mean
        if ratioToMean[i] then healpix_data = reform(hsv.value[*,radInd] / mean(hsv.value[*,radInd])) $
        else healpix_data = reform(hsv.value[*,radInd])
        if plotLog[i] then healpix_data = alog10(healpix_data)

        if valNames[i] eq 'metallicity' then healpix_data /= 0.0127 ; display rescaling

        ; calculate appropriate minmax and clip
        ;if symMinMax[i] then begin
        ;  minVal = -1.0 * max(abs(minmax(healpix_data)))
        ;  maxVal = max(abs(minmax(healpix_data)))
        ;endif else begin
        ;  minVal = min(healpix_data)
        ;  maxVal = max(healpix_data)
        ;endelse
        ;
        ;if mmRound[i] ne 0 then begin
        ;  minVal = round(minVal/mmRound[i])*mmRound[i]
        ;  maxVal = round(maxVal/mmRound[i])*mmRound[i]
        ;endif

        minMaxVal = ranges[i]

        w = where(healpix_data gt minMaxVal[1]*0.99,count)
        if count gt 0 then healpix_data[w] = minMaxVal[1] * 0.99
        w = where(healpix_data lt minMaxVal[0]*0.99,count)
        if count gt 0 then healpix_data[w] = minMaxVal[0] * 0.99
        
        print,minMaxVal
        
        if i eq 0 then $
          plotMollweideProj,healpix_data,rot_ang=rot_ang,title="",bartitle=bartitles[i],pos=pos[i-4],$
            ctName=ctNames[i],minmax=ranges[i]
        if i gt 0 then $
          plotMollweideProj,healpix_data,rot_ang=rot_ang,title="",bartitle=bartitles[i],pos=pos[i-4],$
            /noerase,ctName=ctNames[i],minmax=ranges[i]
  
      endfor
    end_PS, pngResize=60;, /deletePS
  
  endforeach ;radInds

  stop
end

; calcShellInfallFrac(): calculate angular covering fraction of infalling mass flux

function calcShellInfallFrac, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  radFacs     = [0.25,0.5,1.0] ; fractions of the virial radius to save results
  cutSubS     = 1     ; cut satellite substructures out from halo
  minLogMass  = 10.0  ; minimum halo mass
  partType    = 'gas'
  valName     = 'radmassflux'
  threshVal   = 0.0 ; strict inflow
  
  ; check for existence of save
  saveFilename = sP.derivPath + 'shellFrac.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.mm' + str(fix(minLogMass*10)) + '.cut' + str(cutSubS) + '.rad' + $
    str(n_elements(radFacs)) + '.' + partType + '.' + valName + '.sav'
    
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif
  
  ; load subgroup IDs
  gc        = loadGroupCat(sP=sP,/skipIDs)
  priIDs    = gcIDList(gc=gc,select='pri')
  priMasses = codeMassToLogMsun(gc.subgroupMass[priIDs])
  
  w = where(priMasses ge minLogMass,count)
  subgroupIDs = priIDs[w]
  priMasses = priMasses[w]
  
  ; arrays
  r = { fracInfall  : fltarr(n_elements(radFacs),n_elements(subgroupIDs))  ,$
        subgroupIDs : subgroupIDs         ,$
        priMasses   : priMasses           ,$
        count       : count               ,$
        radFacs     : radFacs             ,$
        nRadFacs    : n_elements(radFacs) ,$
        cutSubS     : cutSubS             ,$
        minLogMass  : minLogMass          ,$
        threshVal   : threshVal            }

  ; interpolate all halos (and save)
  ;hsv = haloShellValue(sP=sP,partType=partType,valName=valName,subgroupIDs=subgroupIDs,$
  ;                     cutSubS=cutSubS,radFacs=radFacs)

  print,'calculating...'
  foreach subgroupID,subgroupIDs,k do begin
    if k mod 100 eq 0 then print,k
    
    ; interpolate onto the shell (load)
    hsv = haloShellValue(sP=sP,partType=partType,valName=valName,subgroupIDs=[subgroupID],$
                         cutSubS=cutSubS,radFacs=radFacs)
    
    ; data
    for radInd=0,n_elements(radFacs)-1 do begin
      healpix_data = reform(hsv.value[*,radInd])
      
      ; calculate sky covering fraction with log(rho/mean rho) > threshold
      w = where(healpix_data ge threshVal, countAbove, ncomp=countBelow)
      r.fracInfall[radInd,k] = float(countAbove) / hsv.nPx
      
      ;print,string(subgroupID,format='(i5)')+"  r="+string(radFacs[radInd],format='(f4.2)')+"  "+$
      ;      string(r.fracInfall[radInd,k]*100,format='(f5.2)')+'%'
    endfor
       
  endforeach
  
  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  return, r
end

; plotShellInfallFracComp(): compare arepo/gadget infall covering fraction vs halo mas

pro plotShellInfallFracComp

  sPg = simParams(res=512,run='gadget',redshift=2.0)
  sPa = simParams(res=512,run='tracer',redshift=2.0)
  
  ffG = calcShellInfallFrac(sP=sPg)
  ffA = calcShellInfallFrac(sP=sPa)
  
  ; red/blue as in Vogelsberger+
  colorsA = [getColor24([255,200,200]),getColor24([255,100,100]),getColor24([255,0,0])] ; 128,256,512 AR
  colorsG = [getColor24([200,200,255]),getColor24([100,100,255]),getColor24([0,0,255])] ; 128,256,512 GA
  cInd = 2
  
  ; green/orange higher contrast
  ;colorsA = [getColor24(['91'x,'e5'x,'9b'x]),getColor24(['12'x,'b2'x,'25'x]),getColor24(['03'x,'60'x,'0f'x])]
  ;colorsG = [getColor24(['f6'x,'a1'x,'9c'x]),getColor24(['e6'x,'21'x,'17'x]),getColor24(['7b'x,'0a'x,'04'x])]
  ;cInd = 1
  
  ; bin fractions into halo mass bins and make median lines
  logMassBins=[9.5,10.0,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11.0,$
               11.1,11.25,11.5,11.75,11.9,13.1]
  logMassNBins = n_elements(logMassBins)-1
  logMassBinCen = 0.5 * (logMassBins + shift(logMassBins,-1))
  logMassBinCen = logMassBinCen[0:-2]
  
  medians = { gadget : fltarr(ffG.nRadFacs,logMassNBins) + !values.f_nan ,$
              arepo  : fltarr(ffA.nRadFacs,logMassNBins) + !values.f_nan  }
              
  ; calculate median accretion rate in bins of halo mass
  for i=0,logMassNbins-1 do begin

    w = where(ffG.priMasses gt logMassBins[i] and ffG.priMasses le logMassBins[i+1],count)
    if count gt 0 then for j=0,ffG.nRadFacs-1 do medians.gadget[j,i] = median(ffG.fracInfall[j,w])
      
    w = where(ffA.priMasses gt logMassBins[i] and ffA.priMasses le logMassBins[i+1],count)
    if count gt 0 then for j=0,ffA.nRadFacs-1 do medians.arepo[j,i] = median(ffA.fracInfall[j,w])
    
  endfor  

  ; debug plot (all points)
  start_PS, sPg.plotPath+'fracvsmass.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+str(sPg.res)+'.'+str(sPg.snap)+'.eps'
    ; plot
    xrange = minmax(ffG.priMasses)+[-0.2,0.2]
    yrange = [0.0,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Angular Covering Fraction of Infall",xtitle="log ( Halo Mass ) "+textoidl(" [M_{sun}]")
      
    cgPlot,xrange,[0.5,0.5],line=0,color=cgColor('light gray'),/overplot
    
    for j=0,ffG.nRadFacs-1 do $
      cgPlot,ffG.priMasses,ffG.fracInfall[j,*],psym=4,color=getColor(j),/overplot 
  end_PS

  ; plot
  start_PS, sPg.plotPath+'infallFrac.'+sPg.plotPrefix+'.'+sPa.plotPrefix+'.'+str(sPg.res)+'.'+str(sPg.snap)+'.eps'
    ; plot
    xrange = [10.0,12.5]
    yrange = [0.0,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Angular Covering Fraction of Infall",xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]")
      
    cgPlot,xrange,[0.5,0.5],line=0,color=cgColor('light gray'),/overplot
    
    for j=0,ffG.nRadFacs-1 do $
      cgPlot,logMassBinCen,medians.gadget[j,*],line=j,color=colorsG[cInd],/overplot
    for j=0,ffA.nRadFacs-1 do $
      cgPlot,logMassBinCen,medians.arepo[j,*],line=j,color=colorsA[cInd],/overplot
      
    ; legends
    strings = textoidl('r/r_{vir} = ')+['0.25','0.5','1.0']
    legend,strings,linestyle=indgen(ffG.nRadFacs),linesize=0.25,box=0,/top,/left
    legend,['gadget','arepo'],textcolor=[colorsG[cInd],colorsA[cInd]],box=0,/top,/right
  end_PS

end

; plotHaloShellDensMovie(): 

pro plotHaloShellDensMovie, redshift=redshift, hMassTarget=hMassTarget

  compile_opt idl2, hidden, strictarr, strictarrsubs
    
  ; config
  rot_ang = [0,0] ; [lat,long] center in deg (left,up)
  cutSubS = 0     ; cut satellite substructures out from halo
  
  sP = simParams(res=512,run='gadget',redshift=float(redshift))  
  
  ; select a primary subgroup based on mass
  subgroupID = massTargetToHaloID(hMassTarget,sP=sP)
    
  ; movie configuration
  nFrames  = 450
  radStart = 0.01
  radEnd   = 1.999
  
  ; interpolate onto shells at a set of fixed radii
  print,'interpolating onto shells...'
  hsd_gas = haloShellDensity(sP=sP,partType='gas',subgroupID=subgroupID,cutSubS=cutSubS)
  hsd_dm  = haloShellDensity(sP=sP,partType='dm',subgroupID=subgroupID,cutSubS=cutSubS)  

  ; stepping in radius
  radStep    = (radEnd - radStart) / (nFrames-1)
  frameRadii = radStep * findgen(nFrames) + radStart

  ; normalize each shell by its mean
  for i=0,hsd_gas.nRadFacs-1 do $
    hsd_gas.val_dens[*,i] = alog10(10.0^hsd_gas.val_dens[*,i] / mean(10.0^hsd_gas.val_dens[*,i]))
  for i=0,hsd_dm.nRadFacs-1 do $
    hsd_dm.val_dens[*,i] = alog10(10.0^hsd_dm.val_dens[*,i] / mean(10.0^hsd_dm.val_dens[*,i]))
      
  ; interpolate
  print,'interpolating onto radii...'
  hpShell_gas = fltarr(hsd_gas.nPx,nFrames)
  hpShell_dm  = fltarr(hsd_dm.nPx,nFrames)
  
  for i=0,hsd_gas.nPx-1 do $
    hpShell_gas[i,*] = hermite(hsd_gas.radFacs,hsd_gas.val_dens[i,*],frameRadii)
  for i=0,hsd_dm.nPx-1 do $
    hpShell_dm[i,*] = hermite(hsd_dm.radFacs,hsd_dm.val_dens[i,*],frameRadii)

  print,'rendering frames...'
  for fn=0,nFrames-1 do begin
  ;fn = 100
    print,fn
    radius = frameRadii[fn]

    ; plot 2 vertically
    start_PS, sP.plotPath + 'frames/shell_z'+str(redshift)+'_m'+str(subgroupID)+'_r'+str(fn)+'.eps', xs=8*0.75, ys=8
      !p.multi = [0,2,1]
      
      healpix_data = reform(hpShell_gas[*,fn])
      title = sP.run+" "+str(sP.res)+textoidl("^3")+"  "+textoidl("\rho_{"+hsd_gas.partType+"} (r / r_{vir} = "+$
        string(radius,format='(f5.3)'))+")"
      bartitle = "log ( "+textoidl("\rho / <\rho>")+" )"
      
      plotMollweideProj,healpix_data,rot_ang=rot_ang,title=title,bartitle=bartitle,pos='top'
      
      healpix_data = reform(hpShell_dm[*,fn])
      title = sP.run+" "+str(sP.res)+textoidl("^3")+"  "+textoidl("\rho_{"+hsd_dm.partType+"} (r / r_{vir} = "+$
        string(radius,format='(f5.3)'))+")"
      bartitle = "log ( "+textoidl("\rho / <\rho>")+" )"
      
      plotMollweideProj,healpix_data,rot_ang=rot_ang,title=title,bartitle=bartitle,pos='bottom'
    
      ; halo mass and redshift label
      cgText,0.9,0.05,"log(M) = "+string(hMassTarget,format='(f4.1)'),$
        /normal,alignment=0.5,color=cgColor('dark gray'),charsize=1.0
      cgText,0.9,0.02,"z = "+string(redshift,format='(f3.1)'),$
        /normal,alignment=0.5,color=cgColor('dark gray'),charsize=1.0
    
    end_PS, pngResize=60, /deletePS
  endfor

end

; haloShellAngPowerSpec(): compute and plot the angular power spectrum
;                          e.g. the variation of the spherical harmonic coefficients

pro haloShellAngPowerSpec

  ; config
  radInds      = [6,11,15]  ; pre-saved radFacs
  hMassTargets = [12.5,12.0,11.5,11.0]
  cutSubS  = 1 ; cut satellite substructures out from halo
  redshift = 3
  
  sP = simParams(res=512,run='gadget',redshift=float(redshift))  
 
  ; locate primary subgroupID for mass target
  subgroupIDs = massTargetToHaloID(hMassTargets,sP=sP)
 
  foreach subgroupID,subgroupIDs,m do begin
    ; interpolate onto shells at a set of fixed radii
    hsd_gas = haloShellDensity(sP=sP,partType='gas',subgroupID=subgroupID,cutSubS=cutSubS,/save)
    hsd_dm  = haloShellDensity(sP=sP,partType='dm',subgroupID=subgroupID,cutSubS=cutSubS,/save)
    
    ; choose maximum spherical harmonic order l_max
    ;l_max = fix(3.0 * hsd_gas.nSide - 1)
    l_max = fix(2.0 * hsd_gas.nSide)
  
    ; plot
    start_PS, sP.plotPath + 'powerspec.shell_z'+str(redshift)+'_h'+str(subgroupID)+'.eps'
      
      l_vals = findgen(l_max)
      unit_lambda = 2*!pi / (l_vals + 0.5) ; wavelength on unit sphere
      ang_size = unit_lambda * 180.0 / !pi
      
      first_l = 1
      xrange = [first_l,max(l_vals)]
      yrange = [1e-5,1e0]
      
      !y.margin[1] += 1.0
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,/xlog,/ylog,$
        ytitle=textoidl("l(l+1)C_l/2\pi"),xtitle="l (Spherical Wavenumber)",title=""
        
      ; draw wavelength axis (linear map)
      cgAxis,xs=1,xaxis=1,xtitle="Angular Size [deg]",xrange=[ang_size[first_l],min(ang_size)]
      
      ; function lambdaAxisFunc,index,value
      ;   return,string(2*!pi/(value+0.5),format='(f3.1)')
      ; end
      ;cgAxis,ys=1,yaxis=1,ytitle="Wavelength [ckpc]",ytickformat='lambdaAxisFunc'
      
      ; for this radius, calculate angular power spectrum of gas and dm
      foreach radInd,radInds,k do begin
        print,radInd
        lambda = unit_lambda * hsd_gas.radFacs[radInd] * hsd_gas.rVir ; wavelength in ckpc
        
        healpix_data = alog10(10.0^hsd_gas.val_dens[*,radInd] / mean(10.0^hsd_gas.val_dens[*,radInd]))
        healpix_data = reform(healpix_data)
        ianafast,healpix_data,cl_gas,nlmax=l_max,/cxx,/nested,tmpdir='/n/home07/dnelson/',$
          /silent;,/won,iter_order=2
      
        healpix_data = alog10(10.0^hsd_dm.val_dens[*,radInd] / mean(10.0^hsd_dm.val_dens[*,radInd]))
        healpix_data = reform(healpix_data)
        ianafast,healpix_data,cl_dm,nlmax=l_max,/cxx,/nested,tmpdir='/n/home07/dnelson/',$
          /silent;,/won,iter_order=2
        
        ; plot power spectra
        cgPlot,l_vals,l_vals*(l_vals+1)*cl_gas/2/!pi,line=k,color=getColor(1),/overplot
        cgPlot,l_vals,l_vals*(l_vals+1)*cl_dm/2/!pi,line=k,color=getColor(2),/overplot
      endforeach
      
      ; legends
      legend,[textoidl("r/r_{vir} = ")+string(hsd_gas.radFacs[radInds],format='(f4.2)')],$
        linestyle=indgen(n_elements(radInds)),linesize=0.25,box=0,/top,/right
      legend,['gas','dm'],textcolor=getColor([1,2],/name),box=0,/bottom,/left
      legend,[sP.run+" "+str(sP.res)+textoidl("^3")+" log(M) = "+string(hMassTargets[m],format='(f4.1)')],$
        box=0,/top,/left,textcolor=['light gray']
    end_PS
  endforeach ;subgroupIDs

end
