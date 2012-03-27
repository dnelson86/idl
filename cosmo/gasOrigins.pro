; gasOrigins.pro
; gas accretion project - OLD "gas origins" (UNUSED)
; dnelson mar.2012

; gasOrigins(): from sP.snap load the galaxy catalog and consider the evolution of all the gas 
;               elements snapshot by snapshot backwards in time until either endingSnap (if specified) or
;               back five snapshots (default). at each step save:
;                 1. the radial distance of each from its primary/secondary parent
;                 2. temperature and entropy
;                 3. implicitly, whether the gas is in a primary or secondary subhalo (corresponding to 
;                    the two radial distances being equal or different)              

function gasOrigins, sP=sP, endingSnap=endingSnap

  compile_opt idl2, hidden, strictarr, strictarrsubs
  print,'think more about this function (gadget only anyways?)' & stop
  ; TODO: make this use the tracers if running on an Arepo SP?

  ; if snap specified, run only one snapshot (and/or just return previous results)
  if n_elements(endingSnap) gt 0 then begin
    saveFilename = sP.derivPath + 'gas.origins.'+str(sP.res)+'.'+str(endingSnap)+'-'+str(sP.snap)+'.sav'
    
    ; results exist, return
    if (file_test(saveFilename)) then begin
      restore,saveFilename
      r = {temp_gal:temp_gal,temp_gmem:temp_gmem,entropy_gal:entropy_gal,entropy_gmem:entropy_gmem,$
           indMatch:indMatch}
      return,r
    endif
    
    ; need to compute, set restricted range of snapshots to process
    snapRange = [sP.snap,endingSnap]
  endif else begin
    ; default config
    numSnapsBack = 10
    
    snapRange = [sP.snap,sP.snap-numSnapsBack]
  endelse  
  
  ; load galaxy catalog at starting redshift
  galcat = galaxyCat(sP=sP)
  
  print,'Loaded  ['+str(n_elements(gc.galaxyIDs))+'] ['+str(n_elements(gc.groupmemIDs))+'] from galCat.'
  
  for m=snapRange[0],snapRange[1],-1 do begin
  
    ; set save filename and check for existence
    sP.snap = m
    saveFilename = sP.derivPath + 'gas.origins.'+str(sP.res)+'.'+str(endingSnap)+'-'+str(sP.snap)+'.sav'
    
    if (file_test(saveFilename)) then begin
      print,'Skipping: '+strmid(saveFilename,strlen(sP.derivPath))
      continue
    endif  
  
    ; load gas IDs and match
    ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    
    ; IMPORTANT! rearrange ids_ind to be in the order of galcat.IDs, need this if we want ids[ids_ind], 
    ; temp[ids_ind], etc to be in the same order as the galaxy catalog id list     
    match,galcat.galaxyIDs,ids,galcat_ind,ids_gal_ind,count=countGal,/sort
    ids_gal_ind = ids_gal_ind[sort(galcat_ind)]
    
    match,galcat.groupmemIDs,ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
    ids_gmem_ind = ids_gmem_ind[sort(galcat_ind)]
    
    ids    = !NULL
    galcat_ind = !NULL
    
    ; load u,nelec and calculate temp of gas
    u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
    nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')

    temp_gal  = convertUtoTemp(u[ids_gal_ind], nelec[ids_gal_ind])
    temp_gmem = convertUtoTemp(u[ids_gmem_ind],nelec[ids_gmem_ind])
    
    nelec = !NULL
    
    ; load gas density to calculate entropy
    dens = loadSnapshotSubset(sP=sP,partType='gas',field='density')
    
    entropy_gal  = calcEntropy(u[ids_gal_ind], dens[ids_gal_ind], sP=sP)
    entropy_gmem = calcEntropy(u[ids_gmem_ind],dens[ids_gmem_ind], sP=sP)
    
    u    = !NULL
    dens = !NULL
    
    ; load the galaxy catalog at this redshift and match IDs from the galCat at the target redshift
    galcatCur = galaxyCat(sP=sP)
    
    match,galcat.galaxyIDs,galcatCur.galaxyIDs,galcat_ind_gal,galcatCur_ind_gal,count=countGal,/sort
    match,galcat.galaxyIDs,galcatCur.groupmemIDs,galcat_ind_gal2,galcatCur_ind_gal2,count=countGal2,/sort
    match,galcat.groupmemIDs,galcatCur.groupmemIDs,galcat_ind_gmem,galcatCur_ind_gmem,count=countGmem,/sort
    match,galcat.groupmemIDs,galcatCur.galaxyIDs,galcat_ind_gmem2,galcatCur_ind_gmem2,count=countGmem2,/sort
    
    galcatCur = !NULL
    
    print,'['+str(endingSnap-m)+'] Matched ['+str(countGal)+' + '+str(countGal2)+'] ['+$
          str(countGmem)+' + '+str(countGmem2)+'].'

    ; keep the match indices relating gas elements in the two galaxy catalogs
    ; NOTE: -1 for empty here    
    indMatch = {galcat_ind_gal:galcat_ind_gal,         galcat_ind_gal2:galcat_ind_gal2,         $
                galcatCur_ind_gal:galcatCur_ind_gal,   galcatCur_ind_gal2:galcatCur_ind_gal2,   $
                galcat_ind_gmem:galcat_ind_gmem,       galcat_ind_gmem2:galcat_ind_gmem2,       $
                galcatCur_ind_gmem:galcatCur_ind_gmem, galcatCur_ind_gmem2:galcatCur_ind_gmem2}
    
    ; save
    save,temp_gal,temp_gmem,entropy_gal,entropy_gmem,indMatch,filename=saveFilename
    print,'    Saved: '+strmid(saveFilename,strlen(sP.derivPath))

  endfor

end

; plotGasOrigins(): plot temperature, entropy, and primary/secondary radius for all gas elements
;                   from the galaxy catalog at a target redshift as a function of time backwards

pro plotGasOrigins, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  numSnapsBack = 5           ; how many steps backwards to plot from target redshift
  radBounds    = [0.3,0.325] ; fraction of r_vir in target redshift
  gcType       = 'gmem'      ; plot which component, gal or gmem
  
  colors = ['black','forest green','slate blue','crimson','orange','saddle brown']

  ; subhalo selection config
  sgSelect = 'pri' ; pri,sec,all
  parNorm  = 'sec' ; pri,sec (normalize r_vir by primary or secondary parent)
  
  ; make legend
  legStrs = []
  timeStart = snapNumToAge(sP=sP)
  
  for i=0,numSnapsBack do begin
    timeDiff = (timeStart - snapnumToAge(snap=sP.snap-i)) * 1000 ;Myr
    curStr   = textoidl('\Delta')+"t = "+string(timeDiff,format='(i3)')+" Myr"
    legStrs = [legStrs,curStr]
  endfor
  
  ; plot (1) - overplotted histograms of temp vs number of snapshots back
  plotBase = "gas.origins.temp_"+sgSelect+'_'+parNorm+'_'+str(minNumPart)+'_'+gcType
  plotName = sP.plotPath + plotBase + '_'+str(sP.res)+'_'+str(sP.snap)+'.eps'
  
  start_PS, plotName
    
    ; loop backwards over snapshots
    for m=sP.snap,sP.snap-numSnapsBack,-1 do begin
   
      ; load
      gcTemp = gcSubsetProp(sP=sP,select=sgSelect,oSnap=m,/curTemp)

      ; load galaxy radii catalog to make radial cut at target redshift
      if (m eq sP.snap) then $
        gcRad = gcSubsetProp(sP=sP,select=sgSelect,parNorm=parNorm,oSnap=m,/rVirNorm)

      ; select between galaxy and group member
      if (gcType eq 'gal') then begin
        ; radial selection
        if (m eq sP.snap) then begin
          wGal  = where(gcRad.gal gt radBounds[0] and gcRad.gal le radBounds[1],countGal)
          print,' Gal Radial cut: '+str(countGal)+' of '+str(n_elements(gcRad.gal))
        endif
        
        data = gcTemp.gal[wGal]
      endif
      
      if (gcType eq 'gmem') then begin
        ; radial selection
        if (m eq sP.snap) then begin
          wGmem = where(gcRad.gmem gt radBounds[0] and gcRad.gmem le radBounds[1],countGmem)
          print,' Gmem Radial cut: '+str(countGmem)+' of '+str(n_elements(gcRad.gmem))
        endif
        
        data = gcTemp.gmem[wGmem]
      endif

      ; histogram
      xrange  = [min(data)*0.9,max(data)*1.1]
      binSize = 0.05
      
      hist = histogram(data,binsize=binSize,locations=xpts,min=xrange[0],max=xrange[1])
      
      ; first plot make axes and bounds
      yrange = [10.0,max(hist)*2.0]
      
      if (m eq sP.snap) then $
        fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,$
           title=str(sP.res)+textoidl("^3")+" z = "+string(redshift,format='(f3.1)')+" "+gcType,$
           xtitle="log ( T [K] )",ytitle="Count",/xs,/ys,/ylog
           
      ; plot
      w = where(hist ne 0)
      
      fsc_plot,xpts[w]+binSize/2.0,hist[w],/overplot,line=0,thick=!p.thick,color=fsc_color(colors[sP.snap-m])
               
    endfor
       
    ; legend    
    legend,legStrs,textcolors=colors,box=0,margin=0.25,/right       
               
  end_PS
  
  ; plot (2) - radial distances
  plotBase = "gas.origins.rad_"+sgSelect+'_'+parNorm+'_'+str(minNumPart)+'_'+gcType
  plotName = sP.plotPath + plotBase + '_'+str(sP.res)+'_'+str(sP.snap)+'.eps'
  
  start_PS, plotName
    
    ; loop backwards over snapshots
    for m=sP.snap,sP.snap-numSnapsBack,-1 do begin
      ; load
      gcRad = gcSubsetProp(sP=sP,select=sgSelect,parNorm=parNorm,oSnap=m,/rVirNorm)
      
      ;w = where(gcRad.gal eq -1,count)
      ;print,'['+str(m)+'] Not matched: '+str(count)+' of '+str(n_elements(gcRad.gal))
      
      ; select between galaxy and group member
      if (gcType eq 'gal') then begin
        ; radial selection
        if (m eq sP.snap) then begin
          wGal  = where(gcRad.gal gt radBounds[0] and gcRad.gal le radBounds[1],countGal)
        endif
        
        data = gcRad.gal[wGal]
      endif
      
      if (gcType eq 'gmem') then begin
        ; radial selection
        if (m eq sP.snap) then begin
          wGmem = where(gcRad.gmem gt radBounds[0] and gcRad.gmem le radBounds[1],countGmem)
        endif
        
        data = gcRad.gmem[wGmem]
      endif
      
      ; histogram
      xrange  = [0.0,1.1]
      binSize = 0.025
      
      hist = histogram(data,binsize=binSize,locations=xpts,min=xrange[0],max=xrange[1])

      ; first plot make axes and bounds
      yrange = [10.0,max(hist)*2.0]
      
      if (m eq sP.snap) then $
        fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,$
           title=str(sP.res)+textoidl("^3")+" z = "+string(redshift,format='(f3.1)')+" "+gcType,$
           xtitle="r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}"),ytitle="Count",/xs,/ys,/ylog
           
      ; add drop to zero for first histogram, otherwise restrict to nonzero
      w = where(hist gt 0)
      
      if (m eq sP.snap) then begin
        ; for first overplot a faint histogram of the whole distribution
        hist2 = histogram(gcRad.gal,binsize=binSize,locations=xpts2,min=xrange[0],max=xrange[1])
        
        ; strong vertical lines over radial selection
        xpts = [radBounds[0],radBounds[0],xpts[w] + binSize/2.0,radBounds[1],radBounds[1]]
        hist = [yrange[0],(hist2[w[0]]+hist2[w[0]-1])/2.0,hist2[w],$
                (hist2[w[n_elements(w)-1]]+hist2[w[n_elements(w)-1]+1])/2.0,yrange[0]]

        fsc_plot,xpts2+binSize/2.0,hist2,/overplot,line=0,thick=!p.thick,color=fsc_color('light gray')
        
        ; plot
        fsc_plot,xpts,hist,/overplot,line=0,thick=!p.thick+1.0,color=fsc_color(colors[sP.snap-m])
      endif else begin
        xpts = xpts[w] + binSize/2.0
        hist = hist[w]
        
        ; plot
        fsc_plot,xpts,hist,/overplot,line=0,thick=!p.thick+1.0,color=fsc_color(colors[sP.snap-m])
      endelse
      
    endfor
    
    ; legend    
    legend,legStrs,textcolors=colors,box=0,margin=0.25,/right
    
  end_PS
  
end

; plotGasOriginsTracks(): plot some individual tracks rad vs. time and temp vs. rad using gasOrigins
  
pro plotGasOriginsTracks, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  numSnapsBack = 5           ; how many steps backwards to plot from target redshift
  radBounds    = [0.3,0.302] ; fraction of r_vir in target redshift
  minRadCut    = 0.45        ; gas element must reach this r/r_vir at least to be included
  gcType       = 'gal'      ; plot which component, gal or gmem
  
  colors = ['black','forest green','slate blue','crimson','orange','saddle brown']

  ; subhalo selection config
  sgSelect = 'pri' ; pri,sec,all
  parNorm  = 'sec' ; pri,sec (normalize r_vir by primary or secondary parent)
  
  ; make legend
  legStrs = []
  timeStart = snapNumToAge(sP=sP)
  
  for i=0,numSnapsBack do begin
    timeDiff = (timeStart - snapnumToAge(snap=sP.snap-i)) * 1000 ;Myr
    curStr   = textoidl('\Delta')+"t = "+string(timeDiff,format='(i3)')+" Myr"
    legStrs = [legStrs,curStr]
  endfor  
    
  ; loop backwards over snapshots and make tracks
  for m=sP.snap,sP.snap-numSnapsBack,-1 do begin
    ; load
    gcTemp = gcSubsetProp(sP=sP,select=sgSelect,oSnap=m,/curTemp)
    gcRad  = gcSubsetProp(sP=sP,select=sgSelect,parNorm=parNorm,oSnap=m,/rVirNorm)
    
    ; select between galaxy and group member
    if (gcType eq 'gal') then begin
      ; radial selection
      if (m eq sP.snap) then begin
        wGal  = where(gcRad.gal gt radBounds[0] and gcRad.gal le radBounds[1],countGal)
        xdata = fltarr(numSnapsBack+1,countGal)
        ydata = fltarr(numSnapsBack+1,countGal)
      endif
      
      xdata[sP.snap-m,*] = gcRad.gal[wGal]
      ydata[sP.snap-m,*] = gcTemp.gal[wGal]
    endif
    
    if (gcType eq 'gmem') then begin
      ; radial selection
      if (m eq sP.snap) then begin
        wGmem = where(gcRad.gmem gt radBounds[0] and gcRad.gmem le radBounds[1],countGmem)
        xdata = fltarr(numSnapsBack+1,countGmem)
        ydata = fltarr(numSnapsBack+1,countGmem)
      endif
      
      xdata[sP.snap-m,*] = gcRad.gmem[wGmem]
      ydata[sP.snap-m,*] = gcTemp.gmem[wGmem]
    endif
  endfor
  
  sz = size(xdata)
  
  ; remove -1 in xdata (rad) - not found in prior group
  wBad = where(xdata eq -1,count)
  w2d = array_indices(xdata,wBad)
  
  for i=0,n_elements(wBad)-1 do begin
    xdata[wBad[i]] = xdata[w2d[0,i]-1,w2d[1,i]]
  endfor
  
  ; minimum radial distance reached cut
  wR = lindgen(sz[2])
  if (minRadCut ne 0) then $
    wR = where(max(xdata,dim=1) ge minRadCut,comp=wRcomp)
    
  ; plot (1) - temperature vs. radial distances (all)
  xrange = [0.0,1.1]
  yrange = [4.0,max(ydata)*1.1]  
  
  plotBase = "gas.origins.rt_"+sgSelect+'_'+parNorm+'_'+str(minNumPart)+'_'+gcType
  plotName = sP.plotPath + plotBase + '_'+str(sP.res)+'_'+str(sP.snap)+'.eps'
  
  start_PS, plotName
    
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,$
       title=str(sP.res)+textoidl("^3")+" z = "+string(redshift,format='(f3.1)')+" "+gcType+" (all)",$
       xtitle="r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}"),ytitle="log ( T [K] )",/xs,/ys
    
    ; plot tracks
    foreach i,lindgen(sz[2]) do $
      fsc_plot,xdata[*,i],ydata[*,i],/overplot,psym=-4,symsize=0.4,thick=0.5,color=fsc_color('black')
      
    ; oplot colored markers
    for i=0,sz[1]-1 do $
      fsc_plot,xdata[i,*],ydata[i,*],/overplot,psym=4,symsize=0.4,color=fsc_color(colors[i])
      
    ; oplot dotted indicators for -1
    for i=0,n_elements(wBad)-1 do begin
      fsc_plot,[xdata[w2d[0,i]-1,w2d[1,i]],xdata[w2d[0,i]-1,w2d[1,i]]+0.1],[ydata[wBad[i]],ydata[wBad[i]]],$
      line=1,thick=0.5,/overplot
    endfor
      
    ; legend    
    legend,legStrs,textcolors=colors,box=0,margin=0.25,/right
    
  end_PS    
    
  ; plot (2) - temperature vs. radial distances meeting minRadCut
  plotBase = "gas.origins.rt1_"+sgSelect+'_'+parNorm+'_'+str(minNumPart)+'_'+gcType
  plotName = sP.plotPath + plotBase + '_'+str(sP.res)+'_'+str(sP.snap)+'.eps'
  
  start_PS, plotName
    
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,$
       title=str(sP.res)+textoidl("^3")+" z = "+string(redshift,format='(f3.1)')+" "+gcType+" (minCut)",$
       xtitle="r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}"),ytitle="log ( T [K] )",/xs,/ys
    
    ; plot tracks
    foreach i,wR do $
      fsc_plot,xdata[*,i],ydata[*,i],/overplot,psym=-4,symsize=0.4,thick=0.5,color=fsc_color('black')
      
    ; oplot colored markers
    for i=0,sz[1]-1 do $
      fsc_plot,xdata[i,wR],ydata[i,wR],/overplot,psym=4,symsize=0.4,color=fsc_color(colors[i])
      
    ; oplot dotted indicators for -1
    for i=0,n_elements(wBad)-1 do begin
      if (total(i eq wR) gt 0) then $
        fsc_plot,[xdata[w2d[0,i]-1,w2d[1,i]],xdata[w2d[0,i]-1,w2d[1,i]]+0.1],[ydata[wBad[i]],ydata[wBad[i]]],$
        line=1,thick=0.5,/overplot
    endfor
      
    ; legend    
    legend,legStrs,textcolors=colors,box=0,margin=0.25,/right
    
  end_PS
  
  ; plot (3) - other tracks not meeting minRadCut
  
  plotBase = "gas.origins.rt2_"+sgSelect+'_'+parNorm+'_'+str(minNumPart)+'_'+gcType
  plotName = sP.plotPath + plotBase + '_'+str(sP.res)+'_'+str(sP.snap)+'.eps'
  
  start_PS, plotName
  
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,$
       title=str(sP.res)+textoidl("^3")+" z = "+string(redshift,format='(f3.1)')+" "+gcType+" (failed minCut)",$
       xtitle="r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}"),ytitle="log ( T [K] )",/xs,/ys
    
    ; plot tracks
    foreach i,wRcomp do $
      fsc_plot,xdata[*,i],ydata[*,i],/overplot,psym=-4,symsize=0.4,thick=0.5,color=fsc_color('black')
      
    ; oplot colored markers
    for i=0,sz[1]-1 do $
      fsc_plot,xdata[i,wRcomp],ydata[i,wRcomp],/overplot,psym=4,symsize=0.4,color=fsc_color(colors[i])
      
    ; oplot dotted indicators for -1
    for i=0,n_elements(wBad)-1 do begin
      if (total(i eq wRcomp) gt 0) then $
        fsc_plot,[xdata[w2d[0,i]-1,w2d[1,i]],xdata[w2d[0,i]-1,w2d[1,i]]+0.1],[ydata[wBad[i]],ydata[wBad[i]]],$
        line=1,thick=0.5,/overplot
    endfor
      
    ; legend    
    legend,legStrs,textcolors=colors,box=0,margin=0.25,/right
    
  end_PS
  
  stop
  
end