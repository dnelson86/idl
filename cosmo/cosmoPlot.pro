; cosmoPlot.pro
; gas accretion project - plots
; dnelson mar.2012

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
;
; curTemp=1     : current temperature
; maxPastTemp=1 : maximum previous gas temperature
; tVirNorm=1    : normalize by the virial temperature of the parent halo

pro plotTempRad2DHisto, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; config
  ;massBins = [0.0,100.0] ; no massbins
  massBins = [9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5] ; log(M)
  
  curTemp     = 0 ; either: current temperature vs radius
  maxPastTemp = 1 ; or: maximum previous temperature vs radius
  
   trPopMax   = 0 ; if maxPastTemp: choose one of population min,mean,max to return
   trPopMin   = 0 ; (choice of one requred for allTR=0, and ignored for allTR=1)
   trPopMean  = 0 ; (...)
   
  allTR    = 1     ; plot results for all tracers individually, not collapsed onto their parents
  tVirNorm = 1     ; normalize temperature by virial temp of parent halo
  parNorm  = 'pri' ; pri,sec (normalize r_vir and temp (if doing tVirNorm) by primary or secondary parent)
  sgSelect = 'pri' ; pri,sec,all (subhalo selection config, only plot gas in this category of subhalos)
  
  ; load group catalog just for counts of objects in each mass bin
  gc = loadGroupCat(sP=sP,/skipIDs)
  subgroupMasses = codeMassToLogMsun(gc.subgroupMass)
  gc = !NULL
  
  ; get normalized r/r_vir
  gcRad = gcSubsetProp(sP=sP,select=sgSelect,/rVirNorm,allTR=allTR,parNorm=parNorm)
  
  ; get current or maxPast temperature
  gcTemp = gcSubsetProp(sP=sP,select=sgSelect,curTemp=curTemp,maxPastTemp=maxPastTemp,allTR=allTR,$
                        trPopMin=trPopMin,trPopMax=trPopMax,trPopMean=trPopMean)
  
  if keyword_set(tVirNorm) then begin
    ; calculate temperatures of parents and take (non-log) ratio
    gcVirTemp = gcSubsetProp(sP=sP,select=sgSelect,/virTemp,allTR=allTR,parNorm=parNorm)
    
    gcTemp.gal  = 10.0^gcTemp.gal / gcVirTemp.gal
    gcTemp.gmem = 10.0^gcTemp.gmem / gcVirTemp.gmem
  endif
  
  if ~keyword_set(allTR) then $
    gcMass = gcSubsetProp(sP=sP,select=sgSelect,/curSingleVal,singleValField='mass',allTR=allTR)
  
  ; calculate masses of parents (for mass binning halos only)
  parentMass = gcSubsetProp(sP=sP,select=sgSelect,/parMass,allTR=allTR)
  
  for j=0,n_elements(massBins)-2 do begin
  
    ; plot config
    xrange = alog10([0.01,3.0])
    yrange = [4.0,7.0]
  
    binSizeRad  = 0.04 / (sP.res/128) ;0.04
    binSizeTemp = 0.05 / (sP.res/128) ;0.04
    
    ; preserve number of bins in log(rad) histogram
    nBinsRad_linear = ceil((10.0^xrange[1]-10.0^xrange[0])/binSizeRad)+1
    binSizeRad_log  = (xrange[1]-xrange[0])/(nBinsRad_linear-1)
    
    if (keyword_set(tVirNorm)) then begin
      ; preserve number of bins in tmax/tvir histogram
      nBinsTemp = ((yrange[1]-yrange[0])/binSizeTemp)+1
      yrange = [0.0,1.4]
      binSizeTemp = (yrange[1]-yrange[0])/(nBinsTemp-1) ;adjust
    endif
    
    ; select members of this parent mass bins and r>0<inf
    wGal  = where(parentMass.gal gt massBins[j] and parentMass.gal le massBins[j+1] and $
                  gcRad.gal gt 0.0 and finite(gcRad.gal),count1)
    wGmem = where(parentMass.gmem gt massBins[j] and parentMass.gmem le massBins[j+1] and $
                  gcRad.gmem gt 0.0 and finite(gcRad.gmem),count2)
    
    wGCMassBin = where(subgroupMasses gt massBins[j] and subgroupMasses le massBins[j+1],count_gc)
    
    print,j,count1,count2,count_gc
    
    temp_gal  = gcTemp.gal[wGal]
    temp_gmem = gcTemp.gmem[wGmem]
    
    rad_gal   = alog10( gcRad.gal[wGal] )
    rad_gmem  = alog10( gcRad.gmem[wGmem] )

    if count1 eq 0 or count2 eq 0 then continue ; no halos in this mass bin
    
    ; do mass weighting
    if allTR eq 0 then begin
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
      if sP.trMCPerCell eq 0 then message,'Error: allTR set but no tracers in this run.'
    endelse

    ; 2d histo plot config
    if keyword_set(curTemp) then begin
      plotBase = "temp_"+sgSelect+'_'+parNorm
      ytitle   = "log ( T"+textoidl("_{cur}")+" )"
    endif
    
    if keyword_set(maxPastTemp) then begin
      plotBase = "tmax_"+sgSelect+'_'+parNorm
      ytitle   = "log ( T"+textoidl("_{max}")+" )"
    endif
    
    if keyword_set(tVirNorm) then begin
      plotBase = "tmax_tvirNorm_"+sgSelect+'_'+parNorm
      ytitle   = "T"+textoidl("_{max,cur}")+" / T"+textoidl("_{vir}")
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
    if keyword_set(curTemp) then begin
      plotName = sP.plotPath + 'temp_rad_slices_'+str(sP.res)+'_'+str(sP.snap)+'_'+$
                 sgSelect+'_'+parNorm+'.eps'    
      xtitle   = "log ( T"+textoidl("_{cur}")+ " ) [K]"
    endif
    
    if keyword_set(maxPastTemp) then begin
      plotName = sP.plotPath + 'tmax_rad_slices_'+str(sP.res)+'_'+str(sP.snap)+'_'+$
                 sgSelect+'_'+parNorm+'.eps'    
      xtitle = "log ( T"+textoidl("_{max}")+ " ) [K]"
    endif
    
    if keyword_set(tVirNorm) then begin
      plotName = sP.plotPath + 'tmax_TvirNorm_rad_slices_'+str(sP.res)+'_'+str(sP.snap)+'_'+$
                 sgSelect+'_'+parNorm+'.eps'
      xtitle = "T"+textoidl("_{max,cur}")+" / T"+textoidl("_{vir}")
    endif
    
    if n_elements(massBins) gt 2 then $
      plotName = strmid(plotName,0,strlen(plotName)-4) + '.mbin='+str(j)+'.eps'
    
    xrange = yrange
    yrange = [0.0,1.05]
    
    ; plot vertical slices
    plotVerticalSlices,rad_gal,rad_gmem,temp_gal,temp_gmem,plotName,binSizeTemp*2.0,yrange,xtitle
    
  endfor ;j
  
end
