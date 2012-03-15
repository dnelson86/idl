; cosmoPlot.pro
; gas accretion project - plots
; dnelson jan.2012

; plotGalCatRedshift(): 

pro plotGalCatRedshift, sP=sP

  units = getUnits()
  
  ; config
  colors = ['forest green','slate blue','black']
  
  minNumGasInGal = 100 ; 64, 0=no minimum
  
  ; plot axes
  plotRedshifts = snapNumToRedshift(/all)
  
  xrange = [6.0,0.0]
  yrange = [11.5,13.5]
    
  ; plot
  plotName = sP.plotPath+'galmass_vs_redshift_'+str(sP.res)+'.eps'
  
  if (str(res[0]) eq 'two') then res = [256,128]    
  if (str(res[0]) eq 'all') then res = [512,256,128]    
    
  start_PS, plotName

  fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,$
           xtitle="Redshift",ytitle="log ( total gas mass [M"+textoidl("_{sun}")+"] )",$
           xs=9,/ys,ymargin=[4,3]
  universeage_axis,xrange,yrange
    
  if (n_elements(res) eq 1) then $
    fsc_text,0.06,0.94,/normal,str(res)+"^3",alignment=0.5
  
  ; loop over each resolution
  for j=0,n_elements(res)-1 do begin    
  
    sP.res = res[j]
  
    ; arrays
    mass_tot_gal  = fltarr(sP.snapRange[1])
    mass_tot_gmem = fltarr(sP.snapRange[1])
    
    ; load constant mass factor (gadget)
    ;massfac = loadSnapshotSubset(sP.simPath,snapNum=0,partType='gas',field='mass')
    ;massfac = massfac[0]
    ; TODO FOR UNEQUAL MASSES
    stop
    
    ; load galaxy and group membership catalogs
    for m=sP.groupcatRange[0],sP.groupcatRange[1]-1 do begin
      sP.snap = m
      galcat = galaxyCat(sP=sP)
      
      mass_tot_gal[m]  = n_elements(galcat.galaxyIDs)
      mass_tot_gmem[m] = n_elements(galcat.groupmemIDs)
      
      ; enforce minimum number of particles cut
      if (minNumGasInGal gt 0) then begin
        wGal = where(galcat.galaxyLen ge minNumGasInGal,countGal)
        wGmem = where(galcat.groupmemLen ge minNumGasInGal,countGmem)
        
        if (countGal ne 0) then $
          mass_tot_gal[m] -= total(galcat.galaxyLen[wGal])
        if (countGmem ne 0) then $
          mass_tot_gmem[m] -= total(galcat.groupmemLen[wGmem])
      endif
    endfor
    
    ; convert masses
    mass_tot      = codeMassToLogMsun((mass_tot_gal+mass_tot_gmem)*massfac)
    mass_tot_gal  = codeMassToLogMsun(mass_tot_gal*massfac)
    mass_tot_gmem = codeMassToLogMsun(mass_tot_gmem*massfac)
    
     ; overplot successive resolutions
    overplot = 1 ;1,1,1
    line     = j ;0,1,2
    thick    = !p.thick + 1 - 2*(j gt 0) ;3,1,1

    fsc_plot,plotRedshifts,mass_tot_gal,color=fsc_color(colors[0]),$
             overplot=overplot,line=line,thick=thick
    fsc_plot,plotRedshifts,mass_tot_gmem,color=fsc_color(colors[1]),$
             overplot=overplot,line=line,thick=thick
    fsc_plot,plotRedshifts,mass_tot,color=fsc_color(colors[2]),$
             overplot=overplot,line=line,thick=thick

  endfor ;j

  ; legend
  legend,['gas in all galaxies','bound gas not in galaxies','total'],textcolors=colors,box=0,margin=0.25
  
  end_PS

end

; plotGalCatMassFunction(): 

pro plotGalCatMassFunction, sP=sP

  ; config
  redshifts = [2.0,0.0]
  
  minNumGasInGal = 100 ; 64, 256, 0=no minimum
  
  colors = ['forest green','slate blue','black']

  ; plot
  plotName = sP.plotPath+'gal_massfunction_'+str(sP.res)+'.eps'
  
  xrange = [8.0,11.6]
  yrange = [1,5e3]  
  
  if (str(res[0]) eq 'two') then res = [256,128]    
  if (str(res[0]) eq 'all') then res = [512,256,128]

  start_PS, plotName

  fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,$
           xtitle="log ( galaxy gas mass [M"+textoidl("_{sun}")+"] )",$
           ytitle="Count",/xs,/ys,ymargin=[4,3],/ylog
    
  if (n_elements(res) eq 1) then $
    fsc_text,0.06,0.94,/normal,str(res)+"^3",alignment=0.5
  
  ; loop over each resolution
  for j=0,n_elements(res)-1 do begin
    for m=0,n_elements(redshifts)-1 do begin  

      sP = simParams(res=res[j],run=run,redshift=redshifts[m])
      
      ; load constant mass factor (gadget)
      ;massfac = loadSnapshotSubset(sP.simPath,snapNum=0,partType='gas',field='mass')
      ;massfac = massfac[0]
      ;TODO FOR UNEQUAL MASSES
      stop 
      
      ; load galaxy catalog and select galaxies above particle count cut
      galcat = galaxyCat(sP=sP)
      
      ; enforce minimum number of particles cut and compute galaxy masses
      if (minNumGasInGal gt 0) then begin
        w = where(galcat.galaxyLen ge minNumGasInGal,count)
        
        galMasses = fltarr(count)
        
        for i=0,count-1 do $
          galMasses[i] = galcat.galaxyLen[w[i]] * massfac
      endif else begin
        galMasses = galcat.galaxyLen * massfac
      endelse
      
      galMasses = codeMassToLogMsun(galMasses)
      print,redshifts[m],minmax(galMasses)

       ; overplot successive resolutions
      overplot = 1 ;1,1,1
      line     = j ;0,1,2
      thick    = !p.thick + 1 - 2*(j gt 0) ;3,1,1
      
      bin = 0.2
      
      hist = histogram(galMasses,binsize=bin,locations=xpts,min=xrange[0],max=xrange[1])
      
      w = where(hist ne 0)
      
      fsc_plot,xpts[w]+bin/2.0,hist[w],overplot=overplot,line=line,$
               thick=thick,color=fsc_color(colors[m])
      
    endfor ;m
  endfor ;j
  
  ; legend
  labels = []
  for i=0,n_elements(redshifts)-1 do labels = [labels,"z="+string(redshifts[i],format='(f3.1)')]
  legend,labels,textcolors=colors[0:n_elements(redshifts)-1],box=0,/right
  
  end_PS
      
end

; plotGalCatRadii(): plot distribution of galaxy/groupmem radial distances

pro plotGalCatRadii, sP=sP

  ; config
  colors = ['forest green','slate blue','black']
  
  ; subhalo selection config
  sgSelect = 'pri' ;pri,sec,all
  parNorm  = 'sec' ;pri,sec (normalize r_vir by primary or secondary parent)
  
  ; get normalized r/r_vir
  gcRad = gcSubsetProp(sP=sP,select=sgSelect,parNorm=parNorm,/rVirNorm)
  
  ; plot
  xrange = [-0.1,2.0]
  yrange = [1e1,3e5]
  bin = 0.02
  
  title = str(sP.res) + textoidl("^3") + " z=" + string(sP.redshift,format='(f3.1)')
  xtitle = "r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}")
  ytitle = "Count"
  
  plotName = sP.plotPath + 'galcat_radii_'+str(res)+'_'+str(sP.snap)+'_'+sgSelect+'_'+$
             parNorm+'_'+str(minNumPart)+'.eps'
  
  start_PS, plotName
    w = where(gcRad.gal ge xrange[0] and gcRad.gal le xrange[1])
    plothist,gcRad.gal[w],bin=bin,xrange=xrange,yrange=yrange,xtitle=xtitle,ytitle=ytitle,/ylog,$
             color=fsc_color(colors[0]),title=title,/ys
             
    w = where(gcRad.gmem ge xrange[0] and gcRad.gmem le xrange[1])
    plothist,gcRad.gmem[w],bin=bin,color=fsc_color(colors[1]),/overplot
    
    rad_all = [gcRad.gmem,gcRad.gal]
    
    w = where(rad_all ge xrange[0] and rad_all le xrange[1])
    plothist,rad_all[w],bin=bin,color=fsc_color(colors[2]),line=1,/overplot
             
    ; legend
    legend,['gas in all galaxies','bound gas not in galaxies','total'],$
           textcolors=colors,box=0,margin=0.25,/right
             
  end_PS
end

; plotRhoTempGalCut(): plot the Torrey+ (2010) galaxy cat on the (rho,temp) plane at a given redshift

pro plotRhoTempGalCut, sP=sP
  
  ; config
  galcut_T   = 6.0
  galcut_rho = 0.25
  
  ; subhalo selection config
  sgSelect = 'pri' ; pri,sec,all
  parNorm  = 'sec' ; pri,sec (normalize r_vir by primary or secondary parent)

  ; load (dens,temp) and masses  
  gcDens = gcSubsetProp(sP=sP,select=sgSelect,/curDens)
  gcTemp = gcSubsetProp(sP=sP,select=sgSelect,/curTemp)

  ; 2d histo parameters                
  nbins   = 140*(res/128)
  rMinMax = [-2.0,8.0]
  tMinMax = [3.0,7.0]
  
  dens = [alog10(rhoRatioToCrit(gcDens.gal,redshift=redshift)),$
          alog10(rhoRatioToCrit(gcDens.gmem,redshift=redshift))]
  temp = [gcTemp.gal,gcTemp.gmem]

  ; calculate bin sizes
  binSizeRho  = (rMinMax[1]-rMinMax[0]) / nbins
  binSizeTemp = (tMinMax[1]-tMinMax[0]) / nbins
  
  if not keyword_set(mass) then begin
    ; hist_2d (no weighting)
    h2rt = hist_2d(dens,temp,bin1=binSizeRho,bin2=binSizeTemp,$
                   min1=rMinMax[0],min2=tMinMax[0],max1=rMinMax[1]+binSizeRho,max2=tMinMax[1]+binSizeTemp)
    ; h2rt /= float(max(h2rt)) ;norm s.t. colorbar gives fraction wrt max cell
  endif else begin
    ; hist2d (weighted)
    h2rt = hist2d(dens,temp,mass,$
                  binsize1=binsizeRho,binsize2=binSizeTemp,$
                  min1=rMinMax[0],min2=tMinMax[0],max1=rMinMax[1]+binSizeRho,max2=tMinMax[1]+binSizeTemp)
  endelse
  
  ; take log n_part (~mass_tot)
  b = 1.1
  
  w = where(h2rt ne 0)
  h2rt[w] = alog10(h2rt[w]) / alog10(b)  
  
  plotBase = "galcut_rhot_"+sgSelect+'_'+parNorm+'_'+str(minNumPart)
  plotName = sP.plotPath + plotBase + '_'+str(sP.res)+'_'+str(sP.snap)+'.eps'
  
  start_PS, plotName
  
    ; color table
    loadct, 33, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
      
    tvim,h2rt,pcharsize=!p.charsize,scale=1,clip=[0,100],$;,/c_map
         xtitle="log ("+textoidl("\rho / \rho_{crit,z="+string(redshift,format='(i1)')+"}")+")",ytitle="log (T [K])",$
         stitle="log"+textoidl("_{"+string(b,format='(f3.1)')+"}")+" ( N"+$
                textoidl("_{part}")+" ) or ~log ( M"+textoidl("_{tot}")+" )",$
         barwidth=0.5,lcharsize=!p.charsize-0.5,xrange=rMinMax,yrange=tMinMax
    
    ; scale Torrey+ (2011) galaxy cut to physical density
    scalefac = snapNumToRedshift(sP=sP,/time) ; time flag gives simulation time = scale factor
    a3inv = 1.0 / (scalefac*scalefac*scalefac)
    
    ; draw cut line
    xpts = [1e-9,1e-1] ; comoving density (code units)
    ypts = galcut_T + galcut_rho * alog10(xpts * a3inv) ; cut on physical density  
      
    fsc_plot,alog10(rhoRatioToCrit(xpts,redshift=redshift)),ypts,$
             line=0,thick=!p.thick+1.0,color=fsc_color('light gray'),/overplot
 
  end_PS
           
end

; plotGasOrigins(): plot temperature, entropy, and primary/secondary radius for all gas elements
;                   from the galaxy catalog at a target redshift as a function of time backwards

pro plotGasOrigins, sP=sP

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

; plotGasOriginsTracks():
  
pro plotGasOriginsTracks, sP=sP

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
  
pro plotVerticalSlices,rad_gal,rad_gmem,temp_gal,temp_gmem,plotName,xrange,yrange,xtitle

  ; config
  colors     = ['forest green','slate blue','black']
  
  binSizeRad  = 0.2
  numRadBins  = 6
  
  binSizeTemp = 0.05
  
  start_PS, plotName, xs=8, ys=6
  
  !p.multi = [0,3,2]
  
  !x.margin -= [4.0,1.5]
  !y.margin -= [0.5,0.5]
  
    for i=0,numRadBins-1 do begin
      fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,ytitle="",xtitle=xtitle               
    
      ; select in radial bin
      w = where(rad_gal ge i*binSizeRad and rad_gal lt (i+1)*binSizeRad,count)
      h1_gal  = histogram(temp_gal[w],binsize=binSizeTemp,$
                          min=xrange[0],max=xrange[1],locations=xpts)
                          
      w = where(rad_gmem ge i*binSizeRad and rad_gmem lt (i+1)*binSizeRad,count)                   
      h1_gmem = histogram(temp_gmem[w],binsize=binSizeTemp,$
                          min=xrange[0],max=xrange[1],locations=xpts)
      
      ; move xpts to bin centers and normalize histograms
      xpts += binSizeTemp/2.0
      
      h1_both = h1_gal + h1_gmem
      normfac = 1.0 / max(h1_both)
      
      h1_gal  *= normfac
      h1_gmem *= normfac
      h1_both *= normfac
      
      ; plot
      fsc_plot,xpts,h1_gal,color=fsc_color(colors[0]),thick=!p.thick+2,/overplot
      fsc_plot,xpts,h1_gmem,color=fsc_color(colors[1]),thick=!p.thick+2,/overplot
      fsc_plot,xpts,h1_both,color=fsc_color(colors[2]),line=1,thick=!p.thick+2,/overplot
      
      binStr = string(i*binSizeRad,format='(f3.1)') + "-" + string((i+1)*binSizeRad,format='(f3.1)')
      fsc_text,xrange[1]*0.92,yrange[1]*0.9,binStr,alignment=0.5,charsize=!p.charsize-0.5
    endfor
  
  !p.multi = 0
  
  end_PS
  
end

; plot2DRadHistos(): helper function used by the subsequent 4 plot routines
  
pro plot2DRadHistos,plotBase,sP,h2rt_gal,h2rt_gmem,xrange,yrange,ytitle,$
                    virTempRange=virTempRange,massBinStr=massBinStr

  plotName = sP.plotPath + plotBase + '_rad_both_'+str(sP.res)+'_'+str(sP.snap)+'.eps'
  
  h2rt_both = h2rt_gal + h2rt_gmem               
                        
  ; take log n_part (~mass_tot)
  b = 1.1
  
  w = where(h2rt_gal ne 0)
  h2rt_gal[w] = alog10(h2rt_gal[w]) / alog10(b)
  
  w = where(h2rt_gmem ne 0)
  h2rt_gmem[w] = alog10(h2rt_gmem[w]) / alog10(b)
 
  w = where(h2rt_both ne 0)
  h2rt_both[w] = alog10(h2rt_both[w]) / alog10(b) 
  
  redshift = snapNumToRedshift(sP=sP)
  
  start_PS, plotName
    
    loadct, 33, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
    
    tvim,h2rt_both,pcharsize=!p.charsize,scale=1,clip=[0,100],$;,/c_map
         xtitle="r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}"),$
         ytitle=ytitle,$
         stitle="log"+textoidl("_{"+string(b,format='(f3.1)')+"}")+" ( N"+$
                textoidl("_{part}")+" ) or ~log ( M"+textoidl("_{tot}")+" )",$
         title=str(sP.res)+textoidl("^3")+" z="+string(redshift,format='(f3.1)')+" both",$
         barwidth=0.5,lcharsize=!p.charsize-0.5,$
         xrange=xrange,yrange=yrange;,/rct
            
      if keyword_set(virTempRange) then begin
        fsc_text,xrange[1]*0.82,yrange[1]*0.96,massBinStr,alignment=0.5,color=fsc_color('white')
        fsc_text,xrange[1]*0.85,virTempRange[1]*1.02,"T"+textoidl("_{vir}"),alignment=0.5,color=fsc_color('orange')
        fsc_plot,xrange,[virTempRange[0],virTempRange[0]],line=0,color=fsc_color('orange'),/overplot
        fsc_plot,xrange,[virTempRange[1],virTempRange[1]],line=0,color=fsc_color('orange'),/overplot
      endif
             
  end_PS
  
  plotName = sP.plotPath + plotBase + '_rad_gal_'+str(sP.res)+'_'+str(sP.snap)+'.eps'
  
  start_PS, plotName
    
    loadct, 33, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
    
    tvim,h2rt_gal,pcharsize=!p.charsize,scale=1,clip=[0,100],$;,/c_map
         xtitle="r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}"),$
         ytitle=ytitle,$
         stitle="log"+textoidl("_{"+string(b,format='(f3.1)')+"}")+" ( N"+$
                textoidl("_{part}")+" ) or ~log ( M"+textoidl("_{tot}")+" )",$
         title=str(sP.res)+"^3 z="+string(redshift,format='(f3.1)')+" gal",$
         barwidth=0.5,lcharsize=!p.charsize-0.5,$
         xrange=xrange,yrange=yrange;,/rct
         
      if keyword_set(virTempRange) then begin
        fsc_text,xrange[1]*0.82,yrange[1]*0.96,massBinStr,alignment=0.5,color=fsc_color('white')
        fsc_text,xrange[1]*0.85,virTempRange[1]*1.02,"T"+textoidl("_{vir}"),alignment=0.5,color=fsc_color('orange')
        fsc_plot,xrange,[virTempRange[0],virTempRange[0]],line=0,color=fsc_color('orange'),/overplot
        fsc_plot,xrange,[virTempRange[1],virTempRange[1]],line=0,color=fsc_color('orange'),/overplot
      endif
      
  end_PS
  
  plotName = sP.plotPath + plotBase + '_rad_gmem_'+str(sP.res)+'_'+str(sP.snap)+'.eps'
  
  start_PS, plotName
    
    loadct, 33, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
    
    tvim,h2rt_gmem,pcharsize=!p.charsize,scale=1,clip=[0,100],$;,/c_map
         xtitle="r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}"),$
         ytitle=ytitle,$
         stitle="log"+textoidl("_{"+string(b,format='(f3.1)')+"}")+" ( N"+$
                 textoidl("_{part}")+" ) or ~log ( M"+textoidl("_{tot}")+" )",$
         title=str(sP.res)+"^3 z="+string(redshift,format='(f3.1)')+" gmem",$
         barwidth=0.5,lcharsize=!p.charsize-0.5,$
         xrange=xrange,yrange=yrange;,/rct
         
      if keyword_set(virTempRange) then begin
        fsc_text,xrange[1]*0.82,yrange[1]*0.96,massBinStr,alignment=0.5,color=fsc_color('white')
        fsc_text,xrange[1]*0.85,virTempRange[1]*1.02,"T"+textoidl("_{vir}"),alignment=0.5,color=fsc_color('orange')
        fsc_plot,xrange,[virTempRange[0],virTempRange[0]],line=0,color=fsc_color('orange'),/overplot
        fsc_plot,xrange,[virTempRange[1],virTempRange[1]],line=0,color=fsc_color('orange'),/overplot
      endif
      
  end_PS
  
end

; plotTempRad2DHisto(): plot 2d histogram of gas temperature (log Kelvin) as a function of r_gas/r_vir
;                       as well as vertical slices
;
; curTemp=1     : current temperature
; maxPastTemp=1 : maximum previous gas temperature
; tVirNorm=1    : normalize by the virial temperature of the parent halo

pro plotTempRad2DHisto, sP=sP

  ; config
  massBins = [0.0,100.0] ; no massbins
  ;massBins = [9.5,10.0,10.5,11.0,11.5,12.0]
  
  curTemp     = 1 ; current temperature vs radius
  maxPastTemp = 0 ; maximum previous temperature vs radius
  
  tVirNorm    = 0 ; normalize temperature by virial temp of parent halo
  
  ; subhalo selection config
  sgSelect = 'pri' ; pri,sec,all
  parNorm  = 'sec' ; pri,sec (normalize r_vir by primary or secondary parent)
  
  ; get normalized r/r_vir
  gcRad = gcSubsetProp(sP=sP,select=sgSelect,parNorm=parNorm,/rVirNorm)
  
  ; get current or maxPast temperature
  gcTemp = gcSubsetProp(sP=sP,select=sgSelect,curTemp=curTemp,maxPastTemp=maxPastTemp)
  
  if keyword_set(tVirNorm) then begin
    ; calculate temperatures of parents and take (non-log) ratio
    gcVirTemp = gcSubsetProp(sP=sP,select=sgSelect,/virTemp)
    
    gcTemp.gal  = 10.0^gcTemp.gal / gcVirTemp.gal
    gcTemp.gmem = 10.0^gcTemp.gmem / gcVirTemp.gmem
  endif
  
  ; calculate masses of parents
  parentMass = gcSubsetProp(sP=sP,select=sgSelect,/parMass)
  
  for j=0,n_elements(massBins)-2 do begin
  
    ; plot config
    xrange = [0.0,3.0]
    yrange = [4.0,7.0]
  
    ;binSizeRad  = 0.05 - 0.01*(res/128)
    ;binSizeTemp = 0.05 - 0.01*(res/128)  
    binSizeRad  = 0.04 / (sP.res/128)
    binSizeTemp = 0.04 / (sP.res/128)
    
    if (keyword_set(tVirNorm)) then begin
      yrange = [0.0,1.4]
      binSizeTemp = 0.014 / (res/128)
    endif
  
    ; select members of this parent mass bins  
    wGal  = where(parentMass.gal gt massBins[j] and parentMass.gal le massBins[j+1],count1)
    wGmem = where(parentMass.gmem gt massBins[j] and parentMass.gmem le massBins[j+1],count2)
    
    print,j,count1,count2
    
    temp_gal  = gcTemp.gal[wGal]
    temp_gmem = gcTemp.gmem[wGmem]
    
    rad_gal   = gcRad.gal[wGal]
    rad_gmem  = gcRad.gmem[wGmem]
  
    ; create 2d histo (separate for galaxy/group member and composite) t_cur/t_max ratio
    h2rt_gal  = hist_2d(rad_gal,temp_gal,bin1=binSizeRad,bin2=binSizeTemp,$
                        min1=xrange[0],min2=yrange[0],max1=xrange[1]+binSizeRad,max2=yrange[1]+binSizeTemp)  
    h2rt_gmem = hist_2d(rad_gmem,temp_gmem,bin1=binSizeRad,bin2=binSizeTemp,$
                        min1=xrange[0],min2=yrange[0],max1=xrange[1]+binSizeRad,max2=yrange[1]+binSizeTemp)
       
    ; strip off last empty row/col
    sz = size(h2rt_gal)
    h2rt_gal  = h2rt_gal[0:sz[1]-2,0:sz[2]-2]
    h2rt_gmem = h2rt_gmem[0:sz[1]-2,0:sz[2]-2]
    
    ; 2d histo plot config
    if (keyword_set(curTemp)) then begin
      plotBase = "temp_"+sgSelect+'_'+parNorm
      ytitle   = "log ( T"+textoidl("_{cur}")+" )"
    endif
    
    if (keyword_set(maxPastTemp)) then begin
      plotBase = "tmax_"+sgSelect+'_'+parNorm
      ytitle   = "log ( T"+textoidl("_{max}")+" )"
    endif
    
    if (keyword_set(tVirNorm)) then begin
      plotBase = "tmax_tvirNorm_"+sgSelect+'_'+parNorm
      ytitle   = "T"+textoidl("_{max,cur}")+" / T"+textoidl("_{vir}")
    endif
    
    ; extra config for mass bins
    massBinStr   = !NULL
    virTempRange = !NULL
    
    if (n_elements(massBins) gt 2) then begin
      plotBase = plotBase+"_mbin="+str(j)
      
      massBinStr = string(massBins[j],format='(f4.1)') + '-' + string(massBins[j+1],format='(f4.1)')
      
      massRangeCode = 10.0^[massBins[j],massBins[j+1]] / 1e10
      virTempRange = alog10( codeMassToVirTemp(massRangeCode,redshift) )
    endif
    
    ; plot 2d histo
    plot2DRadHistos,plotBase,sP,h2rt_gal,h2rt_gmem,xrange,yrange,ytitle,$
                    virTempRange=virTempRange,massBinStr=massBinStr    
    
    ; vertical slices plot config
    if (keyword_set(curTemp)) then begin
      plotName = sP.plotPath + 'temp_rad_slices_'+str(sP.res)+'_'+str(sP.snap)+'_'+$
                 sgSelect+'_'+parNorm+'.eps'    
      xtitle   = "log ( T"+textoidl("_{cur}")+ " ) [K]"
    endif
    
    if (keyword_set(maxPastTemp)) then begin
      plotName = sP.plotPath + 'tmax_rad_slices_'+str(sP.res)+'_'+str(sP.snap)+'_'+$
                 sgSelect+'_'+parNorm+'.eps'    
      xtitle = "log ( T"+textoidl("_{max}")+ " ) [K]"
    endif
    
    if (keyword_set(tVirNorm)) then begin
      plotName = sP.plotPath + 'tmax_TvirNorm_rad_slices_'+str(sP.res)+'_'+str(sP.snap)+'_'+$
                 sgSelect+'_'+parNorm+'.eps'
      xtitle = "T"+textoidl("_{max,cur}")+" / T"+textoidl("_{vir}")
    endif
    
    if (n_elements(massBins) gt 2) then $
      plotName = strmid(plotName,0,strlen(plotName)-4) + '.mbin='+str(j)+'.eps'
    
    xrange = yrange
    yrange = [0.0,1.05]
    
    ; plot vertical slices
    plotVerticalSlices,rad_gal,rad_gmem,temp_gal,temp_gmem,plotName,xrange,yrange,xtitle
    
  endfor ;j
  
end
