; cosmoPlot.pro
; gas accretion project - plots
; dnelson nov.2011

pro plotGalCatRedshift, res=res, run=run

  units = getUnits()
  sP    = simParams(res=res,run=run)
  
  ; config
  colors = ['forest green','slate blue','black']
  
  minNumGasInGal = 100 ; 64, 0=no minimum
  
  ; plot axes
  plotRedshifts = snapNumToRedshift(/all)
  
  xrange = [6.0,0.0]
  yrange = [11.5,13.5]
    
  ; plot
  plotName = sP.plotPath+'galmass_vs_redshift_'+str(res)+'.eps'
  
  if (str(res[0]) eq 'two') then res = [256,128]    
  if (str(res[0]) eq 'all') then res = [512,256,128]    
    
  start_PS, plotName

  fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,$
           xtitle="Redshift",ytitle="log ( total gas mass [M"+textoidl("_{sun}")+"] )",xs=9,/ys,ymargin=[4,3]
  universeage_axis,xrange,yrange
    
  if (n_elements(res) eq 1) then $
    fsc_text,0.06,0.94,/normal,str(res)+"^3",alignment=0.5
  
  ; loop over each resolution
  for j=0,n_elements(res)-1 do begin    
  
    sP = simParams(res=res[j],run=run)
  
    ; arrays
    mass_tot_gal  = fltarr(sP.snapRange[1])
    mass_tot_gmem = fltarr(sP.snapRange[1])
    
    ; load constant mass factor (gadget)
    massfac = loadSnapshotSubset(sP.simPath,snapNum=0,partType='gas',field='mass')
    massfac = massfac[0]
    
    ; load galaxy and group membership catalogs
    for m=sP.groupcatRange[0],sP.groupcatRange[1]-1 do begin
      gc = galaxyCat(res=res[j],run=run,snap=m)
      
      mass_tot_gal[m]  = n_elements(gc.galaxyIDs) * massfac
      mass_tot_gmem[m] = n_elements(gc.groupmemIDs) * massfac
      
      ; enforce minimum number of particles cut
      if (minNumGasInGal gt 0) then begin
        wGal = where(gc.galaxyLen ge minNumGasInGal,countGal)
        wGmem = where(gc.groupmemLen ge minNumGasInGal,countGmem)
        
        if (countGal ne 0) then $
          mass_tot_gal[m] -= (total(gc.galaxyLen[wGal]) * massfac)
        if (countGmem ne 0) then $
          mass_tot_gmem[m] -= (total(gc.groupmemLen[wGmem]) * massfac)
      endif
    endfor
    
    ; convert masses
    mass_tot      = codeMassToLogMsun(mass_tot_gal+mass_tot_gmem)
    mass_tot_gal  = codeMassToLogMsun(mass_tot_gal)
    mass_tot_gmem = codeMassToLogMsun(mass_tot_gmem)
    
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

pro plotGalCatMassFunction, res=res, run=run

  sP = simParams(res=res,run=run)
 
  ; config
  redshifts = [2.0,0.0]
  
  minNumGasInGal = 100 ; 64, 256, 0=no minimum
  
  colors = ['forest green','slate blue','black']

  ; plot
  plotName = sP.plotPath+'gal_massfunction_'+str(res)+'.eps'
  
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

      targetSnap = redshiftToSnapnum(redshifts[m])
      sP = simParams(res=res[j],run=run)
      
      ; load constant mass factor (gadget)
      massfac = loadSnapshotSubset(sP.simPath,snapNum=0,partType='gas',field='mass')
      massfac = massfac[0]      
      
      ; load galaxy catalog and select galaxies above particle count cut
      gc = galaxyCat(res=res[j],run=run,snap=targetSnap)
      
      ; enforce minimum number of particles cut and compute galaxy masses
      if (minNumGasInGal gt 0) then begin
        w = where(gc.galaxyLen ge minNumGasInGal,count)
        
        galMasses = fltarr(count)
        
        for i=0,count-1 do $
          galMasses[i] = gc.galaxyLen[w[i]] * massfac
      endif else begin
        galMasses = gc.galaxyLen * massfac
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

pro plotGalCatRadii, res=res, run=run

  sP = simParams(res=res,run=run)
  
  ; config
  redshift = 0.0
  
  colors = ['forest green','slate blue','black']
  
  ; subhalo selection config
  sgSelect = 'pri' ;pri,sec,all
  parNorm  = 'sec' ;pri,sec (normalize r_vir by primary or secondary parent)
  minNumPart = 0
  
  ; get normalized r/r_vir
  sP.snap = redshiftToSnapnum(redshift)
  
  gcRad = gcSubsetProp(sP=sP,sgSelect=sgSelect,parNorm=parNorm,minNumPart=minNumPart,/rVirNorm)
  
  ; plot
  xrange = [-0.1,2.0]
  yrange = [1e1,3e5]
  bin = 0.02
  
  title = str(res) + textoidl("^3") + " z=" + string(redshift,format='(f3.1)')
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

; compTmaxEffEOS(): compare maxiumum temperature reached by gas particles when including time spent
;                   on the effective equation of state vs not

pro compTmaxEffEOS, res=res, run=run

  sP = simParams(res=res,run=run)
  
  ; config
  redshift = 2.0
  
  ; load
  maxt_all = maxTemps(res=res,run=run,redshift=redshift,inclEffEOS=1)
  maxt_noeffeos = maxTemps(res=res,run=run,redshift=redshift,inclEffEOS=0)
  
  ; plot config
  xyrange = [4.0,7.0]
  
  sP.snap = redshiftToSnapnum(redshift)
  plotName = sP.plotPath + 'tmax_effeos_'+str(res)+'_'+str(sP.snap)+'.eps'
  
  start_PS, plotName
  
  fsc_plot,[0],[0],/nodata,xrange=xyrange,yrange=xyrange,$
           xtitle="log ( T"+textoidl("_{max}")+" excluding eff eos )",$
           ytitle="log ( T"+textoidl("_{max}")+" including all time )",/xs,/ys,$
           title=str(res)+textoidl("^3")+" z=30 to z="+string(redshift,format='(i1)')
  
  fsc_plot,maxt_noeffeos.maxTemps_gal,maxt_all.maxTemps_gal,psym=3,$
           color=fsc_color('black'),/overplot
           
  ;fsc_plot,[4.0,7.0],[4.0,7.0],color=fsc_color('orange'),/overplot
           
  end_PS  
  
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
                    blur=blur,virTempRange=virTempRange,massBinStr=massBinStr

  plotName = sP.plotPath + plotBase + '_rad_both_'+str(sP.res)+'_'+str(sP.snap)+'.eps'
  
  h2rt_both = h2rt_gal + h2rt_gmem               
                        
  ; take log n_part (~mass_tot)
  b = 2
  
  w = where(h2rt_gal ne 0)
  h2rt_gal[w] = alog10(h2rt_gal[w]) / alog10(b)
  
  w = where(h2rt_gmem ne 0)
  h2rt_gmem[w] = alog10(h2rt_gmem[w]) / alog10(b)
 
  w = where(h2rt_both ne 0)
  h2rt_both[w] = alog10(h2rt_both[w]) / alog10(b) 
  
  ;gaussian kernel convolution
  if keyword_set(blur) then begin
    h2rt_both = filter_image(h2rt_both,FWHM=[1.2,1.2],/ALL)
    ;h2rt_both = smooth(h2rt_both,3)
    ;gal
    ;gmem
  endif
  
  redshift = snapNumToRedshift(snap=sP.snap)
  
  start_PS, plotName
    
    loadct, 33, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
    
    tvim,h2rt_both,pcharsize=!p.charsize,scale=1,clip=[0,100],$;,/c_map
         xtitle="r"+textoidl("_{gas}")+" / r"+textoidl("_{vir}"),$
         ytitle=ytitle,$
         stitle="log"+textoidl("_{"+str(b)+"}")+" ( N"+textoidl("_{part}")+" ) or ~log ( M"+textoidl("_{tot}")+" )",$
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
         stitle="log"+textoidl("_{"+str(b)+"}")+" ( N"+textoidl("_{part}")+" ) or ~log ( M"+textoidl("_{tot}")+" )",$
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
         stitle="log"+textoidl("_{"+str(b)+"}")+" ( N"+textoidl("_{part}")+" ) or ~log ( M"+textoidl("_{tot}")+" )",$
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

pro plotTempRad2DHisto, res=res, run=run

  sP = simParams(res=res,run=run)
  
  ; config
  redshift   = 2.0
  
  massBins = [0.0,100.0] ;[9.5,10.0,10.5,11.0,11.5,12.0]
  
  curTemp     = 0 ; current temperature vs radius
  maxPastTemp = 1 ; maximum previous temperature vs radius
  
  tVirNorm    = 1 ; normalize temperature by virial temp of parent halo
  
  ; subhalo selection config
  sgSelect = 'pri' ; pri,sec,all
  parNorm  = 'sec' ; pri,sec (normalize r_vir by primary or secondary parent)
  minNumPart = 1   ; enforce minimum 100 gas particles per halo
  
  ; get normalized r/r_vir
  sP.snap = redshiftToSnapnum(redshift)
  
  gcRad = gcSubsetProp(sP=sP,sgSelect=sgSelect,minNumPart=minNumPart,parNorm=parNorm,/rVirNorm)
  
  ; get current or maxPast temperature
  gcTemp = gcSubsetProp(sP=sP,sgSelect=sgSelect,minNumPart=minNumPart,$
                        curTemp=curTemp,maxPastTemp=maxPastTemp)
  
  if keyword_set(tVirNorm) then begin
    ; calculate temperatures of parents and take (non-log) ratio
    gcVirTemp = gcSubsetProp(sP=sP,sgSelect=sgSelect,minNumPart=minNumPart,/virTemp)
    
    gcTemp.gal  = 10.0^gcTemp.gal / gcVirTemp.gal
    gcTemp.gmem = 10.0^gcTemp.gmem / gcVirTemp.gmem
  endif
  
  ; calculate masses of parents
  parentMass = gcSubsetProp(sP=sP,sgSelect=sgSelect,minNumPart=minNumPart,/parMass)
  
  for j=0,n_elements(massBins)-2 do begin
  
    ; plot config
    xrange = [0.0,3.0]
    yrange = [3.5,7.0]
  
    ;binSizeRad  = 0.05 - 0.01*(res/128)
    ;binSizeTemp = 0.05 - 0.01*(res/128)  
    binSizeRad  = 0.04 / (res/128)
    binSizeTemp = 0.04 / (res/128)
    
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

    gcRad  = !NULL
    gcTemp = !NULL
  
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
      plotBase = "temp"
      ytitle   = "log ( T"+textoidl("_{cur}")+" )"
    endif
    
    if (keyword_set(maxPastTemp)) then begin
      plotBase = "tmax"
      ytitle   = "log ( T"+textoidl("_{max}")+" )"
    endif
    
    if (keyword_set(tVirNorm)) then begin
      plotBase = "tmax_tvirNorm"
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
      plotName = sP.plotPath + 'temp_rad_slices_'+str(res)+'_'+str(sP.snap)+'.eps'
      xtitle   = "log ( T"+textoidl("_{cur}")+ " ) [K]"
    endif
    
    if (keyword_set(maxPastTemp)) then begin
      plotName = sP.plotPath + 'tmax_rad_slices_'+str(res)+'_'+str(sP.snap)+'.eps'  
      xtitle = "log ( T"+textoidl("_{max}")+ " ) [K]"
    endif
    
    if (keyword_set(tVirNorm)) then begin
      plotName = sP.plotPath + 'tmax_TvirNorm_rad_slices_'+str(res)+'_'+str(sP.snap)+'.eps'
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
