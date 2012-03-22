; cosmoPlotGalCat.pro
; gas accretion project - plots related to galaxy catalog
; dnelson mar.2012

; plotGalCatRedshift(): plot the total component masses in the catalog as a function of redshift

pro plotGalCatRedshift, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
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

; plotGalCatMassFunction(): plot the galcat mass functions at a few redshifts and resolutions

pro plotGalCatMassFunction, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs

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

; plotGalCatRadii(): plot distribution of galaxy/groupmem radial distances at one snapshot

pro plotGalCatRadii, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs

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
  
  compile_opt idl2, hidden, strictarr, strictarrsubs  
  
  ; config
  galcut_T   = 6.0
  galcut_rho = 0.25
  
  ; subhalo selection config
  sgSelect = 'pri' ; pri,sec,all
  parNorm  = 'sec' ; pri,sec (normalize r_vir by primary or secondary parent)

  ; load (dens,temp) and masses  
  gcDens = gcSubsetProp(sP=sP,select=sgSelect,/curSingleVal,singleValField='dens')
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