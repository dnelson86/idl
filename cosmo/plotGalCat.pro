; plotGalCat.pro
; gas accretion project - plots related to galaxy catalog
; dnelson jul.2013

; plotGalCatRedshift(): plot the total component masses in the catalog as a function of redshift

pro plotGalCatRedshift, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  res    = [128]
  
  ; plot config
  plotRedshifts = snapNumToRedshift(sP=sP,/all)
  
  xrange = [6.0,0.0]
  yrange = [11.5,13.5]
    
  ; plot
  start_PS, sP.plotPath+'galmass_vs_redshift_r'+str(n_elements(sP.res))+'.eps'

  cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,$
           xtitle="Redshift",ytitle="log ( total mass [M"+textoidl("_{sun}")+"] )",$
           xs=9,/ys,ymargin=[4,3]
  universeage_axis,xrange,yrange
    
  if (n_elements(res) eq 1) then $
    cgText,0.06,0.94,/normal,str(res)+"^3",alignment=0.5
  
  ; loop over each resolution
  for j=0,n_elements(res)-1 do begin    
  
    sP = simParams(res=res[j],run=sP.run,redshift=sP.redshift)
    
    ; mass weighting
    if sP.trMCPerCell gt 0 then massWt = sP.trMassConst * units.UnitMass_in_Msun
    if sP.trMCPerCell eq 0 then massWt = sP.targetGasMass * units.UnitMass_in_Msun
    if sP.trMCPerCell lt 0 then message,'error'     
  
    galcat = galaxyCat(sP=sP, /skipSave)
  
    ; arrays for this resolution
    for i=0,n_tags(galcat.types)-1 do $
      mass_tot = mod_struct( mass_tot, (tag_names(galcat.types))[i], fltarr(sP.snapRange[1]) )
    mass_tot = mod_struct( mass_tot, 'TOTAL', fltarr(sP.snapRange[1]) )
    
    ; load galaxy and group membership catalogs
    for m=sP.groupcatRange[0],sP.groupcatRange[1]-1 do begin
      sP.snap = m
      print,m
      galcat = galaxyCat(sP=sP, /skipSave)
      
      for i=0,n_tags(galcat.types)-1 do begin
        w = where(galcat.type eq galcat.types.(i), count)
        if sP.trMCPerCell eq 0 then mass_tot.(i)[m] = count * massWt
        if sP.trMCPerCell gt 0 then mass_tot.(i)[m] = total(galcat.trMC_cc[w]) * massWt
        if sP.trMCPerCell lt 0 then mass_tot.(i)[m] = total(gaclat.trVel_cc[w]) * massWt
      endfor
      
      ; sum for total
      for i=0,n_tags(galcat.types)-1 do $
        mass_tot.total[m] += mass_tot.(i)[m]
      
    endfor
    
    ; convert masses to log msun
    for i=0,n_tags(mass_tot)-1 do $
      mass_tot.(i) = codeMassToLogMsun(mass_tot.(i))
    
     ; overplot successive resolutions
    line     = j ;0,1,2
    thick    = !p.thick + 1 - 2*(j gt 0) ;3,1,1

    for i=0,n_tags(mass_tot)-1 do $
      cgPlot,plotRedshifts,mass_tot.(i),color=cgColor(units.colors[j]),line=line,thick=thick,/overplot

  endfor ;j

  ; legend
  legend,str(res),textcolors=units.colors[0:n_elements(res)-1],/top,/left
  
  end_PS

end

; plotGalCatMassFunction(): plot the galcat mass functions at a few redshifts and resolutions

pro plotGalCatMassFunction, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  redshifts = [2.0,0.0]
  res       = [256,128]

  ; plot config
  xrange = [8.0,11.6]
  yrange = [1,5e3]  
  binSize = 0.2

  ; plot
  start_PS, sP.plotPath+'gal_massfunction_'+str(sP.res)+'.eps'

  cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,$
           xtitle="log ( galaxy gas mass [M"+textoidl("_{sun}")+"] )",$
           ytitle="Count",/xs,/ys,ymargin=[4,3],/ylog
    
  if (n_elements(res) eq 1) then $
    cgText,0.06,0.94,/normal,str(res)+"^3",alignment=0.5
  
  ; loop over each resolution
  for j=0,n_elements(res)-1 do begin
    for m=0,n_elements(redshifts)-1 do begin  

      sP = simParams(res=res[j],run=sP.run,redshift=redshifts[m])
      
      ; mass weighting
      if sP.trMCPerCell gt 0 then massWt = sP.trMassConst * units.UnitMass_in_Msun
      if sP.trMCPerCell eq 0 then massWt = sP.targetGasMass * units.UnitMass_in_Msun
      if sP.trMCPerCell lt 0 then message,'error'    
      
      ; load galaxy catalog and select galaxies above particle count cut
      galcat = galaxyCat(sP=sP, /skipSave)
      
      w = where(galcat.type eq galcat.types.gal, count)
      if sP.trMCPerCell eq 0 then galMasses = count * massWt
      if sP.trMCPerCell gt 0 then galMasses = total(galcat.trMC_cc[w]) * massWt
      if sP.trMCPerCell lt 0 then galMasses = total(gaclat.trVel_cc[w]) * massWt
      
      galMasses = codeMassToLogMsun(galMasses)
      print,redshifts[m],minmax(galMasses)

       ; overplot successive resolutions
      line     = j ;0,1,2
      thick    = !p.thick + 1 - 2*(j gt 0) ;3,1,1

      hist = histogram(galMasses,binsize=binSize,locations=xpts,min=xrange[0],max=xrange[1])
      
      w = where(hist ne 0)
      
      cgPlot,xpts[w]+bin/2.0,hist[w],line=line,thick=thick,color=cgColor(units.colors[m]),/overplot
      
    endfor ;m
  endfor ;j
  
  ; legend
  labels = "z=" + string(redshifts,format='(f3.1)')
  legend,labels,textcolors=units.colors[0:n_elements(redshifts)-1],/top,/right
  
  end_PS
      
end

; plotGalCatRadii(): plot distribution of galaxy/groupmem radial distances at one snapshot

pro plotGalCatRadii, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; get normalized r/r_vir
  gcRad = gcSubsetProp(sP=sP,s/rVirNorm)
  
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
             color=cgColor(units.colors[0]),title=title,/ys
             
    w = where(gcRad.gmem ge xrange[0] and gcRad.gmem le xrange[1])
    plothist,gcRad.gmem[w],bin=bin,color=cgColor(units.colors[1]),/overplot
    
    rad_all = [gcRad.gmem,gcRad.gal]
    
    w = where(rad_all ge xrange[0] and rad_all le xrange[1])
    plothist,rad_all[w],bin=bin,color=cgColor(units.colors[2]),line=1,/overplot
             
    ; legend
    legend,['gas in all galaxies','bound gas not in galaxies','total'],$
           textcolors=units.colors[0:2],box=0,margin=0.25,/right
             
  end_PS
end

; plotRhoTempGalCut(): plot the Torrey+ (2010) galaxy cat on the (rho,temp) plane at a given redshift

pro plotRhoTempGalCut, sP=sP
  
  compile_opt idl2, hidden, strictarr, strictarrsubs  
  
  ; config
  galcut_T   = 6.0
  galcut_rho = 0.25
  
  ; load (dens,temp) and masses  
  gcDens = gcSubsetProp(sP=sP,/curGasVal,curField='dens')
  gcTemp = gcSubsetProp(sP=sP,/curGasVal,curField='temp')

  ; 2d histo parameters                
  nbins   = 140*(res/128)
  rMinMax = [-2.0,8.0]
  tMinMax = [3.0,7.0]
  
  dens = [alog10(rhoRatioToCrit(gcDens.gal,sP=sP)),$
          alog10(rhoRatioToCrit(gcDens.gmem,sP=sP))]
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
  
  plotBase = "galcut_rhot_"+str(minNumPart)
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
      
    cgPlot,alog10(rhoRatioToCrit(xpts,sP=sP)),ypts,$
           line=0,thick=!p.thick+1.0,color=cgColor('light gray'),/overplot
 
  end_PS
           
end