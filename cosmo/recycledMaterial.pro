; recycledMaterial.pro
; feedback project - analysis of recycled material
; dnelson jul.2013

; windCountHisto(): 1d/2d histograms of the trMC windcounter

pro windCountHisto
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  message,'these results seem wierd, at least at 256, check'
  
  ; config
  sP = simParams(res=256,run='feedback',redshift=2.0)
  
  xrange = [9.5,12.5]
  yrange = [-0.5,8.5] ; center on integers, start at zero
  
  binSizeMass = 0.25/2

  hsp = [0.007,0.08]*2 ; mass (0.003 10-12), frac
  nc  = 200 ; number of colors (of 255) to use for background 2d histo
  
  ; load
  parMass = gcSubsetProp(sP=sP,/parMass)
  windc   = gcSubsetProp(sP=sP,/curTracerVal,curField='tracer_windcounter')
  maxt    = gcSubsetProp(sP=sP,/maxPastTemp)

  ; bin fractions into halo mass bins and make median lines
  logMassBins   = [9.5,10.0,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11.0,$
                   11.1,11.25,11.5,11.75,11.9,13.1]
  logMassNBins  = n_elements(logMassBins)-1
  logMassBinCen = 0.5 * (logMassBins + shift(logMassBins,-1))
  logMassBinCen = logMassBinCen[0:-2]
  
  medFracs = { galaxy_zero  : fltarr(logMassNbins) ,$
               halo_zero    : fltarr(logMassNbins) ,$
               galaxy_1plus : fltarr(logMassNbins) ,$
               halo_1plus   : fltarr(logMassNbins) ,$
               galaxy_2plus : fltarr(logMassNbins) ,$
               halo_2plus   : fltarr(logMassNbins)  }
               
  ; calculate median values in bins of halo mass
  for i=0,logMassNbins-1 do begin
   
    galaxyTot = 0
    haloTot   = 0
  
    for j=0,n_tags(parMass)-1 do begin
      ; look in halo mass bin for this type
      w_type = where(parMass.(j) gt logMassBins[i] and parMass.(j) le logMassBins[i+1],count_type)
      
      if count_type eq 0 then continue
      
      ; count windcounter levels
      w = where(windc.(j)[w_type] eq 0,count_zero)
      w = where(windc.(j)[w_type] ge 1,count_1plus)
      w = where(windc.(j)[w_type] ge 2,count_2plus)
      
      ; add to galaxy?
      typeLabel = (tag_names(parMass))[j]
      
      if typeLabel eq 'GAL' or typeLabel eq 'STARS' or typeLabel eq 'BHS' then begin
        medFracs.galaxy_zero[i]  += count_zero
        medFracs.galaxy_1plus[i] += count_1plus
        medFracs.galaxy_2plus[i] += count_2plus
        galaxyTot += count_type
      endif
      
      ; add to halo?
      if typeLabel eq 'GMEM' then begin
        medFracs.halo_zero[i]  += count_zero
        medFracs.halo_1plus[i] += count_1plus
        medFracs.halo_2plus[i] += count_2plus
        haloTot += count_type
      endif
      
    endfor
      
    ; normalize (->mass fractions)
    medFracs.galaxy_zero[i]  /= galaxyTot
    medFracs.galaxy_1plus[i] /= galaxyTot
    medFracs.galaxy_2plus[i] /= galaxyTot
    
    medFracs.halo_zero[i]  /= haloTot
    medFracs.halo_1plus[i] /= haloTot
    medFracs.halo_2plus[i] /= haloTot
  endfor
  
  ; plot (1) - fractions vs halo mass
  colors = ['black','red','blue','forest green','orange']
  plotStr = '.' + str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)
  
  start_PS, sP.plotPath + 'recFracs' + plotStr + '.eps', /big
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=[0.0,1.0],/xs,/ys,$
      xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),ytitle="Mass Fraction"
      
    cgPlot,logMassBinCen,medFracs.galaxy_zero,line=0,color=cgColor(colors[0]),/overplot
    cgPlot,logMassBinCen,medFracs.galaxy_1plus,line=0,color=cgColor(colors[1]),/overplot
    cgPlot,logMassBinCen,medFracs.galaxy_2plus,line=0,color=cgColor(colors[2]),/overplot
    
    cgPlot,logMassBinCen,medFracs.halo_zero,line=1,color=cgColor(colors[0]),/overplot
    cgPlot,logMassBinCen,medFracs.halo_1plus,line=1,color=cgColor(colors[1]),/overplot
    cgPlot,logMassBinCen,medFracs.halo_2plus,line=1,color=cgColor(colors[2]),/overplot
     
    ; legend
    legend,['zero','1+','2+'],textcolors=colors,/top,/left
    legend,['galaxy','halo'],linestyle=[0,1],/top,/right
     
  end_PS
  
  ; plot (2) - 2d histo vs halo mass
  start_PS, sP.plotPath + 'rec2DHisto' + plotStr + '.eps', /big
  
    ; both (gal+stars)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange-[0.05,0.0],/xs,/ys,$
      xtitle="",ytitle="Wind Counter",xtickname=replicate(' ',10),$
      yticks=ceil(max(yrange))+1,ytickv=indgen(ceil(max(yrange))),yminor=0,pos=(sP.pos_2x1)[0]
      
    xx = [parMass.gal,parMass.stars]
    yy = [windc.gal,windc.stars]
      
    f2d = binHisto2D(xx=xx, yy=yy, xmm=xrange, ymm=yrange, xbs=binSizeMass, ybs=1.0)
    oplot2DHistoSq, f2d, hsp=hsp, nc=nc, xrange=xrange, yrange=yrange, /colNorm, /blue
    
    cgText,xrange[0]*1.01,yrange[1]*0.84,"galaxy",alignment=0.0
    
    ; redo plot borders
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange-[0.05,0.0],/xs,/ys,xtitle="",ytitle="",/noerase,$
      yticks=ceil(max(yrange))+1,ytickv=indgen(ceil(max(yrange))),yminor=0,pos=(sP.pos_2x1)[0],$
      xtickname=replicate(' ',10)
    
    ; gmem
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange-[0.05,0.0],/xs,/ys,$
      xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),ytitle="Wind Counter",$
      yticks=ceil(max(yrange))+1,ytickv=indgen(ceil(max(yrange))),yminor=0,pos=(sP.pos_2x1)[1],/noerase
      
    f2d = binHisto2D(xx=parMass.gmem, yy=windc.gmem, xmm=xrange, ymm=yrange, xbs=binSizeMass, ybs=1.0)
    oplot2DHistoSq, f2d, hsp=hsp, nc=nc, xrange=xrange, yrange=yrange, /colNorm, /blue
    
    cgText,xrange[0]*1.01,yrange[1]*0.84,"halo",alignment=0.0
    
    ; redo plot borders
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange-[0.05,0.0],/xs,/ys,xtitle="",ytitle="",/noerase,$
      yticks=ceil(max(yrange))+1,ytickv=indgen(ceil(max(yrange))),yminor=0,pos=(sP.pos_2x1)[1]
    
  end_PS
  
  ; plot (3) - global histogram
  start_PS, sP.plotPath + 'windcountHisto' + plotStr + '.eps', /big
  
    cgPlot,[0],[0],/nodata,xrange=[-0.2,5.2],yrange=[1e-3,1.0],/xs,/ys,/ylog,yminor=0,$
      xtitle=textoidl("windcounter"),ytitle="Fraction"
      
    ; mass bins
    massBins = list( [8.0,13.0], [9.5,10.0], [10.5,11.0], [11.5,12.0] )
    massBinStrs = ['all']
    
    foreach massBin,massBins,i do begin
      w_gal   = where(parMass.gal ge massBin[0]   and parMass.gal lt massBin[1],  count_gal)
      w_stars = where(parMass.stars ge massBin[0] and parMass.stars lt massBin[1],count_stars)
      w_gmem  = where(parMass.gmem ge massBin[0]  and parMass.gmem lt massBin[1], count_gmem)
      
      hh_both = histogram( [windc.gal[w_gal], windc.stars[w_stars]], loc=loc1 )
      if count_gmem gt 0 then hh_gmem = histogram( windc.gmem[w_gmem], loc=loc2 )
      
      cgPlot,loc1,float(hh_both)/total(hh_both),line=0,psym=-4,color=colors[i],/overplot
      if count_gmem gt 0 then cgPlot,loc2,float(hh_gmem)/total(hh_gmem),line=1,psym=-4,color=colors[i],/overplot
      
      if i gt 0 then $
      massBinStrs = [massBinStrs,string(massBin[0],format='(f4.1)')+textoidl(' < log (_{ }M_{h }) < ')+ $
                                 string(massBin[1],format='(f4.1)')]
    endforeach
    
    ; legend
    legend,massBinStrs,textcolor=colors[0:n_elements(massBins)],/top,/right,box=0
  
  end_PS
  
  ; plot (4) - histogram by halo mass bin
  

  
  stop
  
end