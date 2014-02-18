; plotAccTimes.pro
; gas accretion project - past radial history of gas elements (plotting)
; dnelson feb.2014

; helper functions:

function selectAccTimeDelta, sP=sP, at=at, galcat=galcat, $
                             earlierInd=earlierInd, laterInd=laterInd, $
                             gcType=gcType, norm=norm

  units = getUnits()
  
  ; all times are interpolated, so if exactly the same exclude with strict LT
  if gcType eq 'all' then $
    w = where(at.accTime[earlierInd,*] ge 0 and at.accTime[laterInd,*] ge 0 and $
              at.accTime[earlierInd,*] lt at.accTime[laterInd,*],count)
  if gcType eq 'gal' then $
    w = where(at.accTime[earlierInd,*] ge 0 and at.accTime[laterInd,*] ge 0 and $
              at.accTime[earlierInd,*] lt at.accTime[laterInd,*] and $
              (types eq galcat.types.gal or types eq galcat.types.stars),count)
  if gcType eq 'gmem' then $
    w = where(at.accTime[earlierInd,*] ge 0 and at.accTime[laterInd,*] ge 0 and $
              at.accTime[earlierInd,*] lt at.accTime[laterInd,*] and $
              types eq galcat.types.gmem,count)
  if gcType eq 'inter' then $
    w = where(at.accTime[earlierInd,*] ge 0 and at.accTime[laterInd,*] ge 0 and $
              at.accTime[earlierInd,*] lt at.accTime[laterInd,*] and $
              types eq galcat.types.inter,count)
              
  age_earlier = reform( redshiftToAgeFlat( 1.0/at.accTime[earlierInd,w]-1.0 ) )
  age_later   = reform( redshiftToAgeFlat( 1.0/at.accTime[laterInd,w]-1.0 ) )
  
  val = age_later - age_earlier
  
  if count eq 0 then message,'Error'
  
  if keyword_set(norm) then val /= age_earlier
  
  return, {w:w,val:val}
end

function colorMapAccTime, h2, logHist=logHist, byRow=byRow
  ; config
  min = 10.0
  max = 255.0
  
  h2_cmap = h2
  
  ; take values to log and lift all zeros to minimum nonzero
  if keyword_set(logHist) then begin
    w = where(h2_cmap gt 0,count)
    if count gt 0 then h2_cmap[w] = alog10(h2_cmap[w])
  endif
  
  w = where(h2_cmap eq 0,count)
  if count gt 0 then h2_cmap[w] = min(h2_cmap[where(h2_cmap gt 0)])
  print,minmax(h2_cmap)
    
  ; normalization by row?
  if keyword_set(byRow) then begin
    nRows = n_elements(h2_cmap[0,*])
    
    for i=0,nRows-1 do begin
      curRow = reform( h2_cmap[*,i] )
      curRow = fix( (curRow-min(curRow))/(max(curRow)-0.0) * (max-min) + min )
      h2_cmap[*,i] = curRow
    endfor
  endif else begin
    ; linear stretch from 10 (non-white) to 255
    h2_cmap = fix( (h2_cmap-min(h2_cmap))/(max(h2_cmap)-0.0) * (max-min) + min )
  endelse
  
    
  return, h2_cmap
end

function colorMapAccTimeDelta, h2a, h2b, fracMin=fracMin, fracMax=fracMax
  ; config
  min = 10.0
  max = 255.0
  
  if n_elements(fracMin) eq 0 or n_elements(fracMax) eq 0 then message,'Error'
  
  h2d_cmap = h2a*0.0 + 255/2 ; center missing values at white
 
  ; row by row
  nRows = n_elements(h2d_cmap[0,*])

  for i=0,nRows-1 do begin
    ; extract
    aRow = reform( h2a[*,i] )
    bRow = reform( h2b[*,i] )
    
    ; normalize each to their total (now a fraction)
    aTotal = total(aRow,/int)
    bTotal = total(bRow,/int)
    
    if aTotal eq 0 or bTotal eq 0 then continue ; leave output as zero diff
    aRow /= aTotal
    bRow /= bTotal
    
    ; stretch fractional difference between preset bounds
    dRow = alog10( aRow / bRow )
    ;dRow = (aRow - bRow) ;/ bRow
    dRow = (dRow-fracMin)/(fracMax-fracMin) * (max-min) + min
    ;dRow = fix( (dRow-min(dRow))/max(dRow) * (max-min) + min )
    
    ; set any non-finite (e.g. a=0 || b=0) to center value
    w = where( ~finite(dRow), count )
    if count gt 0 then dRow[w] = 255/2
    
    h2d_cmap[*,i] = dRow
  endfor
  
  ; convert to byte, clamp to [0,255]
  dRow = fix(dRow)
  dRow = dRow > 0 < 255
  return, h2d_cmap
end

; plotAccTimeDeltas():

pro plotAccTimeDeltas

  ; config
  redshifts = [2.0] ;[3.0,2.0,1.0,0.0]
  runs      = ['feedback','tracer','gadget']
  res       = 128
  
  ; index selection for difference ([1.0,0.75,0.5,0.25,0.15,0.05,0.01,first0.15,first1.0])
  rVirFacs = ['1.0 rvir','0.75 rvir','0.5 rvir','0.25 rvir','0.15 rvir','0.05 rvir',$
              '0.01 rvir','first 0.15 rvir','first 1.0 rvir']
  earlierInd = -1 ; -1  ; 0
  laterInd   = -2 ; -2  ; 4
  gcType     = 'all' ; all, gal (includes stars), gmem, inter
  
  ; plot config
  yrange     = [1e-4,0.2]
  xrange     = [0.0,5.0]
  xrangeNorm = [7e-3,1.5]
  nBins      = 100
  colorInd   = 1
  
  sP = simParams(res=res,run=runs[0],redshift=redshifts[0])
  
  plotName = 'accTimeDeltas_'+str(res)+'_'+str(earlierInd)+'_'+str(laterInd)+'_'+gcType+'.eps'
  start_PS,sP.plotPath + plotName, xs=3*3.5, ys=4*3.5
  
    pos = plot_pos(rows=4,cols=2,/gap)
    colors = []
    pinfo = {}
    
    cgText,0.5,0.98,textoidl("rVirFacs: t_1 = " + str(rVirFacs[laterInd]) + $
      ", t_2 = " + str(rVirFacs[earlierInd])) + ' (type='+gcType+')',alignment=0.5,/normal
    
    foreach redshift,redshifts,i do begin
    
    xrangeNorm = [7e-3,1.5 + (3-redshift)*1]
    xrange = [0.0,0.7] * redshiftToAgeFlat(redshift)
    
    ; plot (1)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,title="",$
      xtitle=textoidl("( t_{1} - t_{2} ) [Gyr]"),$
      ytitle="Fraction",/xs,/ys,pos=pos[i*2],/noerase,/ylog,yminor=0
    p1 = !P & x1 = !X & y1 = !Y
      
    ; plot (2)
    cgPlot,[0],[0],/nodata,xrange=xrangeNorm,yrange=yrange,title="",$
      xtitle=textoidl("( t_{1} - t_{2} ) / t_{2}"),$
      ytitle="Fraction",/xs,/ys,pos=pos[i*2+1],/noerase,/ylog,yminor=0,/xlog,xminor=0
    p2 = !P & x2 = !X & y2 = !Y
      
	; plot histogram for each run/redshift combination
      foreach run,runs,j do begin
        print, '['+str(j)+'] redshift = '+string(redshift,format='(f4.1)')+' run = '+run
        
        sP_cur = simParams(res=res,run=run,redshift=redshift)
        at     = accretionTimes(sP=sP_cur)
        galcat = galaxyCat(sP=sP_cur)
        
        if sP_cur.trMCPerCell gt 0 then types = ( galcat.type[ replicate_var(galcat.trMC_cc) ] )
        if sP_cur.trMCPerCell eq 0 then types = ( galcat.type )
        if sP_cur.trMCPerCell lt 0 then types = ( galcat.type[ replicate_var(galcat.trVel_cc) ] )
        
        ; add to plot (1)
        delta = selectAccTimeDelta(sP=sP, at=at, galcat=galcat, earlierInd=earlierInd, $
                                 laterInd=laterInd, gcType=gcType)

        !P = p1 & !X = x1 & !Y = y1
        
        ;DEBUG
        maxTemp = gcSubsetProp(sP=sP,/maxPastTemp,/accretionTimeSubset,accMode='all')
        delta2 = gcSubsetProp(sP=sP,accTimeDelta=[earlierInd,laterInd],/accretionTimeSubset,accMode='all')
        stop
        
        ;DEBUG (match to fluidTracks)
        if 0 then begin
          ww = where(delta.val lt -0.4 and delta.val gt -0.5)
          tInd = w[ww[1]]
        
          tracks = tracksFluid(sP=sP_cur)
          mt     = mergerTreeSubset(sP=sP_cur)
        
          gcIDList = mt.gcIndOrigTrMC[tInd]
          print,tInd,gcIDList
        
          haloTvir_t = reverse( mt.hVirTemp[*,gcIDList] )
          haloRvir_t = reverse( mt.hVirRad[*,gcIDList] )
          halo_ages  = reverse( redshiftToAgeFlat(1/mt.times-1) )
          halo_times = reverse( mt.times )
        
          rad_track = tracks.rad[*,tInd] / haloRvir_t
        
          print,rad_track
          print,at.accTime[*,tInd]
          stop
        endif
        ;END DEBUG
        
        binSize = (xrange[1] - xrange[0]) / float(nBins)
	  h = histogram(delta.val, bin=binSize, loc=loc, min=xrange[0], max=xrange[1])
	  h = float(h)/total(h)
        
        colors = [colors,sP_cur.colors[colorInd]]
	  cgPlot,loc+binSize*0.5,h,color=colors[-1],/overplot
        
        ; add to plot (2)
        delta = selectAccTimeDelta(sP=sP, at=at, galcat=galcat, earlierInd=earlierInd, $
                                 laterInd=laterInd, gcType=gcType, /norm)
                                 
        !P = p2 & !X = x2 & !Y = y2
        
        xrNormLog = alog10(xrangeNorm)
        binSize = (xrNormLog[1] - xrNormLog[0]) / float(nBins)
	  h = histogram(alog10(delta.val), bin=binSize, loc=loc, min=xrNormLog[0], max=xrNormLog[1])
	  h = float(h)/total(h)
        
	  cgPlot,10.0^(loc+binSize*0.5),h,color=colors[j],/overplot
	endforeach ; runs
      
      ; legends
      !P = p1 & !X = x1 & !Y = y1
      legend,['z = '+string(redshift,format='(f3.1)'),runs],textcolors=[0L,colors],/top,/right 
      !P = p2 & !X = x2 & !Y = y2
      legend,['z = '+string(redshift,format='(f3.1)'),runs],textcolors=[0L,colors],/top,/left
      
    endforeach ;redshifts
      
  end_PS
  stop
end

; plotAccTimeDeltaVsValMax()

pro plotAccTimeDeltaVsValMax

  ; config
  redshifts = [2.0] ;[3.0,2.0,1.0,0.0]
  runs      = ['feedback','tracer','gadget']
  res       = 128
  
  ; index selection for difference ([1.0,0.75,0.5,0.25,0.15,0.05,0.01,first0.15,first1.0])
  rVirFacs = ['1.0 rvir','0.75 rvir','0.5 rvir','0.25 rvir','0.15 rvir','0.05 rvir',$
              '0.01 rvir','first 0.15 rvir','first 1.0 rvir']
  earlierInd = -1 ; -1  ; 0
  laterInd   = -2 ; -2  ; 4
  gcType     = 'all' ; all, gal (includes stars), gmem, inter
  
  ; plot config
  yrange     = [3.0,8.0]
  yrangeNorm = [0.2,5.0]
  xrange     = [0.0,5.0]
  nBins      = 100
  colorInd   = 1
  
  sP = simParams(res=res,run=runs[0],redshift=redshifts[0])
  
  plotName = 'accTimeDeltasVsValMax_'+str(res)+'_'+str(earlierInd)+'_'+str(laterInd)+'_'+gcType+'.eps'
  start_PS,sP.plotPath + plotName, xs=3*3.5, ys=4*3.5
  
    pos = plot_pos(rows=4,cols=2,/gap)
    colors = []
    pinfo = {}
    
    cgText,0.5,0.98,textoidl("rVirFacs: t_1 = " + str(rVirFacs[laterInd]) + $
      ", t_2 = " + str(rVirFacs[earlierInd])) + ' (type='+gcType+')',alignment=0.5,/normal
    
    foreach redshift,redshifts,i do begin
    
    xrange = [0.0,0.7] * redshiftToAgeFlat(redshift)
    
    ; plot (1)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,title="",$
      xtitle=textoidl("( t_{1} - t_{2} ) [Gyr]"),$
      ytitle="Tmax [log K]",/xs,/ys,pos=pos[i*2],/noerase
    p1 = !P & x1 = !X & y1 = !Y
      
    ; plot (2)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeNorm,title="",$
      xtitle=textoidl("( t_{1} - t_{2} ) / t_{2}"),$
      ytitle="Tmax / Tvir",/xs,/ys,pos=pos[i*2+1],/noerase;,/ylog,yminor=0
    p2 = !P & x2 = !X & y2 = !Y
      
	; plot histogram for each run/redshift combination
      foreach run,runs,j do begin
        print, '['+str(j)+'] redshift = '+string(redshift,format='(f4.1)')+' run = '+run
        
        sP_cur = simParams(res=res,run=run,redshift=redshift)
        at     = accretionTimes(sP=sP_cur)
        galcat = galaxyCat(sP=sP_cur)
        
        if sP_cur.trMCPerCell gt 0 then types = ( galcat.type[ replicate_var(galcat.trMC_cc) ] )
        if sP_cur.trMCPerCell eq 0 then types = ( galcat.type )
        if sP_cur.trMCPerCell lt 0 then types = ( galcat.type[ replicate_var(galcat.trVel_cc) ] )
        
        ; add to plot (1)
        delta = selectAccTimeDelta(sP=sP, at=at, galcat=galcat, earlierInd=earlierInd, $
                                 laterInd=laterInd, gcType=gcType)

        !P = p1 & !X = x1 & !Y = y1
        
        ;DEBUG
        maxVals = maxVals(sP=sP_cur)
        
        colors = [colors,sP_cur.colors[colorInd]]
        
        yy = maxVals.maxTemps[delta.w]
	  cgPlot,delta.val,yy,color=colors[-1],psym=3,/overplot
        
        ; add to plot (2)
        !P = p2 & !X = x2 & !Y = y2
        
        yy = 10.0^maxVals.maxTemps[delta.w] / 10.0^at.accHaloTvir[delta.w]
        cgPlot,delta.val,yy,color=colors[-1],psym=3,/overplot
       ; 
       ; xrNormLog = alog10(xrangeNorm)
       ; binSize = (xrNormLog[1] - xrNormLog[0]) / float(nBins)
	 ; h = histogram(alog10(delta.val), bin=binSize, loc=loc, min=xrNormLog[0], max=xrNormLog[1])
	 ; h = float(h)/total(h)
       ; 
	 ; cgPlot,10.0^(loc+binSize*0.5),h,color=colors[j],/overplot
	endforeach ; runs
      
      ; legends
      !P = p1 & !X = x1 & !Y = y1
      legend,['z = '+string(redshift,format='(f3.1)'),runs],textcolors=[0L,colors],/top,/right 
      !P = p2 & !X = x2 & !Y = y2
      legend,['z = '+string(redshift,format='(f3.1)'),runs],textcolors=[0L,colors],/top,/left
      
    endforeach ;redshifts
      
  end_PS
  stop
end

; plotAccTimeDeltasVsHaloMass()

pro plotAccTimeDeltasVsHaloMass

  ; config
  sP = simParams(res=256,run='feedback',redshift=2.0)
  sP2 = simParams(res=256,run='tracer',redshift=2.0)
  
  ; index selection for difference ([1.0,0.75,0.5,0.25,0.15,0.05,0.01,first0.15,first1.0])
  rVirFacs = ['1.0 rvir','0.75 rvir','0.5 rvir','0.25 rvir','0.15 rvir','0.05 rvir',$
              '0.01 rvir','first 0.15 rvir','first 1.0 rvir']
  earlierInd = -1 ; -1 (0) ; 0 ; 2
  laterInd   = -2 ; -2 (4) ; 2 ; 5
  gcType     = 'all' ; all, gal (includes stars), gmem, inter

  ; plot config
  yrange     = [9.5,12.5]
  binSize_yy = 0.1
  logX       = 1
  norm       = 0
  
  ; load
  at     = accretionTimes(sP=sP)
  galcat = galaxyCat(sP=sP)
    
  if sP.trMCPerCell gt 0 then types = ( galcat.type[ replicate_var(galcat.trMC_cc) ] )
  if sP.trMCPerCell eq 0 then types = ( galcat.type )
  if sP.trMCPerCell lt 0 then types = ( galcat.type[ replicate_var(galcat.trVel_cc) ] )
  
  parentMass = galCatParentProperties(sP=sP,galcat=galcat,trRep=(sP.trMCPerCell ne 0),/mass)
  yrange[1] = max(parentMass) - binSize_yy*0.48
    
  ; plot (1) - 2d histogram
  plotStr = sP.savPrefix + str(sP.res) + '_' + str(earlierInd) + '_' + str(laterInd) + '_' + gcType + $
            "_log=" + str(logX) + "_norm=" + str(norm)
  
  start_PS,sP.plotPath + 'accTimeDeltasVsMass_' + plotStr + '.eps'
  
    cgText,0.5,0.96,textoidl("rVirFacs: t_1 = " + str(rVirFacs[laterInd]) + $
      ", t_2 = " + str(rVirFacs[earlierInd])) + ' (type='+gcType+')',alignment=0.5,/normal
                   
    ; add to plot
    delta = selectAccTimeDelta(at=at, galcat=galcat, earlierInd=earlierInd, $
                             laterInd=laterInd, gcType=gcType, norm=norm)

    yy = parentMass[delta.w]
    
    weights = fltarr(n_elements(delta.val)) + 1.0
    
    if logX eq 1 then begin
      xx = alog10( delta.val )
      xrange = [5e-3,3.0]
      xrangeLog = alog10( xrange )
      binSize_xx = 0.05
      xtickv = [0.01,0.05,0.1,0.2,0.5,1.0,2.0,3.0]
      xtickname = ['0.01','0.05','0.1','0.2','0.5','1','2','3']
    endif else begin
      xx = delta.val
      xrange = [0.0,0.5]
      xrangeLog = xrange
      binSize_xx = 0.01
      xtickv = [0.0,0.1,0.2,0.3,0.4,0.5]
      xtickname = ['0','0.1','0.2','0.3','0.4','0.5']
    endelse
    
    ; if not normalizing, extend (t1-t2) axis as a function of the current redshift
    if norm eq 0 then xrange[1] = 0.7 * redshiftToAgeFlat(sP.redshift)
    
    h2 = hist_nd_weight( transpose( [[xx],[yy]] ), weight=weights, [binSize_xx,binSize_yy], $
        min=[xrangeLog[0]-binSize_xx*0.5,yrange[0]-binSize_yy*0.5],$
        max=[xrangeLog[1]+binSize_xx*0.49,yrange[1]+binSize_yy*0.49])
    
    ; colormap (each halo mass row individually scaled)
    logHist = 0
    h2_cmap = colorMapAccTime(h2,logHist=logHist,/byRow)
    
    ; plot
    pos = [0.14,0.14,0.92,0.9]
    
    loadColorTable, 'helix', /reverse ; data
    tvim,h2_cmap,scale=0,pos=pos,/c_map,/noframe,/noerase
          
    if norm eq 1 then xtitle = textoidl("( t_{1} - t_{2} ) / t_{2}")
    if norm eq 0 then xtitle = textoidl("( t_{1} - t_{2} ) [Gyr]")
    ytitle = textoidl("M_{halo}")
          
    loadColorTable,'bw linear' ; axes/frame
    tvim,h2_cmap,/notv,pcharsize=!p.charsize,scale=0,clip=-1,$
      xtitle=xtitle,ytitle=ytitle,barwidth=0.75,lcharsize=!p.charsize-0.2,$
      xrange=xrange,yrange=yrange,xmargin=2.0,pos=pos,/noerase,$
      xlog=logX,xticks=n_elements(xtickv)-1,xtickv=xtickv,xtickname=xtickname
        
  end_PS
  
  if n_elements(sP2) eq 0 then stop
  
  ; make 2d histogram difference
  at2     = accretionTimes(sP=sP2)
  galcat2 = galaxyCat(sP=sP2)
    
  if sP.trMCPerCell gt 0 then types2 = ( galcat2.type[ replicate_var(galcat2.trMC_cc) ] )
  if sP.trMCPerCell eq 0 then types2 = ( galcat2.type )
  if sP.trMCPerCell lt 0 then types2 = ( galcat2.type[ replicate_var(galcat2.trVel_cc) ] )
  
  parentMass2 = galCatParentProperties(sP=sP2,galcat=galcat2,trRep=(sP2.trMCPerCell ne 0),/mass)
  
  ; plot (2) - 2d difference
  start_PS,sP.plotPath + 'accTimeDeltasVsMassDiff_' + plotStr + '.eps'
    cgText,0.5,0.96,textoidl("rVirFacs: t_1 = " + str(rVirFacs[laterInd]) + $
      ", t_2 = " + str(rVirFacs[earlierInd])) + ' (type='+gcType+')',alignment=0.5,/normal
                       
    ; calculate second histogram (same min/max/binsize/etc)
    delta2 = selectAccTimeDelta(at=at2, galcat=galcat2, earlierInd=earlierInd, $
                              laterInd=laterInd, gcType=gcType, norm=norm)
                              
    yy = parentMass2[delta2.w]
    
    weights = fltarr(n_elements(delta2.val)) + 1.0
    
    if logX eq 1 then xx = alog10( delta2.val ) $
    else xx = delta2.val
  
    h2b = hist_nd_weight( transpose( [[xx],[yy]] ), weight=weights, [binSize_xx,binSize_yy], $
         min=[xrangeLog[0]-binSize_xx*0.5,yrange[0]-binSize_yy*0.5],$
         max=[xrangeLog[1]+binSize_xx*0.49,yrange[1]+binSize_yy*0.49])
    
    ; colormap
    fracMin = -0.6 ; 0.25
    fracMax = 0.6 ; 4.0
    
    h2_cmap = colorMapAccTimeDelta(h2,h2b,fracMin=fracMin,fracMax=fracMax)
    
    ; add difference of two histograms to plot
    pos_plot = [0.14,0.14,0.86,0.9]
    pos_cbar = [0.87,0.14,0.91,0.9]
    ndivs = 5
    
    loadColorTable, 'brewer-purplegreen' ; data
    tvim,h2_cmap,scale=0,pos=pos_plot,/c_map,/noframe,/noerase
          
    loadColorTable,'bw linear' ; axes/frame
    tvim,h2_cmap,/notv,pcharsize=!p.charsize,scale=0,clip=-1,$
      xtitle=xtitle,ytitle=ytitle,barwidth=0.75,lcharsize=!p.charsize-0.2,$
      xrange=xrange,yrange=yrange,xmargin=2.0,pos=pos_plot,/noerase,$
      xlog=logX,xticks=n_elements(xtickv)-1,xtickv=xtickv,xtickname=xtickname
     
   loadColorTable, 'brewer-purplegreen' ; data
   cgColorbar,bottom=1,range=[fracMin,fracMax],position=pos_cbar,$
     /vertical,/right,divisions=ndivs,$ ;,ticknames=ticknames,ncolors=255
     title=textoidl("log ( f_{"+sP.plotPrefix+"} / f_{"+sP2.plotPrefix+"} )")

  end_PS
  
  stop
  
end

; shyPlot(): at a target redshift, make a selection of galaxy star tracers which were in gas cells
;            in the previous timestep (newly formed stars). normalize the time elapsed since the 
;            first 1.0rvir crossing of each such tracer by the age of the universe at the time of
;            that crossing. plot the resulting normalized 'accTime' distribution, comparing
;            TRACER to FEEDBACK at z=0,1,2,3

function shyPlotBin, sP=sP, aMode=aMode, allTypes=allTypes

  if keyword_set(allTypes) then allTypesStr = '.allTypes' else allTypesStr = ''

  ; set saveFilename and check for existence
  saveFilename = sP.derivPath + 'shyPlot.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(sP.snap)+'.mode'+str(aMode)+allTypesStr+'.sav'
  
  if file_test(saveFilename) then begin
    restore, saveFilename
    return, r
  endif

  nSnapsBack = 1
  
  ; load at target redshift
  galcat = galaxyCat(sP=sP)
  
  ; make tracer selection (children of stars in galaxies only)
  if keyword_set(allTypes) then begin
    print,'AllTypes = 1'
    
    ; load all star ids
    star_ids = loadSnapshotSubset(sP=sP, partType='star', field='ids')
    
    ; replicate parent id list for all galcat tracers
    tr_parids = galcat.IDs[ replicate_var(galcat.trMC_cc) ]
    
    ; decide overlap, use those tracers (tr_parids not unique, use value_locate approach)
    sort_inds = calcSort(star_ids)
    star_ids_sorted = star_ids[sort_inds]
        
    par_ind = value_locate(star_ids_sorted,tr_parids) ; indices to star_ids_sorted
    par_ind = sort_inds[par_ind>0] ; indices to star_ids (>0 removes -1 entries, which are removed next line)
    tr_parids_ind = where(star_ids[par_ind] eq tr_parids,count_inPar) ; verify we actually matched the ID

    ; save master list of ids/inds to search for
    tr_inds_orig = par_ind[tr_parids_ind]
    tr_ids_orig  = galcat.trMC_ids[ tr_inds_orig ]
  endif else begin
    print,'AllTypes = 0'
    type = ( galcat.type[ replicate_var(galcat.trMC_cc) ] )
  
    tr_inds_orig = where(type eq galcat.types.stars)
    tr_ids_orig = galcat.trMC_ids[ tr_inds_orig ]
  endelse
  
  ; move back requested number of snapshots in a row, load tracer parent ids for this subset
  tr_ids_final  = []
  tr_inds_final = []

  for i=0,nSnapsBack-1 do begin
    sP.snap -= 1
    tr_ids = loadSnapshotSubset(sP=sP, partType='tracerMC', field='tracerids')
  
    idIndexMap = getIDIndexMap(tr_ids, minid=minid)
    tr_ids = !NULL
    tr_inds = idIndexMap[ tr_ids_orig - minid ]
  
    tr_parids = loadSnapshotSubset(sP=sP, partType='tracerMC', field='parentids')
    tr_parids = tr_parids[ tr_inds ]
  
    ; load gas ids and crossmatch
    gas_ids = loadSnapshotSubset(sP=sP, partType='gas', field='ids')
  
    sort_inds = calcSort(gas_ids)
    gas_ids_sorted = gas_ids[sort_inds]
    gas_ind = value_locate(gas_ids_sorted,tr_parids) ; indices to gas_ids_sorted
    gas_ind = sort_inds[gas_ind>0] ; indices to gas_ids (>0 removes -1 entries, which are removed next line)
    w = where(gas_ids[gas_ind] eq tr_parids,count_inPar) ; verify we actually matched the ID

    tr_ids_final  = [ tr_ids_final, tr_ids_orig[w] ]
    tr_inds_final = [ tr_inds_final, tr_inds_orig[w] ] ; index galcat.trMC_ids
    print,count_inPar
  endfor
  
  ; now just insure uniqueness
  if nSnapsBack gt 1 then begin
    tr_ids_final = tr_ids_final[sort(tr_ids_final)]
    tr_inds_final = tr_inds_final[sort(tr_inds_final)]
    
    w = where( tr_ids_final ne shift(tr_ids_final,1), count)
    if count eq 0 then message,'Error, expect some duplicates.'
    tr_ids_final = tr_ids_final[w]
    
    w = where( tr_inds_final ne shift(tr_inds_final,1), count)
    if count eq 0 then message,'Error, expect some duplicates.'
    tr_inds_final = tr_inds_final[w]
  endif
  
  print,' Found: ['+str(n_elements(tr_ids_final))+'] of ['+str(n_elements(tr_ids_orig))+'] tracers.'

  gas_ids         = !NULL
  star_ids        = !NULL
  gas_ind         = !NULL
  gas_ids_sorted  = !NULL
  star_ids_sorted = !NULL
  tr_parids_ind   = !NULL
  par_ind         = !NULL
  tr_parids       = !NULL
  tr_inds         = !NULL
  sort_inds       = !NULL
  w               = !NULL
  
  ; return to original snapshot
  sP.snap += nSnapsBack
  
  ; aMode=1: load accTimes of first 1.0rvir crossing of main progenitor branch
  if aMode eq 1 then begin
    at = accretionTimes(sP=sP)
    at = reform( at.accTime[-1,tr_inds_final] )
  
    ; restrict to accretion mode (all are +_rec) and valid accTimes
    ;am = accretionMode(sP=sP)
    ;am = reform( am[tr_inds_final] )
    ;if accMode eq 'smooth'   then $
    ;  w_at = where(at ne -1 and (am eq 1 or am eq 11),count)
    ;if accMode eq 'clumpy'   then $
    ;  w_at = where(at ne -1 and (am eq 2 or am eq 3 or am eq 12 or am eq 13),count)
    ;if accMode eq 'stripped' then $
    ;  w_at = where(at ne -1 and (am eq 4 or am eq 14),count)
    ;if accMode eq 'all' then $
    ;  w_at = where(at ne -1,count)
    ;am = !NULL
  
    ; filter out -1 values
    at = at[where(at gt 0.0)]
    
    ; convert scalefac to redshift
    at = 1.0/at-1 ; replace at by at[w_at] to use accMode selection
  endif
  
  ; aMode=2: group membership, use mergerTreeSubset() of the subgroups, walk back through 
  ; the snapshots, update Parent (for each tracer) at each step, calculate current Parent (for each
  ; tracer) at each step, or -2 for not in group catalog, save first mismatch time as 
  ; the accretion time entering the subgroup
  if aMode eq 2 then begin
    mt = mergerTreeSubset(sP=sP)
    
    ; at each snapshot holds the tracked main progenitor parents, for each tracer
    curParInds = mt.gcIndOrigTrMC[ tr_inds_final ]
    
    ; the calculated accretion time
    at = fltarr( n_elements(curParInds) ) - 1.0
    
    origSnap   = sP.snap
    targetSnap = sP.groupCatRange[0]
    
    for m=origSnap,targetSnap+1,-1 do begin
      sP.snap = m
          
      ; at each snapshot holds the actual subhalo parent, for each tracer
      actualParInds = lonarr( n_elements(tr_inds_final) )
          
      ; load group catalog
      gcCur = loadGroupCat(sP=sP,/readIDs)
      
      ; calculate actualParInds
        ; load all tracer IDs and re-locate starting IDs
        tr_ids = loadSnapshotSubset(sP=sP, partType='tracerMC', field='tracerids')
        
        idIndexMap = getIDIndexMap(tr_ids, minid=minid)
        tr_ids = !NULL
        tr_inds = idIndexMap[ tr_ids_final - minid ] ; overwrite previous
        
        ; load parent IDs for this subset (don't care about type)
        tr_parids = loadSnapshotSubset(sP=sP, partType='tracerMC', field='parentids', inds=tr_inds)
        
        ; note: tr_parids are NOT UNIQUE, use a value_locate approach (not match)
          sort_inds = calcSort(gcCur.IDs)
          GC_ids_sorted = gcCur.IDs[sort_inds]
        
          par_ind = value_locate(GC_ids_sorted,tr_parids) ; indices to GC_ids_sorted
          par_ind = sort_inds[par_ind>0] ; indices to gc.IDs (>0 removes -1 entries, which are removed next line)
          tr_parids_ind = where(gcCur.IDs[par_ind] eq tr_parids,count_inPar) ; verify we actually matched the ID
          gcCur_ids_ind = par_ind[tr_parids_ind] ; indices of matched parents
          countMatch = n_elements(tr_parids_ind)
          
          ; for those tracers with non-matching parent IDs (not in group catalog), set to -2
          matchMask = bytarr( n_elements(tr_parids) )
          matchMask[tr_parids_ind] = 1B
          
          w_nonMatch = where(matchMask eq 0B,count)
          
          if count gt 0 then begin
            actualParInds[w_nonMatch] = -2
          endif
          
          ; for those tracers with matching parent IDs
          if countMatch gt 0 then begin
            ; do a value locate of the index (of gcCur.IDs) with group.offset
            ; use whole FOF!
            parCandidates = value_locate( gcCur.groupOffset, gcCur_ids_ind )

            ; calculate delta(index,offset) to check if < group.len, to verify actually in group
            parDelta = gcCur_ids_ind - gcCur.groupOffset[ parCandidates ]
            
            w = where(parDelta ge 0 and parDelta lt gcCur.groupLen[ parCandidates ],$
              countCand,comp=wc,ncomp=countOutsideCand)
            
            if countCand gt 0 then begin
              subgroupParInds = gcCur.groupFirstSub[ parCandidates[w] ]
              actualParInds[tr_parids_ind[w]] = subgroupParInds
            endif
            
            ; those outside their candidates (should not happen since considering full fof groups)
            if countOutsideCand gt 0 then message,'Error, reconsider'
          endif
      
      ; for now just mismatching parInds, set accretionTime
      w = where(curParInds ne actualParInds and curParInds ne -1,count)

      if max(at[w]) gt 0.0 then message,'Error: About to set at for tracers with already set at.'
      if count gt 0 and sP.snap eq origSnap then message,'Error: Bad parents at sP.snap'

      at[w] = snapNumToRedshift(sP=sP)
      curParInds[w] = -1 ; no longer search/set on these tracers
           
      frac = float(count)*100/n_elements(tr_ids_final)
      print,'['+str(m)+'] frac now set = '+string(frac,format='(f4.1)')+'% ('+$
        str(count)+' of '+str(n_elements(tr_ids_final))+')'

      ; load mergerTree and move to Parent
      Parent = mergerTree(sP=sP)
        
      ; change to parent IDs for the next snapshot
      w = where(curParInds ne -1,count)
      if count eq 0 then message,'error'
        
      curParInds[w] = Parent[curParInds[w]]
        
      frac = float(count)*100/n_elements(tr_ids_final)
      print,'['+str(m)+'] frac remaining = '+string(frac,format='(f4.1)')+'% ('+$
        str(count)+' of '+str(n_elements(tr_ids_final))+')'
          
    endfor ; snapshot
    
    ; restrict times to those we found (already in redshift)
    at = at[where(at gt 0.0)]
  endif
  
  if aMode ne 1 and aMode ne 2 then message,'Error'
  
  ; age of universe (in Gyr) at each accTime
  r = { age_at  : redshiftToAgeFlat(at) ,$
        age_cur : redshiftToAgeFlat(sP.redshift) }
        
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        
  
  return, r
  
end

; shyPlot():
;  plot (1): delta time between "halo accretion time" and "star forming time" (Gyr)
;  plot (2): delta time as above, normalized by age of universe at "halo accretion time"
; aMode=1: "halo accretion time" = first 1.0 rvir crossing from accretionTimes()
; aMode=2: "halo accretion time" = group membership
; allTypes=1: all tracers in galaxyCat, otherwise just 'star' (e.g. inside 0.15rvir cut) tracers

pro shyPlot

  ; config
  redshifts = [3.0,2.0,1.0,0.0]
  runs      = ['feedback','tracer']
  res       = 512
  aMode     = 1
  allTypes  = 1
  
  ; plot config
  yrangeRaw  = [0.0,0.2] ;[0.0,0.1]
  yrangeNorm = [0.0,0.3]
  nBins     = 20
  colors    = ['red','blue'] ; one per run
  
  start_PS,'shyPlot_aMode'+str(aMode)+'_allTypes'+str(allTypes)+'_'+str(res)+'b.eps', xs=3*3.5, ys=4*3.5
  
    pos = plot_pos(rows=4,cols=2,/gap)
    
    foreach redshift,redshifts,i do begin
    
    ;xrange = [-0.1, redshiftToAgeFlat(redshift)*1.0]
    xrange = [-0.1, 2.0]
    
    ; plot (1)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeNorm,title="",$
      xtitle=textoidl("( t_{SF} - t_{acc} ) / t_{H}(t_{acc})"),$
      ytitle="Fraction",/xs,/ys,pos=pos[i*2],/noerase
	
	; plot histogram for each run/redshift combination
      foreach run,runs,j do begin
        sP_cur = simParams(res=res,run=run,redshift=redshift)
        bV     = shyPlotBin(sP=sP_cur, aMode=aMode, allTypes=allTypes)
        
        val = (bV.age_cur - bV.age_at) / bV.age_at

        if min(val) lt xrange[0] or max(val) gt xrange[1] then begin
          print,'Warning: normalized val min = ',min(val),' max = ',max(val)
        endif
        
        
        binSize = (xrange[1] - xrange[0]) / float(nBins)
	  h = histogram(val, bin=binSize, loc=loc, min=xrange[0], max=xrange[1])
	  h = float(h)/total(h)
	  
	  cgPlot,loc+binSize*0.5,h,color=cgColor(colors[j]),/overplot
	endforeach
	
	; legend
	legend,['z = '+string(redshift,format='(f3.1)'),runs],textcolors=['black',colors],/top,/right
	
    ; plot (2)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeRaw,title="",$
      xtitle=textoidl("( t_{SF} - t_{acc} ) [Gyr]"),$
      ytitle="Fraction",/xs,/ys,pos=pos[i*2+1],/noerase
	
	; plot histogram for each run/redshift combination
      foreach run,runs,j do begin
        sP_cur = simParams(res=res,run=run,redshift=redshift)
        bV     = shyPlotBin(sP=sP_cur, aMode=aMode, allTypes=allTypes)
        
        val = (bV.age_cur - bV.age_at)
        
        if min(val) lt xrange[0] or max(val) gt xrange[1] then begin
          print,'Warning: raw val min = ',min(val),' max = ',max(val)
        endif
        
        binSize = (xrange[1] - xrange[0]) / float(nBins)
	  h = histogram(val, bin=binSize, loc=loc, min=xrange[0], max=xrange[1])
	  h = float(h)/total(h)
	  
	  cgPlot,loc+binSize*0.5,h,color=cgColor(colors[j]),/overplot
	endforeach
	
	; legend
	legend,['z = '+string(redshift,format='(f3.1)'),runs],textcolors=['black',colors],/top,/right
	
    endforeach ;redshifts
      
  end_PS
  
  stop
  
end