; plotAccTimes.pro
; gas accretion project - past radial history of gas elements (plotting)
; dnelson mar.2014

; helper functions:

function selectAccTimeDelta, sP=sP, at=at, am=am, galcat=galcat, $
                             earlierInd=earlierInd, laterInd=laterInd, $
                             gcType=gcType, accMode=accMode, norm=norm, $
                             parentMass=parentMass, pMassBin=pMassBin

  if keyword_set(pMassBin) and ~keyword_set(parentMass) then message,'Error'
  units = getUnits()
  
  ; mask for accretion mode
  modeMask = accModeInds(at=at, sP=sP, accMode=accMode, /maskAndInds)
  modeMask = modeMask.mask
  
  if gcType ne 'all' then begin
    if sP.trMCPerCell gt 0 then types = ( galcat.type[ replicate_var(galcat.trMC_cc) ] )
    if sP.trMCPerCell eq 0 then types = ( galcat.type )
    if sP.trMCPerCell lt 0 then types = ( galcat.type[ replicate_var(galcat.trVel_cc) ] )
  endif
  
  ; mask for parent halo mass
  massMask = bytarr( n_elements(modeMask) )
  
  if keyword_set(pMassBin) then begin
    w = where(parentMass ge pMassBin[0] and parentMass lt pMassBin[1],count)
    if count eq 0 then message,'Error'
    massMask[w] = 1B
  endif else begin
    massMask += 1B
  endelse
  
  ; all times are interpolated, so if exactly the same exclude with strict LT
  if gcType eq 'all' then $
    w = where(at.accTime[earlierInd,*] ge 0 and at.accTime[laterInd,*] ge 0 and $
              at.accTime[earlierInd,*] lt at.accTime[laterInd,*] and $
              modeMask eq 1B and massMask eq 1B,count)
  if gcType eq 'gal' then $
    w = where(at.accTime[earlierInd,*] ge 0 and at.accTime[laterInd,*] ge 0 and $
              at.accTime[earlierInd,*] lt at.accTime[laterInd,*] and $
              (types eq galcat.types.gal or types eq galcat.types.stars) and $
              modeMask eq 1B and massMask eq 1B,count)
  if gcType eq 'gmem' then $
    w = where(at.accTime[earlierInd,*] ge 0 and at.accTime[laterInd,*] ge 0 and $
              at.accTime[earlierInd,*] lt at.accTime[laterInd,*] and $
              types eq galcat.types.gmem and $
              modeMask eq 1B and massMask eq 1B,count)
  if gcType eq 'inter' then $
    w = where(at.accTime[earlierInd,*] ge 0 and at.accTime[laterInd,*] ge 0 and $
              at.accTime[earlierInd,*] lt at.accTime[laterInd,*] and $
              types eq galcat.types.inter and $
              modeMask eq 1B and massMask eq 1B,count)
         
  age_earlier = reform( redshiftToAgeFlat( 1.0/at.accTime[earlierInd,w]-1.0 ) )
  age_later   = reform( redshiftToAgeFlat( 1.0/at.accTime[laterInd,w]-1.0 ) )
  
  val = age_later - age_earlier
  
  if count eq 0 then message,'Error'
  
  if keyword_set(norm) then val /= age_earlier
  
  return, {w:w,val:val,norm_fac:age_earlier}
end

function colorMapAccTime, h2, logHist=logHist, byRow=byRow, byCol=byCol
  ; config
  min = 10.0
  max = 255.0
  
  h2_cmap = h2
  
  ; take values to log and lift all zeros to minimum nonzero
  if keyword_set(logHist) then begin
    w = where(h2_cmap gt 0,count)
    ;h2_cmap /= float(total(h2_cmap,/int))
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
    
    return, h2_cmap
  endif
  
  ; normalization by column?
  if keyword_set(byCol) then begin
    nRows = n_elements(h2_cmap[*,0])
    
    for i=0,nRows-1 do begin
      curCol = reform( h2_cmap[i,*] )
      curCol = fix( (curCol-min(curCol))/(max(curCol)-0.0) * (max-min) + min )
      h2_cmap[i,*] = curCol
    endfor
    
    return, h2_cmap
  endif
  
  ; else: linear stretch from 10 (non-white) to 255
  h2_cmap = fix( (h2_cmap-min(h2_cmap))/(max(h2_cmap)-0.0) * (max-min) + min )
  
  return, h2_cmap
end

function colorMapAccTimeDelta, h2a, h2b, fracMin=fracMin, fracMax=fracMax, totalNorm=totalNorm
  ; config
  min = 0.0
  max = 255.0
  
  if n_elements(fracMin) eq 0 or n_elements(fracMax) eq 0 then message,'Error'
  
  h2d_cmap = h2a*0.0 + 255/2 ; center missing values at white
 
  ; normalize over entire plane (not by row or column?)
  if n_elements(totalNorm) gt 0 then begin
    h2a /= float(total(h2a,/int))
    h2b /= float(total(h2b,/int))
  endif
 
  ; row by row
  nRows = n_elements(h2d_cmap[0,*])

  for i=0,nRows-1 do begin
    ; extract
    aRow = reform( h2a[*,i] )
    bRow = reform( h2b[*,i] )
    
    ; normalize each to their total (now a fraction)
    if n_elements(totalNorm) eq 0 then begin
      aTotal = total(aRow,/int)
      bTotal = total(bRow,/int)

      if aTotal eq 0 or bTotal eq 0 then begin
        ; set missing data to NaN (or leave at zero difference)
        ;h2d_cmap[*,i] = !values.f_nan
        continue
      endif
    
      aRow /= aTotal
      bRow /= bTotal
    endif
    
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

; plotAccTimeScatter():

pro plotAccTimeScatter

  ; config
  sP = simParams(res=128,run='tracer',redshift=2.0)
  
  ; index selection for difference ([1.0,0.75,0.5,0.25,0.15,0.05,0.01,first0.15,first1.0])
  rVirFacs = ['1.0 rvir','0.75 rvir','0.5 rvir','0.25 rvir','0.15 rvir','0.05 rvir',$
              '0.01 rvir','first 0.15 rvir','first 1.0 rvir']
  earlierInd = -1 ; -1  ; 0
  laterInd   = -2 ; -2  ; 4
  gcType     = 'all' ; all, gal (includes stars), gmem, inter
  accMode    = 'smooth' ; smooth_rec, clumpy_rec, stripped_rec, all
  
  ; plot config
  xyrange      = [0.1,1/(1+sP.redshift)+0.02]
  xyrangeDelta = [0.01, 0.7 * redshiftToAgeFlat(sP.redshift)]
  
  nBins      = 100
  colorInd   = 1
  logHist    = 1 ; log 2d counts in each bin
  
  inds = [0,1,2,3,4,5,7,8]
  indPairs = list( [0,3], [2,4], [-1,-2], [0,5] )
  
  ; load
  at     = accretionTimes(sP=sP)
  am     = accretionMode(sP=sP)
  galcat = galaxyCat(sP=sP)
        
  ; just used for delta.w:
  delta = selectAccTimeDelta(sP=sP, at=at, am=am, galcat=galcat, earlierInd=earlierInd, $
                             laterInd=laterInd, gcType=gcType, accMode=accMode)

  for i=0,n_elements(inds)-1 do begin
    ind0 = inds[i]
    
  ; (A)
  plotStr = sP.plotPrefix+str(sP.res)+'_gc-'+gcType+'_mode-'+accMode+'_log-'+str(logHist)+'_t'+str(ind0)
  
  start_PS,sP.plotPath + 'accTimeScatter_' + plotStr + '_a.eps', xs=10, ys=10
  
    pos = plot_pos(rows=2,cols=2,/gap)
  
    ; plot
    for j=0,3 do begin
      ind1 = inds[j]
      xx = at.accTime[ind0,delta.w]
      yy = at.accTime[ind1,delta.w]
    
      cgPlot,[0],[0],/nodata,xrange=xyrange,yrange=xyrange,title="",$
        xtitle=textoidl("t_{"+str(ind0)+"} ("+rVirFacs[ind0]+")"),$
        ytitle=textoidl("t_{"+str(ind1)+"} ("+rVirFacs[ind1]+")"),$
        /xs,/ys,pos=pos[j],/noerase
      
      cgPlot,[0.12,0.34],[0.12,0.34],line=0,color=cgColor('orange'),/overplot
      cgPlot,xx,yy,psym=3,/overplot ; scatterplot
    endfor
      
  end_PS
  
  ; (B)
  start_PS,sP.plotPath + 'accTimeScatter_' + plotStr + '_b.eps', xs=10, ys=10
  
    pos = plot_pos(rows=2,cols=2,/gap)
  
    ; plot
    for j=4,7 do begin
      ind1 = inds[j]
      xx = at.accTime[ind0,delta.w]
      yy = at.accTime[ind1,delta.w]
    
      cgPlot,[0],[0],/nodata,xrange=xyrange,yrange=xyrange,title="",$
        xtitle=textoidl("t_{"+str(ind0)+"} ("+rVirFacs[ind0]+")"),$
        ytitle=textoidl("t_{"+str(ind1)+"} ("+rVirFacs[ind1]+")"),$
        /xs,/ys,pos=pos[j-4],/noerase
      
      cgPlot,[0.12,0.34],[0.12,0.34],line=0,color=cgColor('orange'),/overplot
      cgPlot,xx,yy,psym=3,/overplot ; scatterplot
    endfor
      
  end_PS
  
  endfor ;ind0
  
  ; "colors": deltas vs eachother
  foreach indPair,indPairs do begin
    
    plotStr = sP.plotPrefix+str(sP.res)+'_gc-'+gcType+'_mode-'+accMode+'_log-'+str(logHist)+$
              '_ind='+str(indPair[0])+'_'+str(indPair[1]);+'_vs='+str(indPairB[0])+'_'+str(indPairB[1])
            
    start_PS,sP.plotPath + 'accTimeScatter_' + plotStr + '.eps', xs=10, ys=10
  
      pos = plot_pos(rows=2,cols=2,/gap)
  
      foreach indPairB,indPairs,j do begin
  
        ; plot
        delta  = selectAccTimeDelta(sP=sP, at=at, galcat=galcat, earlierInd=indPair[0], $
                                   laterInd=indPair[1], gcType=gcType, accMode=accMode)
        deltaB = selectAccTimeDelta(sP=sP, at=at, galcat=galcat, earlierInd=indPairB[0], $
                                   laterInd=indPairB[1], gcType=gcType, accMode=accMode)
                                   
        ;xx = reform( at.accTime[indPair[0],delta.w] ) - reform( at.accTime[indPair[1],delta.w] )
        ;yy = reform( at.accTime[indPairB[0],delta.w] ) - reform( at.accTime[indPairB[1],delta.w] )
        
        cgPlot,[0],[0],/nodata,xrange=xyrangeDelta,yrange=xyrangeDelta,title="",$
          xtitle=textoidl("t_{"+str(indPair[0])+"} - t_{"+str(indPair[1])+"} ("+$
                          rVirFacs[indPair[0]]+" - "+rVirFacs[indPair[1]]+")"),$
          ytitle=textoidl("t_{"+str(indPairB[0])+"} - t_{"+str(indPairB[1])+"} ("+$
                          rVirFacs[indPairB[0]]+" - "+rVirFacs[indPairB[1]]+")"),$
          /xs,/ys,/xlog,xminor=0,/ylog,yminor=0,pos=pos[j],/noerase
      
        cgPlot,[2*xyrangeDelta[0],0.9*xyrangeDelta[1]],[2*xyrangeDelta[0],0.9*xyrangeDelta[1]],$
          line=0,color=cgColor('orange'),/overplot
        cgPlot,delta.val,deltaB.val,psym=3,/overplot ; scatterplot
      
      endforeach ;indPairsB
      
    end_PS
  endforeach ;indPairs
  
end

; plotAccTimeDeltas(): one dimensional histograms, comparing multiple runs at z=0,1,2,3

pro plotAccTimeDeltas

  ; config
  redshifts = [3.0,2.0,1.0,0.0]
  runs      = ['feedback','tracer','gadget']
  res       = 512
  
  ; index selection for difference ([1.0,0.75,0.5,0.25,0.15,0.05,0.01,first0.15,first1.0])
  rVirFacs = ['1.0 rvir','0.75 rvir','0.5 rvir','0.25 rvir','0.15 rvir','0.05 rvir',$
              '0.01 rvir','first 0.15 rvir','first 1.0 rvir']
  earlierInd = -1 ; -1  ; 0
  laterInd   = -2 ; -2  ; 4
  gcType     = 'all' ; all, gal (includes stars), gmem, inter
  accModes   = ['smooth_rec','clumpy_rec','all'] ; smooth_rec, clumpy_rec, stripped_rec, all
  
  ; plot config
  yrange     = [1e-4,0.2]
  xrange     = [0.0,5.0]
  xrangeNorm = [7e-3,1.5]
  nBins      = 100
  colorInd   = 1
  
  foreach accMode,accModes do begin
  
  sP = simParams(res=res,run=runs[0],redshift=redshifts[0])
  
  plotName = 'accTimeDeltas_'+str(res)+'_'+str(earlierInd)+'_'+str(laterInd)+'_gc-'+gcType+$
             '_mode-' + accMode + '.eps'
  start_PS,sP.plotPath + plotName, xs=3*3.5, ys=4*3.5
  
    pos = plot_pos(rows=4,cols=2,/gap)
    colors = []
    pinfo = {}
    
    cgText,0.5,0.98,textoidl("rVirFacs: t_1 = " + str(rVirFacs[laterInd]) + $
      ", t_2 = " + str(rVirFacs[earlierInd])) + ' (type='+gcType+') (mode='+accMode+')',alignment=0.5,/normal
    
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
        
        sP = simParams(res=res,run=run,redshift=redshift)
        at     = accretionTimes(sP=sP)
        am     = accretionMode(sP=sP)
        galcat = galaxyCat(sP=sP)
        
        ; add to plot (1)
        delta = selectAccTimeDelta(sP=sP, at=at, galcat=galcat, earlierInd=earlierInd, $
                                 laterInd=laterInd, gcType=gcType, accMode=accMode)

        !P = p1 & !X = x1 & !Y = y1
        
        ;DEBUG (match to fluidTracks)
        if 0 then begin
          ww = where(delta.val lt -0.4 and delta.val gt -0.5)
          tInd = w[ww[1]]
        
          tracks = tracksFluid(sP=sP)
          mt     = mergerTreeSubset(sP=sP)
        
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
        
        colors = [colors,sP.colors[colorInd]]
	  cgPlot,loc+binSize*0.5,h,color=colors[-1],/overplot
        
        ; add to plot (2)
        delta = selectAccTimeDelta(sP=sP, at=at, am=am, galcat=galcat, earlierInd=earlierInd, $
                                   laterInd=laterInd, gcType=gcType, accMode=accMode, /norm)
                                 
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
  endforeach ;accModes
  
  stop
end

; plotAccTimeDeltaVsValMax(): correlations between time_delta and max_val quantities

pro plotAccTimeDeltaVsValMax

  ; config
  redshift  = 2.0
  runs      = ['feedback','tracer','gadget']
  res       = 128
  
  ; index selection for difference ([1.0,0.75,0.5,0.25,0.15,0.05,0.01,first0.15,first1.0])
  rVirFacs = ['1.0 rvir','0.75 rvir','0.5 rvir','0.25 rvir','0.15 rvir','0.05 rvir',$
              '0.01 rvir','first 0.15 rvir','first 1.0 rvir']
  earlierInd = -1 ; -1  ; 0
  laterInd   = -2 ; -2  ; 4
  gcType     = 'all' ; all, gal (includes stars), gmem, inter
  accMode    = 'smooth' ; smooth_rec, clumpy_rec, stripped_rec, all
  plotVal    = 'temp'
  
  ; plot config
  xrange     = [0.0,0.7] * redshiftToAgeFlat(redshift)
  nBins      = 100
  colorInd   = 1
  logHist    = 1 ; log 2d counts in each bin
  byRow      = 0 ; normalize 2d counts by row
  byCol      = 0 ; normalize 2d counts by column
  massBins   = list( [10.4,10.6], [11.0,11.2], [11.6, 11.8] )
  
  if plotVal eq 'temp' then begin
    yrange     = [3.8,8.2]
    yrangeNorm = [0.01,30.0]
    ytitle     = "Tmax [log K]"
    ytitleNorm = "log ( Tmax / Tvir )"
  endif
  if plotVal eq 'ent' then begin
    yrange     = [5.2,9.6]
    yrangeNorm = [0.05,300.0]
    ytitle     = "S [cgs]"
    ytitleNorm = "log ( S / S200 )"
  endif
  
  sP = simParams(res=res,run=runs[0],redshift=redshift)
  
  plotStr = str(res)+'_'+str(earlierInd)+'_'+str(laterInd)+$
            '_gc-'+gcType+'_mode-'+accMode+'_'+str(logHist)+str(byRow)+str(byCol)+'.eps'
  
  if 0 then begin
  
  start_PS,sP.plotPath + 'accTimeDeltasVs_'+plotVal+'Max_' + plotStr, xs=4*n_elements(runs), ys=8
  
  cgText,0.5,0.98,textoidl("rVirFacs: t_1 = " + str(rVirFacs[laterInd]) + $
    ", t_2 = " + str(rVirFacs[earlierInd])) + ' (type='+gcType+') (mode='+accMode+')',alignment=0.5,/normal
  
  pos = plot_pos(rows=2,cols=n_elements(runs),/gap)
  
  foreach run,runs,j do begin
  
    print, '['+str(j)+'] redshift = '+string(redshift,format='(f4.1)')+' run = '+run
  
    ; plot each run/normalized or not combination
    sP     = simParams(res=res,run=run,redshift=redshift)
    at     = accretionTimes(sP=sP)
    am     = accretionMode(sP=sP)
    galcat = galaxyCat(sP=sP)
        
    ; get accTimeDelta and maximum values
    delta = selectAccTimeDelta(sP=sP, at=at, am=am, galcat=galcat, earlierInd=earlierInd, $
                               laterInd=laterInd, gcType=gcType, accMode=accMode)

    maxVals = maxVals(sP=sP)
    
    parentMass = galCatParentProperties(sP=sP,galcat=galcat,trRep=(sP.trMCPerCell ne 0),/mass)
    parentS200 = codeMassToVirEnt( parentMass, sP=sP, /log ) ; only used for 'ent'

    weights = fltarr(n_elements(delta.val)) + 1.0 ; uniform
    
    ; plot (1): delta_t vs tmax [K]
    if plotVal eq 'temp' then yy = maxVals.maxTemps[delta.w]
    if plotVal eq 'ent'  then yy = maxVals.maxEnt[delta.w]
    
    print,minmax(yy)
    
    binSize_xx = 0.1 / (sP.res/128>1)
    binSize_yy = 0.2 / (sP.res/128>1)
        
    h2 = hist_nd_weight( transpose( [[delta.val],[yy]] ), weight=weights, [binSize_xx,binSize_yy], $
        min=[xrange[0]-binSize_xx*0.50,yrange[0]-binSize_yy*0.50],$
        max=[xrange[1]+binSize_xx*0.49,yrange[1]+binSize_yy*0.49])
    
    ; colormap (each halo mass row individually scaled)
    h2_cmap = colorMapAccTime(h2,logHist=logHist,byRow=byRow,byCol=byCol)
    
    ; plot
    loadColorTable, 'helix', /reverse ; data
    tvim,h2_cmap,scale=0,pos=pos[j],/c_map,/noframe,/noerase
          
    loadColorTable,'bw linear' ; axes/frame
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,title="",$
      xtitle=textoidl("( t_{1} - t_{2} ) [Gyr]"),$
      ytitle=ytitle,/xs,/ys,pos=pos[j],/noerase
    ;cgPlot,delta.val,yy,psym=3,/overplot ; scatterplot
    legend,[sP.simName],textcolors=[sP.colors[colorInd]],/top,/right
        
    ; plot (2): delta_t vs tmax/tvir
    if plotVal eq 'temp' then $
      yy = alog10( 10.0^maxVals.maxTemps[delta.w] / 10.0^at.accHaloTvir[delta.w] )
    if plotVal eq 'ent' then $
      yy = alog10( 10.0^maxVals.maxEnt[delta.w] / 10.0^parentS200[delta.w] )
      
    binSize_xx = 0.1 / (sP.res/128>1)
    binSize_yy = 0.15 / (sP.res/128>1)
        
    h2 = hist_nd_weight( transpose( [[delta.val],[yy]] ), weight=weights, [binSize_xx,binSize_yy], $
        min=[xrange[0]-binSize_xx*0.50,alog10(yrangeNorm[0])-binSize_yy*0.50],$
        max=[xrange[1]+binSize_xx*0.49,alog10(yrangeNorm[1])+binSize_yy*0.49])
    
    ; colormap (each halo mass row individually scaled)
    h2_cmap = colorMapAccTime(h2,logHist=logHist,byRow=byRow,byCol=byCol)
      
    ; plot
    loadColorTable, 'helix', /reverse ; data
    tvim,h2_cmap,scale=0,pos=pos[j+3],/c_map,/noframe,/noerase
          
    loadColorTable,'bw linear' ; axes/frame
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=alog10(yrangeNorm),title="",$
      xtitle=textoidl("( t_{1} - t_{2} ) [Gyr]"),$
      ytitle=ytitleNorm,/xs,/ys,pos=pos[j+3],/noerase;,/ylog,yminor=0
    ;cgPlot,delta.val,yy,psym=3,/overplot ;scatterplot
    legend,[sP.simName],textcolors=[sP.colors[colorInd]],/top,/right  
      
  endforeach ; runs
      
  end_PS
  
  endif ;0
  if 0 then begin
  
  ; plot (2) - different of first two runs
  start_PS,sP.plotPath + 'accTimeDiffVs_'+plotVal+'Max_' + plotStr, xs=12, ys=5
  
    sPs = {}
    hh1 = {} ; tmax
    hh2 = {} ; tmax/tvir
  
    for j=0,1 do begin
      print,runs[j]
      ; load
      sP      = simParams(res=res,run=runs[j],redshift=redshift)
      at      = accretionTimes(sP=sP)
      am      = accretionMode(sP=sP)
      galcat  = galaxyCat(sP=sP)
      maxVals = maxVals(sP=sP)
      delta   = selectAccTimeDelta(sP=sP, at=at, am=am, galcat=galcat, earlierInd=earlierInd, $
                                   laterInd=laterInd, gcType=gcType, accMode=accMode)

      ; histogram
      weights = fltarr(n_elements(delta.val)) + 1.0 ; uniform
      yy = maxVals.maxTemps[delta.w]
    
      binSize_xx = 0.08 / (sP.res/128>1)
      binSize_yy = 0.2 / (sP.res/128>1)
        
      h2 = hist_nd_weight( transpose( [[delta.val],[yy]] ), weight=weights, [binSize_xx,binSize_yy], $
          min=[xrange[0]-binSize_xx*0.50,yrange[0]-binSize_yy*0.50],$
          max=[xrange[1]+binSize_xx*0.49,yrange[1]+binSize_yy*0.49])
          
      hh1 = mod_struct( hh1, 'h_'+str(j), h2 )
      
      ; histogram: tmax/tvirr
      weights = fltarr(n_elements(delta.val)) + 1.0 ; uniform
      yy = alog10( 10.0^maxVals.maxTemps[delta.w] / 10.0^at.accHaloTvir[delta.w] )
	
      binSize_xx = 0.08 / (sP.res/128>1)
      binSize_yy = 0.12 / (sP.res/128>1)
        
      h2 = hist_nd_weight( transpose( [[delta.val],[yy]] ), weight=weights, [binSize_xx,binSize_yy], $
          min=[xrange[0]-binSize_xx*0.50,alog10(yrangeNorm[0])-binSize_yy*0.50],$
          max=[xrange[1]+binSize_xx*0.49,alog10(yrangeNorm[1])+binSize_yy*0.49])
          
      hh2 = mod_struct( hh2, 'h_'+str(j), h2 )
      
      sPs = mod_struct( sPs, 'sP'+str(j), sP )
    endfor
    
    ; different histograms and colormap
    fracMin = -0.6 ; 0.25
    fracMax = 0.6 ; 4.0
    totalNorm = 1
    
    h1_cmap = colorMapAccTimeDelta(hh1.(0),hh1.(1),fracMin=fracMin,fracMax=fracMax,totalNorm=totalNorm)
    h2_cmap = colorMapAccTimeDelta(hh2.(0),hh2.(1),fracMin=fracMin,fracMax=fracMax,totalNorm=totalNorm)

    ; add difference of two histograms to plot
    pos_plot1 = [0.06,0.14,0.44,0.9]
    pos_plot2 = [0.52,0.14,0.90,0.9]
    pos_cbar = [0.91,0.14,0.94,0.9]
    ndivs = 5
    
    ;pos = plot_pos(rows=1,cols=2,/gap)
    
    loadColorTable, 'brewer-redpurple' ; data
    tvim,h1_cmap,scale=0,pos=pos_plot1,/c_map,/noframe,/noerase
          
    loadColorTable,'bw linear' ; axes/frame
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,title="",$
      xtitle=textoidl("( t_{1} - t_{2} ) [Gyr]"),$
      ytitle="Tmax [log K]",/xs,/ys,pos=pos_plot1,/noerase
      
    loadColorTable, 'brewer-redpurple' ; data
    tvim,h2_cmap,scale=0,pos=pos_plot2,/c_map,/noframe,/noerase
          
    loadColorTable,'bw linear' ; axes/frame
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=alog10(yrangeNorm),title="",$
      xtitle=textoidl("( t_{1} - t_{2} ) [Gyr]"),$
      ytitle="log ( Tmax / Tvir )",/xs,/ys,pos=pos_plot2,/noerase
     
   ; colorbar
   loadColorTable, 'brewer-redpurple' ; data
   cgColorbar,bottom=1,range=[fracMin,fracMax],position=pos_cbar,$
     /vertical,/right,divisions=ndivs,$ ;,ticknames=ticknames,ncolors=255
     title=textoidl("log ( f_{"+sPs.(0).plotPrefix+"} / f_{"+sPs.(1).plotPrefix+"} )")
    
  end_PS
  endif ;0
  
  ; plot (3) - first run split into 6 halo mass bins
  ;if 0 then begin
  
  foreach run,runs,j do begin
  
    start_PS,sP.plotPath + 'accTimeDeltasVs_'+sP.plotPrefix + str(sP.res) + "_" +$
      plotVal+'Max_' + plotStr, xs=11, ys=8
  
    cgText,0.5,0.98,textoidl(sP.run + " " + str(sP.res) + " rVirFacs: t_1 = " + str(rVirFacs[laterInd]) + $
      ", t_2 = " + str(rVirFacs[earlierInd])) + ' (type='+gcType+') (mode='+accMode+')',alignment=0.5,/normal
  
    pos = plot_pos(rows=2,cols=3,/gap)
  
    print, '['+run+'] redshift = '+string(redshift,format='(f4.1)')+' run = '+run
  
    ; plot each run/normalized or not combination
    sP     = simParams(res=res,run=run,redshift=redshift)
    at     = accretionTimes(sP=sP)
    am     = accretionMode(sP=sP)
    galcat = galaxyCat(sP=sP)
        
    maxVals = maxVals(sP=sP)
    
    parentMass = galCatParentProperties(sP=sP,galcat=galcat,trRep=(sP.trMCPerCell ne 0),/mass)
    parentS200 = codeMassToVirEnt( parentMass, sP=sP, /log ) ; only used for 'ent'
    
    foreach massBin,massBins,k do begin
    
      ; get accTimeDelta and maximum values
      delta = selectAccTimeDelta(sP=sP, at=at, am=am, galcat=galcat, earlierInd=earlierInd, $
                                 laterInd=laterInd, gcType=gcType, accMode=accMode, $
                                 parentMass=parentMass,pMassBin=massBin)
                                 
      massBinStr = string(massBin[0],format='(f4.1)') + " - " + string(massBin[1],format='(f4.1)')

      weights = fltarr(n_elements(delta.val)) + 1.0 ; uniform
    
      ; plot (1): delta_t vs tmax [K]
      if plotVal eq 'temp' then yy = maxVals.maxTemps[delta.w]
      if plotVal eq 'ent'  then yy = maxVals.maxEnt[delta.w]
    
      print,minmax(yy)
    
      binSize_xx = 0.1 / (sP.res/128>1)
      binSize_yy = 0.2 / (sP.res/128>1)
        
      h2 = hist_nd_weight( transpose( [[delta.val],[yy]] ), weight=weights, [binSize_xx,binSize_yy], $
          min=[xrange[0]-binSize_xx*0.50,yrange[0]-binSize_yy*0.50],$
          max=[xrange[1]+binSize_xx*0.49,yrange[1]+binSize_yy*0.49])
    
      ; colormap (each halo mass row individually scaled)
      h2_cmap = colorMapAccTime(h2,logHist=logHist,byRow=byRow,byCol=byCol)
    
      ; plot
      loadColorTable, 'helix', /reverse ; data
      tvim,h2_cmap,scale=0,pos=pos[k],/c_map,/noframe,/noerase
          
      loadColorTable,'bw linear' ; axes/frame
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,title="",$
        xtitle=textoidl("( t_{1} - t_{2} ) [Gyr]"),$
        ytitle=ytitle,/xs,/ys,pos=pos[k],/noerase
      legend,massBinStr,textcolors=[sP.colors[colorInd]],/top,/right
        
      ; plot (2): delta_t vs tmax/tvir
      if plotVal eq 'temp' then $
        yy = alog10( 10.0^maxVals.maxTemps[delta.w] / 10.0^at.accHaloTvir[delta.w] )
      if plotVal eq 'ent' then $
        yy = alog10( 10.0^maxVals.maxEnt[delta.w] / 10.0^parentS200[delta.w] )
      
      binSize_xx = 0.1 / (sP.res/128>1)
      binSize_yy = 0.15 / (sP.res/128>1)
      
      ; DEBUG:
      z = linspace( min(delta.val), max(delta.val), 1000 )
      zz = replicate( alog10(yrangeNorm[1]) + binSize_yy*0.48, 1000 )
      val = [delta.val,z]
      yy = [yy,zz]
      ; END DEBUG
        
      h2 = hist_nd_weight( transpose( [[val],[yy]] ), weight=weights, [binSize_xx,binSize_yy], $
          min=[xrange[0]-binSize_xx*0.50,alog10(yrangeNorm[0])-binSize_yy*0.50],$
          max=[xrange[1]+binSize_xx*0.49,alog10(yrangeNorm[1])+binSize_yy*0.49], rev=ri)
      stop
      ; colormap (each halo mass row individually scaled)
      h2_cmap = colorMapAccTime(h2,logHist=logHist,byRow=byRow,byCol=byCol)
      
      ; plot
      loadColorTable, 'helix', /reverse ; data
      tvim,h2_cmap,scale=0,pos=pos[k+3],/c_map,/noframe,/noerase
          
      loadColorTable,'bw linear' ; axes/frame
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=alog10(yrangeNorm),title="",$
        xtitle=textoidl("( t_{1} - t_{2} ) [Gyr]"),$
        ytitle=ytitleNorm,/xs,/ys,pos=pos[k+3],/noerase;,/ylog,yminor=0
      legend,massBinStr,textcolors=[sP.colors[colorInd]],/top,/right
    
    endforeach ; massBins
      
  endforeach ; runs
      
  end_PS
  
  ;endif ;0
  
  stop
end

; plotAccTimeDeltasVsHaloMass(): 2d histograms, individual runs and difference histograms between 2 runs

pro plotAccTimeDeltasVsHaloMass

  ; config
  sP = simParams(res=512,run='feedback',redshift=2.0)
  ;sP2 = simParams(res=128,run='tracer',redshift=3.0)
  
  ; index selection for difference ([1.0,0.75,0.5,0.25,0.15,0.05,0.01,first0.15,first1.0])
  rVirFacs = ['1.0 rvir','0.75 rvir','0.5 rvir','0.25 rvir','0.15 rvir','0.05 rvir',$
              '0.01 rvir','first 0.15 rvir','first 1.0 rvir']
  earlierInd = -1 ; -1 (0) ; 0 ; 2
  laterInd   = -2 ; -2 (4) ; 2 ; 5
  gcType     = 'all' ; all, gal (includes stars), gmem, inter
  accMode    = 'smooth' ; smooth_rec, clumpy_rec, stripped_rec, all
  
  ; plot config
  yrange     = [9.5,12.5]
  binSize_yy = 0.1
  logX       = 1
  norm       = 0
  massBins   = list( [9.45,9.55], [9.9,10.1], [10.4,10.6], [10.8,11.2], [11.3,11.7], [11.7,12.3] )
  sK         = 3
  
  ; load
  at     = accretionTimes(sP=sP)
  am     = accretionMode(sP=sP)
  galcat = galaxyCat(sP=sP)
    
  if sP.trMCPerCell gt 0 then types = ( galcat.type[ replicate_var(galcat.trMC_cc) ] )
  if sP.trMCPerCell eq 0 then types = ( galcat.type )
  if sP.trMCPerCell lt 0 then types = ( galcat.type[ replicate_var(galcat.trVel_cc) ] )
  
  parentMass = galCatParentProperties(sP=sP,galcat=galcat,trRep=(sP.trMCPerCell ne 0),/mass)
  yrange[1] = max(parentMass) - binSize_yy*0.48
    
  ; plot (1) - 2d histogram
  plotStr = sP.savPrefix + str(sP.res) + '_' + str(earlierInd) + '_' + str(laterInd) + '_gc-' + gcType + $
            "_mode-" + accMode + "_log=" + str(logX) + "_norm=" + str(norm)
  
  start_PS,sP.plotPath + 'accTimeDeltasVsMass_' + plotStr + '.eps'
  
    cgText,0.5,0.96,textoidl("rVirFacs: t_1 = " + str(rVirFacs[laterInd]) + $
      ", t_2 = " + str(rVirFacs[earlierInd])) + ' (type='+gcType+') (mode='+accMode+')',alignment=0.5,/normal
                   
    ; add to plot
    delta = selectAccTimeDelta(sP=sP, at=at, am=am, galcat=galcat, earlierInd=earlierInd, $
                               laterInd=laterInd, gcType=gcType, accMode=accMode, norm=norm)

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
  
  ; plot (2) - 1d profiles
  set_plot,'ps'
  mbColors = reverse( sampleColorTable('blue-red2', n_elements(massBins), bounds=[0.1,0.9]) )
  
  plotStr = sP.savPrefix + str(sP.res) + '_' + str(earlierInd) + '_' + str(laterInd) + '_gc-' + gcType + $
            "_mode-" + accMode + "_log=" + str(logX) + "_norm=" + str(norm)
  
  start_PS,sP.plotPath + 'accTimeDeltasMassBins_' + plotStr + '.eps', xs=7.5, ys=10.5
    pos = plot_pos(row=2,col=1,/gap)
    
    ; plot A (no normalization)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=[0.001,2.0],pos=pos[0],$
      xtitle=textoidl("( t_{1} - t_{2} ) [Gyr]"),ytitle="Fraction",$
      /xs,/ys,/ylog,yminor=0,xlog=logX,xminor=(logX ne 1),$
      title=textoidl("rVirFacs: t_1 = " + str(rVirFacs[laterInd]) + $
      ", t_2 = " + str(rVirFacs[earlierInd])) + ' (type='+gcType+') (mode='+accMode+')'
    
    ; loop over mass bins
    legendStrs   = []
    legendColors = []
    
    if norm eq 1 then delta.val *= delta.norm_fac ; remove normalization
    xrange = [5e-3,0.8 * redshiftToAgeFlat(sP.redshift)]
    
    foreach massBin,massBins,i do begin
      wMB = where(yy ge massBin[0] and yy lt massBin[1], count)
      print,'['+str(i)+'] '+str(count)
      if count eq 0 then continue
      
      h = histogram( alog10(delta.val[wMB]), binsize=binSize_xx, $
                     min=xrangeLog[0], max=xrangeLog[1], loc=loc )
      
      cgPlot,10.0^(loc+binSize_xx*0.5),smooth(h/float(max(h)),sK),color=mbColors[i],/overplot
      
      legendStrs   = [legendStrs, textoidl('M_{halo} = ' + string(mean(massBin),format='(f4.1)'))]
      legendColors = [legendColors, mbColors[i]]
    endforeach
    
    legend,legendStrs,textcolors=legendColors,/top,/left
    
    ; plot B (normalized)    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=[0.001,2.0],pos=pos[1],/noerase,$
      xtitle=textoidl("( t_{1} - t_{2} ) / t_{2}"),ytitle="Fraction",$
      /xs,/ys,/ylog,yminor=0,xlog=logX,xminor=(logX ne 1),$
      title="(z="+string(sP.redshift,format='(f3.1)')+")"
      
    foreach massBin,massBins,i do begin
      wMB = where(yy ge massBin[0] and yy lt massBin[1], count)
      print,'['+str(i)+'] '+str(count)
      if count eq 0 then continue
      
      h = histogram( alog10(delta.val[wMB]/delta.norm_fac), binsize=binSize_xx, $
                     min=xrangeLog[0], max=xrangeLog[1], loc=loc )
      
      cgPlot,10.0^(loc+binSize_xx*0.5),smooth(h/float(max(h)),sK),color=mbColors[i],/overplot
    endforeach
    
    legend,legendStrs,textcolors=legendColors,/top,/left
    
  end_PS
    
  ; make 2d histogram difference
  ; ----------------------------  
  if n_elements(sP2) eq 0 then stop
  
  at2     = accretionTimes(sP=sP2)
  am2     = accretionMode(sP=sP2)
  galcat2 = galaxyCat(sP=sP2)
    
  if sP2.trMCPerCell gt 0 then types2 = ( galcat2.type[ replicate_var(galcat2.trMC_cc) ] )
  if sP2.trMCPerCell eq 0 then types2 = ( galcat2.type )
  if sP2.trMCPerCell lt 0 then types2 = ( galcat2.type[ replicate_var(galcat2.trVel_cc) ] )
  
  parentMass2 = galCatParentProperties(sP=sP2,galcat=galcat2,trRep=(sP2.trMCPerCell ne 0),/mass)
  
  ; plot (3) - 2d difference
  start_PS,sP.plotPath + 'accTimeDeltasVsMassDiff_' + plotStr + '.eps'
    cgText,0.5,0.96,textoidl("rVirFacs: t_1 = " + str(rVirFacs[laterInd]) + $
      ", t_2 = " + str(rVirFacs[earlierInd])) + ' (type='+gcType+')',alignment=0.5,/normal
                       
    ; calculate second histogram (same min/max/binsize/etc)
    delta2 = selectAccTimeDelta(sP=sP2, at=at2, am=am2, galcat=galcat2, earlierInd=earlierInd, $
                                laterInd=laterInd, gcType=gcType, accMode=accMode, norm=norm)
                              
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
    
    loadColorTable, 'brewer-redpurple' ; data
    tvim,h2_cmap,scale=0,pos=pos_plot,/c_map,/noframe,/noerase
          
    loadColorTable,'bw linear' ; axes/frame
    tvim,h2_cmap,/notv,pcharsize=!p.charsize,scale=0,clip=-1,$
      xtitle=xtitle,ytitle=ytitle,barwidth=0.75,lcharsize=!p.charsize-0.2,$
      xrange=xrange,yrange=yrange,xmargin=2.0,pos=pos_plot,/noerase,$
      xlog=logX,xticks=n_elements(xtickv)-1,xtickv=xtickv,xtickname=xtickname
     
   loadColorTable, 'brewer-redpurple' ; data
   cgColorbar,bottom=1,range=[fracMin,fracMax],position=pos_cbar,$
     /vertical,/right,divisions=ndivs,$ ;,ticknames=ticknames,ncolors=255
     title=textoidl("log ( f_{"+sP.plotPrefix+"} / f_{"+sP2.plotPrefix+"} )")

  end_PS
  
  stop
  
end
