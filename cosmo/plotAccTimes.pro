; plotAccTimes.pro
; gas accretion project - past radial history of gas elements (plotting)
; dnelson aug.2014

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
  if logHist eq 1 then begin
    w = where(h2_cmap gt 0,count)
    if count gt 0 then h2_cmap[w] = alog10(h2_cmap[w])
  endif
  
  w = where(h2_cmap eq 0,count)
  if count gt 0 then h2_cmap[w] = min(h2_cmap[where(h2_cmap gt 0)])
    
  ; normalization by row?
  if byRow eq 1 then begin
    nRows = n_elements(h2_cmap[0,*])
    
    for i=0,nRows-1 do begin
      curRow = reform( h2_cmap[*,i] )
      curRow = fix( (curRow-min(curRow))/(max(curRow)-0.0) * (max-min) + min )
      h2_cmap[*,i] = curRow
    endfor
    
    return, h2_cmap
  endif
  
  ; normalization by column?
  if byCol eq 1 then begin
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

; plotAccTimeDeltas(): global, vs halo mass, and vs valmax

pro plotAccTimeDeltas
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  res      = 512
  runs     = ['feedback','tracer']
  redshift = 2.0
  
  ; index selection for difference ([1.0,0.75,0.5,0.25,0.15,0.05,0.01,first0.15,first1.0])
  rVirFacs = ['1.0 rvir','0.75 rvir','0.5 rvir','0.25 rvir','0.15 rvir','0.05 rvir',$
              '0.01 rvir','first 0.15 rvir','first 1.0 rvir']
  earlierInd = -1 ; -1 (0) ; 0 ; 2
  laterInd   = -2 ; -2 (4) ; 2 ; 5
  gcType     = 'all' ; all, gal (includes stars), gmem, inter
  accMode    = 'smooth' ; smooth_rec, clumpy_rec, stripped_rec, all
  plotVals   = ['temp','ent'] ; 
  
  ; plot config
  yrange     = [9.0,12.0]
  binSize_yy = 0.045
  cInd       = 1
  massBins   = list( [9.0,9.05], [9.45,9.55], [9.9,10.1], $
                     [10.4,10.6], [10.8,11.2], [11.3,11.7], [11.7,12.3] )
  sK         = 3
  plotUnused = 0 ; 0,1
      
  ; plot config
  if redshift eq 2.0 then xrange  = [0.04,4.0] ; delta time in Gyr
  if redshift eq 2.0 then xrangeT = [1.0,3.25] ; age of universe in Gyr
  if redshift eq 1.0 then xrange  = [0.04,6.0]
  if redshift eq 1.0 then xrangeT = [1.0,6.00]
  if redshift eq 0.0 then xrange  = [0.1,14.0]
  if redshift eq 0.0 then xrangeT = [1.0,14.0]
  
  xrangeLog  = alog10( xrange )
  binSize_xx = 0.03
  yrange_pdf = [0.001,0.1]
    
  ; plot config (vs maxvals)
  logHist    = 0 ; log 2d counts in each bin
  byRow      = 0 ; normalize 2d counts by row
  byCol      = 0 ; normalize 2d counts by column
  massBins2  = list( [10.5,10.6], [11.25,11.4], [11.7, 12.0] )
      
  ; load
  foreach run,runs do begin
    sP = simParams(res=res,run=run,redshift=redshift)
    saveFilename = sP.derivPath + 'binnedVals/accTimeDeltas.' + sP.savPrefix + str(sP.res) + $
      '.s' + str(sP.snap) + str(earlierInd) + str(laterInd) + '_' + gcType + '_' + accMode + '.sav'
      
    if ~file_test(saveFilename) then begin
      at     = accretionTimes(sP=sP)
      am     = accretionMode(sP=sP)
      mv     = maxVals(sP=sP)
      galcat = galaxyCat(sP=sP)
        
      types = ( galcat.type[ replicate_var(galcat.trMC_cc) ] )
      
      parentMass = galCatParentProperties(sP=sP,galcat=galcat,/trRep,/mass)
      
      delta = selectAccTimeDelta(sP=sP, at=at, am=am, galcat=galcat, earlierInd=earlierInd, $
                                 laterInd=laterInd, gcType=gcType, accMode=accMode)
                                 
      accHaloTvir = at.accHaloTvir[delta.w]
      parentMass = parentMass[delta.w]
      
      maxVals = { maxTemps    : mv.maxTemps[delta.w]    ,$
                  maxTempTime : mv.maxTempTime[delta.w] ,$
                  maxEnt      : mv.maxEnt[delta.w]      ,$
                  maxDens     : mv.maxDens[delta.w]     ,$
                  atEarly     : at.accTime[earlierInd,delta.w] ,$
                  atLate      : at.accTime[laterInd,delta.w]    }

      save,parentMass,delta,maxVals,accHaloTvir,filename=saveFilename
      print,'Saved TEMP file: ['+saveFilename+']'
    endif else begin
      print,'RESTORE: ['+saveFilename+']'
      restore,saveFilename
    endelse
    
    ; add to keeper
    atd = mod_struct( atd, run, delta )
    pm  = mod_struct( pm, run, parentMass )
    mvs = mod_struct( mvs, run, maxVals )
    htv = mod_struct( htv, run, accHaloTvir )
    sPs = mod_struct( sPs, run, sP )
    
  endforeach ;run
  
  ; plot setup
  simNames  = []
  simColors = []
  plotStrG  = ''
  
  for i=0,n_elements(runs)-1 do begin
    simNames  = [simNames, sPs.(i).simName]
    simColors = [simColors, sPs.(i).colors[cInd]]
    plotStrG  += sPs.(i).savPrefix
  endfor
  
  set_plot,'ps'
  mbColors = reverse( sampleColorTable('blue-red2', n_elements(massBins), bounds=[0.0,1.0]) )
  
  plotStrG += str(sPs.(0).res) + '_' + str(earlierInd) + '_' + str(laterInd) + '_gc-' + gcType + $
              "_mode-" + accMode + '_z' + str(fix(redshift))
   
  if plotUnused then begin
  foreach run,runs,j do begin
    ; plot (1) - 2d histogram
    plotStr = sPs.(j).savPrefix + str(sPs.(j).res) + '_' + str(earlierInd) + '_' + str(laterInd) + $
              '_gc-' + gcType + "_mode-" + accMode
    
    start_PS,sPs.(j).plotPath + 'accTimeDeltasVsMass_' + plotStr + '.eps'
    
      ; add to plot
      yy = alog10( 10.0^(pm.(j)) * units.hubbleParam )
      
      weights = fltarr(n_elements(atd.(j).val)) + 1.0
      
      xx = alog10( atd.(j).val )
      
      h2 = hist_nd_weight( transpose( [[xx],[yy]] ), weight=weights, [binSize_xx,binSize_yy], $
          min=[xrangeLog[0]-binSize_xx*0.5,yrange[0]-binSize_yy*0.5],$
          max=[xrangeLog[1]+binSize_xx*0.49,yrange[1]+binSize_yy*0.49])
      
      ; colormap (each halo mass row individually scaled)
      h2_cmap = colorMapAccTime(h2,logHist=logHist,byRow=1,byCol=0)
      
      h2s = mod_struct( h2s, run, h2 )
      
      ; plot
      pos = [0.14,0.14,0.92,0.9]
      
      loadColorTable, 'helix', /reverse ; data
      tvim,h2_cmap,scale=0,pos=pos,/c_map,/noframe,/noerase
            
      xtitle = textoidl("( t_{gal} - t_{halo} ) [Gyr]")
      ytitle = textoidl("M_{halo}")
            
      loadColorTable,'bw linear' ; axes/frame
      tvim,h2_cmap,/notv,pcharsize=!p.charsize,scale=0,clip=-1,$
        xtitle=xtitle,ytitle=ytitle,barwidth=0.75,lcharsize=!p.charsize,$
        xrange=xrange,yrange=yrange,xmargin=2.0,pos=pos,/noerase,/xlog,xminor=0
          
    end_PS
  
    ; plot (2) - 1d profiles        
    start_PS,sPs.(j).plotPath + 'accTimeDeltasMassBins_' + plotStr + '.eps', /extrabig

      ; (no normalization)
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_pdf,pos=pos,$
        xtitle=textoidl("( t_{gal} - t_{halo} ) [Gyr]"),ytitle="PDF",$
        /xs,/ys,/ylog,yminor=0,/xlog,xminor=0,title=""
      
      yy = alog10( 10.0^(pm.(j)) * units.hubbleParam )
      
      ; loop over mass bins
      legendStrs   = []
      legendColors = []
      
      foreach massBin,massBins,i do begin
        wMB = where(yy ge massBin[0] and yy lt massBin[1], count)
        if count eq 0 then continue
        
        h = histogram( alog10(atd.(j).val[wMB]), binsize=binSize_xx, $
                       min=xrangeLog[0], max=xrangeLog[1], loc=loc )
        
        loc = loc + binSize_xx*0.5
        loc = 10.0^loc
        
        cgPlot,loc,smooth(h/float(total(h)),sK),$
          color=mbColors[i],/overplot
        
        legendStrs   = [legendStrs, textoidl('M_{halo} = ' + string(mean(massBin),format='(f4.1)'))]
        legendColors = [legendColors, mbColors[i]]
      endforeach
      
      legend,legendStrs,textcolors=legendColors,/top,/left
      
    end_PS
    
  endforeach ; j,runs
  
  ; plot (1b) - 2d difference of (delta,halomass)
  start_PS,sPs.(0).plotPath + 'accTimeDeltas2DDiff_' + plotStrG + '.eps'
                
    ; colormap
    fracMin = -0.6 ; 0.25
    fracMax = 0.6 ; 4.0
    
    h2_cmap = colorMapAccTimeDelta(h2s.(0),h2s.(1),fracMin=fracMin,fracMax=fracMax)

    ; add difference of two histograms to plot
    pos_plot = [0.14,0.14,0.84,0.9]
    pos_cbar = [0.85,0.14,0.89,0.9]
    ndivs = 5
    
    loadColorTable, 'brewer-redpurple' ; data
    tvim,h2_cmap,scale=0,pos=pos_plot,/c_map,/noframe,/noerase
          
    loadColorTable,'bw linear' ; axes/frame
    tvim,h2_cmap,/notv,pcharsize=!p.charsize,scale=0,clip=-1,$
      xtitle=xtitle,ytitle=ytitle,barwidth=0.75,lcharsize=!p.charsize,$
      xrange=xrange,yrange=yrange,xmargin=2.0,pos=pos_plot,/noerase,/xlog,xminor=0
     
    loadColorTable, 'brewer-redpurple' ; data
    cgColorbar,bottom=1,range=[fracMin,fracMax],position=pos_cbar,$
      /vertical,/right,divisions=ndivs,$
      title=textoidl("log ( f_{"+sPs.(0).plotPrefix+"} / f_{"+sPs.(1).plotPrefix+"} )")

  end_PS
    
  endif ;plotUnused
  
  ; plot (3) - combined 1d profiles in mass bins
  start_PS,sPs.(0).plotPath + 'accTimeDeltasMassBinsComp_' + plotStrG + '.eps'
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_pdf,pos=pos,$
      xtitle=textoidl("( t_{gal} - t_{halo} ) [Gyr]"),ytitle="PDF",$
      /xs,/ys,/ylog,yminor=0,/xlog,xminor=0
  
    foreach run,runs,j do begin
      ; loop over mass bins
      legendStrs   = []
      legendColors = []
      
      yy = alog10( 10.0^(pm.(j)) * units.hubbleParam )
            
      foreach massBin,massBins,i do begin
        wMB = where(yy ge massBin[0] and yy lt massBin[1], count)
        if count eq 0 then continue
        
        h = histogram( alog10(atd.(j).val[wMB]), binsize=binSize_xx, $
                       min=xrangeLog[0], max=xrangeLog[1], loc=loc )
        
        loc = loc + binSize_xx*0.5
        loc = 10.0^loc
        
        cgPlot,loc,smooth(h/float(total(h)),sK),color=mbColors[i],/overplot,line=j
        
        legendStrs   = [legendStrs, textoidl('M_{halo} = ' + string(mean(massBin),format='(f4.1)'))]
        legendColors = [legendColors, mbColors[i]]
      endforeach
      
    endforeach ; j,runs
    
    legend,legendStrs,textcolors=legendColors,/top,/left
    legend,simNames,textcolor=simColors,linestyle=[0,1],/top,/right
  end_PS
  
  ; plot (4) - combined 1d global profiles
  start_PS,sPs.(0).plotPath + 'accTimeDeltasComp_' + plotStrG + '.eps'
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_pdf,pos=pos,$
      xtitle=textoidl("( t_{gal} - t_{halo} ) [Gyr]"),ytitle="PDF",$
      /xs,/ys,/ylog,yminor=0,/xlog,xminor=0
  
    foreach run,runs,j do begin
    
      h = histogram( alog10(atd.(j).val), binsize=binSize_xx, $
                     min=xrangeLog[0], max=xrangeLog[1], loc=loc )
        
      loc = loc + binSize_xx*0.5
      loc = 10.0^loc
      
      cgPlot,loc,smooth(h/float(total(h)),sK),$
        color=simColors[j],/overplot
      
    endforeach ; j,runs
    
    legend,simNames,textcolor=simColors,/top,/left
  end_PS
  
  ; plot (5) - delta time vs maxvals
  foreach plotVal, plotVals do begin
  
    plotStr = plotStrG+'_'+str(logHist)+str(byRow)+str(byCol)+'.eps'
            
    if plotVal eq 'temp' then begin
      yrange     = [4.0,8.0]
      yrangeNorm = [-1.0,1.0] ;[0.0315,30.0]
      ytitle     = textoidl("T_{max} [_{ }log K_{ }]")
      ytitleNorm = textoidl("log ( T_{max} / T_{vir} )")
    endif
    if plotVal eq 'ent' then begin
      yrange     = [5.2,9.6]
      yrangeNorm = [-0.5,2.5] ;[0.0315,500.0]
      ytitle     = "S [cgs]"
      ytitleNorm = textoidl("log ( S / S_{200} )")
    endif
  
    ; plot (5a)
    start_PS,sP.plotPath + 'accTimeDeltasVs_'+plotVal+'Max_' + plotStrG + '.eps', $
      xs=4.6, ys=3.8*n_elements(runs)
        
    pos = plot_pos(rows=n_elements(runs),cols=1,/gap)
    
    foreach run,runs,j do begin
    
      ; bin
      weights = fltarr(n_elements(atd.(j).val)) + 1.0 ; uniform
      
      parentS200 = codeMassToVirEnt( pm.(j), sP=sPs.(j), /log ) ; only used for 'ent'
      
      xrangeLog2 = alog10( [0.1,2.0] )
      
      xx = alog10( atd.(j).val )
      
      if plotVal eq 'temp' then $
        yy = alog10( 10.0^mvs.(j).maxTemps / 10.0^htv.(j) )
      if plotVal eq 'ent' then $
        yy = alog10( 10.0^mvs.(j).maxEnt / 10.0^parentS200 )
          
      h2 = hist_nd_weight( transpose( [[xx],[yy]] ), weight=weights, [binSize_xx,binSize_yy], $
          min=[xrangeLog2[0]-binSize_xx*0.50,yrangeNorm[0]-binSize_yy*0.50],$
          max=[xrangeLog2[1]+binSize_xx*0.49,yrangeNorm[1]+binSize_yy*0.49])
      
      ; colormap (each halo mass row individually scaled)
      h2_cmap = colorMapAccTime(h2,logHist=logHist,byRow=byRow,byCol=byCol)
        
      ; plot
      loadColorTable, 'helix', /reverse ; data
      tvim,h2_cmap,scale=0,pos=pos[j]+[0.06,0,0,0],/c_map,/noframe,/noerase
            
      loadColorTable,'bw linear' ; axes/frame
      cgPlot,[0],[0],/nodata,xrange=10.0^xrangeLog2,yrange=yrangeNorm,title="",$
        xtitle=textoidl("( t_{gal} - t_{halo} ) [Gyr]"),$
        ytitle=ytitleNorm,/xs,/ys,/xlog,xminor=0,pos=pos[j]+[0.06,0,0,0],/noerase
        
      ; lower envelope line
      xx = [0.3,1.8]
      yy = xx*0.5 - 0.8
      
      cgPlot,xx,yy,line=2,color=cgColor('black'),thick=!p.thick,/overplot
        
      legend,[sPs.(j).simName],textcolors=[sPs.(j).colors[cInd]],/top,/left 
      
    endforeach ; runs
    
    end_PS
    
    ; plot (6) - split into 6 halo mass bins, one plot per run
    start_PS,sP.plotPath + 'accTimeDeltasVs_'+plotVal+'Max_'+'MassBin_'+plotStrG+'.eps', $
      xs=13, ys=7.4
  
    pos = plot_pos(rows=n_elements(runs),cols=3,/gap)
    
    foreach run,runs,j do begin
  
      ; plot each run/normalized or not combination
      foreach massBin,massBins2,k do begin
      
        ; get accTimeDelta and maximum values                                   
        massBinStr = textoidl(string(massBin[0],format='(f4.1)') + " < M_{halo} < " + $
                              string(massBin[1],format='(f4.1)'))

        weights = fltarr(n_elements(atd.(j).val)) + 1.0 ; uniform
      
        ; select
        zz = alog10( 10.0^(pm.(j)) * units.hubbleParam )
            
        wMB = where(zz ge massBin[0] and zz lt massBin[1], count)
        if count eq 0 then continue
        
        ; plot: delta_t vs tmax/tvir
        parentS200 = codeMassToVirEnt( pm.(j), sP=sPs.(j), /log ) ; only used for 'ent'
      
        xx = alog10( atd.(j).val[wMB] )
      
        if plotVal eq 'temp' then $
          yy = alog10( 10.0^mvs.(j).maxTemps[wMB] / (10.0^htv.(j))[wMB] )
        if plotVal eq 'ent' then $
          yy = alog10( 10.0^mvs.(j).maxEnt[wMB] / 10.0^parentS200[wMB] )
        
        h2 = hist_nd_weight( transpose( [[xx],[yy]] ), weight=weights, [binSize_xx,binSize_yy], $
            min=[xrangeLog[0]-binSize_xx*0.50,yrangeNorm[0]-binSize_yy*0.50],$
            max=[xrangeLog[1]+binSize_xx*0.49,yrangeNorm[1]+binSize_yy*0.49], rev=ri)

        ; colormap (each halo mass row individually scaled)
        h2_cmap = colorMapAccTime(h2,logHist=logHist,byRow=byRow,byCol=byCol)
        
        ; plot
        loadColorTable, 'helix', /reverse ; data
        tvim,h2_cmap,scale=0,pos=pos[3*j+k],/c_map,/noframe,/noerase
            
        loadColorTable,'bw linear' ; axes/frame
        cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeNorm,title="",$
          xtitle=textoidl("( t_{gal} - t_{halo} ) [Gyr]"),$
          ytitle=ytitleNorm,/xs,/ys,/xlog,xminor=0,pos=pos[3*j+k],/noerase
          
        legend,massBinStr,textcolors=cgColor('black'),/top,/left
        legend,[sPs.(j).simName],textcolors=[sPs.(j).colors[cInd]],/bottom,/right
      
      endforeach ; massBins2
        
    endforeach ; runs
  
    end_PS
    
    ; plot (7) - difference of first two runs
    start_PS,sPs.(0).plotPath + 'accTimeDiff2DVs_'+plotVal+'Max_' + plotStr, xs=5.6, ys=4.2
    
      hh = {} ; tmax/tvir
    
      foreach run,runs,j do begin
      
        ; histogram: tmax/tvir
        parentS200 = codeMassToVirEnt( pm.(j), sP=sPs.(j), /log ) ; only used for 'ent'
        xx = alog10( atd.(j).val )
        
        if plotVal eq 'temp' then $
          yy = alog10( 10.0^mvs.(j).maxTemps / 10.0^htv.(j) )
        if plotVal eq 'ent' then $
          yy = alog10( 10.0^mvs.(j).maxEnt / 10.0^parentS200 )
    
        h2 = hist_nd_weight( transpose( [[xx],[yy]] ), weight=weights, [binSize_xx,binSize_yy], $
            min=[xrangeLog[0]-binSize_xx*0.50,yrangeNorm[0]-binSize_yy*0.50],$
            max=[xrangeLog[1]+binSize_xx*0.49,yrangeNorm[1]+binSize_yy*0.49])
            
        hh = mod_struct( hh, run, h2 )
        
      endforeach ;runs
      
      ; different histograms and colormap
      fracMin = -0.6 ; 0.25
      fracMax = 0.6 ; 4.0
      totalNorm = 1
      
      hh_cmap = colorMapAccTimeDelta(hh.(0),hh.(1),fracMin=fracMin,fracMax=fracMax,totalNorm=totalNorm)

      ; add difference of two histograms to plot
      pos_plot = [0.16,0.14,0.82,0.9]
      pos_cbar = [0.84,0.14,0.88,0.9]
      ndivs = 5
      
      loadColorTable, 'brewer-redpurple' ; data
      tvim,hh_cmap,scale=0,pos=pos_plot,/c_map,/noframe,/noerase
            
      loadColorTable,'bw linear' ; axes/frame
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeNorm,title="",$
        xtitle=textoidl("( t_{gal} - t_{halo} ) [Gyr]"),ytitle=ytitleNorm,$
        /xs,/ys,/xlog,xminor=0,pos=pos_plot,/noerase
       
      ; colorbar
      loadColorTable, 'brewer-redpurple' ; data
      cgColorbar,bottom=1,range=[fracMin,fracMax],position=pos_cbar,$
        /vertical,/right,divisions=ndivs,$
        title=textoidl("log ( f_{"+sPs.(0).plotPrefix+"} / f_{"+sPs.(1).plotPrefix+"} )")
      
    end_PS
        
  endforeach ; plotVals
  
  ; plot (8) - correlation between TmaxTime and radial crossing times
  start_PS,sP.plotPath + 'accTime_vs_TmaxTime_' + plotStrG + '.eps', $
    xs=4.8, ys=3.8*n_elements(runs)
      
  pos = plot_pos(rows=n_elements(runs),cols=1,/gap)
  
  foreach run,runs,j do begin
  
    ; bin
    weights = fltarr(n_elements(atd.(j).val)) + 1.0 ; uniform
    
    xx = reform( mvs.(j).atEarly ) ; scalefac
    yy = mvs.(j).maxTempTime ; scalefac
    
    xx = redshiftToAgeFlat( 1/xx-1 )
    yy = redshiftToAgeFlat( 1/yy-1 )
    
    h2 = hist_nd_weight( transpose( [[xx],[yy]] ), weight=weights, [binSize_xx,binSize_yy], $
        min=[xrangeT[0]-binSize_xx*0.50,xrangeT[0]-binSize_yy*0.50],$
        max=[xrangeT[1]+binSize_xx*0.49,xrangeT[1]+binSize_yy*0.49])
    
    ; colormap (each halo mass row individually scaled)
    h2_cmap = colorMapAccTime(h2,logHist=logHist,byRow=0,byCol=0)
      
    ; plot (a)
    loadColorTable, 'helix', /reverse ; data
    tvim,h2_cmap,scale=0,pos=pos[j]+[0.06,0,-0.02,0],/c_map,/noframe,/noerase
          
    loadColorTable,'bw linear' ; axes/frame
    cgPlot,[0],[0],/nodata,xrange=xrangeT,yrange=xrangeT,title="",$
      xtitle="",ytitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),$
      /xs,/ys,pos=pos[j]+[0.06,0,-0.02,0],/noerase
      
    ; guiding lines
    nRes = 20
    dynRedshifts = linspace(6.0,sPs.(j).redshift,nRes)
    oneToOne = redshiftToAgeFlat(dynRedshifts)
    cgPlot,oneToOne,oneToOne,color='yellow',line=2,/overplot
    
    ;cgPlot,[oneToOne[1]-1,oneToOne[1]],[oneToOne[1],oneToOne[1]],color=cgColor('orange'),line=1,/overplot
    ;cgPlot,[oneToOne[1],oneToOne[1]],[oneToOne[1]-1,oneToOne[1]],color=cgColor('orange'),line=1,/overplot

    ; dynTime only depends on redshift
    dynTime = fltarr(nRes)
    for q=0,nRes-1 do begin
      sis_dm = sis_profile(1.0,mass=1.0,redshift=dynRedshifts[q])    
      sis_gas = sis_gas_profile(mass_hot=50.0,sis_dm=sis_dm)
      dynTime[q] = sis_gas.dynTime_halo
    endfor
    
    cgPlot,oneToOne,oneToOne+0.5*dynTime,color='orange',line=1,/overplot

    if redshift eq 2.0 then textPos = [2.0,2.4]
    if redshift eq 1.0 then textPos = [3.2,3.8]
    if redshift eq 0.0 then textPos = [5.0,5.0] ; todo
    
    cgText,textPos[0],textPos[1],textoidl("t_{ }(T_{max}) = t_{halo} + t_{dyn }/2"),$
      orientation=40.0,color='orange',alignment=0.5
    
    cgPlot,[0],[0],/nodata,xrange=xrangeT,yrange=xrangeT,title="",$
      xtitle=textoidl("t_{halo} [Gyr]"),ytitle=textoidl("t_{ }( T_{max} ) [Gyr]"),$
      /xs,/ys,pos=pos[j]+[0.06,0,-0.02,0],/noerase
    
    legend,[sPs.(j).simName],textcolors=[sPs.(j).colors[cInd]],/top,/left
    
  endforeach ; runs
  
  end_PS
  
  ; plot (8) - correlation between TmaxTime and radial crossing times
  start_PS,sP.plotPath + 'accTimeRatio_vs_TmaxTime_' + plotStrG + '.eps', $
    xs=4.8, ys=3.8*n_elements(runs)
      
  pos = plot_pos(rows=n_elements(runs),cols=1,/gap)
  
  foreach run,runs,j do begin
  
    ; bin
    weights = fltarr(n_elements(atd.(j).val)) + 1.0 ; uniform
    
    xx = reform( mvs.(j).atEarly ) ; scalefac
    yy = mvs.(j).maxTempTime ; scalefac
    
    xx = redshiftToAgeFlat( 1/xx-1 )
    yy = redshiftToAgeFlat( 1/yy-1 )
    
    yy /= xx ; ratio
    
    binSize_yy2 = 0.01
    xrange2 = [0.8,1.2]
    
    h2 = hist_nd_weight( transpose( [[xx],[yy]] ), weight=weights, [binSize_xx,binSize_yy2], $
        min=[xrangeT[0]-binSize_xx*0.50,xrange2[0]-binSize_yy2*0.50],$
        max=[xrangeT[1]+binSize_xx*0.49,xrange2[1]+binSize_yy2*0.49])
    
    ; colormap (each halo mass row individually scaled)
    h2_cmap = colorMapAccTime(h2,logHist=logHist,byRow=0,byCol=0)
      
    ; plot (a)
    loadColorTable, 'helix', /reverse ; data
    tvim,h2_cmap,scale=0,pos=pos[j]+[0.06,0,-0.02,0],/c_map,/noframe,/noerase
          
    loadColorTable,'bw linear' ; axes/frame
    cgPlot,[0],[0],/nodata,xrange=xrangeT,yrange=xrange2,title="",$
      xtitle="",ytitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),$
      /xs,/ys,pos=pos[j]+[0.06,0,-0.02,0],/noerase
      
    ; guiding lines
    nRes = 20
    dynRedshifts = linspace(6.0,sPs.(j).redshift,nRes)
    oneToOne = redshiftToAgeFlat(dynRedshifts)
    cgPlot,oneToOne,replicate(1.0,nRes),color='yellow',line=2,/overplot
    
    ;cgPlot,[oneToOne[1]-1,oneToOne[1]],[oneToOne[1],oneToOne[1]],color=cgColor('orange'),line=1,/overplot
    ;cgPlot,[oneToOne[1],oneToOne[1]],[oneToOne[1]-1,oneToOne[1]],color=cgColor('orange'),line=1,/overplot

    ; dynTime only depends on redshift
    dynTimeRatio = fltarr(nRes)
    for q=0,nRes-1 do begin
      sis_dm = sis_profile(1.0,mass=1.0,redshift=dynRedshifts[q])    
      sis_gas = sis_gas_profile(mass_hot=50.0,sis_dm=sis_dm)
      dynTimeRatio[q] = sis_gas.dynTime_halo / oneToOne[q]
    endfor
    
    cgPlot,oneToOne,1.0+0.5*dynTimeRatio,color='orange',line=1,/overplot

    textPos = [max(dynRedshifts)*0.5,0.97]
    
    cgText,textPos[0],textPos[1],textoidl("t_{ }(T_{max}) = t_{halo} + t_{dyn }/2"),$
      orientation=40.0,color='orange',alignment=0.5
    
    cgPlot,[0],[0],/nodata,xrange=xrangeT,yrange=xrange2,title="",$
      xtitle=textoidl("t_{halo} [Gyr]"),ytitle=textoidl("t_{ }( T_{max} ) / t_{halo}"),$
      /xs,/ys,pos=pos[j]+[0.06,0,-0.02,0],/noerase
    
    legend,[sPs.(j).simName],textcolors=[sPs.(j).colors[cInd]],/top,/left
    
  endforeach ; runs
  
  end_PS
    
  stop
  
end
