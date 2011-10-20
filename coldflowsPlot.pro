; coldflowsPlot.pro
; cold flows - plots derived from intermediate datafiles
; dnelson oct.2011

; plot TmaxRedshift():

pro plotTmaxRedshift, res=res

  units = getUnits()

  if not keyword_set(res) then begin
    print,'Error: Must specific resolution set.'
    return
  endif
  
  ; config
  dataPath   = '/n/home07/dnelson/coldflows/thermhist.deriv/'
  plotPath   = '/n/home07/dnelson/coldflows/'

  bSH   = 0 ; include background subhalos
  halos = 1 ; find accretion onto main halos (do not set bSH,zWidth,redshiftsCut)
  
  redshifts     = [5.0,4.0,3.0,2.0,1.5,1.0,0.5,0.25,0.0] 
  redshiftNames = ['5','4','3','2','1.5','1.0','0.5','0.25','0']

  ; plot axes
  plotRedshifts = snapNumToRedshift(/all)
  times         = redshiftToAge(plotRedshifts)
  
  xrange = [5.0,0.0]
  yrange = [0.0,1.05]
    
  ; plot
  plotName = plotPath+'maxt_redshift_'+str(res)+'.eps'
  
  if (keyword_set(halos)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.halos.eps'
  if (keyword_set(bSH)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.bSH.eps'
    
  if (str(res[0]) eq 'all') then res = [512,256,128]    
    
  start_PS, plotName

  fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,$
           xtitle="Redshift",ytitle="Normalized Fraction",xs=9,/ys,ymargin=[4,3]
  universeage_axis,xrange,yrange
    
  if (n_elements(res) eq 1) then $
    fsc_text,0.06,0.94,/normal,str(res)+"^3",alignment=0.5
  
  ; calculate histogram
  for m=0,n_elements(redshifts)-1 do begin
    
    ; loop over each resolution
    for j=0,n_elements(res)-1 do begin    
    
      gadgetPath = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res[j])+'_20Mpc/Gadget/output/'
      targetSnap = redshiftToSnapnum(redshifts[m])
      
      ; get list of smoothly accreted gas particle ids
      sgIDs_Acc = findAccretedGas(res=res[j],bSH=bSH,targetSnap=targetSnap,halos=halos)    
      
      ; get time of maximum temperatures
      maxTempSnaps = maxTemperatures(res=res[j], /getSnap)
      maxTempSnaps = maxTempSnaps[sgIDs_Acc]
      maxTempRedshifts = plotRedshifts[maxTempSnaps]
      
      ; safeguard against degenerate size  
      if (n_elements(maxTempRedshifts) eq 1) then maxTempRedshifts = [-1.0,maxTempRedshifts]
      ;maxTempRedshifts = maxTempRedshifts[where(maxTempRedshifts gt 0.05)]
  
       ; overplot successive resolutions
      overplot = 1 ;1,1,1
      line     = j ;0,1,2
      thick    = !p.thick + 1 - 2*(j gt 0) ;3,1,1
  
      bin = 4.0 / sqrt(res[j])
  
      plothist,maxTempRedshifts,bin=bin,/peak,color=fsc_color(units.colors[m]),$
               overplot=overplot,line=line,thick=thick,psym=0

      ; plot legend entry
      if (j eq 0) then begin
        ypos = yrange[1]*0.92 - yrange[1]*0.05*m
        fsc_text,xrange[0]*0.95,ypos,"z = "+redshiftNames[m],$
                 alignment=0.0,color=fsc_color(units.colors[m]),charsize=!p.charsize-0.5
      endif

    endfor ;j
  endfor ;m
  
  end_PS

end

; plotNormScalarProd()

pro plotNormScalarProd, res=res

  ; config
  workingPath = '/n/home07/dnelson/coldflows/thermhist.deriv/'
  plotPath    = '/n/home07/dnelson/coldflows/'  
  
  bSH   = 0 ; include background subhalos
  halos = 1 ; find accretion onto main halos (do not set bSH,zWidth,redshiftsCut)
  
  normRel = 0 ; normalize (within each resolution) the same for hot/cold, which shows better both
              ; their relative dominance and the relative strengths of their peaks near 1.0, but
              ; shows worse the "relative" shape of the 1.0 peak vs the rest
  
  redshifts = [3.0,2.0,1.0,0.0]
  zstrs = 'z='+['3','2','1','0']

  xrange = [-1.25,1.25]
  yrange = [0.0,1.05]

  ; plot
  plotName = plotPath + 'normscalar_'+str(res)+'.eps'
  
  if (keyword_set(halos)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.halos.eps'
  if (keyword_set(bSH)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.bSH.eps'
  if (keyword_set(normRel)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.normRel.eps'
    
  if (str(res[0]) eq 'all') then res = [512,256,128]     
    
  start_PS, plotName
  !p.multi = [0,2,2]
  
  for m=0,n_elements(redshifts)-1 do begin
  
    ; plot frames
    if (m eq 0) then $
      fsc_plot,[0],[0],/nodata,ymargin=[1.0,2.0],xmargin=[7.0,0.0],$
               xrange=xrange,yrange=yrange,/xs,/ys,xtickname=replicate(' ',10),$
               ytickv=[0.25,0.50,0.75,1.0],yticks=3
    if (m eq 1) then $
      fsc_plot,[0],[0],/nodata,ymargin=[1.0,2.0],xmargin=[0.0,7.0],$
               xrange=xrange,yrange=yrange,/xs,/ys,$
               xtickname=replicate(' ',10),ytickname=replicate(' ',10)
    if (m eq 2) then $
      fsc_plot,[0],[0],/nodata,ymargin=[4.0,-1.0],xmargin=[7.0,0.0],$
               xrange=xrange,yrange=yrange,/xs,/ys,$
               ytickv=[0.0,0.25,0.50,0.75,1.0],yticks=5
    if (m eq 3) then $
      fsc_plot,[0],[0],/nodata,ymargin=[4.0,-1.0],xmargin=[0.0,7.0],$
               xrange=xrange,yrange=yrange,/xs,/ys,ytickname=replicate(' ',10)
             
    for j=0,n_elements(res)-1 do begin
  
      nsp = calcNormScalarProd(res=res[j],bSH=bSH,halos=halos,redshift=redshifts[m])

      ; flattten norm scalar products for global statistics
      normProd_Cold = flatten_list(nsp.normProd_Cold)
      normProd_Hot  = flatten_list(nsp.normProd_Hot)
      nsp = !NULL
      
      if (normProd_Cold eq !NULL or normProd_Hot eq !NULL) then begin
        print,'Flattened normProd is NULL, skipping plot (res='+str(res[j])+' '+zstrs[m]+').'
        continue
      endif
      
      ; check for degenerate sizes/entries
      if (n_elements(normProd_Cold) lt 2 or n_elements(normProd_Hot) lt 2) then begin
        print,'Too few elements, skipping plot (res='+str(res[j])+' '+zstrs[m]+').' 
        continue
      endif
      
      uniqHot  = n_elements(normProd_Hot[uniq(normProd_Hot,sort(normProd_Hot))])
      uniqCold = n_elements(normProd_Cold[uniq(normProd_Cold,sort(normProd_Cold))])
      
      if (uniqHot lt 2 or uniqCold lt 2) then begin
        print,'Not enough unique elements, skipping plot (res='+str(res[j])+' '+zstrs[m]+').' 
        continue
      endif
      
       ; overplot successive resolutions
      line     = j ;0,1,2
      thick    = !p.thick + 1 - (j gt 0) ;3,1,1
      
      ; cross-hatching of histograms
      fspacing = 0.15
      fill     = 0
      if (n_elements(res) eq 1) then $
        fill = 1
      
      ; histogram binsize
      bin = 0.1
      ;bin = 1.0 / sqrt(res[j])
      
      ; histogram y-normalization
      if (normRel) then begin
        normfac = max([histogram(normProd_Hot,bin=bin),histogram(normProd_Cold,bin=bin)])
        normfac = 1.0 / normfac
        peak = 0
      endif else begin
        normfac = !NULL
        peak = 1
      endelse

      ; plot histograms
      plothist,normProd_Cold,bin=bin,normfac=normfac,peak=peak,$
               fill=fill,/fline,forientation=-45,fspacing=fspacing,fcolor=fsc_color('blue'),$
               /overplot,color=fsc_color('blue'),line=line,thick=thick
      plothist,normProd_Hot,bin=bin,normfac=normfac,peak=peak,$
               fill=fill,/fline,forientation=45,fspacing=fspacing,fcolor=fsc_color('red'),$
               /overplot,color=fsc_color('red'),line=line,thick=thick
      
      fsc_text,xrange[0]*0.7,yrange[1]*0.8,zstrs[m],alignment=0.5   
       
    endfor ;j
  endfor ;m
  
  ; title and axes labels
  if (n_elements(res) eq 1) then $
    fsc_text,0.5,0.95,str(res)+"^3",alignment=0.5,/normal
    
  fsc_text,0.5,0.03,"Normalized Scalar Product",alignment=0.5,/normal
  fsc_text,0.03,0.5,"Normalized Fraction",alignment=0.5,orientation=90,/normal    
    
  !p.multi = 0
  end_PS

end

; plotAccretionRate(): plot the total smooth accretion rate (decomposed into hot and cold modes) over
;                      the whole simulation time

pro plotAccretionRate, res=res

  units = getUnits()

  if (not keyword_set(res)) then begin
    print,'Error: Must specify resolution.'
    return
  endif
  
  ; config
  dataPath    = '/n/hernquistfs1/dnelson/coldflows/thermhist/'
  plotPath    = '/n/home07/dnelson/coldflows/'
  
  bSH    = 0 ; include background subhalos
  halos  = 0 ; find accretion onto main halos (do not set bSH,zWidth,redshiftsCut)

  yScale  = 1 ; 0=normalized fraction, 1=mass/time/vol (K05)
  kerSize = 5

  ; plot axes
  plotRedshifts = snapNumToRedshift(/all)
  times         = redshiftToAge(plotRedshifts)
  
  xrange = [0.0,max(times)]

  ; plot
  plotName = plotPath+'accretion_rate_'+str(res)+'.eps'
  
  if (keyword_set(halos)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.halos.eps'
  if (keyword_set(bSH)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.bSH.eps'
    
  if (str(res[0]) eq 'all') then res = [512,256,128]        
    
  start_PS, plotName
    
    for j=0,n_elements(res)-1 do begin
      ; get accretion rate (code units mass)
      sa = findAccretionRate(res=res[j],bSH=bSH,halos=halos)
      
      smoothAccHot  = sa.smoothAccHot
      smoothAccCold = sa.smoothAccCold
    
      ; truncate first value (abnormally large)
      ind = min(where(smoothAccCold ne 0.0))
      smoothAccHot[ind]  = 0.01
      smoothAccCold[ind] = 0.01
    
      smoothAccTot  = smoothAccHot + smoothAccCold
    
      ; scale rate to physical units
      if (yScale eq 1) then begin
        ; mass / time / volume
        massFac = (units.UnitMass_in_g / units.Msun_in_g) ; (Msun)
        volFac  = (20.0)^3.0 ;Mpc^3 / h^3
        timeFac = (times - shift(times,1)) * 1e9 ;yr
          
        smoothAccTot  *= massFac / timeFac / volFac
        smoothAccHot  *= massFac / timeFac / volFac
        smoothAccCold *= massFac / timeFac / volFac
    
        yrange = [0.01,max(smoothAccTot)*1.05]
        ytitle = "Accretion Rate [h"+textoidl("^3")+" M  / yr / Mpc"+textoidl("^3")+"]"
      endif else begin
        smoothAccHot  /= max(smoothAccTot)
        smoothAccCold /= max(smoothAccTot)
        smoothAccTot  /= max(smoothAccTot)
        
        yrange = [0.0,1.05]
        ytitle = "Normalized Smooth Accretion Rate"
      endelse
      
      ; frame and text
      if (j eq 0) then begin
        fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/ylog,$
                 xtitle="Time [Gyr]",ytitle=ytitle,xs=9,/ys,ymargin=[4,3]
        redshiftNames = ['30','6','4','3','2','1.5','1.0','0.5','0.25','0']
        redshift_axis,xrange,yrange,/ylog,zTicknames=redshiftNames,/dotted
        
        ; sun symbol
        if (yScale eq 1) then $
          fsc_text,0.065,0.618,/normal,"!9!Z(6E)!X",font=-1,charthick=3,charsize=1.1 ;sunsymbol
        
        ; title
        if (n_elements(res) eq 1) then $
          fsc_text,0.06,0.94,/normal,str(res)+"^3",alignment=0.5
      endif
      
      ; overplot successive resolutions
      line     = j ;0,1,2
      thick    = !p.thick + 1 - (j gt 0) ;3,2,2
  
      fsc_plot,times,smooth(smoothAccCold,kerSize),color=fsc_color('blue'),/overplot,line=line,thick=thick
      fsc_plot,times,smooth(smoothAccHot,kerSize),color=fsc_color('red'),/overplot,line=line,thick=thick
      fsc_plot,times,smooth(smoothAccTot,kerSize),/overplot,line=line,thick=thick
      
    endfor ;j
  
  end_PS
  
end

; plotModeFracVsSubhaloMass(): plot the total smooth accretion rate (decomposed into hot and cold modes)
;                              as a function of total (DM+baryonic) halo (parent or subhalo) mass

pro plotModeFracVsSubhaloMass, res=res

  units = getUnits()

  if (not keyword_set(res)) then begin
    print,'Error: Must specify resolution.'
    return
  endif
  
  ; config
  dataPath    = '/n/hernquistfs1/dnelson/coldflows/thermhist/'
  plotPath    = '/n/home07/dnelson/coldflows/'

  kSP    = 1 ; keres+ 05 spacing
  bSH    = 0 ; include background subhalos
  halos  = 0 ; find accretion onto main halos (do not set bSH,zWidth,redshiftsCut)
  
  redshifts = [3.0,2.0,1.0,0.0]
  zstrs = 'z='+['3','2','1','0']
  
  if (kSP eq 1) then begin
    ; restriction: gas not resolved at keres+ 05 redshift spacing increment before target (single)
    redshiftsCut = [3.25,2.25,1.125,0.125]
  endif else begin
    ; restriction: gas not resolved in any snapshot previous to target
    redshiftsCut = [-1,-1,-1,-1]
  endelse

  ; plot config
  xrange = [8.5,12.5]
  yrange = [-0.08,1.08]
  
  plotsym,0 ;circle
  symsize = 0.3

  ; plot
  plotName = plotPath+'frac_shmass_'+str(res)+'.eps'
  
  if (keyword_set(halos)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.halos.eps'
  if (keyword_set(bSH)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.bSH.eps'
  if (keyword_set(kSP)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.kSP.eps'    
     
  if (str(res[0]) eq 'all') then res = [512,256,128]    
    
  start_PS, plotName
  !p.multi = [0,2,2]

  for m=0,n_elements(redshifts)-1 do begin
    targetSnap    = redshiftToSnapNum(redshifts[m])
    smoothCutSnap = redshiftToSnapnum(redshiftsCut[m])
    
    ; plot borders
    if (m eq 0) then $
      fsc_plot,[0],[0],ymargin=[1.0,2.0],xmargin=[7.0,0.0],xrange=xrange,yrange=yrange,/xs,/ys,$
               xtickname=replicate(' ',10)
    if (m eq 1) then $
      fsc_plot,[0],[0],ymargin=[1.0,2.0],xmargin=[0.0,7.0],xrange=xrange,yrange=yrange,/xs,/ys,$
               xtickname=replicate(' ',10),ytickname=replicate(' ',10)
    if (m eq 2) then $
      fsc_plot,[0],[0],ymargin=[4.0,-1.0],xmargin=[7.0,0.0],xrange=xrange,yrange=yrange,/xs,/ys
    if (m eq 3) then $
      fsc_plot,[0],[0],ymargin=[4.0,-1.0],xmargin=[0.0,7.0],xrange=xrange,yrange=yrange,/xs,/ys,$
               ytickname=replicate(' ',10)
  
    fsc_text,9.1,0.5,zstrs[m],alignment=0.5
    
    for j=0,n_elements(res)-1 do begin
      cf = calcColdFracVsSubhaloMass(res=res[j],halos=halos,bSH=bSH,$
                                     targetSnap=targetSnap,smoothCutSnap=smoothCutSnap)

    ; plot individual systems
    if (n_elements(res) eq 1) then $
      fsc_plot,cf.shmass,cf.coldfrac,psym=8,symsize=symsize,/overplot
      
    ; calculate median line
    massStep = round(100.0/sqrt(n_elements(cf.coldfrac)))/10.0
    ;massStep = 0.2
    massBins = (xrange[1]-xrange[0])/massStep
    massXPts = findgen(massBins)/massBins * (xrange[1]-xrange[0]) + xrange[0] + massStep/2.0
    
    medCold = fltarr(massBins) - 1
    stdCold = fltarr(massBins) - 1
    medCold[0:floor(n_elements(medCold)/2.0)] = 1.0 ;set default value high for first half
    
    for i=0,massBins-1 do begin
      w = where(cf.shmass ge xrange[0]+i*massStep and cf.shmass lt xrange[0]+(i+1)*massStep and $
                cf.coldfrac ne -1,count)
      if (count gt 0) then begin
        medCold[i] = median(cf.coldfrac[w])
        stdCold[i] = stddev(cf.coldfrac[w])
      endif
    endfor
    
     ; overplot successive resolutions
    line     = j ;0,1,2
    thick    = !p.thick + 1 - 2*(j gt 0) ;3,1,1
    smoothSize = 3
    
    ; plot smoothed median line
    w = where(medCold ne -1,count)
    
    if (count lt SmoothSize+1) then continue
      
    fsc_plot,massXPts[w],smooth(medCold[w],smoothSize),color=fsc_color('blue'),line=line,thick=thick,/overplot
    fsc_plot,massXPts[w],smooth(1.0-medCold[w],smoothSize),color=fsc_color('red'),line=line,$
             thick=thick,/overplot

    ; plot error visualization for j=0
    if (j eq 0) then begin
      w = where(stdCold gt 0,count)
      
      if (count lt smoothSize+1) then continue
      
      fsc_plot,massXPts[w],smooth(medCold[w]-stdCold[w],smoothSize),color=fsc_color('skyblue'),$
               line=line,thick=thick,/overplot
      fsc_plot,massXPts[w],smooth(medCold[w]+stdCold[w],smoothSize),color=fsc_color('skyblue'),$
               line=line,thick=thick,/overplot
    endif

    ; fit for medCold=0.5
    w = where(massXPts ge 10.5 and massXPts le 11.5 and medCold ne 1.0 and medCold ne -1.0)
    fit = linfit(massXPts[w],medCold[w])
    fitSHMass = (0.5 - fit[0]) / fit[1]
    
    ; plot best fit coldfrac=0.5
    fsc_plot,[fitSHMass,fitSHMass],[yrange[0],yrange[0]+(yrange[1]-yrange[0])/10.0],$
              line=line,thick=thick,/overplot

    print,'res='+str(res[j])+' cold frac 1/2 best fit subhalo mass = '+str(fitSHMass)

    endfor ;j
  endfor ;m
  
  if (n_elements(res) eq 1) then $
    fsc_text,0.5,0.95,str(res)+"^3",alignment=0.5,/normal
    
  objStr = "Subhalo"
  if (halos) then objStr = "Halo"
  
  fsc_text,0.5,0.05,"log ( Total "+objStr+" Mass )",alignment=0.5,/normal
  fsc_text,0.04,0.5,"Cold Mode Accretion Fraction",alignment=0.5,orientation=90,/normal
  
  !p.multi = 0
  end_PS

end

; plotTmaxHisto(): plot a histogram of maximum past temperature reached by all smoothly accreted gas 
;                  particles in a series of panels at various redshifts

pro plotTmaxHisto, res=res

  units = getUnits()

  if not keyword_set(res) then begin
    print,'Error: Must specific resolution set.'
    return
  endif
  
  ; config 
  dataPath   = '/n/home07/dnelson/coldflows/thermhist.deriv/'
  plotPath   = '/n/home07/dnelson/coldflows/'
  
  kSP    = 0 ; keres+ 05 spacing
  bSH    = 0 ; include background subhalos
  zWidth = 0 ; only if requesting a smoothCutSnap
  halos  = 0 ; find accretion onto main halos (do not set bSH,zWidth,redshiftsCut)
  
  yScale = 0 ; 0=normalized fraction, 1=K05, 2=K09
  
  redshifts = [5.0,4.0,3.0,2.0,1.5,1.0,0.5,0.25,0.0]
  zstrs = 'z='+['5.0','4.0','3.0','2.0','1.5','1.0','0.5','0.25','0.0']  
  
  if (kSP eq 1) then begin
    ; restriction: gas not resolved at keres+ 05 redshift spacing increment before target (single)
    redshiftsCut = [5.5,4.5,3.25,2.25,1.75,1.125,0.625,0.375,0.125]
  endif else begin
    ; restriction: gas not resolved in any snapshot previous to target
    redshiftsCut = [-1,-1,-1,-1,-1,-1,-1,-1,-1]
  endelse
  
  ; set plotName and open
  plotName = plotPath+'maxt_histo_allbins_'+str(res)+'.eps'

  if (keyword_set(halos)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.halos.eps'  
  if (keyword_set(bSH)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.bSH.eps'
  if (keyword_set(kSP) and not keyword_set(zWidth)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.kSP.eps'
  if (keyword_set(zWidth)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.zW.eps'
    
  if (str(res[0]) eq 'all') then res = [512,256,128]
    
  start_PS,plotName,xs=7,ys=6
  !p.multi = [0,3,3]

  ; plot config
  xrange = [3.8,7.2]
  
  if (yScale eq 0) then yrange = [0.0,1.10]
  if (yScale eq 1) then yrange = [0.0,50.0]
  if (yScale eq 2) then yrange = [0.0,20.0]
      
  bin = 0.1
   
  cs = !p.charsize
  pt = !p.thick
  !p.charsize += 1.0
      
  ; loop over each redshift (panel)
  for m=0,n_elements(redshifts)-1 do begin
    targetSnap    = redshiftToSnapnum(redshifts[m])
    smoothCutSnap = redshiftToSnapnum(redshiftsCut[m])
    
    ; loop over each resolution
    for j=0,n_elements(res)-1 do begin
    
      gadgetPath = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res[j])+'_20Mpc/Gadget/output/' 
      
      print,zstrs[m]+' res = ' + str(res[j]) + ' targetSnap = '+str(targetSnap)+' smoothCutSnap = ' + $
            str(smoothCutSnap)
     
      ; get list of smoothly accreted gas particle ids
      sgIDs_Acc = findAccretedGas(res=res[j],bSH=bSH,targetSnap=targetSnap,$
                                  smoothCutSnap=smoothCutSnap,zWidth=zWidth,halos=halos)
      
       ; get maximum temperatures
      maxTemps = maxTemperatures(res=res[j],zMin=redshifts[m])
      maxTemps = maxTemps[sgIDs_Acc]  
    
      if (n_elements(maxTemps) lt 2) then maxTemps = [-1.0,0.0]
      
      ; weighting
      if (yScale ne 0) then begin
        ; calc y-axis mass/time/vol constant
        masses = loadSnapshotSubset(gadgetPath,snapNum=targetSnap,partType='gas',field='mass')
        
        if (min(masses) ne max(masses)) then begin
          print,'ERROR'
          return
        end
        
        ; mass / time / volume
        masses *= (units.UnitMass_in_g / units.Msun_in_g) ; full array (Msun)
        volume = (20.0)^3.0 ;Mpc^3 / h^3
        
        if (keyword_set(kSP)) then $
          time = (redshiftToAge(redshifts[m]) - redshiftToAge(redshiftsCut[m])) * 1e9 ;yr
        if (not keyword_set(kSP)) then $
          time = (snapNumToAge(targetSnap) - snapNumToAge(targetSnap-1)) * 1e9 ;yr
        
        ; (log Tmax)^(-1)
        if (yScale eq 1) then $
          weight = masses / time / volume ;K05
        if (yScale eq 2) then $
          weight = masses / time / volume / maxTemps ;K09
        
        peak = 0
      endif else begin
        weight = !NULL
        peak = 1
      endelse
      
      ; overplot successive resolutions
      overplot = j gt 0 ;0,1,1
      line     = j ;0,1,2
      thick    = !p.thick + 1 - (j gt 0) ;3,2,2

      ; plot histograms
      if (m eq 0) then $
        plothist,maxTemps,bin=bin,peak=peak,weight=weight,xtitle="",ytitle="",$
                 xrange=xrange,yrange=yrange,xtickname=replicate(' ',10),/xs,/ys,$
                 xmargin=[7.0,-3.0],ymargin=[0.0,3.0],overplot=overplot,line=line,$
                 thick=thick;,ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4
      if (m eq 1) then $
        plothist,maxTemps,bin=bin,peak=peak,weight=weight,xtitle="",ytitle="",$
                 xrange=xrange,yrange=yrange,$
                 xtickname=replicate(' ',10),ytickname=replicate(' ',10),/xs,/ys,$
                 xmargin=[3.0,1.0],ymargin=[0.0,3.0],overplot=overplot,line=line,$
                 thick=thick
      if (m eq 2) then $
        plothist,maxTemps,bin=bin,peak=peak,weight=weight,xtitle="",ytitle="",$
                 xrange=xrange,yrange=yrange,$
                 xtickname=replicate(' ',10),ytickname=replicate(' ',10),/xs,/ys,$
                 xmargin=[-1.0,5.0],ymargin=[0.0,3.0],overplot=overplot,line=line,$
                 thick=thick
      if (m eq 3) then $
        plothist,maxTemps,bin=bin,peak=peak,weight=weight,xtitle="",ytitle="",$
                 xrange=xrange,yrange=yrange,xtickname=replicate(' ',10),/xs,/ys,$
                 xmargin=[7.0,-3.0],ymargin=[2.0,0.0],overplot=overplot,line=line,$
                 thick=thick;,ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4
      if (m eq 4) then $
        plothist,maxTemps,bin=bin,peak=peak,weight=weight,xtitle="",ytitle="",$
                 xrange=xrange,yrange=yrange,$
                 xtickname=replicate(' ',10),ytickname=replicate(' ',10),/xs,/ys,$
                 xmargin=[3.0,1.0],ymargin=[2.0,0.0],overplot=overplot,line=line,$
                 thick=thick
      if (m eq 5) then $
        plothist,maxTemps,bin=bin,peak=peak,weight=weight,xtitle="",ytitle="",$
                 xrange=xrange,yrange=yrange,$
                 xtickname=replicate(' ',10),ytickname=replicate(' ',10),/xs,/ys,$
                 xmargin=[-1.0,5.0],ymargin=[2.0,0.0],overplot=overplot,line=line,$
                 thick=thick
      if (m eq 6) then $
        plothist,maxTemps,bin=bin,peak=peak,weight=weight,xtitle="",ytitle="",$
                 xrange=xrange,yrange=yrange,$
                 /xs,/ys,xtickv=[4.0,5.0,6.0,7.0],xtickname=['4','5','6','7'],xticks=3,$
                 xmargin=[7.0,-3.0],ymargin=[5.0,-2.0],overplot=overplot,line=line,$
                 thick=thick
      if (m eq 7) then $
        plothist,maxTemps,bin=bin,peak=peak,weight=weight,xtitle="",ytitle="",$
                 xrange=xrange,yrange=yrange,/xs,/ys,ytickname=replicate(' ',10),$
                 xmargin=[3.0,1.0],ymargin=[5.0,-2.0],xticks=3,$
                 xtickv=[4.0,5.0,6.0,7.0],xtickname=['4','5','6','7'],overplot=overplot,$
                 line=line,thick=thick
      if (m eq 8) then $
        plothist,maxTemps,bin=bin,peak=peak,weight=weight,xtitle="",ytitle="",$
                 xrange=xrange,yrange=yrange,ytickname=replicate(' ',10),/xs,/ys,$
                 xmargin=[-1.0,5.0],ymargin=[5.0,-2.0],$
                 xtickv=[4.0,5.0,6.0,7.0],xtickname=['4','5','6','7'],xticks=3,overplot=overplot,$
                 line=line,thick=thick
      
      if (not overplot) then $
        fsc_text,xrange[1]*0.98,yrange[1]*0.82,zstrs[m],charsize=cs,alignment=1.0,color=fsc_color('orange')
  
    endfor ;j
    
  endfor ;m
  
  !p.charsize = cs
  !p.thick = pt
  
  ; title
  if (n_elements(res) eq 1) then $
    fsc_text,0.5,0.95,str(res)+"^3",alignment=0.5,/normal
    
  ; x-axis title
  fsc_text,0.5,0.04,"log (T"+textoidl("_{max}")+" [K])",alignment=0.5,/normal
  
  ; y-axis title
  if (yScale eq 0) then $
    fsc_text,0.04,0.5,"Normalized Gas Accretion Rate",alignment=0.5,orientation=90,/normal
  if (yScale eq 1) then begin
    fsc_text,0.04,0.5,"Gas Accretion Rate [h"+textoidl("^3")+" M  / yr / Mpc"+$
                      textoidl("^3")+"]",alignment=0.5,orientation=90,/normal
    fsc_text,0.035,0.62,/normal,"!9!Z(6E)!X",font=-1,charthick=3,charsize=1.1 ;sunsymbol
  endif
  if (yScale eq 2) then begin
    fsc_text,0.04,0.5,"Gas Accretion Rate [h"+textoidl("^3")+" M  / yr / Mpc"+$
                      textoidl("^3")+" / log(T"+textoidl("_{max}")+")]",alignment=0.5,orientation=90,/normal
    fsc_text,0.035,0.515,/normal,"!9!Z(6E)!X",font=-1,charthick=3,charsize=1.1 ;sunsymbol
  endif

  !p.multi = 0
  end_PS
  
end

; plotRhoTemp2D(): plot mass-weighted 2d histogram of thermal evolution tracks in (rho,temp) plane

pro plotRhoTemp2D, res=res

  units = getUnits()

  if not keyword_set(res) then begin
    print,'Error: Must specific resolution set.'
    return
  endif

  ; config
  plotPath  = '/n/home07/dnelson/coldflows/'

  bSH   = 0 ; include background subhalos
  halos = 1 ; find accretion onto main halos (do not set bSH,zWidth,redshiftsCut)
  tW    = 0 ; weight histogram bins by time (Gyr)
  
  nbins = 64
  
  redshifts   = [3.0,2.0,1.0,0.0]
  
  ; set plotName and open
  plotName = plotPath+'rhot_2dhisto_'+str(res)+'.nB='+str(nbins)+'.eps'

  if (keyword_set(halos)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.halos.eps'  
  if (keyword_set(bSH)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.bSH.eps' 
  if (keyword_set(tW)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.tW.eps'  
  
  start_PS, plotName, xs=7.0, ys=5
  !p.multi = [0,2,2]  

  ; color table
  loadct, 2, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
  
  xrange = [-2.0,8.0] ;can't change, must redo histo
  yrange = [3.0,7.0]  ;can't change, must redo histo
  clip = [0,100] ;percentage

  for m=0,n_elements(redshifts)-1 do begin
    redshift = redshifts[m]
  
    ; get (rho,temp) 2d histogram
    rt = calcRhoTemp2DHisto(res=res,bSH=bSH,halos=halos,zMin=redshift,nbins=nbins,tW=tW)

    ; convert code->solar mass units
    w = where(rt.h2rt ne 0)
    h2logm = fltarr(nbins,nbins)
    h2logm[w] = alog10( (units.UnitMass_in_g / units.Msun_in_g) * rt.h2rt[w] )

    ; plot
    if (m eq 0) then begin
      tvim,h2logm,clip=clip,position=[0.1,0.53,0.48,0.93],$;,/c_map
         stitle="",barwidth=0.0,lcharsize=!p.charsize,/rct,$
         xrange=xrange,yrange=yrange,xtickname=replicate(' ',10),yticks=3,ytickv=[4,5,6,7]
         ;,nodata=0,rgb_nodata=[1.0,1.0,1.0] ;display zeros as white not black
      fsc_text,xrange[1]*0.84,yrange[1]*0.92,"z=3",alignment=0.5
    endif
    if (m eq 1) then begin
      tvim,h2logm,clip=clip,position=[0.48,0.53,0.86,0.93],$;,/c_map
         stitle="",barwidth=0.0,lcharsize=!p.charsize,/rct,$
         xrange=xrange,yrange=yrange,xtickname=replicate(' ',10),ytickname=replicate(' ',10)
      fsc_text,xrange[1]*0.84,yrange[1]*0.92,"z=2",alignment=0.5
    endif
    if (m eq 2) then begin
      tvim,h2logm,clip=clip,position=[0.1,0.13,0.48,0.53],$;,/c_map
         stitle="",barwidth=0.0,lcharsize=!p.charsize,/rct,$
         xrange=xrange,yrange=yrange,yticks=4,ytickv=[3,4,5,6,7],xticks=6,xtickv=[-2,0,2,4,6,8]
      fsc_text,xrange[1]*0.84,yrange[1]*0.92,"z=1",alignment=0.5
    endif
    if (m eq 3) then begin
      tvim,h2logm,clip=clip,position=[0.48,0.13,0.86,0.53],$;,/c_map
         stitle="",barwidth=0.0,lcharsize=!p.charsize,/rct,$
         xrange=xrange,yrange=yrange,ytickname=replicate(' ',10),xticks=5,xtickv=[0,2,4,6,8]
      fsc_text,xrange[1]*0.84,yrange[1]*0.92,"z=0",alignment=0.5
    endif

  endfor ;m
  
  fsc_text,0.5,0.95,str(res)+"^3",alignment=0.5,/normal
  fsc_text,0.5,0.02,"log ("+textoidl("\rho / \rho_{crit}")+")",alignment=0.5,/normal
  fsc_text,0.04,0.5,"log (T [K])",alignment=0.5,orientation=90,/normal
  fsc_colorbar,position=[0.88,0.13,0.92,0.93],/vertical,/right,/reverse,$
               range=minmax(h2logm)
  
  !p.multi = 0
  end_PS
end

; plotTempTimeTracks(): plot tracks of temperature as a function of time/redshift
; TODO

pro plotTempTimeTracks, res=res

  units = getUnits()

  if not keyword_set(res) then begin
    print,'Error: Must specific resolution set.'
    return
  endif
  
  ; config
  dataPath  = '/n/hernquistfs1/dnelson/coldflows/thermhist/'
  plotPath  = '/n/home07/dnelson/coldflows/'
  
  bSH   = 1 ; include background subhalos
  halos = 0 ; NOT IMPLEMENTED TODO
  
  redshift      = 1.0 ;3,2,1,0
  critLogTemp   = 5.5
  numTracksEach = 20 
 
  ; load maxtemps to decide on tracks to save
  ;TODO
  
  ; find 100 tracks with Tmax<critLogTemp-1.0 (blue) and 100 with Tmax>critLogTemp+1.0 (red)
  wBlue = where(maxTemps lt critLogTemp-1.0, countBlue)
  wRed  = where(maxTemps gt critLogTemp+0.5, countRed)
  
  idBlue = randomu(seed,countBlue)
  idBlue = (sort(idBlue))[0:numTracksEach-1]
  idBlue = wBlue[idBlue]
  
  idRed = randomu(seed,countRed)
  idRed = (sort(idRed))[0:numTracksEach-1]
  idRed = wRed[idRed] 
  
  ; load thermhist data and save specified pid tracks
  ;TODO
  
  ; plot axes
  redshifts = snapNumToRedshift(/all)
  times     = redshiftToAge(redshifts)
  
  xrange = [0,times[redshiftToSnapnum(redshift)]]
  yrange = [1e1,4e6]
  
  ; temp vs. time/redshift tracks
  plotName = plotPath + 'temp.time.tracks.z='+string(redshift,format='(f4.2)')+'.eps'
  
  if (keyword_set(halos)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.halos.eps'
  if (keyword_set(bSH)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.bSH.eps'
    
  start_PS, plotName
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,$
             xtitle="Time [Gyr]",ytitle="Temperature [K]",xs=9,/ys,/ylog,ymargin=[4,3]
    redshift_axis,xrange,yrange,/ylog
    
    ; individual tracks
    for pID=0,numTracksEach-1 do begin
      fsc_plot,times,temp[idBlue[pID],*],/overplot,thick=0.1,color=fsc_color('blue')
      fsc_plot,times,temp[idRed[pID],*],/overplot,thick=0.1,color=fsc_color('red')
    endfor
    
    fsc_text,times[targetSnap]*0.8,3e2,"log "+textoidl("T_{max}")+" < 4.5",$
              alignment=0.5,color=fsc_color('blue')
    fsc_text,times[targetSnap]*0.8,1e2,"log "+textoidl("T_{max}")+" > 6.0",$
              alignment=0.5,color=fsc_color('red')
  end_PS
  
  ; rho/temp plane tracks
  start_PS, plotPath + 'rho.temp.tracks.'+string(redshift,format='(f4.2)')+'.eps'
    fsc_plot,[0],[0],/nodata,xrange=[0.0,4.0],yrange=[3.0,7.0],$
             xtitle="log "+textoidl("\rho / \rho_{crit}"),ytitle="log Temperature [K]",/xs,/ys
            
    plotsym,0
            
    ; individual tracks
    for pID=0,numTracksEach-1 do begin
      ;w = where(density[idBlue[pID],*] ne 0.0 and temp[idBlue[pID],*] ne 0.0, countGood)

      fsc_plot,alog10(density[idBlue[pID],*]/units.rhoCrit),$
               alog10(temp[idBlue[pID],*]),/overplot,$
               thick=0.5,psym=-8,symsize=0.2,color=fsc_color('blue')
      fsc_plot,alog10(density[idRed[pID],*]/units.rhoCrit),$
               alog10(temp[idRed[pID],*]),/overplot,$
               thick=0.5,psym=-8,symsize=0.2,color=fsc_color('red')
    endfor
    
  end_PS

end
