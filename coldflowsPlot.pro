; coldflowsPlot.pro
; cold flows - plots derived from intermediate datafiles
; dnelson nov.2011

; ---------------------------------------------------------------------------------
; ---------------------------------------------------------------------------------
; NOTE: Nothing in this file is still in use and may be out of date and/or buggy.
;       Only still around for ideas when working on cosmo* codebase.
; ---------------------------------------------------------------------------------
; ---------------------------------------------------------------------------------

; plotRedshiftSpacings(): plot output snapshot spacings in redshift for K05, K09, V11

pro plotRedshiftSpacings

  ; config
  plotPath = '/n/home07/dnelson/coldflows/'
  plotName = plotPath+'redshift_spacings.eps'

  xrange = [0.9,5.0]
  yrange = [0.0,1.0]

  start_PS, plotName
  
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,ytickname=replicate(' ',10),$
             xtitle="Redshift",ytitle="Output Snapshot Spacing",/xs,/ys,/xlog,$
             xtickv=[1.0,2.0,3.0,4.0],xticks=3,xtickname=['1','2','3','4']
  
    ; oplot dotted lines
    foreach redshift,list(1.0,2.0,3.0,4.0) do $
      fsc_plot,[redshift,redshift],[0.0,1.0],line=1,thick=!p.thick-1,/overplot
  
    ; K05
    k05_z = findgen(26)+5 ;dz = 1.0, z > 5
    k05_z = [3.5,4.0,4.5,k05_z] ;dz=0.5, 3.5<z<5
    k05_z = [1.0,1.25,1.50,1.75,2.0,2.25,2.50,2.75,3.0,3.25,k05_z] ;dz=0.25, 1<z<3.5
    k05_z = [0.0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,k05_z] ;dz=0.125, 0<z<1
  
    for i=0,n_elements(k05_z)-1 do $
      fsc_plot,[k05_z[i],k05_z[i]],[0.1,0.3],thick=!p.thick+1,color=fsc_color('crimson'),/overplot
    fsc_text,xrange[0]*1.2,0.32,"K05",alignment=0.5,color=fsc_color('crimson')
  
    ; K09 load scalefactors
    pt = {a:0.0}
    pts = loadCSV(0,plotPath+'output.spacings.k09.txt',pt)
    
    k09_z = 1.0 / pts.a - 1.0
    
    for i=0,n_elements(k09_z)-1 do $
      fsc_plot,[k09_z[i],k09_z[i]],[0.4,0.6],thick=!p.thick+1,color=fsc_color('slate blue'),/overplot
    fsc_text,xrange[0]*1.2,0.62,"K09",alignment=0.5,color=fsc_color('slate blue')
    
    ; V11 load redshifts
    v11_z = snapNumToRedshift(/all)
    
    for i=0,n_elements(v11_z)-1 do $
      fsc_plot,[v11_z[i],v11_z[i]],[0.7,0.9],thick=!p.thick+1,color=fsc_color('forest green'),/overplot
    fsc_text,xrange[0]*1.2,0.92,"V11",alignment=0.5,color=fsc_color('forest green')
    
  end_PS

end

; plotHaloMassesVirTemp(): plot halo/subhalo histograms of mass and virial temp

pro plotHaloMassesVirTemp, res=res, run=run, halos=halos

  units = getUnits()
  sP = simParams(res=res,run=run)

  ;redshifts = [4.0,3.0,2.0,1.0,0.0]
  redshifts = [3.0,0.0]
  
  bin = 0.1
  
  massColors = ['forest green','crimson','steel blue','saddle brown','black']
  tempColors = ['forest green','crimson','steel blue','saddle brown','black']
  maxTempColors = ['orange','purple','cyan','saddle brown','purple']
  
  ; config
  plotPath   = '/n/home07/dnelson/coldflows/'
  
  legSp = 0.02
  
  ; plot names
  plotName = sP.plotPath+'halo_masses_tvir_'+str(res)+'.eps'
  
  if (str(res[0]) eq 'two') then res = [256,128]  
  if (str(res[0]) eq 'all') then res = [512,256,128]  
  
  ; masses plot
  start_PS, plotName, xs=6, ys=8
  !p.multi = [0,1,2]
  
  xrange = [7.5,13.0]
  yrange = [1,1e5]
  
  title = ""
  if (n_elements(redshifts) eq 1) then title = "z = "+str(redshifts[0])
  
  fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/ylog,title=title,$
           xtitle="log ( total halo mass )",ytitle="Count",xstyle=1,ystyle=1
  
  for m=0,n_elements(redshifts)-1 do begin
    
    targetSnap = redshiftToSnapnum(redshifts[m])
    
    ; loop over each resolution
    for j=0,n_elements(res)-1 do begin    
    
      gadgetPath = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res[j])+'_20Mpc/Gadget/output/'  
  
      sgTarget = loadSubhaloGroups(sP.simPath,targetSnap,/skipIDs)
      valSGids = getPrimarySubhaloList(sgTarget,halos=halos)

      ;haloMasses = sgTarget.subgroupMass[valSGids] ;a
      ;haloMasses = sgTarget.subgroupMass ;b
      haloMasses = sgTarget.groupMass ;c
      
      haloMasses = alog10( haloMasses * (units.UnitMass_in_g / units.Msun_in_g))
    
      sgTarget = !NULL
      valSGids = !NULL
      
       ; overplot successive resolutions
      overplot = 1 ;1,1,1
      line     = j ;0,1,2
      peak     = 0
      fraction = 0
      thick    = !p.thick + 1 - (j gt 0) ;3,2,2

      plothist, haloMasses,bin=bin,/ylog,peak=peak,overplot=overplot,line=line,$
               fraction=fraction,thick=thick,color=fsc_color(massColors[m])
      
    endfor ;j
  endfor ;m
  
  ; legend
  labels = []
  for i=0,n_elements(redshifts)-1 do labels = [labels,"z="+string(redshifts[i],format='(f3.1)')]
  legend,labels,textcolors=massColors[0:n_elements(redshifts)-1],box=0,/right

  ; virial temps plot
  
  xrange = [2.5,7.5]
  yrange = [1,1e5]

  fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/ylog,$
           xtitle="log ( halo T"+textoidl("_{vir}")+" ) or log ( gas T"+textoidl("_{max}")+" )",$
           ytitle="Count",xstyle=1,ystyle=1
           
  for m=0,n_elements(redshifts)-1 do begin
    
    targetSnap = redshiftToSnapnum(redshifts[m])
    
    ; loop over each resolution
    for j=0,n_elements(res)-1 do begin    
    
      gadgetPath = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res[j])+'_20Mpc/Gadget/output/'  
  
      sgTarget = loadSubhaloGroups(sP.simPath,targetSnap,/skipIDs)
      valSGids = getPrimarySubhaloList(sgTarget,halos=halos)
      
      ;haloMasses = sgTarget.subgroupMass[valSGids] ;a
      ;haloMasses = sgTarget.subgroupMass ;b
      haloMasses = sgTarget.groupMass ;c

      sgTarget = !NULL
      valSGids = !NULL

      haloTvir = alog10(codeMassToVirTemp(haloMasses,redshifts[m])) ;log(K)

       ; overplot successive resolutions
      overplot = 1 ;1,1,1
      line     = j ;0,1,2
      peak     = 0
      thick    = !p.thick + 1 - (j gt 0) ;3,2,2

      plothist, haloTvir,bin=bin,/ylog,peak=peak,overplot=overplot,line=line,$
               thick=thick,color=fsc_color(tempColors[m])

    endfor ;j
  endfor ;m    
  
  ; legend
  labels = []
  textcl = []
  for i=0,n_elements(redshifts)-1 do labels = [labels,"z="+string(redshifts[i],format='(f3.1)')]
  for i=0,n_elements(redshifts)-1 do textcl = [textcl,maxTempColors[i]]
  legend,labels,textcolors=textcl,box=0,margin=0.25,/right
  
  !p.multi = 0
  end_PS
    
end

; plot TmaxRedshift():
;
; halos=1 : find accretion onto main halos (do not set bSH,zWidth,redshiftsCut)

pro plotTmaxRedshift, res=res, halos=halos

  units = getUnits()

  if (not keyword_set(res) or (halos ne 0 and halos ne 1)) then begin
    print,'Error: plotTmaxRedshift: Bad inputs.'
    return
  endif
  
  ; config
  plotPath   = '/n/home07/dnelson/coldflows/'

  bSH   = 0 ; include background subhalos
  
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
;
; halos=1 : find accretion onto main halos (do not set bSH,zWidth,redshiftsCut)

pro plotNormScalarProd, res=res, halos=halos

  if (not keyword_set(res) or (halos ne 0 and halos ne 1)) then begin
    print,'Error: plotNormScalarProd: Bad inputs.'
    return
  endif

  ; config
  workingPath = '/n/home07/dnelson/coldflows/thermhist.deriv/'
  plotPath    = '/n/home07/dnelson/coldflows/'  
  
  bSH = 0 ; include background subhalos
  DM  = 1 ; overplot dark matter calculation as well
  
  normRel = 0 ; 0 = all peaks normalized to 1.0 (strength of peak shown by level of non-peak)
              ; 1 = for each resolution, all histograms are normalized by the same factor such that
              ;     the highest peak reaches 1.0 (relative strengths clear but shapes not as much)
              ; 2 = for each resolution, all histograms normalized such that their center value is 1.0
              ;     (relative strength of peaks and shapes clear)
  
  redshifts = [3.0,2.0,1.0,0.0]
  zstrs = 'z='+['3','2','1','0']

  xrange = [-1.25,1.25]

  if (normRel eq 0 or normRel eq 1) then begin
    ytickv = [0.25,0.50,0.75,1.0]
    yrange = [0.0,1.05]
  endif else begin
    ytickv = [2.0,4.0,6.0,8.0,10.0]
    yrange = [0.0,10.0]
  endelse

  ; plot
  plotName = plotPath + 'normscalar_'+str(res)+'.eps'
  
  if (keyword_set(halos)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.halos.eps'
  if (keyword_set(bSH)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.bSH.eps'
  if (keyword_set(DM)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.DM.eps' 
  if (keyword_set(normRel)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.normRel='+str(normRel)+'.eps'
    
  if (str(res[0]) eq 'all') then res = [512,256,128]     
    
  start_PS, plotName
  !p.multi = [0,2,2]
  
  for m=0,n_elements(redshifts)-1 do begin
             
    for j=0,n_elements(res)-1 do begin
  
      if keyword_set(DM) then begin
        nsp = calcNormScalarProd(res=res[j],bSH=bSH,halos=halos,redshift=redshifts[m],DM=1)
        normProd_DM = flatten_list(nsp.normProd_Cold)
        nsp = !NULL
      endif
  
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
      ;if (n_elements(res) eq 1) then $
      ;  fill = 1
      
      ; histogram binsize
      bin = 0.1
      ;bin = 1.0 / sqrt(res[j])
      
      ; histogram y-normalization
      if (normRel eq 1 or normRel eq 2) then begin
        if (normRel eq 2) then begin
          ; normalize such that all centers are 1.0 and peaks are relative
          normfac = fltarr(3)
          hC = histogram(normProd_Cold,bin=bin)
          normfac[0] = 1.0 / hC[n_elements(hC)/2]
          hH = histogram(normProd_Hot,bin=bin)
          normfac[1] = 1.0 / hH[n_elements(hH)/2]
          
          if keyword_set(DM) then begin
            hDM = histogram(normProd_DM,bin=bin)
            normfac[2] = 1.0 / hDM[n_elements(hDM)/2]
          endif else begin
            hDM = 0.0
          endelse
          
          peak = 0
          ;yrange = [0.0,max([hC*normfac[0],hH*normfac[1],hDM*normfac[2]])*1.05]
          ;ytickv = findgen(5)/5 * round(yrange[1])
          ;ytickv += ytickv[1]
          ;print,yrange,ytickv
        endif ;normRel=2
        
        if (normRel eq 1) then begin
          ; normalize such that highest peak is 1.0 and others are relative
          normfac = max([histogram(normProd_Hot,bin=bin),histogram(normProd_Cold,bin=bin)])
          if keyword_set(DM) then normfac = max([normfac,histogram(normProd_DM,bin=bin)])
          normfac = replicate(1.0 / normfac,3)
          
          peak = 0
        endif ;normRel=1
      endif else begin
        normfac = [0,0,0]
        peak = 1
      endelse
      
      ; plot frames
      if (m eq 0 and j eq 0) then $
        fsc_plot,[0],[0],/nodata,ymargin=[1.0,2.0],xmargin=[7.0,0.0],$
                 xrange=xrange,yrange=yrange,/xs,/ys,xtickname=replicate(' ',10),$
                 ytickv=ytickv,yticks=n_elements(ytickv)-1
      if (m eq 1 and j eq 0) then $
        fsc_plot,[0],[0],/nodata,ymargin=[1.0,2.0],xmargin=[0.0,7.0],$
                 xrange=xrange,yrange=yrange,/xs,/ys,$
                 xtickname=replicate(' ',10),ytickname=replicate(' ',10)
      if (m eq 2 and j eq 0) then $
        fsc_plot,[0],[0],/nodata,ymargin=[4.0,-1.0],xmargin=[7.0,0.0],$
                 xrange=xrange,yrange=yrange,/xs,/ys,$
                 ytickv=[0.0,ytickv],yticks=n_elements(ytickv)
      if (m eq 3 and j eq 0) then $
        fsc_plot,[0],[0],/nodata,ymargin=[4.0,-1.0],xmargin=[0.0,7.0],$
                 xrange=xrange,yrange=yrange,/xs,/ys,ytickname=replicate(' ',10)

      ; plot histograms
      if keyword_set(DM) then $
        plothist,normProd_DM,bin=bin,normfac=normfac[2],peak=peak,$
                 fill=fill,/fline,forientation=30,fspacing=fspacing,fcolor=fsc_color('forest green'),$
                 /overplot,color=fsc_color('forest green'),line=line,thick=thick 
                 
      plothist,normProd_Cold,bin=bin,normfac=normfac[0],peak=peak,$
               fill=fill,/fline,forientation=-45,fspacing=fspacing,fcolor=fsc_color('slate blue'),$
               /overplot,color=fsc_color('slate blue'),line=line,thick=thick
      plothist,normProd_Hot,bin=bin,normfac=normfac[1],peak=peak,$
               fill=fill,/fline,forientation=45,fspacing=fspacing,fcolor=fsc_color('crimson'),$
               /overplot,color=fsc_color('crimson'),line=line,thick=thick
      
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
;
; halos=1 : find accretion onto main halos (do not set bSH,zWidth,redshiftsCut)

pro plotAccretionRate, res=res, halos=halos

  units = getUnits()

  if (not keyword_set(res) or (halos ne 0 and halos ne 1)) then begin
    print,'Error: plotAccretionRate: Bad inputs.'
    return
  endif
  
  ; config
  plotPath    = '/n/home07/dnelson/coldflows/'
  
  bSH     = 0 ; include background subhalos
  
  yScale  = 1 ; 0=normalized fraction, 1=mass/time/vol (K05)
  kerSize = 5 ; smoothing kernel for display

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

; plotTmaxMassHisto()

pro plotTmaxMassHisto, res=res, halos=halos, redshift=redshift

  units = getUnits()

  if (not keyword_set(res) or (halos ne 0 and halos ne 1) or not isnumeric(redshift)) then begin
    print,'Error: plotTmaxMassHisto: Bad inputs.'
    return
  endif
  
  ; config
  plotPath    = '/n/home07/dnelson/coldflows/'

  massBinWidth  = 0.2
  massBinsStart = [9.4,9.7,10.0,10.4,10.7,11.0,11.4,11.7,12.0]

  ; 2d histogram (DO NOT CHANGE - must redo calc)
  massMM = [9.0,13.0]
  tempMM = [4.0,7.0]
  massB  = 0.025
  tempB  = 0.025
  
  massNB = (massMM[1]-massMM[0])/massB
  tempNB = (tempMM[1]-tempMM[0])/tempB

  ; set plotName and open
  plotName = plotPath+'maxt_shmass_histo_'+str(res)+'.z='+string(redshift,format='(f3.1)')+'.eps'

  if (keyword_set(halos)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.halos.eps'  
    
  if (str(res[0]) eq 'two') then res = [256,128]
  if (str(res[0]) eq 'all') then res = [512,256,128]
    
  start_PS,plotName,xs=7,ys=6
  !p.multi = [0,3,3]

  cs = !p.charsize
  pt = !p.thick
  !p.charsize += 1.0

  ; plot config
  xrange = [3.8,7.2]
  yrange = [0.0,1.10]
      
  ; loop over each redshift (panel)
  for m=0,n_elements(massBinsStart)-1 do begin
    targetSnap    = redshiftToSnapnum(redshift)
    smoothCutSnap = !NULL
    
    ; loop over each resolution
    for j=0,n_elements(res)-1 do begin
    
      gadgetPath = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res[j])+'_20Mpc/Gadget/output/'     
      
      print,str(massBinsStart[m])+' res = ' + str(res[j]) + ' targetSnap = '+str(targetSnap)+' smoothCutSnap = ' + $
            str(smoothCutSnap)
     
      ; get 2d histogram
      cf = calcMaxTempVsMassHisto(res=res[j],halos=halos,targetSnap=targetSnap,smoothCutSnap=smoothCutSnap)  
                                     
      shmass = cf.shmass
      tm2d   = cf.tm2d
      
      cf = !NULL
      
      ; binning into mass ranges
      massBinsEnd = massBinsStart[m] + massBinWidth
      
      maxTemps = fltarr(tempNB)
      
      mass_ind1 = (massBinsStart[m] - massMM[0]) / massB
      mass_ind2 = (massBinsEnd - massMM[0]) / massB

      for i=mass_ind1,mass_ind2 do begin
        maxTemps += tm2d[i,*]
      endfor

      ; replicate for histogramming
      if (total(maxTemps) ne 0) then begin
        temps = findgen(tempNB)/tempNB * (tempMM[1]-tempMM[0]) + tempMM[0]
        histoTemps = fltarr(total(maxTemps))
        ind = 0
        for i=0,n_elements(maxTemps)-1 do begin
          len = maxTemps[i]
          if (len eq 0) then continue
          val = temps[i]
          histoTemps[ind:ind+len-1] = val
          ind += len
        endfor
      endif else begin
        histoTemps = [-1,-2]
      endelse
      
      ; overplot successive resolutions
      overplot = j gt 0 ;0,1,1
      line     = j ;0,1,2
      thick    = !p.thick + 1 - (j gt 0) ;3,2,2
      
      peak = 1
      bin  = 0.2
      
      ; 3x3 plot config
      xmargin = [[7.0,-3.0],[3.0,1.0],[-1.0,5.0],$
                 [7.0,-3.0],[3.0,1.0],[-1.0,5.0],$
                 [7.0,-3.0],[3.0,1.0],[-1.0,5.0]]
      ymargin = [[0.0,3.0],[0.0,3.0],[0.0,3.0],$
                 [2.0,0.0],[2.0,0.0],[2.0,0.0],$
                 [5.0,-2.0],[5.0,-2.0],[5.0,-2.0]]

      if (m lt 6) then begin
        xtickname=replicate(' ',10)
        xticks=0
        xtickv=0
      endif else begin
        xtickname=['4','5','6','7']
        xticks=3
        xtickv=[4.0,5.0,6.0,7.0]
      endelse
      
      if (m mod 3 eq 0) then begin ;left column
        ytickname=''
      endif else begin
        ytickname=replicate(' ',10)
      endelse

      ; plot histogram
      plothist,histoTemps,xhist,yhist,bin=bin,peak=peak,weight=weight,xtitle="",ytitle="",$
               xrange=xrange,yrange=yrange,ytickname=ytickname,/xs,/ys,$
               xmargin=xmargin[*,m],ymargin=ymargin[*,m],overplot=overplot,line=line,$
               thick=thick,xtickname=xtickname,xtickv=xtickv,xticks=xticks

      ; virial temperatures overplot
      haloTvir_start = alog10(codeMassToVirTemp(10.0^massBinsStart[m] / 1e10,redshift)) ;log(K)
      haloTvir_end   = alog10(codeMassToVirTemp(10.0^massBinsEnd / 1e10,redshift)) ;log(K)
      
      fsc_plot,[haloTvir_start,haloTvir_start],yrange,line=1,color=fsc_color('red'),/overplot
      fsc_plot,[haloTvir_end,haloTvir_end],yrange,line=1,color=fsc_color('red'),/overplot
      
      if (not overplot) then $
        fsc_text,xrange[1]*0.98,yrange[1]*0.82,$
                 string(massBinsStart[m],format='(f4.1)'),$
                 charsize=cs,alignment=1.0,color=fsc_color('orange')
  
    endfor ;j
    
  endfor ;m
  
  !p.charsize = cs
  !p.thick = pt  
  
  ; title
  if (n_elements(res) eq 1) then resStr = str(res)
  if (n_elements(res) eq 3) then resStr = 'all'
  fsc_text,0.5,0.95,"z = "+string(redshift,format='(f3.1)')+" res = "+resStr+"^3",alignment=0.5,/normal
    
  ; x/y-axis titles
  fsc_text,0.5,0.04,"log ( T"+textoidl("_{max}")+" )",alignment=0.5,/normal
  fsc_text,0.04,0.5,"Normalized Gas Accretion Rate",alignment=0.5,orientation=90,/normal

  !p.multi = 0
  end_PS
  
  ; plot 2d histogram
  for j=0,n_elements(res)-1 do begin
      
      ; load 2d histogram
      cf = calcMaxTempVsMassHisto(res=res[j],halos=halos,targetSnap=targetSnap,smoothCutSnap=smoothCutSnap)
    
      plotName = plotPath+'maxt_shmass_2dhisto_'+str(res[j])+'.z='+string(redshift,format='(f3.1)')+'.eps'   
        
    if (keyword_set(halos)) then $
      plotName = strmid(plotName,0,strlen(plotName)-4) + '.halos.eps'        
      
      ; load masses for weighting
      gadgetPath = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res[j])+'_20Mpc/Gadget/output/'
      masses = loadSnapshotSubset(gadgetPath,snapNum=targetSnap,partType='gas',field='masses')
      masses = masses[0] ;GADGET ONLY

      cf.tm2d = (units.UnitMass_in_g / units.Msun_in_g) * cf.tm2d * masses
      
      ; log(msun)
      ;cf.tm2d[where(cf.tm2d eq 0.0)] = 1.0
      ;cf.tm2d = alog10(cf.tm2d)

      start_PS, plotName
      
        loadct, 33, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
          
        tvim,cf.tm2d,pcharsize=!p.charsize,scale=1,clip=[0,100],$
             title="z = "+string(redshift,format='(f3.1)')+" res = "+str(res[j])+textoidl("^3"),$
             xtitle="log ( halo mass )",ytitle="log ( T"+textoidl("_{max}")+" )",$
             stitle="Total Mass",barwidth=0.8,lcharsize=0.0001,$;lcharsize=!p.charsize-1.5,$
             xrange=massMM,yrange=tempMM;,/rct
             
        tvirs = alog10(codeMassToVirTemp(10.0^massMM / 1e10,redshift))
        fsc_plot,massMM,tvirs,line=0,/overplot
             
       end_PS
   
   endfor ;j
   
end

; plotTransMassVsTcut(): plot the halo transition mass when f_cold=0.5 for varying values of Tcut

pro plotTransMassVsTcut, res=res, halos=halos

  critLogTemp = [5.0,5.5,6.0]
  
  if (not keyword_set(res) or (halos ne 0 and halos ne 1) or not keyword_set(critLogTemp)) then begin
    print,'Error: plotTransMassVsTcut: Bad inputs.'
    return
  endif

  ; config
  plotPath    = '/n/home07/dnelson/coldflows/'

  kSP    = 0 ; keres+ 05 spacing
  bSH    = 0 ; include background subhalos

  redshifts = [3.0,2.0,1.0,0.0]
  zstrs = 'z='+['3','2','1','0']

  colors = ['crimson','steel blue','forest green','cyan']

  ; plot config
  yrange = [10.3,13.3]
  xrange = [4.8,6.2]
  
  ; plot
  plotName = plotPath+'trans_mass_vs_tcut_'+str(res)+'.eps'
  
  if (keyword_set(halos)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.halos.eps'

  if (str(res[0]) eq 'two') then res = [256,128]
  if (str(res[0]) eq 'all') then res = [512,256,128]   

  start_PS, plotName
  
  objStr = "Subhalo"
  if (halos) then objStr = "Halo"
  
  fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
           xtitle="log ( T"+textoidl("_{cut}")+" )",$
           ytitle=textoidl("f_{cold}=1/2")+" Transition "+objStr+" Mass"

  for m=0,n_elements(redshifts)-1 do begin
    targetSnap    = redshiftToSnapNum(redshifts[m])
    smoothCutSnap = !NULL
  
    fsc_text,xrange[0]*1.05,0.5,zstrs[m],alignment=0.5
    
    for j=0,n_elements(res)-1 do begin
    
      fitSHMass = fltarr(n_elements(critLogTemp))
    
      for k=0,n_elements(critLogTemp)-1 do begin

        ; skip z=3 lowest resolutions (no fit)
        if (j gt 0 and m eq 0) then continue
        
        ; load
        cf = calcColdFracVsSubhaloMass(res=res[j],halos=halos,bSH=bSH,critLogTemp=critLogTemp[k],$
                                       targetSnap=targetSnap,smoothCutSnap=smoothCutSnap)
  
        ; calculate median line
        massStep = ceil(100.0/sqrt(n_elements(cf.coldfrac)))/10.0
        ;massStep = 0.2
        massBins = (yrange[1]-yrange[0])/massStep
        massXPts = findgen(massBins)/massBins * (yrange[1]-yrange[0]) + yrange[0] + massStep/2.0
        
        medCold = fltarr(massBins) - 1
        stdCold = fltarr(massBins) - 1
        medCold[0:floor(n_elements(medCold)/2.0)] = 1.0 ;set default value high for first half
        
        for i=0,massBins-1 do begin
          w = where(cf.shmass ge yrange[0]+i*massStep and cf.shmass lt yrange[0]+(i+1)*massStep and $
                    cf.coldfrac ne -1,count)
          if (count gt 0) then begin
            medCold[i] = median(cf.coldfrac[w])
            stdCold[i] = stddev(cf.coldfrac[w])
          endif
        endfor
        
        ; fit for medCold=0.5
        w = where(medCold lt 0.9 and medCold gt 0.1,count)
        if (count lt 2) then continue
        
        fit = linfit(massXPts[w],medCold[w])
        fitSHMass[k] = (0.5 - fit[0]) / fit[1]
        
        print,'res='+str(res[j])+' '+zstrs[m]+' cLT='+str(critLogTemp[k])+' cold frac 0.5 = '+str(fitSHMass)  
      
        if (j eq 0 and m eq 0) then begin
          critLogTempStr = string(critLogTemp[k],format='(f3.1)')
          ;fsc_text,0.5,0.95,"critLogTemp = "+critLogTempStr,alignment=0.5,/normal
        endif
      
      endfor ;k
    
     ; overplot successive resolutions
    line     = j ;0,1,2
    thick    = !p.thick + 1 - 2*(j gt 0) ;3,1,1
    color    = colors[m]
    
    fsc_plot,critLogTemp,fitSHMass,psym=-4,line=line,thick=thick,color=color,/overplot
    
    endfor ;j
  endfor ;m
  
  ; legend
  legend,zstrs,textcolors=colors[0:n_elements(redshifts)-1],box=0
  
  ; Tvir with Mhalo scaling
  xpts = [5.4,5.7]
  const = 0.85
  ypts = const * xpts^(3.0/2.0)
  fsc_plot,xpts,ypts,line=0,thick=!p.thick+1,/overplot
  fsc_text,(xpts[1]+xpts[0])/2+0.1,(ypts[1]+ypts[0])/2-0.25,textoidl("M \propto T^{3/2}"),alignment=0.5
  
  end_PS

  stop
end

; plotModeFracVsSubhaloMass(): plot the total smooth accretion rate (decomposed into hot and cold modes)
;                              as a function of total (DM+baryonic) halo (parent or subhalo) mass
; halos=1 : find accretion onto main halos (do not set bSH,zWidth,redshiftsCut)
;
; critLogTemp='vTC' : "virial temp fraction" instead of "cold mode fraction" by taking critLogTemp equal 
;                     to the virial temperature on a halo by halo basis

pro plotModeFracVsSubhaloMass, res=res, halos=halos, critLogTemp=critLogTemp

  units = getUnits()

  if (not keyword_set(res) or (halos ne 0 and halos ne 1) or not keyword_set(critLogTemp)) then begin
    print,'Error: plotModeFracVsSubhaloMass: Bad inputs.'
    return
  endif
  
  ; config
  plotPath    = '/n/home07/dnelson/coldflows/'

  kSP    = 0 ; keres+ 05 spacing
  bSH    = 0 ; include background subhalos
  
  redshifts = [3.0,2.0,1.0,0.0]
  zstrs = 'z='+['3','2','1','0']
  
  if (kSP eq 1) then begin
    ; restriction: gas not resolved at keres+ 05 redshift spacing increment before target (single)
    redshiftsCut = [3.25,2.25,1.125,0.125]
  endif else begin
    ; restriction: gas not resolved in any snapshot previous to target
    redshiftsCut = [-1,-1,-1,-1]
  endelse
  
  if (str(critLogTemp) eq 'vTC') then begin
    critLogTempStr = 'vTC'
    colors = ['forest green','purple','chartreuse']
  endif else begin
    critLogTempStr = string(critLogTemp,format='(f3.1)')
    colors = ['crimson','steel blue','sky blue']
  endelse

  ; plot config
  xrange = [9.7,13.3]
  yrange = [-0.08,1.08]
  
  plotsym,0 ;circle
  symsize = 0.3

  ; plot
  plotName = plotPath+'frac_shmass_'+str(res)+'.CT='+critLogTempStr+'.eps'
  
  if (keyword_set(halos)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.halos.eps'
  if (keyword_set(bSH)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.bSH.eps'
  if (keyword_set(kSP)) then $
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.kSP.eps'
     
  if (str(res[0]) eq 'two') then res = [256,128]
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
  
    if (str(critLogTemp) eq 'vTC') then $
      fsc_text,xrange[1]*0.95,0.5,zstrs[m],alignment=0.5
    if (str(critLogTemp) ne 'vTC') then $
      fsc_text,xrange[0]*1.05,0.5,zstrs[m],alignment=0.5
    
    for j=0,n_elements(res)-1 do begin
      cf = calcColdFracVsSubhaloMass(res=res[j],halos=halos,bSH=bSH,critLogTemp=critLogTemp,$
                                     targetSnap=targetSnap,smoothCutSnap=smoothCutSnap)

    ; plot individual systems
    if (n_elements(res) eq 1) then $
      fsc_plot,cf.shmass,cf.coldfrac,psym=8,symsize=symsize,/overplot
      
    ; calculate median line
    massStep = ceil(100.0/sqrt(n_elements(cf.coldfrac)))/10.0
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
      
    fsc_plot,massXPts[w],smooth(medCold[w],smoothSize),color=fsc_color(colors[0]),line=line,thick=thick,/overplot
    fsc_plot,massXPts[w],smooth(1.0-medCold[w],smoothSize),color=fsc_color(colors[1]),line=line,$
             thick=thick,/overplot

    ; fit for medCold=0.5
    ;w = where(massXPts ge 11.0 and massXPts le 12.5 and medCold ne 1.0 and medCold ne -1.0)
    w = where(medCold lt 0.9 and medCold gt 0.1,count)
    if (count lt 2) then continue
    
    fit = linfit(massXPts[w],medCold[w])
    fitSHMass = (0.5 - fit[0]) / fit[1]
    
    ; plot best fit coldfrac=0.5
    fsc_plot,[fitSHMass,fitSHMass],[yrange[0],yrange[0]+(yrange[1]-yrange[0])/10.0],$
              line=line,thick=thick,/overplot

    print,'res='+str(res[j])+' '+zstrs[m]+' cold frac 0.5 = '+str(fitSHMass)

    ; plot error visualization for j=0
    if (j eq 0) then begin
      w = where(stdCold gt 0,count)
      
      if (count lt smoothSize+1) then continue
      
      fsc_plot,massXPts[w],smooth(medCold[w]-stdCold[w],smoothSize),color=fsc_color('skyblue'),$
               line=line,thick=thick,/overplot
      fsc_plot,massXPts[w],smooth(medCold[w]+stdCold[w],smoothSize),color=fsc_color('skyblue'),$
               line=line,thick=thick,/overplot
    endif

    endfor ;j
  endfor ;m
  
  if (n_elements(res) eq 1) then $
    fsc_text,0.5,0.95,str(res)+"^3 critLogTemp = "+critLogTempStr,alignment=0.5,/normal
  if (n_elements(res) gt 1) then $
    fsc_text,0.5,0.95,"critLogTemp = "+critLogTempStr,alignment=0.5,/normal
    
  objStr = "Subhalo"
  if (halos) then objStr = "Halo"
  
  fsc_text,0.5,0.05,"log ( Total "+objStr+" Mass )",alignment=0.5,/normal
  
  if (str(critLogTemp) eq 'vTC') then $
    fsc_text,0.04,0.5,"Fraction Never Heated Above T"+textoidl("_{vir}"),alignment=0.5,orientation=90,/normal
  if (str(critLogTemp) ne 'vTC') then $
    fsc_text,0.04,0.5,"Cold Mode Accretion Fraction",alignment=0.5,orientation=90,/normal
  
  !p.multi = 0
  end_PS

end

; plotTmaxHisto(): plot a histogram of maximum past temperature reached by all smoothly accreted gas 
;                  particles in a series of panels at various redshifts
;
; halos=1 : find accretion onto main halos (do not set bSH,zWidth,redshiftsCut)

pro plotTmaxHisto, res=res, halos=halos

  units = getUnits()

  if (not keyword_set(res) or (halos ne 0 and halos ne 1)) then begin
    print,'Error: plotTmaxHisto: Bad inputs.'
    return
  endif
  
  ; config
  plotPath   = '/n/home07/dnelson/coldflows/'
  
  kSP    = 0 ; keres+ 05 spacing
  bSH    = 0 ; include background subhalos
  zWidth = 0 ; only if requesting a smoothCutSnap
  
  yScale = 0 ; 0=normalized fraction, 1=K05, 2=K09
  
  xColor = 'crimson'
  tempColor = 'forest green'
  
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
    
  if (str(res[0]) eq 'two') then res = [256,128]
  if (str(res[0]) eq 'all') then res = [512,256,128]
    
  start_PS,plotName,xs=7,ys=6
  !p.multi = [0,3,3]

  ; plot config
  xrange = [3.8,7.2]
  
  if (yScale eq 0) then yrange = [0.0,1.10]
  if (yScale eq 1) then yrange = [0.0,0.4]
  if (yScale eq 2) then yrange = [0.0,0.4]
      
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

      ; 3x3 plot config
      xmargin = [[7.0,-3.0],[3.0,1.0],[-1.0,5.0],$
                 [7.0,-3.0],[3.0,1.0],[-1.0,5.0],$
                 [7.0,-3.0],[3.0,1.0],[-1.0,5.0]]
      ymargin = [[0.0,3.0],[0.0,3.0],[0.0,3.0],$
                 [2.0,0.0],[2.0,0.0],[2.0,0.0],$
                 [5.0,-2.0],[5.0,-2.0],[5.0,-2.0]]

      if (m lt 6) then begin
        xtickname=replicate(' ',10)
        xticks=0
        xtickv=0
      endif else begin
        xtickname=['4','5','6','7']
        xticks=3
        xtickv=[4.0,5.0,6.0,7.0]
      endelse
      
      if (m mod 3 eq 0) then begin ;left column
        ytickname=''
      endif else begin
        ytickname=replicate(' ',10)
      endelse

      ; plot histogram
      plothist,histoTemps,xhist,yhist,bin=bin,peak=peak,weight=weight,xtitle="",ytitle="",$
               xrange=xrange,yrange=yrange,ytickname=ytickname,/xs,/ys,$
               xmargin=xmargin[*,m],ymargin=ymargin[*,m],overplot=overplot,line=line,$
               thick=thick,xtickname=xtickname,xtickv=xtickv,xticks=xticks
      
      ; scaled up histo for yScale = 1
      if (yScale ge 1 and max(yhist) lt 0.9*yrange[1]) then begin
        fac = floor(0.9*yrange[1] / max(yhist))
        fsc_text,xrange[1]*0.95,yrange[1]*0.65,"x"+str(fac),alignment=0.5,color=fsc_color(xColor),$
                 charsize=!p.charsize-1.5
        plothist,maxTemps,bin=bin,peak=peak,weight=fac*weight,overplot=1,line=line,$
                 thick=thick,color=fsc_color(xColor)
      endif
      
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
;
; halos=1: find accretion onto main halos (do not set bSH,zWidth,redshiftsCut)

pro plotRhoTemp2D, res=res, halos=halos

  units = getUnits()

  if (not keyword_set(res) or (halos ne 0 and halos ne 1)) then begin
    print,'Error: plotRhoTemp2D: Bad inputs.'
    return
  endif
  
  ; config
  plotPath  = '/n/home07/dnelson/coldflows/'

  bSH   = 0 ; include background subhalos
  tW    = 0 ; weight histogram bins by time (Gyr)
  
  nbins = res/2
  
  redshifts = [3.0,2.0,1.0,0.0]
  
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
  
  ; hot vs cold comparison plot (4x2)
  plotName = strmid(plotName,0,strlen(plotName)-4) + '.HC.eps'
  
  start_PS, plotName, xs=7.0, ys=9.0
  !p.multi = [0,4,2]  
  
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
      w = where(rt.h2rt_hot ne 0)
      h2logm_hot = fltarr(nbins,nbins)
      h2logm_hot[w] = alog10( (units.UnitMass_in_g / units.Msun_in_g) * rt.h2rt_hot[w] )
  
      w = where(rt.h2rt_cold ne 0)
      h2logm_cold = fltarr(nbins,nbins)
      h2logm_cold[w] = alog10( (units.UnitMass_in_g / units.Msun_in_g) * rt.h2rt_cold[w] )

      charsize = !p.charsize + 1
      
      xp = [0.1,0.5,0.9]
      yp = [0.1,0.3,0.5,0.7,0.9]

      ; plot
      if (m eq 0) then begin
        tvim,h2logm_hot,clip=clip,position=[xp[0],yp[3],xp[1],yp[4]],$
           stitle="",barwidth=0.0,pcharsize=charsize,/rct,$
           xrange=xrange,yrange=yrange,xtickname=replicate(' ',10),yticks=3,ytickv=[4,5,6,7]
        tvim,h2logm_cold,clip=clip,position=[xp[1],yp[3],xp[2],yp[4]],$
           stitle="",barwidth=0.0,pcharsize=charsize,/rct,$
           xrange=xrange,yrange=yrange,xtickname=replicate(' ',10),$
           ytickname=replicate(' ',10)
        fsc_text,xrange[1]*0.84,yrange[1]*0.92,"z=3",alignment=0.5
      endif
      if (m eq 1) then begin
        tvim,h2logm_hot,clip=clip,position=[xp[0],yp[2],xp[1],yp[3]],$
           stitle="",barwidth=0.0,pcharsize=charsize,/rct,$
           xrange=xrange,yrange=yrange,xtickname=replicate(' ',10),yticks=3,ytickv=[4,5,6,7]
        tvim,h2logm_cold,clip=clip,position=[xp[1],yp[2],xp[2],yp[3]],$
           stitle="",barwidth=0.0,pcharsize=charsize,/rct,$
           xrange=xrange,yrange=yrange,xtickname=replicate(' ',10),ytickname=replicate(' ',10) 
        fsc_text,xrange[1]*0.84,yrange[1]*0.92,"z=2",alignment=0.5
      endif
      if (m eq 2) then begin
        tvim,h2logm_hot,clip=clip,position=[xp[0],yp[1],xp[1],yp[2]],$
           stitle="",barwidth=0.0,pcharsize=charsize,/rct,$
           xrange=xrange,yrange=yrange,yticks=4,ytickv=[3,4,5,6,7],xtickname=replicate(' ',10)
        tvim,h2logm_cold,clip=clip,position=[xp[1],yp[1],xp[2],yp[2]],$
           stitle="",barwidth=0.0,pcharsize=charsize,/rct,$
           xrange=xrange,yrange=yrange,ytickname=replicate(' ',10),xtickname=replicate(' ',10)
        fsc_text,xrange[1]*0.84,yrange[1]*0.92,"z=1",alignment=0.5
      endif
      if (m eq 3) then begin
        tvim,h2logm_hot,clip=clip,position=[xp[0],yp[0],xp[1],yp[1]],$
           stitle="",barwidth=0.0,/rct,pcharsize=charsize,$
           xrange=xrange,yrange=yrange,yticks=3,ytickv=[3,4,5,6],xticks=5,xtickv=[-2,0,2,4,6,8]
        tvim,h2logm_cold,clip=clip,position=[xp[1],yp[0],xp[2],yp[1]],$
           stitle="",barwidth=0.0,pcharsize=charsize,/rct,$
           xrange=xrange,yrange=yrange,ytickname=replicate(' ',10),xticks=5,xtickv=[0,2,4,6,8]
        fsc_text,xrange[1]*0.84,yrange[1]*0.92,"z=0",alignment=0.5
      endif
  
    endfor ;m
    
    ; title
    fsc_text,0.5,0.97,str(res)+"^3",alignment=0.5,/normal
    fsc_text,(xp[1]+xp[0])/2.0,yp[4]+0.02,"Hot",alignment=0.5,/normal
    fsc_text,(xp[2]+xp[1])/2.0,yp[4]+0.02,"Cold",alignment=0.5,/normal
    
    ; axes
    fsc_text,0.5,0.02,"log ("+textoidl("\rho / \rho_{crit}")+")",alignment=0.5,/normal
    fsc_text,0.04,0.5,"log (T [K])",alignment=0.5,orientation=90,/normal
    fsc_colorbar,position=[xp[2]+0.02,yp[1],xp[2]+0.06,yp[3]],/vertical,/right,/reverse,$
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
