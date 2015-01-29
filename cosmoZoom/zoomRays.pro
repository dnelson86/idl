; zoomRays.pro
; 'zoom project' radial rays at different angles
; dnelson jan.2015

; helper function
pro plotColoredRays, rr=rr, co_rr=co_rr, config=config, kBounds=kBounds
  compile_opt idl2, hidden, strictarr, strictarrsubs
  if n_elements(rr) eq 0 or n_elements(co_rr) eq 0 or n_elements(config) eq 0 then message,'Error'
  
  ; config
  thick = 0.1
  
  ; select data
  xx = rr.radFacs
  
  ; colormap
  loadColorTable,config.ctName
  
  co_mm = config.yrange
  co_vv = reform( co_rr.value )
  
  if config.cbarLog eq 1 then begin
    co_mm = alog10(co_mm)
    co_vv = alog10(co_vv)
  endif
  
  co_vv = bytscl( co_vv, min=co_mm[0], max=co_mm[1] )
  
  ; subset of rays?
  if n_elements(kBounds) eq 0 then kBounds = [0,n_elements(rr.value[*,0])-1]
  
  ; plot each radial ray, corresponding to one healpix pixel
  for rayNum=kBounds[0],kBounds[1] do begin
    yy = reform( rr.value[rayNum,*] )
    
    ; select: todo
    ; sort: todo
    
    ; use temperature at each radial sample to color segment, map into [0,255] CT indices
    co = reform(co_vv[rayNum,*])
    
    for rayPt=0,n_elements(yy)-2 do $
      oplot,[xx[rayPt],xx[rayPt+1]],[yy[rayPt],yy[rayPt+1]],color=co[rayPt],line=0,thick=thick
  
  endfor ;rayNum
  
end

; plotZoomRadialRays(): visualize individual/aggregate radial rays

pro plotZoomRadialRays
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  hInds     = [0,1,2,4,5,7,8,9] ;[3,6]
  resLevels = [9,10,11]
  redshift  = 2.0
  newSaves  = 0 ; override existing shell saves
  
  ; binning config
  Nside   = [64] ;[8,16,32,64]  ; healpix resolution parameter, 2=12*4, 8=768, 16~3k, 32~12k, 64~50k
  Nrad    = [400] ;[100,200,400] ; number of radial sampling points
  cutSubS = 2 ; remove substructures? (0=no, 1=satellites, 2=all subfind other than target group)

  ; individual rays (+sorting/selections), all 1d plots, which halo/res/ns/nr?
  hIndiv  = 0   ; which halo to plot individuals from
  rIndiv  = 9   ; which resolution level to plot individuals from
  nsIndiv = 64 ; Nside indiv
  nrIndiv = 400 ; Nrad indiv

  skipRaw   = 1 ; skip plots (1) and (2) showing all rays individually (not a good idea Ns>=32)
  skipSel   = 1 ; skip select testing plots
  skip2D    = 1 ; skip making 2D panel plots
  skipStats = 1 ; skip ray statistics plots
  skipPDFs  = 0 ; skip ray PDF plots
  
  ; plot config
  config = { temp : { valName   : 'temp',$
                      removeLog : 1,$ ; remove log from healpix shell result
                      cbarLog   : 1,$ ; log minmax for colorbar scaling/label
                      cbarDivs  : 6,$ ; number of divisions of colorbar
                      ctName    : 'blue-red2',$
                      yrange    : [1e4,1e7],$ ; Kelvin
                      ytitle    : "T_{gas} [_{ }K_{ }]"},$ 
             vrad : { valName   : 'radvel',$
                      removeLog : 0,$
                      cbarLog   : 0,$
                      cbarDivs  : 9,$ ; 13 for [-400,250], 10 for [-300,200]
                      ctName    : 'blue-red2',$ ;'brewer-brownpurple'
                      yrange    : [-300,150],$ ; km/s, peculiar, including hubble flow
                      ytitle    : "v_{rad} [_{ }km/s_{ }]"} ,$ 
             entr : { valName   : 'entropy',$
                      removeLog : 0,$
                      cbarLog   : 1,$
                      cbarDivs  : 6,$
                      ctName    : 'blue-red2',$ ;'ocR/zeu'
                      yrange    : [1e7,1e10],$ ; cgs (K cm^2)
                      ytitle    : "log S_{gas} [_{ }K cm^{2 }]"} ,$
             dens : { valName   : 'density',$
                      removeLog : 0,$
                      cbarLog   : 1,$
                      cbarDivs  : 6,$
                      ctName    : 'blue-red2',$ ;'helix'
                      yrange    : [1e-6,1.0],$ ; cgs (1/cm^3)
                      ytitle    : "log n_{gas} [_{ }cm^{-3 }]"} ,$
             angm : { valName   : 'angmom',$
                      removeLog : 0,$
                      cbarLog   : 1,$
                      cbarDivs  : 5,$
                      ctName    : 'blue-red2',$
                      yrange    : [1e3,1e5],$ ; specific angular momentum (kpc km/s)
                      ytitle    : "log j_{gas} [_{ }kpc km/s_{ }]"} $
           }
  
  xrange  = [0.0,2.0] ; r/rvir
  xlog    = 0         ; plot radius in log?
  
  ; plot config/parameters (ray analysis)
  yrangePDF  = [1e-4,1e0] ; fractional/cumulative solid angle ray analysis  
  tempFacs   = [2.0,1.5,1.0,0.5,0.2,0.1]
  ytitleAng  = textoidl("Angular Direction (\theta, \phi)")
  
  radNormFac = 0.01 ; solid angle yaxis is: (Omega/4pi) (radNormFac*rvir)^(-1)
  ytitleSA = "Fractional Solid Angle (\Omega / 4\pi) ("+$
             string(radNormFac,format='(f4.2)')+"_{ }r_{vir})^{-1}"
  cIndsPDF  = fix(reverse(linspace(10,225,n_elements(tempFacs))))
  cIndsConv = fix([10,70,180,230]) ; colors for Nside convergence
  nrLines   = [2,1,0] ; linestyles for Nrad convergence
  
  ; analysis parameters
  units = getUnits(redshift=redshift)
  
  haloMass = logMsunToCodeMass(12.0) ; code units, total (dm+gas)
  haloTvir = codeMassToVirTemp(haloMass, redshift=redshift) ; Kelvin
  haloSvir = codeMassToVirEnt(haloMass, redshift=redshift) ; cgs (K cm^2)
  haloRvir = (haloMass * units.G / 100 / units.H_z^2.0)^(1.0/3.0) ; kpc (physical)
  haloVvir = sqrt(units.G * haloMass / haloRvir) ; km/s
  haloJvir = haloVvir * haloRvir ; specific halo angmom (kpc km/s)
  
  print,'haloTvir: '+str(alog10(haloTvir))+' log K'
  print,'haloSvir: '+str(alog10(haloSvir))+' log K cm^2'
  print,'haloRvir: '+str(haloRvir)+' kpc'
  print,'haloVvir: '+str(haloVvir)+' km/s'
  print,'haloJvir: '+str(alog10(haloJvir))+' log kpc km/s'
  
  pdfTvirFacs = [1.0,0.2] ; for tempFracs plots, when showing NsNr convergence or all halos
  ttLines = [3,0] ; line styles for pdfTvirFacs
  
  selections = { sel_0 : { value  : 0.5*haloTvir                       ,$ ; Kelvin (for 12bins)
                           quant  : 'TEMP'                             ,$
                           rMin   : 0.50                               ,$ ; r/rvir
                           rMax   : 2.0                                ,$ ; r/rvir
                           xtitle : "T_{thresh} [_{ }log K_{ }]"       ,$ ; for stats plot
                           xrange : [3.0,6.5]                          ,$ ; log K, for rayStats
                           xstep  : 0.01                               ,$ ; log K, for rayStats
                           xlog   : 1                                  ,$ ; for rayStats
                           vStr   : '0.5T_{vir}'                       ,$ ; for raySelect
                           sTag   : '_TmaxComp'                         } ,$
                 sel_1 : { value  : -0.25*haloVvir                     ,$ ; km/s (for 12bins)
                           quant  : 'VRAD'                             ,$
                           rMin   : 0.3                                ,$ ; r/rvir
                           rMax   : 1.5                                ,$ ; r/rvir
                           xtitle : "v_{rad,thresh} [_{ }km/s_{ }]"    ,$ ; for rayStats
                           xrange : [-400,400]                         ,$ ; km/s, for rayStats
                           xstep  : 1.0                                ,$ ; km/s, for rayStats
                           xlog   : 0                                  ,$ ; for rayStats
                           vStr   : '-0.25v_{vir}'                     ,$ ; for raySelect
                           sTag   : '_vRadMaxComp'                      } ,$
                 sel_2 : { value  : 2.0*haloSvir                       ,$ ; cgs (for 12bins)
                           quant  : 'ENTR'                             ,$
                           rMin   : 0.8                                ,$ ; r/rvir
                           rMax   : 2.0                                ,$ ; r/rvir
                           xtitle : "S_{thresh} [_{ }log K cm^{2 }]"   ,$ ; for rayStats
                           xrange : [6.0,10.0]                         ,$ ; cgs, for rayStats
                           xstep  : 0.01                               ,$ ; cgs, for rayStats
                           xlog   : 1                                  ,$ ; for rayStats
                           vStr   : '2.0S_{vir}'                       ,$ ; for raySelect
                           sTag   : '_EntrMaxComp'                      } ,$
                 sel_3 : { value  : 1.0*haloJvir                       ,$ ; cgs (for 12bins)
                           quant  : 'ANGM'                             ,$
                           rMin   : 1.1                                ,$ ; r/rvir
                           rMax   : 2.0                                ,$ ; r/rvir
                           xtitle : "J_{thresh} [_{ }log kpc km/s_{ }]",$ ; for rayStats
                           xrange : [3.4,5.2]                          ,$ ; cgs, for rayStats
                           xstep  : 0.01                               ,$ ; cgs, for rayStats
                           xlog   : 1                                  ,$ ; for rayStats
                           vStr   : '1.0J_{vir}'                       ,$ ; for raySelect
                           sTag   : '_AngMomMaxComp'                    } }
  
  sorts = { $
    sort_0:{rMin:0.9,rMax:2.0,quant:'TEMP',type:'max',   sTag:'_sortTmax'} ,$
    sort_1:{rMin:0.9,rMax:2.0,quant:'TEMP',type:'radMax',sTag:'_sortRadOfTmax'} ,$
    sort_2:{rMin:0.9,rMax:2.0,quant:'TEMP',type:'maxRTh',sTag:'_sortRmaxTh',valMin:0.5*haloTvir} $
  }
  
  ; for individual rays (always one halo, one resolution)
  nRayBins = 12 ; for visual inspection, how many distinct plots? (40&bin0 for old plot)
                ; note that NESTED pixel scheme is a quadtree with 12 tree structures, thus for 
                ; Nside=8 (768/12=64 per) we can have twelve spatially coherent bins
                ; Nside=16 (3072/12=256 per, or /4=64 per) for 12 or 64 spatially coherent bins
                ; Nside=32 (12288/12/4=256 per, or /4=64 per) for 64 or 256 spatially coherent bins
  line     = 0 ; line style
  cInd     = 1 ; color index
  thick    = 0.1 ; line thickness (for many individual ray lines)
  sK       = 3 ; smoothing kernel (only for analysis plots)
  radLines = [0.15,1.5]
  
  ; for statistics (multiple halos or resolutions)
  ; A (resolutions are different linestyles, all same color and thickness)
  ;rLines   = [1,2,0] ; line style, one per resLevel
  ;cInds    = [1,1,1] ; color index, one per resLevel
  ;rThick   = [1,1,1] ; times !p.thick

  ; B (resolutions are all solid lines, fainter and thinner for lower res)
  rLines  = [0,0,0]
  rThick  = [0.25,0.5,1.0]
  cInds   = [2,1,1]
  
  ; load
  writeNewTempFile = 0
  saveTag = 'h'+strjoin(str(hInds))+'_r'+strjoin(str(resLevels))+$
            '_ns'+strjoin(str(Nside))+'_nr'+strjoin(str(Nrad))
  saveFilename = '/scratch/dnelson/temp_plotZoomRadialRays_'+saveTag+'.sav'

  if ~file_test(saveFilename) or writeNewTempFile eq 1 then begin
  foreach hInd,hInds do begin
    hLoc = {}
        
    foreach resLevel,resLevels do begin
      ; get list of snapshots for this halo in the redshift range
      rLoc = {}
      
      if resLevel eq 11 and total(hInd eq [3,6]) gt 0 then continue
      
      sP = simParams(run='zoom_20Mpc',res=resLevel,hInd=hInd,redshift=redshift)     
      gcInd = zoomTargetHalo(sP=sP)
      
      ; calculate radial rays using nested healpix spheres, for each quantity
      hsv = {}
      
      for i=0,n_tags(config)-1 do begin
        nsStruct = {}
        
        ; for each Nside value (and each Nrad within that)
        foreach Nside_loc,Nside do begin
          nrStruct = {}
          
          foreach Nrad_loc,Nrad do begin
            
            ; r/rvir sampling points
            radFacs = linspace(0.01,2.0,Nrad_loc)
          
            hsv_loc = haloShellValue(sP=sP, partType='gas', valName=config.(i).valName, $
                                     subgroupIDs=[gcInd], Nside=Nside_loc, radFacs=radFacs, $
                                     cutSubS=cutSubS, newSaves=newSaves)
                                     
            if config.(i).removeLog eq 1 then hsv_loc.value = 10.0^hsv_loc.value ; remove log?
          
            nrStruct = mod_struct( nrStruct, 'nrad_'+str(Nrad_loc), hsv_loc )
          endforeach ;Nrad
          
          nsStruct = mod_struct( nsStruct, 'nside_'+str(Nside_loc), nrStruct )
        endforeach ;Nside
        
        hsv = mod_struct( hsv, (tag_names(config))[i], nsStruct )
        
      endfor ;n_tags(config)
      
      ; ray analysis (temporary temp)
      tempFracs = {}
      
      foreach Nside_loc,Nside,i do begin
        nrStruct = {}
        
        foreach Nrad_loc,Nrad,j do begin
        
          ; select, allocate arrays
          hsv_loc = hsv.temp.(i).(j)
        
          rfrac = { countZero  : lonarr(n_elements(tempFacs))       ,$
                    histoOuter : lonarr(hsv_loc.nRadFacs,n_elements(tempFacs))  ,$
                    histoInner : lonarr(hsv_loc.nRadFacs,n_elements(tempFacs))  ,$
                    histoMax   : lonarr(hsv_loc.nRadFacs) }
        
          ; analysis: largest or smallest radii where each ray surpasses Tvir/2
          for m=0,n_elements(tempFacs)-1 do begin
            for k=0,hsv_loc.nPx-1 do begin
              yy = reform( hsv_loc.value[k,*] )
              
              w = where(yy gt tempFacs[m]*haloTvir,count)
              if count eq 0 then begin
                rfrac.countZero[m] += 1
              endif else begin
                rfrac.histoOuter[max(w),m] += 1
                rfrac.histoInner[min(w),m] += 1
              endelse
            endfor ;k
          endfor ;m
            
          ; analysis: radius where each ray reaches its maximum temp
          for k=0,hsv_loc.nPx-1 do begin
            yy = reform( hsv_loc.value[k,*] )
            
            w = where(yy eq max(yy),count)
            rfrac.histoMax[w[0]] += 1
          endfor
        
          nrStruct = mod_struct( nrStruct, 'nrad_'+str(Nrad_loc), rfrac )
        endforeach ;j
      
        tempFracs = mod_struct( tempFracs, 'nside_'+str(Nside_loc), nrStruct )
      endforeach ;i
                
      ; combine all data structures for this halo at this resolution            
      rLoc = mod_struct( rLoc, 'sP', sP )
      rLoc = mod_struct( rLoc, 'hsv', hsv )
      rLoc = mod_struct( rLoc, 'tempFracs', tempFracs )
    
      ; add this resolution level to the halo struct
      hLoc = mod_struct( hLoc, 'L'+str(resLevel), rLoc )     
     
    endforeach ;resLevels
    
    ; add this halo to the save struct
    rr = mod_struct( rr, 'h'+str(hInd), hLoc )
    
  endforeach ;hInds
    save,rr,filename=saveFilename
    print,'Saved: ['+saveFilename+']'
  endif else begin
    restore,saveFilename,/verbose
  endelse

  ; plot setup
  plotStr = 'h'
  hNames  = 'h' + str(hInds)
  hColors = []
  rNames  = 'L' + str(resLevels)
  rLines  = rLines[0:n_elements(resLevels)-1]

  foreach hInd,hInds do plotStr += str(hInd)
  foreach hInd,hInds,i do hColors = [hColors,rr.(i).(0).sP.colors[cInds[-1]]]
  plotStr += '_L'
  foreach resLevel,resLevels do plotStr += str(resLevel)
  plotStr += '_cutSubS-' + str(cutSubS)
  plotStr += '_ns' + str(Nside[-1]) + '_nr' + str(Nrad[-1]) ; maximal Nside,Nrad

  plotStrIndiv = 'h'+str(hIndiv)+'_L'+str(rIndiv)+$
                 '_cutSbS-'+str(cutSubS)+'_ns'+str(nsIndiv)+'_nr'+str(nrIndiv)
  plotStrIndivMax = 'h'+str(hIndiv)+'_L'+str(rIndiv)+$
                    '_cutSbS-'+str(cutSubS)+'_ns'+str(Nside[-1])+'_nr'+str(Nrad[-1])
  plotStrIndiv2 = 'h'+str(hIndiv)+'_L'+str(rIndiv)+'_cutSbS-'+str(cutSubS)
  
  xminor = 2
  if xlog eq 1 then begin
    xminor = 0
    xrange[0] = 0.1
    plotStr += '_xlog'
    plotStrIndiv += '_xlog'
  endif
    
  ; select individual halo                 
  i = where( tag_names(rr) eq 'H'+str(hIndiv), count_i )
  if count_i eq 0 then message,'Error: Bad halo selecting indiv.'
  j = where( tag_names(rr.(i)) eq 'L'+str(rIndiv), count_j )
  if count_j eq 0 then message,'Error: Bad res selecting indiv.'
  k = where( tag_names(rr.(i).(j).hsv.(0)) eq 'NSIDE_'+str(nsIndiv), count_k )
  if count_k eq 0 then message,'Error: Bad ns selecting indiv.'
  m = where( tag_names(rr.(i).(j).hsv.(0).(k)) eq 'NRAD_'+str(nrIndiv), count_m )
  if count_m eq 0 then message,'Error: Bad nr selecting indiv.'
  print,'Individual selection: [h '+str(i)+' L '+str(j)+' ns '+str(k)+' nr '+str(m)+']'
    
  k_max = n_elements(Nside)-1
  m_max = n_elements(Nrad)-1
  nRaysMax = rr.(i).(j).hsv.(0).(k_max).(m_max).nPx
  
  ; plot (1) - individual rays, one plot per quantity which is plotting on the y-axis, vs radius on the x
  ;   for each plot, four panels (4 different colorings)
  nRaysTot = rr.(i).(j).hsv.(0).(k).(m).nPx
  nRaysPerBin = ceil(nRaysTot / float(nRayBins))
  
  print,'nRaysTot: ',nRaysTot,' nRaysPerBin: ',nRaysPerBin
  
  ; sub-divide all rays into some bins, for visual inspection
  if skipRaw ne 1 then begin
  if 0 then begin
  for rayBinNum=0,nRayBins do begin
    kBounds = [rayBinNum * nRaysPerBin, (rayBinNum+1) * nRaysPerBin < nRaysTot-1]
    saveTag = 'bin-'+str(rayBinNum)+'of'+str(nRayBins)
    
    ; on last loop, plot all at once
    if rayBinNum eq nRayBins then begin
      kBounds = [0,nRaysTot-1]
      saveTag = 'bin-all'
    endif
    print,'['+str(rayBinNum)+'] '+str(kBounds[0])+' - '+str(kBounds[1])+' ('+saveTag+')'
    
    saveTag += '_' + plotStrIndiv
      
    foreach yName,tag_names(config) do begin
      
    start_PS, rr.(i).(j).sP.plotPath+'radialRays_'+yName+'_'+saveTag+'.eps', xs=8.0*1.4, ys=5.0*1.4
     
      count = (where( tag_names(config) eq yName ))[0]
      coNames = ['VRAD','ENTR','TEMP','DENS','ANGM']
      
      pos = plot_pos(row=2,col=2,/gap)
      cbarSize = 0.04
      offset = [-0.02,-0.03,-0.0-cbarSize,0]
      
      foreach coName,removeIntersectionFromB(yName,coNames),pInd do begin
      
        if config.(count).cbarLog eq 0 then yminor = 2
        if config.(count).cbarLog eq 1 then yminor = 0
      
        cgPlot,[0],[0],/nodata,xrange=xrange,yrange=config.(count).yrange,/xs,/ys,$
          xlog=xlog,xminor=xminor,ylog=config.(count).cbarLog,yminor=yminor,$
          xtitle=textoidl("r / r_{vir}"),ytitle=textoidl(config.(count).ytitle),$
          /noerase,pos=pos[pInd]+offset

        foreach radLine,radLines do $
          cgPlot,[radLine,radLine],config.(count).yrange+[0.05,-0.05],line=1,color='light gray',/overplot
          
        plotColoredRays, rr=rr.(i).(j).hsv.(count).(k).(m), co=rr.(i).(j).hsv.(co_ind).(k).(m), $
          config=config.(co_ind), kBounds=kBounds
                    
        loadColorTable,config.(co_ind).ctName
        
        posCbar = pos[pInd]+offset
        posCbar[2] += cbarSize*0.8 ; move x of upperright right
        posCbar[0] = posCbar[2] - cbarSize*0.6 ; set width of colorbar
        cgColorbar,pos=posCbar,divisions=config.(co_ind).cbarDivs,/vertical,/right,range=co_mm,$
          title=textoidl(config.(co_ind).ytitle),charsize=!p.charsize*0.7
        
      endforeach ; coNames
        
    end_PS
    
    endforeach ; tag_names(config),yName
  
  endfor ; rayBins
  endif ;0
  
  ; plot (2) - all 12 ray bins for one coloring (one plot per coloring/yaxis combination, 4x3 panels)
  foreach yName,tag_names(config) do begin
  
    count = (where( tag_names(config) eq yName ))[0]
    coNames = ['VRAD','ENTR','TEMP','DENS','ANGM']
    
    pos = plot_pos(row=3,col=4)
    offset = [0.005,0.01,0,0]
    cbarSize = 0.03
    
    foreach coName,removeIntersectionFromB(yName,coNames),pInd do begin
    
      saveTag = str(nRayBins)+'bins_'+yName+'_'+coName+'_'+plotStrIndiv
      start_PS, rr.(i).(j).sP.plotPath+'radialRays_'+saveTag+'.eps', xs=8.0*1.4, ys=5.0*1.4
    
      if config.(count).cbarLog eq 0 then yminor = 2
      if config.(count).cbarLog eq 1 then yminor = 0
    
      co_ind = where( tag_names(config) eq coName )
      loadColorTable,config.(co_ind).ctName
  
      for rayBinNum=0,nRayBins-1 do begin
        kBounds = [rayBinNum * nRaysPerBin, (rayBinNum+1) * nRaysPerBin < nRaysTot-1]
      
        ytitle = ""
        xtitle = ""
        xtickn = ['0','0.5','1','1.5',' ']
        ytickn = ""
        if rayBinNum mod 4 eq 0 then ytitle = textoidl(config.(count).ytitle)
        if rayBinNum mod 4 ne 0 then ytickn = replicate(" ",10)
        if rayBinNum ge 8 then xtitle = textoidl("r / r_{vir}")
        if rayBinNum lt 8 then xtickn = replicate(" ",10)
        if rayBinNum eq nRayBins-1 then xtickn = [xtickn[0:-2],'2']
      
        cgPlot,[0],[0],/nodata,xrange=xrange,yrange=config.(count).yrange,/xs,/ys,$
          xlog=xlog,xminor=xminor,ylog=config.(count).cbarLog,yminor=yminor,$
          xtitle=xtitle,ytitle=ytitle,xtickname=xtickn,ytickname=ytickn,$
          /noerase,pos=pos[rayBinNum]+offset
          
        plotColoredRays, rr=rr.(i).(j).hsv.(count).(k).(m), co_rr=rr.(i).(j).hsv.(co_ind).(k).(m), $
          config=config.(co_ind), kBounds=kBounds
          
      endfor ; rayBins
      
      posCbar = [0.90,0.08,0.90+cbarSize,0.98]+offset
      cgColorbar,pos=posCbar,divisions=config.(co_ind).cbarDivs,/vertical,/right,range=co_mm,$
          title=textoidl(config.(co_ind).ytitle),charsize=!p.charsize*1.0
      
      end_PS
      
    endforeach ; coNames
  endforeach ; tag_names(config),yName
   
  ; plot (2b) - 12 ray bins, binary selections
  for method=0,1 do begin ; method=0,n_tags(selections)-1
    
    ; 0: (select=0 shows Tmax<fac*Tvir rays, select=1 shows Tmax>fac*Tvir rays)
    ; 1: (select=0 shows always fast infall, select=1 shows inflow which slows above X km/s)
    select_ind = (where( tag_names(config) eq selections.(method).quant ))[0]
      
    select_rVals = rr.(i).(j).hsv.(select_ind).(k).(m).radFacs
    select_rInds = where(select_rVals ge selections.(method).rMin and $
                         select_rVals le selections.(method).rMax)
  
    ; loop twice, make series of plots for passing/failing the selection criterion
    for select=0,1 do begin  
    
    foreach yName,tag_names(config) do begin
    
      count = (where( tag_names(config) eq yName ))[0]
      coNames = ['VRAD','ENTR','TEMP','DENS','ANGM']
      
      pos = plot_pos(row=3,col=4)
      offset = [0.005,0.01,0,0]
      cbarSize = 0.03
      
      foreach coName,removeIntersectionFromB(yName,coNames),pInd do begin
      
        saveTag = str(nRayBins)+'bins_'+yName+'_'+coName+'_'+plotStrIndiv+selections.(method).sTag
        start_PS, rr.(i).(j).sP.plotPath+'radialRays_'+saveTag+'_sel'+str(select)+'.eps', $
          xs=8.0*1.4, ys=5.0*1.4
      
        if config.(count).cbarLog eq 0 then yminor = 2
        if config.(count).cbarLog eq 1 then yminor = 0
      
        co_ind = where( tag_names(config) eq coName )
        loadColorTable,config.(co_ind).ctName
    
        for rayBinNum=0,nRayBins-1 do begin
          kBounds = [rayBinNum * nRaysPerBin, (rayBinNum+1) * nRaysPerBin < nRaysTot-1]
        
          ytitle = ""
          xtitle = ""
          xtickn = ['0','0.5','1','1.5',' ']
          ytickn = ""
          if rayBinNum mod 4 eq 0 then ytitle = textoidl(config.(count).ytitle)
          if rayBinNum mod 4 ne 0 then ytickn = replicate(" ",10)
          if rayBinNum ge 8 then xtitle = textoidl("r / r_{vir}")
          if rayBinNum lt 8 then xtickn = replicate(" ",10)
          if rayBinNum eq nRayBins-1 then xtickn = [xtickn[0:-2],'2']
        
          cgPlot,[0],[0],/nodata,xrange=xrange,yrange=config.(count).yrange,/xs,/ys,$
            xlog=xlog,xminor=xminor,ylog=config.(count).cbarLog,yminor=yminor,$
            xtitle=xtitle,ytitle=ytitle,xtickname=xtickn,ytickname=ytickn,$
            /noerase,pos=pos[rayBinNum]+offset
            
          ; select data
          xx = rr.(i).(j).hsv.(count).(k).(m).radFacs
          
          ; colormap
          co_mm = config.(co_ind).yrange
          co_vv = reform( rr.(i).(j).hsv.(co_ind).(k).(m).value )
          
          if config.(co_ind).cbarLog eq 1 then begin
            co_mm = alog10(co_mm)
            co_vv = alog10(co_vv)
          endif
          
          co_vv = bytscl( co_vv, min=co_mm[0], max=co_mm[1] )
          
          ; plot each radial ray, corresponding to one healpix pixel
          for rayNum=kBounds[0],kBounds[1] do begin
            yy = reform( rr.(i).(j).hsv.(count).(k).(m).value[rayNum,*] )
            
            ; ray selection
            yy_select = reform( rr.(i).(j).hsv.(select_ind).(k).(m).value[rayNum,select_rInds] )
            
            if select eq 0 then if max(yy_select) ge selections.(method).value then continue
            if select eq 1 then if max(yy_select) lt selections.(method).value then continue
            
            ; use temperature at each radial sample to color segment, map into [0,255] CT indices
            co = reform(co_vv[rayNum,*])
            
            for rayPt=0,n_elements(yy)-2 do $
              oplot,[xx[rayPt],xx[rayPt+1]],[yy[rayPt],yy[rayPt+1]],color=co[rayPt],line=0,thick=thick
          
          endfor ;rayNum
        
          if method eq 0 and yName eq 'TEMP' then $
            cgPlot,minmax(xx[select_rInds]),replicate(selections.(method).value,2),$
              line=0,thick=!p.thick,color='black',/overplot
        
          if method eq 1 and yName eq 'VRAD' then $
            cgPlot,minmax(xx[select_rInds]),replicate(selections.(method).value,2),$
              line=0,thick=!p.thick,color='black',/overplot
        
        endfor ; rayBins
              
        posCbar = [0.90,0.08,0.90+cbarSize,0.98]+offset
        cgColorbar,pos=posCbar,divisions=config.(co_ind).cbarDivs,/vertical,/right,range=co_mm,$
            title=textoidl(config.(co_ind).ytitle),charsize=!p.charsize*1.0
        
        end_PS
        
      endforeach ; coNames
    endforeach ; tag_names(config),yName
    
    endfor ; select
  endfor ; method
  
  ; plot (2c) - 12 ray bins, sorted
  for method=0,2 do begin

    ; data selection for sort
    select_ind = (where( tag_names(config) eq sorts.(method).quant ))[0]
      
    select_rVals = rr.(i).(j).hsv.(select_ind).(k).(m).radFacs
    select_rInds = where(select_rVals ge sorts.(method).rMin and $
                         select_rVals le sorts.(method).rMax)
    
    select_yy = reform( rr.(i).(j).hsv.(select_ind).(k).(m).value[*,select_rInds] )
    
    ; (i) sort rays based on value maximum at r>selectRmin
    select_min = min( select_yy, select_inds_min, dimension=2 )
    select_max = max( select_yy, select_inds_max, dimension=2 )
    
    if sorts.(method).type eq 'max' then begin
      select_sort = sort(select_max)
    endif
    
    ; (ii) sort rays based on the radius where value maximum (r>selectRmin) is reached
    if sorts.(method).type eq 'radMax' then begin
      select_inds_max_2d = array_indices( select_yy, select_inds_max )
      select_inds_max_2d = select_rInds[ reform(select_inds_max_2d[1,*]) ]
      
      select_sort = sort(select_inds_max_2d)
    endif
    
    ; (iii) sort rays based on the maximum radius where value>fac*quant
    if sorts.(method).type eq 'maxRTh' then begin
      select_rMax_above = lonarr(nRaysTot)
      
      for rayNum=0,nRaysTot-1 do begin
        locMaxRad = where( select_yy[rayNum,*] ge sorts.(method).valMin, count )
        if count gt 0 then select_rMax_above[rayNum] = max( locMaxRad )
      endfor
    
      select_sort = sort(select_rMax_above)
    endif
    
    ; loop for plots
    foreach yName,tag_names(config) do begin
    
      count = (where( tag_names(config) eq yName ))[0]
      coNames = ['VRAD','ENTR','TEMP','DENS','ANGM']
      
      pos = plot_pos(row=3,col=4)
      offset = [0.005,0.01,0,0]
      cbarSize = 0.03
      
      foreach coName,removeIntersectionFromB(yName,coNames),pInd do begin
      
        saveTag = str(nRayBins)+'bins_'+yName+'_'+coName+'_'+plotStrIndiv
        start_PS, rr.(i).(j).sP.plotPath+'radialRays_'+saveTag+sorts.(method).sTag+'.eps', $
          xs=8.0*1.4, ys=5.0*1.4
      
        if config.(count).cbarLog eq 0 then yminor = 2
        if config.(count).cbarLog eq 1 then yminor = 0
      
        co_ind = where( tag_names(config) eq coName )
        loadColorTable,config.(co_ind).ctName
    
        for rayBinNum=0,nRayBins-1 do begin
          kBounds = [rayBinNum * nRaysPerBin, (rayBinNum+1) * nRaysPerBin < nRaysTot-1]
        
          ytitle = ""
          xtitle = ""
          xtickn = ['0','0.5','1','1.5',' ']
          ytickn = ""
          if rayBinNum mod 4 eq 0 then ytitle = textoidl(config.(count).ytitle)
          if rayBinNum mod 4 ne 0 then ytickn = replicate(" ",10)
          if rayBinNum ge 8 then xtitle = textoidl("r / r_{vir}")
          if rayBinNum lt 8 then xtickn = replicate(" ",10)
          if rayBinNum eq nRayBins-1 then xtickn = [xtickn[0:-2],'2']
        
          cgPlot,[0],[0],/nodata,xrange=xrange,yrange=config.(count).yrange,/xs,/ys,$
            xlog=xlog,xminor=xminor,ylog=config.(count).cbarLog,yminor=yminor,$
            xtitle=xtitle,ytitle=ytitle,xtickname=xtickn,ytickname=ytickn,$
            /noerase,pos=pos[rayBinNum]+offset
            
          ; select data
          xx = rr.(i).(j).hsv.(count).(k).(m).radFacs
          
          ; colormap
          co_mm = config.(co_ind).yrange
          co_vv = reform( rr.(i).(j).hsv.(co_ind).(k).(m).value )
          
          if config.(co_ind).cbarLog eq 1 then begin
            co_mm = alog10(co_mm)
            co_vv = alog10(co_vv)
          endif
          
          co_vv = bytscl( co_vv, min=co_mm[0], max=co_mm[1] )
          
          ; plot each radial ray, corresponding to one healpix pixel
          for rayNum=kBounds[0],kBounds[1] do begin
            yy = reform( rr.(i).(j).hsv.(count).(k).(m).value[select_sort[rayNum],*] )
            
            ; use temperature at each radial sample to color segment, map into [0,255] CT indices
            co = reform(co_vv[select_sort[rayNum],*])
            
            for rayPt=0,n_elements(yy)-2 do $
              oplot,[xx[rayPt],xx[rayPt+1]],[yy[rayPt],yy[rayPt+1]],color=co[rayPt],line=0,thick=thick
          
          endfor ;rayNum
        
          ; legends
          if sorts.(method).type eq 'max' then begin
            select_binMin = alog10(min( select_min[select_sort[kBounds[0]:kBounds[1]]] ))
            select_binMax = alog10(max( select_max[select_sort[kBounds[0]:kBounds[1]]] ))
            if sorts.(method).quant eq 'TEMP' then qL = 'T'
            legend,textoidl(qL+'_{max}(r>r_{min}) > '+string(select_binMax,format='(f4.2)')),/top,/right
          endif
          
          if sorts.(method).type eq 'radMax' then begin
            select_binRmin = xx[min( select_inds_max_2d[select_sort[kBounds[0]:kBounds[1]]] )]
            select_binRmax = xx[max( select_inds_max_2d[select_sort[kBounds[0]:kBounds[1]]] )]
            legend,textoidl(string(select_binRmin,format='(f4.2)')+' < r < '+$
                            string(select_binRmax,format='(f4.2)')),/top,/right
          endif
          
          if sorts.(method).type eq 'maxRTh' then begin
            select_binRmin = xx[min( select_Rmax_above[select_sort[kBounds[0]:kBounds[1]]] )]
            select_binRmax = xx[max( select_Rmax_above[select_sort[kBounds[0]:kBounds[1]]] )]
            legend,textoidl(string(select_binRmin,format='(f4.2)')+' < r < '+$
                            string(select_binRmax,format='(f4.2)')),/top,/right
          endif
        
        endfor ; rayBins
              
        posCbar = [0.90,0.08,0.90+cbarSize,0.98]+offset
        cgColorbar,pos=posCbar,divisions=config.(co_ind).cbarDivs,/vertical,/right,range=co_mm,$
            title=textoidl(config.(co_ind).ytitle),charsize=!p.charsize*1.0
        
        end_PS
        
      endforeach ; coNames
    endforeach ; tag_names(config),yName
  endfor ;method
  endif ;skipRaw
  
  ; plot (3a) - combined selections (yes/no side by side), all rays raster, mean on top
  if skipSel eq 0 then begin
  methods = [1,2] ; first splits by column, second splits by row
  yName   = 'TEMP' ; what value to put on the y-axis
  coName  = 'ENTR' ; what value to color by
  xysize  = [8.0,5.0]
  nRandom = 20
  
  start_PS, rr.(i).(j).sP.plotPath + 'raySelect_2x2='+str(methods[0])+str(methods[1])+'_' + $
    plotStrIndiv + '_'+yName+'_'+coName+'.eps', xs=xysize[0]*1.4, ys=xysize[1]*1.4
  
  pos = plot_pos(row=2,col=2,/gap)
  offset = [-0.02,-0.03,-0.06,0]
  
  count = (where( tag_names(config) eq yName ))[0]
  ;coNames = ['VRAD','ENTR','TEMP','DENS','ANGM']
  
  ;foreach coName,removeIntersectionFromB(yName,coNames),pInd do begin
  
  if config.(count).cbarLog eq 0 then yminor = 2
  if config.(count).cbarLog eq 1 then yminor = 0

  co_ind = where( tag_names(config) eq coName )
  loadColorTable,config.(co_ind).ctName
  
  foreach method,methods,method_count do begin
    
    ; selection A
    select_ind = (where( tag_names(config) eq selections.(method).quant ))[0]
      
    select_rVals = rr.(i).(j).hsv.(select_ind).(k).(m).radFacs
    select_rInds = where(select_rVals ge selections.(method).rMin and $
                         select_rVals le selections.(method).rMax)
    ; selection B        
    method2 = methods[ where(methods ne method) ]
    select_ind2 = (where( tag_names(config) eq selections.(method2).quant ))[0]
      
    select_rVals2 = rr.(i).(j).hsv.(select_ind2).(k).(m).radFacs
    select_rInds2 = where(select_rVals ge selections.(method2).rMin and $
                          select_rVals le selections.(method2).rMax)             
  
    ; loop twice, make series of plots for passing/failing the selection criterion
    for select=0,1 do begin
  
      ; select data
      xx = rr.(i).(j).hsv.(count).(k).(m).radFacs
      yy = rr.(i).(j).hsv.(count).(k).(m).value
      yy_select  = reform( rr.(i).(j).hsv.(select_ind).(k).(m).value[*,select_rInds] )
      yy_select2 = reform( rr.(i).(j).hsv.(select_ind2).(k).(m).value[*,select_rInds2] )
      
      ; four permutations
      if method_count eq 0 and select eq 0 then $
        w_select = where( max(yy_select,dim=2) ge selections.(method).value and $
                          max(yy_select2,dim=2) ge selections.(method2).value )
      if method_count eq 0 and select eq 1 then $
        w_select = where( max(yy_select,dim=2) lt selections.(method).value and $
                          max(yy_select2,dim=2) lt selections.(method2).value )
      if method_count eq 1 and select eq 1 then $
        w_select = where( max(yy_select,dim=2) ge selections.(method).value and $
                          max(yy_select2,dim=2) lt selections.(method2).value )
      if method_count eq 1 and select eq 0 then $
        w_select = where( max(yy_select,dim=2) lt selections.(method).value and $
                          max(yy_select2,dim=2) ge selections.(method2).value )
  
      ; bin 2d histo for background
      loc_pos = pos[ method_count*2+select ] + offset
      loc_ar  = ( loc_pos[2]-loc_pos[0] ) / ( loc_pos[3]-loc_pos[1] ) * (xySize[0]/xySize[1])
      
      nBinsX = n_elements(rr.(i).(j).hsv.(count).(k).(m).radFacs)
      nBinsY = round(nBinsX/loc_ar)

      raysRaster = fltarr(nBinsX,nBinsY) - 1.0
      
      yMin = config.(count).yrange[0]
      yMax = config.(count).yrange[1]
      
      if config.(count).cbarLog eq 1 then begin
        yMin = alog10( yMin )
        yMax = alog10( yMax )
      endif
      
      yStep = (yMax-yMin)/nBinsY
      
      for x=0,nBinsX-1 do begin
        yLoc = reform(yy[w_select,x])
        if config.(count).cbarLog eq 1 then yLoc = alog10( yLoc )
      
        for y=0,nBinsY-1 do begin
          ; select based on x-axis and y-axis for this bin
          yMinLoc = yMin + yStep*y
          yMaxLoc = yMin + yStep*(y+1)
          
          w = where(yLoc ge yMinLoc and yLoc lt yMaxLoc,count)
          if count eq 0 then continue
          
          ; mean of color quantity in this bin
          raysRaster[x,y] = mean(rr.(i).(j).hsv.(co_ind).(k).(m).value[w_select[w],x])
        endfor
      endfor
  
      ; colormap 2d histo
      wZero = where(raysRaster eq -1.0,countZero)
      
      co_mm = config.(co_ind).yrange
      if config.(co_ind).cbarLog eq 1 then begin
        co_mm = alog10(co_mm)
        raysRaster = alog10(raysRaster)
      endif
      
      raysRaster = bytscl( raysRaster, min=co_mm[0], max=co_mm[1], top=254 ) + 1
      if countZero gt 0 then raysRaster[wZero] = 0
      
      ; draw 2d histo
      loadColorTable,config.(co_ind).ctName,/white0
      tv, raysRaster, loc_pos[0], loc_pos[1], /normal, xsize=loc_pos[2]-loc_pos[0]
      
      ; panel axes
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=config.(count).yrange,/xs,/ys,$
        xlog=xlog,ylog=config.(count).cbarLog,yminor=yminor,$
        xtitle=textoidl("r / r_{vir}"),ytitle=textoidl(config.(count).ytitle),$
        /noerase,pos=loc_pos
      
      ; draw individual rays
      co_mm = config.(co_ind).yrange
      co_vv = reform( rr.(i).(j).hsv.(co_ind).(k).(m).value )
      
      if config.(co_ind).cbarLog eq 1 then begin
        co_mm = alog10(co_mm)
        co_vv = alog10(co_vv)
      endif
      
      co_vv = bytscl( co_vv, min=co_mm[0], max=co_mm[1] )
      
      ; plot each radial ray, corresponding to one healpix pixel
      rayNums = shuffle(lindgen(n_elements(w_select)),seed=4242) ; shuffle indices of selected rays
      rayNums = rayNums[0:(nRandom<n_elements(w_select))-1] ; sub-select nRandom of them
      rayNums = w_select[rayNums] ; transform to value[] indices
      
      foreach rayNum,rayNums do begin
        yLoc = reform( rr.(i).(j).hsv.(count).(k).(m).value[rayNum,*] )
        
        ; use temperature at each radial sample to color segment, map into [0,255] CT indices
        co = reform(co_vv[rayNum,*])
        
        cgPlot,xx,yLoc,color='black',/overplot,thick=thick,line=0
        ;for rayPt=0,n_elements(yLoc)-2 do $
        ;  oplot,[xx[rayPt],xx[rayPt+1]],[yLoc[rayPt],yLoc[rayPt+1]],color=co[rayPt],line=0,thick=thick
      
      endforeach ;rayNums
      
      ; draw mean
      xx = rr.(i).(j).hsv.(count).(k).(m).radFacs
      yy = mean(rr.(i).(j).hsv.(count).(k).(m).value[w_select,*],dim=1)
      co = mean(co_vv[w_select,*],dim=1)
      
      cgPlot,xx,yy,color='black',line=0,thick=!p.thick+6.0,/overplot
      for rayPt=0,n_elements(yy)-2 do $
        oplot,[xx[rayPt],xx[rayPt+1]],[yy[rayPt],yy[rayPt+1]],color=co[rayPt],line=0,thick=!p.thick+2.0
      
      ; selection radial indicators
      if method_count eq 1 and select eq 0 then begin
        rMin = selections.(method).rMin
        rMax = selections.(method).rMax
        yMin = 10.0^4.1
        yMax = 10.0^4.2
        yMid = mean([yMin,yMax])
        
        cgPlot,[rMin,rMax],[yMid,yMid],line=0,/overplot
        cgPlot,[rMin,rMin],[yMin,yMax],line=0,/overplot
        cgPlot,[rMax,rMax],[yMin,yMax],line=0,/overplot
        
        rMin = selections.(method2).rMin
        rMax = selections.(method2).rMax
        yMin = 10.0^4.2
        yMax = 10.0^4.3
        yMid = mean([yMin,yMax])
        
        cgPlot,[rMin,rMax],[yMid,yMid],line=1,/overplot
        cgPlot,[rMin,rMin],[yMin,yMax],line=0,/overplot
        cgPlot,[rMax,rMax],[yMin,yMax],line=0,/overplot
      endif
      
      ; labels    
      if method_count eq 0 and select eq 0 then $
        selectStr = ['max('+selections.(method).quant + ') > '+textoidl(selections.(method).vStr),$
                     'max('+selections.(method2).quant + ') > '+textoidl(selections.(method2).vStr)]
      if method_count eq 0 and select eq 1 then $
        selectStr = ['max('+selections.(method).quant + ') < '+textoidl(selections.(method).vStr),$
                     'max('+selections.(method2).quant + ') < '+textoidl(selections.(method2).vStr)]
      if method_count eq 1 and select eq 1 then $
        selectStr = ['max('+selections.(method2).quant + ') < '+textoidl(selections.(method2).vStr),$
                     'max('+selections.(method).quant + ') > '+textoidl(selections.(method).vStr)]
      if method_count eq 1 and select eq 0 then $
        selectStr = ['max('+selections.(method2).quant + ') > '+textoidl(selections.(method2).vStr),$
                     'max('+selections.(method).quant + ') < '+textoidl(selections.(method).vStr)]

      legend,selectStr,/top,/right
      
      ; colorbar
      cbarSize = 0.03
      posCbar = [loc_pos[2]+0.01,loc_pos[1],loc_pos[2]+0.01+cbarSize,loc_pos[3]]
      ndiv = config.(co_ind).cbarDivs
      
      cgColorbar,pos=posCbar,divisions=ndiv,range=co_mm,/vertical,/right,$
        title=textoidl(config.(co_ind).ytitle),charsize=!p.charsize*1.0
    
    endfor ;select
    
  endforeach ;method
  end_PS  
  
  ; plot (3b) - combined selections (yes/no side by side), random rays selection, mean on top
  start_PS, rr.(i).(j).sP.plotPath + 'raySelectLines_2x2='+str(methods[0])+str(methods[1])+'_' + $
    plotStrIndiv + '_'+yName+'_'+coName+'.eps', xs=xysize[0]*1.4, ys=xysize[1]*1.4
  
  pos = plot_pos(row=2,col=2,/gap)
  offset = [-0.02,-0.03,-0.06,0]
  
  count = (where( tag_names(config) eq yName ))[0]
  ;coNames = ['VRAD','ENTR','TEMP','DENS','ANGM']
  
  ;foreach coName,removeIntersectionFromB(yName,coNames),pInd do begin
  
  if config.(count).cbarLog eq 0 then yminor = 2
  if config.(count).cbarLog eq 1 then yminor = 0

  co_ind = where( tag_names(config) eq coName )
  loadColorTable,config.(co_ind).ctName
  
  foreach method,methods,method_count do begin
    
    ; selection A
    select_ind = (where( tag_names(config) eq selections.(method).quant ))[0]
      
    select_rVals = rr.(i).(j).hsv.(select_ind).(k).(m).radFacs
    select_rInds = where(select_rVals ge selections.(method).rMin and $
                         select_rVals le selections.(method).rMax)
    ; selection B        
    method2 = methods[ where(methods ne method) ]
    select_ind2 = (where( tag_names(config) eq selections.(method2).quant ))[0]
      
    select_rVals2 = rr.(i).(j).hsv.(select_ind2).(k).(m).radFacs
    select_rInds2 = where(select_rVals ge selections.(method2).rMin and $
                          select_rVals le selections.(method2).rMax)             
  
    ; loop twice, make series of plots for passing/failing the selection criterion
    for select=0,1 do begin
  
      ; select data
      yy_select  = reform( rr.(i).(j).hsv.(select_ind).(k).(m).value[*,select_rInds] )
      yy_select2 = reform( rr.(i).(j).hsv.(select_ind2).(k).(m).value[*,select_rInds2] )
      
      ; four permutations
      if method_count eq 0 and select eq 0 then $
        w_select = where( max(yy_select,dim=2) ge selections.(method).value and $
                          max(yy_select2,dim=2) ge selections.(method2).value )
      if method_count eq 0 and select eq 1 then $
        w_select = where( max(yy_select,dim=2) lt selections.(method).value and $
                          max(yy_select2,dim=2) lt selections.(method2).value )
      if method_count eq 1 and select eq 1 then $
        w_select = where( max(yy_select,dim=2) ge selections.(method).value and $
                          max(yy_select2,dim=2) lt selections.(method2).value )
      if method_count eq 1 and select eq 0 then $
        w_select = where( max(yy_select,dim=2) lt selections.(method).value and $
                          max(yy_select2,dim=2) ge selections.(method2).value )
  
      ; panel axes
      loc_pos = pos[ method_count*2+select ] + offset
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=config.(count).yrange,/xs,/ys,$
        xlog=xlog,ylog=config.(count).cbarLog,yminor=yminor,$
        xtitle=textoidl("r / r_{vir}"),ytitle=textoidl(config.(count).ytitle),$
        /noerase,pos=loc_pos
        
      ; draw individual rays
      co_mm = config.(co_ind).yrange
      co_vv = reform( rr.(i).(j).hsv.(co_ind).(k).(m).value )
      
      if config.(co_ind).cbarLog eq 1 then begin
        co_mm = alog10(co_mm)
        co_vv = alog10(co_vv)
      endif
      
      co_vv = bytscl( co_vv, min=co_mm[0], max=co_mm[1] )
      
      ; plot each radial ray, corresponding to one healpix pixel
      rayNums = shuffle(lindgen(n_elements(w_select)),seed=4242) ; shuffle indices of selected rays
      rayNums = rayNums[0:(nRandom<n_elements(w_select))-1] ; sub-select nRandom of them
      rayNums = w_select[rayNums] ; transform to value[] indices
      
      foreach rayNum,rayNums do begin
        yy = reform( rr.(i).(j).hsv.(count).(k).(m).value[rayNum,*] )
        
        ; use temperature at each radial sample to color segment, map into [0,255] CT indices
        co = reform(co_vv[rayNum,*])
        
        for rayPt=0,n_elements(yy)-2 do $
          oplot,[xx[rayPt],xx[rayPt+1]],[yy[rayPt],yy[rayPt+1]],color=co[rayPt],line=0,thick=thick
      
      endforeach ;rayNums
      
      ; draw mean
      xx = rr.(i).(j).hsv.(count).(k).(m).radFacs
      yy = mean(rr.(i).(j).hsv.(count).(k).(m).value[w_select,*],dim=1)
      co = mean(co_vv[w_select,*],dim=1)
      
      cgPlot,xx,yy,color='black',line=0,thick=!p.thick+6.0,/overplot
      for rayPt=0,n_elements(yy)-2 do $
        oplot,[xx[rayPt],xx[rayPt+1]],[yy[rayPt],yy[rayPt+1]],color=co[rayPt],line=0,thick=!p.thick+2.0
      
      ; selection radial indicators
      if method_count eq 1 and select eq 0 then begin
        rMin = selections.(method).rMin
        rMax = selections.(method).rMax
        yMin = 10.0^4.1
        yMax = 10.0^4.2
        yMid = mean([yMin,yMax])
        
        cgPlot,[rMin,rMax],[yMid,yMid],line=0,/overplot
        cgPlot,[rMin,rMin],[yMin,yMax],line=0,/overplot
        cgPlot,[rMax,rMax],[yMin,yMax],line=0,/overplot
        
        rMin = selections.(method2).rMin
        rMax = selections.(method2).rMax
        yMin = 10.0^4.2
        yMax = 10.0^4.3
        yMid = mean([yMin,yMax])
        
        cgPlot,[rMin,rMax],[yMid,yMid],line=1,/overplot
        cgPlot,[rMin,rMin],[yMin,yMax],line=0,/overplot
        cgPlot,[rMax,rMax],[yMin,yMax],line=0,/overplot
      endif
      
      ; labels    
      if method_count eq 0 and select eq 0 then $
        selectStr = ['max('+selections.(method).quant + ') > '+textoidl(selections.(method).vStr),$
                     'max('+selections.(method2).quant + ') > '+textoidl(selections.(method2).vStr)]
      if method_count eq 0 and select eq 1 then $
        selectStr = ['max('+selections.(method).quant + ') < '+textoidl(selections.(method).vStr),$
                     'max('+selections.(method2).quant + ') < '+textoidl(selections.(method2).vStr)]
      if method_count eq 1 and select eq 1 then $
        selectStr = ['max('+selections.(method2).quant + ') < '+textoidl(selections.(method2).vStr),$
                     'max('+selections.(method).quant + ') > '+textoidl(selections.(method).vStr)]
      if method_count eq 1 and select eq 0 then $
        selectStr = ['max('+selections.(method2).quant + ') > '+textoidl(selections.(method2).vStr),$
                     'max('+selections.(method).quant + ') < '+textoidl(selections.(method).vStr)]

      legend,selectStr,/top,/right
      
      ; colorbar
      cbarSize = 0.03
      posCbar = [loc_pos[2]+0.01,loc_pos[1],loc_pos[2]+0.01+cbarSize,loc_pos[3]]
      ndiv = config.(co_ind).cbarDivs
      
      cgColorbar,pos=posCbar,divisions=ndiv,range=co_mm,/vertical,/right,$
        title=textoidl(config.(co_ind).ytitle),charsize=!p.charsize*1.0
    
    endfor ;select
    
  endforeach ;method
  end_PS
  endif ;skipSel
  
  ; plot (3c) - by hand selections, random rays selection, mean on top
  ddSK = round(nrIndiv/100.0*3)
  ddMinRad = 0.20 ; r/rvir

  cbarSize = 0.015
  
  topSize = 0.1
  yrangeTop = [-0.5,1.5]
  yValsTop = [0.0,0.5,1.0,1.5]
  nTopOther = 5.0 ; show how many other examples in top panel (thin lines)
  
  ddQuant = 'TEMP' ;'ENTR'
  ddLabel = "dT_{gas }/dr" ;"dS_{gas }/dr"
  select_ind_dd = (where( tag_names(config) eq ddQuant ))[0]
  
  laName = 'VRAD' ; what value to add spaced text labels for
  coName = 'ENTR' ; what value to color lines by
  
  if skipSel ne 1 then begin
  for plotVar=0,1 do begin
    if plotVar eq 0 then begin
      dimStr = '3x2'
      xs = xysize[0]*1.8*1.2
      ys = xysize[1]*1.5*1.2
      
      pos = plot_pos(row=2,col=3,/gap)
      offset  = [-0.01,-0.04,-0.07,0.02-topSize]
      
      selMax = 5
    endif
    
    if plotVar eq 1 then begin
      dimStr = '2x2'
      xs = xysize[0]*1.5
      ys = xysize[1]*1.8
      
      pos = plot_pos(row=2,col=2,/gap)
      offset  = [-0.03,-0.04,-0.05,0.01-topSize]
      
      selMax = 3
    endif
  
  start_PS, rr.(i).(j).sP.plotPath + 'raySelectHand_'+dimStr+'_' + $
    plotStrIndiv + '_'+yName+'_'+coName+'.eps', xs=xs, ys=ys
  
    count = (where( tag_names(config) eq yName ))[0]
    
    if config.(count).cbarLog eq 0 then yminor = 2
    if config.(count).cbarLog eq 1 then yminor = 0

    co_ind = where( tag_names(config) eq coName )
    la_ind = where( tag_names(config) eq laName )
    loadColorTable,config.(co_ind).ctName
  
    for select=0,selMax do begin
    
      loc_pos    = pos[select] + offset
      loc_posTop = [loc_pos[0],loc_pos[3],loc_pos[2],loc_pos[3]+topSize]
    
      ; selections
      if select eq 0 then begin
        quant = 'TEMP'        
        value = 2e5
        rMin = 0.2
        rMax = 1.4
        
        rMin2 = 0
        rMax2 = 0
        
        selectStr = 'distant shock'
        protRay = 5
        protRayB = 439 ; 85 dip example
      endif
      
      if select eq 1 then begin
        quant = 'TEMP'        
        value = 3e5
        rMin = 0.2
        rMax = 0.9
        
        value2 = 8e4
        rMin2 = 1.75
        rMax2 = 2.0
        
        selectStr = 'closer shock'
        protRay = 5 ;20
        protRayB = 626 ; dip example
      endif
      
      if select eq 2 then begin
        quant = 'TEMP'
        value = 3e5
        rMin  = 0.2
        rMax  = 0.75   
        
        selectStr = 'gradual heating'
        protRay = 10
      endif
      
      if select eq 3 then begin
        quant = 'ENTR'
        value = 8e8
        rMin  = 0.2
        rMax  = 2.0
        
        selectStr = 'low entropy'
        protRay = 0
      endif
      
      if select eq 4 then begin
        quant = 'TEMP'
        value = 1e5
        rMin  = 0.1
        rMax  = 0.75
        
        quant = 'VRAD'
        value = 0.0
        rMin  = 0.2
        rMax  = 1.0
        
        selectStr = 'big dips'
        protRay = 0
      endif
      
      if select eq 5 then begin
        selectStr = 'remaining'
        protRay = 0
      endif
    
      select_ind = (where( tag_names(config) eq quant ))[0]
      select_rVals = rr.(i).(j).hsv.(select_ind).(k).(m).radFacs
      select_rInds = where(select_rVals ge rMin and select_rVals le rMax)
      select_rInds2 = where(select_rVals ge rMin2 and select_rVals le rMax2)

      ; compute derivatives (-1 means positive deriv for temp increasing towards smaller rad)
      dd = rr.(i).(j).hsv.(select_ind_dd).(k).(m).value
      for d=0,n_elements(dd[*,0])-1 do dd[d,*] = -1.0 * deriv(reform(alog10(dd[d,*])))
      dd = dd/(2.0/nrIndiv)*0.1 ; deriv units: logK/(0.1r/rvir)
    
      ; select data
      yy_select = reform( rr.(i).(j).hsv.(select_ind).(k).(m).value[*,select_rInds] )
      yy_select2 = reform( rr.(i).(j).hsv.(select_ind).(k).(m).value[*,select_rInds2] )
      
      dd_max = max(dd[*,where(select_rVals ge ddMinRad)],dd_maxRad,dim=2)
      dd_min = min(dd[*,where(select_rVals ge ddMinRad)],dd_minRad,dim=2)
      dd_maxRad = select_rVals[reform((array_indices(dd,dd_maxRad))[1,*])] ; convert into r/rvir values
      dd_minRad = select_rVals[reform((array_indices(dd,dd_minRad))[1,*])]
          
      ; selections
      if select eq 0 then $
        w_select = where( smooth(dd_max,ddSK) ge 0.25 and $
                          dd_maxRad ge 1.5 and $
                          dd_maxRad le 2.0 and $
                          min(yy_select,dim=2) ge value, comp=w_remaining )
      if select eq 1 then $
        w_select = where( smooth(dd_max[w_remaining_prev],ddSK) ge 0.25 and $
                          dd_maxRad[w_remaining_prev] ge 0.80 and $
                          dd_maxRad[w_remaining_prev] le 1.5 and $
                          min(yy_select[w_remaining_prev,*],dim=2) ge value and $
                          max(yy_select2[w_remaining_prev,*],dim=2) le value2, comp=w_remaining )
      if select eq 2 then $
        w_select = where( min(yy_select[w_remaining_prev,*],dim=2) ge value and $
                          dd_max[w_remaining_prev] le 0.8, comp=w_remaining )
      if select eq 3 then $
        w_select = where( max(yy_select[w_remaining_prev,*],dim=2) le value, comp=w_remaining)
      if select eq 4 then $
        w_select = where( max(yy_select[w_remaining_prev,*],dim=2) ge value, comp=w_remaining )
      if select eq 5 then $
        w_select = lindgen(n_elements(w_remaining_prev)) ; remaining  
  
      ; map indices into original rays (both for selection and remaining)
      if select gt 0 then begin
        w_select = w_remaining_prev[w_select]
        if select lt 5 then w_remaining = w_remaining_prev[w_remaining]
      endif
      
      w_remaining_prev = w_remaining
    
      ; panel axes
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=config.(count).yrange,xs=5,ys=5,$
        xlog=xlog,ylog=config.(count).cbarLog,yminor=yminor,$
        xtitle="",ytitle="",ytickname=replicate(" ",10),xtickname=replicate(" ",10),$
        /noerase,pos=loc_pos
        
      ; draw individual rays
      co_mm = config.(co_ind).yrange
      co_vv = reform( rr.(i).(j).hsv.(co_ind).(k).(m).value )
      
      if config.(co_ind).cbarLog eq 1 then begin
        co_mm = alog10(co_mm)
        co_vv = alog10(co_vv)
      endif
      
      co_vv = bytscl( co_vv, min=co_mm[0], max=co_mm[1] )
      
      ; plot each radial ray, corresponding to one healpix pixel
      ;rayNums = shuffle(lindgen(n_elements(w_select)),seed=4242) ; shuffle indices of selected rays
      ;;;;rayNums = rayNums[0:(nRandom<n_elements(w_select))-1] ; sub-select nRandom of them
      ;rayNums = w_select[rayNums] ; transform to value[] indices
      rayNums = w_select
      
      foreach rayNum,rayNums do begin
        yy = reform( rr.(i).(j).hsv.(count).(k).(m).value[rayNum,*] )
        
        ; use temperature at each radial sample to color segment, map into [0,255] CT indices
        co = reform(co_vv[rayNum,*])
        
        for rayPt=0,n_elements(yy)-2 do $
          oplot,[xx[rayPt],xx[rayPt+1]],[yy[rayPt],yy[rayPt+1]],color=co[rayPt],line=0,thick=thick
      
      endforeach ;rayNums
      
      ; draw prototypical example(s) thicker
      xx = rr.(i).(j).hsv.(count).(k).(m).radFacs
      yy = rr.(i).(j).hsv.(count).(k).(m).value[rayNums[protRay],*]
      co = co_vv[rayNums[protRay],*]
      
      cgPlot,xx,yy,color='black',line=0,thick=!p.thick+6.0,/overplot
      for rayPt=0,n_elements(yy)-2 do $
        oplot,[xx[rayPt],xx[rayPt+1]],[yy[rayPt],yy[rayPt+1]],color=co[rayPt],line=0,thick=!p.thick+2.0
      
      ; draw axes and ticks on top
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=config.(count).yrange,/xs,/ys,$
        xlog=xlog,ylog=config.(count).cbarLog,yminor=yminor,$
        xtitle=textoidl("r / r_{vir}"),ytitle=textoidl(config.(count).ytitle),$
        /noerase,pos=loc_pos
      
      ; selection radial indicators
      if 0 then begin
        yMin = 10.0^4.1
        yMax = 10.0^4.2
        yMid = mean([yMin,yMax])
        
        cgPlot,[rMin,rMax],[yMid,yMid],line=1,/overplot
        cgPlot,[rMin,rMin],[yMin,yMax],line=0,/overplot
        cgPlot,[rMax,rMax],[yMin,yMax],line=0,/overplot
        
        yMin = 10.0^4.2
        yMax = 10.0^4.28
        yMid = mean([yMin,yMax])
        
        cgPlot,[rMin2,rMax2],[yMid,yMid],line=1,/overplot
        cgPlot,[rMin2,rMin2],[yMin,yMax],line=0,/overplot
        cgPlot,[rMax2,rMax2],[yMin,yMax],line=0,/overplot
      endif ;0
      
      ; labels
      laSpacing = ceil(nrIndiv/8.0)
      la_vv = reform( rr.(i).(j).hsv.(la_ind).(k).(m).value )
      
      for rayPt=laSpacing,n_elements(yy)-laSpacing/2,laSpacing do begin
        xLa = xx[rayPt]
        
        ; offset above local mean y-value
        yLa = max(yy[rayPt-5>0:rayPt+5<n_elements(yy)-1])*1.5 
        
        ; label = mean of selected rays at this rad
        vLa = mean( reform(la_vv[rayNums,rayPt]), dim=1 ) 
        ; label = prototypical ray at this rad
        ;vLa = la_vv[rayNums[protRay],rayPt] 
        ; label = local mean of protRay
        ;vLa = mean( la_vv[rayNums[protRay],rayPt-5>0:rayPt+5<n_elements(yy)-1] )
        
        sLa = str(string(round(vLa/5.0)*5.0,format='(I4)')) ; round to nearest 5 km/s
        
        cgText,xLa,yLa,sLa,color='black',alignment=0.5
      endfor
      
      if select eq 0 or select eq 1 then begin
        ; exception example
        w_ex = where( smooth(dd_max,ddSK) ge 0.1 and $
                          dd_maxRad ge 1.25 and $
                          dd_maxRad le 2.0 and $
                          min(yy_select,dim=2) lt 1e5 )

        yy = rr.(i).(j).hsv.(count).(k).(m).value[protRayB,*]
        cgPlot,xx,yy,color='black',line=0,thick=thick,/overplot
      endif
      
      ; legend
      rayFrac = n_elements(w_select)*100.0/n_elements(yy_select[*,0])
      legend,[textoidl(selectStr),string(rayFrac,format='(f4.1)')+'%'],/top,/right
      
      ; colorbar
      posCbar = [loc_pos[2]+0.01,loc_pos[1],loc_pos[2]+0.01+cbarSize,loc_pos[3]]
      ndiv = config.(co_ind).cbarDivs
      
      cgColorbar,pos=posCbar,divisions=ndiv,range=co_mm,/vertical,/right,$
        title=textoidl(config.(co_ind).ytitle),charsize=!p.charsize*1.0
        
      ; top panel
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeTop,/xs,/ys,xlog=xlog,ylog=0,$
        xtitle="",xtickname=replicate(" ",10),ytitle=textoidl(ddLabel),$
        /noerase,pos=loc_posTop,ytickv=yValsTop,yticks=n_elements(yValsTop)-1
      
      ; draw derivative of prototypical ray
      xx = rr.(i).(j).hsv.(count).(k).(m).radFacs
      yy = dd[rayNums,*]
      ;co = co_vv[rayNums[0],*]
      ww = where(xx ge ddMinRad)
      
      raySkip = round(n_elements(rayNums)/5.0) ; show 5
      
      cgPlot,xx[ww],smooth(reform(yy[protRay,ww]),ddSK),color='black',line=0,thick=!p.thick+1.0,/overplot
      for rayNum=0,n_elements(rayNums)-1,raySkip do $
        cgPlot,xx[ww],reform(yy[rayNum,ww]),color='black',line=0,thick=thick,/overplot

      ; top histogram
      ;posHisto = [loc_posTop[2],loc_posTop[1],loc_posTop[2]+0.05,loc_posTop[3]]
      
      ;dd = rr.(i).(j).hsv.VRAD.(k).(m).value
      ;for d=0,n_elements(dd[*,0])-1 do dd[d,*] = -1.0 * deriv(reform(alog10(dd[d,*])))
      ;dd = dd/(2.0/nrIndiv)*0.1 ; deriv units: logK/(0.1r/rvir)
      
      ;dd_h = dd[*,where(select_rVals ge ddMinRad)]
      ;dd_h = dd_h[w_select,*]
      ;;dd_h = dd_max[w_select]
      ;cgHistoplot,dd_h,binsize=0.1,pos=posHisto,/noerase,/rotate,$
      ;  yrange=yrangeTop,/frequency,xs=9,ys=9,$
      ;  xtickname=replicate(" ",10),ytickname=replicate(" ",10),xtitle=""
      
    endfor ; select
      
  end_PS
  
  endfor ;plotVar  
  endif ;skipSel
  
  ; print stats
  if 0 then begin
  meanStatFracs = fltarr(n_tags(rr),n_tags(rr.(0)),5)
  
  for ii=0,n_tags(rr)-1 do begin
    for jj=0,n_tags(rr.(ii))-1 do begin
      for select=0,4 do begin
      
        ; selections
        if select eq 0 then begin
          quant = 'TEMP'        
          value = 2e5
          rMin = 0.2
          rMax = 1.4
          
          rMin2 = 0
          rMax2 = 0
          
          selectStr = 'distant shock'
          protRay = 5
          protRayB = 439 ; 85 dip example
        endif
        
        if select eq 1 then begin
          quant = 'TEMP'        
          value = 3e5
          rMin = 0.2
          rMax = 0.9
          
          value2 = 8e4
          rMin2 = 1.75
          rMax2 = 2.0
          
          selectStr = 'closer shock'
          protRay = 5 ;20
          protRayB = 626 ; dip example
        endif
        
        if select eq 2 then begin
          quant = 'TEMP'
          value = 3e5
          rMin  = 0.2
          rMax  = 0.75   
          
          selectStr = 'gradual heating'
          protRay = 10
        endif
        
        if select eq 3 then begin
          quant = 'ENTR'
          value = 8e8
          rMin  = 0.2
          rMax  = 2.0
          
          selectStr = 'low entropy'
          protRay = 0
        endif
        
        if select eq 4 then begin
          quant = 'TEMP'
          value = 1e5
          rMin  = 0.1
          rMax  = 0.75
          
          quant = 'VRAD'
          value = 0.0
          rMin  = 0.2
          rMax  = 1.0
          
          selectStr = 'big dips'
          protRay = 0
        endif
        
        if select eq 5 then begin
          selectStr = 'remaining'
          protRay = 0
        endif
      
        select_ind = (where( tag_names(config) eq quant ))[0]
        select_rVals = rr.(ii).(jj).hsv.(select_ind).(k).(m).radFacs
        select_rInds = where(select_rVals ge rMin and select_rVals le rMax)
        select_rInds2 = where(select_rVals ge rMin2 and select_rVals le rMax2)

        ; compute derivatives (-1 means positive deriv for temp increasing towards smaller rad)
        dd = rr.(ii).(jj).hsv.(select_ind_dd).(k).(m).value
        for d=0,n_elements(dd[*,0])-1 do dd[d,*] = -1.0 * deriv(reform(alog10(dd[d,*])))
        dd = dd/(2.0/nrIndiv)*0.1 ; deriv units: logK/(0.1r/rvir)
      
        ; smooth derivative?
        if nrIndiv eq 200 then dd = smooth(dd,3,/nan)
        if nrIndiv eq 400 then dd = smooth(dd,5,/nan)
      
        ; select data
        yy_select = reform( rr.(ii).(jj).hsv.(select_ind).(k).(m).value[*,select_rInds] )
        yy_select2 = reform( rr.(ii).(jj).hsv.(select_ind).(k).(m).value[*,select_rInds2] )
        
        dd_max = max(dd[*,where(select_rVals ge ddMinRad)],dd_maxRad,dim=2)
        dd_min = min(dd[*,where(select_rVals ge ddMinRad)],dd_minRad,dim=2)
        dd_maxRad = select_rVals[reform((array_indices(dd,dd_maxRad))[1,*])] ; convert into r/rvir values
        dd_minRad = select_rVals[reform((array_indices(dd,dd_minRad))[1,*])]
            
        ; selections
        if select eq 0 then $
          w_select = where( smooth(dd_max,ddSK) ge 0.25 and $
                            dd_maxRad ge 1.5 and $
                            dd_maxRad le 2.0 and $
                            min(yy_select,dim=2) ge value, comp=w_remaining )
        if select eq 1 then $
          w_select = where( smooth(dd_max[w_remaining_prev],ddSK) ge 0.25 and $
                            dd_maxRad[w_remaining_prev] ge 0.80 and $
                            dd_maxRad[w_remaining_prev] le 1.5 and $
                            min(yy_select[w_remaining_prev,*],dim=2) ge value and $
                            max(yy_select2[w_remaining_prev,*],dim=2) le value2, comp=w_remaining )
        if select eq 2 then $
          w_select = where( min(yy_select[w_remaining_prev,*],dim=2) ge value and $
                            dd_max[w_remaining_prev] le 0.8, comp=w_remaining )
        if select eq 3 then $
          w_select = where( max(yy_select[w_remaining_prev,*],dim=2) le value, comp=w_remaining)
        if select eq 4 then $
          w_select = where( max(yy_select[w_remaining_prev,*],dim=2) ge value, comp=w_remaining )
        if select eq 5 then $
          w_select = lindgen(n_elements(w_remaining_prev)) ; remaining  

        ; map indices into original rays (both for selection and remaining)
        if select gt 0 then begin
          w_select = w_remaining_prev[w_select]
          if select lt 5 then w_remaining = w_remaining_prev[w_remaining]
        endif
        
        w_remaining_prev = w_remaining
        
        ; fractions
        rayFrac = n_elements(w_select)*100.0/n_elements(yy_select[*,0])
        print,'h'+str(rr.(ii).(jj).sP.hInd)+'L'+str(rr.(ii).(jj).sP.res)+': '+$
          selectStr+' '+string(rayFrac,format='(f4.1)')+'%'
          
        meanStatFracs[ii,jj,select] = rayFrac
      endfor ;select
    endfor ;jj
  endfor ;ii
  
  for select=0,4 do begin
    for jj=0,n_tags(rr.(0))-1 do begin
      xx = reform(meanStatFracs[*,jj,select])
      print,'L'+str(rr.(0).(jj).sP.res)+' mean (select='+str(select)+'):'
      print,' '+string(mean(xx),format='(f4.1)')+' pm '+string(stddev(xx),format='(f4.1)')
    endfor
  endfor
  endif ;0
  
  ; ----------------------------------------------------------------------------------
  
  ; plot (4a) - global 2D plot of all rays, single coloring per plot, one plot per value
  if skip2D ne 1 then begin
  foreach yName,tag_names(config) do begin
    
  start_PS, rr.(i).(j).sP.plotPath + 'radialRays2DAll_'+yName+'_'+plotStrIndiv+'.eps', $
    xs=4.0*1.4, ys=8.0*1.4
   
    ind = (where( tag_names(config) eq yName ))[0]
    
    pos = [0.12,0.06,0.84,0.97]
    cbarSize = 0.03
    
    ; colormap
    loadColorTable,config.(ind).ctName
    
    co_mm = config.(ind).yrange
    co_vv = reform( rr.(i).(j).hsv.(ind).(k).(m).value )
    
    if config.(ind).cbarLog eq 1 then begin
      co_mm = alog10(co_mm)
      co_vv = alog10(co_vv)
    endif
    
    co_vv = bytscl( co_vv, min=co_mm[0], max=co_mm[1] )
    
    ; plot
    co_vv = transpose(co_vv)    
    tvim,co_vv,scale=0,/c_map,pos=pos,/noframe,/noerase
    
    ytickv = lindgen(nRayBins+1)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=[0,1],/xs,/ys,$
      xlog=xlog,xminor=xminor,ylog=0,yminor=4,xtitle=textoidl("r / r_{vir}"),ytitle=ytitleAng,$
      ytickv=ytickv,yticks=n_elements(ytickv)-1,/noerase,pos=pos
    
    ; colorbar
    posCbar = [pos[2]+0.015,pos[1],pos[2]+0.015+cbarSize,pos[3]]
    cgColorbar,pos=posCbar,divisions=config.(ind).cbarDivs,/vertical,/right,range=co_mm,$
      title=textoidl(config.(ind).ytitle),charsize=!p.charsize*1.0
      
  end_PS
  
  endforeach ;tag_names(config),yName
  
  ; plot (4b) - all 2D panels on one plot (5 narrow columns)
  start_PS, rr.(i).(j).sP.plotPath + 'radialRays2DAll_'+plotStrIndiv+'.eps', xs=10.0*1.4, ys=8.0*1.4
   
    posBounds = [0.06,0.06,0.82,0.90]
    cbarSize = 0.03
    
    spacing = 0.05
    indivWidth = (posBounds[2]-posBounds[0]-spacing)/n_tags(config)
   
    foreach yName,tag_names(config),count do begin 
   
      ind = (where( tag_names(config) eq yName ))[0]
      
      pos = posBounds
      pos[0] = pos[0] + count*indivWidth + count*spacing
      pos[2] = pos[0] + indivWidth
      
      ; colormap
      loadColorTable,config.(ind).ctName ;'helix'
      
      co_mm = config.(ind).yrange
      co_vv = reform( rr.(i).(j).hsv.(ind).(k).(m).value )
      
      if config.(ind).cbarLog eq 1 then begin
        co_mm = alog10(co_mm)
        co_vv = alog10(co_vv)
      endif
      
      co_vv = bytscl( co_vv, min=co_mm[0], max=co_mm[1] )
      
      ; plot
      co_vv = transpose(co_vv)      
      tvim,co_vv,scale=0,/c_map,pos=pos,/noframe,/noerase
      
      ytickv = lindgen(nRayBins+1)
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=[0,1],/xs,/ys,$
        xlog=xlog,xminor=xminor,ylog=0,yminor=4,xtitle=textoidl("r / r_{vir}"),ytitle=ytitleAng,$
        ytickv=ytickv,yticks=n_elements(ytickv)-1,/noerase,pos=pos
      
      ; colorbar
      posCbar = [pos[0],pos[3]+0.015,pos[2],pos[3]+0.015+cbarSize]
      ndiv = config.(ind).cbarDivs
      if yName eq 'VRAD' then ndiv = 3
      if yName eq 'ENTR' or yName eq 'TEMP' then ndiv /= 2
      
      cgColorbar,pos=posCbar,divisions=ndiv,range=co_mm,/top,$
        title=textoidl(config.(ind).ytitle),charsize=!p.charsize*1.0
      
    endforeach ;tag_names(config),yName
      
  end_PS
  endif ;skip2D
  
  ; ----------------------------------------------------------------------------------
  
  ; plot (5a) - ray statistics based on selections, one quantity per plot
  if skipStats ne 1 then begin
  for method=0,n_tags(selections)-1 do begin
  
  start_PS, rr.(i).(j).sP.plotPath + 'rayStats'+selections.(method).sTag+'_' + plotStr + '.eps'
  
    cgPlot,[0],[0],/nodata,xrange=selections.(method).xrange,$
      yrange=yrangePDF*[10,1.2],/xs,/ys,xlog=xlog,/ylog,yminor=5,$
      xtitle=textoidl(selections.(method).xtitle),$
      ytitle=textoidl("Fraction of Solid Angle Covered (\Omega / 4\pi)"),/noerase

    ; include all halos
    for ii=0,n_tags(rr)-1 do begin
      jj_L11 = where( tag_names(rr.(ii)) eq 'L11', count_ind )
      if count_ind eq 0 then continue
      
      ; indices
      select_ind = (where( tag_names(config) eq selections.(method).quant ))[0]
        
      select_rVals = rr.(ii).(jj_L11).hsv.(select_ind).(k_max).(m_max).radFacs
      select_rInds = where(select_rVals ge selections.(method).rMin and $
                           select_rVals le selections.(method).rMax)
          
      ; calculate bins
      xrange_loc = selections.(method).xrange
      xBins = ceil( (xrange_loc[1]-xrange_loc[0])/ selections.(method).xstep )
      
      binVals = fltarr(xBins)
      angFrac = fltarr(xBins,6)
      
      ; bin
      yy_select = reform( rr.(ii).(jj_L11).hsv.(select_ind).(k_max).(m_max).value[*,select_rInds] )
      
      yy_max  = max(  yy_select, dimension=2 )
      yy_min  = min(  yy_select, dimension=2 )
      yy_mean = mean( yy_select, dimension=2 )
      
      if selections.(method).xlog eq 1 then begin
        yy_max  = alog10(yy_max)
        yy_min  = alog10(yy_min)
        yy_mean = alog10(yy_mean)
      endif
      
      for bin=0,xBins-1 do begin
        binVals[bin] = xrange_loc[0] + selections.(method).xstep*(bin+0.5)
        
        w = where(yy_max  ge binVals[bin], count_max_ge,  ncomp=count_max_lt)
        w = where(yy_min  ge binVals[bin], count_min_ge,  ncomp=count_min_lt)
        w = where(yy_mean ge binVals[bin], count_mean_ge, ncomp=count_mean_lt) 
        
        angFrac[bin,*] = [count_max_ge,count_max_lt,$
                          count_min_ge,count_min_lt,$
                          count_mean_ge,count_mean_lt]
      endfor
      
      ; normalize
      angFrac /= float(nRaysMax)
      
      ; plot
      co = rr.(ii).(jj_L11).sP.colors[cInds[jj_L11]]
      th = !p.thick * rThick[jj_L11]
          
      cgPlot,binVals,angFrac[*,0],line=0,thick=th,color=co,/overplot
      cgPlot,binVals,angFrac[*,1],line=1,thick=th,color=co,/overplot
      cgPlot,binVals,angFrac[*,2],line=2,thick=th,color=co,/overplot
      cgPlot,binVals,angFrac[*,3],line=3,thick=th,color=co,/overplot
      cgPlot,binVals,angFrac[*,4],line=4,thick=th,color=co,/overplot
      cgPlot,binVals,angFrac[*,5],line=5,thick=th,color=co,/overplot
    
    endfor ;ii
    
    ; legends
    loadColorTable,'bw linear'    
    
    legend,hNames,textcolor=hColors,/top,/right
    
    if selections.(method).quant eq 'TEMP' then qL = 'T'
    if selections.(method).quant eq 'ANGM' then qL = 'J'
    if selections.(method).quant eq 'VRAD' then qL = 'v_{rad}'
    if selections.(method).quant eq 'ENTR' then qL = 'S'
    
    legendStrs = ['max('+qL+') \geq','max('+qL+') <',$
                  'min('+qL+') \geq','min('+qL+') <',$
                  'mean('+qL+') \geq','mean('+qL+') <']
    legend,textoidl(legendStrs),linestyle=[0,1,2,3,4,5],/bottom,/left
    
  end_PS
  
  endfor ;method
  
  ; plot (5b) - ray statistics based on selections, 2x2 panels (four quantities)
  start_PS, rr.(i).(j).sP.plotPath + 'rayStats_2x2_' + plotStr + '.eps', xs=8*1.4, ys=5*1.4
  
  pos = plot_pos(row=2,col=2,/gap)
  offset = [-0.01,-0.03,0.0,0]
  
  for method=0,n_tags(selections)-1 do begin
  
    cgPlot,[0],[0],/nodata,xrange=selections.(method).xrange,$
      yrange=[0.0,1.05],/xs,/ys,xlog=xlog,$;/ylog,yminor=5,$ ;yrangePDF*[10,1.2]
      xtitle=textoidl(selections.(method).xtitle),pos=pos[method]+offset,$
      ytitle=textoidl("Solid Angle Covered (\Omega / 4\pi)"),/noerase

    ; include all halos
    for ii=0,n_tags(rr)-1 do begin
      jj_L11 = where( tag_names(rr.(ii)) eq 'L11', count_ind )
      if count_ind eq 0 then continue
      
      ; indices
      select_ind = (where( tag_names(config) eq selections.(method).quant ))[0]
        
      select_rVals = rr.(ii).(jj_L11).hsv.(select_ind).(k_max).(m_max).radFacs
      select_rInds = where(select_rVals ge selections.(method).rMin and $
                           select_rVals le selections.(method).rMax)
          
      ; calculate bins
      xrange_loc = selections.(method).xrange
      xBins = ceil( (xrange_loc[1]-xrange_loc[0])/ selections.(method).xstep )
      
      binVals = fltarr(xBins)
      angFrac = fltarr(xBins,6)
      if ii eq 0 then meanAngFrac = fltarr(xBins,6)
      
      ; bin
      yy_select = reform( rr.(ii).(jj_L11).hsv.(select_ind).(k_max).(m_max).value[*,select_rInds] )
      
      yy_max  = max(  yy_select, dimension=2 )
      yy_min  = min(  yy_select, dimension=2 )
      yy_mean = mean( yy_select, dimension=2 )
      
      if selections.(method).xlog eq 1 then begin
        yy_max  = alog10(yy_max)
        yy_min  = alog10(yy_min)
        yy_mean = alog10(yy_mean)
      endif
      
      for bin=0,xBins-1 do begin
        binVals[bin] = xrange_loc[0] + selections.(method).xstep*(bin+0.5)
        
        w = where(yy_max  ge binVals[bin], count_max_ge,  ncomp=count_max_lt)
        w = where(yy_min  ge binVals[bin], count_min_ge,  ncomp=count_min_lt)
        w = where(yy_mean ge binVals[bin], count_mean_ge, ncomp=count_mean_lt) 
        
        angFrac[bin,*] = [count_max_ge,count_max_lt,$
                          count_min_ge,count_min_lt,$
                          count_mean_ge,count_mean_lt]
      endfor
      
      ; mean accumulate
      meanAngFrac += angFrac

      ; normalize
      angFrac /= float(nRaysMax)
      
      ; plot
      coth_ind = 0 ;jj_L11 ;0
      co = rr.(ii).(jj_L11).sP.colors[cInds[coth_ind]]
      th = !p.thick * rThick[coth_ind]
          
      for ff=0,5 do $
        cgPlot,binVals,angFrac[*,ff],line=ff,thick=th,color=co,/overplot
    
    endfor ;ii

    ; mean across the individual halos as black lines
    meanAngFrac /= float(nRaysMax*n_tags(rr)) ; normalize

    for ff=0,5 do $
      cgPlot,binVals,meanAngFrac[*,ff],line=ff,thick=!p.thick*rThick[jj_L11]+1.0,color='black',/overplot
    
    ; legends
    loadColorTable,'bw linear'    
    
    if method eq 0 then legend,hNames,textcolor=hColors,pos=[3.3,0.8],/data
    
    legendStrs = ['\geq max','< max',$
                  '\geq min','< min',$
                  '\geq mean','< mean']
    if method eq 2 then legend,textoidl(legendStrs),linestyle=[0,1,2,3,4,5],$
      pos=[6.2,0.75],/data
    
  endfor ;method
    
  end_PS
  endif ;skipStats
      
  ; ----------------------------------------------------------------------------------
      
  ; plot (6a) - ray PDF (fraction of sphere covered) based on temp criterion (indiv halo)
  if skipPDFs eq 0 then begin
  start_PS, rr.(i).(j).sP.plotPath + 'rayRadHisto_TempFacs_' + plotStrIndivMax + '.eps'
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangePDF*[10,0.1],/xs,/ys,xlog=xlog,/ylog,yminor=5,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl(ytitleSA),/noerase

    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrangePDF,line=1,color='light gray',/overplot
      
    loadColorTable,'blue-red2'
      
    xx = rr.(i).(j).hsv.temp.(k_max).(m_max).radFacs
        
    ; normalization
    radBinSize = xrange[1]/Nrad[m_max]
    radNorm = radBinSize/radNormFac
        
    ; plot
    for tt=0,n_elements(tempFacs)-1 do begin
      ; rmax
      yy = rr.(i).(j).tempFracs.(k_max).(m_max).histoOuter[*,tt] / float(nRaysMax)
      yy /= radNorm
      oplot,xx,smooth(yy,sK,/nan),line=0,color=cIndsPDF[tt]
      
      ; rmin
      yy = rr.(i).(j).tempFracs.(k_max).(m_max).histoInner[*,tt] / float(nRaysMax)
      yy /= radNorm
      oplot,xx,smooth(yy,sK,/nan),line=2,color=cIndsPDF[tt]
      
      ; mark zero (at no radius is temp bigger) stats
      yy = rr.(i).(j).tempFracs.(k_max).(m_max).countZero[tt] / float(nRaysMax)
      yy /= radNorm
      ;cgPlots,[0.0],[yy],psym='filled circle',color=cIndsPDF[tt]
      print,'zero fracs: ',tempFacs[tt],yy
    endfor
    
    ; r(Tmax) line
    yy = rr.(i).(j).tempFracs.(k_max).(m_max).histoMax / float(nRaysMax)
    yy /= radNorm
    cgPlot,xx,smooth(yy,sK,/nan),line=1,color='black',/overplot
    
    ; legends
    legendStrs = textoidl("T_{ }/_{ }T_{vir} > "+string(tempFacs,format='(f4.1)'))
    legend,legendStrs,textcolor=cIndsPDF,spacing=1.4*!p.charsize,pos=[1.55,0.08],/data
    
    loadColorTable,'bw linear'
    legend,textoidl(['r_{max} ( T / T_{vir} > x )',$
                     'r_{min} ( T / T_{vir} > x )',$
                     'r ( T_{max} )']),$
      linestyle=[0,2,1],spacing=1.4*!p.charsize,pos=[0.5,yrangePDF[1]*0.07],/data
  
  end_PS
  
  ; plot (6b) - cumulative ray PDF (fraction of sphere covered) based on temp criterion (indiv halo)
  start_PS, rr.(i).(j).sP.plotPath + 'rayRadHisto_TempFacsCumulative_' + plotStrIndivMax + '.eps'
  
    yrange2 = [0.0,1.05]
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange2,/xs,/ys,xlog=xlog,$
      xtitle="",ytitle="",xtickname=replicate(" ",10),ytickname=replicate(" ",10),/noerase

    foreach radLine,[0.5,1.0,1.5] do $
      cgPlot,[radLine,radLine],yrange2,line=1,thick=!p.thick-1.0,color='light gray',/overplot
    foreach radLine,[0.25,0.75,1.25,1.75] do $
      cgPlot,[radLine,radLine],yrange2,line=1,thick=!p.thick-2.0,color='light gray',/overplot
    foreach radLine,[0.2,0.4,0.6,0.8] do $
      cgPlot,xrange,[radLine,radLine],line=1,thick=!p.thick-2.0,color='light gray',/overplot
      
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange2,/xs,/ys,xlog=xlog,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl("Cumulative Solid Angle (\Omega / 4\pi)"),$
      /noerase
      
    loadColorTable,'blue-red2'
      
    xx = rr.(i).(j).hsv.temp.(k_max).(m_max).radFacs
    
    ; plot
    for tt=0,n_elements(tempFacs)-1 do begin
      ; rmax
      yy = rr.(i).(j).tempFracs.(k_max).(m_max).histoOuter[*,tt] / float(nRaysMax)
      oplot,xx,total(yy,/cum,/pres),line=0,color=cIndsPDF[tt] ;,/overplot
      
      ; rmin
      yy = rr.(i).(j).tempFracs.(k_max).(m_max).histoInner[*,tt] / float(nRaysMax)
      oplot,xx,total(yy,/cum,/pres),line=2,color=cIndsPDF[tt] ;,/overplot
      
    endfor
    
    ; r(Tmax) line
    yy = rr.(i).(j).tempFracs.(k_max).(m_max).histoMax / float(nRaysMax)
    cgPlot,xx,total(yy,/cum,/pres),line=1,color='black',/overplot

    ; legends    
    legendStrs = textoidl("T_{ }/_{ }T_{vir} > "+string(tempFacs,format='(f4.1)'))
    legend,legendStrs,textcolor=cIndsPDF,spacing=1.4*!p.charsize,pos=[1.5,0.6],/data
    
    loadColorTable,'bw linear'
    legend,textoidl(['r_{max }|_{ }T / T_{vir} > x',$
                     'r_{min }|_{ }T / T_{vir} > x',$
                     'r_{ }|_{ }T=T_{max}']),$
      linestyle=[0,2,1],spacing=1.4*!p.charsize,pos=[1.3,0.22],/data
  
  end_PS
  
  ; plot (6c) - ray PDF (fraction of sphere covered) based on temp criterion (Nrad,Ns convergence)
  start_PS, rr.(0).(0).sP.plotPath + 'rayRadHisto_TempFacs_NrNsConv_' + plotStrIndiv2 + '.eps',$
    xs=7.5, ys=6.5
  
    pos = [0.14,0.09,0.95,0.96]
    subSize = 0.18
    offset = [0,-0+subSize*2,0,0]  
  
    ; (i) - main panel
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangePDF*[1,0.1],/xs,/ys,xlog=xlog,/ylog,yminor=5,$
      xtitle="",xtickname=replicate(" ",10),ytitle=textoidl(ytitleSA),pos=pos+offset,/noerase

    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrangePDF,line=1,color='light gray',/overplot
      
    loadColorTable,'blue-red2'    
      
    ; plot
    foreach targetFac,pdfTvirFacs do begin
      tt = where(tempFacs eq targetFac)
      
      for kk=0,n_elements(Nside)-1 do begin
        for mm=0,n_elements(Nrad)-1 do begin
          ; normalization
          nRaysLoc = rr.(i).(j).hsv.temp.(kk).(mm).nPx
          radBinSize = xrange[1]/Nrad[mm]
          radNorm = radBinSize/radNormFac
          
          ; rmax
          xx = rr.(i).(j).hsv.temp.(kk).(mm).radFacs
          yy = rr.(i).(j).tempFracs.(kk).(mm).histoOuter[*,tt] / float(nRaysLoc)
          yy /= radNorm
          oplot,xx,smooth(yy,sK,/nan),line=nrLines[mm],color=cIndsConv[kk]
        endfor ;mm
      endfor ;kk
      
    endforeach
    
    ; mark targetFacs
    cgText,0.5,0.03,textoidl('r_{max }|_{ }T / T_{vir} > '+string(pdfTvirFacs[0],format='(f3.1)'))
    cgText,1.45,0.025,textoidl('r_{max }|_{ }T / T_{vir} > '+string(pdfTvirFacs[1],format='(f3.1)'))
      
    ; legends
    legendStrs = textoidl("N_{side} = "+string(Nside,format='(I2)'))
    legend,legendStrs,textcolor=cIndsConv,spacing=1.4*!p.charsize,pos=[0.1,0.001],/data
    
    loadColorTable,'bw linear'
    legend,textoidl('N_{rad} = '+str(Nrad)),$
      linestyle=nrLines[0:n_elements(Nrad)-1],spacing=1.4*!p.charsize,pos=[1.4,0.0007],/data
  
    ; (ii) - two subpanels
    for tt_ind=0,1 do begin
    
      yrangeSub = [0.2,5.0]
      ytickvSub = [0.3,1.0,3.0]
      ytitleSub = "\Omega_{N_{rad}} / \Omega_{"+str(Nrad[-1])+"}" ;"Ratio"
      posSub = pos+offset
      posSub[1] -= subSize*(tt_ind+1) ; move y of lowerleft down
      posSub[3] = posSub[1] + subSize ; set height of subpanel
      
      xtitleSub = ""
      xticknSub = replicate(" ",10)
      if tt_ind eq 1 then xtitleSub = textoidl("r / r_{vir}")
      if tt_ind eq 1 then xticknSub = ''
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeSub,/xs,/ys,xlog=xlog,/ylog,yminor=1,$
        xtitle=xtitleSub,ytitle=textoidl(ytitleSub),ytickv=ytickvSub,yticks=2,$
        xtickname=xticknSub,pos=posSub,/noerase
        
      foreach ytick,[1.0] do $
        cgPlot,xrange,replicate(ytick,2),line=0,color=cgColor('light gray'),thick=!p.thick-1.0,/overplot
    
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeSub,/xs,/ys,xlog=xlog,/ylog,yminor=1,$
        xtitle="",ytitle="",ytickv=ytickvSub,yticks=2,pos=posSub,/noerase,$
        xtickname=replicate(" ",10),ytickname=replicate(" ",10)
    
      loadColorTable,'blue-red2'
    
      tt = where(tempFacs eq pdfTvirFacs[tt_ind])
      
      for kk=0,n_elements(Nside)-1 do begin
      
        ; Nrad=max values
        yy_max = rr.(i).(j).tempFracs.(kk).(m_max).histoOuter[*,tt]
        radNorm = (xrange[1]/Nrad[m_max]/radNormFac)[0]
        yy_max /= radNorm
      
        ; loop over 100,200 for ratios
        for mm=0,n_elements(Nrad)-2 do begin
          ; normalization
          radBinSize = xrange[1]/Nrad[mm]
          radNorm = radBinSize/radNormFac
          
          ; rmax
          xx = rr.(i).(j).hsv.temp.(kk).(mm).radFacs
          yy = rr.(i).(j).tempFracs.(kk).(mm).histoOuter[*,tt]
          yy /= radNorm
          yy /= rebin(yy_max,n_elements(yy)) ; neighbor average yy_400 down to yy number of rad pts
          
          oplot,xx,smooth(yy,sK,/nan),line=nrLines[mm],color=cIndsConv[kk]
        endfor ;mm
      endfor ;kk
      
      legend,textoidl('r_{max }|_{ }T / T_{vir} > '+string(pdfTvirFacs[tt_ind],format='(f3.1)')),$
        top=(tt_ind eq 0),right=(tt_ind eq 0),bottom=(tt_ind eq 1),left=(tt_ind eq 1)
    endfor ;tt_ind
  
  end_PS
  
  ; plot (6d) - ray PDF (fraction of sphere covered) based on temp criterion
  ; all haloes, L11 in main panel, ratios of L9/L11 and L10/L11 in subpanel
  start_PS, rr.(0).(0).sP.plotPath + 'rayRadHisto_TempFacsSub_' + plotStr + '.eps', xs=7.5, ys=8.0
  
    pos = [0.14,0.09,0.95,0.96]
    subSize = 0.30
    offset = [0,-0+subSize,0,0]  
  
    ; (i) - main panel
    yrange2 = [0.0,0.03]
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange2,/xs,/ys,xlog=xlog,$
      xtitle="",xtickname=replicate(" ",10),ytitle=textoidl(ytitleSA),pos=pos+offset,/noerase

    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrangePDF,line=1,color='light gray',/overplot
      
    for ii=0,n_tags(rr)-1 do begin
    
      jj = where( tag_names(rr.(i)) eq 'L11', count_ind )
      if count_ind eq 0 then continue
                
      ; normalization
      radBinSize = xrange[1]/Nrad[m_max]
      radNorm = radBinSize/radNormFac
      
      xx = rr.(ii).(jj).hsv.temp.(k_max).(m_max).radFacs
      co = rr.(ii).(jj).sP.colors[cInds[jj]]
      th = !p.thick * rThick[jj]
      
      ; plot 
      foreach targetFac,pdfTvirFacs,t do begin
        tt = where(tempFacs eq targetFac)
        
        ; rmax
        yy = rr.(ii).(jj).tempFracs.(k_max).(m_max).histoOuter[*,tt] / float(nRaysMax)
        yy /= radNorm
        if targetFac eq 1.0 then yy /= 2.0 ; note! for display on same y-axis
        cgPlot,xx[0:-2],smooth(yy[0:-2],sK+2,/nan),line=ttLines[t],color=co,thick=th,/overplot
        
        ; mark maxima locations?
        if targetFac eq 0.2 then begin
          yy_loc = smooth(yy[0:-2],sK+2,/nan)
          ww = where(yy_loc eq max(yy_loc),count)
          if count ne 1 then print,'WARNING: maxima count: '+str(count)
          ;print,rr.(ii).(0).sP.hInd,xx[ww[0]]
          cgPlot,replicate(xx[ww[0]],2),[0.0195,0.021],line=ttLines[t],color=co,thick=th,/overplot
        endif
      endforeach ;targetFac
    endfor ;n_tags(rr),ii
    
    ; mark targetFacs
    cgText,0.35,0.025,textoidl('r_{max }|_{ }T / T_{vir} > '+string(pdfTvirFacs[0],format='(f3.1)'))
    cgText,1.32,0.018,textoidl('r_{max }|_{ }T / T_{vir} > '+string(pdfTvirFacs[1],format='(f3.1)'))
    
    ; legends
    loadColorTable,'bw linear'    
    
    legend,hNames[0:3],textcolor=hColors[0:3],pos=[1.70,0.0268],/data
    legend,hNames[4:*],textcolor=hColors[4:*],pos=[1.78,0.0268],/data
    legend,rNames[2],linestyle=rLines[2],thick=rThick[2]*!p.thick,/top,/right
    
    ; (ii) - subpanel
    yrangeSub = [0.1,10.0]
    ytickvSub = [0.1,0.2,0.5,1.0,2.0,5.0]
    ytitleSub = "Ln/L11" ;"Ratio"
    posSub = pos+offset
    posSub[1] -= subSize ; move y of lowerleft down
    posSub[3] = posSub[1] + subSize ; set height of subpanel
        
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeSub,/xs,/ys,xlog=xlog,/ylog,yminor=1,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl(ytitleSub),$
      ytickv=ytickvSub,yticks=n_elements(ytickvSub)-1,$
      pos=posSub,/noerase
      
    foreach ytick,[ytickvSub] do $
      cgPlot,xrange,replicate(ytick,2),line=0,color=cgColor('light gray'),thick=!p.thick-1.0,/overplot
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeSub,/xs,/ys,xlog=xlog,/ylog,yminor=1,$
      xtitle="",ytitle="",ytickv=ytickvSub,yticks=n_elements(ytickvSub)-1,pos=posSub,/noerase,$
      xtickname=replicate(" ",10),ytickname=replicate(" ",10)      
  
    yy_stack = fltarr(3,2,nrIndiv)
  
    for ii=0,n_tags(rr)-1 do begin
      jj_L11 = where( tag_names(rr.(i)) eq 'L11', count_ind )
      if count_ind eq 0 then continue
      
      foreach targetFac,pdfTvirFacs,t do begin
        
        tt = where(tempFacs eq targetFac)
        yy_L11 = rr.(ii).(jj_L11).tempFracs.(k_max).(m_max).histoOuter[*,tt]
        
        yy_stack[2,t,*] += reform(yy_L11)
        
        for jj=0,n_tags(rr.(i))-2 do begin

          xx = rr.(ii).(jj).hsv.temp.(k_max).(m_max).radFacs
          co = rr.(ii).(jj).sP.colors[cInds[jj]]
          th = !p.thick * rThick[jj]
          
          ; rmax
          yy = rr.(ii).(jj).tempFracs.(k_max).(m_max).histoOuter[*,tt]
          yy_stack[jj,t,*] += reform(yy)
          
          yy /= float(yy_L11)
          
          if targetFac eq 0.2 then $
            cgPlot,xx,smooth(yy,sK+2,/nan),line=ttLines[t],color=co,thick=th,/overplot

        endfor ;n_tags(rr.(ii)),jj
      
      endforeach ;targetFac
      
    endfor ;n_tags(rr),ii
    
    ; halo-mean
    for jj=0,n_tags(rr.(i))-2 do begin
      foreach targetFac,pdfTvirFacs,t do begin
        yy_loc = reform(yy_stack[jj,t,*])
        yy_norm = reform(yy_stack[2,t,*])
        
        ww = -1
        if targetFac eq 1.0 then ww = where(xx le 1.5)
        if targetFac eq 0.2 then ww = where(xx ge 0.4)
        
        cgPlot,xx[ww],smooth(yy_loc[ww]/yy_norm[ww],sK+2,/nan),$
          line=ttLines[t],color='black',thick=rThick[jj]*!p.thick*3,/overplot
      endforeach
    endfor

    legend,[rNames[0:1]+'/'+rNames[2]],linestyle=rLines[0:1],thick=rThick[0:1]*!p.thick,/top,/right
  
  end_PS
  endif ;skipPDFs
  
  stop

end
