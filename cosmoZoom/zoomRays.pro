; zoomRays.pro
; 'zoom project' radial rays at different angles
; dnelson dec.2014

; plotZoomRadialRays(): visualize individual/aggregate radial rays

pro plotZoomRadialRays
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  hInds     = [0]
  resLevels = [9,10,11]
  redshift  = 2.0
  newSaves  = 0 ; override existing saves
  
  ; binning config
  Nside   = [8,16,32,64]  ; healpix resolution parameter, 2=12*4, 8=768, 16~3k, 32~12k, 64~50k, 128~200k
  Nrad    = [100] ; [200,400] number of radial sampling points
  cutSubS = 2 ; remove substructures? (0=no, 1=satellites, 2=all subfind other than target group)

  ; plots of one halo at one resolution, one binning (Nside/Nrad), which?
  hIndiv  = 0   ; which halo to plot individuals from
  rIndiv  = 11  ; which resolution level to plot individuals from
  nsIndiv = 64  ; Nside indiv
  nrIndiv = 100 ; Nrad indiv
  skipRaw = 1   ; skip plots (1) and (2) showing all rays individually (not a good idea Ns>=32)
  
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
  yrangePDF = [1e-4,1e0] ; fractional/cumulative solid angle ray analysis  
  tempFacs  = [2.0,1.5,1.0,0.5,0.2,0.1]
  ytitleAng = textoidl("Angular Direction (\theta, \phi)")
  
  cIndsPDF  = fix(reverse(linspace(10,225,n_elements(tempFacs))))
  cIndsConv = fix([10,70,180,230]) ; colors for Nside convergence
  nrLines   = [2,1,0] ; linestyles for Nrad convergence
  
  ; analysis parameters
  haloTvirSearch = 1.0 * codeMassToVirTemp( logMsunToCodeMass(12.0), redshift=redshift )
  print,'haloTvirSearch: '+str(alog10(haloTvirSearch))
  
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
  saveFilename = '/n/home07/dnelson/plots/temp_plotZoomRadialRays.sav'
  if ~file_test(saveFilename) or writeNewTempFile eq 1 then begin
  foreach hInd,hInds do begin
    hLoc = {}
        
    foreach resLevel,resLevels do begin
      ; get list of snapshots for this halo in the redshift range
      rLoc = {}
      
      if resLevel eq 11 and total(hInd eq [3,4,6]) gt 0 then continue
      
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
              
              w = where(yy gt tempFacs[m]*haloTvirSearch,count)
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
    
    save,rr,filename=saveFilename
    print,'Saved: ['+saveFilename+']'
  endforeach ;hInds
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
    
  ; plot (1) - individual rays, one plot per quantity which is plotting on the y-axis, vs radius on the x
  ;   for each plot, four panels (4 different colorings)
  nRaysTot = rr.(i).(j).hsv.(0).(k).(m).nPx
  nRaysPerBin = ceil(nRaysTot / float(nRayBins))
  
  print,'nRaysTot: ',nRaysTot,' nRaysPerBin: ',nRaysPerBin
  
  ; sub-divide all rays into some bins, for visual inspection
  if skipRaw ne 1 then begin
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
          
        ; select data
        xx = rr.(i).(j).hsv.(count).(k).(m).radFacs
        
        ; colormap
        co_ind = where( tag_names(config) eq coName )
        loadColorTable,config.(co_ind).ctName
        
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
          
          ; use temperature at each radial sample to color segment, map into [0,255] CT indices
          co = reform(co_vv[rayNum,*])
          
          for rayPt=0,n_elements(yy)-2 do $
            oplot,[xx[rayPt],xx[rayPt+1]],[yy[rayPt],yy[rayPt+1]],color=co[rayPt],line=0,thick=thick
        
        endfor ;rayNum
          
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
          
          ; use temperature at each radial sample to color segment, map into [0,255] CT indices
          co = reform(co_vv[rayNum,*])
          
          for rayPt=0,n_elements(yy)-2 do $
            oplot,[xx[rayPt],xx[rayPt+1]],[yy[rayPt],yy[rayPt+1]],color=co[rayPt],line=0,thick=thick
        
        endfor ;rayNum
      
      endfor ; rayBins
      
      posCbar = [0.90,0.08,0.90+cbarSize,0.98]+offset
      cgColorbar,pos=posCbar,divisions=config.(co_ind).cbarDivs,/vertical,/right,range=co_mm,$
          title=textoidl(config.(co_ind).ytitle),charsize=!p.charsize*1.0
      
      end_PS
      
    endforeach ; coNames
  endforeach ; tag_names(config),yName
  endif ;skipRaw
  
  ; plot (3a) - global 2D plot of all rays, single coloring per plot, one plot per value
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
  
  ; plot (3b) - all values on one plot (5 narrow columns)
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
      
  ; plot (5a) - ray PDF (fraction of sphere covered) based on temp criterion (indiv halo)
  start_PS, rr.(i).(j).sP.plotPath + 'rayRadHisto_TempFacs_' + plotStrIndiv + '.eps'
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangePDF,/xs,/ys,xlog=xlog,/ylog,yminor=5,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl("Fractional Solid Angle (\Omega / 4\pi)"),$
      /noerase

    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrangePDF,line=1,color='light gray',/overplot
      
    loadColorTable,'blue-red2'
      
    xx = rr.(i).(j).hsv.temp.(k).(m).radFacs
        
    ; plot
    for tt=0,n_elements(tempFacs)-1 do begin
      ; rmax
      yy = rr.(i).(j).tempFracs.(k).(m).histoOuter[*,tt] / float(nRaysTot)
      oplot,xx,smooth(yy,sK,/nan),line=0,color=cIndsPDF[tt]
      
      ; rmin
      yy = rr.(i).(j).tempFracs.(k).(m).histoInner[*,tt] / float(nRaysTot)
      oplot,xx,smooth(yy,sK,/nan),line=2,color=cIndsPDF[tt]
      
      ; mark zero (at no radius is temp bigger) stats
      yy = rr.(i).(j).tempFracs.(k).(m).countZero[tt] / float(nRaysTot)
      ;cgPlots,[0.0],[yy],psym='filled circle',color=cIndsPDF[tt]
      print,'zero fracs: ',tempFacs[tt],yy
    endfor
    
    ; r(Tmax) line
    yy = rr.(i).(j).tempFracs.(k).(m).histoMax / float(nRaysTot)
    cgPlot,xx,smooth(yy,sK,/nan),line=1,color='black',/overplot
    
    ; legends
    legendStrs = textoidl("T_{ }/_{ }T_{vir} > "+string(tempFacs,format='(f4.1)'))
    legend,legendStrs,textcolor=cIndsPDF,spacing=1.4*!p.charsize,pos=[1.55,0.8],/data
    
    loadColorTable,'bw linear'
    legend,textoidl(['r_{max} ( T / T_{vir} > x )',$
                     'r_{min} ( T / T_{vir} > x )',$
                     'r ( T_{max} )']),$
      linestyle=[0,2,1],spacing=1.4*!p.charsize,pos=[0.5,yrangePDF[1]*0.7],/data
  
  end_PS
  
  ; plot (5b) - cumulative ray PDF (fraction of sphere covered) based on temp criterion (indiv halo)
  start_PS, rr.(i).(j).sP.plotPath + 'rayRadHisto_TempFacsCumulative_' + plotStrIndiv + '.eps'
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangePDF*[10,1.2],/xs,/ys,xlog=xlog,/ylog,yminor=5,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl("Cumulative Solid Angle (\Omega / 4\pi)"),$
      /noerase

    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrangePDF,line=1,color='light gray',/overplot
      
    loadColorTable,'blue-red2'
      
    xx = rr.(i).(j).hsv.temp.(k).(m).radFacs
    
    ; plot
    for tt=0,n_elements(tempFacs)-1 do begin
      ; rmax
      yy = rr.(i).(j).tempFracs.(k).(m).histoOuter[*,tt] / float(nRaysTot)
      oplot,xx,total(yy,/cum,/pres),line=0,color=cIndsPDF[tt] ;,/overplot
      
      ; rmin
      yy = rr.(i).(j).tempFracs.(k).(m).histoInner[*,tt] / float(nRaysTot)
      oplot,xx,total(yy,/cum,/pres),line=2,color=cIndsPDF[tt] ;,/overplot
      
    endfor
    
    ; r(Tmax) line
    yy = rr.(i).(j).tempFracs.(k).(m).histoMax / float(nRaysTot)
    cgPlot,xx,total(yy,/cum,/pres),line=1,color='black',/overplot

    ; legends    
    legendStrs = textoidl("T_{ }/_{ }T_{vir} > "+string(tempFacs,format='(f4.1)'))
    legend,legendStrs,textcolor=cIndsPDF,spacing=1.4*!p.charsize,pos=[1.5,0.1],/data
    
    loadColorTable,'bw linear'
    legend,textoidl(['r_{max }|_{ }T / T_{vir} > x',$
                     'r_{min }|_{ }T / T_{vir} > x',$
                     'r_{ }|_{ }T=T_{max}']),$
      linestyle=[0,2,1],spacing=1.4*!p.charsize,pos=[1.3,0.005],/data
  
  end_PS
  
  ; plot (5c) - ray PDF (fraction of sphere covered) based on temp criterion (Nrad,Ns convergence)
  start_PS, rr.(0).(0).sP.plotPath + 'rayRadHisto_TempFacs_NrNsConv_' + plotStrIndiv + '.eps'
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangePDF,/xs,/ys,xlog=xlog,/ylog,yminor=5,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl("Fractional Solid Angle (\Omega / 4\pi)"),$
      /noerase

    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrangePDF,line=1,color='light gray',/overplot
      
    loadColorTable,'blue-red2'
      
    xx = rr.(i).(j).hsv.temp.(k).(m).radFacs
    
    ; plot
    foreach targetFac,[1.0,0.2] do begin
      tt = where(tempFacs eq targetFac)
      
      for kk=0,n_elements(Nside)-1 do begin
        for mm=0,n_elements(Nrad)-1 do begin
          ; normalization
          nRaysLoc = rr.(i).(j).hsv.temp.(kk).(mm).nPx
          ; todo rad normalization
          
          ; rmax
          yy = rr.(i).(j).tempFracs.(kk).(mm).histoOuter[*,tt] / float(nRaysLoc)
          oplot,xx,smooth(yy,sK,/nan),line=nrLines[mm],color=cIndsConv[kk]
       
          ; r(Tmax) line
          ;yy = rr.(i).(j).tempFracs.(kk).(mm).histoMax / float(nRaysLoc)
          ;oplot,xx,smooth(yy,sK,/nan),line=nrLines[mm],color=cIndsConv[kk],thick=!p.thick-2.0
        endfor ;mm
      endfor ;kk
      
    endforeach
    
    ; mark targetFacs
    ;cgText,0.1,0.2,textoidl('r_{ }|_{ }T=T_{max}'),color='black'
    cgText,0.5,0.03,textoidl('r_{max }|_{ }T / T_{vir} > 1.0'),color='black'
    cgText,1.2,0.008,textoidl('r_{max }|_{ }T / T_{vir} > 0.2'),color='black'
      
    ; legends
    legendStrs = textoidl("N_{side} = "+string(Nside,format='(I2)'))
    legend,legendStrs,textcolor=cIndsConv,spacing=1.4*!p.charsize,pos=[1.6,0.8],/data
    
    loadColorTable,'bw linear'
    legend,textoidl('N_{rad} = '+str(Nrad)),$
      linestyle=nrLines[0:n_elements(Nrad)-1],spacing=1.4*!p.charsize,pos=[0.7,0.8],/data
  
  end_PS
  
  ; plot (5d) - ray PDF (fraction of sphere covered) based on temp criterion
  ; all haloes, L11 in main panel, ratios of L9/L11 and L10/L11 in subpanel
  start_PS, rr.(0).(0).sP.plotPath + 'rayRadHisto_TempFacsSub_' + plotStr + '.eps', xs=7.5, ys=6.0
  
    pos = [0.14,0.12,0.95,0.96]
    subSize = 0.18
    offset = [0,-0+subSize,0,0]  
  
    ; maximal Nside,Nrad
    ttLines = [1,0]
    
    kk = n_elements(Nside)-1
    mm = n_elements(Nrad)-1
  
    ; main panel
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangePDF,/xs,/ys,xlog=xlog,/ylog,yminor=5,$
      xtitle="",xtickname=replicate(" ",10),ytitle=textoidl("Fractional Solid Angle (\Omega / 4\pi)"),$
      pos=pos+offset,/noerase

    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrangePDF,line=1,color='light gray',/overplot
      
    for ii=0,n_tags(rr)-1 do begin
    
      jj = where( tag_names(rr.(i)) eq 'L11', count_ind )
      if count_ind eq 0 then continue
                
      ; normalization
      nRaysLoc = rr.(ii).(jj).hsv.temp.(kk).(mm).nPx
      ; todo rad normalization
      
      xx = rr.(ii).(jj).hsv.temp.(kk).(mm).radFacs
      co = rr.(ii).(jj).sP.colors[cInds[jj]]
      th = !p.thick * rThick[jj]
      
      ; plot 
      foreach targetFac,[1.0,0.2],t do begin
        tt = where(tempFacs eq targetFac)
        
        ; rmax
        yy = rr.(ii).(jj).tempFracs.(kk).(mm).histoOuter[*,tt] / float(nRaysLoc)
        cgPlot,xx,smooth(yy,sK,/nan),line=ttLines[t],color=co,thick=th,/overplot
      endforeach ;targetFac
      
      ; r(Tmax) line
      ;yy = rr.(ii).(jj).tempFracs.(kk).(mm).histoMax / float(nRaysLoc)
      ;cgPlot,xx,smooth(yy,sK,/nan),lline=lines[jj],color=co,thick=th,/overplot
    
    endfor ;n_tags(rr),ii
    
    ; mark targetFacs
    ;cgText,0.1,0.2,textoidl('r_{ }|_{ }T=T_{max}'),color='black'
    cgText,0.2,0.1,textoidl('r_{max }|_{ }T / T_{vir} > 1.0'),color='black'
    cgText,1.1,0.08,textoidl('r_{max }|_{ }T / T_{vir} > 0.2'),color='black'
    
    ; legends
    loadColorTable,'bw linear'    
    
    legend,hNames,textcolor=hColors,/top,/right   
    legend,[rNames[0:1]+'/'+rNames[2],rNames[2]],linestyle=rLines,thick=rThick*!p.thick,/bottom,/right
  
    ; subpanel
    yrangeSub = [0.05,20.0]
    ytickvSub = [0.1,1.0,10.0]
    ytitleSub = "Ln/L11" ;"Ratio"
    posSub = pos+offset
    posSub[1] -= subSize ; move y of lowerleft down
    posSub[3] = posSub[1] + subSize ; set height of subpanel
        
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeSub,/xs,/ys,xlog=xlog,/ylog,yminor=1,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl(ytitleSub),ytickv=ytickvSub,yticks=2,$
      pos=posSub,/noerase
      
    cgPlot,xrange,[1.0,1.0],line=0,color=cgColor('light gray'),/overplot
  
    for ii=0,n_tags(rr)-1 do begin
      jj_L11 = where( tag_names(rr.(i)) eq 'L11', count_ind )
      if count_ind eq 0 then continue
      
      for jj=0,n_tags(rr.(i))-2 do begin
        ; normalization
        nRaysLoc = rr.(ii).(jj).hsv.temp.(kk).(mm).nPx
        ; todo rad normalization
        
        xx = rr.(ii).(jj).hsv.temp.(kk).(mm).radFacs
        co = rr.(ii).(jj).sP.colors[cInds[jj]]
        th = !p.thick * rThick[jj]
        
        ; plot 
        foreach targetFac,[1.0,0.2],t do begin
          tt = where(tempFacs eq targetFac)
          
          ; rmax
          yy     = rr.(ii).(jj).tempFracs.(kk).(mm).histoOuter[*,tt] / float(nRaysLoc)
          yy_L11 = rr.(ii).(jj_L11).tempFracs.(kk).(mm).histoOuter[*,tt] / float(nRaysLoc)
          
          yy /= yy_L11
          
          cgPlot,xx,smooth(yy,sK,/nan),line=ttLines[t],color=co,thick=th,/overplot
        endforeach ;targetFac
        
        ; r(Tmax) line
        ;yy = rr.(ii).(jj).tempFracs.(kk).(mm).histoMax / float(nRaysLoc)
        ;cgPlot,xx,smooth(yy,sK,/nan),lline=lines[jj],color=co,thick=th,/overplot
    
      endfor ;n_tags(rr.(ii)),jj
    endfor ;n_tags(rr),ii
  
  end_PS
  
  stop

end
