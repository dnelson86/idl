; zoomRays.pro
; 'zoom project' radial rays at different angles
; dnelson nov.2014

; plotZoomRadialRays(): consider individual radial rays evenly distributed in (theta,phi)

pro plotZoomRadialRays
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  hInds     = [0]
  resLevels = [9]
  redshift  = 2.0
  newSaves  = 0 ; override existing saves
  
  ; binning config
  Nside   = 8 ; healpix resolution parameter, 8=768, 16~3k, 32~12k, 64~50k, 128~200k
  radFacs = linspace(0.01,2.0,100) ; r/rvir, sampling points
  cutSubS = 2 ; remove substructures? (0=no, 1=satellites, 2=all subfind other than target group)

  ; plot config
  xrange_rad   = [0.0,2.0]   ; r/rvir
  xrange_log   = 0           ; plot in log? doesn't necessary have to equal binInLog  
  yrange_temp  = [1e4,3.16e6]   ; K
  yrange_dens  = [-28,-24]   ; log cgs (g/cm^3)
  yrange_vrad  = [-350,150]  ; km/s
  yrange_csize = [0.01,2.0]  ; physical kpc
  yrange_entr  = [1e7,1e9]   ; cgs (K cm^2)
  yrange_angm  = [1e3,1e5]   ; kpc km/s
  
  lines    = [0,1,2] ;[1,2,0] ; line style, one per resLevel
  cInds    = [1,1,1] ; color index, one per resLevel
  radLines = [0.15,1.5]
  sK       = 1
  thick    = 0.1
  
  ; load
  foreach hInd,hInds do begin
    hLoc = {}
        
    foreach resLevel,resLevels do begin
      ; get list of snapshots for this halo in the redshift range
      rLoc = {}
      
      if resLevel eq 11 and total(hInd eq [3,4,6]) gt 0 then continue
      
      sP = simParams(run='zoom_20Mpc',res=resLevel,hInd=hInd,redshift=redshift)     
      gcInd = zoomTargetHalo(sP=sP)
      
      ; calculate radial rays using nested healpix spheres
      hsv_temp = haloShellValue(sP=sP, partType='gas', valName='temp', subgroupIDs=[gcInd], $
                                Nside=Nside, radFacs=radFacs, cutSubS=cutSubS, newSaves=newSaves)
      hsv_temp.value = 10.0^hsv_temp.value ; remove log
      
      hsv_vrad = haloShellValue(sP=sP, partType='gas', valName='radvel', subgroupIDs=[gcInd], $
                                Nside=Nside, radFacs=radFacs, cutSubS=cutSubS, newSaves=newSaves)
              
      rLoc = mod_struct( rLoc, 'sP', sP )
      rLoc = mod_struct( rLoc, 'hsv_temp', hsv_temp )
      rLoc = mod_struct( rLoc, 'hsv_vrad', hsv_vrad )
    
      ; add this resolution level to the halo struct
      hLoc = mod_struct( hLoc, 'L'+str(resLevel), rLoc )     
     
    endforeach ;resLevels
    
    ; add this halo to the save struct
    rr = mod_struct( rr, 'h'+str(hInd), hLoc )
    
  endforeach ;hInds
  
  ; plot setup
  plotStr = 'h'
  hNames  = 'h' + str(hInds)
  hColors = []
  rNames  = 'L' + str(resLevels)
  rLines  = lines[0:n_elements(resLevels)-1]

  foreach hInd,hInds do plotStr += str(hInd)
  foreach hInd,hInds,i do hColors = [hColors,rr.(i).(0).sP.colors[cInds[-1]]]
  plotStr += '_L'
  foreach resLevel,resLevels do plotStr += str(resLevel)
  plotStr += '_cutSubS-' + str(cutSubS)

  if xrange_log eq 1 then xrange_rad[0] = 0.01
  if xrange_log eq 1 then plotStr += '_xlog'
  
  ; establish color mapping for temp
  set_plot,'ps'
  mbColors = sampleColorTable('blue-red2', 256, bounds=[0.0,1.0])
  loadColorTable,'bw linear'
  
  pos = [0.12,0.12,0.87,0.9]
  pos_cbar = [0.89,0.12,0.92,0.9]

  ; plot (1) - temperature individual rays
  start_PS, rr.(0).(0).sP.plotPath + 'radialRays_TempIndiv.' + plotStr + '.eps' ;, xs=12.0*1.5, ys=6.0*1.5
  
    ; (a) - gas temperature
    cgPlot,[0],[0],/nodata,xrange=xrange_rad,yrange=yrange_temp,/xs,/ys,xlog=xrange_log,/ylog,yminor=5,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl("T_{gas} [_{ }K_{ }]"),/noerase,pos=pos

    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrange_temp+[0.05,-0.05],line=1,color='light gray',/overplot
      
    for i=0,n_tags(rr)-1 do begin
      for j=0,n_tags(rr.(i))-1 do begin
        xx = rr.(i).(j).hsv_temp.radFacs
        co = rr.(i).(j).sP.colors[cInds[j]]
        
        ; plot each radial ray, corresponding to one healpix pixel
        for k=0,rr.(i).(j).hsv_temp.nPx-1 do begin
          yy = reform( rr.(i).(j).hsv_temp.value[k,*] )
          cgPlot,xx,smooth(yy[*,0],sK),line=lines[j],color=co,/overplot,thick=thick
        endfor
        
      endfor ;n_tags(rr.(i)),j
    endfor ;n_tags(rr),i
    
    legend,hNames,textcolor=hColors,/top,/right
    legend,rNames,linestyle=rLines,/bottom,/left
      
  end_PS
  
  ; plot (1d)
  start_PS, rr.(0).(0).sP.plotPath + 'radialRays_RadHisto.' + plotStr + '.eps' ;, xs=12.0*1.5, ys=6.0*1.5
  
    ; (a) - histo
    yrange2 = [1e-3,1.0]
    cgPlot,[0],[0],/nodata,xrange=xrange_rad,yrange=yrange2,/xs,/ys,xlog=xrange_log,/ylog,yminor=5,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl("PDF"),/noerase,pos=pos

    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrange_temp+[0.05,-0.05],line=1,color='light gray',/overplot
      
    for i=0,n_tags(rr)-1 do begin
      for j=0,n_tags(rr.(i))-1 do begin
        xx = rr.(i).(j).hsv_temp.radFacs
        co = rr.(i).(j).sP.colors[cInds[j]]
        
        countZero  = 0
        histoOuter = lonarr(rr.(i).(j).hsv_temp.nRadFacs)
        histoInner = lonarr(rr.(i).(j).hsv_temp.nRadFacs)
        
        ; analysis: largest or smallest radii where each ray surpasses Tvir/2
        haloTvirSearch = 1.0 * codeMassToVirTemp( rr.(i).(j).sP.targetHaloMass, sP=rr.(i).(j).sP )
        
        for k=0,rr.(i).(j).hsv_temp.nPx-1 do begin
          yy = reform( rr.(i).(j).hsv_temp.value[k,*] )
          w = where(yy gt 1e5,count)
          if count eq 0 then begin
            countZero += 1
          endif else begin
            histoOuter[max(w)] += 1
            histoInner[min(w)] += 1
          endelse
        endfor
        
        cgPlot,xx,histoOuter/float(max(histoOuter)),line=0,color=co,/overplot,thick=thick
        cgPlot,xx,histoInner/float(max(histoInner)),line=2,color=co,/overplot,thick=thick
        
      endfor ;n_tags(rr.(i)),j
    endfor ;n_tags(rr),i
    
    legend,hNames,textcolor=hColors,/top,/right
    legend,['0.5tvir outer','0.5tvir inner'],linestyle=[0,2],/bottom,/left
  
  end_PS
  
  ; plot (1b) - temperature individual rays, colored by vrad
  start_PS, rr.(0).(0).sP.plotPath + 'radialRays_TempIndivColor.' + plotStr + '.eps' ;, xs=12.0*1.5, ys=6.0*1.5
   
    ; (a) - gas temperature
    cgPlot,[0],[0],/nodata,xrange=xrange_rad,yrange=yrange_temp,/xs,/ys,xlog=xrange_log,/ylog,yminor=5,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl("T_{gas} [_{ }K_{ }]"),/noerase,pos=pos

    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrange_temp+[0.05,-0.05],line=1,color='light gray',/overplot
      
    for i=0,n_tags(rr)-1 do begin
      for j=0,n_tags(rr.(i))-1 do begin
        xx = rr.(i).(j).hsv_temp.radFacs
        co_mm = yrange_vrad
        co_vv = reform( rr.(i).(j).hsv_vrad.value )
        co_vv = bytscl( co_vv, min=co_mm[0], max=co_mm[1] )
        
        ; plot each radial ray, corresponding to one healpix pixel
        for k=0,40 do begin ;rr.(i).(j).hsv_temp.nPx-1 do begin
          yy = reform( rr.(i).(j).hsv_temp.value[k,*] )
          
          ; use temperature at each radial sample to color segment, map into [0,255]
          co = mbColors[co_vv[k,*]]
          
          for m=0,n_elements(yy)-2 do begin
            xx0 = xx[m]
            xx1 = xx[m+1]
            yy0 = yy[m]
            yy1 = yy[m+1]
            cgPlot,[xx0,xx1],[yy0,yy1],line=lines[j],color=co[m],/overplot,thick=thick
          endfor ;m
        
        endfor ;k
      endfor ;n_tags(rr.(i)),j
    endfor ;n_tags(rr),i
    
    legend,hNames,textcolor=hColors,/top,/right
    legend,rNames,linestyle=rLines,/bottom,/left
      
    loadColorTable,'blue-red2'
    cgColorbar,pos=pos_cbar,divisions=5,/vertical,/right,range=co_mm,title=textoidl("v_{rad} [_{ }km/s_{ }]")
      
  end_PS
  
  ; plot (2) - vrad individual rays, colored by temperature
  start_PS, rr.(0).(0).sP.plotPath + 'radialRays_VRadIndivColor.' + plotStr + '.eps' ;, xs=12.0*1.5, ys=6.0*1.5
  
    ; (a) - gas radial velocity
    cgPlot,[0],[0],/nodata,xrange=xrange_rad,yrange=yrange_vrad,/xs,/ys,xlog=xrange_log,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl("v_{rad} [_{ }km/s_{ }]"),/noerase,pos=pos

    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrange_vrad+[20,-20],line=1,color='light gray',/overplot
      
    for i=0,n_tags(rr)-1 do begin
      for j=0,n_tags(rr.(i))-1 do begin
        xx = rr.(i).(j).hsv_vrad.radFacs
        co_mm = alog10( yrange_temp )
        co_vv = alog10( reform( rr.(i).(j).hsv_temp.value ) )
        co_vv = bytscl( co_vv, min=co_mm[0], max=co_mm[1] )
        
        if ~array_equal( size(rr.(i).(j).hsv_vrad.value), size(rr.(i).(j).hsv_temp.value) ) then message,'Error'
        
        ; plot each radial ray, corresponding to one healpix pixel
        for k=0,40 do begin ;rr.(i).(j).hsv_temp.nPx-1 do begin
          yy = reform( rr.(i).(j).hsv_vrad.value[k,*] )
          
          ; use temperature at each radial sample to color segment, map into [0,255]
          co = mbColors[co_vv[k,*]]
          
          for m=0,n_elements(yy)-2 do begin
            xx0 = xx[m]
            xx1 = xx[m+1]
            yy0 = yy[m]
            yy1 = yy[m+1]
            cgPlot,[xx0,xx1],[yy0,yy1],line=lines[j],color=co[m],/overplot,thick=thick
          endfor ;m
          
        endfor ;k
      endfor ;n_tags(rr.(i)),j
    endfor ;n_tags(rr),i
    
    legend,hNames,textcolor=hColors,/top,/right
    legend,rNames,linestyle=rLines,/top,/left
      
    loadColorTable,'blue-red2'
    cgColorbar,pos=pos_cbar,divisions=5,/vertical,/right,range=co_mm,title=textoidl("T_{gas} [_{ }log K_{ }]")
      
  end_PS
  stop

end
