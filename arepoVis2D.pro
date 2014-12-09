; arepoVis2D.pro
; visualization for 2D arepo runs
; dnelson jan.2012

; plotVoronoi2D(): plot "voronoi_mesh" output file

pro plotVoronoi2D, fileBase, m, overPlot=overPlot, fillWindow=fillWindow, $
                   zoom=zoom, pos=pos, ids=ids, boxSize=boxSize

  ;config
  nDims   = 2
  color   = fsc_color('dark gray')
  thick   = !p.thick-2.0
  
  ;load
  vor = loadVoronoiMesh(fileBase,m,nDims)
 
  ; start plot
  if (not keyword_set(overPlot)) then begin
  
    if not keyword_set(boxSize) then stop
    
    ; set ranges
    if (not keyword_set(zoom)) then begin
      xrange = [0,boxSize[0]]
      yrange = [0,boxSize[1]]
    endif else begin
      xrange = [boxSize[0]/2.0-boxSize[0]/zoom,boxSize[0]/2.0+boxSize[0]/zoom]
      yrange = [boxSize[1]/2.0-boxSize[1]/zoom,boxSize[1]/2.0+boxSize[1]/zoom]
    endelse
  
    if keyword_set(fillWindow) then $
      plot,[0],[0],/nodata,xrange=xrange,yrange=yrange, $
           xstyle=5,ystyle=5,position=[0.0,0.0,1.0,1.0] ;5=1+4 bitwise, exact+supressed
    if not keyword_set(fillWindow) then $
      plot,[0],[0],/nodata,xrange=xrange,yrange=yrange, $
           xstyle=1,ystyle=1,position=[0.05,0.05,0.95,0.95],charsize=!p.charsize-1.0,$
           ticklen=0.0001,title="snap ["+str(m)+"]"
  endif

  ; for each cell, plot edges of cell
  for j=0L, vor.nGas-1 do begin
    x = transpose(vor.xyzEdges[0, vor.edgeList[vor.nEdgesOffset[j]:vor.nEdgesOffset[j]+vor.nEdges[j]-1]])
    y = transpose(vor.xyzEdges[1, vor.edgeList[vor.nEdgesOffset[j]:vor.nEdgesOffset[j]+vor.nEdges[j]-1]])

    plots,[x,x(0)],[y,y(0)],noclip=0,color=color,thick=thick
  endfor
   
  ; if pos set, overplot positions of gas particles, color code by id
  if (keyword_set(pos) and keyword_set(ids)) then begin
    w1 = where(ids eq -1)
    w2 = where(ids eq -2)
    w3 = where(ids eq -3)
     
    fsc_plot,pos[0,w1],pos[1,w1],psym=4,symsize=0.1,color=fsc_color('crimson'),/overplot
    fsc_plot,pos[0,w2],pos[1,w2],psym=4,symsize=0.1,color=fsc_color('slate blue'),/overplot
    fsc_plot,pos[0,w3],pos[1,w3],psym=4,symsize=0.1,color=fsc_color('forest green'),/overplot
  endif

end

; plotDensityField(): plot "density_field_" or "proj_density_field_" output file

pro plotDensityField, filePath, snaps, axes=axes, writePNG=writePNG, writeJPG=writeJPG, psOut=psOut, $
                      xyScaleFac=xyScaleFac, minMax=minMax, log=log, overPlot=overPlot,$
                      twoAxesSideBySide=twoAxesSideBySide,densTempSideBySide=densTempSideBySide,$
                      sideMinMax=sideMinMax
  
  if (not keyword_set(writePNG) and not keyword_set(writeJPG) and not keyword_set(psOut)) then begin
    print,'ERROR: No output method specified.' & stop
  endif
  
  ; if snaps is an array, treat as list of snapshot numbers and loop over each
  foreach m,snaps do begin

    ; load
    if n_elements(axes) gt 0 then df = loadDensityField(filePath,m,axes=axes[0])
    if n_elements(axes) eq 0 then df = loadDensityField(filePath,m)
    
    ; log scaling of density
    if keyword_set(log) then begin
      w = where(df.dens gt 1e-8,comp=wc)
      df.dens[w] = alog10(df.dens[w])
      df.dens[wc] = min(df.dens[w])
    endif
  
    ; color tables and scaling
    ;loadColorTable, 'blue-red', rgb_table=rgbTable
    loadct,11,rgb_table=rgbTable

    ; override?
    ;w = where(df.dens ge 5.0,count,comp=wc)
    ;if count gt 0 then df.dens[w] = -1.0
    ;print,minmax(df.dens[wc])

    if not keyword_set(minMax) then colorMinMax = minmax(df.dens)
    if keyword_set(minMax) then colorMinMax = minMax
    
    print,' dens minmax:',colorMinMax
  
    ; dens(1)
    colindex = (df.dens - colorMinMax[0])/(colorMinMax[1]-colorMinMax[0])*255.0
    
    ind = where(colindex ge 256.0)
    if ind[0] ne -1 then colindex[ind] = 255.9
    ind = where(colindex lt 0)
    if ind[0] ne -1 then colindex[ind] = 0
    
    colindex = byte(colindex)
    
    ; create image
    sideFac = 1.0
    if keyword_set(twoAxesSideBySide) or keyword_set(densTempSideBySide) then sideFac = 2.0
    pic = bytarr(3,df.nPixelsXY[0]*sideFac,df.nPixelsXY[1])
  
    ; stamp in main density
    pic[0,0:df.nPixelsXY[0]-1,0:df.nPixelsXY[1]-1] = rgbTable[colindex,0]
    pic[1,0:df.nPixelsXY[0]-1,0:df.nPixelsXY[1]-1] = rgbTable[colindex,1]
    pic[2,0:df.nPixelsXY[0]-1,0:df.nPixelsXY[1]-1] = rgbTable[colindex,2]
    
    if keyword_set(twoAxesSideBySide) then begin
      ; load other projection
      df = loadDensityField(filePath,m,axes=axes[1])
    
      ; log scaling of density
      if keyword_set(log) then begin
        w = where(df.dens gt 1e-8,comp=wc)
        df.dens[w] = alog10(df.dens[w])
        df.dens[wc] = min(df.dens[w])
      endif
      
      if not keyword_set(sideMinMax) then colorMinMax = minmax(df.dens)
      if keyword_set(sideMinMax) then colorMinMax = minMax
      
      colindex = (df.dens - colorMinMax[0])/(colorMinMax[1]-colorMinMax[0])*255.0
      
      ind = where(colindex ge 256.0)
      if ind[0] ne -1 then colindex[ind] = 255.9
      ind = where(colindex lt 0)
      if ind[0] ne -1 then colindex[ind] = 0
      
      colindex = byte(colindex)
      
      ; stamp in second density
      pic[0,df.nPixelsXY[0]:2*df.nPixelsXY[0]-1,*] = rgbTable[colindex,0]
      pic[1,df.nPixelsXY[0]:2*df.nPixelsXY[0]-1,*] = rgbTable[colindex,1]
      pic[2,df.nPixelsXY[0]:2*df.nPixelsXY[0]-1,*] = rgbTable[colindex,2]
    endif
    
    if keyword_set(densTempSideBySide) then begin
      ; get temperature from already loaded data
      if not keyword_set(sideMinMax) then colorMinMax = minmax(df.temp)
      if keyword_set(sideMinMax) then colorMinMax = minMax
      
      colindex = (df.temp - colorMinMax[0])/(colorMinMax[1]-colorMinMax[0])*255.0
      
      ind = where(colindex ge 256.0)
      if ind[0] ne -1 then colindex[ind] = 255.9
      ind = where(colindex lt 0)
      if ind[0] ne -1 then colindex[ind] = 0
      
      colindex = byte(colindex)
      
      ; stamp in temperature
      pic[0,df.nPixelsXY[0]:2*df.nPixelsXY[0]-1,*] = rgbTable[colindex,0]
      pic[1,df.nPixelsXY[0]:2*df.nPixelsXY[0]-1,*] = rgbTable[colindex,1]
      pic[2,df.nPixelsXY[0]:2*df.nPixelsXY[0]-1,*] = rgbTable[colindex,2]
    endif
    
    ;rescale image
    if keyword_set(xyScaleFac) then begin
      pic = rebin(pic,3,xyScaleFac*df.nPixelsXY[0],xyScaleFac*df.nPixelsXY[1])
    endif
    
    ;output filename
    if keyword_set(axes) then axesOut = axes
    if not keyword_set(axes) then axesOut = 0
    if n_elements(axes) eq 2 then axesOut = str(axes[0])+'-'+str(axes[1])
    outputFilename = filePath + 'density_' + string(m,format='(i3.3)') + '.' + str(axesOut)
    
    ;write JPG
    if keyword_set(writeJPG) then begin
      if ( writeJPG eq '1' ) then begin
        write_jpeg, outputFilename + ".jpg", pic, true=1, quality=100
      endif else begin
        write_jpeg, writeJPG, pic, true=1, quality=100
      endelse
    endif
      
    ;write PNG
    if keyword_set(writePNG) then begin
      if ( writePNG eq '1' ) then begin
        write_png, outputFilename + ".png", pic
      endif else begin
        write_png, writePNG, pic
      endelse
    endif
    
    ;output PS - assume the ps device is already open
    if keyword_set(psOut) then begin
      tvimage,pic,true=1,position=[0.0,0.0,1.0,1.0],overplot=overPlot
    endif
    
    df = !NULL
  
  endforeach ;m

end

; contourGasSurfDens(): make a contour plot of the gas surface density

pro contourGasSurfDens, filePath=filePath, snapNum=snapNum, zoomSize=zoomSize, gridSize=gridSize

  ; config
  h = loadSnapshotHeader(filePath,snapNum=snapNum)
  
  if keyword_set(zoomSize) then $
    boxSize = [zoomSize,zoomSize,h.boxSize]
  if (not keyword_set(zoomSize)) then $
    boxSize = [h.boxSize,h.boxSize,h.boxSize] ;kpc
    
  boxCen  = [h.boxSize/2.0,h.boxSize/2.0,h.boxSize/2.0] ;kpc
  imgSize = [gridSize,gridSize]  ;px
  
  ; unused
  if 0 then begin

    axis0 = 0 ;x
    axis1 = 1 ;y
    mode  = 1 ;1=col mass, 2=mass-weighted quantity, 3=col density

    ; load positions and densities
    pos  = loadSnapshotSubset(filePath,snapNum=snapNum,partType='gas',field='pos')
    hsml = loadSnapshotSubset(filePath,snapNum=snapNum,partType='gas',field='hsml')
    mass = loadSnapshotSubset(filePath,snapNum=snapNum,partType='gas',field='mass')
  
    ; make spatial subset
    if keyword_set(zoomSize) then begin
      bc = boxCen[0]
      w = where(abs(pos[0,*]-bc) le zoomSize*1.5 and abs(pos[1,*]-bc) le zoomSize*1.5,count)
  
      print,'Found ['+str(count)+'] of ['+str(n_elements(mass))+'] particles inside zoomBox.'
      
      pos = pos[*,w]
      hsml = hsml[w]
      mass = mass[w]
    endif
  
    ; run sph kernel based density projection
    quant = fltarr(n_elements(mass))
    colMassMap = CalcSphMap(pos, hsml, mass, quant, imgSize=imgSize, boxSize=boxSize,$
                            boxCen=boxCen, axis0=axis0, axis1=axis1)
               
  endif
         
            
  ; contour and plot
  xy = findgen(gridSize)/gridSize * zoomSize + h.boxSize/2.0 - zoomSize/2.0

  levels = [1e-3,2e-3,3e-3,4e-3]
 
  fsc_contour,colMassMap,xy,xy,levels=levels,/overplot

end

; scatterPlotPos(): make a simple scatter plot of particle positions for one species

pro scatterPlotPos

  ; config
  basePath = '/n/home07/dnelson/sims.idealized/sims.gasSphere/'
  sP = { simPath  : basePath+'m11_rot0_cool0_wind0/output/' ,$
         plotPath : basePath ,$
         snap     : 0 }
  
  zoomSize = 1000.0 ;kpc
  partType = 'gas'

  axes = [0,2] ;x,z

  ; find number of snapshots and loop
  ;nSnaps = n_elements(file_search(filepath+'snap_*'))
  snapRange = [0,12,1]
  
  for snap=snapRange[0],snapRange[1],snapRange[2] do begin
    ; sizes
    sP.snap = snap
    h = loadSnapshotHeader(sP=sP)
      
    boxCen = h.boxSize/2.0 ;kpc
    
    ; load positions and densities
    pos  = loadSnapshotSubset(sP=sP,partType=partType,field='pos')

    ; make spatial subset
    if keyword_set(zoomSize) then begin
      w = where(abs(pos[axes[0],*]-boxCen) le zoomSize*1.0 and $
                abs(pos[axes[1],*]-boxCen) le zoomSize*1.0,count)
  
      print,'['+str(snap)+'] Found ['+str(count)+'] of ['+str((size(pos))[2])+'] particles inside zoomBox.'
      
      x = pos[axes[0],w]; - boxCen
      y = pos[axes[1],w]; - boxCen
    endif
    
    ; start plot
    start_PS, sP.plotPath + 'scatter_snap='+str(snap)+'_'+str(partType)+'.eps'
    
      xyrange = [boxCen-zoomSize/2.0,boxCen+zoomSize/2.0]
    
      cgPlot,[0],[0],/nodata,xrange=xyrange,yrange=xyrange,xstyle=1,ystyle=1,$
           xtitle="x [kpc]",ytitle="y [kpc]",aspect=1.0,$
           title="particle scatterplot - part ["+str(partType)+"] snap ["+str(snap)+"]"
    
      cgPlot,x,y,psym=3,/overplot
    
    end_PS, pngResize=60, /deletePS
  
  endfor ;i

end

; scatterPlotGasSphere(): make a 3x2 frame for GasSphere (gas and stars, different scales and directions)

pro scatterPlotGasSphere

  ; config
  run = 'm11_rot0.1_cool1_wind1'
  ;run = 'm11_rot0_cool0_wind0'
  
  sP = { simPath : '/n/home07/dnelson/sims.idealized/sims.gasSphere/'+run+'/output/' ,$
         plotPath : '/n/home07/dnelson/data5/frames.gasSphere/' ,$
         snap     : 0 }
  
  pConfigs = { p0 : { size:1000.0, pt:'gas',   axes:[0,2] } ,$
               p1 : { size:100.0,  pt:'gas',   axes:[0,2] } ,$
               p2 : { size:100.0,  pt:'gas',   axes:[0,1] } ,$
               p3 : { size:100.0,  pt:'stars', axes:[0,2] } ,$
               p4 : { size:40.0,   pt:'stars', axes:[0,2] } ,$
               p5 : { size:40.0,   pt:'stars', axes:[0,1] } }

  ; find number of snapshots and loop
  nSnaps = n_elements(file_search(sP.simPath+'snap_*'))
  snapRange = [0,nSnaps-1]
  print,snapRange
  
  for snap=snapRange[0],snapRange[1],1 do begin
    sP.snap = snap
   
    ; load
    h = loadSnapshotHeader(sP=sP)
    boxCen = h.boxSize/2.0 ;kpc
  
    gas_pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
    
    if h.nPartTot[partTypeNum('stars')] gt 0 then $
      stars_pos = loadSnapshotSubset(sP=sP,partType='stars',field='pos')
   
    print,'['+string(snap,format='(I3)')+'] t = '+string(h.time,format='(f5.3)')+$
      ' gas = '+string(h.nPartTot[partTypeNum('gas')],format='(I5)')+$
      ' stars = ' + string(h.nPartTot[partTypeNum('stars')],format='(I5)')
   
    ; already exists?
    if file_test(sP.plotPath + 'gasSphere_'+run+'_'+str(snap)+'.png') then continue
   
    ; start plot
    start_PS, sP.plotPath + 'gasSphere_'+run+'_'+str(snap)+'.eps', xs=16.0, ys=9.0
    
      pos = plot_pos(col=3,row=2,/gap)
      offset = [-0.06, -0.04, 0.02, 0.02]
    
      ; loop over each plot
      for m=0,n_tags(pConfigs)-1 do begin
        ; start plot
        config = pConfigs.(m)
        xyrange = [-config.size/2.0,config.size/2.0]
        
        cgPlot,[0],[0],/nodata,xrange=xyrange,yrange=xyrange,/xs,/ys,$
           xtitle="x [kpc]",ytitle="y [kpc]",aspect=1.0,pos=pos[m]+offset,/noerase
        
        cgText,(pos[m])[0]-0.06,(pos[m])[3]-0.02,$
          config.pt,color=cgColor('orange'),/normal
        
        if config.pt eq 'stars' and h.nPartTot[partTypeNum('stars')] eq 0 then continue
        
        ; spatial subset
        if config.pt eq 'gas'   then loc_pos = gas_pos
        if config.pt eq 'stars' then loc_pos = stars_pos
        
        w = where(abs(loc_pos[config.axes[0],*]-boxCen) le config.size and $
                  abs(loc_pos[config.axes[1],*]-boxCen) le config.size,count)
      
        xx = loc_pos[config.axes[0],w] - boxCen
        yy = loc_pos[config.axes[1],w] - boxCen
        
        ; plot points          
        cgPlot,xx,yy,psym=3,/overplot
           
      endfor ;m
      
      cgText,0.02,0.02,"t = "+string(h.time,format='(f5.3)'),$
        charsize=!p.charsize+0.5,color=cgColor('blue'),/normal
      
    end_PS, pngResize=40, /deletePS
    
  endfor ;snap

end

; plotMeshOverDensity(): plot voronoi over density for a sequence of snapshots (2D only)

pro plotMeshOverDensity

  workingPath = '/n/home07/dnelson/dev.tracer/'
  filePath    = workingPath + 'disk2d.test.VTS/output/'

  ; config
  meshBase   = filePath + "voronoi_mesh_"
  snapBase   = filePath + "snap_"
  
  boxSize    = [6.0,6.0]   ;x,y (code)
  xyScaleFac = 1.0
  
  nSnaps = n_elements(file_search(snapBase+"*"))
  
  for i=0,nSnaps-1 do begin
    print,i
    
    ; plots
    if (file_test(workingPath + 'mesh_'+string(i,format='(I04)')  +".png")) then begin
      print,' skip'
      continue
    endif
    
    ; load gas positions for plotting
    pos = loadSnapshotSubset(filePath,snapNum=i,partType='gas',field='pos')
    ids = loadSnapshotSubset(filePath,snapNum=i,partType='gas',field='ids')

    ;start_PS, workingPath + 'meshdens_'+string(i,format='(I04)')  +".eps", xs=4.0, ys=4.0
    ;  plotDensityField2D, filePath, i, /psOut, xyScaleFac=xyScaleFac
    ;  plotVoronoi2D, meshBase, boxSize, i, /oPlot, /fillWindow
    ;end_PS, pngResize=68, /deletePS
    
    start_PS, workingPath + 'mesh_'+string(i,format='(I04)')  +".eps", xs=4.0, ys=4.0
      plotVoronoi2D, meshBase, boxSize, i
    end_PS, pngResize=78, /deletePS
    
    start_PS, workingPath + 'meshzoom_'+string(i,format='(I04)')  +".eps", xs=4.0, ys=4.0
      plotVoronoi2D, meshBase, boxSize, i, zoom=10, pos=pos, ids=ids
    end_PS, pngResize=78, /deletePS
    
    start_PS, workingPath + 'dens_'+string(i,format='(I04)')  +".eps", xs=4.0, ys=4.0
      plotDensityField2D, filePath, i, /psOut, xyScaleFac=xyScaleFac
    end_PS, pngResize=78, /deletePS

  endfor

end
