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
                        xyScaleFac=xyScaleFac, minMax=minMax, log=log, overPlot=overPlot
  
  if (not keyword_set(writePNG) and not keyword_set(writeJPG) and not keyword_set(psOut)) then begin
    print,'ERROR: No output method specified.'
    stop
  endif
  
  ; if m is an array, treat as list of snapshot numbers and loop over each
  foreach m,snaps do begin

    ; load
    df = loadDensityField(filePath,m,axes=axes)
  
    ; log scaling of density
    if keyword_set(log) then begin
      w = where(df.dens gt 1e-8,comp=wc)
      df.dens[w] = alog10(df.dens[w])
      df.dens[wc] = min(df.dens[w])
    endif
  
    ; color tables and scaling
    loadct, 11, /silent ;11=blue-red
    loadct, 11, rgb_table=rgbTable
    tvlct, r, g, b, /get
    
    if not keyword_set(minMax) then $
      colorMinMax = minmax(df.dens)
    if keyword_set(minMax) then $
      colorMinMax = minMax
    
    print,' dens minmax:',colorMinMax
  
    colindex = (df.dens - colorMinMax[0])/(colorMinMax[1]-colorMinMax[0])*255.0
    
    ind = where(colindex ge 256.0)
    if ind(0) ne -1 then colindex(ind) = 255.9
    
    ind = where(colindex lt 0)
    if ind(0) ne -1 then colindex(ind) = 0
    
    colindex = byte(colindex)
    
    ;create image
    pic=bytarr(3,df.nPixelsXY[0],df.nPixelsXY[1])
  
    pic[0,*,*] = r[colindex]
    pic[1,*,*] = g[colindex]
    pic[2,*,*] = b[colindex]
    
    ;rescale image
    if keyword_set(xyScaleFac) then begin
      pic = rebin(pic,3,xyScaleFac*df.nPixelsXY[0],xyScaleFac*df.nPixelsXY[1])
    endif
  
    ;output filename
    if not keyword_set(axes) then axes = 0
    
    outputFilename = filePath + 'density_' + string(m,format='(i3.3)') + '.' + str(axes)
    
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
    print,'Running sphDensityProjection().'
  
    colMassMap = sphDensityProjection(pos, hsml, mass, imgSize=imgSize, boxSize=boxSize,$
                                      boxCen=boxCen, axis0=axis0, axis1=axis1, mode=mode, periodic=0)
               
  endif
         
            
  ; contour and plot
  xy = findgen(gridSize)/gridSize * zoomSize + h.boxSize/2.0 - zoomSize/2.0

  levels = [1e-3,2e-3,3e-3,4e-3]
 
  fsc_contour,colMassMap,xy,xy,levels=levels,/overplot

end

; scatterPlotPos(): make a simple scatter plot of particle positions for one species

pro scatterPlotPos

  ; config
  workingPath = '/n/home07/dnelson/dev.tracer/'
  filePath    = workingPath + 'col2Sph.gastr.1e5.f1/output/'
  
  zoomSize = 1000.0 ;kpc
  partType = 'gas'

  ; find number of snapshots and loop
  ;nSnaps = n_elements(file_search(filepath+'snap_*'))
  snapRange = [10,40,5]
  
  for snap=snapRange[0],snapRange[1],snapRange[2] do begin
    ; sizes
    h = loadSnapshotHeader(filePath,snapNum=snap)
      
    boxCen = h.boxSize/2.0 ;kpc
    
    ; load positions and densities
    pos  = loadSnapshotSubset(filePath,snapNum=snap,partType=partType,field='pos')
  
    ; make spatial subset
    if keyword_set(zoomSize) then begin
      w = where(abs(pos[0,*]-boxCen) le zoomSize*1.0 and abs(pos[1,*]-boxCen) le zoomSize*1.0,count)
  
      print,'['+str(snap)+'] Found ['+str(count)+'] of ['+str((size(pos))[2])+'] particles inside zoomBox.'
      
      x = pos[0,w]; - boxCen
      y = pos[1,w]; - boxCen
    endif
    
    ; start plot
    start_PS, workingPath + 'scatter_snap='+str(snap)+'_'+str(partType)+'.eps'
    
      xyrange = [boxCen-zoomSize/2.0,boxCen+zoomSize/2.0]
    
      fsc_plot,[0],[0],/nodata,xrange=xyrange,yrange=xyrange,xstyle=1,ystyle=1,$
           xtitle="x [kpc]",ytitle="y [kpc]",aspect=1.0,charsize=!p.charsize-1.0,$
           title="particle scatterplot - part ["+str(partType)+"] snap ["+str(snap)+"]"
    
      fsc_plot,x,y,psym=3,/overplot
    
    end_PS, pngResize=60, /deletePS
  
  endfor ;i

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
