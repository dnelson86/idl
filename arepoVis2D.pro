; arepoVis2D.pro
; dnelson
; dec 2011
;
; visualization for 2D arepo runs

@helper
@arepoLoad

; plotVoronoi2D(): plot "voronoi_mesh"
;

pro plotVoronoi2D, fBase, boxSize, i, oPlot=oPlot, fillWindow=fillWindow, zoom=zoom

  ;config
  nDims   = 2
  color   = fsc_color('dark gray')
  thick   = !p.thick-1.5
  
  ;load
  nLoad = loadVoronoi2D(fBase,nDims,i,nEdges,nEdgesOffset,edgeList,xyzEdges)
  nGas = nLoad[0]

  ; set ranges
  if (not keyword_set(zoom)) then begin
    xrange = [0,boxSize[0]]
    yrange = [0,boxSize[1]]
  endif else begin
    xrange = [boxSize[0]/2.0-boxSize[0]/zoom,boxSize[0]/2.0+boxSize[0]/zoom]
    yrange = [boxSize[1]/2.0-boxSize[1]/zoom,boxSize[1]/2.0+boxSize[1]/zoom]
  endelse

  ;plot
  if keyword_set(fillWindow) then $
    plot,[0],[0],/nodata,xrange=xrange,yrange=yrange, $
         xstyle=5,ystyle=5,position=[0.0,0.0,1.0,1.0],noerase=oPlot ;5=1+4 bitwise, exact+supressed
  if not keyword_set(fillWindow) then $
    plot,[0],[0],/nodata,xrange=xrange,yrange=yrange, $
         xstyle=1,ystyle=1,position=[0.05,0.05,0.95,0.95],charsize=!p.charsize-1.0,noerase=oPlot

  for j=0L, nGas-1 do begin
      x = transpose(xyzEdges[0, edgeList[nEdgesOffset[j]:nEdgesOffset[j]+nEdges[j]-1]])
      y = transpose(xyzEdges[1, edgeList[nEdgesOffset[j]:nEdgesOffset[j]+nEdges[j]-1]])

      plots,[x,x(0)],[y,y(0)],noclip=0,color=color,thick=thick
   endfor

end

; plotDensityField2D(): plot "density_field_"
;

pro plotDensityField2D, fBase, i, writePNG=writePNG, writeJPG=writeJPG, $
                        psOut=psOut, xyScaleFac=xyScaleFac, $
                        dens=dens, colorMinMax=colorMinMax, plotHist=plotHist
  
  ;load
  if not keyword_set(dens) then begin
    nPixelsXY = loadDensityField2D(fBase,i,dens)
  endif else begin
    nPixelsXY = [ (size(dens))[1],(size(dens))[2] ]
  endelse

  ;color tables and scaling
  loadct, 11, /silent
  loadct, 11, rgb_table=rgbTable
  tvlct, r, g, b, /get
  
  if not keyword_set(colorMinMax) then begin
    colorMinMax = minmax(dens)
  endif
  
  print,' dens minmax:',colorMinMax

  if keyword_set(plotHist) then begin
    ext = string(i,format='(i3.3)')
    f = fBase + ext + ".eps"
    w = where(dens gt colorMinMax[0] and dens lt colorMinMax[1])
    start_PS,f
      plothist,dens[w],/auto
    end_PS
  endif

  colindex= (dens - colorMinMax[0])/(colorMinMax[1]-colorMinMax[0])*255.0
  
  ind = where(colindex ge 256.0)
  if ind(0) ne -1 then colindex(ind) = 255.9
  
  ind = where(colindex lt 0)
  if ind(0) ne -1 then colindex(ind) = 0
  
  colindex = byte(colindex)
  
  ;create image
  pic=bytarr(3,nPixelsXY[0],nPixelsXY[1])

  pic[0,*,*] = r[colindex]
  pic[1,*,*] = g[colindex]
  pic[2,*,*] = b[colindex]
  
  ;rescale image
  if keyword_set(xyScaleFac) then begin
    pic = rebin(pic,3,xyScaleFac*nPixelsXY[0],xyScaleFac*nPixelsXY[1])
  endif

  ;write JPG
  if keyword_set(writeJPG) then begin
    if ( writeJPG eq '1' ) then begin
      ext = string(i,format='(i3.3)')
      fName = fBase + ext + ".jpg"
      
      write_jpeg, fName, pic, true=1, quality=100
    endif else begin
      write_jpeg, writeJPG, pic, true=1, quality=100
    endelse
  endif
    
  ;write PNG
  if keyword_set(writePNG) then begin
    if ( writePNG eq '1' ) then begin
      ext = string(i,format='(i3.3)')
      fName = fBase + ext + ".png"
    
      write_png, fName, pic
      ;write_png, ext+".png", rebin(pic,3,3*nPixelsXY[0],3*nPixelsXY[1])
      ;write_png, fName+".2.png", congrid(pic,3,4*nPixelsXY[0],4*nPixelsXY[1],cubic=-0.5)
      ;write_png, fName+".3.png", congrid(pic,3,4*nPixelsXY[0],4*nPixelsXY[1])
      ;write_png, fName+".4.png", rebin(pic,3,4*nPixelsXY[0],4*nPixelsXY[1], /sample)
    endif else begin
      write_png, writePNG, pic
    endelse
  endif
  
  ;output PS
  if keyword_set(psOut) then begin
    ;assume the ps device is already open
    tvimage,pic,true=1,position=[0.0,0.0,1.0,1.0]
  endif

end

;plot voronoi over density
pro plotMeshOverDensity

  workingPath = '/n/home07/dnelson/dev.tracer/'
  filePath    = workingPath + 'disk2d.128/output/'

  ; config
  meshBase   = filePath + "voronoi_mesh_"
  densBase   = filePath + "density_field_"
  snapBase   = filePath + "snap_"
  
  boxSize    = [6.0,6.0]   ;x,y (code)
  xyScaleFac = 1.0
  
  nSnaps = n_elements(file_search(densBase+"*"))
  
  for i=0,nSnaps-1 do begin
    print,i
    
    ; plots
    if (file_test(workingPath + 'mesh_'+string(i,format='(I04)')  +".png")) then begin
      print,' skip'
      ;continue
    endif
    
    ;start_PS, workingPath + 'meshdens_'+string(i,format='(I04)')  +".eps", xs=4.0, ys=4.0
    ;  plotDensityField2D, densBase, i, /psOut, xyScaleFac=xyScaleFac
    ;  plotVoronoi2D, meshBase, boxSize, i, /oPlot, /fillWindow
    ;end_PS, pngResize=68, /deletePS
    
    ;start_PS, workingPath + 'mesh_'+string(i,format='(I04)')  +".eps", xs=4.0, ys=4.0
    ;  plotVoronoi2D, meshBase, boxSize, i
    ;end_PS, pngResize=78, /deletePS
    
    ;start_PS, workingPath + 'meshzoom_'+string(i,format='(I04)')  +".eps", xs=4.0, ys=4.0
    ;  plotVoronoi2D, meshBase, boxSize, i, zoom=10
    ;end_PS, pngResize=78, /deletePS
    
    start_PS, workingPath + 'dens_'+string(i,format='(I04)')  +".eps", xs=4.0, ys=4.0
      plotDensityField2D, densBase, i, /psOut, xyScaleFac=xyScaleFac
    end_PS, pngResize=78, /deletePS

  endfor

end
