; arepoVis2D.pro
; dnelson
; 7/5/10
;
; visualization for 2D arepo runs

@helper
@arepoLoad

; plotVoronoi2D(): plot "voronoi_mesh"
;

pro plotVoronoi2D, fBase, boxSize, i, oPlot=oPlot, fillWindow=fillWindow

  ;config
  nDims   = 2
  color   = fsc_color('dark gray')
  
  ;load
  nLoad = loadVoronoi2D(fBase,nDims,i,nEdges,nEdgesOffset,edgeList,xyzEdges)
  nGas = nLoad[0]

  ;plot
  if keyword_set(fillWindow) then $
    plot,[0],[0],/nodata,xrange=[0,boxSize[0]],yrange=[0,boxSize[1]], $
         xstyle=5,ystyle=5,position=[0.0,0.0,1.0,1.0],noerase=oPlot ;5=1+4 bitwise, exact+supressed
  if not keyword_set(fillWindow) then $
    plot,[0],[0],/nodata,xrange=[0,boxSize[0]],yrange=[0,boxSize[1]], $
         xstyle=1,ystyle=1,noerase=oPlot

  for j=0L, nGas-1 do begin
      x = transpose(xyzEdges[0, edgeList[nEdgesOffset[j]:nEdgesOffset[j]+nEdges[j]-1]])
      y = transpose(xyzEdges[1, edgeList[nEdgesOffset[j]:nEdgesOffset[j]+nEdges[j]-1]])

      plots,[x,x(0)],[y,y(0)],noclip=0,color=color
   endfor

end

; plotDensityField2D(): plot "density_field_"
;

pro plotDensityField2D, fBase, i, writePNG=writePNG, writeJPG=writeJPG, $
                        tvOut=tvOut, psOut=psOut, xyScaleFac=xyScaleFac, $
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
  
  print,'minmax:',colorMinMax

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
  
  ;output to window
  if keyword_set(tvOut) then begin
    if ( !d.x_size ne (size(pic))[2] ) or ( !d.y_size ne (size(pic))[3] ) or ( !d.window eq -1 ) then $
      window,1,xs=(size(pic))[2],ys=(size(pic))[3] ;size(pic) = [ndims,dim1,dim2,dim3,type,nelements]
    tv,pic,true=1
  endif
  
  ;output PS
  if keyword_set(psOut) then begin
    ;assume the ps device is already open
    tvimage,pic,true=1,position=[0.0,0.0,1.0,1.0]
  endif

end

; makeDensityField2D(): no claims this does anything remotely correct
;                       just puts particles in pixel boxes, no SPH smoothing/density recon

function makeDensityField2D, snapBase, i, boxSize, nPixelsXY

  stepX = boxSize[0] / nPixelsXY[0]
  stepY = boxSize[1] / nPixelsXY[1]

  ; load snapshot
  nParticles = loadSnapshotHDF5(snapBase,i,time,pos,vel,id,mass,u,rho)
  
  ; create density grid
  dens = fltarr(nPixelsXY[0],nPixelsXY[1])
  
  for i=0,nPixelsXY[0]-1 do begin
    for j=0,nPixelsXY[1]-1 do begin
      minX = i*stepX
      maxX = minX + stepX
      minY = j*stepY
      maxY = minY + stepY
      
      w = where( (pos[0,*] ge minX) and (pos[0,*] lt maxX) and $
                 (pos[1,*] ge minY) and (pos[1,*] lt maxY), count)
                 
      if (count ne 0) then begin
        rhoAvg = mean(rho[w])
      endif else begin
        rhoAvg = 0.0
      endelse
        
      dens[i,j] = rhoAvg
    
    endfor
  endfor

  return, dens
  
end

;plot voronoi over density
pro plotMeshOverDensity

  filePath    = '/n/home07/dnelson/ArepoTests/KH_Instability_45/output/'
  workingPath = '/n/home07/dnelson/ArepoTests/KH_Instability_45/output2/'

  ; config
  meshBase   = filePath+"voronoi_mesh_"
  densBase   = filePath+"density_field_"
  snapBase   = filePath+"snap_"
  
  boxSize    = [1.0,1.0]   ;x,y (code)
  xyScaleFac = 1.0
            
  ; get global scaling factor for movies
  colorMinMax = [-2.0,5.0]
  
  nSnaps = n_elements(file_search(densBase+"*"))
  
  for i=0,nSnaps-1 do begin
    
    print,i
    
    fileBase = workingPath+'meshdens_'+string(i,format='(I04)')  
  
    ; plot
    PS_Start, FILENAME=fileBase+".eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
              /encapsulated,decomposed=0, xs=4.0, ys=4.0, /inches
  
      plotDensityField2D, densBase, i, /psOut, xyScaleFac=xyScaleFac
      plotVoronoi2D, meshBase, boxSize, i, /oPlot, /fillWindow
     
    PS_End, /PNG, /Delete_PS, Resize=68 ;PNG size=[xs,ys]*300*(resize/100)

  endfor

end
