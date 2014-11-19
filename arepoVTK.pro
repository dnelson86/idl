; arepoVTK.pro
; helper functions to test ArepoRT and ArepoVTK
; dnelson oct.2014

; combineStereoFrames(): combine two frames into side-by-side (and/or do 16->8 bit mapping)

pro combineStereoFrames, mmNum=mmNum

  if n_elements(mmNum) eq 0 then message,'Error: Specify mmNum'
  
  ; file config
  ; A[780,1399], B[1400,1999], C[2000,2599], D[2600,3199], E[3200,3975]
  ;mmNum = [0,3975] ;all
  inStr   = 'A0'
  inPath1 = '/n/home07/dnelson/ArepoVTK/run.subbox0/output/frames_8192/'
  inPath2 = '/n/home07/dnelson/ArepoVTK/run.subbox0/output/frames_360_'+inStr+'/'
  outPath = '/n/home07/dnelson/ArepoVTK/run.subbox0/output/frames_8k8bit/'
  
  fileBase = "frame_"
  fileEnd  = "_16bit.png"
  
  ; image config
  imMinMax = [0.0,0.65]
  imGamma  = 1.0 ; inverse of photoshop
  downto4k = 0
  
  for i=mmNum[0],mmNum[1] do begin
    ; load
    image1 = read_png(inPath1 + fileBase + string(i,format='(I04)') + fileEnd)
    
    if inStr ne 'A0' then begin
      image2 = read_png(inPath2 + fileBase + string(i,format='(I04)') + fileEnd)
      image_out = [[image1],[image2]] ; horizontal stack
    endif else begin
      image_out = image1 ; non-SbS
    endelse
    
    if downto4k eq 1 then image_out = rebin(image_out,3,4096,4096)
    
    ; scaling
    image_out /= 65535.0 ; normalize to [0,1]
    image_out = (image_out - imMinMax[0]) / (imMinMax[1]-imMinMax[0]) ; scale histogram
    image_out = image_out > 0.0 < 1.0 ; clamp
    image_out = image_out^imGamma ; gamma scaling
    image_out = round(image_out * 255.0) ; 8 bit
    image_out = byte(image_out) > 0B < 255B ; clamp
    
    ; write
    outFile = fileBase + string(i,format='(I04)') + ".png"
    write_png,outPath + outFile,image_out
    print,'Wrote: ['+outFile+']'
  endfor
  
  stop

end

; interpMissingFrames(): for corrupt subbox snapshots, create interpolated frames from the two neighbors

pro interpMissingFrames

  ; file config
  ;nums   = [636,1281,1961,2353,2806,2989,3026,3533] ;C1
  nums   = [949] ;[2353,2806,2989,3026,3533]
  inPath = '/n/home07/dnelson/ArepoVTK/run.subbox0/output/frames_8192/'
  
  fileBase = "frame_"
  fileEnd  = "_16bit.png"

  foreach num,nums do begin
    ; load
    image1 = read_png(inPath + fileBase + string(num-1,format='(I04)') + fileEnd)
    image2 = read_png(inPath + fileBase + string(num+1,format='(I04)') + fileEnd)
        
    ; mean combine
    image_out = 0.5 * ( float(image1) + float(image2) )
    image_out = uint(round(image_out)) ; 16 bit
    
    ; write
    outFile = fileBase + string(num,format='(I04)') + fileEnd
    write_png,inPath + outFile,image_out
    print,'Wrote: ['+outFile+']'
  endforeach
  
end

; checkMissingFrames(): check for missing numbered images in sequence

pro checkMissingFrames

  ; config
  mmNum = [0,3975]
  inPath = '/n/home07/dnelson/ArepoVTK/run.subbox0/output/frames_8192/'
  
  fileBase = "frame_"
  fileEnd  = "_16bit.png"
  
  for i=mmNum[0],mmNum[1] do begin
    fileName = inPath + fileBase + string(i,format='(I04)') + fileEnd
    if ~file_test(fileName) then print,fileName
  endfor

end

; plotSpoon3D(): quick point plot of a snapshot with a surface boundary

pro plotSpoon3D
 
  ; config
  axes = [0,2]
  path = '/n/home07/dnelson/sims.idealized/sims.coffee3D/output/'
  snaps = [0,5,10,15,20]
  colors = ['red','blue','green','orange','yellow']
  
  start_PS, path + 'spoon_pos.eps', xs=8, ys=8
  
    cgPlot,[0],[0],/nodata,xrange=[0,1],yrange=[0,1],/xs,/ys,xtitle="axis"+str(axes[0]),$
      ytitle="axis"+str(axes[1]),aspect=1.0
    
    times = []
  
    foreach snap,snaps,k do begin
      fpath = path + 'snap_' + string(snap,format='(I3.3)') + '.hdf5'
      
      x = h5_parse(fpath,/read)
      times = [times, x.header.time._data]
      
      ids = x.parttype0.particleids._data
      pos = x.parttype0.coordinates._data
      dens = x.parttype0.density._data
      mass = double(x.parttype0.masses._data)
      
      ; spoon boundingbox?
      ;print,'mass',total(mass[where(ids ge -1)],/double)
      
      w = where(ids eq -2,count)
      spoon_pos = pos[*,w]
      spoon_bbox = { x:minmax(spoon_pos[0,*]), y:minmax(spoon_pos[1,*]), z:minmax(spoon_pos[2,*]) }
      
      print,'snap = '+str(snap)+' gascount = '+str(count)
      print,' spoon bBox x ['+string(spoon_bbox.x[0],format='(f4.2)')+' '+$
                              string(spoon_bbox.x[1],format='(f4.2)')+']'+$
                        ' y ['+string(spoon_bbox.y[0],format='(f4.2)')+' '+$
                              string(spoon_bbox.y[1],format='(f4.2)')+']'+$
                        ' z ['+string(spoon_bbox.z[0],format='(f4.2)')+' '+$
                              string(spoon_bbox.z[1],format='(f4.2)')+']'
      
      cgPlot,spoon_pos[axes[0],*],spoon_pos[axes[1],*],psym=3,/overplot,color=colors[k]
      
      if k eq 1 then begin
        w = where(ids ge 0 and dens gt 1.5,count)
        print,'fluid',count
        cgPlot,pos[axes[0],w],pos[axes[1],w],psym=3,/overplot,color=cgColor('black')
      endif
    endforeach
    
    legend,'t = '+string(times,format='(f3.1)'),textcolor=colors,/top,/right
  
  end_PS
  
  stop

end

; makeVTKColorTable(): output .tbl discrete color table
; see: http://www.exelisvis.com/docs/LoadingDefaultColorTables.html

pro makeVTKColorTable

  start_PS,'dummy.eps' ; need DEVICE access somewhere
  
  ; load colortable into
  loadColorTable,'ncl/WhiteBlueGreenYellowRed'
  fileName = 'ncl_WhiteBlueGreenYellowRed.tbl'
  
  ; load RGB[3,N]
  tvlct,rgb,/get
  start = 0
  nVals = n_elements(rgb[*,0])

  ; alpha channel
  alpha = fltarr(nVals) + 1.0 ; uniform

  ; write
  openW,lun,fileName,/GET_LUN
  
  ; header
  printf,lun,'# comment'
  printf,lun,str(nVals-start)
  
  for j=start,nVals-1 do begin
    outText = string(rgb[j,0],format='(f5.1)') + " " + $
              string(rgb[j,1],format='(f5.1)') + " " + $
              string(rgb[j,2],format='(f5.1)') + " " + $
              string(alpha[j],format='(f4.2)')
    printf,lun,outText
  endfor
  
  close,lun
  free_lun, lun
    
  end_PS
  spawn,'rm dummy.eps'
    
  stop

end

; makeCosmoCutout(): output HDF5 format spatial cutout around a halo

pro makeVTKCosmoCutout

  ; config
  sP = simParams(res=512,run='tracer',redshift=2.0)
  
  haloID = 304 ;z2.304 z2.301 z2.130 z2.64
  gcID = getMatchedIDs(sPa=sP,sPg=sP,haloID=haloID)
  sizeFac = 3.7 ; times rvir for the bounding box of each cutout

  ; target list
  gc    = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sP)
  rvir  = gc.group_r_crit200[gc.subgroupGrNr[gcID.a]]
  
  ; get subhalo position and size of imaging box
  boxCen     = sgcen[*,gcID.a]
  boxSize    = ceil(sizeFac * rvir / 10.0) * 10.0
  boxSizeImg = [boxSize,boxSize,boxSize] ; cube
  
  ; exclude inner region (galaxy)
  excludeSize = ceil(0.25 * rvir / 10.0) * 10.0

  print,'boxCen',boxCen
  print,'boxSize',boxSize
  print,'excludeSize',excludeSize

  ; make cutout
  createSnapshotCutout,sP=sP,fOut='/n/home07/dnelson/cutout_tracer2.hdf5',$
    cenPos=boxCen,boxSize=boxSizeImg,excludeSize=excludeSize,/includeGas,/convertUtoTemp

end

pro plotScalings

  ; config
  workingPath = "/n/home07/dnelson/vis/ArepoVTK/"
  
  ; strong (const problem size) data (128^3 20mpc)
  strongNumCores        =       [1,     2,     3,     4,     5,     6,     7,     8,     8,    8  ]
  strongNumThreads      = float([8,     16,    24,    32,    40,    48,    56,    64,    128,  256])
  strongCPUTimesInit    =       [54.5,  62.8,  54.6,  54.6,  54.5,  54.6,  54.6,  54.6,  54.7, 54.5]  ;sec
  strongCPUTimesRender  =       [748.0, 451.3, 261.4, 196.5, 163.9, 136.3, 118.1, 102.4, 97.2, 95.4]  ;sec
  strongWallTimesLSF    = float([1117,  1302,  1118,  1115,  1124,  1123,  1127,  1124,  1123, 1118]) ;sec
  
  strongSpeedup    = strongCPUTimesRender[0]/strongCPUTimesRender
  strongEfficiency = strongSpeedup / strongNumCores
  
  ; plot
  start_PS, workingPath + "scaling.strong.eps"
    fsc_plot,[0],[0],/nodata,xrange=[0,9],yrange=[0,9.0],/xs,ys=9,$
             xtitle="Number of Cores",ytitle="Speedup ("+textoidl("T_S/T_P")+")",$
             title="ArepoRT Strong Scaling (pthreads 128^3 20Mpc)",$
             position=[0.15,0.15,0.85,0.85],charsize=!p.charsize-0.5,$
             xticks=[7],xtickv=[1,2,3,4,5,6,7,8]
             
    ; speedup
    fsc_plot, strongNumCores, strongSpeedup, psym=4, /overplot
    fsc_plot,[0.5,8.5],[0.5,8.5],line=1,/overplot
    
    ; efficiency
    ;fsc_plot, strongNumCores, strongEfficiency, psym=2, /overplot
    ;fsc_plot,[1.5,8.5],[1.0,1.0],line=0,thick=!p.thick-0.5,/overplot
    
    ; actual wallclock times for 128^3 20mpc
    fsc_axis,9.0,0.0,0.0,/yaxis,yrange=[850.0,1.0],/ys,charsize=!p.charsize-0.5
    fsc_text,0.92,0.5,"Render Time [sec]",alignment=0.5,orientation=270.0,/normal,charsize=!p.charsize-0.5
  end_PS
  
  
  stop
end


pro makeArepoICs

  ; config
  boxSize = 1.0

  fileName = "/n/home07/dnelson/ArepoVTK/test/Arepo3a.hdf5"
  
  rho_base = 1.1
  u_base   = 5.0
  
  Nx = 4L
  Ny = 4L
  Nz = 4L

  Ntot = Nx*Ny*Nz+1 ;+1 for central point test
  
  Lx = boxSize
  Ly = boxSize
  Lz = boxSize
  
  ;create arrays
  pos  = fltarr(3,Ntot)
  vel  = fltarr(3,Ntot)
  rho  = fltarr(Ntot)
  u    = fltarr(Ntot)
  id   = lindgen(Ntot)+1L
  
  deltaXYZ = [Lx/Nx,Ly/Ny,Lz/Nz]
  
  pid = 0
  
  ;fill
  for i=0L,Nx-1 do begin
    for j=0L,Ny-1 do begin
      for k=0L,Nz-1 do begin
        ;pid = i*(Nx*Ny) + j*Nx + k ;only works for Nx=Ny=Nz
  
        pos[0,pid] = i*deltaXYZ[0]+deltaXYZ[0]/2.0 ;x
        pos[1,pid] = j*deltaXYZ[1]+deltaXYZ[1]/2.0 ;y
        pos[2,pid] = k*deltaXYZ[2]+deltaXYZ[2]/2.0 ;z
 
        ;all points
        rho[pid] = rho_base
        u[pid]   = u_base
        
        print,pid,pos[0,pid],pos[1,pid],pos[2,pid],rho[pid]
        pid++
      endfor
    endfor
  endfor
  
  ; central point
  pos[0,pid] = 0.5
  pos[1,pid] = 0.5
  pos[2,pid] = 0.5
  rho[pid]   = 1.2 * rho_base
  u[pid]     = 2.0 * u_base
  
  print,pid,pos[0,pid],pos[1,pid],pos[2,pid],rho[pid]

  ; Arepo3c: perturb outer layer so we don't have a cubical grid
  w = where( (pos[0,*] lt 0.2 or pos[0,*] gt 0.8) or $
             (pos[1,*] lt 0.2 or pos[1,*] gt 0.8) or $
             (pos[2,*] lt 0.2 or pos[2,*] gt 0.8), count)

  print,'Found ['+str(count)+'] on outer layer.'
  
  seed = 424242L
  
  ;for i=0,2 do begin
  ;  pp = round(5*(randomu(seed,count)-0.5))/50.0 ; -0.04 to +0.04
  ;  pos[i,w] += pp
  ;endfor
  
  ;for i=0,pid-1 do print,i,pos[0,i],pos[1,i],pos[2,i]
  ; end Arepo3c

  ;write
  ;writeICFile,fileName,pos,vel,id,rho,u
  writeICFileHDF5,fileName,boxSize,pos,vel,id,rho,u
  
end

pro makeFloatICs

  boxsize = 1.0

  ; config
  fileName = "~/vis/ArepoVTK/test/test2.txt"
  
  nx = 101
  ny = 101
  nz = 101
  
  dens_base  = 0.2
  dens_over  = 2.0
  dens_under = 0.0
  
  prism  = [0.6,0.9,0.6,0.8,0.0,1.0] ;xmin,xmax,ymin,ymax,zmin,zmax
  sphere = [0.3,0.3,0.3,0.2] ;xcen,ycen,zcen,radius
  
  ntot = long(nx)*ny*nz
  dens = fltarr(nx,ny,nz)
  
  ; calculate point locations
  xStep = boxSize / (nx-1)
  yStep = boxSize / (ny-1)
  zStep = boxSize / (nz-1)

  xPts = findgen(nx) * xStep
  yPts = findgen(ny) * yStep
  zPts = findgen(nz) * zStep

  ; make density structure
  
    ; base constant
    dens += dens_base
    
    for x=0,nx-1 do begin
      for y=0,ny-1 do begin
        for z=0,nz-1 do begin
        
        ; rectangular prism (brick) overdensity
        if(xPts[x] ge prism[0] and xPts[x] lt prism[1] and $
           yPts[y] ge prism[2] and yPts[y] lt prism[3] and $
           zPts[z] ge prism[4] and zPts[z] lt prism[5]) then begin
           
           dens[x,y,z] = dens_over
        endif
        
        ; spherical underdensity
        rdist = sqrt( (xPts[x] - sphere[0])^2.0 + (yPts[y] - sphere[1])^2.0 + (zPts[z] - sphere[2])^2.0 )
        if(rdist le sphere[3]) then begin
          dens[x,y,z] = dens_under
        endif
        
        endfor
      endfor
    endfor

  ; write
  openW, lun, fileName, /GET_LUN
  
    ; header
    printf, lun, str(ntot)
    printf, lun, str(nx)
    printf, lun, str(ny)
    printf, lun, str(nz)
    
    ; data
    ;    x = Clamp(x, 0, nx-1);
    ;    y = Clamp(y, 0, ny-1);
    ;    z = Clamp(z, 0, nz-1);
    ;    return density[z*nx*ny + y*nx + x];
    for z=0,nz-1 do begin
      for y=0,ny-1 do begin
        lineStr = ""
        for x=0,nx-1 do begin
          lineStr = lineStr + string(dens[x,y,z],format='(f4.2)') + " "
        endfor
        printf, lun, strmid(lineStr,0,strlen(lineStr)-1)
      endfor
      ;printf, lun, ""
    endfor
  
  ;close handle
  free_lun, lun

  print,'wrote: ',ntot,' pts to: ',fileName
end
