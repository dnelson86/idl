; arepoVTK.pro
; helper functions to test ArepoRT and ArepoVTK
; dnelson aug.2012

@helper

; makeCosmoCutout(): output HDF5 format spatial cutout around a halo

pro makeVTKCosmoCutout

  ; config
  sP = simParams(res=512,run='tracer',redshift=2.0)
  
  gcID    = 1996 ; NEW512 z2.304 g2342 a1996
  sizeFac = 3.6 ; times rvir for the bounding box of each cutout

  ; target list
  gc    = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sP) 
  
  ; get subhalo position and size of imaging box
  boxCen     = sgcen[*,gcID]
  boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcID]] / 10.0) * 10.0
  boxSizeImg = [boxSize,boxSize,boxSize] ; cube

  print,'boxCen',boxCen
  print,'boxSize',boxSize

  ; make cutout
  createSnapshotCutout,sP=sP,fOut='cutout.hdf5',$
    cenPos=boxCen,boxSize=boxSizeImg,/includeGas,/convertUtoTemp,/verbose

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

; genOrbitConfigs(): generate series of configuration files derived from a template that move the
;                    camera position in a circular orbit about a center (x,y,z) position over a 
;                    specified number of frames, tangent to some axis

pro genOrbitConfigs

  ; orbit config
  xyzCen    = [1123.2,7568.8,16144.2] ;kpc/h in box
  radius    = 500.0 ;kpc

  numFrames = 360 ;360
  tanAxis   = 'y'
  
  ; path config
  path         = '/n/home07/dnelson/vis/ArepoVTK/orbit_test/'
  templateName = path+'config_template.txt'
  filenameBase = 'orbit_test_'
  
  ; generate (x,y,z) camera positions
  cameraPos = fltarr(numFrames,3)
  
  for i=0,numFrames-1 do begin
    if (tanAxis eq 'x') then begin
      cameraPos[i,0] = xyzCen[0]
      cameraPos[i,1] = xyzCen[1] + radius * cos(float(i)/numFrames * 2*!pi)
      cameraPos[i,2] = xyzCen[2] + radius * sin(float(i)/numFrames * 2*!pi)
    endif
    if (tanAxis eq 'y') then begin
      cameraPos[i,0] = xyzCen[0] + radius * cos(float(i)/numFrames * 2*!pi)
      cameraPos[i,1] = xyzCen[1]
      cameraPos[i,2] = xyzCen[2] + radius * sin(float(i)/numFrames * 2*!pi)
    endif
    if (tanAxis eq 'z') then begin
      cameraPos[i,0] = xyzCen[0] + radius * cos(float(i)/numFrames * 2*!pi)
      cameraPos[i,1] = xyzCen[1] + radius * sin(float(i)/numFrames * 2*!pi)
      cameraPos[i,2] = xyzCen[2]
    endif
  endfor
  
  ; load template
  nRows = file_lines(templateName)
  fileText = strarr(nRows)
  
  openR, lun, templateName, /GET_LUN
  
  for i=0,nRows-1 do begin
    tt = ''
    readF,lun,tt
    fileText[i] = tt
  endfor
  
  close, lun
  free_lun, lun
  
  ; write new config files
  for i=0,numFrames-1 do begin
    ; copy template text
    outText = fileText
    
    ; make camera position string
    cameraPosStr = str(string(cameraPos[i,0],format='(f9.3)')) + " " + $
                   str(string(cameraPos[i,1],format='(f9.3)')) + " " + $
                   str(string(cameraPos[i,2],format='(f9.3)'))
    
    ; replace FRAMENAME and CAMERAPOSTUPLE
    strreplace,outText,'FRAMENAME','frame_'+str(i)+'.tga'
    strreplace,outText,'CAMERAPOSTUPLE',cameraPosStr
    
    ; write
    writeName = path + 'config_'+str(i)+'.txt'
    
    openW,lun,writeName,/GET_LUN
    
    for j=0,nRows-1 do begin
      printf,lun,outText[j]
    endfor
    
    close,lun
    free_lun, lun
    
    print,writeName + " -- " + cameraPosStr
  endfor
end

pro makeArepoICs

  ; config
  boxSize = 1.0

  fileName = "/n/home07/dnelson/vis/ArepoVTK/test/Arepo2b.hdf5"
  
  rho_base = 1.1
  u_base   = 5.0
  
  Nx = 2L
  Ny = 2L
  Nz = 2L

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
  
        ;xyz = [pos[0,pid], pos[1,pid], pos[2,pid]]
  
        ;all points
        rho[pid] = rho_base
        u[pid]   = u_base
        
        ;small perturb
        ;if (pid eq 10) then begin
        ;  rho[pid] *= 2.0
        ;endif
        
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
