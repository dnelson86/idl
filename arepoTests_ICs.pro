; arepoTests_ICs.pro
; dnelson
; 03.2011

@helper

; 2D arepo.cuda input
pro gen_arepo_cuda_2D_input

  ; config
  fOut = 'input.txt'
  
  Lx = 1.0
  Ly = 1.0
  
  N = 1024L*1024L
  
  ; setup arrays
  pos = fltarr(3,N)
  
  ; uniformly random positions
  pos[0,*] = randomu(seed,N)
  pos[1,*] = randomu(seed,N)
  
  visFlag = 0
  
  if (visFlag eq 1) then begin
    ; move some points along vertical edges
    numAdd = round(sqrt(N) / 2.0)
    
    ; set x-pos
    pos[0,0:numAdd-1] = 0.0
    pos[0,numAdd:2*numAdd-1] = 1.0
    
    ; set y-pos
    pos[1,0:numAdd-1] = findgen(numAdd) / (numAdd-1) * 0.90
    pos[1,numAdd:2*numAdd-1] = findgen(numAdd) / (numAdd-1) * 0.90
    
    ; compress to y<0.5
    pos[1,*] /= 2.0
  endif

  ; write
  openw,lun,fOut,/GET_LUN
  
  printf,lun,"#"
  printf,lun,"# test input file with "+str(N)+" points in [0,1]x[0,1]"
  printf,lun,"# gen_arepo_cuda_2D_input"
  printf,lun,"#"
  printf,lun,str(N)
  
  for i=0,N-1 do begin
    printf,lun,string(pos[0,i],format='(f12.10)')+" "+string(pos[1,i],format='(f12.10)')
  endfor

  close,lun
  free_lun,lun
  
  print,'wrote ',fOut
  
end

; 2D Kelvin-Helmholtz Instability ICs (only really tested for angle=0)
; BoxSize=1.0 PERIODIC TWODIMS GAMMA=1.4

pro gen_KH_2D_ICs

  ; config
  angle = 0.0 ;CCW from 0 wrt positive x-axis
  fOut  = 'KH_Instability_0.dat'
  
  stripeWidth = 0.5
  
  ga   = 1.4
  P0   = 2.5
  rho1 = 1.0
  rho2 = 2.0
  
  v     = 0.5
  sigma = 0.05/sqrt(2.0)
  w0    = 0.1
  
  Lx = 1.0
  Ly = 1.0
  
  Nx = 64L
  Ny = 64L
  
  ; setup arrays
  pos  = fltarr(3,Nx*Ny)
  vel  = fltarr(3,Nx*Ny)
  dens = fltarr(Nx*Ny)
  u    = fltarr(Nx*Ny)
  id   = lindgen(Nx*Ny)+1L
  
  deltax = Lx/Nx
  deltay = Ly/Ny

  ;angle calc
  x2 = Lx
  y2 = tan(angle*!dtor) * Lx/2.0 + Ly/2.0
  x1 = 0.0
  y1 = Ly/2.0 - tan(angle*!dtor) * Lx/2.0
  
  for i=0L,Nx-1 do begin
   for j=0L,Ny-1 do begin
    pos[0,i+j*Nx] = i*deltax+deltax/2.0  
    pos[1,i+j*Nx] = j*deltay+deltay/2.0    
    pos[2,i+j*Nx] = 0.0
  
    x = pos[0,i+j*Nx]
    y = pos[1,i+j*Nx]
    
    d = abs( (x2-x1)*(y1-y) - (x1-x)*(y2-y1) ) / sqrt( (x2-x1)^2.0 + (y2-y1)^2.0 )

    ;top/bottom
    if (d ge stripeWidth/2.0) then begin   
     dens[i+j*Nx]  = rho1
     P             = P0
     u[i+j*Nx]     = P/rho1/(ga-1.0) 
     vel[0,i+j*Nx] = -v*cos(angle*!dtor)
     vel[1,i+j*Nx] = -v*sin(angle*!dtor)
    endif 
     
    ;middle stripe
    if (d lt stripeWidth/2.0) then begin   
     dens[i+j*Nx]  = rho2
     P             = P0
     u[i+j*Nx]     = P/rho2/(ga-1.0)
     vel[0,i+j*Nx] = +v*cos(angle*!dtor)
     vel[1,i+j*Nx] = +v*sin(angle*!dtor)
    endif 
  
    ;initial eigen perturbation
    vel[0,i+j*Nx] = w0 * sin(4.0*!pi*x*cos(angle*!dtor)) * $
      (exp(-(d+stripeWidth/2.0)^2.0/(2.0*sigma^2.0)) + $
       exp(-(d-stripeWidth/2.0)^2.0/(2.0*sigma^2.0))) * sin(angle*!dtor)
    vel[1,i+j*Nx] = w0 * sin(4.0*!pi*x*cos(angle*!dtor)) * $
      (exp(-(d+stripeWidth/2.0)^2.0/(2.0*sigma^2.0)) + $
       exp(-(d-stripeWidth/2.0)^2.0/(2.0*sigma^2.0))) * cos(angle*!dtor)
  
   endfor
  endfor

  ;write
  gas    = {pos:pos,vel:vel,id:id,mass:mass,u:u}
  
  writeICFile,fOut,part0=gas

end

; 3D Taylor-Sedov (point explosion) blastwave IC
; BoxSize=1.0 PERIODIC GAMMA=1.4

pro gen_TS_BlastWave_3D_ICs

  ;config
  fOut = "TS_BlastWave_3D.dat"
  
  ga      = 1.4
  P1      = 1.0
  rho1    = 1.0
  uCenter = 1e4
  
  Lx = 1.0
  Ly = 1.0
  Lz = 1.0
  
  Nx = 45L
  Ny = 45L
  Nz = 45L
  
  ;create arrays
  pos  = fltarr(3,Nx*Ny*Nz)
  vel  = fltarr(3,Nx*Ny*Nz)
  mass = fltarr(Nx*Ny*Nz)
  u    = fltarr(Nx*Ny*Nz)
  id   = lindgen(Nx*Ny*Nz)+1L
  
  deltaXYZ = [Lx/Nx,Ly/Ny,Lz/Nz]
  
  ;fill
  for i=0L,Nx-1 do begin
    for j=0L,Ny-1 do begin
      for k=0L,Nz-1 do begin
        pid = i*(Nx*Ny) + j*Nx + k ;only works for Nx=Ny=Nz?
  
        pos[0,pid] = i*deltaXYZ[0]+deltaXYZ[0]/2.0 ;x
        pos[1,pid] = j*deltaXYZ[1]+deltaXYZ[1]/2.0 ;y
        pos[2,pid] = k*deltaXYZ[2]+deltaXYZ[2]/2.0 ;z
  
        ;xyz = [pos[0,pid], pos[1,pid], pos[2,pid]]
  
        ;all points
        mass[pid] = rho1*(Lx*Ly*Lz)/(Nx*Ny*Nz) ;only for Cartesian (equal vol cells)
        u[pid]    = P1/rho1/(ga-1.0)
  
        ;center point
        if ( (i eq floor(Nx/2.0)) and (j eq floor(Ny/2.0)) and (k eq floor(Nz/2.0)) ) then begin
          u[pid]=uCenter
      endif
  
      endfor
    endfor
  endfor

  ;write
  gas    = {pos:pos,vel:vel,id:id,mass:mass,u:u}
  
  writeICFile,fOut,part0=gas
  
end

; 3D shocktube with tracer particles
; BoxSize=2.0 PERIODIC LONG_X=10.0 GAMMA=1.4

pro gen_shocktube_tracer_3D_ICs

  ; config
  fOut = "ics.dat"
  
  ga   = 1.4
  
  P1   = 1.0
  rho1 = 1.0
  P2   = 0.1795
  rho2 = 0.25
  
  Lx = 20.0
  Ly = 2.0
  Lz = 2.0
  
  Nx   = 100L
  Ny   = 10L
  Nz   = 10L
  
  tfac = 1L
  
  ; deriv
  if (( round(Nx*Ny*Nz/4) ne (Nx*Ny*Nz/4) ) or ( round(tfac^(1.0/3.0)) ne (tfac^(1.0/3.0)) ) ) then begin
    print,'Error: Choose nicer numbers.'
    return
  endif
  
  nTotGas = Nx*Ny*Nz
  
  deltax = Lx/Nx
  deltay = Ly/Ny
  deltaz = Lz/Nz
  
  ; create gas arrays
  id   = lindgen(nTotGas)+1L
  pos  = fltarr(3,nTotGas)
  vel  = fltarr(3,nTotGas)
  mass = fltarr(nTotGas)
  u    = fltarr(nTotGas)
  
  ; set gas properties
  for i=0L,Nz-1 do begin
    for j=0L,Ny-1 do begin
      for k=0L,Nx-1 do begin
        pid = i*(Nx*Ny) + j*Nx + k

        pos[0,pid] = k*deltax+deltax/2.0
        pos[1,pid] = j*deltay+deltay/2.0
        pos[2,pid] = i*deltaz+deltaz/2.0
    
        ;left  
        if (pos[0,pid] lt Lx/2.0) then begin
          mass[pid] = rho1*(Lx/2.0*Ly*Lz)/(Nx*Ny*Nz/2L)
          u[pid]    = P1/rho1/(ga-1.0)
        endif 
      
        ;right
        if (pos[0,pid] ge Lx/2.0) then begin
          mass[pid] = rho2*(Lx/2.0*Ly*Lz)/(Nx*Ny*Nz/2L) 
          u[pid]    = P2/rho2/(ga-1.0) 
        endif 
      endfor
    endfor
  endfor
  
  ; create tracer arrays
  nTrLeft  = tfac*Nx*Ny*Nz
  nTrRight = tfac*Nx*Ny*Nz/4
  
  nTotTracers = nTrLeft + nTrRight
  
  id2   = lindgen(nTotTracers) + nTotGas + 1L
  pos2  = fltarr(3,nTotTracers)
  vel2  = fltarr(3,nTotTracers)
  mass2 = fltarr(nTotTracers)
  ;u2    = fltarr(nTotTracers) ; u not expected in IC for tracer
  
  ; set tracer properties (right)
  tfac = fix(tfac^(1.0/3.0))
  
  NxTr = tfac*Nx/4
  NyTr = tfac*Ny
  NzTr = tfac*Nz
  
  deltax = Lx/NxTr/2.0
  deltay = Ly/NyTr
  deltaz = Lz/NzTr
  
  prevCount = 0L
  
  for i=0L,NzTr-1 do begin
    for j=0L,NyTr-1 do begin
      for k=0L,NxTr-1 do begin
        pid = i*(NxTr*NyTr) + j*NxTr + k

        pos2[0,pid] = Lx/2.0 + k*deltax+deltax/2.0
        pos2[1,pid] = j*deltay+deltay/2.0
        pos2[2,pid] = i*deltaz+deltaz/2.0
        
        prevCount += 1
      endfor
    endfor
  endfor
  
  ; set tracer properties (left)
  NxTr = tfac*Nx
  NyTr = tfac*Ny
  NzTr = tfac*Nz
  
  deltax = Lx/NxTr/2.0
  
  for i=0L,NzTr-1 do begin
    for j=0L,NyTr-1 do begin
      for k=0L,NxTr-1 do begin
        pid = i*(NxTr*NyTr) + j*NxTr + k + prevCount
        
        pos2[0,pid] = k*deltax+deltax/2.0  
        pos2[1,pid] = j*deltay+deltay/2.0  
        pos2[2,pid] = i*deltaz+deltaz/2.0
      endfor
    endfor
  endfor
  
  ; offset tracers from gas by a small amount
  seed = 45L
  offset_x = ( deltax/4.0 * randomu(seed,1) )[0]
  offset_y = ( deltay/4.0 * randomu(seed,1) )[0]
  offset_z = ( deltaz/4.0 * randomu(seed,1) )[0]

  pos2[0,*] = pos2[0,*] + offset_x
  pos2[1,*] = pos2[1,*] + offset_y
  pos2[2,*] = pos2[2,*] + offset_z
  
  ; quick plot
  start_PS,'ics_xy.eps'
    fsc_plot,pos[0,*],pos[1,*],psym=4,symsize=0.2,xrange=[0.0,Lx],yrange=[0.0,Ly],/xs,/ys,color=fsc_color('forest green')
    fsc_plot,pos2[0,*],pos2[1,*],psym=4,symsize=0.4,color=fsc_color('crimson'),/overplot
  end_PS, pngResize=50, /deletePS
  start_PS,'ics_xz.eps'
    fsc_plot,pos[0,*],pos[2,*],psym=4,symsize=0.2,xrange=[0.0,Lx],yrange=[0.0,Lz],/xs,/ys,color=fsc_color('forest green')
    fsc_plot,pos2[0,*],pos2[2,*],psym=4,symsize=0.4,color=fsc_color('crimson'),/overplot
  end_PS, pngResize=50, /deletePS
  start_PS,'ics_yz.eps'
    fsc_plot,pos[1,*],pos[2,*],psym=4,symsize=0.2,xrange=[0.0,Ly],yrange=[0.0,Lz],/xs,/ys,color=fsc_color('forest green')
    fsc_plot,pos2[1,*],pos2[2,*],psym=4,symsize=0.4,color=fsc_color('crimson'),/overplot
  end_PS, pngResize=50, /deletePS
  
  ; debug
  print, nTotGas+nTotTracers
  print, n_elements(pos(0,*)) + n_elements(pos2(1,*))
  print, n_elements(vel(0,*)) + n_elements(vel2(1,*))
  print, n_elements(id) + n_elements(id2)
  
  ; write
  gas    = {pos:pos,vel:vel,id:id,mass:mass,u:u}
  tracer = {pos:pos2,vel:vel2,id:id2,mass:mass2}
  
  writeICFile,fOut,part0=gas,part2=tracer
  
end

; 2D shocktube with tracer particles
; BoxSize=2.0 PERIODIC TWODIMS LONG_X=10.0 GAMMA=1.4

pro gen_shocktube_tracer_2D_ICs

  ; config
  fOut = "ics.dat"
  
  ga   = 1.4
  
  P1   = 1.0
  rho1 = 1.0
  P2   = 0.1795
  rho2 = 0.25
  
  Lx = 20.0
  Ly = 2.0
  
  Nx   = 100L
  Ny   = 10L
  
  tfac = 1L
  
  ; deriv
  if (( round(Nx*Ny/4) ne (Nx*Ny/4) ) or ( round(sqrt(tfac)) ne (sqrt(tfac)) ) ) then begin
    print,'Error: Choose nicer numbers.'
    return
  endif
  
  nTotGas = Nx*Ny
  
  deltax = Lx/Nx
  deltay = Ly/Ny
  
  ; create gas arrays
  id   = lindgen(nTotGas)+1L
  pos  = fltarr(3,nTotGas)
  vel  = fltarr(3,nTotGas)
  mass = fltarr(nTotGas)
  u    = fltarr(nTotGas)
  
  ; set gas properties
  for i=0L,Nx-1 do begin
    for j=0L,Ny-1 do begin
      pid = i+j*Nx
      
      pos[0,pid] = i*deltax+deltax/2.0  
      pos[1,pid] = j*deltay+deltay/2.0  
      pos[2,pid] = 0.0
  
      ;left  
      if (pos[0,pid] lt Lx/2.0) then begin
        mass[pid] = rho1*(Lx/2.0*Ly)/(Nx*Ny/2L)
        u[pid]    = P1/rho1/(ga-1.0)
      endif 
    
      ;right
      if (pos[0,pid] ge Lx/2.0) then begin
        mass[pid] = rho2*(Lx/2.0*Ly)/(Nx*Ny/2L) 
        u[pid]    = P2/rho2/(ga-1.0) 
      endif 
      
    endfor
  endfor
  
  ; create tracer arrays
  nTrLeft  = tfac*Nx*Ny
  nTrRight = tfac*Nx*Ny/4
  
  nTotTracers = nTrLeft + nTrRight
  
  id2   = lindgen(nTotTracers) + nTotGas + 1L
  pos2  = fltarr(3,nTotTracers)
  vel2  = fltarr(3,nTotTracers)
  mass2 = fltarr(nTotTracers)
  ;u2    = fltarr(nTotTracers) ; u not expected in IC for tracer
  
  ; set tracer properties (right)
  NxTr = sqrt(tfac)*Nx/4
  NyTr = sqrt(tfac)*Ny
  
  deltax = Lx/NxTr/2.0
  deltay = Ly/NyTr
  
  prevCount = 0
  
  for i=0L,NxTr-1 do begin
    for j=0L,NyTr-1 do begin
      pid = i+j*NxTr
      
      pos2[0,pid] = Lx/2.0 + i*deltax+deltax/2.0  
      pos2[1,pid] = j*deltay+deltay/2.0  
      pos2[2,pid] = 0.0
      
      prevCount += 1
      
    endfor
  endfor
  
  ; set tracer properties (left)
  NxTr = sqrt(tfac)*Nx
  NyTr = sqrt(tfac)*Ny
  
  deltax = Lx/NxTr/2.0
  
  for i=0L,NxTr-1 do begin
    for j=0L,NyTr-1 do begin
      pid = prevCount + i+j*NxTr
      
      pos2[0,pid] = i*deltax+deltax/2.0  
      pos2[1,pid] = j*deltay+deltay/2.0  
      pos2[2,pid] = 0.0
      
    endfor
  endfor
  
  ; offset tracers from gas by a small amount
  seed = 45L
  offset_x = ( deltax/4.0 * randomu(seed,1) )[0]
  offset_y = ( deltay/4.0 * randomu(seed,1) )[0]

  pos2[0,*] = pos2[0,*] + offset_x
  pos2[1,*] = pos2[1,*] + offset_y
  
  ; set tracer properties OLD (left)
  ;pos2[0,0:nTrLeft-1] = findgen(nTrLeft)/(nTrLeft+1) * Lx/2.0
  ;pos2[1,0:nTrLeft-1] = findgen(nTrLeft)/(nTrLeft+1) * Ly
  
  ; set tracer properties OLD (right)
  ;pos2[0,nTrLeft:nTotTracers-1] = findgen(nTrRight)/(nTrRight+1) * Lx/2.0 + Lx/2.0
  ;pos2[1,nTrLeft:nTotTracers-1] = findgen(nTrRight)/(nTrRight+1) * Ly
  
  ; quick plot
  start_PS,'ics.eps'
    fsc_plot,pos[0,*],pos[1,*],psym=4,symsize=0.2,xrange=[0.0,Lx],yrange=[0.0,Ly],/xs,/ys,color=fsc_color('forest green')
    fsc_plot,pos2[0,*],pos2[1,*],psym=4,symsize=0.4,color=fsc_color('crimson'),/overplot
  end_PS
  
  ; debug
  print, nTotGas+nTotTracers
  print, n_elements(pos(0,*)) + n_elements(pos2(1,*))
  print, n_elements(vel(0,*)) + n_elements(vel2(1,*))
  print, n_elements(id) + n_elements(id2)
  
  ; write
  gas    = {pos:pos,vel:vel,id:id,mass:mass,u:u}
  tracer = {pos:pos2,vel:vel2,id:id2,mass:mass2}
  
  writeICFile,fOut,part0=gas,part2=tracer
  
end

; 1D shocktube with tracer particles
; BoxSize=2.0
; PERIODIC ONEDIMS LONG_X=10.0 GAMMA=1.4

pro gen_shocktube_tracer_1D_ICs

  ; config
  fOut = "ics.dat"
  
  ga   = 1.4
  
  P1   = 1.0
  rho1 = 1.0
  P2   = 0.1795
  rho2 = 0.25
  
  Lx = 20.0
  
  Nx   = 100L
  
  tfac = 1L
  
  ; deriv
  if ( round(Nx/4) ne (Nx/4) ) then begin
    print,'Error: Choose nicer numbers.'
    return
  endif
  
  nTotGas = Nx
  deltax  = Lx/Nx
  
  ; create gas arrays
  id   = lindgen(nTotGas)+1L
  pos  = fltarr(3,nTotGas)
  vel  = fltarr(3,nTotGas)
  mass = fltarr(nTotGas)
  u    = fltarr(nTotGas)
  
  ; set gas properties
  for i=0L,Nx-1 do begin
      
      pos[0,i] = i*deltax+deltax/2.0  
      pos[1,i] = 0.0
      pos[2,i] = 0.0
  
      ;left  
      if (pos[0,i] lt Lx/2.0) then begin
        mass[i] = rho1*(Lx/2.0)/(Nx/2.0)
        u[i]    = P1/rho1/(ga-1.0)
      endif 
    
      ;right
      if (pos[0,i] ge Lx/2.0) then begin
        mass[i] = rho2*(Lx/2.0)/(Nx/2.0) 
        u[i]    = P2/rho2/(ga-1.0) 
      endif 
      
  endfor
  
  ; create tracer arrays
  nTrLeft  = tfac*Nx
  nTrRight = tfac*Nx/4
  
  nTotTracers = nTrLeft + nTrRight
  
  id2   = lindgen(nTotTracers) + nTotGas + 1L
  pos2  = fltarr(3,nTotTracers)
  vel2  = fltarr(3,nTotTracers)
  mass2 = fltarr(nTotTracers)
  ;u2    = fltarr(nTotTracers) ; u not expected in IC for tracer
  
  ; set tracer properties (right)
  NxTr   = tfac*Nx/4
  deltax = Lx/NxTr/2.0
  
  prevCount = 0
  
  for i=0L,NxTr-1 do begin
      pos2[0,i] = Lx/2.0 + i*deltax+deltax/2.0  
      pos2[1,i] = 0.0
      pos2[2,i] = 0.0
      
      prevCount += 1
  endfor
  
  ; set tracer properties (left)
  NxTr   = tfac*Nx
  deltax = Lx/NxTr/2.0
  
  for i=0L,NxTr-1 do begin
      pos2[0,prevCount + i] = i*deltax+deltax/2.0  
      pos2[1,prevCount + i] = 0.0  
      pos2[2,prevCount + i] = 0.0
  endfor
  
  ; offset tracers from gas by a small amount
  seed = 45L
  offset_x = ( deltax/4.0 * randomu(seed,1) )[0]

  pos2[0,*] = pos2[0,*] + offset_x
  
  ; quick plot
  start_PS,'ics.eps'
    fsc_plot,pos[0,*],pos[1,*],psym=4,symsize=0.2,xrange=[0.0,Lx],yrange=[-1.0,1.0],$
             /xs,/ys,color=fsc_color('forest green')
    fsc_plot,pos2[0,*],pos2[1,*],psym=4,symsize=0.4,color=fsc_color('crimson'),/overplot
  end_PS  
  
  ; write
  gas    = {pos:pos,vel:vel,id:id,mass:mass,u:u}
  tracer = {pos:pos2,vel:vel2,id:id2,mass:mass2}
  stop
  writeICFile,fOut,part0=gas,part2=tracer
  
end

; ---------------------------------------------------------------------
; alternative 1D blastwave (Toro #3):
; ---------------------------------------------------------------------
;  ga=1.4
;  P1=1000.0
;  P2=0.01
;  P3=100.0
;  rho=1.0
;  
;  Lx=1.0
;  
;  Nx=400L
;  
;  pos=fltarr(3,Nx)
;  vel=fltarr(3,Nx)
;  mass=fltarr(Nx)
;  u=fltarr(Nx)
;  id=lindgen(Nx)+1L
;  
;  deltax=Lx/Nx
;  
;  for i=0L,Nx-1 do begin
;    pos(0,i)=i*deltax+deltax/2.0
;    pos(1,i)=0.0
;    pos(2,i)=0.0
;  
;    ;middle
;    mass(i)=rho*Lx/Nx
;    u(i)=P2/rho/(ga-1.0)
;  
;    ;left
;    if (pos(0,i) lt 0.1) then begin
;     u(i)=P1/rho/(ga-1.0)
;    endif
;  
;    ;right
;    if (pos(0,i) gt 0.9) then begin
;     u(i)=P3/rho/(ga-1.0)
;    endif
;  endfor

; ---------------------------------------------------------------------
; original 2D KH:
; ---------------------------------------------------------------------
;  ga=5./3.
;  P0=2.5
;  rho1=1.
;  rho2=2.
;  
;  v=0.5
;  sigma=0.05/sqrt(2.)
;  w0=0.1
;  
;  Lx=1.
;  Ly=1.
;  Nx=512L
;  Ny=512L
;  
;  pos=fltarr(3,Nx*Ny)
;  vel=fltarr(3,Nx*Ny)
;  mass=fltarr(Nx*Ny)
;  u=fltarr(Nx*Ny)
;  id=lindgen(Nx*Ny)+1L
;  
;  deltax=Lx/Nx
;  deltay=Ly/Ny
;  
;  for i=0L,Nx-1 do begin
;   for j=0L,Ny-1 do begin
;    pos(0,i+j*Nx)=i*deltax+deltax/2.0
;    pos(1,i+j*Nx)=j*deltay+deltay/2.0
;    pos(2,i+j*Nx)=0.0
;  
;    x=pos(0,i+j*Nx)
;    y=pos(1,i+j*Nx)
;  
;    ;top/bottom
;    if (abs(y-0.5) ge 0.25) then begin
;     mass(i+j*Nx)=rho1*(Lx*Ly/2.0)/(Nx*Ny/2L)
;     P=P0
;     u(i+j*Nx)=P/rho1/(ga-1.0)
;     vel(0,i+j*Nx)=-v
;    endif
;  
;    ;stripe
;    if (abs(y-0.5) lt 0.25) then begin
;     mass(i+j*Nx)=rho2*(Lx*Ly/2.0)/(Nx*Ny/2L)
;     P=P0
;     u(i+j*Nx)=P/rho2/(ga-1.0)
;     vel(0,i+j*Nx)=+v
;    endif
;  
;    vel(1,i+j*Nx)=w0*sin(4.*!PI*x)*(exp(-(y-0.25)^2./(2.*sigma^2.)) + exp(-(y-0.75)^2./(2.*sigma^2.)))
;  
;   endfor
;  endfor

; ---------------------------------------------------------------------
; original 2D RT:
; BoxSize=0.5 PERIODIC TWODIMS REFLECTIVE_Y GAMMA=1.4 EXTERNALGRAVITY EXTERNALGY=-0.1 LONG_Y=3.0
; ---------------------------------------------------------------------
; 
;  ga=1.4
;  P0=2.5
;  rho1=2.
;  rho2=1.
;  
;  g=-0.1
;  w0=0.0025
;  
;  Lx=0.5
;  Ly=1.5
;  Nx=50L
;  Ny=150L
;  
;  pos=fltarr(3,Nx*Ny)
;  vel=fltarr(3,Nx*Ny)
;  mass=fltarr(Nx*Ny)
;  u=fltarr(Nx*Ny)
;  id=lindgen(Nx*Ny)+1L
;  
;  deltax=Lx/Nx
;  deltay=Ly/Ny
;  
;  for i=0L,Nx-1 do begin
;   for j=0L,Ny-1 do begin
;    pos(0,i+j*Nx)=i*deltax+deltax/2.0
;    pos(1,i+j*Nx)=j*deltay+deltay/2.0
;    pos(2,i+j*Nx)=0.0
;  
;    x=pos(0,i+j*Nx)
;    y=pos(1,i+j*Nx)
;  
;    ;top
;    if (y gt Ly/2.0) then begin
;     mass(i+j*Nx)=rho1*(Lx*Ly/2.0)/(Nx*Ny/2L)
;     P=P0 + g*(y-0.75)*rho1
;     u(i+j*Nx)=P/rho1/(ga-1.0)
;    endif
;  
;    ;bottom
;    if (y le Ly/2.0) then begin
;     mass(i+j*Nx)=rho2*(Lx*Ly/2.0)/(Nx*Ny/2L)
;     P=P0 + g*(y-0.75)*rho2
;     u(i+j*Nx)=P/rho2/(ga-1.0)
;    endif
;  
;    ;velocity perturbation to excite
;    vel(1,i+j*Nx)=w0*(1.-cos(4.*!PI*x))*(1.-cos(4.*!PI*y/3.))
;   endfor
;  endfor

; ---------------------------------------------------------------------
; original 3D Evrard collapse:
; BoxSize=10.0 PERIODIC GAMMA=1.4
; ---------------------------------------------------------------------
; 
;  Lx=10.
;  Ly=10.
;  Lz=10.
;  xc=5.
;  yc=5.
;  zc=5.
;  ga=1.4
;  rho_hot=1.
;  rho_cold=10.
;  P=1.
;  
;  Machnumber=2.7
;  cs=sqrt(ga*P/rho_hot)
;  vext=Machnumber*cs
;  
;  R_cold=1.
;  
;  vol_hot=Lx*Ly*Lz - 4.*!PI/3.*R_cold^3.
;  vol_cold=4.*!PI/3.*R_cold^3.
;  
;  N_cold=2000L
;  mpart=rho_cold*vol_cold/N_cold
;  N_hot=round(rho_hot*vol_hot/mpart)
;  
;  N_total=N_cold + N_hot
;  
;  pos=fltarr(3,N_total)
;  vel=fltarr(3,N_total)
;  mass=fltarr(N_total)
;  u=fltarr(N_total)
;  id=lindgen(N_total)+1L
;  
;  mass+=mpart
;  
;  seed=42L
;  for n=0L, N_hot-1 do begin
;    repeat begin
;     x=Lx*randomu(seed,1)
;     y=Ly*randomu(seed,1)
;     z=Lz*randomu(seed,1)
;     r=sqrt((x-xc)^2. + (y-yc)^2. + (z-zc)^2.)
;    endrep until r gt R_cold
;  
;    pos(0,n)=x
;    pos(1,n)=y
;    pos(2,n)=z
;    u(n)=P/rho_hot/(ga-1.0)
;  
;    vel(0,n)=vext
;  
;  endfor
;  
;  for n=0L, N_cold-1 do begin
;    repeat begin
;     x=xc+R_cold*(randomu(seed,1)-0.5)*2.
;     y=yc+R_cold*(randomu(seed,1)-0.5)*2.
;     z=zc+R_cold*(randomu(seed,1)-0.5)*2.
;     r=sqrt((x-xc)^2. + (y-yc)^2. + (z-zc)^2.)
;    endrep until r lt R_cold
;  
;    pos(0,n+N_hot)=x
;    pos(1,n+N_hot)=y
;    pos(2,n+N_hot)=z
;    u(n+N_hot)=P/rho_cold/(ga-1.0)
;  endfor
