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
  
  N = 10L
  
  ; setup arrays
  pos = fltarr(3,N)
  
  ; uniformly random positions
  pos[0,*] = randomu(seed,N)
  pos[1,*] = randomu(seed,N)

  ; write
  openw,lun,fOut,/GET_LUN
  
  printf,lun,"#"
  printf,lun,"# test input file with "+str(N)+" points in [0,1]x[0,1]"
  printf,lun,"# gen_arepo_cuda_2D_input"
  printf,lun,"#"
  printf,lun,str(N)
  
  for i=0,N-1 do begin
    printf,lun,string(pos[0,i],format='(f6.4)')+" "+string(pos[1,i],format='(f6.4)')
  endfor

  close,lun
  free_lun,lun
  
  print,'wrote ',fOut
  
end

; 2D Kelvin-Helmholtz Instability ICs
pro gen_KH_2D_ICs

  ; config
  angle = 0.0 ;CCW from 0 wrt positive x-axis
  fOut  = 'KH_Instability_0.dat'
  
  stripeWidth = 0.5
  
  ga   = 5.0/3.0
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
  writeICFile,fOut,pos,vel,id,dens,u

end

; 3D Taylor-Sedov (point explosion) blastwave IC
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
  writeICFile,fOut,pos,vel,id,mass,u
  
end
