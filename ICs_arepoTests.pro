; ICs_arepoTests.pro
; initial condition generation - misc idealized tests
; dnelson jan.2012

@helper

; add background grid of specified resolution to gas ICs
function addICBackgroundGrid, gas, boxSize=boxSize, nBackGrid=nBackGrid

  ; if not requested, return un-altered
  if (nBackGrid eq 0) then return, gas
  
  if (n_elements(boxSize) eq 0 or n_elements(gas) eq 0) then stop

  ; config
  massBackGrid   = 1e-20
  uthermBackGrid = 0.0

  backCellSize = boxSize / nBackGrid
  
  xyz_back = findgen(nBackGrid)/nBackGrid * boxSize + backCellSize/2.0
  
  nBackKeep = 0
  pos_back = []
  
  ; find empty background grid cells
  for i=0,nBackGrid-1 do begin
    for j=0,nBackGrid-1 do begin
      for k=0,nBackGrid-1 do begin
        cenBackCell = [xyz_back[i],xyz_back[j],xyz_back[k]]
        min_xyz = cenBackCell - backCellSize/2.0
        max_xyz = cenBackCell + backCellSize/2.0
        
        w = where(gas.pos[0,*] ge min_xyz[0] and gas.pos[0,*] le max_xyz[0] and $
                  gas.pos[1,*] ge min_xyz[1] and gas.pos[1,*] le max_xyz[1] and $
                  gas.pos[2,*] ge min_xyz[2] and gas.pos[2,*] le max_xyz[2], count)
        
        ; keep background point
        if (count eq 0) then begin
          nBackKeep += 1
          pos_back = [[pos_back],[cenBackCell]]
        endif
      endfor
    endfor
  endfor
  
  ; create other arrays
  vel_back  = fltarr(3,nBackKeep)
  mass_back = fltarr(nBackKeep) + massBackGrid
  u_back    = fltarr(nBackKeep) + uthermBackGrid
  id_back   = lindgen(nBackKeep) + max(id) + 1
  
  ; concat background grid and primary cells
  pos  = [[gas.pos],[pos_back]]
  vel  = [[gas.vel],[vel_back]]
  mass = [gas.mass,mass_back]
  u    = [gas.u,u_back]
  id   = [gas.id,id_back]

  r = {pos:pos,vel:vel,mass:mass,u:u,id:id}
  return, r
  
end

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


; original 3D Evrard collapse: BoxSize=10.0 PERIODIC GAMMA=1.4
pro gen_evrard_collapse_3D_ICs

  ; config
  fOut = "evrard_3D_24k.dat"
  
  Lx = 10.0
  Ly = 10.0
  Lz = 10.0
  
  ga       = 5.0/3.0
  rho_hot  = 1.0
  rho_cold = 10.0
  P        = 1.0
  
  Machnumber = 2.7
  
  R_cold = 1.0  
  N_cold = 24000L
  
  ; derived
  xc = Lx/2
  yc = Ly/2
  zc = Lz/2
  
  cs   = sqrt(ga*P/rho_hot)
  vext = Machnumber*cs
  
  vol_hot  = Lx*Ly*Lz - 4.0*!pi/3.0*R_cold^3.0
  vol_cold = 4.0*!pi/3.0*R_cold^3.0
  
  mPart = rho_cold*vol_cold/N_cold
  
  N_hot   = round(rho_hot*vol_hot/mPart)
  N_total = N_cold + N_hot
  
  print,'num cold hot total',N_cold,N_hot,N_total
  
  ; arrays
  pos  = fltarr(3,N_total)
  vel  = fltarr(3,N_total)
  mass = fltarr(N_total)
  u    = fltarr(N_total)
  id   = lindgen(N_total)+1L
  
  mass += mPart
  
  seed = 42L
  
  ; generate hot component
  for n=0L, N_hot-1 do begin
    repeat begin
     x = Lx*randomu(seed,1)
     y = Ly*randomu(seed,1)
     z = Lz*randomu(seed,1)
     r = sqrt((x-xc)^2.0 + (y-yc)^2.0 + (z-zc)^2.0)
    endrep until r gt R_cold
  
    pos[0,n] = x
    pos[1,n] = y
    pos[2,n] = z
    u[n]     = P/rho_hot/(ga-1.0)
    vel[0,n] = vext
  endfor
  
  ; generate cold component
  for n=0L, N_cold-1 do begin
    repeat begin
     x = xc+R_cold*(randomu(seed,1)-0.5)*2.0
     y = yc+R_cold*(randomu(seed,1)-0.5)*2.0
     z = zc+R_cold*(randomu(seed,1)-0.5)*2.0
     r = sqrt((x-xc)^2.0 + (y-yc)^2.0 + (z-zc)^2.0)
    endrep until r lt R_cold
  
    pos[0,n+N_hot] = x
    pos[1,n+N_hot] = y
    pos[2,n+N_hot] = z
    u[n+N_hot]     = P/rho_cold/(ga-1.0)
  endfor
  
  ; add background grid
  gas = {pos:pos,vel:vel,id:id,mass:mass,u:u}
  
  nBackGrid = 0
  
  gasWithBack = addICBackgroundGrid(gas,boxSize=Lx,nBackGrid=nBackGrid)
  
  ; write
  writeICFile,fOut,part0=gasWithBack
  
end

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

