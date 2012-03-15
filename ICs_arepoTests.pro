; ICs_arepoTests.pro
; initial condition generation - misc idealized tests
; dnelson jan.2012

@helper

; addSubboxBackGrid(): add a single nested background grid to a snapshot

pro addSubboxBackGrid

  input_fname  = "snap_subbox_690.hdf5"
  output_fname = "snap_subbox_690"

  subboxSize = 4000.0 ;ckpc
  subboxCen  = [5500,7000,7500] ;ckpc
  nBackGrid  = 16 ;^3, nested twice
  
  ; load input sphere
  h = loadSnapshotHeader(fileName=input_fname,/verbose)

  gas_pos  = loadSnapshotSubset(fileName=input_fname,partType='gas',field='pos')
  gas_mass = loadSnapshotSubset(fileName=input_fname,partType='gas',field='mass')
  gas_vel  = loadSnapshotSubset(fileName=input_fname,partType='gas',field='vel')
  gas_id   = loadSnapshotSubset(fileName=input_fname,partType='gas',field='ids')
  gas_u    = loadSnapshotSubset(fileName=input_fname,partType='gas',field='u')
  
  ; offset template position to box center
  ;gas_pos += boxSize/2.0
  
  ; form gas
  gas = {pos:gas_pos,vel:gas_vel,id:gas_id,u:gas_u,mass:gas_mass}

  ; add nested background grid
  gas = addICBackgroundGrid(gas, boxSize=subboxSize*1.5, boxCen=subboxCen, nBackGrid=nBackGrid)
  gas = addICBackgroundGrid(gas, boxSize=h.boxSize, nBackGrid=nBackGrid)

  ; save
  writeICFile,output_fname,part0=gas,/longIDs
  print,'Suggested target gas mass: '+string(h.masstable[0])

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


; 3D Evrard collapse: BoxSize=10.0 PERIODIC
; 
; this makes a constant density cold sphere surrounded by a constant density hot medium
; which fills the entire box, NOT the density profile rho(r) = M_cold / (2pi R_cold^2 r) for r<R_cold

pro gen_evrard_collapse_3D_ICs

  ; config
  fOut = "evrard_3D_10k_box=50_rhohot=0.001.dat"
  
  Lx = 50.0
  Ly = 50.0
  Lz = 50.0
  
  ga       = 5.0/3.0
  P        = 1.0
  
  rho_hot  = 0.001
  rho_cold = 10.0 ;30.0 -> utherm_cold=0.05
                  ;10.0 -> utherm_cold=0.15
  
  R_cold = 1.0  
  N_cold = 10000L
  
  ; derived
  xc = Lx/2
  yc = Ly/2
  zc = Lz/2
  
  ; don't know what the purpose of this hot medium (x-direction only) velocity was?
  ;Machnumber = 2.7
  ;cs   = sqrt(ga*P/rho_hot)
  ;vext = Machnumber*cs
  ;print,vext
  
  vol_hot  = Lx*Ly*Lz - 4.0*!pi/3.0*R_cold^3.0
  vol_cold = 4.0*!pi/3.0*R_cold^3.0
  
  mPart = rho_cold*vol_cold/N_cold
  
  N_hot   = round(rho_hot*vol_hot/mPart)
  N_total = N_cold + N_hot
  
  print,'num cold hot total',N_cold,N_hot,N_total
  print,'Suggested TargetGasMass: ',mPart
  print,'Suggested MeanVolume: ',vol_hot/N_hot

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
    ;vel[0,n] = vext
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

; 2D uniform medium (for testing domain decomposition, etc)
; BoxSize=2.0 PERIODIC TWODIMS LONG_X=10.0

pro gen_uniform_medium_2D_ICs

  ; config
  fOut = "uniform_2d_20x2.dat"
  
  ga = 2.0
  
  P   = 1.0
  rho = 1.0
 
  Lx  = 20.0
  Ly  = 2.0
  Nx  = long(20)
  Ny  = long(2)

  ; deriv
  uPart    = P/rho/(ga-1.0)
  massPart = rho*(Lx*Ly)/(Nx*Ny)

  nTotGas = Nx*Ny
  
  deltax = Lx/Nx
  deltay = Ly/Ny
  
  ; create gas arrays
  id   = lindgen(nTotGas)+1L
  pos  = fltarr(3,nTotGas)
  vel  = fltarr(3,nTotGas)
  mass = fltarr(nTotGas) + massPart
  u    = fltarr(nTotGas) + uPart
  
  ; set gas properties
  for i=0L,Nx-1 do begin
    for j=0L,Ny-1 do begin
      pid = i+j*Nx
      
      pos[0,pid] = i*deltax+deltax/2.0  
      pos[1,pid] = j*deltay+deltay/2.0  
      pos[2,pid] = 0.0
    endfor
  endfor
  
  ; write gas (partType=0)
  gas    = {pos:pos,vel:vel,id:id,mass:mass,u:u}
  
  writeICFile,fOut,part0=gas

end

; 2D converging flow test
; periodic sinusoidal velocity perturbation, otherwise uniform
; BoxSize=2.0 PERIODIC TWODIMS LONG_X=10.0

pro gen_converging_flow_2D_ICs, gasOnly=gasOnly

  ; config
  fOut = "convFlow_2d_cs5.0_L5_20_ga2.dat"
  
  ga  = 2.0 ;ga=1+2/dof (for monatomic, 2 in 2D, 3 in 3D) (for diatomic, 4 in 2D, 5 in 3D)
  
  P   = 1.0
  rho = 1.0
  
  Lx  = 20.0
  Ly  = 2.0
  Nx  = long(20) ; odd = places one cell exactly centered at x=Lx/2
  Ny  = long(2)
  
  lambda = 5.0 ; wavelength (needs to divide Lx evenly)
  csFac  = 5.0 ; maximum velocity = half the sound speed
  tfac   = 1L  ; 1-to-1 tracer-gas ratio
  
  ; deriv
  uPart    = P/rho/(ga-1.0)
  massPart = rho*(Lx*Ly)/(Nx*Ny)

  print,'Suggested TargetGasMass: ',massPart

  cs = sqrt(ga*P/rho)
  
  nTotGas = Nx*Ny
  
  deltax = Lx/Nx
  deltay = Ly/Ny
  
  ; create gas arrays
  id   = lindgen(nTotGas)+1L
  pos  = fltarr(3,nTotGas)
  vel  = fltarr(3,nTotGas)
  mass = fltarr(nTotGas) + massPart
  u    = fltarr(nTotGas) + uPart
  
  ; set gas properties
  for i=0L,Nx-1 do begin
    for j=0L,Ny-1 do begin
      pid = i+j*Nx
      
      pos[0,pid] = i*deltax+deltax/2.0  
      pos[1,pid] = j*deltay+deltay/2.0  
      pos[2,pid] = 0.0
    
      ; linear velocity (+csFac*cs at x=0, -csFac*cs at x=Lx, zero at x=Lx/2)
      ; vel[0,pid] = (pos[0,pid]/Lx)*(-2.0)*csFac*cs + cs*csFac
      ; sinusoidal velocity (csFac*cs*cos(2pi x/lambda))
      vel[0,pid] = csFac * cs * cos(2*!pi*pos[0,pid]/lambda)
    endfor
  endfor

  if (not keyword_set(gasOnly)) then begin
  
  ; create tracer arrays
  NxTr = tfac*Nx
  NyTr = Ny
  
  nTr = NxTr*NyTr
  
  id2   = lindgen(nTr) + nTotGas + 1L
  pos2  = fltarr(3,nTr)
  vel2  = fltarr(3,nTr)
  mass2 = fltarr(nTr)
  ;u2    = fltarr(nTr) ; u not expected in IC for tracer
  
  ; set tracer properties
  deltax = Lx/NxTr
  deltay = Ly/NyTr
  
  for i=0L,NxTr-1 do begin
    for j=0L,NyTr-1 do begin
      pid = i+j*NxTr
      
      pos2[0,pid] = i*deltax+deltax/2.0  
      pos2[1,pid] = j*deltay+deltay/2.0  
      pos2[2,pid] = 0.0
      
      ; tracers get no velocity from ICs
    endfor
  endfor
  
  ; offset tracers from gas by a small amount
  ;seed = 45L
  ;offset_x = ( deltax/4.0 * randomu(seed,1) )[0]
  ;offset_y = ( deltay/4.0 * randomu(seed,1) )[0]

  ;pos2[0,*] = pos2[0,*] + offset_x
  ;pos2[1,*] = pos2[1,*] + offset_y
  
  ; quick plot
  ;start_PS,'ics.eps'
  ;  fsc_plot,pos[0,*],pos[1,*],psym=4,symsize=0.2,xrange=[0.0,Lx],yrange=[0.0,Ly],/xs,/ys,color=fsc_color('forest green')
  ;  fsc_plot,pos2[0,*],pos2[1,*],psym=4,symsize=0.4,color=fsc_color('crimson'),/overplot
  ;end_PS
  
  ; debug counts
  print, nTotGas+nTr
  print, n_elements(pos(0,*)) + n_elements(pos2(1,*))
  print, n_elements(vel(0,*)) + n_elements(vel2(1,*))
  print, n_elements(id) + n_elements(id2)
  
  ; make struct
  tracer = {pos:pos2,vel:vel2,id:id2,mass:mass2}
  
  endif ;gasOnly
  
  ; write gas (partType=0) and tracers (partType=3)
  gas    = {pos:pos,vel:vel,id:id,mass:mass,u:u}
  
  if (not keyword_set(gasOnly)) then $
    writeICFile,fOut,part0=gas,part3=tracer
  if (keyword_set(gasOnly)) then $
    writeICFile,fOut,part0=gas
  
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

