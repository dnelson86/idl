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

; 2D shocktube with tracer particles
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
  
  tfac = 16L
  
  ; deriv
  if (( round(Nx*Ny/4) ne (Nx*Ny/4) ) or ( round(sqrt(tfac)) ne (sqrt(tfac)) ) ) then begin
    print,'Error: Choose nicer numbers.'
    return
  endif
  
  nTotGas = Nx*Ny
  
  deltax = Lx/Nx
  deltay = Ly/Ny
  
  ; create gas arrays
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
  
  pos2  = fltarr(3,nTotTracers)
  vel2  = fltarr(3,nTotTracers)
  mass2 = fltarr(nTotTracers)
  u2    = fltarr(nTotTracers)
  
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
  
  start_PS,'ics.eps'
    fsc_plot,pos[0,*],pos[1,*],psym=4,symsize=0.2,xrange=[0.0,Lx],yrange=[0.0,Ly],/xs,/ys,color=fsc_color('forest green')
    fsc_plot,pos2[0,*],pos2[1,*],psym=4,symsize=0.4,color=fsc_color('crimson'),/overplot
  end_PS
  
  ; ids
  id   = lindgen(nTotGas + nTotTracers)+1L
  
  ; header
  npart    = lonarr(6) 
  massarr  = dblarr(6)
  npartall = lonarr(6)
  
  npart[0] = nTotGas
  npart[2] = nTotTracers
  npartall[0] = nTotGas
  npartall[2] = nTotTracers
  
  time          = 0.0D
  redshift      = 0.0D
  flag_sfr      = 0L
  flag_feedback = 0L
  bytesleft     = 136
  la            = intarr(bytesleft/2)
  
  ; debug
  print, total(npartall,/int)
  print, n_elements(pos(0,*)) + n_elements(pos2(1,*))
  print, n_elements(vel(0,*)) + n_elements(vel2(1,*))
  print, n_elements(id)
  
  ; write
  openw,1,fOut,/f77_unformatted
  writeu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall,la
  
  writeu,1, pos,pos2
  writeu,1, vel,vel2
  writeu,1, id
  writeu,1, mass, mass2
  writeu,1, u
  
  close,1
  
end