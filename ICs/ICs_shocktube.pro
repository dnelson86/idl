; ICs_shocktube.pro
; initial condition generation for shocktubes/soundwaves
; dnelson jan.2012

@helper

; 3D shocktube with tracer particles
; BoxSize=2.0 PERIODIC REFLECTIVE_X=1 LONG_X=10.0 GAMMA=1.4

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
  
  tfac = 8L
  
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
; BoxSize=2.0 PERIODIC REFLECTIVE_X=1 TWODIMS LONG_X=10.0 GAMMA=1.4

pro gen_shocktube_tracer_2D_ICs, gasOnly=gasOnly

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
  
  if (not keyword_set(gasOnly)) then begin
  
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
  
  tracer = {pos:pos2,vel:vel2,id:id2,mass:mass2}
  
  endif ;gasOnly
  
  ; write
  gas    = {pos:pos,vel:vel,id:id,mass:mass,u:u}
  
  if (not keyword_set(gasOnly)) then $
    writeICFile,fOut,part0=gas,part2=tracer
  if (keyword_set(gasOnly)) then $
    writeICFile,fOut,part0=gas
  
end

; 1D shocktube with tracer particles
; BoxSize=2.0
; PERIODIC REFLECTIVE_X=1 ONEDIMS LONG_X=10.0 GAMMA=1.4

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