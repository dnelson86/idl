; plotTI.pro
; thermal instability IC related
; dnelson 13.march.2011

@helper
@arepoLoad

; TI_1D_IC():
;  - examine Gadget format IC file

pro plotTI_1D_IC, workingPath, fileName

  units = getUnits() 

  ;load
  h = loadSnapshot(fileName,'none',pos,vel,id,mass,u,rho,hsml)

  ; determine xrange
  xStepSize = max(pos[0,*])/n_elements(id)
  xrange = [0.0,max(pos[0,*])+xStepSize] 
  print,xStepSize,xrange
  
  ;plot in physical units?
  vel_cgs  = vel  * units.UnitVelocity_in_cm_per_s        ;cm/s
  u_cgs    = u    * units.UnitEnergy_in_cgs               ;erg
  mass_cgs = mass * units.UnitMass_in_g / units.Msun_in_g ;Msun
  
  PS_Start, FILENAME=workingPath+"plotICs."+fileName+".eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=5, /inches  
  
    !p.multi = [0,1,2]  
    ym = !y.margin
    xm = !x.margin
    !x.margin = [6.0,5.0]
  
    ; plot 1D density/energy profiles
    !y.margin = [0.0,1.0]
    mass *= 10.0e8
    fsc_plot,pos[0,*],mass,    line=0,xrange=xrange,/xs,charsize=1.5, $
         xtitle="",ytitle=textoidl("mass * 10^8 [code]"),xtickname=replicate(' ',10),$
         yrange=[min(mass)*0.9,max(mass)*1.1],/ys
         
    !y.margin = [2.0,0.0]
    fsc_plot,pos[0,*],u,       line=0,xrange=xrange,/xs,charsize=1.5, $
         xtitle="x",ytitle=textoidl("u [code]"),yrange=[min(u)*0.9,max(u)*1.1],/ys
   
    fsc_text,0.75,0.9,fileName,charsize=1.2,alignment=0.5,/normal
   
    !x.margin = xm
    !y.margin = ym      
    !p.multi = 0
  
  PS_End

end

; thermal instability in 1D test ICs
; "initialized with eigenmodes of the instability" - adapted from Piontek & Ostriker (2004)
;
; grid = 1024/4k/8k zones
; box  = 20/100 pc
; T    = 2400 K     initial temperature
; n    = 1 cm^(-3)  initial density
; k    = 0.5/1/2/4/6/8 * k_rho
; w0   = 0.1/0.01

pro gen_TI_1D_ICs_eigen, fOut, Lx, k_fac, w0

  units = getUnits()

  ; config
  ga = 1.667
  mu = 1.22
  
  Nx = 1024L  
  ; w0 = 0.1 ;amplitude of eigenmode perturbation (in percent, 0.1, 0.01, ...)  
  
  ; Lx    = boxsize (0.1=100pc)
  ; k_fac = set eigenmode wavenumber (multiplier of the eigenvalue normalization k_rho)
  
  ; setup ICs
  ;cgs cut from PlotGrowthRateSS02 for: P_k0 = 2400.0,  n_0 = 1.0, gamma = 5.0/3.0, mu = 1.22
  k_rho = 1.38e-19
  k_rho *= units.UnitLength_in_cm ;code
  kvector = k_fac * k_rho

  ; make sure periodic boundary is the same [0,L]
  lambda_rho = 2.0 * !pi / (k_rho*k_fac) ;code
  cycles     = floor(Lx / lambda_rho)
  offset     = (Lx-cycles*lambda_rho) / 2.0 ;center in box
  
  ; debug:
    print,'kvector (code) ',kvector
    print,'lambda_rho ',lambda_rho
    print,'cycles ',cycles
    print,'offset ',offset

  ; uniform initial conditions
  T1   = 2400.0   ;K
  n1   = 1.0      ;cm^(-3)

  ; derived ICs (cgs)
  rho1 = n1 * (mu * units.mass_proton)
  ;temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;
  u1 = T1 * units.boltzmann / (ga-1.0) / (units.mass_proton * mu)
  P1 = u1 * rho1 * (ga-1.0)
  ;P1 = rho1 / (mass_proton*mu) * boltzmann * T1 ;equivalent

  ; convert to code units
  rho1 = rho1 / units.UnitDensity_in_cgs
  P1   = P1   / units.UnitPressure_in_cgs
  
  ; setup arrays
  pos  = fltarr(3,Nx)
  vel  = fltarr(3,Nx)
  dens = fltarr(Nx)
  u    = fltarr(Nx)
  id   = lindgen(Nx)+1L
 
  deltax = Lx/Nx

  ; set each particle
  for i=0L,Nx-1 do begin
    pos(0,i) = i*deltax+deltax/2.0  
    pos(1,i) = 0.0
    pos(2,i) = 0.0

    ; eigenmode density fluctuation
    rho_i = rho1
    if (pos[0,i] ge offset and pos[0,i] le (cycles*lambda_rho+offset)) then $
      rho_i = rho1 + (w0*rho1) + (w0 * rho1) * cos( (pos[0,i]-offset) * kvector + !pi)

    ; spectrum of density fluctuations in [k_min,k_max]
    ;
    
    ; Gaussian density profile (SS02)
    ;

    ; random density fluctuation
    ;rho_i = rho1 * (randomu(seed,1) * (2*w0) + (1-w0))

    ;avg properties  
    dens[i] = rho_i
    u[i]    = P1/rho_i/(ga-1.0)

  endfor


  ;debug output
  print,'minmax dens ',minmax(dens)
  print,'minmax u ',minmax(u)
  
  ;write
  writeICFile,fOut,pos,vel,id,dens,u

end

; gen_tcond_test_ICs():
;  - first pass at some VERY simple ICs to test THERMAL_CONDUCTION
;  - uniform 1D with 1 grid size tophat in middle
;  - write density as mass (MESHRELAX_DENSITY_IN_INPUT)
pro gen_tcond_test_ICs

  units = getUnits()
  
  ; config
  fOut = 'IC.tcond.test.100.dat'

  Nx    = 100L
  Lx    = 0.01   ; 10pc
  
  n0 = 1.0 ;cgs
  
  rho_0 = n0 * units.mass_proton ;cgs
  rho_0_code = rho_0 / units.UnitDensity_in_cgs

  ; setup arrays
  pos  = fltarr(3,Nx)
  vel  = fltarr(3,Nx)
  dens = fltarr(Nx)
  u    = fltarr(Nx)
  id   = lindgen(Nx)+1L
  
  deltax = Lx/Nx

  ; set each particle
  for i=0L,Nx-1 do begin
    pos[0,i] = i*deltax+deltax/2.0  
    pos[1,i] = 0.0
    pos[2,i] = 0.0

    dens[i] = rho_0_code
    u[i]    = 10.0
  endfor

  ; set tophat
  w = floor(Nx/2)
  u[w] *= 2

  ;debug output
  print,'dens ',dens
  print,'u ',u
  
  ;write
  writeICFile,fOut,pos,vel,id,dens,u

end
