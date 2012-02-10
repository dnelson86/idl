; ICs_gasSphere.pro
; initial condition generation - gas cloud in static hernquist potential
; from Mark Vogelsberger's python version
; dnelson jan.2012

@helper

; profile/value functions
function GasRho, r
  COMMON gv
  x0=gas_R0
  x=r/HQ_a
  return, HQ_M/(2*!pi*HQ_a^3.0) * 1/(x+x0) * 1/(x+1)^3.0
end

function HaloRho, r
  COMMON gv
  x=r/HQ_a
  return, HQ_M/(2*!pi*HQ_a^3.0) * 1/x * 1/(x+1)^3.0
end

function Rho, r
  COMMON gv
  return, gas_frac*GasRho(r) + (1.0-gas_frac)*HaloRho(r)
end

function GasMass, r
  COMMON gv
  x0=gas_R0
  x=r/HQ_a
  return, HQ_M * ((1-x0)*x*(x0*(2+3*x)-x)/(1+x)^2.0 + 2*x0^2.0*alog(x0*(1+x)/(x0+x))) / (x0-1)^3.0
end

function HaloMass, r
  COMMON gv
  x=r/HQ_a
  return, HQ_M*x^2.0/(1+x)^2.0
end

function Mass, r
  COMMON gv
  return, gas_frac*GasMass(r) + (1.0-gas_frac)*HaloMass(r)
end

function Sigma_Integrand, r
  COMMON gv
  return, G*Mass(r)*Rho(r)/r^2.0
end

function Sigma, r
  COMMON gv
  ;sqrt(quad(Sigma_Integrand, r, INTERPOL_R_MAX_GAS, epsrel=0.1)[0]/Rho(r))
  val = qpint1d('Sigma_Integrand', r, INTERPOL_R_MAX_GAS, epsrel=0.1)
  test = sqrt(val/Rho(r))
  return,test
end

pro gen_ICs_gasSphere

  COMMON gv,gas_R0,HQ_a,HQ_M,gas_frac,INTERPOL_R_MAX_GAS,G ;globals

  ; INPUT PARAMETERS
  fOut     = 'ics.dat'
  
  N_gas    = 200000
  N_halo   = 200000
  
  gas_frac = 0.1
  halo_fac = 1.0/(1.0-gas_frac)          ; increase halo particles by that factor (for runs with NO_GAS_SELFGRAVITY)

  add_halo = 0                           ; add halo particles?
  add_gas  = 1                           ; add gas particles?

  ; Hernquist profile
  HQ_M = 1e3
  HQ_c = 7.0

  ; bounds
  gas_R0     = 0.01                      ; gas core softening [HQ_a]
  R_min_halo = 1e-5                      ; minimum halo sampling radius [HQ_a]
  R_max_halo = 100.0                     ; maximum halo sampling radius [HQ_a]
  R_min_gas  = 1e-5                      ; minimum gas sampling radius [HQ_a]
  R_max_gas  = 5.0                       ; maximum gas sampling radius [HQ_a]
  R_bins     = 1000                      ; number of interpolation points for function evaluation/inversion

  ; rotation
  Lambda  = 0.0                          ; 0.0 for no rotation
  S_omega = 1.0                          ; 0.0 for rigid body rotation

  ; random seed
  seed = 424242L
  
  units = getUnits()
  G = units.G
  
  ; derived numbers
  gas_mass  = HQ_M * gas_frac / N_gas                     ; per cell
  halo_mass = halo_fac * HQ_M * (1.0-gas_frac) / N_halo   ; per dm particle
  HQ_a      = (G * HQ_M / (100 * 0.1 * 0.1))^(1.0/3.0) / $
              HQ_c * sqrt(2 * (alog(1 + HQ_c) - HQ_c / (1 + HQ_c)))
  
  r_s = HQ_a / sqrt(2*(alog10(1+HQ_c)-HQ_c/(1+HQ_c)))
  r200 = r_s * HQ_c
  
  ; Interpolation parameters
  INTERPOL_BINS       = R_bins
  INTERPOL_R_MIN_HALO = HQ_a*R_min_halo   ; minimum halo sampling radius (halo cut below) 
  INTERPOL_R_MAX_HALO = HQ_a*R_max_halo   ; maximum halo sampling radius (halo cut above) 
  INTERPOL_R_MIN_GAS  = HQ_a*R_min_gas    ; minimum gas sampling radius (gas cut below)     
  INTERPOL_R_MAX_GAS  = HQ_a*R_max_gas    ; maximum gas sampling radius (gas cut above)  

  print,'Interpolation.'
  ; invert function: GasMass^-1 = GasRadius 
  radial_bins = exp(findgen(INTERPOL_BINS)*alog(INTERPOL_R_MAX_GAS/INTERPOL_R_MIN_GAS) / $
                INTERPOL_BINS + alog(INTERPOL_R_MIN_GAS))
  mass_bins_gas = GasMass(radial_bins)
  radius_gas = interpol(radial_bins,mass_bins_gas,randomu(seed,N_gas)*max(mass_bins_gas))

  ; invert function: HaloMass^-1 = HaloRadius
  radial_bins = exp(findgen(INTERPOL_BINS)*alog(INTERPOL_R_MAX_HALO/INTERPOL_R_MIN_HALO) / $
                INTERPOL_BINS + alog(INTERPOL_R_MIN_HALO))
  mass_bins_halo = HaloMass(radial_bins)
  radius_halo = interpol(radial_bins,mass_bins_halo,randomu(seed,N_halo)*max(mass_bins_halo))
  
  ; interpolate sigma and set utherm
  radial_bins = exp(findgen(INTERPOL_BINS)*alog(INTERPOL_R_MAX_GAS/INTERPOL_R_MIN_GAS) / $
                INTERPOL_BINS + alog(INTERPOL_R_MIN_GAS))
  sigma_bins = Sigma(radial_bins)
  sigma_interp = interpol(sigma_bins,radial_bins,radius_gas)
  utherm = 1.5 * sigma_interp^2.0

  print,'Random positions.'
  ; generate random positions (gas)
  phi_gas   = 2.0 * !pi * randomu(seed,N_gas)
  theta_gas = asin(2.0 * randomu(seed,N_gas) - 1.0)
  
  x_gas = radius_gas * cos(theta_gas) * cos(phi_gas)
  y_gas = radius_gas * cos(theta_gas) * sin(phi_gas)
  z_gas = radius_gas * sin(theta_gas)
  
  ; random positions (halo)
  phi_halo   = 2.0 * !pi * randomu(seed,N_halo)
  theta_halo = asin(2.0 * randomu(seed,N_halo) - 1.0)
  
  x_halo = radius_halo * cos(theta_halo) * cos(phi_halo)
  y_halo = radius_halo * cos(theta_halo) * sin(phi_halo)
  z_halo = radius_halo * sin(theta_halo)

  ; rotation
  if (Lambda ne 0.0) then begin
    print,'Velocities.'
    
  endif

  ; prepare for snapshot
  pos_gas  = transpose([[x_gas],[y_gas],[z_gas]])
  vel_gas  = fltarr(3,N_gas)
  id_gas   = lindgen(N_gas)+1
  
  gas = {pos:pos_gas,vel:vel_gas,id:id_gas,u:float(utherm)}
  
  pos_dm  = transpose([[x_halo],[y_halo],[z_halo]])
  vel_dm  = fltarr(3,N_halo)
  id_dm   = lindgen(N_halo)+1
  
  dm = {pos:pos_dm,vel:vel_dm,id:id_dm}
  
  ; write snapshot
  massarr = fltarr(6)
  if (add_gas eq 1)  then massarr[0] = gas_mass
  if (add_halo eq 1) then massarr[1] = halo_mass
  
  if (add_gas eq 1 and add_halo eq 1) then writeICFile,fOut,part0=gas,part1=dm,massarr=massarr
  if (add_gas eq 1 and add_halo eq 0) then writeICFile,fOut,part0=gas,massarr=massarr
  if (add_gas eq 0 and add_halo eq 1) then writeICFile,fOut,part1=dm,massarr=massarr

  print, "HQ_M   = ", HQ_M
  print, "HQ_c   = ", HQ_c
  print, "HQ_a   = ", HQ_a," (r_s = ",r_s," r200 = ",r200,")"
  print, "m_gas  = ", gas_mass
  print, "N_gas  = ", N_gas
  if (add_halo eq 1) then begin
    print, "m_halo = ", halo_mass
    print, "N_halo = ", N_halo
  endif
  stop
end