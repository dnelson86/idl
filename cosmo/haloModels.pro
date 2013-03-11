; haloModels.pro
; theoretical models for DM halos (e.g. NFW, SIS) including their gas
; dnelson feb.2013

; interpolate Sutherland and Dopita (1993) primordial/metal line radiative cooling tables

function interpLambdaSD93, Z=Z, logT=logT, tables=tables, norm=norm

  ; if no input tables, then load them now and return
  if ~keyword_set(tables) then begin
    filePath   = '/n/home07/dnelson/idl/data/sd93/'
    fileNames  = ['mzero','m-30','m-25','m-20','m-15','m-10','m-05','m-00','m+05'] + '.cie'
    dataStruct = {logT:0.0,nelec:0.0,nh:0.0,nt:0.0,logLambdaNet:0.0,logLambdaNorm:0.0,$
                  logU:0.0,logTauCool:0.0,P12:0.0,rho24:0.0,Ci:0.0,muBar:0.0}
    nHeaderLines = 4
                  
    table = loadCSV(nHeaderLines,filePath + fileNames[0],dataStruct,header=header)
    return,table ; just Z=0 for now
  endif
  
  if Z ne 0.0 then message,'Error: Primordial only for now.'
  
  ; interpolate in temperature
  if ~keyword_set(norm) then $
    logLambdaNetInterp = interpol(tables.logLambdaNet,tables.logT,logT,/spline) ; lambda_net
  if keyword_set(norm) then $
    logLambdaNetInterp = interpol(tables.logLambdaNorm,tables.logT,logT,/spline) ; lambda_N
  
  if n_elements(logLambdaNetInterp) eq 1 then return, logLambdaNetInterp[0]
  return,logLambdaNetInterp
end

; haloMAH(): fit of Wechsler+ (2002), given halo mass M0 at redshift z0, return mass at earlier redshift z
;            specify only one z, can specify multiple M0,z0 pairs

function haloMAH, z, M0=M0, z0=z0
  if (n_elements(M0) ne n_elements(z0)) or ~n_elements(M0) or n_elements(z) gt 1 then message,'error'
  
  ; convert redshifts to scale factors
  scalefac  = 1.0 / (1+z)
  scalefac0 = 1.0 / (1+z0)
  
  ; interpolate for ac (scalefactor corresponding to the formation time of the M0 halo)
  ac = interpol([0.45,0.53,0.60,0.75],[10.0,100.0,1000.0,10000.0],M0)
  
  return, M0 * exp( -2*ac * (scalefac0/scalefac - 1) )
end

; nfw (Navarro+)
; --------------

; use Prada+ (2011) analytical to calculate halo concentration given mass and redshift
function nfw_c, mass=mass_codeunits, redshift=redshift

  compile_opt idl2, hidden, strictarr, strictarrsubs
  if ~keyword_set(mass_codeunits) or ~n_elements(redshift) then message,'error'
  message,'not quite right'
  
  omega_L0 = 0.73
  omega_m0 = 0.27
  
  ; find x, D(a)
  scalefac = 1 / (1 + redshift)
  
  x_const = (omega_L0/omega_m0)^(1.0/3.0) ; eqn 13
  x_upperbound = x_const * scalefac
  
  D_a = (5.0/2.0) * (omega_m0/omega_L0)^(1.0/3.0) * sqrt(1+x_upperbound^3.0) / x_upperbound^(3.0/2.0)
  D_a *= qpint1d( 'x^(3.0/2.0) / (1+x^3.0)^(3.0/2.0)', 0.0, x_upperbound, /expr )
  
  ; find B0,B1
  c_min_1393 = 3.681 + (5.033-3.681) * ( atan(6.948*(1.393 - 0.424))/!pi + 0.5 ) ; eqn 19,21
  c_min_x    = 3.681 + (5.033-3.681) * ( atan(6.948*(x_upperbound - 0.424))/!pi + 0.5 )
  
  sigma_min_inv_1393 = 1.047 + (1.646-1.047) * ( atan(7.386*(1.393 - 0.526))/!pi + 0.5 ) ; eqn 20,22
  sigma_min_inv_x    = 1.047 + (1.646-1.047) * ( atan(7.386*(x_upperbound - 0.526))/!pi + 0.5 )
  
  B0 = c_min_x / c_min_1393 ; eqn 18
  B1 = sigma_min_inv_x / sigma_min_inv_1393
  
  ; find sigma(M,a), eqn 23
  y = 1.0 / (mass_codeunits / 100.0)
  sigma_M_a = D_a * ( 16.9 * y^0.41 ) / ( 1 + 1.102 * y^0.2 + 6.22 * y^0.333 )
  
  ; find sigma_prime, calC, concentration
  sigma_prime = B1 * sigma_M_a ; eqn 15
  calC = 2.881 * ( (sigma_prime/1.257)^1.022 + 1.0 ) * exp( 0.06/sigma_prime^2.0 ) ; eqn 16,17
  
  nfw_c = B0 * calC ; eqn 14
  
  return,nfw_c

end

; use Munoz-Cuartas+ (2010) fitting formula for halo concentration (z<2)
function nfw_c_fit, mass=mass_codeunits, redshift=redshift
  
  compile_opt idl2, hidden, strictarr, strictarrsubs
  if ~keyword_set(mass_codeunits) or ~n_elements(redshift) then message,'error'
  
  ; a(z), b(z)
  a_z = 0.029 * redshift - 0.097
  b_z = -110.001 / (redshift + 16.885) + 2469.720 / (redshift + 16.885)^2.0
  
  ; calculate concentration
  mass_msun_hinv = mass_codeunits * 1e10
  log_c = a_z * alog10( mass_msun_hinv ) + b_z
  
  return, 10.0^log_c
end

; use Klypin+ (2011) bolshoi fitting formula for halo concentration (z<5)

function nfw_c_fit_klypin, mass=mass_codeunits, redshift=redshift

  compile_opt idl2, hidden, strictarr, strictarrsubs
  if ~keyword_set(mass_codeunits) or ~n_elements(redshift) then message,'error'
  
  ; setup interpolation
  z_interp  = [0.0,0.5,1.0,2.0,3.0,5.0]
  c0_interp = [9.60,7.08,5.45,3.67,2.83,2.34]
  M0_interp = [1e20,1.5e7,2.5e5,6.8e3,6.3e2,6.6e1] ; code units, 1e10 hinv msun (note first -> zero term)
  
  ; interpolate for c0 and M0 parameters
  c0_z = interpol(c0_interp,z_interp,redshift)
  M0_z = interpol(M0_interp,z_interp,redshift)
  
  ; calculate ceontration
  c = c0_z * (mass_codeunits/100)^(-0.075) * ( 1 + (mass_codeunits/M0_z)^0.26 )
  
  return, c
end

; nfw profile given (M,z)
; inputs: rad_frac in fraction of rvir, mass in code units

function nfw_profile, rad_frac, mass=mass_codeunits, redshift=redshift

  compile_opt idl2, hidden, strictarr, strictarrsubs
  if ~keyword_set(mass_codeunits) or ~n_elements(redshift) then message,'error'
  units = getUnits(redshift=redshift)
  
  ; calculate concentration
  c = nfw_c_fit_klypin(mass=mass_codeunits, redshift=redshift)
  
  ; calculate virial radius (r200), Monaco+ (2007) eqn 5
  r200 = (mass_codeunits * units.G / 100 / units.H_z^2.0)^(1.0/3.0)
  
  ; convert input rad to kpc (code units)
  rad = rad_frac * r200
  
  ; calculate r_s and rho_s
  delta_c = (200.0/3.0) * c^3.0 / ( alog(1+c) - c/(1+c) )
  r_s = r200 / c
  
  x200 = r200 / r_s
  rho_s = mass_codeunits / ( 4 * !pi * r_s^3.0 * ( alog(1+x200) - x200/(1+x200) ) )
  ;rho_s = units.rhoCrit * delta_c ; equivalent, rho_crit(z=0)
  
  ; density profile
  x = rad / r_s
  
  rho_DM = rho_s / ( x*(1+x)^2.0 )
  
  ; enclosed mass (analytic)
  Menc_DM   = 4 * !pi * rho_s * r_s^3.0 * ( alog(1+x) - x/(1+x) )
  Menc_rvir = 100.0 * units.H_z^2.0 * r200^3.0 / units.G
  Menc_rs   = 4 * !pi * rho_s * r_s^3.0 * ( alog(2) - 0.5 )
    
  ; circular velocity
  Vcirc_DM   = sqrt(units.G * Menc_DM / rad)
  Vcirc_rvir = sqrt(units.G * Menc_rvir / r200)
  Vcirc_rs   = sqrt(units.G * Menc_rs / r_s)
  
  ; virial temperature
  Tvir_DM   = float( units.mu * units.mass_proton * Vcirc_DM^2.0   * units.UnitVelocity_in_cm_per_s^2.0 / $
                     (3 * units.boltzmann) )
  Tvir_rvir = float(units.mu * units.mass_proton * Vcirc_rvir^2.0 * units.UnitVelocity_in_cm_per_s^2.0 / $
                     (3 * units.boltzmann) )
  Tvir_rs   = float(units.mu * units.mass_proton * Vcirc_rs^2.0   * units.UnitVelocity_in_cm_per_s^2.0 / $
                     (3 * units.boltzmann) )
  
  ; binding energy
  bindEnergy = units.G * mass_codeunits^2.0 * delta_c^2.0 / r200 * $
               ( (alog(1+c))^2.0 - c*c/(1+c) )
  
  return, {rad:rad, rho_DM:rho_DM, $
           Menc_DM:Menc_DM,   Menc_rvir:Menc_rvir,   Menc_rs:Menc_rs,   $
           Vcirc_DM:Vcirc_DM, Vcirc_rvir:Vcirc_rvir, Vcirc_rs:Vcirc_rs, $
           Tvir_DM:Tvir_DM,   Tvir_rvir:Tvir_rvir,   Tvir_rs:Tvir_rs,   $
           c:c, r_s:r_s, rho_s:rho_s, r200:r200, delta_c:delta_c, $
           mass_codeunits:mass_codeunits, redshift:redshift, bindEnergy:bindEnergy }
  
end

; nfw equilibrium hydrostatic profile, following Suto+ (1998) both iso and polytropic solutions

; suto_model(): calculate profiles given (T0,rho0,n) parameters

function suto_model, T_0, rho_0, n, nfw_dm

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits(redshift=nfw_dm.redshift)
  
  B = ( 4 * !pi * units.Gravity * units.mu * units.mass_proton * $
            units.rhoCrit_z*nfw_dm.delta_c * nfw_dm.r_s^2.0 * units.UnitMass_in_g / units.UnitLength_in_cm ) / $
          ( units.boltzmann * T_0 )
  B_p = B/(n+1)
  
  ; construct density/temperature profiles
  x = (nfw_dm.rad/nfw_dm.r_s)  
  
  ; polytropic
  eps_x = 1 - B_p * ( 1 - alog(1 + x) / x )
  
  temp_gas = T_0   * eps_x
  rho_gas  = rho_0 * eps_x^n
  
  ; isothermal
  rho_gas_iso = rho_0 * exp(-B*( 1 - alog(1 + x) / x ))
  
  return, {B:B, B_p:B_p, temp_gas:temp_gas, rho_gas:rho_gas, rho_gas_iso:rho_gas_iso, $
           T_0:T_0, rho_0:rho_0, n:n}
end

; nfw_rcool_func(): helper function for root finding for r_cool

function nfw_rcool_func, X, params=params
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ;nfw_dm_out = mod_struct(params.nfw_dm,"rad",X) ; replace radii with X
  ;suto = suto_model(params.T_0,params.rho_0,params.n,nfw_dm_out)

  ; dynTime
  xx = X / params.nfw_dm.r_s
  
  ; enclosed mass (analytic)
  Menc   = 4 * !pi * params.nfw_dm.rho_s * params.nfw_dm.r_s^3.0 * ( alog(1+xx) - xx/(1+xx) )
  Vcirc  = sqrt(units.G * Menc / X)
  dynTime = X / (Vcirc * units.kmS_in_kpcYr) / 1e9
  
  ; coolTime
  ;alpha = 3 * units.mu / 2
  alpha = 3.56 ; LU
  
  eps_xx = 1 - params.B_p * ( 1 - alog(1 + xx) / xx )
  temp_gas = params.T_0   * eps_xx
  rho_gas  = params.rho_0 * eps_xx^params.n
  
  lambdaNet = 10.0^(interpLambdaSD93(Z=0.0,logT=alog10(temp_gas),tables=params.tables))
  COF_cgs   = float( alpha * units.mass_proton * units.boltzmann * temp_gas / lambdaNet )
  
  coolTime  = float( COF_cgs / (rho_gas*units.UnitDensity_in_cgs) / units.s_in_Gyr )

  ; find minimum of abs(difference)
  y = coolTime - dynTime

  return, y  
end

; vector "root finder"

function rootFindVector, y1, y2, xx=xx

  w = where(abs(y1-y2) eq min(abs(y1-y2)),count)
  if count gt 0 then w = w[0]
  if count eq 0 or w eq 0 then begin
    x = !values.f_nan
  endif else begin
    x = xx[w]
  endelse
  
  if max(y1) lt min(y2) or min(y1) gt max(y2) then x = !values.f_nan
  
  return,x
end

; nfw_timescales_suto(): calculate timescales vs radius and interesting radii for NFW gas model
  
function nfw_timescales_suto, suto=suto, nfw_dm=nfw_dm, tables=tables
  compile_opt idl2, hidden, strictarr, strictarrsubs
  if ~keyword_set(nfw_dm) or ~keyword_set(suto) or ~keyword_set(tables) then message,'error'
  units = getUnits(redshift=nfw_dm.redshift)
  
  ; timescales: poly
  dynTime = nfw_dm.rad / (nfw_dm.Vcirc_DM * units.kmS_in_kpcYr) / 1e9
  
  ;alpha = 3 * units.mu / 2
  alpha = 3.56 ; LU
  
  lambdaNet = 10.0^(interpLambdaSD93(Z=0.0,logT=alog10(suto.temp_gas),tables=tables))
  COF_cgs   = float( alpha * units.mass_proton * units.boltzmann * suto.temp_gas / lambdaNet )
  
  coolTime  = float( COF_cgs / (suto.rho_gas*units.UnitDensity_in_cgs) / units.s_in_Gyr )
  
  ;hubbleTime = 1/units.H_z ; Gyr
  hubbleTime = redshiftToAgeFlat(nfw_dm.redshift)
  
  ; rcool (where coolTime=dynTime) and rcool (where coolTime=hubbleTime)
  r_cool   = rootFindVector(coolTime, dynTime, xx=nfw_dm.rad)
  r_cool_h = rootFindVector(coolTime, hubbleTime, xx=nfw_dm.rad)
  
  ; try actual root find
  params = { T_0:suto.T_0, rho_0:suto.rho_0, n:suto.n, B_p:suto.B_p, nfw_dm:nfw_dm, tables:tables }
  ;r_zero = zbrent(1.0,1000.0,func_name="nfw_rcool_func",_EXTRA={params:params}) ; kpc
  ;r_zero2 = tnmin("nfw_rcool_func",[nfw_dm.r_s],functargs={params:params},/autoderivative,/quiet)
  
  ; timescales: iso
  lambdaNet = 10.0^(interpLambdaSD93(Z=0.0,logT=alog10(suto.T_0),tables=tables))
  COF_cgs   = float(alpha * units.mass_proton * units.boltzmann * suto.T_0 / lambdaNet)
  
  coolTime_iso = COF_cgs / (suto.rho_gas_iso*units.UnitDensity_in_cgs) / units.s_in_Gyr
  
  r_cool_iso = rootFindVector(coolTime_iso, dynTime, xx=nfw_dm.rad)
  
  return, {dynTime:dynTime, coolTime:coolTime, hubbleTime:hubbleTime, r_cool:r_cool, r_cool_h:r_cool_h, $
           coolTime_iso:coolTime_iso, r_cool_iso:r_cool_iso}
           
end

; nfw_gas_suto(): iteratively solve using total mass and binding energy constraints to fix parameters

function nfw_gas_suto, mass_hot=mass_hot_codeunits, nfw_dm=nfw_dm, tables=tables

  compile_opt idl2, hidden, strictarr, strictarrsubs
  if ~keyword_set(nfw_dm) or ~keyword_set(mass_hot_codeunits) then message,'error'
  units = getUnits(redshift=nfw_dm.redshift)

  iterFac = 0.05
  maxIter = 1000
  errTol  = 1e-5 ; relative
  
  ; parameters
  ;f_shock = 1e-7 ; multiplier
  r_cool  = 0.0 ; inner radius for hot halo gas
  n = 20.0 ;3/2,7
  
  ; initial guess
  T_0   = 2.0 * nfw_dm.Tvir_rs
  rho_0 = 100 * nfw_dm.mass_codeunits / (4/3*!pi*nfw_dm.r200^3.0)
  
  ; iterate to find central temp
  for i=0,maxIter do begin
    suto = suto_model(T_0,rho_0,n,nfw_dm)
    
    ; upper limit of B_p<1 assuming gas extends quite far
    if suto.B_p ge 1.0 then print,'Warning: B_p>1'
    
    ;print,i,B_p,rho_0,T_0
  
    ; cal(I) integrals
    int_rho = qpint1d('x*x*(1-P(0)*(1-alog(1+x)/x))^(P(1))', $
                      r_cool/nfw_dm.r_s, nfw_dm.r200/nfw_dm.r_s, [suto.B_p,n],/expr )
    int_ene = qpint1d('x*x*(1-P(0)*(1-alog(1+x)/x))^(P(1))', $
                      r_cool/nfw_dm.r_s, nfw_dm.r200/nfw_dm.r_s, [suto.B_p,n+1],/expr )

    ; enforce mass constraint -> new rho0
    rho_0_new = mass_hot_codeunits / (4*!pi*nfw_dm.r_s^3.0 * int_rho)
    
    ; enforce energy constraint -> new T0
    ;Eh_H = f_shock * (-0.5) * nfw_dm.bindEnergy * (mass_hot_codeunits/nfw_dm.mass_codeunits)
    ;T_0_new = units.mu * units.mass_proton * (Eh_H*units.UnitEnergy_in_cgs) / $
    ;  (int_ene * 6 * !pi * units.boltzmann * rho_0_new * nfw_dm.r_s^3.0 * units.UnitMass_in_g)
    
    ; or: central temperature is fixed
    T_0_new = T_0
    
    ; check convergence and update
    cur_err1 = abs((T_0_new-T_0)/T_0)
    cur_err2 = abs((rho_0_new-rho_0)/rho_0)
    if cur_err1 lt errTol and cur_err2 lt errTol then break
    
    rho_0 = float( rho_0 + (rho_0_new-rho_0) * iterFac )
    T_0   = float( T_0   + (T_0_new-T_0)     * iterFac )
  
  endfor

  ;print,'final2',alog10(T_0),alog10(rho_0),suto.B_p
  
  ; construct final density/temperature profiles
  suto = suto_model(T_0,rho_0,n,nfw_dm)
  
  ; timescales
  ts = nfw_timescales_suto(suto=suto, nfw_dm=nfw_dm, tables=tables)
  
  return, {B:suto.B, B_p:suto.B_p, temp_gas:suto.temp_gas, rho_gas:suto.rho_gas, $
           T_0:suto.T_0, rho_0:suto.rho_0, n:suto.n ,$
           dynTime:ts.dynTime, coolTime:ts.coolTime, hubbleTime:ts.hubbleTime, $
           r_cool:ts.r_cool, r_cool_h:ts.r_cool_h, r_cool_iso:ts.r_cool_iso, $
           coolTime_iso:ts.coolTime_iso, rho_gas_iso:suto.rho_gas_iso}  
  
end
  
; nfw_gas_fit_func(): fitting from simulations helper

function nfw_gas_fit_func, X, P, B_in=B_in
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; P = [log T_0, log rho_0, n]
  T_0   = P[0]
  rho_0 = P[1]
  n     = P[2]

  ; finish B and compute B_p
  B = B_in / T_0 
  B_p = B/(n+1)
  
  ; polytropic
  X_single = X[0:n_elements(X)/2-1]
  eps_x = 1 - B_p * ( 1 - alog(1 + X_single) / X_single )
  
  temp_gas = T_0   * eps_x
  rho_gas  = rho_0 * eps_x^n
  
  Y = [rho_gas,temp_gas]

  return, Y
end

; nfw_gas_fit(): fit NFW (poly) to radially binned gas density profile from the simulations
; input/output: code units as in simulation (comoving/h)
; note: n->inf is isothermal, and halos (dens,temp) cannot be simultaneously fit well by nfw_poly
;       de-weighting the temperature results in a correct density fit but the temperature is too low

function nfw_gas_fit, rad_frac=rad_frac, dens=dens, temp=temp, x=x, nfw_dm=nfw_dm, tables=tables

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits(redshift=nfw_dm.redshift)
  
  err_y = replicate(1.0,n_elements(rad_frac))
  
  ; fit density/temperature: prepare inputs
  xx    = [rad_frac,rad_frac]
  y     = [dens*1e10,10.0^temp]
  err_y = [err_y/1e2,err_y] ; weight density more?

  ; params: (T_0,rho_0,n)
  p0 = [3*nfw_dm.Tvir_rVir            ,$ ; T_0
        0.1*nfw_dm.rho_dm[0]*1e10     ,$ ; rho_0
        20.0                            ] ; n
  
  parinfo = replicate({limited:[1,1], limits:[0.0,0.0], fixed:0, step:0.0}, 3)

  ; limits
  parinfo[0].limits  = [nfw_dm.Tvir_rVir,10*nfw_dm.Tvir_rVir]
  parinfo[1].limits  = [0.0,nfw_dm.rho_dm[0]*1e10] ; non-negative central density
  parinfo[2].limits  = [1.5,3000.0]
  
  ; fixed at initial guess?
  parinfo[0].fixed = 0
  parinfo[1].fixed = 0
  parinfo[2].fixed = 1
  
  ; stepsizes
  ;parinfo[0].step = 100.0 ; K
  ;parinfo[1].step = 1e-8 ; rho
  ;parinfo[2].step = 0.1 ; polytropic index
  
  ; pre-compute most of B
  B_in = float( ( 4 * !pi * units.Gravity * units.mu * units.mass_proton * $
               units.rhoCrit_z*nfw_dm.delta_c * nfw_dm.r_s^2.0 * $
               units.UnitMass_in_g / units.UnitLength_in_cm ) / units.boltzmann )
  
  p_poly = mpfitfun('nfw_gas_fit_func', xx, y, err_y, p0, functargs={B_in:B_in}, parinfo=parinfo, /quiet)
  f_poly = nfw_gas_fit_func([x,x], p_poly, B_in=B_in)
  
  ; split f_poly into density and temperature
  f_dens = f_poly[0:n_elements(f_poly)/2-1] / 1e10
  f_temp = alog10( f_poly[n_elements(f_poly)/2:n_elements(f_poly)-1] )
  
  ; add in timescales etc
  suto = suto_model(p_poly[0],p_poly[1],p_poly[2],nfw_dm)
  ts = nfw_timescales_suto(suto=suto, nfw_dm=nfw_dm, tables=tables)

  return, { p_poly:p_poly, f_dens:f_dens, f_temp:f_temp, ts:ts }

end

; SIS
; ---

; SIS profile given (M,z)
; inputs: rad_frac in fraction of rvir, mass in code units

function sis_profile, rad_frac, mass=mass_codeunits, redshift=redshift

  compile_opt idl2, hidden, strictarr, strictarrsubs
  if ~keyword_set(mass_codeunits) or ~n_elements(redshift) then message,'error'
  units = getUnits(redshift=redshift)
  
  ; calculate virial radius (r200), Monaco+ (2007) eqn 5
  r200 = (mass_codeunits * units.G / 100 / units.H_z^2.0)^(1.0/3.0)
  ;r200 = (3*mass_codeunits / 800 / units.rhoCrit_z / !pi)^(1.0/3.0) ;equivalent
  
  ; convert input rad to kpc (code units)
  rad = rad_frac * r200
  
  ; density profile
  rho_DM = mass_codeunits / (4 * !pi * rad_frac^2.0 * r200^3.0)
  
  ; enclosed mass
  Menc_DM   = rad * mass_codeunits / r200
  Menc_rVir = mass_codeunits
  
  ; circular velocity (constant)
  Vcirc_DM   = sqrt(units.G * Menc_DM / rad)
  Vcirc_rVir = sqrt(units.G * Menc_rVir / r200)
  
  ; virial temperature
  Tvir_DM   = 35.9 * Vcirc_DM^2.0
  Tvir_rVir = 35.9 * Vcirc_rVir^2.0
  
  return, {rad:rad, r200:r200, rho_DM:rho_DM, mass_codeunits:mass_codeunits, redshift:redshift, $
           Menc_DM:Menc_DM, Menc_rvir:Menc_rvir, $
           Vcirc_DM:Vcirc_DM, Vcirc_rvir:Vcirc_rvir, $
           Tvir_DM:Tvir_DM, Tvir_rVir:tVir_rVir}

end

; gas density and temperature profile
; inputs: mass_hot is the "hot halo gas mass" (e.g. f_b * M_DM - M_op)
;         where M_op = M_cold + M_stars + M_ejecta + ... (other phases already removed from hot, model dependent)
; CROTON/KANG

function sis_gas_profile, mass_hot=mass_hot_codeunits, sis_dm=sis_dm, tables=tables

  compile_opt idl2, hidden, strictarr, strictarrsubs
  if ~keyword_set(mass_hot_codeunits) or ~keyword_set(sis_dm) then message,'error'
  if ~n_elements(tables) then tables = interpLambdaSD93()
  units = getUnits(redshift=sis_dm.redshift)

  ; get cooling rate for T(r) (constant) in cgs (erg/s/cm^3)
  lambdaNet = 10.0^(interpLambdaSD93(Z=0.0,logT=alog10(sis_dm.Tvir_rVir),tables=tables))
  
  ; gas density profile
  rho_0 = mass_hot_codeunits / (4*!pi*sis_dm.r200)
  
  rho_gas = rho_0 / sis_dm.rad^2.0

  ; gas temperature profile (constant/isothermal)
  temp_gas = replicate(sis_dm.Tvir_rVir,n_elements(sis_dm.rad))
  
  ; timescales in Gyr
  ;alpha = 3 * units.mu / 2
  alpha = 3.56 ; LU
  
  COF_cgs  = float(alpha * units.mass_proton * units.boltzmann * sis_dm.Tvir_rVir / lambdaNet)
  
  coolTime   = COF_cgs / (rho_gas*units.UnitDensity_in_cgs) / units.s_in_Gyr
  
  dynTime      = sis_dm.rad / (sis_dm.Vcirc_DM * units.kmS_in_kpcYr) / 1e9
  dynTime_halo = sis_dm.r200 / (sis_dm.Vcirc_rVir * units.kmS_in_kpcYr) / 1e9
  
  hubbleTime = 1/units.H_z ; Gyr
  hubbleTime = redshiftToAgeFlat(sis_dm.redshift)
  
  ; rcool_croton (where coolTime=dynTime_halo)
  r_cool_croton = sqrt(rho_0 * units.UnitDensity_in_cgs * dynTime_halo * units.s_in_Gyr / COF_cgs)
  
  ; rcool_hubble (where coolTime=hubbleTime), used for Kang model
  r_cool_hubble = sqrt(rho_0 * units.UnitDensity_in_cgs * hubbleTime * units.s_in_Gyr / COF_cgs)
  
  ; rcool (where coolTime=dynTime)
  r_cool = rootFindVector(coolTime,dynTime,xx=sis_dm.rad)
    
  return, {rho_gas:rho_gas, temp_gas:temp_gas, coolTime:coolTime, dynTime:dynTime, hubbleTime:hubbleTime, $
           r_cool_croton:r_cool_croton, dynTime_halo:dynTime_halo, $ ; used for Croton model
           r_cool:r_cool, r_cool_hubble:r_cool_hubble, lambdaNet:lambdaNet}
    
end

; fit SIS to radially binned gas density profile from the simulations
; input/output: code units as in simulation (comoving/h)

function sis_gas_fit, rad_frac=rad_frac, dens=dens, temp=temp, x=x

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  err_y = replicate(1.0,n_elements(rad_frac))
  
  ; params: P(0) = mass_hot_codeunits_physical
  ;         P(1) = r200_codeunits_physical
  ;expr = 'P[0] / ( 4 * !pi * x^2.0 * P[1]^3.0 )'
  ;p0 = [1.0,10.0]
  
  ; param: P(0) = rho_0 central density (code units)
  expr = 'P[0] / x^2.0'
  p0 = max(dens)
  
  p_dens = mpfitexpr(expr, rad_frac, dens, err_y, p0, /quiet)
  f_dens = mpevalexpr(expr, x, p_dens)
  
  ; fit temp
  expr = 'P[0]'
  p0 = [5.0]
  
  p_temp = mpfitexpr(expr, rad_frac, temp, err_y, p0, /quiet)
  f_temp = replicate(p_temp[0],n_elements(x))
  
  return, { p_dens:p_dens, f_dens:f_dens, p_temp:p_temp, f_temp:f_temp }

end
