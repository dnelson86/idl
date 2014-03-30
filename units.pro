; units.pro
; units and unit conversions (for Kpc, 10^10 Msun, Gyr unit set)
; dnelson sep.2013
; note: for cosmological+comoving simulations, all masses and lengths are 1/h
;       and all lengths are comoving (multiply by scalefac for physical)
;       and all times derived from Timebase_interval are similar (divide by scalefac for normal)

; getUnits(): return a structure of useful units

function getUnits, redshift=redshift

  units = { units,                                             $
  
            ; units (from parameter file)
            UnitLength_in_cm         : double(3.085678e21)    ,$ ; 1.0 kpc
            UnitMass_in_g            : 1.989*double(10.0)^43  ,$ ; 1.0e10 solar masses
            UnitVelocity_in_cm_per_s : double(1.0e5)          ,$ ; 1 km/sec
            
            ; derived units
            UnitTime_in_s       : 0.0D                        ,$
            UnitDensity_in_cgs  : 0.0D                        ,$
            UnitPressure_in_cgs : 0.0D                        ,$
            UnitEnergy_in_cgs   : 0.0D                        ,$
            UnitTemp_in_cgs     : 0.0D                        ,$
            
            ; non-cgs units
            UnitMass_in_Msun    : 0.0D                        ,$
            UnitTime_in_yr      : 0.0D                        ,$
            
            ; constants
            boltzmann         : double(1.380650e-16)          ,$ ; cgs (erg/K)
            mass_proton       : double(1.672622e-24)          ,$ ; cgs
            hydrogen_massfrac : 0.76                          ,$ ; XH (solar)
            helium_massfrac   : 0.25                          ,$ ; Y (solar)
            mu                : 0.6                           ,$ ; ionized primordial (e.g. hot halo gas)
            HubbleParam       : 0.7                           ,$ ; little h (All.HubbleParam), e.g. H0 in 100 km/s/Mpc
            Gravity           : 6.673e-8                      ,$ ; G in cgs, cm^3/g/s^2
            H0_kmsMpc         : 70.0                          ,$ ; km/s/Mpc
            
            ; derived constants (in code units)
            H0      : 0.0                                     ,$ ; km/s/kpc
            G       : 0.0                                     ,$ ; kpc (km/s)^2 / 1e10 msun
            rhoCrit : 0.0                                     ,$ ; 1e10 msun / kpc^3
            
            ; cosmology parameters (valid for ComparisonProject/Illustris, z=0 values)
            omega_m : 0.27                                    ,$
            omega_L : 0.73                                    ,$
            omega_k : 0.0                                     ,$
            omega_b : 0.044                                   ,$
            
            ; derived cosmology parameters
            f_b : 0.0                                         ,$
            
            ; all previous were z=0, this section for different redshift (if specified) (code units)
            H2_z_fact : !values.f_nan                         ,$
            H_z       : !values.f_nan                         ,$
            rhoCrit_z : !values.f_nan                         ,$
            
            ; color list
            colors : strarr(18)                               ,$
            
            ; unit conversions
            s_in_yr   : 3.155693e7                            ,$
            s_in_Myr  : 3.155693e13                           ,$
            s_in_Gyr  : 3.155693e16                           ,$
            Msun_in_g : 1.98892*double(10.0)^33               ,$
            pc_in_cm  : 3.085680e18                           ,$
            Mpc_in_cm : 3.085680e24                           ,$
            kpc_in_km : 3.085680e16                           ,$
            
            ; derived unit conversions
            kmS_in_kpcYr  : 0.0                               ,$
            kmS_in_kpcGyr : 0.0                                $
      }
      
  ; derived units
  units.UnitTime_in_s       = units.UnitLength_in_cm / units.UnitVelocity_in_cm_per_s
  units.UnitDensity_in_cgs  = units.UnitMass_in_g / units.UnitLength_in_cm^3.0
  units.UnitPressure_in_cgs = units.UnitMass_in_g / units.UnitLength_in_cm / units.UnitTime_in_s^2.0
  units.UnitEnergy_in_cgs   = units.UnitMass_in_g * units.UnitLength_in_cm^2.0 / units.UnitTime_in_s^2.0
  units.UnitTemp_in_cgs     = units.UnitEnergy_in_cgs / units.UnitMass_in_g
  
  ; non-cgs units
  units.UnitMass_in_Msun    = units.UnitMass_in_g / units.Msun_in_g
  units.UnitTime_in_yr      = units.UnitTime_in_s / units.s_in_yr
  
  ; derived constants (in code units)
  units.H0 = units.HubbleParam * 100 * 1e5 / (units.Mpc_in_cm) / $
             units.UnitVelocity_in_cm_per_s * units.UnitLength_in_cm
  units.G  = units.Gravity / units.UnitLength_in_cm^3.0 * units.UnitMass_in_g * units.UnitTime_in_s^2.0
  
  units.rhoCrit = 3.0 * units.H0^2.0 / (8.0*!pi*units.G) ;code, z=0

  ; derived cosmology parameters
  units.f_b = units.omega_b / units.omega_m
  
  ; redshift dependent values (code units)
  if n_elements(redshift) then begin
    units.H2_z_fact = ( units.omega_m*(1+redshift)^3.0 + units.omega_L + units.omega_k*(1+redshift)^2.0 )
    units.H_z       = units.H0 * sqrt(units.H2_z_fact)
    units.rhoCrit_z = units.rhoCrit * units.H2_z_fact
  endif
  
  ; derived unit conversions
  units.kmS_in_kpcYr  = units.s_in_Myr / units.kpc_in_km / 1e6 ; Myr->yr
  units.kmS_in_kpcGyr = units.s_in_Myr / units.kpc_in_km * 1e3 ; Myr->Gyr

  ; color list
  units.colors = ['black','blue','green','red','cyan','magenta','gray','orange', $
                  'brown','chartreuse','violet','papaya','aquamarine', $
                  'firebrick', 'rosy brown', 'gold', 'forest green', 'slate blue']

  return, units
end

; meanmolwt(): from Monaco+ (2007) eqn 14, for hot halo gas

function meanmolwt, Y=Y, Z=Z ; helium fraction (0.25) and metallicity (non-log metal mass/total mass)
  mu = 4.0 / (8 - 5*Y - 6*Z)
  return,mu
end

; codeMassToVirTemp(): convert halo mass (in code units) to virial temperature at specified redshift
; Barkana & Loeb (2001) eqn.26

function codeMassToVirTemp, mass, redshift=redshift, sP=sP, meanmolwt=meanmolwt, log=log

  units = getUnits()
  
  if n_elements(redshift) eq 0 then redshift = sP.redshift
  if redshift eq -1 then begin
    h = loadSnapshotHeader(sP=sP)
    redshift = 1/h.time-1
    print,redshift
  endif
  if redshift eq -1 or redshift lt 0.0 or redshift gt 20.0 then message,'Error: Still do not know redshift.'

  ; mean molecular weight default (valid for ComparisonProject at ionized T>10^4 K)
  if not keyword_set(meanmolwt) then meanmolwt = 0.6

  ; mass to msun
  mass_msun = mass * float(units.UnitMass_in_g / units.Msun_in_g)
  
  little_h = 1.0 ; do not multiply by h since mass_msun is already over h
  
  omega_m_z = units.omega_m * (1+redshift)^3.0 / $
              ( units.omega_m*(1+redshift)^3.0 + units.omega_L + units.omega_k*(1+redshift)^2.0 )
  
  Delta_c = 18*!pi^2 + 82*(omega_m_z-1.0) - 39*(omega_m_z-1.0)^2.0

  Tvir = 1.98e4 * (meanmolwt/0.6) * (mass_msun/1e8*little_h)^(2.0/3.0) * $
         (units.omega_m/omega_m_z * Delta_c / 18.0 / !pi^2.0)^(1.0/3.0) * $
         (1.0 + redshift)/10.0 ;K
         
  if keyword_set(log) then begin
    w = where(Tvir ne 0.0,count)
    if count gt 0 then Tvir[w] = alog10(Tvir[w])
  endif

  return, Tvir
  
end

; logMsunToVirTemp(): convert halo mass (in log msun) to virial temperature at specified redshift

function logMsunToVirTemp, mass, redshift=redshift, sP=sP, meanmolwt=meanmolwt

  units = getUnits()
  print,'Warning: Using this you should have no little h in the log msun.'
  
  if not keyword_set(redshift) then redshift = sP.redshift
  if redshift eq -1 then stop

  ; mean molecular weight default (valid for ComparisonProject at ionized T>10^4 K)
  if not keyword_set(meanmolwt) then meanmolwt = 0.6

  ; mass to msun
  mass_msun = 10.0^mass
  
  ; cosmo
  little_h  = 0.7
  
  omega_m_z = units.omega_m * (1+redshift)^3.0 / $
              ( units.omega_m*(1+redshift)^3.0 + units.omega_L + units.omega_k*(1+redshift)^2.0 )
  
  Delta_c = 18*!pi^2 + 82*(omega_m_z-1.0) - 39*(omega_m_z-1.0)^2.0

  Tvir = 1.98e4 * (meanmolwt/0.6) * (mass_msun/1e8*little_h)^(2.0/3.0) * $
         (units.omega_m/omega_m_z * Delta_c / 18.0 / !pi^2.0)^(1.0/3.0) * $
         (1.0 + redshift)/10.0 ;K
         
  return, Tvir

end

; codeMassToLogMsun(): convert mass in code units to log(msun)

function codeMassToLogMsun, mass
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  mass_msun = mass * float(units.UnitMass_in_g / units.Msun_in_g)
  
  ; log of nonzero
  w = where(mass_msun eq 0.0,count)
  if (count ne 0) then $
    mass_msun[w] = 1.0
  
  return,alog10(mass_msun)
end

function logMsunToCodeMass, mass
  return,10.0^mass / 1e10
end

; codeTempToLogK(): convert temperature in code units (e.g. tracer temp output) to log(kelvin)

function codeTempToLogK, temp
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  temp_k = temp * float(units.UnitTemp_in_cgs)

  ; log of nonzero
  w = where(temp_k eq 0.0,count)
  if (count ne 0) then $
    temp_k[w] = 1.0
  
  return,alog10(temp_k)
end

; codeDensToPhys(): convert comoving->Physical and add little_h factors

function codeDensToPhys, dens, sP=sP, scalefac=scalefac, cgs=cgs
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  if n_elements(scalefac) eq 0 and n_elements(sP) eq 0 then message,'error'
  if n_elements(scalefac) gt 0 and n_elements(sP) gt 0 then message,'error'
  if n_elements(sP) gt 0 then scalefac = 1.0/(1 + sP.redshift)

  dens_phys = dens * units.HubbleParam * units.HubbleParam / scalefac^3.0

  if keyword_set(cgs) then dens_phys *= units.UnitDensity_in_cgs
  
  return, dens_phys
  
end

; convertUtoTemp(): convert u,nelec pair in code units to temperature in Kelvin or log(Kelvin)

function convertUtoTemp, u, nelec, gamma=gamma, hmassfrac=hmassfrac, log=log
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; adiabatic index and hydrogen mass fraction defaults (valid for ComparisonProject)
  if not keyword_set(gamma)     then gamma = 5.0/3.0
  if not keyword_set(hmassfrac) then hmassfrac = units.hydrogen_massfrac
  
  ; calculate mean molecular weight
  meanmolwt = 4.0/(1.0 + 3.0 * hmassfrac + 4.0* hmassfrac * nelec) * units.mass_proton

  ; calculate temperature (K)
  temp = (gamma-1.0) * u / units.boltzmann * units.UnitEnergy_in_cgs / units.UnitMass_in_g * meanmolwt
  
  ; convert to log(K) if requested
  if keyword_set(log) then begin
    w = where(temp eq 0.0,count)
    if (count ne 0) then temp[w] = 1.0
    temp = alog10( temporary(temp) )
  endif
  
  return, float(temp)
end

; convertTempToU(): convert temperature in log(K) to u in code units

function convertTempToU, logTemp, gamma=gamma, hmassfrac=hmassfrac, log=log
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  if max(logTemp) gt 10.0 then message,'Error: input temp probably not in log.'
  
  ; adiabatic index and hydrogen mass fraction defaults (valid for ComparisonProject)
  if not keyword_set(gamma)     then gamma = 5.0/3.0
  if not keyword_set(hmassfrac) then hmassfrac = units.hydrogen_massfrac

  meanmolwt = 0.6 * units.mass_proton ; ionized, T > 10^4 K
  
  ;temp = (gamma-1.0) * u / units.boltzmann * units.UnitEnergy_in_cgs / units.UnitMass_in_g * meanmolwt
  u = 10.0^logTemp * units.boltzmann * units.UnitMass_in_g / (units.UnitEnergy_in_cgs * meanmolwt * (gamma-1.0))
  
  ; convert to log(u) if requested
  if keyword_set(log) then begin
    w = where(u eq 0.0,count)
    if (count ne 0) then u[w] = 1.0
    u = alog10( temporary(u) )
  endif
  
  return, float(u)
  
end

; codeMassToVirEnt(): given a halo mass, return a S200 (e.g. Pvir/rho_200crit^gamma)

function codeMassToVirEnt, mass, redshift=redshift, sP=sP, log=log
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  if n_elements(mass) eq 0 then message,'Error'
  
  if n_elements(redshift) eq 0 then redshift = sP.redshift
  
  virTempLog = codeMassToVirTemp(mass,redshift=redshift,sP=sP,/log)
  virU = convertTempToU(virTempLog)
  r200crit = critRatioToCode(200.0,redshift=redshift)
  
  s200 = calcEntropyCGS(virU, r200crit, sP=sP, log=log)
  
  return, s200

end

; convertCoolingRatetoCGS():

function convertCoolingRatetoCGS, coolrate, h=h
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; default little h
  if not keyword_set(h) then h = 0.7
  
  ; convert code units (du/dt) to erg/s/g (cgs)
  coolrate_cgs = coolrate * units.UnitEnergy_in_cgs * units.UnitTime_in_s^(-1.0) * $
                 units.UnitMass_in_g^(-1.0) * h
                 
  return, coolrate_cgs
end

; convertTracerEntToCGS(): fix cosmological/unit system in TRACER_MC.MaxEnt

function convertTracerEntToCGS, ent, gamma=gamma, log=log, sP=sP

  forward_function snapNumToRedshift
  units = getUnits()
  
  ; adiabatic index default (valid for ComparisonProject)
  if not keyword_set(gamma) then gamma = 5.0/3.0
  
  atime = snapNumToRedshift(sP=sP,/time)
  a3inv = 1.0 / (atime*atime*atime)
  
  ; NOTE: dens=dens*a3inv but in the tracers only converted in dens^gamma not in the pressure
  ; have to make this adjustment in loading tracers
  ; for SFR, for gas and tracers, Pressure = GAMMA_MINUS1 * localSphP[i].Density * localSphP[i].Utherm;
  ; for TRACER_MC, EntMax = SphP.Pressure / pow(SphP.Density * All.cf_a3inv, GAMMA);
  
  ; fix Pressure
  ent *= a3inv * float(units.UnitPressure_in_cgs / units.boltzmann)
  
  ; fix Density
  ent /= float(units.UnitDensity_in_cgs/units.mass_proton)^gamma
  
  ; convert to log(entropy) if requested
  if keyword_set(log) then begin
    w = where(ent eq 0.0,count)
    if (count ne 0) then ent[w] = 1.0
    ent = alog10(ent)
  endif
  
  return, ent ; [K cm^2]
end

; calcEntropyCGS(): calculate entropy as P/rho^gamma (rho is converted from comoving to physical)

function calcEntropyCGS, u, dens, gamma=gamma, log=log, sP=sP

  forward_function snapNumToRedshift
  units = getUnits()
  
  ; adiabatic index default (valid for ComparisonProject)
  if not keyword_set(gamma) then gamma = 5.0/3.0
  
  atime = snapNumToRedshift(sP=sP,/time)
  a3inv = 1.0 / (atime*atime*atime)

  ; cosmological and unit system conversions
  pressure = (gamma-1.0) * u * dens * a3inv * float(units.UnitPressure_in_cgs / units.boltzmann) ; [K/cm^3]
  entropy  = pressure / ( (dens * float(units.UnitDensity_in_cgs/units.mass_proton)*a3inv)^gamma ) ; [K cm^2]
  
  ; convert to log(entropy) if requested
  if keyword_set(log) then begin
    w = where(entropy eq 0.0,count)
    if (count ne 0) then entropy[w] = 1.0
    entropy = alog10(entropy)
  endif  
  
  return, entropy
end

; calcPressureCGS(): calculate pressure as (gamma-1)*u*rho

function calcPressureCGS, u, dens, gamma=gamma, log=log, sP=sP

  forward_function snapNumToRedshift
  units = getUnits()
  
  ; adiabatic index default (valid for ComparisonProject)
  if not keyword_set(gamma) then gamma = 5.0/3.0
  
  atime = snapNumToRedshift(sP=sP,/time)
  a3inv = 1.0 / (atime*atime*atime)
  
  pressure = (gamma-1.0) * u * (dens*a3inv)
  
  ; convert to CGS = 1 barye (ba) = 1 dyn/cm^2 = 0.1 Pa = 0.1 N/m^2 = 0.1 kg/m/s^2
  ; and divide by boltzmann's constant -> [K/cm^3]
  pressure *= float(units.UnitPressure_in_cgs / units.boltzmann)
  
  ; convert to log(pressure) if requested
  if keyword_set(log) then begin
    w = where(pressure eq 0.0,count)
    if (count ne 0) then pressure[w] = 1.0
    pressure = alog10(pressure)
  endif 
  
  return, pressure

end

; rhoRatioToCrit(): normalize density by the critical -baryon- density at some redshift

function rhoRatioToCrit, rho, sP=sP, redshift=redshift

  if n_elements(redshift) eq 0 and n_elements(sP) eq 0 then message,'error'
  if n_elements(redshift) eq 1 and n_elements(sP) eq 1 then message,'error'
  
  if n_elements(redshift) gt 0 then zzz = redshift
  if n_elements(sP) gt 0 then zzz = sP.redshift
  
  units = getUnits(redshift=zzz)
  
  rho_b = units.omega_b * units.rhoCrit_z
  
  return, rho/rho_b

end

; critRatioToCode(): convert a ratio of the critical density at some redshift to a code density

function critRatioToCode, ratioToCrit, redshift=redshift

  units = getUnits(redshift=redshift)
  if n_elements(redshift) eq 0 then message,'specify redshift'

  codeRho = ratioToCrit * units.rhoCrit_z
  
  return,codeRho

end

; critBaryonRatioToCode(): convert a ratio of the critical baryon density at some redshift to a code density

function critBaryonRatioToCode, ratioToCritB, redshift=redshift

  units = getUnits(redshift=redshift)
  if n_elements(redshift) eq 0 then message,'specify redshift'

  codeRho = ratioToCritB * units.omega_b * units.rhoCrit_z
  
  return,codeRho

end

; redshiftToAge(): convert redshift to age of the universe (approximate)
function dtdz, z, lambda0 = lambda0, q0 = q0
  term1 = (1.0d + z)
  term2 = 2.0d * (q0 + lambda0) * z + 1.0d - lambda0
  term3 = (1.0d + z) * (1.0d +z)
  return, 1.0 / (term1 * sqrt(term2 * term3 + lambda0))
end
   
function redshiftToAge, z

  units = getUnits()
  
  ; config
  zform = 1000.0
  H0 = 70.0
  k = 0.0
  Omega_m = 0.27
  Lambda0 = 0.73
  q0 = -0.55
  
  ; arrays
  nz  = N_elements(z)
  age = z * 0.0
  
  ; integrate with qsimp
  for i= 0L, nz-1 do begin
    if (z[i] ge zform) then age_z = 0 else $
        qsimp,'dtdz', z[i], zform, age_z, q0 = q0, lambda0 = lambda0
    age[i] = age_z
  endfor

  return, age * 3.085678e+19 / 3.15567e+7 / H0 / 1e9 ;Gyr

end

; redshiftToAgeFlat(): analytical formula from Peebles, p. 317, eq. 13.2. 

function redshiftToAgeFlat, z

  units = getUnits()
  
  age = 2*asinh(sqrt( (1-units.omega_m)/units.omega_m ) * (1+z)^(-3.0/2.0)) / $
             (units.H0_kmsMpc * 3 * sqrt(1-units.omega_m))
  
  return, age * 3.085678e+19 / 3.15567e+7 / 1e9 ;Gyr

end
