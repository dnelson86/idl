; units.pro
; units and unit conversions (for Kpc, 10^10 Msun, Gyr unit set)
; dnelson mar.2012

; getUnits(): return a structure of useful units

function getUnits

  Hubble  = 1.0      ;H0 in 100km/s/Mpc
  Gravity = 6.673e-8 ;G in cgs, cm^3/g/s^2

  units = { units,                   $
  
            ; units (from parameter file)
            UnitLength_in_cm         : double(3.085678e21)    ,$;  1.0 kpc
            UnitMass_in_g            : 1.989*double(10.0)^43  ,$;  1.0e10 solar masses
            UnitVelocity_in_cm_per_s : double(1.0e5)          ,$;  1 km/sec
            
            ; derived units
            UnitTime_in_s       : 0.0D                        ,$
            UnitDensity_in_cgs  : 0.0D                        ,$
            UnitPressure_in_cgs : 0.0D                        ,$
            UnitEnergy_in_cgs   : 0.0D                        ,$
            UnitTemp_in_cgs     : 0.0D                        ,$
            
            ; non-cgs units
            UnitMass_in_Msun    : 0.0D                        ,$
            
            ; constants
            boltzmann   : double(1.380650e-16)                ,$ ;cgs
            mass_proton : double(1.672622e-24)                ,$ ;cgs
            
            ; derived constants
            H0      : 0.0D                                    ,$
            G       : 0.0D                                    ,$
            rhoCrit : 0.0D                                    ,$
            
            ; color list
            colors : strarr(18)                               ,$
            
            ; unit conversions
            s_in_Myr  : 3.155693e13                           ,$
            s_in_Gyr  : 3.155693e16                           ,$
            Msun_in_g : 1.98892*double(10.0)^33               ,$
            pc_in_cm  : 3.085680e18                           ,$
            Mpc_in_cm : 3.085680e24                           ,$
            kpc_in_km : 3.085680e16                           ,$
            
            ; derived unit conversions
            kmS_in_kpcYr : 0.0                                 $
      }
      
  ; derived units
  units.UnitTime_in_s       = units.UnitLength_in_cm / units.UnitVelocity_in_cm_per_s
  units.UnitDensity_in_cgs  = units.UnitMass_in_g / units.UnitLength_in_cm^3.0
  units.UnitPressure_in_cgs = units.UnitMass_in_g / units.UnitLength_in_cm / units.UnitTime_in_s^2.0
  units.UnitEnergy_in_cgs   = units.UnitMass_in_g * units.UnitLength_in_cm^2.0 / units.UnitTime_in_s^2.0
  units.UnitTemp_in_cgs     = units.UnitEnergy_in_cgs / units.UnitMass_in_g
  
  ; non-cgs units
  units.UnitMass_in_Msun    = units.UnitMass_in_g / units.Msun_in_g
  
  ; derived constants (in code units)
  units.H0 = Hubble * 100 * 1e5 / (units.Mpc_in_cm) / $
             units.UnitVelocity_in_cm_per_s * units.UnitLength_in_cm
  units.G  = Gravity / units.UnitLength_in_cm^3.0 * units.UnitMass_in_g * units.UnitTime_in_s^2.0
  
  units.rhoCrit = 3.0 * units.H0^2.0 / (8.0*!pi*units.G) ;code, z=0

  ; derived unit conversions
  units.kmS_in_kpcYr = units.s_in_Myr / units.kpc_in_km / 1e6 ; Myr->yr

  ; color list
  units.colors = ['black','blue','green','red','cyan','magenta','gray','orange', $
                  'brown','chartreuse','violet','papaya','aquamarine', $
                  'firebrick', 'rosy brown', 'gold', 'forest green', 'slate blue']

  return, units
end

; codeMassToVirTemp(): convert halo mass (in code units) to virial temperature at specified redshift

function codeMassToVirTemp, mass, redshift=redshift, sP=sP, meanmolwt=meanmolwt

  units = getUnits()
  
  if not keyword_set(redshift) then redshift = sP.redshift
  if redshift eq -1 then stop

  ; mean molecular weight default (valid for ComparisonProject at ionized T>10^4 K)
  if not keyword_set(meanmolwt) then meanmolwt = 0.6

  ; mass to msun
  mass_msun = mass * (units.UnitMass_in_g / units.Msun_in_g)
  
  ; cosmo
  omega_m   = 0.27
  omega_L   = 0.73
  omega_k   = 0.0
  little_h  = 0.7
  
  omega_m_z = omega_m * (1+redshift)^3.0 / $
              ( omega_m*(1+redshift)^3.0 + omega_L + omega_k*(1+redshift)^2.0 )
  
  Delta_c = 18*!pi^2 + 82*(omega_m_z-1.0) - 39*(omega_m_z-1.0)^2.0

  Tvir = 1.98e4 * (meanmolwt/0.6) * (mass_msun/1e8*little_h)^(2.0/3.0) * $
         (omega_m/omega_m_z * Delta_c / 18.0 / !pi^2.0)^(1.0/3.0) * $
         (1.0 + redshift)/10.0 ;K

  return, Tvir
  
end

; codeMassToLogMsun(): convert mass in code units to log(msun)

function codeMassToLogMsun, mass

  units = getUnits()
  
  mass_msun = mass * (units.UnitMass_in_g / units.Msun_in_g)
  
  ; log of nonzero
  w = where(mass_msun eq 0.0,count)
  if (count ne 0) then $
    mass_msun[w] = 1.0
  
  return,alog10(mass_msun)
end

; codeTempToLogK(): convert temperature in code units (e.g. tracer temp output) to log(kelvin)

function codeTempToLogK, temp

  units = getUnits()
  
  temp_k = temp * units.UnitTemp_in_cgs

  ; log of nonzero
  w = where(temp_k eq 0.0,count)
  if (count ne 0) then $
    temp_k[w] = 1.0
  
  return,alog10(temp_k)
end

; convertUtoTemp(): convert u,nelec pair in code units to temperature in Kelvin or log(Kelvin)

function convertUtoTemp, u, nelec, gamma=gamma, hmassfrac=hmassfrac, log=log

  units = getUnits()
  
  ; adiabatic index and hydrogen mass fraction defaults (valid for ComparisonProject)
  if not keyword_set(gamma)     then gamma = 5.0/3.0
  if not keyword_set(hmassfrac) then hmassfrac = 0.76
  
  ; calculate mean molecular weight
  meanmolwt = 4.0/(1.0 + 3.0 * hmassfrac + 4.0* hmassfrac * nelec) * units.mass_proton

  ; calculate temperature (K)
  temp = (gamma-1.0) * u / units.boltzmann * units.UnitEnergy_in_cgs / units.UnitMass_in_g * meanmolwt
  
  ; convert to log(K) if requested
  if keyword_set(log) then begin
    w = where(temp eq 0.0,count)
    if (count ne 0) then temp[w] = 1.0
    temp = alog10(temp)
  endif
  
  return, float(temp)
end

; convertCoolingRatetoCGS():

function convertCoolingRatetoCGS, coolrate, h=h

  units = getUnits()
  
  ; default little h
  if not keyword_set(h) then h = 0.7
  
  ; convert code units (du/dt) to erg/s/g (cgs)
  coolrate_cgs = coolrate * units.UnitEnergy_in_cgs * units.UnitTime_in_s^(-1.0) * $
                 units.UnitMass_in_g^(-1.0) * h
                 
  return, coolrate_cgs
end

; calcEntropy(): calculate entropy as P/rho^gamma (rho is converted from comoving to physical)

function calcEntropy, u, dens, gamma=gamma, sP=sP

  forward_function snapNumToRedshift
  
  ; adiabatic index default (valid for ComparisonProject)
  if not keyword_set(gamma)     then gamma = 5.0/3.0
  
  atime = snapNumToRedshift(sP=sP,/time)
  a3inv = 1.0 / (atime*atime*atime)
  stop ; TODO check this (dens=dens*a3inv but in the tracers only the dens^gamma not in the pressure)
  
  pressure = (gamma-1.0) * u * dens
  entropy  = pressure / (dens^gamma)
  
  return, entropy
end

; rhoRatioToCrit(): normalize density by the critical -baryon- density at some redshift

function rhoRatioToCrit, rho, omega_b=omega_b, redshift=redshift

  units = getUnits()
  
  ; cosmo
  omega_m   = 0.27
  omega_L   = 0.73
  omega_k   = 0.0
  
  ; default omega_b (valid for ComparisonProject)
  if not keyword_set(omega_b) then omega_b = 0.044
  
  rho_b = omega_b * units.rhoCrit
  
  ; scale for redshift other than zero
  if keyword_set(redshift) then begin
    H_z_fact = ( omega_m*(1+redshift)^3.0 + omega_L + omega_k*(1+redshift)^2.0 )
    rho_b *= H_z_fact
  endif
  
  return, rho/rho_b

end

; critBaryonRatioToCode(): convert a ratio of the critical baryon density at some redshift to a code density

function critBaryonRatioToCode, ratioToCritB, redshift=redshift

  units = getUnits()

  ; cosmo
  omega_m   = 0.27
  omega_L   = 0.73
  omega_k   = 0.0
  omega_b = 0.044
  
  codeRho = ratioToCritB * omega_b * units.rhoCrit
  
  ; scale for redshift other than zero
  if keyword_set(redshift) then begin
    H_z_fact = ( omega_m*(1+redshift)^3.0 + omega_L + omega_k*(1+redshift)^2.0 )
    codeRho *= H_z_fact
  endif
  
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

  ; config
  H0 = 70.0 ;km/s/Mpc
  Omega_m = 0.27
  
  age = 2*asinh(sqrt( (1-Omega_m)/Omega_m ) * (1+z)^(-3.0/2.0)) / $
             (H0 * 3 * sqrt(1-Omega_m))
  
  return, age * 3.085678e+19 / 3.15567e+7 / 1e9 ;Gyr

end

