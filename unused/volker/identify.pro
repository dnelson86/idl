pro identify

Num      = 4
SimDir   = "/u/dnelson/sims.TNG_method/L25n128_0000/" ; path to simulation directory
;SimDir   = "/hits/universe/IllustrisTNG/L35n540TNG/" ; path to simulation directory


Base     = SimDir + "output/"
Snapbase = "snap"               ; base name of snapshot files



; System of units
UnitLength_in_cm         =    3.085678d21 ;  Mpc/h
UnitMass_in_g            =    1.989d43    ;  1.0e10 Msun/h
UnitVelocity_in_cm_per_s =    1.0d5       ;  1 km/sec

UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s 
UnitDensity_in_cgs= UnitMass_in_g/ UnitLength_in_cm^3
UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ UnitTime_in_s^2
UnitEnergy_in_cgs= UnitMass_in_g * UnitLength_in_cm^2 / UnitTime_in_s^2
GRAVITY    = 6.672d-8
BOLTZMANN  = 1.3806d-16
PROTONMASS = 1.6726d-24
HUBBLE     = 3.2407789d-18   ; in h/sec
 
G=GRAVITY/ UnitLength_in_cm^3 * UnitMass_in_g * UnitTime_in_s^2
H0 = HUBBLE * UnitTime_in_s


Xh=0.76D                        ; mass fraction of hydrogen
HubbleParam  = 0.6774
Omegab       = 0.0486
gamma_gas    = 5.0/3

RhoMean = Omegab * 3 * H0^2 / (8*!PI * G)



get_siminfo, Base, SnapBase, Num, NumGroupFiles, NumSnapFiles, Ngroups, Nsubgroups, NumPartType, Ti, BoxSize
get_group_len_and_offset, Base, Num, NumGroupFiles, Ngroups, GroupLenType, GroupOffsetType

Rho      = get_snap_field(0, 'Density', Base, SnapBase, Num, NumSnapFiles, NumPartType)
U        = get_snap_field(0, 'InternalEnergy', Base, SnapBase, Num, NumSnapFiles, NumPartType)
Nelec    = get_snap_field(0, 'ElectronAbundance', Base, SnapBase, Num, NumSnapFiles, NumPartType)
Sfr      = get_snap_field(0, 'StarFormationRate', Base, SnapBase, Num, NumSnapFiles, NumPartType)
Vel      = get_snap_field(0, 'Velocities', Base, SnapBase, Num, NumSnapFiles, NumPartType)
Mass     = get_snap_field(0, 'Masses', Base, SnapBase, Num, NumSnapFiles, NumPartType)
Bfield   = get_snap_field(0, 'MagneticField', Base, SnapBase, Num, NumSnapFiles, NumPartType)
DivB     = get_snap_field(0, 'MagneticFieldDivergence', Base, SnapBase, Num, NumSnapFiles, NumPartType)
CoolRate = get_snap_field(0, 'GFM_CoolingRate', Base, SnapBase, Num, NumSnapFiles, NumPartType)

Bfield /= (4*!PI) ; back to heavide-lorentz, now have Bc
Vel *= sqrt(ti)   ; now have peculiar velocity

Bv = Bfield(0,*) * vel(0,*)  +  Bfield(1,*) * vel(1,*)  +  Bfield(2,*) * vel(2,*)

PowellTerm = -1.0/ti * DivB * Bv    ; this is dEpsilon/dt, rightmost source term in Eqn (21) in Pakmor & Springel arxiv:1212.1452


rho_cgs = rho / ti^3 * UnitDensity_in_cgs * HubbleParam * HubbleParam;    /* physical density in cgs units */
Nh_cgs = XH * rho_cgs / PROTONMASS
ratefact = Nh_cgs^2 / rho_cgs
dudt_cgs = ratefact * coolrate

dudt = dudt_cgs * UnitTime_in_s / HubbleParam / UnitVelocity_in_cm_per_s^2  ; now in internal units

Heating = ti^2 * rho * dudt      ; convert heating/cooling rate to DEpsilon/dt quantity

print, 'Vel: ',Vel(*,0),Vel(*,1)
print, 'Powell: ', PowellTerm(0:5)
print, 'Heating: ', Heating(0:5)
print, 'dudt: ', dudt_cgs(0:5)

ratio = PowellTerm/Heating
print, 'Ratio: ', ratio(0:5), mean(ratio)
stop

ind = where(sfr gt 0)
if ind(0) ne -1 then Nelec(ind) = 1.2

MeanWeight= 4.0/(3*Xh+1+4*Xh*Nelec) * PROTONMASS
Temp = MeanWeight/BOLTZMANN * (gamma_gas-1) * U * UnitEnergy_in_cgs/ UnitMass_in_g




; Now let's plot the stuff for a random sample of the particles

sel = randomu(seed, n_elements(Temp))
ind = where((sel lt 0.01))

; pick the sample

rho_sel  = Rho(ind)
Temp_sel = Temp(ind)

PowellTerm_sel = PowellTerm(ind)
Heating_sel    = Heating(ind)



;window, xsize=1200, ysize=1200

;plot, rho_sel/rhomean, Temp_sel, psym=3, /xlog,/ylog, charsize=2.0, ytitle = "Temp [K]", xtitle = "rho/<rho>", yrange=[1.0e2, 1.0e8]


; among the sample, let's select the ones where cooling due to
; Powell source term is larger than radiative heating rate

;ind = where((PowellTerm_sel lt 0) and (Heating_sel gt 0) and ((PowellTerm_sel + Heating_sel) lt 0))

;oplot, rho_sel(ind)/rhomean, Temp_Sel(ind), psym=4, color=255


end

