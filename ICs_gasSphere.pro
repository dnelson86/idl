; ICs_gasSphere.pro
; initial condition generation - gas cloud in static hernquist/NFW potential
; dnelson jan.2012

@helper

; setupTwoSphereCollision(): copy gas and tracers from input filename twice, and offset from center

pro setupTwoSphereCollision, gasOnly=gasOnly

  input_fname  = "gasSphere.gasonly.1e4.nfw.M1e2.norot.dat.hdf5"
  output_fname = "col2Sph.gasonly.1e4.nfw.M1e2.norot.dat"

  offsetVec = [200.0,0,0] ;kpc in each direction
  boxSize   = 5000.0 ;kpc
  nBackGrid = 16 ;^3, nested twice
  
  ; load input sphere
  h = loadSnapshotHeader(input_fname,snapNum='none',/verbose)

  gas_pos = loadSnapshotSubset(input_fname,snapNum='none',partType='gas',field='pos')
  gas_vel = loadSnapshotSubset(input_fname,snapNum='none',partType='gas',field='vel')
  gas_u   = loadSnapshotSubset(input_fname,snapNum='none',partType='gas',field='u')

  ; offset template position to box center
  gas_pos += boxSize/2.0
  
  ; decrease thermal energy of gas by a factor to perturb the individual spheres out of equilibrium
  ; and towards collapse (also help prevent slow expansion during impact)
  gas_u *= 0.5

  ; arrays
  nGasFinal = n_elements(gas_ids)*2
  
  ; duplicate, offset positions (right)
  gas_pos_f = [gas_pos[0,*]+offsetVec[0],gas_pos[1,*]+offsetVec[1],gas_pos[2,*]+offsetVec[2]]
  gas_vel_f = [gas_vel]
  gas_u_f   = [gas_u]
  
  ; duplicate, offset positions (left)
  gas_pos_f = [[gas_pos_f],[gas_pos[0,*]-offsetVec[0],gas_pos[1,*]-offsetVec[1],gas_pos[2,*]-offsetVec[2]]]
  gas_vel_f = [[gas_vel_f],[gas_vel]]
  gas_u_f   = [gas_u_f,gas_u]
  
  ; generate matching ids and tracer masses
  gas_id = lindgen(n_elements(gas_u_f)) ;l64indgen()
  
  gas_mass = fltarr(n_elements(gas_id)) + float(h.masstable[0]) ;copy same const gas mass

  ; form gas (partType=0) and tracers (partType=3)
  gas    = {pos:gas_pos_f,vel:gas_vel_f,id:gas_id,u:gas_u_f,mass:gas_mass}

  ; add nested background grid
  boxSizeInner = ceil((max(gas.pos)-min(gas.pos))/100) * 100.0

  gas = addICBackgroundGrid(gas, boxSize=boxSizeInner, boxCen=boxSize/2.0, nBackGrid=nBackGrid)
  gas = addICBackgroundGrid(gas, boxSize=boxSize, nBackGrid=nBackGrid)
  
  if not keyword_set(gasOnly) then begin
    print,'Including tracers!'
    ; repeat all the same with the tracers
    tr_pos  = loadSnapshotSubset(input_fname,snapNum='none',partType='tracer',field='pos')
    tr_vel  = loadSnapshotSubset(input_fname,snapNum='none',partType='tracer',field='vel')
    tr_pos  += boxSize/2.0
    nTrFinal  = n_elements(tr_ids)*2
    
    tr_pos_f = [tr_pos[0,*]+offsetVec[0],tr_pos[1,*]+offsetVec[1],tr_pos[2,*]+offsetVec[2]]
    tr_vel_f = [tr_vel]
    
    tr_pos_f = [[tr_pos_f],[tr_pos[0,*]-offsetVec[0],tr_pos[1,*]-offsetVec[1],tr_pos[2,*]-offsetVec[2]]]
    tr_vel_f = [[tr_vel_f],[tr_vel]]
    
    tr_id  = lindgen(n_elements(tr_pos_f[0,*])) + 100000000L ;l64indgen()
    
    tr_mass  = fltarr(n_elements(tr_id)) 
    
    tracer = {pos:tr_pos_f,vel:tr_vel_f,id:tr_id,mass:tr_mass}  

    writeICFile,output_fname,part0=gas,part3=tracer;,/longIDs
  endif else begin
    print,'Writing gas only!'
    ; save
    writeICFile,output_fname,part0=gas;,/longIDs
  endelse
  
  print,'Suggested target gas mass: '+string(h.masstable[0])

end

; setupSingleSphereIso(): isolated gas sphere setup (add background and center)

pro setupSingleSphereIso

  input_fname  = "gasSphere.gasonly.1e4.nfw.M1e2.norot.dat.hdf5"
  output_fname = "gasSphere.gasonly.1e4.nfw.M1e2.norot.noBack.dat"

  boxSize   = 5000.0 ;kpc
  nBackGrid = 16 ;^3, nested twice
  
  ; load input sphere
  h = loadSnapshotHeader(fileName=input_fname,/verbose)

  gas_pos  = loadSnapshotSubset(fileName=input_fname,partType='gas',field='pos')
  gas_vel  = loadSnapshotSubset(fileName=input_fname,partType='gas',field='vel')
  gas_id   = loadSnapshotSubset(fileName=input_fname,partType='gas',field='ids')
  gas_u    = loadSnapshotSubset(fileName=input_fname,partType='gas',field='u')
  
  ; offset template position to box center
  gas_pos += boxSize/2.0
  
  ; form gas (partType=0) and tracers (partType=3)
  gas_mass = fltarr(n_elements(gas_id)) + float(h.masstable[0]) ;copy same const gas mass
  
  gas = {pos:gas_pos,vel:gas_vel,id:gas_id,u:gas_u,mass:gas_mass}

  ; add nested background grid
  boxSizeInner = ceil((max(gas.pos)-min(gas.pos))/100) * 100.0
  uthermBackGrid = 0 ;1.0e4
  
  ;gas = addICBackgroundGrid(gas, boxSize=boxSizeInner,uthermBackGrid=uthermBackGrid, $
  ;                          boxCen=[boxSize/2.0,boxSize/2.0,boxSize/2.0], nBackGrid=nBackGrid)
  ;gas = addICBackgroundGrid(gas, boxSize=boxSize, nBackGrid=nBackGrid, uthermBackGrid=uthermBackGrid)

  ; save
  writeICFile,output_fname,part0=gas
  print,'Suggested target gas mass: '+string(h.masstable[0])
  print,minmax(gas.pos)

end

; setupFilamentTest(): single cold filament accreting onto an isolated gas sphere

pro setupFilamentTest

  ; cylinder config
  r200 = 162.6 ; kpc nfw
  cyl_length = 100.0  ;kpc
  cyl_radius = 2.0   ;kpc
  cyl_rstart = 1.0*r200 ;kpc (r200 nfw)
  cyl_init_vel = 100.0 ;km/s towards halo
  
  cyl_densratio = 1000 ; specify density as ratio to critical baryon density at z=3
  cyl_dens = critBaryonRatioToCode(cyl_densratio,redshift=3.0)

  cyl_r_res = 10
  cyl_l_res = cyl_r_res * (cyl_length/cyl_radius)
  
  seed = 424242L

  ; gasSphere config
  input_fname  = "gasSphere.gasonly.1e4.nfw.M1e2.norot.dat.hdf5"
  output_fname = "cylTest.1e4.M1e2.norot.c"+str(cyl_r_res)+".d"+str(cyl_densratio)+$
    ".r"+str(fix(cyl_radius))+".withBack.dat"

  boxSize   = 5000.0 ;kpc
  nBackGrid = 16 ;^3, nested twice
  
  ; load input sphere
  h = loadSnapshotHeader(fileName=input_fname,/verbose)

  gas_pos  = loadSnapshotSubset(fileName=input_fname,partType='gas',field='pos')
  gas_vel  = loadSnapshotSubset(fileName=input_fname,partType='gas',field='vel')
  gas_id   = loadSnapshotSubset(fileName=input_fname,partType='gas',field='ids')
  gas_u    = loadSnapshotSubset(fileName=input_fname,partType='gas',field='u')
  
  ; offset template position to box center
  gas_pos += boxSize/2.0
  
  ; form gas for sphere
  gas_mass = fltarr(n_elements(gas_id)) + float(h.masstable[0]) ;copy same const gas mass

  ; cylinder properties
  halo_vol = 4.0/3.0*!pi*r200^3.0
  halo_meandens = total(gas_mass) / halo_vol ;1e10 msun/kpc^3
  
  cyl_vol   = !pi*cyl_radius^2.0*cyl_length
  cyl_npart = cyl_r_res * cyl_l_res
  cyl_masspart = cyl_vol * cyl_dens / cyl_npart

  ; form gas for cylinder (along z-axis)
  rnd_angle = randomu(seed,cyl_npart,/double) * 360.0
  rnd_radii = randomu(seed,cyl_npart,/double) * cyl_radius
  
  cyl_x = rnd_radii * sin(rnd_angle*!dtor)
  cyl_y = rnd_radii * cos(rnd_angle*!dtor)
  cyl_z = randomu(seed,cyl_npart,/double) * cyl_length
  
  cyl_z += cyl_rstart ; displace cylinder to starting radius
  
  cyl_pos = fltarr(3,cyl_npart)
  cyl_pos[0,*] = cyl_x
  cyl_pos[1,*] = cyl_y
  cyl_pos[2,*] = cyl_z
  cyl_pos += boxSize/2.0 ; offset to box center  
  
  ; other cylinder properties
  cyl_vel = fltarr(3,cyl_npart)
  cyl_vel[2,*] = -cyl_init_vel
  cyl_id  = lindgen(cyl_npart) + max(gas_id) + 1L
  cyl_u   = fltarr(cyl_npart) + 200.0 ;~10000K similar to the outskirts of the halo
  cyl_mass = fltarr(cyl_npart) + cyl_masspart

  ; create combined gas struct
  gas = {pos:[[gas_pos],[cyl_pos]],vel:[[gas_vel],[cyl_vel]],id:[gas_id,cyl_id],$
         u:[gas_u,cyl_u],mass:[gas_mass,cyl_mass]}

  ; add nested background grid
  boxSizeInner = ceil((max(gas.pos)-min(gas.pos))/200) * 200.0
  uthermBackGrid = 0 ;1.0e4

  gas = addICBackgroundGrid(gas, boxSize=boxSizeInner, uthermBackGrid=uthermBackGrid, $
                            boxCen=[boxSize/2.0,boxSize/2.0,boxSize/2.0], nBackGrid=nBackGrid)
  gas = addICBackgroundGrid(gas, boxSize=boxSize, nBackGrid=nBackGrid, uthermBackGrid=10*uthermBackGrid)

  ; save
  writeICFile,output_fname,part0=gas
  print,'Cylinder particle count: '+string(cyl_npart)
  print,'t=0 Cylinder id range: '+string(min(cyl_id))+" - "+string(max(cyl_id))
  print,'Suggested target gas mass: '+string(h.masstable[0])

end

; generate hot gas sphere in hydrostatic equilibrium with a Hernquist halo profile
; based on Mark Vogelsberger's python version

pro gen_ICs_gasSphere
  forward_function GasRho,HaloRho,Rho,GasMass,HaloMass,Mass,Sigma_Integrand,Sigma
  COMMON gv,gas_R0,HQ_a,HQ_M,gas_frac,INTERPOL_R_MAX_GAS,G ;globals

  message,'doublecheck before using'

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
