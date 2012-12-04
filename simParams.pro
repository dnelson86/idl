; simParams.pro
; return structure of simulation and analysis parameters with consistent information
; dnelson nov.2012

function simParams, res=res, run=run, redshift=redshift, snap=snap, f=f

  forward_function redshiftToSnapNum

  if not keyword_set(res) or not keyword_set(run) then begin
     print,'Error: simParams: arguments not specified.'
     exit
  endif

  r = {simPath:      '',$    ; root path to simulation snapshots and group catalogs
       arepoPath:    '',$    ; root path to Arepo and param.txt for e.g. projections/fof
       savPrefix:    '',$    ; save prefix for simulation (make unique, e.g. 'G')
       saveTag:      '',$    ; save string: trVel, trMC, or SPH
       plotPrefix:   '',$    ; plot prefix for simulation (make unique, e.g. 'GR')
       boxSize:      0.0,$   ; boxsize of simulation, kpc
       targetGasMass:0.0,$   ; refinement/derefinement target, equal to SPH gas mass in equivalent run
       trMassConst:  0.0,$   ; mass per tracerMC under equal mass assumption (=TargetGasMass/trMCPerCell)
       trMCPerCell:  0,$     ; starting number of monte carlo tracers per cell (copied from f input, 0=none)
       gravSoft:     0.0,$   ; gravitational softening length (ckpc)
       snapRange:    [0,0],$ ; snapshot range of simulation
       groupCatRange:[0,0],$ ; snapshot range of fof/subfind catalogs (subset of above)
       plotPath:     '',$    ; working path to put plots
       derivPath:    '',$    ; path to put derivative (data) files
       snap:         -1,$    ; convenience for passing between functions
       res:          0,$     ; copied from input
       run:          '',$    ; copied from input
       redshift:     -1.0,$  ; copied from input
       $ ; analysis parameters:
       minNumGasPart: 0,$    ; minimum number of gas particles required to use subgroup (unused)
       radcut_rvir:   0.15,$ ; galcat: fraction of rvir as maximum for gal/stars, minimum for gmem (zero to disable)
       radcut_out:    1.5,$  ; galcat: fraction of rvir as maximum for gmem
       galcut_T:      6.0,$  ; galcat: temp coefficient for (rho,temp) galaxy cut
       galcut_rho:    0.25,$ ; galcat: dens coefficient for (rho,temp) galaxy cut
       rVirFacs:      [1.0,0.75,0.5,0.25,0.15,0.05,0.01],$ ; search for accretion times across these fractions of the virial radius
       TcutVals:      [5.3,5.5,5.7],$ ; log(K) for constant threshold comparisons
       TvirVals:      [1.0,0.8,0.4],$ ; T/Tvir coefficients for variable threshold comparisons
       $ ; plotting/vis parameters:
       colorsA:       [getColor24(['00'x,'eb'x,'00'x]),getColor24(['00'x,'bd'x,'00'x]),getColor24(['00'x,'90'x,'00'x])],$ ; green 128,256,512
       colorsG:       [getColor24(['e6'x,'7a'x,'22'x]),getColor24(['b3'x,'5f'x,'1b'x]),getColor24(['80'x,'44'x,'13'x])],$ ; brown 128,256,512
       $ ;colorsA:       [getColor24([255,200,200]),getColor24([255,100,100]),getColor24([255,0,0])],$ ; red 128,256,512 (alternative)
       $ ;colorsG:       [getColor24([200,200,255]),getColor24([100,100,255]),getColor24([0,0,255])],$ ; blue 128,256,512 (alternative)
       radIndHaloAcc: 0,$     ; 1.0 rvir crossing for halo accretion
       radIndGalAcc:  4,$     ; 0.15 rvir crossing for galaxy accretion (or entering rho,temp definition)
       gfmElements:   ['H','He','C','N','O','Ne','Mg','Si','Fe'] ,$
       gfmNumElements: 0, $ ; set to >=1 for GFM runs outputting abundances by metal
       gfmWinds: 0 $ ; set to 1 for GFM_WINDS
      }

  ; copy inputs
  if (isnumeric(res)) then $
    r.res = res
  r.run = run
  if n_elements(redshift) gt 0 then begin
    r.redshift = redshift
    if n_elements(snap) gt 0 then print, 'Warning: snap and redshift both specified!'
  endif
  if n_elements(snap) gt 0 then $
    r.snap = snap
 
  ; sims.feedback (Velocity + f=10 Monte Carlo) 128,256 @ 20Mpc w/ "fiducial winds+BH feedback"
  if (run eq 'feedback') then begin
    r.minNumGasPart  = -1 ; no additional cut
    r.trMCPerCell    = 10
    r.gfmNumElements = 9
    r.gfmWinds       = 1
    
    if res eq 128 or res eq 256 then begin 
      r.boxSize       = 20000.0
      r.snapRange     = [0,139]
      r.groupCatRange = [45,139] ; z6=46, z5=50, z4=55, z3=61, z2=69
    endif else begin
      message,'res error'
    endelse
    
    if res eq 128 then r.targetGasMass = 4.76446157e-03
    if res eq 256 then r.targetGasMass = 5.95556796e-04
    
    r.trMassConst = r.targetGasMass / r.trMCPerCell
    
    if res eq 128 then r.gravSoft = 4.0
    if res eq 256 then r.gravSoft = 2.0
  
    r.simPath    = '/n/home07/dnelson/sims.feedback/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/output/'
    r.arepoPath  = '/n/home07/dnelson/sims.feedback/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/'
    r.savPrefix  = 'F'
    r.saveTag    = 'feMC'
    r.plotPrefix = 'feMC'
    r.plotPath   = '/n/home07/dnelson/coldflows/'
    r.derivPath  = '/n/home07/dnelson/sims.feedback/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/data.files/'
    
    ; if f=-1 use velocity tracers
    if keyword_set(f) then begin
      if f ne -1 then message,'Error.' ; only valid input is -1
      r.trMCPerCell = -1
      r.plotPrefix = 'feVel'
      r.saveTag    = 'feVel'
    endif
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
 
  ; ComparisonProject GADGET 128,256,512 @ 20Mpc, 320,640 @ 40Mpc
  if (run eq 'gadget') or (run eq 'gadget_rad') then begin
    r.minNumGasPart = -1 ; no additional cut
    r.trMCPerCell   = 0  ; none (SPH)
    
    if keyword_set(f) then stop ; shouldn't be specified
    if res ne 128 and res ne 256 and res ne 512 then stop ; only resolutions that exist now
  
    if res eq 128 or res eq 256 or res eq 512 then begin 
      r.boxSize       = 20000.0
      r.snapRange     = [0,313]
      r.groupCatRange = [50,313]
    endif
    if res eq 320 or res eq 640 then begin
      r.boxSize       = 40000.0
      r.snapRange     = [0,0]
      r.groupCatRange = [0,0]
    endif
  
    if res eq 128 then r.targetGasMass = 4.76446157e-03
    if res eq 256 then r.targetGasMass = 5.95556796e-04
    if res eq 512 then r.targetGasMass = 7.44447120e-05
    if res eq 320 then r.targetGasMass = 2.43940430e-03
    if res eq 640 then r.targetGasMass = 3.04925537e-04
    
    if res eq 128 then r.gravSoft = 4.0
    if res eq 256 then r.gravSoft = 2.0
    if res eq 512 then r.gravSoft = 1.0
    if res eq 320 then r.gravSoft = 3.2
    if res eq 640 then r.gravSoft = 1.6
    
    r.simPath    = '/n/home07/dnelson/sims.gadget/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/output/'
    r.arepoPath  = '/n/home07/dnelson/sims.gadget/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/'
    r.savPrefix  = 'G'
    r.saveTag    = 'SPH'
    r.plotPrefix = 'G'
    r.plotPath   = '/n/home07/dnelson/coldflows/'
    r.derivPath  = '/n/home07/dnelson/sims.gadget/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/data.files/'
    
    ; if radially restricted group catalogs, modify derivPath
    if run eq 'gadget_rad' then begin
      r.plotPrefix = r.plotPrefix + 'rad'
      r.derivPath = '/n/home07/dnelson/sims.gadget/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/data.files.rad/'
    endif
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  ; TEMP GADGET 512 old groups and data.files (corresponding to the halo comparison project / for vis)
  if (run eq 'gadgetold') then begin
    r.minNumGasPart = -1 ; no additional cut
    r.trMCPerCell   = 0  ; none (SPH)
    
    if keyword_set(f) then stop ; shouldn't be specified
    if res ne 512 then stop ; only 512
  
    r.boxSize       = 20000.0
    r.snapRange     = [0,313]
    r.groupCatRange = [50,313]
    r.targetGasMass = 7.44447120e-05
    r.gravSoft      = 1.0
    
    r.simPath    = '/n/home07/dnelson/sims.gadget/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_old/output/'
    r.arepoPath  = '/n/home07/dnelson/sims.gadget/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_old/'
    r.savPrefix  = 'G'
    r.saveTag    = 'SPH'
    r.plotPrefix = 'Gold'
    r.plotPath   = '/n/home07/dnelson/coldflows/'
    r.derivPath  = '/n/home07/dnelson/sims.gadget/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_old/data.files/'
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  ; sims.tracers (Velocity + f=10 Monte Carlo) 128,256,512 @ 20Mpc, 320,640 @ 40Mpc
  if (run eq 'tracer') or (run eq 'tracer_rad') or (run eq 'tracer_nouv') then begin
    r.minNumGasPart = -1 ; no additional cut
    r.trMCPerCell   = 10
    
    if res eq 320 or res eq 640 then stop ; don't exist yet
    
    if res eq 128 or res eq 256 or res eq 512 then begin 
      r.boxSize       = 20000.0
      r.snapRange     = [0,313]
      r.groupCatRange = [50,313]
    endif
    if res eq 320 or res eq 640 then begin
      r.boxSize       = 40000.0
      r.snapRange     = [0,0]
      r.groupCatRange = [0,0]
    endif
    
    if res eq 128 then r.targetGasMass = 4.76446157e-03
    if res eq 256 then r.targetGasMass = 5.95556796e-04
    if res eq 512 then r.targetGasMass = 7.44447120e-05
    if res eq 320 then r.targetGasMass = 2.43940430e-03
    if res eq 640 then r.targetGasMass = 3.04925537e-04
    
    r.trMassConst = r.targetGasMass / r.trMCPerCell
    
    if res eq 128 then r.gravSoft = 4.0
    if res eq 256 then r.gravSoft = 2.0
    if res eq 512 then r.gravSoft = 1.0
    if res eq 320 then r.gravSoft = 3.2
    if res eq 640 then r.gravSoft = 1.6
  
    r.simPath    = '/n/home07/dnelson/sims.tracers/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/output/'
    r.arepoPath  = '/n/home07/dnelson/sims.tracers/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/'
    r.savPrefix  = 'N'
    r.saveTag    = 'trMC'
    r.plotPrefix = 'trMC'
    r.plotPath   = '/n/home07/dnelson/coldflows/'
    r.derivPath  = '/n/home07/dnelson/sims.tracers/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/data.files/'
    
    ; if f=-1 use velocity tracers
    if keyword_set(f) then begin
      if f ne -1 then message,'Error.' ; only valid input is -1
      r.trMCPerCell = -1
      r.plotPrefix = 'trVel'
      r.saveTag    = 'trVel'
    endif
    
    ; if radially restricted group catalogs, modify derivPath
    if run eq 'tracer_rad' then begin
      r.plotPrefix = r.plotPrefix + 'rad'
      r.derivPath = '/n/home07/dnelson/sims.tracers/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/data.files.rad/'
    endif
    
    ; if noUV run, modify details and paths
    if run eq 'tracer_nouv' then begin
      if res ne 256 then stop ; only 256 exists
      
      r.snapRange     = [0,139]
      r.groupCatRange = [45,139]
      
      r.simPath    = '/n/home07/dnelson/sims.tracers/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_noUV/output/'
      r.arepoPath  = '/n/home07/dnelson/sims.tracers/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_noUV/'
      r.savPrefix  = 'N'
      r.saveTag    = 'trMC'
      r.plotPrefix = 'trNoUV'
      r.plotPath   = '/n/home07/dnelson/coldflows/'
      r.derivPath  = '/n/home07/dnelson/sims.tracers/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_noUV/data.files/'
    endif
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  ; ComparisonProject arepo (no tracers) 128,256,512 @ 20Mpc
  if (run eq 'arepo') then begin
    r.minNumGasPart = -1 ; no additional cut
    r.trMCPerCell   = 0 ; no actual tracers ;, but flag as non-sph
    
    if keyword_set(f) then stop ; shouldn't be specified
    if res ne 128 and res ne 256 and res ne 512 then stop
    
    r.boxSize       = 20000.0
    r.snapRange     = [0,313]
    r.groupCatRange = [50,313]
    
    if res eq 128 then r.targetGasMass = 4.76446157e-03
    if res eq 256 then r.targetGasMass = 5.95556796e-04
    if res eq 512 then r.targetGasMass = 7.44447120e-05
    
    r.trMassConst = 0.0
    
    if res eq 128 then r.gravSoft = 4.0
    if res eq 256 then r.gravSoft = 2.0
    if res eq 512 then r.gravSoft = 1.0
  
    r.simPath    = '/n/home07/dnelson/sims.tracers/A'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/output/'
    r.arepoPath  = '/n/home07/dnelson/sims.tracers/A'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/'
    r.savPrefix  = 'A'
    r.plotPrefix = 'A'
    r.plotPath   = '/n/home07/dnelson/coldflows/'
    r.derivPath  = '/n/home07/dnelson/sims.tracers/A'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/data.files/'
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  if (run eq 'dev.tracer') then begin
    r.boxSize       = 20000.0
    r.snapRange     = [0,76]
    r.groupCatRange = [25,76]
    r.minNumGasPart = -1 ; no additional cut
    r.trMCPerCell   = -1 ; velocity tracers

    if keyword_set(f) then stop
    if (res ne 128) then stop ; only 128 exists
    
    r.targetGasMass = 4.76446157e-03
    r.trMassConst   = r.targetGasMass
    
    r.gravSoft      = 4.0

    r.simPath    = '/n/scratch2/hernquist_lab/dnelson/dev.tracer/cosmobox.'+str(res)+'_20Mpc/output/'
    r.arepoPath  = '/n/home07/dnelson/dev.tracer/cosmobox.'+str(res)+'_20Mpc/'
    r.savPrefix  = 'D'
    r.plotPrefix = 'D'
    r.plotPath   = '/n/home07/dnelson/dev.tracer/'
    r.derivPath  = '/n/home07/dnelson/dev.tracer/cosmobox.'+str(res)+'_20Mpc/data.files/'
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  if (run eq 'dev.tracer.nonrad') then begin
    r.boxSize       = 20000.0
    r.snapRange     = [0,38]
    r.groupCatRange = [15,38]
    r.minNumGasPart = -1 ; no additional cut
    r.trMCPerCell   = -1 ; velocity tracers
    
    if keyword_set(f) then stop ; shouldn't be set
    if (res ne 128) and (res ne 256) then stop ; only 128 and 256 exist
  
    if (res eq 128) then r.targetGasMass = 4.76446157e-03
    if (res eq 256) then r.targetGasMass = 5.95556796e-04
    
    r.trMassConst   = r.targetGasMass
    
    if (res eq 128) then r.gravSoft = 4.0
    if (res eq 256) then r.gravSoft = 2.0
  
    r.simPath    = '/n/home07/dnelson/dev.tracer/cosmobox.'+str(res)+'_20Mpc_nonrad/output/'
    r.arepoPath  = '/n/home07/dnelson/dev.tracer/cosmobox.'+str(res)+'_20Mpc_nonrad/'
    r.savPrefix  = 'N'
    r.plotPrefix = 'N'
    r.plotPath   = '/n/home07/dnelson/dev.tracer/'
    r.derivPath  = '/n/home07/dnelson/dev.tracer/cosmobox.'+str(res)+'_20Mpc_nonrad/data.files/'
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  if (run eq 'dev.tracerMC.nonrad') then begin ; DEV.ref
    r.boxSize       = 20000.0
    r.snapRange     = [0,38]
    r.groupCatRange = [15,38]
    
    if not keyword_set(f) then stop
    if (res ne 128) and (res ne 256) then stop ; only 128 and 256 exist
    r.trMCPerCell = fix(f)
  
    if (res eq 128) then r.targetGasMass = 4.76446157e-03
    if (res eq 256) then r.targetGasMass = 5.95556796e-04
    
    r.trMassConst = r.targetGasMass / float(r.trMCPerCell)
    
    if (res eq 128) then r.gravSoft = 4.0
    if (res eq 256) then r.gravSoft = 2.0
  
    r.simPath    = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.nonrad.ref/output/'
    r.arepoPath  = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.nonrad.ref/'
    r.savPrefix  = 'N.f'+f
    r.plotPrefix = 'N.f'+f
    r.plotPath   = '/n/home07/dnelson/dev.tracerMC/'
    r.derivPath  = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.nonrad.ref/data.files/'
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  if (run eq 'dev.tracerMC') then begin ; DEV coolSF.GFM.ref
    r.boxSize       = 20000.0
    r.snapRange     = [0,38]
    r.groupCatRange = [15,38]
    
    if not keyword_set(f) then stop
    if (res ne 128) and (res ne 256) then stop ; only 128 and 256 exist
    r.trMCPerCell = fix(f)
  
    if (res eq 128) then r.targetGasMass = 4.76446157e-03
    if (res eq 256) then r.targetGasMass = 5.95556796e-04
    
    r.trMassConst = r.targetGasMass / float(r.trMCPerCell)
    
    if (res eq 128) then r.gravSoft = 4.0
    if (res eq 256) then r.gravSoft = 2.0
  
    r.simPath    = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.coolSF.GFM.ref/output/'
    r.arepoPath  = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.coolSF.GFM.ref/'
    r.savPrefix  = 'T.f'+f
    r.plotPrefix = 'T.f'+f
    r.plotPath   = '/n/home07/dnelson/dev.tracerMC/'
    r.derivPath  = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.coolSF.GFM.ref/data.files/'
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  print,'simParams: ERROR.'
  stop
end
