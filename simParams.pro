; simParams.pro
; return structure of simulation parameters with consistent information
; dnelson jun.2012

function simParams, res=res, run=run, redshift=redshift, snap=snap, f=f

  forward_function redshiftToSnapNum

  if not keyword_set(res) or not keyword_set(run) then begin
     print,'Error: simParams: arguments not specified.'
     exit
  endif

  r = {simPath:      '',$    ; root path to simulation snapshots and group catalogs
       arepoPath:    '',$    ; root path to Arepo and param.txt for e.g. projections/fof
       savPrefix:    '',$    ; save prefix for simulation (make unique, e.g. 'G')
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
       minNumGasPart: 0}     ; minimum number of gas particles required to use subgroup
       
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
    r.plotPrefix = 'Gold'
    r.plotPath   = '/n/home07/dnelson/coldflows/'
    r.derivPath  = '/n/home07/dnelson/sims.gadget/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_old/data.files/'
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  ; sims.tracers (Velocity + f=10 Monte Carlo) 128,256,512 @ 20Mpc, 320,640 @ 40Mpc
  if (run eq 'tracer') or (run eq 'tracer_rad') then begin
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
    r.plotPrefix = 'trMC'
    r.plotPath   = '/n/home07/dnelson/coldflows/'
    r.derivPath  = '/n/home07/dnelson/sims.tracers/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/data.files/'
    
    ; if f=-1 use velocity tracers
    if keyword_set(f) then begin
      if f ne -1 then message,'Error.' ; only valid input is -1
      r.trMCPerCell = -1
      r.plotPrefix = 'trVel'
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
      r.groupCatRange = [45,139] ;z6=45, z5=49, z4=54, z3=60, z2=68
      
      r.simPath    = '/n/home07/dnelson/sims.tracers/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_noUV/output/'
      r.arepoPath  = '/n/home07/dnelson/sims.tracers/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_noUV/'
      r.savPrefix  = 'N'
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
    r.trMCPerCell   = 1 ; no actual tracers, but flag as non-sph
    
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
  
  if (run eq 'dev.tracerMC.SPT') then begin ; shy's bugtest
    r.boxSize       = 20000.0
    r.snapRange     = [1,2]
    r.groupCatRange = [1,2]
    
    if keyword_set(f) then message,'do not set f'
    f = '10'

    if res ne 128 then message,'only res 128 exists'
    r.trMCPerCell = fix(f)
  
    r.targetGasMass = 4.76446157e-03
    r.trMassConst = r.targetGasMass / float(r.trMCPerCell)
    r.gravSoft = 4.0
  
    r.simPath    = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.SPT/output/'
    r.arepoPath  = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.SPT/'
    r.savPrefix  = 'S.f'+f
    r.plotPrefix = 'S.f'+f
    r.plotPath   = '/n/home07/dnelson/dev.tracerMC/'
    r.derivPath  = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.SPT/data.files/'
    
    if (n_elements(redshift) eq 1) then message,'error'
    
    return,r
  endif
  
  print,'simParams: ERROR.'
  stop
end
