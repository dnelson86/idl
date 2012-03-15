; simParams.pro
; return structure of simulation parameters with consistent information
; dnelson mar.2012

function simParams, res=res, run=run, redshift=redshift, snap=snap, f=f

  forward_function redshiftToSnapNum

  if not keyword_set(res) or not keyword_set(run) then begin
     print,'Error: simParams: arguments not specified.'
     exit
  endif

  r = {simPath:      '',$    ; root path to simulation snapshots and group catalogs
       arepoPath:    '',$    ; root path to Arepo and param.txt for e.g. projections/fof
       savPrefix:    '',$    ; save prefix for simulation (make unique, e.g. 'G')
       boxSize:      0.0,$   ; boxsize of simulation, kpc
       trMassConst:  0.0,$   ; mass per tracer under equal mass assumption (=TargetGasMass)
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
 
  if (run eq 'gadget') then begin ; ComparisonProject GADGET 128,256,512
    r.boxSize       = 20000.0
    r.snapRange     = [0,314]
    r.groupCatRange = [50,314]
    r.minNumGasPart = -1 ; no additional cut
    r.trMCPerCell   = 0  ; none (SPH)
    
    if keyword_set(f) then stop ; shouldn't be specified
  
    if (res eq 128) then r.gravSoft = 4.0
    if (res eq 256) then r.gravSoft = 2.0
    if (res eq 512) then r.gravSoft = 1.0
  
    r.simPath   = '/n/home07/dnelson/sims.gadget/'+str(res)+'_20Mpc/output/'
    r.arepoPath = '/n/home07/dnelson/sims.gadget/'+str(res)+'_20Mpc/'
    r.savPrefix = 'G'
    r.plotPath  = '/n/home07/dnelson/coldflows/'
    r.derivPath = '/n/home07/dnelson/sims.gadget/'+str(res)+'_20Mpc/data.files/'
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  if (run eq 'tracer') then begin ; sims.tracers (velocity field) 128,256,512
    r.boxSize       = 20000.0
    r.snapRange     = [0,313]
    r.groupCatRange = [50,313]
    r.minNumGasPart = -1 ; no additional cut
    r.trMCPerCell   = -1 ; velocity tracers
    
    if keyword_set(f) then stop ; shouldn't be specified
    if (res eq 128) or (res eq 256) then stop ; DELETED
    
    if (res eq 128) then r.trMassConst = 4.76446157e-03
    if (res eq 256) then r.trMassConst = 5.95556796e-04
    if (res eq 512) then r.trMassConst = 7.44447120e-05
    
    if (res eq 128) then r.gravSoft = 4.0
    if (res eq 256) then r.gravSoft = 2.0
    if (res eq 512) then r.gravSoft = 1.0
  
    r.simPath   = '/n/home07/dnelson/sims.tracers/'+str(res)+'_20Mpc/output/'
    r.arepoPath = '/n/home07/dnelson/sims.tracers/'+str(res)+'_20Mpc/'
    r.savPrefix = 'T'
    r.plotPath   = '/n/home07/dnelson/coldflows/'
    r.derivPath  = '/n/home07/dnelson/sims.tracers/'+str(res)+'_20Mpc/data.files/'
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  if (run eq 'tracerMC') then begin ; sims.tracersMC (Monte Carlo) 128,256,512
    r.boxSize       = 20000.0
    r.snapRange     = [0,313]
    r.groupCatRange = [50,313]
    r.minNumGasPart = -1 ; no additional cut
    r.trMCPerCell   = 20
    
    if keyword_set(f) then stop ; shouldn't be specified, fixed at 20
    if (res eq 512) then stop ; don't exist yet
    
    if (res eq 128) then r.trMassConst = 4.76446157e-03 / r.trMCPerCell
    if (res eq 256) then r.trMassConst = 5.95556796e-04 / r.trMCPerCell
    if (res eq 512) then r.trMassConst = 7.44447120e-05 / r.trMCPerCell
    
    if (res eq 128) then r.gravSoft = 4.0
    if (res eq 256) then r.gravSoft = 2.0
    if (res eq 512) then r.gravSoft = 1.0
  
    r.simPath   = '/n/home07/dnelson/sims.tracersMC/'+str(res)+'_20Mpc/output/'
    r.arepoPath = '/n/home07/dnelson/sims.tracersMC/'+str(res)+'_20Mpc/'
    r.savPrefix = 'M'
    r.plotPath   = '/n/home07/dnelson/coldflows/'
    r.derivPath  = '/n/home07/dnelson/sims.tracersMC/'+str(res)+'_20Mpc/data.files/'
    
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
    
    r.trMassConst   = 4.76446157e-03
    r.gravSoft      = 4.0

    r.simPath    = '/n/scratch2/hernquist_lab/dnelson/dev.tracer/cosmobox.'+str(res)+'_20Mpc/output/'
    r.arepoPath  = '/n/home07/dnelson/dev.tracer/cosmobox.'+str(res)+'_20Mpc/'
    r.savPrefix  = 'D'
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
    
    if keyword_set(f) then stop
  
    if (res eq 128) then r.trMassConst = 4.76446157e-03
    if (res eq 256) then r.trMassConst = 5.95556796e-04
    
    if (res eq 128) then r.gravSoft = 4.0
    if (res eq 256) then r.gravSoft = 2.0
  
    r.simPath    = '/n/home07/dnelson/dev.tracer/cosmobox.'+str(res)+'_20Mpc_nonrad/output/'
    r.arepoPath  = '/n/home07/dnelson/dev.tracer/cosmobox.'+str(res)+'_20Mpc_nonrad/'
    r.savPrefix  = 'N'
    r.plotPath   = '/n/home07/dnelson/dev.tracer/'
    r.derivPath  = '/n/home07/dnelson/dev.tracer/cosmobox.'+str(res)+'_20Mpc_nonrad/data.files/'
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  if (run eq 'tracerMC.nonrad') then begin ; DEV
    r.boxSize       = 20000.0
    r.snapRange     = [0,38]
    r.groupCatRange = [15,38]
    
    if not keyword_set(f) then stop
    r.trMCPerCell = fix(f)
  
    if (res eq 128) then r.trMassConst = 4.76446157e-03 / float(f)
    if (res eq 256) then r.trMassConst = 5.95556796e-04 / float(f)
    
    if (res eq 128) then r.gravSoft = 4.0
    if (res eq 256) then r.gravSoft = 2.0
  
    r.simPath    = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.nonrad.ref/output/'
    r.arepoPath  = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.nonrad.ref/'
    r.savPrefix  = 'N.f'+f
    r.plotPath   = '/n/home07/dnelson/dev.tracerMC/'
    r.derivPath  = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.nonrad.ref/data.files/'
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  if (run eq 'tracerMC.nonrad.noref') then begin ; DEV
    r.boxSize       = 20000.0
    r.snapRange     = [0,38]
    r.groupCatRange = [15,38]
    
    if not keyword_set(f) then stop
    r.trMCPerCell = fix(f)
    
    if (res eq 128) then r.trMassConst = 4.76446157e-03 / float(f)
    if (res eq 256) then r.trMassConst = 5.95556796e-04 / float(f)
    
    if (res eq 128) then r.gravSoft = 4.0
    if (res eq 256) then r.gravSoft = 2.0
  
    r.simPath    = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.nonrad.noref/output/'
    r.arepoPath  = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.nonrad.noref/'
    r.savPrefix  = 'R.f'+f
    r.plotPath   = '/n/home07/dnelson/dev.tracerMC/'
    r.derivPath  = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.nonrad.noref/data.files/'
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  if (run eq 'tracerMC.noref') then begin ; DEV coolSF.GFM
    r.boxSize       = 20000.0
    r.snapRange     = [0,38]
    r.groupCatRange = [15,38]
    
    if not keyword_set(f) then stop
    r.trMCPerCell = fix(f)
  
    if (res eq 128) then r.trMassConst = 4.76446157e-03 / float(f)
    if (res eq 256) then r.trMassConst = 5.95556796e-04 / float(f)
    
    if (res eq 128) then r.gravSoft = 4.0
    if (res eq 256) then r.gravSoft = 2.0
  
    r.simPath    = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.coolSF.GFM.noref/output/'
    r.arepoPath  = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.coolSF.GFM.noref/'
    r.savPrefix  = 'S.f'+f
    r.plotPath   = '/n/home07/dnelson/dev.tracerMC/'
    r.derivPath  = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.coolSF.GFM.noref/data.files/'
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  if (run eq 'tracerMC.ref') then begin ; DEV coolSF.GFM
    r.boxSize       = 20000.0
    r.snapRange     = [0,38]
    r.groupCatRange = [15,38]
    
    if not keyword_set(f) then stop
    r.trMCPerCell = fix(f)
  
    if (res eq 128) then r.trMassConst = 4.76446157e-03 / float(f)
    if (res eq 256) then r.trMassConst = 5.95556796e-04 / float(f)
    
    if (res eq 128) then r.gravSoft = 4.0
    if (res eq 256) then r.gravSoft = 2.0
  
    r.simPath    = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.coolSF.GFM.ref/output/'
    r.arepoPath  = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.coolSF.GFM.ref/'
    r.savPrefix  = 'T.f'+f
    r.plotPath   = '/n/home07/dnelson/dev.tracerMC/'
    r.derivPath  = '/n/home07/dnelson/dev.tracerMC/cosmobox.'+str(res)+'_20Mpc.f'+str(f)+'.coolSF.GFM.ref/data.files/'
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  print,'simParams: ERROR.'
  stop
end