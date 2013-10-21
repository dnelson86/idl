; simParams.pro
; return structure of simulation and analysis parameters with consistent information
; dnelson sep.2013

function fillZoomParams, r, res=res, hInd=hInd

  ; convert res to levelMax if needed
  r.levelMax = res
  if r.levelMax ge 64 then r.levelMax = alog10(r.levelMax)/alog10(2)
  if r.levelMax-round(r.levelMax) gt 1e-6 then message,'Error: Bad res'
  r.levelMax = round(r.levelMax)
  r.res      = r.levelMax

  if r.levelMin lt 7 or r.levelMax gt 11 or r.levelMax lt r.levelMin then message,'res error'
  if r.levelMax eq 8 then message,'Error: Does not exist.'
  
  if r.levelMin ne r.levelMax then $
    if n_elements(hInd) eq 0 then message,'Error: Must specify hInd to load zoom.'
  
  r.zoomLevel = r.levelMax - r.levelMin
  r.hInd      = hInd
  r.zoomShift = [0,0,0] ; for levelMax=7 (unigrid)

  r.targetGasMass = 4.76446157e-03 ; L7
  r.targetGasMass /= (8^r.zoomLevel) ; factor of eight decrease at each increasing zoom level
  
  r.gravSoft = 4.0 ; L7
  r.gravSoft /= (2^r.zoomLevel) ; factor of two decrease at each increasing zoom level
  
  if r.levelMax eq 9  then r.ids_offset =  10000000L
  if r.levelMax eq 10 then r.ids_offset =  50000000L
  if r.levelMax eq 11 then r.ids_offset = 500000000L
  
  if hInd eq 0 then begin
    ;[ 6] hInd =   95 mass = 11.97 rvir = 239.9 vol =  69.3 pos = [ 7469.41  5330.66  3532.18 ]
    ;ref_center = 0.3938, 0.2466, 0.1794
    ;ref_extent = 0.2450, 0.1850, 0.1778
    r.targetHaloInd  = 95
    r.targetHaloPos  = [7469.41, 5330.66, 3532.18]
    r.targetHaloRvir = 239.9 ; ckpc
    r.targetHaloMass = 11.97 ; log msun
    r.targetRedshift = 2.0
    
    if r.levelMax ge 9 then r.zoomShift = [13,32,41]
  endif
  
  if hInd eq 1 then begin
    ;[ 7] hInd =   98 mass = 11.90 rvir = 218.0 vol =  63.7 zL = 2 rVirFac = 3.8 pos = [ 6994.99 16954.28 16613.29 ]
    ;ref_center = 0.3288, 0.8630, 0.8127
    ;ref_extent = 0.2229, 0.1853, 0.1929
    r.targetHaloInd  = 98
    r.targetHaloPos  = [6994.99, 16954.28, 16613.29]
    r.targetHaloRvir = 218.0 ; ckpc
    r.targetHaloMass = 11.90 ; log msun
    r.targetRedshift = 2.0
    
    if r.levelMax ge 9 then r.zoomShift = [21,-46,-40]
  endif
  
  if hInd eq 2 then begin
    ;[ 9] hInd =  104 mass = 11.82 rvir = 214.8 vol =  76.8 zL = 2 rVirFac = 3.8 pos = [ 4260.38  5453.91  6773.12 ]
    ;boxCenter: [  0.2071 ,  0.2695 ,  0.3588 ]
    ;boxExtent: [  0.2001 ,  0.2151 ,  0.2229 ]
    r.targetHaloInd  = 104
    r.targetHaloPos  = [4260.38, 5453.91, 6773.12]
    r.targetHaloRvir = 214.8 ; ckpc
    r.targetHaloMass = 11.82 ; log msun
    r.targetRedshift = 2.0
    
    if r.levelMax ge 9 then r.zoomShift = [37,29,18]
  endif
  
  if hInd eq 3 then begin
    ;[ 0] hInd =   49 mass = 12.17 rvir = 263.3 vol =  88.2 zL = 2 rVirFac = 3.8 pos = [10805.00  8047.92  4638.30 ]
    ;boxCenter: [  0.5658 ,  0.4138 ,  0.2570 ]
    ;boxExtent: [  0.2306 ,  0.2681 ,  0.1782 ]
    r.targetHaloInd  = 49
    r.targetHaloPos  = [10805.00, 8047.92, 4638.30]
    r.targetHaloRvir = 263.3 ; ckpc
    r.targetHaloMass = 12.17 ; log msun
    r.targetRedshift = 2.0
    
    if r.levelMax ge 9 then r.zoomShift = [-8,11,31]
  endif
  
  ; convert zoomShift to zoomShiftPhys
  r.zoomShiftPhys = r.zoomShift / 2.0^r.levelMin * r.boxSize
  
  if r.targetHaloMass eq 0.0 then message,'Error: Unrecognized zoom hInd.'
  if r.zoomLevel eq 0 then message,'Strange'
  
  return, r
end

function simParams, res=res, run=run, redshift=redshift, snap=snap, hInd=hInd, f=f

  forward_function redshiftToSnapNum

  if not keyword_set(res) or not keyword_set(run) then $
     message,'Error: simParams: arguments not specified.'
  
  run = strlowcase(run)

  r = {simPath:      '',$    ; root path to simulation snapshots and group catalogs
       arepoPath:    '',$    ; root path to Arepo and param.txt for e.g. projections/fof
       savPrefix:    '',$    ; save prefix for simulation (make unique, e.g. 'G')
       saveTag:      '',$    ; save string: trVel, trMC, or SPH
	 simName:      '',$    ; label to add to plot legends (e.g. "GADGET", "AREPO", "FEEDBACK")
       plotPrefix:   '',$    ; plot prefix for simulation (make unique, e.g. 'GR')
       snapRange:    [0,0],$ ; snapshot range of simulation
       groupCatRange:[0,0],$ ; snapshot range of fof/subfind catalogs (subset of above)
       plotPath:     '',$    ; working path to put plots
       derivPath:    '',$    ; path to put derivative (data) files
       snap:         -1,$    ; convenience for passing between functions
       run:          '',$    ; copied from input
       redshift:     -1.0,$  ; copied from input
       $
       $ ; run parameters
       res:          0,$     ; copied from input
       boxSize:      0.0,$   ; boxsize of simulation, kpc
       targetGasMass:0.0,$   ; refinement/derefinement target, equal to SPH gas mass in equivalent run
       gravSoft:     0.0,$   ; gravitational softening length (ckpc)
       subboxCen:    [0.0,0.0,0.0]   ,$ ; subbox0 center
       subboxSize:   [0.0,0.0,0.0]   ,$ ; subbox0 extent (ckpc)
       $
       $ ; zoom runs only
       levelmin:       0,$     ; power of two minimum level parameter (e.g. MUSIC L7=128, L8=256, L9=512, L10=1024)
       levelmax:       0,$     ; power of two maximum level parameter (equals levelmin for non-zoom runs)
       zoomLevel:      0,$     ; levelmax-levelmin
       hInd:           0,$     ; zoom halo index (as in path)
       zoomShift:      [0,0,0]       ,$ ; Music output: "Domain will be shifted by (X, X, X)"
       zoomShiftPhys:  [0.0,0.0,0.0] ,$ ; the domain shift in box units
       targetHaloPos:  [0.0,0.0,0.0] ,$ ; position at targetRedshift in fullbox
       targetHaloInd:  0             ,$ ; hInd (subhalo index) at targetRedshift in fullbox
       targetHaloRvir: 0.0           ,$ ; rvir (ckpc) at targetRedshift
       targetHaloMass: 0.0           ,$ ; mass (logmsun) at targetRedshift
       targetRedshift: 0.0           ,$ ; maximum redshift the halo can be resimulated to
       ids_offset:     0L            ,$ ; IDS_OFFSET configuration parameter
       $
       $ ; tracers
       trMassConst:  0.0,$        ; mass per tracerMC under equal mass assumption (=TargetGasMass/trMCPerCell)
       trMCPerCell:  0,$          ; starting number of monte carlo tracers per cell
       trMCFields:   intarr(13),$ ; which TRACER_MC_STORE_WHAT fields did we save, and in what indices
       trVelPerCell: 0,$          ; starting number of velocity tracers per cell
       $
       $ ; analysis parameters:
       minNumGasPart: 0,$    ; minimum number of gas particles required to use subgroup (unused)
       radcut_rvir:   0.15,$ ; galcat: fraction of rvir as maximum for gal/stars, minimum for gmem (zero to disable)
       radcut_out:    1.5,$  ; galcat: fraction of rvir as maximum for gmem
       galcut_T:      6.0,$  ; galcat: temp coefficient for (rho,temp) galaxy cut
       galcut_rho:    0.25,$ ; galcat: dens coefficient for (rho,temp) galaxy cut
       rVirFacs:      [1.0,0.75,0.5,0.25,0.15,0.05,0.01],$ ; search for accretion times across these fractions of the virial radius
       radIndHaloAcc: 0,$     ; 1.0 rvir crossing for halo accretion
       radIndGalAcc:  4,$     ; 0.15 rvir crossing for galaxy accretion (or entering rho,temp definition)
       atIndMode:     -1,$     ; use first 1.0 rvir crossing to determine mode
       TcutVals:      [5.3,5.5,5.7],$ ; log(K) for constant threshold comparisons
       TvirVals:      [1.0,0.8,0.4],$ ; T/Tvir coefficients for variable threshold comparisons
       $
       $ ; plotting/vis parameters:
       colors:        [0L,0L,0L] ,$ ; color sequence for res 128,256,512
       pos_3x1:       list([0.18,0.67,0.95,0.95],[0.18,0.39,0.95,0.67],[0.18,0.11,0.95,0.39]) ,$
       pos_2x2:       list([0.13,0.5,0.53,0.9],[0.53,0.5,0.93,0.9],[0.13,0.1,0.53,0.5],[0.53,0.1,0.93,0.5]) ,$
       pos_2x1:       list([0.15,0.55,0.95,0.95],[0.15,0.15,0.95,0.55])  ,$
       $
       $ ; GFM and other indications of optional snapshot fields
       gfmElements:   ['H','He','C','N','O','Ne','Mg','Si','Fe'] ,$
       gfmNumElements: 0, $ ; set to >=1 for GFM runs outputting abundances by metal
       gfmBHs:         0, $ ; set to 1 for BLACK_HOLES
       gfmWinds:       0  $ ; set to 1 for GFM_WINDS
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
 
  ; illustris
  if (run eq 'illustris') or (run eq 'illustris_dm') then begin
    if res ne 455 and res ne 910 and res ne 1820 then message,'Error: Invalid res.'
    
    r.minNumGasPart  = -1 ; no additional cut
    r.trVelPerCell   = 0
    r.trMCPerCell    = 1
    r.trMCFields     = [0,1,2,3,4,5,6,7,8,9,10,11,12] ; all (=4096, 13 of 13)
    r.gfmNumElements = 9
    r.gfmWinds       = 1
    r.gfmBHs         = 1
    
    ; TODO: subbox (4 of them)
    ; TODO: SPT / other runs (and some 40Mpc boxes?)
    
    r.boxSize       = 75000.0
    r.snapRange     = [0,135]
    r.groupCatRange = [21,135] ; z6=45, z5=49, z4=54, z3=60, z2=68, z1=85, z0=135
    
    if res eq 455  then r.targetGasMass = 5.66834e-3
    if res eq 910  then r.targetGasMass = 7.08542e-4
    if res eq 1820 then r.targetGasMass = 8.85678e-5
	
    r.trMassConst = r.targetGasMass / r.trMCPerCell
    
    if res eq 455  then r.gravSoft = 4.0 ; comoving until z=1, then fixed at z=1 value after (except DM)
    if res eq 910  then r.gravSoft = 2.0 ; same
    if res eq 1820 then r.gravSoft = 1.0 ; same
    
    if run eq 'illustris' then begin
      r.simPath    = '/n/home07/dnelson/sims.illustris/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_FP/output/'
      r.arepoPath  = '/n/home07/dnelson/sims.illustris/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_FP/'
      r.savPrefix  = 'I'
      r.simName    = 'ILLUSTRIS'
      r.saveTag    = 'il'
      r.plotPrefix = 'il'
      r.plotPath   = '/n/home07/dnelson/plots/'
      r.derivPath  = '/n/home07/dnelson/sims.illustris/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_FP/data.files/'
    endif
    
    if run eq 'illustris_dm' then begin
      r.simPath    = '/n/home07/dnelson/sims.illustris/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_DM/output/'
      r.arepoPath  = '/n/home07/dnelson/sims.illustris/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_DM/'
      r.savPrefix  = 'IDM'
      r.simName    = 'ILLUSTRIS_DM'
      r.saveTag    = 'ilDM'
      r.plotPrefix = 'ilDM'
      r.plotPath   = '/n/home07/dnelson/plots/'
      r.derivPath  = '/n/home07/dnelson/sims.illustris/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_DM/data.files/'
    endif
  
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return, r  
  endif
  
  ; zoom project: DM only single halo zooms (L=7/128 for fullbox)
  if run eq 'zoom_20mpc_dm' then begin
  
    r.minNumGasPart  = -1 ; no additional cut
    r.trVelPerCell   = 0
    r.trMCPerCell    = 0
    r.trMCFields     = replicate(-1,13)
    r.gfmNumElements = 0
    r.gfmWinds       = 0
    r.gfmBHs         = 0
    r.boxSize        = 20000.0
    r.levelMin       = 7 ; uniform box @ 128
    r.levelMax       = 7 ; default
    r.snapRange      = [0,10] ; z99=0, z0=10
    r.groupCatRange  = [2,10]
    r.targetGasMass  = 0.0
    r.trMassConst    = 0.0
    
    if n_elements(hInd) gt 0 then r = fillZoomParams(r,res=res,hInd=hInd)
              
    if r.levelMin ne r.levelMax then $
      pathStr = '128_' + str(fix(r.boxSize/1000)) + 'Mpc_h' + str(hInd) + '_L' + str(r.levelMax) $
    else $
      pathStr = '128_' + str(fix(r.boxSize/1000)) + 'Mpc'
              
    r.simPath    = '/n/home07/dnelson/sims.zooms/'+pathStr+'_dmonly/output/'
    r.arepoPath  = '/n/home07/dnelson/sims.zooms/'+pathStr+'_dmonly/'
    r.savPrefix  = 'Z'
    r.simName    = 'ZOOM_L'+str(r.levelMax)+'_DM'
    r.saveTag    = 'zDmL'+str(r.levelMax)
    r.plotPrefix = 'zDmL'+str(r.levelMax)
    r.plotPath   = '/n/home07/dnelson/plots/'
    r.derivPath  = '/n/home07/dnelson/sims.zooms/'+pathStr+'_dmonly/data.files/'
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return, r
  endif
  
  ; zoom project: DM+gas single halo zooms (now all in 20Mpc box, add boxsize to run label later)
  if run eq 'zoom_20mpc' or run eq 'zoom_20mpc_derefgal' or run eq 'zoom_20mpc_derefgal_nomod' then begin
    if n_elements(hInd) eq 0 then message,'Error: Must specify hInd (halo index) to load zoom.'
    
    r.minNumGasPart  = -1 ; no additional cut
    r.trVelPerCell   = 0
    r.trMCPerCell    = 5
    r.trMCFields     = [0,1,2,3,4,5,6,7,8,9,-1,-1,-1] ; up to and including WIND_COUNTER (=512, 10 of 13)
    r.gfmNumElements = 0
    r.gfmWinds       = 0
    r.gfmBHs         = 0
    r.boxSize        = 20000.0
    r.levelMin       = 7 ; uniform box @ 128
    r.levelMax       = 7 ; default
    r.snapRange      = [0,59]
    r.groupCatRange  = [5,59] ; z6=5, z5=14, z4=21, z3=36, z2=59
    
    if n_elements(hInd) gt 0 then r = fillZoomParams(r,res=res,hInd=hInd)
              
    if r.levelMin ne r.levelMax then $
      pathStr = '128_' + str(fix(r.boxSize/1000)) + 'Mpc_h' + str(hInd) + '_L' + str(r.levelMax) $
    else $
      pathStr = '128_' + str(fix(r.boxSize/1000)) + 'Mpc'
    
    r.trMassConst = r.targetGasMass / r.trMCPerCell
    
    if run eq 'zoom_20mpc_derefgal' then pathStr = pathStr + '_derefgal'
    if run eq 'zoom_20mpc_derefgal_nomod' then pathStr = pathStr + '_derefgal'
    
    dirStr = ''
    if run eq 'zoom_20mpc_derefgal_nomod' then dirStr = '_noMod'
        
    r.simPath    = '/n/home07/dnelson/sims.zooms/'+pathStr+'/output'+dirStr+'/'
    r.arepoPath  = '/n/home07/dnelson/sims.zooms/'+pathStr+'/'
    r.savPrefix  = 'Z'
    r.simName    = 'ZOOM_L'+str(r.levelMax)
    r.saveTag    = 'zL'+str(r.levelMax)
    r.plotPrefix = 'zL'+str(r.levelMax)
    r.plotPath   = '/n/home07/dnelson/plots/'
    r.derivPath  = '/n/home07/dnelson/sims.zooms/'+pathStr+'/data.files'+dirStr+'/'
    
    if run eq 'zoom_20mpc_derefgal' then begin
      r.saveTag = r.saveTag + 'drG'
      r.simName = r.simName + '_drG'
    endif
    if run eq 'zoom_20mpc_derefgal_nomod' then begin
      r.saveTag = r.saveTag + 'drN'
      r.simName = r.simName + '_drN'
    endif
    
    r.colors = [getColor24(['ff'x,'40'x,'40'x]), $ ; red, light to dark (L9/default)
                getColor24(['ff'x,'00'x,'00'x]), $
                getColor24(['a6'x,'00'x,'00'x])]
                
    if r.res eq 10 then $
      r.colors = [getColor24(['ff'x,'5b'x,'00'x]), $ ; red, light to dark (L9/default)
                  getColor24(['ff'x,'7c'x,'00'x]), $
                  getColor24(['ff'x,'9d'x,'00'x])]
                  
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return, r
  endif
 
  ; sims.feedback (Velocity + f=5 Monte Carlo) 128,256,512 @ 20Mpc w/ fiducial Illustris parameters
  if (run eq 'feedback') or (run eq 'feedback_noz') or (run eq 'feedback_nofb') or (run eq 'feedback_nogfm') then begin
    r.minNumGasPart  = -1 ; no additional cut
    r.trVelPerCell   = 0
    r.trMCPerCell    = 5
    r.trMCFields     = [0,1,2,3,4,5,6,7,8,9,-1,-1,-1] ; up to and including WIND_COUNTER (=512, 10 of 13)
    r.gfmNumElements = 9
    r.gfmWinds       = 1
    r.gfmBHs         = 1
    
    if res eq 512 then r.subboxCen  = [5500,7000,7500]
    if res eq 512 then r.subboxSize = [4000,4000,4000]
    
    if res eq 128 or res eq 256 or res eq 512 then begin 
      r.boxSize       = 20000.0
      r.snapRange     = [0,133]
      r.groupCatRange = [5,133] ; z6=5, z5=14, z4=21, z3=36, z2=60, z1=81, z0=130
    endif else begin
      message,'res error'
    endelse
    
    if res eq 128 then r.targetGasMass = 4.76446157e-03
    if res eq 256 then r.targetGasMass = 5.95556796e-04
    if res eq 512 then r.targetGasMass = 7.44447120e-05
	
    r.trMassConst = r.targetGasMass / r.trMCPerCell
    
    if res eq 128 then r.gravSoft = 4.0
    if res eq 256 then r.gravSoft = 2.0
    if res eq 512 then r.gravSoft = 1.0
  
    r.simPath    = '/n/home07/dnelson/sims.feedback/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/output/'
    r.arepoPath  = '/n/home07/dnelson/sims.feedback/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/'
    r.savPrefix  = 'F'
    r.simName    = 'FEEDBACK'
    r.saveTag    = 'feMC'
    r.plotPrefix = 'feMC'
    r.plotPath   = '/n/home07/dnelson/plots/'
    r.derivPath  = '/n/home07/dnelson/sims.feedback/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/data.files/'
    
    r.colors = [getColor24(['3e'x,'41'x,'e8'x]), $ ; blue, light to dark
                getColor24(['15'x,'1a'x,'ac'x]), $
                getColor24(['08'x,'0b'x,'74'x])]
    
    ; if f=-1 use velocity tracers
    if keyword_set(f) then message,'Error' ; no velocity tracers
	
    ; if noZ or noFB runs, modify details and paths
    if run eq 'feedback_noz' then begin
      if res ne 256 then stop ; only 256 exists
      r.trMCFields     = [0,1,2,3,4,5,6,7,8,9,-1,-1,-1]
      
      r.simPath    = '/n/home07/dnelson/sims.feedback/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_noZ/output/'
      r.arepoPath  = '/n/home07/dnelson/sims.feedback/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_noZ/'
      r.saveTag    = 'feNoZ'
      r.simName    = 'AREPO noZ'
      r.plotPrefix = 'feNoZ'
      r.derivPath  = '/n/home07/dnelson/sims.feedback/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_noZ/data.files/'
      
      r.colors = [getColor24(['b2'x,'3a'x,'d4'x]), $ ; purple, light to dark
                  getColor24(['a5'x,'26'x,'c9'x]), $
                  getColor24(['66'x,'12'x,'7e'x])]
    endif
	
    if run eq 'feedback_nofb' then begin
      if res ne 256 then stop ; only 256 exists
      r.trMCFields     = [0,1,2,3,4,5,6,7,8,-1,-1,-1,-1] ; up to and including LST (=256, 9 of 13)
      
      r.gfmWinds   = 0
      r.simPath    = '/n/home07/dnelson/sims.feedback/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_noFB/output/'
      r.arepoPath  = '/n/home07/dnelson/sims.feedback/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_noFB/'
      r.saveTag    = 'feNoFB'
      r.simName    = 'AREPO noFB'
      r.plotPrefix = 'feNoFB'
      r.derivPath  = '/n/home07/dnelson/sims.feedback/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_noFB/data.files/'
      
      r.colors = [getColor24(['ff'x,'40'x,'40'x]), $ ; red, light to dark
                  getColor24(['ff'x,'00'x,'00'x]), $
                  getColor24(['a6'x,'00'x,'00'x])]
    endif
    
    if run eq 'feedback_nogfm' then begin
      if res ne 256 then stop ; only 256 exists
      r.trMCFields     = [0,1,2,3,4,5,6,7,-1,-1,-1,-1,-1] ; up to and including LST (=128, 8 of 13)
      
      r.gfmWinds   = 0
      r.simPath    = '/n/home07/dnelson/sims.feedback/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_noGFM/output/'
      r.arepoPath  = '/n/home07/dnelson/sims.feedback/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_noGFM/'
      r.saveTag    = 'feNoGFM'
      r.simName    = 'AREPO noGFM'
      r.plotPrefix = 'feNoGFM'
      r.derivPath  = '/n/home07/dnelson/sims.feedback/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc_noGFM/data.files/'
      
      r.colors = [getColor24(['00'x,'a7'x,'79'x]), $ ; greenblue, light to dark
                  getColor24(['1f'x,'7d'x,'63'x]), $
                  getColor24(['00'x,'6d'x,'4f'x])]
    endif
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
 
  ; ComparisonProject GADGET 128,256,512 @ 20Mpc
  if (run eq 'gadget') then begin
    r.minNumGasPart = -1 ; no additional cut
    r.trMCPerCell   = 0  ; none (SPH)
    r.trMCFields    = intarr(13)-1 ; none (SPH)
    
    if keyword_set(f) then stop ; shouldn't be specified
    if res ne 128 and res ne 256 and res ne 512 then stop ; only resolutions that exist now
  
    if res eq 128 or res eq 256 or res eq 512 then begin 
      r.boxSize       = 20000.0
      r.snapRange     = [0,314]
      r.groupCatRange = [50,314]
    endif
  
    if res eq 128 then r.targetGasMass = 4.76446157e-03
    if res eq 256 then r.targetGasMass = 5.95556796e-04
    if res eq 512 then r.targetGasMass = 7.44447120e-05
    
    if res eq 128 then r.gravSoft = 4.0
    if res eq 256 then r.gravSoft = 2.0
    if res eq 512 then r.gravSoft = 1.0
    
    r.simPath    = '/n/home07/dnelson/sims.gadget/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/output/'
    r.arepoPath  = '/n/home07/dnelson/sims.gadget/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/'
    r.savPrefix  = 'G'
    r.saveTag    = 'SPH'
    r.simName    = 'GADGET'
    r.plotPrefix = 'G'
    r.plotPath   = '/n/home07/dnelson/plots/'
    r.derivPath  = '/n/home07/dnelson/sims.gadget/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/data.files/'
    
    r.colors = [getColor24(['e6'x,'7a'x,'22'x]),$ ; brown, light to dark
                getColor24(['b3'x,'5f'x,'1b'x]),$
                getColor24(['80'x,'44'x,'13'x])]
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  ; sims.tracers (Velocity + f=10 Monte Carlo) 128,256,512 @ 20Mpc
  if (run eq 'tracer') or (run eq 'tracer_nouv') then begin
    r.minNumGasPart = -1 ; no additional cut
    r.trMCPerCell   = 10
    r.trVelPerCell  = 1
    
    if res eq 128 or res eq 256 then $
      r.trMCFields    = [0,1,5,2,-1,3,4,-1,-1,-1,-1,-1,-1] ; even older code version than tracer.512, indices specified manually in Config.sh
    if res eq 512 then $
      r.trMCFields    = [0,1,4,2,-1,3,5,-1,-1,-1,-1,-1,-1] ; up to and including ENTMAX but in older code, ordering permuted (6 vals, CAREFUL)
    
    if res ne 128 and res ne 256 and res ne 512 then stop
    
    if res eq 128 or res eq 256 or res eq 512 then begin 
      r.boxSize       = 20000.0
      r.snapRange     = [0,314]
      r.groupCatRange = [50,314]
    endif
    
    if res eq 512 then r.subboxCen  = [5500,7000,7500]
    if res eq 512 then r.subboxSize = [4000,4000,4000]
    
    if res eq 128 then r.targetGasMass = 4.76446157e-03
    if res eq 256 then r.targetGasMass = 5.95556796e-04
    if res eq 512 then r.targetGasMass = 7.44447120e-05
    
    r.trMassConst = r.targetGasMass / r.trMCPerCell
    
    if res eq 128 then r.gravSoft = 4.0
    if res eq 256 then r.gravSoft = 2.0
    if res eq 512 then r.gravSoft = 1.0
  
    r.simPath    = '/n/home07/dnelson/sims.tracers/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/output/'
    r.arepoPath  = '/n/home07/dnelson/sims.tracers/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/'
    r.savPrefix  = 'N'
    r.saveTag    = 'trMC'
    r.simName    = 'AREPO'
    r.plotPrefix = 'trMC'
    r.plotPath   = '/n/home07/dnelson/plots/'
    r.derivPath  = '/n/home07/dnelson/sims.tracers/'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/data.files/'
    
    ; originally: 00eb00,00bd00,009000 for gas accretion paper
    r.colors = [getColor24(['00'x,'ab'x,'33'x]), $ ; green, light to dark
                getColor24(['00'x,'7d'x,'23'x]), $
                getColor24(['00'x,'90'x,'13'x])]
    
    ; if f=-1 use velocity tracers
    if keyword_set(f) then begin
      if f ne -1 then message,'Error.' ; only valid input is -1
      r.trMCPerCell = -1
      r.plotPrefix = 'trVel'
      r.saveTag    = 'trVel'
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
	  r.simName    = 'AREPO noUV'
      r.plotPrefix = 'trNoUV'
      r.plotPath   = '/n/home07/dnelson/plots/'
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
    r.trMCFields    = intarr(13)-1 ; none
    
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
	r.saveTag    = 'AR'
	r.simName    = 'AREPO'
    r.plotPrefix = 'A'
    r.plotPath   = '/n/home07/dnelson/plots/'
    r.derivPath  = '/n/home07/dnelson/sims.tracers/A'+str(res)+'_'+str(fix(r.boxSize/1000))+'Mpc/data.files/'
    
    ; if redshift passed in, convert to snapshot number and save
    if (n_elements(redshift) eq 1) then r.snap = redshiftToSnapNum(redshift,sP=r)
    
    return,r
  endif
  
  message,'simParams: ERROR.'

end
