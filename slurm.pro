; slurm.pro
; automation and queue job related
; dnelson sep.2013

; makeArepoFoFBsub(): write slurm job file to invoke FoF/Subfind group finder post-processing

pro makeArepoFoFBsub

  ; config
  sP = simParams(res=512,run='gadget')

  snapRange = [251,314,1]
  ;snaps = [257]

  ; job config
  spawnJobs = 1   ; execute sbatch?
  nSimul    = 4   ; how many maximum to execute simultaneously? 0=all
  nProcs    = 128 ; needed nodes: 512^3 use 128, <=256^3 use 64
  cmdCode   = 3   ; fof/substructure post process
  partName  = "hernquist"
  paramFile = "param_fof.txt"
  
  count = 0
  if nSimul gt 0 then prevJobIDs = lonarr(nSimul)

  for snap=snapRange[0],snapRange[1],snapRange[2] do begin
  ;foreach snap,snaps do begin
  
    ; write bjob file
    jobFileName = sP.plotPath + 'job_fof.slurm'
    
    ; check before overriding
    if file_test(jobFileName) then message,'Error: job_fof.slurm already exists'
    
    print,'['+string(count,format='(I3.3)')+'] snap ['+string(snap,format='(I3.3)')+'] Writing: '+jobFilename
          
    openW, lun, jobFilename, /GET_LUN
    
    ; write header
    printf,lun,'#!/bin/sh'
    printf,lun,'#SBATCH --mail-user dnelson@cfa.harvard.edu'
    printf,lun,'#SBATCH --mail-type=fail'
    printf,lun,'#SBATCH -p ' + partName
    printf,lun,'#SBATCH -J fof_' + str(snap) + ''
    printf,lun,'#SBATCH -o run.out'
    printf,lun,'#SBATCH -e run.err'
    printf,lun,'#SBATCH -t 360 # 6 hours in min'
    printf,lun,'#SBATCH --mem=236800 # 3700MB/core'
    printf,lun,'#SBATCH --exclusive'
    printf,lun,'#SBATCH --ntasks ' + str(nProcs)
    printf,lun,'#SBATCH --workdir ' + sP.arepoPath
    
    ; enforce maximum simultaneous jobs by adding a dependency on a previously subbmited jobID
    if nSimul gt 0 and count ge nSimul then $
      printf,lun,'#SBATCH --dependency=afterany:'+str(prevJobIDs[ count mod nSimul ])

    printf,lun,''
      
    ; write projection commands
    strArray = ['mpirun --mca mpi_leave_pinned 0 -n '+str(nProcs),$
               './Arepo_fof',$
                paramFile,$
                str(cmdCode),$
                str(snap),$
                '>> run_fof_'+str(snap)+'.txt']
    printf,lun,strjoin(strArray,' ')
      
    ; close
    close,lun
    free_lun,lun
      
    ; add to queue if requested
    if (spawnJobs) then begin
      spawn, 'sbatch < ' + sP.plotPath + 'job_fof.slurm', result
      print,'  '+result
      
      ; check for error
      if strpos(result,"sbatch: error") ge 0 then message,'Error detected, paused.'
      
      ; file cleanup
      wait,0.2
      spawn, 'rm ' + sP.plotPath + 'job_fof.slurm'
    endif
    
    if nSimul gt 0 then begin
      resSplit = strsplit(result, " ", /extract)
      prevJobIDs[ count mod nSimul ] = long(resSplit[3])
    endif
    
    count += 1
  
  endfor ;snapRange
  ;endforeach

end

; makeArepoProjBsub(): write bsub file to invoke three different axis aligned projections using the 
;                      Arepo voronoi_makeimage_new() code

pro makeArepoProjBsub

  ; config - path and snapshot
  ;res = 256
  ;run = 'tracerMC.ref'
  ;f   = '1'
  ;subBox = 0 ; search for subbox snapshot versus normal
  
  snapRange = [0,20,1]
  
  ; config - viewsize / object
  ;sP = simParams(res=res,run=run,f=f,snap=snapRange[0])
  sP = { simPath   : '/n/home07/dnelson/dev.tracerMC/gasSphere.cylTest.1e4.norot.nocool.nosg/output/' ,$
         arepoPath : '/n/home07/dnelson/dev.tracerMC/gasSphere.cylTest.1e4.norot.nocool.nosg/',$
         snap      : 0 }
  h = loadSnapshotHeader(sP=sP,subBox=subBox)
  
  ;boxSize = h.boxSize ; full box
  boxSize = 5000.0 ; subbox
  
  ;xyzCen = [1123.20,7568.80,16144.2] ;dusan 512
  xyzCen = [boxSize/2.0,boxSize/2.0,boxSize/2.0]

  sliceWidth  = boxSize ; cube sidelength
  zoomFac     = 10.0   ; final image is boxSize/zoomFac in extent

  ; render config
  spawnJobs   = 1    ; execute bsub?
  nProcs      = 1    ; 128^3=8, 256^3=24, 512^3=256 (minimum for memory)
  ptile       = 1
  dimX        = 1000  ; image dimensions (x pixels)
  dimY        = 1000  ; image dimensions (y pixels)
  
  axesStr = ['0 1 2','0 2 1','1 2 0'] ;xy,xz,yz
  
  ; bbox and projection setup
  cmdCode = 5 ;projection
 
  paramFile = "param.txt"
  
  xMin = xyzCen[0] - sliceWidth/2.0/zoomFac
  xMax = xyzCen[0] + sliceWidth/2.0/zoomFac
  yMin = xyzCen[1] - sliceWidth/2.0/zoomFac
  yMax = xyzCen[1] + sliceWidth/2.0/zoomFac
  zMin = xyzCen[2] - sliceWidth/2.0/zoomFac
  zMax = xyzCen[2] + sliceWidth/2.0/zoomFac
  
  ; integration range for each axis depends on
  axesBB = [[xMin,xMax,yMin,yMax,zMin,zMax],$
            [xMin,xMax,zMin,zmax,yMin,yMax],$
            [yMin,yMax,zMin,zMax,xMin,xMax]]
 
  for snap=snapRange[0],snapRange[1],snapRange[2] do begin
 
    ; write bjob file
    jobFileName = sP.arepoPath+'job_proj.bsub'
    
    ; check before overriding
    if file_test(jobFileName) then message,'Error: job_proj.bsub already exists'
    
    print,'['+str(snap)+'] Writing: '+jobFilename
      
    openW, lun, jobFilename, /GET_LUN
    
    ; write header
    printf,lun,'#!/bin/sh'
    printf,lun,'#BSUB -q keck'
    printf,lun,'#BSUB -J proj_' + str(snap) + ""
    printf,lun,'#BSUB -n ' + str(nProcs)
    printf,lun,'#BSUB -R "rusage[mem=30000] span[ptile=' + str(ptile) + ']"'
    printf,lun,'#BSUB -o run_proj.out'
    printf,lun,'#BSUB -g /dnelson/proj' ; 4 concurrent jobs limit automatically
    printf,lun,'#BSUB -cwd ' + sP.arepoPath
    printf,lun,'#BSUB -a openmpi'
    printf,lun,''
      
    ; write projection commands
    foreach axisStr, axesStr, i do begin 
    
      strArray = ['mpirun.lsf ./Arepo_proj '+paramFile,$
                  str(cmdCode),$
                  str(snap),$
                  str(dimX),str(dimY),$
                  axisStr,$
                  str(axesBB[0,i]),str(axesBB[1,i]),$
                  str(axesBB[2,i]),str(axesBB[3,i]),$
                  str(axesBB[4,i]),str(axesBB[5,i]),$
                  '>> run_proj.txt']
      printf,lun,strjoin(strArray,' ')
      
      ; write file handling commands
      outputFilename = 'proj_density_field_' + string(snap,format='(I3.3)') + '.' + str(i) + '.dat'
      
      printf,lun,''
      printf,lun,'mv ' + sP.simPath + 'proj_density_field_' + string(snap,format='(I3.3)') + ' ' + $
                 sP.simPath + outputFilename
      printf,lun,''
    
    endforeach
      
    ; close
    close,lun
    free_lun,lun
      
    ; add to queue if requested
    if (spawnJobs) then begin
      spawn, 'bsub < ' + sP.arepoPath + 'job_proj.bsub', result
      print,'  '+result
      
      ; file cleanup
      wait,0.5
      spawn, 'rm ' + sP.arepoPath + 'job_proj.bsub'
    endif
  
  endfor ;snapRange

end

; runGasAccPaper(): run all analysis and make all plots (run fof+subfind first)

pro runGasAccPaper

  ; config
  redshifts   = [2.0] ;0.0,1.0,3.0
  resolutions = [128,256,512]
  runs        = ['gadget','tracer'] ;'feedback'
  
  timeWindows = list('all',1000.0,250.0) ; Myr
  accModes    = ['all','smooth','clumpy','stripped']
  
  foreach res,resolutions do begin
    foreach run,runs do begin
      foreach redshift,redshifts do begin
      
        ; simulation parameters
        sP = simParams(res=res,run=run,redshift=redshift)
        snap = sP.snap
        
        ; catalogs at analysis redshift
        x = galaxyCat(sP=sP)
        
        ; primary analysis
        x = maxVals(sP=sP)
        x = mergerTree(sP=sP,makeNum=snap)
        x = mergerTreeSubset(sP=sP)
        x = accretionTimes(sP=sP)
        x = accretionMode(sP=sP)
        
        ; binning for plots
        foreach timeWindow,timeWindows do begin
          foreach accMode,accModes do begin
            x = haloMassBinValues(sP=sP,timeWindow=timeWindow,accMode=accMode)
            x = binValMaxHistos(sP=sP,timeWindow=timeWindow,accMode=accMode)
          endforeach
        endforeach
        
      endforeach
    endforeach
  endforeach
        
  ; line plots
  plotByMode
  plotByMethod
  plotByRes
  plotValMaxHistos
  
  ; visual plots
  scatterMapHalosComp
  mosaicHalosComp
  plotHaloShellValueComp
  
end
