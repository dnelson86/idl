; slurm.pro
; automation and queue job related
; dnelson oct.2013

; runArepoFoF(): write slurm job file to invoke FoF/Subfind group finder post-processing

pro runArepoFoF

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

; runArepoProj(): write job file to invoke a projection, slice, etc using the
;                 Arepo voronoi_makeimage_new() code

function runArepoProj, sP=sP, nPixels=nPixels, axes=axes, xyzCen=xyzCen, xyzSize=xyzSize
  compile_opt idl2, hidden, strictarr, strictarrsubs
  if n_elements(sP) eq 0 or n_elements(xyzCen) ne 3 or n_elements(nPixels) ne 2 or $
     n_elements(axes) ne 3 or n_elements(xyzSize) ne 3 then message,'Error'
  
  ; config
  subBox      = 0    ; unused
  spawnJobs   = 0    ; execute job?
  nCores      = 64   ; number of cores
  cmdCode     = 5    ; projection
  partName    = 'hernquist' ; partition to use
  paramFile   = "param_proj.txt"
  
  xMin = xyzCen[0] - xyzSize[0]*0.5
  xMax = xyzCen[0] + xyzSize[0]*0.5
  yMin = xyzCen[1] - xyzSize[1]*0.5
  yMax = xyzCen[1] + xyzSize[1]*0.5
  zMin = xyzCen[2] - xyzSize[2]*0.5
  zMax = xyzCen[2] + xyzSize[2]*0.5
  
  outputFilename = 'proj.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + $
    '.cen'  + str(round(xyzCen[0]))  + '-' + str(round(xyzCen[1]))  + '-' + str(round(xyzCen[2])) + $
    '.size' + str(round(xyzSize[0])) + '-' + str(round(xyzSize[1])) + '-' + str(round(xyzSize[2])) + $
    '.px'   + str(nPixels[0]) + '-' + str(nPixels[1]) + $
    '.axes' + str(axes[0]) + str(axes[1]) + str(axes[2]) + $
    '.dat'
  
  ; if the requested projection is already done, load it and return
  if file_test(outputFilename) then begin
    openr,1,outputFilename
      ; read header
      nPixelsX = 0L
      nPixelsY = 0L

      readu,1,nPixelsX
      readu,1,nPixelsY
      
      dens = fltarr(nPixelsY, nPixelsX)
      temp = fltarr(nPixelsY, nPixelsX)
      readu,1,dens
      readu,1,temp
    close,1

    dens = transpose(dens)
    temp = transpose(temp)
  
    r = {nPixels:[nPixelsX,nPixelsY],dens:dens,temp:temp}
    return,r
  endif
  
  ; write bjob file
  jobFileName = sP.arepoPath + 'job_proj.slurm'
    
  ; check before overriding
  if file_test(jobFileName) then message,'Error: job_proj.slurm already exists'
    
  print,'['+str(sP.snap)+'] Writing: '+jobFilename
      
  openW, lun, jobFilename, /GET_LUN
    
    ; write header
    printf,lun,'#!/bin/sh'
    printf,lun,'#SBATCH --mail-user dnelson@cfa.harvard.edu'
    printf,lun,'#SBATCH --mail-type=fail'
    printf,lun,'#SBATCH -p ' + partName
    printf,lun,'#SBATCH -J proj_' + str(sP.snap) + ''
    printf,lun,'#SBATCH -o run_proj.txt'
    printf,lun,'#SBATCH -e run_proj.txt'
    printf,lun,'#SBATCH -t 360 # 6 hours in min'
    printf,lun,'#SBATCH --mem=243200 # 3800MB/core'
    printf,lun,'#SBATCH --exclusive'
    printf,lun,'#SBATCH --ntasks ' + str(nCores)
    printf,lun,'#SBATCH --workdir ' + sP.arepoPath
      
    printf,lun,''
      
    ; write projection commands
    strArray = ['mpirun --mca mpi_leave_pinned 0 -n '+str(nCores),$
                './Arepo_proj'    ,$
                paramFile        ,$
                str(cmdCode)     ,$
                str(sP.snap)     ,$
                str(nPixels[0])  ,$
                str(nPixels[1])  ,$ 
                str(axes[0])     ,$
                str(axes[1])     ,$
                str(axes[2])     ,$
                str(xMin),str(xMax),$
                str(yMin),str(yMax),$
                str(zMin),str(zMax),$
                '>> run_proj.txt']
    printf,lun,strjoin(strArray,' ')
      
    ; write file handling commands      
    ;printf,lun,''
    ;printf,lun,'mv ' + sP.simPath + 'proj_density_field_' + string(sP.snap,format='(I3.3)') + ' ' + $
    ;           sP.derivPath + 'binnedVals/' + outputFilename
    ;printf,lun,''
      
  ; close
  close,lun
  free_lun,lun
      
  ; add to queue if requested
  if (spawnJobs) then begin
    spawn, 'sbatch < ' + sP.plotPath + 'job_proj.slurm', result
    print,'  '+result
      
    ; check for error
    if strpos(result,"sbatch: error") ge 0 then message,'Error detected, paused.'
      
    ; file cleanup
    wait,0.2
    spawn, 'rm ' + sP.plotPath + 'job_proj.slurm'
  endif
  
  return,-1

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
