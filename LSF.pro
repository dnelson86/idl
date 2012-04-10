; LSF.pro
; automation and queue job related
; dnelson apr.2012

; makeArepoFoFBsub(): write bsub file to invoke FoF/Subfind group finder post-processing

pro makeArepoFoFBsub

  ; config
  res = 512
  run = 'gadget'
  ;f = '10'
  
  sP = simParams(res=res,run=run)

  snapRange = [314,314,1]

  ; job config
  spawnJobs = 1 ; execute bsub?
  nProcs    = 64 ; needed nodes: 128^3=0.5 (n4tile4 PartAllocFactor=2)
                 ;               256^3=4 (n32tile4 PartAllocFactor=2.5 MaxMemSize=7900)
                 ;               512^3=16 (n64tile4 PartAllocFactor=1.5 though n32tile4 ok until z=2)
  ptile     = 2 ; span[ptile=X]
  cmdCode   = 3 ; fof/substructure post process
  
  excludeHosts = ['hero2402'] ;leave empty otherwise
 
  paramFile = "param_fof.txt"

  for snap=snapRange[0],snapRange[1],snapRange[2] do begin
 
    ; write bjob file
    jobFileName = sP.plotPath + 'job_fof.bsub'
    
    ; check before overriding
    if file_test(jobFileName) then begin
      print,'Error: job_fof.bsub already exists'
      return
    endif
    
    print,'['+str(snap)+'] Writing: '+jobFilename
    
    ; host exclusion string
    selectStr = ''
    if n_elements(excludeHosts) gt 0 then begin
      selectStr = 'select['
      for i=0,n_elements(excludeHosts)-1 do selectStr += 'hname!='+excludeHosts[i]+' '
      selectStr = strmid(selectStr,0,strlen(selectStr)-1) + '] '
    endif
      
    openW, lun, jobFilename, /GET_LUN
    
    ; write header
    printf,lun,'#!/bin/sh'
    printf,lun,'#BSUB -q keck'
    printf,lun,'#BSUB -J fof_' + str(snap) + ''
    printf,lun,'#BSUB -n ' + str(nProcs)
    printf,lun,'#BSUB -R "' + selectStr + 'span[ptile=' + str(ptile) + ']"' ;rusage[mem=31000]
    printf,lun,'#BSUB -x'
    printf,lun,'#BSUB -o run_fof_'+str(snap)+'.out'
    printf,lun,'#BSUB -g /dnelson/fof' ; 4 concurrent jobs limit automatically
    printf,lun,'#BSUB -cwd ' + sP.arepoPath
    printf,lun,'#BSUB -a openmpi'
    printf,lun,''
      
    ; write projection commands
    strArray = ['mpirun.lsf ./Arepo_fof '+paramFile,$
                str(cmdCode),$
                str(snap),$
                '>> run_fof_'+str(snap)+'.txt']
    printf,lun,strjoin(strArray,' ')
      
    ; close
    close,lun
    free_lun,lun
      
    ; add to queue if requested
    if (spawnJobs) then begin
      spawn, 'bsub < ' + sP.plotPath + 'job_fof.bsub', result
      print,'  '+result
      
      ; file cleanup
      wait,0.5
      spawn, 'rm ' + sP.plotPath + 'job_fof.bsub'
    endif
  
  endfor ;snapRange

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
    if file_test(jobFileName) then begin
      print,'Error: job_proj.bsub already exists'
      return
    endif
    
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
