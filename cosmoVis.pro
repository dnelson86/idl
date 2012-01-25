; cosmoVis.pro
; cosmological boxes - 2d visualization
; dnelson jan.2012

; makeArepoFoFBsub(): write bsub file to invoke FoF/Subfind group finder post-processing

pro makeArepoFoFBsub

  ; config
  res = 128
  run = 'dev.tracer.nonrad'
  
  sP = simParams(res=res,run=run)
  
  ;redshift = 3.0
  ;snap     = redshiftToSnapNum(redshift,sP=sP)
  
  snapRange = [0,75,1]
  ;snapRange = [25,63,1]

  ; job config
  spawnJobs = 1 ; execute bsub?
  nProcs    = 8 ; 
  ptile     = 8 ; span[ptile=X]
  cmdCode   = 3 ; fof/substructure post process
 
  paramFile = "param.txt"

  for snap=snapRange[0],snapRange[1],snapRange[2] do begin
 
    ; write bjob file
    jobFileName = sP.plotPath + 'job_fof.bsub'
    
    ; check before overriding
    if file_test(jobFileName) then begin
      print,'Error: job_fof.bsub already exists'
      return
    endif else begin
      print,'['+str(snap)+'] Writing: '+jobFilename
    endelse
      
    openW, lun, jobFilename, /GET_LUN
    
    ; write header
    printf,lun,'#!/bin/sh'
    printf,lun,'#BSUB -q keck'
    printf,lun,'#BSUB -J fof_' + str(snap)
    printf,lun,'#BSUB -n ' + str(nProcs)
    printf,lun,'#BSUB -R "rusage[mem=30000] span[ptile=' + str(ptile) + ']"'
    printf,lun,'#BSUB -x'
    printf,lun,'#BSUB -o run_fof.out'
    printf,lun,'#BSUB -g /dnelson/fof' ; 4 concurrent jobs limit automatically
    printf,lun,'#BSUB -cwd ' + sP.arepoPath
    printf,lun,'#BSUB -a openmpi'
    printf,lun,''
      
    ; write projection commands
    strArray = ['mpirun.lsf ./Arepo_fof '+paramFile,$
                str(cmdCode),$
                str(snap),$
                '>> run_fof.txt']
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
  workingPath = '/n/home07/dnelson/dev.tracer/'
  basePath    = workingPath + 'gasSphere.gasonly.2e5.cooling/'
  snapPath    = workingPath + 'gasSphere.gasonly.2e5.cooling/output/'

  ;redshift = 3.0
  ;snap     = redshiftToSnapNum(redshift,sP=sP)

  ;snapRange = [1,1,1]
  snapRange = [0,11,1]

  ; config - viewsize / object
  h = loadSnapshotHeader(basePath+'output/',snapNum=snapRange[0])
  
  ;xyzCen = [1123.20,7568.80,16144.2] ;dusan 512
  xyzCen = [h.boxSize/2.0,h.boxSize/2.0,h.boxSize/2.0]

  sliceWidth  = h.boxSize ; cube sidelength

  ; render config
  spawnJobs   = 1    ; execute bsub?
  nProcs      = 8    ; 128^3=8, 256^3=24, 512^3=96 (minimum for memory)
  dimX        = 800  ; image dimensions (x pixels)
  dimY        = 800  ; image dimensions (y pixels)
  
  axesStr = ['0 1 2','0 2 1','1 2 0'] ;xy,xz,yz
  
  ; bbox and projection setup
  cmdCode = 5 ;projection
 
  paramFile = "param.txt"
  
  xMin = xyzCen[0] - sliceWidth/2.0
  xMax = xyzCen[0] + sliceWidth/2.0
  yMin = xyzCen[1] - sliceWidth/2.0
  yMax = xyzCen[1] + sliceWidth/2.0
  zMin = xyzCen[2] - sliceWidth/2.0
  zMax = xyzCen[2] + sliceWidth/2.0  
  
  ; integration range for each axis depends on
  axesBB = [[xMin,xMax,yMin,yMax,zMin,zMax],$
            [xMin,xMax,zMin,zmax,yMin,yMax],$
            [yMin,yMax,zMin,zMax,xMin,xMax]]

 
  for snap=snapRange[0],snapRange[1],snapRange[2] do begin
 
    ; write bjob file
    jobFileName = workingPath+'job_proj.bsub'
    
    ; check before overriding
    if file_test(jobFileName) then begin
      print,'Error: job_proj.bsub already exists'
      return
    endif else begin
      print,'['+str(snap)+'] Writing: '+jobFilename
    endelse
      
    openW, lun, jobFilename, /GET_LUN
    
    ; write header
    printf,lun,'#!/bin/sh'
    printf,lun,'#BSUB -q keck'
    printf,lun,'#BSUB -J proj_' + str(snap)
    printf,lun,'#BSUB -n ' + str(nProcs)
    printf,lun,'#BSUB -R "rusage[mem=30000] span[ptile=8]"'
    printf,lun,'#BSUB -o run_proj.out'
    printf,lun,'#BSUB -g /dnelson/proj' ; 4 concurrent jobs limit automatically
    printf,lun,'#BSUB -cwd ' + basePath
    printf,lun,'#BSUB -a openmpi'
    printf,lun,''
      
    ; write projection commands
    foreach axisStr, axesStr, i do begin 
    
      strArray = ['mpirun.lsf ./Arepo '+paramFile,$
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
      printf,lun,'mv ' + snapPath + 'proj_density_field_' + string(snap,format='(I3.3)') + ' ' + $
                 snapPath + outputFilename
      printf,lun,''
    
    endforeach
      
    ; close
    close,lun
    free_lun,lun
      
    ; add to queue if requested
    if (spawnJobs) then begin
      spawn, 'bsub < ' + workingPath + 'job_proj.bsub', result
      print,'  '+result
      
      ; file cleanup
      wait,0.5
      spawn, 'rm ' + workingPath + 'job_proj.bsub'
    endif
  
  endfor ;snapRange

end

; sphMapBox: run sph kernel density projection on whole box
;

pro sphMapBox, res=res, run=run, partType=partType

  sP = simParams(res=res,run=run)
  
  ; config
  redshift = 3.0 ;5.0
  snap     = redshiftToSnapNum(redshift,sP=sP)
  
  nPixels = [800,800] ;px

  zoomFac = 1    ; only in axes, not along projection direction
  nNGB    = 64   ; use CalcHSML for HSML with nNGB
  axes    = [0,1] ; x,y

  ; paths and render config
  h = loadSnapshotHeader(sP.simPath,snapNum=snap)
  
  boxSize = [h.boxSize,h.boxSize,h.boxSize]              ;kpc
  boxCen  = [h.boxSize/2.0,h.boxSize/2.0,h.boxSize/2.0]  ;kpc
  
  foreach k,axes do boxSize[k] /= zoomFac
  
  outFilename = 'sphmap.box_'+str(zoomFac)+'.nNGB='+str(nNGB)+'.snap='+str(snap)+$
                '.box.axis0='+str(axes[0])+'.axis1='+str(axes[0])+'.'+partType

  ; save/restore
  if (file_test(sP.derivPath + outFilename + '.sav')) then begin
    restore,sP.derivPath + outFilename + '.sav',/verbose
  endif else begin

    if (partType eq 'gas') then begin
      mass = loadSnapshotSubset(sP.simPath,snapNum=snap,partType='gas',field='mass',/verbose)
    endif
    
    if (partType eq 'tracer') then begin
       mass_gas = loadSnapshotSubset(sP.simPath,snapNum=snap,partType='gas',field='mass',/verbose)
       mass = replicate(total(mass_gas) / h.nPartTot[3], h.nPartTot[3])
    endif
    
    ; load positions from snapshot
    pos  = loadSnapshotSubset(sP.simPath,snapNum=snap,partType=partType,field='pos',/verbose) 
    
    hsml = calcHSML(pos,ndims=3,nNGB=nNGB,boxSize=boxSize[0])

    ; OR: load HSML from snapshot (only stored for gas)
    ;hsml = loadSnapshotSubset(sP.simPath,snapNum=snap,partType='gas',field='hsml',/verbose)
      
    colMassMap = calcSphMap(pos,hsml,mass,boxSize=boxSize,boxCen=boxCen,nPixels=nPixels,$
                            axes=axes,ndims=3)
              
    save,colMassMap,hsml,filename=sP.derivPath + outFilename + '.sav'
  endelse

  ; rescale
  ;maxVal = max(colMassMap)/2.0
  maxVal = 0.5
  minVal = maxVal / 1e4
  
  print,'min max val: ',minVal,maxVal
  
  w = where(colMassMap eq 0, count, complement=ww)
  if (count ne 0) then colMassMap[w] = min(colMassMap[ww])
  
  colMassMap = colMassMap > minVal < maxVal
  
  colMassMap = alog10(colMassMap)
  colMassMap = (colMassMap-min(colMassMap))*254.0 / (max(colMassMap)-min(colMassMap)) ;0-254
  ;h2d = filter_image(h2d,FWHM=[1.1,1.1],/ALL) ;gaussian kernel convolution
  colMassMap += 1.0 ;1-255

  ; plot
  xMinMax = [boxCen[0]-boxSize[0]/2.0,boxCen[0]+boxSize[0]/2.0]
  yMinMax = [boxCen[1]-boxSize[1]/2.0,boxCen[1]+boxSize[1]/2.0]
  
  start_PS, sP.plotPath + outFilename + '.eps'
    loadct, 4, bottom=1, /silent
    tvim,colMassMap,xrange=xMinMax,yrange=yMinMax,pcharsize=0.0001
  end_PS, pngResize=68, /deletePS
  
end

; sphMapSubhalos: run sph kernel density projection on boxes centered on halos/subhalos

pro sphMapHalos, res=res, run=run, sgIDs=sgIDs

  if (not keyword_set(res)) then begin
    print,'Error: sphMapHalos: Bad inputs.'
    return
  endif

  sP = simParams(res=res,run=run)

  ; config
  boxSize = [200,200,200] ;kpc
  imgSize = [800,800]     ;px
  
  axes = list([0,1]) ;x,y
  mode = 1 ;1=col mass, 2=mass-weighted quantity, 3=col density
  
  redshift = 3.0
  snap     = redshiftToSnapNum(redshift,sP=sP)
  
  ; target list
  sgTarget = loadSubhaloGroups(sP.simPath,snap,/verbose)
    
  if (not keyword_set(sgIDs)) then begin
    ; make list of primary halo IDs (map them all)
    valSGids = sgIDList(sg=sgTarget,select='pri')
  endif else begin
    ; make specified IDs
    valSGids = sgIDs
  endelse

  ; load properties from snapshot
  pos  = loadSnapshotSubset(sP.simPath,snapNum=snap,partType='gas',field='pos',/verbose)
  hsml = loadSnapshotSubset(sP.simPath,snapNum=snap,partType='gas',field='hsml',/verbose)
  mass = loadSnapshotSubset(sP.simPath,snapNum=snap,partType='gas',field='mass',/verbose)
  
  ; loop over all non-background subhalos and image
  foreach sgID, valSGids do begin
  
    foreach axisPair, axes do begin
  
      imgFilename = 'sphmap.G.z='+str(redshift)+'.sgID='+str(sgID)+'.axis0='+$
                    str(axisPair[0])+'.axis1='+str(axisPair[1])+'.res='+str(res)+'.box='+str(boxSize[0])
                    
      if (file_test(sP.plotPath+imgFilename+'.png') or file_test(sP.plotPath+imgFilename+'.eps')) then begin
        print,'Skipping: ' + imgFilename
        continue
      endif     
  
      ; get subhalo position
      boxCen = sgTarget.subgroupPos[*,sgID]
      
      print,'['+string(sgID,format='(I04)')+'] Mapping ['+str(axisPair[0])+' '+$
            str(axisPair[1])+'] with '+str(boxSize[0])+$
            ' kpc box around subhalo center ['+str(boxCen[0])+' '+str(boxCen[1])+' '+str(boxCen[2])+']'
    
      colMassMap = sphDensityProjection(pos, hsml, mass, imgSize=imgSize, boxSize=boxSize,$
                                        boxCen=boxCen, axis0=axisPair[0], axis1=axisPair[1], $
                                        mode=mode, periodic=0, /verbose)
    
      ; rescale
      w = where(colMassMap eq 0, count, complement=ww)
      if (count ne 0) then colMassMap[w] = min(colMassMap[ww])
      
      colMassMap = alog10(colMassMap)
      colMassMap = (colMassMap-min(colMassMap))*254.0 / (max(colMassMap)-min(colMassMap)) ;0-254
      ;h2d = filter_image(h2d,FWHM=[1.1,1.1],/ALL) ;gaussian kernel convolution
      colMassMap += 1.0 ;1-254
    
      ; plot PS and PNG
      xMinMax = [boxCen[0]-boxSize[0]/2.0,boxCen[0]+boxSize[0]/2.0]
      yMinMax = [boxCen[1]-boxSize[1]/2.0,boxCen[1]+boxSize[1]/2.0]
      
      start_PS, sP.plotPath+imgFilename+'.eps'
        loadct, 4, bottom=1, /silent
        tvim,colMassMap,xrange=xMinMax,yrange=yMinMax,pcharsize=0.0001
        fsc_text,0.72,0.05,"z=3 id="+string(sgID,format='(I04)'),alignment=0.5,$
                 color=fsc_color('white'),/normal
      end_PS, pngResize=68, /deletePS
    
    endforeach ;axisPair
  
  endforeach ;valSGids
  
end
