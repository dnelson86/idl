; cosmoVis.pro
; cosmological boxes - 2d visualization
; dnelson feb.2012

; makeArepoFoFBsub(): write bsub file to invoke FoF/Subfind group finder post-processing

pro makeArepoFoFBsub

  ; config
  res = 256
  run = 'gadget'
  ;f = '10'
  
  sP = simParams(res=res,run=run)

  snapRange = [60,99,1]

  ; job config
  spawnJobs = 1 ; execute bsub?
  nProcs    = 32 ; needed nodes: 128^3=0.5 (n4tile4), 256^3=4 (n32tile4)
  ptile     = 4 ; span[ptile=X]
  cmdCode   = 3 ; fof/substructure post process
  
  excludeHosts = ['hero2701','hero1008','hero2405','hero1603'] ;leave empty otherwise
 
  paramFile = "param_fof.txt"

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
    printf,lun,'#BSUB -R "' + selectStr + 'rusage[mem=31000] span[ptile=' + str(ptile) + ']"'
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
  res = 256
  run = 'tracerMC.ref'
  f   = '1'
  subBox = 1 ; search for subbox snapshot versus normal
  
  snapRange = [1,847,1]
  
  ; config - viewsize / object
  sP = simParams(res=res,run=run,f=f,snap=snapRange[0])
  h = loadSnapshotHeader(sP=sP,subBox=subBox)
  
  ;boxSize = h.boxSize ; full box
  boxSize = 2000.0 ; subbox
  
  ;xyzCen = [1123.20,7568.80,16144.2] ;dusan 512
  xyzCen = [boxSize/2.0,boxSize/2.0,boxSize/2.0]

  sliceWidth  = boxSize ; cube sidelength
  zoomFac     = 1.1   ; final image is boxSize/zoomFac in extent

  ; render config
  spawnJobs   = 1    ; execute bsub?
  nProcs      = 1    ; 128^3=8, 256^3=24, 512^3=256 (minimum for memory)
  ptile       = 1
  dimX        = 800  ; image dimensions (x pixels)
  dimY        = 800  ; image dimensions (y pixels)
  
  axesStr = ['0 1 2','0 2 1'] ;['0 1 2','0 2 1','1 2 0'] ;xy,xz,yz
  
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
    endif else begin
      print,'['+str(snap)+'] Writing: '+jobFilename
    endelse
      
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

; sphMapBox: run sph kernel density projection on whole box

pro sphMapBox, res=res, run=run, partType=partType

  ; config
  ;res = 128
  ;run = 'dev.tracer.nocomov'
  ;partType = 'gas'

  ;redshift = 3.0 ;5.0
  ;snap     = redshiftToSnapNum(redshift,sP=sP)
  snap = 19
  
  nPixels = [800,800] ;px

  zoomFac = 1    ; only in axes, not along projection direction
  nNGB    = 64   ; use CalcHSML for HSML with nNGB
  axes    = [0,1] ; x,y

  ; paths and render config
  sP = simParams(res=res,run=run)
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

; sphDensityProjection(): (OLD) make density projection using SPH kernel (inspired by Mark's sphMap)
;                         NOTE: kernel coeffs only valid for 3D!

function sphDensityProjection, pos, hsml, mass, quantity=quantity, imgSize=imgSize, boxSize=boxSize,$
                               boxCen=boxCen, axis0=axis0, axis1=axis1, mode=mode, periodic=periodic,$
                               verbose=verbose

  print,'You should switch this to the calcSphMap C-routine.'
  stop

  ; config
  if not keyword_set(axis0) then axis0 = 0
  if not keyword_set(axis1) then axis1 = 1
  if not keyword_set(verbose) then verbose = 0
  
  if keyword_set(periodic) then begin
    print,'ERROR: PERIODIC not supported.'
    return,0
  endif
  
  if (mode ne 1 and mode ne 2 and mode ne 3) then begin
    print,'ERROR: Unsupported mode='+str(mode)+' parameter.'
    return,0
  endif
  
  ; storage
  p    = dblarr(3)
  pos0 = double(0.0)
  pos1 = double(0.0)
  binnedParticles = 0UL
  
  ; init
  npart = n_elements(hsml)

  grid = fltarr(imgSize[0],imgSize[1])
  
  if keyword_set(quantity) then $
    gridQuantity = fltarr(imgSize[0],imgSize[1])
  
  pxSize = [float(boxSize[0]) / imgSize[0], float(boxSize[1]) / imgSize[1]]
  pxArea = pxSize[0] * pxSize[1]

  if (pxSize[0] lt pxSize[1]) then $
    hMin = 1.001 * pxSize[0] / 2.0
  if (pxSize[0] ge pxSize[1]) then $
    hMin = 1.001 * pxSize[1] / 2.0
    
  hMax = pxSize[0] * 50.0
  
  for part=0, npart-1, 1 do begin
    ; progress report
    if (part mod round(npart/10.0) eq 0 and verbose) then $
      print,'Progress: '+string(100.0*part/npart,format='(I3)')+'%'
      
    ; get particle data
    p[0] = pos[0,part]
    p[1] = pos[1,part]
    p[2] = pos[2,part]
    h    = double(hsml[part])
    v    = double(mass[part])
    
    if keyword_set(quantity) then $
      w    = double(quantity[part])
    
    ; early exit if out of z-bounds
    if (abs(p[3-axis0-axis1] - boxCen[2]) gt boxSize[2] / 2.0) then $
      continue
      
    pos0 = p[axis0] - (boxCen[0] - boxSize[0] / 2.0)
    pos1 = p[axis1] - (boxCen[1] - boxSize[1] / 2.0)
    
    ; clamp hsml
    if (h lt hMin) then h = hMin;
    if (h gt hMax) then h = hMax;
    
    ; early exit if ...
    if (pos0 - 0.0 lt -h or pos1 - 0.0 lt -h or pos0 - boxSize[0] gt h or pos1 - boxSize[1] gt h) then $
      continue
      
    binnedParticles += 1
    
    h2 = h * h;
    
    ; number of pixels covered by particle
    nx = h / pxSize[0] + 1;
    ny = h / pxSize[1] + 1;
    
    ; coordinates of pixel center of particle
    x = (floor(pos0 / pxSize[0]) + 0.5) * pxSize[0]
    y = (floor(pos1 / pxSize[1]) + 0.5) * pxSize[1]
    
    ; normalization constant
    sum = 0.0
    
    for dx = -nx, nx, 1 do begin
      for dy = -ny, ny, 1 do begin
        ; dist of covered pixel from actual position
        xx = x + dx * pxSize[0] - pos0
        yy = y + dy * pxSize[1] - pos1
        r2 = xx*xx + yy*yy
        
        if (r2 < h2) then begin
          ; sph kernel (inlined): sum += _getkernel(h,r2);
          hinv = double(1.0) / h
          u    = sqrt(r2) * hinv
          
          if (u lt 0.5) then begin
            sum += (2.546479089470 + 15.278874536822 * (u - 1.0) * u * u)
          endif else begin
            sum += (5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u))
          endelse
        endif ;r2 < h2
      endfor
    endfor
    
    ; exit if negligible
    if (sum lt 1.0e-10) then $
      continue
      
    ; add contribution to image
    for dx = -nx, nx, 1 do begin
      for dy = -ny, ny, 1 do begin
        ; coordinates of pixel center of covering pixels
        xxx = x + dx * pxSize[0]
        yyy = y + dy * pxSize[1]
        
        ; pixel array indices
        i = floor(xxx / pxSize[0]) ;implicit C cast to int
        j = floor(yyy / pxSize[1]) ;same
        
        if (i ge 0 and i lt imgSize[0] and j ge 0 and j lt imgSize[1]) then begin
          xx = x + dx * pxSize[0] - pos0
          yy = y + dy * pxSize[1] - pos1
          r2 = xx*xx + yy*yy
          
          if (r2 lt h2) then begin
            ; divide by sum for normalization
            ; divide by pixelarea to get column density (optional: /pxArea)
            ; sph kernel (inlined): grid[] += _getkernel(h,r2) * v / sum
            hinv = double(1.0) / h
            u    = sqrt(r2) * hinv
            
            if (u lt 0.5) then begin
              grid[i * imgSize[1] + j] += $
                (2.546479089470 + 15.278874536822 * (u - 1.0) * u * u) * v / sum
              if keyword_set(quantity) then $
                gridQuantity[i * imgSize[1] + j] += $
                  (2.546479089470 + 15.278874536822 * (u - 1.0) * u * u) * v * w / sum
            endif else begin
              grid[i * imgSize[1] + j] += $
                (5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u)) * v / sum
                  if keyword_set(quantity) then $
                  gridQuantity[i * imgSize[1] + j] += $
                  (5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u)) * v * w / sum
            endelse
          
          endif ;r2 < h2
        endif ;i,j
      
      endfor
    endfor

  endfor ;part
  
  if (verbose) then print,'Number of binned particles: ',binnedParticles
  
  if (mode eq 1) then begin
    if (verbose) then print,'Returning: Column Mass Map'
    return,grid
  endif
  if (mode eq 2) then begin
    if (verbose) then print,'Returning: Quantity Mass-Weighted Map'
    return,gridQuantity
  endif
  if (mode eq 3) then begin
    if (verbose) then print,'Returning: Column Density Map'
    for i=0,i lt imgSize[0] do begin
      for j=0,j lt imgSize[1] do begin
        grid[i + imgSize[1] * j] /= pxArea
      endfor
    endfor
    
    return,grid
  endif

end

; sphMapSubhalos: run sph kernel density projection on boxes centered on halos/subhalos

pro sphMapHalos, sP=sP, gcIDs=gcIDs

  ; config
  boxSize = [200,200,200] ;kpc
  imgSize = [800,800]     ;px
  
  axes = list([0,1]) ;x,y
  mode = 1 ;1=col mass, 2=mass-weighted quantity, 3=col density
  
  ; target list
  gcTarget = loadGroupCat(sP=sP,/verbose)
    
  if (not keyword_set(gcIDs)) then begin
    ; make list of primary halo IDs (map them all)
    valGCids = gcIDList(gc=gcTarget,select='pri')
  endif else begin
    ; make specified IDs
    valGCids = gcIDs
  endelse

  ; load properties from snapshot
  pos  = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  hsml = loadSnapshotSubset(sP=sP,partType='gas',field='hsml')
  mass = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
  
  ; loop over all non-background subhalos and image
  foreach gcID, valGCids do begin
  
    foreach axisPair, axes do begin
  
      imgFilename = 'sphmap.G.snap='+str(sP.snap)+'.gcID='+str(gcID)+'.axis0='+$
                    str(axisPair[0])+'.axis1='+str(axisPair[1])+'.res='+str(sP.res)+'.box='+str(boxSize[0])
                    
      if (file_test(sP.plotPath+imgFilename+'.png') or file_test(sP.plotPath+imgFilename+'.eps')) then begin
        print,'Skipping: ' + imgFilename
        continue
      endif     
  
      ; get subhalo position
      boxCen = gcTarget.subgroupPos[*,gcID]
      
      print,'['+string(gcID,format='(I04)')+'] Mapping ['+str(axisPair[0])+' '+$
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
        fsc_text,0.72,0.05,"z=3 id="+string(gcID,format='(I04)'),alignment=0.5,$
                 color=fsc_color('white'),/normal
      end_PS, pngResize=68, /deletePS
    
    endforeach ;axisPair
  
  endforeach ;valGCids
  
end
