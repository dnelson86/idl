; coldflowsVis.pro
; cold flows - 2d visualization
; dnelson oct.2011

; makeArepoProjBsub(): write bsub file to invoke six different axis aligned projections using the 
;                      Arepo voronoi_makeimage_new() code
;                      (output filenames have been made unique with change to Arepo code)

pro makeArepoProjBsub

  ; execute bsub?
  spawnJobs = 1
  paramFile = "paramCut.txt"

  ; object config
  ;xCen = 1123.20 & yCen = 7568.80 & zCen = 16144.2 ;dusan 512
  
  targetRedshift = 3.0
  targetSnap     = redshiftToSnapNum(targetRedshift)
  
  ; render config
  workingPath = '/n/home07/dnelson/coldflows/vis/'
  nProcs      = 8   ; 128^3=8, 256^3=24, 512^3=96
  dimX        = 500  ; image dimensions (x)
  dimY        = 500  ; image dimensions (y)
  
  sliceWidth  = 200.0 ; cube sidelength
  
  ; bbox and projection setup
  cmdCode = 5 ;projection
  
  xMin = xCen - sliceWidth/2.0
  xMax = xCen + sliceWidth/2.0
  yMin = yCen - sliceWidth/2.0
  yMax = yCen + sliceWidth/2.0
  zMin = zCen - sliceWidth/2.0
  zMax = zCen + sliceWidth/2.0  
  
  ; integration range for each axis depends on
  axesBB = [[xMin,xMax,yMin,yMax,zMin,zMax],$
            [xMin,xMax,zMin,zmax,yMin,yMax],$
            [yMin,yMax,zMin,zMax,xMin,xMax]]

  axesStr = ['0 1 2','0 2 1','1 2 0']

  ; write bjob files
  foreach axisStr, axesStr, i do begin  
  
    ; check before overriding
    if (file_test(workingPath+'job'+str(i)+'.bsub')) then begin
      print,'Error: job'+str(i)+'.bsub already exists'
      return
    endif
    
    ; create jobI.bsub
    openW, lun, workingPath+'job'+str(i)+'.bsub', /GET_LUN
    
    ; write header
    printf,lun,'#!/bin/sh'
    printf,lun,'#BSUB -q nancy'
    printf,lun,'#BSUB -J cf_vis_'+str(i)
    printf,lun,'#BSUB -n ' + str(nProcs)
    printf,lun,'#BSUB -R "rusage[mem=30000] span[ptile=8]"'
    printf,lun,'#BSUB -o run.out'
    printf,lun,'#BSUB -e run.err'
    printf,lun,'#BSUB -a openmpi'
    printf,lun,'#BSUB -N'
    printf,lun,''
    
    ; write projection commands
    strArray = ['mpirun.lsf ./Arepo '+paramFile,$
                str(cmdCode),$
                str(targetSnap),$
                str(dimX),str(dimY),$
                axisStr,$
                str(axesBB[0,i]),str(axesBB[1,i]),$
                str(axesBB[2,i]),str(axesBB[3,i]),$
                str(axesBB[4,i]),str(axesBB[5,i]),$
                '> run_'+str(i)+'.txt']
    printf,lun,strjoin(strArray,' ')
    
    ; close
    close,lun
    free_lun,lun
    
  endforeach
  
  ; add to queue if requested
  if (spawnJobs) then for i=0,n_elements(axesStr)-1 do spawn,'bsub < job'+str(i)+'.bsub'

end

; sphMapBox: run sph kernel density projection on whole box

pro sphMapBox, res=res

  if not keyword_set(res) then begin
    print,'Error: Must specific resolution set.'
    return
  endif
  
  ; config
  gadgetPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'
  workingPath  = '/n/home07/dnelson/coldflows/'
  
  boxSize = [10000,10000,10000] ;kpc
  boxCen  = [10000,10000,10000] ;kpc
  imgSize = [500,500]           ;px
  
  axis0 = 0 ;x
  axis1 = 1 ;y
  mode  = 1 ;1=col mass, 2=mass-weighted quantity, 3=col density
  
  targetRedshift = 3.0
  targetSnap     = redshiftToSnapNum(targetRedshift)
  
  imgFilename = 'sphmap.subhalo.z='+str(targetRedshift)+'.box.axis0='+str(axis0)+$
                '.axis1='+str(axis1)+'.res='+str(res)
                
  ; load properties from snapshot
  pos  = loadSnapshotSubset(gadgetPath,snapNum=targetSnap,partType='gas',field='pos',/verbose)
  hsml = loadSnapshotSubset(gadgetPath,snapNum=targetSnap,partType='gas',field='hsml',/verbose)
  mass = loadSnapshotSubset(gadgetPath,snapNum=targetSnap,partType='gas',field='mass',/verbose)

  colMassMap = sphDensityProjection(pos, hsml, mass, imgSize=imgSize, boxSize=boxSize,$
                                    boxCen=boxCen, axis0=axis0, axis1=axis1, mode=mode, periodic=0)

  ; rescale
  maxVal = max(colMassMap)/2.0
  minVal = maxVal / 1e4
  
  w = where(colMassMap eq 0, count, complement=ww)
  if (count ne 0) then colMassMap[w] = min(colMassMap[ww])
  
  w = where(colMassMap gt maxVal, count)
  if (count ne 0) then colMassMap[w] = maxVal
  w = where(colMassMap lt minVal, count)
  if (count ne 0) then colMassMap[w] = minVal
  
  colMassMap = alog10(colMassMap)
  colMassMap = (colMassMap-min(colMassMap))*254.0 / (max(colMassMap)-min(colMassMap)) ;0-250
  ;h2d = filter_image(h2d,FWHM=[1.1,1.1],/ALL) ;gaussian kernel convolution
  colMassMap += 1.0 ;1-251

  ; plot PS and PNG
  xMinMax = [boxCen[0]-boxSize[0]/2.0,boxCen[0]+boxSize[0]/2.0]
  yMinMax = [boxCen[1]-boxSize[1]/2.0,boxCen[1]+boxSize[1]/2.0]
  
  start_PS, workingPath+imgFilename+'.eps'
    loadct, 4, bottom=1, /silent
    tvim,colMassMap,xrange=xMinMax,yrange=yMinMax,pcharsize=0.0001
  end_PS, pngResize=68
  
end

; sphMapSubhalos: run sph kernel density projection on boxes centered on halos/subhalos

pro sphMapHalos, res=res, halos=halos, sgIDs=sgIDs

  if (not keyword_set(res) or (halos ne 0 and halos ne 1)) then begin
    print,'Error: sphMapHalos: Bad inputs.'
    return
  endif

  ; config
  gadgetPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'
  ;workingPath  = '/n/home07/dnelson/coldflows/vis/subhalo_imgs_'+str(res)+'G/'
  workingPath  = '/n/home07/dnelson/coldflows/vis/coldfil/'
  
  boxSize = [200,200,200] ;kpc
  imgSize = [500,500]     ;px
  
  axes = list([0,1]) ;x,y
  mode = 1 ;1=col mass, 2=mass-weighted quantity, 3=col density
  
  targetRedshift = 3.0
  targetSnap     = redshiftToSnapNum(targetRedshift)
  
  sgTarget = loadSubhaloGroups(gadgetPath,targetSnap,/verbose)
    
  if (not keyword_set(sgIDs)) then begin
    ; make list of halos/subhalo IDs (map them all)
    valSGids = getPrimarySubhaloList(sgTarget,halos=halos)
  endif else begin
    ; make specified IDs
    valSGids = sgIDs
  endelse

  ; load properties from snapshot
  pos  = loadSnapshotSubset(gadgetPath,snapNum=targetSnap,partType='gas',field='pos',/verbose)
  hsml = loadSnapshotSubset(gadgetPath,snapNum=targetSnap,partType='gas',field='hsml',/verbose)
  mass = loadSnapshotSubset(gadgetPath,snapNum=targetSnap,partType='gas',field='mass',/verbose)
  
  ; loop over all non-background subhalos and image
  foreach sgID, valSGids do begin
  
    foreach axisPair, axes do begin
  
      imgFilename = 'sphmap.G.z='+str(targetRedshift)+'.sgID='+str(sgID)+'.axis0='+$
                    str(axisPair[0])+'.axis1='+str(axisPair[1])+'.res='+str(res)+'.box='+str(boxSize[0])
                    
      if (file_test(workingPath+imgFilename+'.png') or file_test(workingPath+imgFilename+'.eps')) then begin
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
      
      start_PS, workingPath+imgFilename+'.eps'
        loadct, 4, bottom=1, /silent
        tvim,colMassMap,xrange=xMinMax,yrange=yMinMax,pcharsize=0.0001
        fsc_text,0.72,0.05,"z=3 id="+string(sgID,format='(I04)'),alignment=0.5,$
                 color=fsc_color('white'),/normal
      end_PS, pngResize=68, /deletePS
    
    endforeach ;axisPair
  
  endforeach ;valSGids
  
end