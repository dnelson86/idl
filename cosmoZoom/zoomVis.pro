; zoomVis.pro
; 'zoom project' specialized visualization
; dnelson jan.2015

; zoomMeshAndSliceWithRes(): Voronoi mesh and temperature slice @ L9,L10,L11 for one halo
; NOTE: requires that the snapshots be first post-processed with Arepo (restart-flag=4)
;       with the correct box center and size
; h7L9:  ./Arepo_slice param.txt 4 59 2000 2000 0 1 2 9743.99 10243.99 10269.0 10769.0 10302.2 #d
; h7L10: ./Arepo_slice param.txt 4 59 2000 2000 0 1 2 9736.63 10236.63 10194.3 10694.3 10286.5 #d
; h7L11: ./Arepo_slice param.txt 4 59 2000 2000 0 1 2 9671.24 10171.24 10206.4 10706.4 10220.9 #d
; h7L9:  ./Arepo_slice param.txt 5 59 4000 4000 0 1 2 9743.99 10243.99 10269.0 10769.0 10052.2 10552.2 #d 
; h7L10: ./Arepo_slice param.txt 5 59 4000 4000 0 1 2 9736.63 10236.63 10194.3 10694.3 10036.5 10536.5 #d
; h7L11: ./Arepo_slice param.txt 5 59 4000 4000 0 1 2 9671.24 10171.24 10206.4 10706.4 9970.90 10470.9 #d

; h0L9:  ./Arepo_slice param.txt 4 59 4000 4000 0 1 2 9179.4  9679.4  9992.7 10492.7  9903.73 #d
; h0L10: ./Arepo_slice param.txt 4 59 4000 4000 0 1 2 9192.7  9692.7  9924.0 10424.0  9860.74 #d
; h0L11: ./Arepo_slice param.txt 4 58 4000 4000 0 1 2 9249.0  9749.0 10011.6 10511.6  9775.30 #d(58)
; h0L9:  ./Arepo_slice param.txt 5 59 4000 4000 0 1 2 9179.4  9679.4  9992.7 10492.7  9653.7 10153.7 #d
; h0L10: ./Arepo_slice param.txt 5 59 4000 4000 0 1 2 9192.7  9692.7  9924.0 10424.0  9610.7 10110.7 #d
; h0L11: ./Arepo_slice param.txt 5 58 4000 4000 0 1 2 9249.0  9749.0 10011.6 10511.6  9525.3 10025.3 #d(58)

; h1L9:  ./Arepo_slice param.txt 4 59 4000 4000 0 1 2 10050.3 10550.3  9529.2 10029.2 10350.4 #d
; h1L10: ./Arepo_slice param.txt 4 59 4000 4000 0 1 2 10051.1 10551.1  9670.4 10170.4 10367.0 #d
; h1L11: ./Arepo_slice param.txt 4 59 4000 4000 0 1 2 10244.3 10744.3  9945.7 10445.7 10500.1 #d
; h1L9:  ./Arepo_slice param.txt 5 59 4000 4000 0 1 2 10050.3 10550.3  9529.2 10029.2 10100.4 10600.4 #d
; h1L10: ./Arepo_slice param.txt 5 59 4000 4000 0 1 2 10051.1 10551.1  9670.4 10170.4 10117.0 10617.0 #d
; h1L11: ./Arepo_slice param.txt 5 59 4000 4000 0 1 2 10244.3 10744.3  9945.7 10445.7 10250.1 10750.1 #d

; h8L9:  ./Arepo_slice param.txt 4 59 4000 4000 0 1 2 9501.5 10001.5 10009.8 10509.8 10629.6
; h8L10: ./Arepo_slice param.txt 4 59 4000 4000 0 1 2 9356.9  9856.9 10001.1 10501.1 10622.8
; h8L11: ./Arepo_slice param.txt 4 59 4000 4000 0 1 2 9222.1  9722.1  9849.1 10349.1 10773.8
; h8L9:  ./Arepo_slice param.txt 5 59 4000 4000 0 1 2 9501.5 10001.5 10009.8 10509.8 10379.6 10879.6
; h8L10: ./Arepo_slice param.txt 5 59 4000 4000 0 1 2 9356.9  9856.9 10001.1 10501.1 10372.8 10872.8
; h8L11: ./Arepo_slice param.txt 5 59 4000 4000 0 1 2 9222.1  9722.1  9849.1 10349.1 10523.8 11023.8

pro zoomMeshAndSliceWithRes
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  resLevels = [9,10,11]
  redshift  = 2.0
  hInd      = 0
  boxSize   = 500.0 ;ckpc side-length, but cannot change after post-processing
  fields    = ['density','metal','temp','velocity']
  rVirCircs = [0.15,0.5,1.0] ; times rvir to draw a circle
  boxBounds = [0.0,0.5] ; times boxSize, in both x & y, upper right quadrant
  
  ; load
  foreach resLevel,resLevels do begin
    sP = simParams(run='zoom_20Mpc',res=resLevel,hInd=hInd,redshift=redshift)

    sLoc = {}
    sLoc = mod_struct( sLoc, 'sP', sP )

    ; locate halo
    gc = loadGroupCat(sP=sP,/skipIDs)
    gcInd = zoomTargetHalo(sP=sP, gc=gc)
    gcPos = gc.subgroupPos[*,gcInd]
    
    sLoc = mod_struct( sLoc, 'rVir', gc.group_r_crit200[gcInd] )
    
    print,'L'+str(resLevel)+':', gcPos
    print,'L'+str(resLevel)+':', sP.targetHaloPos + sP.zoomShiftPhys
    print,string(gcPos[0]-0.5*boxSize,format='(f7.1)')+' '+$
          string(gcPos[0]+0.5*boxSize,format='(f7.1)')+' '+$
          string(gcPos[1]-0.5*boxSize,format='(f7.1)')+' '+$
          string(gcPos[1]+0.5*boxSize,format='(f7.1)')+' '+$
          string(gcPos[2]-0.5*boxSize,format='(f7.1)')+' '+$
          string(gcPos[2]+0.5*boxSize,format='(f7.1)')+' '

    ;if resLevel ge 11 then continue
    
    ; faces_list (Voronoi mesh slice): we pre-draw the mesh slice and save it as a raster image
    saveFilename = sP.derivPath + 'mesh_slice_image_' + str(sP.snap) + '_b' + str(fix(boxSize))
    
    if ~file_test(saveFilename+'.png') then begin
      print,'L'+str(resLevel)+': loading mesh...'
      
      flPath = sP.simPath + "faces_list_" + string(sP.snap,format='(I03)') + ".txt"
      readcol,flPath,x0,y0,x1,y1,format='F,F,F,F',/silent
      
      ; shift face vertices to be centered on halo, and restrict to reasonable box
      x0 -= gcPos[0]
      x1 -= gcPos[0]
      y0 -= gcPos[1]
      y1 -= gcPos[1]
      
      w = where((abs(x0) lt boxBounds[1]*boxSize and abs(y0) lt boxBounds[1]*boxSize) or $
                (abs(x1) lt boxBounds[1]*boxSize and abs(y1) lt boxBounds[1]*boxSize),count)
      
      print,'L'+str(resLevel)+': Mesh faces keeping ['+str(count)+'] of ['+str(n_elements(x0))+'].'
      
      ; start image
      start_PS, saveFilename + '.eps', xs=20, ys=20
        cgPlot,[0],[0],/nodata,xrange=[-0.5,0.5]*boxSize,yrange=[-0.5,0.5]*boxSize,xs=5,ys=5,$
          pos=[0,0,1,1]
          
        for k=0,count-1 do begin
          xx = [x0[w[k]],x1[w[k]]]
          yy = [y0[w[k]],y1[w[k]]]
          
          ; plot line segment (automatic clipping to plot box)
          cgPlot,xx,yy,color='black',/overplot
        endfor
        
      end_PS, pngResize=100
      
    endif
    
    ; quantity slices
    foreach field,fields do begin
      rr = loadSliceField(sP.simPath, sP.snap, field)
      sLoc = mod_struct( sLoc, field, rr.field )
    endforeach
    
    ; quantity projections
    rr = loadDensityField(sP.simPath, sP.snap)
    sLoc = mod_struct( sLoc, 'proj_density', rr.dens )
    sLoc = mod_struct( sLoc, 'proj_temp', rr.temp )    
    
    slices = mod_struct( slices, 'res_'+str(resLevel), sLoc )
    
  endforeach
  
  ; positioning
  x0 = 0.05 & x1 = 0.33 & xoff = 0.32
  y0 = 0.08 & y1 = 0.47 & yoff = 0.46
  
  pos = list( [x0+0*xoff, y0+1*yoff, x1+0*xoff, y1+1*yoff] ,$ ; top left
              [x0+1*xoff, y0+1*yoff, x1+1*xoff, y1+1*yoff] ,$ ; top middle
              [x0+2*xoff, y0+1*yoff, x1+2*xoff, y1+1*yoff] ,$ ; top right
              [x0+0*xoff, y0+0*yoff, x1+0*xoff, y1+0*yoff] ,$ ; bottom left
              [x0+1*xoff, y0+0*yoff, x1+1*xoff, y1+0*yoff] ,$ ; bottom middle
              [x0+2*xoff, y0+0*yoff, x1+2*xoff, y1+0*yoff] ) ; bottom right
  
  ; plot
  start_PS, slices.(0).sP.plotPath + 'zoomMeshAndSlices_h'+str(hInd)+'.eps', xs=12.0, ys=8.0

    foreach resLevel,resLevels,i do begin
      j = i < 2
      
      ; ckpc/h -> physical kpc
      scalefac = 1.0/(1+redshift)
      convFac = scalefac / units.HubbleParam
      slices.(j).rVir *= convFac
      xyrange = [0,120] ;xyrange = convFac*boxBounds*boxSize
      
      ; top row: mesh
      posTop = pos[i+0*n_elements(resLevels)]
      cgPlot,[0],[0],/nodata,xrange=xyrange,yrange=xyrange,xs=5,ys=5,$
        pos=posTop,/noerase,aspect=1.0
      
      loadColorTable,'bw linear'
      
      ; load raster image of mesh
      saveFilename = slices.(j).sP.derivPath + 'mesh_slice_image_' + $
                     str(slices.(j).sP.snap) + '_b' + str(fix(boxSize))
      meshRaster = read_png( saveFilename + '.png' )
      meshRaster = reform( meshRaster[0,*,*] )
      
      ; quadrant
      sz = n_elements(meshRaster[0,*])
      meshRaster = meshRaster[sz/2:sz-1,sz/2:sz-1]
      
      tv, meshRaster, posTop[0], posTop[1], /normal, xsize=posTop[2]-posTop[0]
            
      ; circle at virial radius fractions
      foreach rVirCirc,rVirCircs do $
        tvcircle,rVirCirc * slices.(j).rVir,0,0,cgColor('red'),thick=0.8,/data,/quadrant
      
      cgPlot,[0],[0],/nodata,xrange=xyrange,yrange=xyrange,/xs,/ys,$
        pos=posTop,title='h'+str(hInd)+'L'+str(resLevel),/noerase,aspect=1.0
      
      ; bottom row: slice quantity
      posBottom = pos[i+1*n_elements(resLevels)]
      cgPlot,[0],[0],/nodata,xrange=xyrange,yrange=xyrange,xs=5,ys=5,$
        pos=posBottom,/noerase,aspect=1.0
        
      ; colormap (dens slice)
      ;loadColorTable,'helix'
      ;mapMinMax = [-3.0,5.0]
      ;sliceMap = rhoRatioToCrit( slices.(j).density > 1e-20, sP=sP, /log )
      
      ; colormap (temp slice)
      ;loadColorTable,'blue-red2'
      ;mapMinMax = [4.0,5.0]
      ;sliceMap = alog10( slices.(j).temp > 1e-20 )
      
      ; colormap (dens proj)
      ;mapMinMax = [2.0,5.0]
      ;sliceMap = rhoRatioToCrit( slices.(j).proj_density > 1e-20, sP=sP, /log )
      
      ; colormap (temp proj)
      loadColorTable,'blue-red2'
      mapMinMax = [3.4,4.6]
      sliceMap = alog10( slices.(j).proj_temp > 1e-20 )
      
      ; stretch
      sliceMap = (sliceMap-mapMinMax[0])*(255.0) / (mapMinMax[1]-mapMinMax[0])
      sliceMap = fix(sliceMap) > 0 < 255 ; 0-255 
      
      ; quadrant
      sz = n_elements(sliceMap[0,*])
      sliceMap = sliceMap[sz/2:sz-1,sz/2:sz-1]
      
      tv, sliceMap, posBottom[0], posBottom[1], /normal, xsize=posBottom[2]-posBottom[0]

      ; circle at virial radius fractions
      foreach rVirCirc,rVirCircs do $
        tvcircle,rVirCirc * slices.(j).rVir,0,0,cgColor('white'),thick=0.8,/data,/quadrant
      
      cgPlot,[0],[0],/nodata,xrange=xyrange,yrange=xyrange,/xs,/ys,$
        pos=posBottom,/noerase,aspect=1.0
      
    endforeach
  
  end_PS
  
end

; multiMapZoomComp(): compare 3 or 4 physical quantities (using sphMap) of four zooms

pro multiMapZoomComp, doOne=doOne, doTwo=doTwo
  compile_opt idl2, hidden, strictarr, strictarrsubs
  if n_elements(doOne) eq 0 and n_elements(doTwo) eq 0 then message,'Error'
  
  ; config
  redshift = 2.0
  resLevel = [11,11,11,11] ;[9,10,10,10] ;
  hInds    = [2,4,5,9] ;[0,1,7,8] ;[6,3] ;[3,3,5,9] ;

  ; plot config
  sizeFac       = 7.8          ; times rvir (master, must be >=max(sizeFac) below)
  hsmlFac       = 3.0 ; 2.0 in draft0! ; times each cell radius for sph projections
  nPixels       = [1200,1200]  ; px
  xySize        = 3            ; final image is xySize*nPixels[0] high
  axisPair      = [0,1]        ; xy, xz
  barAreaHeight = 0.06         ; fractional
  
  ; use which field and minmax for color mappings?
  if keyword_set(doOne) then begin
    ; gas dens,temp,entr,vrad on the virial scale
    saveTag = 'a'
    
    fields = { $
      field4 : { cF:'overdens', mapMM:[-1.0,4.0], sizeFac:3.0, rVirCircs:[0.15,0.5,1.0] } ,$
      field0 : { cF:'temp',     mapMM:[4.4,6.3],  sizeFac:3.0, rVirCircs:[0.15,0.5,1.0] } ,$
      field1 : { cF:'entropy',  mapMM:[6.5,9.1],  sizeFac:3.0, rVirCircs:[0.15,0.5,1.0] } ,$
      field2 : { cF:'vrad',     mapMM:[-220,220], sizeFac:3.0, rVirCircs:[0.15,0.5,1.0] } }
  endif
  
  if keyword_set(doTwo) then begin
    ; gas dens,stellar dens zoomed in + gas dens,dm dens zoomed out
    saveTag = 'b'
    
    fields = { $
      field4 : { cF:'overdens', mapMM:[1.0,6.0],  sizeFac:1.0, rVirCircs:[0.05,0.15] } ,$
      field0 : { cF:'stardens', mapMM:[2.0,6.0],  sizeFac:1.0, rVirCircs:[0.05,0.15] } ,$
      field1 : { cF:'overdens', mapMM:[-1.0,4.0], sizeFac:5.5, rVirCircs:[1.0,2.0]   } ,$
      field2 : { cF:'dmdens',   mapMM:[1.5,5.5],  sizeFac:5.5, rVirCircs:[1.0,2.0]   } }
  endif
    
  numFields = n_tags(fields)
  numRuns   = n_elements(hInds)
  
  ; load runs
  gcIDs = []
  
  foreach hInd,hInds,i do begin
    sP = mod_struct(sP, 'sP'+str(i), $
                    simParams(run='zoom_20Mpc',res=resLevel[i],hInd=hInd,redshift=redshift))
    gcIDs = [ gcIDs, zoomTargetHalo(sP=sP.(i)) ]
  endforeach
  
  runsStr = 'L' + str(resLevel[0]) + '_' + str(sP.(0).snap) + '_h' + strjoin(str(hInds))
  
  saveFilename = 'zoomMaps_'+runsStr+'_axes'+str(axisPair[0])+str(axisPair[1])+'_'+$
                  saveTag + '_px' + str(nPixels[0]) + '_nF' +str(numFields)+'.eps'
                 
  xs = xySize*numFields
  ys = xySize*numRuns*(1/(1.0-barAreaHeight))
  
  start_PS, sP.(0).plotPath + saveFilename, xs=xs, ys=ys
  
  ; loop over runs
  for i=0,numRuns-1 do begin
  
    print,sP.(i).saveTag + ' ['+str(gcIDs[i])+'] Mapping axes ['+str(axisPair[0])+','+str(axisPair[1])+']'
  
    ; spatial cutouts, increase hsml for vis purposes
    mapCutout = cosmoVisCutout(sP=sP.(i),gcInd=gcIDs[i],sizeFac=sizeFac)
    mapCutout.loc_hsml *= hsmlFac
    
    for j=0,numFields-1 do begin
      ; make boxsize appropriate for map cutout
      mapCutout.boxSizeImg = fields.(j).sizeFac*[mapCutout.haloVirRad,$ ;x
                                                 mapCutout.haloVirRad*nPixels[1]/nPixels[0],$ ;y
                                                 mapCutout.haloVirRad] ;z
  
      ; mapping configuration
      config = {saveFilename:'',nPixels:nPixels,axes:axisPair,fieldMinMax:[0,0],$
                gcID:gcIDs[i],haloMass:mapCutout.haloMass,haloVirRad:mapCutout.haloVirRad,$
                boxCen:[0,0,0],boxSizeImg:mapCutout.boxSizeImg,rVirCircs:fields.(j).rVirCircs,$
                ctNameScat:'',ctNameMap:'',sP:sP.(i),bartype:'',scaleBarLen:100.0,$
                secondCutVal:-1,secondText:'',nbottom:0,secondGt:1,singleColorScale:1,$
                barAreaHeight:barAreaHeight,newBoxSize:fields.(j).sizeFac,$
                colorField:fields.(j).cF,mapMinMax:fields.(j).mapMM}
      if keyword_set(doTwo) then config = mod_struct( config, 'secondScaleBar', 1 )
            
      print,' '+fields.(j).cF
            
      ;tempFilename = sP.(0).plotPath + 'temp_'+saveTag+'_L'+str(resLevel[0])+'_'+str(i)+'_'+str(j)+'.sav'
      ;if file_test(tempFilename) then begin
      ;  restore,tempFilename
      ;endif else begin
        sub = cosmoVisCutoutSub(cutout=mapCutout,mapCutout=mapCutout,config=config)
      ;  save,sub,config,filename=tempFilename
      ;endelse
      
      ; plot
      plotMultiSphmap, map=sub, config=config, row=[i,numRuns], col=[j,numFields]
  
    endfor ; numFields,j
  endfor ; numRuns,i
  
  end_PS, density=ceil(nPixels[0]/xySize), pngResize=100 ;, /deletePS
  
end

; multiMapZoomRotFrames(): compare many quantities and scales of one zoom, and generate frames 
;   rotating around an axis (no time evolution) for a movie
; doSnap: if specified, track through time using tree and plot one frame (no rotation) at each snap
pro multiMapZoomRotFrames, jobNum=jobNum, totJobs=totJobs, do4x2=do4x2, do2x1=do2x1, doSnaps=doSnaps, hInd=hInd
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if n_elements(jobNum) eq 0 or n_elements(totJobs) eq 0 then message,'Error'
  if n_elements(do4x2) eq 0 and n_elements(do2x1) eq 0 then message,'Error'
  
  ; config
  redshift = 2.0
  resLevel = 11
  ;hInd     = 2
  ;if n_elements(hInd) eq 0 then message,'Error'
  outPath  = '/n/home07/dnelson/data5/frames/'

  hsmlFac      = 2.0      ; times each cell radius for sph projections
  axes         = [0,1,2]  ; projection in direction of 3rd axis, rotation about first
  framesPerRot = 600      ; 20 seconds @ 30fps
  rotAxis      = 'y'      ; x or y
  
  ; zoom4x2:
  if keyword_set(do4x2) then begin
    sizeFac       = 7.8           ; times rvir (master, must be sqrt(2)*max(sizeFac) below)
    nPixels       = [960,960]     ; sphMap projection grid size in pixels
    xySize        = [3840,2160]   ; final image size
    barAreaHeight = 0.108 ;14     ; fractional
    numRows       = 2
    fields        = { $
      $ ; top
      field0 : { colorField:'overdens', mapMinMax:[-0.5,4.0], sizeFac:3.0, rVirCircs:[0.15,0.5,1.0] } ,$
      field1 : { colorField:'temp',     mapMinMax:[4.4,6.3],  sizeFac:3.0, rVirCircs:[0.15,0.5,1.0] } ,$
      field2 : { colorField:'entropy',  mapMinMax:[7.5,8.8],  sizeFac:3.0, rVirCircs:[0.15,0.5,1.0] } ,$
      field3 : { colorField:'vrad',     mapMinMax:[-200,200], sizeFac:3.0, rVirCircs:[0.15,0.5,1.0] } ,$
      $ ; bottom
      field4 : { colorField:'overdens', mapMinMax:[1.0,6.0],  sizeFac:1.0, rVirCircs:[0.15] } ,$
      field5 : { colorField:'stardens', mapMinMax:[2.0,6.0],  sizeFac:1.0, rVirCircs:[0.15] } ,$
      field6 : { colorField:'overdens', mapMinMax:[-1.0,4.0], sizeFac:5.5, rVirCircs:[0.5,1.0] } ,$
      field7 : { colorField:'dmdens',   mapMinMax:[1.5,5.5],  sizeFac:5.5, rVirCircs:[0.5,1.0] } }
  endif ;0
  
  ; zoom2x1:
  if keyword_set(do2x1) then begin
    sizeFac       = 4.3           ; times rvir (master, must be sqrt(2)*max(sizeFac) below)
    nPixels       = [1920,2160]   ; sphMap projection grid size in pixels
    xySize        = [3840,2160]   ; final image size
    barAreaHeight = 0             ; fractional
    numRows       = 1
    fields        = { $
      field0 : { colorField:'overdens', mapMinMax:[-0.5,4.5], sizeFac:3.0, rVirCircs:[0.15,0.5,1.0] } ,$
      field1 : { colorField:'temp',     mapMinMax:[4.4,6.3],  sizeFac:3.0, rVirCircs:[0.15,0.5,1.0] }  }
  endif ;0

  ; decide frameMM based on jobNum
  if float(framesPerRot)/totJobs ne round(float(framesPerRot)/totJobs) then $
    message,'Error: framesPerRot not divisable by totJobs.'
    
  framesPerJob = framesPerRot/totJobs
  frameMM = [framesPerJob*jobNum,framesPerJob*(jobNum+1)-1]
  numFields = n_tags(fields)/numRows
  
  ; load
  sP = simParams(run='zoom_20Mpc',res=resLevel,hInd=hInd,redshift=redshift)
  gcID = [zoomTargetHalo(sP=sP)]
  
  ; do rotation at one fixed snapshot?
  if n_elements(doSnaps) eq 0 then begin
    snapNums = [sP.snap]
    print,'Job ['+str(jobNum)+'] of ['+str(totJobs)+'] rotation frames: ',frameMM
  endif
  
  ; do time evolution of tracked halo at one fixed rotation angle?
  if n_elements(doSnaps) ne 0 then begin
    ; time evolution, track using tree
    track = trackHaloPosition(sP=sP, gcID=gcID[0], endSnap=0)
    
    wSnaps = where(track.gcIDs ge 0,countSnaps)
    
    ; override sP.snap, target subgroup ID, frame bracketing
    snapsPerJob = countSnaps/totJobs
    snapInds = wSnaps[snapsPerJob*jobNum:snapsPerJob*(jobNum+1)-1]
    snapNums = track.snaps[snapInds]
    
    print,'Job ['+str(jobNum)+'] of ['+str(totJobs)+'] do snapshots: ',snapNums
    
    frameMM = [0,0]
    gcID = track.gcIDs[snapInds]
  endif
  
  if n_elements(gcID) ne n_elements(snapNums) then message,'Error: Mismatch.'

  foreach snapNum,snapNums,snapCount do begin
  
    sP.snap = snapNum
    sP.redshift = snapNumToRedshift(sP=sP)
  
    mapCutout = cosmoVisCutout(sP=sP,gcInd=gcID[snapCount],sizeFac=sizeFac)
    mapCutout.loc_hsml *= hsmlFac
    
    for frameNum=frameMM[0],frameMM[1] do begin
      rotAngle = (float(frameNum)/framesPerRot*(2*!pi)) ;mod (2*!pi)
      print,'[snap='+string(snapNum,format='(I02)')+'] ['+$
            string(frameNum,format='(I04)')+'] rotAngle = '+string(rotAngle,format='(f6.4)')

      saveFilename = 'zoomMaps_'+str(numFields)+'x'+str(numRows)+'_'+sP.saveTag+'_'+string(sP.snap,format='(I02)')+$
        '_rot-' + rotAxis + '_frame_'+string(frameNum,format='(I04)')+'.eps'
                     
      if file_test(outPath + strmid(saveFilename,0,strlen(saveFilename)-4)+'.png') then begin
        print,' SKIP'
        continue
      endif
                     
      fac = 1.25
      xs = xySize[0]/(300/fac)
      ys = xySize[1]/(300/fac)
      
      start_PS, outPath + saveFilename, xs=xs, ys=ys
      
      ; loop over rows
      for i=0,numRows-1 do begin
        for j=0,numFields-1 do begin
        
          field = fields.(i*numFields+j)
          
          ; make boxsize smaller for map cutout
          mapCutout.boxSizeImg = field.sizeFac*[mapCutout.haloVirRad,$ ;x
                                                mapCutout.haloVirRad*nPixels[1]/nPixels[0],$ ;y
                                                mapCutout.haloVirRad] ;z

          ; mapping configuration
          config = {saveFilename:'',nPixels:nPixels,axes:axes[0:1],fieldMinMax:[0,0],$
                    gcID:gcID[snapCount],haloMass:mapCutout.haloMass,haloVirRad:mapCutout.haloVirRad,$
                    boxCen:[0,0,0],boxSizeImg:mapCutout.boxSizeImg,rVirCircs:field.rVirCircs,$
                    ctNameScat:'',ctNameMap:'',sP:sP,bartype:'',scaleBarLen:200.0,$
                    secondCutVal:-1,secondText:'',nbottom:0,secondGt:1,singleColorScale:1,$
                    barAreaHeight:barAreaHeight,newBoxSize:field.sizeFac,rotAxis:rotAxis,rotAngle:rotAngle,$
                    colorField:field.colorField,mapMinMax:field.mapMinMax}
          
          sub = cosmoVisCutoutSub(cutout=mapCutout,mapCutout=mapCutout,config=config)
          
          ; plot
          plotMultiSphmap, map=sub, config=config, row=[i,numRows], col=[j,numFields]

        endfor ; numFields,j
      endfor ; numRows,i
      
      end_PS, density=300/fac, pngResize=100, /deletePS
      
    endfor ; frameNum
  endforeach ; snapNums
  
end
