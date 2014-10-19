; zoomVis.pro
; 'zoom project' specialized visualization
; dnelson oct.2014

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
      
      ; top row: mesh
      posTop = pos[i+0*n_elements(resLevels)]
      cgPlot,[0],[0],/nodata,xrange=boxBounds*boxSize,yrange=boxBounds*boxSize,xs=5,ys=5,$
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
      
      cgPlot,[0],[0],/nodata,xrange=boxBounds*boxSize,yrange=boxBounds*boxSize,/xs,/ys,$
        pos=posTop,title='h'+str(hInd)+'L'+str(resLevel),/noerase,aspect=1.0
      
      ; bottom row: slice quantity
      posBottom = pos[i+1*n_elements(resLevels)]
      cgPlot,[0],[0],/nodata,xrange=boxBounds*boxSize,yrange=boxBounds*boxSize,xs=5,ys=5,$
        pos=posBottom,/noerase,aspect=1.0
      
      ;loadColorTable,'blue-red2'
      loadColorTable,'helix'
  
      ; colormap (dens slice)
      mapMinMax = [-3.0,5.0]
      sliceMap = rhoRatioToCrit( slices.(j).density > 1e-20, sP=sP, /log )
      
      ; colormap (temp slice)
      ;mapMinMax = [4.0,5.0]
      ;sliceMap = alog10( slices.(j).temp > 1e-20 )
      
      ; colormap (dens proj)
      ;mapMinMax = [2.0,5.0]
      ;sliceMap = rhoRatioToCrit( slices.(j).proj_density > 1e-20, sP=sP, /log )
      
      ; colormap (temp proj)
      ;mapMinMax = [3.4,4.6]
      ;sliceMap = alog10( slices.(j).proj_temp > 1e-20 )
      
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
      
      cgPlot,[0],[0],/nodata,xrange=boxBounds*boxSize,yrange=boxBounds*boxSize,/xs,/ys,$
        pos=posBottom,/noerase,aspect=1.0
      
    endforeach
  
  end_PS
  
end

; multiMapZoomComp(): compare 3 or 4 physical quantities (using sphMap) of four zooms

pro multiMapZoomComp

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  redshift = 2.02 ; NOTE!!!!!!!!!!!!!!!!!
  resLevel = 9
  hInds    = [2,4,5,6] ;[0,1,7,8]

  sizeFacMap  = 3.0            ; times rvir (image width, e.g. 2.0 shows out to the virial radius)
  rVirCircs   = [0.15,0.5,1.0] ; times rvir to draw a circle
  hsmlFac     = 2.0            ; times each cell radius for sph projections
  nPixels     = [1200,1200]    ; px
  xySize      = 3              ; final image is xySize*nPixels[0] high
  axisPair    = [0,1]          ; xy, xz
  barAreaHeight = 0.06         ; fractional
  
  ; use which field and minmax for color mappings?
  fields = { $
    field4 : { colorField : 'overdens',    mapMinMax : [-1.0,4.0] } ,$
    field0 : { colorField : 'temp',        mapMinMax : [4.4,6.2] } ,$
    field1 : { colorField : 'entropy',     mapMinMax : [6.5,8.9] } ,$
    ;field3 : { colorField : 'radmassflux', mapMinMax : [-150,30] } ,$
    field2 : { colorField : 'vrad',        mapMinMax : [-220,220] } };,$
    
  numFields = n_tags(fields)
  numRuns   = n_elements(hInds)
  
  ; load runs
  gcIDs = []
  
  foreach hInd,hInds,i do begin
    sP = mod_struct(sP, 'sP'+str(i), simParams(run='zoom_20Mpc',res=resLevel,hInd=hInd,redshift=redshift))
    gcIDs = [ gcIDs, zoomTargetHalo(sP=sP.(i)) ]
  endforeach
  
  runsStr = sP.(0).saveTag + '_' + str(sP.(0).snap) + '_h' + strjoin(str(hInds))
  
  saveFilename = 'zoomMaps_'+str(sP.(0).res)+'_'+runsStr+'_axes'+str(axisPair[0])+str(axisPair[1])+'_'+$
                  'sf' + str(fix(sizeFacMap*10)) + '_px' + str(nPixels[0]) + '_nF' +str(numFields)+'.eps'
                 
  xs = xySize*numFields
  ys = xySize*numRuns*(1/(1.0-barAreaHeight))
  
  start_PS, sP.(0).plotPath + saveFilename, xs=xs, ys=ys
  
  ; loop over runs
  for i=0,numRuns-1 do begin
  
    print,sP.(i).saveTag + ' ['+str(gcIDs[i])+'] Mapping axes ['+str(axisPair[0])+','+str(axisPair[1])+']'
  
    ; spatial cutouts
    mapCutout = cosmoVisCutout(sP=sP.(i),gcInd=gcIDs[i],sizeFac=sizeFacMap)
    
    ; enhance hsml and make boxsize smaller for map cutout
    mapCutout.loc_hsml *= hsmlFac
    mapCutout.boxSizeImg *= 0.95
    
    for j=0,numFields-1 do begin
  
      config = {saveFilename:'',nPixels:nPixels,axes:axisPair,fieldMinMax:[0,0],$
                gcID:gcIDs[i],haloMass:mapCutout.haloMass,haloVirRad:mapCutout.haloVirRad,$
                boxCen:[0,0,0],boxSizeImg:mapCutout.boxSizeImg,rVirCircs:rVirCircs,$
                ctNameScat:'',ctNameMap:'',sP:sP.(i),bartype:'',scaleBarLen:200.0,$
                secondCutVal:-1,secondText:'',nbottom:0,secondGt:1,singleColorScale:1,$
                barAreaHeight:barAreaHeight,$
                colorField:fields.(j).colorField,mapMinMax:fields.(j).mapMinMax}
            
      print,' '+fields.(j).colorField
            
      sub = cosmoVisCutoutSub(cutout=mapCutout,mapCutout=mapCutout,config=config)
      
      ; plot
      plotMultiSphmap, map=sub, config=config, row=[i,numRuns], col=[j,numFields]
  
    endfor ; numFields,j
  endfor ; numRuns,i
  
  end_PS, density=ceil(nPixels[0]/xySize), pngResize=100 ;, /deletePS
  
end
