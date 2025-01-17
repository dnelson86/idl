; cosmoVisMap.pro
; cosmological boxes - 2d visualization based on sph kernel maps
; dnelson aug.2014

; sphMapBox: run sph kernel density projection on whole box

pro sphMapBox

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  res      = 256
  run      = 'feedback_noFB'
  redshift = 2.0
  partType = 'gas' ;tracervel, tracermc, dm

  ; plot config
  nPixels = [800,800] ;px
  zoomFac = 1    ; only in axes, not along projection direction
  nNGB    = 64   ; use CalcHSML for HSML with nNGB
  axes    = [0,1] ; x,y

  ; paths and render config
  sP = simParams(res=res,run=run,redshift=redshift)
  h = loadSnapshotHeader(sP=sP)
  
  boxSizeImg = [h.boxSize,h.boxSize,h.boxSize]              ;kpc
  boxCen     = [h.boxSize/2.0,h.boxSize/2.0,h.boxSize/2.0]  ;kpc
  
  foreach k,axes do boxSizeImg[k] /= zoomFac
  
  outFilename = 'sphmap.'+sP.saveTag+'.box_'+str(zoomFac)+'.nNGB='+str(nNGB)+'.snap='+str(sP.snap)+$
                '.box.axis0='+str(axes[0])+'.axis1='+str(axes[0])+'.'+partType+'.sav'

  ; save/restore
  if (file_test(sP.derivPath + outFilename)) then begin
    restore,sP.derivPath + outFilename,/verbose
  endif else begin

    ; gas cells or sph particles (pos,hsml stored in snapshots)
    if partType eq 'gas' or partType eq 'stars' or partType eq 'dm' then begin
      mass = loadSnapshotSubset(sP=sP,partType=partType,field='mass')
      pos  = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
      hsml = loadSnapshotSubset(sP=sP,partType=partType,field='hsml')
    endif
    
    ; velocity tracers (pos stored in snapshots, calculate hsmls and use constant eff mass)
    if partType eq 'tracervel' then begin
      mass = replicate(sP.targetGasMass, h.nPartTot[partTypeNum(partType)])
      pos  = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
      hsml = calcHSML(pos,ndims=3,nNGB=nNGB,boxSize=h.boxSize)
    endif
    
    ; monte carlo tracers (use gas pos,hsml and replace mass by num_child_tr*tr_mass_eff)
    if partType eq 'tracermc' then begin
      mass = loadSnapshotSubset(sP=sP,partType='gas',field='numtr')
      mass *= sP.trMassConst
      pos  = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
      hsml = loadSnapshotSubset(sP=sP,partType='gas',field='hsml')
    endif
    
    quant = replicate(1.0, h.nPartTot[partTypeNum(partType)]) ; dummy quant for now
    
    sphMap = calcSphMap(pos,hsml,mass,quant,boxSizeImg=boxSizeImg,boxSizeSim=h.boxSize,boxCen=boxCen,$
                        nPixels=nPixels,axes=axes,ndims=3)
    colMassMap = sphMap.dens_out
    
    save,colMassMap,filename=sP.derivPath + outFilename
  endelse

  ; rescale
  maxVal = max(colMassMap)/2.0
  ;maxVal = 0.5
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
  xMinMax = [boxCen[0]-boxSizeImg[0]/2.0,boxCen[0]+boxSizeImg[0]/2.0]
  yMinMax = [boxCen[1]-boxSizeImg[1]/2.0,boxCen[1]+boxSizeImg[1]/2.0]
  
  start_PS, sP.plotPath + strmid(outFilename,0,strlen(outFilename)-4) + '.eps'
    loadColorTable,'helix'
    tvim,colMassMap,xrange=xMinMax,yrange=yMinMax,pcharsize=0.0001
  end_PS, pngResize=68, /deletePS
  
end

; sphMapSubboxSBS: run sph kernel density projection on subboxes of two runs, side by side

pro sphMapSubboxSBS, jobNum=jobNum, totJobs=totJobs
  compile_opt idl2, hidden, strictarr, strictarrsubs
  if n_elements(jobNum) eq 0 or n_elements(totJobs) eq 0 then message,'Error'

  ; config
  res      = 512
  runs     = ['tracer','feedback']

  subboxMinMax = { x : [3500,7500] ,$
                   y : [5000,9000] ,$
                   z : [5500,9500] }
  valMinMax = [3.5,7.5] ; log K (gas temp)
  
  ; plot config
  nPixels     = [960,960] ;px per side
  zoomFac     = 1    ; only in axes, not along projection direction
  axes        = [0,1] ; x,y
  scaleBarLen = 500.0 ; ckpc
  ctName      = 'blue-red2'
  savePath    = '/n/home07/dnelson/data5/frames/'

  ; setup render config
  boxSizeImg = [ subboxMinMax.x[1]-subboxMinMax.x[0] ,$
                 subboxMinMax.y[1]-subboxMinMax.y[0] ,$
                 subboxMinMax.z[1]-subboxMinMax.z[0] ]
  ;if boxSizeImg[0] ne boxSizeImg[1] or boxSizeImg[0] ne boxSizeImg[2] then message,'Error'
  boxCen     = [mean(subboxMinMax.x), mean(subboxMinMax.y), mean(subboxMinMax.z)]
  
  foreach k,axes do boxSizeImg[k] /= zoomFac
    
  xMinMax = [boxCen[0]-boxSizeImg[0]/2.0,boxCen[0]+boxSizeImg[0]/2.0]
  yMinMax = [boxCen[1]-boxSizeImg[1]/2.0,boxCen[1]+boxSizeImg[1]/2.0]
    
  ; decide on subbox snapshots
  foreach run,runs,i do begin
    sP = simParams(res=res,run=run,redshift=0.0)
    h = loadSnapshotHeader(sP=sP)
    subboxSnapRedshifts = mod_struct( subboxSnapRedshifts, run, snapNumToRedshift(sP=sP,/all,/subBox) )
    w = where( subboxSnapRedshifts.(i) ge 0, count )
    subboxCounts = mod_struct( subboxCounts, run, count )
  endforeach
  
  numFrames = max( [subboxCounts.(0), subboxCounts.(1)] )
  print,'Want to render ['+str(numFrames)+'] total frames.'
  
  ; first run is dominant (has more subboxes)
  if subboxCounts.(0) ne numFrames then message,'Error: rework logic.'
  
  frameSnapNums = mod_struct( frameSnapNums, runs[0], lonarr(numFrames) )
  frameSnapNums = mod_struct( frameSnapNums, runs[1], lonarr(numFrames) )
  
  frameSnapNums.(0) = lindgen(numFrames)
  
  for j=0,numFrames-1 do begin
    curRedshift = (subboxSnapRedshifts.(0))[j]
    w = closest( subboxSnapRedshifts.(1), curRedshift )
    frameSnapNums.(1)[j] = w
  endfor
  
  if max(frameSnapNums.(1)) ge subboxCounts.(1) then message,'Error'
  w = where( frameSnapNums.(1)-shift(frameSnapNums.(1),1) ne 1, count )
  print,'Second run has ['+str(count-1)+'] of ['+str(subboxCounts.(1))+'] frames with non-ideal spacing.'
  
  ; loop over frames
  framesPerJob = ceil(numFrames/totJobs)
  frameMM = [framesPerJob*jobNum > 0,framesPerJob*(jobNum+1)-1 < numFrames]
  print,'Job ['+str(jobNum)+'] of ['+str(totJobs)+'] frames: ',frameMM
  
  for frameNum=frameMM[0],frameMM[1],1 do begin
  ;frameNum = 100 ;100,400,1000,1400
  print,' '+str(frameNum)
  
  maps = {}
  
  ; start plot
  start_PS, savePath + 'frame_'+string(frameNum,format='(I04)')+'.eps', $
    xs=nPixels[0]*2.0005/100.0, ys=nPixels[1]/100.0
  loadColorTable,ctName
  
  ; setup paths and render config
  foreach run,runs,i do begin
    ; render setup
    sP = simParams(res=res,run=run,snap=frameSnapNums.(i)[frameNum])
  
    ; save/restore (TEMPORARY)
    ;outFilename = 'subboxMapTEMP_'+sP.saveTag+'_'+string(frameNum,format='(I04)')
    ;if (file_test(sP.plotPath + outFilename + '.sav')) then begin
    ;  restore,sP.plotPath + outFilename + '.sav',/verbose
    ;endif else begin

      ; load gas quantities
      mass  = loadSnapshotSubset(sP=sP,partType='gas',field='mass',/subBox)
      pos   = loadSnapshotSubset(sP=sP,partType='gas',field='pos',/subBox)
      hsml  = loadSnapshotSubset(sP=sP,partType='gas',field='hsml',/subBox)
      quant = loadSnapshotSubset(sP=sP,partType='gas',field='temp',/subBox)
      quant = 10.0^quant ; take out log
      
      sfr = loadSnapshotSubset(sP=sP,partType='gas',field='sfr',/subBox)
      w = where(sfr gt 0.0,count)
      if count gt 0 then quant[w] = 1000.0 ; set temp of SFing gas to ~ISM temperature (get off eEOS)
      
      sphMap = calcSphMap(pos,hsml,mass,quant,$
                          boxSizeImg=boxSizeImg,boxSizeSim=h.boxSize,boxCen=boxCen,$
                          nPixels=nPixels,axes=axes,ndims=3)
      colMap = sphMap.quant_out
      
    ;  save,colMap,filename=sP.plotPath + outFilename + '.sav'
    ;endelse
    
    ; color scaling
    w = where(colMap eq 0, count, comp=ww)
    if count ne 0 then colMap[w] = min(colMap[ww])
    
    colMap = alog10(colMap)
    colMap = colMap > valMinMax[0] < valMinMax[1]
    colMap = (colMap-valMinMax[0])*255.0 / (valMinMax[1]-valMinMax[0]) ;0-255
    colMap = colMap > 0 < 255
    
    ; add to plot
    pos = [0.5*i,0.0,0.5*(i+1),1.0]
    tvim,colMap,xrange=xMinMax,yrange=yMinMax,pos=pos,/noerase,/c_map,/noframe
    
    ; simname
    cgText,0.5*i+0.01,0.97,sP.simName,alignment=0.0,/normal,color='white'
  endforeach ; runs
    
  ; scale bar
  xpos = [xMinMax[1]*0.92,xMinMax[1]*0.92+scaleBarLen]
  ypos = replicate(yMinMax[1]*0.99,2)
  
  cgText,mean(xpos),ypos*0.99,string(scaleBarLen,format='(i3)')+' kpc',$
    alignment=0.5,color=cgColor('white')
  oplot,xpos,ypos,color=cgColor('white'),thick=!p.thick+0.5
  
  ; redshift
  cgText,xMinMax[1]*0.97,yMinMax[0]*1.02,$
    'z = '+str(string(subboxSnapRedshifts.(0)[frameNum],format='(f5.2)')),$
    alignment=0.5,color=cgColor('white')
  
  ; colorbar
  cbar_pos = [0.4,0.02,0.6,0.055]
  loadColorTable,ctName
  cgColorbar,position=cbar_pos,divisions=0,charsize=0.000001,bottom=0,ticklen=0.00001
  
  cgText,0.41,0.03,string(valMinMax[0],format='(f3.1)'),alignment=0.0,/normal,color='white'
  cgText,0.59,0.03,string(valMinMax[1],format='(f3.1)'),alignment=1.0,/normal,color='white'
  cgText,0.50,0.032,textoidl("T_{gas} [_{ }log K_{ }]"),alignment=0.5,/normal,color='black'
  
  ; finish plot
  end_PS, pngResize=100, density=100, /deletePS
  
  endfor ; frameNum
  
end

; plotScatterAndMap(): plot side by side projection results from CalcSphMap

pro plotScatterAndMap, map=sphmap, scatter=scatter, config=config, column=column
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if ~keyword_set(sphmap) or ~keyword_set(scatter) or $
     ~keyword_set(config) or ~keyword_set(column) then $
    message,'Error: Specify all inputs.'
    
  ; decide column configuration
  curCol  = float(column[0])
  totCols = column[1]
  
  pWidth = 1.0 / totCols
  leftOffset = pWidth * curCol
  
  ; UPPER TWO PANELS - scatter

    xMinMax = [-config.boxSizeScat[0]/2.0,config.boxSizeScat[0]/2.0]
    yMinMax = [-config.boxSizeScat[1]/2.0,config.boxSizeScat[1]/2.0]  
  
    ; plot position (normalized)
    posTop = [leftOffset,0.666,leftOffset+pWidth,1.0]
    posBottom = [leftOffset,0.333,leftOffset+pWidth,0.666]
  
    !p.thick = 2.0
    !p.charsize = 0.8
  
    ; fill with black background
    if curCol eq 0 then $
      cgColorfill,[1,1,0,0,1],[1,0,0,1,1],/normal,color=cgColor('black')
  
    ; color table and establish temperature -> color mapping
    loadColorTable,config.ctNameScat

    ; all gas (left panel)
    cgPlot, /nodata, xMinMax, yMinMax, pos=posTop, xs=5, ys=5, /noerase
    
    ; circle at virial radius
    foreach rVirCirc,config.rVirCircs do $
      tvcircle,rVirCirc * config.haloVirRad,0,0,cgColor('light gray'),thick=0.8,/data
    
    ; particle loop for velocity vector plotting
    nCutoutLeft = n_elements(scatter.pos_left[0,*])
    for i=0L,nCutoutLeft-1,config.interval do $
      oplot,[scatter.pos_left[config.axes[0],i],scatter.pos_left2[config.axes[0],i]],$
             [scatter.pos_left[config.axes[1],i],scatter.pos_left2[config.axes[1],i]],$
             line=0,color=scatter.cinds_left[i],thick=config.lineThick
               
    ; dividing lines
    cgPlot,xMinMax,[yMinMax[0],yMinMax[0]],line=0,thick=1.0,color=cgColor('dark gray'),/overplot
    cgPlot,[xMinMax[0],xMinMax[0]],yMinMax,line=0,thick=1.0,color=cgColor('dark gray'),/overplot
               
    ; cold gas / dark matter (right panel)
    cgPlot, /nodata, xMinMax, yMinMax, pos=posBottom, xs=5, ys=5, /noerase

    foreach rVirCirc,config.rVirCircs do $
      tvcircle,rVirCirc * config.haloVirRad,0,0,cgColor('light gray'),thick=0.8,/data
        
    ; particle loop for velocity vector plotting (cold gas only)
    nCutoutRight = n_elements(scatter.pos_right[0,*])
    for i=0L,nCutoutRight-1,config.interval do $
      oplot,[scatter.pos_right[config.axes[0],i],scatter.pos_right2[config.axes[0],i]],$
             [scatter.pos_right[config.axes[1],i],scatter.pos_right2[config.axes[1],i]],$
             line=0,color=scatter.cinds_right[i],thick=config.lineThick
  
    if curCol eq 0 then begin
      len = 100.0 ;kpc
      
      xpos = [xMinMax[0]*0.9,xMinMax[0]*0.9+len]
      ypos = replicate(yMinMax[1]*0.93,2)
      
      ; skip this for subbox plotting (e.g. scatter and map have same scale)
      if config.haloMass ne -1 then begin
        cgText,mean(xpos),ypos*0.9,string(len,format='(i3)')+' kpc',alignment=0.5,color=cgColor('white')
        loadct,0,/silent
        oplot,xpos,ypos,color=cgColor('white'),thick=!p.thick+0.5
        loadColorTable,config.ctNameScat
      endif
      
    endif
    
    ; dividing lines
    cgPlot,[xMinMax[0],xMinMax[0]],yMinMax,line=0,thick=1.0,color=cgColor('dark gray'),/overplot
    
    ; colorbars
    xthick = !x.thick
    ythick = !y.thick
    !x.thick = 1.0
    !y.thick = 1.0
    
    if config.barType eq '1bar' then begin
      labelText = "log T_{gas} [K]"
      
      ; one colorbar (centered)
      offset = 0.333 ; vertical
      x_left  = 0.5 - pWidth*0.4
      x_right = 0.5 + pWidth*0.4
      pos = [offset+0.01,x_left,offset+0.034,x_right]
      
      loadColorTable,config.ctNameScat
      cgColorbar,position=pos,divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
      
      ; colorbar labels
      cbLabels = str(string([fix(10*config.fieldMinMax[0])/10,fix(10*config.fieldMinMax[1])/10],format='(f3.1)'))
      rAdjust = 0.01*(strlen(cbLabels[1])-1)
      
      cgText,0.5,offset+0.018,textoidl(labelText),alignment=0.5,color=cgColor('black'),/normal
      cgText,x_left + 0.024,offset+0.017,cbLabels[0],alignment=0.5,color=cgColor('black'),/normal
      cgText,x_right -0.008-rAdjust,offset+0.017,cbLabels[1],alignment=0.5,color=cgColor('black'),/normal
    endif
    
    !x.thick = xthick
    !y.thick = ythick

  ; LOWER PANEL - sphmap

  ; box
  xMinMax = [config.boxCen[0]-config.boxSizeImg[0]/2.0,config.boxCen[0]+config.boxSizeImg[0]/2.0]
  yMinMax = [config.boxCen[1]-config.boxSizeImg[1]/2.0,config.boxCen[1]+config.boxSizeImg[1]/2.0]
  
  ; plot position (normalized)
  posBottom = [leftOffset,0.0,leftOffset+pWidth,0.333]
  
  ; output image
  loadColorTable,config.ctNameMap
  
  cgPlot, /nodata, xMinMax, yMinMax, pos=posBottom, xs=5, ys=5, /noerase
  tv, sphmap.quant_out,posBottom[0],posBottom[1],/normal,xsize=pWidth ; mass-weighted quantity

  ; circle at virial radius
  foreach rVirCirc,config.rVirCircs do $
    tvcircle,rVirCirc * config.haloVirRad,0,0,cgColor('white'),thick=0.8,/data    
  
  ; scale bar
  if curCol eq 0 then begin
    xpos = [xMinMax[0]*0.9,xMinMax[0]*0.9+config.scaleBarLen]
    ypos = replicate(yMinMax[1]*0.93,2)
    
    cgText,mean(xpos),ypos*0.9,string(config.scaleBarLen,format='(i3)')+' kpc',$
      alignment=0.5,color=cgColor('white')
    loadct,0,/silent
    oplot,xpos,ypos,color=cgColor('white'),thick=!p.thick+0.5
    loadColorTable,'bw linear'
  endif
  
  ; redshift and halo mass
  if config.haloMass ne -1 then $
    cgText,posTop[2]-0.01,0.666-0.02,"M = "+string(config.haloMass,format='(f4.1)'),alignment=1.0,$
           color=cgColor('white'),/normal
             
  if curCol eq totCols-1 then $
    cgText,0.99,0.666-0.04,"z = "+string(config.sP.redshift,format='(f3.1)'),alignment=1.0,$
           color=cgColor('white'),/normal
  
  ; simulation name
  cgText,mean(posTop[[0,2]]),0.97,config.sP.simName,charsize=!p.charsize+0.5,$
    alignment=0.5,/normal,color=cgColor('white')
          
  ; quantity name
  if curCol eq totCols-1 then begin
    cgText,0.99,0.676,"all gas",charsize=!p.charsize+0.3,alignment=1.0,/normal,color=cgColor('white')
    cgText,0.99,0.343,config.secondText,charsize=!p.charsize+0.3,alignment=1.0,/normal,color=cgColor('white')
    cgText,0.99,0.03,"projected",charsize=!p.charsize+0.3,alignment=1.0,/normal,color=cgColor('white')
    cgText,0.99,0.01,config.colorField,charsize=!p.charsize+0.3,alignment=1.0,/normal,color=cgColor('white')
  endif
    
  ; dividing lines
  cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=1.0,color=cgColor('dark gray'),/overplot
  cgPlot,[xMinMax[0],xMinMax[0]],yMinMax,line=0,thick=1.0,color=cgColor('dark gray'),/overplot
          
end

; plotMultiSphmap(): plot side by side projection results from CalcSphMap

pro plotMultiSphmap, map=sphmap, config=config, row=row, col=col
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if ~keyword_set(sphmap) or ~keyword_set(config) or $
    ~keyword_set(row) or ~keyword_set(col) then $
    message,'Error: Specify all inputs.'
   
  ; decide position configuration
  curRow  = float(row[0])
  totRows = row[1]
  
  curCol  = float(col[0])
  totCols = col[1]
  
  rowHeight = (1.0-config.barAreaHeight) / totRows
  colWidth  = 1.0 / totCols
  
  topOffset  = 1.0 - rowHeight * (curRow+1)
  leftOffset = colWidth * (curCol)
  
  ; different plot types
  zoom2x1 = 0
  zoom4x2 = 0
  zoom4x4 = 0
  if totCols eq 2 and totRows eq 1 then zoom2x1 = 1
  if totCols eq 4 and totRows eq 2 then zoom4x2 = 1
  if totCols eq 4 and totRows eq 4 then zoom4x4 = 1
  
  if zoom4x2 then topOffset -= config.barAreaHeight*0.5
  
  ; box
  xMinMax = [config.boxCen[0]-config.boxSizeImg[0]/2.0,config.boxCen[0]+config.boxSizeImg[0]/2.0]
  yMinMax = [config.boxCen[1]-config.boxSizeImg[1]/2.0,config.boxCen[1]+config.boxSizeImg[1]/2.0]
  
  ; plot position (normalized)
  pos = [leftOffset,topOffset,leftOffset+colWidth,topOffset+rowHeight]

  ; output image
  loadColorTable,config.ctNameMap,black0=(tag_exist(config,'black0'))
  
  cgPlot, /nodata, xMinMax, yMinMax, pos=pos, xs=5, ys=5, /noerase
  tv, sphmap.quant_out,pos[0],pos[1],/normal,xsize=colWidth ; mass-weighted quantity

  ; circle at virial radius
  foreach rVirCirc,config.rVirCircs do $
    tvcircle,rVirCirc * config.haloVirRad,0,0,cgColor('white'),thick=0.8,/data    
  
  ; scale bar(s)
  xfac  = 0.9
  yfac  = 0.86
  yfac2 = 0.93
  
  if zoom4x2 then yfac = 0.88 ; zoom 4x2
  if zoom2x1 then yfac = 0.95
  if zoom2x1 then xfac = 0.95
  if zoom2x1 then yfac2 = 0.96
  
  ; convert ckpc/h -> kpc by adjusting size, leaving label
  scalefac = 1.0/(1+config.sP.redshift)
  convFac = 1.0
  sLabel = 'ckpc'
  
  if zoom4x4 then convFac = scalefac / 0.7
  if zoom4x4 then sLabel = 'kpc'
  
  if curCol eq 0 and curRow eq 0 and ~tag_exist(config,'secondScaleBar') then begin
    xpos = [xMinMax[0]*xfac,xMinMax[0]*xfac+config.scaleBarLen/convFac]
    ypos = replicate(yMinMax[1]*yfac2,2)
    
    if config.sP.zoomLevel ne 99 then begin
      cgText,mean(xpos),ypos*yfac,string(config.scaleBarLen,format='(i3)')+' '+sLabel,$
        alignment=0.5,color=cgColor('white')
      oplot,xpos,ypos,color=cgColor('white'),thick=!p.thick+0.5
    endif
  endif
  
  if zoom4x2 or tag_exist(config,'secondScaleBar') then begin
    scaleBarLen = 0
    if zoom4x2 then targetRow = 1
    if tag_exist(config,'secondScaleBar') then targetRow = 0
    
    if curCol eq 0 and curRow eq targetRow then scaleBarLen = config.scaleBarLen * 0.5
    if curCol eq 2 and curRow eq targetRow then scaleBarLen = config.scaleBarLen * 2.0
    
    if scaleBarLen gt 0 then begin
      xpos = [xMinMax[0]*0.9,xMinMax[0]*0.9+scaleBarLen/convFac]
      ypos = replicate(yMinMax[1]*0.93,2)
      
      cgText,mean(xpos),ypos*yfac,string(scaleBarLen,format='(i3)')+' '+sLabel,$
        alignment=0.5,color=cgColor('white')
      oplot,xpos,ypos,color=cgColor('white'),thick=!p.thick+0.5
    endif
  endif
  
  ; redshift and halo mass
  loadct,0,/silent
  
  if totCols eq 4 then begin
    ; zoom 4x4
    if totRows eq 4 and curRow eq 0 and curCol eq totCols-1 then $
      cgText,pos[2]-0.01,pos[3]-0.02,"z = "+string(config.sP.redshift,format='(f3.1)'),alignment=1.0,$
             color=cgColor('white'),/normal
    ; zoom 4x2
    ;if totRows eq 2 and curRow eq 1 and curCol eq totCols-1 then $
    ;  cgText,pos[2]-0.01,pos[1]+0.025,"z = "+string(config.sP.redshift,format='(f3.1)'),alignment=1.0,$
    ;         color=cgColor('white'),/normal           
  endif else begin
    if zoom2x1 ne 1 then begin
      if (curRow eq 0 and curCol eq totCols-1) or (curRow eq 0 and totCols gt 3) then $
        cgText,pos[2]-0.01,pos[1]+0.015,"M = "+string(config.haloMass,format='(f4.1)'),alignment=1.0,$
               color=cgColor('white'),/normal
               
      if (curRow eq 1 and curCol eq totCols-1) or (curRow eq 1 and totCols gt 3) then $
        cgText,pos[2]-0.01,pos[3]-0.03,"z = "+string(config.sP.redshift,format='(f3.1)'),alignment=1.0,$
               color=cgColor('white'),/normal
    endif
  endelse
  
  ; simulation name
  if curCol eq 0 and ~zoom4x2 and ~zoom2x1 and ~zoom4x4 then $
    cgText,0.03-0.01*(totCols gt 3),mean(pos[[1,3]]),config.sP.simName,charsize=!p.charsize+0.2,$
      alignment=0.5,/normal,orientation=90,color=cgColor('white')
      
  if config.sP.zoomLevel ne 99 then begin
    if curCol eq 0 and zoom4x4 then $
      cgText,0.03-0.01*(totCols gt 3),mean(pos[[1,3]]),$
        'h'+str(config.sP.hIndDisp)+'L'+str(config.sP.levelMax),$
        charsize=!p.charsize+0.2,alignment=0.5,/normal,orientation=90,color=cgColor('white')
      
  endif else begin
    ; iClusters
    if curCol eq 0 then $
      cgText,0.03-0.01*(totCols gt 3),mean(pos[[1,3]]),$
        config.sP.simName+' M='+string(config.haloMass,format='(f4.1)'),$
        charsize=!p.charsize+0.2,alignment=0.5,/normal,orientation=90,color=cgColor('white')
  endelse

  if zoom2x1 then $
    cgText,0.965,0.97,config.sP.simName+'z'+str(fix(config.sP.redshift)),$
      charsize=!p.charsize+0.1,alignment=0.5,/normal,color='white'
      
  ; quantity name
  ;if curRow eq 0 then $
  ;  cgText,mean(pos[[0,2]]),0.97,config.colorField,$
  ;    charsize=!p.charsize,alignment=0.5,/normal,color=cgColor('white')
    
  ; dividing lines
  thick = 5.0
  
  if zoom2x1 ne 1 then begin
    cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=thick,color=getColor24([254,254,254]),/overplot
  
    if curCol gt 0 then $
      cgPlot,[xMinMax[0],xMinMax[0]],yMinMax,line=0,thick=thick,color=getColor24([254,254,254]),/overplot
            
    if curRow eq totRows-1 then $
      cgPlot,xMinMax,[yMinMax[0],yMinMax[0]],line=0,thick=thick,color=getColor24([254,254,254]),/overplot
  endif
          
  ; colorbar(s)
  if curRow eq 0 or (zoom4x2 and curRow eq 1) then begin
    xthick = !x.thick
    ythick = !y.thick
    !x.thick = 1.0
    !y.thick = 1.0
    
    if config.colorField eq 'temp'     then labelText = "log T_{gas} [_{ }K_{ }]"
    if config.colorField eq 'temptvir' then labelText = "log ( T_{gas} / T_{vir} )"
    if config.colorField eq 'entropy'  then labelText = "log S [_{ }K cm^{2 }]"
    if config.colorField eq 'entnorm'  then labelText = "log ( S / S_{200} )"
    if config.colorField eq 'overdens' then labelText = "log ( \rho_{gas} / \rho_{crit,b} )"
    if config.colorField eq 'vrad'     then labelText = "v_{rad} [_{ }km/s_{ }]"
    if config.colorField eq 'vradnorm' then labelText = "v_{rad} / v_{200}"
    if config.colorField eq 'dmdens'   then labelText = "log ( \rho_{dm} / \rho_{crit} )"
    if config.colorField eq 'stardens' then labelText = "log ( \rho_{stars} / \rho_{crit,b} )"
    
    factor = 0.48 ; bar length, 0.5=whole
    height = 0.04 ; colorbar height, fraction of entire figure
    hOffset = 0.4 ; padding between image and start of bar (fraction of height)
    mOffset = 0.008
    yOffset = 0.75
    tOffset = 0.01 ; move colorbar text numbers outside?
    align   = [0.0,1.0]
    
    if zoom4x4 then begin
      height = 0.02
      mOffset = 0.005
      yOffset = 0.9
    endif
    
    if zoom4x2 then begin
      height  = 0.03
      mOffset = 0.007
      yOffset = -0.3 ;0.9
      factor  = 0.34
      tOffset = -0.006
      align   = [1.0,0.0]
    endif
    
    if zoom2x1 then begin
      height  = 0.03
      mOffset = 0.007
      yOffset = -0.3 ;0.9
      factor  = 0.25
      ;tOffset = -0.006
      ;align   = [1.0,0.0]
    endif
    
    x_left  = mean(pos[[0,2]]) - (pos[2]-pos[0])*factor
    x_right = mean(pos[[0,2]]) + (pos[2]-pos[0])*factor
    y_bottom = config.barAreaHeight - height*hOffset - height
    y_top    = config.barAreaHeight - height*hOffset
    
    if zoom4x2 then begin
      if curRow eq 0 then begin ; top scalebars
        y_bottom = pos[3] + height*hOffset ;config.barAreaHeight*0.5
        y_top    = y_bottom + height ; config.barAreaHeight*0.5
      endif
      if curRow eq 1 then begin ; bottom scalebars
        y_bottom -= config.barAreaHeight*0.5
        y_top    -= config.barAreaHeight*0.5
      endif
    endif
    
    if zoom2x1 then begin
      y_bottom = 0.01
      y_top    = y_bottom + height
    endif
    
    ; just one color bar for vs. redshift plot
    if totCols eq 3 then begin
      x_left = 0.25
      x_right = 0.75
    endif
    
    pos = [y_bottom,x_left,y_top,x_right]
    
    loadColorTable,config.ctNameMap,black0=(tag_exist(config,'black0'))
    cgColorbar,position=pos,divisions=0,charsize=0.000001,bottom=0,ticklen=0.00001
    
    ; colorbar labels
    if config.colorField eq 'temp' or config.colorField eq 'entropy' then $
      cbLabels = str(string(config.mapMinMax,format='(f3.1)'))
    if config.colorField eq 'vrad' then $
      cbLabels = str(config.mapMinMax)
    if config.colorField eq 'overdens' or config.colorField eq 'temptvir' or $
       config.colorField eq 'dmdens' or config.colorField eq 'stardens' or $
       config.colorfield eq 'vradnorm' or config.colorField eq 'entnorm' then $
      cbLabels = str(string(config.mapMinMax,format='(f4.1)'))
    
    if zoom4x2 or zoom2x1 then cbLabels = strtrim(cbLabels) ; since outside
    
    y_middle = mean([y_bottom,y_top]) - mOffset
    cgText,mean([x_left,x_right]),y_bottom-height*yOffset,$
      textoidl(labelText),alignment=0.5,color=cgColor('black'),charsize=!p.charsize,/normal
    cgText,x_left  + tOffset,y_middle,cbLabels[0],alignment=align[0],color=cgColor('black'),/normal
    cgText,x_right - tOffset,y_middle,cbLabels[1],alignment=align[1],color=cgColor('black'),/normal
    
    !x.thick = xthick
    !y.thick = ythick
  endif
          
end

; -----------------------------------------------------------------------------

; sphScatterAndMapHaloComp(): compare scatterplots and density/temperature of arepo/gadget in 3x2 panels

pro sphScatterAndMapHaloComp

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  redshift = 2.0
  res      = [9,10,11]
  runs     = ['zoom_20mpc','zoom_20mpc','zoom_20mpc']
    ;['zoom_20mpc','zoom_20mpc_derefgal','zoom_20mpc_derefgal_nomod']
    ;['feedback','feedback_noZ','feedback_noFB']
    ;['gadget','tracer','feedback']  ; two or three supported
                                                
  hInds    = [0,0,0] ; for zooms only
  haloID   = 0 ;zoom.0 z2.304 z2.301 z2.130 z2.64
  
  sizeFacMap  = 3.5       ; times rvir (image width, e.g. 2.0 shows out to the virial radius)
  sizeFacScat = 3.5       ; times rvir
  rVirCircs   = [0.15,0.5,1.0] ; times rvir to draw a circle
  hsmlFac  = 1.5          ; times each cell hsml for sph projections
  nPixels  = [800,800]    ; px
  xySize   = 5            ; final image is xySize*nPixels[0] high
  axisPair = [0,1]        ; xy
  nbottom  = 50           ; lift minimum scatterplot color from zero (to distinguish from black)
  secondGt = 0            ; 1=show greater than cut, 0=show less than cut
  singleColorScale = 1    ; 1=use same color scale for right panel, 0=rescale
  plotAllPts       = 1    ; if 0 skip some points in higher res runs to match visual density of L9
  
  ; use which field and minmax for color mappings? cut value for right panel?
  colorField = 'temp' & fieldMinMax = [4.0,6.2] & mapMinMax = [4.4,6.2] & secondCutVal = 5.0 ; log(K)
  ;colorField = 'entropy' & fieldMinMax = [5.5,8.8] & mapMinMax = [6.5,8.9] & secondCutVal = 7.5 ; log(CGS)
  ;colorField = 'vrad' & fieldMinMax = [-400,400] & mapMinMax = [-350,350] & secondCutVal = -300.0 ; km/s
  ;colorField = 'radmassflux' & fieldMinMax = [-500,50] & mapMinMax = [-150,30] & secondCutVal = -150.0
  ;colorField = 'radmassfluxSA' & fieldMinMax = [-130,50] & mapMinMax = [-100,30] & secondCutVal = -20.0
  
  ; load runs
  runsStr = ''
  gcIDs = []
  foreach run,runs,i do begin
    sP = mod_struct(sP, 'sP'+str(i), simParams(res=res[i],run=run,redshift=redshift,hInd=hInds[i]))
    gcIDs = [ gcIDs, getMatchedIDs(simParams=sP.(i),haloID=haloID) ]
    runsStr += sP.(i).saveTag + '.' + str(sP.(i).snap) + '.h' + str(gcIDs[i]) + '-'
  endforeach
  
  numRuns = n_tags(sP)
  
  saveFilename = 'sc.map.'+str(sP.(0).res)+'.'+runsStr+'axes'+str(axisPair[0])+str(axisPair[1])+'_'+$
                  'sf' + str(fix(sizeFacMap*10)) + '_px' + str(nPixels[0]) + '_all' + $
                  str(plotAllPts) + '_' + colorField+'.eps'
                 
  start_PS, sP.(0).plotPath + saveFilename, xs=xySize*numRuns, ys=xySize*3 
  
  ; loop over runs
  for i=0,numRuns-1 do begin
  
    print,sP.(i).saveTag + ' ['+str(gcIDs[i])+'] Mapping axes ['+str(axisPair[0])+','+str(axisPair[1])+']'
  
    ; spatial cutouts
    mapCutout = cosmoVisCutout(sP=sP.(i),gcInd=gcIDs[i],sizeFac=sizeFacMap)
    cutout    = cosmoVisCutout(sP=sP.(i),gcInd=gcIDs[i],sizeFac=sizeFacScat)

    ; enhance hsml and make boxsize smaller for map cutout
    mapCutout.loc_hsml *= hsmlFac
    mapCutout.boxSizeImg *= 0.95
  
    if plotAllPts eq 1 then begin
      ; config (1): plot all points, reduce line thickness
      interval = 1
      lineThick = 1.5 / 2.0^(sP.(i).zoomLevel-2)
    endif else begin
      ; config (1): plot only every 8th/64th point for L10/L11 (match visual density)
      interval = 8^(sP.(i).zoomLevel-2)
      lineThick = 1.0
    endelse
    
    config = {saveFilename:'',nPixels:nPixels,axes:axisPair,fieldMinMax:fieldMinMax,$
              gcID:gcIDs[i],haloMass:cutout.haloMass,haloVirRad:cutout.haloVirRad,$
              boxCen:[0,0,0],boxSizeImg:mapCutout.boxSizeImg,boxSizeScat:cutout.boxSizeImg,$
              ctNameScat:'helix',ctNameMap:'',sP:sP.(i),bartype:'1bar',scaleBarLen:200.0,$
              lineThick:lineThick,interval:interval,colorField:colorField,secondText:'',$
              nbottom:nbottom,secondGt:secondGt,singleColorScale:singleColorScale,$
              secondCutVal:secondCutVal,mapMinMax:mapMinMax,rVirCircs:rVirCircs}
          
    print,' boxSize ['+str(cutout.boxSizeImg[0])+'] kpc around subhalo center ['+$
          str(cutout.boxCen[0])+' '+str(cutout.boxCen[1])+' '+str(cutout.boxCen[2])+'] interval: '+str(interval)
          
    sub = cosmoVisCutoutSub(cutout=cutout,mapCutout=mapCutout,config=config)
    
    ; plot
    plotScatterAndMap, map=sub, scatter=sub, config=config, column=[i,numRuns]
  
  endfor
  
  end_PS, density=ceil(nPixels[0]/xySize), pngResize=100, /deletePS
  
  stop

end

; multiMapHaloComp(): compare density/temperature/... of two runs in Nx2 panels

pro multiMapHaloComp

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  redshift = 2.0
  res      = [512,512]
  runs     = ['tracer','feedback']

  hInds    = [0,1,2,3] ; for zooms only
  haloID   = 304 ;0 ;zoom.0 z2.304 z2.301 z2.130 z2.64

  sizeFacMap  = 3.5           ; times rvir (image width, e.g. 2.0 shows out to the virial radius)
  rVirCircs   = [0.15,0.5,1.0] ; times rvir to draw a circle
  hsmlFac     = 1.5          ; times each cell hsml for sph projections
  nPixels     = [800,800]    ; px
  xySize      = 3            ; final image is xySize*nPixels[0] high
  axisPair    = [0,1]        ; xy, xz
  barAreaHeight = 0.1        ; fractional
  
  ; use which field and minmax for color mappings?
  fields = { $
    field4 : { colorField : 'overdens',    mapMinMax : [-1.0,4.0] } ,$
    field0 : { colorField : 'temp',        mapMinMax : [4.4,6.2] } ,$
    ;field1 : { colorField : 'entropy',     mapMinMax : [6.5,8.9] } ,$
    ;field3 : { colorField : 'radmassflux', mapMinMax : [-150,30] } ,$
    field2 : { colorField : 'vrad',        mapMinMax : [-180,180] } };,$
    
  numFields = n_tags(fields)
  
  ; load runs
  runsStr = ''
  gcIDs = []
  foreach run,runs,i do begin
    sP = mod_struct(sP, 'sP'+str(i), simParams(res=res[i],run=run,redshift=redshift,hInd=hInds[i]))
    gcIDs = [ gcIDs, getMatchedIDs(simParams=sP.(i),haloID=haloID) ]
    runsStr += sP.(i).saveTag + '.' + str(sP.(i).snap) + '.h' + str(gcIDs[i]) + '-'
  endforeach
  stop
  numRuns = n_tags(sP)
  
  saveFilename = 'sc.map.'+str(sP.(0).res)+'.'+runsStr+'axes'+str(axisPair[0])+str(axisPair[1])+'_'+$
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

; multiRedshiftHaloComp(): compare one field between two runs at multiple redshifts in Nx2 panels

pro multiRedshiftHaloComp

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  res       = 512
  runs      = ['tracer','feedback']
  redshifts = [3.0,2.0,1.0,0.0]

  hInds    = [0,1,2,3] ; for zooms only
  haloID   = 301 ;0 ;zoom.0 z2.304 z2.301 z2.130 z2.64

  sizeFacMap  = 3.5           ; times rvir (image width, e.g. 2.0 shows out to the virial radius)
  rVirCircs   = [0.15,0.5,1.0] ; times rvir to draw a circle
  hsmlFac     = 1.5          ; times each cell hsml for sph projections
  nPixels     = [800,800]    ; px
  xySize      = 3            ; final image is xySize*nPixels[0] high
  axisPair    = [1,2]        ; xy, xz
  barAreaHeight = 0.1        ; fractional
  
  ; use which field and minmax for color mappings?
  ;field = { colorField : 'temp',      mapMinMax : [4.4,6.8]}
  field = { colorField : 'temptvir',  mapMinMax : [-1.0,2.0]}
  ;field = { colorField : 'overdens',  mapMinMax : [-1.0,4.0] }
  ;field = { colorField : 'entropy',   mapMinMax : [6.5,8.9] }
  ;field = { colorField : 'vrad',      mapMinMax : [-180,180] }
    
  ; load runs
  runsStr = ''
  foreach run,runs,i do begin
    sP = simParams(res=res,run=run,redshift=redshifts[0])
    runsStr += sP.saveTag + '-'
  endforeach
  
  numRuns = n_elements(runs)
  numZ = n_elements(redshifts)
  
  saveFilename = 'sc.map.'+str(sP.res)+'.'+runsStr+field.colorField+$
                 '_axes'+str(axisPair[0])+str(axisPair[1])+'_'+$
                 'sf' + str(fix(sizeFacMap*10)) + '_px' + str(nPixels[0]) + '_nZ' +str(numZ)+'.eps'
                 
  xs = xySize*numZ
  ys = xySize*numRuns*(1/(1.0-barAreaHeight))
  
  start_PS, sP.plotPath + saveFilename, xs=xs, ys=ys
  
  ; loop over runs
  for i=0,numRuns-1 do begin
    for j=0,numZ-1 do begin
  
      sP = simParams(res=res,run=runs[i],redshift=redshifts[j])
      
      gcID = getMatchedIDs(simParams=sP,haloID=haloID)
      
      ; spatial cutouts
      mapCutout = cosmoVisCutout(sP=sP,gcInd=gcID,sizeFac=sizeFacMap)
      
      ; enhance hsml and make boxsize smaller for map cutout
      mapCutout.loc_hsml *= hsmlFac
      mapCutout.boxSizeImg *= 0.95
      
      print,sP.saveTag + ' (z='+string(redshifts[j],format='(f3.1)')+') ['+str(gcID)+$
            '] Mapping axes ['+str(axisPair[0])+','+str(axisPair[1])+']'
  
      config = {saveFilename:'',nPixels:nPixels,axes:axisPair,fieldMinMax:[0,0],$
                gcID:gcID,haloMass:mapCutout.haloMass,haloVirRad:mapCutout.haloVirRad,$
                boxCen:[0,0,0],boxSizeImg:mapCutout.boxSizeImg,rVirCircs:rVirCircs,$
                ctNameScat:'',ctNameMap:'',sP:sP,bartype:'',scaleBarLen:200.0,$
                secondCutVal:-1,secondText:'',nbottom:0,secondGt:1,singleColorScale:1,$
                barAreaHeight:barAreaHeight,$
                colorField:field.colorField,mapMinMax:field.mapMinMax}
            
      sub = cosmoVisCutoutSub(cutout=mapCutout,mapCutout=mapCutout,config=config)
      
      ; plot
      plotMultiSphmap, map=sub, config=config, row=[i,numRuns], col=[j,numZ]
  
    endfor ; numZ,j
  endfor ; numRuns,i
  
  end_PS, density=ceil(nPixels[0]/xySize), pngResize=100 ;, /deletePS

end

; sphMapHalosDM: run sph kernel density projection on DM particles with boxes centered on halos

pro sphMapHalosDM, sP=sP, gcIDs=gcIDs

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  sizeFac = 3.5      ; times rvir
  nPixels = [800,800] ; px
  
  ;axes = list([0,1],[0,2],[1,2]) ;xy,xz,yz
  axes = list([0,1])
  
  ; target list
  gc    = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sP) 

  ; load dm positions from snapshot and replicate masses from header
  pos  = loadSnapshotSubset(sP=sP,partType='dm',field='pos')
  
  h = loadSnapshotHeader(sP=sP)
  dmPartMass = float(h.massTable[partTypeNum('dm')])
  mass = replicate(dmPartMass,h.nPartTot[1])
  
  ; load HSMLs and veldisps (~temp), assume ID sorted with no gaps (though DM need not be at the end)
  hsmldir = loadHsmlDir(sP=sP,/readHsml,/readVelDisp,/verbose)
  ids     = loadSnapshotSubset(sP=sP,partType='dm',field='ids')
  hsml    = hsmldir.hsml[ids-1]
  veldisp = hsmldir.veldisp[ids-1]
  
  hsmldir = !NULL
  ids = !NULL
  
  ; loop over all non-background subhalos and image
  foreach gcID, gcIDs do begin
  
    ; get subhalo position and size of imaging box
    haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcID]] ;ckpc
    haloMass   = codeMassToLogMsun(gc.subgroupMass[gcID])
    
    boxCen     = sgcen[*,gcID]
    boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcID]])
    boxSizeImg = [boxSize,boxSize,boxSize] ; cube
  
    foreach axisPair, axes do begin
  
      saveFilename = 'dmMap.'+sP.saveTag+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+$
                     '.axes'+str(axisPair[0])+str(axisPair[1])+'.sav'
                    
      if ~file_test(sP.derivPath+'sphMaps/'+saveFilename) then begin
      
        ; get box center (in terms of specified axes)
        boxCenImg  = [sgcen[axisPair[0],gcID],sgcen[axisPair[1],gcID],sgcen[3-axisPair[0]-axisPair[1],gcID]]
        
        print,'['+string(gcID,format='(i4)')+'] Mapping ['+str(axisPair[0])+' '+$
              str(axisPair[1])+'] with '+str(boxSize[0])+$
              ' kpc box around subhalo center ['+str(boxCen[0])+' '+str(boxCen[1])+' '+str(boxCen[2])+']'

        ; calculate projection using sph kernel
        dmMap = calcSphMap(pos,hsml,mass,veldisp,$
                           boxSizeImg=boxSizeImg,boxSizeSim=sP.boxSize,boxCen=boxCenImg,$
                           nPixels=nPixels,axes=axisPair,ndims=3)
                            
        ; vertical+horizontal flip the C image arrays so they agree with the (x,y) order in scatterMapHalos
        dmMap.dens_out  = reverse(dmMap.dens_out,2)
        dmMap.dens_out  = reverse(dmMap.dens_out,1)
        dmMap.quant_out = reverse(dmMap.quant_out,2)
        dmMap.quant_out = reverse(dmMap.quant_out,1)
        
        ; save
        config = {saveFilename:saveFilename,sizeFac:sizeFac,nPixels:nPixels,axes:axes,boxCen:boxCen,$
                  gcID:gcID,haloMass:haloMass,haloVirRad:haloVirRad,$
                  boxCenImg:boxCenImg,boxSize:boxSize,boxSizeImg:boxSizeImg,sP:sP}
        ;save,dmMap,config,filename=sP.derivPath+'sphMaps/'+saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
      endif else begin
        restore,sP.derivPath +'sphMaps/'+ saveFilename
      endelse
    
      ; plot
      plotSphmapDensQuant, map=dmMap, config=config

    endforeach ;axisPair
  
  endforeach ;gcIDs
  
end
