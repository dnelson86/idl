; cosmoVisMap.pro
; cosmological boxes - 2d visualization based on sph kernel maps
; dnelson sep.2012

; sphMapBox: run sph kernel density projection on whole box

pro sphMapBox, partType=partType

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  res = 128
  run = 'tracernew'
  ;partType = 'tracervel'

  redshift = 0.0
  
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
  
  outFilename = 'sphmap.box_'+str(zoomFac)+'.nNGB='+str(nNGB)+'.snap='+str(sP.snap)+$
                '.box.axis0='+str(axes[0])+'.axis1='+str(axes[0])+'.'+partType

  ; save/restore
  if (file_test(sP.derivPath + outFilename + '.sav')) then begin
    restore,sP.derivPath + outFilename + '.sav',/verbose
  endif else begin

    ; gas cells or sph particles (pos,hsml stored in snapshots)
    if partType eq 'gas' then begin
      mass = loadSnapshotSubset(sP=sP,partType=partType,field='mass',/verbose)
      pos  = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
      hsml = loadSnapshotSubset(sP=sP,partType=partType,field='hsml')
      quant = replicate(1.0,h.nPartTot[0]) ; dummy quant for now
    endif
    
    ; velocity tracers (pos stored in snapshots, calculate hsmls and use constant eff mass)
    if partType eq 'tracervel' then begin
      mass = replicate(sP.targetGasMass, h.nPartTot[2])
      pos  = loadSnapshotSubset(sP=sP,partType=partType,field='pos',/verbose)
      hsml = calcHSML(pos,ndims=3,nNGB=nNGB,boxSize=h.boxSize)
      quant = replicate(1.0,h.nPartTot[2]) ; dummy quant for now
    endif
    
    ; monte carlo tracers (use gas pos,hsml and replace mass by num_child_tr*tr_mass_eff)
    if partType eq 'tracermc' then begin
      mass = loadSnapshotSubset(sP=sP,partType='gas',field='numtr',/verbose)
      mass *= sP.trMassConst
      pos  = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
      hsml = loadSnapshotSubset(sP=sP,partType='gas',field='hsml')
      quant = replicate(1.0,h.nPartTot[0]) ; dummy quant for now
    endif
    
    sphMap = calcSphMap(pos,hsml,mass,quant,boxSizeImg=boxSizeImg,boxSizeSim=h.boxSize,boxCen=boxCen,$
                        nPixels=nPixels,axes=axes,ndims=3)
    colMassMap = sphMap.dens_out
    
    save,colMassMap,filename=sP.derivPath + outFilename + '.sav'
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
  
  start_PS, sP.plotPath + outFilename + '.eps'
    loadct, 4, bottom=1, /silent
    tvim,colMassMap,xrange=xMinMax,yrange=yMinMax,pcharsize=0.0001
  end_PS, pngResize=68, /deletePS
  
end

; plotSphmapDensQuant(): plot side by side projection results from CalcSphMap

pro plotSphmapDensQuant, map=sphmap, config=config, top=top, bottom=bottom

  compile_opt idl2, hidden, strictarr, strictarrsubs

  highmass = 1
  
  if highmass eq 1 then clipFacs = [0.10,400.0,1.05,280.0]
  if highmass eq 0 then clipFacs = [1.5,255.0,1.0,255.0]
  
  ; dens: color map and rescale
  w = where(sphmap.dens_out le 0.0,count,comp=wc)

  sphmap.dens_out = alog10(sphmap.dens_out)
  if highmass eq 1 then begin
    if count gt 0 then sphmap.dens_out[w] = 0.0
    sphmap.dens_out = sphmap.dens_out-clipFacs[0]
  endif
  if highmass eq 0 then begin
    if count gt 0 then sphmap.dens_out[w] = 0.0
    sphmap.dens_out = sphmap.dens_out-clipFacs[0]*min(sphmap.dens_out)
  endif
  
  w = where(sphmap.dens_out lt 0.0,count)
  if count gt 0 then sphmap.dens_out[w] = 0.0
  
  sphmap.dens_out = sphmap.dens_out*clipFacs[1] / (max(sphmap.dens_out)-min(sphmap.dens_out)) < 255 ; 0-255

  ; temp
  w = where(sphmap.quant_out eq 0.0,count,comp=wc)
  if count gt 0 then sphmap.quant_out[w] = min(sphmap.quant_out[wc])
  
  sphmap.quant_out = alog10(sphmap.quant_out)
  sphmap.quant_out = sphmap.quant_out-clipFacs[2]*min(sphmap.quant_out)
  
  w = where(sphmap.quant_out lt 0.0,count)
  if count gt 0 then sphmap.quant_out[w] = 0.0
  
  sphmap.quant_out = sphmap.quant_out*clipFacs[3] / (max(sphmap.quant_out)-min(sphmap.quant_out)) < 255 ;0-255

  ; box
  xMinMax = [config.boxCen[0]-config.boxSizeImg[0]/2.0,config.boxCen[0]+config.boxSizeImg[0]/2.0]
  yMinMax = [config.boxCen[1]-config.boxSizeImg[1]/2.0,config.boxCen[1]+config.boxSizeImg[1]/2.0]
  
  plotPath = '/n/home07/dnelson/coldflows/'
  
  ; plot position (normalized)
  posLeft = [0.0,0.0,0.5,1.0]
  if keyword_set(top) then posLeft = [0.0,0.5,0.5,1.0]
  if keyword_set(bottom) then posLeft = [0.0,0.0,0.5,0.5]
  posRight = posLeft + [0.5,0.0,0.5,0.0]
  
  if ~keyword_set(top) and ~keyword_set(bottom) then $
  start_PS, plotPath + strmid(config.saveFilename,0,strlen(config.saveFilename)-4)+'.eps', xs=8, ys=4
  
    !p.charsize = 0.8
    loadColorTable,'helix'

    ; density
    cgPlot, /nodata, xMinMax, yMinMax, pos=posLeft, xs=5, ys=5, /noerase
    tv, sphmap.dens_out,posLeft[0],posLeft[1],/normal,xsize=0.5
    
    ; circle at virial radius
    tvcircle,config.haloVirRad,0,0,cgColor('white'),thick=0.6,/data
    
    ; scale bar
    if ~keyword_set(bottom) then begin
      len = 200.0 ;kpc
      ;xpos = [config.boxCen[0]-config.boxSizeImg[0]*0.5,config.boxCen[0]-config.boxSizeImg[0]*0.5+len]
      ;ypos = replicate(config.boxCen[1]+config.boxSizeImg[1]*0.48,2)
      
      xpos = [xMinMax[0]*0.9,xMinMax[0]*0.9+len]
      ypos = replicate(yMinMax[1]*0.93,2)
      
      cgText,mean(xpos),ypos*0.9,string(len,format='(i3)')+' kpc',alignment=0.5,color=cgColor('white')
      loadct,0,/silent
      oplot,xpos,ypos,color=cgColor('white'),thick=!p.thick+0.5
      loadColorTable,'helix'
    endif
          
    if keyword_set(bottom) then $
      cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=1.0,color=cgColor('light gray'),/overplot
          
    ; mass weighted temperature
    ;loadColorTable,'blue-red'
    loadct,33,/silent
    cgPlot, /nodata, xMinMax, yMinMax, pos=posRight, xs=5, ys=5, /noerase
    tv, sphmap.quant_out,posRight[0],posRight[1],/normal,xsize=0.5

    ; circle at virial radius
    tvcircle,config.haloVirRad,0,0,cgColor('white'),thick=0.6,/data    
    
    ; redshift and halo mass
    if ~keyword_set(bottom) then begin
      cgText,0.99,0.98,"z = "+string(config.sP.redshift,format='(f3.1)'),alignment=1.0,$
             color=cgColor('white'),/normal
      cgText,0.99,0.96,"M = "+string(config.haloMass,format='(f4.1)'),alignment=1.0,$
             color=cgColor('white'),/normal
    endif
    
    ; simulation name
    if keyword_set(top) then cgText,0.5,0.96,"GADGET",charsize=!p.charsize+0.5,$
      alignment=0.5,/normal,color=cgColor('white')
    if keyword_set(bottom) then cgText,0.5,0.46,"AREPO",charsize=!p.charsize+0.5,$
      alignment=0.5,/normal,color=cgColor('white')    
            
    ; quantity name
    if keyword_set(bottom) then begin
      cgText,0.015,0.02,"density",charsize=!p.charsize+0.3,alignment=0.0,/normal,color=cgColor('white')
      cgText,0.985,0.02,"temperature",charsize=!p.charsize+0.3,alignment=1.0,/normal,color=cgColor('white')
    endif
      
    ; dividing line
    ;cgPlot,[xMinMax[0],xMinMax[0]],yMinMax,line=0,thick=1.0,color=cgColor('light gray'),/overplot
    
    if keyword_set(bottom) then $
      cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=1.0,color=cgColor('light gray'),/overplot
            
    ; colorbars
    !x.thick = 1.0
    !y.thick = 1.0
        
    if config.barType eq '1bar' then begin
      ; one colorbar (centered)
      pos = [0.02,0.35,0.076,0.65]
      fac = 1.0
      if keyword_set(bottom) then fac = 0.5
      pos *= [fac,1,fac,1]
      colorbar,position=pos,divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
      cgText,textoidl("\rho_{gas}")+","+textoidl("T_{gas}"),0.5,0.04*fac,alignment=0.5,/normal,color=cgColor('black')
    endif
            
  if ~keyword_set(top) and ~keyword_set(bottom) then $            
  end_PS, pngResize=60, /deletePS

end

; plotScatterAndSphmap(): plot side by side projection results from CalcSphMap

pro plotScatterAndSphmap, map=sphmap, scatter=scatter, config=config, left=left, right=right

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if ~keyword_set(left) and ~keyword_set(right) then message,'Error: Only for comparison.'
  
  if keyword_set(left) then leftOffset = 0.0
  if keyword_set(right) then leftOffset = 0.5

  ; UPPER TWO PANELS - scatter

    xMinMax = [-config.boxSizeScat[0]/2.0,config.boxSizeScat[0]/2.0]
    yMinMax = [-config.boxSizeScat[1]/2.0,config.boxSizeScat[1]/2.0]  
  
    ; plot position (normalized)
    ;posTop = [leftOffset,0.75,leftOffset+0.5,1.0] ;4x2
    ;posBottom = posTop - [0.0,0.25,0.0,0.25] ;4x2
    
    posTop = [leftOffset,0.666,leftOffset+0.5,1.0]
    posBottom = [leftOffset,0.333,leftOffset+0.5,0.666]
  
    !p.thick = 2.0
    !p.charsize = 0.8
  
    ; fill with black background
    if ~keyword_set(right) then $
      cgColorfill,[1,1,0,0,1],[1,0,0,1,1],/normal,color=cgColor('black')
  
    ; color table and establish temperature -> color mapping
    loadColorTable,'helix';,/reverse
    
    TVLCT, rr, gg, bb, /GET

    newcolors_left = getColor24([[rr[scatter.cinds_left]], [gg[scatter.cinds_left]], [bb[scatter.cinds_left]]])

    ; all gas (left panel)
    cgPlot, /nodata, xMinMax, yMinMax, pos=posTop, xs=5, ys=5, /noerase
    
    ; circle at virial radius
    tvcircle,config.haloVirRad,0,0,cgColor('light gray'),thick=0.8,/data
    
    ; particle loop for velocity vector plotting
    nCutoutLeft = n_elements(scatter.pos_left[0,*])
    for i=0L,nCutoutLeft-1 do $
      oplot,[scatter.pos_left[config.axes[0],i],scatter.pos2_left[config.axes[0],i]],$
             [scatter.pos_left[config.axes[1],i],scatter.pos2_left[config.axes[1],i]],$
             line=0,color=newcolors_left[i]
             
    newcolors_left = !NULL
               
    ; dividing lines
    cgPlot,xMinMax,[yMinMax[0],yMinMax[0]],line=0,thick=1.0,color=cgColor('light gray'),/overplot
    if keyword_set(right) then $
      cgPlot,[xMinMax[0],xMinMax[0]],yMinMax,line=0,thick=1.0,color=cgColor('light gray'),/overplot
               
    ; cold gas / dark matter (right panel)
    cgPlot, /nodata, xMinMax, yMinMax, pos=posBottom, xs=5, ys=5, /noerase

    tvcircle,config.haloVirRad,0,0,cgColor('light gray'),thick=0.8,/data
    
    newcolors_right = getColor24([[rr[scatter.cinds_right]], [gg[scatter.cinds_right]], [bb[scatter.cinds_right]]])
    
    ; particle loop for velocity vector plotting (cold gas only)
    nCutoutRight = n_elements(scatter.pos_right[0,*])
    for i=0L,nCutoutRight-1 do $
      oplot,[scatter.pos_right[config.axes[0],i],scatter.pos2_right[config.axes[0],i]],$
             [scatter.pos_right[config.axes[1],i],scatter.pos2_right[config.axes[1],i]],$
             line=0,color=newcolors_right[i]
             
    newcolors_right = !NULL  
  
    if ~keyword_set(right) then begin
      len = 100.0 ;kpc
      ;xpos = [config.boxCen[0]-config.boxSizeImg[0]*0.5,config.boxCen[0]-config.boxSizeImg[0]*0.5+len]
      ;ypos = replicate(config.boxCen[1]+config.boxSizeImg[1]*0.48,2)
      
      xpos = [xMinMax[0]*0.9,xMinMax[0]*0.9+len]
      ypos = replicate(yMinMax[1]*0.93,2)
      
      cgText,mean(xpos),ypos*0.9,string(len,format='(i3)')+' kpc',alignment=0.5,color=cgColor('white')
      loadct,0,/silent
      oplot,xpos,ypos,color=cgColor('white'),thick=!p.thick+0.5
      loadColorTable,'helix'
    endif
    
    ; dividing lines
    if keyword_set(right) then $
      cgPlot,[xMinMax[0],xMinMax[0]],yMinMax,line=0,thick=1.0,color=cgColor('light gray'),/overplot
    
    ; colorbars
    xthick = !x.thick
    ythick = !y.thick
    !x.thick = 1.0
    !y.thick = 1.0
    
    if config.barType eq '1bar' then begin
      labelText = "log T_{gas} [K]"
      
      ; one colorbar (centered)
      offset = 0.333
      pos = [offset+0.01,0.3,offset+0.034,0.7]
      colorbar,position=pos,divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
      
      ; colorbar labels
      ;cbLabels = str([fix(config.fieldMinMax[0]),fix(config.fieldMinMax[1])]
      cbLabels = str(string([fix(10*config.fieldMinMax[0])/10,fix(10*config.fieldMinMax[1])/10],format='(f3.1)'))
      rAdjust = 0.01*(strlen(cbLabels[1])-1)
      
      cgText,0.5,offset+0.018,textoidl(labelText),alignment=0.5,color=cgColor('black'),/normal
      cgText,0.324,offset+0.017,cbLabels[0],alignment=0.5,color=cgColor('black'),/normal
      cgText,0.692-rAdjust,offset+0.017,cbLabels[1],alignment=0.5,color=cgColor('black'),/normal
    endif
    
    !x.thick = xthick
    !y.thick = ythick

  
  ; LOWER TWO PANELS - sphmap
  
  clipFacs = [0.10,400.0,1.05,280.0]
  
  ; dens: color map and rescale
  w = where(sphmap.dens_out le 0.0,count,comp=wc)

  sphmap.dens_out = alog10(sphmap.dens_out)
  if count gt 0 then sphmap.dens_out[w] = 0.0
  sphmap.dens_out = sphmap.dens_out-clipFacs[0]
  
  w = where(sphmap.dens_out lt 0.0,count)
  if count gt 0 then sphmap.dens_out[w] = 0.0
  
  sphmap.dens_out = sphmap.dens_out*clipFacs[1] / (max(sphmap.dens_out)-min(sphmap.dens_out)) < 255 ; 0-255

  ; temp
  w = where(sphmap.quant_out eq 0.0,count,comp=wc)
  if count gt 0 then sphmap.quant_out[w] = min(sphmap.quant_out[wc])
  
  sphmap.quant_out = alog10(sphmap.quant_out)
  sphmap.quant_out = sphmap.quant_out-clipFacs[2]*min(sphmap.quant_out)
  
  w = where(sphmap.quant_out lt 0.0,count)
  if count gt 0 then sphmap.quant_out[w] = 0.0
  
  sphmap.quant_out = sphmap.quant_out*clipFacs[3] / (max(sphmap.quant_out)-min(sphmap.quant_out)) < 255 ;0-255

  ; box
  xMinMax = [config.boxCen[0]-config.boxSizeImg[0]/2.0,config.boxCen[0]+config.boxSizeImg[0]/2.0]
  yMinMax = [config.boxCen[1]-config.boxSizeImg[1]/2.0,config.boxCen[1]+config.boxSizeImg[1]/2.0]
  
  ; plot position (normalized)
  ;posTop = [leftOffset,0.25,leftOffset+0.5,0.5]
  ;posBottom = posTop - [0.0,0.25,0.0,0.25]
  
  posBottom = [leftOffset,0.0,leftOffset+0.5,0.333]
  
  !p.charsize = 0.8
  loadColorTable,'helix'

  ; density
  ;cgPlot, /nodata, xMinMax, yMinMax, pos=posTop, xs=5, ys=5, /noerase
  ;tv, sphmap.dens_out,posTop[0],posTop[1],/normal,xsize=0.5
  
  ; circle at virial radius
  ;tvcircle,config.haloVirRad,0,0,cgColor('white'),thick=0.8,/data
        
  ; horizontal dividing line
  ;cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=1.0,color=cgColor('light gray'),/overplot
        
  ; mass weighted temperature
  ;loadColorTable,'blue-red'
  loadct,33,/silent
  cgPlot, /nodata, xMinMax, yMinMax, pos=posBottom, xs=5, ys=5, /noerase
  tv, sphmap.quant_out,posBottom[0],posBottom[1],/normal,xsize=0.5

  ; circle at virial radius
  tvcircle,config.haloVirRad,0,0,cgColor('white'),thick=0.8,/data    
  
  ; scale bar
  if ~keyword_set(right) then begin
    len = 200.0 ;kpc
    ;xpos = [config.boxCen[0]-config.boxSizeImg[0]*0.5,config.boxCen[0]-config.boxSizeImg[0]*0.5+len]
    ;ypos = replicate(config.boxCen[1]+config.boxSizeImg[1]*0.48,2)
    
    xpos = [xMinMax[0]*0.9,xMinMax[0]*0.9+len]
    ypos = replicate(yMinMax[1]*0.93,2)
    
    cgText,mean(xpos),ypos*0.9,string(len,format='(i3)')+' kpc',alignment=0.5,color=cgColor('white')
    loadct,0,/silent
    oplot,xpos,ypos,color=cgColor('white'),thick=!p.thick+0.5
    loadColorTable,'helix'
  endif
  
  ; redshift and halo mass
  if ~keyword_set(left) then begin
    cgText,0.99,0.98,"z = "+string(config.sP.redshift,format='(f3.1)'),alignment=1.0,$
           color=cgColor('white'),/normal
    cgText,0.99,0.96,"M = "+string(config.haloMass,format='(f4.1)'),alignment=1.0,$
           color=cgColor('white'),/normal
  endif
  
  ; simulation name
  if keyword_set(left) then cgText,0.25,0.97,"GADGET",charsize=!p.charsize+0.5,$
    alignment=0.5,/normal,color=cgColor('white')
  if keyword_set(right) then cgText,0.75,0.97,"AREPO",charsize=!p.charsize+0.5,$
    alignment=0.5,/normal,color=cgColor('white')    
          
  ; quantity name
  if keyword_set(right) then begin
    cgText,0.99,0.676,"all gas",charsize=!p.charsize+0.3,alignment=1.0,/normal,color=cgColor('white')
    cgText,0.99,0.343,"cold gas",charsize=!p.charsize+0.3,alignment=1.0,/normal,color=cgColor('white')
    cgText,0.99,0.01,"projected temperature",charsize=!p.charsize+0.3,alignment=1.0,/normal,color=cgColor('white')
  endif
    
  ; dividing lines
  cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=1.0,color=cgColor('light gray'),/overplot
  if keyword_set(right) then $
    cgPlot,[xMinMax[0],xMinMax[0]],yMinMax,line=0,thick=1.0,color=cgColor('light gray'),/overplot
          
end

; sphMapHalos: run sph kernel density projection on gas particles/cells with boxes centered on halos

pro sphMapHalos;, sP=sP, gcIDs=gcIDs, coldOnly=coldOnly

  compile_opt idl2, hidden, strictarr, strictarrsubs

  sP = simParams(res=512,run='tracer',redshift=2.0)
  
  haloID = 304 ;z2.304 z2.301 z2.130 z2.64
  gcID = getMatchedIDs(sPa=sP,sPg=sP,haloID=haloID)
  gcIDs = [gcID.a]
  
  if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  sizeFac  = 5.0       ; times rvir
  nPixels  = [800,800] ; px
  axisPair = [0,1]
  
  ; loop over all non-background subhalos and image
  foreach gcID, gcIDs do begin
  
	  cutout = cosmoVisCutout(sP=sP,gcInd=gcID,sizeFac=sizeFac)
    
    saveFilename = 'map.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+$
                   '.axes'+str(axisPair[0])+str(axisPair[1])+'.sav'
               
    ; get box center (we have relative positions centered around origin)
    boxCenImg = [0,0,0]
    
    print,'['+string(gcID,format='(i4)')+'] Mapping ['+str(axisPair[0])+' '+$
          str(axisPair[1])+'] with '+str(cutout.boxSizeImg[0])+$
          ' kpc box around subhalo center ['+str(boxCenImg[0])+' '+str(boxCenImg[1])+' '+str(boxCenImg[2])+']'

    ; calculate projection using sph kernel (boxsize=0 for non-periodic)
    sphmap = calcSphMap(cutout.loc_pos,cutout.loc_hsml,cutout.loc_mass,cutout.loc_temp,$
                        boxSizeImg=cutout.boxSizeImg,boxSizeSim=0,boxCen=boxCenImg,$
                        nPixels=nPixels,axes=axisPair,ndims=3)
  
    ; plot
    config = {saveFilename:saveFilename,sizeFac:sizeFac,nPixels:nPixels,axes:axisPair,$
              gcID:gcID,haloMass:cutout.haloMass,haloVirRad:cutout.haloVirRad,$
              boxCen:boxCenImg,boxSizeImg:cutout.boxSizeImg,sP:sP,bartype:'1bar'}
  
    plotSphmapDensQuant, map=sphmap, config=config
  
  endforeach ;gcIDs
  stop
end

; sphMapHaloComp(): compare density/temperature of arepo/gadget in 2x2 panels

pro sphMapHaloComp

  compile_opt idl2, hidden, strictarr, strictarrsubs

  sPa = simParams(res=512,run='tracer',redshift=2.0)
  sPg = simParams(res=512,run='gadget',redshift=2.0)
  
  ; config
  sizeFac  = 5.0       ; times rvir
  nPixels  = [800,800] ; px
  axisPair = [0,1]
  hsmlFacs = 4*[1.0,2.0] ; arepo gadget, [1.0,2.0] highmass
  
  haloID = 130 ;z2.304 z2.301 z2.130 z2.64
  gcID = getMatchedIDs(sPa=sPa,sPg=sPg,haloID=haloID)
  
  saveFilename = 'map.'+sPa.savPrefix+str(sPa.res)+'.'+str(sPa.snap)+'.h'+str(gcID.a)+$
                 '.'+sPg.savPrefix+str(sPg.res)+'.h'+str(gcID.g)+'.axes'+str(axisPair[0])+str(axisPair[1])+'.eps'
                 
  start_PS, sPa.plotPath + saveFilename, xs=8, ys=8   
  
  ; GADGET
  cutout = cosmoVisCutout(sP=sPg,gcInd=gcID.g,sizeFac=sizeFac)

  ; enhance hsml and make boxsize smaller
  cutout.loc_hsml *= hsmlFacs[0]
  cutout.boxSizeImg *= 0.95
  
  ; get box center (we have relative positions centered around origin)
  boxCenImg = [0,0,0]
  
  print,'['+string(gcID.g,format='(i4)')+'] Mapping ['+str(axisPair[0])+' '+$
        str(axisPair[1])+'] with '+str(cutout.boxSizeImg[0])+$
        ' kpc box around subhalo center ['+str(cutout.boxCen[0])+' '+str(cutout.boxCen[1])+' '+str(cutout.boxCen[2])+']'

  ; calculate projection using sph kernel
  sphmap = calcSphMap(cutout.loc_pos,cutout.loc_hsml,cutout.loc_mass,cutout.loc_temp,$
                      boxSizeImg=cutout.boxSizeImg,boxSizeSim=0,boxCen=boxCenImg,$
                      nPixels=nPixels,axes=axisPair,ndims=3)

  ; plot
  config = {saveFilename:'',sizeFac:sizeFac,nPixels:nPixels,axes:axisPair,$
            gcID:gcID.g,haloMass:cutout.haloMass,haloVirRad:cutout.haloVirRad,$
            boxCen:boxCenImg,boxSizeImg:cutout.boxSizeImg,sP:sPg,bartype:'none'}

  plotSphmapDensQuant, map=sphmap, config=config, /top
  
  ; AREPO
  cutout = cosmoVisCutout(sP=sPa,gcInd=gcID.a,sizeFac=sizeFac)
  
  ; enhance hsml and make boxsize smaller
  cutout.loc_hsml *= hsmlFacs[1]
  cutout.boxSizeImg *= 0.95
 
  ; get box center (we have relative positions centered around origin)
  boxCenImg = [0,0,0]
  
  print,'['+string(gcID.a,format='(i4)')+'] Mapping ['+str(axisPair[0])+' '+$
        str(axisPair[1])+'] with '+str(cutout.boxSizeImg[0])+$
        ' kpc box around subhalo center ['+str(cutout.boxCen[0])+' '+str(cutout.boxCen[1])+' '+str(cutout.boxCen[2])+']'

  ; calculate projection using sph kernel
  sphmap = calcSphMap(cutout.loc_pos,cutout.loc_hsml,cutout.loc_mass,cutout.loc_temp,$
                      boxSizeImg=cutout.boxSizeImg,boxSizeSim=0,boxCen=boxCenImg,$
                      nPixels=nPixels,axes=axisPair,ndims=3)

  ; plot
  config = {saveFilename:'',sizeFac:sizeFac,nPixels:nPixels,axes:axisPair,$
            gcID:gcID.a,haloMass:cutout.haloMass,haloVirRad:cutout.haloVirRad,$
            boxCen:boxCenImg,boxSizeImg:cutout.boxSizeImg,sP:sPa,bartype:'none'}

  plotSphmapDensQuant, map=sphmap, config=config, /bottom
  
  end_PS, pngResize=60;, /deletePS

end

; sphScatterAndMapHaloComp(): compare scatterplots and density/temperature of arepo/gadget in 4x2 panels

pro sphScatterAndMapHaloComp

  compile_opt idl2, hidden, strictarr, strictarrsubs

  sPa = simParams(res=512,run='tracer',redshift=2.0)
  sPg = simParams(res=512,run='gadget',redshift=2.0)
  
  ; config
  sizeFacMap  = 5.0       ; times rvir
  sizeFacScat = 3.5       ; times rvir
  nPixels  = [800,800]    ; px
  axisPair = [0,1]        ; xy
  hsmlFacs = 4*[1.0,2.0]  ; arepo gadget, [1.0,2.0] highmass
  tempCut  = 5.0           ; log K
  fieldMinMax  = [4.0,7.0] ; log K
  
  haloID = 304 ;z2.304 z2.301 z2.130 z2.64
  gcID = getMatchedIDs(sPa=sPa,sPg=sPg,haloID=haloID)
  
  saveFilename = 'sc.map.'+sPa.savPrefix+str(sPa.res)+'.'+str(sPa.snap)+'.h'+str(gcID.a)+$
                 '.'+sPg.savPrefix+str(sPg.res)+'.h'+str(gcID.g)+'.axes'+str(axisPair[0])+str(axisPair[1])+'.eps'
                 
  start_PS, sPa.plotPath + saveFilename, xs=6, ys=9  
  
  ; GADGET
  cutout = cosmoVisCutout(sP=sPg,gcInd=gcID.g,sizeFac=sizeFacMap)

  ; enhance hsml and make boxsize smaller
  cutout.loc_hsml *= hsmlFacs[0]
  cutout.boxSizeImg *= 0.95
  
  ; get box center (we have relative positions centered around origin)
  boxCenImg = [0,0,0]
  
  print,'['+string(gcID.g,format='(i4)')+'] Mapping ['+str(axisPair[0])+' '+$
        str(axisPair[1])+'] with '+str(cutout.boxSizeImg[0])+$
        ' kpc box around subhalo center ['+str(cutout.boxCen[0])+' '+str(cutout.boxCen[1])+' '+str(cutout.boxCen[2])+']'

  ; calculate projection using sph kernel
  sphmap = calcSphMap(cutout.loc_pos,cutout.loc_hsml,cutout.loc_mass,cutout.loc_temp,$
                      boxSizeImg=cutout.boxSizeImg,boxSizeSim=0,boxCen=boxCenImg,$
                      nPixels=nPixels,axes=axisPair,ndims=3)
                      
  cutout2 = cosmoVisCutout(sP=sPg,gcInd=gcID.g,sizeFac=sizeFacScat)

  ; make scatter points  
  colorinds = (cutout2.loc_temp-fieldMinMax[0])*205.0 / (fieldMinMax[1]-fieldMinMax[0]) ;0-205
  colorinds = fix(colorinds + 50.0) > 0 < 255 ;50-255  

  ; second/right panel cutout
  wSecond = where(cutout2.loc_temp le tempCut,nCutoutSecond,comp=wComp)
   
  loc_pos_second   = cutout2.loc_pos[*,wSecond]
  loc_pos2_second  = cutout2.loc_pos2[*,wSecond]
  colorinds_second = colorinds[wSecond]
  
  scatter = {pos_left:cutout2.loc_pos,pos2_left:cutout2.loc_pos2,pos_right:loc_pos_second,$
             pos2_right:loc_pos2_second,cinds_left:colorinds,cinds_right:colorinds_second}  
             
  ; plot
  config = {saveFilename:'',nPixels:nPixels,axes:axisPair,fieldMinMax:fieldMinMax,$
            gcID:gcID.g,haloMass:cutout.haloMass,haloVirRad:cutout.haloVirRad,$
            boxCen:boxCenImg,boxSizeImg:cutout.boxSizeImg,boxSizeScat:cutout2.boxSizeImg,sP:sPg,bartype:'none'}
            
  plotScatterAndSphmap, map=sphmap, scatter=scatter, config=config, /left
  
  ; AREPO
  cutout = cosmoVisCutout(sP=sPa,gcInd=gcID.a,sizeFac=sizeFacMap)
  
  ; enhance hsml and make boxsize smaller
  cutout.loc_hsml *= hsmlFacs[1]
  cutout.boxSizeImg *= 0.95
 
  ; get box center (we have relative positions centered around origin)
  boxCenImg = [0,0,0]
  
  print,'['+string(gcID.a,format='(i4)')+'] Mapping ['+str(axisPair[0])+' '+$
        str(axisPair[1])+'] with '+str(cutout.boxSizeImg[0])+$
        ' kpc box around subhalo center ['+str(cutout.boxCen[0])+' '+str(cutout.boxCen[1])+' '+str(cutout.boxCen[2])+']'

  ; calculate projection using sph kernel
  sphmap = calcSphMap(cutout.loc_pos,cutout.loc_hsml,cutout.loc_mass,cutout.loc_temp,$
                      boxSizeImg=cutout.boxSizeImg,boxSizeSim=0,boxCen=boxCenImg,$
                      nPixels=nPixels,axes=axisPair,ndims=3)

  cutout2 = cosmoVisCutout(sP=sPa,gcInd=gcID.a,sizeFac=sizeFacScat)
         
  ; make scatter points  
  colorinds = (cutout2.loc_temp-fieldMinMax[0])*205.0 / (fieldMinMax[1]-fieldMinMax[0]) ;0-205
  colorinds = fix(colorinds + 50.0) > 0 < 255 ;50-255  

  ; second/right panel cutout
  wSecond = where(cutout2.loc_temp le tempCut,nCutoutSecond,comp=wComp)
   
  loc_pos_second   = cutout2.loc_pos[*,wSecond]
  loc_pos2_second  = cutout2.loc_pos2[*,wSecond]
  colorinds_second = colorinds[wSecond]
  
  scatter = {pos_left:cutout2.loc_pos,pos2_left:cutout2.loc_pos2,pos_right:loc_pos_second,$
             pos2_right:loc_pos2_second,cinds_left:colorinds,cinds_right:colorinds_second}  
       
  ; plot
  config = {saveFilename:'',nPixels:nPixels,axes:axisPair,fieldMinMax:fieldMinMax,$
            gcID:gcID.a,haloMass:cutout.haloMass,haloVirRad:cutout.haloVirRad,$
            boxCen:boxCenImg,boxSizeImg:cutout.boxSizeImg,boxSizeScat:cutout2.boxSizeImg,sP:sPa,bartype:'1bar'}
       
  plotScatterAndSphmap, map=sphmap, scatter=scatter, config=config, /right
  
  end_PS, pngResize=60;, /deletePS

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
  
      saveFilename = 'dmMap.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+$
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
