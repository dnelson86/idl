; cosmoVisStars.pro
; cosmological boxes - 2d visualization of the stars / actual galaxies
; dnelson oct.2013

; plotScatterAndSphmap(): plot side by side projection results from CalcSphMap

pro plotScatterAndSphmap, map=sphmap, scatter=scatter, config=config, $
                          left=left, right=right, one=one, two=two, three=three

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if ~keyword_set(left) and ~keyword_set(right) and ~keyword_set(one) and $
     ~keyword_set(two) and ~keyword_set(three) then message,'error'
  
  ; 2 col, 3 rows
  if keyword_set(left) then leftOffset = 0.0
  if keyword_set(right) then leftOffset = 0.5
  if keyword_set(left) or keyword_set(right) then pWidth = 0.5
  
  ; 3 col, 3 rows
  if keyword_set(one) then leftOffset = 0.0
  if keyword_set(two) then leftOffset = 0.333
  if keyword_set(three) then leftOffset = 0.666
  if keyword_set(one) or keyword_set(two) or keyword_set(three) then pWidth = 0.333

  ; UPPER TWO PANELS - scatter

    xMinMax = [-config.boxSizeScat[0]/2.0,config.boxSizeScat[0]/2.0]
    yMinMax = [-config.boxSizeScat[1]/2.0,config.boxSizeScat[1]/2.0]  
  
    ; plot position (normalized)
    posTop = [leftOffset,0.666,leftOffset+pWidth,1.0]
    posBottom = [leftOffset,0.333,leftOffset+pWidth,0.666]
  
    !p.thick = 2.0
    !p.charsize = 0.8
  
    ; fill with black background
    if ~keyword_set(right) and ~keyword_set(two) and ~keyword_set(three) then $
      cgColorfill,[1,1,0,0,1],[1,0,0,1,1],/normal,color=cgColor('black')
  
    ; color table and establish temperature -> color mapping
    loadColorTable,config.ctNameScat

    ; all gas (left panel)
    cgPlot, /nodata, xMinMax, yMinMax, pos=posTop, xs=5, ys=5, /noerase
    
    ; circle at virial radius
    tvcircle,config.haloVirRad,0,0,cgColor('light gray'),thick=0.8,/data
    
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

    tvcircle,config.haloVirRad,0,0,cgColor('light gray'),thick=0.8,/data
        
    ; particle loop for velocity vector plotting (cold gas only)
    nCutoutRight = n_elements(scatter.pos_right[0,*])
    for i=0L,nCutoutRight-1,config.interval do $
      oplot,[scatter.pos_right[config.axes[0],i],scatter.pos_right2[config.axes[0],i]],$
             [scatter.pos_right[config.axes[1],i],scatter.pos_right2[config.axes[1],i]],$
             line=0,color=scatter.cinds_right[i],thick=config.lineThick
  
    if ~keyword_set(right) and ~keyword_set(two) and ~keyword_set(three) then begin
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
      offset = 0.333
      pos = [offset+0.01,0.3,offset+0.034,0.7]
      loadColorTable,config.ctNameScat
      cgColorbar,position=pos,divisions=0,charsize=0.000001,bottom=50,ticklen=0.00001
      
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
  tvcircle,config.haloVirRad,0,0,cgColor('white'),thick=0.8,/data    
  
  ; scale bar
  if ~keyword_set(right) and ~keyword_set(two) and ~keyword_set(three) then begin
    xpos = [xMinMax[0]*0.9,xMinMax[0]*0.9+config.scaleBarLen]
    ypos = replicate(yMinMax[1]*0.93,2)
    
    cgText,mean(xpos),ypos*0.9,string(config.scaleBarLen,format='(i3)')+' kpc',$
      alignment=0.5,color=cgColor('white')
    loadct,0,/silent
    oplot,xpos,ypos,color=cgColor('white'),thick=!p.thick+0.5
    loadColorTable,'bw linear'
  endif
  
  ; redshift and halo mass
  if ~keyword_set(left) and ~keyword_set(one) and ~keyword_set(two) then begin
    cgText,0.99,0.666-0.02,"z = "+string(config.sP.redshift,format='(f3.1)'),alignment=1.0,$
           color=cgColor('white'),/normal
    if config.haloMass ne -1 then $
      cgText,0.99,0.666-0.04,"M = "+string(config.haloMass,format='(f4.1)'),alignment=1.0,$
             color=cgColor('white'),/normal
  endif
  
  ; simulation name
  cgText,mean(posTop[[0,2]]),0.97,config.sP.simName,charsize=!p.charsize+0.5,$
    alignment=0.5,/normal,color=cgColor('white')
          
  ; quantity name
  if keyword_set(right) or keyword_set(three) then begin
    cgText,0.99,0.676,"all gas",charsize=!p.charsize+0.3,alignment=1.0,/normal,color=cgColor('white')
    cgText,0.99,0.343,config.secondText,charsize=!p.charsize+0.3,alignment=1.0,/normal,color=cgColor('white')
    cgText,0.99,0.03,"projected",charsize=!p.charsize+0.3,alignment=1.0,/normal,color=cgColor('white')
    cgText,0.99,0.01,config.colorField,charsize=!p.charsize+0.3,alignment=1.0,/normal,color=cgColor('white')
  endif
    
  ; dividing lines
  cgPlot,xMinMax,[yMinMax[1],yMinMax[1]],line=0,thick=1.0,color=cgColor('dark gray'),/overplot
  cgPlot,[xMinMax[0],xMinMax[0]],yMinMax,line=0,thick=1.0,color=cgColor('dark gray'),/overplot
          
end

; visStars()

pro visStars

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  sP = simParams(res=9,run='zoom_20mpc',redshift=2.0,hInd=0)
  
  haloID = 0 ;zoom.0 z2.304 z2.301 z2.130 z2.64
  
  nPixels  = [800,800]    ; px
  boxSize  = 20          ; ckpc
  scaleBarLen = 2       ; ckpc
  maxHSMLCutoff = 1.0    ; ckpc
  axes     = [0,1]        ; xy
  hsmlFac  = 2            ; times each cell hsml for sph projections
  nbottom  = 50           ; lift minimum scatterplot color from zero (to distinguish from black)

  ; load halo position
  gcID = getMatchedIDs(haloID=haloID,simParams=sP)
  gc = loadGroupCat(sP=sP,/skipIDs)
  subgroupCen = subgroupPosByMostBoundID(sP=sP)
  
  ; load star properties
  pos  = loadSnapshotSubset(sP=sP,partType='stars',field='pos')
  sft  = loadSnapshotSubset(sP=sP,partType='stars',field='stellarformationtime')
  mass = loadSnapshotSubset(sP=sP,partType='stars',field='mass')
  
  ; (1) restrict stars based on subhalo
  if 0 then begin
    halo_starIDs = gcPIDList(gc=gc, valGCids=[gcID], partType='stars')
  
    ; match to snapshot stars
    ids = loadSnapshotSubset(sP=sP,partType='stars',field='ids')
  
    idIndexMap = getIDIndexMap(ids,minid=minid)
    halo_starIDs_ind = idIndexMap[halo_starIDs-minid]
        
    idIndexMap   = !NULL
    halo_starIDs = !NULL
  
    pos  = pos[*,halo_starIDs_ind]
    sft  = sft[halo_starIDs_ind]
    mass = mass[halo_starIDs_ind]
  endif
  
  ; (2) restrict stars based on spatial box
  dists = periodicDists( reform(subgroupCen[*,gcID]), pos, sP=sP )
  
  inds = where(dists le boxSize)
  
  pos  = pos[*,inds]
  sft  = sft[inds]
  mass = mass[inds]
  
  ; calculate hsml
  hsml = calcHSML(pos, ndims=3, nNGB=32, boxSize=0)
  
  ; make positions relative to halo center
  for i=0,2 do pos[i,*] -= subgroupCen[i,gcID]
  
  ; transform age of stars into SDSS type colors
  ; TODO
  
  ; enforce maximum hsml by setting mass to zero
  w = where(hsml ge maxHSMLCutoff,count)
  if count gt 0 then mass[w] = 0.0
  
  ; do sph map on origin centered positions, scale
  sphmap = calcSphMap(pos,hsml,mass,sft,boxSizeImg=replicate(boxSize,3),boxSizeSim=0,$
                      boxCen=[0,0,0],nPixels=nPixels,axes=axes,ndims=3)
                 
  sphmap.dens_out = alog10( sphmap.dens_out )
  ;mapMinMax = minmax(sphmap.dens_out)
  print,minmax(sphmap.dens_out)
  mapMinMax = [-1.0,0.2] ;[-0.5,1.0]
                 
  sphmap.dens_out = (sphmap.dens_out-mapMinMax[0])*(255.0) / (mapMinMax[1]-mapMinMax[0])
  sphmap.dens_out = fix(sphmap.dens_out) > 0 < 255 ; 0-255    
      
  ; plot setup
  xMinMax = [-boxSize/2.0,boxSize/2.0]
  yMinMax = [-boxSize/2.0,boxSize/2.0]  
    
  ; plot
  start_PS, sP.plotPath + 'visStars_scat.eps', xs=8.0, ys=8.0
  
    cgColorfill,[1,1,0,0,1],[1,0,0,1,1],/normal,color=cgColor('black')
    loadct,0,/silent

    cgPlot, [0], [0], /nodata, xrange=xMinMax, yrange=yMinMax, pos=[0,0,1,1], /xs, /ys, /noerase
    
    ; scatter
    cgPlot,pos[axes[0],*],pos[axes[1],*],psym=3,color=cgColor('red'),/overplot
    
    ; scale bar
    xpos = [xMinMax[0]*0.96,xMinMax[0]*0.96+scaleBarLen]
    ypos = replicate(yMinMax[1]*0.96,2)
    
    cgText,mean(xpos),ypos*0.94,string(scaleBarLen,format='(i3)')+' kpc',$
      alignment=0.5,color=cgColor('white')
    ;loadct,0,/silent
    oplot,xpos,ypos,color=cgColor('white'),thick=!p.thick+0.5
    
  end_PS, pngResize=50, /deletePS
                 
  start_PS, sP.plotPath + 'visStars_map.eps', xs=8.0, ys=8.0 
    
    cgColorfill,[1,1,0,0,1],[1,0,0,1,1],/normal,color=cgColor('black')
    loadct,0,/silent

    cgPlot, [0], [0], /nodata, xrange=xMinMax, yrange=yMinMax, pos=[0,0,1,1], /xs, /ys, /noerase
    
    ; sph
    loadColorTable, 'blue-red2'
    tv, sphmap.dens_out,0.0,0.0,/normal;,xsize=pWidth ; mass-weighted quantity

    ; scale bar
    xpos = [xMinMax[0]*0.96,xMinMax[0]*0.96+scaleBarLen]
    ypos = replicate(yMinMax[1]*0.96,2)
    
    cgText,mean(xpos),ypos*0.94,string(scaleBarLen,format='(i3)')+' kpc',$
      alignment=0.5,color=cgColor('white')
    ;loadct,0,/silent
    oplot,xpos,ypos,color=cgColor('white'),thick=!p.thick+0.5
    
  end_PS, pngResize=50, /deletePS
  
  stop

end
