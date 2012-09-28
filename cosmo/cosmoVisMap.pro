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

pro plotSphmapDensQuant, map=sphmap, config=config

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; color map and rescale
  w = where(sphmap.dens_out eq 0.0,count,comp=wc)
  if count gt 0 then sphmap.dens_out[w] = min(sphmap.dens_out[wc])

  sphmap.dens_out = alog10(sphmap.dens_out)
  sphmap.dens_out = (sphmap.dens_out-min(sphmap.dens_out))*254.0 / $
                    (max(sphmap.dens_out)-min(sphmap.dens_out)) ;0-254
  sphmap.dens_out += 1.0 ;1-255
  
  w = where(sphmap.quant_out eq 0.0,count,comp=wc)
  if count gt 0 then sphmap.quant_out[w] = min(sphmap.quant_out[wc])
  
  sphmap.quant_out = alog10(sphmap.quant_out)
  sphmap.quant_out = (sphmap.quant_out-min(sphmap.quant_out))*254.0 / $
                    (max(sphmap.quant_out)-min(sphmap.quant_out)) ;0-254
  sphmap.quant_out += 1.0 ;1-255

  xMinMax = [config.boxCen[0]-config.boxSizeImg[0]/2.0,config.boxCen[0]+config.boxSizeImg[0]/2.0]
  yMinMax = [config.boxCen[1]-config.boxSizeImg[1]/2.0,config.boxCen[1]+config.boxSizeImg[1]/2.0]
  
  plotPath = '/n/home07/dnelson/data3/HaloComp/'
  
  start_PS, plotPath + strmid(config.saveFilename,0,strlen(config.saveFilename)-4)+'.eps', xs=8, ys=4
  
    !p.charsize = 0.8
    loadColorTable,'helix'
    
    ; density
    cgPlot, /nodata, xMinMax, yMinMax, pos=[0.0,0.0,0.5,1.0], xs=5, ys=5
    tv, sphmap.dens_out,0.0,0.0,/normal,xsize=0.5
    
    ; circle at virial radius
    tvcircle,config.haloVirRad,0,0,cgColor('dark gray'),thick=0.6,/data
    
    ; scale bar
 ;   len = 100.0 ;ckpc
 ;   cgText,mean([config.boxCen[0]-config.boxSizeImg[0]/2.2,$
 ;                config.boxCen[0]-config.boxSizeImg[0]/2.2+len]),$
 ;          config.boxCen[1]+config.boxSizeImg[0]/2.3,$
 ;          string(len,format='(i3)')+' ckpc',alignment=0.5,color=cgColor('black')
 ;   cgPlot,[config.boxCen[0]-config.boxSizeImg[0]/2.2,config.boxCen[0]-config.boxSizeImg[0]/2.2+len],$
 ;          [config.boxCen[1]+config.boxSizeImg[1]/2.1,config.boxCen[1]+config.boxSizeImg[1]/2.1],$
 ;          color=cgColor('black'),/overplot
           
    ; mass weighted temperature
    cgPlot, /nodata, xMinMax, yMinMax, pos=[0.0,0.0,0.5,1.0], xs=5, ys=5, /noerase
    tv, sphmap.quant_out,0.5,0.0,/normal,xsize=0.5
    
    ; redshift and halo mass
 ;   cgText,0.99,0.96,"z = "+string(config.sP.redshift,format='(f3.1)'),alignment=1.0,$
 ;          color=cgColor('black'),/normal
 ;   cgText,0.99,0.92,"M = "+string(config.haloMass,format='(f4.1)'),alignment=1.0,$
 ;          color=cgColor('black'),/normal
             
  end_PS, pngResize=60;, /deletePS

end

; sphMapHalos: run sph kernel density projection on gas particles/cells with boxes centered on halos

pro sphMapHalos, sP=sP, gcIDs=gcIDs, coldOnly=coldOnly

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  sizeFac = 3.5       ; times rvir
  nPixels = [800,800] ; px
  
  axes = list([0,1],[0,2],[1,2]) ;xy,xz,yz
  
  if keyword_set(coldOnly) then coldFlag = '.cold' else coldFlag = ''
  
  ; target list
  gc    = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sP) 

  ; load u,nelec and calculate temperature
  u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
  nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
  temp  = convertUtoTemp(u,nelec)
  u     = !NULL
  nelec = !NULL

  ; load gas positions and masses from snapshot
  pos  = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  mass = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
  
  ; load HSMLs or volumes (convert to sizes)
  if sP.trMCPerCell eq 0 then begin
    hsml = loadSnapshotSubset(sP=sP,partType='gas',field='hsml')
    hsml = 1.0 * temporary(hsml); increase hsml to decrease visualization noise
  endif else begin
    hsml = loadSnapshotSubset(sP=sP,partType='gas',field='vol')
    hsml = (temporary(hsml) * 3.0 / (4*!pi))^(1.0/3.0) ;cellrad [ckpc]
    hsml = 1.75 * temporary(hsml) ; increase hsml to decrease visualization noise
  endelse
  
  ; cold only
  if keyword_set(coldOnly) then begin
    wCold = where(alog10(temp) le 5.0,count)
    print,' Cold cut have remaining ['+str(count)+'] of ['+str(n_elements(temp))+'] particles.'
    
    temp = temp[wCold]
    mass = mass[wCold]
    hsml = hsml[wCold]
    pos  = pos[*,wCold]
  endif

  ; loop over all non-background subhalos and image
  foreach gcID, gcIDs do begin
  
    ; get subhalo position and size of imaging box
    haloVirRad = gc.group_r_crit200[gc.subgroupGrNr[gcID]] ;ckpc
    haloMass   = codeMassToLogMsun(gc.subgroupMass[gcID])
    
    boxCen     = sgcen[*,gcID]
    boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcID]])
    boxSizeImg = [boxSize,boxSize,boxSize] ; cube
  
    foreach axisPair, axes do begin
  
      saveFilename = 'map.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+$
                     '.axes'+str(axisPair[0])+str(axisPair[1])+coldFlag+'.sav'
                    
      if ~file_test(sP.derivPath+'sphMaps/'+saveFilename) then begin
      
        ; get box center (in terms of specified axes)
        boxCenImg  = [sgcen[axisPair[0],gcID],sgcen[axisPair[1],gcID],sgcen[3-axisPair[0]-axisPair[1],gcID]]
        
        print,'['+string(gcID,format='(i4)')+'] Mapping ['+str(axisPair[0])+' '+$
              str(axisPair[1])+'] with '+str(boxSize[0])+$
              ' kpc box around subhalo center ['+str(boxCen[0])+' '+str(boxCen[1])+' '+str(boxCen[2])+']'

        ; calculate projection using sph kernel
        sphmap = calcSphMap(pos,hsml,mass,temp,$
                            boxSizeImg=boxSizeImg,boxSizeSim=sP.boxSize,boxCen=boxCenImg,$
                            nPixels=nPixels,axes=axisPair,ndims=3)
      
        ; save
        config = {saveFilename:saveFilename,sizeFac:sizeFac,nPixels:nPixels,axes:axes,boxCen:boxCen,$
                  gcID:gcID,haloMass:haloMass,haloVirRad:haloVirRad,$
                  boxCenImg:boxCenImg,boxSize:boxSize,boxSizeImg:boxSizeImg,sP:sP}
        ;save,sphmap,config,filename=sP.derivPath+'sphMaps/'+saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
      endif else begin
        restore,sP.derivPath +'sphMaps/'+ saveFilename
      endelse
    
      ; plot
      plotSphmapDensQuant, map=sphmap, config=config

    endforeach ;axisPair
  
  endforeach ;gcIDs
  
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
