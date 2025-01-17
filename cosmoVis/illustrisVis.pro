; illustrisVis.pro
; illustris 1820^3 specialized visualization
; dnelson nov.2013

; findCuboidRemapInds

function findCuboidRemapInds, remapRatio=remapRatio, newBoxSize=newBoxSize, nPixels=nPixels
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; validate remapRatio
  if n_elements(remapRatio) ne 3 then message,'Error'
  remapProduct = remapRatio[0] * remapRatio[1] * remapRatio[2]
  if abs(1.0 - remapProduct) gt 1e-3 then message,'Error: Check L1*L2*L3=1 constraint.'
  
  ; load pre-computed mapping possibilities
  file = '/n/home07/dnelson/idl/CalcBoxRemap/mappings_N7.txt'
  headerLines = 2
  ptStruct = { e1:0.0, e2:0.0, e3:0.0, u11:0,u12:0,u13:0, u21:0,u22:0,u23:0, u31:0,u32:0,u33:0, periodicity:'' }
  
  res = loadCSV(headerLines, file, ptStruct)
  
  ; calculate closest matching edge length set (use abs(xyz) distance metric)
  dists = fltarr(n_elements(res))
  for i=0,n_elements(res)-1 do $
    dists[i] = abs(res[i].e1 - remapRatio[0]) + abs(res[i].e2 - remapRatio[1]) + abs(res[i].e3 - remapRatio[2])
  
  ind = closest(dists, 0.0)
  
  print,'remapRatio: ',remapRatio
  print,'remapProduct: ',remapProduct
  print,'Found distance: ',dists[ind],' with ind',ind
  help,res[ind]
  if dists[ind] ge 0.1 then print,'Warning: Inaccurate remapping matrix chosen.'
  
  remapMatrix = [ res[ind].u11, res[ind].u12, res[ind].u13, $
                  res[ind].u21, res[ind].u22, res[ind].u23 ,$
                  res[ind].u31, res[ind].u32, res[ind].u33 ]
                  
  ; set return
  newBoxSize = [ res[ind].e1, res[ind].e2, res[ind].e3 ]
  
  ; (option 1) calculate new pixel dimensions (enforce aspect ratio of transformation by 
  ; keeping requested width, changing height)
  ;nPixels[1] = round( nPixels[0] * (newBoxSize[1]/newBoxSize[0])  )

  ; (option 2) adjust box size, increase smaller xy dimension to match pixel aspect ratio
  ; and so keep the output image size exactly as requested (black border somewhere)
  print,'ADJUSTING BOXSIZE, KEEPING PIXEL DIMENSIONS CONSTANT'
  if newBoxSize[0] lt newBoxSize[1] then begin
    newBoxSize[0] = newBoxSize[1] * (float(nPixels[0])/nPixels[1])
  endif else begin
    newBoxSize[1] = newBoxSize[0] * (float(nPixels[1])/nPixels[0])
  endelse
                  
  return, remapMatrix

end

; illustrisArepoProj(): handle Arepo-based Voronoi projection map making
;                     for (1) halo-scale cutouts , (2) whole box slices, (3) whole box remaps

function illustrisArepoProj, sP=sP,haloID=haloID,sizeFac=sizeFac,sliceWidth=sliceWidth,$
                             nPixels=nPixels,axes=axes,quantName=quantName
  compile_opt idl2, hidden, strictarr, strictarrsubs

  if keyword_set(wholeBoxSlice) then message,'Error: Whole box slice with voronoi tracing not implemented.'
  if keyword_set(wholeBoxRemap) then message,'Error: Whole box remap with voronoi tracing not implemented.'
  if sliceWidth gt sizeFac then print,'Warning: sliceWidth > sizeFac'
    
  ; config
  if sP.run eq 'scylla' then begin
    filePath = '/n/home07/dnelson/ArepoVTK/scylla/output/' ; proj output dir
  endif else begin
    filePath = '/n/home07/dnelson/ArepoVTK/illustris.fof0/output/'
  endelse
  
  padding    = 40 ; ckpc, to avoid edge effects (~1% of box sidelength)
  quantities = ['density','entropy','metal','temp','velmag',$
                'xray','sz_y','velocity'] ; all have Npx*Npx elems, except velocity with 3*Npx*Npx
  
  if sP.run eq 'scylla' then zStr = '0' else zStr = ''
  
  ; locate subgroup index
  addDMToSave = 0
  gcInd = getMatchedIDs(simParams=sP,haloID=haloID)
    
  fileName = 'vorMap.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.h' + str(gcInd) + '.sf' + str(fix(sizeFac*10)) + '.' + $ 
    str(fix(sliceWidth*10)) + $
    '.px' + str(nPixels) + '.axes' + str(axes[0]) + str(axes[1]) + '_all'  
    
  saveFilename = sP.derivPath + 'binnedVals/' + fileName + '.sav'  
    
  ; if save already exists, immediate load and return requested map
  if file_test(saveFileName) then begin
    ; check existence of temporary files?
    tempFileName = filePath + 'dm_annih_proj_' + str(sP.snap) + '.dat'
    if file_test(tempFileName) then begin
      addDMToSave = 1
      print,'Verify intent to add DM to save, .c to proceed.' & stop
    endif else begin
      restore, saveFilename
      return, vorMaps
    endelse
  endif
    
  ; load metadata
  metaFilePath = sP.derivPath + 'cutouts/cutout.' + sP.savPrefix + str(sP.res) + '.' + $
                 str(sP.snap) + '.h' + str(gcInd) + '.sf' + str(fix(sizeFac*10)) + '_meta.sav'
                 
  restore,metaFilePath
    
  ; file doesn't exist, see if we have .dat files in path to make it (CAREFUL)
  tempFileName = filePath + quantities[0] + '_proj_' + zStr + str(sP.snap) + '.dat'
  
  if ~file_test(tempFileName) then begin
    ; no .dat files exist, need to invoke Arepo to make them 
    
    ; determine xyz bounds
    aMinMax = [padding, 2*boxSize - padding]
    bMinMax = [padding, 2*boxSize - padding]
    cMinMax = [] ; TODO
    message,'TODO'

  endif
  
  ; .dat files exist, combine them by looping over each quantity/file
  if addDMToSave eq 1 then begin
    ; load existing vorMaps structure, and modify quantities
    restore,saveFilename,/verbose
    quantities = ['dm_annih','density','velocity']
  endif else begin
    ; start new vorMaps structure
    vorMaps = {fileName:fileName,boxSize:boxSize,haloRvir:haloRvir}
  endelse
  
  foreach quantName,quantities,k do begin
    fileName = filePath + quantName + '_proj_' + zStr + str(sP.snap) + '.dat'
    if ~file_test(fileName) then message,'Error: Missing output file.'
    
    ; load
    openr,lun,fileName,/get_lun
      ;read header
      nPixelsX = 0L
      nPixelsY = 0L
      readu,lun,nPixelsX,nPixelsY

      ; read field
      if quantName ne 'velocity' then quant = fltarr(nPixelsX, nPixelsY)
      if quantName eq 'velocity' then quant = fltarr(3, nPixelsX, nPixelsY)
      readu,lun,quant
    close,lun
    free_lun,lun
  
    ; units/preprocessing based on quantity
    if quantName eq 'density' or quantName eq 'entropy' or $
       quantName eq 'xray'    or quantName eq 'sz_y'    or $
       quantName eq 'temp'    or quantName eq 'metal'   or $
       quantName eq 'dm_annih' then quant = alog10( quant )
  
    ; add to save structure
    if k eq 0 then begin
      vorMaps = mod_struct( vorMaps, 'nPixelsX', nPixelsX )
      vorMaps = mod_struct( vorMaps, 'nPixelsY', nPixelsY )
    endif
    
    ; only density/"velocity" for DM has special behavior: split first two components as vmag/vdisp
    if addDMToSave eq 1 then begin
    
      if quantName eq 'density' then quantName = 'dm_density'
      
      if quantName eq 'velocity' then begin
        vorMaps = mod_struct( vorMaps, 'dm_vmag',  reform(quant[0,*,*]) )
        vorMaps = mod_struct( vorMaps, 'dm_vdisp', alog10(reform(quant[1,*,*])) )
      endif else begin
        ; density (prepended "dm") or dm_annih
        vorMaps = mod_struct( vorMaps, quantName, quant )
      endelse
      
    endif else begin
      ; normal (gas)
      vorMaps = mod_struct( vorMaps, quantName, quant )
    endelse
  endforeach
    
  ; save and return
  save,vorMaps,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  ; remove .dat files
  foreach quantName,quantities do begin
    fileName = filePath + quantName + '_proj_' + zStr + str(sP.snap) + '.dat'
    cmd = 'rm '+fileName
    print,cmd
    spawn,cmd
  endforeach
  
  return, vorMaps

end

; illustrisMakeMap(): handle CalcSphMap-based map making
;                     for (1) halo-scale cutouts , (2) whole box slices, (3) whole box remaps
                           
function illustrisMakeMap, sP=sP,haloID=haloID,sizeFac=sizeFac,sliceWidth=sliceWidth,remapRatio=remapRatio,$
                           hsmlFac=hsmlFac,nPixels=nPixels,axes=axes,quantName=quantName,$
                           wholeBoxSlice=wholeBoxSlice,wholeBoxRemap=wholeBoxRemap
                           
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; locate subgroup index
  gcInd = getMatchedIDs(simParams=sP,haloID=haloID)
  
  if ~keyword_set(wholeBoxSlice) and ~keyword_set(wholeBoxRemap) then $
    fileName = 'sphMap.' + sP.savPrefix + str(sP.res) + '.' + $
      str(sP.snap) + '.h' + str(gcInd) + '.sf' + str(fix(sizeFac*10)) + '.sw' + $ 
      str(fix(sliceWidth*10)) + '.hf' + str(fix(hsmlFac*100)) + $
      '.px' + str(nPixels) + '.axes' + str(axes[0]) + str(axes[1]) + '_' + quantName
  if keyword_set(wholeBoxSlice) then $
    fileName = 'sliceMapOff1.' + sP.savPrefix + str(sP.res) + '.' + $
      str(sP.snap) + '.h' + str(gcInd) + '.sw' + str(fix(sliceWidth)) + $
      '.hf' + str(fix(hsmlFac*100)) + $
      '.px' + str(nPixels) + '.axes' + str(axes[0]) + str(axes[1]) + str(axes[2]) + '_' + quantName
  if keyword_set(wholeBoxRemap) then $
    fileName = 'boxRemap.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.h' + str(gcInd) + $
      '.rr' + str(fix(remapRatio[0]*100)) + '-' + str(fix(remapRatio[1]*100)) + '-' + str(fix(remapRatio[2]*100)) + $
      '.hf' + str(fix(hsmlFac*100)) + '.px' + str(nPixels[0]) + '-' + str(nPixels[1]) + $
      '.axes' + str(axes[0]) + str(axes[1]) + str(axes[2]) + '_' + quantName ; is actually partName

  saveFilename = sP.derivPath + 'binnedVals/' + fileName + '.sav'
      
  if file_test(saveFilename) then begin
    ; if map already exists, immediate return
    restore,saveFilename
    return, sphMap
  endif

  if ~keyword_set(wholeBoxSlice) and ~keyword_set(wholeBoxRemap) then begin
    ; load
    h = loadSnapshotHeader(sP=sP)
    
    partName  = strsplit(quantName,"_",/extract)
    partName  = partName[0]

    if partName eq 'dm' then $
      mass = float( h.massTable[ partTypeNum('dm') ] )
    if partName ne 'dm' then $
      mass = illustrisVisCutout(sP=sP, gcInd=gcInd, sizeFac=sizeFac, quantName='gas_mass')
      
    pos   = illustrisVisCutout(sP=sP, gcInd=gcInd, sizeFac=sizeFac, quantName=partName+'_pos')
    hsml  = illustrisVisCutout(sP=sP, gcInd=gcInd, sizeFac=sizeFac, quantName=partName+'_hsml')

    if quantName eq 'dm_dens' or quantName eq 'gas_dens' or quantName eq 'stars_dens' then $
      quant = fltarr(n_elements(hsml)) + 1.0 $
    else $
      quant = illustrisVisCutout(sP=sP, gcInd=gcInd, sizeFac=sizeFac, quantName=quantName)
  
    if quantName eq 'gas_xray' then begin
      w = where(quant eq 0.0,count)
      if count gt 0 then mass[w] = 0.0 ; do not consider zero L_X points
    endif
  
    ; load metadata
    metaFilePath = sP.derivPath + 'cutouts/cutout.' + sP.savPrefix + str(sP.res) + '.' + $
                   str(sP.snap) + '.h' + str(gcInd) + '.sf' + str(fix(sizeFac*10)) + '_meta.sav'
                 
    restore,metaFilePath
  
    ; apply hsml factor
    hsml = hsmlFac * temporary(hsml)

    ; image control
    nPixels2D = [nPixels,nPixels]
    boxSizeImg = fltarr(3)
    boxSizeImg[axes[0]] = 0.99 * 2 * boxSize ; avoid edge effects
    boxSizeImg[axes[1]] = 0.99 * 2 * boxSize
     
    projAxis = ( where( histogram(axes,min=0,max=2) eq 0 ) )[0]
    boxSizeImg[projAxis] = sliceWidth * 2 * haloRVir
  
    sphmap = calcSphMap(pos,hsml,mass,quant,$
                        boxSizeImg=boxSizeImg,boxSizeSim=0,boxCen=[0,0,0],$
                        nPixels=nPixels2D,axes=axes,ndims=3)
                        
    ; make save structure
    case quantName of
      'gas_temp'    : localMap  = alog10( sphMap.quant_out )
      'gas_dens'    : localMap  = alog10( sphMap.dens_out )
      'gas_entropy' : localMap  = alog10( sphMap.quant_out )
      'gas_metal'   : localMap = alog10( sphMap.quant_out )
      'gas_xray'    : localMap = alog10( sphMap.quant_out )
      'gas_sz_y'    : localMap = alog10( sphMap.quant_out )
      'gas_nelec'   : localMap = alog10( sphMap.quant_out )
      'gas_sfr'     : localMap = alog10( sphMap.quant_out )
      'dm_dens'     : localMap = alog10( sphMap.dens_out )
      'dm_vdisp'    : localMap = alog10( sphMap.quant_out )
      'dm_annih'    : localMap = alog10( sphMap.quant_out * sphMap.dens_out ) ; remove mass normalization
      else          : localMap = sphMap.quant_out
    endcase
    
    ; consistency with voronoi-projection: flip horizontal, then rotate 90 deg CW
    localMap = rotate(localMap,4) ; transpose
    
    sphMap = {haloRVir:haloRVir,boxSize:boxSize,haloMass:haloMass,fileName:fileName,quant:localMap}
  endif ;~wholeBoxSlice
    
  if keyword_set(wholeBoxRemap) then begin
    ; choose cuboid remapping matrix
    remapMatrix = findCuboidRemapInds(remapRatio=remapRatio,newBoxSize=newBoxSize,nPixels=nPixels)
    xyzCenterInv = 0.5 * newBoxSize ; box center in [0,L1/L2/L3] normalized coordinates
    newBoxSize *= sP.boxSize ; newBoxSize is relative
    
    ; where is the target halo (center in image)
    gc = loadGroupCat(sP=sP,/skipIDs,/skipOffsets)
    haloGrNr = gc.subgroupGrNr[ haloID ]
    xyzCenter = gc.groupPos[*, haloGrNr]
    
    print,'xyzCenter: ',xyzCenter
    print,'BoxSize: ',newBoxSize
    print,'nPixels: ',nPixels

    ; load positions (needed for every map)
    h = loadSnapshotHeader(sP=sP)
    partName = quantName
    print,'Loading pos...'
    
    reportMemory,msg="init"
    
    pos = loadSnapshotSubset(sP=sP,partType=partName,field='pos')
      
    reportMemory,msg="after pos"
      
    ; inverse map new box center into old [0,1]^3 box, and shift the particle positions along
    ; all three axes to place xyzCenter in the map center (only done within CalcSphMap when
    ; periodic BoxSize!=0 is set, and cannot do after transformation since new box is not periodic)
    print,'Shifting...'
    CalcBoxRemap, xyzCenterInv, sP.boxSize, remapMatrix, /inverse
    xyzCenterInv = float( xyzCenterInv * sP.boxSize ) ; [0,sP.boxSize]^3
    
    for j=0,2 do begin
      xyPos = reform( pos[j,*] - xyzCenter[j] + xyzCenterInv[j] )
      correctPeriodicPosVecs, xyPos, boxSize=sP.boxSize
      pos[j,*] = xyPos
    endfor
    xyPos = !NULL
    
    reportMemory,msg="after shift"
      
    ; apply cuboid remapping to positions
    CalcBoxRemap, pos, sP.boxSize, remapMatrix, /skipZ_stride2
    
    reportMemory,msg="after remap"
    
    ; load hsml and mass (needed for every map)
    print,'Loading hsml,mass...'
    
    if partName eq 'gas' then begin
      hsml = loadSnapshotSubset(sP=sP,partType=partName,field='hsml')
      hsml = (temporary(hsml) * 3.0 / (4*!pi))^(1.0/3.0) ;cellrad [ckpc]
      hsml = temporary(hsml) * hsmlFac
    endif
    if partName eq 'dm' or partName eq 'stars' then $
      hsml = loadSnapshotSubset(sP=sP,partType=partName,field='subfind_hsml') * hsmlFac ; ADDED HSMLFAC
      
    reportMemory,msg="after hsml"
      
    if partName eq 'dm' then begin
      mass = float( h.massTable[ partTypeNum('dm') ] )
    endif else begin
      mass = loadSnapshotSubset(sP=sP,partType=partName,field='mass')
    endelse
    
    reportMemory,msg="after mass"
   
    if partName eq 'gas'   then quantities = ['temp','dens','vmag','metal','ent','sz_y','xray']
    if partName eq 'dm'    then quantities = ['dens','vmag','vdisp']
    if partName eq 'stars' then quantities = ['dens']
    
    ; make save structure
    sphMap = {xyzCenter:xyzCenter,xyzCenterInv:xyzCenterInv,nPixels:nPixels,$
              boxSize:newBoxSize,fileName:fileName}
    
    foreach curQuant,quantities do begin
      ; load quantity
      print,'Loading '+curQuant+'...'
      reportMemory,msg="before quant"
      
      quant = !NULL
      
      if curQuant eq 'dens' then $
        quant = fltarr( n_elements(hsml) ) + 1.0
        
      if curQuant eq 'temp' then begin ; includes convertUtoTemp() for memory efficiency
        nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
        nelec = 4.0/(1.0 + 3.0 * units.hydrogen_massfrac + $
                     4.0 * units.hydrogen_massfrac * temporary(nelec)) * float(units.mass_proton) ; mu
        quant = loadSnapshotSubset(sP=sP,partType='gas',field='u') ; u
        quant = (5.0/3.0-1.0) * temporary(quant) * nelec * $
                float( units.UnitEnergy_in_cgs / (units.boltzmann * units.UnitMass_in_g) )
        nelec = !NULL
      endif
      
      if curQuant eq 'vmag' then begin ; quadrature addition in series for memory efficiency
        quant = loadSnapshotSubset(sP=sP, partType=partName,field='velx')^2.0
        quant = temporary(quant) + loadSnapshotSubset(sP=sP, partType=partName,field='vely')^2.0
        quant = temporary(quant) + loadSnapshotSubset(sP=sP, partType=partName,field='velz')^2.0
        quant = reform( sqrt( temporary(quant) ) )
      endif
      
      if curQuant eq 'ent' then begin ; includes CalcEntropyCGS() for memory efficiency
        u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
        quant = loadSnapshotSubset(sP=sP,partType='gas',field='dens')
      
        atime = snapNumToRedshift(sP=sP,/time)
        a3inv = 1.0 / (atime*atime*atime)
        u = (5.0/3.0-1.0) * temporary(u) * quant * a3inv * $
            float(units.UnitPressure_in_cgs / units.boltzmann) ; pressure [K/cm^3]
        quant = u / ( (temporary(quant) * $
                float(units.UnitDensity_in_cgs/units.mass_proton)*a3inv)^(5.0/3.0) ) ; ent [K cm^2]
        u = !NULL
      endif
        
      if curQuant eq 'xray' then begin
        quant = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
        
        quant = 4.0/(1.0 + 3.0 * units.hydrogen_massfrac + $
                     4.0 * units.hydrogen_massfrac * temporary(quant)) * float(units.mass_proton) ; mu
        quant = (5.0/3.0-1.0) * temporary(quant) * $
                float( units.UnitEnergy_in_cgs / (units.boltzmann * units.UnitMass_in_g) ) ; Temp/U
                
        quant = temporary(quant) * loadSnapshotSubset(sP=sP,partType='gas',field='u') ; Temp (K)
        quant = sqrt( temporary(quant) ) ; T^(1/2)
        
        ; zero contribution from gas with T<10^6 K
        w = where(quant le sqrt(1e6),count)
        if count gt 0 then begin
          quant[w] = 0.0
          mass[w] = 0.0 ; zero masses, do xray last so we don't have to restore
          print,'Zero ['+str(count)+'] of ['+str(n_elements(quant))+'] below 1e6 K.'
        endif
    
        ; load SFR, zero contribution from gas with SFR>0
        sfr = loadSnapshotSubset(sP=sP,partType='gas',field='sfr')
        w = where(sfr gt 0.0,count)
        if count gt 0 then begin
          quant[w] = 0.0
          mass[w] = 0.0
          print,'Zero ['+str(count)+'] of ['+str(n_elements(quant))+'] with nonzero SFR.'
        endif
        sfr = !NULL
        w = !NULL
        
        quant = temporary(quant) * loadSnapshotSubset(sP=sP,partType='gas',field='dens')  ; rho^1 * T^(1/2)
        u = loadSnapshotSubset(sP=sP,partType='gas',field='metallicity') ; Z
        
        Y = 0.25
        u = 4.0 / (8 - 5*Y - 6*temporary(u)) ; mu
        
        quant /= (u*u) ; rho/mu^2 * T^(1/2)
        u = !NULL
      endif
      
      if curQuant eq 'sz_y' then begin
        quant = loadSnapshotSubset(sP=sP,partType='gas',field='u') ; u
        quant = temporary(quant) * loadSnapshotSubset(sP=sP,partType='gas',field='nelec') ; u * ne
        metalmu = loadSnapshotSubset(sP=sP,partType='gas',field='metallicity')
        
        Y = 0.25
        metalmu = 4.0 / (8 - 5*Y - 6*temporary(metalmu)) ; mu
        quant = temporary(quant) * metalmu ; u * ne * mu
        metalmu = !NULL
      endif
      
      if n_elements(quant) eq 0 then $ ; normal load
        quant = loadSnapshotSubset(sP=sP,partType=partName,field=curQuant)
        
      reportMemory,msg="after quant"

      ; do sph projection and add contribution to map
      localMap = calcSphMap(pos,hsml,mass,quant,$
                            boxSizeImg=newBoxSize,boxSizeSim=0,boxCen=newBoxSize*0.5,$ 
                            nPixels=nPixels,axes=[0,1],ndims=3)
      
      reportMemory,msg="after map"

      ; add to save structure
      case curQuant of
        'temp'  : localMap = alog10( localMap.quant_out )
        'dens'  : localMap = alog10( localMap.dens_out )
        'ent'   : localMap = alog10( localMap.quant_out )
        'vdisp' : localMap = alog10( localMap.quant_out )
        else    : localMap = localMap.quant_out
      endcase
      
      sphMap = mod_struct( sphMap, partName + '_' + curQuant, localMap )
      
    endforeach
    
  endif
    
  if keyword_set(wholeBoxSlice) then begin
    ; load
    h = loadSnapshotHeader(sP=sP)
    partName = strsplit(quantName,"_",/extract)
    partName = partName[0]
   
    pos  = illustrisSliceCutout( sP=sP, gcInd=gcInd, sliceWidth=sliceWidth, axes=axes, $
                                 quantName=partName+'_pos')
    hsml = illustrisSliceCutout( sP=sP, gcInd=gcInd, sliceWidth=sliceWidth, axes=axes, $
                                 quantName=partName+'_hsml')
                                  
    if partName eq 'gas' or partName eq 'stars' then $
      mass = illustrisSliceCutout( sP=sP, gcInd=gcInd, sliceWidth=sliceWidth, axes=axes, $
                                   quantName=partName+'_mass')
    if partName eq 'dm' then $
      mass = fltarr(n_elements(hsml)) + h.massTable[ partTypeNum('dm') ]
    
    if quantName eq 'gas_dens' or quantName eq 'dm_dens' or quantName eq 'stars_dens' then begin
      quant = fltarr(n_elements(mass)) + 1.0
    endif else begin
      ; mass-weighted quantity
      quant = illustrisSliceCutout( sP=sP, gcInd=gcInd, sliceWidth=sliceWidth, axes=axes, $
                                  quantName=quantName)
    endelse
           
    if quantName eq 'gas_xray' then begin
      w = where(quant eq 0.0,count)
      if count gt 0 then mass[w] = 0.0 ; do not consider zero L_X points
    endif
           
    ; get box center
    metaFilePath = sP.derivPath + 'cutouts/sliceOff1.' + sP.savPrefix + str(sP.res) + '.' + $
      str(sP.snap) + '.h' + str(gcInd) + '.sw' + str(fix(sliceWidth)) + '.axes' + $
      str(axes[0]) + str(axes[1]) + str(axes[2]) + '_meta.sav'
                 
    restore,metaFilePath
           
    ; apply hsml factor (for hsml's derived from cell volumes only)
    if partName eq 'gas' then hsml = hsmlFac * temporary(hsml)

    ; image control
    nPixels2D  = [nPixels,nPixels]
    boxSizeImg = replicate(sP.boxSize,3)
    
    sphMap = calcSphMap(pos,hsml,mass,quant,$ 
                        boxSizeImg=boxSizeImg,boxSizeSim=sP.boxSize,boxCen=xyzCenter,$ 
                        nPixels=nPixels2D,axes=[0,1],ndims=3)
                        
    ; make save structure
    case quantName of
      'gas_dens'   : localMap = alog10( sphMap.dens_out )
      'gas_temp'   : localMap = alog10( sphMap.quant_out )
      'dm_dens'    : localMap = alog10( sphMap.dens_out )
      'stars_dens' : localMap = alog10( sphMap.dens_out )
      'dm_annih'   : localMap = alog10( sphMap.quant_out * sphMap.dens_out ) ; undo mass normalization
      else         : localMap = sphMap.quant_out
    endcase
      
    sphMap = {xyzCenter:xyzCenter,boxSize:sP.boxSize,fileName:fileName,quant:localMap}
      
  endif ; wholeBoxSlice
    
  ; save and return
  save,sphMap,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  return, sphMap
end

; catepillarMakeMap(): handle CalcSphMap-based map making for Grif's catepillar (Aquarius-like) runs
                           
function catepillarMakeMap, sP=sP,hsmlFac=hsmlFac,nPixels=nPixels,axes=axes,quantName=quantName
                           
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  boxSize = 2000.0 ; choose cut  
  
  ; load
  fileName = 'sphMap.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.hf' + str(fix(hsmlFac*100)) + '.bs' + str(fix(boxSize)) + $
    '.px' + str(nPixels) + '.axes' + str(axes[0]) + str(axes[1]) + '_' + quantName

  saveFilename = sP.derivPath + 'binnedVals/' + fileName + '.sav'
      
  if file_test(saveFilename) then begin
    ; if map already exists, immediate return
    restore,saveFilename
    return, sphMap
  endif

  ; load
  h = loadSnapshotHeader(sP=sP)
  
  partNameS  = strsplit(quantName,"_",/extract)
  partName   = partNameS[0]
  rQuantName = partNameS[1]

  mass = float( h.massTable[ partTypeNum('dm') ] )    
  pos  = loadSnapshotSubset(sP=sP,partType='dm',field='pos') ;/doublePrec
  
  pos -= 50.0 ; center (Mpc units)
  pos *= 1000.0 ; kpc

  w = where(abs(pos[0,*]) le boxSize*0.5 and $
            abs(pos[1,*]) le boxSize*0.5 and $
            abs(pos[2,*]) le boxSize*0.5,count)
            
  print,'cut',count,n_elements(pos[0,*])
  
  pos = pos[*,w]
  
  ; CalcHSML: offset into positive [0,boxSize]
  hsmlSaveName = sP.derivPath + 'binnedVals/hsml.sav'
  
  if ~file_test(hsmlSaveName) then begin
    pos += boxSize*0.5
    
    hsml = calcHSML(pos, ndims=3, nNGB=32, boxSize=boxSize)
    ;hsml = loadSnapshotSubset(sP=sP,partType='dm',field='subfind_hsml')

    pos -= boxSize*0.5
    
    save,hsml,filename=hsmlSaveName
    print,'saved hsml'
  endif else begin
    restore,hsmlSaveName,/verbose
  endelse
  
  hsml = hsml[w]
  
  ; load quantity
  if quantName eq 'dm_dens' then $
    quant = fltarr(n_elements(hsml)) + 1.0 $
  else $
    quant = loadSnapshotSubset(sP=sP,partType='dm',field=rQuantName)

  ; apply hsml factor
  hsml = hsmlFac * temporary(hsml)

  ; image control
  nPixels2D = [nPixels,nPixels]
  boxSizeImg = fltarr(3)
  boxSizeImg[axes[0]] = 0.99 * boxSize ; avoid edge effects
  boxSizeImg[axes[1]] = 0.99 * boxSize
   
  projAxis = ( where( histogram(axes,min=0,max=2) eq 0 ) )[0]
  boxSizeImg[projAxis] = 0.99 * boxSize ;sliceWidth * 2 * haloRVir

  sphmap = calcSphMap(pos,hsml,mass,quant,$
                      boxSizeImg=boxSizeImg,boxSizeSim=0,boxCen=[0,0,0],$
                      nPixels=nPixels2D,axes=axes,ndims=3)
                      
  ; make save structure
  case quantName of
    'dm_dens'     : localMap = alog10( sphMap.dens_out )
    'dm_vdisp'    : localMap = alog10( sphMap.quant_out )
    'dm_annih'    : localMap = alog10( sphMap.quant_out * sphMap.dens_out ) ; remove mass normalization
    else          : localMap = sphMap.quant_out
  endcase
  
  sphMap = {boxSize:boxSize,fileName:fileName,quant:localMap}    
    
  ; save and return
  save,sphMap,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  return, sphMap
end

; getColorBarLabel(): helper, shared below

function getColorBarLabel, quantName

  labelText = ''
  
  if quantName eq 'velmag'   or quantName eq 'gas_vmag' or quantName eq 'vmag' then labelText = "|v| [km/s]"
  if quantName eq 'temp'     or quantName eq 'gas_temp' then labelText = "log T_{gas} [K]"
  if quantName eq 'entropy'  or quantName eq 'gas_entropy' then labelText = "log ( S ) [_{ }K cm^{2 }]"
  if quantName eq 'density'  or quantName eq 'gas_dens' then labelText = "log ( \rho_{gas} )"
  if quantName eq 'metal'    or quantName eq 'gas_metal' then labelText = "log Z"
  if quantName eq 'xray'     or quantName eq 'gas_xray' then labelText = "log L_X"
  if quantName eq 'sz_y'     or quantName eq 'gas_sz_y' then labelText = "log y_{SZ}"
  if quantName eq 'vrad'     or quantName eq 'gas_vrad' then labelText = "v_{rad} [km/s]"
  if quantName eq 'dm_dens'  then labelText = "log ( \rho_{DM} )"
  if quantName eq 'dm_density'  then labelText = "log ( \rho_{DM} )"
  if quantName eq 'dm_vmag'  then labelText = "v_{DM} [km/s]"
  if quantName eq 'dm_vdisp' then labelText = "\sigma_{v} (DM) [km/s]"
  if quantName eq 'dm_annih' then labelText = "log L_A (DM)"
  
  if labelText eq '' then message,'Error: Unrecognized quantName for colorbar label.'
  return, labelText

end  

; illustrisProjSingleHalo(): use Arepo-based projection for halo scale images

pro illustrisProjSingleHalo, i_in=i_in
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  sP = simParams(res=1820,run='illustris',redshift=0.0)
  ;sP = simParams(res=1,run='scylla',redshift=2.0)

  haloID    = 0      ; fof number, 0, 1000
  sizeFac   = 1.0    ; boxlength in units of rvir
  axes      = [0,1]  ; 01 02 12 (xy xz yz)
  
  ; plot configuration
  nPixels    = 3000   ; px
  sliceWidth = 1.0    ; depth of box to project through, in units of rvir (maximum is sizeFac)
  scaleBar   = 0.0    ; ckpc, 0 to disable
  colorBar   = 0      ; 0 to disable
  virCircle  = 0      ; 0 to disable
  
  ; I1820 sf20 z=0
  ;pConfigs = { $ ; jjg/physics/visspec ; saga/saga-01
  ;  p0 : {quantName:'temp',      mapMM:[6.8,7.6],   ga:0.6, nB:0, ctName:'blue-red2'} ,$
  ;  p1 : {quantName:'density',   mapMM:[-5.0,-2.0], ga:1.0, nB:0, ctName:'nclR/WhiteBlueGreenYellowRed-dnA'} ,$
  ;  p2 : {quantName:'entropy',   mapMM:[10.5,12.2], ga:1.0, nB:0, ctName:'pmR/f-23-28-3'} ,$
  ;  p3 : {quantName:'metal',     mapMM:[-2.8,-1.7], ga:1.7, nB:1, ctName:'wkp/tubs/nrwc'} ,$ ;pm/f-30-31-32
  ;  p4 : {quantName:'velmag',    mapMM:[250,1100],  ga:1.0, nB:0, ctName:'pm/f-34-35-36'} ,$
  ;  p5 : {quantName:'xray',      mapMM:[15.2,19.8], ga:1.5, nB:10,ctName:'red-temp'} ,$
  ;  p6 : {quantName:'sz_y',      mapMM:[2.5,5.1],   ga:1.6, nB:0, ctName:'ocR/zeu'} ,$
  ;  p7 : {quantName:'dm_density',mapMM:[-5.0,-0.6], ga:1.0, nB:0, ctName:'helix'}  ,$ ; 455
  ;  p8 : {quantName:'dm_vmag',   mapMM:[100,1500],  ga:2.0, nB:0, ctName:'pm/f-34-35-36'} ,$
  ;  p9 : {quantName:'dm_vdisp',  mapMM:[1.6,3.0],   ga:2.0, nB:0, ctName:'esri/events/fire_active_2'} ,$
  ;  p10: {quantName:'dm_annih',  mapMM:[-6.0,1.0],  ga:1.0, nB:0, ctName:'td/DEM_poster'}  $
  ;}
  
  ; I1820 sf10 z=0 (1rvir, 3000px for composite)
  pConfigs = { $
    p0 : {quantName:'temp',      mapMM:[6.8,7.6],   ga:0.6, nB:0, ctName:'blue-red2'} ,$
    p1 : {quantName:'density',   mapMM:[-3.9,-2.0], ga:1.0, nB:0, ctName:'nclR/WhiteBlueGreenYellowRed-dnA'} ,$
    p2 : {quantName:'entropy',   mapMM:[10.0,11.9], ga:1.0, nB:0, ctName:'pmR/f-23-28-3'} ,$
    p3 : {quantName:'metal',     mapMM:[-2.9,-1.7], ga:1.0, nB:0, ctName:'wkp/tubs/nrwc'} ,$ ; wkp/tubs/nrwc
    p4 : {quantName:'velmag',    mapMM:[250,1100],  ga:1.0, nB:0, ctName:'pm/f-34-35-36'} ,$
    p5 : {quantName:'xray',      mapMM:[16.7,19.8], ga:1.0, nB:10,ctName:'red-temp'} ,$
    p6 : {quantName:'sz_y',      mapMM:[3.3,5.0],   ga:1.0, nB:0, ctName:'ocR/zeu'} ,$
    p7 : {quantName:'dm_density',mapMM:[-4.0,-0.6], ga:1.0, nB:0, ctName:'helix'}  ,$ ; 455
    p8 : {quantName:'dm_vmag',   mapMM:[200,1500],  ga:2.0, nB:0, ctName:'pm/f-34-35-36'} ,$
    p9 : {quantName:'dm_vdisp',  mapMM:[100,1060],  ga:1.0, nB:0, ctName:'esri/events/fire_active_2'} ,$ ; REQUIRES UNLOG
    p10: {quantName:'dm_annih',  mapMM:[-6.0,1.0],  ga:1.0, nB:0, ctName:'td/DEM_poster'}  $
  }
  
  ; scylla
  ;pConfigs = { $
  ;  p0 : {quantName:'temp',      mapMM:[6.8,7.6],   ga:0.6, nB:0, ctName:'blue-red2'} ,$
  ;  p1 : {quantName:'density',   mapMM:[-5.9,-3.0], ga:0.8, nB:0, ctName:'nclR/WhiteBlueGreenYellowRed'} ,$
  ;  p2 : {quantName:'entropy',   mapMM:[10.5,12.2], ga:1.0, nB:0, ctName:'pmR/f-23-28-3'} ,$
  ;  p3 : {quantName:'metal',     mapMM:[-2.8,-1.7], ga:1.7, nB:1, ctName:'wkp/tubs/nrwc'} ,$ ;pm/f-30-31-32
  ;  p4 : {quantName:'velmag',    mapMM:[250,1100],  ga:1.0, nB:0, ctName:'pm/f-34-35-36'}  $
  ;}
  
  ; plot which?
  foreach i,i_in do begin
  
  quantName = pConfigs.(i).quantName
  mapMinMax = pConfigs.(i).mapMM
  gamma     = pConfigs.(i).ga
  nBottom   = pConfigs.(i).nB
  mapCtName = pConfigs.(i).ctName
    
  ; calculate projection using either sph kernel or voronoi raytracing via Arepo 
  vorMaps = illustrisArepoProj(sP=sP,haloID=haloID,sizeFac=sizeFac,sliceWidth=sliceWidth,$
                               nPixels=nPixels,axes=axes)
                           
  ; scaling
  quant_ind = (where( tag_names(vorMaps) eq strupcase(quantName) )) [0]
  print,quantName + ' minMax: ',minmax(vorMaps.(quant_ind))
                 
  sphMap2D = (vorMaps.(quant_ind)-mapMinMax[0])*(255.0-nBottom) / (mapMinMax[1]-mapMinMax[0]) > 0
  sphMap2D = fix(sphMap2D + nBottom) < 255 ; nBottom-255    
  
  ; rotate: aligns scylla with maller fig, and aligns vormaps with whole box slices
  sphMap2D = rotate(sphMap2D,4)
  
  ; plot
  xySize = 8 ; keep constant, charsize/pthick relative to pagesize
  density = nPixels / xySize
  
  fileName = vorMaps.fileName + '_' + quantName + '_' + str_replace(mapCtName,"/","-",/global) + '_mm' + $
             str_replace(str(string(mapMinMax[0],format='(f7.1)')),"-","n") + '-' + $
             str_replace(str(string(mapMinMax[1],format='(f7.1)')),"-","n") + $
             '_ga' + str(fix(gamma*10)) + '_nB' + str(nBottom)
             
  start_PS, sP.plotPath + fileName + '.eps', xs=xySize, ys=xySize   
    ; establish axes
    xyMinMax = [-vorMaps.boxSize, vorMaps.boxSize]
    cgPlot, [0], [0], /nodata, xrange=xyMinMax, yrange=xyMinMax, pos=[0,0,1,1], /xs, /ys, /noerase
    
    ; load color table and output image
    loadColorTable, mapCtName, gamma=gamma
    tv, sphMap2D, 0.0, 0.0, /normal

    ; circle at virial radius
    if virCircle gt 0 then $
      tvcircle, vorMaps.haloRVir, 0, 0, cgColor('white'), thick=0.8, /data    
  
    ; scale bar
    if scaleBar gt 0 then begin
      xpos = [xyMinMax[1]*0.96,xyMinMax[1]*0.96-scaleBar]
      ypos = replicate(xyMinMax[1]*0.96,2)
    
      cgText,mean(xpos),ypos*0.94,string(scaleBar,format='(i3)')+' kpc/h',$
        alignment=0.5,color=cgColor('white')
      oplot,xpos,ypos,color=cgColor('white'),thick=!p.thick ;+ 0.5*xySize
    endif
    
    ; color bar
    if colorBar gt 0 then begin
      !x.thick = 1.0
      !y.thick = 1.0
    
      labelText = getColorBarLabel(quantName)
     
      ; calculate position and draw
      offset = 0.01 & height = 0.035 & width = 0.45
      pos = [offset,0.5-width*0.5,offset+height,0.5+width*0.5]
      loadColorTable, mapCtName, gamma=gamma
      cgColorbar,position=pos,divisions=0,charsize=0.000001,bottom=nBottom,ticklen=0.00001
      
      ; colorbar labels
      cbLabels = str( string(mapMinMax,format='(f7.1)') )
      if (mapMinMax[1]/10.0) eq fix(mapMinMax[1]/10.0) then $
        cbLabels = str( string(mapMinMax,format='(I5)') ) ; no decimals
      rAdjust = 0.015*(strlen(cbLabels[1])-1)
      lAdjust = 0
      if mapMinMax[0] lt 0.0 then lAdjust = 0.01
      
      cgText,0.500,offset+height*0.3,textoidl(labelText),alignment=0.5,color=cgColor('black'),/normal
      cgText,0.5-width*0.5+0.03+lAdjust,        offset+height*0.3,cbLabels[0],alignment=0.5,color=cgColor('black'),/normal
      cgText,0.5+width*0.5-0.01-rAdjust,offset+height*0.3,cbLabels[1],alignment=0.5,color=cgColor('black'),/normal
    endif
    
  end_PS, density=density, pngResize=100, /deletePS
  
  endforeach ;i
  
  stop
end

; illustrisVisSingleHalo(): use sphMap-based approach for halo scale images

pro illustrisVisSingleHalo
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  ;sP = simParams(res=1820,run='illustris',redshift=0.0)
  ;sP = simParams(res=11,run='zoom_20mpc',hind=0,redshift=3.0)
  ;sP = simParams(res=1,run='scylla',redshift=2.0)
  sP = simParams(res=1,run='catepillar',redshift=1.3)

  haloID    = 0      ; fof number, 0, 1000
  sizeFac   = 2.1    ; boxlength in units of rvir
  axes      = [0,1]  ; 01 02 12 (xy xz yz)

  ; plot configuration
  nPixels    = 1024   ; px
  sliceWidth = 2.1    ; depth of box to project through, in units of rvir (maximum is sizeFac)
  scaleBar   = 250.0  ; ckpc, 0 to disable
  colorBar   = 0      ; 0 to disable
  
  ; 1820
  ;pConfigs = { $
  ;  p0 : {quantName:'gas_temp',    hF:2.5, mapMM:[5.5,7.6],   ga:1.2, nB:0,  ctName:'blue-red2'} ,$
  ;  p1 : {quantName:'gas_dens',    hF:2.5, mapMM:[-3.9,-1.4], ga:0.8, nB:0,  ctName:'nclR/WhiteBlueGreenYellowRed'} ,$ 
  ;  p2 : {quantName:'gas_entropy', hF:2.5, mapMM:[8.5,10.2],  ga:1.0, nB:0,  ctName:'pmR/f-23-28-3'} ,$
  ;  p3 : {quantName:'gas_metal',   hF:2.5, mapMM:[-2.8,-2.0], ga:1.0, nB:0,  ctName:'wkp/tubs/nrwc'} ,$ ;pm/f-30-31-32
  ;  p4 : {quantName:'gas_vrad',    hF:2.5, mapMM:[-350,350],  ga:0.5, nB:0,  ctName:'pm/f-34-35-36'} ,$
  ;  p5 : {quantName:'gas_xray',    hF:2.5, mapMM:[-3.8,-0.8], ga:1.0, nB:20, ctName:'red-temp'} ,$
  ;  p6 : {quantName:'gas_sz_y',    hF:2.5, mapMM:[4.7,6.4],   ga:1.8, nB:40, ctName:'ocR/zeu'} ,$ ; 1820
  ;  p7 : {quantName:'gas_nelec',   hF:2.5, mapMM:[0.06,0.07], ga:1.0, nB:0,  ctName:'rainbow'} ,$
  ;  p8 : {quantName:'gas_sfr',     hF:2.5, mapMM:[-9.0,-1.2], ga:1.0, nB:10, ctName:'rainbow'} ,$
  ;  p9 : {quantName:'dm_dens',     hF:0.5, mapMM:[-5.0,0.2],  ga:1.0, nB:0,  ctName:'helix'}  ,$ ; 455
  ;  p10: {quantName:'dm_vmag',     hF:0.5, mapMM:[300,1500],  ga:2.0, nB:0,  ctName:'pm/f-34-35-36'} ,$
  ;  p11: {quantName:'dm_vdisp',    hF:0.5, mapMM:[1.9,3.0],   ga:2.0, nB:0,  ctName:'esri/events/fire_active_2'} ,$
  ;  p12: {quantName:'dm_annih',    hF:1.0, mapMM:[-9.0,-1.0], ga:1.0, nB:0,  ctName:'td/DEM_poster'}  $ ; BAD for 1820, ok for 455?
  ;}
  
  ; zoom_L11
  pConfigs = { $
    p0 : {quantName:'gas_temp',    hF:2.5, mapMM:[4.4,6.2],   ga:1.0, nB:0,  ctName:'blue-red2'} ,$ ;
    p1 : {quantName:'gas_dens',    hF:2.5, mapMM:[-5.3,-2.5], ga:0.8, nB:0,  ctName:'nclR/WhiteBlueGreenYellowRed'} ,$
    p2 : {quantName:'gas_entropy', hF:2.5, mapMM:[5.0,8.5],   ga:1.0, nB:0,  ctName:'pmR/f-23-28-3'} ,$
    p3 : {quantName:'gas_vrad',    hF:2.5, mapMM:[-400,400],  ga:1.0, nB:0,  ctName:'brewer-redpurple'} ,$
    p4 : {quantName:'gas_xray',    hF:2.5, mapMM:[-4.5,-1.2], ga:1.0, nB:20, ctName:'red-temp'} ,$ 
    p5 : {quantName:'gas_sz_y',    hF:2.5, mapMM:[2.6,4.6],   ga:1.8, nB:40, ctName:'ocR/zeu'} ,$
    p6 : {quantName:'dm_dens',     hF:0.5, mapMM:[-5.5,-1.0],  ga:1.0, nB:0,  ctName:'dnelson/dmdens'}  ,$ ; helix
    p7 : {quantName:'dm_vmag',     hF:0.5, mapMM:[300,1500],  ga:2.0, nB:0,  ctName:'pm/f-34-35-36'} ,$
    p8 : {quantName:'dm_vdisp',    hF:0.5, mapMM:[1.9,3.0],   ga:2.0, nB:0,  ctName:'esri/events/fire_active_2'} ,$
    p9 : {quantName:'dm_annih',    hF:1.0, mapMM:[-9.0,-1.0], ga:1.0, nB:0,  ctName:'td/DEM_poster'}  $ 
  }
  
  ; plot which?
  foreach i,[6] do begin
  quantName = pConfigs.(i).quantName
  mapMinMax = pConfigs.(i).mapMM
  gamma     = pConfigs.(i).ga
  nBottom   = pConfigs.(i).nB
  mapCtName = pConfigs.(i).ctName
  hsmlFac   = pConfigs.(i).hF
  
  ; calculate projection using mass-weighted sph kernel
  if sliceWidth gt sizeFac or sliceWidth le 0.0 then message,'Error'
  
  ; illustris:
  ;sphMap = illustrisMakeMap(sP=sP,haloID=haloID,sizeFac=sizeFac,sliceWidth=sliceWidth,$
  ;                          hsmlFac=hsmlFac,nPixels=nPixels,axes=axes,quantName=quantName)
                            
  ; catepillar:
  sphMap = catepillarMakeMap(sP=sP,hsmlFac=hsmlFac,nPixels=nPixels,axes=axes,quantName=quantName)

  ; scaling
  w = where(~finite(sphMap.quant),count,comp=wc)
  if count gt 0 then begin
    print,'non finite: ',count,' of ',n_elements(sphMap.quant)
   sphMap.quant[w] = 1.0 * min( sphMap.quant[wc] )
  endif
  
  print,quantName+' minMax: ',minmax(sphMap.quant)
                 
  sphMap2D = (sphMap.quant-mapMinMax[0])*(255.0-nBottom) / (mapMinMax[1]-mapMinMax[0]) > 0
  sphMap2D = fix(sphMap2D + nBottom) < 255 ; nBottom-255    
  
  ; rotate?
  ;print,'ROTATING (scylla)'
  ;sphMap2D = rotate(sphMap2D,4)
  
  ; plot
  xySize = 8 ; keep constant, charsize/pthick relative to pagesize
  density = nPixels / xySize
  
  fileName = sphMap.fileName + '_' + str_replace(mapCtName,"/","-",/global) + '_mm' + $
             str_replace(str(string(mapMinMax[0],format='(f7.1)')),"-","n") + '-' + $
             str_replace(str(string(mapMinMax[1],format='(f7.1)')),"-","n") + $
             '_ga' + str(fix(gamma*10)) + '_nB' + str(nBottom)
             
  start_PS, sP.plotPath + fileName + '.eps', xs=xySize, ys=xySize
  
    ; establish axes
    xyMinMax = [-sphMap.boxSize, sphMap.boxSize]
    cgPlot, [0], [0], /nodata, xrange=xyMinMax, yrange=xyMinMax, pos=[0,0,1,1], /xs, /ys, /noerase
    
    ; load color table
    loadColorTable, mapCtName, gamma=gamma
    
    ; output image
    tv, sphMap2D, 0.0, 0.0, /normal

    ; circle at virial radius
    ;tvcircle, sphMap.haloRVir, 0, 0, cgColor('white'), thick=0.8, /data    
  
    ; contours?
    ;cgContour, smooth(sphMap2D,ceil(nPixels/100)), levels=[140,180,220,240], label=0, /onImage
  
    ; scale bar
    if scaleBar gt 0 then begin
      xpos = [xyMinMax[1]*0.96,xyMinMax[1]*0.96-scaleBar/0.7]
      ypos = replicate(xyMinMax[1]*0.96,2)
    
      cgText,mean(xpos),ypos*0.94,string(scaleBar,format='(i3)')+' ckpc',$
        alignment=0.5,color=cgColor('white')
      oplot,xpos,ypos,color=cgColor('white'),thick=!p.thick+0.5
    endif
    
    ; color bar
    if colorBar gt 0 then begin
      !x.thick = 1.0
      !y.thick = 1.0
    
      labelText = getColorBarLabel(quantName)
       
      ; calculate position and draw
      offset = 0.01 & height = 0.035 & width = 0.45
      pos = [offset,0.5-width*0.5,offset+height,0.5+width*0.5]
      loadColorTable,mapCtName,gamma=gamma
      cgColorbar,position=pos,divisions=0,charsize=0.000001,bottom=nBottom,ticklen=0.00001
      
      ; colorbar labels
      cbLabels = str( string(mapMinMax,format='(f7.1)') )
      if (mapMinMax[1]/10.0) eq fix(mapMinMax[1]/10.0) then $
        cbLabels = str( string(mapMinMax,format='(I5)') ) ; no decimals
      rAdjust = 0.015*(strlen(cbLabels[1])-1)
      lAdjust = 0
      if mapMinMax[0] lt 0.0 then lAdjust = 0.01
      
      cgText,0.500,offset+height*0.3,textoidl(labelText),alignment=0.5,color=cgColor('black'),/normal
      cgText,0.5-width*0.5+0.03+lAdjust,        offset+height*0.3,cbLabels[0],alignment=0.5,color=cgColor('black'),/normal
      cgText,0.5+width*0.5-0.01-rAdjust,offset+height*0.3,cbLabels[1],alignment=0.5,color=cgColor('black'),/normal
    endif
    
  end_PS, density=density, pngResize=100, /deletePS

  endforeach ;i

end

; illustrisVisSlice():

pro illustrisVisSlice, i_in=i_in
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  sP = simParams(res=1820,run='illustris',redshift=0.0)
  
  haloID     = 0       ; fof number, 0, 1000
  sliceWidth = 7500    ; slice extent in
  axes       = [0,1,2] ; image in xy, slice in z
  
  ; plot configuration
  nPixels    = 8192  ; px
  reduce     = 1.0    ; reduce final output image size (1.0 to disable)
  scaleBar   = 0      ; cMpc, 0 to disable
  
  pConfigs = { $ ; tailored to 1820, 16k px
    p0 : {qName:'gas_temp',    hF:2.5, mapMM:[3.0,7.2],    ga:1.2, nB:0, ctName:'blue-red2'} ,$
    p1 : {qName:'gas_dens',    hF:1.0, mapMM:[-4.7,-1.6],  ga:1.3, nB:0, ctName:'nclR/WhiteBlueGreenYellowRed-dnA'} ,$ ;  ; helix
    p2 : {qName:'gas_entropy', hF:2.5, mapMM:[8.0,11.2],   ga:2.0, nB:0, ctName:'pm/f-23-28-3-dnB' } ,$ ;
    p3 : {qName:'gas_metal',   hF:2.5, mapMM:[-5.0,-2.0],  ga:1.0, nB:5, ctName:'wkp/tubs/nrwc'} ,$ ; pm/f-30-31-32
    p4 : {qName:'gas_vmag',    hF:2.5, mapMM:[50.0,960.0], ga:1.0, nB:0, ctName:'pm/f-34-35-36'} ,$
    p5 : {qName:'gas_xray',    hF:2.5, mapMM:[-7.6,-2.5],  ga:1.5, nB:0, ctName:'red-temp'} ,$
    p6 : {qName:'gas_sz_y',    hF:2.5, mapMM:[2.5,5.6],    ga:0.5, nB:0, ctName:'ocR/zeu'} ,$
    p7 : {qName:'dm_dens',     hF:1.0, mapMM:[-2.8,0.0],   ga:1.0, nB:0, ctName:'dnelson/dmdens'} ,$
    p8 : {qName:'dm_vmag',     hF:1.0, mapMM:[50.0,960.0], ga:1.0, nB:0, ctName:'pm/f-34-35-36'} ,$
    p9 : {qName:'dm_vdisp',    hF:1.0, mapMM:[0.6,3.0],    ga:1.0, nB:0, ctName:'esri/events/fire_active_2'} ,$
    p10: {qName:'dm_annih',    hF:1.0, mapMM:[-12.0,-1.2], ga:1.0, nB:0, ctName:'td/DEM_poster'} ,$
    p11: {qName:'stars_dens',  hF:2.5, mapMM:[-10.5,-2.0], ga:3.0, nB:0, ctName:'brewerC-cool'} $
  }
  
  ;pConfigs = { $ ; tailored to 1820, 2048 px
  ;  p1 : {qName:'gas_dens',    hF:2.5, mapMM:[-2.9,0.1],   ga:1.3, nB:0, ctName:'nclR/WhiteBlueGreenYellowRed-dnA'} ,$ ;  ; helix
  ;  p7 : {qName:'dm_dens',     hF:1.0, mapMM:[-2.2,1.6],   ga:0.7, nB:0, ctName:'helix'} $
  ;}
  
  ; plot which?
  foreach i,i_in do begin
  
    quantName = pConfigs.(i).qName
    mapMinMax = pConfigs.(i).mapMM
    gamma     = pConfigs.(i).ga
    nBottom   = pConfigs.(i).nB
    hsmlFac   = pConfigs.(i).hF
    ctName    = pConfigs.(i).ctName
  
  ; calculate projection using either sph kernel or voronoi raytracing via Arepo
  sphMap = illustrisMakeMap(sP=sP,haloID=haloID,sizeFac=sizeFac,sliceWidth=sliceWidth,$
                            hsmlFac=hsmlFac,nPixels=nPixels,axes=axes,quantName=quantName,$
                            wholeBoxSlice=1)
  
  ; scaling
  if quantName eq 'dm_vdisp' or quantName eq 'gas_entropy' or quantName eq 'gas_metal' or $
     quantName eq 'gas_xray' or quantName eq 'gas_sz_y' then sphMap.quant = alog10( sphMap.quant )
  
  w = where(~finite(sphMap.quant),count,comp=wc)
  if count gt 0 then begin
    print,'non finite: ',count,' of ',n_elements(sphMap.quant)
   sphMap.quant[w] = 1.0 * min( sphMap.quant[wc] )
  endif
  
  print,quantName + ' minMax: ',minmax(sphMap.quant)

  sphMap2D = (sphMap.quant-mapMinMax[0])*(255.0-nBottom) / (mapMinMax[1]-mapMinMax[0]) > 0
  sphMap2D = fix(sphMap2D + nBottom) < 255 ; nBottom-255    
  
  ; plot
  fileName = sphMap.fileName + '_' + str_replace(ctName,"/","-",/global) + '_mm' + $
             str_replace(str(string(mapMinMax[0],format='(f7.1)')),"-","n") + '-' + $
             str_replace(str(string(mapMinMax[1],format='(f7.1)')),"-","n") + $
             '_ga' + str(fix(gamma*10)) + '_nB' + str(nBottom)
             
  xySize = 16
  density = nPixels / xySize * reduce ; 512 for 8192px, 90 for 1440px
  
  start_PS, sP.plotPath + fileName + '.eps', xs=xySize, ys=xySize
  
    ; establish axes
    xyMinMax = [0.0, sphMap.boxSize]
    cgPlot, [0], [0], /nodata, xrange=xyMinMax, yrange=xyMinMax, pos=[0,0,1,1], /xs, /ys, /noerase
    
    ; load color table
    loadColorTable, ctName, gamma=gamma
    
    ; output image
    tv, sphMap2D, 0.0, 0.0, /normal

    ; scale bar
    if scaleBar gt 0 then begin
      xpos = [xyMinMax[1]*0.98,xyMinMax[1]*0.98-scaleBar*1000.0]
      ypos = replicate(xyMinMax[1]*0.98,2)
    
      cgText,mean(xpos),ypos*0.97,string(scaleBar,format='(i3)')+' Mpc/h',$
        alignment=0.5,color=cgColor('white')
      oplot,xpos,ypos,color=cgColor('white'),thick=!p.thick+0.5
    endif
        
  end_PS, density=density, pngResize=100, /deletePS
  
  endforeach ; pConfigs
end

; illustrisVisBox(): sphMap entire box with a CuboidRemap applied

pro illustrisVisBox, i_in=i_in, res=res, nPixels=nPixels
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  sP = simParams(res=res,run='illustris',snap=123)
  ;sP = simParams(res=256,run='feedback',redshift=0.0)
  
  haloID     = 0         ; fof number, defines center of image
  axes       = [0,1,2]   ; image in xy, slice in z
  hsmlFac    = 10.0      ; increase arepo 'hsml' to decrease visualization noise
  scaleBar   = 50.0      ; cMpc, 0 to disable
  nBottom    = 0
  
  ;nPixels     = [1024,1024]
  remapRatio  = [2.44,2.44,0.168]
  ;nPixels    = [1920,1080]*5       ; (width,height), should have the same aspect ratio as remapRatio[0,1]
  ;remapRatio = [5.0,2.8125,0.0711]  ; must satisfy L1*L2*L3=1 constraint, last entry gives fractional z-width

  if res eq 1820 then $
  pConfigs = { $
    p0 : {pName:'gas',   qName:'temp',  mapMM:[2.9,7.6],    ga:1.2, ctName:'blue-red2'} ,$ ; h5/dkbluered
    p1 : {pName:'gas',   qName:'dens',  mapMM:[-6.0,1.1],   ga:1.0, ctName:'nclR/WhiteBlueGreenYellowRed'} ,$
    p2 : {pName:'gas',   qName:'ent',   mapMM:[8.0,11.0],   ga:2.0, ctName:'pm/f-23-28-3-dnB'} ,$
    p3 : {pName:'gas',   qName:'metal', mapMM:[-5.0,-1.5],  ga:1.0, ctName:'wkp/tubs/nrwc'} ,$
    p4 : {pName:'gas',   qName:'vmag',  mapMM:[50.0,960.0], ga:1.0, ctName:'pm/f-34-35-36'} ,$
    p5 : {pName:'gas',   qName:'xray',  mapMM:[-6.0,-1.0],  ga:1.5, ctName:'red-temp'} ,$
    p6 : {pName:'gas',   qName:'sz_y',  mapMM:[2.0,5.8],    ga:1.2, ctName:'ocR/zeu'} ,$
    p7 : {pName:'dm',    qName:'dens',  mapMM:[-3.5,2.6],   ga:0.6, ctName:'helix'} ,$ ;1820
    p8 : {pName:'dm',    qName:'vmag',  mapMM:[50.0,960.0], ga:1.0, ctName:'pm/f-34-35-36'} ,$
    p9 : {pName:'dm',    qName:'vdisp', mapMM:[1.1,2.9],    ga:1.0, ctName:'esri/events/fire_active_2'} ,$
    p10: {pName:'stars', qName:'dens',  mapMM:[0.1,4.5],    ga:1.5, ctName:'brewerC-cool'} $ ;1820
  }
  
  if res eq 455 then $
  pConfigs = { $
    p0 : {pName:'gas',   qName:'temp',  mapMM:[2.9,7.6],    ga:1.2, ctName:'blue-red2'} ,$ ; h5/dkbluered
    p1 : {pName:'gas',   qName:'dens',  mapMM:[-5.0,1.1],   ga:1.0, ctName:'nclR/WhiteBlueGreenYellowRed'} ,$
    p2 : {pName:'gas',   qName:'ent',   mapMM:[7.5,11.0],   ga:2.0, ctName:'pm/f-23-28-3-dnB'} ,$
    p3 : {pName:'gas',   qName:'metal', mapMM:[-6.0,-2.0],  ga:1.0, ctName:'wkp/tubs/nrwc'} ,$
    p4 : {pName:'gas',   qName:'vmag',  mapMM:[50.0,960.0], ga:1.0, ctName:'pm/f-34-35-36'} ,$
    p5 : {pName:'gas',   qName:'xray',  mapMM:[-6.5,-1.5],  ga:1.5, ctName:'red-temp'} ,$ ;pm/pm3d11
    p6 : {pName:'gas',   qName:'sz_y',  mapMM:[1.5,5.5],    ga:1.2, ctName:'ocR/zeu'} ,$
    p7 : {pName:'dm',    qName:'dens',  mapMM:[-2.5,2.1],   ga:0.6, ctName:'helix'} ,$  ;455
    p8 : {pName:'dm',    qName:'vmag',  mapMM:[50.0,960.0], ga:1.0, ctName:'pm/f-34-35-36'} ,$
    p9 : {pName:'dm',    qName:'vdisp', mapMM:[1.1,2.9],    ga:1.0, ctName:'esri/events/fire_active_2'} ,$
    p10: {pName:'stars', qName:'dens',  mapMM:[-12.0,-2.0], ga:1.5, ctName:'brewerC-cool'} $
  }
  
  ; plot which?
  foreach i,i_in do begin
  partName  = pConfigs.(i).pName
  quantName = pConfigs.(i).qName
  mapMinMax = pConfigs.(i).mapMM
  gamma     = pConfigs.(i).ga
  mapCtName = pConfigs.(i).ctName
  
  ; calculate projection using either sph kernel or voronoi raytracing via Arepo
  boxMaps = illustrisMakeMap(sP=sP,haloID=haloID,remapRatio=remapRatio,$
                             hsmlFac=hsmlFac,nPixels=nPixels,axes=axes,quantName=partName,$
                             wholeBoxRemap=1)

  ; scaling
  quant_ind = (where( tag_names(boxMaps) eq strupcase(partName) + '_' + strupcase(quantName) )) [0]
  if quantName eq 'metal' or quantName eq 'xray' or quantName eq 'sz_y' then $
    boxMaps.(quant_ind) = alog10( boxMaps.(quant_ind) )
    
  w = where(~finite(boxMaps.(quant_ind)),count,comp=wc)
  if count gt 0 then begin
    print,'non finite: ',count,' of ',n_elements(boxMaps.(quant_ind))
    boxMaps.(quant_ind)[w] = 1.0 * min( boxMaps.(quant_ind)[wc] )
  endif
    
  print,quantName + ' minMax: ',minmax(boxMaps.(quant_ind))

  sphMap2D = (boxMaps.(quant_ind)-mapMinMax[0])*(255.0-nBottom) / (mapMinMax[1]-mapMinMax[0]) > 0
  sphMap2D = fix(sphMap2D + nBottom) < 255 ; nBottom-255    

  ; plot
  xySize = boxMaps.nPixels * 0.1 ; 0.1 = inverse of density
  fileName = boxMaps.fileName + '_' + quantName + '_' + str_replace(mapCtName,"/","-",/global) + '_mm' + $
             str_replace(str(string(mapMinMax[0],format='(f7.1)')),"-","n") + '-' + $
             str_replace(str(string(mapMinMax[1],format='(f7.1)')),"-","n") + $
             '_ga' + str(fix(gamma*10))
             
  start_PS, sP.plotPath + fileName + '.eps', xs=xySize[0], ys=xySize[1]
  
    ; establish axes
    xyMinMax = boxMaps.boxSize[0:1]
    cgPlot, [0], [0], /nodata, xrange=[0,xyMinMax[0]], yrange=[0,xyMinMax[1]], $
      pos=[0,0,1,1], /xs, /ys, /noerase
    
    ; load color table
    loadColorTable, mapCtName, gamma=gamma
    
    ; output image
    tv, sphMap2D, 0.0, 0.0, /normal

    ; scale bar
    if scaleBar gt 0 then begin
      loadColorTable,'bw linear'
      xpos = [xyMinMax[0]*0.98,xyMinMax[0]*0.98-scaleBar*1000.0]
      ypos = replicate(xyMinMax[1]*0.97,2)
    
      cgText,mean(xpos),mean(ypos)*0.97,string(scaleBar,format='(i3)')+' Mpc/h',$
        alignment=0.5,color=cgColor('white'),charsize=!p.charsize * 15 * (nPixels[0]/1920)
      oplot,xpos,ypos,color=cgColor('white'),thick=!p.thick * 20 * (nPixels[0]/1920)
    endif
        
  end_PS, density=10, pngResize=100, /deletePS
  
  endforeach ;i
  
end
