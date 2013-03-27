; filamentSearch.pro
; locating filamentary structures/streams in halos (inc. using healpix spheres)
; dnelson apr.2012

; angularFoF(): angular friends of friends algorithm on a healpix map mask

function angularFoF, healpix_mask, discrete=discrete, verbose=verbose

  if ~keyword_set(discrete) then message,'Error: Only discrete (touching pixels) implemented.'
  
  fragCutoff = n_elements(healpix_mask) / 1500 ; ~30, minimum number of pixels to save group
  
  ; group membership mask
  if max(healpix_mask) eq 0 then message,'Error: Empty mask.'
  fofMemMask = intarr(n_elements(healpix_mask))
  
  Nside = npix2nside(n_elements(healpix_mask))
  if Nside eq -1 then message,'Error: Bad pixel count in mask'
  
  nextGroupID = 1
  
  ; find starting orphans
  orphanInds = where(healpix_mask eq 1B,orphanCount)
    
  ; process each orphan pixel once
  for i=0L,orphanCount-1 do begin
    pxInd = orphanInds[i]
    
    ; locate immediately neighboring pixels in mask
    neighbours_nest,Nside,pxInd,neighborInds,nNeighbors
    neighborInds = neighborInds[where(healpix_mask[neighborInds] eq 1B,nnCount)]
    
    ; candidate groupID as minimum of all neighbor+self groupIDs, or next available if all zero
    candID = fofMemMask[pxInd]
    if nnCount gt 0 then candID = min([candID,fofMemMask[neighborInds]])
    
    if candID eq 0 then begin
      candID = nextGroupID
      nextGroupID += 1
    endif
    
    ; assign groupID
    if nnCount gt 0 then fofMemMask[neighborInds] = candID
    fofMemMask[pxInd] = candID
    
  endfor
  
  ; determine number of fof groups and merge
  nGroups = max(fofMemMask)
  collisionFlag = 1
  count = 0L
  
  if keyword_set(verbose) then $
    print,'Found ['+str(nGroups)+'] total groups, need to merge...'

  ; loop while there are potentially touching groups
  while collisionFlag do begin
    count += 1
    if keyword_set(verbose) then print,count
    collisionFlag = 0
    
    ; process each grouped (originally orphaned) pixel
    for i=0L,orphanCount-1 do begin
      pxInd = orphanInds[i]
      
      ; locate immediately neighboring pixels and their groupIDs
      neighbours_nest,Nside,pxInd,neighborInds,nNeighbors
      neighborInds = neighborInds[where(healpix_mask[neighborInds] eq 1B,nnCount)]
      
      colGroupIDs = fofMemMask[pxInd]
      if nnCount gt 0 then colGroupIDs = [colGroupIDs,fofMemMask[neighborInds]]
      
      ; if there is more than one unique value, flag collision and replace all with minimum
      if nuniq(colGroupIDs) gt 1 then begin
        collisionFlag = 1
        candID = min(colGroupIDs)
        fofMemMask[pxInd] = candID
        if nnCount gt 0 then fofMemMask[neighborInds] = candID
      endif
    endfor
  endwhile
  
  ; now we have the minimal set, compress group IDs so they are [1,...,nGroups]
  oldGroupIDs = uniqvals(fofMemMask[where(healpix_mask eq 1B)])
  nextGroup = 1
  
  foreach oldGID,oldGroupIDs do begin
    w = where(fofMemMask eq oldGID,count)
    if count eq 0 then message,'Error'

    ; discard small fragments
    if count lt fragCutoff then begin
      fofMemMask[w] = 0
    endif else begin
      ; group size large enough, keep
      fofMemMask[w] = nextGroup
      nextGroup += 1
    endelse
  endforeach
  
  if keyword_set(verbose) then $
    print,'Found ['+str(max(fofMemMask))+'] final groups.'
  
  return, fofMemMask  
end

; haloFilamentCrossSec():

function haloFilamentCrossSec, sP=sP, subgroupID=subgroupID, verbose=verbose, radIndOffset=radIndOffset

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config  
  threshOverDensAll  = 1.5 ; minimum rho/meanrho for candidate filament pixel
  threshOverDensRvir = 2.0

  widthShells = 5        ; odd, how many shells to enforce overdensity over, centered at rvir
  wrapTol     = !pi/20.0 ; radians, unwrapping lat,long at edges
  querySize   = !pi/4.0  ; radians, cutout size around each filament to calculate radii to
  circCutoff  = 0.1  ; exclude highly irregular shaped cross sections
  thetaTol    = 0.05 ; cannot be too close to N/S poles for phi(long) unwrap
  
  ; load the density shells
  hsv = haloShellValue(sP=sP,partType='gas',valName='density',subgroupID=subgroupID,/cutSubS)
  
  ; transform each shell to the ratio of the shell median
  for radInd=0,hsv.nRadFacs-1 do begin
    hsv.value[*,radInd] = reform(hsv.value[*,radInd] / mean(hsv.value[*,radInd]))
  endfor
  
  ; locate r=rvir radInd
  rvirInd = where(hsv.radFacs eq 1.0,count)
  if count ne 1 then message,'Error: Failed to locate rvir shell.'
  
  minInd = rvirInd[0] - floor(widthShells/2.0)
  maxInd = rvirInd[0] + floor(widthShells/2.0)
  
  ; mask filament pixel candidates based on overdensity requirement across radii range
  odMask = bytarr(hsv.nPx) + 1B

  for radInd=minInd,maxInd do begin
    w = where(hsv.value[*,radInd] lt threshOverDensAll,count)
    if count gt 0 then odMask[w] = 0B
  endfor
  
  w = where(odMask eq 1B,count)
  if keyword_set(verbose) then $
    print,'['+str(count)+'] of ['+str(hsv.nPx)+'] global filament pixel candidates.'

  ; find filaments on one shell
  radInd = rvirInd[0] + radIndOffset
 
  cenPos     = [] ; density weighted centroid (x,y,z)
  cenPosLL   = [] ; theta,phi = colat,long
  effCirSize = [] ; effective area of each filament (rad^2)
  cirMeasure = [] ; circularity measure [0,1]  
  
  filPxInds = []
  pxNums    = []
  ckpcDists = []
  
  ; further remove from mask more stringent overdensity requirement at this shell
  w = where(hsv.value[*,radInd] lt threshOverDensRvir,count)
  if count gt 0 then odMask[w] = 0B
  
  if keyword_set(verbose) then begin
    start_PS,'testmask2_r'+str(radInd)+'.eps'
      plotMollweideProj,float(odMask),title="",bartitle="",minmax=[0,1]
    end_PS,/deletePS,pngResize=40
  endif
  
  w = where(odMask eq 1B,count)
  if keyword_set(verbose) then print,'['+str(count)+'] of ['+str(hsv.nPx)+'] filament pixel candidates.'
  if count eq 0 then message,'Error: Too stringent local cut to mask.'

  ; run angular friends of friends
  angfof = angularFoF(odMask,/discrete)
  nGroupsCand = max(angfof)
    
  ; determine a overdensity weighted center position of each filament on the shell
  for i=0,nGroupsCand-1 do begin
    pxInds = where(angfof eq i+1)
    pix2ang_nest, hsv.nSide, pxInds, pxTheta, pxPhi
    pxWts = hsv.value[pxInds,radInd]
    
    if max(pxTheta) gt !pi-thetaTol or min(pxTheta) lt thetaTol then continue ; pole skip
    
    ; detect wrap (theta in [0,pi], phi in [0,2pi]) and correct
    if min(pxTheta) lt wrapTol and max(pxTheta) gt !pi-wrapTol then begin
      w = where(pxTheta gt !pi*0.4 and pxTheta lt !pi*0.6,count)
      if count ne 0 then continue ; something odd
      
      ; move lower points up by pi
      wTheta = where(pxTheta lt !pi*0.5,countMoveTheta)
      pxTheta[wTheta] += !pi
    endif
    
    if min(pxPhi) lt wrapTol and max(pxPhi) gt 2*!pi-wrapTol then begin
      w = where(pxPhi gt 2*!pi*0.4 and pxPhi lt 2*!pi*0.6,count)
      if count ne 0 then continue ; something odd
      
      ; move lower points up by 2pi
      wPhi = where(pxPhi lt 2*!pi*0.5,countMovePhi)
      pxPhi[wPhi] += 2*!pi
    endif
    
    ; weighted mean for centroid
    wtSum = total(pxWts)
    cenTheta = total(pxTheta * pxWts) / wtSum
    cenPhi = total(pxPhi * pxWts) / wtSum
    
    ; maximum extent in both axes for effective circular size and circularity measure      
    maxExtent = gcDist([max(pxTheta),max(pxPhi)],[min(pxTheta),min(pxPhi)]) ;[lat,long]
    curECS = 2*!pi*(1-cos(maxExtent/2.0))
    curCM  = (4*!pi/hsv.nPx*n_elements(pxInds)) / curECS
    
    if keyword_set(verbose) then print,radInd,i,maxExtent,curECS,curCM
    
    if curCM lt circCutoff then continue ; circ measure skip

    ; move center back into range if necessary
    if cenTheta gt !pi then cenTheta -= !pi
    if cenPhi gt 2*!pi then cenPhi -= 2*!pi
    
    ; convert centroid to pixel to vector
    ang2pix_nest, hsv.nSide, cenTheta, cenPhi, cenPxInd
    pix2vec_nest, hsv.nSide, cenPxInd, cenVec3
    
    ; query circular disc around centroid pixel and recalculate center position (without mask bias)
    query_disc, hsv.nSide, cenVec3, maxExtent/2.0, qdPxInds, /nested
    
    pix2ang_nest, hsv.nSide, qdPxInds, pxTheta, pxPhi
    pxWts = hsv.value[qdPxInds,radInd]
    
    ; detect wrap (theta in [0,pi], phi in [0,2pi]) and correct
    if min(pxTheta) lt wrapTol and max(pxTheta) gt !pi-wrapTol then begin
      w = where(pxTheta gt !pi*0.4 and pxTheta lt !pi*0.6,count)
      if count ne 0 then continue ; something odd
      
      ; move lower points up by pi
      wTheta = where(pxTheta lt !pi*0.5,countMoveTheta)
      pxTheta[wTheta] += !pi
    endif
    
    if min(pxPhi) lt wrapTol and max(pxPhi) gt 2*!pi-wrapTol then begin
      w = where(pxPhi gt 2*!pi*0.4 and pxPhi lt 2*!pi*0.6,count)
      if count ne 0 then continue ; something odd
      
      ; move lower points up by 2pi
      wPhi = where(pxPhi lt 2*!pi*0.5,countMovePhi)
      pxPhi[wPhi] += 2*!pi
    endif
    
    ; weighted mean for centroid
    wtSum = total(pxWts)
    cenTheta2 = total(pxTheta * pxWts) / wtSum
    cenPhi2 = total(pxPhi * pxWts) / wtSum
    
    ; store center
    ang2pix_nest, hsv.nSide, cenTheta2, cenPhi2, cenPxInd
    pix2vec_nest, hsv.nSIde, cenPxInd, vec3xyz

    ; query a new selection around this center
    query_disc, hsv.nSide, vec3xyz, querySize, qdPxInds, /nested

    ; calculate distance to each pixel on the sphere in ckpc
    pix2ang_nest, hsv.nSide, qdPxInds, pxTheta, pxPhi
    
    ; detect wrap (theta in [0,pi], phi in [0,2pi]) and correct
    if min(pxTheta) lt wrapTol and max(pxTheta) gt !pi-wrapTol then begin
      w = where(pxTheta gt !pi*0.4 and pxTheta lt !pi*0.6,count)
      if count ne 0 then continue ; something odd
      
      ; move lower points up by pi
      wTheta = where(pxTheta lt !pi*0.5,countMoveTheta)
      pxTheta[wTheta] += !pi
    endif
    
    if min(pxPhi) lt wrapTol and max(pxPhi) gt 2*!pi-wrapTol then begin
      w = where(pxPhi gt 2*!pi*0.4 and pxPhi lt 2*!pi*0.6,count)
      if count ne 0 then continue ; something odd
      
      ; move lower points up by 2pi
      wPhi = where(pxPhi lt 2*!pi*0.5,countMovePhi)
      pxPhi[wPhi] += 2*!pi
    endif
    
    angDists = fltarr(n_elements(qdPxInds))
    for j=0,n_elements(qdPxInds)-1 do $
      angDists[j] = gcDist([cenTheta2,cenPhi2],[pxTheta[j],pxPhi[j]])
    
    ; store values
    effCirSize = [effCirSize, curECS]
    cirMeasure = [cirMeasure, curCM]
    cenPos = [[cenPos], [reform(float(vec3xyz) * hsv.radFacs[radInd] * hsv.rVir)]]
    cenPosLL = [[cenPosLL], [cenTheta2,cenPhi2]]
    
    filPxInds = [filPxInds, qdPxInds]
    pxNums    = [pxNums, n_elements(qdPxInds)]
    ckpcDists = [ckpcDists, angDists * hsv.rVir]
  endfor ; nGroups
  
  if keyword_set(verbose) then begin
    start_PS,'testmask2_r'+str(radInd)+'_fof.eps'
      plotMollweideProj,float(angfof),title="",bartitle="",minmax=[0,max(angfof)]
    end_PS,/deletePS,pngResize=40
  endif

  r = {}
  if n_elements(pxNums) gt 0 then $
  r = { nFilaments:n_elements(pxNums), cenPos:cenPos, cenPosLL:cenPosLL, effCirSize:effCirSize, $
        cirMeasure:cirMeasure, filPxInds:filPxInds, pxNums:pxNums, ckpcDists:ckpcDists, $
        radInd:radInd}

  return,r
end

; haloFilamentSearch(): apply filament search algorithm

function haloFilamentSearch, sP=sP, subgroupID=subgroupID, verbose=verbose

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config  
  threshOverDensAll  = 1.0 ;1.5 ; minimum rho/meanrho for candidate filament pixel
  threshOverDensRvir = 1.5 ;2.0
  
  widthShells      = 5        ; odd, how many shells to enforce overdensity over, centered at rvir
  minNumDetections = 3       ; number of detection points per filament required to fit & keep
  wrapTol          = !pi/20.0 ; radians, unwrapping lat,long at edges
  maxNGroups       = 10       ; max filament detections per shell
  circMeasCutoff   = 0.2      ; circularity measure minimum for acceptance
  maxExtentCutoff  = !pi/4.0 ; radians, maximum size of filament cross section to accept
  filDistTol       = !pi/40.0 ; radians, matching filament centroids between shells
  
  ; load the density shells
  hsv = haloShellValue(sP=sP,partType='gas',valName='density',subgroupID=subgroupID,/cutSubS)
  
  ; transform each shell to the ratio of the shell median
  for radInd=0,hsv.nRadFacs-1 do begin
    hsv.value[*,radInd] = reform(hsv.value[*,radInd] / mean(hsv.value[*,radInd]))
  endfor
  
  ; locate r=rvir radInd
  rvirInd = where(hsv.radFacs eq 1.0,count)
  if count ne 1 then message,'Error: Failed to locate rvir shell.'
  
  minInd = rvirInd[0] - floor(widthShells/2.0)
  maxInd = rvirInd[0] + floor(widthShells/2.0)
  
  ; mask filament pixel candidates based on overdensity requirement across radii range
  odMask = bytarr(hsv.nPx) + 1B

  for radInd=minInd,maxInd do begin
    w = where(hsv.value[*,radInd] lt threshOverDensAll,count)
    if count gt 0 then odMask[w] = 0B
  endfor
  
  w = where(odMask eq 1B,count)
  if keyword_set(verbose) then $
    print,'['+str(count)+'] of ['+str(hsv.nPx)+'] global filament pixel candidates.'
  
  if keyword_set(verbose) then begin
    start_PS,'testmask.eps'
      plotMollweideProj,float(odMask),title="",bartitle="",minmax=[0,1]
    end_PS,/deletePS,pngResize=40
  endif
  
  ; arrays
  radMasks   = bytarr(widthShells,hsv.nPx)
  angfof     = intarr(widthShells,hsv.nPx)
  nGroups    = intarr(widthShells)
  cenPos     = fltarr(widthShells,maxNGroups,3) ; density weighted centroid (x,y,z)
  cenPosLL   = fltarr(widthShells,maxNGroups,2) ; theta,phi = colat,long
  effCirSize = fltarr(widthShells,maxNGroups) ; effective area of each filament (rad^2)
  cirMeasure = fltarr(widthShells,maxNGroups) ; circularity measure [0,1]
  
  ; 2d source location at each plane
  k = 0
  for radInd=minInd,maxInd do begin
    radMasks[k,*] = odMask  
  
    ; further remove from mask more stringent overdensity requirement at this shell
    w = where(hsv.value[*,radInd] lt threshOverDensRvir,count)
    if count gt 0 then radMasks[k,w] = 0B
    
    if keyword_set(verbose) then begin
      start_PS,'testmask_r'+str(radInd)+'.eps'
        plotMollweideProj,float(reform(radMasks[k,*])),title="",bartitle="",minmax=[0,1]
      end_PS,/deletePS,pngResize=40
    endif
    
    w = where(radMasks[k,*] eq 1B,count)
    if keyword_set(verbose) then print,'['+str(count)+'] of ['+str(hsv.nPx)+'] filament pixel candidates.'
    if count eq 0 then message,'Error: Too stringent local cut to mask.'

    ; run angular friends of friends
    angfof[k,*] = angularFoF(reform(radMasks[k,*]),/discrete)
    nGroupsCand = max(angfof[k,*])
      
    ; determine a overdensity weighted center position of each filament on the shell
    nextGroup = 0
    for i=0,nGroupsCand-1 do begin
      pxInds = where(reform(angfof[k,*]) eq i+1)
      pix2ang_nest, hsv.nSide, pxInds, pxTheta, pxPhi
      pxWts = hsv.value[pxInds,radInd]
      
      ; detect wrap (theta in [0,pi], phi in [0,2pi]) and correct
      if min(pxTheta) lt wrapTol and max(pxTheta) gt !pi-wrapTol then begin
        w = where(pxTheta gt !pi*0.4 and pxTheta lt !pi*0.6,count)
        if count ne 0 then message,'Error: Something odd in theta unwrap.'
        
        ; move lower points up by pi
        wTheta = where(pxTheta lt !pi*0.5,countMoveTheta)
        pxTheta[wTheta] += !pi
      endif
      
      if min(pxPhi) lt wrapTol and max(pxPhi) gt 2*!pi-wrapTol then begin
        w = where(pxPhi gt 2*!pi*0.4 and pxPhi lt 2*!pi*0.6,count)
        if count ne 0 then message,'Error: Something odd in phi unwrap.'
        
        ; move lower points up by 2pi
        wPhi = where(pxPhi lt 2*!pi*0.5,countMovePhi)
        pxPhi[wPhi] += 2*!pi
      endif
      
      ; weighted mean for centroid
      wtSum = total(pxWts)
      cenTheta = total(pxTheta * pxWts) / wtSum
      cenPhi = total(pxPhi * pxWts) / wtSum
      
      ; maximum extent in both axes for effective circular size and circularity measure      
      maxExtent = gcDist([max(pxTheta),max(pxPhi)],[min(pxTheta),min(pxPhi)]) ;[lat,long]

      effCirSize[k,nextGroup] = 2*!pi*(1-cos(maxExtent/2.0))
      cirMeasure[k,nextGroup] = (4*!pi/hsv.nPx*n_elements(pxInds)) / effCirSize[k,nextGroup]
      
      if keyword_set(verbose) then print,radInd,i,maxExtent,cirMeasure[k,nextGroup]
      
      ; move center back into range if necessary
      if cenTheta gt !pi then cenTheta -= !pi
      if cenPhi gt 2*!pi then cenPhi -= 2*!pi
      
      ; convert centroid to pixel to vector
      ang2pix_nest, hsv.nSide, cenTheta, cenPhi, cenPxInd
      pix2vec_nest, hsv.nSide, cenPxInd, cenVec3
      
      ; query circular disc around centroid pixel and recalculate center position (without mask bias)
      query_disc, hsv.nSide, cenVec3, maxExtent/2.0, qdPxInds, /nested
      
      pix2ang_nest, hsv.nSide, qdPxInds, pxTheta, pxPhi
      pxWts = hsv.value[qdPxInds,radInd]
      
      ; detect wrap (theta in [0,pi], phi in [0,2pi]) and correct
      if min(pxTheta) lt wrapTol and max(pxTheta) gt !pi-wrapTol then begin
        w = where(pxTheta gt !pi*0.4 and pxTheta lt !pi*0.6,count)
        if count ne 0 then message,'Error: Something odd in theta unwrap.'
        
        ; move lower points up by pi
        wTheta = where(pxTheta lt !pi*0.5,countMoveTheta)
        pxTheta[wTheta] += !pi
      endif
      
      if min(pxPhi) lt wrapTol and max(pxPhi) gt 2*!pi-wrapTol then begin
        w = where(pxPhi gt 2*!pi*0.4 and pxPhi lt 2*!pi*0.6,count)
        if count ne 0 then message,'Error: Something odd in phi unwrap.'
        
        ; move lower points up by 2pi
        wPhi = where(pxPhi lt 2*!pi*0.5,countMovePhi)
        pxPhi[wPhi] += 2*!pi
      endif
      
      ; weighted mean for centroid
      wtSum = total(pxWts)
      cenTheta2 = total(pxTheta * pxWts) / wtSum
      cenPhi2 = total(pxPhi * pxWts) / wtSum
      
      ; store x,y,z of center position
      ang2pix_nest, hsv.nSide, cenTheta2, cenPhi2, cenPxInd
      pix2vec_nest, hsv.nSIde, cenPxInd, vec3xyz
      
      cenPos[k,nextGroup,*] = reform(float(vec3xyz) * hsv.radFacs[radInd] * hsv.rVir) ; xyz (ckpc) halo centered
      cenPosLL[k,nextGroup,0] = cenTheta2
      cenPosLL[k,nextGroup,1] = cenPhi2

      ; check we satisfied the circularity minimum
      if cirMeasure[k,nextGroup] ge circMeasCutoff and maxExtent lt maxExtentCutoff then begin
        nextGroup += 1
        nGroups[k] += 1
        if nextGroup gt maxNGroups-1 then message,'Error: Exceeded maxNGroups.'
      endif
    endfor ; nGroups
    
    if keyword_set(verbose) then begin
      start_PS,'testmask_r'+str(radInd)+'_fof.eps'
        plotMollweideProj,float(reform(angfof[k,*])),title="",bartitle="",minmax=[0,max(angfof[k,*])]
      end_PS,/deletePS,pngResize=40
    endif
    
    k += 1
  endfor ; radInds
  
  ; if any nGroups are zero then we failed to locate any filaments
  if min(nGroups) eq 0 then begin
    print,'Warning: Failed to find any coherent filaments.'
    return,0
  endif
  
  ; matching arrays
  nFilaments   = nGroups[floor(widthShells/2.0)] ; at rvir
  filMatchInds = intarr(max(nGroups),widthShells) - 1
  filMinDists  = fltarr(max(nGroups),widthShells)
  
  for i=0,nFilaments-1 do filMatchInds[i,floor(widthShells/2.0)] = i ; at rvir filaments match themselves
  
  ; crossmatch sources between the planes: from rvir out
  for radInd=rvirInd[0],maxInd-1 do begin
    k = radInd - rVirInd[0] + floor(widthShells/2.0)

    ; loop over each filament
    for i=0,nFilaments-1 do begin
      curCenPos = cenPosLL[k,i,*]
      
      ; look to next shell outwards and find closest center match
      minDist = 999.9
      minFilInd  = -1
      for j=0,nGroups[k+1]-1 do begin
        otherCenPos = cenPosLL[k+1,j,*]
        otherDist = gcDist(curCenPos,otherCenPos)

        ; new closest?
        if otherDist lt minDist and otherDist lt filDistTol then begin
          if keyword_set(verbose) then print,k,i,j,otherDist
          
          ; has this k+1 candidate already been assigned to a previous filament at level k?
          if total(filMatchInds[*,k+1] eq j) gt 0 then begin
            ; if our current minDist is smaller than the previously used minDist, steal it
            if otherDist lt filMinDists[filMatchInds[i,k],k+1] then begin
              minDist = otherDist
              minFilInd = j
              
              ; set previous holder to unmatched
              w = where(filMatchInds[*,k+1] eq j,count)
              if count ne 1 then message,'Error2'
              filMatchInds[w,k+1] = -1
            endif
          endif else begin
            ; new unique assignment
            minDist = otherDist
            minFilInd = j
          endelse
        endif ; otherDist
      endfor ; nGroups
      
      ; save match
      if keyword_set(verbose) then if minFilInd eq -1 then print,'Warning: Match fail: ',k,k+1,i
      
      if filMatchInds[i,k] ne -1 then filMatchInds[filMatchInds[i,k],k+1] = minFilInd
      if filMatchInds[i,k] eq -1 then message,'Hmm error'
    endfor ; nFilaments
  endfor ; radInds
  
  ; crossmatch sources between the planes: from rvir in
  for radInd=rvirInd[0],minInd+1,-1 do begin
    k = radInd - rVirInd[0] + floor(widthShells/2.0)

    ; loop over each filament
    for i=0,nFilaments-1 do begin
      curCenPos = cenPosLL[k,i,*]
      
      if filMatchInds[i,k] eq -1 then begin
        print,'Warning: Skipping orphan filament at this level.',k,i
        continue
      endif
      
      ; look to next shell outwards and find closest center match
      minDist = 999.9
      minFilInd  = -1
      for j=0,nGroups[k-1]-1 do begin
        otherCenPos = cenPosLL[k-1,j,*]
        otherDist = gcDist(curCenPos,otherCenPos)

        ; new closest?
        if otherDist lt minDist and otherDist lt filDistTol then begin
          if keyword_set(verbose) then print,k,i,j,otherDist
          
          ; has this k-1 candidate already been assigned to a previous filament at level k?
          if total(filMatchInds[*,k-1] eq j) gt 0 then begin
            ; if our current minDist is smaller than the previously used minDist, steal it
            if otherDist lt filMinDists[filMatchInds[i,k],k-1] then begin
              minDist = otherDist
              minFilInd = j
              
              ; set previous holder to unmatched
              w = where(filMatchInds[*,k-1] eq j,count)
              if count ne 1 then message,'Error3'
              filMatchInds[w,k-1] = -1
            endif
          endif else begin
            ; new unique assignment
            minDist = otherDist
            minFilInd = j
          endelse
        endif ; otherDist
      endfor ; nGroups
      
      ; save match
      if keyword_set(verbose) then if minFilInd eq -1 then print,'Warning: Match fail: ',k,k-1,i
      
      if filMatchInds[i,k] ne -1 then begin
        filMatchInds[filMatchInds[i,k],k-1] = minFilInd
        filMinDists[filMatchInds[i,k],k-1] = minDist
      endif
      
    endfor ; nFilaments
  endfor ; radInds
  
  filIntPts    = [] ; best fit intersection point for each filament
  filUnitVecs  = [] ; best fit unit vector along the direction for each filament
  
  ; for each filament, collect the shell intersection points and calculate the ODR fit
  nFound = 0
  for i=0,nFilaments-1 do begin
    xpts = [] & ypts = [] & zpts = []
    
    for k=0,widthShells-1 do begin
      ; if this filament wasn't linked to this shell, don't contribute coordinates to fit
      if filMatchInds[i,k] ne -1 then begin
        xpts = [xpts,cenPos[k,reform(filMatchInds[i,k]),0]]
        ypts = [ypts,cenPos[k,reform(filMatchInds[i,k]),1]]
        zpts = [zpts,cenPos[k,reform(filMatchInds[i,k]),2]]
      endif
    endfor
    
    ; keep filament if we found N or more shell intersections
    if n_elements(xpts) ge minNumDetections then begin
      ; orthogonal distance regression -> best fit line
      odrFit = fitODRPts3D(xpts,ypts,zpts)
      
      filIntPts = [[filIntPts],[odrFit.centroid]]
      filUnitVecs = [[filUnitVecs],[odrFit.line_unitvec]]
      nFound += 1
    endif
  endfor
  
  r = { nFilaments:nFound, filIntPts:filIntPts, filUnitVecs:filUnitVecs, sP:sP, subgroupID:subgroupID, $
        threshOverDensAll:threshOverDensAll, threshOverDensRvir:threshOverDensRvir, $
        widthShells:widthShells, minNumDetections:minNumDetections, wrapTol:wrapTol, $
        maxNGroups:maxNGroups, circMeasCutoff:circMeasCutoff, maxExtentCutoff:maxExtentCutoff, $
        filDistTol:filDistTol, filMatchInds:filMatchInds, filMinDists:filMinDists, $
        cenPosLL:cenPosLL, cirMeasure:cirMeasure}

  return,r

end

; haloRefineFilaments(): refine angularfof based filament directions with a overdensity maximization search

function haloRefineFilaments, sP=sP, hfs=hfs

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  cylRadiusSQ = 40.0*40.0  ; ckpc, pre-square to avoid sqrt
  radRange    = [0.75,1.5] ; r/rvir of points to consider
  uvStepSize  = 0.01
  uvStepRange = 0.2

  ; load gas positions and densities
  pos_gas  = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  dens_gas = loadSnapshotSubset(sP=sP,partType='gas',field='dens')
  u_gas    = loadSnapshotSubset(sP=sP,partType='gas',field='u')
  
  ; spatial subset
  gc = loadGroupCat(sP=sP,/readIDs)
  sgpos  = subgroupPosByMostBoundID(sP=sP)
  cenPos = sgpos[*,hfs.subgroupID]
  rVir   = gc.group_r_crit200[gc.subgroupGrNr[hfs.subgroupID]]
  
  rad_gas = periodicDists(cenPos,pos_gas,sP=sP)
  
  wRadCut = where(rad_gas le radRange[1]*rVir and rad_gas gt radRange[0]*rVir,sCount)

  rad_gas  = rad_gas[wRadCut] / rVir
  pos_gas  = pos_gas[*,wRadCut]
  dens_gas = dens_gas[wRadCut]
  u_gas    = u_gas[wRadCut]

  ; remove substructures
;  ids_gas = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
;  ids_gas = ids_gas[wRadCut]
      
  ; make a list of satellites of this halo and their particle ids
;  nSubs    = gc.groupNSubs[gc.subgroupGrNr[hfs.subgroupID]]
;  firstSub = gc.groupFirstSub[gc.subgroupGrNr[hfs.subgroupID]]
;  satGCids = lindgen(nSubs-1) + firstSub + 1
;  satPIDs = gcPIDList(gc=gc,select='sec',valGCids=satGCids,partType='gas')
      
  ; remove the intersection of (satPIDs,loc_ids) from posval
;  calcMatch,satPIDs,ids_gas,sat_ind,ids_ind,count=count
;  sat_ind = !NULL
      
;  all = bytarr(n_elements(ids_gas))
;  if count gt 0 then all[ids_ind] = 1B
;  wSubSComp = where(all eq 0B, ncomp)
      
;  print,'Substructures cut ['+str(count)+'] of ['+str(n_elements(ids_gas))+'] have left: '+str(ncomp)
   
;  ids_ind = !NULL
;  ids_gas = !NULL
  wRadCut = !NULL
  sgpos = !NULL
  gc = !NULL

  ; take non-substructure cut and make positios relative to halo center
  ;if ncomp gt 0 then begin
  ;  rad_gas  = rad_gas[*,wSubSComp]
  ;  pos_gas  = pos_gas[*,wSubSComp]
  ;  dens_gas = dens_gas[wSubSComp]
  ;  u_gas    = u_gas[wSubSComp]
  ;endif
  
  pos_gas[0,*] -= cenPos[0]
  pos_gas[1,*] -= cenPos[1]
  pos_gas[2,*] -= cenPos[2]
    
  ; arrays
  nUvSteps = fix(uvStepRange / uvStepSize)
  nTrials = long(long(nUvSteps) * nUvSteps * nUvSteps)
  
  newUnitVecs = fltarr(3,hfs.nFilaments)
  
  ; loop over each filament
  for i=0,hfs.nFilaments-1 do begin
    print,i
    ; generate steps on filUnitVec
    trialUnitVecs0 = findgen(nUvSteps)/(nUvSteps-1) * uvStepRange - uvStepRange/2.0 + hfs.filUnitVecs[0,i]
    trialUnitVecs1 = findgen(nUvSteps)/(nUvSteps-1) * uvStepRange - uvStepRange/2.0 + hfs.filUnitVecs[1,i]
    trialUnitVecs2 = findgen(nUvSteps)/(nUvSteps-1) * uvStepRange - uvStepRange/2.0 + hfs.filUnitVecs[2,i]
    
    densFacs = fltarr(nTrials)
    maxInd = 0 & maxInd0 = 0 & maxInd1 = 0 & maxInd2 = 0
    
    ; loop over all trials
    for j=0UL,nTrials-1 do begin
      ind0 = fix(j mod nUvSteps)
      ind1 = fix(j/nUvSteps mod nUvSteps)
      ind2 = fix(j/(nUvSteps*nUvSteps))
      
      ;print,'['+string(j,format='(i3)')+' / '+str(nTrials)+'] '+str(ind0)+' '+str(ind1)+' '+str(ind2)
      
      ; calculate two line points x1,x2
      x1 = hfs.filIntPts[*,i] - 100.0*[trialUnitVecs0[ind0],trialUnitVecs1[ind1],trialUnitVecs2[ind2]]
      x2 = hfs.filIntPts[*,i] + 100.0*[trialUnitVecs0[ind0],trialUnitVecs1[ind1],trialUnitVecs2[ind2]]
      
      ; parametric solution and minimum distance to line (gas)  
      n21  = (x2[0]-x1[0])^2.0 + (x2[1]-x1[1])^2.0 + (x2[2]-x1[2])^2.0 ;n21=transpose(x1-x2)#(x2-x1)
      n10  = reform( (x1[0]-pos_gas[0,*])^2.0 + (x1[1]-pos_gas[1,*])^2.0 + (x1[2]-pos_gas[2,*])^2.0 )
      dotp = reform( (x1[0]-pos_gas[0,*])*(x2[0]-x1[0]) + (x1[1]-pos_gas[1,*])*(x2[1]-x1[1]) + $
                     (x1[2]-pos_gas[2,*])*(x2[2]-x1[2]) )
    
      t_gas = -1.0 * dotp / n21
      d_gas = ( n10 * n21 - dotp^2.0 ) / n21
      ;d_gas = sqrt(d_gas)
      
      ; select gas particles near the filament line
      w = where(d_gas le cylRadiusSQ,count)
      
      ; calculate "goodness of fit" function as mean(density/temperature) to bias towards
      ; low temperature mass overdensities with cylindrical geometry
      if count gt 0 then begin
        densFacs[j] = mean(dens_gas[w]/u_gas[w])
        
        ; if greater than the baseline mark as new best unitvecs
        if densFacs[j] gt densFacs[0] then begin
          maxInd = j & maxInd0 = ind0 & maxInd1 = ind1 & maxInd2 = ind2
        endif
      endif
    endfor
    
    ; if we found a maxInd other than the starting unit vector, add as best
    if maxInd eq 0 then print,'Warning: No better match?'
    
    newUnitVecs[0,i] = trialUnitVecs0[maxInd0]
    newUnitVecs[1,i] = trialUnitVecs1[maxInd1]
    newUnitVecs[2,i] = trialUnitVecs2[maxInd2]
    
  endfor

  hfs.filUnitVecs = newUnitVecs
  return,hfs
end

; makeFilamentProfile(): run filament search and create radial and cross sectional profiles (individual)

function makeFilamentProfile, sP=sP, subgroupID=subgroupID

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config  
  cylMaxRadius = 100.0
  cylRadRange  = [0.1,2.0]
  
  saveFilename = sP.derivPath+'hFil/hFil.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(subgroupID)+'.sav'
 
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif
  
  ; run filament search and refine procedure
  hfs = haloFilamentSearch(sP=sP,subgroupID=subgroupID)
  hfs = haloRefineFilaments(sP=sP,hfs=hfs)

  ; temp: make string
  uvStr = '['
  posStr = '['
  for i=0,hfs.nFilaments-1 do begin
    uvStr = uvStr + '['+string(hfs.filUnitVecs[0,i],format='(f6.3)')+','+$
                        string(hfs.filUnitVecs[1,i],format='(f6.3)')+','+$
                        string(hfs.filUnitVecs[2,i],format='(f6.3)')+']'
    posStr = posStr + '['+string(hfs.filIntPts[0,i],format='(f7.2)')+','+$
                          string(hfs.filIntPts[1,i],format='(f7.2)')+','+$
                          string(hfs.filIntPts[2,i],format='(f7.2)')+']'
                          
    if i ne hfs.nFilaments-1 then uvStr = uvStr + ','                   
    if i ne hfs.nFilaments-1 then posStr = posStr + ','
  endfor
  
  print,uvStr+']'
  print,posStr+']'
  
  ; load halo properties
  sgpos  = subgroupPosByMostBoundID(sP=sP)
  cenPos = sgpos[*,hfs.subgroupID]
  sgpos  = !NULL
  
  gc   = loadGroupCat(sP=sP,/readIDs)
  rVir = gc.group_r_crit200[gc.subgroupGrNr[hfs.subgroupID]]
  
  ; load gas positions and take radial subset
  pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  rad = periodicDists(cenPos,pos,sP=sP)
  wRadCut = where(rad le cylRadRange[1]*1.1*rVir and rad gt cylRadRange[0]*0.5*rVir,sCount)
  rad = rad[wRadCut] / rVir
  
  pos = pos[*,wRadCut]
  
  ; remove substructures
    ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    ids = ids[wRadCut]
    
    ; make a list of satellites of this halo and their particle ids
    satPIDs = gcPIDList(gc=gc,select='sec',partType='gas')
    
    ; remove the intersection of (satPIDs,ids)
    calcMatch,satPIDs,ids,sat_ind,ids_ind,count=count
    sat_ind = !NULL
    
    all = bytarr(n_elements(ids))
    if count gt 0 then all[ids_ind] = 1B
    wSubSComp = where(all eq 0B, ncomp)
    
    print,'Substructures cut ['+str(count)+'] of ['+str(n_elements(ids))+'] have left: '+str(ncomp)

    if ncomp gt 0 then begin
      wRadCut = wRadCut[wSubSComp]
      rad = rad[wSubSComp]
      pos = pos[*,wSubSComp]
    endif
    
    ids_ind = !NULL & ids = !NULL & satPIDs = !NULL & wSubSComp = !NULL
 
  ; load other gas properties and make same subset
  dens  = loadSnapshotSubset(sP=sP,partType='gas',field='dens')
  dens  = dens[wRadCut]
  u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
  u     = u[wRadCut]
  nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
  nelec = nelec[wRadCut]
  temp  = convertUtoTemp(u,nelec,/log)
  nelec = !NULL
  entr  = calcEntropyCGS(u,dens,sP=sP,/log)
  pres  = calcPressureCGS(u,dens,sP=sP,/log)
  u     = !NULL
  
  ; make velocity products
    vel = loadSnapshotSubset(sP=sP,partType='gas',field='vel')
    vel = vel[*,wRadCut]
  
    ; ----- vrad ------
    gVel = gc.subgroupVel[*,subgroupID]
    vel[0,*] = reform(vel[0,*] - gVel[0])
    vel[1,*] = reform(vel[1,*] - gVel[1])
    vel[2,*] = reform(vel[2,*] - gVel[2])
    
    ; make normalized position vector wrt halo center = vec(r) / ||r|| where r from particle to center
    ; means that radvel<0 is inflow and radvel>0 is outflow
    rnorm0 = reform(pos[0,*] - cenPos[0])
    rnorm1 = reform(pos[1,*] - cenPos[1])
    rnorm2 = reform(pos[2,*] - cenPos[2])
    
    correctPeriodicDistVecs, rnorm0, sP=sP
    correctPeriodicDistVecs, rnorm1, sP=sP
    correctPeriodicDistVecs, rnorm2, sP=sP
    
    ; denominator and do divide
    rnorm = sqrt(rnorm0*rnorm0 + rnorm1*rnorm1 + rnorm2*rnorm2)

    rnorm0 /= rnorm
    rnorm1 /= rnorm
    rnorm2 /= rnorm
    
    ; dot(vel,rnorm) gives the magnitude of the projection of vel onto vec(r)
    vrad = reform(vel[0,*]*rnorm0 + vel[1,*]*rnorm1 + vel[2,*]*rnorm2)
  
    ; ----- angm ------
    rnorm0 = reform(cenPos[0] - pos[0,*])
    rnorm1 = reform(cenPos[1] - pos[1,*])
    rnorm2 = reform(cenPos[2] - pos[2,*])
    
    correctPeriodicDistVecs, rnorm0, sP=sP
    correctPeriodicDistVecs, rnorm1, sP=sP
    correctPeriodicDistVecs, rnorm2, sP=sP
    
    ; angular momentum magnitude
    angm = fltarr(3,ncomp)
    angm[0,*] = rnorm1 * vel[2,*] - rnorm2 * vel[1,*]
    angm[1,*] = rnorm2 * vel[0,*] - rnorm0 * vel[2,*]
    angm[2,*] = rnorm0 * vel[1,*] - rnorm1 * vel[0,*]
    
    ; magnitude of specific angular momentum = rvec x vel
    angm = reform(sqrt(angm[0,*]*angm[0,*] + angm[1,*]*angm[1,*] + angm[2,*]*angm[2,*]))  
  
    rnorm0 = !NULL & rnorm1 = !NULL & rnorm2 = !NULL & rnorm = !NULL & vel = !NULL
  
  ; make positions relative to halo center
  pos[0,*] -= cenPos[0]
  pos[1,*] -= cenPos[1]
  pos[2,*] -= cenPos[2]  
  
  ; arrays
  fil_dist = []
  fil_rad  = []
  fil_dens = []
  fil_temp = []
  fil_entr = []
  fil_pres = []
  fil_vrad = []
  fil_angm = []
  fil_num  = []

  ; transform gas positions into a (distance from filament, radius from halo center) coord system
  for i=0,hfs.nFilaments-1 do begin
    ; calculate two line points x1,x2
    x1 = hfs.filIntPts[*,i] - 100.0*[hfs.filUnitVecs[0,i],hfs.filUnitVecs[1,i],hfs.filUnitVecs[2,i]]
    x2 = hfs.filIntPts[*,i] + 100.0*[hfs.filUnitVecs[0,i],hfs.filUnitVecs[1,i],hfs.filUnitVecs[2,i]]
    
    start_PS,'test_xy'+str(i)+'.eps',xs=6,ys=6
      cgplot,pos[0,*],pos[1,*],psym=3,xtitle="x pos",ytitle="y pos"
      cgplot,[x1[0],x2[0]],[x1[1],x2[1]],line=0,color=cgColor('red'),/overplot
    end_PS
    
    ; parametric solution and minimum distance to line (gas)  
    n21  = (x2[0]-x1[0])^2.0 + (x2[1]-x1[1])^2.0 + (x2[2]-x1[2])^2.0 ;n21=transpose(x1-x2)#(x2-x1)
    n10  = reform( (x1[0]-pos[0,*])^2.0 + (x1[1]-pos[1,*])^2.0 + (x1[2]-pos[2,*])^2.0 )
    dotp = reform( (x1[0]-pos[0,*])*(x2[0]-x1[0]) + (x1[1]-pos[1,*])*(x2[1]-x1[1]) + $
                   (x1[2]-pos[2,*])*(x2[2]-x1[2]) )
  
    t_gas = -1.0 * dotp / n21
    d_gas = ( n10 * n21 - dotp^2.0 ) / n21
    d_gas = sqrt(d_gas)

    ; take the subset near the filament
    w = where(d_gas le cylMaxRadius,count)
    print,i,count
    if count eq 0 then message,'Error: No gas found near filament.'
    
    fil_dist = [fil_dist,d_gas[w]]
    fil_rad  = [fil_rad,rad[w]]
    fil_dens = [fil_dens,dens[w]]
    fil_temp = [fil_temp,temp[w]]
    fil_entr = [fil_entr,entr[w]]
    fil_pres = [fil_pres,pres[w]]
    fil_vrad = [fil_vrad,vrad[w]]
    fil_angm = [fil_angm,angm[w]]
    
    fil_num = [fil_num,count]
  endfor
  
  r = { hfs:hfs, fil_dist:fil_dist, fil_rad:fil_rad, fil_dens:fil_dens, fil_temp:fil_temp, $
        fil_entr:fil_entr, fil_pres:fil_pres, fil_vrad:fil_vrad, fil_angm:fil_angm, fil_num:fil_num, $
        cylRadRange:cylRadRange, cylMaxRadius:cylMaxRadius }
  
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))

  return,r
end