; externalC.pro
; wrappers for external C-routines
; dnelson nov.2012

; calcHSML(): use CalcHSML external C-routine for the tree, neighbor searchings, and smoothing lengths

function calcHSML, Pos, ndims=ndims, nNGB=nNGB, boxSize=boxSize

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if (ndims ne 1 and ndims ne 2 and ndims ne 3) then stop

  ; prepare inputs
  npos = (size(pos))[2]

  NumPart = long(npos)
  Mass    = fltarr(npos)+1.0 ;dummy
  
  DesNumNgb    = long(nNGB) ; number of neighbors to use
  DesNumNgbDev = long(0)
  boxSize      = float(boxSize)
  HsmlGuess    = float(1.0)
  Softening    = float(1.0)
  
  hsml_out = fltarr(NumPart)
  
  ; call CalcHSML
  libName = '/n/home07/dnelson/idl/CalcHSML/CalcHSML_'+str(ndims)+'D.so'
  ret = Call_External(libName, 'CalcHSML', $
                      NumPart,Pos,Mass,DesNumNgb,DesNumNgbDev,boxSize,HsmlGuess,Softening,hsml_out, $
                      /CDECL)
   
  return, hsml_out
                     
end

; calcHSMLds(): use CalcHSMLds external C-routine to calculate the smoothing length needed to locate
;               a given number of neighbors with positions Pos around each position in SearchPos
;               where the two are generally different (e.g. the size of a tophat filter)

function calcHSMLds, Pos, SearchPos, ndims=ndims, nNGB=nNGB, boxSize=boxSize

  compile_opt idl2, hidden, strictarr, strictarrsubs

  npos = size(Pos)
  nsrc = size(SearchPos)

  if ndims ne 1 and ndims ne 2 and ndims ne 3 then message,'Error: Need ndims=1,2,3.'
  if npos[0] ne 2 or npos[1] ne 3 then message,'Error: Point position array shape.'
  if nsrc[0] ne 2 or nsrc[1] ne 3 then message,'Error: Search position array shape.'
  if npos[2] lt nNGB then message,'Error: Point count too low for nNGB.'

  ; prepare inputs
  NumPart   = long(npos[2])
  NumSearch = long(nsrc[2])
  
  DesNumNgb    = long(nNGB)     ; number of neighbors to use
  DesNumNgbDev = long(0)        ; deviation allowed
  boxSize      = float(boxSize) ; use zero for non-periodic search
  
  hsml_out = fltarr(NumSearch)
  
  ; make sure point arrays are float triples since we direct cast now
  Pos = float(Pos)
  SearchPos = float(SearchPos)
  
  ; call CalcHSMLds
  libName = '/n/home07/dnelson/idl/CalcHSMLds/CalcHSMLds_'+str(ndims)+'D.so'
  ret = Call_External(libName, 'CalcHSMLds', $
                      NumPart,Pos,NumSearch,SearchPos,DesNumNgb,DesNumNgbDev,boxSize,hsml_out, $
                      /CDECL)
                      
  return, hsml_out              
end

; calcTHVal(): use CalcTHVal external C-routine to calculate the tophat kernel estimated value of a 
;              function specified for each particle/cell over a given number of neighbors with positions 
;              Pos around each position in SearchPos, where the two are generally different (as in CalcHSMLds)

function calcTHVal, PosVal, SearchPos, thMode=thMode, ndims=ndims, nNGB=nNGB, boxSize=boxSize, $
                    weighted=weighted

  compile_opt idl2, hidden, strictarr, strictarrsubs

  npos = size(PosVal)
  nsrc = size(SearchPos)

  if ndims ne 1 and ndims ne 2 and ndims ne 3 then message,'Error: Need ndims=1,2,3.'
  if thMode ne 1 and thMode ne 2 and thMode ne 3 then message,'Error: Need thMode=1,2,3.'
  if npos[0] ne 2 then message,'Error: Point position array shape.'
  if nsrc[0] ne 2 or nsrc[1] ne 3 then message,'Error: Search position array shape.'
  if npos[2] lt nNGB then message,'Error: Point count too low for nNGB.'
  
  if ~keyword_set(weighted) and npos[1] ne 4 then message,'Error: Posval should be 4xN.'
  if keyword_set(weighted) and npos[1] ne 5 then message,'Error: Posval(wt) should be 5xN.'

  ; prepare inputs
  NumPart   = long(npos[2])
  NumSearch = long(nsrc[2])
  
  DesNumNgb    = long(nNGB)     ; number of neighbors to use
  boxSize      = float(boxSize) ; use zero for non-periodic search
  thMode       = fix(thMode)    ; 1=mean, 2=total, 3=total/volume (density)
  
  val_out = fltarr(NumSearch)
  
  ; make sure point arrays are 32bit float since we direct cast now
  PosVal    = float(PosVal)
  SearchPos = float(SearchPos)
  
  if ~keyword_set(weighted) then begin
    ; call CalcTHVal
    libName = '/n/home07/dnelson/idl/CalcTHVal/CalcTHVal_'+str(ndims)+'D.so'
    ret = Call_External(libName, 'CalcTHVal', $
                        NumPart,PosVal,NumSearch,SearchPos,DesNumNgb,thMode,boxSize,val_out, $
                        /CDECL)
  endif else begin
    ; call CalcTHValWt where posval=posvalwt
    libName = '/n/home07/dnelson/idl/CalcTHValWt/CalcTHValWt_'+str(ndims)+'D.so'
    ret = Call_External(libName, 'CalcTHValWt', $
                        NumPart,PosVal,NumSearch,SearchPos,DesNumNgb,thMode,boxSize,val_out, $
                        /CDECL)
  endelse
                      
  return, val_out              
end

; calcNN(): use CalcNN external C-routine for the tree and neighbor search
;           return the index of Pos_SrcTargs closest to each of Pos_SrcOrigs

function calcNN, Pos_SrcTargs, Pos_SrcOrigs, boxSize=boxSize, ndims=ndims

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; prepare inputs
  n_srcTargs = long( n_elements(Pos_SrcTargs[0,*]) )
  n_srcOrigs = long( n_elements(Pos_SrcOrigs[0,*]) )
  
  boxSize = float(boxSize)
  
  ; make sure floats for direct cast
  Pos_SrcTargs = float(Pos_SrcTargs)
  Pos_SrcOrigs = float(Pos_SrcOrigs)
  
  ; prepare return
  ind_out = lonarr(n_srcOrigs)
  
  ; call CalcNN
  libName = '/n/home07/dnelson/idl/CalcNN/CalcNN_'+str(ndims)+'D.so'
  ret = Call_External(libName, 'CalcNN', $
                      n_srcTargs,n_srcOrigs,Pos_SrcTargs,Pos_SrcOrigs,boxSize,ind_out, $
                      /CDECL)
    
  return, ind_out
end

; calcSphMap(): use CalcSphMap external C-routine to simultaneously calculate a map of projected 
;               density and some other mass weighted quantity (e.g. temperature) with the sph 
;               spline kernel (for non-periodic set boxsize=0)

function calcSphMap, pos, hsml, mass, quant, axes=axes, boxSizeImg=boxSizeImg, boxSizeSim=boxSizeSim, $
                     boxCen=boxCen, nPixels=nPixels, ndims=ndims

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; prepare inputs
  NumPart = (size(pos))[2]
  
  boxSizeImg = float(boxSizeImg)
  boxSizeSim = float(boxSizeSim)
  boxCen     = float(boxCen)
  nPixels    = fix(nPixels)
  axes       = fix(axes)
  
  ; check inputs
  if (n_elements(boxSizeImg) ne 3) then stop
  if (n_elements(boxSizeSim) ne 1) then stop
  if (n_elements(boxCen)  ne 3) then stop
  if (n_elements(nPixels) ne 2) then stop
  if (n_elements(axes)    ne 2) then stop
  
  if (size(pos))[0] ne 2 or (size(pos))[1] ne 3 then stop
  if (size(hsml))[0] ne 1 then stop
  if (size(mass))[0] ne 1 then stop
  if (size(quant))[0] ne 1 then stop
  
  if (size(pos))[2] ne (size(hsml))[1] or $
     (size(pos))[2] ne (size(mass))[1] or $
     (size(pos))[2] ne (size(quant))[1] then stop
  
  if (axes[0] ne 0 and axes[0] ne 1 and axes[0] ne 2) then stop
  if (axes[1] ne 0 and axes[1] ne 1 and axes[1] ne 2) then stop
  
  ; we direct cast so ensure everything is the right size
  pos   = float(pos)
  hsml  = float(hsml)
  mass  = float(hsml)
  quant = float(quant)
  
  ; make return
  dens_out  = fltarr(nPixels[0],nPixels[1])
  quant_out = fltarr(nPixels[0],nPixels[1])

  ; call CalcSphMap
  libName = '/n/home07/dnelson/idl/CalcSphMap/CalcSphMap_'+str(ndims)+'D.so'
  ret = Call_External(libName, 'CalcSphMap', $
                      NumPart,pos,hsml,mass,quant,dens_out,quant_out,$
                      boxSizeImg[0],boxSizeImg[1],boxSizeImg[2],boxSizeSim,$
                      boxCen[0],boxCen[1],boxCen[2],axes[0],axes[1],nPixels[0],nPixels[1],$
                      /CDECL)

  return, { dens_out:dens_out, quant_out:quant_out }
end

; calcGridData(): given points on a plane with coordinates (xx,yy) interpolate to a regular grid
; using the sph kernel (optional weights mm and quantities to be weighted qq)

function calcGridData, xx=xx, yy=yy, mm=mm, qq=qq, nPixels=nPixels, $
                       xMinMax=xmm, yMinMax=ymm, logY=logY, logX=logX, colNorm=colNorm

  compile_opt idl2, hidden, strictarr, strictarrsubs

  nNGB = 64
  
  ; copy to not override inputs
  xxx = xx
  yyy = yy
  xMinMax = xmm
  yMinMax = ymm
  
  ; take log(x)?
  if keyword_set(logX) then begin
    w = where(finite(alog10(xx)))
    yy = yy[w]
    xxx = alog10(xx[w])
    xMinMax = alog10(xmm)
  endif
  
  ; take log(y)?
  if keyword_set(logY) then begin
    w = where(finite(alog10(yy)))
    xxx = xxx[w]
    yyy = alog10(yy[w])
    yMinMax = alog10(ymm)
  endif
  
  ; prepare inputs
  NumPart = n_elements(xxx)
  
  boxSizeImg = [xMinMax[1]-xMinMax[0], yMinMax[1]-yMinMax[0], 1.0]
  boxSizeSim = 0.0
  boxCen     = [mean(xMinMax), mean(yMinMax), 0.0]
  nPixels    = fix(nPixels)
  axes       = fix([0,1])
  
  ; get hsml
  zPts = replicate(0.0, NumPart)
  pos = float(transpose([[xxx],[yyy],[zPts]]))
  
  hsml = fltarr(NumPart)
  hsml = calcHSML(pos,ndims=2,nNGB=nNGB,boxSize=0.0)
  
  ; optional weights and quantity to take weighted average of (otherwise just density)  
  if keyword_set(mm) then begin
    mass = float(mm)
  endif else begin
    mass  = replicate(1.0, NumPart)
  endelse
  
  if keyword_set(qq) then begin
    quant = float(qq)
  endif else begin
    quant = replicate(1.0, NumPart)
  endelse
  
  ; make return
  dens_out  = fltarr(nPixels[0],nPixels[1])
  quant_out = fltarr(nPixels[0],nPixels[1])

  ; call CalcSphMap
  libName = '/n/home07/dnelson/idl/CalcSphMap/CalcSphMap_2D.so'
  ret = Call_External(libName, 'CalcSphMap', $
                      NumPart,pos,hsml,mass,quant,dens_out,quant_out,$
                      boxSizeImg[0],boxSizeImg[1],boxSizeImg[2],boxSizeSim,$
                      boxCen[0],boxCen[1],boxCen[2],axes[0],axes[1],nPixels[0],nPixels[1],$
                      /CDECL)

  ; normalize column by column?
  if keyword_set(colNorm) then begin
    for i=0,nPixels[0]-1 do begin
      dens_out[i,*] /= max(dens_out[i,*],/nan)
      quant_out[i,*] /= max(quant_out[i,*],/nan)
    endfor
  endif
                      
  ; return
  xPts = linspace(xMinMax[0],xMinMax[1],nPixels[0])
  yPts = linspace(yMinMax[0],yMinMax[1],nPixels[1])
  
  xxPts = cmreplicate(xPts,nPixels[1])
  yyPts = transpose(cmreplicate(yPts,nPixels[0]))
  
  return, { dens_out:dens_out, quant_out:quant_out, xPts:xPts, yPts:yPts, xxPts:xxPts, yyPts:yyPts }

end

; calcCoolTime(): use the primordial cooling network (KWH) to calculate the cooling times of gas

function calcCoolTime, u, rho, nelec, flag=flag, scalefac=scalefac

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if keyword_set(scalefac) then begin
    if scalefac lt 0.0 or scalefac gt 1.0 then message,'error scalefac'
    scalefac = float(scalefac)
  endif else begin
    print,'Warning: Using no scalefactor (comoving integration effectively disabled).'
    scalefac = float(0.0)
  endelse
  
  if n_elements(flag) eq 0 then flag = 0
  
  if n_elements(u) ne n_elements(rho) or n_elements(u) ne n_elements(nelec) or n_elements(u) eq 0 then $
    message,'Error: Bad input array sizes.'

  ; prepare inputs
  npts  = long(n_elements(u))
  u     = float(u)
  rho   = float(rho)
  nelec = float(nelec)
  flag  = long(flag)

  cooltime_out = fltarr(npts)
  
  ; call CalcCoolTime
  libName = '/n/home07/dnelson/idl/CalcCoolTime/CalcCoolTime.so'
  ret = Call_External(libName, 'CalcCoolTime',npts,u,rho,nelec,cooltime_out,scalefac,flag,/CDECL)
   
  return, cooltime_out
                     
end

; estimateDensityTophat(): spatial density estimator for an input position tuple array and
;                          CONSTANT mass per particle, by using HSMLs as tophat filter sizes
; pos_search : if specified, estimate density at this point set not at the positions of the particles

function estimateDensityTophat, pos, pos_search=pos_search, mass=mass, $
                                ndims=ndims, nNGB=nNGB, boxSize=boxSize

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if ~keyword_set(nNGB) or ~keyword_set(ndims) then message,'Error: Must specify nNGB and ndims.'
  if ~keyword_set(mass) or n_elements(mass) ne 1 then message,'Error: Expected one constant mass.'
  
  ; calculate smoothing lengths
  if ~keyword_set(pos_search) then $
    hsml_out = calcHSML(pos,ndims=ndims,nNGB=nNGB,boxSize=boxSize)
  if keyword_set(pos_search) then $
    hsml_out = calcHSMLds(pos,pos_search,ndims=ndims,nNGB=nNGB,boxSize=boxSize)
  
  ; convert smoothing lengths to densities
  if (ndims eq 1) then hsml_out = mass * nNGB / (2.0 * temporary(hsml_out))
  if (ndims eq 2) then hsml_out = mass * nNGB / (!pi * temporary(hsml_out)^2.0)
  if (ndims eq 3) then hsml_out = mass * nNGB / (4.0*!pi/3.0 * temporary(hsml_out)^3.0)
  
  return,hsml_out
end

; calcSort(): identical behavior as idl sort(), i.e. return permutation indices to sort array 
;   in ascending order, unless inPlace=1, in which case the array elements are rearranged 
;   and the return is 1. uses the FHTR multi-threaded qsort like implementation

function calcSort, array, inPlace=inPlace

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; sanity checks
  if n_elements(array) gt 2e9 then message,'Error: Need to use 64 bit indices.'
  if size(array,/tname) ne 'LONG' and size(array,/tname) ne 'LONG64' then $
    message,'Error: Unsupported array variable type.'
  
  ; prepare input
  NumData = long(n_elements(array))
  
  if size(array,/tname) eq 'LONG'   then soName = 'int32'
  if size(array,/tname) eq 'LONG64' then soName = 'int64'
  
  ; inPlace? default method=2 is the FHTR (method=1 is glibc qsort)
  if ~keyword_set(inPlace) then begin
    method = 2L
    inds = lindgen(NumData)
  endif
  
  if keyword_set(inPlace) then begin
    method = 12L
    inds = [0]
  endif
  
  ret = Call_External('/n/home07/dnelson/idl/CalcMatch/CalcMatch_'+soName+'.so','CalcSort',$
                      NumData,array,inds,method,/CDECL)
                      
  if keyword_set(inPlace) then return,ret
  if ~keyword_set(inPlace) then return,inds
end

; calcMatch(): for two vectors A,B which contain only UNIQUE elements, find the intersection
;   where the return is such that A[ind1] = B[ind2]
;   note that like match.pro A[ind1] does not sample A in its (unsorted) order
; noSortA=1 : skip sort of A, and ind1 contains the permutation indices, but note 
;             that these ind1 will be overwritten by the subset match with B

pro calcMatch, A, B, ind1, ind2, count=count, noSortA=noSortA

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; sanity checks
  numA = long(n_elements(A))
  numB = long(n_elements(B))
  
  if numA gt 2e9 or numB gt 2e9 then message,'Error: Move to 64bit.'
  if size(A,/type) ne size(B,/type) then message,'Error: A and B have different var types.'
  if size(A,/tname) ne 'LONG' and size(A,/tname) ne 'LONG64' then $
    message,'Error: Unsupported variable type in A,B.'
    
  ; prepare inputs
  method = 1L

  ; prepare outputs
  if ~keyword_set(noSortA) then begin
    ind1 = lindgen(numA)
  endif else begin
    method = 11L ; skip A sort in external
    if n_elements(ind1) ne numA then message,'Error: Input ind1 size mismatch with A.'
  endelse
  
  ind2 = lindgen(numB)
  
  count = -1L
  
  ; pick library based on 32bit or 64bit IDs
  if size(A,/tname) eq 'LONG'   then soName = 'int32'
  if size(A,/tname) eq 'LONG64' then soName = 'int64'
    
  ret = Call_External('/n/home07/dnelson/idl/CalcMatch/CalcMatch_'+soName+'.so', 'CalcMatch', $
                      method,numA,numB,A,B,ind1,ind2,count,/CDECL)

  ; take index subsets
  if count eq -1L then message,'Error: Count unchanged.'
  
  if count gt 0 then begin
    if count lt numA then ind1 = ind1[0:count-1]
    if count lt numB then ind2 = ind2[0:count-1]
  endif else begin
    ind1 = -1
    ind2 = -1
  endelse

end

; calcMatchBlock(): gives the functionality of calcMatch() above but only the first array A is
; required. the second array B is loaded block by block such that only a small portion of it 
; need be stored in memory at any one time. the matches from each block are combined. the 
; array A is sorted once at the beginning and end, and B only once at the end, although nBlock 
; subsets of B are also sorted as intermediate steps, leading to at least 100% extra work.
; sP,partType,field : specify how to load B
; nBlocks : default to 10

pro calcMatchBlock, A, ind1, ind2, sP=sP, partType=partType, field=field, count=count, nBlocks=nBlocks

  compile_opt idl2, hidden, strictarr, strictarrsubs
    
  if ~keyword_set(sP) or n_elements(partType) eq 0 or ~keyword_set(field) then message,'Error.'
  
  numA = n_elements(A)
  
  ; determine blocking size
  if ~keyword_set(nBlocks) then nBlocks = 10
  
  nBlocks = fix(nBlocks)

  h = loadSnapshotHeader(sP=sP)
  numB = h.nPartTot[partTypeNum(partType)]
  
  blockSize = ceil(numB / nBlocks)
  
  ; make index return arrays
  ind1 = lonarr(numA) - 1
  ind2 = lonarr(numB) - 1
  
  ; make a copy of A and sort it (need sorted A and its indices)
  AA_sortinds = calcSort(A)
  if ~array_equal(AA_sortinds,sort(A)) then message,'Error: calcSort failure.'
  
  ; loop over blocks
  totCount = 0LL
  B_offset = 0LL
  
  for i=0,nBlocks-1 do begin
    ; load part of B
    indicesStart = i*blockSize
    indicesEnd   = (i+1)*blockSize-1
    
    if indicesStart lt 0 or indicesEnd ge numB then message,'Error: Bad indices.'
    if i eq nBlocks-1 then indicesEnd = numB - 1 ; make sure we read all the remaining
    
    indRange = [indicesStart, indicesEnd]
    
    B_subset = loadSnapshotSubset(sP=sP,partType=partType,field=field,indRange=indRange)
    
    ; run CalcMatch on B_subset, request we skip the sort of A
    sub_ind1 = AA_sortinds ; sort indices are input to CalcMatch, overwritten on output
    CalcMatch, A, B_subset, sub_ind1, sub_ind2, count=sub_count, /noSortA
    
    ;print,i,indicesStart,indicesEnd,sub_count
    B_offset += n_elements(B_subset)
    
    ; DEBUG
    ;match,A,B_subset,check_ind1,check_ind2,count=check_count
    ;if sub_count ne check_count then message,'Error.'
    ;if ~array_equal(sub_ind1,check_ind1) then message,'Error: ind1.'
    ;if ~array_equal(sub_ind2,check_ind2) then message,'Error: ind2.'  
    
    if sub_count eq 0 then continue
    
    ; stamp
    ind1[totCount : (totCount + sub_count - 1)] = sub_ind1
    ind2[totCount : (totCount + sub_count - 1)] = sub_ind2 + (B_offset-n_elements(B_subset))
    
    totCount += sub_count
  endfor
  
  AA_sortinds = !NULL
  B_subset    = !NULL
  sub_ind1    = !NULL
  sub_ind2    = !NULL
  
  if indicesEnd ne numB-1 then message,'Error: Failed to load all of B.'
  
  ; sanity check
  w = where(ind1 eq -1,count)
  if count gt 0 and min(w) lt totCount then message,'Error: Bad ind1.'
  w = where(ind2 eq -1,count)
  if count gt 0 and min(w) lt totCount then message,'Error: Bad ind2.'
  
  ; take overall subset
  if totCount gt 0 then begin
    if totCount lt numA then ind1 = ind1[0:totCount-1]
    if totCount lt numB then ind2 = ind2[0:totCount-1]
    
    ; now finished, while we have A loaded fix the ordering of ind1 with a global sort
    ind1 = ind1[ calcSort(A[ind1]) ]
  
    ; delete A and load full B
    A = !NULL
    B = loadSnapshotSubset(sP=sP,partType=partType,field=field)
    
    ; use the global sort of the matched elements of B to rearrange ind2
    ind2 = ind2[ calcSort(B[ind2]) ]
    B = !NULL
  endif else begin
    ; no matches
    ind1 = -1
    ind2 = -1
    A = !NULL ; for consistency
  endelse
  
  ; DEBUG: check the whole calculation
  ;match,A,B,check_ind1,check_ind2,count=check_count
  ;if totCount ne check_count then message,'Error.'
  ;if ~array_equal(ind1,check_ind1) then message,'Error: ind1.'
  ;if ~array_equal(ind2,check_ind2) then message,'Error: ind2.'  
   
  count = totCount
  
end

; mymatch(): produces identical output of match.pro, but no concatenated array

pro mymatch, A, B, ind1, ind2, count=count, retUnsorted=retUnsorted

  ; sort A and B if not already
  sort_A = sort(A)
  sort_B = sort(B)
  AA = A[sort_A]
  BB = B[sort_B]
  
  i = 0LL
  j = 0LL
  
  numA = n_elements(AA)
  numB = n_elements(BB)
  
  if numA gt 2e9 or numB gt 2e9 then message,'Error: Move to 64bit.'
  
  ind1 = lonarr(numA)
  ind2 = lonarr(numB)
  
  offset = 0LL
  
  while i lt numA and j lt numB do begin
    if AA[i] eq BB[j] then begin
      if keyword_set(retUnsorted) then begin
        ind1[offset] = i
        ind2[offset] = j
      endif else begin
        ; matched element
        ind1[offset] = sort_A[i]
        ind2[offset] = sort_B[j]
      endelse
      i += 1 & j += 1 & offset += 1
    endif else if AA[i] lt BB[j] then begin
      ; element mismatch, forward in A
      i += 1
    endif else begin
      ; element mismatch, forward in B
      j += 1
    endelse
  endwhile
  
  if offset gt 0 then begin
    ind1 = ind1[0:offset-1]
    ind2 = ind2[0:offset-1]
  endif else begin
    ind1 = -1
    ind2 = -1
  endelse
  
  count = offset
  
end

; testMatchBlock()

pro testMatchBlock

  ; generate some data
  num = 50
  A = shuffle([indgen(num)+30,indgen(num)+100])
  B = shuffle([indgen(num)+1000,indgen(num)+80])
  
  ; split B and match halves
  num = n_elements(B)
  
  b1 = b[0:num/2-1]
  b2 = b[num/2:num-1]
  match,a,b1,inda1,indb1,count=count1
  match,a,b2,inda2,indb2,count=count2
  
  count = count1+count2 ; number matched
  
  ; test
  mymatch,a,b1,inda1_test,indb1_test,count=count3,/retUnsorted
  mymatch,a,b2,inda2_test,indb2_test,count=count4,/retUnsorted
  
  inda_test = [inda1_test,inda2_test]
  indb_test = [indb1_test,indb2_test + n_elements(b1)]
  
  ; rearrange inda_test
  inda_fixed1 = (sort(a))[inda1_test]
  inda_fixed2 = (sort(a))[inda2_test]
  
  inda_fixed = [inda_fixed1,inda_fixed2]
  inda_fixed = inda_fixed[ sort(a[inda_fixed]) ]
  
  ; rearrange indb_test by the global sort of b
  indb_fixed1 = (sort(b1))[indb1_test]
  indb_fixed2 = (sort(b2))[indb2_test] + n_elements(b1)
  
  indb_fixed = [indb_fixed1,indb_fixed2] ; correct, just in wrong order
  indb_fixed = indb_fixed[ sort(b[indb_fixed]) ]
  
  ; do it again with a loop, have only full A to start
  ind1 = lonarr(n_elements(a))
  ind2 = lonarr(n_elements(b))
  
  offset = 0L
  b_offset = 0L
  
  for i=0,1 do begin
    ; load b subset
    if i eq 0 then b_subset = b[0:num/2-1]
    if i eq 1 then b_subset = b[num/2:num-1]
    
    ; match
    mymatch,a,b_subset,inda_sub,indb_sub,count=count_sub;,/retUnsorted
    b_offset += n_elements(b_subset)
    
    if count_sub eq 0 then continue
    
    ind1_sub = inda_sub
    ind2_sub = indb_sub + (b_offset-n_elements(b_subset))
    
    ind1[offset : offset+count_sub-1] = ind1_sub
    ind2[offset : offset+count_sub-1] = ind2_sub
    
    offset += count_sub
  endfor
  
  ; take matched subset
  if offset gt 0 then begin
    ind1 = ind1[0:offset-1]
    ind2 = ind2[0:offset-1]
  endif else begin
    ind1 = -1
    ind2 = -1
  endelse
  
  ; now finished, fix ind1 while we have A
  ind1 = ind1[ sort(a[ind1]) ]
  ; delete A and load full B, fix ind2
  ind2 = ind2[ sort(b[ind2]) ]
  
  ; match full for comparison
  match,a,b,Ia,Ib,count=count5

  print,'new:'
  print,array_equal(inda_fixed,Ia)
  print,array_equal(indb_fixed,Ib)
  print,'loop:'
  print,array_equal(ind1,Ia)
  print,array_equal(ind2,Ib)
  print,'counts:'
  print,count5,count1+count2,count3+count4
  
  stop
  
  ; OLD
  ; concatenate subindex arrays
  second_b_offset = count1
  
  inda       = [inda1,inda2]
  indb       = [indb1,indb2 + second_b_offset]
  
  ; now here we could delete full A and load full B, sort full B
  b_sorted    = b[sort(b)]
  b_subsorted = [b1(sort(b1)),b2(sort(b2))]
  
  b_subsortinds = [sort(b1),sort(b2) + second_b_offset]
  b_sortinds = sort(b)
  
  ; allocate arrays for our final solution
  inda_fixed = lonarr(num)
  indb_fixed = lonarr(num)
  
  ; walk through b_subsorted (and so inda,indb)
  if 0 then begin
  for i=0,count-1 do begin
    ; where in full b_sorted is this b_subsorted element?
    new_index = ( where(b_sorted eq b_subsorted[i]) )[0]
    new_index2 = b_subsorted[i] ; test, wrong
    
    ; move this inda index to the index corresponding to b_sorted (instead of b_subsorted)
    inda_fixed[new_index] = inda[i]
    indb_fixed[new_index] = indb[i]
    print,new_index,new_index2
  endfor
  endif ;0
  
end

; testMatching():

pro testMatching
  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function simParams, loadSnapshotSubset
  
  ; run each sort a few times and store mean runtime  
  nTests    = 4
  meanTimes = fltarr(nTests)
  
  nIter    = 1
  runTimes = fltarr(nIter)
  
  ; start memory monitor and timer
  print,'Loading...'
  
  mem = memory(/current) / 1024.0^3 ; GB
  time = systime(/sec) ; start timer (microsecond resolution)
  
  ; load some 64 bit data
  partType = 'gas'
  field    = 'ids'
  sP = simParams(res=128,run='tracer',redshift=2.0)
  ids_1 = loadSnapshotSubset(sP=sP,partType=partType,field=field)
  
  sP.snap -= 1
  ids_2 = loadSnapshotSubset(sP=sP,partType=partType,field=field)
  
  NumData = long(n_elements(ids_1))  
 
  ; report time and memory
  maxMem = memory(/highwater) / 1024.0^3 - mem
  curMem = memory(/current) / 1024.0^3
  print,'Load took: ['+string(systime(/sec)-time,format='(f5.1)')+'] sec, ['+$
    string(maxMem,format='(f5.2)')+' GB] max, ['+string(curMem,format='(f5.2)'),' GB] cur'
 
  ; run normal IDL match
  ;if 0 then begin
  print,'Running IDL match...' 
  
  for i=0,nIter-1 do begin  
    time = systime(/sec)
  
    match,ids_1,ids_2,ind1_a,ind2_a,count=count_a
  
    maxMem = memory(/highwater) / 1024.0^3 - mem
    curMem = memory(/current) / 1024.0^3
    runTimes[i] = systime(/sec) - time
    print,'['+str(i)+'] IDL match took: ['+string(runTimes[i],format='(f5.1)')+'] sec, ['+$
      string(maxMem,format='(f5.2)')+' GB] max, ['+string(curMem,format='(f5.2)'),' GB] cur'
  endfor
  
  meanTimes[0] = mean(runTimes)

  ; run pure IDL mymatch
  print,'Running pure IDL mymatch...'
  
  for i=0,nIter-1 do begin  
    time = systime(/sec)
  
    mymatch,ids_1,ids_2,ind1_b,ind2_b,count=count_b
  
    maxMem = memory(/highwater) / 1024.0^3 - mem
    curMem = memory(/current) / 1024.0^3
    runTimes[i] = systime(/sec) - time
    print,'['+str(i)+'] IDL mymatch took: ['+string(runTimes[i],format='(f5.1)')+'] sec, ['+$
      string(maxMem,format='(f5.2)')+' GB] max, ['+string(curMem,format='(f5.2)'),' GB] cur'
  endfor
  
  meanTimes[1] = mean(runTimes)
  ;endif ;0
  
  ; run CalcMatch one-pass
  print,'Running CalcMatch one-pass...'
  
  for i=0,nIter-1 do begin  
    time = systime(/sec)
    ids_2 = loadSnapshotSubset(sP=sP,partType=partType,field=field)
    CalcMatch,ids_1,ids_2,ind1_c,ind2_c,count=count_c
  
    maxMem = memory(/highwater) / 1024.0^3 - mem
    curMem = memory(/current) / 1024.0^3
    extMem = (n_elements(ids_1) + n_elements(ids_2)) * 4 / 1024.0^3
    runTimes[i] = systime(/sec) - time
    print,'['+str(i)+'] CalcMatch took: ['+string(runTimes[i],format='(f5.1)')+'] sec, ['+$
      string(maxMem,format='(f5.2)')+' GB] max, ['+string(extMem,format='(f5.2)')+' GB] ext, ['+$
      string(curMem,format='(f5.2)'),' GB] cur'
  endfor
  
  meanTimes[2] = mean(runTimes)
  
  ; run CalcMatch blocked
  print,'Running CalcMatch blocked...'
  
  for i=0,nIter-1 do begin  
    time = systime(/sec)
  
    CalcMatchBlock,ids_1,ind1_d,ind2_d,sP=sP,partType=partType,field=field,count=count_d
  
    maxMem = memory(/highwater) / 1024.0^3 - mem
    curMem = memory(/current) / 1024.0^3
    extMem = (n_elements(ids_1) + n_elements(ids_2)/10) * 4 / 1024.0^3
    runTimes[i] = systime(/sec) - time
    print,'['+str(i)+'] CalcMatchBlock took: ['+string(runTimes[i],format='(f5.1)')+'] sec, ['+$
      string(maxMem,format='(f5.2)')+' GB] max, ['+string(extMem,format='(f5.2)')+' GB] ext, ['+$
      string(curMem,format='(f5.2)'),' GB] cur'
  endfor
  
  meanTimes[3] = mean(runTimes)
  
  ; verify
  strings = ['IDL match  ','IDL mymatch','CalcMatch  ','CalcBlock  ']
  foreach sName,strings,i do print,sName+' : ['+string(meanTimes[i],format='(f5.1)')+'] sec.'
  
  if count_a ne count_b then message,'Error.'
  if count_a ne count_c then message,'Error.'
  if count_a ne count_d then message,'Error.'
  if ~array_equal(ind1_a,ind1_b) then message,'Error 1b.'
  if ~array_equal(ind2_a,ind2_b) then message,'Error 2b.'
  if ~array_equal(ind1_a,ind1_c) then message,'Error 1c.'
  if ~array_equal(ind2_a,ind2_c) then message,'Error 2c.'
  if ~array_equal(ind1_a,ind1_d) then message,'Error 1d.'
  if ~array_equal(ind2_a,ind2_d) then message,'Error 2d.'
  stop
  
end

pro monteCarloTestMatching

  nIter = 10
  
  ; load full data set
  partType = 'tracermc'
  field    = 'tracerids'
  sP = simParams(res=128,run='tracer',redshift=2.0)
  ids_1_full = loadSnapshotSubset(sP=sP,partType=partType,field=field)
  
  sP.snap -= 1
  ids_2_full = loadSnapshotSubset(sP=sP,partType=partType,field=field)
  
  ; loop
  for i=0,nIter-1 do begin
    ; make subsets, shuffle
    startInd = round( randomu(seed,1) * 0.5 * n_elements(ids_1) ) > 0
    endInd   = round( randomu(seed,1) * 2.0 * n_elements(ids_1) ) $
               > startInd + 1000 < n_elements(ids_1)-1   
    
    ids_1 = shuffle(ids_1_full[startInd:endInd])
  
    startInd = round( randomu(seed,1) * 0.5 * n_elements(ids_1) ) > 0
    endInd   = round( randomu(seed,1) * 2.0 * n_elements(ids_1) ) $
               > startInd + 1000 < n_elements(ids_1)-1  

    ids_2 = shuffle(ids_2_full[startInd:endInd])       
  
    ; run CalcMatch, CalcMatchBlock
    match,ids_1,ids_2,ind1_a,ind2_a,count=count_a
    
    CalcMatch,ids_1,ids_2,ind1_c,ind2_c,count=count_c
    
    ; verify
    if count_a ne count_c then message,'Error.'
    if ~array_equal(ind1_a,ind1_c) then message,'Error 1c.'
    if ~array_equal(ind2_a,ind2_c) then message,'Error 2c.'
    
    ; run CalcMatch with ids_2_full, CalcMatchBlock
    CalcMatch,ids_1,ids_2_full,ind1_a,ind2_a,count=count_a
    
    CalcMatchBlock,ids_1,ind1_c,ind2_c,sP=sP,partType=partType,field=field,count=count_c
    
    ; verify
    if count_a ne count_c then message,'Error.'
    if ~array_equal(ind1_a,ind1_c) then message,'Error 1c.'
    if ~array_equal(ind2_a,ind2_c) then message,'Error 2c.'
    
    print,i,startInd,endInd,'Passed.'
  endfor

end

; testSorting(): testing sort algorithms, IDL vs various C implementations, single vs multi threaded

pro testSorting
  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function simParams, loadSnapshotSubset
    
  ; run each sort a few times and store mean runtime
  nTests    = 4
  meanTimes = fltarr(nTests)
  
  nIter    = 3
  runTimes = fltarr(nIter)

  ; disable IDL threadpool
  ;CPU, tpool_nthreads = 1
  
  ; start memory monitor and timer
  print,'Loading...'
  
  mem = memory(/current) / 1024.0^3 ; GB
  time = systime(/sec) ; start timer (microsecond resolution)
  
  ; load some 64 bit data
  sP = simParams(res=256,run='tracer',redshift=2.0)
  ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
  
  NumData = long(n_elements(ids))  
  
  ; report time and memory
  maxMem = memory(/highwater) / 1024.0^3 - mem
  curMem = memory(/current) / 1024.0^3
  print,'Load took: ['+string(systime(/sec)-time,format='(f5.1)')+'] sec, ['+$
    string(maxMem,format='(f5.2)')+' GB] max, ['+string(curMem,format='(f5.2)'),' GB] cur'
  
  ; run normal IDL sort
  print,'Running IDL sort...'
  
  inds_0 = lonarr(NumData)
  
  for i=0,nIter-1 do begin
    time = systime(/sec) ; start timer
    inds_0 = sort(ids)
  
    maxMem = memory(/highwater) / 1024.0^3 - mem
    curMem = memory(/current) / 1024.0^3
    runTimes[i] = systime(/sec) - time
    print,'['+str(i)+'] IDL sort took: ['+string(runTimes[i],format='(f5.1)')+'] sec, ['+$
      string(maxMem,format='(f5.2)')+' GB] max, ['+string(curMem,format='(f5.2)'),' GB] cur'
  endfor
  
  meanTimes[0] = mean(runTimes)
  
  ; run external qsort (glibc, single threaded)
  print,'Running glibc qsort (single threaded)...'
  method = 1L

  inds_1 = lindgen(NumData)
  
  for i=0,nIter-1 do begin
    time = systime(/sec) ; start timer
  
    ret = Call_External('/n/home07/dnelson/idl/CalcMatch/CalcMatch_int64.so','CalcSort',$
                        NumData,ids,inds_1,method,/CDECL)
  
    maxMem = memory(/highwater) / 1024.0^3 - mem
    extMem = 0 ;(bytesMyInt * NumData + bytesMyInd * NumData) / 1024.0^3
    curMem = memory(/current) / 1024.0^3
    runTimes[i] = systime(/sec) - time
    print,'['+str(i)+'] glibc qsort took: ['+string(runTimes[i],format='(f5.1)')+'] sec, ['+$
      string(maxMem,format='(f5.2)')+' GB] max, ['+string(extMem,format='(f5.2)')+' GB] ext, ['+$
      string(curMem,format='(f5.2)'),' GB] cur'
  endfor
  
  meanTimes[1] = mean(runTimes)
  
  ; run external qsort (threaded, fhtr generic)
  print,'Running fhtr qsort (multi-threaded)...'
  method = 2L
  
  inds_2 = lindgen(NumData)
  
  for i=0,nIter-1 do begin
    time = systime(/sec) ; start timer
  
    ret = Call_External('/n/home07/dnelson/idl/CalcMatch/CalcMatch_int64.so','CalcSort',$
                        NumData,ids,inds_2,method,/CDECL)
  
    maxMem = memory(/highwater) / 1024.0^3 - mem
    extMem = 0 ;(bytesMyInt * NumData + bytesMyInd * NumData) / 1024.0^3
    curMem = memory(/current) / 1024.0^3
    runTimes[i] = systime(/sec) - time
    print,'['+str(i)+'] fhtr-generic qsort took: ['+string(runTimes[i],format='(f5.1)')+'] sec, ['+$
      string(maxMem,format='(f5.2)')+' GB] max, ['+string(extMem,format='(f5.2)')+' GB] ext, ['+$
      string(curMem,format='(f5.2)'),' GB] cur'
  endfor
  
  meanTimes[2] = mean(runTimes)
  
  ; run inplace
  print,'Running fhtr qsort inplace (multi-threaded)...'
  method = 12L
  inds_3 = [0]
  
  for i=0,nIter-1 do begin
    dataIn = ids
    time = systime(/sec) ; start timer
  
    ret = Call_External('/n/home07/dnelson/idl/CalcMatch/CalcMatch_int64.so','CalcSort',$
                        NumData,dataIn,inds_3,method,/CDECL)
  
    maxMem = memory(/highwater) / 1024.0^3 - mem
    extMem = 0 ; O(logn) at most since inplace
    curMem = memory(/current) / 1024.0^3
    runTimes[i] = systime(/sec) - time
    print,'['+str(i)+'] fhtr-inplace qsort took: ['+string(runTimes[i],format='(f5.1)')+'] sec, ['+$
      string(maxMem,format='(f5.2)')+' GB] max, ['+string(extMem,format='(f5.2)')+' GB] ext, ['+$
      string(curMem,format='(f5.2)'),' GB] cur'
  endfor
  
  meanTimes[3] = mean(runTimes)
  
  ; verify sorts agree
  strings = ['IDL   ','glibc ','fhtr1 ','fhtrIP']
  foreach sName,strings,i do print,sName+' sort: ['+string(meanTimes[i],format='(f5.1)')+'] sec.'
  
  if ~array_equal(inds_0,inds_1) then message,'Error: Fail 1.'
  if ~array_equal(inds_0,inds_2) then message,'Error: Fail 2.'
  stop
end
