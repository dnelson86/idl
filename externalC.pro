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