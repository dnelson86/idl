; tracers.pro
; dev for tracer particles
; dnelson jan.2012

function interpNN, pos_src, pos_dest, val_src

    ; something fancy (linear or cubic?)
    ;dens_interp = interpol(dens_gas,pos_gas,pos_tracer)
    ;res = dens_interp - dens_tracer
    
    ; nearest neighbor interpolation by indices
    ;inds = round(pos_gas / 20.0 * n_elements(dens_tracer))
    
    ; nearest neighbor interpolation (manual)
    val_interp = fltarr(n_elements(pos_src))
    
    for i=0,n_elements(pos_src)-1 do begin
      ind = where(abs(pos_dest-pos_src[i]) eq min(abs(pos_dest-pos_src[i])),count)
      
      if (count eq 0) then begin
        print,'WARNING'
        stop
      endif
      val_interp[i] = val_src[ind[0]]
    endfor
    
  return, val_interp
end

; calcHSML(): use CalcHSML external C-routine for the tree, neighbor searchings, and smoothing lengths

function calcHSML, Pos, ndims=ndims, nNGB=nNGB, boxSize=boxSize

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

; calcNN(): use CalcNN external C-routine for the tree and neighbor search
;           return the index of Pos_SrcTargs closest to each of Pos_SrcOrigs

function calcNN, Pos_SrcTargs, Pos_SrcOrigs, boxSize=boxSize, ndims=ndims

  ; prepare inputs
  n_srcTargs = long( n_elements(Pos_SrcTargs[0,*]) )
  n_srcOrigs = long( n_elements(Pos_SrcOrigs[0,*]) )
  
  boxSize   = float(boxSize)
  
  ind_out = lonarr(n_srcOrigs)
  
  ; call CalcNN
  libName = '/n/home07/dnelson/idl/CalcNN/CalcNN_'+str(ndims)+'D.so'
  ret = Call_External(libName, 'CalcNN', $
                      n_srcTargs,n_srcOrigs,Pos_SrcTargs,Pos_SrcOrigs,boxSize,ind_out, $
                      /CDECL)
    
  return, ind_out
  
end

; calcSphMap(): use CalcSphMap external C-routine to calculate a map of projected densities
;               with the sph spline kernel

function calcSphMap, pos, hsml, mass, axes=axes, boxSize=boxSize, boxCen=boxCen, $
                     nPixels=nPixels, ndims=ndims

  ; prepare inputs
  npos = (size(pos))[2]
  
  boxSize  = float(boxSize)
  boxCen   = float(boxCen)
  nPixels  = fix(nPixels)
  axes     = fix(axes)
  
  mode     = fix(1)
  periodic = fix(1)
  
  ; check inputs
  if (n_elements(boxSize) ne 3) then stop
  if (n_elements(boxCen)  ne 3) then stop
  if (n_elements(nPixels) ne 2) then stop
  if (n_elements(axes)    ne 2) then stop
  
  if (mode ne 1 and mode ne 2 and mode ne 3) then stop
  if (periodic ne 0 and periodic ne 1) then stop
  
  if (axes[0] ne 0 and axes[0] ne 1 and axes[0] ne 2) then stop
  if (axes[1] ne 0 and axes[1] ne 1 and axes[1] ne 2) then stop
  
  ; make return
  grid_out = fltarr(nPixels[0],nPixels[1])
  
  ; call CalcSphMap
  libName = '/n/home07/dnelson/idl/CalcSphMap/CalcSphMap_'+str(ndims)+'D.so'
  ret = Call_External(libName, 'CalcSphMap', $
                      npos,pos,hsml,mass,grid_out,boxSize[0],boxSize[1],boxSize[2],$
                      boxCen[0],boxCen[1],boxCen[2],axes[0],axes[1],nPixels[0],nPixels[1],$
                      mode,periodic,/CDECL)

  return, grid_out
  
end

; estimateDensityTophat(): spatial density estimator for an input position tuple array and
;                          mass array for some particle type
;   do density calculation by calculating smoothing lengths for all the particles
;   

function estimateDensityTophat, pos, mass=mass, ndims=ndims, nNGB=nNGB, boxSize=boxSize

  if (not keyword_set(nNGB) or not keyword_set(ndims)) then stop

  ; calculate smoothing lengths
  hsml_out = calcHSML(pos,ndims=ndims,nNGB=nNGB,boxSize=boxSize)
                      
  ; estimate densities on eval_pos using hsml
  if (ndims eq 1) then $
    eval_dens = nNGB / (2.0 * hsml_out)
    
  if (ndims eq 2) then $
    eval_dens = nNGB / (!pi * hsml_out^2.0)
    
  if (ndims eq 3) then $
    eval_dens = nNGB / (4.0*!pi/3.0 * hsml_out^3.0)
  
  ; add in mass
  if keyword_set(mass) then $
    eval_dens *= mass[0]
  
  return,eval_dens

end
