; tracers.pro
; dnelson
; jan 2012
;
; dev for tracer particles

@helper
@cosmoUtil
@cosmoLoad
@arepoVis2D

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

; estimateDensityTophat(): spatial density estimator for an input position tuple array and
;                          mass array for some particle type
;   do density calculation by calculating smoothing lengths for all the particles
;   use CalcHSML external C-routine for the tree, neighbor searchings, and smoothing lengths

function estimateDensityTophat, pos, mass, ndims=ndims, nNGB=nNGB

  if (not keyword_set(nNGB) or not keyword_set(ndims)) then begin
    print,'Error: Must specify both nNGB and ndims.'
    return,0
  endif

  ; config
  sz_pos = size(pos)
  npos   = sz_pos[2]
  
  ; arrays
  eval_dens = fltarr(npos)

  ; prepare inputs
  NumPart = long(npos)
  ;Pos     = pos
  ;Mass    = mass
  
  DesNumNgb    = long(nNGB) ; number of neighbors to use
  DesNumNgbDev = long(0)
  BoxSize      = 0.0
  HsmlGuess    = float(1.0)
  Softening    = float(1.0)
  
  hsml_out = fltarr(NumPart)
  
  ; call CalcHSML
  libName = '/n/home07/dnelson/idl/CalcHSML/CalcHSML_'+str(ndims)+'D.so'
  ret = Call_External(libName, 'CalcHSML', $
                      NumPart,Pos,Mass,DesNumNgb,DesNumNgbDev,BoxSize,HsmlGuess,Softening,hsml_out, $
                      /CDECL)

  ; estimate densities on eval_pos using hsml
  if (ndims eq 1) then $
    eval_dens = DesNumNgb / hsml_out
    
  if (ndims eq 2) then $
    eval_dens = DesNumNgb / (!pi * hsml_out^2.0)
    
  if (ndims eq 3) then $
    eval_dens = DesNumNgb / (4.0*!pi/3.0 * hsml_out^3.0)
  
  ; add in mass
  eval_dens *= mass[0]
  
  return,eval_dens

end

@tracersDisks
@tracersShocktube