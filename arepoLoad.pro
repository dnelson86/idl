; arepoLoad.pro
; dnelson
; 7/5/10
;
; loading functions for various Arepo outputs

; loadVoronoi2D(): load "voronoi_mesh_" output
;
; inputs:  fBase, nDims, i
; outputs: [nGas,nEl,nDt] returned
;          nEdges,nEdgesOffset,edgeList,xyzEdges are overwritten with data

function loadVoronoi2D, fBase, nDims, i, $
                        nEdges, nEdgesOffset, edgeList, xyzEdges

  ;set filename
  if (str(i) eq 'none') then begin
    f = fBase
  endif else begin  
    ext = string(i,format='(i3.3)')
    f = fBase + ext
  endelse

  openr,1,f
    ;read header
    nGas   = 0L
    nEl    = 0L
    Ndt    = 0L ;nEdgePoints
    
    readu,1,nGas,nEl,Ndt
    ;print,i,nGas,nEl,Ndt
    
    ;read mesh
    nEdges       = lonarr(nGas)
    nEdgesOffset = lonarr(nGas)
    edgeList     = lonarr(nEl)
    xyzEdges     = fltarr(nDims,Ndt)
    
    readu,1,nEdges
    readu,1,NedgesOffset
    readu,1,edgeList
    readu,1,xyzEdges
    
  close,1

  return,[nGas,nEl,nDt]
  
end

; loadDensityField2D():
;

function loadDensityField2D, fBase, i, dens

  ;set filename
  ext = string(i,format='(i3.3)')
  f = fBase + ext

  openr,1,f
    ;read header
    nPixelsX = 0L
    nPixelsY = 0L
    
    readu,1,nPixelsX
    readu,1,nPixelsY

    ;read density field
    dens= fltarr(nPixelsY, nPixelsX)
    readu,1,dens
    
  close,1

  dens = transpose(dens)

  return,[nPixelsX,nPixelsY]

end

; getDensityMinMax(): return minimum and maximum for the whole sequence of density fields for output display scaling in a movie
;

function getDensityMinMax, densBase, nSnapshots

  rMinMax = fltarr(2)

  for i=0,nSnapshots,1 do begin 
  
    nPixelsXY = loadDensityField2D(densBase,i,dens)
    tMinMax = minmax(dens)
    
    rMinMax[0] = min([rMinMax[0],tMinMax[0]])
    rMinMax[1] = max([rMinMax[1],tMinMax[1]])
  
  endfor
  
  return, rMinMax

end

; loadSnapshot():
;
; inputs:  fBase, i
; outputs: h=header returned
;          pos,vel,id,mass,u,rho,hsml overwritten with data

function loadSnapshot, fBase, i, $
                       pos, vel, id, mass, u, rho, hsml

  ;set filename
  if (str(i) eq 'none') then begin
    f = fBase
  endif else begin  
    ext = string(i,format='(i3.3)')
    f = fBase + ext
  endelse
  
  ; check HDF5 format
  if h5f_is_hdf5(f) then begin
    print, 'snap is HDF5.'
    return,0
  endif

  ; header structure
  bytesLeft     = 136
  h = { header,                           $
        nPartTot:     lonarr(6),          $
        massTable:    dblarr(6),          $
        time:         0.0D,               $
        redshift:     0.0D,               $
        flagSFR:      0L,                 $
        flagFeedback: 0L,                 $
        npartall:     lonarr(6),          $
        la:           intarr(bytesLeft/2) $
      }
  
  openr,1,f,/f77_unformatted
    ;read header
    readu,1,h

    nGas = h.nPartTot[0]
    
    ;read blocks for all gas particles
    pos = fltarr(3,nGas)
    vel = fltarr(3,nGas)
    id  = lonarr(nGas)
    readu,1, pos
    readu,1, vel
    readu,1, id
    
    ; masses for variable mass particles (usually gas+stars)
    ind=where((h.nPartTot gt 0) and (h.massTable eq 0))
    
    if ind[0] ne -1 then begin
      nMass = total(h.nPartTot[ind])
      mass = fltarr(nMass)
      readu,1,mass
    endif
    
    ; internal energy per unit mass
    u=fltarr(nGas)
    readu,1,u
    
    ; comoving gas density
    if (str(i) ne 'none') then begin ;rho is added in snapshots but doesn't exist in ICs
      rho=fltarr(nGas)
      readu,1,rho
    
      if h.flagSFR gt 0 then begin
        ; gas electron abundance relative to hydrogen
        Nelec=fltarr(nGas)
        readu,1,nElec
        ; neutral hydrogen abundance relative to hydrogen
        NH0=fltarr(nGas)
        readu,1,NH0
      endif
      
      ; SPH smoothing length
      hsml=fltarr(nGas)
      readu,1,hsml  
   endif                  
  
  close,1
  
  return, h
end

; loadSnapshotHDF5(): GAS ONLY
;
; inputs:  fBase, i
; outputs: h=header returned
;          pos,vel,id,mass,u,rho overwritten with data

function loadSnapshotHDF5, fBase, i, $
                           pos, vel, id, mass, u, rho, vol, hsml

  ; set filename
  if (i ne 'none') then begin
    ext = string(i,format='(i3.3)')
    f = fBase + ext + '.hdf5'
  endif else begin
    f = fBase
  endelse
  
  ; check HDF5 format
  if not h5f_is_hdf5(f) then begin
    print, 'error: snap is not HDF5.'
    return,0
  endif
    
  ; recursively load whole hdf5 file into a structure
  s = h5_parse(f,/read_data)
  
  ; header structure
  h = { headerHDF                                                      ,$
        nPartThisFile       : s.header.numPart_ThisFile._DATA          ,$
        nPartTot            : s.header.numPart_Total._DATA             ,$
        nPartTotHighword    : s.header.numPart_Total_Highword._DATA    ,$
        massTable           : s.header.massTable._DATA                 ,$
        time                : s.header.time._DATA                      ,$
        redshift            : s.header.redshift._DATA                  ,$
        boxSize             : s.header.boxSize._DATA                   ,$
        numFilesPerSnapshot : s.header.numFilesPerSnapshot._DATA       ,$
        Omega0              : s.header.Omega0._DATA                    ,$
        OmegaLambda         : s.header.OmegaLambda._DATA               ,$
        hubbleParam         : s.header.hubbleParam._DATA               ,$
        flagSFR             : s.header.flag_SFR._DATA                  ,$
        flagCooling         : s.header.flag_Cooling._DATA              ,$
        flagStellarAge      : s.header.flag_StellarAge._DATA           ,$
        flagMetals          : s.header.flag_Metals._DATA               ,$
        flagFeedback        : s.header.flag_Feedback._DATA             ,$
        flagDoublePrecision : s.header.flag_DoublePrecision._DATA      $
        ;compositionVecLen   : s.header.Composition_Vector_Length._DATA  $ ;not in Gadget3 output 
      }
  
  ; particle type 0 (gas)
  pos  = s.parttype0.coordinates._DATA
  vel  = s.parttype0.velocities._DATA
  id   = s.parttype0.particleIDs._DATA
  mass = s.parttype0.masses._DATA
  u    = s.parttype0.internalEnergy._DATA
  rho  = s.parttype0.density._DATA
  vol  = s.parttype0.volume._DATA
  hsml = s.parttype0.smoothingLength._DATA
  if (total(tag_names(s.parttype0) eq 'VOLUME') ge 1.0) then $
    vol  = s.parttype0.volume._DATA

  return, h
end
