; arepoLoad.pro
; loading/writing functions for various Arepo outputs, snapshots, ICs
; dnelson feb.2012

; loadVoronoiMesh(): load "voronoi_mesh_" output

function loadVoronoiMesh, fileBase, m, nDims

  ;set filename
  if (str(m) eq 'none') then begin
    f = fileBase
  endif else begin  
    ext = string(m,format='(i3.3)')
    f = fileBase + ext
  endelse

  if (not file_test(f)) then begin
    print, 'ERROR: voronoi mesh [' + str(m) + '] at ' + fileBase + ' does not exist!'
    stop
  endif
  
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

  r = {nGas:nGas,nEl:nEl,nDt:nDt,nEdges:nEdges,$
       nEdgesOffset:nEdgesOffset,edgeList:edgeList,xyzEdges:xyzEdges}
  return,r
  
end

; loadDensityField(): load "density_field_" or "proj_density_field" output file

function loadDensityField, fileBase, m, axes=axes

  ; set filename
  if (str(m) eq 'none') then begin
    ext = ''
    f = fileBase
  endif else begin
    ext = string(m,format='(I3.3)')
    f = fileBase + ext
  endelse
  
  ; check for subdirectories
  if file_test(fileBase + 'proj_density_field/',/directory) then $
    fileBase = fileBase + 'proj_density_field/'
   
  if file_test(fileBase+'density_field_'+ext) then $
    f = fileBase + 'density_field_' + ext
  if file_test(fileBase+'proj_density_field_'+ext) then $
    f = fileBase + 'proj_density_field_' + ext

  if (not file_test(f)) then begin
    ; check for multiple density fields (projections along different axes)
    f = fileBase + 'proj_density_field_' + ext + '.' + str(axes) + '.dat'
      
    if (not file_test(f)) then begin
      print,'ERROR: snap ['+str(m)+'] axes ['+str(axes)+'] at ' + fileBase + ' not found!'
      stop
    endif
  endif

  ; load
  openr,1,f
    ;read header
    nPixelsX = 0L
    nPixelsY = 0L
    
    readu,1,nPixelsX
    readu,1,nPixelsY

    ;read density field
    dens = fltarr(nPixelsY, nPixelsX)
    readu,1,dens
    
    ;read temperature field
    temp = fltarr(nPixelsY, nPixelsX)
    ;readu,1,temp
    
  close,1

  dens = transpose(dens)
  temp = transpose(temp)
  
  r = {nPixelsXY:[nPixelsX,nPixelsY],dens:dens,temp:temp}
  return,r

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

; loadSnapshotOld(): load an old (non-HDF5) snapshot

function loadSnapshotOld, fileBase, m
  
  ;set filename
  if (str(m) eq 'none') then begin
    f = fileBase
  endif else begin  
    ext = string(m,format='(i3.3)')
    f = fileBase + ext
  endelse
  
  ; check HDF5 format
  if h5f_is_hdf5(f) then message, 'Error: snap is HDF5.'

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
      
  openr,lun,f,/get_lun,/f77_unformatted
    ; read header
    readu, lun, h

    nGas = h.nPartTot[0]
    nTot = total(h.nPartTot,/int)
    
    ; masses for variable mass particles (usually gas+stars)
    ind = where((h.nPartTot gt 0) and (h.massTable eq 0), count)
    
    nMass = 0
    if count gt 0 then nMass = total(h.nPartTot[ind],/int)
    
    ; create structure for fields, all particles
    pos  = fltarr(3,nTot)
    vel  = fltarr(3,nTot)
    ids  = lonarr(nTot)
    if nMass gt 0 then mass = fltarr(nMass)
    if nGas  gt 0 then u    =  fltarr(nGas)
            
    ; read gas values that are always present
    readu, lun, pos
    readu, lun, vel
    readu, lun, ids
    if nMass gt 0 then readu, lun, mass
    if nGas  gt 0 then readu, lun, u
    
    ; fields that are written to snapshots but do not exist in ICs
    if ~EOF(lun) then begin
    
      rho   = fltarr(nGas)
      nelec = fltarr(nGas)
      nh0   = fltarr(nGas)
      hsml  = fltarr(nGas)
               
      readu, lun, rho ; density (code units)
    
      if h.flagSFR gt 0 then begin
        readu, lun, nelec ; gas electron abundance relative to hydrogen
        readu, lun, nh0   ; neutral hydrogen abundance relative to hydrogen
      endif
      
      readu, lun, hsml ; SPH smoothing length
    endif
    
  close,lun
  
  ; split by particle type
  offset = 0L
  offsetMass = 0L
  
  for i=0,n_elements(h.nPartTot)-1 do begin
    if h.nPartTot[i] eq 0 then continue
    
    partType = { pos : pos[*,offset : offset + h.nPartTot[i]-1] ,$
                 vel : vel[*,offset : offset + h.nPartTot[i]-1] ,$
                 ids : ids[  offset : offset + h.nPartTot[i]-1] }
         
    ; add masses if massTable[partType] is nonzero
    if h.massTable[i] eq 0 then begin
      partType = mod_struct( partType, 'mass', mass[offsetMass : offsetMass + h.nPartTot[i]-1] )
      offsetMass += h.nPartTot[i]
    endif
    
    ; if gas, add extra fields
    if i eq 0 then begin
      if n_elements(u)     gt 0 then partType = mod_struct( partType, 'u', u )
      if n_elements(rho)   gt 0 then partType = mod_struct( partType, 'rho', rho )
      if n_elements(nelec) gt 0 then partType = mod_struct( partType, 'nelec', nelec )
      if n_elements(nh0)   gt 0 then partType = mod_struct( partType, 'nh0', nh0 )
      if n_elements(hsml)  gt 0 then partType = mod_struct( partType, 'hsml', hsml )
    endif
    
    r = mod_struct( r, 'PARTTYPE'+str(i), partType )
    offset += h.nPartTot[i]
  endfor
  
  r = mod_struct( r, 'Header', h )
  
  return, r
end

; writeICFile(): write old Gadget format IC file with gas particles and tracers
;                each partX struct should contain {id,pos,vel,mass,u} in the usual format

pro writeICFile, fOut, part0=part0, part1=part1, part2=part2, part3=part3, massarr=massarr, $
                 longIDs=longIDs, doublePrecision=doublePrecision

  ; arrays
  pos  = []
  vel  = []
  id   = []
  mass = []
  u    = []

  ; type checking (type code 4 = FLOAT precision, 3 = LONG)
  if keyword_set(massarr) then $
    if ( (size(massarr))[2] ne 4 ) then print,'WARNING: massarr type.'
    
  valTypeCode = 4 ; FLOAT
  idTypeCode  = 3 ; LONG32
  
  if keyword_set(doublePrecision) then valTypeCode = 5 ;DOUBLE
  if keyword_set(longIDs) then idTypeCode = 14 ;LONG64

  ; create header
  npart    = lonarr(6)  
  if not keyword_set(massarr) then massarr  = dblarr(6)
  npartall = lonarr(6)

  ; add to particle counts and concat arrays
  if keyword_set(part0) then begin
    ; GAS
    npart[0]    = n_elements(part0.id)
    npartall[0] = n_elements(part0.id)
    
    ; check typing
    if (size(part0.pos))[3] ne valTypeCode then print,'WARNING: part0 pos type.'
    if (size(part0.vel))[3] ne valTypeCode then print,'WARNING: part0 vel type.'
    if (size(part0.id))[2] ne idTypeCode then print,'WARNING: part0 id type.'
    if (size(part0.u))[2] ne valTypeCode then print,'WARNING: part0 u type.'    
    
    pos  = [[pos], [part0.pos]]
    vel  = [[vel], [part0.vel]]
    id   = [id,    part0.id]
    if (massarr[0] eq 0.0) then begin ; if massTable[partType]=0 then expect mass block in ICs
      if (size(part0.mass))[2] ne valTypeCode then print,'WARNING: part0 mass type.'
      mass = [mass,  part0.mass]
    endif    
    u    = [u,     part0.u]
  endif

  if keyword_set(part1) then begin
    ; DM
    npart[1]    = n_elements(part1.id)
    npartall[1] = n_elements(part1.id)
    
    pos  = [[pos], [part1.pos]]
    vel  = [[vel], [part1.vel]]
    id   = [id,    part1.id]
    mass = [mass,  part1.mass]
    u    = [u,     part1.u]
  endif
  
  if keyword_set(part2) then begin
    print,'Update the script.'
    stop
    ; TRACER (old partType)
  endif
  
  if keyword_set(part3) then begin
    ; TRACER
    npart[3]    = n_elements(part3.id)
    npartall[3] = n_elements(part3.id)
    
    ; check typing
    if (size(part3.pos))[3] ne valTypeCode then print,'WARNING: part3 pos type.'
    if (size(part3.vel))[3] ne valTypeCode then print,'WARNING: part3 vel type.'
    if (size(part3.id))[2] ne idTypeCode then print,'WARNING: part3 id type.'
    
    pos  = [[pos], [part3.pos]]
    vel  = [[vel], [part3.vel]]
    id   = [id,    part3.id]
    ; mass not expected in input for tracer (after my change to not output mass in snapshots)
    ; note: this is different in different versions of Arepo right now (true in gasSphere, convFlow)
    ;if (massarr[3] eq 0.0) then begin ; if massTable[partType]=3 then expect mass block in ICs
    ;  if (size(part3.mass))[2] ne valTypeCode then print,'WARNING: part3 mass type.'
    ;  mass = [mass,  part3.mass]
    ;endif
    ;u    = [u,     part3.u] ; u not expected in input for tracer
  endif

  ; double precision?
  if keyword_set(doublePrecision) then begin
    ; force double
    pos = double(pos)
    vel = double(vel)
    if (n_elements(mass) gt 0) then mass = double(mass)
    u = double(u)
  endif else begin
    ; force single
    pos = float(pos)
    vel = float(vel)
    if (n_elements(mass) gt 0) then mass = float(mass)
    u = float(u)
  endelse
  
  ; long (64bit) ids?
  if keyword_set(longIDs) then begin
    id = long64(id)
  endif else begin
    id = long(id)
  endelse

  ; header
  time          = 0.0D
  redshift      = 0.0D
  flag_sfr      = 0L
  flag_feedback = 0L
  
  bytesleft = 136
  la        = intarr(bytesleft/2)

  ; write IC file
  openw,1,fOut,/f77_unformatted
  writeu,1,npart,double(massarr),time,redshift,flag_sfr,flag_feedback,npartall,la

  writeu,1, pos
  writeu,1, vel
  writeu,1, id
  if (n_elements(mass) gt 0) then $
    writeu,1, mass    
  writeu,1, u
  close,1
  
  print,'wrote ',fOut
  print,massarr
  print,n_elements(pos[0,*]),n_elements(vel[0,*]),n_elements(id),n_elements(mass),n_elements(u)
end

; writeICFileHDF5(): GAS ONLY
pro writeICFileHDF5, fOut, boxSize, pos, vel, id, massOrDens, u

  ; load hdf5 template
  templatePath = '/n/home07/dnelson/make.ics/ArepoTemplate.hdf5'
  
  s = h5_parse(templatePath, /read)

  ; modify base
  s._NAME    = fOut
  s._FILE    = fOut
  s._COMMENT = "dnelson IC gen"
  
  ; modify HEADER
  s.HEADER._FILE                      = fOut
  s.HEADER.NUMPART_THISFILE._DATA     = [n_elements(id),0,0,0,0,0]
  s.HEADER.NUMPART_TOTAL._DATA        = [n_elements(id),0,0,0,0,0]
  s.HEADER.MASSTABLE._DATA            = [0.0,0.0,0.0,0.0,0.0,0.0]
  s.HEADER.TIME._DATA                 = 0.0
  s.HEADER.REDSHIFT._DATA             = 0.0
  s.HEADER.BOXSIZE._DATA              = boxSize
  s.HEADER.NUMFILESPERSNAPSHOT._DATA  = 1
  s.HEADER.OMEGA0._DATA               = 0.0
  s.HEADER.OMEGALAMBDA._DATA          = 0.0
  s.HEADER.HUBBLEPARAM._DATA          = 1.0
  s.HEADER.FLAG_SFR._DATA             = 0
  s.HEADER.FLAG_COOLING._DATA         = 0  
  s.HEADER.FLAG_STELLARAGE._DATA      = 0
  s.HEADER.FLAG_METALS._DATA          = 0
  s.HEADER.FLAG_FEEDBACK._DATA        = 0
  s.HEADER.FLAG_DOUBLEPRECISION._DATA = 0
  
  s.HEADER.COMPOSITION_VECTOR_LENGTH._DATA = 0 ;?

  ; these are IO blocks for snapshots but not for ICs
  s1 = mod_struct(s.PARTTYPE0,'DENSITY',/delete)
  s1 = mod_struct(s1,'SMOOTHINGLENGTH',/delete)
  s1 = mod_struct(s1,'VOLUME',/delete)
  s = mod_struct(s,'PARTTYPE0',s1)

  ;s.PARTTYPE0.DENSITY._DATA[*]         = 0.0
  ;s.PARTTYPE0.SMOOTHINGLENGTH._DATA[*] = 0.0
  ;s.PARTTYPE0.VOLUME._DATA[*]          = 0.0
  
  ; modify data parameters
  s.PARTTYPE0._FILE                = fOut
  s.PARTTYPE0.COORDINATES._FILE    = fOut
  s.PARTTYPE0.VELOCITIES._FILE     = fOut
  s.PARTTYPE0.PARTICLEIDS._FILE    = fOut
  s.PARTTYPE0.INTERNALENERGY._FILE = fOut
    
  ; modify data
  s.PARTTYPE0.COORDINATES._DIMENSIONS    = [3,n_elements(id)]
  s.PARTTYPE0.COORDINATES._NELEMENTS     = n_elements(pos)
  s1 = mod_struct(s.PARTTYPE0.COORDINATES,'_DATA',pos) ;change _DATA size
  s2 = mod_struct(s.PARTTYPE0,'COORDINATES',s1) ;update PARTTYPE0 with child

  s.PARTTYPE0.VELOCITIES._DIMENSIONS     = [3,n_elements(id)]
  s.PARTTYPE0.VELOCITIES._NELEMENTS      = n_elements(vel)
  s1 = mod_struct(s.PARTTYPE0.VELOCITIES,'_DATA',vel)
  s2 = mod_struct(s2,'VELOCITIES',s1)
  
  s.PARTTYPE0.PARTICLEIDS._DIMENSIONS    = [n_elements(id)]
  s.PARTTYPE0.PARTICLEIDS._NELEMENTS     = n_elements(id)
  s1 = mod_struct(s.PARTTYPE0.PARTICLEIDS,'_DATA',id)
  s2 = mod_struct(s2,'PARTICLEIDS',s1)

  s.PARTTYPE0.INTERNALENERGY._DIMENSIONS = [n_elements(u)]
  s.PARTTYPE0.INTERNALENERGY._NELEMENTS  = n_elements(u)
  s1 = mod_struct(s.PARTTYPE0.INTERNALENERGY,'_DATA',u)
  s2 = mod_struct(s2,'INTERNALENERGY',s1)
  
  s.PARTTYPE0.MASSES._DIMENSIONS = [n_elements(massOrDens)]
  s.PARTTYPE0.MASSES._NELEMENTS  = n_elements(massOrDens)
  s1 = mod_struct(s.PARTTYPE0.MASSES,'_DATA',massOrDens)
  s2 = mod_struct(s2,'MASSES',s1)
  
  s = mod_struct(s,'PARTTYPE0',s2) ;import new PARTTYPE0 structure

  ; output
  h5_create, fOut, s

end

; addICBackgroundGrid(): add background grid of specified resolution nBackGrid^3 of size boxSize 
;                        centered at [0,0,0] (or [boxCen,boxCen,boxCen] if specified) to gas ICs 
;                        (only add background cells that would be empty)

function addICBackgroundGrid, gas, boxSize=boxSize, boxCen=boxCen, nBackGrid=nBackGrid, $
                              massBackGrid=massBackGrid, uthermBackGrid=uthermBackGrid

  ; if not requested, return un-altered
  if (nBackGrid eq 0) then return, gas
  
  if (n_elements(boxSize) eq 0 or n_elements(gas) eq 0) then stop

  ; config
  if ~keyword_set(massBackGrid) then massBackGrid   = 1e-20
  if ~keyword_set(uthermBackGrid) then uthermBackGrid = 0.0

  backCellSize = boxSize / nBackGrid
  
  x_back = findgen(nBackGrid)/nBackGrid * boxSize + backCellSize/2.0
  y_back = x_back
  z_back = x_back
  
  if keyword_set(boxCen) then begin
    x_back += boxCen[0] - boxSize/2.0
    y_back += boxCen[1] - boxSize/2.0
    z_back += boxCen[2] - boxSize/2.0
  endif

  nBackKeep = 0
  pos_back = []
  
  ; find empty background grid cells
  for i=0,nBackGrid-1 do begin
    for j=0,nBackGrid-1 do begin
      for k=0,nBackGrid-1 do begin
        cenBackCell = [x_back[i],y_back[j],z_back[k]]
        min_xyz = cenBackCell - backCellSize/2.0
        max_xyz = cenBackCell + backCellSize/2.0
        
        w = where(gas.pos[0,*] ge min_xyz[0] and gas.pos[0,*] le max_xyz[0] and $
                  gas.pos[1,*] ge min_xyz[1] and gas.pos[1,*] le max_xyz[1] and $
                  gas.pos[2,*] ge min_xyz[2] and gas.pos[2,*] le max_xyz[2], count)
        
        ; keep background point
        if (count eq 0) then begin
          nBackKeep += 1
          pos_back = [[pos_back],[cenBackCell]]
        endif
      endfor
    endfor
  endfor
  
  print,'Added ['+str(nBackKeep)+'] background cells of '+str(nBackGrid)+$
        '^3 ('+str(nBackGrid^3)+') inside boxSize = '+string(boxSize)+' suggest meanVolume = '+$
        string(backCellSize^3.0)
  
  ; create other arrays
  vel_back  = fltarr(3,nBackKeep)
  mass_back = fltarr(nBackKeep) + massBackGrid
  u_back    = fltarr(nBackKeep) + uthermBackGrid
  id_back   = lindgen(nBackKeep) + max(gas.id) + 1
  
  ; concat background grid and primary cells
  pos  = [[gas.pos],[pos_back]]
  vel  = [[gas.vel],[vel_back]]
  mass = [gas.mass,mass_back]
  u    = [gas.u,u_back]
  id   = [gas.id,id_back]

  r = {pos:pos,vel:vel,mass:mass,u:u,id:id}
  return, r
  
end
