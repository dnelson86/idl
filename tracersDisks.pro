; tracersDisks.pro
; dnelson
; jan 2012
;
; dev for tracer particles related to disk (2d/3d) tests

; fixBoundaryDisk(): load one of Diego's circumstellar 2d disk ICs and strip out the boundary cells
;                    which have duplicate, negative IDs, and add a background grid in the central 
;                    hole and outer empty regions

pro fixBoundaryDisk, addBackGrid=addBackGrid

  in_file = 'ics.dat'
  out_file = 'ics.dat.new'
  
  ; load
  h = loadSnapshotHeader(in_file,snapNum='none')
  print,'Read: '+str(in_file)+' ('+str(h.nPartTot[0])+' total type0 particles)'
  
  pos    = loadSnapshotSubset(in_file,snapNum='none',partType='gas',field='pos')
  u      = loadSnapshotSubset(in_file,snapNum='none',partType='gas',field='u')
  masses = loadSnapshotSubset(in_file,snapNum='none',partType='gas',field='masses')
  ids    = loadSnapshotSubset(in_file,snapNum='none',partType='gas',field='ids')
  vels   = loadSnapshotSubset(in_file,snapNum='none',partType='gas',field='vel')
  
  ; make selection
  w = where(ids ge 0,ncomp=ncomp)
  
  pos    = pos[*,w]
  u      = u[w]
  masses = masses[w]
  ids    = ids[w]
  vels   = vels[*,w]
  
  print,'Removed ['+str(ncomp)+'] boundary particles.'
  
  ; add background grid if requested
  if keyword_set(addBackGrid) then begin
    nBackGrid = 32 ;squared
    
    backIDStart = max(ids)+1
    
    ; arrays
    back_pos    = fltarr(3,nBackGrid*nBackGrid)
    back_u      = fltarr(nBackGrid*nBackGrid)
    back_masses = fltarr(nBackGrid*nBackGrid)
    back_ids    = lonarr(nBackGrid*nBackGrid)
    back_vels   = fltarr(3,nBackGrid*nBackGrid)
    
    ; fill in properties
    deltax = h.boxSize/nBackGrid
    deltay = h.boxSize/nBackGrid
    
    for i=0L,nBackGrid-1 do begin
      for j=0L,nBackGrid-1 do begin
        pid = i+j*nBackGrid
        
        back_pos[0,pid] = i*deltax+deltax/2.0  
        back_pos[1,pid] = j*deltay+deltay/2.0  
        back_pos[2,pid] = 0.0
        
        back_ids[pid]    = backIDStart + pid 
        ;back_vels[*,pid] = 0.0
      endfor
    endfor    
    
    ; radii of all particles
    boxCen = h.boxSize/2.0
    
    rad      = reform(sqrt((pos[0,*]-boxCen)*(pos[0,*]-boxCen) + $
                           (pos[1,*]-boxCen)*(pos[1,*]-boxCen)))
    back_rad = reform(sqrt((back_pos[0,*]-boxCen)*(back_pos[0,*]-boxCen) + $
                           (back_pos[1,*]-boxCen)*(back_pos[1,*]-boxCen)))
    
    ; assign continuous properties for mass and u from the disk
    ; inner hole
    back_w = where(back_rad lt min(rad)*0.99)
    orig_w = where(rad lt min(rad)*1.01)
    
    back_u[back_w] += mean(u[orig_w])
    back_masses[back_w] += mean(masses[orig_w])
    
    ; outer border
    back_w = where(back_rad gt min(rad)*1.01)
    orig_w = where(rad gt max(rad)*0.99)

    back_u[back_w] += mean(u[orig_w])
    back_masses[back_w] += mean(masses[orig_w])
    
    ; select only those outside of disk region
    w = where(back_rad lt min(rad)*0.99 or back_rad gt max(rad)*1.01)
    
    back_pos    = back_pos[*,w]
    back_u      = back_u[w]
    back_masses = back_masses[w]
    back_ids    = back_ids[w]
    back_vels   = back_vels[*,w]
    
    ; concat disk and background
    pos    = [[pos],[back_pos]]
    u      = [u,back_u]
    masses = [masses,back_masses]
    ids    = [ids,back_ids]
    vels   = [[vels],[back_vels]]

    print,'Added ['+str(n_elements(back_ids))+'] background points from '+str(nBackGrid)+'^2 grid.'

  endif

  ; save new IC
  print,'Writing: ',out_file,' (kept '+str(n_elements(ids))+' particles)'

  writeICFileHDF5, out_file, h.boxSize, pos, vels, ids, masses, u

end