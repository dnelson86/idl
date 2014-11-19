; illustrisVisCutout.pro
; illustris 1820^3 specialized visualization (cutouts/preprocessing)
; dnelson nov.2013

; makeHDF5SnapFromCutout(): make Arepo/ArepoVTK format snapshot out of cutout files for ray-tracing/vis

pro wrapHDF5Write, name, data, loc_id, attribute=attribute, dataset=dataset
  compile_opt idl2, hidden, strictarr, strictarrsubs
  if ~keyword_set(attribute) and ~keyword_set(dataset) then message,'Error'
  
  ; get data type and space, needed to create the dataset
  datatype_id  = H5T_IDL_CREATE(data)
  dataspace_id = H5S_CREATE_SIMPLE(size(data,/DIMENSIONS))
  
  if keyword_set(dataset) then begin
    ; create dataset in the output file
    dataset_id = H5D_CREATE(loc_id,name,datatype_id,dataspace_id)
    H5D_WRITE,dataset_id,data
    H5D_CLOSE,dataset_id
  endif
  
  if keyword_set(attribute) then begin
    ; create attribute in the output file
    attribute_id = H5A_CREATE(loc_id,name,datatype_id,dataspace_id)
    H5A_WRITE,attribute_id,data
    H5A_CLOSE,attribute_id
  endif
  
  ; close open identifiers
  H5S_CLOSE,dataspace_id
  H5T_CLOSE,datatype_id
  
end

pro makeHDF5SnapFromCutout
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  ;sP = simParams(res=1820,run='illustris',snap=135)
  sP = simParams(res=1,run='scylla',redshift=2.0)
  haloID  = 0
  sizeFac = 9.0
  dmAsGas = 0 ; if 1, write DM positions/mass as GAS PartType0 (Voronoi DM estimator)
  
  ; load
  gcInd = getMatchedIDs(haloID=haloID,simParams=sP)

  saveFileBase = sP.derivPath + 'cutouts/cutout.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.h' + str(gcInd) + '.sf' + str(fix(sizeFac*10))
    
  restore,saveFileBase + '_meta.sav'
    
  if dmAsGas gt 0 then begin
    print,'Writing DM as Gas!'
    ; replace nCutout with number of DM in cutout
    restore,saveFileBase + '_dm_pos.sav'
    size = size(quant)
    nCutout = ulong(size[2])
  endif
    
  NumPart = ulonarr(6)
  NumPart[0] = ulong(nCutout)
    
  ; open HDF5 and write header
  if dmAsGas eq 0 then dmTag = '' else dmTag = 'DM_'
  
  saveFileHdf5 = sP.derivPath + 'cutouts/cutout.' + sP.savPrefix + str(sP.res) + $
    '.h' + str(gcInd) + '.sf' + str(fix(sizeFac*10)) + '_' + dmTag + str(sP.snap) + '.hdf5'
    
  fid = H5F_CREATE(saveFileHdf5)
  hid = H5G_CREATE(fid,'Header')

    wrapHDF5Write, 'NumPart_ThisFile',          NumPart,                     hid, /attribute
    wrapHDF5Write, 'NumPart_Total',             NumPart,                     hid, /attribute
    wrapHDF5Write, 'NumPart_Total_HighWord',    ulonarr(6),                  hid, /attribute
    wrapHDF5Write, 'MassTable',                 dblarr(6),                   hid, /attribute
    wrapHDF5Write, 'Time',                      dblarr(1),                   hid, /attribute
    wrapHDF5Write, 'Redshift',                  dblarr(1),                   hid, /attribute
    wrapHDF5Write, 'BoxSize',                   [double(2*boxSize)],         hid, /attribute
    wrapHDF5Write, 'NumFilesPerSnapshot',       [long(1)],                   hid, /attribute
    wrapHDF5Write, 'Omega0',                    [double(units.Omega_M)],     hid, /attribute
    wrapHDF5Write, 'OmegaLambda',               [double(units.Omega_L)],     hid, /attribute
    wrapHDF5Write, 'HubbleParam',               [double(units.hubbleParam)], hid, /attribute
    wrapHDF5Write, 'Flag_Sfr',                  lonarr(1),                   hid, /attribute
    wrapHDF5Write, 'Flag_Cooling',              lonarr(1),                   hid, /attribute
    wrapHDF5Write, 'Flag_StellarAge',           lonarr(1),                   hid, /attribute
    wrapHDF5Write, 'Flag_Metals',               lonarr(1),                   hid, /attribute
    wrapHDF5Write, 'Flag_Feedback',             lonarr(1),                   hid, /attribute
    wrapHDF5Write, 'Flag_DoublePrecision',      lonarr(1),                   hid, /attribute
    wrapHDF5Write, 'Composition_vector_length', lonarr(1),                   hid, /attribute
    
  H5G_close,hid

  ; make group for PartType0 (skip: Density, SmoothingLength, Volume, ...)
  gid = H5G_CREATE(fid,'PartType0')
  
    if dmAsGas eq 0 then begin
      ; write gas
      restore,saveFileBase + '_gas_pos.sav'
      
      ; centered halo is at (0,0), shift to boxcenter (boxsize) of [0,2*boxsize]
      quant += boxSize
      w = where(quant lt 0.0 or quant gt 2*boxSize,count)
      if count gt 0 then message,'Check potential points outside box'
    
      if sP.res eq 1820 and sP.snap eq 135 and sizeFac eq 2.0 then begin
        print,'WARNING: ADDING NOISE AROUND DEGEN LOCATION.'
        dPos = [406.8,3167.6,21.3]
        dists = periodicDists(dPos, quant, sP=sP)
        w = where(dists le 10.0,count)
        print,'Count: ',count
        dists = !NULL
        seed = 42424244L
        quant[0,w] += 0.01 * randomu(seed,count,/uniform) ; seed turns into 6x6 state matrix
        quant[1,w] += 0.01 * randomu(seed,count,/uniform) ; pass it back to continue sequence
        quant[2,w] += 0.01 * randomu(seed,count,/uniform)
      endif
      
      if sP.run eq 'scylla' and sP.redshift eq 2.0 then begin
        print,'WARNING: ADDING SCYLLA NOISE AROUND DEGEN LOCATION.'
        posList = list([2447.3,3306.8,1789.5],[1711.6,1708.07,1710.5],[1192.3,3192.7,3250.1],$
                       [1709.9,1711.1,1710.4],[2451.4,3266.4,1855.7])
        foreach dPos,posList,k do begin
          dists = periodicDists(dPos, quant, sP=sP)
          w = where(dists le 2.5,count)
          print,k,' Count: ',count
          dists = !NULL
          seed = 42424244L
          quant[0,w] += 0.01 * randomu(seed,count,/uniform) ; seed turns into 6x6 state matrix
          quant[1,w] += 0.01 * randomu(seed,count,/uniform) ; pass it back to continue sequence
          quant[2,w] += 0.01 * randomu(seed,count,/uniform)
        endforeach
      endif
    
      wrapHDF5Write, 'Coordinates', quant, gid, /dataset
    
      restore,saveFileBase + '_gas_temp.sav'
      wrapHDF5Write, 'InternalEnergy', quant, gid, /dataset
      
      restore,saveFileBase + '_gas_density.sav'
      wrapHDF5Write, 'Density', quant, gid, /dataset
    
      restore,saveFileBase + '_gas_mass.sav'
      wrapHDF5Write, 'Masses', quant, gid, /dataset
    
      restore,saveFileBase + '_gas_sfr.sav'
      wrapHDF5Write, 'StarFormationRate', quant, gid, /dataset
    
      if file_test(saveFileBase + '_gas_metal.sav') then begin
        restore,saveFileBase + '_gas_metal.sav'
        wrapHDF5Write, 'Metallicity', quant, gid, /dataset
      endif
    
      restore,saveFileBase + '_gas_nelec.sav'
      wrapHDF5Write, 'ElectronAbundance', quant, gid, /dataset
    
      restore,saveFileBase + '_gas_vel.sav'
      wrapHDF5Write, 'Velocities', quant, gid, /dataset
    endif else begin
      ; write DM AS GAS!
      h = loadSnapshotHeader(sP=sP)
      
      ;restore,saveFileBase + '_dm_pos.sav' ; already loaded
      quant += boxSize
      w = where(quant lt 0.0 or quant gt 2*boxSize,count)
      if count gt 0 then message,'Check potential points outside box'
    
      wrapHDF5Write, 'Coordinates', quant, gid, /dataset
      
      ; write constant masses out
      quant = fltarr(nCutout) + h.massTable[ partTypeNum('dm') ]
      wrapHDF5Write, 'Masses', quant, gid, /dataset
      
      ; write vmag,vdisp as vel[0,1] components
      vel = fltarr(3,nCutout)
      restore,saveFileBase + '_dm_vmag.sav'
      vel[0,*] = quant
      restore,saveFileBase + '_dm_vdisp.sav'
      vel[1,*] = quant
      
      wrapHDF5Write, 'Velocities', vel, gid, /dataset
    endelse
    
    quant = lindgen(nCutout) + 1L
    wrapHDF5Write, 'ParticleIDs', quant, gid, /dataset
    
  ; close
  H5G_CLOSE,gid
  H5F_CLOSE,fid
  
  print,'Saved: '+saveFileHdf5
  
end

; illustrisVisCutout(): do spatial cubic cutout from the box, one quantity at a time for 1820^3

function illustrisVisCutout, sP=sP, gcInd=gcInd, sizeFac=sizeFac, quantName=quantName
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  if ~keyword_set(sP) or n_elements(gcInd) ne 1 or n_elements(sizeFac) eq 0 then message,'Error'
  
  saveFileBase = sP.derivPath + 'cutouts/cutout.' + sP.savPrefix + str(sP.res) + '.' + $
      str(sP.snap) + '.h' + str(gcInd) + '.sf' + str(fix(sizeFac*10)) + '_'
  
  ; if one quantity requested, load and return
  if n_elements(quantName) gt 0 then begin
    saveFilename = saveFileBase + quantName + '.sav'
      
    ; dosen't exist? make the cutout now
    if ~file_test(saveFilename) then dummy = illustrisVisCutout(sP=sP, gcInd=gcInd, sizeFac=sizeFac)
    
    restore,saveFilename
    return,quant
  endif
  
  ; proceed with cutout
  gc = loadGroupCat(sP=sP,/skipIDs,/skipOffsets,/verbose)
  
  ; halo properties
  haloGrNr = gc.subgroupGrNr[gcInd]
  haloRVir = gc.group_r_crit200[ haloGrNr ] ;ckpc
  haloMass = codeMassToLogMsun(gc.subgroupMass[gcInd])
  haloM200 = codeMassToLogMsun(gc.group_m_crit200[ haloGrNr ])
  haloV200 = sqrt(units.G * gc.subgroupMass[gcInd] / haloRVir )
  haloPos  = gc.groupPos[*, haloGrNr]
  
  boxSize  = ceil(sizeFac * haloRVir / 10.0) * 10.0
  
  print,'Halo gcInd',gcInd
  print,'Halo GrNr ',haloGrNr
  print,'Halo rVir ',haloRvir
  print,'Halo Mass ',haloMass
  print,'Halo pos  ',haloPos
  print,'boxSize   ',boxSize

  ; if 'inds' does not exist, make the spatial subselection now
  saveFilename = saveFileBase + 'gas_inds.sav'
  if ~file_test(saveFilename) then begin
    print,'Determining indices...'
    ; load all gas positions
    h   = loadSnapshotHeader(sP=sP)
    pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
    
    ; account for periodic box
    for i=0,2 do begin
      xyzDist = pos[i,*] - haloPos[i]
      correctPeriodicDistVecs, xyzDist, sP=sP
      pos[i,*] = xyzDist
    endfor
    
    xyzDist = !NULL
    
    ; local (cube) cutout
    inds = where(abs(pos[0,*]) le boxSize and $
                 abs(pos[1,*]) le boxSize and $
                 abs(pos[2,*]) le boxSize,nCutout)
                 
    fracCutout = float(nCutout) / h.nPartTot[ partTypeNum('gas') ] * 100
    print,'Found ['+str(nCutout)+'] ('+string(fracCutout,format='(f4.1)')+'%), shuffling...'
    
    ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
    iseed = 424242L
    sort_inds = sort(randomu(iseed,n_elements(inds)))
    inds = inds[sort_inds]
    sort_inds = !NULL
    
    ; save the indices
    save,inds,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
    
    ; save the metadata
    saveFilename = saveFileBase + 'meta.sav'
    save,haloRVir,haloMass,haloM200,haloV200,haloPos,boxSize,nCutout,fracCutout,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
    
    ; since we also have them, save the relative positions
    quant = pos[*,inds]
    pos   = !NULL
    
    saveFilename = saveFileBase + 'gas_pos.sav'
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif else begin
    restore,saveFilename
  endelse
  
  if n_elements(inds) eq 0 then message,'Error'
  
  ; vrad
  saveFilename = saveFileBase + 'gas_vrad.sav'
  if ~file_test(saveFilename) then begin
    if n_elements(quant)/3 ne n_elements(inds) then $ ; lazy load
      quant = illustrisVisCutout(sP=sP, gcInd=gcInd, sizeFac=sizeFac, quantName='pos')

    ; load velocities
    vel = loadSnapshotSubset(sP=sP,partType='gas',field='vel',inds=inds)
    
    ; calculate norm of radial velocity vector
    rad = reform(quant[0,*]^2.0 + quant[1,*]^2.0 + quant[2,*]^2.0) ; quant=pos

    sgcen_vel = gc.subgroupVel[*, gcInd]
    
    quant = reform( (vel[0,*]-sgcen_vel[0]) * quant[0,*] + $
                    (vel[1,*]-sgcen_vel[1]) * quant[1,*] + $
                    (vel[2,*]-sgcen_vel[2]) * quant[2,*]) / sqrt(rad)
                       
    rad = !NULL
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; vel (xyz)
  saveFilename = saveFileBase + 'gas_vel.sav'
  if ~file_test(saveFilename) then begin
    if n_elements(vel)/3 ne n_elements(inds) then $ ; lazy load
      quant = loadSnapshotSubset(sP=sP,partType='gas',field='vel',inds=inds) $
    else $
      quant = temporary(vel)
      
    vel = !NULL
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; temp
  saveFilename = saveFileBase + 'gas_temp.sav'
  if ~file_test(saveFilename) then begin
    u     = loadSnapshotSubset(sP=sP,partType='gas',field='u',inds=inds)
    nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec',inds=inds)
    
    quant = convertUtoTemp(u,nelec)

    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; nelec
  saveFilename = saveFileBase + 'gas_nelec.sav'
  if ~file_test(saveFilename) then begin
    if n_elements(nelec) ne n_elements(inds) then $ ; lazy load
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec',inds=inds)

    quant = nelec
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; sfr
  saveFilename = saveFileBase + 'gas_sfr.sav'
  if ~file_test(saveFilename) then begin
    quant = loadSnapshotSubset(sP=sP,partType='gas',field='sfr',inds=inds)
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; metallicity
  saveFilename = saveFileBase + 'gas_metal.sav'
  if ~file_test(saveFilename) and snapshotFieldExists(sP=sP, field='gfm_metallicity') then begin
    quant = loadSnapshotSubset(sP=sP,partType='gas',field='metallicity',inds=inds)
        
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; xray
  saveFilename = saveFileBase + 'gas_xray.sav'
  if ~file_test(saveFilename) then begin
    if n_elements(u) ne n_elements(inds) then $ ; lazy load
      u   = loadSnapshotSubset(sP=sP,partType='gas',field='u',inds=inds)
    if n_elements(nelec) ne n_elements(inds) then $ ; lazy load
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec',inds=inds)
    
    quant = sqrt( convertUtoTemp(u,nelec) ) ; T^(1/2)
    
    ; zero contribution from gas with T<10^6 K
    w = where(quant le sqrt(1e6),count)
    if count gt 0 then begin
      quant[w] = 0.0
      print,'Zero ['+str(count)+'] of ['+str(n_elements(quant))+'] below 1e6 K.'
    endif
    
    ; load SFR, zero contribution from gas with SFR>0
    sfr = loadSnapshotSubset(sP=sP,partType='gas',field='sfr',inds=inds)
    w = where(sfr gt 0.0,count)
    if count gt 0 then begin
      quant[w] = 0.0
      print,'Zero ['+str(count)+'] of ['+str(n_elements(quant))+'] with nonzero SFR.'
    endif
    sfr = !NULL
    w = !NULL
    
    quant *= loadSnapshotSubset(sP=sP,partType='gas',field='dens',inds=inds) ; rho^1 * T^(1/2)
    
    if snapshotFieldExists(sP=sP, field='gfm_metallicity') then begin
      metalmu = loadSnapshotSubset(sP=sP,partType='gas',field='metallicity',inds=inds)
      metalmu = meanmolwt(Y=0.25,Z=metalmu) ; mu
    endif else begin
      metalmu = fltarr( n_elements(quant) ) + meanmolwt(Y=0.25,Z=0.1)
    endelse
    
    quant /= (metalmu*metalmu) ; mu^(-2) * rho * T^(1/2)
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; SZ y-parameter
  saveFilename = saveFileBase + 'gas_sz_y.sav'
  if ~file_test(saveFilename) then begin    
    if n_elements(nelec) ne n_elements(inds) then $ ; lazy load
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec',inds=inds)
    if n_elements(u) ne n_elements(inds) then $ ; lazy load
      u   = loadSnapshotSubset(sP=sP,partType='gas',field='u',inds=inds)
      
    if n_elements(metalmu) ne n_elements(inds) then begin ; lazy load
      if snapshotFieldExists(sP=sP, field='gfm_metallicity') then begin
        metalmu = loadSnapshotSubset(sP=sP,partType='gas',field='metallicity',inds=inds)
        metalmu = meanmolwt(Y=0.25,Z=metalmu) ; mu
      endif else begin
        metalmu = fltarr( n_elements(nelec) ) + meanmolwt(Y=0.25,Z=0.1)
      endelse
    endif
      
    quant = u * nelec * metalmu
      
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  nelec = !NULL
  metalmu = !NULL
  
  ; entropy
  saveFilename = saveFileBase + 'gas_entropy.sav'
  if ~file_test(saveFilename) then begin
    if n_elements(dens) ne n_elements(inds) then $ ; lazy load
      dens = loadSnapshotSubset(sP=sP,partType='gas',field='dens',inds=inds)
    if n_elements(u) ne n_elements(inds) then $ ; lazy load
      u   = loadSnapshotSubset(sP=sP,partType='gas',field='u',inds=inds)
      
    quant = calcEntropyCGS(u,dens,sP=sP)
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  u = !NULL
  
  ; dens
  saveFilename = saveFileBase + 'gas_density.sav'
  if ~file_test(saveFilename) then begin
    if n_elements(dens) ne n_elements(inds) then $ ; lazy load
      dens = loadSnapshotSubset(sP=sP,partType='gas',field='dens',inds=inds)
    
    quant = dens

    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
    
  dens = !NULL
  
  ; volumes (hsml)
  saveFilename = saveFileBase + 'gas_hsml.sav'
  if ~file_test(saveFilename) then begin
    ; load volumes (convert to sizes)
    quant = loadSnapshotSubset(sP=sP,partType='gas',field='vol',inds=inds)
    quant = (temporary(quant) * 3.0 / (4*!pi))^(1.0/3.0) ;cellrad [ckpc]
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; mass
  saveFilename = saveFileBase + 'gas_mass.sav'
  if ~file_test(saveFilename) then begin
    quant = loadSnapshotSubset(sP=sP,partType='gas',field='mass',inds=inds)
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  if ~snapshotFieldExists(sP=sP,field='SubfindHsml') then begin
    print,'WARNING: Did not find SubfindHsml in snapshot, skipping DM!'
    return,0
  endif
  
  ; make DM selection and save positions
  saveFilename = saveFileBase + 'dm_inds.sav'
  if ~file_test(saveFilename) then begin
    print,'Determining DM indices...'
    ; load all dm positions
    h   = loadSnapshotHeader(sP=sP)
    pos = loadSnapshotSubset(sP=sP,partType='dm',field='pos')
        
    ; account for periodic box
    for i=0,2 do begin
      xyzDist = pos[i,*] - haloPos[i]
      correctPeriodicDistVecs, xyzDist, sP=sP
      pos[i,*] = xyzDist
    endfor
    
    xyzDist = !NULL
    
    ; local (cube) cutout
    inds = where(abs(pos[0,*]) le boxSize and $
                 abs(pos[1,*]) le boxSize and $
                 abs(pos[2,*]) le boxSize,nCutout)
                 
    fracCutout = float(nCutout) / h.nPartTot[ partTypeNum('dm') ] * 100
    print,'Found ['+str(nCutout)+'] DM ('+string(fracCutout,format='(f4.1)')+'%), shuffling...'
    
    ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
    iseed = 424242L
    sort_inds = sort(randomu(iseed,n_elements(inds)))
    inds = inds[sort_inds]
    sort_inds = !NULL
    
    ; save the indices
    save,inds,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
    
    ; since we also have them, save the relative positions
    quant = pos[*,inds]
    pos   = !NULL
    
    saveFilename = saveFileBase + 'dm_pos.sav'
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
                
  endif else begin
    restore,saveFilename
  endelse
  
  ; dm_hsml
  saveFilename = saveFileBase + 'dm_hsml.sav'
  if ~file_test(saveFilename) then begin
    quant = loadSnapshotSubset(sP=sP,partType='dm',field='subfind_hsml',inds=inds)
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; dm_vmag
  saveFilename = saveFileBase + 'dm_vmag.sav'
  if ~file_test(saveFilename) then begin
    quant = loadSnapshotSubset(sP=sP,partType='dm',field='vel',inds=inds)
    
    ; calculate magnitude of radial velocity vector
    quant = reform( sqrt(quant[0,*]^2.0 + quant[1,*]^2.0 + quant[2,*]^2.0) )
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; dm_vdisp
  saveFilename = saveFileBase + 'dm_vdisp.sav'
  if ~file_test(saveFilename) then begin
    quant = loadSnapshotSubset(sP=sP,partType='dm',field='veldisp',inds=inds)
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; dm_annih
  saveFilename = saveFileBase + 'dm_annih.sav'
  if ~file_test(saveFilename) then begin
    quant = loadSnapshotSubset(sP=sP,partType='dm',field='subfind_dens',inds=inds)
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  print,'Done.'

end

; illustrisSliceCutout(): axes[3 int] = [image_x, image_y, slice_axis]

function illustrisSliceCutout, sP=sP, gcInd=gcInd, sliceWidth=sliceWidth, axes=axes, quantName=quantName
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  if ~keyword_set(sP) or n_elements(gcInd) ne 1 or n_elements(sliceWidth) ne 1 or $
    n_elements(axes) ne 3 then message,'Error'
  
  saveFileBase = sP.derivPath + 'cutouts/sliceOff1.' + sP.savPrefix + str(sP.res) + '.' + $
      str(sP.snap) + '.h' + str(gcInd) + '.sw' + str(fix(sliceWidth)) + '.axes' + $
      str(axes[0]) + str(axes[1]) + str(axes[2]) + '_'
  
  ; if one quantity requested, load and return
  if n_elements(quantName) gt 0 then begin
    saveFilename = saveFileBase + quantName + '.sav'
      
    ; dosen't exist? make the cutout now
    if ~file_test(saveFilename) then $
      dummy = illustrisSliceCutout(sP=sP, gcInd=gcInd, sliceWidth=sliceWidth, axes=axes)
    
    restore,saveFilename
    return,quant
  endif
  
  ; proceed with cutout
  h = loadSnapshotHeader(sP=sP)
  
  ; if 'inds' does not exist, make the spatial subselection now
  saveFilename = saveFileBase + 'gas_inds.sav'
  if ~file_test(saveFilename) then begin
    ; load group catalog for center
    print,'WARNING: Using z=0 fof position to center slice!'
    redshiftOld = sP.redshift
    snapOld     = sP.snap
    sP.redshift = 0.0
    sP.snap = redshiftToSnapNum(sP.redshift, sP=sP)
    
    gc = loadGroupCat(sP=sP,/skipIDs,/skipOffsets,/verbose)
      
    ; restore sP
    sP.snap     = snapOld
    sP.redshift = redshiftOld
    
    ; continue
    haloGrNr  = gc.subgroupGrNr[gcInd]
    xyzCenter = gc.groupPos[*, haloGrNr]
    gc = !NULL
  
    print,'Halo GrNr ',haloGrNr
    print,'Halo pos  ',xyzCenter
    
    ; offset slice
    xyzCenter[axes[2]] += sliceWidth*2.0
    print,'OFFSET NEW pos ',xyzCenter
    
    print,'Determining gas indices...'
    ; load all gas positions
    pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
        
    ; slice (zmin-zmax) restriction
    inds = where( abs(pos[axes[2],*] - xyzCenter[axes[2]]) le sliceWidth,nCutout)
                 
    fracCutout = float(nCutout) / h.nPartTot[ partTypeNum('gas') ] * 100
    print,'Found ['+str(nCutout)+'] ('+string(fracCutout,format='(f4.1)')+'%), shuffling...'
    
    ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
    iseed = 424242L
    sort_inds = sort(randomu(iseed,n_elements(inds)))
    inds = inds[sort_inds]
    sort_inds = !NULL
    
    ; save the indices
    save,inds,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
    
    ; save the metadata
    saveFilename = saveFileBase + 'meta.sav'
    save,xyzCenter,nCutout,fracCutout,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif else begin
    restore,saveFilename
    restore,saveFileBase + 'meta.sav'
  endelse
  
  if n_elements(inds) eq 0 then message,'Error'
  
  ; gas_pos
  saveFilename = saveFileBase + 'gas_pos.sav'
  if ~file_test(saveFilename) then begin
    if n_elements(pos)/3 ne h.nPartTot[partTypeNum('gas')] then $ ; lazy load
      pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  
    ; save the positions (orthogonal to slice axis only)
    quant = fltarr( 2, n_elements(inds) )
    for i=0,1 do quant[i,*] = reform( pos[ axes[i], inds ] )
    pos = !NULL
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; gas_vmag
  saveFilename = saveFileBase + 'gas_vmag.sav'
  if ~file_test(saveFilename) then begin
    quant = loadSnapshotSubset(sP=sP,partType='gas',field='vel',inds=inds)

    ; calculate magnitude of radial velocity vector
    quant = reform( sqrt(quant[0,*]^2.0 + quant[1,*]^2.0 + quant[2,*]^2.0) )
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; gas_temp
  saveFilename = saveFileBase + 'gas_temp.sav'
  if ~file_test(saveFilename) then begin
    u     = loadSnapshotSubset(sP=sP,partType='gas',field='u',inds=inds)
    nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec',inds=inds)
    
    quant = convertUtoTemp(u,nelec)
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; gas_xray
  saveFilename = saveFileBase + 'gas_xray.sav'
  if ~file_test(saveFilename) then begin
    if n_elements(u) ne n_elements(inds) then $ ; lazy load
      u   = loadSnapshotSubset(sP=sP,partType='gas',field='u',inds=inds)
    if n_elements(nelec) ne n_elements(inds) then $ ; lazy load
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec',inds=inds)
    
    quant = sqrt( convertUtoTemp(u,nelec) ) ; T^(1/2)
    
    ; zero contribution from gas with T<10^6 K
    w = where(quant le sqrt(1e6),count)
    if count gt 0 then begin
      quant[w] = 0.0
      print,'Zero ['+str(count)+'] of ['+str(n_elements(quant))+'] below 1e6 K.'
    endif
    
    ; load SFR, zero contribution from gas with SFR>0
    sfr = loadSnapshotSubset(sP=sP,partType='gas',field='sfr',inds=inds)
    w = where(sfr gt 0.0,count)
    if count gt 0 then begin
      quant[w] = 0.0
      print,'Zero ['+str(count)+'] of ['+str(n_elements(quant))+'] with nonzero SFR.'
    endif
    sfr = !NULL
    w = !NULL
    
    quant *= loadSnapshotSubset(sP=sP,partType='gas',field='dens',inds=inds) ; rho^1 * T^(1/2)
    
    metalmu = loadSnapshotSubset(sP=sP,partType='gas',field='metallicity',inds=inds)
    metalmu = meanmolwt(Y=0.25,Z=metalmu) ; mu
    quant /= (metalmu*metalmu) ; rho/mu^2 * T^(1/2)
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; gas SZ y-parameter
  saveFilename = saveFileBase + 'gas_sz_y.sav'
  if ~file_test(saveFilename) then begin
    if n_elements(metalmu) ne n_elements(inds) then begin ; lazy load
      metalmu = loadSnapshotSubset(sP=sP,partType='gas',field='metallicity',inds=inds)
      metalmu = meanmolwt(Y=0.25,Z=metalmu) ; mu
    endif
    if n_elements(nelec) ne n_elements(inds) then $ ; lazy load
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec',inds=inds)
    if n_elements(u) ne n_elements(inds) then $ ; lazy load
      u   = loadSnapshotSubset(sP=sP,partType='gas',field='u',inds=inds)
      
    quant = u * nelec * metalmu
      
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  nelec = !NULL
  metalmu = !NULL
  
  ; gas_entropy
  saveFilename = saveFileBase + 'gas_entropy.sav'
  if ~file_test(saveFilename) then begin
    dens  = loadSnapshotSubset(sP=sP,partType='gas',field='dens',inds=inds)
    if n_elements(u) ne n_elements(inds) then $ ; lazy load
      u   = loadSnapshotSubset(sP=sP,partType='gas',field='u',inds=inds)
      
    quant = calcEntropyCGS(u,dens,sP=sP)
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  u     = !NULL
  dens  = !NULL
  
  ; gas_metallicity
  saveFilename = saveFileBase + 'gas_metal.sav'
  if ~file_test(saveFilename) then begin
    quant = loadSnapshotSubset(sP=sP,partType='gas',field='metallicity',inds=inds)
    
    ; convert to log(metallicity) for positive values, otherwise set GFM_MIN_METAL = -20 (log)
    ;w = where(quant gt 0.0,count,comp=wc,ncomp=ncomp)
    ;if count gt 0 then quant[w] = alog10(quant[w])
    ;if ncomp gt 0 then quant[wc] = -20.0
    ;w  = !NULL
    ;wc = !NULL
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; gas_volumes (hsml)
  saveFilename = saveFileBase + 'gas_hsml_vol.sav'
  if ~file_test(saveFilename) then begin
    ; load volumes (convert to sizes)
    quant = loadSnapshotSubset(sP=sP,partType='gas',field='vol',inds=inds)
    quant = (temporary(quant) * 3.0 / (4*!pi))^(1.0/3.0) ;cellrad [ckpc]
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; gas_hsml (subfind N-neighbor based)
  saveFilename = saveFileBase + 'gas_hsml.sav'
  if ~file_test(saveFilename) then begin
    quant = loadSnapshotSubset(sP=sP,partType='gas',field='subfind_hsml',inds=inds)
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; gas_mass
  saveFilename = saveFileBase + 'gas_mass.sav'
  if ~file_test(saveFilename) then begin
    quant = loadSnapshotSubset(sP=sP,partType='gas',field='mass',inds=inds)
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  quant = !NULL
  inds  = !NULL
  
  ; make DM selection and save positions
  saveFilename = saveFileBase + 'dm_inds.sav'
  if ~file_test(saveFilename) then begin
    print,'Determining DM indices...'
    ; load all dm positions
    h   = loadSnapshotHeader(sP=sP)
    pos = loadSnapshotSubset(sP=sP,partType='dm',field='pos')
        
    ; slice (zmin-zmax) restriction
    inds = where( abs(pos[axes[2],*] - xyzCenter[axes[2]]) le sliceWidth,nCutout)
                 
    fracCutout = float(nCutout) / h.nPartTot[ partTypeNum('dm') ] * 100
    print,'Found ['+str(nCutout)+'] ('+string(fracCutout,format='(f4.1)')+'%), shuffling...'
    
    ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
    iseed = 424242L
    sort_inds = sort(randomu(iseed,n_elements(inds)))
    inds = inds[sort_inds]
    sort_inds = !NULL
    
    ; save the indices
    save,inds,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif else begin
    restore,saveFilename
  endelse
  
  ; dm_pos
  saveFilename = saveFileBase + 'dm_pos.sav'
  if ~file_test(saveFilename) then begin
    if n_elements(pos)/3 ne h.nPartTot[partTypeNum('dm')] then $ ; lazy load
      pos = loadSnapshotSubset(sP=sP,partType='dm',field='pos')
  
    ; save the positions (orthogonal to slice axis only)
    quant = fltarr( 2, n_elements(inds) )
    for i=0,1 do quant[i,*] = reform( pos[ axes[i], inds ] )
    pos = !NULL
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; dm_hsml
  saveFilename = saveFileBase + 'dm_hsml.sav'
  if ~file_test(saveFilename) then begin
    quant = loadSnapshotSubset(sP=sP,partType='dm',field='subfind_hsml',inds=inds)
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; dm_vel
  saveFilename = saveFileBase + 'dm_vmag.sav'
  if ~file_test(saveFilename) then begin
    quant = loadSnapshotSubset(sP=sP,partType='dm',field='vel',inds=inds)
    
    ; calculate magnitude of radial velocity vector
    quant = reform( sqrt(quant[0,*]^2.0 + quant[1,*]^2.0 + quant[2,*]^2.0) )
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; dm_vdisp
  saveFilename = saveFileBase + 'dm_vdisp.sav'
  if ~file_test(saveFilename) then begin
    quant = loadSnapshotSubset(sP=sP,partType='dm',field='veldisp',inds=inds)
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; dm_annih
  saveFilename = saveFileBase + 'dm_annih.sav'
  if ~file_test(saveFilename) then begin
    quant = loadSnapshotSubset(sP=sP,partType='dm',field='subfind_dens',inds=inds)
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; make STARS selection and save positions
  saveFilename = saveFileBase + 'stars_inds.sav'
  if ~file_test(saveFilename) then begin
    print,'Determining STARS indices...'
    ; load all dm positions
    h   = loadSnapshotHeader(sP=sP)
    pos = loadSnapshotSubset(sP=sP,partType='stars',field='pos')
        
    ; slice (zmin-zmax) restriction
    inds = where( abs(pos[axes[2],*] - xyzCenter[axes[2]]) le sliceWidth,nCutout)
                 
    fracCutout = float(nCutout) / h.nPartTot[ partTypeNum('stars') ] * 100
    print,'Found ['+str(nCutout)+'] ('+string(fracCutout,format='(f4.1)')+'%), shuffling...'
    
    ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
    iseed = 424242L
    sort_inds = sort(randomu(iseed,n_elements(inds)))
    inds = inds[sort_inds]
    sort_inds = !NULL
    
    ; save the indices
    save,inds,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif else begin
    restore,saveFilename
  endelse
  
  ; stars_pos
  saveFilename = saveFileBase + 'stars_pos.sav'
  if ~file_test(saveFilename) then begin
    if n_elements(pos)/3 ne h.nPartTot[partTypeNum('stars')] then $ ; lazy load
      pos = loadSnapshotSubset(sP=sP,partType='stars',field='pos')
  
    ; save the positions (orthogonal to slice axis only)
    quant = fltarr( 2, n_elements(inds) )
    for i=0,1 do quant[i,*] = reform( pos[ axes[i], inds ] )
    pos = !NULL
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif
  
  ; stars_mass
  saveFilename = saveFileBase + 'stars_mass.sav'
  if ~file_test(saveFilename) then begin
    quant = loadSnapshotSubset(sP=sP,partType='stars',field='mass',inds=inds)
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif  
  
  ; stars_hsml
  saveFilename = saveFileBase + 'stars_hsml.sav'
  if ~file_test(saveFilename) then begin
    quant = loadSnapshotSubset(sP=sP,partType='stars',field='subfind_hsml',inds=inds)
    
    save,quant,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  endif  
  
  print,'Done.'

end
