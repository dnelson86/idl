; sphere.pro
; gas accretion project - interpolation of quantities onto healpix spheres
; dnelson dec.2014

; sphereXYZCoords(): return a HEALPix subdivision at Nside resolution parameter scaled out to radius

function sphereXYZCoords, Nside=Nside, radius=radius, center=center

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; init HEALPix
  init_healpix
  
  w = where(!HEALPIX.Nside eq Nside,count)
  if (count ne 1) then message,'ERROR: Bad Nside.'
  
  ; set parameters
  Npix = nside2npix(Nside)
  
  ; get (x,y,z) positions
  pxIDs = lindgen(Npix)
  
  pix2vec_nest, Nside, pxIDs, vec
  
  ; rescale to radius
  if keyword_set(radius) then vec *= radius
  
  ; tolerance for values near zero
  epsTol = 1e-15
  
  w = where(abs(vec) lt epsTol,count)
  if (count gt 0) then vec[w] = 0.0
  
  ; transpose to shape [3,N]
  vec = transpose(vec)
  
  ; re-position to center
  if keyword_set(center) then begin
    vec[0,*] += center[0]
    vec[1,*] += center[1]
    vec[2,*] += center[2]
  endif
  
  return, vec
end

; haloShellValue(): for a given snapshot and subgroupID, evaluate some property/function value (valName)
;                   of a specified particle type on a series of radial shells (radFacs)
;
; cutSubS = cut substructures (satellite subgroups) out before estimating densities

function haloShellValue, sP=sP, partType=partType, valName=valName, subgroupIDs=subgroupIDs, $
                         Nside=Nside, radFacs=radFacs, cutSubS=cutSubS, newSaves=newSaves
                         
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits(redshift=sP.redshift)
  
  if partType ne 'dm' and partType ne 'gas' then message,'Error: Check. Thinkt through for stars.'
  
  ; config
  nNGB = 20  ; neighbor search in CalcTHVal
  if n_elements(cutSubS) eq 0 then cutSubS = 0
  
  ; healpix resolution parameter, 8=768, 16~3k, 32~12k, 64~50k, 128~200k, 256~750k, 512~3M
  if ~keyword_set(Nside) then Nside = 64
  
  ; r/r_vir list of shells to compute
  if ~keyword_set(radFacs) then $
    radFacs = [0.01,0.05,0.1,0.25,0.5,0.75,0.9,1.0,1.1,1.25,1.5,1.75,2.0]
    
  padFac = 2.0*max(radFacs) ; times r_vir maximum search radius (cubic cutout0)
    
  if ~keyword_set(valName) then message,'Error: Must specify valName'
    
  ; if one subgroupID, check for existence of a save
  csTag = ''
  if cutSubS eq 1 then csTag = '.cutSubS'
  if cutSubS eq 2 then csTag = '.cutSubS2'
  
  saveFilename = sP.derivPath+'hShells/snap_'+str(sP.snap)+'/hShells.'+partType+'.'+valName+'.'+$
                 sP.savPrefix+str(sP.res)+csTag+'.ns'+str(Nside)+'.'+str(sP.snap)+'.h'+$
                 str(subgroupIDs[0])+'.'+str(n_elements(radFacs)) + '.sav'
 
  if n_elements(subgroupIDs) eq 1 then begin
    if file_test(saveFilename) and ~keyword_set(newSaves) then begin
      restore,saveFilename
      return,r
    endif
  endif
  
  ; if multiple subgroupIDs, check for existence of all saves
  workFlag = 0
  
  saveFilenames = sP.derivPath+'hShells/snap_'+str(sP.snap)+'/hShells.'+partType+'.'+valName+'.'+$
                  sP.savPrefix+str(sP.res)+csTag+'.ns'+str(Nside)+'.'+str(sP.snap)+'.h'+$
                  str(subgroupIDs)+'.'+str(n_elements(radFacs)) + '.sav'
                  
  foreach fname, saveFilenames do begin
    if ~file_test(fname) then workFlag = 1
  endforeach
  if workFlag eq 0 and ~keyword_set(newSaves) then begin
    print,valName+': All already done ('+str(n_elements(subgroupIDs))+' subgroups), returning.'
    return,[]
  endif
  
  file_mkdir, sP.derivPath + 'hShells/snap_' + str(sP.snap) ; if it does not already exist
  
  Npx      = nside2npix(Nside)
  nRadFacs = n_elements(radFacs)
  
  ; "radial mass flux" as density * radial velocity (area normalization / r^2 sphere factor omitted)
  ; instead of estimating this value on each particle and then tophat smoothing, we use the
  ; estimates of each avaiable already on each point on the sphere
  if valName eq 'radmassflux' then begin
    hsv_dens   = haloShellValue(sP=sP,partType=partType,valName='density',subgroupIDs=subgroupIDs,$
                                Nside=Nside,radFacs=radFacs,cutSubS=cutSubS)
    hsv_radvel = haloShellValue(sP=sP,partType=partType,valName='radvel',subgroupIDs=subgroupIDs,$
                                Nside=Nside,radFacs=radFacs,cutSubS=cutSubS)
                                
    if n_elements(subgroupIDs) gt 1 then begin
      print,'radmassflux: Computed and saved, to return request only one halo.'
      return,[]
    endif
  
    if ~tag_exist(hsv_dens,'valName') then return,hsv_dens ; flag=-1
  
    hsv_dens.valName = 'radmassflux'
    hsv_dens.value = hsv_dens.value * units.UnitMass_in_Msun * (hsv_radvel.value * units.kmS_in_kpcYr)
    hsv_dens.value *= 1e6 ; Msun / kpc^2 / year  -->  Msun / kpc^2 / Myr
    return, hsv_dens 
  endif
  
  ; load group catalog and primary list
  h  = loadSnapshotHeader(sP=sP)
  if keyword_set(cutSubS) then gc = loadGroupCat(sP=sP,/verbose,/readIDs) $
  else gc = loadGroupCat(sP=sP,/verbose,/skipIDs)
  
  priSGIDs = gcIDList(gc=gc,select='pri')
  sgpos  = subgroupPosByMostBoundID(sP=sP)
  
  ; if this is an Arepo run, we need the cell masses to use as weights
  if sP.trMCPerCell ne 0 then begin
    if partType eq 'dm' then $
      mass = replicate(h.massTable[partTypeNum('dm')],h.nPartTot[partTypeNum('dm')]) $
    else $
      mass = loadSnapshotSubset(sP=sP,partType=partType,field='mass')
  endif
  
  ; load particle positions
  pos = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
  
  ; load data value requested
  if valName eq 'temp' then begin
    u = loadSnapshotSubset(sP=sP,partType=partType,field='u')
    nelec = loadSnapshotSubset(sP=sP,partType=partType,field='nelec')
    
    temp_pres_ent = convertUtoTemp(u,nelec)

    u = !NULL
    nelec = !NULL
  endif
  
  if valName eq 'pressure' then begin
    u = loadSnapshotSubset(sP=sP,partType=partType,field='u')
    dens = loadSnapshotSubset(sP=sP,partType=partType,field='dens')
    
    temp_pres_ent = calcPressureCGS(u,dens,sP=sP)
    
    u = !NULL
    dens = !NULL
  endif
  
  if valName eq 'entropy' then begin
    u = loadSnapshotSubset(sP=sP,partType=partType,field='u')
    dens = loadSnapshotSubset(sP=sP,partType=partType,field='dens')
    
    temp_pres_ent = calcEntropyCGS(u,dens,sP=sP)

    u = !NULL
    dens = !NULL
  endif
  
  if valName eq 'metallicity' then $
    metal_sfr = loadSnapshotSubset(sP=sP,partType=partType,field='metallicity')
  
  if valName eq 'sfr' then $
    metal_sfr = loadSnapshotSubset(sP=sP,partType=partType,field='sfr')
  
  if valName eq 'radvel' then $
    vel = loadSnapshotSubset(sP=sP,partType=partType,field='vel')
  
  if valName eq 'density' then begin
    if partType eq 'dm' then $
      mass = replicate(h.massTable[partTypeNum('dm')],h.nPartTot[partTypeNum('dm')])
      
    if partType eq 'gas' then $
      dens = loadSnapshotSubset(sP=sP,partType=partType,field='density')
  endif
    
  if valName eq 'angmom' then begin
    pos = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
    vel = loadSnapshotSubset(sP=sP,partType=partType,field='vel')
  endif
  
  if valName eq 'radvel' or valName eq 'angmom' then begin
    ; convert particle velocities to proper
    scalefac = 1.0 / (1.0 + sP.redshift)
    if scalefac le 0.0 or scalefac gt 1.0 then message,'Error'
    vel *= sqrt(scalefac)
  endif
  
  if valName eq 'veldisp' then begin
    veldisp = loadHsmlDir(sP=sP,partType=partType,/readVelDisp,/verbose)
  endif
  
  ; if cutting substructures, load ids and make secondary list
  if keyword_set(cutSubS) then begin
    ids = loadSnapshotSubset(sP=sP,partType=partType,field='ids')
    satPIDs = gcPIDList(gc=gc,select='sec',partType=partType)
  endif
  
  ; loop over each requested halo
  foreach subgroupID,subgroupIDs do begin
    ; skip if this halo exists
    saveFilename = sP.derivPath+'hShells/snap_'+str(sP.snap)+'/hShells.'+partType+'.'+valName+'.'+$
                   sP.savPrefix+str(sP.res)+csTag+'.ns'+str(Nside)+'.'+str(sP.snap)+'.h'+$
                   str(subgroupID)+'.'+str(n_elements(radFacs)) + '.sav'

    if file_test(saveFilename) and ~keyword_set(newSaves) then begin
      print,'Skipping: ',subgroupID
      continue
    endif

    ; find halo position and virial radius
    cenPos = sgpos[*,subgroupID]
    rVir   = gc.group_r_crit200[gc.subgroupGrNr[subgroupID]]
    
    ; verify that the requested subgroupID is a primary subgroup
    w = where(priSGIDs eq subgroupID,countMatch)
    if ~countMatch then message,'Error: Only know how to do this for primary subgroups for now.'
    
    ; take conservative subset of points using periodic distances
    rad = periodicDists(cenPos,pos,sP=sP)
    wRadCut = where(rad le padFac*rVir,sCount)

    if sCount eq 0 then begin
      r = {flag:-1,valName:valName,value:0}
      save,r,filename=saveFilename
      print,'Warning: No positions found near specified radius:'
      print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))+' (FLAG=-1)'
      continue
    endif
    
    rad = !NULL
    
    loc_pos = pos[*,wRadCut]
    
    ; create local data value and make posval = [[pos],[val]]
    if valName eq 'temp' or valName eq 'pressure' or valName eq 'entropy' then begin     
      value = reform(temp_pres_ent[wRadCut],[1,n_elements(wRadCut)]) ; make 1xN vector
      
      posval = [loc_pos,value]
      thMode = 1 ; mean
    endif
    
    if valName eq 'metallicity' or valName eq 'sfr' then begin
      value = reform(metal_sfr[wRadCut],[1,n_elements(wRadCut)]) ; make 1xN vector
      
      posval = [loc_pos,value]
      thMode = 1 ; mean
    endif
    
    if valName eq 'radvel' then begin
      loc_vel = vel[*,wRadCut]
      
      ; make normalized position vector wrt halo center = vec(r) / ||r|| where r from particle to center
      ; means that radvel<0 is inflow and radvel>0 is outflow
      rnorm0 = reform(loc_pos[0,*] - cenPos[0])
      rnorm1 = reform(loc_pos[1,*] - cenPos[1])
      rnorm2 = reform(loc_pos[2,*] - cenPos[2])
      
      correctPeriodicDistVecs, rnorm0, sP=sP
      correctPeriodicDistVecs, rnorm1, sP=sP
      correctPeriodicDistVecs, rnorm2, sP=sP
      
      ; denominator
      rnorm = sqrt(rnorm0*rnorm0 + rnorm1*rnorm1 + rnorm2*rnorm2)
      
      ; make velocities relative to bulk halo motion
      gVel = gc.subgroupVel[*,subgroupID] ; already peculiar (no scalefac correction needed)
      
      loc_vel[0,*] = reform(loc_vel[0,*] - gVel[0])
      loc_vel[1,*] = reform(loc_vel[1,*] - gVel[1])
      loc_vel[2,*] = reform(loc_vel[2,*] - gVel[2])

      ; dot(vel,rnorm) gives the magnitude of the projection of vel onto vec(r)
      value = (loc_vel[0,*]*rnorm0 + loc_vel[1,*]*rnorm1 + loc_vel[2,*]*rnorm2) / rnorm ; 1xN
      
      ; add hubble flow
      value += rnorm * units.H_z
      
      posval = [loc_pos,value]
      thMode = 1 ; mean
    endif
    
    if valName eq 'angmom' then begin
      ; mean magnitude of specific angular momentum = rvec x vel
      rvec0 = reform(cenPos[0] - loc_pos[0,*])
      rvec1 = reform(cenPos[1] - loc_pos[1,*])
      rvec2 = reform(cenPos[2] - loc_pos[2,*])
      
      correctPeriodicDistVecs, rvec0, sP=sP
      correctPeriodicDistVecs, rvec1, sP=sP
      correctPeriodicDistVecs, rvec2, sP=sP
      
      rvec0 /= units.HubbleParam ; remove little h from Coordinates for angmom
      rvec1 /= units.HubbleParam
      rvec2 /= units.HubbleParam
      
      rvec0 *= scalefac ; convert to peculiar for angmom
      rvec1 *= scalefac
      rvec2 *= scalefac
      
      loc_vel = vel[*,wRadCut]
      
      ; make velocities relative to bulk halo motion
      gVel = gc.subgroupVel[*,subgroupID] ; already peculiar (no scalefac correction needed)
      loc_vel[0,*] = reform(loc_vel[0,*] - gVel[0])
      loc_vel[1,*] = reform(loc_vel[1,*] - gVel[1])
      loc_vel[2,*] = reform(loc_vel[2,*] - gVel[2])
      
      ; angular momentum magnitude
      jvec = fltarr(3,sCount)
      jvec[0,*] = rvec1 * loc_vel[2,*] - rvec2 * loc_vel[1,*]
      jvec[1,*] = rvec2 * loc_vel[0,*] - rvec0 * loc_vel[2,*]
      jvec[2,*] = rvec0 * loc_vel[1,*] - rvec1 * loc_vel[0,*]
      jnorm = reform( sqrt(jvec[0,*]*jvec[0,*] + jvec[1,*]*jvec[1,*] + jvec[2,*]*jvec[2,*]) )
      
      posval = [loc_pos,jnorm]
      thMode = 1 ; mean
    endif
    
    if valName eq 'veldisp' then begin
      value = reform(veldisp[wRadCut],[1,n_elements(wRadCut)]) ; 1xN
      
      posval = [loc_pos,value]
      thMode = 1 ; mean
    endif
    
    if valName eq 'density' then begin
      ; gas: mean (take average of stored gas cell densities falling within kernel)
      if partType eq 'gas' then begin
        value = reform(dens[wRadCut],[1,n_elements(wRadCut)]) ; 1xN
        thMode = 1
      endif
      
      ; dm: total/volume (use the adaptive kernel as a spatial density estimator)
      if partType eq 'dm' then begin
        value = reform(mass[wRadCut],[1,n_elements(wRadCut)]) ; 1xN
        thMode = 3
      endif
      
      posval = [loc_pos,value]
    endif
    
    if n_elements(posval) eq 0 then message,'Error: Unrecognized value name.'
    
    ; if cutting substructure, match secondary ids against local ids
    if cutSubS eq 1 then begin
      loc_ids = ids[wRadCut]
      
      ; remove the intersection of (satPIDs,loc_ids) from posval
      calcMatch,satPIDs,loc_ids,sat_ind,ids_ind,count=count
      sat_ind = !NULL
      
      all = bytarr(n_elements(loc_ids))
      if count gt 0 then all[ids_ind] = 1B
      wSubSComp = where(all eq 0B, ncomp)
      
      ;print,'Substructures cut ['+str(count)+'] of ['+str(n_elements(loc_ids))+'] have left: '+str(ncomp)
      
      ids_ind = !NULL
      loc_ids = !NULL

      if ncomp gt 0 then posval = posval[*,wSubSComp]
    endif
    
    ; remove substructures (all particles in subfind groups other than target)? load ids and make secondary list
    if cutSubS eq 2 then begin
      ; verify our target ID is a primary (otherwise...)
      priGCIDs = gcIDList(gc=gc,select='pri')
      if total(subgroupID eq priGCIDs) eq 0 then message,'Error'
      
      loc_ids = ids[wRadCut]
      
      ; get ids of all gas cells in all substructures other than the target
      ; (this differs from ids in the target, as it includes cells not in any subfind group)
      myPIDs  = gcPIDList(gc=gc,valGCids=[subgroupID],partType='gas')
      subPIDs = gcPIDList(gc=gc,select='all',partType='gas')
      remPIDs = removeIntersectionFromB(myPIDs,subPIDs)
      if n_elements(remPIDs) ne n_elements(subPIDs)-n_elements(myPIDs) then message,'Error'
      
      ; remove the intersection of (satPIDs,ids) from wRadCut
      calcMatch,remPIDs,loc_ids,ind1,ind2,count=count
      
      all = bytarr(n_elements(loc_ids))
      if count gt 0 then all[ind2] = 1B ; flag substructure gas cells, then select against
      wSubSComp = where(all eq 0B, ncomp)
      
      if ncomp gt 0 then posval = posval[*,wSubSComp]
    endif
  
    ; allocate save structure
    r = { Nside      : Nside                 ,$
          Npx        : Npx                   ,$
          subgroupID : subgroupID            ,$
          sP         : sP                    ,$
          rVir       : rVir                  ,$
          cenPos     : cenPos                ,$
          partType   : partType              ,$
          radFacs    : radFacs               ,$
          padFac     : padFac                ,$
          sCount     : sCount                ,$
          nNGB       : nNGB                  ,$
          nRadFacs   : nRadFacs              ,$
          valName    : valName               ,$
          cutSubS    : cutSubS               ,$
          flag       : 1                     ,$ ; 1=ok, -1=bad
          value      : fltarr(Npx,nRadFacs)   }
  
    ; check if we have very few points
    npos = size(posval)
    if npos[2] le nNGB then begin
      r.flag = -1
      save,r,filename=saveFilename
      print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))+' (FLAG=-1)'
      continue
    endif
  
    sphereXYZ = fltarr(3,Npx*nRadFacs)
  
    ; loop over all requested shells and generate all the sphere points
    for i=0,nRadFacs-1 do begin
      radius = radFacs[i] * rVir ;kpc
    
      ; get sphere (x,y,z) positions
      locSphereXYZ = sphereXYZCoords(Nside=Nside,radius=radius,center=cenPos)
  
      ; periodic wrap any sphere points that landed outside the box (periodic ok in CalcHSMLds)
      w = where(locSphereXYZ lt 0.0,count)
      if count gt 0 then locSphereXYZ[w] += sP.boxSize
      w = where(locSphereXYZ gt sP.boxSize,count)
      if count gt 0 then locSphereXYZ[w] -= sP.boxSize
      
      ; store
      sphereXYZ[*,i*Npx:(i+1)*Npx-1] = locSphereXYZ
    endfor
  
    ; calculate tophat density estimate of all points on all spheres (one tree build)
    if sP.trMCPerCell eq 0 then begin
      r.value = CalcTHVal(posval,sphereXYZ,ndims=3,nNGB=nNGB,thMode=thMode,boxSize=sP.boxSize)
    endif else begin
      ; if this is an Arepo run, make the mass subset for weighting and do tophat estimate
      weights = reform(mass[wRadCut],[1,n_elements(wRadCut)])
      if keyword_set(cutSubS) then weights = weights[0,wSubSComp]
      posvalwt = [posval,weights]
      r.value = CalcTHVal(posvalwt,sphereXYZ,ndims=3,nNGB=nNGB,thMode=thMode,boxSize=sP.boxSize,/weighted)
    endelse
    
    ; units, log
    if valName eq 'temp' then r.value = alog10(r.value)
      
    if valName eq 'density' then begin
      ; density unit conversion (from code_mass/code_volume)
      r.value = codeDensToPhys(r.value, sP=sP, /cgs)
      r.value /= units.mass_proton ;g/cm^3 -> cm^(-3) number density
    endif
      
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
    posval = !NULL
    value  = !NULL
    sphereXYZ = !NULL
  
  endforeach
  
  if n_elements(subgroupIDs) eq 1 then return,r
end

; uniformHaloMassSelection(): select primary halo indices with uniform density in mass
;   e.g. counteract the steep mass function

function uniformHaloMassSelection, sP=sP, gc=gc, minMass=minMass, binSize=binSize, midMass=midMass, $
                                   priMasses=priMasses, verbose=verbose, massBinCens=massBinCens
                                  
  if ~keyword_set(sP) or ~keyword_set(minMass) or ~keyword_set(binSize) or ~keyword_set(midMass) then $
    message,'Error: Must specify inputs'
  if ~keyword_set(gc) then gc = loadGroupCat(sP=sP,/skipIDs)
  
  priIDs    = gcIDList(gc=gc,select='pri')
  parMasses = codeMassToLogMsun(gc.subgroupMass[priIDs])
  
  w = where(parMasses ge midMass,countAbove)
  subgroupIDs = priIDs[w]
  priMasses = parMasses[w]
  
  print,'Above ['+string(midMass,format='(f4.1)')+'] picked = '+str(countAbove)
  
  ; make selection by bin below midMass
  seed    = 42424242L
  nPerBin = 40 ; make constant across different 512 runs
  
  w = where(parMasses ge minMass and parMasses lt midMass,count)
  nBins = floor((midMass - minMass) / binSize)
  
  for i=0,nBins-1 do begin
    binMin = midMass - (i+1)*binSize
    binMax = binMin + binSize
    
    wBin = where(parMasses ge binMin and parMasses lt binMax,count)
    if count eq 0 then continue
    
    ; if specific massBin points specified, only include these
    if n_elements(massBinCens) gt 0 then begin
      if min( abs(binMin-massBinCens) ) gt 0.01 then continue
    endif
    
    ; highest mass bin, store number and accept all
    if i eq 0 then begin
      wSelect = wBin
      ;nPerBin = count
    endif else begin
      ; shuffle indices and pick the first nPerBin
      wBin = shuffle(wBin, seed=seed)
      wSelect = wBin[0:( (nPerBin-1) < (n_elements(wBin)-1) )]
    endelse
    
    ; add selection to keeper arrays
    subgroupIDs = [subgroupIDs, priIDs[wSelect]]
    priMasses   = [priMasses, parMasses[wSelect]]
    
    if keyword_set(verbose) then $
      print,' ['+string(i,format='(I2)')+'] binMinMax: '+string(binMin,format='(f4.1)')+' to '+$
        string(binMax,format='(f4.1)')+' [pick '+str(n_elements(wSelect))+' of '+str(count)+']'
  endfor
  
  if keyword_set(verbose) then $
    print,'Below ['+string(midMass,format='(f4.1)')+'] picked = '+str(n_elements(priMasses)-countAbove)
    
  return, subgroupIDs
end

; calcShellInfallFrac(): calculate angular covering fraction of infalling mass flux

function calcShellInfallFrac, sP=sP, massBinCens=massBinCens, blindSearch=blindSearch
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  radFacs      = [0.25,0.5,0.75,1.0,1.5] ; fractions of the virial radius to save results
  cutSubS      = 1     ; cut satellite substructures out from halo
  threshValAbs = 20.0  ; e.g. greater than +/- 20 km/s vrad
  threshValRel = 0.1   ; e.g. greater than +/- 10% of circular velocity of NFW halo fit at that rad
  
  ; config - halo selection
  allAboveLogMsun = 11.0 ; every halo above this value
  minLogMsun      = 9.0  ; minimum halo mass
  belowBinSize    = 0.1  ; below allAboveLogMsun, randomly get equal number per mass bin
  
  ; check for existence of save
  massStr = 'm' + string(allAboveLogMsun,format='(f4.1)') + '-' + string(minLogMsun,format='(f3.1)') + $
            '-' + string(belowBinSize,format='(f3.1)')
  if n_elements(massBinCens) gt 0 then massStr += '.mBC' + str(n_elements(massBinCens))
  
  saveFilename = sP.derivPath + 'hShells/shellFrac.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.cut' + str(cutSubS) + '.rad' + $
    str(n_elements(radFacs)) + '.' + massStr + '.infallFrac.sav'
    
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif else begin
    ; if we are searching for snapshots which have been analyzed, return a false here
    if keyword_set(blindSearch) then return,[]
  endelse
  
  ; make halo selection load subgroup IDs
  subgroupIDs = uniformHaloMassSelection(sP=sP, minMass=minLogMsun, binSize=belowBinSize, $
                                         midMass=allAboveLogMsun, priMasses=priMasses, $
                                         massBinCens=massBinCens, /verbose)
                                         
  ; interpolate all halos (and save)
  hsv = haloShellValue(sP=sP,partType='gas',valName='radmassflux',subgroupIDs=subgroupIDs,$
                       cutSubS=cutSubS,radFacs=radFacs)  
  
  ; save array
  r = { fracInfall_noTh   : fltarr(n_elements(radFacs),n_elements(subgroupIDs)) + !values.f_nan ,$
        fracInfall_relTh  : fltarr(n_elements(radFacs),n_elements(subgroupIDs)) + !values.f_nan ,$
        fracInfall_absTh  : fltarr(n_elements(radFacs),n_elements(subgroupIDs)) + !values.f_nan ,$
        fracOutflow_noTh  : fltarr(n_elements(radFacs),n_elements(subgroupIDs)) + !values.f_nan ,$ ; 1-infall by def
        fracOutflow_relTh : fltarr(n_elements(radFacs),n_elements(subgroupIDs)) + !values.f_nan ,$
        fracOutflow_absTh : fltarr(n_elements(radFacs),n_elements(subgroupIDs)) + !values.f_nan ,$
        massFluxIn        : fltarr(n_elements(radFacs),n_elements(subgroupIDs)) + !values.f_nan ,$
        massFluxOut       : fltarr(n_elements(radFacs),n_elements(subgroupIDs)) + !values.f_nan ,$
        subgroupIDs       : subgroupIDs         ,$
        priMasses         : priMasses           ,$
        radFacs           : radFacs             ,$
        nRadFacs          : n_elements(radFacs) ,$
        cutSubS           : cutSubS             ,$
        threshValRel      : threshValRel        ,$
        threshValAbs      : threshValAbs         }

  print,'calculating...'
  foreach subgroupID,subgroupIDs,k do begin
    if k mod 10 eq 0 then print,k
    
    ; interpolate onto the shell (load)
    hsv_vrad = haloShellValue(sP=sP,partType='gas',valName='radvel',subgroupIDs=[subgroupID],$
                              cutSubS=cutSubS,radFacs=radFacs)
    hsv = haloShellValue(sP=sP,partType='gas',valName='radmassflux',subgroupIDs=[subgroupID],$
                         cutSubS=cutSubS,radFacs=radFacs)
                         
    ; check for skipped halo
    if tag_exist(hsv,'flag') then begin
      if hsv.flag lt 0 then continue
    endif
    
    ; area per healpixel in code units (cKpc^2)
    pixel_areas = (4 * !pi * (radFacs*hsv.rVir)^2.0) / hsv.nPx
                         
    ; get nfw model for this halo mass
    nfw = nfw_profile(radFacs,mass=logMsunToCodeMass(priMasses[k]),redshift=sP.redshift)

    ; loop over radii
    for radInd=0,n_elements(radFacs)-1 do begin
      healpix_data = reform(hsv.value[*,radInd])
      healpix_data_vrad = reform(hsv_vrad.value[*,radInd])
      
      ; fracInfall/Outflow
      w = where(healpix_data lt 0.0, count, ncomp=ncomp)
      
      r.fracInfall_noTh[radInd,k]  = float(count) / hsv.nPx
      r.fracOutflow_noTh[radInd,k] = float(ncomp) / hsv.nPx
      
      w = where(healpix_data_vrad lt -threshValAbs, count, ncomp=ncomp)
      
      r.fracInfall_absTh[radInd,k]  = float(count) / hsv.nPx
      r.fracOutflow_absTh[radInd,k] = float(ncomp) / hsv.nPx
      
      w = where(healpix_data_vrad lt -threshValRel*nfw.vCirc_DM[radInd], count, ncomp=ncomp)
      
      r.fracInfall_relTh[radInd,k]  = float(count) / hsv.nPx
      r.fracOutflow_relTh[radInd,k] = float(ncomp) / hsv.nPx
      
      ; massFluxIn/Out (multiply by pixel area, Msun/h / kpc^2 / Myr -> Msun/h/Myr)
      w = where(healpix_data lt 0.0,count,comp=wc,ncomp=ncomp)
      
      if count gt 0 then r.massFluxIn[radInd,k]  = total( healpix_data[w] * pixel_areas[radInd] )
      if ncomp gt 0 then r.massFluxOut[radInd,k] = total( healpix_data[wc] * pixel_areas[radInd] )
    endfor
       
  endforeach
  
  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  return, r
end

; calcShellValuesAcrossRedshifts()

function calcShellValuesAcrossRedshifts, run=run, doSnap=doSnap, getParams=getParams
  compile_opt idl2, hidden, strictarr, strictarrsubs
  if n_elements(run) eq 0 then message,'Error'
  
  ; config
  res         = 512
  runs        = ['feedback','tracer','gadget']
  snapFac     = 2 ; 2=every other snap, 3=every third snap, etc
  maxRedshift = 6.0 ; don't go back further than this
  minRedshift = 0.0 ; don't go forward further than this
  massBinCens = [9.0,9.5,10.0,10.5,11.0,11.5,12.0]
  
  ; make authoritative list of target snapshots
  sP = simParams(res=512,run='feedback',redshift=0.0) ; master
  ss = snapNumToRedshift(sP=sP,/all)
  
  targetSnaps = where(ss ge minRedshift and ss lt maxRedshift)
  targetSnaps = targetSnaps[0:-1:snapFac]
  
  snapRedshifts = ss[targetSnaps]
  
  ; for this one particular run, match the authoritative list to a list for this run
  sP = simParams(res=res,run=run,snap=snap)
  targetSnaps = redshiftToSnapNum(snapRedshifts,sP=sP)
   
  ; returning parameters? exit now
  if keyword_set(getParams) then return, {targetSnaps:targetSnaps, massBinCens:massBinCens}
  
  if n_elements(doSnap) gt 0 then begin
    ; process just one snapshot?
    print,'Running snapshot ['+str(doSnap)+'] only!'
       
    sP = simParams(res=res,run=run,snap=doSnap)
    x = calcShellInfallFrac(sP=sP,massBinCens=massBinCens)

  endif else begin
    ; process all selected snapshots (will skip already processed snaps with calcShellInfallFrac save)
    print,'Running ['+str(n_elements(targetSnaps))+'] snapshots:',targetSnaps

    foreach snap,targetSnaps do begin
      print,' -- [' + str(snap) + '] --'
      sP = simParams(res=res,run=run,snap=snap)
      x = calcShellInfallFrac(sP=sP,massBinCens=massBinCens)
    endforeach
  endelse
  
end
