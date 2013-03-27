; sphere.pro
; gas accretion project - interpolation of quantities onto healpix spheres
; dnelson mar.2013

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
                         Nside=Nside, radFacs=radFacs, cutSubS=cutSubS
                         
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  nNGB = 32  ; neighbor search in CalcTHVal
  
  ; healpix resolution parameter, 8=768, 16~3k, 32~12k, 64~50k, 128~200k, 256~750k, 512~3M
  if ~keyword_set(Nside) then Nside = 64
  
  ; r/r_vir list of shells to compute
  if ~keyword_set(radFacs) then $
    radFacs = [0.01,0.05,0.1,0.25,0.5,0.75,0.9,1.0,1.1,1.25,1.5,1.75,2.0]
    
  padFac = 2.0*max(radFacs) ; times r_vir maximum search radius  
    
  if ~keyword_set(valName) then message,'Error: Must specify valName'
    
  ; if one subgroupID, check for existence of a save
  if keyword_set(cutSubS) then csTag = '.cutSubS' else csTag = ''
  
  saveFilename = sP.derivPath+'hShells/hShells.'+partType+'.'+valName+'.'+sP.savPrefix+str(sP.res)+csTag+$
                 '.ns'+str(Nside)+'.'+str(sP.snap)+'.h'+str(subgroupIDs[0])+'.'+str(n_elements(radFacs)) + '.sav'
 
  if n_elements(subgroupIDs) eq 1 then begin
    if file_test(saveFilename) then begin
      restore,saveFilename
      return,r
    endif
  endif
  
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
                                
    if n_elements(subgroupIDs) gt 1 then message,'Computed and saved, to return radmassflux request only one halo.'                       
    hsv_dens.valName = 'radmassflux'
    hsv_dens.value = alog10(hsv_dens.value * units.UnitMass_in_Msun) * (hsv_radvel.value * units.kmS_in_kpcYr)
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
  if sP.trMCPerCell ne 0 and valName ne 'density' then begin
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
    
    temp_pres_ent = convertUtoTemp(u,nelec,/log)

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
      mass = replicate(h.massTable[partTypeNum('dm')],h.nPartTot[partTypeNum('dm')]) $
    else $
      mass = loadSnapshotSubset(sP=sP,partType=partType,field='mass')
  endif
    
  if valName eq 'angmom' then begin
    pos = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
    vel = loadSnapshotSubset(sP=sP,partType=partType,field='vel')
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
    saveFilename = sP.derivPath+'hShells/hShells.'+partType+'.'+valName+'.'+sP.savPrefix+str(sP.res)+csTag+$
                   '.ns'+str(Nside)+'.'+str(sP.snap)+'.h'+str(subgroupID)+'.'+str(n_elements(radFacs)) + '.sav'

    if file_test(saveFilename) then begin
      print,'Skipping: ',subgroupID
      continue
    endif

    ; find halo position and virial radius
    cenPos = sgpos[*,subgroupID]
    rVir   = gc.group_r_crit200[gc.subgroupGrNr[subgroupID]]
    
    ; verify that the requested subgroupID is a primary subgroup
    ;w = where(priSGIDs eq subgroupID,countMatch)
    ;if ~countMatch then message,'Error: Only know how to do this for primary subgroups for now.'
    
    ; take conservative subset of points using periodic distances
    rad = periodicDists(cenPos,pos,sP=sP)
    
    wRadCut = where(rad le padFac*rVir,sCount)
    if ~sCount then message,'Error: No positions found near specified radius.'
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
      
      ; make velocities relative to bulk halo motion
      gVel = gc.subgroupVel[*,subgroupID]
      loc_vel[0,*] = reform(loc_vel[0,*] - gVel[0])
      loc_vel[1,*] = reform(loc_vel[1,*] - gVel[1])
      loc_vel[2,*] = reform(loc_vel[2,*] - gVel[2])
      
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
      
      ; dot(vel,rnorm) gives the magnitude of the projection of vel onto vec(r)
      value = (loc_vel[0,*]*rnorm0 + loc_vel[1,*]*rnorm1 + loc_vel[2,*]*rnorm2) / rnorm; 1xN
      
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
      
      loc_vel = vel[*,wRadCut]
      
      ; make velocities relative to bulk halo motion
      gVel = gc.subgroupVel[*,subgroupID]
      loc_vel[0,*] = reform(loc_vel[0,*] - gVel[0])
      loc_vel[1,*] = reform(loc_vel[1,*] - gVel[1])
      loc_vel[2,*] = reform(loc_vel[2,*] - gVel[2])
      
      ; angular momentum magnitude
      jvec = fltarr(3,sCount)
      jvec[0,*] = rvec1 * loc_vel[2,*] - rvec2 * loc_vel[1,*]
      jvec[1,*] = rvec2 * loc_vel[0,*] - rvec0 * loc_vel[2,*]
      jvec[2,*] = rvec0 * loc_vel[1,*] - rvec1 * loc_vel[0,*]
      jnorm = sqrt(jvec[0,*]*jvec[0,*] + jvec[1,*]*jvec[1,*] + jvec[2,*]*jvec[2,*])
      
      posval = [loc_pos,jnorm]
      thMode = 1 ; mean
    endif
    
    if valName eq 'veldisp' then begin
      value = reform(veldisp[wRadCut],[1,n_elements(wRadCut)]) ; 1xN
      
      posval = [loc_pos,value]
      thMode = 1 ; mean
    endif
    
    if valName eq 'density' then begin
      value = reform(mass[wRadCut],[1,n_elements(wRadCut)]) ; 1xN
      
      posval = [loc_pos,value]
      thMode = 3 ; total/volume
    endif
    
    if n_elements(posval) eq 0 then message,'Error: Unrecognized value name.'
    
    ; if cutting substructure, match secondary ids against local ids
    if keyword_set(cutSubS) then begin
      loc_ids = ids[wRadCut]
      
      ; make a list of satellites of this halo and their particle ids
      ;nSubs    = gc.groupNSubs[gc.subgroupGrNr[subgroupID]]
      ;firstSub = gc.groupFirstSub[gc.subgroupGrNr[subgroupID]]
      ;satGCids = lindgen(nSubs-1) + firstSub + 1
      ;satPIDs = gcPIDList(gc=gc,select='sec',valGCids=satGCids,partType=partType)
      
      ; remove the intersection of (satPIDs,loc_ids) from posval
      calcMatch,satPIDs,loc_ids,sat_ind,ids_ind,count=count
      sat_ind = !NULL
      
      all = bytarr(n_elements(loc_ids))
      if count gt 0 then all[ids_ind] = 1B
      wSubSComp = where(all eq 0B, ncomp)
      
      print,'Substructures cut ['+str(count)+'] of ['+str(n_elements(loc_ids))+'] have left: '+str(ncomp)
      
      ids_ind = !NULL
      loc_ids = !NULL

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
          value      : fltarr(Npx,nRadFacs)   }
  
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
      
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
    posval = !NULL
    value  = !NULL
    sphereXYZ = !NULL
  
  endforeach
  
  if n_elements(subgroupIDs) eq 1 then return,r
end

; calcShellInfallFrac(): calculate angular covering fraction of infalling mass flux

function calcShellInfallFrac, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  radFacs     = [0.25,0.5,1.0] ; fractions of the virial radius to save results
  cutSubS     = 1     ; cut satellite substructures out from halo
  minLogMass  = 10.0  ; minimum halo mass
  partType    = 'gas'
  valName     = 'radmassflux'
  threshVal   = 0.0 ; strict inflow
  
  ; check for existence of save
  saveFilename = sP.derivPath + 'shellFrac.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.mm' + str(fix(minLogMass*10)) + '.cut' + str(cutSubS) + '.rad' + $
    str(n_elements(radFacs)) + '.' + partType + '.' + valName + '.sav'
    
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif
  
  ; load subgroup IDs
  gc        = loadGroupCat(sP=sP,/skipIDs)
  priIDs    = gcIDList(gc=gc,select='pri')
  priMasses = codeMassToLogMsun(gc.subgroupMass[priIDs])
  
  w = where(priMasses ge minLogMass,count)
  subgroupIDs = priIDs[w]
  priMasses = priMasses[w]
  
  ; arrays
  r = { fracInfall  : fltarr(n_elements(radFacs),n_elements(subgroupIDs))  ,$
        subgroupIDs : subgroupIDs         ,$
        priMasses   : priMasses           ,$
        count       : count               ,$
        radFacs     : radFacs             ,$
        nRadFacs    : n_elements(radFacs) ,$
        cutSubS     : cutSubS             ,$
        minLogMass  : minLogMass          ,$
        threshVal   : threshVal            }

  ; interpolate all halos (and save)
  ;hsv = haloShellValue(sP=sP,partType=partType,valName=valName,subgroupIDs=subgroupIDs,$
  ;                     cutSubS=cutSubS,radFacs=radFacs)

  print,'calculating...'
  foreach subgroupID,subgroupIDs,k do begin
    if k mod 100 eq 0 then print,k
    
    ; interpolate onto the shell (load)
    hsv = haloShellValue(sP=sP,partType=partType,valName=valName,subgroupIDs=[subgroupID],$
                         cutSubS=cutSubS,radFacs=radFacs)
    
    ; data
    for radInd=0,n_elements(radFacs)-1 do begin
      healpix_data = reform(hsv.value[*,radInd])
      
      ; calculate sky covering fraction with log(rho/mean rho) > threshold
      w = where(healpix_data ge threshVal, countAbove, ncomp=countBelow)
      r.fracInfall[radInd,k] = float(countAbove) / hsv.nPx
      
      ;print,string(subgroupID,format='(i5)')+"  r="+string(radFacs[radInd],format='(f4.2)')+"  "+$
      ;      string(r.fracInfall[radInd,k]*100,format='(f5.2)')+'%'
    endfor
       
  endforeach
  
  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  return, r
end
