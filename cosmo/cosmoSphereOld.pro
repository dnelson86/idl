; cosmoSphereOld.pro
; gas accretion project - interpolation and visualization of quantities onto spheres (OLD UNUSED stuff)
; dnelson apr.2012

; det4x4(): take 4x4 matrix determinant where vi are row vectors (Laplace expansion)

function det4x4, v0, v1, v2, v3

  det = v0[3]*v1[2]*v2[1]*v3[0] - v0[2]*v1[3]*v2[1]*v3[0] - $
        v0[3]*v1[1]*v2[2]*v3[0] + v0[1]*v1[3]*v2[2]*v3[0] + $
                                                            $
        v0[2]*v1[1]*v2[3]*v3[0] - v0[1]*v1[2]*v2[3]*v3[0] - $
        v0[3]*v1[2]*v2[0]*v3[1] + v0[2]*v1[3]*v2[0]*v3[1] + $
                                                            $
        v0[3]*v1[0]*v2[2]*v3[1] - v0[0]*v1[3]*v2[2]*v3[1] - $
        v0[2]*v1[0]*v2[3]*v3[1] + v0[0]*v1[2]*v2[3]*v3[1] + $
                                                            $
        v0[3]*v1[1]*v2[0]*v3[2] - v0[1]*v1[3]*v2[0]*v3[2] - $
        v0[3]*v1[0]*v2[1]*v3[2] + v0[0]*v1[3]*v2[1]*v3[2] + $
                                                            $
        v0[1]*v1[0]*v2[3]*v3[2] - v0[0]*v1[1]*v2[3]*v3[2] - $
        v0[2]*v1[1]*v2[0]*v3[3] + v0[1]*v1[2]*v2[0]*v3[3] + $
                                                            $
        v0[2]*v1[0]*v2[1]*v3[3] - v0[0]*v1[2]*v2[1]*v3[3] - $
        v0[1]*v1[0]*v2[2]*v3[3] + v0[0]*v1[1]*v2[2]*v3[3]
        
  return, det
end

; barycenCoords(): is the pt[x,y,z] inside the tetra defined by the four vertices of (x,y,z)?
;                  if yes return the barycentric coordinates (b0+b1+b2+b3=1), else return -1

function barycenCoords, x, y, z, pt

 ;   ; T matrix
 ;   T = [ [x[0]-x[3], x[1]-x[3], x[2]-x[3]] , $
 ;         [y[0]-y[3], y[1]-y[3], y[2]-y[3]] , $
 ;         [z[0]-z[3], z[1]-z[3], z[2]-z[3]] ]
 ;   ; lambda 1,2,3,4
 ;   bc123 = invert(T) ## (pt-[x[3],y[3],z[3]])
 ;   bc123b = cramer(T,pt-[x[3],y[3],z[3]])

 ;   b4 = 1 - total(bc123)
 ;   b = [bc123[0],bc123[1],bc123[2],b4]
 ;
 ;   ; is point inside? return barycen coords
 ;   inside = (total(b gt 0 and b lt 1) eq 4)
 ;   
 ;   if (inside eq 1) then return, b
 ;   return, -1 ; point not inside tetra
    
  v0 = [x[0],y[0],z[1],1]
  v1 = [x[1],y[1],z[1],1]
  v2 = [x[2],y[2],z[2],1]
  v3 = [x[3],y[3],z[3],1]
  p0 = [pt[0],pt[1],pt[2],1]
  
  det0 = det4x4(v0,v1,v2,v3)
  det1 = det4x4(p0,v1,v2,v3)
  det2 = det4x4(v0,p0,v2,v3)
  det3 = det4x4(v0,v1,p0,v3)
  det4 = det4x4(v0,v1,v2,p0)
  
  ; if det0==0 then the tetra was degenerate in the first place
  if (det0 lt 0.0) then begin
    if ((det1 lt 0.0) and (det2 lt 0.0) and (det3 lt 0.0) and (det4 lt 0.0)) then $
      return, [det1/det0,det2/det0,det3/det0,det4/det0]
  endif else begin ;det0>0
    if ((det1 gt 0.0) and (det2 gt 0.0) and (det3 gt 0.0) and (det4 gt 0.0)) then $
      return, [det1/det0,det2/det0,det3/det0,det4/det0]
  endelse
  
  return,-1 ; point not inside tetra
end

function sphereInterp, pos, val, sphereXYZ, verbose=verbose

  ; debugging
  ;seed = 424242L
  ;pos = (randomu(seed,3,10000)-0.5) * 3.2
  ;rad = reform( sqrt(pos[0,*]^2.0 + pos[1,*]^2.0 + pos[2,*]^2.0) )
  ;print,minmax(rad)
  ;val = 100.0/rad^1.5
  ;print,minmax(val)

  ; alternate (modified shepard's method)
  start = systime(/sec)
  
  val_interp = grid3(pos[0,*],pos[1,*],pos[2,*],val,sphereXYZ[0,*],sphereXYZ[1,*],sphereXYZ[2,*])
  
  if keyword_set(verbose) then print,'ShepardMethod took: ['+str(systime(/sec)-start)+'] seconds.'

  return, reform(val_interp)

  ; delaunay triangulate input points using quickhull
  if 0 then begin
    start = systime(/sec)
    qhull, pos, dt, /delaunay, bounds=bounds
    if keyword_set(verbose) then print,'Triangulation took: ['+str(systime(/sec)-start)+'] seconds.'
    
    ; array for interpolated values
    val_interp = fltarr(n_elements(sphereXYZ[0,*]))
    
    start = systime(/sec)
    for i=0,n_elements(sphereXYZ[0,*])-1 do begin
      ; interpolate onto a 2x2x2 regular grid of (x,y,z) points near the target
      pt = sphereXYZ[*,i]
      res = qgrid3(pos, val, dt, delta=0.01,dimension=2,start=pt)
      
      ; save interpolated value
      val_interp[i] = mean(res)
    endfor
    if keyword_set(verbose) then print,'Interpolation took: ['+str(systime(/sec)-start)+'] seconds.'
  endif ;0
  
  ; code below for barycentric interpolation (not finished)
  if 0 then begin
    par_tt = lonarr(n_elements(sphereXYZ[0,*])) ;parent tetra indices
    
    ind_pos = calcNN(pos,sphereXYZ,boxSize=0,ndims=3) ;non-periodic
    
    ; debug: verify distances (IDL!)
    ind_pos2 = lonarr(n_elements(sphereXYZ[0,*]))
    for i=0,n_elements(sphereXYZ[0,*])-1 do begin
      pt = sphereXYZ[*,i]
      dists = sqrt( (pos[0,*]-pt[0])^2.0 + (pos[1,*]-pt[1])^2.0 + (pos[2,*]-pt[2])^2.0 )
      w = where(abs(dists) eq min(abs(dists)),count)
      if (count ne 1) then stop
      ind_pos2[i] = w[0]
    endfor
    if ( total(abs(ind_pos-ind_pos2)) ne 0 ) then print,'ERROR: DEBUG FAILED.'
    
    ; find all tetra containing this vertex and compute the barycentric coordinates of each
    for i=0,n_elements(sphereXYZ[0,*])-1 do begin
      pt = sphereXYZ[*,i]
      w = where(dt eq ind_pos[i],count)
      print,pt,count
      
      tt = (array_indices(dt,w))[1,*] ; tetra indices
      
      for j=0,n_elements(tt)-1 do begin
        ; get vertices of tetra and use double precision to help vs edge/vertex degen
        x = double(pos[0,dt[*,tt[j]]])
        y = double(pos[1,dt[*,tt[j]]])
        z = double(pos[2,dt[*,tt[j]]])
        
        b = barycenCoords(x,y,z,pt)
  
        if (n_elements(b) gt 1) then begin
          if (par_tt[i] ne 0) then begin print,'already have parent!' & stop & endif
          val_interp[i] = b[0] * val[dt[0,tt[j]]] + b[1] * val[dt[1,tt[j]]] + $
                          b[2] * val[dt[2,tt[j]]] + b[3] + val[dt[3,tt[j]]]
          par_tt[i] = tt[j]
          print,val_interp[i],par_tt[i]
        endif
      endfor
      
      ; check we found a parent
      if (par_tt[i] eq 0) then begin print,'no parent found!' & stop & endif
  
    endfor
  endif ;0
  
  ; debug
  ;start_PS,'test2.eps'
  ;  plotsym,0,/fill
  ;  fsc_plot,pos[0,*],pos[1,*],xtitle="x",ytitle="y",aspect=1.0,psym=8,symsize=0.5
  ;  for i=0,n_elements(dt[0,*])-1 do begin
  ;    for j=0,3 do begin
  ;      inds = [dt[j,i],dt[(j+1) mod 3,i],dt[(j+2) mod 3,i]]
  ;      fsc_plot,pos[0,inds],pos[1,inds],line=0,/overplot,thick=!p.thick-2.5
  ;    endfor
  ;  endfor
  ;end_PS
  
  stop
end

; interpHaloShell(): for a given snapshot and haloID, interpolation of the density field at a 
;                    specified fraction of the virial radius

function interpHaloShell, sP=sP, haloID=haloID, radFac=radFac, Nside=Nside
  
  ; check for existence of a save
  saveFilename = sP.derivPath + sP.savPrefix + str(sP.res) + '.hShell.nS=' + str(Nside) + '.snap=' + $
                 str(sP.snap) + '.' + str(haloID) + '.' + str(radFac*100) + '.sav'
                 
  if file_test(saveFilename) then begin
    restore,saveFilename,/verbose
  endif else begin
    ; load catalog and find halo position
    h  = loadSnapshotHeader(sP=sP)
    gc = loadGroupCat(sP=sP,/verbose,/skipIDs)
    
    haloPos = gc.groupPos[*,haloID]
    radius  = radFac * gc.group_r_crit200[haloID] ;kpc
    
    ; load gas positions and densities
    pos  = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
    dens = loadSnapshotSubset(sP=sP,partType='gas',field='dens')
      
    ; correction if shell is near or crosses periodic box edge (add ghosts)
    ; TODO
    if (total(haloPos-2*radius lt 0.0 or haloPos+2*radius gt h.boxSize) ne 0) then $
      message,'Warning: Halo near box edge. Periodic not implemented.'
    
    ; take conservative subset of points
    rad = reform( sqrt( (pos[0,*]-haloPos[0])*(pos[0,*]-haloPos[0]) + $
                        (pos[1,*]-haloPos[1])*(pos[1,*]-haloPos[1]) + $
                        (pos[2,*]-haloPos[2])*(pos[2,*]-haloPos[2]) ) )
    
    minMax = [0.1*radius,2.0*radius]
    if (minMax[0] lt 10.0) then minMax[0] = 0.0
    if (minMax[1] lt gc.group_r_crit200[haloID]) then minMax[1] = gc.group_r_crit200[haloID]
    
    w = where(rad ge minMax[0] and rad le minMax[1],count)
    if ~count then message,'Error: No positions found near specified radius.'

    pos  = pos[*,w]
    dens = alog10(rhoRatioToCrit(dens[w],redshift=sP.redshift))
    
    rad = !NULL
    
    ; get sphere (x,y,z) positions
    sphereXYZ = sphereXYZCoords(Nside=Nside,radius=radius,center=haloPos)
  
    ; rescale all positions to [0,0,0]-[1,1,1] by dividing by boxSize (optimize interp accuracy)
    pos /= h.boxSize
    sphereXYZ /= h.boxSize
  
    ; interpolate
    val_interp = sphereInterp(pos,dens,sphereXYZ,/verbose)

    ; save
    ;save,val_interp,Nside,radFac,haloID,count,filename=saveFilename
  endelse
  
  return, val_interp
end

; plotHaloShell(): test mollview plot of the sphereInterp results using the healpix points (OLD)

pro plotHaloShell
  message,'old'
  ; config
  Nside  = 64 ;256~750k, 512~3M
  radFac = 0.25 ; times r_vir
  haloID = 0
  
  sP = simParams(res=128,run='tracer',redshift=3.0)

  ; interpolate onto the shell
  val_interp = interpHaloShell(sP=sP,haloID=haloID,radFac=radFac,Nside=Nside)
  
  ; plot
  start_PS, sP.plotPath + 'shell_r'+str(fix(radFac*100))+'.eps'
    rot_ang = [60,-45] ;deg (left,up)
    
    title = sP.run+" "+str(sP.res)+textoidl("^3")+textoidl("\rho_{gas}  -  r / r_{vir} = "+$
      string(radFac,format='(f4.2)'))
    bartitle = textoidl("\rho_{gas}")
    
    plotMollweideProj,val_interp,rot_ang=rot_ang,title=title,bartitle=bartitle
  
  end_PS
  stop
end

