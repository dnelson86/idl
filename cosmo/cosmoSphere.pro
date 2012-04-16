; cosmoSphere.pro
; gas accretion project - interpolation of quantities onto healpix spheres
; dnelson apr.2012

; sphereXYZCoords(): return a HEALPix subdivision at Nside resolution parameter scaled out to radius

function sphereXYZCoords, Nside=Nside, radius=radius, center=center

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; init HEALPix
  init_healpix
  
  w = where(!HEALPIX.Nside eq Nside,count)
  if (count ne 1) then begin print,'ERROR: Bad Nside.' & stop & endif
  
  ; set parameters
  Npix = nside2npix(Nside)
  
  ; get (x,y,z) positions
  pxIDs = lindgen(Npix)
  
  pix2vec_nest, Nside, pxIDs, vec
  
  ; rescale to radius
  if keyword_set(radius) then $
    vec *= radius
  
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
  
  ; debug
  ;ptRads = fltarr(Npix)
  ;for i=0,Npix-1 do ptRads[i] = sqrt(vec[i,0]^2.0 + vec[i,1]^2.0 + vec[i,2]^2.0)
  ;start_PS,'test.eps'
  ;  fsc_plot,vec[*,0],vec[*,1],xtitle="x",ytitle="y",aspect=1.0,psym=3
  ;end_PS
  
  return, vec
end

; plotMollweideProj(): plot the all-sky data in mollweide projection using some of the healpix tools
; 
; rot_ang = [lat,long]
; minmax : constrain minmax used for color table
; pos = 'bottom'/'top' : tile x2 in the vertical direction
; pos = 'ul'/'ur'/'ll'/'lr' : tile x2 in both horizontal and vertical directions
; noerase : set for all plots after the first for a compound plot

pro plotMollweideProj, data, rot_ang=rot_ang, minmax=minmax, ctName=ctName, $
                       title=title, bartitle=bartitle, pos=pos, noerase=noerase, bigbar=bigbar

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  pxsize   = 1600 ; raster resolution of mollweide map
  charSize = 1.0  ; charsize factor
  du_dv = 2.0     ; aspect ratio
  fudge = 1.02    ; spare some space around the Mollweide egg
  
  colorBar = 1    ; display color bar
  if ~keyword_set(bartitle) then colorBar = 0
  if ~keyword_set(ctName) then ctName = 'helix'
  
  ; init
  defsysv, '!healpix', exists = exists
  if (exists ne 1) then init_healpix
  
  @viewcom ; define common
  data_plot = 0 ; empty common array
  loadsky ; cgis package routine, define rotation matrices
  
  ;fits = ""       ; output healpix format fits file
  pix_type = 'N'  ; nested pixel configuration
  
  ; parse data
  ; ----------
  dim1 = n_elements(data[*,0])
  dim2 = n_elements(data[0,*])

  npix = (size(data))[1]
  nside = npix2nside(npix)
  
  ; setup rotation
  if n_elements(rot_ang) gt 0 then rot_ang = ([rot_ang,0.,0.])[0:2] else rot_ang = [0., 0., 0.]
  eul_mat = euler_matrix_new(rot_ang[0], -rot_ang[1], rot_ang[2], /Deg, /ZYX)
  
  ; minmax specified?
  if keyword_set(minmax) then begin
    min_set = minmax[0]
    max_set = minmax[1]
  endif
  
  ; convert data to mollweide projection
  ; ------------------------------------
  bad_data = !healpix.bad_value
  
  !P.BACKGROUND = 1               ; white background
  !P.COLOR = 0                    ; black foreground  
  
  if DEFINED(pxsize) then xsize = long(pxsize>200) else xsize = 800L
  ysize = xsize/2L
  n_uv  = xsize*ysize
  small_file = (n_uv GT npix) && ~keyword_set(fits)
  
  if (small_file) then begin
    ; file smaller than final map, make costly operation on the file
    ; initial data is destroyed and replaced by color
    data_plot = data
    ; color observed pixels
    find_min_max_valid, data[*,0], mindata, maxdata, valid=Obs, bad_data= 0.9 * bad_data
    data = COLOR_MAP(data, mindata, maxdata, Obs, color_bar = color_bar, mode=mode_col, $
                     minset = min_set, maxset = max_set, /silent )
  
    if defined(Obs) then Obs = 0
    Tmin = mindata & Tmax = maxdata
    planmap = MAKE_ARRAY(/BYTE, xsize, ysize, 1, Value = !P.BACKGROUND) ; white
  endif else begin ; large file
    print,'Warning: using large file strategy'
    planmap = MAKE_ARRAY(/FLOAT, xsize, ysize, 1, Value = bad_data) 
    plan_off = 0L
  endelse  

  ; generate the (u,v) position on the mollweide map
  xll= 0 & xur =  xsize-1
  yll= 0 & yur =  ysize-1
  xc = 0.5*(xll+xur) & dx = (xur - xc)
  yc = 0.5*(yll+yur) & dy = (yur - yc)
  
  yband = LONG(5.e5 / FLOAT(xsize))
  for ystart = 0, ysize - 1, yband do begin 
    yend   = (ystart + yband - 1) < (ysize - 1)
    nband = yend - ystart + 1
    u = FINDGEN(xsize)     # REPLICATE(1,nband)
    v = REPLICATE(1,xsize) # (FINDGEN(nband) + ystart)
    u =  du_dv*(u - xc)/(dx/fudge)   ; in [-2,2]*fudge
    v =        (v - yc)/(dy/fudge)   ; in [-1,1] * fudge
  
    ; for each point on the mollweide map find the corresponding position vector on the sphere
    ellipse  = where( (u^2/4. + v^2) LE 1. , nellipse)
    if (~small_file) then begin
      off_ellipse = where( (u^2/4. + v^2) GT 1. , noff_ell)
      if (noff_ell NE 0) then plan_off = [plan_off, ystart*xsize+off_ellipse]
    endif
    if (nellipse gt 0) then begin
      u1 =  u[ellipse]
      v1 =  v[ellipse]
      u = 0 & v = 0
      s1 =  sqrt( (1-v1)*(1+v1) )
      a1 =  asin(v1)
  
      z = 2.0/!PI * ( a1 + v1*s1)
      phi = -1.0 * !Pi/2.0 * u1/s1 ; lon in [-pi,pi], the minus sign is here to fit astro convention
      sz = sqrt( (1.0 - z)*(1. + z) )
      vector = [[sz * cos(phi)], [sz * sin(phi)], [z]]
      u1 = 0 & v1 = 0 & s1 = 0 & a1 = 0 & z = 0 & phi = 0 & sz = 0
  
      ; rotation
      if n_elements(rot_ang) gt 0 then vector = vector # eul_mat
      ; convert position on the sphere into pixel number and project the corresponding data value on the map
      case pix_type of
        'R' : VEC2PIX_RING, nside, vector, id_pix ; Healpix ring
        'N' : VEC2PIX_NEST, nside, vector, id_pix ; Healpix nest
        else : print,'error on pix_type'
      endcase
      if (small_file) then begin ; (data is already rescaled and color coded)
        planmap[ystart*xsize+ellipse] = sample_sparse_array(data,id_pix,in_pix=pixel_list,default=2B) ; temperature
      endif else begin ; (large file : do the projection first)
        planmap[ystart*xsize+ellipse] = sample_sparse_array(data,id_pix,in_pix=pixel_list,default=!healpix.bad_value) ; temperature
      endelse
    endif
    ellipse = 0 & id_pix = 0
  endfor
  
  ; write fits if requested
  ;if keyword_set(fits) then begin 
  ;  reso_arcmin = 60.d0 * 360.d0/(xsize-1) * fudge
  ;  reso_arcmin *=  sqrt(8.d0) / !dpi ; WCS convention, ellipse surface is 4Pi
  ;  proj2fits,planmap,fits,projection='MOLL', $
  ;            rot=rot_ang,coord='G',reso=reso_arcmin,unit=sunits,min=mindata,max=maxdata
  ;endif
  
  if (small_file) then begin
      data = 0 & pol_data = 0
  endif else begin
    ; file larger than final map, make costly coloring operation on the Mollweide map
    data_plot = temporary(data)
    pol_data = 0
    find_min_max_valid, planmap, mindata, maxdata, valid= Obs, bad_data = 0.9 * bad_data
    ; same for truecolors=1 and false colors:
    planmap = COLOR_MAP(planmap, mindata, maxdata, Obs, color_bar = color_bar, mode=mode_col, $
                        minset = min_set, maxset = max_set, /silent)
    planmap[plan_off+n_uv] = !p.background ; white
    Obs = 0 & plan_off = 0
    Tmin = mindata & Tmax = maxdata
  endelse
  
  ; setup plot
  ; ----------
  xsize = (size(planmap))[1]
  ysize = (size(planmap))[2]

  xc = 0.5*(xsize-1) & delta_x = (xsize-1 - xc)
  yc = 0.5*(ysize-1) & delta_y = (ysize-1 - yc)
  
  ; x and y range of egg
  umin = - du_dv * fudge & umax = du_dv * fudge
  vmin = - fudge         & vmax =         fudge
  
  ; position of the egg in the final window: single plot?
  if ~keyword_set(pos) then begin
    w_xll = 0.0 & w_xur = 1.0 & w_yll = 0.1 & w_yur = 0.9
    x_title = 0.5 & y_title = 0.90
    pxsize = 1.0
    cbar_units = 3.0
  endif else begin
    ; configure for multiple plots on the same page
    if pos eq 'top' then begin
      w_xll = 0.0 & w_xur = 1.0 & w_yll = 0.58 & w_yur = 0.98
      x_title = 0.5 & y_title = 0.97
      pxsize = 1.0
      cbar_units = 2.0
    endif
    if pos eq 'bottom' then begin
      w_xll = 0.0 & w_xur = 1.0 & w_yll = 0.08 & w_yur = 0.48
      x_title = 0.5 & y_title = 0.47
      pxsize = 1.0
      cbar_units = 2.0
    endif
    
    ; 2x2 with space for an independent color bar for each
    if pos eq 'ul' then begin
      w_xll = 0.0 & w_xur = 0.5 & w_yll = 0.54 & w_yur = 0.94
      x_title = 0.5 & y_title = 0.95
      pxsize = 0.5
      cbar_units = 3.0
    endif
    if pos eq 'ur' then begin
      w_xll = 0.5 & w_xur = 1.0 & w_yll = 0.54 & w_yur = 0.94
      x_title = 0.5 & y_title = 0.95
      pxsize = 0.5
      cbar_units = 3.0
    endif
    if pos eq 'll' then begin
      w_xll = 0.0 & w_xur = 0.5 & w_yll = 0.07 & w_yur = 0.47
      x_title = 0.5 & y_title = 0.95
      pxsize = 0.5
      cbar_units = 3.0
    endif
    if pos eq 'lr' then begin
      w_xll = 0.5 & w_xur = 1.0 & w_yll = 0.07 & w_yur = 0.47
      x_title = 0.5 & y_title = 0.95
      pxsize = 0.5
      cbar_units = 3.0
    endif
    
    ; 2x2 compact tiling with room for the colorbar at the bottom
    if pos eq 'ul_nb' then begin
      w_xll = 0.0 & w_xur = 0.5 & w_yll = 0.54 & w_yur = 0.94
      x_title = 0.5 & y_title = 0.95
      pxsize = 0.5
      cbar_units = 3.0
    endif
    if pos eq 'ur_nb' then begin
      w_xll = 0.5 & w_xur = 1.0 & w_yll = 0.54 & w_yur = 0.94
      x_title = 0.5 & y_title = 0.95
      pxsize = 0.5
      cbar_units = 3.0
    endif
    if pos eq 'll_nb' then begin
      w_xll = 0.0 & w_xur = 0.5 & w_yll = 0.12 & w_yur = 0.52
      x_title = 0.5 & y_title = 0.95
      pxsize = 0.5
      cbar_units = 3.0
    endif
    if pos eq 'lr_nb' then begin
      w_xll = 0.5 & w_xur = 1.0 & w_yll = 0.12 & w_yur = 0.52
      x_title = 0.5 & y_title = 0.95
      pxsize = 0.5
      cbar_units = 3.0
    endif
  endelse
  
  w_dx = w_xur - w_xll
  w_dy = w_yur - w_yll
  w_dx_dy = w_dx / w_dy  ; 1/0.8
  
  ; color bar, position, dimension
  cbar_dx = 1.0/3.0
  cbar_dy = 1.0/70.0
  
  if keyword_set(pos) then begin
    ; for 2x2 tiling reduce the size of the individual colorbars below each egg
    if pos eq 'ul' or pos eq 'ur' or pos eq 'll' or pos eq 'lr' then begin
      cbar_dx /= 1.5
      cbar_dy /= 1.5
    endif
  endif 
  
  cbar_xll = (w_xll+w_xur)/2.0 - cbar_dx/2.0
  cbar_xur = (w_xll+w_xur)/2.0 + cbar_dx/2.0
  cbar_yur = w_yll - cbar_dy
  cbar_yll = cbar_yur - cbar_dy
  
  cbar_units_y = cbar_yll-cbar_dy*cbar_units
  
  ; "big bar" override with bottom centered colorbar (assume no individual colorbars for spacing)
  if keyword_set(bigbar) then begin
    cbar_dx = 1.0/2.0
    cbar_dy = 1.0/50.0
    cbar_xll = 0.5 - cbar_dx/2.0
    cbar_xur = 0.5 + cbar_dx/2.0
    cbar_yur = 0.1 - cbar_dy
    cbar_yll = cbar_yur - cbar_dy
    
    cbar_units_y = cbar_yll-cbar_dy*2
  endif
  
  ; load color table
  loadColorTable,ctName
  tvlct,red,green,blue,/get
  
  ; set up some specific color definitions: reserve first colors for Black, White and Neutral grey
  idx_black = 0B & idx_white = 1B   & idx_grey = 2B   & idx_bwg = [idx_black, idx_white, idx_grey]
  col_black = 0B & col_white = 255B & col_grey = 175B & col_bwg = [col_black, col_white, col_grey]
  red  [idx_bwg] = col_bwg
  green[idx_bwg] = col_bwg
  blue [idx_bwg] = col_bwg
  TVLCT,red,green,blue

  ; make the plot
  cgPlot, /nodata, [umin,umax], [vmin,vmax], pos=[w_xll,w_yll,w_xur,w_yur], xs=5, ys=5, noerase=noerase
  tv, planmap,w_xll,w_yll,/normal,xsize=pxsize

  ; title
  if keyword_set(title) then cgText, x_title, y_title ,title, align=0.5, /normal  
  
  ;  the color bar
  if keyword_set(colorBar) then begin
    color_bar_out = BYTE(CONGRID(color_bar,xsize*cbar_dx)) # REPLICATE(1.,(ysize*cbar_dy*w_dx_dy)>1)
    ;back = REPLICATE(BYTE(!P.BACKGROUND),xsize,(ysize*cbar_dy*w_dx_dy)>1)
    ;back(xsize*cbar_xll,0) = color_bar_out
    TV,color_bar_out,cbar_xll,cbar_yll,/normal,xsize = cbar_dx

    strmin = str(string(Tmin,format='(f6.2)'))
    strmax = str(string(Tmax,format='(f6.2)'))
    
    ; different integer formatting for diverging large bars (non-normalized e.g. radvel)
    if Tmin lt 0 and Tmax gt 100 then begin
      strmin = str(string(Tmin,format='(i5)'))
      strmax = str(string(Tmax,format='(i5)'))
    endif
    
    cgText,cbar_xll,cbar_yll,strmin+'  ',ALIGN=1.0,/normal,charsize=charSize
    cgText,cbar_xur,cbar_yll,'  '+strmax,ALIGN=0.0,/normal,charsize=charSize
    cgText,mean([cbar_xll,cbar_xur]),cbar_units_y,bartitle,alignment=0.5,/normal,charsize=charSize
  endif
  
end

; haloShellDensity(): for a given snapshot and subgroupID, evaluate the underlying density distribution 
;                     of a specified particle type on a series of radial shells
;
; cutSubS = cut substructures (satellite subgroups) out before estimating densities

function haloShellDensity, sP=sP, partType=partType, subgroupID=subgroupID, $
                           Nside=Nside, radFacs=radFacs, save=save, cutSubS=cutSubS
  
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  nNGB   = 32  ; neighbor search in CalcHSMLds
  padFac = 4.0 ; times r_vir maximum search radius
  
  ; healpix resolution parameter, 8=768, 16~3k, 32~12k, 64~50k, 128~200k, 256~750k, 512~3M
  if ~keyword_set(Nside) then Nside = 64
  
  ; r/r_vir list of shells to compute
  if ~keyword_set(radFac) then $
    radFacs = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0]
    
  ; by construction this assumes constant mass particles, bad idea for arepo
  if sP.trMCPerCell ne 0 then message,'Error: Use haloShellValue for arepo density estimates.'
    
  ; check for existence of a save
  csTag = ''
  if keyword_set(cutSubS) then csTag = '.cutSubS'
  saveFilename = sP.derivPath+'hShells/hShells.'+sP.savPrefix+str(sP.res)+'.'+partType+csTag+'.ns'+$
                 str(Nside)+'.'+str(sP.snap)+'.'+str(subgroupID)+'.'+str(n_elements(radFacs)) + '.sav'
                 
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif
  
  Npx      = nside2npix(Nside)
  nRadFacs = n_elements(radFacs)
  
  if padFac lt 1.9*max(radFacs) then message,'Error: Suggest increasing padFac.'
  
  ; load group catalog and find halo position and virial radius
  h  = loadSnapshotHeader(sP=sP)
  gc = loadGroupCat(sP=sP,/verbose,/readIDs)
  
  sgpos  = subgroupPosByMostBoundID(sP=sP)
  cenPos = sgpos[*,subgroupID]
  rVir   = gc.group_r_crit200[gc.subgroupGrNr[subgroupID]]
  
  ; verify that the requested subgroupID is a primary subgroup
  priSGIDs = gcIDList(gc=gc,select='pri')
  w = where(priSGIDs eq subgroupID,countMatch)
  if ~countMatch then message,'Error: Only know how to do this for primary subgroups for now.'
  
  ; load particle positions and decide constant mass
  pos = loadSnapshotSubset(sP=sP,partType=partType,field='pos')

  if partType eq 'gas'   then massPart = sP.targetGasMass
  if partType eq 'dm'    then massPart = h.massTable[partTypeNum(partType)]
  if partType eq 'trmc'  then massPart = sP.trMassConst
  if partType eq 'trvel' then massPart = sP.targetGasMass

  ; take conservative subset of points using periodic distances
  rad = periodicDists(cenPos,pos,sP=sP)
  
  w = where(rad le padFac*rVir,sCount)
  if ~sCount then message,'Error: No positions found near specified radius.'

  pos = pos[*,w]
  
  ; if cutting substructure, load particle ids and make same radial cut
  if keyword_set(cutSubS) then begin
    ids = loadSnapshotSubset(sP=sP,partType=partType,field='ids')
    ids = ids[w]
    
    ; make a list of satellites of this halo
    nSubs    = gc.groupNSubs[gc.subgroupGrNr[subgroupID]]
    firstSub = gc.groupFirstSub[gc.subgroupGrNr[subgroupID]]
    
    if firstSub ne subgroupID then message,'Warning: firstSub is not subgroupID'

    satGCids = indgen(nSubs-1) + firstSub + 1

    ; make a list of member particle ids of these satellites for the requested particle type
    satPIDs = gcPIDList(gc=gc,select='secondary',valGCids=satGCids,partType=partType)
    gc = !NULL
    
    ; remove the intersection of (satPIDs,ids) from pos
    match,satPIDs,ids,sat_ind,ids_ind,count=count,/sort
    sat_ind = !NULL
    satPIDs = !NULL
    
    all = bytarr(n_elements(ids))
    if count gt 0 then all[ids_ind] = 1B
    w = where(all eq 0B, ncomp)
    
    ids_ind = !NULL
    ids     = !NULL
    
    print,'Substructures cut ['+str(count)+'] of ['+str(n_elements(ids))+'] have left: '+str(ncomp)
    if ncomp gt 0 then pos = pos[*,w]
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
        val_dens   : fltarr(Npx,nRadFacs)   }

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
  r.val_dens = alog10(estimateDensityTophat(pos,pos_search=sphereXYZ,mass=massPart,$
                                            ndims=3,nNGB=nNGB,boxSize=sP.boxSize))
  
  if keyword_set(save) then begin
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
  endif
  
  return, r
end

; haloShellValue(): for a given snapshot and subgroupID, evaluate some property/function value (valName)
;                   of a specified particle type on a series of radial shells (radFacs)
;
; cutSubS = cut substructures (satellite subgroups) out before estimating densities

function haloShellValue, sP=sP, partType=partType, valName=valName, subgroupID=subgroupID, $
                         Nside=Nside, radFacs=radFacs, save=save, cutSubS=cutSubS
                         
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  nNGB   = 32  ; neighbor search in CalcTHVal
  padFac = 4.0 ; times r_vir maximum search radius
  
  ; healpix resolution parameter, 8=768, 16~3k, 32~12k, 64~50k, 128~200k, 256~750k, 512~3M
  if ~keyword_set(Nside) then Nside = 64
  
  ; r/r_vir list of shells to compute
  if ~keyword_set(radFac) then $
    radFacs = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0]
    
  if ~keyword_set(valName) then message,'Error: Must specify valName'
    
  ; check for existence of a save
  csTag = ''
  if keyword_set(cutSubS) then csTag = '.cutSubS'
  
  saveFilename = sP.derivPath+'hShells/hShells.'+valName+'.'+sP.savPrefix+str(sP.res)+'.'+partType+csTag+$
                 '.ns'+str(Nside)+'.'+str(sP.snap)+'.'+str(subgroupID)+'.'+str(n_elements(radFacs)) + '.sav'
                 
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif
  
  Npx      = nside2npix(Nside)
  nRadFacs = n_elements(radFacs)
  
  if padFac lt 1.9*max(radFacs) then message,'Error: Suggest increasing padFac.'
  
  ; "radial mass flux" as density * radial velocity (area normalization / r^2 sphere factor omitted)
  ; instead of estimating this value on each particle and then tophat smoothing, we use the
  ; estimates of each avaiable already on each point on the sphere
  if valName eq 'radmassflux' then begin
    hsv_dens   = haloShellValue(sP=sP,partType=partType,valName='density',subgroupID=subgroupID,$
                                Nside=Nside,radFacs=radFacs,cutSubS=cutSubS,save=save)
    hsv_radvel = haloShellValue(sP=sP,partType=partType,valName='radvel',subgroupID=subgroupID,$
                                Nside=Nside,radFacs=radFacs,cutSubS=cutSubS,save=save)
                                
    hsv_dens.valName = 'radmassflux'
    hsv_dens.value = alog10(hsv_dens.value * units.UnitMass_in_Msun) * (hsv_radvel.value * units.kmS_in_kpcYr)
    return, hsv_dens ; Msun / kpc^2 * year
  endif
  
  ; load group catalog and find halo position and virial radius
  h  = loadSnapshotHeader(sP=sP)
  if keyword_set(cutSubS) then gc = loadGroupCat(sP=sP,/verbose,/readIDs) $
  else gc = loadGroupCat(sP=sP,/verbose,/skipIDs)
  
  sgpos  = subgroupPosByMostBoundID(sP=sP)
  cenPos = sgpos[*,subgroupID]
  rVir   = gc.group_r_crit200[gc.subgroupGrNr[subgroupID]]
  
  ; verify that the requested subgroupID is a primary subgroup
  priSGIDs = gcIDList(gc=gc,select='pri')
  w = where(priSGIDs eq subgroupID,countMatch)
  if ~countMatch then message,'Error: Only know how to do this for primary subgroups for now.'
  
  ; load particle positions
  pos = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
  
  ; take conservative subset of points using periodic distances
  rad = periodicDists(cenPos,pos,sP=sP)
  
  wRadCut = where(rad le padFac*rVir,sCount)
  if ~sCount then message,'Error: No positions found near specified radius.'
  rad = !NULL
  
  pos = pos[*,wRadCut]
  
  ; load data value requested and make posval = [[pos],[val]]
  if valName eq 'temp' then begin
    u = loadSnapshotSubset(sP=sP,partType=partType,field='u')
    u = u[wRadCut]
    
    nelec = loadSnapshotSubset(sP=sP,partType=partType,field='nelec')
    nelec = nelec[wRadCut]
    
    value = reform(convertUtoTemp(u,nelec,/log),[1,n_elements(wRadCut)]) ; make 1xN vector
    
    posval = [pos,value]
    thMode = 1 ; mean
  endif
  
  if valName eq 'radvel' then begin
    vel = loadSnapshotSubset(sP=sP,partType=partType,field='vel')
    vel = vel[*,wRadCut]
    
    ; make normalized position vector wrt halo center = vec(r) / ||r||
    rnorm0 = reform(cenPos[0] - pos[0,*])
    rnorm1 = reform(cenPos[1] - pos[1,*])
    rnorm2 = reform(cenPos[2] - pos[2,*])
    
    correctPeriodicDistVecs, rnorm0, sP=sP
    correctPeriodicDistVecs, rnorm1, sP=sP
    correctPeriodicDistVecs, rnorm2, sP=sP
    
    ; denominator and do divide
    rnorm = sqrt(rnorm0*rnorm0 + rnorm1*rnorm1 + rnorm2*rnorm2)

    rnorm0 /= rnorm
    rnorm1 /= rnorm
    rnorm2 /= rnorm
    rnorm = !NULL
    
    ; dot(vel,rnorm) gives the magnitude of the projection of vel onto vec(r)
    value = vel[0,*]*rnorm0 + vel[1,*]*rnorm1 + vel[2,*]*rnorm2 ; 1xN
    rnorm0 = !NULL
    rnorm1 = !NULL
    rnorm2 = !NULL
    
    posval = [pos,value]
    thMode = 1 ; mean
  endif
  
  if valName eq 'density' then begin
    mass = loadSnapshotSubset(sP=sP,partType=partType,field='mass')
    mass = reform(mass[wRadCut],[1,n_elements(wRadCut)]) ; 1xN
    
    posval = [pos,mass]
    thMode = 3 ; total/volume
  endif
  
  if n_elements(posval) eq 0 then message,'Error: Unrecognized value name.'
  
  ; if cutting substructure, load particle ids and make same radial cut
  if keyword_set(cutSubS) then begin
    ids = loadSnapshotSubset(sP=sP,partType=partType,field='ids')
    ids = ids[wRadCut]
    
    ; make a list of satellites of this halo
    nSubs    = gc.groupNSubs[gc.subgroupGrNr[subgroupID]]
    firstSub = gc.groupFirstSub[gc.subgroupGrNr[subgroupID]]
    
    if firstSub ne subgroupID then message,'Warning: firstSub is not subgroupID'

    satGCids = indgen(nSubs-1) + firstSub + 1

    ; make a list of member particle ids of these satellites for the requested particle type
    satPIDs = gcPIDList(gc=gc,select='secondary',valGCids=satGCids,partType=partType)
    gc = !NULL
    
    ; remove the intersection of (satPIDs,ids) from pos
    match,satPIDs,ids,sat_ind,ids_ind,count=count,/sort
    sat_ind = !NULL
    satPIDs = !NULL
    
    all = bytarr(n_elements(ids))
    if count gt 0 then all[ids_ind] = 1B
    wSubSComp = where(all eq 0B, ncomp)
    
    ids_ind = !NULL
    ids     = !NULL
    
    print,'Substructures cut ['+str(count)+'] of ['+str(n_elements(ids))+'] have left: '+str(ncomp)
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
  r.value = CalcTHVal(posval,sphereXYZ,ndims=3,nNGB=nNGB,thMode=thMode,boxSize=sP.boxSize)

  if keyword_set(save) then begin
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
  endif
  
  return, r
end
