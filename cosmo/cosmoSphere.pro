; cosmoSphere.pro
; gas accretion project - interpolation and visualization of quantities onto spheres
; dnelson feb.2012

; sphereXYZCoords(): return a HEALPix subdivision at Nside resolution parameter scaled out to radius

function sphereXYZCoords, Nside=Nside, radius=radius, center=center

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

pro plotMollweideProj, data, rot_ang=rot_ang, minmax=minmax, $
                       title=title, bartitle=bartitle, pos=pos, noerase=noerase, bigbar=bigbar

  ; config
  pxsize   = 1600 ; raster resolution of mollweide map
  charSize = 1.0  ; charsize factor
  du_dv = 2.0     ; aspect ratio
  fudge = 1.02    ; spare some space around the Mollweide egg
  
  colorBar = 1    ; display color bar
  if ~keyword_set(bartitle) then colorBar = 0
  
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
  if n_elements(rot_ang) gt 0 then rot_ang = ([rot_ang,0.,0.])(0:2) else rot_ang = [0., 0., 0.]
  eul_mat = euler_matrix_new(rot_ang(0), -rot_ang(1), rot_ang(2), /Deg, /ZYX)
  
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
      u1 =  u(ellipse)
      v1 =  v(ellipse)
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
  xsize = (size(planmap))(1)
  ysize = (size(planmap))(2)

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
  
  ; "big bar" override with bottom centered colorbar
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
  loadColorTable,'helix'
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

    strmin = str(STRING(Tmin,format='(f6.2)'))
    strmax = str(STRING(Tmax,format='(f6.2)'))
    cgText,cbar_xll,cbar_yll,strmin+'  ',ALIGN=1.0,/normal,charsize=charSize
    cgText,cbar_xur,cbar_yll,'  '+strmax,ALIGN=0.0,/normal,charsize=charSize
    cgText,mean([cbar_xll,cbar_xur]),cbar_units_y,bartitle,alignment=0.5,/normal,charsize=charSize
  endif
  
end

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

; haloShellDensity(): for a given snapshot and subgroupID, evaluate the underlying density distribution 
;                     of a specified particle type on a series of radial shells

function haloShellDensity, sP=sP, partType=partType, subgroupID=subgroupID, $
                           Nside=Nside, radFacs=radFacs, save=save
  
  ; config
  nNGB   = 32  ; neighbor search in CalcHSMLds
  padFac = 4.0 ; times r_vir maximum search radius
  
  ; healpix resolution parameter, 8=768, 16~3k, 32~12k, 64~50k, 128~200k, 256~750k, 512~3M
  if ~keyword_set(Nside) then Nside = 64
  
  ; r/r_vir list of shells to compute
  if ~keyword_set(radFac) then $
    radFacs = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0]
    
  ; check for existence of a save
  if keyword_set(save) then begin
    saveFilename = sP.derivPath + 'hShells.'+sP.savPrefix+str(sP.res)+'.'+partType+'.ns'+str(Nside)+'.'+ $
                   str(sP.snap) + '.' + str(subgroupID) + '.' + str(n_elements(radFacs)) + '.sav'
                   
    if file_test(saveFilename) then begin
      restore,saveFilename
      return,r
    endif
  endif
  
  Npx      = nside2npix(Nside)
  nRadFacs = n_elements(radFacs)
  
  if padFac lt 1.9*max(radFacs) then message,'Error: Suggest increasing padFac.'
  
  ; load group catalog and find halo position and virial radius
  h  = loadSnapshotHeader(sP=sP)
  gc = loadGroupCat(sP=sP,/verbose,/skipIDs)
  
  sgpos  = subgroupPosByMostBoundID(sP=sP)
  cenPos = sgpos[*,subgroupID]
  rVir   = gc.group_r_crit200[gc.subgroupGrNr[subgroupID]]
  
  ; verify that the requested subgroupID is a primary subgroup
  priSGIDs = gcIDList(gc=gc,select='pri')
  w = where(priSGIDs eq subgroupID,countMatch)
  if ~countMatch then message,'Error: Only know how to do this for primary subgroups for now.'
  
  ; load particle positions and decide constant mass
  pos = loadSnapshotSubset(sP=sP,partType=partType,field='pos')

  if partType eq 'gas'   then massPart = sP.targetGasMass ; slightly worrying for Arepo
  if partType eq 'dm'    then massPart = h.massTable[partTypeNum(partType)]
  if partType eq 'trmc'  then massPart = sP.trMassConst
  if partType eq 'trvel' then massPart = sP.targetGasMass

  ; take conservative subset of points using periodic distances
  rad = periodicDists(cenPos,pos,sP=sP)
  
  w = where(rad le padFac*rVir,sCount)
  if ~sCount then message,'Error: No positions found near specified radius.'

  pos = pos[*,w]

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

; plotHaloShellDensComp(): compare different mass halos at one radial shell and one redshift

pro plotHaloShellDensComp

  ; config
  redshift = 1
  partType = 'gas'
  rot_ang = [0,0] ;[60,-45] ;[lat,long] center in deg (left,up)
  
  sP = simParams(res=512,run='gadget',redshift=float(redshift))

  ; get IDs of mass targets
  hMassTargets = [12.5,12.0,11.5,11.0]
  radInd = 6 ; pre-saved radFacs
  minmax = [-1.0,5.0] ; log (rho/rho_b_crit_z)

  gc = loadGroupCat(sP=sP,/skipIDs)
  priSGIDs = gcIDList(gc=gc,select='pri')
  hMasses = codeMassToLogMsun(gc.subgroupMass[priSGIDs])

  hInds = value_locate(hMasses,hMassTargets)
  w = where(hInds eq -1,count)
  if count gt 0 then hInds[w] = 0 ; largest mass halo available

  subgroupIDs = priSGIDs[hInds]

  ; plot
  start_PS, sP.plotPath + 'shell_rcomp_z'+str(redshift)+'_r'+str(radInd)+'_'+partType+'.eps', xs=6*1.5, ys=6
    
    pos = ['ul_nb','ur_nb','ll_nb','lr_nb']
    xtpos = [0.07,0.57,0.07,0.57]
    ytpos = [0.55,0.55,0.13,0.13]
    
    bartitle = "log ( "+textoidl("\rho / \rho_{crit,b,z="+string(sP.redshift,format='(i1)')+"}")+" )"
  
    for i=0,3 do begin
      ; interpolate onto the shell
      hsd_gas = haloShellDensity(sP=sP,partType=partType,subgroupID=subgroupIDs[i],/save)

      ; convert densities into ratios to the critical baryon density at this redshift
      hsd_gas.val_dens = alog10(rhoRatioToCrit(10.0^hsd_gas.val_dens,redshift=sP.redshift)) 

      healpix_data = reform(hsd_gas.val_dens[*,radInd])
      
      title = sP.run+" "+str(sP.res)+textoidl("^3")+"  z = "+string(sP.redshift,format='(f3.1)')+" "+$
              textoidl("\rho_{"+hsd_gas.partType+"} (r / r_{vir} = "+string(hsd_gas.radFacs[radInd],format='(f4.2)'))+")"

      if i eq 0 then $
        plotMollweideProj,healpix_data,rot_ang=rot_ang,title=title,bartitle=bartitle,$
          minmax=minmax,/bigbar,pos=pos[i]
      if i gt 0 then $
        plotMollweideProj,healpix_data,rot_ang=rot_ang,title="",bartitle="",$
          minmax=minmax,pos=pos[i],/noerase

      cgText,xtpos[i],ytpos[i],"M = "+string(hMasses[hInds[i]],format='(f4.1)'),$
        /normal,alignment=0.5,color=cgColor('forest green'),charsize=1.0
    endfor
  end_PS, pngResize=60, /deletePS
  
end

; plotHaloShellDensMovie(): 

pro plotHaloShellDensMovie, redshift=redshift, hMassTarget=hMassTarget
  
  ; config
  rot_ang = [0,0] ;[60,-45] ;[lat,long] center in deg (left,up)
  
  sP = simParams(res=512,run='gadget',redshift=float(redshift))  
  
  ; select a primary subgroup based on mass
  gc = loadGroupCat(sP=sP,/skipIDs)
  priSGIDs = gcIDList(gc=gc,select='pri')
  hMasses = codeMassToLogMsun(gc.subgroupMass[priSGIDs])
  
  hInd = value_locate(hMasses,hMassTarget)
  if hInd eq -1 then hInd = 0 ; largest mass halo
  
  subgroupID = priSGIDs[hInd]
  print,'selected halo ind ['+str(hInd)+'] sgID ['+str(subgroupID)+'] mass = '+string(hMasses[hInd],format='(f5.2)')
  
  ; movie configuration
  nFrames  = 450
  radStart = 0.01
  radEnd   = 1.999
  
  ; interpolate onto shells at a set of fixed radii
  print,'interpolating onto shells...'
  hsd_gas = haloShellDensity(sP=sP,partType='gas',subgroupID=subgroupID,/save)
  hsd_dm  = haloShellDensity(sP=sP,partType='dm',subgroupID=subgroupID,/save)  
  
  ; stepping in radius
  radStep    = (radEnd - radStart) / (nFrames-1)
  frameRadii = radStep * findgen(nFrames) + radStart

  ; interpolate
  print,'interpolating onto radii...'
  hpShell_gas = fltarr(hsd_gas.nPx,nFrames)
  hpShell_dm  = fltarr(hsd_dm.nPx,nFrames)
  
  for i=0,hsd_gas.nPx-1 do $
    hpShell_gas[i,*] = hermite(hsd_gas.radFacs,hsd_gas.val_dens[i,*],frameRadii)
  for i=0,hsd_dm.nPx-1 do $
    hpShell_dm[i,*] = hermite(hsd_dm.radFacs,hsd_dm.val_dens[i,*],frameRadii)

  print,'rendering frames...'
  for fn=0,nFrames-1 do begin
  ;fn = 100
    print,fn
    radius = frameRadii[fn]

    ; plot 2 vertically
    start_PS, sP.plotPath + 'frames/shell_z'+str(redshift)+'_m'+str(subgroupID)+'_r'+str(fn)+'.eps', xs=8*0.75, ys=8
      !p.multi = [0,2,1]
      
      healpix_data = reform(hpShell_gas[*,fn])
      title = sP.run+" "+str(sP.res)+textoidl("^3")+"  "+textoidl("\rho_{"+hsd_gas.partType+"} (r / r_{vir} = "+$
        string(radius,format='(f5.3)'))+")"
      bartitle = "log ( "+textoidl("\rho / \rho_{crit,b,z="+string(sP.redshift,format='(i1)')+"}")+" )"
      
      plotMollweideProj,healpix_data,rot_ang=rot_ang,title=title,bartitle=bartitle,pos='top'
      
      healpix_data = reform(hpShell_dm[*,fn])
      title = sP.run+" "+str(sP.res)+textoidl("^3")+"  "+textoidl("\rho_{"+hsd_dm.partType+"} (r / r_{vir} = "+$
        string(radius,format='(f5.3)'))+")"
      bartitle = "log ( "+textoidl("\rho / \rho_{crit,z="+string(sP.redshift,format='(i1)')+"}")+" )"
      
      plotMollweideProj,healpix_data,rot_ang=rot_ang,title=title,bartitle=bartitle,pos='bottom'
    
      ; halo mass and redshift label
      cgText,0.9,0.06,"log(M) = "+string(hMasses[hInd],format='(f4.1)'),$
        /normal,alignment=0.5,color=cgColor('forest green'),charsize=1.0
      cgText,0.9,0.03,"z = "+string(redshift,format='(f3.1)'),$
        /normal,alignment=0.5,color=cgColor('forest green'),charsize=1.0
    
    end_PS, pngResize=60, /deletePS
  endfor

end

; haloShellAngPowerSpec(): compute and plot the angular power spectrum
;                          e.g. the variation of the spherical harmonic coefficients

pro haloShellAngPowerSpec

  ; config
  Nside      = 256
  radFacs    = [0.5]
  subgroupID = 0
  
  sP = simParams(res=512,run='gadget',redshift=3.0)  
 
  ; interpolate onto shells at a set of fixed radii
  hsd_gas = haloShellDensity(sP=sP,partType='gas',subgroupID=subgroupID,Nside=Nside,radFacs=radFacs)
  hsd_dm  = haloShellDensity(sP=sP,partType='dm',subgroupID=subgroupID,Nside=Nside,radFacs=radFacs)
  
  ; choose maximum spherical harmonic order l_max
  ;l_max = fix(3.0 * hsd_gas.nSide - 1)
  l_max = fix(2.0 * hsd_gas.nSide)
  
  ; for this radius, calculate angular power spectrum 
  healpix_data = reform(hsd_gas.val_dens[*,0])
  ianafast,healpix_data,cl_gas,nlmax=l_max,/cxx,/nested,tmpdir='/n/home07/dnelson/',$
    /silent;,/won,iter_order=2

  healpix_data = reform(hsd_dm.val_dens[*,0])
  ianafast,healpix_data,cl_dm,nlmax=l_max,/cxx,/nested,tmpdir='/n/home07/dnelson/',$
    /silent;,/won,iter_order=2
    
  ; plot
  start_PS, sP.plotPath + 'powerspec.shell.eps'
    
    l_vals = findgen(l_max) + 2 ;?    
    
    xrange = [1,max(l_vals)]
    yrange = [0,10]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/xlog,$
      ytitle=textoidl("l(l+1)C_l/2\pi"),xtitle=textoidl("l"),$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+")"
      

    cgPlot,l_vals,l_vals*(l_vals+1)*cl_gas/2/!pi,psym=-8,color=getColor(1),/overplot
    cgPlot,l_vals,l_vals*(l_vals+1)*cl_dm/2/!pi,psym=-8,color=getColor(2),/overplot
      
  end_PS
stop

end

