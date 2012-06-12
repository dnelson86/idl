; cosmoSphere.pro
; gas accretion project - interpolation of quantities onto healpix spheres
; dnelson apr.2012

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
    
    ; 2x2 same as above but no title
    if pos eq 'ul_nt' then begin
      w_xll = 0.0 & w_xur = 0.5 & w_yll = 0.6 & w_yur = 1.0
      x_title = 0.5 & y_title = 0.95
      pxsize = 0.5
      cbar_units = 3.0
    endif
    if pos eq 'ur_nt' then begin
      w_xll = 0.5 & w_xur = 1.0 & w_yll = 0.6 & w_yur = 1.0
      x_title = 0.5 & y_title = 0.95
      pxsize = 0.5
      cbar_units = 3.0
    endif
    if pos eq 'll_nt' then begin
      w_xll = 0.0 & w_xur = 0.5 & w_yll = 0.10 & w_yur = 0.50
      x_title = 0.5 & y_title = 0.95
      pxsize = 0.5
      cbar_units = 3.0
    endif
    if pos eq 'lr_nt' then begin
      w_xll = 0.5 & w_xur = 1.0 & w_yll = 0.10 & w_yur = 0.50
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
  cbar_dy = 1.0/60.0 ;70
  
  if keyword_set(pos) then begin
    ; for 2x2 tiling reduce the size of the individual colorbars below each egg
    if pos eq 'ul' or pos eq 'ur' or pos eq 'll' or pos eq 'lr' then begin
      cbar_dx /= 1.5
      cbar_dy /= 1.0
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

    strmin = str(string(Tmin,format='(f5.1)'))
    strmax = str(string(Tmax,format='(f5.1)'))
    
    ; different integer formatting for diverging large bars (non-normalized e.g. radvel)
    if Tmin lt 0 and Tmax gt 100 then begin
      strmin = str(string(Tmin,format='(i5)'))
      strmax = str(string(Tmax,format='(i5)'))
    endif
    
    cgText,cbar_xll,cbar_yll,strmin+'   ',ALIGN=1.0,/normal,charsize=charSize
    cgText,cbar_xur,cbar_yll,'   '+strmax,ALIGN=0.0,/normal,charsize=charSize
    cgText,mean([cbar_xll,cbar_xur]),cbar_units_y,textoidl(bartitle),alignment=0.5,/normal,charsize=charSize
  endif
  
end


; angularFoF(): angular friends of friends algorithm on a healpix map mask

function angularFoF, healpix_mask, discrete=discrete, verbose=verbose

  if ~keyword_set(discrete) then message,'Error: Only discrete (touching pixels) implemented.'
  
  fragCutoff = n_elements(healpix_mask) / 1500 ; ~30, minimum number of pixels to save group
  
  ; group membership mask
  if max(healpix_mask) eq 0 then message,'Error: Empty mask.'
  fofMemMask = intarr(n_elements(healpix_mask))
  
  Nside = npix2nside(n_elements(healpix_mask))
  if Nside eq -1 then message,'Error: Bad pixel count in mask'
  
  nextGroupID = 1
  
  ; find starting orphans
  orphanInds = where(healpix_mask eq 1B,orphanCount)
    
  ; process each orphan pixel once
  for i=0L,orphanCount-1 do begin
    pxInd = orphanInds[i]
    
    ; locate immediately neighboring pixels in mask
    neighbours_nest,Nside,pxInd,neighborInds,nNeighbors
    neighborInds = neighborInds[where(healpix_mask[neighborInds] eq 1B,nnCount)]
    
    ; candidate groupID as minimum of all neighbor+self groupIDs, or next available if all zero
    candID = fofMemMask[pxInd]
    if nnCount gt 0 then candID = min([candID,fofMemMask[neighborInds]])
    
    if candID eq 0 then begin
      candID = nextGroupID
      nextGroupID += 1
    endif
    
    ; assign groupID
    if nnCount gt 0 then fofMemMask[neighborInds] = candID
    fofMemMask[pxInd] = candID
    
  endfor
  
  ; determine number of fof groups and merge
  nGroups = max(fofMemMask)
  collisionFlag = 1
  count = 0L
  
  if keyword_set(verbose) then $
    print,'Found ['+str(nGroups)+'] total groups, need to merge...'

  ; loop while there are potentially touching groups
  while collisionFlag do begin
    count += 1
    if keyword_set(verbose) then print,count
    collisionFlag = 0
    
    ; process each grouped (originally orphaned) pixel
    for i=0L,orphanCount-1 do begin
      pxInd = orphanInds[i]
      
      ; locate immediately neighboring pixels and their groupIDs
      neighbours_nest,Nside,pxInd,neighborInds,nNeighbors
      neighborInds = neighborInds[where(healpix_mask[neighborInds] eq 1B,nnCount)]
      
      colGroupIDs = fofMemMask[pxInd]
      if nnCount gt 0 then colGroupIDs = [colGroupIDs,fofMemMask[neighborInds]]
      
      ; if there is more than one unique value, flag collision and replace all with minimum
      if nuniq(colGroupIDs) gt 1 then begin
        collisionFlag = 1
        candID = min(colGroupIDs)
        fofMemMask[pxInd] = candID
        if nnCount gt 0 then fofMemMask[neighborInds] = candID
      endif
    endfor
  endwhile
  
  ; now we have the minimal set, compress group IDs so they are [1,...,nGroups]
  oldGroupIDs = uniqvals(fofMemMask[where(healpix_mask eq 1B)])
  nextGroup = 1
  
  foreach oldGID,oldGroupIDs do begin
    w = where(fofMemMask eq oldGID,count)
    if count eq 0 then message,'Error'

    ; discard small fragments
    if count lt fragCutoff then begin
      fofMemMask[w] = 0
    endif else begin
      ; group size large enough, keep
      fofMemMask[w] = nextGroup
      nextGroup += 1
    endelse
  endforeach
  
  if keyword_set(verbose) then $
    print,'Found ['+str(max(fofMemMask))+'] final groups.'
  
  return, fofMemMask  
end

; haloFilamentCrossSec():

function haloFilamentCrossSec, sP=sP, subgroupID=subgroupID, verbose=verbose, radIndOffset=radIndOffset

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config  
  threshOverDensAll  = 1.5 ; minimum rho/meanrho for candidate filament pixel
  threshOverDensRvir = 2.0

  widthShells = 5        ; odd, how many shells to enforce overdensity over, centered at rvir
  wrapTol     = !pi/20.0 ; radians, unwrapping lat,long at edges
  querySize   = !pi/4.0  ; radians, cutout size around each filament to calculate radii to
  circCutoff  = 0.1  ; exclude highly irregular shaped cross sections
  thetaTol    = 0.05 ; cannot be too close to N/S poles for phi(long) unwrap
  
  ; load the density shells
  hsv = haloShellValue(sP=sP,partType='gas',valName='density',subgroupID=subgroupID,/cutSubS)
  
  ; transform each shell to the ratio of the shell median
  for radInd=0,hsv.nRadFacs-1 do begin
    hsv.value[*,radInd] = reform(hsv.value[*,radInd] / mean(hsv.value[*,radInd]))
  endfor
  
  ; locate r=rvir radInd
  rvirInd = where(hsv.radFacs eq 1.0,count)
  if count ne 1 then message,'Error: Failed to locate rvir shell.'
  
  minInd = rvirInd[0] - floor(widthShells/2.0)
  maxInd = rvirInd[0] + floor(widthShells/2.0)
  
  ; mask filament pixel candidates based on overdensity requirement across radii range
  odMask = bytarr(hsv.nPx) + 1B

  for radInd=minInd,maxInd do begin
    w = where(hsv.value[*,radInd] lt threshOverDensAll,count)
    if count gt 0 then odMask[w] = 0B
  endfor
  
  w = where(odMask eq 1B,count)
  if keyword_set(verbose) then $
    print,'['+str(count)+'] of ['+str(hsv.nPx)+'] global filament pixel candidates.'

  ; find filaments on one shell
  radInd = rvirInd[0] + radIndOffset
 
  cenPos     = [] ; density weighted centroid (x,y,z)
  cenPosLL   = [] ; theta,phi = colat,long
  effCirSize = [] ; effective area of each filament (rad^2)
  cirMeasure = [] ; circularity measure [0,1]  
  
  filPxInds = []
  pxNums    = []
  ckpcDists = []
  
  ; further remove from mask more stringent overdensity requirement at this shell
  w = where(hsv.value[*,radInd] lt threshOverDensRvir,count)
  if count gt 0 then odMask[w] = 0B
  
  if keyword_set(verbose) then begin
    start_PS,'testmask2_r'+str(radInd)+'.eps'
      plotMollweideProj,float(odMask),title="",bartitle="",minmax=[0,1]
    end_PS,/deletePS,pngResize=40
  endif
  
  w = where(odMask eq 1B,count)
  if keyword_set(verbose) then print,'['+str(count)+'] of ['+str(hsv.nPx)+'] filament pixel candidates.'
  if count eq 0 then message,'Error: Too stringent local cut to mask.'

  ; run angular friends of friends
  angfof = angularFoF(odMask,/discrete)
  nGroupsCand = max(angfof)
    
  ; determine a overdensity weighted center position of each filament on the shell
  for i=0,nGroupsCand-1 do begin
    pxInds = where(angfof eq i+1)
    pix2ang_nest, hsv.nSide, pxInds, pxTheta, pxPhi
    pxWts = hsv.value[pxInds,radInd]
    
    if max(pxTheta) gt !pi-thetaTol or min(pxTheta) lt thetaTol then continue ; pole skip
    
    ; detect wrap (theta in [0,pi], phi in [0,2pi]) and correct
    if min(pxTheta) lt wrapTol and max(pxTheta) gt !pi-wrapTol then begin
      w = where(pxTheta gt !pi*0.4 and pxTheta lt !pi*0.6,count)
      if count ne 0 then continue ; something odd
      
      ; move lower points up by pi
      wTheta = where(pxTheta lt !pi*0.5,countMoveTheta)
      pxTheta[wTheta] += !pi
    endif
    
    if min(pxPhi) lt wrapTol and max(pxPhi) gt 2*!pi-wrapTol then begin
      w = where(pxPhi gt 2*!pi*0.4 and pxPhi lt 2*!pi*0.6,count)
      if count ne 0 then continue ; something odd
      
      ; move lower points up by 2pi
      wPhi = where(pxPhi lt 2*!pi*0.5,countMovePhi)
      pxPhi[wPhi] += 2*!pi
    endif
    
    ; weighted mean for centroid
    wtSum = total(pxWts)
    cenTheta = total(pxTheta * pxWts) / wtSum
    cenPhi = total(pxPhi * pxWts) / wtSum
    
    ; maximum extent in both axes for effective circular size and circularity measure      
    maxExtent = gcDist([max(pxTheta),max(pxPhi)],[min(pxTheta),min(pxPhi)]) ;[lat,long]
    curECS = 2*!pi*(1-cos(maxExtent/2.0))
    curCM  = (4*!pi/hsv.nPx*n_elements(pxInds)) / curECS
    
    if keyword_set(verbose) then print,radInd,i,maxExtent,curECS,curCM
    
    if curCM lt circCutoff then continue ; circ measure skip

    ; move center back into range if necessary
    if cenTheta gt !pi then cenTheta -= !pi
    if cenPhi gt 2*!pi then cenPhi -= 2*!pi
    
    ; convert centroid to pixel to vector
    ang2pix_nest, hsv.nSide, cenTheta, cenPhi, cenPxInd
    pix2vec_nest, hsv.nSide, cenPxInd, cenVec3
    
    ; query circular disc around centroid pixel and recalculate center position (without mask bias)
    query_disc, hsv.nSide, cenVec3, maxExtent/2.0, qdPxInds, /nested
    
    pix2ang_nest, hsv.nSide, qdPxInds, pxTheta, pxPhi
    pxWts = hsv.value[qdPxInds,radInd]
    
    ; detect wrap (theta in [0,pi], phi in [0,2pi]) and correct
    if min(pxTheta) lt wrapTol and max(pxTheta) gt !pi-wrapTol then begin
      w = where(pxTheta gt !pi*0.4 and pxTheta lt !pi*0.6,count)
      if count ne 0 then continue ; something odd
      
      ; move lower points up by pi
      wTheta = where(pxTheta lt !pi*0.5,countMoveTheta)
      pxTheta[wTheta] += !pi
    endif
    
    if min(pxPhi) lt wrapTol and max(pxPhi) gt 2*!pi-wrapTol then begin
      w = where(pxPhi gt 2*!pi*0.4 and pxPhi lt 2*!pi*0.6,count)
      if count ne 0 then continue ; something odd
      
      ; move lower points up by 2pi
      wPhi = where(pxPhi lt 2*!pi*0.5,countMovePhi)
      pxPhi[wPhi] += 2*!pi
    endif
    
    ; weighted mean for centroid
    wtSum = total(pxWts)
    cenTheta2 = total(pxTheta * pxWts) / wtSum
    cenPhi2 = total(pxPhi * pxWts) / wtSum
    
    ; store center
    ang2pix_nest, hsv.nSide, cenTheta2, cenPhi2, cenPxInd
    pix2vec_nest, hsv.nSIde, cenPxInd, vec3xyz

    ; query a new selection around this center
    query_disc, hsv.nSide, vec3xyz, querySize, qdPxInds, /nested

    ; calculate distance to each pixel on the sphere in ckpc
    pix2ang_nest, hsv.nSide, qdPxInds, pxTheta, pxPhi
    
    ; detect wrap (theta in [0,pi], phi in [0,2pi]) and correct
    if min(pxTheta) lt wrapTol and max(pxTheta) gt !pi-wrapTol then begin
      w = where(pxTheta gt !pi*0.4 and pxTheta lt !pi*0.6,count)
      if count ne 0 then continue ; something odd
      
      ; move lower points up by pi
      wTheta = where(pxTheta lt !pi*0.5,countMoveTheta)
      pxTheta[wTheta] += !pi
    endif
    
    if min(pxPhi) lt wrapTol and max(pxPhi) gt 2*!pi-wrapTol then begin
      w = where(pxPhi gt 2*!pi*0.4 and pxPhi lt 2*!pi*0.6,count)
      if count ne 0 then continue ; something odd
      
      ; move lower points up by 2pi
      wPhi = where(pxPhi lt 2*!pi*0.5,countMovePhi)
      pxPhi[wPhi] += 2*!pi
    endif
    
    angDists = fltarr(n_elements(qdPxInds))
    for j=0,n_elements(qdPxInds)-1 do $
      angDists[j] = gcDist([cenTheta2,cenPhi2],[pxTheta[j],pxPhi[j]])
    
    ; store values
    effCirSize = [effCirSize, curECS]
    cirMeasure = [cirMeasure, curCM]
    cenPos = [[cenPos], [reform(float(vec3xyz) * hsv.radFacs[radInd] * hsv.rVir)]]
    cenPosLL = [[cenPosLL], [cenTheta2,cenPhi2]]
    
    filPxInds = [filPxInds, qdPxInds]
    pxNums    = [pxNums, n_elements(qdPxInds)]
    ckpcDists = [ckpcDists, angDists * hsv.rVir]
  endfor ; nGroups
  
  if keyword_set(verbose) then begin
    start_PS,'testmask2_r'+str(radInd)+'_fof.eps'
      plotMollweideProj,float(angfof),title="",bartitle="",minmax=[0,max(angfof)]
    end_PS,/deletePS,pngResize=40
  endif

  r = {}
  if n_elements(pxNums) gt 0 then $
  r = { nFilaments:n_elements(pxNums), cenPos:cenPos, cenPosLL:cenPosLL, effCirSize:effCirSize, $
        cirMeasure:cirMeasure, filPxInds:filPxInds, pxNums:pxNums, ckpcDists:ckpcDists, $
        radInd:radInd}

  return,r
end

; haloFilamentSearch(): apply filament search algorithm

function haloFilamentSearch, sP=sP, subgroupID=subgroupID, verbose=verbose

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config  
  threshOverDensAll  = 1.0 ;1.5 ; minimum rho/meanrho for candidate filament pixel
  threshOverDensRvir = 1.5 ;2.0
  
  widthShells      = 5        ; odd, how many shells to enforce overdensity over, centered at rvir
  minNumDetections = 3       ; number of detection points per filament required to fit & keep
  wrapTol          = !pi/20.0 ; radians, unwrapping lat,long at edges
  maxNGroups       = 10       ; max filament detections per shell
  circMeasCutoff   = 0.2      ; circularity measure minimum for acceptance
  maxExtentCutoff  = !pi/4.0 ; radians, maximum size of filament cross section to accept
  filDistTol       = !pi/40.0 ; radians, matching filament centroids between shells
  
  ; load the density shells
  hsv = haloShellValue(sP=sP,partType='gas',valName='density',subgroupID=subgroupID,/cutSubS)
  
  ; transform each shell to the ratio of the shell median
  for radInd=0,hsv.nRadFacs-1 do begin
    hsv.value[*,radInd] = reform(hsv.value[*,radInd] / mean(hsv.value[*,radInd]))
  endfor
  
  ; locate r=rvir radInd
  rvirInd = where(hsv.radFacs eq 1.0,count)
  if count ne 1 then message,'Error: Failed to locate rvir shell.'
  
  minInd = rvirInd[0] - floor(widthShells/2.0)
  maxInd = rvirInd[0] + floor(widthShells/2.0)
  
  ; mask filament pixel candidates based on overdensity requirement across radii range
  odMask = bytarr(hsv.nPx) + 1B

  for radInd=minInd,maxInd do begin
    w = where(hsv.value[*,radInd] lt threshOverDensAll,count)
    if count gt 0 then odMask[w] = 0B
  endfor
  
  w = where(odMask eq 1B,count)
  if keyword_set(verbose) then $
    print,'['+str(count)+'] of ['+str(hsv.nPx)+'] global filament pixel candidates.'
  
  if keyword_set(verbose) then begin
    start_PS,'testmask.eps'
      plotMollweideProj,float(odMask),title="",bartitle="",minmax=[0,1]
    end_PS,/deletePS,pngResize=40
  endif
  
  ; arrays
  radMasks   = bytarr(widthShells,hsv.nPx)
  angfof     = intarr(widthShells,hsv.nPx)
  nGroups    = intarr(widthShells)
  cenPos     = fltarr(widthShells,maxNGroups,3) ; density weighted centroid (x,y,z)
  cenPosLL   = fltarr(widthShells,maxNGroups,2) ; theta,phi = colat,long
  effCirSize = fltarr(widthShells,maxNGroups) ; effective area of each filament (rad^2)
  cirMeasure = fltarr(widthShells,maxNGroups) ; circularity measure [0,1]
  
  ; 2d source location at each plane
  k = 0
  for radInd=minInd,maxInd do begin
    radMasks[k,*] = odMask  
  
    ; further remove from mask more stringent overdensity requirement at this shell
    w = where(hsv.value[*,radInd] lt threshOverDensRvir,count)
    if count gt 0 then radMasks[k,w] = 0B
    
    if keyword_set(verbose) then begin
      start_PS,'testmask_r'+str(radInd)+'.eps'
        plotMollweideProj,float(reform(radMasks[k,*])),title="",bartitle="",minmax=[0,1]
      end_PS,/deletePS,pngResize=40
    endif
    
    w = where(radMasks[k,*] eq 1B,count)
    if keyword_set(verbose) then print,'['+str(count)+'] of ['+str(hsv.nPx)+'] filament pixel candidates.'
    if count eq 0 then message,'Error: Too stringent local cut to mask.'

    ; run angular friends of friends
    angfof[k,*] = angularFoF(reform(radMasks[k,*]),/discrete)
    nGroupsCand = max(angfof[k,*])
      
    ; determine a overdensity weighted center position of each filament on the shell
    nextGroup = 0
    for i=0,nGroupsCand-1 do begin
      pxInds = where(reform(angfof[k,*]) eq i+1)
      pix2ang_nest, hsv.nSide, pxInds, pxTheta, pxPhi
      pxWts = hsv.value[pxInds,radInd]
      
      ; detect wrap (theta in [0,pi], phi in [0,2pi]) and correct
      if min(pxTheta) lt wrapTol and max(pxTheta) gt !pi-wrapTol then begin
        w = where(pxTheta gt !pi*0.4 and pxTheta lt !pi*0.6,count)
        if count ne 0 then message,'Error: Something odd in theta unwrap.'
        
        ; move lower points up by pi
        wTheta = where(pxTheta lt !pi*0.5,countMoveTheta)
        pxTheta[wTheta] += !pi
      endif
      
      if min(pxPhi) lt wrapTol and max(pxPhi) gt 2*!pi-wrapTol then begin
        w = where(pxPhi gt 2*!pi*0.4 and pxPhi lt 2*!pi*0.6,count)
        if count ne 0 then message,'Error: Something odd in phi unwrap.'
        
        ; move lower points up by 2pi
        wPhi = where(pxPhi lt 2*!pi*0.5,countMovePhi)
        pxPhi[wPhi] += 2*!pi
      endif
      
      ; weighted mean for centroid
      wtSum = total(pxWts)
      cenTheta = total(pxTheta * pxWts) / wtSum
      cenPhi = total(pxPhi * pxWts) / wtSum
      
      ; maximum extent in both axes for effective circular size and circularity measure      
      maxExtent = gcDist([max(pxTheta),max(pxPhi)],[min(pxTheta),min(pxPhi)]) ;[lat,long]

      effCirSize[k,nextGroup] = 2*!pi*(1-cos(maxExtent/2.0))
      cirMeasure[k,nextGroup] = (4*!pi/hsv.nPx*n_elements(pxInds)) / effCirSize[k,nextGroup]
      
      if keyword_set(verbose) then print,radInd,i,maxExtent,cirMeasure[k,nextGroup]
      
      ; move center back into range if necessary
      if cenTheta gt !pi then cenTheta -= !pi
      if cenPhi gt 2*!pi then cenPhi -= 2*!pi
      
      ; convert centroid to pixel to vector
      ang2pix_nest, hsv.nSide, cenTheta, cenPhi, cenPxInd
      pix2vec_nest, hsv.nSide, cenPxInd, cenVec3
      
      ; query circular disc around centroid pixel and recalculate center position (without mask bias)
      query_disc, hsv.nSide, cenVec3, maxExtent/2.0, qdPxInds, /nested
      
      pix2ang_nest, hsv.nSide, qdPxInds, pxTheta, pxPhi
      pxWts = hsv.value[qdPxInds,radInd]
      
      ; detect wrap (theta in [0,pi], phi in [0,2pi]) and correct
      if min(pxTheta) lt wrapTol and max(pxTheta) gt !pi-wrapTol then begin
        w = where(pxTheta gt !pi*0.4 and pxTheta lt !pi*0.6,count)
        if count ne 0 then message,'Error: Something odd in theta unwrap.'
        
        ; move lower points up by pi
        wTheta = where(pxTheta lt !pi*0.5,countMoveTheta)
        pxTheta[wTheta] += !pi
      endif
      
      if min(pxPhi) lt wrapTol and max(pxPhi) gt 2*!pi-wrapTol then begin
        w = where(pxPhi gt 2*!pi*0.4 and pxPhi lt 2*!pi*0.6,count)
        if count ne 0 then message,'Error: Something odd in phi unwrap.'
        
        ; move lower points up by 2pi
        wPhi = where(pxPhi lt 2*!pi*0.5,countMovePhi)
        pxPhi[wPhi] += 2*!pi
      endif
      
      ; weighted mean for centroid
      wtSum = total(pxWts)
      cenTheta2 = total(pxTheta * pxWts) / wtSum
      cenPhi2 = total(pxPhi * pxWts) / wtSum
      
      ; store x,y,z of center position
      ang2pix_nest, hsv.nSide, cenTheta2, cenPhi2, cenPxInd
      pix2vec_nest, hsv.nSIde, cenPxInd, vec3xyz
      
      cenPos[k,nextGroup,*] = reform(float(vec3xyz) * hsv.radFacs[radInd] * hsv.rVir) ; xyz (ckpc) halo centered
      cenPosLL[k,nextGroup,0] = cenTheta2
      cenPosLL[k,nextGroup,1] = cenPhi2

      ; check we satisfied the circularity minimum
      if cirMeasure[k,nextGroup] ge circMeasCutoff and maxExtent lt maxExtentCutoff then begin
        nextGroup += 1
        nGroups[k] += 1
        if nextGroup gt maxNGroups-1 then message,'Error: Exceeded maxNGroups.'
      endif
    endfor ; nGroups
    
    if keyword_set(verbose) then begin
      start_PS,'testmask_r'+str(radInd)+'_fof.eps'
        plotMollweideProj,float(reform(angfof[k,*])),title="",bartitle="",minmax=[0,max(angfof[k,*])]
      end_PS,/deletePS,pngResize=40
    endif
    
    k += 1
  endfor ; radInds
  
  ; if any nGroups are zero then we failed to locate any filaments
  if min(nGroups) eq 0 then begin
    print,'Warning: Failed to find any coherent filaments.'
    return,0
  endif
  
  ; matching arrays
  nFilaments   = nGroups[floor(widthShells/2.0)] ; at rvir
  filMatchInds = intarr(max(nGroups),widthShells) - 1
  filMinDists  = fltarr(max(nGroups),widthShells)
  
  for i=0,nFilaments-1 do filMatchInds[i,floor(widthShells/2.0)] = i ; at rvir filaments match themselves
  
  ; crossmatch sources between the planes: from rvir out
  for radInd=rvirInd[0],maxInd-1 do begin
    k = radInd - rVirInd[0] + floor(widthShells/2.0)

    ; loop over each filament
    for i=0,nFilaments-1 do begin
      curCenPos = cenPosLL[k,i,*]
      
      ; look to next shell outwards and find closest center match
      minDist = 999.9
      minFilInd  = -1
      for j=0,nGroups[k+1]-1 do begin
        otherCenPos = cenPosLL[k+1,j,*]
        otherDist = gcDist(curCenPos,otherCenPos)

        ; new closest?
        if otherDist lt minDist and otherDist lt filDistTol then begin
          if keyword_set(verbose) then print,k,i,j,otherDist
          
          ; has this k+1 candidate already been assigned to a previous filament at level k?
          if total(filMatchInds[*,k+1] eq j) gt 0 then begin
            ; if our current minDist is smaller than the previously used minDist, steal it
            if otherDist lt filMinDists[filMatchInds[i,k],k+1] then begin
              minDist = otherDist
              minFilInd = j
              
              ; set previous holder to unmatched
              w = where(filMatchInds[*,k+1] eq j,count)
              if count ne 1 then message,'Error2'
              filMatchInds[w,k+1] = -1
            endif
          endif else begin
            ; new unique assignment
            minDist = otherDist
            minFilInd = j
          endelse
        endif ; otherDist
      endfor ; nGroups
      
      ; save match
      if keyword_set(verbose) then if minFilInd eq -1 then print,'Warning: Match fail: ',k,k+1,i
      
      if filMatchInds[i,k] ne -1 then filMatchInds[filMatchInds[i,k],k+1] = minFilInd
      if filMatchInds[i,k] eq -1 then message,'Hmm error'
    endfor ; nFilaments
  endfor ; radInds
  
  ; crossmatch sources between the planes: from rvir in
  for radInd=rvirInd[0],minInd+1,-1 do begin
    k = radInd - rVirInd[0] + floor(widthShells/2.0)

    ; loop over each filament
    for i=0,nFilaments-1 do begin
      curCenPos = cenPosLL[k,i,*]
      
      if filMatchInds[i,k] eq -1 then begin
        print,'Warning: Skipping orphan filament at this level.',k,i
        continue
      endif
      
      ; look to next shell outwards and find closest center match
      minDist = 999.9
      minFilInd  = -1
      for j=0,nGroups[k-1]-1 do begin
        otherCenPos = cenPosLL[k-1,j,*]
        otherDist = gcDist(curCenPos,otherCenPos)

        ; new closest?
        if otherDist lt minDist and otherDist lt filDistTol then begin
          if keyword_set(verbose) then print,k,i,j,otherDist
          
          ; has this k-1 candidate already been assigned to a previous filament at level k?
          if total(filMatchInds[*,k-1] eq j) gt 0 then begin
            ; if our current minDist is smaller than the previously used minDist, steal it
            if otherDist lt filMinDists[filMatchInds[i,k],k-1] then begin
              minDist = otherDist
              minFilInd = j
              
              ; set previous holder to unmatched
              w = where(filMatchInds[*,k-1] eq j,count)
              if count ne 1 then message,'Error3'
              filMatchInds[w,k-1] = -1
            endif
          endif else begin
            ; new unique assignment
            minDist = otherDist
            minFilInd = j
          endelse
        endif ; otherDist
      endfor ; nGroups
      
      ; save match
      if keyword_set(verbose) then if minFilInd eq -1 then print,'Warning: Match fail: ',k,k-1,i
      
      if filMatchInds[i,k] ne -1 then begin
        filMatchInds[filMatchInds[i,k],k-1] = minFilInd
        filMinDists[filMatchInds[i,k],k-1] = minDist
      endif
      
    endfor ; nFilaments
  endfor ; radInds
  
  filIntPts    = [] ; best fit intersection point for each filament
  filUnitVecs  = [] ; best fit unit vector along the direction for each filament
  
  ; for each filament, collect the shell intersection points and calculate the ODR fit
  nFound = 0
  for i=0,nFilaments-1 do begin
    xpts = [] & ypts = [] & zpts = []
    
    for k=0,widthShells-1 do begin
      ; if this filament wasn't linked to this shell, don't contribute coordinates to fit
      if filMatchInds[i,k] ne -1 then begin
        xpts = [xpts,cenPos[k,reform(filMatchInds[i,k]),0]]
        ypts = [ypts,cenPos[k,reform(filMatchInds[i,k]),1]]
        zpts = [zpts,cenPos[k,reform(filMatchInds[i,k]),2]]
      endif
    endfor
    
    ; keep filament if we found N or more shell intersections
    if n_elements(xpts) ge minNumDetections then begin
      ; orthogonal distance regression -> best fit line
      odrFit = fitODRPts3D(xpts,ypts,zpts)
      
      filIntPts = [[filIntPts],[odrFit.centroid]]
      filUnitVecs = [[filUnitVecs],[odrFit.line_unitvec]]
      nFound += 1
    endif
  endfor
  
  r = { nFilaments:nFound, filIntPts:filIntPts, filUnitVecs:filUnitVecs, sP:sP, subgroupID:subgroupID, $
        threshOverDensAll:threshOverDensAll, threshOverDensRvir:threshOverDensRvir, $
        widthShells:widthShells, minNumDetections:minNumDetections, wrapTol:wrapTol, $
        maxNGroups:maxNGroups, circMeasCutoff:circMeasCutoff, maxExtentCutoff:maxExtentCutoff, $
        filDistTol:filDistTol, filMatchInds:filMatchInds, filMinDists:filMinDists, $
        cenPosLL:cenPosLL, cirMeasure:cirMeasure}

  return,r

end

; haloRefineFilaments(): refine angularfof based filament directions with a overdensity maximization search

function haloRefineFilaments, sP=sP, hfs=hfs

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  cylRadiusSQ = 40.0*40.0  ; ckpc, pre-square to avoid sqrt
  radRange    = [0.75,1.5] ; r/rvir of points to consider
  uvStepSize  = 0.01
  uvStepRange = 0.2

  ; load gas positions and densities
  pos_gas  = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  dens_gas = loadSnapshotSubset(sP=sP,partType='gas',field='dens')
  u_gas    = loadSnapshotSubset(sP=sP,partType='gas',field='u')
  
  ; spatial subset
  gc = loadGroupCat(sP=sP,/readIDs)
  sgpos  = subgroupPosByMostBoundID(sP=sP)
  cenPos = sgpos[*,hfs.subgroupID]
  rVir   = gc.group_r_crit200[gc.subgroupGrNr[hfs.subgroupID]]
  
  rad_gas = periodicDists(cenPos,pos_gas,sP=sP)
  
  wRadCut = where(rad_gas le radRange[1]*rVir and rad_gas gt radRange[0]*rVir,sCount)

  rad_gas  = rad_gas[wRadCut] / rVir
  pos_gas  = pos_gas[*,wRadCut]
  dens_gas = dens_gas[wRadCut]
  u_gas    = u_gas[wRadCut]

  ; remove substructures
;  ids_gas = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
;  ids_gas = ids_gas[wRadCut]
      
  ; make a list of satellites of this halo and their particle ids
;  nSubs    = gc.groupNSubs[gc.subgroupGrNr[hfs.subgroupID]]
;  firstSub = gc.groupFirstSub[gc.subgroupGrNr[hfs.subgroupID]]
;  satGCids = lindgen(nSubs-1) + firstSub + 1
;  satPIDs = gcPIDList(gc=gc,select='sec',valGCids=satGCids,partType='gas')
      
  ; remove the intersection of (satPIDs,loc_ids) from posval
;  match,satPIDs,ids_gas,sat_ind,ids_ind,count=count,/sort
;  sat_ind = !NULL
      
;  all = bytarr(n_elements(ids_gas))
;  if count gt 0 then all[ids_ind] = 1B
;  wSubSComp = where(all eq 0B, ncomp)
      
;  print,'Substructures cut ['+str(count)+'] of ['+str(n_elements(ids_gas))+'] have left: '+str(ncomp)
   
;  ids_ind = !NULL
;  ids_gas = !NULL
  wRadCut = !NULL
  sgpos = !NULL
  gc = !NULL

  ; take non-substructure cut and make positios relative to halo center
  ;if ncomp gt 0 then begin
  ;  rad_gas  = rad_gas[*,wSubSComp]
  ;  pos_gas  = pos_gas[*,wSubSComp]
  ;  dens_gas = dens_gas[wSubSComp]
  ;  u_gas    = u_gas[wSubSComp]
  ;endif
  
  pos_gas[0,*] -= cenPos[0]
  pos_gas[1,*] -= cenPos[1]
  pos_gas[2,*] -= cenPos[2]
    
  ; arrays
  nUvSteps = fix(uvStepRange / uvStepSize)
  nTrials = long(long(nUvSteps) * nUvSteps * nUvSteps)
  
  newUnitVecs = fltarr(3,hfs.nFilaments)
  
  ; loop over each filament
  for i=0,hfs.nFilaments-1 do begin
    print,i
    ; generate steps on filUnitVec
    trialUnitVecs0 = findgen(nUvSteps)/(nUvSteps-1) * uvStepRange - uvStepRange/2.0 + hfs.filUnitVecs[0,i]
    trialUnitVecs1 = findgen(nUvSteps)/(nUvSteps-1) * uvStepRange - uvStepRange/2.0 + hfs.filUnitVecs[1,i]
    trialUnitVecs2 = findgen(nUvSteps)/(nUvSteps-1) * uvStepRange - uvStepRange/2.0 + hfs.filUnitVecs[2,i]
    
    densFacs = fltarr(nTrials)
    maxInd = 0 & maxInd0 = 0 & maxInd1 = 0 & maxInd2 = 0
    
    ; loop over all trials
    for j=0UL,nTrials-1 do begin
      ind0 = fix(j mod nUvSteps)
      ind1 = fix(j/nUvSteps mod nUvSteps)
      ind2 = fix(j/(nUvSteps*nUvSteps))
      
      ;print,'['+string(j,format='(i3)')+' / '+str(nTrials)+'] '+str(ind0)+' '+str(ind1)+' '+str(ind2)
      
      ; calculate two line points x1,x2
      x1 = hfs.filIntPts[*,i] - 100.0*[trialUnitVecs0[ind0],trialUnitVecs1[ind1],trialUnitVecs2[ind2]]
      x2 = hfs.filIntPts[*,i] + 100.0*[trialUnitVecs0[ind0],trialUnitVecs1[ind1],trialUnitVecs2[ind2]]
      
      ; parametric solution and minimum distance to line (gas)  
      n21  = (x2[0]-x1[0])^2.0 + (x2[1]-x1[1])^2.0 + (x2[2]-x1[2])^2.0 ;n21=transpose(x1-x2)#(x2-x1)
      n10  = reform( (x1[0]-pos_gas[0,*])^2.0 + (x1[1]-pos_gas[1,*])^2.0 + (x1[2]-pos_gas[2,*])^2.0 )
      dotp = reform( (x1[0]-pos_gas[0,*])*(x2[0]-x1[0]) + (x1[1]-pos_gas[1,*])*(x2[1]-x1[1]) + $
                     (x1[2]-pos_gas[2,*])*(x2[2]-x1[2]) )
    
      t_gas = -1.0 * dotp / n21
      d_gas = ( n10 * n21 - dotp^2.0 ) / n21
      ;d_gas = sqrt(d_gas)
      
      ; select gas particles near the filament line
      w = where(d_gas le cylRadiusSQ,count)
      
      ; calculate "goodness of fit" function as mean(density/temperature) to bias towards
      ; low temperature mass overdensities with cylindrical geometry
      if count gt 0 then begin
        densFacs[j] = mean(dens_gas[w]/u_gas[w])
        
        ; if greater than the baseline mark as new best unitvecs
        if densFacs[j] gt densFacs[0] then begin
          maxInd = j & maxInd0 = ind0 & maxInd1 = ind1 & maxInd2 = ind2
        endif
      endif
    endfor
    
    ; if we found a maxInd other than the starting unit vector, add as best
    if maxInd eq 0 then print,'Warning: No better match?'
    
    newUnitVecs[0,i] = trialUnitVecs0[maxInd0]
    newUnitVecs[1,i] = trialUnitVecs1[maxInd1]
    newUnitVecs[2,i] = trialUnitVecs2[maxInd2]
    
  endfor

  hfs.filUnitVecs = newUnitVecs
  return,hfs
end

; makeFilamentProfile(): run filament search and create radial and cross sectional profiles (individual)

function makeFilamentProfile, sP=sP, subgroupID=subgroupID

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config  
  cylMaxRadius = 100.0
  cylRadRange  = [0.1,2.0]
  
  saveFilename = sP.derivPath+'hFil/hFil.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(subgroupID)+'.sav'
 
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif
  
  ; run filament search and refine procedure
  hfs = haloFilamentSearch(sP=sP,subgroupID=subgroupID)
  hfs = haloRefineFilaments(sP=sP,hfs=hfs)

  ; temp: make string
  uvStr = '['
  posStr = '['
  for i=0,hfs.nFilaments-1 do begin
    uvStr = uvStr + '['+string(hfs.filUnitVecs[0,i],format='(f6.3)')+','+$
                        string(hfs.filUnitVecs[1,i],format='(f6.3)')+','+$
                        string(hfs.filUnitVecs[2,i],format='(f6.3)')+']'
    posStr = posStr + '['+string(hfs.filIntPts[0,i],format='(f7.2)')+','+$
                          string(hfs.filIntPts[1,i],format='(f7.2)')+','+$
                          string(hfs.filIntPts[2,i],format='(f7.2)')+']'
                          
    if i ne hfs.nFilaments-1 then uvStr = uvStr + ','                   
    if i ne hfs.nFilaments-1 then posStr = posStr + ','
  endfor
  
  print,uvStr+']'
  print,posStr+']'
  
  ; load halo properties
  sgpos  = subgroupPosByMostBoundID(sP=sP)
  cenPos = sgpos[*,hfs.subgroupID]
  sgpos  = !NULL
  
  gc   = loadGroupCat(sP=sP,/readIDs)
  rVir = gc.group_r_crit200[gc.subgroupGrNr[hfs.subgroupID]]
  
  ; load gas positions and take radial subset
  pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  rad = periodicDists(cenPos,pos,sP=sP)
  wRadCut = where(rad le cylRadRange[1]*1.1*rVir and rad gt cylRadRange[0]*0.5*rVir,sCount)
  rad = rad[wRadCut] / rVir
  
  pos = pos[*,wRadCut]
  
  ; remove substructures
    ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    ids = ids[wRadCut]
    
    ; make a list of satellites of this halo and their particle ids
    satPIDs = gcPIDList(gc=gc,select='sec',partType='gas')
    
    ; remove the intersection of (satPIDs,ids)
    match,satPIDs,ids,sat_ind,ids_ind,count=count,/sort
    sat_ind = !NULL
    
    all = bytarr(n_elements(ids))
    if count gt 0 then all[ids_ind] = 1B
    wSubSComp = where(all eq 0B, ncomp)
    
    print,'Substructures cut ['+str(count)+'] of ['+str(n_elements(ids))+'] have left: '+str(ncomp)

    if ncomp gt 0 then begin
      wRadCut = wRadCut[wSubSComp]
      rad = rad[wSubSComp]
      pos = pos[*,wSubSComp]
    endif
    
    ids_ind = !NULL & ids = !NULL & satPIDs = !NULL & wSubSComp = !NULL
 
  ; load other gas properties and make same subset
  dens  = loadSnapshotSubset(sP=sP,partType='gas',field='dens')
  dens  = dens[wRadCut]
  u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
  u     = u[wRadCut]
  nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
  nelec = nelec[wRadCut]
  temp  = convertUtoTemp(u,nelec,/log)
  nelec = !NULL
  entr  = calcEntropyCGS(u,dens,sP=sP,/log)
  pres  = calcPressureCGS(u,dens,sP=sP,/log)
  u     = !NULL
  
  ; make velocity products
    vel = loadSnapshotSubset(sP=sP,partType='gas',field='vel')
    vel = vel[*,wRadCut]
  
    ; ----- vrad ------
    gVel = gc.subgroupVel[*,subgroupID]
    vel[0,*] = reform(vel[0,*] - gVel[0])
    vel[1,*] = reform(vel[1,*] - gVel[1])
    vel[2,*] = reform(vel[2,*] - gVel[2])
    
    ; make normalized position vector wrt halo center = vec(r) / ||r|| where r from particle to center
    ; means that radvel<0 is inflow and radvel>0 is outflow
    rnorm0 = reform(pos[0,*] - cenPos[0])
    rnorm1 = reform(pos[1,*] - cenPos[1])
    rnorm2 = reform(pos[2,*] - cenPos[2])
    
    correctPeriodicDistVecs, rnorm0, sP=sP
    correctPeriodicDistVecs, rnorm1, sP=sP
    correctPeriodicDistVecs, rnorm2, sP=sP
    
    ; denominator and do divide
    rnorm = sqrt(rnorm0*rnorm0 + rnorm1*rnorm1 + rnorm2*rnorm2)

    rnorm0 /= rnorm
    rnorm1 /= rnorm
    rnorm2 /= rnorm
    
    ; dot(vel,rnorm) gives the magnitude of the projection of vel onto vec(r)
    vrad = reform(vel[0,*]*rnorm0 + vel[1,*]*rnorm1 + vel[2,*]*rnorm2)
  
    ; ----- angm ------
    rnorm0 = reform(cenPos[0] - pos[0,*])
    rnorm1 = reform(cenPos[1] - pos[1,*])
    rnorm2 = reform(cenPos[2] - pos[2,*])
    
    correctPeriodicDistVecs, rnorm0, sP=sP
    correctPeriodicDistVecs, rnorm1, sP=sP
    correctPeriodicDistVecs, rnorm2, sP=sP
    
    ; angular momentum magnitude
    angm = fltarr(3,ncomp)
    angm[0,*] = rnorm1 * vel[2,*] - rnorm2 * vel[1,*]
    angm[1,*] = rnorm2 * vel[0,*] - rnorm0 * vel[2,*]
    angm[2,*] = rnorm0 * vel[1,*] - rnorm1 * vel[0,*]
    
    ; magnitude of specific angular momentum = rvec x vel
    angm = reform(sqrt(angm[0,*]*angm[0,*] + angm[1,*]*angm[1,*] + angm[2,*]*angm[2,*]))  
  
    rnorm0 = !NULL & rnorm1 = !NULL & rnorm2 = !NULL & rnorm = !NULL & vel = !NULL
  
  ; make positions relative to halo center
  pos[0,*] -= cenPos[0]
  pos[1,*] -= cenPos[1]
  pos[2,*] -= cenPos[2]  
  
  ; arrays
  fil_dist = []
  fil_rad  = []
  fil_dens = []
  fil_temp = []
  fil_entr = []
  fil_pres = []
  fil_vrad = []
  fil_angm = []
  fil_num  = []

  ; transform gas positions into a (distance from filament, radius from halo center) coord system
  for i=0,hfs.nFilaments-1 do begin
    ; calculate two line points x1,x2
    x1 = hfs.filIntPts[*,i] - 100.0*[hfs.filUnitVecs[0,i],hfs.filUnitVecs[1,i],hfs.filUnitVecs[2,i]]
    x2 = hfs.filIntPts[*,i] + 100.0*[hfs.filUnitVecs[0,i],hfs.filUnitVecs[1,i],hfs.filUnitVecs[2,i]]
    
    start_PS,'test_xy'+str(i)+'.eps',xs=6,ys=6
      cgplot,pos[0,*],pos[1,*],psym=3,xtitle="x pos",ytitle="y pos"
      cgplot,[x1[0],x2[0]],[x1[1],x2[1]],line=0,color=cgColor('red'),/overplot
    end_PS
    
    ; parametric solution and minimum distance to line (gas)  
    n21  = (x2[0]-x1[0])^2.0 + (x2[1]-x1[1])^2.0 + (x2[2]-x1[2])^2.0 ;n21=transpose(x1-x2)#(x2-x1)
    n10  = reform( (x1[0]-pos[0,*])^2.0 + (x1[1]-pos[1,*])^2.0 + (x1[2]-pos[2,*])^2.0 )
    dotp = reform( (x1[0]-pos[0,*])*(x2[0]-x1[0]) + (x1[1]-pos[1,*])*(x2[1]-x1[1]) + $
                   (x1[2]-pos[2,*])*(x2[2]-x1[2]) )
  
    t_gas = -1.0 * dotp / n21
    d_gas = ( n10 * n21 - dotp^2.0 ) / n21
    d_gas = sqrt(d_gas)

    ; take the subset near the filament
    w = where(d_gas le cylMaxRadius,count)
    print,i,count
    if count eq 0 then message,'Error: No gas found near filament.'
    
    fil_dist = [fil_dist,d_gas[w]]
    fil_rad  = [fil_rad,rad[w]]
    fil_dens = [fil_dens,dens[w]]
    fil_temp = [fil_temp,temp[w]]
    fil_entr = [fil_entr,entr[w]]
    fil_pres = [fil_pres,pres[w]]
    fil_vrad = [fil_vrad,vrad[w]]
    fil_angm = [fil_angm,angm[w]]
    
    fil_num = [fil_num,count]
  endfor
  
  r = { hfs:hfs, fil_dist:fil_dist, fil_rad:fil_rad, fil_dens:fil_dens, fil_temp:fil_temp, $
        fil_entr:fil_entr, fil_pres:fil_pres, fil_vrad:fil_vrad, fil_angm:fil_angm, fil_num:fil_num, $
        cylRadRange:cylRadRange, cylMaxRadius:cylMaxRadius }
  
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))

  return,r
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
  nNGB   = 32  ; neighbor search in CalcTHVal
  padFac = 4.0 ; times r_vir maximum search radius
  
  ; healpix resolution parameter, 8=768, 16~3k, 32~12k, 64~50k, 128~200k, 256~750k, 512~3M
  if ~keyword_set(Nside) then Nside = 64
  
  ; r/r_vir list of shells to compute
  if ~keyword_set(radFacs) then $
    radFacs = [0.01,0.05,0.1,0.25,0.5,0.75,0.9,1.0,1.1,1.25,1.5,1.75,2.0]
    
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
  
  if padFac lt 1.9*max(radFacs) then message,'Error: Suggest increasing padFac.'
  
  ; "radial mass flux" as density * radial velocity (area normalization / r^2 sphere factor omitted)
  ; instead of estimating this value on each particle and then tophat smoothing, we use the
  ; estimates of each avaiable already on each point on the sphere
  if valName eq 'radmassflux' then begin
    if n_elements(subgroupIDs) gt 1 then message,'This makes little sense right now.'
    hsv_dens   = haloShellValue(sP=sP,partType=partType,valName='density',subgroupIDs=subgroupIDs,$
                                Nside=Nside,radFacs=radFacs,cutSubS=cutSubS)
    hsv_radvel = haloShellValue(sP=sP,partType=partType,valName='radvel',subgroupIDs=subgroupIDs,$
                                Nside=Nside,radFacs=radFacs,cutSubS=cutSubS)
                                
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
      
      ; denominator and do divide
      rnorm = sqrt(rnorm0*rnorm0 + rnorm1*rnorm1 + rnorm2*rnorm2)
  
      rnorm0 /= rnorm
      rnorm1 /= rnorm
      rnorm2 /= rnorm
      
      ; dot(vel,rnorm) gives the magnitude of the projection of vel onto vec(r)
      value = loc_vel[0,*]*rnorm0 + loc_vel[1,*]*rnorm1 + loc_vel[2,*]*rnorm2 ; 1xN
      
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
      match,satPIDs,loc_ids,sat_ind,ids_ind,count=count,/sort
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
    
    saveFilename = sP.derivPath+'hShells/hShells.'+partType+'.'+valName+'.'+sP.savPrefix+str(sP.res)+csTag+$
                   '.ns'+str(Nside)+'.'+str(sP.snap)+'.h'+str(subgroupID)+'.'+str(n_elements(radFacs)) + '.sav'
  
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
    posval = !NULL
    value  = !NULL
    sphereXYZ = !NULL
  
  endforeach
  
  if n_elements(subgroupIDs) eq 1 then return,r
end
