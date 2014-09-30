; spherePlot.pro
; gas accretion project - visualization/plotting of quantities onto healpix spheres
; dnelson sep.2014

; plotMollweideProj(): plot the all-sky data in mollweide projection using some of the healpix tools
; 
; rot_ang = [lat,long]
; minmax : constrain minmax used for color table
; pos = 'bottom'/'top' : tile x2 in the vertical direction
; pos = 'ul'/'ur'/'ll'/'lr' : tile x2 in both horizontal and vertical directions
; noerase : set for all plots after the first for a compound plot

pro plotMollweideProj, data, rot_ang=rot_ang, minmax=minmax, ctName=ctName, $
                       title=title, bartitle=bartitle, pos=pos, noerase=noerase, bigbar=bigbar, $
                       inds=inds

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
    
    ; 3x2 compact tiling with room for 3 colorbars and GADGET/AREPO title
    if pos eq 'ul_3' then begin
      w_xll = 0.0 & w_xur = 0.5 & w_yll = 0.70 & w_yur = 0.96
      x_title = 0.25 & y_title = 0.96
      pxsize = 0.5
      cbar_units = 1.5
    endif
    if pos eq 'ur_3' then begin
      w_xll = 0.5 & w_xur = 1.0 & w_yll = 0.70 & w_yur = 0.96
      x_title = 0.75 & y_title = 0.96
      pxsize = 0.5
      cbar_units = 1.5
    endif
    if pos eq 'cl_3' then begin
      w_xll = 0.0 & w_xur = 0.5 & w_yll = 0.39 & w_yur = 0.65
      x_title = 0.5 & y_title = 0.96
      pxsize = 0.5
      cbar_units = 1.5
    endif
    if pos eq 'cr_3' then begin
      w_xll = 0.5 & w_xur = 1.0 & w_yll = 0.39 & w_yur = 0.65
      x_title = 0.5 & y_title = 0.96
      pxsize = 0.5
      cbar_units = 1.5
    endif
    if pos eq 'll_3' then begin
      w_xll = 0.0 & w_xur = 0.5 & w_yll = 0.08 & w_yur = 0.34
      x_title = 0.5 & y_title = 0.96
      pxsize = 0.5
      cbar_units = 1.5
    endif
    if pos eq 'lr_3' then begin
      w_xll = 0.5 & w_xur = 1.0 & w_yll = 0.08 & w_yur = 0.34
      x_title = 0.5 & y_title = 0.96
      pxsize = 0.5
      cbar_units = 1.5
    endif
    
    ; 3x3 compact tiling, individual colorbars and titles
    if pos eq '3x3' then begin
      if ~keyword_set(inds) then message,'Error: Specify indices pair with 3x3.'

      w_xll = ( [0.04,0.36,0.68] )[ inds[0] ]
      w_xur = ( [0.34,0.66,0.98] )[ inds[0] ]
      w_yll = ( [0.09,0.42,0.75] )[ inds[1] ]
      w_yur = ( [0.40,0.73,1.03] )[ inds[1] ]
      x_title = w_xur - 0.05
      y_title = w_yur - 0.05
      pxsize = 0.30
      cbar_units = 1.3
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
    
    ; for 3x2 center colorbars for each row
    if pos eq 'ul_3' or pos eq 'cl_3' or pos eq 'll_3' then begin
      cbar_dx = 0.5
      cbar_dy = 1.0/70.0
    endif
    
    ; for 3x3 reduce
    if pos eq '3x3' then begin
      cbar_dx /= 1.5
      cbar_dy /= 1.0
    endif
  endif 
  
  cbar_xll = (w_xll+w_xur)/2.0 - cbar_dx/2.0
  cbar_xur = (w_xll+w_xur)/2.0 + cbar_dx/2.0
  cbar_yur = w_yll - cbar_dy
  cbar_yll = cbar_yur - cbar_dy
  
  cbar_units_y = cbar_yll-cbar_dy*cbar_units
  
  ; for 3x2 center colorbar for each row
  if pos eq 'ul_3' or pos eq 'cl_3' or pos eq 'll_3' then begin
    cbar_xll = 0.5 - cbar_dx/2.0
    cbar_xur = 0.5 + cbar_dx/2.0
  endif
  
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

; plotHaloShellValComp(): compare shell density for four different mass halos at one redshift

pro plotHaloShellValComp

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  partType = 'dm'
  valName  = 'density'
  radInds  = [4,5,7]    ; pre-saved radFacs
  minmax   = [-0.6,2.0] ; log (rho/mean rho)
  rot_ang  = [0,0]      ; [60,-45] ;[lat,long] center in deg (left,up)
  cutSubS  = 1          ; cut satellite substructures out from halo
  takeLog  = 1          ; log quantity
  
  hMassTargets = [12.0,11.5,11.0,10.5]
  
  sP = simParams(res=512,run='tracer',redshift=2.0)

  ; get IDs of mass targets
  subgroupIDs  = massTargetToHaloID(hMassTargets,sP=sP)

  ; plot
  foreach radInd,radInds do begin
    print,radInd
    
    if cutSubS then csTag = '.cutSubS' else csTag = ''
    start_PS, sP.plotPath+'shell_rcomp.' + sP.savPrefix + str(sP.res) + '.' + $
      str(sP.snap) + '.h'+str(radInd)+'_'+partType+csTag+'.eps', xs=6*1.5, ys=6
      
      pos = ['ul_nb','ur_nb','ll_nb','lr_nb']
      xtpos = [0.06,0.56,0.06,0.56]
      ytpos = [0.55,0.55,0.13,0.13]
      
      bartitle = "log ( "+textoidl("\rho / <\rho>")+" )"
        
      for i=0,3 do begin
        ; interpolate onto the shell
        hsd = haloShellValue(sP=sP,partType=partType,valName=valName,subgroupID=subgroupIDs[i],cutSubS=cutSubS)
  
        ; convert densities into ratios to the mean
        healpix_data = hsd.value[*,radInd] / mean(hsd.value[*,radInd])
        healpix_data = reform(healpix_data)
        if takeLog then healpix_data = alog10(healpix_data)
        
        title = sP.run+" "+str(sP.res)+textoidl("^3")+"  z = "+string(sP.redshift,format='(f3.1)')+" "+$
                textoidl("\rho_{"+hsd.partType+"} (r / r_{vir} = "+string(hsd.radFacs[radInd],format='(f4.2)'))+")"
  
        if i eq 0 then $
          plotMollweideProj,healpix_data,rot_ang=rot_ang,title=title,bartitle=bartitle,$
            minmax=minmax,/bigbar,pos=pos[i]
        if i gt 0 then $
          plotMollweideProj,healpix_data,rot_ang=rot_ang,title="",bartitle="",$
            minmax=minmax,pos=pos[i],/noerase
  
        cgText,xtpos[i],ytpos[i],"M = "+string(hMassTargets[i],format='(f4.1)'),$
          /normal,alignment=0.5,color=cgColor('dark gray'),charsize=1.0
      endfor
    end_PS, pngResize=60, /deletePS

  endforeach ;radInds
  
end
  
; plotShellDifferences(): compare temperature and density of inflow/outflow between two runs

pro plotShellDifferences
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  redshift = 2
  res      = 512
  runs     = ['feedback','tracer']
  haloIDs  = [304,314,315,316]        ; [304 314 315 316] z2.304 z2.301 z2.130 z2.64
  radInds  = [5] ;[3,4,5,7]  ; pre-saved radFacs (3=0.25, 4=0.5, 5=0.75, 7=rvir)
  cutSubS  = 1          ; cut satellite substructures out from halo
  
  vRadThresh = 50 ;km/s

  colors = ['red','blue'] ; inflow,outflow
  
  ; load
  foreach run,runs,m do begin
    sP = simParams(res=res,run=run,redshift=redshift)
    sPs = mod_struct( sPs, 'sP'+str(m)+'_'+run, sP )
  endforeach
  
  ; start plot
  start_PS, sPs.(0).plotPath+'shell_diff_'+$
    sPs.(0).savPrefix+'_'+sPs.(1).savPrefix+str(sPs.(0).res)+'_z'+$
    str(redshift)+'_h4.eps', /big

    pos = plot_pos(row=2,col=2,/gap)
    
    foreach haloID,haloIDs,i do begin
    
    ; load
    foreach run,runs,m do begin
      sP = simParams(res=res,run=run,redshift=redshift)
      gcID = getMatchedIDs(simParam=sP,haloID=haloID)
      
      temp = mod_struct( temp, run, $
        haloShellValue(sP=sP,partType='gas',valName='temp',subgroupID=gcID,cutSubS=cutSubS) )
      dens = mod_struct( dens, run, $
        haloShellValue(sP=sP,partType='gas',valName='density',subgroupID=gcID,cutSubS=cutSubS) )
      vrad = mod_struct( vrad, run, $
        haloShellValue(sP=sP,partType='gas',valName='radvel',subgroupID=gcID,cutSubS=cutSubS) )
      rmf  = mod_struct( rmf, run, $
        haloShellValue(sP=sP,partType='gas',valName='radmassflux',subgroupID=gcID,cutSubS=cutSubS) )
      
      sPs = mod_struct( sPs, 'sP'+str(m)+'_'+run, sP )
    endforeach
    
    ;healpix_data = reform(hsv.value[*,radInd] / mean(hsv.value[*,radInd])) 
    
    foreach radInd,radInds,j do begin
      ; select inflow and outflow
      w_in0  = where(vrad.(0).value[*,radInd] lt -vRadThresh,count_in)
      w_out0 = where(vrad.(0).value[*,radInd] gt vRadThresh,count_out)
      w_in1  = where(vrad.(1).value[*,radInd] lt -vRadThresh,count_in)
      w_out1 = where(vrad.(1).value[*,radInd] gt vRadThresh,count_out)
      
      ; histogram temperature difference
      if count_in ge 2 then begin
        temp_loc0 = reform( temp.(0).value[w_in0,radInd] )
        temp_loc1 = reform( temp.(1).value[w_in1,radInd] )
        
        ;plothist, temp_loc0-temp_loc1, /auto, color=cgColor(colors[0]), $
        ;  xtitle="Gas Temp (FB-noFB)",ytitle="N"
          
        zz = temp_loc0-temp_loc1
        w = where(zz ge 0.0,count)
        print,'fraction of inflow hotter in FB run: '+str(float(count)/n_elements(zz)*100)+'%'
      endif
      
      if count_out ge 2 then begin
        temp_loc0 = reform( temp.(0).value[w_out0,radInd] )
        temp_loc1 = reform( temp.(1).value[w_out1,radInd] )
        
        ;plothist, temp_loc0-temp_loc1, /auto, color=cgColor(colors[1]), /overplot
      endif
      
      ; histogram both
      plothist, temp.(0).value[w_in0,radInd], /auto, color=cgColor(colors[0]), $
        xtitle="Gas Temp of Inflow [log K]",ytitle="N", xrange=[3.9,6.4], /xs, pos=pos[i], /noerase
      plothist, temp.(1).value[w_in1,radInd], /auto, color=cgColor(colors[1]), /overplot
      
      ; legend
      cgPlot,[0,0],[0,1000],color=cgColor('light gray'),/overplot
      ;legend,['inflow','outflow'],textcolor=colors,/top,/right
      legend,['FB','noFB'],textcolor=colors,/top,/left
      ;legend,[string(temp.(0).radFacs[radInd],format='(f5.2)')+' rvir',$
      ;        'vradthresh = '+string(vRadThresh,format='(f5.1)')],/top,/left
          
    endforeach ;radInds,j
    
    endforeach ;haloIDs,i
    
  end_PS

  stop

end
  
; plot2Halo3ValComp(): compare two matched halos with three values each (2 column)

pro plot2Halo3ValComp
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  redshift = 2
  
  sPg = simParams(res=512,run='feedback',redshift=float(redshift))
  sPa = simParams(res=512,run='tracer',redshift=float(redshift))
  
  radInd   = 7        ; pre-saved radFacs (3=0.25, 4=0.5, 5=0.75, 7=rvir)
  rot_ang  = [0,45]   ; [lat,long] center in deg (left,up)
  cutSubS  = 1        ; cut satellite substructures out from halo

  haloID = 304 ;z2[304 314 315 316] ;z2.304 z2.301 z2.130 z2.64
  gcID = getMatchedIDs(sPa=sPa,sPg=sPg,haloID=haloID)

  partTypes = ['gas','gas','gas']
  valNames  = ['temp','density','radmassflux']
  ctNames   = ['blue-red2','helix','brewer-browngreen']
  bartitles = ["T_{gas} [_{ }log K_{ }]",$
               "log ( \rho / <\rho> )",$
               "Radial Mass Flux [_{ }log M_{sun} kpc^{-2} Myr^{-1 }]"]

  ; quantity bounds
  ranges = list([4.4,6.2],[-1.0,1.5],[-300,300]) ;[-2.5,2.5]

  ratioToMean = [0,1,0] ; plot value/mean(value) ratio
  plotLog     = [0,1,0] ; plot log(value)
  
  pos = ['ul_3','cl_3','ll_3','ur_3','cr_3','lr_3']

  if cutSubS then csTag = '.cutSubS' else csTag = ''
  
  ; start plot
  start_PS, sPg.plotPath+'shell_valcomp_'+sPg.savPrefix+str(sPg.res)+'_'+sPa.savPrefix+str(sPa.res)+'_z'+$
    str(redshift)+'_h'+str(haloID)+'_r'+str(radInd)+csTag+'.eps', xs=6*1.5, ys=9.5

    for i=0,2 do begin
      ; RUN 1
      ; interpolate onto the shell
      hsv = haloShellValue(sP=sPg,partType=partTypes[i],valName=valNames[i],$
                           subgroupID=gcID.G,cutSubS=cutSubS)

      ; convert values into ratios to the mean
      if ratioToMean[i] then healpix_data = reform(hsv.value[*,radInd] / mean(hsv.value[*,radInd])) $
      else healpix_data = reform(hsv.value[*,radInd])
            
      if plotLog[i] then begin
        w = where(healpix_data gt 0.0,count,comp=wc,ncomp=ncomp)
        if count gt 0 then healpix_data[w]  = alog10(healpix_data[w])
        if ncomp gt 0 then healpix_data[wc] = -alog10(-healpix_data[wc])
      endif
      
      minMaxVal = ranges[i]
      print,'run1: ',minmax(healpix_data),minMaxVal
      
      w = where(healpix_data gt minMaxVal[1]*0.99,count)
      if count gt 0 then healpix_data[w] = minMaxVal[1] * 0.99
      w = where(healpix_data lt minMaxVal[0]*0.99,count)
      if count gt 0 then healpix_data[w] = minMaxVal[0] * 0.99
      
      if i eq 0 then $
        plotMollweideProj,healpix_data,rot_ang=rot_ang,bartitle=bartitles[i],pos=pos[i],$
          ctName=ctNames[i],minmax=ranges[i],title=sPg.simName
      if i gt 0 then $
        plotMollweideProj,healpix_data,rot_ang=rot_ang,title="",bartitle=bartitles[i],pos=pos[i],$
          /noerase,ctName=ctNames[i],minmax=ranges[i]
          
      healpix_data_old = healpix_data
          
      ; RUN 2
      ; interpolate onto the shell
      hsv = haloShellValue(sP=sPa,partType=partTypes[i],valName=valNames[i],$
                           subgroupID=gcID.A,cutSubS=cutSubS)

      ; convert values into ratios to the mean
      if ratioToMean[i] then healpix_data = reform(hsv.value[*,radInd] / mean(hsv.value[*,radInd])) $
      else healpix_data = reform(hsv.value[*,radInd])

      if plotLog[i] then begin
        w = where(healpix_data lt 0.0,count)
        healpix_data = alog10(abs(healpix_data))
        healpix_data[w] = -healpix_data[w]
      endif
      
      minMaxVal = ranges[i]
      print,'run2: ',minmax(healpix_data),minMaxVal
            
      ;w = where(healpix_data gt minMaxVal[1]*0.99,count)
      ;if count gt 0 then healpix_data[w] = minMaxVal[1] * 0.99
      ;w = where(healpix_data lt minMaxVal[0]*0.99,count)
      ;if count gt 0 then healpix_data[w] = minMaxVal[0] * 0.99
      
      if i eq 0 then $
        plotMollweideProj,healpix_data,rot_ang=rot_ang,bartitle="",pos=pos[i+3],$
          /noerase,ctName=ctNames[i],minmax=ranges[i],title=sPa.simName
      if i gt 0 then $
        plotMollweideProj,healpix_data,rot_ang=rot_ang,title="",bartitle="",pos=pos[i+3],$
          /noerase,ctName=ctNames[i],minmax=ranges[i]
    endfor
    
  end_PS, pngResize=60;, /deletePS
  
end

; plot2Halo3ValComp(): compare three matched halos with three values each (3x3)

pro plot3Halo3ValComp

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  res      = 512
  redshift = 2.0
  
  sPs = mod_struct( sPs, 'sP0', simParams(res=res,run='feedback',redshift=redshift) )
  sPs = mod_struct( sPs, 'sP1', simParams(res=res,run='tracer',redshift=redshift) )
  sPs = mod_struct( sPs, 'sP2', simParams(res=res,run='gadget',redshift=redshift) )
  
  radInd    = 7        ; pre-saved radFacs (3=0.25, 4=0.5, 7=rvir)
  rot_ang   = [0,45]   ; [lat,long] center in deg (left,up)
  cutSubS   = 1        ; cut satellite substructures out from halo

  haloID    = 304 ;z2.304 z2.301 z2.130 z2.64
  partTypes = ['gas','gas','gas']
  valNames  = ['density','radvel','radmassflux']
  ctNames   = ['helix','brewer-redblue','brewer-browngreen']
  bartitles = ["log ( \rho / <\rho> )",$
               "V_{rad} [_{ }km/s_{ }]",$
               "log ( dM_{rad }/dt ) [_{ }M_{sun} kpc^{-2} Myr^{-1 }]"]

  ; quantity bounds
  ranges      = list([-0.5,1.5],[-350,350],[-3.0,3.0])
  ratioToMean = [1,0,0] ; plot value/mean(value) ratio
  plotLog     = [1,0,1] ; plot log(value)

  ; find matching subgroup IDs
  subgroupIDs = lonarr( n_tags(sPs) )
  for i=0,n_tags(sPs)-1 do subgroupIDs[i] = getMatchedIDs(simParams=sPs.(i),haloID=haloID)
  
  if cutSubS then csTag = '.cutSubS' else csTag = ''
  
  ; start plot
  plotStr = 'shell_valcomp_'
  for i=0,n_tags(sPs)-1 do plotStr += sPs.(i).savPrefix+str(sPs.(i).res)+'_'
  plotStr += 'z'+string(redshift,format='(f3.1)')+'_h'+str(haloID)+'_r'+str(radInd)+csTag
            
  start_PS, sPs.(0).plotPath + plotStr + '.eps', xs=9*1.5, ys=9.5

    for i=0,2 do begin
      for j=0,2 do begin
      ; shell cutout
      hsv = haloShellValue(sP=sPs.(j),partType=partTypes[i],valName=valNames[i],$
                           subgroupIDs=subgroupIDs[j],cutSubS=cutSubS)

      ; convert values into ratios to the mean
      if ratioToMean[i] then healpix_data = reform(hsv.value[*,radInd] / mean(hsv.value[*,radInd])) $
      else healpix_data = reform(hsv.value[*,radInd])
      
      if plotLog[i] then begin
        w = where(healpix_data gt 0.0,count,comp=wc,ncomp=ncomp)
        if count gt 0 then healpix_data[w]  = alog10(healpix_data[w])
        if ncomp gt 0 then healpix_data[wc] = -alog10(-healpix_data[wc])
      endif

      minMaxVal = ranges[i]
      print,minmax(healpix_data),minMaxVal

      w = where(healpix_data gt minMaxVal[1]*0.99,count)
      if count gt 0 then healpix_data[w] = minMaxVal[1] * 0.99
      w = where(healpix_data lt minMaxVal[0]*0.99,count)
      if count gt 0 then healpix_data[w] = minMaxVal[0] * 0.99

      plotMollweideProj,healpix_data,rot_ang=rot_ang,bartitle=bartitles[i],pos='3x3',$
        /noerase,ctName=ctNames[i],minmax=ranges[i],title="",inds=[i,j]
          
      endfor ;j
    endfor ;i
    
    for j=0,2 do begin
      ypos = ( [0.20,0.53,0.86] )[j]
      cgText,0.026,ypos,sPs.(j).simName,orientation=90,alignment=0.5,color=cgColor('orange'),/normal
    endfor
    
  end_PS, pngResize=60, /deletePS

end

; plotHaloShellValueComp(): compare four different particle fields for one halo at one redshift

pro plotHaloShellValueComp

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  redshift = 2
  
  sP = simParams(res=512,run='gadget',redshift=float(redshift))
  
  subgroupIDs = [981] ;z2.304 g2342 a2132
                      ;z2.301 g2289 a2034
                      ;z2.314 g981 a927
  
  radInds  = [7]        ; pre-saved radFacs (3=0.25, 4=0.5, 7=rvir)
  rot_ang  = [0,45]      ; [lat,long] center in deg (left,up)
  cutSubS  = 1          ; cut satellite substructures out from halo
  
  ; value config
  partTypes = ['gas','gas','gas','gas',$
               'gas','gas','gas','gas'] ;'dm','dm'
  valNames  = ['temp','density','pressure','radmassflux',$
               'radvel','entropy','metallicity','angmom'] ;dm 'density','veldisp'
  ctNames   = ['helix','helix','helix','brewer-redblue',$
               'brewer-redpurple','helix','helix','helix']
  
  bartitles = ["T_{gas} [_{ }log K_{ }]",$
               "log ( \rho / <\rho> )",$
               "log ( P / k_B ) [_{ }K cm^{-3 }]",$
               "Radial Mass Flux [_{ }M_{sun} kpc^{-2} Myr^{-1 }]",$
               "v_{rad} [_{ }km s^{-1 }]",$
               "log ( Entropy ) [_{ }K cm^{2 }]",$
               ;"log ( \rho_{DM} / <\rho_{DM}> )",$
               ;"\sigma_{vel,DM} [_{ } km/s_{ }]"]
               "Metallicity / 0.0127",$
               "log ( Angular Momentum ) [kpc km/s]"]
               
  ; rvir
  ranges = list([4.3,6.5],[-1.0,1.5],[2.0,3.0],[-1.0,1.0],$
                [-400,400],[6.0,9.0],[0.0,1.0],[3.8,5.0]) ;dm[-0.5,1.0],[0.0,200.0]
        
  ; 0.5/0.25 rvir
  ;ranges = list([4.3,7.0],[-1.0,1.5],[3.0,4.5],[-1.5,1.5],$
  ;              [-400,400],[6.0,9.0],[0.0,2.0],[3.5,4.5]) ;dm[-0.5,1.0],[100.0,300.0]
               
  ratioToMean = [0,1,0,0,0,0,0,0] ; plot value/mean(value) ratio
  plotLog     = [0,1,1,0,0,1,0,1] ; plot log(value)
  
  pos = ['ul_nt','ur_nt','ll_nt','lr_nt']
  
  ; plot
  foreach radInd,radInds do begin
    print,radInd
    
    if cutSubS then csTag = '.cutSubS' else csTag = ''
    start_PS, sP.plotPath+'shell_valcomp_'+sP.savPrefix+str(sP.res)+'_z'+$
      str(redshift)+'_h'+str(subgroupIDs[0])+'_r'+str(radInd)+csTag+'.eps', xs=6*1.5, ys=6

      for i=0,3 do begin
        ; interpolate onto the shell
        hsv = haloShellValue(sP=sP,partType=partTypes[i],valName=valNames[i],$
                             subgroupID=subgroupIDs[0],cutSubS=cutSubS)

        ; convert values into ratios to the mean
        if ratioToMean[i] then healpix_data = reform(hsv.value[*,radInd] / mean(hsv.value[*,radInd])) $
        else healpix_data = reform(hsv.value[*,radInd])
        
        if plotLog[i] then begin
          w = where(healpix_data gt 0.0,count,comp=wc,ncomp=ncomp)
          if count gt 0 then healpix_data[w]  = alog10(healpix_data[w])
          if ncomp gt 0 then healpix_data[wc] = -alog10(-healpix_data[wc])
        endif
        
        minMaxVal = ranges[i]

        w = where(healpix_data gt minMaxVal[1]*0.99,count)
        if count gt 0 then healpix_data[w] = minMaxVal[1] * 0.99
        w = where(healpix_data lt minMaxVal[0]*0.99,count)
        if count gt 0 then healpix_data[w] = minMaxVal[0] * 0.99
        
        print,minMaxVal
  
        if i eq 0 then $
          plotMollweideProj,healpix_data,rot_ang=rot_ang,title="",bartitle=bartitles[i],pos=pos[i],$
            ctName=ctNames[i],minmax=ranges[i]
        if i gt 0 then $
          plotMollweideProj,healpix_data,rot_ang=rot_ang,title="",bartitle=bartitles[i],pos=pos[i],$
            /noerase,ctName=ctNames[i],minmax=ranges[i]
      endfor
    end_PS, pngResize=60;, /deletePS
  
    if cutSubS then csTag = '.cutSubS' else csTag = ''
    start_PS, sP.plotPath+'shell_valcomp2_'+sP.savPrefix+str(sP.res)+'_z'+$
      str(redshift)+'_h'+str(subgroupIDs[0])+'_r'+str(radInd)+csTag+'.eps', xs=6*1.5, ys=6
      
      for i=4,7 do begin
        ; interpolate onto the shell
        hsv = haloShellValue(sP=sP,partType=partTypes[i],valName=valNames[i],$
                             subgroupID=subgroupIDs[0],cutSubS=cutSubS)

        ; convert values into ratios to the mean
        if ratioToMean[i] then healpix_data = reform(hsv.value[*,radInd] / mean(hsv.value[*,radInd])) $
        else healpix_data = reform(hsv.value[*,radInd])

        if plotLog[i] then begin
          w = where(healpix_data gt 0.0,count,comp=wc,ncomp=ncomp)
          if count gt 0 then healpix_data[w]  = alog10(healpix_data[w])
          if ncomp gt 0 then healpix_data[wc] = -alog10(-healpix_data[wc])
        endif

        if valNames[i] eq 'metallicity' then healpix_data /= 0.0127 ; display rescaling

        minMaxVal = ranges[i]

        w = where(healpix_data gt minMaxVal[1]*0.99,count)
        if count gt 0 then healpix_data[w] = minMaxVal[1] * 0.99
        w = where(healpix_data lt minMaxVal[0]*0.99,count)
        if count gt 0 then healpix_data[w] = minMaxVal[0] * 0.99
        
        print,minMaxVal
        
        if i eq 0 then $
          plotMollweideProj,healpix_data,rot_ang=rot_ang,title="",bartitle=bartitles[i],pos=pos[i-4],$
            ctName=ctNames[i],minmax=ranges[i]
        if i gt 0 then $
          plotMollweideProj,healpix_data,rot_ang=rot_ang,title="",bartitle=bartitles[i],pos=pos[i-4],$
            /noerase,ctName=ctNames[i],minmax=ranges[i]
  
      endfor
    end_PS, pngResize=60;, /deletePS
  
  endforeach ;radInds

  stop
end
