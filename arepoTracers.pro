; arepoTracers.pro
; dnelson
; dec 2011
;
; dev for tracer particles

@helper
@cosmoUtil
@cosmoLoad

; plotShocktube():

pro plotShocktube, snap=snap, resize=resize

  print,'Running: snap ['+str(snap)+']'
  
  if (keyword_set(resize)) then deletePS = 1
  if (not keyword_set(resize)) then deletePS = 0

  units = getUnits()
  
  ; config
  tfac = '1.line'
  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath    = workingPath + 'shocktube.tfac'+str(tfac)+'/output/'
  
  colors = ['forest green','crimson']
  
  ; load
  h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
  
  pos_gas    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='x')  
  dens_gas   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='rho')
  u_gas      = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='u')
  
  pos_tracer  = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='x')
  dens_tracer = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='tracer_rho')
  temp_tracer = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='tracer_temp')
  
  ; sort on x-position
  sort_ind_gas    = sort(pos_gas)
  sort_ind_tracer = sort(pos_tracer)
  
  pos_gas  = pos_gas[sort_ind_gas]
  dens_gas = dens_gas[sort_ind_gas]
  u_gas    = u_gas[sort_ind_gas]
  
  pos_tracer  = pos_tracer[sort_ind_tracer]
  dens_tracer = dens_tracer[sort_ind_tracer]
  temp_tracer = temp_tracer[sort_ind_tracer]
  
  ; calculate gas temperature (all units=1)
  meanmolwt = 1.0 * units.mass_proton
  gamma     = 1.4
  
  temp_gas = (gamma-1.0) * u_gas / units.boltzmann * meanmolwt
  
  ; debugging: load y,z
    y_gas    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='y')
    z_gas    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='z')
    
    y_tracer = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='y')
    
    y_gas = y_gas[sort_ind_gas]
    z_gas = z_gas[sort_ind_gas]
    
    y_tracer = y_tracer[sort_ind_tracer]
    
    ; collapse gas onto x-coordinate, check for any differences vs y or z
    num_per = 10
    
    num_gas_new = n_elements(pos_gas) / num_per
    pos_gas_new = fltarr(3,num_gas_new)
    
    for i=0,num_gas_new-1 do begin
      pos_gas_new[0,i] = mean(pos_gas[i*num_per:(i+1)*num_per-1])
      pos_gas_new[1,i] = mean(y_gas[i*num_per:(i+1)*num_per-1])
      pos_gas_new[2,i] = mean(z_gas[i*num_per:(i+1)*num_per-1])
    endfor
  
  ; interpolate tracer density+temp to positions of gas particles
    ;dens_interp = interpol(dens_gas,pos_gas,pos_tracer) ; something fancy
    ;res = dens_interp - dens_tracer
    
    ; nearest neighbor interpolation
    ;inds = round(pos_gas / 20.0 * n_elements(dens_tracer))
    dens_interp = fltarr(n_elements(pos_gas))
    temp_interp = fltarr(n_elements(pos_gas))
    
    for i=0,n_elements(pos_gas)-1 do begin
      ind = where(abs(pos_tracer-pos_gas[i]) eq min(abs(pos_tracer-pos_gas[i])),count)
      if (count eq 0) then begin
        print,'WARNING'
        stop
      endif
      dens_interp[i] = dens_tracer[ind[0]]
      temp_interp[i] = temp_tracer[ind[0]]
    endfor    
  
  ; plot (0) - scatterplot of (x,y) gas and tracer positions
  start_PS, workingPath + 'shocktube_2d_spos_'+string(snap,format='(I3.3)')+'.eps', xs=10.5, ys=2
  
    xrange = [0.0,20.0]
    yrange = [0.0,2.0]
    
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="x Position [Code]",ytitle="y",$
              title="snap="+str(snap)+" time="+string(h.time,format='(f5.2)')+"",$
              charsize=!p.charsize-0.6,position=[0.06,0.2,0.97,0.86]
    
    fsc_plot,pos_tracer,y_tracer,psym=4,symsize=0.2,/overplot,color=fsc_color(colors[1])
    fsc_plot,pos_gas,y_gas,psym=4,symsize=0.2,/overplot,color=fsc_color(colors[0])
  
  end_PS, pngResize=resize, deletePS=deletePS
  
  ; plot (1) - density comparison
  start_PS, workingPath + 'shocktube_2d_dens_'+string(snap,format='(I3.3)')+'.eps'
  
    yrange = [0.05,max(dens_gas)*1.1]
  
    ; plot
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="",ytitle="Density [Code]",$
              title="snap="+str(snap)+" time="+string(h.time,format='(f5.2)')+"",$
              position=[0.25,0.3,0.95,0.9],xtickname=replicate(' ',10)
    
    psym = -4
    line = 0
    thick = 0.5
    
    ; plot tracer
    fsc_plot, pos_tracer,dens_tracer,psym=psym,line=line,thick=thick,$
              symsize=0.2,/overplot,color=fsc_color(colors[1])
    
    ; plot gas
    fsc_plot, pos_gas,dens_gas,psym=psym,line=line,thick=thick,$
              symsize=0.4,/overplot,color=fsc_color(colors[0])

    ; legend
    legend,['gas','tracer'],textcolors=colors,box=0,margin=0.25,/right
    
    ; residual (difference)
    res = abs(dens_interp - dens_gas) / dens_gas
    print,minmax(res)
    
    yrange = minmax(res)*[0.9,1.1]
    ;yrange=[-1e-6,1e-6]
    
    fsc_plot,[0],[0],ytitle="Frac Res",xtitle="x Position [Code]",xrange=xrange,yrange=yrange,/xs,/ys,$
              xticklen=0.05,/nodata,position=[0.25,0.15,0.95,0.3],/noerase,yticks=2
    
    fsc_plot,pos_gas,res,psym=-4,line=line,thick=thick,symsize=0.2,/overplot   
  
  end_PS, pngResize=resize, deletePS=deletePS
  
  ; plot (2) - temperature comparison
  start_PS, workingPath + 'shocktube_2d_temp_'+string(snap,format='(I3.3)')+'.eps'
  
    yrange = [-7.88,-8.07]
  
    ; plot
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="",ytitle="Log ( Temp [Code] )",$
              title="snap="+str(snap)+" time="+string(h.time,format='(f5.2)')+"",$
              position=[0.25,0.3,0.95,0.9],xtickname=replicate(' ',10)
  
    ; plot tracer
    fsc_plot, pos_tracer,alog10(temp_tracer),psym=psym,line=line,thick=thick,$
              symsize=0.2,/overplot,color=fsc_color(colors[1])  
  
    ; plot gas
    fsc_plot, pos_gas,alog10(temp_gas),psym=psym,line=line,thick=thick,$
              symsize=0.4,/overplot,color=fsc_color(colors[0])  
  
    ; legend
    legend,['gas','tracer'],textcolors=colors,box=0,margin=0.25,/right
    
    ; residual (difference)
    res = abs(temp_interp - temp_gas) / abs(temp_gas)
    print,minmax(res)
    yrange = minmax(res)*[0.8,1.2]
    
    fsc_plot,[0],[0],ytitle="Res [in Log]",xtitle="x Position [Code]",xrange=xrange,yrange=yrange,/xs,/ys,$
              xticklen=0.05,/nodata,position=[0.25,0.15,0.95,0.3],/noerase,yticks=2,/ylog
  
    fsc_plot,pos_gas,res,psym=-4,line=line,thick=thick,symsize=0.2,/overplot  
  
  end_PS, pngResize=resize, deletePS=deletePS



  ; load vertical velocity field
  vy_gas    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='vely')
  vy_tracer = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='vely')
  
  vy_gas    = vy_gas[sort_ind_gas]
  vy_tracer = vy_tracer[sort_ind_tracer]
  
  ; plot (3) - v_y
  start_PS, workingPath + 'shocktube_2d_vy_'+string(snap,format='(I3.3)')+'.eps'
  
    yrange = [-1e-2,1e-2]
  
    ; plot
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="x Position [Code]",ytitle="v_y [Code]",$
              title="snap="+str(snap)+" time="+string(h.time,format='(f5.2)')+""
  
    ; plot gas
    fsc_plot, pos_gas,vy_gas,psym=4,line=line,thick=thick,$
              symsize=0.4,/overplot,color=fsc_color(colors[0])  
  
    ; plot tracer
    fsc_plot, pos_tracer,vy_tracer,psym=4,line=line,thick=thick,$
              symsize=0.2,/overplot,color=fsc_color(colors[1])  
  
    ; legend
    legend,['gas','tracer'],textcolors=colors,box=0,margin=0.25,/right
    
  end_PS, pngResize=resize, deletePS=deletePS
  
  ; load vertical velocity field
  vz_gas    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='velz')
  vz_tracer = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='velz')
  
  vz_gas    = vz_gas[sort_ind_gas]
  vz_tracer = vz_tracer[sort_ind_tracer]
  
  ; plot (3) - v_y
  start_PS, workingPath + 'shocktube_2d_vz_'+string(snap,format='(I3.3)')+'.eps'
  
    yrange = [-1e-2,1e-2]
  
    ; plot
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="x Position [Code]",ytitle="v_z [Code]",$
              title="snap="+str(snap)+" time="+string(h.time,format='(f5.2)')+""
  
    ; plot gas
    fsc_plot, pos_gas,vz_gas,psym=4,line=line,thick=thick,$
              symsize=0.4,/overplot,color=fsc_color(colors[0])  
  
    ; plot tracer
    fsc_plot, pos_tracer,vz_tracer,psym=4,line=line,thick=thick,$
              symsize=0.2,/overplot,color=fsc_color(colors[1])  
  
    ; legend
    legend,['gas','tracer'],textcolors=colors,box=0,margin=0.25,/right
    
  end_PS, pngResize=resize, deletePS=deletePS
  stop
end

; plotShocktubeMovie()

pro plotShocktubeMovie

  start = 0
  stop  = 40
  
  resize = 40
  
  for i=start,stop do begin
    plotShocktube,snap=i,resize=resize
  endfor

end

; estimateDensityTophat(): 2D or 3D number density estimator for a given particle type using a square
;                          or cubic filter function (ideally should do with HSML)

function estimateDensityTophat, pos, mass, eval_pos, ndims=ndims, $
                                filterSize=filterSize, nNeighbors=nNeighbors

  if ((keyword_set(filterSize) and keyword_set(nNeighbors)) or not keyword_set(ndims)) then begin
    print,'Error: Only specify one of filterSize or nNeighbors, and ndims.'
    return,0
  endif

  ; config
  ;boxSize = [20.0,2.0,0.0]
  sz_eval = size(eval_pos)
  neval = sz_eval[2]
  
  ;sz_pos = size(pos)
  ;npos   = sz_pos[2]
  
  ; arrays
  eval_dens = fltarr(neval)
  
  ; do density calculation with a constant filter size
  if keyword_set(filterSize) then begin
    invAreaVol = 1.0 / filterSize^ndims
    
    for i=0,neval-1 do begin
   
      ;vecs = fltarr(3,npos)
      
      ; 1d - select neighbors in filter (can use in 2d, e.g. ignore y-coordinate)
      if (ndims eq 1) then $
        dists = sqrt( (pos[0,*]-eval_pos[0,i])^2.0 )
        
      ; 2d - select neighbors in filter
      if (ndims eq 2) then $
        dists = sqrt( (pos[0,*]-eval_pos[0,i])^2.0 + $
                      (pos[1,*]-eval_pos[1,i])^2.0 )
      
      ; 3d - select neighbors in filter
      if (ndims eq 3) then $
        dists = sqrt( (pos[0,*]-eval_pos[0,i])^2.0 + $
                      (pos[1,*]-eval_pos[1,i])^2.0 + $
                      (pos[2,*]-eval_pos[2,i])^2.0 )
      
      ; correct distance vectors for periodic B.C.
      ;for i=0,2 do begin
      ;  w = where(vecs[i,*] gt boxSize[i]*0.5,count)
      ;  if (count ne 0) then $
      ;    vecs[i,w] = vecs[i,w] - boxSize[i]
      ;    
      ;  w = where(vecs[i,*] lt -boxSize[i]*0.5,count)
      ;  if (count ne 0) then $
      ;    vecs[i,w] = boxSize[i] + vecs[i,w]
      ;endfor
      ;
      ;dists = sqrt(vecs[0,*] * vecs[0,*] + vecs[1,*] * vecs[1,*] + vecs[2,*] * vecs[2,*])
      
      ; calculate density = number/area_or_vol
      w = where(dists le filterSize,count)
      
      if (count gt 0) then $
        eval_dens[i] = total(mass[w]) * invAreaVol

    endfor
  endif
  
  ; do density calculation with a constant number of nearest neighbors
  if keyword_set(nNeighbors) then begin
    print,'ERROR implement'
    return,0
  endif
  
  return,eval_dens

end

; compMassDist(): compare mass distribution of tracer particles to the underlying density field of gas

pro compMassDist, snap=snap

  print,'Running: snap ['+str(snap)+']'

  units = getUnits()
  
  ; config
  tfac = 4
  ndims = 1
  
  filterSize = 0.6 ; constant tophat filter size
  nNeighbors = 0   ; constant number of nearest neighbors
  
  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath    = workingPath + 'shocktube.tfac'+str(tfac)+'/output/'
  
  colors = ['forest green','crimson']
  
  ; load
  h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
  
  pos_gas    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pos')  
  dens_gas   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='rho')
  mass_gas =  loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='mass')
  
  pos_tracer  = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='pos')
  mass_tracer = replicate(total(mass_gas) / h.nPartTot[2], h.nPartTot[2])

  ; evaluate tracer density on gas particle positions
  eval_dens = estimateDensityTophat(pos_tracer,mass_tracer,pos_gas,ndims=ndims,$
                                    filterSize=filterSize,nNeighbors=nNeighbors)
  
  ; sort on gas x-position
  sort_ind_gas    = sort(pos_gas[0,*])
  
  x_gas  = reform(pos_gas[0,sort_ind_gas])
  
  dens_gas    = dens_gas[sort_ind_gas]
  dens_tracer = eval_dens[sort_ind_gas]
  
  ; some factor I don't understand offset in density for ndims=1
  if (ndims eq 1) then dens_tracer /= 4.0
  
  ; plot comparison
  start_PS, workingPath + 'shocktube_2d_massdist_'+string(snap,format='(I3.3)')+'.eps'
  
    xrange = [0.0,20.0]
    yrange = [0.05,max(dens_tracer)*1.1]
  
    ; plot
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="",ytitle="Density [Code]",$
              title="snap="+str(snap)+" time="+string(h.time,format='(f5.2)')+$
              " (filterSize="+string(filterSize,format='(f3.1)')+" tFac="+str(tfac)+")",$
              position=[0.25,0.3,0.95,0.9],xtickname=replicate(' ',10)
              
    psym = 4
    line = 0
    thick = 0.5
    
    ; plot gas
    fsc_plot, x_gas,dens_gas,psym=psym,line=line,thick=thick,$
              symsize=0.6,/overplot,color=fsc_color(colors[0])    
    
    ; plot tracer
    fsc_plot, x_gas,dens_tracer,psym=psym,line=line,thick=thick,$
              symsize=0.3,/overplot,color=fsc_color(colors[1])
    
    ; legend
    legend,['gas','tracer'],textcolors=colors,box=0,margin=0.25,/right
    
    ; residual (difference)
    res = dens_tracer - dens_gas
    print,minmax(res)
    
    yrange = minmax(res)
    ;yrange=[-1e-6,1e-6]
    
    fsc_plot,[0],[0],ytitle="Residual",xtitle="Position [Code]",xrange=xrange,yrange=yrange,/xs,/ys,$
              xticklen=0.05,/nodata,position=[0.25,0.15,0.95,0.3],/noerase,yticks=2
    
    fsc_plot,x_gas,res,psym=4,line=line,thick=thick,symsize=0.2,/overplot  
              
  end_PS
  
  stop
  
end
