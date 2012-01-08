; arepoTracers.pro
; dnelson
; jan 2012
;
; dev for tracer particles

@helper
@cosmoUtil
@cosmoLoad

function interpNN, pos_src, pos_dest, val_src

    ; something fancy (linear or cubic?)
    ;dens_interp = interpol(dens_gas,pos_gas,pos_tracer)
    ;res = dens_interp - dens_tracer
    
    ; nearest neighbor interpolation by indices
    ;inds = round(pos_gas / 20.0 * n_elements(dens_tracer))
    
    ; nearest neighbor interpolation (manual)
    val_interp = fltarr(n_elements(pos_src))
    
    for i=0,n_elements(pos_src)-1 do begin
      ind = where(abs(pos_dest-pos_src[i]) eq min(abs(pos_dest-pos_src[i])),count)
      
      if (count eq 0) then begin
        print,'WARNING'
        stop
      endif
      val_interp[i] = val_src[ind[0]]
    endfor
    
  return, val_interp
end

; plotShocktube():

pro plotShocktube, snap=snap

  print,'Running: snap ['+str(snap)+']'

  units = getUnits()
  
  ; config
  tfac = '1'
  workingPath = '/n/home07/dnelson/dev.tracer/'
  ;snapPath    = workingPath + 'st2d.tfac'+str(tfac)+'.1storder/output/'
  snapPath = workingPath + 'st2d.autogen/output/'
  plotBase    = 'st2d_autogen_'
  
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
  dens_interp = interpNN(pos_gas,pos_tracer,dens_tracer)
  temp_interp = interpNN(pos_gas,pos_tracer,temp_tracer)
  
  ; plot (0) - scatterplot of (x,y) gas and tracer positions
  start_PS, workingPath + plotBase + 'spos_'+string(snap,format='(I3.3)')+'.tfac='+str(tfac)+'.eps',$
            xs=10.5, ys=2
  
    xrange = [0.0,20.0]
    yrange = [0.0,2.0]
    
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="Position [Code]",ytitle="y",$
              title="snap="+str(snap)+" time="+string(h.time,format='(f5.2)')+"",$
              charsize=!p.charsize-0.6,position=[0.06,0.2,0.97,0.86]
    
    fsc_plot,pos_gas,y_gas,psym=4,symsize=0.2,/overplot,color=fsc_color(colors[0])
    fsc_plot,pos_tracer,y_tracer,psym=4,symsize=0.2,/overplot,color=fsc_color(colors[1])
  
  end_PS
  
  ; plot (1) - density comparison
  start_PS, workingPath + plotBase + 'dens_'+string(snap,format='(I3.3)')+'.tfac='+str(tfac)+'.eps'
  
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
    w = where(res ne 0.0)
    print,'residual minmax: ',minmax(res[w])
    
    ;yrange = alog10(minmax(res[w]))
    yrange = [-7.5,-0.5]
    
    fsc_plot,[0],[0],ytitle="log(Frac Error)",xtitle="Position [Code]",xrange=xrange,yrange=yrange,/xs,/ys,$
              xticklen=0.05,/nodata,position=[0.25,0.15,0.95,0.3],/noerase,yticks=2
    
    fsc_plot,xrange,[-2.0,-2.0],line=0,color=fsc_color('gray'),thick=!p.thick-1.0,/overplot
    fsc_plot,pos_gas,alog10(res[w]),psym=-4,line=line,thick=thick,symsize=0.2,/overplot   
  
  end_PS
  
  ; plot (2) - temperature comparison
  start_PS, workingPath + plotBase + 'temp_'+string(snap,format='(I3.3)')+'.tfac='+str(tfac)+'.eps'
  
    yrange = [-7.88,-8.07]
  
    ; plot
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="",ytitle="log(Temp) [Code]",$
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
    w = where(res ne 0.0)
    
    yrange = alog10(minmax(res[w]))

    fsc_plot,[0],[0],ytitle="log(Frac Error)",xtitle="Position [Code]",xrange=xrange,yrange=yrange,/xs,/ys,$
              xticklen=0.05,/nodata,position=[0.25,0.15,0.95,0.3],/noerase,yticks=2
  
    fsc_plot,xrange,[-2.0,-2.0],line=0,color=fsc_color('gray'),thick=!p.thick-1.0,/overplot
    fsc_plot,pos_gas,alog10(res[w]),psym=-4,line=line,thick=thick,symsize=0.2,/overplot  
  
  end_PS

  ; load vertical velocity field
  vy_gas    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='vely')
  vy_tracer = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='vely')
  
  vy_gas    = vy_gas[sort_ind_gas]
  vy_tracer = vy_tracer[sort_ind_tracer]
  
  ; plot (3) - v_y
  start_PS, workingPath + plotBase + 'vy_'+string(snap,format='(I3.3)')+'.tfac='+str(tfac)+'.eps'
  
    yrange = [-1e-4,1e-4]
  
    ; plot
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="Position [Code]",ytitle="v_y [Code]",$
              title="snap="+str(snap)+" time="+string(h.time,format='(f5.2)')+""
  
    ; plot gas
    fsc_plot, pos_gas,vy_gas,psym=4,line=line,thick=thick,$
              symsize=0.4,/overplot,color=fsc_color(colors[0])  
  
    ; plot tracer
    fsc_plot, pos_tracer,vy_tracer,psym=4,line=line,thick=thick,$
              symsize=0.2,/overplot,color=fsc_color(colors[1])  
  
    ; legend
    legend,['gas','tracer'],textcolors=colors,box=0,margin=0.25,/right
    
  end_PS
  
  ; load vertical velocity field
  vz_gas    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='velz')
  vz_tracer = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='velz')
  
  vz_gas    = vz_gas[sort_ind_gas]
  vz_tracer = vz_tracer[sort_ind_tracer]
  
  ; plot (3) - v_z
  start_PS, workingPath + plotBase + 'vz_'+string(snap,format='(I3.3)')+'.tfac='+str(tfac)+'.eps'
  
    yrange = [-1e-4,1e-4]
  
    ; plot
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="Position [Code]",ytitle="v_z [Code]",$
              title="snap="+str(snap)+" time="+string(h.time,format='(f5.2)')+""
  
    ; plot gas
    fsc_plot, pos_gas,vz_gas,psym=4,line=line,thick=thick,$
              symsize=0.4,/overplot,color=fsc_color(colors[0])  
  
    ; plot tracer
    fsc_plot, pos_tracer,vz_tracer,psym=4,line=line,thick=thick,$
              symsize=0.2,/overplot,color=fsc_color(colors[1])  
  
    ; legend
    legend,['gas','tracer'],textcolors=colors,box=0,margin=0.25,/right
    
  end_PS
  
end

; estimateDensityTophat(): spatial density estimator for an input position tuple array and
;                          mass array for some particle type
;   do density calculation by calculating smoothing lengths for all the particles
;   use CalcHSML external C-routine for the tree, neighbor searchings, and smoothing lengths

function estimateDensityTophat, pos, mass, ndims=ndims, nNGB=nNGB

  if (not keyword_set(nNGB) or not keyword_set(ndims)) then begin
    print,'Error: Must specify both nNGB and ndims.'
    return,0
  endif

  ; config
  sz_pos = size(pos)
  npos   = sz_pos[2]
  
  ; arrays
  eval_dens = fltarr(npos)

  ; prepare inputs
  NumPart = long(npos)
  ;Pos     = pos
  ;Mass    = mass
  
  DesNumNgb    = long(nNGB) ; number of neighbors to use
  DesNumNgbDev = long(0)
  BoxSize      = 0.0
  HsmlGuess    = float(1.0)
  Softening    = float(1.0)
  
  hsml_out = fltarr(NumPart)
  
  ; call CalcHSML
  libName = '/n/home07/dnelson/idl/CalcHSML/CalcHSML_'+str(ndims)+'D.so'
  ret = Call_External(libName, 'CalcHSML', $
                      NumPart,Pos,Mass,DesNumNgb,DesNumNgbDev,BoxSize,HsmlGuess,Softening,hsml_out, $
                      /CDECL)

  ; estimate densities on eval_pos using hsml
  if (ndims eq 1) then $
    eval_dens = DesNumNgb / hsml_out
    
  if (ndims eq 2) then $
    eval_dens = DesNumNgb / (!pi * hsml_out^2.0)
    
  if (ndims eq 3) then $
    eval_dens = DesNumNgb / (4.0*!pi/3.0 * hsml_out^3.0)
  
  ; add in mass
  eval_dens *= mass[0]
  
  return,eval_dens

end

; compMassDist(): compare mass distribution of tracer particles to the underlying density field of gas

pro compMassDist, snap=snap

  print,'Running: snap ['+str(snap)+']'

  units = getUnits()
  
  ; config
  tfac = '16'
  ndims = 2
  nNGB  = 128 ; number of neighbors to use with external CalcHSML routine
             ; which finds smoothing lengths for tophat density estimate
  
  xrange = [0.0,20.0]  
  
  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath    = workingPath + 'st2d.tfac'+str(tfac)+'.1storder/output/'
  ;snapPath = workingPath + 'st2d.autogen/output/'
  plotBase    = 'st2d_1st_'
  
  colors = ['forest green','crimson']
  
  ; load
  h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
  
  pos_gas    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pos')  
  dens_gas   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='rho')
  mass_gas =  loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='mass')
  
  pos_tracer  = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='pos')
  mass_tracer = replicate(total(mass_gas) / h.nPartTot[2], h.nPartTot[2])

  ; evaluate tracer density
  eval_dens = estimateDensityTophat(pos_tracer,mass_tracer,ndims=ndims,nNGB=nNGB)
  
  ; sort on gas x-position
  sort_ind_gas    = sort(pos_gas[0,*])
  sort_ind_tracer = sort(pos_tracer[0,*])
  
  x_gas    = reform(pos_gas[0,sort_ind_gas])
  x_tracer = reform(pos_tracer[0,sort_ind_tracer])
  y_tracer = reform(pos_tracer[1,sort_ind_tracer])
  z_tracer = reform(pos_tracer[2,sort_ind_tracer])
  
  dens_gas    = dens_gas[sort_ind_gas]
  dens_tracer = eval_dens[sort_ind_tracer]

  ; tracer selection to avoid periodic edges (CalcHSML not doing periodic correctly for LONG_X)
  wBounds = [0.6,1.4]
  if (ndims eq 2) then $
    wTracer = where(y_tracer ge wBounds[0] and y_tracer le wBounds[1],count)
  if (ndims eq 3) then $
    wTracer = where(y_tracer ge wBounds[0] and y_tracer le wBounds[1] and $
                    z_tracer ge wBounds[0] and z_tracer le wBounds[1],count)
  
  ; interpolate tracer density to positions of gas particles
  dens_interp = interpNN(x_gas,x_tracer[wTracer],dens_tracer[wTracer])  

  ; plot comparison
  start_PS, workingPath + plotBase + 'massdist_'+string(snap,format='(I3.3)')+'.tfac='+str(tfac)+'.ngb='+$
            str(nNGB)+'.eps'

    yrange = [0.05,max(dens_tracer)*1.1]
  
    ; plot
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="",ytitle="Density [Code]",$
              title="snap="+str(snap)+" time="+string(h.time,format='(f5.2)')+$
              " (nNGB="+str(nNGB)+" tFac="+str(tfac)+")",$
              position=[0.25,0.4,0.95,0.9],xtickname=replicate(' ',10)
              
    psym = 4
    line = 0
    thick = 0.5
    
    ; plot gas
    fsc_plot, x_gas,dens_gas,psym=-psym,line=line,thick=thick,$
              symsize=0.6,/overplot,color=fsc_color(colors[0])    
    
    ; plot tracer (avoiding edges with wTracer)   
    fsc_plot, x_tracer[wTracer],dens_tracer[wTracer],psym=psym,line=line,thick=thick,$
              symsize=0.2,/overplot,color=fsc_color(colors[1])
    
    ; legend
    legend,['gas','tracer'],textcolors=colors,box=0,margin=0.25,/right
    
    ; residual (difference)
    res = abs(dens_interp - dens_gas) / dens_gas
    w = where(res ne 0.0 and x_gas gt 1.0 and x_gas lt 19.0)
    print,'residual minmax: ',minmax(res[w])
    
    ;yrange = alog10(minmax(res[w]))
    yrange = [-5.0,-0.5]
    
    fsc_plot,[0],[0],ytitle="log(Frac Error)",xtitle="Position [Code]",xrange=xrange,yrange=yrange,/xs,/ys,$
              xticklen=0.05,/nodata,position=[0.25,0.15,0.95,0.4],/noerase,yticks=3
    
    fsc_plot,xrange,[-2.0,-2.0],line=0,color=fsc_color('gray'),thick=!p.thick-1.0,/overplot
    
    fsc_plot,x_gas[w],alog10(res[w]),psym=-4,line=line,thick=thick,symsize=0.2,/overplot  
              
  end_PS, pngResize=40
  
  stop
  
end

; compMassDistTwoRuns():  compare mass distribution of tracer particles to the underlying density field of gas

pro compMassDistTwoRuns, snap=snap

  print,'Running: snap ['+str(snap)+']'

  units = getUnits()
  
  ; config
  ndims = 2
  nNGB  = 32 ; number of neighbors to use with external CalcHSML routine
             ; which finds smoothing lengths for tophat density estimate
  
  xrange = [0.0,20.0]  
  
  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath1   = workingPath + 'st2d.tfac1.1storder/output/'
  snapPath2   = workingPath + 'st2d.tfac1.2ndorder/output/'
  
  ;legend = ['original','16x smaller timestep']
  ;legend = ['global timestep','variable timestep']
  legend = ['first order','second order']
  
  colors = ['forest green','crimson']
  
  ; load (1)
  h = loadSnapshotHeader(snapPath1,snapNum=snap,/verbose)
  
  pos_gas1    = loadSnapshotSubset(snapPath1,snapNum=snap,partType='gas',field='pos')  
  dens_gas1   = loadSnapshotSubset(snapPath1,snapNum=snap,partType='gas',field='rho')
  mass_gas1 =  loadSnapshotSubset(snapPath1,snapNum=snap,partType='gas',field='mass')
  
  pos_tracer1  = loadSnapshotSubset(snapPath1,snapNum=snap,partType='tracer',field='pos')
  mass_tracer1 = replicate(total(mass_gas1) / h.nPartTot[2], h.nPartTot[2])

  ; evaluate tracer density
  eval_dens1 = estimateDensityTophat(pos_tracer1,mass_tracer1,ndims=ndims,nNGB=nNGB)
  
  ; sort on gas x-position
  sort_ind_gas    = sort(pos_gas1[0,*])
  sort_ind_tracer = sort(pos_tracer1[0,*])
  
  x_gas1    = reform(pos_gas1[0,sort_ind_gas])
  x_tracer1 = reform(pos_tracer1[0,sort_ind_tracer])
  y_tracer1 = reform(pos_tracer1[1,sort_ind_tracer])
  
  dens_gas1    = dens_gas1[sort_ind_gas]
  dens_tracer1 = eval_dens1[sort_ind_tracer]
  
  ; load (2)
  h = loadSnapshotHeader(snapPath2,snapNum=snap,/verbose)
  
  pos_gas2  = loadSnapshotSubset(snapPath2,snapNum=snap,partType='gas',field='pos')  
  dens_gas2 = loadSnapshotSubset(snapPath2,snapNum=snap,partType='gas',field='rho')
  mass_gas2 = loadSnapshotSubset(snapPath2,snapNum=snap,partType='gas',field='mass')
  
  pos_tracer2  = loadSnapshotSubset(snapPath2,snapNum=snap,partType='tracer',field='pos')
  mass_tracer2 = replicate(total(mass_gas2) / h.nPartTot[2], h.nPartTot[2])

  ; evaluate tracer density
  eval_dens2 = estimateDensityTophat(pos_tracer2,mass_tracer2,ndims=ndims,nNGB=nNGB)

  ; sort on gas x-position
  sort_ind_gas    = sort(pos_gas2[0,*])
  sort_ind_tracer = sort(pos_tracer2[0,*])
  
  x_gas2    = reform(pos_gas2[0,sort_ind_gas])
  x_tracer2 = reform(pos_tracer2[0,sort_ind_tracer])
  y_tracer2 = reform(pos_tracer2[1,sort_ind_tracer])
  
  dens_gas2    = dens_gas2[sort_ind_gas]
  dens_tracer2 = eval_dens2[sort_ind_tracer]
  
  ; gas and tracer density differences
  gas_dens_diff    = abs(dens_gas1 - dens_gas2) / dens_gas1
  tracer_dens_diff = abs(dens_tracer1 - dens_tracer2) / dens_tracer1
  
  ; tracer selection to avoid vertical edges (CalcHSML not doing periodic correctly)
  wTracer1 = where(y_tracer1 ge 0.6 and y_tracer1 le 1.4,count1)
  wTracer2 = where(y_tracer2 ge 0.6 and y_tracer2 le 1.4,count2)
  wTracer  = where(y_tracer1 ge 0.6 and y_tracer2 ge 0.6 and $
                   y_tracer1 lt 1.4 and y_tracer2 lt 1.4 and tracer_dens_diff ne 0.0,count3)

  ; plot comparison (gas)
  start_PS, workingPath + 'st2d_twosims_gas_'+string(snap,format='(I3.3)')+'.'+str(nNGB)+'.eps'
  
    psym = 4
    line = 0
    thick = 0.5
  
    ; main plot
    yrange = [min(dens_gas1)*0.9,max(dens_gas1)*1.1]
    
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="",ytitle="Gas Rho [Code]",$
              title="snap="+str(snap)+" time="+string(h.time,format='(f5.2)')+"",$
              position=[0.25,0.4,0.95,0.9],xtickname=replicate(' ',10)

    fsc_plot, x_gas1,dens_gas1,psym=-psym,line=line,thick=thick,$
              symsize=0.4,/overplot,color=fsc_color(colors[0])   
    fsc_plot, x_gas2,dens_gas2,psym=-psym,line=line,thick=thick,$
              symsize=0.4,/overplot,color=fsc_color(colors[1])    
    
    legend,legend,textcolors=colors,box=0,margin=0.25,/right
    
    ; residuals
    w = where(gas_dens_diff ne 0.0)
    ;yrange = alog10([min(gas_dens_diff[w])*0.9,max(gas_dens_diff[w])*1.1])
    yrange = [-7.0,-1.0]
    
    fsc_plot,[0],[0],ytitle="log(Frac Error)",xtitle="Position [Code]",xrange=xrange,yrange=yrange,/xs,/ys,$
              xticklen=0.05,/nodata,position=[0.25,0.15,0.95,0.4],/noerase,yticks=3
    
    fsc_plot,xrange,[-2.0,-2.0],line=0,color=fsc_color('gray'),thick=!p.thick-1.0,/overplot
    fsc_plot,x_gas1[w],alog10(gas_dens_diff[w]),psym=-4,line=line,thick=thick,symsize=0.2,/overplot  
    
  end_PS, pngResize=40
  
  ; plot comparison (tracer)
  start_PS, workingPath + 'st2d_twosims_tracer_'+string(snap,format='(I3.3)')+'.'+str(nNGB)+'.eps'
  
    ; main plot
    yrange = [min(dens_tracer1)*0.9,max(dens_tracer1)*1.1]
    
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="",ytitle="Tracer Density",$
              title="snap="+str(snap)+" time="+string(h.time,format='(f5.2)')+" (nNGB="+str(nNGB)+")",$
              position=[0.25,0.4,0.95,0.9],xtickname=replicate(' ',10)
    
    fsc_plot, x_tracer1[wTracer1],dens_tracer1[wTracer1],psym=psym,line=line,thick=thick,$
              symsize=0.2,/overplot,color=fsc_color(colors[0])
              
    fsc_plot, x_tracer2[wTracer2],dens_tracer2[wTracer2],psym=psym,line=line,thick=thick,$
              symsize=0.2,/overplot,color=fsc_color(colors[1])
              
    legend,legend,textcolors=colors,box=0,margin=0.25,/right
    
    ; residuals
    ;yrange = alog10([min(tracer_dens_diff[wTracer])*0.9,max(tracer_dens_diff[wTracer])*1.1])
    yrange = [-5.0,-1.0]
    
    fsc_plot,[0],[0],ytitle="log(Frac Error)",xtitle="Position [Code]",xrange=xrange,yrange=yrange,/xs,/ys,$
              xticklen=0.05,/nodata,position=[0.25,0.15,0.95,0.4],/noerase,yticks=4
              
    fsc_plot,xrange,[-2.0,-2.0],line=0,color=fsc_color('gray'),thick=!p.thick-1.0,/overplot    
    fsc_plot,x_tracer1[wTracer],alog10(tracer_dens_diff[wTracer]),$
             psym=-4,line=line,thick=thick,symsize=0.2,/overplot  
    
  end_PS, pngResize=40

end

; compGasExact(): compare gas properties between runs (which should be identical)

pro compGasExact, snap=snap

  print,'Running: snap ['+str(snap)+']'
  
  workingPath = '/n/home07/dnelson/dev.tracer/'
  
  ;snapPath1   = workingPath + 'st2d.tfac1.2ndorder/output/'

  ;snapPath1   = workingPath + 'st2d.tfac1.1storder/output.timebins/'
  ;snapPath1   = workingPath + 'st2d.tfac1.1storder/output.direct/'
  snapPath1   = workingPath + 'st2d.tfac1.1storder/output/'
  
  snapPath2   = workingPath + 'st2d.gasonly/output/'
  
  ; load (1)
  x_gas1    = loadSnapshotSubset(snapPath1,snapNum=snap,partType='gas',field='x')
  dens_gas1 = loadSnapshotSubset(snapPath1,snapNum=snap,partType='gas',field='rho')
  ids_gas1  = loadSnapshotSubset(snapPath1,snapNum=snap,partType='gas',field='ids')
  
  ; load (2)
  x_gas2    = loadSnapshotSubset(snapPath2,snapNum=snap,partType='gas',field='x')
  dens_gas2 = loadSnapshotSubset(snapPath2,snapNum=snap,partType='gas',field='rho')
  ids_gas2  = loadSnapshotSubset(snapPath2,snapNum=snap,partType='gas',field='ids')
  
  ; sort on id
  sort1 = sort(ids_gas1)
  sort2 = sort(ids_gas2)
  
  dens_gas_sorted1 = dens_gas1[sort1]
  dens_gas_sorted2 = dens_gas2[sort2]
  
  diff = dens_gas_sorted1 - dens_gas_sorted2
  
  ; where is it different
  w = where(diff gt 0.001,count)
  
  x_gas_sorted1 = x_gas1[sort1]
  x_gas_sorted2 = x_gas2[sort2]
  
  print,'gas dens total diff: ',total(diff),count
  if (count gt 0) then $
    print,'positions: ',x_gas_sorted1[w]
  
  ; quick plot
  start_PS,'test2.eps'
  
    fsc_plot,x_gas1,dens_gas1,psym=4,symsize=0.6,color=fsc_color('crimson'),xrange=[0.0,20.0]
    fsc_plot,x_gas2,dens_gas2,psym=4,symsize=0.6,color=fsc_color('slate blue'),/overplot
  
  end_PS
  
  ;start_PS,'test3.eps'
  ;  plothist,diff,/auto
  ;end_PS
  stop
  
end

; compCurMaxTemp(): compare profiles of current to past maximum tracer temperature

pro compCurMaxTemp, snap=snap

  print,'Running: snap ['+str(snap)+']'

  units = getUnits()
  
  ; config
  ndims = 2
  
  xrange = [0.0,20.0]  
  
  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath    = workingPath + 'test.maxtemp/output/'
  plotBase    = 'testmaxt_'
  
  colors = ['forest green','crimson']
  
  ; load
  h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
  
  pos_tracer     = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='pos')
  curtemp_tracer = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='tracer_temp')
  maxtemp_tracer = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='tracer_maxtemp')

  ; sort on tracer x-position
  sort_ind_tracer = sort(pos_tracer[0,*])
  
  x_tracer = reform(pos_tracer[0,sort_ind_tracer])
  y_tracer = reform(pos_tracer[1,sort_ind_tracer])
  z_tracer = reform(pos_tracer[2,sort_ind_tracer])
  
  curtemp_tracer = curtemp_tracer[sort_ind_tracer]
  maxtemp_tracer = maxtemp_tracer[sort_ind_tracer]
  
  ; tracer selection to avoid periodic edges (CalcHSML not doing periodic correctly for LONG_X)
  wBounds = [0.4,1.6]
  if (ndims eq 2) then $
    wTracer = where(y_tracer ge wBounds[0] and y_tracer le wBounds[1],count)
  if (ndims eq 3) then $
    wTracer = where(y_tracer ge wBounds[0] and y_tracer le wBounds[1] and $
                    z_tracer ge wBounds[0] and z_tracer le wBounds[1],count)
  
  ; plot comparison
  start_PS, workingPath + plotBase +string(snap,format='(I3.3)')+'.eps'

    yrange = alog10(minmax(curtemp_tracer))*[1.005,0.995]
  
    ; plot
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="",ytitle="log( Tracer Temp )",$
              title="snap="+str(snap)+" time="+string(h.time,format='(f5.2)'),$
              position=[0.25,0.4,0.95,0.9],xtickname=replicate(' ',10)
              
    psym = 4
    line = 0
    thick = 0.5
    
    ; plot tracer (avoiding edges with wTracer)   
    fsc_plot, x_tracer[wTracer],alog10(curtemp_tracer[wTracer]),psym=psym,line=line,thick=thick,$
              symsize=0.2,/overplot,color=fsc_color(colors[0])
    fsc_plot, x_tracer[wTracer],alog10(maxtemp_tracer[wTracer]),psym=psym,line=line,thick=thick,$
              symsize=0.2,/overplot,color=fsc_color(colors[1])         
    
    ; legend
    legend,['current','past max'],textcolors=colors,box=0,margin=0.25,/right
    
    ; residual (difference)
    res = abs(curtemp_tracer - maxtemp_tracer) / curtemp_tracer
    w = where(res ne 0.0)
    
    yrange = alog10(minmax(res[w]))
    if (snap eq 0) then yrange = [-3.0,-1.0] ; fix degeneracy when no res
    
    fsc_plot,[0],[0],ytitle="log(Frac Diff)",xtitle="Position [Code]",xrange=xrange,yrange=yrange,/xs,/ys,$
              xticklen=0.05,/nodata,position=[0.25,0.15,0.95,0.4],/noerase,yticks=3
    
    fsc_plot,xrange,[-2.0,-2.0],line=0,color=fsc_color('gray'),thick=!p.thick-1.0,/overplot
    
    if (snap ne 0) then $
      fsc_plot,x_tracer[w],alog10(res[w]),psym=-4,line=line,thick=thick,symsize=0.2,/overplot  
              
  end_PS, pngResize=50, /deletePS

end
