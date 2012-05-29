; tracersShocktube.pro
; dev for tracer particles (shocktube and converging flow tests)
; dnelson feb.2012

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

pro plotShocktube

  snaps = indgen(1001)
  
  foreach snap,snaps do begin

    print,'Running: snap ['+str(snap)+']'
    
    ; config
    workingPath = '/n/home07/dnelson/dev.tracer/'
    snapPath    = workingPath + 'convFlow.cs0.02.L5.1e2.f1/output.ga53/'
    plotBase    = 'convFlow.cs0.02.L5.1e2.f1'
    
    gamma = 5.0/3.0
    ndims = 2
    nNGB  = 32 ; number of neighbors to use with external CalcHSML routine
               ; which finds smoothing lengths for tophat density estimate
    
    xrange = [0.0,20.0]
    
    ; load
    h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
    
    pos_gas    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='x')
    dens_gas   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='rho')
    u_gas      = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='u')
    pres_gas   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pressure')
    mass_gas   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='mass')

    entr_gas   = pres_gas / dens_gas^gamma
    
    pos_tracer  = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='x')
    
    y_gas    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='y')
    y_tracer = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='y')
    
    ; sort on x-position
    sort_ind_gas    = sort(pos_gas)
    sort_ind_tracer = sort(pos_tracer)
    
    pos_gas  = pos_gas[sort_ind_gas]
    dens_gas = dens_gas[sort_ind_gas]
    u_gas    = u_gas[sort_ind_gas]
    pres_gas = pres_gas[sort_ind_gas]
    mass_gas = mass_gas[sort_ind_gas]
    entr_gas = entr_gas[sort_ind_gas]
    
    pos_tracer  = pos_tracer[sort_ind_tracer]
    
    y_gas    = y_gas[sort_ind_gas]
    y_tracer = y_tracer[sort_ind_tracer]
    
    ; load horizontal velocity field
    vx_gas    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='velx')
    vx_tracer = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='velx')
    
    vx_gas    = vx_gas[sort_ind_gas]
    vx_tracer = vx_tracer[sort_ind_tracer]    
    
    ; load vertical velocity field
    vy_gas    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='vely')
    vy_tracer = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='vely')
    
    vy_gas    = vy_gas[sort_ind_gas]
    vy_tracer = vy_tracer[sort_ind_tracer]    
    
    ; total gas energy
    tot_e = total(0.5 * (vx_gas^2.0 + vy_gas^2.0) * double(mass_gas) + mass_gas * u_gas)
    print,tot_e
    
    ; evaluate tracer density
    pos_tracer3d  = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='pos')
    mass_tracer = replicate(total(mass_gas) / h.nPartTot[3], h.nPartTot[3])
  
    eval_dens = estimateDensityTophat(pos_tracer3d,mass=mass_tracer,ndims=ndims,nNGB=nNGB,boxSize=0)
    eval_dens = eval_dens[sort_ind_tracer]
    pos_tracer3d = pos_tracer3d[*,sort_ind_tracer]
    
    ; subselect tophat densities to avoid top/bottom box boundaries
    wTrDens = where(abs(pos_tracer3d[1,*]-1.1) lt 0.05,count)

    ; debugging: load y,z
    ;  z_gas    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='z')
    ;  y_gas = y_gas[sort_ind_gas]
    ;  z_gas = z_gas[sort_ind_gas]
    ;  
    ;  y_tracer = y_tracer[sort_ind_tracer]
    ;  
    ;  ; collapse gas onto x-coordinate, check for any differences vs y or z
    ;  num_per = 10
    ;  
    ;  num_gas_new = n_elements(pos_gas) / num_per
    ;  pos_gas_new = fltarr(3,num_gas_new)
    ;  
    ;  for i=0,num_gas_new-1 do begin
    ;    pos_gas_new[0,i] = mean(pos_gas[i*num_per:(i+1)*num_per-1])
    ;    pos_gas_new[1,i] = mean(y_gas[i*num_per:(i+1)*num_per-1])
    ;    pos_gas_new[2,i] = mean(z_gas[i*num_per:(i+1)*num_per-1])
    ;  endfor
    
    ; plot - fluid and tracer quantities and positions
    start_PS, workingPath + plotBase + '_gastr_'+string(snap,format='(I04)')+'.eps',ys=6.0
    
      !p.multi = [0,1,3]
    
      psym    = 0
      symsize = 0.4
    
      ; density comparison
      yrange = [0.8,1.6] ;cs0.5=[0.5,2.5] cs2=[0.0,4.5] cs0.02=[0.8,1.6]
      !y.margin = [-6.0,2.0]
      fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
                xtitle="",ytitle="Density / Energy",$
                title="snap="+str(snap)+" time="+string(h.time,format='(f5.2)')+""
      
      fsc_plot, pos_gas,dens_gas,psym=psym,line=0,symsize=symsize,/overplot,color=getColor(0)
      fsc_plot, pos_tracer3d[0,wTrDens],eval_dens[wTrDens],$
        psym=8,symsize=symsize,/overplot,color=getColor(3)
      
      fsc_plot, pos_gas,u_gas,psym=psym,line=0,symsize=symsize,/overplot,color=getColor(4) 
      fsc_plot, pos_gas,entr_gas,psym=psym,line=1,symsize=symsize,/overplot,color=getColor(5)
  
      ; legend
      strings = ['gas dens','tracer dens','gas u','gas '+textoidl('P/\rho^\gamma')]
      legend,strings,textcolors=getColor([0,3,4,5],/name),$
        box=0,margin=0.25,/right,charsize=!p.charsize-0.4
  
      ; vx/vy comparison
      yrange = [-0.025,0.025]  ;cs0.5=[-3.4,-3.4] cs2=[-3.4,-3.4 cs0.02=[-0.025,0.025]
      !y.margin = [-8.0,8.0]
      fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
                xtitle="",ytitle="Velocity"
                
      fsc_plot, pos_gas,vx_gas,psym=psym,line=0,symsize=symsize,/overplot,color=getColor(0) 
      fsc_plot, pos_tracer,vx_tracer,psym=8,symsize=symsize,/overplot,color=getColor(1)
  
      fsc_plot, pos_gas,vy_gas,psym=psym,line=0,symsize=symsize,/overplot,color=getColor(2) 
      fsc_plot, pos_tracer,vy_tracer,psym=8,symsize=symsize,/overplot,color=getColor(3)
  
      ; legend
      strings = ['gas vx','tracer vx','gas vy','tracer vy']
      legend,strings,textcolors=getColor([0,1,2,3],/name),$
        box=0,margin=0.25,/right,charsize=!p.charsize-0.4
    
      ; scatterplot positions
      yrange = [0.0,2.0]
      !y.margin = [4.0,10.0]
      fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
                xtitle="Position",ytitle="y";,position=[0.06,0.2,0.97,0.86]
      
      plotsym,0
      fsc_plot,pos_gas,y_gas,psym=8,symsize=0.5,/overplot,color=getColor(1)
      plotsym,0,/fill
      fsc_plot,pos_tracer,y_tracer,psym=8,symsize=0.2,/overplot,color=getColor(3)
    
    end_PS, pngResize=50, /deletePS
    
  endforeach
 
end

; plotShocktubeMeshDens(): same as above but with mesh and projected density output files

pro plotShocktubeMeshDens

  ; config
  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath    = workingPath + 'convFlow.cs5.L5.1e2.f1/output/'
  plotBase    = 'convFlow.cs5.L5.1e2.f1'
  
  yrange1 = [0.0,8.5] ;cs0.5=[0.5,2.5] cs1.5/2=[0.0,4.5] cs0.02=[0.9,1.1] cs5=[0.0,8.5]
  yrange2 = [-3.4,3.4]  ;cs0.5=[-3.4,-3.4] cs1.5/2=[-3.4,-3.4] cs0.02=[-0.025,0.025]
  
  gamma = 2.0
  ndims = 2
  nNGB_2D  = 128 ; number of neighbors to use with external CalcHSML routine
                 ; which finds smoothing lengths for tophat density estimate
  nNGB_1D  = 32
  trPartType = 3
  
  xrange = [0.0,20.0]

  snaps = indgen(1001)
  
  ; equivalent tracer mass
  sP = { simPath:snapPath, snap:0 }
  h = loadSnapshotHeader(sP=sP)
  
  mass_gas    = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
  trMassConst = total(mass_gas) / h.nPartTot[3]
  
  foreach snap,snaps do begin

    print,'Running: snap ['+str(snap)+']'
    sP.snap = snap
    
    ; load
    h = loadSnapshotHeader(sP=sP,/verbose)
    
    pos_gas    = loadSnapshotSubset(sP=sP,partType='gas',field='x')
    dens_gas   = loadSnapshotSubset(sP=sP,partType='gas',field='rho')
    u_gas      = loadSnapshotSubset(sP=sP,partType='gas',field='u')
    pres_gas   = loadSnapshotSubset(sP=sP,partType='gas',field='pressure')
    ;mass_gas   = loadSnapshotSubset(sP=sP,partType='gas',field='mass')

    entr_gas   = pres_gas / dens_gas^gamma
    
    pos_tracer  = loadSnapshotSubset(sP=sP,partType=trPartType,field='x')
    
    y_gas    = loadSnapshotSubset(sP=sP,partType='gas',field='y')
    y_tracer = loadSnapshotSubset(sP=sP,partType=trPartType,field='y')
    
    ; sort on x-position
    sort_ind_gas    = sort(pos_gas)
    sort_ind_tracer = sort(pos_tracer)
    
    pos_gas  = pos_gas[sort_ind_gas]
    dens_gas = dens_gas[sort_ind_gas]
    u_gas    = u_gas[sort_ind_gas]
    pres_gas = pres_gas[sort_ind_gas]
    ;mass_gas = mass_gas[sort_ind_gas]
    entr_gas = entr_gas[sort_ind_gas]
    
    pos_tracer  = pos_tracer[sort_ind_tracer]
    
    y_gas    = y_gas[sort_ind_gas]
    y_tracer = y_tracer[sort_ind_tracer]
    
    ; load horizontal velocity field
    vx_gas    = loadSnapshotSubset(sP=sP,partType='gas',field='velx')
    vx_tracer = loadSnapshotSubset(sP=sP,partType=trPartType,field='velx')
    
    vx_gas    = vx_gas[sort_ind_gas]
    vx_tracer = vx_tracer[sort_ind_tracer]    
    
    ; load vertical velocity field
    vy_gas    = loadSnapshotSubset(sP=sP,partType='gas',field='vely')
    vy_tracer = loadSnapshotSubset(sP=sP,partType=trPartType,field='vely')
    
    vy_gas    = vy_gas[sort_ind_gas]
    vy_tracer = vy_tracer[sort_ind_tracer]    
    
    ; total gas energy
    ;tot_e = total(0.5 * (vx_gas^2.0 + vy_gas^2.0) * double(mass_gas) + mass_gas * u_gas)
    
    ; evaluate tracer density
    pos_tracer3d  = loadSnapshotSubset(sP=sP,partType=trPartType,field='pos')
    
    eval_dens = estimateDensityTophat(pos_tracer3d,mass=trMassConst,ndims=ndims,nNGB=nNGB_2D,boxSize=0)
    eval_dens = eval_dens[sort_ind_tracer]
    pos_tracer3d = pos_tracer3d[*,sort_ind_tracer]
    
    ; subselect tophat densities to avoid top/bottom box boundaries
    wTrDens = where(abs(pos_tracer3d[1,*]-1.1) lt 0.05,count)
    
    ; re-calculate tracer density using one-dimensional CalcHSML
    eval_dens_1d = estimateDensityTophat(pos_tracer3d[*,wTrDens],mass=trMassConst,$
                                         ndims=1,nNGB=nNGB_1D,boxSize=0)

    ; plot - fluid and tracer quantities and positions
    start_PS, workingPath + plotBase + '_gastr_'+string(snap,format='(I04)')+'.eps',xs=10.0,ys=5.0
    
      !p.multi = [0,1,4]
      !p.charsize += 0.4
    
      psym    = 0
      symsize = 0.4
    
      ; density comparison
      !y.margin = [-6.0,2.0]
      fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange1,/xs,/ys,$
                xtitle="",ytitle="Density / Energy",$
                title="",xtickname=replicate(' ',10),yticklen=0.005
      
      fsc_plot, pos_gas,dens_gas,psym=psym,line=0,symsize=symsize,/overplot,color=getColor(0)
      fsc_plot, pos_tracer3d[0,wTrDens],eval_dens[wTrDens],psym=8,symsize=symsize,/overplot,color=getColor(3)
      fsc_plot, pos_tracer3d[0,wTrDens],5*eval_dens_1d,psym=8,symsize=symsize,/overplot,color=getColor(7)
      
      fsc_plot, pos_gas,u_gas,psym=psym,line=0,symsize=symsize,/overplot,color=getColor(4) 
      fsc_plot, pos_gas,entr_gas,psym=psym,line=1,symsize=symsize,/overplot,color=getColor(5)
  
      ; legend
      strings = ['gas '+textoidl('\rho'),'tracer '+textoidl('\rho_{2D}'),'tracer '+textoidl('\rho_{1D}'),$
                 'gas u','gas '+textoidl('P/\rho^\gamma')]
      legend,strings,textcolors=getColor([0,3,7,4,5],/name),$
        box=0,position=[17.5,8.2],charsize=!p.charsize-0.6
  
      ; vx/vy comparison
      !y.margin = [-8.0,6.0]
      fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange2,/xs,/ys,$
                xtitle="",ytitle="Velocity",xtickname=replicate(' ',10),yticklen=0.005
                
      fsc_plot, pos_gas,vx_gas,psym=psym,line=0,symsize=symsize,/overplot,color=getColor(0) 
      fsc_plot, pos_tracer,vx_tracer,psym=8,symsize=symsize,/overplot,color=getColor(1)
  
      ;fsc_plot, pos_gas,vy_gas,psym=psym,line=0,symsize=symsize,/overplot,color=getColor(2) 
      ;fsc_plot, pos_tracer,vy_tracer,psym=8,symsize=symsize,/overplot,color=getColor(3)
  
      ; legend
      strings = ['gas '+textoidl('v_x'),'tracer '+textoidl('v_x')];,'gas vy','tracer vy']
      legend,strings,textcolors=getColor([0,1],/name),$
        box=0,position=[17.5,3.1],charsize=!p.charsize-0.6
    
      ; plot mesh with tracers over
      yrange = [0.0,2.0]
      !y.margin = [-2.0,8.0]
      fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
                xtitle="",ytitle="y",xtickname=replicate(' ',10),ytickv=[0.5,1.0,1.5],yticks=2,yticklen=0.005
     
      plotVoronoi2D, snapPath+'voronoi_mesh_', snap, /overPlot
       
      ;plotsym,0
      ;fsc_plot,pos_gas,y_gas,psym=8,symsize=0.5,/overplot,color=getColor(1)
      plotsym,0,/fill
      fsc_plot,pos_tracer,y_tracer,psym=8,symsize=0.2,/overplot,color=getColor(3)
    
      ; projected density
      !y.margin = [4.0,2.0]
      fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
               xtitle="Position",ytitle="y",ytickv=[0.5,1.0,1.5],yticks=2,yticklen=0.005
               
      plotDensityField, snapPath,snap,/psOut,/overPlot,minMax=yrange1
    
    end_PS, pngResize=50
    
  endforeach
 
end

; compWaveOverDens(): compare maximum tracer overdensity vs time for different methods

pro compWaveOverDens

  ; config
  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath    = workingPath + 'convFlow.cs0.02.L5.1e2.f1/output/'
  plotBase    = 'convFlow.cs0.02.L5.1e2.f1'
  
  yrange = [0.5,2.0] ;cs5=[0,16] cs0.02/cs0.5=[0.5,2.0] cs1.5=[0.5,4.0] cs5f4 or cs5_1e3=[0.0,30.0]
  
  nNGB_1D = [8,16,32,64]
  nNGB_2D = [32,64,128,256]

  snaps = indgen(1001) ;indgen(50)*20 ;
  
  ; arrays
  trOverGas_1D = fltarr(n_elements(nNGB_1D),n_elements(snaps))
  trOverGas_2D = fltarr(n_elements(nNGB_2D),n_elements(snaps))
  
  times = fltarr(n_elements(snaps))
  
  ; equivalent tracer mass
  h = loadSnapshotHeader(snapPath,snapNum=0)
  
  mass_gas    = loadSnapshotSubset(snapPath,snapNum=0,partType='gas',field='mass')
  trMassConst = total(mass_gas) / h.nPartTot[3]
  
  foreach snap,snaps,j do begin

    print,'Running: snap ['+str(snap)+']'
    
    ; load
    h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
    
    times[j] = h.time
    
    pos_gas    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pos')
    dens_gas   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='rho')
    pos_tracer = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='pos')
    
    ; select gas row and sort on x-position
    wGas = where(abs(pos_gas[1,*]-1.1) lt 0.05,count)
    
    pos_gas = pos_gas[0,wGas]
    pos_gas = pos_gas[sort(pos_gas)]
    
    dens_gas = dens_gas[wGas]
    dens_gas = dens_gas[sort(pos_gas)]
    
    ; select tracer row and prepare sort on x-position
    pos_tracer = pos_tracer[*,sort(pos_tracer[0,*])]
    wTr = where(abs(pos_tracer[1,*]-1.1) lt 0.05,count)

    ; evaluate tracer density (2D)
    ndims = 2
    foreach nNGB,nNGB_2D,k do begin
      ; tophat evaluation
      eval_dens = estimateDensityTophat(pos_tracer,mass=trMassConst,ndims=ndims,nNGB=nNGB,boxSize=0)

      eval_dens = eval_dens[wTr]
      
      ; linear interpolation onto gas points
      eval_interp = interpol(eval_dens,pos_tracer[0,wTr],pos_gas)
      
      ; save maximum overdensity
      ratio = eval_interp / dens_gas
      trOverGas_2d[k,j] = max(ratio)
    endforeach

    ; zero y-coordinates
    pos_tracer[1,*] = 0.0

    ; evaluate tracer density (1D)
    ndims = 1
    foreach nNGB,nNGB_1D,k do begin
      ; tophat evaluation
      eval_dens_1d = estimateDensityTophat(pos_tracer[*,wTr],mass=trMassConst,$
                                           ndims=ndims,nNGB=nNGB,boxSize=0)
      eval_dens_1d = 5.0 * eval_dens_1d[wTr] ;5.0 = Ny/Ly (convert tracer rho to 2D, or gas rho to 1D)
      
      ; linear interpolation onto gas points
      eval_interp = interpol(eval_dens_1d,pos_tracer[0,wTr],pos_gas)
      
      ; save maximum overdensity
      ratio = eval_interp / dens_gas
      trOverGas_1d[k,j] = max(ratio)
    endforeach
    
  endforeach

  ; plot
  start_PS, workingPath + plotBase + '_overDens_'+string(snap,format='(I04)')+'.eps',ys=5.0
    
    fsc_plot, [0],[0],/nodata,xrange=minmax(times),yrange=yrange,/xs,/ys,$
              xtitle="Time",ytitle="Maximum Tracer/Gas Overdensity"
    
    foreach nNGB,nNGB_2D,k do begin
      fsc_plot,times,trOverGas_2d[k,*],line=0,color=getColor(k),/overplot
    endforeach
    
    foreach nNGB,nNGB_1D,k do begin
      fsc_plot,times,trOverGas_1d[k,*],line=0,color=getColor(k+n_elements(nNGB_2D)),/overplot
    endforeach

    ; legend
    colors  = getColor(indgen(n_elements(nNGB_1D)+n_elements(nNGB_2D)),/name)
    strings = ['2D nNGB = '+str(nNGB_2D),'1D nNGB = '+str(nNGB_1D)]
    legend,strings,textcolors=colors,box=0,margin=0.25,/right,charsize=!p.charsize-0.3

  end_PS
stop
end

; compMassDist(): compare mass distribution of tracer particles to the underlying density field of gas

pro compMassDist, snap=snap

  print,'Running: snap ['+str(snap)+']'

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
  
  pos_gas  = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pos')  
  dens_gas = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='rho')
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
