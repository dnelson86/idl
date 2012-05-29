; tracersMC_2D.pro
; dev for tracer particles (2d hydro tests)
; dnelson feb.2012

; ssF(): symsize function

function ssF, nT, base

  minSS = 0.1 ;minimum size
  maxSS = 1.0 ;maximum size

  ss = 0.4 ;base size
  
  if (nT lt base) then ss -= abs(nT-base)/(base/5.0) > minSS
  if (nT gt base) then ss += abs(nT-base)/(base/5.0) < maxSS
  
  return,ss

end

; plotConvFlowMeshDens(): same as above but with mesh and projected density output files

pro plotConvFlowMeshDens, cs=cs, fn=fn, mesh=mesh

  ; config
  ;cs = '5.0'
  ;fn = '100'
  ;mesh = 'fixed'

  ; paths
  workingPath = '/n/home07/dnelson/dev.tracerMC/'
  snapPath    = workingPath + 'convFlow.cs'+cs+'.1e2.f'+fn+'.'+mesh+'/output/'
  ;plotBase    = 'convFlow.cs'+cs+'.1e2.f'+fn+'.'+mesh+'/frames/convFlow.cs'+cs+'.1e2.f'+fn+'.'+mesh
  plotBase = 'convFlow.cs'+cs+'.1e2.f'+fn+'.'+mesh
  
  snaps = indgen(n_elements(file_search(snapPath+'snap_*.hdf5')))
  ;snaps = [400,500,600,700,800,900,1000,1024] ;indgen(1001) ;[10,100,500]  
  
  ;snapPath = workingPath + 'uniform.2d.1e2.fixed/output/'
  ;plotBase = 'frames.uniform.f1/uniform.2d.1e2.f1'
  
  TracerMCPerCell = float(fn)
  
  dispFac = 50/TracerMCPerCell>1.0 ; divisor of difference from rho=1 for TR density estimates  
  
  if (cs eq '0.02') then begin
    yrange1 = [0.95,1.05]
    yrange2 = [-0.02,0.02]
  endif
  if (cs eq '1.5') then begin
    yrange1 = [0.0,4.5]
    yrange2 = [-3.0,3.0]
  endif
  if (cs eq '5.0') then begin
    yrange1 = [0.0,4.5] ;cs0.5=[0.5,2.5] cs1.5/2=[0.0,4.5] cs0.02=[0.95,1.05] cs5=[0.0,8.5]
    yrange2 = [-3.4,3.4]  ;cs0.5=[-3.4,-3.4] cs1.5/2=[-3.4,-3.4] cs0.02=[-0.02,0.02]
  endif
  
  gamma = 2.0
  ndims = 2
  
  xrange = [0.0,20.0]
  
  ; equivalent tracer mass
  print,snapPath
  sP = { simPath:snapPath, snap:0 }
  h = loadSnapshotHeader(sP=sP)
  
  mass_gas    = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
  trMassConst = total(mass_gas) / TracerMCPerCell
  trMassConst /= (20.0*2.0) ;divide by vol -> dens
  
  foreach snap,snaps do begin

    print,'Running: snap ['+str(snap)+']'
    sP.snap = snap
    
    plotName = workingPath + plotBase + '_gastr_'+string(snap,format='(I04)')
    if (file_test(plotName+'.png')) then begin
      print,' skip'
      continue
    endif
    
    ; load
    h = loadSnapshotHeader(sP=sP,/verbose)
    
    pos_gas    = loadSnapshotSubset(sP=sP,partType='gas',field='x')
    y_gas      = loadSnapshotSubset(sP=sP,partType='gas',field='y')
    
    dens_gas   = loadSnapshotSubset(sP=sP,partType='gas',field='rho')
    u_gas      = loadSnapshotSubset(sP=sP,partType='gas',field='u')
    pres_gas   = loadSnapshotSubset(sP=sP,partType='gas',field='pressure')
    numtr_gas  = loadSnapshotSubset(sP=sP,partType='gas',field='numtr')

    entr_gas   = pres_gas / dens_gas^gamma
    
    ; sort on x-position
    sort_ind_gas    = sort(pos_gas)
    
    pos_gas  = pos_gas[sort_ind_gas]
    y_gas    = y_gas[sort_ind_gas]
    
    dens_gas  = dens_gas[sort_ind_gas]
    u_gas     = u_gas[sort_ind_gas]
    pres_gas  = pres_gas[sort_ind_gas]
    numtr_gas = numtr_gas[sort_ind_gas]
    
    entr_gas = entr_gas[sort_ind_gas]
  
    ; load horizontal velocity field
    vx_gas    = loadSnapshotSubset(sP=sP,partType='gas',field='velx')
    vx_gas    = vx_gas[sort_ind_gas]  

    ; total gas energy
    ;tot_e = total(0.5 * (vx_gas^2.0 + vy_gas^2.0) * double(mass_gas) + mass_gas * u_gas)
 
    ; scatterplot a density estimae for each of the rows separately
    dens_tr = numtr_gas / TracerMCPerCell

    ; for differing densities, adjust so that the difference is smaller (visible)
    w = where(numtr_gas gt TracerMCPerCell,count)
    if (count gt 0) then dens_tr[w] = (dens_tr[w]-1.0)/dispFac + 1.0
    w = where(numtr_gas lt TracerMCPerCell,count)
    if (count gt 0) then dens_tr[w] = 1.0 - (1.0-dens_tr[w])/dispFac

    ; collapse gas onto x-coordinate and mean tracer numbers
    num_per = 10 ; y-direction tiling
    
    num_gas_new = n_elements(pos_gas) / num_per
    numtr_new   = fltarr(num_gas_new)
    pos_gas_new = fltarr(num_gas_new)
    pos_gas_err = fltarr(num_gas_new)
    
    for i=0,num_gas_new-1 do begin
      pos_gas_new[i] = mean(pos_gas[i*num_per:(i+1)*num_per-1])
      pos_gas_err[i] = stddev(pos_gas[i*num_per:(i+1)*num_per-1])
      numtr_new[i]   = total(numtr_gas[i*num_per:(i+1)*num_per-1])
    endfor
    
    dens_tr_new = numtr_new / TracerMCPerCell / float(num_per)
    
    ; plot - fluid and tracer quantities and positions
    start_PS,plotName+'.eps',xs=10.0,ys=5.0
    
      !p.multi = [0,1,4]
      !p.charsize += 0.4
    
      psym    = 0
      symsize = 0.4
    
      ; density comparison
      !y.margin = [-6.0,2.0]
      fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange1,/xs,/ys,$
                xtitle="",ytitle="Density / Energy",$
                title="",xtickname=replicate(' ',10),yticklen=0.005
      
      fsc_plot, pos_gas,u_gas,psym=psym,line=0,symsize=symsize,/overplot,color=getColor(4) 
      fsc_plot, pos_gas,entr_gas,psym=psym,line=1,symsize=symsize,/overplot,color=getColor(5)      
      
      fsc_plot, pos_gas,dens_gas,psym=psym,line=0,symsize=symsize,/overplot,color=getColor(0)
      fsc_plot, pos_gas,dens_tr,psym=8,symsize=symsize,/overplot,color=getColor(3)
      fsc_plot, pos_gas_new,dens_tr_new,psym=8,symsize=symsize,/overplot,color=getColor(1)
  
      ; legend
      strings = ['gas '+textoidl('\rho'),'tracer indiv '+textoidl('\rho_'),'tracer mean '+textoidl('\rho_'),$
                 'gas u','gas '+textoidl('P/\rho^\gamma')]
      legend,strings,textcolors=getColor([0,3,1,4,5],/name),$
        box=0,margin=0.25,/right,charsize=!p.charsize-0.6
  
      ; vx/vy comparison
      !y.margin = [-8.0,6.0]
      fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange2,/xs,/ys,$
                xtitle="",ytitle="Velocity",xtickname=replicate(' ',10),yticklen=0.005
                
      fsc_plot, pos_gas,vx_gas,psym=psym,line=0,symsize=symsize,/overplot,color=getColor(0) 
  
      ; legend
      strings = ['gas '+textoidl('v_x')];,'gas vy','tracer vy']
      legend,strings,textcolors=getColor([0],/name),box=0,margin=0.25,/right,charsize=!p.charsize-0.6
    
      ; plot mesh with tracers over
      yrange = [0.0,2.0]
      !y.margin = [-2.0,8.0]
      fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
                xtitle="",ytitle="y",xtickname=replicate(' ',10),ytickv=[0.5,1.0,1.5],yticks=2,yticklen=0.005
     
      plotVoronoi2D, snapPath+'voronoi_mesh_', snap, /overPlot
       
      ; dot of size proportional to number of tracers
      for i=0,n_elements(pos_gas)-1 do begin
        if (numtr_gas[i] gt 0) then begin
          symsize = numtr_gas[i]/TracerMCPerCell * 0.2 ;ssF(numtr_gas[i],TracerMCPerCell)
          fsc_plot,[pos_gas[i]],[y_gas[i]],psym=8,symsize=symsize,/overplot,color=getColor(3)
        endif
      endfor
      
      ; projected density
      !y.margin = [4.0,2.0]
      fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
               xtitle="Position",ytitle="y",ytickv=[0.5,1.0,1.5],yticks=2,yticklen=0.005
               
      plotDensityField, snapPath,snap,/psOut,/overPlot,minMax=yrange1
    
    end_PS, pngResize=50
    
  endforeach
 
end

; plotTotNumTracersWithTime():

pro plotTotNumTracersWithTime

  ; config
  workingPath = '/n/home07/dnelson/dev.tracerMC/'
  snapPath    = workingPath + 'convFlow.cs0.02.1e2.f1.fixed/output/'
  plotBase    = 'convFlow.cs0.02.1e2.f1.fixed'
  
  snaps = indgen(101)
  
  ; arrays
  numTracers = []
  times = []
  
  sP = { simPath: snapPath, snap: 0 }
  
  foreach snap,snaps do begin

    print,'Running: snap ['+str(snap)+']'
    sP.snap = snap
    
    ; load
    h = loadSnapshotHeader(sP=sP,/verbose)
    
    ntr = loadSnapshotSubset(sP=sP,partType='gas',field='numtr')
    
    ; store
    numTracers = [numTracers,total(ntr)]
    print,total(ntr)
    times = [times,h.time]

  endforeach
  
  ; plot
  start_PS, workingPath + plotBase + '_numtr.eps'
  
    xrange = minmax(times)
    yrange = [900,1100]
    yrange = [4000,4100]
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="Time",ytitle="Total Number of Tracers"

    fsc_plot,times,numTracers,psym=-8,/overplot
    
  end_PS
end

; plotTotErrWithTime(): use reduced chi-squared to estimate tracer density error

pro plotTotErrorWithTime

  ; config
  fnSet = ['1','100','1e4'] ;1e4,1e5
  cs    = '0.02'
  mesh  = 'fixed'
  
  workingPath = '/n/home07/dnelson/dev.tracerMC/'
  plotBase    = 'convFlow.error'
  
  strings = []
  
  sP = { simPath: '', snap: 0 }

  ; start plot
  start_PS, workingPath + plotBase + '_cs'+cs+'_'+mesh+'.eps', xs=10.0, ys=6.0
  
    xrange = [0.0,10.0]
    yrange = [1e-7,1.0] ;fixed=[1e-7,1.0] moving=[5e-9,3e-2]
    
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
              xtitle="Time",ytitle=textoidl("\Sigma ( \rho_{tr} - \rho_{gas} )^2 / N")

    foreach fn,fnSet,k do begin
    
      snapPath = workingPath + 'convFlow.cs'+cs+'.1e2.f'+fn+'.'+mesh+'/output/'
      print,'Running: f'+fn+' cs'+cs+' '+mesh
      sP.simPath = snapPath
      
      ; determine number of snapshots
      snaps = indgen(n_elements(file_search(snapPath+'snap_*.hdf5')))
      ;snaps = indgen(90)
      
      ; arrays
      chiSq1 = [] ;unstacked (f=fn)
      chiSq2 = [] ;stacked (f=10fn)
      times = []
      
      foreach snap,snaps do begin
    
        if (snap mod 100 eq 0) then print,' snap ['+str(snap)+']'
        
        ; load
        sP.snap = snap
        h = loadSnapshotHeader(sP=sP,/verbose)
        dP = h.flagDoublePrecision
        
        ntr_gas  = loadSnapshotSubset(sP=sP,doublePrec=dP,partType='gas',field='numtr')  
        pos_gas  = loadSnapshotSubset(sP=sP,doublePrec=dP,partType='gas',field='x')
        dens_gas = loadSnapshotSubset(sP=sP,doublePrec=dP,partType='gas',field='rho')
        y_gas    = loadSnapshotSubset(sP=sP,doublePrec=dP,partType='gas',field='y')

        ; sort on x-position
        sort_ind_gas = sort(pos_gas)
        pos_gas  = pos_gas[sort_ind_gas]
        dens_gas = dens_gas[sort_ind_gas]
        ntr_gas  = long(ntr_gas[sort_ind_gas])
        y_gas    = y_gas[sort_ind_gas]
        
        ; select one row for unstacked
        wRow = where(abs(y_gas-1.1) lt 0.01,wRowCount)
        
        ; collapse gas onto x-coordinate and mean tracer numbers
        num_per = 10 ; y-direction tiling
        
        num_gas_new  = n_elements(pos_gas) / num_per
        ntr_new      = lonarr(num_gas_new)
        pos_gas_new  = fltarr(num_gas_new)
        dens_gas_new = fltarr(num_gas_new)
        
        for i=0,num_gas_new-1 do begin
          pos_gas_new[i]  = mean(pos_gas[i*num_per:(i+1)*num_per-1])
          dens_gas_new[i] = mean(dens_gas[i*num_per:(i+1)*num_per-1])
          ntr_new[i]      = total(ntr_gas[i*num_per:(i+1)*num_per-1])
        endfor
        
        ; tracer density estimates
        dens_tr     = ntr_gas / double(fn)
        dens_tr_new = ntr_new / double(fn) / double(num_per)        
        
        ; calculate reduced chi^2 and store
        ; TODO normalize by stddev(tracer density)?
        chiSq1Indiv = total( (dens_tr[wRow]-dens_gas[wRow])^2.0 ) / (wRowCount+1.0)
        chiSq2Indiv = total( (dens_tr_new-dens_gas_new)^2.0 ) / (num_gas_new+1.0)

        chiSq1 = [chiSq1,chiSq1Indiv]
        chiSq2 = [chiSq2,chiSq2Indiv]
        times = [times,h.time]
    
      endforeach
      
      ; overplot for this sim
      fsc_plot,times,chiSq1,line=0,color=getColor(2*k+0),/overplot
      fsc_plot,times,chiSq2,line=0,color=getColor(2*k+1),/overplot
      strings = [strings,textoidl('10^'+string(alog10(fn*1.0),format='(i1)')) +' '+mesh+' mesh',$
                         textoidl('10^'+string(alog10(fn*10.0),format='(i1)'))+' '+mesh+' mesh']
    
    endforeach
 
    ; legend
    legend,strings,textcolors=getColor(indgen(n_elements(fnSet)*2),/name),$
        box=0,margin=0.25,/right,/bottom,charsize=!p.charsize-0.3
        
  ; end plot
  end_PS
stop
end
