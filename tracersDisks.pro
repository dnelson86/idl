; tracersDisks.pro
; dnelson
; jan 2012
;
; dev for tracer particles related to disk (2d/3d) tests

; fixBoundaryDisk(): load one of Diego's circumstellar 2d disk ICs and strip out the boundary cells
;                    which have duplicate, negative IDs, and add a background grid in the central 
;                    hole and outer empty regions

pro fixBoundaryDisk, addBackGrid=addBackGrid

  in_file = 'ics.dat.orig'
  out_file = 'ics.dat.new'
  
  ; master switch - use constant continuation for center of disk or vacuum
  useCenVacuum = 0
  
  ; central parameters
  Central_Dens = 6.3662e-4
  
  csnd = 0.05 * sqrt(0.4)
  gamma = 5.0/3.0
  Central_U    = csnd^2.0 / (gamma*(gamma-1.0))
  
  Vacuum_Dens = 1e-7
  Vacuum_U    = 1e-7
  
  ; load
  h = loadSnapshotHeader(in_file,snapNum='none')
  print,'Read: '+str(in_file)+' ('+str(h.nPartTot[0])+' total type0 particles)'
  
  pos    = loadSnapshotSubset(in_file,snapNum='none',partType='gas',field='pos')
  u      = loadSnapshotSubset(in_file,snapNum='none',partType='gas',field='u')
  masses = loadSnapshotSubset(in_file,snapNum='none',partType='gas',field='masses')
  ids    = loadSnapshotSubset(in_file,snapNum='none',partType='gas',field='ids')
  vels   = loadSnapshotSubset(in_file,snapNum='none',partType='gas',field='vel')

  ; make selection
  w = where(ids ge 0,ncomp=ncomp)
  
  pos    = pos[*,w]
  u      = u[w]
  masses = masses[w]
  ids    = ids[w]
  vels   = vels[*,w]
  
  print,'Removed ['+str(ncomp)+'] boundary particles.'
  
  ; add background grid if requested
  if keyword_set(addBackGrid) then begin
    nBackGrid = 32 ;squared
    
    backIDStart = max(ids)+1
    
    ; arrays
    back_pos    = fltarr(3,nBackGrid*nBackGrid)
    back_u      = fltarr(nBackGrid*nBackGrid)
    back_masses = fltarr(nBackGrid*nBackGrid)
    back_ids    = lonarr(nBackGrid*nBackGrid)
    back_vels   = fltarr(3,nBackGrid*nBackGrid)
    
    ; fill in properties
    deltax = h.boxSize/nBackGrid
    deltay = h.boxSize/nBackGrid
    
    for i=0L,nBackGrid-1 do begin
      for j=0L,nBackGrid-1 do begin
        pid = i+j*nBackGrid
        
        back_pos[0,pid] = i*deltax+deltax/2.0 + randomu(seed)/100.0 
        back_pos[1,pid] = j*deltay+deltay/2.0 + randomu(seed)/100.0  
        back_pos[2,pid] = 0.0
        
        back_ids[pid]    = backIDStart + pid 
        ;back_vels[*,pid] = 0.0
      endfor
    endfor    
    
    ; radii of all particles
    boxCen = h.boxSize/2.0
    
    rad      = reform(sqrt((pos[0,*]-boxCen)*(pos[0,*]-boxCen) + $
                           (pos[1,*]-boxCen)*(pos[1,*]-boxCen)))
    back_rad = reform(sqrt((back_pos[0,*]-boxCen)*(back_pos[0,*]-boxCen) + $
                           (back_pos[1,*]-boxCen)*(back_pos[1,*]-boxCen)))
    
    ; assign continuous properties for mass and u from the disk for inner hole
    back_w = where(back_rad lt min(rad)*0.99)
    
    if (useCenVacuum eq 0) then begin
      back_u[back_w]      += Central_U
      back_masses[back_w] += Central_Dens
      print,'  Setting center to continuous values.'
    endif else begin
      back_u[back_w]      += Vacuum_U
      back_masses[back_w] += Vacuum_Dens
      print,'  Setting center to vacuum.'
    endelse

    ; outer region as vacuum
    back_w = where(back_rad gt min(rad)*1.01)

    back_u[back_w]      += Vacuum_U
    back_masses[back_w] += Vacuum_Dens
    
    ; select only those outside of disk region
    w = where(back_rad lt min(rad)*0.99 or back_rad gt max(rad)*1.01)
    
    back_pos    = back_pos[*,w]
    back_u      = back_u[w]
    back_masses = back_masses[w]
    back_ids    = back_ids[w]
    back_vels   = back_vels[*,w]
    
    ; concat disk and background
    pos    = [[pos],[back_pos]]
    u      = [u,back_u]
    masses = [masses,back_masses]
    ids    = [ids,back_ids]
    vels   = [[vels],[back_vels]]

    print,'Added ['+str(n_elements(back_ids))+'] background points from '+str(nBackGrid)+'^2 grid.'

  endif

  ; save new IC
  print,'Writing: ',out_file,' (kept '+str(n_elements(ids))+' particles)'

  writeICFileHDF5, out_file, h.boxSize, pos, vels, ids, masses, u

end

; diskGasRadProfiles(): calculate radial profiles of density, temperature, velocity of gas

pro diskGasRadProfiles

  ; config
  res = '2'
  
  numBins = 64
  
  snapRange = [0,9]
  snapStep = 2
  
  radMinMax = [0.0,30.0]
  
  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath = workingPath + 'disk3d.'+res+'/output/'
  plotBase    = 'disk3d_'+res+'_'

  units = getUnits() ;colors
  
  ; arrays
  nSnaps = (snapRange[1]-snapRange[0]+1) / snapStep
  
  radDens    = fltarr(nSnaps,numBins)
  radNumDens = fltarr(nSnaps,numBins)
  radTemp    = fltarr(nSnaps,numBins)  
  
  ; radial bins
  radBins = findgen(numBins+1)/numBins * radMinMax[1] + radMinMax[0]
  radBinC = radBins + 0.5 * (radMinMax[1]-radMinMax[0])/numBins  
  
  ; load
  k = 0
  for snap=snapRange[0],snapRange[1],snapStep do begin
    print,'snap: ',str(snap)
    h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
    
    pos    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pos')
    dens   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='rho')
    u      = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='u')
  
    ; calculate radii of particles
    boxCen = h.boxSize / 2.0
    
    rad = reform(sqrt((pos[0,*]-boxCen)*(pos[0,*]-boxCen) + $
                      (pos[1,*]-boxCen)*(pos[1,*]-boxCen)))
                      
    ; do binning
    for j=0, numBins-1,1 do begin
      w = where(rad ge radBins[j] and rad lt radBins[j+1], count)
    
      if (count ne 0.0) then begin
        radDens[k,j]    = mean(dens[w])
        radNumDens[k,j] = count
        radTemp[k,j]    = mean(u[w])
      endif
      
    endfor
    
    k += 1
    
  endfor ;snap
  
  ; plots
  start_PS, workingPath + plotBase + 'radDens_'+string(snapRange[0],format='(I3.3)')+'_'+$
            string(snapRange[1],format='(I3.3)')+'.eps'
  
    ;yrange = [max(radDens)*0.95,max(radDens)*1.01] ;2d 0-20
    ;yrange = [0.0005,0.0008] ;VTS
    yrange = [1e-6,max(radDens)*1.05] ;3d
    
    fsc_plot, [0],[0],/nodata,xrange=radMinMax+[0.0,0.2],yrange=yrange,/xs,/ys,$
              xtitle="Radius",ytitle="Gas Density",$
              title="snapStart="+str(snapRange[0])+" snapEnd="+str(snapRange[1]),$
              charsize=!p.charsize-0.5,/ylog
    
    ; inner and outer boundaries
    fsc_plot,[0.4,0.4],yrange,color=fsc_color('light gray'),/overplot
    fsc_plot,[2.5,2.5],yrange,color=fsc_color('light gray'),/overplot
    
    ; density profiles for successive snapshots
    k = 0
    legendStrs = []
    legendColors = []
    
    for snap=snapRange[0],snapRange[1],snapStep do begin
      w = where(radDens[k,*] gt 1e-5,count)
      
      if (count ne 0) then $
        fsc_plot,radBinC[w],radDens[k,w],psym=-4,symsize=0.2,/overplot,$
                 color=fsc_color(units.colors[k mod 20]),thick=!p.thick+0.5
      
      legendStrs = [legendStrs,'snap = '+str(snap)]
      legendColors = [legendColors,units.colors[k mod 20]]      
      k += 1
    endfor
    
    ; snapshot legend
    legend,legendStrs,textcolors=legendColors,box=0,margin=0.25,/right,charsize=!p.charsize-0.5

  
  end_PS

end

; getTracerPos(): return time series of tracer positions for a given set of target starting radii
;                 works for 2d and 3d (calculates radii using only x,y)

function getTracerPos, snapPath=snapPath, snapRange=snapRange, snapStep=snapStep, radTargets=radTargets, $
                       verbose=verbose

  ; arrays
  numTracers = n_elements(radTargets)
  nSnaps     = (snapRange[1]-snapRange[0]+1) / snapStep
  
  trPos     = dblarr(nSnaps,numTracers,2)
  trRad     = dblarr(nSnaps,numTracers)
  
  idTargets = lonarr(numTracers)
  times     = fltarr(nSnaps)
  
  ; verify spacing ok
  nSnaps2 = (float(snapRange[1])-snapRange[0]+1.0) / float(snapStep)
  if (nSnaps ne nSnaps2) then begin
    print,'Error: Spacing must evenly divide snapshot range.'
    return,0
  endif
  
  ; loop over requested snapshots
  k = 0
  
  for snap=snapRange[0],snapRange[1],snapStep do begin
    if keyword_set(verbose) then $
      print,'snap: ',str(snap)
    
    ; load header and store time
    h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
    times[k] = h.time  
      
    x = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='x')
    y = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='y')
    
    ids = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='ids')
    
    ; select N tracers near radTarget for first snapshot
    if (snap eq snapRange[0]) then begin      
      boxCen = h.boxSize / 2.0
      
      ; find targets
      rad = reform(sqrt((x-boxCen)*(x-boxCen) + (y-boxCen)*(y-boxCen)))
      
      for j=0,numTracers-1 do begin
        w = where(abs(rad-radTargets[j]) eq min(abs(rad-radTargets[j])),count)
      
        ; only one and get ID
        if (count gt 1) then w = w[0]
        idTargets[j] = ids[w]
        
        if keyword_set(verbose) then $
          print,'selected id='+str(idTargets[j])+' rad='+str(rad[w])
      endfor
      
    endif
    
    ; locate indices of target IDs and save positions
    for j=0,numTracers-1 do begin
      ind = where(ids eq idTargets[j])
      
      if keyword_set(verbose) then $
        print,' ['+str(j)+'] refound id at ind='+str(ind)

      trPos[k,j,0] = x[ind]
      trPos[k,j,1] = y[ind]
      trRad[k,j]   = sqrt( (x[ind]-boxCen)*(x[ind]-boxCen) + (y[ind]-boxCen)*(y[ind]-boxCen) )
    endfor
  
    k += 1
    
  endfor ;snap
  
  r = {trPos:trPos,trRad:trRad,idTargets:idTargets,times:times,boxSize:h.boxSize}
  return, r

end

; diskTracerTraj2D(): plot a few individual tracer trajectories over the disk mesh

pro diskTracerTraj2D

  res = '128'
  ts  = 'GC' ;GC/VTS
  
  snapRange = [0,99]
  snapStep  = 10
  
  radTargets = [1.0,1.1,1.5,2.0]

  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath = workingPath + 'disk2d.'+res+'.'+ts+'/output/'
  plotBase    = 'disk2d_'+ts+'_'+res+'_'
  
  units = getUnits() ;colors
  
  ; get tracer positions
  tp = getTracerPos(snapPath=snapPath, snapRange=snapRange, snapStep=snapStep, $
                       radTargets=radTargets, /verbose)

  ; plots
  start_PS, workingPath + plotBase + 'trTraj_'+string(snapRange[0],format='(I3.3)')+'_'+$
            string(snapRange[1],format='(I3.3)')+'.eps', xs=6, ys=6
  
    xyrange = [0,tp.boxSize]
              
    fsc_plot,[0],[0],/nodata,xrange=xyrange,yrange=xyrange, $
         xstyle=1,ystyle=1,position=[0.08,0.08,0.92,0.92],charsize=!p.charsize-0.5,$
         title="tracer trajectory - "+ts+res+" snapStart="+str(snapRange[0])+" snapEnd="+str(snapRange[1])

    plotVoronoi2D, snapPath+"voronoi_mesh_", [tp.boxSize,tp.boxSize], snapRange[1], /oPlot

    for i=0,n_elements(radTargets)-1 do begin
      xPts = tp.trPos[*,i,0]
      yPts = tp.trPos[*,i,1]
      
      plotsym,0,/fill
      fsc_plot, xPts, yPts, psym=-8, symsize=0.4, color=fsc_color(units.colors[i+1 mod 20]), /overplot
    endfor
  
  end_PS, pngResize=50, /deletePS
  
end

; diskTracerTraj3D(): plot a few individual tracer trajectories over contour of the disk gas density

pro diskTracerTraj3D
 
  run = '2'

  snapRange = [0,10]
  snapStep  = 1
  
  radTargets = [2.0]

  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath = workingPath + 'disk3d.'+run+'/output/'
  plotBase    = 'disk3d_'+run+'_'
  
  units = getUnits() ;colors
  
  ; get tracer positions
  tp = getTracerPos(snapPath=snapPath, snapRange=snapRange, snapStep=snapStep, $
                       radTargets=radTargets, /verbose)
            
  ; plots
  start_PS, workingPath + plotBase + 'trTraj_'+string(snapRange[0],format='(I3.3)')+'_'+$
            string(snapRange[1],format='(I3.3)')+'.eps', xs=6, ys=6
  
    fac = 50.0
    xyrange = [tp.boxSize/2.0 - tp.boxSize/fac,tp.boxSize/2.0 + tp.boxSize/fac]
              
    fsc_plot,[0],[0],/nodata,xrange=xyrange,yrange=xyrange, $
         xstyle=1,ystyle=1,position=[0.1,0.1,0.9,0.9],charsize=!p.charsize-0.5,$
         title="tracer trajectory - disk3d."+run+" snapStart="+str(snapRange[0])+" snapEnd="+str(snapRange[1]),$
         xtitle="x [kpc]",ytitle="y [kpc]"

    ;plotVoronoi2D, snapPath+"voronoi_mesh_", [tp.boxSize,tp.boxSize], snapRange[1], /oPlot

    for i=0,n_elements(radTargets)-1 do begin
      xPts = tp.trPos[*,i,0]
      yPts = tp.trPos[*,i,1]
      
      plotsym,0,/fill
      fsc_plot, xPts, yPts, psym=-8, symsize=0.4, color=fsc_color(units.colors[i+1 mod 20]), /overplot
    endfor
  
  end_PS, pngResize=50, /deletePS
                       
  stop  
  
end

; diskTracerRadiusSingle(): plot evolution of tracer radius with time

pro diskTracerRadiusSingle

  res = '128'
  ts  = 'GC' ;GC/VTS
  
  radTargets = [1.0,1.1,1.5,2.0]

  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath = workingPath + 'disk2d.'+res+'.'+ts+'/output/'
  plotBase    = 'disk2d_'+ts+'_'+res+'_'
  
  ; find number of available snapshots
  nSnaps = n_elements(file_search(snapPath+'snap_*'))
  nSnaps = fix(floor(nSnaps*10.0)/10.0)
  
  snapRange = [0,nSnaps-1]
  snapStep  = 1
  
  units = getUnits() ;colors
  
  ; get tracer positions, radii, times
  tp = getTracerPos(snapPath=snapPath, snapRange=snapRange, snapStep=snapStep, $
                       radTargets=radTargets)
  
  ; normalize times to multiples of the rotation period at R=1 (2pi)
  tp.times /= (2*!pi)
  
  ; normalize radii to their starting radii
  for i=0,n_elements(radTargets)-1 do $
    tp.trRad[*,i] /= tp.trRad[0,i]
  
  ; plot
  start_PS, workingPath + plotBase + 'trRad_'+string(snapRange[0],format='(I3.3)')+'_'+$
            string(snapRange[1],format='(I3.3)')+'.eps'
  
    yrange = [0.9905,1.010]
    
    fsc_plot,[0],[0],/nodata,xrange=[0,round(max(tp.times))],yrange=yrange, $
         xstyle=1,ystyle=1,$
         title="tracer radial evolution - "+ts+"."+res,$
         xtitle="t / P",ytitle="R / R"+textoidl("_0")

    ; loop over tracers and plot
    legendStrs = []
    legendColors = []
    
    for i=0,n_elements(radTargets)-1 do begin
      fsc_plot,tp.times,tp.trRad[*,i],line=0,color=fsc_color(units.colors[i+1 mod 20]), /overplot
      
      legendStrs = [legendStrs,'R'+textoidl("_0")+' = '+string(radTargets[i],format='(f3.1)')]
      legendColors = [legendColors,units.colors[i+1 mod 20]]    
    endfor
    
    ; R0 legend
    legend,legendStrs,textcolors=legendColors,box=0,margin=0.25,/left,charsize=!p.charsize-0.5
  
  end_PS;, pngResize=50

end

; diskTracerRadiusEveryTS(): compare moving vs. static mesh with EveryTimestep outputs

pro diskTracerRadiusEveryTS

  radTargets = [1.0,1.1,1.5,2.0]

  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath1 = workingPath + 'disk2d.everyts.movingmesh/output/'
  snapPath2 = workingPath + 'disk2d.everyts.staticmesh/output/'
  plotBase    = 'disk2d_everyTS_'
  
  units = getUnits() ;colors  
  
  ; find number of available snapshots (1)
  nSnaps1 = n_elements(file_search(snapPath1+'snap_*'))
  nSnaps2 = n_elements(file_search(snapPath2+'snap_*'))
  nSnaps1 = fix(floor(nSnaps1*10.0)/10.0)
  nSnaps2 = fix(floor(nSnaps2*10.0)/10.0)
  
  nSnaps1 = 80
  nSnaps2 = 200
  
  snapRange1 = [0,nSnaps1-1]
  snapRange2 = [0,nSnaps2-1]
  snapStep  = 1
  
  ; get tracer positions, radii, times (1)
  tp1 = getTracerPos(snapPath=snapPath1, snapRange=snapRange1, snapStep=snapStep, $
                       radTargets=radTargets,/verbose)
                       
  ; get tracer positions, radii, times (2)
  tp2 = getTracerPos(snapPath=snapPath2, snapRange=snapRange2, snapStep=snapStep, $
                       radTargets=radTargets,/verbose)
                       
  ; normalize times to multiples of the rotation period at R=1 (2pi)
  tp1.times /= (2*!pi)
  tp2.times /= (2*!pi)
  
  ; normalize radii to their starting radii
  for i=0,n_elements(radTargets)-1 do begin
    tp1.trRad[*,i] /= tp1.trRad[0,i]
    tp2.trRad[*,i] /= tp2.trRad[0,i]
  endfor
  
  ; start plot
  start_PS, workingPath + plotBase + 'trRad.eps', xs=7, ys=6
  
    xrange = [0,max([tp1.times,tp2.times])]
    yrange = [0.999,1.001]
    
    !p.multi = [0,1,2]
    
    ; plot (1)
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,$
         title="tracer radial evolution - moving mesh",$
         xtitle="t / P",ytitle="R / R"+textoidl("_0"),charsize=!p.charsize-0.5

    ; loop over tracers and plot
    legendStrs = []
    legendColors = []
    
    for i=0,n_elements(radTargets)-1 do begin
      fsc_plot,tp1.times,tp1.trRad[*,i],line=0,color=fsc_color(units.colors[i+1 mod 20]), /overplot
      
      legendStrs = [legendStrs,'R'+textoidl("_0")+' = '+string(radTargets[i],format='(f3.1)')]
      legendColors = [legendColors,units.colors[i+1 mod 20]]    
    endfor
    
    ; R0 legend
    legend,legendStrs,textcolors=legendColors,box=0,margin=0.25,/left,charsize=!p.charsize-0.5
  
    ; plot (2)
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,$
         title="tracer radial evolution - static mesh",$
         xtitle="t / P",ytitle="R / R"+textoidl("_0"),charsize=!p.charsize-0.5
         
    ; loop over tracers and plot
    for i=0,n_elements(radTargets)-1 do $
      fsc_plot,tp2.times,tp2.trRad[*,i],line=0,color=fsc_color(units.colors[i+1 mod 20]), /overplot
       
  end_PS;, pngResize=50
stop
end

; diskTracerRadius3x3(): plot diskTracerRadiusSingle for each res=128,256,512 and GC,VTS

pro diskTracerRadius3x3

  cc = [['GC','128'], ['VTS','128'], $
        ['GC','256'], ['VTS','256'], $
        ['GC','512'], ['VTS','512']]
  
  radTargets = double([1.0,1.1,1.5,2.0])
  
  units = getUnits() ;colors

  ; start plot
  workingPath = '/n/home07/dnelson/dev.tracer/'
  start_PS, workingPath + 'trRadMult.eps', xs=8, ys=10

  xrange = [0,4.0]
  ;yrange = [0.990,1.010]
  yrange = [0.998,1.002]
  
  !p.multi = [0,2,3]

  ; loop over requested configurations
  nConfigs = (size(cc))[2]
  
  for i=0,nConfigs-1 do begin
    ts  = cc[0,i]
    res = cc[1,i]
  
    snapPath = workingPath + 'disk2d.'+res+'.'+ts+'/output/'
    
    ; find number of available snapshots
    nSnaps = n_elements(file_search(snapPath+'snap_*'))
    nSnaps = fix(floor(nSnaps/10.0)*10.0)
    
    snapRange = [0,nSnaps-1]
    snapStep  = 1
    
    ; get tracer positions, radii, times
    if (nSnaps ge 10) then begin
      print,res+' '+ts,snapRange
      
      tp = getTracerPos(snapPath=snapPath, snapRange=snapRange, snapStep=snapStep, $
                        radTargets=radTargets)
      
      ; normalize times to multiples of the rotation period at R=1 (2pi)
      tp.times /= (2*!pi)
      
      ; normalize radii to their starting radii
      for j=0,n_elements(radTargets)-1 do $
        tp.trRad[*,j] /= tp.trRad[0,j]
      
      ; plot this config
      ;xrange = [0,round(max(tp.times))]
      fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,$
           xtitle="t / P",ytitle="R / R"+textoidl("_0")
         
      legendStrs = []
      legendColors = []
         
      for j=0,n_elements(radTargets)-1 do begin
        fsc_plot,tp.times,tp.trRad[*,j],line=0,color=fsc_color(units.colors[j+1 mod 20]),/overplot
        
        legendStrs = [legendStrs,'R'+textoidl("_0")+' = '+string(radTargets[j],format='(f3.1)')]
        legendColors = [legendColors,units.colors[j+1 mod 20]]    
      endfor
      
      fsc_text,xrange[1]*0.2,yrange[0]+0.001,res+ts,alignment=0.5,charsize=!p.charsize-0.5,/data
      
      ; R0 legend
      if (i eq 0) then $
      legend,legendStrs,textcolors=legendColors,box=0,margin=0.1,/left,charsize=!p.charsize-0.7
    
    endif else begin
      ; not enough data yet, empty plot
      fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,$
           xtitle="t / P",ytitle="R / R"+textoidl("_0")
           
      fsc_text,xrange[1]*0.2,yrange[0]+0.001,res+ts,alignment=0.5,charsize=!p.charsize-0.5,/data
    endelse
    
  endfor
  
  ; finish plot
  end_PS
  
  stop

end
