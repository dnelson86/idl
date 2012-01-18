; tracersDisks.pro
; dev for tracer particles related to disk (2d/3d) tests
; dnelson jan.2012

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
  dim = '3d'
  run = 'lowres.rigid'
  
  numBins = 64
  
  snapRange = [0,129]
  snapStep = 10
  
  if (dim eq '2d') then radMinMax = [0.0,3.0]
  if (dim eq '3d') then radMinMax = [0.0,25.0]
  
  ; paths
  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath = workingPath + 'disk'+dim+'.'+run+'/output/'
  plotBase    = 'disk'+dim+'_'+run+'_'

  units = getUnits() ;colors
  
  ; number of snapshots and verify spacing ok
  nSnaps = (snapRange[1]-snapRange[0]+1) / snapStep
  
  nSnaps2 = (float(snapRange[1])-snapRange[0]+1.0) / float(snapStep)
  if (nSnaps ne nSnaps2) then begin
    print,'Error: Spacing must evenly divide snapshot range.'
    return
  endif
  
  ; arrays
  radDens    = fltarr(nSnaps,numBins)
  radNumDens = fltarr(nSnaps,numBins)
  radTemp    = fltarr(nSnaps,numBins)
  
  times = fltarr(nSnaps)  
  
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
        times[k]        = h.time
      endif
      
    endfor
    
    k += 1
  endfor ;snap
  
  ; plots
  start_PS, workingPath + plotBase + 'radDens_'+string(snapRange[0],format='(I3.3)')+'_'+$
            string(snapRange[1],format='(I3.3)')+'.eps'
  
    if (dim eq '2d') then yrange = [5e-4,max(radDens)*1.01]
    if (dim eq '3d') then yrange = [1e-6,max(radDens)*1.1]
    
    if (dim eq '2d') then ylog = 0
    if (dim eq '3d') then ylog = 1
    
    fsc_plot, [0],[0],/nodata,xrange=radMinMax+[0.0,0.2],yrange=yrange,/xs,/ys,$
              xtitle="radius",ytitle="gas density [code]",$
              title="snapStart="+str(snapRange[0])+" snapEnd="+str(snapRange[1]),ylog=ylog
    
    ; inner and outer boundaries
    if (dim eq '2d') then begin
      fsc_plot,[0.4,0.4],yrange,color=fsc_color('light gray'),/overplot
      fsc_plot,[2.5,2.5],yrange,color=fsc_color('light gray'),/overplot
    endif
    
    ; density profiles for successive snapshots
    k = 0
    legendStrs = []
    legendColors = []
    
    for snap=snapRange[0],snapRange[1],snapStep do begin
      w = where(radDens[k,*] gt 1e-6,count)
      
      if (count ne 0) then $
        fsc_plot,radBinC[w],radDens[k,w],psym=-4,symsize=0.2,/overplot,$
                 color=fsc_color(units.colors[k mod 20]),thick=!p.thick+0.5
      
      if (dim eq '2d') then legendStrs = [legendStrs,'t/P = '+string(snap/10.0,format='(I2)')]
      if (dim eq '3d') then legendStrs = [legendStrs,'t = '+string(times[k],format='(f4.1)')+' Gyr']
      legendColors = [legendColors,units.colors[k mod 20]]      
      k += 1
    endfor
    
    ; snapshot legend
    legend,legendStrs,textcolors=legendColors,box=0,margin=0.25,/right,charsize=!p.charsize-0.5
  
  end_PS

end

; getTracerPos(): return time series of tracer positions for a given set of target starting radii
;                 works for 2d and 3d (calculates radii using only x,y)
;
; radTargets : select a handful of tracers closest to the specified starting radii
; radRange   : select all tracers within a radius range

function getTracerPos, snapPath=snapPath, snapRange=snapRange, snapStep=snapStep, $
                       radTargets=radTargets, radRange=radRange, verbose=verbose

  useMatch = 0 ;use match for ID location instead of where loop

  ; arrays
  nSnaps  = (snapRange[1]-snapRange[0]+1) / snapStep
  times   = fltarr(nSnaps)
  
  ; verify spacing ok
  nSnaps2 = (float(snapRange[1])-snapRange[0]+1.0) / float(snapStep)
  if (nSnaps ne nSnaps2) then begin
    print,'Error: Spacing must evenly divide snapshot range.'
    return,0
  endif
  
  ; loop over requested snapshots
  k = 0
  
  for snap=snapRange[0],snapRange[1],snapStep do begin
    if keyword_set(verbose) then if (snap mod verbose eq 0) then $
      print,'snap: ',str(snap)
    
    ; load header and store time
    h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
    times[k] = h.time  
      
    x = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='x')
    y = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='y')
    z = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='z')
    
    ids = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='ids')
    
    ; make tracer selection on first snapshot
    if (snap eq snapRange[0]) then begin      
      boxCen = h.boxSize / 2.0

      rad = reform(sqrt((x-boxCen)*(x-boxCen) + (y-boxCen)*(y-boxCen)))
      
      ; RADTARGETS = select N tracers near radTargets
      if keyword_set(radTargets) then begin
        ; arrays
        numTracers = n_elements(radTargets)
        trPos     = dblarr(nSnaps,numTracers,3)
        trRad     = dblarr(nSnaps,numTracers)
        
        idTargets = lonarr(numTracers)
        
        ; find IDs of targets
        for j=0,numTracers-1 do begin
          w = where(abs(rad-radTargets[j]) eq min(abs(rad-radTargets[j])),count)
        
          ; only one and get ID
          if (count gt 1) then w = w[0]
          idTargets[j] = ids[w]
          
          if keyword_set(verbose) then $
            print,'selected id='+str(idTargets[j])+' rad='+str(rad[w])
        endfor
      endif
      
      ; RADRANGE = select all tracers within radRange
      if keyword_set(radRange) then begin
        ; make selection
        w = where(rad ge radRange[0] and rad le radRange[1],count)
        
        if (count eq 0) then begin
          print,'Error: No tracers in specified radRange.'
          return,0
        endif else begin
          print,'Found ['+str(count)+'] tracers in radRange.'
        endelse
        
        ; arrays
        numTracers = count
        trPos     = dblarr(nSnaps,numTracers,3)
        trRad     = dblarr(nSnaps,numTracers)
        
        idTargets = lonarr(numTracers)
        
        ; find IDs of targets
        idTargets = ids[w]
      endif
      
    endif ;snap0
    
    ; for all subsequent snapshots - locate indices of target IDs and save positions
    if (useMatch eq 1) then begin
      ; use match instead
      match,ids,idTargets,ids_ind,idTargets_ind,count=count
      
      if (count ne numTracers) then begin
        print,'Error: Failed to match all targets.'
        return,0
      endif
      
    trPos[k,*,0] = x[ids_ind]
    trPos[k,*,1] = y[ids_ind]
    trPos[k,*,2] = z[ids_ind]
    trRad[k,*]   = sqrt( (x[ids_ind]-boxCen)*(x[ids_ind]-boxCen) + $
                         (y[ids_ind]-boxCen)*(y[ids_ind]-boxCen) )      
    endif else begin
      ; locate IDs using where loop
      for j=0,numTracers-1 do begin
        ind = where(ids eq idTargets[j])
        
        if keyword_set(verbose) then if (snap mod verbose eq 0) then $
          print,' ['+str(j)+'] refound id at ind='+str(ind)
    
        trPos[k,j,0] = x[ind]
        trPos[k,j,1] = y[ind]
        trPos[k,j,2] = z[ind]
        trRad[k,j]   = sqrt( (x[ind]-boxCen)*(x[ind]-boxCen) + (y[ind]-boxCen)*(y[ind]-boxCen) )
      endfor
    endelse
      
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
 
  run = 'lowres.rigid'
  
  radTargets = [2.1,5.1,10.1,15.1]

  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath = workingPath + 'disk3d.'+run+'/output/'
  plotBase    = 'disk3d_'+run+'_'
  
  units = getUnits() ;colors  
  
  ; find number of available snapshots
  nSnaps = n_elements(file_search(snapPath+'snap_*'))
  nSnaps = fix(floor(nSnaps/1.0)*1.0)
  
  snapRange = [0,nSnaps-1]
  snapStep  = 1
  
  ; get tracer positions
  tp = getTracerPos(snapPath=snapPath, snapRange=snapRange, snapStep=snapStep, $
                       radTargets=radTargets, verbose=20)
            
  ; plots
  start_PS, workingPath + plotBase + 'trTraj_'+string(snapRange[0],format='(I3.3)')+'_'+$
            string(snapRange[1],format='(I3.3)')+'.eps'
  
    ;multiplot,[2,1],/square,xgap=0.04,mTitle="tracer trajectory - disk3d."+run+" 1Gyr"
    fsc_text,0.5,0.8,"tracer trajectory - disk3d."+run+" 1Gyr",/normal,alignment=0.5
    
    ; PLOT ONE - (x,y)
    fac = 50.0
    xyrange = [tp.boxSize/2.0 - tp.boxSize/fac,tp.boxSize/2.0 + tp.boxSize/fac]
              
    fsc_plot,[0],[0],/nodata,xrange=xyrange,yrange=xyrange,xstyle=1,ystyle=1,charsize=!p.charsize-0.5,$
         xtitle="x [kpc]",ytitle="y [kpc]",position=[0.1,0.2,0.45,0.75],/noerase

    ; contour density field
    ;contourGasSurfDens, filePath=snapPath, snapNum=snapRange[1], zoomSize=10.0, gridSize=64

    ; overplot tracer trajectories
    for i=0,n_elements(radTargets)-1 do begin
      xPts = tp.trPos[*,i,0]
      yPts = tp.trPos[*,i,1]
      
      ; underplot circles at targetR
      tvcircle,radTargets[i],tp.boxSize/2.0,tp.boxSize/2.0,color=fsc_color('light gray'),/data
      
      ; circles around markers at the start and final time
      ;tvcircle,0.4,xPts[0],yPts[0],color=fsc_color('dark gray'),/data
      tvcircle,0.4,xPts[n_elements(xPts)-1],yPts[n_elements(yPts)-1],color=fsc_color('dark gray'),/data
      
      plotsym,0,/fill
      fsc_plot, xPts, yPts, psym=8, symsize=0.3, color=fsc_color(units.colors[i+1 mod 20]), /overplot
      
    endfor
    
    ; PLOT TWO (r,z)
    ;multiplot,/doyaxis
    rrange = [0.0,17.0]
    zrange = [-0.8,0.8]
    
    fsc_plot,[0],[0],/nodata,xrange=rrange,yrange=zrange,xstyle=1,ystyle=1,charsize=!p.charsize-0.5,$
         xtitle="radius [kpc]",ytitle="z [kpc]",position=[0.55,0.2,0.95,0.75],/noerase
         
    for i=0,n_elements(radTargets)-1 do begin
      xPts = tp.trPos[*,i,0]
      yPts = tp.trPos[*,i,1]
      zPts = tp.trPos[*,i,2] - tp.boxSize/2.0
      
      bc = tp.boxSize/2.0
      rPts = sqrt((xPts-bc)*(xPts-bc) + (yPts-bc)*(yPts-bc))
      
      ; underplot line at targetR
      fsc_plot,[radTargets[i],radTargets[i]],zrange,color=fsc_color('light gray'),/overplot
      
      ; circles around markers at the start and final time
      tvcircle,0.2,rPts[0],zPts[0],color=fsc_color('dark gray'),/data
      tvcircle,0.2,rPts[n_elements(rPts)-1],zPts[n_elements(zPts)-1],color=fsc_color('dark gray'),/data  
      
      plotsym,0,/fill
      fsc_plot, rPts, zPts, psym=8, symsize=0.3, color=fsc_color(units.colors[i+1 mod 20]), /overplot
      
    endfor
  
  end_PS, pngResize=50
  
end

; diskTracerRadSingle2D(): plot evolution of tracer radius with time

pro diskTracerRadSingle2D

  run = '128.GC'
  
  radTargets = [1.0,1.1,1.5,2.0]   ;2d
  ;radTargets = [1.000,1.005,1.010,1.015]
  
  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath = workingPath + 'disk2d.'+run+'/output/'
  plotBase    = 'disk2d.'+run+'_'
  
  ; find number of available snapshots
  spacing = 10.0
  
  nSnaps = n_elements(file_search(snapPath+'snap_*'))
  nSnaps = fix(floor(nSnaps/spacing)*spacing)
  
  snapRange = [0,nSnaps-1]
  snapStep  = spacing

  units = getUnits() ;colors
  
  ; get tracer positions, radii, times
  tp = getTracerPos(snapPath=snapPath, snapRange=snapRange, snapStep=snapStep, $
                       radTargets=radTargets, verbose=10)
  
  ; normalize times to multiples of the rotation period at R=1 (2pi) for 2D only
  tp.times /= (2*!pi)
  
  ; normalize radii to their starting radii
  for i=0,n_elements(radTargets)-1 do $
    tp.trRad[*,i] /= tp.trRad[0,i]
  
  ; plot
  start_PS, workingPath + plotBase + 'trRad_'+string(snapRange[0],format='(I3.3)')+'_'+$
            string(snapRange[1],format='(I3.3)')+'.eps'
  
    xrange = [0,ceil(max(tp.times))]
    yrange = [0.983,1.017]
        
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,$
         title="tracer radial evolution - disk2d "+run,$
         xtitle="t / P",ytitle="R / R"+textoidl("_0")

    fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot,thick=!p.thick-0.5

    ; loop over tracers and plot
    legendStrs = []
    legendColors = []
    
    for i=0,n_elements(radTargets)-1 do begin
      fsc_plot,tp.times,tp.trRad[*,i],line=0,color=fsc_color(units.colors[i+1 mod 20]), /overplot
      
      legendStrs = [legendStrs,'R'+textoidl("_0")+' = '+string(radTargets[i],format='(f3.1)')]
      legendColors = [legendColors,units.colors[i+1 mod 20]]    
    endfor
    
    ; R0 legend
    legend,legendStrs,textcolors=legendColors,box=0,margin=0.25,/left,charsize=!p.charsize-0.4
  
  end_PS;, pngResize=50

end

; diskTracerRadSingle3D(): plot evolution of tracer radius with time (3D)

pro diskTracerRadSingle3D

  run   = 'lowres.rigid'
  
  radTargets = [2.1,5.1,10.1,15.1]
  
  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath = workingPath + 'disk3d.'+run+'/output/'
  plotBase    = 'disk3d.'+run+'_'
  
  ; find number of available snapshots
  spacing = 1.0
  
  nSnaps = n_elements(file_search(snapPath+'snap_*'))
  nSnaps = fix(floor(nSnaps/spacing)*spacing)
  
  snapRange = [0,nSnaps-1]
  snapStep  = spacing

  units = getUnits() ;colors
  
  ; get tracer positions, radii, times
  tp = getTracerPos(snapPath=snapPath, snapRange=snapRange, snapStep=snapStep, $
                       radTargets=radTargets, verbose=20)
  
  ; normalize radii to their starting radii
  for i=0,n_elements(radTargets)-1 do $
    tp.trRad[*,i] /= tp.trRad[0,i]
  
  ; plot
  start_PS, workingPath + plotBase + 'trRad_'+string(snapRange[0],format='(I3.3)')+'_'+$
            string(snapRange[1],format='(I3.3)')+'.eps'
  
    xrange = [0,max(tp.times)]
    yrange = [0.95,1.05]
    
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,$
         title="tracer radial evolution - disk3d "+run,$
         xtitle="time [Gyr]",ytitle="R / R"+textoidl("_0")+""

    fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot,thick=!p.thick-0.5

    ; loop over tracers and plot
    legendStrs = []
    legendColors = []
    plotsym,0,/fill
    
    for i=0,n_elements(radTargets)-1 do begin
      fsc_plot,tp.times,tp.trRad[*,i],psym=-8,symsize=0.4,$
               color=fsc_color(units.colors[i+1 mod 20]),/overplot
      
      legendStrs = [legendStrs,'R'+textoidl("_0")+' = '+string(radTargets[i],format='(f4.1)')]
      legendColors = [legendColors,units.colors[i+1 mod 20]]    
    endfor
    
    ; R0 legend
    legend,legendStrs,textcolors=legendColors,box=0,margin=0.25,/right,charsize=!p.charsize-0.4
  
  end_PS;, pngResize=50

end

; diskTracerRadiusCompTwo(): compare moving vs. static mesh with EveryTimestep outputs

pro diskTracerRadiusCompTwo

  radTargets = [1.0,1.1,1.5,2.0]

  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath1 = workingPath + 'disk2d.128.GC/output/'
  snapPath2 = workingPath + 'disk2d.128.GC.novelgrad/output/'
  plotBase    = 'disk2d_novelgrad_'
  
  units = getUnits() ;colors  
  
  ; find number of available snapshots (1)
  nSnaps1 = n_elements(file_search(snapPath1+'snap_*'))
  nSnaps2 = n_elements(file_search(snapPath2+'snap_*'))
  nSnaps1 = fix(floor(nSnaps1/10.0)*10.0)
  nSnaps2 = fix(floor(nSnaps2/10.0)*10.0)
  
  snapRange1 = [0,nSnaps1-1]
  snapRange2 = [0,nSnaps2-1]
  snapStep  = 1
  
  ; get tracer positions, radii, times (1)
  tp1 = getTracerPos(snapPath=snapPath1, snapRange=snapRange1, snapStep=snapStep, $
                       radTargets=radTargets,verbose=20)
                       
  ; get tracer positions, radii, times (2)
  tp2 = getTracerPos(snapPath=snapPath2, snapRange=snapRange2, snapStep=snapStep, $
                       radTargets=radTargets,verbose=20)
                       
  ; normalize times to multiples of the rotation period at R=1 (2pi)
  tp1.times /= (2*!pi)
  tp2.times /= (2*!pi)
  
  ; normalize radii to their starting radii
  for i=0,n_elements(radTargets)-1 do begin
    tp1.trRad[*,i] /= tp1.trRad[0,i]
    tp2.trRad[*,i] /= tp2.trRad[0,i]
  endfor
  
  ; start plot
  start_PS, workingPath + plotBase + 'trRad_'+string(snapRange2[0],format='(I3.3)')+'_'+$
            string(snapRange2[1],format='(I4.4)')+'.eps', xs=7, ys=6
  
    xrange = [0.0,min([max(tp1.times),max(tp2.times)])]
    yrange = [0.98,1.02]
    
    !p.multi = [0,1,2]
    
    ; plot (1)
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,$
         title="tracer radial evolution - 128.GC.velgrad",$
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
    yrange = [0.88,1.12]
    
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,$
         title="tracer radial evolution - 128.GC.novelgrad",$
         xtitle="t / P",ytitle="R / R"+textoidl("_0"),charsize=!p.charsize-0.5
         
    ; loop over tracers and plot
    for i=0,n_elements(radTargets)-1 do $
      fsc_plot,tp2.times,tp2.trRad[*,i],line=0,color=fsc_color(units.colors[i+1 mod 20]), /overplot
       
  end_PS

end

; diskTracerRadiusEveryTS(): three different timescale views of an "everyTS" output

pro diskTracerRadiusEveryTS

  radTargets = [1.0,1.1,1.5,2.0]
  ;radTargets = [1.000,1.005,1.010,1.015]

  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath = workingPath + 'disk2d.everyTS.novelgrad/output/'
  plotBase    = 'disk2d_everyTS_novelgrad_'
  
  ; find number of available snapshots
  nSnaps = n_elements(file_search(snapPath+'snap_*'))
  nSnaps = fix(floor(nSnaps/10.0)*10.0)
  
  snapRange = [0,nSnaps-1]
  snapStep  = 1

  units = getUnits() ;colors
  
  ; get tracer positions, radii, times
  tp = getTracerPos(snapPath=snapPath, snapRange=snapRange, snapStep=snapStep, $
                       radTargets=radTargets, verbose=100)
  
  ; normalize times to multiples of the rotation period at R=1 (2pi)
  tp.times /= (2*!pi)
  
  ; normalize radii to their starting radii
  for i=0,n_elements(radTargets)-1 do $
    tp.trRad[*,i] /= tp.trRad[0,i]
  
  ; plot
  start_PS, workingPath + plotBase + 'trRad_'+string(snapRange[0],format='(I3.3)')+'_'+$
            string(snapRange[1],format='(I4.4)')+'.eps', xs=7, ys=8
  
    !p.multi = [0,1,3]
    
    ; PLOT ONE
    xrange = [0,0.02]
    yrange = [0.99999,1.00001]
    
    !x.margin += [5.0,3.0]
    
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,$
         title="tracer radial evolution - GC.128.everyTS",$
         xtitle="t / P",ytitle="R / R"+textoidl("_0")

    fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot,thick=!p.thick-0.5

    ; loop over tracers and plot
    legendStrs = []
    legendColors = []
    
    for i=0,n_elements(radTargets)-1 do begin
      fsc_plot,tp.times,tp.trRad[*,i],line=0,color=fsc_color(units.colors[i+1 mod 20]), /overplot
      
      legendStrs = [legendStrs,'R'+textoidl("_0")+' = '+string(radTargets[i],format='(f3.1)')]
      legendColors = [legendColors,units.colors[i+1 mod 20]]    
    endfor
    
    ; R0 legend
    legend,legendStrs,textcolors=legendColors,box=0,margin=0.25,/left,charsize=!p.charsize-0.7
  
    ; PLOT TWO
    xrange = [0,0.2]
    yrange = [0.999,1.001]
    
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,$
         xtitle="t / P",ytitle="R / R"+textoidl("_0")
    
    for i=0,n_elements(radTargets)-1 do $
      fsc_plot,tp.times,tp.trRad[*,i],line=0,color=fsc_color(units.colors[i+1 mod 20]), /overplot
      
    ; PLOT THREE
    xrange = [0,2.0]
    yrange = [0.995,1.015]
    
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,$
         xtitle="t / P",ytitle="R / R"+textoidl("_0")
    
    for i=0,n_elements(radTargets)-1 do $
      fsc_plot,tp.times,tp.trRad[*,i],line=0,color=fsc_color(units.colors[i+1 mod 20]), /overplot
  
  end_PS;, pngResize=50

end

; diskTracerRadius3(): plot diskTracerRadiusSingle for each res=128,256,512

pro diskTracerRadius3

  cc = [['GC','128'], ['GC','256'], ['GC','512']]
  labels = ['128x384','256x768','512x1536']
  
  radTargets = double([1.0,1.1,1.5,2.0])
  
  units = getUnits() ;colors

  ; start plot
  workingPath = '/n/home07/dnelson/dev.tracer/'
  start_PS, workingPath + 'trRad_resolution.eps', xs=7, ys=8

  xrange = [0,7.0]
  yrange = [0.998,1.002]
  
  !p.multi = [0,1,3]

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
      fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,$
           xtitle="t / P",ytitle="R / R"+textoidl("_0")
         
      legendStrs = []
      legendColors = []
         
      for j=0,n_elements(radTargets)-1 do begin
        fsc_plot,tp.times,tp.trRad[*,j],line=0,color=fsc_color(units.colors[j+1 mod 20]),/overplot
        
        legendStrs = [legendStrs,'R'+textoidl("_0")+' = '+string(radTargets[j],format='(f3.1)')]
        legendColors = [legendColors,units.colors[j+1 mod 20]]    
      endfor
      
      fsc_text,xrange[1]*0.15,yrange[0]+0.001,labels[i],alignment=0.5,charsize=!p.charsize-0.5,/data
      
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
  
end

; diskTracersDisp(): plot dispersion of radii (histogram) of a large population of tracers

pro diskTracersDisp

  res = '128'
  ts  = 'GC.novelgrad'
  
  radRanges = [[1.0,1.1],[1.5,1.6],[2.0,2.1]]
  
  timeTargets = [10.0,20.0,30.0,40.0] ;t/P

  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath = workingPath + 'disk2d.'+res+'.'+ts+'/output/'
  plotBase    = 'disk2d_'+ts+'_'+res+'_'
  
  numRanges = (size(radRanges))[2]
  
  ; find number of available snapshots
  spacing = 2.0
  nSnaps = n_elements(file_search(snapPath+'snap_*'))
  nSnaps = fix(floor(nSnaps/spacing)*spacing)
  
  snapRange = [0,nSnaps-1]
  snapStep  = spacing
  
  units = getUnits() ;colors
    
  ; plot time series
  ;if 0 then begin
  start_PS, workingPath + plotBase + 'trDisp_'+string(snapRange[0],format='(I3.3)')+'_'+$
            string(snapRange[1],format='(I3.3)')+'.eps', xs=7, ys=9
    
    !p.multi = [0,1,numRanges]
    
    for j=0,numRanges-1 do begin
      ; get tracer positions, radii, times
      tp = getTracerPos(snapPath=snapPath, snapRange=snapRange, snapStep=snapStep, radRange=radRanges[*,j])
      
      numTracers = n_elements(tp.idTargets)
      
      ; normalize times to multiples of the rotation period at R=1 (2pi)
      tp.times /= (2*!pi)
      
      ; normalize radii to their starting radii
      for i=0,numTracers-1 do $
        tp.trRad[*,i] /= tp.trRad[0,i]
    
      ; make plot
      xrange = [0,ceil(max(tp.times))]
      yrange = [0.8,1.2] ;[0.96,1.04] for 128.GC
    
      fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,$
           title="tracer radial evolution - "+ts+"."+res,$
           xtitle="t / P",ytitle="R / R"+textoidl("_0")
  
      fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot,thick=!p.thick-0.5
  
      ; loop over tracers and plot    
      for i=0L,numTracers-1 do $
        fsc_plot,tp.times,tp.trRad[*,i],line=0,thick=!p.thick-2.5,color=fsc_color('black'), /overplot
  
      ; write text for radRange
      legStr = string(radRanges[0,j],format='(f3.1)') + " < R"+textoidl("_0")+" < " + $
               string(radRanges[1,j],format='(f3.1)') + " ("+str(numTracers)+" tracers)"
      fsc_text,xrange[1]*0.3,yrange[0]+(yrange[1]-yrange[0])/10.0,$
               legStr,alignment=0.5,/data,charsize=!p.charsize-0.4
             
    endfor ;j

  end_PS, pngResize=50, /deletePS
  ;endif ;0
  
  ; plot 4 histograms of R/R0 at t/P=10,20,30,40
  start_PS, workingPath + plotBase + 'trDispHist_'+string(snapRange[0],format='(I3.3)')+'_'+$
            string(snapRange[1],format='(I3.3)')+'.eps', xs=7, ys=9
    
    !p.multi = [0,1,numRanges]
    
    for j=0,numRanges-1 do begin
      ; get tracer positions, radii, times
      tp = getTracerPos(snapPath=snapPath, snapRange=snapRange, snapStep=snapStep, radRange=radRanges[*,j])
      
      numTracers = n_elements(tp.idTargets)
      
      ; normalize times to multiples of the rotation period at R=1 (2pi)
      tp.times /= (2*!pi)
      
      ; normalize radii to their starting radii
      for i=0,numTracers-1 do $
        tp.trRad[*,i] /= tp.trRad[0,i]
        
      ; start plot
      yrange = [0,400]   ;[0,200] for 128.GC
      xrange = [0.9,1.1] ;[0.98,1.02] for 128.GC
      
      fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,$
           title="tracer radial evolution - "+ts+"."+res+" ("+string(radRanges[0,j],format='(f3.1)') + $
                 " < R"+textoidl("_0")+" < " + string(radRanges[1,j],format='(f3.1)')+")",$
           xtitle="R / R"+textoidl("_0"),ytitle="N"
           
      fsc_plot,[1.0,1.0],yrange,line=0,color=fsc_color('light gray'),/overplot,thick=!p.thick-0.5     
           
      ; loop over each time target an overplot histogram
      legendStrs = []
      legendColors = []
      
      for k=0,n_elements(timeTargets)-1 do begin
        ; find snapshot in time
        ind = where(abs(tp.times-timeTargets[k]) eq min(abs(tp.times-timeTargets[k])),count)
        if (count gt 1) then ind = ind[0]
        
        print,j,k,ind
        
        ; histogram radii
        rad = tp.trRad[ind,*]
        
        plothist,rad,/auto,/overplot,color=fsc_color(units.colors[k+1])
            
        legendStrs = [legendStrs,'t/P = '+string(timeTargets[k],format='(I2)')]
        legendColors = [legendColors,units.colors[k+1 mod 20]]    
      endfor ;k
      
      ; legend
      legend,legendStrs,textcolors=legendColors,box=0,margin=0.1,/left,charsize=!p.charsize-0.5
    endfor ;j
    
  end_PS, pngResize=50
  
  stop
end

