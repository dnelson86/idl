; arepoSphSym.pro
; snapshot analysis for 3D spherically symmetric cases
; dnelson feb.2012

pro gasRadialProfiles

  ; paths
  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath    = workingPath + 'evrardOrig.10k/output/'
  plotBase    = 'evrardOrig.10k'
  
  nBins = 30 ;40
  
  yRangeRho    = [0.1,500.0]
  yRangeRatio  = [0.01,3.0]
  
  yRangeNum      = [1,1e5]
  yRangeRatioNum = [0.1,2.0]
  
  yRangeVel    = [-60,60]
  yRangeUtherm = [0.1,1000.0]  
  
  ; snapshot selection
  snaps = [0,1,2,3,4,5,6]
  snaps = [0,6,7,8,9,10,11]
  snaps = [0,11,12,13,14,15,16]
  snaps = [0,17,18,19,20,21,22]
  snaps = [0,23,24,25,26,27,28]
  snaps = [0,29,30,31,32,33,34]
  snaps = [0,35,36,37,38,39,40]
  nSnaps = n_elements(snaps)
  
  h = loadSnapshotHeader(snapPath,snapNum=snaps[0])
 
  radMinMax = alog10([0.1,h.boxSize/2.0])
  xyzCen    = [h.boxSize,h.boxSize,h.boxSize]/2.0

  ; arrays
  radDens_gas   = fltarr(nSnaps,nBins)
  radNum_gas    = fltarr(nSnaps,nBins)
  radVelR_gas   = fltarr(nSnaps,nBins)
  radUtherm_gas = fltarr(nSnaps,nBins)
  
  radNum_tr     = fltarr(nSnaps,nBins)
  radDens_tr    = fltarr(nSnaps,nBins)
  
  times = fltarr(nSnaps)

  ; tracer mass
  mass   = loadSnapshotSubset(snapPath,snapNum=0,partType='gas',field='mass')
  trMassConst = total(mass) / n_elements(mass)

  foreach snap,snaps,k do begin
  
    ; gas load and radii
    h = loadSnapshotHeader(snapPath,snapNum=snap)
    print,snap,h.nparttot[0],h.nparttot[3]
    
    pos    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pos')
    vel    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='vel')
    mass   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='mass')
    u      = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='u')
    
    pos_tr = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='pos')
    
    rad_vec = fltarr(3,n_elements(mass))
    rad_vec[0,*] = pos[0,*] - xyzCen[0]
    rad_vec[1,*] = pos[1,*] - xyzCen[1]
    rad_vec[2,*] = pos[2,*] - xyzCen[2]
    
    rad = reform(sqrt(rad_vec[0,*]^2.0 + rad_vec[1,*]^2.0 + rad_vec[2,*]^2.0))
    rad_tr = reform(sqrt((pos_tr[0,*]-xyzCen[0])^2.0 + $
                         (pos_tr[1,*]-xyzCen[1])^2.0 + $
                         (pos_tr[2,*]-xyzCen[2])^2.0))
                         
    vel_r = fltarr(3,n_elements(mass))
    vel_r[0,*] = vel[0,*] * rad_vec[0,*]
    vel_r[1,*] = vel[1,*] * rad_vec[1,*]
    vel_r[2,*] = vel[2,*] * rad_vec[2,*]
    vel_r /= rad ;norm
   
    radBins = [logspace(radMinMax[0],radMinMax[1],nBins+1)]
    midBins = [logspace(radMinMax[0],radMinMax[1],nBins+1,/mid)]

    ; do binning
    for j=0,nbins-1,1 do begin
      ; shell volume normalization
      vol = 4*!pi/3 * (radBins[j+1]^3.0 - radBins[j]^3.0)
      
      w    = where(rad ge radBins[j] and rad lt radBins[j+1], count)
      w_tr = where(rad_tr ge radBins[j] and rad_tr lt radBins[j+1], count_tr)
      
      if (count ne 0) then radDens_gas[k,j] = total(mass[w]) / vol
      if (count ne 0) then radNum_gas[k,j]  = count
      
      if (count ne 0) then radVelR_gas[k,j] = mean(vel_r[w])
      if (count ne 0) then radUtherm_gas[k,j] = mean(u[w])
      
      if (count_tr ne 0) then radNum_tr[k,j]  = count_tr
      if (count_tr ne 0) then radDens_tr[k,j] = trMassConst * count_tr / vol
    endfor
    
    print,total(radNum_gas[k,*]),total(radNum_tr[k,*])
    
    times[k] = h.time
  endforeach
  
  ; start plot
  start_PS, workingPath + plotBase + '_rad.'+str(min(snaps))+'-'+str(max(snaps))+'.eps', /big
  
  !p.multi = [0,2,3]
  !p.charsize -= 0.3

    ; legend
    strings = "t = "+string(times,format='(f4.2)')
    colors  = getColor(indgen(nSnaps),/name)

    xrange  = 10.0^radMinMax
  
    ; density
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yRangeRho,/ys,/xs, $
         xtitle="radius",ytitle=textoidl("density"),/xlog,/ylog
         
    foreach snap,snaps,k do begin
      fsc_plot,midBins,radDens_gas[k,*],line=0,color=getColor(k),/overplot
      fsc_plot,midBins,radDens_tr[k,*],line=1,color=getColor(k),/overplot
    endforeach
         
    ; tr/gas density ratio
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yRangeRatio,/ys,/xs, $
         xtitle="radius",ytitle=textoidl("\rho_{tr} / \rho_{gas}"),/xlog,/ylog
         
    foreach snap,snaps,k do $
      fsc_plot,midBins,radDens_tr[k,*]/radDens_gas[k,*],line=0,color=getColor(k),/overplot
         
    legend,strings,textcolors=colors,box=0,/bottom,/right,charsize=!p.charsize-0.3    
         
    ; number counts
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yRangeNum,/ys,/xs, $
         xtitle="radius",ytitle=textoidl("number counts"),/xlog,/ylog
         
    foreach snap,snaps,k do begin
      fsc_plot,midBins,radNum_gas[k,*],line=0,color=getColor(k),/overplot
      fsc_plot,midBins,radNum_tr[k,*],line=1,color=getColor(k),/overplot
    endforeach
    
    legend,strings,textcolors=colors,box=0,/bottom,/right,charsize=!p.charsize-0.3  
    
    ; tr/gas number ratio
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yRangeRatioNum,/ys,/xs, $
         xtitle="radius",ytitle=textoidl("N_{tr} / N_{gas}"),/xlog,/ylog
         
    foreach snap,snaps,k do $
      fsc_plot,midBins,radNum_tr[k,*]/radNum_gas[k,*],line=0,color=getColor(k),/overplot
         
    ; vel_r
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yRangeVel,/ys,/xs, $
         xtitle="radius",ytitle=textoidl("radial velocity"),/xlog
         
    foreach snap,snaps,k do $
      fsc_plot,midBins,radVelR_gas[k,*],line=0,color=getColor(k),/overplot
         
    ; utherm
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yRangeUtherm,/ys,/xs, $
         xtitle="radius",ytitle="utherm",/xlog,/ylog
         
    foreach snap,snaps,k do $
      fsc_plot,midBins,radUtherm_gas[k,*],line=0,color=getColor(k),/overplot
     
  
  !p.multi = 0
  end_PS
    
  stop
end