; tracersShocktube.pro
; dev for tracer particles (equil gas spheres and colliding gas spheres)
; dnelson feb.2012

; findAspectRatio():

function findAspectRatio, pos, r200, boxCen, gas_mass=gas_mass, targetN=targetN

  rf    = 10    ; fraction of r200 for cutting off-axis
  rfTar = 1.5   ; fraction of r200 to locate
  bs    = 10.0  ; kpc
  maxR  = 190.0 ; kpc
  
  ; DIRECTION ONE
  w = where(abs(pos[2,*]-boxCen) le r200/rf and abs(pos[1,*]-boxCen) le r200/rf,count)
  
  rad = reform(sqrt((pos[0,w]-boxCen)*(pos[0,w]-boxCen) + $
                    (pos[1,w]-boxCen)*(pos[1,w]-boxCen) + $
                    (pos[2,w]-boxCen)*(pos[2,w]-boxCen)))
  
  if keyword_set(gas_mass) then $
    h = hist1d(rad,gas_mass[w],binsize=bs,obin=loc,binedge=-1) ;weight by mass
  if not keyword_set(gas_mass) then $
    h = histogram(rad,binsize=bs,locations=loc) ;number density only
    
  loc += bs/2.0 
  
  rpts = findgen(1001)/1000 * maxR
  nInterp = interpol(h,loc,rpts)
        
  if (targetN eq 0) then begin
    ; locate target number density
    w = where(abs(rpts-r200/rfTar) eq min(abs(rpts-r200/rfTar)),count)
    if (count gt 1) then w = w[0]
    targetN = nInterp[w[0]]
    print,'new target ',targetN
  endif
  
  ; save radius of target number density
  w = where(abs(nInterp-targetN) eq min(abs(nInterp-targetN)),count)
  if (count gt 1) then w = w[0]
  sizeX = rpts[w[0]]

  ; DIRECTION TWO
  w = where(abs(pos[0,*]-boxCen) le r200/rf and abs(pos[1,*]-boxCen) le r200/rf,count)
  
  rad = reform(sqrt((pos[0,w]-boxCen)*(pos[0,w]-boxCen) + $
                    (pos[1,w]-boxCen)*(pos[1,w]-boxCen) + $
                    (pos[2,w]-boxCen)*(pos[2,w]-boxCen)))
 
  if keyword_set(gas_mass) then $
    h = hist1d(rad,gas_mass[w],binsize=bs,obin=loc,binedge=-1) ;weight by mass
  if not keyword_set(gas_mass) then $
    h = histogram(rad,binsize=bs,locations=loc) ;number density only
    
  loc += bs/2.0 
  
  rpts = findgen(1001)/1000 * maxR
  nInterp = interpol(h,loc,rpts)
        
  ; locate target number density and save radius
  w = where(abs(nInterp-targetN) eq min(abs(nInterp-targetN)),count)
  if (count gt 1) then w = w[0]
  sizeZ = rpts[w[0]]
  
  aspectRatio = sizeZ / sizeX
  ;stop
  return, aspectRatio
 
end

; gasSphereRadHistogram(): calculate radial profiles of gas/tracer properties

function gasSphereRadHistogram, snapPath=snapPath, snaps=snaps

  if n_elements(snapPath) eq 0 or n_elements(snaps) eq 0 then stop

  nSnaps = n_elements(snaps)

  ; config
  ;nbins = 30
  
  ;radMinMax = alog10([1.0,300.0])
  radMinMax = alog10([1.0,500.0])
  
  ;CUSTOM BINS
  nbins = 10
  radBins = [0.0,5.0,12.0,25.0,45.0,75.0,100.0,130.0,180.0,300.0,500.0]
  midBins = [2.5,8.5,18.5,35.0,60.0,87.5,115.0,155.0,240.0,400.0]   
  
  ; arrays
  radDens    = fltarr(nSnaps,nbins)
  radDensTR  = fltarr(nSnaps,nbins)
  
  numGas     = fltarr(nSnaps,nbins)
  numTR      = fltarr(nSnaps,nbins)
  
  times = fltarr(nSnaps)
  
  ; aspect ratio arrays
  targetNgas = 0.0
  targetNtr  = 0.0
  
  gasAR   = fltarr(nSnaps)
  trAR    = fltarr(nSnaps)
  
  ; radial bins
  ;radBins = 10.0^( findgen(nbins+1)/nbins*(radMinMax[1]-radMinMax[0]) + radMinMax[0] )
  ;midBins = 10.0^( (findgen(nbins)+0.5)/nbins*(radMinMax[1]-radMinMax[0]) + radMinMax[0] )
  
  ; tracer mass (t=0)
  h    = loadSnapshotHeader(snapPath,snapNum=0,/verbose)
  mass = loadSnapshotSubset(snapPath,snapNum=0,partType='gas',field='mass')
  trMassConst = total(mass) / h.nPartTot[3]
  ;trMassConst = 0.0001

  boxCen = h.boxSize / 2.0
  
  ; load
  foreach snap,snaps,k do begin
    print,'snap: ',str(snap)
    
    h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
    
    pos    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pos')
    mass   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='mass')
    ;dens   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='rho')

    ; calculate radii of particles
    rad = reform(sqrt((pos[0,*]-boxCen)*(pos[0,*]-boxCen) + $
                      (pos[1,*]-boxCen)*(pos[1,*]-boxCen) + $
                      (pos[2,*]-boxCen)*(pos[2,*]-boxCen)))

    ; load tracer positions and calculate radii
    pos_tr = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='pos')
    
    rad_tr = reform(sqrt((pos_tr[0,*]-boxCen)*(pos_tr[0,*]-boxCen) + $
                         (pos_tr[1,*]-boxCen)*(pos_tr[1,*]-boxCen) + $
                         (pos_tr[2,*]-boxCen)*(pos_tr[2,*]-boxCen)))

    ; calculate aspect ratios
    ;gasAR[k] = findAspectRatio(pos,r200,boxCen,targetN=targetNgas,gas_mass=mass)
    ;trAR[k]  = findAspectRatio(pos_tr,r200,boxCen,targetN=targetNtr)
               
    ; do binning
    for j=0, nbins-1,1 do begin
      ; shell volume normalization
      vol = 4*!pi/3 * (radBins[j+1]^3.0 - radBins[j]^3.0)
      
      w = where(rad ge radBins[j] and rad lt radBins[j+1], count)
      w_tr = where(rad_tr ge radBins[j] and rad_tr lt radBins[j+1], count_tr)
      
      if (count ne 0) then begin
        radDens[k,j] = total(mass[w]) / vol * 1e10
        numGas[k,j]  = count
      endif
      
      if (count_tr ne 0) then begin
        radDensTR[k,j] = count_tr * trMassConst / vol * 1e10
        numTR[k,j]     = count_tr
      endif
      
      times[k] = h.time

    endfor
    
  endforeach ;snap 
  
  r = {radMinMax:radMinMax,radBins:radBins,midBins:midBins,nbins:nbins,$
       radDens:radDens,radDensTR:radDensTR,numGas:numGas,numTR:numTR,$
       times:times,gasAR:gasAR,trAR:trAR,trMassConst:trMassConst}
  return,r 
  
end

; gasSphereRadProfiles(): calculate radial profiles of density, temperature, velocity of gas

pro gasSphereRadProfiles

  ;r_s  = 327.8 ;kpc (hernquist)
  r200 = 162.6 ;kpc (nfw)
  
  ; paths
  workingPath = '/n/home07/dnelson/dev.tracerMC/'
  ;snapPath    = workingPath + 'gasSphere.gasonly.1e4.norot.nocool.nosg/output/'
  ;plotBase    = 'gasSphere.gasonly.1e4.norot.nocool.nosg'
  
  ;snapPath    = workingPath + 'col2Sph.gasonly.1e4.norot.nocool.nosg.f100/output/'
  ;plotBase    = 'col2Sph.gasonly.1e4.norot.nocool.nosg.f100'
  
  snapPath = workingPath + 'evrard.gasonly.3D.1e4.f1/output/'
  plotBase = 'evrard.gasonly.3D.1e4.f1'
  
  TracerMCPerCell = 1
  
  ;snaps =  [0,1,2,3,4,5,10,15,20,25,50,75,100,125,150,175,200] ;indgen(200)
  snaps = indgen(18)
  
  ; calculate radial/box profiles
  ;gsp = gasSphereRadHistogram(snapPath=snapPath,snaps=snaps)
  gsp = gasSphereRadHistogramMC(snapPath=snapPath,snaps=snaps,TracerMCPerCell=TracerMCPerCell)

  ; plots
  start_PS, workingPath + plotBase + '_radDens_'+string(min(snaps),format='(I3.3)')+'_'+$
            string(max(snaps),format='(I3.3)')+'.eps'
  
    xrange = 10.0^gsp.radMinMax
    yrange = [min(gsp.radDens[where(gsp.radDens ne 0)])*0.5,max([gsp.radDens,gsp.radDensTR])*2.0]
    
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="",ytitle=textoidl("\rho_{gas,tr}(r)"),$ ; [M_{sun} kpc^{-3}]
              title=plotBase,/ylog,$
              position=[0.18,0.35,0.9,0.9],xtickname=replicate(' ',10)
    
    ; density profiles for successive snapshots   
    foreach snap,snaps,k do begin
        fsc_plot,gsp.midBins,gsp.radDens[k,*],  line=0,/overplot,color=getColor(k);,thick=!p.thick+0.5
        fsc_plot,gsp.midBins,gsp.radDensTR[k,*],line=1,/overplot,color=getColor(k);,thick=!p.thick+0.5 
    endforeach
    
    ; snapshot legend
    strings = 't = '+string(gsp.times,format='(f4.2)')+' Gyr'
    colors  = getColor(indgen(n_elements(snaps)),/name)
    legend,strings,textcolors=colors,box=0,/top,/right,charsize=!p.charsize-0.6
    
    legend,['gas','tracer'],linestyle=[0,1],textcolors=['black','black'],$
           box=0,/left,/top,charsize=!p.charsize-0.4
           
    ; residual plot
    yrange = [0.0,2.0]
    fsc_plot,[0],[0],/nodata,/noerase,xrange=xrange,yrange=yrange,/xs,/ys,$
             xtitle="radius [kpc]",$
             ytitle=textoidl("\rho_{tr} / \rho_{gas}"),position=[0.18,0.15,0.9,0.35],$
             ytickv=[0.5,1.0,1.5],yticks=2
             
    fsc_plot,xrange,[0.5,0.5],line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[1.5,1.5],line=1,color=fsc_color('light gray'),/overplot
    
    foreach snap,snaps,k do begin
      fsc_plot,gsp.midBins,gsp.radDensTR[k,*]/gsp.radDens[k,*],$
        psym=8,symsize=0.7,color=getColor(k),/overplot
    endforeach

  end_PS

  ; plot two - overdensity with time for different midbins
  start_PS, workingPath + plotBase + '_overDens_'+string(min(snaps),format='(I3.3)')+'_'+$
            string(max(snaps),format='(I3.3)')+'.eps'
  
    ratio = gsp.radDensTR/gsp.radDens
    ratio = ratio[where(finite(ratio))]
    yrange = [min(ratio)*0.5,max(ratio)*1.05]>[0.1,1.2]<[0.8,4.0]
    ;yrange = [0.95,1.05]
  
    fsc_plot, [0],[0],/nodata,xrange=minmax(gsp.times),yrange=yrange,/xs,/ys,$
              xtitle="time [Gyr]",ytitle=textoidl("\rho_{tr} / \rho_{gas}"),$ ; [M_{sun} kpc^{-3}]
              title=plotBase
  
    fsc_plot, minmax(gsp.times), [1.0,1.0], line=0, color=fsc_color('light gray'),/overplot
  
    ;bins = [0,1,2,3,5] ;gasSphere or col2Sph
    bins = [0,2,4,6,8,10] ;evrard
    strings = []
    
    foreach bin,bins,k do begin
      ratio = gsp.radDensTR[*,bin] / gsp.radDens[*,bin]
      
      fsc_plot,gsp.times,ratio,line=0,color=getColor(k),/overplot
      strings = [strings,'<r> = '+string(gsp.midBins[bin],format='(f5.1)')+' kpc']
    endforeach
    
    ; legend
    legend,strings,textcolors=getColor(indgen(n_elements(bins)),/name),$
      box=0,/right,/top,charsize=!p.charsize-0.3
    
  end_PS
  stop
  ; plot three - individual snapshots
  for k=0,n_elements(snaps)-1 do begin
    print,gsp.times[k]
  start_PS, workingPath + plotBase + 'radDens_'+str(k)+'_'+string(min(snaps),format='(I3.3)')+'_'+$
            string(max(snaps),format='(I3.3)')+'.eps'
            
    xrange = 10.0^gsp.radMinMax
    yrange = [min(gsp.radDens[where(gsp.radDens ne 0)])*0.5,max([gsp.radDens,gsp.radDensTR])*2.0]
    
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="",ytitle=textoidl("\rho_{gas,tr}(r)"),$ ; [M_{sun} kpc^{-3}]
              title=plotBase,/ylog,$
              position=[0.18,0.35,0.9,0.9],xtickname=replicate(' ',10)
    
    ; density profile
    fsc_plot,gsp.midBins,gsp.radDens[k,*],  line=0,/overplot,color=getColor(k);,thick=!p.thick+0.5
    fsc_plot,gsp.midBins,gsp.radDensTR[k,*],line=1,/overplot,color=getColor(k);,thick=!p.thick+0.5 

    ; snapshot legend
    strings = 't = '+string(gsp.times,format='(f4.2)')+' Gyr'
    colors  = getColor(indgen(n_elements(snaps)),/name)
    legend,strings,textcolors=colors,box=0,/top,/right,charsize=!p.charsize-0.6
    
    legend,['gas','tracer'],linestyle=[0,1],textcolors=['black','black'],$
           box=0,/left,/top,charsize=!p.charsize-0.4
           
    ; residual plot
    yrange = [0.0,2.0]
    fsc_plot,[0],[0],/nodata,/noerase,xrange=xrange,yrange=yrange,/xs,/ys,$
             xtitle="radius [kpc]",$
             ytitle=textoidl("\rho_{tr} / \rho_{gas}"),position=[0.18,0.15,0.9,0.35],$
             ytickv=[0.5,1.0,1.5],yticks=2
             
    fsc_plot,xrange,[0.5,0.5],line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[1.5,1.5],line=1,color=fsc_color('light gray'),/overplot
    
    fsc_plot,gsp.midBins,gsp.radDensTR[k,*]/gsp.radDens[k,*],$
      psym=8,symsize=0.7,color=getColor(k),/overplot
            
  end_PS, pngResize=50, /deletePS
  endfor ;k
stop
end

; gasSphereBoxHistogram(): calculate box profiles of gas/tracer properties

function gasSphereBoxHistogram, snapPath=snapPath, snaps=snaps

  if n_elements(snapPath) eq 0 or n_elements(snaps) eq 0 then stop

  nSnaps = n_elements(snaps)

  ; config
  nbins = 20
  
  boxSize = 5.0 ;kpc
  
  ; binning (along x axis)
  radMinMax = alog10([0.00001,(nbins-1)*boxSize])
  midBins = linspace(0.0,(nbins-1)*boxSize,nbins)
  
  axes = [0,1,2] ;move along first, always within boxSize of second and third
  
  ; arrays
  radDens    = fltarr(nSnaps,nbins)
  radDensTR  = fltarr(nSnaps,nbins)
  
  numGas     = fltarr(nSnaps,nbins)
  numTR      = fltarr(nSnaps,nbins)
  
  times = fltarr(nSnaps)
  
  ; tracer mass (t=0)
  h    = loadSnapshotHeader(snapPath,snapNum=0,/verbose)
  mass = loadSnapshotSubset(snapPath,snapNum=0,partType='gas',field='mass')
  trMassConst = total(mass) / h.nPartTot[3]
  ;trMassConst = 0.0001

  boxCen = h.boxSize / 2.0
  
  ; load
  foreach snap,snaps,k do begin
    print,'snap: ',str(snap)
    
    h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
    
    pos    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pos')
    mass   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='mass')
    
    ; load tracer positions and calculate radii
    pos_tr = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='pos')

    ; shift positions relative to big box center
    pos    -= boxCen
    pos_tr -= boxCen

    ; pre-select on axes[1,2]
    w = where(abs(pos[axes[1],*]) le boxSize*0.5 and abs(pos[axes[2],*]-0.0) le boxSize*0.5, count)
    pos = pos[axes[0],w]
    mass = mass[w]
    
    w = where(abs(pos_tr[axes[1],*]) le boxSize*0.5 and abs(pos_tr[axes[2],*]-0.0) le boxSize*0.5, count_tr)
    pos_tr = pos_tr[axes[0],w]

    ; do binning
    for j=0, nbins-1,1 do begin
      ; box volume normalization
      vol = boxSize^3.0
      
      ; subselect within this box
      w = where( abs(pos-midBins[j]) le boxSize*0.5,count)
                 
      w = where( abs(pos_tr-midBins[j]) le boxSize*0.5,count_tr)
      
      if (count ne 0) then begin
        radDens[k,j] = total(mass[w]) / vol * 1e10
        numGas[k,j]  = count
      endif
      
      if (count_tr ne 0) then begin
        radDensTR[k,j] = count_tr * trMassConst / vol * 1e10
        numTR[k,j]     = count_tr
      endif
      
      times[k] = h.time

    endfor
    
  endforeach ;snap 

  r = {radMinMax:radMinMax,midBins:midBins,nbins:nbins,$
       radDens:radDens,radDensTR:radDensTR,numGas:numGas,numTR:numTR,$
       times:times,trMassConst:trMassConst,boxSize:boxSize}
  return,r 
  
end

; gasSphereBoxProfiles(): calculate box-radial profiles of density, temperature, velocity of gas

pro gasSphereBoxProfiles

  ; paths
  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath    = workingPath + 'col2Sph.gastr.1e5.f1/output/'
  plotBase    = 'col2Sph.1e5.f1_box'
  
  snaps = indgen(201)
  
  ; calculate radial/box profiles
  gsp = gasSphereBoxHistogram(snapPath=snapPath,snaps=snaps)

  ; plots
  start_PS, workingPath + plotBase + '_radDens_'+string(min(snaps),format='(I3.3)')+'_'+$
            string(max(snaps),format='(I3.3)')+'.eps'
  
    xrange = 10.0^gsp.radMinMax
    yrange = [min(gsp.radDens[where(gsp.radDens ne 0)])*0.5,max([gsp.radDens,gsp.radDensTR])*2.0]
    ;yrange = [1e4,5e7]
    
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="",ytitle=textoidl("\rho_{gas,tr}(r)"),$ ; [M_{sun} kpc^{-3}]
              title=plotBase,/ylog,$
              position=[0.18,0.35,0.9,0.9],xtickname=replicate(' ',10)
    
    ; density profiles for successive snapshots   
    foreach snap,snaps,k do begin
        fsc_plot,gsp.midBins,gsp.radDens[k,*],  line=0,/overplot,color=getColor(k);,thick=!p.thick+0.5
        fsc_plot,gsp.midBins,gsp.radDensTR[k,*],line=1,/overplot,color=getColor(k);,thick=!p.thick+0.5 
    endforeach
    
    ; snapshot legend
    strings = 't = '+string(gsp.times,format='(f4.2)')+' Gyr'
    colors  = getColor(indgen(n_elements(snaps)),/name)
    legend,strings,textcolors=colors,box=0,/top,/right,charsize=!p.charsize-0.6
    
    legend,['gas','tracer'],linestyle=[0,1],textcolors=['black','black'],$
           box=0,/left,/top,charsize=!p.charsize-0.4
           
    ; residual plot
    yrange = [0.0,2.0]
    fsc_plot,[0],[0],/nodata,/noerase,xrange=xrange,yrange=yrange,/xs,/ys,$
             xtitle="radius [kpc]",$
             ytitle=textoidl("\rho_{tr} / \rho_{gas}"),position=[0.18,0.15,0.9,0.35],$
             ytickv=[0.5,1.0,1.5],yticks=2
             
    fsc_plot,xrange,[0.5,0.5],line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[1.5,1.5],line=1,color=fsc_color('light gray'),/overplot
    
    foreach snap,snaps,k do begin
      fsc_plot,gsp.midBins,gsp.radDensTR[k,*]/gsp.radDens[k,*],$
        psym=8,symsize=0.7,color=getColor(k),/overplot
    endforeach

  end_PS

  ; plot two - overdensity with time for different midbins
  start_PS, workingPath + plotBase + '_overDens_'+string(min(snaps),format='(I3.3)')+'_'+$
            string(max(snaps),format='(I3.3)')+'.eps'
  
    ;yrange = [min(radDensTR/radDens),max(radDensTR/radDens)]
    yrange = [0.5,2.0]
  
    fsc_plot, [0],[0],/nodata,xrange=minmax(gsp.times),yrange=yrange,/xs,/ys,$
              xtitle="time [Gyr]",ytitle=textoidl("\rho_{tr} / \rho_{gas}"),$ ; [M_{sun} kpc^{-3}]
              title=plotBase
  
    fsc_plot, minmax(gsp.times), [1.0,1.0], line=0, color=fsc_color('light gray'),/overplot
  
    bins = [0,1,2,3,5,7]
    strings = []
    
    foreach bin,bins,k do begin
      ratio = gsp.radDensTR[*,bin] / gsp.radDens[*,bin]
      
      fsc_plot,gsp.times,ratio,line=0,color=getColor(k),/overplot
      strings = [strings,'<r> = '+string(gsp.midBins[bin]+gsp.boxSize/2.0,format='(f5.1)')+' kpc']
    endforeach
    
    ; legend
    legend,strings,textcolors=getColor(indgen(n_elements(bins)),/name),$
      box=0,/right,/top,charsize=!p.charsize-0.3
    
  end_PS
  
end

; gasSphereRunsProf(): overplot tracer density profiles and ratios for different runs at one time

pro gasSphereRunsProf

  ; config
  nbins = 20
  r200  = 162.6 ;kpc
  radMinMax = alog10([1.0,200.0])

  ; paths
  workingPath = '/n/home07/dnelson/dev.tracer/'
  plotBase    = 'gasSphere.gastr.res'
  
  ; snapshot selection
  snap = 4
  ;runs = '2e5.f'+['1','4','10','0']
  runs = 'res.'+['2e5','1e5','5e4','1e4','5e3']
  
  ; arrays
  nRuns = n_elements(runs)
  
  radDensGas = fltarr(nRuns,nbins)
  radDensTR  = fltarr(nRuns,nbins)
  
  ; radial bins
  radBins = [0.0,           logspace(radMinMax[0],radMinMax[1],nbins)]
  midBins = [radBins[0]/2.0,logspace(radMinMax[0],radMinMax[1],nbins,/mid)]
  
  ; load
  foreach run,runs,k do begin
    print,'run: ',str(run)
    
    snapPath    = workingPath + 'gasSphere.gastr.'+run+'/output/'
    
    ; tracer mass (t=0)
    h    = loadSnapshotHeader(snapPath,snapNum=0,/verbose)
    mass = loadSnapshotSubset(snapPath,snapNum=0,partType='gas',field='mass')
    
    trMassConst = total(mass) / h.nPartTot[3]
    
    ; load
    h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
    
    boxCen = h.boxSize / 2.0
    
    pos    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pos')
    mass   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='mass')

    ; calculate radii of particles
    rad_gas = reform(sqrt((pos[0,*]-boxCen)*(pos[0,*]-boxCen) + $
                          (pos[1,*]-boxCen)*(pos[1,*]-boxCen) + $
                          (pos[2,*]-boxCen)*(pos[2,*]-boxCen)))

    ; load tracer positions and calculate radii
    pos_tr = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='pos')
    
    rad_tr = reform(sqrt((pos_tr[0,*]-boxCen)*(pos_tr[0,*]-boxCen) + $
                         (pos_tr[1,*]-boxCen)*(pos_tr[1,*]-boxCen) + $
                         (pos_tr[2,*]-boxCen)*(pos_tr[2,*]-boxCen)))

    ; do binning
    for j=0, nbins-1 do begin
      ; shell volume normalization
      vol = 4*!pi/3 * (radBins[j+1]^3.0 - radBins[j]^3.0)
      
      w    = where(rad_gas ge radBins[j] and rad_gas lt radBins[j+1], count_gas)
      w_tr = where(rad_tr  ge radBins[j] and rad_tr  lt radBins[j+1], count_tr)
      
      if (count_gas ne 0) then $
        radDensGas[k,j] = total(mass[w]) / vol * 1e10
      if (count_tr ne 0) then $
        radDensTR[k,j]  = count_tr * trMassConst / vol * 1e10

    endfor
    
  endforeach ;snap

  ; plots
  start_PS, workingPath + plotBase + '_snap='+str(snap)+'.radDens.eps'
  
    xrange = 10.0^radMinMax
    yrange = [1e7,max([radDensTR])*1.5]
    yrange = [1e7,1e8]
    
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="",ytitle=textoidl("\rho_{tr}(r) * R^2 [M_{sun} kpc^{-1}]"),$
              title=plotBase + textoidl(" \Lambda=0.1")+" t = "+string(h.time,format='(f3.1)')+" Gyr",$
              /ylog,/xlog,position=[0.18,0.35,0.9,0.9],xtickname=replicate(' ',10)
    
    ; density profiles for successive snapshots  
    fsc_plot,[r200,r200],yrange,line=2,color=fsc_color('light gray'),/overplot   
    fsc_text,r200*0.75,yrange[0]*2,textoidl("r_{200}"),alignment=0.0,color=fsc_color('light gray'),$
             charsize=!p.charsize-0.4
   
    ; gas profile only for first run
    fsc_plot,midBins,radDensGas[0,*]*midBins^2.0,line=1,/overplot,color=fsc_color('black')
    
    foreach run,runs,k do begin
      fsc_plot,midBins,radDensTR[k,*]*midBins^2.0,line=0,/overplot,color=getColor(k)
    endforeach
    
    ; snapshot legend
    ;strings = ['1=','4=','10=','1/4='] + textoidl('N_{tr} / N_{gas}')
    strings = ['gas 2e5',runs] ;+ textoidl('N_{gas}')
    lines = [1,intarr(Nruns)]
    colors  = getColor([0,indgen(nRuns)],/name)
    
    legend,strings,textcolors=colors,linestyle=lines,$
      box=0,/top,/right,charsize=!p.charsize-0.3,linesize=0.3
           
    ; residual plot
    yrange = [0.6,1.4]
    fsc_plot,[0],[0],/nodata,/noerase,xrange=xrange,yrange=yrange,/xs,/ys,/xlog,$
             xtitle="radius [kpc]",$
             ytitle=textoidl("\rho_{tr,res} / \rho_{tr,2e5}"),position=[0.18,0.15,0.9,0.35],$
             ytickv=[0.8,1.0,1.2],yticks=2
  
    fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot    
    fsc_plot,xrange,[0.8,0.8],line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[1.2,1.2],line=1,color=fsc_color('light gray'),/overplot
    
    fsc_plot,[r200,r200],yrange,line=2,color=fsc_color('light gray'),/overplot
    
    ; just interpolate both onto a set of radii
    nbins = 100
    res_pts = 10.0^( findgen(nbins+1)/nbins * (alog10(r200)-alog10(1)) + alog10(1) )
    
    tr_res_f1 = interpol(radDensTR[0,*],midBins,res_pts)
    
    for k=1,n_elements(runs)-1 do begin
      ;gas_res = interpol(radDensGas[k,*],midBins,res_pts)
      tr_res  = interpol(radDensTR[k,*],midBins,res_pts)
  
      fsc_plot,res_pts,tr_res/tr_res_f1,line=0,color=getColor(k),/overplot
    endfor
    
  end_PS
  stop
end

; make a frame of (xy) (yz) (xz) projected density also scatterplots of gas and tracer positions (3x3)

pro col2SphFrame3x3

  ; config
  workingPath = '/n/home07/dnelson/dev.tracerMC/'
  snapPath    = workingPath + 'col2Sph.gasonly.1e4.norot.nocool.nosg.f1/output/'
  plotBase    = 'col2Sph.gasonly.1e4.norot.nocool.nosg.f1'

  snaps = indgen(n_elements(file_search(snapPath+'snap_*.hdf5')))
  
  foreach snap,snaps do begin

    print,'Running: snap ['+str(snap)+']'
    
    ; load
    h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
    
    pos_gas = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pos')
    pos_tr  = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='pos')
    
    ; plot extent
    zoomFac = 10.0
    xyzrange = [h.boxSize/2.0 - h.boxSize/2.0/zoomFac,h.boxSize/2.0 + h.boxSize/2.0/zoomFac]
    
    ; plot - fluid and tracer quantities and positions
    start_PS, workingPath + plotBase + '_'+string(snap,format='(I04)')+'.eps',xs=7.0,ys=7.0
    
      psym = 3
      minMax = [-6.0,-3.0]
      
      ; top row (projected density images)
      fsc_plot, [0],[0],/nodata,xrange=xyzrange,yrange=xyzrange,/xs,/ys,$
                position=[0.05,0.65,0.35,0.95],$
                xtickname=replicate(' ',10),ytickname=replicate(' ',10) ;UL
      
        plotDensityField,snapPath,snap,axes=0,/psOut,/overPlot,/log,minMax=minMax
      
      fsc_plot, [0],[0],/nodata,xrange=xyzrange,yrange=xyzrange,/xs,/ys,$
                position=[0.35,0.65,0.65,0.95],/noerase,$
                xtickname=replicate(' ',10),ytickname=replicate(' ',10) ;UC
      
        plotDensityField,snapPath,snap,axes=1,/psOut,/overPlot,/log,minMax=minMax
  
      fsc_plot, [0],[0],/nodata,xrange=xyzrange,yrange=xyzrange,/xs,/ys,$
                position=[0.65,0.65,0.95,0.95],/noerase,$
                xtickname=replicate(' ',10),ytickname=replicate(' ',10) ;UR
      
        plotDensityField,snapPath,snap,axes=2,/psOut,/overPlot,/log,minMax=minMax
  
      ; middle row (gas position scatterplots)
      fsc_plot, [0],[0],/nodata,xrange=xyzrange,yrange=xyzrange,/xs,/ys,$
                position=[0.05,0.35,0.35,0.65],/noerase,$
                xtickname=replicate(' ',10),ytickname=replicate(' ',10) ;UL
      
        fsc_plot, pos_gas[0,*],pos_gas[1,*],psym=psym,/overplot,color=getColor(0)
      
      fsc_plot, [0],[0],/nodata,xrange=xyzrange,yrange=xyzrange,/xs,/ys,$
                position=[0.35,0.35,0.65,0.65],/noerase,$
                xtickname=replicate(' ',10),ytickname=replicate(' ',10) ;UC
      
        fsc_plot, pos_gas[0,*],pos_gas[2,*],psym=psym,/overplot,color=getColor(0)
  
      fsc_plot, [0],[0],/nodata,xrange=xyzrange,yrange=xyzrange,/xs,/ys,$
                position=[0.65,0.35,0.95,0.65],/noerase,$
                xtickname=replicate(' ',10),ytickname=replicate(' ',10) ;UR
      
        fsc_plot, pos_gas[1,*],pos_gas[2,*],psym=psym,/overplot,color=getColor(0)
      
      ; bottom row (tracer position scatterplots)
      fsc_plot, [0],[0],/nodata,xrange=xyzrange,yrange=xyzrange,/xs,/ys,$
                position=[0.05,0.05,0.35,0.35],/noerase,$
                xtickname=replicate(' ',10),ytickname=replicate(' ',10) ;LL
      
        fsc_plot, pos_tr[0,*],pos_tr[1,*],psym=psym,/overplot,color=getColor(0)
      
      fsc_plot, [0],[0],/nodata,xrange=xyzrange,yrange=xyzrange,/xs,/ys,$
                position=[0.35,0.05,0.65,0.35],/noerase,$
                xtickname=replicate(' ',10),ytickname=replicate(' ',10) ;LL
      
        fsc_plot, pos_tr[0,*],pos_tr[2,*],psym=psym,/overplot,color=getColor(0)
  
      fsc_plot, [0],[0],/nodata,xrange=xyzrange,yrange=xyzrange,/xs,/ys,$
                position=[0.65,0.05,0.95,0.35],/noerase,$
                xtickname=replicate(' ',10),ytickname=replicate(' ',10) ;LR
      
        fsc_plot, pos_tr[1,*],pos_tr[2,*],psym=psym,/overplot,color=getColor(0)
      
      ; labels
      fsc_text,0.07,0.07,"tracers",charsize=!p.charsize-0.2,color=getColor(1),/normal
      fsc_text,0.07,0.37,"gas",charsize=!p.charsize-0.2,color=getColor(2),/normal
      
      ; title
      fsc_text,0.5,0.96,"snap="+str(snap)+" time="+string(h.time,format='(f5.2)'),alignment=0.5,$
               charsize=!p.charsize-0.2,/normal
    
    end_PS, pngResize=50, /deletePS
    
  endforeach

end

; make a frame of (xy) projected density also scatterplots of gas and tracer positions (2x1)

pro col2SphFrame2x1

  ; config
  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath    = workingPath + 'col2Sph.gastr.1e5.f1/output/'
  plotBase    = 'col2Sph'

  snaps = [200]
  
  foreach snap,snaps do begin

    print,'Running: snap ['+str(snap)+']'
    
    ; load
    h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
    
    pos_gas = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pos')
    pos_tr  = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='pos')
    
    ; plot extent
    zoomFac = 10.0
    xyzrange = [h.boxSize/2.0 - h.boxSize/2.0/zoomFac,h.boxSize/2.0 + h.boxSize/2.0/zoomFac]
    
    ; plot - fluid and tracer quantities and positions
    start_PS, workingPath + plotBase + '_2x1_'+string(snap,format='(I04)')+'.eps',xs=12.0,ys=6.0
    
      psym = 3
      minMax = [-6.0,-3.0]
      
      ; left (projected density image)
      fsc_plot, [0],[0],/nodata,xrange=xyzrange,yrange=xyzrange,/xs,/ys,$
                position=[0.0,0.0,0.5,1.0],$
                xtickname=replicate(' ',10),ytickname=replicate(' ',10)
      
        plotDensityField,snapPath,snap,axes=0,/psOut,/overPlot,/log,minMax=minMax
      
      ; right (scatterplot gas over tracers)
      fsc_plot, [0],[0],/nodata,xrange=xyzrange,yrange=xyzrange,/xs,/ys,$
                position=[0.5,0.0,1.0,1.0],/noerase,$
                xtickname=replicate(' ',10),ytickname=replicate(' ',10)
      
        fsc_plot, pos_tr[0,*], pos_tr[1,*], psym=psym,/overplot,color=getColor(1)
        fsc_plot, pos_gas[0,*],pos_gas[1,*],psym=psym,/overplot,color=getColor(0)

        
      ; title
      fsc_text,0.91,0.05,"snap="+string(snap,format='(i03)')+" time="+string(h.time,format='(f5.2)'),$
        alignment=0.5,charsize=!p.charsize-0.3,/normal,color=fsc_color('orange')
    
    end_PS, pngResize=50, /deletePS
    
  endforeach

end
