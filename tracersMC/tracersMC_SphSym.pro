; tracersMC_SphSym.pro
; dev for MC tracer particles (spherically symmetric setups)
; dnelson feb.2012

; plotTracerRadiiWithTime():

pro plotTracerRadiiWithTime

  ; config
  workingPath = '/n/home07/dnelson/dev.tracerMC/'
  snapPath    = workingPath + 'gasSphere.gasonly.1e4.L0.2.cool.nosg.f1/output/'
  plotBase    = 'gasSphere.gasonly.1e4.L0.2.cool.nosg.f1'
  
  snaps = indgen(n_elements(file_search(snapPath+'snap_*.hdf5')))
  ;snaps = indgen(10)
  
  ; tracer selection
  nTr = 10
  tracerIDBase = 1000000000LL
  
  seed = 4243L ;4243L
  tracerIDs = tracerIDBase + round(randomu(seed,nTr)*1e4*1)
  
  ; arrays
  trRad = fltarr(nTr,n_elements(snaps))
  times = fltarr(n_elements(snaps))
  
  ; loop over each snapshot
  foreach snap,snaps,k do begin

    print,'Running: snap ['+str(snap)+']'
    
    ; load
    h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
    
    pos_gas = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pos')
    id_gas  = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='ids')
    
    trid_tr  = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='tracerID')
    parid_tr = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='parentID')

    ; locate target tracerIDs, match parentIDs to gas IDs
    for i=0,n_elements(tracerIDs)-1 do begin
      w = where(trid_tr eq tracerIDs[i],count)
      if (count ne 1) then stop
      
      parentID = parid_tr[w[0]]
      w = where(id_gas eq parentID,count)
      if (count ne 1) then stop
      
      ; calculate radii and save
      parentPos = pos_gas[*,w[0]]

      trRad[i,k] = sqrt( (parentPos[0]-h.boxSize*0.5)^2.0 + $
                         (parentPos[1]-h.boxSize*0.5)^2.0 + $
                         (parentPos[2]-h.boxSize*0.5)^2.0 )
    endfor
    
    times[k] = h.time
  endforeach
  
  ; plot
  start_PS, workingPath + plotBase + '_'+str(min(snaps))+'-'+str(max(snaps))+'.eps'
  
    xrange = minmax(times)
    yrange = [0.0,max(trRad)*1.2]
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="Time [Gyr]",ytitle="Radii of Tracers [Kpc]"

    for i=0,nTr-1 do $
      fsc_plot,times,trRad[i,*],line=0,/overplot,color=getColor(i)
    
  end_PS
  
end

; gasSphereRadHistogramMC(): calculate radial profiles of gas/tracer (MC) properties

function gasSphereRadHistogramMC, snapPath=snapPath, snaps=snaps, TracerMCPerCell=TracerMCPerCell

  if n_elements(snapPath) eq 0 or n_elements(snaps) eq 0 or n_elements(TracerMCPerCell) eq 0 then stop

  nSnaps = n_elements(snaps)

  ; config (CUSTOM BINS)
  ;radMinMax = alog10([1.0,160.0])
  ;nbins = 7
  ;radBins = [0.0,5.0,12.0,25.0,45.0,75.0,100.0,160.0]
  ;midBins = [2.5,8.5,18.5,35.0,60.0,87.5,130.0]   
  
  radMinMax = alog10([0.1,1.2])
  nbins = 12
  radBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2]
  midBins = [0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15]
  
  ; arrays
  radDens    = fltarr(nSnaps,nbins)
  radDensTR  = fltarr(nSnaps,nbins)
  
  numGas     = fltarr(nSnaps,nbins)
  numTR      = fltarr(nSnaps,nbins)
  
  times = fltarr(nSnaps)
  
  ; equivalent tracer mass (exclude background gas cells)
  h = loadSnapshotHeader(snapPath,snapNum=0)
  mass_gas = loadSnapshotSubset(snapPath,snapNum=0,partType='gas',field='mass')
  w = where(mass_gas gt 0.01*mean(mass_gas))
  trMassConst = mean(mass_gas[w]) / TracerMCPerCell
  print,'trMassConst ',trMassConst
  
  boxCen = h.boxSize / 2.0
  
  ; load
  foreach snap,snaps,k do begin
    print,'snap: ',str(snap)
    
    h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
    
    pos    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pos')
    mass   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='mass')
    num_tr = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='numtr')

    ; calculate radii of particles
    rad = reform(sqrt((pos[0,*]-boxCen)*(pos[0,*]-boxCen) + $
                      (pos[1,*]-boxCen)*(pos[1,*]-boxCen) + $
                      (pos[2,*]-boxCen)*(pos[2,*]-boxCen)))

    ; do binning
    for j=0, nbins-1,1 do begin
      ; shell volume normalization
      vol = 4*!pi/3 * (radBins[j+1]^3.0 - radBins[j]^3.0)
      
      w = where(rad ge radBins[j] and rad lt radBins[j+1], count)
      
      if (count ne 0) then begin
        radDens[k,j] = total(mass[w]) / vol * 1e10
        numGas[k,j]  = count
        
        ; count child tracers
        count_tr = total(num_tr[w])

        radDensTR[k,j] = count_tr * trMassConst / vol * 1e10
        numTR[k,j]     = count_tr
      endif
      
      times[k] = h.time

    endfor
    
  endforeach ;snap 
  
  r = {radMinMax:radMinMax,radBins:radBins,midBins:midBins,nbins:nbins,$
       radDens:radDens,radDensTR:radDensTR,numGas:numGas,numTR:numTR,$
       times:times,trMassConst:trMassConst}
  return,r 
  
end

; gasSphereRadProfiles(): calculate radial profiles of density, temperature, velocity of gas

pro gasSphereRadProfiles

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
end

; make a frame of (xy) projected density also scatterplots of gas positions w/ tracer indicators (2x1)

pro col2SphFrame2x1

  ; config
  workingPath = '/n/home07/dnelson/dev.tracerMC/'
  snapPath    = workingPath + 'col2Sph.gasonly.1e4.norot.nocool.nosg.f100/output/'
  plotBase    = 'col2Sph.norot.nocool.nosg.f100'

  trMCPerCell = 100

  snaps = indgen(n_elements(file_search(snapPath+'snap_*.hdf5')))
  ;snaps = [10]
 
  foreach snap,snaps do begin

    print,'Running: snap ['+str(snap)+']'
    
    ; load
    h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
    
    pos_gas = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pos')
    num_tr  = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='numtr')
    
    ; plot extent
    zoomFac = 10.0
    xyzrange = [h.boxSize/2.0 - h.boxSize/2.0/zoomFac,h.boxSize/2.0 + h.boxSize/2.0/zoomFac]
    
    ; plot - fluid and tracer quantities and positions
    start_PS, workingPath + plotBase + '_2x1_'+string(snap,format='(I04)')+'.eps',xs=12.0,ys=6.0
    
      psym = 8 ; circles good for ~10k, use dots (3) for >=100k
      minMax = [-6.0,-3.0]
      
      ; left (projected density image)
      fsc_plot, [0],[0],/nodata,xrange=xyzrange,yrange=xyzrange,/xs,/ys,$
                position=[0.0,0.0,0.5,1.0],$
                xtickname=replicate(' ',10),ytickname=replicate(' ',10)
      
        plotDensityField,snapPath,snap,axes=0,/psOut,/overPlot,/log,minMax=minMax
      
      ; right (scatterplot all gas)
      fsc_plot, [0],[0],/nodata,xrange=xyzrange,yrange=xyzrange,/xs,/ys,$
                position=[0.5,0.0,1.0,1.0],/noerase,$
                xtickname=replicate(' ',10),ytickname=replicate(' ',10)
      
        fsc_plot, pos_gas[0,*],pos_gas[1,*],psym=psym,symsize=0.2,/overplot,color=fsc_color('black')

        ; mark gas with no tracers, less than starting, more than starting
        w = where(num_tr eq 0,count)
        if (count gt 0) then $
          fsc_plot, pos_gas[0,w],pos_gas[1,w],psym=psym,symsize=0.2,/overplot,color=fsc_color('red')
          
        w = where(num_tr eq trMCPerCell,count)
        if (count gt 0) then $
          fsc_plot, pos_gas[0,w],pos_gas[1,w],psym=psym,symsize=0.2,/overplot,color=fsc_color('orange')
          
        if (trMCPerCell gt 1) then begin
          w = where(num_tr gt 0 and num_tr lt trMCPerCell,count)
          if (count gt 0) then $
            fsc_plot, pos_gas[0,w],pos_gas[1,w],psym=psym,symsize=0.2,/overplot,color=fsc_color('blue')
        endif
          
        w = where(num_tr gt trMCPerCell,count)
        if (count gt 0) then $
          fsc_plot, pos_gas[0,w],pos_gas[1,w],psym=psym,symsize=0.2,/overplot,color=fsc_color('green')
        
      ; title
      fsc_text,0.91,0.05,"snap="+string(snap,format='(i03)')+" time="+string(h.time,format='(f5.2)'),$
        alignment=0.5,charsize=!p.charsize-0.3,/normal,color=fsc_color('orange')
    
    end_PS, pngResize=50, /deletePS
    
  endforeach

end
