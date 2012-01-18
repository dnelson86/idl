; tracersCosmo.pro
; dev for tracer particles related to cosmological boxes
; dnelson jan.2012

; cosmoTracerParents(): return indices (or optionally IDs) of parent gas cells of all tracers

function cosmoTracerParents, sP=sP, getInds=getInds, getIDs=getIDs
 
  if (not keyword_set(getInds) and not keyword_set(getIDs)) then begin
    print,'Error: No return specified.'
    stop
  endif
  
  ; check for save file
  saveFilename = sP.derivPath + 'trPar_' + sP.savPrefix + str(sP.res) + '_' + str(sP.snap) + '.sav'
  
  if file_test(saveFilename) then begin
    restore, saveFilename
  endif else begin
    ; load
    gas_pos = loadSnapshotSubset(sP.simPath, snapNum=sP.snap, field='pos', partType='gas')
    tr_pos  = loadSnapshotSubset(sP.simPath, snapNum=sP.snap, field='pos', partType='tracer')
    
    nTr  = (size(tr_pos))[2]
    nGas = (size(gas_pos))[2]
    
    par_ind  = lonarr(nTr)
    
    ; IDL (DEBUG!) - calculate nearest neighbor to each tracer position
    if 0 then begin 
    par_ind2 = lonarr(nTr)
      for i=0,nTr-1 do begin
        if (i mod 10 eq 0) then print,string(ceil(float(i)),format='(I3)')+'% done.'
        
        ; calculate distances to all gas particles
        dx = reform(gas_pos[0,*] - tr_pos[0,i])
        dy = reform(gas_pos[1,*] - tr_pos[1,i])
        dz = reform(gas_pos[2,*] - tr_pos[2,i])
        
        ; correct for periodic B.C.
        correctPeriodicDistVecs, dx, sP=sP
        correctPeriodicDistVecs, dy, sP=sP
        correctPeriodicDistVecs, dz, sP=sP
        
        ; find minimum distance
        dists = sqrt( dx*dx + dy*dy + dz*dz )
        w = where(dists eq min(dists),count)
        
        if (count ne 1) then begin
          print,'WARNING'
          stop
        endif
        
        par_ind2[i] = w[0]
      endfor
    endif
    
    ; EXTERNAL C - calculate nearest neighbor to each tracer position
    par_ind = calcNN(gas_pos,tr_pos,boxSize=sP.boxSize,ndims=3)

    ; save
    save,par_ind,filename=saveFilename
    print,'Saved: ',saveFilename
  endelse
  
  ; if IDs requested, load gas IDs and do crossmatch
  if keyword_set(getIDs) then begin
    gas_ids = loadSnapshotSubset(sP.simPath, snapNum=sP.snap, field='ids', partType='gas')
    
    ; return parent gas ids
    return,gas_ids[par_ind]
    
  endif else begin
    ; return parent gas inds
    return, par_ind
  endelse

end

; getCosmoTracerPos(): return time series of tracer positions for a given number of tracers nearest
;                      to a given ending position
;
; targetPos  : select N tracers closest to the specified ending position
; trackDegen : track all tracers with degenerate ending position

function getCosmoTracerPos, snapPath=snapPath, snapRange=snapRange, $
                            targetPos=targetPos, trackDegen=trackDegen, verbose=verbose

  useMatch     = 0    ; use match for ID location instead of where loop
  distDegenTol = 1e-2 ; distance tolerance for degenerate positions

  ; arrays
  nSnaps  = (snapRange[0]-snapRange[1]+1) / snapRange[2]
  times   = fltarr(nSnaps)
  
  ; verify spacing ok
  nSnaps2 = (float(snapRange[0])-snapRange[1]+1.0) / float(snapRange[2])
  if (nSnaps ne nSnaps2) then begin
    print,'Error: Spacing must evenly divide snapshot range.'
    return,0
  endif
  
  ; loop over requested snapshots
  k = 0
  
  for snap=snapRange[0],snapRange[1],-snapRange[2] do begin
    if keyword_set(verbose) then if (snap mod verbose eq 0) then $
      print,'snap: ',str(snap)
    
    ; load header and store time
    h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
    times[k] = h.time  
      
    pos = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='pos')
    ids = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='ids')
    
    ; make tracer selection on first snapshot
    if (snap eq snapRange[0]) then begin      

      ; NORMAL = select nearest tracer to each targetPos
      if not keyword_set(trackDegen) then begin
        ; arrays
        numTracers = (size(targetPos))[0]
        trPos      = dblarr(nSnaps,numTracers,3)
        trDist     = dblarr(nSnaps,numTracers)
        
        idTargets  = lonarr(numTracers)
        
        ; find IDs of targets
        for j=0,numTracers-1 do begin
          dists = abs( sqrt( (pos[0,*]-targetPos[0])*(pos[0,*]-targetPos[0]) + $
                             (pos[1,*]-targetPos[1])*(pos[1,*]-targetPos[1]) + $
                             (pos[2,*]-targetPos[2])*(pos[2,*]-targetPos[2]) ) )
          w = where(dists eq min(dists),count)
        
          ; only one and get ID
          if (count gt 1) then w = w[0]
          idTargets[j] = ids[w]
          
          if keyword_set(verbose) then $
            print,'selected id='+str(idTargets[j])+' rad='+str(rad[w])
        endfor
      endif
      
      ; TRACKDEGEN = select tracers near (<1e-6) targetPos
      if keyword_set(trackDegen) then begin
        ; make selection
        dists = abs( sqrt( (pos[0,*]-targetPos[0])*(pos[0,*]-targetPos[0]) + $
                           (pos[1,*]-targetPos[1])*(pos[1,*]-targetPos[1]) + $
                           (pos[2,*]-targetPos[2])*(pos[2,*]-targetPos[2]) ) )
        w = where(dists le distDegenTol,count)
        
        if (count eq 0) then begin
          print,'Error: No tracers near specified targetPos.'
          stop
        endif else begin
          print,'Found ['+str(count)+'] tracers near targetPos.'
        endelse
        
        ; arrays
        numTracers = count
        trPos      = dblarr(nSnaps,numTracers,3)
        ;trDist     = dblarr(nSnaps,numTracers)
        
        idTargets  = lonarr(numTracers)
        
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
      
    trPos[k,*,0] = pos[0,ids_ind]
    trPos[k,*,1] = pos[1,ids_ind]
    trPos[k,*,2] = pos[2,ids_ind]
    ;trDist[k,*]  = sqrt( (x[ids_ind]-boxCen)*(x[ids_ind]-boxCen) + $
    ;                     (y[ids_ind]-boxCen)*(y[ids_ind]-boxCen) )      
    endif else begin
      ; locate IDs using where loop
      for j=0,numTracers-1 do begin
        ind = where(ids eq idTargets[j])
        
        if keyword_set(verbose) then if (snap mod verbose eq 0) then $
          print,' ['+str(j)+'] refound id at ind='+str(ind)
    
        trPos[k,j,0] = pos[0,ind]
        trPos[k,j,1] = pos[1,ind]
        trPos[k,j,2] = pos[2,ind]
        ;trDist[k,j]   = sqrt( (x[ind]-boxCen)*(x[ind]-boxCen) + (y[ind]-boxCen)*(y[ind]-boxCen) )
      endfor
    endelse
      
    k += 1
    
  endfor ;snap
  
  r = {trPos:trPos,idTargets:idTargets,times:times,boxSize:h.boxSize} ;trDist
  return, r

end

; cosmoTracerTraj(): plot a few individual tracer trajectories over contour of the disk gas density

pro cosmoTracerTraj
 
  ; config
  res = 128
  run = 'dev.tracer'
  
  ;targetPos = [2614.26,123.019,3410.30] ;ckpc ;orig @snap95
  targetPos = [1731.88,3787.28,3495.91]  ;ckpc, _test2 @snap6

  ; paths and snapshots
  sP    = simParams(res=res,run=run)
  units = getUnits() ;colors
  
  snapRange = [sP.snapRange[1],0,1]
  plotBase  = 'cosmoRad_'+run+'_'
  
  ; get tracer positions
  saveFilename = sP.derivPath + plotBase + 'trTraj_'+string(snapRange[1],format='(I3.3)')+'_'+$
            string(snapRange[0],format='(I3.3)')+'.sav'
  
  if file_test(saveFilename) then begin
    restore,saveFilename,/verbose
  endif else begin
    tp = getCosmoTracerPos(snapPath=sP.simPath, snapRange=snapRange, targetPos=targetPos, $
                           /trackDegen, verbose=1)
    save,tp,filename=saveFilename
  endelse
            
  ; plots
  start_PS, sP.plotPath + plotBase + 'trTraj_'+string(snapRange[1],format='(I3.3)')+'_'+$
            string(snapRange[0],format='(I3.3)')+'.eps'
  
    fsc_text,0.5,0.8,"tracer trajectory - cosmo "+run+"",/normal,alignment=0.5
    
    ; PLOT ONE - (x,y)
    boxsize = 30.0 ;ckpc
    xrange = minmax(tp.trPos[*,*,0])+[-boxsize,boxsize] > [0,0]
    yrange = minmax(tp.trPos[*,*,1])+[-boxsize,boxsize] > [0,0]
    zrange = minmax(tp.trPos[*,*,2])+[-boxsize,boxsize] > [0,0]
    
    ;boxsize = 10.0 ;ckpc
    ;xrange = [tp.trPos[n_elements(tp.times)-1,0,0] - boxsize, $
    ;          tp.trPos[n_elements(tp.times)-1,0,0] + boxsize]
    ;yrange = [tp.trPos[n_elements(tp.times)-1,0,1] - boxsize, $
    ;          tp.trPos[n_elements(tp.times)-1,0,1] + boxsize]
    ;zrange = [tp.trPos[n_elements(tp.times)-1,0,2] - boxsize, $
    ;          tp.trPos[n_elements(tp.times)-1,0,2] + boxsize]
              
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,charsize=!p.charsize-0.5,$
         xtitle="x [ckpc]",ytitle="y [ckpc]",position=[0.1,0.2,0.45,0.75],/noerase

    ; overplot tracer trajectories
    for i=0,n_elements(tp.idTargets)-1 do begin
      xPts = tp.trPos[*,i,0]
      yPts = tp.trPos[*,i,1]
      
      ; underplot circles at targetPos
      tvcircle,boxsize/20.0,targetPos[0],targetPos[1],color=fsc_color('light gray'),/data
      
      plotsym,0,/fill
      fsc_plot, xPts, yPts, psym=-8, symsize=0.3, color=getColor(i), /overplot
      
    endfor
    
    ; PLOT TWO (x,z)
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=zrange,xstyle=1,ystyle=1,charsize=!p.charsize-0.5,$
         xtitle="x [ckpc]",ytitle="z [ckpc]",position=[0.55,0.2,0.95,0.75],/noerase
         
    for i=0,n_elements(tp.idTargets)-1 do begin
      xPts = tp.trPos[*,i,0]
      zPts = tp.trPos[*,i,2]
      
      ; underplot circles at targetPos
      tvcircle,boxsize/20.0,targetPos[0],targetPos[2],color=fsc_color('light gray'),/data
      
      plotsym,0,/fill
      fsc_plot, xPts, zPts, psym=-8, symsize=0.3, color=getColor(i), /overplot
      
      ; legend
      zStrs = "z="+string(1.0/tp.times-1.0,format='(f4.1)')
      legend,zStrs,textcolors=['black'],/right,$
             margin=0.25,charsize=!p.charsize-0.5,box=0
      
    endfor
  
  end_PS, pngResize=50

end

; cosmoTracerGasDensComp(): compare spatially estimated tracer density to parent cell gas density

pro cosmoTracerGasDensComp, redshift=redshift, nNGB=nNGB

  ; config
  res = 128
  run = 'dev.tracer'
  
  sP    = simParams(res=res,run=run,redshift=redshift)
  units = getUnits() ;colors

  ; load
  h = loadSnapshotHeader(sP.simPath, snapNum=sP.snap)

  ; tracer mass estimate (uniform density ICs assumption)
  mass_gas = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='mass',/verbose)
  tr_mass = total(mass_gas) / h.nPartTot[3]
  mass_gas = !NULL
  
  ; estimate tracer densities
  saveFilename = sP.derivPath + sP.savPrefix + str(res) + '.trDens.nNGB=' + str(nNGB) + '.snap=' + $
                 str(sP.snap) + '.sav'
                 
  if file_test(saveFilename) then begin
    restore,saveFilename,/verbose
  endif else begin
    tr_pos = loadSnapshotSubset(sP.simPath, snapNum=sP.snap, field='pos', partType='tracer')
    tr_dens = estimateDensityTophat(tr_pos,mass=tr_mass,ndims=3,nNGB=nNGB,boxSize=h.boxSize)
    tr_pos = !NULL
    save,tr_dens,filename=saveFilename
  endelse
  
  ; load tracer parents and gas densities
  tr_par_ind = cosmoTracerParents(sP=sP, /getInds)
  
  gas_dens = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='density',/verbose)
  
  ; reorder gas densities to correspond to tr_dens parents
  gas_dens = gas_dens[tr_par_ind]

  ; plot
  if 0 then begin
  start_PS, sP.plotPath+sP.savPrefix+str(res)+'.densComp1.nNGB='+str(nNGB)+'.snap='+str(sP.snap)+'.eps'
  
    xrange = [1e-12,1e-8]
    yrange = xrange
    
    w = where(tr_dens ge xrange[0] and tr_dens le xrange[1])
    
    ;ypts = (tr_dens - gas_dens) / gas_dens
    ypts = gas_dens
  
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,charsize=!p.charsize-0.5,$
         xtitle="tracer density",ytitle="parent gas cell density"
         
    fsc_plot,tr_dens[w],ypts[w],psym=3,/overplot
         
  end_PS, pngResize=50, /deletePS

  start_PS, sP.plotPath+sP.savPrefix+str(res)+'.densComp2.nNGB='+str(nNGB)+'.snap='+str(sP.snap)+'.eps'
  
    xrange = [1e-12,1e-7]
    yrange = [0.0,5.0]
    
    w = where(tr_dens ge xrange[0] and tr_dens le xrange[1])
    
    ypts = tr_dens / gas_dens
  
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,charsize=!p.charsize-0.5,$
         xtitle="gas density",ytitle="tracer density / parent gas cell density"
         
    fsc_plot,gas_dens[w],ypts[w],psym=3,/overplot
         
  end_PS, pngResize=50, /deletePS
  
  ; ratio to critical at z=3 and alog
  gas_dens2 = rhoRatioToCrit(gas_dens,redshift=3.0)
  gas_dens2 = alog10(gas_dens2)
  
  start_PS, sP.plotPath+sP.savPrefix+str(res)+'.densComp3.nNGB='+str(nNGB)+'.snap='+str(sP.snap)+'.eps'
  
    xrange = [-2.0,8.0]
    yrange = [0.0,8.0]
    
    w = where(gas_dens2 ge xrange[0] and gas_dens2 le xrange[1])
    
    ypts = tr_dens / gas_dens
  
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,charsize=!p.charsize-0.5,$
         xtitle="gas density / baryon crit density at z=3",ytitle="tracer dens / parent gas cell dens"
         
    fsc_plot,gas_dens2[w],ypts[w],psym=3,/overplot
         
  end_PS, pngResize=50, /deletePS
  
  endif ;0
  
  start_PS, sP.plotPath+sP.savPrefix+str(res)+'.densComp4.nNGB='+str(nNGB)+'.snap='+str(sP.snap)+'.eps'
  
    xrange = [0.0,5.0]
    yrange = [1e2,1e5]
    
    w = where(tr_dens ge xrange[0] and tr_dens le xrange[1])
    
    ypts = tr_dens / gas_dens
    ypts = ypts[where(ypts le xrange[1])]
  
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,charsize=!p.charsize-0.5,$
         xtitle="tracer density / parent cell gas density",ytitle="N tracers",/ylog
         
    fsc_plot,[1.0,1.0],yrange,color=fsc_color('light gray'),/overplot
    
    binsize = 0.02
    
    hist = histogram(ypts,binsize=binsize,min=xrange[0],max=xrange[1],loc=loc)
    
    fsc_plot,loc+binsize/2.0,hist,line=0,/overplot
         
  end_PS
  
end

; cosmoTracerGasDensCompNGB(): compare spatially estimated tracer density to parent cell gas density
;                              as a function of nNGB used to calculate tracer densities

pro cosmoTracerGasDensCompNGB, redshift=redshift

  if not keyword_set(redshift) then stop

  ; config
  res = 128
  run = 'dev.tracer'
  
  nNGBs = [32,64,128,256]
  
  binsize = 0.02
  
  sP    = simParams(res=res,run=run,redshift=redshift)
  units = getUnits() ;colors

  ; load
  h = loadSnapshotHeader(sP.simPath, snapNum=sP.snap)

  ; tracer mass estimate (uniform density ICs assumption)
  mass_gas = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='mass',/verbose)
  tr_mass = total(mass_gas) / h.nPartTot[3]
  mass_gas = !NULL
  
  ; load tracer parents and gas densities
  tr_par_ind = cosmoTracerParents(sP=sP, /getInds)
  
  gas_dens = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='density',/verbose)  
  
  ; reorder gas densities to correspond to tr_dens parents
  gas_dens = gas_dens[tr_par_ind]  
  
  ; start plot
  start_PS, sP.plotPath+sP.savPrefix+str(res)+'.densCompNGB.snap='+str(sP.snap)+'.eps'
  
    xrange = [0.0,5.0]
    yrange = [1e2,1e5]
  
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,charsize=!p.charsize-0.5,$
         xtitle="tracer density / parent cell gas density",ytitle="N tracers",/ylog,$
         title=run+" "+str(res)+" z="+string(redshift,format='(f3.1)')
        
    fsc_plot,[1.0,1.0],yrange,color=fsc_color('light gray'),/overplot
         
    legendColors = []
    legendStrs   = []
  
    foreach nNGB,nNGBs,i do begin

      ; estimate tracer densities
      saveFilename = sP.derivPath + sP.savPrefix + str(res) + '.trDens.nNGB=' + str(nNGB) + '.snap=' + $
                     str(sP.snap) + '.sav'
                     
      if file_test(saveFilename) then begin
        restore,saveFilename,/verbose
      endif else begin
        tr_pos = loadSnapshotSubset(sP.simPath, snapNum=sP.snap, field='pos', partType='tracer')
        tr_dens = estimateDensityTophat(tr_pos,mass=tr_mass,ndims=3,nNGB=nNGB,boxSize=h.boxSize)
        tr_pos = !NULL
        save,tr_dens,filename=saveFilename
      endelse
      
      ypts = tr_dens / gas_dens
      ypts = ypts[where(ypts le xrange[1])]
      
      ; plot
      hist = histogram(ypts,binsize=binsize,min=xrange[0],max=xrange[1],loc=loc)
      
      fsc_plot,loc+binsize/2.0,hist,line=0,/overplot,color=getColor(i)
      
      legendStrs   = [legendStrs,'nNGB = '+str(nNGB)]
      legendColors = [legendColors,getColor(i,/name)]
      
    endforeach
    
    ; legend
    legend,legendStrs,textcolors=legendColors,/right,margin=0.25,charsize=!p.charsize-0.5,box=0
  
  ; end plot
  end_PS

end

; cosmoTracerGasDensCompRedshift(): compare spatially estimated tracer density to parent cell gas density
;                                   as a function of redshift

pro cosmoTracerGasDensCompRedshift

  ; config
  res = 128
  run = 'dev.tracer'
  
  nNGB = 32
  
  redshifts = [3.0,4.0,6.0,10.0,22.0]
  
  binsize = 0.02
  
  units = getUnits() ;colors
  
  ; start plot
  start_PS, sP.plotPath+sP.savPrefix+str(res)+'.densCompNGB.snap='+str(sP.snap)+'.eps'
  
    xrange = [0.0,5.0]
    yrange = [1e2,1e5]
  
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,charsize=!p.charsize-0.5,$
         xtitle="tracer density / parent cell gas density",ytitle="N tracers",/ylog,$
         title=run+" "+str(res)+" z="+string(redshift,format='(f3.1)')
        
    fsc_plot,[1.0,1.0],yrange,color=fsc_color('light gray'),/overplot
         
    legendColors = []
    legendStrs   = []  
  
  ; loop over redshifts
  foreach redshift,redshifts,i do begin
  
    sP = simParams(res=res,run=run,redshift=redshift)
    
    ; load
    h = loadSnapshotHeader(sP.simPath, snapNum=sP.snap)
  
    ; tracer mass estimate (uniform density ICs assumption)
    mass_gas = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='mass',/verbose)
    tr_mass = total(mass_gas) / h.nPartTot[3]
    mass_gas = !NULL
    
    ; load tracer parents and gas densities
    tr_par_ind = cosmoTracerParents(sP=sP, /getInds)
    
    gas_dens = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='density',/verbose)  
    
    ; reorder gas densities to correspond to tr_dens parents
    gas_dens = gas_dens[tr_par_ind]  
  

    ; estimate tracer densities
    saveFilename = sP.derivPath + sP.savPrefix + str(res) + '.trDens.nNGB=' + str(nNGB) + '.snap=' + $
                   str(sP.snap) + '.sav'
                   
    if file_test(saveFilename) then begin
      restore,saveFilename,/verbose
    endif else begin
      tr_pos = loadSnapshotSubset(sP.simPath, snapNum=sP.snap, field='pos', partType='tracer')
      tr_dens = estimateDensityTophat(tr_pos,mass=tr_mass,ndims=3,nNGB=nNGB,boxSize=h.boxSize)
      tr_pos = !NULL
      save,tr_dens,filename=saveFilename
    endelse
      
    ypts = tr_dens / gas_dens
    ypts = ypts[where(ypts le xrange[1])]
    
    ; plot
    hist = histogram(ypts,binsize=binsize,min=xrange[0],max=xrange[1],loc=loc)
    
    fsc_plot,loc+binsize/2.0,hist,line=0,/overplot,color=getColor(i)
    
    legendStrs   = [legendStrs,'z = '+string(redshift,format='(f4.1)')]
    legendColors = [legendColors,getColor(i,/name)]

  endforeach
  
  ; end plot
  legend,legendStrs,textcolors=legendColors,/right,margin=0.25,charsize=!p.charsize-0.5,box=0
  end_PS

end
