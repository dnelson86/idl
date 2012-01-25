; tracersCosmo.pro
; dev for tracer particles related to cosmological boxes
; dnelson jan.2012

; getTracerSpatialDens(): wrapper to do tophat density calculation for tracers and save result

function getTracerSpatialDens, sP=sP, nNGB=nNGB

  saveFilename = sP.derivPath + sP.savPrefix + str(sP.res) + '.trDens.nNGB=' + str(nNGB) + '.snap=' + $
                 str(sP.snap) + '.sav'
                 
  if file_test(saveFilename) then begin
    restore,saveFilename,/verbose
  endif else begin
    ; load snapshot header for boxSize
    h = loadSnapshotHeader(sP.simPath, snapNum=sP.snap)
    
    ;tr_mass = total(mass_gas) / h.nPartTot[3] ; should be mass gas at t=0
    
    ; load positions and calculate densities via HSMLs
    tr_pos  = loadSnapshotSubset(sP.simPath, snapNum=sP.snap, field='pos', partType='tracer')
    tr_dens = estimateDensityTophat(tr_pos,mass=sP.trMassConst,ndims=3,nNGB=nNGB,boxSize=h.boxSize)
    tr_pos  = !NULL
    
    ; save
    save,tr_dens,filename=saveFilename
  endelse
  
  return,tr_dens
  
end

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

; cosmoCompMassFunctions(): compare gas, tracer, and DM FoF mass functions

pro cosmoCompMassFunctions

  ; config
  res      = 128
  run      = 'dev.tracer.noref'
  redshift = 1.0
  
  ; load group catalog
  sP    = simParams(res=res,run=run,redshift=redshift)
  units = getUnits()

  h  = loadSnapshotHeader(sP.simPath, snapNum=sP.snap)
  gc = loadGroupCat(sP.simPath,sP.snap,/verbose)
  
  ; load gas masses, calculate dm and tr masses
  gas_mass = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='mass')
  dm_mass  = h.massTable[1]  
  
  ; calculate halo masses
  hm_gas  = reform(gc.groupMassType[0,*]) * units.UnitMass_in_Msun
  hm_dm   = reform(gc.groupMassType[1,*]) * units.UnitMass_in_Msun
  hm_star = reform(gc.groupMassType[4,*]) * units.UnitMass_in_Msun
  hm_tr   = reform(gc.groupLenType[3,*]) * sP.trMassConst * units.UnitMass_in_Msun
  ;hm_dm2  = reform(gc.groupLenType[1,*]) * dm_mass * units.UnitMass_in_Msun ;debug check
  
  ; nonzero only
  hm_gas  = hm_gas[where(hm_gas ne 0)]
  hm_dm   = hm_dm[where(hm_dm ne 0)]
  hm_star = hm_star[where(hm_star ne 0)]
  hm_tr   = hm_tr[where(hm_tr ne 0)]
  hm_bar  = [hm_gas,hm_star] ;baryonic
  
  print,'Found: ['+str(n_elements(hm_gas))+'] gas, ['+str(n_elements(hm_dm))+'] dm, ['+$
        str(n_elements(hm_star))+'] stars, ['+str(n_elements(hm_tr))+'] tracers.'
  
  ; sort ascending
  hm_gas  = hm_gas[sort(hm_gas)]
  hm_dm   = hm_dm[sort(hm_dm)]
  hm_star = hm_star[sort(hm_star)]
  hm_tr   = hm_tr[sort(hm_tr)]
  hm_bar  = hm_bar[sort(hm_bar)]
  
  ; y-vals (cumulative number count) and normalize by box volume
  y_gas  = reverse(indgen(n_elements(hm_gas)) + 1)   / (h.boxSize/1000)^3.0 ;Mpc
  y_dm   = reverse(indgen(n_elements(hm_dm)) + 1)    / (h.boxSize/1000)^3.0
  y_star = reverse(indgen(n_elements(hm_star)) + 1)  / (h.boxSize/1000)^3.0
  y_tr   = reverse(indgen(n_elements(hm_tr)) + 1)    / (h.boxSize/1000)^3.0
  y_bar  = reverse(indgen(n_elements(hm_bar)) + 1)   / (h.boxSize/1000)^3.0
  
  ; plot
  start_PS,sP.plotPath+sP.savPrefix+str(res)+'.massFuncs.snap='+str(sP.snap)+'.eps'
    xrange = [2e7,4e12]
    yrange = [1e-4,2e0]
    
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,/xlog,$
         xtitle="",xtickname=replicate(' ',10),$
         ytitle="number ("+textoidl("\geq")+" M) [h"+textoidl("^3")+" Mpc"+textoidl("^{-3}")+"]",$
         title=str(res)+textoidl("^3")+" z="+string(redshift,format='(f3.1)')+" (FoF Mass Functions)",$
         position=[0.18,0.35,0.9,0.9]
         
    fsc_plot,hm_dm,y_dm,line=0,/overplot,color=getColor(1)
    fsc_plot,hm_star,y_star,line=0,/overplot,color=getColor(2)
    fsc_plot,hm_gas,y_gas,line=0,/overplot,color=getColor(3)
    fsc_plot,hm_tr,y_tr,line=0,/overplot,color=getColor(7)
    fsc_plot,hm_bar,y_bar,line=0,/overplot,color=getColor(8)
    
    ; legend
    legend,['dm','star','gas','tracer','gas+stars'],$
           textcolors=getColor([1,2,3,7,8],/name),/bottom,/left,box=0,margin=0.1
           
    ; gas/tr residual plot
    yrange = [0.66,1.34]
    fsc_plot,[0],[0],/nodata,/noerase,xrange=xrange,yrange=yrange,/xs,/ys,/xlog,$
             xtitle="mass [h"+textoidl("^{-1}")+" M"+textoidl("_{sun}")+"]",$
             ytitle="ratio",ytickv=[0.8,1.0,1.2],yticks=2,$
             position=[0.18,0.15,0.9,0.35]
             
    ; just interpolate both onto a set of masses then compare
    nbins = 100
    res_pts = 10.0^( findgen(nbins+1)/nbins * (10.9-7.5) + 7.5 )
    gas_res = interpol(y_gas,hm_gas,res_pts)
    tr_res  = interpol(y_tr,hm_tr,res_pts)
    
    ; plot
    fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[1.1,1.1],line=2,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[0.9,0.9],line=2,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[1.2,1.2],line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[0.8,0.8],line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,res_pts,tr_res/gas_res,line=0,color=getColor(3),/overplot
    
    ; do the same for the baryonic MF instead of just gas
    res_pts = 10.0^( findgen(nbins+1)/nbins * (11.1-7.5) + 7.5 )
    tr_res  = interpol(y_tr,hm_tr,res_pts)
    bar_res = interpol(y_bar,hm_bar,res_pts)
    fsc_plot,res_pts,tr_res/bar_res,line=0,color=getColor(8),/overplot
    
    ; legend
    legend,['tr/gas','tr/bar'],textcolors=getColor([3,8],/name),/right,/top,box=0
    
  end_PS
  
  stop  
  
end

; cosmoCompRadProfiles(): compare gas, tracer, and DM radial profiles of halos

pro cosmoCompRadProfiles, massBin=massBin

  if not keyword_set(massBin) then stop

  ; config
  res      = 128
  run      = 'dev.tracer.nonrad'
  redshift = 1.0
  
  ; load group catalog
  sP    = simParams(res=res,run=run,redshift=redshift)
  units = getUnits() ;colors

  h  = loadSnapshotHeader(sP.simPath, snapNum=sP.snap)
  gc = loadGroupCat(sP.simPath,sP.snap,/verbose)

  ; halo selection (manual)
  ;haloIDs = [0]
  
  ;plotTitle = str(res)+textoidl("^3")+" z="+string(redshift,format='(f3.1)')+" - haloID="+str(haloIDs[0])
  ;saveTag   = 'halo='+str(haloIDs[0])

  ; halo selection (mass bin)
  haloIDs = where(gc.groupMass ge massBin[0] and gc.groupMass lt massBin[1],countMassBin)
  print,'Found ['+str(countMassBin)+'] halos in mass bin.'

  plotTitle = str(res)+textoidl("^3")+" z="+string(redshift,format='(f3.1)')+" ("+$
              string(massBin[0],format='(f4.1)')+" < log(M) < "+string(massBin[1],format='(f4.1)')+")"
  saveTag   = 'massbin='+string(massBin[0],format='(f4.1)')+"-"+string(massBin[1],format='(f4.1)')

  ; setup binning
  nbins  = 20
  minmax = alog10([0.02,1.5]) ; ratio to r200
 
  rho_dm    = fltarr(nbins)
  rho_gas   = fltarr(nbins)
  rho_stars = fltarr(nbins)
  rho_tr    = fltarr(nbins)
  
  radBins = 10.0^( findgen(nbins+1)/nbins*(minmax[1]-minmax[0]) + minmax[0] )
  midBins = 10.0^( (findgen(nbins)+0.5)/nbins*(minmax[1]-minmax[0]) + minmax[0] )

  ; load gas,tr,dm,star positions
  pos_gas   = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='pos')
  pos_tr    = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='tracer',field='pos')
  pos_dm    = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='dm',field='pos')
  ;pos_stars = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='stars',field='pos')
  
  gas_mass   = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='mass')
  ;stars_mass = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='stars',field='mass')
  dm_mass    = h.massTable[1]  

  ; locate halo
  foreach haloID,haloIDs do begin
  
    haloPos = gc.groupPos[*,haloID]
    haloRad = gc.group_r_crit200[haloID]
    
    ; calculate radii and make radial cut
    rad = reform( sqrt( (pos_gas[0,*]-haloPos[0])*(pos_gas[0,*]-haloPos[0]) + $
                        (pos_gas[1,*]-haloPos[1])*(pos_gas[1,*]-haloPos[1]) + $
                        (pos_gas[2,*]-haloPos[2])*(pos_gas[2,*]-haloPos[2]) ) )

    gas_ind = where(rad le haloRad*2.0,gas_count)
    rad_gas = rad[gas_ind]/haloRad
    
    rad = reform( sqrt( (pos_tr[0,*]-haloPos[0])*(pos_tr[0,*]-haloPos[0]) + $
                        (pos_tr[1,*]-haloPos[1])*(pos_tr[1,*]-haloPos[1]) + $
                        (pos_tr[2,*]-haloPos[2])*(pos_tr[2,*]-haloPos[2]) ) )
    
    tr_ind = where(rad le haloRad*2.0,tr_count)
    rad_tr = rad[tr_ind]/haloRad
    
    rad = reform( sqrt( (pos_dm[0,*]-haloPos[0])*(pos_dm[0,*]-haloPos[0]) + $
                        (pos_dm[1,*]-haloPos[1])*(pos_dm[1,*]-haloPos[1]) + $
                        (pos_dm[2,*]-haloPos[2])*(pos_dm[2,*]-haloPos[2]) ) )
    
    dm_ind = where(rad le haloRad*2.0,dm_count)
    rad_dm = rad[dm_ind]/haloRad
    
    ;rad = reform( sqrt( (pos_stars[0,*]-haloPos[0])*(pos_stars[0,*]-haloPos[0]) + $
    ;                    (pos_stars[1,*]-haloPos[1])*(pos_stars[1,*]-haloPos[1]) + $
    ;                    (pos_stars[2,*]-haloPos[2])*(pos_stars[2,*]-haloPos[2]) ) )
    ;
    ;stars_ind = where(rad le haloRad*2.0,stars_count)
    ;rad_stars = rad[stars_ind]/haloRad
    stars_count = 0
    
    ; check for degenerate radii
    w = where(finite(rad_gas),comp=wc,ncomp=ncomp)
    if (ncomp gt 0) then rad_gas[wc] = 0.0
    w = where(finite(rad_tr),comp=wc,ncomp=ncomp)
    if (ncomp gt 0) then rad_tr[wc] = 0.0
    w = where(finite(rad_dm),comp=wc,ncomp=ncomp)
    if (ncomp gt 0) then rad_dm[wc] = 0.0
    ;w = where(finite(rad_stars),comp=wc,ncomp=ncomp)
    ;if (ncomp gt 0) then rad_stars[wc] = 0.0
    
    print,'(' + str(haloID) + ') Found ['+str(gas_count)+'] gas ['+str(tr_count)+'] tracer ['+$
          str(dm_count)+'] dm ['+str(stars_count)+'] stars inside cut.'
                     
    ; subselect gas,stars masses
    gas_mass_sub   = gas_mass[gas_ind]
    ;stars_mass_sub = stars_mass[stars_ind]
    
    ; do binning 
    for i=0,nbins-1 do begin
      ; shell volume normalization
      vol = 4*!pi/3 * (radBins[i+1]^3.0 - radBins[i]^3.0) * haloRad^3.0
      
      w = where(rad_gas gt radBins[i] and rad_gas le radBins[i+1],count)
      if (count gt 0) then rho_gas[i] += total(gas_mass_sub[w]) / vol
      
      ;w = where(rad_gas gt radBins[i] and rad_stars le radBins[i+1],count)
      ;if (count gt 0) then rho_stars[i] += total(stars_mass_sub[w]) / vol
      
      w = where(rad_dm gt radBins[i] and rad_dm le radBins[i+1],count)
      if (count gt 0) then rho_dm[i] += dm_mass * count / vol
      
      w = where(rad_tr gt radBins[i] and rad_tr le radBins[i+1],count)
      if (count gt 0) then rho_tr[i] += sP.trMassConst * count / vol
    endfor
  
  endforeach
  
  ; normalize stacked profiles by number of halos
  rho_gas   /= n_elements(haloIDs)
  rho_tr    /= n_elements(haloIDs)
  rho_dm    /= n_elements(haloIDs)
  ;rho_stars /= n_elements(haloIDs)

  ; plot
  start_PS,sP.plotPath+sP.savPrefix+str(res)+'.'+saveTag+'.radProfiles.snap='+str(sP.snap)+'.eps'
    xrange = 10.0^minmax
    yrange = [1e-8,1e-2]
    
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
         xtitle="r / r"+textoidl("_{vir}"),ytitle="density [h"+textoidl("^2")+" M"+textoidl("_{sun}")+$
         " kpc"+textoidl("^{-3}")+"]",title=plotTitle,/ylog,/xlog
         
    fsc_plot,midBins,rho_gas,line=0,/overplot,color=getColor(1)
    fsc_plot,midBins,rho_dm,line=0,/overplot,color=getColor(2)
    fsc_plot,midBins,rho_tr,line=0,/overplot,color=getColor(3)
    ;fsc_plot,midBins,rho_stars,line=0,/overplot,color=getColor(7)
    ;fsc_plot,midBins,rho_gas+rho_stars,line=0,/overplot,color=getColor(8)
    
    ; legend
    legend,['gas','dm','tracer','stars','gas+stars','('+str(n_elements(haloIDs))+' halos)'],$
           textcolors=[getColor([1,2,3,7,8],/name),'black'],$
           /right, box=0, margin=0.25
  end_PS
  
  stop
end

; cosmoTracerParentHisto():

pro cosmoTracerParentHisto

  ; config
  res = 128
  run = 'dev.tracer.noref'
  
  redshifts = [5.0,3.0,2.0,1.0]
  nNGB = 32
  units = getUnits() ;colors
  
  ; start plot
  sP = simParams(res=res,run=run,redshift=redshift)
  start_PS, sP.plotPath+sP.savPrefix+str(res)+'.parHisto.nNGB='+str(nNGB)+'.eps'
  
  num = 20
  xrange = [0,num]
  yrange = [1e1,4e6]
  
  fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
       xtitle="number of tracers in a gas cell",ytitle="N gas cells",$
       title=str(res)+textoidl("^3")+" "+run,/ylog  
  
  legendColors = []
  legendStrs   = []
  
  foreach redshift,redshifts,j do begin
  
    sP    = simParams(res=res,run=run,redshift=redshift)
    
    ; load
    h = loadSnapshotHeader(sP.simPath, snapNum=sP.snap)
  
    ; load tracer parents and sort
    tr_par_ind = cosmoTracerParents(sP=sP, /getInds)
    tr_par_ind = tr_par_ind[sort(tr_par_ind)]
    ;tr_par_ind = [0,1,2,3,4,4,5,6,6,6,7,9,15,16,16]
    
    par_histo = lonarr(10000)
    
    parCount = 1
    
    for i=1,n_elements(tr_par_ind)-1 do begin
      if (tr_par_ind[i] eq tr_par_ind[i-1]) then begin
        ; this parent matches the last
        parCount += 1
      endif else begin
        ; this parent different than the last
        par_histo[parCount] += 1
        
        ; reset counter
        parCount = 0
      endelse
    endfor
    
    ; maximum number of parents
    w = where(par_histo ne 0)
    print,'z='+str(redshift)+' max number of parents: ',max(w)+1

    ; overplot
    plotsym,0,/fill
    fsc_plot,indgen(num)+1,par_histo[0:num-1],psym=-8,thick=!p.thick+1.0,$
             color=getColor(j),/overplot
    
    legendStrs   = [legendStrs,'z = '+string(redshift,format='(f4.1)')]
    legendColors = [legendColors,getColor(j,/name)]
    
  endforeach
  
  ; legend
  legend,legendStrs,textcolors=legendColors,/right,/top,box=0
  
  ; end plot
  end_PS
  
  ; x-log plot
  sP = simParams(res=res,run=run,redshift=redshift)
  start_PS, sP.plotPath+sP.savPrefix+str(res)+'.parHisto2.nNGB='+str(nNGB)+'.eps'
  
  num = 500
  xrange = [1,num]
  yrange = [1,4e6]
  
  fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
       xtitle="number of tracers in a gas cell",ytitle="N gas cells",$
       title=str(res)+textoidl("^3")+" "+run,/ylog,/xlog
  
  legendColors = []
  legendStrs   = []
  
  foreach redshift,redshifts,j do begin
  
    sP    = simParams(res=res,run=run,redshift=redshift)
    
    ; load
    h = loadSnapshotHeader(sP.simPath, snapNum=sP.snap)
  
    ; load tracer parents and sort
    tr_par_ind = cosmoTracerParents(sP=sP, /getInds)
    tr_par_ind = tr_par_ind[sort(tr_par_ind)]
    
    par_histo = lonarr(10000)
    
    parCount = 1
    
    for i=1,n_elements(tr_par_ind)-1 do begin
      if (tr_par_ind[i] eq tr_par_ind[i-1]) then begin
        ; this parent matches the last
        parCount += 1
      endif else begin
        ; this parent different than the last
        par_histo[parCount] += 1
        
        ; reset counter
        parCount = 0
      endelse
    endfor
    
    ; maximum number of parents
    w = where(par_histo ne 0)
    print,'z='+str(redshift)+' max number of parents: ',max(w)+1

    ; overplot
    fsc_plot,indgen(num)+1,par_histo[0:num-1],line=0,thick=!p.thick+1.0,$
             color=getColor(j),/overplot
    
    legendStrs   = [legendStrs,'z = '+string(redshift,format='(f4.1)')]
    legendColors = [legendColors,getColor(j,/name)]
    
  endforeach
  
  ; legend
  legend,legendStrs,textcolors=legendColors,/right,/top,box=0
  
  ; end plot
  end_PS
  
  stop
end

; cosmoTracerGasDensComp(): compare spatially estimated tracer density to parent cell gas density

pro cosmoTracerGasDensComp, redshift=redshift, nNGB=nNGB

  ; config
  res = 128
  run = 'dev.tracer.noref'
  
  sP    = simParams(res=res,run=run,redshift=redshift)
  units = getUnits() ;colors

  ; load
  h = loadSnapshotHeader(sP.simPath, snapNum=sP.snap)
  
  ; get tracer and gas densities
  tr_dens  = getTracerSpatialDens(sP=sP, nNGB=nNGB)
  gas_dens = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='density',/verbose)
  
  ; load tracer parents and gas densities
  tr_par_ind = cosmoTracerParents(sP=sP, /getInds)
  
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

pro cosmoTracerGasDensCompNGB

  ; config
  res = 128
  run = 'dev.tracer.nonrad'
  
  redshift = 3.0
  
  nNGBs = [32,64,128,256]
  
  binsize = 0.02
  
  sP    = simParams(res=res,run=run,redshift=redshift)
  units = getUnits() ;colors

  ; load
  h = loadSnapshotHeader(sP.simPath, snapNum=sP.snap)
  
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
      tr_dens  = getTracerSpatialDens(sP=sP, nNGB=nNGB)
      
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
  run = 'dev.tracer.noref'
  
  nNGB = 32
  
  redshifts = [1.0,3.0,6.0,10.0,30.0]
  
  binsize = 0.02
  
  sP    = simParams(res=res,run=run)
  units = getUnits() ;colors
  
  ; start plot
  start_PS, sP.plotPath+sP.savPrefix+str(res)+'.densCompRedshift.nNGB='+str(nNGB)+'.eps'
  
    xrange = [0.0,5.0]
    yrange = [1e2,5e5]
  
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,charsize=!p.charsize-0.5,$
         xtitle="tracer density / parent cell gas density",ytitle="N tracers",/ylog,$
         title=run+" "+str(res)+" nNGB="+str(nNGB)
        
    fsc_plot,[1.0,1.0],yrange,color=fsc_color('light gray'),/overplot
         
    legendColors = []
    legendStrs   = []  
  
  ; loop over redshifts
  foreach redshift,redshifts,i do begin
  
    sP = simParams(res=res,run=run,redshift=redshift)

    ; load
    h = loadSnapshotHeader(sP.simPath, snapNum=sP.snap)
  
    ; load tracer parents and gas densities
    tr_par_ind = cosmoTracerParents(sP=sP, /getInds)
    gas_dens   = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='density',/verbose)  
    
    ; reorder gas densities to correspond to tr_dens parents
    gas_dens = gas_dens[tr_par_ind]  
  
    ; estimate tracer densities
    tr_dens  = getTracerSpatialDens(sP=sP, nNGB=nNGB)
      
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
stop
end
