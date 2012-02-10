; tracersCosmo.pro
; dev for tracer particles related to cosmological boxes
; dnelson feb.2012

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
; numTracers=1 : track N tracers closest to targetPos
; maxDist=1    : track all tracers within minDist of targetPos

function getCosmoTracerPos, snapPath=snapPath, snapRange=snapRange, numTracers=numTracers, $
                            targetPos=targetPos, maxDist=maxDist, verbose=verbose

  useMatch     = 1    ; use match for ID location instead of where loop

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
    
    ; skip over missing snapshots
    if (n_tags(h) eq 0) then begin
      if keyword_set(verbose) then if (snap mod verbose eq 0) then $
        print,' skipped'
      continue
    endif
    times[k] = h.time  
      
    pos = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='pos')
    ids = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='ids')
    
    ; make tracer selection on first snapshot
    if (snap eq snapRange[0]) then begin      

      ; make selection
      dists = reform(abs(sqrt( (pos[0,*]-targetPos[0])*(pos[0,*]-targetPos[0]) + $
                               (pos[1,*]-targetPos[1])*(pos[1,*]-targetPos[1]) + $
                               (pos[2,*]-targetPos[2])*(pos[2,*]-targetPos[2]) )))

      ; NUMTRACERS = select numTracers nearest to each targetPos
      if keyword_set(numTracers) then begin
       
        w = (sort(dists))[0:numTracers-1]
        
        ; arrays
        trPos      = dblarr(nSnaps,numTracers,3)
        idTargets  = lonarr(numTracers)
        
        ; find IDs of targets
        idTargets = ids[w]
        
        print,'Minimum and maximum of located tracers: ',min(dists[w]),max(dists[w])
      endif
      
      ; MAXDIST = select tracers near (<distDegenTol) targetPos
      if keyword_set(maxDist) then begin

        w = where(dists le maxDist,count)
        
        if (count eq 0) then begin
          print,'Error: No tracers within maxDist of specified targetPos.'
          stop
        endif else begin
          print,'Found ['+str(count)+'] tracers near targetPos.'
        endelse
        
        ; arrays
        numTracers = count
        trPos      = dblarr(nSnaps,numTracers,3)
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
    endif else begin
      ; locate IDs using where loop
      for j=0,numTracers-1 do begin
        ind = where(ids eq idTargets[j])
        
        if keyword_set(verbose) then if (snap mod verbose eq 0) then $
          print,' ['+str(j)+'] refound id at ind='+str(ind)
    
        trPos[k,j,0] = pos[0,ind]
        trPos[k,j,1] = pos[1,ind]
        trPos[k,j,2] = pos[2,ind]
      endfor
    endelse
      
    k += 1
    
  endfor ;snap
  
  w = where(times ne 0)

  r = {trPos:trPos[w,*,*],idTargets:idTargets,times:times[w]} ;trDist
  return, r

end

; cosmoTracerTrajPretty(): make interesting image

pro cosmoTracerTrajPretty

  ; config
  res = 128
  run = 'dev.tracer'
  redshift = 2.0
  
  numTracers = 2000
  
  ; load group catalog
  sP    = simParams(res=res,run=run,redshift=redshift)
  units = getUnits()

  h  = loadSnapshotHeader(sP.simPath, snapNum=sP.snap)
  gc = loadGroupCat(sP.simPath,sP.snap,/verbose)
  
  haloIDs = [1]
  
  ; paths and snapshots  
  snapRange = [sP.snap,sP.groupCatRange[0],1]
  plotBase  = 'cosmoPretty_'+run+'_'
  
  ; start plot
  start_PS, sP.plotPath + plotBase + string(snapRange[1],format='(I3.3)')+'_'+$
            string(snapRange[0],format='(I3.3)')+'.eps', xs=7.0, ys=7.0

  ;xrange = [0,sP.boxSize]
  ;yrange = [0,sP.boxSize]
  boxsize = 500.0
  axes = [1,2]
  xrange = [gc.groupPos[axes[0],1]-boxsize,gc.groupPos[axes[0],1]+boxsize]
  yrange = [gc.groupPos[axes[1],1]-boxsize,gc.groupPos[axes[1],1]+boxsize]

  fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=4,ystyle=4,charsize=!p.charsize-0.5,$
       xtitle="",ytitle="",xtickname=replicate(' ',10),ytickname=replicate(' ',10),$
       xticklen=0.00001,yticklen=0.00001,position=[0.0,0.0,1.0,1.0]
         
  foreach haloID,haloIDs,k do begin
    targetPos = gc.groupPos[*,haloID]
    
    ; get tracer positions
    saveFilename = sP.derivPath + 'trTraj_B_'+sP.savPrefix + str(sP.res) + "_" + $
              string(snapRange[1],format='(I3.3)')+'_'+string(snapRange[0],format='(I3.3)')+'_'+$
              str(haloID)+'.sav'
    
    if file_test(saveFilename) then begin
      restore,saveFilename,/verbose
    endif else begin
      tp = getCosmoTracerPos(snapPath=sP.simPath, snapRange=snapRange, targetPos=targetPos, $
                             numTracers=numTracers, verbose=1)
      save,tp,filename=saveFilename
    endelse
    
    ; add to plot
    for i=0,n_elements(tp.idTargets)-1 do begin
      xPts = tp.trPos[*,i,axes[0]]
      yPts = tp.trPos[*,i,axes[1]]
      
      fsc_plot, xPts, yPts, line=0, color=fsc_color('black'), thick=!p.thick-4.0, /overplot
    endfor
    
  endforeach
  
  ; end plot
  end_PS, pngResize=50, /deletePS
  
end

; cosmoTracerTraj(): plot a few individual tracer trajectories over contour of the disk gas density

pro cosmoTracerTraj
 
  ; config
  res = 128
  run = 'dev.tracer.nonrad'
  redshift = 1.0
  
  maxDist = 20.0 ;ckpc
  
  ; load group catalog
  sP    = simParams(res=res,run=run,redshift=redshift)
  units = getUnits()

  h  = loadSnapshotHeader(sP.simPath, snapNum=sP.snap)
  gc = loadGroupCat(sP.simPath,sP.snap,/verbose)
  
  ; select target position
  ;targetPos = [2614.26,123.019,3410.30] ;ckpc

  haloID = 0
  targetPos = gc.groupPos[*,haloID]
 
  snapRange = [sP.snap,sP.groupCatRange[0],1]
  
  ; get tracer positions
  saveFilename = sP.derivPath + 'trTraj_'+sP.savPrefix + str(sP.res) + ".haloID="+str(haloID)+"." + $
            string(snapRange[1],format='(I3.3)')+'-'+string(snapRange[0],format='(I3.3)')+$
            "."+str(maxDist*100)+'.sav'
  
  if file_test(saveFilename) then begin
    restore,saveFilename,/verbose
  endif else begin
    tp = getCosmoTracerPos(snapPath=sP.simPath, snapRange=snapRange, targetPos=targetPos, $
                           maxDist=maxDist, verbose=1)
    save,tp,filename=saveFilename
  endelse
           
  ; cKpc -> cMpc
  tp.trPos /= 1000.0
  targetPos /= 1000.0
            
  ; plots
  start_PS, sP.plotPath + sP.savPrefix + str(res) + '.trTraj_'+string(snapRange[1],format='(I3.3)')+'_'+$
            string(snapRange[0],format='(I3.3)')+'.eps'
  
    fsc_text,0.5,0.8,"tracer trajectory - cosmo "+run+" (6.0 < z < "+string(redshift,format='(f3.1)')+")",$
             /normal,alignment=0.5
    
    ; PLOT ONE - (x,y)
    boxsize = 4.2 ;cMpc
    ;xrange = minmax(tp.trPos[*,*,0])+[-boxsize,boxsize] > [0,0]
    ;yrange = minmax(tp.trPos[*,*,1])+[-boxsize,boxsize] > [0,0]
    ;zrange = minmax(tp.trPos[*,*,2])+[-boxsize,boxsize] > [0,0]
    
    xrange = [targetPos[0]-boxsize/2.0,targetPos[0]+boxsize/2.0] > [0,0] < [sP.boxSize,sP.boxSize]
    yrange = [targetPos[1]-boxsize/2.0,targetPos[1]+boxsize/2.0] > [0,0] < [sP.boxSize,sP.boxSize]
    zrange = [targetPos[2]-boxsize/2.0,targetPos[2]+boxsize/2.0] > [0,0] < [sP.boxSize,sP.boxSize]
              
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,charsize=!p.charsize-0.5,$
         xtitle="x [cMpc]",ytitle="y [cMpc]",position=[0.1,0.2,0.45,0.75],/noerase

    plotsym,0,/fill

    ; overplot tracer trajectories
    for i=0,n_elements(tp.idTargets)-1 do begin
      xPts = tp.trPos[*,i,0]
      yPts = tp.trPos[*,i,1]
      
      ;psym=-8,symsize=0.3
      fsc_plot, xPts, yPts, line=0, color=getColor(i), thick=!p.thick-3.0, /overplot
    endfor
    
    ; underplot circles at targetPos
    tvcircle,(xrange[1]-xrange[0])/80.0,targetPos[0],targetPos[1],color=fsc_color('black'),/data,/fill
    
    ; PLOT TWO (x,z)
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=zrange,xstyle=1,ystyle=1,charsize=!p.charsize-0.5,$
         xtitle="x [cMpc]",ytitle="z [cMpc]",position=[0.55,0.2,0.95,0.75],/noerase
         
    for i=0,n_elements(tp.idTargets)-1 do begin
      xPts = tp.trPos[*,i,0]
      zPts = tp.trPos[*,i,2]
      
      ;psym=-8,symsize=0.3
      fsc_plot, xPts, zPts, line=0, color=getColor(i), thick=!p.thick-3.0, /overplot
      
      ; legend
      ;zStrs = "z="+string(1.0/tp.times-1.0,format='(f4.1)')
      ;legend,zStrs,textcolors=['black'],/right,margin=0.25,charsize=!p.charsize-0.5,box=0
    endfor
  
    ; underplot circles at targetPos
    tvcircle,(xrange[1]-xrange[0])/80.0,targetPos[0],targetPos[2],color=fsc_color('black'),/data,/fill
  
  end_PS, pngResize=50, /deletePS
stop
end

; cosmoCompMassFunctions(): compare gas, tracer, and DM FoF mass functions

pro cosmoCompMassFunctions

  ; config
  res      = 256
  run      = 'tracer'
  redshift = 3.0
  
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
    xrange = [2e7,max(hm_dm)]
    yrange = [1e-4,1e0]
    
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
    yrange = [0.65,1.35]
    fsc_plot,[0],[0],/nodata,/noerase,xrange=xrange,yrange=yrange,/xs,/ys,/xlog,$
             xtitle="mass [h"+textoidl("^{-1}")+" M"+textoidl("_{sun}")+"]",$
             ytitle="ratio",ytickv=[0.8,1.0,1.2],yticks=2,$
             position=[0.18,0.15,0.9,0.35]
             
    ; just interpolate both onto a set of masses then compare
    nbins = 100
    res_pts = 10.0^( findgen(nbins+1)/nbins * (11.5-7.5) + 7.5 )
    gas_res = interpol(y_gas,hm_gas,res_pts)
    tr_res  = interpol(y_tr,hm_tr,res_pts)
    
    ; plot
    fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[1.2,1.2],line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[0.8,0.8],line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,res_pts,tr_res/gas_res,line=0,color=getColor(3),/overplot
    
    ; do the same for the baryonic MF instead of just gas
    res_pts = 10.0^( findgen(nbins+1)/nbins * (11.5-7.5) + 7.5 )
    tr_res  = interpol(y_tr,hm_tr,res_pts)
    bar_res = interpol(y_bar,hm_bar,res_pts)
    fsc_plot,res_pts,tr_res/bar_res,line=0,color=getColor(8),/overplot
    
    ; legend
    legend,['tr/gas','tr/baryon'],textcolors=getColor([3,8],/name),/right,/top,box=0
    
  end_PS
  
  stop  
  
end

;cosmoCompAxisProfiles(): compare gas, tracer, and DM profiles across a box axis

pro cosmoCompAxisProfiles, redshift=redshift

  if not keyword_set(redshift) then stop
  
  ; config
  resSet   = [128]
  run      = 'dev.tracer.nonrad'
  
  ax = 0 ;x
  range = [5000.0,5500.0] ;ckpc
  
  ; loop over resolutions
  foreach res,resSet,k do begin
  
    ; load snapshot info
    sP = simParams(res=res,run=run,redshift=redshift)
    h  = loadSnapshotHeader(sP.simPath, snapNum=sP.snap)
  
    ; setup binning
    nBins  = 100
    
    rho_dm    = fltarr(nbins)
    rho_gas   = fltarr(nbins)
    ;rho_stars = fltarr(nbins)
    rho_tr    = fltarr(nbins)
    
    xBins   = linspace(range[0],range[1],nBins)
    binSize = (range[1]-range[0])/nBins
    midBins = xBins + binSize/2.0
    
    ; load gas,tr,dm,star positions
    pos_gas   = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='pos')
    pos_tr    = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='tracer',field='pos')
    pos_dm    = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='dm',field='pos')
    ;pos_stars = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='stars',field='pos')
    
    gas_mass   = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='mass')
    ;stars_mass = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='stars',field='mass')
    dm_mass    = h.massTable[1]  
  
    ; do binning
    rho_gas   = hist1d(pos_gas[ax,*],gas_mass,binsize=binSize,min=range[0],max=range[1],obin=obin)
    ;rho_stars = hist1d(pos_stars[ax,*],stars_mass,binsize=binSize,min=0.0,max=h.boxSize,obin=obin)
    rho_dm    = hist1d(pos_dm[ax,*],binsize=binSize,min=range[0],max=range[1],obin=obin)
    rho_tr    = hist1d(pos_tr[ax,*],binsize=binSize,min=range[0],max=range[1],obin=obin)
    
    rho_gas = rho_gas[0:n_elements(rho_gas)-2]
    ;rho_stars = rho_stars[0:n_elements(rho_stars)-2]
    rho_dm = rho_dm[0:n_elements(rho_dm)-2]
    rho_tr = rho_tr[0:n_elements(rho_tr)-2]
    
    ; slab volume normalization and mass->Msun
    vol = h.boxSize * h.boxSize * binSize ;kpc^3
    
    rho_gas   = rho_gas / vol * 1e10
    ;rho_stars = rho_stars / vol * 1e10
    rho_dm    = rho_dm * dm_mass / vol * 1e10
    rho_tr    = rho_tr * sP.trMassConst / vol * 1e10
    
    ; normalize DM density down for better plotting
    ;rho_dm /= 4.0

    ; multi plot config
    line  = k
  
    ; start plot
    if (k eq 0) then begin
      start_PS,sP.plotPath+sP.savPrefix+str(n_elements(resSet))+'.axisProf.'+str(ax)+$
               '.snap='+str(sP.snap)+'.eps'

      xrange = range
      yrange = [min([rho_gas[where(rho_gas ne 0)],$
                     rho_dm[where(rho_dm ne 0)]])/1.2,$
                     max([rho_dm,rho_gas])*1.2]
      
      fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
           xtitle="x [kpc]",ytitle="density [h"+textoidl("^2")+" M"+textoidl("_{sun}")+$
           " kpc"+textoidl("^{-3}")+"]",$
           title="non-rad z="+string(redshift,format='(f3.1)')+" box density profile",/ylog
      
      ; plot gas dens
      fsc_plot,midBins,rho_gas,psym=-8,symsize=0.7,/overplot,color=getColor(1)
    endif else begin
      fsc_plot,midBins,rho_gas,line=line,/overplot,color=getColor(1)
    endelse
    
    ; plot other densities
    fsc_plot,midBins,rho_dm,line=line,/overplot,color=getColor(2)
    fsc_plot,midBins,rho_tr,line=line,/overplot,color=getColor(3)
    ;fsc_plot,midBins,rho_stars,line=0,/overplot,color=getColor(7)
    ;fsc_plot,midBins,rho_gas+rho_stars,line=0,/overplot,color=getColor(8)
             
  endforeach ;resSet
  
  ; legend
  strs = ['gas ','dm ',' tracer ','gas ','dm ',' tracer '] + $
         [textoidl('128^3'),textoidl('128^3'),textoidl('128^3'), $
          textoidl('256^3'),textoidl('256^3'),textoidl('256^3')]
  colors = getColor([1,2,3,1,2,3],/name)
  styles = [1,1,1,0,0,0]
  legend, strs, textcolors=colors, linestyle=styles, $
    /right, box=0, margin=0.25, linesize=0.25, charsize=!p.charsize-0.4
  
  end_PS;, pngResize=50
  stop
end

; cosmoStackGroupsRad(): stack radial properties of groups in mass bins and save

function cosmoStackGroupsRad, sP=sP, massBinLog=massBinLog, hIDs=hIDs, stars=stars

  ; load snapshot info and group catalog
  h  = loadSnapshotHeader(sP.simPath, snapNum=sP.snap)
  gc = loadGroupCat(sP.simPath,sP.snap)

  ;nbins = 20
  nbins = 10 ;CUSTOM

  if (n_elements(hIDs) gt 0) then begin
    haloIDs = hIDs ; to prevent overwriting input
    
    ; halo selection (manual)
    ;print,'Using ['+str(n_elements(haloIDs))+'] specified halo IDs.'
    haloRadii = gc.group_r_crit200[haloIDs]
    
    ; setup binning
    logminmax = alog10([1.0,500.0])
    
    saveTag = 'halo='+str(haloIDs[0])+'.'+str(n_elements(haloIDs))+$
              '.snap='+str(sP.snap)+'.nBins='+str(nbins)
    
  endif else begin
    ; halo selection (mass bin)
    massBin = 10.0^massBinLog / 1e10
    haloIDs = where(gc.groupMass ge massBin[0] and gc.groupMass lt massBin[1],count)
    print,'Found ['+str(n_elements(haloIDs))+'] halos in mass bin.'
    if (count eq 0) then stop
    
    ; setup binning
    haloRadii = gc.group_r_crit200[haloIDs]
    
    ;logminmax    = [0.01,2.0] * minmax(haloRadii) > [0.1,100.0]
    ;logminmax[0] = floor(logminmax[0]*10)/10.0
    ;logminmax[1] = ceil(logminmax[1]/100)*100.0
    ;logminmax    = alog10(logminmax)
    logminmax = alog10([1.0,500.0])
    
    saveTag   = 'massbin='+string(massBinLog[0],format='(f4.1)')+"-"+string(massBinLog[1],format='(f4.1)')+$
                '.snap='+str(sP.snap)+'.nBins='+str(nbins)
  endelse

  ; save/restore
  saveFilename = sP.derivPath + sP.savPrefix + str(sP.res) + '.stackRad.' + saveTag + '.sav'
                 
  if file_test(saveFilename) then begin
    restore,saveFilename
  endif else begin 
    ; arrays
    rho_dm    = fltarr(nbins)
    rho_gas   = fltarr(nbins)
    size_gas  = fltarr(nbins)
    num_gas   = fltarr(nbins)
    rho_tr    = fltarr(nbins)
    rho_stars = fltarr(nbins)
    
    ;radBins = [0.0,           logspace(logminmax[0],logminmax[1],nbins)]
    ;midBins = [radBins[0]/2.0,logspace(logminmax[0],logminmax[1],nbins,/mid)]
    ;TODO CUSTOM BINS
    radBins = [0.0,5.0,12.0,25.0,45.0,75.0,100.0,130.0,180.0,300.0,500.0]
    midBins = [2.5,8.5,18.5,35.0,60.0,87.5,115.0,155.0,240.0,400.0]

    ; load gas,tr,dm,star positions
    pos_gas   = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='pos')
    pos_tr    = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='tracer',field='pos')
    pos_dm    = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='dm',field='pos')
    if keyword_set(stars) then $
      pos_stars = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='stars',field='pos')
    
    gas_mass   = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='mass')
    if keyword_set(stars) then $
      stars_mass = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='stars',field='mass')
    dm_mass    = h.massTable[1]
    
    gas_size    = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='vol')
    gas_size    = (gas_size * 3.0 / (4*!pi))^(1.0/3.0) ;cellrad [ckpc]
  
    ; find alternative halo centers via iterative CM fitting
    ;iterCM = groupCenterPosByIterativeCM(sP=sP,gc=gc,haloIDs=haloIDs)
    ; find alternative halo centers via most bound particle ID
    idMBCM = groupCenterPosByMostBoundID(sP=sP)
  
    ; locate halo
    foreach haloID,haloIDs,j do begin
    
      ;haloPos = gc.groupPos[*,haloID] ;use FoF centers
      ;haloPos = iterCM.iterDM[*,j] ;use iterative DM CoM centers
      ;haloPos = gc.subgroupCM[*,gc.groupFirstSub[haloID]] ;use subfind CoM centers
      if (gc.subgroupGrNr[gc.groupFirstSub[haloID]] ne haloID) then begin
        print,'WARNING'
        
        stop
      endif
      haloPos = idMBCM[*,gc.groupFirstSub[haloID]] ;use most bound particle id centers
      
      haloRad = gc.group_r_crit200[haloID]
  
      ; calculate radii and make radial cut
      rad     = periodicRadialDists(haloPos,pos_gas,sP=sP)
      gas_ind = where(rad le 2*10.0^logminmax[1],gas_count)
      rad_gas = rad[gas_ind]
  
      rad    = periodicRadialDists(haloPos,pos_tr,sP=sP)
      tr_ind = where(rad le 2*10.0^logminmax[1],tr_count)
      rad_tr = rad[tr_ind]
  
      rad    = periodicRadialDists(haloPos,pos_dm,sP=sP)
      dm_ind = where(rad le 2*10.0^logminmax[1],dm_count)
      rad_dm = rad[dm_ind]
  
      if keyword_set(stars) then begin
        rad       = periodicRadialDists(haloPos,pos_stars,sP=sP)
        stars_ind = where(rad le 2*10.0^logminmax[1],stars_count)
        rad_stars = rad[stars_ind]
      endif else begin
        stars_count = 0
      endelse
      
      print,'(' + str(haloID) + ') Found ['+str(gas_count)+'] gas ['+str(tr_count)+'] tracer ['+$
            str(dm_count)+'] dm ['+str(stars_count)+'] stars inside cut.'
                       
      ; subselect gas,stars masses
      gas_mass_sub   = gas_mass[gas_ind]
      gas_size_sub   = gas_size[gas_ind]
      if keyword_set(stars) then $
        stars_mass_sub = stars_mass[stars_ind]
      
      ; do binning 
      for i=0,nbins-1 do begin
        ; gas
        w1 = where(rad_gas gt radBins[i] and rad_gas le radBins[i+1],count1)
        if (count1 gt 0) then begin
          rho_gas[i]  += total(gas_mass_sub[w1])
          size_gas[i] += total(gas_size_sub[w1])
          num_gas[i]  += count1
        endif
        
        ; stars
        if keyword_set(stars) then begin
          w2 = where(rad_stars gt radBins[i] and rad_stars le radBins[i+1],count2)
          if (count2 gt 0) then $
            rho_stars[i] += total(stars_mass_sub[w2])
        endif
        
        ; dm
        w3 = where(rad_dm gt radBins[i] and rad_dm le radBins[i+1],count3)
        if (count3 gt 0) then $
          rho_dm[i] += dm_mass * count3
        
        ; tracers
        w4 = where(rad_tr gt radBins[i] and rad_tr le radBins[i+1],count4)
        if (count4 gt 0) then $
          rho_tr[i] += sP.trMassConst * count4
      endfor
    
    endforeach

    ; normalize stacked profiles
    for i=0,nbins-1 do begin
      ; shell volume normalization, average over number of halos, and convert to Msun
      vol = 4*!pi/3 * (radBins[i+1]^3.0 - radBins[i]^3.0) ;kpc^3
      
      rho_gas[i]   = rho_gas[i]   / vol / n_elements(haloIDs) * 1e10
      rho_tr[i]    = rho_tr[i]    / vol / n_elements(haloIDs) * 1e10
      rho_dm[i]    = rho_dm[i]    / vol / n_elements(haloIDs) * 1e10
      rho_stars[i] = rho_stars[i] / vol / n_elements(haloIDs) * 1e10
      
      ; normalize average cell size
      if (num_gas[i] gt 0) then size_gas[i]  = size_gas[i] / num_gas[i]
    endfor
    
    num_gas = num_gas / n_elements(haloIDs)  
    
    r = {rho_gas:rho_gas,size_gas:size_gas,num_gas:num_gas,$
         rho_tr:rho_tr,rho_dm:rho_dm,rho_stars:rho_stars,$
         haloIDs:haloIDs,haloRadii:haloRadii,$
         logminmax:logminmax,radBins:radBins,midBins:midBins,$
         saveTag:saveTag}

    ; save
    save,r,filename=saveFilename
  endelse
  
  return, r

end

; cosmoCompareHaloCenters():

pro cosmoCompareHaloCenters

  res = 256
  run = 'dev.tracer.nonrad'
  redshift = 3.0
  
  sP = simParams(res=res,run=run,redshift=redshift)
  
  gc   = loadGroupCat(sP.simPath,sP.snap)
  gCen = groupCenterPosByMostBoundID(sP=sP)
  
  binSize = 0.25
  min = 0.0
  max = 20.0
  
  num = 20000
  
  ; group catalog consistency checks
  print,'nGroups nSubgroups',gc.nGroupsTot,gc.nSubgroupsTot
  for i=0,gc.nGroupsTot-1 do begin
    fsInd = gc.groupFirstSub[i]
    grNr  = gc.subgroupGrNr[fsInd]
    nSubs = gc.groupNSubs[i]
    if (grNr ne i) then begin & print,'ERROR ',i,fsInd,grNr & stop & endif
  endfor
  
  for i=0,gc.nGroupsTot-1 do begin
    nSubs = gc.groupNSubs[i]
    fsInd = gc.groupFirstSub[i]
    grNrs = gc.subgroupGrNr[fsInd:fsInd+nSubs]
    for j=0,n_elements(grNrs)-1 do begin
      if (grNrs[j] ne i) then stop
    endfor
  endfor
  
  ; centers
  cen_fof  = gc.groupPos[*,0:num]
  cen_sfcm = gc.subgroupCM[*,gc.groupFirstSub[0:num]]
  cen_mb   = gCen[*,gc.groupFirstSub[0:num]]
  
  ; distances
  dist_fof_sfcm = sqrt( (cen_fof[0,*]-cen_sfcm[0,*])^2.0 + $
                        (cen_fof[1,*]-cen_sfcm[1,*])^2.0 + $
                        (cen_fof[2,*]-cen_sfcm[2,*])^2.0 )
  
  dist_fof_mb = sqrt( (cen_fof[0,*]-cen_mb[0,*])^2.0 + $
                      (cen_fof[1,*]-cen_mb[1,*])^2.0 + $
                      (cen_fof[2,*]-cen_mb[2,*])^2.0 )
                      
  dist_sfcm_mb = sqrt( (cen_sfcm[0,*]-cen_mb[0,*])^2.0 + $
                       (cen_sfcm[1,*]-cen_mb[1,*])^2.0 + $
                       (cen_sfcm[2,*]-cen_mb[2,*])^2.0 )
  
  ; plot histos
  hist_fof_sfcm = histogram(reform(dist_fof_sfcm),binsize=binSize,min=min,max=max,loc=loc1)
  hist_fof_mb   = histogram(reform(dist_fof_mb),binsize=binSize,min=min,max=max,loc=loc2)
  hist_sfcm_mb  = histogram(reform(dist_sfcm_mb),binsize=binSize,min=min,max=max,loc=loc3)
  
  start_PS,'hist_fof_sfcm.eps'
    fsc_plot,loc1,hist_fof_sfcm,xtitle="dist [kpc]",ytitle="N"
  end_PS
  
  start_PS,'hist_fof_mb.eps'
    fsc_plot,loc2,hist_fof_mb,xtitle="dist [kpc]",ytitle="N"
  end_PS
  
  start_PS,'hist_sfcm_mb.eps'
    fsc_plot,loc3,hist_sfcm_mb,xtitle="dist [kpc]",ytitle="N"
  end_PS
  
  stop

end

; cosmoCompRadProfiles(): compare gas, tracer, and DM radial profiles of halos

pro cosmoCompRadProfiles, massBinLog=massBinLog, haloIDs=haloIDs, redshift=redshift

  if (n_elements(massBinLog) eq 0 and n_elements(haloIDs) eq 0) then stop
  if not keyword_set(redshift) then stop

  ; config
  resSet = [256,128] ;[256,128]
  run    = 'dev.tracer.nonrad'
  stars  = 0
  
  ; run match
  ;sP1 = simParams(res=resSet[1],run=run,redshift=redshift)
  sP2 = simParams(res=resSet[0],run=run,redshift=redshift)
  ;match = findMatchedHalos(sP1=sP1,sP2=sP2)
  
  ; start plot book
  ;start_PS,sP1.plotPath+sP1.savPrefix+'.book.snap='+str(sP1.snap)+'.radProfiles.ps',eps=0,xs=7,ys=9
  massBinsTag = 'massbin='+string(massBinLog[0],format='(f4.1)')+"-"+string(massBinLog[1],format='(f4.1)')
  start_PS,sP2.plotPath+sP2.savPrefix+'.'+massBinsTag+'.snap='+str(sP2.snap)+'.radProfiles.eps',xs=7,ys=9
  
  ;for i=0,288,2 do begin ;288 matched at z=1
  ;  print,''
  ;  print,i
  ;  hIDs = [ match.matchedInds[match.wMatch[i]], match.wMatch[i] ] ;sP1,sP2 (res[1],res[0]) (128,256)
  
    !p.multi = [0,1,2]
  
    ; plot 1
    ; ------
    foreach res,resSet,k do begin
    ;foreach run,runSet,k do begin
    
      ; load group catalog
      sP = simParams(res=res,run=run,redshift=redshift)
      gc = loadGroupCat(sP.simPath,sP.snap)
    
      ; load radially stacked results
      rs = cosmoStackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=haloIDs,stars=stars)
      
      ; use matched halos between resolutions
      ;rs = cosmoStackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=hIDs[k],stars=stars)
  
      ; check positions
      ;sfInd = gc.groupFirstSub[hIDs[k]]
      ;print,gc.groupPos[*,hIDs[k]]
      ;print,gc.subgroupCM[*,sfInd]
      ;print,'fof sfcm dist, match dist ',$
      ;  periodicRadialDists(gc.groupPos[*,hIDs[k]],gc.subgroupCM[*,sfInd],sP=sP),match.posDiffs[hIDs[1]]
    
      ; multi plot config
      line  = k
    
      ; start plot
      if (k eq 0) then begin
        ;start_PS,sP.plotPath+sP.savPrefix+'.res'+str(n_elements(resSet))+'.'+rs.saveTag+'.radProfiles.eps'
  
        plotStr = "nonRad"
        ;plotStr = str(res)+textoidl('^3')
        
        plotTitle = plotStr+" z="+string(redshift,format='(f3.1)')+" ("+$
                    string(massBinLog[0],format='(f4.1)')+" < log(M) < "+$
                    string(massBinLog[1],format='(f4.1)')+") " + str(n_elements(rs.haloIDs)) + " halos"
  
        ;plotTitle = plotStr+" z="+string(redshift,format='(f3.1)')+" halo ID="+str(haloIDs[0])+$
        ;            " log(M)="+string(alog10(gc.groupmass[haloIDs[0]]*1e10),format='(f5.2)')
              
        ;plotTitle = plotStr+" z="+string(redshift,format='(f3.1)')+" matched"+$
        ;            " log(M)="+string(alog10(gc.groupMass[match.wMatch[i]]*1e10),format='(f5.2)')+" dist="+$
        ;            string(match.posDiffs[match.wMatch[i]],format='(f4.1)')
  
        rs.logminmax = alog10([1.0,500.0])
        xrange = 10.0^rs.logminmax
        ;yrange = [min([rho_gas[where(rho_gas ne 0)],rho_tr,rho_dm,rho_stars[where(rho_stars ne 0)]])/2,$
        ;          max([rho_dm,rho_gas,rho_tr,rho_stars])*2]
        yrange = [min([rs.rho_gas[where(rs.rho_gas ne 0)],rs.rho_dm[where(rs.rho_dm ne 0)]])/2,$
                  max([rs.rho_dm,rs.rho_gas])*2]
         
        fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
             xtitle="r [ckpc]",$
             ytitle="mass density [h"+textoidl("^2")+" M"+textoidl("_{sun}")+$
             " ckpc"+textoidl("^{-3}")+"]",title=plotTitle,/ylog,/xlog
        
        ; r200 lines
        r200 = minmax(gc.group_r_crit200[rs.haloIDs])
  
        fsc_plot,[r200[0],r200[0]],yrange*[1.1,0.98],line=2,color=fsc_color('light gray'),/overplot
        fsc_plot,[r200[1],r200[1]],yrange*[1.1,0.98],line=2,color=fsc_color('light gray'),/overplot
        fsc_text,mean(r200)*0.95,yrange[0]*2,textoidl("r_{200}"),alignment=0.5,$
          charsize=!p.charsize-0.3,color=fsc_color('light gray')
        
        ; plot gas dens
        fsc_plot,rs.midBins,psym=-8,rs.rho_gas,/overplot,color=getColor(1)
      endif else begin
        fsc_plot,rs.midBins,rs.rho_gas,line=line,/overplot,color=getColor(1)
      endelse
      
      ; plot other densities
      fsc_plot,rs.midBins,rs.rho_dm,line=line,/overplot,color=getColor(2)
      fsc_plot,rs.midBins,rs.rho_tr,line=line,/overplot,color=getColor(3)
      
      if (stars eq 1) then begin
        fsc_plot,rs.midBins,rs.rho_stars,line=line,/overplot,color=getColor(7)
        fsc_plot,rs.midBins,rs.rho_gas+rho_stars,line=line,/overplot,color=getColor(8)
      endif
      
      ; softening lines
      soft = sP.gravSoft * [1.0,2.8]
      
     fsc_plot,[soft[0],soft[0]],[yrange[0],yrange[0]*10],line=line,/overplot,color=fsc_color('dark gray')
     fsc_plot,[soft[1],soft[1]],[yrange[0],yrange[0]*10],line=line,/overplot,color=fsc_color('dark gray') 
     
      if (k eq 0) then begin
        fsc_text,soft[0]*0.8,yrange[0]*5,textoidl('\epsilon_{grav}'),color=fsc_color('dark gray'),$
          charsize=!p.charsize-0.5,alignment=0.5
        fsc_text,soft[1]*0.75,yrange[0]*5,textoidl('2.8\epsilon_{grav}'),color=fsc_color('dark gray'),$
          charsize=!p.charsize-0.5,alignment=0.5
      endif
               
    endforeach ;resSet
    
    ; legend (two res one run)
    strs = ['gas ','dm ','tracer ','gas ','dm ','tracer '] + $
           [textoidl('128^3'),textoidl('128^3'),textoidl('128^3'), $
            textoidl('256^3'),textoidl('256^3'),textoidl('256^3')]
    colors = getColor([1,2,3,1,2,3],/name)
    styles = [1,1,1,0,0,0]
    
    ; legend (one res one run)
    ;strs = ['gas ','dm ','tracer ']+$
    ;       [textoidl(str(res)+'^3'),textoidl(str(res)+'^3'),textoidl(str(res)+'^3')]
    ;colors = getColor([1,2,3],/name)
    ;styles = [0,0,0]
    
    ; legend (one res two run)
    ;strs = ['gas non-rad','dm non-rad','tracer non-rad','gas no-PM','dm no-PM','tracer no-PM']
    ;colors = getColor([1,2,3,1,2,3],/name)
    ;styles = [1,1,1,0,0,0]
    
    legend, strs, textcolors=colors, linestyle=styles, $
      /top, /right, box=0, margin=0.25, linesize=0.25, charsize=!p.charsize-0.2
    
    ; plot 2
    ; ------
    foreach res,resSet,k do begin
    
      ; load group catalog
      sP = simParams(res=res,run=run,redshift=redshift)
    
      ; load radially stacked results
      rs = cosmoStackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=haloIDs,stars=stars)
      
      ; load radially stacked results (matched)
      ;rs = cosmoStackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=hIDs[k],stars=stars)
  
      ; reconstruct tracer number density since we didn't save it
      ;nbins  = 20
      ;logminmax = alog10([1.0,500.0])
      ;radBins = [0.0, logspace(logminmax[0],logminmax[1],nbins)]
      nbins = 10
      radBins = [0.0,5.0,12.0,25.0,45.0,75.0,100.0,130.0,180.0,300.0,500.0]
      
      num_tr = fltarr(nbins)
      
      for j=0,nbins-1 do begin
        vol = 4*!pi/3 * (radBins[j+1]^3.0 - radBins[j]^3.0) ;kpc^3
        num_tr[j] = rs.rho_tr[j] * vol * 1 / 1e10 / sP.trMassConst ;=count4
      endfor
      
      ; start plot
      line = k
      
      if (k eq 0) then begin
        rs.logminmax = alog10([1.0,500.0])
        xrange = 10.0^rs.logminmax
        yrange = [0.8,max([rs.num_gas,rs.size_gas])*2]
        
        plotTitle = "gas size and number counts"
        
        fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
             xtitle="r [ckpc]",$
             ytitle="<r"+textoidl("_{gas cell}")+"> [ckpc]  or  N"+textoidl("_{gas}")+"  or  N"+textoidl("_{tr}"),$
             title=plotTitle,/ylog,/xlog
        
        fsc_plot,xrange,[100.0,100.0],line=1,color=fsc_color('light gray'),/overplot
        fsc_plot,xrange,[10.0,10.0],line=1,color=fsc_color('light gray'),/overplot
        
        ; r200 lines
        r200 = minmax(rs.haloRadii)
  
        fsc_plot,[r200[0],r200[0]],yrange*[1.1,0.98],line=2,color=fsc_color('light gray'),/overplot
        fsc_plot,[r200[1],r200[1]],yrange*[1.1,0.98],line=2,color=fsc_color('light gray'),/overplot
        fsc_text,mean(r200)*0.95,yrange[0]*2,textoidl("r_{200}"),alignment=0.5,$
          charsize=!p.charsize-0.3,color=fsc_color('light gray')
        
        ; plot gas dens
        fsc_plot,rs.midBins,rs.num_gas,psym=-8,/overplot,color=getColor(1)
        fsc_plot,rs.midBins,num_tr,psym=-8,/overplot,color=getColor(3)
      endif else begin
        fsc_plot,rs.midBins,rs.num_gas,line=line,/overplot,color=getColor(1)
        fsc_plot,rs.midBins,num_tr,line=line,/overplot,color=getColor(3)
      endelse
      
      ; plot other densities
      fsc_plot,rs.midBins,rs.size_gas,line=line,/overplot,color=getColor(2)
      
      ; softening lines
      soft = sP.gravSoft * [1.0,2.8]
      
     fsc_plot,[soft[0],soft[0]],[yrange[0],yrange[0]*6],line=line,/overplot,color=fsc_color('dark gray')
     fsc_plot,[soft[1],soft[1]],[yrange[0],yrange[0]*6],line=line,/overplot,color=fsc_color('dark gray') 
     
      if (k eq 0) then begin
        fsc_text,soft[0]*0.8,yrange[0]*2,textoidl('\epsilon_{grav}'),color=fsc_color('dark gray'),$
          charsize=!p.charsize-0.5,alignment=0.5
        fsc_text,soft[1]*0.75,yrange[0]*2,textoidl('2.8\epsilon_{grav}'),color=fsc_color('dark gray'),$
          charsize=!p.charsize-0.5,alignment=0.5
      endif
               
      ; legend
      strs = [textoidl('128^3'),textoidl('128^3'),textoidl('128^3'), $
              textoidl('256^3'),textoidl('256^3'),textoidl('256^3')]+$
              [' N'+textoidl("_{gas}"),' <r'+textoidl("_{gas cell}")+'>',' N'+textoidl("_{tr}"),$
              ' N'+textoidl("_{gas}"),' <r'+textoidl("_{gas cell}")+'>',' N'+textoidl("_{tr}")]
      colors = getColor([1,2,3,1,2,3],/name)
      styles = [1,1,1,0,0,0]
      legend, strs, textcolors=colors, linestyle=styles, $
        /top, /left, box=0, margin=0.25, linesize=0.25, charsize=!p.charsize-0.3, spacing=!p.charsize+0.2
    
    endforeach ;resSet
    
    ;end_PS;, pngResize=50
    
  ;endfor ;match
  
  end_PS ;book
stop
  return
end

; cosmoDiffRadProfiles(): calculate difference between gas and tracer radial profiles at several radii
;                         overplot different resolutions and radii as a function of halo mass bin
;                         plot for multiple resolutions

pro cosmoDiffRadProfiles

  resSet = [256,128]
  run    = 'dev.tracer.nonrad'
  stars  = 0
  
  redshifts = [4.0,3.0,2.0,1.0]
  massBins  = [[11.5,12.0],[11.0,11.5],[10.5,11.0],[10.0,10.5]]
  bins      = [4,9,12] ;~2%,10%,25% of r200

  ; arrays
  ratio_trgas = fltarr(n_elements(bins)+1,n_elements(redshifts),$
                       n_elements(resSet),n_elements(massBins[0,*]))
  
  ; load data
  foreach redshift,redshifts,i do begin
    foreach res,resSet,j do begin
      for k=0,n_elements(massBins[0,*])-1 do begin
        massBinLog = massBins[*,k]
        print,'massBinLog: ',massBinLog
  
        ; load group catalog
        sP = simParams(res=res,run=run,redshift=redshift)
        gc = loadGroupCat(sP.simPath,sP.snap)
      
        ; load radially stacked results
        rs = cosmoStackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=haloIDs,stars=stars)
    
        ; save ratios at specified bins
        foreach bin,bins,m do begin
          ratio_trgas[m,i,j,k] = rs.rho_tr[bin] / rs.rho_gas[bin]
        endforeach
        
        ; save mean ratio
        w = where(rs.rho_tr ne 0 and rs.rho_gas ne 0)
        ratio_trgas[m,i,j,k] = mean(rs.rho_tr[w] / rs.rho_gas[w])
  
      endfor ;massBin
    endforeach ;res
  endforeach ;redshift
  
  ; plot once for each massBin
  for k=0,n_elements(massBins[0,*])-1 do begin
    start_PS,sP.plotPath+sP.savPrefix+'.massbin='+string(massBins[0,k],format='(f4.1)')+"-"+$
      string(massBins[1,k],format='(f4.1)')+'.radDiff.eps'
  
    xrange = [6.0,0.0]
    yrange = [0.0,5.0]
    
    plotTitle = "non-rad ("+string(massBins[0,k],format='(f4.1)')+$
      " < log("+textoidl("M_{tot}")+") < "+string(massBins[1,k],format='(f4.1)')+")"
    
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
         xtitle="redshift",ytitle=textoidl("\rho_{tr} / \rho_{gas}"),title=plotTitle
  
    fsc_plot,xrange,[1.0,1.0],line=2,color=fsc_color('light gray'),/overplot
  
    binStrs = ['r/r'+textoidl('_{vir}'),'r/r'+textoidl('_{vir}'),'r/r'+textoidl('_{vir}'),''] + $
              ['~0.02','~0.10','~0.25','mean']
    strings = []
    lines   = []
    colors  = []
  
    foreach res,resSet,j do begin
      for m=0,n_elements(bins) do begin
        ratios = ratio_trgas[m,*,j,k] ;vs redshift
        fsc_plot,redshifts,ratios,psym=-8,line=j,color=getColor(m),/overplot
        
        strings = [strings,str(res)+textoidl('^3')+' '+binStrs[m]]
        lines   = [lines,j]
        colors  = [colors,getColor(m,/name)]
      endfor ;m
    endforeach
    
    ; legend
    legend,strings,linestyle=lines,textcolors=colors,box=0,/left,/top,$
      charsize=!p.charsize-0.2,linesize=0.25
  
    end_PS
  endfor
  
  ;endforeach ;redshift
  stop
  
end

; tracerParentHisto():

pro tracerParentHisto

  ; config
  res = 128 ;not used for non-cosmo
  run = 'dev.tracer.gassphere'
  
  ;redshifts = [5.0,3.0,2.0,1.0] ;redshifts or snap numbers for non-cosmo
  redshifts = [0,1,2,4,6,8,10]
  nNGB = 32
  units = getUnits() ;colors
  
  ; start plot
  sP = simParams(res=res,run=run,redshift=redshifts[0])
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
    
    ;legendStrs   = [legendStrs,'z = '+string(redshift,format='(f4.1)')]
    legendStrs   = [legendStrs,'snap = '+str(redshift)]
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
    
    ;legendStrs   = [legendStrs,'z = '+string(redshift,format='(f4.1)')]
    legendStrs   = [legendStrs,'snap = '+str(redshift)]
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
