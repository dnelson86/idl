; tracersCosmo.pro
; dev for tracer particles related to cosmological boxes
; dnelson feb.2012

; getTracerVelSpatialDens(): wrapper to do tophat density calculation for tracers and save result

function getTracerVelSpatialDens, sP=sP, nNGB=nNGB

  saveFilename = sP.derivPath + sP.savPrefix + str(sP.res) + '.trDens.nNGB=' + str(nNGB) + '.snap=' + $
                 str(sP.snap) + '.sav'
                 
  if file_test(saveFilename) then begin
    restore,saveFilename,/verbose
  endif else begin
    ; load snapshot header for boxSize
    h = loadSnapshotHeader(sP=sP)
    
    ;tr_mass = total(mass_gas) / h.nPartTot[3] ; should be mass gas at t=0
    
    ; load positions and calculate densities via HSMLs
    tr_pos  = loadSnapshotSubset(sP=sP, field='pos', partType='tracerVel')
    tr_dens = estimateDensityTophat(tr_pos,mass=sP.trMassConst,ndims=3,nNGB=nNGB,boxSize=h.boxSize)
    tr_pos  = !NULL
    
    ; save
    save,tr_dens,filename=saveFilename
  endelse
  
  return,tr_dens
  
end

; cosmoTracerVelParents(): return indices (or optionally IDs) of parent gas cells of all tracers

function cosmoTracerVelParents, sP=sP, getInds=getInds, getIDs=getIDs
 
  if (not keyword_set(getInds) and not keyword_set(getIDs)) then stop
  
  ; check for save file
  saveFilename = sP.derivPath + 'trPar_' + sP.savPrefix + str(sP.res) + '_' + str(sP.snap) + '.sav'
  
  if file_test(saveFilename) then begin
    restore, saveFilename
  endif else begin
    ; load
    gas_pos = loadSnapshotSubset(sP=sP, field='pos', partType='gas')
    tr_pos  = loadSnapshotSubset(sP=sP, field='pos', partType='tracerVel')
    
    nTr  = (size(tr_pos))[2]
    nGas = (size(gas_pos))[2]
    
    par_ind  = lonarr(nTr)
    
    ; EXTERNAL C - calculate nearest neighbor to each tracer position
    par_ind = calcNN(gas_pos,tr_pos,boxSize=sP.boxSize,ndims=3)

    ; DEBUG: for a subset, verify nearest neighbors to each tracer position
    nVerify = 100
    indVerify = floor(randomu(seed,nVerify)*nTr)

    for i=0,nVerify-1 do begin
      ; find minimum distance
      dists = periodicDists(tr_pos[*,indVerify[i]],gas_pos,sP=sP)
      w = where(dists eq min(dists),count)
      
      if (count ne 1) then message,'ERROR: More than one mindist in NN debug search!'
      
      ; fail?
      if (w[0] ne par_ind[indVerify[i]]) then begin
        print,'ERROR: Nearest neighbor verification failed.',i,indVerify[i],w[0],par_ind[indVerify[i]]
        stop
      endif
    endfor

    ; save
    save,par_ind,filename=saveFilename
    print,'Saved: ',saveFilename
  endelse
  
  ; if IDs requested, load gas IDs and do crossmatch
  if keyword_set(getIDs) then begin
    gas_ids = loadSnapshotSubset(sP=sP, field='ids', partType='gas')
    
    ; return parent gas ids
    return,gas_ids[par_ind]
    
  endif else begin
    ; return parent gas inds
    return, par_ind
  endelse

end

; cosmoTracerVelChildren(): return indices (or optionally IDs) of child tracer particles of 
;                        specified gas cells (by indices gasInds)
;
; child_counts and child_inds : pass both to skip load and reverse histo (for e.g. loops over more 
; than one gas cell)
; child_counts: pass alone to return number of child tracers per gas cell as an optional output

function cosmoTracerVelChildren, sP=sP, getInds=getInds, getIDs=getIDs, $
                                 gasInds=gasInds, gasIDs=gasIDs, $ ; input: gas cells to search
                                 child_counts=child_counts, child_inds=child_inds

  if (n_elements(gasInds) eq 0 and n_elements(gasIDs) eq 0) then stop
  if (not keyword_set(getInds) and not keyword_set(getIDs)) then stop

  if (n_elements(gasInds) eq 0) then begin
    ; convert input gas IDs into indices
    if n_elements(gasIDs) eq 0 then stop ; IDs required if indices not specified
    gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    idsIndMap = getIDIndexMap(gas_ids,minid=minid)
    gas_ids = !NULL
    
    gasInds = idsIndMap[gasIDs-minid]
    idsIndMap = !NULL
  endif

  ; load if required
  if (n_elements(child_counts) eq 0 or n_elements(child_inds) eq 0) then begin
    ; get tracer parent indices
    tr_par_ind = cosmoTracerVelParents(sP=sP,/getInds)
  
    ; reverse histogram
    child_counts = histogram(tr_par_ind,min=0,rev=child_inds)
  endif  
  
  ; find gas cells with at least one child tracer
  ;numChildren = child_inds[gasInds+1] - child_inds[gasInds] ; equivalent to below
  child_counts = child_counts[gasInds]
  
  w = where(child_counts gt 0,count)
  if (count eq 0) then return, []
  
  ; add all children tracer indices to keeper array
  if total(child_counts,/pres) gt 2e9 then stop ; change tr_inds to lon64arr
  tr_inds = ulonarr(total(child_counts,/pres))
  start = 0L
  
  foreach gasInd,gasInds[w],i do begin
    tr_inds[start:start+child_counts[w[i]]-1] = child_inds[child_inds[gasInd]:child_inds[gasInd+1]-1]
    start += child_counts[w[i]]
  endforeach
  
  ; check for 32 bit long overflow
  if (min(tr_inds) lt 0) then stop
  
  ; reduce memory of return
  if max(child_counts) lt 32767 then child_counts = uint(child_counts)

  ; debug: (slower loop with concat)
  ;tr_inds2 = []
  ;foreach gasInd,gasInds do begin
  ;  ; if number of children is nonzero, add tracer indices to keeper array
  ;  if (child_inds[gasInd+1]-1-child_inds[gasInd] ge 0) then $
  ;    tr_inds2 = [tr_inds2,child_inds[child_inds[gasInd]:child_inds[gasInd+1]-1]]
  ;endforeach
  ;if (total(tr_inds-tr_inds2) ne 0) then stop
  
  ;print,'found ['+str(n_elements(tr_inds))+'] matching tracer children.'
  
  ; if IDs requested, load tracer IDs and do crossmatch
  if keyword_set(getIDs) then begin
    tr_ids = loadSnapshotSubset(sP=sP, field='ids', partType='tracerVel')
    
    ; return children tracer ids
    return,tr_ids[tr_inds]
    
  endif else begin
    ; return children tracer indices
    return, tr_inds
  endelse
end

; cosmoTracerVelMasses(): return array of masses for the tracers by dividing the total mass of each
;                         parent gas cell among all its child tracers (UNUSED)

function cosmoTracerVelMasses, sP=sP
 
  ; check for save file
  saveFilename = sP.derivPath + 'trMass_' + sP.savPrefix + str(sP.res) + '_' + str(sP.snap) + '.sav'
  
  if file_test(saveFilename) then begin
    restore, saveFilename
  endif else begin
    ; load gas mass and tracer parent indices
    gas_mass = loadSnapshotSubset(sP=sP, field='mass', partType='gas')
    tr_par_ind = cosmoTracerVelParents(sP=sP,/getInds)
    
    ; reverse histogram
    child_counts = histogram(tr_par_ind,rev=child_inds)    
    
    ; construct mass array
    tr_mass = fltarr(n_elements(tr_par_ind))
    
    for i=0,n_elements(gas_mass)-1 do begin
      cInds = cosmoTracerVelChildren(sP=sP,gasInds=[i],/getInds,$
                                     child_counts=child_counts,child_inds=child_inds)
      if (n_elements(cInds) gt 0) then $
        tr_mass[cInds] = gas_mass[i] / n_elements(cInds)
    endfor
    
    ; DEBUG: for a subset, verify the tracers have the correct masses
    nVerify = 100
    indVerify = floor(randomu(seed,nVerify)*n_elements(tr_par_ind))

    for i=0,nVerify-1 do begin
      ; find gas parent
      par_ind = tr_par_ind[indVerify[i]]
      
      ; find all tracers with this gas parent
      w = where(tr_par_ind eq par_ind,count)
      if (count eq 0) then message,'ERROR: trVel parent verify failed!'
      
      ; calculate mass using same subdivision scheme
      mass = gas_mass[par_ind] / count
      
      ; fail?
      epsTol = 1e-6
      if (abs(mass-tr_mass[indVerify[i]]) gt epsTol) then begin
        print,'ERROR: Tracer mass verification failed.',i,indVerify[i],mass,tr_mass[indVerify[i]]
        stop
      endif
    endfor

    ; save
    save,tr_mass,filename=saveFilename
    print,'Saved: ',saveFilename
  endelse
  
  ; debug
  ;start_PS,'tr_mass_hist.eps'
  ;  plothist,tr_mass/sP.trMassConst,/auto,/ylog,xrange=[0.0,2.02],yrange=[1e1,5e4],$
  ;    xtitle="tracer mass / gas target mass",ytitle="N",$
  ;    title="M_tr = parent cell M_gas / N_tr children"
  ;end_PS

  return, tr_mass
end

; getCosmoTracerVelPos(): return time series of tracer positions for a given number of tracers nearest
;                      to a given ending position
;
; numTracers=1 : track N tracers closest to targetPos
; maxDist=1    : track all tracers within minDist of targetPos

function getCosmoTracerVelPos, sP=sP, snapRange=snapRange, numTracers=numTracers, $
                               targetPos=targetPos, maxDist=maxDist, verbose=verbose

  useMatch = 1    ; use match for ID location instead of where loop

  ; arrays
  nSnaps  = (snapRange[0]-snapRange[1]+1) / snapRange[2]
  times   = fltarr(nSnaps)
  
  ; verify spacing ok
  nSnaps2 = (float(snapRange[0])-snapRange[1]+1.0) / float(snapRange[2])
  if (nSnaps ne nSnaps2) then message,'Error: Spacing must evenly divide snapshot range.'

  ; loop over requested snapshots
  k = 0
  
  for snap=snapRange[0],snapRange[1],-snapRange[2] do begin
    sP.snap = snap
    if keyword_set(verbose) then if (snap mod verbose eq 0) then $
      print,'snap: ',str(sP.snap)
    
    ; load header and store time
    h = loadSnapshotHeader(sP=sP,/verbose)
    
    ; skip over missing snapshots
    if (n_tags(h) eq 0) then begin
      if keyword_set(verbose) then if (snap mod verbose eq 0) then $
        print,' skipped'
      continue
    endif
    times[k] = h.time  
      
    pos = loadSnapshotSubset(sP=sP,partType='tracerVel',field='pos')
    ids = loadSnapshotSubset(sP=sP,partType='tracerVel',field='ids')
    
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
        
        if (count eq 0) then message,'Error: No tracers within maxDist of specified targetPos.'
        print,'Found ['+str(count)+'] tracers near targetPos.'
        
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
      
      if (count ne numTracers) then message,'Error: Failed to match all targets.'
      
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

; cosmoTracerVelTrajPretty(): make interesting image

pro cosmoTracerVelTrajPretty, numHalos

  ; config
  res = 128
  run = 'dev.tracer.nonrad'
  redshift = 2.0
  
  numTracers = 500
  axes = [1,2]  
  
  ; load group catalog
  sP    = simParams(res=res,run=run,redshift=redshift)
  units = getUnits()

  h  = loadSnapshotHeader(sP=sP)
  gc = loadGroupCat(sP=sP,/verbose)
  
  ; halo selection (manual)
  ;haloIDs = indgen(20)
  
  ; halo selection (~evenly spaced on image plane)
  haloIDs = []
  
  stepSize = h.boxSize / numHalos
  for i=0,numHalos-1 do begin
    for j=0,numHalos-1 do begin
      hDists = sqrt( ((stepSize*i)-gc.groupPos[axes[0],*])^2.0 + ((stepSize*j)-gc.groupPos[axes[1],*])^2.0 )
      ; select quasi random (most massive) from within a box of stepSize x stepSize and avoid edges
      w = where(hDists le stepSize/2.0 and $
                (gc.groupPos[axes[0],*] gt h.boxSize*0.05) and $
                (gc.groupPos[axes[0],*] lt h.boxSize*0.95) and $
                (gc.groupPos[axes[1],*] gt h.boxSize*0.05) and $
                (gc.groupPos[axes[1],*] lt h.boxSize*0.95),count)
      if (count gt 0) then haloIDs = [haloIDs,w[0]]
    endfor
  endfor
  
  ; paths and snapshots  
  snapRange = [sP.snap,sP.groupCatRange[0],1]
  plotBase  = 'cosmoPretty_'+run+'_'
  
  ; start plot
  start_PS, sP.plotPath + plotBase + string(snapRange[1],format='(I3.3)')+'_'+$
            string(snapRange[0],format='(I3.3)')+'_h'+str(n_elements(haloIDs))+'.eps', xs=7.0, ys=7.0

  ; plot range
  xrange = [0,sP.boxSize]
  yrange = [0,sP.boxSize]

  ;boxsize = 500.0
  ;xrange = [gc.groupPos[axes[0],1]-boxsize,gc.groupPos[axes[0],1]+boxsize]
  ;yrange = [gc.groupPos[axes[1],1]-boxsize,gc.groupPos[axes[1],1]+boxsize]

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
      tp = getCosmoTracerVelPos(sP=sP, snapRange=snapRange, targetPos=targetPos, $
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

; cosmoTracerVelTraj(): plot a few individual tracer trajectories over contour of the disk gas density

pro cosmoTracerVelTraj
 
  ; config
  res = 128
  run = 'dev.tracer.nonrad'
  redshift = 1.0
  
  maxDist = 20.0 ;ckpc
  
  ; load group catalog
  sP    = simParams(res=res,run=run,redshift=redshift)
  units = getUnits()

  h  = loadSnapshotHeader(sP=sP)
  gc = loadGroupCat(sP=sP,/verbose)
  
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
    tp = getCosmoTracerVelPos(sP=sP, snapRange=snapRange, targetPos=targetPos, $
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

;cosmoCompAxisProfiles(): compare gas, tracer, and DM profiles across a box axis

pro cosmoCompAxisProfiles
  
  ; config
  resSet   = [128]
  run      = 'dev.tracer.nonrad'
  redshift = 1.0
  
  ax = 0 ;x
  range = [5000.0,5500.0] ;ckpc
  
  ; loop over resolutions
  foreach res,resSet,k do begin
  
    ; load snapshot info
    sP = simParams(res=res,run=run,redshift=redshift)
    ;sP = simParams(res=res,run=run)
    ;sP.snap = 10
    h  = loadSnapshotHeader(sP=sP)
  
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
    pos_gas   = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
    pos_tr    = loadSnapshotSubset(sP=sP,partType='tracerVel',field='pos')
    pos_dm    = loadSnapshotSubset(sP=sP,partType='dm',field='pos')
    ;pos_stars = loadSnapshotSubset(sP=sP,partType='stars',field='pos')
    
    gas_mass   = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
    ;stars_mass = loadSnapshotSubset(sP=sP,partType='stars',field='mass')
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
    rho_dm /= 4.0

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

; tracerVelParentOffsetHisto(): offset of tracers from cell center, split into different histograms
;                            based on the number of child tracers in each gas cell

pro tracerVelParentOffsetHisto

  ; config
  resSet = [256,128]
  run = 'dev.tracer.nonrad'
  
  redshift = 3.0
  
  binSize = 0.05
  minMax  = [0.0,2.0]
  
  foreach res,resSet,k do begin
    ; load simulation parameters
    sP = simParams(res=res,run=run,redshift=redshift)
    
    ; start plot
    line=k
    
    if (k eq 0) then begin
      start_PS, sP.plotPath+sP.savPrefix+'.parOffset.snap='+str(sP.snap)+'.eps'
      fsc_plot,[0],[0],/nodata,xrange=minMax,yrange=[1e-5,1.0],/xs,/ys,$
           xtitle="tracer offset from parent gas center / parent gas cell radius",$
           ytitle="fraction of total number of tracers",$
           title=run+' z='+string(redshift,format='(f3.1)'),/ylog
    endif
    
    ; load tracer parents and reverse histogram
    tr_par_ind = cosmoTracerVelParents(sP=sP,/getInds)
    child_counts = histogram(tr_par_ind,rev=child_inds)    
    
    ; get distance (for each tracer) from parent
    gas_pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
    tr_pos  = loadSnapshotSubset(sP=sP,partType='tracerVel',field='pos')
    
    par_offset = periodicDists(tr_pos,gas_pos[*,tr_par_ind],sP=sP)
    
    gas_pos = !NULL
    tr_pos  = !NULL
    
    ; get gas cell sizes
    gas_size = loadSnapshotSubset(sP=sP,partType='gas',field='vol')
    gas_size = (gas_size * 3.0 / (4*!pi))^(1.0/3.0) ;cellrad [ckpc]
    
    ; normalize by r_cell
    par_offset /= gas_size[tr_par_ind]
    
    gas_size = !NULL
  
    ; all
    h = histogram(par_offset,binsize=binSize,min=minMax[0],max=minMax[1],loc=loc)
    
    h /= float(n_elements(tr_par_ind))
    
    fsc_plot,loc,h,line=line,color=getColor(0),/overplot
  
    ; find tracers in gas cells with N_children=1 and plot histogram
    gasInds = where(child_counts eq 1,count)
    
    cInds = cosmoTracerVelChildren(sP=sP,gasInds=gasInds,/getInds,$
                                   child_counts=child_counts,child_inds=child_inds)
    print,'N=1 ',n_elements(cInds),median(par_offset[cInds])
    
    h = histogram(par_offset[cInds],binsize=binSize,min=minMax[0],max=minMax[1],loc=loc)
    
    h /= float(n_elements(tr_par_ind))
    
    fsc_plot,loc,h,line=line,color=getColor(1),/overplot
    
    ; find tracers in gas cells with N_children>1<5 and plot histogram
    gasInds = where(child_counts gt 1 and child_counts lt 5,count)
    
    cInds = cosmoTracerVelChildren(sP=sP,gasInds=gasInds,/getInds,$
                                   child_counts=child_counts,child_inds=child_inds)
    print,'N>1<5 ',n_elements(cInds),median(par_offset[cInds])
    
    h = histogram(par_offset[cInds],binsize=binSize,min=minMax[0],max=minMax[1],loc=loc)
    
    h /= float(n_elements(tr_par_ind))
    
    fsc_plot,loc,h,line=line,color=getColor(2),/overplot
    
    ; find tracers in gas cells with N_children>1<5 and plot histogram
    gasInds = where(child_counts ge 5,count)
    
    cInds = cosmoTracerVelChildren(sP=sP,gasInds=gasInds,/getInds,$
                                   child_counts=child_counts,child_inds=child_inds)
    print,'N>=5 ',n_elements(cInds),median(par_offset[cInds])
    
    h = histogram(par_offset[cInds],binsize=binSize,min=minMax[0],max=minMax[1],loc=loc)
    
    h /= float(n_elements(tr_par_ind))
    
    fsc_plot,loc,h,line=line,color=getColor(3),/overplot
    
  endforeach ;res
  
  ; legend
  strs = ['all ',$
          'gas with N'+textoidl("_{tr}=1")+' children ',$
          'gas with N'+textoidl("_{tr}>1<5")+' children ',$
          'gas with N'+textoidl("_{tr}>=5")+' children ',$
          'all ',$
          'N'+textoidl("_{tr}=1")+' ',$
          'N'+textoidl("_{tr}>1<5")+' ',$
          'N'+textoidl("_{tr}>=5")+' '] + $
         [textoidl('128^3'),textoidl('128^3'),textoidl('128^3'),textoidl('128^3'), $
          textoidl('256^3'),textoidl('256^3'),textoidl('256^3'),textoidl('256^3')]
  colors = getColor([0,1,2,3,0,1,2,3],/name)
  styles = [1,1,1,1,0,0,0,0]
  
  legend, strs, textcolors=colors, linestyle=styles, $
    /top, /right, box=0, margin=0.25, linesize=0.25, charsize=!p.charsize-0.2
      
  ; end plot
  end_PS
stop
end

; tracerVelParentHisto(): number of tracers per gas cell

pro tracerVelParentHisto

  ; config
  res = 128 ;not used for non-cosmo
  run = 'dev.tracer.nonrad'
  
  redshifts = [5.0,3.0,2.0,1.0] ;redshifts or snap numbers for non-cosmo
  ;redshifts = [0,1,2,4,6,8,10]
  
  ; start plot
  sP = simParams(res=res,run=run,redshift=redshifts[0])
  start_PS, sP.plotPath+sP.savPrefix+str(res)+'.parHisto.eps'
  
  xrange = [-1,20]
  yrange = [1e1,4e6]
  
  fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
       xtitle="number of tracers in a gas cell",ytitle="N gas cells",$
       title=str(res)+textoidl("^3")+" "+run,/ylog  
  
  legendColors = []
  legendStrs   = []
  
  foreach redshift,redshifts,j do begin
  
    sP = simParams(res=res,run=run,redshift=redshift)
    
    ; load
    h = loadSnapshotHeader(sP=sP)
  
    ; load tracer parents and sort
    tr_par_ind = cosmoTracerVelParents(sP=sP, /getInds)
    
    par_histo = histogram(tr_par_ind)
    par_histo = histogram(par_histo,loc=loc)
    
    ; maximum number of parents
    print,'z='+str(redshift)+' max number of parents: ',max(loc)

    ; overplot
    plotsym,0,/fill
    fsc_plot,loc,par_histo,psym=-8,thick=!p.thick+1.0,color=getColor(j),/overplot
    
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
  start_PS, sP.plotPath+sP.savPrefix+str(res)+'.parHisto2.eps'
  
  xrange = [1,100]
  yrange = [1,4e6]
  
  fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
       xtitle="number of tracers in a gas cell",ytitle="N gas cells",$
       title=str(res)+textoidl("^3")+" "+run,/ylog,/xlog
  
  legendColors = []
  legendStrs   = []
  
  foreach redshift,redshifts,j do begin
  
    sP    = simParams(res=res,run=run,redshift=redshift)
    
    ; load
    h = loadSnapshotHeader(sP=sP)
  
    ; load tracer parents and sort
    tr_par_ind = cosmoTracerVelParents(sP=sP, /getInds)
    
    par_histo = histogram(tr_par_ind)
    par_histo = histogram(par_histo,loc=loc)
    
    ; maximum number of parents
    print,'z='+str(redshift)+' max number of parents: ',max(loc)

    ; overplot
    fsc_plot,loc,par_histo,line=0,thick=!p.thick+1.0,color=getColor(j),/overplot
    
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

; cosmoTracerVel_GasDensComp(): compare spatially estimated tracer density to parent cell gas density

pro cosmoTracerVel_GasDensComp

  ; config
  res = 128
  redshift = 1.0
  run = 'dev.tracer.nocomov'
  
  nNGB = 32
  
  sP    = simParams(res=res,run=run) ;,redshift=redshift
  sP.snap = 10
  units = getUnits() ;colors

  ; load
  h = loadSnapshotHeader(sP=sP)
  
  ; get tracer and gas densities
  tr_dens  = getTracerVelSpatialDens(sP=sP, nNGB=nNGB)
  gas_dens = loadSnapshotSubset(sP=sP,partType='gas',field='density',/verbose)
  
  ; load tracer parents and gas densities
  tr_par_ind = cosmoTracerVelParents(sP=sP, /getInds)
  
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

; cosmoTracerVel_GasDensCompNGB(): compare spatially estimated tracer density to parent cell gas density
;                              as a function of nNGB used to calculate tracer densities

pro cosmoTracerVel_GasDensCompNGB

  ; config
  res = 128
  run = 'dev.tracer.nonrad'
  
  redshift = 3.0
  
  nNGBs = [32,64,128,256]
  
  binsize = 0.02
  
  sP    = simParams(res=res,run=run,redshift=redshift)
  units = getUnits() ;colors

  ; load
  h = loadSnapshotHeader(sP=sP)
  
  ; load tracer parents and gas densities
  tr_par_ind = cosmoTracerVelParents(sP=sP, /getInds)
  
  gas_dens = loadSnapshotSubset(sP=sP,partType='gas',field='density',/verbose)  
  
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
      tr_dens  = getTracerVelSpatialDens(sP=sP, nNGB=nNGB)
      
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

; cosmoTracerVel_GasDensCompRedshift(): compare spatially estimated tracer density to parent cell gas density
;                                   as a function of redshift

pro cosmoTracerVel_GasDensCompRedshift

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
    h = loadSnapshotHeader(sP=sP)
  
    ; load tracer parents and gas densities
    tr_par_ind = cosmoTracerVelParents(sP=sP, /getInds)
    gas_dens   = loadSnapshotSubset(sP=sP,partType='gas',field='density',/verbose)  
    
    ; reorder gas densities to correspond to tr_dens parents
    gas_dens = gas_dens[tr_par_ind]  
  
    ; estimate tracer densities
    tr_dens  = getTracerVelSpatialDens(sP=sP, nNGB=nNGB)
      
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
