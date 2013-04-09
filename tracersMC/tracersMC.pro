; tracersMC.pro
; dev for MC tracer particles (spherically symmetric setups)
; dnelson nov.2012

; checkSnapshotIntegrity(): check the tracer output in a snapshot makes sense

pro checkSnapshotIntegrity

  ; config
  sP = simParams(res=256,run='debora_test',redshift=2.0)
  
  snaps = lindgen(28)
  subBox = 0 ; subBox snapshot numbering or main snapshot numbering

  foreach snap,snaps do begin
    ; load
    sP.snap = snap
    print,'loading ['+str(sP.snap)+']...'
    
    h = loadSnapshotHeader(sP=sP,subBox=subBox)
    
    ; load parent IDs
    gas_ids    = loadSnapshotSubset(sP=sP,partType='gas',field='ids',subBox=subBox)
    
    star_ids = []
    if (h.nPartTot[4] gt 0) then $
      star_ids   = loadSnapshotSubset(sP=sP,partType='stars',field='ids',subBox=subBox)
      
    bh_ids = []
    if (h.nPartTot[5] gt 0) then $
      bh_ids   = loadSnapshotSubset(sP=sP,partType='bh',field='ids',subBox=subBox)

    ; get tracer children counts
    gas_numtr  = loadSnapshotSubset(sP=sP,partType='gas',field='numtr',subBox=subBox)
    
    star_numtr = [0]
    if (h.nPartTot[4] gt 0) then $
      star_numtr = loadSnapshotSubset(sP=sP,partType='stars',field='numtr',subBox=subBox)
    
    bh_numtr = [0]
    if (h.nPartTot[5] gt 0) then $
      bh_numtr = loadSnapshotSubset(sP=sP,partType='bh',field='numtr',subBox=subBox)
    
    tr_ids    = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids',subBox=subBox)
    tr_parids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentids',subBox=subBox)
    
    print,' load done, checking...'

    ; check for strange IDs
    w = where(tr_ids eq 0,count_zero)
    if (count_zero ne 0) then stop ;bunch with ID=0 are all attached to a star (must be a bug)
    w = where(tr_parids lt 0 or tr_parids gt max([gas_ids,star_ids]),count_badpar)
    if (count_badpar ne 0) then stop  
    
    ; check tracer IDs are unique
    nunique = nuniq(tr_ids)
    print,nunique,n_elements(tr_ids)
    if (nunique ne n_elements(tr_ids)) then stop
    
    ; check sum of num_tr equals total number of tracers output (gas+stars)
    tot_numtr = total(gas_numtr) + total(star_numtr) + total(bh_numtr)
    print,tot_numtr,n_elements(tr_ids)
    if (tot_numtr ne n_elements(tr_ids)) then stop
    
    ; check number of gas cells with tracers equals number of matched parent IDs vs gas IDs
    tr_parids_uniq = tr_parids[uniq(tr_parids,sort(tr_parids))]
    
    match,gas_ids,tr_parids_uniq,ind_gas,ind_tr,count=count_gas
    w = where(gas_numtr gt 0,count_gas_children)
    if (count_gas ne count_gas_children) then stop
    
    ; check number of stars with tracer children equals number of matched parent IDs vs star IDs
    if (h.nPartTot[4] gt 0) then begin
      match,star_ids,tr_parids_uniq,ind_star,ind_tr,count=count_star
      w = where(star_numtr gt 0,count_star_children)
      if (count_star ne count_star_children) then stop
    endif
    
    ; bh
    if (h.nPartTot[5] gt 0) then begin
      match,bh_ids,tr_parids_uniq,ind_bh,ind_tr,count=count_bh
      w = where(bh_numtr gt 0,count_bh_children)
      if (count_bh ne count_bh_children) then stop
    endif
    
    ; check gas and star IDs don't intersect
    if (h.nPartTot[4] gt 0) then begin
      match,gas_ids,star_ids,ind_gas,ind_star,count=count_int
      if (count_int gt 0) then stop
    endif  
    
    ; check gas/star parents don't collide
    match,gas_ids, tr_parids_uniq,gas_ind, loc_ind_gas, count=count_gas
    
    if (h.nPartTot[4] gt 0) then begin
      match,star_ids,tr_parids_uniq,star_ind,loc_ind_star,count=count_star
      match,loc_ind_gas,loc_ind_star,ind1,ind2,count=count_collide
      if (count_collide ne 0) then stop
    endif
    
    ; check gas/bh parents don't collide
    if (h.nPartTot[5] gt 0) then begin
      match,bh_ids,tr_parids_uniq,bh_ind,loc_ind_bh,count=count_bh
      match,loc_ind_gas,loc_ind_bh,ind1,ind2,count=count_collide
      if (count_collide ne 0) then stop
    endif
    
    ; check star/bh parents don't collide
    if (h.nPartTot[4] gt 0 and h.nPartTot[5] gt 0) then begin
      match,loc_ind_star,loc_ind_bh,ind1,ind2,count=count_collide
      if (count_collide ne 0) then stop
    endif
    
    ; for each gas cell, verify its number of tracer children equals the number of tracers pointing to it
    trhist = histogram(tr_parids,min=0,rev=child_inds)
    
    for i=0,count_gas-1 do begin
      ;if (i mod round(count_gas/10.0) eq 0) then print,float(i)/count_gas*100.0
      gas_cur_id = gas_ids[gas_ind[i]]
      gas_cur_numtr = gas_numtr[gas_ind[i]]
      num_children = trhist[gas_cur_id]
      if (num_children ne gas_cur_numtr) then stop
    endfor
    
    ; for each star particle, verify its number of tracer children equals the number of tracers pointing to it
    if (h.nPartTot[4] gt 0) then begin    
      for i=0,count_star-1 do begin
        ;if (i mod round(count_star/10.0) eq 0) then print,float(i)/count_star*100.0
        star_cur_id = star_ids[star_ind[i]]
        star_cur_numtr = star_numtr[star_ind[i]]
        num_children = trhist[star_cur_id]
        if (num_children ne star_cur_numtr) then stop
      endfor
    endif
    
    ; bh
    if (h.nPartTot[5] gt 0) then begin    
      for i=0,count_bh-1 do begin
        ;if (i mod round(count_star/10.0) eq 0) then print,float(i)/count_star*100.0
        bh_cur_id = bh_ids[bh_ind[i]]
        bh_cur_numtr = bh_numtr[bh_ind[i]]
        num_children = trhist[bh_cur_id]
        if (num_children ne bh_cur_numtr) then stop
      endfor
    endif
    
    print,' passed.'
  endforeach ;snaps
end

; cosmoTracerChildren(): return indices (or optionally IDs) of child tracer particles of 
;                        specified gas cells/stars (by indices gasInds or ids gasIDs or ids starIDs)
; note: for MC tracers should not mix gas and stellar searches, do separately

function cosmoTracerChildren, sP=sP, getInds=getInds, getIDs=getIDs, $
                              gasInds=gasInds, gasIDs=gasIDs, starIDs=starIDs, $ ; input: gas cells/stars to search
                              tr_parids=tr_parids, $ ; optional input, if already loaded
                              child_counts=child_counts ; optional output
                              
  compile_opt idl2, hidden, strictarr, strictarrsubs               
  useExternalLowMem = 0 ; roughly 50% slower, but half the peak memory usage
                      
  if (n_elements(gasInds) eq 0 and n_elements(gasIDs) eq 0 and n_elements(starIDs) eq 0) then message,'Input required.'
  if (not keyword_set(getInds) and not keyword_set(getIDs)) then message,'Output type required.'
  if (n_elements(gasIDs) gt 0 and n_elements(starIDs) gt 0) then message,'Either gas or stars.'
  
  if n_elements(gasIDs) eq 0 and n_elements(starIDs) eq 0 then begin
    ; convert input gas indices into IDs
    if n_elements(gasInds) eq 0 then stop ; gas indices required if IDs not specified (no starInds support)
    gasIDs = loadSnapshotSubset(sP=sP,partType='gas',field='ids',inds=gasInds)
  endif
  
  pt = 'gas'
  
  ; if we input stars, switch particle type and override gasIDs with starIDs
  if n_elements(starIDs) gt 0 then begin
    pt = 'stars'
    gasIDs = starIDs
    starIDs = !NULL
  endif
  
  ; get tracer parent IDs
  if ~keyword_set(tr_parids) then $
    tr_parids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentids')
  
  if useExternalLowMem eq 0 then begin
  
    ; (option 1) use reverse histogram approach
    prevMax = max(tr_parids)
    child_counts = histogram(tr_parids,min=0,max=prevMax+1,rev=child_inds,/L64)

    if min(child_inds) lt 0 then message,'Error: Corrupt RI.'
    if max(tr_parids) ne prevMax then message,'Error: Corrupted histo input.'
  
    ; number of child tracers for each requested parent
    child_counts = child_counts[gasIDs]
  
    ; reduce memory usage if possible
    if max(child_counts) lt 32767 then child_counts = long(child_counts)

    ; find gas cells with at least one child tracer
    w = where(child_counts gt 0,count)
    if (count eq 0) then return, []
  
    ; add all children tracer indices to keeper array
    tr_inds = lon64arr(total(child_counts,/int))
    start = 0LL
  
    foreach gasID,gasIDs[w],i do begin
      tr_inds[start:start+child_counts[w[i]]-1] = child_inds[child_inds[gasID]:child_inds[gasID+1]-1]
      start += child_counts[w[i]]
    endforeach
  
  endif else begin
  
    ; (option 2) use external routine
    tr_inds = calcMatchDupe(gasIDs,tr_parids,dupe_counts=child_counts,count=count)
    
  endelse
  
  ; check for 32 bit long overflow
  if min(tr_inds) lt 0 or min(child_counts) lt 0 then message,'Error: Likely overflow.'
  
  ; DEBUG: sanity check on child counts
  ;gas_ids = loadSnapshotSubset(sP=sP,partType=pt,field='ids')
  ;placeMap = getIDIndexMap(gas_ids,minid=minid)
  ;gas_ids = !NULL
  ;num_tr = loadSnapshotSubset(sP=sP,partType=pt,field='numtr',inds=placeMap[gasIDs-minid])
  ;placeMap = !NULL
  ;if not array_equal(child_counts,num_tr) then message,'Error: Tracer child count mismatch.'
  
  ; DEBUG: sanity check (slower loop with concat)
  ;tr_inds2 = []
  ;foreach gasID,gasIDs do begin
  ;  ; if number of children is nonzero, add tracer indices to keeper array
  ;  if (child_inds[gasID+1]-1-child_inds[gasID] ge 0) then $
  ;    tr_inds2 = [tr_inds2,child_inds[child_inds[gasID]:child_inds[gasID+1]-1]]
  ;endforeach
  ;if not array_equal(tr_inds,tr_inds2) then stop
  
  if keyword_set(getIDs) then begin
    ; if tracer IDs requested, load tracer IDs and return subset
    return, loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids',inds=tr_inds)
  endif else begin
    ; return children tracer indices
    return, tr_inds
  endelse
  
end

; checkTracerDistInGal()

pro checkTracerDistInGal

  sP = simParams(res=512,run='tracer',redshift=2.0)
  gc = galaxyCat(sP=sP)
  
  ; load gas IDs and numtr
  ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
  numtr = loadSnapshotSubset(sP=sP,partType='gas',field='numtr')
  
  match,ids,gc.galaxyIDs,ids_ind,gc_ind,count=countMatch
  
  numtr_gal = numtr[ids_ind]
  
  start_PS,'test.eps'
    plothist,numtr_gal,/auto
  end_PS
  
  hist = histogram(numtr_gal,min=0,max=30,loc=loc)
  
  res = gaussfit(loc,hist,A,nterms=3)
  
  print,'gaussian center: ',A[1]
  print,'gaussian width: ',A[2]
  
  stop

end

; tracerMCParentHisto(): number of tracers per gas cell

pro tracerMCParentHisto

  ; config
  res = 128
  run = 'tracerMC.nonrad'
  f   = '1'
  
  redshifts = [30.0,6.0,3.0,1.0]
  
  sP = simParams(res=res,run=run,redshift=redshifts[0],f=f)
  
  start_PS, sP.plotPath+'_'+str(min(redshifts))+'-'+str(max(redshifts))+'.parHisto.eps'
  
  h = loadSnapshotHeader(sP=sP)
  xrange = [-1,sP.trMCPerCell*5]
  yrange = [1,total(h.nPartTot[0])]
  
  fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
       xtitle="number of tracers in a gas cell",ytitle="N gas cells",title=sP.run
  
  fsc_plot,[sP.trMCPerCell,sP.trMCPerCell],yrange,line=0,color=fsc_color("light gray"),/overplot
  
  legendStrs = []
  
  foreach snap,snaps,j do begin
  
    ; load
    sP.snap = snap
    h      = loadSnapshotHeader(sP=sP)
    num_tr = loadSnapshotSubset(sP=sP,partType='gas',field='numtr')
    
    par_histo = histogram(num_tr,loc=loc)

    print,'t ='+string(h.time,format='(f4.1)')+' Gyr (max number of parents: '+str(max(loc))+')'

    ; overplot
    fsc_plot,loc,par_histo,line=0,color=getColor(j),/overplot
    
    ; normalize by mass?
    mass = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
    num_tr /= (mass/mean(mass))
    
    par_histo = histogram(num_tr,loc=loc)
    
    ; overplot
    fsc_plot,loc,par_histo,line=1,color=getColor(j),/overplot
    
    legendStrs   = [legendStrs,'z = '+string(1/h.time-1,format='(f4.1)')]
    ;legendStrs   = [legendStrs,'t ='+string(h.time,format='(f4.1)')+' Gyr']
    
  endforeach
  
  ; legend
  colors=[getColor(indgen(n_elements(snaps)),/name),'light gray']
  legendStrs = [legendStrs,'(dotted = '+textoidl('N_{tr}\cdot m/m_{avg}')+')']
  legend,legendStrs,textcolors=colors,/right,/top,box=0
  
  ; end plot
  end_PS
  stop
end

; tracerMCSpatialDisp(): plot the spatial dispersion of tracer positions which were initially
;                        coincident (children of the same gas cell) as a function of time

pro tracerMCSpatialDisp

  ; config
  redshiftStart = 30.0
  redshifts     = [6.0,3.0,1.0]
  
  sP = simParams(res=256,run='tracerMC.nonrad',redshift=redshiftStart,f='10')
  
  ; analysis (1)
  if 0 then begin
  ; load starting snapshot, build tracer list
  gas_ids   = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
  gas_numtr = loadSnapshotSubset(sP=sP,partType='gas',field='numtr')

  child_tr_ids = cosmoTracerChildren(sP=sP, /getIDs, gasIDs=gas_ids)

  nGasCells = n_elements(gas_ids)
  gas_ids = !NULL
  
  startSnap = redshiftToSnapnum(redshiftStart,sP=sP)
  
  ; arrays
  medianDispAll = fltarr(nGasCells,n_elements(redshifts))
  meanDispAll   = fltarr(nGasCells,n_elements(redshifts))
  maxDispAll    = fltarr(nGasCells,n_elements(redshifts))

  ; for each subsequent redshift
  foreach redshift,redshifts,k do begin
  
    sP.snap = redshiftToSnapNum(redshift,sP=sP)
  
    ; save/restore
    saveFilename = sP.derivPath + sP.savPrefix + '.' + str(sP.res) + '.spatialDisp.snap=' + $
                   str(startSnap) + '-' + str(sP.snap) + '.sav'
                   
    if file_test(saveFilename) then begin
      restore,saveFilename
      medianDispAll[*,k] = medianDisp
      meanDispAll[*,k]   = meanDisp
      maxDispAll[*,k]    = maxDisp
      continue
    endif
    
    medianDisp = fltarr(nGasCells)
    meanDisp   = fltarr(nGasCells)
    maxDisp    = fltarr(nGasCells)
    
    ; load gas ids and find current tracer children
    gas_ids   = loadSnapshotSubset(sP=sP,partType='gas',field='ids')

    child_tr_inds = cosmoTracerChildren(sP=sP, /getInds, gasIDs=gas_ids)
    gas_ids = !NULL
    
    ; load gas positions and use match to transform to tracer positions
    gas_pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
    
    tr_pos = gas_pos[*,child_tr_inds]
    gas_pos = !NULL
    child_tr_inds = !NULL

    ; relocate all tracers by initial "coincident" groups
    tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
    tr_indmap = getIDIndexMap(tr_ids,minid=tr_minid)
    tr_ids = !NULL
    
    ; loop over each group
    start = 0L
    
    for i=0,nGasCells-1 do begin
      if (i mod round(nGasCells/10.0) eq 0) then print,float(i)/nGasCells*100.0
      
      if (gas_numtr[i] gt 1) then begin
        ; use offset to find current tracer positions
        orig_tr_ids = child_tr_ids[start:start+gas_numtr[i]-1]
        cur_tr_inds = tr_indmap[orig_tr_ids-tr_minid]
        cur_tr_pos  = tr_pos[*,cur_tr_inds]
        
        ; calculate median,mean,max of all pairwise distances between tracers in the group
        ;dists = distance_measure(cur_tr_pos) ; nice but not periodic
        dists = periodicPairwiseDists(cur_tr_pos,sP=sP)
        
        medianDisp[i] = median(dists)
        meanDisp[i]   = mean(dists)
        maxDisp[i]    = max(dists)
       
        start += gas_numtr[i]
      endif
    endfor
    
    ; save
    save,medianDisp,meanDisp,maxDisp,filename=saveFilename
    
    medianDispAll[*,k] = medianDisp
    meanDispAll[*,k]   = meanDisp
    maxDispAll[*,k]    = maxDisp
    
  endforeach
  endif
  ; analysis (2)
  ; load parent IDs (and sort by tracer ID for later) to find fraction in same gas cell
  tr_ids    = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
  tr_parids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentids')
  tr_parids_orig = tr_parids[sort(tr_ids)]
  
  ; for each subsequent redshift
  foreach redshift,redshifts,k do begin
  
    sP.snap = redshiftToSnapNum(redshift,sP=sP)
  
    ; relocate the parents of all the tracers
    tr_ids    = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
    tr_parids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentids')
    tr_parids = tr_parids[sort(tr_ids)]

    ; find matching parents with original parent list
    w = where(tr_parids_orig eq tr_parids,countOrig)
    
    ; sanity checks
    if n_elements(tr_parids) ne n_elements(tr_parids_orig) then stop
    
    ; output
    print,redshift,float(countOrig)/n_elements(tr_parids_orig)
    
  endforeach
  stop
  ; set zero distances to very small
  w = where(medianDispAll eq 0,count)
  if (count gt 0) then medianDispAll[w] = 0.1
  w = where(meanDispAll eq 0,count)
  if (count gt 0) then meanDispAll[w] = 0.1
  w = where(maxDispAll eq 0,count)
  if (count gt 0) then maxDispAll[w] = 0.1  
  
  ; plot (1) - kpc separation
  start_PS, sP.plotPath+sP.run+'_'+str(sP.res)+'.f'+str(sP.trMCPerCell)+'.trDisp.eps'
  
    xrange = [1.0,20000.0] ;kpc
    yrange = [0.9,10.0^ceil(max(alog10(medianDispAll)))]
    bin = 0.04
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/xlog,/ylog,$
         xtitle="Separation [kpc]",ytitle="N",$
         title=sP.run+" "+str(sP.res)+textoidl('^3')+" f="+str(sP.trMCPerCell)+" "+$
         textoidl("z_{init}=")+string(redshiftStart,format='(I2)')
    
    strings = []
    
    foreach redshift,redshifts,k do begin
      ; histogram
      h_median  = histogram(alog10(medianDispAll[*,k]),min=-1.0,bin=bin,loc=loc_median)
      h_mean    = histogram(alog10(meanDispAll[*,k]),min=-1.0,bin=bin,loc=loc_mean)
      h_max     = histogram(alog10(maxDispAll[*,k]),min=-1.0,bin=bin,loc=loc_max)

      ; plot
      fsc_plot,10.0^loc_median,h_median,line=0,/overplot,color=getColor(k)
      fsc_plot,10.0^loc_mean,h_mean,line=1,/overplot,color=getColor(k)
      fsc_plot,10.0^loc_max,h_max,line=2,/overplot,color=getColor(k)
      
      frac = string(h_mean[0] / total(h_mean) * 100.0,format='(I2)')
      
      strings = [strings,'z = '+string(redshift,format='(f3.1)')+' ('+frac+'% at zero)']
    endforeach
    
    ; median,mean,max legend
    legend,['median','mean','max'],textcolors=getColor([0,0,0],/name),linestyle=[0,1,2],$
      /right,/bottom,box=0,margin=0.25,linesize=0.25
    legend,strings,textcolors=getColor(indgen(n_elements(redshifts)),/name),$
      /left,/top,box=0,margin=0.25
  end_PS
  
  ; plot (2) - normalize separation by mean gas cell size at that redshift
  start_PS, sP.plotPath+sP.run+'_'+str(sP.res)+'.f'+str(sP.trMCPerCell)+'.trDisp2.eps'
  
    xrange = [0.05,200.0] ;kpc
    yrange = [0.9,1e5]
    bin = 0.04
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/xlog,/ylog,$
         xtitle="Separation / Mean Gas Cell Size",ytitle="N",$
         title=sP.run+" "+str(sP.res)+textoidl('^3')+" f="+str(sP.trMCPerCell)+" "+$
         textoidl("z_{init}=")+string(redshiftStart,format='(I2)')
    
    strings = []
    
    foreach redshift,redshifts,k do begin
      sP.snap = redshiftToSnapNum(redshift,sP=sP)
      
      ; load gas cell volumes
      gas_size = loadSnapshotSubset(sP=sP,partType='gas',field='vol')
      gas_size = (gas_size * 3.0 / (4*!pi))^(1.0/3.0) ;cellrad [ckpc]
      gas_size = mean(gas_size)
      
      ; histogram
      h_median  = histogram(alog10(medianDispAll[*,k]),min=-1.0,bin=bin,loc=loc_median)
      h_mean    = histogram(alog10(meanDispAll[*,k]),min=-1.0,bin=bin,loc=loc_mean)
      h_max     = histogram(alog10(maxDispAll[*,k]),min=-1.0,bin=bin,loc=loc_max)

      ; plot
      fsc_plot,10.0^loc_median/gas_size,h_median,line=0,/overplot,color=getColor(k)
      fsc_plot,10.0^loc_mean/gas_size,h_mean,line=1,/overplot,color=getColor(k)
      fsc_plot,10.0^loc_max/gas_size,h_max,line=2,/overplot,color=getColor(k)
      
      frac = string(h_mean[0] / total(h_mean) * 100.0,format='(I2)')
      
      strings = [strings,'z = '+string(redshift,format='(f3.1)')+' ('+frac+'% at zero)']
    endforeach
    
    ; median,mean,max legend
    legend,['median','mean','max'],textcolors=getColor([0,0,0],/name),linestyle=[0,1,2],$
      /right,/bottom,box=0,margin=0.25,linesize=0.25
    legend,strings,textcolors=getColor(indgen(n_elements(redshifts)),/name),$
      /left,/top,box=0,margin=0.25
  end_PS
  
stop
end

; tracerMCGasStarBalance(): plot the balance between gas and star tracers with time

pro tracerMCGasStarBalance

  ; config
  workingPath = '/n/home07/dnelson/dev.tracerMC/'
  ;snapPath    = workingPath + 'col2Sph.gasonly.1e4.norot.cool.nosg.SF.GFM.f10.noref/output/'
  ;plotBase    = 'col2Sph.gasonly.1e4.cool.nosg.SF.GFM.f10.noref'
  snapPath    = workingPath + 'cosmobox.128_20Mpc.f1.coolSF.GFM.noref/output/'
  plotBase    = 'cosmobox.128_20Mpc.f1.coolSF.GFM.noref'
  
  trMCPerCell = 1
  
  snaps = indgen(n_elements(file_search(snapPath+'snap_*.hdf5')))
  
  ; arrays
  times        = fltarr(n_elements(snaps))
  num_gastr    = fltarr(n_elements(snaps))
  num_startr   = fltarr(n_elements(snaps))
  tot_gasmass  = fltarr(n_elements(snaps))
  tot_starmass = fltarr(n_elements(snaps))
  tot_sfr      = fltarr(n_elements(snaps))
  
  ; loop over each snapshot
  foreach snap,snaps,j do begin
  
    ; load
    print,snap
    h = loadSnapshotHeader(snapPath,snapNum=snap)
    
    num_tr_gas   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='numtr')
    num_tr_stars = loadSnapshotSubset(snapPath,snapNum=snap,partType='stars',field='numtr')
    mass_gas     = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='mass')
    mass_stars   = loadSnapshotSubset(snapPath,snapNum=snap,partType='stars',field='mass')
    sfr          = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='sfr')
    
    times[j] = h.time
    
    num_gastr[j]    = total(num_tr_gas)
    num_startr[j]   = total(num_tr_stars)
    tot_gasmass[j]  = total(mass_gas)
    tot_starmass[j] = total(mass_stars)
    tot_sfr[j]      = total(sfr)
    
  endforeach  
  
  ; start plot
  start_PS, plotBase+'_'+str(min(snaps))+'-'+str(max(snaps))+'.balance.eps',xs=8,ys=7
  
    ; convert times to redshift?
    times = 1/times-1
    
    ;xrange = minmax(times)
    xrange = [10.0,1.0]
    ;yrange = [1,max(num_gastr+num_startr)*2.0]
    yrange = [0.9,1.1]
    psym = 0
    
    !p.multi = [0,1,2]
    !p.charsize -= 0.4
    
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
         xtitle="Redshift",ytitle="Tracer Counts / Masses",title=plotBase
  
    ;fsc_plot,times,num_gastr,psym=psym,color=getColor(1),/overplot
    ;fsc_plot,times,num_startr,psym=psym,color=getColor(2),/overplot
    ;fsc_plot,times,tot_gasmass*1e2,psym=psym,color=getColor(3),/overplot
    ;fsc_plot,times,tot_starmass*1e2,psym=psym,color=getColor(4),/overplot
    ;fsc_plot,times,tot_sfr*1e2,psym=psym,color=getColor(5),/overplot
    ;fsc_plot,times,(num_gastr+num_startr),psym=psym,color=getColor(0),/overplot
    
    fsc_plot,times,num_gastr/(num_gastr+num_startr),psym=psym,color=getColor(1),/overplot
    fsc_plot,times,1.0+num_startr/(num_gastr+num_startr),psym=psym,color=getColor(2),/overplot
    
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
         xtitle="Redshift",ytitle="Tracer Counts / Masses",title=plotBase    
    
    fsc_plot,times,tot_gasmass/(tot_gasmass+tot_starmass),psym=psym,color=getColor(3),/overplot
    fsc_plot,times,1.0+tot_starmass/(tot_gasmass+tot_starmass),psym=psym,color=getColor(4),/overplot
    
    ; legend
    ;strings = ['N'+textoidl('_{tr}')+' in gas','N'+textoidl('_{tr}')+' in stars',$
    ;           '100x M'+textoidl('_{gas,tot}'),'100x M'+textoidl('_{stars,tot}'),$
    ;           'total instant. SFR','total N'+textoidl('_{tr}')+'']
    strings = [textoidl('N_{tr,gas} / N_{tr,tot}'),textoidl('1 + N_{tr,stars} / N_{tr,tot}'),$
               textoidl('M_{gas} / M_{bary}'),textoidl('1 + M_{stars} / M_{bary}')]
    colors=getColor([1,2,3,4],/name)
    legend,strings,textcolors=colors,/left,/bottom,box=0
  
  ; end plot
  end_PS
  stop
  
end

; plotTracerMCTiming(): plot scalings for tracerMC+tracerVEL

pro plotTracerMCTiming

  ; data (n=1 MPI task, no domain decompositions)
  trVel_trMC_perCell  = [1,     5,     20,    50]
  trVel_trMC_runtimes = ([366.6, 384.0, 404.0, 462.8] + [364.7, 374.7, 403.3, 464.8])/2
  
  noVel_trMC_perCell  = [1,     5,     20,     50]
  noVel_trMC_runtimes = ([346.0, 354.2, 384.9,  445.0] + [345.5, 355.2, 384.2, 446.0])/2
  
  noVel_noMC_runtime  = (342.0 + 343.3)/2
  trVel_noMC_runtime  = (359.6 + 349.8)/2
  
  ; plot
  start_PS, '/n/home07/dnelson/dev.tracerMC/trMC.timings.eps'
    
    xrange = [-2,52]
    yrange = [0.98,1.4]
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
         xtitle="TracerMCPerCell",ytitle="CPU Time / No Tracer CPU Time",$
         title="tracer timings (MPI n=1) evrard 1e4 gas (1e4 trVel)",$
         xticks=3,xtickv=[1,5,20,50],xtickname=['1','5','20','50']
         
    fsc_plot,xrange,[1.0,1.0],line=1,color=fsc_color('light gray'),/overplot
         
    fsc_plot,trVel_trMC_perCell,trVel_trMC_runtimes/noVel_noMC_runtime,psym=-8,color=getColor(1),/overplot
    fsc_plot,noVel_trMC_perCell,noVel_trMC_runtimes/noVel_noMC_runtime,psym=-8,color=getColor(2),/overplot
    fsc_plot,[0.0],[trVel_noMC_runtime]/noVel_noMC_runtime,psym=8,color=getColor(3),/overplot
    fsc_plot,[0.0],[noVel_noMC_runtime]/noVel_noMC_runtime,psym=8,color=getColor(4),/overplot
    
    ; legend
    strings = ['trVel.trMC','noVel.trMC','trVel.noMC','noVel.noMC']
    legend,strings,textcolor=getColor([1,2,3,4],/name),/top,/left,box=0
         
  end_PS
  
  ; data (n=2 MPI tasks, domain decomposition every timestep ~200 times)
  trVel_trMC_perCell  = [1,     5,     20,    50]
  trVel_trMC_runtimes = [474.6, 482.1, 520.6, 610.1] + [476.6, 486.9, 520.13, 597.8]
  
  noVel_trMC_perCell  = [1,     5,     20,     50]
  noVel_trMC_runtimes = [446.0, 456.4, 494.4, 580.4] + [450.6, 458.4, 494.3, 573.4]
  
  noVel_noMC_runtime  = (445.4 + 444.4)/2
  trVel_noMC_runtime  = (465.1 + 467.7)/2
  
  ; plot
  start_PS, '/n/home07/dnelson/dev.tracerMC/trMC.timings2.eps'
    
    xrange = [-2,52]
    yrange = [0.95,3.0]
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
         xtitle="TracerMCPerCell",ytitle="CPU Time / No Tracer CPU Time",$
         title="tracer timings (MPI n=2, domain decomp everyTS)",$
         xticks=3,xtickv=[1,5,20,50],xtickname=['1','5','20','50']
         
    fsc_plot,xrange,[1.0,1.0],line=1,color=fsc_color('light gray'),/overplot
         
    fsc_plot,trVel_trMC_perCell,trVel_trMC_runtimes/noVel_noMC_runtime,psym=-8,color=getColor(1),/overplot
    fsc_plot,noVel_trMC_perCell,noVel_trMC_runtimes/noVel_noMC_runtime,psym=-8,color=getColor(2),/overplot
    fsc_plot,[0.0],[trVel_noMC_runtime]/noVel_noMC_runtime,psym=8,color=getColor(3),/overplot
    fsc_plot,[0.0],[noVel_noMC_runtime]/noVel_noMC_runtime,psym=8,color=getColor(4),/overplot
    
    ; legend
    strings = ['trVel.trMC','noVel.trMC','trVel.noMC','noVel.noMC']
    legend,strings,textcolor=getColor([1,2,3,4],/name),/top,/left,box=0
         
  end_PS

end
