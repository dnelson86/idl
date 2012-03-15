; tracersMC.pro
; dev for MC tracer particles (spherically symmetric setups)
; dnelson mar.2012

; checkSnapshotIntegrity(): check the tracer output in a snapshot makes sense

pro checkSnapshotIntegrity

  ; config
  sP = simParams(res=128,run='tracerMC.ref',f='1')
  
  snaps = indgen(350)
  subBox = 1 ; subBox snapshot numbering or main snapshot numbering

  foreach snap,snaps do begin
    ; load
    sP.snap = snap
    print,'loading ['+str(sP.snap)+']...
    
    h = loadSnapshotHeader(sP=sP,subBox=subBox)
    
    gas_ids    = loadSnapshotSubset(sP=sP,partType='gas',field='ids',subBox=subBox)
    star_ids   = loadSnapshotSubset(sP=sP,partType='stars',field='ids',subBox=subBox)
    gas_numtr  = loadSnapshotSubset(sP=sP,partType='gas',field='numtr',subBox=subBox)
    star_numtr = loadSnapshotSubset(sP=sP,partType='stars',field='numtr',subBox=subBox)
    
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
    tot_numtr = total(gas_numtr) + total(star_numtr)
    print,tot_numtr,n_elements(tr_ids)
    if (tot_numtr ne n_elements(tr_ids)) then stop
    
    ; check number of gas cells with tracers equals number of matched parent IDs vs gas IDs
    match,gas_ids,tr_parids,ind_gas,ind_tr,count=count_gas
    w = where(gas_numtr gt 0,count_gas_children)
    print,count_gas,count_gas_children
    if (count_gas ne count_gas_children) then stop
    
    ; check number of stars with tracer children equals number of matched parent IDs vs star IDs
    if (h.nPartTot[4] gt 0) then begin
      match,star_ids,tr_parids,ind_star,ind_tr,count=count_star
      w = where(star_numtr gt 0,count_star_children)
      print,count_star,count_star_children
      if (count_star ne count_star_children) then stop
    endif
    
    ; check gas and star IDs don't intersect
    if (h.nPartTot[4] gt 0) then begin
      match,gas_ids,star_ids,ind_gas,ind_star,count=count_int
      if (count_int gt 0) then stop
    endif  
    
    ; check gas/star parents don't collide
    match,gas_ids, tr_parids,gas_ind, loc_ind_gas, count=count_gas
    
    if (h.nPartTot[4] gt 0) then begin
      match,star_ids,tr_parids,star_ind,loc_ind_star,count=count_star
      match,loc_ind_gas,loc_ind_star,ind1,ind2,count=count_collide
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
    
    print,' passed.'
  endforeach ;snaps
end

; cosmoTracerChildren(): return indices (or optionally IDs) of child tracer particles of 
;                        specified gas cells (by indices gasInds or ids gasIDs)

function cosmoTracerChildren, sP=sP, getInds=getInds, getIDs=getIDs, $
                              gasInds=gasInds, gasIDs=gasIDs, $ ; input: gas cells to search
                              child_counts=child_counts ; optional output
                              
  if (n_elements(gasInds) eq 0 and n_elements(gasIDs) eq 0) then stop
  if (not keyword_set(getInds) and not keyword_set(getIDs)) then stop
  
  if (n_elements(gasIDs) eq 0) then begin
    ; convert input gas indices into IDs
    if n_elements(gasInds) eq 0 then stop ; indices required if IDs not specified
    gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')    
    gasIDs = gas_ids[gasInds]
    gas_ids = !NULL
  endif
  
  ; get tracer parent IDs
  tr_parids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentids')
  
  ; reverse histogram
  child_counts = histogram(tr_parids,min=0,rev=child_inds)
  
  ; find gas cells with at least one child tracer
  child_counts = child_counts[gasIDs]

  ; DEBUG: sanity check on child counts
  ;num_tr = loadSnapshotSubset(sP.simPath,snapNum=sP.snap,partType='gas',field='numtr')
  ;num_tr = num_tr[gasInds]
  ;if not array_equal(child_counts,num_tr) then stop

  w = where(child_counts gt 0,count)
  if (count eq 0) then return, []
  
  ; add all children tracer indices to keeper array
  if total(child_counts,/pres) gt 2e9 then stop ; change tr_inds to lon64arr
  tr_inds = ulonarr(total(child_counts,/pres))
  start = 0LL
  
  foreach gasID,gasIDs[w],i do begin
    tr_inds[start:start+child_counts[w[i]]-1] = child_inds[child_inds[gasID]:child_inds[gasID+1]-1]
    start += child_counts[w[i]]
  endforeach
  
  ; check for 32 bit long overflow
  if (min(tr_inds) lt 0) then stop
  
  ; reduce memory of return
  if max(child_counts) lt 32767 then child_counts = uint(child_counts)
  
  ; DEBUG: sanity check (slower loop with concat)
  ;tr_inds2 = []
  ;foreach gasID,gasIDs do begin
  ;  ; if number of children is nonzero, add tracer indices to keeper array
  ;  if (child_inds[gasID+1]-1-child_inds[gasID] ge 0) then $
  ;    tr_inds2 = [tr_inds2,child_inds[child_inds[gasID]:child_inds[gasID+1]-1]]
  ;endforeach
  ;if not array_equal(tr_inds,tr_inds2) then stop
  
  print,'found ['+str(n_elements(tr_inds))+'] matching tracer children.'

  ; if tracer IDs requested, load tracer IDs and do crossmatch
  if keyword_set(getIDs) then begin
    tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
    
    ; return children tracer ids
    return,tr_ids[tr_inds]
    
  endif else begin
    ; return children tracer indices
    return, tr_inds
  endelse
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
