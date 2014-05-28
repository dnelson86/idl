; testing.pro
; anything temporary
; dnelson may.2014

pro sarahCheck
  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  res          = 1820 ;455,910
  run          = 'illustris' ;feedback,illustris
  targetSGInds = [1772,1739] ;,1501,8422]
  targetSnaps  = [65,64] ;,63,62]
  
  ; load
  childTrIDs = {}

  foreach snap,targetSnaps,k do begin
    print,'k snap: ',k,snap
    start_time = systime(/seconds)
    
    sP = simParams(res=res,run=run,snap=snap)
    gc = loadGroupCat(sP=sP,/readIDs)
  
    ; find gas IDs in target subgroup
    gasIDs  = gcPIDList(gc=gc, valGCids=[targetSGInds[k]], partType='gas')
    print,size(gasIDs,/tname)
    
    if run eq 'illustris' then begin
      ; check with direct IDs from group ordered snapshots
      gcIndsType = lindgen(gc.subgrouplentype[0,targetSGInds[k]]) + $
                   gc.snapOffsets.subgroupType[0,targetSGInds[k]]
    
      gasIDsCheck = loadSnapshotSubset(sP=sP,partType='gas',field='ids',inds=gcIndsType)
      print,'ID check passed: ',array_equal(gasIDs,gasIDsCheck)
      ;message,'Above should fail (check gc.snapOffsets)'
    endif else begin
      all_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
      calcMatch,all_ids,gasIDs,ind1,ind2,count=countMatch
      if countMatch ne n_elements(gasIDs) then message,'Error'
      
      gcIndsType = ind1[sort(ind2)]
    endelse
    
    ; find child tracers of those gas cells, and save them
    local_childTrIDs2  = cosmoTracerChildren(sP=sP,gasIDs=gasIDs,/getIDs,/useTwo)
    local_childTrIDs   = cosmoTracerChildren(sP=sP,gasIDs=gasIDs,/getIDs,/useOne)
    print,'one: ',size(local_childTrIDs,/tname)
    print,'two: ',size(local_childTrIDs2,/tname)
    print,'trids equal: ',array_equal(local_childTrIDs,local_childTrIDs2)
    
    numTr = loadSnapshotSubset(sP=sP,partType='gas',field='numTr',inds=gcIndsType)
    print,'numTr comp',total(numTr,/int),n_elements(local_childTrIDs)
    
    childTrIDs = mod_struct( childTrIDs, 'snap'+str(snap), local_childTrIDs )
  
    ; print number tracers found and timing up to this point
    print,"snap "+str(snap)+": Found ["+str(n_elements(local_childTrIDs))+"] tracers (gas parents) in subgroup "+$
          str(targetSGInds[k])+" (time: "+str(systime(/seconds)-start_time)+" sec)"
  
  endforeach
  
  ; compare tracer IDs from first snapshot previous snapshots
  for i=0,n_elements(targetSGInds)-2 do begin
    calcMatch, childTrIDs.(i), childTrIDs.(i+1), inds0, inds1, count=countMatch
    match, childTrIDs.(i), childTrIDs.(i+1), inds0, inds1, count=countMatch2
    print,"["+str(i)+" "+str(i+1)+" ] Matched ["+str(countMatch)+" / "+str(countMatch2)+$
          "] tracers (gas-gas) between the two groups."
  endfor
  
  stop

end


; mcIllustrisCheck():

pro mcIllustrisCheck

  ; config
  ;sP = simParams(res=128,run='feedback',redshift=2.0)
  sP = simParams(res=1820,run='illustris',redshift=2.0)
  
  ; DEBUG
  ;ids_gas = loadSnapshotSubset(sP=sP,partType='gas',field='ids',/verbose)
  ;print,minmax(ids_gas)
  ;gcPIDs = ids_gas[0:100]
  ;calcMatch,gcPIDs,ids_gas,gc_ind,ids_ind,count=countMatch
  ; END DEBUG
  
  gc = galaxyCat(sP=sP)
  print,'ids minmax: ',minmax(gc.ids)
  
  ; run
  print,'runOne:'
  x1 = cosmoTracerChildren( sP=sP, /getInds, gasIDs=gc.ids, child_counts=galcat_cc1, /useOne)
  
  print,'runTwo:'
  x2 = cosmoTracerChildren( sP=sP, /getInds, gasIDs=gc.ids, child_counts=galcat_cc2, /useTwo)
  
  print,'inds equal: ',array_equal(x1,x2)
  print,'cc   equal: ',array_equal(galcat_cc1,galcat_cc2)
  stop
end

; check1820galcat

pro check1820galcat

  sP1 = simParams(res=1820,run='illustris',redshift=7.0)
  sP2 = simParams(res=256,run='feedback',redshift=3.0)
 
  gc1 = galaxyCat(sP=sP1)
  gc2 = galaxyCat(sP=sP2)
  
  rvir1 = galCatParentProperties(sP=sP1, galcat=gc1, /rVirNorm)
  rvir2 = galCatParentProperties(sP=sP2, galcat=gc2, /rVirNorm)
  
  start_PS,'test_rad.eps'
    w = where(gc1.rad gt 0.0 and finite(gc1.rad))
    plothist,alog10( gc1.rad[w] ),/auto,color=cgColor('blue'),xtitle="log (rad)"
    
    w = where(gc2.rad gt 0.0 and finite(gc2.rad))
    plothist,alog10( gc2.rad[w] ),/auto,/overplot,color=cgColor('red')
    
    legend,['illustris','feedback'],textcolor=['blue','red'],/top,/left
  end_PS
  
  start_PS,'test_radnorm.eps'
    w = where(rvir1 gt 0.0 and finite(rvir1))
    plothist,alog10( rvir1[w] ),/auto,color=cgColor('blue'),xtitle="log (rad)"
    
    w = where(rvir2 gt 0.0 and finite(rvir2))
    plothist,alog10( rvir2[w] ),/auto,/overplot,color=cgColor('red')
    
    legend,['illustris','feedback'],textcolor=['blue','red'],/top,/left
  end_PS
  
  start_PS,'test_len1.eps'
    w = where(gc1.len gt 0.0)
    plothist,alog10( gc1.len[w] ),/auto,color=cgColor('blue'),/ylog,xtitle="log (length)"
    
    legend,['illustris','feedback'],textcolor=['blue','red'],/top,/left
  end_PS
  
    start_PS,'test_len2.eps'
    w = where(gc2.len gt 0.0)
    plothist,alog10( gc2.len[w] ),/auto,color=cgColor('red'),/ylog,xtitle="log (length)"
    
    legend,['illustris','feedback'],textcolor=['blue','red'],/top,/left
  end_PS  
  


  stop
end

; check1820b

pro check1820b

  sP = simParams(res=1820,run='illustris',redshift=2.0)
  h = loadSnapshotHeader(sP=sP)
  
  ; load ids of particles in all primary subfind groups
  gc = loadGroupCat(sP=sP,/readIDs)
  gcPIDs = gcPIDList(gc=gc,select='pri',partType='gas')
  
  print,'minmax gcPIDs: ',minmax(gcPIDs)
  print,'num gcPIDs: ',n_elements(gcPIDs)

  gcIDs = gcIDList(gc=gc,select='pri')
  
  print,'total pri gas: ',total( gc.subgroupLenType[partTypeNum('gas'),gcIDs],/int )

  ; load gas ids and match to catalog
  ids_gas = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
  print,'minmax ids_gas: ',minmax(ids_gas)
  print,'num ids_gas: ',n_elements(ids_gas)
  
  ; indexing check
  if 1 then begin
    num_test = h.nPartTot[partTypeNum('gas')]/1000
    test_inds = [l64indgen(num_test)*100,l64indgen(num_test)*500]
  
    ids_test = loadSnapshotSubset(sP=sP,partType='gas',field='ids',inds=test_inds,/verbose)
    print,'first: ',array_equal(ids_test,ids_gas[test_inds])
  
    test_inds = shuffle(test_inds)
  
    ids_test = loadSnapshotSubset(sP=sP,partType='gas',field='ids',inds=test_inds)
    print,'shuffled: ',array_equal(ids_test,ids_gas[test_inds])
  endif
  
  print,'header check: ',n_elements(ids_gas) eq h.nPartTot[ partTypeNum('gas') ]
  
  ; verify matching
  calcMatch,gcPIDs,ids_gas,gc_ind,ids_ind,count=countMatch
  
  ;match,gcPIDs,ids_gas,gc_ind2,ids_ind2,count=countMatch2
  ;print,'match1: ',countMatch eq countMatch2
  ;print,'match2: ',array_equal(gc_ind,gc_ind2)
  ;print,'match3: ',array_equal(ids_ind,ids_ind2)
  ;gc_ind2 = !NULL & ids_ind2 = !NULL
  
  print,'match check1: ',array_equal( gcPIDs[gc_ind], ids_gas[ids_ind] )
  print,'match check2: ',countMatch eq n_elements(gc_ind)
  print,'match check3: ',countMatch eq n_elements(ids_ind)
  print,'match check4: ',min(gc_ind) ge 0 and max(gc_ind) lt n_elements(gcPIDs)
  print,'match check5: ',min(ids_ind) ge 0 and max(ids_ind) lt n_elements(ids_gas)

  ; verify sorting
  ids_ind1 = ids_ind[calcSort(gc_ind)] ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs
  ids_ind2 = ids_ind[sort(gc_ind)]
  
  print,'sort: ',array_equal(ids_ind1,ids_ind2)
  
  ids_ind2 = !NULL
  
  stop

end

; checkAccModeSum

pro checkAccModeSum

  ; config
  sP = simParams(res=128,run='tracer',redshift=2.0)
  accModes = ['smooth','clumpy','stripped']
  
  ; load
  galcat = galaxyCat(sP=sP)
  mt = mergerTreeSubset(sP=sP)
  am = accretionMode(sP=sP)
  at = accretionTimes(sP=sP)

  w1 = where(am gt 0)
  w2 = where(at.accTime[sP.atIndMode,*] ge 0)
  print,array_equal(w1,w2)
  
  ; index check
  tot_inds = 0L
  inds_all = accModeInds(sP=sP,mt=mt,accMode='all')
  for i=0,n_tags(inds_all)-1 do begin
    tot_inds += n_elements(inds_all.(i))
    print,'all',(tag_names(inds_all))[i],minmax(am[inds_all.(i)])
  endfor
  print,tot_inds,n_elements(am)
  
  tot_inds_all_modes = 0L
  foreach accMode,accModes do begin
    tot_inds_mode = 0L
    inds_mode = accModeInds(sP=sP,mt=mt,accMode=accMode)
    for i=0,n_tags(inds_mode)-1 do begin
      tot_inds_mode += n_elements(inds_mode.(i))
      print,accMode,(tag_names(inds_mode))[i],minmax(am[inds_mode.(i)])
    endfor
    tot_inds_all_modes += tot_inds_mode
  endforeach
  
  w=where(am eq 0,count)
  print,tot_inds_all_modes+count,tot_inds
  
  ; check against an atS quantity
  accTvir = gcSubsetProp(sP=sP,/accTvir,/accretionTimeSubset,accMode='all')
  
  tot_quant = 0L
  for i=0,n_tags(accTvir)-1 do tot_quant += n_elements(accTvir.(i))
  
  ; load IDs based on each mode, then match back against galaxy catalog, 
  ; then subset am and make sure they agree
  foreach accMode,[accModes,'all'] do begin
    elemIDs = gcSubsetProp(sP=sP,/elemIDs,/accretionTimeSubset,accMode=accMode)
  
    for i=0,n_tags(elemIDs)-1 do begin
      calcMatch,elemIDs.(i),galcat.trMC_ids,ind1,ind2,count=countMatch
      if countMatch ne n_elements(elemIDs.(i)) then message,'Error'
      am_local = am[ind2]
      print,accMode,(tag_names(elemIDs))[i],minmax(am_local)
    endfor
  endforeach ;accModes
  
  stop

end

; checkGroupOrderedOffsets():

pro checkGroupOrderedOffsets

  ; config
  sP = simParams(res=455,run='illustris',redshift=2.0)
  pt = 'gas'
  groupID = 7000 ;7000, 8000
  
  ; load
  ptNum = partTypeNum(pt)
  gc = loadGroupCat(sP=sP,/skipIDs)
  subID = gc.groupFirstSub[groupID]
  
  print,'Picked subid ['+str(subID)+'] of ('+str(gc.groupNsubs[groupID])+') in this group.'
  
  ; pick one halo, make its index list into the snapshot
  ind_min = gc.snapOffsets.subgroupType[ptNum,subID]
  ind_max = ind_min + gc.subgroupLenType[ptNum,subID] - 1
  
  gr_min = gc.snapOffsets.groupType[ptNum,groupID]
  gr_max = gr_min + gc.groupLenType[ptNum,groupID] - 1
  
  print,'Reading indices: ',ind_min,ind_max
  print,'Group indices: ',gr_min,gr_max
  
  pos = loadSnapshotSubset(sP=sP,partType=pt,field='pos',indRange=[ind_min,ind_max])
  print,'In subhalo:'
  print,minmax(pos[0,*])
  print,minmax(pos[1,*])
  print,minmax(pos[2,*])
  
  pos = loadSnapshotSubset(sP=sP,partType=pt,field='pos',indRange=[gr_min,gr_max])
  print,'In group:'
  print,minmax(pos[0,*])
  print,minmax(pos[1,*])
  print,minmax(pos[2,*])
  
  pos = loadSnapshotSubset(sP=sP,partType=pt,field='pos',indRange=[gr_max+1,gr_max+10])
  print,'Boundary (1):'
  print,minmax(pos[0,*])
  print,minmax(pos[1,*])
  print,minmax(pos[2,*])
  
  pos = loadSnapshotSubset(sP=sP,partType=pt,field='pos',indRange=[gr_min-10,gr_min-1])
  print,'Boundary (2):'
  print,minmax(pos[0,*])
  print,minmax(pos[1,*])
  print,minmax(pos[2,*])
  stop


end


