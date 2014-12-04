; testing.pro
; anything temporary
; dnelson dec.2014

pro actualTmaxColdHot
  ; feedback.paper referee question
  foreach redshift,[0.0,2.0,5.0] do begin
    print,'z = ',redshift
    sP = simParams(res=512,run='feedback',redshift=redshift)
    
    ; load
    maxTemp = gcSubsetProp(sP=sP,/maxPastTemp,/accretionTimeSubset,accMode='smooth')
    accTvir = gcSubsetProp(sP=sP,/accTvir,/accretionTimeSubset,accMode='smooth')
    parMass = gcSubsetProp(sP=sP,/parMass,/accretionTimeSubset,accMode='smooth')
    
    ; combine
    maxTemp = [maxTemp.gal,maxTemp.gmem,maxTemp.inter,maxTemp.stars,maxTemp.bhs]
    accTvir = [accTvir.gal,accTvir.gmem,accTvir.inter,accTvir.stars,accTvir.bhs]
    parMass = [parMass.gal,parMass.gmem,parMass.inter,parMass.stars,parMass.bhs]
    parMass = alog10(10.0^parMass / 0.7) ; log msun
    
    ; select parent halo mass
    w = where(parMass ge 11.3 and parMass lt 11.4)
    
    maxTemp = maxTemp[w]
    accTvir = accTvir[w]
    
    ; select hot/cold
    w_cold = where(maxTemp/accTvir le 1.0,count_cold)
    w_hot  = where(maxTemp/accTvir gt 1.0,count_hot)
    
    print,' cold: ',count_cold,' hot: ',count_hot
    print,' mean cold temp: ',mean(maxTemp[w_cold])
    print,' mean hot temp: ',mean(maxTemp[w_hot])
  
  endforeach
  
  stop
  
  
end

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


; checkStarIDs(): make sure SPH star ids turn once into gas IDs moving backwards in time
; note: in Arepo runs spawned (not converted) stars will have new IDs with no progenitor gas cell info

pro checkStarIDs

  sP = simParams(res=128,run='gadget',redshift=2.0)
  
  ; load all star particle IDs
  ids = loadSnapshotSubset(sP=sP,partType='stars',field='ids')
  ids_sort = sort(ids)
  mask = intarr(n_elements(ids))
   
  for m=sP.snap,0,-1 do begin
    sP.snap = m
    
    ; load gas ids and match
    gas_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    match,gas_ids,ids,gas_ind,star_ind,count=count1
    star_ind = star_ind[ids_sort]
    
    if count1 gt 0 then mask[star_ind] += 1
    
    ; load star ids and match
    star_ids = loadSnapshotSubset(sP=sP,partType='stars',field='ids')
    match,star_ids,ids,star_ind_cur,star_ind_orig,count=count2
    
    if count1+count2 ne n_elements(ids) then message,'Did not find all original star IDs.'
    
    ; for those stars that are still stars, make sure we have never seen them as gas
    star_ind_orig = star_ind_orig[ids_sort]
    if max(mask[star_ind_orig]) gt 0 then message,'Error: Flip gas back to star.'
    
    print,m,float(count2)/(count1+count2),float(count1)/(count1+count2)
  endfor

end

; groupCatZeroSFRCheck():

pro groupCatZeroSFRCheck

  sP = simParams(res=512,run='tracer',snap=i)

  for i=sP.groupCatRange[0],sP.groupCatRange[1] do begin
    sP.snap = i
    gc = loadGroupCat(sP=sP,/skipIDs)
    print,i,max(gc.subgroupSFR)
  endfor

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


; checkMusicICsVel

pro checkMusicICsVel

  ; config
  file1 = '/n/home07/dnelson/sims.zooms/128_20Mpc_dmonly/ics.hdf5'
  file2 = '/n/hernquistfs1/mvogelsberger/ComparisonProject/512_20Mpc/ICs_DM/output/ICs.' ;0-63
  
  ; load
  x = h5_parse(file1,/read)
  
  y = loadSnapshotOld(file2+'0','none')
  
  ; pos
  print,minmax(x.parttype1.coordinates._data)
  print,minmax(y.parttype1.pos)
  
  ; ids
  print,minmax(x.parttype1.particleids._data)
  print,minmax(y.parttype1.ids)
  
  ; vel
  print,'vel:'
  print,minmax(x.parttype1.velocities._data)
  print,minmax(y.parttype1.vel)
  
  ; load all old ICs
  print,'loop:'
  vel = fltarr(3,y.header.nPartAll[1])
  count = 0
  
  for i=0,63 do begin
    y = loadSnapshotOld(file2+str(i),'none')
    for j=0,2 do vel[j,count:count+y.header.nPartTot[1]-1] = y.parttype1.vel[j,*]
    count += y.header.nPartTot[1]
    print,i,minmax(y.parttype1.vel)
  endfor
  
  print,minmax(vel)
  
  stop



end

; process16bitimg

pro process16bitimg

  ; config
  fname_in  = "frame_8k180_2000_16bit.png"
  fname_out = "test.png"
  min_val   = 0.02
  max_val   = 0.46
  gamma     = 1.30
  sat_val   = 1.2
  
  ; load
  image = read_png(fname_in)
  
  ; convert into float, clip at [min_val,max_val] and scale into [0.0,1.0]
  image_new = float(image) / 65535.0 ;[0,1]
  image_new = (image_new - min_val) / (max_val-min_val)
  image_new = image_new > 0.0 < 1.0
  
  ; apply gamma scaling
  ; todo
  ; apply saturation
  ; todo
  
  ; convert into [0,255] and write
  image_out = byte(round(image_new * 255.0)) > 0 < 255
  
  write_png, fname_out, image_out
  stop

end
