; testing.pro
; anything temporary
; dnelson feb.2014

pro sarahCheck
  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  res          = 910 ;455,910
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
      message,'Above should fail (check gc.snapOffsets)'
    endif else begin
      all_ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
      calcMatch,all_ids,gasIDs,ind1,ind2,count=countMatch
      if countMatch ne n_elements(gasIDs) then message,'Error'
      
      gcIndsType = ind1[sort(ind2)]
    endelse
    
    ; find child tracers of those gas cells, and save them
    local_childTrIDs2   = cosmoTracerChildren(sP=sP,gasIDs=gasIDs,/getIDs)
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