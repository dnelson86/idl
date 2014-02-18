; testing.pro
; anything temporary
; dnelson feb.2014

pro sarahCheck
  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; 910 result:
  ; 0          65
  ;snap 65: Found [3081] tracers in subgroup 1772 (time: 406.11550 sec)
  ; 1          64
  ;snap 64: Found [818] tracers in subgroup 1739 (time: 362.33149 sec)  
  ;[0 to 1] Matched [0] tracers between the two groups.
  
  ; 910:
  ;         0          65
  ;    58.9394
  ;    305.654      30224.5      21226.9
  ;      9174       13420           0           0        6302           1
  ;         1          64
  ;    62.9272
  ;    285.191      30454.6      21424.4
  ;      9847       14551           0           0        5725           1
  
  
  ; config
  res          = 910
  targetSGInds = [1772,1739] ;,1501,8422]
  targetSnaps  = [65,64] ;,63,62]
  
  ; load
  childTrIDs = {}

  
  foreach snap,targetSnaps,k do begin
    print,k,snap
    start_time = systime(/seconds)
    
    sP = simParams(res=res,run='illustris',snap=snap)
    gc = loadGroupCat(sP=sP,/readIDs)
  
    ; find gas IDs in target subgroup
    gasIDs  = gcPIDList(gc=gc, valGCids=[targetSGInds[k]], partType='gas')
    starIDs = gcPIDList(gc=gc, valGCids=[targetSGInds[k]], partType='stars')
    
    ; find child tracers of those gas cells, and save them
    local_childTrIDs_gas   = cosmoTracerChildren(sP=sP,gasIDs=gasIDs,/getIDs)
    local_childTrIDs_stars = cosmoTracerChildren(sP=sP,starIDs=starIDs,/getIDs)
    
    childTrIDs = mod_struct( childTrIDs, 'snap'+str(snap), {gas:local_childTrIDs_gas, stars:local_childTrIDs_stars} )
  
    ; print number tracers found and timing up to this point
    print,"snap "+str(snap)+": Found ["+str(n_elements(local_childTrIDs_gas))+" "+str(n_elements(local_childTrIDs_stars))+$
          "] tracers in subgroup "+str(targetSGInds[k])+" (time: "+str(systime(/seconds)-start_time)+" sec)"
  
  endforeach
  
  ; compare tracer IDs from first snapshot previous snapshots
  for i=0,n_elements(targetSGInds)-2 do begin
    calcMatch, childTrIDs.(i).gas, childTrIDs.(i+1).gas, inds0, inds1, count=countMatch
    print,"["+str(i)+" gas "+str(i+1)+" gas] Matched ["+str(countMatch)+"] tracers between the two groups."
    
    calcMatch, childTrIDs.(i).gas, childTrIDs.(i+1).stars, inds0, inds1, count=countMatch
    print,"["+str(i)+" gas "+str(i+1)+" str] Matched ["+str(countMatch)+"] tracers between the two groups."

    calcMatch, childTrIDs.(i).stars, childTrIDs.(i+1).gas, inds0, inds1, count=countMatch
    print,"["+str(i)+" str "+str(i+1)+" gas] Matched ["+str(countMatch)+"] tracers between the two groups."   

    calcMatch, childTrIDs.(i).stars, childTrIDs.(i+1).stars, inds0, inds1, count=countMatch
    print,"["+str(i)+" str "+str(i+1)+" str] Matched ["+str(countMatch)+"] tracers between the two groups."    
  endfor
  
  ; 
  
  stop

end