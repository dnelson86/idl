; old.pro
; unused routines (cut down on number loaded by .r helper)
; dnelson mar.2013

; mymatch(): produces identical output of match.pro, but no concatenated array

pro mymatch, A, B, ind1, ind2, count=count, retUnsorted=retUnsorted

  ; sort A and B if not already
  sort_A = sort(A)
  sort_B = sort(B)
  AA = A[sort_A]
  BB = B[sort_B]
  
  i = 0LL
  j = 0LL
  
  numA = n_elements(AA)
  numB = n_elements(BB)
  
  if numA gt 2e9 or numB gt 2e9 then message,'Error: Move to 64bit.'
  
  ind1 = lonarr(numA)
  ind2 = lonarr(numB)
  
  offset = 0LL
  
  while i lt numA and j lt numB do begin
    if AA[i] eq BB[j] then begin
      if keyword_set(retUnsorted) then begin
        ind1[offset] = i
        ind2[offset] = j
      endif else begin
        ; matched element
        ind1[offset] = sort_A[i]
        ind2[offset] = sort_B[j]
      endelse
      i += 1 & j += 1 & offset += 1
    endif else if AA[i] lt BB[j] then begin
      ; element mismatch, forward in A
      i += 1
    endif else begin
      ; element mismatch, forward in B
      j += 1
    endelse
  endwhile
  
  if offset gt 0 then begin
    ind1 = ind1[0:offset-1]
    ind2 = ind2[0:offset-1]
  endif else begin
    ind1 = -1
    ind2 = -1
  endelse
  
  count = offset
  
end

; mymatchdupe(): produces identical output of match.pro, but no concatenated array

function mymatchdupe, A, B, ind=ind, count=count

  ; sort A and B if not already
  sort_A = sort(A)
  sort_B = sort(B)
  AA = A[sort_A]
  BB = B[sort_B]
  
  i = 0LL
  j = 0LL
  
  numA = n_elements(AA)
  numB = n_elements(BB)
  
  if numA gt 2e9 or numB gt 2e9 then message,'Error: Move to 64bit.'
  
  ind1 = lonarr(numA)
  ind2 = lonarr(numB)-1
  
  offsetTable = lonarr(numA+1)
  
  offset = 0LL
  count = 0
  
  while i lt numA and j lt numB do begin
    if AA[i] eq BB[j] then begin
      ; matched element
      ind1[sort_A[i]] += 1 ; child_counts in unsorted order
      ind2[offset] = sort_B[j]
      
      j += 1
      offset += 1
      count += 1 ; children of this parent
    endif else if AA[i] lt BB[j] then begin
      ; element mismatch, forward in A
      i += 1
      offsetTable[i] = offsetTable[i-1] + count
      count = 0
    endif else begin
      ; element mismatch, forward in B
      j += 1
    endelse
  endwhile
  
  count = offset
  
  ; rearrange ind2 based on the sorting of A (using offsetTable)
  ind3 = lonarr(count)
  
  xx = lonarr(numA)
  xx[sort_A] = lindgen(numA)
  
  print,xx[0:9]
  
  offset = 0LL
  
  for i=0,numA-1 do begin
    ind = xx[i]
    numChildren = ind1[i]
  
    if numChildren eq 0 then continue
    ind3[offset : offset + numChildren - 1] = ind2[offsetTable[ind] : offsetTable[ind] + numChildren - 1]
    offset += numChildren
  endfor
  
  ind = ind1 ; child_counts
  return, ind3 ; tr_inds
  
end


; testMatchBlock()

pro testMatchBlock

  ; generate some data
  num = 50
  A = shuffle([indgen(num)+30,indgen(num)+100])
  B = shuffle([indgen(num)+1000,indgen(num)+80])
  
  ; split B and match halves
  num = n_elements(B)
  
  b1 = b[0:num/2-1]
  b2 = b[num/2:num-1]
  match,a,b1,inda1,indb1,count=count1
  match,a,b2,inda2,indb2,count=count2
  
  count = count1+count2 ; number matched
  
  ; test
  mymatch,a,b1,inda1_test,indb1_test,count=count3,/retUnsorted
  mymatch,a,b2,inda2_test,indb2_test,count=count4,/retUnsorted
  
  inda_test = [inda1_test,inda2_test]
  indb_test = [indb1_test,indb2_test + n_elements(b1)]
  
  ; rearrange inda_test
  inda_fixed1 = (sort(a))[inda1_test]
  inda_fixed2 = (sort(a))[inda2_test]
  
  inda_fixed = [inda_fixed1,inda_fixed2]
  inda_fixed = inda_fixed[ sort(a[inda_fixed]) ]
  
  ; rearrange indb_test by the global sort of b
  indb_fixed1 = (sort(b1))[indb1_test]
  indb_fixed2 = (sort(b2))[indb2_test] + n_elements(b1)
  
  indb_fixed = [indb_fixed1,indb_fixed2] ; correct, just in wrong order
  indb_fixed = indb_fixed[ sort(b[indb_fixed]) ]
  
  ; do it again with a loop, have only full A to start
  ind1 = lonarr(n_elements(a))
  ind2 = lonarr(n_elements(b))
  
  offset = 0L
  b_offset = 0L
  
  for i=0,1 do begin
    ; load b subset
    if i eq 0 then b_subset = b[0:num/2-1]
    if i eq 1 then b_subset = b[num/2:num-1]
    
    ; match
    mymatch,a,b_subset,inda_sub,indb_sub,count=count_sub;,/retUnsorted
    b_offset += n_elements(b_subset)
    
    if count_sub eq 0 then continue
    
    ind1_sub = inda_sub
    ind2_sub = indb_sub + (b_offset-n_elements(b_subset))
    
    ind1[offset : offset+count_sub-1] = ind1_sub
    ind2[offset : offset+count_sub-1] = ind2_sub
    
    offset += count_sub
  endfor
  
  ; take matched subset
  if offset gt 0 then begin
    ind1 = ind1[0:offset-1]
    ind2 = ind2[0:offset-1]
  endif else begin
    ind1 = -1
    ind2 = -1
  endelse
  
  ; now finished, fix ind1 while we have A
  ind1 = ind1[ sort(a[ind1]) ]
  ; delete A and load full B, fix ind2
  ind2 = ind2[ sort(b[ind2]) ]
  
  ; match full for comparison
  match,a,b,Ia,Ib,count=count5

  print,'new:'
  print,array_equal(inda_fixed,Ia)
  print,array_equal(indb_fixed,Ib)
  print,'loop:'
  print,array_equal(ind1,Ia)
  print,array_equal(ind2,Ib)
  print,'counts:'
  print,count5,count1+count2,count3+count4
  
  stop
  
  ; OLD
  ; concatenate subindex arrays
  second_b_offset = count1
  
  inda       = [inda1,inda2]
  indb       = [indb1,indb2 + second_b_offset]
  
  ; now here we could delete full A and load full B, sort full B
  b_sorted    = b[sort(b)]
  b_subsorted = [b1(sort(b1)),b2(sort(b2))]
  
  b_subsortinds = [sort(b1),sort(b2) + second_b_offset]
  b_sortinds = sort(b)
  
  ; allocate arrays for our final solution
  inda_fixed = lonarr(num)
  indb_fixed = lonarr(num)
  
  ; walk through b_subsorted (and so inda,indb)
  if 0 then begin
  for i=0,count-1 do begin
    ; where in full b_sorted is this b_subsorted element?
    new_index = ( where(b_sorted eq b_subsorted[i]) )[0]
    new_index2 = b_subsorted[i] ; test, wrong
    
    ; move this inda index to the index corresponding to b_sorted (instead of b_subsorted)
    inda_fixed[new_index] = inda[i]
    indb_fixed[new_index] = indb[i]
    print,new_index,new_index2
  endfor
  endif ;0
  
end

; testMatching():

pro testMatching
  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function simParams, loadSnapshotSubset
  
  ; run each sort a few times and store mean runtime  
  nTests    = 4
  meanTimes = fltarr(nTests)
  
  nIter    = 1
  runTimes = fltarr(nIter)
  
  ; start memory monitor and timer
  print,'Loading...'
  
  mem = memory(/current) / 1024.0^3 ; GB
  time = systime(/sec) ; start timer (microsecond resolution)
  
  ; load some 64 bit data
  partType = 'gas'
  field    = 'ids'
  sP = simParams(res=128,run='tracer',redshift=2.0)
  ids_1 = loadSnapshotSubset(sP=sP,partType=partType,field=field)
  
  sP.snap -= 1
  ids_2 = loadSnapshotSubset(sP=sP,partType=partType,field=field)
  
  NumData = long(n_elements(ids_1))  
 
  ; report time and memory
  maxMem = memory(/highwater) / 1024.0^3 - mem
  curMem = memory(/current) / 1024.0^3
  print,'Load took: ['+string(systime(/sec)-time,format='(f5.1)')+'] sec, ['+$
    string(maxMem,format='(f5.2)')+' GB] max, ['+string(curMem,format='(f5.2)'),' GB] cur'
 
  ; run normal IDL match
  print,'Running IDL match...' 
  
  for i=0,nIter-1 do begin  
    time = systime(/sec)
  
    match,ids_1,ids_2,ind1_a,ind2_a,count=count_a
  
    maxMem = memory(/highwater) / 1024.0^3 - mem
    curMem = memory(/current) / 1024.0^3
    runTimes[i] = systime(/sec) - time
    print,'['+str(i)+'] IDL match took: ['+string(runTimes[i],format='(f5.1)')+'] sec, ['+$
      string(maxMem,format='(f5.2)')+' GB] max, ['+string(curMem,format='(f5.2)'),' GB] cur'
  endfor
  
  meanTimes[0] = mean(runTimes)

  ; run pure IDL mymatch
  print,'Running pure IDL mymatch...'
  
  for i=0,nIter-1 do begin  
    time = systime(/sec)
  
    mymatch,ids_1,ids_2,ind1_b,ind2_b,count=count_b
  
    maxMem = memory(/highwater) / 1024.0^3 - mem
    curMem = memory(/current) / 1024.0^3
    runTimes[i] = systime(/sec) - time
    print,'['+str(i)+'] IDL mymatch took: ['+string(runTimes[i],format='(f5.1)')+'] sec, ['+$
      string(maxMem,format='(f5.2)')+' GB] max, ['+string(curMem,format='(f5.2)'),' GB] cur'
  endfor
  
  meanTimes[1] = mean(runTimes)
  
  ; run CalcMatch one-pass
  print,'Running CalcMatch one-pass...'
  
  for i=0,nIter-1 do begin  
    time = systime(/sec)
    ids_2 = loadSnapshotSubset(sP=sP,partType=partType,field=field)
    CalcMatch,ids_1,ids_2,ind1_c,ind2_c,count=count_c
  
    maxMem = memory(/highwater) / 1024.0^3 - mem
    curMem = memory(/current) / 1024.0^3
    extMem = (n_elements(ids_1) + n_elements(ids_2)) * 4 / 1024.0^3
    runTimes[i] = systime(/sec) - time
    print,'['+str(i)+'] CalcMatch took: ['+string(runTimes[i],format='(f5.1)')+'] sec, ['+$
      string(maxMem,format='(f5.2)')+' GB] max, ['+string(extMem,format='(f5.2)')+' GB] ext, ['+$
      string(curMem,format='(f5.2)'),' GB] cur'
  endfor
  
  meanTimes[2] = mean(runTimes)
  
  ; run CalcMatch blocked
  print,'Running CalcMatch blocked...'
  
  for i=0,nIter-1 do begin  
    time = systime(/sec)
  
    CalcMatchBlock,ids_1,ind1_d,ind2_d,sP=sP,partType=partType,field=field,count=count_d
  
    maxMem = memory(/highwater) / 1024.0^3 - mem
    curMem = memory(/current) / 1024.0^3
    extMem = (n_elements(ids_1) + n_elements(ids_2)/10) * 4 / 1024.0^3
    runTimes[i] = systime(/sec) - time
    print,'['+str(i)+'] CalcMatchBlock took: ['+string(runTimes[i],format='(f5.1)')+'] sec, ['+$
      string(maxMem,format='(f5.2)')+' GB] max, ['+string(extMem,format='(f5.2)')+' GB] ext, ['+$
      string(curMem,format='(f5.2)'),' GB] cur'
  endfor
  
  meanTimes[3] = mean(runTimes)
  
  ; verify
  strings = ['IDL match  ','IDL mymatch','CalcMatch  ','CalcBlock  ']
  foreach sName,strings,i do print,sName+' : ['+string(meanTimes[i],format='(f5.1)')+'] sec.'
  
  if count_a ne count_b then message,'Error.'
  if count_a ne count_c then message,'Error.'
  if count_a ne count_d then message,'Error.'
  if ~array_equal(ind1_a,ind1_b) then message,'Error 1b.'
  if ~array_equal(ind2_a,ind2_b) then message,'Error 2b.'
  if ~array_equal(ind1_a,ind1_c) then message,'Error 1c.'
  if ~array_equal(ind2_a,ind2_c) then message,'Error 2c.'
  if ~array_equal(ind1_a,ind1_d) then message,'Error 1d.'
  if ~array_equal(ind2_a,ind2_d) then message,'Error 2d.'
  stop
  
end

pro monteCarloTestMatching

  nIter = 10
  
  ; load full data set
  partType = 'tracermc'
  field    = 'tracerids'
  sP = simParams(res=128,run='tracer',redshift=2.0)
  ids_1_full = loadSnapshotSubset(sP=sP,partType=partType,field=field)
  
  sP.snap -= 1
  ids_2_full = loadSnapshotSubset(sP=sP,partType=partType,field=field)
  
  ; loop
  for i=0,nIter-1 do begin
    ; make subsets, shuffle
    startInd = round( randomu(seed,1) * 0.5 * n_elements(ids_1) ) > 0
    endInd   = round( randomu(seed,1) * 2.0 * n_elements(ids_1) ) $
               > startInd + 1000 < n_elements(ids_1)-1   
    
    ids_1 = shuffle(ids_1_full[startInd:endInd])
  
    startInd = round( randomu(seed,1) * 0.5 * n_elements(ids_1) ) > 0
    endInd   = round( randomu(seed,1) * 2.0 * n_elements(ids_1) ) $
               > startInd + 1000 < n_elements(ids_1)-1  

    ids_2 = shuffle(ids_2_full[startInd:endInd])       
  
    ; run CalcMatch, CalcMatchBlock
    match,ids_1,ids_2,ind1_a,ind2_a,count=count_a
    
    CalcMatch,ids_1,ids_2,ind1_c,ind2_c,count=count_c
    
    ; verify
    if count_a ne count_c then message,'Error.'
    if ~array_equal(ind1_a,ind1_c) then message,'Error 1c.'
    if ~array_equal(ind2_a,ind2_c) then message,'Error 2c.'
    
    ; run CalcMatch with ids_2_full, CalcMatchBlock
    CalcMatch,ids_1,ids_2_full,ind1_a,ind2_a,count=count_a
    
    CalcMatchBlock,ids_1,ind1_c,ind2_c,sP=sP,partType=partType,field=field,count=count_c
    
    ; verify
    if count_a ne count_c then message,'Error.'
    if ~array_equal(ind1_a,ind1_c) then message,'Error 1c.'
    if ~array_equal(ind2_a,ind2_c) then message,'Error 2c.'
    
    print,i,startInd,endInd,'Passed.'
  endfor

end

; testSorting(): testing sort algorithms, IDL vs various C implementations, single vs multi threaded

pro testSorting
  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function simParams, loadSnapshotSubset
    
  ; run each sort a few times and store mean runtime
  nTests    = 4
  meanTimes = fltarr(nTests)
  
  nIter    = 3
  runTimes = fltarr(nIter)

  ; disable IDL threadpool
  ;CPU, tpool_nthreads = 1
  
  ; start memory monitor and timer
  print,'Loading...'
  
  mem = memory(/current) / 1024.0^3 ; GB
  time = systime(/sec) ; start timer (microsecond resolution)
  
  ; load some 64 bit data
  sP = simParams(res=256,run='tracer',redshift=2.0)
  ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
  
  NumData = long(n_elements(ids))  
  
  ; report time and memory
  maxMem = memory(/highwater) / 1024.0^3 - mem
  curMem = memory(/current) / 1024.0^3
  print,'Load took: ['+string(systime(/sec)-time,format='(f5.1)')+'] sec, ['+$
    string(maxMem,format='(f5.2)')+' GB] max, ['+string(curMem,format='(f5.2)'),' GB] cur'
  
  ; run normal IDL sort
  print,'Running IDL sort...'
  
  inds_0 = lonarr(NumData)
  
  for i=0,nIter-1 do begin
    time = systime(/sec) ; start timer
    inds_0 = sort(ids)
  
    maxMem = memory(/highwater) / 1024.0^3 - mem
    curMem = memory(/current) / 1024.0^3
    runTimes[i] = systime(/sec) - time
    print,'['+str(i)+'] IDL sort took: ['+string(runTimes[i],format='(f5.1)')+'] sec, ['+$
      string(maxMem,format='(f5.2)')+' GB] max, ['+string(curMem,format='(f5.2)'),' GB] cur'
  endfor
  
  meanTimes[0] = mean(runTimes)
  
  ; run external qsort (glibc, single threaded)
  print,'Running glibc qsort (single threaded)...'
  method = 1L

  inds_1 = lindgen(NumData)
  
  for i=0,nIter-1 do begin
    time = systime(/sec) ; start timer
  
    ret = Call_External('/n/home07/dnelson/idl/CalcMatch/CalcMatch_int64.so','CalcSort',$
                        NumData,ids,inds_1,method,/CDECL)
  
    maxMem = memory(/highwater) / 1024.0^3 - mem
    extMem = 0 ;(bytesMyInt * NumData + bytesMyInd * NumData) / 1024.0^3
    curMem = memory(/current) / 1024.0^3
    runTimes[i] = systime(/sec) - time
    print,'['+str(i)+'] glibc qsort took: ['+string(runTimes[i],format='(f5.1)')+'] sec, ['+$
      string(maxMem,format='(f5.2)')+' GB] max, ['+string(extMem,format='(f5.2)')+' GB] ext, ['+$
      string(curMem,format='(f5.2)'),' GB] cur'
  endfor
  
  meanTimes[1] = mean(runTimes)
  
  ; run external qsort (threaded, fhtr generic)
  print,'Running fhtr qsort (multi-threaded)...'
  method = 2L
  
  inds_2 = lindgen(NumData)
  
  for i=0,nIter-1 do begin
    time = systime(/sec) ; start timer
  
    ret = Call_External('/n/home07/dnelson/idl/CalcMatch/CalcMatch_int64.so','CalcSort',$
                        NumData,ids,inds_2,method,/CDECL)
  
    maxMem = memory(/highwater) / 1024.0^3 - mem
    extMem = 0 ;(bytesMyInt * NumData + bytesMyInd * NumData) / 1024.0^3
    curMem = memory(/current) / 1024.0^3
    runTimes[i] = systime(/sec) - time
    print,'['+str(i)+'] fhtr-generic qsort took: ['+string(runTimes[i],format='(f5.1)')+'] sec, ['+$
      string(maxMem,format='(f5.2)')+' GB] max, ['+string(extMem,format='(f5.2)')+' GB] ext, ['+$
      string(curMem,format='(f5.2)'),' GB] cur'
  endfor
  
  meanTimes[2] = mean(runTimes)
  
  ; run inplace
  print,'Running fhtr qsort inplace (multi-threaded)...'
  method = 12L
  inds_3 = [0]
  
  for i=0,nIter-1 do begin
    dataIn = ids
    time = systime(/sec) ; start timer
  
    ret = Call_External('/n/home07/dnelson/idl/CalcMatch/CalcMatch_int64.so','CalcSort',$
                        NumData,dataIn,inds_3,method,/CDECL)
  
    maxMem = memory(/highwater) / 1024.0^3 - mem
    extMem = 0 ; O(logn) at most since inplace
    curMem = memory(/current) / 1024.0^3
    runTimes[i] = systime(/sec) - time
    print,'['+str(i)+'] fhtr-inplace qsort took: ['+string(runTimes[i],format='(f5.1)')+'] sec, ['+$
      string(maxMem,format='(f5.2)')+' GB] max, ['+string(extMem,format='(f5.2)')+' GB] ext, ['+$
      string(curMem,format='(f5.2)'),' GB] cur'
  endfor
  
  meanTimes[3] = mean(runTimes)
  
  ; verify sorts agree
  strings = ['IDL   ','glibc ','fhtr1 ','fhtrIP']
  foreach sName,strings,i do print,sName+' sort: ['+string(meanTimes[i],format='(f5.1)')+'] sec.'
  
  if ~array_equal(inds_0,inds_1) then message,'Error: Fail 1.'
  if ~array_equal(inds_0,inds_2) then message,'Error: Fail 2.'
  stop
end

; sphDensityProjection(): (OLD) make density projection using SPH kernel (inspired by Mark's sphMap)
;                         NOTE: kernel coeffs only valid for 3D!

function sphDensityProjection, pos, hsml, mass, quantity=quantity, imgSize=imgSize, boxSize=boxSize,$
                               boxCen=boxCen, axis0=axis0, axis1=axis1, mode=mode, periodic=periodic,$
                               verbose=verbose

  print,'You should switch this to the calcSphMap C-routine.'
  stop

  ; config
  if not keyword_set(axis0) then axis0 = 0
  if not keyword_set(axis1) then axis1 = 1
  if not keyword_set(verbose) then verbose = 0
  
  if keyword_set(periodic) then begin
    print,'ERROR: PERIODIC not supported.'
    return,0
  endif
  
  if (mode ne 1 and mode ne 2 and mode ne 3) then begin
    print,'ERROR: Unsupported mode='+str(mode)+' parameter.'
    return,0
  endif
  
  ; storage
  p    = dblarr(3)
  pos0 = double(0.0)
  pos1 = double(0.0)
  binnedParticles = 0UL
  
  ; init
  npart = n_elements(hsml)

  grid = fltarr(imgSize[0],imgSize[1])
  
  if keyword_set(quantity) then $
    gridQuantity = fltarr(imgSize[0],imgSize[1])
  
  pxSize = [float(boxSize[0]) / imgSize[0], float(boxSize[1]) / imgSize[1]]
  pxArea = pxSize[0] * pxSize[1]

  if (pxSize[0] lt pxSize[1]) then $
    hMin = 1.001 * pxSize[0] / 2.0
  if (pxSize[0] ge pxSize[1]) then $
    hMin = 1.001 * pxSize[1] / 2.0
    
  hMax = pxSize[0] * 50.0
  
  for part=0, npart-1, 1 do begin
    ; progress report
    if (part mod round(npart/10.0) eq 0 and verbose) then $
      print,'Progress: '+string(100.0*part/npart,format='(I3)')+'%'
      
    ; get particle data
    p[0] = pos[0,part]
    p[1] = pos[1,part]
    p[2] = pos[2,part]
    h    = double(hsml[part])
    v    = double(mass[part])
    
    if keyword_set(quantity) then $
      w    = double(quantity[part])
    
    ; early exit if out of z-bounds
    if (abs(p[3-axis0-axis1] - boxCen[2]) gt boxSize[2] / 2.0) then $
      continue
      
    pos0 = p[axis0] - (boxCen[0] - boxSize[0] / 2.0)
    pos1 = p[axis1] - (boxCen[1] - boxSize[1] / 2.0)
    
    ; clamp hsml
    if (h lt hMin) then h = hMin;
    if (h gt hMax) then h = hMax;
    
    ; early exit if ...
    if (pos0 - 0.0 lt -h or pos1 - 0.0 lt -h or pos0 - boxSize[0] gt h or pos1 - boxSize[1] gt h) then $
      continue
      
    binnedParticles += 1
    
    h2 = h * h;
    
    ; number of pixels covered by particle
    nx = h / pxSize[0] + 1;
    ny = h / pxSize[1] + 1;
    
    ; coordinates of pixel center of particle
    x = (floor(pos0 / pxSize[0]) + 0.5) * pxSize[0]
    y = (floor(pos1 / pxSize[1]) + 0.5) * pxSize[1]
    
    ; normalization constant
    sum = 0.0
    
    for dx = -nx, nx, 1 do begin
      for dy = -ny, ny, 1 do begin
        ; dist of covered pixel from actual position
        xx = x + dx * pxSize[0] - pos0
        yy = y + dy * pxSize[1] - pos1
        r2 = xx*xx + yy*yy
        
        if (r2 < h2) then begin
          ; sph kernel (inlined): sum += _getkernel(h,r2);
          hinv = double(1.0) / h
          u    = sqrt(r2) * hinv
          
          if (u lt 0.5) then begin
            sum += (2.546479089470 + 15.278874536822 * (u - 1.0) * u * u)
          endif else begin
            sum += (5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u))
          endelse
        endif ;r2 < h2
      endfor
    endfor
    
    ; exit if negligible
    if (sum lt 1.0e-10) then $
      continue
      
    ; add contribution to image
    for dx = -nx, nx, 1 do begin
      for dy = -ny, ny, 1 do begin
        ; coordinates of pixel center of covering pixels
        xxx = x + dx * pxSize[0]
        yyy = y + dy * pxSize[1]
        
        ; pixel array indices
        i = floor(xxx / pxSize[0]) ;implicit C cast to int
        j = floor(yyy / pxSize[1]) ;same
        
        if (i ge 0 and i lt imgSize[0] and j ge 0 and j lt imgSize[1]) then begin
          xx = x + dx * pxSize[0] - pos0
          yy = y + dy * pxSize[1] - pos1
          r2 = xx*xx + yy*yy
          
          if (r2 lt h2) then begin
            ; divide by sum for normalization
            ; divide by pixelarea to get column density (optional: /pxArea)
            ; sph kernel (inlined): grid[] += _getkernel(h,r2) * v / sum
            hinv = double(1.0) / h
            u    = sqrt(r2) * hinv
            
            if (u lt 0.5) then begin
              grid[i * imgSize[1] + j] += $
                (2.546479089470 + 15.278874536822 * (u - 1.0) * u * u) * v / sum
              if keyword_set(quantity) then $
                gridQuantity[i * imgSize[1] + j] += $
                  (2.546479089470 + 15.278874536822 * (u - 1.0) * u * u) * v * w / sum
            endif else begin
              grid[i * imgSize[1] + j] += $
                (5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u)) * v / sum
                  if keyword_set(quantity) then $
                  gridQuantity[i * imgSize[1] + j] += $
                  (5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u)) * v * w / sum
            endelse
          
          endif ;r2 < h2
        endif ;i,j
      
      endfor
    endfor

  endfor ;part
  
  if (verbose) then print,'Number of binned particles: ',binnedParticles
  
  if (mode eq 1) then begin
    if (verbose) then print,'Returning: Column Mass Map'
    return,grid
  endif
  if (mode eq 2) then begin
    if (verbose) then print,'Returning: Quantity Mass-Weighted Map'
    return,gridQuantity
  endif
  if (mode eq 3) then begin
    if (verbose) then print,'Returning: Column Density Map'
    for i=0,i lt imgSize[0] do begin
      for j=0,j lt imgSize[1] do begin
        grid[i + imgSize[1] * j] /= pxArea
      endfor
    endfor
    
    return,grid
  endif

end