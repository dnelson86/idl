; helper.pro
; helper functions (NOTE: all my IDL routines loaded at bottom of this file)
; dnelson feb.2013

; one line utility functions
; --------------------------

function str, tString
  compile_opt idl2, hidden, strictarr, strictarrsubs
  return, strcompress(string(tString),/remove_all)
end

function isnumeric, input
  compile_opt idl2, hidden, strictarr, strictarrsubs
  on_ioerror, false
  test = double(input)
  return, 1
  false: return, 0
end

function getColor, i, name=name
  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function getUnits
  units = getUnits()
  ind = (i) mod (n_elements(units.colors)-1)
  
  if keyword_set(name) then return,units.colors[ind]
  return,cgColor(units.colors[ind])
end

function getColor24, color
    s = size(color)
    if s[0] eq 1 then begin
       if s[1] ne 3 then message, 'Error: Color must be a size 3 vector.'
       return, color[0] + (color[1] * 2L^8) + (color[2] * 2L^16)
    endif else begin
       if s[2] gt 3 then message, 'Error: Color must be an Nx3 array.'
       return, color[*,0] + (color[*,1] * 2L^8) + (color[*,2] * 2L^16)
    endelse
end

function linspace, a, b, N
  compile_opt idl2, hidden, strictarr, strictarrsubs
  vals = findgen(N) / (N-1.0) * (b-a) + a
  return, vals
end

function logspace, a, b, N, mid=mid
  compile_opt idl2, hidden, strictarr, strictarrsubs
  vals = findgen(N) / (N-1.0) * (b-a) + a
  
  ; return mid-bin points instead
  if keyword_set(mid) then $
    vals = (findgen(N-1)+0.5) / (N-1.0) * (b-a) + a
  
  vals = 10.0^vals
  
  return, vals
end

function nuniq, arr
  compile_opt idl2, hidden, strictarr, strictarrsubs
  return, n_elements(uniq(arr,sort(arr)))
end

function uniqvals, arr
  compile_opt idl2, hidden, strictarr, strictarrsubs
  return, arr[uniq(arr,sort(arr))]
end

function shuffle, array, seed=seed
  compile_opt idl2, hidden, strictarr, strictarrsubs
  if n_elements(seed) ne 0 then iseed=seed
  return,array[sort(randomu(iseed,n_elements(array)))]
end

; general algorithms
; ------------------

; replicate_var(): given a number of children for each parent, replicate the parent indices for each child

function replicate_var, child_counts

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  parent_inds = ulonarr(total(child_counts,/int))
  
  ; loop approach
  
  offset = 0UL
  for i=0UL,n_elements(child_counts)-1 do begin
    if child_counts[i] gt 0 then $
      parent_inds[ offset : offset+child_counts[i]-1 ] = replicate(i,child_counts[i])
    offset += child_counts[i]
  endfor
  
  ; vector approach:
  
  ;if n_elements(child_counts) le 1 then message,'error'
  ;cc = total(child_counts,/int,/cum)
  ;; remove last element if duplicate
  ;if cc[n_elements(cc)-1] eq cc[n_elements(cc)-2] then cc = cc[0:n_elements(cc)-2]
  ;      
  ;; mark first child of each parent
  ;parent_inds[cc[0:n_elements(cc)-2]] = 1
  ;
  ;; compensate for zero child counts
  ;w = where(child_counts eq 0,count)
  ;if count gt 0 then begin
  ; 
  ;  ; need to handle multiple zeros on end
  ;  ; TODO
  ;  if w[-1] eq n_elements(child_counts)-1 then w = w[0:-2]
 ; 
 ;   ; possible multiple zeros sequentially
 ;   nummod = histogram(cc[w],loc=locmod)
 ;   parent_inds[locmod] += nummod
 ;   
 ; endif
 ; 
 ; ; scan sum
 ; parent_inds = total(parent_inds,/cum,/int)
 
  return,parent_inds
end

; pSplit(): divide work for embarassingly parallel problems

function pSplit, arr, split=split ; split=[numProcs,curProc] with curProc zero indexed

  ; no split, return whole job load to caller
  if n_elements(split) ne 2 then return,arr

  ; split arr into split[0] segments and return split[1] segment
  splitSize = fix(n_elements(arr) / split[0])
  arrSplit  = arr[split[1]*splitSize:(split[1]+1)*splitSize-1]
  
  ; for last split, make sure it takes any leftovers
  if split[0]-1 eq split[1] then arrSplit = arr[split[1]*splitSize:*]
  
  return,arrSplit
end

; fitODRPts3D(): Orthogonal Distance Regression method to fit a line/plane through points in 3D

function fitODRPts3D, x, y, z

  if n_elements(x) ne n_elements(y) or n_elements(x) ne n_elements(z) then message,'Error'
  
  data = transpose([[x],[y],[z]])
  centroid = total(data,2) / n_elements(x)
  
  data[0,*] -= centroid[0]
  data[1,*] -= centroid[1]
  data[2,*] -= centroid[2]
  
  SVDC, data, W, U, V
  
  smallest_singVal = min(W,ind)
  plane_normal = reform(V[ind,*])
  
  largest_singVal = max(W,ind)
  line_unitVec = reform(V[ind,*])
  
  r = { centroid:centroid, line_unitvec:line_unitvec, plane_normal:plane_normal }
  return,r
end

; gcDist(): "great circle" (i.e. arc) distance between two lat,long points on the surface of a sphere

function gcDist, latlong1, latlong2 ; in radians
  d = acos( sin(latlong1[0]) * sin(latlong2[0]) + $
            cos(latlong1[0]) * cos(latlong2[0]) * cos(latlong1[1]-latlong2[1]) )
  ; A mathematically equivalent formula, which is less subject to rounding error for short distances is:
  ; d=2*asin(sqrt((sin((lat1-lat2)/2))^2 + cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2))^2))
  return,d
end

; flatten_list

function flatten_list, list
  compile_opt idl2, hidden, strictarr, strictarrsubs
  if (list.Count() eq 0) then return,[0]
  
  arr = []
  for i=0ULL,list.Count()-1 do $
    arr = [arr,list[i]]

  return,arr
  
end

; removeIntersectionFromB(): return a modified version of B with all those elements also found in
;                            A (the collision/intersection) removed

function removeIntersectionFromB, A, B, union=union

  compile_opt idl2, hidden, strictarr, strictarrsubs

    match, A, B, A_ind, B_ind, count=count, /sort
    
    A_ind = !NULL ;unused
    
    if (count gt 0) then begin
      ; remove B[B_ind] using complement
      all = bytarr(n_elements(B))
      if (B_ind[0] ne -1L) then all[B_ind] = 1B
      w = where(all eq 0B, ncomp)
    
      if (ncomp ne n_elements(B)-count) then $
        message,'removeIntersectionFromB: ERROR!'
      
      ; set union return
      if (keyword_set(union)) then union=B[B_ind]
      
      return, B[w]
    endif else begin
      print,'Warning: removeIntersectionFromB returning unmodified.'
      return, B
    endelse
end

; getIDIndexMapSparse(): return an array which maps ID->indices within dense, disjoint subsets which are 
;                        allowed to be sparse in the entire ID range. within each subset i of size binsize
;                        array[ID-minids[i]+offset[i]] is the index of the original array ids where ID 
;                        is found (assumes no duplicate IDs)

function getIDIndexMapSparse, ids, map=map

  compile_opt idl2, hidden, strictarr, strictarrsubs

  minid = min(ids)
  maxid = max(ids)

  if n_elements(ids) gt 2e9 then message,'Error: Going to overrun arr.'

  ; if map passed in, use it to return the indices of the ids instead of making a new map
  if keyword_set(map) then begin
    inds = lonarr(n_elements(ids))
    totCount = 0ULL
    
    for i=0,map.nBlocks-1 do begin
      ; process ids in this block
      w = where(ids ge map.blockStarts[i] and ids le map.blockEnds[i],count)
      
      if count gt 0 then inds[w] = map.map[ids[w]-map.blockStarts[i]+map.blockOffsets[i]]
      totCount += count
    endfor
    
    if totCount ne n_elements(ids) then message,'Error: Failed to map all IDs to indices.'
    return,inds
  endif

  ; make broad histogram of ids to identify sparsity
  binsize = round(sqrt(max(ids))/10000.0) * 10000 > 10 < 1000000 ; adjust to range
  print,binsize
  hist = histogram(ids,binsize=binsize,min=0,loc=loc)
  occBlocks = where(hist ne 0,count)
  
  if count eq 0 then message,'Error.'
  
  nBlocks = 0
  blockStarts = []
  blockLens   = []
  blockEnds   = []
  prevBlock = ulong64(0)
  count = 0
  
  ; search for contiguous id blocks
  for i=0,n_elements(occBlocks)-1 do begin
    count += 1
    
    ; jump detected (zero block length ge binsize)
    if occBlocks[i] ne (prevBlock+1) then begin
      nBlocks += 1
      blockStarts = [blockStarts,loc[occBlocks[i]]]
      blockLens   = [blockLens,count]
      blockEnds   = [blockEnds,loc[prevBlock]+count*binsize-1]
      if i ne n_elements(occBlocks)-1 then count = 0
    endif
    
    prevBlock = occBlocks[i]
  endfor
  
  ; remove first (zero) block length and add final
  if nBlocks gt 1 then begin
    count += 1
    blockLens = [blockLens[1:*],count]
    blockEnds = [blockEnds[1:*],loc[prevBlock]+count*binsize-1]
  endif
  
  if nBlocks le 0 then message,'Error: No dense blocks found.'
  
  ; make offset array into mapping
  blockOffsets = ulonarr(nBlocks)
  for i=0,nBlocks-2 do begin
    blockOffsets[i+1] = blockOffsets[i] + blockLens[i]*binsize
  endfor
  
  ; create mapping array inside return structure
  r = { map          : ulonarr(total(blockLens)*binsize) ,$
        blockStarts  : blockStarts                       ,$
        blockLens    : blockLens                         ,$
        blockEnds    : blockEnds                         ,$
        blockOffsets : blockOffsets                      ,$
        nBlocks      : nBlocks                            }
  
  ; fill mapping for each subset
  totCount = 0ULL
  
  for i=0,nBlocks-1 do begin
    ; what ids in this block?
    w = where(ids ge blockStarts[i] and ids le blockEnds[i],count)
    if count eq 0 then message,'Error'
    
    for j=0ULL,n_elements(w)-1 do r.map[ids[w[j]]-blockStarts[i]+blockOffsets[i]] = w[j] ;j + blockOffsets[i]
    totCount += count
  endfor
  
  if totCount ne n_elements(ids) then message,'Error: Did not map all IDs.'
  if nuniq(r.map) ne n_elements(ids) then message,'Error: Non-unique mapping.' ; DEBUG ONLY (slow)
  
  ; DEBUG: for a subset, verify mapping
  nVerify = 100 < n_elements(ids)
  indVerify = floor(randomu(seed,nVerify)*n_elements(ids))
  indTest = getIDIndexMap(ids[indVerify],map=r)
  
  for i=0,nVerify-1 do begin
    w = where(ids eq ids[indVerify[i]],count)
    ;print,i,w[0],indTest[i]
    if count ne 1 or w[0] ne indTest[i] then message,'Error. Bad mapping.'
  endfor
  
  return,r
  
end
  
; getIDIndexMap(): return an array of size max(ids)-min(ids) such that array[ID-min(ids)] is the 
;                     index of the original array ids where ID is found (assumes a one to one mapping, 
;                     not repeated indices as in the case of parentIDs for tracers)

function getIDIndexMap, ids, minid=minid

  compile_opt idl2, hidden, strictarr, strictarrsubs

  minid = min(ids)
  maxid = max(ids)

  if n_elements(ids) gt 2e9 then message,'Error: Going to overrun arr.'
  
  ; C-style loop approach (good for sparse IDs)
  arr = ulonarr(maxid-minid+1)
  for i=0ULL,n_elements(ids)-1L do arr[ids[i]-minid] = i
  
  ; looped where approach (never a good idea)
  ;arr = l64indgen(maxid-minid+1)
  ;for i=minid,maxid do begin
  ;  w = where(ids eq i,count)
  ;  if (count gt 0) then arr[i] = w[0]
  ;endfor

  ; reverse histogram approach (good for dense ID sampling, maybe better by factor of ~2)
  ;arr = l64indgen(maxid-minid+1)
  ;h = histogram(ids,rev=rev,omin=omin)
  ;for i=0L,n_elements(h)-1 do if (rev[i+1] gt rev[i]) then arr[i] = rev[rev[i]:rev[i+1]-1]

  return, arr
end

; basic IO
; --------

; loadColorTable(): load a custom or builtin IDL color table into the display
;                   if rgb_table specified do not load into display, just return the table

pro loadColorTable, ctName, bottom=bottom, rgb_table=rgb_table, reverse=reverse

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; cubehelix CT implementation
  if ctName eq 'helix' then begin
    ; helix parameters
    start = 2.0  ; color, 1=r,2=g,3=b,0.5=purple
    rots  = 1.5  ; color rotations, typically -1.5 to 1.5
    hue   = 1.0  ; hue intensity scaling, typically 0 (BW) to 1
    gamma = 0.75 ; gamma intensity expontent (1=normal/linear,<1 emphasize low vals,>1 emphasize high vals)
    
    ; calculate R,G,B based on helix parameters
    nlev = 256
    RGB = fltarr(nlev,3)
    
    fracs = findgen(nlev)/float(nlev-1) * 1.0; scale [0,1]
    phi   = 2*!pi*(start/3.0 + 1.0 + rots*fracs)
    fracs = fracs^float(gamma)
    amplt = hue * fracs * (1-fracs)/2.0
    
    RGB[*,0] = fracs + amplt * (-0.14861*cos(phi) + 1.78277*sin(phi)) ;R
    RGB[*,1] = fracs + amplt * (-0.29227*cos(phi) - 0.90649*sin(phi)) ;G
    RGB[*,2] = fracs + amplt * (+1.97294*cos(phi)) ;B
    
    ; clip to [0,1] and expand to [0,255]
    w = where(RGB lt 0.0,count)
    if count gt 0 then RGB[w] = 0.0
    w = where(RGB gt 1.0,count)
    if count gt 0 then RGB[w] = 1.0
    
    RGB *= 255.0 ; cubehelix creates RGB in [0,1] but we use [0,255] for IDL display devices
    RGB = fix(round(RGB)) ; rounded INT
    
    ; reverse (light=low to dark=high) if requested
    if keyword_set(reverse) then begin
      RGB = reverse(RGB)
      RGB[0,*] = [0,0,0]
    endif
    
    ; fill rgb_table if requested, or if not, load RGB into display
    if arg_present(rgb_table) then rgb_table = RGB

    if not arg_present(rgb_table) then begin
      ; start above zero on the CT?
      cbot = 0
      if n_elements(bottom) gt 0 then cbot = bottom > 0 < nlev-1
      
      tvlct, RGB, cbot
    endif
    return
  endif
  
  ; otherwise, load normal IDL table
  if ctName eq 'bw linear'          then loadct,0,bottom=bottom,/silent
  if ctName eq 'green-white linear' then loadct,8,bottom=bottom,/silent
  if ctName eq 'green-white exp'    then loadct,9,bottom=bottom,/silent
  if ctName eq 'blue-red'           then loadct,11,bottom=bottom,/silent
  if ctName eq 'plasma'             then loadct,32,bottom=bottom,/silent
  if ctName eq 'blue-red2'          then loadct,33,bottom=bottom,/silent
  
  ; brewer (nice diverging options, both ends dark, light in middle) (17-26)
  if ctName eq 'brewer-brownpurple' then cgLoadct,18,ncolors=255,bottom=bottom,/brewer
  if ctName eq 'brewer-browngreen'  then cgLoadct,19,ncolors=255,bottom=bottom,/brewer
  if ctName eq 'brewer-purplegreen' then cgLoadct,20,ncolors=255,bottom=bottom,/brewer
  if ctName eq 'brewer-pinkgreen'   then cgLoadct,21,ncolors=255,bottom=bottom,/brewer
  if ctName eq 'brewer-redblue'     then cgLoadct,22,ncolors=255,bottom=bottom,/brewer
  if ctName eq 'brewer-redpurple'   then cgLoadct,25,ncolors=255,bottom=bottom,/brewer
  if ctName eq 'brewer-redgreen'    then cgLoadct,26,ncolors=255,bottom=bottom,/brewer
 
 ; brewer (convergent, light->dark)
 if ctName eq 'brewerc-yellowblue'      then cgLoadct,1,ncolors=255,bottom=bottom,/brewer
 if ctName eq 'brewerc-whitegreen'      then cgLoadct,3,ncolors=255,bottom=bottom,/brewer
 if ctName eq 'brewerc-purplebluegreen' then cgLoadct,4,ncolors=255,bottom=bottom,/brewer
 if ctName eq 'brewerc-redpurple'       then cgLoadct,7,ncolors=255,bottom=bottom,/brewer
 if ctName eq 'brewerc-orangered'       then cgLoadct,9,ncolors=255,bottom=bottom,/brewer
 if ctName eq 'brewerc-blues'           then cgLoadct,13,ncolors=255,bottom=bottom,/brewer
 if ctName eq 'brewerc-greens'          then cgLoadct,14,ncolors=255,bottom=bottom,/brewer
 if ctName eq 'brewerc-reds'            then cgLoadct,16,ncolors=255,bottom=bottom,/brewer
 
  ; reverse normal/brewer if requested
  if keyword_set(reverse) then begin
    tvlct,r,g,b,/get
    tvlct,reverse(r),reverse(g),reverse(b)
  endif
 
  ; brewer (sequential) (0-16)

end

; loadCSV(): load all the lines of a textfile with column template ptStruct
;            skip headerLines at the beginning and put their contents as a string into header

function loadCSV, headerLines, fileName, ptStruct, header=header;, format=format

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ;prepare data containers
  nRows = File_Lines(fileName)
  if (nRows - headerLines ne 0) then $
    pts  = replicate(ptStruct,nRows - headerLines)
  if (headerLines ne 0) then $
    header = strarr(headerLines)
  
  ;open and read file
  openR, lun, fileName, /GET_LUN
  
  if (headerLines ne 0) then $
    readF, lun, header
  if (nRows - headerLines ne 0) then $
    readF, lun, pts
  
  ;close handle
  free_lun, lun

  if (nRows - headerLines ne 0) then $
    return, pts
  
end

; loadBinary(): right now just reads Stars_X.bin (first float indicates how many pts)

function loadBinary, fileName, ptStruct

  compile_opt idl2, hidden, strictarr, strictarrsubs

  openR, lun, fileName, /GET_LUN
    nPts = 0UL
    readU, lun, nPts ; header
    pts  = replicate(ptStruct,nPts) ; replicate
    readU, lun, pts ; fill
  close, lun
  free_lun, lun

  return, pts
end

; loadBinarySequence(): right now just reads Stars_X_Y where X=num, Y=node

function loadBinarySequence, fileBase, ptStruct

  compile_opt idl2, hidden, strictarr, strictarrsubs

  fileNames = file_search(fileBase+'*')
  pCount = 0UL
  tempPt = replicate(ptStruct,1)
  
  ; open first file
  openR, lun, fileNames[0], /GET_LUN
    ;header
    nPts = 0UL
    readU, lun, nPts
    ;sanity check
    if (nPts eq 0) then begin
      print,'WARNING! nPts in Stars0 eq 0, hardcoding to 10mil.'
      nPts = 10000000UL
    endif
    ;replicate
    pts = replicate(ptStruct,nPts)
    ;fill from first file
    while( not EOF(lun) ) do begin
      readU, lun, tempPt
      pts[pCount] = tempPt
      pCount = pCount + 1
    endwhile
  close, lun
  free_lun, lun
  
  ; loop over remaining files
  for i=1,n_elements(fileNames)-1 do begin
    openR, lun, fileNames[i], /GET_LUN
    
    while( not EOF(lun) ) do begin
      readU, lun, tempPt
      pts[pCount] = tempPt
      pCount = pCount + 1
    endwhile
    
    close, lun
    free_lun, lun
  endfor
  
  return, pts
end

; postscript output
; -----------------

pro start_PS, filename, xs=xs, ys=ys, eps=eps, big=big, extrabig=extrabig

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if not keyword_set(xs) then xs=7.5
  if not keyword_set(ys) then ys=5.0
  if n_elements(eps) eq 0 then eps=1
  
  ; make the page bigger
  if n_elements(big) eq 1 then begin
    xs *= 1.2 ;9.0
    ys *= 1.2 ;6.0
  endif
  if n_elements(extrabig) eq 1 then begin
    xs *= 1.4 ;10.5
    ys *= 1.4 ;7
  endif

  PS_Start, FILENAME=filename, /nomatch, /quiet, bits_per_pixel=8, color=1, $
            encapsulated=eps, decomposed=0, xs=xs, ys=ys, /inches, font=0;, $
            ;/dejavusans ;requires idl8.2 
            ;font=1,tt_font='Helvetica'
 
  !p.charsize  = 1.4
  ;!p.charthick = 1.4
  !p.thick    = 5.0
  
  ; make axis lines thicker
  !x.thick += 1.0
  !y.thick += 1.0
  !z.thick += 1.0
  
  ; reduce number of minor ticks to clean up plots
  !x.minor = 2
  !y.minor = 2
  
  ; make default custom psym (8) a filled circle
  plotsym,0,/fill
  !p.symsize = 0.7
            
end

pro end_PS, pngResize=pngResize, deletePS=deletePS, im_options=im_options

  compile_opt idl2, hidden, strictarr, strictarrsubs

 ;PNG size=[xs,ys]*300*(resize/100)

  if not keyword_set(pngResize) then $
    PS_End
    
  if keyword_set(pngResize) then $
    PS_End, /PNG, Delete_PS=deletePS, im_options=im_options, Resize=pngResize, /nofix;, /showcmd

end

; save_eps(): idl 8.x compatible EPS save

pro save_eps, p, plotName, width=width, height=height, savePDF=savePDF, savePNG=savePNG

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if not keyword_set(width)  then width=7.5
  if not keyword_set(height) then height=5.0

  p->save, plotName, page_size=[width,height]
  
  ; other save formats
  if (keyword_set(savePDF)) then begin
    plotName = strmid(plotName,0,strlen(plotName)-4) + '.pdf'
    p->save, plotName, page_size=[width,height]
  endif

end

; multithreading
; --------------

; reportCPUThreads(): return info from CPU structure for threading pool

pro reportCPUThreads
  compile_opt idl2, hidden, strictarr, strictarrsubs
  help,/memory
  print,'!CPU.HW_NCPU = ' + str(!CPU.HW_NCPU) + ' TPool_NThreads = ' + str(!CPU.TPOOL_NTHREADS)
end

; testBridgePro

pro testBridgePro
  compile_opt idl2, hidden, strictarr, strictarrsubs
  ; retrieve from $MAIN$
  ;redshift = SCOPE_VARFETCH("redshift", LEVEL=1)
  ;res      = SCOPE_VARFETCH("res", LEVEL=1)

  workingPath  = '/n/home07/dnelson/coldflows/'
  
  result = redshift*2.0
  
  wait,2.0
  
  save,result,filename=workingPath+"test."+str(redshift)+".sav"

end

; runBridge():

pro runBridge, res=res
  compile_opt idl2, hidden, strictarr, strictarrsubs
  reportCPUThreads
  
  redshifts = [3.0,2.0,1.0,0.0]
 
  start_time = systime(/seconds)
  
  ; launch children
  oB = objarr(n_elements(redshifts))
  for i=0,n_elements(redshifts)-1 do begin
    oB[i] = obj_new('IDL_IDLBridge') ;OUTPUT='' send child output to tty
    oB[i]->execute, ".r coldflows"
    oB[i]->SetVar, "redshift", redshifts[i]
    oB[i]->SetVar, "res", res
    oB[i]->execute, "testBridgePro", /NOWAIT ; asynchronous
  endfor
  
  ; wait for children to finish and cleanup
  for i=0,n_elements(redshifts)-1 do $
    while (oB[i]->Status() ne 0) do wait,0.1
  obj_destroy,oB
  
  print,"Elapsed time: "+str(systime(/seconds)-start_time)+" sec"

end

; load routines for use
; ---------------------
@externalC
@units
@haloModels
@simParams
@cosmoUtil
@groupCat
@cosmoLoad

@mergerTree
@galaxyCat
@accretionMode
@accretionTimes
@accretionVel
@maxTemps
@timeScales
@timeScalesPlot

@accretionTraj
@accretionTrajVis
@cosmoVis
@cosmoVisMap
@cosmoSphere
@cosmoOverDens
@haloCompProj
@galaxyCatVis

@plotGalCat
@plotMaxTemps
@plotSphere
@binVsHaloMass
@plotVsHaloMass
@plotRadProfiles

@tracersVel_Cosmo
@tracersVel_Halos
;@tracersVel_Disks
;@tracersVel_2D
;@tracersVel_SphSym

@tracersMC
@tracersMC_Halos
;@tracersMC_2D
;@tracersMC_SphSym

@arepoLoad
@arepoVis2D
@arepoSphSym

@filamentTest
@LSF
