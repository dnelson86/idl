; helper.pro
; helper functions (NOTE: all my IDL routines loaded at bottom of this file)
; dnelson jan.2014

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

function mylog10, x ; leave zeros as zero, instead of converting to NaN
  w = where(x eq 0,count)
  if count gt 0 then x[w] = 1.0
  x = alog10(x)
  return, x
end

function closest, arr, x
  w = where( abs(arr-x) eq min(abs(arr-x)), count )
  if count eq 0 then message,'error'
  if count ge 1 then w = w[0]
  return,w
end

function percentiles, arr, percentiles=percentiles

  if n_elements(percentiles) eq 0 then percentiles = [0.25, 0.50, 0.75]
   
  if min(percentiles) lt 0.0 or max(percentiles) gt 1.0 then message,'Error'

  index = sort(arr)
  ind = findgen(n_elements(arr)+1)/float(n_elements(arr))
  arr_ind = value_locate(ind, percentiles)
  result = arr[index[arr_ind]]

  return, result
   
end

pro reportMemory, msg=msg
  if ~keyword_set(msg) then msg = ''
  mem = memory()
  high_GB = string( float(mem[3]) / (1024L^3), format='(f6.2)')
  cur_GB  = string( float(mem[0]) / (1024L^3), format='(f6.2)')
  print,msg+': memory highwater: ['+high_GB+'] GB, current: ['+cur_GB+'] GB. ('+str(mem[1])+' '+str(mem[2])+')'
end

; general algorithms
; ------------------

; replicate_var(): given a number of children for each parent, replicate the parent indices for each child
; subset_inds=1 : still need to walk the full child_counts, but only want the parent indices of a subset

function replicate_var, child_counts, subset_inds=subset_inds

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; loop approach
  offset = 0UL
  
  if ~keyword_set(subset_inds) then begin
    ; full
    parent_inds = ulonarr(total(child_counts,/int))
    
    for i=0UL,n_elements(child_counts)-1 do begin
      if child_counts[i] gt 0 then $
        parent_inds[ offset : offset+child_counts[i]-1 ] = replicate(i,child_counts[i])
      offset += child_counts[i]
    endfor
    
    return,parent_inds
    
  endif
  
  if keyword_set(subset_inds) then begin
    ;subset
    totChildren = total(child_counts[subset_inds],/int)
    
    ; we also return the child index array (i.e. which children belong to the subset_inds parents)
    r = { parent_inds : ulonarr(totChildren) ,$
          child_inds  : ulonarr(totChildren)  }
    
    offset_sub = 0UL
    subset_mask = bytarr(n_elements(child_counts))
    subset_mask[subset_inds] = 1B
    
    for i=0UL,n_elements(child_counts)-1 do begin
      if subset_mask[i] eq 1B and child_counts[i] gt 0 then begin
        r.parent_inds[ offset_sub : offset_sub+child_counts[i]-1 ] = replicate(i,child_counts[i])
        r.child_inds[ offset_sub : offset_sub+child_counts[i]-1 ] = lindgen(child_counts[i]) + offset
        offset_sub += child_counts[i]
      endif
      offset += child_counts[i]
    endfor
    
    return, r
    
  endif
  
  ; vector approach (full, uncomplete):
  
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
  ;  message,'TODO'
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

; intersection(): return intersection of A and B
; xor=1: only values in either A or B, but not both, are returned

function intersection, A, B, ret_xor=ret_xor

  array = [A, B]
  array = array(sort(array))

  if keyword_set(ret_xor) then begin
    samp1 = intarr(n_elements(array))
    samp2 = samp1
    
    i1 = where(array ne shift(array, -1),count)
    if count gt 0 then samp1[i1]=1
    i2 = where(array ne shift(array,  1),count)
    if count gt 0 then samp2[i2]=1
    
    indices = where(samp1 eq samp2, count)
  
  endif else begin
    indices = where(array eq shift(array, -1), count)
  endelse
  
  if count gt 0 then return, array[indices]
  return, -1

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

  if n_elements(ids) gt 2e9 then message,'Error: Going to overrun arr, switch to long64.'
  
  ; Paul's direct indexing approach (fastest)  
  arr = lonarr(maxid-minid+1)
  arr[ids-minid] = lindgen(n_elements(ids))
  
  ; C-style loop approach (good for sparse IDs)
  ;arr = ulonarr(maxid-minid+1)
  ;for i=0ULL,n_elements(ids)-1L do arr[ids[i]-minid] = i
  
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

; plotting related
; ----------------

; fitRadProfile(): helper function to fit median radial profile

function fitRadProfile, radii=radii, vals=vals, range=range, radBins=radBins

  r = { binSize    : (range[1]-range[0])/radBins          ,$
        binCen     : linspace(range[0],range[1],radBins)  ,$
        radMean    : fltarr(radBins)                      ,$
        radMedian  : fltarr(radBins)                      ,$
        radStddev  : fltarr(radBins)                      ,$
        radNum     : lonarr(radBins)                       }
        
  for i=0,radBins-1 do begin
    binEdges = range[0] + [i,i+1]*r.binSize
    w = where(radii ge binEdges[0] and radii lt binEdges[1],count)
    if count gt 0 then begin
      r.radMean[i]   = mean(vals[w])
      r.radMedian[i] = median(vals[w])
      r.radStddev[i] = stddev(vals[w])
      r.radNum[i]    = count
    endif
  endfor
         
  return,r

end

; binHisto2D()

function binHisto2D, xx=xx, yy=yy, wt=wt, xmm=xmm, ymm=ymm, xbs=xbs, ybs=ybs

  compile_opt idl2, hidden, strictarr, strictarrsubs
  if ~keyword_set(xx) or ~keyword_set(yy) or ~keyword_set(xmm) or ~keyword_set(ymm) or $
     ~keyword_set(xbs) or ~keyword_set(ybs) then message,'error'
     
  if ~keyword_set(wt) then wt = replicate(1.0,n_elements(xx))
  
  nXBins = ceil((xmm[1]-xmm[0])/xbs)
  nYBins = ceil((ymm[1]-ymm[0])/ybs)
  
  binCenX = linspace(xmm[0],xmm[1]-xbs,nXBins) + xbs/2
  binCenY = linspace(ymm[0],ymm[1]-ybs,nYBins) + ybs/2
  
  ; return array
  r = { nXBins:nXBins, nYBins:nYBins, binCenX:binCenX, binCenY:binCenY, $
        binSizeX:xbs, binSizeY:ybs, $
        h2:fltarr(nXBins,nYBins)  }
        
  ; bin manually for mass weighting
  for i=0,nXBins-1 do begin
    xBin = [binCenX[i]-xbs/2,binCenX[i]+xbs/2]
  
    for j=0,nYBins-1 do begin
      yBin = [binCenY[j]-ybs/2,binCenY[j]+ybs/2]
      
      w = where(xx ge xBin[0] and xx lt xBin[1] and yy ge yBin[0] and yy lt yBin[1],count)
      
      if count gt 0 then r.h2[i,j] = total(wt[w])
    
    endfor
  endfor

  return, r

end

; oplot2DHistoSq(): plot a 2d histogram in the separated squares style
; hsp: 2d array, gap size (in data units) along x and y directions
; nc: highest number of the color table to use (e.g. for 0-255 normally white-dark, exclude the darkest colors)
; nonZero=1: replace all zero counts by lowest nonzero value
; logX/Y=1: indicates the y-axis of the plot is log
; colNorm=1: normalize each column independently (i.e. for xaxis=halo mass to offset mass function)
    
pro oplot2DHistoSq, ct2d, hsp=hsp, nc=nc, xRange=xRange, yRange=yRange, $
  nonZero=nonZero, logY=logY, logX=logX, colNorm=colNorm, gray=gray, green=green, ctName=ctName

  if ~keyword_set(hsp) or ~keyword_set(nc) then message,'error'

  ; store current table and load for 2d histo
  tvlct, rr, gg, bb, /get
  
  if keyword_set(gray)   then loadColorTable,'bw linear', /reverse
  if keyword_set(green)  then loadColorTable,'green-white linear', /reverse
  if keyword_set(ctName) then loadColorTable, ctName
    
  ; process data (stretching data values)
  w = where(ct2d.h2 eq 0.0,count,comp=wc)
    
  if keyword_set(nonZero) then begin
    ct2d.h2 = ct2d.h2*1e10
    if count gt 0 then ct2d.h2[w] = min(ct2d.h2,/nan)
  endif
  
  ; normalize column by column
  if keyword_set(colNorm) then begin
    for i=0,ct2d.nXBins-1 do begin
      colTot = max(ct2d.h2[i,*],/nan)
      ct2d.h2[i,*] /= colTot
    endfor
  endif
    
  ; colorscale range
  fieldMinMax = [min(ct2d.h2[wc],/nan),max(ct2d.h2,/nan)]
    
  for i=0,ct2d.nXBins-1 do begin
    x = [ct2d.binCenX[i]-ct2d.binSizeX/2+hsp[0], ct2d.binCenX[i]+ct2d.binSizeX/2-hsp[0], $ ; ll, lr
         ct2d.binCenX[i]+ct2d.binSizeX/2-hsp[0], ct2d.binCenX[i]-ct2d.binSizeX/2+hsp[0]]   ; ur, ul
           
    if keyword_set(logX) then x = 10.0^x
           
    if min(x) lt xRange[0] or max(x) gt xRange[1] then continue
             
    for j=0,ct2d.nYBins-1 do begin
      y = [ct2d.binCenY[j]-ct2d.binSizeY/2+hsp[1], ct2d.binCenY[j]-ct2d.binSizeY/2+hsp[1], $ ; ll, lr
           ct2d.binCenY[j]+ct2d.binSizeY/2-hsp[1], ct2d.binCenY[j]+ct2d.binSizeY/2-hsp[1]]   ; ur, ul
             
      if keyword_set(logY) then y = 10.0^y
             
      if min(y) lt yRange[0] or max(y) gt yRange[1] then continue
             
      ; determine color and make polygon
      colorind = (ct2d.h2[i,j]-fieldMinMax[0])*nc / (fieldMinMax[1]-fieldMinMax[0]) ;0-nc
      colorind = fix(colorind + 0.0) > 0 < 255 ;0-nc
      
      if colorind lt ceil(0.01*nc) then continue ; skip marginal bins
        
      cgPolygon,x,y,color=colorind,/fill
    endfor
  endfor
  
  ; restore original CT
  tvlct, rr, gg, bb
  
end

; oplotBand(): fill band between two lines with color

pro oplotBand, x, y1, y2, color=color, yrange=yrange

  if n_elements(x) ne n_elements(y2) or n_elements(x) ne n_elements(y2) then $
    message,'Error'
    
  xp  = x
  y1p = y1 ;smooth(y1,3)
  y2p = y2 ;smooth(y2,3)
  
  xx=[xp,reverse(xp)]
  yy=[y1p,reverse(y2p)]

  ymin = min(yrange)
  ymax = max(yrange)
  yy=(ymin > yy) < ymax
  ;xx=(xmin > xx) < xmax
  
  cgColorFill,xx,yy,color=color ;polyfill,xxx,yyy,_extra=_extra

end

; basic IO
; --------

; loadColorTable(): load a custom or builtin IDL color table into the display
;                   if rgb_table specified do not load into display, just return the table

pro loadColorTable, ctName, bottom=bottom, rgb_table=rgb_table, reverse=reverse, gamma=GA_in, $
  GA_double2=GA_double2

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; sabotage CT at beginning to make sure we load something
  tvlct,fltarr(256,3)
  
  ; cubehelix CT implementation
  if ctName eq 'helix' then begin
    ; helix parameters
    ;start = 2.0  ; color, 1=r,2=g,3=b,0.5=purple
    ;rots  = 1.5  ; color rotations, typically -1.5 to 1.5
    ;hue   = 1.0  ; hue intensity scaling, typically 0 (BW) to 1
    ;gamma = 0.75 ; gamma intensity expontent (1=normal/linear,<1 emphasize low vals,>1 emphasize high vals)
    
    ; illustris dm_dens:
    start = 2.5
    rots = 1.3
    hue = 1.0
    gamma = 1.0
    
    if n_elements(GA_in) gt 0 then gamma = GA_in
    
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
    
    ; fill rgb_table if requested, or if not, load RGB into display
    if arg_present(rgb_table) then rgb_table = RGB

    ; start above zero on the CT?
    cbot = 0
    if n_elements(bottom) gt 0 then cbot = bottom > 0 < nlev-1
      
    tvlct, RGB, cbot
    
  endif
  
  ; normal IDL tables
  if ctName eq 'bw linear'          then loadct,0,bottom=bottom,/silent
  if ctName eq 'red-temp'           then loadct,3,bottom=bottom,/silent
  if ctName eq 'prism'              then loadct,6,bottom=bottom,/silent
  if ctName eq 'green-white linear' then loadct,8,bottom=bottom,/silent
  if ctName eq 'green-white exp'    then loadct,9,bottom=bottom,/silent
  if ctName eq 'blue-red'           then loadct,11,bottom=bottom,/silent
  if ctName eq 'rainbow'            then loadct,13,bottom=bottom,/silent
  if ctName eq 'purple-red-stripes' then loadct,23,bottom=bottom,/silent
  if ctName eq 'hardcandy'          then loadct,28,bottom=bottom,/silent
  if ctName eq 'nature'             then loadct,29,bottom=bottom,/silent
  if ctName eq 'plasma'             then loadct,32,bottom=bottom,/silent
  if ctName eq 'blue-red2'          then loadct,33,bottom=bottom,/silent
  if ctName eq 'waves'              then loadct,37,bottom=bottom,/silent
  
  ; brewer (nice diverging options, both ends dark, light in middle) (17-26)
  if ctName eq 'brewer-brownpurple' then cgLoadct,18,ncolors=256,bottom=bottom,/brewer
  if ctName eq 'brewer-browngreen'  then cgLoadct,19,ncolors=256,bottom=bottom,/brewer
  if ctName eq 'brewer-purplegreen' then cgLoadct,20,ncolors=256,bottom=bottom,/brewer
  if ctName eq 'brewer-pinkgreen'   then cgLoadct,21,ncolors=256,bottom=bottom,/brewer
  if ctName eq 'brewer-redblue'     then cgLoadct,22,ncolors=256,bottom=bottom,/brewer
  if ctName eq 'brewer-redpurple'   then cgLoadct,25,ncolors=256,bottom=bottom,/brewer
  if ctName eq 'brewer-redgreen'    then cgLoadct,26,ncolors=256,bottom=bottom,/brewer
 
 ; brewer (convergent, light->dark)
 if ctName eq 'brewerC-yellowgreen'     then cgLoadct,0,ncolors=256,bottom=bottom,/brewer
 if ctName eq 'brewerC-yellowblue'      then cgLoadct,1,ncolors=256,bottom=bottom,/brewer
 if ctName eq 'brewerC-yellowgreenblue' then cgLoadct,2,ncolors=256,bottom=bottom,/brewer
 if ctName eq 'brewerC-whitegreen'      then cgLoadct,3,ncolors=256,bottom=bottom,/brewer
 if ctName eq 'brewerC-purplebluegreen' then cgLoadct,4,ncolors=256,bottom=bottom,/brewer
 if ctName eq 'brewerC-purpleblue'      then cgLoadct,5,ncolors=256,bottom=bottom,/brewer
 if ctName eq 'brewerC-whitebluepurple' then cgLoadct,6,ncolors=256,bottom=bottom,/brewer
 if ctName eq 'brewerC-redpurple'       then cgLoadct,7,ncolors=256,bottom=bottom,/brewer
 if ctName eq 'brewerC-whiteredpink'    then cgLoadct,8,ncolors=256,bottom=bottom,/brewer
 if ctName eq 'brewerC-orangered'       then cgLoadct,9,ncolors=256,bottom=bottom,/brewer
 if ctName eq 'brewerC-blues'           then cgLoadct,13,ncolors=256,bottom=bottom,/brewer
 if ctName eq 'brewerC-greens'          then cgLoadct,14,ncolors=256,bottom=bottom,/brewer
 if ctName eq 'brewerC-reds'            then cgLoadct,16,ncolors=256,bottom=bottom,/brewer
 if ctName eq 'brewerC-cubehelix'       then cgLoadCt,28,ncolors=256,bottom=bottom,/brewer
 if ctName eq 'brewerC-cool'            then cgLoadCt,29,ncolors=256,bottom=bottom,/brewer
 
 ; reversed brewer (convergent, dark->light)
 if ctName eq 'brewerR-yellowgreen'     then cgLoadct,0,ncolors=256,bottom=bottom,/brewer,/reverse
 if ctName eq 'brewerR-yellowblue'      then cgLoadct,1,ncolors=256,bottom=bottom,/brewer,/reverse
 if ctName eq 'brewerR-yellowgreenblue' then cgLoadct,2,ncolors=256,bottom=bottom,/brewer,/reverse
 if ctName eq 'brewerR-whitegreen'      then cgLoadct,3,ncolors=256,bottom=bottom,/brewer,/reverse
 if ctName eq 'brewerR-purplebluegreen' then cgLoadct,4,ncolors=256,bottom=bottom,/brewer,/reverse
 if ctName eq 'brewerR-purpleblue'      then cgLoadct,5,ncolors=256,bottom=bottom,/brewer,/reverse
 if ctName eq 'brewerR-whitebluepurple' then cgLoadct,6,ncolors=256,bottom=bottom,/brewer,/reverse
 if ctName eq 'brewerR-redpurple'       then cgLoadct,7,ncolors=256,bottom=bottom,/brewer,/reverse
 if ctName eq 'brewerR-whiteredpink'    then cgLoadct,8,ncolors=256,bottom=bottom,/brewer,/reverse
 if ctName eq 'brewerR-orangered'       then cgLoadct,9,ncolors=256,bottom=bottom,/brewer,/reverse
 if ctName eq 'brewerR-blues'           then cgLoadct,13,ncolors=256,bottom=bottom,/brewer,/reverse
 if ctName eq 'brewerR-greens'          then cgLoadct,14,ncolors=256,bottom=bottom,/brewer,/reverse
 if ctName eq 'brewerR-reds'            then cgLoadct,16,ncolors=256,bottom=bottom,/brewer,/reverse
 if ctName eq 'brewerR-cubehelix'       then cgLoadCt,28,ncolors=256,bottom=bottom,/brewer,/reverse
 if ctName eq 'brewerR-cool'            then cgLoadCt,29,ncolors=256,bottom=bottom,/brewer,/reverse
 
 ; matplotlib
 if ctName eq 'spectral'                then vis_loadct,15,/mpl
 
 ; cpt-city (name specifications must include directory)
 if strpos(ctName,"/") ge 1 and strpos(ctName,"R/") eq -1 then vis_loadct,cpt_filename=ctName
 if strpos(ctName,"R/") ge 1 then vis_loadct,cpt_filename=str_replace(ctName,"R/","/"),/reverse
 
 ; make sure we loaded something
 tvlct,rgb,/get
 if total(rgb,/int) eq 0 then message,'Error: Probably unrecognized CT name.'
 
 ; reverse if requested
 if keyword_set(reverse) then begin
   tvlct,r,g,b,/get
   tvlct,reverse(r),reverse(g),reverse(b)
 endif
  
 ; gamma scaling?
 if keyword_set(GA_in) and ctName ne 'helix' then begin
   tvlct,r,g,b,/get
   nElem = n_elements(r)
   
   r = rebin(r,nElem*100)
   g = rebin(g,nElem*100)
   b = rebin(b,nElem*100)
   
   n = n_elements(r)
   s = long(n*((findgen(n)/n)^GA_in))
   r = r[s] & g = g[s] & b = b[s]
   
   r = rebin(r,nElem)
   g = rebin(g,nElem)
   b = rebin(b,nElem)
   
   tvlct,r,g,b
 endif
 
 if keyword_set(GA_double2) then begin
   ; experimental
   tvlct,r,g,b,/get
   n = n_elements(r)
   
   GA_try = findgen(n)/(n+1) * GA_double2 - 0.5*GA_double2 + 1.0
   GA_try = reverse(GA_try)
   print,GA_try
   
   ss = long(n*((findgen(n)/n)^GA_try))
   r = r[ss] & g = g[ss] & b = b[ss]
   tvlct,r,g,b
 endif

end

; sampleColorTable(): grab a sequence of colors, evenly spaced, from a given colortable

function sampleColorTable, ctname, num, bounds=bounds
  compile_opt idl2, hidden, strictarr, strictarrsubs
  colors = lonarr(num)
  
  set_plot,'ps'
  if ~keyword_set(bounds) then bounds = [0.0,1.0]
  ; get current CT
  tvlct,r,g,b,/get
  
  ; load new colortable
  loadColorTable, ctname
  tvlct,rr,gg,bb,/get
  
  ; determine indices, for each convert color to 24-bit and save
  inds = round( linspace(bounds[0],bounds[1],num) * n_elements(rr) )
  inds = inds > 0 < (n_elements(rr)-1)
  print,inds
  
  foreach ind,inds,i do $
    colors[i] = getColor24( [ rr[ind],gg[ind],bb[ind] ] )

  ; restore CT and return
  tvlct,r,g,b
  
  return,colors

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

; h5_write(): wraps the already simple h5_create

pro h5_write, array, dataset_name, filename

  dataset = { _Name: dataset_name, _TYPE: 'Dataset', _DATA: array }
  h5_create, filename, dataset

end

; mvbakrestarts(): wipe out Arepo restart files with backups

pro mvbakrestarts
  ;message,'remove for safety'
  files = file_search("bak-restart.*")

  foreach file,files,k do begin
    newfname = strmid(file,4,strlen(file)-4)
    cmd = "mv "+file+" "+newfname
    print,cmd
    spawn,cmd
    wait,0.1
  endforeach

end


; postscript output
; -----------------

pro start_PS, filename, xs=xs, ys=ys, eps=eps, big=big, extrabig=extrabig, huge=huge, small=small

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if not keyword_set(xs) then xs=7.5
  if not keyword_set(ys) then ys=5.0
  if n_elements(eps) eq 0 then eps=1
  
  ; make the page bigger
  if n_elements(small) eq 1 then begin
    xs *= 0.9
    ys *= 0.9
  endif
  if n_elements(big) eq 1 then begin
    xs *= 1.2 ;9.0
    ys *= 1.2 ;6.0
  endif
  if n_elements(extrabig) eq 1 then begin
    xs *= 1.4 ;10.5
    ys *= 1.4 ;7
  endif
  if n_elements(huge) eq 1 then begin
    xs *= 2.0 ; 15
    ys *= 2.0 ; 10
  endif

  PS_Start, FILENAME=filename, /nomatch, /quiet, bits_per_pixel=8, color=1, $
            encapsulated=eps, decomposed=0, xs=xs, ys=ys, /inches, font=0, $
            /dejavusans ;requires idl8.2 
            ;font=1,tt_font='Helvetica'
 
  !p.charsize  = 1.1 ; 1.4 for gas accretion paper
  ;!p.charthick = 1.4
  !p.thick    = 4.0 ; 5.0 for gas accretion paper
  
  ; make axis lines thicker
  !x.thick += 1.0
  !y.thick += 1.0
  !z.thick += 1.0
  
  ; reduce number of minor ticks to clean up plots
  !x.minor = 2
  !y.minor = 2
  
  ; make default custom psym (8) a filled circle
  plotsym,0,/fill
  !p.symsize = 0.6
            
end

pro end_PS, pngResize=pngResize, density=density, deletePS=deletePS, im_options=im_options

  compile_opt idl2, hidden, strictarr, strictarrsubs

 ;PNG size=[xs,ys]*300*(resize/100)

  if not keyword_set(pngResize) then $
    PS_End
    
  if keyword_set(pngResize) then $
    PS_End, /PNG, Delete_PS=deletePS, im_options=im_options, Resize=pngResize, density=density, /nofix;, /showcmd

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

; plot_pos(): get some pre-defined plot positioning numbers

function plot_pos, rows=rows, cols=cols, total=total, gap=gap

  compile_opt idl2, hidden, strictarr, strictarrsubs
  if ~keyword_set(total) and ~keyword_set(rows) and ~keyword_set(cols) then $
    message,'Error: Must specify at least one.'
  if keyword_set(total) and (keyword_set(rows) or keyword_set(cols)) then $
    message,'Error: Not both.'
  if (keyword_set(rows) and ~keyword_set(cols)) or (~keyword_set(rows) and keyword_set(cols)) then $
    message,'Error: rows and cols must go together.'
    
  ; only total number of plots set? make a reasonable row/col arrangement
  if keyword_set(total) then begin
    ; one row
    if total eq 1 or total eq 2 or total eq 3 or total eq 5 then begin
      rows = 1 & cols = total
    endif
    
    ; two rows
    if total eq 4 or total eq 6 then begin
      rows = 2 & cols = total/2
    endif
    
    ; three rows
    if total eq 9 or total eq 15 then begin
      rows = 3 & cols = total/3
    endif
    
    ; four rows
    if total eq 8 or total eq 12 then begin
      rows = 4 & cols = total/4
    endif
  endif
  
  if keyword_set(gap) then begin
    ; gap: spacing between all plots such that all axes and axes labels can be drawn
    if rows eq 1 and cols eq 2 then begin
      ; 1x2 (xs=10.5, ys=3.5)
      x0 = 0.10 & x1 = 0.45 & xoff = 0.47
      y0 = 0.16 & y1 = 0.93
                
      pos = list( [x0,y0,x1,y1] ,$ ; left
                  [x0+xoff,y0,x1+xoff,y1] ) ; right
    endif
    
    if rows eq 2 and cols eq 1 then begin
      ; 2x1
      x0 = 0.12 & x1 = 0.96
      y0 = 0.08 & y1 = 0.56 & ysize = 0.40
                
      pos = list( [x0,y1,x1,y1+ysize] ,$ ; top
                  [x0,y0,x1,y0+ysize] ) ; bottom
    endif
    
    if rows eq 1 and cols eq 3 then begin ;xs=9, ys=6
      ; 1x3
      x0 = 0.08 & x1 = 0.34 & xoff = 0.32
      y0 = 0.12 & y1 = 0.94
      
      pos = list( [x0+0*xoff, y0, x1+0*xoff, y1] ,$ ; left
                  [x0+1*xoff, y0, x1+1*xoff, y1] ,$ ; middle
                  [x0+2*xoff, y0, x1+2*xoff, y1] )  ; right
      
    endif
    
    if rows eq 2 and cols eq 2 then begin
      ; 2x2
      x0 = 0.10 & x1 = 0.46 & xoff = 0.49
      y0 = 0.12 & y1 = 0.47 & yoff = 0.48
                
      pos = list( [x0,y0+yoff,x1,y1+yoff] ,$ ; UL
                  [x0+xoff,y0+yoff,x1+xoff,y1+yoff] ,$ ; UR
                  [x0,y0,x1,y1] ,$ ; LL
                  [x0+xoff,y0,x1+xoff,y1] ) ; LR
    endif
    
    if rows eq 4 and cols eq 2 then begin ;xs=9, ys=12
      ; 4x2
      x0 = 0.10 & x1 = 0.46 & xoff = 0.49
      y0 = 0.06 & y1 = 0.24 & yoff = 0.24
      
      pos = list( [x0,y0+3*yoff,x1,y1+3*yoff] ,$ ; top left
                  [x0+xoff,y0+3*yoff,x1+xoff,y1+3*yoff] ,$ ; top right
                  [x0,y0+2*yoff,x1,y1+2*yoff] ,$ ; second row left
                  [x0+xoff,y0+2*yoff,x1+xoff,y1+2*yoff] ,$ ; second row right
                  [x0,y0+1*yoff,x1,y1+1*yoff] ,$ ; third row left
                  [x0+xoff,y0+1*yoff,x1+xoff,y1+1*yoff] ,$ ; third row right
                  [x0,y0,x1,y1] ,$ ; bottom left
                  [x0+xoff,y0,x1+xoff,y1] ) ; bottom right
      
    endif
    
    if rows eq 2 and cols eq 3 then begin ;xs=9, ys=12
      ; 2x3
      x0 = 0.08 & x1 = 0.34 & xoff = 0.32
      y0 = 0.12 & y1 = 0.47 & yoff = 0.48
      
      pos = list( [x0+0*xoff, y0+1*yoff, x1+0*xoff, y1+1*yoff] ,$ ; top left
                  [x0+1*xoff, y0+1*yoff, x1+1*xoff, y1+1*yoff] ,$ ; top middle
                  [x0+2*xoff, y0+1*yoff, x1+2*xoff, y1+1*yoff] ,$ ; top right
                  [x0+0*xoff, y0+0*yoff, x1+0*xoff, y1+0*yoff] ,$ ; bottom left
                  [x0+1*xoff, y0+0*yoff, x1+1*xoff, y1+0*yoff] ,$ ; bottom middle
                  [x0+2*xoff, y0+0*yoff, x1+2*xoff, y1+0*yoff] ) ; bottom right
      
    endif
    
    
    
    if rows eq 1 and cols eq 3 then begin ;xs=9, ys=4 (2d maxtemp histo for 3 runs)
    
    endif
  endif else begin
    ; no gap: boxgrid configuration as in gas accretion paper plots
    if rows eq 4 and cols eq 1 then begin ;xs=6, ys=12
      x0 = 0.14 & x1 = 0.94
      y0 = 0.05 & y1 = 0.27 & yoff = 0.22
      
      pos = list( [x0,y0+3*yoff,x1,y1+3*yoff] ,$ ; top
                  [x0,y0+2*yoff,x1,y1+2*yoff] ,$ ; second row
                  [x0,y0+1*yoff,x1,y1+1*yoff] ,$ ; third row
                  [x0,y0,x1,y1]                ) ; bottom
    endif
    
    if rows eq 2 and cols eq 3 then begin ;/big, e.g. plotValMaxHistos
      x0 = 0.15 & x1 = 0.42 & x2 = 0.69 & x3 = 0.96
      y0 = 0.15 & y1 = 0.55 & y2 = 0.95
      pos = list( [x0,y1,x1,y2] ,$ ; ul
                  [x1,y1,x2,y2] ,$ ; uc
                  [x2,y1,x3,y2] ,$ ; ur
                  [x0,y0,x1,y1] ,$ ; ll
                  [x1,y0,x2,y1] ,$ ; lc
                  [x2,y0,x3,y1] )  ; lr
    endif
  endelse

  if n_elements(pos) eq 0 then message,'Error: Requested plot config not yet set.'
  return, pos
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

  workingPath  = '/n/home07/dnelson/plots/'
  
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
    oB[i]->execute, ".r helper"
    oB[i]->SetVar, "redshift", redshifts[i]
    oB[i]->SetVar, "res", res
    oB[i]->execute, "testBridgePro", /NOWAIT ; asynchronous
  endfor
  
  ; wait for children to finish and cleanup (DESTROY IN REVERSE ORDER see CR64611)
  for i=0,n_elements(redshifts)-1,-1 do $
    while (oB[i]->Status() ne 0) do wait,0.1
  obj_destroy,oB
  
  print,"Elapsed time: "+str(systime(/seconds)-start_time)+" sec"

end

; load routines for use
; ---------------------
@externalC
@units
@simParams
@cosmoUtil
@groupCat
@cosmoLoad

@haloModels
@haloModelsPlot

@mergerTree
@galaxyCat
@galaxyHaloCat
@accretionMode
@accretionTimes
@accretionFlags
;@accretionVel
@maxVals
@timeScales
@timeScalesPlot
@growthRates

;@accretionRates
@accretionTraj
;@accretionTrajVis
@cosmoVis
@cosmoVisMap
@cosmoVisMovie
@cosmoVisStars
;@cosmoVisTry
;@cosmoOverDens
;@haloCompProj

@illustrisVisCutout
@illustrisVis

@sphere
@spherePlot
@spherePlotStat
@spherePowerSpec
;@filamentSearch
;@filamentSearchPlot

@plotGalCat
@plotMaxVals
@plotAccTimes
@binVsHaloMass
@plotVsHaloMass
@shyPlot
@plotRadProfiles

; new feedback project analysis
@plotVsRedshift
@recycledMaterial
@tracksFluid

; zooms:
@ICs_cosmoZoom
@zoomVis
@zoomEvo

@tracersVel_Cosmo
;@tracersVel_Halos
;@tracersVel_Disks
;@tracersVel_2D
;@tracersVel_SphSym

@tracersMC
;@tracersMC_Halos
;@tracersMC_2D
;@tracersMC_SphSym

;@arepoLoad
;@arepoVis2D
;@arepoSphSym

@slurm
@testing
