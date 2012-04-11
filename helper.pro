; helper.pro
; helper functions
; dnelson apr.2012
;
; NOTE: all my IDL routines loaded at bottom of this file

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
  return,fsc_color(units.colors[ind])
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

function shuffle, array, seed=seed
  compile_opt idl2, hidden, strictarr, strictarrsubs
  if n_elements(seed) ne 0 then iseed=seed
  return,array[sort(randomu(iseed,n_elements(array)))]
end

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

; partTypeNum(): convert a string description of a particle type to its numeric value

function partTypeNum, PT

  compile_opt idl2, hidden, strictarr, strictarrsubs

  partType = PT ; so we don't change the input
  if not isnumeric(partType) then partType = strlowcase(str(partType))

  if (strcmp(partType,'gas')       or strcmp(partType,'hydro'))      then partType = 0
  if (strcmp(partType,'dm')        or strcmp(partType,'darkmatter')) then partType = 1
  if (strcmp(partType,'tracervel') or strcmp(partType,'tracersvel')) then partType = 2
  if (strcmp(partType,'tracermc')  or strcmp(partType,'tracersmc'))  then partType = 3
  if (strcmp(partType,'stars')     or strcmp(partType,'star'))       then partType = 4
  
  if (strcmp(partType,'tracer') or strcmp(partType,'tracers')) then $
    message,'ERROR: Please specify which type of tracers!'
  
  if not isnumeric(partType) then $
    message,'ERROR: Unrecognized partType!'
  
  if (partType lt 0 or partType gt 4) then $
    message,'ERROR: partType = ' + str(partType) + ' out of bounds!'

  return, partType
end

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
    gamma = 1.0  ; gamma intensity expontent (1=normal/linear,<1 emphasize low vals,>1 emphasize high vals)
    
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
    if keyword_set(reverse) then RGB = reverse(RGB)

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

end

; basic IO
; --------

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
  
    ; header
    nPts = 0UL
    readU, lun, nPts

    ; replicate
    pts  = replicate(ptStruct,nPts)
    
    ; fill
    readU, lun, pts
  
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

pro start_PS, filename, xs=xs, ys=ys, eps=eps, big=big

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if not keyword_set(xs) then xs=7.5
  if not keyword_set(ys) then ys=5.0
  if n_elements(eps) eq 0 then eps=1
  
  ; make the page bigger
  if n_elements(big) eq 1 then begin
    xs *= 1.2 ;9.0
    ys *= 1.2 ;6.0
  endif 

  PS_Start, FILENAME=filename, /nomatch, /quiet, bits_per_pixel=8, color=1, $
            encapsulated=eps, decomposed=0, xs=xs, ys=ys, /inches, font=0;, tt_font='Helvetica' ;3/2  
 
  !p.charsize  = 1.4
  ;!p.charthick = 1.4
  !p.thick    = 5.0
  
  !x.thick += 1.0
  !y.thick += 1.0
  !z.thick += 1.0
  
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
    PS_End, /PNG, Delete_PS=deletePS, im_options=im_options, Resize=pngResize, /showcmd

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

; getIDIndexMap(): return an array of size max(ids)-min(ids) such that array[ID-min(ids)] is the 
;                  index of the original array ids where ID is found (assumes a one to one mapping, 
;                  not repeated indices as in the case of parentIDs for tracers)

function getIDIndexMap, ids, minid=minid

  compile_opt idl2, hidden, strictarr, strictarrsubs

  minid = long(min(ids))
  maxid = long(max(ids))
  
  if (maxid-minid) gt 2e9 then stop ; should change arr to lon64arr

  ; looped where approach (never a good idea)
  ;arr = l64indgen(maxid-minid+1)
  ;for i=minid,maxid do begin
  ;  w = where(ids eq i,count)
  ;  if (count gt 0) then arr[i] = w[0]
  ;endfor

  ; C-style loop approach (good for sparse IDs)
  arr = ulonarr(maxid-minid+1)
  for i=0UL,n_elements(ids)-1L do arr[ids[i]-minid] = i

  ; reverse histogram approach (good for dense ID sampling, maybe better by factor of ~2)
  ;arr = l64indgen(maxid-minid+1)
  ;h = histogram(ids,rev=rev,omin=omin)
  ;for i=0L,n_elements(h)-1 do if (rev[i+1] gt rev[i]) then arr[i] = rev[rev[i]:rev[i+1]-1]

  return, arr
end

; external C-routine interfaces
; -----------------------------

; calcHSML(): use CalcHSML external C-routine for the tree, neighbor searchings, and smoothing lengths

function calcHSML, Pos, ndims=ndims, nNGB=nNGB, boxSize=boxSize

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if (ndims ne 1 and ndims ne 2 and ndims ne 3) then stop

  ; prepare inputs
  npos = (size(pos))[2]

  NumPart = long(npos)
  Mass    = fltarr(npos)+1.0 ;dummy
  
  DesNumNgb    = long(nNGB) ; number of neighbors to use
  DesNumNgbDev = long(0)
  boxSize      = float(boxSize)
  HsmlGuess    = float(1.0)
  Softening    = float(1.0)
  
  hsml_out = fltarr(NumPart)
  
  ; call CalcHSML
  libName = '/n/home07/dnelson/idl/CalcHSML/CalcHSML_'+str(ndims)+'D.so'
  ret = Call_External(libName, 'CalcHSML', $
                      NumPart,Pos,Mass,DesNumNgb,DesNumNgbDev,boxSize,HsmlGuess,Softening,hsml_out, $
                      /CDECL)
   
  return, hsml_out
                     
end

; calcHSMLds(): use CalcHSMLds external C-routine to calculate the smoothing length needed to locate
;               a given number of neighbors with positions Pos around each position in SearchPos
;               where the two are generally different (e.g. the size of a tophat filter)

function calcHSMLds, Pos, SearchPos, ndims=ndims, nNGB=nNGB, boxSize=boxSize

  compile_opt idl2, hidden, strictarr, strictarrsubs

  npos = size(Pos)
  nsrc = size(SearchPos)

  if ndims ne 1 and ndims ne 2 and ndims ne 3 then message,'Error: Need ndims=1,2,3.'
  if npos[0] ne 2 or npos[1] ne 3 then message,'Error: Point position array shape.'
  if nsrc[0] ne 2 or nsrc[1] ne 3 then message,'Error: Search position array shape.'
  if npos[2] lt nNGB then message,'Error: Point count too low for nNGB.'

  ; prepare inputs
  NumPart   = long(npos[2])
  NumSearch = long(nsrc[2])
  
  DesNumNgb    = long(nNGB)     ; number of neighbors to use
  DesNumNgbDev = long(0)        ; deviation allowed
  boxSize      = float(boxSize) ; use zero for non-periodic search
  
  hsml_out = fltarr(NumSearch)
  
  ; make sure point arrays are float triples since we direct cast now
  Pos = float(Pos)
  SearchPos = float(SearchPos)
  
  ; call CalcHSML
  libName = '/n/home07/dnelson/idl/CalcHSMLds/CalcHSMLds_'+str(ndims)+'D.so'
  ret = Call_External(libName, 'CalcHSMLds', $
                      NumPart,Pos,NumSearch,SearchPos,DesNumNgb,DesNumNgbDev,boxSize,hsml_out, $
                      /CDECL)
                      
  return, hsml_out              
end

; calcNN(): use CalcNN external C-routine for the tree and neighbor search
;           return the index of Pos_SrcTargs closest to each of Pos_SrcOrigs

function calcNN, Pos_SrcTargs, Pos_SrcOrigs, boxSize=boxSize, ndims=ndims

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; prepare inputs
  n_srcTargs = long( n_elements(Pos_SrcTargs[0,*]) )
  n_srcOrigs = long( n_elements(Pos_SrcOrigs[0,*]) )
  
  boxSize   = float(boxSize)
  
  ind_out = lonarr(n_srcOrigs)
  
  ; call CalcNN
  libName = '/n/home07/dnelson/idl/CalcNN/CalcNN_'+str(ndims)+'D.so'
  ret = Call_External(libName, 'CalcNN', $
                      n_srcTargs,n_srcOrigs,Pos_SrcTargs,Pos_SrcOrigs,boxSize,ind_out, $
                      /CDECL)
    
  return, ind_out
end

; calcSphMap(): use CalcSphMap external C-routine to calculate a map of projected densities
;               with the sph spline kernel

function calcSphMap, pos, hsml, mass, axes=axes, boxSize=boxSize, boxCen=boxCen, $
                     nPixels=nPixels, ndims=ndims

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; prepare inputs
  npos = (size(pos))[2]
  
  boxSize  = float(boxSize)
  boxCen   = float(boxCen)
  nPixels  = fix(nPixels)
  axes     = fix(axes)
  
  mode     = fix(1)
  periodic = fix(1)
  
  ; check inputs
  if (n_elements(boxSize) ne 3) then stop
  if (n_elements(boxCen)  ne 3) then stop
  if (n_elements(nPixels) ne 2) then stop
  if (n_elements(axes)    ne 2) then stop
  
  if (mode ne 1 and mode ne 2 and mode ne 3) then stop
  if (periodic ne 0 and periodic ne 1) then stop
  
  if (axes[0] ne 0 and axes[0] ne 1 and axes[0] ne 2) then stop
  if (axes[1] ne 0 and axes[1] ne 1 and axes[1] ne 2) then stop
  
  ; make return
  grid_out = fltarr(nPixels[0],nPixels[1])
  
  ; call CalcSphMap
  libName = '/n/home07/dnelson/idl/CalcSphMap/CalcSphMap_'+str(ndims)+'D.so'
  ret = Call_External(libName, 'CalcSphMap', $
                      npos,pos,hsml,mass,grid_out,boxSize[0],boxSize[1],boxSize[2],$
                      boxCen[0],boxCen[1],boxCen[2],axes[0],axes[1],nPixels[0],nPixels[1],$
                      mode,periodic,/CDECL)

  return, grid_out
end

; estimateDensityTophat(): spatial density estimator for an input position tuple array and
;                          CONSTANT mass per particle, by using HSMLs as tophat filter sizes
; pos_search : if specified, estimate density at this point set not at the positions of the particles

function estimateDensityTophat, pos, pos_search=pos_search, mass=mass, $
                                ndims=ndims, nNGB=nNGB, boxSize=boxSize

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if ~keyword_set(nNGB) or ~keyword_set(ndims) then message,'Error: Must specify nNGB and ndims.'
  if ~keyword_set(mass) or n_elements(mass) ne 1 then message,'Error: Expected one constant mass.'
  
  ; calculate smoothing lengths
  if ~keyword_set(pos_search) then $
    hsml_out = calcHSML(pos,ndims=ndims,nNGB=nNGB,boxSize=boxSize)
  if keyword_set(pos_search) then $
    hsml_out = calcHSMLds(pos,pos_search,ndims=ndims,nNGB=nNGB,boxSize=boxSize)
  
  ; convert smoothing lengths to densities
  if (ndims eq 1) then hsml_out = mass * nNGB / (2.0 * temporary(hsml_out))
  if (ndims eq 2) then hsml_out = mass * nNGB / (!pi * temporary(hsml_out)^2.0)
  if (ndims eq 3) then hsml_out = mass * nNGB / (4.0*!pi/3.0 * temporary(hsml_out)^3.0)
  
  return,hsml_out
end

; load routines for use
; ---------------------
@units
@simParams
@cosmoUtil
@cosmoLoad
@LSF

@mergerTree
@cosmoAnalysis
@cosmoHist
@accretionTraj
@accretionTrajVis
@cosmoVis
@cosmoSphere
@cosmoSpherePlot
@cosmoPlotGalCat
@cosmoPlot

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

