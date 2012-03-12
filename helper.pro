; helper.pro
; helper functions
; dnelson mar.2012
;
; NOTE: all my IDL routines loaded at bottom of this file

; getUnits(): return a structure of useful units

function getUnits

  Hubble  = 1.0      ;H0 in 100km/s/Mpc
  Gravity = 6.673e-8 ;G in cgs, cm^3/g/s^2

  units = { units,                   $
  
            ; units (from parameter file)
            UnitLength_in_cm         : double(3.085678e21)    ,$;  1.0 kpc
            UnitMass_in_g            : 1.989*double(10.0)^43  ,$;  1.0e10 solar masses
            UnitVelocity_in_cm_per_s : double(1.0e5)          ,$;  1 km/sec
            
            ; derived units
            UnitTime_in_s       : 0.0D                        ,$
            UnitDensity_in_cgs  : 0.0D                        ,$
            UnitPressure_in_cgs : 0.0D                        ,$
            UnitEnergy_in_cgs   : 0.0D                        ,$
            
            ; non-cgs units
            UnitMass_in_Msun    : 0.0D                        ,$
            
            ; constants
            boltzmann   : double(1.38066e-16)                 ,$ ;cgs
            mass_proton : double(1.6727e-24)                  ,$ ;cgs
            
            ; derived constants
            H0      : 0.0D                                    ,$
            G       : 0.0D                                    ,$
            rhoCrit : 0.0D                                    ,$
            
            ; color list
            colors : strarr(17)                               ,$
            
            ; unit conversions
            s_in_Myr  : 3.156e13                              ,$
            Msun_in_g : 1.989*double(10.0)^33                 ,$
            pc_in_cm  : 3.0868e18                             ,$
            Mpc_in_cm : 3.0868e24                             ,$
            kpc_in_km : 3.0856e16                              $
      }
      
  ; derived units
  units.UnitTime_in_s       = units.UnitLength_in_cm / units.UnitVelocity_in_cm_per_s
  units.UnitDensity_in_cgs  = units.UnitMass_in_g / units.UnitLength_in_cm^3.0
  units.UnitPressure_in_cgs = units.UnitMass_in_g / units.UnitLength_in_cm / units.UnitTime_in_s^2.0
  units.UnitEnergy_in_cgs   = units.UnitMass_in_g * units.UnitLength_in_cm^2.0 / units.UnitTime_in_s^2.0
  
  ; non-cgs units
  units.UnitMass_in_Msun = units.UnitMass_in_g / units.Msun_in_g
  
  ; derived constants (in code units)
  units.H0 = Hubble * 100 * 1e5 / (units.Mpc_in_cm) / $
             units.UnitVelocity_in_cm_per_s * units.UnitLength_in_cm
  units.G  = Gravity / units.UnitLength_in_cm^3.0 * units.UnitMass_in_g * units.UnitTime_in_s^2.0
  
  units.rhoCrit = 3.0 * units.H0^2.0 / (8.0*!pi*units.G) ;code, z=0

  ; color list
  units.colors = ['black','blue','green','red','cyan','magenta','gray','orange', $
                  'brown','chartreuse','violet','papaya','aquamarine', $
                  'firebrick', 'rosy brown', 'gold', 'olive']

  return, units
end

; loadCSV()

function loadCSV, headerLines, fileName, ptStruct, header=header;, format=format

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

; loadBinary()
; right now just reads Stars_X.bin
; first float indicates how many pts

function loadBinary, fileName, ptStruct

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

; loadBinarySequence()
; right now just reads Stars_X_Y where X=num, Y=node

function loadBinarySequence, fileBase, ptStruct

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

; one line utility functions
; --------------------------

function str, tString
  return, strcompress(string(tString),/remove_all)
end

function isnumeric, input
  on_ioerror, false
  test = double(input)
  return, 1
  false: return, 0
end

function getColor, i, name=name
  units = getUnits()
  ind = (i) mod (n_elements(units.colors)-1)
  
  if keyword_set(name) then return,units.colors[ind]
  return,fsc_color(units.colors[ind])
end

function linspace, a, b, N
  vals = findgen(N) / (N-1.0) * (b-a) + a
  return, vals
end

function logspace, a, b, N, mid=mid
  vals = findgen(N) / (N-1.0) * (b-a) + a
  
  ; return mid-bin points instead
  if keyword_set(mid) then $
    vals = (findgen(N-1)+0.5) / (N-1.0) * (b-a) + a
  
  vals = 10.0^vals
  
  return, vals
end

function nuniq, arr
  return, n_elements(uniq(arr,sort(arr)))
end

function shuffle, array, seed=seed
  if n_elements(seed) ne 0 then iseed=seed
  return,array[sort(randomu(iseed,n_elements(array)))]
end

; postscript output
; -----------------

pro start_PS, filename, xs=xs, ys=ys, eps=eps, big=big

  if not keyword_set(xs) then xs=7.5
  if not keyword_set(ys) then ys=5.0
  if n_elements(eps) eq 0 then eps=1
  
  ; make the page bigger
  if n_elements(big) eq 1 then begin
    xs *= 1.2 ;9.0
    ys *= 1.2 ;6.0
  endif 

  PS_Start, FILENAME=filename, /nomatch, /quiet, bits_per_pixel=8, color=1, $
            encapsulated=eps, decomposed=0, xs=xs, ys=ys, /inches, font=0, tt_font='Times' ;3/2  
 
  !p.charsize = 1.4
  !p.thick    = 5.0
  
  !x.thick += 2.0
  !y.thick += 2.0
  !z.thick += 2.0
  
  ; make default custom psym (8) a filled circle
  plotsym,0,/fill
  !p.symsize = 0.7
            
end

pro end_PS, pngResize=pngResize, deletePS=deletePS

 ;PNG size=[xs,ys]*300*(resize/100)

  if not keyword_set(pngResize) then $
    PS_End
    
  if keyword_set(pngResize) then $
    PS_End, /PNG, Delete_PS=deletePS, Resize=pngResize

end

; save_eps(): idl 8.x compatible EPS save

pro save_eps, p, plotName, width=width, height=height, savePDF=savePDF, savePNG=savePNG

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
  print,'!CPU.HW_NCPU = ' + str(!CPU.HW_NCPU) + ' TPool_NThreads = ' + str(!CPU.TPOOL_NTHREADS)
end

; testBridgePro

pro testBridgePro

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
  
  if (list.Count() eq 0) then return,[0]
  
  arr = []
  for i=0ULL,list.Count()-1 do $
    arr = [arr,list[i]]

  return,arr
  
end

; removeIntersectionFromB(): return a modified version of B with all those elements also found in
;                            A (the collision/intersection) removed

function removeIntersectionFromB, A, B, union=union

    match, A, B, A_ind, B_ind, count=count
    
    A_ind = !NULL ;unused
    
    if (count gt 0) then begin
      ; remove B[B_ind] using complement
      all = bytarr(n_elements(B))
      if (B_ind[0] ne -1L) then all[B_ind] = 1B
      w = where(all eq 0B, ncomp)
    
      if (ncomp ne n_elements(B)-count) then begin
        print,'removeIntersectionFromB: ERROR ',ncomp,n_elements(B),count
        return,0
      endif
      
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

  minid = long(min(ids))
  maxid = long(max(ids))

  ; looped where approach (never a good idea)
  ;arr = l64indgen(maxid-minid+1)
  ;for i=minid,maxid do begin
  ;  w = where(ids eq i,count)
  ;  if (count gt 0) then arr[i] = w[0]
  ;endfor

  ; C-style loop approach (good for sparse IDs)
  arr = lon64arr(maxid-minid+1)
  for i=0L,n_elements(ids)-1L do arr[ids[i]-minid] = i

  ; reverse histogram approach (good for dense ID sampling, maybe better by factor of ~2)
  ;arr = l64indgen(maxid-minid+1)
  ;h = histogram(ids,rev=rev,omin=omin)
  ;for i=0L,n_elements(h)-1 do if (rev[i+1] gt rev[i]) then arr[i] = rev[rev[i]:rev[i+1]-1]

  return, arr
end

; external C-routine interfaces
; -----------------------------
; 
; calcHSML(): use CalcHSML external C-routine for the tree, neighbor searchings, and smoothing lengths

function calcHSML, Pos, ndims=ndims, nNGB=nNGB, boxSize=boxSize

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

; calcNN(): use CalcNN external C-routine for the tree and neighbor search
;           return the index of Pos_SrcTargs closest to each of Pos_SrcOrigs

function calcNN, Pos_SrcTargs, Pos_SrcOrigs, boxSize=boxSize, ndims=ndims

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
;                          mass array for some particle type
;   do density calculation by calculating smoothing lengths for all the particles
;   

function estimateDensityTophat, pos, mass=mass, ndims=ndims, nNGB=nNGB, boxSize=boxSize

  if (not keyword_set(nNGB) or not keyword_set(ndims)) then stop

  ; calculate smoothing lengths
  hsml_out = calcHSML(pos,ndims=ndims,nNGB=nNGB,boxSize=boxSize)
                      
  ; estimate densities on eval_pos using hsml
  if (ndims eq 1) then $
    eval_dens = nNGB / (2.0 * hsml_out)
    
  if (ndims eq 2) then $
    eval_dens = nNGB / (!pi * hsml_out^2.0)
    
  if (ndims eq 3) then $
    eval_dens = nNGB / (4.0*!pi/3.0 * hsml_out^3.0)
  
  ; add in mass
  if keyword_set(mass) then $
    eval_dens *= mass[0]
  
  return,eval_dens

end

; load routines for use
; ---------------------
@simParams
@cosmoUtil
@cosmoLoad

@tracersMC
@tracersMC_2D
@tracersMC_SphSym
@tracersMC_Halos

;@tracers
;@tracersCosmo
;@tracersCosmoHalos
;@tracersDisks
;@tracersShocktube
;@tracersSpheres

@cosmoVis
@cosmoSphere
;@cosmoPlot
;@cosmoAnalysis

@arepoLoad
@arepoVis2D
@arepoSphSym
