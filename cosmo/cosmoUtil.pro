; cosmoUtil.pro
; cosmological simulations - utility functions
; dnelson nov.2012

; redshiftToSnapNum(): convert redshift to the nearest snapshot number

function redshiftToSnapNum, redshiftList, sP=sP, verbose=verbose

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if not keyword_set(verbose) then verbose = 0
  
  saveFileName = sP.derivPath + sP.savPrefix + '_snapnum.redshift.sav'

  if not (file_test(saveFileName)) then begin

    nSnaps = 400 ; maximum
    
    redshifts = fltarr(nSnaps) - 1
    times     = fltarr(nSnaps) - 1
    
    for m=0,nSnaps-1 do begin
      ; format filename
      ext = str(string(m,format='(I3.3)'))
      f = sP.simPath + 'snapdir_' + ext + '/snap_' + ext + '.0.hdf5'
      
      ; single file per group catalog
      if (not file_test(f)) then $
        f = sP.simPath + 'snap_' + ext + '.hdf5'
        
      ; single groupordered file per group catalog
      if (not file_test(f)) then $
        f = sP.simPath + 'snap-groupordered_' + ext + '.hdf5'
      
      ; if file doesn't exist yet, skip (e.g. every other file deleted)
      if (not file_test(f)) then continue
    
      ; load hdf5 header and save time+redshift
      fileID   = h5f_open(f)
      headerID = h5g_open(fileID,"Header")
      
      redshifts[m] = h5a_read(h5a_open_name(headerID,"Redshift"))
      times[m]     = h5a_read(h5a_open_name(headerID,"Time"))
      
      h5g_close, headerID
      h5f_close, fileID
    endfor
  
    ; save/restore
    save,nSnaps,redshifts,times,filename=saveFileName
  endif else begin
    restore,saveFileName
  endelse

  ; for -1 return !NULL
  if (total(redshiftList) eq -1) then return,!NULL
  
  ; return array
  snapNum = intarr(n_elements(redshiftList))
  
  foreach redshift,redshiftList,i do begin

    dists = abs(redshifts - redshift)
    w = where(dists eq min(dists),count)
    if count eq 2 then w = w[0]
    if count eq 0 then stop
    if min(dists) gt 0.1 then print,'Warning! Snapshot selected with redshift error = ',min(dists)
    snapNum[i] = w[0]
    
    if (verbose) then $
      print,'Found nearest snapshot to z = ' + str(redshift) + ' (num = ' + str(snapNum[i]) + ')'
  endforeach

  if (n_elements(snapNum) eq 1) then snapNum = snapNum[0]

  return,snapNum
end

; mergerTreeChild(): construct a Child array given a Parent array
; 
; childPrev : if specified, compose the two mappings such that the returned Child points not at the
;             immediate children of Parent but to a prior child of Parent (say, the first)

function mergerTreeChild, Parent, ChildPrev=ChildPrev

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  Child = lonarr(max(Parent)+1)-1
  for i=0UL,n_elements(Parent)-1L do if Parent[i] ne -1 then Child[Parent[i]] = i

  ; compose with previous Child mapping?
  if n_elements(ChildPrev) gt 0 then begin
    ChildNew = lonarr(n_elements(Child))-1
    
    ; for each immediate Child, replace Cur->Child pointer by Cur->Child->PrevChild pointer
    for i=0UL,n_elements(Child)-1 do begin
      ; discard pointers to children greater than the maximum in ChildPrev, which can happen when
      ; the Child pointed to is not part of the halo set being tracked (the composed mapping is
      ; undefined for this Child)
      if Child[i] lt n_elements(ChildPrev) then $
        if Child[i] ne -1 then ChildNew[i] = ChildPrev[Child[i]] ; only follow valid Child links
    endfor
    
    return, ChildNew
  endif

  return, Child
end

; correctPeriodicDistVecs(): enforce periodic B.C. for distance vecotrs (effectively component by 
;                            component), input vecs in format fltarr[3,n]

pro correctPeriodicDistVecs, vecs, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  w = where(vecs gt sP.boxSize*0.5,count)
  if (count ne 0) then $
    vecs[w] = vecs[w] - sP.boxSize
    
  w = where(vecs lt -sP.boxSize*0.5,count)
  if (count ne 0) then $
    vecs[w] = sP.boxSize + vecs[w]

end

; periodicDists(): calculate distances correctly taking into account periodic B.C.
; 
; if pt is one point: distance from pt to all vecs
; if pt is several points: distance from each pt to each vec (must have same number of points)

function periodicDists, pt, vecs, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  nDimsPt = (size(pt))[0]

  if ( ((size(vecs))[0] ne 1 and (size(vecs))[0] ne 2) or $
       (size(vecs))[1] ne 3) then message,'Error: vecs not in expected shape'
  if (nDimsPt ne 1 and nDimsPt ne 2) then message,'Error: something is wrong'
  if not keyword_set(sP) then message,'Error: need sP for boxSize'

  ; distances from one point to a vector of other points [3,n]
  if (nDimsPt eq 1) then begin
    xDist = vecs[0,*]-pt[0]
    yDist = vecs[1,*]-pt[1]
    zDist = vecs[2,*]-pt[2]
  endif
  
  ; distances from a vector of points [3,n] to another vector of other points [3,n]
  if (nDimsPt eq 2) then begin
    xDist = vecs[0,*]-pt[0,*]
    yDist = vecs[1,*]-pt[1,*]
    zDist = vecs[2,*]-pt[2,*]
  endif
  
  correctPeriodicDistVecs, xDist, sP=sP
  correctPeriodicDistVecs, yDist, sP=sP
  correctPeriodicDistVecs, zDist, sP=sP
  
  dists = reform( sqrt( xDist*xDist + yDist*yDist + zDist*zDist ) )
  
  return, dists

end

; periodicPairwiseDists(): calculate pairwise distances between all 3D points, correctly
;                          taking into account periodic B.C.

function periodicPairwiseDists, pts, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; construct indices for vector computation of pairwise distances
  dims = size(pts,/dimensions)
  
  num = dims[1]*(dims[1]-1)/2
  
  ii = 0L
  index0 = lindgen(dims[1]-1) + 1
  index1 = lonarr(num,/nozero)
  index2 = lonarr(num,/nozero)
  
  for i=0,dims[1]-2 do begin
    n1 = dims[1] - (i+1)
    index1[ii:ii+n1-1] = i
    index2[ii] = index0[0:n1-1] + i
    ii += n1
  endfor
  
  ; component wise difference
  xDist = pts[0,index1] - pts[0,index2]
  yDist = pts[1,index1] - pts[1,index2]
  zDist = pts[2,index1] - pts[2,index2]
  
  ; account for periodic distance function
  correctPeriodicDistVecs, xDist, sP=sP
  correctPeriodicDistVecs, yDist, sP=sP
  correctPeriodicDistVecs, zDist, sP=sP
  
  dists = reform( sqrt( xDist*xDist + yDist*yDist + zDist*zDist ) )

  return,dists

end

; snapNumToRedshift(): convert snapshot number to redshift or time (scale factor)

function snapNumToRedshift, time=time, all=all, sP=sP, snap=snap

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if not keyword_set(sP) then message,'Error: Need sP to convert snapshot number to redshift!'

  saveFileName = sP.derivPath + sP.savPrefix + '_snapnum.redshift.sav'

  if not file_test(saveFileName) then message,'Error: Failed to find snapnum save file.'
  if n_elements(snap) eq 0 then snap = sP.snap
  if snap[0] eq -1 then message,'Error: sP.snap used but was undefined.'
  
  ; restore
  restore,saveFilename

  if (not keyword_set(time)) then begin
    if (keyword_set(all)) then return,redshifts
    
    if (snap[0] ge 0 and snap[0] lt n_elements(redshifts)) then $
      return,redshifts[snap]
  endif
  
  ; keyword_set(time)
  if (keyword_set(all)) then return,times
      
  if (snap[0] ge 0 and snap[0] lt n_elements(redshifts)) then $
    return,times[snap]

end

; snapNumToAge(): convert snapshot number to approximate age of the universe at that point

function snapNumToAge, snap=snap, sP=sP

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if not keyword_set(snap) then snap = sP.snap
  z = snapNumToRedshift(snap=snap,sP=sP)
  return, redshiftToAge(z)
end

function snapNumToAgeFlat, snap=snap, sP=sP
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if not keyword_set(snap) then snap = sP.snap
  z = snapNumToRedshift(snap=snap,sP=sP)
  return, redshiftToAgeFlat(z)
end
  
; partTypeNum(): convert a string description of a particle type to its numeric value

function partTypeNum, PT

  compile_opt idl2, hidden, strictarr, strictarrsubs

  partType = PT ; so we don't change the input
  if not isnumeric(partType) then partType = strlowcase(str(partType))

  if (strcmp(partType,'gas')        or strcmp(partType,'hydro'))      then partType = 0
  if (strcmp(partType,'dm')         or strcmp(partType,'darkmatter')) then partType = 1
  if (strcmp(partType,'tracervel')  or strcmp(partType,'tracersvel')) then partType = 2
  if (strcmp(partType,'tracermc')   or strcmp(partType,'tracersmc'))  then partType = 3
  if (strcmp(partType,'stars')      or strcmp(partType,'star'))       then partType = 4 ;StellarFormationTime>0
  if (strcmp(partType,'wind')       or strcmp(partType,'windstars'))  then partType = 4 ;StellarFormationTime<0
  if (strcmp(partType,'blackholes') or strcmp(partType,'bh') $
                                    or strcmp(partType,'bhs'))        then partType = 5
  
  if (strcmp(partType,'tracer') or strcmp(partType,'tracers')) then $
    message,'ERROR: Please specify which type of tracers!'
  
  if not isnumeric(partType) then $
    message,'ERROR: Unrecognized partType!'
  
  if (partType lt 0 or partType gt 5) then $
    message,'ERROR: partType = ' + str(partType) + ' out of bounds!'

  return, partType
end  
  
; rhoTHisto(): make mass-weighted density-temperature 2d histogram

function rhoTHisto, dens_in, temp_in, mass=mass, nbins=nbins, plot=plot

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; config
  if not keyword_set(nbins) then nbins = 20
  dens = alog10(rhoRatioToCrit(dens_in))
  temp = alog10(temp_in)
  rMinMax = [-2.0,8.0]
  tMinMax = [3.0,7.0]
  
  ; calculate bin sizes
  binSizeRho  = (rMinMax[1]-rMinMax[0]) / nbins
  binSizeTemp = (tMinMax[1]-tMinMax[0]) / nbins
  
  ; if not mass weighting assign equal weights
  if not keyword_set(mass) then mass = replicate(1.0,n_elements(dens))
  
  ; use hist_nd with centered bins and boundary fixes
  h2rt = hist_nd_weight(transpose([[dens],[temp]]),weight=mass,[binSizeRho,binSizeTemp],$
                        min=[rMinMax[0]-binSizeRho*0.5,tMinMax[0]-binSizeTemp*0.5],$
                        max=[rMinMax[1]+binSizeRho*0.49,tMinMax[1]+binSizeTemp*0.49])
  
  ; plot
  if keyword_set(plot) then begin
    ; color table
    loadct, 2, bottom=1, /silent ;bw linear=0, white-green exp=9 (33=blue-red)
      
    tvim,h2rt,pcharsize=!p.charsize-1.0,scale=1,clip=[10,100],$;,/c_map
         xtitle="log ("+textoidl("\rho / \rho_{crit}")+")",ytitle="log (T [K])",$
         stitle="Total Mass ("+textoidl("M_{sun}")+")",barwidth=0.5,lcharsize=!p.charsize-1.5,$
         xrange=[-2.0,8.0],yrange=[3.0,7.0],$;xrange=rMinMax,yrange=tMinMax,$
         /rct;,nodata=0,rgb_nodata=[1.0,1.0,1.0] ;display zeros as white not black
  endif
  
  return,h2rt

end

; redshift_axis(): draw redshift axis

pro redshift_axis, xRange, yRange, ylog=ylog, sP=sP, dotted=dotted, zTicknames=zTicknames

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  logFac = 1.0
  if keyword_set(ylog) then logFac = 10.0
  
  yTickLen = (yRange[1]-yRange[0]) / 40.0 * logFac
  yTextOff = (yRange[1]-yRange[0]) / 50.0 * logFac

  if (not keyword_set(zTicknames)) then $
    zTicknames = ['30','6','4','3','2','1','0.5','0.25','0']
  nZ = n_elements(zTicknames)
  
  zXPos = fltarr(nZ)
  
  ; plot "maximum" redshift label
  if (xRange[0] le 0.0) then $
    fsc_text,0.0,yRange[1]+yTextOff,zTicknames[0],alignment=0.5
  
  ; skip z=30 (highest) at t=0
  for i=1,nZ-1 do begin
    if (sP.snap ne -1) then begin ;x-axis in snapshot number
      zXPos[i] = redshiftToSnapNum(float(zTicknames[i]),sP=sP)
    endif else begin ;x-axis in time [gyr]
      zXPos[i] = redshiftToAge(float(zTicknames[i]))
    endelse
    
    ; plot tick mark and label if inside plotrange
    if (zXPos[i] ge xRange[0] and zXPos[i] le xRange[1]) then begin
      fsc_plot,[zXPos[i],zXPos[i]],[yRange[1],yRange[1]-yTickLen],/overplot
      fsc_text,zXPos[i],yRange[1]+yTextOff,zTicknames[i],alignment=0.5
    endif
    
    ; plot vertical dotted line at each redshift mark if requested
    if keyword_set(dotted) then $
      fsc_plot,[zXPos[i],zXPos[i]],yRange,line=dotted,/overplot,thick=!p.thick-0.5
  endfor
  
  fsc_plot,xRange,[yRange[1],yRange[1]],/overplot
  fsc_text,0.5,0.94,"Redshift",/normal

end

; universeage_axis(): draw age of universe (elapsed) axis

pro universeage_axis, xRange, yRange, ylog=ylog, dotted=dotted

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  logFac = 1.0
  if keyword_set(ylog) then logFac = 10.0
  
  yTickLen = (yRange[1]-yRange[0]) / 40.0 * logFac
  yTextOff = (yRange[1]-yRange[0]) / 50.0 * logFac

  if (not keyword_set(zTicknames)) then begin
    zTicknames = ['0.1','1','1.5','2','3','4','5','6','8','13.8'] ;10=0.33
    zRedshifts = [30,5.7,4.21,3.25,2.23,1.65,1.265,0.98,0.595,0.0]
  endif
  nZ = n_elements(zTicknames)
  
  zXPos = fltarr(nZ)
  
  ; plot "maximum" redshift label
  if (xRange[0] le 0.0) then $
    fsc_text,0.0,yRange[1]+yTextOff,zTicknames[0],alignment=0.5
  
  ; skip t=0 (highest) at z=inf (?)
  for i=1,nZ-1 do begin

    zXPos[i] = zRedshifts[i]
    
    ; plot tick mark and label if inside plotrange
    if (zXPos[i] le xRange[0] and zXPos[i] ge xRange[1]) then begin
      fsc_plot,[zXPos[i],zXPos[i]],[yRange[1],yRange[1]-yTickLen],/overplot
      fsc_text,zXPos[i],yRange[1]+yTextOff,zTicknames[i],alignment=0.5
    endif
    
    ; plot vertical dotted line at each redshift mark if requested
    if keyword_set(dotted) then $
      fsc_plot,[zXPos[i],zXPos[i]],yRange,line=dotted,/overplot,thick=!p.thick-0.5
  endfor
  
  fsc_plot,xRange,[yRange[1],yRange[1]],/overplot
  fsc_text,(xRange[0]-xRange[1])/2.0,yRange[1]+yTextOff*5.0,"Time [Gyr]",alignment=0.5

end

; exportParticlesAscii(): export some part of a snapshot as text

pro exportParticlesAscii

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function simParams, loadSnapshotSubset

  partType = 'gas'
  fileName = '128fb.gas.txt'
  sP = simParams(res=128,run='feedback',redshift=0.0)
  
  pos  = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
  dens = loadSnapshotSubset(sP=sP,partType=partType,field='dens')
  u    = loadSnapshotSubset(sP=sP,partType=partType,field='u')
  
  print,n_elements(u)
  stop
  
  outBuf = fltarr(5,n_elements(u))
  outBuf[0,*] = pos[0,*]
  outBuf[1,*] = pos[1,*]
  outBuf[2,*] = pos[2,*]
  outBuf[3,*] = dens / mean(dens)
  outBuf[4,*] = u
  
  openW,lun,fileName,/get_lun
  
    ;for i=0,n_elements(u)-1 do $
      printf,lun,outBuf,format='(f10.4,1x,f10.4,1x,f10.4,1x,f14.8,1x,f10.1)'
  
  close,lun
  free_lun,lun

end
