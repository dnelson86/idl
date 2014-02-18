; cosmoUtil.pro
; cosmological simulations - utility functions
; dnelson aug.2013

; redshiftToSnapNum(): convert redshift to the nearest snapshot number

function redshiftToSnapNum, redshiftList, sP=sP, verbose=verbose, subBox=subBox

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if not keyword_set(verbose) then verbose = 0
  if ~keyword_set(sP) then message,'Error: Must input sP.'
  if n_elements(redshiftList) eq 0 then redshiftList = [sP.redshift]
  
  if keyword_set(subBox) then sbstr = 'subbox0_' else sbstr = ''
  if keyword_set(subBox) then sbst2 = 'subbox0/' else sbst2 = ''
  
  saveFileName = sP.derivPath + sP.savPrefix + '_' + sbstr + 'snapnum.redshift.sav'

  if not (file_test(saveFileName)) then begin

    nSnaps = 400 ; maximum
    if keyword_set(subBox) then nSnaps *= 10
    
    redshifts = fltarr(nSnaps) - 1
    times     = fltarr(nSnaps) - 1
    
    for m=0,nSnaps-1 do begin
      ; format snapshot number
      if (m le 999) then $
        ext = string(m,format='(I3.3)')
      if (m gt 999) then $
        ext = string(m,format='(I4.4)')
        
      ; format filename
      f = sP.simPath + sbst2 + 'snapdir_' + sbstr + ext + '/snap_' + sbstr + ext + '.0.hdf5'
      
      ; single file per group catalog
      if (not file_test(f)) then $
        f = sP.simPath + sbst2 + 'snap_' + sbstr + ext + '.hdf5'
        
      ; single groupordered file per group catalog
      if (not file_test(f)) then $
        f = sP.simPath + sbst2 + 'snap-groupordered_' + ext + '.hdf5'
      
      ; if file doesn't exist yet, skip (e.g. every other file deleted)
      if (not file_test(f)) then continue
      
      ; corrupt/non-HDF5?
      if ~h5f_is_hdf5(f) then begin
        print,'CORRUPT, SKIPPING: ',f
        continue
      endif
    
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

; gasMassesFromIDs(): return individual gas cell/particle masses given input ID list
  
function gasMassesFromIDs, search_ids, sP=sP

  h = loadSnapshotHeader(sP=sP)

  if sP.trMCPerCell eq 0 then begin
    ; SPH case: all particles have constant mass
    massPerPart = sP.targetGasMass
    masses = replicate(massPerPart,h.nPartTot[partTypeNum('gas')])
  endif else begin
    ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    
    idIndexMap = getIDIndexMap(ids,minid=minid)
    ids = !NULL
    ids_ind = idIndexMap[search_ids-minid]
    idIndexMap = !NULL
    
    masses = loadSnapshotSubset(sP=sP,partType='gas',field='mass',inds=ids_ind)
  endelse
  
  return, masses
  
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

; correctPeriodicDistVecs(): enforce periodic B.C. for distance vectors (effectively component by 
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

; correctPeriodicPosVecs(): enforce periodic B.C. for positions (add boxSize to any negative
;                           points, subtract boxSize from any points outside boxSize)

pro correctPeriodicPosVecs, vecs, boxSize=boxSize, sP=sP
  compile_opt idl2, hidden, strictarr, strictarrsubs
  if n_elements(boxSize) eq 0 then boxSize = sP.boxSize
  
  w = where(vecs lt 0.0, count)
  if count ne 0 then vecs[w] += boxSize
  
  w = where(vecs gt boxSize,count)
  if count ne 0 then vecs[w] -= boxSize

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

function snapNumToRedshift, time=time, all=all, sP=sP, snap=snap, subBox=subBox

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if not keyword_set(sP) then message,'Error: Need sP to convert snapshot number to redshift!'

  sbstr = ''
  if keyword_set(subBox) then sbstr = 'subbox0_'  
  
  saveFileName = sP.derivPath + sP.savPrefix + '_' + sbstr + 'snapnum.redshift.sav'

  if not file_test(saveFileName) then dummy = redshiftToSnapNum(sP=sP, subBox=subBox)
  if n_elements(snap) eq 0 then snap = sP.snap
  if snap[0] eq -1 and ~keyword_set(all) then message,'Error: sP.snap used but was undefined.'
  
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

  ; non-zoom, with or without GFM
  if (strcmp(partType,'gas')        or strcmp(partType,'hydro'))      then partType = 0
  if (strcmp(partType,'dm')         or strcmp(partType,'darkmatter')) then partType = 1
  if (strcmp(partType,'tracervel')  or strcmp(partType,'tracersvel')) then partType = 2
  if (strcmp(partType,'tracermc')   or strcmp(partType,'tracersmc'))  then partType = 3
  if (strcmp(partType,'stars')      or strcmp(partType,'star'))       then partType = 4 ;StellarFormationTime>0
  if (strcmp(partType,'wind')       or strcmp(partType,'windstars'))  then partType = 4 ;StellarFormationTime<0
  if (strcmp(partType,'blackholes') or strcmp(partType,'bh') $
                                    or strcmp(partType,'bhs'))        then partType = 5
  
  ; zoom
  if (strcmp(partType,'highres_dm')  or strcmp(partType,'dm_highres')) then partType = 1
  if (strcmp(partType,'lowres_dm')   or strcmp(partType,'dm_lowres'))  then partType = 2
  if (strcmp(partType,'coarse_dm')   or strcmp(partType,'dm_coarse'))  then partType = 2
  
  if (strcmp(partType,'tracer') or strcmp(partType,'tracers')) then $
    message,'ERROR: Please specify which type of tracers!'
  
  if not isnumeric(partType) then $
    message,'ERROR: Unrecognized partType!'
  
  if (partType lt 0 or partType gt 5) then $
    message,'ERROR: partType = ' + str(partType) + ' out of bounds!'

  return, partType
end  

; plotRedshiftSpacings(): compare different runs

pro plotRedshiftSpacings

  ; config
  sP = mod_struct( sP, 'sP0', simParams(res=128,run='gadget') )
  sP = mod_struct( sP, 'sP1', simParams(res=128,run='tracer') )
  sP = mod_struct( sP, 'sP2', simParams(res=128,run='feedback') )
  ;sP = mod_struct( sP, 'sP3', simParams(res=512,run='zoom_20mpc',hind=0) )
  sP = mod_struct( sP, 'sP4', simParams(res=910,run='illustris') )

  xrange = [0.0,14.0]
  yrange = [0.5,n_tags(sP)+0.5]
  
  runNames = []
  for i=0,n_tags(sP)-1 do runNames = [runNames,sP.(i).run]

  start_PS, sP.(0).plotPath + 'redshift_spacing.eps', xs=16, ys=6
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,yticks=n_tags(sP)-1,$
      ytickv=lindgen(n_tags(sP))+1,ytickname=[runNames],yminor=1,$
      xtitle="Age of Universe [Gyr]",ytitle="Output Snapshot Spacing",xs=9,/ys,$
      pos=[0.12,0.12,0.96,0.94]
  
    ; loop over each run
    for i=0,n_tags(sP)-1 do begin
      zVals = snapNumToRedshift(sP=sP.(i),/all)
      zVals = redshiftToAgeFlat(zVals)
      
      yLoc = (i+1) + [-0.4,0.4]
      for j=0,n_elements(zVals)-1 do $
        cgPlot,[zVals[j],zVals[j]],yLoc,color=sP.(i).colors[0],thick=!p.thick-2,/overplot
    endfor
    
    ; redshift axis
    sP.(0).snap = 0 ; indicate xaxis in Gyr not snapshot number
    redshift_axis, xrange, yrange, sP=sP.(0)
   
  end_PS

  stop
  
end
  
; rhoTHisto(): make mass-weighted density-temperature 2d histogram

pro rhoTHisto

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  nbins = 40  
  sP = simParams(res=256,run='feedback_noFB',redshift=2.0)
  
  rMinMax = [-2.0,8.0]
  tMinMax = [3.0,7.0]
  
  ; load
  u_in    = loadSnapshotSUbset(sP=sP,partType='gas',field='u')
  ne_in   = loadSnapshotSUbset(sP=sP,partType='gas',field='nelec')
  temp_in = convertUtoTemp(u_in,ne_in)
  
  u_in  = !NULL
  ne_in = !NULL
  
  dens_in = loadSnapshotSubset(sP=sP,partType='gas',field='dens')
  mass    = loadSnapshotSUbset(sP=sP,partType='gas',field='mass')
  
  ; prepare data
  dens = alog10(rhoRatioToCrit(dens_in,sP=sP))
  temp = alog10(temp_in)
  
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
  start_PS,'test2.eps'
    ; color table
    loadColorTable,'helix'
      
    tvim,h2rt,scale=1,clip=[10,100],$;,/c_map
         xtitle="log ("+textoidl("\rho / \rho_{crit}")+")",ytitle="log (T [K])",$
         stitle="Total Mass ("+textoidl("M_{sun}")+")",barwidth=0.5,lcharsize=!p.charsize-1.0,$
         xrange=rMinMax,yrange=tMinMax,/rct;,nodata=0,rgb_nodata=[1.0,1.0,1.0] ;display zeros as white not black
  end_PS

end

; redshift_axis(): draw redshift axis

pro redshift_axis, xRange, yRange, ylog=ylog, sP=sP, dotted=dotted, zTicknames=zTicknames

  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  logFac = 1.0
  if keyword_set(ylog) then logFac = 10.0
  
  yTickLen = (yRange[1]-yRange[0]) / 40.0 * logFac
  yTextOff = (yRange[1]-yRange[0]) / 50.0 * logFac

  if (not keyword_set(zTicknames)) then $
    zTicknames = ['30','6','5','4','3','2','1','0.5','0.25','0']
  nZ = n_elements(zTicknames)
  
  zXPos = fltarr(nZ)
  
  ; plot horizontal line along top xaxis
  cgPlot,xRange,[yRange[1],yRange[1]]-0.0005,/overplot
  ;cgText,0.5,0.94,"Redshift",/normal  
  
  ; plot "maximum" redshift label
  if (xRange[0] le 0.0) then $
    cgText,0.0,yRange[1]+yTextOff,zTicknames[0],alignment=0.5
  
  ; skip z=30 (highest) at t=0
  for i=1,nZ-1 do begin
    if sP.snap eq -1 then begin ;x-axis in snapshot number
      zXPos[i] = redshiftToSnapNum(float(zTicknames[i]),sP=sP)
    endif else begin ;x-axis in time [gyr]
      zXPos[i] = redshiftToAgeFlat(float(zTicknames[i]))
    endelse
    
    ; plot tick mark and label if inside plotrange
    if (zXPos[i] ge xRange[0] and zXPos[i] le xRange[1]) then begin
      cgPlot,[zXPos[i],zXPos[i]],[yRange[1],yRange[1]-yTickLen],/overplot
      cgText,zXPos[i],yRange[1]+yTextOff,zTicknames[i],alignment=0.5
    endif
    
    ; plot vertical dotted line at each redshift mark if requested
    if keyword_set(dotted) then $
      cgPlot,[zXPos[i],zXPos[i]],yRange,line=dotted,/overplot,thick=!p.thick-0.5
  endfor
  
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
    cgText,0.0,yRange[1]+yTextOff,zTicknames[0],alignment=0.5
  
  ; skip t=0 (highest) at z=inf (?)
  for i=1,nZ-1 do begin

    zXPos[i] = zRedshifts[i]
    
    ; plot tick mark and label if inside plotrange
    if (zXPos[i] le xRange[0] and zXPos[i] ge xRange[1]) then begin
      cgPlot,[zXPos[i],zXPos[i]],[yRange[1],yRange[1]-yTickLen],/overplot
      cgText,zXPos[i],yRange[1]+yTextOff,zTicknames[i],alignment=0.5
    endif
    
    ; plot vertical dotted line at each redshift mark if requested
    if keyword_set(dotted) then $
      cgPlot,[zXPos[i],zXPos[i]],yRange,line=dotted,/overplot,thick=!p.thick-0.5
  endfor
  
  cgPlot,xRange,[yRange[1],yRange[1]],/overplot
  cgText,(xRange[0]-xRange[1])/2.0,yRange[1]+yTextOff*5.0,"Time [Gyr]",alignment=0.5

end

; exportParticlesAscii(): export some part of a snapshot as text

pro exportParticlesAscii

  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function simParams, loadSnapshotSubset

  partType = 'gas'
  fileName = '512fb.gas.txt'
  sP = simParams(res=512,run='tracer',redshift=0.0)
  
  pos  = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
  dens = loadSnapshotSubset(sP=sP,partType=partType,field='dens')
  u    = loadSnapshotSubset(sP=sP,partType=partType,field='u')
  
  print,n_elements(u)
  
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

; exportGroupCatAscii(): export some part of a group catalog as text (for MySQL import)
;
;   % mysql --local-infile -u <username> -p <DatabaseName>
;
;   USE subgroups;
;   LOAD DATA LOCAL INFILE 'out.txt'
;   INTO TABLE snap_135
;   COLUMNS TERMINATED BY ' ' LINES TERMINATED BY '\n'
;   (id,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,mass,r200,masstype_0,masstype_1,masstype_4,
;    masstype_5,spin_x,spin_y,spin_z,veldisp,velmax,velmaxrad,halfmassrad_0,halfmassrad_1,
;    halfmassrad_4,grnr,sfr,phot_U,phot_B,phot_V,phot_K,phot_G,phot_R,phot_I,phot_Z,fits_dirnr,
;    primary_flag,image_x,image_y)
;   set point_xy = PointFromText(CONCAT('POINT(',image_x,' ',image_y,')'));
;
;   then add SPATIAL INDEX on point_xy

pro exportGroupCatAscii
  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function simParams, loadGroupCat
  
  ; config
  sP = simParams(res=1820,run='illustris',redshift=0.0)
  
  fileName   = 'out.txt'
  minNumPart = 100
  zWidth     = 0 ; 15000.0 ; cMpc/h centered on fof0 group position (or 0 to disable)
  
  ; load group catalog and select subgroups for export
  gc = loadGroupCat(sP=sP,/skipIDs,/skipOffsets)

  zBounds = [-0.1, sP.boxSize + 0.1]
  if zWidth gt 0 then begin
    zBounds[0] = gc.groupPos[2,0] - zWidth*0.5
    zBounds[1] = gc.groupPos[2,0] + zWidth*0.5
  endif
  
  ; make image (fof0 centered) coordinates
  image_x = reform( gc.subgroupPos[0,*] )
  image_y = reform( gc.subgroupPos[1,*] )
  image_x = image_x - gc.subgroupPos[0,0] + sP.boxSize*0.5
  image_y = image_y - gc.subgroupPos[1,0] + sP.boxSize*0.5

  correctPeriodicPosVecs, image_x, sP=sP
  correctPeriodicPosVecs, image_y, sP=sP

  ; make export selection
  w = where(gc.subgroupLen ge minNumPart and $
            gc.subgroupPos[2,*] ge zBounds[0] and $
            gc.subgroupPos[2,*] le zBounds[1],count)
  print,'Exporting ['+str(count)+'] SUBGROUPS.'

  ; load paul's additional array
  pt_directory = '/n/ghernquist/ptorrey/Illustris/SubfindParsedSnapshots/Illustris-1/snapshot_135/directory.txt'
  readcol,pt_directory,pt_subnr,pt_dirnr,format='L,L'
  
  pt_full_dirnr_array = intarr(gc.nSubgroupsTot) - 1   ; set all to -1
  pt_full_dirnr_array[pt_subnr] = pt_dirnr             ; set directory elements

  ; fill output buffer
  outBuf = strarr(count)
  
  for i=0,count-1 do begin
    if i mod round(count*0.1) eq 0 then print,i,count
    
    gInd = gc.SubgroupGrNr[w[i]]
    pri_flag = gc.groupFirstSub[gInd] eq w[i]
    pt_full_dirnr = pt_full_dirnr_array[w[i]]
    
    outBuf[i] = str(string(w[i],                     format='(I8)'))    + " " + $ ; ID
                str(string(gc.SubgroupPos[0,w[i]],   format='(f12.6)')) + " " + $ ; pos_x
                str(string(gc.SubgroupPos[1,w[i]],   format='(f12.6)')) + " " + $ ; pos_y
                str(string(gc.SubgroupPos[2,w[i]],   format='(f12.6)')) + " " + $ ; pos_z
                str(string(gc.SubgroupVel[0,w[i]],   format='(f12.6)')) + " " + $ ; vel_x
                str(string(gc.SubgroupVel[1,w[i]],   format='(f12.6)')) + " " + $ ; vel_y
                str(string(gc.SubgroupVel[2,w[i]],   format='(f12.6)')) + " " + $ ; vel_z
                str(string(gc.SubgroupMass[w[i]],    format='(f12.6)')) + " " + $ ; mass
                str(string(gc.group_r_crit200[gInd], format='(f8.2)'))  + " " + $ ; r200
                str(string(gc.SubgroupMassType[0,w[i]],            format='(f12.6)')) + " " + $
                str(string(gc.SubgroupMassType[1,w[i]],            format='(f12.6)')) + " " + $
                str(string(gc.SubgroupMassType[4,w[i]],            format='(f12.6)')) + " " + $
                str(string(gc.SubgroupMassType[5,w[i]],            format='(f12.6)')) + " " + $
		    str(string(gc.SubgroupSpin[0,w[i]],                format='(f12.6)')) + " " + $
                str(string(gc.SubgroupSpin[1,w[i]],                format='(f12.6)')) + " " + $
                str(string(gc.SubgroupSpin[2,w[i]],                format='(f12.6)')) + " " + $
		    str(string(gc.SubgroupVelDisp[w[i]],               format='(f12.6)')) + " " + $
		    str(string(gc.SubgroupVmax[w[i]],                  format='(f12.6)')) + " " + $
		    str(string(gc.SubgroupVmaxRad[w[i]],               format='(f12.6)')) + " " + $
		    str(string(gc.SubgroupHalfMassRadType[0,w[i]], 	 format='(f12.6)')) + " " + $
                str(string(gc.SubgroupHalfMassRadType[1,w[i]], 	 format='(f12.6)')) + " " + $
                str(string(gc.SubgroupHalfMassRadType[4,w[i]], 	 format='(f12.6)')) + " " + $
		    str(string(gc.SubgroupGrnr[w[i]],                  format='(I8)'))    + " " + $
		    str(string(gc.SubgroupSFR[w[i]],                   format='(f12.6)')) + " " + $
		    str(string(gc.SubgroupStellarPhotometrics[0,w[i]], format='(f12.6)')) + " " + $
                str(string(gc.SubgroupStellarPhotometrics[1,w[i]], format='(f12.6)')) + " " + $
                str(string(gc.SubgroupStellarPhotometrics[2,w[i]], format='(f12.6)')) + " " + $
                str(string(gc.SubgroupStellarPhotometrics[3,w[i]], format='(f12.6)')) + " " + $
                str(string(gc.SubgroupStellarPhotometrics[4,w[i]], format='(f12.6)')) + " " + $
                str(string(gc.SubgroupStellarPhotometrics[5,w[i]], format='(f12.6)')) + " " + $
                str(string(gc.SubgroupStellarPhotometrics[6,w[i]], format='(f12.6)')) + " " + $
                str(string(gc.SubgroupStellarPhotometrics[7,w[i]], format='(f12.6)')) + " " + $
		    str(string(pt_full_dirnr,	                         format='(I4)'))    + " " + $
                str(string(pri_flag,                               format='(I1)'))    + " " + $ ; pri_flag
                str(string(image_x[w[i]],                          format='(f12.6)')) + " " + $ ; image_x
                str(string(image_y[w[i]],                          format='(f12.6)'))           ; image_y
                
  endfor
  
  ; open ascii, write and close
  openW,lun,fileName,/get_lun
  
    printf,lun,outBuf
  
  close,lun
  free_lun,lun

  stop

end
