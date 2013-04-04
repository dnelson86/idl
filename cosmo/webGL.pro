; webGL.pro
; cosmological boxes - exporters, etc for webgl demos
; dnelson apr.2013

; vorMeshExport(): export the VORONOI_MESHOUTPUT data files

pro vorMeshExport

  ; config
  fileNameOut = "vorMesh_Arepo3b.dat"
  fileBaseIn = "voronoi_mesh_0_"
  bboxSize = 1000.0 ; don't change
  
  ; Arepo2b/3b/other test boxes
  x_minmax = [0,1]
  y_minmax = [0,1]
  z_minmax = [0,1]
  
  ; DIEGO
  ;x_minmax = [990,1020]
  ;y_minmax = [990,1020]
  ;z_minmax = [990,1020]
  
  ; 512z2h304
  ;x_minmax = [16465,17120]
  ;y_minmax = [00545,01200]
  ;z_minmax = [14580,15230]
  
  fname = { coords  : fileBaseIn + "coordinates.dat" ,$
            inds    : fileBaseIn + "indices.dat"     ,$
	      normals : fileBaseIn + "normals.dat"      }
	    
  ; load vertices (coordinates)
  openr,lun,fname.coords,/get_lun
    nVertsTot = 0L
    readU, lun, nVertsTot ; header
    pts = fltarr(3,nVertsTot) ; replicate
    readU, lun, pts ; fill
  close,lun
  free_lun,lun
  
  print,'Loaded ['+str(nVertsTot)+'] total vertices (min='+str(min(pts))+' max='+str(max(pts))+').'
  
  ; arrays for face info
  nFacesRead   = 0L
  nEntriesRead = 0L
  nIndsRead    = 0L
  
  maxNumFaces = 30000000L
  faceIDs     = lonarr(maxNumFaces)
  faceLengths = intarr(maxNumFaces) - 1
  faceInds    = lonarr(50*maxNumFaces)
  
  ; load inds (faces)
  openr,lun,fname.inds,/get_lun
    nEntriesTot = 0L
    readU, lun, nEntriesTot ; header
    
    ; loop over faces one by one
    facesLeftToRead = 1L
    while nEntriesRead lt nEntriesTot do begin
      ; read the length (number of verts) for this face
      face_id = 0L
      face_len = 0L
      readU, lun, face_id
      readU, lun, face_len
      if face_len lt 3 then message,'Error: Face has less than 3 vertices.'
      
      ; read the indices (into the vertex array) and store
      face_inds = lonarr(face_len)
      readU, lun, face_inds
      
      ; read the 3 normals (in other file)
      ;face_normals = fltarr(3)
      ;readU, lun, face_normals
      
      ; store face data
      faceIDs[nFacesRead] = face_id
      faceLengths[nFacesRead] = face_len
      faceInds[nIndsRead : nIndsRead + face_len - 1] = face_inds
      
      nFacesRead += 1
      nEntriesRead += face_len + 2
      nIndsRead += face_len
      
      if nFacesRead mod 100000L eq 0 then print,face_id,face_len,nFacesRead,nEntriesRead
      if nFacesRead ge maxNumFaces then message,'Error: Too many faces.'
    endwhile
  close,lun
  free_lun,lun
  
  ; verify counts
  if nEntriesTot ne nEntriesRead then message,'Error: Failed to read correct number of total face entries.'
  print,'Loaded ['+str(nFacesRead)+'] total faces, ['+str(nEntriesRead)+'] total entries.'
  
  ; take subset of valid faces
  faceIDs     = faceIDs[0:nFacesRead-1]
  faceLengths = faceLengths[0:nFacesRead-1]
  faceInds    = faceInds[0:nIndsRead-1]
  
  ; make offset table
  faceOffsets = [0,total(faceLengths,/cum,/int)]
  
  ; load normals (skip)
  
  ; compute face centroids
  offset = 0L
  faceCentroids = fltarr(3,nFacesRead)
  
  for i=0L,nFacesRead-1 do begin
    vert_inds = faceInds[offset : offset + faceLengths[i] - 1]
    
    faceCentroids[0,i] = mean(pts[0,vert_inds])
    faceCentroids[1,i] = mean(pts[1,vert_inds])
    faceCentroids[2,i] = mean(pts[2,vert_inds])
    
    offset += faceLengths[i]
  endfor
            
  ; make a spatial subset?
  if n_elements(x_minmax) gt 0 then begin
    w = where(faceCentroids[0,*] ge x_minmax[0] and faceCentroids[0,*] le x_minmax[1] and $
              faceCentroids[1,*] ge y_minmax[0] and faceCentroids[1,*] le y_minmax[1] and $
              faceCentroids[2,*] ge z_minmax[0] and faceCentroids[2,*] le z_minmax[1],count)  

    print,'Found ['+str(count)+'] faces in spatial subset.'
            
    faceCentroids = faceCentroids[*,w]
    faceLengths   = faceLengths[w]
    faceOffsets   = faceOffsets[w]              
  endif
  
  ; construct dataout1 array
  dataout1 = fix(faceLengths)
  
  if min(faceLengths) lt 0 or max(faceLengths) gt 20 then message,'Error: Strange face lengths.'
  if min(faceInds) lt 0 then print,'WARNING: A face references the bounding tetra.'
  
  ; compress points into [-1000,1000] range (maintain aspect ratio if not cube)
  xyzr = float([x_minmax[1]-x_minmax[0], y_minmax[1]-y_minmax[0], z_minmax[1]-z_minmax[0]])
  largest = max(xyzr,largest_ind)
  aspect_mod = xyzr / xyzr[largest_ind]
  
  xyzmin = [x_minmax[0],y_minmax[0],z_minmax[0]]
  
  for j=0,2 do begin
    pts[j,*] = (pts[j,*] - xyzmin[j]) / xyzr[j] * aspect_mod[j] * bboxSize*2 - bboxSize
    faceCentroids[j,*] = (faceCentroids[j,*] - xyzmin[j]) / xyzr[j] * aspect_mod[j] * bboxSize*2 - bboxSize
  endfor
  
  ; construct dataout2 array
  dataout2 = fltarr( (total(faceLengths,/int)) * 3 + n_elements(faceCentroids) )
  offset = 0L
  
  for i=0,count-1 do begin
    ; for this face, first send centroid x,y,z
    dataout2[offset+0] = faceCentroids[0,i]
    dataout2[offset+1] = faceCentroids[1,i]
    dataout2[offset+2] = faceCentroids[2,i]
    offset += 3
    
    ; get indices of vertices
    vert_inds = faceInds[faceOffsets[i] : faceOffsets[i] + faceLengths[i] - 1]
    
    ; then send each vertex x,y,z
    for j=0,faceLengths[i]-1 do begin
      dataout2[offset+0] = pts[0,vert_inds[j]]
      dataout2[offset+1] = pts[1,vert_inds[j]]
      dataout2[offset+2] = pts[2,vert_inds[j]]
      offset += 3
    endfor
  endfor
  
  if n_elements(dataout2) ne offset then message,'Error: Failed to fill dataout2 completely.'

  ; header information
  simType  = 2

  ;dataout2 -= 10000.0 ; [-10000,10000]
  ;dataout2 /= 10.0 ; [-1000,1000]
  ;print,'COMPRESSING COSMO'
  
  ; write binary file
  outStruct = { fileType  : fix(1)              ,$ ; 0=debug, 1=no normals
                simType   : fix(simType)        ,$ ; 1=gadget, 2=arepo
                bboxSize  : float(bboxSize)     ,$
                nFaces    : long(count)         ,$
                hVertices : long(nVertsTot)     ,$ ; not actually used for anything as vertices are send explicitly
                data1     : dataout1            ,$
                data2     : dataout2             $
              }
		
  ; write file
  openw,lun,fileNameOut,/get_lun
  writeu,lun,outStruct
  close,lun
  free_lun,lun
  
  stop

end

; trTrajExport(): export velocity tracer tracks with time relative to tracked halo position

pro trTrajExport

  ; config
  sP = simParams(res=128,run='tracer',redshift=1.0)
  minSnap = redshiftToSnapNum(5.0,sP=sP)
  minLogMass = 11.5
  
  ; select halo
  mt = mergerTreeSubset(sP=sP)
  gc = loadGroupCat(sP=sP,/skipIDs)
  galcat = galaxyCat(sP=sP)
  
  w = where(codeMassToLogMsun(gc.subgroupMass[mt.galcatIDList]) gt minLogMass and $
            mt.hMinSnap lt minSnap and mt.hMinSnap ne -1,count)
	    
  mtID = w[0]
  gcID = mt.galcatIDList[mtID] ;5
  
  ; find all velocity tracer children in galaxy at this redshift
  gasInds = galcatINDList(galcat=galcat,gcIDList=[gcID])
  trIDs = cosmoTracerVelChildren(sP=sP,gasIDs=galcat.galaxyIDs[gasInds.gal],/getIDs)
  
  ; override hMinSnap to take less, and/or number of tracers
  ;mt.hMinSnap[mtID] = 200
  ;trIDs = trIDs[floor(randomu(424242L,20)*n_elements(trIDs))]
  
  ; header information
  nTracers = n_elements(trIDs)
  nTimes   = sP.snap - mt.hMinSnap[mtID] + 1
  fileName = 'trTraj.' + sP.plotPrefix + '.' + str(sP.res) + '.' + str(sP.snap) + '-' + $
             str(mt.hMinSnap[mtID]) + '.h' + str(gcID) + '.ntr' + str(nTracers) + '.dat'
  simType  = 2 ; arepo
  
  ; allocate arrays
  times    = fltarr(nTimes)
  virRads  = mt.hVirRad[0:(sP.snap-mt.hMinSnap[mtID]),mtID]
  virTemps = mt.hVirTemp[0:(sP.snap-mt.hMinSnap[mtID]),mtID]
  
  pos  = fltarr(nTracers,3,nTimes)
  temp = fltarr(nTracers,nTimes)
  
  ; loop back as far as halo is tracked
  snapRange = [sP.snap,mt.hMinSnap[mtID]]
  
  for m=snapRange[0],snapRange[1],-1 do begin
    ; load snapshot data
    sP.snap = m & print,m
    curInd = snapRange[0] - m
    
    h = loadSnapshotHeader(sP=sP)
    times[curInd] = h.time
    
    ; load tracer ids and match
    loc_ids = loadSnapshotSubset(sP=sP,partType='tracerVel',field='ids')
    calcMatch,loc_ids,trIDs,loc_ids_ind,trIDs_ind,count=countMatch
    loc_ids_ind = loc_ids_ind[calcSort(trIDs_ind)]
    
    if countMatch ne nTracers then message,'Error: Failed to find all tracers'
    
    ; load positions and temperatures
    loc_pos = loadSnapshotSubset(sP=sP,partType='tracerVel',field='pos')
    for i=0,2 do pos[*,i,curInd] = loc_pos[i,loc_ids_ind]
    for i=0,2 do pos[*,i,curInd] -= reform(replicate(mt.hPosSm[curInd,i,mtID],nTracers))
    
    loc_temp = loadSnapshotSubset(sP=sP,partType='tracerVel',field='tracer_maxtemp')
    temp[*,curInd] = codeTempToLogK(loc_temp[loc_ids_ind]) ; vel tracer output in unit system
    
  endfor
  
  ; prepare output data
  fltperpt = 3+1
  dataout = fltarr(nTracers*nTimes*fltperpt)
  
  count = 0L
  for i=0L,nTracers-1 do begin
    xpos = reform(pos[i,0,*])
    ypos = reform(pos[i,1,*])
    zpos = reform(pos[i,2,*])
    tt   = reform(temp[i,*])
    chunk = transpose([[xpos],[ypos],[zpos],[tt]]) ;x0y0z0t0x1y1...
    chunk = reform(chunk,n_elements(chunk)) ; 1d
    dataout[count:count+n_elements(chunk)-1] = chunk
    count += n_elements(chunk)
  endfor
  
  if count ne nTracers*nTimes*fltperpt then message,'error'
    
  ; write binary file
  outStruct = { fileType : fix(1)              ,$ ; 0=debug, 1=pos,temp
                simType  : fix(simType)        ,$ ; 1=gadget, 2=arepo
                nTracers : long(nTracers)      ,$
                nTimes   : long(nTimes)        ,$
                hInd     : long(gcID)          ,$ ; subgroup id
                times    : float(times)        ,$ ; scale factor
                virRads  : float(virRads)      ,$ ; ckpc
                virTemps : float(virTemps)     ,$ ; log K
                data     : dataout $
              }
		
  ; write file
  openw,lun,fileName,/get_lun
  writeu,lun,outStruct
  close,lun  
  
  w = where(pos eq 0.0,count1)
  w = where(temp eq 0.0,count2)
  
  print,count1,count2
  
  stop

end

; trTrajTest(): debug export for tracer trajectory demo

pro trTrajTest

  ; config
  seed = 424242L
  fileName = 'trTraj.test10.dat'
  simType  = 2
  nTracers = 10
  nTimes   = 100
  gcID     = 333
  
  ; time arrays
  times    = findgen(nTimes)/(nTimes-1)*0.7+0.1
  virRads  = findgen(nTimes)*2+100.0
  virTemps = findgen(nTimes)/50+5.0

  ; data arrays
  pos = fltarr(nTracers,3,nTimes)
  temp = fltarr(nTracers,nTimes)
  
  ; fill positions
  ms = randomn(424242L,3,nTracers)
  w = where(ms ge 0.0,comp=wc)
  ms[w] = 1.0 & ms[wc] = -1.0
  
  for i=0,nTracers-1 do begin
    pos[i,0,*] = ms[0,i] * findgen(nTimes)/(nTimes-1) * 500.0 + 10.0
    pos[i,1,*] = ms[1,i] * findgen(nTimes)/(nTimes-1) * 500.0 + 20.0
    pos[i,2,*] = ms[2,i] * findgen(nTimes)/(nTimes-1) * 500.0 + 30.0
    
    ; random wiggle
    pos[i,0,*] += randomu(seed,nTimes)*50
    pos[i,1,*] += randomu(seed,nTimes)*50
    pos[i,2,*] += randomu(seed,nTimes)*50
  endfor
  
  ; fill temps
  for i=0,nTracers-1 do begin
    temp[i,*] = findgen(nTimes)/(nTimes-1) * 3.0 + 4.0 ; 4-7
  endfor
  
  ; prepare output data
  fltperpt = 3+1
  dataout = fltarr(nTracers*nTimes*fltperpt)
  
  count = 0
  for i=0,nTracers-1 do begin
    xpos = reform(pos[i,0,*])
    ypos = reform(pos[i,1,*])
    zpos = reform(pos[i,2,*])
    tt   = reform(temp[i,*])
    chunk = transpose([[xpos],[ypos],[zpos],[tt]]) ;x0y0z0t0x1y1...
    chunk = reform(chunk,n_elements(chunk)) ; 1d
    dataout[count:count+n_elements(chunk)-1] = chunk
    count += n_elements(chunk)
  endfor
  
  if count ne nTracers*nTimes*fltperpt then message,'error'
    
  ; write binary file
  outStruct = { fileType : fix(0)              ,$ ; 0=debug
                simType  : fix(simType)        ,$ ; 1=gadget, 2=arepo
                nTracers : long(nTracers)      ,$
                nTimes   : long(nTimes)        ,$
                hInd     : long(gcID)          ,$ ; subgroup id
                times    : float(times)        ,$ ; scale factor
                virRads  : float(virRads)      ,$ ; ckpc
                virTemps : float(virTemps)     ,$ ; log K
                data     : dataout $
              }
		
  ; write file
  openw,lun,fileName,/get_lun
  writeu,lun,outStruct
  close,lun

end

; webglCutouts(): make halo centered cutouts for the webGL app

pro webglCutouts;, sP=sP, gcIDs=gcIDs

  compile_opt idl2, hidden, strictarr, strictarrsubs

  sP = simParams(res=512,run='gadget',redshift=2.0)
  gcIDs = [2342] ;g2342 a2132 (z2.304)
  ;if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  sizeFac = 5.0 ; times rvir for the bounding box of each cutout
  
  ; target list
  gc    = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sP) 

  ; load gas properties
  u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
  nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
  pos   = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  dens  = loadSnapshotSubset(sP=sP,partType='gas',field='dens')
  hsml  = loadSnapshotSubset(sP=sP,partType='gas',field='hsml')
  vel   = loadSnapshotSubset(sP=sP,partType='gas',field='vel')

  ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
  print,'shuffling...'
  iseed = 424242L
  sort_inds = sort(randomu(iseed,n_elements(u)))
  
  u     = u[sort_inds]
  nelec = nelec[sort_inds]
  dens  = dens[sort_inds]
  hsml  = hsml[sort_inds]
  pos   = pos[*,sort_inds]
  vel   = vel[*,sort_inds]
  
  sort_inds = !NULL
  
  if sP.trMCPerCell eq 0 then simType = 1 ;gadget
  if sP.trMCPerCell ne 0 then simType = 2 ;arepo
  
  print,'saving...'
  ; loop over all requested halos
  foreach gcID, gcIDs do begin
  
    ; get subhalo position and size of imaging box
    boxCen     = sgcen[*,gcID]
    boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcID]] / 10.0) * 10.0
    boxSizeImg = [boxSize,boxSize,boxSize] ; cube
  
    ; make conservative cutout greater than boxsize accounting for periodic (do cube not sphere)
    xDist = pos[0,*] - boxCen[0]
    yDist = pos[1,*] - boxCen[1]
    zDist = pos[2,*] - boxCen[2]
    
    correctPeriodicDistVecs, xDist, sP=sP
    correctPeriodicDistVecs, yDist, sP=sP
    correctPeriodicDistVecs, zDist, sP=sP
    
    rvir = gc.group_r_crit200[gc.subgroupGrNr[gcID]]
    virtemp = codeMassToVirTemp(gc.subgroupmass[gcID],sP=sP)
    rad = reform(sqrt(xDist*xDist + yDist*yDist + zDist*zDist)) / rvir[0]
  
    ; local (cube) cutout
    wCut = where(abs(xDist) le 0.5*boxSize and abs(yDist) le 0.5*boxSize and $
                 abs(zDist) le 0.5*boxSize,nCutout)
    print,gcID,nCutout
           
    ; field subsets   
    loc_u     = u[wCut]
    loc_nelec = nelec[wCut]
    loc_dens  = alog10(dens[wCut] * 1e10) ; log msun/ckpc^3 (i.e. comoving density)
    loc_hsml  = hsml[wCut]
    
    ; derived
    loc_temp = alog10(convertUtoTemp(loc_u,loc_nelec)) ; log T (K)
    loc_pres = alog10(calcPressureCGS(loc_u,loc_dens,sP=sP)) ; log P/k (K/cm^3)
    loc_entr = alog10(calcEntropyCGS(loc_u,loc_dens,sP=sP)) ; log S (K cm^2)
    
    ; relative positions
    loc_pos  = fltarr(3,nCutout)
    loc_pos[0,*] = xDist[wCut] ; delta
    loc_pos[1,*] = yDist[wCut]
    loc_pos[2,*] = zDist[wCut]
    
    xDist = !NULL
    yDist = !NULL
    zDist = !NULL
    
    ; velocities: calculate radial velocity relative to bulk halo motion
    loc_vel = vel[*,wCut]
    
    gVel = gc.subgroupVel[*,gcID]
    loc_vel[0,*] = reform(loc_vel[0,*] - gVel[0])
    loc_vel[1,*] = reform(loc_vel[1,*] - gVel[1])
    loc_vel[2,*] = reform(loc_vel[2,*] - gVel[2])
    
    ; make normalized position vector wrt halo center = vec(r) / ||r|| where r from particle to center
    ; means that radvel<0 is inflow and radvel>0 is outflow
    rnorm0 = reform(loc_pos[0,*] - boxCen[0])
    rnorm1 = reform(loc_pos[1,*] - boxCen[1])
    rnorm2 = reform(loc_pos[2,*] - boxCen[2])
    
    correctPeriodicDistVecs, rnorm0, sP=sP
    correctPeriodicDistVecs, rnorm1, sP=sP
    correctPeriodicDistVecs, rnorm2, sP=sP
    
    ; denominator and do divide
    rnorm = sqrt(rnorm0*rnorm0 + rnorm1*rnorm1 + rnorm2*rnorm2)

    rnorm0 /= rnorm
    rnorm1 /= rnorm
    rnorm2 /= rnorm
    
    ; dot(vel,rnorm) gives the magnitude of the projection of vel onto vec(r)
    loc_radvel = loc_vel[0,*]*rnorm0 + loc_vel[1,*]*rnorm1 + loc_vel[2,*]*rnorm2 ; 1xN

    ; velocities: create endpoint for each position point for the velocity vector line
    ;loc_pos2 = fltarr(3,nCutout)
    ;loc_pos2[0,*] = loc_pos[0,*] + loc_vel[0,*]*velVecFac
    ;loc_pos2[1,*] = loc_pos[1,*] + loc_vel[1,*]*velVecFac
    ;loc_pos2[2,*] = loc_pos[2,*] + loc_vel[2,*]*velVecFac
  
    loc_hsml = reform(loc_hsml,[1,nCutout]) ; make 1xN vector
    loc_temp = reform(loc_temp,[1,nCutout])
    loc_dens = reform(loc_dens,[1,nCutout])
    loc_pres = reform(loc_pres,[1,nCutout])
    loc_entr = reform(loc_entr,[1,nCutout])

    ; prepare output data
    dataout = [loc_pos,loc_hsml,loc_temp,loc_dens,loc_pres,loc_entr,loc_radvel]
    
    ; write binary file
    fileName = 'cutout.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+'.dat'
    
    outStruct = { fileType : fix(1)              ,$ ; 0=gas[x,y,z,hsml,temp,dens]
                                                  $ ; 1=gas[x,y,z,hsml,temp,dens,pres,entr,radvel]
                  simType  : fix(simType)        ,$ ; 1=gadget, 2=arepo
                  nPts     : long(nCutout)       ,$
                  redshift : float(sP.redshift)  ,$
                  hInd     : long(gcID)          ,$ ; subgroup id
                  hVirRad  : float(rvir)         ,$ ; ckpc
                  hVirTemp : float(virtemp)      ,$ ; log K
                  sizeFac  : float(sizeFac)      ,$ ; bounding box size
                  data     : reform(dataout,n_elements(dataout)) $
                }
                
    ; write file
    openw,lun,fileName,/get_lun
    writeu,lun,outStruct
    close,lun
  endforeach
end
