; spiralsLoad.pro
; snapshot/data loading
; dnelson may.2011

; loadMGParameterFile(): load simName.hdf5.parameters file and fill h2 struct

function loadMGSimParams, filePath, h2=h2

  ; check file exists
  if not (file_test(filePath)) then return,0
  
  headerLines = 0
  lines = loadCSV(headerLines,filePath,'')
  
  ; verify format
  if (lines[n_elements(lines)-1] ne 'done') then return,0

  ; parse lines into structure
  nameVals = replicate({name:'',val:''},n_elements(lines)-1)
  
  for i=0,n_elements(lines)-2 do begin
    line = strsplit(strcompress(lines[i])," ",/extract)
    nameVals[i].name = line[0]
    nameVals[i].val  = line[1]
  endfor

  ; create new sP struct
  sP = mrd_struct(nameVals.name,nameVals.val,1)
  
  ; fill old h2 struct
  if keyword_set(h2) then begin
    h2.c    = sP.c
    h2.v200 = sP.v200
    h2.f_d  = sP.md
    h2.f_R  = sP.sig_r
    h2.m200 = sP.m200
    h2.a    = sP.rh
    h2.h    = sP.rd
    ;lc
  endif
  
  return,sP
end

; loadSimParams():

function loadSimParams, simName

  ; header2 (sim parameters)
  h2 = { c:      0.0,     $
         v200:   0.0,     $
         f_d:    0.0,     $
         f_R:    0.0,     $ ;RDF
         m200:   0.0,     $ ;=M_DM, MDG output
         a:      0.0,     $ ;halo, MDG output
         h:      0.0,     $ ;disk, MDG output
         lc:     0.0      $ ;lamdbda crit at 2h
       }
   
  ; load from MakeGalaxy output parameter file if available
  mgParamFile = "/n/home07/dnelson/spirals/lambdacrit/sims/ICs/" + $
                simName + "/" + simName + ".hdf5.parameters"
  if file_test(mgParamFile) then begin
    sP = loadMGSimParams(mgParamFile,h2=h2)
    return,h2
  endif
       
  if (strpos(simName,'LC_10m_disk_2') ne -1) then begin
    h2.c    = 9.8
    h2.v200 = 157.0
    h2.f_d  = 0.019
    h2.f_R  = 1.15
    h2.m200 = 89.9826
    h2.a    = 27.4893
    h2.h    = 3.10
    h2.lc   = 3.50
  endif
  if (strpos(simName,'LC_10m_disk_4') ne -1) then begin
    h2.c    = 8.2
    h2.v200 = 202.0
    h2.f_d  = 0.021
    h2.f_R  = 1.35
    h2.m200 = 191.652
    h2.a    = 40.1453
    h2.h    = 3.93
    h2.lc   = 5.50
  endif
  if (strpos(simName,'LC_10m_disk_6') ne -1) then begin
    h2.c    = 6.0
    h2.v200 = 210.0
    h2.f_d  = 0.023
    h2.f_R  = 1.50
    h2.m200 = 215.337
    h2.a    = 51.6477
    h2.h    = 4.29
    h2.lc   = 7.50
  endif
  if (strpos(simName,'LC_10m_disk_8') ne -1 or strpos(simName,'LC_100m_disk_8') ne -1 or $
      strpos(simName,'LC_50m_disk_8') ne -1 or strpos(simName,'LC_5m_disk_8')   ne -1 or $
      strpos(simName,'LC_1m_disk_8')  ne -1 or $
      strpos(simName,'LCLH_1e6_disk_8') ne -1 or strpos(simName,'LCLH_1e5_disk_8') ne -1 or $
      strpos(simName,'LCMGE_1m_disk_8') ne -1 or strpos(simName,'LCMGM_1m_disk_8') ne -1) then begin
    h2.c    = 5.6
    h2.v200 = 243.0
    h2.f_d  = 0.025
    h2.f_R  = 1.60
    h2.m200 = 333.640
    h2.a    = 62.5395
    h2.h    = 4.62
    h2.lc   = 9.50
  endif
  
  if (strpos(simName,'bar_1m_a') ne -1 or strpos(simName,'bar_1m_b') ne -1) then begin
    h2.c    = 6.0
    h2.v200 = 240.0
    h2.f_d  = 0.04
    h2.f_R  = 1.0 ;NOTE, LAMBDA CHANGED TO 0.05 FROM 0.033
    h2.m200 = 321.435
    h2.a    = 59.0259
    h2.h    = 3.60
    h2.lc   = 0.0
  endif  
  
  if (strpos(simName,'CT_10mil_') ne -1) then begin ;nopert,sigma0,sigma0.1,sigma0.33,sigma1
    h2.c    = 9.0
    h2.v200 = 160.0
    h2.f_d  = 0.02
    h2.f_R  = 1.0
    h2.m200 = 95.2401
    h2.a    = 29.7754
    h2.h    = 3.13
    h2.lc   = 0.0
  endif  
  
  if (strpos(simName,'m101.elena') ne -1) then begin ;has bulge
    h2.c    = 9.0
    h2.v200 = 160.0
    h2.f_d  = 0.04
    h2.f_R  = 1.0
    h2.m200 = 95.2401
    h2.a    = 29.7754
    h2.h    = 2.76
    h2.lc   = 0.0
  endif  

  if (h2.c eq 0) then $
    print,' LOADSIMPARAMS FAILED WARNING'

  return,h2
end

; loadSnapshotHDF5(): load HDF5 (format=3) snapshot
;
; inputs:  filePath, s_pos, s_vel, s_id, c_pos, c_vel, c_id
; outputs: complete struct (not called directly, only from loadSnapshot)

function loadSnapshotHDF5, filePath, s_pos, s_vel, s_id, c_pos, c_vel, c_id

  ; recursively load whole hdf5 file into a structure
  s = h5_parse(filePath,/read_data)

  ; header structure
  h = { headerHDF                                                      ,$
        nPartThisFile       : s.header.numPart_ThisFile._DATA          ,$
        nPartTot            : s.header.numPart_Total._DATA             ,$
        nPartTotHighword    : s.header.numPart_Total_Highword._DATA    ,$
        massTable           : s.header.massTable._DATA                 ,$
        time                : s.header.time._DATA                      ,$
        redshift            : s.header.redshift._DATA                  ,$
        boxSize             : s.header.boxSize._DATA                   ,$
        numFilesPerSnapshot : s.header.numFilesPerSnapshot._DATA       ,$
        Omega0              : s.header.Omega0._DATA                    ,$
        OmegaLambda         : s.header.OmegaLambda._DATA               ,$
        hubbleParam         : s.header.hubbleParam._DATA               ,$
        flagSFR             : s.header.flag_SFR._DATA                  ,$
        flagCooling         : s.header.flag_Cooling._DATA              ,$
        flagStellarAge      : s.header.flag_StellarAge._DATA           ,$
        flagMetals          : s.header.flag_Metals._DATA               ,$
        flagFeedback        : s.header.flag_Feedback._DATA             ,$
        flagDoublePrecision : s.header.flag_DoublePrecision._DATA      $
        ;compositionVecLen   : s.header.Composition_Vector_Length._DATA  $ ;not in Gadget3 output 
        ;flagICInfo          : s.header.flag_IC_Info._DATA               $ ;?
      }
  
  ; particle type 0 (gas)
  ; particle type 1 (halo)

  ; particle type 2 (disk stars)
  if (total(tag_names(s) eq 'PARTTYPE2') ge 1.0) then begin
    s_pos  = s.parttype2.coordinates._DATA
    s_vel  = s.parttype2.velocities._DATA
    s_id   = s.parttype2.particleIDs._DATA
  endif
  
  ; particle type 3 (MCs if no bulge/bulge)
  if (total(tag_names(s) eq 'PARTTYPE3') ge 1.0) then begin
    if (h.nPartTot[4] eq 0) then begin
      ; no bulge
      c_pos  = s.parttype3.coordinates._DATA
      c_vel  = s.parttype3.velocities._DATA
      c_id   = s.parttype3.particleIDs._DATA
    endif
  endif
  
  ; particle type 4 ("stars"/MCs if bulge)
  ; particle type 5 (boundary)

  return, h
end

; loadSnapshot():
;
; inputs:  fBase, i
; outputs: h=header returned
;          pos,vel,id,mass,u,rho,hsml overwritten with data

function loadSnapshot, fBase, i, $
                       s_pos, s_vel, s_id, $
                       c_pos, c_vel, c_id

  ; set filename
  if (str(i) eq 'none') then begin
    f = fBase
  endif else begin  
    ext = string(i,format='(i3.3)')
    f = fBase + ext
  endelse
  
  ; check existance
  if not file_test(f) then begin
  
    ; check for HDF5 format
    if file_test(f+".hdf5") then begin
      if h5f_is_hdf5(f+".hdf5") then begin
        h = loadSnapshotHDF5(f+".hdf5",s_pos,s_vel,s_id,c_pos,c_vel,c_id)
        return,h
      endif
    endif
    
    print, 'ERROR: snapshot ' + fBase + str(i) + ' does not exist!'
    return,0
  endif
  
  ; header structure
  bytesLeft     = 136
  h  = { header,                          $
        nPartTot:     lonarr(6),          $
        massTable:    dblarr(6),          $
        time:         0.0D,               $
        redshift:     0.0D,               $
        flagSFR:      0L,                 $
        flagFeedback: 0L,                 $
        npartall:     lonarr(6),          $
        la:           intarr(bytesLeft/2) $
      }
  
  openr,1,f,/f77_unformatted
    ;read header
    readu,1,h

    nStars  = h.nPartTot[2] ;stars, parttype=2
    
    ; no bulge
    if (h.nPartTot[4] eq 0) then begin
      nBulge  = 0
      nClouds = h.nPartTot[3] ;clouds (usually bulge), parttype=3
      nTot    = nStars + nClouds
    endif
    
    ; bulge
    if (h.nPartTot[4] ne 0) then begin
      nBulge  = h.nPartTot[3]
      nClouds = h.nPartTot[4]
      nTot    = nStars + nBulge + nClouds
    endif
    
    ;create arrays
    pos = fltarr(3,nTot)
    vel = fltarr(3,nTot)
    id  = lonarr(nTot)
    
    s_pos = fltarr(3,nStars+nBulge)
    s_vel = fltarr(3,nStars+nBulge)
    s_id  = lonarr(nStars+nBulge)
    
    if (nClouds ne 0) then begin
      c_pos = fltarr(3,nClouds)
      c_vel = fltarr(3,nClouds)
      c_id  = lonarr(nClouds)
    endif
    
    ;read blocks for all particles
    readu,1, pos
    readu,1, vel
    readu,1, id

    ; masses for variable mass particles (none for now)
    ind=where((h.nPartTot gt 0) and (h.massTable eq 0))
    
    if ind[0] ne -1 then begin
      nMass = total(h.nPartTot[ind])
      mass = fltarr(nMass)
      readu,1,mass
    endif      
    
    ;no u, rho, hsml (sph gas only)
    
    ; split star/bulge and cloud particles
    s_pos[0,*] = pos[0,0 : nStars+nBulge-1]
    s_pos[1,*] = pos[1,0 : nStars+nBulge-1]    
    s_pos[2,*] = pos[2,0 : nStars+nBulge-1]
    
    s_vel[0,*] = vel[0,0 : nStars+nBulge-1]    
    s_vel[1,*] = vel[1,0 : nStars+nBulge-1] 
    s_vel[2,*] = vel[2,0 : nStars+nBulge-1]    
    
    s_id = id[0 : nStars-1] 
    
    if (nClouds ne 0) then begin
      c_pos[0,*] = pos[0,nStars+nBulge : nStars+nBulge+nClouds-1]
      c_pos[1,*] = pos[1,nStars+nBulge : nStars+nBulge+nClouds-1]    
      c_pos[2,*] = pos[2,nStars+nBulge : nStars+nBulge+nClouds-1]
    
      c_vel[0,*] = vel[0,nStars+nBulge : nStars+nBulge+nClouds-1]    
      c_vel[1,*] = vel[1,nStars+nBulge : nStars+nBulge+nClouds-1] 
      c_vel[2,*] = vel[2,nStars+nBulge : nStars+nBulge+nClouds-1]    
    
      c_id = id[nStars+nBulge : nStars+nBulge+nClouds-1] 
    endif
    
  close,1

  return, h
end

; loadMCTracer: OLD-UNUSED
;
; inputs:  path, type ("old"/"new")
; outputs: t,x,y,z of tracer MC
;

function loadMCTracer, logPath, logType

  ; data structs     
  pOld = { name: '', $
           time: 0.0, $
           mass: 0.0, $
           x:    0.0, $
           y:    0.0, $
           z:    0.0, $
           unk:  0 $
         }
         
   pNew = { time: 0.0, $
            die:  0.0, $
            x:    0.0, $
            y:    0.0, $
            z:    0.0  $
           }

  ; load file and reformat
  if (logType eq "old") then begin
    pts = loadCSV(0, logPath, '')
    
    p = replicate(pOld,n_elements(pts))
    
    for i=0,n_elements(pts)-1 do begin
      tStr = strsplit(pts[i]," ",/extract)
      p[i].name = tStr[0]
      p[i].time = tStr[1]
      p[i].mass = tStr[2]
      p[i].x    = tStr[3]
      p[i].y    = tStr[4]
      p[i].z    = tStr[5]
      p[i].unk  = tStr[6]
    endfor
  endif
  if (logType eq "new") then begin
    pts = loadCSV(0, logPath, '')
    
    p = replicate(pNew,n_elements(pts))
    
    for i=0,n_elements(pts)-1 do begin
      tStr = strsplit(pts[i]," ",/extract)
      
      p[i].time = tStr[5]
      p[i].die  = strmid(tStr[7],0,strlen(tStr[7])-1)
      p[i].x    = strmid(tStr[10],1)
      p[i].y    = tStr[11]
      p[i].z    = strmid(tStr[12],0,strlen(tStr[12])-1)
    endfor
  endif
  
  return, p
end

; loadMCs
;
; load output file created by outputMolClouds()
; complete details on MCs only with fine time resolution

function loadMCs, filePath, i, header

  nHeaderLines = 1
  filename = filePath + "output/MCs_" + str(i) + ".txt"
  
  ptStruct = { MC_index: 0, $
               task:     0, $
               ind:      0, $
               TimeBorn: 0.0, $
               LifeTime: 0.0, $
               mass:     0.0, $
               x:        0.0, $
               y:        0.0, $
               z:        0.0, $
               velx:     0.0, $
               vely:     0.0, $
               velz:     0.0}
               
  ;if not (file_test(filename)) then return,0
  
  pts = loadCSV(nHeaderLines, filename, ptStruct, header=header)

  header = strsplit(header, " ",/extract)

  return, pts
end

; loadStars
;
; load output file created by outputStarPositions()
; (x,y,z) float triplets for all stars -> 12 bytes per star
; header=NumPart (4 byte float)

function loadStars, filePath, i

  ptStruct = { x:        0.0, $
               y:        0.0, $
               z:        0.0  $
               }  

  ; load Stars_x.bin (after catstars)
  filename = filePath + "output/Stars_" + str(i) + ".bin"
  
  if file_test(filename) eq 1 then begin
    pts = loadBinary(filename, ptStruct)
    return, pts
  endif

  ; if not, load Stars_x_* (before catstars)
  fileBase = filePath + "output/Stars_" + str(i) + "_"
  
  if file_test(fileBase+'0') eq 1 then begin
    pts = loadBinarySequence(fileBase, ptStruct)
    return, pts
  endif
  
  return,0
end

; loadFourierSaved (OLD-UNUSED)

function loadFourierSaved, filePath, i

  fileName = filePath + "output/fourierSaved_" + str(i) + ".txt"
  
  ;templates
  header = { nModes:0.0, nRadialBins:0.0, time: 0.0 }
             
  POINTS = 0
  count  = 0

  ;open file
  openR, lun, fileName, /GET_LUN

  ;read header
  readF, lun, header

  ;setup arrays
  tempA = fltarr(header.nRadialBins)
  A = fltarr(header.nRadialBins, header.nModes)
  r = fltarr(header.nRadialBins)
  
  ;read radial bins
  readF, lun, r
  
  ;read fourier amplitudes
  for i=0,header.nModes-1 do begin
    readF, lun, tempA
    A[*,i] = tempA
  endfor

  ;close handle
  close, lun
  free_lun, lun
  
  ;return struct
  f = { nModes:header.nModes, nRadialBins:header.nRadialBins, time:header.time, A:A, r:r}
  return,f

end

; loadMDGCurve
;
; load curve.txt output from MakeDiskGalaxy (vel curve and toomre's Q)
; schema:
; ------------------------------------------
; cc
; v200
; lambda
; md (jd?)
; gasfraction (?)
; diskheight
; bulgesize (?)
; 
; v200 (again?)
; H
; 
; POINTS
; R[POINTS]
; VC2 comp_Dphi_R(R,0) [POINTS]
; VC2 comp_Dphi_R_disk_tree(R,0) [POINTS]
; VC2 comp_DPhi_R_halo(R,0) [POINTS]
; VC2 comp_Dphi_R_bulge(R,0) [POINTS]
; 
; count
; list_R[count]
; Q[count]
; ------------------------------------------

function loadMDGCurve, filePath, vcs, qs

  ;templates
  header = { c:0.0, v200:0.0, lambda: 0.0, md: 0.0, gasFraction: 0.0, diskHeight: 0.0, $
             bulgeSize: 0.0, empty1: '', v200b: 0.0, H: 0.0, empty2: ''}
             
  POINTS = 0
  count  = 0

  ;open file
  openR, lun, filePath, /GET_LUN

  ;read header
  readF, lun, header
  
  ;read VC header
  readF, lun, POINTS
  
  ;create VC arrays
  vcs = {R: fltarr(POINTS), VC2_1: fltarr(POINTS), VC2_2: fltarr(POINTS), $
                            VC2_3: fltarr(POINTS), VC2_4: fltarr(POINTS)}
  
  ;read VC
  readF,lun,vcs
  
  ;read Q header
  readF,lun,''
  readF,lun,count
  
  ;create Q arrays
  qs = {list_R: fltarr(count), Q: fltarr(count)}

  ;read Q
  readF,lun,qs
  
  ;close handle
  close, lun
  free_lun, lun

  ; Sigma0 = (M_DISK) / (2 * PI * H * H);
  ; Q[i] = RadialDispersionFactor * sqrt(VelDispRz_disk[i][0]) * 
  ;        sqrt(epi_kappa2[i]) / (3.36 * G * Sigma0 * exp(-list_R[i] / H));

  ; #define  RSIZE 512
  ; #define  ZSIZE 512

  ; Baselen = 0.001 * H;
  ; LL = 10.0*R200;
  ; list_R[i] = exp(log(Baselen) + (i - 1) * (log(LL) - log(Baselen)) / (RSIZE - 1));
  ; list_z[i] = exp(log(Baselen) + (i - 1) * (log(LL) - log(Baselen)) / (ZSIZE - 1));

  return, header
end

; loadVelDispDump

;  fd = fopen("deldisp.dat", "w");
;  si = RSIZE + 1;
;  fwrite(&si, sizeof(int), 1, fd);
;  si = ZSIZE + 1;
;  fwrite(&si, sizeof(int), 1, fd);
;  fwrite(list_R, sizeof(double), RSIZE + 1, fd);
;  fwrite(list_z, sizeof(double), ZSIZE + 1, fd);
;  fwrite(&VelDispRz_disk[0][0], sizeof(double), (ZSIZE + 1) * (RSIZE + 1), fd);
;  fwrite(&VelDispPhi_disk[0][0], sizeof(double), (ZSIZE + 1) * (RSIZE + 1), fd);
;  fwrite(&VelStreamPhi_disk[0][0], sizeof(double), (ZSIZE + 1) * (RSIZE + 1), fd);
;  fwrite(&epi_gamma2[0], sizeof(double), RSIZE + 1, fd);
;  fwrite(&epi_kappa2[0], sizeof(double), RSIZE + 1, fd);
;  fclose(fd);

function loadVelDispDump, filePath, velDisp

  header = {rSize:0l, zSize:0l}

  ;open file
  openR, lun, filePath, /GET_LUN

  ;read header
  readU, lun, header

  ;create arrays
  velDisp = {list_R:            dblarr(header.rSize)               ,$
             list_Z:            dblarr(header.zSize)               ,$
             velDispRz_disk:    dblarr(header.zSize,header.rSize)  ,$
             velDispPhi_disk:   dblarr(header.zSize,header.rSize)  ,$
             velStreamPhi_disk: dblarr(header.zSize,header.rSize)  ,$
             epi_gamma2:        dblarr(header.rSize)               ,$
             epi_kappa2:        dblarr(header.rSize)                $
            }
  
  ;read VC
  readU,lun,velDisp

  ;close handle
  close, lun
  free_lun, lun

  return, header
end

; convertG2K():
; convert Gadget/Arepo snapshot to Krakatoa PRT format

pro convertG2K

  basePath    = '/n/home07/dnelson/spirals/cloudbias/RUN_test/output/snap_'
  ;basePath    = '/n/home07/dnelson/spirals/lambdacrit/sims/LC_10m_disk_6/output/snap_'
  prtFilename = '/n/home07/dnelson/test.100k.prt'
  
  m = 5
  
  nChannels = 1 ;pos, fltarr(3)

  ;load snapshot
  h = loadSnapshot(basePath,m,s_pos,s_vel,s_id,c_pos,c_vel,c_id)

  nPart = n_elements(s_id)
  print,'dumping nPart: ',nPart

  ;open PRT1 (uncompressed header)
  ;-------------------------------
  openw,1,prtFilename+".1"
  
  ;magic number
  writeu,1,'C0'XB ;192
  writeu,1,'PRT'
  writeu,1,'0D'XB ;\r
  writeu,1,'0A'XB ;\n
  writeu,1,'1A'XB ;26
  writeu,1,'0A'XB ;\n

  ;header length
  writeu,1,56L
  
  ;ident string
  writeu,1,'Extensible Particle Format'
  writeu,1,[0,0,0] ;null termination
  
  ;version number
  writeu,1,1L
  
  ;particle count
  writeu,1,long64(nPart)

  ;reserved
  writeu,1,4L
  
  ;number of channels
  writeu,1,long(nChannels)
  
    ;channel 1
    ;---------
    
    ;channel entry length = 44
    writeu,1,long(44)
    
    ;channel name
    writeu,1,'Position'
    writeu,1,bytarr(32-strlen('Position')) ;pad to 32 bytes
    
    ;data type: float32
    writeu,1,long(4)
    
    ;arity: 3
    writeu,1,long(3)
    
    ;data offset: 0
    writeu,1,long(0)
    
    ;channel 2
    ;---------
    ;empty


  ;open PRT2 (compressed data)
  ;---------------------------
  close,1
  openw,1,prtFilename+".2",/compress
  
  writeu,1,s_pos ;[3,nPart] = x1,y1,z1,x2,y2,z2,...
  
  ;close PRT
  close,1

  ;open PRT3 (data uncompressed debug)
  ;----------------------------------
  openw,1,prtFilename+".3"

  writeu,1,s_pos

  close,1
  
  ;open CSV (debug)
  ;---------------
  openw,1,prtFilename+".csv"
  for i=0L,nPart-1 do begin
    lineStr = strcompress(string(s_pos[0,i],format='(f20.16)'),/remove_all) + "," + $
              strcompress(string(s_pos[1,i],format='(f20.16)'),/remove_all) + "," + $
              strcompress(string(s_pos[2,i],format='(f20.16)'),/remove_all)
    printf,1,lineStr
  endfor
  close,1

  ;cat PRT
  cmd = 'cat '+prtFilename+".1 "+prtFilename+".2 > "+prtFilename
  spawn,cmd
  ;cmd = 'rm '+prtFilename+".1"
  ;spawn,cmd
  ;cmd = 'rm '+prtFilename+".2"
  ;spawn,cmd

end
