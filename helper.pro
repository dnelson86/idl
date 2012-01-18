; helper.pro
; helper functions
; dnelson feb.2011
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
  units.UnitTime_in_s = units.UnitLength_in_cm / units.UnitVelocity_in_cm_per_s
  units.UnitDensity_in_cgs  = units.UnitMass_in_g / units.UnitLength_in_cm^3.0
  units.UnitPressure_in_cgs = units.UnitMass_in_g / units.UnitLength_in_cm / units.UnitTime_in_s^2.0
  units.UnitEnergy_in_cgs   = units.UnitMass_in_g * units.UnitLength_in_cm^2.0 / units.UnitTime_in_s^2.0
  
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

; writeICFile(): write old Gadget format IC file with gas particles and tracers
;                each partX struct should contain {id,pos,vel,mass,u} in the usual format

pro writeICFile, fOut, part0=part0, part1=part1, part2=part2

  ; arrays
  pos  = []
  vel  = []
  id   = []
  mass = []
  u    = []

  ; create header
  npart    = lonarr(6)  
  massarr  = dblarr(6)
  npartall = lonarr(6)

  ; add to particle counts and concat arrays
  if keyword_set(part0) then begin
    ; GAS
    npart(0)    = n_elements(part0.id)
    npartall(0) = n_elements(part0.id)
    
    pos  = [[pos], [part0.pos]]
    vel  = [[vel], [part0.vel]]
    id   = [id,    part0.id]
    mass = [mass,  part0.mass]
    u    = [u,     part0.u]
  endif
  
  if keyword_set(part1) then begin
    npart(1)    = n_elements(part1.id)
    npartall(1) = n_elements(part1.id)
    
    pos  = [[pos], [part1.pos]]
    vel  = [[vel], [part1.vel]]
    id   = [id,    part1.id]
    mass = [mass,  part1.mass]
    u    = [u,     part1.u]
  endif
  
  if keyword_set(part2) then begin
    ; TRACER
    npart(2)    = n_elements(part2.id)
    npartall(2) = n_elements(part2.id)
    
    pos  = [[pos], [part2.pos]]
    vel  = [[vel], [part2.vel]]
    id   = [id,    part2.id]
    mass = [mass,  part2.mass]
    ;u    = [u,     part2.u] ; u not expected in input for tracer
  endif

  time          = 0.0D
  redshift      = 0.0D
  flag_sfr      = 0L
  flag_feedback = 0L
  
  bytesleft = 136
  la        = intarr(bytesleft/2)

  ; write IC file
  openw,1,fOut,/f77_unformatted
  writeu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall,la

  writeu,1, pos
  writeu,1, vel
  writeu,1, id
  writeu,1, mass
  writeu,1, u                       ; internal energy per unit mass
  close,1
  
  print,'wrote ',fOut

end

; writeICFileHDF5(): GAS ONLY
pro writeICFileHDF5, fOut, boxSize, pos, vel, id, massOrDens, u

  ; load hdf5 template
  templatePath = '/n/home07/dnelson/make.ics/ArepoTemplate.hdf5'
  
  s = h5_parse(templatePath, /read)

  ; modify base
  s._NAME    = fOut
  s._FILE    = fOut
  s._COMMENT = "dnelson IC gen"
  
  ; modify HEADER
  s.HEADER._FILE                      = fOut
  s.HEADER.NUMPART_THISFILE._DATA     = [n_elements(id),0,0,0,0,0]
  s.HEADER.NUMPART_TOTAL._DATA        = [n_elements(id),0,0,0,0,0]
  s.HEADER.MASSTABLE._DATA            = [0.0,0.0,0.0,0.0,0.0,0.0]
  s.HEADER.TIME._DATA                 = 0.0
  s.HEADER.REDSHIFT._DATA             = 0.0
  s.HEADER.BOXSIZE._DATA              = boxSize
  s.HEADER.NUMFILESPERSNAPSHOT._DATA  = 1
  s.HEADER.OMEGA0._DATA               = 0.0
  s.HEADER.OMEGALAMBDA._DATA          = 0.0
  s.HEADER.HUBBLEPARAM._DATA          = 1.0
  s.HEADER.FLAG_SFR._DATA             = 0
  s.HEADER.FLAG_COOLING._DATA         = 0  
  s.HEADER.FLAG_STELLARAGE._DATA      = 0
  s.HEADER.FLAG_METALS._DATA          = 0
  s.HEADER.FLAG_FEEDBACK._DATA        = 0
  s.HEADER.FLAG_DOUBLEPRECISION._DATA = 0
  
  s.HEADER.COMPOSITION_VECTOR_LENGTH._DATA = 0 ;?

  ; for some reason these are expected IO blocks even for ICs
  s1 = mod_struct(s.PARTTYPE0,'DENSITY',/delete)
  s1 = mod_struct(s1,'SMOOTHINGLENGTH',/delete)
  s1 = mod_struct(s1,'VOLUME',/delete)
  s = mod_struct(s,'PARTTYPE0',s1)

  ;s.PARTTYPE0.DENSITY._DATA[*]         = 0.0
  ;s.PARTTYPE0.SMOOTHINGLENGTH._DATA[*] = 0.0
  ;s.PARTTYPE0.VOLUME._DATA[*]          = 0.0
  
  ; modify data parameters
  s.PARTTYPE0._FILE                = fOut
  s.PARTTYPE0.COORDINATES._FILE    = fOut
  s.PARTTYPE0.VELOCITIES._FILE     = fOut
  s.PARTTYPE0.PARTICLEIDS._FILE    = fOut
  s.PARTTYPE0.INTERNALENERGY._FILE = fOut
    
  ; modify data
  s.PARTTYPE0.COORDINATES._DIMENSIONS    = [3,n_elements(id)]
  s.PARTTYPE0.COORDINATES._NELEMENTS     = n_elements(pos)
  s1 = mod_struct(s.PARTTYPE0.COORDINATES,'_DATA',pos) ;change _DATA size
  s2 = mod_struct(s.PARTTYPE0,'COORDINATES',s1) ;update PARTTYPE0 with child

  s.PARTTYPE0.VELOCITIES._DIMENSIONS     = [3,n_elements(id)]
  s.PARTTYPE0.VELOCITIES._NELEMENTS      = n_elements(vel)
  s1 = mod_struct(s.PARTTYPE0.VELOCITIES,'_DATA',vel)
  s2 = mod_struct(s2,'VELOCITIES',s1)
  
  s.PARTTYPE0.PARTICLEIDS._DIMENSIONS    = [n_elements(id)]
  s.PARTTYPE0.PARTICLEIDS._NELEMENTS     = n_elements(id)
  s1 = mod_struct(s.PARTTYPE0.PARTICLEIDS,'_DATA',id)
  s2 = mod_struct(s2,'PARTICLEIDS',s1)

  s.PARTTYPE0.INTERNALENERGY._DIMENSIONS = [n_elements(u)]
  s.PARTTYPE0.INTERNALENERGY._NELEMENTS  = n_elements(u)
  s1 = mod_struct(s.PARTTYPE0.INTERNALENERGY,'_DATA',u)
  s2 = mod_struct(s2,'INTERNALENERGY',s1)
  
  s.PARTTYPE0.MASSES._DIMENSIONS = [n_elements(massOrDens)]
  s.PARTTYPE0.MASSES._NELEMENTS  = n_elements(massOrDens)
  s1 = mod_struct(s.PARTTYPE0.MASSES,'_DATA',massOrDens)
  s2 = mod_struct(s2,'MASSES',s1)
  
  s = mod_struct(s,'PARTTYPE0',s2) ;import new PARTTYPE0 structure

  ; output
  h5_create, fOut, s

end

; one line utility functions: str(), isnumeric()

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

; startPS, endPS: my version

pro start_PS, filename, xs=xs, ys=ys

  if not keyword_set(xs) then xs=7.5
  if not keyword_set(ys) then ys=5.0

  PS_Start, FILENAME=filename, /nomatch, /quiet, bits_per_pixel=8, color=1, $
            /encapsulated, decomposed=0, xs=xs, ys=ys, /inches, font=0, tt_font='Times' ;3/2  
 
  !p.charsize = 1.5
  !p.thick    = 3.0          
            
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

; load routines for use
; ---------------------
@cosmoUtil
@cosmoLoad

@tracers
@tracersCosmo
;@tracersDisks
;@tracersShocktube

@cosmoVis
;@cosmoPlot
;@cosmoAnalysis

@arepoLoad
@arepoVis2D
;@arepoSphSym
