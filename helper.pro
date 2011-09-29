; helper.pro
; helper functions
; dnelson feb.2011

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
            colors : strarr(21)                               ,$
            
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
  
  units.rhoCrit = 3.0 * units.H0^2.0 / (8.0*!pi*units.G) ;code

  ; color list
  units.colors = ['black','blue','green','red','cyan','magenta','gray','orange', $
                  'brown','black','sienna','chartreuse','violet','papaya','yellow','aquamarine', $
                  'firebrick', 'rosy brown', 'gold', 'olive', 'saddle brown']

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

; writeICFile():
pro writeICFile, fOut, pos, vel, id, mass, u

  ; create header
  npart    = lonarr(6)  
  massarr  = dblarr(6)
  npartall = lonarr(6)

  npart(0)    = n_elements(id)
  npartall(0) = n_elements(id)

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
  
  ; remove non-IC data (for some reason these are expected IO blocks even for ICs, just zero them)
  ;s1 = mod_struct(s.PARTTYPE0,'DENSITY',/delete)
  ;s1 = mod_struct(s1,'SMOOTHINGLENGTH',/delete)
  ;s1 = mod_struct(s1,'VOLUME',/delete)
  ;s = mod_struct(s,'PARTTYPE0',s1)
  
  s.PARTTYPE0.DENSITY._DATA[*]         = 0.0
  s.PARTTYPE0.SMOOTHINGLENGTH._DATA[*] = 0.0
  s.PARTTYPE0.VOLUME._DATA[*]          = 0.0
  
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

; str()

function str, tString
  return, strcompress(string(tString),/remove_all)
end

function isnumeric, input
  on_ioerror, false
  test = double(input)
  return, 1
  false: return, 0
end

; startPS, endPS: my version

pro start_PS, filename, xs=xs, ys=ys

  if not keyword_set(xs) then xs=7.5
  if not keyword_set(ys) then ys=5.0

  PS_Start, FILENAME=filename, /nomatch, /quiet, bits_per_pixel=8, color=1, $
            /encapsulated, decomposed=0, xs=xs, ys=ys, /inches, font=0, tt_font='Times' ;3/2  
            
  ;!p.charsize = 1.2
  !p.thick    = 3.0          
            
end

pro end_PS, pngResize=pngResize, deletePS=deletePS

 ;PNG size=[xs,ys]*300*(resize/100)

  if not keyword_set(pngResize) then $
    PS_End
    
  if keyword_set(pngResize) then $
    PS_End, /PNG, Delete_PS=deletePS, Resize=pngResize

end