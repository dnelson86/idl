; cosmoExport.pro
; cosmological simulations - utility export functions (ascii)
; dnelson jul.2014

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
;
;   -----
;   to add a new column to the existing DB, where a text file 'in.txt' contains one value 
;   for each subhalo (linenumber=subhaloID+1), do:
;   
;   USE subgroups;
;   CREATE TEMPORARY TABLE foo (id_in INT AUTO_INCREMENT,value SMALLINT,PRIMARY KEY (id_in));
;   LOAD DATA LOCAL INFILE 'in.txt' INTO TABLE foo 
;     COLUMNS TERMINATED BY ' ' LINES TERMINATED BY '\n' (value);
;   UPDATE foo SET id_in=id_in-1 WHERE 1;
;   UPDATE snap_135 INNER JOIN foo ON snap_135.id=foo.id_in SET snap_135.value=foo.value;
;   DROP TEMPORARY TABLE foo;

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

; exportBlackholesAscii(): export all the blackhole info as text (for MySQL import)
;
;   % mysql --local-infile -u <username> -p <DatabaseName>
;
;   USE subgroups;
;   LOAD DATA LOCAL INFILE 'out.txt'
;   INTO TABLE blackholes_135
;   COLUMNS TERMINATED BY ' ' LINES TERMINATED BY '\n'
;   (id,particle_id,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,mass,
;    hosthalomass,cumegyinjection_qm,cummassgrowth_qm,bh_density,bh_mass,
;    bh_mass_bubbles,bh_mass_ini,bh_mdot,bh_pressure,bh_u,bh_progs,
;    subhalo_parent_id,image_x,image_y)
;   set point_xy = PointFromText(CONCAT('POINT(',image_x,' ',image_y,')'));
;
;   then add SPATIAL INDEX on point_xy:
;       ALTER TABLE `blackholes_135` ADD SPATIAL INDEX (point_xy);
;   then add FOREIGN KEY CONSTRAINT (NULL=True) on subhalo_parent_id -> (Subhalo135.id):
;       ALTER TABLE `blackholes_135` ADD FOREIGN KEY (`subhalo_parent_id`) 
;       REFERENCES `snap_135`(`id`) ON UPDATE NO ACTION ON DELETE NO ACTION;


pro exportBlackholesAscii
  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function simParams, loadGroupCat
  
  ; config
  sP         = simParams(res=1820,run='illustris',redshift=0.0)
  fileName   = 'out.txt'
  fields     = ['pos', 'vel', 'ids', 'mass', 'hosthalomass', $
                'cumegyinjection_qm', 'cummassgrowth_qm', 'bh_density', 'bh_mass', $
                'bh_mass_bubbles', 'bh_mass_ini', 'bh_mdot', 'bh_pressure', 'bh_progs', 'bh_u']
  
  ; load
  gc = loadGroupCat(sP=sP,/skipIDs)
  
  data = {}
  foreach field,fields do begin
    print, field
    load = loadSnapshotSubset(sP=sP,partType='bhs',field=field)
    data = mod_struct( data, field, load )
  endforeach
      
  ; make image (fof0 centered) coordinates
  image_x = reform( data.pos[0,*] )
  image_y = reform( data.pos[1,*] )
  image_x = image_x - gc.subgroupPos[0,0] + sP.boxSize*0.5
  image_y = image_y - gc.subgroupPos[1,0] + sP.boxSize*0.5

  correctPeriodicPosVecs, image_x, sP=sP
  correctPeriodicPosVecs, image_y, sP=sP
  
  ; create parent subhalo IDs for each blackhole
  bh_parents = gcRepParentIDs(sP=sP, gc=gc, partType='bh', /subhaloIDs)
  
  ; fill output buffer
  print,'Exporting ['+str(count)+'] BLACKHOLES.'
  
  outBuf = strarr(count)
  
  for i=0,count-1 do begin
    if i mod round(count*0.1) eq 0 then print,i,count
    
    if bh_parents[i] ge 0 then subgroupNumber = string(bh_parents[i], format='(I8)')
    if bh_parents[i] lt 0 then subgroupNumber = "\N      " ; MySQL foreign key null flag
    
    outBuf[i] = str(string(i,                 format='(I7)'))    + " " + $ ; index
                str(string(data.ids[i],       format='(I20)'))   + " " + $ ; ID
                str(string(data.pos[0,i],     format='(f12.6)')) + " " + $ ; pos_x
                str(string(data.pos[1,i],     format='(f12.6)')) + " " + $ ; pos_y
                str(string(data.pos[2,i],     format='(f12.6)')) + " " + $ ; pos_z
                str(string(data.vel[0,i],     format='(f12.6)')) + " " + $ ; vel_x
                str(string(data.vel[1,i],     format='(f12.6)')) + " " + $ ; vel_y
                str(string(data.vel[2,i],     format='(f12.6)')) + " " + $ ; vel_z
                str(string(data.mass[i],               format='(f14.8)')) + " " + $ ; mass
                str(string(data.hosthalomass[i],       format='(f14.6)')) + " " + $ ; hosthalomass
                str(string(alog10(data.cumegyinjection_qm[i]), format='(f12.6)')) + " " + $ ; cumegyinjection_qm
                str(string(alog10(data.cummassgrowth_qm[i]),   format='(f12.6)')) + " " + $ ; cummassgrowth_qm
                str(string(data.bh_density[i],         format='(f14.8)')) + " " + $ ; bh_density
                str(string(data.bh_mass[i],            format='(f14.8)')) + " " + $ ; bh_mass
                str(string(data.bh_mass_bubbles[i],    format='(f14.8)')) + " " + $ ; bh_mass_bubbles
                str(string(data.bh_mass_ini[i],        format='(f14.6)')) + " " + $ ; bh_mass_ini
                str(string(data.bh_mdot[i],            format='(f14.8)')) + " " + $ ; bh_mdot
                str(string(data.bh_pressure[i],        format='(f14.6)')) + " " + $ ; bh_pressure
                str(string(data.bh_u[i],               format='(f14.6)')) + " " + $ ; bh_u
                str(string(fix(data.bh_progs[i]),      format='(I4)'))    + " " + $ ; bh_progs
		    str(subgroupNumber)                                       + " " + $ ; parent SG number
                str(string(image_x[i],                 format='(f12.6)')) + " " + $ ; image_x
                str(string(image_y[i],                 format='(f12.6)'))           ; image_y
                
  endfor
  
  ; open ascii, write and close
  openW,lun,fileName,/get_lun
  
    printf,lun,outBuf
  
  close,lun
  free_lun,lun

  stop

end