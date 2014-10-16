; zoomEvo.pro
; 'zoom project' evolution of halo properties with redshift
; dnelson oct.2014

pro tempDo1

  hInds = [3,4,5,6,9]
  resLevels = [9,10]
  
  foreach hInd,hInds do begin
    foreach resLevel,resLevels do begin
      sP = simParams(res=resLevel,hInd=hInd,run='zoom_20mpc',redshift=2.0)
      x = mergerTree(sP=sP,makeNum=60)
    endforeach
  endforeach

end

pro tempDo2

  hInds = [7]
  resLevels = [9,10,11]
  
  foreach hInd,hInds do begin
    foreach resLevel,resLevels do begin
      sP = simParams(res=resLevel,hInd=hInd,run='zoom_20mpc',redshift=2.0)
      x = mergerTree(sP=sP,makeNum=60)
    endforeach
  endforeach

end

pro tempDo3

  hInds = [8]
  resLevels = [9,10,11]
  
  foreach hInd,hInds do begin
    foreach resLevel,resLevels do begin
      sP = simParams(res=resLevel,hInd=hInd,run='zoom_20mpc',redshift=2.0)
      x = mergerTree(sP=sP,makeNum=60)
    endforeach
  endforeach

end

; plotZoomMassEvo(): plot time evolution of various masses

pro plotZoomMassEvo
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  hInds     = [0,1,7,8] ; TODO: h2 has serious problems (cannot make L9 galcat @ snap59, maybe stop at snap 58?)
  ; similar problem with h3L10
  resLevels = [9,10,11]
  newSaves  = 0 ; override existing saves
  
  ; plot config
  xrange_z    = [1.85,6.15] ; redshift
  yrange_mass = [10.0,12.5]   ; log halo mass
  
  lines   = [1,2,0] ; one per resLevel
  cInd    = 1       ; color index
  psym    = 0       ; plot symbol
  symsize = 0.9     ; symbol size
  
  ; load
  foreach hInd,hInds do begin
    hLoc = {}
        
    foreach resLevel,resLevels do begin
      ; get list of snapshots for this halo in the redshift range
      rLoc = {}
      
      if resLevel eq 11 and total(hInd eq [2,3,4,5,6,9]) gt 0 then continue
      
      sP = simParams(run='zoom_20Mpc',res=resLevel,hInd=hInd,redshift=2.0)     
      saveFilename = sP.derivPath + 'binnedVals/massEvo.' + sP.savPrefix + str(sP.res) + '.sav'
      
      if ~file_test(saveFilename) or newSaves eq 1 then begin
                  
        ; halo tracking
        tarInd = zoomTargetHalo(sP=sP)
        mts    = mergerTreeSubset(sP=sP,/verbose)
        nSnaps = n_elements(mts.snaps)
        
        ; allocate arrays
        mass = { total       : fltarr(nSnaps)        ,$ ; primary subfind, all types
                 total_gas   : fltarr(nSnaps)        ,$ ; primary subfind, gas
                 total_stars : fltarr(nSnaps)        ,$ ; primary subfind, stars
                 halo_gas    : fltarr(nSnaps)         } ; primary subfind+radial cuts, gas
        rvir = fltarr(nSnaps)
                
        if max(mts.hInd[*,tarInd]) gt 0 then print,'Warning: Target halo becomes not most massive in the past.'
        if mts.hMinSnap[tarInd] ne -1 then message,'Error: Target halo not tracked all the way.'
        
        ; loop over snaps, calculate mass in each, save
        foreach snap,mts.snaps,m do begin
          ; load group catalog
          sP.snap = snap
          gc = loadGroupCat(sP=sP,/skipIDs)
          gcInd = mts.hInd[m,tarInd] ; for this snapshot
          
          rvir[m] = gc.group_r_crit200[gc.subgroupGrNr[gcInd]]
          
          print,' h'+str(hInd)+'L'+str(resLevel)+': ['+string(snap,format='(I2)')+'] gcInd = ' + str(gcInd)
          
          mass.total[m]       = gc.subgroupMass[gcInd]
          mass.total_gas[m]   = gc.subgroupMassType[partTypeNum('gas'),gcInd]
          mass.total_stars[m] = gc.subgroupMassType[partTypeNum('stars'),gcInd]
          
          ; load positions of gas, make radial cut for halo selection
          ; TODO
          
          ; calculate masses of halo gas, based on temp/vrad cuts
          mass.halo_gas[m] = 0.0 ; TODO
          
        endforeach ;snaps
        
        ; units and little h factors
        mass.total       = codeMassToLogMsun( mass.total / units.hubbleParam )
        mass.total_gas   = codeMassToLogMsun( mass.total_gas / units.hubbleParam )
        mass.total_stars = codeMassToLogMsun( mass.total_stars / units.hubbleParam )
        mass.halo_gas    = codeMassToLogMsun( mass.halo_gas / units.hubbleParam )
        
        rvir /= units.hubbleParam ; still comoving
        
        ; save all values for this halo+resLevel combination
        rLoc = mod_struct( rLoc, 'sP', sP )
        rLoc = mod_struct( rLoc, 'redshifts', 1/mts.times-1 )
        rLoc = mod_struct( rLoc, 'mass', mass )     
       
        save,rLoc,filename=saveFilename
        print,'Saved TEMP file: ['+saveFilename+']'
      endif else begin
        print,'RESTORE: ['+saveFilename+']'
        restore,saveFilename
      endelse
    
      ; add this resolution level to the halo struct
      hLoc = mod_struct( hLoc, 'L'+str(resLevel), rLoc )     
     
    endforeach ;resLevels
    
    ; add this halo to the save struct
    evo = mod_struct( evo, 'h'+str(hInd), hLoc )
    
  endforeach ;hInds
  
  ; plot setup
  plotStr = 'h'
  hNames  = 'h' + str(hInds)
  hColors = []
  rNames  = 'L' + str(resLevels)
  rLines  = lines[0:n_elements(resLevels)-1]

  foreach hInd,hInds do plotStr += str(hInd)
  foreach hInd,hInds,i do hColors = [hColors,evo.(i).(0).sP.colors[cInd]]
  plotStr += '_L'
  foreach resLevel,resLevels do plotStr += str(resLevel)

  ; plot (1) - total mass vs redshift, all halos, all resolutions (single panel)
  start_PS, evo.(0).(0).sP.plotPath + 'massEvoTest.' + plotStr + '.eps'
  
    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_mass,/xs,/ys,$
      xtitle="Redshift",ytitle=textoidl("M_{halo} [_{ }log M_{sun }]")
    
    for i=0,n_tags(evo)-1 do begin
      for j=0,n_tags(evo.(i))-1 do begin
        xx = evo.(i).(j).redshifts
        yy = evo.(i).(j).mass.total
        co = evo.(i).(j).sP.colors[cInd]
        
        cgPlot,xx,yy,psym=psym,symsize=symsize,line=lines[j],color=co,/overplot
      endfor ;n_tags(evo.(i)),j
    endfor ;n_tags(evo),i
    
    legend,hNames,textcolor=hColors,/top,/right
    legend,rNames,linestyle=rLines,/top,/left
      
      
  end_PS
  
  ; plot (2)
  stop

end
