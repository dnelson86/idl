; zoomEvo.pro
; 'zoom project' evolution of halo properties with redshift
; dnelson nov.2014

; plotZoomMassEvo(): plot time evolution of various masses

pro plotZoomMassEvo
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  hInds     = [0,1,7,8] ; TODO: h2L9 and h3L10 cannot make galcat @ snap59, maybe stop at snap 58?)
  resLevels = [9,10,11]
  newSaves  = 0 ; override existing saves
  
  ; plot config
  xrange_z          = [2.0,6.0]    ; redshift
  yrange_mass       = [10.0,12.1]  ; log msun, total halo mass
  yrange_mass_gas   = [9.5,10.8]   ; log msun
  yrange_mass_stars = [9.0,11.5]   ; log msun
  yrange_rad        = [100,350]    ; ckpc
  yrange_logT       = [5.4,6.2]    ; log K
  
  lines   = [1,2,0] ; line style, one per resLevel
  cInds   = [1,1,1] ; color index, one per resLevel
  psym    = 8       ; plot symbol
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
                 
        rVir = fltarr(nSnaps)
        tVir = fltarr(nSnaps)
                
        if max(mts.hInd[*,tarInd]) gt 0 then print,'Warning: Target halo becomes not most massive in the past.'
        if mts.hMinSnap[tarInd] ne -1 then message,'Error: Target halo not tracked all the way.'
        
        ; loop over snaps, calculate mass in each, save
        foreach snap,mts.snaps,m do begin
          ; load group catalog
          sP.snap = snap
          sP.redshift = 1/mts.times[m]-1
          
          gc = loadGroupCat(sP=sP,/skipIDs)
          gcInd = mts.hInd[m,tarInd] ; for this snapshot
          
          rVir[m] = gc.group_r_crit200[gc.subgroupGrNr[gcInd]]
          tVir[m] = codeMassToVirTemp( gc.subgroupMass[gcInd], redshift=sP.redshift, /log )
          
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
        
        rVir /= units.hubbleParam ; still comoving
        
        ; save all values for this halo+resLevel combination
        rLoc = mod_struct( rLoc, 'sP', sP )
        rLoc = mod_struct( rLoc, 'redshifts', 1/mts.times-1 )
        rLoc = mod_struct( rLoc, 'mass', mass )     
        rLoc = mod_struct( rLoc, 'rVir', rVir )
        rLoc = mod_struct( rLoc, 'tVir', tVir )
       
        save,rLoc,filename=saveFilename
        print,'Saved file: ['+saveFilename+']'
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
  foreach hInd,hInds,i do hColors = [hColors,evo.(i).(0).sP.colors[cInds[-1]]]
  plotStr += '_L'
  foreach resLevel,resLevels do plotStr += str(resLevel)

  ; plot (1) - total mass vs redshift, all halos, all resolutions (single panel)
  start_PS, evo.(0).(0).sP.plotPath + 'zoomMassEvo_total.' + plotStr + '.eps'
  
    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_mass,/xs,/ys,$
      xtitle="Redshift",ytitle=textoidl("M_{halo} [_{ }log M_{sun }]")
    
    for i=0,n_tags(evo)-1 do begin
      for j=0,n_tags(evo.(i))-1 do begin
        xx = evo.(i).(j).redshifts
        yy = evo.(i).(j).mass.total
        co = evo.(i).(j).sP.colors[cInds[j]]
        
        cgPlot,xx,yy,line=lines[j],color=co,/overplot
      endfor ;n_tags(evo.(i)),j
      
      ; endpoint symbol for highest resolution
      cgPlot,xx[0],yy[0],psym=psym,symsize=symsize,color=co,/overplot
    endfor ;n_tags(evo),i
    
    legend,hNames,textcolor=hColors,/bottom,/left
    legend,rNames,linestyle=rLines,/top,/right
      
  end_PS
  
  ; plot (2)
  start_PS, evo.(0).(0).sP.plotPath + 'zoomMassEvo_3x2.' + plotStr + '.eps', xs=12.0, ys=6.0
  
    pos = plot_pos(col=3,row=2,/gap)
    offset = [0,0,-0.02,0]
    
    ; (a) - total mass  
    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_mass,/xs,/ys,$
      xtitle="Redshift",ytitle=textoidl("M_{sub,total} [_{ }log M_{sun }]"),pos=pos[0]+offset,/noerase
    
    for i=0,n_tags(evo)-1 do begin
      for j=0,n_tags(evo.(i))-1 do begin
        xx = evo.(i).(j).redshifts
        yy = evo.(i).(j).mass.total
        co = evo.(i).(j).sP.colors[cInds[j]]
        
        cgPlot,xx,yy,line=lines[j],color=co,/overplot
      endfor ;n_tags(evo.(i)),j
      
      ; endpoint symbol for highest resolution
      cgPlot,xx[0],yy[0],psym=psym,symsize=symsize,color=co,/overplot
    endfor ;n_tags(evo),i
    
    legend,hNames,textcolor=hColors,/bottom,/left
      
    ; (b) - total gas mass
    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_mass_gas,/xs,/ys,$
      xtitle="Redshift",ytitle=textoidl("M_{sub,gas} [_{ }log M_{sun }]"),pos=pos[1]+offset,/noerase
    
    for i=0,n_tags(evo)-1 do begin
      for j=0,n_tags(evo.(i))-1 do begin
        xx = evo.(i).(j).redshifts
        yy = evo.(i).(j).mass.total_gas
        co = evo.(i).(j).sP.colors[cInds[j]]
        
        cgPlot,xx,yy,line=lines[j],color=co,/overplot
      endfor ;n_tags(evo.(i)),j
      
      ; endpoint symbol for highest resolution
      cgPlot,xx[0],yy[0],psym=psym,symsize=symsize,color=co,/overplot
    endfor ;n_tags(evo),i
    
    legend,rNames,linestyle=rLines,/top,/right
    
    ; (c) - total stellar mass
    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_mass_stars,/xs,/ys,$
      xtitle="Redshift",ytitle=textoidl("M_{sub,stars} [_{ }log M_{sun }]"),pos=pos[2]+offset,/noerase
    
    for i=0,n_tags(evo)-1 do begin
      for j=0,n_tags(evo.(i))-1 do begin
        xx = evo.(i).(j).redshifts
        yy = evo.(i).(j).mass.total_stars
        co = evo.(i).(j).sP.colors[cInds[j]]
        
        cgPlot,xx,yy,line=lines[j],color=co,/overplot
      endfor ;n_tags(evo.(i)),j
      
      ; endpoint symbol for highest resolution
      cgPlot,xx[0],yy[0],psym=psym,symsize=symsize,color=co,/overplot
    endfor ;n_tags(evo),i
    
    ; (d) - virial radius, comoving
    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_rad,/xs,/ys,$
      xtitle="Redshift",ytitle=textoidl("R_{vir} [_{ }ckpc_{ }]"),pos=pos[3]+offset,/noerase
    
    for i=0,n_tags(evo)-1 do begin
      for j=0,n_tags(evo.(i))-1 do begin
        xx = evo.(i).(j).redshifts
        yy = evo.(i).(j).rVir
        co = evo.(i).(j).sP.colors[cInds[j]]
        
        cgPlot,xx,yy,line=lines[j],color=co,/overplot
      endfor ;n_tags(evo.(i)),j
      
      ; endpoint symbol for highest resolution
      cgPlot,xx[0],yy[0],psym=psym,symsize=symsize,color=co,/overplot
    endfor ;n_tags(evo),i
    
    ; (e) - virial temperature, log K
    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_logT,/xs,/ys,$
      xtitle="Redshift",ytitle=textoidl("T_{vir} [_{ }log K_{ }]"),pos=pos[4]+offset,/noerase
    
    for i=0,n_tags(evo)-1 do begin
      for j=0,n_tags(evo.(i))-1 do begin
        xx = evo.(i).(j).redshifts
        yy = evo.(i).(j).tVir
        co = evo.(i).(j).sP.colors[cInds[j]]
        
        cgPlot,xx,yy,line=lines[j],color=co,/overplot
      endfor ;n_tags(evo.(i)),j
      
      ; endpoint symbol for highest resolution
      cgPlot,xx[0],yy[0],psym=psym,symsize=symsize,color=co,/overplot
    endfor ;n_tags(evo),i
    
    ; (f)
    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_logT,/xs,/ys,$
      xtitle="Redshift",ytitle=textoidl(""),pos=pos[5]+offset,/noerase
      
  end_PS
  
  stop

end
