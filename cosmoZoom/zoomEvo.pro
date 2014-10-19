; zoomEvo.pro
; 'zoom project' evolution of halo properties with redshift
; dnelson oct.2014

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

; plotZoomRadialProfiles() - 1d with all halos on same plot, and 2d with each halo in separate panel

pro plotZoomRadialProfiles
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  hInds     = [0,1,7,8]
  resLevels = [9,10,11]
  redshift  = 2.0
  newSaves  = 0 ; override existing saves
  
  ; binning config
  nBaseBins = 100 ; for base resolution level, increased at higher resolutions
  rBinMax   = 2.0 ; r/rvir
  cutSubS   = 1 ; remove substructures?
  
  ; plot config
  xrange_rad   = [0.0,2.0]   ; r/rvir
  xrange_log   = 0           ; plot in log? doesn't necessary have to equal binInLog  
  yrange_temp  = [5e4,1e7]   ; K
  yrange_dens  = [-28,-24]   ; log cgs (g/cm^3)
  yrange_vrad  = [-200,100]  ; km/s
  yrange_csize = [0.01,2.0]  ; physical kpc
  yrange_entr  = [1e7,1e9]   ; cgs (K cm^2)
  yrange_angm  = [1e3,1e5]   ; kpc km/s
  
  lines    = [1,2,0] ; line style, one per resLevel
  cInds    = [1,1,1] ; color index, one per resLevel
  radLines = [0.15,1.5]
  sK       = 1
  
  ; load
  foreach hInd,hInds do begin
    hLoc = {}
        
    foreach resLevel,resLevels do begin
      ; get list of snapshots for this halo in the redshift range
      rLoc = {}
      
      if resLevel eq 11 and total(hInd eq [2,3,4,5,6,9]) gt 0 then continue
      
      sP = simParams(run='zoom_20Mpc',res=resLevel,hInd=hInd,redshift=redshift)     
      saveFilename = sP.derivPath + 'binnedVals/radProfiles.h'+str(hInd)+'.' + sP.savPrefix + str(sP.res) + $
                     '.cutSubS-' + str(cutSubS) + '.sav'
      
      if ~file_test(saveFilename) or newSaves eq 1 then begin
                  
        ; adjust binning based on resolution (increase as cube root of number of particles)
        nRadBins = nBaseBins * 2.0^(sP.zoomLevel-2)
                  
        ; allocate arrays (4 per quantity = [mean,lower quartile,median,upper quartile])
        gas = { rBinCen : fltarr(nRadBins)        ,$ ; r/rvir
                temp    : fltarr(nRadBins,4)      ,$ ; gas temp, K
                dens    : fltarr(nRadBins,4)      ,$ ; gas density, cgs
                vrad    : fltarr(nRadBins,4)      ,$ ; radial velocity, km/s
                csize   : fltarr(nRadBins,4)      ,$ ; cellsize, physical kpc
                entr    : fltarr(nRadBins,4)      ,$ ; entropy, cgs
                angm    : fltarr(nRadBins,4)       } ; specific angular momentum, kpc km/s
                
        ; load group catalog, locate halo
        gc = loadGroupCat(sP=sP,readIDs=(cutSubS eq 1),skipIDs=(cutSubS eq 0))
        gcInd = zoomTargetHalo(sP=sP)
        
        rVir = gc.group_r_crit200[gc.subgroupGrNr[gcInd]]
        tVir = codeMassToVirTemp( gc.subgroupMass[gcInd], redshift=sP.redshift )
        pPos = gc.subgroupPos[*,gcInd]
        
        print,' h'+str(hInd)+'L'+str(resLevel)+': ['+string(sP.snap,format='(I2)')+'] gcInd = ' + str(gcInd) + $
              ' rVir = ' + string(rVir,format='(f5.1)') + ' tVir = ' + string(alog10(tVir),format='(f3.1)')
        
        ; load positions of gas, make radial cut for halo selection
        pos = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
          
        rad_pri  = periodicDists(pPos,pos,sP=sP)
        rad_pri /= rVir
        
        wRad = where(rad_pri le rBinMax*1.5,countRad)
        print,'  within '+string(rBinMax*1.5,format='(f3.1)')+' rvir, have ['+str(countRad)+$
              '] of ['+str(n_elements(rad_pri))+']'
              
        ; remove substructures? load ids and make secondary list
        if cutSubS eq 1 then begin
          ; verify our target ID is a primary (otherwise...)
          priGCIDs = gcIDList(gc=gc,select='pri')
          if total(gcInd eq priGCIDs) eq 0 then message,'Error'
          
          ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids',inds=wRad)
          
          ; get ids of all gas cells in all substructures other than the target
          ; (this differs from ids in the target, as it includes cells not in any subfind group)
          myPIDs  = gcPIDList(gc=gc,valGCids=[gcInd],partType='gas')
          subPIDs = gcPIDList(gc=gc,select='all',partType='gas')
          remPIDs = removeIntersectionFromB(myPIDs,subPIDs)
          if n_elements(remPIDs) ne n_elements(subPIDs)-n_elements(myPIDs) then message,'Error'
          
          ; remove the intersection of (satPIDs,ids) from wRad
          calcMatch,remPIDs,ids,ind1,ind2,count=count
          
          all = bytarr(countRad)
          if count gt 0 then all[ind2] = 1B ; flag substructure gas cells, then select against
          
          ; update wRad and countRad
          wRad = wRad[ where(all eq 0B, countRad) ]
          
          print,'  substructures ['+str(count)+'] removed! now have ['+str(countRad)+'] in radius.'
        endif
          
        pos     = pos[*,wRad]
        rad_pri = rad_pri[wRad]
          
        ; load further gas properties, handle little h and unit conversions
        scalefac = 1.0 / (1.0 + sP.redshift)
        
        u     = loadSnapshotSubset(sP=sP,partType='gas',field='u',inds=wRad)
        nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec',inds=wRad)
        dens  = loadSnapshotSubset(sP=sP,partType='gas',field='density',inds=wRad)
        vel   = loadSnapshotSubset(sP=sP,partType='gas',field='vel',inds=wRad)
        csize = loadSnapshotSubset(sP=sP,partType='gas',field='cellsize',inds=wRad)
        sfr   = loadSnapshotSubset(sP=sP,partType='gas',field='sfr',inds=wRad)
    
        temp  = convertUtoTemp(u, nelec)
        entr  = calcEntropyCGS(u, dens, sP=sP)
        dens  = codeDensToPhys(dens, sP=sP, /cgs)
        csize = csize * (scalefac^3.0 / units.HubbleParam^3.0)
        
        w = where(sfr gt 0.0,countSfr) ; truncate eEOS temperature down to minimum of h0L10@z=2 (3000K)
        if countSfr gt 0 then temp[w] = 3000.0
        
        ; vrad calculation: replace coordinates by relative coordinates (radial vectors) to direct parent
        for i=0,2 do begin
          pos_rel = reform(pos[i,*] - pPos[i,gcInd])
          correctPeriodicDistVecs, pos_rel, sP=sP
          pos[i,*] = pos_rel
        endfor
      
        vel *= sqrt(scalefac)

        ; subgroupVel already peculiar (no scalefac correction needed)
        vrad = ((vel[0,*] - gc.subgroupVel[0,gcInd]) * pos[0,*] + $ 
                (vel[1,*] - gc.subgroupVel[1,gcInd]) * pos[1,*] + $
                (vel[2,*] - gc.subgroupVel[2,gcInd]) * pos[2,*]) $
                    / (rad_pri*rVir)
                
        ; angm calculation: mean magnitude of specific angular momentum = rvec x vel
        ; make velocities relative to bulk halo motion
        vel[0,*] = reform(vel[0,*] - gc.subgroupVel[0,gcInd])
        vel[1,*] = reform(vel[1,*] - gc.subgroupVel[1,gcInd])
        vel[2,*] = reform(vel[2,*] - gc.subgroupVel[2,gcInd])
        
        ; angular momentum magnitude
        jvec = fltarr(3,countRad)
        jvec[0,*] = pos[1,*] * vel[2,*] - pos[2,*] * vel[1,*]
        jvec[1,*] = pos[2,*] * vel[0,*] - pos[0,*] * vel[2,*]
        jvec[2,*] = pos[0,*] * vel[1,*] - pos[1,*] * vel[0,*]
        jnorm = reform( sqrt(jvec[0,*]*jvec[0,*] + jvec[1,*]*jvec[1,*] + jvec[2,*]*jvec[2,*]) )
                
        ; setup binning
        rBinMin = 0.0
        rSize = (rBinMax-rBinMin)/nRadBins
        
        ; binning
        for i=0,nRadBins-1 do begin
          rLocMin = rBinMin + rSize*i
          rLocMax = rBinMin + rSize*(i+1)
          w = where(rad_pri ge rLocMin and rad_pri lt rLocMax,count)
          
          gas.rBinCen[i] = mean([rLocMin,rLocMax])
          
          print,'  ['+string(i,format='(I3)')+'] rLocMin = '+string(rLocMin,format='(f5.3)')+$
                ' rLocMax = '+string(rLocMax,format='(f5.3)')+' rBinCen = '+string(gas.rBinCen[i],format='(f5.3)')+$
                ' ['+str(count)+']'
          
          if count eq 0 then continue
          
          gas.temp[i,0]   = mean( temp[w] )
          gas.temp[i,1:3] = percentiles( temp[w] )
          gas.dens[i,0]   = mean( dens[w] )
          gas.dens[i,1:3] = percentiles( dens[w] )
          gas.vrad[i,0]   = mean( vrad[w] )
          gas.vrad[i,1:3] = percentiles( vrad[w] )
          
          gas.csize[i,0]   = mean( csize[w] )
          gas.csize[i,1:3] = percentiles( csize[w] )
          gas.entr[i,0]    = mean( entr[w] )
          gas.entr[i,1:3]  = percentiles( entr[w] )
          gas.angm[i,0]    = mean( jnorm[w] )
          gas.angm[i,1:3]  = percentiles( jnorm[w] )
        endfor
        
        ; save all values for this halo+resLevel combination
        halo = { tVir : tVir, rVir : rVir, pPos : pPos }
        
        rLoc = mod_struct( rLoc, 'sP', sP )
        rLoc = mod_struct( rLoc, 'gas', gas )
        rLoc = mod_struct( rLoc, 'halo', halo )
       
        save,rLoc,filename=saveFilename
        print,'Saved file: ['+strmid(saveFilename,strlen(sP.derivPath))+']'
      endif else begin
        print,'RESTORE: ['+strmid(saveFilename,strlen(sP.derivPath))+']'
        restore,saveFilename
      endelse
    
      ; add this resolution level to the halo struct
      hLoc = mod_struct( hLoc, 'L'+str(resLevel), rLoc )     
     
    endforeach ;resLevels
    
    ; add this halo to the save struct
    rp = mod_struct( rp, 'h'+str(hInd), hLoc )
    
  endforeach ;hInds
  
  ; plot setup
  plotStr = 'h'
  hNames  = 'h' + str(hInds)
  hColors = []
  rNames  = 'L' + str(resLevels)
  rLines  = lines[0:n_elements(resLevels)-1]

  foreach hInd,hInds do plotStr += str(hInd)
  foreach hInd,hInds,i do hColors = [hColors,rp.(i).(0).sP.colors[cInds[-1]]]
  plotStr += '_L'
  foreach resLevel,resLevels do plotStr += str(resLevel)
  plotStr += '_cutSubS-' + str(cutSubS)

  if xrange_log eq 1 then xrange_rad[0] = 0.01
  if xrange_log eq 1 then plotStr += '_xlog'
  
  ; plot (1) - radial profiles of all quantities, one per panel, all halos all resolutions
  start_PS, rp.(0).(0).sP.plotPath + 'zoomRadProfiles_1D.' + plotStr + '.eps', xs=12.0*1.4, ys=6.0*1.4
  
    pos = plot_pos(col=3,row=2,/gap)
    offset = [-0.03,-0.02,-0.02,0]
  
    ; (a) - gas temperature
    cgPlot,[0],[0],/nodata,xrange=xrange_rad,yrange=yrange_temp,/xs,/ys,xlog=xrange_log,/ylog,yminor=5,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl("T_{gas} [_{ }K_{ }]"),pos=pos[0]+offset,/noerase
    
    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrange_temp+[1e4,1e6],line=1,color='light gray',/overplot
    
    for i=0,n_tags(rp)-1 do begin
      for j=0,n_tags(rp.(i))-1 do begin
        xx = rp.(i).(j).gas.rBinCen
        yy = rp.(i).(j).gas.temp
        co = rp.(i).(j).sP.colors[cInds[j]]
                
        cgPlot,xx,smooth(yy[*,0],sK),line=lines[j],color=co,/overplot
      endfor ;n_tags(rp.(i)),j
      
      ; highest resolution, quartile band
      ;oplotBand,xx,yy[*,1],yy[*,3],color='light gray',yrange=yrange_temp
    endfor ;n_tags(rp),i
        
    legend,rNames,linestyle=rLines,/top,/right
    
    ; (b) - gas dens
    cgPlot,[0],[0],/nodata,xrange=xrange_rad,yrange=yrange_dens,/xs,/ys,xlog=xrange_log,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl("log \rho_{gas} [_{ }g/cm^{3 }]"),pos=pos[1]+offset,/noerase
      
    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrange_dens+[0.05,-0.05],line=1,color='light gray',/overplot
      
    for i=0,n_tags(rp)-1 do begin
      for j=0,n_tags(rp.(i))-1 do begin
        xx = rp.(i).(j).gas.rBinCen
        yy = alog10( rp.(i).(j).gas.dens )
        co = rp.(i).(j).sP.colors[cInds[j]]
                
        cgPlot,xx,smooth(yy[*,0],sK),line=lines[j],color=co,/overplot
      endfor ;n_tags(rp.(i)),j
    endfor ;n_tags(rp),i
    
    legend,hNames,textcolor=hColors,/top,/right
      
    ; (c) - gas vrad
    cgPlot,[0],[0],/nodata,xrange=xrange_rad,yrange=yrange_vrad,/xs,/ys,xlog=xrange_log,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl("v_{rad} [_{ }km/s_{ }]"),pos=pos[2]+offset,/noerase
      
    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrange_vrad+[0.05,-0.05],line=1,color='light gray',/overplot
      
    for i=0,n_tags(rp)-1 do begin
      for j=0,n_tags(rp.(i))-1 do begin
        xx = rp.(i).(j).gas.rBinCen
        yy = rp.(i).(j).gas.vrad
        co = rp.(i).(j).sP.colors[cInds[j]]
                
        cgPlot,xx,smooth(yy[*,0],sK),line=lines[j],color=co,/overplot
      endfor ;n_tags(rp.(i)),j
    endfor ;n_tags(rp),i
    
    ; (d) - gas cell size
    cgPlot,[0],[0],/nodata,xrange=xrange_rad,yrange=yrange_csize,/xs,/ys,xlog=xrange_log,/ylog,yminor=5,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl("r_{cell} [_{ }kpc_{ }]"),pos=pos[3]+offset,/noerase
      
    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrange_csize+[0.05,-0.05],line=1,color='light gray',/overplot
      
    for i=0,n_tags(rp)-1 do begin
      for j=0,n_tags(rp.(i))-1 do begin
        xx = rp.(i).(j).gas.rBinCen
        yy = rp.(i).(j).gas.csize
        co = rp.(i).(j).sP.colors[cInds[j]]
                
        cgPlot,xx,smooth(yy[*,0],sK),line=lines[j],color=co,/overplot
      endfor ;n_tags(rp.(i)),j
    endfor ;n_tags(rp),i
    
    ; (e) - entropy
    cgPlot,[0],[0],/nodata,xrange=xrange_rad,yrange=yrange_entr,/xs,/ys,xlog=xrange_log,/ylog,yminor=5,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl("log S_{gas} [_{ }K cm^{2 }]"),pos=pos[4]+offset,/noerase
      
    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrange_entr+[0.05,-0.05],line=1,color='light gray',/overplot
      
    for i=0,n_tags(rp)-1 do begin
      for j=0,n_tags(rp.(i))-1 do begin
        xx = rp.(i).(j).gas.rBinCen
        yy = rp.(i).(j).gas.entr
        co = rp.(i).(j).sP.colors[cInds[j]]
                
        cgPlot,xx,smooth(yy[*,0],sK),line=lines[j],color=co,/overplot
      endfor ;n_tags(rp.(i)),j
    endfor ;n_tags(rp),i
    
    ; (f) - ang mom
    cgPlot,[0],[0],/nodata,xrange=xrange_rad,yrange=yrange_angm,/xs,/ys,xlog=xrange_log,/ylog,yminor=5,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl("J_{gas} [_{ }kpc km/s_{ }]"),pos=pos[5]+offset,/noerase
      
    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrange_angm+[0.05,-0.05],line=1,color='light gray',/overplot
      
    for i=0,n_tags(rp)-1 do begin
      for j=0,n_tags(rp.(i))-1 do begin
        xx = rp.(i).(j).gas.rBinCen
        yy = rp.(i).(j).gas.angm
        co = rp.(i).(j).sP.colors[cInds[j]]
                
        cgPlot,xx,smooth(yy[*,0],sK),line=lines[j],color=co,/overplot
      endfor ;n_tags(rp.(i)),j
    endfor ;n_tags(rp),i
      
  end_PS

  stop

end

; plotZoomRadialRays(): consider individual radial rays evenly distributed in (theta,phi)

pro plotZoomRadialRays
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  hInds     = [0]
  resLevels = [9]
  redshift  = 2.0
  newSaves  = 0 ; override existing saves
  
  ; binning config
  Nside   = 8 ; healpix resolution parameter, 8=768, 16~3k, 32~12k, 64~50k, 128~200k
  radFacs = linspace(0.01,2.0,100) ; r/rvir, sampling points
  cutSubS = 2 ; remove substructures? (0=no, 1=satellites, 2=all subfind other than target group)

  ; plot config
  xrange_rad   = [0.0,2.0]   ; r/rvir
  xrange_log   = 0           ; plot in log? doesn't necessary have to equal binInLog  
  yrange_temp  = [1e4,1e7]   ; K
  yrange_dens  = [-28,-24]   ; log cgs (g/cm^3)
  yrange_vrad  = [-350,250]  ; km/s
  yrange_csize = [0.01,2.0]  ; physical kpc
  yrange_entr  = [1e7,1e9]   ; cgs (K cm^2)
  yrange_angm  = [1e3,1e5]   ; kpc km/s
  
  lines    = [0,1,2] ;[1,2,0] ; line style, one per resLevel
  cInds    = [1,1,1] ; color index, one per resLevel
  radLines = [0.15,1.5]
  sK       = 1
  thick    = 0.1
  
  ; load
  foreach hInd,hInds do begin
    hLoc = {}
        
    foreach resLevel,resLevels do begin
      ; get list of snapshots for this halo in the redshift range
      rLoc = {}
      
      if resLevel eq 11 and total(hInd eq [2,3,4,5,6,9]) gt 0 then continue
      
      sP = simParams(run='zoom_20Mpc',res=resLevel,hInd=hInd,redshift=redshift)     
      gcInd = zoomTargetHalo(sP=sP)
      
      ; calculate radial rays using nested healpix spheres
      hsv_temp = haloShellValue(sP=sP, partType='gas', valName='temp', subgroupIDs=[gcInd], $
                                Nside=Nside, radFacs=radFacs, cutSubS=cutSubS, newSaves=newSaves)
      hsv_temp.value = 10.0^hsv_temp.value ; remove log
      
      hsv_vrad = haloShellValue(sP=sP, partType='gas', valName='radvel', subgroupIDs=[gcInd], $
                                Nside=Nside, radFacs=radFacs, cutSubS=cutSubS, newSaves=newSaves)
              
      rLoc = mod_struct( rLoc, 'sP', sP )
      rLoc = mod_struct( rLoc, 'hsv_temp', hsv_temp )
      rLoc = mod_struct( rLoc, 'hsv_vrad', hsv_vrad )
    
      ; add this resolution level to the halo struct
      hLoc = mod_struct( hLoc, 'L'+str(resLevel), rLoc )     
     
    endforeach ;resLevels
    
    ; add this halo to the save struct
    rr = mod_struct( rr, 'h'+str(hInd), hLoc )
    
  endforeach ;hInds
  
  ; plot setup
  plotStr = 'h'
  hNames  = 'h' + str(hInds)
  hColors = []
  rNames  = 'L' + str(resLevels)
  rLines  = lines[0:n_elements(resLevels)-1]

  foreach hInd,hInds do plotStr += str(hInd)
  foreach hInd,hInds,i do hColors = [hColors,rr.(i).(0).sP.colors[cInds[-1]]]
  plotStr += '_L'
  foreach resLevel,resLevels do plotStr += str(resLevel)
  plotStr += '_cutSubS-' + str(cutSubS)

  if xrange_log eq 1 then xrange_rad[0] = 0.01
  if xrange_log eq 1 then plotStr += '_xlog'
  
  ; establish color mapping for temp
  set_plot,'ps'
  mbColors = reverse( sampleColorTable('blue-red2', 256, bounds=[0.0,1.0]) )
  loadColorTable,'bw linear'
  
  pos = [0.12,0.12,0.87,0.9]
  pos_cbar = [0.89,0.12,0.92,0.9]
  
  ; plot (1) - temperature individual rays
  start_PS, rr.(0).(0).sP.plotPath + 'radialRays_TempIndiv.' + plotStr + '.eps' ;, xs=12.0*1.5, ys=6.0*1.5
  
    ; (a) - gas temperature
    cgPlot,[0],[0],/nodata,xrange=xrange_rad,yrange=yrange_temp,/xs,/ys,xlog=xrange_log,/ylog,yminor=5,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl("T_{gas} [_{ }K_{ }]"),/noerase,pos=pos

    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrange_temp+[0.05,-0.05],line=1,color='light gray',/overplot
      
    for i=0,n_tags(rr)-1 do begin
      for j=0,n_tags(rr.(i))-1 do begin
        xx = rr.(i).(j).hsv_temp.radFacs
        co = rr.(i).(j).sP.colors[cInds[j]]
        
        ; plot each radial ray, corresponding to one healpix pixel
        for k=0,rr.(i).(j).hsv_temp.nPx-1 do begin
          yy = reform( rr.(i).(j).hsv_temp.value[k,*] )
          cgPlot,xx,smooth(yy[*,0],sK),line=lines[j],color=co,/overplot,thick=thick
        endfor
        
      endfor ;n_tags(rr.(i)),j
    endfor ;n_tags(rr),i
    
    legend,hNames,textcolor=hColors,/top,/right
    legend,rNames,linestyle=rLines,/bottom,/left
      
  end_PS
  
  ; plot (1b) - temperature individual rays, colored by vrad
  start_PS, rr.(0).(0).sP.plotPath + 'radialRays_TempIndivColor.' + plotStr + '.eps' ;, xs=12.0*1.5, ys=6.0*1.5
   
    ; (a) - gas temperature
    cgPlot,[0],[0],/nodata,xrange=xrange_rad,yrange=yrange_temp,/xs,/ys,xlog=xrange_log,/ylog,yminor=5,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl("T_{gas} [_{ }K_{ }]"),/noerase,pos=pos

    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrange_temp+[0.05,-0.05],line=1,color='light gray',/overplot
      
    for i=0,n_tags(rr)-1 do begin
      for j=0,n_tags(rr.(i))-1 do begin
        xx = rr.(i).(j).hsv_temp.radFacs
        co_mm = yrange_vrad
        co_vv = reform( rr.(i).(j).hsv_vrad.value )
        co_vv = bytscl( co_vv, min=co_mm[0], max=co_mm[1] )
        
        ; plot each radial ray, corresponding to one healpix pixel
        for k=0,rr.(i).(j).hsv_temp.nPx-1 do begin
          yy = reform( rr.(i).(j).hsv_temp.value[k,*] )
          
          ; use temperature at each radial sample to color segment, map into [0,255]
          co = mbColors[co_vv[k,*]]
          
          for m=0,n_elements(yy)-2 do begin
            xx0 = xx[m]
            xx1 = xx[m+1]
            yy0 = yy[m]
            yy1 = yy[m+1]
            cgPlot,[xx0,xx1],[yy0,yy1],line=lines[j],color=co[m],/overplot,thick=thick
          endfor ;m
        
        endfor ;k
      endfor ;n_tags(rr.(i)),j
    endfor ;n_tags(rr),i
    
    legend,hNames,textcolor=hColors,/top,/right
    legend,rNames,linestyle=rLines,/bottom,/left
      
    loadColorTable,'blue-red2'
    cgColorbar,pos=pos_cbar,divisions=6,/vertical,/right,range=co_mm,title=textoidl("v_{rad} [_{ }km/s_{ }]")
      
  end_PS
  
  ; plot (2) - vrad individual rays, colored by temperature
  start_PS, rr.(0).(0).sP.plotPath + 'radialRays_VRadIndivColor.' + plotStr + '.eps' ;, xs=12.0*1.5, ys=6.0*1.5
  
    ; (a) - gas radial velocity
    cgPlot,[0],[0],/nodata,xrange=xrange_rad,yrange=yrange_vrad,/xs,/ys,xlog=xrange_log,$
      xtitle=textoidl("r / r_{vir}"),ytitle=textoidl("v_{rad} [_{ }km/s_{ }]"),/noerase,pos=pos

    foreach radLine,radLines do $
      cgPlot,[radLine,radLine],yrange_vrad+[20,-20],line=1,color='light gray',/overplot
      
    for i=0,n_tags(rr)-1 do begin
      for j=0,n_tags(rr.(i))-1 do begin
        xx = rr.(i).(j).hsv_vrad.radFacs
        co_mm = alog10( yrange_temp )
        co_vv = alog10( reform( rr.(i).(j).hsv_temp.value ) )
        co_vv = bytscl( co_vv, min=co_mm[0], max=co_mm[1] )
        
        if ~array_equal( size(rr.(i).(j).hsv_vrad.value), size(rr.(i).(j).hsv_temp.value) ) then message,'Error'
        
        ; plot each radial ray, corresponding to one healpix pixel
        for k=0,rr.(i).(j).hsv_temp.nPx-1 do begin
          yy = reform( rr.(i).(j).hsv_vrad.value[k,*] )
          
          ; use temperature at each radial sample to color segment, map into [0,255]
          co = mbColors[co_vv[k,*]]
          
          for m=0,n_elements(yy)-2 do begin
            xx0 = xx[m]
            xx1 = xx[m+1]
            yy0 = yy[m]
            yy1 = yy[m+1]
            cgPlot,[xx0,xx1],[yy0,yy1],line=lines[j],color=co[m],/overplot,thick=thick
          endfor ;m
          
        endfor ;k
      endfor ;n_tags(rr.(i)),j
    endfor ;n_tags(rr),i
    
    legend,hNames,textcolor=hColors,/top,/right
    legend,rNames,linestyle=rLines,/top,/left
      
    loadColorTable,'blue-red2'
    cgColorbar,pos=pos_cbar,divisions=6,/vertical,/right,range=co_mm,title=textoidl("T_{gas} [_{ }log K_{ }]")
      
  end_PS
  stop

end