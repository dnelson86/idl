; zoomEvo.pro
; 'zoom project' evolution of halo properties with redshift
; dnelson dec.2014

; plotZoomMassEvo(): plot time evolution of various masses

pro plotZoomMassEvo
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  hInds     = [0,1,5,7,8,9] ; TODO: h2L9 and h3L10 cannot make galcat @ snap59, maybe stop at snap 58?)
  resLevels = [9,10,11]
  newSaves  = 0 ; override existing saves
  
  ; plot config
  xrange = [2.0,6.0]    ; redshift
  xtitle = "Redshift"
  xlog   = 1
  xtickv = [2,3,4,5,6]
  
  ytitle = { mass_tot   : "M_{sub,total} [_{ }log M_{sun }]" ,$ ; MUST align with vals struct
             mass_gas   : "M_{sub,gas} [_{ }log M_{sun }]"   ,$
             mass_stars : "M_{sub,stars} [_{ }log M_{sun }]" ,$
             mass_gas2  : "M_{gas,<r} [_{ }log M_{sun }]"    ,$
             halo_rvir  : "R_{vir} [_{ }ckpc_{ }]"           ,$
             halo_tvir  : "T_{vir} [_{ }log K_{ }]"           }
  yrange = { mass_tot   : [10.0,12.1] ,$ ; log msun, total halo mass
             mass_gas   : [9.5,10.8]  ,$ ; log msun, total subhalo gas mass
             mass_stars : [9.0,11.5]  ,$ ; log msun, total subhalo stellar mass
             mass_gas2  : [9.5,10.8]  ,$ ; log msun, gas mass within some radius
             halo_rvir  : [100,350]   ,$ ; ckpc
             halo_tvir  : [5.5,6.201]  } ; log K
  ylog   = { mass_tot     : 1 ,$ ; are values in log
             mass_gas     : 1 ,$
             mass_stars   : 1 ,$
             mass_gas2    : 1 ,$
             halo_rvir    : 0 ,$
             halo_tvir    : 1  }
  
  psym    = 8       ; plot symbol for ending value
  symsize = 0.9     ; symbol size for ending value
  
  yrangeSub = [0.1,10.0]
  ytickvSub = [0.2,1.0,5.0]
  ytitleSub = "Ln/L11" ;"Ratio"
  
  ; A (resolutions are different linestyles, all same color and thickness)
  ;lines    = [1,2,0] ; line style, one per resLevel
  ;cInds    = [1,1,1] ; color index, one per resLevel
  ;thick    = [1,1,1] ; times !p.thick

  ; B (resolutions are all solid lines, fainter and thinner for lower res)
  lines   = [0,0,0]
  thick   = [0.25,0.5,1.0]
  cInds   = [2,1,1]
  
  ; load
  foreach hInd,hInds do begin
    hLoc = {}
        
    foreach resLevel,resLevels do begin
      ; get list of snapshots for this halo in the redshift range
      rLoc = {}
      
      ;if resLevel eq 11 and total(hInd eq [2,3,4,5,6,9]) gt 0 then continue
      
      sP = simParams(run='zoom_20Mpc',res=resLevel,hInd=hInd,redshift=2.0)     
      saveFilename = sP.derivPath + 'binnedVals/massEvo.' + sP.savPrefix + str(sP.res) + '.sav'
      
      if ~file_test(saveFilename) or newSaves eq 1 then begin
                  
        ; halo tracking
        tarInd = zoomTargetHalo(sP=sP)
        mts    = mergerTreeSubset(sP=sP,/verbose)
        nSnaps = n_elements(mts.snaps)
        
        ; allocate arrays
        vals = { mass_tot      : fltarr(nSnaps) ,$ ; primary subfind, all types
                 mass_gas      : fltarr(nSnaps) ,$ ; primary subfind, gas
                 mass_stars    : fltarr(nSnaps) ,$ ; primary subfind, stars
                 mass_halo_gas : fltarr(nSnaps) ,$ ; primary subfind+radial cuts, gas
                 rVir          : fltarr(nSnaps) ,$ ; virial radius, comoving kpc
                 tVir          : fltarr(nSnaps)  } ; virial temperature, log K
                
        if max(mts.hInd[*,tarInd]) gt 0 then print,'Warning: Target halo becomes not most massive in the past.'
        if mts.hMinSnap[tarInd] ne -1 then message,'Error: Target halo not tracked all the way.'
        
        ; loop over snaps, calculate mass in each, save
        foreach snap,mts.snaps,m do begin
          ; load group catalog
          sP.snap = snap
          sP.redshift = 1/mts.times[m]-1
          
          gc = loadGroupCat(sP=sP,/skipIDs)
          gcInd = mts.hInd[m,tarInd] ; for this snapshot
          
          vals.rVir[m] = gc.group_r_crit200[gc.subgroupGrNr[gcInd]]
          vals.tVir[m] = codeMassToVirTemp( gc.subgroupMass[gcInd], redshift=sP.redshift, /log )
          
          print,' h'+str(hInd)+'L'+str(resLevel)+': ['+string(snap,format='(I2)')+'] gcInd = ' + str(gcInd)
          
          vals.mass_tot[m]       = gc.subgroupMass[gcInd]
          vals.mass_gas[m]   = gc.subgroupMassType[partTypeNum('gas'),gcInd]
          vals.mass_stars[m] = gc.subgroupMassType[partTypeNum('stars'),gcInd]
          
          ; load positions of gas, make radial cut for halo selection
          ; TODO
          
          ; calculate masses of halo gas, based on temp/vrad cuts
          vals.mass_halo_gas[m] = 0.0 ; TODO
          
        endforeach ;snaps
        
        ; units and little h factors
        vals.mass_tot      = codeMassToLogMsun( vals.mass_tot / units.hubbleParam )
        vals.mass_gas      = codeMassToLogMsun( vals.mass_gas / units.hubbleParam )
        vals.mass_stars    = codeMassToLogMsun( vals.mass_stars / units.hubbleParam )
        vals.mass_halo_gas = codeMassToLogMsun( vals.mass_halo_gas / units.hubbleParam )
        
        vals.rVir /= units.hubbleParam ; still comoving
        
        ; save all values for this halo+resLevel combination
        rLoc = mod_struct( rLoc, 'sP', sP )
        rLoc = mod_struct( rLoc, 'redshifts', 1/mts.times-1 )
        rLoc = mod_struct( rLoc, 'vals', vals )
       
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
  start_PS, evo.(0).(0).sP.plotPath + 'zoomMassEvo_total_' + plotStr + '.eps'
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange.mass_tot,/xs,/ys,$
      xlog=xlog,xminor=2,xtickv=xtickv,xticks=n_elements(xtickv)-1,$
      xtitle=xtitle,ytitle=textoidl(ytitle.mass_tot)
    
    for i=0,n_tags(evo)-1 do begin
      for j=0,n_tags(evo.(i))-1 do begin
        xx = evo.(i).(j).redshifts
        yy = evo.(i).(j).vals.mass_tot
        co = evo.(i).(j).sP.colors[cInds[j]]
        
        cgPlot,xx,yy,line=lines[j],color=co,thick=!p.thick*thick[j],/overplot
      endfor ;n_tags(evo.(i)),j
      
      ; endpoint symbol for highest resolution
      cgPlot,xx[0],yy[0],psym=psym,symsize=symsize,color=co,/overplot
    endfor ;n_tags(evo),i
    
    legend,hNames,textcolor=hColors,/bottom,/left
    legend,rNames,linestyle=rLines,/top,/right
      
  end_PS
  
  ; plot (2) - five mass components, vs redshift (3x2 panels)
  start_PS, evo.(0).(0).sP.plotPath + 'zoomMassEvo_3x2_' + plotStr + '.eps', xs=12.0, ys=6.0
  
    pos = plot_pos(col=3,row=2,/gap)
    offset = [0,0,-0.02,0]
    
    plotTags = ['MASS_TOT','MASS_GAS','MASS_STARS','HALO_RVIR','HALO_TVIR']
    
    ; loop over each mass component (separate panel for each)
    for m=0,n_elements(plotTags)-1 do begin
      ; find index of this mass component
      count = where( tag_names(ytitle) eq plotTags[m] )
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange.(count),/xs,/ys,$
        xlog=xlog,xminor=2,xtickv=xtickv,xticks=n_elements(xtickv)-1,$
        xtitle=xtitle,ytitle=textoidl(ytitle.(count)),pos=pos[m]+offset,/noerase
      
      for i=0,n_tags(evo)-1 do begin
        for j=0,n_tags(evo.(i))-1 do begin
          xx = evo.(i).(j).redshifts
          yy = evo.(i).(j).vals.(count)          
          co = evo.(i).(j).sP.colors[cInds[j]]
          
          cgPlot,xx,yy,line=lines[j],color=co,thick=!p.thick*thick[j],/overplot
        endfor ;n_tags(evo.(i)),j
        
        ; endpoint symbol for highest resolution
        cgPlot,xx[0],yy[0],psym=psym,symsize=symsize,color=co,/overplot
      endfor ;n_tags(evo),i
      
      if m eq 0 then legend,hNames,textcolor=hColors,/bottom,/left
      if m eq 1 then legend,rNames,linestyle=rLines,/top,/right
    endfor
      
  end_PS
  
  ; plot (3) - three mass components and virial temp (2x2 panels with subpanels for resolution ratios)
  start_PS, evo.(0).(0).sP.plotPath + 'zoomMassEvoSubPanels_' + plotStr + '.eps', xs=8.0*1.4, ys=6.0*1.4
  
    pos = plot_pos(col=2,row=2,/gap)
    subSize = 0.1
    offset = [-0.03,-0.04+subSize,0.02,0]
      
    plotTags = ['MASS_TOT','MASS_GAS','MASS_STARS','HALO_TVIR']
      
    ; (i) - main panel
    for m=0,n_elements(plotTags)-1 do begin
      ; find index of this mass component
      count = where( tag_names(ytitle) eq plotTags[m] )
      print,plotTags[m],count
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange.(count),/xs,/ys,$
        xlog=xlog,xminor=2,xtickv=xtickv,xticks=n_elements(xtickv)-1,$
        xtitle="",xtickname=replicate(" ",10),ytitle=textoidl(ytitle.(count)),$
        pos=pos[m]+offset,/noerase
      
      for i=0,n_tags(evo)-1 do begin
        j = where( tag_names(evo.(i)) eq 'L11', count_ind )
        if count_ind eq 0 then continue
        ;for j=0,n_tags(evo.(i))-1 do begin
          xx = evo.(i).(j).redshifts
          yy = evo.(i).(j).vals.(count)          
          co = evo.(i).(j).sP.colors[cInds[j]]
          
          cgPlot,xx,yy,line=lines[j],color=co,thick=!p.thick*thick[j],/overplot
        ;endfor ;n_tags(evo.(i)),j
        
        ; endpoint symbol for highest resolution
        cgPlot,xx[0],yy[0],psym=psym,symsize=symsize,color=co,/overplot
      endfor ;n_tags(evo),i
      
      if m eq 0 then legend,hNames,textcolor=hColors,/bottom,/left
      if m eq 1 then legend,[rNames[0:1]+'/'+rNames[2],rNames[2]],$
        linestyle=rLines,thick=thick*!p.thick,/top,/right
      
      ; (ii) - subpanel
      posSub = pos[m]+offset
      posSub[1] -= subSize ; move y of lowerleft down
      posSub[3] = posSub[1] + subSize ; set height of subpanel
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeSub,/xs,/ys,ylog=1,yminor=1,$
        xlog=xlog,xminor=2,xtickv=xtickv,xticks=n_elements(xtickv)-1,$
        xtitle=textoidl(xtitle),ytitle=textoidl(ytitleSub),ytickv=ytickvSub,yticks=2,pos=posSub,/noerase
      
      cgPlot,xrange,[1.0,1.0],line=0,color=cgColor('light gray'),/overplot
      
      for i=0,n_tags(evo)-1 do begin
        ind_L11 = where( tag_names(evo.(i)) eq 'L11', count_ind )
        if count_ind eq 0 then continue
        
        yy_L11 = evo.(i).(ind_L11).vals.(count)
        if ylog.(count) eq 1 then yy_L11 = 10.0^yy_L11 ; remove log
        
        for j=0,n_tags(evo.(i))-2 do begin
          xx = evo.(i).(j).redshifts
          
          yy = evo.(i).(j).vals.(count)      
          if ylog.(count) eq 1 then yy = 10.0^yy ; remove log
          yy /= yy_L11 ; resolution ratio
          
          co = evo.(i).(j).sP.colors[cInds[j]]
          
          cgPlot,xx,yy,line=lines[j],color=co,thick=!p.thick*thick[j],/overplot
        endfor ;n_tags(evo.(i)),j
        
        ; endpoint symbol for highest resolution
        cgPlot,xx[0],yy[0],psym=psym,symsize=symsize,color=co,/overplot
      endfor ;n_tags(evo),i
      
    endfor ;m
    
  end_PS
  
  stop

end
