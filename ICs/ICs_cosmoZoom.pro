; ICs_cosmoZoom.pro
; zoom project - helper for generating zoom ICs with MUSIC (O. Hahn)
; dnelson sep.2013

; checkBHDerefGal(): debugging

pro checkBHDerefGal

  ; how similar are the snapshots?
  snapNum = 0
  sP_old = simParams(run='zoom_20Mpc',res=9,hInd=0,snap=snapNum)
  sP_new = simParams(run='zoom_20Mpc_derefgal_nomod',res=9,hInd=0,snap=snapNum)
  
  h_old = loadSnapshotHeader(sP=sP_old)
  h_new = loadSnapshotHeader(sP=sP_new)
  
  print,h_old.nPartTot
  print,h_new.nPartTot
  ; random number differences due to BHs being there? (new has more gas cells, less stars)
  
  stop
  
end

pro checkZoomDerefineGal

  ; config
  sP = simParams(run='zoom_20Mpc_derefgal',res=9,hInd=0,redshift=2.2)
  sP.snap = 48
  fileName = '/n/home07/dnelson/sims.zooms/128_20Mpc_h0_L9_derefgal/output/bhderef.txt'
  ptStruct = {line:''}
  headerLines = 1
  
  outer_mass_fac  = 1.0
  inner_mass_fac  = 64.0
  outer_deref_rad = 10.0
  colors = ['red','blue','orange']
  
  ; load
  x = loadCSV(headerLines, fileName, ptStruct)
  x = x.line
  
  ; count proper lines (some are malformed)
  count = 0L
  
  for i=0,n_elements(x)-1 do begin
    startPos = strpos(x[i], ':', 0) + 2
    line = strmid(x[i], startPos)
    line = strsplit(line,' ',/extract)
    if n_elements(line) ne 13 and n_elements(line) ne 14 then continue
    
    count += 1
  endfor
  
  print,n_elements(x),count,float(count)/n_elements(x)*100.0
  
  createFlag = intarr(count) ; 1=create at this time, 0=reposition
  time       = fltarr(count)
  bh_id      = lonarr(count)
  timebin    = intarr(count)
  minpot     = fltarr(count)
  redshift   = fltarr(count)
  xyz        = fltarr(3,count)
  
  ; parse
  count = 0L
  
  for i=0,n_elements(x)-1 do begin
    ; cut off task number and BH_DEREFINE_GAL:
    startPos = strpos(x[i], ':', 0) + 2
    line = strmid(x[i], startPos)
    
    line = strsplit(line,' ',/extract)
    if n_elements(line) ne 13 and n_elements(line) ne 14 then continue
    
    if line[0] eq 'BH_Create:' then begin
      ; create
      createFlag[count] = 1
      time[count]     = float(line[2])
      bh_id[count]    = long(line[5])
      timebin[count]  = 0
      minpot[count]   = 0.0
      xyz[0,count]    = line[10]
      xyz[1,count]    = line[11]
      xyz[2,count]    = line[12]
    endif else begin
      ; reposition
      time[count]     = float(line[2])
      bh_id[count]    = long(line[5])
      timebin[count]  = fix(line[7])
      minpot[count]   = float(line[9])
      xyz[0,count]    = line[11]
      xyz[1,count]    = line[12]
      xyz[2,count]    = line[13]
    endelse
    
    count += 1
  endfor
  
  ; stats
  redshift = 1.0/time-1.0
  
  w = where(createFlag eq 0,comp=wc)
  print,'Found ['+str(nuniq(bh_id[w]))+'] unique BH IDs.'
  print,'Time: ',minmax(time)
  print,'Timebins: ',minmax(timebin[w])
  print,'XYZ: ',minmax(xyz)
  print,'Minpot: ',minmax(minpot[w])
  
  ids_repo = bh_id[w]
  uniq_ids = ids_repo[uniq(ids_repo,sort(ids_repo))]
  
  ; plot
  xrange = minmax(time)*[0.99,1.01]
  
  start_PS,sP.plotPath+'derefGal_pos_'+sP.saveTag+'.eps', /big
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=minmax(xyz)*[0.99,1.01]/1e4,/xs,/ys,$
      ytitle="XYZ Coordinate",xtitle="Scale Factor"
      
    for i=0,2 do begin
      cgPlot,time,xyz[i,*]/1e4,psym=3,color=cgColor(colors[i]),/overplot
      cgPlot,time[wc],xyz[i,wc]/1e4,psym=4,color=cgColor(colors[i]),/overplot
    endfor
    
    legend,['x','y','z'],textcolor=colors,/top,/left
  end_PS
  
  start_PS,sP.plotPath+'derefGal_minpot_'+sP.saveTag+'.eps', /big
    w = where(minpot ne 0.0)
    yrange = minmax(minpot[w])*[1.01,0.99]/1e4
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Minimum Potential Value",xtitle="Scale Factor"
      
    cgPlot,time,minpot/1e4,psym=3,color=cgColor('blue'),/overplot
  end_PS
  
  ; load group catalog, find most massive center
  gc = loadGroupCat(sP=sP,/skipIDs)
  
  hInd = 0
  rmax = 20.0 ; ckpc
  
  grNr = gc.subgroupGrNr[hInd]
  groupPos = gc.groupPos[*,grNr]
  print,'Fof GrNr: '+str(grNr)+' Pos ['+$
    string(groupPos[0],format='(f7.1)')+' '+$
    string(groupPos[1],format='(f7.1)')+' '+$
    string(groupPos[2],format='(f7.1)')+']'
  print,'SF hInd: '+str(hInd)+' Pos ['+$
    string(gc.subgroupPos[0,hInd],format='(f7.1)')+' '+$
    string(gc.subgroupPos[0,hInd],format='(f7.1)')+' '+$
    string(gc.subgroupPos[0,hInd],format='(f7.1)')+']'
  print,'SF hInd: '+str(hInd)+' CM ['+$
    string(gc.subgroupCM[0,hInd],format='(f7.1)')+' '+$
    string(gc.subgroupCM[0,hInd],format='(f7.1)')+' '+$
    string(gc.subgroupCM[0,hInd],format='(f7.1)')+']'
    
  ; load BH positions
  h = loadSnapshotHeader(sP=sP)
  numBHs = h.nPartTot[partTypeNum('bh')]
  
  pos_bh = loadSnapshotSubset(sP=sP,partType='bh',field='pos')
  nuniq_pos = nuniq( pos_bh[0,*]+pos_bh[1,*]+pos_bh[2,*] )
  
  print,'In snapshot z='+string(snapNumToRedshift(sP=sP),format='(f4.2)')+' have ['+str(numBHs)+'] BHs with ['+$
    str(nuniq_pos)+'] unique positions (snap='+str(sP.snap)+')'
  print,pos_bh
  
  start_PS,sP.plotPath+'derefGal_bhPos_'+sP.saveTag+'.eps', xs=8.0, ys=8.0
    xyrange = minmax(pos_bh)*[0.99,1.01]/1e4
    cgPlot,[0],[0],/nodata,xrange=xyrange,yrange=xyrange,xtitle="X",ytitle="Y",/xs,/ys,aspect=1.0
    
    ; most massive fof/subfind position
    cgPlot,groupPos[0]/1e4,groupPos[1]/1e4,psym='filled circle',$
      symsize=2.0,color=cgColor('red'),/overplot
    cgPlot,gc.subgroupPos[0]/1e4,gc.subgroupPos[1]/1e4,psym='filled circle',$
      symsize=2.0,color=cgColor('blue'),/overplot
    cgPlot,gc.subgroupCM[0]/1e4,gc.subgroupCM[1]/1e4,psym='filled circle',$
      symsize=2.0,color=cgColor('orange'),/overplot
    
    ; bh positions
    cgPlot,pos_bh[0,*]/1e4,pos_bh[1,*]/1e4,psym=4,color=cgColor('black'),/overplot
  end_PS
    
  ; load gas positions (make into radial distances) and masses
  pos  = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  mass = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
  ids  = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
  
  ; histogram of gas cell masses relative to target
  start_PS,sP.plotPath+'derefGal_massHisto_'+sP.saveTag+'.eps'
    plothist,alog10(mass/sP.targetGasMass),/auto,xrange=[-1,3.3],$
      xtitle="log (Gas Cell Mass / targetGasMass)",ytitle="N",/ylog,yminor=0
  end_PS
  
  ; scatterplot of >100 target mass in (x,y) plane
  start_PS,sP.plotPath+'derefGal_highMassPos_'+sP.saveTag+'.eps', xs=8, ys=8
    zoomFac = 5.0
    axes = [0,1]
    
    xrange = groupPos[axes[0]] + [-sP.boxSize/zoomFac,sP.boxSize/zoomFac]
    yrange = groupPos[axes[1]] + [-sP.boxSize/zoomFac,sP.boxSize/zoomFac]
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      xtitle="axis"+str(axes[0]),ytitle="axis"+str(axes[1]),aspect=1.0
    
    xx = mass/sP.targetGasMass
    cuts = [5, 10, 20]
    
    w = where(xx ge cuts[0] and xx le cuts[1] and ids le sP.ids_offset)
    cgPlot,pos[axes[0],w],pos[axes[1],w],psym=4,color=cgColor('blue'),/overplot
    w = where(xx ge cuts[1] and xx le cuts[2] and ids le sP.ids_offset)
    cgPlot,pos[axes[0],w],pos[axes[1],w],psym=4,color=cgColor('red'),/overplot
    w = where(xx ge cuts[2] and ids le sP.ids_offset)
    cgPlot,pos[axes[0],w],pos[axes[1],w],psym=4,color=cgColor('orange'),/overplot
    
    cgPlot,groupPos[axes[0]],groupPos[axes[1]],psym='open circle',symsize=4.0,/overplot
  end_PS
  
  ; calculate radial distances
  dists = periodicDists(groupPos,pos,sP=sP)
  
  w = where(dists le rmax,count)
  if count eq 0 then message,'Error'
  
  dists = dists[w]
  mass  = mass[w]
  
  ; plot expected mass(r) profile for gas cells vs scatter of actual
  start_PS,sP.plotPath+'derefGal_massRadProfile_'+sP.saveTag+'.eps'
    cgPlot,[0],[0],/nodata,xrange=[0,rmax],yrange=[0.3,300.0],/xs,/ys,/ylog,yminor=0,$
      xtitle="Radial Distance [ckpc]",ytitle="Gas Cell Mass / targetGasMass"
      
    !p.thick -= 1
    cgPlot,dists,mass/sP.targetGasMass,psym=4,/overplot
    
    ; theory (1)
    inner_slope = 0
    r = linspace(0.0,rmax,1000)
    t = r / outer_deref_rad
    localTargetFac = outer_mass_fac * (3*t*t - 2*t*t*t) + $ ; h00
                     inner_mass_fac * (1 - 3*t*t + 2*t*t*t) + $ ;h01
                     ( (t*t*t)-2*t*t+t ) * inner_slope ; h10 (tangent at t=0)
                     
    ;localTargetFac = outer_mass_fac / (t*t) * (1.0-t)
    
    ; t>1 set ltF=1
    w = where(t gt 1.0,count)
    if count gt 0 then localTargetFac[w] = 1.0
    
    ;cgPlot,r,localTargetFac,line=0,thick=!p.thick+2,color=cgColor('blue'),/overplot
    
    ; theory (2)
    KERNEL_COEFF_1 = 2.546479089470
    KERNEL_COEFF_1_INV = 0.3926990926
    KERNEL_COEFF_2 = 15.278874536822
    KERNEL_COEFF_3 = 5.092958178941
    
    u = r / outer_deref_rad
    
    localTargetFac = fltarr(n_elements(u))
    w = where(u lt 0.5)
    localTargetFac[w] = KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u[w] - 1.0) * u[w] * u[w]
    w = where(u gt 0.5 and u le 1.0)
    localTargetFac[w] = KERNEL_COEFF_3 * (1.0 - u[w]) * (1.0 - u[w]) * (1.0 - u[w])
    w = where(u gt 1.0)
    localTargetFac[w] = 0.0
    
    localTargetFac *= (inner_mass_fac-outer_mass_fac) * KERNEL_COEFF_1_INV
    localTargetFac += outer_mass_fac
    
    cgPlot,r,localTargetFac,line=0,thick=!p.thick+2,color=cgColor('orange'),/overplot
    
    ltfHalf = 0.5 * localTargetFac
    ltfDbl  = 2.0 * localTargetFac
    
    cgPlot,r,ltfHalf,line=2,thick=!p.thick+2,color=cgColor('orange'),/overplot
    cgPlot,r,ltfDbl,line=2,thick=!p.thick+2,color=cgColor('orange'),/overplot
    
    cgPlot,[1,8.5],inner_mass_fac*[1,1],line=1,color=cgColor('orange'),/overplot
    cgText,8,inner_mass_fac*1.2,str(fix(inner_mass_fac)),color=cgColor('orange'),alignment=0.5
  end_PS
  
  ; load volumes (-> cell radii) and histogram
  hsml = loadSnapshotSubset(sP=sP,partType='gas',field='vol')
  hsml = (temporary(hsml) * 3.0 / (4*!pi))^(1.0/3.0) ;cellrad [ckpc]
  
  start_PS,sP.plotPath+'derefGal_cellSizeHisto_'+sP.saveTag+'.eps'
    plothist,alog10(hsml),/auto,xrange=[-2,3.0],$
      xtitle="log (Gas Cell Size [ckpc])",ytitle="N",/ylog,yminor=0
    legend,['minimum: '+string(min(hsml)*1000,format='(f5.1)')+' pc'],/top,/left
  end_PS
  
  stop
  
end

; zoomSelectHaloSimple(): select a single halo at z=0 based on desired mass and determine bounding region in ICs
;   by matching all IDs from the primary subgroup (subfind selection, or spatial) in one step
; NOTE: assumes snapshot zero is the ICs at z=zstart

pro zoomSelectHaloSimple

  ; config
  targetHaloMass = 12.1 ; log msun at z=target
  targetRedshift = 1.66 ; resimulate until this redshift only
  zoomLevel      = 2    ; levelmax-levelmin
  
  ; plot/output config
  colors         = ['red','blue','orange'] ;x,y,z
  binSize        = 40.0 ; kpc
  axes           = [0,1] ; for box
  posFac         = 1000.0 ; kpc->Mpc
  
  sP = simParams(res=128,run='zoom_20Mpc_dm',redshift=targetRedshift)

  ; Onorbe+ (2013) suggested R_tb / rVirFac (requires >~4000 particles in halo in fullbox
  ; we just about satisfy this, ~2500-3500 for a 10^12 halo @ 20Mpc/128 or 40Mpc/256
  rVirFac = (1.5 * zoomLevel + 1)
  print,'Using rVirFac: '+str(rVirFac)+' for zoomLevel='+str(zoomLevel)
  
  ; load z=target
  ; -------------
  gc = loadGroupCat(sP=sP,/readIDs)
  
  gcIDs    = gcIDList(gc=gc,select='pri')
  gcMasses = codeMassToLogMsun( gc.subgroupMass[ gcIDs ] )
  
  ; find closest to target mass
  w = closest(gcMasses, targetHaloMass)
  gcID = gcIDs[ w[0] ]
  gc_rVir = gc.group_r_crit200[ gc.subgroupGrNr[ gcID ] ]
  
  ; Power+ (2003) suggested highres softening
  M_vir = gc.group_m_crit200[ gc.subgroupGrNr[ gcID ] ]
  gc_ids_all = gcPIDList(gc=gc,valGCids=gcID, partType='dm')
  
  highres_part_mass = (M_vir/n_elements(gc_ids_all)) / 8.0^zoomLevel
  N_vir = M_vir / highres_part_mass
  
  soft = 4.0 * gc_rVir / sqrt( N_vir )
  print,'Power soft: ',soft
  print,'Normal soft: ',20000.0/(128*2.0^zoomlevel)/40.0
  print,'Num DM in pri subfind: '+str(n_elements(gc_ids_all))
  
  ; load z=0 dm_ids and dm_pos
  ztarget_ids = loadSnapshotSubset(sP=sP,partType='dm',field='ids')
  ztarget_pos = loadSnapshotSubset(sP=sP,partType='dm',field='pos')

  ; selection type:
  if rVirFac eq 0.0 then begin
    ; find all DM ids in parent fof group
    gc_ids = gcPIDList(gc=gc, valGCids=gcID, partType='dm')
    
    ; crossmatch ids
    idIndexMap = getIDIndexMap(ztarget_ids,minid=minid)
    ztarget_ind = idIndexMap[gc_ids-minid]
  endif else begin
    ; calculate radial distance of each DM particle to halo center
    rad = periodicDists( gc.subgroupPos[*,gcID],ztarget_pos,sP=sP)
    
    ; select indices and ids
    ztarget_ind = where(rad le rvirFac * gc_rVir, count)
    if count eq 0 then message,'Error'
    print,'Num DM in rVirFac*rVir: '+str(count)
    
    ; record all DM ids in this selection
    gc_ids = ztarget_ids[ztarget_ind]
  endelse
  
  ztarget_ids = !NULL
  
  ; position subset for particle selection (and rectangular bBox)
  ztarget_pos_halo = ztarget_pos[*,ztarget_ind]
  ztarget_boundBox = { x : minmax(ztarget_pos_halo[0,*]), $
                       y : minmax(ztarget_pos_halo[1,*]), $
                       z : minmax(ztarget_pos_halo[2,*]) }
  
  print,'Selected hInd ['+str(gcID)+'] at z=0 (#'+str(w[0])+' most massive) with ['+str(n_elements(ztarget_ind))+'] dm particles.'
  
  print,' z='+string(sP.redshift,format='(f4.1)')+' pos: [ '+$
    string(gc.subgroupPos[0,gcID],format='(f7.1)')+' , '+$
    string(gc.subgroupPos[1,gcID],format='(f7.1)')+' , '+$
    string(gc.subgroupPos[2,gcID],format='(f7.1)')+' ]'
  print,' z='+string(sP.redshift,format='(f4.1)')+' boundBox:   [ '+$
    string(ztarget_boundBox.x[0],format='(i5)')+' - '+string(ztarget_boundBox.x[1],format='(i5)')+' , '+$
    string(ztarget_boundBox.y[0],format='(i5)')+' - '+string(ztarget_boundBox.y[1],format='(i5)')+' , '+$
    string(ztarget_boundBox.z[0],format='(i5)')+' - '+string(ztarget_boundBox.z[1],format='(i5)')+' ]'
         
  ; load z=zstart dm_ids
  ; --------------------
  sP.snap = sP.snapRange[0]
  sP.redshift = snapNumToRedshift(sP=sP)
  
  zstart_ids = loadSnapshotSubset(sP=sP,partType='dm',field='ids')

  ; crossmatch
  idIndexMap = getIDIndexMap(zstart_ids,minid=minid)
  zstart_ind = idIndexMap[gc_ids-minid]
  zstart_ids = !NULL
  
  ; load z=zstart dm_pos
  zstart_pos = loadSnapshotSubset(sP=sP,partType='dm',field='pos')
  zstart_pos_halo = zstart_pos[*,zstart_ind]
  
  ; center of mass and bounding box
  massCenter = { x : mean(zstart_pos_halo[0,*]),   y : mean(zstart_pos_halo[1,*]),   z : mean(zstart_pos_halo[2,*]) }
  boundBox   = { x : minmax(zstart_pos_halo[0,*]), y : minmax(zstart_pos_halo[1,*]), z : minmax(zstart_pos_halo[2,*]) }
  
  ; box center and extent
  boxCenter = { x : mean(boundBox.x), y : mean(boundBox.y), z : mean(boundBox.z) }
  boxExtent = { x : boundBox.x[1]-boundBox.x[0], y : boundBox.y[1]-boundBox.y[0], z : boundBox.z[1]-boundBox.z[0] }  
  
  ; check if bounds wrap periodic edge? if so, cannot place highres region on this halo
  for i=0,2 do if abs(boundBox.(i)[1]-boundBox.(i)[0]) gt sP.boxSize/2 then message,'Error: Wrap'
  
  print,' z='+string(sP.redshift,format='(f4.1)')+' massCenter: [ '+$
    string(massCenter.x,format='(f7.1)')+' , '+$
    string(massCenter.y,format='(f7.1)')+' , '+$
    string(massCenter.z,format='(f7.1)')+' ]'
  print,' z='+string(sP.redshift,format='(f4.1)')+' boundBox:   [ '+$
    string(boundBox.x[0],format='(i5)')+' - '+string(boundBox.x[1],format='(i5)')+' , '+$
    string(boundBox.y[0],format='(i5)')+' - '+string(boundBox.y[1],format='(i5)')+' , '+$
    string(boundBox.z[0],format='(i5)')+' - '+string(boundBox.z[1],format='(i5)')+' ]'
  
  ; MUSIC: ref_center and ref_extent fields (normalized to 1.0)
  print,' z='+string(sP.redshift,format='(f4.1)')+' boxCenter: [ '+$
    string(boxCenter.x/sP.boxSize,format='(f7.4)')+' , '+$
    string(boxCenter.y/sP.boxSize,format='(f7.4)')+' , '+$
    string(boxCenter.z/sP.boxSize,format='(f7.4)')+' ]'
  print,' z='+string(sP.redshift,format='(f4.1)')+' boxExtent: [ '+$
    string(boxExtent.x/sP.boxSize,format='(f7.4)')+' , '+$
    string(boxExtent.y/sP.boxSize,format='(f7.4)')+' , '+$
    string(boxExtent.z/sP.boxSize,format='(f7.4)')+' ]' 
  
  ; plot z=target relative positions
  start_PS,'zoomSelectHaloSimple.relPos.ztarget.eps'
    x_rel = reform( ztarget_pos_halo[0,*]-gc.subgroupPos[0,gcID] )
    y_rel = reform( ztarget_pos_halo[1,*]-gc.subgroupPos[1,gcID] )
    z_rel = reform( ztarget_pos_halo[2,*]-gc.subgroupPos[2,gcID] )
    
    xrange = minmax( [x_rel,y_rel,z_rel] ) * 1.1
    yrange = [1, n_elements(ztarget_ind)/4.0]
    plothist, x_rel, bin=binSize, color=cgColor(colors[0]), xrange=xrange, yrange=yrange, $
      /xs, /ys, /ylog, yminor=0, xtitle="Distance [kpc]", ytitle=textoidl("N_{DM}")
    plothist, y_rel, bin=binSize, color=cgColor(colors[1]), /overplot
    plothist, z_rel, bin=binSize, color=cgColor(colors[2]), /overplot
    legend,['x','y','z'],textcolors=colors,/top,/right
  end_PS
  
  ; plot z=zstart relative positions
  start_PS,'zoomSelectHaloSimple.relPos.zstart.eps'
    binFac = round( (boundBox.x[1]-boundBox.x[0]) / (ztarget_boundBox.x[1]-ztarget_boundBox.x[0]) )

    x_rel = reform( zstart_pos_halo[0,*]-massCenter.x )
    y_rel = reform( zstart_pos_halo[1,*]-massCenter.y )
    z_rel = reform( zstart_pos_halo[2,*]-massCenter.z )
    
    xrange = minmax( [x_rel,y_rel,z_rel] ) * 1.1
    yrange = [1, n_elements(zstart_ind)/5.0]
    plothist, x_rel, bin=binSize*binFac, color=cgColor(colors[0]), xrange=xrange, yrange=yrange, $
      /xs, /ys, /ylog, yminor=0, xtitle="Distance [kpc]", ytitle=textoidl("N_{DM}")
    plothist, y_rel, bin=binSize*binfac, color=cgColor(colors[1]), /overplot
    plothist, z_rel, bin=binSize*binfac, color=cgColor(colors[2]), /overplot
    legend,['x','y','z'],textcolors=colors,/top,/right
  end_PS
  
  if 0 then begin
  ; plot z=target scatterbox
  start_PS,'zoomSelectHaloSimple.box.ztarget.eps', xs=8, ys=8
    pTitle = textoidl(str(sP.res)+"^3")+" "+str(round(sP.boxSize/posFac))+"Mpc hInd="+str(gcID)+" z="+$
             string(targetRedshift,format='(f3.1)')
    cgPlot, [0], [0], /nodata, xrange=[0,sP.boxSize/posFac], yrange=[0,sP.boxSize/posFac], /xs, /ys, $
      xtitle="axis0 [Mpc]", ytitle="axis1 [Mpc]", title=pTitle
    
    plots, ztarget_pos[axes[0],*]/posFac, ztarget_pos[axes[1],*]/posFac, psym=3 ;black
    plots, ztarget_pos[axes[0],ztarget_ind]/posFac, ztarget_pos[axes[1],ztarget_ind]/posFac, psym=3, color=cgColor('red')
    
    cgPlot, ztarget_boundBox.(axes[0])[[0,1,1,0,0]]/posFac, ztarget_boundBox.(axes[1])[[0,0,1,1,0]]/posFac, $
      line=0, color=cgColor('orange'), /overplot
  end_PS, pngResize=40, /deletePS
  
  ; plot z=zstart scatterbox
  start_PS,'zoomSelectHaloSimple.box.zstart.eps', xs=8, ys=8
    pTitle = textoidl(str(sP.res)+"^3")+" "+str(round(sP.boxSize/posFac))+"Mpc hInd="+str(gcID)+" z=start"
    cgPlot, [0], [0], /nodata, xrange=[0,sP.boxSize/posFac], yrange=[0,sP.boxSize/posFac], /xs, /ys, $
      xtitle="axis0 [Mpc]", ytitle="axis1 [Mpc]", title=pTitle
    
    plots, zstart_pos[axes[0],*]/posFac,          zstart_pos[axes[1],*]/posFac,      psym=3 ;black
    plots, zstart_pos[axes[0],zstart_ind]/posFac, zstart_pos[axes[1],zstart_ind]/posFac, psym=3, color=cgColor('red')
    
    cgPlot, boundBox.(axes[0])[[0,1,1,0,0]]/posFac, boundBox.(axes[1])[[0,0,1,1,0]]/posFac, $
      line=0, color=cgColor('orange'), /overplot
  end_PS, pngResize=40, /deletePS
  endif ;0
  
  stop

end

; orderLagrangianVolumes(): take a selection of halos at z=target based on desired mass and rank the volumes of 
;   their bounding lagrangian regions at z=start (choose the smallest/smaller to resimulate)
; NOTE: assumes snapshot zero is the ICs at z=zstart

pro orderLagrangianVolumes

  ; config
  haloMassRange  = [11.8,12.2] ; log msun
  targetRedshift = 2.0 ; resimulate until this redshift only
  zoomLevel      = 4    ; levelmax-levelmin
  
  sP = simParams(res=128,run='zoom_20Mpc_dm',redshift=targetRedshift)

  ; Onorbe+ (2013) suggested R_tb / rVirFac (requires >~4000 particles in halo in fullbox
  ; we just about satisfy this, ~2500-3500 for a 10^12 halo @ 20Mpc/128 or 40Mpc/256
  rVirFac = (-0.1 * zoomLevel + 4) ; shy's formula
  print,'Using rVirFac: '+str(rVirFac)+' for zoomLevel='+str(zoomLevel)
  
  ; load z=target
  ; -------------
  gc = loadGroupCat(sP=sP,/readIDs)
  
  gcIDs    = gcIDList(gc=gc,select='pri')
  gcMasses = codeMassToLogMsun( gc.subgroupMass[ gcIDs ] )
  
  gcInds = where(gcMasses ge haloMassRange[0] and gcMasses lt haloMassRange[1], count)
  w = where(gcMasses gt haloMassRange[1], count2)
  
  gcIDs    = gcIDs[ gcInds ]
  gcMasses = gcMasses[ gcInds ]
  
  print,'Found ['+str(count)+'] halos in mass range, ['+str(count2)+'] more massive exist.'
  
  ; load z=0 dm_ids and dm_pos
  ztarget_ids = loadSnapshotSubset(sP=sP,partType='dm',field='ids')
  ztarget_pos = loadSnapshotSubset(sP=sP,partType='dm',field='pos')

  ; loop over each halo
  halos = {}
  foreach gcID,gcIDs,k do begin
    if k mod round(n_elements(gcIDs)/5.0) eq 0 then print,string(float(k)/n_elements(gcIDs)*100,format='(f4.1)')+'%'
    
    ; calculate radial distance of each DM particle to halo center
    gc_rVir = gc.group_r_crit200[ gc.subgroupGrNr[ gcID ] ]
    rad = periodicDists( gc.subgroupPos[*,gcID],ztarget_pos,sP=sP)
    
    ; select indices and ids
    ztarget_ind = where(rad le rVirFac * gc_rVir, count)
    if count eq 0 then message,'Error'
    
    ; add entry to halos
    halo = { hInd : gcID                     ,$
             pos  : gc.subgroupPos[*,gcID]   ,$
             rvir : gc_rVir                  ,$
             mass : gcMasses[k]              ,$
             ids  : ztarget_ids[ztarget_ind]  }
             
    halos = mod_struct( halos, 'hInd'+str(gcID), halo )
             
  endforeach
  
  ; load z=zstart dm_ids
  ; --------------------
  sP.snap = sP.snapRange[0]
  sP.redshift = snapNumToRedshift(sP=sP)
  
  zstart_ids = loadSnapshotSubset(sP=sP,partType='dm',field='ids')
  zstart_pos = loadSnapshotSubset(sP=sP,partType='dm',field='pos')
  
  ; make id->index map
  idIndexMap = getIDIndexMap(zstart_ids,minid=minid)
  
  ; loop over each halo
  haloStats = {}
  volumes   = []
  
  print,''
  print,'Candidate halos at z='+string(targetRedshift,format='(f3.1)')
  print,''
  
  for k=0,n_tags(halos)-1 do begin
    halo = halos.(k)
    
    zstart_ind = idIndexMap[ halo.ids - minid ]
    zstart_pos_halo = zstart_pos[*,zstart_ind]
    
    ; center of mass, bounding box, and volume
    massCenter = { x : mean(zstart_pos_halo[0,*]),   y : mean(zstart_pos_halo[1,*]),   z : mean(zstart_pos_halo[2,*]) }
    boundBox   = { x : minmax(zstart_pos_halo[0,*]), y : minmax(zstart_pos_halo[1,*]), z : minmax(zstart_pos_halo[2,*]) }
    volume     = (boundBox.x[1]-boundBox.x[0]) * (boundBox.y[1]-boundBox.y[0]) * (boundBox.z[1]-boundBox.z[0])
    volume    /= (1000.0^3) ;kpc^3 -> Mpc^3
    
    ; box center and extent
    boxCenter = { x : mean(boundBox.x), y : mean(boundBox.y), z : mean(boundBox.z) }
    boxExtent = { x : boundBox.x[1]-boundBox.x[0], y : boundBox.y[1]-boundBox.y[0], z : boundBox.z[1]-boundBox.z[0] }  
    
    ; check if bounds wrap periodic edge? if so, cannot place highres region on this halo
    wrap = 0 ; ok
    for i=0,2 do $
      if abs(boundBox.(i)[1]-boundBox.(i)[0]) gt sP.boxSize/2 then wrap = 1 ; bad
      
    ; add entry to haloStats
    haloStat = { massCenter:massCenter, boundBox:boundBox, volume:volume, $
                 boxCenter:boxCenter, boxExtent:boxExtent, wrap:wrap }
                 
    haloStats = mod_struct(haloStats, 'hInd'+str(halo.hInd), haloStat )
    volumes = [volumes,volume]
    
    ; if we don't wrap, print out this halo
    if wrap eq 0 then $
      print,'['+string(k,format='(i2)')+']'+$
            ' hInd = '+string(halo.hInd,format='(i4)')+$
            ' mass = '+string(halo.mass,format='(f5.2)')+$
            ' rvir = '+string(halo.rVir,format='(f5.1)')+$
            ' vol = '+ string(haloStat.volume,format='(f5.1)')+$
            ' pos = ['+string(halo.pos[0],format='(f8.2)')+' '+$
                       string(halo.pos[1],format='(f8.2)')+' '+$
                       string(halo.pos[2],format='(f8.2)')+' ]'
    
  endfor
 
  ; sort by ascending volume, then print out the fields required for MUSIC ICs
  sort_inds = sort(volumes)
  
  print,''
  print,'Sorted by ascending Lagrangian volume, non-wrapping final candidates in mass range:'
  print,'(Showing MUSIC fractional boxCenter and boxExtent parameters at z='+str(fix(sP.redshift))+')'
  print,''
  
  foreach ind,sort_inds do begin
    if haloStats.(ind).wrap eq 1 then continue
    
    halo = halos.(ind)
    haloStat = haloStats.(ind)
    
    print,'['+string(ind,format='(i2)')+']'+$
          ' hInd = '+string(halo.hInd,format='(i4)')+$
          ' mass = '+string(halo.mass,format='(f5.2)')+$
          ' rvir = '+string(halo.rVir,format='(f5.1)')+$
          ' vol = '+ string(haloStat.volume,format='(f5.1)')+$
          ' zL = '+str(zoomLevel)+$
          ' rVirFac = '+string(rVirFac,format='(f3.1)')+$
          ' pos = ['+string(halo.pos[0],format='(f8.2)')+' '+$
                     string(halo.pos[1],format='(f8.2)')+' '+$
                     string(halo.pos[2],format='(f8.2)')+' ]'
    print,'      boxCenter: [ '+$
      string(haloStat.boxCenter.x/sP.boxSize,format='(f7.4)')+' , '+$
      string(haloStat.boxCenter.y/sP.boxSize,format='(f7.4)')+' , '+$
      string(haloStat.boxCenter.z/sP.boxSize,format='(f7.4)')+' ]'
    print,'      boxExtent: [ '+$
      string(haloStat.boxExtent.x/sP.boxSize,format='(f7.4)')+' , '+$
      string(haloStat.boxExtent.y/sP.boxSize,format='(f7.4)')+' , '+$
      string(haloStat.boxExtent.z/sP.boxSize,format='(f7.4)')+' ]' 
  endforeach
  
  
  stop
end

; zoomTargetHalo(): locate our target halo in a zoom snapshot

function zoomTargetHalo, sP=sP, gc=gc

  if n_elements(sP) eq 0 then message,'Error: Must supply sP.'
  if n_elements(gc) eq 0 then gc = loadGroupCat(sP=sP,/skipIDs)

  expectedHaloPos = sP.targetHaloPos + sP.zoomShiftPhys
  
  massDiffs = abs( sP.targetHaloMass - codeMassToLogMsun(gc.subgroupMass) )
  hInd = ( where(massDiffs eq min(massDiffs),count) )[0]
  if count ne 1 then message,'Error'
  
  posDiff = expectedHaloPos - gc.subgroupPos[*,hInd]
  foundMass = codeMassToLogMsun( gc.subgroupMass[hInd] )
  foundRvir = gc.group_r_crit200[ gc.subgroupGrNr[hInd] ]
  distNorm = sqrt( posDiff[0]^2 + posDiff[1]^2 + posDiff[2]^2 ) / foundRvir
  
  print,' zoom hInd: '+str(hInd)
  print,' targetMass ['+str(sP.targetHaloMass)+'] foundMass ['+str(foundMass)+']'
  print,' targetRvir ['+str(sP.targetHaloRvir)+'] foundRvir ['+str(foundRvir)+']'
  print,' posDiff: ',posDiff
  print,' distNorm: ',distNorm
  
  return, hInd
end

; checkDMContamination(): look for intrusion of lowres DM particles into the Lagrangian region of the target halo

pro checkDMContamination

  ; config
  sPs = mod_struct( sPs, 'sP0', simParams(run='zoom_20Mpc',res=9,hInd=0,redshift=2.0) )
  sPs = mod_struct( sPs, 'sP1', simParams(run='zoom_20Mpc',res=9,hInd=1,redshift=2.0) )
  ;sPs = mod_struct( sPs, 'sP2', simParams(run='zoom_20Mpc_dm',res=9,hInd=0,redshift=2.0) )
  ;sPs = mod_struct( sPs, 'sP3', simParams(run='zoom_20Mpc_dm',res=10,hInd=0,redshift=2.0) )
  ;sPs = mod_struct( sPs, 'sP4', simParams(run='zoom_20Mpc_dm',res=11,hInd=0,redshift=2.0) )
  
  colors = ['red','blue','green','orange']
  nn     = 99
  
  start_PS,sPs.(0).plotPath + 'radialContamination_DM.eps'
    cgPlot,[0],[0],/nodata,xrange=[1,3.5],yrange=[0.9,100],/xs,/ys,/ylog,yminor=0,$
      xtitle=textoidl("r / r_{vir}"),ytitle="Cumulative Number of Lowres DM Within Radius"
  
  strings = []
  
  for i=0,n_tags(sPs)-1 do begin

    ; load
    sP = sPs.(i)
    print,'RUN: '+sP.run+' [res='+str(sP.res)+']'
    gc   = loadGroupCat(sP=sP,/skipIDs)
    hInd = zoomTargetHalo(sP=sP, gc=gc) ; which subhalo is our target halo?
  
    ; group catalog check
    nCoarseSubhalo = gc.subgroupLenType[ partTypeNum('lowres_dm'),hInd ]
    nCoarseFof     = gc.groupLenType[ partTypeNum('lowres_dm'), gc.subgroupGrNr[hInd] ]
  
    print,' nCoarseSubhalo ['+str(nCoarseSubhalo)+'] nCoarseFof ['+str(nCoarseFof)+']'
  
    ; spatial check: load lowres DM positions
    rVirFac = (1.5 * sP.zoomLevel + 1)
    pos = loadSnapshotSubset(sP=sP,partType='lowres_dm',field='pos')
  
    dists = periodicDists(gc.subgroupPos[*,hInd],pos,sP=sP)
    dists /= gc.group_r_crit200[ gc.subgroupGrNr[hInd] ] ; divide by found rvir
  
    w = where(dists le rVirFac,nRvirFac)
    print,' Spatial within rVirFac='+str(rVirFac)+' have ['+str(nRvirFac)+'] coarse DM.'
    w = where(dists le 1.5,n15rvir)
    print,' Spatial within rVirFac=1.5 have ['+str(n15rvir)+'] coarse DM.'
  
    if n15rvir eq 0 and nCoarseSubhalo eq 0 and nCoarseFof eq 0 then print,' - PASSED -' else print,' - FAILED -'
  
    ; out to what radius are we uncontaminated?
    dists = dists[sort(dists)]
    print,' Zero lowres DM out to radius of ['+string(min(dists),format='(f4.2)')+' rvir] ('+$
      string(min(dists)/rVirFac,format='(f4.2)')+')'
  
    ; add cumulative histogram to plot
    strings = [strings, str(sP.simName)+' h'+str(sP.hInd)+' ('+str(2^sP.res)+')']
    cgPlot,dists[0:nn],lindgen(nn),psym=-4,/overplot,color=cgColor(colors[i])
    
  endfor
  
  legend,strings,textcolor=colors,/top,/left
  
  end_PS
  
  stop

end

; checkGasContamination(): look for the intrusion of MC tracers starting in lowres gas cells

pro checkGasContamination

  ; config
  sP = simParams(run='zoom_20Mpc',res=9,hInd=1,redshift=2.0)

  ; load
  h = loadSnapshotHeader(sP=sP)

  tr_ids    = loadSnapshotSubset(sP=sP, partType='tracerMC', field='tracerids')
  tr_parids = loadSnapshotSubset(sP=sP, partType='tracerMC', field='parentids')

  ; trIDs > IDS_OFFSET*trMCPerCell are in lowres
  ; trIDs < IDS_OFFSET*trMCPerCell are in highres
  w_tr_highres = where(tr_ids lt sP.ids_offset * sP.trMCPerCell, comp=w_tr_lowres)
  
  tr_ids_highres    = tr_ids[ w_tr_highres ]
  tr_ids_lowres     = tr_ids[ w_tr_lowres ]
  tr_parids_highres = tr_parids[ w_tr_highres ]
  tr_parids_lowres  = tr_parids[ w_tr_lowres ]
  
  ; sanity check
  calcMatch,tr_ids,tr_ids_highres,ind1,ind_high,count=count1
  calcMatch,tr_ids,tr_ids_lowres,ind2,ind_low,count=count2
  
  if count1+count2 ne n_elements(tr_ids) then message,'Fail'
  if intersection(tr_ids_highres,tr_ids_lowres) ne -1 then message,'Fail'
  if intersection(ind1,ind2) ne -1 then message,'Fail'
  
  mask = bytarr(n_elements(tr_ids))
  mask[ind1] = 1B
  mask[ind2] = 1B
  w = where(mask eq 0B,count)
  if count gt 0 then message,'Fail'
  
  if count1 ne h.nPartTot[ partTypeNum('highres_dm') ] * sP.trMCPerCell then message,'Fail'
  if count2 ne h.nPartTot[ partTypeNum('lowres_dm') ] * sP.trMCPerCell then message,'Fail'
  
  ; tracer IDs are ok, load gas cells by group catalog
  gc   = loadGroupCat(sP=sP,/readIDs)
  hInd = zoomTargetHalo(sP=sP, gc=gc) ; which subhalo is our target halo?
  
  ; primary subgroup: all=gas+stars
  par_ids = gcPIDList(gc=gc,valGCids=[hInd],partType='all')
  par_ids = par_ids[sort(par_ids)]
  
  ; crossmatch lowres (e.g. starting in lowres gas) tracer parentids to target halo ids
  par_ind = value_locate(par_ids,tr_parids_lowres) > 0
  w = where(par_ids[par_ind] eq tr_parids_lowres,count_inSG) ; verify we actually matched the ID
  
  ; repeat for full fof group, use all types (no tracers have DM parents)
  grNr = gc.subgroupGrNr[hInd]
  par_ids = lindgen(gc.groupLen[grNr]) + gc.groupOffset[grNr]
  par_ids = gc.ids[par_ids]
  par_ids = par_ids[sort(par_ids)]
  
  par_ind = value_locate(par_ids,tr_parids_lowres) > 0
  w = where(par_ids[par_ind] eq tr_parids_lowres,count_inFoF) ; verify we actually matched the ID
  
  print,'Found lowres tracers: ['+str(count_inSG)+'] in subgroup and ['+str(count_inFoF)+'] in fof.'
  
  ; want: dists to each lowres tracer, first change tr_parids to tr_parinds then calculate par dists
  dists = fltarr(n_elements(tr_parids_lowres))
  mask  = intarr(n_elements(tr_parids_lowres))
  
  foreach parType,['gas','stars'] do begin
    ; load
    par_ids = loadSnapshotSubset(sP=sP, partType=parType, field='ids')
    par_pos = loadSnapshotSubset(sP=sP, partType=parType, field='pos')
  
    ; crossmatch
    sort_inds = calcSort(par_ids)
    par_ids_sorted = par_ids[sort_inds]
    par_ind = value_locate(par_ids_sorted,tr_parids_lowres)
    par_ind = sort_inds[par_ind>0]
    w = where(par_ids[par_ind] eq tr_parids_lowres,count) ; indices to tr_parids_lowres, dists

    if count eq 0 then continue ; none, skip this parent type
    
    tr_parids_inPar = par_ind[w] ; indices to par_ids, par_pos, par_dists
  
    ; calculate par dists (all)
    par_dists = periodicDists(gc.subgroupPos[*,hInd],par_pos,sP=sP)
    par_dists /= gc.group_r_crit200[ gc.subgroupGrNr[hInd] ] ; divide by found rvir
  
    dists[w] = par_dists[tr_parids_inPar] ; stamp matched parents in
    mask[w] += 1 ; debug
  endforeach
  
  w = where(mask ne 1,count)
  if count gt 0 then message,'Calculating lowres tracer dists FAIL.'
  
  ; out to what radius are we uncontaminated?
  dists = dists[sort(dists)]
  rVirFac = (1.5 * sP.zoomLevel + 1)
  print,' Zero lowres TR out to radius of ['+string(min(dists),format='(f4.2)')+' rvir] ('+$
    string(min(dists)/rVirFac,format='(f4.2)')+')'
  
  ; plot
  start_PS,sP.plotPath + 'radialContamination_Gas.eps'
    cgPlot,[0],[0],/nodata,xrange=[1,3.5],yrange=[0.9,100],/xs,/ys,/ylog,yminor=0,$
      xtitle=textoidl("r / r_{vir}"),ytitle="Cumulative Number of Lowres TR Within Radius"
      
    nn = 99
    cgPlot,dists[0:nn],lindgen(nn),psym=-4,/overplot,color=cgColor('red')
    
    legend,[str(sP.simName)+' h'+str(sP.hInd)+' ('+str(2^sP.res)+')'],textcolor=['red'],/top,/left
  end_PS
  
  stop

end

; checkMusicHDF5ICs(): verify output of the arepo hdf5 plugin for MUSIC agrees with the gadget output plugin
;  note: gas and highres DM are binary identical, but coarse DM is different at the level of roundoff
;  e.g. boxSize_kpc*1e-7 which is ~10pc for 20/40Mpc boxes, due to saving the intermediate values 
;  as float not double for combining the coarse components
;  note: using OMP_NUM_THREADS>1 also changes all the outputs at the level of roundoff

@arepoLoad

pro checkMusicHDF5ICs

  forward_function loadSnapshotOld

  ;fpath = '/n/home07/dnelson/sims.zooms/ICs/128_20Mpc_dm/output/' ; pass OMP=1
  ;fpath = '/n/home07/dnelson/sims.zooms/ICs/128_20Mpc_h0/output_L8_dmonly/' ; pass OMP=1
  fpath = '/n/home07/dnelson/sims.zooms/ICs/128_20Mpc_h0/output_L9_dmonly/' ; pass OMP=1
  ;fpath = '/n/home07/dnelson/sims.zooms/ICs/128_20Mpc_h0/output_L8/' ; pass OMP=1
  ;fpath = '/n/home07/dnelson/sims.zooms/ICs/128_20Mpc_h0/output_L9/' ; pass OMP=1
  
  old = loadSnapshotOld(fpath + 'ics_gadget.dat', 'none')
  new = h5_parse(fpath + 'ics_arepo.hdf5', /read)
  
  ; check counts and masstable
  old_npart = old.header.nparttot
  new_npart = new.header.numpart_total._data
  
  old_massTable = old.header.massTable
  new_massTable = new.header.massTable._data
  diff_massTable = max(abs(old_massTable-new_massTable))
  
  print,'npart',array_equal(old_npart,new_npart)
  print,'masstable diff: ',diff_massTable
  
  ; compare gas,highresDM,coarseDM
  foreach partType,[0,1,2] do begin
    tagInd_old = where( (tag_names(old)) eq 'PARTTYPE'+str(partType), count1)
    tagInd_new = where( (tag_names(new)) eq 'PARTTYPE'+str(partType), count2)
    
    if count1 eq 0 or count2 eq 0 then begin
      print,'Skip: '+str(partType)
      continue
    endif
    
    old_pt_pos = old.(tagInd_old).pos
    new_pt_pos = new.(tagInd_new).coordinates._data
    old_pt_vel = old.(tagInd_old).vel
    new_pt_vel = new.(tagInd_new).velocities._data
    old_pt_ids = old.(tagInd_old).ids
    new_pt_ids = new.(tagInd_new).particleids._data
  
    ; old is reordered due to temporary file process
    old_sort = sort(old_pt_ids)
    old_pt_pos = old_pt_pos[*,old_sort]
    old_pt_vel = old_pt_vel[*,old_sort]
    old_pt_ids = old_pt_ids[old_sort]
  
    print,'PT'+str(partType)+' count',n_elements(old_pt_pos[0,*]),n_elements(new_pt_pos[0,*])
    print,'PT'+str(partType)+' pos',array_equal(old_pt_pos,new_pt_pos)
    print,'PT'+str(partType)+' vel',array_equal(old_pt_vel,new_pt_vel)
    print,'PT'+str(partType)+' ids',array_equal(old_pt_ids,new_pt_ids)
  
    ; masses exist?
    if tag_exist(new.(tagInd_new), 'masses') then begin
      old_pt_mass = old.(tagInd_old).mass
      new_pt_mass = new.(tagInd_new).masses._data
      old_pt_mass = old_pt_mass[old_sort]
      print,'PT'+str(partType)+' mass',array_equal(old_pt_mass,new_pt_mass)
    endif
  
    ; position/velocity difference?
    vel_diff = old_pt_vel - new_pt_vel
	pos_diff = old_pt_pos - new_pt_pos
    vel_diff = max(abs(vel_diff))
	pos_diff = max(abs(pos_diff))
    print,'PT'+str(partType)+': vel diff magnitude: ',vel_diff
	print,'PT'+str(partType)+': pos diff magnitude: ',pos_diff

  endforeach
  
  stop
end

; plotCAMB(): plot transfer functions and powerspectrum output from CAMB, check realspace TF from MUSIC

pro plotCAMB

  f1 = '/n/home07/dnelson/make.ics/camb.run/zoom_'
  f2 = '/n/home07/dnelson/sims.zooms/ICs/256_40Mpc_dm/transfer_real_total.txt'
  
  ; matter power spectrum
  headerLines = 0
  template    = { k:0.0, p:0.0 }
  
  mps = loadCSV(headerLines, f1 + 'matterpower.dat', template)
  
  ; transfer function ( wavenumber, CDM, baryon, photon, massless neutrino, massive neutrino, total )
  template = { k:0.0, cdm:0.0, b:0.0, g:0.0, r:0.0, nu:0.0, tot:0.0 }
  tf = loadCSV(headerLines, f1 + 'transfer_out.dat', template)
  
  template = { r:0.0, tf:0.0, err:0.0 }
  tf_real = loadCSV(headerLines, f2, template)
  
  ; plot
  start_PS,'camb_matterpower.eps'
    cgPlot, [0], [0], /nodata, /xlog, /ylog, xminor=0, yminor=0, xtitle="k (h/Mpc)", ytitle="P(k)", $
      xrange=minmax(mps.k), yrange=minmax(mps.p)
      
    cgPlot, mps.k, mps.p, /overplot, color=cgColor('red')
  end_PS
  
  start_PS,'camb_transferfunc.eps'
    cgPlot, [0], [0], /nodata, /xlog, /ylog, xminor=0, yminor=0, xtitle="k/h", ytitle="T(k)", $
      xrange=minmax(tf.k), yrange=minmax(tf.tot)
    cgPlot, tf.k, tf.tot * tf.k^2, /overplot, color=cgColor('black'), line=0
    cgPlot, tf.k, tf.cdm * tf.k^2, /overplot, color=cgColor('red'), line=0
    cgPlot, tf.k, tf.b   * tf.k^2, /overplot, color=cgColor('blue'), line=0
    
    legend,['CDM','baryon'],textcolor=['red','blue'],/top,/right
  end_PS
  
  start_PS,'camb_transferreal.eps'
    cgPlot, [0], [0], /nodata, /xlog, /ylog, xminor=0, yminor=0, xtitle="r (Mpc/h)", ytitle="T(r)", $
      xrange=minmax(tf_real.r), yrange=[1e-22,max(tf_real.tf)], /xs, /ys
    cgPlot, tf_real.r, tf_real.tf, /overplot, color=cgColor('black'), line=0
    cgPlot, [40.0,40.0], [1e-20,1e5], line=1, color=cgColor('red'), /overplot ; 40Mpc
    cgPlot, [10.0/1000^2, 10.0/1000^2], [1e-20,1e5], line=1, color=cgColor('green'), /overplot ;10pc
  end_PS
  stop

end
