; ICs_filamentTest.pro
; idealized filament impacting hot halo atmosphere test
; dnelson mar.2012

; filamentPlot

pro filamentPlot
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  basePath = '/n/home07/dnelson/sims.idealized/'
  d = '1000'
  
  sPa = { simPath   : basePath+'cylTest.1e4.M1e2.norot.c10.d'+d+'.r2.arepo/output/' ,$
          plotPath  : basePath ,$
          derivPath : basePath+'cylTest.1e4.M1e2.norot.c10.d'+d+'.r2.arepo/data.files/',$
          savPrefix : 'F',$
          res       : '1',$
          boxSize   : 5000.0,$
          snap      : 0 }
  sPg = { simPath   : basePath+'cylTest.1e4.M1e2.norot.c10.d'+d+'.r2.gadget/output/' ,$
          plotPath  : basePath ,$
          derivPath : basePath+'cylTest.1e4.M1e2.norot.c10.d'+d+'.r2.gadget/data.files/',$
          savPrefix : 'F',$
          res       : '1',$
          boxSize   : 5000.0,$
          snap      : 0 }
  
  haloTvir = alog10(codeMassToVirTemp(100.0,redshift=2.0))
  haloRvir = 162.6 ;kpc
         
  minFilID = 10001L ; for cylTest.1e4.c10.r2
  maxFilID = 15000L ; for cylTest.1e4.c10.r2
  
  nGasFil = maxFilID - minFilID + 1
  
  snapRange = [0,80]
  nSnaps = snapRange[1]-snapRange[0] + 1
  
  ; save/restore
  saveFilename = sPa.derivPath + 'filsave_'+str(snapRange[0])+'_'+str(snapRange[1])+$
    '-'+str(minFilID)+'_'+str(maxFilID)+'.sav'
    
  if file_test(saveFilename) then begin
    restore,saveFilename,/verbose
  endif else begin
  
    ; arepo: load snap=0 and select filament gas cells
    gas_ids = loadSnapshotSubset(sP=sPa,partType='gas',field='ids')
    
    w_gas = where(gas_ids ge minFilID and gas_ids le maxFilID,count_gas)
    gas_ids = gas_ids[w_gas]
    
    ; find all MC/vel tracer children
    trids_mc  = cosmoTracerChildren(   sP=sPa,/getIDs,gasIDs=gas_ids,child_counts=child_counts_mc)
    trids_vel = cosmoTracerVelChildren(sP=sPa,/getIDs,gasInds=w_gas, child_counts=child_counts_vel)
  
    ; arrays
    times = fltarr(nSnaps)
    
    rad  = { ga  : fltarr(nSnaps,nGasFil)               ,$
             mc  : fltarr(nSnaps,n_elements(trids_mc))  ,$
             vel : fltarr(nSnaps,n_elements(trids_vel))  }
    tmax = { ga  : fltarr(nSnaps,nGasFil)               ,$
             mc  : fltarr(nSnaps,n_elements(trids_mc))  ,$
             vel : fltarr(nSnaps,n_elements(trids_vel))  }    
            
    ; loop through snapshots and record tracer details
    k = 0
    for m=snapRange[0],snapRange[1],1 do begin
      sPa.snap = m
      sPg.snap = m
      print,m
      ; load header and save time
      h = loadSnapshotHeader(sP=sPa)
      times[k] = h.time
      
      ; load tracer IDs
      loc_trids_mc  = loadSnapshotSubset(sP=sPa,partType='tracerMC',field='tracerids')
      loc_trids_vel = loadSnapshotSubset(sP=sPa,partType='tracerVel',field='ids')
      
      ; match starting to local tracers
      match,trids_mc,loc_trids_mc,ind_global,ind_loc_mc,count=count_mc
      match,trids_vel,loc_trids_vel,ind_global,ind_loc_vel,count=count_vel
      
      ; load MC tracer parents, gas IDs and match
      parids_mc = loadSnapshotSubset(sP=sPa,partType='tracerMC',field='parentid')
      gas_ids   = loadSnapshotSubset(sP=sPa,partType='gas',field='ids')
      
      parids_mc = parids_mc[ind_loc_mc] ; subselect gas IDs for parents only
      
      gasIDMap = getIDIndexMap(gas_ids,minid=minid)
      gas_ids  = !NULL
      
      parinds_mc = gasIDMap[parids_mc-minid]
      gasIDMap = !NULL
      
      ; load gas positions and transform to tracerMC positions
      gas_pos = loadSnapshotSubset(sP=sPa,partType='gas',field='pos')
      gas_rad = periodicDists([sPa.boxSize/2,sPa.boxSize/2,sPa.boxSize/2],gas_pos,sP=sPa)
      gas_pos = !NULL
      
      rad.mc[k,*] = gas_rad[parinds_mc]
      gas_rad = !NULL
      
      ; load VEL tracer positions
      trvel_pos = loadSnapshotSubset(sP=sPa,partType='tracerVel',field='pos')
      trvel_rad = periodicDists([sPa.boxSize/2,sPa.boxSize/2,sPa.boxSize/2],trvel_pos,sP=sPa)
      trvel_pos = !NULL
      
      rad.vel[k,*] = trvel_rad[ind_loc_vel]
      trvel_rad = !NULL
      
      ; load tmax for both tracer types and save
      tmax_mc  = loadSnapshotSubset(sP=sPa,partType='tracerMC',field='tracer_maxtemp')
      tmax_vel = loadSnapshotSubset(sP=sPa,partType='tracerVel',field='tracer_maxtemp')
      
      tmax_mc  = tmax_mc[ind_loc_mc]
      tmax_vel = tmax_vel[ind_loc_vel]
      
      ; replace tmax if new maxima
      tmax.mc[k,*]  = tmax_mc
      tmax.vel[k,*] = tmax_vel
      
      tmax_mc  = !NULL
      tmax_vel = !NULL
      
      ; gadget: load gas particle ids and positions and save filament gas radii
      gas_ids = loadSnapshotSubset(sP=sPg,partType='gas',field='ids')
      
      w_gas_ga = where(gas_ids ge minFilID and gas_ids le maxFilID,count_gas)
      gas_ids = !NULL
      
      gas_pos = loadSnapshotSubset(sP=sPg,partType='gas',field='pos')
      gas_rad = periodicDists([sPa.boxSize/2,sPa.boxSize/2,sPa.boxSize/2],gas_pos,sP=sPa)
      gas_pos = !NULL
      
      rad.ga[k,*] = gas_rad[w_gas_ga]
      
      ; gadget: load gas utherm,ne and save temps if new maxima
      gas_u  = loadSnapshotSubset(sP=sPg,partType='gas',field='u')
      gas_ne = loadSnapshotSubset(sP=sPg,partType='gas',field='ne')
      
      gas_temp = convertUtoTemp(gas_u,gas_ne)
      gas_temp = alog10(gas_temp[w_gas_ga])
      
      tmax.ga[k,*] = gas_temp ; actually current temp for sph
      
      k += 1
    endfor
    
    ; tracer output still in unit system. convert temperatures from code to log(K)
    tmax.mc  = codeTempToLogK(tmax.mc)
    tmax.vel = codeTempToLogK(tmax.vel)
    
    ; make a global maximum for mc,vel,ga over all snapshots
    tmax_all = { mc  : max(tmax.mc,dimension=1)  ,$
                 vel : max(tmax.vel,dimension=1) ,$
                 ga  : max(tmax.ga,dimension=1)   }
               
    ; save
    save,times,tmax,rad,tmax_all,snapRange,sPa,sPg,filename=saveFilename
  endelse ; file_test
  
  binSizeRad  = 0.2
  binSizeTemp = 0.02
  
  ; plot (1) - radial distribution of tracers
  start_PS, sPa.plotPath + 'cyl_rad_d'+d+'.eps'
    ;!p.thick += 1
    xrange = [0,1.9]
    yrange = [5e-4,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="fraction",xtitle=textoidl("R_{tr} / R_{vir}")+""
    cgPlot,[1,1],yrange,line=0,color=fsc_color('light gray'),/overplot
    cgPlot,[1,1]+100.0/haloRvir,yrange,line=0,color=fsc_color('light gray'),/overplot
    
    for k=0,nSnaps-1,1 do begin
      ; histogram tracerMC and tracerVEL radii normalized by virial radius of halo
      h_mc  = histogram(rad.mc[k,*]/haloRvir,min=0.0,binsize=binSizeRad,loc=loc_mc)
      h_vel = histogram(rad.vel[k,*]/haloRvir,min=0.0,binsize=binSizeRad,loc=loc_vel)
      h_ga  = histogram(rad.ga[k,*]/haloRvir,min=0.0,binsize=binSizeRad,loc=loc_ga)
      
      ; plot fractional histograms
      cgPlot,loc_mc,float(h_mc)/total(h_mc),line=0,/overplot,color=getColor(k)
      cgPlot,loc_vel,float(h_vel)/total(h_vel),line=1,/overplot,color=getColor(k)
      cgPlot,loc_ga,float(h_ga)/total(h_ga),line=2,/overplot,color=getColor(k)
    endfor
    
    ; legend
    legend,string(times,format='(f4.2)')+' Gyr',textcolors=getColor(indgen(nSnaps),/name),box=0,/left,/top
    legend,['MC','VEL','SPH'],linestyle=[0,1,2],box=0,linesize=0.25,/right,/top
    
  end_PS
  
  ; plot (2) - tmax histograms at final snapshot
  start_PS, sPa.plotPath + 'cyl_tmax_final_d'+d+'.eps'
    ;!p.thick += 1
    xrange = [0.5,1.5]
    yrange = [5e-4,1.0]
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,$
      ytitle="fraction",xtitle=textoidl("log(T_{max,tr}) / log(T_{vir,halo})")+""
    cgPlot,[1,1],yrange,line=0,color=fsc_color('light gray'),/overplot
    
    ; histogram tracerMC and tracerVEL radii normalized by virial radius of halo
    h_mc  = histogram(tmax_all.mc[*]/haloTvir,binsize=binSizeTemp,loc=loc_mc)
    h_vel = histogram(tmax_all.vel[*]/haloTvir,binsize=binSizeTemp,loc=loc_vel)
    h_ga  = histogram(tmax_all.ga[*]/haloTvir,binsize=binSizeTemp,loc=loc_ga)
    
    ; plot fractional histograms
    cgPlot,loc_mc,float(h_mc)/total(h_mc),line=0,/overplot,color=getColor(0)
    cgPlot,loc_vel,float(h_vel)/total(h_vel),line=1,/overplot,color=getColor(0)
    cgPlot,loc_ga,float(h_ga)/total(h_ga),line=2,/overplot,color=getColor(0)

    ; legend
    legend,['MC','VEL','SPH'],textcolors=getColor([0,0,0],/name),linestyle=[0,1,2],$
      linesize=0.25,box=0,/right,/top
    
  end_PS

stop
end

; filGetPositions(): load positions(time) and track through snapshots

function filGetPositions, sPa=sPa, sPg=sPg, snapRange=snapRange, minFilID=minFilID, maxFilID=maxFilID

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; deriv
  nGasFil   = maxFilID - minFilID + 1
  gasFilIDs = lindgen(nGasFil) + minFilID
  
  nSnaps = snapRange[1]-snapRange[0] + 1

  ; save/restore
  saveFilename = sPa.derivPath + 'filpos_'+str(snapRange[0])+'_'+str(snapRange[1])+$
    '-'+str(minFilID)+'_'+str(maxFilID)+'.sav'
    
  if file_test(saveFilename) then begin
    restore,saveFilename,/verbose
    return,pos
  endif
  
  ; arepo: load snap=0 and select filament gas cells
  gas_ids = loadSnapshotSubset(sP=sPa,partType='gas',field='ids')
  
  w_gas = where(gas_ids ge minFilID and gas_ids le maxFilID,count_gas)
  gas_ids = gas_ids[w_gas]
  
  ; find all MC/vel tracer children
  trids_mc  = cosmoTracerChildren(   sP=sPa,/getIDs,gasIDs=gas_ids,child_counts=child_counts_mc)
  trids_vel = cosmoTracerVelChildren(sP=sPa,/getIDs,gasInds=w_gas, child_counts=child_counts_vel)
  gas_ids = !NULL
  
  ; arrays
  pos  = { ga    : fltarr(nSnaps,3,nGasFil)               ,$
           mc    : fltarr(nSnaps,3,n_elements(trids_mc))  ,$
           vel   : fltarr(nSnaps,3,n_elements(trids_vel)) ,$
           times_gadget : fltarr(nSnaps)                  ,$
           times_arepo  : fltarr(nSnaps)                   }
           
  ; loop through snapshots and record tracer details
  k = 0
  for m=snapRange[0],snapRange[1],1 do begin
    sPa.snap = m
    sPg.snap = m
    print,m
    ; load header and save time
    h = loadSnapshotHeader(sP=sPa)
    pos.times_arepo[k] = h.time
    
    ; load tracer IDs
    loc_trids_mc  = loadSnapshotSubset(sP=sPa,partType='tracerMC',field='tracerids')
    loc_trids_vel = loadSnapshotSubset(sP=sPa,partType='tracerVel',field='ids')
    
    ; match starting to local tracers
    match,trids_mc,loc_trids_mc,ind_global_mc,ind_loc_mc,count=count_mc
    match,trids_vel,loc_trids_vel,ind_global_vel,ind_loc_vel,count=count_vel
    if count_mc ne n_elements(trids_mc) or count_vel ne n_elements(trids_vel) then message,'errortr'
    
    ind_loc_mc  = ind_loc_mc[sort(ind_global_mc)]
    ind_loc_vel = ind_loc_vel[sort(ind_global_vel)]
    
    ; load MC tracer parents, gas IDs and match
    parids_mc = loadSnapshotSubset(sP=sPa,partType='tracerMC',field='parentid')
    gas_ids   = loadSnapshotSubset(sP=sPa,partType='gas',field='ids')
    
    parids_mc = parids_mc[ind_loc_mc] ; subselect gas IDs for parents only
    
    gasIDMap = getIDIndexMap(gas_ids,minid=minid)
    gas_ids  = !NULL
    
    parinds_mc = gasIDMap[parids_mc-minid]
    gasIDMap = !NULL
    
    ; load gas positions and transform to tracerMC positions
    gas_pos = loadSnapshotSubset(sP=sPa,partType='gas',field='pos')
    gas_pos = gas_pos[*,parinds_mc]
    
    pos.mc[k,0,*] = gas_pos[0,*]
    pos.mc[k,1,*] = gas_pos[1,*]
    pos.mc[k,2,*] = gas_pos[2,*]
    
    ; load VEL tracer positions
    trvel_pos = loadSnapshotSubset(sP=sPa,partType='tracerVel',field='pos')
    trvel_pos = trvel_pos[*,ind_loc_vel]
    
    pos.vel[k,0,*] = trvel_pos[0,*]
    pos.vel[k,1,*] = trvel_pos[1,*]
    pos.vel[k,2,*] = trvel_pos[2,*]
    
    ; load SPH particle positions
    h = loadSnapshotHeader(sP=sPg)
    pos.times_gadget[k] = h.time
    
    gas_ids = loadSnapshotSubset(sP=sPg,partType='gas',field='ids')
    match,gas_ids,gasFilIDs,ind_global_ga,ind_loc_ga,count=count_ga
    if count_ga ne n_elements(gasFilIDs) then message,'errorga'

    ;ind_loc_ga = ind_loc_ga[sort(ind_global_ga)]
    ind_place = sort(ind_global_ga)
    gas_ids = !NULL

    gas_pos = loadSnapshotSubset(sP=sPg,partType='gas',field='pos')
    
    pos.ga[k,0,*] = gas_pos[0,ind_global_ga]
    pos.ga[k,1,*] = gas_pos[1,ind_global_ga]
    pos.ga[k,2,*] = gas_pos[2,ind_global_ga]
    
    gas_pos = !NULL
    
    k += 1
  endfor ;m
  
  ; save
  save,pos,snapRange,sPa,sPg,filename=saveFilename

  return,pos

end

; filInterpImage(): load positions(time) and make interpolated images of tracer/sph positions

pro filInterpImage

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  basePath = '/n/home07/dnelson/sims.idealized/'
  d = '500'
  
  sPa = { simPath   : basePath+'cylTest.1e4.M1e2.norot.c10.d'+d+'.r2.arepo/output/' ,$
          plotPath  : basePath+'frames.cylTest2/',$
          derivPath : basePath+'cylTest.1e4.M1e2.norot.c10.d'+d+'.r2.arepo/data.files/',$
          savPrefix : 'F',$
          res       : '1',$
          boxSize   : 5000.0,$
          snap      : 0 }
  sPg = { simPath   : basePath+'cylTest.1e4.M1e2.norot.c10.d'+d+'.r2.gadget/output/' ,$
          plotPath  : basePath ,$
          derivPath : basePath+'cylTest.1e4.M1e2.norot.c10.d'+d+'.r2.gadget/data.files/',$
          savPrefix : 'F',$
          res       : '1',$
          boxSize   : 5000.0,$
          snap      : 0 }
  
  haloTvir = alog10(codeMassToVirTemp(100.0,redshift=2.0))
  haloRvir = 162.6 ;kpc
         
  minFilID = 10001L ; for cylTest.1e4
  maxFilID = 15000L ; for cylTest.1e4
  
  snapRange = [0,72] ;--- 80@2Gyr

  ; get gas/tracer positions with time
  pos = filGetPositions(sPa=sPa,sPg=sPg,snapRange=snapRange,minFilID=minFilID,maxFilID=maxFilID)
  
  ; count
  nTrMC  = n_elements(pos.mc[0,0,*])
  nTrVel = n_elements(pos.vel[0,0,*])
  nSph   = n_elements(pos.ga[0,0,*])

  ; view configuration
  haloRvir = 162.6 ;kpc
  boxCen   = [sPa.boxSize/2.0, sPa.boxSize/2.0]; + [0,haloRvir];xz
  zoomSize = 330.0 ;kpc
  axes     = [0,2]
  
  ; movie configuration
  nFrames = 361 ;401  ;--- 80@2Gyr
  timeStart = 0.0 + 1e-4 ;Gyr
  timeEnd   = 1.8 - 1e-4 ;Gyr  ;--- 80@2Gyr
  
  if timeStart lt min(pos.times_arepo) then message,'Error: Requested timeStart beyond snapshot range.'
  if timeEnd gt max(pos.times_arepo) then message,'Error: Requested timeEnd beyond snapshot range.'
  
  ; trail length and color sequence
  trailLenTime = 0.05 ;Gyr

  ; deriv
  timeStep = (timeEnd - timeStart) / (nFrames-1)
  frameTimes = timeStep * findgen(nFrames)
  
  trailLenFrame = fix(trailLenTime / timeStep)
  
  ; setup fading trail
  rgb = fix(findgen(trailLenFrame)/(trailLenFrame-1) * 254)
  trailColors = lonarr(trailLenFrame)
  for i=0,trailLenFrame-1 do $
    trailColors[i] = rgb[i] + rgb[i]*(256L) + rgb[i]*(256L)^2

  trailColors = reverse(trailColors)
  
  ; spline interpolation of axes[0],axes[1] coordinates
  velPos = fltarr(nFrames,2,nTrVel)
  gaPos  = fltarr(nFrames,2,nSph)
  
  print,'interp vel...'
  for i=0,nTrVel-1 do begin
    velPos[*,0,i] = hermite(pos.times_arepo,pos.vel[*,axes[0],i],frameTimes)
    velPos[*,1,i] = hermite(pos.times_arepo,pos.vel[*,axes[1],i],frameTimes)
  endfor
  
  print,'interp sph...'
  for i=0,nSph-1 do begin
    gaPos[*,0,i] = hermite(pos.times_gadget,pos.ga[*,axes[0],i],frameTimes)
    gaPos[*,1,i] = hermite(pos.times_gadget,pos.ga[*,axes[1],i],frameTimes)
  endfor

  ; move tracers left a boxquarter, move sph right a boxquarter
  velPos[*,axes[0],*] -= zoomSize/4.0
  gaPos[*,axes[0],*]  += zoomSize/4.0
  
  print,'rendering frames...'
  for fn=0,nFrames-1 do begin
  ;fn = 300
    ; determine time of this frame and bracketing snapshots
    time = timeStart + timeStep * fn
    
    ; determine earliest frame number to draw trail to
    fnTrail = fn - trailLenFrame + 1 < fn > 0
    numTrail = fn - fnTrail + 1
    
    print,' ['+string(fn,format='(i3)')+'] '+string(time*1000,format='(i4)')+' Myr ( '+str(fn)+' - '+str(fnTrail)+' )'
    
    ; plot
    start_PS, sPa.plotPath + 'frame_'+str(fn)+'.eps', xs=8, ys=4
    
      !p.thick = 1.0
      
      xrange = [boxCen[0]-zoomSize/2.0,boxCen[0]+zoomSize/2.0] ;tile x2
      yrange = [boxCen[1],boxCen[1]+zoomSize/2.0] ; top half only
    
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xstyle=5,ystyle=5,$
             xtitle="",ytitle="",title="",position=[0,0,1,1],xtickname=replicate(' ',10),$
             ytickname=replicate(' ',10)
           
      ; left: draw rvir and rvir/10
      tvcircle,haloRvir/10,sPa.boxSize/2.0-zoomSize/4.0,sPa.boxSize/2.0,$
        cgColor('green'),thick=0.5,/data
      tvcircle,haloRvir/100,sPa.boxSize/2.0-zoomSize/4.0,sPa.boxSize/2.0,$
        cgColor('green'),thick=0.5,/data
      tvcircle,haloRvir,sPa.boxSize/2.0-zoomSize/4.0,sPa.boxSize/2.0,$
        cgColor('green'),thick=0.5,/data
        
      tvcircle,haloRvir/10,sPa.boxSize/2.0+zoomSize/4.0,sPa.boxSize/2.0,$
        cgColor('green'),thick=0.5,/data
      tvcircle,haloRvir/100,sPa.boxSize/2.0+zoomSize/4.0,sPa.boxSize/2.0,$
        cgColor('green'),thick=0.5,/data
      tvcircle,haloRvir,sPa.boxSize/2.0+zoomSize/4.0,sPa.boxSize/2.0,$
        cgColor('green'),thick=0.5,/data 

      ; to plot the trail, plots each particle path separately
      for i=0,nTrVel-1 do $
        for j=fnTrail,fn-1 do $
          plots,velPos[j:(j+1),0,i],velPos[j:(j+1),1,i],color=trailColors[j-fnTrail]
    
      for i=0,nSph-1 do $
        for j=fnTrail,fn-1 do $
          plots,gaPos[j:(j+1),0,i],gaPos[j:(j+1),1,i],color=trailColors[j-fnTrail]
    
      ; draw time ticker and names
      strTime = string(time*1000,format='(i4)')+' Myr'
      cgText,0.5,0.03,strTime,alignment=0.5,color=cgColor('purple'),/normal,charsize=!p.charsize-0.4
      cgText,0.05,0.03,"Arepo",alignment=0.5,color=cgColor('dark orchid'),/normal,charsize=!p.charsize-0.4
      cgText,0.95,0.03,"Gadget",alignment=0.5,color=cgColor('dark orchid'),/normal,charsize=!p.charsize-0.4
      
    end_PS, pngResize=60, im_options='-negate', /deletePS

  endfor

end
