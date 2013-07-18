; tracksFluid.pro
; feedback project - time series analysis of history of gas element properties
; dnelson jun.2013

; tracksFluid(): record tracks in (temp,ent,dens,rad) space back in time for each snapshot for each gas element/tracer

function tracksFluid, sP=sP

  forward_function cosmoTracerChildren, cosmoTracerVelParents
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; set minimum snapshot (maxmimum redshift)
  zStart = 6.0
  minSnap = redshiftToSnapnum(zStart,sP=sP)

  ; set maximum snapshot (minimum redshift)
  maxSnap = sP.snap
  
  snapRange = [minSnap,maxSnap]
  nSnaps = maxSnap - minSnap + 1
  
  ; load mtS for radial tracking
  mt = mergerTreeSubset(sP=sP,/verbose)
  
  ; set saveFilename and check for existence
  saveFilename = sP.derivPath + 'tracksFluid.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(minSnap)+'-'+str(maxSnap)+'.sav'
                     
  if file_test(saveFilename) then begin
    restore, saveFilename
    return, rtr
  endif

  resFilename = sP.derivPath + 'tracksFluid.'+sP.saveTag+'.restart.'+sP.savPrefix+str(sP.res)+'.'+$
                        str(minSnap)+'-'+str(maxSnap)+'.sav'
  
  ; load galaxy/group member catalogs at zMin for gas ids to search for
  galcat = galaxyCat(sP=sP) ;sP.snap is still at zMin
  
  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
    message,'todo'
  endif
  
  ; MONTE CARLO TRACERS CASE - for each original gas cell, determine some statistics of its
  ; population of tracers and an estimate for the dispersion in those statistics
  ; ------------------------
  if sP.trMCPerCell gt 0 then begin
  
    print,'Calculating new tracksFluid using ( TracerMC ) res = '+str(sP.res)+$
      ' in range ['+str(minSnap)+'-'+str(maxSnap)+'].'
      
    if ~file_test(resFilename) then begin ; no restart   
      
      ; locate tracer children in mtS at starting snapshot
      gcIndOrigTr = mergerTreeRepParentIDs(mt=mt,galcat=galcat,sP=sP,/compactMtS,trids=galcat_trids)
      
      galcat = !NULL ; not used past this point
      
      rr = { snaps : lonarr(nSnaps), times : fltarr(nSnaps), redshifts : fltarr(nSnaps) }
      
      ; store the main arrays for all tracers as structures so we can write them directly
      rtr  = { temp : fltarr(nSnaps,n_elements(galcat_trids))  ,$
               ent  : fltarr(nSnaps,n_elements(galcat_trids))  ,$
               dens : fltarr(nSnaps,n_elements(galcat_trids))  ,$
               rad  : fltarr(nSnaps,n_elements(galcat_trids))  ,$
               flag : intarr(nSnaps,n_elements(galcat_trids))  , rr : rr }
             
      ; for determining flags
      lastTime = 1.0/(1+zStart)
      tr_wc_last = intarr(n_elements(galcat_trids))
      
    endif else begin
      ; restart
      if ~keyword_set(restart) then message,'Error: Restart file exists but restart not requested.'
      restore,resFilename,/verbose
      snapRange[0] = m
    endelse 

    for m=snapRange[0],snapRange[1],1 do begin
      sP.snap = m
      print,m
      
      h = loadSnapshotHeader(sP=sP)
      
      rr.snaps    [ m-minSnap ] = m
      rr.redshifts[ m-minSnap ] = 1/h.time-1
      rr.times    [ m-minSnap ] = redshiftToAgeFlat(1/h.time-1)
  
      ; save restart?
      if m mod 10 eq 0 and m gt minSnap and keyword_set(restart) then begin
        print,' --- Writing restart! ---'
        save,rtr,galcat_trids,gcIndOrigTr,mts_trids,$
          tr_wc_last,lastTime,rr,m,filename=resFilename
        print,' --- Done! ---'
      endif
      
      ; load tracer ids and match to child ids from zMin
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
       
      idIndexMap = getIDIndexMap(tr_ids,minid=minid)
          
      trids_ind  = idIndexMap[galcat_trids-minid]
      
      idIndexMap = !NULL
      tr_ids     = !NULL
      galcat_ind = !NULL
      
      ; rad: load parentids
      tr_parids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentid')
      tr_parids = tr_parids[trids_ind]
      
      ; --- for each each possible parent particle type, match child tracers and save times ---
      parPartTypes = ['gas','stars']
      
      foreach partType,parPartTypes do begin
        ; if we are outside [mt.minSnap,mt.maxSnap] cannot record rad, so skip
        if m lt mt.minSnap or m gt mt.maxSnap then continue
        
        ; load parent IDs and convert tracer parent IDs -> indices (for those tracers with this partType parent now), also positions
        par_ids = loadSnapshotSubset(sP=sP,partType=partType,field='ids')
        par_pos = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
        
        ; note: tr_parids_gal,gmem,stars are NOT UNIQUE, use a value_locate approach (not match)
        sort_inds = calcSort(par_ids)
        par_ids_sorted = par_ids[sort_inds]
        
        par_ind = value_locate(par_ids_sorted,tr_parids) ; indices to par_ids_sorted
        par_ind = sort_inds[par_ind>0] ; indices to par_ids (>0 removes -1 entries, which are removed next line)
        galcat_gal_ind_inPar = where(par_ids[par_ind] eq tr_parids,count_inPar) ; verify we actually matched the ID
        
        if countGal_inPar gt 0 then begin
          tr_pos_type = par_pos[ *,par_ind[galcat_ind_inPar] ]
          par_ind = !NULL
          
          ; calculate current distance of parent from smoothed halo center position for galaxy members
          rad_pri  = periodicDists( $
            reform(mt.hPos[mt.maxSnap-m,*,gcIndOrigTr[galcat_ind_inPar]]),tr_pos_type,sP=sP)
          
          ; store current radius
          rtr.rad[m-minSnap,galcat_ind_inPar] = rad_pri
        endif
        
        ; free some memory for next load
        par_ids_sorted = !NULL
        tr_pos_type    = !NULL
        sort_inds = !NULL
        par_pos   = !NULL
        rad_pri   = !NULL
        
      endforeach ; partType
        
      ; tracer maximum temperature
      tr_maxtemp = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxtemp')
      tr_maxtemp = mylog10(tr_maxtemp) ; convert temp to log
       
      rtr.temp[m-minSnap,*] = tr_maxtemp[trids_ind]
        
      tr_maxtemp = !NULL
       
      ; tracer maximum entropy, convert to log(cgs)
      tr_maxent = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxent')
      tr_maxent = convertTracerEntToCGS(tr_maxent,/log,sP=sP)
        
      rtr.ent[m-minSnap,*]   = tr_maxent[trids_ind] ; zero if in star, or in wind
        
      tr_maxent = !NULL
        
      ; tracer maximum density
      tr_maxdens = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxdens')
      
      ; convert densities into log(rho ratio to crit)
      w = where(tr_maxdens ne 0,comp=wc,ncomp=ncomp)
      tr_maxdens[w] = alog10( rhoRatioToCrit(10.0^tr_maxdens[w], sP=sP) )
      if ncomp gt 0 then tr_maxdens[wc] = -10.0 ; very small
  
      rtr.dens[m-minSnap,*] = tr_maxdens[trids_ind] ; zero if still in star, or in wind
        
      tr_maxdens = !NULL
        
      if sP.gfmWinds eq 1 then begin
        
        ; note: all fluid quantities will be zero if the tracer was in a star (or wind) particle for the 
        ; entire duration between the two snapshots
      
        ; check if last_star_time falls within this new snapshot interval
        tr_lst = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_laststartime')
                 
        flags = intarr(n_elements(tr_lst))
      
        ; (a) LST = All.TimeMax+1 (2) if currently in star, or All.TimeMax+2 (3) if currently in wind, at this snapshot
        ; (b) LST = All.Time (0<LST<1) if in gas and last moved from a star (stellar mass loss)
        ; (c) LST = -All.Time (-1<LST<0) if in gas and last moved from a wind (recoupling)
        AllTimeMax = 1.0
      
        ; FLAG = 1 (in star between this snap and previous), FLAG = 2 (in wind between this snap and previous)
        ; FLAG += 100 (wind counter increased between this snap and previous)
      
        ; case (a)
        w = where(tr_lst gt AllTimeMax+0.99 and tr_lst le AllTimeMax+1.01,count)
        if count gt 0 then flags[w] += 1
        w = where(tr_lst gt AllTimeMax+1.99 and tr_lst le AllTimeMax+2.01,count)
        if count gt 0 then flags[w] += 2
      
        ; case (b) - if bracketed by previous and current snaps
        w = where(tr_lst ge lastTime and tr_lst le h.time,count)
        if count gt 0 then flags[w] += 1
      
        ; case (c)
        w = where(tr_lst le -lastTime and tr_lst ge -h.time,count)
        if count gt 0 then flags[w] += 2
      
        ; store flags for gal,gmem,stars and update lastTime
        rtr.flag[m-minSnap,*]   = flags[trids_ind]
      
        lastTime = h.time
      
        ; check if wind_counter increased
        tr_lst = !NULL
      
        tr_wc = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_windcounter')
        tr_wc = tr_wc[trids_ind]
                
        w = where(tr_wc gt tr_wc_last, count)
        if count gt 0 then rtr.flag[m-minSnap,w] += 100
                
        tr_wc_last = tr_wc
      
      endif ; sP.gfmWinds
      
      ; SAVE?
      if sP.snap eq maxSnap then begin
        
        ; (1) full tracer information (galaxy members) - set savefilename
        saveFilename = sP.derivPath + 'tracksFluid.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(maxSnap)+'.sav'

        rtr.rr = rr
        save,rtr,filename=saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        
      endif ;save
    endfor ;m  
  endif
  
  ; VELOCITY TRACERS case - will be similar to above since there could be multiple
  ; ---------------------
  if sP.trMCPerCell eq -1 then begin
    message,'todo'
  endif
  
end

; plotFluidTracks():

pro plotTracksSingle, tracks=tracks, val_ind=val_ind, inds=inds, colors=colors, $
                      xlines=xlines, yrange=yrange ; vertical lines along the x-axes
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  windLineLen   = 0.05 ; fraction of plot height
  windLineThick = 6.0
  xlineStyle    = 1
  
  ; loop over each requested individual track
  for i=0,n_elements(inds)-1 do begin
    yy = (tracks.(val_ind))[*,inds[i]]
    ff = tracks.flag[*,inds[i]]
  
    ; non-zero values
    w = where(yy ne 0,count,comp=wc,ncomp=ncomp)
    cgPlot,tracks.rr.times[w],yy[w],psym='open diamond',color=colors[i],/overplot
      
    ; interpolate over zero values
    tt = interpol(yy[w],tracks.rr.times[w],tracks.rr.times)
      
    ; zero values
    if ncomp gt 0 then $
      cgPlot,tracks.rr.times[wc],tt[wc],psym='X',color=colors[i],/overplot
      
    ; mark windcounter increases with vertical lines
    w = where(ff ge 100,count)
    if count eq 0 or ~keyword_set(yrange) then continue
    
    wc_x = tracks.rr.times[w]
      
    foreach wc_xpos,wc_x do begin
      cgPlot,[wc_xpos,wc_xpos],[yrange[1],yrange[1]-windLineLen*(yrange[1]-yrange[0])],$
        line=0,thick=!p.thick+windLineThick,color=colors[i],/overplot
      cgPlot,[wc_xpos,wc_xpos],[yrange[0],yrange[0]+windLineLen*(yrange[1]-yrange[0])],$
        line=0,thick=!p.thick+windLineThick,color=colors[i],/overplot
    endforeach
    
    ; separate xlines requested? plot the one corresponding to this track
    if keyword_set(xlines) then begin
      cgPlot,[xlines[i],xlines[i]],[yrange[1],yrange[1]-windLineLen*1.3*(yrange[1]-yrange[0])],$
        line=xlineStyle,thick=!p.thick+windLineThick,color=colors[i],/overplot
      cgPlot,[xlines[i],xlines[i]],[yrange[0],yrange[0]+windLineLen*1.3*(yrange[1]-yrange[0])],$
        line=xlineStyle,thick=!p.thick+windLineThick,color=colors[i],/overplot
    endif
    
  endfor
  
end  

pro plotFluidTracks
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; config
  sP = simParams(res=256,run='feedback',redshift=0.0)
  ;sgSelect = 'pri'
  gcIDList = getMatchedIDs(simParams=sP,haloID=304)

  ; tracks config
  seed  = 424242L
  nRand = 3
  randType = 'hot' ; all,hot,cold
  accMode  = 'smooth'
  
  windLineLen = 0.05 ; fraction of plot height
  
  endTime = ceil(redshiftToAgeFlat(sP.redshift))
  timeRange = [0.75,endTime] ; z=6-z with some padding, in age of universe
  tempRange = [3.0,7.5]  ; log K
  entRange  = [3.5,9.5]  ; log cgs
  densRange = [-0.5,6.0] ; log rho ratio to crit
  radRange  = [0.5,3.0]  ; log kpc

  ; load
  tracks = gcSubsetProp(sP=sP,gcIDList=gcIDList,select=sgSelect,/tracksFluid)
  ids    = gcSubsetProp(sP=sP,gcIDList=gcIDList,select=sgSelect,/elemIDs)
  
  ; calculate maximum temps and times for galaxy
  maxtemp = max(tracks.gal.temp, maxind, dim=1)
  maxtemp_times =  tracks.rr.times[ reform((array_indices(tracks.gal.temp, maxind))[0,*]) ]
                    
  maxind = !NULL
  
  gc = loadGroupCat(sP=sP,/skipIDs)
  mt = mergerTreeSubset(sP=sP)
  
  ; if we want to choose based on accretion mode, just crossmatch the tracer IDs
  ;ids_mode = gcSubsetProp(sP=sP,select='pri',accMode=accMode,/accretionTimeSubset,/elemIDs)
  ;calcMatch,ids,ids_mode,ind1,ind2,count=countMatch
  mode_mask = intarr(n_elements(ids))+1
  ;mode_mask[ind1] = 1
  
  ; find mt index of this halo, and rvir and tvir history
  mtInd = ( where(mt.galcatIDList eq gcIDList,count) )[0]
  if count eq 0 then message,'Error: haloID not in mtS?'
  
  haloTvir_t = reverse( mt.hVirTemp[*,mtInd] )
  haloRvir_t = mylog10( reverse( mt.hVirRad[*,mtInd] ) )
  halo_t     = reverse( redshiftToAgeFlat(1/mt.times-1) )
  
  ; make a rough hot/cold selection in addition to the mode restriction
  haloTvir = codeMassToVirTemp(gc.subgroupMass[gcIDList],sP=sP,/log)

  if randType eq 'all' then $
    w_select = where(mode_mask eq 1, count_select)
  if randType eq 'hot' then $
    w_select  = where(10.0^maxtemp / 10.0^haloTvir ge 2.0 and mode_mask eq 1,count_select)
  if randType eq 'cold' then $
    w_select = where(10.0^maxtemp / 10.0^haloTvir le 0.5 and mode_mask eq 1,count_select)
    
  ; make random indice selection
  nRand = nRand < count_select ; limit for low res
  rnd_inds = floor(randomu(seed,nRand) * count_select)
  rnd_inds = w_select[ rnd_inds[ uniq(rnd_inds, sort(rnd_inds)) ] ] 
  
  ; convert rad into log(r)
  tracks.rad   = mylog10( tracks.rad )
  
  ; plot (0) - composite 4x1 (rad,ent,dens,temp) vs. time
  plotStr = '.'+randType+'-'+accMode+'.'+str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)
  
  colors = ['red','blue','forest green']
  
  start_PS, sP.plotPath + 'tracks.4panel' + plotStr + '.eps', xs=6, ys=12
  
    pos = plot_pos(rows=4,cols=1)
    
    ; top: radius
    cgPlot,[0],[0],/nodata,xrange=timeRange,yrange=radRange,xs=9,/ys,$
      xtitle="",ytitle="log( R ) [kpc]",pos=pos[0],xtickname=replicate(' ',10)
      
    plotTracksSingle, tracks=tracks, val_ind=3, inds=rnd_inds, colors=colors, $
                      yrange=radRange, xlines=maxtemp_times[rnd_inds]
      
    redshift_axis, timeRange, radRange, sP=sP
      
    ; plot halo r_vir(t)
    w = where(haloRvir_t gt 0)
    cgPlot,halo_t[w],haloRvir_t[w],line=2,/overplot
      
    ; second row: entropy
    cgPlot,[0],[0],/nodata,xrange=timeRange,yrange=entRange,/xs,/ys,$
      xtitle="",ytitle=textoidl("log( S ) [K cm^2]"),xtickname=replicate(' ',10),pos=pos[1],/noerase
      
    plotTracksSingle, tracks=tracks, val_ind=1, inds=rnd_inds, colors=colors, $
                      yrange=entRange, xlines=maxtemp_times[rnd_inds]
      
    ; third row: density
    cgPlot,[0],[0],/nodata,xrange=timeRange,yrange=densRange,/xs,/ys,$
      xtitle="",ytitle=textoidl("log( \rho / \rho_{crit,b,z} )"),$
      xtickname=replicate(' ',10),pos=pos[2],/noerase
      
    plotTracksSingle, tracks=tracks, val_ind=2, inds=rnd_inds, colors=colors, $
                      yrange=densRange, xlines=maxtemp_times[rnd_inds]
      
    ; bottom: temp
    cgPlot,[0],[0],/nodata,xrange=timeRange,yrange=tempRange,xs=9,/ys,$
      xtitle=textoidl("t_{age} [Gyr]"),ytitle="log( T ) [K]",pos=pos[3],/noerase
      
    plotTracksSingle, tracks=tracks, val_ind=0, inds=rnd_inds, colors=colors, $
                      yrange=tempRange, xlines=maxtemp_times[rnd_inds]
    
    ; plot halo T_vir(t)
    w = where(haloTvir_t gt 0)
    cgPlot,halo_t[w],haloTvir_t[w],line=2,/overplot
    
    ; labels
    cgText,mean( (pos[0])[0:2:2] ), (pos[0])[3]+0.028, "Redshift", alignment=0.5, /normal

  end_PS
  
  ; plot (6)
  start_PS, sP.plotPath + 'track.tempdens' + plotStr + '.eps', /big
  
    cgPlot,[0],[0],/nodata,xrange=densRange,yrange=tempRange,/xs,/ys,$
      xtitle="Density",ytitle="Temperature"
      
    for i=0,n_elements(rnd_inds)-1 do begin
      xx = tracks.gal.dens[*,rnd_inds[i]]
      yy = tracks.gal.temp[*,rnd_inds[i]]
      ff = tracks.gal.flag[*,rnd_inds[i]]
  
      ; non-zero values
      w = where(xx ne 0 and yy ne 0,count,comp=wc,ncomp=ncomp)
      cgPlot,xx[w],yy[w],psym='open diamond',color=units.colors[i],/overplot
      
      ; interpolate over zero values
      xx_i = interpol(xx[w],tracks.rr.times[w],tracks.rr.times)
      yy_i = interpol(yy[w],tracks.rr.times[w],tracks.rr.times)
      
      ; zero values with crosses
      if ncomp gt 0 then $
        cgPlot,xx_i[wc],yy_i[wc],psym='X',color=units.colors[i],/overplot
      
      ; progression line
      cgPlot,xx_i,yy_i,line=0,thick=1.2,color=units.colors[i],/overplot
      
      ; mark windcounter increases with discrete symbols
      w = where(ff ge 100,count)
      if count gt 0 then $
        cgPlot,xx_i[w],yy_i[w],psym='open circle',color=units.colors[i],symsize=2.0,/overplot
  
      ; mark final spot
      cgPlot,xx[n_elements(xx)-1],yy[n_elements(xx)-1],psym='filled square',color=units.colors[i],/overplot
  
    endfor
      
  end_PS  

  stop

end

; meanFluidTracks()

pro meanFluidTracks
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; config
  sP = simParams(res=256,run='tracer',redshift=2.0)
  sgSelect = 'pri'

  sK       = 3
  hcSplit  = [0.5,1.5] ; required tmax/tviracc ratio above or below which is hot or cold
  accMode  = 'all' ;['all','smooth','clumpy','stripped','recycled']
  
  timeRange = [0.75,3.5] ; z=6-2 with some padding, in age of universe
  tempRange = [3.0,6.5]  ; log K
  entRange  = [4.5,8.5]  ; log cgs
  densRange = [-0.5,4.5] ; log rho ratio to crit
  radRange  = [0.5,4.0]  ; log kpc

  ; load
  tracks  = gcSubsetProp(sP=sP,select=sgSelect,accMode=accMode,/accretionTimeSubset,/tracksFluid)
  ids     = gcSubsetProp(sP=sP,select=sgSelect,accMode=accMode,/accretionTimeSubset,/elemIDs)
  maxtemp = gcSubsetProp(sP=sP,select=sgSelect,accMode=accMode,/accretionTimeSubset,/maxPastTemp)
  tviracc = gcSubsetProp(sP=sP,select=sgSelect,accMode=accMode,/accretionTimeSubset,/accTvir)

  ; make a hot/cold selection in addition to the mode restriction
  w_hot = { gal   : where(10.0^maxtemp.gal / 10.0^tviracc.gal ge hcSplit[1])     ,$
            gmem  : where(10.0^maxtemp.gmem / 10.0^tviracc.gmem ge hcSplit[1])   ,$
            stars : where(10.0^maxtemp.stars / 10.0^tviracc.stars ge hcSplit[1])  }
            
  w_cold = { gal   : where(10.0^maxtemp.gal / 10.0^tviracc.gal le hcSplit[0])     ,$
             gmem  : where(10.0^maxtemp.gmem / 10.0^tviracc.gmem le hcSplit[0])   ,$
             stars : where(10.0^maxtemp.stars / 10.0^tviracc.stars le hcSplit[0])  } 
  
  ; convert rad into log(r)
  tracks.gal.rad   = mylog10( tracks.gal.rad )
  tracks.gmem.rad  = mylog10( tracks.gmem.rad )
  tracks.stars.rad = mylog10( tracks.stars.rad )
  
  ; replace all zero's in tracks with NaN's so that they do not contribute to mean
  tracks.gal.temp[where(tracks.gal.temp eq 0.0)]     = !values.f_nan
  tracks.gmem.temp[where(tracks.gmem.temp eq 0.0)]   = !values.f_nan
  tracks.stars.temp[where(tracks.stars.temp eq 0.0)] = !values.f_nan
  
  mean_temp = { gal_cold   : mean( tracks.gal.temp[*,w_cold.gal], dim=2, /nan )     ,$
                gal_hot    : mean( tracks.gal.temp[*,w_hot.gal],  dim=2, /nan )     ,$
                gmem_cold  : mean( tracks.gmem.temp[*,w_cold.gmem], dim=2, /nan )   ,$
                gmem_hot   : mean( tracks.gmem.temp[*,w_hot.gmem], dim=2, /nan )    ,$
                stars_cold : mean( tracks.stars.temp[*,w_cold.stars], dim=2, /nan ) ,$
                stars_hot  : mean( tracks.stars.temp[*,w_hot.stars], dim=2, /nan )   }
                
  mean_dens = { gal_cold   : mean( tracks.gal.dens[*,w_cold.gal], dim=2, /nan )     ,$
                gal_hot    : mean( tracks.gal.dens[*,w_hot.gal],  dim=2, /nan )     ,$
                gmem_cold  : mean( tracks.gmem.dens[*,w_cold.gmem], dim=2, /nan )   ,$
                gmem_hot   : mean( tracks.gmem.dens[*,w_hot.gmem], dim=2, /nan )    ,$
                stars_cold : mean( tracks.stars.dens[*,w_cold.stars], dim=2, /nan ) ,$
                stars_hot  : mean( tracks.stars.dens[*,w_hot.stars], dim=2, /nan )   }
                
  mean_ent  = { gal_cold   : mean( tracks.gal.ent[*,w_cold.gal], dim=2, /nan )     ,$
                gal_hot    : mean( tracks.gal.ent[*,w_hot.gal],  dim=2, /nan )     ,$
                gmem_cold  : mean( tracks.gmem.ent[*,w_cold.gmem], dim=2, /nan )   ,$
                gmem_hot   : mean( tracks.gmem.ent[*,w_hot.gmem], dim=2, /nan )    ,$
                stars_cold : mean( tracks.stars.ent[*,w_cold.stars], dim=2, /nan ) ,$
                stars_hot  : mean( tracks.stars.ent[*,w_hot.stars], dim=2, /nan )   }        

  mean_rad  = { gal_cold   : mean( tracks.gal.rad[*,w_cold.gal], dim=2, /nan )     ,$
                gal_hot    : mean( tracks.gal.rad[*,w_hot.gal],  dim=2, /nan )     ,$
                gmem_cold  : mean( tracks.gmem.rad[*,w_cold.gmem], dim=2, /nan )   ,$
                gmem_hot   : mean( tracks.gmem.rad[*,w_hot.gmem], dim=2, /nan )    ,$
                stars_cold : mean( tracks.stars.rad[*,w_cold.stars], dim=2, /nan ) ,$
                stars_hot  : mean( tracks.stars.rad[*,w_hot.stars], dim=2, /nan )   }              
  
  ; plot (0) - composite 4x1 (rad,ent,dens,temp) vs. time
  plotStr = '.'+accMode+'.'+str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)
  
  colors = ['black','red','blue','forest green']
  
  start_PS, sP.plotPath + 'tracksMean.4panel' + plotStr + '.eps', xs=6, ys=12
  
    pos = plot_pos(rows=4,cols=1)
    
    ; top: radius
    cgPlot,[0],[0],/nodata,xrange=timeRange,yrange=radRange,xs=9,/ys,$
      xtitle="",ytitle="log( R ) [kpc]",pos=pos[0],xtickname=replicate(' ',10)
      
    cgPlot,tracks.rr.times[1:*],smooth(mean_rad.gal_cold[1:*],sK),line=0,color=colors[2],/overplot
    cgPlot,tracks.rr.times[1:*],smooth(mean_rad.gal_hot[1:*],sK),line=0,color=colors[1],/overplot
    
    cgPlot,tracks.rr.times[1:*],smooth(mean_rad.gmem_cold[1:*],sK),line=1,color=colors[2],/overplot
    cgPlot,tracks.rr.times[1:*],smooth(mean_rad.gmem_hot[1:*],sK),line=1,color=colors[1],/overplot
    
    redshift_axis, timeRange, radRange, sP=sP
      
    ; second row: entropy
    cgPlot,[0],[0],/nodata,xrange=timeRange,yrange=entRange,/xs,/ys,$
      xtitle="",ytitle=textoidl("log( S ) [K cm^2]"),xtickname=replicate(' ',10),pos=pos[1],/noerase
      
    cgPlot,tracks.rr.times,smooth(mean_ent.gal_cold,sK),line=0,color=colors[2],/overplot
    cgPlot,tracks.rr.times,smooth(mean_ent.gal_hot,sK),line=0,color=colors[1],/overplot
    
    cgPlot,tracks.rr.times,smooth(mean_ent.gmem_cold,sK),line=1,color=colors[2],/overplot
    cgPlot,tracks.rr.times,smooth(mean_ent.gmem_hot,sK),line=1,color=colors[1],/overplot
    
    ; third row: density
    cgPlot,[0],[0],/nodata,xrange=timeRange,yrange=densRange,/xs,/ys,$
      xtitle="",ytitle=textoidl("log( \rho / \rho_{crit,b,z} )"),$
      xtickname=replicate(' ',10),pos=pos[2],/noerase
      
    cgPlot,tracks.rr.times,smooth(mean_dens.gal_cold,sK),line=0,color=colors[2],/overplot
    cgPlot,tracks.rr.times,smooth(mean_dens.gal_hot,sK),line=0,color=colors[1],/overplot
    
    cgPlot,tracks.rr.times,smooth(mean_dens.gmem_cold,sK),line=1,color=colors[2],/overplot
    cgPlot,tracks.rr.times,smooth(mean_dens.gmem_hot,sK),line=1,color=colors[1],/overplot
    
    ; bottom: temp
    cgPlot,[0],[0],/nodata,xrange=timeRange,yrange=tempRange,xs=9,/ys,$
      xtitle=textoidl("t_{age} [Gyr]"),ytitle="log( T ) [K]",pos=pos[3],/noerase
      
    cgPlot,tracks.rr.times,smooth(mean_temp.gal_cold,sK),line=0,color=colors[2],/overplot
    cgPlot,tracks.rr.times,smooth(mean_temp.gal_hot,sK),line=0,color=colors[1],/overplot
    
    cgPlot,tracks.rr.times,smooth(mean_temp.gmem_cold,sK),line=1,color=colors[2],/overplot
    cgPlot,tracks.rr.times,smooth(mean_temp.gmem_hot,sK),line=1,color=colors[1],/overplot
    
    ; labels
    cgText,mean( (pos[0])[0:2:2] ), (pos[0])[3]+0.028, "Redshift", alignment=0.5, /normal

  end_PS
  
  stop
end