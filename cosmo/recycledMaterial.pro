; recycledMaterial.pro
; feedback project - analysis of recycled material
; dnelson jun.2013

; tracksFluid(): record tracks in (temp,ent,dens) space back in time for each snapshot for each gas element/tracer

function tracksFluid, sP=sP, loadAllTrGal=loadAllTrGal, loadAllTrGmem=loadAllTrGmem, loadAllTrStars=loadAllTrStars

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
  
  ; set saveFilename and check for existence
  if keyword_set(loadAllTrGal) or keyword_set(loadAllTrGmem) or keyword_set(loadAllTrStars) then begin
    if ((keyword_set(loadAllTrGas) or keyword_set(loadAllTrGmem) or keyword_set(loadAllTrStars)) and $
      sP.trMCPerCell eq 0) then message,'Cannot load all tracers for SPH type.'
    
    ; load the results for all MC tracers (galaxy)
    if keyword_set(loadAllTrGal) then begin
      saveFilename = sP.derivPath + 'tracksFluid.'+sP.saveTag+'.gal.'+sP.savPrefix+str(sP.res)+'.'+$
                     str(minSnap)+'-'+str(maxSnap)+'.sav'
                     
      if not file_test(saveFilename) then message,'Error: Specified maxTemps gal not found!'
      restore, saveFilename
      return, rtr_gal
    endif
    
    ; load the results for all MC tracers (groupmem)
    if keyword_set(loadAllTrGmem) then begin
        saveFilename = sP.derivPath + 'tracksFluid.'+sP.saveTag+'.gmem.'+sP.savPrefix+str(sP.res)+'.'+$
                        str(minSnap)+'-'+str(maxSnap)+'.sav'
                     
      if not file_test(saveFilename) then message,'Error: Specified maxTemps gmem not found!'
      restore, saveFilename
      return, rtr_gmem
    endif
    
    ; load the results for all MC tracers (stars)
    if keyword_set(loadAllTrStars) then begin
        saveFilename = sP.derivPath + 'tracksFluid.'+sP.saveTag+'.stars.'+sP.savPrefix+str(sP.res)+'.'+$
                        str(minSnap)+'-'+str(maxSnap)+'.sav'
                     
      if not file_test(saveFilename) then message,'Error: Specified maxTemps stars not found!'
      restore, saveFilename
      return, rtr_stars
    endif
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
      
      ; locate tracer children (indices) of gas id subsets
      tr_parids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='parentids')
      galcat_gal_trids   = cosmoTracerChildren(sP=sP, /getInds, gasIDs=galcat.galaxyIDs, tr_parids=tr_parids)
      galcat_gmem_trids  = cosmoTracerChildren(sP=sP, /getInds, gasIDs=galcat.groupmemIDs, tr_parids=tr_parids)
      galcat_stars_trids = cosmoTracerChildren(sP=sP, /getInds, starIDs=galcat.stellarIDs, tr_parids=tr_parids)
      tr_parids = !NULL
      
      ; convert tracer children indices to tracer IDs at this zMin
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
      
      galcat_gal_trids   = tr_ids[galcat_gal_trids]
      galcat_gmem_trids  = tr_ids[galcat_gmem_trids]
      galcat_stars_trids = tr_ids[galcat_stars_trids]
  
      tr_ids    = !NULL
      galcat    = !NULL ; not used past this point
      
      rr = { snaps : lonarr(nSnaps), times : fltarr(nSnaps), redshifts : fltarr(nSnaps) }
      
      ; store the main arrays for all tracers as structures so we can write them directly
      rtr_gal  = { temp : fltarr(nSnaps,n_elements(galcat_gal_trids))  ,$
                   ent  : fltarr(nSnaps,n_elements(galcat_gal_trids))  ,$
                   dens : fltarr(nSnaps,n_elements(galcat_gal_trids))  ,$
                   flag : intarr(nSnaps,n_elements(galcat_gal_trids))  , rr : rr }
             
      rtr_gmem = { temp : fltarr(nSnaps,n_elements(galcat_gmem_trids))  ,$
                   ent  : fltarr(nSnaps,n_elements(galcat_gmem_trids))  ,$
                   dens : fltarr(nSnaps,n_elements(galcat_gmem_trids))  ,$
                   flag : intarr(nSnaps,n_elements(galcat_gmem_trids))  , rr : rr  }
      
      rtr_stars = { temp : fltarr(nSnaps,n_elements(galcat_stars_trids))  ,$
                    ent  : fltarr(nSnaps,n_elements(galcat_stars_trids))  ,$
                    dens : fltarr(nSnaps,n_elements(galcat_stars_trids))  ,$
                    flag : intarr(nSnaps,n_elements(galcat_stars_trids))  , rr : rr }     
                    
      ; for determining flags
      lastTime = 1.0/(1+zStart)
      tr_wc_last = { gal   : intarr(n_elements(galcat_gal_trids))   ,$
                     gmem  : intarr(n_elements(galcat_gmem_trids))  ,$
                     stars : intarr(n_elements(galcat_stars_trids))  }
      
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
        save,rtr_gal,rtr_gmem,rtr_stars,galcat_gal_trids,galcat_gmem_trids,galcat_stars_trids,$
          tr_wc_last,lastTime,rr,m,filename=resFilename
        print,' --- Done! ---'
      endif
      
      ; load tracer ids and match to child ids from zMin
      tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
       
      idIndexMap = getIDIndexMap(tr_ids,minid=minid)
          
      trids_gal_ind   = idIndexMap[galcat_gal_trids-minid]
      trids_gmem_ind  = idIndexMap[galcat_gmem_trids-minid]
      trids_stars_ind = idIndexMap[galcat_stars_trids-minid]
      idIndexMap = !NULL
        
      tr_ids     = !NULL
      galcat_ind = !NULL
        
      ; tracer maximum temperature
      tr_maxtemp = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxtemp')
      tr_maxtemp = mylog10(tr_maxtemp) ; convert temp to log
       
      rtr_gal.temp[m-minSnap,*]   = tr_maxtemp[trids_gal_ind]
      rtr_gmem.temp[m-minSnap,*]  = tr_maxtemp[trids_gmem_ind]
      rtr_stars.temp[m-minSnap,*] = tr_maxtemp[trids_stars_ind]
        
      tr_maxtemp = !NULL
       
      ; tracer maximum entropy, convert to log(cgs)
      tr_maxent = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxent')
      
      tr_maxent = convertTracerEntToCGS(tr_maxent,/log,sP=sP)
        
      rtr_gal.ent[m-minSnap,*]   = tr_maxent[trids_gal_ind]
      rtr_gmem.ent[m-minSnap,*]  = tr_maxent[trids_gmem_ind]
      rtr_stars.ent[m-minSnap,*] = tr_maxent[trids_stars_ind] ; zero if still in star, or in wind
        
      tr_maxent = !NULL
        
      ; tracer maximum density
      tr_maxdens = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_maxdens')
      tr_maxdens = mylog10(tr_maxdens) ; convert dens to log
      
      rtr_gal.dens[m-minSnap,*]   = tr_maxdens[trids_gal_ind]
      rtr_gmem.dens[m-minSnap,*]  = tr_maxdens[trids_gmem_ind]
      rtr_stars.dens[m-minSnap,*] = tr_maxdens[trids_stars_ind] ; zero if still in star, or in wind
        
      tr_maxdens = !NULL
        
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
      rtr_gal.flag[m-minSnap,*]   = flags[trids_gal_ind]
      rtr_gmem.flag[m-minSnap,*]  = flags[trids_gmem_ind]
      rtr_stars.flag[m-minSnap,*] = flags[trids_stars_ind]
      
      lastTime = h.time
      
      ; check if wind_counter increased
      tr_lst = !NULL
      
      tr_wc = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracer_windcounter')
      
      tr_wc = { gal   : tr_wc[trids_gal_ind]  ,$
                gmem  : tr_wc[trids_gmem_ind] ,$
                stars : tr_wc[trids_stars_ind] }
                
      w = where(tr_wc.gal gt tr_wc_last.gal, count)
      if count gt 0 then rtr_gal.flag[m-minSnap,w] += 100
      w = where(tr_wc.gmem gt tr_wc_last.gmem, count)
      if count gt 0 then rtr_gmem.flag[m-minSnap,w] += 100
      w = where(tr_wc.stars gt tr_wc_last.stars, count)
      if count gt 0 then rtr_stars.flag[m-minSnap,w] += 100
                
      tr_wc_last = tr_wc
      
      ; SAVE?
      if sP.snap eq maxSnap then begin
        
        ; (1) full tracer information (galaxy members) - set savefilename
        saveFilename = sP.derivPath + 'tracksFluid.'+sP.saveTag+'.gal.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(maxSnap)+'.sav'

        rtr_gal.rr = rr
        save,rtr_gal,filename=saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        
        ; (2) full tracer information (group members) - set savefilename
        saveFilename = sP.derivPath + 'tracksFluid.'+sP.saveTag+'.gmem.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(maxSnap)+'.sav'

        rtr_gmem.rr = rr
        save,rtr_gmem,filename=saveFilename
        print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        
        ; (2) full tracer information (group members) - set savefilename
        saveFilename = sP.derivPath + 'tracksFluid.'+sP.saveTag+'.stars.'+sP.savPrefix+str(sP.res)+'.'+$
                       str(minSnap)+'-'+str(maxSnap)+'.sav'

        rtr_stars.rr = rr
        save,rtr_stars,filename=saveFilename
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

pro windCountHisto
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()

  ; config
  sP = simParams(res=128,run='feedback',redshift=2.0)
  ;sgSelect = 'pri'
  gcIDList = getMatchedIDs(simParams=sP,haloID=314)
  
  xrange = [9.5,12.5]
  yrange = [-0.5,8.5] ; center on integers, start at zero
  
  binSizeMass = 0.25/2

  hsp = [0.007,0.08]*2 ; mass (0.003 10-12), frac
  nc  = 200 ; number of colors (of 255) to use for background 2d histo
  
  ; tracks config
  seed  = 424241L
  nRand = 10
  randType = 'hot' ; all,hot,cold
  
  timeRange = [0.75,3.5] ; z=6-2 with some padding, in age of universe
  tempRange = [3.0,8.0]  ; log K
  entRange  = [3.5,9.0] ; log cgs
  densRange = [-1.0,6.0] ; log rho ratio to crit
  
  ; load windcounter and parent halo mass
  parMass = gcSubsetProp(sP=sP,gcIDList=gcIDList,select=sgSelect,/parMass)
  windc   = gcSubsetProp(sP=sP,gcIDList=gcIDList,select=sgSelect,/curTracerVal,singleTracerField='tracer_windcounter')
  
  ; more load
  gc = loadGroupCat(sP=sP,/skipIDs)
  haloTvir = codeMassToVirTemp(gc.subgroupMass[gcIDList],sP=sP,/log)
  
  tracks = gcSubsetProp(sP=sP,gcIDList=gcIDList,select=sgSelect,/tracksFluid)
  maxt   = gcSubsetProp(sP=sP,gcIDList=gcIDList,select=sgSelect,/maxPastTemp)
  
  ; make a rough hot/cold selection and pick random indices
  if randType eq 'all' then begin
    rnd_inds = floor(randomu(seed,nRand) * n_elements(tracks.gal.temp[0,*]))
    rnd_inds = rnd_inds[ uniq(rnd_inds, sort(rnd_inds)) ]
  endif
  
  if randType eq 'hot' then begin
    w_hot  = where(10.0^maxt.gal / 10.0^haloTvir ge 2.0,count_hot)
    rnd_inds = floor(randomu(seed,nRand) * count_hot)
    rnd_inds = w_hot[ rnd_inds[ uniq(rnd_inds, sort(rnd_inds)) ] ]
  endif
  
  if randType eq 'cold' then begin
    w_cold = where(10.0^maxt.gal / 10.0^haloTvir le 0.5,count_cold)
    rnd_inds = floor(randomu(seed,nRand) * count_cold)
    rnd_inds = w_cold[ rnd_inds[ uniq(rnd_inds, sort(rnd_inds)) ] ]
  endif
  
  ; convert densities into log(rho ratio to crit)
  w = where(tracks.gal.dens ne 0,comp=wc,ncomp=ncomp)
  tracks.gal.dens[w] = alog10( rhoRatioToCrit(10.0^tracks.gal.dens[w], sP=sP) )
  if ncomp gt 0 then tracks.gal.dens[wc] = -10.0 ; very small, off left edge
  
  w = where(tracks.gmem.dens ne 0,comp=wc,ncomp=ncomp)
  tracks.gmem.dens[w] = alog10( rhoRatioToCrit(10.0^tracks.gmem.dens[w], sP=sP) )
  if ncomp gt 0 then tracks.gmem.dens[wc] = -10.0
  
  w = where(tracks.stars.dens ne 0,comp=wc,ncomp=ncomp)
  tracks.stars.dens[w] = alog10( rhoRatioToCrit(10.0^tracks.stars.dens[w], sP=sP) )
  if ncomp gt 0 then tracks.stars.dens[wc] = -10.0

  ; bin fractions into halo mass bins and make median lines
  logMassBins   = [9.5,10.0,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11.0,$
                   11.1,11.25,11.5,11.75,11.9,13.1]
  logMassNBins  = n_elements(logMassBins)-1
  logMassBinCen = 0.5 * (logMassBins + shift(logMassBins,-1))
  logMassBinCen = logMassBinCen[0:-2]

  medFracs = { both_zero  : fltarr(logMassNbins) + !values.f_nan ,$
               gmem_zero  : fltarr(logMassNbins) + !values.f_nan ,$
               both_1plus : fltarr(logMassNbins) + !values.f_nan ,$
               gmem_1plus : fltarr(logMassNbins) + !values.f_nan ,$
               both_2plus : fltarr(logMassNbins) + !values.f_nan ,$
               gmem_2plus : fltarr(logMassNbins) + !values.f_nan  }
               
  ; calculate median values in bins of halo mass
  for i=0,logMassNbins-1 do begin
    w_gal   = where(parMass.gal gt logMassBins[i] and parMass.gal le logMassBins[i+1],count_gal)
    w_gmem  = where(parMass.gmem gt logMassBins[i] and parMass.gmem le logMassBins[i+1],count_gmem)
    w_stars = where(parMass.stars gt logMassBins[i] and parMass.stars le logMassBins[i+1],count_stars)
    
    ; count both
    w = where(windc.gal[w_gal] eq 0,count1)
    w = where(windc.stars[w_stars] eq 0,count2)
    medFracs.both_zero[i] = count1 + count2
    
    w = where(windc.gal[w_gal] ge 1,count1)
    w = where(windc.stars[w_stars] ge 1,count2)
    medFracs.both_1plus[i] = count1 + count2
    
    w = where(windc.gal[w_gal] ge 2,count1)
    w = where(windc.stars[w_stars] ge 2,count2)
    medFracs.both_2plus[i] = count1 + count2
    
    ; count gmem
    w = where(windc.gmem[w_gmem] eq 0,count1)
    medFracs.gmem_zero[i] = count1
    
    w = where(windc.gmem[w_gmem] ge 1,count1)
    medFracs.gmem_1plus[i] = count1
    
    w = where(windc.gmem[w_gmem] ge 2,count1)
    medFracs.gmem_2plus[i] = count1
    
    ; normalize (->mass fractions)
    medFracs.both_zero[i]  /= (count_gal+count_stars)
    medFracs.both_1plus[i] /= (count_gal+count_stars)
    medFracs.both_2plus[i] /= (count_gal+count_stars)
    
    medFracs.gmem_zero[i]  /= (count_gmem)
    medFracs.gmem_1plus[i] /= (count_gmem)
    medFracs.gmem_2plus[i] /= (count_gmem)
  endfor
  
  ; plot (1) - fractions vs halo mass
  colors = ['red','blue','forest green']
  
  start_PS, sP.plotPath + 'recFracs.'+str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=[0.0,1.0],/xs,/ys,$
      xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),ytitle="Mass Fraction"
      
    cgPlot,logMassBinCen,medFracs.both_zero,line=0,color=cgColor(colors[0]),/overplot
    cgPlot,logMassBinCen,medFracs.both_1plus,line=0,color=cgColor(colors[1]),/overplot
    cgPlot,logMassBinCen,medFracs.both_2plus,line=0,color=cgColor(colors[2]),/overplot
    
    cgPlot,logMassBinCen,medFracs.gmem_zero,line=1,color=cgColor(colors[0]),/overplot
    cgPlot,logMassBinCen,medFracs.gmem_1plus,line=1,color=cgColor(colors[1]),/overplot
    cgPlot,logMassBinCen,medFracs.gmem_2plus,line=1,color=cgColor(colors[2]),/overplot
     
    ; legend
    legend,['zero','1+','2+'],textcolors=colors,/top,/left
    legend,['central galaxy','halo atmosphere'],linestyle=[0,1],/top,/right
     
  end_PS
  
  ; plot (2) - 2d histo
  start_PS, sP.plotPath + 'rec2DHisto.'+str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
  
    ; both (gal+stars)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange-[0.05,0.0],/xs,/ys,$
      xtitle="",ytitle="Wind Counter",xtickname=replicate(' ',10),$
      yticks=ceil(max(yrange))+1,ytickv=indgen(ceil(max(yrange))),yminor=0,pos=(sP.pos_2x1)[0]
      
    xx = [parMass.gal,parMass.stars]
    yy = [windc.gal,windc.stars]
      
    f2d = binHisto2D(xx=xx, yy=yy, xmm=xrange, ymm=yrange, xbs=binSizeMass, ybs=1.0)
    oplot2DHistoSq, f2d, hsp=hsp, nc=nc, xrange=xrange, yrange=yrange, /colNorm, /blue
    
    cgText,xrange[0]*1.01,yrange[1]*0.84,"Central Galaxy",alignment=0.0
    
    ; redo plot borders
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange-[0.05,0.0],/xs,/ys,xtitle="",ytitle="",/noerase,$
      yticks=ceil(max(yrange))+1,ytickv=indgen(ceil(max(yrange))),yminor=0,pos=(sP.pos_2x1)[0],$
      xtickname=replicate(' ',10)
    
    ; gmem
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange-[0.05,0.0],/xs,/ys,$
      xtitle=textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]"),ytitle="Wind Counter",$
      yticks=ceil(max(yrange))+1,ytickv=indgen(ceil(max(yrange))),yminor=0,pos=(sP.pos_2x1)[1],/noerase
      
    f2d = binHisto2D(xx=parMass.gmem, yy=windc.gmem, xmm=xrange, ymm=yrange, xbs=binSizeMass, ybs=1.0)
    oplot2DHistoSq, f2d, hsp=hsp, nc=nc, xrange=xrange, yrange=yrange, /colNorm, /blue
    
    cgText,xrange[0]*1.01,yrange[1]*0.84,"Halo Atmosphere",alignment=0.0
    
    ; redo plot borders
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange-[0.05,0.0],/xs,/ys,xtitle="",ytitle="",/noerase,$
      yticks=ceil(max(yrange))+1,ytickv=indgen(ceil(max(yrange))),yminor=0,pos=(sP.pos_2x1)[1]
    
  end_PS
  
  ; plot (3) - tracks  
  start_PS, sP.plotPath + 'track.temp.'+randType+'.'+str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
  
    cgPlot,[0],[0],/nodata,xrange=timeRange,yrange=tempRange,xs=9,/ys,$
      xtitle="Age of Universe [Gyr]",ytitle="Temp"
      
    redshift_axis, timeRange, tempRange, sP=sP;, dotted=dotted, zTicknames=zTicknames
    
    for i=0,n_elements(rnd_inds)-1 do begin
    yy = tracks.gal.temp[*,rnd_inds[i]]
    ff = tracks.gal.flag[*,rnd_inds[i]]
  
    ; non-zero values
    w = where(yy ne 0,count,comp=wc,ncomp=ncomp)
    cgPlot,tracks.rr.times[w],yy[w],psym='open diamond',color=units.colors[i],/overplot
      
    ; interpolate over zero values
    tt = interpol(yy[w],tracks.rr.times[w],tracks.rr.times)
      
    ; zero values
    if ncomp gt 0 then $
      cgPlot,tracks.rr.times[wc],tt[wc],psym='X',color=units.colors[i],/overplot
      
    ; mark windcounter increases with discrete symbols / vertical lines
    w = where(ff ge 100,count)
    ;if count gt 0 then $
    ;  cgPlot,tracks.rr.times[w],tt[w],psym='filled circle',color=units.colors[i],/overplot
    if count gt 0 then for j=0,count-1 do $
      cgPlot,replicate( tracks.rr.times[w[j]]+0.005*randomu(seed,1) ,2),tempRange*[1.1,0.9],$
        line=0,thick=1.0,color=units.colors[i],/overplot
  
    ; mark maximum past temp
    cgPlot,timeRange[0]*[1.2,1.6],replicate(maxt.gal[rnd_inds[i]], 2),line=2,color=units.colors[i],/overplot
  
    ; debug
    w = where(ff ge 100,count)
    if windc.gal[rnd_inds[i]] ne count then print,'Warning'
  
    endfor
      
  end_PS
  
  ; plot (4)
  start_PS, sP.plotPath + 'track.ent.'+randType+'.'+str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
  
    cgPlot,[0],[0],/nodata,xrange=timeRange,yrange=entRange,xs=9,/ys,$
      xtitle="Age of Universe [Gyr]",ytitle="Entropy"
      
    redshift_axis, timeRange, entRange, sP=sP;, dotted=dotted, zTicknames=zTicknames
    
    for i=0,n_elements(rnd_inds)-1 do begin
    yy = tracks.gal.ent[*,rnd_inds[i]]
    ff = tracks.gal.flag[*,rnd_inds[i]]
  
    ; non-zero values
    w = where(yy ne 0,count,comp=wc,ncomp=ncomp)
    cgPlot,tracks.rr.times[w],yy[w],psym='open diamond',color=units.colors[i],/overplot
      
    ; interpolate over zero values
    tt = interpol(yy[w],tracks.rr.times[w],tracks.rr.times)
      
    ; zero values
    if ncomp gt 0 then $
      cgPlot,tracks.rr.times[wc],tt[wc],psym='X',color=units.colors[i],/overplot
      
    ; mark windcounter increases with discrete symbols / vertical lines
    w = where(ff ge 100,count)
    ;if count gt 0 then $
    ;  cgPlot,tracks.rr.times[w],tt[w],psym='filled circle',color=units.colors[i],/overplot
    if count gt 0 then for j=0,count-1 do $
      cgPlot,replicate( tracks.rr.times[w[j]]+0.005*randomu(seed,1) ,2),entRange*[1.1,0.9],$
        line=0,thick=1.0,color=units.colors[i],/overplot
  
    ; debug
    w = where(ff ge 100,count)
    if windc.gal[rnd_inds[i]] ne count then print,'Warning'
  
    endfor
      
  end_PS
  
  ; plot (5)
  start_PS, sP.plotPath + 'track.dens.'+randType+'.'+str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
  
    cgPlot,[0],[0],/nodata,xrange=timeRange,yrange=densRange,xs=9,/ys,$
      xtitle="Age of Universe [Gyr]",ytitle="Density"
      
    redshift_axis, timeRange, densRange, sP=sP;, dotted=dotted, zTicknames=zTicknames
    
    for i=0,n_elements(rnd_inds)-1 do begin
    yy = tracks.gal.dens[*,rnd_inds[i]]
    ff = tracks.gal.flag[*,rnd_inds[i]]
  
    ; non-zero values
    w = where(yy ne 0,count,comp=wc,ncomp=ncomp)
    cgPlot,tracks.rr.times[w],yy[w],psym='open diamond',color=units.colors[i],/overplot
      
    ; interpolate over zero values
    tt = interpol(yy[w],tracks.rr.times[w],tracks.rr.times)
      
    ; zero values
    if ncomp gt 0 then $
      cgPlot,tracks.rr.times[wc],tt[wc],psym='X',color=units.colors[i],/overplot
      
    ; mark windcounter increases with discrete symbols / vertical lines
    w = where(ff ge 100,count)
    ;if count gt 0 then $
    ;  cgPlot,tracks.rr.times[w],tt[w],psym='filled circle',color=units.colors[i],/overplot
    if count gt 0 then for j=0,count-1 do $
      cgPlot,replicate( tracks.rr.times[w[j]]+0.005*randomu(seed,1) ,2),densRange*[0.9,0.9],$
        line=0,thick=1.0,color=units.colors[i],/overplot
  
    ; debug
    w = where(ff ge 100,count)
    if windc.gal[rnd_inds[i]] ne count then print,'Warning'
  
    endfor
      
  end_PS 

  ; plot (6)
  start_PS, sP.plotPath + 'track.tempdens.'+randType+'.'+str(sP.res)+'.'+sP.plotPrefix+'.'+str(sP.snap)+'.eps', /big
  
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