; tracersMC_Halos.pro
; dev for MC tracer particles (group catalogs in cosmological boxes, radial profiles, etc)
; dnelson mar.2012

; cosmoCompMassFunctions(): compare gas, tracer, and DM FoF mass functions

function cosmoCompMassFunctions, sP=sP, noPlot=noPlot

  ; config
  ;res      = 128
  ;f        = '10'
  ;run      = 'tracerMC.nonrad' ;ref
  ;redshift = 1.0
  
  ; load group catalog
  ;sP    = simParams(res=res,run=run,redshift=redshift,f=f)
  units = getUnits()

  h  = loadSnapshotHeader(sP=sP)
  gc = loadGroupCat(sP=sP,/verbose,/readIDs,/getSortedIDs)
  
  dm_mass  = h.massTable[1]  
  
  ; calculate halo masses from group catalog
  hm_gas  = reform(gc.groupMassType[0,*]) * units.UnitMass_in_Msun
  hm_dm   = reform(gc.groupMassType[1,*]) * units.UnitMass_in_Msun
  hm_star = reform(gc.groupMassType[4,*]) * units.UnitMass_in_Msun
  
  hm_tr_gas  = fltarr(n_elements(hm_gas))
  hm_tr_star = fltarr(n_elements(hm_gas))
  
  ; load gas and star IDs and make ID->index maps
  gas_ids  = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
  gas_indmap = getIDIndexMap(gas_ids,minid=gas_minid)
  gas_ids = !NULL
  
  star_ids = loadSnapshotSubset(sP=sP,partType='star',field='ids')
  star_indmap = getIDIndexMap(star_ids,minid=star_minid)
  star_ids = !NULL
    
  ; load tracer counts
  gas_numtr  = loadSnapshotSubset(sP=sP,partType='gas',field='numtr')
  star_numtr = loadSnapshotSubset(sP=sP,partType='star',field='numtr')

  ; for each fof group, record the total mass in children tracers of member particles of both types
  for i=0,n_elements(hm_gas)-1 do begin
    if (gc.groupLenType[0,i] gt 0) then begin
      ; calculate gas tracer mass in halos using gas indices
      gas_ids  = gc.IDsSorted[gc.groupOffsetType[0,i]:gc.groupOffsetType[0,i]+gc.groupLenType[0,i]-1]
      gas_inds = gas_indmap[gas_ids-gas_minid]
      numChildrenGas = total(gas_numtr[gas_inds])
      
      hm_tr_gas[i] = numChildrenGas * sP.trMassConst * units.UnitMass_in_Msun
    endif
    
    if (gc.groupLenType[4,i] gt 0) then begin
      ; repeat for star tracers
      star_ids  = gc.IDsSorted[gc.groupOffsetType[4,i]:gc.groupOffsetType[4,i]+gc.groupLenType[4,i]-1]
      star_inds = star_indmap[star_ids-star_minid]
      numChildrenStar = total(star_numtr[star_inds])
      
      hm_tr_star[i] = numChildrenStar * sP.trMassConst * units.UnitMass_in_Msun
    endif
  endfor
  
  ; plot (1) - scatterplot of hm_gas vs hm_tr_gas, hm_star vs hm_tr_star, and hm_bar vs hm_tr_bar
  if not keyword_set(noPlot) then begin
  start_PS,sP.plotPath+sP.savPrefix+'.'+str(sP.res)+'.gasScatter.snap='+str(sP.snap)+'.eps'
    xrange = [5e8,5e12]
    yrange = [0.0,2.0]
    
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/xlog,$
         xtitle="Halo Mass [DM h"+textoidl("^{-1}")+" M"+textoidl("_{sun}")+"]",$
         ytitle="Gas Tracer Mass / Gas Mass",$
         title=str(sP.res)+textoidl("^3")+" "+sP.run+" (f="+str(sP.trMCPerCell)+$
               ") z="+string(sP.redshift,format='(f3.1)')+" (FoF Catalog)"
         
    fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot
    fsc_plot,hm_dm,hm_tr_gas/hm_gas,psym=8,symsize=0.4,/overplot,color=getColor(1)
    
  end_PS
    
  if (h.nPartTot[4] gt 0) then begin
    start_PS,sP.plotPath+sP.savPrefix+'.'+str(sP.res)+'.starScatter.snap='+str(sP.snap)+'.eps'
      xrange = [5e8,5e12]
      yrange = [0.0,2.0]
      
      fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/xlog,$
           xtitle="Halo Mass [DM h"+textoidl("^{-1}")+" M"+textoidl("_{sun}")+"]",$
           ytitle="Star Tracer Mass / Star Mass",$
           title=str(sP.res)+textoidl("^3")+" "+sP.run+" (f="+str(sP.trMCPerCell)+$
                 ") z="+string(sP.redshift,format='(f3.1)')+" (FoF Catalog)"
           
      fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot
      fsc_plot,hm_dm,hm_tr_star/hm_star,psym=8,symsize=0.4,/overplot,color=getColor(1)
      
      ; legend
      strings = ['tr mass/star mass', 'mean','stddev']
      colors  = getColor([1,3,4],/name)
      legend,strings,textcolors=colors,/right,/top,box=0    
    end_PS
  endif

  ; plot (2) - standard mass functions
    
    ; sort ascending
    hm_gas  = hm_gas[sort(hm_gas)]
    hm_dm   = hm_dm[sort(hm_dm)]
    hm_star = hm_star[sort(hm_star)]
    
    hm_tr_gas  = hm_tr_gas[sort(hm_tr_gas)]
    hm_tr_star = hm_tr_star[sort(hm_tr_star)]
    
    ; y-vals (cumulative number count) and normalize by box volume
    y_gas    = reverse(indgen(n_elements(hm_gas)) + 1)    / (h.boxSize/1000)^3.0 ;Mpc
    y_dm     = reverse(indgen(n_elements(hm_dm)) + 1)     / (h.boxSize/1000)^3.0    
    y_tr_gas = reverse(indgen(n_elements(hm_tr_gas)) + 1) / (h.boxSize/1000)^3.0
    
    y_star    = reverse(indgen(n_elements(hm_star)) + 1)    / (h.boxSize/1000)^3.0
    y_tr_star = reverse(indgen(n_elements(hm_tr_star)) + 1) / (h.boxSize/1000)^3.0
    
  ; plot
  start_PS,sP.plotPath+sP.savPrefix+'.'+str(sP.res)+'.massFuncs.snap='+str(sP.snap)+'.eps'
    xrange = [1e8,max(hm_dm)*1.2]
    yrange = [1e-4,1e0]
    
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,/xlog,$
         xtitle="",xtickname=replicate(' ',10),$
         ytitle="number ("+textoidl("\geq")+" M) [h"+textoidl("^3")+" Mpc"+textoidl("^{-3}")+"]",$
         title=str(sP.res)+textoidl("^3")+" "+sP.run+" (f="+str(sP.trMCPerCell)+$
               ") z="+string(sP.redshift,format='(f3.1)')+" (FoF catalog)",$
         position=[0.18,0.35,0.9,0.9]
         
    fsc_plot,hm_dm,y_dm,line=0,/overplot,color=getColor(1)
    fsc_plot,hm_gas,y_gas,line=0,/overplot,color=getColor(3)
    fsc_plot,hm_star,y_star,line=0,/overplot,color=getColor(2)
    
    fsc_plot,hm_tr_gas,y_tr_gas,line=0,/overplot,color=getColor(7)
    fsc_plot,hm_tr_star,y_tr_star,line=0,/overplot,color=getColor(4)
    
    ; legend
    if (h.nPartTot[4] gt 0) then $
      legend,['dm','gas','star','gas tr','star tr'],$
             textcolors=getColor([1,3,2,7,4],/name),/bottom,/left,box=0,margin=0.1
    if (h.nPartTot[4] eq 0) then $
      legend,['dm','gas','gas tr'],$
             textcolors=getColor([1,3,7],/name),/bottom,/left,box=0,margin=0.1
               
    ; gas/tr residual plot
    yrange = [0.3,1.7]
    fsc_plot,[0],[0],/nodata,/noerase,xrange=xrange,yrange=yrange,/xs,/ys,/xlog,$
             xtitle="Mass [h"+textoidl("^{-1}")+" M"+textoidl("_{sun}")+"]",$
             ytitle="ratio",ytickv=[0.5,1.0,1.5],yticks=2,$
             position=[0.18,0.15,0.9,0.35]
             
    ; just interpolate both onto a set of masses then compare
    nbins = 100
    res_pts = 10.0^( findgen(nbins+1)/nbins * (11.8-8.0) + 8.0 )
    gas_res = interpol(y_gas,hm_gas,res_pts)
    tr_gas_res = interpol(y_tr_gas,hm_tr_gas,res_pts)
    
    ; plot
    fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[1.5,1.5],line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[0.5,0.5],line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,res_pts,tr_gas_res/gas_res,line=0,color=getColor(3),/overplot
    
    ; do the same for the stars
    res_pts = 10.0^( findgen(nbins+1)/nbins * (11.8-8.0) + 8.0 )
    star_res = interpol(y_star,hm_star,res_pts)
    tr_star_res = interpol(y_tr_star,hm_tr_star,res_pts)
    fsc_plot,res_pts,tr_star_res/star_res,line=0,color=getColor(8),/overplot
    
    ; legend
    if (h.nPartTot[4] gt 0) then $
      legend,['gas','star'],textcolors=getColor([3,8],/name),/right,/top,box=0
    
  end_PS
  endif ; noPlot
  
  r = {hm_dm:hm_dm,hm_gas:hm_gas,hm_star:hm_star,hm_tr_gas:hm_tr_gas,hm_tr_star:hm_tr_star}
  return, r
  
end

; cosmoCompMassFunctionsMulti()

pro cosmoCompMassFunctionsMulti

  ; config (MC)
  resSet      = [128,128,128,256]
  fSet        = ['1','10','100','10']
  runSet      = ['tracerMC.nonrad','tracerMC.nonrad','tracerMC.nonrad','tracerMC.nonrad']
  redshiftSet = [3.0,3.0,3.0,3.0]
  
  ; config (VEL)
  resSetVel = [128,256]
  runSetVel = ['dev.tracer.nonrad','dev.tracer.nonrad']
  redshiftSetVel = [3.0,3.0]
  
  strings = []
  
  for k=0,n_elements(resSet)-1 do begin
  
    ; load group catalog
    sP  = simParams(res=resSet[k],run=runSet[k],redshift=redshiftSet[k],f=fSet[k])
    hmf = cosmoCompMassFunctions(sP=sP,/noPlot)

    ; plot (1) mean and dispersion binned into halo mass bins
  
    ; histogram and construct mean and dispersion
    hist = histogram(alog10(hmf.hm_dm),binsize=0.2,loc=loc,rev=rev)
    
    gas_mean  = fltarr(n_elements(loc))
    gas_disp  = fltarr(n_elements(loc))
    star_mean = fltarr(n_elements(loc))
    star_disp = fltarr(n_elements(loc))
    
    for i=0,n_elements(loc)-1 do begin
      if (rev[i] ne rev[i+1]) then begin
        bin_inds = rev[rev[i]:rev[i+1]-1]
        gas_mean[i]  = mean(hmf.hm_tr_gas[bin_inds]/hmf.hm_gas[bin_inds])
        gas_disp[i]  = stddev(hmf.hm_tr_gas[bin_inds]/hmf.hm_gas[bin_inds])
        star_mean[i] = mean(hmf.hm_tr_star[bin_inds]/hmf.hm_star[bin_inds])
        star_disp[i] = stddev(hmf.hm_tr_star[bin_inds]/hmf.hm_star[bin_inds])
      endif
    endfor

    ; start plot on first loop
    if (k eq 0) then begin
      start_PS,sP.plotPath + 'cosmo.multi.gasScatter.mean.snap='+str(sP.snap)+'.eps'
      xrange = [1e8,1e12]
      yrange = [0.8,1.2]
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/xlog,$
           xtitle="Halo Mass [DM h"+textoidl("^{-1}")+" M"+textoidl("_{sun}")+"]",$
           ytitle="Tracer Mass / Baryon Mass",$
           title="Tracer Halo Mass Agreement z="+string(redshiftSet[0],format='(f3.1)')+" (FoF)"

      cgPlot,xrange*[1.3,0.75],[1.0,1.0],line=0,/overplot,color=fsc_color('light gray')
      
      cgPlot,10.0^loc,gas_mean,line=0,/overplot,color=getColor(k)
      ;cgPlot,10.0^loc,star_mean,line=0,/overplot,color=getColor(2*k+1)
      
    endif else begin
      ; on subsequent loops, overplot
      cgPlot,10.0^loc,gas_mean,line=0,/overplot,color=getColor(k)

    endelse ;k=0
    
    strings = [strings,str(resSet[k])+textoidl('^3')+' f='+fSet[k]]
  
  endfor ;nRunsMC
  
  for k=0,n_elements(resSetVel)-1 do begin
  
    ; load group catalog
    sP  = simParams(res=resSetVel[k],run=runSetVel[k],redshift=redshiftSetVel[k])
    hmf = cosmoTracerVel_CompMassFunctions(sP=sP,/noPlot)

    ; plot (1) mean and dispersion binned into halo mass bins

    ; histogram and construct mean and dispersion
    hist = histogram(alog10(hmf.hm_dm),binsize=0.2,loc=loc,rev=rev)
    
    gas_mean = fltarr(n_elements(loc))
    gas_disp = fltarr(n_elements(loc))
    bar_mean = fltarr(n_elements(loc))
    bar_disp = fltarr(n_elements(loc))
    
    for i=0,n_elements(loc)-1 do begin
      if (rev[i] ne rev[i+1]) then begin
        bin_inds = rev[rev[i]:rev[i+1]-1]
        gas_mean[i] = mean(hmf.hm_tr[bin_inds]/hmf.hm_gas[bin_inds])
        gas_disp[i] = stddev(hmf.hm_tr[bin_inds]/hmf.hm_gas[bin_inds])
        bar_mean[i] = mean(hmf.hm_tr[bin_inds]/hmf.hm_bar[bin_inds])
        bar_disp[i] = stddev(hmf.hm_tr[bin_inds]/hmf.hm_bar[bin_inds])
      endif
    endfor

    ; on subsequent loops, overplot
    cgPlot,10.0^loc,gas_mean,line=1,/overplot,color=getColor(n_elements(resSet)+k)

    strings = [strings,str(resSetVel[k])+textoidl('^3')+' tracerVel']
  
  endfor ;nRunsVel

  ; legend
  colors  = getColor(indgen(n_elements(resSet)+n_elements(resSetVel)),/name)
  legend,strings,textcolors=colors,/right,/bottom,box=0
  
  ; plot (2) dispersion
  for k=0,n_elements(resSet)-1 do begin
  
    ; load group catalog
    sP  = simParams(res=resSet[k],run=runSet[k],redshift=redshiftSet[k],f=fSet[k])
    hmf = cosmoCompMassFunctions(sP=sP,/noPlot)

    ; plot (1) mean and dispersion binned into halo mass bins
  
    ; histogram and construct mean and dispersion
    hist = histogram(alog10(hmf.hm_dm),binsize=0.2,loc=loc,rev=rev)
    
    gas_mean  = fltarr(n_elements(loc))
    gas_disp  = fltarr(n_elements(loc))
    star_mean = fltarr(n_elements(loc))
    star_disp = fltarr(n_elements(loc))
    
    for i=0,n_elements(loc)-1 do begin
      if (rev[i] ne rev[i+1]) then begin
        bin_inds = rev[rev[i]:rev[i+1]-1]
        gas_mean[i]  = mean(hmf.hm_tr_gas[bin_inds]/hmf.hm_gas[bin_inds])
        gas_disp[i]  = stddev(hmf.hm_tr_gas[bin_inds]/hmf.hm_gas[bin_inds])
        star_mean[i] = mean(hmf.hm_tr_star[bin_inds]/hmf.hm_star[bin_inds])
        star_disp[i] = stddev(hmf.hm_tr_star[bin_inds]/hmf.hm_star[bin_inds])
      endif
    endfor
  
    ; start plot on first loop
    if (k eq 0) then begin
      start_PS,sP.plotPath + 'cosmo.multi.gasScatter.disp.snap='+str(sP.snap)+'.eps'
      xrange = [1e8,1e12]
      yrange = [0.001,0.5]
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/xlog,/ylog,$
           xtitle="Halo Mass [DM h"+textoidl("^{-1}")+" M"+textoidl("_{sun}")+"]",$
           ytitle="Tracer Dispersion"
      
      cgPlot,10.0^loc,gas_disp,line=0,/overplot,color=getColor(k)
      ;cgPlot,10.0^loc,star_disp,line=0,/overplot,color=getColor(2*k+1)
      print,linfit(loc,alog10(gas_disp))
      
    endif else begin
      cgPlot,10.0^loc,gas_disp,line=0,/overplot,color=getColor(k)
      print,linfit(loc,alog10(gas_disp))
    endelse ;k=0  

  endfor ;nRunsMC
  
  for k=0,n_elements(resSetVel)-1 do begin
  
    ; load group catalog
    sP  = simParams(res=resSetVel[k],run=runSetVel[k],redshift=redshiftSetVel[k])
    hmf = cosmoTracerVel_CompMassFunctions(sP=sP,/noPlot)

    ; plot (1) mean and dispersion binned into halo mass bins
  
    ; histogram and construct mean and dispersion
    hist = histogram(alog10(hmf.hm_dm),binsize=0.2,loc=loc,rev=rev)
    
    gas_mean = fltarr(n_elements(loc))
    gas_disp = fltarr(n_elements(loc))
    bar_mean = fltarr(n_elements(loc))
    bar_disp = fltarr(n_elements(loc))
    
    for i=0,n_elements(loc)-1 do begin
      if (rev[i] ne rev[i+1]) then begin
        bin_inds = rev[rev[i]:rev[i+1]-1]
        gas_mean[i]  = mean(hmf.hm_tr[bin_inds]/hmf.hm_gas[bin_inds])
        gas_disp[i]  = stddev(hmf.hm_tr[bin_inds]/hmf.hm_gas[bin_inds])
        bar_mean[i] = mean(hmf.hm_tr[bin_inds]/hmf.hm_bar[bin_inds])
        bar_disp[i] = stddev(hmf.hm_tr[bin_inds]/hmf.hm_bar[bin_inds])
      endif
    endfor

    ; on subsequent loops, overplot
    cgPlot,10.0^loc,gas_disp,line=1,/overplot,color=getColor(n_elements(resSet)+k)
    print,linfit(loc,alog10(gas_disp))
  endfor ;nRunsVel
  
  ; legend
  colors  = getColor(indgen(n_elements(resSet)+n_elements(resSetVel)),/name)
  legend,strings,textcolors=colors,/bottom,/left,box=0  
  
  end_PS
stop
end

; cosmoStackGroupsRad(): stack radial properties of groups in mass bins and save

function cosmoStackGroupsRad, sP=sP, massBinLog=massBinLog, hIDs=hIDs, stars=stars

  ; load snapshot info and group catalog
  h  = loadSnapshotHeader(sP=sP)
  gc = loadGroupCat(sP=sP,/skipIDs)

  nbins = 10 ;CUSTOM

  if (n_elements(hIDs) gt 0) then begin
    ; halo selection (manual)
    haloIDs = hIDs ; to prevent overwriting input
    
    ; setup binning
    haloRadii = gc.group_r_crit200[haloIDs]
    logminmax = alog10([1.0,500.0])
    
    saveTag = 'halo='+str(haloIDs[0])+'.'+str(n_elements(haloIDs))+$
              '.snap='+str(sP.snap)+'.nBins='+str(nbins)
    
  endif else begin
    ; halo selection (mass bin)
    massBin = 10.0^massBinLog / 1e10
    haloIDs = where(gc.groupMass ge massBin[0] and gc.groupMass lt massBin[1],count)
    print,'Found ['+str(n_elements(haloIDs))+'] halos in mass bin.'
    if (count eq 0) then stop
    
    ; setup binning
    haloRadii = gc.group_r_crit200[haloIDs]
    logminmax = alog10([1.0,500.0])
    
    saveTag   = 'massbin='+string(massBinLog[0],format='(f4.1)')+"-"+string(massBinLog[1],format='(f4.1)')+$
                '.snap='+str(sP.snap)+'.nBins='+str(nbins)
  endelse
 
  ; save/restore
  saveFilename = sP.derivPath + sP.savPrefix + str(sP.res) + '.stackRad.' + saveTag + '.sav'
                 
  if file_test(saveFilename) then begin
    restore,saveFilename
  endif else begin 
    ; arrays
    rho_dm    = fltarr(nbins)
    
    rho_gas      = fltarr(nbins)
    size_gas     = fltarr(nbins)
    num_gas      = fltarr(nbins)
    
    num_gas_notr = fltarr(nbins)
    num_gas_motr = fltarr(nbins)
    
    rho_tr    = fltarr(nbins)
    num_tr    = fltarr(nbins)
    
    rho_stars    = fltarr(nbins)
    num_stars    = fltarr(nbins)
    rho_tr_stars = fltarr(nbins)
    num_tr_stars = fltarr(nbins)
    
    num_stars_notr = fltarr(nbins)
    num_stars_motr = fltarr(nbins)
    
    ;radBins = [0.0,           logspace(logminmax[0],logminmax[1],nbins)]
    ;midBins = [radBins[0]/2.0,logspace(logminmax[0],logminmax[1],nbins,/mid)]
    ;TODO CUSTOM BINS
    radBins = [0.0,5.0,12.0,25.0,45.0,75.0,100.0,130.0,180.0,300.0,500.0]
    midBins = [2.5,8.5,18.5,35.0,60.0,87.5,115.0,155.0,240.0,400.0]

    ; load positions
    pos_gas   = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
    pos_dm    = loadSnapshotSubset(sP=sP,partType='dm',field='pos')
    if keyword_set(stars) then $
      pos_stars = loadSnapshotSubset(sP=sP,partType='stars',field='pos')
    
    ; load masses
    gas_mass   = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
    if keyword_set(stars) then $
      stars_mass = loadSnapshotSubset(sP=sP,partType='stars',field='mass')
    dm_mass    = h.massTable[1]
    
    ; load tracer counts
    gas_numtr = loadSnapshotSubset(sP=sP,partType='gas',field='numtr')
    if keyword_set(stars) then $
      stars_numtr = loadSnapshotSubset(sP=sP,partType='stars',field='numtr')
    
    gas_size    = loadSnapshotSubset(sP=sP,partType='gas',field='vol')
    gas_size    = (gas_size * 3.0 / (4*!pi))^(1.0/3.0) ;cellrad [ckpc]
    
    ; find alternative halo centers via iterative CM fitting
    ;iterCM = groupCenterPosByIterativeCM(sP=sP,gc=gc,haloIDs=haloIDs)
    ; use subgroup most bound ID centers (cannot without actually matching FOF catalogs)
    idMBCM = subgroupPosByMostBoundID(sP=sP)
  
    haloCount = 0
  
    ; locate halo
    foreach haloID,haloIDs,j do begin
    
      ;haloPos = gc.groupPos[*,haloID] ;use FoF centers
      ;haloPos = iterCM.iterDM[*,j] ;use iterative DM CoM centers
      ;if (haloPos[0,j] eq -1) then begin & print,' skip' & continue & endif
    
      ; take most bound particle ID position and r200_crit
      ; skip if group has no subgroup to obtain this value from
      gcInd = gcPriChildInd(gc=gc,haloID=haloID)
      if (gcInd eq -1) then begin & print,' skip' & continue & endif
      haloPos = idMBCM[*,gcInd]
  
      ; calculate radii and make radial cut
      rad     = periodicDists(haloPos,pos_gas,sP=sP)
      gas_ind = where(rad le 2*10.0^logminmax[1],gas_count)
      rad_gas = rad[gas_ind]
      
      tr_count = total(gas_numtr[gas_ind])
      
      rad    = periodicDists(haloPos,pos_dm,sP=sP)
      dm_ind = where(rad le 2*10.0^logminmax[1],dm_count)
      rad_dm = rad[dm_ind]
  
      if keyword_set(stars) then begin
        rad       = periodicDists(haloPos,pos_stars,sP=sP)
        stars_ind = where(rad le 2*10.0^logminmax[1],stars_count)
        rad_stars = rad[stars_ind]
        tr_count_stars = total(stars_numtr[stars_ind])
      endif else begin
        stars_count = 0
        tr_count_stars = 0
      endelse
      
      ; subselect gas,stars masses
      gas_mass_sub   = gas_mass[gas_ind]
      gas_numtr_sub  = gas_numtr[gas_ind]
      gas_size_sub   = gas_size[gas_ind]
      if keyword_set(stars) then begin
        stars_mass_sub = stars_mass[stars_ind]
        stars_numtr_sub = stars_numtr[stars_ind]
      endif   
      
      print,'(' + str(haloID) + ') Found ['+str(gas_count)+'] gas ['+str(tr_count)+'] tracer ['+$
            str(dm_count)+'] dm ['+str(stars_count)+'] stars ['+str(tr_count_stars)+'] star tracers, inside cut.'
                       
      ; do binning 
      for i=0,nbins-1 do begin
        ; gas
        w1 = where(rad_gas gt radBins[i] and rad_gas le radBins[i+1],count1)
        if (count1 gt 0) then begin
          rho_gas[i]  += total(gas_mass_sub[w1])
          size_gas[i] += total(gas_size_sub[w1])
          num_gas[i]  += count1

          ; tracers
          rho_tr[i] += sP.trMassConst * total(gas_numtr_sub[w1])
          num_tr[i] += total(gas_numtr_sub[w1])
          
          ; number of gas with no tracers / more than 2x starting tracers
          w = where(gas_numtr_sub[w1] eq 0,numZeroCount)
          w = where(gas_numtr_sub[w1] gt 2.0*sP.trMCPerCell,numMultCount)
          
          num_gas_notr[i] += numZeroCount
          num_gas_motr[i] += numMultCount
        endif
        
        ; stars
        if keyword_set(stars) then begin
          w2 = where(rad_stars gt radBins[i] and rad_stars le radBins[i+1],count2)
          if (count2 gt 0) then begin
            rho_stars[i] += total(stars_mass_sub[w2])
            num_stars[i] += count2
            
            ; tracers
            rho_tr_stars[i] += sP.trMassConst * total(stars_numtr_sub[w2])
            num_tr_stars[i] += total(stars_numtr_sub[w2])
            
            ; number of stars with no tracers / more than 2x starting tracers
            w = where(stars_numtr_sub[w2] eq 0,numZeroCount2)
            w = where(stars_numtr_sub[w2] gt 2.0*sP.trMCPerCell,numMultCount2)
            
            num_stars_notr[i] += numZeroCount2
            num_stars_motr[i] += numMultCount2
          endif
        endif
        
        ; dm
        w3 = where(rad_dm gt radBins[i] and rad_dm le radBins[i+1],count3)
        if (count3 gt 0) then $
          rho_dm[i] += dm_mass * count3
        
      endfor
      
      haloCount += 1
    
    endforeach

    ; normalize stacked profiles
    for i=0,nbins-1 do begin
      ; shell volume normalization, average over number of halos, and convert to Msun
      vol = 4*!pi/3 * (radBins[i+1]^3.0 - radBins[i]^3.0) ;kpc^3
      
      rho_gas[i]   = rho_gas[i]   / vol / haloCount * 1e10
      rho_tr[i]    = rho_tr[i]    / vol / haloCount * 1e10
      rho_dm[i]    = rho_dm[i]    / vol / haloCount * 1e10
      
      rho_stars[i]    = rho_stars[i] / vol / haloCount * 1e10
      rho_tr_stars[i] = rho_tr_stars[i] / vol / haloCount * 1e10
      
      ; normalize average cell size
      if (num_gas[i] gt 0) then size_gas[i]  = size_gas[i] / num_gas[i]
    endfor
    
    num_gas        = num_gas / haloCount ;avg
    num_stars      = num_stars / haloCount
    num_tr         = num_tr / haloCount ;gas
    num_tr_stars   = num_tr_stars / haloCount
    
    num_gas_notr   = num_gas_notr / haloCount
    num_gas_motr   = num_gas_motr / haloCount
    num_stars_notr = num_stars_notr / haloCount
    num_stars_motr = num_stars_motr / haloCount
    
    r = {rho_gas:rho_gas,size_gas:size_gas,num_gas:num_gas,$
         num_gas_notr:num_gas_notr,num_gas_motr:num_gas_motr,$
         rho_tr:rho_tr,num_tr:num_tr,$
         rho_dm:rho_dm,$
         rho_stars:rho_stars,num_tr_stars:num_tr_stars,rho_tr_stars:rho_tr_stars,$
         num_stars:num_stars,num_stars_notr:num_stars_notr,num_stars_motr:num_stars_motr,$
         haloIDs:haloIDs,haloRadii:haloRadii,haloCount:haloCount,$
         logminmax:logminmax,radBins:radBins,midBins:midBins,$
         saveTag:saveTag}

    ; save
    save,r,filename=saveFilename
  endelse
  
  return, r

end

; cosmoCompRadProfiles(): compare gas, tracer, and DM radial profiles of halos

pro cosmoCompRadProfiles, massBinLog=massBinLog, haloIDs=haloIDs, redshift=redshift, snap=snap

  if (n_elements(massBinLog) eq 0 and n_elements(haloIDs) eq 0) then stop
  if ~keyword_set(redshift) and ~keyword_set(snap) then stop

  ; config
  resSet = [128] ;[256,128]
  ;f      = '1'
  run    = 'dev.tracerMC.SPT' ;coolSF.GFM
  stars  = 1
  
  ; A. run match
  ;sP1 = simParams(res=resSet[1],run=run,redshift=redshift)
  sP2 = simParams(res=resSet[0],run=run,redshift=redshift,snap=snap) ;for res comparisons
  ;match = findMatchedHalos(sP1=sP1,sP2=sP2)
  
  ; B. single run only
  ;sP2 = simParams(res=resSet[0],run=run,redshift=redshift,snap=snap) ;f=f
  
  ; A. start plot book
  start_PS,sP2.plotPath+sP2.savPrefix+'.book.snap='+str(sP2.snap)+'.radProfiles.ps',eps=0,xs=7,ys=9
  
  ; B. start single (stacked) plot
  ;massBinsTag = 'massbin='+string(massBinLog[0],format='(f4.1)')+"-"+string(massBinLog[1],format='(f4.1)')
  ;start_PS,sP2.plotPath+sP2.savPrefix+'.'+massBinsTag+'.snap='+str(sP2.snap)+'.radProfiles.eps',xs=7,ys=9
  
  ;for i=0,288,2 do begin ;288 matched at z=1
  for i=0,n_elements(haloIDs)-1 do begin ;non-matched, just lots of single halos
  ;  hIDs = [ match.matchedInds[match.wMatch[i]], match.wMatch[i] ] ;sP1,sP2 (res[1],res[0]) (128,256)
  
    !p.multi = [0,1,2]
  
    ; plot 1
    ; ------
    foreach res,resSet,k do begin
    ;foreach run,runSet,k do begin
    
      ; load group catalog
      sP = simParams(res=res,run=run,redshift=redshift,snap=snap) ;f=f
      gc = loadGroupCat(sP=sP,/skipIDs)
    
      ; A. load radially stacked results
      ;rs = cosmoStackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=haloIDs,stars=stars)
      
      ; B. use individual halos
      rs = cosmoStackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=haloIDs[i],stars=stars)
      
      ; C. use matched halos between resolutions
      ;rs = cosmoStackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=hIDs[k],stars=stars)
  
      ; multi plot config
      line  = k
    
      ; start plot
      if (k eq 0) then begin
        ;plotStr = "nonRad"
        plotStr = str(res)+textoidl('^3')
        
        ;plotTitle = plotStr+" z="+string(sP.redshift,format='(f3.1)')+" ("+$
        ;            string(massBinLog[0],format='(f4.1)')+" < log(M) < "+$
        ;            string(massBinLog[1],format='(f4.1)')+") " + str(n_elements(rs.haloIDs)) + " halos"
  
        plotTitle = plotStr+" z="+string(sP.redshift,format='(f3.1)')+" halo ID="+str(haloIDs[i])+$
                    " log(M)="+string(alog10(gc.groupmass[haloIDs[i]]*1e10),format='(f5.2)')
              
        ;plotTitle = plotStr+" z="+string(sP.redshift,format='(f3.1)')+" matched"+$
        ;            " log(M)="+string(alog10(gc.groupMass[match.wMatch[i]]*1e10),format='(f5.2)')+" dist="+$
        ;            string(match.posDiffs[match.wMatch[i]],format='(f4.1)')
  
        rs.logminmax = alog10([1.0,500.0])
        xrange = 10.0^rs.logminmax
        ;yrange = [min([rho_gas[where(rho_gas ne 0)],rho_tr,rho_dm,rho_stars[where(rho_stars ne 0)]])/2,$
        ;          max([rho_dm,rho_gas,rho_tr,rho_stars])*2]
        yrange = [min([rs.rho_gas[where(rs.rho_gas ne 0)],rs.rho_dm[where(rs.rho_dm ne 0)]])/2,$
                  max([rs.rho_dm,rs.rho_gas])*2]
         
        fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
             xtitle="r [ckpc]",$
             ytitle="mass density [h"+textoidl("^2")+" M"+textoidl("_{sun}")+$
             " ckpc"+textoidl("^{-3}")+"]",title=plotTitle,/ylog,/xlog
        
        ; r200 lines
        r200 = minmax(gc.group_r_crit200[rs.haloIDs])
  
        fsc_plot,[r200[0],r200[0]],yrange*[1.1,0.98],line=2,color=fsc_color('light gray'),/overplot
        fsc_plot,[r200[1],r200[1]],yrange*[1.1,0.98],line=2,color=fsc_color('light gray'),/overplot
        fsc_text,mean(r200)*0.95,yrange[0]*2,textoidl("r_{200}"),alignment=0.5,$
          charsize=!p.charsize-0.3,color=fsc_color('light gray')
        
        ; plot gas dens
        fsc_plot,rs.midBins,psym=-8,rs.rho_gas,/overplot,color=getColor(1)
      endif else begin
        fsc_plot,rs.midBins,rs.rho_gas,line=line,/overplot,color=getColor(1)
      endelse
      
      ; plot other densities
      fsc_plot,rs.midBins,rs.rho_dm,line=line,/overplot,color=getColor(2)
      fsc_plot,rs.midBins,rs.rho_tr,line=line,/overplot,color=getColor(3)
      
      if (stars eq 1) then begin
        fsc_plot,rs.midBins,rs.rho_stars,line=line,/overplot,color=getColor(7)
        fsc_plot,rs.midBins,rs.rho_tr_stars,line=line,/overplot,color=getColor(8)
      endif
      
      ; softening lines
      soft = sP.gravSoft * [1.0,2.8]
      
     fsc_plot,[soft[0],soft[0]],[yrange[0],yrange[0]*10],line=line,/overplot,color=fsc_color('dark gray')
     fsc_plot,[soft[1],soft[1]],[yrange[0],yrange[0]*10],line=line,/overplot,color=fsc_color('dark gray') 
     
      if (k eq 0) then begin
        fsc_text,soft[0]*0.8,yrange[0]*5,textoidl('\epsilon_{grav}'),color=fsc_color('dark gray'),$
          charsize=!p.charsize-0.5,alignment=0.5
        fsc_text,soft[1]*0.75,yrange[0]*5,textoidl('2.8\epsilon_{grav}'),color=fsc_color('dark gray'),$
          charsize=!p.charsize-0.5,alignment=0.5
      endif
               
    endforeach ;resSet
    
    ; A. legend (two res one run)
    ;strs = ['gas ','dm ','tracer ','gas ','dm ','tracer '] + $
    ;       [textoidl('128^3'),textoidl('128^3'),textoidl('128^3'), $
    ;        textoidl('256^3'),textoidl('256^3'),textoidl('256^3')]
    ;colors = getColor([1,2,3,1,2,3],/name)
    ;styles = [1,1,1,0,0,0]
    
    ; B. legend (one res one run)
    strs = ['gas ','dm ','gas tracer ']+$
           [textoidl(str(res)+'^3'),textoidl(str(res)+'^3'),textoidl(str(res)+'^3')]
    colors = getColor([1,2,3],/name)
    styles = [0,0,0]
    
    ; C. legend (one res two run)
    ;strs = ['gas non-rad','dm non-rad','tracer non-rad','gas no-PM','dm no-PM','tracer no-PM']
    ;colors = getColor([1,2,3,1,2,3],/name)
    ;styles = [1,1,1,0,0,0]
    
    if (stars eq 1) then begin
      strs = [strs,'stars '+textoidl(str(res)+'^3'),'star tracer '+textoidl(str(res)+'^3')]
      colors = [colors,getColor([7,8],/name)]
      styles = [styles,0,0]
    endif
    
    legend, strs, textcolors=colors, linestyle=styles, $
      /top, /right, box=0, margin=0.25, linesize=0.25, charsize=!p.charsize-0.2
    
    ; plot 2
    ; ------
    foreach res,resSet,k do begin
    
      ; load group catalog
      sP = simParams(res=res,run=run,redshift=redshift,snap=snap) ;,f=f
    
      ; A. load radially stacked results
      ;rs = cosmoStackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=haloIDs,stars=stars)
      
      ; B. use individual halos
      rs = cosmoStackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=haloIDs[i],stars=stars)
      
      ; C. load radially stacked results (matched)
      ;rs = cosmoStackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=hIDs[k],stars=stars)
  
      ; start plot
      line = k
      
      if (k eq 0) then begin
        rs.logminmax = alog10([1.0,500.0])
        xrange = 10.0^rs.logminmax
        yrange = [0.8,max([rs.num_gas,rs.size_gas])*2]
        
        plotTitle = "gas size and number counts"
        
        fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
             xtitle="r [ckpc]",$
             ytitle="<r"+textoidl("_{gas cell}")+"> [ckpc] / <N"+textoidl("_{gas}>")+" / <N"+textoidl("_{tr}>"),$
             title=plotTitle,/ylog,/xlog
        
        fsc_plot,xrange,[100.0,100.0],line=1,color=fsc_color('light gray'),/overplot
        fsc_plot,xrange,[10.0,10.0],line=1,color=fsc_color('light gray'),/overplot
        
        ; r200 lines
        r200 = minmax(rs.haloRadii)
  
        fsc_plot,[r200[0],r200[0]],yrange*[1.1,0.98],line=2,color=fsc_color('light gray'),/overplot
        fsc_plot,[r200[1],r200[1]],yrange*[1.1,0.98],line=2,color=fsc_color('light gray'),/overplot
        fsc_text,mean(r200)*0.95,yrange[0]*2,textoidl("r_{200}"),alignment=0.5,$
          charsize=!p.charsize-0.3,color=fsc_color('light gray')
        
        ; plot gas/star dens and number of tracers
        fsc_plot,rs.midBins,rs.num_gas,psym=-8,/overplot,color=getColor(1)
        fsc_plot,rs.midBins,rs.num_tr,psym=-8,/overplot,color=getColor(3)
        fsc_plot,rs.midBins,rs.num_stars,psym=-8,/overplot,color=getColor(4)
        fsc_plot,rs.midBins,rs.num_tr_stars,psym=-8,/overplot,color=getColor(6)
      endif else begin
        fsc_plot,rs.midBins,rs.num_gas,line=line,/overplot,color=getColor(1)
        fsc_plot,rs.midBins,rs.num_tr,line=line,/overplot,color=getColor(3)
        fsc_plot,rs.midBins,rs.num_stars,line=line,/overplot,color=getColor(4)
        fsc_plot,rs.midBins,rs.num_tr_stars,line=line,/overplot,color=getColor(6)
      endelse
      
      ; plot other quantities
      fsc_plot,rs.midBins,rs.size_gas,line=line,/overplot,color=getColor(2)
      fsc_plot,rs.midBins,rs.num_gas_notr/rs.num_gas*100,line=line,/overplot,color=getColor(5)
      fsc_plot,rs.midBins,rs.num_gas_motr/rs.num_gas*100,line=line,/overplot,color=getColor(7)
      
      if (stars eq 1) then begin
        fsc_plot,rs.midBins,rs.num_stars_notr/rs.num_stars*100,line=line,/overplot,color=getColor(8)
        fsc_plot,rs.midBins,rs.num_stars_motr/rs.num_stars*100,line=line,/overplot,color=getColor(12)
      endif
      
      ; softening lines
      soft = sP.gravSoft * [1.0,2.8]
      
     fsc_plot,[soft[0],soft[0]],[yrange[0],yrange[0]*6],line=line,/overplot,color=fsc_color('dark gray')
     fsc_plot,[soft[1],soft[1]],[yrange[0],yrange[0]*6],line=line,/overplot,color=fsc_color('dark gray') 
     
      if (k eq 0) then begin
        fsc_text,soft[0]*0.8,yrange[0]*2,textoidl('\epsilon_{grav}'),color=fsc_color('dark gray'),$
          charsize=!p.charsize-0.5,alignment=0.5
        fsc_text,soft[1]*0.75,yrange[0]*2,textoidl('2.8\epsilon_{grav}'),color=fsc_color('dark gray'),$
          charsize=!p.charsize-0.5,alignment=0.5
      endif
               
      ; legend
      
      ; A. two res sets
      ;strs = [textoidl('128^3'),textoidl('128^3'),textoidl('128^3'), $
      ;        textoidl('256^3'),textoidl('256^3'),textoidl('256^3')]+$
      ;        [' <N'+textoidl("_{gas}>"),' <r'+textoidl("_{gas cell}")+'>',' <N'+textoidl("_{tr}>"),$
      ;         ' <N'+textoidl("_{gas}>"),' <r'+textoidl("_{gas cell}")+'>',' <N'+textoidl("_{tr}>")]
      ;colors = getColor([1,2,3,1,2,3],/name)
      ;styles = [1,1,1,0,0,0]
      
      ; B. one res set with stars
      strs = ['<N'+textoidl("_{gas}>"),'<r'+textoidl("_{gas cell}")+'>','<N'+textoidl("_{tr,gas}>"),$
               '<N'+textoidl("_{stars}>"),'<N'+textoidl("_{tr,stars}>")]
      colors = getColor([1,2,3,4,6],/name)
      styles = [0,0,0,0,0]
      legend, strs, textcolors=colors, $ ;linestyle=styles, (only for multiple res)
        /top, /left, box=0, margin=0.0, linesize=0.25, charsize=!p.charsize-0.3, spacing=!p.charsize+0.2
    
      ; A. two res sets
      ;strs = ['fraction gas w/ N'+textoidl("_{tr}")+'=0 ','fraction gas w/ N'+textoidl("_{tr}")+'>2x ',$
      ;         'fraction gas w/ N'+textoidl("_{tr}")+'=0 ','fraction gas w/ N'+textoidl("_{tr}")+'>2x ']+$
      ;        [textoidl('128^3'),textoidl('128^3'), $
      ;        textoidl('256^3'),textoidl('256^3')]
      ;colors = getColor([5,7,5,7],/name)
      ;styles = [1,1,0,0]
      
      ; B. one res set with stars
      strs = ['fraction gas w/ N'+textoidl("_{tr}")+'=0 ','fraction gas w/ N'+textoidl("_{tr}")+'>2x ',$
               'fraction stars w/ N'+textoidl("_{tr}")+'=0 ','fraction stars w/ N'+textoidl("_{tr}")+'>2x ']
      colors = getColor([5,7,8,12],/name)
      styles = [0,0,0,0]
      
      legend, strs, textcolors=colors, $ ;linestyle=styles, (only for multiple res)
        /bottom, /right, box=0, margin=0.1, linesize=0.25, charsize=!p.charsize-0.3, spacing=!p.charsize+0.2
    
    endforeach ;resSet
    
  endfor ;match or individuals
  
  end_PS

end
