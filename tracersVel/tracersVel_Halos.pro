; tracersCosmoHalos.pro
; dev for tracer particles related to cosmological boxes and group catalogs / halos
; dnelson feb.2012

; cosmoTracerVel_CompMassFunctions(): compare gas, tracer, and DM FoF mass functions

function cosmoTracerVel_CompMassFunctions, sP=sP, noPlot=noPlot

  ; config
  ;res      = 256
  ;run      = 'dev.tracer.nonrad'
  ;redshift = 3.0
  trPT = partTypeNum('tracerVEL')
  
  ; load group catalog
  ;sP    = simParams(res=res,run=run,redshift=redshift)
  units = getUnits()

  h  = loadSnapshotHeader(sP=sP)
  gc = loadGroupCat(sP=sP,/skipIDs)
  
  ; load gas masses, calculate dm and tr masses
  gas_mass = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
  dm_mass  = h.massTable[1]  
  
  ; calculate halo masses
  hm_gas  = reform(gc.groupMassType[0,*]) * units.UnitMass_in_Msun
  hm_dm   = reform(gc.groupMassType[1,*]) * units.UnitMass_in_Msun
  hm_star = reform(gc.groupMassType[4,*]) * units.UnitMass_in_Msun
  hm_tr   = reform(gc.groupLenType[trPT,*]) * sP.targetGasMass * units.UnitMass_in_Msun
  
  hm_bar  = hm_gas + hm_star
  
  print,'Found: ['+str(n_elements(hm_gas))+'] gas, ['+str(n_elements(hm_dm))+'] dm, ['+$
        str(n_elements(hm_star))+'] stars, ['+str(n_elements(hm_tr))+'] tracers.'
        
  ; plot (1) - scatterplot of hm_gas vs hm_tr, hm_bar vs hm_tr
  if not keyword_set(noPlot) then begin
  start_PS,sP.plotPath+sP.savPrefix+'.'+str(sP.res)+'.gasScatter.snap='+str(sP.snap)+'.eps'
    xrange = [1e8,1e12]
    yrange = [0.0,2.0]
    
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/xlog,$
         xtitle="Halo Mass [DM h"+textoidl("^{-1}")+" M"+textoidl("_{sun}")+"]",$
         ytitle="Tracer Mass / Gas Mass",$
         title=str(sP.res)+textoidl("^3")+" "+sP.run+" z="+$
               string(sP.redshift,format='(f3.1)')+" (FoF Catalog)"
         
    fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot
    fsc_plot,hm_dm,hm_tr/hm_gas,psym=8,symsize=0.4,/overplot,color=getColor(1)
    
  end_PS
    
  if (h.nPartTot[4] gt 0) then begin
    start_PS,sP.plotPath+sP.savPrefix+'.'+str(sP.res)+'starScatter.snap='+str(sP.snap)+'.eps'
      xrange = [1e8,1e12]
      yrange = [0.0,2.0]
      
      fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/xlog,$
           xtitle="Halo Mass [DM h"+textoidl("^{-1}")+" M"+textoidl("_{sun}")+"]",$
           ytitle="Tracer Mass / Total Baryonic Mass",$
           title=str(sP.res)+textoidl("^3")+" "+sP.run+" z="+$
                 string(sP.redshift,format='(f3.1)')+" (FoF Catalog)"
           
      fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot
      fsc_plot,hm_dm,hm_tr/hm_bar,psym=8,symsize=0.4,/overplot,color=getColor(1)
      
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
    hm_bar  = hm_bar[sort(hm_bar)]
    
    hm_tr  = hm_tr[sort(hm_tr)]
    
    ; y-vals (cumulative number count) and normalize by box volume
    y_gas  = reverse(indgen(n_elements(hm_gas)) + 1)    / (h.boxSize/1000)^3.0 ;Mpc
    y_dm   = reverse(indgen(n_elements(hm_dm)) + 1)     / (h.boxSize/1000)^3.0    
    y_tr   = reverse(indgen(n_elements(hm_tr)) + 1) / (h.boxSize/1000)^3.0
    y_star = reverse(indgen(n_elements(hm_star)) + 1)    / (h.boxSize/1000)^3.0
    y_bar  = reverse(indgen(n_elements(hm_bar)) + 1) / (h.boxSize/1000)^3.0
    
  ; plot
  start_PS,sP.plotPath+sP.savPrefix+'.'+str(sP.res)+'.massFuncs.snap='+str(sP.snap)+'.eps'
    xrange = [1e8,max(hm_dm)*1.2]
    yrange = [1e-4,1e0]
    
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,/xlog,$
         xtitle="",xtickname=replicate(' ',10),$
         ytitle="number ("+textoidl("\geq")+" M) [h"+textoidl("^3")+" Mpc"+textoidl("^{-3}")+"]",$
         title=str(sP.res)+textoidl("^3")+" "+sP.run+" z="+$
               string(sP.redshift,format='(f3.1)')+" (FoF catalog)",$
         position=[0.18,0.35,0.9,0.9]
         
    fsc_plot,hm_dm,y_dm,line=0,/overplot,color=getColor(1)
    fsc_plot,hm_gas,y_gas,line=0,/overplot,color=getColor(3)
    fsc_plot,hm_tr,y_tr,line=0,/overplot,color=getColor(7)
    
    if (h.nPartTot[4] gt 0) then begin
      fsc_plot,hm_star,y_star,line=0,/overplot,color=getColor(2)
      fsc_plot,hm_bar,y_bar,line=0,/overplot,color=getColor(4)
    endif
    
    ; legend
    if (h.nPartTot[4] gt 0) then $
      legend,['dm','gas','star','tracer','baryons'],$
             textcolors=getColor([1,3,2,7,4],/name),/bottom,/left,box=0,margin=0.1
    if (h.nPartTot[4] eq 0) then $
      legend,['dm','gas','tracer'],$
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
    tr_res = interpol(y_tr,hm_tr,res_pts)
    
    ; plot
    fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[1.5,1.5],line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[0.5,0.5],line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,res_pts,tr_res/gas_res,line=0,color=getColor(3),/overplot
    
    ; do the same for the total baryons
    res_pts = 10.0^( findgen(nbins+1)/nbins * (11.8-8.0) + 8.0 )
    bar_res = interpol(y_bar,hm_bar,res_pts)
    tr_res = interpol(y_tr,hm_tr,res_pts)
    fsc_plot,res_pts,tr_res/bar_res,line=0,color=getColor(8),/overplot
    
    ; legend
    if (h.nPartTot[4] gt 0) then $
      legend,['gas','gas+stars'],textcolors=getColor([3,8],/name),/right,/top,box=0
    
  end_PS
  endif ; noPlot

  r = {hm_dm:hm_dm,hm_gas:hm_gas,hm_star:hm_star,hm_bar:hm_bar,hm_tr:hm_tr}
  return, r
  
end

; cosmoTracerVel_StackGroupsRad(): stack radial properties of groups in mass bins and save

function cosmoTracerVel_StackGroupsRad, sP=sP, massBinLog=massBinLog, hIDs=hIDs, stars=stars

  ; load snapshot info and group catalog
  h  = loadSnapshotHeader(sP=sP)
  gc = loadGroupCat(sP=sP,/skipIDs)

  ;nbins = 20
  nbins = 10 ;CUSTOM

  if (n_elements(hIDs) gt 0) then begin
    haloIDs = hIDs ; to prevent overwriting input
    
    ; halo selection (manual)
    ;print,'Using ['+str(n_elements(haloIDs))+'] specified halo IDs.'
    haloRadii = gc.group_r_crit200[haloIDs]
    ;haloRadii = fltarr(n_elements(haloIDs))
    ;if tag_exist(gc,'group_r_crit200') then stop ; catch
    
    ; setup binning
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
    
    ;logminmax    = [0.01,2.0] * minmax(haloRadii) > [0.1,100.0]
    ;logminmax[0] = floor(logminmax[0]*10)/10.0
    ;logminmax[1] = ceil(logminmax[1]/100)*100.0
    ;logminmax    = alog10(logminmax)
    logminmax = alog10([1.0,500.0])
    
    saveTag   = 'massbin='+string(massBinLog[0],format='(f4.1)')+"-"+string(massBinLog[1],format='(f4.1)')+$
                '.snap='+str(sP.snap)+'.nBins='+str(nbins)
  endelse
  
  ; save/restore
  saveFilename = sP.derivPath + sP.savPrefix + str(sP.res) + '.stackRadVel.' + saveTag + '.sav'
                 
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
    rho_stars = fltarr(nbins)
    
    ;radBins = [0.0,           logspace(logminmax[0],logminmax[1],nbins)]
    ;midBins = [radBins[0]/2.0,logspace(logminmax[0],logminmax[1],nbins,/mid)]
    ;TODO CUSTOM BINS
    radBins = [0.0,5.0,12.0,25.0,45.0,75.0,100.0,130.0,180.0,300.0,500.0]
    midBins = [2.5,8.5,18.5,35.0,60.0,87.5,115.0,155.0,240.0,400.0]

    ; load positions
    pos_gas   = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
    pos_tr    = loadSnapshotSubset(sP=sP,partType='tracerVel',field='pos')
    pos_dm    = loadSnapshotSubset(sP=sP,partType='dm',field='pos')
    if keyword_set(stars) then $
      pos_stars = loadSnapshotSubset(sP=sP,partType='stars',field='pos')
    
    ; load masses
    gas_mass   = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
    if keyword_set(stars) then $
      stars_mass = loadSnapshotSubset(sP=sP,partType='stars',field='mass')
    dm_mass    = h.massTable[1]
    
    gas_size    = loadSnapshotSubset(sP=sP,partType='gas',field='vol')
    gas_size    = (gas_size * 3.0 / (4*!pi))^(1.0/3.0) ;cellrad [ckpc]
  
    ; load tracer parents and reverse histogram
    tr_par_ind = cosmoTracerVelParents(sP=sP,/getInds)
    child_counts = histogram(tr_par_ind,rev=child_inds)
    
    ; find alternative halo centers via iterative CM fitting
    ;iterCM = groupCenterPosByIterativeCM(sP=sP,gc=gc,haloIDs=haloIDs)
    ; find alternative halo centers via most bound particle ID
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
  
      rad    = periodicDists(haloPos,pos_tr,sP=sP)
      tr_ind = where(rad le 2*10.0^logminmax[1],tr_count)
      rad_tr = rad[tr_ind]
  
      rad    = periodicDists(haloPos,pos_dm,sP=sP)
      dm_ind = where(rad le 2*10.0^logminmax[1],dm_count)
      rad_dm = rad[dm_ind]
  
      if keyword_set(stars) then begin
        rad       = periodicDists(haloPos,pos_stars,sP=sP)
        stars_ind = where(rad le 2*10.0^logminmax[1],stars_count)
        rad_stars = rad[stars_ind]
      endif else begin
        stars_count = 0
      endelse
      
      print,'(' + str(haloID) + ') Found ['+str(gas_count)+'] gas ['+str(tr_count)+'] tracer ['+$
            str(dm_count)+'] dm ['+str(stars_count)+'] stars inside cut.'
                       
      ; subselect gas,stars masses
      gas_mass_sub   = gas_mass[gas_ind]
      gas_size_sub   = gas_size[gas_ind]
      if keyword_set(stars) then $
        stars_mass_sub = stars_mass[stars_ind]
      
      ; do binning 
      for i=0,nbins-1 do begin
        ; gas
        w1 = where(rad_gas gt radBins[i] and rad_gas le radBins[i+1],count1)
        if (count1 gt 0) then begin
          rho_gas[i]  += total(gas_mass_sub[w1])
          size_gas[i] += total(gas_size_sub[w1])
          num_gas[i]  += count1
        endif
        
        ; number of gas with no tracers
        globalGasInds = gas_ind[w1]
        
        numChildren = child_inds[globalGasInds+1]-child_inds[globalGasInds]
                                    
        w = where(numChildren eq 0,numZeroCount)
        w = where(numChildren gt 1,numMultCount)
        num_gas_notr[i] += numZeroCount
        num_gas_motr[i] += numMultCount

        ; stars
        if keyword_set(stars) then begin
          w2 = where(rad_stars gt radBins[i] and rad_stars le radBins[i+1],count2)
          if (count2 gt 0) then $
            rho_stars[i] += total(stars_mass_sub[w2])
        endif
        
        ; dm
        w3 = where(rad_dm gt radBins[i] and rad_dm le radBins[i+1],count3)
        if (count3 gt 0) then $
          rho_dm[i] += dm_mass * count3
        
        ; tracers
        w4 = where(rad_tr gt radBins[i] and rad_tr le radBins[i+1],count4)
        if (count4 gt 0) then begin
          if keyword_set(trMass) then begin
            rho_tr[i] += total(tr_mass_sub[w4])
          endif else begin
            rho_tr[i] += sP.targetGasMass * count4
          endelse
          num_tr[i] += count4
        endif
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
      rho_stars[i] = rho_stars[i] / vol / haloCount * 1e10
      
      ; normalize average cell size
      if (num_gas[i] gt 0) then size_gas[i]  = size_gas[i] / num_gas[i]
    endfor
    
    num_gas      = num_gas / haloCount ;avg
    num_tr       = num_tr / haloCount
    num_gas_notr = num_gas_notr / haloCount
    num_gas_motr = num_gas_motr / haloCount
    
    r = {rho_gas:rho_gas,size_gas:size_gas,num_gas:num_gas,$
         num_gas_notr:num_gas_notr,num_gas_motr:num_gas_motr,$
         rho_tr:rho_tr,num_tr:num_tr,rho_dm:rho_dm,rho_stars:rho_stars,$
         haloIDs:haloIDs,haloRadii:haloRadii,haloCount:haloCount,$
         logminmax:logminmax,radBins:radBins,midBins:midBins,$
         saveTag:saveTag}

    ; save
    save,r,filename=saveFilename
  endelse
  
  return, r

end

; cosmoCompRadProfilesVel(): compare gas, tracer, and DM radial profiles of halos

pro cosmoCompRadProfilesVel, massBinLog=massBinLog, haloIDs=haloIDs, redshift=redshift

  if (n_elements(massBinLog) eq 0 and n_elements(haloIDs) eq 0) then stop
  if not keyword_set(redshift) then stop

  ; config
  resSet = [128] ;[256,128]
  run    = 'tracernew'
  stars  = 1
  
  ; run match
  ;sP1 = simParams(res=resSet[1],run=run,redshift=redshift)
  sP2 = simParams(res=resSet[0],run=run,redshift=redshift) ;for res comparisons
  ;match = findMatchedHalos(sP1=sP1,sP2=sP2)
  
  ; start plot book
  ;start_PS,sP1.plotPath+sP1.savPrefix+'.book.snap='+str(sP1.snap)+'.radProfiles.ps',eps=0,xs=7,ys=9
  massBinsTag = 'massbin='+string(massBinLog[0],format='(f4.1)')+"-"+string(massBinLog[1],format='(f4.1)')
  
  start_PS,sP2.plotPath+sP2.savPrefix+'.trVel.'+massBinsTag+'.snap='+str(sP2.snap)+'.radProfiles.eps',xs=7,ys=9
  
  ;for i=0,288,2 do begin ;288 matched at z=1
  ;for i=0,n_elements(haloIDs)-1 do begin ;non-matched, just lots of single halos
  ;  print,''
  ;  print,i
  ;  hIDs = [ match.matchedInds[match.wMatch[i]], match.wMatch[i] ] ;sP1,sP2 (res[1],res[0]) (128,256)
  
    !p.multi = [0,1,2]
  
    ; plot 1
    ; ------
    foreach res,resSet,k do begin
    ;foreach run,runSet,k do begin
    
      ; load group catalog
      sP = simParams(res=res,run=run,redshift=redshift)
      ;sP = sP1
      gc = loadGroupCat(sP=sP,/skipIDs)
    
      ; load radially stacked results
      rs = cosmoTracerVel_StackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=haloIDs,stars=stars)
      
      ; use individual halos
      ;rs = cosmoTracerVel_StackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=haloIDs[i],stars=stars)
      
      ; use matched halos between resolutions
      ;rs = cosmoTracerVel_StackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=hIDs[k],trMass=trMass,stars=stars)
  
      ; check positions
      ;sfInd = gc.groupFirstSub[hIDs[k]]
      ;print,gc.groupPos[*,hIDs[k]]
      ;print,gc.subgroupCM[*,sfInd]
      ;print,'fof sfcm dist, match dist ',$
      ;  periodicDists(gc.groupPos[*,hIDs[k]],gc.subgroupCM[*,sfInd],sP=sP),match.posDiffs[hIDs[1]]
    
      ; multi plot config
      line  = k
    
      ; start plot
      if (k eq 0) then begin
        ;plotStr = "nonRad"
        plotStr = str(res)+textoidl('^3')
        
        plotTitle = plotStr+" z="+string(redshift,format='(f3.1)')+" ("+$
                    string(massBinLog[0],format='(f4.1)')+" < log(M) < "+$
                    string(massBinLog[1],format='(f4.1)')+") " + str(n_elements(rs.haloIDs)) + " halos"
  
        ;plotTitle = plotStr+" z="+string(redshift,format='(f3.1)')+" halo ID="+str(haloIDs[i])+$
        ;            " log(M)="+string(alog10(gc.groupmass[haloIDs[i]]*1e10),format='(f5.2)')
              
        ;plotTitle = plotStr+" z="+string(redshift,format='(f3.1)')+" matched"+$
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
        ;if tag_exist(gc,'group_r_crit200') then stop ; catch
  
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
        fsc_plot,rs.midBins,rs.rho_gas+rs.rho_stars,line=line,/overplot,color=getColor(8)
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
    
    ; legend (two res one run)
    ;strs = ['gas ','dm ','tracer ','gas ','dm ','tracer '] + $
    ;       [textoidl('128^3'),textoidl('128^3'),textoidl('128^3'), $
    ;        textoidl('256^3'),textoidl('256^3'),textoidl('256^3')]
    ;colors = getColor([1,2,3,1,2,3],/name)
    ;styles = [1,1,1,0,0,0]
    
    ; legend (one res one run)
    strs = ['gas ','dm ','tracer ']+$
           [textoidl(str(res)+'^3'),textoidl(str(res)+'^3'),textoidl(str(res)+'^3')]
    colors = getColor([1,2,3],/name)
    styles = [0,0,0]
    
    ; legend (one res two run)
    ;strs = ['gas non-rad','dm non-rad','tracer non-rad','gas no-PM','dm no-PM','tracer no-PM']
    ;colors = getColor([1,2,3,1,2,3],/name)
    ;styles = [1,1,1,0,0,0]
    
    legend, strs, textcolors=colors, linestyle=styles, $
      /top, /right, box=0, margin=0.25, linesize=0.25, charsize=!p.charsize-0.2
    
    ; plot 2
    ; ------
    foreach res,resSet,k do begin
    
      ; load group catalog
      ;sP = sP1
      sP = simParams(res=res,run=run,redshift=redshift)
    
      ; load radially stacked results
      rs = cosmoTracerVel_StackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=haloIDs,stars=stars)
      
      ; use individual halos
      ;rs = cosmoTracerVel_StackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=haloIDs[i],stars=stars)
      
      ; load radially stacked results (matched)
      ;rs = cosmoTracerVel_StackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=hIDs[k],stars=stars)
  
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
        
        ; plot gas dens
        fsc_plot,rs.midBins,rs.num_gas,psym=-8,/overplot,color=getColor(1)
        fsc_plot,rs.midBins,rs.num_tr,psym=-8,/overplot,color=getColor(3)
      endif else begin
        fsc_plot,rs.midBins,rs.num_gas,line=line,/overplot,color=getColor(1)
        fsc_plot,rs.midBins,rs.num_tr,line=line,/overplot,color=getColor(3)
      endelse
      
      ; plot other densities
      fsc_plot,rs.midBins,rs.size_gas,line=line,/overplot,color=getColor(2)
      fsc_plot,rs.midBins,rs.num_gas_notr/rs.num_gas*100,line=line,/overplot,color=getColor(5)
      fsc_plot,rs.midBins,rs.num_gas_motr/rs.num_gas*100,line=line,/overplot,color=getColor(7)
      
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
      strs = [textoidl('128^3'),textoidl('128^3'),textoidl('128^3'), $
              textoidl('256^3'),textoidl('256^3'),textoidl('256^3')]+$
              [' <N'+textoidl("_{gas}>"),' <r'+textoidl("_{gas cell}")+'>',' <N'+textoidl("_{tr}>"),$
               ' <N'+textoidl("_{gas}>"),' <r'+textoidl("_{gas cell}")+'>',' <N'+textoidl("_{tr}>")]
      colors = getColor([1,2,3,1,2,3],/name)
      styles = [1,1,1,0,0,0]
      legend, strs, textcolors=colors, linestyle=styles, $
        /top, /left, box=0, margin=0.25, linesize=0.25, charsize=!p.charsize-0.3, spacing=!p.charsize+0.2
    
      strs = ['fraction gas w/ N'+textoidl("_{tr}")+'=0 ','fraction gas w/ N'+textoidl("_{tr}")+'>1 ',$
               'fraction gas w/ N'+textoidl("_{tr}")+'=0 ','fraction gas w/ N'+textoidl("_{tr}")+'>1 ']+$
              [textoidl('128^3'),textoidl('128^3'), $
              textoidl('256^3'),textoidl('256^3')]
              
      colors = getColor([5,7,5,7],/name)
      styles = [1,1,0,0]
      legend, strs, textcolors=colors, linestyle=styles, $
        /bottom, /right, box=0, margin=0.25, linesize=0.25, charsize=!p.charsize-0.3, spacing=!p.charsize+0.2
    
    endforeach ;resSet
    
  ;endfor ;match or individuals
  
  ;end_PS ;book

  return
end

; cosmoTracerVel_DiffRadProfiles(): calculate difference between gas and tracer radial profiles at several radii
;                         overplot different resolutions and radii as a function of halo mass bin
;                         plot for multiple resolutions

pro cosmoTracerVel_DiffRadProfiles

  resSet = [256,128]
  run    = 'dev.tracer.nonrad'
  stars  = 0
  
  redshifts = [4.0,3.0,2.0,1.0]
  massBins  = [[11.5,12.0],[11.0,11.5],[10.5,11.0],[10.0,10.5]]
  bins      = [4,9,12] ;~2%,10%,25% of r200

  ; arrays
  ratio_trgas = fltarr(n_elements(bins)+1,n_elements(redshifts),$
                       n_elements(resSet),n_elements(massBins[0,*]))
  
  ; load data
  foreach redshift,redshifts,i do begin
    foreach res,resSet,j do begin
      for k=0,n_elements(massBins[0,*])-1 do begin
        massBinLog = massBins[*,k]
        print,'massBinLog: ',massBinLog
  
        ; load group catalog
        sP = simParams(res=res,run=run,redshift=redshift)
        gc = loadGroupCat(sP=sP,/skipIDs)
      
        ; load radially stacked results
        rs = cosmoTracerVel_StackGroupsRad(sP=sP,massBinLog=massBinLog,hIDs=haloIDs,stars=stars)
    
        ; save ratios at specified bins
        foreach bin,bins,m do begin
          ratio_trgas[m,i,j,k] = rs.rho_tr[bin] / rs.rho_gas[bin]
        endforeach
        
        ; save mean ratio
        w = where(rs.rho_tr ne 0 and rs.rho_gas ne 0)
        ratio_trgas[m,i,j,k] = mean(rs.rho_tr[w] / rs.rho_gas[w])
  
      endfor ;massBin
    endforeach ;res
  endforeach ;redshift
  
  ; plot once for each massBin
  for k=0,n_elements(massBins[0,*])-1 do begin
    start_PS,sP.plotPath+sP.savPrefix+'.massbin='+string(massBins[0,k],format='(f4.1)')+"-"+$
      string(massBins[1,k],format='(f4.1)')+'.radDiff.eps'
  
    xrange = [6.0,0.0]
    yrange = [0.0,5.0]
    
    plotTitle = "non-rad ("+string(massBins[0,k],format='(f4.1)')+$
      " < log("+textoidl("M_{tot}")+") < "+string(massBins[1,k],format='(f4.1)')+")"
    
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
         xtitle="redshift",ytitle=textoidl("\rho_{tr} / \rho_{gas}"),title=plotTitle
  
    fsc_plot,xrange,[1.0,1.0],line=2,color=fsc_color('light gray'),/overplot
  
    binStrs = ['r/r'+textoidl('_{vir}'),'r/r'+textoidl('_{vir}'),'r/r'+textoidl('_{vir}'),''] + $
              ['~0.02','~0.10','~0.25','mean']
    strings = []
    lines   = []
    colors  = []
  
    foreach res,resSet,j do begin
      for m=0,n_elements(bins) do begin
        ratios = ratio_trgas[m,*,j,k] ;vs redshift
        fsc_plot,redshifts,ratios,psym=-8,line=j,color=getColor(m),/overplot
        
        strings = [strings,str(res)+textoidl('^3')+' '+binStrs[m]]
        lines   = [lines,j]
        colors  = [colors,getColor(m,/name)]
      endfor ;m
    endforeach
    
    ; legend
    legend,strings,linestyle=lines,textcolors=colors,box=0,/left,/top,$
      charsize=!p.charsize-0.2,linesize=0.25
  
    end_PS
  endfor
  
  ;endforeach ;redshift
  
end
