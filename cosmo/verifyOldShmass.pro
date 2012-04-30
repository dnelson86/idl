; verifyOldShmass.pro
; gas accretion project - remake old gadget plot of cold fraction vs halo mass that transitions
;                         strongly at M~11.5 for Tcut~5.5 in agreement with Keres+ 05 etc.
; dnelson apr.2012

; using galaxycat
; ---------------

function maxTempOld, sP=sP, eosok=eosok

  ; set saveFilename and check for existence
  maxSnap = sP.snap-1
  minSnap = 0
  
  eosokTag = ''
  if keyword_set(eosok) then eosoktag = '.eosOK'
  
  saveFilename = sP.derivPath + 'maxTemp.SPH.old'+eosoktag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(maxSnap)+'-'+str(minSnap)+'.sav'  
  
  if file_test(saveFilename) then begin
    print,'Loading: '+strmid(saveFilename,strlen(sp.derivPath))
    restore, saveFilename
    return, r
  endif

  ; load galaxy/group member catalogs at zMin for gas ids to search for
  galcat = galaxyCat(sP=sP)

  ; store the main arrays as a structure so we can write them directly
  r = {maxTemps_gal       : fltarr(n_elements(galcat.galaxyIDs))   ,$
       maxTemps_gmem      : fltarr(n_elements(galcat.groupmemIDs))  }

  for m=minSnap,maxSnap,1 do begin
    sP.snap = m
    ; load local group catalog
    h = loadSnapshotHeader(sP=sP)
    
    ; load gas ids and match to catalog
    ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    
    ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs, need this if we want ids[ids_ind], 
    ; temp[ids_ind], etc to be in the same order as the group catalog id list    
    match,galcat.galaxyIDs,ids,galcat_ind,ids_gal_ind,count=countGal,/sort
    ids_gal_ind = ids_gal_ind[sort(galcat_ind)]
    
    match,galcat.groupmemIDs,ids,galcat_ind,ids_gmem_ind,count=countGmem,/sort
    ids_gmem_ind = ids_gmem_ind[sort(galcat_ind)]
    
    ids        = !NULL
    galcat_ind = !NULL
    
    ; load u,nelec to calculate temperatures
    u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
    nelec = loadSnapshotSubset(sP=sP,partType='gas',field='ne')
    
    u_gal  = u[ids_gal_ind]
    u_gmem = u[ids_gmem_ind]
    u = !NULL
    
    nelec_gal  = nelec[ids_gal_ind]
    nelec_gmem = nelec[ids_gmem_ind]
    nelec = !NULL
    
    temp_gal  = convertUtoTemp(u_gal,nelec_gal,/log)
    temp_gmem = convertUtoTemp(u_gmem,nelec_gmem,/log)
    
    if keyword_set(eosok) then begin
      ; replace existing values (effective EOS not considered, all temps included)
      w1 = where(temp_gal gt r.maxTemps_gal,count1)
      if (count1 gt 0) then r.maxTemps_gal[w1]    = temp_gal[w1]
      
      w2 = where(temp_gmem gt r.maxTemps_gmem,count2)
      if (count2 gt 0) then r.maxTemps_gmem[w2]    = temp_gmem[w2]
    
    endif else begin
      ; load gas SFR and select off the effective EOS (SFR must be zero)
      sfr = loadSnapshotSubset(sP=sP,partType='gas',field='sfr')
      
      sfr_gal  = sfr[ids_gal_ind]
      sfr_gmem = sfr[ids_gmem_ind]
      
      sfr = !NULL
      
      ; replace existing values if current snapshot has higher temps (enforce off effective EOS)
      w1 = where(temp_gal gt r.maxTemps_gal and sfr_gal eq 0.0,count1)
      if (count1 gt 0) then r.maxTemps_gal[w1]    = temp_gal[w1]
      
      ; replace existing values if current snapshot has higher temps (group members)
      w2 = where(temp_gmem gt r.maxTemps_gmem and sfr_gmem eq 0.0,count2)
      if (count2 gt 0) then r.maxTemps_gmem[w2]    = temp_gmem[w2]
    endelse
    
    print,' ['+string(m,format='(i3)')+'] gal: '+string(count1,format='(i6)')+$
      ' gmem: '+string(count2,format='(i6)')
    
  endfor

  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))

  return,r
end

function accModeOld, sP=sP

  ; set saveFilename and check for existence
  maxSnap = sP.snap-1
  minSnap = sP.groupCatRange[0]
  
  saveFilename = sP.derivPath + 'accMode.SPH.old.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(maxSnap)+'-'+str(minSnap)+'.sav'  
  
  if file_test(saveFilename) then begin
    restore, saveFilename
    return, r
  endif

  ; load galaxy/group member catalogs at zMin for gas ids to search for
  galcat = galaxyCat(sP=sP)

  ; store the main arrays as a structure so we can write them directly
  r = {accMode_gal       : intarr(n_elements(galcat.galaxyIDs))+1  ,$
       accMode_gmem      : intarr(n_elements(galcat.groupmemIDs))+1  }

  for m=maxSnap,minSnap,-1 do begin
    sP.snap = m
    ; load local group catalog
    h = loadSnapshotHeader(sP=sP)
    gc = loadGroupCat(sP=sP,/readIDs)
    
    ; old method for comparison: mark as non-smooth if we find in the groupcat at any preivous time
    ;curPIDs = gcPIDList(gc=gc,select='all',partType='gas')
    curPIDs = gc.IDs ; including fuzz (old style)
    
    ; gal
    match,galcat.galaxyIDs,curPIDs,galcat_ind,gc_ind,count=count,/sort
    if count gt 0 then r.accMode_gal[galcat_ind] = 2 ;clumpy

    ; gmem
    match,galcat.groupmemIDs,curPIDs,galcat_ind,gc_ind,count=count2,/sort
    if count2 gt 0 then r.accMode_gmem[galcat_ind] = 2 ;clumpy
    
    ; verbose
    w = where(r.accMode_gal eq 1,count_smooth_gal)
    w = where(r.accMode_gmem eq 1,count_smooth_gmem)
    w = where(r.accMode_gal eq 2,count_clumpy_gal)
    w = where(r.accMode_gmem eq 2,count_clumpy_gmem)
    
    print,' ['+string(m,format='(i3)')+'] gal smooth: '+string(count_smooth_gal,format='(i5)')+$
      ' gal clumpy: '+string(count_clumpy_gal,format='(i6)')+' gmem smooth: '+$
      string(count_smooth_gmem,format='(i5)')+' gmem clumpy: '+$
      string(count_clumpy_gmem,format='(i6)')
        
  endfor

  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))

  return,r
end

pro remakeOldShmass

  ; config
  sP = simParams(res=256,run='gadget',redshift=2.0)
  
  minNum  = 6
  xrange  = [10.0,13.0]
  yrange  = [0.0,1.1]
  
  logMassBinSize = 0.1
  
  TcutVal = 5.5 ; log(K)
  eosok   = 1 ; 0=exclude maxtemp on effEOS, 1=all temps included
  
  ; select on accretion mode by accMode=1 (smooth)
  am = accModeOld(sP=sP)

  gal_w  = where(am.accMode_gal eq 1,count_gal)
  gmem_w = where(am.accMode_gmem eq 1,count_gmem)
  
  ; load galaxycat and maxtemps for all
  galcat  = galaxyCat(sP=sP)
  maxTemp = maxTempOld(sP=sP, eosok=eosok)

  ; load group cat for subgroup masses
  gc = loadGroupCat(sP=sP,/skipIDs)
  gcMasses = codeMassToLogMsun(gc.subgroupMass)

  ; structures to store results (Tmax)
  coldCount  = fltarr(n_elements(gcMasses))
  totalCount = fltarr(n_elements(gcMasses))
  
  ; loop over all subgroups of some type
  gcIDList = gcIDList(gc=gc,select='pri')
  
  foreach gcId,gcIDList do begin
    if galcat.galaxyLen[gcId] gt 0 then begin
      ; list of indices of galaxy gas particles in this subgroup
      loc_inds_gal = lindgen(galcat.galaxyLen[gcId]) + galcat.galaxyOff[gcId]
      
      ; select only smooth
      w = where(am.accMode_gal[loc_inds_gal] eq 1,count)
      
      if count gt 0 then begin
        loc_inds_gal = loc_inds_gal[w]
        loc_maxt_gal = maxTemp.maxTemps_gal[loc_inds_gal]
      
        ; count fraction Tmax below each constant temperature threshold
        w = where(loc_maxt_gal le TcutVal,count_below)
        coldCount[gcId]  = count_below
        totalCount[gcId] = n_elements(loc_inds_gal)
      endif
    endif
    
    if galcat.groupmemLen[gcId] gt 0 then begin
      ; list of indices of groupmem gas particles in this subgroup
      loc_inds_gmem = lindgen(galcat.groupmemLen[gcId]) + galcat.groupmemOff[gcId]
      
      ; select only smooth
      w = where(am.accMode_gmem[loc_inds_gmem] eq 1,count2)
      
      if count2 gt 0 then begin      
        loc_inds_gmem = loc_inds_gmem[w]
        loc_maxt_gmem = maxTemp.maxTemps_gmem[loc_inds_gmem]
      
        ; count fraction Tmax below each constant temperature threshold
        w = where(loc_maxt_gmem le TcutVal,count_below2)
        coldCount[gcId]  += count_below2
        totalCount[gcId] += n_elements(loc_inds_gmem)
      endif
      
      print,[galcat.galaxyIDs[loc_inds_gal],galcat.groupmemIDs[loc_inds_gmem]]
      print,[loc_maxt_gal,loc_maxt_gmem]
      print,count+count2
      stop
    endif
    
  endforeach
  
  ; compute cold fraction by halo
  w = where(totalCount ge minNum)
  coldFrac = coldCount[w] / totalCount[w]
  gcMasses = gcMasses[w]
  
  print,totalCount[0:10]
  print,totalCount[gcIDList[0:10]]
  print,gcIDList[0:10]
  stop
  
  ; calculate median line
  massStep = ceil(100.0/sqrt(n_elements(coldFrac)))/10.0
  ;massStep = 0.2
  massBins = (xrange[1]-xrange[0])/massStep
  massXPts = findgen(massBins)/massBins * (xrange[1]-xrange[0]) + xrange[0] + massStep/2.0
  
  medCold = fltarr(massBins) - 1
  stdCold = fltarr(massBins) - 1
  medCold[0:floor(n_elements(medCold)/2.0)] = 1.0 ;set default value high for first half
  
  for i=0,massBins-1 do begin
    w = where(gcMasses ge xrange[0]+i*massStep and gcMasses lt xrange[0]+(i+1)*massStep,count)
    if (count gt 0) then begin
      medCold[i] = median(coldFrac[w])
      stdCold[i] = stddev(coldFrac[w])
    endif
  endfor  
  
  ; plot
  start_PS, sP.plotPath + 'coldFrac.old_'+str(sP.res)+'.eos'+str(eosok)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("log(M_{halo})")+"",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+") OLD"
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; plot raw points
    psym = 4
    symsize = 0.5
    
    cgPlot,gcMasses,coldFrac,psym=psym,symsize=symsize,thick=1.0,/overplot
    
    ; plot median lines
    cgPlot,massXpts,medCold,color=getColor(3),line=0,/overplot
    
  end_PS
  stop
end

; without galaxycat
; -----------------

function accModeOldNoGC, sP=sP

  ; set saveFilename and check for existence
  maxSnap = sP.snap-1
  minSnap = sP.groupCatRange[0]
  
  saveFilename = sP.derivPath + 'accMode.SPH.old.noGC.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(maxSnap)+'-'+str(minSnap)+'.sav'  
  
  if file_test(saveFilename) then begin
    restore, saveFilename
    return, accMode
  endif

  ; load group catalog at zMin for gas ids to search for
  gcOrig = loadGroupCat(sP=sP,/readIDs)
  
  origIDs = gcOrig.IDs ; noGC1: old method w/o galaxycat (includes fuzz)
  ;origIDs = gcPIDList(gc=gcOrig,select='all',partType='all')
  
  ; store the main arrays so we can write them directly
  accMode = intarr(n_elements(gcOrig.IDs))+1

  for m=maxSnap,minSnap,-1 do begin
    sP.snap = m
    ; load local group catalog
    h = loadSnapshotHeader(sP=sP)
    gc = loadGroupCat(sP=sP,/readIDs)
    
    ; old method for comparison: mark as non-smooth if we find in the groupcat at any preivous time
    ;curPIDs = gcPIDList(gc=gc,select='all',partType='gas')
    curPIDs = gc.IDs ; including fuzz (old style)
    
    match,origIDs,curPIDs,gcorig_ind,gccur_ind,count=count,/sort
    if count gt 0 then accMode[gcorig_ind] = 2 ;clumpy

    ; verbose
    w = where(accMode eq 1,count_smooth)
    w = where(accMode eq 2,count_clumpy)
    
    print,' ['+string(m,format='(i3)')+'] smooth: '+string(count_smooth,format='(i7)')+$
      ' clumpy: '+string(count_clumpy,format='(i7)')
        
  endfor

  ; save
  save,accMode,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))

  return,accMode
end

function maxTempOldNoGC, sP=sP, eosok=eosok

  ; set saveFilename and check for existence
  maxSnap = sP.snap-1
  minSnap = 0
  
  eosokTag = ''
  if keyword_set(eosok) then eosoktag = '.eosOK'
  
  saveFilename = sP.derivPath + 'maxTemp.SPH.old.noGC2'+eosoktag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(maxSnap)+'-'+str(minSnap)+'.sav'  
  
  if file_test(saveFilename) then begin
    print,'Loading: '+strmid(saveFilename,strlen(sp.derivPath))
    restore, saveFilename
    return, maxTemp
  endif

  ; load group catalog at zMin for gas ids to search for
  gcOrig = loadGroupCat(sP=sP,/readIDs)

  ;origIDs = gcOrig.IDs ; noGC1: old method w/o galaxycat (includes fuzz)
  origIDs = gcPIDList(gc=gcOrig,select='all',partType='all') ; noGC2

  ; store the main arrays so we can write them directly
  maxTemp = fltarr(n_elements(origIDs))

  for m=minSnap,maxSnap,1 do begin
    sP.snap = m
    ; load local group catalog
    h = loadSnapshotHeader(sP=sP)
    
    ; load gas ids and match to catalog
    ids = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    
    ; IMPORTANT! rearrange ids_ind to be in the order of gcPIDs, need this if we want ids[ids_ind], 
    ; temp[ids_ind], etc to be in the same order as the group catalog id list    
    match,origIDs,ids,gcorig_ind,ids_ind,count=count,/sort
    ids_ind = ids_ind[sort(gcorig_ind)]
    
    ids        = !NULL
    galcat_ind = !NULL
    
    ; load u,nelec to calculate temperatures
    u = loadSnapshotSubset(sP=sP,partType='gas',field='u')
    u = u[ids_ind]
    
    nelec = loadSnapshotSubset(sP=sP,partType='gas',field='ne')
    nelec = nelec[ids_ind]
    
    temps = convertUtoTemp(u,nelec,/log)
    
    if keyword_set(eosok) then begin
      ; replace existing values (effective EOS not considered, all temps included)
      w1 = where(temps gt maxTemp,count1)
      if (count1 gt 0) then maxTemp[w1] = temps[w1]
    
    endif else begin
      ; load gas SFR and select off the effective EOS (SFR must be zero)
      sfr = loadSnapshotSubset(sP=sP,partType='gas',field='sfr')
      sfr = sfr[ids_ind]
      
      ; replace existing values if current snapshot has higher temps (enforce off effective EOS)
      w1 = where(temps gt maxTemp and sfr eq 0.0,count1)
      if (count1 gt 0) then maxTemp[w1] = temps[w1]
    endelse
    
    print,' ['+string(m,format='(i3)')+'] new maxt: '+string(count1,format='(i6)')
    
  endfor

  ; save
  save,maxTemp,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))

  return,maxTemp
end

pro remakeOldShmassNoGC

  ; config
  sP = simParams(res=256,run='gadget',redshift=2.0)
  
  minNum  = 6
  xrange  = [10.0,13.0]
  yrange  = [0.0,1.1]
  
  logMassBinSize = 0.1
  
  TcutVal = 5.5 ; log(K)
  eosok   = 1 ; 0=exclude maxtemp on effEOS, 1=all temps included
  partType = partTypeNum('gas')
  
  ; select on accretion mode by accMode=1 (smooth)
  accMode  = accModeOldNoGC(sP=sP)
  smooth_w = where(accMode eq 1,count_smooth)
  
  ; load maxtemps for all
  maxTemp = maxTempOldNoGC(sP=sP, eosok=eosok)

  ; load group cat for subgroup masses
  gc = loadGroupCat(sP=sP,/readIDs)
  gcMasses = codeMassToLogMsun(gc.subgroupMass)

  ; structures to store results (Tmax)
  coldCount  = fltarr(n_elements(gcMasses))
  totalCount = fltarr(n_elements(gcMasses))
  
  ; loop over all subgroups of some type
  gcIDList = gcIDList(gc=gc,select='pri')
  
  foreach gcId,gcIDList do begin
    if gc.subgroupLenType[partType,gcId] gt 0 then begin
      ; list of indices of gas particles in this subgroup
      loc_inds = lindgen(gc.subgroupLenType[partType,gcId]) + gc.subgroupOffsetType[partType,gcId]

      ; select only smooth
      w = where(accMode[loc_inds] eq 1,count)

      if count gt 0 then begin
        loc_inds = loc_inds[w]
        loc_maxt = maxTemp[loc_inds]
      
        ; count fraction Tmax below each constant temperature threshold
        w = where(loc_maxt le TcutVal,count_below)
        coldCount[gcId]  = count_below
        totalCount[gcId] = n_elements(loc_inds)
      endif
      
      print,gc.IDs[loc_inds]
      print,loc_maxt
      print,count_below
      stop
    endif
  endforeach
  
  ; compute cold fraction by halo
  w = where(totalCount ge minNum)
  coldFrac = coldCount[w] / totalCount[w]
  gcMasses = gcMasses[w]
  
  ; calculate median line
  massStep = ceil(100.0/sqrt(n_elements(coldFrac)))/10.0
  ;massStep = 0.2
  massBins = (xrange[1]-xrange[0])/massStep
  massXPts = findgen(massBins)/massBins * (xrange[1]-xrange[0]) + xrange[0] + massStep/2.0
  
  medCold = fltarr(massBins) - 1
  stdCold = fltarr(massBins) - 1
  medCold[0:floor(n_elements(medCold)/2.0)] = 1.0 ;set default value high for first half
  
  for i=0,massBins-1 do begin
    w = where(gcMasses ge xrange[0]+i*massStep and gcMasses lt xrange[0]+(i+1)*massStep,count)
    if (count gt 0) then begin
      medCold[i] = median(coldFrac[w])
      stdCold[i] = stddev(coldFrac[w])
    endif
  endfor  
  
  ; plot
  start_PS, sP.plotPath + 'coldFrac.old.noGC_'+str(sP.res)+'.eos'+str(eosok)+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
      ytitle="Cold Fraction",xtitle=textoidl("log(M_{halo})")+"",$
      title=str(sP.res)+textoidl("^3")+" "+sP.run+" (z="+string(sP.redshift,format='(f3.1)')+") OLD"
    cgPlot,xrange,[0.5,0.5],line=0,color=fsc_color('light gray'),/overplot
    
    ; plot raw points
    psym = 4
    symsize = 0.5
    
    cgPlot,gcMasses,coldFrac,psym=psym,symsize=symsize,thick=1.0,/overplot
    
    ; plot median lines
    cgPlot,massXpts,medCold,color=getColor(3),line=0,/overplot
    
  end_PS
  stop
end

pro compareMaxTemps

  sP = simParams(res=128,run='gadget',redshift=2.0)
  
  gcId = 0
  
  ; load galaxycat and maxtemps for all
  galcat  = galaxyCat(sP=sP)
  maxTemp = maxTempOld(sP=sP, eosok=1)
  
  ; get IDs and maxTemps for all particles in gcId
  inds_gal  = lindgen(galcat.galaxyLen[gcId]) + galcat.galaxyOff[gcId]
  inds_gmem = lindgen(galcat.groupmemLen[gcId]) + galcat.groupmemOff[gcId]

  maxt_gal  = maxTemp.maxTemps_gal[inds_gal]
  maxt_gmem = maxTemp.maxTemps_gmem[inds_gmem]
  
  IDs1  = [galcat.galaxyIDs[inds_gal],galcat.groupmemIDs[inds_gmem]]
  maxt1 = [maxt_gal,maxt_gmem] 
  
  ; load maxtemps for all groupcat IDs
  gc = loadGroupCat(sP=sP,/readIDs)
  maxTemp2 = maxTempOldNoGC(sP=sP, eosok=1)
  
  ; get IDs and maxTemps for all particles in gcId
  inds = lindgen(gc.subgroupLenType[0,gcId]) + gc.subgroupOffsetType[0,gcId]
  
  IDs2  = gc.IDs[inds]
  maxt2 = maxTemp2[inds]
  
  ; match
  match,IDs1,IDs2,ind1,ind2,count=count
  print,count,n_elements(IDs1),n_elements(IDs2)
  
  print,min(maxt1),min(maxt2),max(maxt1),max(maxt2)
  
  ; rearrange them both to be in sorted ID order
  IDs1 = IDs1[ind1]
  IDs2 = IDs2[ind2]
  maxt1 = maxt1[ind1]
  maxt2 = maxt2[ind2]
  
  stop
  ; load at some snapshot and do a ~maxTemps procedure
  snapStart = 180
  snapEnd   = 188
  
  ; maxTemps for the two different kinds
  maxTemp = fltarr(n_elements(inds))
  r = { gal : fltarr(n_elements(inds_gal)), gmem : fltarr(n_elements(inds_gmem)) }
  
  gal_IDs  = galcat.galaxyIDs[inds_gal]
  gmem_IDs = galcat.groupmemIDs[inds_gmem]

  for m=snapStart,snapEnd,1 do begin
    sP.snap = m
  
    ids   = loadSnapshotSubset(sP=sP,partType='gas',field='ids')
    u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
    nelec = loadSnapshotSubset(sP=sP,partType='gas',field='ne')
   
    ; do galaxycat
    match,gal_IDs,ids,gal_ind,ids1_gal_ind,count=countgal,/sort
    ids1_gal_ind = ids1_gal_ind[sort(gal_ind)]
    
    match,gmem_IDs,ids,gmem_ind,ids1_gmem_ind,count=countgmem,/sort
    ids1_gmem_ind = ids1_gmem_ind[sort(gmem_ind)]
    
    temps_gal  = convertUtoTemp(u[ids1_gal_ind],nelec[ids1_gal_ind],/log)
    temps_gmem = convertUtoTemp(u[ids1_gmem_ind],nelec[ids1_gmem_ind],/log)
    
    wgal = where(temps_gal gt r.gal,count1gal)
    if (count1gal gt 0) then r.gal[wgal] = temps_gal[wgal]
    
    wgmem = where(temps_gmem gt r.gmem,count1gmem)
    if (count1gmem gt 0) then r.gmem[wgmem] = temps_gmem[wgmem]
    
    ; do non-galaxycat
    match,IDs2,ids,ids2_loc_ind,ids2_ind,count=count2,/sort
    ids2_ind = ids2_ind[sort(ids2_loc_ind)]
  
    temps2 = convertUtoTemp(u[ids2_ind],nelec[ids2_ind],/log)
    
    w2 = where(temps2 gt maxTemp,count2)
    if (count2 gt 0) then maxTemp[w2] = temps2[w2]
    
    print,m,min(maxTemp),min([r.gal,r.gmem]),max(maxTemp),max([r.gal,r.gmem])
  endfor
  
  stop
  
end
