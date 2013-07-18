; plotMaxTemps.pro
; gas accretion project - plots related to maximum past temperature of the gas
; dnelson jun.2013

pro checkOldNew

  ; load accTimes
  restore,'/n/home07/dnelson/data3/sims.gadget/128_20Mpc/data.files/accTimes.SPH.G128.189-51.sav'
  at_new = r
  restore,'/n/home07/dnelson/data3/sims.gadget/128_20Mpc/data.files.old/accModeTimesOld/accTimesAda.SPH.G128.189-51.sav'
  at_old = r

  ; load accModes
  restore,'/n/home07/dnelson/data3/sims.gadget/128_20Mpc/data.files/accMode.SPH.G128.189-51.sav'
  accMode_new = accMode
  accMode = !NULL
  restore,'/n/home07/dnelson/data3/sims.gadget/128_20Mpc/data.files.old/accModeTimesOld/accMode.SPH.G128.189-51.sav'
  accMode_old = r
  
  ; load mergerTreeSubset
  restore,'/n/home07/dnelson/data3/sims.gadget/128_20Mpc/data.files/mTreeAdaSub.G128.189-51.sK3.sav'
  mt_new = r
  restore,'/n/home07/dnelson/data3/sims.gadget/128_20Mpc/data.files.old/mTreeAdaSub.G128.189-51.sK3.sav'
  mt_old = r
  
  ; lost maxtemps
  restore,'/n/home07/dnelson/data3/sims.gadget/128_20Mpc/data.files/maxVals.SPH.G128.0-188.sav'
  restore,'/n/home07/dnelson/data3/sims.gadget/128_20Mpc/data.files.old/maxtemp.SPH.G128.0-188.sav'
  
  ; load galaxyCatalogs
  restore,'/n/home07/dnelson/data3/sims.gadget/128_20Mpc/data.files.old/galcat.G128.189.sav'
  restore,'/n/home07/dnelson/data3/sims.gadget/128_20Mpc/data.files.old/groupmemcat.G128.189.sav'
  restore,'/n/home07/dnelson/data3/sims.gadget/128_20Mpc/data.files.old/starcat.G128.189.sav'
  galcatOld = { galaxyIDs:galaxyids, groupmemIDs:groupmemids, stellarIDs:stellarIDs ,$
              galaxyOff:galaxyoff, groupmemOff:groupmemoff, stellarOff:stellaroff ,$
              galaxyLen:galaxylen, groupmemLen:groupmemlen, stellarLen:stellarlen }
  
  sP = simParams(res=128,run='gadget',redshift=2.0)
  gcInds = galcatINDList(sP=sP,/gcSep)
  galcatNew = galaxyCat(sP=sP)
  
  colors = ['black','red','blue','green','orange'] ; gal,gmem,stars
  xrange = [4.0,7.0]
  peak   = 1.0
  
  ; maxtemp
  start_PS,'/n/home07/dnelson/old_maxtemp.eps'
    plothist,r.maxtemps_gal,/auto,peak=peak,color=cgColor(colors[0]),xrange=xrange
    plothist,r.maxtemps_gmem,/auto,peak=peak,color=cgColor(colors[1]),/overplot
    plothist,r.maxtemps_stars,/auto,peak=peak,color=cgColor(colors[2]),/overplot
    legend,['gal','gmem','stars'],textcolor=colors[0:2],/top,/right
  end_PS

  start_PS,'/n/home07/dnelson/new_maxtemp.eps'
    plothist,rtr.maxtemps[gcInds.gal],/auto,peak=peak,color=cgColor(colors[0]),xrange=xrange
    plothist,rtr.maxtemps[gcInds.gmem],/auto,peak=peak,color=cgColor(colors[1]),/overplot
    plothist,rtr.maxtemps[gcInds.stars],/auto,peak=peak,color=cgColor(colors[2]),/overplot
    plothist,rtr.maxtemps[gcInds.inter],/auto,peak=peak,color=cgColor(colors[3]),/overplot
    legend,['gal','gmem','stars','inter'],textcolor=colors[0:3],/top,/right
  end_PS
  
  ; check if galaxyCatalogs match?
  ids_gal_old = galaxyids
  ids_gal_new = galcatNew.ids[gcInds.gal]
  
  ids_gal_old = ids_gal_old[sort(ids_gal_old)]
  ids_gal_new = ids_gal_new[sort(ids_gal_new)]
  
  if ~array_equal(ids_gal_new,ids_gal_old) then message,'Fail'
  
  ids_gmem_old = groupmemids
  ids_gmem_new = galcatNew.ids[gcInds.gmem]
  
  ids_gmem_old = ids_gmem_old[sort(ids_gmem_old)]
  ids_gmem_new = ids_gmem_new[sort(ids_gmem_new)]
  
  if ~array_equal(ids_gmem_new,ids_gmem_old) then message,'Fail2'
  
  ; find mtS intersection: GAL
  mts_type = galcatNew.type[mt_new.galcatSub]
  
  mts_galcatSub_gal  = where(mts_type eq 1)
  mts_galcat_gal     = mt_new.galcatSub[mts_galcatSub_gal]
  
  new_gal_ids = galcatNew.ids[mts_galcat_gal]
  old_gal_ids = galcatOld.galaxyids[mt_old.galcatsub.gal]
  
  calcMatch,new_gal_ids,old_gal_ids,ind_new,ind_old,count=countMatch
  if countMatch ne n_elements(new_gal_ids) then message,'Fail'
  
  mts_galcat_gal = mts_galcat_gal[ind_new[sort(ind_old)]] ; rearrange new mts_X_gal to match old order
  mts_galcatSub_gal = mts_galcatSub_gal[ind_new[sort(ind_old)]]
  
  if ~array_equal(galcatNew.ids[mts_galcat_gal],old_gal_ids) then message,'Fail3'
  
  ; find mtS intersection: GMEM
  mts_galcatSub_gmem  = where(mts_type eq 2)
  mts_galcat_gmem     = mt_new.galcatSub[mts_galcatSub_gmem]
  
  new_gmem_ids = galcatNew.ids[mts_galcat_gmem]
  old_gmem_ids = galcatOld.groupmemids[mt_old.galcatsub.gmem]
  
  calcMatch,new_gmem_ids,old_gmem_ids,ind_new,ind_old,count=countMatch
  if countMatch ne n_elements(new_gmem_ids) then message,'Fail'
  
  mts_galcat_gmem = mts_galcat_gmem[ind_new[sort(ind_old)]] ; rearrange new mts_X_gal to match old order
  mts_galcatSub_gmem = mts_galcatSub_gmem[ind_new[sort(ind_old)]]
  
  if ~array_equal(galcatNew.ids[mts_galcat_gmem],old_gmem_ids) then message,'Fail4'
  
  ; accTime
  for i=0,6 do begin
    at_old_gal  = reform( at_old.accTime_gal[i,*] )
    at_old_gal2 = reform( at_old.accTimeRT_gal[*] )
    at_new_gal  = reform( at_new.accTime[i, mts_galcatSub_gal] )
    at_new_gal2 = reform( at_new.accTimeRT[mts_galcatSub_gal] )
  
    start_PS,'/n/home07/dnelson/acctime_'+str(i)+'.eps'
      xrange = [0.0,0.5]
      plothist,at_old_gal,/auto,peak=peak,color=cgColor(colors[0]),xrange=xrange
      plothist,at_new_gal,/auto,peak=peak,color=cgColor(colors[1]),/overplot
    
      plothist,at_old_gal2,/auto,peak=peak,color=cgColor(colors[2]),/overplot
      plothist,at_new_gal2,/auto,peak=peak,color=cgColor(colors[3]),/overplot
      
      legend,['old rad','new rad','old RT','new RT'],textcolor=colors[0:3],/top,/right
    end_PS
  endfor
  
  
  ; compare accTime_new inds=0,6,7
  start_PS,'/n/home07/dnelson/acctime_new_067.eps'
    at_new_gal  = reform( at_new.accTime[0, mts_galcatSub_gal] )
    plothist,at_new_gal,/auto,color=cgColor(colors[0]),xrange=xrange
    
    at_new_gal  = reform( at_new.accTime[6, mts_galcatSub_gal] )
    plothist,at_new_gal,/auto,color=cgColor(colors[1]),/overplot
    
    at_new_gal  = reform( at_new.accTime[7, mts_galcatSub_gal] )
    plothist,at_new_gal,/auto,color=cgColor(colors[2]),/overplot
    legend,['ind0','ind6','ind7'],textcolor=colors[0:2],/top,/right
  end_PS
    
  ; accTime_new: all ind7 (earliest rvir crossing) should be earlier than ind0 (latest rvir crossing)
  w = where(at_new.accTime[7,*] gt at_new.accTime[0,*],count)
    
  if count gt 0 then print,'EARLIEST RVIR CROSSING BROKEN'
    
  ; accHaloTvir
  start_PS,'/n/home07/dnelson/acctime_tvir_gmem.eps'
    xrange = [3.0,8.0]
    plothist, [at_old.accHaloTvir_gmem], /auto,peak=peak,color=cgColor(colors[2]),xrange=xrange
    plothist, at_new.accHaloTvir[[mts_galcatSub_gmem]], /auto,peak=peak,color=cgColor(colors[3]),/overplot
  end_PS
  start_PS,'/n/home07/dnelson/acctime_tvir_gal.eps'
    plothist, at_old.accHaloTvir_gal, /auto,peak=peak,color=cgColor(colors[0]),xrange=xrange
    plothist, at_new.accHaloTvir[mts_galcatSub_gal], /auto,peak=peak,color=cgColor(colors[1]),/overplot
  end_PS
  
  ; number of non-found accTimes matches?
  for i=0,6 do begin
    w = where(at_old.accTime_gal[i,*] ne -1,count_old)
    w = where(at_new.accTime[i,mts_galcatSub_gal] ne -1,count_new)
  
    print,'i='+str(i) + ' old = '+str(count_old) + ' new = ' + str(count_new)
  endfor
  
  w = where(at_new.accTime[7,mts_galcatSub_gal] ne -1,count_new)
  print,'i='+str(7) + ' new = '+str(count_new)
  
  w = where(at_old.accTimeRT_gal ne -1,count_old)
  w = where(at_new.accTimeRT[mts_galcatSub_gal] ne -1,count_new)
  
  print,'RT  old = '+str(count_old)+' new = '+str(count_new)
  
  ; mt and mtS seem to agree
  ; accTime(rad) and accTimeRT both also seem to agree
  ; accMode seems to more or less agree, or at least does not seem that it can account for the
  ; log(T)>6 peak that appears in "new gal smooth"
  ; suspect again this is an indexing error when choosing a mode, i.e. when going all the way
  ; to the end of the index subset hierarchy
  
  ; at_new has (non -1 times) 314297 (.new accTime) or 314838 (old accTime, also searching star parents) ?
  
  ; TODO!
  ; seems prudent to wipe all this out, make all indices the same
  
  ; accMode inds for mapping
  sPatIndMode = 0
  
  w_at_new = where(at_new.accTime[sPatIndMode,*] ne -1)
  ats_type = galcatNew.type[mt_new.galcatSub[w_at_new]]
  
  ats_modeSub_gal   = where(ats_type eq 1)
  ats_galcatSub_gal = w_at_new[ats_modeSub_gal]
  ats_galcat_gal    = mt_new.galcatSub[ats_galcatSub_gal]
  
  ; accMode
  start_PS,'/n/home07/dnelson/accmode.eps'
    xrange = [0,6]
    plothist,accMode_new[ats_modeSub_gal],binsize=1,color=cgColor(colors[0]),xrange=xrange
    h = histogram(accMode_new[ats_modeSub_gal])
    print,'newMode ',h
    print,'newMode ',float(h)/total(h)
    
    plothist,accMode_old.accMode_gal,binsize=1,color=cgColor(colors[1]),/overplot
    h = histogram(accMode_old.accMode_gal)
    print,'oldMode ',h
    print,'oldMode ',float(h)/total(h)
    
    legend,['new gal','old gal'],textcolor=colors[0:1],/top,/right
  end_PS
  
  ; select smooth, plot tmax distribution for smooth only
  smModeInd = 1
  
  w_smooth_old = where(accMode_old.accMode_gal eq smModeInd,count_smooth_old)
  w_smooth_new = where(accMode_new[ats_modeSub_gal] eq smModeInd,count_smooth_new)
  
  print,'smooth: old = '+str(count_smooth_old)+' new = '+str(count_smooth_new)
  
  w_at_old = where(at_old.accTime_gal[0,*] ne -1)
  old_tmax_gal_smooth = r.maxtemps_gal[mt_old.galcatSub.gal[w_at_old[w_smooth_old]]]
  new_tmax_gal_smooth_bad = rtr.maxtemps[mt_new.galcatSub[w_at_new[w_smooth_new]]]
  
  start_PS,'/n/home07/dnelson/tmax_smooth.eps'
    xrange = [4.0,8.0]
    plothist, old_tmax_gal_smooth, /auto, peak=peak, color=cgColor(colors[0]), xrange=xrange
    plothist, new_tmax_gal_smooth_bad, /auto, peak=peak, color=cgColor(colors[1]), /overplot
    legend,['old gal smooth','new gal smooth (bad)'],textcolor=colors[0:1],/top,/right
  end_PS
  
  start_PS,'/n/home07/dnelson/tmax_hierarchy_old.eps'
    plothist,r.maxtemps_gal,/auto,color=cgColor(colors[0]),xrange=xrange
    plothist,r.maxtemps_gal[mt_old.galcatSub.gal],/auto,color=cgColor(colors[1]),/overplot
    plothist,r.maxtemps_gal[mt_old.galcatSub.gal[w_at_old]],/auto,color=cgColor(colors[2]),/overplot
    plothist,r.maxtemps_gal[mt_old.galcatSub.gal[w_at_old[w_smooth_old]]],/auto,color=cgColor(colors[3]),/overplot
    plothist,r.maxtemps_gal[mt_old.galcatSub.gal[w_at_old[w_smooth_old]]],/auto,color=cgColor(colors[3]),/overplot
    legend,['all maxtemps','only those in mtS','those with non neg1 accTimes','smooth'],$
      textcolor=colors[0:3],/top,/right
  end_PS
  
  start_PS,'/n/home07/dnelson/tmax_hierarchy_new.eps'
    ww = where(galcatNew.type eq 1)
    plothist,rtr.maxtemps[ww],/auto,color=cgColor(colors[0]),xrange=xrange
    plothist,rtr.maxtemps[mts_galcat_gal],/auto,color=cgColor(colors[1]),/overplot
    
    ; those with non neg1 accTimes
    w_at_new_gal = intersection( w_at_new, mts_galcatSub_gal )
    plothist,rtr.maxtemps[w_at_new_gal],/auto,color=cgColor(colors[2]),/overplot
    
    ; good smooth indexing
    w_smooth_new2 = where(accMode_new eq smModeInd,count_smooth_new2)
    w_at_new_smooth_gal = intersection( w_at_new[w_smooth_new2], mts_galcatSub_gal)
    plothist,rtr.maxtemps[w_at_new_smooth_gal],/auto,color=cgColor(colors[3]),/overplot
    
    ; bad smooth indexing
    plothist,rtr.maxtemps[mt_new.galcatSub[w_at_new[w_smooth_new]]],/auto,line=2,color=cgColor(colors[4]),/overplot
    ; why is this bad? cannot subset w_at_new/accMode with w_smooth_new, since
    ; w_smooth_new must subset ats_modeSub_gal first
    
    ; corrected smooth indexing
    plothist,rtr.maxtemps[mt_new.galcatSub[w_at_new[ats_modeSub_gal[w_smooth_new]]]],/auto,line=0,color=cgColor(colors[4]),/overplot

    ; BUT: why is the good and the corrected different? can't see a problem with either
    
    legend,['all maxtemps','only those in mtS','those with non neg1 accTimes',$
            'gal smooth (good)','gal smooth (solid ok, dotted bad)'],$
      textcolor=colors[0:4],/top,/right
  end_PS
  
  ; do 2 methods agree on galsmooth indices?
  ind1 = w_at_new_smooth_gal
  ind2 = mt_new.galcatSub[w_at_new[ats_modeSub_gal[w_smooth_new]]]
  
  ; check accModeInds
  ;accModeInds_all    = accModeInds(sP=sP,at=at_new,accMode='all')
  accModeInds_smooth = accModeInds(sP=sP,at=at_new,accMode='smooth')
  
  local_accModeInds_smooth_gal = w_at_new[ats_modeSub_gal[w_smooth_new]]
  
  if ~array_equal(accModeInds_smooth.gal, local_accModeInds_smooth_gal) then message,'Fail'
  
  ; check gcSubsetProp
  gsp_maxtemp_all_new = gcSubsetProp(sP=sP,select='pri',/accretionTimeSubset,/maxPastTemp)
  
  if ~array_equal(rtr.maxtemps[mt_new.galcatSub[w_at_new[ats_modeSub_gal]]],$
    gsp_maxtemp_all_new.gal) then message,'Fail'
    
  gsp_maxtemp_smooth_new = gcSubsetProp(sP=sP,select='pri',accMode='smooth',/accretionTimeSubset,/maxPastTemp)
  
  if ~array_equal(rtr.maxtemps[mt_new.galcatSub[w_at_new[ats_modeSub_gal[w_smooth_new]]]],$
    gsp_maxtemp_smooth_new.gal) then message,'Fail'
  
  ; compare the saved, binned histograms
  restore,'/n/home07/dnelson/data3/sims.gadget/128_20Mpc/data.files.old/binnedVals.old/' + $
    'binTemp.SPH.G128.189.mb7.pri.smooth_tw1000.00.sav'
  bintemp_old = r
  restore,'/n/home07/dnelson/data3/sims.gadget/128_20Mpc/data.files/binnedVals/' + $
    'binMaxVals.SPH.G128.189.mb7.smooth_tw1000.00.sav'
  binmaxvals_new = r
    
  ; Q: why is nonneg-1 at for atIndMode=-1 (7) essentially every sph particle?
  ; earliest rvir crossing for sph is definitely broken (see acctime_new_067.eps)
  
  for i=0,n_elements(bintemp_old.massbins)-2 do begin
    old_yy = bintemp_old.hGalTmax[i,*]
    new_yy = binmaxvals_new.gal.Tmax[i,*]
    
    old_xx = bintemp_old.binLocTemp
    new_xx = binmaxvals_new.params.binLoc.tmax
    
    start_PS,'bintemp_'+str(i)+'.eps'
      cgPlot,old_xx,old_yy,psym=-4,xrange=xrange,color=cgColor(colors[0])
      cgPlot,new_xx,new_yy,psym=-4,/overplot,color=cgColor(colors[1])
      legend,['old','new'],textcolor=colors[0:1],/top,/right
    end_PS
  endfor
    
  stop
end

; binValMaxHistos()

function binValMaxHistos, sP=sP, accMode=accMode, timeWindow=TW

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  if sP.gfmWinds eq 1 and (accMode eq 'smooth' or accMode eq 'stripped' or accMode eq 'clumpy') then begin
    print,'Switching [' + accMode + '] to [' + accMode + '_rec] for GFM run.'
    accMode = accMode + '_rec'
  endif
  
  ; config (1D)
  sgSelect   = 'pri' ; only option for atS
  massBins   = [9.0,9.5,10.0,10.5,11.0,11.5,12.0] ; log(M)
  mmMass     = minmax(massBins)
  nMassBins  = n_elements(massBins)-1
  
  ; config (2D)
  binSizeMass = 0.10 / (sP.res/128)
  mmMass2D    = [9.0,12.0]-[0.0,0.0001] ; log(mhalo [msun])
  nMassBins2D = ceil( (mmMass2D[1]-mmMass2D[0]) / binSizeMass )
  
  ; value binning config
  nValBins = ( mmMass[1] - mmMass[0] ) / binSizeMass
  
  mm = { TmaxTviracc : [-2.2,1.2]  - [0.0,0.0001] ,$ ; log(tmax/tvir)
         Tmax        : [3.8,7.2]   - [0.0,0.0001] ,$ ; log(tmax [K])
         EntMax      : [4.0,10.0]  - [0.0,0.0001] ,$ ; log(cgs)
         DensMax     : [-10.0,0.0] - [0.0,0.0001] ,$ ; log(code)
         MachMax     : [0.0,100.0] - [0.0,0.0001]  }; unitless
  
  for i=0,n_tags(mm)-1 do $
    binSize = mod_struct( binSize, (tag_names(mm))[i], ceil( (mm.(i)[1] - mm.(i)[0])/nValBins*100 )/100.0 ) ; 0.01 round
  binSize.MachMax = ceil(binSize.MachMax) ; 1.0 round
  
  for i=0,n_tags(mm)-1 do begin
    nBins  = mod_struct( nBins, (tag_names(mm))[i], ceil( (mm.(i)[1] - mm.(i)[0])/binSize.(i) ) )
    binLoc = mod_struct( binLoc, (tag_names(mm))[i], fltarr( nBins.(i) ) )
  endfor
  
  ; current time
  h = loadSnapshotHeader(sP=sP)
  curtime = 1/h.time - 1 ; redshift
  curtime = redshiftToAgeFlat(curtime)*1e9 ; yr
  
  snapTimes = redshiftToAgeFlat(snapNumToRedshift(/all,sP=sP))*1e9 ; yr
  
  ; time window to consider accretion over
  if ~keyword_set(TW) then message,'time window required (in Myr)'
  
  if str(TW) eq 'all' then begin
    timeWindow = curtime - redshiftToAgeFlat(6.0)*1e9 ; go back to z=6 (in yr)
  endif else begin
    timeWindow = TW * 1e6 ; convert input Myr to yr
  endelse

  ; check if save exists
  saveFilename = sP.derivPath + 'binnedVals/binMaxVals.' + sP.saveTag + '.' + sP.savPrefix + str(sP.res) + '.' + $
    str(sP.snap) + '.mb' + str(n_elements(massBins)) + '.' + accMode + '_tw' + str(TW) + '.sav'
  
  ; results exist, return
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif  

  ; debug:
  ;gc = galaxyCat(sP=sP)
  ;mt = mergerTreeSubset(sP=sP)
  ;gcInds = galcatINDList(sP=sP,/gcSep)
  
  ; load accretion times (need wAm for accMode selection) for timewindow
  at = accretionTimes(sP=sP)
  wAm = accModeInds(at=at,accMode=accMode,sP=sP)

  ; do the timewindow restriction immediately
  ; ---

    ; GALAXY (1): accretion defined as (rho,temp) joining time or 0.15rvir crossing time (most recent)
    loc_atime = reform(at.accTimeRT[wAm.gal])
    
    r_crossing_time = reform(at.accTime[sP.radIndGalAcc,wAm.gal])
    w = where(r_crossing_time gt loc_atime,count)
    if count gt 0 then loc_atime[w] = r_crossing_time[w]
    
    loc_atime = 1/loc_atime - 1 ; redshift
    loc_atime = redshiftToAgeFlat(loc_atime)*1e9 ; yr
    
    ; make a count of those falling in the time window
    w_gal = where(curtime - loc_atime le timeWindow,nloc)

    ; GMEM (2)
    loc_atime = reform(at.accTime[sP.radIndHaloAcc,wAm.gmem])
    loc_atime = 1/loc_atime - 1 ; redshift
    loc_atime = redshiftToAgeFlat(loc_atime)*1e9 ; yr
    
    ; make a count of those falling in the time window
    w_gmem = where(curtime - loc_atime le timeWindow,nloc)
    
    ; INTER (3) - look at outflow time within this TW
    loc_atime = reform(at.outTimeRT[wAm.inter])
    
    r_crossing_time = reform(at.outTime[sP.radIndGalAcc,wAm.inter])
    w = where(r_crossing_time gt loc_atime,count)
    if count gt 0 then loc_atime[w] = r_crossing_time[w]
    
    loc_atime = 1/loc_atime - 1 ; redshift
    loc_atime = redshiftToAgeFlat(loc_atime)*1e9 ; yr
    
    ; make a count of those falling in the time window
    w_inter = where(curtime - loc_atime le timeWindow,nloc)
    
    ; STARS (4)
    loc_atime = reform(at.accTimeRT[wAm.stars])
    
    r_crossing_time = reform(at.accTime[sP.radIndGalAcc,wAm.stars])
    w = where(r_crossing_time gt loc_atime,count)
    if count gt 0 then loc_atime[w] = r_crossing_time[w]
    
    ; convert from scale factor to age of the universe
    loc_atime = 1/loc_atime - 1 ; redshift
    loc_atime = redshiftToAgeFlat(loc_atime)*1e9 ; yr
    
    ; make a count of those falling in the time window
    w_stars = where(curtime - loc_atime le timeWindow,nloc)

    w_tw = { gal:w_gal, gmem:w_gmem, inter:w_inter, stars:w_stars }
    
    ; BHs
    if sP.gfmBHs ne 0 then begin
      loc_atime = reform(at.accTimeRT[wAm.bhs])
    
      r_crossing_time = reform(at.accTime[sP.radIndGalAcc,wAm.bhs])
      w = where(r_crossing_time gt loc_atime,count)
      if count gt 0 then loc_atime[w] = r_crossing_time[w]
    
      loc_atime = 1/loc_atime - 1 ; redshift
      loc_atime = redshiftToAgeFlat(loc_atime)*1e9 ; yr
    
      ; make a count of those falling in the time window
      w_bhs = where(curtime - loc_atime le timeWindow,nloc)
      w_tw = { gal:w_gal, gmem:w_gmem, inter:w_inter, stars:w_stars, bhs:w_bhs }
    endif

  w   = !NULL
  at  = !NULL
  ;wAm = !NULL
  loc_atime       = !NULL
  r_crossing_time = !NULL
  
  ; binTemp, Ent, Dens, MachNum: load temps and do TW subsets
  gcSP = gcSubsetProp(sP=sP,select=sgSelect,/accTvir,/accretionTimeSubset,accMode=accMode)
  for i=0,n_tags(gcSP)-1 do accTvir = mod_struct( accTvir, (tag_names(gcSP))[i], gcSP.(i)[w_tw.(i)] )
    
  gcSP = gcSubsetProp(sP=sP,select=sgSelect,/maxPastTemp,/accretionTimeSubset,accMode=accMode)  
  for i=0,n_tags(gcSP)-1 do maxTemp = mod_struct( maxTemp, (tag_names(gcSP))[i], gcSP.(i)[w_tw.(i)] )

  ; DEBUG:
  ;start_PS,'debug.eps'
  ;  plothist,maxTemp.gal,/auto
  ;end_PS
  ;stop
  
  gcSP = gcSubsetProp(sP=sP,select=sgSelect,/maxPastEnt,/accretionTimeSubset,accMode=accMode)  
  for i=0,n_tags(gcSP)-1 do maxEnt = mod_struct( maxEnt, (tag_names(gcSP))[i], gcSP.(i)[w_tw.(i)] )
  
  gcSP = gcSubsetProp(sP=sP,select=sgSelect,/maxPastDens,/accretionTimeSubset,accMode=accMode)  
  for i=0,n_tags(gcSP)-1 do maxDens = mod_struct( maxDens, (tag_names(gcSP))[i], gcSP.(i)[w_tw.(i)] )
  
  ; no mach num for SPH
  if sP.trMCPerCell ne 0 then begin
    gcSP = gcSubsetProp(sP=sP,select=sgSelect,/maxPastMach,/accretionTimeSubset,accMode=accMode)  
    for i=0,n_tags(gcSP)-1 do maxMach = mod_struct( maxMach, (tag_names(gcSP))[i], gcSP.(i)[w_tw.(i)] )
  endif else begin
    for i=0,n_tags(gcSP)-1 do maxMach = mod_struct( maxMach, (tag_names(gcSP))[i], fltarr(n_elements(w_tw.(i)))-1 )
  endelse

  ; load parent halo masses so we can make halo massbins
  gcSP = gcSubsetProp(sP=sP,select=sgSelect,/parMass,/accretionTimeSubset,accMode=accMode)
  for i=0,n_tags(gcSP)-1 do parentMass = mod_struct( parentMass, (tag_names(gcSP))[i], gcSP.(i)[w_tw.(i)] )
  
  ; check: if no particles of a type are in w_tw, then replace their entries by a [0,0] array
  ; to allow histogram() to work, since they are now single numbers indexed by -1
  for i=0,n_tags(gcSP)-1 do begin
    if n_elements(w_tw.(i)) ne 1 then continue
    if w_tw.(i) ne -1 then message,'Strange, just one single valid entry?'
    
    accTvir = mod_struct( accTvir, (tag_names(gcSP))[i], [100,100.1] )
    maxTemp = mod_struct( maxTemp, (tag_names(gcSP))[i], [0,0.1] )
    maxEnt  = mod_struct( maxEnt, (tag_names(gcSP))[i], [0,0.1] )
    maxDens = mod_struct( maxDens, (tag_names(gcSP))[i], [-1,-1.1] )
    maxMach = mod_struct( maxMach, (tag_names(gcSP))[i], [-1,-1.1] )
    parentMass = mod_struct( parentMass, (tag_names(gcSP))[i], [0,0.1] )
  endfor

  gcSP = !NULL
  w_tw = !NULL

  ; uniform weighting by mass for 2D histograms
  if sP.trMCPerCell gt 0 then massWt = sP.trMassConst * units.UnitMass_in_Msun
  if sP.trMCPerCell eq 0 then massWt = sP.targetGasMass * units.UnitMass_in_Msun
  if sP.trMCPerCell eq -1 then message,'error'  
  
  ; make return structures
  typeLabels = ['gal','gmem','inter','stars','bhs','allgal','total']
  
  ; which galaxyCat types contribute to the "allGal" and "total"?
  allgalInds = [0,3,4] ; gal,stars,bhs
  totalInds  = [0,1,3,4] ; gal,gmem,stars,bhs (not inter)
  
  allgalInd = ( where( typeLabels eq 'allgal' ) )[0]
  totalInd  = ( where( typeLabels eq 'total' ) )[0]
  
  rr = { binLoc      : binLoc      ,$
         binSize     : binSize     ,$
         mm          : mm          ,$
         massBins    : massBins    ,$ ; 1D
         binSizeMass : binSizeMass ,$ ; 2D
         mmMass2D    : mmMass2D     } ; 2D
       
  ; 1D, vs halo mass and global, and 2D
  for i=0,n_tags(mm)-1 do $
    template = mod_struct( template, (tag_names(mm))[i], fltarr(nMassBins, nBins.(i)) )
  for i=0,n_tags(mm)-1 do $
    template = mod_struct( template, 'Global' + (tag_names(mm))[i], fltarr( nBins.(i)) )
  for i=0,n_tags(mm)-1 do $
    template = mod_struct( template, 'h2_' + (tag_names(mm))[i], fltarr(nMassBins2D, nBins.(i)) )
    
  for i=0,n_elements(typeLabels)-1 do $
    r = mod_struct( r, typeLabels[i], template)
        
  ; loop over halo mass bins
  for j=0,n_elements(massBins)-2 do begin
  
    ; for each type (gal,gmem,stars,inter,bhs) histogram within mass bin
    foreach tInd,totalInds do begin
      if sP.gfmBHs eq 0 and typeLabels[tInd] eq 'bhs' then continue
      if (tag_names(parentMass))[tInd] ne strupcase(typeLabels[tInd]) then message,'Careful'
      
      w_type = where(parentMass.(tInd) gt massBins[j] and parentMass.(tInd) le massBins[j+1],count)
    
      if count eq 0 then continue
    
      ; binTemp: log(tmax/tviracc)
      vals = [10.0^maxTemp.(tInd)[w_type]/10.0^accTvir.(tInd)[w_type]]
      r.(tInd).TmaxTviracc[j,*] = histogram(alog10(vals),$
        binsize=binsize.TmaxTviracc,min=mm.TmaxTviracc[0],max=mm.TmaxTviracc[1])
    
      ; binTemp: log(tmax)
      r.(tInd).Tmax[j,*] = histogram(maxTemp.(tInd)[w_type],binsize=binsize.Tmax,min=mm.Tmax[0],max=mm.Tmax[1])
      
      if tInd eq 0 then begin
      start_PS,'debug_'+str(j)+'.eps'
        plothist,maxTemp.(tInd)[w_type],/auto
      end_PS
      endif
      
      ; binEnt, binDens, binMach
      r.(tInd).EntMax[j,*]  = histogram(maxEnt.(tInd)[w_type],$
        binsize=binsize.EntMax,min=mm.EntMax[0],max=mm.EntMax[1])
      r.(tInd).DensMax[j,*] = histogram(alog10(maxDens.(tInd)[w_type]),$
        binsize=binsize.DensMax,min=mm.DensMax[0],max=mm.DensMax[1])
      r.(tInd).MachMax[j,*] = histogram(maxMach.(tInd)[w_type],$
        binsize=binsize.MachMax,min=mm.MachMax[0],max=mm.MachMax[1])
    endforeach
    
  endfor
  
  ; GLOBAL quantities: for each type (gal,gmem,stars,inter,bhs) histogram across all halo masses
  foreach tInd,totalInds do begin
    if sP.gfmBHs eq 0 and typeLabels[tInd] eq 'bhs' then continue
    
    ;binTemp: global tmax/tviracc and global tmax
    vals = [10.0^maxTemp.(tInd)/10.0^accTvir.(tInd)]
    r.(tInd).GlobalTmaxTviracc[*] = histogram(alog10(vals),$
      binsize=binsize.TmaxTviracc,min=mm.TmaxTviracc[0],max=mm.TmaxTviracc[1],loc=loc)
    rr.binLoc.TmaxTviracc = loc + binSize.TmaxTviracc*0.5 ; save ratio midbins
    
    r.(tInd).GlobalTmax[*] = histogram(maxTemp.(tInd),$
      binsize=binsize.Tmax,min=mm.Tmax[0],max=mm.Tmax[1],loc=loc)
    rr.binLoc.Tmax = loc + binSize.Tmax*0.5 ; save temp [K] midbins
    
    ; binEnt, binDens, binMach
    r.(tInd).GlobalEntMax[*] = histogram(maxEnt.(tInd),$
      binsize=binsize.EntMax,min=mm.EntMax[0],max=mm.EntMax[1],loc=loc)
    rr.binLoc.EntMax = loc + binsize.EntMax*0.5
    
    r.(tInd).GlobalDensMax[*] = histogram(alog10(maxDens.(tInd)),$
      binsize=binsize.DensMax,min=mm.DensMax[0],max=mm.DensMax[1],loc=loc)
    rr.binLoc.DensMax = loc + binsize.DensMax*0.5
    
    r.(tInd).GlobalMachMax[*] = histogram(maxMach.(tInd),$
      binsize=binsize.MachMax,min=mm.MachMax[0],max=mm.MachMax[1],loc=loc)
    rr.binLoc.MachMax = loc + binsize.MachMax*0.5
    
    ; 2D histograms: temp
    vals = alog10( 10.0^maxTemp.(tInd) / 10.0^accTvir.(tInd) )
    r.(tInd).h2_TmaxTviracc = hist_nd( transpose([[parentMass.(tInd)],[vals]]), [binSizeMass,binSize.TmaxTviracc], $
                                min=[mmMass2D[0],mm.TmaxTviracc[0]], max=[mmMass2D[1],mm.TmaxTviracc[1]] ) * massWt
                                   
    r.(tInd).h2_Tmax = hist_nd( transpose([[parentMass.(tInd)],[maxTemp.(tInd)]]), [binSizeMass,binSize.Tmax], $
                                min=[mmMass2D[0],mm.Tmax[0]], max=[mmMass2D[1],mm.Tmax[1]] ) * massWt
 
    ; 2D histograms: ent, dens, mach
    r.(tInd).h2_EntMax = $
      hist_nd( transpose([[parentMass.(tInd)],[maxEnt.(tInd)]]), [binSizeMass,binSize.EntMax], $
               min=[mmMass2D[0],mm.EntMax[0]], max=[mmMass2D[1],mm.EntMax[1]] ) * massWt
    r.(tInd).h2_DensMax = $
      hist_nd( transpose([[parentMass.(tInd)],[alog10(maxDens.(tInd))]]), [binSizeMass,binSize.DensMax], $
               min=[mmMass2D[0],mm.DensMax[0]], max=[mmMass2D[1],mm.DensMax[1]] ) * massWt
    r.(tInd).h2_MachMax = $
      hist_nd( transpose([[parentMass.(tInd)],[maxMach.(tInd)]]), [binSizeMass,binSize.MachMax], $
               min=[mmMass2D[0],mm.MachMax[0]], max=[mmMass2D[1],mm.MachMax[1]] ) * massWt
  endforeach
  
  ; do allGal and total
  for k=0,n_tags(r.(0))-1 do begin
    foreach q,allgalInds do r.(allgalInd).(k)[*] += r.(q).(k)[*]
    foreach q,totalInds  do r.(totalInd).(k)[*] += r.(q).(k)[*]
  endfor
  
  ; add general configuration parameters to save struct
  r = mod_struct( r, 'params', rr )
  
  ; DEBUG:
  ;start_PS,'debug2.eps'
  ;  cgPlot,rr.binLoc.Tmax,r.gal.GlobalTmax,psym=4
  ;end_PS
  ;stop
  
  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  return, r

end

; plotValMaxHistos(); plot (1) the previous max temp normalized by tviracc for arepo vs. gadget, gal vs. 
;                     gmem, (2) same but unnormalized by tviracc, (3) global not binned by halo mass but
;                     normalized by each parent tviracc, (4) same global without normalization

pro plotValMaxHistos

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  runs       = ['gadget','tracer','feedback'] ;['feedback','feedback_noZ','feedback_noFB']
  redshift   = 2.0
  res        = 256
  timeWindow = 1000.0 ; Myr
  accMode    = 'smooth' ;'all','smooth','clumpy','stripped','recycled'
  tagNames   = ['Tmax'] ;['TmaxTviracc','Tmax','EntMax','DensMax','MachMax'] ; plot each quantity
  
  ; plot config
  lines   = [0,1] ; gal,gmem
  cInd    = 1 ; color index for simName labels
  
  ; 2D plot config
  exp     = 0.5   ; gamma exponent for non-linear color scaling
  ndivs   = 5     ; number of divisions on colorbar   
  Tc_val  = 5.5   ; log(K) for constant temp line
  lines2D = [1,2] ; Tc,Tvir line styles
  colors  = ['black','black'] ; Tc,Tvir line colors
  
  if isnumeric(timeWindow) then twStr = str(fix(timeWindow))
  if ~isnumeric(timeWindow) then twStr = "-"+timeWindow
     
  ; load
  foreach run,runs,i do begin
    sP  = mod_struct(sP, 'sP'+str(i), simParams(res=res,run=run,redshift=redshift))
    bV  = mod_struct(bV, 'bV'+str(i), binValMaxHistos(sP=sP.(i),accMode=accMode,timeWindow=timeWindow))
  endforeach
    
  ; strings
  plotStr   = ''
  simNames  = []
  simColors = []
  
  foreach run,runs,i do begin
    plotStr   = plotStr + sP.(i).plotPrefix + '.'
    simNames  = [simNames, sP.(i).simName]
    simColors = [simColors, sP.(i).colors[cInd]]
  endforeach

  plotStr += str(res) + '_' + str(sP.(0).snap) + '_tw' + twStr + '_am-' + accMode
    
  ; plot config
  pos     = plot_pos(rows=2, cols=3)
  yrange  = [6e-4,2.0]
  pParams = { TmaxTviracc : {xrange:[-2.2,1.4], xtickv:[-2,-1,0,1],   ylabel : "log ( T_{max} / T_{vir,acc} )"} ,$
              Tmax        : {xrange:[3.0,8.0],  xtickv:[4,5,6,7],     ylabel : "log ( T_{max} )"}               ,$
              EntMax      : {xrange:[5.5,9.5],  xtickv:[6,7,8,9],     ylabel : "log ( S_{max} ) [K cm^{2 }]"}   ,$
              DensMax     : {xrange:[-10,0],    xtickv:[-8,-6,-4,-2], ylabel : "log ( \rho_{max} )"}            ,$
              MachMax     : {xrange:[0,100],    xtickv:[10,20,50,80], ylabel : "M_{max}"}                        }     
              
  foreach tagName,tagNames do begin

    bVind   = ( where(tag_names(bV.(0).(0)) eq strupcase(tagName)) )[0]
    bVind2d = ( where(tag_names(bV.(0).(0)) eq 'H2_'+strupcase(tagName)) )[0]
    pPind   = ( where(tag_names(pParams) eq strupcase(tagName)) )[0]
    if bVind eq -1 or pPind eq -1 or bvInd2d eq -1 then message,'Error'
    print,tagName,bVind,pPind,bVind2d
    
    ; plot (1) - 3x2 mass bins separated out and each panel with all runs, gal vs. gmem (VIR)
    start_PS, sP.(0).plotPath + tagName + '_3x2_allGal_gmem.'+plotStr+'.eps', /big
      !p.thick += 1
      
      xrange = pParams.(pPind).xrange
      xtickv = pParams.(pPind).xtickv
      ylabel = pParams.(pPind).ylabel
      
      for j=0,n_elements(bV.(0).params.massBins)-2 do begin
        
        if j eq 0 or j eq 3 then ytickname = '' else ytickname = replicate(' ',10)
        if j ge 3 then xtickname = '' else xtickname = replicate(' ',10)
        if j gt 0 then noerase = 1 else noerase = 0
        
        cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,pos=pos[j],$
          ytitle="",xtitle="",xticks=n_,xtickv=xtickv,xtickname=xtickname,ytickname=ytickname,noerase=noerase       
        
        cgPlot,[0,0],[8e-4,0.25],line=2,color=cgColor('black'),thick=!p.thick-0.0,/overplot
        
        for i=0,n_tags(sP)-1 do begin
          ; gal
          hist = bV.(i).allGal.(bVind)[j,*]
          cgPlot,bV.(i).params.binLoc.(bvInd),float(hist)/total(hist),line=lines[0],color=sP.(i).colors[cInd],/overplot

          ; gmem
          hist = bV.(i).gmem.(bVind)[j,*]
          cgPlot,bV.(i).params.binLoc.(bvInd),float(hist)/total(hist),line=lines[1],color=sP.(i).colors[cInd],/overplot
        endfor
        
        ; legends
        massBinStr = string(bV.(0).params.massBins[j],format='(f4.1)') + ' < log(M) < ' + $
                     string(bV.(0).params.massBins[j+1],format='(f4.1)')
  
        cgText,mean(xrange),yrange[1]*0.4,massBinStr,charsize=!p.charsize-0.0,alignment=0.5
            
        if j eq 0 then legend,simNames,textcolors=simColors,box=0,position=[xrange[0],0.5]
        if j eq 3 then legend,['galaxy','halo'],linestyle=[0,1],color=cgColor('dark gray'),$
          textcolors=['dark gray','dark gray'],linesize=0.25,box=0,position=[xrange[0],0.5]
      endfor
      
      ; axis labels
      cgText,0.05,mean([ (pos[0])[3], (pos[3])[1] ]),"Gas Mass Fraction",alignment=0.5,orientation=90.0,/normal
      cgText,mean([ (pos[0])[0], (pos[2])[3] ]),0.05,textoidl(ylabel),alignment=0.5,/normal
      
    end_PS
    
    ; plot (2) - 2D histogram
    if n_elements(runs) eq 2 then $
      pos2D = list( [0.07,0.18,0.44,0.94]  ,$ ; run1
                    [0.47,0.18,0.84,0.94] ,$ ; run2
                    [0.89,0.18,0.93,0.94]   ) ; cbar
    if n_elements(runs) eq 3 then $
      pos2D = list( [0.06,0.18,0.32,0.94] ,$ ; run1
                    [0.39,0.18,0.65,0.94] ,$ ; run2
                    [0.72,0.18,0.98,0.94]  ) ; run3
    
    start_PS, sP.(0).plotPath + tagName + '2D_allgal.'+plotStr+'.eps', xs=13.0, ys=4.0
          
      ; color range (same for all three panels, NOT used)
      h2all = []
      for i=0,n_tags(sP)-1 do h2all = [h2all,bV.(i).allGal.(bVind2d)]
      crange = minmax(h2all[where(h2all ne 0.0)]) ;* [2.0,0.7]
      print,crange
      
      ; loop over each run
      for i=0,n_tags(sP)-1 do begin
      
        xtickv   = [9,10,11,12]
        xrange2D = bV.(i).params.mmMass2D
        yrange2D = bV.(i).params.mm.(bVind)
        xtitle   = textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]")
        ytitle   = textoidl( pParams.(pPind).ylabel )
      
        ; plot 2d histo
        h2mt   = bV.(i).allGal.(bVind2d)
        
        loadColorTable, 'helix', /reverse
        tvim,h2mt^exp,scale=0,clip=-1,xtitle=xtitle,ytitle=ytitle,xrange=xrange2D,yrange=yrange2D,$
           xticks=n_elements(xtickv)-1,xtickv=xtickv,position=pos2D[i],noerase=(i gt 0),range=crange^exp
           
        ; temp lines and legend
        ;tvir_vals = alog10(10.0^Tc_val / codeMassToVirTemp(10.0^xrange2D/units.UnitMass_in_Msun,sP=sP.(i)))
        ;cgPlot,xrange2D,tvir_vals,line=lines2D[0],color=cgColor(colors[1]),/overplot
        
        ; simName
        legend,[sP.(i).simName],textcolor=[sP.(i).colors[cInd]],/bottom,/left
           
      endfor
           
      ; colorbar (NOTE: range is just representative of the last 2d histo)
      if n_elements(runs) eq 2 then begin
        barvals = findgen(ndivs+1)/ndivs*(max(h2mt^exp)-min(h2mt^exp)) + min(h2mt^exp)
        ticknames = textoidl(str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)')))
        ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
      
        loadColorTable, 'helix', /reverse
        cgColorbar,bottom=1,range=minmax(h2mt),position=pos2D[2],/vertical,/right,title="",$
           divisions=ndivs,ticknames=ticknames,ncolors=255
           
        cgText,textoidl("M_{tot} [_{ }log h^{-1} M_{sun }]"),0.98,(0.18+0.94)/2,alignment=0.5,orientation=90,/normal
      endif
                
    end_PS
        
  endforeach ; tagNames
  
  stop
end
