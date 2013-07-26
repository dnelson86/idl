; checkOldNew.pro
; dnelson jul.2013

; check the old results (pre r99, when galaxyCat was still split into components and no inter existed)
; with new results (after galaxyCat.ids unified, but prior to eliminating the indice subsetting, where
; these results now reside in /data.files.newer/)

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