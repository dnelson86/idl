; coldflowsLoad.pro
; cold flows - main
; dnelson sep.2011

@helper
@coldflowsLoad
@coldflowsUtil
@coldflowsVis

; findThermalHistories(): for a given subhalo selection at a given redshift, load all previous
;                         snapshots and record the density+temperature for each gas particle in 
;                         that halo for all prior redshifts
;
; smoothAccretionOnly=1:  only include gas particles from the target subhalo that were not a 
;                         member of any other subfind group at the previous output time

pro findThermalHistories, res=res

  if not keyword_set(res) then begin
    print,'Error: Must specific resolution set.'
    return
  endif
  
  ; config
  gadgetPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'
  workingPath  = '/n/hernquistfs1/dnelson/coldflows/thermhist/'
  
  smoothAccretionOnly = 1
  
  redshiftBins = [6.0,5.0,4.0,3.0,2.0,1.5,1.0,0.5,0.25,0.0]

  ; do complete analysis for all previous time for each redshift bin of "smooth accretion"
  for j=0,n_elements(redshiftBins)-2 do begin
    targetRedshift = redshiftBins[j+1]
    endingRedshift = redshiftBins[j]
    
    ; make halo selection at target redshift
    targetSnap = redshiftToSnapNum(targetRedshift)
    endingSnap = redshiftToSnapNum(endingRedshift)
    
    print,''
    print,'Runing: ',targetRedshift,endingRedshift,targetSnap,endingSnap
    print,''
  
    saveFilename = workingPath + 'thermhist.'+str(res)+'.m='+str(targetSnap)+'.to.'+str(endingSnap)+'.sav'
    
    ; load subhalo group catalogs
    sg    = loadSubhaloGroups(gadgetPath,targetSnap)
    
    ; make list of sgIDs excluding background (id=0) subgroups
    valSGids   = []
    prevGrNr   = -1
    
    for i=0,n_elements(sg.subgroupLen)-1 do begin
      if (sg.subgroupGrnr[i] ne prevGrNr) then begin
        prevGrNr = sg.subgroupGrnr[i]
      endif else begin
        valSGids = [valSGids,i]
      endelse
    endfor
  
    if (keyword_set(smoothAccretionOnly)) then begin
      sgEnd = loadSubhaloGroups(gadgetPath,endingSnap)
      saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.sAO.sav'
      print,'Smooth Accretion Only!'
    endif else begin
      print,'Smooth Accretion DISABLED!'
    endelse
  
    if not (file_test(saveFileName)) then begin  
    
      ; load particle counts from header of first snapshot
      h = loadSnapshotHeader(gadgetPath, snapNum=0)
      
      ; load gas ids from targetSnap to restrict to gas
      gas_ids = loadSnapshotSubset(gadgetPath,snapNum=targetSnap,partType='gas',field='ids',/verbose)
      
      ; for all subhalos
      countSGIDs_All = 0
      sgIDs_All = []
      
      ; loop over each valid subhalo
      foreach sgID, valSGids do begin
      
        ; select subgroup  
        nPart = sg.subGroupLen[sgID]
        sgIDs = sg.subGroupIDs[sg.subGroupOffset[sgID] : sg.subGroupOffset[sgID] + sg.subGroupLen[sgID] - 1]
        
        sgIDs_All = [sgIDs_All, sgIDs]
        countSGIDs_All += n_elements(sgIDs)
      endforeach
      
      ; restrict to only gas
      match, gas_ids, sgIDs_All, gas_ids_ind, sgIDs_ind, count=count_gas
      
      if (count_gas gt 0) then begin
          print,'At z='+str(targetRedshift)+' found ['+str(n_elements(sgIDs_All))+'] particles (all subhalos)'+$
                ' after gas particle cut have ['+str(count_gas)+'] left, lost '+str(n_elements(sgIDs_All)-count_gas)+'.'
                
          ; keep only sgIDs[sgIDs_ind]
          sgIDs_All = sgIDs_All[sgIDs_ind]
      endif
      
      ; enforce smooth accretion only (if requested)
      if (keyword_set(smoothAccretionOnly)) then begin
        match, sgEnd.subgroupIDs, sgIDs_All, snap_ind, sgids_loc_ind, count=count_sm
        
        if (count_sm gt 0) then begin
          ; remove sgIDs[sgids_loc_ind] and modify count
          all = bytarr(n_elements(sgIDs_All))
          if (sgids_loc_ind[0] ne -1L) then all[sgids_loc_ind] = 1B
          w = where(all eq 0B, ncomp)
  
          if (ncomp ne n_elements(sgIDs_All)-count_sm) then begin
            print,'ERROR',ncomp,n_elements(sgIDs_All),count_sm
            return
          endif
          
          print,'At z='+str(targetRedshift)+' found ['+str(n_elements(sgIDs_All))+'] gas particles (all subhalos)'+$
                ' after smooth accretion cut have ['+str(ncomp)+'] left, lost '+str(count_sm)+'.'
          
          sgIDs_All = sgIDs_All[w]
        endif
      endif
      
      countSGIDs_All = n_elements(sgIDs_All)
  
      ; allocate struct
      temp    = fltarr(countSGIDs_All,targetSnap+1)
      density = fltarr(countSGIDs_All,targetSnap+1)
      masses  = fltarr(countSGIDs_All,targetSnap+1)
      
      ; loop from target snapshot through all prior snapshots
      for m=targetSnap,0,-1 do begin
        ; load ids and make particle selection
        ids = loadSnapshotSubset(gadgetPath,snapNum=m,partType='gas',field='ids')
        match, ids, sgIDs_All, ids_ind, sgIDs_ind, count=count
        ids = 0 ;free
        
        ; load u and nelec and calculate temperature
        u     = loadSnapshotSubset(gadgetPath,snapNum=m,partType='gas',field='u')
        u     = u[ids_ind]
        nelec = loadSnapshotSubset(gadgetPath,snapNum=m,partType='gas',field='ne')
        nelec = nelec[ids_ind]
        
        t = convertUtoTemp(u,nelec)
        
        ; load density and masses
        rho   = loadSnapshotSubset(gadgetPath,snapNum=m,partType='gas',field='density')
        rho   = rho[ids_ind]
        mass  = loadSnapshotSubset(gadgetPath,snapNum=m,partType='gas',field='mass')
        mass  = mass[ids_ind]
        
        ; density or temperature selection?
        
        ; save properties
        temp[sgIDs_ind,m]    = t
        density[sgIDs_ind,m] = rho
        masses[sgIDs_ind,m]  = mass
        
      endfor
    
      ; save/restore
      save,temp,density,masses,sgIDs_All,h,filename=saveFileName
    endif else begin
      restore,saveFileName
    endelse 
  
  endfor ; redshiftBins
  
  print,''
  print,'Done.'
end

; plot TmaxRedshift():

pro plotTmaxRedshift, res=res

  units = getUnits()

  if not keyword_set(res) then begin
    print,'Error: Must specific resolution set.'
    return
  endif
  
  ; config
  dataPath  = '/n/hernquistfs1/dnelson/coldflows/thermhist/'
  plotPath     = '/n/home07/dnelson/coldflows/'

  smoothAccretionOnly = 1
  
  redshiftBins = [6.0,5.0,4.0,3.0,2.0,1.5,1.0,0.5,0.25,0.0]  

  ; plot axes
  redshifts = snapNumToRedshift(/all)
  times     = redshiftToAge(redshifts)
  
  xrange = [0.0,max(times)]
  yrange = [0.0,1.05]
 
  ; pdf
  pdfRes = 1000.0
  pdfWidth = 20.0 * max(times) / pdfRes
  pdfPts = findgen(pdfRes)/pdfRes * max(times)
  pdf = fltarr(n_elements(redshiftBins)-1,pdfRes)
    
  ; plot
  start_PS,plotPath+'maxt_redshift_'+str(res)+'.eps';,xs=7,ys=6

  fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,$
           xtitle="Time [Gyr]",ytitle="Normalized PDF",xs=9,/ys,ymargin=[4,3]
  redshift_axis,xrange,yrange,/dotted
  fsc_text,xrange[1]*0.9,yrange[1]*0.9,str(res)+"^3",alignment=0.5
  
  for j=0,n_elements(redshiftBins)-2 do begin
    targetRedshift = redshiftBins[j+1]
    endingRedshift = redshiftBins[j]
    
    ; make halo selection at target redshift
    targetSnap = redshiftToSnapNum(targetRedshift)
    endingSnap = redshiftToSnapNum(endingRedshift)
  
    saveFilename = dataPath + 'thermhist.'+str(res)+'.m='+str(targetSnap)+'.to.'+str(endingSnap)+'.sav'
    
    if (keyword_set(smoothAccretionOnly)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.sAO.sav'
    
    if not (file_test(saveFileName)) then begin
      print,'WARNING: Cannot load savefile: '+saveFilename
      continue
    endif
  
    restore, saveFilename
    print,'Working: '+saveFilename
    
    ; calc max temp times
    maxTempTimes = fltarr(n_elements(sgIDs_All))
    for i=0,n_elements(sgIDs_All)-1 do begin
      w = where(max(temp[i,*]) eq temp[i,*],count)
      if (count ne 1) then begin
        if (count ne 2) then begin
          print,'ERROR'
          stop
        endif
        print,'WARNING: Two snaps with same max temps',w[0],temp[i,w[0]],w[1],temp[i,w[1]]
        w=w[0]
      endif
      maxTempTimes[i] = times[w]
    endfor
    
    ; calculate PDF
    for i=0,pdfRes-2 do begin
      w = where(maxTempTimes ge pdfPts[i]-pdfWidth/2.0 and $
                maxTempTimes lt pdfPts[i+1]+pdfWidth/2.0, count)
      pdf[j,i] = count
    endfor
 
  endfor ;j
  
  ; normalize
  pdf /= max(pdf)  
  
  redshiftNames = ['6','5','4','3','2','1.5','1.0','0.5','0.25','0']
  for j=0,n_elements(redshiftBins)-2 do begin
    ; plot pdf
    fsc_plot,pdfPts,pdf[j,*],/overplot,color=fsc_color(units.colors[j])
    
    ; plot legend entry
    ypos = yrange[1]*0.8 - yrange[1]*0.05*j
    ;fsc_plot,[xrange[1]*0.9,xrange[1]*0.92],[ypos,ypos],line=0,color=fsc_color(units.colors[j]),/overplot
    fsc_text,xrange[1]*0.97,ypos,redshiftNames[j+1]+" < z < "+redshiftNames[j],$
             alignment=1.0,color=fsc_color(units.colors[j]),charsize=!p.charsize-0.5
  endfor
  
  end_PS

end

; plotModeFracVsSubhaloMass()

pro plotModeFracVsSubhaloMass, res=res

  units = getUnits()

  if not keyword_set(res) then begin
    print,'Error: Must specific resolution set.'
    return
  endif
  
  ; config
  gadgetPath  = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'
  dataPath    = '/n/hernquistfs1/dnelson/coldflows/thermhist/'
  plotPath    = '/n/home07/dnelson/coldflows/'

  smoothAccretionOnly = 1
  critLogTemp = 5.5
  
  minNumGasPart = 10 ;smoothly accreted gas per subhalo
  
  redshiftBins = [6.0,5.0,4.0,3.0,2.0,1.5,1.0,0.5,0.25,0.0]  
  binVals      = [2,3,5,8] ;z=3,2,1,0 halo selection

  ; plot
  start_PS,plotPath+'frac_shmass_'+str(res)+'.eps'
  !p.multi = [0,2,2]

  for j=0,n_elements(binVals)-1 do begin
    targetRedshift = redshiftBins[binVals[j]+1]
    endingRedshift = redshiftBins[binVals[j]]
    
    ; make halo selection at target redshift
    targetSnap = redshiftToSnapNum(targetRedshift)
    endingSnap = redshiftToSnapNum(endingRedshift)
  
    saveFilename = dataPath + 'thermhist.'+str(res)+'.m='+str(targetSnap)+'.to.'+str(endingSnap)+'.sav'
    
    if (keyword_set(smoothAccretionOnly)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.sAO.sav'
    
    if not (file_test(saveFileName)) then begin
      print,'WARNING: Cannot load savefile: '+saveFilename
      continue
    endif
    
    restore,saveFilename
    print,'Working: '+saveFilename
    
    ; calc max temp
    maxTemps = fltarr(n_elements(sgIDs_All))
    for i=0,n_elements(sgIDs_All)-1 do begin
      maxTemps[i] = alog10(max(temp[i,*]))
    endfor    
    
    ; load subhalo catalog
    sg    = loadSubhaloGroups(gadgetPath,targetSnap)
    
    ; make list of sgIDs excluding background (id=0) subgroups
    valSGids   = []
    prevGrNr   = -1
    
    for i=0,n_elements(sg.subgroupLen)-1 do begin
      if (sg.subgroupGrnr[i] ne prevGrNr) then begin
        prevGrNr = sg.subgroupGrnr[i]
      endif else begin
        valSGids = [valSGids,i]
      endelse
    endfor
    
    ; subhalo masses (x-axis)
    shmass = alog10( (units.UnitMass_in_g / units.Msun_in_g) * sg.subgroupMass[valSGids] )
    
    ; subhalo cold fractions (y-axis)
    coldfrac = fltarr(n_elements(valSGids))
    
    for i=0,n_elements(valSGids)-1 do begin
      sgIDs = sg.subGroupIDs[sg.subGroupOffset[valSGids[i]] : $
                             sg.subGroupOffset[valSGids[i]] + sg.subGroupLen[valSGids[i]] - 1]

      match,sgIDs,sgIDs_All,sh_ind,sgIDs_ind,count=count
      
      if (count eq 0 or count lt minNumGasPart) then begin
        ;print,'Warning: halo sgid='+str(i)+' had no matches in sgIDs_All'
        coldfrac[i] = -0.1
        continue
      endif
      
      wCold = where(maxTemps[sgIDs_ind] lt critLogTemp, count_cold)
      if (count_cold ne 0) then begin
        coldfrac[i] = float(count_cold) / count
      endif
    endfor

    ; plot config
    xrange = [8.5,12.5]
    yrange = [-0.08,1.08]
    
    plotsym,0 ;circle
    symsize = 0.3
    
    if (j eq 0) then begin
      fsc_plot,shmass,coldfrac,psym=8,ymargin=[1.0,2.0],xmargin=[7.0,0.0],$
               xrange=xrange,yrange=yrange,/xs,/ys,xtickname=replicate(' ',10),symsize=symsize
      fsc_text,xrange[1]*0.96,yrange[1]*0.72,"z=3",alignment=0.5,color=fsc_color('orange')
    endif
    if (j eq 1) then begin
      fsc_plot,shmass,coldfrac,psym=8,ymargin=[1.0,2.0],xmargin=[0.0,7.0],$
               xrange=xrange,yrange=yrange,/xs,/ys,$
               xtickname=replicate(' ',10),ytickname=replicate(' ',10),symsize=symsize
      fsc_text,xrange[1]*0.96,yrange[1]*0.72,"z=2",alignment=0.5,color=fsc_color('orange')
    endif
    if (j eq 2) then begin
      fsc_plot,shmass,coldfrac,psym=8,ymargin=[4.0,-1.0],xmargin=[7.0,0.0],$
               xrange=xrange,yrange=yrange,/xs,/ys,symsize=symsize
      fsc_text,xrange[1]*0.96,yrange[1]*0.72,"z=1",alignment=0.5,color=fsc_color('orange')
    endif
    if (j eq 3) then begin
      fsc_plot,shmass,coldfrac,psym=8,ymargin=[4.0,-1.0],xmargin=[0.0,7.0],$
               xrange=xrange,yrange=yrange,/xs,/ys,ytickname=replicate(' ',10),$
               symsize=symsize
      fsc_text,xrange[1]*0.96,yrange[1]*0.72,"z=0",alignment=0.5,color=fsc_color('orange')
    endif

    ; median line
    ;massStep = round(100.0/sqrt(n_elements(valSGids)))/10.0
    massStep = 0.2
    massBins = (xrange[1]-xrange[0])/massStep
    massXPts = findgen(massBins)/massBins * (xrange[1]-xrange[0]) + xrange[0] + massStep/2.0
    
    medCold = fltarr(massBins)
    medCold[0:floor(n_elements(medCold)/2.0)] = 1.0 ;set default value high for first half
    
    for i=0,massBins-1 do begin
      w = where(shmass ge xrange[0]+i*massStep and shmass lt xrange[0]+(i+1)*massStep and $
                coldfrac ne -0.1,count)
      if (count gt 0) then begin
        medCold[i] = mean(coldfrac[w])
      endif
    endfor
    
    ; plot smoothed median line
    fsc_plot,massXPts,smooth(medCold,3),color=fsc_color('blue'),line=0,/overplot
    fsc_plot,massXPts,smooth(1.0-medCold,3),color=fsc_color('red'),line=0,thick=!p.thick-0.5,/overplot
    
    ; fit for medCold=0.5
    w = where(massXPts ge 10.5 and massXPts le 11.5 and medCold ne 1.0)
    fit = linfit(massXPts[w],medCold[w])
    fitSHMass = (0.5 - fit[0]) / fit[1]

    print,'cold frac 1/2 best fit subhalo mass = '+str(fitSHMass)

  endfor ;j
  
  fsc_text,0.5,0.95,str(res)+"^3",alignment=0.5,/normal
  fsc_text,0.5,0.05,"log ( Subhalo Mass )",alignment=0.5,/normal
  fsc_text,0.04,0.5,"Normalized Fraction",alignment=0.5,orientation=90,/normal
  
  !p.multi = 0
  end_PS

end

; plotTmaxHisto():

pro plotTmaxHisto, res=res

  if not keyword_set(res) then begin
    print,'Error: Must specific resolution set.'
    return
  endif
  
  ; config
  dataPath  = '/n/hernquistfs1/dnelson/coldflows/thermhist/'
  plotPath     = '/n/home07/dnelson/coldflows/'
  
  smoothAccretionOnly = 1
  
  redshiftBins = [6.0,5.0,4.0,3.0,2.0,1.5,1.0,0.5,0.25,0.0]
  
  ; maxt_histo_allbins
  start_PS,plotPath+'maxt_histo_allbins_'+str(res)+'.eps',xs=7,ys=6
  !p.multi = [0,3,3]
  
  for j=0,n_elements(redshiftBins)-2 do begin
    targetRedshift = redshiftBins[j+1]
    endingRedshift = redshiftBins[j]
    
    ; make halo selection at target redshift
    targetSnap = redshiftToSnapNum(targetRedshift)
    endingSnap = redshiftToSnapNum(endingRedshift)
  
    saveFilename = dataPath + 'thermhist.'+str(res)+'.m='+str(targetSnap)+'.to.'+str(endingSnap)+'.sav'
    
    if (keyword_set(smoothAccretionOnly)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.sAO.sav'
    
    if not (file_test(saveFileName)) then begin
      print,'WARNING: Cannot load savefile: '+saveFilename
      continue
    endif
  
    restore, saveFilename
    print,'Working: '+saveFilename
    
    ; calc max temp
    maxTemps = fltarr(n_elements(sgIDs_All))
    for i=0,n_elements(sgIDs_All)-1 do begin
      maxTemps[i] = alog10(max(temp[i,*]))
    endfor
  
    ; plot
    xrange = [3.8,7.2]
    yrange = [0.0,1.05]
    
    bin = 0.1
 
    cs = !p.charsize
    pt = !p.thick
    !p.charsize += 1.0
    !p.thick += 1.0
    
    zstrs = 'z='+['5.0','4.0','3.0','2.0','1.5','1.0','0.5','0.25','0.0']
 
    if (j eq 0) then $
      plothist,maxTemps,bin=bin,xtitle="",ytitle="",/peak,$
               xrange=xrange,yrange=yrange,xtickname=replicate(' ',10),/xs,/ys,$
               xmargin=[7.0,-3.0],ymargin=[0.0,3.0],$
               ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4
    if (j eq 1) then $
      plothist,maxTemps,bin=bin,xtitle="",ytitle="",/peak,$
               xrange=xrange,yrange=yrange,$
               xtickname=replicate(' ',10),ytickname=replicate(' ',10),/xs,/ys,$
               xmargin=[3.0,1.0],ymargin=[0.0,3.0]
    if (j eq 2) then $
      plothist,maxTemps,bin=bin,xtitle="",ytitle="",/peak,$
               xrange=xrange,yrange=yrange,$
               xtickname=replicate(' ',10),ytickname=replicate(' ',10),/xs,/ys,$
               xmargin=[-1.0,5.0],ymargin=[0.0,3.0]
    if (j eq 3) then $
      plothist,maxTemps,bin=bin,xtitle="",ytitle="",/peak,$
               xrange=xrange,yrange=yrange,xtickname=replicate(' ',10),/xs,/ys,$
               xmargin=[7.0,-3.0],ymargin=[2.0,0.0],$
               ytickv=[0.2,0.4,0.6,0.8,1.0],yticks=4
    if (j eq 4) then $
      plothist,maxTemps,bin=bin,xtitle="",ytitle="",/peak,$
               xrange=xrange,yrange=yrange,$
               xtickname=replicate(' ',10),ytickname=replicate(' ',10),/xs,/ys,$
               xmargin=[3.0,1.0],ymargin=[2.0,0.0]
    if (j eq 5) then $
      plothist,maxTemps,bin=bin,xtitle="",ytitle="",/peak,$
               xrange=xrange,yrange=yrange,$
               xtickname=replicate(' ',10),ytickname=replicate(' ',10),/xs,/ys,$
               xmargin=[-1.0,5.0],ymargin=[2.0,0.0]   
    if (j eq 6) then $
      plothist,maxTemps,bin=bin,xtitle="",ytitle="",/peak,$
               xrange=xrange,yrange=yrange,$
               /xs,/ys,xtickv=[4.0,5.0,6.0,7.0],xtickname=['4','5','6','7'],xticks=3,$
               xmargin=[7.0,-3.0],ymargin=[5.0,-2.0]
    if (j eq 7) then $
      plothist,maxTemps,bin=bin,xtitle="",ytitle="",/peak,$
               xrange=xrange,yrange=yrange,/xs,/ys,ytickname=replicate(' ',10),$
               xmargin=[3.0,1.0],ymargin=[5.0,-2.0],$
               xtickv=[4.0,5.0,6.0,7.0],xtickname=['4','5','6','7'],xticks=3
    if (j eq 8) then $
      plothist,maxTemps,bin=bin,xtitle="",ytitle="",/peak,$
               xrange=xrange,yrange=yrange,ytickname=replicate(' ',10),/xs,/ys,$
               xmargin=[-1.0,5.0],ymargin=[5.0,-2.0],$
               xtickv=[4.0,5.0,6.0,7.0],xtickname=['4','5','6','7'],xticks=3
               
    !p.charsize = cs
    !p.thick = pt
    
    if (j lt 6) then begin & textx = 5.8 & texty = 0.76 & endif
    if (j ge 6) then begin & textx = 5.1 & texty = 0.10 & endif
    fsc_text,textx,texty,zstrs[j],color=fsc_color('orange')

  endfor ;j
  
  fsc_text,0.5,0.95,str(res)+"^3",alignment=0.5,/normal
  fsc_text,0.5,0.04,"log (max T [K])",alignment=0.5,/normal
  fsc_text,0.04,0.5,"Normalized Fraction",alignment=0.5,orientation=90,/normal

  !p.multi = 0
  end_PS
  
end

; plotTempTracks()

pro plotTempTracks, res=res

  if not keyword_set(res) then begin
    print,'Error: Must specific resolution set.'
    return
  endif
  
  ; config
  dataPath  = '/n/hernquistfs1/dnelson/coldflows/thermhist/'
  plotPath  = '/n/home07/dnelson/coldflows/'
  
  smoothAccretionOnly = 1
  
  redshiftBins = [6.0,5.0,4.0,3.0,2.0,1.5,1.0,0.5,0.25,0.0]
 
  ; load data
  j = 5 ;5,3,2
  targetRedshift = redshiftBins[j+1]
  endingRedshift = redshiftBins[j]
  
  targetSnap = redshiftToSnapNum(targetRedshift)
  endingSnap = redshiftToSnapNum(endingRedshift)

  saveFilename = dataPath + 'thermhist.'+str(res)+'.m='+str(targetSnap)+'.to.'+str(endingSnap)+'.sav'
  
  if (keyword_set(smoothAccretionOnly)) then $
  saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.sAO.sav'
  
  if not (file_test(saveFileName)) then begin
    print,'Error: Cannot load savefile: '+saveFilename
    return
  endif
  
  restore, saveFilename

  ; calc max temp
  maxTemps = fltarr(n_elements(sgIDs_All))
  for i=0,n_elements(sgIDs_All)-1 do begin
    maxTemps[i] = alog10(max(temp[i,*]))
  endfor

  ; find 100 tracks with Tmax<critLogTemp-1.0 (blue) and 100 with Tmax>critLogTemp+1.0 (red)
  units = getUnits()
  
  critLogTemp = 5.5
  
  numTracksEach = 20
  
  ; select
  wBlue = where(maxTemps lt critLogTemp-1.0, countBlue)
  wRed  = where(maxTemps gt critLogTemp+0.5, countRed)
  
  idBlue = randomu(seed,countBlue)
  idBlue = (sort(idBlue))[0:numTracksEach-1]
  idBlue = wBlue[idBlue]
  
  idRed = randomu(seed,countRed)
  idRed = (sort(idRed))[0:numTracksEach-1]
  idRed = wRed[idRed] 
  
  ;if 0 then begin
  
  ; plot axes
  redshifts = snapNumToRedshift(/all)
  times     = redshiftToAge(redshifts)
  
  xrange = [0,times[targetSnap]] ;targetSnap-1
  yrange = [1e1,4e6]
  
  ; temp vs. time/redshift tracks
  start_PS, plotPath + 'temp.time.tracks.z='+string(redshiftBins[j+1],format='(f4.2)')+'.eps'
    fsc_plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,$
             xtitle="Time [Gyr]",ytitle="Temperature [K]",xs=9,/ys,/ylog,ymargin=[4,3]
    redshift_axis,xrange,yrange,/ylog
    
    ; individual tracks
    for pID=0,numTracksEach-1 do begin
      fsc_plot,times,temp[idBlue[pID],*],/overplot,thick=0.1,color=fsc_color('blue')
      fsc_plot,times,temp[idRed[pID],*],/overplot,thick=0.1,color=fsc_color('red')
    endfor
    
    fsc_text,times[targetSnap]*0.8,3e2,"log "+textoidl("T_{max}")+" < 4.5",$
              alignment=0.5,color=fsc_color('blue')
    fsc_text,times[targetSnap]*0.8,1e2,"log "+textoidl("T_{max}")+" > 6.0",$
              alignment=0.5,color=fsc_color('red')
  end_PS
  
  ;endif ;0
  
  ; rho/temp plane tracks
  start_PS, plotPath + 'rho.temp.tracks.'+string(redshiftBins[j+1],format='(f4.2)')+'.eps'
    fsc_plot,[0],[0],/nodata,xrange=[0.0,4.0],yrange=[3.0,7.0],$
             xtitle="log "+textoidl("\rho / \rho_{crit}"),ytitle="log Temperature [K]",/xs,/ys
            
    plotsym,0
            
    ; individual tracks
    for pID=0,numTracksEach-1 do begin
      ;w = where(density[idBlue[pID],*] ne 0.0 and temp[idBlue[pID],*] ne 0.0, countGood)

      fsc_plot,alog10(density[idBlue[pID],*]/units.rhoCrit),$
               alog10(temp[idBlue[pID],*]),/overplot,$
               thick=0.5,psym=-8,symsize=0.2,color=fsc_color('blue')
      fsc_plot,alog10(density[idRed[pID],*]/units.rhoCrit),$
               alog10(temp[idRed[pID],*]),/overplot,$
               thick=0.5,psym=-8,symsize=0.2,color=fsc_color('red')
    endfor
    
  end_PS

  stop

end

; selectFilament():

pro selectFilament

  units = getUnits()

  ; config
  targetSnap  = 189
  subhaloPath = '/n/hernquistfs1/mvogelsberger/ComparisonProject/512_20Mpc/Arepo_ENERGY/output/'
  workingPath = '/n/home07/dnelson/coldflows/vis/cutout_512A/'
  
  ; filament description
  x1 = [1023.2,7468.8,16090.0] ;start
  x2 = [1110.0,7540.0,16130.0] ;end  
  
  filRad   = 10.0      ; kpc
  filBound = [0.0,1.0] ; t
  
  ; 200kpc and 400kpc bboxes
  bbox    = [1023.20,1223.20,7468.80,7668.80,16044.2,16244.2]/1e5 ;Mpc
  bboxBig = [923.2,1323.2,7368.8,7768.8,15944.2,16344.2]/1e5 ;Mpc
  
  ; load subhalo cat and select in bbox
  sh = loadSubhaloGroups(subhaloPath,targetSnap,/verbose)
  
  wGR = where(sh.groupPos[0,*] gt bboxBig[0]*1e5 and sh.groupPos[0,*] lt bboxBig[1]*1e5 and $
              sh.groupPos[1,*] gt bboxBig[2]*1e5 and sh.groupPos[1,*] lt bboxBig[3]*1e5 and $
              sh.groupPos[2,*] gt bboxBig[4]*1e5 and sh.groupPos[2,*] lt bboxBig[5]*1e5, countGR)
  wSH = where(sh.subgroupPos[0,*] gt bboxBig[0]*1e5 and sh.subgroupPos[0,*] lt bboxBig[1]*1e5 and $
              sh.subgroupPos[1,*] gt bboxBig[2]*1e5 and sh.subgroupPos[1,*] lt bboxBig[3]*1e5 and $
              sh.subgroupPos[2,*] gt bboxBig[4]*1e5 and sh.subgroupPos[2,*] lt bboxBig[5]*1e5, countSH)
  print,'Found ['+str(countGR)+'] groups and ['+str(countSH)+'] subgroups in big bbox.'

  ; temp: 3d bbox pts for plotting
  bboxPts = [[bbox[0],bbox[2],bbox[4]],$
             [bbox[0],bbox[3],bbox[4]],$
             [bbox[1],bbox[3],bbox[4]],$
             [bbox[1],bbox[2],bbox[4]],$
             [bbox[0],bbox[2],bbox[4]],$
             [bbox[0],bbox[2],bbox[5]],$
             [bbox[0],bbox[3],bbox[5]],$
             [bbox[1],bbox[3],bbox[5]],$
             [bbox[1],bbox[2],bbox[5]],$
             [bbox[0],bbox[2],bbox[5]]]
  
  ; load particle positions
  x0_gas   = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='gas',field='pos',/verbose)
  x0_stars = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='stars',field='pos',/verbose) 
  x0_dm    = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='dm',field='pos',/verbose)
  
  ; parametric solution and minimum distance to line (gas)  
  n21  = (x2[0]-x1[0])^2.0 + (x2[1]-x1[1])^2.0 + (x2[2]-x1[2])^2.0 ;n21=transpose(x1-x2)#(x2-x1)
  n10  = reform( (x1[0]-x0_gas[0,*])^2.0 + (x1[1]-x0_gas[1,*])^2.0 + (x1[2]-x0_gas[2,*])^2.0 )
  dotp = reform( (x1[0]-x0_gas[0,*])*(x2[0]-x1[0]) + (x1[1]-x0_gas[1,*])*(x2[1]-x1[1]) + $
                 (x1[2]-x0_gas[2,*])*(x2[2]-x1[2]) )

  t_gas = -1.0 * dotp / n21
  d2    = ( n10 * n21 - dotp^2.0 ) / n21
  d_gas = sqrt(d2)
  
  ; (stars)
  n21  = (x2[0]-x1[0])^2.0 + (x2[1]-x1[1])^2.0 + (x2[2]-x1[2])^2.0 ;n21=transpose(x1-x2)#(x2-x1)
  n10  = reform( (x1[0]-x0_stars[0,*])^2.0 + (x1[1]-x0_stars[1,*])^2.0 + (x1[2]-x0_stars[2,*])^2.0 )
  dotp = reform( (x1[0]-x0_stars[0,*])*(x2[0]-x1[0]) + (x1[1]-x0_stars[1,*])*(x2[1]-x1[1]) + $
                 (x1[2]-x0_stars[2,*])*(x2[2]-x1[2]) )

  t_stars = -1.0 * dotp / n21
  d2      = ( n10 * n21 - dotp^2.0 ) / n21
  d_stars = sqrt(d2)
  
  ; (dm)
  n21  = (x2[0]-x1[0])^2.0 + (x2[1]-x1[1])^2.0 + (x2[2]-x1[2])^2.0 ;n21=transpose(x1-x2)#(x2-x1)
  n10  = reform( (x1[0]-x0_dm[0,*])^2.0 + (x1[1]-x0_dm[1,*])^2.0 + (x1[2]-x0_dm[2,*])^2.0 )
  dotp = reform( (x1[0]-x0_dm[0,*])*(x2[0]-x1[0]) + (x1[1]-x0_dm[1,*])*(x2[1]-x1[1]) + $
                 (x1[2]-x0_dm[2,*])*(x2[2]-x1[2]) )

  t_dm = -1.0 * dotp / n21
  d2   = ( n10 * n21 - dotp^2.0 ) / n21
  d_dm = sqrt(d2)
  
  ; make selection
  wGas   = where(t_gas   gt filBound[0] and d_gas   lt filRad and t_gas   lt filBound[1], countGas)
  wStars = where(t_stars gt filBound[0] and d_stars lt filRad and t_stars lt filBound[1], countStars)
  wDM    = where(t_dm    gt filBound[0] and d_dm    lt filRad and t_dm    lt filBound[1], countDM)
  
  print,'Filament selection: '+str(countGas)+ ' gas, '+str(countStars)+' stars, '+str(countDM)+' DM.'
  
  ; plot dist from line vs t
  start_PS,workingPath+'dist_vs_t.eps',xs=7,ys=10
    ; gas
    fsc_plot,[0],[0],/nodata,xrange=[0,40],yrange=[-1.4,2.0],$
             xtitle="Distance from Line",ytitle="Parameter t (Dist along Line)",/xs,/ys,$
             title="Gas",position=[0.1,0.1,0.3,0.9],charsize=!p.charsize-1
    fsc_plot,d_gas,t_gas,psym=3,symsize=1.0,/overplot
    fsc_plot,d_gas[wGas],t_gas[wGas],psym=3,symsize=2.0,color=fsc_color('green'),/overplot

    ; stars
    fsc_plot,[0],[0],/nodata,xrange=[0,40],yrange=[-1.4,2.0],$
             xtitle="Distance from Line",ytitle="Parameter t (Dist along Line)",/xs,/ys,$
             title="Stars",position=[0.4,0.1,0.6,0.9],charsize=!p.charsize-1,/noerase
    fsc_plot,d_stars,t_stars,psym=3,symsize=1.0,/overplot
    fsc_plot,d_stars[wStars],t_stars[wStars],psym=3,symsize=2.0,color=fsc_color('green'),/overplot
    
    ; dm
    fsc_plot,[0],[0],/nodata,xrange=[0,40],yrange=[-1.4,2.0],$
             xtitle="Distance from Line",ytitle="Parameter t (Dist along Line)",/xs,/ys,$
             title="DM",position=[0.7,0.1,0.9,0.9],charsize=!p.charsize-1,/noerase
    fsc_plot,d_dm,t_dm,psym=3,symsize=1.0,/overplot
    fsc_plot,d_dm[wDM],t_dm[wDM],psym=3,symsize=2.0,color=fsc_color('green'),/overplot
  end_PS
  
  ; plot selection (gas)
  start_PS, workingPath+'gas_select_xyz.eps', xs=8.0, ys=8.0
    !p.multi = [0,2,2]
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[0],bboxBig[1]],yrange=[bboxBig[2],bboxBig[3]],$
               xtitle="x [Mpc]",ytitle="y [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_gas[0,*]/1e5,x0_gas[1,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_gas[0,wGas]/1e5,x0_gas[1,wGas]/1e5,psym=3,symsize=3.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[0],x2[0]]/1e5,[x1[1],x2[1]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      ;for i=0,countGR-1 do begin
      ;  print,sh.grouppos[0,wGR[i]]/1e5,sh.groupPos[1,wGR[i]]/1e5,sh.group_R_Mean200[wGR[i]]/1e5
      ;  tvcircle,sh.group_r_mean200[wGR[i]]/1e5,sh.grouppos[0,wGR[i]]/1e5,sh.grouppos[1,wGR[i]]/1e5,$
      ;           color=fsc_color('blue'),/data,thick=!p.thick-1.0
      ;endfor
      for i=0,countSH-1 do begin
        ;print,sh.subgroupMass[wSH[i]]/5e4,sh.subgroupPos[0,wSH[i]]/1e5,sh.subgroupPos[1,wSH[i]]/1e5
        tvcircle,sh.subgroupMass[wSH[i]]/5e4,sh.subgroupPos[0,wSH[i]]/1e5,sh.subgroupPos[1,wSH[i]]/1e5,$
                 color=fsc_color('purple'),/data,thick=!p.thick-1.0
      endfor
      
      fsc_plot,[bbox[0],bbox[1],bbox[1],bbox[0],bbox[0]],$
               [bbox[2],bbox[2],bbox[3],bbox[3],bbox[2]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
      
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[2],bboxBig[3]],yrange=[bboxBig[4],bboxBig[5]],$
               xtitle="y [Mpc]",ytitle="z [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_gas[1,*]/1e5,x0_gas[2,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_gas[1,wGas]/1e5,x0_gas[2,wGas]/1e5,psym=3,symsize=1.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[1],x2[1]]/1e5,[x1[2],x2[2]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      for i=0,countSH-1 do begin
        tvcircle,sh.subgroupMass[wSH[i]]/5e4,sh.subgroupPos[1,wSH[i]]/1e5,sh.subgroupPos[2,wSH[i]]/1e5,$
                 color=fsc_color('purple'),/data,thick=!p.thick-1.0
      endfor
      
      fsc_plot,[bbox[2],bbox[3],bbox[3],bbox[2],bbox[2]],$
               [bbox[4],bbox[4],bbox[5],bbox[5],bbox[4]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
      
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[0],bboxBig[1]],yrange=[bboxBig[4],bboxBig[5]],$
               xtitle="x [Mpc]",ytitle="z [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_gas[0,*]/1e5,x0_gas[2,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_gas[0,wGas]/1e5,x0_gas[2,wGas]/1e5,psym=3,symsize=1.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[0],x2[0]]/1e5,[x1[2],x2[2]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      for i=0,countSH-1 do begin
        tvcircle,sh.subgroupMass[wSH[i]]/5e4,sh.subgroupPos[0,wSH[i]]/1e5,sh.subgroupPos[2,wSH[i]]/1e5,$
                 color=fsc_color('purple'),/data,thick=!p.thick-1.0
      endfor
      
      fsc_plot,[bbox[0],bbox[1],bbox[1],bbox[0],bbox[0]],$
               [bbox[4],bbox[4],bbox[5],bbox[5],bbox[4]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
               
    !p.multi = 0
  end_PS
  
  ; (stars)
  start_PS, workingPath+'stars_select_xyz.eps', xs=8.0, ys=8.0
    !p.multi = [0,2,2]
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[0],bboxBig[1]],yrange=[bboxBig[2],bboxBig[3]],$
               xtitle="x [Mpc]",ytitle="y [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_stars[0,*]/1e5,x0_stars[1,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_stars[0,wStars]/1e5,x0_stars[1,wStars]/1e5,psym=3,symsize=3.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[0],x2[0]]/1e5,[x1[1],x2[1]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      fsc_plot,[bbox[0],bbox[1],bbox[1],bbox[0],bbox[0]],$
               [bbox[2],bbox[2],bbox[3],bbox[3],bbox[2]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
      
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[2],bboxBig[3]],yrange=[bboxBig[4],bboxBig[5]],$
               xtitle="y [Mpc]",ytitle="z [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_stars[1,*]/1e5,x0_stars[2,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_stars[1,wStars]/1e5,x0_stars[2,wStars]/1e5,psym=3,symsize=1.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[1],x2[1]]/1e5,[x1[2],x2[2]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      fsc_plot,[bbox[2],bbox[3],bbox[3],bbox[2],bbox[2]],$
               [bbox[4],bbox[4],bbox[5],bbox[5],bbox[4]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
      
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[0],bboxBig[1]],yrange=[bboxBig[4],bboxBig[5]],$
               xtitle="x [Mpc]",ytitle="z [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_stars[0,*]/1e5,x0_stars[2,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_stars[0,wStars]/1e5,x0_stars[2,wStars]/1e5,psym=3,symsize=1.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[0],x2[0]]/1e5,[x1[2],x2[2]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      fsc_plot,[bbox[0],bbox[1],bbox[1],bbox[0],bbox[0]],$
               [bbox[4],bbox[4],bbox[5],bbox[5],bbox[4]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
               
    !p.multi = 0
  end_PS
  
  ; (dm)
  start_PS, workingPath+'dm_select_xyz.eps', xs=8.0, ys=8.0
    !p.multi = [0,2,2]
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[0],bboxBig[1]],yrange=[bboxBig[2],bboxBig[3]],$
               xtitle="x [Mpc]",ytitle="y [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_dm[0,*]/1e5,x0_dm[1,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_dm[0,wDM]/1e5,x0_dm[1,wDM]/1e5,psym=3,symsize=3.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[0],x2[0]]/1e5,[x1[1],x2[1]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      fsc_plot,[bbox[0],bbox[1],bbox[1],bbox[0],bbox[0]],$
               [bbox[2],bbox[2],bbox[3],bbox[3],bbox[2]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
      
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[2],bboxBig[3]],yrange=[bboxBig[4],bboxBig[5]],$
               xtitle="y [Mpc]",ytitle="z [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_dm[1,*]/1e5,x0_dm[2,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_dm[1,wDM]/1e5,x0_dm[2,wDM]/1e5,psym=3,symsize=1.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[1],x2[1]]/1e5,[x1[2],x2[2]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      fsc_plot,[bbox[2],bbox[3],bbox[3],bbox[2],bbox[2]],$
               [bbox[4],bbox[4],bbox[5],bbox[5],bbox[4]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
      
      fsc_plot,[0],[0],/nodata,xrange=[bboxBig[0],bboxBig[1]],yrange=[bboxBig[4],bboxBig[5]],$
               xtitle="x [Mpc]",ytitle="z [Mpc]",/xs,/ys,charsize=1.0
      fsc_plot,x0_dm[0,*]/1e5,x0_dm[2,*]/1e5,psym=3,symsize=1.0,/overplot
      fsc_plot,x0_dm[0,wDM]/1e5,x0_dm[2,wDM]/1e5,psym=3,symsize=1.0,color=fsc_color('green'),/overplot
      fsc_plot,[x1[0],x2[0]]/1e5,[x1[2],x2[2]]/1e5,line=0,/overplot,color=fsc_color('red')
      
      fsc_plot,[bbox[0],bbox[1],bbox[1],bbox[0],bbox[0]],$
               [bbox[4],bbox[4],bbox[5],bbox[5],bbox[4]],$
               line=0,color=fsc_color('orange'),thick=!p.thick-1,/overplot
               
    !p.multi = 0
  end_PS
  
  ; load gas properties
  dens = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='gas',field='density')
  u     = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='gas',field='u')
  nelec = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='gas',field='ne')
  mass = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='gas',field='masses')
  
  dens = dens[wGas]
  u     = u[wGas]
  nelec = nelec[wGas]
  mass = mass[wGas]

  temp = convertUtoTemp(u,nelec)
  entropy = calcEntropy(u,dens)
  
  ; plot gas properties
  start_PS, workingPath+'gas_properties.eps'
    !p.multi = [0,2,2]
       plothist,alog10(rhoRatioToCrit(dens)),/auto,$
                xtitle="log ("+textoidl("\rho / \rho_{crit}")+")",ytitle="N",$
                title="Gas Density",charsize=!p.charsize-1
       
       plothist,alog10(temp),/auto,xtitle="log (T [K])",ytitle="N",$
                title="Gas Temperature",charsize=!p.charsize-1
                
       h2rt = rhoTHisto(dens,temp,mass=mass,/plot)
       
       plothist,alog10(entropy),/auto,xtitle="log (Entropy) [Code]",ytitle="N",$
                title="Gas Entropy",charsize=!p.charsize-1
    !p.multi = 0
  end_PS
  
  ; velocities - find distribution of cos(theta) between line and individual vel vectors
  vels = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='gas',field='vel')
  vels = vels[*,wGas]

  line = x2-x1

  dotp  = reform( line[0]*vels[0,*] + line[1]*vels[1,*] + line[2]*vels[2,*] )
  norm1 = line[0]*line[0] + line[1]*line[1] + line[2]*line[2]
  norm2 = reform( (vels[0,*]*vels[0,*] + vels[1,*]*vels[1,*] + vels[2,*]*vels[2,*]) )

  cos_theta = dotp / sqrt(norm1*norm2)
  theta_deg_gas = acos(cos_theta) * 180.0/!pi
  
  ; (dm)
  vels = loadSnapshotSubset(workingPath,snapNum=targetSnap,partType='dm',field='vel')
  vels = vels[*,wDM]

  line = x2-x1

  dotp  = reform( line[0]*vels[0,*] + line[1]*vels[1,*] + line[2]*vels[2,*] )
  norm1 = line[0]*line[0] + line[1]*line[1] + line[2]*line[2]
  norm2 = reform( (vels[0,*]*vels[0,*] + vels[1,*]*vels[1,*] + vels[2,*]*vels[2,*]) )

  cos_theta = dotp / sqrt(norm1*norm2)
  theta_deg_dm = acos(cos_theta) * 180.0/!pi

  start_PS, workingPath+'select_vel_dir.eps'
      plothist,theta_deg_gas,/auto,xtitle="Relative Angle between Part. Vel and Filament [deg]",ytitle="N",$
               xrange=[-2.0,180.0],yrange=[0,1.05],/peak,/ys,/xs
      plothist,theta_deg_dm,/auto,/overplot,/peak,color=fsc_color('orange')
      fsc_text,150.0,0.9,"Gas",color=fsc_color('black')
      fsc_text,150.0,0.84,"DM",color=fsc_color('orange')
  end_PS
  
  
  
  stop
end

