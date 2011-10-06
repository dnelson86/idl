; coldflowsPlot.pro
; cold flows - plots derived from intermediate datafiles
; dnelson oct.2011

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
  includeBackgroundSH = 1
  
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
    if (keyword_set(includeBackgroundSH)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.bSH.sav'
    
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
  includeBackgroundSH = 1
  
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
    if (keyword_set(includeBackgroundSH)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.bSH.sav'
    
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
  includeBackgroundSH = 1
  
  massWeighted = 0
  
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
    if (keyword_set(includeBackgroundSH)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.bSH.sav'
    
    if not (file_test(saveFileName)) then begin
      print,'WARNING: Cannot load savefile: '+saveFilename
      continue
    endif
  
    restore, saveFilename,/verbose
    print,'Working: '+saveFilename

    ; calc max temp
    maxTemps = fltarr(n_elements(sgIDs_All))
    for i=0,n_elements(sgIDs_All)-1 do begin
      maxTemps[i] = alog10(max(temp[i,*]))
    endfor
    
    ; calc total masses
    totMasses = fltarr(n_elements(sgIDs_All))
    if (massWeighted eq 1) then begin
      for i=0,n_elements(sgIDs_All)-1 do begin
        totMasses[i] = 0
      endfor
    endif
  stop
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
  includeBackgroundSH = 1
  
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
  if (keyword_set(includeBackgroundSH)) then $
  saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.bSH.sav'
  
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

end
