; spiralsLC.pro
; "lambda crit" project specific (correlations and multi-sim comparisons)
; dnelson may.2011

; compAmpGrowth(): compare growth of structure (without perturbers) at different res

pro compAmpGrowth

  units = getUnits()

  ;simNames   = ['LC_1m_disk_8_nomc','LC_5m_disk_8_nomc','LC_10m_disk_8_nomc',$
  ;              'LC_50m_disk_8_nomc','LC_100m_disk_8_nomc',$
  ;              'LC_1m_disk_8_nomc_e','LC_5m_disk_8_nomc_e','LC_50m_disk_8_nomc_e']
  ;simNames = ['LC_10m_disk_2_nomc','LC_10m_disk_4_nomc','LC_10m_disk_6_nomc','LC_10m_disk_8_nomc']
  simNames = ['LC_10m_disk_2','LC_10m_disk_2_nomc','LC_10m_disk_2_mcsoff']

  ; path config
  fileBase    = '/n/home07/dnelson/spirals/lambdacrit/sims/'
  filePaths   = fileBase + simNames + '/'
  workingPath = '/n/home07/dnelson/spirals/lambdacrit/'

  nComps = n_elements(filePaths)
  
  ; config
  nRad    = 360.0 ;number for azimuthal histogram
  widthR  = 0.5   ;kpc  
  numSLs  = 1.0   ;times scale length
    
  ; compute overdensity at targetR for all sims
  for i=0,nComps-1 do begin

    h2 = loadSimParams(simNames[i])
    targetR = numSLs * h2.h

    saveFileName = filePaths[i] + 'ampGrowth.sav'
    
    if not (file_test(saveFileName)) then begin
  
      nSnaps = max([n_elements(file_search(filePaths[i]+"output/Stars_*.bin")),$
                    n_elements(file_search(filePaths[i]+"output/Stars_*_0"))])
   
      ; setup arrays
      annulusA   = fltarr(nSnaps)
      times      = fltarr(nSnaps)
   
      for m=0,nSnaps-1 do begin
      
        if not (file_test(filePaths[i]+"output/MCs_"+str(m+1)+".txt")) then continue
      
        ; load
        xyz      = loadStars(filePaths[i],(m+1))
        mcs      = loadMCs(filePaths[i],(m+1),h)
        times[m] = h[1] * 1000.0 ;Myr
      
        ; select annulus
        r = reform(sqrt(xyz.x^2.0 + xyz.y^2.0))
        w = where(r ge (targetR-widthR/2.0) and r lt (targetR+widthR/2.0),count)
        
        if (count ne 0) then begin
          arctan,xyz[w].x,xyz[w].y,theta_rad,theta_deg ;[0,2pi]
          nRad = round(count / 1000.0)
          hist1d = histogram(theta_rad,nbins=nRad,min=0.0,max=2*!pi)
          annulusA[m] = max(hist1d) / mean(hist1d)
          print,i,m,times[m],annulusA[m],count,count/nRad,1.0/sqrt(count/nRad)
        endif
      
      endfor ;m

      ; save/restore
      save,numSLs,nRad,nSnaps,widthR,annulusA,times,filename=saveFileName
      print,'saved ',saveFileName
    endif else begin
      restore,saveFileName
    endelse
  
  endfor ;i

  ; plot
  colors = ['blue','forest green','black','red','orange','blue','forest green','red']

  start_PS, workingPath+"compAmpGrowth.eps"
  
    fsc_plot,[0],[0],ytitle="Maximum Surface Overdensity (at h)",xtitle="Time [Myr]",$
             xrange=[0,3000.0],/xs,yrange=[1.0,2.0],/ys,/nodata
             
    for i=0,nComps-1 do begin
    
      ; restore save
      saveFileName = filePaths[i] + 'ampGrowth.sav'
      restore,saveFileName
      
      w = where(times ne 0,count)
      
      if (count eq 0) then continue
     
      ;parse simName
      ;shortName = strmid(simNames[i],3,strpos(simNames[i],"_",3)-3) ;LC8
      ;shortName = strmid(simNames[i],7,6) ;LCs
      shortName = strmid(simNames[i],7) ;LC2 3Gyr

      line  = 0
      if (strpos(simNames[i],"_e") ne -1) then begin
        line  = 1
        shortName = shortName + "_e"
      endif
          
      ;smoothed
      fsc_plot, times[w], smooth(annulusA[w],nSnaps/75.0), line=line, thick=3.0-line, $
                color=fsc_color(colors[i]), /overplot
                
      ; legend
      legendXStart  = 150.0 ;0.75 (right side)
      legendXLength = 100.0
      legendYStart  = 1.9
      legendYSpace  = 0.05
      
      fsc_plot,[legendXStart,legendXStart+legendXLength],$
               [legendYStart-legendYSpace*i,legendYStart-legendYSpace*i],$
               color=fsc_color(colors[i]),line=line,thick=2.0,/overplot
      fsc_text,legendXStart+legendXLength*1.2,legendYStart-(legendYSpace/5)-legendYSpace*i,$
               shortName,charsize=1.0,alignment=0.0,/data
                  
    endfor
  
  end_PS

end

; compFourierWithTime

function compFourierWithTime, simNames, plot=plot

  units = getUnits()

  ; config
  fileBase = '/n/home07/dnelson/spirals/lambdacrit/sims/'
  ;simNames = ['LC_10m_disk_2','LC_10m_disk_4','LC_10m_disk_6','LC_10m_disk_8']
  filePaths = fileBase + simNames + '/'
  workingPath = '/n/home07/dnelson/spirals/lambdacrit/'

  nComps = n_elements(filePaths)
  
  nModes      = 8
  numSLs      = 2.0 ; times scale length, bar=0.5, spiral=2.0
  
  ; find number of saved fourier
  nSnaps  = 1
  targetR = fltarr(nComps)
  
  for i=0,nComps-1 do begin
    nSnaps = max([nSnaps,n_elements(file_search(filePaths[i]+"output2/fourierSaved_*.sav"))])
    
    h2 = loadSimParams(filePaths[i])
    targetR[i] = numSLs * h2.h
  endfor
  
  ; setup arrays
  fourierA   = fltarr(nComps,nSnaps,nModes+1)
  fourierMax = fltarr(nComps,nSnaps)
  times      = fltarr(nComps,nSnaps) 

  ; load pre-computed fourier results
  for i=0,nComps-1 do begin
    nSnaps = n_elements(file_search(filePaths[i]+"output2/fourierSaved_*.sav"))
  
    for j=0,nSnaps-1,1 do begin
    
     fileName = filePaths[i] + "output2/fourierSaved_" + str(j+1) + ".sav"
     restore,fileName   
    
      ; find annulus closest to target radius
      w = where(abs(f.R - targetR[i]) eq min(abs(f.R - targetR[i])))
      A = reform(f.A[w,*])
      
      ; pick maximum mode
      w2 = where(A[1:*] eq max(A[1:*]))

      ; save mode amplitudes
      fourierA[i,j,*] = A
      
      ; save max mode number and time
      fourierMax[i,j] = w2[0]+1
      times[i,j] = time ;Myr

    endfor
  endfor

  ; plot
  if keyword_set(plot) then begin
  start_PS, workingPath+"fourierTimeMax.eps"
    fsc_plot,[0],[0],ytitle="Maximum Fourier Mode",xtitle="Time [Myr]",$
             xrange=[0,round(max(times)/100)*100],/xs,yrange=[0.5,nModes+0.5],/ys,/nodata
             
    ; legend
    legendXStart  = 100.0
    legendXLength = 100.0
    legendYStart  = 1.5
    legendYSpace  = 0.3
             
    for i=0,nComps-1 do begin
      w = where(times[i,*] ne 0)
     
      ;original
      ;fsc_plot, times[i,w], fourierMax[i,w]+0.01*i, line=0, $
      ;          color=fsc_color(units.colors[i]), thick=!p.thick-2, /overplot
                
      ;smoothed
      fsc_plot, times[i,w], smooth(fourierMax[i,w]+0.01*i,nSnaps/50), line=0, $
                color=fsc_color(units.colors[i]), /overplot
                
      ;legend
      shortName = strmid(simNames[i],7)
      fsc_plot,[legendXStart,legendXStart+legendXLength],$
               [legendYStart-legendYSpace*i,legendYStart-legendYSpace*i],$
               color=fsc_color(units.colors[i]),line=line,thick=2.0,/overplot
      fsc_text,legendXStart+legendXLength*1.2,legendYStart-(legendYSpace/5)-legendYSpace*i,$
               shortName,charsize=1.0,alignment=0.0,/data
    endfor
  end_PS
  
  ymax = max(fourierA[*,*,1:*]/fourierA[*,*,0]) * 1.1
  ymax = 0.03
  print,ymax
  
  start_PS, workingPath+"fourierTimeModes.eps"
    fsc_plot,[0],[0],ytitle="Fourier Mode Amplitude",xtitle="Time [Myr]",$
             xrange=[0,round(max(times)/100)*100],/xs,yrange=[0,ymax],/ys,/nodata
             
    for i=0,nComps-1 do begin
      w = where(times[i,*] ne 0)
      for j=1,nModes do begin
        fsc_plot, times[i,w], fourierA[i,w,j]/fourierA[i,*,0], color=fsc_color(units.colors[j]), $
                  line=i, /overplot, thick=!p.thick-1.5 ;/max(fourierA[i,*,*])
      endfor
    endfor
             
    ; legend
    legendXStart = 0.05 ;0.75 (right side)
    legendYStart = ymax * 0.98
    legendYSpace = ( ymax )/24.0

    for j=1,nModes do begin
      fsc_plot,[max(times)*(legendXStart),max(times)*(legendXStart+0.08)],$
               [legendYStart-legendYSpace*j,legendYStart-legendYSpace*j],$
               color=fsc_color(units.colors[j]),/overplot
      fsc_text,max(times)*(legendXStart+0.14),legendYStart-(legendYSpace/5.0)-legendYSpace*j,$
               "m = "+str(j),charsize=!p.charsize-0.5,alignment=0.5,/data
    endfor    
  end_PS
  endif ;plot
  
  r = {times:times,fourierA:fourierA,fourierMax:fourierMax}
  
  return,r
end

; LCcorrs(): temp plot everything vs lambda crit

pro LCcorrs

  units = getUnits()
  
  ; config
  workingPath = '/n/home07/dnelson/spirals/lambdacrit/'
  fileBase    = '/n/home07/dnelson/spirals/lambdacrit/sims/'
  simNames    = ['LC_10m_disk_2','LC_10m_disk_4','LC_10m_disk_6','LC_10m_disk_8']
  
  filePaths = fileBase + simNames + '/'
  
  nComps = n_elements(filePaths)
  
  ; sizing
  resN = 2
  
  ; arrays
  barAngMean    = fltarr(nComps)
  barAngMeanErr = fltarr(nComps)
  barFrac       = fltarr(nComps)
  
  rtwModelParams = fltarr(resN,nComps)
  rtwModelErrs   = fltarr(resN,nComps)
  
  binnedTW    = fltarr(nComps)
  binnedTWErr = fltarr(nComps)
  
  fourierMax  = fltarr(nComps)
  
  lambdaCrit  = fltarr(nComps)
  diskSL      = fltarr(nComps)
  
  ; load
  for i=0,nComps-1 do begin
  
    ;IC info
    h2 = loadSimParams(simNames[i])
    
    lambdaCrit[i] = h2.lc
    diskSL[i]     = h2.h
    
    ;bar
    r = barSpeed(simNames[i])
    
    w = where(abs(r.radii - 2.0*h2.h) eq min(abs(r.radii - 2.0*h2.h)),count)
    
    if (count eq 0) then $
      print,'error1'
      
    barAngMean[i]    = r.angmean[w]
    barAngMeanErr[i] = r.angmeanerr[w]
    barFrac[i]       = r.frac[w]
    
    ;RTW model parameters
    restore,filePaths[i]+"pSS.30.sav"

    rtwModelParams[0,i] = mean(omega_p_rad_mpfit[0,*])
    rtwModelParams[1,i] = mean(omega_p_rad_mpfit[1,*])
    rtwModelErrs[0,i] = stddev(omega_p_rad_mpfit[0,*])
    rtwModelErrs[1,i] = stddev(omega_p_rad_mpfit[1,*])
    
    ;binned TW at ~2h (second bin box)
    omegap1 = mean(omega_p_const_rad[2,*])
    omegap2 = mean(omega_p_const_rad[8-2-1,*])
    
    omegae1 = stddev(omega_p_const_rad[2,*])
    omegae2 = stddev(omega_p_const_rad[8-2-1,*])
    
    if (omegap1 le 0.0 or omegap2 le 0.0) then $
      print,'error2'
      
    binnedTW[i]    = (omegap1 + omegap2) / 2.0
    binnedTWErr[i] = (omegae1 + omegae2) / 2.0

    ;fourier
    r = compFourierWithTime([simNames[i]])
    
    fourierMax[i] = median(r.fourierMax)

  endfor
  
  ; plot
  start_PS, workingPath+"LCcorrs.eps"
    !p.multi = [0,2,2]
    
    ; upper left
    fsc_plot,lambdaCrit,binnedTW,xtitle="Lambda Crit (2h) [kpc]",ytitle="Binned TW Omega_P(2h)",$
             psym=16,charsize=1.0,yrange=[min(binnedTW)*0.8,max(binnedTW)*1.2],/ys
             
    for i=0,n_elements(lambdaCrit)-1 do begin
      fsc_plot,[lambdaCrit[i],lambdaCrit[i]],[binnedTW[i]-binnedTWErr[i],binnedTW[i]+binnedTWErr[i]],$
               line=0,/overplot
    endfor
    
    ; upper right
    fsc_plot,lambdaCrit,fourierMax,xtitle="Lambda Crit (2h) [kpc]",ytitle="Max Fourier Mode (w/ time)",$
             psym=16,charsize=1.0,yrange=[2,8],/ys
             
    fsc_plot,lambdaCrit,(2*!pi*diskSL)/lambdaCrit,psym=15,color=fsc_color('orange'),/overplot
  
    ; lower left
    fsc_plot,lambdaCrit,rtwModelParams[0,*],xtitle="Lambda Crit (2h) [kpc]",ytitle="RTW Plaw Amp / 50x Slope",$
             psym=16,charsize=1.0
             
    for i=0,n_elements(lambdaCrit)-1 do begin
      fsc_plot,[lambdaCrit[i],lambdaCrit[i]],$
               [rtwModelParams[0,i]-rtwModelErrs[0,i],rtwModelParams[0,i]+rtwModelErrs[0,i]],$
               line=0,/overplot
    endfor
    
    fsc_plot,lambdaCrit,rtwModelParams[1,*]*50,psym=15,/overplot,color=fsc_color('green')
  
    for i=0,n_elements(lambdaCrit)-1 do begin
      fsc_plot,[lambdaCrit[i],lambdaCrit[i]],$
               [rtwModelParams[1,i]-rtwModelErrs[1,i],rtwModelParams[1,i]+rtwModelErrs[1,i]]*50,$
               line=0,/overplot,color=fsc_color('green')
    endfor
  
    ; lower right
    fsc_plot,lambdaCrit,barAngMean,xtitle="Lambda Crit (2h) [kpc]",ytitle="Direct Overdensity Omega_P(2h)",$
             psym=16,charsize=1.0,yrange=[min(barAngMean)*0.8,max(barAngMean)*1.2],/ys

    for i=0,n_elements(lambdaCrit)-1 do begin
      fsc_plot,[lambdaCrit[i],lambdaCrit[i]],$
               [barAngMean[i]-barAngMeanErr[i],barAngMean[i]+barAngMeanErr[i]],$
               line=0,/overplot
    endfor
  
    !p.multi = 0
  end_PS

  stop
end

; compPatSpeedWithTime()

pro compPatSpeedWithTime
  
  units = getUnits()
  
  ; config
  workingPath = '/n/home07/dnelson/spirals/lambdacrit/'
  fileBase    = '/n/home07/dnelson/spirals/lambdacrit/sims/'
  simNames    = ['LC_10m_disk_2','LC_10m_disk_4','LC_10m_disk_6','LC_10m_disk_8']
  
  filePaths = fileBase + simNames
  
  nComps = n_elements(filePaths)
  
  ; plot
  start_PS, workingPath+"compPS.eps"
            
    fsc_plot,[0.0,1.0],[0.0,50.0],/nodata,$
             xtitle="time [Gyr]",ytitle="globally constant pattern speed [km/s/kpc]"       
            
    for i=0,nComps-1 do begin
  
      restore,filePaths[i]+"pSS.30.sav"
  
      fsc_plot,times,omega_p[0,*],line=i,/overplot,color=fsc_color(units.colors[1]) 
      fsc_plot,times,omega_p[1,*],line=i,/overplot,color=fsc_color(units.colors[2])
      fsc_plot,times,omega_p[2,*],line=i,/overplot,color=fsc_color(units.colors[3]) 
      fsc_plot,times,omega_p[3,*],line=i,/overplot,color=fsc_color(units.colors[4])  
 
      fsc_plot,[0.85,0.89],[7.0-2*i,7.0-2*i],line=i,/overplot
      fsc_text,0.90,6.5-2*i,"LC-"+str((i+1)*2),charsize=1.2
  
    endfor

  end_PS
 
end
