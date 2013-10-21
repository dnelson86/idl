; plotMaxTemps.pro
; gas accretion project - plots related to maximum past temperature of the gas
; dnelson jul.2013

; binValMaxHistos()

function binValMaxHistos, sP=sP, accMode=accMode, timeWindow=timeWindow

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  if sP.gfmWinds eq 1 and (accMode eq 'smooth' or accMode eq 'stripped' or accMode eq 'clumpy') then begin
    print,'Switching [' + accMode + '] to [' + accMode + '_rec] for GFM run.'
    accMode = accMode + '_rec'
  endif
  
  ; config (1D)
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
  if ~keyword_set(timeWindow) then message,'time window required (in Myr)'
  
  ; check if save exists
  saveFilename = sP.derivPath + 'binnedVals/binMaxVals.' + sP.saveTag + '.' + sP.savPrefix + str(sP.res) + '.' + str(sP.snap) + '.' + accMode + '_tw' + str(timeWindow) + '.sav'
  
  ; results exist, return
  if file_test(saveFilename) then begin
    restore,saveFilename
    return,r
  endif  

  ; load
  accTvir = gcSubsetProp(sP=sP,timeWindow=timeWindow,/accTvir,/accretionTimeSubset,accMode=accMode)
  maxTemp = gcSubsetProp(sP=sP,timeWindow=timeWindow,/maxPastTemp,/accretionTimeSubset,accMode=accMode)
  maxEnt  = gcSubsetProp(sP=sP,timeWindow=timeWindow,/maxPastEnt,/accretionTimeSubset,accMode=accMode)
  maxDens = gcSubsetProp(sP=sP,timeWindow=timeWindow,/maxPastDens,/accretionTimeSubset,accMode=accMode)
  maxMach = gcSubsetProp(sP=sP,timeWindow=timeWindow,/maxPastMach,/accretionTimeSubset,accMode=accMode)

  ; load parent halo masses so we can make halo massbins
  parentMass = gcSubsetProp(sP=sP,timeWindow=timeWindow,/parMass,/accretionTimeSubset,accMode=accMode)

  ; check: if no particles of a type are in w_tw, then replace their entries by a [0,0] array
  ; to allow histogram() to work, since they are now single numbers indexed by -1
  for i=0,n_tags(accTvir)-1 do begin
    if n_elements(accTvir.(i)) ne 1 then continue
    ;if accTvir.(i) ne -1 then message,'Strange, just one single valid entry? (sometimes with bhs in 128)'
    accTvir = mod_struct( accTvir, (tag_names(accTvir))[i], [100,100.1] )
    maxTemp = mod_struct( maxTemp, (tag_names(accTvir))[i], [0,0.1] )
    maxEnt  = mod_struct( maxEnt, (tag_names(accTvir))[i], [0,0.1] )
    maxDens = mod_struct( maxDens, (tag_names(accTvir))[i], [-1,-1.1] )
    maxMach = mod_struct( maxMach, (tag_names(accTvir))[i], [-1,-1.1] )
    parentMass = mod_struct( parentMass, (tag_names(accTvir))[i], [0,0.1] )
  endfor

  ; uniform weighting by mass for 2D histograms
  if sP.trMCPerCell gt 0 then massWt = sP.trMassConst * units.UnitMass_in_Msun
  if sP.trMCPerCell eq 0 then massWt = sP.targetGasMass * units.UnitMass_in_Msun
  if sP.trMCPerCell lt 0 then message,'error'  
  
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
  
  ; save
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sP.derivPath))
  
  return, r

end

; plotValMaxHistos(); plot (1) the previous max temp normalized by tviracc for arepo vs. gadget, gal vs. 
;                     gmem, (2) same but unnormalized by tviracc, (3) global not binned by halo mass but
;                     normalized by each parent tviracc, (4) same global without normalization

pro plotValMaxHistos, redshift=redshift, res=res

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  runs       = ['gadget','tracer','feedback'] ;['feedback','feedback_noZ','feedback_noFB']
  ;redshift   = 2.0
  ;res        = 256
  timeWindow = 500.0 ; Myr
  accMode    = 'smooth' ;'all','smooth','clumpy','stripped','recycled'
  tagNames   = ['Tmax','TmaxTviracc','EntMax','DensMax','MachMax'] ; plot each quantity
  
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
  yrangeG = [1e-4,0.3]
  pParams = { TmaxTviracc : {xrange:[-2.2,1.4], xtickv:[-2,-1,0,1],   ylabel : "log ( T_{max} / T_{vir,acc} )"} ,$
              Tmax        : {xrange:[3.0,8.0],  xtickv:[4,5,6,7],     ylabel : "log ( T_{max} )"}               ,$
              EntMax      : {xrange:[4.5,10.0], xtickv:[5,6,7,8,9],   ylabel : "log ( S_{max} ) [K cm^{2 }]"}   ,$
              DensMax     : {xrange:[-10,0],    xtickv:[-8,-6,-4,-2], ylabel : "log ( \rho_{max} )"}            ,$
              MachMax     : {xrange:[0,100],    xtickv:[10,30,50,80], ylabel : "M_{max}"}                        }     
              
  foreach tagName,tagNames do begin

    bVind   = ( where(tag_names(bV.(0).(0)) eq strupcase(tagName)) )[0]
    bVindG  = ( where(tag_names(bV.(0).(0)) eq 'GLOBAL'+strupcase(tagName)) )[0]
    bVind2d = ( where(tag_names(bV.(0).(0)) eq 'H2_'+strupcase(tagName)) )[0]
    pPind   = ( where(tag_names(pParams) eq strupcase(tagName)) )[0]
    if bVind eq -1 or pPind eq -1 or bvInd2d eq -1 or bVindG eq -1 then message,'Error'
    print,tagName,bVind,pPind,bVind2d,bVindG
    
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
          ytitle="",xtitle="",xticks=n_elements(xtickv)-1,xtickv=xtickv,$
          xtickname=xtickname,ytickname=ytickname,noerase=noerase       
        
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
      ;h2all = []
      ;for i=0,n_tags(sP)-1 do h2all = [h2all,bV.(i).allGal.(bVind2d)]
      ;crange = minmax(h2all[where(h2all ne 0.0)]) ;* [2.0,0.7]
      ;print,crange
      
      ; loop over each run
      for i=0,n_tags(sP)-1 do begin
      
        xtickv   = [9,10,11,12]
        xrange2D = bV.(i).params.mmMass2D
        yrange2D = bV.(i).params.mm.(bVind)
        xtitle   = textoidl("M_{halo} [_{ }log h^{-1} M_{sun }]")
        ytitle   = textoidl( pParams.(pPind).ylabel )
      
        ; plot 2d histo
        h2mt   = bV.(i).allGal.(bVind2d)
        crange = minmax(h2mt[where(h2mt ne 0.0)])
        
        loadColorTable, 'helix', /reverse
        tvim,h2mt^exp,scale=0,clip=-1,xtitle=xtitle,ytitle=ytitle,xrange=xrange2D,yrange=yrange2D,$
           xticks=n_elements(xtickv)-1,xtickv=xtickv,position=pos2D[i],noerase=(i gt 0),range=crange^exp
           
        ; simName
        legend,[sP.(i).simName],textcolor=[sP.(i).colors[cInd]],/bottom,/left
           
      endfor
           
      ; colorbar (NOTE: range is just representative of the last 2d histo, unless global crange)
      if n_elements(runs) eq 2 then begin
        barvals = findgen(ndivs+1)/ndivs*(max(h2mt^exp)-min(h2mt^exp)) + min(h2mt^exp)
        ticknames = textoidl(str(string(round(alog10(barvals^(1/exp))*10.0)/10.0,format='(f4.1)')))
        ticknames = ['0',ticknames[1:n_elements(ticknames)-1]]
      
        loadColorTable, 'helix', /reverse
        cgColorbar,bottom=1,range=minmax(h2mt),position=pos2D[2],/vertical,/right,title="",$
           divisions=ndivs,ticknames=ticknames,ncolors=255
           
        cgText,textoidl("M_{tot} [_{ }log h^{-1} M_{sun }]"),0.98,(0.18+0.94)/2,alignment=0.5,orientation=90,/normal
      endif
      
      ; plot (3) - global
      start_PS, sP.(0).plotPath + tagName + '_global_allGal_gmem.'+plotStr+'.eps', /big
        !p.thick += 1
        
        xrange = pParams.(pPind).xrange
        xtickv = pParams.(pPind).xtickv
        ylabel = pParams.(pPind).ylabel
        
        cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeG,/xs,/ys,/ylog,yminor=0,$
          ytitle="Gas Mass Fraction",xtitle=textoidl(ylabel),xticks=n_elements(xtickv)-1,xtickv=xtickv
        
        cgPlot,[0,0],[2e-4,0.14],line=2,color=cgColor('black'),thick=!p.thick-0.0,/overplot
        
        for i=0,n_tags(sP)-1 do begin
          ; gal
          hist = bV.(i).allGal.(bVindG)[*]
          cgPlot,bV.(i).params.binLoc.(bvInd),float(hist)/total(hist),line=lines[0],color=sP.(i).colors[cInd],/overplot

          ; gmem
          hist = bV.(i).gmem.(bVindG)[*]
          cgPlot,bV.(i).params.binLoc.(bvInd),float(hist)/total(hist),line=lines[1],color=sP.(i).colors[cInd],/overplot
        endfor
        
        ; legend
        legend,simNames,textcolors=simColors,box=0,/top,/left
        legend,['galaxy','halo'],linestyle=[0,1],color=cgColor('dark gray'),$
          textcolors=['dark gray','dark gray'],linesize=0.25,box=0,/top,/right
      end_PS
                
    end_PS
        
  endforeach ; tagNames
  
end

; plotValMaxHistosByMode(); as above, but separated into modes (line plots only)

pro plotValMaxHistosByMode, redshift=redshift, res=res

  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  runs       = ['tracer']; ,'gadget'];,'feedback'] ;['feedback','feedback_noZ','feedback_noFB']
  ;redshift   = 2.0
  ;res        = 256
  timeWindow = 500.0 ; Myr
  accModes   = ['all','smooth','clumpy','stripped','recycled']
  tagNames   = ['Tmax','TmaxTviracc','EntMax','DensMax','MachMax'] ; plot each quantity
  
  ; plot config
  lines   = [0,1,2,3,4] ; for each accMode
  cInd    = 1 ; color index for simName labels
  
  if isnumeric(timeWindow) then twStr = str(fix(timeWindow))
  if ~isnumeric(timeWindow) then twStr = "-"+timeWindow
     
  ; load
  foreach run,runs,i do begin
    sP  = mod_struct(sP, 'sP'+str(i), simParams(res=res,run=run,redshift=redshift))
    bV_mode = {}
    
    ; make for all the accretion modes
    foreach accMode,accModes,j do begin
      if accMode eq 'recycled' and sP.(i).gfmWinds eq 0 then continue ; skip recycled for nonWinds     
      bV_mode = mod_struct(bV_mode, 'mode'+str(j), $
        binValMaxHistos(sP=sP.(i),accMode=accMode,timeWindow=timeWindow))
    endforeach    
    
    ; put this mode collection into mbv, once per run
    bV = mod_struct(bV, 'bV'+str(i), bV_mode)
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

  plotStr += str(res) + '_' + str(sP.(0).snap) + '_tw' + twStr
    
  ; plot config
  pos     = plot_pos(rows=2, cols=3)
  yrange  = [6e-4,2.0]
  yrangeG = [1e-4,0.3]
  pParams = { TmaxTviracc : {xrange:[-2.2,1.4], xtickv:[-2,-1,0,1],   ylabel : "log ( T_{max} / T_{vir,acc} )"} ,$
              Tmax        : {xrange:[3.0,8.0],  xtickv:[4,5,6,7],     ylabel : "log ( T_{max} )"}               ,$
              EntMax      : {xrange:[4.5,10.0], xtickv:[5,6,7,8,9],   ylabel : "log ( S_{max} ) [K cm^{2 }]"}   ,$
              DensMax     : {xrange:[-10,0],    xtickv:[-8,-6,-4,-2], ylabel : "log ( \rho_{max} )"}            ,$
              MachMax     : {xrange:[0,100],    xtickv:[10,30,50,80], ylabel : "M_{max}"}                        }     
              
  foreach tagName,tagNames do begin

    bVind   = ( where(tag_names(bV.(0).(0).(0)) eq strupcase(tagName)) )[0]
    bVindG  = ( where(tag_names(bV.(0).(0).(0)) eq 'GLOBAL'+strupcase(tagName)) )[0]
    pPind   = ( where(tag_names(pParams) eq strupcase(tagName)) )[0]
    if bVind eq -1 or pPind eq -1 or bVindG eq -1 then message,'Error'
    print,tagName,bVind,pPind,bVindG
    
    ; plot (1) - 3x2 mass bins separated out and each panel with all runs, allGal
    start_PS, sP.(0).plotPath + tagName + '_byMode_3x2_allGal.'+plotStr+'.eps', /big
      !p.thick += 1
      
      xrange = pParams.(pPind).xrange
      xtickv = pParams.(pPind).xtickv
      ylabel = pParams.(pPind).ylabel
      
      for j=0,n_elements(bV.(0).(0).params.massBins)-2 do begin
        
        if j eq 0 or j eq 3 then ytickname = '' else ytickname = replicate(' ',10)
        if j ge 3 then xtickname = '' else xtickname = replicate(' ',10)
        if j gt 0 then noerase = 1 else noerase = 0
        
        cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,/ylog,yminor=0,pos=pos[j],$
          ytitle="",xtitle="",xticks=n_elements(xtickv)-1,xtickv=xtickv,$
          xtickname=xtickname,ytickname=ytickname,noerase=noerase       
        
        cgPlot,[0,0],[8e-4,0.25],line=2,color=cgColor('black'),thick=!p.thick-0.0,/overplot
        
        for i=0,n_tags(sP)-1 do begin
          for k=0,n_tags(bV.(i))-1 do begin
            ; gal or gmem
            hist = bV.(i).(k).allGal.(bVind)[j,*]
            cgPlot,bV.(i).(k).params.binLoc.(bvInd),float(hist)/total(hist),line=lines[k],color=sP.(i).colors[cInd],/overplot

            ; gmem
            ;hist = bV.(i).(k).gmem.(bVind)[j,*]
            ;cgPlot,bV.(i).(k).params.binLoc.(bvInd),float(hist)/total(hist),line=lines[k],color=sP.(i).colors[cInd],/overplot
          endfor
        endfor
        
        ; legends
        massBinStr = string(bV.(0).(0).params.massBins[j],format='(f4.1)') + ' < log(M) < ' + $
                     string(bV.(0).(0).params.massBins[j+1],format='(f4.1)')
  
        cgText,mean(xrange),yrange[1]*0.4,massBinStr,charsize=!p.charsize-0.0,alignment=0.5
            
        if j eq 0 then legend,simNames,textcolors=simColors,box=0,position=[xrange[0],0.5]
        if j eq 3 then legend,accModes,linestyle=lines,color=cgColor('dark gray'),$
          textcolors=replicate('dark gray',n_elements(accModes)),linesize=0.25,box=0,position=[xrange[0],0.5]
      endfor
      
      ; axis labels
      cgText,0.05,mean([ (pos[0])[3], (pos[3])[1] ]),"Gas Mass Fraction",alignment=0.5,orientation=90.0,/normal
      cgText,mean([ (pos[0])[0], (pos[2])[3] ]),0.05,textoidl(ylabel),alignment=0.5,/normal
      
    end_PS
    
    ; plot (2) - global
    start_PS, sP.(0).plotPath + tagName + '_byMode_global_allGal_gmem.'+plotStr+'.eps', /big
      !p.thick += 1
        
      xrange = pParams.(pPind).xrange
      xtickv = pParams.(pPind).xtickv
      ylabel = pParams.(pPind).ylabel
        
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeG,/xs,/ys,/ylog,yminor=0,$
        ytitle="Gas Mass Fraction",xtitle=textoidl(ylabel),xticks=n_elements(xtickv)-1,xtickv=xtickv
        
      cgPlot,[0,0],[2e-4,0.14],line=2,color=cgColor('black'),thick=!p.thick-0.0,/overplot
      
      for i=0,n_tags(sP)-1 do begin
        for k=0,n_tags(bV.(i))-1 do begin
          ; gal
          hist = bV.(i).(k).allGal.(bVindG)[*]
          cgPlot,bV.(i).(k).params.binLoc.(bvInd),float(hist)/total(hist),line=lines[k],color=sP.(i).colors[cInd],/overplot
 
          ; gmem
          ;hist = bV.(i).(k).gmem.(bVindG)[*]
          ;cgPlot,bV.(i).(k).params.binLoc.(bvInd),float(hist)/total(hist),line=lines[1],color=sP.(i).colors[cInd],/overplot
        endfor
      endfor
      
      ; legend
      legend,simNames,textcolors=simColors,/top,/left
      legend,accModes,linestyle=lines,color=cgColor('dark gray'),$
        textcolors=replicate('dark gray',n_elements(accModes)),linesize=0.35,position=[0.80,0.92],/normal
    end_PS
      
  endforeach ; tagNames
  
end
