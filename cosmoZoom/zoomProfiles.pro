; zoomProfiles.pro
; 'zoom project' 1D and 2D radial profiles, phase diagrams, PDFs of various quantities
; dnelson nov.2014

; zoomRadialBin(): do radial binning (1d and 2d) for multiple particle types

function zoomRadialBin, sP=sP, gc=gc, partType=partType, halo=halo, mode=mode
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits(redshift=sP.redshift)
  
  ; mode
  if n_elements(mode) eq 0 then mode = 0
  
  ; setup binning
  binSize = {}
  binCen  = {}
  
  for i=0,n_tags(halo.binMinMax)-1 do begin
    bsLoc = ( halo.binMinMax.(i)[1] - halo.binMinMax.(i)[0] ) / halo.nBins
    binSize = mod_struct( binSize, (tag_names(halo.binMinMax))[i], bsLoc )
    binCen  = mod_struct( binCen, (tag_names(halo.binMinMax))[i], fltarr(halo.nBins) )
  endfor
  
  saveFilename = sP.derivPath + 'cutouts/zoomRadialBin.gc'+str(halo.gcInd)+'.'+sP.savPrefix+str(sP.res) + $
                 '.cutSubS-' + str(halo.cutSubS) + '_'+partType+'_r'+str(fix(halo.binMinMax.rad[1]*100))+'.sav'

  if ~file_test(saveFilename) then begin
    ; load positions of gas, make radial cut for halo selection
    pos = loadSnapshotSubset(sP=sP,partType=partType,field='pos')
      
    rad_pri  = periodicDists(halo.pPos,pos,sP=sP)
    rad_pri /= halo.rVir
    
    wRad = where(rad_pri le halo.binMinMax.rad[1]*1.5,countRad)
    print,' ['+partType+'] within '+string(halo.binMinMax.rad[1]*1.5,format='(f3.1)')+' rvir, have ['+$
      str(countRad)+'] of ['+str(n_elements(rad_pri))+']'
          
    ; remove substructures? load ids and make secondary list
    if halo.cutSubS eq 1 then begin
      ; verify our target ID is a primary (otherwise...)
      priGCIDs = gcIDList(gc=gc,select='pri')
      if total(halo.gcInd eq priGCIDs) eq 0 then message,'Error'
      
      ids = loadSnapshotSubset(sP=sP,partType=partType,field='ids',inds=wRad)
      
      ; get ids of all gas cells in all substructures other than the target
      ; (this differs from ids in the target, as it includes cells not in any subfind group)
      myPIDs  = gcPIDList(gc=gc,valGCids=[halo.gcInd],partType=partType)
      subPIDs = gcPIDList(gc=gc,select='all',partType=partType)
      remPIDs = removeIntersectionFromB(myPIDs,subPIDs)
      if n_elements(remPIDs) ne n_elements(subPIDs)-n_elements(myPIDs) then message,'Error'
      
      ; remove the intersection of (satPIDs,ids) from wRad
      calcMatch,remPIDs,ids,ind1,ind2,count=count
      
      all = bytarr(countRad)
      if count gt 0 then all[ind2] = 1B ; flag substructure gas cells, then select against
      
      ; update wRad and countRad
      wRad = wRad[ where(all eq 0B, countRad) ]
      
      print,'  substructures ['+str(count)+'] removed! now have ['+str(countRad)+'] in radius.'
    endif
      
    pos     = pos[*,wRad]
    rad_pri = rad_pri[wRad]
      
    ; load further properties, handle little h and unit conversions
    scalefac = 1.0 / (1.0 + sP.redshift)
    
    if partType eq 'gas' then begin
      u     = loadSnapshotSubset(sP=sP,partType='gas',field='u',inds=wRad)
      nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec',inds=wRad)
      dens  = loadSnapshotSubset(sP=sP,partType='gas',field='density',inds=wRad)
      nh0   = loadSnapshotSubset(sP=sP,partType='gas',field='neutralhydrogenabundance',inds=wRad)
      mass  = loadSnapshotSubset(sP=sP,partType='gas',field='mass',inds=wRad)
      csize = loadSnapshotSubset(sP=sP,partType='gas',field='cellsize',inds=wRad)
      sfr   = loadSnapshotSubset(sP=sP,partType='gas',field='sfr',inds=wRad)

      temp  = convertUtoTemp(u, nelec)
      
      ; modifiers for star forming gas
      w = where(sfr gt 0.0,countSfr) 
      if countSfr gt 0 then temp[w] = 3000.0  ; set eEOS temperature down to minimum of h0L10@z=2 (3000K)
      if countSfr gt 0 then nh0[w] = 1.0 ; set neutral fraction=1.0 for eEOS (avoid effective temp)
      
      mass  = mass * 1e10 / units.hubbleParam ; 10^10 msun/h -> msun
      entr  = calcEntropyCGS(u, dens, sP=sP)
      csize = csize * (scalefac / units.HubbleParam) ; ckpc/h -> physical
      nh0   = nH0ToPhys(nh0, dens, sP=sP, /cgs) ; cgs mass density of H0 (divide by m_p for number dens)
    endif
    
    if partType eq 'dm' then begin
      h = loadSnapshotHeader(sP=sP)
      dmPartMass = h.massTable[ partTypeNum('dm') ]
      
      ; constant mass based tophat density estimator
      dens = estimateDensityTophat(pos, mass=dmPartMass, ndims=3, nNGB=32, boxSize=sP.boxSize)
    endif
    
    if partType eq 'stars' then begin
      thMode = 3 ; 3=total/volume (density)
      mass   = loadSnapshotSubset(sP=sP,partType='stars',field='mass',inds=wRad)
      mass   = reform(mass, [1,n_elements(mass)])
      posVal = [pos,mass] 
      
      ; quantity estimator in tophat neighbor search
      dens = calcTHVal(posVal, pos, thMode=thMode, ndims=3, nNGB=32, boxSize=sP.boxSize)
    endif
    
    ; density unit conversion (from code_mass/code_volume)
    dens = codeDensToPhys(dens, sP=sP, /cgs)
    dens /= units.mass_proton ;g/cm^3 -> cm^(-3) number density
    
    ; vrad calculation: replace coordinates by relative coordinates (radial vectors) to direct parent
    for i=0,2 do begin
      pos_rel = reform(pos[i,*] - halo.pPos[i])
      correctPeriodicDistVecs, pos_rel, sP=sP
      pos[i,*] = pos_rel
    endfor
    
    vel = loadSnapshotSubset(sP=sP,partType=partType,field='vel',inds=wRad)
    vel *= sqrt(scalefac)

    ; subgroupVel already peculiar (no scalefac correction needed)
    vrad = ((vel[0,*] - halo.pVel[0]) * pos[0,*] + $ 
            (vel[1,*] - halo.pVel[1]) * pos[1,*] + $
            (vel[2,*] - halo.pVel[2]) * pos[2,*]) $
                / (rad_pri*halo.rVir)
            
    vrad += (rad_pri*halo.rVir) * units.H_z ; add in hubble expansion (neglect mass growth term)
            
    ; angm calculation: mean magnitude of specific angular momentum = rvec x vel
    ; make velocities relative to bulk halo motion
    vel[0,*] = reform(vel[0,*] - halo.pVel[0])
    vel[1,*] = reform(vel[1,*] - halo.pVel[1])
    vel[2,*] = reform(vel[2,*] - halo.pVel[2])
    
    ; angular momentum magnitude
    jvec = fltarr(3,countRad)
    jvec[0,*] = pos[1,*] * vel[2,*] - pos[2,*] * vel[1,*]
    jvec[1,*] = pos[2,*] * vel[0,*] - pos[0,*] * vel[2,*]
    jvec[2,*] = pos[0,*] * vel[1,*] - pos[1,*] * vel[0,*]
    jnorm = reform( sqrt(jvec[0,*]*jvec[0,*] + jvec[1,*]*jvec[1,*] + jvec[2,*]*jvec[2,*]) )
          
    if partType eq 'gas' then $
      save,rad_pri,dens,vrad,jnorm,temp,mass,entr,csize,nh0,filename=saveFilename
    if partType eq 'dm' then $
      save,rad_pri,dens,vrad,jnorm,filename=saveFilename
    if partType eq 'stars' then $
      save,rad_pri,dens,vrad,jnorm,mass,filename=saveFilename
      
    print,'Saved file: ['+strmid(saveFilename,strlen(sP.derivPath))+']'
  endif else begin
    print,'RESTORE: ['+strmid(saveFilename,strlen(sP.derivPath))+']'
    restore,saveFilename
  endelse
  
  ; final mass modifications
  if partType eq 'stars' then begin
    mass *= 1e10 / units.hubbleParam
  endif
  if partType eq 'dm' then begin
    h = loadSnapshotHeader(sP=sP)
    dmPartMass = h.massTable[ partTypeNum('dm') ]
    mass = replicate( dmPartMass * 1e10 / units.hubbleParam, n_elements(rad_pri) )
  endif
  
  ; stuff into struct for automated binning vs. all parameters
  if partType eq 'gas' then $
    data = {rad:rad_pri,dens:dens,angm:jnorm,vrad:vrad,temp:temp,csize:csize,entr:entr,nh0:nh0}
  if partType eq 'dm' then $
    data = {rad:rad_pri,dens:dens,angm:jnorm,vrad:vrad}
  if partType eq 'stars' then $
    data = {rad:rad_pri,dens:dens,angm:jnorm,vrad:vrad}
         
  gasOnlyFields = ['TEMP','NH0','CSIZE','ENTR']
         
  ; 1D and 2D histograms (all gas with r/rvir<2) binning
  ; ----------------------------------------------------
  if mode eq 0 then begin
    
    foreach binFieldName,tag_names(halo.binMinMax) do begin
      ; make struct (4 per quantity = [mean,lower quartile,median,upper quartile])
      r = {}
      
      foreach fName,tag_names(halo.binMinMax) do begin
        if total(fName eq gasOnlyFields gt 0) and partType ne 'gas' then continue
        r = mod_struct( r, fName, fltarr(halo.nBins,4) )
        r = mod_struct( r, fName+'_2D', fltarr(halo.nBins,halo.nBins) )
      endforeach
      
      ; indices
      if total(binFieldName eq gasOnlyFields gt 0) and partType ne 'gas' then continue
      tD1 = where( strcmp(tag_names(data),binFieldName) eq 1 )
      tR1 = where( strcmp(tag_names(r),binFieldName) eq 1 ) ; unused
      tH1 = where( strcmp(tag_names(halo.binMinMax),binFieldName) eq 1 ) ; also binSize,binCen 
      
      print,' '+binFieldName
      
      for i=0,halo.nBins-1 do begin
        ; first dimension: boundaries
        locMin = 10.0^( halo.binMinMax.(tH1)[0] + binSize.(tH1)*i )
        locMax = 10.0^( halo.binMinMax.(tH1)[0] + binSize.(tH1)*(i+1) )
        if binFieldName eq 'RAD' or binFieldName eq 'VRAD' then begin ; only linear ones
          locMin = halo.binMinMax.(tH1)[0] + binSize.(tH1)*i
          locMax = halo.binMinMax.(tH1)[0] + binSize.(tH1)*(i+1)
        endif
        
        w = where(data.(tD1) ge locMin and data.(tD1) lt locMax,count)
        
        binCen.(tH1)[i] = mean([locMin,locMax]) ; unused
        
        if count eq 0 then continue
        
        foreach fieldName,tag_names(halo.binMinMax) do begin
          if total(fieldName eq gasOnlyFields gt 0) and partType ne 'gas' then continue

          ; indices
          tD2   = where( strcmp(tag_names(data),fieldName) eq 1 )
          tR_1D = where( strcmp(tag_names(r),fieldName) eq 1 )
          tR_2D = where( strcmp(tag_names(r),fieldName+'_2D') eq 1 )
          tH2   = where( strcmp(tag_names(halo.binMinMax),fieldName) eq 1 ) ; also binSize,binCen 
          
          ; 1D binning
          r.(tR_1D)[i,0]   = mean( data.(tD2)[w] )
          r.(tR_1D)[i,1:3] = percentiles( data.(tD2)[w] )
          
          ; 2D binning
          for j=0,halo.nBins-1 do begin
            ; second dimension: boundaries
            locMin = 10.0^( halo.binMinMax.(tH2)[0] + binSize.(tH2)*j )
            locMax = 10.0^( halo.binMinMax.(tH2)[0] + binSize.(tH2)*(j+1) )
            if fieldName eq 'RAD' or fieldName eq 'VRAD' then begin ; only linear ones
              locMin = halo.binMinMax.(tH2)[0] + binSize.(tH2)*j
              locMax = halo.binMinMax.(tH2)[0] + binSize.(tH2)*(j+1)
            endif
            
            wr = where(data.(tD2)[w] ge locMin and data.(tD2)[w] lt locMax,count_r)

            ;binCen.(tH2)[j] = mean([locMin,locMax]) ; unused
            ;if i eq 92 and j eq 92 and binFieldName eq 'TEMP' and fieldName eq 'RAD' then stop
            if count_r eq 0 then continue
            
            r.(tR_2D)[i,j] = total( mass[w[wr]] )
          endfor ;nBins,j
        endforeach ;fieldName
      endfor ;nBins,i
      
      r = mod_struct( r, 'binCen', binCen )
      r = mod_struct( r, 'binSize', binSize )
      
      rr = mod_struct( rr, binFieldName, r )
    endforeach ; binFieldNames
    
    ;if partType eq 'gas' then begin ; DEBUG
    ;  rr = mod_struct( rr, 'data', data )
    ;endif ; DEBUG END
  endif ; mode=0
  
  ; 1d histograms (PDFs) in restricted radial ranges
  ; ------------------------------------------------
  if mode eq 1 then begin
    foreach binFieldName,tag_names(halo.binMinMax) do begin
      ; make struct (4 per quantity = [mean,lower quartile,median,upper quartile])
      r = {}
      
      foreach fName,tag_names(halo.binMinMax) do begin
        if total(fName eq gasOnlyFields gt 0) and partType ne 'gas' then continue
        r = mod_struct( r, fName, fltarr(halo.nBins,4) )
      endforeach
      
      ; indices
      if total(binFieldName eq gasOnlyFields gt 0) and partType ne 'gas' then continue
      tD1 = where( strcmp(tag_names(data),binFieldName) eq 1 )
      tR1 = where( strcmp(tag_names(r),binFieldName) eq 1 ) ; unused
      tH1 = where( strcmp(tag_names(halo.binMinMax),binFieldName) eq 1 ) ; also binSize,binCen 
      
      print,' '+binFieldName
      
      for i=0,halo.nBins-1 do begin
        ; first dimension: boundaries
        locMin = 10.0^( halo.binMinMax.(tH1)[0] + binSize.(tH1)*i )
        locMax = 10.0^( halo.binMinMax.(tH1)[0] + binSize.(tH1)*(i+1) )
        if binFieldName eq 'RAD' or binFieldName eq 'VRAD' then begin ; only linear ones
          locMin = halo.binMinMax.(tH1)[0] + binSize.(tH1)*i
          locMax = halo.binMinMax.(tH1)[0] + binSize.(tH1)*(i+1)
        endif
        
        w = where(data.(tD1) ge locMin and data.(tD1) lt locMax,count)
        
        binCen.(tH1)[i] = mean([locMin,locMax]) ; unused
        
        if count eq 0 then continue
        
        foreach fieldName,tag_names(halo.binMinMax) do begin
          if total(fieldName eq gasOnlyFields gt 0) and partType ne 'gas' then continue

          ; indices
          tD2   = where( strcmp(tag_names(data),fieldName) eq 1 )
          tR_1D = where( strcmp(tag_names(r),fieldName) eq 1 )
          
          ; 1D binning
          r.(tR_1D)[i,0]   = mean( data.(tD2)[w] )
          r.(tR_1D)[i,1:3] = percentiles( data.(tD2)[w] )
          
        endforeach ;fieldName
      endfor ;nBins,i
      
      r = mod_struct( r, 'binCen', binCen )
      r = mod_struct( r, 'binSize', binSize )
      
      rr = mod_struct( rr, binFieldName, r )
    endforeach ; binFieldNames
  endif
  
  return, rr
        
end

; plotZoomRadialProfiles() - 1d with all halos on same plot, and 2d with each halo in separate panel

pro plotZoomRadialProfiles ;, input=input
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  hInds     = [0,1,7,8] ;[0,1,2,3,4,5,6,7,8,9] [0,1,7,8]
  resLevels = [9,10,11] ;[9,10,11]
  redshift  = 2.0
  newSaves  = 0 ; override existing saves
  
  ; binning config
  nBaseBins = 100 ; for base resolution level, increased at higher resolutions
  cutSubS   = 1   ; remove substructures?
  
  ; 2D plot config (binMinMax also plot bounds for 2D histos, must rebin if you change)
  binMinMax = { rad   : [0.0,2.0]   ,$
                temp  : [4.0,7.0]   ,$
                dens  : [-6.0,0.0]  ,$
                vrad  : [-400,250]  ,$
                csize : [-1.0,0.75] ,$
                entr  : [5.5,9.5]   ,$
                angm  : [3.0,5.0]   ,$
                nh0   : [-28,-24]    }
                
  ctName2D      = 'ncl/WhiteBlueGreenYellowRed' ;-dnA ; 2d histogram
  lineColor2D   = ['black','black'] ; mean/median line color
  line2D        = [3,0]             ; mean/median line style
  thick2D       = [4.0,4.0]         ; mean/median line thickness
  
  h2dMinMax   = [0,8e7] ; linear color scaling of total mass per pixel (comment out to scale each 
                        ; 2D histogram independently from its min to max)
  logHist     = 0       ; take log of 2d histogram pixel values?
  byRow       = 0       ; independently normalize 2d histograms row by row?
  byCol       = 1       ; independently normalize 2d histograms col by col?

  matrixPlots = ['RAD','TEMP','DENS','VRAD','ENTR','ANGM'] ; correlation matrix quantities
  indivPlotsX = ['RAD']   ; single-panel 2d histos?
  indivPlotsY = ['CSIZE'] ; single-panel 2d histos?
  
  ; plot config
  pConfig = { rad   : { label:"r / r_{vir}",                 range:[0.0,2.0],  log:0 } ,$
              temp  : { label:"log T_{gas} [_{ }K_{ }]",     range:[4.5,6.5],  log:1 } ,$
              dens  : { label:"log n_{gas} [_{ }cm^{-3 }]",  range:[-4,-1],    log:1 } ,$
              vrad  : { label:"v_{rad} [_{ }km/s_{ }]",      range:[-200,100], log:0 } ,$
              csize : { label:"log r_{cell} [_{ }kpc_{ }]",  range:[-1.0,1.0], log:1 } ,$
              entr  : { label:"log S_{gas} [_{ }K cm^{2 }]", range:[7.0,9.0],  log:1 } ,$
              angm  : { label:"j_{gas} [_{ }kpc km/s_{ }]",  range:[3.0,5.0],  log:1 }  }
              
  pConfig_dm = { rad   : { label:"r / r_{vir}",                range:[0.0,2.0],  log:0 } ,$
                 dens  : { label:"log n_{DM} [_{ }cm^{-3 }]",  range:[-6.0,1.0], log:1 } ,$
                 vrad  : { label:"v_{rad,DM} [_{ }km/s_{ }]",  range:[-100,200], log:0 } ,$
                 angm  : { label:"j_{DM} [_{ }kpc km/s_{ }]",  range:[3.0,5.0],  log:1 }  }
                 
  pConfig_stars = { rad   : { label:"r / r_{vir}",                   range:[0.0,2.0],  log:0 } ,$
                    dens  : { label:"log n_{stars} [_{ }cm^{-3 }]",  range:[-6.0,1.0], log:1 } ,$
                    vrad  : { label:"v_{rad,stars} [_{ }km/s_{ }]",  range:[-100,200], log:0 } ,$
                    angm  : { label:"j_{stars} [_{ }kpc km/s_{ }]",  range:[3.0,5.0],  log:1 }  }
         
  radLines = [0.15,1.5]
  sK       = 1
  
  ; A (resolutions are different linestyles, all same color and thickness)
  ;lines    = [1,2,0] ; line style, one per resLevel
  ;cInds    = [1,1,1] ; color index, one per resLevel
  ;thick    = [1,1,1] ; times !p.thick

  ; B (resolutions are all solid lines, fainter and thinner for lower res)
  lines   = [0,0,0]
  thick   = [0.25,0.5,1.0]
  cInds   = [2,1,1]
  
  ; load  
  foreach hInd,hInds do begin
    hLoc = {}
        
    foreach resLevel,resLevels do begin
      ; get list of snapshots for this halo in the redshift range
      rLoc = {}
      
      if resLevel eq 11 and total(hInd eq [2,3,4,5,6,9]) gt 0 then continue
      
      sP = simParams(run='zoom_20Mpc',res=resLevel,hInd=hInd,redshift=redshift)     
      saveFilename = sP.derivPath + 'binnedVals/radProfiles.h'+str(hInd)+'.'+sP.savPrefix+str(sP.res) + $
                     '.cutSubS-' + str(cutSubS) + '.sav'
      
      if ~file_test(saveFilename) or newSaves eq 1 then begin
                  
        ; adjust binning based on resolution (increase as cube root of number of particles)
        nBins = nBaseBins * 2.0^(sP.zoomLevel-2)
        
        ; load group catalog, locate halo
        gc = loadGroupCat(sP=sP,readIDs=(cutSubS eq 1),skipIDs=(cutSubS eq 0))
        gcInd = zoomTargetHalo(sP=sP)
        
        rVir = gc.group_r_crit200[gc.subgroupGrNr[gcInd]]
        tVir = codeMassToVirTemp( gc.subgroupMass[gcInd], redshift=sP.redshift )
        pPos = gc.subgroupPos[*,gcInd]
        pVel = gc.subgroupVel[*,gcInd]
        
        print,' h'+str(hInd)+'L'+str(resLevel)+': ['+string(sP.snap,format='(I2)')+'] gcInd = ' + $
              str(gcInd) + ' rVir = ' + string(rVir,format='(f5.1)') + $
              ' tVir = ' + string(alog10(tVir),format='(f3.1)')
        
        ; radial binning, 1D and 2D, all quantity combinations
        halo = { tVir:tVir, rVir:rVir, pPos:pPos, pVel:pVel, gcInd:gcInd, $
                 nBins:nBins, binMinMax:binMinMax, cutSubS:cutSubS }
        
        gas   = zoomRadialBin(sP=sP, gc=gc, partType='gas', halo=halo)
        dm    = zoomRadialBin(sP=sP, gc=gc, partType='dm', halo=halo)
        stars = zoomRadialBin(sP=sP, gc=gc, partType='stars', halo=halo)
        
        ; save all values for this halo+resLevel combination
        rLoc = mod_struct( rLoc, 'sP', sP )
        rLoc = mod_struct( rLoc, 'gas', gas )
        rLoc = mod_struct( rLoc, 'dm',  dm )
        rLoc = mod_struct( rLoc, 'stars', stars )
        rLoc = mod_struct( rLoc, 'halo', halo )
       
        save,rLoc,filename=saveFilename
        print,'Saved file: ['+strmid(saveFilename,strlen(sP.derivPath))+']'
      endif else begin
        print,'RESTORE: ['+strmid(saveFilename,strlen(sP.derivPath))+']'
        restore,saveFilename
      endelse
    
      ; add this resolution level to the halo struct
      hLoc = mod_struct( hLoc, 'L'+str(resLevel), rLoc )     
     
    endforeach ;resLevels
    
    ; add this halo to the save struct
    rp = mod_struct( rp, 'h'+str(hInd), hLoc )
    
  endforeach ;hInds

  ; plot setup
  plotStr = 'h'
  hNames  = 'h' + str(hInds)
  hColors = []
  rNames  = 'L' + str(resLevels)
  rLines  = lines[0:n_elements(resLevels)-1]
  rThick  = thick[0:n_elements(resLevels)-1]

  foreach hInd,hInds do plotStr += str(hInd)
  foreach hInd,hInds,i do hColors = [hColors,rp.(i).(0).sP.colors[cInds[-1]]]
  plotStr += '_L'
  foreach resLevel,resLevels do plotStr += str(resLevel)
  plotStr += '_cutSubS-' + str(cutSubS)
  
  !except = 0
  
  ; plot (1) - 1D radial profiles of all quantities, one per panel, all halos all resolutions
  foreach fName,tag_names(pConfig) do begin
  
    qInd = where( tag_names(rp.(0).(0).gas) eq fName )
    bInd = where( tag_names(rp.(0).(0).gas.(qInd).binCen) eq fName )
    
    pInd = where( tag_names(pConfig) eq fName )
    xtitle = pConfig.(pInd).label
    xrange = pConfig.(pInd).range
    
    start_PS, rp.(0).(0).sP.plotPath + 'zoomProfiles_'+fName+'_Gas_1D_' + plotStr + '.eps', $
      xs=12.0*1.4, ys=6.0*1.4
    
      pos = plot_pos(col=3,row=2,/gap)
      offset = [-0.03,-0.02,-0.02,0]
      
      count = 0
      
      foreach curField,tag_names(pConfig) do begin
        if curField eq fName then continue
        
        yInd = where( tag_names(pConfig) eq curField )
        zInd = where( tag_names(rp.(0).(0).gas.(qInd)) eq curField )
        yrange = pConfig.(yInd).range
        ytitle = pConfig.(yInd).label
        
        ; plot
        cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xlog=0,$
          xtitle=textoidl(xtitle),ytitle=textoidl(ytitle),pos=pos[count]+offset,/noerase
          
        foreach radLine,radLines do $
          cgPlot,[radLine,radLine],yrange+[0.05,-0.05],line=1,color='light gray',/overplot
          
        for i=0,n_tags(rp)-1 do begin
          for j=0,n_tags(rp.(i))-1 do begin
            xx = rp.(i).(j).gas.(qInd).binCen.(bInd)
            if pConfig.(pInd).log eq 1 then xx = alog10( xx )
            yy = rp.(i).(j).gas.(qInd).(zInd)
            if pConfig.(yInd).log eq 1 then yy = alog10( yy )
            
            co = rp.(i).(j).sP.colors[cInds[j]]
            th = !p.thick * thick[j]
                    
            cgPlot,xx,smooth(yy[*,0],sK),line=lines[j],color=co,thick=th,/overplot
          endfor ;n_tags(rp.(i)),j
        endfor ;n_tags(rp),i
        
        if count eq 0 then legend,hNames,textcolor=hColors,/top,/right
        if count eq 1 then legend,rNames,linestyle=rLines,thick=rThick*!p.thick,/top,/right
        
        count += 1
      endforeach ; tag_names(pConfig)
        
    end_PS
  
  endforeach ; tag_names(pConfig)
  
  stop
  
  ; plot (2) - 1D radial profiles of DM/stars quantities, one per panel, all halos all resolutions
  foreach fName,tag_names(pConfig_dm) do begin
  
    start_PS, rp.(0).(0).sP.plotPath + 'zoomProfiles_'+fName+'_DMStars_1D_' + plotStr + '.eps', $
      xs=12.0*1.4, ys=6.0*1.4
    
      pos = plot_pos(col=3,row=2,/gap)
      offset = [-0.03,-0.02,-0.02,0]
      
      count = 0
      
      ; (abc) - upper three panels are DM
      qInd = where( tag_names(rp.(0).(0).dm) eq fName )
      bInd = where( tag_names(rp.(0).(0).dm.(qInd).binCen) eq fName )
      
      pInd = where( tag_names(pConfig_dm) eq fName )
      xtitle = pConfig_dm.(pInd).label
      xrange = pConfig_dm.(pInd).range
    
      foreach curField,tag_names(pConfig_dm) do begin
        if curField eq fName then continue
        
        yInd = where( tag_names(pConfig_dm) eq curField )
        zInd = where( tag_names(rp.(0).(0).dm.(qInd)) eq curField )
        yrange = pConfig_dm.(yInd).range
        ytitle = pConfig_dm.(yInd).label
        
        ; plot
        cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xlog=0,$
          xtitle=textoidl(xtitle),ytitle=textoidl(ytitle),pos=pos[count]+offset,/noerase
          
        foreach radLine,radLines do $
          cgPlot,[radLine,radLine],yrange+[0.05,-0.05],line=1,color='light gray',/overplot
          
        for i=0,n_tags(rp)-1 do begin
          for j=0,n_tags(rp.(i))-1 do begin
            xx = rp.(i).(j).dm.(qInd).binCen.(bInd)
            if pConfig_dm.(pInd).log eq 1 then xx = alog10( xx )
            yy = rp.(i).(j).dm.(qInd).(zInd)
            if pConfig_dm.(yInd).log eq 1 then yy = alog10( yy )
            
            co = rp.(i).(j).sP.colors[cInds[j]]
            th = !p.thick * thick[j]
                    
            cgPlot,xx,smooth(yy[*,0],sK),line=lines[j],color=co,thick=th,/overplot
          endfor ;n_tags(rp.(i)),j
        endfor ;n_tags(rp),i
        
        count += 1
      endforeach ; tag_names(pConfig_dm)
      
      ; (def) - lower three panels are stellar
      qInd = where( tag_names(rp.(0).(0).stars) eq fName )
      bInd = where( tag_names(rp.(0).(0).stars.(qInd).binCen) eq fName )
      
      pInd = where( tag_names(pConfig_stars) eq fName )
      xtitle = pConfig_stars.(pInd).label
      xrange = pConfig_stars.(pInd).range
    
      foreach curField,tag_names(pConfig_stars) do begin
        if curField eq fName then continue
        
        yInd = where( tag_names(pConfig_stars) eq curField )
        zInd = where( tag_names(rp.(0).(0).stars.(qInd)) eq curField )
        yrange = pConfig_stars.(yInd).range
        ytitle = pConfig_stars.(yInd).label
        
        ; plot
        cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xlog=0,$
          xtitle=textoidl(xtitle),ytitle=textoidl(ytitle),pos=pos[count]+offset,/noerase
          
        foreach radLine,radLines do $
          cgPlot,[radLine,radLine],yrange+[0.05,-0.05],line=1,color='light gray',/overplot
          
        for i=0,n_tags(rp)-1 do begin
          for j=0,n_tags(rp.(i))-1 do begin
            xx = rp.(i).(j).stars.(qInd).binCen.(bInd)
            if pConfig_stars.(pInd).log eq 1 then xx = alog10( xx )
            yy = rp.(i).(j).stars.(qInd).(zInd)
            if pConfig_stars.(yInd).log eq 1 then yy = alog10( yy )
            
            co = rp.(i).(j).sP.colors[cInds[j]]
            th = !p.thick * thick[j]
                    
            cgPlot,xx,smooth(yy[*,0],sK),line=lines[j],color=co,thick=th,/overplot
          endfor ;n_tags(rp.(i)),j
        endfor ;n_tags(rp),i
        
        count += 1
      endforeach ; tag_names(pConfig_dm)
        
    end_PS
  
  endforeach ; tag_names(pConfig_dm)
  
  ; plot (3) - 2d histograms of all quantities (3x2), one per panel, one halo one resolution
  for i=0,n_tags(rp)-1 do begin
    for j=0,n_tags(rp.(i))-1 do begin
    
      plotStr2 = (tag_names(rp))[i] + (tag_names(rp.(i)))[j]
      plotStr2 += '_cutSubS-' + str(cutSubS)
    
      foreach fName,tag_names(pConfig) do begin
      
        qInd = where( tag_names(rp.(0).(0).gas) eq fName )
        bInd = where( tag_names(rp.(0).(0).gas.(qInd).binCen) eq fName )
        mInd1 = where( tag_names(binMinMax) eq fName )
        
        pInd = where( tag_names(pConfig) eq fName )
        xtitle = pConfig.(pInd).label
        xrange = binMinMax.(mInd1)
        
        start_PS, rp.(0).(0).sP.plotPath + 'zoomProfiles_'+fName+'_Gas_2D_' + plotStr2 + '.eps', $
          xs=12.0*1.4, ys=6.0*1.4
        
          pos = plot_pos(col=3,row=2,/gap)
          offset = [-0.03,-0.02,-0.02,0]
          
          count = 0
          
          foreach curField,tag_names(pConfig) do begin
            if curField eq fName then continue
            
            yInd = where( tag_names(pConfig) eq curField )
            mInd2 = where( tag_names(binMinMax) eq curField )
            zInd_1D = where( tag_names(rp.(0).(0).gas.(qInd)) eq curField )
            zInd_2D = where( tag_names(rp.(0).(0).gas.(qInd)) eq curField+'_2D' )
            yrange = binMinMax.(mInd2)
            ytitle = pConfig.(yInd).label
            
            ; plot
            h2d = rp.(i).(j).gas.(qInd).(zInd_2D)
            h2d = colorMapAccTime(h2d, logHist=logHist, byRow=byRow, byCol=byCol, min=0, range=h2dMinMax)
            loadColorTable, ctName2D
            tvim,h2d,scale=0,/c_map,pos=pos[count]+offset,/noframe,/noerase
            
            cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xlog=0,$
              xtitle=textoidl(xtitle),ytitle=textoidl(ytitle),pos=pos[count]+offset,/noerase
               
            foreach radLine,radLines do $
              cgPlot,[radLine,radLine],yrange+[0.05,-0.05],line=1,color='light gray',/overplot
            
            ; mean/median
            xx = rp.(i).(j).gas.(qInd).binCen.(bInd)
            if pConfig.(pInd).log eq 1 then xx = alog10( xx )
            yy = rp.(i).(j).gas.(qInd).(zInd_1D)
            if pConfig.(yInd).log eq 1 then yy = alog10( yy )
            
            cgPlot,xx,smooth(yy[*,0],sK),line=line2D[0],color=lineColor2D[0],thick=thick2D[0],/overplot
            cgPlot,xx,smooth(yy[*,2],sK),line=line2D[1],color=lineColor2D[1],thick=thick2D[1],/overplot
                    
            count += 1
          endforeach ; tag_names(pConfig)
            
        end_PS
        
        ; plot (3b) - individual plots for some quantities        
        foreach curField,indivPlotsY do begin
          if curField eq fName then continue
          if total(fName eq indivPlotsX) eq 0 then continue
          
          start_PS, rp.(0).(0).sP.plotPath + 'zoomProfiles_'+fName+'-'+curField+$
            '_Gas_2D_' + plotStr2 + '.eps' ;, xs=12.0*1.4, ys=6.0*1.4
          
            pos = [0.13,0.14,0.92,0.92]
          
            yInd = where( tag_names(pConfig) eq curField )
            mInd2 = where( tag_names(binMinMax) eq curField )
            zInd_1D = where( tag_names(rp.(0).(0).gas.(qInd)) eq curField )
            zInd_2D = where( tag_names(rp.(0).(0).gas.(qInd)) eq curField+'_2D' )
            yrange = binMinMax.(mInd2)
            ytitle = pConfig.(yInd).label
            
            ; plot
            h2d = rp.(i).(j).gas.(qInd).(zInd_2D)
            h2d = colorMapAccTime(h2d, logHist=logHist, byRow=byRow, byCol=byCol, min=0, range=h2dMinMax)
            loadColorTable, ctName2D
            tvim,h2d,scale=0,/c_map,pos=pos,/noframe,/noerase
            
            cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xlog=0,$
              xtitle=textoidl(xtitle),ytitle=textoidl(ytitle),pos=pos,/noerase
               
            foreach radLine,radLines do $
              cgPlot,[radLine,radLine],yrange+[0.05,-0.05],line=1,color='light gray',/overplot
            
            ; mean/median
            xx = rp.(i).(j).gas.(qInd).binCen.(bInd)
            if pConfig.(pInd).log eq 1 then xx = alog10( xx )
            yy = rp.(i).(j).gas.(qInd).(zInd_1D)
            if pConfig.(yInd).log eq 1 then yy = alog10( yy )
            
            cgPlot,xx,smooth(yy[*,0],sK),line=line2D[0],color=lineColor2D[0],thick=thick2D[0],/overplot
            cgPlot,xx,smooth(yy[*,2],sK),line=line2D[1],color=lineColor2D[1],thick=thick2D[1],/overplot
                  
          end_PS ;(3b)
          
        endforeach ; indivPlotsY (3b)
      
      endforeach ; tag_names(pConfig)
  
    endfor ;j
  endfor ;i
    
  ; plot (4) - 2d histograms of all DM/stars, one per panel, one halo one resolution
  for i=0,n_tags(rp)-1 do begin
    for j=0,n_tags(rp.(i))-1 do begin
    
      plotStr2 = (tag_names(rp))[i] + (tag_names(rp.(i)))[j]
      plotStr2 += '_cutSubS-' + str(cutSubS)
          
      foreach fName,tag_names(pConfig_dm) do begin
          
        start_PS, rp.(0).(0).sP.plotPath + 'zoomProfiles_'+fName+'_DMStars_2D_' + plotStr2 + '.eps', $
          xs=12.0*1.4, ys=6.0*1.4
        
          pos = plot_pos(col=3,row=2,/gap)
          offset = [-0.03,-0.02,-0.02,0]
          
          count = 0
          
          ; (abc) - upper three panels are DM
          qInd = where( tag_names(rp.(0).(0).dm) eq fName )
          bInd = where( tag_names(rp.(0).(0).dm.(qInd).binCen) eq fName )
          mInd1 = where( tag_names(binMinMax) eq fName )
          
          pInd = where( tag_names(pConfig_dm) eq fName )
          xtitle = pConfig_dm.(pInd).label
          xrange = binMinMax.(mInd1)
      
          foreach curField,tag_names(pConfig_dm) do begin
            if curField eq fName then continue
            
            yInd = where( tag_names(pConfig_dm) eq curField )
            mInd2 = where( tag_names(binMinMax) eq curField )
            zInd_1D = where( tag_names(rp.(0).(0).dm.(qInd)) eq curField )
            zInd_2D = where( tag_names(rp.(0).(0).dm.(qInd)) eq curField+'_2D' )
            yrange = binMinMax.(mInd2)
            ytitle = pConfig_dm.(yInd).label
            
            ; plot
            h2d = rp.(i).(j).dm.(qInd).(zInd_2D)
            h2d = colorMapAccTime(h2d, logHist=logHist, byRow=byRow, byCol=byCol, min=0, range=h2dMinMax)
            loadColorTable, ctName2D
            tvim,h2d,scale=0,/c_map,pos=pos[count]+offset,/noframe,/noerase
            
            cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xlog=0,$
              xtitle=textoidl(xtitle),ytitle=textoidl(ytitle),pos=pos[count]+offset,/noerase
               
            foreach radLine,radLines do $
              cgPlot,[radLine,radLine],yrange+[0.05,-0.05],line=1,color='light gray',/overplot
            
            ; mean/median
            xx = rp.(i).(j).dm.(qInd).binCen.(bInd)
            if pConfig_dm.(pInd).log eq 1 then xx = alog10( xx )
            yy = rp.(i).(j).dm.(qInd).(zInd_1D)
            if pConfig_dm.(yInd).log eq 1 then yy = alog10( yy )
            
            cgPlot,xx,smooth(yy[*,0],sK),line=line2D[0],color=lineColor2D[0],thick=thick2D[0],/overplot
            cgPlot,xx,smooth(yy[*,2],sK),line=line2D[1],color=lineColor2D[1],thick=thick2D[1],/overplot
                    
            count += 1
          endforeach ; tag_names(pConfig_dm)
        
          ; (def) - lower three panels are stars
          qInd = where( tag_names(rp.(0).(0).stars) eq fName )
          bInd = where( tag_names(rp.(0).(0).stars.(qInd).binCen) eq fName )
          mInd1 = where( tag_names(binMinMax) eq fName )
          
          pInd = where( tag_names(pConfig_stars) eq fName )
          xtitle = pConfig_stars.(pInd).label
          xrange = binMinMax.(mInd1)
      
          foreach curField,tag_names(pConfig_dm) do begin
            if curField eq fName then continue
            
            yInd = where( tag_names(pConfig_stars) eq curField )
            mInd2 = where( tag_names(binMinMax) eq curField )
            zInd_1D = where( tag_names(rp.(0).(0).stars.(qInd)) eq curField )
            zInd_2D = where( tag_names(rp.(0).(0).stars.(qInd)) eq curField+'_2D' )
            yrange = binMinMax.(mInd2)
            ytitle = pConfig_stars.(yInd).label
            
            ; plot
            h2d = rp.(i).(j).stars.(qInd).(zInd_2D)
            h2d = colorMapAccTime(h2d, logHist=logHist, byRow=byRow, byCol=byCol, min=0, range=h2dMinMax)
            loadColorTable, ctName2D
            tvim,h2d,scale=0,/c_map,pos=pos[count]+offset,/noframe,/noerase
            
            cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xlog=0,$
              xtitle=textoidl(xtitle),ytitle=textoidl(ytitle),pos=pos[count]+offset,/noerase
               
            foreach radLine,radLines do $
              cgPlot,[radLine,radLine],yrange+[0.05,-0.05],line=1,color='light gray',/overplot
            
            ; mean/median
            xx = rp.(i).(j).stars.(qInd).binCen.(bInd)
            if pConfig_stars.(pInd).log eq 1 then xx = alog10( xx )
            yy = rp.(i).(j).stars.(qInd).(zInd_1D)
            if pConfig_stars.(yInd).log eq 1 then yy = alog10( yy )
            
            cgPlot,xx,smooth(yy[*,0],sK),line=line2D[0],color=lineColor2D[0],thick=thick2D[0],/overplot
            cgPlot,xx,smooth(yy[*,2],sK),line=line2D[1],color=lineColor2D[1],thick=thick2D[1],/overplot
                    
            count += 1
          endforeach ; tag_names(pConfig_dm)
        
        end_PS
      
      endforeach ; tag_names(pConfig_dm)
      
    endfor ;j
  endfor ;i
  
  ; plot (5) - correlation matrix (for each halo, at each res)
  for i=0,n_tags(rp)-1 do begin
    for j=0,n_tags(rp.(i))-1 do begin
    
      plotStr2 = (tag_names(rp))[i] + (tag_names(rp.(i)))[j]
      plotStr2 += '_cutSubS-' + str(cutSubS)
    
      start_PS, rp.(0).(0).sP.plotPath + 'zoomMatrix_Gas_2D_' + plotStr2 + '.eps', $
        xs=14.0*1.0, ys=10.0*1.0
    
        foreach fName,matrixPlots,colNum do begin ; columns
        
          qInd = where( tag_names(rp.(0).(0).gas) eq fName )
          bInd = where( tag_names(rp.(0).(0).gas.(qInd).binCen) eq fName )
          mInd1 = where( tag_names(binMinMax) eq fName )
          
          pInd = where( tag_names(pConfig) eq fName )
          xtitle = pConfig.(pInd).label
          xrange = binMinMax.(mInd1)
          xtickf = ""
          xtickv = []
          xticks = 0
          
          if fName eq 'TEMP' then xtickv = [4.5,5.0,5.5,6.0,6.5]
          if fName eq 'VRAD' then xtickv = [-300,-150,0,150] ; make sure inside 2D bounds
          
          if n_elements(xtickv) gt 0 then xticks = n_elements(xtickv)-1
                  
          foreach curField,matrixPlots,rowNum do begin ; rows
            
            yInd = where( tag_names(pConfig) eq curField )
            mInd2 = where( tag_names(binMinMax) eq curField )
            zInd_1D = where( tag_names(rp.(0).(0).gas.(qInd)) eq curField )
            zInd_2D = where( tag_names(rp.(0).(0).gas.(qInd)) eq curField+'_2D' )
            yrange = binMinMax.(mInd2)
            ytitle = pConfig.(yInd).label
            ytickf = ""
            ytickv = []
            yticks = 0
            
            if curField eq 'TEMP' then ytickv = [4.5,5.0,5.5,6.0,6.5]
            if curField eq 'VRAD' then ytickv = [-300,-150,0,150]
            
            if n_elements(ytickv) gt 0 then yticks = n_elements(ytickv)-1
            
            ; calculate position and labelling
            x0 = 0.08 & x1 = 0.96
            y0 = 0.08 & y1 = 0.96
            xs = (x1-x0)/n_elements(matrixPlots)
            ys = (y1-y0)/n_elements(matrixPlots)
            
            pos = [x0+xs*colNum, y0+ys*rowNum, x0+xs*(colNum+1), y0+ys*(rowNum+1)]
            if rowNum gt 0 then xtitle=""
            if rowNum gt 0 then xtickf="(A1)"
            if colNum gt 0 then ytitle=""
            if colNum gt 0 then ytickf="(A1)"
           
            ; plot
            h2d = rp.(i).(j).gas.(qInd).(zInd_2D)
            h2d = colorMapAccTime(h2d, logHist=logHist, byRow=byRow, byCol=byCol, min=0, range=h2dMinMax)
            loadColorTable, ctName2D
            tvim,h2d,scale=0,/c_map,pos=pos,/noframe,/noerase
            
            cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xlog=0,$
              xtitle=textoidl(xtitle),ytitle=textoidl(ytitle),$
              xtickformat=xtickf,ytickformat=ytickf,$
              xtickv=xtickv,xticks=xticks,ytickv=ytickv,yticks=yticks,pos=pos,/noerase
               
            ; mean/median
            xx = rp.(i).(j).gas.(qInd).binCen.(bInd)
            if pConfig.(pInd).log eq 1 then xx = alog10( xx )
            yy = rp.(i).(j).gas.(qInd).(zInd_1D)
            if pConfig.(yInd).log eq 1 then yy = alog10( yy )
            
            ;cgPlot,xx,smooth(yy[*,0],sK),line=line2D[0],color=lineColor2D[0],thick=thick2D[0],/overplot
            cgPlot,xx,smooth(yy[*,2],sK),line=line2D[1],color=lineColor2D[1],thick=thick2D[1],/overplot
                    
          endforeach ; matrixPlots,rowNum
        endforeach ; matrixPlots,colNum
        
        end_PS
    endfor ;j
  endfor ;i
  
  ; plot (5b) - correlation matrix (stacked at each res)
  foreach stackResLevel,resLevels do begin
  
    plotStr2 = 'h'
    foreach hInd,hInds do plotStr2 += str(hInd)
    plotStr2 += '_L' + str(stackResLevel) + 'stack' + '_cutSubS-' + str(cutSubS)
    
      start_PS, rp.(0).(0).sP.plotPath + 'zoomMatrix_Gas_2D_' + plotStr2 + '.eps', $
        xs=14.0*1.0, ys=10.0*1.0
    
        foreach fName,matrixPlots,colNum do begin ; columns
        
          qInd = where( tag_names(rp.(0).(0).gas) eq fName )
          bInd = where( tag_names(rp.(0).(0).gas.(qInd).binCen) eq fName )
          mInd1 = where( tag_names(binMinMax) eq fName )
          
          pInd = where( tag_names(pConfig) eq fName )
          xtitle = pConfig.(pInd).label
          xrange = binMinMax.(mInd1)
          xtickf = ""
          xtickv = []
          xticks = 0
          
          if fName eq 'TEMP' then xtickv = [4.5,5.0,5.5,6.0,6.5]
          if fName eq 'VRAD' then xtickv = [-300,-150,0,150] ; make sure inside 2D bounds
          
          if n_elements(xtickv) gt 0 then xticks = n_elements(xtickv)-1
                  
          foreach curField,matrixPlots,rowNum do begin ; rows
          
            yInd = where( tag_names(pConfig) eq curField )
            mInd2 = where( tag_names(binMinMax) eq curField )
            zInd_1D = where( tag_names(rp.(0).(0).gas.(qInd)) eq curField )
            zInd_2D = where( tag_names(rp.(0).(0).gas.(qInd)) eq curField+'_2D' )
            yrange = binMinMax.(mInd2)
            ytitle = pConfig.(yInd).label
            ytickf = ""
            ytickv = []
            yticks = 0
            
            if curField eq 'TEMP' then ytickv = [4.5,5.0,5.5,6.0,6.5]
            if curField eq 'VRAD' then ytickv = [-300,-150,0,150]
            
            if n_elements(ytickv) gt 0 then yticks = n_elements(ytickv)-1
            
            ; calculate position and labelling
            x0 = 0.05 & x1 = 0.98
            y0 = 0.06 & y1 = 0.98
            xs = (x1-x0)/n_elements(matrixPlots)
            ys = (y1-y0)/n_elements(matrixPlots)
            
            pos = [x0+xs*colNum, y0+ys*rowNum, x0+xs*(colNum+1), y0+ys*(rowNum+1)]
            if rowNum gt 0 then xtitle=""
            if rowNum gt 0 then xtickf="(A1)"
            if colNum gt 0 then ytitle=""
            if colNum gt 0 then ytickf="(A1)"
           
            ; stack
            j = where( tag_names(rp.(0)) eq 'L'+str(stackResLevel), count )
            h2d = rp.(0).(j).gas.(qInd).(zInd_2D) * 0.0
            
            for i=0,n_tags(rp)-1 do begin
              j = where( tag_names(rp.(i)) eq 'L'+str(stackResLevel), count )
              if count eq 0 then message,'Error (just continue)'
              h2d += rp.(i).(j).gas.(qInd).(zInd_2D)
            endfor
            
            ; plot (skip diagonal for stacked)
            if curField ne fName then begin
              h2d = colorMapAccTime(h2d, logHist=logHist, byRow=byRow, byCol=byCol, min=0, range=h2dMinMax)
              loadColorTable, ctName2D
              tvim,h2d,scale=0,/c_map,pos=pos,/noframe,/noerase
            endif
            
            cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xlog=0,$
              xtitle=textoidl(xtitle),ytitle=textoidl(ytitle),$
              xtickformat=xtickf,ytickformat=ytickf,$
              xtickv=xtickv,xticks=xticks,ytickv=ytickv,yticks=yticks,pos=pos,/noerase
               
            ; mean/median
            for i=0,n_tags(rp)-1 do begin
              j = where( tag_names(rp.(i)) eq 'L'+str(stackResLevel), count )
              
              xx = rp.(i).(j).gas.(qInd).binCen.(bInd)
              if pConfig.(pInd).log eq 1 then xx = alog10( xx )
              yy = rp.(i).(j).gas.(qInd).(zInd_1D)
              if pConfig.(yInd).log eq 1 then yy = alog10( yy )
              
              ;cgPlot,xx,smooth(yy[*,2],sK),line=line2D[1],color=lineColor2D[1],thick=thick2D[1],/overplot
            endfor
            
          endforeach ; matrixPlots,rowNum
        endforeach ; matrixPlots,colNum
        
        end_PS
    
  endforeach ; resLevels stack
  
  ; plot (6) - stacked 2d histograms of all quantities, one per panel, all halos at each res
  foreach stackResLevel,resLevels do begin
  
    plotStr2 = 'h'
    foreach hInd,hInds do plotStr2 += str(hInd)
    plotStr2 += '_L' + str(stackResLevel) + 'stack' + '_cutSubS-' + str(cutSubS)
  
    foreach fName,tag_names(pConfig) do begin
    
      qInd = where( tag_names(rp.(0).(0).gas) eq fName )
      bInd = where( tag_names(rp.(0).(0).gas.(qInd).binCen) eq fName )
      mInd1 = where( tag_names(binMinMax) eq fName )
      
      pInd = where( tag_names(pConfig) eq fName )
      xtitle = pConfig.(pInd).label
      xrange = binMinMax.(mInd1)
      
      start_PS, rp.(0).(0).sP.plotPath + 'zoomProfiles_'+fName+'_Gas_2D_' + plotStr2 + '.eps', $
        xs=12.0*1.4, ys=6.0*1.4
      
        pos = plot_pos(col=3,row=2,/gap)
        offset = [-0.03,-0.02,-0.02,0]
        
        counter = 0
        
        foreach curField,tag_names(pConfig) do begin
          if curField eq fName then continue
          
          yInd = where( tag_names(pConfig) eq curField )
          mInd2 = where( tag_names(binMinMax) eq curField )
          zInd_1D = where( tag_names(rp.(0).(0).gas.(qInd)) eq curField )
          zInd_2D = where( tag_names(rp.(0).(0).gas.(qInd)) eq curField+'_2D' )
          yrange = binMinMax.(mInd2)
          ytitle = pConfig.(yInd).label
          
          ; stack
          j = where( tag_names(rp.(0)) eq 'L'+str(stackResLevel), count )
          h2d = rp.(0).(j).gas.(qInd).(zInd_2D) * 0.0
          
          for i=0,n_tags(rp)-1 do begin
            j = where( tag_names(rp.(i)) eq 'L'+str(stackResLevel), count )
            if count eq 0 then message,'Error (just continue)'
            h2d += rp.(i).(j).gas.(qInd).(zInd_2D)
          endfor
          
          ; plot
          h2d = colorMapAccTime(h2d, logHist=logHist, byRow=byRow, byCol=byCol, min=0, range=h2dMinMax)
          loadColorTable, ctName2D
          tvim,h2d,scale=0,/c_map,pos=pos[counter]+offset,/noframe,/noerase
          
          cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xlog=0,$
            xtitle=textoidl(xtitle),ytitle=textoidl(ytitle),pos=pos[counter]+offset,/noerase
             
          foreach radLine,radLines do $
            cgPlot,[radLine,radLine],yrange+[0.05,-0.05],line=1,color='light gray',/overplot
          
          ; mean/median of each run
          for i=0,n_tags(rp)-1 do begin
            j = where( tag_names(rp.(i)) eq 'L'+str(stackResLevel), count )
            
            xx = rp.(i).(j).gas.(qInd).binCen.(bInd)
            if pConfig.(pInd).log eq 1 then xx = alog10( xx )
            yy = rp.(i).(j).gas.(qInd).(zInd_1D)
            if pConfig.(yInd).log eq 1 then yy = alog10( yy )
            
            ;cgPlot,xx,smooth(yy[*,0],sK),line=line2D[0],color=lineColor2D[0],thick=thick2D[0],/overplot
            cgPlot,xx,smooth(yy[*,2],sK),line=line2D[1],color=lineColor2D[1],thick=thick2D[1],/overplot
          endfor
          
          counter += 1
        endforeach ; tag_names(pConfig)
          
      end_PS
      
      ; plot (6b) - individual plots for some quantities        
      foreach curField,indivPlotsY do begin
        if curField eq fName then continue
        if total(fName eq indivPlotsX) eq 0 then continue
        
        start_PS, rp.(0).(0).sP.plotPath + 'zoomProfiles_'+fName+'-'+curField+$
          '_Gas_2D_' + plotStr2 + '.eps' ;, xs=12.0*1.4, ys=6.0*1.4
          
          pos = [0.13,0.14,0.92,0.92]
          
          yInd = where( tag_names(pConfig) eq curField )
          mInd2 = where( tag_names(binMinMax) eq curField )
          zInd_1D = where( tag_names(rp.(0).(0).gas.(qInd)) eq curField )
          zInd_2D = where( tag_names(rp.(0).(0).gas.(qInd)) eq curField+'_2D' )
          yrange = binMinMax.(mInd2)
          ytitle = pConfig.(yInd).label
          
          ; stack
          j = where( tag_names(rp.(0)) eq 'L'+str(stackResLevel), count )
          h2d = rp.(0).(j).gas.(qInd).(zInd_2D) * 0.0
          
          for i=0,n_tags(rp)-1 do begin
            j = where( tag_names(rp.(i)) eq 'L'+str(stackResLevel), count )
            if count eq 0 then message,'Error (just continue)'
            h2d += rp.(i).(j).gas.(qInd).(zInd_2D)
          endfor
          
          ; plot
          h2d = colorMapAccTime(h2d, logHist=logHist, byRow=byRow, byCol=byCol, min=0, range=h2dMinMax)
          loadColorTable, ctName2D
          tvim,h2d,scale=0,/c_map,pos=pos,/noframe,/noerase
          
          cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xlog=0,$
            xtitle=textoidl(xtitle),ytitle=textoidl(ytitle),pos=pos,/noerase
             
          foreach radLine,radLines do $
            cgPlot,[radLine,radLine],yrange+[0.05,-0.05],line=1,color='light gray',/overplot
          
          ; mean/median of each run
          for i=0,n_tags(rp)-1 do begin
            j = where( tag_names(rp.(i)) eq 'L'+str(stackResLevel), count )
            
            xx = rp.(i).(j).gas.(qInd).binCen.(bInd)
            if pConfig.(pInd).log eq 1 then xx = alog10( xx )
            yy = rp.(i).(j).gas.(qInd).(zInd_1D)
            if pConfig.(yInd).log eq 1 then yy = alog10( yy )
            
            ;cgPlot,xx,smooth(yy[*,0],sK),line=line2D[0],color=lineColor2D[0],thick=thick2D[0],/overplot
            cgPlot,xx,smooth(yy[*,2],sK),line=line2D[1],color=lineColor2D[1],thick=thick2D[1],/overplot
          endfor
                    
        end_PS ; (5b)
      endforeach ; indivPlotsY (5b)
    
    endforeach ; tag_names(pConfig)
  
  endforeach ; resLevels stack
    
  !except = 1
  stop

end

; plotZoomSlicePDFs() - PDFs of gas dens,temp,etc in a few radial bins

pro plotZoomSlicePDFs ;, input=input
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  hInds     = [0] ;[0,1,2,3,4,5,6,7,8,9] [0,1,7,8]
  resLevels = [9] ;[9,10,11]
  redshift  = 2.0
  newSaves  = 0 ; override existing saves
  cutSubS   = 1   ; remove substructures?
  
  ; 2D plot config (binMinMax also plot bounds for 2D histos, must rebin if you change)
  binMinMax = { rad   : [0.0,2.0]   ,$
                temp  : [4.0,7.0]   ,$
                dens  : [-6.0,0.0]  ,$
                vrad  : [-400,250]  ,$
                csize : [-1.0,0.75] ,$
                entr  : [5.5,9.5]   ,$
                angm  : [3.0,5.0]   ,$
                nh0   : [-28,-24]    }
                
  ctName2D      = 'ncl/WhiteBlueGreenYellowRed' ;-dnA ; 2d histogram
  lineColor2D   = ['black','black'] ; mean/median line color
  line2D        = [3,0]             ; mean/median line style
  thick2D       = [4.0,4.0]         ; mean/median line thickness
  
  h2dMinMax   = [0,8e7] ; linear color scaling of total mass per pixel (comment out to scale each 
                        ; 2D histogram independently from its min to max)
  logHist     = 0       ; take log of 2d histogram pixel values?
  byRow       = 0       ; independently normalize 2d histograms row by row?
  byCol       = 1       ; independently normalize 2d histograms col by col?

  matrixPlots = ['RAD','TEMP','DENS','VRAD','ENTR','ANGM'] ; correlation matrix quantities
  indivPlotsX = ['RAD']   ; single-panel 2d histos?
  indivPlotsY = ['CSIZE'] ; single-panel 2d histos?
  
  ; plot config
  pConfig = { rad   : { label:"r / r_{vir}",                 range:[0.0,2.0],  log:0 } ,$
              temp  : { label:"log T_{gas} [_{ }K_{ }]",     range:[4.5,6.5],  log:1 } ,$
              dens  : { label:"log n_{gas} [_{ }cm^{-3 }]",  range:[-4,-1],    log:1 } ,$
              vrad  : { label:"v_{rad} [_{ }km/s_{ }]",      range:[-200,100], log:0 } ,$
              csize : { label:"log r_{cell} [_{ }kpc_{ }]",  range:[-1.0,1.0], log:1 } ,$
              entr  : { label:"log S_{gas} [_{ }K cm^{2 }]", range:[7.0,9.0],  log:1 } ,$
              angm  : { label:"j_{gas} [_{ }kpc km/s_{ }]",  range:[3.0,5.0],  log:1 }  }
              
  radLines = [0.15,1.5]
  sK       = 1
  
  ; A (resolutions are different linestyles, all same color and thickness)
  lines    = [1,2,0] ; line style, one per resLevel
  cInds    = [1,1,1] ; color index, one per resLevel
  thick    = [1,1,1] ; times !p.thick

  ; B (resolutions are all solid lines, fainter and thinner for lower res)
  ;lines   = [0,0,0]
  ;thick   = [0.25,0.5,1.0]
  ;cInds   = [2,1,0]
  
  ; load  
  foreach hInd,hInds do begin
    hLoc = {}
        
    foreach resLevel,resLevels do begin
      ; get list of snapshots for this halo in the redshift range
      rLoc = {}
      
      if resLevel eq 11 and total(hInd eq [2,3,4,5,6,9]) gt 0 then continue
      
      sP = simParams(run='zoom_20Mpc',res=resLevel,hInd=hInd,redshift=redshift)     
      saveFilename = sP.derivPath + 'binnedVals/radProfiles.h'+str(hInd)+'.'+sP.savPrefix+str(sP.res) + $
                     '.cutSubS-' + str(cutSubS) + '.sav'
      
      if ~file_test(saveFilename) or newSaves eq 1 then begin
                  
        ; adjust binning based on resolution (increase as cube root of number of particles)
        nBins = nBaseBins * 2.0^(sP.zoomLevel-2)
        
        ; load group catalog, locate halo
        gc = loadGroupCat(sP=sP,readIDs=(cutSubS eq 1),skipIDs=(cutSubS eq 0))
        gcInd = zoomTargetHalo(sP=sP)
        
        rVir = gc.group_r_crit200[gc.subgroupGrNr[gcInd]]
        tVir = codeMassToVirTemp( gc.subgroupMass[gcInd], redshift=sP.redshift )
        pPos = gc.subgroupPos[*,gcInd]
        pVel = gc.subgroupVel[*,gcInd]
        
        print,' h'+str(hInd)+'L'+str(resLevel)+': ['+string(sP.snap,format='(I2)')+'] gcInd = ' + $
              str(gcInd) + ' rVir = ' + string(rVir,format='(f5.1)') + $
              ' tVir = ' + string(alog10(tVir),format='(f3.1)')
        
        ; radial binning, 1D and 2D, all quantity combinations
        halo = { tVir:tVir, rVir:rVir, pPos:pPos, pVel:pVel, gcInd:gcInd, $
                 nBins:nBins, binMinMax:binMinMax, cutSubS:cutSubS }
        
        gas   = zoomRadialBin(sP=sP, gc=gc, partType='gas', halo=halo)
        dm    = zoomRadialBin(sP=sP, gc=gc, partType='dm', halo=halo)
        stars = zoomRadialBin(sP=sP, gc=gc, partType='stars', halo=halo)
        
        ; save all values for this halo+resLevel combination
        rLoc = mod_struct( rLoc, 'sP', sP )
        rLoc = mod_struct( rLoc, 'gas', gas )
        rLoc = mod_struct( rLoc, 'dm',  dm )
        rLoc = mod_struct( rLoc, 'stars', stars )
        rLoc = mod_struct( rLoc, 'halo', halo )
       
        save,rLoc,filename=saveFilename
        print,'Saved file: ['+strmid(saveFilename,strlen(sP.derivPath))+']'
      endif else begin
        print,'RESTORE: ['+strmid(saveFilename,strlen(sP.derivPath))+']'
        restore,saveFilename
      endelse
    
      ; add this resolution level to the halo struct
      hLoc = mod_struct( hLoc, 'L'+str(resLevel), rLoc )     
     
    endforeach ;resLevels
    
    ; add this halo to the save struct
    rp = mod_struct( rp, 'h'+str(hInd), hLoc )
    
  endforeach ;hInds
  
  ; TODO
  
end