; zoomProfiles.pro
; 'zoom project' 1D and 2D radial profiles, phase diagrams
; dnelson jan.2015

; zoomRadialBin(): do radial binning (1d and 2d) for multiple particle types

function zoomRadialBin, sP=sP, gc=gc, partType=partType, halo=halo, mode=mode
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits(redshift=sP.redshift)
  
  if n_elements(mode) eq 0 then message,'Error: Specify mode.'
  
  ; setup binning
  binSize = {}
  binCen  = {}
  
  for i=0,n_tags(halo.binMinMax)-1 do begin
    bsLoc = ( halo.binMinMax.(i)[1] - halo.binMinMax.(i)[0] ) / float(halo.nBins)
    binSize = mod_struct( binSize, (tag_names(halo.binMinMax))[i], bsLoc )
    binCen  = mod_struct( binCen, (tag_names(halo.binMinMax))[i], fltarr(halo.nBins) )
  endfor
  
  saveFilename = sP.derivPath + 'cutouts/zoomRadialBin.gc'+str(halo.gcInd)+'.'+sP.savPrefix+str(sP.res) + $
                 '.cutSubS-' + str(halo.cutSubS) + '_'+partType+'_r'+str(fix(halo.binMinMax.rad[1]*100))+$
                 '.sav'

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

    ; halo.pVel already peculiar (no scalefac correction needed)
    ; scalefac in pos cancels with scalefac in rad_pri*rVir (for comoving conversion)
    vrad = ((vel[0,*] - halo.pVel[0]) * pos[0,*] + $ 
            (vel[1,*] - halo.pVel[1]) * pos[1,*] + $
            (vel[2,*] - halo.pVel[2]) * pos[2,*]) $
                / (rad_pri*halo.rVir)
            
    vrad += (rad_pri*halo.rVir) * units.H_z ; add in hubble expansion (neglect mass growth term)
    vrad = reform(vrad)
    
    ; angm calculation: mean magnitude of specific angular momentum = rvec x vel
    ; make velocities relative to bulk halo motion
    vel[0,*] = reform(vel[0,*] - halo.pVel[0])
    vel[1,*] = reform(vel[1,*] - halo.pVel[1])
    vel[2,*] = reform(vel[2,*] - halo.pVel[2])
    
    ; angular momentum magnitude
    jvec = fltarr(3,countRad)
    pos /= units.HubbleParam ; remove little h from Coordinates for angmom
    pos *= scalefac ; convert to peculiar for angmom
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
  
    ; radial restriction?
    if tag_exist(halo,'radRange') then begin
      if halo.radRange eq 0 then restrictLocRad = [0.0,2.0] ; no restriction
      if halo.radRange eq 1 then restrictLocRad = [0.2,2.0]
      if halo.radRange eq 2 then restrictLocRad = [0.2,1.0]
      if halo.radRange eq 3 then restrictLocRad = [1.1,1.25]
      
      w = where(data.rad ge restrictLocRad[0] and data.rad le restrictLocRad[1],countLocRad)
      print,'Restrict: radRange='+str(halo.radRange)+' found ['+str(countLocRad)+$
            '] of ['+str(n_elements(data.rad))+']'
            
      if partType eq 'gas' then $
        data = {rad:data.rad[w],dens:data.dens[w],angm:data.angm[w],vrad:data.vrad[w],$
                temp:data.temp[w],csize:data.csize[w],entr:data.entr[w],nh0:data.nh0[w]}
      if partType eq 'dm' then $
        data = {rad:data.rad[w],dens:data.dens[w],angm:data.angm[w],vrad:data.vrad[w]}
      if partType eq 'stars' then $
        data = {rad:data.rad[w],dens:data.dens[w],angm:data.angm[w],vrad:data.vrad[w]}
    endif
    
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
    foreach radRange,halo.radRanges,m do begin
      ; make struct (6 per quantity = [mean,lower quartile,median,upper quartile,count,mass frac])
      r = {}
      
      foreach fName,tag_names(halo.binMinMax) do begin
        if total(fName eq gasOnlyFields gt 0) and partType ne 'gas' then continue
        r = mod_struct( r, fName, fltarr(halo.nBins,6) )
      endforeach
      
      ; radial restriction
      w_rad = where(data.rad ge radRange[0] and data.rad lt radRange[1],count_rad)
      if count_rad eq 0 then message,'Error'
      
      ; loop over each quantity
      foreach fieldName,tag_names(halo.binMinMax) do begin
        if total(fieldName eq gasOnlyFields gt 0) and partType ne 'gas' then continue
        
        ; indices
        tD = (where( strcmp(tag_names(data),fieldName) eq 1 ))[0]
        tH = (where( strcmp(tag_names(halo.binMinMax),fieldName) eq 1 ))[0]
        tR = (where( strcmp(tag_names(r),fieldName) eq 1 ))[0]
        
        ; local data
        data_loc = data.(tD)[w_rad]
        mass_loc = mass[w_rad]
        mass_tot = total( mass_loc )
        
        print,' '+fieldName
      
        for i=0,halo.nBins-1 do begin
          ; bin boundaries
          locMin = 10.0^( halo.binMinMax.(tH)[0] + binSize.(tH)*i )
          locMax = 10.0^( halo.binMinMax.(tH)[0] + binSize.(tH)*(i+1) )

          if fieldName eq 'RAD' or fieldName eq 'VRAD' then begin ; only linear ones
            locMin = halo.binMinMax.(tH)[0] + binSize.(tH)*i
            locMax = halo.binMinMax.(tH)[0] + binSize.(tH)*(i+1)
          endif
          
          ; binning
          w = where(data_loc ge locMin and data_loc lt locMax,count)
          
          binCen.(tH)[i] = mean([locMin,locMax])
          
          if count eq 0 then continue
            
          r.(tR)[i,0]   = mean( data_loc[w] )
          r.(tR)[i,1:3] = percentiles( data_loc[w] )
          r.(tR)[i,4]   = count
          r.(tR)[i,5]   = total( mass_loc[w] ) / mass_tot
          
        endfor ;nBins,i
      
      endforeach ;fieldName
      
      r = mod_struct( r, 'binCen', binCen )
      r = mod_struct( r, 'binSize', binSize )
      r = mod_struct( r, 'radRange', radRange )
      
      rr = mod_struct( rr, 'range'+str(m), r )
    endforeach ; binFieldNames
  endif
  
  ; just print some stats
  if mode eq 2 then begin
    if partType ne 'gas' then return,0
    print,'h'+str(sP.hInd)+' L'+str(sP.levelMax)
    
    ; total gas mass within 1.0 rvir
    ww = where(data.rad ge 0.15 and data.rad lt 1.5,count)
    print,total( mass[ww] )/1e10
    return,total( mass[ww] )/1e10
  endif
  
  return, rr
        
end

; plotZoomRadialProfiles() - 1d with all halos on same plot, and 2d with each halo in separate panel

pro plotZoomRadialProfiles ;, input=input
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  hInds     = [0,1,7,8,2,4,5,9] ;re-ordered such that hIndDisp are in order ;[3,6]
  resLevels = [9,10,11] ;[9,10,11]
  redshift  = 2.0
  newSaves  = 0 ; override existing saves
  
  ;hInds = [ hInds[input] ] ; CAREFUL
  
  ; binning config
  nBaseBins = 100 ; for base resolution level, increased at higher resolutions
  cutSubS   = 1   ; remove substructures?
  radRange  = 0   ; 0=[0,2], 1=[0.2,2], 2=[0.2,1.0], 3=[1.1,1.25] rvir
  
  ; 2D plot config (binMinMax also plot bounds for 2D histos, must rebin if you change)
  binMinMax = { rad   : [0.0,2.0]   ,$
                temp  : [3.47,7.0]  ,$ ; 3000K star forming gas = 3.477 (want it to show up)
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
  pConfig = { rad   : { label:"r / r_{vir}",                    range:[0.0,2.0],  log:0 } ,$
              temp  : { label:"log T_{gas} [_{ }K_{ }]",        range:[4.5,6.5],  log:1 } ,$
              dens  : { label:"log n_{gas} [_{ }cm^{-3 }]",     range:[-5,-2],    log:1 } ,$
              vrad  : { label:"v_{rad} [_{ }km/s_{ }]",         range:[-200,100], log:0 } ,$
              csize : { label:"log r_{cell} [_{ }kpc_{ }]",     range:[-1.0,0.5], log:1 } ,$
              entr  : { label:"log S_{gas} [_{ }K cm^{2 }]",    range:[7.5,9.0],  log:1 } ,$
              angm  : { label:"log j_{gas} [_{ }kpc km/s_{ }]", range:[3.0,4.5],  log:1 }  }
              
  pConfig_dm = { rad   : { label:"r / r_{vir}",                   range:[0.0,2.0],  log:0 } ,$
                 dens  : { label:"log n_{DM} [_{ }cm^{-3 }]",     range:[-6.0,1.0], log:1 } ,$
                 vrad  : { label:"v_{rad,DM} [_{ }km/s_{ }]",     range:[-100,200], log:0 } ,$
                 angm  : { label:"log j_{DM} [_{ }kpc km/s_{ }]", range:[3.0,5.0],  log:1 }  }
                 
  pConfig_stars = { rad   : { label:"r / r_{vir}",                      range:[0.0,2.0],  log:0 } ,$
                    dens  : { label:"log n_{stars} [_{ }cm^{-3 }]",     range:[-6.0,1.0], log:1 } ,$
                    vrad  : { label:"v_{rad,stars} [_{ }km/s_{ }]",     range:[-100,200], log:0 } ,$
                    angm  : { label:"log j_{stars} [_{ }kpc km/s_{ }]", range:[3.0,5.0],  log:1 }  }
         
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
  
  ; test
  masses = fltarr(8,3)
  foreach hInd,hInds,i do begin
    foreach resLevel,resLevels,j do begin
      sP  = simParams(run='zoom_20Mpc',res=resLevel,hInd=hInd,redshift=redshift)
      gc = {}
      halo = {binMinMax:binMinMax,nBins:nBaseBins,gcInd:0,cutSubS:cutSubS}
      gas = zoomRadialBin(sP=sP, gc=gc, partType='gas', halo=halo, mode=2)
      masses[i,j] = gas
    endforeach
  endforeach    
      
  stop
  
  ; load  
  foreach hInd,hInds do begin
    hLoc = {}
        
    foreach resLevel,resLevels do begin
      ; get list of snapshots for this halo in the redshift range
      rLoc = {}
      
      if resLevel eq 11 and total(hInd eq [3,6]) gt 0 then continue
      
      rrTag = ''
      if radRange gt 0 then rrTag = '_rr'+str(radRange)      
      
      sP = simParams(run='zoom_20Mpc',res=resLevel,hInd=hInd,redshift=redshift)     
      saveFilename = sP.derivPath + 'binnedVals/radProfiles.h'+str(hInd)+'.'+sP.savPrefix+str(sP.res) + $
                     '.cutSubS-' + str(cutSubS) + rrTag + '.sav'
      
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
        mode = 0
        halo = { tVir:tVir, rVir:rVir, pPos:pPos, pVel:pVel, gcInd:gcInd, $
                 nBins:nBins, binMinMax:binMinMax, cutSubS:cutSubS, radRange:radRange }
        
        gas   = zoomRadialBin(sP=sP, gc=gc, partType='gas', halo=halo, mode=mode)
        dm    = zoomRadialBin(sP=sP, gc=gc, partType='dm', halo=halo, mode=mode)
        stars = zoomRadialBin(sP=sP, gc=gc, partType='stars', halo=halo, mode=mode)

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
  hColors = []
  hNames  = []
  rNames  = 'L' + str(resLevels)
  rLines  = lines[0:n_elements(resLevels)-1]
  rThick  = thick[0:n_elements(resLevels)-1]

  foreach hInd,hInds,i do begin
    plotStr += str(hInd)
    ;hColors = [hColors,rp.(i).(0).sP.colors[cInds[-1]]] ; old colors
    ;hNames  = [hNames,'h'+str(rp.(i).(0).sP.hIndDisp)] ; need to remake saves for new sP
    
    ; hack: grab new sP for hColor and hName, replace stored color!
    sP = simParams(run='zoom_20Mpc',res=resLevels[0],hInd=hInd,redshift=redshift)
    hColors = [hColors,sP.colors[cInds[-1]]]
    hNames  = [hNames,'h'+str(sP.hIndDisp)]
    
    foreach resLevel,resLevels,j do $
      rp.(i).(j).sP.colors = sP.colors
    
  endforeach
  plotStr += '_L'
  foreach resLevel,resLevels do plotStr += str(resLevel)
  plotStr += '_cutSubS-' + str(cutSubS)
  
  rrTag = ''
  if radRange gt 0 then rrTag = '_rr'+str(radRange)
  
  !except = 0
  
  ; plot (1) - 1D radial (median) profiles of all quantities, one per panel, all halos all resolutions
  foreach fName,['RAD'] do begin ;tag_names(pConfig) do begin
  
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
                    
            cgPlot,xx,smooth(yy[*,2],sK),line=lines[j],color=co,thick=th,/overplot
          endfor ;n_tags(rp.(i)),j
        endfor ;n_tags(rp),i
        
        if count eq 0 then legend,hNames,textcolor=hColors,/top,/right
        if count eq 1 then legend,rNames,linestyle=rLines,thick=rThick*!p.thick,/top,/right
        
        count += 1
      endforeach ; tag_names(pConfig)
        
    end_PS
    
    ; plot (1b) - lower resolutions split into subpanels
    start_PS, rp.(0).(0).sP.plotPath + 'zoomProfilesSubPanels_'+fName+'_Gas_1D_' + plotStr + '.eps', $
      xs=8.0*1.4, ys=9.0*1.4
    
      pos = plot_pos(col=2,row=3,/gap)
      subSize = 0.1
      offset = [-0.02,0+subSize,0.02,0]
      
      count = 0
      
      foreach curField,tag_names(pConfig) do begin
        if curField eq fName then continue
        
        yInd = where( tag_names(pConfig) eq curField )
        zInd = where( tag_names(rp.(0).(0).gas.(qInd)) eq curField )
        yrange = pConfig.(yInd).range
        ytitle = pConfig.(yInd).label
        
        ; plot (i) - main panels
        if curField eq 'CSIZE' then begin
          cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xlog=0,$
            xtitle="",xtickname=replicate(" ",10),ytitle=textoidl(ytitle),pos=pos[count]+offset,$
            ytickv=[-1.0,-0.5,0.0,0.5],yticks=3,/noerase

        endif else begin
          if curField eq 'ENTR' then begin
            cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xlog=0,$
              xtitle="",xtickname=replicate(" ",10),ytitle=textoidl(ytitle),pos=pos[count]+offset,$
              ytickv=[7.5,8.0,8.5,9.0],yticks=3,/noerase
          endif else begin
            if curField eq 'ANGM' then begin
            cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xlog=0,$
              xtitle="",xtickname=replicate(" ",10),ytitle=textoidl(ytitle),pos=pos[count]+offset,$
              ytickv=[3.0,3.5,4.0,4.5],yticks=3,/noerase
            endif else begin
              cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xlog=0,$
                xtitle="",xtickname=replicate(" ",10),ytitle=textoidl(ytitle),pos=pos[count]+offset,/noerase
            endelse
          endelse
        endelse
        
          
        foreach radLine,radLines do $
          cgPlot,[radLine,radLine],yrange+[0.05,-0.05],line=1,color='light gray',/overplot
          
        for i=0,n_tags(rp)-1 do begin
        
          j = where( tag_names(rp.(i)) eq 'L11', count_ind )
          if count_ind eq 0 then continue
          
          ;for j=0,n_tags(rp.(i))-1 do begin
            xx = rp.(i).(j).gas.(qInd).binCen.(bInd)
            if pConfig.(pInd).log eq 1 then xx = alog10( xx )
            yy = rp.(i).(j).gas.(qInd).(zInd)
            if pConfig.(yInd).log eq 1 then yy = alog10( yy )
            
            co = rp.(i).(j).sP.colors[cInds[j]]
            th = !p.thick * thick[j]
                    
            cgPlot,xx,smooth(yy[*,2],sK+2,/nan),line=lines[j],color=co,thick=th,/overplot
          ;endfor ;n_tags(rp.(i)),j
        endfor ;n_tags(rp),i
        
        ; some powerlaw scalings
        if fName eq 'RAD' then begin
        
          mass_codeunits = 10^(12.0) / 1e10 / 0.7
          units = getUnits()
        
          if curField eq 'TEMP' then begin
            
            ; e.g. suto 1998
            xpl = linspace(0.01,1.25,100)
            nfw = nfw_profile(xpl, mass=mass_codeunits, redshift=redshift)
            yy = alog10( 1e6*(nfw.menc_dm/nfw.rad) )
            cgPlot,xpl,yy,line=1,/overplot
            cgText,0.9,6.2,textoidl("M_{enc}(r)/r"),/data
            
            ; r^(-p)
            xpl = linspace(0.01,1.245,100)
            yy = alog10( 4e5*xpl^(-0.5) )
            cgPlot,xpl,yy,line=2,thick=!p.thick+2,/overplot
            cgText,0.4,5.5,textoidl('r^{-1/2}'),/data
          endif
          
          if curField eq 'DENS' then begin
            
            ; r^(-p)
            xpl = linspace(0.1,1.0,100)
            yy = alog10( 0.0002*xpl^(-2.5) )
            cgPlot,xpl,yy,line=1,thick=!p.thick+2,/overplot
            cgText,0.6,-2.8,textoidl('r^{-5/2}'),/data
            
            ; suto isothermal
            xpl = linspace(0.05,1.49,100)
            nfw_dm  = nfw_profile(xpl, mass=mass_codeunits, redshift=redshift)
            xsuto = (nfw_dm.rad/nfw_dm.r_s)  
  
            m1_rscaling = alog(1+xsuto) - xsuto/(1+xsuto)
            ;B_rscaling = (nfw_dm.menc_dm/nfw_dm.rad) / m1_rscaling ; Tvir(r)/m(r)
            B_rscaling = (7.0+1.0)*1.0
            rho_gas_iso = 0.008 * exp(-B_rscaling*( 1 - alog(1 + xsuto) / xsuto ))
            
            yy = alog10( rho_gas_iso )
            cgPlot,xpl,yy,line=2,thick=!p.thick+2,/overplot
            cgText,0.6,-4.5,'Suto',/data
          endif
          
          if curField eq 'VRAD' then begin
            xpl = linspace(1.0,2.0,100)
            
            ; r^(-p)
            ;yy = -10.0*xpl^(-1.5)
            ;cgPlot,xpl,yy,line=1,thick=!p.thick+2,/overplot
            ;cgText,1.5,30,textoidl('r^{-3/2}'),/data
          endif
          
          if curField eq 'CSIZE' then begin
            xpl = linspace(0.2,1.4,100)
            
            ; r^(-p)
            yy = alog10( 1.0*xpl^(5.0/6.0) )
            cgPlot,xpl,yy,line=2,thick=!p.thick+2,/overplot
            cgText,0.6,-0.4,textoidl('r^{5/6}'),/data
          endif
          
          if curField eq 'ENTR' then begin
            xpl = linspace(0.1,1.0,100)
            yy = alog10( 7e8*xpl^(1.1) )
            cgPlot,xpl,yy,line=2,thick=!p.thick+2,/overplot
            cgText,0.6,8.8,textoidl('r^{1.1}'),/data
            
            xpl = linspace(1.20,1.95,100)
            yy = alog10( 9e8*xpl^(-2.5) )
            cgPlot,xpl,yy,line=1,thick=!p.thick+2,/overplot
            cgText,1.3,8.8,textoidl('r^{-5/2}'),/data
          endif
          
          if curField eq 'ANGM' then begin
            ; sqrt(Menc/r)
            xpl = linspace(0.18,1.4,100)
            nfw = nfw_profile(xpl, mass=mass_codeunits, redshift=redshift)
            yy = alog10( 1.2e2*sqrt(nfw.menc_dm/nfw.rad)*nfw.rad )
            cgPlot,xpl,yy,line=2,/overplot
            cgText,0.6,4.2,textoidl("r \times v_c(r)"),/data
           
            ; r^0.5
            xpl = linspace(0.2,1.4,100)
            yy = alog10( 5e3*xpl^(0.5) )
            cgPlot,xpl,yy,line=1,thick=!p.thick+2,/overplot
            cgText,0.8,3.4,textoidl('r^{1/2}'),/data
          endif
        
        endif
        
        ;if count eq 1 then legend,rNames,linestyle=rLines,thick=rThick*!p.thick,/top,/right
        loadColorTable,'bw linear'
        if count eq 3 then legend,rNames[2],linestyle=rLines[2],thick=rThick[2]*!p.thick,/top,/right
        if count eq 1 then legend,hNames[0:3],textcolor=hColors[0:3],pos=[1.64,-2.1],/data
        if count eq 1 then legend,hNames[4:-1],textcolor=hColors[4:-1],pos=[1.76,-2.1],/data
        
        ; plot (ii) - subpanels, ratios L9/L11 and L10/L11
        yrangeSub = [0.1,10.0]
        ytickvSub = [0.2,0.5,1.0,2.0,5.0]
        ytitleSub = "Ln/L11" ;"Ratio"
        posSub = pos[count]+offset
        posSub[1] -= subSize ; move y of lowerleft down
        posSub[3] = posSub[1] + subSize ; set height of subpanel
        
        cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeSub,xs=5,ys=5,xlog=0,ylog=1,yminor=1,$
          xtitle="",ytitle="",$
          ytickv=ytickvSub,yticks=n_elements(ytickvSub)-1,pos=posSub,/noerase
          
        foreach yLine,[1.0,2.0,0.5] do $
          cgPlot,xrange,[yLine,yLine],line=0,color=cgColor('light gray'),/overplot
        
        cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeSub,/xs,/ys,xlog=0,ylog=1,yminor=1,$
          xtitle=textoidl(xtitle),ytitle=textoidl(ytitleSub),$
          ytickv=ytickvSub,yticks=n_elements(ytickvSub)-1,pos=posSub,/noerase
        
        yy_stack = fltarr(nBaseBins,3)
        
        for i=0,n_tags(rp)-1 do begin
          ind_L11 = where( tag_names(rp.(i)) eq 'L11', count_ind )
          if count_ind eq 0 then continue
          
          yy_L11 = rp.(i).(ind_L11).gas.(qInd).(zInd)
          yy_L11 = rebin(reform(yy_L11[*,2]),nBaseBins)
          
          yy_stack[*,2] += yy_L11
          
          for j=0,n_tags(rp.(i))-2 do begin
            xx = rp.(i).(j).gas.(qInd).binCen.(bInd)
            xx = rebin(xx,nBaseBins)
            if pConfig.(pInd).log eq 1 then xx = alog10( xx )
            
            yy = rp.(i).(j).gas.(qInd).(zInd)
            yy = rebin(reform(yy[*,2]),nBaseBins)
            
            yy_stack[*,j] += yy
            yy /= yy_L11 ; plot ratio to highest resolution level
            
            co = rp.(i).(j).sP.colors[cInds[j]]
            th = !p.thick * thick[j]
                    
            cgPlot,xx,smooth(yy,sK+2,/nan),line=lines[j],color=co,thick=th,psym=psym,/overplot
          endfor ;n_tags(rp.(i)),j
        endfor ;n_tags(rp),i
        
        if count eq 3 then legend,[rNames[0:1]+'/'+rNames[2]],$
          linestyle=rLines[0:1],thick=rThick[0:1]*!p.thick,/bottom,/right
        
        ; stacked
        for j=0,n_tags(rp.(0))-2 do begin
          yy_stack[*,j] /= reform(yy_stack[*,2])
          cgPlot,xx,yy_stack[*,j],line=lines[j],thick=fix(!p.thick*thick[j]*3),color='black',/overplot
          print,curField,j
          print,reform(yy_stack[*,j])
          ww = where(xx ge 0.15 and xx le 1.5)
          print,mean(yy_stack[*,j]),mean(yy_stack[ww,j])
        endfor
          
        count += 1
      endforeach ; tag_names(pConfig)
        
    end_PS
  
  endforeach ; tag_names(pConfig)
  stop
  
  if 0 then begin
  ; plot (2) - 1D radial (median) profiles of DM/stars quantities, one per panel, all halos all resolutions
  foreach fName,['RAD'] do begin ;tag_names(pConfig_dm) do begin
  
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
                    
            cgPlot,xx,smooth(yy[*,2],sK),line=lines[j],color=co,thick=th,/overplot
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
                    
            cgPlot,xx,smooth(yy[*,2],sK),line=lines[j],color=co,thick=th,/overplot
          endfor ;n_tags(rp.(i)),j
        endfor ;n_tags(rp),i
        
        count += 1
      endforeach ; tag_names(pConfig_dm)
        
    end_PS
  
  endforeach ; tag_names(pConfig_dm)
  stop
  
  ; plot (3) - 2d histograms of all quantities (3x2), one per panel, one halo one resolution
  for i=0,n_tags(rp)-1 do begin
    for j=0,n_tags(rp.(i))-1 do begin
    
      plotStr2 = (tag_names(rp))[i] + (tag_names(rp.(i)))[j]
      plotStr2 += '_cutSubS-' + str(cutSubS) + rrTag
    
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
      plotStr2 += '_cutSubS-' + str(cutSubS) + rrTag
          
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
  
  endif ;0
  
  ; plot (5) - correlation matrix (for each halo, at each res)
  for i=0,n_tags(rp)-1 do begin
    for j=0,n_tags(rp.(i))-1 do begin
    
      plotStr2 = (tag_names(rp))[i] + (tag_names(rp.(i)))[j]
      plotStr2 += '_cutSubS-' + str(cutSubS) + rrTag
    
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
          
          if fName eq 'TEMP' then xtickv = [4.0,4.5,5.0,5.5,6.0,6.5]
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
            
            if curField eq 'TEMP' then ytickv = [4.0,4.5,5.0,5.5,6.0,6.5]
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
    plotStr2 += '_L' + str(stackResLevel) + 'stack' + '_cutSubS-' + str(cutSubS) + rrTag
    
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
          
          if fName eq 'TEMP' then xtickv = [4.0,4.5,5.0,5.5,6.0,6.5]
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
            
            if curField eq 'TEMP' then ytickv = [4.0,4.5,5.0,5.5,6.0,6.5]
            if curField eq 'VRAD' then ytickv = [-300,-150,0,150]
            
            if n_elements(ytickv) gt 0 then yticks = n_elements(ytickv)-1
            
            ; calculate position and labelling
            x0 = 0.06 & x1 = 0.94
            y0 = 0.06 & y1 = 0.94
            xs = (x1-x0)/n_elements(matrixPlots)
            ys = (y1-y0)/n_elements(matrixPlots)
            
            pos = [x0+xs*colNum, y0+ys*rowNum, x0+xs*(colNum+1), y0+ys*(rowNum+1)]
            xtitle_p = xtitle
            ytitle_p = ytitle
            if rowNum gt 0 then xtitle_p=""
            if rowNum gt 0 then xtickf="(A1)"
            if colNum gt 0 then ytitle_p=""
            if colNum gt 0 then ytickf="(A1)"
           
            ; stack
            j = where( tag_names(rp.(0)) eq 'L'+str(stackResLevel), count )
            h2d = rp.(0).(j).gas.(qInd).(zInd_2D) * 0.0
            
            for i=0,n_tags(rp)-1 do begin
              j = where( tag_names(rp.(i)) eq 'L'+str(stackResLevel), count )
              if count eq 0 then message,'Error (just continue)'
              h2d += rp.(i).(j).gas.(qInd).(zInd_2D)
            endfor
            
            ; for TEMP row, colormap without regard to bottom bin (eEOS 3000K gas)
            if curField eq 'TEMP' then begin
              print,'col='+fName+' row='+curField
              
              for cc=0,n_elements(h2d[0,*])-1 do $
                h2d[cc,11:25] += h2d[cc,0] / 30.0 ; spread it out over some bins
              h2d[*,0] = 0.0
            endif
            
            ; plot (skip diagonal for stacked)
            if curField ne fName then begin
              h2d = colorMapAccTime(h2d, logHist=logHist, byRow=byRow, byCol=byCol, min=0, range=h2dMinMax)
              loadColorTable, ctName2D
              tvim,h2d,scale=0,/c_map,pos=pos,/noframe,/noerase
            endif
            
            ;if curField eq 'TEMP' then stop
            
            cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,xlog=0,$
              xtitle=textoidl(xtitle_p),ytitle=textoidl(ytitle_p),$
              xtickformat=xtickf,ytickformat=ytickf,$
              xtickv=xtickv,xticks=xticks,ytickv=ytickv,yticks=yticks,pos=pos,/noerase
               
            ; add duplicate axes on right and top of figure
            if rowNum eq n_elements(matrixPlots)-1 then begin
              cgAxis,xaxis=1,xrange=xrange,xlog=0,/xs,title=textoidl(xtitle),$
                xtickv=xtickv,xticks=xticks;,xtickformat=xtickf
              ;cgText,mean([pos[0],pos[2]]),0.98,textoidl(xtitle),alignment=0.5
            endif
            if colNum eq n_elements(matrixPlots)-1 then begin
              cgAxis,yaxis=1,yrange=yrange,/ys,title=textoidl(ytitle),$
                ytickv=ytickv,yticks=yticks;,ytickformat=ytickf
              ;cgText,0.98,mean([pos[1],pos[3]]),"test"
            endif
               
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
  stop
  
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
