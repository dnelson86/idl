; haloCompProj.pro
; cosmological boxes - halo comparison project webpages/helpers
; dnelson sep.2012

; webglCutouts(): make halo centered cutouts for the webGL app

pro webglCutouts;, sP=sP, gcIDs=gcIDs

  compile_opt idl2, hidden, strictarr, strictarrsubs

  sP = simParams(res=512,run='gadget',redshift=2.0)
  gcIDs = [2342] ;g2342 a2132 (z2.304)
  ;if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  sizeFac = 5.0 ; times rvir for the bounding box of each cutout
  
  ; target list
  gc    = loadGroupCat(sP=sP,/skipIDs,/verbose)
  sgcen = subgroupPosByMostBoundID(sP=sP) 

  ; load gas properties
  u     = loadSnapshotSubset(sP=sP,partType='gas',field='u')
  nelec = loadSnapshotSubset(sP=sP,partType='gas',field='nelec')
  pos   = loadSnapshotSubset(sP=sP,partType='gas',field='pos')
  dens  = loadSnapshotSubset(sP=sP,partType='gas',field='dens')
  hsml  = loadSnapshotSubset(sP=sP,partType='gas',field='hsml')
  vel   = loadSnapshotSubset(sP=sP,partType='gas',field='vel')

  ; randomly shuffle the points (break the peano ordering to avoid "square" visualization artifacts)
  print,'shuffling...'
  iseed = 424242L
  sort_inds = sort(randomu(iseed,n_elements(u)))
  
  u     = u[sort_inds]
  nelec = nelec[sort_inds]
  dens  = dens[sort_inds]
  hsml  = hsml[sort_inds]
  pos   = pos[*,sort_inds]
  vel   = vel[*,sort_inds]
  
  sort_inds = !NULL
  
  if sP.trMCPerCell eq 0 then simType = 1 ;gadget
  if sP.trMCPerCell ne 0 then simType = 2 ;arepo
  
  print,'saving...'
  ; loop over all requested halos
  foreach gcID, gcIDs do begin
  
    ; get subhalo position and size of imaging box
    boxCen     = sgcen[*,gcID]
    boxSize    = ceil(sizeFac * gc.group_r_crit200[gc.subgroupGrNr[gcID]] / 10.0) * 10.0
    boxSizeImg = [boxSize,boxSize,boxSize] ; cube
  
    ; make conservative cutout greater than boxsize accounting for periodic (do cube not sphere)
    xDist = pos[0,*] - boxCen[0]
    yDist = pos[1,*] - boxCen[1]
    zDist = pos[2,*] - boxCen[2]
    
    correctPeriodicDistVecs, xDist, sP=sP
    correctPeriodicDistVecs, yDist, sP=sP
    correctPeriodicDistVecs, zDist, sP=sP
    
    rvir = gc.group_r_crit200[gc.subgroupGrNr[gcID]]
    virtemp = codeMassToVirTemp(gc.subgroupmass[gcID],sP=sP)
    rad = reform(sqrt(xDist*xDist + yDist*yDist + zDist*zDist)) / rvir[0]
  
    ; local (cube) cutout
    wCut = where(abs(xDist) le 0.5*boxSize and abs(yDist) le 0.5*boxSize and $
                 abs(zDist) le 0.5*boxSize,nCutout)
    print,gcID,nCutout
           
    ; field subsets   
    loc_u     = u[wCut]
    loc_nelec = nelec[wCut]
    loc_dens  = alog10(dens[wCut] * 1e10) ; log msun/ckpc^3 (i.e. comoving density)
    loc_hsml  = hsml[wCut]
    
    ; derived
    loc_temp = alog10(convertUtoTemp(loc_u,loc_nelec)) ; log T (K)
    loc_pres = alog10(calcPressureCGS(loc_u,loc_dens,sP=sP)) ; log P/k (K/cm^3)
    loc_entr = alog10(calcEntropyCGS(loc_u,loc_dens,sP=sP)) ; log S (K cm^2)
    
    ; relative positions
    loc_pos  = fltarr(3,nCutout)
    loc_pos[0,*] = xDist[wCut] ; delta
    loc_pos[1,*] = yDist[wCut]
    loc_pos[2,*] = zDist[wCut]
    
    xDist = !NULL
    yDist = !NULL
    zDist = !NULL
    
    ; velocities: calculate radial velocity relative to bulk halo motion
    loc_vel = vel[*,wCut]
    
    gVel = gc.subgroupVel[*,gcID]
    loc_vel[0,*] = reform(loc_vel[0,*] - gVel[0])
    loc_vel[1,*] = reform(loc_vel[1,*] - gVel[1])
    loc_vel[2,*] = reform(loc_vel[2,*] - gVel[2])
    
    ; make normalized position vector wrt halo center = vec(r) / ||r|| where r from particle to center
    ; means that radvel<0 is inflow and radvel>0 is outflow
    rnorm0 = reform(loc_pos[0,*] - boxCen[0])
    rnorm1 = reform(loc_pos[1,*] - boxCen[1])
    rnorm2 = reform(loc_pos[2,*] - boxCen[2])
    
    correctPeriodicDistVecs, rnorm0, sP=sP
    correctPeriodicDistVecs, rnorm1, sP=sP
    correctPeriodicDistVecs, rnorm2, sP=sP
    
    ; denominator and do divide
    rnorm = sqrt(rnorm0*rnorm0 + rnorm1*rnorm1 + rnorm2*rnorm2)

    rnorm0 /= rnorm
    rnorm1 /= rnorm
    rnorm2 /= rnorm
    
    ; dot(vel,rnorm) gives the magnitude of the projection of vel onto vec(r)
    loc_radvel = loc_vel[0,*]*rnorm0 + loc_vel[1,*]*rnorm1 + loc_vel[2,*]*rnorm2 ; 1xN

    ; velocities: create endpoint for each position point for the velocity vector line
    ;loc_pos2 = fltarr(3,nCutout)
    ;loc_pos2[0,*] = loc_pos[0,*] + loc_vel[0,*]*velVecFac
    ;loc_pos2[1,*] = loc_pos[1,*] + loc_vel[1,*]*velVecFac
    ;loc_pos2[2,*] = loc_pos[2,*] + loc_vel[2,*]*velVecFac
  
    loc_hsml = reform(loc_hsml,[1,nCutout]) ; make 1xN vector
    loc_temp = reform(loc_temp,[1,nCutout])
    loc_dens = reform(loc_dens,[1,nCutout])
    loc_pres = reform(loc_pres,[1,nCutout])
    loc_entr = reform(loc_entr,[1,nCutout])

    ; prepare output data
    dataout = [loc_pos,loc_hsml,loc_temp,loc_dens,loc_pres,loc_entr,loc_radvel]
    
    ; write binary file
    fileName = 'cutout.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+'.dat'
    
    outStruct = { fileType : fix(1)              ,$ ; 0=gas[x,y,z,hsml,temp,dens]
                                                  $ ; 1=gas[x,y,z,hsml,temp,dens,pres,entr,radvel]
                  simType  : fix(simType)        ,$ ; 1=gadget, 2=arepo
                  nPts     : long(nCutout)       ,$
                  redshift : float(sP.redshift)  ,$
                  hInd     : long(gcID)          ,$ ; subgroup id
                  hVirRad  : float(rvir)         ,$ ; ckpc
                  hVirTemp : float(virtemp)      ,$ ; log K
                  sizeFac  : float(sizeFac)      ,$ ; bounding box size
                  data     : reform(dataout,n_elements(dataout)) $
                }
                
    ; write file
    openw,lun,fileName,/get_lun
    writeu,lun,outStruct
    close,lun
  endforeach
end

; makeHaloComparisonImages(): create a mosaic of halo comparison images between gadget and arepo

pro makeHaloComparisonImages, select=select, redshift=redshift

  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; config
  res = 512  
  minLogMsun = 11.25
  maxLogMsun = 13.60
  
  sPa = simParams(res=res,run='arepo',redshift=redshift)
  sPg = simParams(res=res,run='gadget',redshift=redshift)
  
  ; gadget group catalog
  if select eq 'gadget' then begin
    gcg = loadGroupCat(sP=sPg,/skipIDs)
    priGIDs = gcIDList(gc=gcg,select='pri')
    priGMasses = codeMassToLogMsun(gcg.subgroupMass[priGIDs])
    w = where(priGMasses ge minLogMsun and priGMasses le maxLogMsun,countG)
    gcIDs_gadget = priGIDs[w]
    
    print,'Mapping ['+str(countG)+'] gadget halos above minLogMsun.'
    ;sphMapHalos,sP=sPg,gcIDs=gcIDs_gadget;,/coldOnly
    ;scatterMapHalos,sP=sPg,gcIDs=gcIDs_gadget
    ;scatterMapHalosGasDM,sP=sPg,gcIDs=gcIDs_gadget

    valNames  = ['temp','density','pressure','radvel','entropy','metallicity','angmom']
    foreach valName,valNames do $
      hsv = haloShellValue(sP=sPg,partType='gas',valName=valName,subgroupIDs=gcIDs_gadget,/cutSubS)

  endif
  
  ; arepo group catalog
  if select eq 'arepo' then begin
    gca = loadGroupCat(sP=sPa,/skipIDs)
    priAIDs = gcIDList(gc=gca,select='pri')
    priAMasses = codeMassToLogMsun(gca.subgroupMass[priAIDs])
    w = where(priAMasses ge minLogMsun and priAMasses le maxLogMsun,countA)
    gcIDs_arepo = priAIDs[w]
    
    print,'Mapping ['+str(countA)+'] arepo halos above minLogMsun.'
    ;sphMapHalos,sP=sPa,gcIDs=gcIDs_arepo;,/coldOnly
    ;scatterMapHalos,sP=sPa,gcIDs=gcIDs_arepo
    ;scatterMapHalosGasDM,sP=sPa,gcIDs=gcIDs_arepo

    valNames  = ['temp','density','pressure','radvel','entropy','metallicity','angmom']
    foreach valName,valNames do $
      hsv = haloShellValue(sP=sPa,partType='gas',valName=valName,subgroupIDs=gcIDs_arepo,/cutSubS)

  endif
  
  print,'done.'

end

pro makeHaloComparisonPage
  
  compile_opt idl2, hidden, strictarr, strictarrsubs  
  
  ; config
  res = 512
  redshifts = [0.0,1.0,2.0,3.0]
  
  massBins = [10.75,11.00,11.25,11.50,11.75,12.00,13.60]
  
  axesPairs = ['01','02','12']
  axesNames = ['xy','xz','yz']
  
  redshift = 2.0
  
  haloCounter = 0UL
  
  ; load group catalogs and find halo IDs
  sPa = simParams(res=res,run='arepo',redshift=redshift)
  sPg = simParams(res=res,run='gadget',redshift=redshift)
  
  gcg = loadGroupCat(sP=sPg,/skipIDs)
  gca = loadGroupCat(sP=sPa,/skipIDs)
  
  sgceng = subgroupPosByMostBoundID(sP=sPg)
  sgcena = subgroupPosByMostBoundID(sP=sPa)
  
  priGIDs = gcIDList(gc=gcg,select='pri')
  priGMasses = codeMassToLogMsun(gcg.subgroupMass[priGIDs])   
  priAIDs = gcIDList(gc=gca,select='pri')
  priAMasses = codeMassToLogMsun(gca.subgroupMass[priAIDs])
  
  ; match only for webpage output
  match = findMatchedHalos(sP1=sPa,sP2=sPg)  

  ; loop over mass bins
  for j=0,n_elements(massBins)-2 do begin
    print,massBins[j],massBins[j+1]
    
    ; open file and write header with links
    outputName = "HaloComp_z"+string(redshift,format='(i1)')+"mb"+string(massBins[j],format='(f5.2)')+".htm"
    openw,lun,outputName,/get_lun
    
    printf,lun,"<div class='compHead'><table>"
    printf,lun,"  <tr><th>Redshifts</th><th colspan='"+str(n_elements(massBins)-1)+"'>Log(Msun) Bins</th></tr>"
    foreach red,redshifts,k do begin
      printf,lun,"  <tr><td>"+string(red,format='(f3.1)')+"</td>"
      for mb=0,n_elements(massBins)-2 do begin
        printf,lun,"    <td><a href='HaloComp_z"+string(red,format='(i1)')+"mb"+$
                   string(massBins[mb],format='(f5.2)')+".htm'>"+string(massBins[mb],format='(f5.2)')+$
                   " - "+string(massBins[mb+1],format='(f5.2)')+"</a></td>"
      endfor
      printf,lun,"  </tr>"
    endforeach
    printf,lun,"</table></div>"
    
    printf,lun,""
    printf,lun,"<div class='compBody'>"
    
    ; redshift header
    printf,lun,""
    printf,lun,"<h2>z = "+string(redshift,format='(f3.1)')+"</h2>"
    printf,lun,""
  
    ; halo IDs in this mass bin
    wG = where(priGMasses ge massBins[j] and priGMasses lt massBins[j+1],countG)
    if countG gt 0 then gcIDs_gadget = priGIDs[wG]
    wA = where(priAMasses ge massBins[j] and priAMasses lt massBins[j+1],countA)
    if countA gt 0 then gcIDs_arepo = priAIDs[wA]
    
    ; decide which arepo halos have matching gadget halos
    if countA gt 0 then begin
      matchedInds = match.matchedInds[wA]
      wMatch = where(matchedInds ne -1,countMatch,comp=wNoMatch,ncomp=countNoMatchA)
      
      ; find subset of gadget halos that were not matched
      match,match.matchedInds,wG,inds_match,inds_wG,count=countNoMatchG,/sort
      if countNoMatchG gt 0 then begin
        wGMatched = wG[inds_wG]
        wGNoMatch = removeIntersectionFromB(wGMatched,wG)
        countNoMatchG = n_elements(wGNoMatch)
      endif
    endif else begin
      ; if arepo has no halos in this massbin, there can be no matches
      countMatch = 0
      countNoMatchA = countA
      countNoMatchG = countG
      wGNoMatch = wG
    endelse

    ; mass bin header
    printf,lun,""
    printf,lun,"<a name='z"+str(k)+"_massbin"+str(j)+"'></a>"
    printf,lun,"<h3>"+string(massBins[j],format='(f4.1)')+" < log(M) < "+$
                      string(massBins[j+1],format='(f4.1)')+"</h3>"
    printf,lun,""    
    printf,lun,"<table>"
    printf,lun," <tr><th>Info</th><th>Gadget</th><th>Arepo</th></tr>"
    
    ; matched
    ; -------
    for i=0,countMatch-1 do begin
      haloCounter += 1
      
      ; make stats string
      gcIndA = priAIDs[wA[wMatch[i]]]
      gcIndG = priGIDs[matchedInds[wMatch[i]]]
      
      haloMassA = codeMassToLogMsun(gca.subgroupMass[gcIndA])
      haloMassG = codeMassToLogMsun(gcg.subgroupMass[gcIndG])
      
      virRadA  = gca.group_r_crit200[gca.subgroupGrNr[gcIndA]]
      virRadG  = gcg.group_r_crit200[gcg.subgroupGrNr[gcIndG]]
      virTempA = alog10(codeMassToVirTemp(gca.subgroupMass[gcIndA],sP=sPa))
      virTempG = alog10(codeMassToVirTemp(gcg.subgroupMass[gcIndG],sP=sPg))
      
      statsString = "<span class='haloid'>ID: z"+string(redshift,format='(i1)')+"."+$
                    str(haloCounter)+"</span><br><br> "+$
                    "Gadget:<br> log(M) = "+string(haloMassG,format='(f5.2)')+" <br> "+$
                    " xyz = "+string(sgceng[0,gcIndG],format='(f7.1)')+" "+$
                              string(sgceng[1,gcIndG],format='(f7.1)')+" "+$
                              string(sgceng[2,gcIndG],format='(f7.1)')+" <br> "+$
                    " virRad = "+string(virRadG,format='(i3)')+" ckpc<br> "+$
                    " virTemp = "+string(virTempG,format='(f3.1)')+" <br><br> "+$
                    "Arepo:<br> log(M) = "+string(haloMassA,format='(f5.2)')+" <br> "+$
                    " xyz = "+string(sgcena[0,gcIndA],format='(f7.1)')+" "+$
                              string(sgcena[1,gcIndA],format='(f7.1)')+" "+$
                              string(sgcena[2,gcIndA],format='(f7.1)')+" <br> "+$
                    " virRad = "+string(virRadA,format='(i3)')+" ckpc<br> "+$
                    " virTemp = "+string(virTempA,format='(f3.1)')+" <br> "
                              
      ; make gadget image thumbnail/link string
      gadgetString = ''
      
      foreach axisPair,axesPairs,m do begin        
        ; sphMap
        fpath = "sphMaps.gadget.z"+string(redshift,format='(f3.1)')+"/"
        fname = "map.G"+str(res)+"."+str(sPg.snap)+".h"+str(gcIndG)+".axes"+axisPair+".png"
        gadgetString += axesNames[m]+": <a title='Gadget' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a>"
                        
        ; scatterMap
        fpath = "scatterMaps.gadget.z"+string(redshift,format='(f3.1)')+"/"
        fname = "scatter.G"+str(res)+"."+str(sPg.snap)+".h"+str(gcIndG)+".axes"+axisPair+".png"
        gadgetString += " <a title='Gadget' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a>"
        
        ; gasdm
        fpath = "gasdm.gadget.z"+string(redshift,format='(f3.1)')+"/"
        fname = "gasdm.G"+str(res)+"."+str(sPg.snap)+".h"+str(gcIndG)+".axes"+axisPair+".png"
        gadgetString += " <a title='Gadget' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a><br>"
      endforeach
      
      gadgetString = strmid(gadgetString,0,strlen(gadgetString)-4)
      
      ; make arepo image thumbnail/link string
      arepoString = ''
      
      foreach axisPair,axesPairs,m do begin        
        ; sphMap
        fpath = "sphMaps.arepo.z"+string(redshift,format='(f3.1)')+"/"
        fname = "map.A"+str(res)+"."+str(sPa.snap)+".h"+str(gcIndA)+".axes"+axisPair+".png"
        arepoString += axesNames[m]+": <a title='Arepo' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a>"
                        
        ; scatterMap
        fpath = "scatterMaps.arepo.z"+string(redshift,format='(f3.1)')+"/"
        fname = "scatter.A"+str(res)+"."+str(sPa.snap)+".h"+str(gcIndA)+".axes"+axisPair+".png"
        arepoString += " <a title='Arepo' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a>"
       
        ; gasdm
        fpath = "gasdm.arepo.z"+string(redshift,format='(f3.1)')+"/"
        fname = "gasdm.A"+str(res)+"."+str(sPa.snap)+".h"+str(gcIndA)+".axes"+axisPair+".png"
        arepoString += " <a title='Arepo' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a><br>"
      endforeach
      
      arepoString = strmid(arepoString,0,strlen(arepoString)-4)
      
      printf,lun," <tr class='row'>"
      printf,lun,"  <td class='stcell'>"+statsString+"</td>"
      printf,lun,"  <td class='gacell'>"+gadgetString+"</td>"
      printf,lun,"  <td class='arcell'>"+arepoString+"</td>"
      printf,lun," </tr>"
    endfor
    
    ; arepo (unmatched)
    ; -----------------
    for i=0,countNoMatchA-1 do begin
      haloCounter += 1
      
      ; make stats string
      gcIndA = priAIDs[wA[wNoMatch[i]]]
      
      haloMassA = codeMassToLogMsun(gca.subgroupMass[gcIndA])
      virRadA  = gca.group_r_crit200[gca.subgroupGrNr[gcIndA]]
      virTempA = alog10(codeMassToVirTemp(gca.subgroupMass[gcIndA],sP=sPa))
      
      statsString = "<span class='haloid'>ID: z"+string(redshift,format='(i1)')+"."+$
                    str(haloCounter)+"</span><br><br> "+$
                    "Gadget:<br> <i>Unmatched.</i> <br><br> "+$
                    "Arepo:<br> log(M) = "+string(haloMassA,format='(f5.2)')+" <br> "+$
                    " xyz = "+string(sgcena[0,gcIndA],format='(f7.1)')+" "+$
                              string(sgcena[1,gcIndA],format='(f7.1)')+" "+$
                              string(sgcena[2,gcIndA],format='(f7.1)')+" <br> "+$
                    " virRad = "+string(virRadA,format='(i3)')+" ckpc<br> "+$
                    " virTemp = "+string(virTempA,format='(f3.1)')+" <br> "
                       
      ; make arepo image thumbnail/link string
      arepoString = ''
      
      foreach axisPair,axesPairs,m do begin        
        ; sphMap
        fpath = "sphMaps.arepo.z"+string(redshift,format='(f3.1)')+"/"
        fname = "map.A"+str(res)+"."+str(sPa.snap)+".h"+str(gcIndA)+".axes"+axisPair+".png"
        arepoString += axesNames[m]+": <a title='Arepo' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a>"
                        
        ; scatterMap
        fpath = "scatterMaps.arepo.z"+string(redshift,format='(f3.1)')+"/"
        fname = "scatter.A"+str(res)+"."+str(sPa.snap)+".h"+str(gcIndA)+".axes"+axisPair+".png"
        arepoString += " <a title='Arepo' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a>"
        
        ; gasdm
        fpath = "gasdm.arepo.z"+string(redshift,format='(f3.1)')+"/"
        fname = "gasdm.A"+str(res)+"."+str(sPa.snap)+".h"+str(gcIndA)+".axes"+axisPair+".png"
        arepoString += " <a title='Arepo' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a><br>"
      endforeach
      
      arepoString = strmid(arepoString,0,strlen(arepoString)-4)
      
      printf,lun," <tr class='row'>"
      printf,lun,"  <td class='stcell'>"+statsString+"</td>"
      printf,lun,"  <td class='gacell'><i>(Unmatched)</i></td>"
      printf,lun,"  <td class='arcell'>"+arepoString+"</td>"
      printf,lun," </tr>"
    endfor
    
    ; gadget (unmatched)
    ; -----------------
    for i=0,countNoMatchG-1 do begin
      ; check: was this gadget halo matched by an arepo halo in some other mass bin? if so just skip
      ;gcIndG = priGIDs[wGNoMatch[i]]
      ;w = where(match.matchedInds eq gcIndG,countOMM)
      ;if countOMM gt 0 then continue
      
      haloCounter += 1
      
      ; make stats string
      haloMassG = codeMassToLogMsun(gcg.subgroupMass[gcIndG])
      virRadG   = gcg.group_r_crit200[gcg.subgroupGrNr[gcIndG]]
      virTempG  = alog10(codeMassToVirTemp(gcg.subgroupMass[gcIndG],sP=sPg))
      
      statsString = "<span class='haloid'>ID: z"+string(redshift,format='(i1)')+"."+$
                    str(haloCounter)+"</span><br><br> "+$
                    "Gadget:<br> log(M) = "+string(haloMassG,format='(f5.2)')+" <br> "+$
                    " xyz = "+string(sgceng[0,gcIndG],format='(f7.1)')+" "+$
                              string(sgceng[1,gcIndG],format='(f7.1)')+" "+$
                              string(sgceng[2,gcIndG],format='(f7.1)')+" <br> "+$
                    " virRad = "+string(virRadG,format='(i3)')+" ckpc<br> "+$
                    " virTemp = "+string(virTempG,format='(f3.1)')+" <br><br> "+$
                    "Arepo:<br> <i>Unmatched.</i>  <br> "
                       
      ; make arepo image thumbnail/link string
      gadgetString = ''
      
      foreach axisPair,axesPairs,m do begin        
        ; sphMap
        fpath = "sphMaps.gadget.z"+string(redshift,format='(f3.1)')+"/"
        fname = "map.G"+str(res)+"."+str(sPg.snap)+".h"+str(gcIndG)+".axes"+axisPair+".png"
        gadgetString += axesNames[m]+": <a title='Gadget' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a>"
                        
        ; scatterMap
        fpath = "scatterMaps.gadget.z"+string(redshift,format='(f3.1)')+"/"
        fname = "scatter.G"+str(res)+"."+str(sPg.snap)+".h"+str(gcIndG)+".axes"+axisPair+".png"
        gadgetString += " <a title='Gadget' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a>"
        
        ; gasdm
        fpath = "gasdm.gadget.z"+string(redshift,format='(f3.1)')+"/"
        fname = "gasdm.G"+str(res)+"."+str(sPg.snap)+".h"+str(gcIndG)+".axes"+axisPair+".png"
        gadgetString += " <a title='Gadget' href='"+fpath+fname+"'><img src='"+fpath+"thumbnails/"+fname+".jpg'></a><br>"
      endforeach
      
      gadgetString = strmid(gadgetString,0,strlen(gadgetString)-4)
      
      printf,lun," <tr class='row'>"
      printf,lun,"  <td class='stcell'>"+statsString+"</td>"
      printf,lun,"  <td class='gacell'>"+gadgetString+"</td>"
      printf,lun,"  <td class='arcell'><i>(Unmatched)</i></td>"
      printf,lun," </tr>"
    endfor
    
    printf,lun,"</table>"
    printf,lun,""
    printf,lun,"</div>"
    printf,lun,""
    printf,lun,"<p id='footer'>The end.</p>"
    printf,lun,"</body>"
    printf,lun,"</html>"
     
    ; close file
    free_lun,lun
    
    ; if header.txt exists, cat the two together
    if file_test("header.txt") then begin
      print,' added header'
      spawn,'cat header.txt '+outputName+' > newpage.htm' ;make new
      spawn,'mv newpage.htm '+outputName ;overwrite original
    endif
        
  endfor ;massBins
  
end
