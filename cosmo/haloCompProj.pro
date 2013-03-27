; haloCompProj.pro
; cosmological boxes - halo comparison project webpages/helpers
; dnelson sep.2012

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
      calcMatch,match.matchedInds,wG,inds_match,inds_wG,count=countNoMatchG
      
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
