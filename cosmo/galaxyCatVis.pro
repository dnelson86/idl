; galaxyCatVis.pro
; 2d visualization (using cosmoVis.pro) specialized for the galaxyCat
; dnelson apr.2013

; scatterMapHalos: plot colored scatter plots with velocity vectors on boxes centered on halos
; showGmem/showAll (use galaxy catalog and display gmem only, or all in spatial subset)

pro scatterMapHalos, showGmem=showGmem, showAll=showAll

  compile_opt idl2, hidden, strictarr, strictarrsubs
  if ~keyword_set(showGmem) and ~keyword_set(showAll) then message,'Choose one.'
  
  sP = simParams(res=512,run='feedback',redshift=2.0)
  units = getUnits()

  haloID = 304 ;z2.304 z2.301 z2.130 z2.64
  gcID = getMatchedIDs(sPa=sP,sPg=sP,haloID=haloID)
  gcIDs = [gcID.a]

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  singleColorScale = 0 ; 1=use same color scale for right panel, 0=rescale
  secondGt         = 0 ; 1=show greater than cut, 0=show less than cut
  axes             = list([0,1],[0,2],[1,2]) ;xy,xz,yz
  nbottom          = 50
  sizeFac          = 3.5 ; times rvir
  ctName           = 'helix'
  
  ; use which field and cut value for right panel?
  ;secondField = 'temp'     & secondCutVal = 5.0
  ;secondField = 'entropy'  & secondCutVal = 7.5
  secondField = 'metal'     & secondCutVal = -1.0
  ;secondField = 'vrad'      & secondCutVal = -200.0
  ;secondField = 'vradnorm'  & secondCutVal = 2.0
  ;secondField = 'coolTime'  & secondCutVal = 1.0 ; gmem only
  ;secondField = 'dynTime'   & secondCutVal = 0.4 ; gmem only
  ;secondField = 'timeRatio' & secondCutVal = 1.0 ; gmem only
  
  ; use which field and minmax for color mapping?
  ;colorField = 'temp'     & fieldMinMax  = [4.0,7.0]
  ;colorField = 'entropy'  & fieldMinMax  = [6.5,8.5] ; log(CGS)
  colorField = 'metal'     & fieldMinMax = [-4.0,1.0] ; log(Z/Zsun)
  ;colorField = 'vrad'      & fieldMinMax = [-400,400] ; km/s
  ;colorField = 'vradnorm'  & fieldMinMax = [0.0,4.0] ; vrad/v200
  ;colorField  = 'coolTime' & fieldMinMax = [0.0,5.0] ; gmem only
  ;colorField = 'dynTime'   & fieldMinMax = [0.0,0.8] ; gmem only
  ;colorField = 'timeRatio' & fieldMinMax = [0.0,4.0] ; gmem only

  ; pre-make cutouts (multiple or single)
  ;if keyword_set(showGmem) then cutout = cosmoVisCutout(sP=sP,gcInd=gcIDs,/selectGmem)
  ;if keyword_set(showAll)  then cutout = cosmoVisCutout(sP=sP,gcInd=gcIDs,sizeFac=sizeFac)

  ; target list
  gc    = loadGroupCat(sP=sP,/skipIDs)
  sgcen = subgroupPosByMostBoundID(sP=sP) 
  
  print,'rendering...'
  ; loop over all requested halos and image
  foreach gcID, gcIDs do begin
    
    ; load cutout
    if keyword_set(showGmem) then cutout = cosmoVisCutout(sP=sP,gcInd=gcID,/selectGmem)
    if keyword_set(showAll)  then cutout = cosmoVisCutout(sP=sP,gcInd=gcID,sizeFac=sizeFac)

    ; create color index mapping
    if colorField eq 'temp'      then fieldVal = cutout.loc_temp
    if colorField eq 'entropy'   then fieldVal = cutout.loc_ent
    if colorField eq 'density'   then fieldVal = cutout.loc_dens
    if colorField eq 'metal'     then fieldVal = cutout.loc_metal
    if colorField eq 'vrad'      then fieldVal = reform(cutout.loc_vrad)
    if colorField eq 'vradnorm'  then fieldVal = reform(cutout.loc_vrad)/cutout.haloV200
    if colorField eq 'coolTime'  then fieldVal = cutout.loc_coolTime
    if colorField eq 'dynTime'   then fieldVal = cutout.loc_dynTime
    if colorField eq 'timeRatio' then fieldVal = cutout.loc_coolTime / cutout.loc_dynTime
    
    colorinds = (fieldVal-fieldMinMax[0])*(255.0-nbottom) / (fieldMinMax[1]-fieldMinMax[0])
    colorinds = fix(colorinds + nbottom) > 0 < 255 ;nbottom-255  
  
    ; second/right panel cutout
    wSecond = where(fieldVal le secondCutVal,nCutoutSecond,comp=wComp)
      
    ; show gas above this cut value (instead of below)?
    if secondGt then wSecond = wComp
    
    loc_pos_second   = cutout.loc_pos[*,wSecond]
    loc_pos2_second  = cutout.loc_pos2[*,wSecond]
    colorinds_second = colorinds[wSecond]
    
    ; use instead a differently scaled color mapping for the second panel?
    if singleColorScale eq 0 then begin
      if secondGt eq 1 then $
        colorinds_second = (fieldVal-secondCutVal)*(255.0-nbottom) / (fieldMinMax[1]-secondCutVal)
      if secondGt eq 0 then $
        colorinds_second = (fieldVal-fieldMinMax[0])*(255.0-nbottom) / (secondCutVal-fieldMinMax[0])
      
      colorinds_second = fix(colorinds_second[wSecond] + nbottom) > 0 < 255 ; nbottom-255
    endif

    ; make a plot for each requested projection direction
    foreach axisPair, axes do begin
           
      ; get box center (in terms of specified axes)
      boxCenImg  = [sgcen[axisPair[0],gcID],sgcen[axisPair[1],gcID],sgcen[3-axisPair[0]-axisPair[1],gcID]]
      
      print,'['+string(gcID,format='(i4)')+'] Mapping ['+str(axisPair[0])+' '+$
            str(axisPair[1])+'] with '+str(cutout.boxSizeImg[0])+$
            ' kpc box around center ['+str(boxCenImg[0])+' '+str(boxCenImg[1])+' '+str(boxCenImg[2])+']'

      pFilename = 'scatter.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+$
                     '.axes'+str(axisPair[0])+str(axisPair[1])+'-'+$
                     colorField+'-'+secondField+'-'+str(secondGt)+'sCS'+str(singleColorScale)+'.eps'

      config = {boxSizeImg:cutout.boxSizeImg,plotFilename:pFilename,haloVirRad:cutout.haloVirRad,$
                haloMass:cutout.haloMass,axisPair:axisPair,sP:sP,$
                colorField:colorField,fieldMinMax:fieldMinMax,secondCutVal:secondCutVal,$
                secondGt:secondGt,nbottom:nbottom,barMM:fieldMinMax,ctName:ctName,barType:'2bar'}
      
      ; plot
      plotScatterComp,cutout.loc_pos,cutout.loc_pos2,loc_pos_second,loc_pos2_second,$
        colorinds,colorinds_second,config=config            

    endforeach ;axisPair

  endforeach ;gcIDs
  
  stop

end

; scatterMap4Panels(): four slices of some field for one halo
; showGmem/showAll (use galaxy catalog and display gmem only, or all in spatial subset)

pro scatterMap4Panels, showGmem=showGmem, showAll=showAll

  compile_opt idl2, hidden, strictarr, strictarrsubs

  sP = simParams(res=512,run='tracer',redshift=2.0)
  
  haloID = 304 ;z2.304 z2.301 z2.130 z2.64
  gcID = getMatchedIDs(sPa=sP,sPg=sP,haloID=haloID)
  gcIDs = [gcID.a]

  compile_opt idl2, hidden, strictarr, strictarrsubs

  if ~keyword_set(gcIDs) then message,'Error: Must specify gcIDs.'

  ; config
  singleColorScale = 1 ; 1=use same color scale for all panels, 0=rescale
  secondGt         = 0 ; 1=show greater than cut, 0=show less than cut
  axes             = list([0,1]) ;xy,xz,yz
  sizeFac          = 1.5 ; times rvir
  
  ; use which field and minmax for color mapping?
  ; temp [log K]
  ;colorField  = 'temp'
  ;fieldMinMax = [4.0,7.0]
  ;fieldRanges = list([4.5,5.0],[5.0,5.5],[5.5,6.0],[6.0,7.0])
  ;subtitles   = ['< 5.0','5.0 - 5.5','5.5 - 6.0','> 6.0']
  ;nbottom     = 50
  ;ctName      = 'helix' 
  
  ; vrad [km/s]
  colorField  = 'vrad'
  fieldMinMax = [-400,400]
  fieldRanges = list([-400,400],[-50,50],[-400,-300],[100,400])
  subtitles   = ['all','zero','inflow','outflow']
  nbottom     = 0
  ctName      = 'brewer-redblue'
  
  ; vrad norm [vrad/v200]
  ;colorField  = 'vradnorm'
  ;fieldMinMax = [0.0,4.0]
  ;fieldRanges = list([-5,3],[-1,1],[-5,-2],[2,3])
  ;subtitles   = ['all','zero','inflow,'outflow']
  ;nbottom     = 0
  ;ctName      = 'brewer-redblue'
  
  ; make cutouts (multiple or single)
  ;if keyword_set(showAll)  then cutout = cosmoVisCutout(sP=sP,gcInd=gcIDs,/selectGmem)
  ;if keyword_set(showAll)  then cutout = cosmoVisCutout(sP=sP,gcInd=gcIDs,sizeFac=sizeFac)

  ; target list
  gc    = loadGroupCat(sP=sP,/skipIDs)
  sgcen = subgroupPosByMostBoundID(sP=sP) 
  
  print,'rendering...'
  ; loop over all requested halos and image
  foreach gcID, gcIDs do begin
  
    ; load cutout
    if keyword_set(showAll)  then cutout = cosmoVisCutout(sP=sP,gcInd=gcID,/selectGmem)
    if keyword_set(showAll)  then cutout = cosmoVisCutout(sP=sP,gcInd=gcID,sizeFac=sizeFac)

    ; create color index mapping
    if colorField eq 'temp'      then fieldVal = cutout.loc_temp
    if colorField eq 'vrad'      then fieldVal = reform(cutout.loc_vrad)
    if colorField eq 'vradnorm'  then fieldVal = reform(cutout.loc_vrad)/cutout.haloV200
    
    ; make a plot for each requested projection direction
    foreach axisPair, axes do begin    
    
      pFilename = 'scatter4.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.h'+str(gcID)+$
                     '.axes'+str(axisPair[0])+str(axisPair[1])+'-'+$
                     colorField+'-'+str(secondGt)+'sCS'+str(singleColorScale)+'.eps'
      
      config = {boxSizeImg:cutout.boxSizeImg*sizeFac,plotFilename:pFilename,haloVirRad:cutout.haloVirRad,$
                haloMass:cutout.haloMass,axisPair:axisPair,sP:sP,barMM:fieldMinMax,$
                colorField:colorField,fieldMinMax:fieldMinMax,$
                ctName:ctName,nbottom:nbottom,barType:'1bar'}

      if ~singleColorScale then config.barType = ''
                
      start_PS, sP.plotPath + config.plotFilename, xs=8, ys=8
                
      ; cutouts and plot (k=0 top, k=1 bottom)
      for k=0,1 do begin
        wCutout = where(fieldVal ge fieldRanges[k*2+0,0] and fieldVal lt fieldRanges[k*2+0,1],nCutout)
          
        loc_pos_left = cutout.loc_pos[*,wCutout]
        loc_pos2_left = cutout.loc_pos2[*,wCutout]
        
        if singleColorScale then $
          colorinds = (fieldVal[wCutout]-fieldMinMax[0])*(255.0-nbottom) / $
                      (fieldMinMax[1]-fieldMinMax[0])
        if ~singleColorScale then $
          colorinds = (fieldVal[wCutout]-fieldRanges[k*2+0,0])*(255.0-nbottom) / $
                      (fieldRanges[k*2+0,1]-fieldRanges[k*2+0,0])
          
        colorinds_left = fix(colorinds + nbottom) > 0 < 255 ;nbottom-255  
        
        wCutout = where(fieldVal ge fieldRanges[k*2+1,0] and fieldVal lt fieldRanges[k*2+1,1],nCutout)
          
        loc_pos_right = cutout.loc_pos[*,wCutout]
        loc_pos2_right = cutout.loc_pos2[*,wCutout]
        
        if singleColorScale then $
          colorinds = (fieldVal[wCutout]-fieldMinMax[0])*(255.0-nbottom) / $
                      (fieldMinMax[1]-fieldMinMax[0])
        if ~singleColorScale then $
          colorinds = (fieldVal[wCutout]-fieldRanges[k*2+1,0])*(255.0-nbottom) / $
                      (fieldRanges[k*2+1,1]-fieldRanges[k*2+1,0])
          
        colorinds_right = fix(colorinds + nbottom) > 0 < 255 ;nbottom-255
        
        
        plotScatterComp,loc_pos_left,loc_pos2_left,loc_pos_right,loc_pos2_right,$
          colorinds_left,colorinds_right,config=config,$
          top=(k eq 0),bottom=(k eq 1),subtitle=subtitles[(k*2+0):(k*2+1)]
        
      endfor
      
      end_PS, pngResize=60  
  
    endforeach ;axisPair  
  
  endforeach ;gcIDs  
  stop
end