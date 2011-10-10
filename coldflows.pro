; coldflows.pro
; cold flows - main
; dnelson oct.2011

@helper
@coldflowsUtil
@coldflowsLoad
@coldflowsVis

; getSubgroupIDList(): get list of gas particle ids belonging to subgroups at specified snapshot
;                      endingSnap (including or excluding background subhalos) used as an
;                      exclusion list for deciding accreted gas
;
; bSH=1:  allow accretion from the background subhalo of each FOF group as not
;         part of any 'merger' related accretion and thus in the 'smooth' mode

function getSubgroupIDList, res=res, endingSnap=endingSnap, bSH=bSH

  if not keyword_set(res) or not keyword_set(endingSnap) then begin
     print,'Error: getSubgroupIDList() arguments not specified.'
     return,0
  endif

  gadgetPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'

  ; load subfind cat
  sgEnd = loadSubhaloGroups(gadgetPath,endingSnap)
  sgEnd_subgroupIDs = sgEnd.subgroupIDs
  
  if (keyword_set(bSH)) then begin
    ; make new sgEnd_subgroupIDs list specifically including only particles from primary
    ; (non-background) subhalos, since this is the exclusion list
    valSGids_end  = getPrimarySubhaloList(sgEnd)
    sgEnd_subgroupIDs = []
    
    foreach sgID, valSGids_end do begin
    
      ; select subgroup  
      nPart = sgEnd.subGroupLen[sgID]
      sgIDs = sgEnd.subGroupIDs[sgEnd.subGroupOffset[sgID] : sgEnd.subGroupOffset[sgID] + $
                                sgEnd.subGroupLen[sgID] - 1]
      
      sgEnd_subgroupIDs = [sgEnd_subgroupIDs, sgIDs]
    endforeach
    
  endif

  return, sgEnd_subgroupIDs

end

; findAccretedGas(): find gas particle ids from subhalos in targetSnap that were either:
;                    (a) not members of any subhalo at one specified earlier snapshot (smoothCutSnap)
;                    (b) not members of any subhalo for all previous snapshots (smoothCutSnap not set)
;                    [background subhalos optionally included or not through bSH flag]

function findAccretedGas, res=res, bSH=bSH, targetSnap=targetSnap, sgTarget=sgTarget, $
                          smoothCutSnap=smoothCutSnap

  if not keyword_set(res) or not keyword_set(targetSnap) or not keyword_set(sgTarget) then begin
     print,'Error: findAccretedGas() arguments not specified.'
     return,0
  endif
  
  ; config
  gadgetPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'  
  workingPath  = '/n/home07/dnelson/coldflows/thermhist.deriv/'  
  
  ; set range for enforcing smooth accretion
  if keyword_set(smoothCutSnap) then begin
      smoothStart = smoothCutSnap
      smoothEnd   = smoothCutSnap
  endif else begin
      smoothStart = targetSnap - 1
      smoothEnd   = 50 ; first group catalogs at z=6
  endelse
  
  ; set saveFilename and check existence
  saveFilename = workingPath + 'accreted.gas.'+str(res)+'.target='+str(targetSnap)+'.'+$
                  str(smoothStart)+'-'+str(smoothEnd)+'.sav'
                  
  if (keyword_set(bSH)) then $
    saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.bSH.sav'
                  
  if (file_test(saveFilename)) then begin
    restore,saveFilename
    return,sgIDs_All
  endif

  ; make list of sgIDs excluding background (id=0) subgroups
  valSGids = getPrimarySubhaloList(sgTarget)  
  
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
    nPart = sgTarget.subGroupLen[sgID]
    sgIDs = sgTarget.subGroupIDs[sgTarget.subGroupOffset[sgID] : sgTarget.subGroupOffset[sgID] + $
                                                                 sgTarget.subGroupLen[sgID] - 1]
    
    sgIDs_All = [sgIDs_All, sgIDs]
    countSGIDs_All += n_elements(sgIDs)
  endforeach  
  
  ; restrict to only gas
  match, gas_ids, sgIDs_All, gas_ids_ind, sgIDs_ind, count=count_gas
  
  if (count_gas gt 0) then begin
      print,'Found ['+str(n_elements(sgIDs_All))+'] particles (all subhalos)'+$
            ' after gas particle cut have ['+str(count_gas)+'] left, lost '+str(n_elements(sgIDs_All)-count_gas)+'.'
            
      ; keep only sgIDs[sgIDs_ind]
      sgIDs_All = sgIDs_All[sgIDs_ind]
  endif
  
  ; enforce smooth accretion
  print,'Running smooth accretion cut over snapshots ['+str(smoothStart)+'-'+str(smoothEnd)+'].'
  for m=smoothStart,smoothEnd,-1 do begin
    ; make (exclusion) list of particles belonging to subhalos
    sgEnd_subgroupIDs = getSubgroupIDList(res=res,endingSnap=m, bSH=bSH)  
    
    match, sgEnd_subgroupIDs, sgIDs_All, snap_ind, sgids_loc_ind, count=count_sm
    
    if (count_sm gt 0) then begin
      ; remove sgIDs[sgids_loc_ind] using complement
      all = bytarr(n_elements(sgIDs_All))
      if (sgids_loc_ind[0] ne -1L) then all[sgids_loc_ind] = 1B
      w = where(all eq 0B, ncomp)
  
      if (ncomp ne n_elements(sgIDs_All)-count_sm) then begin
        print,'ERROR',ncomp,n_elements(sgIDs_All),count_sm
        return,0
      endif
      
      sgIDs_All = sgIDs_All[w]
    endif
  endfor ;m
    
  print,'After smooth accretion cut have ['+str(ncomp)+'] left, lost '+str(count_gas-ncomp)+'.'
  
  save,sgIDs_All,filename=saveFilename
  
  return, sgIDs_All
  
end

; saveThermalHistories(): record the density+temperature for all gas particles at all snapshots

pro saveThermalHistories, res=res
  
  if not keyword_set(res) then begin
    print,'Error: Must specify resolution.'
    return
  endif
  
  ; config
  gadgetPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'
  workingPath  = '/n/hernquistfs1/dnelson/coldflows/thermhist/'
  saveFilebase = workingPath + 'thermhist.gas.'+str(res)+'_'
  
  nSnaps = n_elements(file_search(gadgetPath+'snapdir_*')) - 1 ;exclude final duplicate
  
  ; load first snapshot header
  h = loadSnapshotHeader(gadgetPath,snapNum=0)
  nGas = h.nPartTot[0]
  
  print,'Saving ['+str(nGas)+'] gas particles (rho,T) over ['+str(nSnaps)+'] snapshots.'
  
  ; loop from target snapshot through all prior snapshots
  for m=0,nSnaps-1,1 do begin
  
    ; check we haven't already done this
    saveFilename=saveFilebase + str(m) + '.sav'
    
    if (file_test(saveFilename)) then begin
      print,' Skipping: ' + saveFilename
      continue
    endif
  
    ; arrays
    temp    = fltarr(nGas+1)
    density = fltarr(nGas+1)
    
    ; setup such that first array index = particle ID, so index=0 is undefined
    temp[0]    = -1.0
    density[0] = -1.0
  
    if (m mod 10 eq 0) then print, ' '+str(m)
    ; load u and nelec and calculate temperature
    u     = loadSnapshotSubset(gadgetPath,snapNum=m,partType='gas',field='u')
    nelec = loadSnapshotSubset(gadgetPath,snapNum=m,partType='gas',field='ne')
    
    t = convertUtoTemp(u,nelec)
    
    u     = !NULL
    nelec = !NULL
    
    ; load density
    rho   = loadSnapshotSubset(gadgetPath,snapNum=m,partType='gas',field='density')
    
    ; load ids for sorted insert
    ids = loadSnapshotSubset(gadgetPath,snapNum=m,partType='gas',field='ids')
    
    ; save properties
    temp[ids]    = t
    density[ids] = rho

    rho = !NULL
    t   = !NULL
    ids = !NULL
    
    ; save
    save,temp,density,nGas,nSnaps,filename=saveFileName
    
  endfor

end

; maxTemperatures(): find maximum temperature for all gas particles in redshift range [zMin,zMax]

function maxTemperatures, res=res, bSH=bSH, zMin=zMin, zMax=zMax, getSnap=getSnap

  if not keyword_set(zMin) then zMin = 0.0
  if not keyword_set(zMax) then zMax = 30.0
  
  ; config
  dataPath     = '/n/hernquistfs1/dnelson/coldflows/thermhist/'
  workingPath  = '/n/home07/dnelson/coldflows/thermhist.deriv/'

  ; restrict temp array to redshift range
  maxSnap = redShiftToSnapnum(zMin)
  minSnap = redShiftToSnapnum(zMax)

  ; set saveFilenames and check for existence
  saveFilename1 = workingPath + 'maxtemp.'+str(res)+'.'+str(minSnap)+'-'+str(maxSnap)+'.sav'
  saveFilename2 = workingPath + 'maxtempsnap.'+str(res)+'.'+str(minSnap)+'-'+str(maxSnap)+'.sav'
  
  if (keyword_set(getSnap) and file_test(saveFilename2)) then begin
    restore, saveFilename2
    return, maxTempSnap
  endif
  if (not keyword_set(getSnap) and file_test(saveFilename1)) then begin
    restore, saveFilename1
    return, maxTemps
  endif
  
  print,'Calculating new maxtemp for res = '+str(res)+' in range ['+str(minSnap)+'-'+str(maxSnap)+'].'
  
  ; load first for nGas and make array
  restore,dataPath + 'thermhist.gas.'+str(res)+'_0.sav'
  maxTemps     = fltarr(nGas+1)
  maxTempSnap  = intarr(nGas+1)
  maxTemps[0]    = -1.0
  maxTempSnap[0] = -1.0
  
  ; load thermal history
  for m=minSnap,maxSnap,1 do begin
    thFilename = dataPath + 'thermhist.gas.'+str(res)+'_'+str(m)+'.sav'
    
    restore,thFilename
    
    ; take log
    w = where(temp eq 0,count)
    temp[w] = 1.0
    
    temp = alog10(temp)
    
    ; calc max temp
    w = where(temp gt maxTemps,count)
    if (count gt 0) then begin
      maxTemps[w]    = temp[w]
      maxTempSnap[w] = m
    endif
  
  endfor ;m

  ; save for future lookups
  if (not file_test(saveFilename1)) then $
    save,maxTemps,filename=saveFilename1
  if (not file_test(saveFilename2)) then $
    save,maxTempSnap,filename=saveFilename2

  if (keyword_set(getSnap)) then $
    return,maxTempSnap
  if (not keyword_set(getSnap)) then $
    return,maxTemps

end

; calcNormScalarProd(): calculate normalized scalar products of velocity vectors of individual
;                       gas particles prior to accretion, towards the subhalo center

pro calcNormScalarProd, res=res

  ; config
  gadgetPath   = '/n/hernquistfs1/mvogelsberger/ComparisonProject/'+str(res)+'_20Mpc/Gadget/output/'
  workingPath  = '/n/home07/dnelson/coldflows/thermhist.deriv/'
  
  bSH = 0 ; include background subhalos
  
  critLogTemp = 5.5

  redshifts = [3.0,2.0,1.0,0.0]
  targetSnaps = redshiftToSnapNum(redshifts)
  endingSnaps = targetSnaps - 1 ;immediate preceeding

  for m=0,n_elements(targetSnaps)-1 do begin
    targetSnap = targetSnaps[m]
    endingSnap = endingSnaps[m]

    ; set saveFilename
    saveFilename = workingPath + 'normscalar.'+str(res)+'.'+str(targetSnap)+'-'+str(endingSnap)+'.sav'
    
    if (keyword_set(bSH)) then $
      saveFilename = strmid(saveFilename,0,strlen(saveFilename)-4) + '.bSH.sav'
  
    if (file_test(saveFilename)) then begin
      print,'Skipping: '+saveFilename
      continue
    endif
  
    ; load subhalo group catalogs
    sgTarget = loadSubhaloGroups(gadgetPath,targetSnap) 
    
    ; get list of smoothly accreted gas particle ids
    sgIDs_Acc = findAccretedGas(res=res,bSH=bSH,targetSnap=targetSnap,sgTarget=sgTarget)
  
    ; get maximum temperatures
    maxTemps = maxTemperatures(res=res,bSH=bSH,zMin=redshifts[m])
    maxTemps = maxTemps[sgIDs_Acc]
  
    ; arrays
    normProd_Cold = list()
    normProd_Hot  = list()
    
    parentIDs      = intarr(n_elements(sgIDs_Acc))
    radialVectors  = fltarr(3,n_elements(sgIDs_Acc))
    
    ; load positions and particle ids
    ids = loadSnapshotSubset(gadgetPath,snapNum=endingSnap,partType='gas',field='ids')
    pos = loadSnapshotSubset(gadgetPath,snapNum=endingSnap,partType='gas',field='pos')
    
    ; loop over each accreted gas particle
    foreach pid,sgIDs_Acc,i do begin
      ; find parent subhalo ID
      pid_ind = where(sgTarget.subgroupIDs eq pid, count)
      
      if (count ne 1) then begin
        print,'ERROR'
        return
      endif
      
      parentID = where(pid_ind[0] - float(sgTarget.subgroupOffset) gt 0, count)
      
      if (count gt 0) then begin
        parentIDs[i] = parentID[count-1]
      endif else begin
        parentIDs[i] = 0
        print,'Warning: ParentID=0'
      endelse
      
      ; find parent subhalo (x,y,z) position
      parentCen = sgTarget.subgroupPos[*,parentIDs[i]]
  
      ; save radial vector
      w = where(ids eq pid,count)
      
      if (count ne 1) then begin
        print,'ERROR2'
        return
      endif
      
      pid_ind = w[0]
      radialVectors[*,i] = pos[*,pid_ind] - parentCen
      
    endforeach
  
    ; find unique parent IDs
    uniqParentIDs = parentIDs[uniq(parentIDs,sort(parentIDs))]
    
    print,'Processing [' + str(n_elements(uniqParentIDs)) + '] unique parent subhalos:'
  
    ; loop over each unique parent ID
    foreach parentID,uniqParentIDs,k do begin
      if (k mod 100 eq 0) then print,' '+str(k)
      ; find all COLD accreted gas particles of this parent subhalo
      w = where(parentIDs eq parentID and maxTemps lt critLogTemp,count)
      
      if (count gt 1) then begin
        ; COLD: pairwise (double loop) compute: normalized scalar product and save
        for i=0,count-1 do begin
          for j=i+1,count-1 do begin
            dotp = radialVectors[*,i] ## transpose(radialVectors[*,j])
            norm = sqrt(radialVectors[*,i] ## transpose(radialVectors[*,i])) * $
                   sqrt(radialVectors[*,j] ## transpose(radialVectors[*,j]))
            normProd_Cold->add, (dotp/norm)
          endfor
        endfor
      endif
      
      ; find all HOT accreted gas particles of this parent subhalo
      w = where(parentIDs eq parentID and maxTemps ge critLogTemp,count)
      
      if (count gt 1) then begin
        ; HOT: pairwise product
        for i=0,count-1 do begin
          for j=i+1,count-1 do begin
            dotp = radialVectors[*,i] ## transpose(radialVectors[*,j])
            norm = sqrt(radialVectors[*,i] ## transpose(radialVectors[*,i])) * $
                   sqrt(radialVectors[*,j] ## transpose(radialVectors[*,j]))
            normProd_Hot->add, (dotp/norm)
          endfor
        endfor
      endif
      
    endforeach
    
    ; save results
    normProd_Hot  = normProd_Hot.toArray()
    normProd_Cold = normProd_Cold.toArray()
    save,sgIDs_Acc,normProd_Hot,normProd_Cold,parentIDs,$
         radialVectors,targetSnap,endingSnap,filename=saveFilename
    
    print,'Saved: '+saveFilename
  endfor ; targetSnaps
  
end

; selectFilament(): beginning of an algorithm to select a filamentary structure incident on a 
;                   specified subhalo [x,y,z] position based on a column-like spatial structure
;                   for all properties falling inside selection, plot various quantities

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

@coldflowsPlot
