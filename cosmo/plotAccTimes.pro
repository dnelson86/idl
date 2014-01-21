; plotAccTimes.pro
; gas accretion project - past radial history of gas elements (plotting)
; dnelson jan.2014

; plotAccTimeDeltas():

pro plotAccTimeDeltas

  ; config
  redshifts = [3.0,2.0,1.0,0.0]
  runs      = ['feedback','tracer']
  res       = 128
  
  ; index selection for difference ([1.0,0.75,0.5,0.25,0.15,0.05,0.01,first0.15,first1.0])
  earlierInd = -1 
  laterInd   = -2

  ; plot config
  yrange     = [0.0,1.0]
  xrange     = [0.0,5.0] ;[-5.0,5.0]
  xrangeNorm = [-1.0,1.5]
  nBins      = 20
  colorInd   = 1
  
  sP = simParams(res=res,run=runs[0],redshift=redshifts[0])
  
  start_PS,sP.plotPath+'accTimeDeltas_'+str(res)+'.eps', xs=3*3.5, ys=4*3.5
  
    pos = plot_pos(rows=4,cols=2,/gap)
    colors = []
    
    foreach redshift,redshifts,i do begin
    
    xrange = [-1,1] * redshiftToAgeFlat(redshift)
    
    ; plot (1)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,title="",$
      xtitle=textoidl("( t_{1} - t_{2} ) [Gyr]"),$
      ytitle="Fraction",/xs,/ys,pos=pos[i*2],/noerase
	
	; plot histogram for each run/redshift combination
      foreach run,runs,j do begin
        print,redshift,run
        sP_cur = simParams(res=res,run=run,redshift=redshift)
        at     = accretionTimes(sP=sP_cur)
        
        w = where(at.accTime[earlierInd,*] ge 0 and at.accTime[laterInd,*] ge 0,count)
        if count eq 0 then message,'Error'
        
        age_earlier = reform( redshiftToAgeFlat( 1.0/at.accTime[earlierInd,w]-1.0 ) )
        age_later   = reform( redshiftToAgeFlat( 1.0/at.accTime[laterInd,w]-1.0 ) )
        
        val = age_later - age_earlier

        if min(val) lt xrange[0] or max(val) gt xrange[1] then begin
          print,'Warning: val min = ',min(val),' max = ',max(val)
        endif
        
        binSize = (xrange[1] - xrange[0]) / float(nBins)
	  h = histogram(val, bin=binSize, loc=loc, min=xrange[0], max=xrange[1])
	  h = float(h)/total(h)
	  
        colors = [colors,sP_cur.colors[colorInd]]
	  cgPlot,loc+binSize*0.5,h,color=colors[-1],/overplot
	endforeach
	
	; legend
	legend,['z = '+string(redshift,format='(f3.1)'),runs],textcolors=[0L,colors],/top,/right
	
    ; plot (2)
    cgPlot,[0],[0],/nodata,xrange=xrangeNorm,yrange=yrange,title="",$
      xtitle=textoidl("( t_{1} - t_{2} ) / t_{1}"),$
      ytitle="Fraction",/xs,/ys,pos=pos[i*2+1],/noerase
	
	; plot histogram for each run/redshift combination
      foreach run,runs,j do begin
        print,redshift,run
        sP_cur = simParams(res=res,run=run,redshift=redshift)
        at     = accretionTimes(sP=sP_cur)
        
        w = where(at.accTime[earlierInd,*] ge 0 and at.accTime[laterInd,*] ge 0,count)
        if count eq 0 then message,'Error'
        
        age_earlier = reform( redshiftToAgeFlat( 1.0/at.accTime[earlierInd,w]-1.0 ) )
        age_later   = reform( redshiftToAgeFlat( 1.0/at.accTime[laterInd,w]-1.0 ) )
        
        val = age_later - age_earlier
        val /= age_earlier
        
        if min(val) lt xrange[0] or max(val) gt xrange[1] then begin
          print,'Warning: norm val min = ',min(val),' max = ',max(val)
        endif
        
        binSize = (xrange[1] - xrange[0]) / float(nBins)
	  h = histogram(val, bin=binSize, loc=loc, min=xrange[0], max=xrange[1])
	  h = float(h)/total(h)
        
	  cgPlot,loc+binSize*0.5,h,color=colors[j],/overplot
	endforeach
	
	; legend
	legend,['z = '+string(redshift,format='(f3.1)'),runs],textcolors=[0L,colors],/top,/right
	
    endforeach ;redshifts
      
  end_PS
  stop
end

; shyPlot(): at a target redshift, make a selection of galaxy star tracers which were in gas cells
;            in the previous timestep (newly formed stars). normalize the time elapsed since the 
;            first 1.0rvir crossing of each such tracer by the age of the universe at the time of
;            that crossing. plot the resulting normalized 'accTime' distribution, comparing
;            TRACER to FEEDBACK at z=0,1,2,3

function shyPlotBin, sP=sP, aMode=aMode, allTypes=allTypes

  if keyword_set(allTypes) then allTypesStr = '.allTypes' else allTypesStr = ''

  ; set saveFilename and check for existence
  saveFilename = sP.derivPath + 'shyPlot.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(sP.snap)+'.mode'+str(aMode)+allTypesStr+'.sav'
  
  if file_test(saveFilename) then begin
    restore, saveFilename
    return, r
  endif

  nSnapsBack = 1
  
  ; load at target redshift
  galcat = galaxyCat(sP=sP)
  
  ; make tracer selection (children of stars in galaxies only)
  if keyword_set(allTypes) then begin
    print,'AllTypes = 1'
    
    ; load all star ids
    star_ids = loadSnapshotSubset(sP=sP, partType='star', field='ids')
    
    ; replicate parent id list for all galcat tracers
    tr_parids = galcat.IDs[ replicate_var(galcat.trMC_cc) ]
    
    ; decide overlap, use those tracers (tr_parids not unique, use value_locate approach)
    sort_inds = calcSort(star_ids)
    star_ids_sorted = star_ids[sort_inds]
        
    par_ind = value_locate(star_ids_sorted,tr_parids) ; indices to star_ids_sorted
    par_ind = sort_inds[par_ind>0] ; indices to star_ids (>0 removes -1 entries, which are removed next line)
    tr_parids_ind = where(star_ids[par_ind] eq tr_parids,count_inPar) ; verify we actually matched the ID

    ; save master list of ids/inds to search for
    tr_inds_orig = par_ind[tr_parids_ind]
    tr_ids_orig  = galcat.trMC_ids[ tr_inds_orig ]
  endif else begin
    print,'AllTypes = 0'
    type = ( galcat.type[ replicate_var(galcat.trMC_cc) ] )
  
    tr_inds_orig = where(type eq galcat.types.stars)
    tr_ids_orig = galcat.trMC_ids[ tr_inds_orig ]
  endelse
  
  ; move back requested number of snapshots in a row, load tracer parent ids for this subset
  tr_ids_final  = []
  tr_inds_final = []

  for i=0,nSnapsBack-1 do begin
    sP.snap -= 1
    tr_ids = loadSnapshotSubset(sP=sP, partType='tracerMC', field='tracerids')
  
    idIndexMap = getIDIndexMap(tr_ids, minid=minid)
    tr_ids = !NULL
    tr_inds = idIndexMap[ tr_ids_orig - minid ]
  
    tr_parids = loadSnapshotSubset(sP=sP, partType='tracerMC', field='parentids')
    tr_parids = tr_parids[ tr_inds ]
  
    ; load gas ids and crossmatch
    gas_ids = loadSnapshotSubset(sP=sP, partType='gas', field='ids')
  
    sort_inds = calcSort(gas_ids)
    gas_ids_sorted = gas_ids[sort_inds]
    gas_ind = value_locate(gas_ids_sorted,tr_parids) ; indices to gas_ids_sorted
    gas_ind = sort_inds[gas_ind>0] ; indices to gas_ids (>0 removes -1 entries, which are removed next line)
    w = where(gas_ids[gas_ind] eq tr_parids,count_inPar) ; verify we actually matched the ID

    tr_ids_final  = [ tr_ids_final, tr_ids_orig[w] ]
    tr_inds_final = [ tr_inds_final, tr_inds_orig[w] ] ; index galcat.trMC_ids
    print,count_inPar
  endfor
  
  ; now just insure uniqueness
  if nSnapsBack gt 1 then begin
    tr_ids_final = tr_ids_final[sort(tr_ids_final)]
    tr_inds_final = tr_inds_final[sort(tr_inds_final)]
    
    w = where( tr_ids_final ne shift(tr_ids_final,1), count)
    if count eq 0 then message,'Error, expect some duplicates.'
    tr_ids_final = tr_ids_final[w]
    
    w = where( tr_inds_final ne shift(tr_inds_final,1), count)
    if count eq 0 then message,'Error, expect some duplicates.'
    tr_inds_final = tr_inds_final[w]
  endif
  
  print,' Found: ['+str(n_elements(tr_ids_final))+'] of ['+str(n_elements(tr_ids_orig))+'] tracers.'

  gas_ids         = !NULL
  star_ids        = !NULL
  gas_ind         = !NULL
  gas_ids_sorted  = !NULL
  star_ids_sorted = !NULL
  tr_parids_ind   = !NULL
  par_ind         = !NULL
  tr_parids       = !NULL
  tr_inds         = !NULL
  sort_inds       = !NULL
  w               = !NULL
  
  ; return to original snapshot
  sP.snap += nSnapsBack
  
  ; aMode=1: load accTimes of first 1.0rvir crossing of main progenitor branch
  if aMode eq 1 then begin
    at = accretionTimes(sP=sP)
    at = reform( at.accTime[-1,tr_inds_final] )
  
    ; restrict to accretion mode (all are +_rec) and valid accTimes
    ;am = accretionMode(sP=sP)
    ;am = reform( am[tr_inds_final] )
    ;if accMode eq 'smooth'   then $
    ;  w_at = where(at ne -1 and (am eq 1 or am eq 11),count)
    ;if accMode eq 'clumpy'   then $
    ;  w_at = where(at ne -1 and (am eq 2 or am eq 3 or am eq 12 or am eq 13),count)
    ;if accMode eq 'stripped' then $
    ;  w_at = where(at ne -1 and (am eq 4 or am eq 14),count)
    ;if accMode eq 'all' then $
    ;  w_at = where(at ne -1,count)
    ;am = !NULL
  
    ; filter out -1 values
    at = at[where(at gt 0.0)]
    
    ; convert scalefac to redshift
    at = 1.0/at-1 ; replace at by at[w_at] to use accMode selection
  endif
  
  ; aMode=2: group membership, use mergerTreeSubset() of the subgroups, walk back through 
  ; the snapshots, update Parent (for each tracer) at each step, calculate current Parent (for each
  ; tracer) at each step, or -2 for not in group catalog, save first mismatch time as 
  ; the accretion time entering the subgroup
  if aMode eq 2 then begin
    mt = mergerTreeSubset(sP=sP)
    
    ; at each snapshot holds the tracked main progenitor parents, for each tracer
    curParInds = mt.gcIndOrigTrMC[ tr_inds_final ]
    
    ; the calculated accretion time
    at = fltarr( n_elements(curParInds) ) - 1.0
    
    origSnap   = sP.snap
    targetSnap = sP.groupCatRange[0]
    
    for m=origSnap,targetSnap+1,-1 do begin
      sP.snap = m
          
      ; at each snapshot holds the actual subhalo parent, for each tracer
      actualParInds = lonarr( n_elements(tr_inds_final) )
          
      ; load group catalog
      gcCur = loadGroupCat(sP=sP,/readIDs)
      
      ; calculate actualParInds
        ; load all tracer IDs and re-locate starting IDs
        tr_ids = loadSnapshotSubset(sP=sP, partType='tracerMC', field='tracerids')
        
        idIndexMap = getIDIndexMap(tr_ids, minid=minid)
        tr_ids = !NULL
        tr_inds = idIndexMap[ tr_ids_final - minid ] ; overwrite previous
        
        ; load parent IDs for this subset (don't care about type)
        tr_parids = loadSnapshotSubset(sP=sP, partType='tracerMC', field='parentids', inds=tr_inds)
        
        ; note: tr_parids are NOT UNIQUE, use a value_locate approach (not match)
          sort_inds = calcSort(gcCur.IDs)
          GC_ids_sorted = gcCur.IDs[sort_inds]
        
          par_ind = value_locate(GC_ids_sorted,tr_parids) ; indices to GC_ids_sorted
          par_ind = sort_inds[par_ind>0] ; indices to gc.IDs (>0 removes -1 entries, which are removed next line)
          tr_parids_ind = where(gcCur.IDs[par_ind] eq tr_parids,count_inPar) ; verify we actually matched the ID
          gcCur_ids_ind = par_ind[tr_parids_ind] ; indices of matched parents
          countMatch = n_elements(tr_parids_ind)
          
          ; for those tracers with non-matching parent IDs (not in group catalog), set to -2
          matchMask = bytarr( n_elements(tr_parids) )
          matchMask[tr_parids_ind] = 1B
          
          w_nonMatch = where(matchMask eq 0B,count)
          
          if count gt 0 then begin
            actualParInds[w_nonMatch] = -2
          endif
          
          ; for those tracers with matching parent IDs
          if countMatch gt 0 then begin
            ; do a value locate of the index (of gcCur.IDs) with group.offset
            ; use whole FOF!
            parCandidates = value_locate( gcCur.groupOffset, gcCur_ids_ind )

            ; calculate delta(index,offset) to check if < group.len, to verify actually in group
            parDelta = gcCur_ids_ind - gcCur.groupOffset[ parCandidates ]
            
            w = where(parDelta ge 0 and parDelta lt gcCur.groupLen[ parCandidates ],$
              countCand,comp=wc,ncomp=countOutsideCand)
            
            if countCand gt 0 then begin
              subgroupParInds = gcCur.groupFirstSub[ parCandidates[w] ]
              actualParInds[tr_parids_ind[w]] = subgroupParInds
            endif
            
            ; those outside their candidates (should not happen since considering full fof groups)
            if countOutsideCand gt 0 then message,'Error, reconsider'
          endif
      
      ; for now just mismatching parInds, set accretionTime
      w = where(curParInds ne actualParInds and curParInds ne -1,count)

      if max(at[w]) gt 0.0 then message,'Error: About to set at for tracers with already set at.'
      if count gt 0 and sP.snap eq origSnap then message,'Error: Bad parents at sP.snap'

      at[w] = snapNumToRedshift(sP=sP)
      curParInds[w] = -1 ; no longer search/set on these tracers
           
      frac = float(count)*100/n_elements(tr_ids_final)
      print,'['+str(m)+'] frac now set = '+string(frac,format='(f4.1)')+'% ('+$
        str(count)+' of '+str(n_elements(tr_ids_final))+')'

      ; load mergerTree and move to Parent
      Parent = mergerTree(sP=sP)
        
      ; change to parent IDs for the next snapshot
      w = where(curParInds ne -1,count)
      if count eq 0 then message,'error'
        
      curParInds[w] = Parent[curParInds[w]]
        
      frac = float(count)*100/n_elements(tr_ids_final)
      print,'['+str(m)+'] frac remaining = '+string(frac,format='(f4.1)')+'% ('+$
        str(count)+' of '+str(n_elements(tr_ids_final))+')'
          
    endfor ; snapshot
    
    ; restrict times to those we found (already in redshift)
    at = at[where(at gt 0.0)]
  endif
  
  if aMode ne 1 and aMode ne 2 then message,'Error'
  
  ; age of universe (in Gyr) at each accTime
  r = { age_at  : redshiftToAgeFlat(at) ,$
        age_cur : redshiftToAgeFlat(sP.redshift) }
        
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
        
  
  return, r
  
end

; shyPlot():
;  plot (1): delta time between "halo accretion time" and "star forming time" (Gyr)
;  plot (2): delta time as above, normalized by age of universe at "halo accretion time"
; aMode=1: "halo accretion time" = first 1.0 rvir crossing from accretionTimes()
; aMode=2: "halo accretion time" = group membership
; allTypes=1: all tracers in galaxyCat, otherwise just 'star' (e.g. inside 0.15rvir cut) tracers

pro shyPlot

  ; config
  redshifts = [3.0,2.0,1.0,0.0]
  runs      = ['feedback','tracer']
  res       = 512
  aMode     = 1
  allTypes  = 1
  
  ; plot config
  yrangeRaw  = [0.0,0.2] ;[0.0,0.1]
  yrangeNorm = [0.0,0.3]
  nBins     = 20
  colors    = ['red','blue'] ; one per run
  
  start_PS,'shyPlot_aMode'+str(aMode)+'_allTypes'+str(allTypes)+'_'+str(res)+'b.eps', xs=3*3.5, ys=4*3.5
  
    pos = plot_pos(rows=4,cols=2,/gap)
    
    foreach redshift,redshifts,i do begin
    
    ;xrange = [-0.1, redshiftToAgeFlat(redshift)*1.0]
    xrange = [-0.1, 2.0]
    
    ; plot (1)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeNorm,title="",$
      xtitle=textoidl("( t_{SF} - t_{acc} ) / t_{H}(t_{acc})"),$
      ytitle="Fraction",/xs,/ys,pos=pos[i*2],/noerase
	
	; plot histogram for each run/redshift combination
      foreach run,runs,j do begin
        sP_cur = simParams(res=res,run=run,redshift=redshift)
        bV     = shyPlotBin(sP=sP_cur, aMode=aMode, allTypes=allTypes)
        
        val = (bV.age_cur - bV.age_at) / bV.age_at

        if min(val) lt xrange[0] or max(val) gt xrange[1] then begin
          print,'Warning: normalized val min = ',min(val),' max = ',max(val)
        endif
        
        
        binSize = (xrange[1] - xrange[0]) / float(nBins)
	  h = histogram(val, bin=binSize, loc=loc, min=xrange[0], max=xrange[1])
	  h = float(h)/total(h)
	  
	  cgPlot,loc+binSize*0.5,h,color=cgColor(colors[j]),/overplot
	endforeach
	
	; legend
	legend,['z = '+string(redshift,format='(f3.1)'),runs],textcolors=['black',colors],/top,/right
	
    ; plot (2)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeRaw,title="",$
      xtitle=textoidl("( t_{SF} - t_{acc} ) [Gyr]"),$
      ytitle="Fraction",/xs,/ys,pos=pos[i*2+1],/noerase
	
	; plot histogram for each run/redshift combination
      foreach run,runs,j do begin
        sP_cur = simParams(res=res,run=run,redshift=redshift)
        bV     = shyPlotBin(sP=sP_cur, aMode=aMode, allTypes=allTypes)
        
        val = (bV.age_cur - bV.age_at)
        
        if min(val) lt xrange[0] or max(val) gt xrange[1] then begin
          print,'Warning: raw val min = ',min(val),' max = ',max(val)
        endif
        
        binSize = (xrange[1] - xrange[0]) / float(nBins)
	  h = histogram(val, bin=binSize, loc=loc, min=xrange[0], max=xrange[1])
	  h = float(h)/total(h)
	  
	  cgPlot,loc+binSize*0.5,h,color=cgColor(colors[j]),/overplot
	endforeach
	
	; legend
	legend,['z = '+string(redshift,format='(f3.1)'),runs],textcolors=['black',colors],/top,/right
	
    endforeach ;redshifts
      
  end_PS
  
  stop
  
end