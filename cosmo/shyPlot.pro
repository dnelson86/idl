; shyPlot.pro
; for shy's main illustris paper
; dnelson apr.2014

; galcatIndsParPartType(): return all trMC in galcat with parents of a particular particle type
;   e.g. handle mixed inter
    
function galcatIndsParPartType, sP=sP, galcat=galcat, partType=partType
  compile_opt idl2, hidden, strictarr, strictarrsubs
  ; load all parent type ids
  par_ids = loadSnapshotSubset(sP=sP, partType=partType, field='ids')
    
  ; replicate parent id list for all galcat tracers
  tr_parids = galcat.IDs[ replicate_var(galcat.trMC_cc) ]
    
  ; decide overlap, use those tracers (tr_parids not unique, use value_locate approach)
  sort_inds = calcSort(par_ids)
  par_ids_sorted = par_ids[sort_inds]
        
  par_ind = value_locate(par_ids_sorted,tr_parids) ; indices to par_ids_sorted
  par_ind = sort_inds[par_ind>0] ; indices to par_ids (>0 removes -1 entries, which are removed next line)
  tr_inds = where(par_ids[par_ind] eq tr_parids,count_inPar) ; verify we actually matched the ID

  ; master list of ids/inds to search for
  return, tr_inds
end

; shyPlotBin(): at a target redshift, make a selection of galaxy star tracers which were in gas cells
;            in the previous timestep (newly formed stars). normalize the time elapsed since the 
;            first 1.0rvir crossing of each such tracer by the age of the universe at the time of
;            that crossing. plot the resulting normalized 'accTime' distribution, comparing
;            TRACER to FEEDBACK at z=0,1,2,3

function shyPlotBin, sP=sP, aMode=aMode, allTypes=allTypes, blindSearch=blindSearch
  compile_opt idl2, hidden, strictarr, strictarrsubs
  if keyword_set(allTypes) then allTypesStr = '.allTypes' else allTypesStr = ''
  if aMode ne 1 and aMode ne 2 then message,'Error'
  if allTypes ne 0 and allTypes ne 1 then message,'Error'
  
  ; set saveFilename and check for existence
  saveFilename = sP.derivPath + 'shyPlot/shyPlot.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(sP.snap)+'.mode'+str(aMode)+allTypesStr+'.sav'
  
  addParMassFlag = 0
  
  if file_test(saveFilename) then begin
    restore, saveFilename
   
    ; return regardless
    return,r
 
    ; return only if parMass already is calculated, otherwise calculate and save it now
    if tag_exist(r, 'parMass') then return,r
    addParMassFlag = 1 & print,'Request parMass add: '+saveFilename
    
  endif else begin
    ; file does not exist, and we're just searching to see what is available, so return empty
    if keyword_set(blindSearch) then return, {flag:-1}
  endelse
  
  nSnapsBack = 1
  
  ; load at target redshift
  galcat = galaxyCat(sP=sP)
  
  ; make tracer selection (children of stars in galaxies only)
  if keyword_set(allTypes) then begin
    print,'AllTypes = 1'
    
    tr_inds_orig = galcatIndsParPartType(sP=sP, galcat=galcat, partType='star')
    tr_ids_orig = galcat.trMC_ids[ tr_inds_orig ]
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
    
    ; missing parMasses field?
    parMass = galCatParentProperties(sP=sP,/mass,/trRep)
    parMass = parMass[tr_inds_final[ where(at gt 0.0) ]]
    
    if addParMassFlag eq 1 then begin
      r = mod_struct( r, 'parMass', parMass )
      save,r,filename=saveFilename
      print,'Saved (added parMass): '+strmid(saveFilename,strlen(sp.derivPath))
      return, r
    endif
  
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

      if count gt 0 and max(at[w]) gt 0.0 then message,'Error: About to set at for tracers with already set at.'
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
    
    sP.snap = origSnap
    
    ; request to return parent masses of aMode selection?
    parMass = galCatParentProperties(sP=sP,/mass,/trRep)
    parMass = parMass[tr_inds_final[ where(at gt 0.0) ]]
    
    if addParMassFlag eq 1 then begin
      r = mod_struct( r, 'parMass', parMass )
      save,r,filename=saveFilename
      print,'Saved (added parMass): '+strmid(saveFilename,strlen(sp.derivPath))
      return, r
    endif
    
    ; restrict times to those we found (already in redshift)
    at = at[where(at gt 0.0)]
  endif
  
  ; age of universe (in Gyr) at each accTime
  r = { age_at  : redshiftToAgeFlat(at) ,$
        age_cur : redshiftToAgeFlat(sP.redshift),$
        parMass : parMass }
        
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
  
  return, r
  
end

; shyInvPlotAccretedTracers(): find tracers which traverse 1.0rvir between sP.snap and next
  
function shyInvPlotAccretedTracers, sP=sP, parMass=parMass
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  cutRadius = 1.0 ; times rvir
  
  ; load at target redshift
  galcat = galaxyCat(sP=sP, /skipSave)
  rad = galCatParentProperties(sP=sP,galcat=galcat,/rVirNorm,/trRep)
  
  ; initial tracer selection  
  tr_inds_orig = galcatIndsParPartType(sP=sP, galcat=galcat, partType='gas')
  w = where( rad[tr_inds_orig] ge cutRadius, count )
  if count eq 0 then message,'Error'

  tr_ids_orig = galcat.trMC_ids[ tr_inds_orig[w] ]
  
  ; increment snapshot by one and load (memory only)
  sP.snap += 1
  
  galcat = galaxyCat(sP=sP, /skipSave)
  rad = galCatParentProperties(sP=sP,galcat=galcat,/rVirNorm,/trRep)
  
  w = where( rad lt cutRadius, count )
  if count eq 0 then message,'Error'
  
  tr_ids_final = galcat.trMC_ids[ w ]
  
  ; get tracer ID selection satisfying 1.0rvir traversal criterion
  tr_ids_final = intersection( tr_ids_final, tr_ids_orig )
  
  ; parent halo masses requested?
  if keyword_set(parMass) then begin
    parMass = galCatParentProperties(sP=sP,galcat=galcat,/mass,/trRep)
    
    calcMatch,tr_ids_final,galcat.trMC_ids,ind1,ind2,count=countMatch
    ind2 = ind2[sort(ind1)]
    
    parMass = parMass[ind2]
  endif
  
  ; restore to original snap
  sP.snap -= 1
  
  return, tr_ids_final
end

; shyInvPlotBin(): for all tracers in gas parents inside galCat at target redshift, and which 
;   are outside 1.0 rvir but within 1.0 rvir in the immediately next snapshot (lower redshift), 
;   the fraction of these which reside in star parents at z=0

function shyInvPlotBin, sP=sP
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; config
  endingRedshift = 0.0 ; when to compare to star child tracers

  ; set saveFilename and check for existence
  saveFilename = sP.derivPath + 'shyPlot/shyInvPlot.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+str(sP.snap)+'.sav'

  if file_test(saveFilename) then begin
    restore, saveFilename
    return, r
  endif

  ; get accreted tracers
  parMass = 1
  tr_ids_final = shyInvPlotAccretedTracers(sP=sP,parMass=parMass)
  
  ; load z=0 all star_ids, get all tracer children
  sP.redshift = endingRedshift
  sP.snap     = redshiftToSnapNum(sP=sP)
  
  star_ids = loadSnapshotSubset(sP=sP, partType='star', field='ids')
  tr_children_stars = cosmoTracerChildren(sP=sP,/getIDs,starIDs=star_ids)
  
  ; intersect final candidates with redshift zero star parent tracers
  tr_ids_stars  = intersection( tr_ids_final, tr_children_stars )
  frac_to_stars = float( n_elements(tr_ids_stars) ) / n_elements(tr_ids_final)
  
  r = { tr_ids_stars:tr_ids_stars, frac_to_stars:frac_to_stars, parMass:parMass }
  
  save,r,filename=saveFilename
  print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
  
  return, r

end

; shyFullInvBin():

function shyFullInvBin, sP=sP, getResult=getResult, incWind=incWind
  compile_opt idl2, hidden, strictarr, strictarrsubs
  if ~keyword_set(sP) then message,'Error'
  
  ; config
  redshifts = [5.0,4.0,3.5,3.0,2.5,2.0,1.5,1.0,0.75,0.5,0.25]

  ; set target save points and restart filename
  if keyword_set(incWind) then windStr = '.wWind' else windStr = ''
  saveSnaps = redshiftToSnapNum( redshifts, sP=sP )

  restartFilename = sP.derivPath + 'shyPlot/shyFullInvRestart.'+sP.saveTag+'.'+$
                    sP.savPrefix+str(sP.res)+windStr+'.sav'
  
  ; result requested? return or throw error for lack of existence
  if keyword_set(getResult) then begin
    saveFilename = sP.derivPath + 'shyPlot/shyFullInvPlot.'+sP.saveTag+'.'+$
                   sP.savPrefix+str(sP.res)+'.' + str(sP.snap) + windStr + '.sav'
  
    if file_test(saveFilename) then begin
      restore, saveFilename
      return, r
    endif else begin
      print,'Skip: save not found ['+sP.run+' '+str(sP.snap)+']!.'
      return, []
    endelse
  endif
  
  snapRange = reverse( sP.groupCatRange )
  snapStep  = -1
  
  ; restart file exists? if so, load it now
  if file_test(restartFilename) then begin
    restore, restartFilename, /verbose
    snapRange[0] = m
  endif else begin
    ; allocate (lastStarTime[tr_ID-minid] indexed)
    tr_ids = loadSnapshotSubset(sP=sP,partType='tracerMC',field='tracerids')
    
    maxid = max(tr_ids)
    minid = min(tr_ids)
    print,'Allocate for min='+str(minid)+' max='+str(maxid)+'  '+windStr
    tr_ids = !NULL
    
    lastStarTime = fltarr( maxid-minid+1 ) - 1
  endelse
  
  for m=snapRange[0],snapRange[1],snapStep do begin
    sP.snap = m & print,m
    sP.redshift = snapNumToRedshift(sP=sP)
      
    ; save restart?
    if m mod 10 eq 0 and m lt snapRange[0] then begin
      print,' --- Writing restart! ---'
      save,lastStarTime,m,maxid,minid,filename=restartFilename
      print,' --- Done! ---'
    endif
      
    ; get all star children tracers
    star_ids = loadSnapshotSubset(sP=sP, partType='star', field='ids')
    
    ; restrict to real stars only?
    if keyword_set(noWind) then begin
      star_sft = loadSnapshotSubset(sP=sP, partType='stars', field='gfm_sftime')
      if n_elements(star_ids) ne n_elements(star_sft) then message,'Error'
      
      w = where(star_sft ge 0.0,count)
      if count eq 0 then message,'Error'
      star_ids = star_ids[w]
    endif
    
    tr_children_stars = cosmoTracerChildren(sP=sP,/getIDs,starIDs=star_ids)

    lastStarTime[ tr_children_stars-minid ] = 1.0/(1+sP.redshift)
    
    if min(tr_children_stars) lt minid or max(tr_children_stars) gt maxid then message,'Error'
    
    ; are we at a save point? if so, calculate the instantaneously accreting tracers
    if total( saveSnaps eq sP.snap ) gt 0 then begin
      print,'Save point reached: ['+str(sP.snap)+'] (z='+string(sP.redshift,format='(f4.1)')+')'
      
      saveFilename = sP.derivPath + 'shyPlot/shyFullInvPlot.'+sP.saveTag+'.'+$
                     sP.savPrefix+str(sP.res)+'.' + str(sP.snap) + windStr + '.sav'
                     
      ; get accreted tracers
      parMass = 1
      tr_ids_acc = shyInvPlotAccretedTracers( sP=sP, parMass=parMass )
      
      ; which of those have been already in stars by this point?
      w = where( lastStarTime[tr_ids_acc-minid] gt 0.0, count )
      if count eq 0 then message,'Error'
      
      ; calculate save values
      frac_to_stars = float(count) / n_elements(tr_ids_acc)
      age_sf        = redshiftToAgeFlat( 1.0/lastStarTime[tr_ids_acc[w]-minid]-1.0 )
      age_at        = redshiftToAgeFlat( sP.redshift )
      parMass       = parMass[w]
      
      r = { age_at:age_at, age_sf:age_sf, frac_to_stars:frac_to_stars, $
            parMass:parMass, tr_ids_acc:tr_ids_acc }

      save,r,filename=saveFilename
      print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    endif
  endfor
    
    
end

; shyInvPlot(): plot the fraction of tracers which have just -now- (instantaneously, from one snapshot to the next)
;   accreted onto the halo (i.e. traversed 1.0 rvir) in GAS parents, which at z=0 now reside in STAR parents

pro shyInvPlot, doRun=doRun, doParMass=doParMass

  ; config
  res       = 512
  redshifts = [5.0,4.0,3.5,3.0,2.5,2.0,1.5,1.0,0.75,0.5,0.25]
  runs      = ['feedback','tracer']
  massBins  = list( [10.0,10.5], [10.5,11.0], [11.0,11.5], [11.5,12.0] )
  noParMass = 1
  
  if n_elements(doParMass) gt 0 then noParMass = doParMass
  if keyword_set(doRun) then runs = doRun
  
  ; plot config
  lines  = [1,0,2]    ; linestyle
  syms   = [-4,-6]    ; psym
  cInd   = 1          ; color index
  xrange = [-0.1,5.1] ; redshift
  yrange = [0.0,1.0]  ; fraction
  
  ; load
  sP = simParams(res=res,run=runs[0],redshift=0.0)
  if noParMass eq 1 then pmStr = '.noParMass' else pmStr = ''
  saveFilename = sP.derivPath + 'shyPlot/shyInvPlotRedshiftBin' + pmStr + '.sav'
  
  if file_test(saveFilename) then begin
    restore, saveFilename
  endif else begin
    foreach run,runs do begin
      fracs_old    = []
      fracs_full   = []
      fracs_fullNW = []
      
      mb_fracs_old    = {}
      mb_fracs_full   = {}
      mb_fracs_fullNW = {}
      foreach massBin,massBins,k do begin
        mb_fracs_old  = mod_struct(mb_fracs_old, 'mb'+str(k),fltarr(n_elements(redshifts))+!values.f_nan)
        mb_fracs_full = mod_struct(mb_fracs_full,'mb'+str(k),fltarr(n_elements(redshifts))+!values.f_nan)
        mb_fracs_fullNW = mod_struct(mb_fracs_fullNW,'mb'+str(k),fltarr(n_elements(redshifts))+!values.f_nan)
      endforeach
      
      foreach redshift,redshifts,j do begin
        print,run,redshift
        sP = simParams(res=res,run=run,redshift=redshift)
        
        ; only z=0 parents
        iB = shyInvPlotBin(sP=sP)
      
        fracs_old = [fracs_old, iB.frac_to_stars]
        
        ; full snapshot walk (all parttype4)
        iBF = shyFullInvBin(sP=sP,/incWind,/getResult)
        
        if n_elements(iBF) gt 0 then fracs_full = [fracs_full, iBf.frac_to_stars]
        if n_elements(iBF) eq 0 then fracs_full = [fracs_full, !values.f_nan]
        
        ; full snapshot walk (noWind)
        iBFnW = shyFullInvBin(sP=sP,/getResult)
        
        if n_elements(iBFnW) gt 0 then fracs_fullNW = [fracs_fullNW, iBfnW.frac_to_stars]
        if n_elements(iBFnW) eq 0 then fracs_fullNW = [fracs_fullNW, !values.f_nan]
        
        ; get the missing information we need to compute by massbin
        if noParMass eq 0 then begin
        parMass = 1
        tr_ids_acc = shyInvPlotAccretedTracers( sP=sP, parMass=parMass )
        
        ; now also bin by parent mass
        foreach massBin,massBins,k do begin
          ; only z=0
          w = where(iB.parMass ge (massBins[k])[0] and iB.parMass lt (massBins[k])[1],count_all)
          tr_int = intersection( tr_ids_acc[w], iB.tr_ids_stars )
          if count_all gt 0 then mb_fracs_old.(k)[j] = float( n_elements(tr_int) ) / count_all
          
          ; full
          if n_elements(iBF) gt 0 then begin
            w = where(iBF.parMass ge (massBins[k])[0] and iBF.parMass lt (massBins[k])[1],count_1)
            w = where(parMass ge (massBins[k])[0] and parMass lt (massBins[k])[1],count_all)
            if count_all gt 0 then mb_fracs_full.(k)[j] = float(count_1) / count_all
          endif
          
          ; full (noWind)
          if n_elements(iBFnW) gt 0 then begin
            w = where(iBFnW.parMass ge (massBins[k])[0] and iBFnW.parMass lt (massBins[k])[1],count_1)
            w = where(parMass ge (massBins[k])[0] and parMass lt (massBins[k])[1],count_all)
            if count_all gt 0 then mb_fracs_fullNW.(k)[j] = float(count_1) / count_all
          endif
        endforeach
        endif ; noParMass
          
      endforeach
    
      stampStruct = {fracs_old:fracs_old, fracs_full:fracs_full, fracs_fullNW:fracs_fullNW, $
                     mb_fracs_old:mb_fracs_old, mb_fracs_full:mb_fracs_full, $
                     mb_fracs_fullNW:mb_fracs_fullNW}
      runFracs = mod_struct( runFracs, run, stampStruct )
    endforeach
    
    save,runFracs,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
  endelse
  
  ; plot (1) - all halo mass
  start_PS,sP.plotPath + 'shyInvPlotRedshift_'+str(res)+'.eps'
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,$
      title="",xtitle="Redshift",ytitle="Fraction",/xs,/ys
    
    cgPlot,xrange+[0.2,-0.2],[0.5,0.5],line=0,color=cgColor('light gray'),/overplot
      
    legendNames  = []
    legendColors = []
    
    foreach run,runs,i do begin
      sP = simParams(res=res,run=run)
      
      legendNames  = [legendNames,sP.simName]
      legendColors = [legendColors,sP.colors[cInd]]
      
      cgPlot, redshifts, runFracs.(i).fracs_old, $
        line=lines[0], color=sP.colors[cInd], /overplot
      cgPlot, redshifts, runFracs.(i).fracs_full, $
        psym=syms[0], line=lines[1], color=sP.colors[cInd], /overplot
      cgPlot, redshifts, runFracs.(i).fracs_fullNW, $
        psym=syms[1], line=lines[2], color=sP.colors[cInd], /overplot
    endforeach
	
    ; legend
    legend,legendNames,textcolors=legendColors,/top,/left
    legend,['only z=0 parent check','full walk to z=0 (all PT4)','full walk to z=0 (no wind)'],$
      psym=[0,syms[0],syms[1]],linestyle=lines,/bottom,/right ;psym=syms,
      
  end_PS
  
  ; plot (2) - binned into 2x2 by halo mass
  start_PS, sP.plotPath + 'shyInvPlotRedshift_massBin_'+str(res)+'.eps', /big
    pos = plot_pos(col=2, row=2, /gap)
    
    for k=0,3 do begin
      ; plot
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,pos=pos[k],noerase=(k gt 0),$
        title="",xtitle="Redshift",ytitle="Fraction",/xs,/ys
      
      cgPlot,xrange+[0.2,-0.2],[0.5,0.5],line=0,color=cgColor('light gray'),/overplot
    
      foreach run,runs,i do begin
        sP = simParams(res=res,run=run)
      
        cgPlot, redshifts, runFracs.(i).mb_fracs_old.(k), $
          line=lines[0], color=sP.colors[cInd], /overplot
        cgPlot, redshifts, runFracs.(i).mb_fracs_full.(k), $
          psym=syms[0], line=lines[1], color=sP.colors[cInd], /overplot
        cgPlot, redshifts, runFracs.(i).mb_fracs_fullNW.(k), $
          psym=syms[1], line=lines[2], color=sP.colors[cInd], /overplot
      endforeach
    
      ; legends
      legend,legendNames,textcolors=legendColors,/top,/left
      massBinStr = string((massBins[k])[0],format='(f4.1)')+' < M < '+$
                   string((massBins[k])[1],format='(f4.1)')
      legend,massBinStr,/bottom,/right
      
    endfor ;k
  end_PS
  
  ; write out text
  openw,lun,'dumpInv.txt',width=400,/get_lun
  printf,lun,'redshifts: ',redshifts
  foreach run,tag_names(runFracs),i do begin
    printf,lun,run+' fracs_z0:         ',runFracs.(i).fracs_old
    printf,lun,run+' fracs_fullWWind:  ',runFracs.(i).fracs_full
    printf,lun,run+' fracs_fullNoWind: ',runFracs.(i).fracs_fullNW
  endforeach
  
  foreach run,tag_names(runFracs),i do begin
    foreach massBin,tag_names(runFracs.(i).mb_fracs_old),k do begin
      printf,lun,run+' massBin ('+string( (massBins[k])[0], format='(f4.1)' )+' - '+$
                                  string( (massBins[k])[1], format='(f4.1)' )+') fracs_z0:   ',$
            runFracs.(i).mb_fracs_old.(k)
      printf,lun,run+' massBin ('+string( (massBins[k])[0], format='(f4.1)' )+' - '+$
                                  string( (massBins[k])[1], format='(f4.1)' )+') fracs_fullWWind:   ',$
            runFracs.(i).mb_fracs_full.(k)
      printf,lun,run+' massBin ('+string( (massBins[k])[0], format='(f4.1)' )+' - '+$
                                  string( (massBins[k])[1], format='(f4.1)' )+') fracs_fullNoWind:   ',$
            runFracs.(i).mb_fracs_fullNW.(k)
    endforeach
  endforeach
  
  close,lun
  free_lun,lun
  print,'Wrote [dumpInv.txt].'
  
  stop

end

; shyPlot():
;  plot (1): delta time between "halo accretion time" and "star forming time" (Gyr)
;  plot (2): delta time as above, normalized by age of universe at "halo accretion time"
; aMode=1: "halo accretion time" = first 1.0 rvir crossing from accretionTimes()
; aMode=2: "halo accretion time" = group membership
; allTypes=1: all tracers in galaxyCat, otherwise just 'star' (e.g. inside 0.15rvir cut) tracers

pro shyPlot, redshifts=redshifts, allTypes=allTypes, modeIn=modeIn, doRun=doRun

  ; config
  ;redshifts = [3.0,2.0,1.0,0.0] ;c
  ;redshifts = [5.0,4.0,3.5,2.5] ;a
  ;redshifts = [1.5,0.75,0.5,0.25] ;b
  runs      = ['feedback','tracer']
  res       = 512
  aMode     = 2
  ;allTypes  = 0
  
  if keyword_set(modeIn) then aMode = modeIn
  if keyword_set(doRun)  then runs  = doRun
  
  ; plot config
  yrangeRaw  = [0.0,0.2] ;[0.0,0.1]
  yrangeNorm = [0.0,0.3]
  nBins     = 20
  nBins2    = 20
  colors    = ['red','blue'] ; one per run
  
  sP = simParams(res=res,run=runs[0])
  
  start_PS,sP.plotPath+'shyPlot_aMode'+str(aMode)+'_allTypes'+str(allTypes)+'_'+str(res)+'b.eps', xs=3*3.5, ys=4*3.5
  
    pos = plot_pos(rows=4,cols=2,/gap)
    
    foreach redshift,redshifts,i do begin
    
    xrange = [-0.1, redshiftToAgeFlat(redshift)*1.0]
    xrange2 = [-0.1, 2.0]
    
    ; plot (1)
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeNorm,title="",$
      xtitle=textoidl("( t_{SF} - t_{acc} ) / t_{H}(t_{acc})"),$
      ytitle="Fraction",/xs,/ys,pos=pos[i*2],/noerase
	
	; plot histogram for each run/redshift combination
      foreach run,runs,j do begin
        sP = simParams(res=res,run=run,redshift=redshift)
        bV = shyPlotBin(sP=sP, aMode=aMode, allTypes=allTypes)
        
        val = (bV.age_cur - bV.age_at) / bV.age_at
        
        binSize = (xrange[1] - xrange[0]) / float(nBins)
	  h = histogram(val, bin=binSize, loc=loc, min=xrange[0], max=xrange[1])
	  h = float(h)/total(h)
	  
	  cgPlot,loc+binSize*0.5,h,color=cgColor(colors[j]),/overplot
	endforeach
	
	; legend
	legend,['z = '+string(redshift,format='(f3.1)'),runs],textcolors=['black',colors],/top,/right
	
    ; plot (2)
    cgPlot,[0],[0],/nodata,xrange=xrange2,yrange=yrangeRaw,title="",$
      xtitle=textoidl("( t_{SF} - t_{acc} ) / t_{H}(t_{acc})"),$
      ytitle="Fraction",/xs,/ys,pos=pos[i*2+1],/noerase
	
	; plot histogram for each run/redshift combination
      foreach run,runs,j do begin
        sP = simParams(res=res,run=run,redshift=redshift)
        bV = shyPlotBin(sP=sP, aMode=aMode, allTypes=allTypes)
        
        val = (bV.age_cur - bV.age_at) / bV.age_at
        
        binSize = (xrange2[1] - xrange2[0]) / float(nBins2)
	  h = histogram(val, bin=binSize, loc=loc, min=xrange2[0], max=xrange2[1])
	  h = float(h)/total(h)
	  
	  cgPlot,loc+binSize*0.5,h,color=cgColor(colors[j]),/overplot
	endforeach
	
	; legend
	legend,['z = '+string(redshift,format='(f3.1)'),runs],textcolors=['black',colors],/top,/right
	
    endforeach ;redshifts
      
  end_PS
  
end

; shyPlotRedshift(): plot mean vs redshift

pro shyPlotRedshift, doRun=doRun

  ; config
  runs      = ['feedback','tracer']
  res       = 512
  aModes    = [1,2]
  allTypes  = [0,1]
  redshifts = [5.0,4.0,3.5,3.0,2.5,2.0,1.5,1.0,0.75,0.5,0.25,0.0]
  massBins  = list( [10.0,10.5], [10.5,11.0], [11.0,11.5], [11.5,12.0] )
  
  if keyword_set(doRun) then runs = doRun
  
  ; plot config
  xrange    = [-0.1,5.1] ; redshift
  yrange    = [0.0,2.5] ; delta_t/t_H
  cInd      = 1
  lines     = [0,1,2,3,4]
  syms      = [-4,-5,-6,-7,-3]
  skipInv   = 0 ; 0=include inverse measurements, 1=do not plot
  
  ; set saveFilename and check for existence
  sP = simParams(res=res,run=runs[0],redshift=0.0)
  saveFilename = sP.derivPath + 'shyPlot/shyPlotRedshiftBin.sav'
  
  if file_test(saveFilename) then begin
    restore, saveFilename
  endif else begin
  
    ; load
    foreach run,runs,j do begin
      print,run
      runSet  = {}
      
      foreach aMode,aModes do begin
      foreach allType,allTypes do begin
      
        val_mean   = fltarr(n_elements(redshifts)) + !values.f_nan
        val_median = val_mean

        mb_mean   = {}
        mb_median = {}
        foreach massBin,massBins,k do begin
          mb_mean   = mod_struct( mb_mean,   "mb"+str(k), val_mean )
          mb_median = mod_struct( mb_median, "mb"+str(k), val_median )
        endforeach
        
        foreach redshift,redshifts,i do begin
          sP = simParams(res=res,run=run,redshift=redshift)

          ; if doing single run, we are trying to fill in points, so make sure save exists, if not, make  it
          if keyword_set(doRun) then begin
            aa = shyPlotBin(sP=sP, aMode=aMode, allTypes=allType)
            continue            
          endif
     
          ; NORMAL (backwards in time)
          bV = shyPlotBin(sP=sP, aMode=aMode, allTypes=allType, /blindSearch)
      
          if tag_exist(bV,'flag') then continue
          print,run,redshift,aMode,allType
          
          ; process data, add to linear array
          val = (bV.age_cur - bV.age_at) / bV.age_at
          
          val_mean[i]   = mean(val)
          val_median[i] = median(val)
          
          ; also do mean/medians by mass bin
          if tag_exist(bV,"parMass") then begin
            if n_elements(bV.parMass) ne n_elements(bV.age_at) then message,'Error'
          
            foreach massBin,massBins,k do begin
              w = where(bV.parMass ge (massBins[k])[0] and bV.parMass lt (massBins[k])[1],count)
              if count gt 0 then mb_mean.(k)[i]   = mean(val[w])
              if count gt 0 then mb_median.(k)[i] = median(val[w])
            endforeach
          endif
          
        endforeach ;redshifts
        
        ; add this (aMode,allType) combination to runSet
        sName  = 'aM'+str(aMode)+'_aT'+str(allType)
        runSet = mod_struct( runSet, sName, { aMode:aMode, allType:allType, $
                                              val_mean:val_mean, val_median:val_median,$
                                              mb_mean:mb_mean, mb_median:mb_median } )
        
      endforeach ;allTypes
      endforeach ;allModes
      
      ; INVERSE (forwards in time - no aMode,allType permutations)
      val_mean   = fltarr(n_elements(redshifts)) + !values.f_nan
      val_median = val_mean

      mb_mean   = {}
      mb_median = {}
      foreach massBin,massBins,k do begin
        mb_mean   = mod_struct( mb_mean,   "mb"+str(k), val_mean )
        mb_median = mod_struct( mb_median, "mb"+str(k), val_median )
      endforeach
        
      foreach redshift,redshifts,i do begin
        sP = simParams(res=res,run=run,redshift=redshift)
        
        bVF = shyFullInvBin(sP=sP, /getResult)
        
        if n_elements(bVF) eq 0 then continue
        
        ; process data, add to linear array
        val = (bVF.age_sf - bVF.age_at) / bVF.age_at
          
        val_mean[i]   = mean(val)
        val_median[i] = median(val)
         
        ; also do mean/medians by mass bin
        if tag_exist(bVF,"parMass") then begin          
          foreach massBin,massBins,k do begin
            w = where(bVF.parMass ge (massBins[k])[0] and bVF.parMass lt (massBins[k])[1],count)
            if count gt 0 then mb_mean.(k)[i]   = mean(val[w])
            if count gt 0 then mb_median.(k)[i] = median(val[w])
          endforeach
        endif
      endforeach
      
      runSet = mod_struct( runSet, 'inv', {val_mean:val_mean, val_median:val_median,$
                                           mb_mean:mb_mean, mb_median:mb_median} )

      ; add runSet to keeper struct
      r = mod_struct( r, run, runSet )    

    endforeach ;runs
  
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
  
  endelse ; file_test(saveFilename)
  
  ; plot (1) - all mass bins, mean
  start_PS,sP.plotPath + 'shyPlotRedshift_'+str(res)+'_mean.eps'
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,title="(mean)",$
      xtitle="Redshift",ytitle=textoidl("( t_{SF} - t_{acc} ) / t_{H}(t_{acc})"),/xs,/ys
    
    cgPlot,xrange+[0.2,-0.2],[0.8,0.8],line=0,color=cgColor('light gray'),/overplot
      
    ; plot statistic for each save file that exists
    legendStrs   = []
    legendLines  = []
    legendSym    = []
    
    legendNames  = []
    legendColors = []
    
    foreach run,runs,i do begin
      sP = simParams(res=res,run=run)
      
      lineInd = 0
      legendNames  = [legendNames,sP.simName]
      legendColors = [legendColors,sP.colors[cInd]]
      
      for j=0,n_tags(r.(i))-1-skipInv do begin
        ; mean
        cgPlot,redshifts,r.(i).(j).val_mean,psym=syms[lineInd],line=lines[lineInd],color=sP.colors[cInd],/overplot
        
        ; add to legend
        if i eq 0 then begin
          if j ne n_tags(r.(i))-1 then $
            legendStrs   = [legendStrs,'aMode='+str(r.(i).(j).aMode)+' allType='+str(r.(i).(j).allType)+'']
          if j eq n_tags(r.(i))-1 then $
            legendStrs   = [legendStrs,'(inverse)']
          legendLines  = [legendLines,lines[lineInd+0]]
          legendSym    = [legendSym,syms[lineInd]]
        endif
        
        lineInd += 1

      endfor
    endforeach
	
    ; legend
    legend,legendStrs,linestyle=legendLines,psym=legendSym,/top,/right
    legend,legendNames,textcolors=legendColors,/bottom,/left
      
  end_PS
  
  ; plot (2) - all mass bins, median
  start_PS,sP.plotPath + 'shyPlotRedshift_'+str(res)+'_median.eps'
  
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,title="(median)",$
      xtitle="Redshift",ytitle=textoidl("( t_{SF} - t_{acc} ) / t_{H}(t_{acc})"),/xs,/ys
    
    cgPlot,xrange+[0.2,-0.2],[0.8,0.8],line=0,color=cgColor('light gray'),/overplot
      
    ; plot statistic for each save file that exists
    foreach run,runs,i do begin
      sP = simParams(res=res,run=run)
      
      lineInd = 0
      
      for j=0,n_tags(r.(i))-1-skipInv do begin
        ; median
        cgPlot,redshifts,r.(i).(j).val_median,psym=syms[lineInd],line=lines[lineInd],color=sP.colors[cInd],/overplot
        
        lineInd += 1

      endfor
    endforeach
	
    ; legend
    legend,legendStrs,linestyle=legendLines,psym=legendSym,/top,/right
    legend,legendNames,textcolors=legendColors,/bottom,/left
      
  end_PS
  
  ; plot (3) - 2x2 for four massbins, mean
  start_PS, sP.plotPath + 'shyPlotRedshift_massBin_'+str(res)+'_mean.eps', /extrabig
    pos = plot_pos(col=2, row=2, /gap)
    
    for k=0,3 do begin
      ; plot
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,pos=pos[k],noerase=(k gt 0),$
        title="",xtitle="Redshift",ytitle=textoidl("( t_{SF} - t_{acc} ) / t_{H}(t_{acc})"),/xs,/ys
      
      cgPlot,xrange+[0.2,-0.2],[0.8,0.8],line=0,color=cgColor('light gray'),/overplot
    
      foreach run,runs,i do begin
        sP = simParams(res=res,run=run)
      
        lineInd = 0
      
        for j=0,n_tags(r.(i))-1-skipInv do begin
          ; mean
          cgPlot,redshifts,r.(i).(j).mb_mean.(k),psym=syms[lineInd],line=lines[lineInd],color=sP.colors[cInd],/overplot
          lineInd += 1
        endfor
      endforeach
    
      ; legends
      if k eq 1 then $
        legend,legendNames,textcolors=legendColors,/top,/right
      massBinStr = string((massBins[k])[0],format='(f4.1)')+' < M < '+$
                   string((massBins[k])[1],format='(f4.1)')
      legend,massBinStr,textcolors=[cgColor('orange')],/bottom,/left
      if k eq 0 then $
        legend,legendStrs,linestyle=legendLines,psym=legendSym,/top,/right
      
    endfor ;k
  end_PS
  
  ; write out text
  openw,lun,'dump.txt',width=400,/get_lun
  printf,lun,'redshifts: ',redshifts
  foreach run,tag_names(r),i do begin
    for j=0,n_tags(r.(i))-1-skipInv do begin
      method = ( tag_names(r.(i)) )[j]
      printf,lun,run+' '+method+' mean:   ',r.(i).(j).val_mean
      printf,lun,run+' '+method+' median: ',r.(i).(j).val_median
    endfor
  endforeach
  
  foreach run,tag_names(r),i do begin
    for j=0,n_tags(r.(i))-1-skipInv do begin
      method = ( tag_names(r.(i)) )[j]
      foreach massBin,tag_names(r.(i).(j).mb_mean),k do begin
        printf,lun,run+' '+method+' massBin ('+string( (massBins[k])[0], format='(f4.1)' )+' - '+$
                                          string( (massBins[k])[1], format='(f4.1)' )+') mean:   ',$
              r.(i).(j).mb_mean.(k)
        printf,lun,run+' '+method+' massBin ('+string( (massBins[k])[0], format='(f4.1)' )+' - '+$
                                          string( (massBins[k])[1], format='(f4.1)' )+') median: ',$
              r.(i).(j).mb_median.(k)
      endforeach
    endfor
  endforeach
  
  close,lun
  free_lun,lun
  print,'Wrote [dump.txt].'
  
  stop
  
end

; shyPlot2D(): plot 2d histogram

pro shyPlot2D

  ; config
  runs      = ['feedback','tracer']
  res       = 512
  redshifts = [5.0,4.0,3.5,3.0,2.5,2.0,1.5,1.0,0.75,0.5,0.25]
  massBins  = list( [10.0,10.5], [10.5,11.0], [11.0,11.5], [11.5,12.0] )
  atBinSize = 0.3 ; Gyr
  atBinMM   = [0.0,14.0]
  
  ; plot config
  xrange    = [0.0,14.0] ; age of universe
  yrange    = [0.0,14.0] ; age of universe
  yrangeF   = [1e-3,0.4] ; fraction
  cInd      = 1
  sK        = 7
  
  ; set saveFilename and check for existence
  nBins   = floor(( atBinMM[1] - atBinMM[0] ) / atBinSize) + 1
  binCens = linspace( atBinMM[0], atBinMM[1], nBins )
  
  sP = simParams(res=res,run=runs[0],redshift=0.0)
  saveFilename = sP.derivPath + 'shyPlot/shyPlotRedshiftBin2D.nB'+str(nBins)+'.sav'
  
  if file_test(saveFilename) then begin
    restore, saveFilename
  endif else begin
  
    ; load
    foreach run,runs,j do begin
      ; INVERSE (forwards in time)
      val_mean = fltarr(n_elements(redshifts),nBins) + !values.f_nan
      val_tot  = lonarr(n_elements(redshifts)) + !values.f_nan
      val_acc  = fltarr(n_elements(redshifts)) + !values.f_nan

      mb_mean = {}
      foreach massBin,massBins,k do begin
        mb_mean = mod_struct( mb_mean,   "mb"+str(k), val_mean )
      endforeach
        
      foreach redshift,redshifts,i do begin
        print,run,redshift
        sP = simParams(res=res,run=run,redshift=redshift)
        
        bVF = shyFullInvBin(sP=sP, /getResult)
        
        if n_elements(bVF) eq 0 then continue
        
        ; process data, add to linear array
        h = histogram( bVF.age_sf, binSize=atBinSize, min=atBinMM[0], max=atBinMM[1], loc=loc )

        val_mean[i,*] = h
        val_tot[i]    = n_elements(bVF.age_sf)
        val_acc[i]    = bVF.age_at
         
        ; also do mean/medians by mass bin
        if tag_exist(bVF,"parMass") then begin          
          foreach massBin,massBins,k do begin
            w = where(bVF.parMass ge (massBins[k])[0] and bVF.parMass lt (massBins[k])[1],count)
            
            if count eq 0 then continue
            
            h = histogram( bVF.age_sf[w], binSize=atBinSize, min=atBinMM[0], max=atBinMM[1], loc=loc )
            mb_mean.(k)[i,*] = h
          endforeach
        endif
      endforeach
      
      ; add runSet to keeper struct
      r = mod_struct(r,run,{val_mean:val_mean,mb_mean:mb_mean,val_tot:val_tot,val_acc:val_acc,loc:loc})    

    endforeach ;runs
  
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
  
  endelse ; file_test(saveFilename)
  
  ; plot init
  set_plot,'ps'
  zColors = reverse( sampleColorTable('blue-red2', n_elements(redshifts), bounds=[0.1,0.9]) )
  
  ; plot (1) - all mass bins, 1d profiles from each redshift
  start_PS,sP.plotPath + 'shyPlot2DRedshift_'+str(res)+'_nB'+str(nBins)+'_profiles.eps', ys=10.0, xs=7.5
    pos = plot_pos(row=2,col=1,/gap)
    
    ; RUN 1
    i = 0
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeF,title="",/ylog,yminor=0,pos=pos[i],$
      xtitle=textoidl("t_{age} [Gyr] (SF)"),ytitle=textoidl("Fraction"),/xs,/ys
    
    cgPlot,xrange+[0.2,-0.2],[0.5,0.5],line=0,color=cgColor('light gray'),/overplot
      
    sP = simParams(res=res,run=runs[i])
      
    for j=0,n_elements(redshifts)-1 do begin
      ; mean
      yy = r.(i).val_mean[j,*] / float(r.(i).val_tot[j])
        
      cgPlot,binCens,smooth(yy,sK),line=0,color=zColors[j],/overplot
    endfor
	
    ; legend
    legend,[sP.simName],textcolors=[sP.colors[cInd]],/top,/right
    legend,textoidl("z_{acc} = ")+string(redshifts[0:-1:2],format='(f4.2)'),textcolors=zColors[0:-1:2],/top,/left
    
    ; RUN 2
    i = 1
    cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrangeF,title="",/ylog,yminor=0,pos=pos[i],$
      xtitle=textoidl("t_{age} [Gyr] (SF)"),ytitle=textoidl("Fraction"),/xs,/ys,/noerase
    
    cgPlot,xrange+[0.2,-0.2],[0.5,0.5],line=0,color=cgColor('light gray'),/overplot
    
    sP = simParams(res=res,run=runs[i])
      
    for j=0,n_elements(redshifts)-1 do begin
      ; mean
      yy = r.(i).val_mean[j,*] / float(r.(i).val_tot[j])
        
      cgPlot,binCens,smooth(yy,sK),line=0,color=zColors[j],/overplot
    endfor
    
    ; legend
    legend,[sP.simName],textcolors=[sP.colors[cInd]],/top,/right
      
  end_PS
  
  ; plot (2) - 2d histogram
  start_PS,sP.plotPath + 'shyPlot2DRedshift_'+str(res)+'_nB'+str(nBins)+'.eps', xs=7.5, ys=10.0
    pos = plot_pos(row=2,col=1,/gap)
    
    for i=0,1 do begin
      sP = simParams(res=res,run=runs[i])
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,title="",pos=pos[i],noerase=(i gt 0),$
        xtitle=textoidl("t_{age} [Gyr] (ACC)"),$
        ytitle=textoidl("t_{age} [Gyr] (SF)"),/xs,/ys
      
      cgPlot,xrange,yrange,line=2,color=cgColor('light gray'),/overplot
      
      ; form 2d histogram
      f2d = { h2       : smooth(r.(i).val_mean,[1,3]) ,$
              nXBins   : n_elements(r.(i).val_acc) ,$
              binSizeX : 0.3                       ,$
              binCenX  : r.(i).val_acc             ,$
              nYBins   : nBins                     ,$
              binSizeY : atBinSize*1.0             ,$
              binCenY  : binCens                    }
    
      ; plot 2d histogram
      hsp = [0.04,0.04]
      nc = 230
      oplot2DHistoSq, f2d, hsp=hsp, nc=nc, xrange=xrange, yrange=yrange, ctName='brewerC-purplebluegreen', /colNorm
      legend,[sP.simName],textcolors=[sP.colors[cInd]],/bottom,/right
      
    endfor ;i
  end_PS
  
  ;NN = 1000000000
  ;;NNN = 100
  ;zz = ptrarr(NNN)
  ;for i=0,NNN-1 do begin
  ;  zz[i] = ptr_new( lindgen(NN) )
  ;  reportMemory, msg=str(i)
  ;endfor
  
  ; write out text
  openw,lun,'dump2D.txt',width=400,/get_lun
  printf,lun,'redshifts: ',redshifts
  printf,lun,''
  printf,lun,'t_SF_binCens: ',binCens
  printf,lun,''
  foreach run,tag_names(r),i do begin
    for j=0,n_elements(redshifts)-1 do begin
      printf,lun,run+' t_ACC='+str(r.(i).val_acc[j])+':   '
      printf,lun,r.(i).val_mean[j,*]
      printf,lun,''
    endfor
  endforeach
  close,lun
  free_lun,lun
  print,'Wrote [dump.txt].'
  
  
  stop
  
end

; shyMassRatio():

pro shyMassRatio

  ; config
  sP = simParams(res=256,run='feedback',redshift=2.0)


  ; load at target redshift
  galcat = galaxyCat(sP=sP)
  
  ; ----------------------------------------
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
  ;tr_inds_orig1 = par_ind[tr_parids_ind]
  tr_inds_orig1 = tr_parids_ind
  tr_ids_orig1  = galcat.trMC_ids[ tr_inds_orig1 ]
    
  ; ----------------------------------------
  print,'AllTypes = 0'
  type = ( galcat.type[ replicate_var(galcat.trMC_cc) ] )
  
  tr_inds_orig0 = where(type eq galcat.types.stars)
  tr_ids_orig0 = galcat.trMC_ids[ tr_inds_orig0 ]

  ; -----------------------------------------
  help,tr_ids_orig0
  help,tr_ids_orig1
  calcMatch,tr_ids_orig0,tr_ids_orig1,ind1,ind2,count=countMatch
  print,'Matched between 0 and 1: ',countMatch
  
  ; get all tracer children of all stars
  tr_children_stars = cosmoTracerChildren(sP=sP,/getIDs,starIDs=star_ids)
  
  calcMatch,tr_children_stars,tr_ids_orig0,ind1,ind2,count=countMatch
  print,'Matched between all and 0: ',countMatch
  
  calcMatch,tr_children_stars,tr_ids_orig1,ind1,ind2,count=countMatch
  print,'Matched between all and 1: ',countMatch
  stop
  
end
