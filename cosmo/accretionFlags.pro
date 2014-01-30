; accretionFlags.pro
; feedback project - flagging method for defining accreted particles/tracers within a timeWindow
; dnelson dec.2013

; accretionFlags(): testing new procedure for measuring accretion rates based on group membership at a preivous snapshot

function accretionFlags, sP=sP, targetSnap=targetSnap

  forward_function cosmoTracerChildren, cosmoTracerVelParents
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  if ~keyword_set(sP) or n_elements(targetSnap) eq 0 then message,'Error'

  ; set saveFilename and check for existence
  saveFilename = sP.derivPath + 'accFlags.'+sP.saveTag+'.'+sP.savPrefix+str(sP.res)+'.'+$
                 str(sP.snap)+'-'+str(targetSnap)+'.sav'
  
  if file_test(saveFilename) then begin
    restore, saveFilename
    return, r
  endif
  
  ; first, walk back through the merger tree and find primary subhalos with good parent histories
  origSnap = sP.snap
  galcat   = galaxyCat(sP=sP)
  mt       = mergerTreeSubset(sP=sP)
  
  ; ALL CASES
  ; ---------
  hInds = lindgen(n_elements(mt.hMinSnap))
    
  ; walk the mergerTree, follow the Parent's back to the targetSnap, such that the output
  ; hInds then contains the mapping of the original indices into their parents at targetSnap
  for m=origSnap,targetSnap+1,-1 do begin
    sP.snap = m
           
    ; load mergerTree and move to Parent
    Parent = mergerTree(sP=sP)
        
    ; change to parent IDs for the next snapshot
    w = where(hInds ne -1,count)
    if count eq 0 then message,'error'
        
    hInds[w] = Parent[hInds[w]]
        
    frac = float(count)*100/n_elements(mt.hMinSnap)
    print,'['+str(m)+'] frac remaining = '+string(frac,format='(f4.1)')+'% ('+$
      str(count)+' of '+str(n_elements(mt.hMinSnap))+')'
          
    ; sanity check no IDs are -1 (we should only be processing the mergerTreeSubset)
    w = where(hInds lt 0,count)        
    if count gt 0 then $
      if min( mt.hMinSnap[w] ) lt sP.snap then message,'Error: Bad parent.'
      
  endfor ; snapshot
    
  ; load/create the galaxy catalog at this previous snapshot
  sP.snap  = targetSnap
  galcatTarget = galaxyCat(sP=sP)
  mtTarget     = mergerTreeSubset(sP=sP)
  
  ; make inverse hInds mtS mapping (e.g. hIndsInv[targetSnapParInd] = origSnapParInd
  hIndsInv = lonarr(galcatTarget.nGroups) - 1
  w = where(hInds ne -1,count)
  if count eq 0 then message,'Error'
  hIndsInv[hInds[w]] = w
  
  ; NO TRACERS CASE - track the particles themselves back in time (SPH)
  ; ---------------
  if sP.trMCPerCell eq 0 then begin
  
    print,'Calculating new accretion flags using ( SPH ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'

    nTotSearch = n_elements(galcat.ids)
    nTotSearchTarg = n_elements(galcatTarget.ids)
    
    ; map the parent index list to a new parentIndex array at targetSnap
    origParIDs = hInds[ mt.gcIndOrig ]    
    
    ; make association of galaxycat IDs between the two galaxy catalogs
    match,galcat.ids,galcatTarget.ids,gcOrig_inds,gcTarget_inds,count=countMatch
    sort_inds = sort(gcOrig_inds)
    gcOrig_inds = gcOrig_inds[ sort_inds ]
    gcTarget_inds = gcTarget_inds[ sort_inds ]
    sort_inds = !NULL
  
    frac = float(countMatch)*100 / nTotSearch
    print,' Matched fraction: '+string(frac,format='(f4.1)')+'% ('+$
      str(countMatch)+' of '+str(nTotSearch)+')'
      
    ; mapping of members from galcat.ids to galcatTarget.ids
    gcMapInds = lonarr(nTotSearch) - 1
    gcMapInds[gcOrig_inds] = gcTarget_inds
    
    ; mapping of members from galcatTarget.ids to galcat.ids
    gcMapIndsInv = lonarr(nTotSearchTarg) - 1
    gcMapIndsInv[gcTarget_inds] = gcOrig_inds
    
    ; replicate inverse hInds mapping for each targetSnap member, and make mask for -matched- members
    origParIDsTarg = hIndsInv[ mtTarget.gcIndOrig ]
    
    targMatchedMask = bytarr(nTotSearchTarg)
    targMatchedMask[gcTarget_inds] = 1B
    
    ; transform parentIndices at origSnap, to be ordered as one origSnapParent for each targetSnap member
    gcParIndRep = galcatRepParentIDs(galcat=galcat)
    
    galcatIndOrig = lonarr(nTotSearchTarg) - 2 ; default -2, see below
    
    w = where(gcMapIndsInv ne -1,count)
    if count eq 0 then message,'Error'
    galcatIndOrig[w] = gcParIndRep[ gcMapIndsInv[w] ]
    
    gcParIndRep = !NULL
    gcOrig_inds = !NULL
    gcTarget_inds = !NULL
  
    ; transform parentIndices at targetSnap, to be ordered as one targetSnapParent for each origSnap member
    gcTargParIndRep = galCatRepParentIDs(galcat=galcatTarget)
    
    galcatTargetIndOrig = lonarr(nTotSearch) - 2 ; default -2, to prevent match to -1 below
    
    w = where(gcMapInds ne -1,count)
    if count eq 0 then message,'Error'
    galcatTargetIndOrig[w] = gcTargParIndRep[ gcMapInds[w] ]
    
    gcTargParIndRep = !NULL
    
    ; transform the targetSnap member types to be ordered as one type for each 
    ; origSnap member, with type -2 for unmatched
    gc_targ_type = intarr(nTotSearch) - 2
    
    w = where(gcMapInds ne -1,count)
    if count eq 0 then message,'Error'
    gc_targ_type[w] = galcatTarget.type[ gcMapInds[w] ]
  
    ; start: store the main arrays as a structure so we can write them directly for all tracers
    r = { accFlag           : intarr(nTotSearch)-2        ,$
          accFlagTarg       : intarr(nTotSearchTarg)-2    ,$
          origParIDsTarg    : origParIDsTarg              ,$
          targMatchedMask   : targMatchedMask            }
          
    targMatchedMask = !NULL
    origParIDsTarg = !NULL
          
    gcMaskDebug = intarr(nTotSearch) ; DEBUG
    gcMaskDebugTarg = intarr(nTotSearchTarg) ; DEBUG
    
    ; which of these tracers are also in galcatTarget?
    wMap = where(gcMapInds ge 0,countMap,comp=wNoMap,ncomp=countNoMap)
    
    if countMap eq 0 or countNoMap eq 0 then message,'Error: Unlikely, and below will fail.'
    ; (1a). those which did not map:
      ; flag (1) any gal/stars members from origSnap but not in the progenitor at targetSnap
      gc_types_nomap = galcat.type[ wNoMap ]
        
      w = where(gc_types_nomap eq galcat.types.gal or $
                gc_types_nomap eq galcat.types.stars, count_galAcc_noMap)
                  
      if count_galAcc_noMap eq 0 then message,'Error'
      r.accFlag[wNoMap[w]] = 1
      gcMaskDebug[wNoMap[w]] += 1
      
      ; flag (2) any gmem tracers from origSnap but not in the progenitor at targetSnap
      ; as net halo inflow (tracers in inter at origSnap do not contribute to anything)
      w = where(gc_types_nomap eq galcat.types.gmem, count_gmemAcc_noMap)
        
      if count_gmemAcc_noMap eq 0 then message,'Error'
      r.accFlag[wNoMap[w]] = 2
      gcMaskDebug[wNoMap[w]] += 1
    
    ; (1b). those which did map:
    
      ; (2a). flag (-1) any tracers of halos with failed mergerTree connectivity, do not include these halos
      w = where(origParIDs[wMap] eq -1,count,comp=wTracked)
      frac = float(count)*100 / nTotSearch
      print,' Failed connectivity: '+string(frac,format='(f4.1)')+'% ('+$
      str(count)+' of '+str(nTotSearch)+')'
    
      r.accFlag[wMap[w]] = -1
      gcMaskDebug[wMap[w]] += 1
      
      ; (2b):
      wTracked = wMap[wTracked] ; index galcat.trMC_ids directly
      gc_types_orig  = galcat.type[ wTracked ]
      gc_types_targ  = gc_targ_type[ wTracked ]
      gc_parind_orig = origParIDs[ wTracked ]
      gc_parind_targ = galcatTargetIndOrig[ wTracked ]
      
        ; (3a1): flag (1) any gal/stars members from origSnap in the progenitor 
        ; but in the inter/gmem type, as net galaxy inflow
        w = where((gc_types_orig eq galcat.types.gal   or $
                   gc_types_orig eq galcat.types.stars ) and $
                  (gc_types_targ eq galcat.types.inter or $
                   gc_types_targ eq galcat.types.gmem)   and $
                  gc_parind_orig eq gc_parind_targ, count_galAcc_sameHalo)
                  
        if count_galAcc_sameHalo eq 0 then message,'Error'
        r.accFlag[wTracked[w]] = 1
        gcMaskDebug[wTracked[w]] += 1
      
        ; (3b): in different halo, flag (1) or (2) based on gal/halo type at origSnap
        w = where((gc_types_orig eq galcat.types.gal   or $
                   gc_types_orig eq galcat.types.stars ) and $
                  gc_parind_orig ne gc_parind_targ, count_galAcc_diffHalo)
                  
        if count_galAcc_diffHalo eq 0 then message,'Error'
        r.accFlag[wTracked[w]] = 1
        gcMaskDebug[wTracked[w]] += 1
        
        w = where(gc_types_orig eq galcat.types.gmem and $
                  gc_parind_orig ne gc_parind_targ, count_haloAcc_diffHalo)
                  
        if count_haloAcc_diffHalo eq 0 then message,'Error'
        r.accFlag[wTracked[w]] = 2
        gcMaskDebug[wTracked[w]] += 1
      
    ; fill in the cases where we have no contribution as a verification
        ; (1a) inter tracers not in galcatTarget do not contribute to anything
        w = where(gc_types_nomap eq galcat.types.inter, count)
        
        if count eq 0 then message,'Error'
        r.accFlag[wNoMap[w]] = 0
        gcMaskDebug[wNoMap[w]] += 1
      
        ; (3a2+3a4): in same halo, was inter/gmem in origSnap and in targetSnap
        w = where((gc_types_orig eq galcat.types.inter   or $
                   gc_types_orig eq galcat.types.gmem ) and $
                  gc_parind_orig eq gc_parind_targ, count)
                  
        if count eq 0 then message,'Error'
        r.accFlag[wTracked[w]] = 0
        gcMaskDebug[wTracked[w]] += 1
        
        ; (3a3): in same halo, was gal/stars/bh in origSnap and in targetSnap
        w = where((gc_types_orig eq galcat.types.gal   or $
                   gc_types_orig eq galcat.types.stars ) and $
                  (gc_types_targ eq galcat.types.gal or $
                   gc_types_targ eq galcat.types.stars ) and $
                  gc_parind_orig eq gc_parind_targ, count)
                  
        if count eq 0 then message,'Error'
        r.accFlag[wTracked[w]] = 0
        gcMaskDebug[wTracked[w]] += 1
        
        ; (3b): in different halo, inter in origSnap
        w = where(gc_types_orig eq galcat.types.inter and $
                  gc_parind_orig ne gc_parind_targ, count)
                  
        if count eq 0 then message,'Error'
        r.accFlag[wTracked[w]] = 0
        gcMaskDebug[wTracked[w]] += 1
        
    ; consider the material in the (tracked) progenitor of each halo at targetSnap
    ; which does not appear in the origSnap galaxyCat
      ; (4a): galaxy outflow
      w = where(r.targMatchedMask eq 0B and $
                (galcatTarget.type eq galcat.types.gal or $
                 galcatTarget.type eq galcat.types.stars),count)
      if count eq 0 then message,'Error'

      r.accFlagTarg[w] = 1
      gcMaskDebugTarg[w] += 1
      
      ; (4b): halo outflow
      w = where(r.targMatchedMask eq 0B and $
                galcatTarget.type eq galcat.types.gmem,count)
      if count eq 0 then message,'Error'
      
      r.accFlagTarg[w] = 2
      gcMaskDebugTarg[w] += 1
      
      ; (4c): inter, not unused, just flag
      w = where(r.targMatchedMask eq 0B and $
                galcatTarget.type eq galcat.types.inter,count)
      if count eq 0 then message,'Error'
      
      r.accFlagTarg[w] = 0
      gcMaskDebugTarg[w] += 1
      
    ; consider the material in the (tracked) progenitor of each halo at targetSnap
    ; which is in the origSnap galaxyCat, but in a different descendent
      ; (4b): galaxy outflow
      w = where(r.targMatchedMask eq 1B and $
                (galcatTarget.type eq galcat.types.gal or $
                 galcatTarget.type eq galcat.types.stars) and $
                 galcatIndOrig ne r.origParIDsTarg,count)
      if count eq 0 then message,'Error'

      r.accFlagTarg[w] = 1
      gcMaskDebugTarg[w] += 1
      
      ; (5b): halo outflow
      w = where(r.targMatchedMask eq 1B and $
                galcatTarget.type eq galcat.types.gmem and $
                galcatIndOrig ne r.origParIDsTarg,count)
      if count eq 0 then message,'Error'

      r.accFlagTarg[w] = 2
      gcMaskDebugTarg[w] += 1
      
      ; (5c): 
      w = where(r.targMatchedMask eq 1B and $
                galcatTarget.type eq galcat.types.inter and $
                galcatIndOrig ne r.origParIDsTarg,count)
      if count eq 0 then message,'Error'

      r.accFlagTarg[w] = 0
      gcMaskDebugTarg[w] += 1
      
    ; matched material in targetSnap in the descendent (just mark, does not contribute)
      ; galaxy/halo/inter
      w = where(r.targMatchedMask eq 1B and $
                galcatIndOrig eq r.origParIDsTarg,count)
      if count eq 0 then message,'Error'
      
      r.accFlagTarg[w] = 0
      gcMaskDebugTarg[w] += 1
        
    ; ambiguities:
    ; 1. what about material that is in the galaxy at targetSnap but not at origSnap?
    ;    this is somehow net outflow, do we subtract it or not? (NO)
    ; 2. what about material in gmem at origSnap and in gal/stars/bh/inter at targetSnap?
    ;    this is somehow net inflow into the halo, from within? (do not include as halo acc)
    
    w = where(gcMaskDebug ne 1,count)
    if count gt 0 then message,'Error: Some members not considered.'
     
    w = where(gcMaskDebugTarg ne 1,count)
    if count gt 0 then message,'Error: Some targ members not considered.'
     
    w = where(r.accFlag lt -1 or r.accFlag gt 2,count)
    if count gt 0 then message,'Error: Strange accFlag values.'
    
    w = where(r.accFlagTarg lt -1 or r.accFlag gt 2,count)
    if count gt 0 then message,'Error: Strange accFlagTarg values.'
        
    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
  endif
  
  ; MONTE CARLO TRACERS CASE - for all galaxy catalog members, track back all child tracers
  ; ------------------------
  if sP.trMCPerCell gt 0 then begin
                  
    print,'Calculating new accretion flags using ( TracerMC ) res = '+str(sP.res)+$
      ' in range ['+str(mt.minSnap)+'-'+str(mt.maxSnap)+'].'
      
    nTotSearch = n_elements(galcat.trMC_ids)
    nTotSearchTarg = n_elements(galcatTarget.trMC_ids)
    
    ; map the replicated parent index list to a new parentIndex array at targetSnap
    origParIDs = hInds[ mt.gcIndOrigTrMC ]
           
    ; replicate gc.type for trMC
    trMC_orig_type = ( galcat.type[ replicate_var(galcat.trMC_cc) ] )
    
    ; if we are dealing with a run without BHs, add a type of -1 that will never match
    galcat_types = galcat.types
    if sP.gfmBHs eq 0 then $
      galcat_types = mod_struct( galcat_types, 'bhs', -1 )
              
    ; make association of tracer IDs between the two galaxy catalogs
    match,galcat.trMC_ids,galcatTarget.trMC_ids,trOrig_inds,trTarget_inds,count=countMatch
    sort_inds = sort(trOrig_inds)
    trOrig_inds = trOrig_inds[sort_inds]
    trTarget_inds = trTarget_inds[sort_inds]
    sort_inds = !NULL
    
    frac = float(countMatch)*100 / nTotSearch
    print,' Matched fraction: '+string(frac,format='(f4.1)')+'% ('+$
      str(countMatch)+' of '+str(nTotSearch)+')'
      
    ; mapping of tracers from galcat.trMC to galcatTarget.trMC
    trMapInds = lonarr(nTotSearch) - 1
    trMapInds[trOrig_inds] = trTarget_inds
    
    ; mapping of tracers from galcatTarget.trMC to galcat.trMC
    trMapIndsInv = lonarr(nTotSearchTarg) - 1
    trMapIndsInv[trTarget_inds] = trOrig_inds
    
    ; replicate to each trMC at targetSnap, and make mask for -matched- between the two galaxyCats
    ; we will use these unmatched trMCs @ targetSnap to mark outflow to parent indices @ origSnap
    origParIDsTarg = hIndsInv[ mtTarget.gcIndOrigTrMC ]
    
    targTrMatchedMask = bytarr(nTotSearchTarg)
    targTrMatchedMask[trTarget_inds] = 1B
    
    ; make replicated tracer parent indices at origSnap, then transform this to be ordered as
    ; one origSnapParent for each targetSnap tracer
    trParIndRep = galCatRepParentIDs(galcat=galcat,child_counts=galcat.trMC_cc)
    
    galcatIndOrigTrMC = lonarr(nTotSearchTarg) - 2 ; default -2, see below
    
    w = where(trMapIndsInv ne -1,count)
    if count eq 0 then message,'Error'
    galcatIndOrigTrMC[w] = trParIndRep[ trMapIndsInv[w] ]    
    
    trParIndRep = !NULL
    trOrig_inds = !NULL
    trTarget_inds = !NULL
    
    ; make replicated tracer parent indices at targetSnap, then transform this to be ordered as
    ; one targetSnapParent for each origSnap tracer
    trTargParIndRep = galCatRepParentIDs(galcat=galcatTarget,child_counts=galcatTarget.trMC_cc)
    
    galcatTargetIndOrigTrMC = lonarr(nTotSearch) - 2 ; default -2, to prevent match to -1 below
    
    w = where(trMapInds ne -1,count)
    if count eq 0 then message,'Error'
    galcatTargetIndOrigTrMC[w] = trTargParIndRep[ trMapInds[w] ]
    
    trTargParIndRep = !NULL
    
    ; replicate gcTarget.type for trMC, then transform this to be ordered as one type for each
    ; origSnap tracer, with type -2 for unmatched
    trTargTypeRep = ( galcatTarget.type[ replicate_var(galcatTarget.trMC_cc) ] )
    
    trMC_targ_type = intarr(nTotSearch) - 2
    
    w = where(trMapInds ne -1,count)
    if count eq 0 then message,'Error'
    trMC_targ_type[w] = trTargTypeRep[ trMapInds[w] ]

    ; start: store the main arrays as a structure so we can write them directly for all tracers
    r = { accFlag           : intarr(nTotSearch)-2        ,$
          accFlagTarg       : intarr(nTotSearchTarg)-2    ,$
          origParIDsTarg    : origParIDsTarg              ,$
          targTrMatchedMask : targTrMatchedMask            }
          
    targTrMatchedMask = !NULL
    origParIDsTarg = !NULL
    
    trMaskDebug = intarr(nTotSearch) ; DEBUG
    trMaskDebugTarg = intarr(nTotSearchTarg) ; DEBUG
    
    ; which of these tracers are also in galcatTarget?
    wMap = where(trMapInds ge 0,countMap,comp=wNoMap,ncomp=countNoMap)
    
    if countMap eq 0 or countNoMap eq 0 then message,'Error: Unlikely, and below will fail.'
    ; (1a). those which did not map:
      ; flag (1) any gal/stars/bh tracers from origSnap but not in the progenitor at targetSnap
      tr_types_nomap = trMC_orig_type[ wNoMap ]
        
      w = where(tr_types_nomap eq galcat_types.gal or $
                tr_types_nomap eq galcat_types.stars or $
                tr_types_nomap eq galcat_types.bhs, count_galAcc_noMap)
                  
      if count_galAcc_noMap eq 0 then message,'Error'
      r.accFlag[wNoMap[w]] = 1
      trMaskDebug[wNoMap[w]] += 1
      
      ; flag (2) any gmem tracers from origSnap but not in the progenitor at targetSnap
      ; as net halo inflow (tracers in inter at origSnap do not contribute to anything)
      w = where(tr_types_nomap eq galcat_types.gmem, count_gmemAcc_noMap)
        
      if count_gmemAcc_noMap eq 0 then message,'Error'
      r.accFlag[wNoMap[w]] = 2
      trMaskDebug[wNoMap[w]] += 1
    
    ; (1b). those which did map:
    
      ; (2a). flag (-1) any tracers of halos with failed mergerTree connectivity, do not include these halos
      w = where(origParIDs[wMap] eq -1,count,comp=wTracked)
      frac = float(count)*100 / nTotSearch
      print,' Failed connectivity: '+string(frac,format='(f4.1)')+'% ('+$
      str(count)+' of '+str(nTotSearch)+')'
    
      r.accFlag[wMap[w]] = -1
      trMaskDebug[wMap[w]] += 1
      
      ; (2b):
      wTracked = wMap[wTracked] ; index galcat.trMC_ids directly
      tr_types_orig  = trMC_orig_type[ wTracked ]
      tr_types_targ  = trMC_targ_type[ wTracked ]
      tr_parind_orig = origParIDs[wTracked]
      tr_parind_targ = galcatTargetIndOrigTrMC[wTracked]
      
        ; (3a1): flag (1) any gal/stars/bh tracers from origSnap in the progenitor 
        ; but in the inter/gmem type, as net galaxy inflow
        w = where((tr_types_orig eq galcat_types.gal   or $
                   tr_types_orig eq galcat_types.stars or $
                   tr_types_orig eq galcat_types.bhs   ) and $
                  (tr_types_targ eq galcat_types.inter or $
                   tr_types_targ eq galcat_types.gmem)   and $
                  tr_parind_orig eq tr_parind_targ, count_galAcc_sameHalo)
                  
        if count_galAcc_sameHalo eq 0 then message,'Error'
        r.accFlag[wTracked[w]] = 1
        trMaskDebug[wTracked[w]] += 1
      
        ; (3b): in different halo, flag (1) or (2) based on gal/halo type at origSnap
        w = where((tr_types_orig eq galcat_types.gal   or $
                   tr_types_orig eq galcat_types.stars or $
                   tr_types_orig eq galcat_types.bhs   ) and $
                  tr_parind_orig ne tr_parind_targ, count_galAcc_diffHalo)
                  
        if count_galAcc_diffHalo eq 0 then message,'Error'
        r.accFlag[wTracked[w]] = 1
        trMaskDebug[wTracked[w]] += 1
        
        w = where(tr_types_orig eq galcat_types.gmem and $
                  tr_parind_orig ne tr_parind_targ, count_haloAcc_diffHalo)
                  
        if count_haloAcc_diffHalo eq 0 then message,'Error'
        r.accFlag[wTracked[w]] = 2
        trMaskDebug[wTracked[w]] += 1
      
    ; fill in the cases where we have no contribution as a verification
        ; (1a) inter tracers not in galcatTarget do not contribute to anything
        w = where(tr_types_nomap eq galcat_types.inter, count)
        
        if count eq 0 then message,'Error'
        r.accFlag[wNoMap[w]] = 0
        trMaskDebug[wNoMap[w]] += 1
      
        ; (3a2+3a4): in same halo, was inter/gmem in origSnap and in targetSnap
        w = where((tr_types_orig eq galcat_types.inter   or $
                   tr_types_orig eq galcat_types.gmem ) and $
                  tr_parind_orig eq tr_parind_targ, count)
                  
        if count eq 0 then message,'Error'
        r.accFlag[wTracked[w]] = 0
        trMaskDebug[wTracked[w]] += 1
        
        ; (3a3): in samehalo, was gal/stars/bh in origSnap and in targetSnap
        w = where((tr_types_orig eq galcat_types.gal   or $
                   tr_types_orig eq galcat_types.stars or $
                   tr_types_orig eq galcat_types.bhs   ) and $
                  (tr_types_targ eq galcat_types.gal or $
                   tr_types_targ eq galcat_types.stars or $
                   tr_types_targ eq galcat_types.bhs   ) and $
                  tr_parind_orig eq tr_parind_targ, count)
                  
        if count eq 0 then message,'Error'
        r.accFlag[wTracked[w]] = 0
        trMaskDebug[wTracked[w]] += 1
        
        ; (3b): in different halo, inter in origSnap
        w = where(tr_types_orig eq galcat_types.inter and $
                  tr_parind_orig ne tr_parind_targ, count)
                  
        if count eq 0 then message,'Error'
        r.accFlag[wTracked[w]] = 0
        trMaskDebug[wTracked[w]] += 1
        
    ; consider the material in the (tracked) progenitor of each halo at targetSnap
    ; which does not appear in the origSnap galaxyCat
      ; (4a): galaxy outflow
      w = where(r.targTrMatchedMask eq 0B and $
                (trTargTypeRep eq galcat_types.gal or $
                 trTargTypeRep eq galcat_types.stars or $
                 trTargTypeRep eq galcat_types.bhs),count)
      if count eq 0 then message,'Error'

      r.accFlagTarg[w] = 1
      trMaskDebugTarg[w] += 1
      
      ; (4b): halo outflow
      w = where(r.targTrMatchedMask eq 0B and $
                trTargTypeRep eq galcat_types.gmem,count)
      if count eq 0 then message,'Error'
      
      r.accFlagTarg[w] = 2
      trMaskDebugTarg[w] += 1
      
      ; (4c): inter, not unused, just flag
      w = where(r.targTrMatchedMask eq 0B and $
                trTargTypeRep eq galcat_types.inter,count)
      if count eq 0 then message,'Error'
      
      r.accFlagTarg[w] = 0
      trMaskDebugTarg[w] += 1
      
    ; consider the material in the (tracked) progenitor of each halo at targetSnap
    ; which is in the origSnap galaxyCat, but in a different descendent
      ; (4b): galaxy outflow
      w = where(r.targTrMatchedMask eq 1B and $
                (trTargTypeRep eq galcat_types.gal or $
                 trTargTypeRep eq galcat_types.stars or $
                 trTargTypeRep eq galcat_types.bhs) and $
                 galcatIndOrigTrMC ne r.origParIDsTarg,count)
      if count eq 0 then message,'Error'

      r.accFlagTarg[w] = 1
      trMaskDebugTarg[w] += 1
      
      ; (5b): halo outflow
      w = where(r.targTrMatchedMask eq 1B and $
                trTargTypeRep eq galcat_types.gmem and $
                galcatIndOrigTrMC ne r.origParIDsTarg,count)
      if count eq 0 then message,'Error'

      r.accFlagTarg[w] = 2
      trMaskDebugTarg[w] += 1
      
      ; (5c): 
      w = where(r.targTrMatchedMask eq 1B and $
                trTargTypeRep eq galcat_types.inter and $
                galcatIndOrigTrMC ne r.origParIDsTarg,count)
      if count eq 0 then message,'Error'

      r.accFlagTarg[w] = 0
      trMaskDebugTarg[w] += 1
      
    ; matched material in targetSnap in the descendent (just mark, does not contribute)
      ; galaxy/halo/inter
      w = where(r.targTrMatchedMask eq 1B and $
                galcatIndOrigTrMC eq r.origParIDsTarg,count)
      if count eq 0 then message,'Error'
      
      r.accFlagTarg[w] = 0
      trMaskDebugTarg[w] += 1
        
    ; ambiguities:
    ; 1. what about material in gmem at origSnap and in gal/stars/bh/inter at targetSnap?
    ;    this is somehow net inflow into the halo, from within? (do not include as halo acc)
    
    w = where(trMaskDebug ne 1,count)
    if count gt 0 then message,'Error: Some tracers not considered.'
          
    w = where(trMaskDebugTarg ne 1,count)
    if count gt 0 then message,'Error: Some targ tracers not considered.'
          
    w = where(r.accFlag lt -1 or r.accFlag gt 2,count)
    if count gt 0 then message,'Error: Strange accFlag values.'
    
    w = where(r.accFlagTarg lt -1 or r.accFlag gt 2,count)
    if count gt 0 then message,'Error: Strange accFlagTarg values.'

    ; save
    save,r,filename=saveFilename
    print,'Saved: '+strmid(saveFilename,strlen(sp.derivPath))
    
  endif
  
  ; VELOCITY TRACERS case - will be similar to above since there could be multiple
  ; ---------------------
  if sP.trMCPerCell eq -1 then begin
    message,'TODO'                
  endif
  
  sP.snap = origSnap ; restore sP.snap
  return,r

end
