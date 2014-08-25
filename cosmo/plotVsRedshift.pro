; plotVsRedshift.pro
; feedback - plots skipping tconst/tvircur/tviracc definitions and timewindow variations
;   in favor of redshift evolution and redshift panels
; dnelson jul.2014

; plotsVsRedshiftBin(): bin and save values for below plotting
    
pro plotsVsRedshiftBin, tempFilename, runs, redshifts, res, $
                             timeWindow, tVirInd, accModes, accRateModel
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  foreach run,runs,i do begin
    
    sP_z  = {}
    mbv_z = {}
    nModes = 0
  
    ; loop over each redshift
    foreach redshift,redshifts,j do begin
    
      sP_loc = simParams(res=res,run=run,redshift=redshift)
      sP_loc.accRateModel = accRateModel
      
      sP_mode  = {}
      mbv_mode = {}
      
      if run eq 'gadget' and redshift ne 2.0 then continue ; only want z=2 for gadget
    
      ; load halo masses to re-bin median lines of some quantities
      gc = loadGroupCat(sP=sP_loc,/skipIDs)
      gcMasses = codeMassToLogMsun(gc.subgroupMass)
      gcW = where(gc.subgroupMass ge logMsunToCodeMass(8.9)*units.hubbleParam,count_gcW)
    
      ; loop over each accretion mode
      foreach accMode,accModes,k do begin
        if accMode eq 'recycled' and sP_loc.gfmWinds eq 0 then continue ; skip recycled for nonWinds 

        mbv_loc = haloMassBinValues(sP=sP_loc,accMode=accMode,timeWindow=timeWindow,/search)
      
        if n_tags(mbv_loc) eq 0 then begin
          print,run,' ',redshift,' ',accMode,' (SKIP)'
          continue
        endif
        
        ; OLD: empty structure to hold ratios
        numMassBins = n_elements(mbv_loc.logMassBinCen)
        
        mode_ratio   = { gal_hot   : fltarr(numMassBins) + !values.f_nan ,$
                         gal_cold  : fltarr(numMassBins) + !values.f_nan ,$
                         halo_hot  : fltarr(numMassBins) + !values.f_nan ,$
                         halo_cold : fltarr(numMassBins) + !values.f_nan  } 
        specific_net = { gal_hot   : fltarr(numMassBins) + !values.f_nan ,$
                         gal_cold  : fltarr(numMassBins) + !values.f_nan ,$
                         halo_hot  : fltarr(numMassBins) + !values.f_nan ,$
                         halo_cold : fltarr(numMassBins) + !values.f_nan  }
        coldfrac_net = { gal  : fltarr(numMassBins) + !values.f_nan ,$
                         halo : fltarr(numMassBins) + !values.f_nan  }
                         
        ; NEW: empty structures to hold galaxy_acc,galaxy_out,halo_acc,halo_out
        aoTemplate = { hot   : fltarr(numMassBins,3) + !values.f_nan ,$
                       cold  : fltarr(numMassBins,3) + !values.f_nan ,$
                       total : fltarr(numMassBins,3) + !values.f_nan}
                         
        galaxyRates = { $
          galaxy_acc    : aoTemplate ,$
          galaxy_out    : aoTemplate ,$
          halo_acc      : aoTemplate ,$
          halo_out      : aoTemplate ,$
          yvals_gal_out : aoTemplate ,$ ; (#2+#3) galaxy only
          yvals_gal_acc : aoTemplate ,$ ; (#1+#2) galaxy only
          yvals_halo    : aoTemplate ,$ ; halo only (halo_haloacc+galaxy_haloacc)
          abcd_gal_out  : aoTemplate ,$ ; A+B+C = 3+2
          abcd_gal_acc  : aoTemplate ,$ ; A+C+D = 1+4
          abcd_gal_net  : aoTemplate  $ ; B-D = 1+4-2-3
        }
        
        ; NEW: empty structures to hold updated mode ratios and fractions
        modeRatiosNet    = { galaxy : aoTemplate ,$
                             halo   : aoTemplate  }
        coldFractionsNet = { galaxy : fltarr(numMassBins,3) + !values.f_nan ,$
                             halo   : fltarr(numMassBins,3) + !values.f_nan  }
        coldFractionsAcc = { galaxy : fltarr(numMassBins,3) + !values.f_nan ,$
                             halo   : fltarr(numMassBins,3) + !values.f_nan  }     

        ; ABCD: mode ratios and fractions
        modeRatiosABCD = { galaxy : aoTemplate ,$
                           halo   : aoTemplate  }
        coldFractionsNetABCD = { galaxy : fltarr(numMassBins,3) + !values.f_nan ,$
                                 halo   : fltarr(numMassBins,3) + !values.f_nan  }
        coldFractionsAccABCD = { galaxy : fltarr(numMassBins,3) + !values.f_nan ,$
                                 halo   : fltarr(numMassBins,3) + !values.f_nan  }  

        ; OLD: empty structures to hold (acc/out/net) percentiles
        galaxyPerc = { $
          acc : { hot  : fltarr(numMassBins,3) + !values.f_nan ,$
                  cold : fltarr(numMassBins,3) + !values.f_nan  } ,$
          out : { hot  : fltarr(numMassBins,3) + !values.f_nan ,$
                  cold : fltarr(numMassBins,3) + !values.f_nan  } ,$  
          net : { hot  : fltarr(numMassBins,3) + !values.f_nan ,$
                  cold : fltarr(numMassBins,3) + !values.f_nan  } $
        }
        
        ; extract fields we care about
        mbv_locAll = mbv_loc
        mbv_loc = {galaxyMedian   : mbv_loc.galaxyMedian  ,$
                   haloMedian     : mbv_loc.haloMedian    ,$
                   logMassBinCen  : mbv_loc.logMassBinCen ,$
                   mode_ratio     : mode_ratio            ,$ ; old
                   specific_net   : specific_net          ,$ ; old
                   coldfrac_net   : coldfrac_net          ,$ ; old
                   galaxyRates    : galaxyRates           ,$ ; in/out of galaxy
                   haloRates      : galaxyRates           ,$ ; in/out of halo
                   modeRatiosNet  : modeRatiosNet         ,$ ; ratios to all
                   coldFractionsNet : coldFractionsNet    ,$ ; ratios of cold/(hot+cold)
                   coldFractionsAcc : coldFractionsAcc    ,$
                   galaxyPerc     : galaxyPerc            ,$ ; old
                   haloPerc       : galaxyPerc             } ; old
                   
        ; eta_w comparison (mode='all' only)
        if k eq 0 then begin
          mbv_loc = mod_struct( mbv_loc, 'eta_w_binned',   fltarr(numMassBins,3) )
          mbv_loc = mod_struct( mbv_loc, 'eta_w_xvals',    fltarr(count_gcW) + !values.f_nan )
          mbv_loc = mod_struct( mbv_loc, 'eta_w_unbinned', fltarr(count_gcW) + !values.f_nan )
          
          if redshift le 1.0 then begin
            mbv_loc = mod_struct( mbv_loc, 'eta_sfr2_net', fltarr(count_gcW) + !values.f_nan )
            mbv_loc = mod_struct( mbv_loc, 'eta_fit',      fltarr(2) )
          endif
          
          ; ABCD approach
          mbv_loc = mod_struct( mbv_loc, 'eta_abcd_binned',    fltarr(numMassBins,3) )
          mbv_loc = mod_struct( mbv_loc, 'eta_abcd2_binned',   fltarr(numMassBins,3) )
          mbv_loc = mod_struct( mbv_loc, 'eta_abcd_unbinned',  fltarr(count_gcW) + !values.f_nan )
          mbv_loc = mod_struct( mbv_loc, 'eta_abcd2_unbinned', fltarr(count_gcW) + !values.f_nan )
          
          mbv_loc = mod_struct( mbv_loc, 'eta_sfr2_net_abcd', fltarr(count_gcW) + !values.f_nan )
          mbv_loc = mod_struct( mbv_loc, 'eta_sfr2_net_abcd2', fltarr(count_gcW) + !values.f_nan )
        endif
        
        print,run,' ',redshift,' ',accMode
        sP_mode  = mod_struct(sP_mode, 'mode'+str(k)+'_'+accMode, sP_loc)
        mbv_mode = mod_struct(mbv_mode, 'mode'+str(k)+'_'+accMode, mbv_loc)
        mbv_modeAll = mod_struct(mbv_modeAll, 'mode'+str(k)+'_'+accMode, mbv_locAll) ; not saved
      endforeach ; accModes,k
    
      if nModes eq 0 then nModes = n_tags(mbv_mode) ; set
      
      if n_tags(sP_mode) eq 0 or n_tags(sP_mode) lt nModes then begin
        print,'SKIP REDSHIFT: ',redshift,' ',run,' (missing mode)'
        continue
      endif
      
      ; truncate gcMasses to be the same as from mbv
      maxHist = n_elements(mbv_locAll.galaxy.netRate.total.hot.num)
      gcMasses = gcMasses[0:maxHist-1]
      
      for k=0,nModes-1 do begin
        
        ; old
        ; ---
        ratio_gal_hot   = reform(mbv_modeAll.(k).galaxy.netRate.total.hot.tVirAcc[tVirInd,*] / $
                                 mbv_modeAll.(0).galaxy.netRate.total.hot.tVirAcc[tVirInd,*])
        ratio_gal_cold  = reform(mbv_modeAll.(k).galaxy.netRate.total.cold.tVirAcc[tVirInd,*] / $
                                 mbv_modeAll.(0).galaxy.netRate.total.cold.tVirAcc[tVirInd,*])
        ratio_halo_hot  = reform(mbv_modeAll.(k).halo.netRate.total.hot.tVirAcc[tVirInd,*] / $
                                 mbv_modeAll.(0).halo.netRate.total.hot.tVirAcc[tVirInd,*])
        ratio_halo_cold = reform(mbv_modeAll.(k).halo.netRate.total.cold.tVirAcc[tVirInd,*] / $
                                 mbv_modeAll.(0).halo.netRate.total.cold.tVirAcc[tVirInd,*])

        gal_hot   = reform(mbv_modeAll.(k).galaxy.netRate.total.hot.tVirAcc[tVirInd,*]) / $
                    10.0^gcMasses * 1e9
        gal_cold  = reform(mbv_modeAll.(k).galaxy.netRate.total.cold.tVirAcc[tVirInd,*]) / $
                    10.0^gcMasses * 1e9
        halo_hot  = reform(mbv_modeAll.(k).halo.netRate.total.hot.tVirAcc[tVirInd,*]) / $
                    10.0^gcMasses * 1e9
        halo_cold = reform(mbv_modeAll.(k).halo.netRate.total.cold.tVirAcc[tVirInd,*]) / $
                    10.0^gcMasses * 1e9
                    
        ; calculate galaxy_galacc,galaxy_galout,halo_galacc,halo_galout:
        ; ---------
        ; galaxy_galacc (#1): rate into galaxy for material now in the galaxy
        ; since any material that left is not included, this is a net inflow rate
        galaxy_galacc_hot = $
          mbv_modeAll.(k).galaxy.accRate.gal.hot.tVirAcc[0,*] + $
          mbv_modeAll.(k).galaxy.accRate.stars.hot.tVirAcc[0,*]
        galaxy_galacc_cold = $
          mbv_modeAll.(k).galaxy.accRate.gal.cold.tVirAcc[0,*] + $
          mbv_modeAll.(k).galaxy.accRate.stars.cold.tVirAcc[0,*]
          
        if sP_loc.gfmWinds then begin
          galaxy_galacc_hot  += mbv_modeAll.(k).galaxy.accRate.bhs.hot.tVirAcc[0,*]
          galaxy_galacc_cold += mbv_modeAll.(k).galaxy.accRate.bhs.cold.tVirAcc[0,*]
        endif
            
        ; halo_galacc (#4): rate into galaxy for material now in the halo
        ; if it went in, but is outside, then it must have come back out
        ; (it is counted already in halo_out)
        ; the outflow rate of material which entered the galaxy within the last TW
        halo_galacc_hot = $
          mbv_modeAll.(k).galaxy.accRate.inter.hot.tVirAcc[0,*] + $
          mbv_modeAll.(k).galaxy.accRate.gmem.hot.tVirAcc[0,*]
          
        halo_galacc_cold = $
          mbv_modeAll.(k).galaxy.accRate.inter.cold.tVirAcc[0,*] + $
          mbv_modeAll.(k).galaxy.accRate.gmem.cold.tVirAcc[0,*]
          
        ; halo_galout (#2): rate out of galaxy for material now in the halo
        ; should be the sum of halo_acc and the outflow of material which entered
        ;   the galaxy more than TW ago
        ; halo_galout-halo_galacc=outflow from old material
        halo_galout_hot = $
          mbv_modeAll.(k).galaxy.outRate.inter.hot.tVirAcc[0,*] + $
          mbv_modeAll.(k).galaxy.outRate.gmem.hot.tVirAcc[0,*]

        halo_galout_cold = $
          mbv_modeAll.(k).galaxy.outRate.inter.cold.tVirAcc[0,*] + $
          mbv_modeAll.(k).galaxy.outRate.gmem.cold.tVirAcc[0,*]
          
        ; galaxy_galout (#3): rate out of galaxy for material now in the galaxy
        ; if it went out, but is inside, then it must have come back
        ; (its accretion is counted already in galaxy_acc)
        ; galaxy_galout/halo_galout = ratio of material that gas exited and returned within TW
        galaxy_galout_hot = $
          mbv_modeAll.(k).galaxy.outRate.gal.hot.tVirAcc[0,*] + $
          mbv_modeAll.(k).galaxy.outRate.stars.hot.tVirAcc[0,*]
          
        galaxy_galout_cold = $
          mbv_modeAll.(k).galaxy.outRate.gal.cold.tVirAcc[0,*] + $
          mbv_modeAll.(k).galaxy.outRate.stars.cold.tVirAcc[0,*]
          
        if sP_loc.gfmWinds then begin
          galaxy_galout_hot  += mbv_modeAll.(k).galaxy.outRate.bhs.hot.tVirAcc[0,*]
          galaxy_galout_cold += mbv_modeAll.(k).galaxy.outRate.bhs.cold.tVirAcc[0,*]
        endif
        
        ; calculate galaxy_haloacc,galaxy_haloout,halo_haloacc,halo_haloout:
        ; ---------
        ; galaxy_haloacc: rate into halo for material now in the galaxy
        galaxy_haloacc_hot = $
          mbv_modeAll.(k).halo.accRate.gal.hot.tVirAcc[0,*] + $
          mbv_modeAll.(k).halo.accRate.stars.hot.tVirAcc[0,*]
        galaxy_haloacc_cold = $
          mbv_modeAll.(k).halo.accRate.gal.cold.tVirAcc[0,*] + $
          mbv_modeAll.(k).halo.accRate.stars.cold.tVirAcc[0,*]
          
        if sP_loc.gfmWinds then begin
          galaxy_haloacc_hot  += mbv_modeAll.(k).halo.accRate.bhs.hot.tVirAcc[0,*]
          galaxy_haloacc_cold += mbv_modeAll.(k).halo.accRate.bhs.cold.tVirAcc[0,*]
        endif
            
        ; halo_haloacc: rate into halo for material now in the halo
        ; if it went in, but isn't in the halo any longer, it is either in the galaxy or 
        ; outside the halo. so this is the net accretion rate into the halo, -excluding- 
        ; material that makes it into the galaxy
        halo_haloacc_hot = $
          mbv_modeAll.(k).halo.accRate.inter.hot.tVirAcc[0,*] + $
          mbv_modeAll.(k).halo.accRate.gmem.hot.tVirAcc[0,*]
          
        halo_haloacc_cold = $
          mbv_modeAll.(k).halo.accRate.inter.cold.tVirAcc[0,*] + $
          mbv_modeAll.(k).halo.accRate.gmem.cold.tVirAcc[0,*]
          
        ; halo_haloout: rate out of halo for material now in the halo
        ; if it went out, but is now inside, then is must have come back (counted in halo_haloacc)
        ; halo_haloout/(halo_haloacc+halo_haloout)=fraction of material in halo that has
        ; recycled in the past TW
        halo_haloout_hot = $
          mbv_modeAll.(k).halo.outRate.inter.hot.tVirAcc[0,*] + $
          mbv_modeAll.(k).halo.outRate.gmem.hot.tVirAcc[0,*]

        halo_haloout_cold = $
          mbv_modeAll.(k).halo.outRate.inter.cold.tVirAcc[0,*] + $
          mbv_modeAll.(k).halo.outRate.gmem.cold.tVirAcc[0,*]
          
        ; galaxy_haloout: rate out of halo for material now in the galaxy
        ; if it went out, but is inside, then it must have come back
        ; (its accretion is counted already in galaxy_haloacc)
        ; this should be quite small
        galaxy_haloout_hot = $
          mbv_modeAll.(k).halo.outRate.gal.hot.tVirAcc[0,*] + $
          mbv_modeAll.(k).halo.outRate.stars.hot.tVirAcc[0,*]
          
        galaxy_haloout_cold = $
          mbv_modeAll.(k).halo.outRate.gal.cold.tVirAcc[0,*] + $
          mbv_modeAll.(k).halo.outRate.stars.cold.tVirAcc[0,*]
          
        if sP_loc.gfmWinds then begin
          galaxy_haloout_hot  += mbv_modeAll.(k).halo.outRate.bhs.hot.tVirAcc[0,*]
          galaxy_haloout_cold += mbv_modeAll.(k).halo.outRate.bhs.cold.tVirAcc[0,*]
        endif
        
        ; updated mode ratios (for net rates)
        ; galaxy net = galaxy_galacc_hot+galaxy_galacc_cold
        ; halo net = halo_haloacc_hot+halo_haloacc_cold + galaxy_haloacc_hot+galaxy_haloacc_cold
        ; --------------
        
        ; save mode0 for mode ratios if k=0, otherwise compute mode ratios
        if k eq 0 then begin
          mode0_saved = { galaxy_galacc_hot  : galaxy_galacc_hot  ,$
                          galaxy_galacc_cold : galaxy_galacc_cold ,$
                          galaxy_galout_hot  : galaxy_galout_hot  ,$
                          galaxy_galout_cold : galaxy_galout_cold ,$
                          halo_galacc_hot  : halo_galacc_hot  ,$
                          halo_galacc_cold : halo_galacc_cold ,$
                          halo_galout_hot  : halo_galout_hot  ,$
                          halo_galout_cold : halo_galout_cold ,$
                          galaxy_haloacc_hot  : galaxy_haloacc_hot  ,$
                          galaxy_haloacc_cold : galaxy_haloacc_cold ,$
                          galaxy_haloout_hot  : galaxy_haloout_hot  ,$
                          galaxy_haloout_cold : galaxy_haloout_cold ,$
                          halo_haloacc_hot  : halo_haloacc_hot  ,$
                          halo_haloacc_cold : halo_haloacc_cold ,$
                          halo_haloout_hot  : halo_haloout_hot  ,$
                          halo_haloout_cold : halo_haloout_cold  }
        endif else begin
          netratio_gal_hot   = galaxy_galacc_hot / mode0_saved.galaxy_galacc_hot
          netratio_gal_cold  = galaxy_galacc_cold / mode0_saved.galaxy_galacc_cold
          netratio_gal_total = (galaxy_galacc_hot + galaxy_galacc_cold) / $
                               (mode0_saved.galaxy_galacc_hot + mode0_saved.galaxy_galacc_cold)
                               
          netratio_halo_hot   = (halo_haloacc_hot + galaxy_haloacc_hot) / $
                                (mode0_saved.halo_haloacc_hot + mode0_saved.galaxy_haloacc_hot)
          netratio_halo_cold  = (halo_haloacc_cold + galaxy_haloacc_cold) / $
                                (mode0_saved.halo_haloacc_cold + mode0_saved.galaxy_haloacc_cold)
          netratio_halo_total = (halo_haloacc_hot + galaxy_haloacc_hot + $
                                 halo_haloacc_cold + galaxy_haloacc_cold) / $
                                (mode0_saved.halo_haloacc_hot + mode0_saved.galaxy_haloacc_hot + $
                                 mode0_saved.halo_haloacc_cold + mode0_saved.galaxy_haloacc_cold)

          ; ABCD method (TODO)
          abcd_ratio_gal_hot = 0
        endelse
                
        ; wind_eta comparison (for mode='all' only)
        ; --------
        if k eq 0 then begin
        
          ; OLD: definitions
          eta_out = galaxy_galout_hot + galaxy_galout_cold + halo_galout_hot + halo_galout_cold ;#2+#3
          eta_net = galaxy_galacc_hot + galaxy_galacc_cold ;#1
          
          ; calculate correction factor
          if run eq 'tracer' and redshift gt 1.0 then begin
            ; for tracer z>1, group catalog has missing SFRs, just use z=1 values
            fit = [0.19466,-0.043346]
          endif else begin
            ; normal
            sfr1 = gc.subgroupSfr[0:maxHist-1]
            sfr2 = gc.subgroupSfrInRad[0:maxHist-1]
            
            yy = alog10(sfr2/eta_net)
            w=where(gcMasses ge 10.0 and gcMasses le 13.0 and yy ge -4.0 and yy lt 2.0)
            fit = linfit(gcMasses[w],yy[w])
          endelse
          
          ; un-binned
          eta_w_unbinned = eta_out / eta_net
          
          ; apply correction (using linear fit)
          sfr_acc_correction_factor = 10.0^(gcMasses * fit[1] + fit[0])
          eta_w_unbinned /= sfr_acc_correction_factor
          
          ; save unbinned values (only those in halo mass >8.9 visible in plots)
          xvals = codeMassToLogMsun( logMsunToCodeMass( gcMasses ) / units.hubbleParam )
                      
          mbv_mode.(k).eta_w_xvals    = xvals[gcW]
          mbv_mode.(k).eta_w_unbinned = eta_w_unbinned[gcW]
          
          if redshift le 1.0 then begin
            mbv_mode.(k).eta_sfr2_net = yy[gcW]
            mbv_mode.(k).eta_fit      = fit
          endif
          
          ; NEW: ABCD approach
          eta_out = galaxy_galout_hot + galaxy_galout_cold + halo_galout_hot + halo_galout_cold ;#2+#3
          eta_acc = galaxy_galacc_hot + galaxy_galacc_cold + halo_galacc_hot + halo_galacc_cold ;#1+#4
          eta_net = galaxy_galacc_hot + galaxy_galacc_cold + halo_galacc_hot + halo_galacc_cold $
                  - (halo_galout_hot + halo_galout_cold + galaxy_galout_hot + galaxy_galout_cold) 
                  ;#1+#4-#2-#3
          
          eta_abcd_unbinned = eta_out / eta_net
          eta_abcd2_unbinned = eta_out / eta_acc
          
          mbv_mode.(k).eta_abcd_unbinned  = eta_abcd_unbinned[gcW]
          mbv_mode.(k).eta_abcd2_unbinned = eta_abcd2_unbinned[gcW]
          
          if redshift le 1.0 then begin
            mbv_mode.(k).eta_sfr2_net_abcd = (alog10(sfr2/eta_net))[gcW]
            mbv_mode.(k).eta_sfr2_net_abcd2 = (alog10(sfr2/eta_acc))[gcW]
          endif
        endif ;wind_eta
            
        ; rebin in halo mass
        for m=0,n_elements(mbv_loc.logMassBinCen)-1 do begin
          ; restrict medians to primary only (num>0)
          w = where(gcMasses gt mbv_modeAll.(k).logMassBins[m] and $
                    gcMasses le mbv_modeAll.(k).logMassBins[m+1] and $
                    (mbv_modeAll.(0).galaxy.coldFrac.total.num+$
                     mbv_modeAll.(0).halo.coldFrac.total.num) gt 0,count)
                    
          if count eq 0 then continue
          if max(w) ge maxHist then message,'Error'
          
          ; (1) OLD: mode_ratio: bin ratios of modes vs 'all' mode
          mbv_mode.(k).mode_ratio.gal_hot[m]   = median( ratio_gal_hot[w] )
          mbv_mode.(k).mode_ratio.gal_cold[m]  = median( ratio_gal_cold[w] )
          mbv_mode.(k).mode_ratio.halo_hot[m]  = median( ratio_halo_hot[w] )
          mbv_mode.(k).mode_ratio.halo_cold[m] = median( ratio_halo_cold[w] )
          
          ; (2) OLD: specific_net and coldfrac_net: specific net rates/coldfracs, galaxy by galaxy (yr -> Gyr)
          mbv_mode.(k).specific_net.gal_hot[m]   = median( gal_hot[w] )
          mbv_mode.(k).specific_net.gal_cold[m]  = median( gal_cold[w] )
          mbv_mode.(k).specific_net.halo_hot[m]  = median( halo_hot[w] )
          mbv_mode.(k).specific_net.halo_cold[m] = median( halo_cold[w] )
          
          mbv_mode.(k).coldfrac_net.gal[m]  = median( gal_cold[w] / (gal_hot[w]+gal_cold[w]) )
          mbv_mode.(k).coldfrac_net.halo[m] = median( halo_cold[w] / (halo_hot[w]+halo_cold[w]) )
          
          ; (3) NEW: galaxy_acc,galaxy_out,halo_acc,halo_out
          mbv_mode.(k).galaxyRates.galaxy_acc.hot[m,*] = percentiles( galaxy_galacc_hot[w] )
          mbv_mode.(k).galaxyRates.galaxy_out.hot[m,*] = percentiles( galaxy_galout_hot[w] )
          mbv_mode.(k).galaxyRates.halo_acc.hot[m,*]   = percentiles( halo_galacc_hot[w] )
          mbv_mode.(k).galaxyRates.halo_out.hot[m,*]   = percentiles( halo_galout_hot[w] )
          
          mbv_mode.(k).galaxyRates.galaxy_acc.cold[m,*] = percentiles( galaxy_galacc_cold[w] )
          mbv_mode.(k).galaxyRates.galaxy_out.cold[m,*] = percentiles( galaxy_galout_cold[w] )
          mbv_mode.(k).galaxyRates.halo_acc.cold[m,*]   = percentiles( halo_galacc_cold[w] )
          mbv_mode.(k).galaxyRates.halo_out.cold[m,*]   = percentiles( halo_galout_cold[w] )
          
          mbv_mode.(k).galaxyRates.galaxy_acc.total[m,*] = $
            percentiles( galaxy_galacc_hot[w]+galaxy_galacc_cold[w] )
          mbv_mode.(k).galaxyRates.galaxy_out.total[m,*] = $
            percentiles( galaxy_galout_hot[w]+galaxy_galout_cold[w] )
          mbv_mode.(k).galaxyRates.halo_acc.total[m,*]   = $
            percentiles( halo_galacc_hot[w]+halo_galacc_cold[w] )
          mbv_mode.(k).galaxyRates.halo_out.total[m,*]   = $
            percentiles( halo_galout_hot[w]+halo_galout_cold[w] )
          
          ; medians of sums
          mbv_mode.(k).galaxyRates.yvals_gal_out.hot[m,*] = $
            percentiles( halo_galout_hot[w] + galaxy_galout_hot[w] )
          mbv_mode.(k).galaxyRates.yvals_gal_acc.hot[m,*] = $
            percentiles( galaxy_galacc_hot[w] + halo_galout_hot[w] )

          mbv_mode.(k).galaxyRates.yvals_gal_out.cold[m,*] = $
            percentiles( halo_galout_cold[w] + galaxy_galout_cold[w] )
          mbv_mode.(k).galaxyRates.yvals_gal_acc.cold[m,*] = $
            percentiles( galaxy_galacc_cold[w] + halo_galout_cold[w] )
            
          mbv_mode.(k).galaxyRates.yvals_gal_out.total[m,*] = $
            percentiles( halo_galout_hot[w] + galaxy_galout_hot[w] + $
                         halo_galout_cold[w] + galaxy_galout_cold[w] )
          mbv_mode.(k).galaxyRates.yvals_gal_acc.total[m,*] = $
            percentiles( galaxy_galacc_hot[w] + halo_galout_hot[w] + $
                         galaxy_galacc_cold[w] + halo_galout_cold[w] )
                         
          ; (3b) NEW: for halo
          mbv_mode.(k).haloRates.galaxy_acc.hot[m,*] = percentiles( galaxy_haloacc_hot[w] )
          mbv_mode.(k).haloRates.galaxy_out.hot[m,*] = percentiles( galaxy_haloout_hot[w] )
          mbv_mode.(k).haloRates.halo_acc.hot[m,*]   = percentiles( halo_haloacc_hot[w] )
          mbv_mode.(k).haloRates.halo_out.hot[m,*]   = percentiles( halo_haloout_hot[w] )
          
          mbv_mode.(k).haloRates.galaxy_acc.cold[m,*] = percentiles( galaxy_haloacc_cold[w] )
          mbv_mode.(k).haloRates.galaxy_out.cold[m,*] = percentiles( galaxy_haloout_cold[w] )
          mbv_mode.(k).haloRates.halo_acc.cold[m,*]   = percentiles( halo_haloacc_cold[w] )
          mbv_mode.(k).haloRates.halo_out.cold[m,*]   = percentiles( halo_haloout_cold[w] )
          
          mbv_mode.(k).haloRates.galaxy_acc.total[m,*] = $
            percentiles( galaxy_haloacc_hot[w]+galaxy_haloacc_cold[w] )
          mbv_mode.(k).haloRates.galaxy_out.total[m,*] = $
            percentiles( galaxy_haloout_hot[w]+galaxy_haloout_cold[w] )
          mbv_mode.(k).haloRates.halo_acc.total[m,*]   = $
            percentiles( halo_haloacc_hot[w]+halo_haloacc_cold[w] )
          mbv_mode.(k).haloRates.halo_out.total[m,*]   = $
            percentiles( halo_haloout_hot[w]+halo_haloout_cold[w] )
          
          ; medians of sums (net acc rate into halo = galaxy_haloacc + halo_haloacc)
          mbv_mode.(k).haloRates.yvals_halo.hot[m,*] = $
            percentiles( galaxy_haloacc_hot[w] + halo_haloacc_hot[w] )

          mbv_mode.(k).haloRates.yvals_halo.cold[m,*] = $
            percentiles( galaxy_haloacc_cold[w] + halo_haloacc_cold[w] )
            
          mbv_mode.(k).haloRates.yvals_halo.total[m,*] = $
            percentiles( galaxy_haloacc_hot[w]  + halo_haloacc_hot[w] + $
                         galaxy_haloacc_cold[w] + halo_haloacc_cold[w] )
          
          ; (3c) NEW: mode ratios
          if k gt 0 then begin
            mbv_mode.(k).modeRatiosNet.galaxy.hot[m,*]   = percentiles( netratio_gal_hot[w] )
            mbv_mode.(k).modeRatiosNet.galaxy.cold[m,*]  = percentiles( netratio_gal_cold[w] )
            mbv_mode.(k).modeRatiosNet.galaxy.total[m,*] = percentiles( netratio_gal_total[w] )
            
            mbv_mode.(k).modeRatiosNet.halo.hot[m,*]   = percentiles( netratio_halo_hot[w] )
            mbv_mode.(k).modeRatiosNet.halo.cold[m,*]  = percentiles( netratio_halo_cold[w] )
            mbv_mode.(k).modeRatiosNet.halo.total[m,*] = percentiles( netratio_halo_total[w] )
          endif
          
          ; (3d) NEW: fractions
          mbv_mode.(k).coldFractionsNet.galaxy[m,*] = $ ; of net inflow
            percentiles( galaxy_galacc_cold[w] / (galaxy_galacc_hot[w]+galaxy_galacc_cold[w]) )
            
          mbv_mode.(k).coldFractionsNet.halo[m,*] = $
            percentiles( (halo_haloacc_cold[w] + galaxy_haloacc_cold[w]) / $
                         (halo_haloacc_cold[w] + halo_haloacc_cold[w] + $
                          halo_haloacc_hot[w]  + halo_haloacc_hot[w]) )
                          
          mbv_mode.(k).coldFractionsAcc.galaxy[m,*] = $ ; of raw inflow
            percentiles( (galaxy_galacc_cold[w] + halo_galout_cold[w]) / $
                         (galaxy_galacc_cold[w] + halo_galout_cold[w] + $
                          galaxy_galacc_hot[w]  + halo_galout_hot[w]) )
            
          mbv_mode.(k).coldFractionsAcc.halo[m,*] = $
            percentiles( (halo_haloacc_cold[w] + galaxy_haloacc_cold[w]) / $
                         (halo_haloacc_cold[w] + galaxy_haloacc_cold[w] + $
                          halo_haloacc_hot[w]  + galaxy_haloacc_hot[w]) )  
          
          ; (4) ABCD rates
          mbv_mode.(k).galaxyRates.abcd_gal_out.hot[m,*] = $ ; A+B+C = 3+2
            percentiles( galaxy_galout_hot[w] + halo_galout_hot[w] )
          mbv_mode.(k).galaxyRates.abcd_gal_acc.hot[m,*] = $ ; A+C+D = 1+4
            percentiles( galaxy_galacc_hot[w] + halo_galacc_hot[w] )
          mbv_mode.(k).galaxyRates.abcd_gal_net.hot[m,*] = $ ; B-D = 1+4-2-3
            percentiles( galaxy_galacc_hot[w] + halo_galacc_hot[w] - $
                         galaxy_galout_hot[w] - halo_galout_hot[w])

          mbv_mode.(k).galaxyRates.abcd_gal_out.cold[m,*] = $ 
            percentiles( galaxy_galout_cold[w] + halo_galout_cold[w] )
          mbv_mode.(k).galaxyRates.abcd_gal_acc.cold[m,*] = $
            percentiles( galaxy_galacc_cold[w] + halo_galacc_cold[w] )
          mbv_mode.(k).galaxyRates.abcd_gal_net.cold[m,*] = $
            percentiles( galaxy_galacc_cold[w] + halo_galacc_cold[w] - $
                         galaxy_galout_cold[w] - halo_galout_cold[w])
            
          mbv_mode.(k).galaxyRates.abcd_gal_out.total[m,*] = $ 
            percentiles( galaxy_galout_hot[w]  + halo_galout_hot[w] + $
                         galaxy_galout_cold[w] + halo_galout_cold[w] )
          mbv_mode.(k).galaxyRates.abcd_gal_acc.total[m,*] = $
            percentiles( galaxy_galacc_hot[w]  + halo_galacc_hot[w] + $
                         galaxy_galacc_cold[w] + halo_galacc_cold[w] )
          mbv_mode.(k).galaxyRates.abcd_gal_net.total[m,*] = $
            percentiles( (galaxy_galacc_hot[w]  + halo_galacc_hot[w] + $
                          galaxy_galacc_cold[w] + halo_galacc_cold[w]) - $
                         (galaxy_galout_hot[w]  + halo_galout_hot[w] + $
                          galaxy_galout_cold[w] + halo_galout_cold[w]) )
          
          ; (5) OLD: rates (acc/out/net) with percentiles
          mbv_mode.(k).galaxyPerc.acc.hot[m,*] = $
            percentiles( mbv_modeAll.(k).galaxy.accRate.total.hot.tVirAcc[tVirInd,w] )
          mbv_mode.(k).galaxyPerc.out.hot[m,*] = $
            percentiles( mbv_modeAll.(k).galaxy.outRate.total.hot.tVirAcc[tVirInd,w] )
          mbv_mode.(k).galaxyPerc.net.hot[m,*] = $
            percentiles( mbv_modeAll.(k).galaxy.netRate.total.hot.tVirAcc[tVirInd,w] )
          mbv_mode.(k).galaxyPerc.acc.cold[m,*] = $
            percentiles( mbv_modeAll.(k).galaxy.accRate.total.cold.tVirAcc[tVirInd,w] )
          mbv_mode.(k).galaxyPerc.out.cold[m,*] = $
            percentiles( mbv_modeAll.(k).galaxy.outRate.total.cold.tVirAcc[tVirInd,w] )
          mbv_mode.(k).galaxyPerc.net.cold[m,*] = $
            percentiles( mbv_modeAll.(k).galaxy.netRate.total.cold.tVirAcc[tVirInd,w] ) 
            
          mbv_mode.(k).haloPerc.acc.hot[m,*] = $
            percentiles( mbv_modeAll.(k).halo.accRate.total.hot.tVirAcc[tVirInd,w] )
          mbv_mode.(k).haloPerc.out.hot[m,*] = $
            percentiles( mbv_modeAll.(k).halo.outRate.total.hot.tVirAcc[tVirInd,w] )
          mbv_mode.(k).haloPerc.net.hot[m,*] = $
            percentiles( mbv_modeAll.(k).halo.netRate.total.hot.tVirAcc[tVirInd,w] )
          mbv_mode.(k).haloPerc.acc.cold[m,*] = $
            percentiles( mbv_modeAll.(k).halo.accRate.total.cold.tVirAcc[tVirInd,w] )
          mbv_mode.(k).haloPerc.out.cold[m,*] = $
            percentiles( mbv_modeAll.(k).halo.outRate.total.cold.tVirAcc[tVirInd,w] )
          mbv_mode.(k).haloPerc.net.cold[m,*] = $
            percentiles( mbv_modeAll.(k).halo.netRate.total.cold.tVirAcc[tVirInd,w] ) 
            
          ; wind_eta comparison (mode='all' only)
          if k eq 0 then begin
            ; definitions
            eta_out = galaxy_galout_hot + galaxy_galout_cold + halo_galout_hot + halo_galout_cold ;#2+#3
            eta_net = galaxy_galacc_hot + galaxy_galacc_cold ;#1
            
            mbv_mode.(k).eta_w_binned[m,*] = percentiles( eta_out[w] / eta_net[w] )

            ; apply correction (using linear fit)
            sfr_acc_correction_factor = 10.0^(mbv_loc.logMassBinCen[m] * fit[1] + fit[0])
            mbv_mode.(k).eta_w_binned[m,*] /= sfr_acc_correction_factor
            
            ; ABCD approach
            eta_out = galaxy_galout_hot + galaxy_galout_cold + halo_galout_hot + halo_galout_cold ;#2+#3
            eta_acc = galaxy_galacc_hot + galaxy_galacc_cold + halo_galacc_hot + halo_galacc_cold ;#1+#4
            eta_net = galaxy_galacc_hot + galaxy_galacc_cold + halo_galacc_hot + halo_galacc_cold $
                    - (halo_galout_hot + halo_galout_cold + galaxy_galout_hot + galaxy_galout_cold) 
                    ;#1+#4-#2-#3
                   
            mbv_mode.(k).eta_abcd_binned[m,*]  = percentiles( eta_out[w] / eta_net[w] )
            mbv_mode.(k).eta_abcd2_binned[m,*] = percentiles( eta_out[w] / eta_acc[w] )
            
          endif
          
        endfor ; m (mass bins)
        
      endfor ; k (modes)
             
      sP_z  = mod_struct(sP_z, 'redshift'+str(j)+'_'+str(round(redshift*100)), sP_mode)
      mbv_z = mod_struct(mbv_z, 'redshift'+str(j)+'_'+str(round(redshift*100)), mbv_mode)
    endforeach ; redshifts,j
    
    ; put this mode collection into mbv, once per run, and sP collection also
    mbv = mod_struct(mbv, 'mbv'+str(i)+'_'+run, mbv_z)
    sP  = mod_struct(sP, 'sP'+str(i)+'_'+run, sP_z)
  endforeach ;runs,i

  save,mbv,sP,filename=tempFilename
  print,'Saved TEMP file: ['+tempFilename+']'
  
end

; plotsVsRedshift(): plot fractional rate contributions by mode, as a function of redshift
;  and also at 4 discrete redshifts, as a function of halo mass

pro plotsVsRedshift
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  
  ; config
  runs       = ['feedback','tracer']
  redshifts  = [5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.75, 0.5, 0.25, 0.0]
  res        = 512
  timeWindow = 500.0    ; consider accretion over this past time range (Myr)
                        ; 250.0 500.0 1000.0 "all" "tVir_tIGM" or "tVir_tIGM_bin"
  tVirInd    = 0        ; use Tmax/Tviracc=1 to separate hot vs. cold
  accModes   = ['all','smooth','clumpy','stripped','recycled']
  accRateModel = 0      ; 0=prim+rec, 2=prim
  
  ; plot config
  massInds = [21]        ; index logMassBinCen[], for plot vs redshift
  zInds    = [4,6,8,12]  ; index redshifts[], for 2x2 panels vs halo mass
  
  plotUnused   = 0     ; 0=only make plots in paper, 1=make all
  
  xrange_z     = [-0.15,5.15] ; redshift
  xrange_halo  = [9.0,12.0]   ; log halo mass
  yrange_frac  = [0.01,2.0]   ; log fraction
  yrange_frac2 = [0.0,1.0]    ; non-log fraction
  yrange_spec  = [-4.0,-1.0]  ; log specific net rate
  yrange_rate  = [1e-3,2e1]   ; net rate
  
  lines   = [4,0,1,2,3]  ; for each accmode ratio to all
  linesHC = [0,1,2]      ; hot,cold or gal,halo (3rd for total/combined)
  linesIO = [1,2,0]      ; in,out,net
  sK      = 3            ; smoothing kernel size
  psyms   = [-16,-9]      ; symbols for runs (closed circles=FB, none=noFB)
  cts     = ['brewerC-blues','brewerC-greens']
  cInd    = 1     ; color index into sP.colors
  symsize = 0.9   ; symbol size
  
  py_min  = 0.25  ; gray band min in fraction
  py_max  = 1.0   ; gray band max in fraction
  py_min2 = 0.25  ; 2nd gray band min in fraction
  py_max2 = 0.5   ; 2nd gray band max in fraction
  px_cl   = 'light gray'  ; mass bin
  py_cl   = [240,240,240] ; fraction, lighter
  py_cl2  = [230,230,230] ; fraction, medium
  px_cl3  = [200,200,200] ; darker
  
  ; load
  if accModes[0] ne 'all' then message,'Not going to work.'
  if runs[0] ne 'feedback' then message,'Maybe a problem.'
  
  tempFilename = '/n/home07/dnelson/plots/temp_plotsVsRedshift_'+str(n_elements(redshifts))+$
                 '_tw'+str(round(timeWindow))+'_arm'+str(accRateModel)+'.sav'
                 
  if ~file_test(tempFilename) then $
    plotsVsRedshiftBin, tempFilename, runs, redshifts, res, timeWindow, tVirInd, accModes, accRateModel
    
  restore,tempFilename,/verbose

  wind_model = wind_model(simPath=sP.(0).(0).(0).simPath, $
                          redshifts=redshifts, $
                          sMass=mbv.(0).(0).(0).logMassBinCen[massInds])
  
  ; strings
  if isnumeric(timeWindow) then twStr = str(fix(timeWindow))
  if ~isnumeric(timeWindow) then twStr = "-"+timeWindow
  
  plotStr   = ''
  simNames  = []
  simColors = []
  simNames2  = []
  simColors2 = []
  
  for i=0,n_elements(runs)-1 do begin
    if sP.(i).(0).(0).run eq 'gadget' then continue
    
    plotStr   = plotStr + sP.(i).(0).(0).plotPrefix + '.'
    simNames  = [simNames, sP.(i).(0).(0).simName]
    simColors = [simColors, sP.(i).(0).(0).colors[cInd]]
    simNames2  = [simNames2, sP.(i).(0).(0).simName, sP.(i).(0).(0).simName]
    simColors2 = [simColors2, sampleColorTable(cts[i],2,bounds=[0.4,0.9])]
  endfor

  if accRateModel eq 0 then dMdt_top = "p+r"
  if accRateModel eq 2 then dMdt_top = "prim"
  
  plotStr  += str(res) + '_' + str(sP.(0).(0).(0).snap) + '_tw' + twStr + $
              '_arm' + str(accRateModel) 
  plotStrMB = plotStr + '_mB-' + strjoin(str(massInds),",")
  print,'Plotting mass bin centered on (log msun/h): '+str(mbv.(0).(0).(0).logMassBinCen[massInds])
  print,'Plotting redshifts: '+strjoin(redshifts[zInds]," ")
             
  minMass = mean(mbv.(0).(0).(0).logMassBinCen[massInds-1:massInds])
  maxMass = mean(mbv.(0).(0).(0).logMassBinCen[massInds:massInds+1])
  px_min = codeMassToLogMsun( logMsunToCodeMass( minMass ) / units.hubbleParam )
  px_max = codeMassToLogMsun( logMsunToCodeMass( maxMass ) / units.hubbleParam )
  print,'Plotting massBin (log msun no h units): ',px_min,px_max
             
  ; plots vs redshift (all accModes together, @ one halo mass bin)
  ; ---------------------------------------
  pos_single = [0.17,0.13,0.94,0.87] ; single plot with second x-axis above

  ; plot (1) - fractional rates, gal hot/cold (all modes at once)
  if plotUnused then begin
  start_PS, sP.(0).(0).(0).plotPath + 'rateFracsVsRedshift.galaxyOld.' + plotStrMB + '-allModes.eps'
  
    yrange2 = yrange_frac * [1.0,0.7]

    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=5,ys=5,/ylog,/noerase,pos=pos_single
      
    cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
      [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
    cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
      [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
        
    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=9,/ys,/ylog,yminor=0,$
      ytitle=textoidl("(dM/dt)_{mode} / (dM/dt)_{total}"),xtitle="Redshift",/noerase,pos=pos_single
    universeage_axis, xrange_z, yrange2, /ylog
        
    ; loop over runs
    for i=0,n_tags(sP)-1 do begin
      ; loop over mass bins (just one now)
      massStrs = []
      foreach massInd,massInds,k do begin
        
        massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
          
        ; loop over each accMode
        for j=1,n_tags(sP.(i).(0))-1 do begin
        
          xvals           = []
          yvals_gal_hot   = []
          yvals_gal_cold  = []
          yvals_halo_hot  = []
          yvals_halo_cold = []
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals = [xvals, sP.(i).(m).(j).redshift]
            yvals_gal_hot   = [yvals_gal_hot,   mbv.(i).(m).(j).mode_ratio.gal_hot[massInd]]
            yvals_gal_cold  = [yvals_gal_cold,  mbv.(i).(m).(j).mode_ratio.gal_cold[massInd]]
            yvals_halo_hot  = [yvals_halo_hot,  mbv.(i).(m).(j).mode_ratio.halo_hot[massInd]]
            yvals_halo_cold = [yvals_halo_cold, mbv.(i).(m).(j).mode_ratio.halo_cold[massInd]]
          endfor ;m
        
          ; smooth and set zeros to a plottable value
          yvals_gal_hot   = smooth(yvals_gal_hot,sK,/nan) < 1.0
          yvals_gal_cold  = smooth(yvals_gal_cold,sK,/nan) < 1.0
          yvals_halo_hot  = smooth(yvals_halo_hot,sK,/nan) < 1.0
          yvals_halo_cold = smooth(yvals_halo_cold,sK,/nan) < 1.0
      
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal_hot,color=simColors2[i*2+0],$
            line=lines[j],psym=-9,symsize=symsize,/overplot ;psym=psyms[i],
          cgPlot,xvals,yvals_gal_cold,color=simColors2[i*2+1],$
            line=lines[j],psym=-16,symsize=symsize,/overplot ;psym=psyms[i],
          
        endfor ; accModes 
      endforeach ; massBins
    endfor ; runs
      
    ;legend,textoidl("M_{halo} = "+string(massStrs,format='(f4.1)')),/top,/right
    legend,accModes[1:*],linestyle=lines[1:*],pos=[1.1,yrange2[0]*3.4]
      
    tvirStrs = textoidl(["T_{max} > T_{vir,acc}","T_{max} < T_{vir,acc}",$
                         "T_{max} > T_{vir,acc}","T_{max} < T_{vir,acc}"])
    legend,tvirStrs+" ("+simNames2+")",textcolor=simColors2,/bottom,/right,spacing=1.5
  
  end_PS
  endif ;plotUnused
  
  ; plot (1b) - fractional rates, gal hot/cold (all modes at once, updated)
  start_PS, sP.(0).(0).(0).plotPath + 'rateFracsVsRedshift.galaxy.' + plotStrMB + '-allModes.eps'
  
    yrange2 = yrange_frac * [1.0,0.7]

    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=5,ys=5,/ylog,/noerase,pos=pos_single
      
    cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
      [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
    cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
      [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
        
    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=9,/ys,/ylog,yminor=0,$
      ytitle=textoidl("(dM/dt)_{mode} / (dM/dt)_{total}"),xtitle="Redshift",/noerase,pos=pos_single
    universeage_axis, xrange_z, yrange2, /ylog
        
    ; loop over runs
    for i=0,n_tags(sP)-1 do begin
      ; loop over mass bins (just one now)
      massStrs = []
      foreach massInd,massInds,k do begin
        
        massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
          
        ; loop over each accMode
        for j=1,n_tags(sP.(i).(0))-1 do begin
        
          xvals            = fltarr(n_elements(redshifts))
          yvals_gal_hot    = fltarr(n_elements(redshifts),3)
          yvals_gal_cold   = fltarr(n_elements(redshifts),3)
          yvals_gal_total  = fltarr(n_elements(redshifts),3)
          yvals_halo_hot   = fltarr(n_elements(redshifts),3)
          yvals_halo_cold  = fltarr(n_elements(redshifts),3)
          yvals_halo_total = fltarr(n_elements(redshifts),3)
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals[m] = sP.(i).(m).(j).redshift
            
            yvals_gal_hot[m,*]    = mbv.(i).(m).(j).modeRatiosNet.galaxy.hot[massInd,*]
            yvals_gal_cold[m,*]   = mbv.(i).(m).(j).modeRatiosNet.galaxy.cold[massInd,*]
            yvals_gal_total[m,*]  = mbv.(i).(m).(j).modeRatiosNet.galaxy.total[massInd,*]
           
            yvals_halo_hot[m,*]   = mbv.(i).(m).(j).modeRatiosNet.halo.hot[massInd,*]
            yvals_halo_cold[m,*]  = mbv.(i).(m).(j).modeRatiosNet.halo.cold[massInd,*]
            yvals_halo_total[m,*] = mbv.(i).(m).(j).modeRatiosNet.halo.total[massInd,*]
          endfor ;m
        
          ; smooth and set zeros to a plottable value
          for m=0,2 do begin
            yvals_gal_hot[*,m]   = smooth(yvals_gal_hot[*,m],sK)
            yvals_gal_cold[*,m]  = smooth(yvals_gal_cold[*,m],sK)
            yvals_gal_total[*,m] = smooth(yvals_gal_total[*,m],sK)
            
            yvals_halo_hot[*,m]   = smooth(yvals_halo_hot[*,m],sK)
            yvals_halo_cold[*,m]  = smooth(yvals_halo_cold[*,m],sK)
            yvals_halo_total[*,m] = smooth(yvals_halo_total[*,m],sK)
          endfor
        
          ; smooth and set zeros to a plottable value
          yvals_gal_hot   = smooth(yvals_gal_hot,sK,/nan) < 1.0
          yvals_gal_cold  = smooth(yvals_gal_cold,sK,/nan) < 1.0
          yvals_halo_hot  = smooth(yvals_halo_hot,sK,/nan) < 1.0
          yvals_halo_cold = smooth(yvals_halo_cold,sK,/nan) < 1.0
      
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal_hot[*,1],color=simColors2[i*2+0],$
            line=lines[j],psym=-9,symsize=symsize,/overplot
          cgPlot,xvals,yvals_gal_cold[*,1],color=simColors2[i*2+1],$
            line=lines[j],psym=-16,symsize=symsize,/overplot
           
          ;cgPlot,xvals,yvals_gal_total[*,1],color=simColors2[i*2+1],$
          ;  line=lines[j],psym=-16,symsize=symsize,/overplot
          
        endfor ; accModes 
      endforeach ; massBins
    endfor ; runs
      
    ;legend,textoidl("M_{halo} = "+string(massStrs,format='(f4.1)')),/top,/right
    legend,accModes[1:*],linestyle=lines[1:*],pos=[1.1,yrange2[0]*3.4]
      
    tvirStrs = textoidl(["T_{max} > T_{vir,acc}","T_{max} < T_{vir,acc}",$
                         "T_{max} > T_{vir,acc}","T_{max} < T_{vir,acc}"])
    legend,tvirStrs+" ("+simNames2+")",textcolor=simColors2,/bottom,/right,spacing=1.5
  
  end_PS
    
  ; plot (2) - fractional rates, gal hot/cold (2x2 panels, one per mode)
  if plotUnused then begin
  start_PS, sP.(0).(0).(0).plotPath + 'rateFracsVsRedshiftPanels.galaxyOld.' + plotStrMB + '-allModes.eps', /extrabig
  
    pos = plot_pos(rows=2,cols=2,/gap)
    
    ; loop over each accMode
    for j=1,n_tags(sP.(0).(0))-1 do begin
    
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_frac,xs=5,ys=5,/ylog,/noerase,$
        pos=pos[j-1]-[0,0.02,0,0.04]
        
      cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
        [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
      cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
        [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
        
      barStr = str(round(py_max2*100))+'% - '+str(round(py_max*100))+'%'
      if j eq 3 then $
        cgText, xrange_z[1]-1.0, py_min+0.39, barStr, alignment=0.5, color=cgColor('gray')
        
      barStr = str(round(py_min*100))+'% - '+str(round(py_max2*100))+'% '
      if j eq 3 then $
        cgText, xrange_z[1]-1.0, py_min+0.06, barStr, alignment=0.5, color=cgColor('gray')
        
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_frac,xs=9,/ys,/ylog,yminor=0,$
        ytitle=textoidl("(dM/dt)^{"+dMdt_top+"}_{"+accModes[j]+"} / (dM/dt)^{"+dMdt_top+"}_{total}"),$
        xtitle="Redshift",/noerase,$
        pos=pos[j-1]-[0,0.02,0,0.04]
      universeage_axis, xrange_z, yrange_frac, /ylog, spaceFac=1.5
    
      ; loop over runs
      for i=0,n_tags(sP)-1 do begin
        ; accMode does not exist for this run? then skip
        if j ge n_tags(sP.(i).(0)) then continue
        
        ; loop over mass bins (just one now)
        massStrs = []
        foreach massInd,massInds,k do begin
        
          massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
          
          xvals           = []
          yvals_gal_hot   = []
          yvals_gal_cold  = []
          yvals_halo_hot  = []
          yvals_halo_cold = []
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals = [xvals, sP.(i).(m).(j).redshift]
            yvals_gal_hot   = [yvals_gal_hot,   mbv.(i).(m).(j).mode_ratio.gal_hot[massInd]]
            yvals_gal_cold  = [yvals_gal_cold,  mbv.(i).(m).(j).mode_ratio.gal_cold[massInd]]
            yvals_halo_hot  = [yvals_halo_hot,  mbv.(i).(m).(j).mode_ratio.halo_hot[massInd]]
            yvals_halo_cold = [yvals_halo_cold, mbv.(i).(m).(j).mode_ratio.halo_cold[massInd]]
          endfor ;m
        
          ; smooth and set zeros to a plottable value
          yvals_gal_hot   = smooth(yvals_gal_hot,sK) < 1.0
          yvals_gal_cold  = smooth(yvals_gal_cold,sK) < 1.0
          yvals_halo_hot  = smooth(yvals_halo_hot,sK) < 1.0
          yvals_halo_cold = smooth(yvals_halo_cold,sK) < 1.0
      
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal_hot,color=simColors[i],psym=psyms[i],symsize=symsize,line=linesHC[0],/overplot
          cgPlot,xvals,yvals_gal_cold,color=simColors[i],psym=psyms[i],symsize=symsize,line=linesHC[1],/overplot
          
        endforeach ; massBins
      endfor ; runs
    
      if j eq 4 then $
      legend,textoidl(["T_{max} > T_{vir,acc}","T_{max} < T_{vir,acc}"]),$
        linestyle=linesHC[0:1],color=cgColor('black'),textcolor=cgColor('black'),/bottom,/right
      if j eq 2 then $
      legend,simNames,textcolor=simColors,/bottom,/left
    
    endfor ; accModes
  
  end_PS
  endif ;plotUnused
  
  ; plot (2b) - fractional rates, gal hot/cold (2x2 panels, one per mode, updated)
  start_PS, sP.(0).(0).(0).plotPath + 'rateFracsVsRedshiftPanels.galaxy.' + plotStrMB + '-allModes.eps', /extrabig
  
    pos = plot_pos(rows=2,cols=2,/gap)
    
    ; loop over each accMode
    for j=1,n_tags(sP.(0).(0))-1 do begin
    
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_frac,xs=5,ys=5,/ylog,/noerase,$
        pos=pos[j-1]-[0,0.02,0,0.04]
        
      cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
        [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
      cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
        [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
        
      barStr = str(round(py_max2*100))+'% - '+str(round(py_max*100))+'%'
      if j eq 3 then $
        cgText, xrange_z[1]-1.0, py_min+0.39, barStr, alignment=0.5, color=cgColor('gray')
        
      barStr = str(round(py_min*100))+'% - '+str(round(py_max2*100))+'% '
      if j eq 3 then $
        cgText, xrange_z[1]-1.0, py_min+0.06, barStr, alignment=0.5, color=cgColor('gray')
        
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_frac,xs=9,/ys,/ylog,yminor=0,$
        ytitle=textoidl("(dM/dt)^{"+dMdt_top+"}_{"+accModes[j]+"} / (dM/dt)^{"+dMdt_top+"}_{total}"),$
        xtitle="Redshift",/noerase,$
        pos=pos[j-1]-[0,0.02,0,0.04]
      universeage_axis, xrange_z, yrange_frac, /ylog, spaceFac=1.5
    
      ; loop over runs
      for i=0,n_tags(sP)-1 do begin
        ; accMode does not exist for this run? then skip
        if j ge n_tags(sP.(i).(0)) then continue
        
        ; loop over mass bins (just one now)
        massStrs = []
        foreach massInd,massInds,k do begin
        
          massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
          
          xvals            = fltarr(n_elements(redshifts))
          yvals_gal_hot    = fltarr(n_elements(redshifts),3)
          yvals_gal_cold   = fltarr(n_elements(redshifts),3)
          yvals_gal_total  = fltarr(n_elements(redshifts),3)
          yvals_halo_hot   = fltarr(n_elements(redshifts),3)
          yvals_halo_cold  = fltarr(n_elements(redshifts),3)
          yvals_halo_total = fltarr(n_elements(redshifts),3)
          
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals[m] = sP.(i).(m).(j).redshift
            
            yvals_gal_hot[m,*]    = mbv.(i).(m).(j).modeRatiosNet.galaxy.hot[massInd,*]
            yvals_gal_cold[m,*]   = mbv.(i).(m).(j).modeRatiosNet.galaxy.cold[massInd,*]
            yvals_gal_total[m,*]  = mbv.(i).(m).(j).modeRatiosNet.galaxy.total[massInd,*]
           
            yvals_halo_hot[m,*]   = mbv.(i).(m).(j).modeRatiosNet.halo.hot[massInd,*]
            yvals_halo_cold[m,*]  = mbv.(i).(m).(j).modeRatiosNet.halo.cold[massInd,*]
            yvals_halo_total[m,*] = mbv.(i).(m).(j).modeRatiosNet.halo.total[massInd,*]
          endfor ;m
        
          ; smooth and set zeros to a plottable value
          for m=0,2 do begin
            yvals_gal_hot[*,m]   = smooth(yvals_gal_hot[*,m],sK)
            yvals_gal_cold[*,m]  = smooth(yvals_gal_cold[*,m],sK)
            yvals_gal_total[*,m] = smooth(yvals_gal_total[*,m],sK)
            
            yvals_halo_hot[*,m]   = smooth(yvals_halo_hot[*,m],sK)
            yvals_halo_cold[*,m]  = smooth(yvals_halo_cold[*,m],sK)
            yvals_halo_total[*,m] = smooth(yvals_halo_total[*,m],sK)
          endfor
      
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal_hot[*,1],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesHC[0],/overplot
          cgPlot,xvals,yvals_gal_cold[*,1],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesHC[1],/overplot
          
          ;cgPlot,xvals,yvals_gal_total[*,1],$
          ;  color=simColors[i],psym=psyms[i],symsize=symsize,line=linesHC[2],/overplot
          
        endforeach ; massBins
      endfor ; runs
    
      if j eq 4 then $
      legend,textoidl(["T_{max} > T_{vir,acc}","T_{max} < T_{vir,acc}"]),$
        linestyle=linesHC[0:1],color=cgColor('black'),textcolor=cgColor('black'),/bottom,/right
      if j eq 2 then $
      legend,simNames,textcolor=simColors,/bottom,/left
    
    endfor ; accModes
  
  end_PS  
  
  ; plot (3) - coldfrac_net vs redshift, galaxy only, all accModes on plot
  if plotUnused then begin
  start_PS, sP.(0).(0).(0).plotPath + 'coldFracNetVsRedshift.galaxyOld.' + plotStrMB + '-allModes.eps', /small
  
    yrange2 = yrange_frac * [1.0,0.7]

    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=5,ys=5,/ylog,/noerase,pos=pos_single
      
    cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
      [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
    cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
      [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
        
    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=9,/ys,/ylog,yminor=0,$
      ytitle=textoidl("(dM/dt)^{"+dMdt_top+"}_{T_{max}<T_{vir,acc}} / (dM/dt)^{"+dMdt_top+$
             "}_{total,mode }"),$
      xtitle="Redshift",/noerase,pos=pos_single
    universeage_axis, xrange_z, yrange2, /ylog
        
    ; loop over runs
    for i=0,n_tags(sP)-1 do begin
      ; loop over mass bins (just one now)
      massStrs = []
      foreach massInd,massInds,k do begin
        
        massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
          
        ; loop over each accMode (skip all)
        for j=1,n_tags(sP.(i).(0))-1 do begin
        
          xvals      = []
          yvals_gal  = []
          yvals_halo = []
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals = [xvals, sP.(i).(m).(j).redshift]
            
            yvals_gal  = [yvals_gal,  mbv.(i).(m).(j).coldfrac_net.gal[massInd]]
            yvals_halo = [yvals_halo, mbv.(i).(m).(j).coldfrac_net.halo[massInd]]
          endfor ;m
        
          ; smooth and set zeros to a plottable value
          yvals_gal  = smooth(yvals_gal,sK,/nan) < 1.0
          yvals_halo = smooth(yvals_halo,sK,/nan) < 1.0
      
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal,color=simColors[i],psym=psyms[i],line=lines[j],symsize=symsize,/overplot
          
        endfor ; accModes 
      endforeach ; massBins
    endfor ; runs
      
    ;legend,textoidl("M_{halo} = "+string(massStrs,format='(f4.1)')),/top,/right
    legend,accModes[1:*],linestyle=lines[1:*],linesize=0.6,/bottom,/right
    legend,simNames,textcolor=simColors,/top,/left
  
  end_PS
  endif ;plotUnused
  
  ; plot (3b) - coldfrac_net vs redshift, galaxy only, all accModes on plot (updated)
  start_PS, sP.(0).(0).(0).plotPath + 'coldFracNetVsRedshift.galaxy.' + plotStrMB + '-allModes.eps', /small
  
    yrange2 = yrange_frac * [1.0,0.7]

    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=5,ys=5,/ylog,/noerase,pos=pos_single
      
    cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
      [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
    cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
      [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
        
    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=9,/ys,/ylog,yminor=0,$
      ytitle=textoidl("(dM/dt)^{"+dMdt_top+"}_{T_{max}<T_{vir,acc}} / (dM/dt)^{"+dMdt_top+$
             "}_{total,mode }"),$
      xtitle="Redshift",/noerase,pos=pos_single
    universeage_axis, xrange_z, yrange2, /ylog
        
    ; loop over runs
    for i=0,n_tags(sP)-1 do begin
      ; loop over mass bins (just one now)
      massStrs = []
      foreach massInd,massInds,k do begin
        
        massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
          
        ; loop over each accMode (skip all)
        for j=1,n_tags(sP.(i).(0))-1 do begin
        
          xvals      = fltarr(n_elements(redshifts))
          yvals_gal  = fltarr(n_elements(redshifts),3)
          yvals_halo = fltarr(n_elements(redshifts),3)
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals[m] = sP.(i).(m).(j).redshift
            
            yvals_gal[m,*]  = mbv.(i).(m).(j).coldFractionsNet.galaxy[massInd,*]
            yvals_halo[m,*] = mbv.(i).(m).(j).coldFractionsNet.halo[massInd,*]
          endfor ;m
        
          ; smooth and set zeros to a plottable value
          for m=0,2 do begin
            yvals_gal[*,m]  = smooth(yvals_gal[*,m],sK,/nan)
            yvals_halo[*,m] = smooth(yvals_halo[*,m],sK,/nan)
          endfor
      
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal[*,1],color=simColors[i],psym=psyms[i],line=lines[j],symsize=symsize,/overplot
          
          ;cgPlot,xvals,yvals_halo[*,1],color=simColors[i],psym=psyms[i],line=lines[j],symsize=symsize,/overplot
          
        endfor ; accModes 
      endforeach ; massBins
    endfor ; runs
      
    legend,accModes[1:*],linestyle=lines[1:*],linesize=0.6,/bottom,/right
    legend,simNames,textcolor=simColors,/top,/left
  
  end_PS
  
  ; plot (3c) - coldfrac_acc vs redshift, galaxy only, all accModes on plot (updated)
  if plotUnused then begin
  start_PS, sP.(0).(0).(0).plotPath + 'coldFracAccVsRedshift.galaxy.' + plotStrMB + '-allModes.eps', /small
  
    yrange2 = yrange_frac * [1.0,0.7]

    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=5,ys=5,/ylog,/noerase,pos=pos_single
      
    cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
      [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
    cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
      [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
        
    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=9,/ys,/ylog,yminor=0,$
      ytitle=textoidl("(dM/dt)^{"+dMdt_top+"}_{T_{max}<T_{vir,acc}} / (dM/dt)^{"+dMdt_top+$
             "}_{total,mode }"),$
      xtitle="Redshift",/noerase,pos=pos_single
    universeage_axis, xrange_z, yrange2, /ylog
        
    ; loop over runs
    for i=0,n_tags(sP)-1 do begin
      ; loop over mass bins (just one now)
      massStrs = []
      foreach massInd,massInds,k do begin
        
        massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
          
        ; loop over each accMode (skip all)
        for j=1,n_tags(sP.(i).(0))-1 do begin
        
          xvals      = fltarr(n_elements(redshifts))
          yvals_gal  = fltarr(n_elements(redshifts),3)
          yvals_halo = fltarr(n_elements(redshifts),3)
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals[m] = sP.(i).(m).(j).redshift
            
            yvals_gal[m,*]  = mbv.(i).(m).(j).coldFractionsAcc.galaxy[massInd,*]
            yvals_halo[m,*] = mbv.(i).(m).(j).coldFractionsAcc.halo[massInd,*]
          endfor ;m
        
          ; smooth and set zeros to a plottable value
          for m=0,2 do begin
            yvals_gal[*,m]  = smooth(yvals_gal[*,m],sK,/nan)
            yvals_halo[*,m] = smooth(yvals_halo[*,m],sK,/nan)
          endfor
      
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal[*,1],color=simColors[i],psym=psyms[i],line=lines[j],symsize=symsize,/overplot
          
          ;cgPlot,xvals,yvals_halo[*,1],color=simColors[i],psym=psyms[i],line=lines[j],symsize=symsize,/overplot
          
        endfor ; accModes 
      endforeach ; massBins
    endfor ; runs
      
    legend,accModes[1:*],linestyle=lines[1:*],linesize=0.6,/bottom,/right
    legend,simNames,textcolor=simColors,/top,/left
  
  end_PS
  endif ;plotUnused

  ; plots vs redshift (one plot per accMode, @ one halo mass bin)
  ; --------------------------------------------------------
  temp = {}
  
  for j=0,n_tags(sP.(1).(0))-1 do begin
    plotStrAM = plotStrMB + '_am-'+accModes[j] 
    
    if plotUnused then begin
    ; plot (4) - net rates vs redshift, galaxy hot/cold
    start_PS, sP.(0).(0).(0).plotPath + 'netRateVsRedshift.galaxyOld.' + plotStrAM + '.eps', /small
    
      yrange2 = yrange_rate
      if accModes[j] eq 'all' then yrange2 *= 4
      
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=9,/ys,/ylog,yminor=0,$
        ytitle=textoidl("dM^{"+dMdt_top+"}_{gas,"+accModes[j]+" }/dt [_{ }M_{sun }/_{ }yr_{ }]"),$
        xtitle="Redshift",/noerase,pos=pos_single
        
      universeage_axis, xrange_z, yrange2, /ylog
      
      ; loop over runs
      temp_sub = {}
      
      for i=0,n_tags(sP)-1 do begin
        ; loop over mass bins (just one now)
        massStrs = []
        foreach massInd,massInds,k do begin
        
          massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
        
          xvals           = []
          yvals_gal_hot   = []
          yvals_gal_cold  = []
          yvals_halo_hot  = []
          yvals_halo_cold = []
          yvals_gal_tot   = []
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals = [xvals, sP.(i).(m).(j).redshift]
            yvals_gal_hot   = [yvals_gal_hot, $
              mbv.(i).(m).(j).galaxyMedian.netRate.total.hot.tVirAcc[tVirInd,massInd]]
            yvals_gal_cold  = [yvals_gal_cold, $
              mbv.(i).(m).(j).galaxyMedian.netRate.total.cold.tVirAcc[tVirInd,massInd]]
            yvals_halo_hot  = [yvals_halo_hot, $
              mbv.(i).(m).(j).haloMedian.netRate.total.hot.tVirAcc[tVirInd,massInd]]
            yvals_halo_cold = [yvals_halo_cold, $
              mbv.(i).(m).(j).haloMedian.netRate.total.cold.tVirAcc[tVirInd,massInd]]
            yvals_gal_tot = [yvals_gal_tot, yvals_gal_hot[-1]+yvals_gal_cold[-1]]
          endfor ;m
                          
          ; smooth and set zeros to a plottable value
          yvals_gal_hot   = smooth(yvals_gal_hot/units.hubbleParam ,sK,/nan)
          yvals_gal_cold  = smooth(yvals_gal_cold/units.hubbleParam ,sK,/nan)
          yvals_halo_hot  = smooth(yvals_halo_hot/units.hubbleParam ,sK,/nan)
          yvals_halo_cold = smooth(yvals_halo_cold/units.hubbleParam ,sK,/nan)
          yvals_gal_tot   = smooth(yvals_gal_tot/units.hubbleParam ,sK,/nan) ;> yrange2[0]/10
        
          ; keep values for this run for subplot
          save_struct = { gal_hot:yvals_gal_hot, gal_cold:yvals_gal_cold, gal_tot:yvals_gal_tot,z:xvals }
          temp_sub = mod_struct( temp_sub, sP.(i).(0).(0).run, save_struct )
        
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal_hot,color=simColors[i],psym=psyms[i],symsize=symsize,line=linesHC[0],/overplot
          cgPlot,xvals,yvals_gal_cold,color=simColors[i],psym=psyms[i],symsize=symsize,line=linesHC[1],/overplot
          ;cgPlot,xvals,yvals_gal_tot,color=simColors[i],psym=psyms[i],symsize=symsize,line=2,/overplot
          
          ;cgPlot,xvals,yvals_halo_hot,color=simColors[i],psym=psyms[i],line=linesHC[0],/overplot
          ;cgPlot,xvals,yvals_halo_cold,color=simColors[i],psym=psyms[i],line=linesHC[1],/overplot
        endforeach ; massBins
        
      endfor ;i
      
      ;legend,textoidl("M_{halo} = "+string(massStrs,format='(f4.1)')),/top,/right
      
      legend,textoidl(["T_{max} > T_{vir,acc}","T_{max} < T_{vir,acc}"]),$
        linestyle=linesHC[0:1],color=cgColor('black'),textcolor=cgColor('black'),/top,/left
        
      if accRateModel eq 2 then $
        legend,simNames,textcolor=simColors,pos=[xrange_z[1]-0.18,0.4],/right
      if accRateModel eq 0 then $
        legend,simNames,textcolor=simColors,/top,/right
        
      ; subplot
      pos_sub = [0.60,0.23,0.89,0.47]
      yrange_sub = [0.0,2.5]
      if accRateModel eq 0 then yrange_sub = [-1.0,1.0]
      
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_sub,/xs,/ys,$
        ytitle=textoidl("log_{ }(_{ }noFB_{ }/_{ }FB_{ })"),$
        xtitle="Redshift",/noerase,pos=pos_sub,charsize=!p.charsize-0.3

      if n_tags(temp_sub) ne 2 then print,'WARNING, check'
      
      if n_tags(temp_sub) eq 2 then begin
        if accRateModel eq 0 then $
          cgPlot,xrange_z,[0.0,0.0],line=0,/overplot,color=cgColor('gray')
        
        yy = smooth( (temp_sub.tracer.gal_hot>0.0) / $
                     (temp_sub.feedback.gal_hot>0.0) ,sK,/nan)

        cgPlot,xvals,alog10(yy),line=linesHC[0],/overplot,color=cgColor('black')
              
        yy = smooth( (temp_sub.tracer.gal_cold>0.0) / $
                     (temp_sub.feedback.gal_cold>0.0) ,sK,/nan)

        cgPlot,xvals,alog10(yy),line=linesHC[1],/overplot,color=cgColor('black')
      endif
      
    end_PS
    endif ;plotUnused
    
    ; plot (4b) - net rates vs redshift, galaxy hot/cold (alternative definitions)
    start_PS, sP.(0).(0).(0).plotPath + 'netRateVsRedshift.galaxy.' + plotStrAM + '.eps', /small
    
      ; plot config
      yrange2 = yrange_rate*3
      if accModes[j] eq 'all' then yrange2 *= 4
      if accRateModel eq 2 then yrange2 /= 4
      
      if accRateModel eq 2 and accModes[j] eq 'smooth'   then yrange2 = [0.001,20.0]
      if accRateModel eq 2 and accModes[j] eq 'all'      then yrange2 = [0.005,20.0]
      if accRateModel eq 2 and accModes[j] eq 'stripped' then yrange2 = [0.0005,2.0]
      if accRateModel eq 2 and accModes[j] eq 'clumpy'   then yrange2 = [0.003,20.0]
      
      if accRateModel eq 0 and accModes[j] eq 'smooth'   then yrange2 = [0.01,20.0]
      if accRateModel eq 0 and accModes[j] eq 'all'      then yrange2 = [0.1,50.0]
      if accRateModel eq 0 and accModes[j] eq 'stripped' then yrange2 = [0.005,5.0]
      if accRateModel eq 0 and accModes[j] eq 'clumpy'   then yrange2 = [0.05,20.0]
      
      ; subplot config
      pos_sub = [0.58,0.23,0.89,0.48]
      yrange_sub = [0.0,5.0]
      if accRateModel eq 0 then yrange_sub = [0.0,4.0]
      
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=9,/ys,/ylog,yminor=0,$
        ytitle=textoidl("dM^{"+dMdt_top+"}_{gas,"+accModes[j]+" }/dt [_{ }M_{sun }/_{ }yr_{ }]"),$
        xtitle="Redshift",/noerase,pos=pos_single
        
      universeage_axis, xrange_z, yrange2, /ylog
      
      ; loop over runs
      temp_sub = {}
      
      for i=0,n_tags(sP)-1 do begin
        ; loop over mass bins (just one now)
        massStrs = []
        foreach massInd,massInds,k do begin
        
          massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
        
          xvals           = fltarr(n_elements(redshifts))
          yvals_gal_hot   = fltarr(n_elements(redshifts),3)
          yvals_gal_cold  = fltarr(n_elements(redshifts),3)
          yvals_gal_tot   = fltarr(n_elements(redshifts),3)
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals[m] = sP.(i).(m).(j).redshift
            
            yvals_gal_hot[m,*]  = reform( mbv.(i).(m).(j).galaxyRates.galaxy_acc.hot[massInd,*] )
            yvals_gal_cold[m,*] = reform( mbv.(i).(m).(j).galaxyRates.galaxy_acc.cold[massInd,*] )
            yvals_gal_tot[m,*]  = reform( mbv.(i).(m).(j).galaxyRates.galaxy_acc.total[massInd,*] )
          endfor ;m
                          
          ; smooth and set zeros to a plottable value
          for m=0,2 do begin
            yvals_gal_hot[*,m]  = smooth(yvals_gal_hot[*,m]/units.hubbleParam,sK,/nan)
            yvals_gal_cold[*,m] = smooth(yvals_gal_cold[*,m]/units.hubbleParam,sK,/nan)
            yvals_gal_tot[*,m]  = smooth(yvals_gal_tot[*,m]/units.hubbleParam,sK,/nan)
          endfor
        
          ; keep values for this run for subplot
          save_struct = { gal_hot:yvals_gal_hot, gal_cold:yvals_gal_cold, gal_tot:yvals_gal_tot,z:xvals }
          temp_sub = mod_struct( temp_sub, sP.(i).(0).(0).run, save_struct )
        
          ; 25th-75th
          ;if i eq 0 then $
          ;oplotBand,xvals,yvals_gal_hot[*,0],yvals_gal_hot[*,2],color='light gray',yrange=yrange2
        
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal_hot[*,1],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesHC[0],/overplot
          cgPlot,xvals,yvals_gal_cold[*,1],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesHC[1],/overplot
          ;cgPlot,xvals,yvals_gal_tot[*,1],color=simColors[i],psym=psyms[i],symsize=symsize,line=2,/overplot
          
        endforeach ; massBins
        
      endfor ;i
      
      ;legend,textoidl("M_{halo} = "+string(massStrs,format='(f4.1)')),/top,/right
      
      legend,textoidl(["T_{max} > T_{vir,acc}","T_{max} < T_{vir,acc}"]),$
        linestyle=linesHC[0:1],color=cgColor('black'),textcolor=cgColor('black'),/top,/left
        
      if accRateModel eq 2 then $
        legend,simNames,textcolor=simColors,pos=[pos_sub[2]+0.02,pos_sub[3]+0.1],/right,/normal
      if accRateModel eq 0 then $
        legend,simNames,textcolor=simColors,/top,/right
        
      ; subplot
      cgColorFill,pos_sub[[0,2,2,0,0]],pos_sub[[1,1,3,3,1]],/normal,color='WT2'
      
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_sub,/xs,/ys,$
        ytitle=textoidl("noFB_{ }/_{ }FB_{ }"),$
        xtitle="Redshift",/noerase,pos=pos_sub,charsize=!p.charsize-0.3

      if n_tags(temp_sub) ne 2 then print,'WARNING, check'
      
      if n_tags(temp_sub) eq 2 then begin
        cgPlot,xrange_z+[0.3,-0.3],[1.0,1.0],line=0,/overplot,color=cgColor('light gray')
        
        yy = smooth( temp_sub.tracer.gal_hot[*,1] / temp_sub.feedback.gal_hot[*,1] ,sK,/nan)
        cgPlot,xvals,yy,line=linesHC[0],/overplot,color=cgColor('black')
              
        yy = smooth( temp_sub.tracer.gal_cold[*,1] / temp_sub.feedback.gal_cold[*,1] ,sK,/nan)
        cgPlot,xvals,yy,line=linesHC[1],/overplot,color=cgColor('black')
      endif
      
    end_PS
    
    ; plot (4c) - net rates vs redshift, galaxy hot/cold (abcd method)
    start_PS, sP.(0).(0).(0).plotPath + 'netRateVsRedshift.galaxy_abcd.' + plotStrAM + '.eps', /small
    
      ; plot config
      yrange2 = yrange_rate*3
      if accModes[j] eq 'all' then yrange2 *= 4
      if accRateModel eq 2 then yrange2 /= 4
      
      if accRateModel eq 2 and accModes[j] eq 'smooth'   then yrange2 = [0.0005,15.0]
      if accRateModel eq 2 and accModes[j] eq 'all'      then yrange2 = [0.005,20.0]
      if accRateModel eq 2 and accModes[j] eq 'stripped' then yrange2 = [0.0005,2.0]
      if accRateModel eq 2 and accModes[j] eq 'clumpy'   then yrange2 = [0.003,20.0]
      
      if accRateModel eq 0 and accModes[j] eq 'smooth'   then yrange2 = [0.001,15.0]
      if accRateModel eq 0 and accModes[j] eq 'all'      then yrange2 = [0.1,50.0]
      if accRateModel eq 0 and accModes[j] eq 'stripped' then yrange2 = [0.005,5.0]
      if accRateModel eq 0 and accModes[j] eq 'clumpy'   then yrange2 = [0.05,20.0]
      
      ; subplot config
      pos_sub = [0.58,0.23,0.89,0.48]
      if accRateModel eq 2 then yrange_sub = [1.0,200.0]
      if accRateModel eq 0 then yrange_sub = [0.0,8.0]
      
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=9,/ys,/ylog,yminor=0,$
        ytitle=textoidl("dM^{"+dMdt_top+"}_{gas,"+accModes[j]+" }/dt [_{ }M_{sun }/_{ }yr_{ }]"),$
        xtitle="Redshift",/noerase,pos=pos_single
        
      universeage_axis, xrange_z, yrange2, /ylog
      
      ; loop over runs
      temp_sub = {}
      
      for i=0,n_tags(sP)-1 do begin
        ; loop over mass bins (just one now)
        massStrs = []
        foreach massInd,massInds,k do begin
        
          massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
        
          xvals           = fltarr(n_elements(redshifts))
          yvals_gal_hot   = fltarr(n_elements(redshifts),3)
          yvals_gal_cold  = fltarr(n_elements(redshifts),3)
          yvals_gal_tot   = fltarr(n_elements(redshifts),3)
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals[m] = sP.(i).(m).(j).redshift
            
            yvals_gal_hot[m,*]  = reform( mbv.(i).(m).(j).galaxyRates.abcd_gal_net.hot[massInd,*] )
            yvals_gal_cold[m,*] = reform( mbv.(i).(m).(j).galaxyRates.abcd_gal_net.cold[massInd,*] )
            yvals_gal_tot[m,*]  = reform( mbv.(i).(m).(j).galaxyRates.abcd_gal_net.total[massInd,*] )
          endfor ;m
                          
          ; smooth and set zeros to a plottable value
          for m=0,2 do begin
            yvals_gal_hot[*,m]  = smooth(yvals_gal_hot[*,m]  / units.hubbleParam,sK,/nan)
            yvals_gal_cold[*,m] = smooth(yvals_gal_cold[*,m] / units.hubbleParam,sK,/nan)
            yvals_gal_tot[*,m]  = smooth(yvals_gal_tot[*,m]  / units.hubbleParam,sK,/nan)
          endfor
        
          ; keep values for this run for subplot
          save_struct = { gal_hot:yvals_gal_hot, gal_cold:yvals_gal_cold, gal_tot:yvals_gal_tot,z:xvals }
          temp_sub = mod_struct( temp_sub, sP.(i).(0).(0).run, save_struct )
        
          ; 25th-75th
          ;if i eq 0 then $
          ;oplotBand,xvals,yvals_gal_hot[*,0],yvals_gal_hot[*,2],color='light gray',yrange=yrange2
        
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal_hot[*,1],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesHC[0],/overplot
          cgPlot,xvals,yvals_gal_cold[*,1],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesHC[1],/overplot
          ;cgPlot,xvals,yvals_gal_tot[*,1],color=simColors[i],psym=psyms[i],symsize=symsize,line=2,/overplot
          
        endforeach ; massBins
        
      endfor ;i
      
      ;legend,textoidl("M_{halo} = "+string(massStrs,format='(f4.1)')),/top,/right
      
      legend,textoidl(["T_{max} > T_{vir,acc}","T_{max} < T_{vir,acc}"]),$
        linestyle=linesHC[0:1],color=cgColor('black'),textcolor=cgColor('black'),/top,/left
        
      if accRateModel eq 2 then $
        legend,simNames,textcolor=simColors,pos=[pos_sub[2]+0.02,pos_sub[3]+0.1],/right,/normal
      if accRateModel eq 0 then $
        legend,simNames,textcolor=simColors,/top,/right
        
      ; subplot
      cgColorFill,pos_sub[[0,2,2,0,0]],pos_sub[[1,1,3,3,1]],/normal,color='WT2'
      
      if accRateModel eq 2 then yminor = 0
      if accRateModel eq 0 then yminor = 2
      
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_sub,/xs,/ys,$
        ylog=(accRateModel eq 2),yminor=yminor,$
        ytitle=textoidl("noFB_{ }/_{ }FB_{ }"),$
        xtitle="Redshift",/noerase,pos=pos_sub,charsize=!p.charsize-0.3

      if n_tags(temp_sub) ne 2 then print,'WARNING, check'
      
      if n_tags(temp_sub) eq 2 then begin
        if accRateModel eq 0 then yval = 1.0
        if accRateModel eq 2 then yval = 10.0
        
        cgPlot,xrange_z+[0.3,-0.3],[yval,yval],line=0,/overplot,color=cgColor('light gray')
        
        yy = smooth( temp_sub.tracer.gal_hot[*,1] / temp_sub.feedback.gal_hot[*,1] ,sK,/nan)
        cgPlot,xvals,yy,line=linesHC[0],/overplot,color=cgColor('black')
              
        yy = smooth( temp_sub.tracer.gal_cold[*,1] / temp_sub.feedback.gal_cold[*,1] ,sK,/nan)
        cgPlot,xvals,yy,line=linesHC[1],/overplot,color=cgColor('black')
      endif
      
    end_PS
    
    ; plot (5) - net/in/out rates vs redshift, galaxy hot/cold
    if plotUnused then begin
    start_PS, sP.(0).(0).(0).plotPath + 'accOutNetRatesVsRedshift.galaxyOld.' + plotStrAM + '.eps', /small
    
      yrange2 = yrange_rate*[10,2]
      if accModes[j] eq 'all' then yrange2 *= 4
      
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=9,/ys,/ylog,yminor=0,$
        ytitle=textoidl("dM^{"+dMdt_top+"}_{gas,"+accModes[j]+" }/dt [_{ }M_{sun }/_{ }yr_{ }]"),$
        xtitle="Redshift",/noerase,pos=pos_single
        
      universeage_axis, xrange_z, yrange2, /ylog
        
      ; loop over runs      
      for i=0,n_tags(sP)-1 do begin
        ; loop over mass bins (just one now)
        massStrs = []
        foreach massInd,massInds,k do begin
        
          massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
        
          xvals          = []
          yvals_gal_acc  = []
          yvals_gal_out  = []
          yvals_gal_net  = []
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals = [xvals, sP.(i).(m).(j).redshift]
            yvals_gal_acc = [yvals_gal_acc, $
              mbv.(i).(m).(j).galaxyMedian.accRate.total.hot.tVirAcc[tVirInd,massInd] + $
              mbv.(i).(m).(j).galaxyMedian.accRate.total.cold.tVirAcc[tVirInd,massInd]]
            yvals_gal_out = [yvals_gal_out, $
              mbv.(i).(m).(j).galaxyMedian.outRate.total.hot.tVirAcc[tVirInd,massInd] + $
              mbv.(i).(m).(j).galaxyMedian.outRate.total.cold.tVirAcc[tVirInd,massInd]]
            yvals_gal_net = [yvals_gal_net, $
              mbv.(i).(m).(j).galaxyMedian.netRate.total.hot.tVirAcc[tVirInd,massInd] + $
              mbv.(i).(m).(j).galaxyMedian.netRate.total.cold.tVirAcc[tVirInd,massInd]]
          endfor ;m
        
          ; smooth and set zeros to a plottable value
          yvals_gal_acc = smooth(yvals_gal_acc/units.hubbleParam,sK,/nan)
          yvals_gal_out = smooth(yvals_gal_out/units.hubbleParam,sK,/nan)
          yvals_gal_net = smooth(yvals_gal_net/units.hubbleParam,sK,/nan)
          
          temp = mod_struct( temp, 'ratio_gal_'+sP.(i).(0).(0).run+'_'+accModes[j], (yvals_gal_out/yvals_gal_acc) )
        
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal_acc,color=simColors[i],psym=psyms[i],symsize=symsize,line=linesIO[0],/overplot
          cgPlot,xvals,yvals_gal_out,color=simColors[i],psym=psyms[i],symsize=symsize,line=linesIO[1],/overplot
          cgPlot,xvals,yvals_gal_net,color=simColors[i],psym=psyms[i],symsize=symsize,line=linesIO[2],/overplot
          
        endforeach ; massBins
        
      endfor ;i
      
      ;legend,textoidl("M_{halo} = "+string(massStrs,format='(f4.1)')),/top,/right
      legend,textoidl(["acc","out","net"]),$
        linestyle=linesIO,color=cgColor('black'),textcolor=cgColor('black'),/bottom,/right
      legend,simNames,textcolor=simColors,/top,/left
    
    end_PS
    endif ;plotUnused
    
    ; plot (5b) - net/in/out rates vs redshift, galaxy hot/cold (alternative definitions)
    start_PS, sP.(0).(0).(0).plotPath + 'accOutNetRatesVsRedshift.galaxy.' + plotStrAM + '.eps', /small
    
      yrange2 = yrange_rate*[20,0.7]
      if accModes[j] eq 'all' then yrange2 *= [4,3.8]
      if accRateModel eq 0 then yrange2 *= [20,1.5]
      
      if accRateModel eq 2 and accModes[j] eq 'smooth'   then yrange2 = [0.05,20.0]
      if accRateModel eq 2 and accModes[j] eq 'all'      then yrange2 = [0.1,40.0]
      if accRateModel eq 2 and accModes[j] eq 'stripped' then yrange2 = [0.01,5.0]
      if accRateModel eq 2 and accModes[j] eq 'clumpy'   then yrange2 = [0.1,40.0]
      
      if accRateModel eq 0 and accModes[j] eq 'smooth'   then yrange2 = [0.8,20.0]
      if accRateModel eq 0 and accModes[j] eq 'all'      then yrange2 = [1.0,100.0]
      if accRateModel eq 0 and accModes[j] eq 'stripped' then yrange2 = [0.5,10.0]
      if accRateModel eq 0 and accModes[j] eq 'clumpy'   then yrange2 = [0.5,50.0]
      
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=9,/ys,/ylog,yminor=0,$
        ytitle=textoidl("dM^{"+dMdt_top+"}_{gas,"+accModes[j]+" }/dt [_{ }M_{sun }/_{ }yr_{ }]"),$
        xtitle="Redshift",/noerase,pos=pos_single
        
      universeage_axis, xrange_z, yrange2, /ylog
        
      ; loop over runs      
      for i=n_tags(sP)-1,0,-1 do begin
        ;print,sP.(i).(0).(0).run + ' ' + accModes[j]
        
        ; loop over mass bins (just one now)
        massStrs = []
        foreach massInd,massInds,k do begin
        
          massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
        
          xvals          = fltarr(n_elements(redshifts))
          yvals_gal_acc  = fltarr(n_elements(redshifts),3)
          yvals_gal_out  = fltarr(n_elements(redshifts),3)
          yvals_gal_net  = fltarr(n_elements(redshifts),3)
           
          yvals_gal_acc2  = fltarr(n_elements(redshifts),3)
          yvals_gal_out2  = fltarr(n_elements(redshifts),3)
           
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals[m] = sP.(i).(m).(j).redshift
            
            halo_acc   = reform( mbv.(i).(m).(j).galaxyRates.halo_acc.total[massInd,*] )
            galaxy_acc = reform( mbv.(i).(m).(j).galaxyRates.galaxy_acc.total[massInd,*] )
            halo_out   = reform( mbv.(i).(m).(j).galaxyRates.halo_out.total[massInd,*] )
            galaxy_out = reform( mbv.(i).(m).(j).galaxyRates.galaxy_out.total[massInd,*] )
            
            yv_gal_out = reform( mbv.(i).(m).(j).galaxyRates.yvals_gal_out.total[massInd,*] )
            yv_gal_acc = reform( mbv.(i).(m).(j).galaxyRates.yvals_gal_acc.total[massInd,*] )
            
            ; NET = sum of accRate of material still inside the galaxy (#1)
            yvals_gal_net[m,*] = galaxy_acc
              
            ; OUT = sum of outRate of material now outside the galaxy (#2+#3)
            yvals_gal_out[m,*] = yv_gal_out ; halo_out + galaxy_out ;(but add prior to median)
            
            ; ACC = sum of net and out (#1+#2+#3)
            ;yvals_gal_acc[m,*] = yvals_gal_net[-1] + yvals_gal_out[-1]
            ; ACC = sum of net and out, but #3 is a strict subset of #1 (so #1+#2)
            yvals_gal_acc[m,*] = yv_gal_acc ; galaxy_acc + halo_out ;(but add prior to median)
            ; ACC = sum of accRate_galaxy and accRate_halo (#1+#4)
            ;yvals_gal_acc[m,*] = galaxy_acc + halo_acc
            
            ;print,galaxy_acc[1]+halo_acc[1]-halo_out[1]-galaxy_out[1],galaxy_acc[1]
            
          endfor ;m
        
          ; smooth and set zeros to a plottable value
          for m=0,2 do begin
            yvals_gal_acc[*,m] = smooth(yvals_gal_acc[*,m]/units.hubbleParam,3,/nan)
            yvals_gal_out[*,m] = smooth(yvals_gal_out[*,m]/units.hubbleParam,3,/nan)
            yvals_gal_net[*,m] = smooth(yvals_gal_net[*,m]/units.hubbleParam,3,/nan)
          endfor
          
          temp = mod_struct( temp, 'ratio_gal_'+sP.(i).(0).(0).run+'_'+accModes[j], (yvals_gal_out/yvals_gal_acc) )
        
          ; plot this combination of: run,massbin (vs redshift)
          
          ; 25th-75th percentile (color band)
          if i eq 0 then tColor = 'light gray'
          if i eq 1 then tColor = 'light gray'
          if i eq 1 then $
          oplotBand,xvals,yvals_gal_net[*,0],yvals_gal_net[*,2],color=tColor,yrange=yrange2
          
          ; 25th-75th percentile (adjacent lines of lesser thickness)
          ;cgPlot,xvals,yvals_gal_acc[*,0],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesIO[0],thick=!p.thick-1.0,/overplot
          ;cgPlot,xvals,yvals_gal_acc[*,2],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesIO[0],thick=!p.thick-1.0,/overplot
          
          ; median
          cgPlot,xvals,yvals_gal_acc[*,1],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesIO[0],/overplot
          cgPlot,xvals,yvals_gal_out[*,1],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesIO[1],/overplot
          cgPlot,xvals,yvals_gal_net[*,1],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesIO[2],/overplot
          
          ; 25th-75th percentile (errorbars)
          ;for m=0,n_elements(xvals)-1 do begin
          ;  cgPlot,[xvals[m],xvals[m]],[yvals_gal_net[m,0],yvals_gal_net[m,2]],$
          ;    color=simColors[i],line=0,/overplot
          ;endfor
          
          ; redshift scaling
          if accRateModel eq 2 and accModes[j] eq 'smooth' then begin
            xx = [0.5,2.5]
            yy = xx^1.0 * 4.8
            cgPlot,xx,yy,line=0,color='gray',/overplot
            cgText,xx[0],yy[0]+1.3,textoidl("\propto z^1"),color='gray'
            
            xx = [1.0,3.5]
            yy = xx^0.5 * 0.3
            cgPlot,xx,yy,line=0,color='gray',/overplot
            cgText,xx[0],yy[0]-0.1,textoidl("\propto z^{1/2}"),color='gray'
          endif
          
        endforeach ; massBins
        
      endfor ;i
      
      ;legend,textoidl("M_{halo} = "+string(massStrs,format='(f4.1)')),/top,/right
      legend,textoidl(["acc","out","net"]),$
        linestyle=linesIO,color=cgColor('black'),textcolor=cgColor('black'),/bottom,/right
      legend,simNames,textcolor=simColors,/top,/left
    
    end_PS
    
    ; plot (5c) - net/in/out rates vs redshift, galaxy hot/cold (abcd)
    start_PS, sP.(0).(0).(0).plotPath + 'accOutNetRatesVsRedshift.galaxy_abcd.' + plotStrAM + '.eps', /small
    
      yrange2 = yrange_rate*[20,0.7]
      if accModes[j] eq 'all' then yrange2 *= [4,3.8]
      if accRateModel eq 0 then yrange2 *= [20,1.5]
      
      if accRateModel eq 2 and accModes[j] eq 'smooth'   then yrange2 = [0.05,20.0]
      if accRateModel eq 2 and accModes[j] eq 'all'      then yrange2 = [0.1,40.0]
      if accRateModel eq 2 and accModes[j] eq 'stripped' then yrange2 = [0.01,5.0]
      if accRateModel eq 2 and accModes[j] eq 'clumpy'   then yrange2 = [0.1,40.0]
      
      if accRateModel eq 0 and accModes[j] eq 'smooth'   then yrange2 = [0.8,20.0]
      if accRateModel eq 0 and accModes[j] eq 'all'      then yrange2 = [1.0,100.0]
      if accRateModel eq 0 and accModes[j] eq 'stripped' then yrange2 = [0.5,10.0]
      if accRateModel eq 0 and accModes[j] eq 'clumpy'   then yrange2 = [0.5,50.0]
      
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange2,xs=9,/ys,/ylog,yminor=0,$
        ytitle=textoidl("dM^{"+dMdt_top+"}_{gas,"+accModes[j]+" }/dt [_{ }M_{sun }/_{ }yr_{ }]"),$
        xtitle="Redshift",/noerase,pos=pos_single
        
      universeage_axis, xrange_z, yrange2, /ylog
        
      ; loop over runs      
      for i=n_tags(sP)-1,0,-1 do begin
        ;print,sP.(i).(0).(0).run + ' ' + accModes[j]
        
        ; loop over mass bins (just one now)
        massStrs = []
        foreach massInd,massInds,k do begin
        
          massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
        
          xvals          = fltarr(n_elements(redshifts))
          yvals_gal_acc  = fltarr(n_elements(redshifts),3)
          yvals_gal_out  = fltarr(n_elements(redshifts),3)
          yvals_gal_net  = fltarr(n_elements(redshifts),3)
           
          yvals_gal_acc2  = fltarr(n_elements(redshifts),3)
          yvals_gal_out2  = fltarr(n_elements(redshifts),3)
           
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals[m] = sP.(i).(m).(j).redshift
                        
            ; NET = #1+#4 - (#2+#3)
            yvals_gal_net[m,*] = reform( mbv.(i).(m).(j).galaxyRates.abcd_gal_net.total[massInd,*] )
              
            ; OUT = #3+#2
            yvals_gal_out[m,*] = reform( mbv.(i).(m).(j).galaxyRates.abcd_gal_out.total[massInd,*] )
            
            ; ACC = #1+#4
            yvals_gal_acc[m,*] = reform( mbv.(i).(m).(j).galaxyRates.abcd_gal_acc.total[massInd,*] )
            
          endfor ;m
        
          ; smooth and set zeros to a plottable value
          for m=0,2 do begin
            yvals_gal_acc[*,m] = smooth(yvals_gal_acc[*,m]/units.hubbleParam,3,/nan)
            yvals_gal_out[*,m] = smooth(yvals_gal_out[*,m]/units.hubbleParam,3,/nan)
            yvals_gal_net[*,m] = smooth(yvals_gal_net[*,m]/units.hubbleParam,3,/nan)
          endfor
          
          ; plot this combination of: run,massbin (vs redshift)
          
          ; 25th-75th percentile (color band)
          if i eq 0 then tColor = 'light gray'
          if i eq 1 then tColor = 'light gray'
          if i eq 1 then $
          oplotBand,xvals,yvals_gal_net[*,0],yvals_gal_net[*,2],color=tColor,yrange=yrange2
          
          ; 25th-75th percentile (adjacent lines of lesser thickness)
          ;cgPlot,xvals,yvals_gal_acc[*,0],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesIO[0],thick=!p.thick-1.0,/overplot
          ;cgPlot,xvals,yvals_gal_acc[*,2],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesIO[0],thick=!p.thick-1.0,/overplot
          
          ; median
          cgPlot,xvals,yvals_gal_acc[*,1],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesIO[0],/overplot
          cgPlot,xvals,yvals_gal_out[*,1],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesIO[1],/overplot
          cgPlot,xvals,yvals_gal_net[*,1],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesIO[2],/overplot
          
          ; 25th-75th percentile (errorbars)
          ;for m=0,n_elements(xvals)-1 do begin
          ;  cgPlot,[xvals[m],xvals[m]],[yvals_gal_net[m,0],yvals_gal_net[m,2]],$
          ;    color=simColors[i],line=0,/overplot
          ;endfor
          
          ; redshift scaling
          if accRateModel eq 2 and accModes[j] eq 'smooth' then begin
            xx = [0.5,2.5]
            yy = xx^1.0 * 4.8
            cgPlot,xx,yy,line=0,color='gray',/overplot
            cgText,xx[0],yy[0]+1.3,textoidl("\propto z^1"),color='gray'
            
            xx = [1.0,3.5]
            yy = xx^0.5 * 0.3
            cgPlot,xx,yy,line=0,color='gray',/overplot
            cgText,xx[0],yy[0]-0.1,textoidl("\propto z^{1/2}"),color='gray'
          endif
          
        endforeach ; massBins
        
      endfor ;i
      
      ;legend,textoidl("M_{halo} = "+string(massStrs,format='(f4.1)')),/top,/right
      legend,textoidl(["acc","out","net"]),$
        linestyle=linesIO,color=cgColor('black'),textcolor=cgColor('black'),/bottom,/right
      legend,simNames,textcolor=simColors,/top,/left
    
    end_PS
    
    ; plot (6) - net cold fractions vs redshift, galaxy hot/cold
    if plotUnused then begin
    start_PS, sP.(0).(0).(0).plotPath + 'coldFracNetVsRedshift.galaxy-haloOld.' + plotStrAM + '.eps', /small
    
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_frac,xs=5,ys=5,/ylog,/noerase,$
        pos=pos_single
        
      cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
        [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
      cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
        [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
    
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_frac,xs=9,/ys,/ylog,yminor=0,$
        ytitle=textoidl("(dM/dt)^{"+dMdt_top+"}_{T_{max}<T_{vir,acc}} / (dM/dt)^{"+dMdt_top+"}_{total,"+$
               accModes[j]+" }"),$
        xtitle="Redshift",/noerase,pos=pos_single
        
      universeage_axis, xrange_z, yrange_frac, /ylog
        
      ; loop over runs      
      for i=0,n_tags(sP)-1 do begin
        ; loop over mass bins (just one now)
        massStrs = []
        foreach massInd,massInds,k do begin
        
          massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
        
          xvals      = []
          yvals_gal  = []
          yvals_halo = []
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals = [xvals, sP.(i).(m).(j).redshift]
            
            yvals_gal  = [yvals_gal,  mbv.(i).(m).(j).coldfrac_net.gal[massInd]]
            yvals_halo = [yvals_halo, mbv.(i).(m).(j).coldfrac_net.halo[massInd]]
          endfor ;m
        
          ; smooth and set zeros to a plottable value
          yvals_gal  = smooth(yvals_gal,sK,/nan)
          yvals_halo = smooth(yvals_halo,sK,/nan)
        
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal,color=simColors[i],psym=psyms[i],symsize=symsize,line=linesHC[0],/overplot
          cgPlot,xvals,yvals_halo,color=simColors[i],psym=psyms[i],symsize=symsize,line=linesHC[1],/overplot
          
        endforeach ; massBins
        
      endfor ;i
      
      ;legend,textoidl("M_{halo} = "+string(massStrs,format='(f4.1)')),/top,/right
      legend,['galaxy','halo'],linestyle=linesHC[0:1],$
        color=cgColor('black'),textcolors=['black','black'],linesize=0.4,box=0,/bottom,/right
      legend,simNames,textcolor=simColors,/top,/right
  
    end_PS
    endif ;plotUnused
    
    ; plot (6b) - net cold fractions vs redshift, galaxy hot/cold (updated)
    start_PS, sP.(0).(0).(0).plotPath + 'coldFracNetVsRedshift.galaxy-halo.' + plotStrAM + '.eps', /small
    
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_frac,xs=5,ys=5,/ylog,/noerase,$
        pos=pos_single
        
      cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
        [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
      cgPolygon,[xrange_z[0],xrange_z[1],xrange_z[1],xrange_z[0],xrange_z[0]],$
        [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
    
      cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_frac,xs=9,/ys,/ylog,yminor=0,$
        ytitle=textoidl("(dM/dt)^{"+dMdt_top+"}_{T_{max}<T_{vir,acc}} / (dM/dt)^{"+dMdt_top+"}_{total,"+$
               accModes[j]+" }"),$
        xtitle="Redshift",/noerase,pos=pos_single
        
      universeage_axis, xrange_z, yrange_frac, /ylog
        
      ; loop over runs      
      for i=0,n_tags(sP)-1 do begin
        ; loop over mass bins (just one now)
        massStrs = []
        foreach massInd,massInds,k do begin
        
          massStrs = [massStrs, mbv.(i).(0).(0).logMassBinCen[massInd]]
        
          xvals      = fltarr(n_elements(redshifts))
          yvals_gal  = fltarr(n_elements(redshifts),3)
          yvals_halo = fltarr(n_elements(redshifts),3)
        
          ; loop over redshifts
          for m=0,n_tags(sP.(i))-1 do begin
            xvals[m] = sP.(i).(m).(j).redshift
            
            yvals_gal[m,*]  = mbv.(i).(m).(j).coldFractionsNet.galaxy[massInd,*]
            yvals_halo[m,*] = mbv.(i).(m).(j).coldFractionsNet.halo[massInd,*]
          endfor ;m
        
          ; smooth and set zeros to a plottable value
          for m=0,2 do begin
            yvals_gal[*,m]  = smooth(yvals_gal[*,m],sK,/nan)
            yvals_halo[*,m] = smooth(yvals_halo[*,m],sK,/nan)
          endfor
        
          ; plot this combination of: run,massbin (vs redshift)
          cgPlot,xvals,yvals_gal[*,1],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesHC[0],/overplot
          cgPlot,xvals,yvals_halo[*,1],color=simColors[i],psym=psyms[i],symsize=symsize,line=linesHC[1],/overplot
          
          ; 25th-75th percentile (errorbars)
          ;for m=0,n_elements(xvals)-1 do begin
          ;  tFac = 1.0
          ;  xOffset = 0.01 * ( (i eq 0)*2 - 1.0 )
          ;  
          ;  cgPlot,[xvals[m],xvals[m]]+xOffset,[yvals_gal[m,0],yvals_gal[m,2]],$
          ;    color=simColors[i],line=linesHC[0],thick=!p.thick-tFac,/overplot
          ;  cgPlot,[xvals[m],xvals[m]]+xOffset,[yvals_halo[m,0],yvals_halo[m,2]],$
          ;    color=simColors[i],line=linesHC[1],thick=!p.thick-tFac,/overplot
          ;endfor
          
        endforeach ; massBins
        
      endfor ;i
      
      ;legend,textoidl("M_{halo} = "+string(massStrs,format='(f4.1)')),/top,/right
      legend,['galaxy','halo'],linestyle=linesHC[0:1],$
        color=cgColor('black'),textcolors=['black','black'],linesize=0.4,box=0,/bottom,/right
      legend,simNames,textcolor=simColors,/top,/right
  
    end_PS

  endfor ;accModes,j
     
  ; plots vs halo mass (all accModes together, one plot per redshift)
  ; ------------------------------------------------------------
  foreach zInd,zInds,k do begin
  
    ; plot (7) - fraction of accretion by mode, galaxy hot/cold
    if plotUnused then begin
    start_PS, sP.(0).(0).(0).plotPath + $
      'rateFracs4Redshifts.zInd-'+str(zInd)+'.galaxy.'+plotStr+'-allModes.eps', /extrabig
      
      pos = plot_pos(rows=2,cols=2,/gap)
      
      ; loop over accretion modes
      for m=1,n_elements(accModes)-1 do begin
        
        cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_frac,xs=5,ys=5,/ylog,/noerase,$
          pos=pos[m-1]-[0,0.02,0,0.04]
          
        cgPolygon,[xrange_halo[0],xrange_halo[1],xrange_halo[1],xrange_halo[0],xrange_halo[0]],$
          [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
        
        cgPolygon,[xrange_halo[0],xrange_halo[1],xrange_halo[1],xrange_halo[0],xrange_halo[0]],$
          [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
          
        cgPlot,[px_min,px_min],[yrange_frac[0],yrange_frac[1]],line=2,/overplot,color=cgColor(px_cl3)
        cgPlot,[px_max,px_max],[yrange_frac[0],yrange_frac[1]],line=2,/overplot,color=cgColor(px_cl3)
        
        barStr = str(round(py_max2*100))+'% - '+str(round(py_max*100))+'%'
        if m eq 4 then $
          cgText, xrange_halo[0]+0.6, py_min+0.39, barStr, alignment=0.5, color=cgColor('gray')
          
        barStr = str(round(py_min*100))+'% - '+str(round(py_max2*100))+'% '
        if m eq 4 then $
          cgText, xrange_halo[0]+0.6, py_min+0.06, barStr, alignment=0.5, color=cgColor('gray')
        
        cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_frac,/xs,/ys,/ylog,yminor=0,$
          ytitle=textoidl("(dM/dt)^{"+dMdt_top+"}_{"+accModes[m]+"} / (dM/dt)^{"+dMdt_top+"}_{total}"),$
          xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),$
          /noerase,pos=pos[m-1]-[0,0.02,0,0.04]
      
        for i=0,n_tags(sP)-1 do begin
          if m ge n_tags(mbv.(i).(zInd)) then continue ; mode does not exist for this run
          
          ; gal hot
          xvals = alog10( 10.0^mbv.(i).(zind).(0).logMassBinCen / units.hubbleParam )
          yvals = smooth( mbv.(i).(zInd).(m).mode_ratio.gal_hot, sK, /nan )
          cgPlot,xvals,yvals,color=sP.(i).(zind).(0).colors[cInd],line=linesHC[0],/overplot
          
          ; gal cold
          xvals = alog10( 10.0^mbv.(i).(zind).(0).logMassBinCen / units.hubbleParam )
          yvals = smooth( mbv.(i).(zInd).(m).mode_ratio.gal_cold, sK, /nan )
          cgPlot,xvals,yvals,color=sP.(i).(zind).(0).colors[cInd],line=linesHC[1],/overplot
            
        endfor
        
        ; panel-specific legends
        if m eq 0 then $
          legend,textoidl('z_{ }=_{ }'+string(redshifts[zind],format='(f3.1)')),/top,/right
        if m eq 2 then $
          legend,textoidl(["T_{max} > T_{vir,acc}","T_{max} < T_{vir,acc}"]),$
            linestyle=linesHC[0:1],color=cgColor('black'),textcolor=cgColor('black'),/bottom,/left
        if m eq 3 then $
          legend,simNames,textcolors=simColors,box=0,/top,/left
          
      endfor ;m

    end_PS
    endif ;plotUnused
  endforeach
       
  ; plots vs halo mass (all 2x2 redshifts together, one plot per accMode)
  ; ------------------------------------------------------------
  for j=0,n_tags(sP.(1).(0))-1 do begin
    pos = plot_pos(total=n_elements(zInds),/gap) ; plot positioning (3x2, 2x2, or 1x2 with gaps)
    
    plotStrAM = plotStr + '_am-'+accModes[j]
  
    ; plot (8) - net cold fractions, galaxy/halo
    if plotUnused then begin
    start_PS, sP.(0).(0).(0).plotPath + 'coldFracNet4Redshifts.galaxy-haloOld.' + plotStrAM + '.eps', /extrabig
      
      foreach zInd,zInds,k do begin
        cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_frac,xs=5,ys=5,$
          /noerase,pos=pos[k],/ylog,yminor=0 ; suppress all axes, only establish coordinate system
      
        cgPolygon,[xrange_halo[0],xrange_halo[1],xrange_halo[1],xrange_halo[0],xrange_halo[0]],$
          [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
        cgPolygon,[xrange_halo[0],xrange_halo[1],xrange_halo[1],xrange_halo[0],xrange_halo[0]],$
          [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
          
        cgPlot,[px_min,px_min],[yrange_frac[0],yrange_frac[1]],line=2,/overplot,color=cgColor(px_cl3)
        cgPlot,[px_max,px_max],[yrange_frac[0],yrange_frac[1]],line=2,/overplot,color=cgColor(px_cl3)
          
        cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_frac,/xs,/ys,/ylog,yminor=0,$
          ytitle=textoidl("(dM/dt)^{"+dMdt_top+"}_{T_{max}<T_{vir,acc}} / (dM/dt)^{"+dMdt_top+$
                 "}_{total,"+accModes[j]+" }"),$
          xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),/noerase,pos=pos[k]
      
        for i=0,n_tags(sP)-1 do begin
          xvals = alog10( 10.0^mbv.(i).(zInd).(j).logMassBinCen / units.hubbleParam )
          
          yvals = smooth(mbv.(i).(zInd).(j).coldfrac_net.gal,sK,/nan) < 1.0
          cgPlot,xvals,yvals,color=simColors[i],line=linesHC[0],/overplot ; gal
          
          yvals = smooth(mbv.(i).(zInd).(j).coldfrac_net.halo,sK,/nan) < 1.0
          cgPlot,xvals,yvals,color=simColors[i],line=linesHC[1],/overplot ; halo
        endfor
        
        ; add gadget line at z=2
        if k eq 1 and 0 then begin
          gaColor = sP.sP2_gadget.(0).(0).colors[cInd]
          
          xvals = alog10( 10.0^mbv.mbv2_Gadget.redshift6_200.(j).logMassBinCen / units.hubbleParam )
          
          yvals = smooth(mbv.mbv2_Gadget.redshift6_200.(j).coldfrac_net.gal,sK,/nan) < 1.0
          cgPlot,xvals,yvals,color=gaColor,line=linesHC[0],/overplot ; gal
          
          yvals = smooth(mbv.mbv2_Gadget.redshift6_200.(j).coldfrac_net.halo,sK,/nan) < 1.0
          cgPlot,xvals,yvals,color=gaColor,line=linesHC[1],/overplot ; halo
          
          legend,['GADGET'],textcolors=[gaColor],box=0,/top,/left
        endif
        
        ; redshift legend
        legend,textoidl('z_{ }=_{ }'+string(redshifts[zInd],format='(f3.1)')),$
          pos=[xrange_halo[1]-0.6,yrange_frac[1]-0.25]
          
        ; panel-specific legends
        if k eq 2 then $
          legend,['galaxy','halo'],linestyle=linesHC[0:1],$
          color=cgColor('black'),textcolors=['black','black'],linesize=0.4,box=0,/top,/left
          
        if k eq 3 then $
          legend,simNames,textcolors=simColors,box=0,/top,/left
      endforeach

    end_PS
    endif ;plotUnused
    
    ; plot (8b) - net cold fractions, galaxy/halo (updated)
    start_PS, sP.(0).(0).(0).plotPath + 'coldFracNet4Redshifts.galaxy-halo.' + plotStrAM + '.eps', /extrabig
      
      foreach zInd,zInds,k do begin
        cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_frac,xs=5,ys=5,$
          /noerase,pos=pos[k],/ylog,yminor=0 ; suppress all axes, only establish coordinate system
      
        cgPolygon,[xrange_halo[0],xrange_halo[1],xrange_halo[1],xrange_halo[0],xrange_halo[0]],$
          [py_min,py_min,py_max,py_max,py_min],color=cgColor(py_cl),/fill
        cgPolygon,[xrange_halo[0],xrange_halo[1],xrange_halo[1],xrange_halo[0],xrange_halo[0]],$
          [py_min2,py_min2,py_max2,py_max2,py_min2],color=cgColor(py_cl2),/fill
          
        cgPlot,[px_min,px_min],[yrange_frac[0],yrange_frac[1]],line=2,/overplot,color=cgColor(px_cl3)
        cgPlot,[px_max,px_max],[yrange_frac[0],yrange_frac[1]],line=2,/overplot,color=cgColor(px_cl3)
          
        cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange_frac,/xs,/ys,/ylog,yminor=0,$
          ytitle=textoidl("(dM/dt)^{"+dMdt_top+"}_{T_{max}<T_{vir,acc}} / (dM/dt)^{"+dMdt_top+$
                 "}_{total,"+accModes[j]+" }"),$
          xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),/noerase,pos=pos[k]
      
        for i=0,n_tags(sP)-1 do begin
          xvals = alog10( 10.0^mbv.(i).(zInd).(j).logMassBinCen / units.hubbleParam )
          
          yvals = smooth(mbv.(i).(zInd).(j).coldFractionsNet.galaxy,sK,/nan) ;< 1.0
          cgPlot,xvals,yvals[*,1],color=simColors[i],line=linesHC[0],/overplot ; gal
          
          yvals = smooth(mbv.(i).(zInd).(j).coldFractionsNet.halo,sK,/nan) ;< 1.0
          cgPlot,xvals,yvals[*,1],color=simColors[i],line=linesHC[1],/overplot ; halo
        endfor
                
        ; redshift legend
        legend,textoidl('z_{ }=_{ }'+string(redshifts[zInd],format='(f3.1)')),$
          pos=[xrange_halo[1]-0.6,yrange_frac[1]-0.25]
          
        ; panel-specific legends
        if k eq 2 then $
          legend,['galaxy','halo'],linestyle=linesHC[0:1],$
          color=cgColor('black'),textcolors=['black','black'],linesize=0.4,box=0,/top,/left
          
        if k eq 3 then $
          legend,simNames,textcolors=simColors,box=0,/top,/left
      endforeach

    end_PS
    
    ; plot (9) - specific net accretion rates, galaxy hot/cold
    if plotUnused then begin
    start_PS, sP.(0).(0).(0).plotPath + 'netSpecificRate4Redshifts.galaxy.' + plotStrAM + '.eps', /extrabig
      
      foreach zInd,zInds,k do begin
        yrange2 = yrange_spec - 0.5*k
        
        cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange2,xs=5,ys=5,$
          /noerase,pos=pos[k] ; suppress all axes, only establish coordinate system
      
        cgPolygon,[px_min,px_max,px_max,px_min,px_min],$
          [yrange2[0],yrange2[0],yrange2[1],yrange2[1],yrange2[0]],$
          color=cgColor(px_cl),/fill
          
        cgPlot,[0],[0],/nodata,xrange=xrange_halo,yrange=yrange2,/xs,/ys,$
          ytitle=textoidl("dM^{"+dMdt_top+"}_{gas,"+accModes[j]+$
                 "}/dt  M_{halo}^{-1} [_{ }log Gyr^{-1 }]"),$
          xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),/noerase,pos=pos[k]
      
        for i=0,n_tags(sP)-1 do begin
          xvals = alog10( 10.0^mbv.(i).(zInd).(j).logMassBinCen / units.hubbleParam )
          
          yvals = mbv.(i).(zInd).(j).specific_net.gal_hot
          yvals = alog10( smooth(yvals,sK,/nan) )
          cgPlot,xvals,yvals,color=simColors[i],line=linesHC[0],/overplot ; galaxy, hot
          
          yvals = mbv.(i).(zInd).(j).specific_net.gal_cold
          yvals = alog10( smooth(yvals,sK,/nan) )
          cgPlot,xvals,yvals,color=simColors[i],line=linesHC[1],/overplot ; galaxy, cold
        endfor
        
        ; redshift legend
        legend,textoidl('z_{ }=_{ }'+string(redshifts[zind],format='(f3.1)')),/top,/right
        
        ; panel-specific legends
        if k eq 0 then $
          legend,textoidl(['T_{max} > T_{vir,acc}','T_{max} < T_{vir,acc}']),$
          linestyle=linesHC[0:1],color=cgColor('black'),textcolors=replicate('black',2),box=0,/top,/left
          
        if k eq 3 then $
          legend,simNames,textcolors=simColors,box=0,/top,/left
      endforeach

    end_PS
    endif ;plotUnused
    
  endfor ;accModes,j

  ; wind comparison plots
  ; ---------------------
  if accRateModel eq 0 or plotUnused then begin
  ; plot (10) - vs halo mass at one redshift
  foreach zInd,zInds do begin
  
    sMasses = codeMassToLogMsun( logMsunToCodeMass( mbv.(0).(0).(0).logMassBinCen ) / units.hubbleParam )
    
    start_PS,sP.(0).(0).(0).plotPath + 'wind_scalings_' + plotStr + '_zInd'+str(zInd)+'.eps', $
      xs=10.0, ys=4.0
    
      ; plot config
      xrange = xrange_halo + [0.0,1.0]
      tColor = [150,150,150]
      bColor = [150,50,50]
      yrange2 = [1,10000]
      yrange3 = [0.6,1.05]
      
      legendStrs   = [sP.(0).(0).(0).simName,$
                      sP.(1).(0).(0).simName,$
                      'Illustris \eta_w',$
                      'Illustris v_w [km/s]']
      legendColors = [sP.(0).(0).(0).colors[cInd],sP.(1).(0).(0).colors[cInd],$
                      getColor24(tColor),getColor24(bColor)]

      plot_pos = plot_pos(rows=1,cols=2,/gap)
          
      ; plot (a) - eta_w
      cgPlot,[0],[0],/nodata,xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),$
        ytitle=textoidl("\eta_w or v_w"),$
        xrange=xrange,/xs,/ylog,yminor=0,yrange=yrange2,pos=plot_pos[0],/noerase
        
      cgPlot,wind_model.vs_M.(zInd).halo_mass,wind_model.vs_M.(zInd).wind_eta,$
        line=0,thick=!p.thick*4,color=cgColor(tColor),/overplot
        
      cgPlot,[px_min,px_min],[yrange2[0],yrange2[1]],line=2,/overplot,color=cgColor(px_cl3)
      cgPlot,[px_max,px_max],[yrange2[0],yrange2[1]],line=2,/overplot,color=cgColor(px_cl3)
      
      cgPlot,wind_model.vs_M.(zInd).halo_mass,wind_model.vs_M.(zInd).wind_vel,$
        line=0,thick=!p.thick*2,color=cgColor(bColor),/overplot
      
      foreach run,runs,i do begin
        cgPlot,mbv.(i).(zInd).mode0_all.eta_w_xvals,mbv.(i).(zInd).mode0_all.eta_w_unbinned,$
          psym=3,color=sP.(i).(0).(0).colors[cInd],/overplot
        
        w = where(mbv.(i).(zInd).mode0_all.eta_w_binned eq 0.0,count)
        if count gt 0 then mbv.(i).(zInd).mode0_all.eta_w_binned[w] = !values.f_nan
        
        cgPlot,sMasses,smooth(mbv.(i).(zInd).mode0_all.eta_w_binned[*,1],sK,/nan),$
          line=0,color=sP.(i).(0).(0).colors[cInd],/overplot
          
        for j=0,2 do $
          cgPlot,sMasses,smooth(mbv.(i).(zInd).mode0_all.eta_w_binned[*,j],sK,/nan),$
            line=2,color=sP.(i).(0).(0).colors[cInd],/overplot ; 25th/75th percentiles
      endforeach

      legend,textoidl(legendStrs),textcolors=legendColors,/top,/right  
      
      ; plot (b) - eta/(1+eta)
      cgPlot,[0],[0],/nodata,$
        xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),ytitle=textoidl("\eta_w / (1+\eta_w)"),$
        xrange=xrange,/xs,yrange=yrange3,/ys,pos=plot_pos[1],/noerase
        
      cgPlot,xrange,[1.0,1.0],line=0,color=cgColor('light gray'),/overplot
      
      w_eta = wind_model.vs_M.(zInd).wind_eta
      cgPlot,wind_model.vs_M.(zInd).halo_mass,w_eta/(1+w_eta),$
        line=0,thick=!p.thick*4,color=cgColor(tColor),/overplot
      
      cgPlot,[px_min,px_min],[yrange3[0],yrange3[1]],line=2,/overplot,color=cgColor(px_cl3)
      cgPlot,[px_max,px_max],[yrange3[0],yrange3[1]],line=2,/overplot,color=cgColor(px_cl3)
      
      foreach run,runs,i do begin
        yy = mbv.(i).(zInd).mode0_all.eta_w_unbinned
        cgPlot,mbv.(i).(zInd).mode0_all.eta_w_xvals,yy/(1+yy),$
          psym=3,color=sP.(i).(0).(0).colors[cInd],/overplot
        
        yy = mbv.(i).(zInd).mode0_all.eta_w_binned[*,1]
        cgPlot,sMasses,smooth(yy/(1+yy),sK,/nan),$
          line=0,color=sP.(i).(0).(0).colors[cInd],/overplot
          
        for j=0,2 do begin
          yy = mbv.(i).(zInd).mode0_all.eta_w_binned[*,j]
          cgPlot,sMasses,smooth(yy/(1+yy),sK,/nan),$
            line=2,color=sP.(i).(0).(0).colors[cInd],/overplot ; 25th/75th percentiles
        endfor  
      endforeach
           
      legend,['z='+string(redshifts[zInd],format='(f3.1)')],/top,/right
           
    end_PS
    
    ; plot (10b) - abcd method
    start_PS,sP.(0).(0).(0).plotPath + 'wind_scalings_abcd_' + plotStr + '_zInd'+str(zInd)+'.eps', $
      xs=10.0, ys=4.0
    
      ; plot config
      xrange = xrange_halo + [0.0,1.0]
      tColor = [150,150,150]
      bColor = [150,50,50]
      yrange2 = [1,10000]
      yrange3 = [0.6,1.05]
      
      legendStrs   = [sP.(0).(0).(0).simName,$
                      sP.(1).(0).(0).simName,$
                      'Illustris \eta_w',$
                      'Illustris v_w [km/s]']
      legendColors = [sP.(0).(0).(0).colors[cInd],sP.(1).(0).(0).colors[cInd],$
                      getColor24(tColor),getColor24(bColor)]

      plot_pos = plot_pos(rows=1,cols=2,/gap)
          
      ; plot (a) - eta_w
      cgPlot,[0],[0],/nodata,xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),$
        ytitle=textoidl("\eta_w or v_w"),$
        xrange=xrange,/xs,/ylog,yminor=0,yrange=yrange2,pos=plot_pos[0],/noerase
        
      cgPlot,wind_model.vs_M.(zInd).halo_mass,wind_model.vs_M.(zInd).wind_eta,$
        line=0,thick=!p.thick*4,color=cgColor(tColor),/overplot
        
      cgPlot,[px_min,px_min],[yrange2[0],yrange2[1]],line=2,/overplot,color=cgColor(px_cl3)
      cgPlot,[px_max,px_max],[yrange2[0],yrange2[1]],line=2,/overplot,color=cgColor(px_cl3)
      
      cgPlot,wind_model.vs_M.(zInd).halo_mass,wind_model.vs_M.(zInd).wind_vel,$
        line=0,thick=!p.thick*2,color=cgColor(bColor),/overplot
      
      foreach run,runs,i do begin
        cgPlot,mbv.(i).(zInd).mode0_all.eta_w_xvals,mbv.(i).(zInd).mode0_all.eta_abcd_unbinned,$
          psym=3,color=sP.(i).(0).(0).colors[cInd],/overplot
        
        w = where(mbv.(i).(zInd).mode0_all.eta_abcd_binned eq 0.0,count)
        if count gt 0 then mbv.(i).(zInd).mode0_all.eta_abcd_binned[w] = !values.f_nan
        
        cgPlot,sMasses,smooth(mbv.(i).(zInd).mode0_all.eta_abcd_binned[*,1],sK,/nan),$
          line=0,color=sP.(i).(0).(0).colors[cInd],/overplot
          
        for j=0,2 do $
          cgPlot,sMasses,smooth(mbv.(i).(zInd).mode0_all.eta_abcd_binned[*,j],sK,/nan),$
            line=2,color=sP.(i).(0).(0).colors[cInd],/overplot ; 25th/75th percentiles
      endforeach

      legend,textoidl(legendStrs),textcolors=legendColors,/top,/right  
      
      ; plot (b) - eta/(1+eta)
      cgPlot,[0],[0],/nodata,$
        xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),ytitle=textoidl("\eta_w / (1+\eta_w)"),$
        xrange=xrange,/xs,yrange=yrange3,/ys,pos=plot_pos[1],/noerase
        
      cgPlot,xrange,[1.0,1.0],line=0,color=cgColor('light gray'),/overplot
      
      w_eta = wind_model.vs_M.(zInd).wind_eta
      cgPlot,wind_model.vs_M.(zInd).halo_mass,w_eta/(1+w_eta),$
        line=0,thick=!p.thick*4,color=cgColor(tColor),/overplot
      
      cgPlot,[px_min,px_min],[yrange3[0],yrange3[1]],line=2,/overplot,color=cgColor(px_cl3)
      cgPlot,[px_max,px_max],[yrange3[0],yrange3[1]],line=2,/overplot,color=cgColor(px_cl3)
      
      foreach run,runs,i do begin
        yy = mbv.(i).(zInd).mode0_all.eta_abcd_unbinned
        cgPlot,mbv.(i).(zInd).mode0_all.eta_w_xvals,yy/(1+yy),$
          psym=3,color=sP.(i).(0).(0).colors[cInd],/overplot
        
        yy = mbv.(i).(zInd).mode0_all.eta_abcd_binned[*,1]
        cgPlot,sMasses,smooth(yy/(1+yy),sK,/nan),$
          line=0,color=sP.(i).(0).(0).colors[cInd],/overplot
          
        for j=0,2 do begin
          yy = mbv.(i).(zInd).mode0_all.eta_abcd_binned[*,j]
          cgPlot,sMasses,smooth(yy/(1+yy),sK,/nan),$
            line=2,color=sP.(i).(0).(0).colors[cInd],/overplot ; 25th/75th percentiles
        endfor  
      endforeach
           
      legend,['z='+string(redshifts[zInd],format='(f3.1)')],/top,/right
           
    end_PS
  
    ; plot (10c) - abcd2 method
    start_PS,sP.(0).(0).(0).plotPath + 'wind_scalings_abcd2_' + plotStr + '_zInd'+str(zInd)+'.eps', $
      xs=10.0, ys=4.0
    
      ; plot config
      xrange = xrange_halo + [0.0,1.0]
      tColor = [150,150,150]
      bColor = [150,50,50]
      yrange2 = [1,10000]
      yrange3 = [0.6,1.05]
      
      legendStrs   = [sP.(0).(0).(0).simName,$
                      sP.(1).(0).(0).simName,$
                      'Illustris \eta_w',$
                      'Illustris v_w [km/s]']
      legendColors = [sP.(0).(0).(0).colors[cInd],sP.(1).(0).(0).colors[cInd],$
                      getColor24(tColor),getColor24(bColor)]

      plot_pos = plot_pos(rows=1,cols=2,/gap)
          
      ; plot (a) - eta_w
      cgPlot,[0],[0],/nodata,xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),$
        ytitle=textoidl("\eta_w or v_w"),$
        xrange=xrange,/xs,/ylog,yminor=0,yrange=yrange2,pos=plot_pos[0],/noerase
        
      cgPlot,wind_model.vs_M.(zInd).halo_mass,wind_model.vs_M.(zInd).wind_eta,$
        line=0,thick=!p.thick*4,color=cgColor(tColor),/overplot
        
      cgPlot,[px_min,px_min],[yrange2[0],yrange2[1]],line=2,/overplot,color=cgColor(px_cl3)
      cgPlot,[px_max,px_max],[yrange2[0],yrange2[1]],line=2,/overplot,color=cgColor(px_cl3)
      
      cgPlot,wind_model.vs_M.(zInd).halo_mass,wind_model.vs_M.(zInd).wind_vel,$
        line=0,thick=!p.thick*2,color=cgColor(bColor),/overplot
      
      foreach run,runs,i do begin
        cgPlot,mbv.(i).(zInd).mode0_all.eta_w_xvals,mbv.(i).(zInd).mode0_all.eta_abcd2_unbinned,$
          psym=3,color=sP.(i).(0).(0).colors[cInd],/overplot
        
        w = where(mbv.(i).(zInd).mode0_all.eta_abcd2_binned eq 0.0,count)
        if count gt 0 then mbv.(i).(zInd).mode0_all.eta_abcd2_binned[w] = !values.f_nan
        
        cgPlot,sMasses,smooth(mbv.(i).(zInd).mode0_all.eta_abcd2_binned[*,1],sK,/nan),$
          line=0,color=sP.(i).(0).(0).colors[cInd],/overplot
          
        for j=0,2 do $
          cgPlot,sMasses,smooth(mbv.(i).(zInd).mode0_all.eta_abcd2_binned[*,j],sK,/nan),$
            line=2,color=sP.(i).(0).(0).colors[cInd],/overplot ; 25th/75th percentiles
      endforeach

      legend,textoidl(legendStrs),textcolors=legendColors,/top,/right  
      
      ; plot (b) - eta/(1+eta)
      cgPlot,[0],[0],/nodata,$
        xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),ytitle=textoidl("\eta_w / (1+\eta_w)"),$
        xrange=xrange,/xs,yrange=yrange3,/ys,pos=plot_pos[1],/noerase
        
      cgPlot,xrange,[1.0,1.0],line=0,color=cgColor('light gray'),/overplot
      
      w_eta = wind_model.vs_M.(zInd).wind_eta
      cgPlot,wind_model.vs_M.(zInd).halo_mass,w_eta/(1+w_eta),$
        line=0,thick=!p.thick*4,color=cgColor(tColor),/overplot
      
      cgPlot,[px_min,px_min],[yrange3[0],yrange3[1]],line=2,/overplot,color=cgColor(px_cl3)
      cgPlot,[px_max,px_max],[yrange3[0],yrange3[1]],line=2,/overplot,color=cgColor(px_cl3)
      
      foreach run,runs,i do begin
        yy = mbv.(i).(zInd).mode0_all.eta_abcd2_unbinned
        cgPlot,mbv.(i).(zInd).mode0_all.eta_w_xvals,yy/(1+yy),$
          psym=3,color=sP.(i).(0).(0).colors[cInd],/overplot
        
        yy = mbv.(i).(zInd).mode0_all.eta_abcd2_binned[*,1]
        cgPlot,sMasses,smooth(yy/(1+yy),sK,/nan),$
          line=0,color=sP.(i).(0).(0).colors[cInd],/overplot
          
        for j=0,2 do begin
          yy = mbv.(i).(zInd).mode0_all.eta_abcd2_binned[*,j]
          cgPlot,sMasses,smooth(yy/(1+yy),sK,/nan),$
            line=2,color=sP.(i).(0).(0).colors[cInd],/overplot ; 25th/75th percentiles
        endfor  
      endforeach
           
      legend,['z='+string(redshifts[zInd],format='(f3.1)')],/top,/right
           
    end_PS
  
  endforeach ;zInds
  
  ; plot (11)
  start_PS,sP.(0).(0).(0).plotPath + 'wind_scalings_z_' + plotStrMB + '.eps', xs=10.0, ys=4.0
   
    plot_pos = plot_pos(row=1,col=2,/gap)
    
    if massInd ge 20 then yrange2=[1,100]
    if massInd lt 20 then yrange2=[1,1000]
    yrange3 = [0.5,1.05]
        
    ; plot (a)
    cgPlot,[0],[0],/nodata,$
      xtitle="Redshift",ytitle=textoidl("\eta_w"),$
      xrange=xrange_z,/xs,yrange=yrange2,/ys,/ylog,yminor=0,pos=plot_pos[0]
      
    cgPlot,wind_model.vs_z.z_vals,wind_model.vs_z.wind_eta,color=cgColor(tColor),thick=!p.thick*4,/overplot
            
    foreach run,runs,i do begin
      ; collect in redshift
      zz_eta_w = fltarr(n_elements(redshifts),3)
      foreach redshift,redshifts,j do zz_eta_w[j,*] = mbv.(i).(j).mode0_all.eta_w_binned[massInd,*]
      
      cgPlot,redshifts,zz_eta_w[*,1],psym=-16,color=sP.(i).(0).(0).colors[cInd],/overplot
      cgPlot,redshifts,zz_eta_w[*,0],line=2,color=sP.(i).(0).(0).colors[cInd],/overplot ;25th
      cgPlot,redshifts,zz_eta_w[*,2],line=2,color=sP.(i).(0).(0).colors[cInd],/overplot ;75th
      
      ; un-binned
      foreach redshift,redshifts,j do begin
        w = where(mbv.(i).(j).mode0_all.eta_w_xvals ge px_min and $
                  mbv.(i).(j).mode0_all.eta_w_xvals lt px_max,count)
        if count eq 0 then continue
        
        zz_eta_w_unbinned = mbv.(i).(j).mode0_all.eta_w_unbinned[w]
        cgPlot,replicate(redshift,count),zz_eta_w_unbinned,$
          psym=3,color=sP.(i).(0).(0).colors[cInd],/overplot
      endforeach
    endforeach
    
    legend,textoidl([string(px_min,format='(f4.1)')+' < M_{halo} < '+string(px_max,format='(f4.1)')]),$
      /top,/right
    
    ; plot (b)
    cgPlot,[0],[0],/nodata,$
      xtitle="Redshift",ytitle=textoidl("\eta_w / (1+\eta_w)"),$
      xrange=xrange_z,/xs,yrange=yrange3,/ys,pos=plot_pos[1],/noerase
      
    cgPlot,[0.0,5.0],[1.0,1.0],line=0,color=cgColor('light gray'),/overplot
    cgPlot,wind_model.vs_z.z_vals,wind_model.vs_z.wind_eta/(1+wind_model.vs_z.wind_eta),$
      color=cgColor(tColor),thick=!p.thick*4,/overplot
    
    foreach run,runs,i do begin
      ; collect in redshift
      zz_eta_w = fltarr(n_elements(redshifts),3)
      foreach redshift,redshifts,j do zz_eta_w[j,*] = mbv.(i).(j).mode0_all.eta_w_binned[massInd,*]
      
      cgPlot,redshifts,zz_eta_w[*,1]/(1+zz_eta_w[*,1]),$
        line=0,psym=-16,color=sP.(i).(0).(0).colors[cInd],/overplot ; median
      cgPlot,redshifts,zz_eta_w[*,0]/(1+zz_eta_w[*,0]),$
        line=2,color=sP.(i).(0).(0).colors[cInd],/overplot ;25th
      cgPlot,redshifts,zz_eta_w[*,2]/(1+zz_eta_w[*,2]),$
        line=2,color=sP.(i).(0).(0).colors[cInd],/overplot ;75th
      
      ; un-binned
      foreach redshift,redshifts,j do begin
        w = where(mbv.(i).(j).mode0_all.eta_w_xvals ge px_min and $
                  mbv.(i).(j).mode0_all.eta_w_xvals lt px_max,count)
        if count eq 0 then continue
        
        zz_eta_w_unbinned = mbv.(i).(j).mode0_all.eta_w_unbinned[w]
        cgPlot,replicate(redshift,count),zz_eta_w_unbinned/(1+zz_eta_w_unbinned),$
          psym=3,color=sP.(i).(0).(0).colors[cInd],/overplot
      endforeach
    endforeach
      
    legend,textoidl(legendStrs[0:-2]),textcolors=legendColors[0:-2],/bottom,/left
  
  end_PS
  
  ; plot (12)
  start_PS,sP.(0).(0).(0).plotPath + 'wind_scalings_z_' + plotStrMB + '_zInd'+str(zInd) + '.eps', xs=10.0, ys=4.0
   
    zInd = where(redshifts eq 0.0)
    plot_pos = plot_pos(row=1,col=2,/gap)
    
    ; plot (a): left of plot10
    sMasses = codeMassToLogMsun( logMsunToCodeMass( mbv.(0).(0).(0).logMassBinCen ) / units.hubbleParam )
        
    xrange = xrange_halo + [0.0,1.0]
    tColor = [150,150,150]
    bColor = [150,50,50]
    yrange2 = [1,10000]
    yrange3 = [0.6,1.05]
    
    legendStrs   = [sP.(0).(0).(0).simName,$
                    sP.(1).(0).(0).simName,$
                    'Illustris \eta_w',$
                    'Illustris v_w [km/s]']
    legendColors = [sP.(0).(0).(0).colors[cInd],sP.(1).(0).(0).colors[cInd],$
                    getColor24(tColor),getColor24(bColor)]

    plot_pos = plot_pos(rows=1,cols=2,/gap)
        
    cgPlot,[0],[0],/nodata,xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),$
      ytitle=textoidl("\eta_w"),$
      xrange=xrange,/xs,/ylog,yminor=0,yrange=yrange2,ys=9,pos=plot_pos[0],/noerase
      
    cgPlot,wind_model.vs_M.(zInd).halo_mass,wind_model.vs_M.(zInd).wind_eta,$
      line=0,thick=!p.thick*4,color=cgColor(tColor),/overplot
      
    cgPlot,[px_min,px_min],[yrange2[0],yrange2[1]],line=2,/overplot,color=cgColor(px_cl3)
    cgPlot,[px_max,px_max],[yrange2[0],yrange2[1]],line=2,/overplot,color=cgColor(px_cl3)
       
    foreach run,runs,i do begin
      cgPlot,mbv.(i).(zInd).mode0_all.eta_w_xvals,mbv.(i).(zInd).mode0_all.eta_w_unbinned,$
        psym=3,color=sP.(i).(0).(0).colors[cInd],/overplot
      
      w = where(mbv.(i).(zInd).mode0_all.eta_w_binned eq 0.0,count)
      if count gt 0 then mbv.(i).(zInd).mode0_all.eta_w_binned[w] = !values.f_nan
      
      cgPlot,sMasses,smooth(mbv.(i).(zInd).mode0_all.eta_w_binned[*,1],sK,/nan),$
        line=0,color=sP.(i).(0).(0).colors[cInd],/overplot
        
      for j=0,2 do $
        cgPlot,sMasses,smooth(mbv.(i).(zInd).mode0_all.eta_w_binned[*,j],sK,/nan),$
          line=2,color=sP.(i).(0).(0).colors[cInd],/overplot ; 25th/75th percentiles
    endforeach
    
    axis,yaxis=1,/ylog,yrange=[10,1000],ytitle=textoidl("v_w [km/s]"),yminor=0,/save;,$
      ;yticks=2,ytickv=[10,100,1000],ytickn=textoidl(['10^1','10^2','10^3'])
    cgPlot,wind_model.vs_M.(zInd).halo_mass,wind_model.vs_M.(zInd).wind_vel,$
      line=0,thick=!p.thick*2,color=cgColor(bColor),/overplot

    legend,textoidl(legendStrs[2:3]),textcolors=legendColors[2:3],/top,/left   
    legend,['z='+string(redshifts[zInd],format='(f3.1)')],/top,/right
    
    ; plot (b): left of plot11
    if massInd ge 20 then yrange2=[1,100]
    if massInd lt 20 then yrange2=[1,1000]
    yrange3 = [0.5,1.05]
    
    cgPlot,[0],[0],/nodata,xtitle="Redshift",ytitle=textoidl("\eta_w"),$
      xrange=xrange_z,/xs,yrange=yrange2,/ys,/ylog,yminor=0,pos=plot_pos[1]+[0.04,0,0.04,0],/noerase
      
    cgPlot,wind_model.vs_z.z_vals,wind_model.vs_z.wind_eta,color=cgColor(tColor),thick=!p.thick*4,/overplot
            
    foreach run,runs,i do begin
      ; collect in redshift
      zz_eta_w = fltarr(n_elements(redshifts),3)
      foreach redshift,redshifts,j do zz_eta_w[j,*] = mbv.(i).(j).mode0_all.eta_w_binned[massInd,*]
      
      cgPlot,redshifts,zz_eta_w[*,1],psym=-16,color=sP.(i).(0).(0).colors[cInd],/overplot
      cgPlot,redshifts,zz_eta_w[*,0],line=2,color=sP.(i).(0).(0).colors[cInd],/overplot ;25th
      cgPlot,redshifts,zz_eta_w[*,2],line=2,color=sP.(i).(0).(0).colors[cInd],/overplot ;75th
      
      ; un-binned
      foreach redshift,redshifts,j do begin
        w = where(mbv.(i).(j).mode0_all.eta_w_xvals ge px_min and $
                  mbv.(i).(j).mode0_all.eta_w_xvals lt px_max,count)
        if count eq 0 then continue
        
        zz_eta_w_unbinned = mbv.(i).(j).mode0_all.eta_w_unbinned[w]
        cgPlot,replicate(redshift,count),zz_eta_w_unbinned,$
          psym=3,color=sP.(i).(0).(0).colors[cInd],/overplot
      endforeach
    endforeach
    
    legend,textoidl([string(px_min,format='(f4.1)')+' < M_{halo} < '+string(px_max,format='(f4.1)')]),$
      /top,/left
    legend,textoidl(legendStrs[0:1]),textcolors=legendColors[0:1],/top,/right 
      
  end_PS
  
  ; plot (13) - correction factor given the SFR/NetAcc ratio
  z_targets = [0.0,1.0]
  
  foreach z_target,z_targets do begin
    start_PS,sP.(0).(0).(0).plotPath + 'wind_correction_' + plotStr + '_z'+str(round(z_target*100))+'.eps', xs=10.5, ys=3.5
   
      zInd = closest(redshifts,z_target)
   
      ; plot config
      xrange  = xrange_halo + [0.0,1.0]
      yrange2 = [-4.0,4.0]
      
      plot_pos = plot_pos(rows=1,cols=2,/gap)
          
      foreach run,runs,i do begin
        cgPlot,[0],[0],/nodata,xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),$
          ytitle=textoidl("log( SFR / GalaxyNetAccRate )"),$
          xrange=xrange,/xs,yrange=yrange2,pos=plot_pos[i],/noerase
          
        cgPlot,xrange,[0.0,0.0],line=0,color=cgColor('gray'),/overplot
          
        cgPlot,[px_min,px_min],[yrange2[0],yrange2[1]],line=2,/overplot,color=cgColor(px_cl3)
        cgPlot,[px_max,px_max],[yrange2[0],yrange2[1]],line=2,/overplot,color=cgColor(px_cl3)
      
        cgPlot,mbv.(i).(zInd).mode0_all.eta_w_xvals,mbv.(i).(zInd).mode0_all.eta_sfr2_net,$
          psym=4,color=cgColor('black'),/overplot
        
        fit = mbv.(i).(zInd).mode0_all.eta_fit
        cgPlot,xrange,xrange*fit[1]+fit[0],line=0,color=sP.(i).(zInd).mode0_all.colors[cInd],/overplot
        
        legend,[sP.(i).(0).(0).simName,'z='+string(z_target,format='(f4.1)')],$
          textcolor=[sP.(i).(0).(0).colors[cInd],0],/top,/right
      endforeach
  
    end_PS
    
    ; plot (13b) - abcd method
    start_PS,sP.(0).(0).(0).plotPath + 'wind_correction_abcd_' + plotStr + '_z'+str(round(z_target*100))+'.eps', xs=10.5, ys=3.5
   
      zInd = closest(redshifts,z_target)
   
      ; plot config
      xrange  = xrange_halo + [0.0,1.0]
      yrange2 = [-4.0,4.0]
      
      plot_pos = plot_pos(rows=1,cols=2,/gap)
          
      foreach run,runs,i do begin
        cgPlot,[0],[0],/nodata,xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),$
          ytitle=textoidl("log( SFR / GalaxyNetAccRateAbcd )"),$
          xrange=xrange,/xs,yrange=yrange2,pos=plot_pos[i],/noerase
          
        cgPlot,xrange,[0.0,0.0],line=0,color=cgColor('gray'),/overplot
          
        cgPlot,[px_min,px_min],[yrange2[0],yrange2[1]],line=2,/overplot,color=cgColor(px_cl3)
        cgPlot,[px_max,px_max],[yrange2[0],yrange2[1]],line=2,/overplot,color=cgColor(px_cl3)
      
        cgPlot,mbv.(i).(zInd).mode0_all.eta_w_xvals,mbv.(i).(zInd).mode0_all.eta_sfr2_net_abcd,$
          psym=4,color=cgColor('black'),/overplot
        
        fit = mbv.(i).(zInd).mode0_all.eta_fit
        cgPlot,xrange,xrange*fit[1]+fit[0],line=0,color=sP.(i).(zInd).mode0_all.colors[cInd],/overplot
        
        legend,[sP.(i).(0).(0).simName,'z='+string(z_target,format='(f4.1)')],$
          textcolor=[sP.(i).(0).(0).colors[cInd],0],/top,/right
      endforeach
  
    end_PS
    
    ; plot (13c) - abcd2 method
    start_PS,sP.(0).(0).(0).plotPath + 'wind_correction_abcd2_' + plotStr + '_z'+str(round(z_target*100))+'.eps', xs=10.5, ys=3.5
   
      zInd = closest(redshifts,z_target)
   
      ; plot config
      xrange  = xrange_halo + [0.0,1.0]
      yrange2 = [-4.0,4.0]
      
      plot_pos = plot_pos(rows=1,cols=2,/gap)
          
      foreach run,runs,i do begin
        cgPlot,[0],[0],/nodata,xtitle=textoidl("M_{halo} [_{ }log M_{sun }]"),$
          ytitle=textoidl("log( SFR / GalaxyNetAccRateAbcd2 )"),$
          xrange=xrange,/xs,yrange=yrange2,pos=plot_pos[i],/noerase
          
        cgPlot,xrange,[0.0,0.0],line=0,color=cgColor('gray'),/overplot
          
        cgPlot,[px_min,px_min],[yrange2[0],yrange2[1]],line=2,/overplot,color=cgColor(px_cl3)
        cgPlot,[px_max,px_max],[yrange2[0],yrange2[1]],line=2,/overplot,color=cgColor(px_cl3)
      
        cgPlot,mbv.(i).(zInd).mode0_all.eta_w_xvals,mbv.(i).(zInd).mode0_all.eta_sfr2_net_abcd2,$
          psym=4,color=cgColor('black'),/overplot
        
        fit = mbv.(i).(zInd).mode0_all.eta_fit
        cgPlot,xrange,xrange*fit[1]+fit[0],line=0,color=sP.(i).(zInd).mode0_all.colors[cInd],/overplot
        
        legend,[sP.(i).(0).(0).simName,'z='+string(z_target,format='(f4.1)')],$
          textcolor=[sP.(i).(0).(0).colors[cInd],0],/top,/right
      endforeach
  
    end_PS
    
  endforeach ; z_target
  endif ;acr=0 or plotUnused
  
  stop       
end

; plotMaxValsVsRedshift(): TODO

pro plotMaxValsVsRedshift
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
  message,'revamp'
  
  ; config
  runs       = ['feedback','tracer'] ;,'gadget']
  redshifts  = [3.0,2.0,1.0,0.0]
  res        = 512
  accMode    = 'smooth' ; accretion mode: all, smooth, clumpy, stripped, recycled
  timeWindow = 500.0    ; consider accretion over this past time range (Myr)
                        ; 250.0 500.0 1000.0 "all" "tVir_tIGM" or "tVir_tIGM_bin"
  tVirInd    = 0        ; use Tmax/Tviracc=1 to separate hot vs. cold
  
  ; max histograms config
  massBinInd = 4        ; plot 11.0<logM<11.5 halos for ValMax histos
  tagNames   = ['Tmax','TmaxTviracc','EntMax','DensMax','MachMax'] ; plot each quantity
  pParams = { TmaxTviracc : {xrange:[-2.2,1.4], xtickv:[-2,-1,0,1],   xlabel : "log ( T_{max} / T_{vir,acc} )"} ,$
              Tmax        : {xrange:[3.5,7.5],  xtickv:[4,5,6,7],     xlabel : "log ( T_{max} )"}               ,$
              EntMax      : {xrange:[4.5,10.0], xtickv:[5,6,7,8,9],   xlabel : "log ( S_{max} ) [K cm^{2 }]"}   ,$
              DensMax     : {xrange:[-10,0],    xtickv:[-8,-6,-4,-2], xlabel : "log ( \rho_{max} )"}            ,$
              MachMax     : {xrange:[0,100],    xtickv:[10,30,50,80], xlabel : "M_{max}"} }
  
  ; plot config
  lines   = [0,1]       ; gal/both,gmem
  sK      = 3           ; smoothing kernel size  
  cInd    = 1           ; color index

  xrange_halo = [9.0,12.0]
  yrange_hist = [6e-4,2e-1]
  
  ; load
  tempFilename = '/n/home07/dnelson/temp_plotRatesMaxValsAt4Redshifts.sav'
  if ~file_test(tempFilename) then begin
  
  foreach run,runs,i do begin
    sP_z = {} & mbv_z = {} & bth_z = {}
    
    ; make for all the redshifts
    foreach redshift,redshifts,j do begin
      sP_loc = simParams(res=res,run=run,redshift=redshift)
      mbv_loc = haloMassBinValues(sP=sP_loc,accMode=accMode,timeWindow=timeWindow)

      ; load group cat for subgroup masses
      gc = loadGroupCat(sP=sP_loc,/skipIDs)
      gcMasses = codeMassToLogMsun(gc.subgroupMass)
      
      ; calculate specific net rates, galaxy by galaxy (yr -> Gyr)
      maxHist = n_elements(mbv_loc.galaxy.netRate.total.hot.num)
      gcMasses = gcMasses[0:maxHist-1]
        
      gal_hot   = reform(mbv_loc.galaxy.netRate.total.hot.tVirAcc[tVirInd,*]) / 10.0^gcMasses * 1e9
      gal_cold  = reform(mbv_loc.galaxy.netRate.total.cold.tVirAcc[tVirInd,*]) / 10.0^gcMasses * 1e9
      halo_hot  = reform(mbv_loc.halo.netRate.total.hot.tVirAcc[tVirInd,*]) / 10.0^gcMasses * 1e9
      halo_cold = reform(mbv_loc.halo.netRate.total.cold.tVirAcc[tVirInd,*]) / 10.0^gcMasses * 1e9

      specific_net = { gal_hot   : fltarr(n_elements(mbv_loc.logMassBinCen)) + !values.f_nan ,$
                       gal_cold  : fltarr(n_elements(mbv_loc.logMassBinCen)) + !values.f_nan ,$
                       halo_hot  : fltarr(n_elements(mbv_loc.logMassBinCen)) + !values.f_nan ,$
                       halo_cold : fltarr(n_elements(mbv_loc.logMassBinCen)) + !values.f_nan  }
      coldfrac_net = { gal  : fltarr(n_elements(mbv_loc.logMassBinCen)) + !values.f_nan ,$
                       halo : fltarr(n_elements(mbv_loc.logMassBinCen)) + !values.f_nan  }
      
      for k=0,n_elements(mbv_loc.logMassBinCen)-1 do begin
        ; restrict medians to primary only (num>0)
        w = where(gcMasses gt mbv_loc.logMassBins[k] and gcMasses le mbv_loc.logMassBins[k+1] and $
                  (mbv_loc.galaxy.coldFrac.total.num+mbv_loc.halo.coldFrac.total.num) gt 0,count)
                  
        if count eq 0 then continue
        if max(w) ge maxHist then message,'Error'
        
        specific_net.gal_hot[k]   = median( gal_hot[w] )
        specific_net.gal_cold[k]  = median( gal_cold[w] )
        specific_net.halo_hot[k]  = median( halo_hot[w] )
        specific_net.halo_cold[k] = median( halo_cold[w] )
        
        coldfrac_net.gal[k]  = median( gal_cold[w] / (gal_hot[w]+gal_cold[w]) )
        coldfrac_net.halo[k] = median( halo_cold[w] / (halo_hot[w]+halo_cold[w]) )
      endfor ;k
      
      mbv_loc = { galaxyMedian  : mbv_loc.galaxyMedian  ,$
                  haloMedian    : mbv_loc.haloMedian    ,$
                  logMassBinCen : mbv_loc.logMassBinCen ,$
                  specific_net  : specific_net          ,$
                  coldfrac_net  : coldfrac_net           }
      
      ; add to keeper array for this redshift
      sP_z  = mod_struct(sP_z, 'redshift'+str(j), sP_loc)
      mbv_z = mod_struct(mbv_z, 'redshift'+str(j), mbv_loc)
      ;bth_z = mod_struct(bth_z, 'redshift'+str(j), $
      ;  binValMaxHistos(sP=sP_z.(j),accMode=accMode,timeWindow=timeWindow))
    endforeach
    
    ; put this mode collection into mbv, once per run, and sP collection also
    mbv = mod_struct(mbv, 'mbv'+str(i), mbv_z)
    ;bth = mod_struct(bth, 'bth'+str(i), bth_z)
    sP  = mod_struct(sP, 'sP'+str(i), sP_z)
  endforeach
  
    save,mbv,sP,filename=tempFilename
    print,'Saved TEMP file: ['+tempFilename+']'
  endif else begin
    restore,tempFilename,/verbose
  endelse
  
  ; strings
  if isnumeric(timeWindow) then twStr = str(fix(timeWindow))
  if ~isnumeric(timeWindow) then twStr = "-"+timeWindow
  
  plotStr   = ''
  simNames  = []
  simColors = []
  
  foreach run,runs,i do begin
    plotStr   = plotStr + sP.(i).(0).plotPrefix + '.'
    simNames  = [simNames, sP.(i).(0).simName]
    simColors = [simColors, sP.(i).(0).colors[cInd]]
  endforeach

  plotStr += str(res) + '_' + str(sP.(0).(0).snap) + '_tw' + twStr + '_am-' + accMode
  
  ; maximum val (temp,temp/tvir,ent,dens,mach) histos (allgal,gmem) (one halo mass bin)
  foreach tagName,tagNames do begin
    bVind   = ( where(tag_names(bth.(0).(0).allGal) eq strupcase(tagName)) )[0]
    pPind   = ( where(tag_names(pParams) eq strupcase(tagName)) )[0]
    
    xrange = pParams.(pPind).xrange
    xtickv = pParams.(pPind).xtickv
    
    start_PS, sP.(0).(0).plotPath + 'maxHistosRedshift_' + tagName + '.massBin' + $
              str(massBinInd) + '.' + plotStr + '.eps', /big
    
    for zind=0,3 do begin
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange_hist,/xs,/ys,/ylog,yminor=0,$
        ytitle=textoidl("Gas Mass Fraction"),$
        xtitle=textoidl( pParams.(pPind).xlabel ),$
        /noerase,pos=pos[zind],xticks=n_elements(xtickv)-1,xtickv=xtickv
    
      for i=0,n_tags(sP)-2 do begin
      
        ; gal
        xvals = bth.(i).(zind).params.binLoc.(bVind)
        yvals = float( bth.(i).(zind).allGal.(bVind)[massBinInd,*] ) / $
                total( bth.(i).(zind).allGal.(bVind)[j,*] )
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[0],/overplot ; allgal
        
        ; gmem
        yvals = float( bth.(i).(zind).gmem.(bVind)[massBinInd,*] ) / $
                total( bth.(i).(zind).gmem.(bVind)[j,*] )
        cgPlot,xvals,yvals,color=sP.(i).(zind).colors[cInd],line=lines[1],/overplot ; gmem
      endfor
      
      ; redshift legend
      legend,textoidl('z_{ }=_{ }'+string(redshifts[zind],format='(f3.1)')),/bottom,/left,$
        charsize=!p.charsize-0.2
      
      cgPlot,[0,0],[9e-4,0.13],line=2,color=cgColor('black'),/overplot
      
      ; panel-specific legends
      if zind eq 2 then $
        legend,['galaxy','halo'],linestyle=lines,color=cgColor('dark gray'),$
        textcolors=replicate('dark gray',2),box=0,/top,/left,charsize=!p.charsize-0.2
        
      if zind eq 3 then $
        legend,simNames,textcolors=simColors,box=0,/top,/left,charsize=!p.charsize-0.2
    endfor

    end_PS
  endforeach ;tagName
   
  stop
end

; plotEvolIGMTemp(): plot evolving temperature of the IGM gas (outside all FOFs) vs redshift

pro plotEvolIGMTemp

  ; config
  sP_FB   = simParams(res=512,run='feedback')
  sP_noFB = simParams(res=512,run='tracer')
  
  igmEvol_FB   = evolIGMTemp(sP=sP_FB)
  igmEvol_noFB = evolIGMTemp(sP=sP_noFB)
  
  lines    = [0,2,1,3]    ; mean,median,min,max
  cInd     = 1            ; color index
  xrange_z = [-0.15,6.15] ; redshift
  yrange_T = [1.5,9.8]    ; temp (log K)
  yrange_P = [1e-7,1e-1]  ; pdf
  zDists   = [5.0,3.0,2.0,1.0,0.5,0.0] ; 6 redshifts to plot full distributions
  
  ; plot (1) - mean/median evolution of IGM gas temperature with redshift
  pos_single = [0.14,0.11,0.94,0.87] ; single plot with second x-axis above
  
  plotStr = sP_FB.savPrefix + str(sP_FB.res) + "_" + sP_noFB.savPrefix + str(sP_noFB.res)
  
  start_PS, sP_FB.plotPath + 'igmTemp.comp.'+plotStr+'.eps'
    
    cgPlot,[0],[0],/nodata,xrange=xrange_z,yrange=yrange_T,xs=9,/ys,$
      ytitle="IGM Gas Temperature [ log K ]",xtitle="Redshift",pos=pos_single
    
    universeage_axis, xrange_z, yrange_T
    
    ; with feedback
    cgPlot,igmEvol_FB.redshifts,igmEvol_FB.meanTemp,$
      color=sP_FB.colors[cInd],line=lines[0],/overplot
    cgPlot,igmEvol_FB.redshifts,igmEvol_FB.medianTemp,$
      color=sP_FB.colors[cInd],line=lines[1],/overplot
    cgPlot,igmEvol_FB.redshifts,igmEvol_FB.minTemp,$
      color=sP_FB.colors[cInd],line=lines[2],/overplot
    cgPlot,igmEvol_FB.redshifts,igmEvol_FB.maxTemp,$
      color=sP_FB.colors[cInd],line=lines[3],/overplot

    ; without feedback
    cgPlot,igmEvol_noFB.redshifts,igmEvol_noFB.meanTemp,$
      color=sP_noFB.colors[cInd],line=lines[0],/overplot
    cgPlot,igmEvol_noFB.redshifts,igmEvol_noFB.medianTemp,$
      color=sP_noFB.colors[cInd],line=lines[1],/overplot
    cgPlot,igmEvol_noFB.redshifts,igmEvol_noFB.minTemp,$
      color=sP_noFB.colors[cInd],line=lines[2],/overplot
    cgPlot,igmEvol_noFB.redshifts,igmEvol_noFB.maxTemp,$
      color=sP_noFB.colors[cInd],line=lines[3],/overplot
    
    ; legends
    legend,[sP_FB.simName,sP_noFB.simName],textcolors=[sP_FB.colors[1],sP_noFB.colors[1]],$
      pos=[xrange_z[0]+2.6,yrange_T[1]-0.2]
      
    legend,['mean','median','min','max'],linestyle=lines,/top,/right
    
  end_PS
  
  ; plot (2) - distributions
  start_PS, sP_FB.plotPath + 'igmTemp.compdist.'+plotStr+'.eps', xs=10*1.4, ys=6*1.4
  
    pos = plot_pos(total=n_elements(zDists),/gap)
    
    for i=0,n_elements(zDists)-1 do begin
     
      ; locate index
      curInd_FB = closest( igmEvol_FB.redshifts, zDists[i] )
      curInd_noFB = closest( igmEvol_FB.redshifts, zDists[i] )
    
      ; plot
      cgPlot,[0],[0],/nodata,xrange=yrange_T,yrange=yrange_P,xs=1,ys=1,/ylog,yminor=0,$
        ytitle="PDF",xtitle="IGM Gas Temperature [ log K ]",pos=pos[i],/noerase
  
      yy = igmEvol_FB.tempDist[curInd_FB,*]
      cgPlot,igmEvol_FB.tempBinCen,yy/total(yy),color=sP_FB.colors[cInd],/overplot
      
      yy = igmEvol_noFB.tempDist[curInd_noFB,*]
      cgPlot,igmEvol_noFB.tempBinCen,yy/total(yy),color=sP_noFB.colors[cInd],/overplot
      
      legend,["z="+string(zDists[i],format='(f3.1)')],/top,/left
  
      if i eq 2 then $
        legend,[sP_FB.simName,sP_noFB.simName],textcolors=[sP_FB.colors[1],sP_noFB.colors[1]],/top,/right
    endfor
  
  end_PS
  
  stop

end
