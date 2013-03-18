; haloModelsPlot.pro
; plotting: theoretical models for DM halos (e.g. NFW, SIS) including their gas
;           comparison with SAM models
; dnelson mar.2013

; crotonTest

pro crotonTest
 
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()
 
  ; config
  zBounds = [0.0,7.0]
  zSteps  = 100
  
  M0 = 100.0 ; 10^12 Msun
  
  redshifts = linspace(zBounds[1],zBounds[0],zSteps)
  times = redshiftToAgeFlat(redshifts) * 1e9 ; yr
  
  ; load Sutherland & Dopita (1993) cooling tables
  tables = interpLambdaSD93()
  
  ; arrays: halo properties
  M_DM    = fltarr(zSteps) ; Msun, total halo mass
  M_b     = fltarr(zSteps) ; Msun, total baryonic halo mass
  R_vir   = fltarr(zSteps) ; km
  V_c     = fltarr(zSteps) ; km/s
  T_vir   = fltarr(zSteps) ; K
  dM_b    = fltarr(zSteps) ; Msun/yr
  
  ; Croton model
  M_cold     = fltarr(zSteps) ; Msun, e.g. the only destination for hot gas with no SF/ejecta
  M_hot      = fltarr(zSteps) ; Msun, e.g. Mgas0 (hot gas available to cool)
  dM_cool    = fltarr(zSteps) ; Msun/yr
  r_cool     = fltarr(zSteps) ; km
  rapid_mode = bytarr(zSteps) ; 0=no, 1=yes
  
  ; Kang model
  M_cold_kang     = fltarr(zSteps) ; Msun, e.g. the only destination for hot gas with no SF/ejecta
  M_hot_kang      = fltarr(zSteps) ; Msun, e.g. Mgas0 (hot gas available to cool)
  dM_cool_kang    = fltarr(zSteps) ; Msun/yr
  r_cool_kang     = fltarr(zSteps) ; km
  rapid_mode_kang = bytarr(zSteps) ; 0=no, 1=yes
  
  ; timestep zero
  M_b[0] = units.f_b * haloMAH(redshifts[0],M0=M0,z0=zBounds[0])
  
  ; evolve in time
  ; --------------
  for i=1,zSteps-1 do begin
    ; halo properties
    delta_t = ( times[i] - times[i-1] )
    
    ; new halo baryonic mass
    M_DM[i]  = haloMAH(redshifts[i],M0=M0,z0=zBounds[0])
    M_b[i]   = M_DM[i] * units.f_b
    
    ; rate of growth of halo baryons
    dM_b[i] = (M_b[i] - M_b[i-1]) / delta_t ; 1e10 Msun/yr
    
    ; SIS DM halo
    sis_dm   = sis_profile(1.0, mass=M_DM[i], redshift=redshifts[i])
    R_vir[i] = sis_dm.r200
    
    ; Croton model
    ; ------------
      M_hot[i] = M_DM[i] * units.f_b - M_cold[i-1]
      
      if M_hot[i] gt 0.0 then begin
      
        sis_gas = sis_gas_profile(mass_hot=M_hot[i], sis_dm=sis_dm, tables=tables)
        r_cool[i] = sis_gas.r_cool_croton
    
        ; compare rcool to rvir, "rapid mode" if rcool>rvir
        if sis_gas.r_cool_croton ge sis_dm.r200 then rapid_mode[i] = 1B
      
        ; set cold gas accretion rate based on mode (1e10 Msun/yr)
        if rapid_mode[i] then begin
          ; rapid mode
          dM_cool[i] = M_hot[i] / delta_t ;(sis_gas.dynTime*1e9)
        endif else begin
          ; slow mode
          dM_cool[i] = 0.5 * M_hot[i] * sis_gas.r_cool_croton / (sis_dm.r200 * sis_gas.dynTime_halo*1e9)
        endelse
      
        if dM_cool[i] * delta_t gt M_hot[i] then begin
          print,'dM_cool * dt > M_hot',i
          dM_cool[i] = M_hot[i] / delta_t
        endif
      
        ; update cold gas
        M_cold[i] = M_cold[i-1] + dM_cool[i] * delta_t
      
      endif else begin
        print,'warning M_hot zero',i
      endelse
      
    ; Kang model
    ; ----------
      M_hot_kang[i] = M_DM[i] * units.f_b - M_cold_kang[i-1]
      
      if M_hot_kang[i] le 0.0 then message,'error'
      
      sis_gas = sis_gas_profile(mass_hot=M_hot_kang[i], sis_dm=sis_dm, tables=tables)
      r_cool_kang[i] = sis_gas.r_cool_hubble
      
      ; compare rcool_h to rvir, "rapid mode" if rcool_h>rvir
      if sis_gas.r_cool_hubble ge sis_dm.r200 then rapid_mode_kang[i] = 1B
      
      ; set cold gas accretion rate based on mode
      kang_timenorm = times[i] ;sis_gas.hubbleTime*1e9 ;
      
      if rapid_mode_kang[i] then begin
        dM_cool_kang[i] = M_hot_kang[i] / kang_timenorm
      endif else begin
        dM_cool_kang[i] = 0.5 * M_hot_kang[i] * sis_gas.r_cool_hubble / (sis_dm.r200 * kang_timenorm)
      endelse
      
      if dM_cool_kang[i] * delta_t gt M_hot_kang[i] then message,'error'
      
      ; update cold gas
      M_cold_kang[i] = M_cold_kang[i-1] + dM_cool_kang[i] * delta_t
      
  endfor
  
  ; plot
  xrange = [1,7]
  
  start_PS,'halo_SAM_vs_redshift.eps', ys=8, xs=6
    ; halo baryons
    cgPlot,1+redshifts,dM_b*1e10,xtitle='',ytitle=textoidl('dM/dt [M_{sun}/year]'),$
      /xlog,xrange=xrange,/xs,xminor=0,/ylog,yrange=[0.06,60.0],/ys,yminor=0,$
      xtickv=[1,2,3,4,5,6,7],xtickname=replicate(' ',10),xticks=6,position=(sP.pos_3x1)[0]
      
    ; croton model, dM/dt cold split into rapid and slow
    cgPlot,1+redshifts,dM_cool*1e10,/overplot,color=cgColor('red')
    
    w_rapid = where(rapid_mode eq 1B,count)
    cgPlot,1+redshifts[w_rapid],dM_cool[w_rapid]*1e10,psym=4,/overplot,color=cgColor('red')
    
    ; kang model
    cgPlot,1+redshifts,dM_cool_kang*1e10,/overplot,color=cgColor('blue')
    
    w_rapid = where(rapid_mode_kang eq 1B,count)
    cgPlot,1+redshifts[w_rapid],dM_cool_kang[w_rapid]*1e10,psym=4,/overplot,color=cgColor('blue')

    ; halo baryons
    cgPlot,1+redshifts,M_b*1e10,xtitle='',ytitle=textoidl('M_{cold} [M_{sun}]'),$
      /xlog,xrange=[1,7],/xs,xminor=0,/ylog,yrange=[1.5e8,2e11],/ys,yminor=0,$
      xtickv=[1,2,3,4,5,6,7],xtickname=replicate(' ',10),xticks=6,/noerase,position=(sP.pos_3x1)[1]
      
    ; croton model, dM/dt cold split into rapid and slow
    cgPlot,1+redshifts,M_cold*1e10,/overplot,color=cgColor('red')
    
    w_rapid = where(rapid_mode eq 1B,count)
    cgPlot,1+redshifts[w_rapid],M_cold[w_rapid]*1e10,psym=4,/overplot,color=cgColor('red')
    
    ; kang model
    cgPlot,1+redshifts,M_cold_kang*1e10,/overplot,color=cgColor('blue')
    
    w_rapid = where(rapid_mode_kang eq 1B,count)
    cgPlot,1+redshifts[w_rapid],M_cold_kang[w_rapid]*1e10,psym=4,/overplot,color=cgColor('blue')
    
    legend,['croton','kang'],textcolors=['red','blue'],box=0,/top,/right

    cgPlot,[0],[0],/nodata,xtitle='1+z',ytitle='Radius [kpc]',$
      /xlog,xrange=[1,7],/xs,xminor=0,yrange=[0,150],/ys,$
      xtickv=[1,2,3,4,5,6,7],xtickname=['1','2','3','4','5','6','7'],xticks=6,/noerase,position=(sP.pos_3x1)[2]
      
    cgPlot,1+redshifts,R_vir,/overplot,color=cgColor('black')
    cgPlot,1+redshifts,R_cool,/overplot,color=cgColor('red')
    cgPlot,1+redshifts,r_cool_kang,/overplot,color=cgColor('blue')
    legend,textoidl(['R_{vir}','R_{cool,croton}','R_{cool,kang}']),textcolors=['black','red','blue'],box=0,$
      position=[2.5,140]
  end_PS
  
end
  
; compVsMass(): at one redshift, compute models as a function of halo mass
; calculate critical halo mass crossover for rcool<>rvir
  
pro compVsMass
  
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()  
  
  ; config
  redshift = 2.0
  M0_range = [10.0,13.0] ; log(msun) at z=0
  
  hRes = 200
  hMasses = 10.0^(linspace(M0_range[0],M0_range[1],hRes))/1e10
  
  rPts = linspace(0.01,10.0,5000)
  
  ; load Sutherland & Dopita (1993) cooling tables
  tables = interpLambdaSD93()
  
  ; arrays
  hMass_z = fltarr(hRes)
  virRad  = fltarr(hRes)
  virTemp = fltarr(hRes)
  
  sis_coolRad = fltarr(hRes)
  nfw_coolRad = fltarr(hRes)
  nfw_coolRad_iso = fltarr(hRes)
  
  for i=0,hRes-1 do begin
    ; calculate halo masses / hot gas mass, at target redshift
    hMass_z[i] = haloMAH(redshift,M0=hMasses[i],z0=0.0)
    m_gas0 = hMass_z[i] * units.f_b
    print,codeMassToLogMsun(hMass_z[i])
    
    ; SIS DM halo and gas profiles
    sis_dm  = sis_profile(rPts, mass=hMass_z[i], redshift=redshift)
    sis_gas = sis_gas_profile(mass_hot=m_gas0, sis_dm=sis_dm, tables=tables)

    nfw_dm  = nfw_profile(rPts, mass=hMass_z[i], redshift=redshift)
    nfw_gas = nfw_gas_suto(mass_hot=m_gas0, nfw_dm=nfw_dm, tables=tables)
    
    ; radii, code units (kpc)
    virRad[i]  = sis_dm.r200
    virTemp[i] = sis_dm.Tvir_rVir
    
    sis_coolRad[i] = sis_gas.r_cool
    nfw_coolRad[i] = nfw_gas.r_cool
    nfw_coolRad_iso[i] = nfw_gas.r_cool_iso
  endfor
  
  ; plot
  start_PS,'halo_rads_vs_mass.eps'
    cgPlot,[0],[0],/nodata,xtitle='Halo Mass at z=2',ytitle='Radius [kpc]',$
      xrange=[8.5,12],/xs,yrange=[0,220],/ys
      
    ; mark tvir=10^4 (end of cooling tables)
    w = min(where(alog10(virTemp) ge 4.0))
    cgPlot,codeMassToLogMsun([hMass_z[w],hMass_z[w]]),[50,200],/overplot,line=1
      
    ; rvir and rcool
    cgPlot,codeMassToLogMsun(hMass_z),virRad,/overplot,color=cgColor('black'),line=2
    
    cgPlot,codeMassToLogMsun(hMass_z),sis_coolRad,/overplot,color=cgColor('orange'),line=0
    cgPlot,codeMassToLogMsun(hMass_z),nfw_coolRad,/overplot,color=cgColor('red'),line=0
    cgPlot,codeMassToLogMsun(hMass_z),nfw_coolRad_iso,/overplot,color=cgColor('blue'),line=0
    
    legend,textoidl(['R_{vir}','SIS R_{cool}','NFW poly R_{cool}','NFW iso R_{cool}']),$
      textcolors=['black','orange','red','blue'],linestyle=[2,0,0,0],box=0,/top,/right,linesize=0.25
  end_PS
  
  stop
  
end
  
; compProfiles(): for one redshift and halo mass, compare different theoretical profiles
;                 density,temp,timescales vs radius
  
pro compProfiles
  
  compile_opt idl2, hidden, strictarr, strictarrsubs
  units = getUnits()    
  
  ; config
  redshift = 2.0
  hMass    = 10.0 ; code units at redshift
  
  ; load Sutherland & Dopita (1993) cooling tables
  tables = interpLambdaSD93()
  
  ; arrays
  rRes = 2000
  rPts = linspace(0.01,5.0,rRes) ; r/rvir
  
  ; DM and gas profiles
  sis_dm  = sis_profile(rPts, mass=hMass, redshift=redshift)
  sis_gas = sis_gas_profile(mass_hot=hMass*units.f_b, sis_dm=sis_dm, tables=tables)
  
  nfw_dm  = nfw_profile(rPts, mass=hMass, redshift=redshift)
  nfw_gas = nfw_gas_suto(mass_hot=hMass*units.f_b, nfw_dm=nfw_dm, tables=tables)

  ; plot
  xrange = [0.01,1.5]
  pos = list([0.18,0.67,0.95,0.95] ,$
             [0.18,0.39,0.95,0.67] ,$
             [0.18,0.11,0.95,0.39])
  
  start_PS,'halo_denstempts_vs_rad.eps', ys=8, xs=6
    cgPlot,[0],[0],/nodata,xtitle='',ytitle=textoidl('\rho(r) / \rho_{crit,z=2}'),$
      xrange=xrange,/xs,/xlog,xminor=0,yrange=10.0^[1.1,8],/ys,/ylog,yminor=0,position = pos[0],xtickname=replicate(' ',10)
      
    ; DM profile
    cgPlot,rPts,rhoRatioToCrit(sis_dm.rho_DM,redshift=redshift),/overplot,color=cgColor('gray')
    cgPlot,rPts,rhoRatioToCrit(nfw_dm.rho_DM,redshift=redshift),/overplot,color=cgColor('black')
    
    ; r_s line for NFW
    cgPlot,[nfw_dm.r_s/nfw_dm.r200,nfw_dm.r_s/nfw_dm.r200],[1.5,3],line=2,/overplot
    
    ; gas profile
    cgPlot,rPts,rhoRatioToCrit(sis_gas.rho_gas,redshift=redshift),/overplot,line=2,color=cgColor('orange')
    cgPlot,rPts,rhoRatioToCrit(nfw_gas.rho_gas,redshift=redshift),/overplot,line=2,color=cgColor('red')
    cgPlot,rPts,rhoRatioToCrit(nfw_gas.rho_gas_iso,redshift=redshift),/overplot,line=2,color=cgColor('blue')
    
    legend,['SIS DM','NFW DM'],textcolors=['gray','black'],box=0,/right,/top
    
    ; temp
    cgPlot,[0],[0],/nodata,xtitle="",xtickname=replicate(' ',10),ytitle=textoidl('T_{gas} [log K]'),$
       xrange=xrange,/xs,/xlog,xminor=0,yrange=[5,7],/ys,/noerase,position = pos[1]
      
    cgPlot,rPts,alog10(sis_gas.temp_gas),/overplot,color=cgColor('orange')
    cgPlot,rPts,alog10(nfw_gas.temp_gas),/overplot,color=cgColor('red')
    cgPlot,rPts,alog10(replicate(nfw_gas.T_0,n_elements(rPts))),/overplot,color=cgColor('blue')
    
    legend,['SIS','NFW poly','NFW iso'],textcolors=['orange','red','blue'],box=0,/right,/top
    
    ; timescales
    cgPlot,[0],[0],/nodata,xtitle=textoidl(' r / r_{vir} '),ytitle='Timescale [Gyr]',$
      xrange=xrange,/xs,yrange=[0.001,5.0],/ys,/ylog,/xlog,xminor=0,yminor=0,/noerase,position = pos[2]
      
    cgPlot,rPts,sis_gas.coolTime,/overplot,color=cgColor('orange'),line=0
    cgPlot,rPts,sis_gas.dynTime,/overplot,color=cgColor('orange'),line=2
    cgPlot,rPts,nfw_gas.coolTime,/overplot,color=cgColor('red'),line=0
    cgPlot,rPts,nfw_gas.dynTime,/overplot,color=cgColor('red'),line=2
    cgPlot,rPts,nfw_gas.coolTime_iso,/overplot,color=cgColor('blue'),line=0
    cgPlot,rPts,nfw_gas.dynTime*0.97,/overplot,color=cgColor('blue'),line=2 ; visual offset
    
    ; verify r_cool
    y_rcool = [0.003,0.03]
    cgPlot,[sis_gas.r_cool,sis_gas.r_cool]/nfw_dm.r200,y_rcool,line=1,color=cgColor('orange'),/overplot
    cgPlot,[nfw_gas.r_cool,nfw_gas.r_cool]/nfw_dm.r200,y_rcool,line=1,color=cgColor('red'),/overplot
    cgPlot,[nfw_gas.r_cool_iso,nfw_gas.r_cool_iso]/nfw_dm.r200,y_rcool,line=1,color=cgColor('blue'),/overplot
    
    legend,textoidl(['\tau_{cool}','\tau_{dyn}']),linestyle=[0,2],box=0,/top,/left,linesize=0.25
    
  end_PS
  
  stop
end
