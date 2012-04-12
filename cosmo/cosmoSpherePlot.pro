
; cosmoSpherePlot.pro
; gas accretion project - visualization/plotting of quantities onto healpix spheres
; dnelson apr.2012

; plotHaloShellDensComp(): compare different mass halos at one radial shell and one redshift

pro plotHaloShellDensComp

  ; config
  redshift = 2
  partType = 'gas'
  radInds  = [6,11,15]  ; pre-saved radFacs
  minmax   = [-0.6,2.0] ; log (rho/mean rho)
  rot_ang  = [0,0]      ; [60,-45] ;[lat,long] center in deg (left,up)
  cutSubS  = 1          ; cut satellite substructures out from halo
  
  sP = simParams(res=512,run='gadget',redshift=float(redshift))

  ; get IDs of mass targets
  hMassTargets = [12.5,12.0,11.5,11.0]
  subgroupIDs  = massTargetToHaloID(hMassTargets,sP=sP)

  ; plot
  foreach radInd,radInds do begin
    print,radInd
    
    if cutSubS then csTag = '.cutSubS' else csTag = ''
    start_PS, sP.plotPath+'shell_rcomp_z'+str(redshift)+'_r'+str(radInd)+'_'+partType+csTag+'.eps', xs=6*1.5, ys=6
      
      pos = ['ul_nb','ur_nb','ll_nb','lr_nb']
      xtpos = [0.06,0.56,0.06,0.56]
      ytpos = [0.55,0.55,0.13,0.13]
      
      bartitle = "log ( "+textoidl("\rho / <\rho>")+" )"
        
      for i=0,3 do begin
        ; interpolate onto the shell
        hsd = haloShellDensity(sP=sP,partType=partType,subgroupID=subgroupIDs[i],cutSubS=cutSubS,/save)
  
        ; convert densities into ratios to the mean
        healpix_data = alog10(10.0^hsd.val_dens[*,radInd] / mean(10.0^hsd.val_dens[*,radInd]))
        healpix_data = reform(healpix_data)
        
        title = sP.run+" "+str(sP.res)+textoidl("^3")+"  z = "+string(sP.redshift,format='(f3.1)')+" "+$
                textoidl("\rho_{"+hsd.partType+"} (r / r_{vir} = "+string(hsd.radFacs[radInd],format='(f4.2)'))+")"
  
        if i eq 0 then $
          plotMollweideProj,healpix_data,rot_ang=rot_ang,title=title,bartitle=bartitle,$
            minmax=minmax,/bigbar,pos=pos[i]
        if i gt 0 then $
          plotMollweideProj,healpix_data,rot_ang=rot_ang,title="",bartitle="",$
            minmax=minmax,pos=pos[i],/noerase
  
        cgText,xtpos[i],ytpos[i],"M = "+string(hMassTargets[i],format='(f4.1)'),$
          /normal,alignment=0.5,color=cgColor('dark gray'),charsize=1.0
      endfor
    end_PS, pngResize=60, /deletePS
  
  endforeach ;radInds
  
end

; plotHaloShellDensMovie(): 

pro plotHaloShellDensMovie, redshift=redshift, hMassTarget=hMassTarget
  
  ; config
  rot_ang = [0,0] ; [lat,long] center in deg (left,up)
  cutSubS = 0     ; cut satellite substructures out from halo
  
  sP = simParams(res=512,run='gadget',redshift=float(redshift))  
  
  ; select a primary subgroup based on mass
  subgroupID = massTargetToHaloID(hMassTarget,sP=sP)
    
  ; movie configuration
  nFrames  = 450
  radStart = 0.01
  radEnd   = 1.999
  
  ; interpolate onto shells at a set of fixed radii
  print,'interpolating onto shells...'
  hsd_gas = haloShellDensity(sP=sP,partType='gas',subgroupID=subgroupID,cutSubS=cutSubS,/save)
  hsd_dm  = haloShellDensity(sP=sP,partType='dm',subgroupID=subgroupID,cutSubS=cutSubS,/save)  

  ; stepping in radius
  radStep    = (radEnd - radStart) / (nFrames-1)
  frameRadii = radStep * findgen(nFrames) + radStart

  ; normalize each shell by its mean
  for i=0,hsd_gas.nRadFacs-1 do $
    hsd_gas.val_dens[*,i] = alog10(10.0^hsd_gas.val_dens[*,i] / mean(10.0^hsd_gas.val_dens[*,i]))
  for i=0,hsd_dm.nRadFacs-1 do $
    hsd_dm.val_dens[*,i] = alog10(10.0^hsd_dm.val_dens[*,i] / mean(10.0^hsd_dm.val_dens[*,i]))
      
  ; interpolate
  print,'interpolating onto radii...'
  hpShell_gas = fltarr(hsd_gas.nPx,nFrames)
  hpShell_dm  = fltarr(hsd_dm.nPx,nFrames)
  
  for i=0,hsd_gas.nPx-1 do $
    hpShell_gas[i,*] = hermite(hsd_gas.radFacs,hsd_gas.val_dens[i,*],frameRadii)
  for i=0,hsd_dm.nPx-1 do $
    hpShell_dm[i,*] = hermite(hsd_dm.radFacs,hsd_dm.val_dens[i,*],frameRadii)

  print,'rendering frames...'
  for fn=0,nFrames-1 do begin
  ;fn = 100
    print,fn
    radius = frameRadii[fn]

    ; plot 2 vertically
    start_PS, sP.plotPath + 'frames/shell_z'+str(redshift)+'_m'+str(subgroupID)+'_r'+str(fn)+'.eps', xs=8*0.75, ys=8
      !p.multi = [0,2,1]
      
      healpix_data = reform(hpShell_gas[*,fn])
      title = sP.run+" "+str(sP.res)+textoidl("^3")+"  "+textoidl("\rho_{"+hsd_gas.partType+"} (r / r_{vir} = "+$
        string(radius,format='(f5.3)'))+")"
      bartitle = "log ( "+textoidl("\rho / <\rho>")+" )"
      
      plotMollweideProj,healpix_data,rot_ang=rot_ang,title=title,bartitle=bartitle,pos='top'
      
      healpix_data = reform(hpShell_dm[*,fn])
      title = sP.run+" "+str(sP.res)+textoidl("^3")+"  "+textoidl("\rho_{"+hsd_dm.partType+"} (r / r_{vir} = "+$
        string(radius,format='(f5.3)'))+")"
      bartitle = "log ( "+textoidl("\rho / <\rho>")+" )"
      
      plotMollweideProj,healpix_data,rot_ang=rot_ang,title=title,bartitle=bartitle,pos='bottom'
    
      ; halo mass and redshift label
      cgText,0.9,0.05,"log(M) = "+string(hMassTarget,format='(f4.1)'),$
        /normal,alignment=0.5,color=cgColor('dark gray'),charsize=1.0
      cgText,0.9,0.02,"z = "+string(redshift,format='(f3.1)'),$
        /normal,alignment=0.5,color=cgColor('dark gray'),charsize=1.0
    
    end_PS, pngResize=60, /deletePS
  endfor

end

; haloShellAngPowerSpec(): compute and plot the angular power spectrum
;                          e.g. the variation of the spherical harmonic coefficients

pro haloShellAngPowerSpec

  ; config
  radInds      = [6,11,15]  ; pre-saved radFacs
  hMassTargets = [12.5,12.0,11.5,11.0]
  cutSubS  = 1 ; cut satellite substructures out from halo
  redshift = 3
  
  sP = simParams(res=512,run='gadget',redshift=float(redshift))  
 
  ; locate primary subgroupID for mass target
  subgroupIDs = massTargetToHaloID(hMassTargets,sP=sP)
 
  foreach subgroupID,subgroupIDs,m do begin
    ; interpolate onto shells at a set of fixed radii
    hsd_gas = haloShellDensity(sP=sP,partType='gas',subgroupID=subgroupID,cutSubS=cutSubS,/save)
    hsd_dm  = haloShellDensity(sP=sP,partType='dm',subgroupID=subgroupID,cutSubS=cutSubS,/save)
    
    ; choose maximum spherical harmonic order l_max
    ;l_max = fix(3.0 * hsd_gas.nSide - 1)
    l_max = fix(2.0 * hsd_gas.nSide)
  
    ; plot
    start_PS, sP.plotPath + 'powerspec.shell_z'+str(redshift)+'_h'+str(subgroupID)+'.eps'
      
      l_vals = findgen(l_max)
      unit_lambda = 2*!pi / (l_vals + 0.5) ; wavelength on unit sphere
      ang_size = unit_lambda * 180.0 / !pi
      
      first_l = 1
      xrange = [first_l,max(l_vals)]
      yrange = [1e-5,1e0]
      
      !y.margin[1] += 1.0
      
      cgPlot,[0],[0],/nodata,xrange=xrange,yrange=yrange,xs=9,ys=1,/xlog,/ylog,$
        ytitle=textoidl("l(l+1)C_l/2\pi"),xtitle="l (Spherical Wavenumber)",title=""
        
      ; draw wavelength axis (linear map)
      cgAxis,xs=1,xaxis=1,xtitle="Angular Size [deg]",xrange=[ang_size[first_l],min(ang_size)]
      
      ; function lambdaAxisFunc,index,value
      ;   return,string(2*!pi/(value+0.5),format='(f3.1)')
      ; end
      ;cgAxis,ys=1,yaxis=1,ytitle="Wavelength [ckpc]",ytickformat='lambdaAxisFunc'
      
      ; for this radius, calculate angular power spectrum of gas and dm
      foreach radInd,radInds,k do begin
        print,radInd
        lambda = unit_lambda * hsd_gas.radFacs[radInd] * hsd_gas.rVir ; wavelength in ckpc
        
        healpix_data = alog10(10.0^hsd_gas.val_dens[*,radInd] / mean(10.0^hsd_gas.val_dens[*,radInd]))
        healpix_data = reform(healpix_data)
        ianafast,healpix_data,cl_gas,nlmax=l_max,/cxx,/nested,tmpdir='/n/home07/dnelson/',$
          /silent;,/won,iter_order=2
      
        healpix_data = alog10(10.0^hsd_dm.val_dens[*,radInd] / mean(10.0^hsd_dm.val_dens[*,radInd]))
        healpix_data = reform(healpix_data)
        ianafast,healpix_data,cl_dm,nlmax=l_max,/cxx,/nested,tmpdir='/n/home07/dnelson/',$
          /silent;,/won,iter_order=2
        
        ; plot power spectra
        cgPlot,l_vals,l_vals*(l_vals+1)*cl_gas/2/!pi,line=k,color=getColor(1),/overplot
        cgPlot,l_vals,l_vals*(l_vals+1)*cl_dm/2/!pi,line=k,color=getColor(2),/overplot
      endforeach
      
      ; legends
      legend,[textoidl("r/r_{vir} = ")+string(hsd_gas.radFacs[radInds],format='(f4.2)')],$
        linestyle=indgen(n_elements(radInds)),linesize=0.25,box=0,/top,/right
      legend,['gas','dm'],textcolor=getColor([1,2],/name),box=0,/bottom,/left
      legend,[sP.run+" "+str(sP.res)+textoidl("^3")+" log(M) = "+string(hMassTargets[m],format='(f4.1)')],$
        box=0,/top,/left,textcolor=['light gray']
    end_PS
  endforeach ;subgroupIDs

end
