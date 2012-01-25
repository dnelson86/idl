; tracers.pro
; dev for tracer particles
; dnelson jan.2012

function interpNN, pos_src, pos_dest, val_src

    ; something fancy (linear or cubic?)
    ;dens_interp = interpol(dens_gas,pos_gas,pos_tracer)
    ;res = dens_interp - dens_tracer
    
    ; nearest neighbor interpolation by indices
    ;inds = round(pos_gas / 20.0 * n_elements(dens_tracer))
    
    ; nearest neighbor interpolation (manual)
    val_interp = fltarr(n_elements(pos_src))
    
    for i=0,n_elements(pos_src)-1 do begin
      ind = where(abs(pos_dest-pos_src[i]) eq min(abs(pos_dest-pos_src[i])),count)
      
      if (count eq 0) then begin
        print,'WARNING'
        stop
      endif
      val_interp[i] = val_src[ind[0]]
    endfor
    
  return, val_interp
end

; calcHSML(): use CalcHSML external C-routine for the tree, neighbor searchings, and smoothing lengths

function calcHSML, Pos, ndims=ndims, nNGB=nNGB, boxSize=boxSize

  ; prepare inputs
  npos = (size(pos))[2]

  NumPart = long(npos)
  Mass    = fltarr(npos)+1.0 ;dummy
  
  DesNumNgb    = long(nNGB) ; number of neighbors to use
  DesNumNgbDev = long(0)
  boxSize      = float(boxSize)
  HsmlGuess    = float(1.0)
  Softening    = float(1.0)
  
  hsml_out = fltarr(NumPart)
  
  ; call CalcHSML
  libName = '/n/home07/dnelson/idl/CalcHSML/CalcHSML_'+str(ndims)+'D.so'
  ret = Call_External(libName, 'CalcHSML', $
                      NumPart,Pos,Mass,DesNumNgb,DesNumNgbDev,boxSize,HsmlGuess,Softening,hsml_out, $
                      /CDECL)
   
  return, hsml_out
                     
end

; calcNN(): use CalcNN external C-routine for the tree and neighbor search
;           return the index of Pos_SrcTargs closest to each of Pos_SrcOrigs

function calcNN, Pos_SrcTargs, Pos_SrcOrigs, boxSize=boxSize, ndims=ndims

  ; prepare inputs
  n_srcTargs = long( (size(Pos_SrcTargs))[2] )
  n_srcOrigs = long( (size(Pos_SrcOrigs))[2] )
  
  boxSize   = float(boxSize)
  
  ind_out = lonarr(n_srcOrigs)
  
  ; call CalcNN
  libName = '/n/home07/dnelson/idl/CalcNN/CalcNN_'+str(ndims)+'D.so'
  ret = Call_External(libName, 'CalcNN', $
                      n_srcTargs,n_srcOrigs,Pos_SrcTargs,Pos_SrcOrigs,boxSize,ind_out, $
                      /CDECL)
    
  return, ind_out
  
end

; calcSphMap(): use CalcSphMap external C-routine to calculate a map of projected densities
;               with the sph spline kernel

function calcSphMap, pos, hsml, mass, axes=axes, boxSize=boxSize, boxCen=boxCen, $
                     nPixels=nPixels, ndims=ndims

  ; prepare inputs
  npos = (size(pos))[2]
  
  boxSize  = float(boxSize)
  boxCen   = float(boxCen)
  nPixels  = fix(nPixels)
  axes     = fix(axes)
  
  mode     = fix(1)
  periodic = fix(1)
  
  ; check inputs
  if (n_elements(boxSize) ne 3) then stop
  if (n_elements(boxCen)  ne 3) then stop
  if (n_elements(nPixels) ne 2) then stop
  if (n_elements(axes)    ne 2) then stop
  
  if (mode ne 1 and mode ne 2 and mode ne 3) then stop
  if (periodic ne 0 and periodic ne 1) then stop
  
  if (axes[0] ne 0 and axes[0] ne 1 and axes[0] ne 2) then stop
  if (axes[1] ne 0 and axes[1] ne 1 and axes[1] ne 2) then stop
  
  ; make return
  grid_out = fltarr(nPixels[0],nPixels[1])
  
  ; call CalcSphMap
  libName = '/n/home07/dnelson/idl/CalcSphMap/CalcSphMap_'+str(ndims)+'D.so'
  ret = Call_External(libName, 'CalcSphMap', $
                      npos,pos,hsml,mass,grid_out,boxSize[0],boxSize[1],boxSize[2],$
                      boxCen[0],boxCen[1],boxCen[2],axes[0],axes[1],nPixels[0],nPixels[1],$
                      mode,periodic,/CDECL)

  return, grid_out
  
end

; estimateDensityTophat(): spatial density estimator for an input position tuple array and
;                          mass array for some particle type
;   do density calculation by calculating smoothing lengths for all the particles
;   

function estimateDensityTophat, pos, mass=mass, ndims=ndims, nNGB=nNGB, boxSize=boxSize

  if (not keyword_set(nNGB) or not keyword_set(ndims)) then begin
    print,'Error: Must specify both nNGB and ndims.'
    return,0
  endif

  ; calculate smoothing lengths
  hsml_out = calcHSML(pos,ndims=ndims,nNGB=nNGB,boxSize=boxSize)
                      
  ; estimate densities on eval_pos using hsml
  if (ndims eq 1) then $
    eval_dens = nNGB / hsml_out
    
  if (ndims eq 2) then $
    eval_dens = nNGB / (!pi * hsml_out^2.0)
    
  if (ndims eq 3) then $
    eval_dens = nNGB / (4.0*!pi/3.0 * hsml_out^3.0)
  
  ; add in mass
  if keyword_set(mass) then $
    eval_dens *= mass[0]
  
  return,eval_dens

end

; gasSphereRadProfiles(): calculate radial profiles of density, temperature, velocity of gas

pro gasSphereRadProfiles

  ; config
  nbins = 30
  
  radMinMax = alog10([1.0,600.0])
  
  r_s = 327.8 ;kpc
  r_s = 152.1 ;kpc M100
  
  ; paths
  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath    = workingPath + 'gasSphere.gasonly.2e5.cooling/output/'
  plotBase    = 'gasSphere.2e5.cooling.'
  
  ; snapshot selection
  ;snaps = [0,5,10,22]
  snaps = [0,1,2,3,4,5,6,7,8,9,10,11]
  nSnaps = n_elements(snaps)
  
  ; arrays
  radDens    = fltarr(nSnaps,nbins)
  radDens2   = fltarr(nSnaps,nbins)
  radDensTR  = fltarr(nSnaps,nbins)
  ;radNumDens = fltarr(nSnaps,nbins)
  ;radTemp    = fltarr(nSnaps,nbins)
  
  times = fltarr(nSnaps)  
  
  ; radial bins
  radBins = 10.0^( findgen(nbins+1)/nbins*(radMinMax[1]-radMinMax[0]) + radMinMax[0] )
  midBins = 10.0^( (findgen(nbins)+0.5)/nbins*(radMinMax[1]-radMinMax[0]) + radMinMax[0] )
  
  ; load
  foreach snap,snaps,k do begin
    print,'snap: ',str(snap)
    h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
    
    pos    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pos')
    dens   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='rho')
    mass   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='mass')
    ;u      = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='u')

    trMassConst = total(mass) / h.nPartTot[3]

    ; calculate radii of particles
    boxCen = h.boxSize / 2.0
    
    rad = reform(sqrt((pos[0,*]-boxCen)*(pos[0,*]-boxCen) + $
                      (pos[1,*]-boxCen)*(pos[1,*]-boxCen) + $
                      (pos[2,*]-boxCen)*(pos[2,*]-boxCen)))
                      
    ; load tracer positions and calculate radii
    pos_tr = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='pos')
    
    rad_tr = reform(sqrt((pos_tr[0,*]-boxCen)*(pos_tr[0,*]-boxCen) + $
                         (pos_tr[1,*]-boxCen)*(pos_tr[1,*]-boxCen) + $
                         (pos_tr[2,*]-boxCen)*(pos_tr[2,*]-boxCen)))
                                    
    ; do binning
    for j=0, nbins-1,1 do begin
      ; shell volume normalization
      vol = 4*!pi/3 * (radBins[j+1]^3.0 - radBins[j]^3.0)
      
      w = where(rad ge radBins[j] and rad lt radBins[j+1], count)
      w_tr = where(rad_tr ge radBins[j] and rad_tr lt radBins[j+1], count_tr)
    
      if (count ne 0.0) then begin
        radDens2[k,j]   = mean(dens[w]) * 1e10
        radDens[k,j]    = total(mass[w]) / vol * 1e10
        radDensTR[k,j]  = count_tr * trMassConst / vol * 1e10
        ;radNumDens[k,j] = count
        ;radTemp[k,j]    = mean(u[w])
        times[k]        = h.time
      endif
      
    endfor
    
  endforeach ;snap
  
  ; plots
  start_PS, workingPath + plotBase + 'radDens_'+string(min(snaps),format='(I3.3)')+'_'+$
            string(max(snaps),format='(I3.3)')+'.eps'
  
    xrange = 10.0^radMinMax
    yrange = [5e0,max(radDens)*1.5]
    
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="",ytitle=textoidl("\rho_{gas,tr}(r) [M_{sun} kpc^{-3}]"),$
              title=plotBase + " radial density profiles",/ylog,/xlog,$
              position=[0.18,0.35,0.9,0.9],xtickname=replicate(' ',10)
    
    ; density profiles for successive snapshots
    legendStrs = []
    legendColors = []
    
    fsc_plot,[r_s,r_s],yrange,line=0,color=fsc_color('dark gray'),/overplot   
    fsc_text,0.81,0.39,textoidl("r_s"),/normal,alignment=0.0,color=fsc_color('dark gray'),$
             charsize=!p.charsize-0.5
   
    foreach snap,snaps,k do begin
      w = where(radDens[k,*] gt 1e-8,count)
      
      if (count ne 0) then begin
        fsc_plot,midBins[w],radDens[k,w],  line=0,/overplot,color=getColor(k);,thick=!p.thick+0.5
        ;fsc_plot,midBins[w],radDens2[k,w], line=2,/overplot,color=getColor(k);,thick=!p.thick+0.5
        fsc_plot,midBins[w],radDensTR[k,w],line=1,/overplot,color=getColor(k);,thick=!p.thick+0.5
      endif
      
      legendStrs = [legendStrs,'t = '+string(times[k]*1000.0,format='(I3)')+' Myr']
      legendColors = [legendColors,getColor(k,/name)]      
    endforeach
    
    ; snapshot legend
    legend,legendStrs,textcolors=legendColors,box=0,position=[80.0,yrange[1]*0.6],charsize=!p.charsize-0.5
    legend,['gas','tracer'],linestyle=[0,1],textcolors=['black','black'],$
           box=0,/left,/bottom,charsize=!p.charsize-0.5
           
    ; residual plot
    yrange = [0.3,1.7]
    fsc_plot,[0],[0],/nodata,/noerase,xrange=xrange,yrange=yrange,/xs,/ys,/xlog,$
             xtitle="radius [kpc]",$
             ytitle=textoidl("\rho_{tr} / \rho_{gas}"),ytickv=[0.5,1.0,1.5],yticks=2,$
             position=[0.18,0.15,0.9,0.35]
  
    fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[1.1,1.1],line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[0.9,0.9],line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[1.5,1.5],line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[0.5,0.5],line=1,color=fsc_color('light gray'),/overplot
    
    fsc_plot,[r_s,r_s],yrange,line=0,color=fsc_color('dark gray'),/overplot
    
    ; just interpolate both onto a set of radii
    nbins = 100
    res_pts = 10.0^( findgen(nbins+1)/nbins * (alog10(400)-alog10(1)) + alog10(1) )
    
    foreach snap,snaps,k do begin
      gas_res = interpol(radDens[k,*],midBins,res_pts)
      tr_res  = interpol(radDensTR[k,*],midBins,res_pts)
  
      fsc_plot,res_pts,tr_res/gas_res,line=0,color=getColor(k),/overplot
    endforeach
    
  end_PS
stop
end
