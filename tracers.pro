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
  ;n_srcTargs = long( (size(Pos_SrcTargs))[2] )
  ;n_srcOrigs = long( (size(Pos_SrcOrigs))[2] )
  n_srcTargs = long( n_elements(Pos_SrcTargs[0,*]) )
  n_srcOrigs = long( n_elements(Pos_SrcOrigs[0,*]) )
  
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

; findAspectRatio():

function findAspectRatio, pos, r200, boxCen, gas_mass=gas_mass, targetN=targetN

  rf    = 10    ; fraction of r200 for cutting off-axis
  rfTar = 1.5   ; fraction of r200 to locate
  bs    = 10.0  ; kpc
  maxR  = 190.0 ; kpc
  
  ; DIRECTION ONE
  w = where(abs(pos[2,*]-boxCen) le r200/rf and abs(pos[1,*]-boxCen) le r200/rf,count)
  
  rad = reform(sqrt((pos[0,w]-boxCen)*(pos[0,w]-boxCen) + $
                    (pos[1,w]-boxCen)*(pos[1,w]-boxCen) + $
                    (pos[2,w]-boxCen)*(pos[2,w]-boxCen)))
  
  if keyword_set(gas_mass) then $
    h = hist1d(rad,gas_mass[w],binsize=bs,obin=loc,binedge=-1) ;weight by mass
  if not keyword_set(gas_mass) then $
    h = histogram(rad,binsize=bs,locations=loc) ;number density only
    
  loc += bs/2.0 
  
  rpts = findgen(1001)/1000 * maxR
  nInterp = interpol(h,loc,rpts)
        
  if (targetN eq 0) then begin
    ; locate target number density
    w = where(abs(rpts-r200/rfTar) eq min(abs(rpts-r200/rfTar)),count)
    if (count gt 1) then w = w[0]
    targetN = nInterp[w[0]]
    print,'new target ',targetN
  endif
  
  ; save radius of target number density
  w = where(abs(nInterp-targetN) eq min(abs(nInterp-targetN)),count)
  if (count gt 1) then w = w[0]
  sizeX = rpts[w[0]]

  ; DIRECTION TWO
  w = where(abs(pos[0,*]-boxCen) le r200/rf and abs(pos[1,*]-boxCen) le r200/rf,count)
  
  rad = reform(sqrt((pos[0,w]-boxCen)*(pos[0,w]-boxCen) + $
                    (pos[1,w]-boxCen)*(pos[1,w]-boxCen) + $
                    (pos[2,w]-boxCen)*(pos[2,w]-boxCen)))
 
  if keyword_set(gas_mass) then $
    h = hist1d(rad,gas_mass[w],binsize=bs,obin=loc,binedge=-1) ;weight by mass
  if not keyword_set(gas_mass) then $
    h = histogram(rad,binsize=bs,locations=loc) ;number density only
    
  loc += bs/2.0 
  
  rpts = findgen(1001)/1000 * maxR
  nInterp = interpol(h,loc,rpts)
        
  ; locate target number density and save radius
  w = where(abs(nInterp-targetN) eq min(abs(nInterp-targetN)),count)
  if (count gt 1) then w = w[0]
  sizeZ = rpts[w[0]]
  
  aspectRatio = sizeZ / sizeX
  ;stop
  return, aspectRatio
 
end

; gasSphereRadProfiles(): calculate radial profiles of density, temperature, velocity of gas

pro gasSphereRadProfiles

  ; config
  nbins = 20
  
  ;radMinMax = alog10([1.0,300.0])
  radMinMax = alog10([0.2,10.0])
  
  ;r_s  = 327.8 ;kpc (hernquist)
  ;r200 = 162.6 ;kpc (nfw)
  
  ; paths
  workingPath = '/n/home07/dnelson/dev.tracer/'
  snapPath    = workingPath + 'evrardOrig.50k/output/'
  plotBase    = 'evrardOrig.50k'
  
  ; snapshot selection
  snaps = [0,1,2,3,4]
  nSnaps = n_elements(snaps)
  
  ; arrays
  radDens    = fltarr(nSnaps,nbins)
  radDensTR  = fltarr(nSnaps,nbins)
  radDensST = fltarr(nSnaps,nbins)
  ;radDens2   = fltarr(nSnaps,nbins)
  ;radTemp    = fltarr(nSnaps,nbins)
  
  times = fltarr(nSnaps)
  
  ; aspect ratio arrays
  targetNgas = 0.0
  targetNtr  = 0.0
  
  gasAR   = fltarr(nSnaps)
  trAR    = fltarr(nSnaps)
  
  ; radial bins
  radBins = 10.0^( findgen(nbins+1)/nbins*(radMinMax[1]-radMinMax[0]) + radMinMax[0] )
  midBins = 10.0^( (findgen(nbins)+0.5)/nbins*(radMinMax[1]-radMinMax[0]) + radMinMax[0] )
  
  ; tracer mass (t=0)
  h    = loadSnapshotHeader(snapPath,snapNum=0,/verbose)
  mass = loadSnapshotSubset(snapPath,snapNum=0,partType='gas',field='mass')
  
  trMassConst = total(mass) / h.nPartTot[3]
  
  boxCen = h.boxSize / 2.0
  
  ; load
  foreach snap,snaps,k do begin
    print,'snap: ',str(snap)
    
    h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
    
    pos    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pos')
    mass   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='mass')
    ;dens   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='rho')

    ; calculate radii of particles
    rad = reform(sqrt((pos[0,*]-boxCen)*(pos[0,*]-boxCen) + $
                      (pos[1,*]-boxCen)*(pos[1,*]-boxCen) + $
                      (pos[2,*]-boxCen)*(pos[2,*]-boxCen)))

    ; load tracer positions and calculate radii
    pos_tr = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='pos')
    cumsfr_tr = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='tracer_cumsfr')
    
    rad_tr = reform(sqrt((pos_tr[0,*]-boxCen)*(pos_tr[0,*]-boxCen) + $
                         (pos_tr[1,*]-boxCen)*(pos_tr[1,*]-boxCen) + $
                         (pos_tr[2,*]-boxCen)*(pos_tr[2,*]-boxCen)))

    ; calculate aspect ratios
    ;gasAR[k] = findAspectRatio(pos,r200,boxCen,targetN=targetNgas,gas_mass=mass)
    ;trAR[k]  = findAspectRatio(pos_tr,r200,boxCen,targetN=targetNtr)
               
    ; load star positions, calculate radii
    if (h.nPartTot[4] gt 0) then begin
      pos_stars  = loadSnapshotSubset(snapPath,snapNum=snap,partType='stars',field='pos')
      mass_stars = loadSnapshotSubset(snapPath,snapNum=snap,partType='stars',field='mass')
      
      print,total(mass),total(mass_stars),total(mass)+total(mass_stars),$
            n_elements(mass),n_elements(mass_stars),n_elements(mass)+n_elements(mass_stars),$
            n_elements(cumsfr_tr)*trMassConst
      
      if (n_elements(pos_stars) gt 1) then $
        rad_st = reform(sqrt((pos_stars[0,*]-boxCen)*(pos_stars[0,*]-boxCen) + $
                             (pos_stars[1,*]-boxCen)*(pos_stars[1,*]-boxCen) + $
                             (pos_stars[2,*]-boxCen)*(pos_stars[2,*]-boxCen)))        
    endif      
    ; do binning
    for j=0, nbins-1,1 do begin
      ; shell volume normalization
      vol = 4*!pi/3 * (radBins[j+1]^3.0 - radBins[j]^3.0)
      
      w = where(rad ge radBins[j] and rad lt radBins[j+1], count)
      w_tr = where(rad_tr ge radBins[j] and rad_tr lt radBins[j+1], count_tr)
      
      ;if (n_elements(rad_st) gt 0) then $
      ;  w_st = where(rad_st ge radBins[j] and rad_st lt radBins[j+1], count_st)
    
      ;radDens2[k,j]  = mean(dens[w]) * 1e10
      if (count ne 0) then $
        radDens[k,j]    = total(mass[w]) / vol * 1e10
      if (count_tr ne 0) then $
        radDensTR[k,j]  = count_tr * trMassConst / vol * 1e10
      ;if (n_elements(w_st) gt 0) then if (count_st gt 0) then $
      ;  radDensST[k,j]  = total(mass_stars[w_st]) / vol * 1e10

      times[k] = h.time

    endfor
    
  endforeach ;snap

  ; plots
  start_PS, workingPath + plotBase + '_radDens_'+string(min(snaps),format='(I3.3)')+'_'+$
            string(max(snaps),format='(I3.3)')+'.eps'
  
    xrange = 10.0^radMinMax
    ;yrange = [min(radDens[where(radDens ne 0)])*0.5,max(radDens)*2.0]
    yrange = [1e9,5e11]
    
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="",ytitle=textoidl("\rho_{gas,tr}(r)"),$ ; [M_{sun} kpc^{-3}]
              title=plotBase,/ylog,/xlog,$
              position=[0.18,0.35,0.9,0.9],xtickname=replicate(' ',10)
    
    ; density profiles for successive snapshots
    legendStrs = []
    legendColors = []
    
    ;fsc_plot,[r200,r200],yrange,line=2,color=fsc_color('light gray'),/overplot   
    ;fsc_text,0.78,0.39,textoidl("r_{200}"),/normal,alignment=0.0,color=fsc_color('light gray'),$
    ;         charsize=!p.charsize-0.5
   
    foreach snap,snaps,k do begin
      ;if (count ne 0) then begin
        fsc_plot,midBins,radDens[k,*],  line=0,/overplot,color=getColor(k);,thick=!p.thick+0.5
        ;fsc_plot,midBins,radDens2[k,*], line=2,/overplot,color=getColor(k);,thick=!p.thick+0.5
        fsc_plot,midBins,radDensTR[k,*],line=1,/overplot,color=getColor(k);,thick=!p.thick+0.5
      ;endif
      
      ;legendStrs = [legendStrs,'t = '+string(times[k],format='(f3.1)')+' Gyr (h'+textoidl("_{gas}")+" = "+$
      ;              string(gasAR[k],format='(f4.2)')+" h"+textoidl("_{tr}")+" = "+$
      ;              string(trAR[k],format='(f4.2)')+")"]
      legendStrs = [legendStrs,'t = '+string(times[k],format='(f3.1)')]
      legendColors = [legendColors,getColor(k,/name)]      
    endforeach
    
    ; snapshot legend
    legend,legendStrs,textcolors=legendColors,box=0,/bottom,/left,charsize=!p.charsize-0.6
    legend,['gas','tracer'],linestyle=[0,1],textcolors=['black','black'],$
           box=0,/right,/top,charsize=!p.charsize-0.4
           
    ; residual plot
    yrange = [0.6,1.4]
    ;yrange = [-0.2,1.0]
    fsc_plot,[0],[0],/nodata,/noerase,xrange=xrange,yrange=yrange,/xs,/ys,/xlog,$
             xtitle="radius",$
             ytitle=textoidl("\rho_{tr} / \rho_{gas}"),position=[0.18,0.15,0.9,0.35],$
             ytickv=[0.8,1.0,1.2],yticks=2
             ;ytickv=[alog10(1.0),alog10(2.0),alog10(5.0)],$
             ;ytickname=['1','2','5'],yticks=2
  
    fsc_plot,xrange,alog10([1.0,1.0]),line=0,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,alog10([2.0,2.0]),line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,alog10([5.0,5.0]),line=1,color=fsc_color('light gray'),/overplot
    
    fsc_plot,xrange,[0.8,0.8],line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[1.2,1.2],line=1,color=fsc_color('light gray'),/overplot
    
    ;fsc_plot,[r200,r200],yrange,line=2,color=fsc_color('light gray'),/overplot
    
    ; just interpolate both onto a set of radii
    nbins = 100
    ;minmax = alog10([1.0,280.0])
    minmax = alog10([0.1,10.0])
    res_pts = 10.0^( findgen(nbins+1)/nbins * (minmax[1]-minmax[0]) + minmax[0] )
    
    foreach snap,snaps,k do begin
      gas_res = interpol(radDens[k,*],midBins,res_pts)
      tr_res  = interpol(radDensTR[k,*],midBins,res_pts)
      ;st_res  = interpol(radDensST[k,*],midBins,res_pts)
  
      fsc_plot,res_pts,tr_res/gas_res,line=0,color=getColor(k),/overplot
    endforeach

  end_PS
  stop
  ; plot two - individual snapshots
  for k=0,n_elements(snaps)-1 do begin
    print,times[k]
  start_PS, workingPath + plotBase + 'radDens_'+str(k)+'_'+string(min(snaps),format='(I3.3)')+'_'+$
            string(max(snaps),format='(I3.3)')+'.eps'
            
    xrange = 10.0^radMinMax
    yrange = [5e1,5e8]
    ;yrange = [1e-6,1e7]
    
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="radius [kpc]",ytitle=textoidl("\rho_{gas,tr}(r) [M_{sun} kpc^{-3}]"),$
              title=plotBase + textoidl(" \Lambda=0.1")+"",/ylog,/xlog
    
    ; density profiles for successive snapshots
    fsc_plot,[r200,r200],yrange,line=2,color=fsc_color('light gray'),/overplot   
    fsc_text,0.82,0.24,textoidl("r_{200}"),/normal,alignment=0.0,color=fsc_color('light gray'),$
             charsize=!p.charsize-0.5
   
    ; interpolate onto a set of radii
    nbins = 100
    res_pts = 10.0^( findgen(nbins+1)/nbins * (alog10(280)-alog10(1)) + alog10(1) )
    
    gas_res = interpol(radDens[k,*],midBins,res_pts)
    tr_res  = interpol(radDensTR[k,*],midBins,res_pts)
    st_res  = interpol(radDensST[k,*],midBins,res_pts)
    
    if (count ne 0) then begin
      fsc_plot,res_pts,gas_res,line=0,/overplot,color=getColor(1);,thick=!p.thick+0.5
      fsc_plot,res_pts,st_res, line=0,/overplot,color=getColor(2);,thick=!p.thick+0.5
      fsc_plot,res_pts,tr_res, line=0,/overplot,color=getColor(3);,thick=!p.thick+0.5
    endif
      
    ; legend
    legend,['gas','stars','tracers'],textcolors=getColor([1,2,3],/name),$
           box=0,/top,/right,charsize=!p.charsize-0.2
    legend,['t = '+string(times[k],format='(f3.1)')+' Gyr (h'+textoidl("_{gas}")+" = "+$
                    string(gasAR[k],format='(f4.2)')+" h"+textoidl("_{tr}")+" = "+$
                    string(trAR[k],format='(f4.2)')+")"],$
           textcolors=['black'],box=0,/bottom,/left,charsize=!p.charsize-0.5
            
  end_PS, pngResize=50, /deletePS
  endfor ;k
stop
end

; gasSphereRunsProf(): overplot tracer density profiles and ratios for different runs at one time

pro gasSphereRunsProf

  ; config
  nbins = 20
  r200  = 162.6 ;kpc
  radMinMax = alog10([1.0,200.0])

  ; paths
  workingPath = '/n/home07/dnelson/dev.tracer/'
  plotBase    = 'gasSphere.gastr.res'
  
  ; snapshot selection
  snap = 4
  ;runs = '2e5.f'+['1','4','10','0']
  runs = 'res.'+['2e5','1e5','5e4','1e4','5e3']
  
  ; arrays
  nRuns = n_elements(runs)
  
  radDensGas = fltarr(nRuns,nbins)
  radDensTR  = fltarr(nRuns,nbins)
  
  ; radial bins
  radBins = [0.0,           logspace(radMinMax[0],radMinMax[1],nbins)]
  midBins = [radBins[0]/2.0,logspace(radMinMax[0],radMinMax[1],nbins,/mid)]
  
  ; load
  foreach run,runs,k do begin
    print,'run: ',str(run)
    
    snapPath    = workingPath + 'gasSphere.gastr.'+run+'/output/'
    
    ; tracer mass (t=0)
    h    = loadSnapshotHeader(snapPath,snapNum=0,/verbose)
    mass = loadSnapshotSubset(snapPath,snapNum=0,partType='gas',field='mass')
    
    trMassConst = total(mass) / h.nPartTot[3]
    
    ; load
    h = loadSnapshotHeader(snapPath,snapNum=snap,/verbose)
    
    boxCen = h.boxSize / 2.0
    
    pos    = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='pos')
    mass   = loadSnapshotSubset(snapPath,snapNum=snap,partType='gas',field='mass')

    ; calculate radii of particles
    rad_gas = reform(sqrt((pos[0,*]-boxCen)*(pos[0,*]-boxCen) + $
                          (pos[1,*]-boxCen)*(pos[1,*]-boxCen) + $
                          (pos[2,*]-boxCen)*(pos[2,*]-boxCen)))

    ; load tracer positions and calculate radii
    pos_tr = loadSnapshotSubset(snapPath,snapNum=snap,partType='tracer',field='pos')
    
    rad_tr = reform(sqrt((pos_tr[0,*]-boxCen)*(pos_tr[0,*]-boxCen) + $
                         (pos_tr[1,*]-boxCen)*(pos_tr[1,*]-boxCen) + $
                         (pos_tr[2,*]-boxCen)*(pos_tr[2,*]-boxCen)))

    ; do binning
    for j=0, nbins-1 do begin
      ; shell volume normalization
      vol = 4*!pi/3 * (radBins[j+1]^3.0 - radBins[j]^3.0)
      
      w    = where(rad_gas ge radBins[j] and rad_gas lt radBins[j+1], count_gas)
      w_tr = where(rad_tr  ge radBins[j] and rad_tr  lt radBins[j+1], count_tr)
      
      if (count_gas ne 0) then $
        radDensGas[k,j] = total(mass[w]) / vol * 1e10
      if (count_tr ne 0) then $
        radDensTR[k,j]  = count_tr * trMassConst / vol * 1e10

    endfor
    
  endforeach ;snap

  ; plots
  start_PS, workingPath + plotBase + '_snap='+str(snap)+'.radDens.eps'
  
    xrange = 10.0^radMinMax
    yrange = [1e7,max([radDensTR])*1.5]
    yrange = [1e7,1e8]
    
    fsc_plot, [0],[0],/nodata,xrange=xrange,yrange=yrange,/xs,/ys,$
              xtitle="",ytitle=textoidl("\rho_{tr}(r) * R^2 [M_{sun} kpc^{-1}]"),$
              title=plotBase + textoidl(" \Lambda=0.1")+" t = "+string(h.time,format='(f3.1)')+" Gyr",$
              /ylog,/xlog,position=[0.18,0.35,0.9,0.9],xtickname=replicate(' ',10)
    
    ; density profiles for successive snapshots  
    fsc_plot,[r200,r200],yrange,line=2,color=fsc_color('light gray'),/overplot   
    fsc_text,r200*0.75,yrange[0]*2,textoidl("r_{200}"),alignment=0.0,color=fsc_color('light gray'),$
             charsize=!p.charsize-0.4
   
    ; gas profile only for first run
    fsc_plot,midBins,radDensGas[0,*]*midBins^2.0,line=1,/overplot,color=fsc_color('black')
    
    foreach run,runs,k do begin
      fsc_plot,midBins,radDensTR[k,*]*midBins^2.0,line=0,/overplot,color=getColor(k)
    endforeach
    
    ; snapshot legend
    ;strings = ['1=','4=','10=','1/4='] + textoidl('N_{tr} / N_{gas}')
    strings = ['gas 2e5',runs] ;+ textoidl('N_{gas}')
    lines = [1,intarr(Nruns)]
    colors  = getColor([0,indgen(nRuns)],/name)
    
    legend,strings,textcolors=colors,linestyle=lines,$
      box=0,/top,/right,charsize=!p.charsize-0.3,linesize=0.3
           
    ; residual plot
    yrange = [0.6,1.4]
    fsc_plot,[0],[0],/nodata,/noerase,xrange=xrange,yrange=yrange,/xs,/ys,/xlog,$
             xtitle="radius [kpc]",$
             ytitle=textoidl("\rho_{tr,res} / \rho_{tr,2e5}"),position=[0.18,0.15,0.9,0.35],$
             ytickv=[0.8,1.0,1.2],yticks=2
  
    fsc_plot,xrange,[1.0,1.0],line=0,color=fsc_color('light gray'),/overplot    
    fsc_plot,xrange,[0.8,0.8],line=1,color=fsc_color('light gray'),/overplot
    fsc_plot,xrange,[1.2,1.2],line=1,color=fsc_color('light gray'),/overplot
    
    fsc_plot,[r200,r200],yrange,line=2,color=fsc_color('light gray'),/overplot
    
    ; just interpolate both onto a set of radii
    nbins = 100
    res_pts = 10.0^( findgen(nbins+1)/nbins * (alog10(r200)-alog10(1)) + alog10(1) )
    
    tr_res_f1 = interpol(radDensTR[0,*],midBins,res_pts)
    
    for k=1,n_elements(runs)-1 do begin
      ;gas_res = interpol(radDensGas[k,*],midBins,res_pts)
      tr_res  = interpol(radDensTR[k,*],midBins,res_pts)
  
      fsc_plot,res_pts,tr_res/tr_res_f1,line=0,color=getColor(k),/overplot
    endfor
    
  end_PS
  stop
end