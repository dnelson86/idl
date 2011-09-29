; arepoSphSym.pro
; dnelson
; 7/7/10
;
; snapshot analysis for spherically symmetric cases

@arepoLoad

; temp TS blastwave
pro temp2D

  ; config - TS_BlastWave.2D
  fBase="c:\zStuff\IDL.work\Default\TS_BlastWave.2D\snap_"
  xrange     = [0,1.0]
  yrange_vel = [-5.0,5.0]
  yrange_rho = [0.0,5.0]
  xyCen      = [0.5,0.5]
  
  for i=0,15,1 do begin
  ;i=2
  
    h = loadSnapshotHDF5(fBase,i,pos,vel,id,mass,u,rho,hsml)

    r = reform( sqrt( (pos[0,*]-xyCen[0])^2.0 + (pos[1,*]-xyCen[1])^2.0 ) )
    rrange = [0,max(r)]

    !p.multi = [0,2,1]
    plot,r,vel[0,*],psym=2,xrange=rrange,yrange=yrange_vel,/ys,/xs,charsize=1.3, $
         xtitle="r",ytitle=textoidl("v")
    plot,r,rho,     psym=2,xrange=rrange,yrange=yrange_rho,/ys,/xs,charsize=1.3, $
         xtitle="r",ytitle=textoidl("\rho")
    
    print,'time=',h.time
    !p.multi = 0
    
    ;; check for moving mesh
    ;;plot,pos[0,*],pos[1,*],psym=3,xrange=xrange,yrange=yrange,/ys,/xs
    
    ;stop
    wait,0.1
  endfor
  
end


pro temp3D

  ; config - TS_BlastWave.2D
  fBase     = "c:\zStuff\IDL.work\Default\TS_BlastWave.3D\snap_"
  yRangeVel = [-5.0,5.0]
  yRangeRho = [0.0,5.0]
  xyzCenter = [0.5,0.5,0.5]
  nBins     = 100
  
  velBin  = fltarr(nBins)
  massBin = fltarr(nBins)
  uBin    = fltarr(nBins)
  rhoBin  = fltarr(nBins)

  for i=0,20,1 do begin
  ;i=2
  
    h = loadSnapshot(fBase,i,pos,vel,id,mass,u,rho,hsml)

    ; calculate radius of each point
    r = reform( sqrt( (pos[0,*]-xyzCenter[0])^2.0 + $
                      (pos[1,*]-xyzCenter[1])^2.0 + $
                      (pos[2,*]-xyzCenter[2])^2.0 ) )
    rRange = [0,max(r)]
    
    ; bin in r
    dr = ( rRange[1] - rRange[0] ) / nBins
    rPts = findgen(nBins)/nBins * ( rRange[1] - rRange[0] ) + dr/2.0
    
    for j=0,nBins-1 do begin
      rMin = j*dr
      rMax = rMin + dr
      w = where( (r ge rMin) and (r lt rMax), count)
      
      if ( count ne 0 ) then begin
        velBin[j] = mean( sqrt( vel[0,w]^2.0 + vel[1,w]^2.0 + vel[2,w]^2.0 ) )
        rhoBin[j] = mean( rho[w] )
      endif
    endfor

    ; plot
    !p.multi = [0,2,1]
    plot,rPts,velBin,psym=2,xrange=rRange,yrange=yRangeVel,/ys,/xs,charsize=1.3, $
         xtitle="r",ytitle=textoidl("v")
    plot,rPts,rhoBin,psym=2,xrange=rRange,yrange=yRangeRho,/ys,/xs,charsize=1.3, $
         xtitle="r",ytitle=textoidl("\rho")
    
    print,'time=',h.time
    !p.multi = 0
    
    ;stop
    wait,0.1
  endfor

end