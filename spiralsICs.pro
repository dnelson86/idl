; spiralsICs.pro
; analysis of ICs and analytic results
; dnelson feb.2011

; comparisons
; -----------

pro compTwoICs;, f1, f2

  basePath = '/n/home07/dnelson/make.ics/CompareTest/'
  f1 = 'LCold/LCold.dat'
  f2 = 'LCMGM/LCMGM.dat'

  vlims=[[-400,400],[-5,5],$
         [-400,400],[-5,5],$
         [-4.0,4.0],[-0.1,0.1]]

  h1 = loadSnapshot(basePath+f1,'none',s_pos1,s_vel1,s_id1)
  h2 = loadSnapshot(basePath+f2,'none',s_pos2,s_vel2,s_id2)
  
  for i=0,2 do begin
    print,'total diff pos ',i,total(abs(s_pos1[i,*]-s_pos2[i,*]))
    print,'total diff vel ',i,total(abs(s_vel1[i,*]-s_vel2[i,*]))   
  endfor

  for i=0,2 do begin
  
    start_PS, basePath+'pos_'+str(i)+'.eps'
      fsc_plot,s_pos1[i,*],s_pos2[i,*],psym=3,xtitle=f1+" pos"+str(i),ytitle=f2+" pos"+str(i)
      fsc_plot,[-100.0,100.0],[-100.0,100.0],line=1,/overplot
    end_PS
    
    start_PS, basePath+'vel_'+str(i)+'.eps'
      !p.multi = [0,1,2]
      
      fsc_plot,s_vel1[i,*],s_vel2[i,*],psym=3,xtitle=f1+" vel"+str(i),ytitle=f2+" vel"+str(i),$
               xrange=vlims[*,i*2],yrange=vlims[*,i*2],/xs,/ys,charsize=0.6
      
      fsc_plot,s_vel1[i,*],s_vel2[i,*],psym=4,xtitle=f1+" vel"+str(i),ytitle=f2+" vel"+str(i),$
               xrange=vlims[*,i*2+1],yrange=vlims[*,i*2+1],/xs,/ys,charsize=0.6
      fsc_plot,[-300.0,300.0],[-300.0,300.0],line=1,/overplot
      
      !p.multi = 0
    end_PS
  
  endfor

end

; exponential disk
; ----------------

function func_LambdaCrit_ExpDisk, R, R_H

  ;disk fraction
  f_d = 0.02

  ; = 2pi * R^3/R_H^2 * f * ( (1-f)(1+R/R_H) + f(R/R_H)^2 )^(-1)
  lambdaCrit = 2 * !pi * R^3.0 / R_H^2.0 
  lambdaCrit *= exp(-R/R_H)*(1 + R/R_H) 
  lambdaCrit /= ((1 - exp(-R/R_H)*(1 + R/R_H)) * (1+R/R_H) + exp(-R/R_H)*(1 + R/R_H)*(R/R_H)^2.0)

  lambdaCrit *= f_d

  return, lambdaCrit
end

function func_Xparam_ExpDisk, R, R_H, m

  k_crit = 2 * !pi / func_LambdaCrit_ExpDisk(R,R_H)

  Xparam = k_crit * R / m
  
  return,Xparam
end

; plotExpDisk()

pro plotExpDisk

  ;config
  minR = 0.01     ;kpc
  maxR = 15.5    ;kpc
  resR = 100.0
  R_H  = [3.13, 2.54, 1.0, 3.72, 6.0] ;iterative result from MakeGalaxyDisk
  modes = [1,2,3,4]
  
  R = findgen(resR)/resR * (maxR-minR) + minR

  ; plot lambda crit
  PS_Start, FILENAME="~/spirals/lambdaCrit.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=5.0, /inches

    fsc_plot,[0],[0],xtitle="r [kpc]",ytitle=textoidl("\lambda_{crit}")+"(R,R"+textoidl("_d")+") [kpc]", $
             charsize=1.2,yrange=[0.0,0.6],/ys,xrange=[0,15],/xs,/nodata
    fsc_plot,R,func_LambdaCrit_ExpDisk(R,R_H[0]),/overplot,line=0
    fsc_plot,R,func_LambdaCrit_ExpDisk(R,R_H[1]),/overplot,line=0,color=fsc_color('sky blue')
    fsc_plot,R,func_LambdaCrit_ExpDisk(R,R_H[2]),/overplot,line=1
    fsc_plot,R,func_LambdaCrit_ExpDisk(R,R_H[3]),/overplot,line=0,color=fsc_color('salmon')
    fsc_plot,R,func_LambdaCrit_ExpDisk(R,R_H[4]),/overplot,line=1   
    
    fsc_text,11,0.43,'3.13',/data
    fsc_text,10,0.33,'2.54',/data
    fsc_text,8,0.04,'1.0',/data
    fsc_text,8,0.50,'3.72',/data
    fsc_text,4,0.37,'6.0',/data

  PS_End
  
  ; plot X parameter for our R_H
  PS_Start, FILENAME="~/spirals/xParam.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=5.0, /inches

    fsc_plot,R,func_Xparam_ExpDisk(R,R_H[0],modes[0]),xtitle="r [kpc]", $
             ytitle="X"+textoidl("_m")+"(R,R"+textoidl("_d")+")",charsize=1.2, $
             yrange=[0,200],/ys,xrange=[0,15],/xs
    fsc_plot,R,func_Xparam_ExpDisk(R,R_H[0],modes[1]),/overplot
    fsc_plot,R,func_Xparam_ExpDisk(R,R_H[0],modes[2]),/overplot
    fsc_plot,R,func_Xparam_ExpDisk(R,R_H[0],modes[3]),/overplot
    
    fsc_plot,R,func_Xparam_ExpDisk(R,R_H[1],modes[0]),/overplot,line=0,color=fsc_color('sky blue')
    fsc_plot,R,func_Xparam_ExpDisk(R,R_H[1],modes[1]),/overplot,line=1,color=fsc_color('sky blue')
    fsc_plot,R,func_Xparam_ExpDisk(R,R_H[1],modes[2]),/overplot,line=0,color=fsc_color('sky blue')
    fsc_plot,R,func_Xparam_ExpDisk(R,R_H[1],modes[3]),/overplot,line=1,color=fsc_color('sky blue')
    
    fsc_plot,R,func_Xparam_ExpDisk(R,R_H[3],modes[0]),/overplot,line=0,color=fsc_color('salmon')
    fsc_plot,R,func_Xparam_ExpDisk(R,R_H[3],modes[1]),/overplot,line=1,color=fsc_color('salmon')
    fsc_plot,R,func_Xparam_ExpDisk(R,R_H[3],modes[2]),/overplot,line=0,color=fsc_color('salmon')
    fsc_plot,R,func_Xparam_ExpDisk(R,R_H[3],modes[3]),/overplot,line=1,color=fsc_color('salmon')    
    
    fsc_text,3,90,'m=1',/data,alignment=1
    fsc_text,3,48,'2',/data,alignment=1
    fsc_text,3,32,'3',/data,alignment=1
    fsc_text,3,11,'4',/data,alignment=1

  PS_End
  
end 


; plotQAndRotCurve()

pro plotQAndRotCurve

  units = getUnits()

  ;load
  fileBase = '/n/home07/dnelson/spirals/lambdacrit/ICs/'
  fileNames = ['LC_10m_disk_2/','LC_10m_disk_4/','LC_10m_disk_6/','LC_10m_disk_8/']
  ;fileNames = ['bar_1m_a/']
  
  colors = ['black','orange','sky blue','olive']
  
  ; plot toomre Q
  
  start_PS, "~/spirals/lambdacrit/toomreQAndRotCurve.eps", xs=6.0, ys=7.0
    
    ym = !y.margin
    xm = !x.margin
    !y.margin = [0]
    
    fsc_plot,[0],[0],xtitle="",ytitle="Stellar Q",position=[0.15,0.6,0.9,0.9],$
              xrange=[0,5.0],/xs,yrange=[0.8,2.5],/ys,/nodata,xtickname=replicate(' ',10)
              
    for i=0,n_elements(fileNames)-1 do begin
      h1 = loadMDGCurve(fileBase+fileNames[i]+"curve.txt",vcs,qs)
      h2 = loadSimParams(fileNames[i])
      fsc_plot,qs.list_R/h2.h,qs.Q,/overplot,color=fsc_color(colors[i])
    endfor
    
    fsc_plot,[0,15],[1,1],line=1,/overplot
    fsc_plot,[0,15],[1.2,1.2],line=1,/overplot
 
    
    !y.margin = [4,0]
    
    fsc_plot,[0],[0],xtitle="Radius / Disk Scale Length",ytitle="Rotation Velocity [km/s]", $
              xrange=[0,5.0],/xs,yrange=[10,max(sqrt(vcs.VC2_1))*1.1],/ys,/nodata,$
              position=[0.15,0.15,0.9,0.6],/noerase
    
    for i=0,n_elements(fileNames)-1 do begin
      h1 = loadMDGCurve(fileBase+fileNames[i]+"curve.txt",vcs,qs)
      h2 = loadSimParams(fileNames[i])
      fsc_plot,vcs.R/h2.h,sqrt(vcs.VC2_1),/overplot,color=fsc_color(colors[i]),line=0
      fsc_plot,vcs.R/h2.h,sqrt(vcs.VC2_2),/overplot,color=fsc_color(colors[i]),line=1
      fsc_plot,vcs.R/h2.h,sqrt(vcs.VC2_3),/overplot,color=fsc_color(colors[i]),line=2
      
      fsc_text,0.5,250-20*i,"LC-"+str((i+1)*2),alignment=0.5,color=fsc_color(colors[i]),charsize=1.2,/data
    endfor
    
    !y.margin = ym
    !x.margin = xm
     
  end_PS
  
end

; exponential disk + hernquist DM profile
; ---------------------------------------

function func_kappa2_ExpHern, G, r, M_DM, M_star, h, a

  y = r / (2.0*h)

  I0 = BESELI(y,0,/double)
  I1 = BESELI(y,1,/double)
  K0 = BESELK(y,0,/double)
  K1 = BESELK(y,1,/double)
  
  kappa2 = G*M_DM*(3*a + r) / ( r*(a+r)^3.0 ) + G*M_star / (2.0 * h^4.0) * $
           ( (4*h*I0 + r*I1)*K0 - (r*I0 + 2*h*I1)*K1 )
  
  return, kappa2
end

function func_omega2_ExpHern, G, r, M_DM, M_star, h, a

  y = r / (2.0*h)

  I0 = BESELI(y,0,/double)
  I1 = BESELI(y,1,/double)
  K0 = BESELK(y,0,/double)
  K1 = BESELK(y,1,/double)
  
  omega2 = G*M_DM / (r*(a+r)^2.0) + G*M_star / (2.0 * h^3.0) * (I0 * K0 - I1 * K1)
  
  return, omega2
end

function func_LambdaCrit_ExpHern, G, r, M_DM, M_star, h, a
  sigmar = M_star / (2*!pi*h^2.0) * exp(-r / h)
  kappa2 = func_kappa2_ExpHern(G,r,M_DM,M_star,h,a)
  lambdaCrit = 4 * !pi^2.0 * G * sigmar / kappa2 

  return, lambdaCrit
end

function func_Xparam_ExpHern, G, r, M_DM, M_star, h, a, m
  k_crit = 2 * !pi / func_LambdaCrit_ExpHern(G,r,M_DM,M_star,h,a)
  Xparam = k_crit * r / m
  
  return,Xparam
end

function func_veldispR_ExpHern, G, r, M_DM, M_star, h, a, z0, thick=thick, thin=thin

  veldispR2_thin  = (G*M_star) / (2*h^2.0) * exp(-r / h)
  veldispR2_thick = (G*M_star*z0) / (h^2.0) * exp(-r / h) + (2.773 * G * M_DM * z0^2.0) / (r*(r+a)^2.0)

  ;choose thin or thick
  if keyword_set(thick) then begin
    veldispR = veldispR2_thick
  endif
  if keyword_set(thin) then begin
    veldispR = veldispR2_thin
  endif
    
  return, velDispR
end

function func_Q_ExpHern, G, r, M_DM, M_star, h, a, z0, RDF, thick=thick,thin=thin

  ; MakeDiskGalaxy
  ;RDF = 1.0

  veldispR = func_veldispR_ExpHern(G,r,M_DM,M_star,h,a,z0,thick=thick,thin=thin)
  sigmar   = M_star / (2*!pi*h^2.0) * exp(-r / h)
  kappa2   = func_kappa2_ExpHern(G,r,M_DM,M_star,h,a)
  
  Q = RDF * sqrt(veldispR) * sqrt(kappa2) / (3.36 * G * sigmar)

  return, Q
end

pro plotExpHern

  workingPath = "c:\zStuff\IDL\Default\"

  ;config
  minR = 0.01    ;kpc
  maxR = 15.5    ;kpc
  resR = 200.0
  
  R    = findgen(resR)/resR * (maxR-minR) + minR
  Rlog = 2.0^R / max(2.0^R) * (maxR-minR) 
  
  M_DM   = 1e12                         ;mass of Hernquist halo, Msun
  M_star = 1e10                         ;mass of stellar disk, Msun
  R_H    = [3.13, 2.54, 3.72, 1.6, 6.0] ;iterative result from MakeGalaxyDisk, kpc
  a      = 100.0                        ;Hernquist scale length, kpc
  z0     = 0.2                          ;disc vertical scale height, fraction of R_H
  RDF    = 1.0
  modes  = [1,2,3]  
 
  ;convert to cgs
  kpc_in_m   = 3.08568e19
  Msun_in_kg = 1.989e30
  yr_in_s    = 31556926
  
  G = 4.30135 ;in kpc^3/Msun/s^2
  
  colors = fsc_color(['black','sky blue','orange','salmon','dark green'])
  
  ; plot epicyclic frequency
  PS_Start, FILENAME=workingPath+"kappa2_ExpHern.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=5.0, /inches
            
    fsc_plot,[0],[0],xtitle="r [kpc]",ytitle=textoidl("\kappa^2 [1/s^2]"),charsize=1.2, $
              xrange=[0.1,15],/xlog,/xs,yrange=[1e9,1e11],/ys,/nodata

    for i=0,n_elements(R_H)-1 do begin
      fsc_plot,Rlog,func_kappa2_ExpHern(G,Rlog,M_DM,M_star,R_H[i],a),/overplot,line=0,color=colors[i]
    endfor
            
  PS_End
  
  ; plot lambda crit
  PS_Start, FILENAME=workingPath+"lambdaCrit_ExpHern.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=5.0, /inches
            
    fsc_plot,[0],[0],xtitle="r [kpc]",ytitle=textoidl("critical lambda [kpc]"),charsize=1.2, $
              xrange=[0,15],/xs,yrange=[0,3],/ys,/nodata
              
    for i=0,n_elements(R_H)-1 do begin
      fsc_plot,R,func_LambdaCrit_ExpHern(G,R,M_DM,M_star,R_H[i],a),/overplot,line=0,color=colors[i]
    endfor
    
  PS_End
  
  ; plot X parameter for our R_H
  PS_Start, FILENAME=workingPath+"XParam_ExpHern.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=5.0, /inches

    fsc_plot,[0],[0],xtitle="r [kpc]", $
             ytitle="X"+textoidl("_m")+"(R,R"+textoidl("_d")+")",charsize=1.2, $
             yrange=[0,60],/ys,xrange=[0,15],/xs,/nodata
             
    for i=0,n_elements(R_H)-3 do begin
      fsc_plot,R,func_Xparam_ExpHern(G,R,M_DM,M_star,R_H[i],a,modes[0]),line=0,/overplot,color=colors[i]
      fsc_plot,R,func_Xparam_ExpHern(G,R,M_DM,M_star,R_H[i],a,modes[1]),line=1,/overplot,color=colors[i]
      fsc_plot,R,func_Xparam_ExpHern(G,R,M_DM,M_star,R_H[i],a,modes[2]),line=2,/overplot,color=colors[i]
    endfor

    fsc_text,11,50,'m=1',/data,alignment=1
    fsc_text,11,25,'2',/data,alignment=1
    fsc_text,11,10,'3',/data,alignment=1

  PS_End
  
  ; plot Q parameter (need to choose thin or thick disc!)
  PS_Start, FILENAME=workingPath+"Q_ExpHern.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=5.0, /inches

    fsc_plot,[0],[0],xtitle="r [kpc]",ytitle=textoidl("Toomre's Q"),charsize=1.2, $
              xrange=[0,15],/xs,yrange=[0.5,6.0],/ys,/nodata
              
    for i=0,n_elements(R_H)-1 do begin
      fsc_plot,R,func_Q_ExpHern(G,R,M_DM,M_star,R_H[i],a,z0,RDF,/thick),/overplot,line=0,color=colors[i]
      fsc_plot,R,func_Q_ExpHern(G,R,M_DM,M_star,R_H[i],a,z0,RDF,/thin),/overplot,line=1,color=colors[i]
    endfor    
    
    fsc_plot,[0,maxR],[1.0,1.0],line=2,/overplot
  
  PS_End
  
end

; Hernquist / NFW and parameter mapping
; -------------------------------------

; f(c) of the concentration parameter
function func_fc, c
  fc = c * (0.5 - 0.5 / (1+c)^2.0 - alog(1+c) / (1+c)) / ( alog(1+c)-c / (1+c) )^2.0
  return,fc
end

; hernquist scale length given NFW scale length and c
function func_a, r_s, c
  a = r_s * sqrt( 2*(alog(1+c) - c/(1+c)) )
  return,a
end

; disk scale length - "first guess"
function func_h0, lambda, c, r200
  h0 = sqrt(2)/2.0 * lambda / func_fc(c) * r200
  return,h0
end

; mapInputParams(): change from MakeDiskGalaxy inputs to theory model parameter set
function mapInputParams, c, v200, f_d, f_z

  ; get G, H0
  units = getUnits()
  
  ; make sure floats
  c    = float(c)
  v200 = float(v200)
  
  ; constant input params  
  lambda = 0.033 ; angular momentum parameter of halo
  j_d    = 0.02  ; angular momentum of disk (fraction of halo J)

  ; calculate M_DM, M_star, a
  M200   = v200^3.0 / (10*units.G*units.H0)
  r200   = v200 / (10*units.H0)
  r_s    = r200 / c
  M_star = f_d * M200
  M_DM   = M200 - M_star ;note M_star << M200
  a      = func_a(r_s,c)
  
  ; calculate first guesses for h, z0 (z0->z0/2 mapping!)
  h  = func_h0(lambda,c,r200)
  z0 = f_z * h / 2.0

  pS = {M_DM:M_DM,M_star:M_star,h:h,a:a,z0:z0}

  return, pS
end

; compLambdaCrit

pro compLambdaCrit

  ; set G to code units
  units = getUnits()
  G = units.G

  f_z  = 0.1 ; constant

  ; radial calc config
  minR = 0.01    ; kpc (code units)
  maxR = 30.0    ; kpc (code units)
  resR = 300.0
  R    = findgen(resR)/resR * (maxR-minR) + minR

  ; config 
  simNames = ['LC_10m_disk_2','LC_10m_disk_4','LC_10m_disk_6','LC_10m_disk_8']

  fileBase    = '/n/home07/dnelson/spirals/lambdacrit/ICs/'
  filePathsVD = fileBase + simNames + '/deldisp.dat'
  workingPath = '/n/home07/dnelson/spirals/lambdacrit/'
  
  nComps = n_elements(filePathsVD)

  ; plot
  start_PS, workingPath+"compLambdaCrit.eps"
    fsc_plot,[0],[0],ytitle="Lambda Crit [kpc]",xtitle="Radius / Disk Scale Length",charsize=1.2, $
             xrange=[0,5],/xs,yrange=[0.0,10.0],/ys,/nodata
    
    ; legend
    legendXStart  = 4.1
    legendXLength = 0.25
    legendYStart  = 9.4
    legendYSpace  = 0.4
    
    for i=0,nComps-1 do begin
      ; load sim params
      h2 = loadSimParams(simNames[i])
      shortName = strmid(simNames[i],7,6)
      
      ; map to pS
      pS = mapInputParams(h2.c,h2.v200,h2.f_d,f_z)
      
      ; theoretical
      sigmaR  = ps.M_star / (2*!pi*h2.h^2.0) * exp(-1.0 * R / h2.h)
      tKappa2 = func_kappa2_ExpHern(G,R,pS.M_DM,pS.M_star,h2.h,pS.a)  
      
      ; save lambda_crit
      lambda_crit = 4 * !pi^2.0 * G * sigmaR / tKappa2
    
      print,simNames[i] + " - " + units.colors[i],max(lambda_crit)
    
      ; calculate directly from ICs
      h3 = loadVelDispDump(filePathsVD[i],velDisp)
      
      sigmaR_meas = ps.M_star / (2*!pi*h2.h^2.0) * exp(-1.0 * velDisp.list_R / h2.h)
      
      lambda_crit2 = 4 * !pi^2.0 * G * sigmaR_meas / velDisp.epi_kappa2
      
      ; plot theoretical
      fsc_plot,R/h2.h,lambda_crit,line=0,color=fsc_color(units.colors[i]),/overplot
      
      ; plot using kappa2 measured from ICs
      fsc_plot,velDisp.list_R/h2.h,lambda_crit2,line=1,color=fsc_color(units.colors[i]),/overplot
      
      ; plot theoretical, normalized as lambda_crit/h
      fsc_plot,R/h2.h,lambda_crit/h2.h,line=1,color=fsc_color(units.colors[i]),/overplot,thick=!p.thick-2
      
      ; plot legend
      fsc_plot,[legendXStart,legendXStart+legendXLength],$
               [legendYStart-legendYSpace*i,legendYStart-legendYSpace*i],$
               color=fsc_color(units.colors[i]),line=line,thick=2.0,/overplot
      fsc_text,legendXStart+legendXLength*1.2,legendYStart-(legendYSpace/5)-legendYSpace*i,$
               shortName,charsize=1.0,alignment=0.0,/data          
      
    endfor
    
  end_PS
  
  ; plot
  start_PS, workingPath+"compSurfDens.eps"
    fsc_plot,[0],[0],ytitle="Surface Density [Msun/pc^2]",xtitle="Radius / Disk Scale Length",charsize=1.2, $
             xrange=[0.01,5],/xs,yrange=[1.0,800.0],/ys,/nodata
  
    ; legend
    legendXStart  = 4.1
    legendXLength = 0.25
    legendYStart  = 580
    legendYSpace  = 50
    
    for i=0,nComps-1 do begin
      ; load sim params
      h2 = loadSimParams(simNames[i])
      shortName = strmid(simNames[i],7,6)
      
      ; map to pS
      pS = mapInputParams(h2.c,h2.v200,h2.f_d,f_z)
      
      ; theoretical
      sigmaR  = pS.M_star * 1e10 / (2*!pi*(h2.h*1000)^2.0) * exp(-1.0 * R / h2.h) ;msun/pc^2
      
      ; plot theoretical
      fsc_plot,R/h2.h,sigmaR,line=0,color=fsc_color(units.colors[i]),/overplot
        
      ; plot legend
      fsc_plot,[legendXStart,legendXStart+legendXLength],$
               [legendYStart-legendYSpace*i,legendYStart-legendYSpace*i],$
               color=fsc_color(units.colors[i]),line=line,thick=2.0,/overplot
      fsc_text,legendXStart+legendXLength*1.2,legendYStart-(legendYSpace/5)-legendYSpace*i,$
               shortName,charsize=1.0,alignment=0.0,/data        
      
    endfor
        
  end_PS
  
  stop
end

; debug: check our model agreement with MakeGalaxyDisk outputs

pro matchActual

  workingPath = "c:\zStuff\IDL\Default\"
  filePathQ   = 'c:\zStuff\IDL\Default\10mil.curve.txt'
  filePathVD  = 'c:\zStuff\IDL\Default\10mil.deldisp.dat'
  
  ; input params
  c    = 9.0
  v200 = 160.0 ; km/s (code units)
  f_d  = 0.02
  f_z  = 0.1
  RDF  = 1.0
  
  ; radial calc config
  minR = 0.01    ; kpc (code units)
  maxR = 15.0    ; kpc (code units)
  resR = 100.0
  R    = findgen(resR)/resR * (maxR-minR) + minR
  Rlog = 2.0^R / max(2.0^R) * (maxR-minR)  
  
  ; set G to code units
  units = getUnits()
  G = units.G
  
  ; theoretical Q(r), kappa2, veldispR
  pS = mapInputParams(c,v200,f_d,f_z)
  
  ; NOTE: MakeNewDisk uses zeta(z)=(2/(exp(z/z0)+exp(-z/z0)))^2=sech^2(z/z0)
  ;  ---> MakeNewDisk z0 = 2z0 of our zeta(z)=sech^2(z/2z0)
  
  theoryQ        = func_Q_ExpHern(G,Rlog,pS.M_DM,pS.M_star,pS.h,pS.a,pS.z0,RDF,/thick)  
  theoryKappa2   = func_kappa2_ExpHern(G,Rlog,pS.M_DM,pS.M_star,pS.h,pS.a)
  theoryVelDispR = func_velDispR_ExpHern(G,Rlog,pS.M_DM,pS.M_star,pS.h,pS.a,pS.z0,/thick)

  ; load actual Q(r), kappa2, velDisp
  h1 = loadMDGCurve(filePathQ,vcs,qs)
  h2 = loadVelDispDump(filePathVD,velDisp)
  ;velDispR(r) = velDisp.velDispRz_Disk[0,*]
  
  ; make Q(r) with theory kappa2 + velDisp
  sigmaR      = ps.M_star / (2*!pi*pS.h^2.0) * exp(-1.0 * velDisp.list_R / pS.h)
  tKappa2     = func_kappa2_ExpHern(G,velDisp.list_R,pS.M_DM,pS.M_star,pS.h,pS.a)
  
  halfQ  = 1.0 * sqrt(velDisp.velDispRz_Disk[0,*]) * sqrt(tKappa2) / (3.36*G*sigmaR)
  
  ; plotting
  start_PS, workingPath+"matchActual_kappa2.eps"
  
    fsc_plot,[0],[0],ytitle="Epicyclic Frequency",xtitle="",charsize=1.2, $
             xrange=[0.1,15],/xs,/xlog,/ylog,yrange=[2e2,2e5],/ys,/nodata,position=[0.1,0.3,0.9,0.9], $
             xtickname=replicate(' ',10)
  
    fsc_plot,Rlog,theoryKappa2,line=0,/overplot
    fsc_plot,velDisp.list_R,velDisp.epi_kappa2,line=1,/overplot
    
    fsc_text,1.2,4e4,"Theoretical",charsize=1.2,alignment=0.5,/data
    fsc_text,0.7,1e4,"MakeDiskGalaxy",charsize=1.2,alignment=0.5,/data
    
    fsc_plot,[0],[0],ytitle="Fractional Residual",xtitle="Radius [kpc]",charsize=1.2,$
             xrange=[0.1,15],/xs,/xlog,xticklen=0.05, $
             yrange=[1e-5,3e-2],/ys,/ylog,/nodata,position=[0.1,0.1,0.9,0.3],/noerase
    
    ; residuals
    res = abs(theoryKappa2 - interpol(velDisp.epi_kappa2,velDisp.list_R,Rlog))
    res /= velDisp.epi_kappa2 ;norm
    
    fsc_plot,Rlog,res,/overplot       
    
  end_PS
  
  start_PS, workingPath+"matchActual_velDispR.eps"
  
    fsc_plot,[0],[0],ytitle="Radial Velocity Dispersion ["+textoidl("km^2/s^2")+"]",xtitle="",charsize=1.2, $
             xrange=[0.1,15],/xs,/xlog,yrange=[0,5e3],/ys,/nodata,position=[0.1,0.3,0.9,0.9], $
             xtickname=replicate(' ',10)
  
    fsc_plot,Rlog,theoryVelDispR,line=0,/overplot ;thick
    fsc_plot,velDisp.list_R,velDisp.velDispRz_Disk[0,*],line=1,/overplot
    
    fsc_text,0.6,2000,"Theoretical",charsize=1.2,alignment=0.5,/data
    fsc_text,0.4,900,"MakeDiskGalaxy",charsize=1.2,alignment=0.5,/data
    
    fsc_plot,[0],[0],ytitle="Fractional Residual",xtitle="Radius [kpc]",charsize=1.2,$
             xrange=[0.1,15],/xs,/xlog,$
             ytickname=textoidl(['10^{-3}','10^{-2}','10^{-1}','10^0']),xticklen=0.05, $
             yrange=[0.001,3.0],/ys,/ylog,/nodata,position=[0.1,0.1,0.9,0.3],/noerase
    
    ; residuals
    res = theoryVelDispR - interpol(velDisp.velDispRz_Disk[0,*],velDisp.list_R,Rlog)
    res /= velDisp.velDispRz_Disk[0,*] ;norm
    
    fsc_plot,Rlog,res,/overplot
    
  end_PS
  
  start_PS, workingPath+"matchActual_Q.eps"
            
    fsc_plot,[0],[0],ytitle="Toomre's Stellar Q(r) Parameter",xtitle="",charsize=1.2, $
             xrange=[0.1,15],/xs,/xlog,yrange=[0.5,5.0],/ys,/nodata,position=[0.1,0.3,0.9,0.9], $
             xtickname=replicate(' ',10)
             
    fsc_plot,Rlog,theoryQ,line=0,/overplot ;thick    
    ;fsc_plot,velDisp.list_R,halfQ,line=2,color=fsc_color('orange'),/overplot ;theory kappa2 + data veldispR
    fsc_plot,qs.list_R,qs.Q,line=1,/overplot
    
    fsc_text,0.5,2.5,"Theoretical",charsize=1.2,alignment=0.5,/data
    fsc_text,0.4,1.3,"MakeDiskGalaxy",charsize=1.2,alignment=0.5,/data
    
    fsc_plot,[0],[0],ytitle="Fractional Residual",xtitle="Radius [kpc]",charsize=1.2,$
             xrange=[0.1,15],/xs,/xlog,$
             ytickname=textoidl(['10^{-3}','10^{-2}','10^{-1}','10^0']),xticklen=0.05, $
             yrange=[0.001,3.0],/ys,/ylog,/nodata,position=[0.1,0.1,0.9,0.3],/noerase
    
    ; residuals
    res = abs(theoryQ - interpol(qs.Q,qs.list_R,Rlog))
    res /= qs.Q ;norm
    
    fsc_plot,Rlog,res,/overplot    
    
  end_PS
  
end

; model making
; ------------
; for LC1-6: c=[6.0,12.0], v200=[140.0,220.0],f_d=[0.017,0.023],RDF=[1.0,1.5]
;            resP=[61,81,7,11] ;[13,17,7,6] ;[25,33,7,11] ;[61,81,7,11]

; inputPlane()
pro inputPlane

  workingPath = "/n/home07/dnelson/spirals/lambdacrit/"

  ; set G to code units
  units = getUnits()
  G = units.G

  ; how many scale lengths to find lambda crit at?
  numSLs = 2

  ; fixed inputs
  f_z = 0.1 ; DiskHeight (fraction of H)

  ; inputs range (code units)
  c    = [4.0,12.0]    ; concentration parameter
  v200 = [140.0,280.0] ; km/s (code units)
  f_d  = [0.017,0.025] ; mass of disk (fraction of halo mass)
  RDF  = [1.0,1.5]     ; RadialDispersionFactor (mult factor sigma_r^2 to sigma_z^2)
  
  origParams  = [9.0,160.0,0.02,1.0]
  
  ; radial calc config
  minR = 0.01    ; kpc (code units)
  maxR = 15.0    ; kpc (code units)
  resR = 100.0
  R    = findgen(resR)/resR * (maxR-minR) + minR

  ; generate curves over variable parametera
  resP = [81,141,9,11]
  
  c_vals    = findgen(resP[0])/(resP[0]-1) * (c[1] - c[0]) + c[0]
  v200_vals = findgen(resP[1])/(resP[1]-1) * (v200[1] - v200[0]) + v200[0]
  f_d_vals  = findgen(resP[2])/(resP[2]-1) * (f_d[1] - f_d[0]) + f_d[0]
  RDF_vals  = findgen(resP[3])/(resP[3]-1) * (RDF[1] - RDF[0]) + RDF[0]
  
  print,'c_vals: ',c_vals
  print,'v200_vals: ',v200_vals
  print,'f_d_vals: ',f_d_vals
  print,'RDF_vals: ',RDF_vals
  
  minQs       = fltarr(resP[0],resP[1],resP[2],resP[3]) ;4 params
  lambdaCrits = fltarr(resP[0],resP[1],resP[2],resP[3]) ;4 params
  
  ; loop over parameter space and calculate minimum Q and lambda crit values
  for i=0,resP[0]-1 do begin
    for j=0,resP[1]-1 do begin
      for k=0,resP[2]-1 do begin
        for l=0,resP[3]-1 do begin
    
          ; transform input parameters to model parameters
          pS = mapInputParams(c_vals[i],v200_vals[j],f_d_vals[k],f_z)
          
          tempQ = func_Q_ExpHern(G,R,pS.M_DM,pS.M_star,pS.h,pS.a,pS.z0,RDF_vals[l],/thick)
          minQs[i,j,k,l] = min( tempQ )
          
          tempLC = func_LambdaCrit_ExpHern(G,R,pS.M_DM,pS.M_star,pS.h,pS.a)
          diffLC = abs(tempLC - numSLs * pS.h)
          w = where(diffLC eq min(diffLC),count)
          
          if (count) then $
            lambdaCrits[i,j,k,l] = tempLC[w]
            
        endfor
      endfor
    endfor
  endfor
  
  ; pick 5 fiducial models at Q~targetQ
  targetQ     = 1.2
  targetCL    = [2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5]
  
  searchFac   = 3.0 ; distance expansion about target vals for other candidates (1.0=disable)
  weightFac   = 0.2 ; parameter weighting in distance search, [0,1]
  
  fModels = fltarr(n_elements(targetCL),6) ;[Q,CL,4params]

  for i=0, n_elements(targetCL)-1 do begin
  
    ; choose model closest to target Q and lambda_crit
    dist2DFromCLQ = (targetCL[i] - lambdaCrits)^2.0 + (minQs - targetQ)^2.0
    
    w = where(dist2DFromCLQ le min(dist2DFromCLQ)*searchFac^2.0, count)
    
    ;debug:
    ;wT = where(dist2DFromCLQ eq min(dist2DFromCLQ))
    ;t4 = array_indices(minQs,wT)
    ;print,'unweighted: ',minQs[wT],lambdaCrits[wT],$
    ;      c_vals[t4[0]],v200_vals[t4[1]],f_d_vals[t4[2]],f_z_vals[t4[3]]
    print,'weighted search among: ',count
    
    ; include weighting on distance from mean of each parameter to get the most reasonable
    effWeights     = fltarr(4,count)
    effWeightDists = fltarr(count)

    for j=0l,count-1 do begin
      ;normalized (fractional) error in each parameter
      inds4d = array_indices(minQs,w[j])
      effWeights[0,j] = (origParams[0]-c_vals[inds4d[0]])^2.0 / origParams[0]^2.0
      effWeights[1,j] = (origParams[1]-v200_vals[inds4d[1]])^2.0 / origParams[1]^2.0
      effWeights[2,j] = (origParams[2]-f_d_vals[inds4d[2]])^2.0 / origParams[2]^2.0
      effWeights[3,j] = (origParams[3]-RDF_vals[inds4d[3]])^2.0 / origParams[3]^2.0
      
      effWeightDists[j] = weightFac * total(effWeights[*,j])      
    endfor
  
    ; choose new effectively closest model
    effDist6D = sqrt(dist2DFromCLQ[w] + effWeightDists)  
    
    w2 = where(effDist6D eq min(effDist6D), count2)
    inds4d = array_indices(minQs,w[w2])
    
    ;debug:
    ;print,'weighted: ',minQs[w[w2]],lambdaCrits[w[w2]],$
    ;      c_vals[inds4d[0]],v200_vals[inds4d[1]],f_d_vals[inds4d[2]],f_z_vals[inds4d[3]]

    fModels[i,0:1]  = [minQs[w[w2]],lambdaCrits[w[w2]]]
    fModels[i,2:5]  = [c_vals[inds4d[0]],v200_vals[inds4d[1]],f_d_vals[inds4d[2]],RDF_vals[inds4d[3]]]

  endfor

  ; plotting
  PS_Start, FILENAME=workingPath+"inputPlane.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=5.0, /inches
            
    yrange = [0.5,8.0]
            
    fsc_plot,[0],[0],xtitle="Minimum of Q within 15 kpc",$
             ytitle="Critical Wavelength at "+str(numSLs)+" Scale Lengths [kpc]",charsize=1.2, $
             xrange=[0.5,2.0],yrange=yrange,/ys,/nodata 
           
    fsc_plot,[1.1,1.1],yrange,line=0,color=fsc_color('gray'),/overplot
    fsc_plot,[1.3,1.3],yrange,line=0,color=fsc_color('gray'),/overplot           
    
    colors1 = ['blue','brown','orange','cyan','deep pink','red','orchid','green']
    
    ; plot as lines (numP-1 dim)         
    for i=0,resP[0]-1 do begin
      for j=0,resP[1]-1 do begin
        for k=0,resP[2]-1 do begin
          fsc_plot,minQs[i,j,k,*],lambdaCrits[i,j,k,*],psym=14,symsize=0.3,$
                   color=fsc_color(colors1[k mod n_elements(colors1)]),/overplot
        endfor
      endfor
    endfor
    
    ; overplot fModels
    fsc_plot,fModels[*,0],fModels[*,1],psym=14,symsize=0.6,color=fsc_color('black'),/overplot
    
  PS_End 

  ; dump parameters of fModels for table
  for i=0,n_elements(targetCL)-1 do begin
    print,string(format='(i," & ",f4.2," & ",f4.2," & ",f5.2," & ",f5.1," & ",f5.3," & ",f5.3," \\")',$
                 i,fModels[i,*])
  endfor

end

; radiallyBinStars():
function radiallyBinStars, snapPath, iSnap
  
  ; config
  rMin = 0.0
  rMax = 30.0
  rRes = 100
  
  ; bins
  rPts = findgen(rRes)/rRes * (rMax - rMin) + rMin
  rLog = 1.5^rPts / max(1.5^rPts) * (rMax-rMin)
  
  rS = {rCenters:fltarr(rRes-1),$
        areas:fltarr(rRes-1),$
        counts:lonarr(rRes-1),$
        massPerStar:0.0} 

  ; load ICs
  h = loadSnapshot(snapPath,iSnap,s_pos,s_vel,s_id,c_pos,c_vel,c_id)

  r = reform(sqrt(s_pos[0,*]^2.0 + s_pos[1,*]^2.0))

  for i=0,rRes-2 do begin
    w = where(r ge rPts[i] and r lt rPts[i+1], count)
    
    rS.rCenters[i] = (rPts[i]+rPts[i+1])/2.0
    rS.areas[i]    = !pi*rPts[i+1]^2.0 - !pi*rPts[i]^2.0
    rS.counts[i]   = count
  endfor

  rS.massPerStar = h.massTable[2]
  return,rS
end

; compare6(): plot comparison of disk and DM profiles for our "6" proposed models
pro compare6

  units = getUnits()

  ; config
  snapPath = '/n/home07/dnelson/spirals/lambdacrit/ICs/'
  workingPath = '/n/home07/dnelson/spirals/lambdacrit/'

  fileNames = ['MW_disk_c9_10mil.dat',$
               'LC_disk_1.dat','LC_disk_2.dat','LC_disk_3.dat',$
               'LC_disk_4.dat','LC_disk_5.dat','LC_disk_6.dat']
               
  R_H  = [29.7754, 26.1946, 28.7315, 36.2805, 37.6057, 44.5400, 27.2245]
  M_DM = [95.2401, 72.8103, 105.022, 111.042, 138.107, 138.107, 160.748]
               
  colors = ['black','blue','brown','orange','cyan','deep pink','red','orchid','green']
  
  ; start disk profiles plot
  PS_Start, FILENAME=workingPath+"compare6.disk.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=7.0, /inches
    
    fsc_plot,[0],[0],xtitle="radius [kpc]",ytitle="disk surface density [Msun/pc^2]",$
             charsize=1.5,/nodata,xrange=[1.0,20],/xs,yrange=[0,300],/ys,/xlog
  
  for i=0,n_elements(fileNames)-1 do begin
    ; load IC and bin radially
    rS = radiallyBinStars(snapPath+fileNames[i],'none')
    
    ; convert counts to surface density (Msun/pc^2)
    areas_pc2     = rS.areas * 1e6
    mass_msun     = rS.massPerStar * units.UnitMass_in_g / units.Msun_in_g
    
    dens = rS.counts * mass_msun / areas_pc2
    
    ; add to disk profiles plot
    fsc_plot,rS.rCenters,dens,line=0,color=fsc_color(colors[i]),/overplot
    
    ; add to legend
    fsc_plot,[10.0,12.0],[290-10*i,290-10*i],line=0,color=fsc_color(colors[i]),/overplot
    fsc_text,13.0,287-10*i,str(i),charsize=1.2,alignment=0.5
    
    ; debug:
    print,'i mass_per_star dens_minmax = ',i,rS.massPerStar,minmax(dens)
  endfor

  ; end disk profiles plot
  PS_End

  ; halo profiles plot
  PS_Start, FILENAME=workingPath+"compare6.halo.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=7.0, /inches
            
    fsc_plot,[0],[0],xtitle="radius [kpc]",ytitle="halo density [10^10 Msun/kpc^3]",$
             charsize=1.5,/nodata,xrange=[1.0,20],/xs,yrange=[0,0.03],/ys,/xlog

  for i=0,n_elements(fileNames)-1 do begin
    ; r points
    rMin = 0.1
    rMax = 30.0
    rRes = 100
  
    rPts = findgen(rRes)/rRes * (rMax - rMin) + rMin
  
    ; hernquist profile
    dens = M_DM[i] / (2*!pi) * R_H[i] / (rPts*(rPts+R_H[i])^3.0)
    print,minmax(dens)
    
    ; add to plot
    fsc_plot,rPts,dens,line=0,color=fsc_color(colors[i]),/overplot

  endfor

  PS_End

  stop
end

; epicycles
; ---------

pro plotEpi

  ; config
  M  = 1
  G  = 1
  Lz = 0.1
  Rg = 4.0
  
  A      = 1.0
  Psi    = 0.0
  Theta0 = 0.0

  ; epicyclic frequency
  ratios = [1.0,1.01,1.1,1.26,1.4,1.5,1.74,1.75,1.99]
  nRatios = n_elements(ratios)
  
  omegaG = Lz / Rg^2.0
  ;kappa  = -2.0 * G * M / Rg^3.0 + 3.0 * Lz^2.0 / Rg^4.0
  kappa = omegaG*ratios ;kepler potential
  
  P_G     = 2*!pi/omegaG
  P_kappa = 2*!pi/kappa
  
  ; time
  tRes = 2000
  tMin = 0.0
  tMax = P_G*10
  tPts = findgen(tRes)/tRes * (tMax-tMin) + tMin
    
  x = fltarr(tRes,nRatios)
  y = fltarr(tRes,nRatios)
  rTheta_x = fltarr(tRes,nRatios)
  rTheta_y = fltarr(tRes,nRatios)
  
  for i=0,nRatios-1 do begin
    ; x,y cartesian displacements
    x[*,i] = A * cos(kappa[i]*tPts + Psi)
    y[*,i] = -2.0 * omegaG * A / kappa[i] * sin(kappa[i]*tPts + Psi)
  
    ; r,theta global cylindrical
    r     = x[*,i] + Rg
    theta = (y[*,i]/Rg) + omegaG*tPts + Theta0
    
    ; convert to x,y for plotting
    rTheta_x[*,i]   = r * cos(theta)
    rTheta_y[*,i]   = r * sin(theta)
  endfor
  
  ; plot x,y
  PS_Start, FILENAME="~/spirals/epiXY.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=7.0, /inches
    
    xyR = [-2.5,2.5]
    fsc_plot,[0],[0],xtitle="x",ytitle="y",charsize=1.5,/nodata,xrange=xyR,/xs,yrange=xyR,/ys
    fsc_plot,x[*,0],y[*,0],/overplot

  PS_End
  
  ; plot r,theta
  PS_Start, FILENAME="~/spirals/epiRTheta.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=8.0, /inches

    !p.multi = [0,3,3]
    xm = !x.margin
    ym = !y.margin
    !x.margin = [0,0]
    !y.margin = [0,0]
    
    xyR = [-6,6]
    
    for i=0,nRatios-1 do begin
      fsc_plot,[0],[0],charsize=1.5,/nodata,xrange=xyR,/xs,yrange=xyR,/ys, $
        xtickname=replicate(' ',10),ytickname=replicate(' ',10)
      fsc_plot,rTheta_x[*,i],rTheta_y[*,i],/overplot
      
      ; overplot circular orbit at Rg
      Rg_x = Rg * cos(findgen(100)/100 * 2 * !pi)
      Rg_y = Rg * sin(findgen(100)/100 * 2 * !pi)
      
      fsc_plot,Rg_x,Rg_y,/overplot,line=1
      
      ; overplot kappa/omegaG ratio
      fsc_text,3.5,-5,string(ratios[i],format='(f4.2)'),/data,charsize=1.5
    endfor
    
    !x.margin = xm
    !y.margin = ym
    !p.multi = [0]

  PS_End
  
  stop
end

; findICBins():
; find equal stellar number bins given ICs and check Poisson noise level

pro findICBins

  ;config
  ;ICsPath = '/n/home07/dnelson/spirals/perturbers/ICs/MW_disk_and_mc_c9_pro_test.dat'
  ICsPath = '/n/home07/dnelson/spirals/perturbers/ICs/MW_disk_and_mc_c9_10mil_pro.dat'

  ;100,000 stars
  ;nRadialBins  = 25
  ;nAngularBins = 18
  
  ;10 million stars
  nRadialBins  = 50
  nAngularBins = 90
  
  minMaxRadius = [0.0, 50.0]
  
  ;make angular bins
  angularBins = fltarr(nAngularBins+1)
  for m=0,nAngularBins do begin
    angularBins[m] = 360.0 / nAngularBins * m;
  endfor
  
  ;find poisson
  radialBins = findPoissonBins(ICsPath, nRadialBins, nAngularBins, minMaxRadius)
  
  ;check poisson
  checkPoissonBins, ICsPath, radialBins, angularBins
end
