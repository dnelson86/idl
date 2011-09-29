; TI_analysis.pro
; thermal instability related / Arepo test results
; dnelson 13.march.2011

@helper
@arepoLoad
@arepoVis2D
@TI_plot

;+
; NAME:
;        CUBIC
;
; PURPOSE:
;        Solve for the roots of a cubic polynomial OR
;        find the value(s) of a cubic polynomial (INVERSE)
;
; CATEGORY:
;        Math.
;
; CALLING SEQUENCE:
;
;        Result = CUBIC( X, Coeff )
;
; INPUTS:
;        X:   For the forward transform, it is the VALUE of the cubic
;             polynomial. If the INVERSE keyword is set then it is
;             the ARGUMENT of the cubic polynomial.
;
;        Coeff:    Array of coefficients of the cubic polynomial, fltarr(4).
;
; KEYWORD PARAMETERS:
;
;        INVERSE:  Set this keyword to determine the value of the cubic
;             polynomial, (0=Default).
;
; OUTPUTS:
;        Returns the ARGUMENT of a cubic polynomial for the forward
;        transform or the VALUE of a cubic polynomial for the inverse
;        transform (if the INVERSE keyword is set).
;
; MODIFICATION HISTORY:
;        Written by:    Han Wen, October 1994.
;
;-
function SIGN, A, B

         NP   = N_PARAMS()
         npt  = n_elements( A )

         if NP eq 1 then begin
              n    = A
              m    = replicate( 1.d0,npt )
         endif else begin
              n    = B
              m    = A
         endelse

         if npt eq 1 then begin
              if n ge 0 then $
                   return, m(0) $
              else           $
                   return,-m(0)
         endif else begin
              signs   = A*0
              here_ge = WHERE( n ge 0, nge )
              here_lt = WHERE( n lt 0, nlt )
              if nge gt 0 then signs( here_ge ) =  m( here_ge )
              if nlt gt 0 then signs( here_lt ) = -m( here_lt )
              return, signs
         endelse

end

function CUBIC, X, Coeff, INVERSE=Inverse

         a    = Coeff(0)
         b    = Coeff(1)
         c    = Coeff(2)
         d    = Coeff(3)

         if keyword_set( INVERSE ) then begin

              h=X
              t=a+b*h+c*h^2+d*h^3

              return, t

         endif else begin

              t=X

              s1= ATAN(SQRT(3)*(27*a*d^2-9*b*c*d+2*c^3-27*d^2*t)/(9*d*(c^2 $
                       -3*b*d)^(1.5)*SQRT((27*a^2*d^2-2*a*(9*b*c*d $
                       -2*c^3+27*d^2*t)+4*b^3*d-b^2*c^2+18*b*c*d*t $
                       -t*(4*c^3-27*d^2*t))/(3*b*d-c^2)^3)) $
                      )/3

              r1= SQRT(3)*SQRT(c^2-3*b*d)*COS(s1)/(3*d)
              r2= SQRT(c^2-3*b*d)*SIN(s1)/(3*d)
              r3= -c/(3*d)

              h1= 2*SQRT(c^2-3*b*d)*SIN(s1)/(3*ABS(d))+r3
              h2=-SIGN(d)*(r1 + r2) + r3
              h3= SIGN(d)*(r1 - r2) + r3

              h = [[h1],[h2],[h3]]
              h = REFORM(h)

              return, h

         endelse
end


function cuberoot,cc
;+
; NAME:
; CUBEROOT
; PURPOSE:
; Real roots of a cubic equation. Complex roots set to -1.0e30.
; Called by HALFAGAUSS.
; CALLING SEQUENCE:
; roots = CUBEROOT(cc)
; INPUT:
; cc = 4 element vector, giving coefficients of cubic polynomial, 
;   [c0,c1,c2,c3], where one seeks the roots of
;   c0 + c1*x + c2*x^2 + c3*x^3
; OUTPUT:
; Function returns a vector of 3 roots. If only one root is real, then 
; that becomes the 1st element.
; EXAMPLE:
; Find the roots of the equation
;   3.2 + 4.4*x -0.5*x^2 -x^3 = 0
; IDL> x = cuberoot([3.2,4.4,-0.5,-1])
; 
; will return a 3 element vector with the real roots 
;   -1.9228, 2.1846, -0.7618
; REVISION HISTORY:
; Henry Freudenreich, 1995
;-

  on_error,2
  unreal=-1.0e30
  
  a1=cc(2)/cc(3)
  a2=cc(1)/cc(3)
  a3=cc(0)/cc(3)
  
  q=(a1^2-3.*a2)/9.            &  r=(2.*a1^3-9.*a1*a2+27.*a3)/54.
  
  if r^2 lt q^3 then begin
  ;  3 real roots.
     theta=acos(r/q^1.5)
     x1=-2.*sqrt(q)*cos(theta/3.)-a1/3.
     x2=-2.*sqrt(q)*cos((theta+6.28319)/3.)-a1/3.
     x3=-2.*sqrt(q)*cos((theta-6.28319)/3.)-a1/3.
  endif else begin
  ;  Get the one real root:
     a=-r/abs(r) * (abs(r)+sqrt(r^2-q^3))^.33333
     if a eq 0. then b=0. else b=q/a
     x1=(a+b)-a1/3.
     x2=unreal
     x3=unreal
  endelse 
  
  roots=[x1,x2,x3]
  return,roots
end




; plotEquilCurve(): plot GK10 P-rho equil curve using all data & non-constant heating w/ SS02 curve

pro plotEquilCurve

  mass_hydrogen = double(1.67e-24) ;g

  ; load GK10 (all)
  fileName = 'C:\zStuff\IDL.work\Default\cooling\gk10\dl.Z=1.0_UV=1.res'
  headerLines = 0
  ptStruct = {n_b:0.0,  T:0.0, col_b:0.0, Z:0.0, f_HI:0.0, f_HII:0.0, f_H2:0.0, U_MW:0.0, $
              D_MW:0.0, rho_over:0.0, fc:0.0, fh:0.0}
              
  gk = loadCSV(headerLines,fileName,ptStruct)
  
  ; load GK10 (tabulated)
  fileName = 'C:\zStuff\IDL.work\Default\cooling\gk10\dl.Z=1.0_UV=1.res.txt'
  headerLines = 1 ;nLines
  ptStruct = {logT:0.0, logLambda:0.0, logHeat:0.0, ne_to_nh:0.0}
              
  gkTab = loadCSV(headerLines,fileName,ptStruct)
  lambdaGK = 10.0^(gkTab.logLambda) / mass_hydrogen^2.0
  
  ; load SS02
  
  tMin = 1.0 ;log K
  tMax = 5.0 ;log K
  tRes = 100
  
  tPts = findgen(tRes)/tRes * (tMax - tMin) + tMin
  tPts = 10.0^(tPts)
  
  lambdaPts = LambdaArrSS02(tPts)
  
  ; calculate GK10 (all) P-rho equil curve
  
  equilRhoGK = gk.n_b * 10.0^(gk.fh) / 10.0^(gk.fc) ;cm^(-3)
  
  P_k_GK = gk.T * equilRhoGK ;K cm^(-3)
  
  ; calculate the SS02 P-rho equilibrium cooling curve
  
  fHeat = 0.015 ;constant heating function, erg/s/g
  
  equilRhoPts = fHeat / lambdaPts ; 0 = fHL = equilRhoPts * Lambda - fHeat ;erg/s/g
  equilRhoPts /= mass_hydrogen ;cm^(-3)
  
  P_k = tPts * equilRhoPts ;K cm^(-3)  
  
  ; calculate GK10 (tabulated) P-rho equil curve
  Tmax = 7.1 ; beyond this we are using GS07 cooling functions with no heating func
  
  w = where(gkTab.logT le Tmax)
  
  equilRhoGKTab = 10.0^(gkTab[w].logHeat) / (10.0^(gkTab[w].logLambda)) ;cm^(-3)
  
  ;equilRhoGKTab = 100*fHeat / (10.0^(gkTab[w].logLambda)) ;cm^(-3) g^(-1)
  ;equilRhoGKtab *= mass_hydrogen ;cm^(-3)
  
  P_k_GKTab = 10.0^(gkTab[w].logT) * equilRhoGKTab ;K cm^(-3)

  ; plot P/k vs rho
  
  ;set_plot,'PS'
  ;device,filename='ss02gk10.eps',bits_per_pixel=8,color=1,/landscape,/encapsulated,/decomposed  
  
    plot,[0],[0],xtitle="Log "+textoidl('\rho')+" [cm"+textoidl('^{-3}')+"]", $
         ytitle="Log P/k [K cm"+textoidl('^{-3}')+"]", $
         thick=1.1,charsize=1.5,/xstyle,/ystyle,/nodata,xrange=[-4,3], $
         yrange=[2,6]
    
    oplot,alog10(equilRhoGK),alog10(P_k_GK),psym=3,color=fsc_color('orange')
    oplot,alog10(equilRhoGKTab),alog10(P_k_GKtab),psym=0,thick=2,color=fsc_color('red')
    oplot,alog10(equilRhoPts),alog10(P_k),psym=0,thick=2,color=fsc_color('green') 

  ;device,/close_file
  ;set_plot,'WIN'

  stop

end

; plotGrowthRateGK10(): TI growth rate using GK10 cooling function

pro plotGrowthRateGK10

  ; init
  Rgas          = 8.3314e7 ;gas constant, erg/K/mol
  mass_hydrogen = 1.67e-24 ;g
  boltzmann     = 1.38e-16 ;erg/K

  ; config
  P_k0  = 2400.0           ;K/cm^3
  n_0   = 1.0              ;cm^(-3)
  
  ;fHeat = 0.015   ;constant heating function, erg/s/g
  
  kappa = 0.0     ;conduction, erg/cm/K/s
  gamma = 5.0/3.0 ;specific heat ratio
  mu    = 1.22    ;mean molecular weight  
  
  T_0   = P_k0 / n_0                  ;K
  P_0   = P_k0 * boltzmann            ;ba
  rho_0 = n_0 * (mu * mass_hydrogen)  ;g/cm^3
  
  csnd  = (gamma * P_0 / rho_0)^(0.5) ;cm/s
  
  ; load GK10
  fileNameGK10 = 'C:\zStuff\IDL.work\Default\cooling\gk10\dl.Z=1.0_UV=1.res.txt'
  headerLines = 1 ;nLines
  ptStruct = {logT:0.0, logLambda:0.0, logHeat:0.0, ne_to_nh:0.0}
              
  gk = loadCSV(headerLines,fileNameGK10,ptStruct)
  
  ; compute partial derivatives of heatloss function (3pt Lagrangian method)
  pdXSpread = [0.98,0.99,1.0,1.01,1.02]
  
  fHL_X = double(T_0) * pdXSpread
  fHL_L = 10.0^(interpol(gk.logLambda, gk.logT, alog10(fHL_X)))
  fHL_H = 10.0^(interpol(gk.logHeat,   gk.logT, alog10(fHL_X)))
  fHL_Y = (fHL_L * n_0 - fHL_H) / mass_hydrogen
  
  pdHL_T   = ( deriv(fHL_X,fHL_Y) )[2]

  fHL_X = double(rho_0) * pdXSpread
  fHL_L = 10.0^(interpol(gk.logLambda, gk.logT, alog10(T_0)))
  fHL_H = 10.0^(interpol(gk.logHeat, gk.logT, alog10(T_0)))
  fHL_Y = ( (fHL_X / mass_hydrogen) * fHL_L - fHL_H ) / mass_hydrogen
  
  pdHL_rho = ( deriv(fHL_X,fHL_Y) )[2] 

  ; compute wavenumbers [1/cm]
  k_T   = (mu * (gamma - 1.0) * pdHL_T) / (Rgas * csnd)
  k_rho = (mu * (gamma - 1.0) * rho_0 * pdHL_rho) / (Rgas * csnd * T_0)
  k_K   = (Rgas * csnd * rho_0) / (mu * (gamma - 1.0) * kappa)

  ;k_rho = 2.0 * k_T   ;a = 1/2
  ;k_K   = 100 * k_rho ;beta=0.1 (coefficient=1/beta)

  ;compute kappa=0 k->inf limit of growth rate
  n_1 = -1.0 * (gamma - 1) * mu * (T_0 * pdHL_T - rho_0 * pdHL_rho) / (gamma * Rgas * T_0)
  n_1_norm = n_1 / (k_rho * csnd)

  ; wavenumber points
  kMin = 0.01*k_rho ;norm
  kMax = 6*k_rho   ;norm
  kRes = 100
  
  kPts = findgen(kRes)/kRes * (kMax - kMin) + kMin
  
  ; solve cubic characteristic equation
  roots = fltarr(3,n_elements(kPts))
  
  c0 = csnd^3.0 * kPts^2.0 / gamma * (k_T - k_rho + kPts^2.0/k_K)
  c1 = csnd^2.0 * kPts^2.0
  c2 = csnd * (k_T + kPts^2.0/k_K)
  c3 = 1.0

  for i=0,n_elements(kPts)-1 do begin
    roots[*,i] = cuberoot([c0[i],c1[i],c2[i],c3])
    ;roots[*,i] = cubic(0.0,[c0[i],c1[i],c2[i],c3])
  endfor
  
  ; plot normalized
  plotX = kPts / k_rho
  
  plotY = fltarr(n_elements(kPts)) ;handle cuberoot() return
  for i=0,n_elements(kPts)-1 do begin
    if ( roots[0,i] > 0 ) then begin
      plotY[i] = roots[0,i] / (k_rho * csnd)
    endif
    if ( roots[0,i] < 0 and roots[1,i] > 0 ) then begin
      plotY[i] = roots[1,i] / (k_rho * csnd)
    endif
  endfor

    plot,[0],[0],xtitle="Normalized Wavenumber k / k"+textoidl('_\rho')+"", $
         ytitle="Normalized Growth Rate n / (k"+textoidl('_\rho')+" c"+textoidl('_s')+")", $
         thick=1.1,charsize=1.5,/xstyle,/ystyle,/nodata,xrange=minmax(plotX), $
         yrange=[0,max(plotY)*1.1]
    
    oplot,plotX,plotY,psym=3;,color=fsc_color('green')
    oplot,[minmax(plotX)],[n_1_norm,n_1_norm]*1.01,line=2
  
  stop
end

; eigenSS02(): compute eigenvalues and eigenvectors for SS02

pro eigenSS02

  ;---------- (common with plotGrowthRateSS02)

    ;init
    Rgas          = 8.3314e7 ;gas constant, erg/K/mol
    mass_hydrogen = 1.67e-24 ;g
    boltzmann     = 1.38e-16 ;erg/K
    cm_in_pc      = 3.086e18
  
    ; config
    P_k0  = 2400.0          ;K/cm^3
    n_0   = 1.0             ;cm^(-3)
    
    fHeat = 0.015   ;constant heating function, erg/s/g
    
    kappa = 0.0     ;conduction, erg/cm/K/s
    ;kappa = 7.48e6
    gamma = 5.0/3.0 ;specific heat ratio
    mu    = 1.22     ;mean molecular weight  
    
    T_0   = P_k0 / n_0                  ;K
    P_0   = P_k0 * boltzmann            ;ba
    rho_0 = n_0 * (mu * mass_hydrogen)  ;g/cm^3
    
    csnd  = (gamma * P_0 / rho_0)^(0.5) ;cm/s
  
    ; compute partial derivatives of heatloss function (3pt Lagrangian method)
    ;pdXSpread = [0.98,0.99,1.0,1.01,1.02]
    pdXSpread = [0.80,0.90,1.0,1.10,1.20]
    fHL_X = double(T_0) * pdXSpread
    fHL_Y = rho_0 * LambdaArrSS02(fHL_X) - fHeat
    
    pdHL_T   = ( deriv(fHL_X,fHL_Y) )[2]
    
    fHL_X = double(rho_0) * pdXSpread
    fHL_Y = fHL_X * LambdaArrSS02(T_0) - fHeat
    
    pdHL_rho = ( deriv(fHL_X,fHL_Y) )[2] 
    ;pdHL_rho = LambdaArrSS02(T_0) ;analytically, equivalent
  
    ; compute wavenumbers [1/cm]
    k_T   = (mu * (gamma - 1.0) * pdHL_T) / (Rgas * csnd)
    k_rho = (mu * (gamma - 1.0) * rho_0 * pdHL_rho) / (Rgas * csnd * T_0)
    k_K   = (Rgas * csnd * rho_0) / (mu * (gamma - 1.0) * kappa)
  
  ;--------
  
  ; choose k
  kvec = 1.0 * k_rho
  i_const = complex(0,1)
  
  ; solve cubic characteristic equation for n
  roots = fltarr(3)
  
  c0 = csnd^3.0 * kvec^2.0 / gamma * (k_T - k_rho + kvec^2.0/k_K)
  c1 = csnd^2.0 * kvec^2.0
  c2 = csnd * (k_T + kvec^2.0/k_K)
  c3 = 1.0

  roots = cuberoot([c0,c1,c2,c3])
  n = roots[0]
  
  ; create coefficient matrix (4x4)
  A = [[ n,                                                  0,                   0,                           i_const*rho_0 ], $
       [ 0,                                                  i_const*kvec^2,      0,                           rho_0*n       ], $
       [ rho_0*pdHL_rho - (gamma)/(gamma-1.0)*(P_0/rho_0)*n, (1.0/(gamma-1.0))*n, rho_0*pdHL_T + kappa*kvec^2, 0             ], $
       [ -1.0/rho_0,                                         1.0/P_0,             -1.0/T_0,                    0             ]]
  
  ; re-compute eigenvalues
  hes = la_elmhes(A, QZ, permute_result=PR, scale_result=SR)
  evals = la_hqr(hes, QZ, permute_result=PR)
  
  schur_T = hes
  schur_v = QZ
  
  ; compute eigenvectors
  evecs = la_eigenvec(schur_T, schur_v, /double, permute_result=PR, scale_result=SR)
  
    stop
  
end

; plotGrowthRateSS02(): TI growth rate using SS02 cooling function

pro plotGrowthRateSS02

  workingPath     = '/n/home07/dnelson/compTI.w0/'

  ;init
  Rgas          = 8.3314e7 ;gas constant, erg/K/mol
  mass_hydrogen = 1.67e-24 ;g
  boltzmann     = 1.38e-16 ;erg/K
  cm_in_pc      = 3.086e18

  ; config
  P_k0  = 2400.0          ;K/cm^3
  n_0   = 1.0             ;cm^(-3)
  
  fHeat = 0.015   ;constant heating function, erg/s/g ;0.015 previously
  
  kappa = [0.0,1.0e6,1.0e7,2.5e7,1.0e8]     ;conduction, erg/cm/K/s
  gamma = 1.667 ;specific heat ratio
  mu    = 1.22     ;mean molecular weight  
  
  T_0   = P_k0 / n_0                  ;K
  P_0   = P_k0 * boltzmann            ;ba (erg/cm^3)
  rho_0 = n_0 * (mu * mass_hydrogen)  ;g/cm^3
  
  csnd  = sqrt(gamma * P_0 / rho_0) ;cm/s

  ; compute partial derivatives of heatloss function (3pt Lagrangian method)
  ;pdXSpread = [0.98,0.99,1.0,1.01,1.02]
  pdXSpread = [0.80,0.90,1.0,1.10,1.20]
  fHL_X = double(T_0) * pdXSpread
  fHL_Y = rho_0 * LambdaArrSS02(fHL_X) - fHeat
  
  pdHL_T   = ( deriv(fHL_X,fHL_Y) )[2]

  fHL_X = double(rho_0) * pdXSpread
  fHL_Y = fHL_X * LambdaArrSS02(T_0) - fHeat
  
  pdHL_rho = ( deriv(fHL_X,fHL_Y) )[2] 
  ;pdHL_rho = LambdaArrSS02(T_0) ;analytically, equivalent

  ; compute wavenumbers [1/cm]
  k_T   = (mu * (gamma - 1.0) * pdHL_T) / (Rgas * csnd)
  k_rho = (mu * (gamma - 1.0) * rho_0 * pdHL_rho) / (Rgas * csnd * T_0)
  k_K   = (Rgas * csnd * rho_0) / (mu * (gamma - 1.0) * kappa)

  ;debug:
    print,'csnd: ',csnd
    print,'k_T: ',k_T
    print,'k_rho: ',k_rho
    print,'k_K: ',k_K
    ;k_rho = 2.0 * k_T    ;a = 1/2
    ;k_K   = 10.0 * k_rho ;beta=0.2
  
  ;compute kappa=0 k->inf limit of growth rate
  n_1 = -1.0 * (gamma - 1) * mu * (T_0 * pdHL_T - rho_0 * pdHL_rho) / (gamma * Rgas * T_0)
  n_1_norm = n_1 / (k_rho * csnd)

  ; wavenumber points
  kMin = 0.0
  kMax = 6.0
  kRes = 100
  
  kPts = findgen(kRes)/kRes * (kMax - kMin) + kMin
  kLog = 2.0^kPts / max(2.0^kPts) 
  kLog = (kLog - min(kLog)) * (kMax-kMin)  
  
  ; solve cubic characteristic equation
  plotY = fltarr(n_elements(kappa),kRes) ;handle cuberoot() return

  ; loop over conduction values
  for j=0,n_elements(kappa)-1 do begin
    roots = fltarr(3,kRes)
  
    ;convert range to cgs (with respect to k_rho)
    kPts = kLog * k_rho
    
    c0 = csnd^3.0 * kPts^2.0 / gamma * (k_T - k_rho + kPts^2.0/k_K[j])
    c1 = csnd^2.0 * kPts^2.0
    c2 = csnd * (k_T + kPts^2.0/k_K[j])
    c3 = 1.0
  
    for i=0,n_elements(kPts)-1 do begin
      roots[*,i] = cuberoot([c0[i],c1[i],c2[i],c3])
      ;roots[*,i] = cubic(0.0,[c0[i],c1[i],c2[i],c3])
  
      if ( roots[0,i] > 0 ) then begin
        plotY[j,i] = roots[0,i]              ;s^(-1)
      endif
      if ( roots[0,i] < 0 and roots[1,i] > 0 ) then begin
        plotY[j,i] = roots[1,i]              ;s^(-1)
      endif
    endfor
    
  endfor
  
  ; plot normalized
  plotX = kPts / k_rho
  ;plotX = 2.0 * !pi / kPts / cm_in_pc   ;pc
  plotY /= (k_rho * csnd)
  
  PS_Start, FILENAME=workingPath+"growthRate.SS02.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=5, /inches
              
    ;plot in cgs
    ;plot,[0],[0],xtitle="Wavelength "+textoidl('\lambda')+" [pc]", $
    ;     ytitle="Growth Rate n [s"+textoidl('^{-1}')+"]", $
    ;     thick=1.1,charsize=1.5,/xstyle,/xlog,/ystyle,/nodata,xrange=[0.1,max(plotX)], $
    ;     yrange=[0,max(plotY)*1.1]
    ;oplot,plotX,plotY,psym=3;,color=fsc_color('green') 
         
    ;plot normalized
    fsc_plot,[0],[0],xtitle="Normalized Wavenumber k / k"+textoidl('_\rho')+"", $
         ytitle="Normalized Growth Rate n / (k"+textoidl('_\rho')+" c"+textoidl('_s')+")", $
         thick=1.1,charsize=1.5,/xstyle,/ystyle,/nodata,xrange=[0,max(plotX)], $
         yrange=[0,max(plotY)*1.2]
    fsc_plot,plotX,plotY[0,*],line=0,/overplot;,color=fsc_color('green')
    
    for i=0,n_elements(kappa)-1 do begin
      fsc_plot,plotX,plotY[i,*],line=1,/overplot
    endfor

    fsc_plot,[minmax(plotX)],[n_1_norm,n_1_norm]*1.01,line=2,/overplot

  PS_End
  
end

; plotGrowthRateKI02(): TI growth rate using KI02 cooling function

pro plotGrowthRateKI02

  ;init
  Rgas          = 8.3314e7 ;gas constant, erg/K/mol
  mass_hydrogen = 1.67e-24 ;g
  boltzmann     = 1.38e-16 ;erg/K
  cm_in_pc      = 3.086e18

  ; config
  ;P_k0  = 756.0          ;K/cm^3
  ;n_0   = 2.93             ;cm^(-3)
  
  P_k0  = 2400.0          ;K/cm^3
  n_0   = 1.0             ;cm^(-3)
  
  fHeat = 0.012   ;constant heating function, erg/s/g
  
  kappa = 0.0     ;conduction, erg/cm/K/s
  ;kappa = 7.48e6
  gamma = 5.0/3.0 ;specific heat ratio
  mu    = 1.22     ;mean molecular weight  
  
  T_0   = P_k0 / n_0                  ;K
  P_0   = P_k0 * boltzmann            ;ba
  rho_0 = n_0 * (mu * mass_hydrogen)  ;g/cm^3
  
  csnd  = (gamma * P_0 / rho_0)^(0.5) ;cm/s

  ; compute partial derivatives of heatloss function (3pt Lagrangian method)
  ;pdXSpread = [0.98,0.99,1.0,1.01,1.02]
  pdXSpread = [0.80,0.90,1.0,1.10,1.20]
  fHL_X = double(T_0) * pdXSpread
  fHL_Y = rho_0 * KI02(fHL_X)/mass_hydrogen/mass_hydrogen - fHeat
  
  pdHL_T   = ( deriv(fHL_X,fHL_Y) )[2]
  
  fHL_X = double(rho_0) * pdXSpread
  fHL_Y = fHL_X * KI02(T_0)/mass_hydrogen/mass_hydrogen - fHeat
  
  pdHL_rho = ( deriv(fHL_X,fHL_Y) )[2] 
  ;pdHL_rho = LambdaArrSS02(T_0) ;analytically, equivalent

  ; compute wavenumbers [1/cm]
  k_T   = (mu * (gamma - 1.0) * pdHL_T) / (Rgas * csnd)
  k_rho = (mu * (gamma - 1.0) * rho_0 * pdHL_rho) / (Rgas * csnd * T_0)
  k_K   = (Rgas * csnd * rho_0) / (mu * (gamma - 1.0) * kappa)
stop
  ;k_rho = 2.0 * k_T   ;a = 1/2
  ;k_K   = 10.0 * k_rho ;beta=0.2
  
  ;compute kappa=0 k->inf limit of growth rate
  n_1 = -1.0 * (gamma - 1) * mu * (T_0 * pdHL_T - rho_0 * pdHL_rho) / (gamma * Rgas * T_0)
  n_1_norm = n_1 / (k_rho * csnd)

  ; wavenumber points
  ;lambdaMin = alog10(0.1*cm_in_pc) ;cm
  ;lambdaMax = alog10(1e4*cm_in_pc) ;cm
  ;lambdaRes = 100
  
  ;lambdaPts = findgen(lambdaRes)/lambdaRes * (lambdaMax - lambdaMin) + lambdaMin
  ;lambdaPts = 10.0^(lambdaPts)
  
  kMin = 0.01*k_rho ;norm
  kMax = 3*k_rho   ;norm
  kRes = 100
  
  kPts = findgen(kRes)/kRes * (kMax - kMin) + kMin
  
  ; solve cubic characteristic equation
  roots = fltarr(3,n_elements(kPts))
  
  c0 = csnd^3.0 * kPts^2.0 / gamma * (k_T - k_rho + kPts^2.0/k_K)
  c1 = csnd^2.0 * kPts^2.0
  c2 = csnd * (k_T + kPts^2.0/k_K)
  c3 = 1.0

  for i=0,n_elements(kPts)-1 do begin
    roots[*,i] = cuberoot([c0[i],c1[i],c2[i],c3])
    ;roots[*,i] = cubic(0.0,[c0[i],c1[i],c2[i],c3])
  endfor
  
  ;plot in cgs  
    ;plotX = 2.0 * !pi / kPts / cm_in_pc   ;pc
  
  plotY = fltarr(n_elements(kPts)) ;handle cuberoot() return
  for i=0,n_elements(kPts)-1 do begin
    if ( roots[0,i] > 0 ) then begin
      plotY[i] = roots[0,i]              ;s^(-1)
    endif
    if ( roots[0,i] < 0 and roots[1,i] > 0 ) then begin
      plotY[i] = roots[1,i]              ;s^(-1)
    endif
  endfor
  
  ; plot normalized
  plotX = kPts / k_rho
  plotY = roots[0,*] / (k_rho * csnd)
  
  PS_Start, FILENAME=workingPath+"growthRateKI02.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=5, /inches
            
    ;plot in cgs
    ;plot,[0],[0],xtitle="Wavelength "+textoidl('\lambda')+" [pc]", $
    ;     ytitle="Growth Rate n [s"+textoidl('^{-1}')+"]", $
    ;     thick=1.1,charsize=1.5,/xstyle,/xlog,/ystyle,/nodata,xrange=[0.1,max(plotX)], $
    ;     yrange=[0,max(plotY)*1.1]
    ;oplot,plotX,plotY,psym=3;,color=fsc_color('green') 
         
    ;plot normalized
    fsc_plot,[0],[0],xtitle="Normalized Wavenumber k / k"+textoidl('_\rho')+"", $
         ytitle="Normalized Growth Rate n / (k"+textoidl('_\rho')+" c"+textoidl('_s')+")", $
         thick=1.1,charsize=1.5,/xstyle,/ystyle,/nodata,xrange=minmax(plotX), $
         yrange=[0,max(plotY)*1.2]     
    fsc_plot,plotX,plotY,psym=4,/overplot;,color=fsc_color('green') 
    fsc_plot,[minmax(plotX)],[n_1_norm,n_1_norm]*1.01,line=2,/overplot

  PS_End
  
end

; compTI():
;  - ArepoStatic / ArepoMoving / Gadget TI growth comparison

pro compTI

  ;config
  units = getUnits()
  
  xrange = [0.0,0.1]
  
  workingPath     = '/n/home07/dnelson/compTI/compTI.w01.'
  arepoStaticPath = '/n/home07/dnelson/compTI/sim.arepo.moving.w01.eigen05/output/snap_'
  arepoMovingPath = '/n/home07/dnelson/compTI/sim.arepo.moving.w01/output/snap_' 
  gadgetPath = '/n/home07/dnelson/compTI/sim.gadget.w01/output/snap_'  
  
  ;load one time and plot density comparison
  snapNums = [25,185,200]
  
  PS_Start, FILENAME=workingPath+"threepanel.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=5.0, /inches      
  
  xm = !x.margin
  ym = !y.margin
  !y.margin = [2.0,4.0]
  !x.margin = [7.0,2.0]
  !p.multi = [0,3,3]
  
    sn = 20
    h1 = loadSnapshotHDF5(arepoStaticPath,sn,pos1,vel1,id1,mass1,u1,rho1,hsml1)
    h2 = loadSnapshotHDF5(arepoMovingPath,sn,pos2,vel2,id2,mass2,u2,rho2,hsml2)
    h3 = loadSnapshotHDF5(gadgetPath,sn,pos3,vel3,id3,mass3,u3,rho3,hsml3)
    
    fsc_plot,pos1[0,*],vel1[0,*],line=0,xrange=xrange,/xs,charsize=1.5, $
         xtitle="",ytitle="vel, rho, u",xtickname=replicate(' ',10)
    fsc_plot,pos1[0,*],rho1,     line=0,xrange=xrange,/xs,charsize=1.5, $
         xtitle="",ytitle="",xtickname=replicate(' ',10)
    fsc_plot,pos1[0,*],u1,       line=0,xrange=xrange,/xs,charsize=1.5, $
         xtitle="",ytitle="",xtickname=replicate(' ',10)
         
    fsc_plot,pos2[0,*],vel2[0,*],line=0,xrange=xrange,/xs,charsize=1.5, $
         xtitle="",ytitle="vel, rho, u"
    fsc_plot,pos2[0,*],rho2,     line=0,xrange=xrange,/xs,charsize=1.5, $
         xtitle="",ytitle=""
    fsc_plot,pos2[0,*],u2,       line=0,xrange=xrange,/xs,charsize=1.5, $
         xtitle="",ytitle=""
         
    fsc_plot,pos3[0,*],vel3[0,*],line=0,xrange=xrange,/xs,charsize=1.5, $
         xtitle="x",ytitle="vel, rho, u"
    fsc_plot,pos3[0,*],rho3,     line=0,xrange=xrange,/xs,charsize=1.5, $
         xtitle="x",ytitle=""
    fsc_plot,pos3[0,*],u3,       line=0,xrange=xrange,/xs,charsize=1.5, $
         xtitle="x",ytitle=""
        
    fsc_text,0.1,0.95,"t="+string(h1.time*1000,format='(f4.1)'),/normal,alignment=0.5,charsize=1.5  
         
  !p.multi = 0
  !x.margin = xm
  !y.margin = ym
    
  PS_End  
  
  PS_Start, FILENAME=workingPath+"dens.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=8.0, /inches    
    
  xm = !x.margin
  ym = !y.margin
  !y.margin = [1.0,1.0]
  !x.margin = [4.0,2.0]
  !p.multi = [0,1,n_elements(snapNums)]  
    
  for i=0,n_elements(snapNums)-1 do begin  
  
    h1 = loadSnapshotHDF5(arepoStaticPath,snapNums[i],pos1,vel1,id1,mass1,u1,rho1,hsml1)
    h2 = loadSnapshotHDF5(arepoMovingPath,snapNums[i],pos2,vel2,id2,mass2,u2,rho2,hsml2)  
    h3 = loadSnapshotHDF5(gadgetPath,snapNums[i],pos3,vel3,id3,mass3,u3,rho3,hsml3)
    
    ;subset of cell positions to plot
    select = indgen(204)*5
    pos1_x = pos1[0,[select]]
    pos1_y = fltarr(n_elements(select))+0.0011
    pos2_x = pos2[0,[select]]
    pos2_y = fltarr(n_elements(select))+0.0012  
    pos3_x = pos3[0,[select]]
    pos3_y = fltarr(n_elements(select))+0.0013      
        
    fsc_plot,pos1[0,*],rho1,line=0,xrange=xrange,/xs,/ylog,charsize=1.5, $
         xtitle="",ytitle="",xtickname=replicate(' ',10)
    fsc_plot,pos2[0,*],rho2*1.0,line=0,color=fsc_color('red'),/overplot
    fsc_plot,pos3[0,*],rho3*1.0,line=0,color=fsc_color('orange'),/overplot
    
    fsc_plot,pos1_x,pos1_y,psym=9,symsize=0.4,/overplot ;3=dot,9=open circle
    fsc_plot,pos2_x,pos2_y,psym=9,symsize=0.4,/overplot,color=fsc_color('red')
    fsc_plot,pos3_x,pos3_y,psym=9,symsize=0.4,/overplot,color=fsc_color('orange') 
    
    fsc_text,0.93,0.95-0.335*i,"t="+string(h1.time*1000,format='(f4.1)'),/normal,alignment=0.5,charsize=1.5  
            
  endfor
  
  !p.multi = 0
  !x.margin = xm
  !y.margin = ym
            
  PS_End  
  
  ;load time sequence and comparison growth rates

end

; compGrowthRates():

pro compGrowthRates

  ;config
  units = getUnits()

  linStart = 10 ;snapshot (~1 Myr, let initial bumps die down)
  linStop  = 80 ;snapshot (~4 Myr)
  
  colors = ['blue','brown','orange','cyan','deep pink','red','orchid','green']
  
  ; cut from PlotGrowthRateSS02 for: P_k0 = 2400.0,  n_0 = 1.0, gamma = 1.667, mu = 1.22
    csnd  = 520565.0
    k_rho = 1.3817e-19

  workingPath = '/n/home07/dnelson/compTI.w0/'
  simPostfix  = '/output/snap_'
  
  simNames = ['sim.e2.20pc.w01','sim.e2.20pc.w005','sim.e2.20pc.w002','sim.e2.20pc.w001',$
              'sim.e2.35pc.w012','sim.e2.35pc.w01','sim.e2.35pc.w008',$
              'sim.e2.35pc.w005','sim.e2.35pc.w002','sim.e2.35pc.w001',$
              'sim.e2.50pc.w01','sim.e2.50pc.w005','sim.e2.50pc.w002','sim.e2.50pc.w001']

  eigenVals = [2.0, 2.0, 2.0, 2.0,$
               2.0, 2.0, 2.0,$
               2.0, 2.0, 2.0,$
               2.0, 2.0, 2.0, 2.0]
               
  w0Vals   = [0.1, 0.05, 0.02, 0.01,$
              0.12,0.10, 0.08, $
              0.05,0.02, 0.01, $
              0.1, 0.05, 0.02, 0.01]
               
  boxSizes = [20.0, 20.0, 20.0, 20.0,$
              35.0, 35.0, 35.0, $
              35.0, 35.0, 35.0, $
              50.0, 50.0, 50.0, 50.0]
          
  ; grouping config    
  mf = [0,4,6,4]

  growthRate2 = fltarr(n_elements(simNames))

  for i=0,n_elements(simNames)-1 do begin
    ; load snapshots
    fBase = workingPath + simNames[i] + simPostfix
    nSnaps = n_elements(file_search(fBase+"*"))
  
    maxRho = fltarr(nSnaps)
    times  = fltarr(nSnaps)
  
    ; save max density with time
    for j=linStart,linStop do begin
      h = loadSnapshotHDF5(fBase,j,pos,vel,id,mass,u,rho,hsml)

      maxRho[j] = max(rho)
      times[j]  = h.time
    endfor
    
    ; compute growth rate
    deltaTime = (times[linStop] - times[linStart]) * units.UnitTime_in_s ;s
    deltaDens2 = alog(maxRho[linStop]*units.UnitDensity_in_cgs) - alog(maxRho[linStart]*units.UnitDensity_in_cgs)
    
    growthRate2[i] = deltaDens2 / deltaTime
    growthRate2[i] /= (csnd * k_rho) ;normalize
    
    print,simNames[i]+'     = ',growthRate2[i]
  
  endfor
  
  ; plot
  ;PS_Start, FILENAME=workingPath+"compGrowthRates.eigen.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
  ;          /encapsulated,decomposed=0, xs=8.0, ys=8.0, /inches  
            
  ;  fsc_plot,[0],[0],line=0,xrange=[min(eigenVals)-1,max(eigenVals)+1],/xs,$
  ;           yrange=[min(growthRate2)*0.95,max(growthRate2)*1.05],$
  ;           charsize=1.5,xtitle="eigenmode k / k_rho",ytitle="growth rate (norm)"
    
  ;  for i=0,n_elements(mf)-2 do begin
  ;    ind1=total(mf[0:i])
  ;    ind2=total(mf[0:(i+1)])-1
  ;    fsc_plot,eigenVals[ind1:ind2],growthRate2[ind1:ind2],$
  ;             line=0,color=fsc_color(colors[i]),/overplot
  ;    fsc_plot,eigenVals[ind1:ind2],growthRate2[ind1:ind2],psym=4,/overplot
  ;  endfor
    
  ;PS_End
  
  PS_Start, FILENAME=workingPath+"compGrowthRates.w0.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=8.0, /inches  
            
    fsc_plot,[0],[0],line=0,xrange=[min(w0Vals)*0.9,max(w0Vals)*1.1],/xs,$
             yrange=[min(growthRate2)*0.95,max(growthRate2)*1.05],$
             charsize=1.5,xtitle="perturbation strength (frac)",ytitle="growth rate (norm)"
    
    for i=0,n_elements(mf)-2 do begin
      ind1=total(mf[0:i])
      ind2=total(mf[0:(i+1)])-1
      print,ind1,ind2
      fsc_plot,w0Vals[ind1:ind2],growthRate2[ind1:ind2],$
               line=0,color=fsc_color(colors[i]),/overplot
      fsc_plot,w0Vals[ind1:ind2],growthRate2[ind1:ind2],psym=4,/overplot
    endfor
    
  PS_End 
            
  stop
end

; TI_1D_Movie()

pro TI_1D_Movie

  units = getUnits()

  ; config
  simName = 'sim.eigen2.1n2.nocool'
  xrange  = [0,0.035]
  
  workingPath = '/n/home07/dnelson/dev.thermcond/'
  fBase       = workingPath + simName + '/output/snap_'
  fBaseOut    = workingPath + simName + '/output2/frame_'

  nSnaps = n_elements(file_search(fBase+"*"))

  ; plotting snap by snap, establish const limits
  yrange_vel = [1e8,0]
  yrange_rho = [1e8,0]
  yrange_u   = [1e8,0]
  
  for i=0,nSnaps-1,1 do begin ;for i=2,...
    h = loadSnapshotHDF5(fBase,i,pos,vel,id,mass,u,rho,hsml)
    yrange_vel = [min([yrange_vel[0],min(vel[0,*])]), $
                  max([yrange_vel[1],max(vel[0,*])])]
    yrange_rho = [min([yrange_rho[0],min(rho)]), $
                  max([yrange_rho[1],max(rho)])]
    yrange_u   = [min([yrange_u[0],min(u)]), $
                  max([yrange_u[1],max(u)])]                  
  endfor

  yrange_vel *= [0.5,1.5]
  yrange_rho *= [0.8,1.2]
  yrange_u   *= [0.5,1.5]

  ;plot in physical units?
  yrange_rho *= units.UnitDensity_in_cgs / units.mass_proton 
  yrange_vel *= units.UnitVelocity_in_cm_per_s / 1e5
  ;yrange_u *= units.UnitEnergy_in_cgs ;erg (now code)

  ;load initial
  h0 = loadSnapshotHDF5(fBase,0,pos0,vel0,id0,mass0,u0,rho0,hsml0)
  rho0_cgs = rho0 * units.UnitDensity_in_cgs / units.mass_proton ;cm^(-3)
  vel0_cgs = vel0 * units.UnitVelocity_in_cm_per_s / 1e5 ;cm/s
  u0_cgs   = u0 ; code   * units.UnitEnergy_in_cgs ;erg

  for i=0,nSnaps-1,1 do begin
    h = loadSnapshotHDF5(fBase,i,pos,vel,id,mass,u,rho,hsml)

    ;plot in physical units?
    rho_cgs = rho * units.UnitDensity_in_cgs / units.mass_proton ;cm^(-3)
    vel_cgs = vel * units.UnitVelocity_in_cm_per_s / 1e5 ;km/s
    u_cgs   = u ; code   * units.UnitEnergy_in_cgs ;erg
    
    ;SPH output should be sorted for plotting
    ;sortkeys     = sort(pos[0,*])
    ;rho_cgs      = rho_cgs[sortkeys]
    ;vel_cgs[0,*] = vel_cgs[0,sortkeys]
    ;u_cgs        = u_cgs[sortkeys]
    ;pos[0,*]     = pos[0,sortkeys]    
    
    PS_Start, FILENAME=fBaseOut+string(i,format='(I04)')+".eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
              /encapsulated,decomposed=0, xs=6.0, ys=8.0, /inches    
    
      !p.multi = [0,1,3]
      xm = !x.margin
      ym = !y.margin
      ;!x.margin = [0,0]
      
      curTitle = simName + " t = " + string(h.time*1000,format='(f5.2)')
          
      ; plot 1D density/energy profiles
      !y.margin = [0,3]
      fsc_plot,pos[0,*],vel_cgs[0,*],line=0,xrange=xrange,yrange=yrange_vel,/ys,/xs,charsize=1.5, $
           xtitle="",ytitle=textoidl("velocity [km/s]"),xtickname=replicate(' ',10),title=curTitle
      fsc_plot,pos0[0,*],vel0_cgs[0,*],line=1,/overplot
      
      !y.margin = [0,0]
      fsc_plot,pos[0,*],rho_cgs,     line=0,xrange=xrange,yrange=yrange_rho,/ys,/xs,charsize=1.5, $
           xtitle="",ytitle=textoidl("number density [1/cm^3]"),xtickname=replicate(' ',10),/ylog
      fsc_plot,pos0[0,*],rho0_cgs,line=1,/overplot
      
      !y.margin = [3,0]
      fsc_plot,pos[0,*],u_cgs,       line=0,xrange=xrange,yrange=yrange_u,/ys,/xs,charsize=1.5, $
           xtitle="x",ytitle=textoidl("internal energy [code]"),/ylog
      fsc_plot,pos0[0,*],u0_cgs,line=1,/overplot
      
      !x.margin = xm
      !y.margin = ym
      !p.multi = 0     
    
    PS_End, /PNG, Resize=40, /Delete_PS ;PNG size=[xs,ys]*300*(resize/100)
    
    vel2    = vel[0,*]*vel[0,*] + vel[1,*]+vel[1,*] + vel[2,*]*vel[2,*]
    KEtoU   = total(0.5*vel2) / total(u)
 
    print,'[i] time max_rho total_mass ke/u',i,h.time,max(rho),total(mass),KEtoU
  endfor

end

; TI_1D():
;  - run analysis on snapshot sequence

pro TI_1D

  units = getUnits()

  ; config
  simName = 'sim.e2.35pc.w005'
  
  workingPath = '/n/home07/dnelson/compTI.w0/'
  fBase       = workingPath + simName + '/output/snap_'
  fBaseOut    = workingPath + simName + '.'
    
  ; find number of snapshots
  nSnaps = n_elements(file_search(fBase+"*"))

  ; cut from PlotGrowthRateSS02 for: P_k0 = 2400.0,  n_0 = 1.0, gamma = 1.667, mu = 1.22
    csnd  = 520565.0
    k_rho = 1.3817e-19

  ; quantities
  avgVel  = fltarr(nSnaps)
  avgRho  = fltarr(nSnaps)
  totMass = fltarr(nSnaps)
  totU    = fltarr(nSnaps)
  totKE   = fltarr(nSnaps)
  
  KEtoU   = fltarr(nSnaps)
  
  maxRho  = fltarr(nSnaps)
  times   = fltarr(nSnaps)
  
  for i=0,nSnaps-1,1 do begin
    h = loadSnapshotHDF5(fBase,i,pos,vel,id,mass,u,rho,hsml)
    
    ; compute quantities for this timestep
    avgVel[i]  = mean(vel[0,*])
    avgRho[i]  = mean(rho)
    totMass[i] = total(mass)
    totU[i]    = total(u)
    
    vel2       = vel[0,*]*vel[0,*] + vel[1,*]+vel[1,*] + vel[2,*]*vel[2,*]
    totKE[i]   = total(0.5*vel2)
    KEtoU[i]   = totKE[i] / totU[i]
 
    maxRho[i]  = max(rho)
    times[i]   = h.time

    print,'[i] max_rho total_mass ke/u',i,maxRho[i],totMass[i],KEtoU[i]
  endfor

  ; convert code time to Myr
  plotX = times * units.UnitTime_in_s / units.s_in_Myr

  ; timeRange = [0.1,max(plotX)*1.5] ;log
  timeRange = [0,10]

  ; plot quantities over time  
  PS_Start, FILENAME=fBaseOut+"runValsMeas.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=5, /inches

    xm = !x.margin
    ym = !y.margin
    ;!x.margin = [4,2]
    ;!y.margin = [0,0]
    !p.multi = [0,2,2]
    
    fsc_plot,plotX,totKE,psym=-3,xrange=timeRange,yrange=[1e-3,max(totKE)*1.5], $
         /ys,/xs,/ylog,charsize=1.2,xtitle="", $
         ytitle=textoidl("avg kinetic energy")      
    fsc_plot,plotX,maxRho,xrange=timeRange,yrange=[min(maxRho)*0.9,max(maxRho)*1.1], $
         /ys,/xs,/ylog,charsize=1.2,xtitle="", $
         ytitle=textoidl("max density");,psym=-3

    fsc_plot,plotX,KEtoU,psym=-3,yrange=[1e-8,max(KEtoU)*1.5],xrange=timeRange, $
         /ys,/xs,/ylog,charsize=1.2, $
         xtitle="time [Myr]",ytitle=textoidl("kinetic to thermal energy")
    fsc_plot,plotX,totU,psym=-3,xrange=timeRange,yrange=[min(totU)*0.9,max(totU)*1.1], $
         /ys,/xs,/ylog,charsize=1.2, $
         xtitle="time [Myr]",ytitle=textoidl("total thermal energy") 
            
    !p.multi = 0  
    !x.margin = xm
    !y.margin = ym
  
  PS_End
  
  ; calculate a simple slope in the linear regime
  linStart = 40 ;snapshot (~2 Myr, let initial bumps die down)
  linStop = 100 ;snapshot (~10 Myr)
  
  deltaTime = (times[linStop] - times[linStart]) * units.UnitTime_in_s ;s
  ;deltaDens = (maxRho[linStop] - maxRho[linStart]) * units.UnitDensity_in_cgs
  deltaDens2 = alog(maxRho[linStop]*units.UnitDensity_in_cgs) - alog(maxRho[linStart]*units.UnitDensity_in_cgs)
  
  ;growthRate = deltaDens / deltaTime ;cgs
  ;growthRate /= (csnd * k_rho) ;normalize
  
  ;print,'growthRate = ',growthRate
  
  growthRate2 = deltaDens2 / deltaTime
  growthRate2 /= (csnd * k_rho) ;normalize
  
  print,'growthRate (log) = ',growthRate2
  
  ; calculate logarithmic change of density with time
  nLogDens = 5
  densRate  = fltarr(n_elements(times)-nLogDens)
  densRate2 = fltarr(n_elements(times)-nLogDens)
  plotX     = fltarr(n_elements(times)-nLogDens)
  
  for i=1,n_elements(times)-nLogDens-1,1 do begin
  
    yDens    = alog( maxRho[i:i+nLogDens-1] * units.UnitDensity_in_cgs )
    xDens    = times[i:i+nLogDens-1] * units.UnitTime_in_s
  
    densRate[i]  = (yDens[nLogDens-1] - yDens[0]) / (xDens[nLogDens-1] - xDens[0])
    densRate2[i] = (deriv(xDens[0:nLogDens-1], yDens[0:nLogDens-1]))[floor(nLogDens/2)]
    plotX[i]     = mean(xDens) / units.s_in_Myr
  
    ; normalize growth rate
    densRate[i]  /= (csnd * k_rho)
    densRate2[i] /= (csnd * k_rho)

  endfor
  
  w = where(plotX ge timeRange[0] and plotX le timeRange[1])
  
  PS_Start, FILENAME=fBaseOut+"growthRateMeas.eps", /nomatch, /quiet, font=1, bits_per_pixel=8,color=1, $
            /encapsulated,decomposed=0, xs=8.0, ys=5, /inches
  
    fsc_plot,[0],[0],psym=-3,xrange=timeRange,yrange=[min(densRate2[w]*1.3),max(densRate2[w]*1.1)], $
       /ys,/xs,charsize=1.2, $
       xtitle="time [Myr]",ytitle="normalized Growth Rate n / (k"+textoidl('_\rho')+" c"+textoidl('_s')+")";,/xlog;,/ylog,
    ;fsc_plot,plotX,densRate,line=0,/overplot
    fsc_plot,plotX,densRate2,line=0,/overplot

  PS_End

  print,'done.'
end

pro TI_2D

  meshBase   = "c:\zStuff\IDL.work\Default\ctTests.SS02\TI_2D.HIM.512\voronoi_mesh_"
  densBase   = "c:\zStuff\IDL.work\Default\ctTests.SS02\TI_2D.HIM.512\density_field_"
  snapBase   = "c:\zStuff\IDL.work\Default\ctTests.SS02\TI_2D.HIM.512\snap_"
  
  boxSize    = [1.0,1.0]   ;x,y (sim units)
  xyScaleFac = 1.0  
  
  ; find number of snapshots
  nSnaps = n_elements(file_search(densBase+"*"))
  
  ;colorMinMax = getDensityMinMax(densBase, nSnaps) / 2.0
  colorMinMax = [-5.0,15.0]

  for i=0,nSnaps,1 do begin 
    ;h = loadSnapshot(snapBase,i,pos,vel,id,mass,u,rho,hsml)
    sz = loadDensityField2D(densBase, i, dens)
    print,i,minmax(dens*100)
    ;plotDensityField2D, densBase, i, /writeJPG, xyScaleFac=xyScaleFac, colorMinMax=colorMinMax
    plotDensityField2D, densBase, i, /tvOut, xyScaleFac=xyScaleFac, colorMinMax=colorMinMax, dens=dens*100
    ;stop
  endfor

end
