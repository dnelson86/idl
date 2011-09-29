;TI_plot.pro
;plot cooling functions
;dnelson jun.2001

@helper

; ------------------------------------------------------------------------------
; analytic fits
; ------------------------------------------------------------------------------

  ; Dalgarno & McCray (1972):
  ; electronic collisions with singly charged positive ions (C+, Si+, Fe+)
  
; e + C+(2P1/2) -> e + C+(2P3/2)
function L_C1, T
  L_e_C = 7.9*10.0^(-20.0) * T^(-0.5) * exp(-92.0/T)

  return, L_e_C
end

; e + Si+(2P1/2) -> e + Si+(2P3/2)
function L_Si1, T
  L_e_Si = 1.9*10.0^(-18.0) * T^(-0.5) * exp(-413.0/T)

  return, L_e_Si
end

; e + Fe+(a6D9/2) -> e + Fe+(a6D7/2)
; and
; e + Fe+(a6D9/2) -> e + Fe+(a6D5/2)
function L_Fe1, T
  L_e_Fe = 1.1*10.0^(-18.0) * T^(-0.5) * ( exp(-554.0/T) + 1.3*exp(-961.0/T) )
  
  return, L_e_Fe
end

; e + 0(3P2) -> e + O(3P1,0)
function L_O1, T
  L_e_O = 1.74*10.0^(-24.0) * T^(0.5) * ( (1-7.6*T^(-0.5)) * exp(-228.0/T) + $
                                      0.38*(1-7.7*T^(-0.5)) * exp(-326.0/T) )
  
  return, L_e_O
end

; e + O(3P) -> e + O(1D) metastable
function L_O2, T
  L_e_O = 9.4*10.0^(-23.0) * T^(0.5) * exp(-22700.0/T)
  
  return, L_e_O
end

; Koyama & Inutsaka (2002)
; analytical fit to KI02 [erg cm^3/s]
function KI02, T

  L = 10D7 * exp(-114800.0 / (T + 1000)) 
  L += 14.0 * sqrt(T) * exp(-92.0 / T)
  L *= 2.0D-26

  return, L
end

; Sanchez-Salcedo & Vazquez-Semadeni (2002) analytic cooling function
; [Lambda] = [erg / s / g^2 / cm^3 ]

function LambdaArrSS02, T

  ;cooling function
  L = fltarr(n_elements(T))
  
  ;for each temperature regime fill L with appropriate fit
  w = where(T lt 10, count)
  if (count ne 0.0) then $
    L[w] = 0.0
  
  w = where((10 le T) and (T lt 141), count)
  if (count ne 0.0) then $
    L[w] = (3.42 * 1e16) * T[w]^(2.12)
    
  w = where((141 le T) and (T lt 313), count)
  if (count ne 0.0) then $
    L[w] = (9.10 * 1e18) * T[w]^(1.0)
    
  w = where((313 le T) and (T lt 6102), count)
  if (count ne 0.0) then $
    L[w] = (1.11 * 1e20) * T[w]^(0.56)
    
  w = where((6102 le T) and (T le 10L^5), count)
  if (count ne 0.0) then $
    L[w] = (2.00 * 1e8) * T[w]^(3.67)
    
  w = where(10L^5 lt T, count)
  if (count ne 0.0) then $
    L[w] = 0.0

  ; return single value if requested
  if (n_elements(T) eq 1) then $
    return, L[0]
    
  return, L

end

; inverse function (return T given Lambda)

function LambdaInvSS02, L
  
  ;piecewise break points
  L10    = LambdaArrSS02(10)
  L141   = LambdaArrSS02(141)
  L313   = LambdaArrSS02(313)
  L6102  = LambdaArrSS02(6102)
  L10to5 = LambdaArrSS02(10L^5)
  
  ;cooling function
  T = fltarr(n_elements(L))
  
  ;for each temperature regime fill L with appropriate fit
  w = where(L lt L10, count)
  if (count ne 0.0) then $
    T[w] = 0.0
  
  w = where((L10 le L) and (L lt L141), count)
  if (count ne 0.0) then $
    T[w] = (L[w] / (3.42 * 1e16))^(1/2.12)

  w = where((L141 le L) and (L lt L313), count)
  if (count ne 0.0) then $
    T[w] = (L[w] / (9.10 * 1e18))^(1/1.00)
  

  w = where((L313 le L) and (L lt L6102), count)
  if (count ne 0.0) then $
    T[w] = (L[w] / (1.11 * 1e20))^(1/0.56)

  w = where((L6102 le L) and (L lt L10to5), count)
  if (count ne 0.0) then $
    T[w] = (L[w] / (2.00 * 1e8))^(1/3.67)

  w = where(L10to5 le L, count)
  if (count ne 0.0) then $
    T[w] = 0.0

  ; return single value if requested
  if (n_elements(L) eq 1) then $
    return, T[0]
    
  return, T

end

;Rosen & Bregman (1995)
;single return version
function RB95, T

  if (T lt 300) then $
    L = 0.0
    
  if ((300 le T) and (T lt 2000)) then $
    L = (2.2380 * 10L^(-32.0)) * T^(2.0)
    
  if ((2000 le T) and (T lt 8000)) then $
    L = (1.0012 * 10L^(-30.0)) * T^(1.5)
    
  if ((8000 le T) and (T lt 10L^(5.0))) then $
    L = (4.6240 * 10L^(-36.0)) * T^(2.867)
    
  if ((10L^(5.0) le T) and (T lt 4.0*10L^(7.0))) then $
    L = (1.7800 * 10L^(-18.0)) * T^(-0.65)
    
  if (4.0*10L^(7.0) le T) then $
    L = (3.2217 * 10L^(-27.0)) * T^(0.5)

  return, L
end

;array version
function RB95arr, T
  
  ;cooling function
  L = fltarr(n_elements(T))
  
  ;for each temperature regime fill L with appropriate fit
  w = where(T lt 300, count)
  if (count ne 0.0) then $
    L[w] = 0.0
  
  w = where((300 le T) and (T lt 2000), count)
  if (count ne 0.0) then $
    L[w] = (2.2380 * 10L^(-32.0)) * T[w]^(2.0)
    
  w = where((2000 le T) and (T lt 8000), count)
  if (count ne 0.0) then $
    L[w] = (1.0012 * 10L^(-30.0)) * T[w]^(1.5)
    
  w = where((8000 le T) and (T lt 10L^(5.0)), count)
  if (count ne 0.0) then $
    L[w] = (4.6240 * 10L^(-36.0)) * T[w]^(2.867)
    
  w = where((10L^(5.0) le T) and (T lt 4.0*10L^(7.0)), count)
  if (count ne 0.0) then $
    L[w] = (1.7800 * 10L^(-18.0)) * T[w]^(-0.65)
    
  w = where(4.0*10L^(7.0) le T, count)
  if (count ne 0.0) then $
    L[w] = (3.2217 * 10L^(-27.0)) * T[w]^(0.5)

  return, L
end

; ------------------------------------------------------------------------------
; create interpolated tables from GK10
; use GS07 data for 7<logT<8
; ------------------------------------------------------------------------------
pro makeTables

  fileName1 = 'C:\zStuff\IDL.work\Default\cooling\gk10\dl.Z=1.0_UV=1.res'
  fileName2 = 'C:\zStuff\IDL.work\Default\cooling\gk10\dl.Z=0.1_UV=1.res'
  fileName3 = 'C:\zStuff\IDL.work\Default\cooling\gk10\dl.Z=0.01_UV=1.res'
  fileName4 = 'C:\zStuff\IDL.work\Default\cooling\gk10\dl.Z=0.001_UV=1.res'

  ;note for Z=0.1 have T=7.2 7.3 points (delete 2 lines from GS07 in output)
  fileNameGK10 = fileName1

  ;load GK10 data
  tMin = 1e1
  tMax = 1e8

  headerLines = 0

  ptStruct = {n_b:0.0,  T:0.0, col_b:0.0, Z:0.0, f_HI:0.0, f_HII:0.0, f_H2:0.0, U_MW:0.0, $
              D_MW:0.0, rho_over:0.0, fc:0.0, fh:0.0}
              
  z1 = loadCSV(headerLines,fileNameGK10,ptStruct)

  ;load GS07 data
  zVals = [0.001,0.01,0.1,1.0,2.0]
  headerLines = 23
  fileNameGS07 = 'C:\zStuff\IDL.work\Default\cooling\gs07\tab13.txt'
  ptStruct2 = {T:0.0, Z1:0.0, Z2:0.0, Z3:0.0, Z4:0.0, Z5:0.0}
  
  fc = loadCSV(headerLines,fileNameGS07,ptStruct2)
  
  ;loglinear interpolation
  z1.T   = alog10(z1.T)
  fc.T   = alog10(fc.T)
              
  ;averaging/interpolation bins
  nbLogBins = [-4.0,-2.0,0.0,2.0,4.0]
  tLogBins  = findgen(62)*0.1 + 1.95
  
  ;bin for table
  binnedCool = fltarr(n_elements(tLogBins)-1)
  binnedHeat = fltarr(n_elements(tLogBins)-1)
  binnedT    = fltarr(n_elements(tLogBins)-1)
  binnedNe   = fltarr(n_elements(tLogBins)-1)
  
  ;open output file
  openW, lunW, fileNameGK10+'.txt', /GET_LUN
  
  ;bin GK10 in T (only) for now   
  for j=0,n_elements(tLogBins)-2 do begin
  
    tPts = where( (z1.T ge tLogBins[j]) and (z1.T lt tLogBins[j+1]), count)
                  
    if ( count ne 0 ) then begin
      binnedT[j]    = mean( [tLogBins[j],tLogBins[j+1]] )
      binnedCool[j] = mean( z1[tPts].fc / (0.76*z1[tPts].n_b)^1.0 ) ;normalized by n_H
      ;binnedHeat[j] = mean( z1[tPts].fh / (0.76*z1[tPts].n_b)^1.0 ) ;normalized by n_H
      binnedHeat[j] = mean( z1[tPts].fh  )
      binnedNe[j]   = mean( 0.88 * z1[tPts].f_HII / 0.76 )      ;fraction, normalized by n_H
      
      tVar = variance( alog10(z1[tPts].fc / (0.76*z1[tPts].n_b)) )
      
      ;write line
      lineStr = string(binnedT[j],   format="(f5.3)") + ' ' + $
                string(alog10(binnedCool[j]),format="(f8.4)") + ' ' + $
                string(alog10(binnedHeat[j]),format="(f8.4)") + ' ' + $
                string(binnedNe[j],  format="(e9.2)")
      print,string(j)+' ['+string(count)+'] ('+string(tVar)+') '+lineStr
      printF, lunW, lineStr
    endif
      
  endfor
  
  ;bin GS07 in T for high temp
  ;w = where(tLogBins eq 7.15)
  w = [52]
  for j=w[0], n_elements(tLogBins)-2 do begin
  
    tPts = where( (fc.T ge tLogBins[j]) and (fc.T lt tLogBins[j+1]), count)
                  
    if ( count ne 0 ) then begin
      binnedT[j]    = mean( [tLogBins[j],tLogBins[j+1]] )
      binnedCool[j] = mean( fc[tPts].Z1 )
      binnedNe[j]   = 0.88 / 0.76 ;frac, normed by n_H, for X=0.76 Y=0.24 (Z~=0)
      
      tVar = variance( alog10(fc[tPts].Z1) )
      
      ;write line
      lineStr = string(binnedT[j],   format="(f5.3)") + ' ' + $
                string(alog10(binnedCool[j]),format="(f8.4)") + ' ' + $
                '000.0000' + ' ' + $ ;temp binnedHeat
                string(binnedNe[j],  format="(e9.2)")
      print,string(j)+' ['+string(count)+'] ('+string(tVar)+') '+lineStr
      printF, lunW, lineStr
    endif
  endfor
  
  ;close output file
  close, lunW

  set_plot,'PS'
  device,filename=fileNameGK10+'.eps',bits_per_pixel=8,color=1,/landscape,/encapsulated,/decomposed
         
  ;plot binned results
    plot,[0],[0],xtitle="Temperature [K]", $
         ytitle="Cooling Efficiency [erg cm"+textoidl('^3')+" s"+textoidl('^{-1}')+"]", $
         thick=1.1,charsize=1.5,/xlog,/ylog,/xstyle,/ystyle,/nodata,xrange=[tMin,tMax], $
         yrange=[10.0^(-32.0),10.0^(-20.0)]
    
    oplot,10.0^(z1.T),(z1.fc/(0.76*z1.n_b)),psym=3,color=fsc_color('green')
    oplot,10.0^(z1.T),(z1.fh/(0.76*z1.n_b)),psym=3,color=fsc_color('yellow')    
    oplot,10.0^(binnedT),(binnedCool),thick=3.0,color=fsc_color('orange')
    oplot,10.0^(binnedT),(binnedHeat),thick=3.0,color=fsc_color('red')
    oplot,10.0^(binnedT),(binnedCool-binnedHeat),thick=3.0,color=fsc_color('black')
    
    w = where( (binnedCool-binnedHeat) lt 0, count)
    if (count) then begin
      oplot,10.0^(binnedT),-1.0*(binnedCool[w]-binnedHeat[w]),thick=3.0,color=fsc_color('black'),line=2
    endif
 
  
  device,/close_file
  set_plot,'WIN'
  
  print,'done.'
  stop
  
end


; ------------------------------------------------------------------------------
; visualization
; ------------------------------------------------------------------------------

;for comparison
pro plotMulti

  tMin = 1e1
  tMax = 1e8

  ;GS07
  headerLines = 23
  fileName = 'C:\zStuff\IDL.work\Default\cooling\gs07\tab13.txt'
  ptStruct = {T:0.0, Z1:0.0, Z2:0.0, Z3:0.0, Z4:0.0, Z5:0.0}
  
  fc = loadCSV(headerLines,fileName,ptStruct)
  
  ;GK10
  headerLines = 0
  fileName1 = 'C:\zStuff\IDL.work\Default\cooling\gk10\dl.Z=1.0_UV=1.res'
  fileName2 = 'C:\zStuff\IDL.work\Default\cooling\gk10\dl.Z=0.01_UV=1.res'
  ptStruct = {n_b:0.0,  T:0.0, C3:0.0, C4:0.0, C5:0.0, C6:0.0, C7:0.0, C8:0.0, $
              C9:0.0, C10:0.0, fc:0.0, fh:0.0}
              
  z1 = loadCSV(headerLines,fileName1,ptStruct)
  z2 = loadCSV(headerLines,fileName2,ptStruct)
  
  ;RB95
  tRes = 1000
  tPts    = findgen(tRes)/tRes * tMax + tMin
  tPtsLog = 10.0^( findgen(tRes)/tRes * alog10(tMax) + alog10(tMin) )
  
  ;open ps
  set_plot,'PS'
  device,filename='C:\zStuff\IDL.work\Default\cooling\multi.eps',bits_per_pixel=8,color=1, $
         /landscape,/encapsulated,/decomposed
  
    ;setup plot
    plot,[0],[0],xtitle="Temperature [K]", $
         ytitle="Cooling Efficiency [erg cm"+textoidl('^3')+" s"+textoidl('^{-1}')+"]", $
         thick=1.1,charsize=1.5,/xlog,/ylog,/xstyle,/ystyle,/nodata,xrange=[tMin,tMax], $
         yrange=[10.0^(-27.0),10.0^(-20.0)]
         
    ;GK10
    ;oplot,z1.T,((z1.fh + z1.fc)/z1.n_b),psym=3,color=fsc_color('green')
    ;oplot,z2.T,((z2.fh + z2.fc)/z2.n_b),psym=3,color=fsc_color('orange')
    
    oplot,z1.T,((z1.fc)/(0.76*z1.n_b) - z1.fh),psym=3,color=fsc_color('orange')
    ;oplot,z1.T,((z1.fh) + (z1.fc)/z1.n_b),psym=3,color=fsc_color('green')
    oplot,z1.T,(z1.fc/(0.76*z1.n_b)),psym=3,color=fsc_color('red')
    
    ;GS07
    oplot,fc.T,fc.Z1,line=0,thick=1.5,color=fsc_color('blue')
    oplot,fc.T,fc.Z2,line=1,thick=1.5,color=fsc_color('blue')
    oplot,fc.T,fc.Z3,line=2,thick=1.5,color=fsc_color('blue')
    oplot,fc.T,fc.Z4,line=3,thick=1.5,color=fsc_color('blue')
    oplot,fc.T,fc.Z5,line=4,thick=1.5,color=fsc_color('blue')
    
    ;RB95
    ;oplot,tPtsLog,RB95arr(tPtsLog),thick=2.0,line=0,color=fsc_color('red')
    
    ;KI02
    ;w = where(tPtsLog le 1e5)
    ;oplot,tPtsLog[w],KI02(tPtsLog[w]),thick=2.0,line=0,color=fsc_color('purple')
    
  device,/close_file
  set_plot,'WIN'
    
  stop
end


;Gnedin & Kravstov (2010b)
;cooling function L(T,Z,n_b,U) full temperature range
;metallicity dependence by filename, each file contains n_b,T,heating,cooling
;only U=1 (Milky Way value) for now
;
;Each file contains 12 columns:
;
;1. number density of baryons in cm^{-3} - multiply by 0.76 to get the numberdensity of hydrogen atoms.
;2. temperature in K
;11. cooling rate in erg/s - divide by n_b to get the cooling function per baryon
;12. heating rate in erg/s

; header:
;
;  "baryon number density (cm^{-3}",
;  "temperature (K)",
;  "baryon column density (cm^{-2})",
;  "gas metallicity (solar units)",
;  "HI fraction",
;  "HII fraction",
;  "H_2 fraction",
;  "radiation field in MW units (U_{MW})",
;  "dust-to-gas ratio in MW units (D_{MW}",
;  "density in this cell over the average density in neighboring cells",
;  "cooling rate (erg/s)",
;  "heating rate (erg/s)"

pro plotGK10

  headerLines = 0
  fileName1 = 'C:\zStuff\IDL.work\Default\cooling\gk10\dl.Z=1.0_UV=1.res'
  ;fileName2 = 'C:\zStuff\IDL.work\Default\cooling\gk10\dl.Z=0.01_UV=1.res'
  ;fileName3 = 'C:\zStuff\IDL.work\Default\cooling\gk10\dl.Z=0.001_UV=1.res'
  ptStruct = {n_b:0.0,  T:0.0, C3:0.0, C4:0.0, C5:0.0, C6:0.0, C7:0.0, C8:0.0, $
              C9:0.0, C10:0.0, fc:0.0, fh:0.0}

  z1 = loadCSV(headerLines,fileName1,ptStruct)
  ;z2 = loadCSV(headerLines,fileName2,ptStruct)
  ;z3 = loadCSV(headerLines,fileName3,ptStruct)

  ;plot,gk10.n_b,gk10.T,xtitle="Baryon Number Density [cm^-3]",ytitle="Temperature [K]", $
  ;     thick=1.1,charsize=1.5,psym=3,/xlog,/ylog,symsize=1.5
  set_plot,'PS'
  device,filename='C:\zStuff\IDL.work\Default\cooling\gk10.eps',bits_per_pixel=8,color=1,/landscape,/encapsulated,/decomposed
    plot,z1.T,((z1.fh + z1.fc)/z1.n_b),xtitle="Temperature [K]",ytitle="Cooling Efficiency [erg cm"+textoidl('^3')+" s"+textoidl('^{-1}')+"]", $
         thick=1.1,charsize=1.5,/xlog,/ylog,xrange=[1e1,1e8],yrange=[10L^(-30.0),10L^(-20.0)],/xstyle,/ystyle,/nodata
    oplot,z1.T,(z1.fc/z1.n_b),psym=3,color=fsc_color('green')
    oplot,z1.T,z1.fh,psym=3,color=fsc_color('sky blue')
    ;oplot,z1.T,((z1.fh*1e3 - z1.fc)/z1.n_b),psym=3,color=fsc_color('orange')
    ;oplot,z3.T,(z3.fc/z3.n_b),psym=3,color=fsc_color('green')
  device,/close_file
  set_plot,'WIN'
       
  stop
end

; plotGK10n_b_bin()
;

pro plotGK10n_b_bin

  ;load data
  tMin = 1e1
  tMax = 1e8

  headerLines = 0
  fileName1 = 'C:\zStuff\IDL.work\Default\cooling\gk10\dl.Z=1.0_UV=1.res'
  ptStruct = {n_b:0.0,  T:0.0, col_b:0.0, Z:0.0, f_HI:0.0, f_HII:0.0, f_H2:0.0, U_MW:0.0, $
              D_MW:0.0, rho_over:0.0, fc:0.0, fh:0.0}
              
  z1 = loadCSV(headerLines,fileName1,ptStruct)

  ;loglinear interpolation
  z1.n_b = alog10(z1.n_b)
  z1.T   = alog10(z1.T)

  ;histogram n_b
  histoplot,z1.n_b,axiscolorname='black',thick=1.1,ytitle='',nbins=50,charsize=1.5, $
            xtitle="Log (Baryonic Number Density) [cm"+textoidl("^{-3}")+"]",/fillpolygon, $
            /frequency


  ;plot temp hist by n_b bin
  set_plot,'PS'
  device,filename='C:\zStuff\IDL.work\Default\cooling\nbBinTempHist.eps',bits_per_pixel=8,color=1, $
         /landscape,/encapsulated,/decomposed
         
    !p.multi = [0,2,2,0,0]
    for i=0,n_elements(nbLogBins)-2 do begin
    
      binMin = nbLogBins[i]
      binMax = nbLogBins[i+1]
      
      tPts = where( (z1.n_b ge binMin) and (z1.n_b lt binMax) )
      histoplot,z1[tPts].T,axiscolorname='black',thick=1.1,ytitle='',nbins=50,charsize=1.2, $
              xtitle=string(binMin,format="(f4.1)")+" < log(n"+textoidl("_b")+") < "+ $
              string(binMax,format="(f4.1)"),/fillpolygon,/frequency
    
    endfor
    !p.multi = 0
  
  device,/close_file
  set_plot,'WIN'
  
  ;plot lambda vs T color coded by n_b bin  
  set_plot,'PS'
  device,filename='C:\zStuff\IDL.work\Default\cooling\nbBin.eps',bits_per_pixel=8,color=1, $
         /landscape,/encapsulated,/decomposed
         
    plot,[0],[0],xtitle="Temperature [K]", $
         ytitle="Cooling Efficiency [erg cm"+textoidl('^3')+" s"+textoidl('^{-1}')+"]", $
         thick=1.1,charsize=1.5,/xlog,/ylog,/xstyle,/ystyle,/nodata,xrange=[tMin,tMax], $
         yrange=[10.0^(-27.0),10.0^(-20.0)]
         
    oplot,10.0^(z1.T),(z1.fc/10.0^(z1.n_b)),psym=3,color=fsc_color('black')
    colorNames = ['orange','green','sky blue','red']
    for i=0,n_elements(nbLogBins)-2 do begin
  
      binMin = nbLogBins[i]
      binMax = nbLogBins[i+1]
    
      tPts = where( (z1.n_b ge binMin) and (z1.n_b lt binMax) )
      print,minmax(z1[tPts].T)
      oplot,10.0^(z1[tPts].T),(z1[tPts].fc/10.0^(z1[tPts].n_b)),psym=3,color=fsc_color(colorNames[i])
  
    endfor
       
  device,/close_file
  set_plot,'WIN'
  
  stop

end

;Gnat & Sternberg (2007)
;metallicity dependent cooling function (T > 1e4 regime)
;similar to (updated version of) Sutherland & Dopita (1993)
pro plotGS07

  temps = [18250.0,106000.0,1.047e6,1.096e7,9.12e7]
  zVals = [0.001,0.01,0.1,1.0,2.0]

  ;load cooling function
  headerLines = 23
  fileName = 'C:\zStuff\IDL.work\Default\cooling\gs07\tab13.txt'
  ptStruct = {T:0.0, Z1:0.0, Z2:0.0, Z3:0.0, Z4:0.0, Z5:0.0}
  
  fc = loadCSV(headerLines,fileName,ptStruct)
  
  ;load ionization fractions

  ;open ps
  set_plot,'PS'
  device,filename='C:\zStuff\IDL.work\Default\cooling\gs07.eps',bits_per_pixel=8,color=1,/landscape,/encapsulated,/decomposed
  
    ;plot cooling functions vs T for each metallicity
    plot,fc.T,fc.Z1,xtitle="Temperature [K]", ytitle="Cooling Efficiency [erg cm"+textoidl('^3')+" s"+textoidl('^{-1}')+"]", $
         thick=1.1,charsize=1.5,/xlog,/ylog,/xstyle,/nodata,xrange=[1e4,1e8];,yrange=[10.0^(-27.0),10.0^(-22.0)]
    oplot,fc.T,fc.Z1,line=0
    oplot,fc.T,fc.Z2,line=1
    oplot,fc.T,fc.Z3,line=2
    oplot,fc.T,fc.Z4,line=3
    oplot,fc.T,fc.Z5,line=4
    
    xyouts,10L^6,10L^(-23.5),'Z=0.001',charsize=1.5,alignment=0.5
    xyouts,10L^6,10L^(-23.05),'0.01',charsize=1.5,alignment=0.5
    xyouts,10L^6,10L^(-22.6),'0.1',charsize=1.5,alignment=0.5
    xyouts,10L^6,10L^(-21.75),'1.0',charsize=1.5,alignment=0.5
    xyouts,10L^6,10L^(-21.4),'2.0',charsize=1.5,alignment=0.5
    
    for i=0,n_elements(temps)-1 do begin
      oplot,[temps[i],temps[i]],[1e-24,1e-21],line=2
    endfor
  
  device,/close_file
  set_plot,'WIN'
  
  set_plot,'PS'
  device,filename='C:\zStuff\IDL.work\Default\cooling\gs07_Z.eps',bits_per_pixel=8,color=1,/landscape,/encapsulated,/decomposed
    
    ;plot a few cooling functions vs Z for a few T
    plot,[0],[0],xtitle="Metallicity [Z/Zsun]",ytitle="Cooling Efficiency [erg cm"+textoidl('^3')+" s"+textoidl('^{-1}')+"]", $
         thick=1.1,charsize=1.5,/xlog,/ylog,xrange=[0.0005,3.0],/xstyle,yrange=[min(fc.Z1),max(fc.Z5)]
         
    for i=0,n_elements(temps)-1 do begin
      w = where(fc.T eq temps[i])
      L = [fc[w].Z1,fc[w].Z2,fc[w].Z3,fc[w].Z4,fc[w].Z5]
      oplot,zVals,L,line=i
      oplot,zVals,L,psym=4
      xyouts,zVals[1]*0.3,L[0]*0.8,string(temps[i],format="(e9.2)"),/data,alignment=1.0,charsize=1.0
    endfor
    
  device,/close_file
  set_plot,'WIN'
  
  stop
end

;Koyama & Inutsaka (2002) & Sanchez-Salcedo (2002)
; - analytic fit "reproduces relevant features" of Koyama & Inutsaka (2000)
; - plaw fit to Wolfire (95)
pro plotKI02SS02

  mass_hydrogen = 1.67e-24 ;g

  ;temp points
  tRes = 1000
  tMin = 1e1
  tMax = 1e4

  tPts    = findgen(tRes)/tRes * tMax + tMin
  tPtsLog = 10.0^( findgen(tRes)/tRes * alog10(tMax) + alog10(tMin) )

  plot,tPtsLog,KI02(tPtsLog),xtitle="Temperature [K]", ytitle="Cooling Efficiency [erg cm^3 / s]", $
       thick=1.1,charsize=1.5,/xlog,/ylog,/xstyle,xrange=[1e1,1e8],yrange=[10L^(-27.0),10L^(-16.0)]
       
  oplot,tPtsLog,(LambdaArrSS02(tPtsLog) * (mass_hydrogen)*(mass_hydrogen)),line=2
       
  stop
end

;Rosen & Bregman (1995)
;approximation to the radiative cooling function of Dalgarno & McCray (1972)
;and Raymond, Cox, & Smith (1976)
pro plotRB95

  ;temp points
  tRes = 1000
  tMin = 1e1
  tMax = 1e8

  tPts    = findgen(tRes)/tRes * tMax + tMin
  tPtsLog = 10.0^( findgen(tRes)/tRes * alog10(tMax) + alog10(tMin) )
  
  plot,tPtsLog,RB95arr(tPtsLog),xtitle="Temperature [K]", ytitle="Cooling Efficiency [erg cm^3 / s]", $
       thick=1.1,charsize=1.5,/xlog,/ylog,/xstyle,yrange=[10L^(-27.0),10L^(-20.0)]
  oplot,[1e4,1e4],[10L^(-28.0),10L^(-20.0)],line=2
  stop
end

;Dalgarno & McCray (1972)
;individual cooling functions by species (T<1e4 regime)
pro plotDC72

  ;abundances
  A_O  = 4.40*10.0^(-4.0)
  A_C  = 3.75*10.0^(-4.0)
  A_N  = 8.70*10.0^(-5.0)
  A_Si = 3.20*10.0^(-5.0)
  A_Fe = 3.20*10.0^(-5.0)

  ;temp points
  tRes = 1000
  tMin = 1e1
  tMax = 1e4

  tPts    = findgen(tRes)/tRes * tMax + tMin
  tPtsLog = 10.0^( findgen(tRes)/tRes * alog10(tMax) + alog10(tMin) )
  
  ;psopen,'twoslit.singlebulb.eps',/encapsulated,/landscape
    plot,tPtsLog,L_C1(tPtsLog)*A_C,xtitle="Temperature [K]", ytitle="Abundance * Cooling Efficiency [erg cm^3 / s]", $
         yrange=[10.0^(-27.0),10.0^(-22.0)],thick=1.1,charsize=1.5,/xlog,/ylog
         
    oplot,tPtsLog,L_Si1(tPtsLog)*A_Si,line=1,thick=1.1
    oplot,tPtsLog,L_Fe1(tPtsLog)*A_Fe,line=2,thick=1.1
    oplot,tPtsLog,(L_O1(tPtsLog)+L_O2(tPtsLog))*A_O,line=3,thick=1.1
  ;psclose

  stop
end