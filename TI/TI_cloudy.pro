; TI_cloudy.pro
; dnelson
; 22/9/10
;
; visualization for grids of Cloudy outputs

@helper

; plotSingleCloudyCooling()
pro plotSingleCloudyCooling, fBase, runID

  fileName = fBase + 'run' + str(runID) + '.dat'
  tMinMax  = [1e1,1e8]
  hcMinMax = [1e-30,1e-20]
  
  ; load
  ; ----
  
  headerLines = 15
  ptStruct = {T:0.0, heat:0.0, cool:0.0, mu:0.0}
  
  cOut = loadCSV(headerLines,fileName,ptStruct)
  
  ; plot
  ; ----

  ;set_plot,'PS'
  ;device,filename='cCool_'+runID+'.eps',bits_per_pixel=8,color=1,/landscape,/encapsulated,/decomposed


    plot,[0],[0],xtitle="Temperature [K]", $
         ytitle="Cooling Efficiency [erg cm"+textoidl('^3')+" s"+textoidl('^{-1}')+"]", $
         thick=1.1,charsize=1.5,/xlog,/ylog,/xstyle,/ystyle,/nodata,xrange=tMinMax, $
         yrange=hcMinMax
         
    oplot,cOut.T,cOut.cool,thick=3.0,color=fsc_color('orange')
    oplot,cOut.T,cOut.heat,thick=3.0,color=fsc_color('red')
    oplot,cOut.T,(cOut.cool - cOut.heat),thick=3.0,color=fsc_color('white')

  ;device,/close_file
  ;set_plot,'WIN'
  
  stop

end

; -----

pro temp1

  ;fBase   = "c:\zStuff\IDL.work\Default\cloudy\testism.1\cooling_run_"
  fBase   = "c:\zStuff\IDL.work\Default\cloudy\testism.2\cooling_run_"
  
  plotSingleCloudyCooling, fBase, 1

end