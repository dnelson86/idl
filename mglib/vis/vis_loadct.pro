; :Params:
;    table : in, optional, type=long
;       table number, 0-40 if using default color table file, 0-34 for Brewer
;       color tables, 0-6 for the Yorick/Gist color tables, or 0-15 for the
;       matplotlib color tables

;
; :Keywords:
;    file : in, optional, type=string, default=colors.tbl
;       filename of color table file; this is present to make `VIS_LOADCT`
;       completely implement `LOADCT`'s interface, it would normally not be used
;    brewer : in, optional, type=boolean
;       set to use the Brewer color tables
;    gmt : in, optional, type=boolean
;       set to use the GMT color tables
;    mpl : in, optional, type=boolean
;       set to use the matplotlib color tables
;    gist : in, optional, type=boolean
;       set to use the Gist/Yorick color tables
;    chaco : in, optional, type=boolean
;       set to use the Chaco color tables
;    vis : in, optional, type=boolean
;       set to use the VIS library color tables
;    rgb_table : out, optional, type="lonarr(ncolors, 2)"
;       set to a named variable to retrieve the color table
;    reverse : in, optional, type=boolean
;       set to reverse color table
;    show_tables : in, optional, type=boolean
;       set to print a listing of the available color tables
;    cpt_filename : in, optional, type=string
;       filename of `.cpt` file to load a color table from; the `.cpt` 
;       filename extension is optional; the filename given can be absolute, 
;       relative from the current working directory, or relative from the 
;       `cpt-city` directory in the VIS library; see `cptcity_catalog.idldoc`
;       for a listing of the `.cpt` files provided with the VIS library
;    _ref_extra : in, out, optional, type=keyword
;       keywords to `LOADCT`
;-
function vis_src_root
  return,'/n/home07/dnelson/idl/mglib/vis/'
end

pro vis_loadct, table, file=file, $
                brewer=brewer, gmt=gmt, mpl=mpl, gist=gist, chaco=chaco, vis=vis, $
                rgb_table=rgbTable, $
                reverse=reverse, $
                show_tables=showtables, $
                cpt_filename=cptFilename, $
                _ref_extra=e
  compile_opt strictarr
  on_error, 2
  
  if (n_elements(cptFilename) gt 0L) then begin
    if (strmid(cptFilename, 3, /reverse_offset) ne '.cpt') then begin
      _cptFilename = cptFilename + '.cpt'
    endif else _cptFilename = cptFilename
    
    if (~file_test(_cptFilename)) then begin
      _cptFilename = filepath(_cptFilename, $
                              subdir=['cpt-city'], $
                              root=vis_src_root())
    endif
    
    if (~file_test(_cptFilename)) then begin
      message, '.cpt file not found, ' + cptFilename
    endif
    
    rgb = vis_cpt2ct(_cptFilename, name=ctnames)
    
    if (keyword_set(showTables)) then begin
      print, ctnames
      
      return
    endif
    
    if (keyword_set(reverse)) then rgb = reverse(rgb, 1)
    
    if (arg_present(rgbTable)) then begin
      rgbTable = rgb
    endif else begin
      tvlct, rgb
    endelse
    
    return
  endif
  
  case 1 of
    keyword_set(brewer): ctfilename = filepath('brewer.tbl', root=vis_src_root())
    keyword_set(gmt): ctfilename = filepath('gmt.tbl', root=vis_src_root())
    keyword_set(mpl): ctfilename = filepath('mpl.tbl', root=vis_src_root())
    keyword_set(gist): ctfilename = filepath('gist.tbl', root=vis_src_root())
    keyword_set(chaco): ctfilename = filepath('chaco.tbl', root=vis_src_root())
    keyword_set(vis): ctfilename = filepath('vis.tbl', root=vis_src_root())
    n_elements(file) gt 0L: ctfilename = file
    else:
  endcase
  
  if (keyword_set(showTables)) then begin
    loadct, get_names=ctnames, file=ctfilename
    print, ctnames
    return
  endif
  
  ; search for color table name if it is specified as a string
  if (size(table, /type) eq 7L) then begin
    loadct, get_names=ctnames, file=ctfilename
    ind = where(stregex(ctnames, table, /boolean, /fold_case), count)
    if (count gt 0L) then begin
      _table = ind[0]
    endif else begin
      message, 'Error: Color table not found'
      return
    endelse
  endif else begin
    if (n_elements(table) gt 0L) then _table = table
  endelse
  
  loadct, _table, rgb_table=rgbTable, file=ctfilename, _strict_extra=e
  
  if (keyword_set(reverse)) then begin
    rgbTable = reverse(rgbTable, 1)
  endif 
  
  if (~arg_present(rgbTable)) then tvlct, rgbTable
end

