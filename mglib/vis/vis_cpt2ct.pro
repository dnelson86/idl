; docformat = 'rst'

;+
; Converts a GMT color table file to an RGB color table.
;
; :Examples:
;    For example, if the GMT color tables files were stored in the `cpt` 
;    directory::
;
;       IDL> tvlct, vis_cpt2ct('cpt/GMT_relief')
; 
; :Returns:
;    bytarr(256, 3)
;
; :Params:
;    filename : in, required, type=string
;       filename of GMT color table file, i.e., .cpt file
;
; :Keywords:
;    name : out, optional, type=string
;       color table name
;-
function vis_cpt2ct, filename, name=name
  compile_opt strictarr

  nlines = file_lines(filename)
  lines = strarr(nlines)
  openr, lun, filename, /get_lun
  readf, lun, lines
  free_lun, lun
  
  colorModel = 'RGB'
  matches = ''
  i = 0
  while (matches[0] eq '' && strmid(lines[i], 0, 1) eq '#') do begin
    matches = stregex(lines[i++], '^# COLOR_MODEL = ([+HSVRGB]*)', /subexpr, /extract)
    if (matches[0] ne '') then colorModel = matches[1]
  endwhile
  
  dataLines = stregex(lines, '^[[:space:]]*[-[:digit:]]', /boolean)
  dataLinesInd = where(dataLines eq 1B, count)

  data = fltarr(8, count)
  for i = 0L, count - 1L do begin
    data[*, i] = float((strsplit(lines[dataLinesInd[i]], /extract))[0:7])
  endfor
  
  minValue = data[0, 0]
  maxValue = data[4, count - 1]
  
  cutoffs = value_locate([reform(data[0, *]), data[4, count - 1L]], $
                         (maxValue - minValue) * findgen(256) / 255. + minValue)
  
  ncolors = histogram(cutoffs < (count - 1L))

  result = bytarr(256, 3)
  pos = 0L
  for i = 0L, count - 1L do begin
    c1 = reform(data[1:3, i])
    c2 = reform(data[5:7, i])

    if (ncolors[i] gt 0) then begin
      colors = congrid(transpose([[c1], [c2]]), ncolors[i], 3, /interp, /minus_one)

      if (colorModel eq 'HSV' || colorModel eq '+HSV') then begin
        color_convert, colors[*, 0], colors[*, 1], colors[*, 2], r, g, b, /hsv_rgb
        colors[*, 0] = r
        colors[*, 1] = g
        colors[*, 2] = b
      endif
    
      result[pos, 0] = colors
    endif
    
    pos += ncolors[i]
  endfor
  
  name = stregex(file_basename(filename), 'GMT_([_[:alnum:]]*).cpt', /subexpr, /extract)
  if (name[0] eq '') then name = stregex(file_basename(filename), '([_[:alnum:]]*).cpt', /subexpr, /extract)
  name = name[1]
 
  return, result
end

end
