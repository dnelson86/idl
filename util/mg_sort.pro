; docformat = 'rst'

;+
; An alternative to IDL's SORT function that employs a 
; `radix-sort algorithm <http://en.wikipedia.org/wiki/Radix_sort>`. This 
; algorithm is stable but slower than the 
; `quicksort algorithm <http://en.wikipedia.org/wiki/Quicksort>` used in IDL's
; `SORT` function.
;
; :History:
; 	Derived from code by Atle Borsholm, VIS, 2012
;-


;+
; Handle special case of sorting string arrays.
;
; Note that this is very slow compared to `SORT` in most cases. However, it is 
; a stable sort.
;
; :Returns:
;    `lonarr(n_elements(data))`
;
; :Params:
;    data : in, required, type=strarr
;       array of strings to sort
;-
function mg_sort_strings, data
  compile_opt idl2, logical_predicate

  len = strlen(data)
  maxlen = max(len)
  sorted = data
  radix = 256 ; fixed to byte for each character
  indices = lindgen(n_elements(data))
  
  for i = maxlen - 1L, 0L, -1L do begin
    ; use a nice feature of IDL: byte('') = [0b]
    digit = byte(strmid(sorted, i, 1L))
    h = histogram(digit, reverse_indices=ri)
    sorted = sorted[ri[radix + 1L:*]]
    indices = indices[ri[radix + 1L:*]]
  endfor
  
  return, indices
end


;+
; Radix sort input data.
;
; :Returns:
;    `lonarr(n_elements(data))`
;
; :Params:
;    data : in, required, type="int, long, float, double, or string array"
;       data to sort
;
; :Keywords:
;    radix : in, optional, type=int, default=256
;       radix to use for sort
;-
function mg_sort, data, radix=radix
  compile_opt idl2, logical_predicate
   
   ; default radix if not specified
   _radix = n_elements(radix) eq 0L ? 256 : radix
     
  ; support signed ints
  case size(data, /type) of
    2: sorted = uint(data) + ishft(1us, 15) ; +/- are equivalent here
    3: sorted = ulong(data) + ishft(1ul, 31)
    4: begin
         sorted = ulong(data, 0, n_elements(data))
         sorted xor= ((sorted and ishft(1ul, 31)) ne 0) * (not ishft(1ul, 31))
         sorted += ishft(1ul, 31)
       end
    5: begin
         sorted = ulong64(data, 0, n_elements(data))
         ; IEEE does not use 2's complement for negative numbers so, flip the 
         ; bits for the negative numbers with this trick this may not be 100% 
         ; IEEE compliant sorting, but works well for finite numbers at least
         sorted xor= ((sorted and ishft(1ull, 63)) ne 0) * (not ishft(1ull, 63))
         sorted += ishft(1ull, 63)
       end
    7: return, mg_sort_strings(data)
    14: sorted = ulong64(data) + ishft(1ull, 63)
    else: sorted = data
  endcase

  indices = lindgen(n_elements(data))   
  mx = max(sorted, min=mn)
   
  ; implement a slight speed improvement for small data ranges
  rng = mx - mn
  if (rng lt mn / _radix) then begin
    sorted -= mn
    mx -= mn
  endif
   
  factor = 1ull
  while (mx gt 0) do begin
    mx /= _radix
    rem = sorted / factor
    digit = rem mod _radix
    factor = factor * _radix
    h = histogram(digit, min=0, max=_radix-1, binsize=1, reverse_indices=ri)
    ind = ri[_radix + 1:*]
    sorted = sorted[ind]
    indices = indices[ind]
  endwhile
   
  return, indices
end
