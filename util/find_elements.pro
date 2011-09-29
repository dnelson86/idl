;   FIND_ELEMENTS.PRO   5-01-95
;
; Function to return the array indices for the elements of one array
; that are contained in another array (see Limitations below!). If
; a third argument is supplied then these elements are adjusted to
; the values contained in the third argument, and returned in the
; first argument.
;
; Limitation:  It is assumed that there are a relatively small number
;   of elements in the second (and third) arrays...since each element
;   must be processed separately, this routine will be very slow if
;   these arrays have many elements. 
;
; Error returns:
;
;   Returns -2 if one of the arguments is inappropriate, and -1 if
;   none of the values in FIND are found in ARRAY.

FUNCTION find_elements, array, find, adjust

if (n_params() gt 2) then begin
	if (n_elements(find) ne n_elements(adjust)) then begin
		message, 'Array arguments two and three must be the same size', $
				/continue
		return, -2
	endif
	array_copy = array(*)                          ; Use copy of original array
	adjust_array = 1 
endif else begin
	adjust_array = 0
endelse

flag = 0                                           ; Flag elements found

for i = 0, n_elements(find) - 1 do begin
	found = where(array(*) eq find(i))
	if (found(0) ne -1) then begin
		if (adjust_array eq 1) then array_copy(found) = adjust(i)
		if (flag eq 0) then found_total = found $
			else found_total = [found_total, found]
		flag = 1
	endif
endfor

if (adjust_array eq 1 and flag eq 1) then $        ; Adjust elements of Array
	array(found_total) = array_copy(found_total)   ;  if 3rd arg specified

if (n_elements(found_total) gt 0) then $
	return, found_total $
else $
	return, -1
END
