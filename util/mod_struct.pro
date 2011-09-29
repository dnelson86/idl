;+
; NAME:
;       MOD_STRUCT
;
; PURPOSE:
;       Modify an existing IDL structure: change, add or delete a tag
;
; CATEGORY:
;       General IDL tools.
;
; CALLING SEQUENCE:
;       newstr = mod_struct(oldstr, tagnam, tagval [, positi=0] [, /delete])
;
; REQUIRED INPUTS:
;       oldstr      the old structure to be modified
;       tagnam      name of the tag to be modified or added
;
; OPTIONAL INPUTS:
;       tagval      value of the tag to be modified or added
;
; KEYWORDS:
;       positi      position of the new or modified tag in the new struct
;       delete      delete the tag instead of adding or modifying it
;
; OUTPUTS:
;       newstr      a new structure with modified tags
;
; DESCRIPTION AND EXAMPLE:
;       This function modifies the IDL structure provided in OLDSTR by adding 
;       the new structure tag with the name TAGNAM and the value TAGVAL at 
;       the position specified by POSITI. If the tag already exists, then 
;       the old tag is overwritten, i.e. this function can be used to modify 
;       the data type of a tag or to move it within the structure. If no 
;       position is specified, then the tag is added at the end of the 
;       structure. For POSITI < 0 the tag is removed from the struct if 
;       previously present. The same result is obtained when the keyword 
;       DELETE is set. In this case TAGVAL does not have to be specified.
;
; CALLED BY:
;       midi_chpdat
;       midi_mskfit
;
; CALLING:
;       none
;
; MODIFICATION HISTORY:
;       2010-11-02 Written by Konrad R. W. Tristram
;       2010-11-04 K. Tristram: Return the input struct in case of input error.
;
FUNCTION MOD_STRUCT, OLDSTR, TAGNAM, TAGVAL, POSITI=POSITI, DELETE=DELETE

; CHECK IF INPUT PARAMETERS ARE OK
;-------------------------------------------------------------------------------
if (size(oldstr, /type) ne 8) or (size(tagnam, /type) ne 7) then begin
	message, 'Error in the parameters provided.', /continue
	message, 'Usage: newstr = mod_struct(oldstr, tagnam, tagval ' + $
	         '[, positi=0], [, /delete])', /continue
	message, '    oldstr - structure to be modified', /continue
	message, '    tagnam - name of the tag to be added', /continue
	message, '    help, tagval - value of the tag to be added', /continue
	message, '    positi - position of the new tag', /continue
	message, '    delete - delete the tag instead of adding it', /continue
	return, oldstr
endif

; INITIALISE THE COMMAND STRING AND THE SEPARATOR STRING
;-------------------------------------------------------------------------------
comman = 'newstr = create_struct('
sepstr = ''

; GET THE OLD TAG NAMES AND THE LOCATION OF THE SPECIFIED TAG
;-------------------------------------------------------------------------------
oldtag = tag_names(oldstr)
oldpos = where(oldtag eq strupcase(tagnam))

; IF THE TAG IS NEW THEN ADD IT AT THE END (UNLESS POSITI IS DEFINED)
;-------------------------------------------------------------------------------
if oldpos eq -1 then oldpos = n_tags(oldstr)

; SET THE NEW POSITION OF THE SPECIFIED TAG
;-------------------------------------------------------------------------------
if n_elements(positi) eq 0 then newpos = oldpos else newpos = positi
if keyword_set(delete) then newpos = -1

; LOOP OVER THE TAGS IN THE STRUCTURE
;-------------------------------------------------------------------------------
for i=0,n_tags(oldstr)-1 do begin
	; IF AT THE CORRECT POSITION, THEN INSERT THE NEW TAG
	;-----------------------------------------------------------------------
	if i eq newpos then begin
		comman += sepstr+'tagnam, tagval'
		sepstr = ', '
	endif
	; ADD ALL OTHER TAGS TO THE COMMAND
	;-----------------------------------------------------------------------
	if i ne oldpos then begin
		comman += sepstr+'oldtag[' +strtrim(i, 2)+'], ' + $
	                         'oldstr.('+strtrim(i, 2)+')'
		sepstr = ', '
	endif
endfor

; IF THE NEW TAG IS SUPPOSED TO BE AT THE END OF THE STRUCT THEN ADD IT NOW
;-------------------------------------------------------------------------------
if newpos ge i then comman += sepstr+'tagnam, tagval'

; FINISH THE COMMAND STRING AND EXECUTE THE COMMAND
;-------------------------------------------------------------------------------
comman += ')'
tmpvar = execute(comman)

; RETURN THE RESULT
;-------------------------------------------------------------------------------
return, newstr

END
