pro print_to_width, str
;+
; Name
; ----
;   PRINT_TO_WIDTH
;
; Purpose
; -------
;   Prints text to the terminal or IDL workbench. Strings that
;   are longer than the width of the window will be wrapped to
;   the next line without splitting words. This is different
;   from IDL's ``print`` procedure, which will split words.
;
; Calling Sequence
; ----------------
;   ::
;
;       print_to_width, str
;
; Input
; -----
;   ``str`` : string scalar
;       The string that is to be printed.
;
; Output
; ------
;   The input string is printed to the screen, without splitting words.
;
; Modification History
; --------------------
;   - 2022/05/12: Created (Keith Doore)
;-
 On_error, 2
 Compile_opt idl2

; Error handling
 if n_elements(str) eq 0 then message, 'Variable is undefined: STR.'
 if size(str, /type) ne 7 then message, 'STR must be of type string.'
 if size(str, /n_dim) ne 0 then message, 'STR must be a scalar.'


; Get the terminal width
 term_width = (TERMINAL_SIZE())[0]


; Check if str is longer than the width
 if strlen(str) gt term_width then begin
 ; If it is longer, split into lines that are just less than the width
   words = STRSPLIT(str, /EXTRACT)+' '
   substr = ''
   substr_arr = !null
   for i=0, n_elements(words)-1 do begin
     if (strlen(substr)+strlen(words[i])) gt term_width then begin
       substr_arr = [substr_arr, substr]
       substr = ''
     endif
     substr = substr+words[i]
     if i eq n_elements(words)-1 then substr_arr = [substr_arr, substr]
   endfor

   print, substr_arr
 endif else begin
   ; If not, print to screen
   print, str
 endelse

end