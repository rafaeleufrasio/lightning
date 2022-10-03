;+
; Name
; ----
;   LIGHTNING_BUILD
;
; Purpose
; -------
;   Generates a build file containing all of the compiled 
;   Lightning functions and procedures. Additionally, the
;   lightning_dir system variable is added to the user's startup
;   file, or one is created if it does not exist.
;
; Calling Sequence
; ----------------
;   ``@lightning_build``
;
; Output
; ------
;   A save file (``lightning.sav``) containing all of the compiled
;   functions and procedures required to run Lightning.
;
; Modification History
; --------------------
;   - 2022/03/17: Created (Erik B. Monson)
;   - 2022/05/23: Added lightning_dir system variable to startup file (Keith Doore)
;   - 2022/06/24: Added documentation (Keith Doore)
;   - 2022/06/27: Copied default lightning_configuration.pro to top directory for user (Keith Doore)
;   - 2022/08/11: Updated path for lightning_compile_all (Keith Doore)
;   - 2022/08/23: Fix issue where a blank line would be added to top of existing startup file (Keith Doore)
;-
 On_error, 2
 compile_opt idl2


; Get lightning directory from current directory, since
;  one needs to be in the lightning directory to build it.
 cd, current=lightning_dir
 ; Include directory marker, since cd does not
 lightning_dir = lightning_dir+'/'

; Get startup file
 start_file = pref_get('IDL_STARTUP')

 file_array = !null
 ; If startup file exists, read each line into file_array
 line = ''
 if start_file ne '' then begin $
   openr, lun, start_file, /get_lun &$
   while ~eof(lun) do begin $
     readf, lun, line &$
     file_array = [file_array, line]
 ; Else set to blank string
 if n_elements(file_array) eq 0 then file_array = '' else free_lun, lun

 ; Check if line defining !lightning_dir system variable is already present
 bool = stregex(file_array, "^defsysv, ['"+'"]'+"!lightning_dir['"+'"]'+", .*, 1$", /boolean)

 ; If not present, add line to file_array as final entry (line)
 if total(bool) eq 0 then file_array=[file_array, '', "defsysv, '!lightning_dir', '"+lightning_dir+"', 1"]

 ; If no startup file exists, create one
 if start_file eq '' then start_file = file_search('~/', /expand_tilde)+'/.idl/startup.pro'
 ; Write edited or created startup file
 openw, lun, start_file, /get_lun
 for i=0,n_elements(file_array)-1 do printf, lun, file_array[i]
 free_lun, lun
 
 pref_set, 'IDL_STARTUP', start_file, /commit


; Copy lighting_configuration.pro to top level directory
 file_copy, lightning_dir+'pro/configuration/lightning_configure_defaults.pro', $
            lightning_dir+'lightning_configure.pro',  /overwrite


; Start from scratch before compiling
 .reset


; Compile functions and procedures
 @pro/build/lightning_compile_all


; Save the compiled functions and procedures to `lightning.sav`
 save, filename='lightning.sav', /routines, description='Lightning'

 exit