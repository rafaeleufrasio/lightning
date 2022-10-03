pro lightning, input_file, resume=resume, interactive=interactive
;+
; Name
; ----
;   LIGHTNING
;
; Purpose
; -------
;   Calls lower level Lightning procedures and functions to
;
;   1) generate, check, and save the configuration file;
;   2) load, check, and save the SED data;
;   3) call the Lightning fitting procedure to fit each SED using the specified configuration;
;   4) and post-process the fitting results.
;
;   Each call to the Lightning fitting procedure is optimized to run
;   in parallel depending on the allowed CPU usage.
;
; Calling Sequence
; ----------------
;   ::
;
;       lightning, input_file [, /resume, /interactive]
;
; Input
; -----
;   ``input_file`` : string scalar
;       The name (including path) to the file containing the SED fluxes
;       and distances (or redshifts) in a data table. If the file is a
;       FITS file, the table must be in the first extension. (See
;       :ref:`input-formats-label` for details, required contents and
;       format of the data table.)
;
; Optional Inputs
; ---------------
;   ``resume`` : flag
;       If set, then Lightning will resume running where it left off assuming
;       some issue caused it to stop mid-run. This means Lightning would check
;       what SEDs have been fit and would only fit those that were not completed.
;   ``interactive`` : flag
;       If set, then Lightning will give prompts to the command line allowing
;       for a user to set the Lightning configuration interactively, rather
;       than having to manually edit the configuration file.
;
; Modification History
; --------------------
;   - 2022/05/04: Created (Keith Doore)
;   - 2022/05/18: Rearranged input and configuration to allow for cosmology in input (Keith Doore)
;   - 2022/06/06: Moved tabulated prior check to its own procedure (Keith Doore)
;   - 2022/08/01: Added ``/silent`` to ``mrdfits`` (Keith Doore)
;   - 2022/08/02: Added ASCII input to FITS file conversion (Keith Doore)
;   - 2022/08/18: Added progress printing (Keith Doore)
;   - 2022/09/01: Swapped ``input_dir`` to absolute path, and passed ``config`` to ``lightning_input`` (Erik B. Monson)
;   - 2022/09/14: Updates to allow fitting with X-ray fluxes (Erik B. Monson)
;   - 2022/09/16: Added check for syntax errors in user edited ``lightning_configure.pro`` (Keith Doore)
;   - 2022/09/26: Fixed ``input_dir`` to absolute path bug (Keith Doore)
;-
 On_error, 2
 compile_opt idl2


; Check if using an IDL virtual machine and allow for required inputs
 vm = lmgr(/vm)
 ; If an IDL virtual machine is used, only interactive mode can be used
 if vm then interactive = 1
 ; Input file
 if vm then begin
   error = 1
   input_message = 'Please specify the file (including path) to the FITS file '+$
                   'containing the input data table (e.g., /<file_path>/<file_name>.fits).'
   print_to_width, '======================='
   while error do begin
     print_to_width, ''
     print_to_width, input_message

     temp = ''
     read, temp, prompt='INPUT_FILE: ', form='(A0)'
     temp = strtrim(temp, 2)

     if ~(file_test(temp)) then begin
       input_message = 'Given file does not exist. Please confirm the file and path are '+$
                       'correct and try again.'
     endif else error = 0
   endwhile
   input_file = temp
 endif
 ; Allow for resuming
 if vm then begin
   error = 1
   resume_message = 'Would you like Lightning to resume running where it left '+$
                    'off, assuming some issue caused it to stop mid-run? Please '+$
                    'indicate yes or no. (Default: no)'
   print_to_width, '======================='
   while error do begin
     print_to_width, ''
     print_to_width, resume_message

     temp = ''
     read,temp, prompt='RESUME: ',form='(A0)'
     temp = strtrim(temp, 2)
     if temp eq '' then temp = 'NO'

     if ~(strupcase(temp) eq 'YES' or strupcase(temp) eq 'NO') then begin
       resume_message = 'Please type yes or no.'
     endif else error = 0
   endwhile
   if strupcase(temp) eq 'YES' then temp = 1 else temp = 0
   resume = temp
 endif
 ; Lightning directory
 if vm then begin
   ; If using a virtual machine, create system variable for Lighting directory,
   ;   since startup files are not executed. The current directory is the directory
   ;   containing the .sav file (i.e., the Lightning directory).
   cd, current=lightning_dir
   ; Include directory marker, since cd does not
   defsysv, '!lightning_dir', lightning_dir+'/', 1
 endif

; Get time for progress printing
 t0 = systime(/sec)

; Error check input_file
 if n_elements(input_file) eq 0 then message, 'Variable is undefined: INPUT_FILE.'
 if size(input_file, /type) ne 7 then message, 'INPUT_FILE must be of type string.'
 if size(input_file, /n_dim) ne 0 then message, 'INPUT_FILE must be a scalar.'
 if ~file_test(input_file) then message, 'No input data file exists. '+$
                                        'Please confirm specified file name leads to the correct file.'

; Check Lightning directory
 if n_elements(!lightning_dir) eq 0 then message, $
                  '!LIGHTNING_DIR system variable is not defined. '+$
   'Please confirm Lightning package has been built and try again.'
 if size(!lightning_dir, /type) ne 7 then message, $
                     '!LIGHTNING_DIR system variable is not of type string. '+$
   'Please confirm Lightning package has been built correctly and try again.'
 if size(!lightning_dir, /n_dim) ne 0 then message, $
                          '!LIGHTNING_DIR system variable must be a scalar. '+$
   'Please confirm Lightning package has been built correctly and try again.'
 if ~file_test(!lightning_dir+'lightning.sav') then message, $
     '!LIGHTNING_DIR system variable does not lead to Lightning build file. '+$
   'Please confirm Lightning package has been built correctly and try again.'



; Restore the compiled functions and procedures,
;   and generate constants system variable
 restore, !lightning_dir+'lightning.sav'
 lightning_constants


; Check if a configuration .sav or .pro file exists in input file directory.
;   If not, create a configuration structure.
 tconfig0 = systime(/sec)
 input_dir = file_search(file_dirname(input_file), /fully_qualify_path) + '/'
 input_dir = input_dir[0]
 config_sav_exists = file_test(input_dir+'lightning_output/lightning_configure.sav')
 config_pro_exists = file_test(input_dir+'lightning_configure.pro')

 if keyword_set(resume) then begin
   if ~config_sav_exists then message, 'The Lightning configuration .sav file was not '+$
       'found within the lightning_output subdirectory of the input file directory. '+$
       'Please do not specify the RESUME keyword unless the .sav file exists.'
   restore, input_dir+'lightning_output/lightning_configure.sav'
 endif else if keyword_set(interactive) then begin
   if config_sav_exists then begin
     restore, input_dir+'lightning_output/lightning_configure.sav'
     config = lightning_configure_interactive(config_edit=config)
   endif else begin
     config = lightning_configure_interactive()
   endelse
 endif else begin
   if config_pro_exists then begin
    ; Recompile the configuration structure procedure file, since it was potentially
    ;   edited from the default
     cd, current=cwd
     cd, input_dir
     resolve_routine, 'lightning_configure', /is_function
     if strupcase(strtrim(!error_state.msg, 2)) eq 'SYNTAX ERROR.' then begin
       cd, cwd
       message, 'A syntax error was found in the edited lightning_configure.pro file. '+$
                'See above error message for general location of syntax error in the file.'
     endif
     cd, cwd
   endif else begin
     message, 'Using the default Lightning configuration.', /informational
   endelse
   config = lightning_configure()
 endelse

; Check the configuration structure for errors and save it to the input file directory
 config = lightning_configure_check(config)
 file_mkdir, input_dir+'lightning_output/'
 save, config, filename=input_dir+'lightning_output/lightning_configure.sav'
 ; Copy the default configuration structure procedure file to the lightning_output directory
 if ~config_pro_exists then begin
   file_copy, !lightning_dir+'pro/configuration/lightning_configure_defaults.pro', $
              input_dir+'lightning_configure.pro'
 endif
 if config.PRINT_PROGRESS then print, 'Configuration free of errors. Scan complete in '+$
                                      string(systime(/sec)-tconfig0, f='(F0.2)')+' s. Current total run time: '+$
                                      string(systime(/sec)-t0, f='(F0.2)')+' s.'


; Read in the input file and save each SED to a .sav file.
 ; Check if file is an ascii or fits file. If not ascii, convert to fits.
 ;   Easier to check for FITS extension rather than plethora of ascii table extensions.
 tinput0 = systime(/sec)
 void = where(strupcase(strsplit(file_basename(input_file), '.', /extract)) eq 'FITS', Nfits, /null)
 if Nfits gt 0 then input_file_fits = input_file else ascii_input_to_fits, input_file, input_file_fits=input_file_fits

 lightning_input, input_file_fits, _EXTRA=config
 sed_files = file_search(input_dir+'lightning_output/input_sav_files/lightning_input_*.sav')
 sed_id = strmid(file_basename(sed_files, '.sav'), strlen('lightning_input_'))

 ; If resuming, remove the fit SED id and file from list.
 if keyword_set(resume) then begin
   fit_sed_files = file_search(input_dir+'lightning_output/output_sav_files/lightning_output_*.sav', $
                               count=nout_files)
   ; Extract SED IDs from output files
   if nout_files gt 0 then begin
     fit_sed_id = strmid(file_basename(fit_sed_files, '.sav'), strlen('lightning_output_'))

     incompleted = intarr(n_elements(sed_id))
     for i=0, n_elements(sed_id)-1 do if total(sed_id[i] eq fit_sed_id) eq 0 then $
       incompleted[i] = 1
     sed_id = sed_id[where(incompleted, /null)]
     sed_files = sed_files[where(incompleted, /null)]
     if n_elements(sed_files) eq 0 then $
       message, 'All SEDs have already been successfully fit. Please do not specify the '+$
                'RESUME keyword unless there are some unfit SEDs.'
   endif
 endif
 if config.PRINT_PROGRESS then print, 'Input file read in and free of errors. Process complete in '+$
                                      string(systime(/sec)-tinput0, f='(F0.2)')+' s. Current total run time: '+$
                                      string(systime(/sec)-t0, f='(F0.2)')+' s.'


; Check if xray-emission is set in config, and if set, then Xray data must exist
 if config.XRAY_EMISSION then begin
   ; Only need to check one file as all SED data files will or will not have xray data
   restore, sed_files[0]
   xray_tag_count = total(tag_names(sed_data) eq 'NET_COUNTS') + total(tag_names(sed_data) eq 'XRAY_LNU_OBS')
   if (xray_tag_count eq 0) then $
     message, 'The X-ray emission model was selected to be used in the Lightning configuration. '+$
       'However, no X-ray data was found in the input file. '+$
       'Please either include X-ray data in the input file or deselect the X-ray emission model.'
 endif



; Check if any priors are tabulated and confirm that they have the correct format for Lightning.
;   This is outside lightning_configure_check because each SED can have unique tabulated priors.
 tabulated_prior_check, config, sed_id



; Check if running the Doore21 attenuation. If so, check if the correct Lbol_abs tables are generated.
;   If not, generate them, since all SEDs will share the same SFH bins.
 tdooretable0 = systime(/sec)
 if strupcase(config.ATTEN_CURVE) eq 'DOORE21' then begin
   data = mrdfits(input_file_fits, 1, /silent)
   data_tags = tag_names(data)
   if total(strupcase(data_tags) eq 'REDSHIFT') eq 1 then begin
     redshifts = data.redshift
     ; Converting to string assures same format as extracted redshifts from file names
     redshifts = string(redshifts[uniq(redshifts, sort(redshifts))], f='(f0.6)')
   endif else redshifts = '0.000000'
   data = !null

   Lbol_abs_tables_exist = file_test(input_dir+'lightning_output/doore21_Lbol_abs_table/redshift_0.000000.fits')
   if Lbol_abs_tables_exist then begin
     Lbol_abs_table_files = file_search(input_dir+'lightning_output/doore21_Lbol_abs_table/redshift_*.fits')
     Lbol_abs_redshifts = strmid(file_basename(Lbol_abs_table_files, '.fits'), strlen('redshift_'))

     missing = intarr(n_elements(redshifts))
     for i=0, n_elements(redshifts)-1 do if total(redshifts[i] eq Lbol_abs_redshifts) eq 0 then $
       missing[i] = 1

     Lbol_abs_table0 = mrdfits(input_dir+'lightning_output/doore21_Lbol_abs_table/redshift_0.000000.fits', 1, /silent)
     steps_bounds = config.STEPS_BOUNDS < galage(0.d, 1e3, /silent, _extra=config.COSMOLOGY)
     steps_bounds = steps_bounds[uniq(steps_bounds, sort(steps_bounds))]

     incorrect_missing = ~((n_elements(Lbol_abs_table0.STEPS_BOUNDS) eq n_elements(steps_bounds)) and $
                           (total(Lbol_abs_table0.STEPS_BOUNDS eq steps_bounds) eq n_elements(steps_bounds))) or $
                         total(missing) ne 0
     if incorrect_missing then begin
       if config.PRINT_PROGRESS then print, 'Generating Doore+2021 Lbol_abs tables.'
       generate_doore21_lbol_abs_table, input_dir, double(redshifts), config
       if config.PRINT_PROGRESS then print, 'Table generation complete in '+string(systime(/sec)-tdooretable0, f='(F0.2)')+$
                                            ' s. Current total run time: '+string(systime(/sec)-t0, f='(F0.2)')+' s.'
     endif
   endif else begin
     if config.PRINT_PROGRESS then print, 'Generating Doore+2021 Lbol_abs tables.'
     generate_doore21_lbol_abs_table, input_dir, double(redshifts), config
     if config.PRINT_PROGRESS then print, 'Table generation complete in '+string(systime(/sec)-tdooretable0, f='(F0.2)')+$
                                          ' s. Current total run time: '+string(systime(/sec)-t0, f='(F0.2)')+' s.'
   endelse
 endif


; Make the output directory for the fitting results
 file_mkdir, input_dir+'lightning_output/output_sav_files/'


; Run Lightning in parallel utilizing the IDL_IDL bridge, where the number of bridges
;   is given by the number of specified number of CPUs
 ; If config.MAX_CPUS exceeds the number of CPUs on the machine, then all CPUs will be used.
 num_bridges = config.MAX_CPUS < !CPU.HW_NCPU

 ; Virtual machine cannot run IDL_IDL bridge. So, set number of bridges to 1
 if vm then num_bridges = 1

 ; Determine the number of groups of bridges that will need to be run
 if (n_elements(sed_files) mod float(num_bridges)) ne 0 then num_groups = fix(n_elements(sed_files)/float(num_bridges))+1
 if (n_elements(sed_files) mod float(num_bridges)) eq 0 then num_groups = fix(n_elements(sed_files)/float(num_bridges))

 if config.PRINT_PROGRESS then print, 'Beginning SED fitting. Current total run time: '+$
                                      string(systime(/sec)-t0, f='(F0.2)')+' s.'

 tfit0 = systime(/sec)
 for j=0,(num_groups-1) do begin
   ; If a virtual machine (1 CPU) is used, then do not (no need to) use IDL_IDL bridge
   if num_bridges eq 1 then begin
     tfit = systime(/sec)

     restore, sed_files[j]
     lightning_fit, input_dir, sed_data, config

     if config.PRINT_PROGRESS then begin
       tnow = systime(/sec)
       t_elapsed =  tnow - tfit0
       avg_speed = (j + 1) / t_elapsed
       t_remaining = (num_groups - (j + 1)) / avg_speed
       print, 'SED '+string(j+1, f='(I0)')+' of '+string(num_groups, f='(I0)')+' fit in '+$
              string(tnow-tfit, f='(F0.2)')+' s. Estimated fitting time remaining: '+$
              string(t_remaining, f='(F0.1)')+' s.'
     endif
   endif else begin

     ;Last group of bridges may be smaller group. Set to include this possibility
     if j eq (num_groups-1) then num_bridges_temp = n_elements(sed_files)-(j*num_bridges) else $
                                 num_bridges_temp = num_bridges

     ;Below for loop creates and executes a group of bridges
     bridge = objarr(num_bridges_temp)
     for i=0, (num_bridges_temp-1) do begin
       bridge[i]=IDL_IDLbridge()
       ; The child process will not have access to data, compiled routines, or system variables.
       bridge[i]->EXECUTE, "defsysv, '!lightning_dir', '"+!lightning_dir+"', 1"
       bridge[i]->EXECUTE, "restore, '"+!lightning_dir+"lightning.sav'"
       bridge[i]->EXECUTE, 'lightning_constants'
       bridge[i]->SetVar,  'input_dir', input_dir
       bridge[i]->EXECUTE, "restore, '"+input_dir+"lightning_output/lightning_configure.sav'"
       bridge[i]->EXECUTE, "restore, '"+sed_files[(i+(j*num_bridges))]+"'"
       bridge[i]->EXECUTE, 'lightning_fit, input_dir, sed_data, config', /nowait
     endfor

     ;Monitor progress of bridges. When status is all 0 then all bridges have finished running.
     status = replicate(1, num_bridges_temp)
     while max(status) gt 0 do begin
       wait, 1
       for i=0,(num_bridges_temp-1) do begin
         status[i] = bridge[i]->Status()
       endfor
     endwhile

     ;Destroy all bridges once all are complete so next set runs on new bridges and they do not accumulate.
     for i=0,(num_bridges_temp-1) do begin
       bridge[i]->Cleanup
     endfor
     ;Reset bridge in case number of bridges changes
     bridge=!null

     ; If progress printing, print progress bar.
     if config.PRINT_PROGRESS then lightning_print_progress, j, num_groups, tfit0, funcname='LIGHTNING_FIT'
   endelse
 endfor



; ======= Post processing =========
 if config.PRINT_PROGRESS then print, 'Fitting complete in '+string(systime(/sec)-tfit0, f='(F0.2)')+$
                                      ' s. Beginning post-processing. Current total run time: '+$
                                      string(systime(/sec)-t0, f='(F0.2)')+' s.'

 ; Post-process the SED fits for all galaxies not just those if resuming
 sed_files = file_search(input_dir+'lightning_output/input_sav_files/lightning_input_*.sav')
 sed_id = strmid(file_basename(sed_files, '.sav'), strlen('lightning_input_'))

 tpostprocess0 = systime(/sec)
 lightning_postprocessing, input_dir, config, sed_id

 if config.PRINT_PROGRESS then print, 'Post-processing complete in '+string(systime(/sec)-tpostprocess0, f='(F0.2)')+$
                                      ' s. Total Lightning run time: '+string(systime(/sec)-t0, f='(F0.2)')+' s.'


end