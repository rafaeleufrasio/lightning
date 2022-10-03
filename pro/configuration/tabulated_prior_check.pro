pro tabulated_prior_check, config, sed_id
;+
; Name
; ----
;   TABULATED_PRIOR_CHECK
;
; Purpose
; -------
;   Checks that all tabulated priors specified in the Lightning 
;   configuration structure are properly formatted. The correct format
;   requires that each specified tabulated prior is found within the 
;   tabulated prior files for each SED, has gridded values in ascending 
;   order, and non-negative probabilities for each grid value.
;
; Calling Sequence
; ----------------
;   ::
;
;       tabulated_prior_check, config, sed_id
;
; Inputs
; ------
;   ``config`` : structure
;       A Lightning configuration structure. (See
;       ``lightning_configure_defaults.pro`` for details and contents.)
;   ``sed_id`` : string array(Nsed)
;       The ID of each SED.
;
; Output
; ------
;   An error message stating if any errors are found within the tabulated prior structures.
;   Error messages are relatively specific and should help users resolve any issues.
;
; Modification History
; --------------------
;   - 2022/06/06: Created (Keith Doore)
;   - 2022/08/01: Added ``/silent`` to ``mrdfits`` (Keith Doore)
;   - 2022/08/02: Added check to ensure initialization range is within tabulated range (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(config) eq 0 then message, 'Variable is not defined: CONFIG.'
 if size(config, /type) ne 8 then message, 'CONFIG is not of type structure.'

 if n_elements(sed_id) eq 0 then message, 'Variable is not defined: SED_ID.'
 if size(sed_id, /type) ne 7 then message, 'SED_ID is not of type string.'
 if size(sed_id, /n_dim) gt 1 then message, 'SED_ID must be a scalar or 1D array.'


; Check tabulated priors
; Make list of all priors that are tabulated and the specified directories
 tabulated_priors = !null
 init_range = !null
 config_hash = orderedhash(config, /extract)
 foreach prior_substruct, config_hash, param_name do begin
   if size(prior_substruct, /type) eq 11 then begin
     if prior_substruct.haskey('PRIOR') then begin
       if total(strupcase(prior_substruct['PRIOR']) eq 'TABULATED') gt 0 then begin
         if n_elements(prior_substruct['PRIOR']) gt 1 then begin
           for i=0, n_elements(prior_substruct['PRIOR'])-1 do begin 
             tabulated_priors = [[tabulated_priors], [param_name+'_'+strtrim(string(i+1,f='(I0)')), $
                                                      prior_substruct['PRIOR_ARG', i]]]
             init_range = [[init_range], [reform(config_hash[param_name, 'INITIALIZATION_RANGE', i, *])]]
           endfor
         endif else begin
           tabulated_priors = [[tabulated_priors], [param_name, prior_substruct['PRIOR_ARG']]]
           init_range = [[init_range], [config_hash[param_name, 'INITIALIZATION_RANGE', *]]]
         endelse
       endif
     endif
   endif
 endforeach


; All parameters with tabulated priors must be in the same file for a given galaxy.
 if n_elements(tabulated_priors) gt 0 then begin
   prior_dir = tabulated_priors[1, *]
   if n_elements(uniq(prior_dir, sort(prior_dir))) ne 1 then $
     message, 'The directories specified for all tabulated priors must be the same. '+$
              'Please check that the specified directories for all parameters with '+$
              'tabulated priors are identical.'

   ; Check for proper format of data in each prior file.
   for i=0, n_elements(sed_id)-1 do begin
     if file_test(tabulated_priors[1,0]+'tabulated_prior_'+sed_id[i]+'.fits') then begin
       tabulated_file = mrdfits(tabulated_priors[1,0]+'tabulated_prior_'+sed_id[i]+'.fits', 1, /silent)
       prior_tags = tag_names(tabulated_file)
       expected_prior_tags = [reform(tabulated_priors[0,*])+'_values', $
                              reform(tabulated_priors[0,*])+'_pdf']
       prior_tags = prior_tags[sort(prior_tags)]
       expected_prior_tags = expected_prior_tags[sort(expected_prior_tags)]
       if total(strupcase(prior_tags) eq strupcase(expected_prior_tags)) ne n_elements(expected_prior_tags) then $
         message, 'The tabulated prior .fits file for SED_ID '+sed_id[i]+' does '+$
                  'not have all the required parameters tabulated. Please '+$
                  'check that all parameters with tabulated priors are in the file '+$
                  'with the correct tag names.'

       ; Convert to hash as to index by name.
       tabulated_hash = orderedhash(tabulated_file, /extract)
       ; Check that each parameter has the same number of elements in values and pdf tags,
       ;   the values tag is in ascending order, and pdf tag is non-negative.
       foreach element, reform(tabulated_priors[0,*]), j do begin
         ele = strupcase(element)
         if size(tabulated_hash[ele+'_VALUES'], /type) lt 2 or size(tabulated_hash[ele+'_VALUES'], /type) gt 5 or $
            size(tabulated_hash[ele+'_PDF'], /type)    lt 2 or size(tabulated_hash[ele+'_PDF'], /type)    gt 5 then $
           message, 'The tabulated '+ele+' prior for SED_ID '+sed_id[i]+' is of incorrect data type.'

         if size(tabulated_hash[ele+'_VALUES'], /n_dim) ne 1 or size(tabulated_hash[ele+'_PDF'], /n_dim) ne 1 then $
           message, 'The tabulated '+ele+' prior for SED_ID '+sed_id[i]+' must be a 1D array for '+$
                    'both the PDF and corresponding values.'

         if n_elements(tabulated_hash[ele+'_VALUES']) ne n_elements(tabulated_hash[ele+'_PDF']) then $
           message, 'The tabulated '+ele+' prior for SED_ID '+sed_id[i]+' does not have '+$
                    'the same number of elements for the PDF and corresponding values.'

         if total(sort(tabulated_hash[ele+'_VALUES']) ne lindgen(n_elements(tabulated_hash[ele+'_VALUES']))) ne 0 then $
           message, 'The tabulated '+ele+' prior for SED_ID '+sed_id[i]+' must have its '+$
                    'grid of values in ascending order.'

         if total(tabulated_hash[ele+'_PDF'] lt 0) gt 0 then $
           message, 'The tabulated '+ele+' prior for SED_ID '+sed_id[i]+' must have non-negative '+$
                    'values for the PDF.'

         if init_range[0, j] lt min(tabulated_hash[ele+'_VALUES']) then $
           message, "Error found in the configuration structure. INITIALIZATION_RANGE tag of "+$
                    ele+" substructure must have a min bound greater than or equal to "+$
                    "the minimum gridded bound of the tabulated prior for SED_ID "+sed_id[i]+"."
         if init_range[1, j] gt max(tabulated_hash[ele+'_VALUES']) then $
           message, "Error found in the configuration structure. INITIALIZATION_RANGE tag of "+$
                    ele+" substructure must have a max bound less than or equal to "+$
                    "the maximum gridded bound of the tabulated prior for SED_ID "+sed_id[i]+"."

         in_initialize_range = where(tabulated_hash[ele+'_VALUES'] ge init_range[0, j] and $
                                     tabulated_hash[ele+'_VALUES'] le init_range[1, j], /null)
         if n_elements(in_initialize_range) le 1 then $
           message, "Error found in the configuration structure. INITIALIZATION_RANGE tag of "+$
                    ele+" substructure must have a range that includes multiple gridded values "+$
                    "of the tabulated prior for SED_ID "+sed_id[i]+"."
       endforeach

     endif else begin
       message, 'The tabulated prior .fits file was not found within '+$
                'the specified directory for SED_ID '+sed_id[i]+'. Please '+$
                'check that the tabulated prior .fits file exists for the SED.'
     endelse 
   endfor
 endif

end