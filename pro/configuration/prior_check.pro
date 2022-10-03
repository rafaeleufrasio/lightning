pro prior_check, prior_struct, param_name, npriors, param_range, prior_options, $
                 prior_options_narg
;+
; Name
; ----
;   PRIOR_CHECK
;
; Purpose
; -------
;   Checks the Lightning configuration prior substructures to make
;   sure all inputs are correctly formatted and can be understood by
;   Lightning. If they are not, an error message is printed.
;
; Calling Sequence
; ----------------
;   ::
;
;       prior_check, prior_struct, param_name, npriors, $
;                    param_range, prior_options, prior_options_narg
;
; Inputs
; ------
;   ``prior_struct`` : structure
;       The prior structure for a parameter as generated from the Lightning 
;       configuration. The full description of the structure is as follows:
;
;       ====================     ======================================     ==============================
;       TAG                      TYPE                                       DESCRIPTION
;       ====================     ======================================     ==============================
;       PRIOR                    string(Npriors)                            Chosen prior distribution type
;       PRIOR_ARG                int/float/double/string(Npriors, Narg)     Distribution shape arguments
;       INITIALIZATION_RANGE     int/float/double(Npriors, 2)               Parameter initialization range
;       ====================     ======================================     ==============================
;
;   ``param_name`` : string scalar
;       The name of the parameter associated with the prior structure.
;   ``npriors`` : int, float, or double scalar
;       The number of priors that should be specified for the parameter.
;   ``param_range`` : int, float, or double array(2)
;       The allowed range of the parameter contained within the prior 
;       substructure, given in terms of ``[min, max]``.
;   ``prior_options`` : string array(Nprior_options)
;       The allowed names of the types of priors that the parameter can be
;       set (e.g., ``'uniform'`` or ``'normal'``).
;   ``prior_options_narg`` : int, float, or double array(Nprior_options)
;       The number of arguments associated with each with each prior option.
;
; Output
; ------
;   An error message is output if any errors are found within the prior substructures.
;   Error messages are relatively specific and should help users resolve any issues.
;
; Modification History
; --------------------
;   - 2022/04/26: Created (Keith Doore)
;   - 2022/05/19: Added check to ensure tabulated directory exists (Keith Doore)
;   - 2022/06/02: Added additional checks to confirm parameters do not result in priors with
;     machine over/underflow (Keith Doore)
;   - 2022/06/27: Renamed inputs for consistent naming scheme (Keith Doore)
;   - 2022/08/01: Added check for initialization range (Keith Doore)
;   - 2022/08/01: Fixed bug if input multiple priors and not all had same ``Narg`` (Keith Doore)
;-
 On_error, 2
 Compile_opt idl2

; No error handling of inputs is included since there should be no user interaction 
;  with this function.
 tags = tag_names(prior_struct)
 param_name = strupcase(param_name)
 base_err = 'Error found in the configuration structure. '


; Check PRIOR tag
 if npriors eq 1 then begin
   prior_size = 'scalar.'
   prior_dim = 0
   prior_element = ''
 endif else if npriors gt 1 then begin
   prior_size = '1-D array.'
   prior_dim = 1
   prior_element = 'Each element in the '
 endif

 if n_elements(where(strupcase(tags) eq 'PRIOR', /null)) eq 0 then $
   message, base_err+'PRIOR tag is missing for '+param_name+' substructure.'
 
 if size(prior_struct.PRIOR, /type) ne 7 then $
   message, base_err+'PRIOR tag of '+param_name+' substructure must be of type string.'

 if n_elements(prior_struct.PRIOR) ne npriors then $
   message, base_err+'PRIOR tag of '+param_name+' substructure must have '+$
            string(npriors, f='(I0)')+' element(s).'

 if size(prior_struct.PRIOR, /n_dim) ne prior_dim then $
   message, base_err+'PRIOR tag of '+param_name+' substructure must be a '+prior_size

 foreach i, prior_struct.PRIOR do if total(strupcase(i) eq strupcase(prior_options)) eq 0 then $
   message, base_err+prior_element+'PRIOR tag of '+param_name+' substructure must '+$
            "be set to one of the following: '"+strjoin(prior_options, "', '")+"'."

 if npriors gt 1 then $
   if total(total(strupcase(prior_struct.PRIOR) eq 'TABULATED') eq [0, npriors]) eq 0 then $
     message, base_err+"All or no elements of PRIOR tag of "+param_name+" substructure must be 'tabulated'."



; Check PRIOR_ARG tag
 if n_elements(where(strupcase(tags) eq 'PRIOR_ARG', /null)) eq 0 then $
   message, base_err+'PRIOR_ARG tag is missing for '+param_name+' substructure.'

 ; Confirm that any potential extra dimensions (brackets) are removed, or
 ;   any missing brackets are added.
 prior_struct.PRIOR_ARG = reform(prior_struct.PRIOR_ARG)

 Narg = intarr(npriors)
 for i=0, npriors-1 do begin
   ind = where(strupcase(prior_struct.PRIOR[i]) eq strupcase(prior_options), /null)
   Narg[i] = prior_options_narg[ind]
 endfor


 if strupcase(prior_struct.PRIOR[0]) eq 'TABULATED' then begin

   if size(prior_struct.PRIOR_ARG, /type) ne 7 then $
     message, base_err+"PRIOR_ARG tag of "+param_name+" substructure must be of type string"+$
              " if 'tabulated' prior is used."

   if size(prior_struct.PRIOR_ARG, /n_dim) ne 1 then $
     message, base_err+"PRIOR_ARG tag of "+param_name+" substructure must be a 1-D array "+$
              " if 'tabulated' prior is used."

   file_path = file_test(prior_struct.PRIOR_ARG, /directory)
   if total(file_path) ne n_elements(file_path) then $
     message, base_err+"PRIOR_ARG tag of "+param_name+" substructure has 'tabulated' prior that"+$
              "does not lead to any existing directory."

 endif else begin

   if size(prior_struct.PRIOR_ARG, /type) lt 2 or size(prior_struct.PRIOR_ARG, /type) gt 5 then $ 
     message, base_err+'PRIOR_ARG tag of '+param_name+' substructure must be of type int, float, or double.'

   if npriors eq 1 then begin

     if size(prior_struct.PRIOR_ARG, /n_dim) ne 1 then $
       message, base_err+"PRIOR_ARG tag of "+param_name+" substructure must be a 1-D array."

   endif else begin

     if total(strupcase(prior_struct.PRIOR) eq 'FIXED') eq npriors then begin
       if size(prior_struct.PRIOR_ARG, /n_dim) ne 1 then $
         message, base_err+"PRIOR_ARG tag of "+param_name+" substructure must be a 1-D array if "+$
                  "all priors are 'fixed'."
     endif else begin
       if size(prior_struct.PRIOR_ARG, /n_dim) ne 2 then $
         message, base_err+"PRIOR_ARG tag of "+param_name+" substructure must be a 2-D array with "+$
                  "current prior choices."
     endelse
   endelse
 endelse

 if npriors eq 1 then begin
   if size(prior_struct.PRIOR_ARG, /dim) ne Narg then $
      message, base_err+'PRIOR_ARG tag of '+param_name+' substructure must have Narg elements.'
 endif else begin
   if (size(prior_struct.PRIOR_ARG, /dim))[0] ne npriors then $
     message, base_err+'PRIOR_ARG tag of '+param_name+' substructure must have a first dimension'+$
              ' with '+string(npriors, f='(I0)')+' elements.'

   if max(Narg) gt 1 then $
     if (size(prior_struct.PRIOR_ARG, /dim))[1] ne max(Narg) then $
       message, base_err+'PRIOR_ARG tag of '+param_name+' substructure must have a second'+$
                ' dimension with the largest Narg value of elements.'
 endelse


 if total(strupcase(prior_struct.PRIOR) eq 'TABULATED') ne npriors then begin
   for i=0, Npriors-1 do begin

     if npriors eq 1 then prior_arg = prior_struct.PRIOR_ARG else prior_arg = reform((prior_struct.PRIOR_ARG)[i, 0:Narg[i]-1])

     if Narg[i] eq 1 then begin
       if prior_arg lt param_range[0] or prior_arg gt param_range[1] then $
         message, base_err+"PRIOR_ARG tag of "+param_name+" substructure must have a 'fixed' prior "+$
                  "value within allowed parameter range."
     endif
     if Narg[i] ge 2 then begin
       if prior_arg[0] eq prior_arg[1] then $
         message, base_err+'PRIOR_ARG tag of '+param_name+' substructure must have unique min'+$
                  " and max bounds. If intending to fix the parameter, use 'fixed' prior intstead."+$
                  " Additionally, the bounds may be too small for machine precision to differentiate"+$
                  " between them."
       if prior_arg[0] gt prior_arg[1] then $
         message, base_err+'PRIOR_ARG tag of '+param_name+' substructure has incorrect ordering of '+$
                  'min and max bounds.'

       if prior_arg[0] lt param_range[0] or prior_arg[0] gt param_range[1] then $
         message, base_err+"PRIOR_ARG tag of "+param_name+" substructure must have a min bound on"+$
                  " the prior within the allowed parameter range."
       if prior_arg[1] lt param_range[0] or prior_arg[1] gt param_range[1] then $
         message, base_err+"PRIOR_ARG tag of "+param_name+" substructure must have a max bound on"+$
                  " the prior within the allowed parameter range."

       ; Check confirms that uniform distributions will be finite when computed
       if ~finite(alog(1.d/(prior_arg[1]-prior_arg[0]))) then $
         message, base_err+'PRIOR_ARG tag of '+param_name+' substructure has min'+$
                  " and max bounds that result in machine precision errors. Please specify"+$
                  " values that are finite and well within machine precision range."
     endif
     if Narg[i] eq 4 then begin
       if prior_arg[3] eq 0 then $
         message, base_err+"PRIOR_ARG tag of "+param_name+" substructure must have a positive "+$
                  "standard deviation of 'normal' prior. If intending to fix the "+$
                  "parameter, use 'fixed' prior intstead."
       if prior_arg[3] lt 0 then $
         message, base_err+"PRIOR_ARG tag of "+param_name+" substructure must have a positive "+$
                  "standard deviation of 'normal' prior."

       if prior_arg[3] gt 10.d*(prior_arg[1]-prior_arg[0]) then $
         message, base_err+"PRIOR_ARG tag of "+param_name+" substructure must have a "+$
                  "standard deviation comparable in size to the effective distribution range. "+$
                  "Current standard deviation results in effective uniform distribution. Please "+$
                  "use an 'uniform' prior instead or reduce standard deviation value."

       ; Check confirms that normal distributions will be finite when computed.
       ;   Very large or small sigmas can result in machine precision errors
       if ~finite(2.d/(prior_arg[3]*sqrt(2.d*!dpi))) or ~finite(-0.5d/prior_arg[3]^2.d) then $
         message, base_err+"PRIOR_ARG tag of "+param_name+" substructure has a "+$
                  "standard deviation that results in machine precision errors. Please specify"+$
                  " a smaller or larger value that will accommodate machine precision range."

       ; Issues with error function can occur when in the tails of the distribution or when
       ;   min and max bounds are close to each other.
       erf_cal = erf((prior_arg[1]-prior_arg[2])/prior_arg[3]/sqrt(2.d))-erf((prior_arg[0]-prior_arg[2])/prior_arg[3]/sqrt(2.d))
       if erf_cal eq 0 then $
         message, base_err+"PRIOR_ARG tag of "+param_name+" substructure results in a "+$
                  "prior with machine precision errors. This is likely due to being in a single tail "+$
                  "of the normal distribution or having min and max bounds that are very close together. "+$
                  "Please specify bounds that span a nice range of the distribution or adjust the "+$
                  "peak value of the distribution accordingly."

       ; Issues can occur with the values for the parameter when the squared difference from
       ;   mu is too large or small.
       if ~finite((prior_arg[0]-prior_arg[2])^2.d) or ~finite((prior_arg[1]-prior_arg[2])^2.d) then $
         message, base_err+"PRIOR_ARG tag of "+param_name+" substructure results in a "+$
                  "prior with machine precision errors. This is likely due to values of the parameter "+$
                  "exceeding machine precision when squared. Please adjust the prior bounds to be well "+$
                  "within machine precision when squared."

       ; Check for any final machine precision errors that may occur from combined normal dist parameters
       lnprob_min_bound = -0.5d/prior_arg[3]^2.d * (prior_arg[0]-prior_arg[2])^2.d + $
                          alog(erf_cal * 2.d/(prior_arg[3]*sqrt(2.d*!dpi)))
       lnprob_max_bound = -0.5d/prior_arg[3]^2.d * (prior_arg[1]-prior_arg[2])^2.d + $
                          alog(erf_cal * 2.d/(prior_arg[3]*sqrt(2.d*!dpi)))
       if ~finite(lnprob_min_bound) or ~finite(lnprob_max_bound) then $
         message, base_err+"PRIOR_ARG tag of "+param_name+" substructure results in a "+$
                  "prior with machine precision errors. This is a final catch-all error. "+$
                  "Resolving the error will require user inspection of prior arguments."
     endif
   endfor
 endif


 
; Check INITIALIZATION_RANGE tag
 if n_elements(where(strupcase(tags) eq 'INITIALIZATION_RANGE', /null)) eq 0 then $
   message, base_err+'INITIALIZATION_RANGE tag is missing for '+param_name+' substructure.'

 if size(prior_struct.INITIALIZATION_RANGE, /type) lt 2 or size(prior_struct.INITIALIZATION_RANGE, /type) gt 5 then $ 
   message, base_err+'INITIALIZATION_RANGE tag of '+param_name+' substructure must be of type int, float, or double.'

 if npriors eq 1 then begin
   if size(prior_struct.INITIALIZATION_RANGE, /n_dim) ne 1 then $
     message, base_err+"INITIALIZATION_RANGE tag of "+param_name+" substructure must be a 1-D array."
 endif else begin
   if size(prior_struct.INITIALIZATION_RANGE, /n_dim) ne 2 then $
     message, base_err+"INITIALIZATION_RANGE tag of "+param_name+" substructure must be a 2-D array."
 endelse


 for i=0, Npriors-1 do begin
   if npriors eq 1 then init_range = prior_struct.INITIALIZATION_RANGE else $
                        init_range = reform((prior_struct.INITIALIZATION_RANGE)[i, *])
   if npriors eq 1 then prior_arg = prior_struct.PRIOR_ARG else prior_arg = reform((prior_struct.PRIOR_ARG)[i, 0:Narg[i]-1])

   if Narg[i] ge 2 then begin
     if init_range[0] gt init_range[1] then $
       message, base_err+'INITIALIZATION_RANGE tag of '+param_name+' substructure has '+$
                  'incorrect ordering of min and max bounds.'

     if init_range[0] lt prior_arg[0] then $
       message, base_err+"INITIALIZATION_RANGE tag of "+param_name+" substructure must have a min bound"+$
                " greater than or equal to the minimum bound of the prior."
     if init_range[1] gt prior_arg[1] then $
       message, base_err+"INITIALIZATION_RANGE tag of "+param_name+" substructure must have a max bound"+$
                " less than or equal to the maximum bound of the prior."

     ; Check confirms that uniform distributions will be finite when computed
     if ~finite(alog(1.d/(init_range[1]-init_range[0]))) then $
       message, base_err+'INITIALIZATION_RANGE tag of '+param_name+' substructure has min'+$
                " and max bounds that result in machine precision errors. Please specify"+$
                " values that are finite and well within machine precision range."

   endif
   if Narg[i] eq 4 then begin
     ; Issues with error function can occur when in the tails of the distribution or when
     ;   min and max bounds are close to each other.
     erf_cal = erf((init_range[1]-prior_arg[2])/prior_arg[3]/sqrt(2.d))-erf((init_range[0]-prior_arg[2])/prior_arg[3]/sqrt(2.d))
     if erf_cal eq 0 then $
       message, base_err+"INITIALIZATION_RANGE tag of "+param_name+" substructure results in a randomly sampled "+$
                "prior distribution with machine precision errors. This is likely due to being in a single tail "+$
                "of the normal distribution or having min and max initialization bounds that are very close together. "+$
                "Please specify a initialization range that span a nice range of the distribution."

     ; Issues can occur with the values for the parameter when the squared difference from
     ;   mu is too large or small.
     if ~finite((init_range[0]-prior_arg[2])^2.d) or ~finite((init_range[1]-prior_arg[2])^2.d) then $
       message, base_err+"INITIALIZATION_RANGE tag of "+param_name+" substructure results in a randomly sampled "+$
                "prior with machine precision errors. This is likely due to values of the parameter "+$
                "exceeding machine precision when squared. Please adjust the initialization bounds bounds to be well "+$
                "within machine precision when squared."

     ; Check for any final machine precision errors that may occur from combined normal dist parameters
     lnprob_min_bound = -0.5d/prior_arg[3]^2.d * (init_range[0]-prior_arg[2])^2.d + $
                        alog(erf_cal * 2.d/(prior_arg[3]*sqrt(2.d*!dpi)))
     lnprob_max_bound = -0.5d/prior_arg[3]^2.d * (init_range[1]-prior_arg[2])^2.d + $
                        alog(erf_cal * 2.d/(prior_arg[3]*sqrt(2.d*!dpi)))
     if ~finite(lnprob_min_bound) or ~finite(lnprob_max_bound) then $
       message, base_err+"INITIALIZATION_RANGE tag of "+param_name+" substructure results in a randomly sampled "+$
                "prior with machine precision errors. This is a final catch-all error. "+$
                "Resolving the error will require user inspection of initialization range arguments."
   
   endif
   if strupcase(prior_struct.PRIOR[i]) eq 'TABULATED' then begin
    ;NOTE - Initialization_range is checked to be within min and max bounds of tabulated prior in tabulated_prior_check.pro
     if init_range[0] gt init_range[1] then $
       message, base_err+'INITIALIZATION_RANGE tag of '+param_name+' substructure has '+$
                  'incorrect ordering of min and max bounds.'

     if init_range[0] lt param_range[0] or init_range[0] gt param_range[1] then $
       message, base_err+"INITIALIZATION_RANGE tag of "+param_name+" substructure must have a min bound"+$
                " within the allowed parameter range."
     if init_range[1] lt param_range[0] or init_range[1] gt param_range[1] then $
       message, base_err+"INITIALIZATION_RANGE tag of "+param_name+" substructure must have a max bound"+$
                " within the allowed parameter range."
   endif
 endfor

end