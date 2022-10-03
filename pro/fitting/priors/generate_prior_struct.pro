function generate_prior_struct, config, sed_id, config_nopriors=config_nopriors
;+
; Name
; ----
;   GENERATE_PRIOR_STRUCT
;
; Purpose
; -------
;   Generates the priors structure to be used when computing the
;   prior probability of each parameter.
;
; Calling Sequence
; ----------------
;   ::
;
;       priors = generate_prior_struct(config, sed_id [, config_nopriors=config_nopriors])
;
; Inputs
; ------
;   ``config`` : structure
;       A Lightning configuration structure. (See
;       ``lightning_configure_defaults.pro`` for details and contents.)
;   ``sed_id`` : string scalar
;       A unique SED identifier
;
; Output
; ------
;   ``priors`` : structure
;       This structure includes the hyper-parameters needed to compute the prior
;       probability for each parameter in the chosen models. It also includes
;       the indices indicating what type of prior each parameter has.
;       The full description of the structure is as follows:
; 
;       ================     =================     ================================================================================
;       TAG                  TYPE                  DESCRIPTION
;       ================     =================     ================================================================================
;       PARAMETER_NAME       string(Nparam)        Name of each parameter
;       IDCS.FIXED           int(...)              Indices of parameters that are fixed (size varies)
;       IDCS.VARIABLE        int(...)              Indices of parameters that are not fixed (size varies)
;       IDCS.UNIFORM         int(...)              Indices of parameters that have uniform priors (size varies)
;       IDCS.NORMAL          int(...)              Indices of parameters that have normal priors (size varies)
;       IDCS.ANALYTICAL      int(...)              Indices of parameters that have uniform or normal priors (size varies)
;       IDCS.TABULATED       int(...)              Indices of parameters that have tabulated priors (size varies)
;       INITIALIZE_RANGE     double(Nparam, 2)     The lower and upper bounds for the initialization range of the fitting algorithm
;       FIXED                double(Nparam)        Values of fixed parameters (``NaN`` if not fixed)
;       MIN_BOUND            double(Nparam)        Values of the minimum bounds for the priors
;       MAX_BOUND            double(Nparam)        Values of the maximum bounds for the priors
;       ANALYTICAL.WIDTH     double(Nparam)        Values of the width for the normal priors (:math:`-0.5/\sigma^2`)
;       ANALYTICAL.MU        double(Nparam)        Values of the peak value for the normal priors
;       ANALYTICAL.CONST     double(Nparam)        Values of the normalization constant for the uniform and normal priors
;       TABULATED            structure             A structure containing the tabulated values for the tabulated priors
;       ================     =================     ================================================================================
;
; Optional Output
; ---------------
;   ``config_nopriors`` : structure
;       The input Lightning configuration structure edited to remove the prior 
;       substructures for each parameter.
;
; Modification History
; --------------------
;   - 2022/06/06: Created (Keith Doore)
;   - 2022/06/20: Fixed issue with tabulated prior indexing of hash keys (Keith Doore)
;   - 2022/07/05: Allowed for other parameters with multiple priors besides ``PSI`` (Keith Doore)
;   - 2022/07/05: Removed parameter substructures from ``config`` after extracting values (Keith Doore)
;   - 2022/08/01: Added ``/silent`` to ``mrdfits`` (Keith Doore)
;   - 2022/08/02: Added ``initialization_range`` as ``initialize_range`` to prior structure (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(config) eq 0 then message, 'Variable is undefined: CONFIG.'
 if size(config, /type) ne 8 then message, 'CONFIG is not of type structure.'

 if n_elements(sed_id) eq 0 then message, 'Variable is undefined: SED_ID.'
 if size(sed_id, /type) ne 7 then message, 'SED_ID is not of type string.'
 if size(sed_id, /n_dim) ne 0 then message, 'SED_ID must be a scalar.'



; Make list of all parameter names and their prior types.
;   If a parameter has multiple priors (e.g., PSI), separate them out for easier 
;   extraction, will not change config
 config_hash = orderedhash(config, /extract)

 foreach prior_substruct, config_hash, param_name do begin
   if size(prior_substruct, /type) eq 11 then begin
     if prior_substruct.haskey('PRIOR') then begin
       if n_elements(prior_substruct['PRIOR']) gt 1 then begin
         for i=0, n_elements(prior_substruct['PRIOR'])-1 do begin
           config_hash[param_name+'_'+strtrim(string(i+1,f='(I0)')), 'PRIOR'] = config_hash[param_name, 'PRIOR', i]
           config_hash[param_name+'_'+strtrim(string(i+1,f='(I0)')), 'PRIOR_ARG'] = $
                 reform(config_hash[param_name, 'PRIOR_ARG', i, *])
           config_hash[param_name+'_'+strtrim(string(i+1,f='(I0)')), 'INITIALIZATION_RANGE'] = $
                 reform(config_hash[param_name, 'INITIALIZATION_RANGE', i, *])
         endfor
         config_hash.remove, param_name
       endif
     endif
   endif
 endforeach


; Extract parameter name and prior type
 prior_name_type = !null
 foreach prior_substruct, config_hash, param_name do begin
   if size(prior_substruct, /type) eq 11 then begin
     if prior_substruct.haskey('PRIOR') then begin
       prior_name_type = [[prior_name_type], [param_name, prior_substruct['PRIOR']]]
     endif
   endif
 endforeach
 Nparam = n_elements(prior_name_type[0, *])
 prior_name_type = strupcase(prior_name_type)



; Find indices of for each prior type
 variable_param = where(prior_name_type[1, *] ne 'FIXED', comp=fixed_param, /null)
 tabulate_param = where(prior_name_type[1, *] eq 'TABULATED', /null)
 analytic_param = where(prior_name_type[1, *] eq 'UNIFORM' or prior_name_type[1, *] eq 'NORMAL', /null)
 normal_param   = where(prior_name_type[1, *] eq 'NORMAL', /null)
 uniform_param  = where(prior_name_type[1, *] eq 'UNIFORM', /null)



; Get initialization range
 initialization_range = dblarr(Nparam, 2)
 for i=0, Nparam-1 do initialization_range[i, *] = config_hash[prior_name_type[0, i], 'INITIALIZATION_RANGE']
 ; Set range for fixed parameters to fixed values
 if n_elements(fixed_param) gt 0 then begin
   for i=0, n_elements(fixed_param)-1 do $
     initialization_range[fixed_param[i], *] = config_hash[prior_name_type[0, fixed_param[i]], 'PRIOR_ARG']
 endif


; Get fixed values
 fixed_values = replicate(!values.D_NaN, Nparam)
 if n_elements(fixed_param) gt 0 then begin
   for i=0, n_elements(fixed_param)-1 do $
     fixed_values[fixed_param[i]] = config_hash[prior_name_type[0, fixed_param[i]], 'PRIOR_ARG']
 endif



; Generate analytical prior coefficients for easy computation in `lightning_priors.pro`
;   Parameters in order are: a -> minimum bound, b -> maximum bound, mu -> expected value (mean), sig -> standard deviation
 analytic_values = dblarr(Nparam, 4)
 if n_elements(analytic_param) gt 0 then begin
   for i=0, n_elements(analytic_param)-1 do begin
     if prior_name_type[1, analytic_param[i]] eq 'UNIFORM' then temp=1 else temp=3
     analytic_values[analytic_param[i], 0:temp] = config_hash[prior_name_type[0, analytic_param[i]], 'PRIOR_ARG', 0:temp]
   endfor
 endif
 ; lnprob = width *(x-mu)^2 + const
 ; Uniform  -> width = 0.d0
 ;             mu = 0.d
 ;             const = alog(1.d0/(b-a))
 ; Gaussian -> width = -0.5d/sig^2.d
 ;             mu = mu
 ;             const = alog(2.d/(sig*sqrt(2.d*!dpi)) / (erf((b-mu)/sig/sqrt(2.d))-erf((a-mu)/sig/sqrt(2.d))))
 width_arr = dblarr(Nparam)
 mu_arr = dblarr(Nparam)
 const_arr = dblarr(Nparam)
 if n_elements(uniform_param) gt 0 then begin
   const_arr[uniform_param] = alog(1.d0/(analytic_values[uniform_param, 1]-analytic_values[uniform_param, 0]))
 endif
 if n_elements(normal_param) gt 0 then begin
   width_arr[normal_param] = -0.5d/analytic_values[normal_param, 3]^2.d
   mu_arr[normal_param] = analytic_values[normal_param, 2]
   const_arr[normal_param] = alog(2.d/(analytic_values[normal_param, 3]*sqrt(2.d*!dpi)) / $
     (erf((analytic_values[normal_param, 1]-analytic_values[normal_param, 2])/analytic_values[normal_param, 3]/sqrt(2.d)) - $
      erf((analytic_values[normal_param, 0]-analytic_values[normal_param, 2])/analytic_values[normal_param, 3]/sqrt(2.d))))
 endif
 ; There should be no NaNs or Infinities due to checking config, but double check to confirm
 if total(~finite(width_arr)) ne 0 or total(~finite(mu_arr)) ne 0 or total(~finite(const_arr)) ne 0 then $
   message, 'Non-finite values have occurred in the prior for SED_ID '+sed_id+'. Please update '+$
            'the prior shape arguments to prevent non-finite values.'



; Generate minimum and maximum bounds for each parameter.
;   For tabulated priors, this is the first and last gridded values.
 min_bounds_arr = dblarr(Nparam)
 max_bounds_arr = dblarr(Nparam)
 if n_elements(fixed_param) gt 0 then begin
   min_bounds_arr[fixed_param] = fixed_values[fixed_param]
   max_bounds_arr[fixed_param] = fixed_values[fixed_param]
 endif
 if n_elements(analytic_param) gt 0 then begin
   min_bounds_arr[analytic_param] = analytic_values[analytic_param, 0]
   max_bounds_arr[analytic_param] = analytic_values[analytic_param, 1]
 endif


; Generate tabulated prior structure for easy computation in `lightning_priors.pro`.
;   The structure consists of substructures for each parameter with tags of pdf and
;   values. Parameters without tabulated priors will just have a NaN within each tag.
;   Note that values have already been checked for errors in `tabulated_prior_check.pro`.
 tabulated_struct = orderedhash()
 if n_elements(tabulate_param) gt 0 then begin
   tabulated_file = mrdfits(config_hash[prior_name_type[0, tabulate_param[0]], 'PRIOR_ARG']+$
                            'tabulated_prior_'+sed_id+'.fits', 1, /silent)
   tabulated_hash = orderedhash(tabulated_file, /extract)
 endif
 for i=0, Nparam-1 do begin
   if total(i eq tabulate_param) eq 1 then begin
     tabulated_struct[prior_name_type[0, i], 'VALUES'] = tabulated_hash[strupcase(prior_name_type[0, i])+'_VALUES']
     tabulated_struct[prior_name_type[0, i], 'PDF'] = tabulated_hash[strupcase(prior_name_type[0, i])+'_PDF']
     min_bounds_arr[i] = tabulated_hash[strupcase(prior_name_type[0, i])+'_VALUES', 0]
     max_bounds_arr[i] = tabulated_hash[strupcase(prior_name_type[0, i])+'_VALUES', -1]
   endif else begin
     tabulated_struct[prior_name_type[0, i], 'VALUES'] = !values.D_NaN
     tabulated_struct[prior_name_type[0, i], 'PDF'] = !values.D_NaN
   endelse
 endfor
 tabulated_struct = tabulated_struct.tostruct(/recursive)


; Edit the Lightning configuration structure to remove the prior substructures for each parameter.
 config_hash.Remove, reform(prior_name_type[0, *])
 config_nopriors = config_hash.ToStruct(/recursive)



; Create the prior structure
;   Indices are not guaranteed to have a value (may be null), if they are empty
;   replace with NaN
 if n_elements(fixed_param   ) eq 0 then fixed_param    = !values.D_NaN
 if n_elements(variable_param) eq 0 then variable_param = !values.D_NaN
 if n_elements(uniform_param ) eq 0 then uniform_param  = !values.D_NaN
 if n_elements(normal_param  ) eq 0 then normal_param   = !values.D_NaN
 if n_elements(analytic_param) eq 0 then analytic_param = !values.D_NaN
 if n_elements(tabulate_param) eq 0 then tabulate_param = !values.D_NaN

 priors = {parameter_name: reform(prior_name_type[0, *]) ,$ 
           idcs: {fixed: fixed_param                     ,$
                  variable: variable_param               ,$
                  uniform: uniform_param                 ,$
                  normal: normal_param                   ,$
                  analytical: analytic_param             ,$
                  tabulated: tabulate_param               $
                  }                                      ,$
           initialize_range: initialization_range        ,$
           fixed: fixed_values                           ,$
           min_bound: min_bounds_arr                     ,$
           max_bound: max_bounds_arr                     ,$
           analytical: {width: width_arr                 ,$
                        mu: mu_arr                       ,$
                        const: const_arr                  $
                        }                                ,$
           tabulated: tabulated_struct                    $
           }


 return, priors

end