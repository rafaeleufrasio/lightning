function prior_interactive, param_name, param_descript, param_range, default_prior, $
                            default_prior_args, default_initialization, prior_options, $
                            edit=edit, prior=prior, initialization_range=initialization_range
;+
; Name
; ----
;   PRIOR_INTERACTIVE
;
; Purpose
; -------
;   Generates the prior type, corresponding shape arguments, and
;   initialization range for a prior structure in the Lightning 
;   configuration using interactive prompts. It additionally checks
;   the prompt inputs for errors and ensures proper format.
;
; Calling Sequence
; ----------------
;   ::
;
;       prior_arg = prior_interactive(param_name, param_descript, param_range, default_prior, $
;                                     default_prior_args, default_initialization, prior_options [, $
;                                     /edit, prior=prior, initialization_range=initialization_range])
;
; Inputs
; ------
;   ``param_name`` : string scalar
;       The name of the parameter associated with the prior.
;   ``param_descript`` : string scalar
;       The description of the parameter.
;   ``param_range`` : int, float, or double array(2)
;       The allowed range of the parameter, given in terms of ``[min, max]``.
;   ``default_prior`` : string scalar
;       The default prior type choice for the parameter.
;   ``default_prior_args`` : int, float, or double array(5)
;       The default shape arguments associated with each prior type, given
;       in order of ``[fixed_value, min_bound, max_bound, normal_peak, normal_stddev]``.
;   ``default_initialization`` : int, float, or double array(2)
;       The default initialization range of the parameter indicating the minimum
;       and maximum bounds for the random initialization of the fitting algorithm.
;   ``prior_options`` : string array(Nprior_options)
;       The allowed names of the types of priors that the parameter can be
;       set (e.g., ``'uniform'`` or ``'normal'``).
;
; Optional Input
; --------------
;   ``edit`` : flag
;       If set, then the configuration is being edited and previous values should
;       be given rather than defaults.
;
; Output
; ------
;   ``prior_arg`` : int, float, double, or string array(Narg)
;       The distribution shape argument array whose size depends
;       on the chosen distribution type.
;
; Optional Outputs
; ----------------
;   ``prior`` : string scalar
;       The chosen prior distribution type.
;   ``initialization_range`` : double array(2)
;       The chosen initialization range of the parameter.
;
; Modification History
; --------------------
;   - 2022/04/26: Created (Keith Doore)
;   - 2022/05/05: Allowed for previous configuration to be printed with ``edit`` keyword (Keith Doore)
;   - 2022/05/19: Allowed for linear prior inputs to be converted to log-space (Keith Doore)
;   - 2022/05/31: Removed log priors (Keith Doore)
;   - 2022/06/27: Updated documentation (Keith Doore)
;   - 2022/06/27: Updated variable names (Keith Doore)
;   - 2022/08/01: Added initialization range (Keith Doore)
;-
 Compile_opt idl2
 On_error,2

 prior_arg = !null
 initialization_range = dblarr(2)
 if keyword_set(edit) then default_str = '(Previous configuration: ' else default_str = '(Default: '


; ========= Prior type ========
 error = 1
 prior_message = 'Please specify the type of prior to use on '+param_name+', '+$
                 param_descript+'. Current options: '+strjoin(prior_options, ', ')+$
                 '. '+default_str+default_prior+')'
 print_to_width, '======================='
 while error do begin
   print_to_width, ''
   print_to_width, prior_message
 
   temp = ''
   read, temp, prompt=strupcase(param_name)+': ', form='(A0)'
   if temp eq '' then temp = strupcase(default_prior)
   temp = strtrim(temp, 2) 
 
   if total(strupcase(temp) eq strupcase(prior_options)) eq 0 then begin
     prior_message = 'Please specify one of the current options: '+strjoin(prior_options, ', ')+'.'
   endif else error = 0 
 endwhile
 prior = temp

 param_range_str = strtrim(string(param_range[0]), 2)+' to '+strtrim(string(param_range[1]), 2)
 param_defaults_str = strtrim(string(default_prior_args), 2)
 initial_defaults_str = strtrim(string(default_initialization), 2)

 
; ========= Prior shape arguments ========
  ; ========= Fixed Prior ========
 if strupcase(prior) eq 'FIXED' then begin
   if strupcase(default_prior)  eq 'FIXED' and keyword_set(edit) then $
     default_str = '(Previous configuration: ' else default_str = '(Default: '
   error = 1
   prior_arg_message = 'Please specify the value at which to fix the '+param_name+$
                       ' prior within the allowed range of '+param_range_str+'. '+$
                       default_str+param_defaults_str[0]+')'
   print_to_width, '======================='
   while error do begin
     print_to_width, ''
     print_to_width, prior_arg_message
   
     temp = ''
     read, temp, prompt='FIXED_'+strupcase(param_name)+': ', form='(A0)'
     if temp eq '' then temp = param_defaults_str[0]
     temp = double(temp)

     if temp lt param_range[0] or temp gt param_range[1] then begin
       prior_arg_message = 'Please specify a value within the allowed range of '+param_range_str+'.'
     endif else error = 0 
   endwhile
   prior_arg = [prior_arg, temp]
 
  ; ========= Tabulated Prior ========
 endif else if strupcase(prior) eq 'TABULATED' then begin
   error = 1
   prior_arg_message = 'Please specify the path to the tabulated prior file '+$
                       'for '+param_name+' (e.g., /<prior_path>/).'
   print_to_width, '======================='
   while error do begin
     print_to_width, ''
     print_to_width, prior_arg_message
   
     temp = ''
     read, temp, prompt='TABULATED_'+strupcase(param_name)+': ', form='(A0)'
   
     file_path = file_test(temp, /directory)
     if ~file_path then begin
       path_message = 'Given path does not lead to any existing directory. Please '+$
                      'specify the correct path.'
     endif else error = 0 
   endwhile
   prior_arg = [prior_arg, temp]
 
  ; ========= Min/max bounds of uniform and normal priors ========
 endif else if total(strupcase(prior) eq ['UNIFORM', 'NORMAL']) eq 1 then begin
   if total(strupcase(default_prior) eq ['UNIFORM', 'NORMAL']) eq 1 and keyword_set(edit) then $
     default_str = '(Previous configuration: ' else default_str = '(Default: '
   error = 1
   prior_arg_message = 'Please specify the prior minimum bound for '+strupcase(param_name)+$
                       ' within the allowed range of '+param_range_str+'. '+$
                       default_str+param_defaults_str[1]+')'
   print_to_width, '======================='
   while error do begin
     print_to_width, ''
     print_to_width, prior_arg_message
   
     temp = ''
     read, temp, prompt='MIN_BOUND_'+strupcase(param_name)+': ', form='(A0)'
     if temp eq '' then temp = param_defaults_str[1]
     temp = double(temp)
   
     if temp lt param_range[0] or temp gt param_range[1] then begin
       prior_arg_message = 'Please specify a value within the allowed range of '+param_range_str+'.'
     endif else error = 0 
   endwhile
   prior_arg = [prior_arg, temp]
 
   error = 1
   prior_arg_message = 'Please specify the prior maximum bound for '+strupcase(param_name)+$
                       ' within the allowed range of '+param_range_str+'. '+$
                       default_str+param_defaults_str[2]+')'
   print_to_width, '======================='
   while error do begin
     print_to_width, ''
     print_to_width, prior_arg_message
   
     temp = ''
     read, temp, prompt='MAX_BOUND_'+strupcase(param_name)+': ', form='(A0)'
     if temp eq '' then temp = param_defaults_str[2]
     temp = double(temp)

     if temp le prior_arg[0] or temp gt param_range[1] then begin
       prior_arg_message = 'Please specify a maximum bound larger than the minimum bound and'+$
                           ' within the allowed range of '+param_range_str+'.'
     endif else error = 0 
   endwhile
   prior_arg = [prior_arg, temp]
 
  ; ========= Expected value and stddev of normal priors ========
   if total(strupcase(prior) eq ['NORMAL', 'LOG-NORMAL']) eq 1 then begin
     if total(strupcase(default_prior) eq ['NORMAL', 'LOG-NORMAL']) eq 1 and keyword_set(edit) then $
       default_str = '(Previous configuration: ' else default_str = '(Default: '
     error = 1
     prior_arg_message = 'Please specify the peak value of the normal prior for '+$
                         strupcase(param_name)+'. '+default_str+param_defaults_str[3]+')'
     print_to_width, '======================='
     print_to_width, ''
     print_to_width, prior_arg_message
     
     temp = ''
     read, temp, prompt='PEAK_'+strupcase(param_name)+': ', form='(A0)'
     if temp eq '' then temp = param_defaults_str[3]
     temp = double(temp)
     
     prior_arg = [prior_arg, temp]
 
     prior_arg_message = 'Please specify the standard deviation of the normal prior for '+$
                          strupcase(param_name)+'. '+default_str+param_defaults_str[4]+')'
     print_to_width, '======================='
     while error do begin
       print_to_width, ''
       print_to_width, prior_arg_message
     
       temp = ''
       read, temp, prompt='STDDEV_'+strupcase(param_name)+': ', form='(A0)'
       if temp eq '' then temp = param_defaults_str[4]
       temp = double(temp)
     
       case 1 of
         temp le 0: $
           prior_arg_message = 'Please specify a positive value.'

         temp gt 10.d*(prior_arg[1]-prior_arg[0]): $
           prior_arg_message = 'Standard deviation is comparable in size to the effective distribution range, '+$
                               'resulting in an effective uniform distribution. Please reduce the '+$
                               'standard deviation value.'

         else: error = 0 
       endcase
     endwhile
     prior_arg = [prior_arg, temp]
   endif
 endif


; ========= Initialization bounds ========
 if total(strupcase(prior) eq ['UNIFORM', 'NORMAL', 'TABULATED']) eq 1 then begin
   if total(strupcase(default_prior) eq ['UNIFORM', 'NORMAL', 'TABULATED']) eq 1 and keyword_set(edit) then $
     default_str = '(Previous configuration: ' else default_str = '(Default: '
   error = 1
   prior_arg_message = 'Please specify the initialization minimum bound for '+strupcase(param_name)+$
                       ' within the allowed range of '+param_range_str+'. '+$
                       default_str+initial_defaults_str[0]+')'
   print_to_width, '======================='
   while error do begin
     print_to_width, ''
     print_to_width, prior_arg_message
   
     temp = ''
     read, temp, prompt='INITIAL_MIN_BOUND_'+strupcase(param_name)+': ', form='(A0)'
     if temp eq '' then temp = initial_defaults_str[0]
     temp = double(temp)
   
     if temp lt param_range[0] or temp gt param_range[1] then begin
       prior_arg_message = 'Please specify a value within the allowed range of '+param_range_str+'.'
     endif else error = 0 
   endwhile
   initialization_range[0] = temp
 
   error = 1
   prior_arg_message = 'Please specify the initialization  maximum bound for '+strupcase(param_name)+$
                       ' within the allowed range of '+param_range_str+'. '+$
                       default_str+initial_defaults_str[1]+')'
   print_to_width, '======================='
   while error do begin
     print_to_width, ''
     print_to_width, prior_arg_message
   
     temp = ''
     read, temp, prompt='INITIAL_MAX_BOUND_'+strupcase(param_name)+': ', form='(A0)'
     if temp eq '' then temp = initial_defaults_str[1]
     temp = double(temp)

     if temp lt initialization_range[0] or temp gt param_range[1] then begin
       prior_arg_message = 'Please specify a initialization maximum bound greater than or equal to the'+$
                           ' minimum bound and within the allowed range of '+param_range_str+'.'
     endif else error = 0 
   endwhile
   initialization_range[1] = temp
 endif


 return, prior_arg

end