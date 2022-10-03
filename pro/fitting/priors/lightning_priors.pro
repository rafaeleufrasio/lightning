function lightning_priors, parameters, priors, error_check=error_check, lnprior_param=lnprior_param, $
                           in_bounds_idcs=in_bounds_idcs, out_bounds_idcs=out_bounds_idcs
;+
; Name
; ----
;   LIGHTNING_PRIORS
;
; Purpose
; -------
;   Calculates the log probability of each prior using the input prior
;   shape arguments. The log probability of each parameter is then 
;   summed to give the total prior log probability.
;
; Calling Sequence
; ----------------
;   ::
;
;       lnprior = lightning_priors(parameters, priors [, /error_check, lnprior_param=lnprior_param, $
;                                  in_bounds_idcs=in_bounds_idcs, out_bounds_idcs=out_bounds_idcs])
;
; Inputs
; ------
;   ``parameters`` : int, float, or double array(Nparam, Nmodels)
;       Parameters of the model(s). The actual parameters contained in this
;       array depend on the chosen model(s) during configuration.
;   ``priors`` : structure
;       A structure containing the prior hyper-parameters. (See
;       ``generate_prior_struct.pro`` for details and contents.)
;
; Optional Input
; --------------
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``lnprior`` : double array(Nmodels)
;       The total prior log probability of each model(s).
;
; Optional Outputs
; ----------------
;   ``lnprior_param`` : double array(Nparam, Nmodels)
;       The prior log probability for each parameter of each model.
;   ``in_bounds_idcs`` : int array(Nmodels)
;       The model indices that indicate if the given model(s) had all parameters
;       within the specified bounds.
;   ``out_bounds_idcs`` : int array(Nmodels)
;       The model indices that indicate if the given model(s) had one or more
;       parameters outside the specified bounds.
;
; Notes
; -----
;   `Uniform distribution equation <https://en.wikipedia.org/wiki/Continuous_uniform_distribution>`_ and
;   `Truncated normal distribution equation <https://en.wikipedia.org/wiki/Truncated_normal_distribution>`_
;
; Modification History
; --------------------
;   - 2022/06/02: Created (Keith Doore)
;-
 On_error, 2
 Compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(parameters) eq 0 then message, 'Variable is undefined: PARAMETERS.'
   if size(parameters, /type) lt 2 or size(parameters, /type) gt 5 then $
     message, 'PARAMETERS is not of type int, float, or double.'
   if size(parameters, /n_dim) gt 2 then $
     message, 'PARAMETERS must be a scalar, 1-D, or 2-D array.'
  
   if n_elements(priors) eq 0 then message, 'Variable is undefined: PRIORS.'
   if size(priors, /type) ne 8 then message, 'PRIORS is not of type structure.'
 endif


 param_size = size(parameters, /dim)
 Nparam = param_size[0]
 if size(parameters, /n_dim) le 1 then Nmodels = 1 else Nmodels = param_size[1]

; Check that all parameters for a given model are within the bounds, if not set
;   probability to 0
 lnprob_bounds = dblarr(Nparam, Nmodels)
 max_bound = rebin(priors.max_bound, Nparam, Nmodels)
 min_bound = rebin(priors.min_bound, Nparam, Nmodels)
 out_bounds_mask = parameters gt temporary(max_bound) or parameters lt temporary(min_bound)
 lnprob_bounds[where(out_bounds_mask, /null)] = alog(0.d0)
 in_bounds_idcs = where(total(out_bounds_mask, 1) eq 0, comp=out_bounds_idcs, /null)


; Compute analytical prior probability
;   lnprob = width *(x-mu)^2 + const
 width = rebin(priors.analytical.width, Nparam, Nmodels)
 mu = rebin(priors.analytical.mu, Nparam, Nmodels)
 const = rebin(priors.analytical.const, Nparam, Nmodels)
 lnprob_analytical = (temporary(width) * (parameters - temporary(mu))^2.d + temporary(const))


; Compute tabulated prior probability
 lnprob_tabulated = dblarr(Nparam, Nmodels)
 if total(finite(priors.idcs.tabulated)) gt 0 then begin
   for i=0, n_elements(priors.idcs.tabulated)-1 do begin
     tabulated_prior = priors.tabulated.((priors.idcs.tabulated)[i])

     ; Need '> 0' to prevent negative values due to any potential interpolation error
     lnprob_tabulated[(priors.idcs.tabulated)[i], *] = alog(interpol(tabulated_prior.pdf, $
                      tabulated_prior.values, parameters[(priors.idcs.tabulated)[i], *]) > 0)
   endfor
 endif


 lnprior_param =  lnprob_bounds + lnprob_analytical + lnprob_tabulated
; Set total log probability to a very negative number rather than -Infinity, since 
;   -Infinity-(-Infinity) = -NaN
 return, total(lnprior_param, 1) > (-1.d140)

end