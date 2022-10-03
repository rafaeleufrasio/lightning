function prior_sampled_initialization, priors, Nmodels, sigma_start=sigma_start
;+
; Name
; ----
;   PRIOR_SAMPLED_INITIALIZATION
;
; Purpose
; -------
;   Generates the initial parameter values from the prior, restricted to the
;   initialization range, for sampling with the MCMC or MPFIT algorithms.
;
; Calling Sequence
; ----------------
;   ::
;
;       parameter_start = prior_sampled_initialization(priors, Nmodels [, sigma_start=sigma_start])
;
; Inputs
; ------
;   ``priors`` : structure
;        A structure containing the prior hyper-parameters and initialization
;        range. (See ``generate_prior_struct.pro`` for details and contents.)
;   ``Nmodels`` : int, float, or double scalar
;        The number of model parameter sets to initialize.
;
; Output
; ------
;   ``parameter_start`` : double array(Nparam, Nmodels)
;       Initial parameter values for the model(s) as sampled from
;       each prior distribution within initialization range.
;
; Optional Output
; ---------------
;   ``sigma_start`` : double array(Nparam)
;       Initial proposal distribution sigma values for each parameter
;       as determined from half of the 16th and 84th percentile width
;       of each prior distribution.
;
; Modification History
; --------------------
;   - 2022/06/16: Created (Keith Doore)
;   - 2022/08/02: Updated to use initialization range to initialize within prior (Keith Doore)
;-
 On_error, 2
 Compile_opt idl2

; Error handling
 if n_elements(priors) eq 0 then message, 'Variable is undefined: PRIORS.'
 if size(priors, /type) ne 8 then message, 'PRIORS is not of type structure.'

 if n_elements(Nmodels) eq 0 then message, 'Variable is undefined: NMODELS.'
 if size(Nmodels, /type) lt 2 or size(Nmodels, /type) gt 5 then $
   message, 'NMODELS is not of type int, float, or double.'
 if size(Nmodels, /n_dim) ne 0 then message, 'NMODELS must be a scalar.'
 if min(Nmodels) le 0 then message, 'NMODELS must be a positive value.'


; Generate initial values
 parameter_start = dblarr(n_elements(priors.parameter_name), Nmodels)
 sigma_start = dblarr(n_elements(priors.parameter_name))


; Set fixed parameters to user specified value. Leave sigma set to 0.
 if total(finite(priors.idcs.fixed)) gt 0 then begin
   parameter_start[priors.idcs.fixed, *] = rebin((priors.fixed)[priors.idcs.fixed], n_elements(priors.idcs.fixed), Nmodels)
 endif

; Generate random sampling from initialization range within uniform priors.
;   Note: Should not need to check for infinity as prior_check.pro requires finite values.
 if total(finite(priors.idcs.uniform)) gt 0 then begin
   uniform_width = rebin(priors.INITIALIZE_RANGE[priors.idcs.uniform, 1]-priors.INITIALIZE_RANGE[priors.idcs.uniform, 0], $
                            n_elements(priors.idcs.uniform), Nmodels)
   uniform_min = rebin(priors.INITIALIZE_RANGE[priors.idcs.uniform, 0], n_elements(priors.idcs.uniform), Nmodels)
   parameter_start[priors.idcs.uniform, *] = randomu(seed, n_elements(priors.idcs.uniform), Nmodels) * $
                                             uniform_width + uniform_min

   ; Half the width between the 16th and 84th percentile of a uniform distribution is 0.34x the full width.
   sigma_start[priors.idcs.uniform] = 0.34d * (priors.max_bound-priors.min_bound)[priors.idcs.uniform]
 endif


; Generate random sampling from initialization range within normal priors.
; For normal distribution, have to interpolate from CDF since inverse erf does not exist in IDL.
;   This should not be an issue since we only do this once at the start.
 if total(finite(priors.idcs.normal)) gt 0 then begin
   for i=0, n_elements(priors.idcs.normal)-1 do begin
     normal_min = priors.INITIALIZE_RANGE[(priors.idcs.normal)[i], 0]
     normal_max = priors.INITIALIZE_RANGE[(priors.idcs.normal)[i], 1]
     normal_grid  = normal_min + dindgen(1e3) * (normal_max - normal_min)/(1.d3-1)

     width = (priors.analytical.width)[(priors.idcs.normal)[i]]
     mu = (priors.analytical.mu)[(priors.idcs.normal)[i]]
     ; Recompute constant with initialization range
     const = alog(2.d/sqrt(-1.d/width*!dpi) / (erf((normal_max-mu)/sqrt(-1.d/width)) - erf((normal_min-mu)/sqrt(-1.d/width))))
     normal_pdf_grid = exp(width * (normal_grid - mu)^2.d + const)

     normal_cdf_grid = total(normal_pdf_grid, /cum) / max(total(normal_pdf_grid, /cum))
     parameter_start[(priors.idcs.normal)[i], *] = interpol(normal_grid, normal_cdf_grid, randomu(seed, Nmodels))
   endfor
   ; Use sigma input by user for normal distribution (width = -0.5/sig^2 ->  sig = sqrt(-2*width))
   sigma_start[priors.idcs.normal] = sqrt(-2.d*(priors.analytical.width)[priors.idcs.normal])
 endif


; Generate random sampling from initialization range within tabulated priors.
; For tabulated distribution, have to interpolate since values are tabulated.
 if total(finite(priors.idcs.tabulated)) gt 0 then begin
   for i=0, n_elements(priors.idcs.tabulated)-1 do begin
     tabulated_prior = priors.tabulated.((priors.idcs.tabulated)[i])
     in_initialize_range = where(tabulated_prior.values ge priors.INITIALIZE_RANGE[(priors.idcs.tabulated)[i], 0] and $
                                 tabulated_prior.values le priors.INITIALIZE_RANGE[(priors.idcs.tabulated)[i], 1], /null)
     tabulated_cdf = total(tabulated_prior.pdf[in_initialize_range], /cum) / max(total(tabulated_prior.pdf[in_initialize_range], /cum))
     parameter_start[(priors.idcs.tabulated)[i], *] = $
       interpol(tabulated_prior.values[in_initialize_range], tabulated_cdf, randomu(seed, Nmodels))

     ; Calculate 16th and 84th percentiles and find half the difference for sigma
     percentile_values = interpol(tabulated_prior.values[in_initialize_range], tabulated_cdf, [0.16d, 0.84d])
     sigma_start[(priors.idcs.tabulated)[i]] = (percentile_values[1] - percentile_values[0])/2.d
   endfor
 endif


 return, parameter_start

end