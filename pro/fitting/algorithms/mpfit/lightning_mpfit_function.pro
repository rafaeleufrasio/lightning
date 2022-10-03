function lightning_mpfit_function, parameters, sed_data=sed_data, config_nopriors=config_nopriors, $
                                   models=models, priors=priors, error_check=error_check
;+
; Name
; ----
;   LIGHTNING_MPFIT_FUNCTION
;
; Purpose
; -------
;   Calculates the deviates to be minimized by MPFIT given the data
;   and a set of parameters.
;
; Calling Sequence
; ----------------
;   ::
;
;       deviates = lightning_mpfit_function(parameters [, sed_data = , config_nopriors = , $
;                                           models = , priors = , /error_check])
;
; Inputs
; ------
;   ``parameters`` : int, float, or double array(Nparam)
;       Parameters of the model. The actual parameters contained in this
;       array depend on the chosen model during configuration.
;
; Optional Inputs
; ---------------
;   ``sed_data`` : structure
;       A structure containing the SED luminosities and uncertainties, filter
;       labels, distances, redshifts, and optional X-ray data. (See 
;       ``lightning_input.pro`` for details and contents.)
;   ``config_nopriors`` : structure
;       A Lightning configuration structure edited to remove the prior 
;       substructures for each parameter. (See ``lightning_configure_defaults.pro``
;       for details and contents.)
;   ``models`` : structure
;       A structure containing each model structure (stellar, dust, AGN, 
;       X-ray) as a substructure. (See ``lightning_models.pro`` for details
;       and contents.)
;   ``priors`` : structure
;       A structure containing the prior hyper-parameters. (See
;       ``generate_prior_struct.pro`` for details and contents.)
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``deviates`` : float or double array(Nfilters + Nxray)
;       The deviates of the model for the specified parameters.
;
; Modification History
; --------------------
;   - 2022/08/15: Created (Keith Doore)
;   - 2022/09/21: Updated to allow for X-ray fluxes (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(parameters) eq 0 then message, 'Variable is undefined: PARAMETERS.'
   if size(parameters, /type) lt 2 or size(parameters, /type) gt 5 then $
     message, 'PARAMETERS must be of type int, float, or double.'
   if size(parameters, /n_dim) gt 1 then $
     message, 'PARAMETERS must be a scalar or 1-D array.'

   if n_elements(sed_data) eq 0 then message, 'Variable is undefined: SED_DATA.'
   if size(sed_data, /type) ne 8 then message, 'SED_DATA must be of type structure.'

   if n_elements(config_nopriors) eq 0 then message, 'Variable is undefined: CONFIG_NOPRIORS.'
   if size(config_nopriors, /type) ne 8 then message, 'CONFIG_NOPRIORS must be of type structure.'

   if n_elements(models) eq 0 then message, 'Variable is undefined: MODELS.'
   if size(models, /type) ne 8 then message, 'MODELS must be of type structure.'

   if n_elements(priors) eq 0 then message, 'Variable is undefined: PRIORS.'
   if size(priors, /type) ne 8 then message, 'PRIORS must be of type structure.'
 endif


; Call lightning_model_lnprob to generate lnu_mod and xray_counts_mod
 lnprob = lightning_model_lnprob(sed_data, parameters, config_nopriors, models, priors, $
                                 Lnu_mod=Lnu_mod, xray_counts_mod=xray_counts_mod, $
                                 Lnu_xray_mod=Lnu_xray_mod)
 ; Remove potential extra padded dimension
 Lnu_mod = reform(Lnu_mod)
 if n_elements(xray_counts_mod) gt 0 then xray_counts_mod = reform(xray_counts_mod)
 if n_elements(Lnu_xray_mod) gt 0 then Lnu_xray_mod = reform(Lnu_xray_mod)


; Calculate deviates
;   Models will alway be non-NaN as MPFIT keeps values within bound limits
 ; If Lnu_unc = 0 for a band, sigma_total should remain 0 or be NaN, since that band is to be ignored.
 ;   unc_buffer will be -NaN if Lnu_unc = 0 or 1 if Lnu_unc ne 0.
 unc_buffer = sed_data.Lnu_unc/sed_data.Lnu_unc
 sigma_total = sqrt(sed_data.Lnu_unc^2 + (config_nopriors.MODEL_UNC * Lnu_mod * unc_buffer)^2)
 UVIR_deviates = (sed_data.Lnu_obs - Lnu_mod)/sigma_total
 ; Exclude from the deviates, measurements with uncertainty eq 0 (i.e., deviate is NaN)
 UVIR_deviates = UVIR_deviates[where(finite(UVIR_deviates), /null)]


 ; Generate X-ray model counts and compute X-ray deviates statistic
 if config_nopriors.XRAY_EMISSION then begin

   case strupcase(config_nopriors.XRAY_UNIT) of
     'COUNTS': begin
         xray_counts = sed_data.NET_COUNTS
         xray_counts_unc = sed_data.NET_COUNTS_UNC

         ; Calculate X-ray chi2 contribution.
         ;   Rather than attempting to use one of the Poisson likelihoods (which are not
         ;   appropriate for use with background-subtracted data), we use a chi2
         ;   with the X-ray counts and the selected or user supplied count uncertainty.
         xray_deviates = (xray_counts_mod - xray_counts) / xray_counts_unc
     end
     'FLUX': begin
         Lnu_obs_xray = sed_data.XRAY_LNU_OBS
         Lnu_unc_xray = sed_data.XRAY_LNU_UNC
         xray_unc_buffer = Lnu_unc_xray/Lnu_unc_xray

         Lnu_unc_xray_total = sqrt(Lnu_unc_xray^2 + (config_nopriors.MODEL_UNC * Lnu_xray_mod * xray_unc_buffer)^2)
         xray_deviates = (Lnu_xray_mod - Lnu_obs_xray) / Lnu_unc_xray_total
     end
   endcase
 endif else begin
   xray_deviates = !null
 endelse

; To implement the effect of normal priors, we use regularization. (param - mu)^2/sigma
 if total(finite(priors.idcs.normal)) gt 0 then begin
   regularize = (parameters - priors.ANALYTICAL.MU)^2.d/sqrt(-0.5d/priors.ANALYTICAL.WIDTH)
   regularize = total(regularize[priors.IDCS.NORMAL])
 endif else regularize = 0.d


 deviates = [UVIR_deviates, xray_deviates]

 return, sqrt(deviates^2 + regularize)

end