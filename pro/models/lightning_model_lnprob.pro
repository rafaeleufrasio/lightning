function lightning_model_lnprob, sed_data, parameters, config_nopriors, models, priors, $
                                 negative=negative, error_check=error_check, $
                                 UVIR_chi2=UVIR_chi2, xray_chi2=xray_chi2,$
                                 lnlike=lnlike, lnprior=lnprior, Lnu_mod=Lnu_mod, $
                                 xray_counts_mod=xray_counts_mod, Lnu_xray_mod=Lnu_xray_mod
;+
; Name
; ----
;   LIGHTNING_MODEL_LNPROB
;
; Purpose
; -------
;   Calculates the natural log of the model(s) probability
;   (i.e. ln(prior * likelihood)) given the data for a set
;   of parameters.
;
; Calling Sequence
; ----------------
;   ::
;
;       lnprob = lightning_model_lnprob(sed_data, parameters, config_nopriors, models, priors [, $
;                                       /negative, /error_check, UVIR_chi2=UVIR_chi2, $
;                                       xray_chi2=xray_chi2, lnlike=lnlike, lnprior=lnprior, $
;                                       Lnu_mod=Lnu_mod, xray_counts_mod=xray_counts_mod, $
;                                       Lnu_xray_mod=Lnu_xray_mod])
;
; Inputs
; ------
;   ``sed_data`` : structure
;       A structure containing the SED luminosities and uncertainties, filter
;       labels, distances, redshifts, and optional X-ray data. (See
;       ``lightning_input.pro`` for details and contents.)
;   ``parameters`` : int, float, or double array(Nparam, Nmodels)
;       Parameters of the model(s). The actual parameters contained in this
;       array depend on the chosen model(s) during configuration.
;   ``config_nopriors`` : structure
;       A Lightning configuration structure edited to remove the prior
;       substructures for each parameter. (See ``lightning_configure_defaults.pro``
;       for details and contents.)
;   ``models`` : structure
;       A structure containing each model structure (stellar, dust, AGN,
;       X-ray) as a substructure. (See ``lightning_models.pro`` for details
;       and contents.)
;   ``priors`` : structure
;        A structure containing the prior hyper-parameters. (See
;        ``generate_prior_struct.pro`` for details and contents.)
;
; Optional Inputs
; ---------------
;   ``negative`` : flag
;       If set, the function returns the negative of the log probability
;       (for minimization purposes).
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``lnprob`` : double array(Nmodels)
;       The log probability of the model(s) for the specified parameters.
;
; Optional Outputs
; ----------------
;   ``UVIR_chi2`` : double array(Nmodels)
;       The :math:`\chi^2` of the UV-to-IR portion of the model(s).
;   ``xray_chi2`` : double array(Nmodels)
;       The X-ray statistic of the model(s).
;   ``lnlike`` : double array(Nmodels)
;       The total likelihood log probability of the model(s).
;   ``lnprior`` : double array(Nmodels)
;       The total prior log probability of each model(s).
;   ``Lnu_mod`` : double array(Nfilters, Nmodels)
;       The model Lnu for each filter of each model generated from the
;       specified parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   ``xray_counts_mod`` : double array(Nxray, Nmodels)
;       The model X-ray counts for each band of each model generated from
;       the specified parameters :math:`[{\rm counts}]`.
;   ``Lnu_xray_mod`` : double array(Nxray, Nmodels)
;       The X-ray model Lnu for each band of each model generated from the
;       specified parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;
; Modification History
; --------------------
;   - 2022/02/01: Created (Erik Monson)
;   - 2022/06/07: Major update to include new implementation (e.g., prior, config, etc.) (Keith Doore)
;   - 2022/07/05: Renamed config to ``config_nopriors`` to reflect removed parameters (Keith Doore)
;   - 2022/07/25: Fixed bug where in bounds indexing for ``AGN_mass`` and ``AGN_logmdot`` in qsosed ``L2500`` computation was missing (Keith Doore)
;   - 2022/09/01: Added handling for user-supplied X-ray count uncertainties (Erik B. Monson)
;   - 2022/09/14: Updates to allow fitting with X-ray fluxes (Erik B. Monson)
;-
 On_error, 2
 Compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(sed_data) eq 0 then message, 'Variable is undefined: SED_DATA.'
   if size(sed_data, /type) ne 8 then message, 'SED_DATA is not of type structure.'

   if n_elements(parameters) eq 0 then message, 'Variable is undefined: PARAMETERS.'
   if size(parameters, /type) lt 2 or size(parameters, /type) gt 5 then $
     message, 'PARAMETERS is not of type int, float, or double.'
   if size(parameters, /n_dim) gt 2 then $
     message, 'PARAMETERS must be a scalar, 1-D, or 2-D array.'

   if n_elements(config_nopriors) eq 0 then message, 'Variable is undefined: CONFIG_NOPRIORS.'
   if size(config_nopriors, /type) ne 8 then message, 'CONFIG_NOPRIORS is not of type structure.'

   if n_elements(models) eq 0 then message, 'Variable is undefined: MODELS.'
   if size(models, /type) ne 8 then message, 'MODELS is not of type structure.'

   if n_elements(priors) eq 0 then message, 'Variable is undefined: PRIORS.'
   if size(priors, /type) ne 8 then message, 'PRIORS is not of type structure.'
 endif


; Generate log prior probability
 lnprior = lightning_priors(parameters, priors, in_bounds_idcs=in_bounds_idcs, out_bounds_idcs=out_bounds_idcs)
 count_in_bounds = n_elements(in_bounds_idcs)


; Generate UV-IR Lnu for model
 Nfilters = n_elements(sed_data.Lnu_obs)
 if size(parameters, /n_dim) le 1 then Nmodels = 1 else Nmodels = (size(parameters, /dim))[1]
 Lnu_obs = rebin(sed_data.Lnu_obs, Nfilters, Nmodels)
 Lnu_unc = rebin(sed_data.Lnu_unc, Nfilters, Nmodels)
 Lnu_mod = replicate(!values.D_NaN, Nfilters, Nmodels)
 if (count_in_bounds gt 0) then begin

   ; Get the qsosed L2500 if need be, so that the AGN model luminosity is correct.
   L2500 = !null
   if config_nopriors.XRAY_EMISSION and (strupcase(config_nopriors.AGN_MODEL) eq 'SKIRTOR') then begin
     if strupcase(config_nopriors.XRAY_AGN_MODEL) eq 'QSOSED' then begin
       ; Need to make `where` index a scalar, since can not index by two different size arrays.
       agn_mass    = reform(parameters[(where(priors.parameter_name eq 'AGN_MASS', /null))[0], in_bounds_idcs])
       agn_logmdot = reform(parameters[(where(priors.parameter_name eq 'AGN_LOGMDOT', /null))[0], in_bounds_idcs])

       L2500 = qsosed_L2500(models.xray_models, agn_mass, agn_logmdot)
     endif
   endif

   Lnu_mod[*, in_bounds_idcs] = lightning_model_lnu(parameters[*, in_bounds_idcs], priors.parameter_name, models, $
                                                    Lbol_AGN_model=Lbol_AGN_model, L2500=L2500, _extra=config_nopriors)
 endif


; Calculate chi2
 ; If Lnu_unc = 0 for a band, sigma_total should remain 0 or be NaN, since that band is to be ignored.
 ;   unc_buffer will be -NaN if Lnu_unc = 0 or 1 if Lnu_unc ne 0.
 unc_buffer = Lnu_unc/Lnu_unc
 sigma_total = sqrt(Lnu_unc^2 + (config_nopriors.MODEL_UNC * Lnu_mod * unc_buffer)^2)

 ; NaNs occur from non-computed models or sigma_total = 0
 UVIR_chi2 = total((Lnu_obs - Lnu_mod)^2 / sigma_total^2, 1, /NaN)
 ; Set chi2 of non-computed models to NaN, since chi2 can’t be computed
 UVIR_chi2[out_bounds_idcs] = !values.D_NaN


; Generate X-ray model counts and compute X-ray statistic
 if config_nopriors.XRAY_EMISSION then begin

   case strupcase(config_nopriors.XRAY_UNIT) of
       'COUNTS': begin
           Nxray = n_elements(sed_data.NET_COUNTS)
           xray_counts = rebin(sed_data.NET_COUNTS, Nxray, Nmodels)
           xray_counts_unc = rebin(sed_data.NET_COUNTS_UNC, Nxray, Nmodels)

           xray_counts_mod = replicate(!values.D_NaN, Nxray, Nmodels)
           if (count_in_bounds gt 0) then begin
             xray_counts_mod[*, in_bounds_idcs] = lightning_model_counts(parameters[*, in_bounds_idcs], priors.parameter_name, models, $
                                                                         Lbol_AGN_model=Lbol_AGN_model, L2500=L2500, $
                                                                         _extra=config_nopriors)
           endif

           ; Calculate X-ray chi2 contribution.
           ;   Rather than attempting to use one of the Poisson likelihoods (which are not
           ;   appropriate for use with background-subtracted data), we use a chi2
           ;   with the X-ray counts and the selected or user supplied count uncertainty.
           xray_chi2 = total((xray_counts_mod - xray_counts)^2 / xray_counts_unc^2, 1, /NaN)
       end
       'FLUX': begin
           Nxray = n_elements(sed_data.XRAY_LNU_OBS)
           Lnu_obs_xray = rebin(sed_data.XRAY_LNU_OBS, Nxray, Nmodels)
           Lnu_unc_xray = rebin(sed_data.XRAY_LNU_UNC, Nxray, Nmodels)
           xray_unc_buffer = Lnu_unc_xray/Lnu_unc_xray


           Lnu_xray_mod = replicate(!values.D_NaN, Nxray, Nmodels)
           if (count_in_bounds gt 0) then begin
             Lnu_xray_mod[*, in_bounds_idcs] = lightning_model_lnuxray(parameters[*, in_bounds_idcs], priors.parameter_name, models, $
                                                                       Lbol_AGN_model=Lbol_AGN_model, L2500=L2500, $
                                                                       _extra=config_nopriors)
           endif

           Lnu_unc_xray_total = sqrt(Lnu_unc_xray^2 + (config_nopriors.MODEL_UNC * Lnu_xray_mod * xray_unc_buffer)^2)
           xray_chi2 = total((Lnu_xray_mod - Lnu_obs_xray)^2 / Lnu_unc_xray_total^2, 1, /NaN)
       end
   endcase

   xray_chi2[out_bounds_idcs] = !values.D_NaN

 endif else begin
   xray_chi2 = 0.0
 endelse


 chi2 = UVIR_chi2 + xray_chi2
 lnlike = -0.5 * chi2

; Calculate total log probability.
;   Ignore NaNs from likelihood, due to model not being computed.
;   Prior values is hugely negative if out of bounds, so rejection is guaranteed.
 lnprob = total([[lnlike], [lnprior]], 2, /NaN)

 if keyword_set(negative) then begin
   return, -1 * lnprob
 endif else begin
   return, lnprob
 endelse

end
