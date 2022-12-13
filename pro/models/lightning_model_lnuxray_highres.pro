function lightning_model_lnuxray_highres, parameters, parameter_name, models, Lbol_AGN_model=Lbol_AGN_model, $
                                     L2500=L2500, xray_agn_model=xray_agn_model, agn_model=agn_model,  $
                                     error_check=error_check, get_wave=get_wave, wave=wave, $
                                     Lnu_xray_stellar_highres=Lnu_xray_stellar_highres, $
                                     Lnu_xray_unabs_stellar_highres=Lnu_xray_unabs_stellar_highres, $
                                     Lnu_xray_AGN_highres=Lnu_xray_AGN_highres, $
                                     Lnu_xray_unabs_AGN_highres=Lnu_xray_unabs_AGN_highres, Lnu_xray_mod=Lnu_xray_mod
;+
; Name
; ----
;   LIGHTNING_MODEL_LNUXRAY_HIGHRES
;
; Purpose
; -------
;   Generates the observed-frame luminosity density of a given lightning
;   X-ray model for a given set (or sets) of parameters.
;
; Calling Sequence
; ----------------
;   ::
;
;       Lnu_xray_highres = lightning_model_lnuxray_highres(parameters, parameter_name, models [, Lbol_AGN_model = ,$
;                                                     L2500 = , xray_agn_model = , agn_model = , /error_check, $
;                                                     /get_wave, wave=wave, Lnu_xray_stellar_highres=Lnu_xray_stellar_highres, $
;                                                     Lnu_xray_unabs_stellar_highres=Lnu_xray_unabs_stellar_highres, $
;                                                     Lnu_xray_AGN_highres=Lnu_xray_AGN_highres, $
;                                                     Lnu_xray_unabs_AGN_highres=Lnu_xray_unabs_AGN_highres, $
;                                                     Lnu_xray_mod=Lnu_xray_mod])
;
; Inputs
; ------
;   ``parameters`` : int, float, or double array(Nparam, Nmodels)
;       Parameters of the model(s). The actual parameters contained in this
;       array depend on the chosen model(s) during configuration.
;   ``parameter_name`` : string array(Nparam)
;       The names associated with the parameters of the model(s) given in the
;       same order as ``parameters``.
;   ``models`` : structure
;       A structure containing each model structure (stellar, dust, AGN,
;       X-ray) as a substructure. (See ``lightning_models.pro`` for details
;       and contents.)
;
; Optional Inputs
; ---------------
;   ``Lbol_AGN_model`` : double array(Nmodels)
;       Bolometric luminosity of the AGN model :math:`[L_\odot]`. (Required if using a power
;       law X-ray AGN emission model.)
;   ``L2500`` : double array(Nmodels)
;       The rest-frame 2500 Angstrom monochromatic luminosity shifted to the observed
;       frame :math:`[L_\odot\ {\rm Hz}^{-1}]`. (Required if using a power law X-ray
;       AGN emission model.)
;   ``agn_model`` : string scalar
;       The UV-to-IR AGN emission model to use. Current options are: ``'SKIRTOR'`` or ``'NONE'``.
;       (Default = ``'NONE'``)
;   ``xray_agn_model`` : string scalar
;       The X-ray AGN emission model to use. Current options are: ``'PLAW'``, ``'QSOSED'``, ``'NONE'``.
;       (Default = ``'QSOSED'``)
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;   ``wave`` : int, float, or double array(Nwave)
;       The common wavelength grid of the spectra :math:`[\mu \rm m]`.
;   ``get_wave`` : flag
;       If set, the function returns a value of ``-1`` after generating the default wavelength grid.
;       This is useful when one wants to get the default wavelength grid without computing
;       any high resolution spectra.
;
; Output
; ------
;   ``Lnu_xray_highres`` : double array(Nwave, Nmodels)
;       The total, high resolution X-ray spectral model, after absorption (intrinsic and
;       Galactic), calculated on a common wavelength grid for each set of parameters
;       :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;
; Optional Outputs
; ----------------
;   ``wave`` : double array(Nwave)
;       The common wavelength grid of the spectra :math:`[\mu \rm m]`.
;   ``Lnu_xray_stellar_highres`` : double array(Nwave, Nmodels)
;       Stellar emission component of ``Lnu_xray_highres`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   ``Lnu_xray_unabs_stellar_highres`` : double array(Nwave, Nmodels)
;       Intrinsic high resolution stellar X-ray spectral model calculated on a
;       common wavelength grid for each set of parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   ``Lnu_xray_AGN_highres`` : double array(Nwave, Nmodels)
;       AGN emission component of ``Lnu_xray_highres`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   ``Lnu_xray_unabs_AGN_highres`` : double array(Nwave, Nmodels)
;       Intrinsic high resolution AGN X-ray spectral model calculated on a common
;       wavelength grid for each set of parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   ``Lnu_xray_mod`` : double array(Nbands)
;       The mean luminosity calculated from ``Lnu_xray_highres`` using each of the X-ray
;       bandpasses :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;
; Modification History
; --------------------
;   - 2022/07/25: Created (Erik B. Monson)
;   - 2022/07/25: Renamed variables to match common naming scheme (Keith Doore)
;   - 2022/07/25: Added error handling with ``error_check`` keyword (Keith Doore)
;   - 2022/07/27: Added fixed observed-frame wavelength grid (Keith Doore)
;   - 2022/08/11: Renamed function to match function naming scheme (Keith Doore)
;   - 2022/09/08: Allow user-specified arbitary wavelength grid (Erik B. Monson)
;   - 2022/11/02: Galactic NH is now in units of 1e20 cm-2 (Erik B. Monson)
;   - 2022/12/13: Fixed bug in ``wave`` error check (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(parameters) eq 0 then message, 'Variable is undefined: PARAMETERS.'
   if size(parameters, /type) lt 2 or size(parameters, /type) gt 5 then $
     message, 'PARAMETERS must be of type int, float, or double.'
   if size(parameters, /n_dim) gt 2 then $
     message, 'PARAMETERS must be a scalar, 1-D, or 2-D array.'

   if n_elements(parameter_name) eq 0 then message, 'Variable is undefined: PARAMETER_NAME.'
   if size(parameter_name, /type) ne 7 then $
     message, 'PARAMETER_NAME must be of type string.'
   if size(reform(parameter_name), /n_dim) ne 1 then $
     message, 'PARAMETER_NAME must be a scalar or 1-D array.'

   if n_elements(models) eq 0 then message, 'Variable is undefined: MODELS.'
   if size(models, /type) ne 8 then message, 'MODELS must be of type structure.'

   if n_elements(agn_model) ne 0 then begin
     if size(agn_model, /type) ne 7 then message, 'AGN_MODEL must be of type string.'
     if size(agn_model, /dim) ne 0 then message, 'AGN_MODEL must be a scalar.'
     if total(strupcase(agn_model) eq ['SKIRTOR', 'NONE']) ne 1 then $
       message, "AGN_MODEL must be set to either 'SKIRTOR' or 'NONE'."
   endif

   if n_elements(xray_agn_model) ne 0 then begin
     if size(xray_agn_model, /type) ne 7 then message, 'XRAY_AGN_MODEL must be of type string.'
     if size(xray_agn_model, /dim) ne 0 then message, 'XRAY_AGN_MODEL must be a scalar.'
     if total(strupcase(xray_agn_model) eq ['PLAW', 'QSOSED', 'NONE']) ne 1 then $
       message, "XRAY_AGN_MODEL must be set to either 'PLAW', 'QSOSED', or 'NONE'."

     if strupcase(xray_agn_model) eq 'PLAW' then begin
       if n_elements(Lbol_AGN_model) eq 0 then $
         message, 'LBOL_AGN_MODEL must be specified if using power law XRAY_AGN_MODEL.'
       if n_elements(L2500) eq 0 then $
         message, 'L2500 must be specified if using power law XRAY_AGN_MODEL.'
     endif
   endif

   if n_elements(wave) ne 0 then begin
     if size(wave, /type) lt 2 or size(wave, /type) gt 5 then $
         message, 'WAVE must be of type int, float, or double.'
     if size(wave, /n_dim) gt 1 then message, 'WAVE must be a scalar or 1-D array.'
     if min(wave) lt 0 then $
         message, 'WAVE must only contain positive values.'
   endif
 endif

 if n_elements(agn_model) eq 0 then agn_model = 'NONE'
 if n_elements(xray_agn_model) eq 0 then xray_agn_model = 'QSOSED'

 if n_elements(wave) eq 0 then begin
     wave_logmin = -5.d
     wave_logmax = -2.d
     Nwave = 200.d
     wave = 10 ^ (wave_logmin + (wave_logmax - wave_logmin) / (Nwave - 1) * dindgen(Nwave))
 endif else begin
     Nwave = n_elements(wave)
 endelse

 if keyword_set(get_wave) then return, -1
 if size(parameters, /n_dim) le 1 then Nmodels = 1 else Nmodels = (size(parameters, /dim))[1]

 Nwave_mod = n_elements(models.xray_models.WAVE_REST)
 redshift = models.xray_models.REDSHIFT


; Extract parameters
 ; Transpose parameters array to eliminate need to reform to remove padded dimension after indexing.
 parameters_transposed = transpose(parameters)
 psi = transpose(parameters_transposed[*, where(strmatch(parameter_name, 'PSI_*'), /null)])
 tauv = parameters_transposed[*, where(parameter_name eq 'TAUV' or $
                                 parameter_name eq 'TAUV_DIFF', /null)]
 nh = parameters_transposed[*, where(parameter_name eq 'NH')]


 ; Handle absorption
 nh_milkyway = models.xray_models.GALACTIC_NH
 exp_neg_tau_milkyway = rebin(reform(models.xray_models.EXP_NEG_TAU_XRAY_MW, Nwave_mod, 1), Nwave_mod, Nmodels)^$
                        replicate(nh_milkyway, Nwave_mod, Nmodels)

 exp_neg_tau_nucleus = rebin(reform(models.xray_models.EXP_NEG_TAU_XRAY, Nwave_mod, 1), Nwave_mod, Nmodels)^$
                       rebin(reform(nh, 1, Nmodels), Nwave_mod, Nmodels)

 case strupcase(agn_model) of
   'SKIRTOR': begin
        nh_galaxy = 22.4 * (-2.5 * alog10(exp(-1 * tauV)))
        exp_neg_tau_galaxy = rebin(reform(models.xray_models.EXP_NEG_TAU_XRAY, Nwave_mod, 1), Nwave_mod, Nmodels)^$
                             rebin(reform(nh_galaxy, 1, Nmodels), Nwave_mod, Nmodels)
      end

   'NONE': exp_neg_tau_galaxy = exp_neg_tau_nucleus
 endcase


 case strupcase(agn_model) of
   'SKIRTOR': begin
        case strupcase(xray_agn_model) of
          'QSOSED': begin
               AGN_mass = parameters_transposed[*, where(parameter_name eq 'AGN_MASS', /null)]
               AGN_logmdot = parameters_transposed[*, where(parameter_name eq 'AGN_LOGMDOT', /null)]

               ; Countrate is ignored here. Lnu is output by reference.
               unabs_countrate_AGN = qsosed_spectrum(models.xray_models, AGN_mass, AGN_logmdot, $
                                                     lnu=Lnu_xray_unabs_AGN)
               Lnu_xray_AGN = exp_neg_tau_milkyway * exp_neg_tau_nucleus * Lnu_xray_unabs_AGN

               ; Interpolate AGN spectrum to common wavelength grid
               Lnu_xray_unabs_AGN_highres = dblarr(Nwave, Nmodels)
               Lnu_xray_AGN_highres = dblarr(Nwave, Nmodels)
               for i=0, Nmodels-1 do begin
                 Lnu_xray_unabs_AGN_highres[*, i] = interpol(Lnu_xray_unabs_AGN[*, i], models.xray_models.WAVE_OBS, wave)
                 Lnu_xray_AGN_highres[*, i] = interpol(Lnu_xray_AGN[*, i], models.xray_models.WAVE_OBS, wave)
               endfor
                Lnu_xray_unabs_AGN_highres[where(Lnu_xray_unabs_AGN_highres lt 0, /null)] = 0.d0
                Lnu_xray_AGN_highres[where(Lnu_xray_AGN_highres lt 0, /null)] = 0.d0
             end

          'PLAW': begin
               log_l_agn = parameters_transposed[*, where(parameter_name eq 'LOG_L_AGN', /null)]
               L2kev_agn = L2keV_LR17((10.d0 ^ log_l_agn) * L2500 / Lbol_AGN_model / (1 + redshift))

               Lnu_xray_unabs_AGN = (1 + redshift) * $
                                            rebin(reform(L2kev_agn, 1, Nmodels), Nwave_mod, Nmodels) *$
                                            rebin(reform(models.xray_models.LNU_AGN, Nwave_mod, 1), Nwave_mod, Nmodels)
               Lnu_xray_AGN = exp_neg_tau_milkyway * exp_neg_tau_nucleus * Lnu_xray_unabs_AGN

               ; Interpolate AGN spectrum to common wavelength grid
               Lnu_xray_unabs_AGN_highres = dblarr(Nwave, Nmodels)
               Lnu_xray_AGN_highres = dblarr(Nwave, Nmodels)
               for i=0, Nmodels-1 do begin
                 Lnu_xray_unabs_AGN_highres[*, i] = interpol(Lnu_xray_unabs_AGN[*, i], models.xray_models.WAVE_OBS, wave)
                 Lnu_xray_AGN_highres[*, i] = interpol(Lnu_xray_AGN[*, i], models.xray_models.WAVE_OBS, wave)
               endfor
                Lnu_xray_unabs_AGN_highres[where(Lnu_xray_unabs_AGN_highres lt 0, /null)] = 0.d0
                Lnu_xray_AGN_highres[where(Lnu_xray_AGN_highres lt 0, /null)] = 0.d0
             end

           'NONE': begin
               Lnu_xray_unabs_AGN_highres = dblarr(Nwave, Nmodels)
               Lnu_xray_AGN_highres = dblarr(Nwave, Nmodels)
             end
        endcase
      end

   'NONE': begin
       Lnu_xray_unabs_AGN_highres = dblarr(Nwave, Nmodels)
       Lnu_xray_AGN_highres = dblarr(Nwave, Nmodels)
     end
 endcase

 ; Get the normalizations for the X-rays from the stellar
 ; population. These are rest-frame 2-10 keV Lnu
 Lnu_XRB = gilbertson22_LX_tau(models.stellar_models, psi)
 Lnu_LMXB = rebin(reform(Lnu_XRB[0,*], 1, Nmodels), Nwave_mod, Nmodels)
 Lnu_HMXB = rebin(reform(Lnu_XRB[1,*], 1, Nmodels), Nwave_mod, Nmodels)

 Lnu_xray_unabs_stellar = (1 + redshift) * (Lnu_LMXB + Lnu_HMXB) * $
                                  rebin(reform(models.xray_models.LNU_XRB, Nwave_mod, 1), Nwave_mod, Nmodels)
 Lnu_xray_stellar = exp_neg_tau_milkyway * exp_neg_tau_galaxy * Lnu_xray_unabs_stellar


 ; Interpolate all highres models to common wavelength grid
 Lnu_xray_stellar_highres       = dblarr(Nwave, Nmodels)
 Lnu_xray_unabs_stellar_highres = dblarr(Nwave, Nmodels)
 for i=0, Nmodels-1 do begin
     Lnu_xray_stellar_highres[*, i]       = interpol(Lnu_xray_stellar[*, i], models.xray_models.WAVE_OBS, wave)
     Lnu_xray_unabs_stellar_highres[*, i] = interpol(Lnu_xray_unabs_stellar[*, i], models.xray_models.WAVE_OBS, wave)
 endfor
 Lnu_xray_stellar_highres[where(Lnu_xray_stellar_highres lt 0, /null)] = 0.d0
 Lnu_xray_unabs_stellar_highres[where(Lnu_xray_unabs_stellar_highres lt 0, /null)] = 0.d0


 Lnu_xray_highres = Lnu_xray_AGN_highres + Lnu_xray_stellar_highres


 ; Get the average Lnu in each bandpass if we want it
 if arg_present(Lnu_xray_mod) then begin
     xray_bandpass = models.xray_models.XRAY_BANDPASS
     Nxray = n_elements(models.xray_models.XRAY_EXPOSURE)

     Lnu_xray_mod = dblarr(Nxray, Nmodels)

     hc = (!lightning_cgs.hplanck / !lightning_cgs.keV) * (1e4 * !lightning_cgs.clight) ; h * c in keV * micron
     E_obs = hc / wave

     for j=0, Nmodels - 1 do begin
         for i=0, Nxray - 1 do begin
             filter = 0.d0 * E_obs
             filter[where(E_obs ge xray_bandpass[0, i] and E_obs le xray_bandpass[1, i], /null)] = 1.d0
             filter = filter / trapez(filter, wave)

             Lnu_xray_mod[i, j] = trapez(filter * Lnu_xray_highres[*, j], wave)
         endfor
     endfor
 endif

 return, Lnu_xray_highres

end
