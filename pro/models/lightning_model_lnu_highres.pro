function lightning_model_lnu_highres, parameters, parameter_name, models, ssp=ssp, sfh=sfh, $
                                      dust_model=dust_model, agn_model=agn_model, energy_balance=energy_balance, $
                                      error_check=error_check, get_wave=get_wave, wave=wave, $
                                      Lnu_stellar_highres=Lnu_stellar_highres, Lnu_unred_stellar_highres=Lnu_unred_stellar_highres, $
                                      Lnu_AGN_highres=Lnu_AGN_highres, Lnu_unred_AGN_highres=Lnu_unred_AGN_highres, $
                                      Lnu_dust_highres=Lnu_dust_highres, LTIR=LTIR, Lbol_AGN_model=Lbol_AGN_model, $
                                      L2500=L2500, _extra=_extra
;+
; Name
; ----
;   LIGHTNING_MODEL_LNU_HIGHRES
;
; Purpose
; -------
;   Generates the observed-frame luminosity spectrum of a given lightning model
;   for a given set (or sets) of parameters.
;
; Calling Sequence
; ----------------
;   ::
;
;       Lnu_mod_highres = lightning_model_lnu_highres(parameters, parameter_name, models [, ssp = , $
;                                                     sfh = , dust_model = , agn_model = , /energy_balance, $
;                                                     /error_check, /get_wave, Lnu_stellar_highres=Lnu_stellar_highres, $
;                                                     wave=wave, Lnu_unred_stellar_highres=Lnu_unred_stellar_highres, $
;                                                     Lnu_AGN_highres=Lnu_AGN_highres, Lnu_unred_AGN_highres=Lnu_unred_AGN_highres, $
;                                                     Lnu_dust_highres=Lnu_dust_highres, LTIR=LTIR, $
;                                                     Lbol_AGN_model=Lbol_AGN_model, L2500=L2500, _extra = ])
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
;   ``ssp`` : string scalar
;       The simple stellar population (SSP) models to use for the stellar
;       population. Current options  are: ``'PEGASE'`` or ``'NONE'``.
;       (Default = ``'PEGASE'``)
;   ``sfh`` : string scalar
;       The type of SFH to assume if the model is to include a stellar
;       population. Current options are: ``'NON-PARAMETRIC'``.
;       (Default = ``'NON-PARAMETRIC'``)
;   ``dust_model`` : string scalar
;       The dust emission model to use. Current options are: ``'DL07'`` or ``'NONE'``.
;       (Default = ``'DL07'``)
;   ``agn_model`` : string scalar
;       The UV-to-IR AGN emission model to use. Current options are: ``'SKIRTOR'`` or ``'NONE'``.
;       (Default = ``'NONE'``)
;   ``energy_balance`` : flag
;       If set, the total integrated IR luminosity (normalization) of the dust emission
;       is tied to the total absorbed stellar (and, if set, AGN) emission.
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;   ``wave`` : int, float, or double array(Nwave)
;       The common wavelength grid of the spectra :math:`[\mu \rm m]`.
;   ``get_wave`` : flag
;       If set, the function returns a value of ``-1`` after generating the default wavelength grid.
;       This is useful when one wants to get the default wavelength grid without computing
;       any high resolution spectra.
;   ``L2500`` : double array(Nmodels)
;       The rest-frame 2500 Angstrom monochromatic luminosity shifted to the observed
;       frame :math:`[L_\odot\ {\rm Hz}^{-1}]`. (Used as input if using a QSOSED X-ray
;       AGN model.)
;   ``_extra`` : structure
;       Additional optional inputs that are passed to ``binned_stellar_spectrum.pro`` and
;       ``skirtor_spectrum.pro``.
;
; Output
; ------
;   ``Lnu_mod_highres`` : double array(Nwave, Nmodels)
;       The model luminosity spectrum of each model generated on a common wavelength
;       grid from the specified parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;
; Optional Outputs
; ----------------
;   ``wave`` : double array(Nwave)
;       The common wavelength grid of the spectra :math:`[\mu \rm m]`.
;   ``Lnu_stellar_highres`` : double array(Nwave, Nmodels)
;       Stellar emission component of ``Lnu_mod_highres`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   ``Lnu_unred_stellar_highres`` : double array(Nfilters, Nmodels)
;       The unattenuated stellar emission of each model generated on a common wavelength
;       grid from the specified parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   ``Lnu_AGN_highres`` : double array(Nwave, Nmodels)
;       AGN emission component of ``Lnu_mod_highres`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   ``Lnu_unred_AGN_highres`` : double array(Nwave, Nmodels)
;       The unattenuated AGN emission of each model generated on a common wavelength
;       grid from the specified parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   ``Lnu_dust_highres`` : double array(Nwave, Nmodels)
;       Dust emission component of ``Lnu_mod_highres`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   ``LTIR`` : double array(Nmodels)
;       The total integrated IR luminosity (i.e., the bolometric luminosity)
;       of the dust model for each set of model parameters :math:`[L_\odot]`.
;   ``Lbol_AGN_model`` : double array(Nmodels)
;       Bolometric luminosity of the AGN model :math:`[L_\odot]`.
;   ``L2500`` : double array(Nmodels)
;       The rest-frame 2500 Angstrom monochromatic luminosity shifted to the observed
;       frame for the current AGN model or qsosed X-ray AGN model if input
;       :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;
; Modification History
; --------------------
;   - 2022/04/14: Created (Erik B. Monson)
;   - 2022/07/06: Major update to match ``lightning_model_lnu`` (Keith Doore)
;   - 2022/07/22: Added optional output of unreddened ``Lnu_stellar_highres`` (Keith Doore)
;   - 2022/07/22: Added optional output of the common wavelength grid (Keith Doore)
;   - 2022/07/22: Added optional output of unreddened ``Lnu_AGN_highres`` (Keith Doore)
;   - 2022/07/27: Renamed unreddened Lnus to prevent ambiguous keywords (Keith Doore)
;   - 2022/09/08: Allow user-specified arbitary wavelength grid (Erik B. Monson)
;   - 2022/10/25: Renamed SPS to SSP (Keith Doore)
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

   if n_elements(ssp) ne 0 then begin
     if size(ssp, /type) ne 7 then message, 'SSP must be of type string.'
     if size(ssp, /dim) ne 0 then message, 'SSP must be a scalar.'
     if total(strupcase(ssp) eq ['PEGASE', 'NONE']) ne 1 then $
       message, "SSP must be set to either 'PEGASE' or 'NONE'."
   endif

   if n_elements(sfh) ne 0 then begin
     if size(sfh, /type) ne 7 then message, 'SFH must be of type string.'
     if size(sfh, /dim) ne 0 then message, 'SFH must be a scalar.'
     if total(strupcase(sfh) eq ['NON-PARAMETRIC']) ne 1 then $
       message, "SFH must be set to either 'NON-PARAMETRIC'."
   endif

   if n_elements(dust_model) ne 0 then begin
     if size(dust_model, /type) ne 7 then message, 'DUST_MODEL must be of type string.'
     if size(dust_model, /dim) ne 0 then message, 'DUST_MODEL must be a scalar.'
     if total(strupcase(dust_model) eq ['DL07', 'NONE']) ne 1 then $
       message, "DUST_MODEL must be set to either 'DL07' or 'NONE'."
   endif

   if n_elements(agn_model) ne 0 then begin
     if size(agn_model, /type) ne 7 then message, 'AGN_MODEL must be of type string.'
     if size(agn_model, /dim) ne 0 then message, 'AGN_MODEL must be a scalar.'
     if total(strupcase(agn_model) eq ['SKIRTOR', 'NONE']) ne 1 then $
       message, "AGN_MODEL must be set to either 'SKIRTOR' or 'NONE'."
   endif

   if n_elements(L2500) ne 0 then begin
     if size(L2500, /type) lt 2 or size(L2500, /type) gt 5 then $
       message, 'L2500 must be of type int, float, or double.'
     if size(L2500, /dim) gt 1 then message, 'L2500 must be a scalar or 1-D array.'
     if min(L2500) lt 0 then $
       message, "L2500 must only contain positive values."
   endif

   if n_elements(wave) ne 0 then begin
     if size(wave, /type) lt 2 or size(wave, /type) gt 5 then $
         message, 'WAVE must be of type int, float, or double.'
     if size(wave, /n_dim) gt 1 then message, 'WAVE must be a scalar or 1-D array.'
     if min(wave) lt 0 then $
         message, 'WAVE must only contain positive values.'
   endif
 endif

 if n_elements(ssp) eq 0 then ssp = 'PEGASE'
 if n_elements(sfh) eq 0 then sfh = 'NON-PARAMETRIC'
 if n_elements(dust_model) eq 0 then dust_model = 'DL07'
 if n_elements(agn_model) eq 0 then agn_model = 'NONE'


 if n_elements(wave) eq 0 then begin
     ; Generate a common wavelength grid
     ;   Make a logarithmically spaced grid from just shy of the Lyman break to ~1585 microns
     wave_logmin = -2.d
     wave_logmax = 3.2d
     Nwave = 1000.d
     wave = 10 ^ (wave_logmin + (wave_logmax - wave_logmin) / (Nwave - 1) * dindgen(Nwave))
 endif else begin
     Nwave = n_elements(wave)
 endelse

 if keyword_set(get_wave) then return, -1


 if size(parameters, /n_dim) le 1 then Nmodels = 1 else Nmodels = (size(parameters, /dim))[1]

 ; Transpose parameters array to eliminate need to reform to remove padded dimension after indexing.
 parameters_transposed = transpose(parameters)


; Generate attenuated stellar emission
 case strupcase(ssp) of
   'PEGASE': begin
        case strupcase(sfh) of
          'NON-PARAMETRIC': begin
               psi = transpose(parameters_transposed[*, where(strmatch(parameter_name, 'PSI_*'), /null)])

               tauv_diff = parameters_transposed[*, where(parameter_name eq 'TAUV' or $
                                               parameter_name eq 'TAUV_DIFF', /null)]
               delta     = parameters_transposed[*, where(parameter_name eq 'DELTA', /null)]
               tauv_bc   = parameters_transposed[*, where(parameter_name eq 'TAUV_BC', /null)]

               taub_f  = parameters_transposed[*, where(parameter_name eq 'TAUB_F', /null), *]
               f_clump = parameters_transposed[*, where(parameter_name eq 'F_CLUMP', /null), *]
               cosi    = parameters_transposed[*, where(parameter_name eq 'COSI', /null), *]
               b_to_d  = parameters_transposed[*, where(parameter_name eq 'B_TO_D', /null), *]

               Lnu_stellar = binned_stellar_spectrum(models.stellar_models, psi, $
                                                     tauV_DIFF=tauV_DIFF, delta=delta, tauV_BC=tauV_BC, $
                                                     tauB_f=tauB_f, F_clump=F_clump, cosi=cosi, b_to_d=b_to_d, $
                                                     atten_models=models.atten_models, Lbol_abs_stellar=Lbol_abs_stellar, $
                                                     Lnu_spec_unred_stellar=Lnu_unred_stellar, _extra=_extra)

               ; Interpolate stellar spectrum to common wavelength grid
               Lnu_stellar_highres = dblarr(Nwave, Nmodels)
               Lnu_unred_stellar_highres = dblarr(Nwave, Nmodels)
               for i=0, Nmodels-1 do begin
                   Lnu_stellar_highres[*, i] = interpol(Lnu_stellar[*, i], models.stellar_models.WAVE_OBS, wave)
                   Lnu_unred_stellar_highres[*, i] = interpol(Lnu_unred_stellar[*, i], models.stellar_models.WAVE_OBS, wave)
               endfor
               Lnu_stellar_highres[where(Lnu_stellar_highres lt 0, /null)] = 0.d0
               Lnu_unred_stellar_highres[where(Lnu_unred_stellar_highres lt 0, /null)] = 0.d0
             end

          'PARAMETRIC': begin
              Lnu_stellar_highres = dblarr(Nwave, Nmodels)
              Lnu_unred_stellar_highres = dblarr(Nwave, Nmodels)
              Lbol_abs_stellar = dblarr(Nmodels)
            end
        endcase
      end

   'NONE': begin
       Lnu_stellar_highres = dblarr(Nwave, Nmodels)
       Lnu_unred_stellar_highres = dblarr(Nwave, Nmodels)
       Lbol_abs_stellar = dblarr(Nmodels)
     end
 endcase




; Generate AGN emission
 case strupcase(agn_model) of
   'SKIRTOR': begin
        tau97    = parameters_transposed[*, where(parameter_name eq 'TAU97', /null)]
        agn_cosi = parameters_transposed[*, where(parameter_name eq 'AGN_COSI', /null)]

        tauv_diff = parameters_transposed[*, where(parameter_name eq 'TAUV' or $
                                        parameter_name eq 'TAUV_DIFF', /null)]
        delta     = parameters_transposed[*, where(parameter_name eq 'DELTA', /null)]

        ; Since the models are gridded in terms of inclination angle, transform back from cosi
        i_agn = 180.d0 / !dpi * acos(agn_cosi)

        Lnu_spec_AGN = skirtor_spectrum(models.agn_models, tau97=tau97, i_agn=i_agn, tauV_diff=tauV_diff,$
                                        delta=delta, Lbol_abs_AGN=Lbol_abs_SKIRTOR, Lbol_AGN=Lbol_SKIRTOR, $
                                        L2500=L2500_SKIRTOR, Lnu_spec_unred_AGN=Lnu_spec_unred_AGN, _extra=_extra)

        ; Scale AGN emission
        ;   Check if L2500 is input from QSOSED model. If so, use it for scaling
        if n_elements(L2500) ne 0 then begin
          scale_factor = L2500/L2500_SKIRTOR
        endif else begin
          Lbol_AGN = 10.d0^parameters_transposed[*, where(parameter_name eq 'LOG_L_AGN', /null)]
          scale_factor = Lbol_AGN/Lbol_SKIRTOR
          L2500 = L2500_SKIRTOR
        endelse
        Lbol_AGN_model = Lbol_SKIRTOR

        Lnu_AGN = rebin(reform(scale_factor, 1, Nmodels), n_elements(models.agn_models.WAVE_OBS), Nmodels) * Lnu_spec_AGN
        Lnu_unred_AGN = rebin(reform(scale_factor, 1, Nmodels), n_elements(models.agn_models.WAVE_OBS), Nmodels) * Lnu_spec_unred_AGN
        Lbol_abs_AGN = scale_factor * Lbol_abs_SKIRTOR

        ; Interpolate AGN spectrum to common wavelength grid
        Lnu_agn_highres = dblarr(Nwave, Nmodels)
        Lnu_unred_AGN_highres = dblarr(Nwave, Nmodels)
        for i=0, Nmodels-1 do begin
            Lnu_agn_highres[*, i] = interpol(Lnu_AGN[*, i], models.agn_models.WAVE_OBS, wave)
            Lnu_unred_AGN_highres[*, i] = interpol(Lnu_unred_AGN[*, i], models.agn_models.WAVE_OBS, wave)
        endfor
        Lnu_agn_highres[where(Lnu_agn_highres lt 0, /null)] = 0.d0
        Lnu_unred_AGN_highres[where(Lnu_unred_AGN_highres lt 0, /null)] = 0.d0
      end

   'NONE': begin
        Lnu_agn_highres = dblarr(Nwave, Nmodels)
        Lnu_unred_AGN_highres = dblarr(Nwave, Nmodels)
        Lbol_abs_AGN = dblarr(Nmodels)
        Lbol_AGN_model = dblarr(Nmodels)
        L2500 = dblarr(Nmodels)
      end
 endcase



; Generate dust emission
 case strupcase(dust_model) of
   'DL07': begin
        alpha = parameters_transposed[*, where(parameter_name eq 'ALPHA', /null)]
        umin  = parameters_transposed[*, where(parameter_name eq 'UMIN', /null)]
        umax  = parameters_transposed[*, where(parameter_name eq 'UMAX', /null)]
        gam   = parameters_transposed[*, where(parameter_name eq 'GAMMA', /null)]
        qPAH  = parameters_transposed[*, where(parameter_name eq 'QPAH', /null)]

        Lnu_spec_dust = dl07_spectrum(models.dust_models, alpha=alpha, umin=umin, umax=umax, gam=gam, $
                                      qPAH=qPAH, LTIR=LTIR_model)

        ; Scale dust emission
        if energy_balance then begin
          LTIR = Lbol_abs_stellar + Lbol_abs_AGN
        endif else begin
          LTIR = parameters_transposed[*, where(parameter_name eq 'LTIR', /null)]
        endelse
        Lnu_dust = rebin(reform(LTIR/LTIR_model, 1, Nmodels), n_elements(models.dust_models.wave_obs), Nmodels) * Lnu_spec_dust

        ; Interpolate dust spectrum to common wavelength grid
        Lnu_dust_highres = dblarr(Nwave, Nmodels)
        for i=0, Nmodels-1 do begin
            Lnu_dust_highres[*, i] = interpol(Lnu_dust[*, i], models.dust_models.WAVE_OBS, wave)
        endfor
        Lnu_dust_highres[where(Lnu_dust_highres lt 0, /null)] = 0.d0
      end

   'NONE': begin
        Lnu_dust_highres = dblarr(Nwave, Nmodels)
        LTIR = dblarr(Nmodels)
      end
 endcase



 return, Lnu_stellar_highres + Lnu_agn_highres + Lnu_dust_highres

end
