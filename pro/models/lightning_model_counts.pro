function lightning_model_counts, parameters, parameter_name, models, Lbol_AGN_model=Lbol_AGN_model, $
                                 L2500=L2500, agn_model=agn_model, xray_agn_model=xray_agn_model, $
                                 error_check=error_check, counts_xrb=counts_xrb, counts_agn=counts_agn
                                 
;+
; Name
; ----
;   LIGHTNING_MODEL_COUNTS
;
; Purpose
; -------
;   Generates the X-ray detector counts (folded through the instrumental 
;   response) of a lightning model under a specified set(s) of parameters.
;
; Calling Sequence
; ----------------
;   ::
;
;       counts = lightning_model_counts(parameters, parameter_name, models [, Lbol_AGN_model = , $
;                                       L2500 = , agn_model = , xray_agn_model = , $
;                                       /error_check, counts_xrb=counts_xrb, counts_agn=counts_agn])
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
;
; Output
; ------
;   ``counts`` : double array(Nxray, Nmodels)
;       The total instrumental counts produced by the X-ray model under the set
;       observing conditions :math:`[{\rm counts}]`.
;
; Optional Outputs
; ----------------
;   ``counts_xrb`` : double array(Nxray, Nmodels)
;       The instrumental counts produced by the X-ray binary model :math:`[{\rm counts}]`.
;   ``counts_agn`` : double array(Nxray, Nmodels)
;       The instrumental counts produced by the X-ray AGN model :math:`[{\rm counts}]`.
;
; Modification History
; --------------------
;   - 2022/06/07: Removed AGN covering factor from model (Erik B. Monson)
;   - 2022/06/20: Major update to include new implementation (e.g., prior, config, etc.) (Keith Doore)
;   - 2022/06/20: Replaced ``!cv`` with ``!lightning_cgs`` (Keith Doore)
;   - 2022/06/20: Updated documentation (Keith Doore)
;   - 2022/06/20: Added error handling (Keith Doore)
;   - 2022/06/20: Added ``error_check`` keyword to do error handling (Keith Doore)
;   - 2022/07/06: Removed ``config`` and replaced ``config`` tag calls with inputs (Keith Doore)
;   - 2022/07/06: Transposed parameters array to eliminate need to reform it after indexing. (Keith Doore)
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
 endif

 if n_elements(agn_model) eq 0 then agn_model = 'NONE'
 if n_elements(xray_agn_model) eq 0 then xray_agn_model = 'QSOSED'


; Generate model counts
 redshift = models.xray_models.REDSHIFT
 Nxray = n_elements(models.xray_models.XRAY_EXPOSURE)
 Nwave = n_elements(models.xray_models.WAVE_OBS)
 if size(parameters, /n_dim) le 1 then Nmodels = 1 else Nmodels = (size(parameters, /dim))[1]

 xray_bandpass_nu = models.xray_models.XRAY_BANDPASS / (!lightning_cgs.hplanck / !lightning_cgs.keV)
 nu_obs = !lightning_cgs.clight * 1e4 / models.xray_models.WAVE_OBS


 ; Transpose parameters array to eliminate need to reform to remove padded dimension after indexing.
 parameters_transposed = transpose(parameters)


 psi = transpose(parameters_transposed[*, where(strmatch(parameter_name, 'PSI_*'), /null)])
 Lnu_XRB = gilbertson22_LX_tau(models.stellar_models, psi)
 Lnu_LMXB = Lnu_XRB[0,*]
 Lnu_HMXB = Lnu_XRB[1,*]

 ; In Woody's parameterization, the hot gas is basically included in the HMXBs?
 Lnu_LMXB = rebin(reform(Lnu_LMXB, 1, Nmodels), Nxray, Nmodels)
 Lnu_HMXB = rebin(reform(Lnu_HMXB, 1, Nmodels), Nxray, Nmodels)
 ; Lnu_hotgas = rebin(reform(Lnu_hotgas, 1, Nmodels), Nxray, Nparallel)

 counts_xrb_norm = dblarr(Nxray, Nmodels)
 ; counts_hotgas = dblarr(Nxray, Nmodels)
 nH = parameters_transposed[*, where(parameter_name eq 'NH', /null)]
 nH = rebin(reform(nH, 1, Nmodels), n_elements(models.xray_models.WAVE_OBS), Nmodels)
 ; If we're using the agn model, we'll treat the stellar column density as a constant scaled from the
 ; galaxy's AV. Otherwise, the nH parameter will control only the absorption of the stellar population.
 case strupcase(agn_model) of
   'SKIRTOR': begin
           ; tauV can be either tauV from Calzetti model or tauV_diff from modified Calzetti
           tauv = parameters_transposed[*, where(parameter_name eq 'TAUV' or $
                                      parameter_name eq 'TAUV_DIFF', /null)]
           nH_stellar = 22.4 * (-2.5 * alog10(exp(-1 * tauV)))
           nH_stellar = rebin(reform(nH_stellar, 1, Nmodels), n_elements(models.xray_models.WAVE_OBS), Nmodels)
      end
   'NONE': nH_stellar = nH
 endcase

 for i=0, Nxray - 1 do begin
   for j=0, Nmodels - 1 do begin
     counts_xrb_norm[i, j] = models.xray_models.XRAY_EXPOSURE[i] * trap_int(reverse(nu_obs), $
                             reverse(models.xray_models.EXP_NEG_TAU_XRAY^(nH_stellar[*, j]) * models.xray_models.XRB_MODEL), $
                             XRANGE=xray_bandpass_nu[*, i])
     ; counts_hotgas[i, j] = models.xray_models.XRAY_EXPOSURE[i] * trap_int(reverse(nu_X_obs), $
     ;                              reverse(models.xray_models.EXP_NEG_TAU_XRAY^(nH[*, j]) * models.xray_models.GAS_MODEL), $
     ;                              XRANGE=xray_box_nu[*, i])
   endfor
 endfor

 ; Note that since we integrate against nu_obs, the model counts_* have dimensions [counts / Lnu_obs].
 ; Since the Lnu_LMXB etc. are Lnu_rest, we multiply by (1 + z) to get the proper number of counts.
 counts_xrb = (1 + redshift) * (Lnu_LMXB * counts_xrb_norm + Lnu_HMXB * counts_xrb_norm)


 case strupcase(agn_model) of
   'SKIRTOR': begin
      counts_AGN_norm = dblarr(Nxray, Nmodels)

      case strupcase(xray_agn_model) of
        'PLAW': begin
            ; The L_2500A we get from the template should be normalized by the bolometric
            ; luminosity of the template and then multiplied by the current AGN luminosity parameter.
            ; The L_2500A is also an observed-frame Lnu, like all the Lnu that come out of our models, so we convert it
            ; to rest-frame before plugging it into the LR17 relation.
            ; Here, this means monochromatic 2 keV luminosity in Lnu_obs units.
            ; Redshift does not cancel, since L2keV_LR17() has multiplicative factor in front of the log
            log_l_agn = parameters_transposed[*, where(parameter_name eq 'LOG_L_AGN', /null)]
            agn_xray_norm = L2keV_LR17((10.d0 ^ log_l_agn) * L2500 / Lbol_AGN_model / (1 + redshift))
            agn_xray_norm = (1 + redshift) * rebin(reform(agn_xray_norm, 1, Nmodels), Nxray, Nmodels)
            agn_xray_model = models.xray_models.AGN_MODEL
            agn_xray_model = rebin(agn_xray_model, Nwave, Nmodels)
          end

        'QSOSED': begin
            ; In the case of the qsosed model we normalize the counts directly by L2500, so we don't need the LR17
            ; relationship. Note that we also don't multiply by (1 + z) because L_2500A is observed frame Lnu as returned
            ; by the AGN model.
            ; L_2500A_norm = (10.d0 ^ log_l_agn) * agn_L_2500A / agn_L_bol
            ; L_2500A_norm = rebin(reform(L_2500A_norm, 1, Nparallel), Nxray, Nparallel)
            ; counts_new = counts_new + (L_2500A_norm * counts_AGN)
   
            ; ; Here, this means monochromatic 2500 A luminosity in Lnu_obs units
            ; agn_xray_norm = (10.d0 ^ log_l_agn) * L2500 / Lbol_AGN_model
            ; agn_xray_norm = rebin(reform(agn_xray_norm, 1, Nmodels), Nxray, Nmodels)
            agn_mass = parameters_transposed[*, where(parameter_name eq 'AGN_MASS', /null)]
            agn_logmdot = parameters_transposed[*, where(parameter_name eq 'AGN_LOGMDOT', /null)]
            agn_xray_norm = replicate(1, Nxray, Nmodels)
            agn_xray_model = qsosed_spectrum(models.xray_models, agn_mass, agn_logmdot)
          end

        'NONE': begin
            agn_xray_norm = dblarr(Nxray, Nmodels)
            agn_xray_model = dblarr(Nwave, Nmodels)
          end
      endcase

      for i=0, Nxray - 1 do begin
        for j=0, Nmodels - 1 do begin
          counts_AGN_norm[i, j] = models.xray_models.XRAY_EXPOSURE[i] * trap_int(reverse(nu_obs), $
                                  reverse(models.xray_models.EXP_NEG_TAU_XRAY^(nH[*, j]) * $
                                  agn_xray_model[*, j]), XRANGE=xray_bandpass_nu[*, i])
        endfor
      endfor

      ; X-ray anisotropy is implemented as L_X = cosi * L_X(0)
      ; if (not keyword_set(agn_isotropy)) then begin
      ;     agn_xray_norm = rebin(reform(parameters[13,*], 1, Nmodels), Nxray, Nmodels) * agn_xray_norm
      ; endif else begin
      ;     ; cos(60 degrees) = 0.5
      ;     agn_xray_norm = 0.5 * agn_xray_norm
      ; endelse
   
      ; TODO: figure out how to do AGN X-ray anisotropy with the new qsosed model. It's backwards now?
      ;       who knows.

      counts_agn = agn_xray_norm * counts_AGN_norm
    end   

   'NONE': counts_agn = dblarr(Nxray, Nmodels)
 endcase


 return, counts_xrb + counts_agn

end