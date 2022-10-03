function skirtor_sed, skirtor, tau97=tau97, i_agn=i_agn, tauV_diff=tauV_diff, delta=delta, $
                      atten_curve=atten_curve, uv_bump=uv_bump, error_check=error_check, $
                      xray_isotropic=xray_isotropic, mean_Lnu_unred_AGN=mean_Lnu_unred_AGN, $
                      Lbol_abs_AGN=Lbol_abs_AGN, Lbol_AGN=Lbol_AGN, L2500=L2500
;+
; Name
; ----
;   SKIRTOR_SED
;
; Purpose
; -------
;   Generates the AGN model SED for a given set of AGN model parameters and
;   filters. The AGN model is based on the SKIRTOR models (Stalevski et al. 2016),
;   which include the parameters of :math:`\tau_{9.7}` and :math:`i_{\rm AGN}`.
;
; Calling Sequence
; ----------------
;   ::
;
;       mean_Lnu_AGN = skirtor_sed(skirtor [, tau97 = , i_agn = , tauV_diff = , delta = , $
;                                  atten_curve = , /uv_bump, /error_check, /xray_isotropic, $
;                                  mean_Lnu_unred_AGN=mean_Lnu_unred_AGN, Lbol_abs_AGN=Lbol_abs_AGN, $
;                                  Lbol_AGN=Lbol_AGN, L2500=L2500])
;
; Input
; -----
;   ``skirtor`` : structure
;       A structure containing the spectra, SEDs, and AGN parameters for the
;       SKIRTOR model. (See ``skirtor_models.pro`` for details and contents.)
;
; Optional Inputs
; ---------------
;   ``tau97`` : int, float, or double array(Nmodels)
;       Edge-on optical depth of the AGN torus at 9.7 um. (Default = ``3.0``)
;   ``i_agn`` : int, float, or double array(Nmodels)
;       Inclination from the polar axis to the line of sight :math:`[{\rm degrees}]`.
;       (Default = ``0.0``)
;   ``tauV_diff`` : int, float, or double array(Nmodel)
;       The V-band optical depth of the diffuse dust. If this input is set,
;       the model will be attenuated.
;   ``delta`` : int, float, or double array(Nmodels)
;       The power law value to change the attenuation curve slope.
;   ``atten_curve`` : string scalar
;       The name of the attenuation curve to apply to the AGN models. Current
;       options are ``'CALZETTI00'`` and ``'CALZETTI_MOD'``.
;   ``uv_bump`` : flag
;       If set, then a 2175 Angstrom UV bump feature will be added to the
;       attenuation curve.
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;   ``xray_isotropic`` : flag
;       Deprecated (?). If set, the output ``L2500`` is assumed to be isotropic from the
;       face-on model ``L2500``. Otherwise, it is assumed to be a function of viewing angle.
;
; Output
; ------
;   ``mean_Lnu_AGN`` : double array(Nfilters, Nmodels)
;       The mean luminosity of each filter and set of AGN model parameters
;       :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;
; Optional Outputs
; ----------------
;   ``mean_Lnu_unred_AGN`` : double array(Nfilters, Nmodels)
;       The unattenuated mean luminosity of each filter and set of AGN model
;       parameters :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   ``Lbol_abs_AGN``:  double array(Nmodels)
;       The bolometric luminosity of the attenuated AGN emission for each
;       set of model parameters, integrated for inclination variation :math:`[L_\odot]`.
;   ``Lbol_AGN`` : double array(Nmodels)
;       The bolometric luminosity of the AGN model for each set of model parameters,
;       integrated for inclination variation :math:`[L_\odot]`.
;   ``L2500`` : double array(Nmodels)
;       The rest-frame 2500 Angstrom monochromatic luminosity of the intrinsic accretion
;       disk spectrum shifted to the observed frame :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;
; Note
; ----
;   Defaults will only be set if the optional ``error_check`` input is set.
;
; Reference
; ---------
;   `Stalevski, M., Ricci, C., Ueda, Y., et al. 2016, MNRAS, 458, 2288 <https://ui.adsabs.harvard.edu/abs/2016MNRAS.458.2288S/abstract>`_
;
; Modification History
; --------------------
;   - 2021/03/03: Created (Erik B. Monson)
;   - 2021/03/19: Added optional output for rest-frame ``L2500`` (Erik B. Monson)
;   - 2021/09/14: Interpolation now proceeds in log(Lnu) (Erik B. Monson)
;   - 2022/03/15: Documentation update (Erik B. Monson)
;   - 2022/06/09: Changed parameter names to standardized format and updated documentation (Keith Doore)
;   - 2022/06/09: Added error handling (Keith Doore)
;   - 2022/06/09: Added ``error_check`` keyword to do error handling (Keith Doore)
;   - 2022/06/15: Replaced ``!cv`` with ``!lightning_cgs`` (Keith Doore)
;   - 2022/06/30: Updated documentation and renamed variables (Keith Doore)
;   - 2022/06/30: Removed ``config`` and replaced with the two accessed configuration values (Keith Doore)
;   - 2022/06/30: Fixed potential bug of ``atten_curve`` not being set (Keith Doore)
;   - 2022/07/22: Added optional output of unreddened mean Lnu (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(skirtor) eq 0 then message, 'Variable is undefined: SKIRTOR.'
   if size(skirtor, /type) ne 8 then message, 'SKIRTOR must be of type structure.'

   if n_elements(tau97) ne 0 then begin
     if size(tau97, /type) lt 2 or size(tau97,/type) gt 5 then $
              message, 'TAU97 must be of type int, float, or double.'
     if size(reform(tau97), /n_dim) ne 1 then message, 'TAU97 must be a scalar or 1-D array.'
     if min(tau97) lt 3 or max(tau97) gt 11 then message, 'TAU97 must contain values between 3 and 11.'
     Nmodels = n_elements(tau97)
   endif

   if n_elements(i_agn) ne 0 then begin
     if size(i_agn, /type) lt 2 or size(i_agn,/type) gt 5 then $
              message, 'I_AGN must be of type int, float, or double.'
     if size(reform(i_agn), /n_dim) ne 1 then message, 'I_AGN must be a scalar or 1-D array.'
     if min(i_agn) lt 0 or max(i_agn) gt 90 then message,'I_AGN must contain values between 0 and 90.'
     if n_elements(Nmodels) eq 0 then Nmodels = n_elements(i_agn)
   endif

   if n_elements(Nmodels) ne 0 then begin
     if n_elements(tau97) eq 0 then tau97 = replicate(3.d, Nmodels)
     if n_elements(i_agn) eq 0 then i_agn = replicate(0.d, Nmodels)
   endif else begin
     Nmodels = 1
     tau97 = 3.d
     i_agn = 0.d
   endelse

   if n_elements(tau97) ne n_elements(i_agn) then $
      message,'TAU97 and I_AGN must have the same size.'


   if n_elements(tauV_diff) ne 0 then begin
     if size(tauV_diff, /type) lt 2 or size(tauV_diff,/type) gt 5 then $
              message, 'TAUV_DIFF must be of type int, float, or double.'
     if size(reform(tauV_diff), /n_dim) ne 1 then message, 'TAUV_DIFF must be a scalar or 1-D array.'
     if min(tauV_diff) lt 0 then message,'TAUV_DIFF must only contain non-negative values.'
     if n_elements(tauV_diff) ne Nmodels then message,'TAUV_DIFF must have the same size as TAU97 and I_AGN.'
   endif

   if n_elements(delta) ne 0 then begin
     if size(delta, /type) lt 2 or size(delta,/type) gt 5 then $
              message, 'DELTA must be of type int, float, or double.'
     if size(reform(delta), /n_dim) ne 1 then message, 'DELTA must be a scalar or 1-D array.'
     if n_elements(delta) ne Nmodels then message,'DELTA must have the same size as TAU97 and I_AGN.'
   endif

   if n_elements(atten_curve) ne 0 then begin
     if size(atten_curve, /type) ne 7 then $
       message, 'ATTEN_CURVE must be of type string.'
     if size(atten_curve, /n_dim) ne 0 then message, 'ATTEN_CURVE must be a scalar.'
     if total(strupcase(atten_curve) eq ['CALZETTI00', 'CALZETTI_MOD']) ne 1 then $
       message, "ATTEN_CURVE must be set to either 'CALZETTI00' or 'CALZETTI_MOD'."
   endif
 endif

 if n_elements(atten_curve) eq 0 and n_elements(tauV_diff) ne 0 then begin
   if n_elements(delta) ne 0 then atten_curve = 'CALZETTI_MOD' else atten_curve = 'CALZETTI00'
 endif

; Generate model Lnu
 Nmodels = n_elements(tau97)

 tau_grid = [3,5,7,9,11]
 N_tau_grid = n_elements(tau_grid)
 tau_grid_min = min(tau_grid)
 tau_grid_max = max(tau_grid)

 inc_grid = [0,10,20,30,40,50,60,70,80,90]
 cosi_grid = cos(!DPI * inc_grid / 180.d0)
 N_inc_grid = n_elements(inc_grid)
 inc_idcs = indgen(N_inc_grid)
 ;inc_grid_min = min(inc_grid)
 ;inc_grid_max = max(inc_grid)


 Nfilters = n_elements(skirtor.FILTER_LABELS)
 mean_Lnu_AGN = dblarr(Nfilters, Nmodels)
 Lbol_abs_AGN = dblarr(Nmodels)
 Lbol_AGN = dblarr(Nmodels)
 L2500 = dblarr(Nmodels)


 nu_obs = (1.d4 * !lightning_cgs.clight) / skirtor.WAVE_OBS

 ; Normalize array coordinates for input inclination and optical depth

 ; This is only correct for an evenly spaced grid.
 tau_norm = (n_tau_grid - 1) * (tau97 - tau_grid_min) / (tau_grid_max - tau_grid_min)

 ; To get the correct indices for cosi (which is super not evenly spaced)
 ; we need to do another (cheap) interpolation beforehand.
 cosi_agn = cos(!DPI * i_agn / 180.d0)
 inc_norm = interpol(inc_idcs, cosi_grid, cosi_agn)

 Lnu_spec = 10.d0 ^ interpolate(alog10(skirtor.LNU_TOTAL), inc_norm, tau_norm)

 if n_elements(tauV_DIFF) ne 0 then begin
   case strupcase(atten_curve) of
     'CALZETTI00': begin
          exp_neg_tau    = calzetti00_atten(skirtor.WAVE_REST,   $
                                            tauV_DIFF=tauV_DIFF,    $
                                            delta=dblarr(Nmodels),  $
                                            tauV_BC=dblarr(Nmodels))
        end

     'CALZETTI_MOD': begin
          exp_neg_tau    = calzetti00_atten(skirtor.WAVE_REST,    $
                                            tauV_DIFF=tauV_DIFF,     $
                                            delta=delta,             $
                                            tauV_BC=dblarr(Nmodels), $
                                            uv_bump=uv_bump)
        end
   endcase

   ; Much of the emission from the AGN accretion disk comes out at UV energies higher than
   ; the lyman limit. We should make sure the ISM is opaque to this. Although this is ionizing
   ; flux, so really it should get put into emission lines instead of dust heating, but we don't
   ; have a great way of doing that. I'll investigate that later. For now, this.
   ;LyC_idcs = where(rebin(skirtor.WAVE_REST, n_elements(skirtor.WAVE_REST), Nmodels) lt 0.0912, /null)
   LyC_idcs = where(skirtor.WAVE_REST lt 0.0912, /null)
   exp_neg_tau[LyC_idcs, *] = 0.d0
 endif else begin
   exp_neg_tau = replicate(1.d0, n_elements(skirtor.WAVE_REST), Nmodels)
 endelse

 L2500 = replicate(skirtor.L2500A_ZERO, Nmodels)

 ; Interpolate the cosi-integrated spectrum so that we can get the bolometric luminosity
 ; (which should be the same between all the models, but isn't quite) and the
 ; attenuated bolometric luminosity
 Lnu_spec_inc_integrated = 10.d0 ^ interpolate(alog10(skirtor.LNU_INTEGRATED), inc_norm, tau_norm)

 ; Generate unattenuated emission
 mean_Lnu_unred_AGN = 10.d0 ^ interpolate(alog10(skirtor.MEAN_LNU), inc_norm, tau_norm)

 ; Loop
 for nm=0, Nmodels - 1 do begin

   Lbol_AGN[nm] = abs(trap_int(nu_obs, Lnu_spec_inc_integrated[*, nm]))
   if n_elements(tauV_DIFF) ne 0 then $
     Lbol_abs_AGN[nm] = abs(trap_int(reverse(nu_obs), reverse((1 - exp_neg_tau[*, nm]) * Lnu_spec_inc_integrated[*, nm]), $
                              XRANGE=reverse([min(nu_obs[LyC_idcs]), min(nu_obs)]))) else Lbol_abs_AGN[nm] = 0.d0

   for f=0, Nfilters-1 do begin
     ; Transmission is negative due to nu_obs being reversed (big to small)
     transmission = reform(skirtor.FILTERS[f, *]) / trapez(reform(skirtor.FILTERS[f, *]), nu_obs)
     ; mean_Lnu_AGN is positive due to nu_obs being reversed (big to small) and transmission being negative
     mean_Lnu_AGN[f, nm] = trapez(transmission * exp_neg_tau[*, nm] * Lnu_spec[*, nm], nu_obs)
   endfor

 endfor

 return, mean_Lnu_AGN

end
