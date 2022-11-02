function xrb_xagn_models, xray_bandpass, xray_exposure=xray_exposure,$
                          arf_E_lo=arf_E_lo, arf_E_hi=arf_E_hi, arf_specresp=arf_specresp, $
                          redshift=redshift, lumin_dist=lumin_dist, galactic_nH=galactic_nH, $
                          xray_abs_model=xray_abs_model, xray_agn_model=xray_agn_model, xray_unit=xray_unit,$
                          error_check=error_check
;+
; Name
; ----
;   XRB_XAGN_MODELS
;
; Purpose
; -------
;   Generates the spectra, counts, and X-ray parameters for given a set of
;   X-ray bandpasses and ARF for an X-ray binary population (XRB). The
;   spectra and counts can include or not include Galactic absorption
;   and X-ray AGN emission.
;
; Calling Sequence
; ----------------
;   ::
;
;       xray_models = xrb_xagn_models(xray_bandpass [, xray_exposure = , arf_E_lo = , arf_E_hi = , $
;                                     arf_specresp = , redshift = , lumin_dist = , $
;                                     galactic_nH = , xray_abs_model = , $
;                                     xray_agn_model = , xray_unit = , /error_check])
;
; Input
; -----
;   ``xray_bandpass`` : int, float or double array(2, Nxray)
;       The bandpasses for the X-ray spectrum. The first column should be the lower
;       energy bound, and the second should be the upper energy bound :math:`[{\rm keV}]`.
;
; Optional Inputs
; ---------------
;   ``xray_exposure`` : int, float or double array(Nxray)
;       The exposure time of the observations, one per band :math:`[{\rm s}]`.
;       Required to generate model count-rate spectrum.
;   ``arf_E_lo`` : float or double array(Nchannels)
;       Lower energy bounds of each channel in the ARF :math:`[{\rm keV}]`.
;       Required to generate model count-rate spectrum.
;   ``arf_E_hi`` : float or double array(Nchannels)
;       Upper energy bounds of each channel in the ARF :math:`[{\rm keV}]`.
;       Required to generate model count-rate spectrum.
;   ``arf_specresp`` : float or double array(Nchannels)
;       The spectral response of the ARF at each channel :math:`[{\rm cm}^2]`.
;       Required to generate model count-rate spectrum.
;   ``redshift`` : int, float, or double scalar
;       The redshift of the model. (Default = ``0.0``)
;   ``lumin_dist`` : int, float, double scalar
;       The luminosity distance of the model :math:`[{\rm Mpc}]`. (Default = ``10``)
;   ``galactic_nH`` : int, float, or double scalar
;       Galactic, i.e. Milky Way, neutral Hydrogen column density along the line of sight
;       :math:`[10^{20}\ {\rm cm}^{-2}]`. (Default = ``0``)
;   ``xray_abs_model`` : string scalar
;       The name of the X-ray absorption model to apply to the X-ray emission. Current options
;       are: ``'TBABS-WILM'``, ``'TBABS-ANGR'``, ``'ATTEN'``, and ``'NONE'``.
;       (Default = ``'TBABS-WILM'``)
;   ``xray_agn_model`` : string scalar
;       The X-ray AGN emission model to use. Current options are: ``'PLAW'``, ``'QSOSED'``, ``'NONE'``.
;       (Default = ``'QSOSED'``)
;   ``xray_unit`` : string scalar
;       The type of X-ray data to use in fitting. Current options are ``'COUNTS'`` and ``'FLUX'``.
;       (Default = ``'COUNTS'``)
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``xray_models`` : structure
;       A structure including the the XRB, hot gas, and X-ray AGN models in terms
;       of luminosities and count rates given the input ARF.
;       The full description of the structure is as follows:
;
;       ======================     ===========================     ================================================================================================================================================
;       TAG                        TYPE                            DESCRIPTION
;       ======================     ===========================     ================================================================================================================================================
;       XRAY_BANDPASS              double(2, Nxray)                X-ray bandpasses :math:`[{\rm keV}]`, same as input
;       XRAY_EXPOSURE [1]_         double(Nxray)                   X-ray exposure time per band :math:`[{\rm s}]`, same as input
;       WAVE_REST                  double(Nwave)                   Rest-frame wavelength of the models :math:`[\mu \rm m]`
;       WAVE_OBS                   double(Nwave)                   Observed-frame wavelength of the models :math:`[\mu \rm m]`
;       EXP_NEG_TAU_XRAY           double(Nwave)                   Rest-frame absorption curve normalized to :math:`{\rm nH} = 10^{20}\ {\rm cm}^2`
;       EXP_NEG_TAU_XRAY_MW        double(Nwave)                   Observed-frame absorption curve normalized to :math:`{\rm nH} = 10^{20}\ {\rm cm}^2`
;       LNU_AGN [2]_               double(Nwave, Nmass, Nmdot)     AGN model normalized to 1 :math:`L_\odot\ {\rm Hz}^{-1}` at rest-frame 2 keV :math:`[L_\odot\ {\rm Hz}^{-1}]`
;       LNU_XRB                    double(Nwave)                   XRB model normalized to 1 :math:`L_\odot\ {\rm Hz}^{-1}` over rest-frame 2-10 keV bandpass :math:`[L_\odot\ {\rm Hz}^{-1}]`
;       LNU_GAS                    double(Nwave)                   Hot gas model normalized to 1 :math:`L_\odot\ {\rm Hz}^{-1}` over rest-frame 0.5-2 keV bandpass :math:`[L_\odot\ {\rm Hz}^{-1}]`
;       AGN_MODEL [1]_ [2]_        double(Nwave, Nmass, Nmdot)     Instrumental count rate density produced by normalized AGN model :math:`[{\rm counts\ s^{-1}\ Hz^{-1}}]`
;       XRB_MODEL [1]_             double(Nwave)                   Instrumental count rate density produced by normalized XRB model :math:`[{\rm counts\ s^{-1}\ Hz^{-1}}]`
;       GAS_MODEL [1]_             double(Nwave)                   Instrumental count rate density produced by normalized hot gas model :math:`[{\rm counts\ s^{-1}\ Hz^{-1}}]`
;       XRB_COUNTS_LNU [1]_        double(Nxray)                   Instrumental counts produced by normalized XRB model integrated over the supplied bandpass (Galactic absorption only) :math:`[{\rm counts}]`
;       GAS_COUNTS_LNU [1]_        double(Nxray)                   Instrumental counts produced by normalized hot gas model integrated over the supplied bandpass (Galactic absorption only) :math:`[{\rm counts}]`
;       GALACTIC_NH                double                          Galactic neutral hydrogen column density :math:`[10^{20}\ {\rm cm}^{-2}]`, same as input
;       L2500 [3]_                 double(Nmass, Nmdot)            2500 Angstrom intrinsic luminosity grid from qsosed model, shifted to the observed frame :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;       AGN_MASS [3]_              double(Nmass, Nmdot)            Supermassive black hole mass grid from qsosed model :math:`[M_\odot]`
;       AGN_LOGMDOT [3]_           double(Nmass, Nmdot)            Log10 of SMBH accretion rate grid from qsosed model, normalized by the Eddington rate
;       REDSHIFT                   double                          Redshift of the model, same as input
;       XRAY_ABS_MODEL             string                          Name of the X-ray absorption model, same as input
;       XRAY_AGN_MODEL             string                          Name of the AGN model, same as input
;       ======================     ===========================     ================================================================================================================================================
;
;   .. [1] ``XRAY_EXPOSURE``, ``AGN_MODEL``, ``XRB_MODEL``, ``GAS_MODEL``, ``XRB_COUNTS_LNU``, and ``GAS_COUNTS_LNU``
;      are set to ``NaN`` if the ``xray_exposure`` and ``arf_*`` keywords are not specified.
;   .. [2] ``LNU_AGN`` and ``AGN_MODEL`` will only have the second and third dimensions if
;      using the qsosed model. Additionally, the qsosed model is not normalized to 1 :math:`L_\odot\ {\rm Hz}^{-1}`.
;   .. [3] ``L2500``, ``AGN_MASS``, and ``AGN_LOGMDOT`` only have values if using the qsosed
;      model. Otherwise, they will be ``NaN``.
;
; Modification History
; --------------------
;   - 2021/09/21: Created (Erik B. Monson).
;   - 2021/10/11: Added X-ray absorption (Erik B. Monson).
;   - 2022/03/16: Moved to separate file and documentation improved (Erik B. Monson).
;   - 2022/04/18: qsosed option now loads entire grid (Erik B. Monson)
;   - 2022/06/22: Major update to include new implementation (Keith Doore)
;   - 2022/07/07: Changed several variable names and updated documentation (Keith Doore)
;   - 2022/09/14: Updates to allow fitting with X-ray fluxes (Erik B. Monson)
;   - 2022/11/02: Galactic NH is now in units of 1e20 cm-2 (Erik B. Monson)
;-
 On_error, 2
 compile_opt idl2

 ; Error handling
 if keyword_set(error_check) then begin
   if n_elements(xray_bandpass) eq 0 then message, 'Variable is undefined: XRAY_BANDPASS.'
   if size(xray_bandpass, /type) lt 2 or size(xray_bandpass, /type) gt 5 then $
     message, 'XRAY_BANDPASS must be of type int, float, or double.'
   if size(reform(xray_bandpass), /n_dim) lt 1 or size(reform(xray_bandpass), /n_dim) gt 2 then $
     message, 'XRAY_BANDPASS must be a 1-D or 2-D array.'
   if min(xray_bandpass) le 0 then message, 'XRAY_BANDPASS must only contain positive values.'
   if size(xray_bandpass, /n_dim) eq 1 then Nxray = 1 else Nxray = (size(xray_bandpass, /dim))[1]

   if n_elements(xray_unit) ne 0 then begin
     if size(xray_unit, /type) ne 7 then message, 'XRAY_UNIT must be of type string.'
     if size(xray_unit, /dim) ne 0 then message, 'XRAY_UNIT must be a scalar.'
     if total(strupcase(xray_unit) eq ['COUNTS', 'FLUX']) ne 1 then $
       message, "XRAY_UNIT must be set to either 'COUNTS' or 'FLUX'."
   endif else xray_unit = 'COUNTS'

   ; The following inputs are only necessary if we need the model to predict
   ;    the number of observed counts.
   if strupcase(xray_unit) eq 'COUNTS' then begin

       if n_elements(xray_exposure) eq 0 then message, 'Variable is undefined: XRAY_EXPOSURE.'
       if size(xray_exposure, /type) lt 2 or size(xray_exposure, /type) gt 5 then $
         message, 'XRAY_EXPOSURE must be of type int, float, or double.'
       if size(reform(xray_exposure), /n_dim) ne 1 then $
         message, 'XRAY_EXPOSURE must be a 1-D array.'
       if min(xray_exposure) le 0 then message, 'XRAY_EXPOSURE must only contain positive values.'
       if n_elements(xray_exposure) ne Nxray then $
         message, 'XRAY_EXPOSURE must have the same number of elements as the second dimension of XRAY_BANDPASS'

       if n_elements(arf_E_lo) eq 0 then message, 'Variable is undefined: ARF_E_LO.'
       if size(arf_E_lo, /type) lt 2 or size(arf_E_lo, /type) gt 5 then $
         message, 'ARF_E_LO must be of type int, float, or double.'
       if size(reform(arf_E_lo), /n_dim) ne 1 then $
         message, 'ARF_E_LO must be a 1-D array.'
       if min(arf_E_lo) le 0 then message, 'ARF_E_LO must only contain positive values.'

       if n_elements(arf_E_hi) eq 0 then message, 'Variable is undefined: ARF_E_HI.'
       if size(arf_E_hi, /type) lt 2 or size(arf_E_hi, /type) gt 5 then $
         message, 'ARF_E_HI must be of type int, float, or double.'
       if size(reform(arf_E_hi), /n_dim) ne 1 then $
         message, 'ARF_E_HI must be a 1-D array.'
       if min(arf_E_hi) le 0 then message, 'ARF_E_HI must only contain positive values.'
       if n_elements(arf_E_hi) ne n_elements(arf_E_lo) then $
         message, 'ARF_E_HI must have the same number of elements as ARF_E_LO.'

       if n_elements(arf_specresp) eq 0 then message, 'Variable is undefined: ARF_SPECRESP.'
       if size(arf_specresp, /type) lt 2 or size(arf_specresp, /type) gt 5 then $
         message, 'ARF_SPECRESP must be of type int, float, or double.'
       if size(reform(arf_specresp), /n_dim) ne 1 then $
         message, 'ARF_SPECRESP must be a 1-D array.'
       if min(arf_specresp) le 0 then message, 'ARF_SPECRESP must only contain positive values.'
       if n_elements(arf_specresp) ne n_elements(arf_E_lo) then $
         message, 'ARF_SPECRESP must have the same number of elements as ARF_E_LO.'
   endif

   if n_elements(redshift) ne 0 then begin
     if size(redshift, /type) lt 2 or size(redshift,/type) gt 5 then $
              message, 'REDSHIFT must be of type int, float, or double.'
     if size(redshift, /n_dim) ne 0 then message, 'REDSHIFT must be a scalar.'
     if redshift lt 0 then message, 'REDSHIFT must be a non-negative value.'
   endif

   if n_elements(lumin_dist) ne 0 then begin
     if size(lumin_dist, /type) lt 2 or size(lumin_dist,/type) gt 5 then $
              message, 'LUMIN_DIST must be of type int, float, or double.'
     if size(lumin_dist, /n_dim) ne 0 then message, 'LUMIN_DIST must be a scalar.'
     if lumin_dist le 0 then message, 'LUMIN_DIST must be a positive value.'
   endif

   if n_elements(galactic_nH) ne 0 then begin
     if size(galactic_nH, /type) lt 2 or size(galactic_nH,/type) gt 5 then $
              message, 'GALACTIC_NH must be of type int, float, or double.'
     if size(galactic_nH, /n_dim) ne 0 then message, 'GALACTIC_NH must be a scalar.'
     if galactic_nH lt 0 then message, 'GALACTIC_NH must be a non-negative value.'
   endif

   if n_elements(xray_abs_model) ne 0 then begin
     if size(xray_abs_model, /type) ne 7 then message, 'XRAY_ABS_MODEL must be of type string.'
     if size(xray_abs_model, /n_dim) ne 0 then message, 'XRAY_ABS_MODEL must be a scalar.'
     if total(strupcase(xray_abs_model) eq ['TBABS-WILM', 'TBABS-ANGR', 'ATTEN', 'NONE']) ne 1 then $
       message, "XRAY_ABS_MODEL must be set to either 'TBABS-WILM', 'TBABS-ANGR', 'ATTEN', or 'NONE'."
   endif

   if n_elements(xray_agn_model) ne 0 then begin
     if size(xray_agn_model, /type) ne 7 then message, 'XRAY_AGN_MODEL must be of type string.'
     if size(xray_agn_model, /n_dim) ne 0 then message, 'XRAY_AGN_MODEL must be a scalar.'
     if total(strupcase(xray_agn_model) eq ['PLAW', 'QSOSED', 'NONE']) ne 1 then $
       message, "XRAY_AGN_MODEL must be set to either 'PLAW', 'QSOSED', or 'NONE'."
   endif
 endif

 if size(xray_bandpass, /n_dim) eq 1 then Nxray = 1 else Nxray = (size(xray_bandpass, /dim))[1]
 if n_elements(redshift) eq 0 then redshift = 0.d
 if n_elements(lumin_dist) eq 0 then lumin_dist = 10.d
 if n_elements(galactic_nH) eq 0 then galactic_nH = 0.d
 if n_elements(xray_abs_model) eq 0 then xray_abs_model = 'TBABS-WILM'
 if n_elements(xray_agn_model) eq 0 then xray_agn_model = 'QSOSED'



; Compute the x-ray models
 wave_rest = 10.d0^dindgen(185, start=-6, increment=0.02)
 Nwave = n_elements(wave_rest)
 ; Shift the rest wavelength to the observed frame
 wave_obs = (1 + redshift) * wave_rest


 ; Construct the underlying X-ray models.
 ; These give Lnu / Lnu_band, where band is defined by the normbox keyword; for the AGN
 ; we use a narrow band around 2 keV. Lnu_band is in Lsun Hz-1.
 ; We might want two spectral shapes for LMXBs and HMXBs with different photon indices,
 ; but for now we only use one.
 Lnu_x_xrb = xray_plaw_expcut(wave_rest, plaw_gamma=1.8, E_cut=100.0, normbox=[2.d0, 10.d0])
 Lnu_x_hotgas = xray_plaw_expcut(wave_rest, plaw_gamma=1.0, E_cut=1.0, normbox=[0.5d0, 2.d0])


 ; Load absorption curve, if any
 exp_neg_tau_xray = load_xray_abs(wave_rest, xray_abs_model)
 ; ...and shift it to the observed frame, for Galactic absorption.
 exp_neg_tau_xray_MW = interpol(exp_neg_tau_xray, wave_rest, wave_obs)
 ; Interpolation issue at z >~ 1.0: soft X-ray absorption extrapolated to negative values
 negative_idcs = where(exp_neg_tau_xray_MW lt 0.d0, count_negative, /null)
 if (count_negative gt 0) then exp_neg_tau_xray_MW[negative_idcs] = 0.d0

 if strupcase(xray_unit) eq 'COUNTS' then begin
     arf_interp = interpolate_arf(arf_E_lo, arf_E_hi, arf_specresp, wave_obs)
 endif else begin
     arf_interp = !values.D_NAN + dblarr(Nwave)
     xray_exposure = !values.D_NAN + dblarr(Nxray)
 endelse


 ; We need the luminosity distance in cm to convert Lnu to fnu to counts s-1 Hz-1
 DL = lumin_dist
 if (DL lt 10) then DL = 10.d ; We can't really handle very small distances. 10 Mpc is, I think, like z = 0.0023
 DL = DL * 1.d6 * !lightning_cgs.pc
 LtoF = 1.d0 / (4.d0 * !dpi * DL^2)
 ; Planck constant * speed of light is the conversion between wavelength and photon energy.
 hc = (!lightning_cgs.hplanck / !lightning_cgs.keV) * (1e4 * !lightning_cgs.clight) ; h * c in keV * micron
 ; Convert from energy units to photons to counts (recall that the ARF is counts photon-1 cm2)
 ; We want the photon energy in Lsun s for convenience later.
 photon_energy_obs = hc / wave_obs * (!lightning_cgs.keV / !lightning_cgs.Lsun)

 ; 1 / (4 pi DL^2) is a *very* small number in cm-2
 ; as is the energy of an X-ray photon in Lsun s, so we combine them
 ; here to avoid underflowing our models to 0
 LtoF_over_photon_energy_obs = LtoF / photon_energy_obs


 ; Convert the X-ray models to counts s-1 Hz-1 / L_band, where L_band is as above.
 ;  If we've got our X-ray data in fluxes, this will all just come out as NaNs
 model_xrb_count_rate = arf_interp * (exp_neg_tau_xray_MW ^ (galactic_nH) * Lnu_x_xrb * LtoF_over_photon_energy_obs)
 model_hotgas_count_rate = arf_interp * (exp_neg_tau_xray_MW ^ (galactic_nH) * Lnu_x_hotgas * LtoF_over_photon_energy_obs)


 ; Integrate the XRB and hot gas models to get counts / L_band for each bin.
 ; AGN model must be integrated inside the loops.
 xray_bandpass_nu = xray_bandpass / (!lightning_cgs.hplanck / !lightning_cgs.keV)
 nu_obs = (1e4 * !lightning_cgs.clight) / wave_obs

 xrb_counts = !values.D_NAN + dblarr(Nxray)
 hotgas_counts = !values.D_NAN +  dblarr(Nxray)

 if strupcase(xray_unit) eq 'COUNTS' then begin
   for i=0, Nxray - 1 do begin
     xrb_counts[i] = xray_exposure[i] * trap_int(reverse(nu_obs), reverse(model_xrb_count_rate), $
                                                     XRANGE=xray_bandpass_nu[*, i])
     hotgas_counts[i] = xray_exposure[i] * trap_int(reverse(nu_obs), reverse(model_hotgas_count_rate), $
                                                        XRANGE=xray_bandpass_nu[*, i])
   endfor
 endif



 case strupcase(XRAY_AGN_MODEL) of
   'QSOSED': model_agn_count_rate = qsosed_models(redshift, wave_rest, arf_interp, exp_neg_tau_xray_MW, galactic_nH, $
                                                  lumin_dist, Lnu_x_agn=Lnu_x_agn, L2500=L2500, $
                                                  agn_mass=agn_mass, agn_logmdot=agn_logmdot)

   'PLAW': begin
       Lnu_x_agn = xray_plaw_expcut(wave_rest, plaw_gamma=1.8, E_cut=300.0)
       model_agn_count_rate = arf_interp * (exp_neg_tau_xray_MW ^ (galactic_nH) * Lnu_x_agn * LtoF_over_photon_energy_obs)
       L2500 = !values.D_NaN
       agn_mass = !values.D_NaN
       agn_logmdot = !values.D_NaN
     end

   'NONE': begin
       Lnu_x_agn = !values.D_NaN
       model_agn_count_rate = !values.D_NaN
       L2500 = !values.D_NaN
       agn_mass = !values.D_NaN
       agn_logmdot = !values.D_NaN
     end
 endcase


 xray_models = {XRAY_BANDPASS: xray_bandpass,             $
                XRAY_EXPOSURE: xray_exposure,             $
                WAVE_REST: wave_rest,                     $
                WAVE_OBS: wave_obs,                       $
                EXP_NEG_TAU_XRAY: exp_neg_tau_xray,       $
                EXP_NEG_TAU_XRAY_MW: exp_neg_tau_xray_MW, $
                LNU_AGN: Lnu_x_agn,                       $
                LNU_XRB: Lnu_x_xrb,                       $
                LNU_GAS: Lnu_x_hotgas,                    $
                AGN_MODEL: model_agn_count_rate,          $
                XRB_MODEL: model_xrb_count_rate,          $
                GAS_MODEL: model_hotgas_count_rate,       $
                XRB_COUNTS_LNU: xrb_counts,               $
                GAS_COUNTS_LNU: hotgas_counts,            $
                GALACTIC_NH: galactic_nH,                 $
                L2500: L2500,                             $
                AGN_MASS: agn_mass,                       $
                AGN_LOGMDOT: agn_logmdot,                 $
                REDSHIFT: redshift,                       $
                XRAY_ABS_MODEL: xray_abs_model,       $
                XRAY_AGN_MODEL: xray_agn_model}


 return, xray_models

end
