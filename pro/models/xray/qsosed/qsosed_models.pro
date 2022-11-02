function qsosed_models, redshift, wave_rest, arf_interp, exp_neg_tau_xray_MW, galactic_nH, $
                        lumin_dist, Lnu_x_agn=Lnu_x_agn, L2500=L2500, $
                        agn_mass=agn_mass, agn_logmdot=agn_logmdot

;+
; Name
; ----
;   QSOSED_MODELS
;
; Purpose
; -------
;   Generates the spectra and X-ray parameters for the QSOSED model as
;   described in Kubota & Done (2018). The spectra can include or not
;   include Galactic absorption.
;
; Calling Sequence
; ----------------
;   ::
;
;       model_agn_count_rate = qsosed_models(redshift, wave_rest, arf_interp, exp_neg_tau_xray_MW, $
;                                   galactic_nH, lumin_dist [, Lnu_x_agn=Lnu_x_agn, $
;                                   L2500=L2500, agn_mass=agn_mass, agn_logmdot=agn_logmdot])
;
; Inputs
; ------
;   ``redshift`` : int, float, or double scalar
;       The redshift of the model.
;   ``wave_rest`` : int, float or double array(Nwave)
;       The rest-frame wavelengths at which to interpolate the qsosed models
;       :math:`[\mu \rm m]`.
;   ``arf_interp`` : double array(Nwave)
;       A grid of ARF values interpolated from the input ARF file :math:`[{\rm cm}^2]`.
;   ``exp_neg_tau_xray_MW`` : float or double array(Nwave)
;       The attenuation from the Milky Way to be applied to the X-ray data in terms
;       of :math:`e^{-\tau}`.
;   ``galactic_nH`` : int, float, or double scalar
;       Galactic, i.e. Milky Way, neutral Hydrogen column density along the line
;       of sight :math:`[10^{20}\ {\rm cm}^2]`.
;   ``lumin_dist`` : int, float, double scalar
;       The luminosity distance of the model :math:`[{\rm Mpc}]`.
;
; Output
; ------
;   ``model_agn_count_rate`` : double array(Nwave, Nmass, Nmdot)
;       Instrumental count rate density produced by qsosed model grid
;       :math:`[{\rm counts\ s^{-1}\ Hz^{-1}}]`.
;
; Optional Outputs
; ----------------
;   ``Lnu_x_agn`` : double array(Nwave, Nmass, Nmdot)
;       The qsosed model luminosity :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   ``L2500`` : double array(Nmass, Nmdot)
;       The rest-frame 2500 Angstrom monochromatic luminosity shifted to the
;       observed frame :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   ``agn_mass`` : double array(Nmass, Nmdot)
;       The SMBH mass grid :math:`[M_\odot]`.
;   ``agn_logmdot`` : double array(Nmass, Nmdot)
;       The Log10 of SMBH accretion rate grid, normalized by the Eddington rate.
;
; Reference
; ---------
;   `Kubota, A., & Done, C. 2018, MNRAS, 480, 1247 <https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.1247K/abstract>`_
;
; Modification History
; --------------------
;   - 2021/09/21: Created (Erik B. Monson).
;   - 2022/06/22: Moved to separate file and documentation improved (Keith Doore).
;   - 2022/06/22: Major update to include new implementation (Keith Doore)
;   - 2022/11/02: Galactic NH is now in units of 1e20 cm-2 (Erik B. Monson)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(redshift) eq 0 then message, 'Variable is undefined: REDSHIFT.'
 if size(redshift, /type) lt 2 or size(redshift, /type) gt 5 then $
   message, 'REDSHIFT must be of type int, float, or double.'
 if size(redshift, /n_dim) ne 0  then $
   message, 'REDSHIFT must be a scalar.'
 if min(redshift) lt 0 then message, 'REDSHIFT must only contain non-negative values.'

 if n_elements(wave_rest) eq 0 then message, 'Variable is undefined: WAVE_REST.'
 if size(wave_rest, /type) lt 2 or size(wave_rest, /type) gt 5 then $
   message, 'WAVE_REST must be of type int, float, or double.'
 if size(reform(wave_rest), /n_dim) ne 1  then $
   message, 'WAVE_REST must be a scalar or 1-D array.'
 if min(wave_rest) le 0 then message, 'WAVE_REST must only contain positive values.'

 if n_elements(arf_interp) eq 0 then message, 'Variable is undefined: ARF_INTERP.'
 if size(arf_interp, /type) lt 2 or size(arf_interp, /type) gt 5 then $
   message, 'ARF_INTERP must be of type int, float, or double.'
 if size(reform(arf_interp), /n_dim) ne 1  then $
   message, 'ARF_INTERP must be a scalar or 1-D array.'
 if n_elements(arf_interp) ne n_elements(wave_rest) then $
   message, 'ARF_INTERP must have the same number of elements as WAVE_REST.'

 if n_elements(exp_neg_tau_xray_MW) eq 0 then message, 'Variable is undefined: EXP_NEG_TAU_XRAY_MW.'
 if size(exp_neg_tau_xray_MW, /type) lt 2 or size(exp_neg_tau_xray_MW, /type) gt 5 then $
   message, 'EXP_NEG_TAU_XRAY_MW must be of type int, float, or double.'
 if size(reform(exp_neg_tau_xray_MW), /n_dim) ne 1  then $
   message, 'EXP_NEG_TAU_XRAY_MW must be a scalar or 1-D array.'
 if n_elements(exp_neg_tau_xray_MW) ne n_elements(wave_rest) then $
   message, 'EXP_NEG_TAU_XRAY_MW must have the same number of elements as WAVE_REST.'

 if n_elements(galactic_nH) eq 0 then message, 'Variable is undefined: GALACTIC_NH.'
 if size(galactic_nH, /type) lt 2 or size(galactic_nH, /type) gt 5 then $
   message, 'GALACTIC_NH must be of type int, float, or double.'
 if size(galactic_nH, /n_dim) ne 0  then $
   message, 'GALACTIC_NH must be a scalar.'
 if min(galactic_nH) lt 0 then message, 'GALACTIC_NH must only contain non-negative values.'

 if n_elements(lumin_dist) eq 0 then message, 'Variable is undefined: LUMIN_DIST.'
 if size(lumin_dist, /type) lt 2 or size(lumin_dist, /type) gt 5 then $
   message, 'LUMIN_DIST must be of type int, float, or double.'
 if size(lumin_dist, /n_dim) ne 0  then $
   message, 'LUMIN_DIST must be a scalar.'
 if lumin_dist le 0 then $
   message, 'LUMIN_DIST must be a positive value.'



; Generate QSOSED xray models
 Nwave = n_elements(wave_rest)

 qsosed = mrdfits(!lightning_dir + 'models/xray/qsosed/qsosed.fits.gz', 1, qsosed_hdr, /SILENT)

 mass_card_idx = (where(stregex(qsosed_hdr, 'N_MASS') ne -1, /NULL))[0]
 Nmass = fix(strtrim((strsplit(qsosed_hdr[mass_card_idx], '=', /EXTRACT))[1], 2))

 mdot_card_idx = (where(stregex(qsosed_hdr, 'N_MDOT') ne -1, /NULL))[0]
 Nmdot = fix(strtrim((strsplit(qsosed_hdr[mdot_card_idx], '=', /EXTRACT))[1], 2))

 Nenergy = (size(qsosed.LNU, /DIMENSIONS))[0]
 ; Rest-frame energy and wavelength
 qsosed_energy = (qsosed.E_MID)[*, 0]
 qsosed_wave = (qsosed.WAVE_MID)[*, 0]

 L2500 = qsosed.L2500
 qsosed_gamma = qsosed.GAMMA_EFFECTIVE

 ; All this transposing is to get us to M begin the first dimension and
 ; log(mdot) the second. Array indexing is the other way around from python,
 ; where these models were generated in Sherpa.
 qsosed_matr = transpose(reform(transpose(qsosed.LNU), Nmdot, Nmass, 250), [2,1,0])
 M_matr = transpose(reform(transpose(qsosed.mass), Nmdot, Nmass))
 logmdot_matr = transpose(reform(transpose(qsosed.LOG_MDOT), Nmdot, Nmass))
 L2500_matr = transpose(reform(transpose(qsosed.L2500), Nmdot, Nmass))

 ; Interpolate to our normal wavelength grid
 qsosed_matr_interp = dblarr(Nwave, Nmass, Nmdot)
 for i=0, Nmass-1 do begin
   for j=0, Nmdot-1 do begin
     qsosed_matr_interp[*, i, j] = interpol(reform(qsosed_matr[*, i, j]), qsosed_wave, wave_rest)
     out_of_bounds = (wave_rest lt min(qsosed_wave)) or (wave_rest gt max(qsosed_wave))
     zero_idcs = where((out_of_bounds or (qsosed_matr_interp[*, i, j] lt 0)), /NULL)
     qsosed_matr_interp[zero_idcs, i, j] = 0
   endfor
 endfor


 ; We need the luminosity distance in cm to convert Lnu to fnu to counts s-1 Hz-1
 DL = lumin_dist
 if (DL lt 10) then DL = 10 ; We can't really handle very small distances. 10 Mpc is, I think, like z = 0.0023
 DL = DL * 1.d6 * !lightning_cgs.pc
 LtoF = 1.d0 / (4.d0 * !dpi * DL^2)
 ; Planck constant * speed of light is the conversion between wavelength and photon energy.
 hc = (!lightning_cgs.hplanck / !lightning_cgs.keV) * (1e4 * !lightning_cgs.clight) ; h * c in keV * micron
 ; Convert from energy units to photons to counts (recall that the ARF is counts photon-1 cm2)
 ; We want the photon energy in Lsun s for convenience later.
 wave_X_obs = (1 + redshift) * wave_rest
 photon_energy_obs = hc / wave_X_obs * (!lightning_cgs.keV / !lightning_cgs.Lsun)

 ; 1 / (4 pi DL^2) is a *very* small number in cm-2
 ; as is the energy of an X-ray photon in Lsun s, so we combine them
 ; here to avoid underflowing our models to 0
 LtoF_over_photon_energy_obs = LtoF / photon_energy_obs


 ; Inside the loops as we sample parameter space we'll need to further multiply the AGN model by
 ; exptau_xray ^ (NH_AGN), changing the shape of the model.
 model_agn_count_rate = rebin(reform(arf_interp * (exp_neg_tau_xray_MW ^ (galactic_nH)), Nwave, 1, 1), $
                          Nwave, Nmass, Nmdot) * qsosed_matr_interp * $
                    rebin(reform(LtoF_over_photon_energy_obs, Nwave, 1, 1), Nwave, Nmass, Nmdot)


; Compute additional outputs, which need to be in the observed-frame
 Lnu_x_agn = (1 + redshift) * qsosed_matr_interp
 L2500 = (1 + redshift) * L2500_matr
 agn_mass = M_matr
 agn_logmdot = logmdot_matr


 return, (1 + redshift) * model_agn_count_rate

end
