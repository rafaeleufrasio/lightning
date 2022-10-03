function gilbertson22_LX_tau, stellar_models, psi, error_check=error_check, $
                              Lnu_LMXB_steps=Lnu_LMXB_steps, Lnu_HMXB_steps=Lnu_HMXB_steps
;+
; Name
; ----
;   GILBERTSON22_LX_TAU
;
; Purpose
; -------
;   Produces the 2-10 keV LMXB and HMXB luminosities appropriate for the stellar 
;   population(s), according to the :math:`L_X/M` parametrizations with stellar age 
;   (that is, :math:`\tau = \log(\rm age)`) from Gilbertson et al. (2022).
;
; Calling Sequence
; ----------------
;   ::
;
;       Lnu_XRB = gilbertson22_LX_tau(stellar_models, psi [, /error_check, $
;                        Lnu_LMXB_steps=Lnu_LMXB_steps, Lnu_HMXB_steps=Lnu_HMXB_steps])
;
; Inputs
; ------
;   ``stellar_models`` : structure
;       A structure containing the spectra, SEDs, and stellar parameters for the 
;       non-parametric stellar model. (See ``binned_stellar_models.pro`` for
;       details and contents.)
;   ``psi`` : int, float, or double array(Nsteps, Nmodels)
;       The non-parametric SFH coefficients :math:`[M_\odot\ {\rm yr}^{-1}]`.
;
; Optional Input
; --------------
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``Lnu_XRB`` : double array(2, Nmodels) 
;       The total LMXB (first column) and HMXB (second column) 2-10 keV
;       luminosities :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;
; Optional Outputs
; ----------------
;   ``Lnu_LMXB_steps`` : double array(Nsteps, Nmodels)
;       2-10 keV Lnu from LMXBs for each stellar age bin :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   ``Lnu_HMXB_steps`` : double array(Nsteps, Nmodels)
;       2-10 keV Lnu from HMXBs for each stellar age bin :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;
; Notes
; -----
;   The parametrizations of :math:`\gamma = L_X / M` as a function of stellar age are as follows:
;
;       :math:`\gamma(t_{\rm age}) = \gamma_{\rm HMXB}(t_{\rm age}) + \gamma_{\rm LMXB}(t_{\rm age})
;       = [a_0(\tau - b_0)^2 + c_0] + [a_1(\tau - b_1)^2 + c_1]`
;
;   where :math:`\tau = \log_{10}(t_{\rm age})`. Gilbertson et al. (2022) has:
;
;       - :math:`a_0 = -0.24 \ [{\rm erg s^{-1}\ M_\odot^{-1}}\ \log({\rm yr})^{-1}]`
;       - :math:`b_0 = 5.23 \ [\log({\rm yr})^{-1}]`
;       - :math:`c_0 = 32.54 \ [{\rm erg s^{-1}\ M_\odot^{-1}}]`
;       - :math:`a_1 = -1.21 \ [{\rm erg s^{-1}\ M_\odot^{-1}}\ \log({\rm yr})^{-1}]`
;       - :math:`b_1 = 9.32 \ [\log({\rm yr})^{-1}]`
;       - :math:`c_1 = 29.09 \ [{\rm erg s^{-1}\ M_\odot^{-1}}]`
;
; Reference
; ---------
;   `Gilbertson, W., Lehmer, B. D., Doore, K., et al. 2022, ApJ, 926, 28 <https://ui.adsabs.harvard.edu/abs/2022ApJ...926...28G/abstract>`_
;
; Modification History
; --------------------
;   - 2022/04/21: Created (Erik B. Monson)
;   - 2022/06/20: Renamed variables to common naming scheme (Keith Doore)
;   - 2022/06/20: Replaced ``!cv`` with ``!lightning_cgs`` (Keith Doore)
;   - 2022/06/20: Updated documentation (Keith Doore)
;   - 2022/06/20: Added error handling (Keith Doore)
;   - 2022/06/20: Added ``error_check`` keyword to do error handling (Keith Doore)
;   - 2022/07/12: Removed reforming of ``psi`` if ``Nmodels=1`` as it is not necessary (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(stellar_models) eq 0 then message, 'Variable is undefined: STELLAR_MODELS.'
   if size(stellar_models, /type) ne 8 then message, 'STELLAR_MODELS must be of type structure.'
  
   if n_elements(psi) eq 0 then message, 'Variable is not defined: PSI.'
   if size(psi, /type) lt 2 or size(psi, /type) gt 5 then $
     message, 'PSI must be of type int, float, or double.'
   if size(psi, /n_dim) gt 2 then $
     message, 'PSI must be a scalar, 1-D, or 2-D array.'
 endif


 Nsteps = n_elements(stellar_models.bounds) - 1
 if size(psi, /n_dim) lt 2 then Nmodels = 1 else Nmodels = (size(psi, /dim))[1]


 ; We need the current stellar mass in each bin to calculate the HMXB and LMXB luminosities.
 steps_Mstar = rebin(stellar_models.MSTAR, Nsteps, Nmodels) * psi
 ; The parametrizations are a function of stellar age
 avg_age_per_bin = 0.5 * (stellar_models.BOUNDS[0:-2] + stellar_models.BOUNDS[1:-1])
 tau_age = alog10(rebin(avg_age_per_bin, Nsteps, Nmodels))

 ; These are rest frame 2-10 keV luminosities.
 dnu_2_10 = 8.d0 * !lightning_cgs.keV / !lightning_cgs.hplanck

 ; log(2-10 keV Lx / Mstar) converted to log(Lnu) in log(Lsun Hz-1) based on current psi
 Lnu_LMXB_steps = 10.d0^((-1.21 * (tau_age - 9.32) ^ 2 + 29.09) + alog10(steps_Mstar) - alog10(dnu_2_10) - alog10(!lightning_cgs.Lsun))
 Lnu_LMXB = total(Lnu_LMXB_steps, 1)

 ; log(2-10 keV Lx / Mstar) converted to log(Lnu) in log(Lsun Hz-1) based on current psi
 Lnu_HMXB_steps = 10.d0^((-0.24 * (tau_age - 5.23) ^ 2 + 32.54) + alog10(steps_Mstar) - alog10(dnu_2_10) - alog10(!lightning_cgs.Lsun))
 Lnu_HMXB = total(Lnu_HMXB_steps, 1)

 return, transpose([[Lnu_LMXB], [Lnu_HMXB]])

end