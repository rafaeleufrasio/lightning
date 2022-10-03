function qsosed_L2500, xray_models, agn_mass, agn_logmdot, error_check=error_check
;+
; Name
; ----
;   QSOSED_L2500
;
; Purpose
; -------
;   Generates the qsosed :math:`L_{2500}` (monochromatic luminosity at 2500 Angstroms)
;   corresponding to given values of :math:`M` and :math:`\log(\dot{M})` using the
;   models from Kubota & Done (2018).
;
; Calling Sequence
; ----------------
;   ::
;
;       L2500 = qsosed_L2500(xray_models, agn_mass, agn_logmdot [, /error_check])
;
; Inputs
; ------
;   ``xray_models`` : structure
;       A structure containing the spectra, counts, and X-ray model parameters.
;       (See ``xrb_xagn_models.pro`` for details and contents.)
;   ``agn_mass`` : int, float, or double array(Nmodels)
;       Supermassive black hole mass :math:`[M_\odot]`.
;   ``agn_logmdot`` : int, float, or double array(Nmodels)
;       Log10 of SMBH accretion rate, normalized by the Eddington rate.
;
; Optional Input
; --------------
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``L2500`` : double array(Nmodels)
;       The interpolated rest-frame 2500 Angstrom monochromatic luminosity of 
;       shifted to the observed frame :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;
; Reference
; ---------
;   `Kubota, A., & Done, C. 2018, MNRAS, 480, 1247 <https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.1247K/abstract>`_
;
; Modification History
; --------------------
;   - 2022/04/18: Created (Erik B. Monson)
;   - 2022/06/10: Added error handling (Keith Doore)
;   - 2022/06/10: Added ``error_check`` keyword to do error handling (Keith Doore)
;   - 2022/06/10: Made ``mass`` and ``logmdot`` required inputs (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(xray_models) eq 0 then message, 'Variable is undefined: XRAY_MODELS.'
   if size(xray_models, /type) ne 8 then message, 'XRAY_MODELS must be of type structure.'
  
   if n_elements(agn_mass) eq 0 then message, 'Variable is undefined: AGN_MASS.'
   if size(agn_mass, /type) lt 2 or size(agn_mass,/type) gt 5 then $
            message, 'AGN_MASS must be of type int, float, or double.'
   if size(reform(agn_mass), /n_dim) ne 1 then message, 'AGN_MASS must be a scalar or 1-D array.'
   if min(agn_mass) lt 1.d5 or max(agn_mass) gt 1.d10 then message, 'AGN_MASS must contain values between 1e5 and 1e10.'

   if n_elements(agn_logmdot) eq 0 then message, 'Variable is undefined: AGN_LOGMDOT.'
   if size(agn_logmdot, /type) lt 2 or size(agn_logmdot,/type) gt 5 then $
            message, 'AGN_LOGMDOT must be of type int, float, or double.'
   if size(reform(agn_logmdot), /n_dim) ne 1 then message, 'AGN_LOGMDOT must be a scalar or 1-D array.'
   if min(agn_logmdot) lt -1.5d or max(agn_logmdot) gt 0.3d then message, 'AGN_LOGMDOT  must contain values between -1.5 and 0.3.'
   if n_elements(agn_mass) ne n_elements(agn_logmdot) then $
            message, 'AGN_MASS and AGN_LOGMDOT must have the same number of elements.'
 endif

 Nmodels = n_elements(agn_mass)

 mass_min = min(xray_models.AGN_MASS)
 mass_max = max(xray_models.AGN_MASS)
 Nmass = (size(xray_models.AGN_MASS, /DIMENSIONS))[0]

 logmdot_min = min(xray_models.AGN_LOGMDOT)
 logmdot_max = max(xray_models.AGN_LOGMDOT)
 Nmdot = (size(xray_models.AGN_LOGMDOT, /DIMENSIONS))[1]

 M_interp_norm = (Nmass - 1) * (alog10(agn_mass) - alog10(mass_min)) / (alog10(mass_max) - alog10(mass_min))
 logmdot_interp_norm = (Nmdot - 1) * (agn_logmdot - logmdot_min) / (logmdot_max - logmdot_min)

 L2500 = interpolate(xray_models.L2500, M_interp_norm, logmdot_interp_norm)

 return, L2500

end