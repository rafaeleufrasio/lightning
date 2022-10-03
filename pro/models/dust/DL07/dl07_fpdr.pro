function dl07_fpdr, alpha, umin, umax, gam, error_check=error_check
;+
; Name
; ----
;   DL07_FPDR
;
; Purpose
; -------
;   Computes fPDR (fraction of total luminosity from the photodissociation
;   regions) according to Eq. 4 of Aniano et al. (2020). This calculation
;   assumes the parameters are defined as in Draine & Li (2007).
;
; Calling Sequence
; ----------------
;   ::
;
;       fPDR = dl07_fpdr(alpha, umin, umax, gam [, /error_check])
;
; Inputs
; ------
;   ``alpha`` : int, float, or double array(Nmodels)
;       The exponent of the power-law distribution of heating starlight
;       intensities between ``umin`` and ``umax``.
;   ``umin`` : int, float, or double array(Nmodels)
;       The minimum radiation field intensity of the diffuse ISM radiation
;       field heating the dust.
;   ``umax`` : int, float, or double array(Nmodels)
;       The maximum radiation field intensity of the power-law distribution
;       of heating starlight intensities.
;   ``gam`` : int, float, or double array(Nmodels)
;       The fraction of the dust mass exposed to the power-law distribution of
;       starlight intensities.
;
; Optional Input
; --------------
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``fPDR`` : double array(Nmodel, ...)
;       The fraction of the total dust luminosity that is radiated by
;       dust in regions where ``U > 1e2`` (i.e., fraction of total luminosity
;       from the photodissociation regions).
;
; References
; ----------
;   - `Draine, B. T., & Li, A. 2007, ApJ, 657, 810 <https://ui.adsabs.harvard.edu/abs/2007ApJ...657..810D/abstract>`_
;   - `Aniano, G., Draine, B. T., Hunt, L. K., et al. 2020, ApJ, 889, 150 <https://ui.adsabs.harvard.edu/abs/2020ApJ...889..150A/abstract>`_
;
; Modification History
; --------------------
;   - 2020/08/19: Created (Rafael T. Eufrasio)
;   - 2022/03/15: Added error handling (Keith Doore)
;   - 2022/03/18: Updated documentation (Keith Doore)
;   - 2022/04/12: Allowed for inputs to have degenerate dimensions (Keith Doore)
;   - 2022/04/12: Allowed integer inputs (Keith Doore)
;   - 2022/07/12: Renamed reformed ``gam`` input as to not change input value (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(alpha) eq 0 then message, 'Variable is undefined: ALPHA.'
   if size(alpha, /type) lt 2 or size(alpha,/type) gt 5 then message, 'ALPHA must be of type int, float, or double.'
   if size(reform(alpha), /n_dim) lt 1 then message, 'ALPHA must be an array.'
  
   if n_elements(umin) eq 0 then message, 'Variable is undefined: UMIN.'
   if size(umin, /type) lt 2 or size(umin,/type) gt 5 then message, 'UMIN must be of type int, float, or double.'
   if size(reform(umin), /n_dim) lt 1 then message, 'UMIN must be an array.'
  
   if n_elements(umax) eq 0 then message, 'Variable is undefined: UMAX.'
   if size(umax, /type) lt 2 or size(umax,/type) gt 5 then message, 'UMAX must be of type int, float, or double.'
   if size(reform(umax), /n_dim) lt 1 then message, 'UMAX must be an array.'
   if max(umin/umax) gt 1 then message, 'UMAX must only contain values larger than UMIN.'
  
   if n_elements(gam) eq 0 then message, 'Variable is undefined: GAM.'
   if size(gam, /type) lt 2 or size(gam,/type) gt 5 then message, 'GAM must be of type int, float, or double.'
   if size(reform(gam), /n_dim) lt 1 then message, 'GAM must be an array.'
   if min(gam) lt 0 or max(gam) gt 1 then message, 'GAM must contain values between 0 and 1.'
  
   if n_elements(alpha) ne n_elements(umin) or n_elements(alpha) ne n_elements(umax) or $
      n_elements(alpha) ne n_elements(gam) then $
      message,'ALPHA, UMIN, UMAX, and GAM must have the same size.'
 endif

; Compute fPDR
 hund_umax = 1.d2 / reform(umax)
 umin_umax = reform(umin / umax)
 alpha_minus_one = reform(alpha) - 1.d0
 two_minus_alpha = 2.d0 - reform(alpha)
 gamm = reform(gam)

 num = dblarr(size(gamm, /dim))
 den = dblarr(size(gamm, /dim))

 alpha_eq2 = where(alpha eq 2.d0, Neq2, /null)
 alpha_ne2 = where(alpha ne 2.d0, Nne2, /null)

 ;for alpha = 2
 IF Neq2 ne 0 then begin
   num[alpha_eq2] = -gamm[alpha_eq2]*alog(hund_umax[alpha_eq2])
   den[alpha_eq2] = (1d0 - gamm[alpha_eq2])*(1d0 - Umin_Umax[alpha_eq2]) - gamm[alpha_eq2]*alog(Umin_Umax[alpha_eq2])
 endif

 ;for alpha != 2
 if Nne2 ne 0 then begin
   num[alpha_ne2] = gamm[alpha_ne2] * (1d0 - (hund_umax[alpha_ne2])^two_minus_alpha[alpha_ne2])
   den[alpha_ne2] = (1d0-gamm[alpha_ne2]) * two_minus_alpha[alpha_ne2]/alpha_minus_one[alpha_ne2]   $
                               * (Umin_Umax[alpha_ne2]^two_minus_alpha[alpha_ne2])       $
                               * (1d0 - umin_umax[alpha_ne2]^alpha_minus_one[alpha_ne2]) $
                + gamm[alpha_ne2]    * (1d0 - umin_umax[alpha_ne2]^two_minus_alpha[alpha_ne2])
 endif

 fPDR = num/den
 fPDR[where(alpha_minus_one eq 0, /null)] = 1.d0

 if n_elements(fPDR) eq 1 then return, fPDR[0] else $
   return, fPDR

end
