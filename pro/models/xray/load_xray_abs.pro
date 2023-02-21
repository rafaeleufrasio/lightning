function load_xray_abs, wave, xray_abs_model
;+
; Name
; ----
;   LOAD_XRAY_ABS
;
; Purpose
; -------
;   Loads the specified X-ray absorption model,and interpolates
;   it to the input wavelength grid.
;
; Calling Sequence
; ----------------
;   ::
;
;       exp_neg_tau_xray = load_xray_abs(wave, xray_abs_model)
;
; Inputs
; ------
;   ``wave`` : int, float or double array(Nwave)
;       A grid of wavelengths to be interpolated from the read in
;       X-ray absorption model :math:`[\mu \rm m]`.
;   ``xray_abs_model`` : string scalar
;       The name of the X-ray absorption model to apply to the X-ray emission.
;       Current options are ``'TBABS-WILM'``, ``'ATTEN'``,
;       and ``'NONE'``.
;
; Output
; ------
;   ``exp_neg_tau_xray`` : float or double array(Nwave)
;       The absorption in terms of :math:`e^{-\tau}` interpolated from the read
;       in X-ray absorption model.
;
; Modification History
; --------------------
;   - 2021/09/21: Created (Erik B. Monson).
;   - 2022/06/21: Added documentation (Keith Doore)
;   - 2022/06/21: Converted to using ``config`` from ``curve`` keyword (Keith Doore)
;   - 2022/07/01: Converted from ``config`` to ``xray_abs_model`` keyword since it was the only accessed value from ``config`` (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(wave) eq 0 then message, 'Variable is undefined: WAVE.'
 if size(wave, /type) lt 2 or size(wave, /type) gt 5 then $
   message, 'WAVE must be of type int, float, or double.'
 if size(reform(wave), /n_dim) ne 1  then $
   message, 'WAVE must be a scalar or 1-D array.'
 if min(wave) le 0 then message, 'WAVE must only contain positive values.'

 if n_elements(xray_abs_model) eq 0 then message, 'Variable is undefined: XRAY_ABS_MODEL.'
 if size(xray_abs_model, /type) ne 7 then message, 'XRAY_ABS_MODEL must be of type string.'
 if size(xray_abs_model, /n_dim) ne 0 then message, 'XRAY_ABS_MODEL must be a scalar.'
 if total(strupcase(xray_abs_model) eq ['TBABS-WILM', 'ATTEN', 'NONE']) ne 1 then $
   message, "XRAY_ABS_MODEL must be set to either 'TBABS-WILM', 'ATTEN', or 'NONE'."


; Read in specified absorption file
 case strupcase(xray_abs_model) of
   'TBABS-WILM': curve_path = !lightning_dir + 'models/xray/absorption/tbabs_wilm.txt'
   ;'TBABS-ANGR': curve_path = !lightning_dir + 'models/xray/absorption/tbabs_angr.txt'
   'ATTEN':      curve_path = !lightning_dir + 'models/xray/absorption/atten.txt'
   'NONE':       return, replicate(1.d0, n_elements(wave))
 endcase

 readcol, curve_path, wave_tab, exp_neg_tau_tab, comment='#', format='D,x,D', /silent

 ; The curves are dumped at the same wavelength grid as we use for
 ; the xray_models.pro. So, this interpolation should not really be
 ; necessary, but it's also relatively quick and we'll only do it once.
 exp_neg_tau_xray = interpol(exp_neg_tau_tab, wave_tab, wave)

 return, exp_neg_tau_xray

end
