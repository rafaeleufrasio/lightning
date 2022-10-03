function xray_plaw_expcut, wave, plaw_gamma=plaw_gamma, E_cut=E_cut, normbox=normbox
;+
; Name
; ----
;   XRAY_PLAW_EXPCUT
;
; Purpose
; -------
;   Generates an X-ray power law spectrum with a high energy exponential 
;   cutoff (i.e., :math:`L_{\nu} \propto F_{\nu} \propto E^{-\gamma + 1}
;   e^{-E / E_{\rm cut}}`). The spectrum is normalized by default to the
;   monochromatic flux at 2 keV, computed in a 0.5 keV (~100 Angstrom) wide
;   box filter.
;
; Calling Sequence
; ----------------
;   ::
;
;       plaw_spec = xray_plaw_expcut(wave [, plaw_gamma = , E_cut = , normbox = ])
;
; Input
; -----
;   ``wave`` : int, float or double array(Nwave)
;       A grid of wavelengths at which evaluate the power law spectrum
;       :math:`[\mu \rm m]`.
;
; Optional Inputs
; ---------------
;   ``plaw_gamma`` : int, float or double scalar
;       The photon index of the power law. (Default = ``1.8``)
;   ``E_cut`` : int, float or double scalar
;       The high energy exponential cutoff value :math:`[{\rm keV}]`. (Default = ``300``)
;   ``normbox`` : int, float or double array(2)
;       The lower and upper energies of a box in which to normalize
;       the spectrum :math:`[{\rm keV}]`. (Default = ``[1.75, 2.25]``)
;
; Output
; ------
;   ``plaw_spec`` : double array(Nwave)
;       The normalized power law spectrum.
;
; Modification History
; --------------------
;   - 2021/03/19: Created (Erik B. Monson).
;   - 2022/06/22: Added documentation (Keith Doore)
;   - 2022/06/22: Replaced ``!cv`` with ``!lightning_cgs`` (Keith Doore)
;   - 2022/06/22: Updated documentation and renamed variables (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(wave) eq 0 then message, 'Variable is undefined: WAVE.'
 if size(wave, /type) lt 2 or size(wave, /type) gt 5 then $
   message, 'WAVE must be of type int, float, or double.'
 if size(reform(wave), /n_dim) ne 1 then $
   message, 'WAVE must be a scalar or 1-D array.'
 if min(wave) le 0 then message, 'WAVE must only contain positive values.'

 if n_elements(plaw_gamma) ne 0 then begin
   if size(plaw_gamma, /type) lt 2 or size(plaw_gamma, /type) gt 5 then $
     message, 'PLAW_GAMMA must be of type int, float, or double.'
   if size(plaw_gamma, /n_dim) ne 0 then $
     message, 'PLAW_GAMMA must be a scalar.'
 endif else plaw_gamma = 1.8d0

 if n_elements(E_cut) ne 0 then begin
   if size(E_cut, /type) lt 2 or size(E_cut, /type) gt 5 then $
     message, 'E_CUT must be of type int, float, or double.'
   if size(E_cut, /n_dim) ne 0 then $
     message, 'E_CUT must be a scalar.'
   if min(E_cut) le 0 then message, 'E_CUT must be a positive value.'
 endif else E_cut = 300.d0

 if n_elements(normbox) ne 0 then begin
   if size(normbox, /type) lt 2 or size(normbox, /type) gt 5 then $
     message, 'NORMBOX must be of type int, float, or double.'
   if size(normbox, /n_dim) ne 1 then $
     message, 'NORMBOX must be a 1-D array.'
   if min(normbox) le 0 then message, 'NORMBOX must only contain positive values.'
 endif else normbox = [1.75d0, 2.25d0]

 ; hc in keV * micron
 hc = (!lightning_cgs.hplanck / !lightning_cgs.keV) * (1e4 * !lightning_cgs.clight)

 if (hc / max(normbox) lt min(wave)) or (hc / min(normbox) gt max(wave)) then $
     message, 'NORMBOX values extend outside of the range of the model (1e-6 to 1e-3 micron).'


; Compute the power law spectrum
 E = (hc / wave)
 plaw_spec = (E ^ (1.0 - plaw_gamma)) * exp(-1.0 * E / E_cut)


 normbox_minwave = hc / max(normbox)
 normbox_maxwave = hc / min(normbox)
 filter = dblarr(n_elements(wave))
 mask = ((wave ge normbox_minwave) and (wave le normbox_maxwave))
 filter[where(mask, /null)] = 1
 filter = filter / trap_int(wave, filter)


 normalization = trap_int(wave, plaw_spec * filter)
 plaw_spec = plaw_spec / normalization


 return, plaw_spec

end