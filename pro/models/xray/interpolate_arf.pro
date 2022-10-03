function interpolate_arf, E_lo, E_hi, specresp, wave
;+
; Name
; ----
;   INTERPOLATE_ARF
;
; Purpose
; -------
;   Gets the values of the specified Auxiliary Response Function (ARF)
;   at the specified wavelengths.
;
; Calling Sequence
; ----------------
;   ::
;
;       specresp_interp = interpolate_arf(E_lo, E_hi, specresp, wave)
;
; Inputs
; ------
;   ``E_lo`` : int, float or double array(Nchannels)
;       Lower energy bounds of each channel in the ARF :math:`[{\rm keV}]`
;   ``E_hi`` : int, float or double array(Nchannels)
;       Upper energy bounds of each channel in the ARF :math:`[{\rm keV}]`
;   ``specresp`` : int, float or double array(Nchannels)
;       The spectral response of the ARF at each channel :math:`[{\rm cm}^2]`
;   ``wave`` : int, float or double array(Nwave)
;       A grid of wavelengths at which to interpolate the ARF :math:`[\mu \rm m]`.
;
; Output
; ------
;   ``specresp_interp`` : double array(Nwave)
;       A grid of ARF values interpolated from the input ARF file at each
;       specified wavelength :math:`[{\rm cm}^2]`.
;
; Modification History
; --------------------
;   - 2021/08/26: Created (Erik B. Monson).
;   - 2022/06/22: Added documentation (Keith Doore)
;   - 2022/06/22: Replaced ``!cv`` with ``!lightning_cgs`` (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(E_lo) eq 0 then message, 'Variable is undefined: E_LO.'
 if size(E_lo, /type) lt 2 or size(E_lo, /type) gt 5 then $
   message, 'E_LO must be of type int, float, or double.'
 if size(reform(E_lo), /n_dim) ne 1  then $
   message, 'E_LO must be a scalar or 1-D array.'
 if min(E_lo) le 0 then message, 'E_LO must only contain positive values.'

 if n_elements(E_hi) eq 0 then message, 'Variable is undefined: E_HI.'
 if size(E_hi, /type) lt 2 or size(E_hi, /type) gt 5 then $
   message, 'E_HI must be of type int, float, or double.'
 if size(reform(E_hi), /n_dim) ne 1  then $
   message, 'E_HI must be a scalar or 1-D array.'
 if min(E_hi) le 0 then message, 'E_HI must only contain positive values.'
 if n_elements(E_hi) ne n_elements(E_lo) then message, 'E_HI must have the same number of elements as E_LO.'

 if n_elements(specresp) eq 0 then message, 'Variable is undefined: SPECRESP.'
 if size(specresp, /type) lt 2 or size(specresp, /type) gt 5 then $
   message, 'SPECRESP must be of type int, float, or double.'
 if size(reform(specresp), /n_dim) ne 1  then $
   message, 'SPECRESP must be a 1-D array.'

 if n_elements(wave) eq 0 then message, 'Variable is undefined: WAVE.'
 if size(wave, /type) lt 2 or size(wave, /type) gt 5 then $
   message, 'WAVE must be of type int, float, or double.'
 if size(reform(wave), /n_dim) ne 1  then $
   message, 'WAVE must be a scalar or 1-D array.'
 if min(wave) le 0 then message, 'WAVE must only contain positive values.'


; Interpolate ARF
 ; hc in keV * micron
 hc = (!lightning_cgs.hplanck / !lightning_cgs.keV) * (1e4 * !lightning_cgs.clight)

 wave_lo = hc / E_hi
 wave_hi = hc / E_lo

 specresp_interp = dblarr(n_elements(wave))

; Find where input wavelength is in ARF energy bin range
 for i=0, n_elements(wave) - 1 do begin
   idx = where((wave_lo lt wave[i]) and (wave_hi ge wave[i]), /null, num)
   if num ne 0 then specresp_interp[i] = specresp[idx]
 endfor

 return, specresp_interp

end
