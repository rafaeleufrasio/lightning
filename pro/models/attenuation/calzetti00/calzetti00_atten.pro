function calzetti00_atten, wave, tauV_diff=tauV_diff, delta=delta, $
                           tauV_BC=tauV_BC, uv_bump=uv_bump, error_check=error_check
;+
; Name
; ----
;   CALZETTI00_ATTEN
;
; Purpose
; -------
;   Generates attenuation values at the input wavelengths using the 
;   Calzetti et al. (2000) attenuation curve or its modified version.
;   Besides the normalization in the standard Calzetti curve, the 
;   modified version includes a variable slope as described in 
;   Noll et al. (2009), an optional 2175 Angstrom bump feature specified
;   in Kriek & Conroy (2013), and birth cloud attenuation as described in
;   Eufrasio et al. (2017).
;
; Calling Sequence
; ----------------
;   ::
;
;       exp_neg_tau = calzetti00_atten(wave [, tauV_diff = , delta = , $
;                                      tauV_BC = , /uv_bump, /error_check])
;
; Input
; -----
;   ``wave`` : int, float, or double array(Nwave)
;       The wavelength at which to determine the attenuation :math:`[\mu \rm m]`.
;
; Optional Inputs
; ---------------
;   ``tauV_diff`` : int, float, or double array(Nmodels)
;       The V-band optical depth of the diffuse dust. (Default = ``1.d0``)
;   ``delta`` : int, float, or double array(Nmodels)
;       The power law value to change the attenuation curve slope.
;       (Default = ``0.d0``)
;   ``tauV_BC`` : int, float, or double array(Nmodels)
;       The V-band optical depth of the birth cloud. (Default = ``0.d0``)
;   ``uv_bump`` : flag
;       If set, then a 2175 Angstrom UV bump feature will be added to the 
;       attenuation curve.
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``exp_neg_tau`` : double array(Nwave, Nmodels)
;       The attenuation in terms of :math:`e^{-\tau}` (i.e., :math:`e^{-\tau} = 10^{-0.4 A_{\lambda}}`).
;
; Note
; ----
;   Defaults will only be set if the optional ``error_check`` input is set.
;
; References
; ----------
;   - `Calzetti, D., Armus, L., Bohlin, R. C., et al. 2000, ApJ, 533, 682 <https://ui.adsabs.harvard.edu/abs/2000ApJ...533..682C/abstract>`_
;   - `Noll, S., Burgarella, D., Giovannoli, E., et al. 2009, A&A, 507, 1793 <https://ui.adsabs.harvard.edu/abs/2009A%26A...507.1793N/abstract>`_
;   - `Kriek, M., & Conroy, C. 2013, ApJ, 775, L16 <https://ui.adsabs.harvard.edu/abs/2013ApJ...775L..16K/abstract>`_
;   - `Eufrasio, R. T., Lehmer, B. D., Zezas, A., et al. 2017, ApJ, 851, 10 <https://ui.adsabs.harvard.edu/abs/2017ApJ...851...10E/abstract>`_
;
; Modification History
; --------------------
;   - 2016/05/01: Created (Rafael T. Eufrasio)
;   - 2020/05/06: Added ability to run with pure Calzetti curve (Keith Doore)
;   - 2022/03/15: Added proper error handling (Keith Doore)
;   - 2022/03/15: Renamed variables to standard format (Keith Doore)
;   - 2022/03/15: Cleaned up method for separating pure and modified Calzetti (Keith Doore)
;   - 2022/03/18: Updated documentation (Keith Doore)
;   - 2022/04/07: Allowed for inputs to have degenerate dimensions (Keith Doore)
;   - 2022/04/07: Allowed for inputs to be scalars (Keith Doore)
;   - 2022/04/07: Allowed integer inputs (Keith Doore)
;   - 2022/04/07: Made optical depths positive input values (Keith Doore)
;   - 2022/04/07: Allowed for array of only one optional input to be given (Keith Doore)
;   - 2022/06/09: Changed ``no_bump`` keyword to ``uv_bump`` keyword and adjusted logic accordingly (Keith Doore)
;   - 2022/06/09: Added ``error_check`` keyword to do error handling (Keith Doore)
;   - 2022/11/04: Corrected extrapolation below 0.12um to use method from Noll+2009 (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(wave) eq 0 then message, 'Variable is undefined: WAVE.'
   if size(wave, /type) lt 2 or size(wave,/type) gt 5 then message, 'WAVE must be of type int, float, or double.'
   if size(reform(wave), /n_dim) ne 1 then message, 'WAVE must be a scalar or 1-D array.'
   if min(wave) le 0 then message,'WAVE must only contain positive values.'

   if n_elements(tauV_diff) ne 0 then begin
     if size(tauV_diff, /type) lt 2 or size(tauV_diff,/type) gt 5 then $
              message, 'TAUV_DIFF must be of type int, float, or double.'
     if size(reform(tauV_diff), /n_dim) ne 1 then message, 'TAUV_DIFF must be a scalar or 1-D array.'
     if min(tauV_diff) lt 0 then message, 'TAUV_DIFF must only contain non-negative values.'
     Nmodels = n_elements(tauV_diff)
   endif

   if n_elements(delta) ne 0 then begin
     if size(delta, /type) lt 2 or size(delta,/type) gt 5 then $
              message, 'DELTA must be of type int, float, or double.'
     if size(reform(delta), /n_dim) ne 1 then message, 'DELTA must be a scalar or 1-D array.'
     if n_elements(Nmodels) eq 0 then Nmodels = n_elements(delta)
   endif

   if n_elements(tauV_BC) ne 0 then begin
     if size(tauV_BC, /type) lt 2 or size(tauV_BC,/type) gt 5 then $
              message, 'TAUV_BC must be of type int, float, or double.'
     if size(reform(tauV_BC), /n_dim) ne 1 then message, 'TAUV_BC must be a scalar or 1-D array.'
     if min(tauV_BC) lt 0 then message, 'TAUV_BC must only contain non-negative values.'
     if n_elements(Nmodels) eq 0 then Nmodels = n_elements(tauV_BC)
   endif

   if n_elements(Nmodels) ne 0 then begin
     if n_elements(tauV_diff) eq 0 then tauV_diff = replicate(1.d0,Nmodels)
     if n_elements(delta    ) eq 0 then delta     = replicate(0.d0,Nmodels)
     if n_elements(tauV_BC  ) eq 0 then tauV_BC   = replicate(0.d0,Nmodels)
   endif else begin
     Nmodels = 1
     tauV_diff = 1.d0
     delta     = 0.d0
     tauV_BC   = 0.d0
   endelse

   if n_elements(tauV_diff) ne n_elements(delta) or n_elements(tauV_diff) ne n_elements(tauV_BC) then $
      message,'TAUV_DIFF, DELTA, and TAUV_BC must have the same size.'

 endif
 Nmodels = n_elements(tauV_diff)


; Derive exp_neg_tau at input wavelength
 wavelength=reform(wave)
 Nwave = n_elements(wavelength)

; UV bump feature
 Dlam=dblarr(nwave, Nmodels)
 if keyword_set(uv_bump) then begin
   FWHM_bump = 0.0350d0 ; 350 Angs
   lambda0_bump = 0.2175d0   ; 2175 Angs
   Dlam = (rebin(reform(0.85d0-1.9d0*delta, 1, Nmodels), nwave, Nmodels)) / $
          (rebin((((wavelength^2 - lambda0_bump^2)/(wavelength*FWHM_bump))^2 + 1.d0), nwave, Nmodels))
 endif


;Calzetti extinction curve originally derived from 0.0912 to 2.2 microns
 klam = dblarr(Nwave)
 RV = 4.05d0 ;Calzetti+2000

 w1 = where((wavelength GE 0.6300) AND (wavelength LE 2.2000), nw1, /null)
 w2 = where((wavelength GE 0.1200) AND (wavelength LT 0.6300), nw2, /null)
 inverse_wave  = 1.d0/wavelength  ;Wavelength in inverse microns
 if (nw1 ne 0) then klam[w1] = 2.659d0*(-1.857d0 + 1.040d0*inverse_wave[w1]) + RV
 if (nw2 ne 0) then klam[w2] = 2.659d0*(poly(inverse_wave[w2], [-2.156, 1.509d0, -0.198d0, 0.011d0])) + RV

 ; Linearly extrapolate from 0.12 to 0.0912 as stated in Noll+2009
 ;slope = 2.659d0*(-1.509/0.12^2+ 2*0.198/0.12^3 - 3*0.011/0.12^4) = -92.449445
 ;value = 2.659d0*(poly(1/0.12, [-2.156, 1.509d0, -0.198d0, 0.011d0])) + RV = 12.119377
 w3 = where((wavelength GE 0.0912) AND (wavelength LT 0.1200), nw3, /null)
 slope = -92.449445d
 value =  12.119377d
 if (nw3 ne 0) then klam[w3] = slope * wavelength[w3] + (value - slope * 0.1200)

 ; Combine bump feature, variable slope, and Calzetti curve
 ; If uv_bump is not set and delta is 0, then this results in original Calzetti curve
 flam = (rebin(klam, nwave, Nmodels) + temporary(Dlam))/4.05 * $
        (rebin(wavelength, nwave, Nmodels)/0.55)^(rebin(reform(delta, 1, Nmodels), nwave, Nmodels))
 tau_lam_diff = temporary(flam)*rebin(reform(tauV_diff, 1, Nmodels), nwave, Nmodels)

 ; birth cloud attenuation
 tau_lam_BC = rebin(1.d0/(wavelength/0.55d0), nwave, Nmodels)*rebin(reform(tauV_BC, 1, Nmodels), nwave, Nmodels)

 ; Combine Calzetti component and birth cloud component
 exp_neg_tau = exp(-1.d0*(temporary(tau_lam_diff) + temporary(tau_lam_BC)))

 if Nmodels eq 1 and Nwave eq 1 then return, exp_neg_tau[0] else $
   return, exp_neg_tau

end