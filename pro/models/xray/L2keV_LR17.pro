function L2keV_LR17, L2500, error_check=error_check
;+
; Name
; ----
;   L2KEV_LR17
;
; Purpose
; -------
;   Calculates the the 2 keV monochromatic luminosity of an AGN 
;   given its intrinsic monochromatic luminosity at 2500
;   Angstroms using the Lusso & Risaliti (2017) relationship:
;   
;       :math:`\log(L_{\rm 2keV}) = 0.633 (\log(L_{2500A}) - 25) - 1.959 + 25` 
;       :math:`\log(L_{\rm 2keV}) = 0.633 \log(L_{2500A}) + 7.216`
;   
;   where both monochromatic luminosities are in :math:`\rm ergs\ s^{-1}\ Hz^{-1}`.
;   
; Calling Sequence
; ----------------
;   ::
;   
;       L2keV = L2keV_LR17(L2500 [, /error_check])
;   
; Input
; -----
;   ``L2500`` : float or double array(...)
;       The rest-frame 2500 Angstrom monochromatic luminosity
;       :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   
; Optional Input
; --------------
;   ``error_check`` : flag
;       If set, the input is checked for errors. Otherwise, the input is
;       assumed to be of correct format.
;   
; Output
; ------
;   ``L2keV`` : double array(...)
;       The 2 keV monochromatic luminosity :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   
; Notes
; -----
;   The relation quoted above is from Kubota & Done (2018) and differs
;   slightly from that in Lusso & Risaliti (2017), though it is still
;   attributed to that paper. We also assume 
;   :math:`v_{FWHM} = 2000\ {\rm km\ s}^{-1}` as measured from the Mg II line.
;   
; References
; ----------
;   - `Lusso, E., & Risaliti, G. 2017, A&A, 602, A79 <https://ui.adsabs.harvard.edu/abs/2017A%26A...602A..79L/abstract>`_
;   - `Kubota, A., & Done, C. 2018, MNRAS, 480, 1247 <https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.1247K/abstract>`_
;   
; Modification History
; --------------------
;   - 2021/07/06: Created (Erik B. Monson)
;   - 2022/03/16: Moved to separate file and documentation improved (Erik B. Monson).
;   - 2022/06/20: Replaced ``!cv`` with ``!lightning_cgs`` (Keith Doore)
;   - 2022/06/20: Updated documentation (Keith Doore)
;   - 2022/06/20: Added error handling (Keith Doore)
;   - 2022/06/20: Added ``error_check`` keyword to do error handling (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(L2500) eq 0 then message, 'Variable is undefined: L2500.'
   if size(L2500, /type) lt 2 or size(L2500, /type) gt 5 then $
     message, 'L2500 must be of type int, float, or double.'
   if min(L2500) lt 0 then message, 'L2500 must only contain non-negative values.'
 endif

 ; We convert to erg s-1 Hz-1
 log_L2keV = 0.633 * alog10(L2500 * !lightning_cgs.Lsun) + 7.216

 ; And convert back to Lsun Hz-1
 return, (10.d0^log_L2keV) / !lightning_cgs.Lsun

end