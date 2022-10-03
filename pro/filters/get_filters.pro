pro get_filters, filter_labels, nu, filters, Lnu=Lnu, error_check=error_check, $
                 plot_filters=plot_filters, mean_nu=mean_nu, sigma_nu=sigma_nu, $
                 mean_Lnu=mean_Lnu, wave=wave, Llambda=Llambda, mean_wave=mean_wave, $
                 sigma_wave=sigma_wave, mean_Llambda=mean_Llambda
                 
;+
; Name
; ----
;   GET_FILTERS
;
; Purpose
; -------
;   Generates a filter transmission array for the input filter labels.
;   The specified filters are read in from each tabulated filter profile
;   and interpolated onto the input frequency array.
;
; Calling Sequence
; ----------------
;   ::
;
;       get_filters, filter_labels, nu, filters [, Lnu = ,   $
;                    /error_check, /plot_filters, mean_nu=mean_nu, $
;                    sigma_nu=sigma_nu, mean_Lnu=mean_Lnu, $
;                    wave=wave, Llambda=Llambda, mean_wave=mean_wave, $
;                    sigma_wave=sigma_wave, mean_Llambda=mean_Llambda]
;
; Inputs
; ------
;   ``filter_labels`` : string array(Nfilters)
;       The names of the filters.
;   ``nu`` : int, float, or double array(Nnu)
;       The frequencies for interpolation :math:`[{\rm Hz}]`.
;
; Optional Inputs
; ---------------
;   ``Lnu`` : int, float, or double array(NLnu, Nnu)
;       The luminosity densities at each element of ``nu`` to be used in
;       the calculation of the optional outputs. The luminosity units can
;       be user chosen and will dictate the optional output luminosity units.
;       However, overall unit must be per Hertz (i.e., :math:`[L\ {\rm Hz}^{-1}]`).
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;   ``plot_filters`` : flag
;       If set, a plot will be generated showing the normalized filter
;       transmissions as a function of wavelength.
;
; Output
; ------
;   ``filters`` : double array(Nfilters, Nnu)
;       The tabulated filter transmission functions interpolated and normalized
;       onto the input frequency array.
;
; Optional Outputs
; ----------------
;   ``mean_nu`` : double array(NLnu, Nfilters)
;       The mean frequency of each filter as computed from each set of ``Lnu``
;       :math:`[\rm Hz]`.
;   ``sigma_nu`` : double array(NLnu, Nfilters)
;       The standard deviation width of each filter as computed
;       from each set of ``Lnu`` :math:`[\rm Hz]`.
;   ``mean_Lnu`` : double array(NLnu, Nfilters)
;       The mean luminosity density of each filter in the same units as
;       ``Lnu`` (:math:`[L\ {\rm Hz}^{-1}]`) as computed from each set of ``Lnu``.
;   ``wave`` : double array(Nnu)
;       The wavelengths converted from ``nu`` :math:`[\mu \rm m]`.
;   ``Llambda`` : double array(NLnu, Nnu)
;       The luminosity density converted from ``Lnu`` in the same luminosity
;       units as ``Lnu`` per micron (i.e. :math:`[L\ {\mu \rm m}^{-1}]`).
;   ``mean_wave`` : double array(NLnu, Nfilters)
;       The mean wavelength of each filter as computed from each set of ``Lnu``
;       :math:`[\mu \rm m]`.
;   ``sigma_wave`` : double array(NLnu, Nfilters)
;       The standard deviation width of each filter as computed from each
;       set of ``Lnu`` :math:`[\mu \rm m]`.
;   ``mean_Llambda`` : double array(NLnu, Nfilters)
;       The mean luminosity density of each filter in the same units as ``Llambda``
;       (:math:`[L\ {\mu \rm m}^{-1}]`) as computed from each set of ``Llambda``.
;
; Notes
; -----
;   The optional outputs will only be computed if the optional input ``Lnu`` is
;   specified.
;
; Modification History
; --------------------
;   - 2016/05/01: Created (Rafael T. Eufrasio)
;   - 2019/04/28: Allowed for tabulated filter profile to exceed input ``nu`` range (Rafael T. Eufrasio)
;   - 2022/03/02: Added error handling (Keith Doore)
;   - 2022/03/02: Changed list of filters and directories from case statement to json file (Keith Doore)
;   - 2022/03/02: Minor bug fixes (Keith Doore)
;   - 2022/03/18: Updated documentation (Keith Doore)
;   - 2022/03/22: Updated to used plot function vs plot procedure (Keith Doore)
;   - 2022/04/12: Allowed for ``nu`` to have degenerate dimensions (Keith Doore)
;   - 2022/04/12: Allowed for ``filter_labels`` to be scalars (Keith Doore)
;   - 2022/04/12: Allowed integer inputs (Keith Doore)
;   - 2022/04/12: Made the output ``filters`` frequency normalized (Keith Doore)
;   - 2022/05/16: Removed ``filters_dir`` and changed to use ``!lightning_dir`` system variable (Keith Doore)
;   - 2022/05/16: Changed ``!cv`` constants call to ``!lightning_cgs`` call (Keith Doore)
;   - 2022/05/17: Added ``error_check`` keyword to do error handling (Keith Doore)
;   - 2022/06/08: Added ``strtrim`` to ``filter_labels`` to remove any extra blank padding (Keith Doore)
;   - 2022/06/29: Updated documentation and changed ``L_wave`` to ``Llambda`` (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(filter_labels) eq 0 then message, 'Variable is undefined: FILTER_LABELS.'
   if size(filter_labels, /type) ne 7 then message, 'FILTER_LABELS must be of type string.'
   if size(reform(filter_labels), /n_dim) ne 1 then message, 'FILTER_LABELS must be a scalar or 1-D array.'
  
   if n_elements(nu) eq 0 then message, 'Variable is undefined: NU.'
   if size(nu, /type) lt 2 or size(nu,/type) gt 5 then message, 'NU must be of type int, float, or double.'
   if size(reform(nu), /n_dim) ne 1 then message, 'NU must be a 1-D array.'
   if n_elements(nu) eq 1 then message, 'NU must have more than one element.'
   if min(nu) le 0 then message, 'NU must only contain positive values.'
  
   if n_elements(Lnu) ne 0 then begin
     if size(Lnu, /type) lt 2 or size(Lnu, /type) gt 5 then message, 'LNU must be of type int, float, or double.'
     if size(Lnu, /n_dim) ne 2 then message, 'LNU must be a 2-D array.'
     if (size(lnu, /dim))[1] ne n_elements(nu) then message, 'The second dimension of LNU must have the '+$
                                                             'same number of elements as NU.'
   endif
 endif
 
 if n_elements(Lnu) eq 0 then analysis = 0 else analysis = 1


; Read in Filter files and interpolate
 Nfilters = n_elements(filter_labels)
 Nwave = n_elements(nu)

 Filters = dblarr(Nfilters, Nwave)

 wave = 1.d4*!lightning_cgs.clight/reform(nu)

 if analysis then begin
   Llambda = 1.d4*!lightning_cgs.clight*Lnu/(rebin(reform(wave, 1,n_elements(wave)), size(Lnu, /dim)))^2
   Nlnu = (size(Lnu, /dim))[0]
   mean_nu = dblarr(Nlnu, Nfilters)
   mean_Lnu = dblarr(Nlnu, Nfilters)
   mean_nu2 = dblarr(Nlnu, Nfilters)
   mean_Llambda = dblarr(Nlnu, Nfilters)
   mean_wave = dblarr(Nlnu, Nfilters)
   mean_wave2 = dblarr(Nlnu, Nfilters)
 endif

 json = JSON_PARSE(!lightning_dir+'filters/filters.json')

 for f=0,(Nfilters-1) do begin
   trans_curve = !lightning_dir+'filters/'+json[strtrim(filter_labels[f], 2)]

   ; Read wavelength in microns and normalized trasmission
   readcol, trans_curve, wraw, Kraw, form='D,D', /silent

   ; Combine read wavelength and input wavelength vector, keeping only unique values
   wcomb = [wraw, wave]
   wcomb = wcomb[uniq(wcomb, sort(wcomb))]


   ; Check where input and combined wavelengths exceed read wavelength range
   offwraw       = where((wave  lt min(wraw)) or (wave  gt max(wraw)), /null, comp=inwraw)
   wcomb_offwraw = where((wcomb lt min(wraw)) or (wcomb gt max(wraw)), /null, comp=wcomb_inwraw)

   ; Interpolate normalized trasmission to combined wavelength grid and normalize area
   Kcomb = dblarr(n_elements(wcomb))
   Kcomb[wcomb_inwraw] = interpol(Kraw, wraw, wcomb[wcomb_inwraw])
   Kcomb /= trap_int(wcomb, Kcomb)

   ; Convolution is only possible if read wavelength is fully included within input wavelengths
   ; If not possible, set the filter interpolation to NaN
   possible = ((min(wraw) ge min(wave)) and (max(wraw) le max(wave)))
   if possible then begin

     ; If input wavelength grid is too coarse, the full range of the read wavelength may be
     ;   between two consecutive points of the input wavelength grid.
     ; If this happens, for instance with PEGASE2 models convolved with FIR filters,
     ;   set the filter interpolation to NaN.
     if n_elements(inwraw) ne 0 then begin
       wcomb_inwave = dblarr(Nwave)
       for kk=0,(Nwave-1) do wcomb_inwave[kk] = where(wcomb eq wave[kk])
       Kwave = Kcomb[wcomb_inwave]

       ; If all elements of Filter are zero, set the filter interpolation to NaN.
       if max(Kwave) eq 0 then begin
         Filters[f,*] = !values.D_NaN
       endif else begin
         ; Normalize to frequency grid (absolute value for potentially reversed nu)
         Kwave /= abs(trap_int(nu, Kwave))
         Filters[f,*] = (Kwave)
       endelse

     endif else begin
       Filters[f,*] = !values.D_NaN
     endelse

   endif else begin
     Filters[f,*] = !values.D_NaN
   endelse

   if analysis then begin
     ; For each set of LNU, generate the corresponding set of optional outputs
     for t=0,(Nlnu-1) do begin
       nu_wcomb = 1.d4*!lightning_cgs.clight/wcomb
       Lnu_wcomb = interpol(reform(Lnu[t,*]), wave, wcomb)
       filter_Lnu_wcomb = Kcomb*Lnu_wcomb
       fnorm = trap_int(nu_wcomb, Kcomb)
       mean_Lnu[t,f] = trap_int(nu_wcomb, filter_Lnu_wcomb)/temporary(fnorm)
       fnorm = trap_int(nu_wcomb, filter_Lnu_wcomb)
       mean_nu[t,f]  = trap_int(nu_wcomb, nu_wcomb  *filter_Lnu_wcomb)/fnorm
       mean_nu2[t,f] = trap_int(nu_wcomb, nu_wcomb^2*filter_Lnu_wcomb)/temporary(fnorm)

       Llambda_wcomb = 1.d4*!lightning_cgs.clight*Lnu_wcomb/wcomb^2
       filter_Llambda_wcomb = Kcomb*Llambda_wcomb
       fnorm = trap_int(wcomb, Kcomb)
       mean_Llambda[t,f]  = trap_int(wcomb, filter_Llambda_wcomb)/temporary(fnorm)
       fnorm = trap_int(wcomb, filter_Llambda_wcomb)
       mean_wave [t,f] = trap_int(wcomb, wcomb  *filter_Llambda_wcomb)/          fnorm
       mean_wave2[t,f] = trap_int(wcomb, wcomb^2*filter_Llambda_wcomb)/temporary(fnorm)
     endfor
   endif
 endfor

 if analysis then begin
   sigma_nu   = sqrt(mean_nu2 - mean_nu^2)
   sigma_wave = sqrt(mean_wave2 - mean_wave^2)
 endif

 if keyword_set(plot_filters) then begin
   ; Plot normalized Filters
   max_filt=max(Filters, dim=2)
   x=plot(wave, Filters[0,*]/max_filt[0], xrange=minmax(wave), yrange=[0,1.1], /xlog, $
          xtitle='Wavelength [$\mu$m]', ytitle='Unity Normalized Transmission', thick=2, $
          vert_color=BYTSCL(alog10(wave[where(Filters[0,*] eq max_filt[0])]), $
          min=alog10(min(wave)), max=alog10(max(wave))), rgb_table=34)
   for i=0,(Nfilters-1) do x=plot(wave, Filters[i,*]/max(Filters[i,*]), /over, rgb_table=34, thick=2, $
           vert_color=BYTSCL(alog10(wave[where(Filters[i,*] eq max_filt[i])]), min=alog10(min(wave)), $
           max=alog10(max(wave))))
 endif

 return

end
