function binned_stellar_spectrum, stellar_models, psi, atten_curve=atten_curve, tauV_diff=tauV_diff, $
                                  delta=delta, tauV_BC=tauV_BC, tauB_f=tauB_f, F_clump=F_clump, cosi=cosi, $
                                  b_to_d=b_to_d, rold0_ages=rold0_ages, atten_models=atten_models, $
                                  uv_bump=uv_bump, error_check=error_check, steps_Lnu_spec=steps_Lnu_spec, $
                                  Lnu_spec_unred_stellar=Lnu_spec_unred_stellar, Lbol_abs_stellar=Lbol_abs_stellar, $
                                  steps_Lbol_abs=steps_Lbol_abs
;+
; Name
; ----
;   BINNED_STELLAR_SPECTRUM
;
; Purpose
; -------
;   Generates the attenuated stellar model spectra for the non-parametric SFH.
;   The output model spectra depend on the age bins and attenuation.
;
; Calling Sequence
; ----------------
;   ::
;
;       Lnu_spec_stellar = binned_stellar_spectrum(stellar_models, psi [, atten_curve = , $
;                                       tauV_diff = , delta = , tauV_BC = , tauB_f = , F_clump = , $
;                                       cosi = , b_to_d = , rold0_ages = , atten_models = , /uv_bump, $
;                                       /error_check, steps_Lnu_spec=steps_Lnu_spec, $
;                                       Lnu_spec_unred_stellar=Lnu_spec_unred_stellar,$
;                                       Lbol_abs_stellar=Lbol_abs_stellar, steps_Lbol_abs=steps_Lbol_abs])
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
; Optional Inputs
; ---------------
;   ``atten_curve`` : string scalar
;       The name of the attenuation curve to apply to the stellar models. Current
;       options are ``'CALZETTI00'``, ``'CALZETTI_MOD'``, and ``'DOORE21'``.
;       (Default = ``'CALZETTI00'``)
;   ``tauV_diff`` : int, float, or double array(Nmodels)
;       The V-band optical depth of the diffuse dust for the Calzetti attenuation.
;       (Default = ``1.0``)
;   ``delta`` : int, float, or double array(Nmodels)
;       The power law value to change the attenuation curve slope for the Calzetti
;       attenuation. (Default = ``0.d0``)
;   ``tauV_BC`` : int, float, or double array(Nmodels)
;       The V-band optical depth of the birth cloud for the Calzetti attenuation.
;       (Default = ``0.0``)
;   ``tauB_f`` : int, float, or double array(Nmodels)
;       The face-on optical depth in the B-band for the Doore attenuation.
;       (Default = ``1.0``)
;   ``F_clump`` : int, float, or double array(Nmodels)
;       The clumpiness factor F for the Doore attenuation. (Default = ``0.0``)
;   ``cosi`` : int, float, or double array(Nmodels)
;       The inclination of the galactic disk in terms of cos(i) for the Doore
;       attenuation. (Default = ``1.d0``)
;   ``b_to_d`` : int, float, or double array(Nmodels)
;       The bulge-to-disk ratio for the Doore attenuation. (Default = ``0.0``)
;   ``rold0_ages`` : int, float, or double array(Nsteps)
;       The binary parameter ``rold0``, designating each SFH age bin as part of
;       the young or old population when using the Doore attenuation. A value
;       of ``0`` for the corresponding age bin considers it to be part of the young
;       population, and a value of ``1`` considers it to be part of the old
;       populations. (Default = Ages < 500 Myr are ``0`` else ``1``)
;   ``atten_models`` : structure
;       A structure containing the preloaded files for the Doore attenuation.
;   ``uv_bump`` : flag
;       If set, then a 2175 Angstrom UV bump feature will be added to the 
;       attenuation curve.
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``Lnu_spec_stellar`` : double array(Nwave, Nmodels)
;       The luminosity spectrum for each set of stellar model parameters at
;       each wavelength :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;
; Optional Outputs
; ----------------
;   ``steps_Lnu_spec`` : double array(Nwave, Nsteps, Nmodels)
;       The luminosity spectrum for each SFH age bin and set of stellar model 
;       parameters at each wavelength :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   ``Lnu_spec_unred_stellar`` : double array(Nwave, Nmodels)
;       The unattenuated luminosity spectrum for each set of stellar model
;       parameters at each wavelength :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   ``Lbol_abs_stellar`` : double array(Nmodels)
;       The bolometric luminosity of the attenuated stellar emission for each
;       set of stellar model parameters :math:`[L_\odot]`.
;   ``steps_Lbol_abs`` : double array(Nsteps, Nmodels)
;       The bolometric luminosity of the attenuated stellar emission for
;       each SFH age bin and set of stellar model parameters :math:`[L_\odot]`.
;
; Note
; ----
;   Defaults will only be set if the optional ``error_check`` input is set.
;
; Modification History
; --------------------
;   - 2022/04/13: Created (Keith Doore)
;   - 2022/06/30: Updated to match format of ``binned_stellar_sed.pro`` (Keith Doore)
;   - 2022/07/07: Change name of ``sfh_coeff`` to ``psi`` (Keith Doore)
;   - 2022/07/22: Added optional output of unreddened Lnu spectrum (Keith Doore)
;   - 2022/11/15: Fixed bug with variable redshift and ``rold0_ages`` (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(stellar_models) eq 0 then message, 'Variable is undefined: STELLAR_MODELS.'
   if size(stellar_models, /type) ne 8 then message, 'STELLAR_MODELS is not of type structure.'

   if n_elements(psi) eq 0 then message, 'Variable is undefined: PSI.'
   if size(psi, /type) lt 2 or size(psi, /type) gt 5 then $
     message, 'PSI must be of type int, float, or double.'
   if size(psi, /n_dim) gt 2 then $
     message, 'PSI must be a scalar, 1-D, or 2-D array.'
   Nsteps = (size(psi, /dim))[0]
   if size(psi, /n_dim) eq 2 then Nmodels = (size(psi, /dim))[1] else Nmodels = 1

   if n_elements(atten_curve) ne 0 then begin
     if size(atten_curve, /type) ne 7 then $
              message, 'ATTEN_CURVE must be of type string.'
     if size(atten_curve, /n_dim) ne 0 then message, 'ATTEN_CURVE must be a scalar.'
     if total(strupcase(atten_curve) eq ['CALZETTI00', 'CALZETTI_MOD', 'DOORE21']) ne 1 then $
       message, "ATTEN_CURVE must be set to either 'CALZETTI00', 'CALZETTI_MOD', or 'DOORE21'."
   endif else atten_curve = 'CALZETTI00'

   case strupcase(atten_curve) of
     'CALZETTI00': begin
          if n_elements(tauV_diff) ne 0 then begin
            if size(tauV_diff, /type) lt 2 or size(tauV_diff,/type) gt 5 then $
                     message, 'TAUV_DIFF must be of type int, float, or double.'
            if size(reform(tauV_diff), /n_dim) ne 1 then message, 'TAUV_DIFF must be a scalar or 1-D array.'
            if min(tauV_diff) lt 0 then message, 'TAUV_DIFF must only contain non-negative values.'
          endif else tauV_diff = replicate(1.0d, Nmodels)
  
          if n_elements(tauV_diff) ne Nmodels then $
            message, 'TAUV_DIFF must contain the same number of elements as the second dimension of PSI.'
        end

     'CALZETTI_MOD': begin
          if n_elements(tauV_diff) ne 0 then begin
            if size(tauV_diff, /type) lt 2 or size(tauV_diff,/type) gt 5 then $
                     message, 'TAUV_DIFF must be of type int, float, or double.'
            if size(reform(tauV_diff), /n_dim) ne 1 then message, 'TAUV_DIFF must be a scalar or 1-D array.'
            if min(tauV_diff) lt 0 then message, 'TAUV_DIFF must only contain non-negative values.'
          endif else tauV_diff = replicate(1.0d, Nmodels)

          if n_elements(delta) ne 0 then begin
            if size(delta, /type) lt 2 or size(delta,/type) gt 5 then $
                     message, 'DELTA must be of type int, float, or double.'
            if size(reform(delta), /n_dim) ne 1 then message, 'DELTA must be a scalar or 1-D array.'
          endif else delta = replicate(0.0d, Nmodels)

          if n_elements(tauV_BC) ne 0 then begin
            if size(tauV_BC, /type) lt 2 or size(tauV_BC,/type) gt 5 then $
                     message, 'TAUV_BC must be of type int, float, or double.'
            if size(reform(tauV_BC), /n_dim) ne 1 then message, 'TAUV_BC must be a scalar or 1-D array.'
            if min(tauV_BC) lt 0 then message, 'TAUV_BC must only contain non-negative values.'
          endif else tauV_BC = replicate(0.0d, Nmodels)

          if n_elements(tauV_diff) ne Nmodels or n_elements(delta) ne Nmodels or n_elements(tauV_BC) ne Nmodels then $
            message, 'TAUV_DIFF, DELTA, and TAUV_BC must contain the same number of elements as the second dimension of PSI.'
        end

     'DOORE21': begin
          if n_elements(tauB_f) ne 0 then begin
            if size(tauB_f, /type) lt 2 or size(tauB_f,/type) gt 5 then $
                     message, 'TAUB_F must be of type int, float, or double.'
            if size(reform(tauB_f), /n_dim) ne 1 then message, 'TAUB_F must be a scalar or 1-D array.'
            if min(tauB_f) lt 0 or max(tauB_f) gt 8 then message, 'TAUB_F must contain values between 0 and 8.'
          endif else tauB_f = replicate(1.0d, Nmodels)

          if n_elements(F_clump) ne 0 then begin
            if size(F_clump, /type) lt 2 or size(F_clump,/type) gt 5 then $
                     message, 'F_CLUMP must be of type int, float, or double.'
            if size(reform(F_clump), /n_dim) ne 1 then message, 'F_CLUMP must be a scalar or 1-D array.'
            if min(F_clump) lt 0 or max(F_clump) gt 0.61 then message, 'F_CLUMP must contain values between 0 and 0.61.'
          endif else F_clump = replicate(0.0d, Nmodels)

          if n_elements(cosi) ne 0 then begin
            if size(cosi, /type) lt 2 or size(cosi,/type) gt 5 then $
                     message, 'COSI must be of type int, float, or double.'
            if size(reform(cosi), /n_dim) ne 1 then message, 'COSI must be a scalar or 1-D array.'
            if min(cosi) lt 0 or max(cosi) gt 1 then message, 'COSI must contain values between 0 and 1.'
          endif else cosi = replicate(1.0d, Nmodels)

          if n_elements(b_to_d) ne 0 then begin
            if size(b_to_d, /type) lt 2 or size(b_to_d,/type) gt 5 then $
                     message, 'B_TO_D must be of type int, float, or double.'
            if size(reform(b_to_d), /n_dim) ne 1 then message, 'B_TO_D must be a scalar or 1-D array.'
            if min(b_to_d) lt 0 then message, 'B_TO_D must only contain non-negative values.'
          endif else b_to_d = replicate(0.0d, Nmodels)

          if n_elements(tauB_f) ne Nmodels or n_elements(F_clump) ne Nmodels or $
             n_elements(cosi) ne Nmodels   or n_elements(b_to_d) ne Nmodels then $
            message, 'TAUB_F, F_CLUMP, COSI, and B_TO_D must contain the same number of elements as the second dimension of PSI.'

          if n_elements(rold0_ages) ne 0 then begin
            if size(rold0_ages, /type) lt 2 or size(rold0_ages,/type) gt 5 then $
                     message, 'ROLD0_AGES must be of type int, float, or double.'
            if size(rold0_ages, /n_dim) gt 1 then message, 'ROLD0_AGES must be a scalar or 1-D array.'
            if n_elements(rold0_ages) ne Nsteps then $
              message, 'ROLD0_AGES must contain the same number of elements as the first dimension of PSI.'
            if total(rold0_ages eq 0 or rold0_ages eq 1) ne Nsteps then message, 'ROLD0_AGES must only contain values of 0 or 1.'
          endif else begin
            rold0_ages = dblarr(Nsteps)
            rold0_ages[where((stellar_models.bounds)[1:*] gt 5e8, /null)]= 1
          endelse

          if n_elements(atten_models) ne 0 then begin
            if size(atten_models, /type) ne 8 then message, 'ATTEN_MODELS is not of type structure.'
            if total(tag_names(atten_models) eq 'DOORE21_LBOL_ABS_TABLE') ne 1 then $
              message, 'ATTEN_MODELS must contain DOORE21_LBOL_ABS_TABLE tag if using the Doore21 attenuation curves.'
          endif
        end
   endcase  
 endif


 Nsteps = n_elements(stellar_models.bounds) - 1
 wave_rest = stellar_models.wave_rest ; restframe wavelength
 wave_obs  = stellar_models.wave_obs  ; observed wavelength
 Nwave = n_elements(wave_rest)
 nu_rest = 1.d4 * !lightning_cgs.clight / wave_rest
 nu_obs  = 1.d4 * !lightning_cgs.clight / wave_obs
 redshift = stellar_models.redshift


; Generate the attenuation in terms of exp(-tau)
 case strupcase(atten_curve) of
   'CALZETTI00': begin
        Nmodels = n_elements(tauV_diff)
        exp_neg_tau    = calzetti00_atten(wave_rest,              $
                                          tauV_diff=tauV_diff,    $
                                          delta=dblarr(Nmodels),   $
                                          tauV_BC=dblarr(Nmodels))
      end

   'CALZETTI_MOD': begin
        Nmodels = n_elements(tauV_diff)
        exp_neg_tau_BC = calzetti00_atten(wave_rest,              $
                                          tauV_diff=tauV_diff,    $
                                          delta=delta,            $
                                          tauV_BC=tauV_BC,        $
                                          uv_bump=uv_bump)
        exp_neg_tau    = calzetti00_atten(wave_rest,              $
                                          tauV_diff=tauV_diff,    $
                                          delta=delta,            $
                                          tauV_BC=dblarr(Nmodels), $
                                          uv_bump=uv_bump)
      end

   'DOORE21': begin
        Nmodels = n_elements(tauB_f)
        rold0_y = dblarr(Nmodels)
        exp_neg_tau_young = doore21_atten(wave_rest,       $
                                          tauB_f=tauB_f,   $
                                          rold0=rold0_y,   $
                                          b_to_d=b_to_d,   $
                                          F_clump=F_clump, $
                                          cosi=cosi,       $
                                          _extra=atten_models)
        rold0_o = replicate(1.d0, Nmodels)
        exp_neg_tau_old   = doore21_atten(wave_rest,       $
                                          tauB_f=tauB_f,   $
                                          rold0=rold0_o,   $
                                          b_to_d=b_to_d,   $
                                          F_clump=F_clump, $
                                          cosi=cosi,       $
                                          _extra=atten_models)
      end
 endcase

 steps_Lnu_red = dblarr(Nwave, Nsteps, Nmodels)
 steps_Lbol_abs = dblarr(Nsteps, Nmodels)

; Calculate the mean Lnu of each filter for each model and bolometric luminosity of absorbed emission
 for st=0,(Nsteps-1) do begin
   steps_Lnu_unred = rebin(stellar_models.lnu[*, st], Nwave, Nmodels)

   ; Attenuate stellar models
   case strupcase(atten_curve) of
     'CALZETTI00': steps_Lnu_red[*, st, *] = exp_neg_tau * steps_Lnu_unred

     'CALZETTI_MOD': begin
          if st eq 0 then steps_Lnu_red[*, st, *] = exp_neg_tau_BC * steps_Lnu_unred else $
                          steps_Lnu_red[*, st, *] = exp_neg_tau * steps_Lnu_unred
        end

     'DOORE21': begin
          if rold0_ages[st] eq 0 then steps_Lnu_red[*, st, *] = exp_neg_tau_young * steps_Lnu_unred
          if rold0_ages[st] eq 1 then steps_Lnu_red[*, st, *] = exp_neg_tau_old * steps_Lnu_unred
        end
   endcase

   ;Lbol = -1.d0*trapez(steps_Lnu_unred, nu_obs)   ; does not include nebular emission absorption
   Lbol = (stellar_models.Lbol[st])[0]
   for nm=0,(Nmodels-1) do begin
     if strupcase(atten_curve) ne 'DOORE21' then begin
       ; gives negative value from integral bc nu_obs is reversed (big to small)
       steps_Lbol_abs[st, nm] = Lbol + trapez(steps_Lnu_red[*, st, nm], nu_obs)
     endif
   endfor
 endfor

 if strupcase(atten_curve) eq 'DOORE21' then begin
   steps_Lbol_abs = doore21_interp_lbol_abs_table(tauB_f, F_clump, b_to_d, rold0_ages[0:(Nsteps-1)], atten_models.doore21_Lbol_abs_table)
 endif

 ; Exclude NaNs in case a filter included NaNs.
 steps_Lnu_spec = temporary(steps_Lnu_red)
 Lnu_stellar_spec = total(steps_Lnu_spec * rebin(reform(psi, 1, Nsteps, Nmodels), Nwave, Nsteps, Nmodels), 2, /nan)
 Lnu_spec_unred_stellar = total(rebin(reform(stellar_models.Lnu, Nwave, Nsteps, 1), Nwave, Nsteps, Nmodels) * $
                                rebin(reform(psi, 1, Nsteps, Nmodels), Nwave, Nsteps, Nmodels), 2, /NaN)
 Lbol_abs_stellar = total(steps_Lbol_abs * psi, 1, /nan)

 return, Lnu_stellar_spec

end
