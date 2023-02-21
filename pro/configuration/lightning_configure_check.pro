function lightning_configure_check, config
;+
; Name
; ----
;   LIGHTNING_CONFIGURE_CHECK
;
; Purpose
; -------
;   Checks the Lightning configuration structure to make sure
;   all inputs are correctly formatted and can be understood by
;   Lightning. All unused parameters are removed from the
;   structure.
;
; Calling Sequence
; ----------------
;   ::
;
;       config_checked = lightning_configure_check(config)
;
; Input
; -----
;   ``config`` : structure
;       A Lightning configuration structure. (See
;       ``lightning_configure_defaults.pro`` for details and contents.)
;
; Output
; ------
;   ``config_checked`` : structure
;       The error checked Lightning configuration structure with
;       unused tags removed.
;
; Modification History
; --------------------
;   - 2022/04/26: Created (Keith Doore)
;   - 2022/05/03: Updated method for removing key from hash (Keith Doore)
;   - 2022/05/17: Added ability to check chosen cosmology (Keith Doore)
;   - 2022/05/20: Removed loguniform and lognormal prior options (Keith Doore)
;   - 2022/06/16: Added check for number of parallel chains if using affine MCMC (Keith Doore)
;   - 2022/06/29: Added conversion of ``nebular_extinction`` to ``no_nebular_extinction`` (Keith Doore)
;   - 2022/06/29: Added conversion of ``emission_lines`` to ``no_emission_lines`` (Keith Doore)
;   - 2022/07/08: Added check to make sure ``energy_balance`` is ``0`` if only using a dust model (Keith Doore)
;   - 2022/08/16: Added MPFIT keywords (Keith Doore)
;   - 2022/08/18: Added checks to ensure MCMC post-processing values are reasonable for number of trials (Keith Doore)
;   - 2022/09/01: Replaced ``XRAY_STAT`` with ``XRAY_UNC`` (Erik B. Monson)
;   - 2022/09/14: Updates to allow fitting with X-ray fluxes (Erik B. Monson)
;   - 2022/09/14: Fixed issue with MCMC post-processing size check for affine MCMC (Keith Doore)
;   - 2022/10/24: Added option to choose stranded walker deviation value for affine MCMC (Keith Doore)
;   - 2022/10/25: Renamed SPS to SSP (Keith Doore)
;   - 2022/12/13: Prevented ``XRAY_UNC`` from begin checked if ``XRAY_UNIT='FLUX'`` (Keith Doore)
;   - 2023/01/31: Added check of added ``OUTPUT_FILENAME`` option (Keith Doore)
;   - 2023/02/17: Added metallicity of Z=0.01 (Keith Doore)
;-
 On_error, 2
 Compile_opt idl2

; Confirm input is a structure
 if n_elements(config) eq 0 then message, 'Variable is undefined: CONFIG.'
 if size(config, /type) ne 8 then message, 'CONFIG must be a structure.'
 tags=tag_names(config)

; Create a hash from the configuration structure for parameter removal
 config_hash = orderedhash(config)

; Current options for Lightning parameters (makes for easier updating if more options are added).
 prior_options = ['fixed', 'uniform', 'normal', 'tabulated']
 prior_options_narg = [1, 2, 4, 1]
 ssp_options = ['PEGASE', 'none']
 imf_options = ['Kroupa01']
 zmetal_options = [0.001, 0.004, 0.008, 0.01, 0.02, 0.05, 0.1]
 sfh_options = ['NON-PARAMETRIC']
 atten_curve_options = ['CALZETTI00', 'CALZETTI_MOD', 'DOORE21']
 dust_emission_options = ['DL07', 'none']
 agn_emission_options = ['SKIRTOR', 'none']
 xray_unit_options = ['COUNTS', 'FLUX']
 xray_unc_options = ['SQRT', 'GEHRELS', 'USER']
 xray_abs_options = ['TBABS-WILM', 'ATTEN']
 xray_agn_options = ['QSOSED', 'PLAW', 'none']
 ;fit_algorithm_options = ['GRID', 'INVERSION', 'MCMC-ADAPTIVE', 'MCMC-AFFINE', 'MPFIT']
 fit_algorithm_options = ['MCMC-ADAPTIVE', 'MCMC-AFFINE', 'MPFIT']

; Check tags and print error message if error is found.
 base_err = 'Error found in the configuration structure. '


; ===================    Core Parameter Check    =============================================
 if n_elements(where(strupcase(tags) eq 'OUTPUT_FILENAME', /null)) eq 0 then $
   message, base_err+'OUTPUT_FILENAME tag is missing.'
 if size(config.OUTPUT_FILENAME, /type) ne 7 then $
   message, base_err+'OUTPUT_FILENAME tag must be of type string.'
 if size(config.OUTPUT_FILENAME, /n_dim) ne 0 then $
   message, base_err+'OUTPUT_FILENAME tag must be a scalar.'

 if n_elements(where(strupcase(tags) eq 'PRINT_PROGRESS', /null)) eq 0 then $
   message, base_err+'PRINT_PROGRESS tag is missing.'
 if size(config.PRINT_PROGRESS, /type) lt 2 or size(config.PRINT_PROGRESS, /type) gt 3 then $
   message, base_err+'PRINT_PROGRESS tag must be of type int.'
 if size(config.PRINT_PROGRESS, /n_dim) ne 0 then $
   message, base_err+'PRINT_PROGRESS tag must be a scalar.'
 if ~(config.PRINT_PROGRESS eq 0 or config.PRINT_PROGRESS eq 1) then $
   message, base_err+'PRINT_PROGRESS tag must be a flag with a value of 0 or 1.'

 if n_elements(where(strupcase(tags) eq 'MAX_CPUS', /null)) eq 0 then $
   message, base_err+'MAX_CPUS tag is missing.'
 if size(config.MAX_CPUS, /type) lt 2 or size(config.MAX_CPUS, /type) gt 3 then $
   message, base_err+'MAX_CPUS tag must be of type int.'
 if size(config.MAX_CPUS, /n_dim) ne 0 then $
   message, base_err+'MAX_CPUS tag must be a scalar.'
 if config.MAX_CPUS le 0 then $
   message, base_err+'MAX_CPUS tag must be a positive value.'

 if n_elements(where(strupcase(tags) eq 'COSMOLOGY', /null)) eq 0 then $
   message, base_err+'COSMOLOGY tag is missing.'
 if size(config.COSMOLOGY, /type) ne 8 then $
   message, base_err+'COSMOLOGY tag must be of type structure.'
 cosmo_tags = strupcase(tag_names(config.COSMOLOGY))
 if total(cosmo_tags[sort(cosmo_tags)] eq ['H0','K','LAMBDA0','OMEGA_M','Q0']) ne 5 then $
   message, base_err+'H0, OMEGA_M, LAMBDA0, Q0, or K tags are missing from the COSMOLOGY structure.'
 if size(config.COSMOLOGY.H0, /type) lt 2 or size(config.COSMOLOGY.H0, /type) gt 5 then $
   message, base_err+'H0 tag in COSMOLOGY structure must be of type int, float, or double.'
 if size(config.COSMOLOGY.OMEGA_M, /type) lt 2 or size(config.COSMOLOGY.OMEGA_M, /type) gt 5 then $
   message, base_err+'OMEGA_M tag in COSMOLOGY structure must be of type int, float, or double.'
 if size(config.COSMOLOGY.LAMBDA0, /type) lt 2 or size(config.COSMOLOGY.LAMBDA0, /type) gt 5 then $
   message, base_err+'LAMBDA0 tag in COSMOLOGY structure must be of type int, float, or double.'
 if size(config.COSMOLOGY.Q0, /type) lt 2 or size(config.COSMOLOGY.Q0, /type) gt 5 then $
   message, base_err+'Q0 tag in COSMOLOGY structure must be of type int, float, or double.'
 if size(config.COSMOLOGY.K, /type) lt 2 or size(config.COSMOLOGY.K, /type) gt 5 then $
   message, base_err+'K tag in COSMOLOGY structure must be of type int, float, or double.'
 if size(config.COSMOLOGY.H0, /n_dim) ne 0 then $
   message, base_err+'H0 tag in COSMOLOGY structure must be a scalar.'
 if size(config.COSMOLOGY.OMEGA_M, /n_dim) ne 0 then $
   message, base_err+'OMEGA_M tag in COSMOLOGY structure must be a scalar.'
 if size(config.COSMOLOGY.LAMBDA0, /n_dim) ne 0 then $
   message, base_err+'LAMBDA0 tag in COSMOLOGY structure must be a scalar.'
 if size(config.COSMOLOGY.Q0, /n_dim) ne 0 then $
   message, base_err+'Q0 tag in COSMOLOGY structure must be a scalar.'
 if size(config.COSMOLOGY.K, /n_dim) ne 0 then $
   message, base_err+'K tag in COSMOLOGY structure must be a scalar.'
 if config.COSMOLOGY.H0 le 0 then $
   message, base_err+'H0 tag in COSMOLOGY structure must be a positive value.'
 if config.COSMOLOGY.OMEGA_M lt 0 then $
   message, base_err+'OMEGA_M tag in COSMOLOGY structure must be a non-negative value.'

 if n_elements(where(strupcase(tags) eq 'ENERGY_BALANCE', /null)) eq 0 then $
   message, base_err+'ENERGY_BALANCE tag is missing.'
 if size(config.ENERGY_BALANCE, /type) lt 2 or size(config.ENERGY_BALANCE, /type) gt 3 then $
   message, base_err+'ENERGY_BALANCE tag must be of type int.'
 if size(config.ENERGY_BALANCE, /n_dim) ne 0 then $
   message, base_err+'ENERGY_BALANCE tag must be a scalar.'
 if ~(config.ENERGY_BALANCE eq 0 or config.ENERGY_BALANCE eq 1) then $
   message, base_err+'ENERGY_BALANCE tag must be a flag with a value of 0 or 1.'

 if n_elements(where(strupcase(tags) eq 'MODEL_UNC', /null)) eq 0 then $
   message, base_err+'MODEL_UNC tag is missing.'
 if size(config.MODEL_UNC, /type) lt 2 or size(config.MODEL_UNC, /type) gt 5 then $
   message, base_err+'MODEL_UNC tag must be of type int, float, or double.'
 if size(config.MODEL_UNC, /n_dim) ne 0 then $
   message, base_err+'MODEL_UNC tag must be a scalar.'
 if config.MODEL_UNC lt 0 then $
   message, base_err+'MODEL_UNC tag must be a be a non-negative value.'




; ===================    Stellar Emission Parameter Check   =============================================
 if n_elements(where(strupcase(tags) eq 'SSP', /null)) eq 0 then $
   message, base_err+'SSP tag is missing.'
 if size(config.SSP, /type) ne 7 then $
   message, base_err+'SSP tag must be of type string.'
 if size(config.SSP, /n_dim) ne 0 then $
   message, base_err+'SSP tag must be a scalar.'
 if total(strupcase(config.SSP) eq strupcase(ssp_options)) eq 0 then $
   message, base_err+"SSP tag must be set to one of the following: '"+$
            strjoin(ssp_options, "', '")+"'."

     case strupcase(config.SSP) of
     ; ===================    Stellar Population Parameter Check   =============================================
       'PEGASE': begin
          if n_elements(where(strupcase(tags) eq 'IMF', /null)) eq 0 then $
            message, base_err+'IMF tag is missing.'
          if size(config.IMF, /type) ne 7 then $
            message, base_err+'IMF tag must be of type string.'
          if size(config.IMF, /n_dim) ne 0 then $
            message, base_err+'IMF tag must be a scalar.'
          if total(strupcase(config.IMF) eq strupcase(imf_options)) eq 0 then $
            message, base_err+"IMF tag must be set to one of the following: '"+$
                     strjoin(imf_options, "', '")+"'."

          if n_elements(where(strupcase(tags) eq 'ZMETAL', /null)) eq 0 then $
            message, base_err+'ZMETAL tag is missing.'
          if size(config.ZMETAL, /type) lt 4 or size(config.ZMETAL, /type) gt 5 then $
            message, base_err+'ZMETAL tag must be of type float or double.'
          if size(config.ZMETAL, /n_dim) ne 0 then $
            message, base_err+'ZMETAL tag must be a scalar.'
          if total(float(config.ZMETAL) eq ZMETAL_options) eq 0 then $
            message, base_err+'ZMETAL tag must be set to one of the following values: '+$
                     strjoin(string(zmetal_options, f='(F0)'), ', ')+'.'

          if n_elements(where(strupcase(tags) eq 'EMISSION_LINES', /null)) eq 0 then $
            message, base_err+'EMISSION_LINES tag is missing.'
          if size(config.EMISSION_LINES, /type) lt 2 or size(config.EMISSION_LINES, /type) gt 3 then $
            message, base_err+'EMISSION_LINES tag must be of type int.'
          if size(config.EMISSION_LINES, /n_dim) ne 0 then $
            message, base_err+'EMISSION_LINES tag must be a scalar.'
          if ~(config.EMISSION_LINES eq 0 or config.EMISSION_LINES eq 1) then $
            message, base_err+'EMISSION_LINES tag must be a flag with a value of 0 or 1.'
          ; Convert to NO_EMISSION_LINES, such that in functions the default flag is to include them.
          config_hash['NO_EMISSION_LINES'] = ~config.EMISSION_LINES

          if n_elements(where(strupcase(tags) eq 'NEBULAR_EXTINCTION', /null)) eq 0 then $
            message, base_err+'NEBULAR_EXTINCTION tag is missing.'
          if size(config.NEBULAR_EXTINCTION, /type) lt 2 or size(config.NEBULAR_EXTINCTION, /type) gt 3 then $
            message, base_err+'NEBULAR_EXTINCTION tag must be of type int.'
          if size(config.NEBULAR_EXTINCTION, /n_dim) ne 0 then $
            message, base_err+'NEBULAR_EXTINCTION tag must be a scalar.'
          if ~(config.NEBULAR_EXTINCTION eq 0 or config.NEBULAR_EXTINCTION eq 1) then $
            message, base_err+'NEBULAR_EXTINCTION tag must be a flag with a value of 0 or 1.'
          ; Convert to NO_NEBULAR_EXTINCTION, such that in functions the default flag is to include them.
          config_hash['NO_NEBULAR_EXTINCTION'] = ~config.NEBULAR_EXTINCTION

          if n_elements(where(strupcase(tags) eq 'SFH', /null)) eq 0 then $
            message, base_err+'SFH tag is missing.'
          if size(config.SFH, /type) ne 7 then $
            message, base_err+'SFH tag must be of type string.'
          if size(config.SFH, /n_dim) ne 0 then $
            message, base_err+'SFH tag must be a scalar.'
          if total(strupcase(config.SFH) eq strupcase(sfh_options)) eq 0 then $
            message, base_err+"SFH tag must be set to one of the following: '"+$
                     strjoin(sfh_options, "', '")+"'."


          case strupcase(config.SFH) of
          ; ===================    Non-Parametric SFH Parameter Check   =============================================
            'NON-PARAMETRIC' : begin
                if n_elements(where(strupcase(tags) eq 'STEPS_BOUNDS', /null)) eq 0 then $
                  message, base_err+'STEPS_BOUNDS tag is missing.'
                if size(config.STEPS_BOUNDS, /type) lt 2 or size(config.STEPS_BOUNDS, /type) gt 5 then $
                  message, base_err+'STEPS_BOUNDS tag must be of type int, float, or double.'
                if size(config.STEPS_BOUNDS, /n_dim) ne 1 then $
                  message, base_err+'STEPS_BOUNDS tag must be a 1-D array.'
                if total(config.STEPS_BOUNDS lt 0) gt 0 then $
                  message, base_err+'STEPS_BOUNDS tag must only contain non-negative values.'
                if total(sort(config.STEPS_BOUNDS) ne indgen(n_elements(config.STEPS_BOUNDS))) ne 0 then $
                  message, base_err+'STEPS_BOUNDS tag must contain values in ascending order.'
                if n_elements(uniq(config.STEPS_BOUNDS, sort(config.STEPS_BOUNDS))) ne $
                   n_elements(config.STEPS_BOUNDS) then $
                  message, base_err+'STEPS_BOUNDS tag must only contain unique values.'
                Nsteps = n_elements(config.STEPS_BOUNDS) - 1

                if n_elements(where(strupcase(tags) eq 'DTIME_SF', /null)) eq 0 then $
                  message, base_err+'DTIME_SF tag is missing.'
                if size(config.DTIME_SF, /type) lt 2 or size(config.DTIME_SF, /type) gt 5 then $
                  message, base_err+'DTIME_SF tag must be of type int, float, or double.'
                if size(config.DTIME_SF, /n_dim) ne 0 then $
                  message, base_err+'DTIME_SF tag must be a scalar.'
                if config.DTIME_SF le 0 then $
                  message, base_err+'DTIME_SF tag must be a positive value.'

                if n_elements(where(strupcase(tags) eq 'PSI', /null)) eq 0 then $
                  message, base_err+'PSI tag is missing.'
                if size(config.PSI, /type) ne 8 then $
                  message, base_err+'PSI tag must be of type structure.'
                psi_range = [0.0d, !values.d_infinity]
                prior_check, config.PSI, 'PSI', Nsteps, psi_range, prior_options, prior_options_narg
              end
          endcase
        end

       'NONE': begin
          unused_tags = ['IMF', 'ZMETAL', 'EMISSION_LINES', 'NEBULAR_EXTINCTION', 'SFH']
          unused_tags = [unused_tags, 'STEPS_BOUNDS','DTIME_SF', 'PSI']

          ; Only remove tag if it is in the structure, otherwise an error will be thrown if
          ;   trying to remove not existent tag
          remove_tag = intarr(n_elements(unused_tags))
          for i=0, n_elements(unused_tags)-1 do remove_tag[i] = total(strmatch(strupcase(tags), unused_tags[i])) eq 1
          unused_tags = unused_tags[where(remove_tag, /null)]
          if n_elements(unused_tags) ne 0 then config_hash.remove, unused_tags
        end
     endcase



; ===================    Dust Attenuation Parameter Check   =============================================
 if n_elements(where(strupcase(tags) eq 'ATTEN_CURVE', /null)) eq 0 then $
   message, base_err+'ATTEN_CURVE tag is missing.'
 if size(config.ATTEN_CURVE, /type) ne 7 then $
   message, base_err+'ATTEN_CURVE tag must be of type string.'
 if size(config.ATTEN_CURVE, /n_dim) ne 0 then $
   message, base_err+'ATTEN_CURVE tag must be a scalar.'
 ; Doore21 attenuation requires stellar emission
 if strupcase(config.SSP) eq 'NONE' then $
   atten_curve_options = ['CALZETTI00', 'CALZETTI_MOD']
 if total(strupcase(config.ATTEN_CURVE) eq strupcase(atten_curve_options)) eq 0 then $
   message, base_err+"ATTEN_CURVE tag must be set to one of the following: '"+$
            strjoin(atten_curve_options, "', '")+"'."


     case strupcase(config.ATTEN_CURVE) of
     ; ===================    Calzetti00 Attenuation Parameter Check   ==========================================
       'CALZETTI00' : begin
           if n_elements(where(strupcase(tags) eq 'TAUV', /null)) eq 0 then $
             message, base_err+'TAUV tag is missing.'
           if size(config.TAUV, /type) ne 8 then $
             message, base_err+'TAUV tag must be of type structure.'
           tauv_range = [0.0d, !values.d_infinity]
           prior_check, config.TAUV, 'TAUV', 1, tauv_range, prior_options, prior_options_narg

           unused_tags = ['TAUV_DIFF', 'DELTA', 'TAUV_BC', 'UV_BUMP', $
                          'TAUB_F', 'F_CLUMP', 'COSI', 'B_TO_D', 'ROLD0_AGES']
           remove_tag = intarr(n_elements(unused_tags))
           for i=0, n_elements(unused_tags)-1 do remove_tag[i] = total(strmatch(strupcase(tags), unused_tags[i])) eq 1
           unused_tags = unused_tags[where(remove_tag, /null)]
           if n_elements(unused_tags) ne 0 then config_hash.remove, unused_tags
         end


     ; ===================    Modified Calzetti Attenuation Parameter Check   ===================================
       'CALZETTI_MOD': begin
           if n_elements(where(strupcase(tags) eq 'TAUV_DIFF', /null)) eq 0 then $
             message, base_err+'TAUV_DIFF tag is missing.'
           if size(config.TAUV_DIFF, /type) ne 8 then $
             message, base_err+'TAUV_DIFF tag must be of type structure.'
           tauv_diff_range = [0.0d, !values.d_infinity]
           prior_check, config.TAUV_DIFF, 'TAUV_DIFF', 1, tauv_diff_range, prior_options, $
                        prior_options_narg

           if n_elements(where(strupcase(tags) eq 'DELTA', /null)) eq 0 then $
             message, base_err+'DELTA tag is missing.'
           if size(config.DELTA, /type) ne 8 then $
             message, base_err+'DELTA tag must be of type structure.'
           delta_range = [-1.0d*!values.d_infinity, !values.d_infinity]
           prior_check, config.DELTA, 'DELTA', 1, delta_range, prior_options, prior_options_narg

           if n_elements(where(strupcase(tags) eq 'TAUV_BC', /null)) eq 0 then $
             message, base_err+'TAUV_BC tag is missing.'
           if size(config.TAUV_BC, /type) ne 8 then $
             message, base_err+'TAUV_BC tag must be of type structure.'
           tauv_bc_range = [0.0d, !values.d_infinity]
           prior_check, config.TAUV_BC, 'TAUV_BC', 1, tauv_bc_range, prior_options, prior_options_narg

           if n_elements(where(strupcase(tags) eq 'UV_BUMP', /null)) eq 0 then $
             message, base_err+'UV_BUMP tag is missing.'
           if size(config.UV_BUMP, /type) lt 2 or size(config.UV_BUMP, /type) gt 3 then $
             message, base_err+'UV_BUMP tag must be of type int.'
           if size(config.UV_BUMP, /n_dim) ne 0 then $
             message, base_err+'UV_BUMP tag must be a scalar.'
           if ~(config.UV_BUMP eq 0 or config.UV_BUMP eq 1) then $
             message, base_err+'UV_BUMP tag must be a flag with a value of 0 or 1.'

           unused_tags = ['TAUV', $
                          'TAUB_F', 'F_CLUMP', 'COSI', 'B_TO_D', 'ROLD0_AGES']
           remove_tag = intarr(n_elements(unused_tags))
           for i=0, n_elements(unused_tags)-1 do remove_tag[i] = total(strmatch(strupcase(tags), unused_tags[i])) eq 1
           unused_tags = unused_tags[where(remove_tag, /null)]
           if n_elements(unused_tags) ne 0 then config_hash.remove, unused_tags
         end


     ; ===================    Doore21 Attenuation Parameter Check   ======================================
       'DOORE21': begin
           if n_elements(where(strupcase(tags) eq 'TAUB_F', /null)) eq 0 then $
             message, base_err+'TAUB_F tag is missing.'
           if size(config.TAUB_F, /type) ne 8 then $
             message, base_err+'TAUB_F tag must be of type structure.'
           taub_f_range = [0.0d, 8.d]
           prior_check, config.TAUB_F, 'TAUB_F', 1, taub_f_range, prior_options, prior_options_narg

           if n_elements(where(strupcase(tags) eq 'F_CLUMP', /null)) eq 0 then $
             message, base_err+'F_CLUMP tag is missing.'
           if size(config.F_CLUMP, /type) ne 8 then $
             message, base_err+'F_CLUMP tag must be of type structure.'
           f_clump_range = [0.0d, 0.61d]
           prior_check, config.F_CLUMP, 'F_CLUMP', 1, f_clump_range, prior_options, prior_options_narg

           if n_elements(where(strupcase(tags) eq 'COSI', /null)) eq 0 then $
             message, base_err+'COSI tag is missing.'
           if size(config.COSI, /type) ne 8 then $
             message, base_err+'COSI tag must be of type structure.'
           cosi_range = [0.0d, 1.0d]
           prior_check, config.COSI, 'COSI', 1, cosi_range, prior_options, prior_options_narg

           if n_elements(where(strupcase(tags) eq 'B_TO_D', /null)) eq 0 then $
             message, base_err+'B_TO_D tag is missing.'
           if size(config.B_TO_D, /type) ne 8 then $
             message, base_err+'B_TO_D tag must be of type structure.'
           b_to_d_range = [0.0d, !values.d_infinity]
           prior_check, config.B_TO_D, 'B_TO_D', 1, b_to_d_range, prior_options, prior_options_narg

           if n_elements(where(strupcase(tags) eq 'ROLD0_AGES', /null)) eq 0 then $
             message, base_err+'ROLD0_AGES tag is missing.'
           if size(config.ROLD0_AGES, /type) lt 2 or size(config.ROLD0_AGES, /type) gt 5 then $
             message, base_err+'ROLD0_AGES tag must be of type int, float, or double.'
           if size(config.ROLD0_AGES, /n_dim) ne 1 then $
             message, base_err+'ROLD0_AGES tag must be a 1-D array.'
           if n_elements(config.ROLD0_AGES) ne Nsteps then $
             message, base_err+'ROLD0_AGES tag must have Nsteps ('+string(Nsteps, f='(I0)')+$
                      ') number of elements.'
           if total(config.ROLD0_AGES eq 0 or config.ROLD0_AGES eq 1) ne Nsteps then $
             message, base_err+'ROLD0_AGES tag must only contain values of 0 or 1.'

           unused_tags = ['TAUV', $
                          'TAUV_DIFF', 'DELTA', 'TAUV_BC', 'UV_BUMP']
           remove_tag = intarr(n_elements(unused_tags))
           for i=0, n_elements(unused_tags)-1 do remove_tag[i] = total(strmatch(strupcase(tags), unused_tags[i])) eq 1
           unused_tags = unused_tags[where(remove_tag, /null)]
           if n_elements(unused_tags) ne 0 then config_hash.remove, unused_tags
         end
     endcase




; ===================    Dust Emission Parameter Check   =============================================
 if n_elements(where(strupcase(tags) eq 'DUST_MODEL', /null)) eq 0 then $
   message, base_err+'DUST_MODEL tag is missing.'
 if size(config.DUST_MODEL, /type) ne 7 then $
   message, base_err+'DUST_MODEL tag must be of type string.'
 if size(config.DUST_MODEL, /n_dim) ne 0 then $
   message, base_err+'DUST_MODEL tag must be a scalar.'
 if total(strupcase(config.DUST_MODEL) eq strupcase(dust_emission_options)) eq 0 then $
   message, base_err+"DUST_MODEL tag must be set to one of the following: '"+$
            strjoin(dust_emission_options, "', '")+"'."


     case strupcase(config.DUST_MODEL) of
     ; ===================    DL07 Dust Emission Parameter Check   ======================================
       'DL07' : begin
           if n_elements(where(strupcase(tags) eq 'UMIN', /null)) eq 0 then $
             message, base_err+'UMIN tag is missing.'
           if size(config.UMIN, /type) ne 8 then $
             message, base_err+'UMIN tag must be of type structure.'
           umin_range = [0.10d, 25.d]
           prior_check, config.UMIN, 'UMIN', 1, umin_range, prior_options, prior_options_narg

           if n_elements(where(strupcase(tags) eq 'UMAX', /null)) eq 0 then $
             message, base_err+'UMAX tag is missing.'
           if size(config.UMAX, /type) ne 8 then $
             message, base_err+'UMAX tag must be of type structure.'
           umax_range = [1.d3, 3.d5]
           prior_check, config.UMAX, 'UMAX', 1, umax_range, prior_options, prior_options_narg

           if n_elements(where(strupcase(tags) eq 'ALPHA', /null)) eq 0 then $
             message, base_err+'ALPHA tag is missing.'
           if size(config.ALPHA, /type) ne 8 then $
             message, base_err+'ALPHA tag must be of type structure.'
           alpha_range = [-10.d, 4.0d]
           prior_check, config.ALPHA, 'ALPHA', 1, alpha_range, prior_options, prior_options_narg

           if n_elements(where(strupcase(tags) eq 'GAMMA', /null)) eq 0 then $
             message, base_err+'GAMMA tag is missing.'
           if size(config.GAMMA, /type) ne 8 then $
             message, base_err+'GAMMA tag must be of type structure.'
           gamma_range = [0.0d, 1.0d]
           prior_check, config.GAMMA, 'GAMMA', 1, gamma_range, prior_options, prior_options_narg

           if n_elements(where(strupcase(tags) eq 'QPAH', /null)) eq 0 then $
             message, base_err+'QPAH tag is missing.'
           if size(config.QPAH, /type) ne 8 then $
             message, base_err+'QPAH tag must be of type structure.'
           qpah_range = [4.7d-3, 4.58d-2]
           prior_check, config.QPAH, 'QPAH', 1, qpah_range, prior_options, prior_options_narg

           if ~(config.ENERGY_BALANCE) then begin
             if n_elements(where(strupcase(tags) eq 'LTIR', /null)) eq 0 then $
               message, base_err+'LTIR tag is missing.'
             if size(config.LTIR, /type) ne 8 then $
               message, base_err+'LTIR tag must be of type structure.'
             ltir_range = [0.0d, !values.d_infinity]
             prior_check, config.LTIR, 'LTIR', 1, ltir_range, prior_options, prior_options_narg
           endif else begin
             unused_tags = ['LTIR']
             remove_tag = total(strmatch(strupcase(tags), unused_tags[0])) eq 1
             unused_tags = unused_tags[where(remove_tag, /null)]
             if n_elements(unused_tags) ne 0 then config_hash.remove, unused_tags
           endelse
         end


       'NONE' : begin
           unused_tags = ['UMIN', 'UMAX', 'ALPHA', 'GAMMA', 'QPAH', 'LTIR']
           remove_tag = intarr(n_elements(unused_tags))
           for i=0, n_elements(unused_tags)-1 do remove_tag[i] = total(strmatch(strupcase(tags), unused_tags[i])) eq 1
           unused_tags = unused_tags[where(remove_tag, /null)]
           if n_elements(unused_tags) ne 0 then config_hash.remove, unused_tags
        end
     endcase




; ===================    Xray Emission Parameter Check   =============================================
 if n_elements(where(strupcase(tags) eq 'XRAY_EMISSION', /null)) eq 0 then $
   message, base_err+'XRAY_EMISSION tag is missing.'
 if size(config.XRAY_EMISSION, /type) lt 2 or size(config.XRAY_EMISSION, /type) gt 3 then $
   message, base_err+'XRAY_EMISSION tag must be of type int.'
 if size(config.XRAY_EMISSION, /n_dim) ne 0 then $
   message, base_err+'XRAY_EMISSION tag must be a scalar.'
 if ~(config.XRAY_EMISSION eq 0 or config.XRAY_EMISSION eq 1) then $
   message, base_err+'XRAY_EMISSION tag must be a flag with a value of 0 or 1.'
 if config.SSP eq 'NONE' and config.XRAY_EMISSION eq 1 then $
   message, base_err+'XRAY_EMISSION tag must be set to 0 if there is no stellar emission.'

 if config.XRAY_EMISSION then begin
     if n_elements(where(strupcase(tags) eq 'XRAY_UNIT', /null)) eq 0 then $
       message, base_err+'XRAY_UNIT tag is missing.'
     if size(config.XRAY_UNIT, /type) ne 7 then $
       message, base_err+'XRAY_UNIT tag must be of type string.'
     if size(config.XRAY_UNIT, /n_dim) ne 0 then $
       message, base_err+'XRAY_UNIT tag must be a scalar.'
     if total(strupcase(config.XRAY_UNIT) eq strupcase(xray_unit_options)) eq 0 then $
       message, base_err+"XRAY_UNIT tag must be set to one of the following: '"+$
                strjoin(xray_unit_options, "', '")+"'."

     if strupcase(config.XRAY_UNIT) eq 'COUNTS' then begin
       if n_elements(where(strupcase(tags) eq 'XRAY_UNC', /null)) eq 0 then $
         message, base_err+'XRAY_UNC tag is missing.'
       if size(config.XRAY_UNC, /type) ne 7 then $
         message, base_err+'XRAY_UNC tag must be of type string.'
       if size(config.XRAY_UNC, /n_dim) ne 0 then $
         message, base_err+'XRAY_UNC tag must be a scalar.'
       if total(strupcase(config.XRAY_UNC) eq strupcase(xray_unc_options)) eq 0 then $
         message, base_err+"XRAY_UNC tag must be set to one of the following: '"+$
                  strjoin(xray_unc_options, "', '")+"'."
     endif

     if n_elements(where(strupcase(tags) eq 'XRAY_ABS_MODEL', /null)) eq 0 then $
       message, base_err+'XRAY_ABS_MODEL tag is missing.'
     if size(config.XRAY_ABS_MODEL, /type) ne 7 then $
       message, base_err+'XRAY_ABS_MODEL tag must be of type string.'
     if size(config.XRAY_ABS_MODEL, /n_dim) ne 0 then $
       message, base_err+'XRAY_ABS_MODEL tag must be a scalar.'
     if total(strupcase(config.XRAY_ABS_MODEL) eq strupcase(xray_abs_options)) eq 0 then $
       message, base_err+"XRAY_ABS_MODEL tag must be set to one of the following: '"+$
                strjoin(xray_abs_options, "', '")+"'."

     if n_elements(where(strupcase(tags) eq 'NH', /null)) eq 0 then $
       message, base_err+'NH tag is missing.'
     if size(config.NH, /type) ne 8 then $
       message, base_err+'NH tag must be of type structure.'
     nh_range = [1.d-4, 1.d5]
     prior_check, config.NH, 'NH', 1, nh_range, prior_options, prior_options_narg


     if n_elements(where(strupcase(tags) eq 'XRAY_AGN_MODEL', /null)) eq 0 then $
       message, base_err+'XRAY_AGN_MODEL tag is missing.'
     if size(config.XRAY_AGN_MODEL, /type) ne 7 then $
       message, base_err+'XRAY_AGN_MODEL tag must be of type string.'
     if size(config.XRAY_AGN_MODEL, /n_dim) ne 0 then $
       message, base_err+'XRAY_AGN_MODEL tag must be a scalar.'
     if total(strupcase(config.XRAY_AGN_MODEL) eq strupcase(xray_agn_options)) eq 0 then $
       message, base_err+"XRAY_AGN_MODEL tag must be set to one of the following: '"+$
                strjoin(xray_agn_options, "', '")+"'."


     case strupcase(config.XRAY_AGN_MODEL) of
     ; ===================    SKIRTOR AGN Parameter Check   ======================================
       'QSOSED' : begin
           if n_elements(where(strupcase(tags) eq 'AGN_MASS', /null)) eq 0 then $
             message, base_err+'AGN_MASS tag is missing.'
           if size(config.AGN_MASS, /type) ne 8 then $
             message, base_err+'AGN_MASS tag must be of type structure.'
           agn_mass_range = [1.d5, 1.d10]
           prior_check, config.AGN_MASS, 'AGN_MASS', 1, agn_mass_range, prior_options, $
                        prior_options_narg

           if n_elements(where(strupcase(tags) eq 'AGN_LOGMDOT', /null)) eq 0 then $
             message, base_err+'AGN_LOGMDOT tag is missing.'
           if size(config.AGN_LOGMDOT, /type) ne 8 then $
             message, base_err+'AGN_LOGMDOT tag must be of type structure.'
           agn_logmdot_range = [-1.5d, 0.3d]
           prior_check, config.AGN_LOGMDOT, 'AGN_LOGMDOT', 1, agn_logmdot_range, prior_options, $
                        prior_options_narg
         end


       'PLAW' : begin
           unused_tags = ['AGN_MASS', 'AGN_LOGMDOT']
           remove_tag = intarr(n_elements(unused_tags))
           for i=0, n_elements(unused_tags)-1 do remove_tag[i] = total(strmatch(strupcase(tags), unused_tags[i])) eq 1
           unused_tags = unused_tags[where(remove_tag, /null)]
           if n_elements(unused_tags) ne 0 then config_hash.remove, unused_tags
         end

       'NONE' : begin
           unused_tags = ['AGN_MASS', 'AGN_LOGMDOT']
           remove_tag = intarr(n_elements(unused_tags))
           for i=0, n_elements(unused_tags)-1 do remove_tag[i] = total(strmatch(strupcase(tags), unused_tags[i])) eq 1
           unused_tags = unused_tags[where(remove_tag, /null)]
           if n_elements(unused_tags) ne 0 then config_hash.remove, unused_tags
         end
     endcase
 endif else begin
     unused_tags = ['XRAY_UNIT', 'XRAY_UNC', 'XRAY_ABS_MODEL', 'NH', 'XRAY_AGN_MODEL', $
                    'AGN_MASS', 'AGN_LOGMDOT']
     remove_tag = intarr(n_elements(unused_tags))
     for i=0, n_elements(unused_tags)-1 do remove_tag[i] = total(strmatch(strupcase(tags), unused_tags[i])) eq 1
     unused_tags = unused_tags[where(remove_tag, /null)]
     if n_elements(unused_tags) ne 0 then config_hash.remove, unused_tags
 endelse



; ===================    AGN Emission Parameter Check   =============================================
 if n_elements(where(strupcase(tags) eq 'AGN_MODEL', /null)) eq 0 then $
   message, base_err+'AGN_MODEL tag is missing.'
 if size(config.AGN_MODEL, /type) ne 7 then $
   message, base_err+'AGN_MODEL tag must be of type string.'
 if size(config.AGN_MODEL, /n_dim) ne 0 then $
   message, base_err+'AGN_MODEL tag must be a scalar.'
 ; Check if DOORE21 attenuation is used, since it is not compatible
 if strupcase(config.ATTEN_CURVE) eq 'DOORE21' then begin
   agn_emission_options = ['none']
   agn_message = "AGN_MODEL tag must be set to 'NONE', since it is incompatible with the 'DOORE21' "+$
                 'attenuation curve.'
 endif else begin
   agn_message = "AGN_MODEL tag must be set to one of the following: '"+$
                  strjoin(agn_emission_options, "', '")+"'."
 endelse
 if total(strupcase(config.AGN_MODEL) eq strupcase(agn_emission_options)) eq 0 then $
   message, base_err+agn_message


     case strupcase(config.AGN_MODEL) of
     ; ===================    SKIRTOR AGN Parameter Check   ======================================
       'SKIRTOR' : begin
           ; Complicated logic due to 'XRAY_AGN_MODEL' tag not existing if 'XRAY_EMISSION' tag is 0
           if config.XRAY_EMISSION then begin
             if strupcase(config.XRAY_AGN_MODEL) ne 'QSOSED' then begin
               if n_elements(where(strupcase(tags) eq 'LOG_L_AGN', /null)) eq 0 then $
                 message, base_err+'LOG_L_AGN tag is missing.'
               if size(config.LOG_L_AGN, /type) ne 8 then $
                 message, base_err+'LOG_L_AGN tag must be of type structure.'
               log_l_agn_range = [0.0d, 20.d0]
               prior_check, config.LOG_L_AGN, 'LOG_L_AGN', 1, log_l_agn_range, prior_options, $
                            prior_options_narg
             endif else begin
               unused_tags = ['LOG_L_AGN']
               remove_tag = total(strmatch(strupcase(tags), unused_tags[0])) eq 1
               unused_tags = unused_tags[where(remove_tag, /null)]
               if n_elements(unused_tags) ne 0 then config_hash.remove, unused_tags
             endelse
           endif else begin
             if n_elements(where(strupcase(tags) eq 'LOG_L_AGN', /null)) eq 0 then $
               message, base_err+'LOG_L_AGN tag is missing.'
             if size(config.LOG_L_AGN, /type) ne 8 then $
               message, base_err+'LOG_L_AGN tag must be of type structure.'
             log_l_agn_range = [0.0d, 20.d0]
             prior_check, config.LOG_L_AGN, 'LOG_L_AGN', 1, log_l_agn_range, prior_options, $
                          prior_options_narg
           endelse

           if n_elements(where(strupcase(tags) eq 'TAU97', /null)) eq 0 then $
             message, base_err+'TAU97 tag is missing.'
           if size(config.TAU97, /type) ne 8 then $
             message, base_err+'TAU97 tag must be of type structure.'
           tau97_range = [3.0d, 11.d0]
           prior_check, config.TAU97, 'TAU97', 1, tau97_range, prior_options, prior_options_narg

           if n_elements(where(strupcase(tags) eq 'AGN_COSI', /null)) eq 0 then $
             message, base_err+'AGN_COSI tag is missing.'
           if size(config.AGN_COSI, /type) ne 8 then $
             message, base_err+'AGN_COSI tag must be of type structure.'
           agn_cosi_range = [0.0d, 1.0d]
           prior_check, config.AGN_COSI, 'AGN_COSI', 1, agn_cosi_range, prior_options, prior_options_narg
         end


       'NONE' : begin
           unused_tags = ['LOG_L_AGN', 'TAU97', 'AGN_COSI']
           remove_tag = intarr(n_elements(unused_tags))
           for i=0, n_elements(unused_tags)-1 do remove_tag[i] = total(strmatch(strupcase(tags), unused_tags[i])) eq 1
           unused_tags = unused_tags[where(remove_tag, /null)]
           if n_elements(unused_tags) ne 0 then config_hash.remove, unused_tags
         end
     endcase



; ====== Additional check to make sure AGN and/or stellar models are set if using a dust model and energy balance
 if strupcase(config.AGN_MODEL) eq 'NONE' and strupcase(config.SSP) eq 'NONE' and $
    strupcase(config.DUST_MODEL) ne 'NONE' then $
   if config.ENERGY_BALANCE then $
     message, base_err+'ENERGY_BALANCE tag cannot be set if using a dust emission model '+$
              'with no stellar or AGN emission models. Please unset the ENERGY_BALANCE tag '+$
              'when only fitting a dust emission model.'



; ===================    Fitting Algorithm Parameter Check   =============================================
 if n_elements(where(strupcase(tags) eq 'METHOD', /null)) eq 0 then $
   message, base_err+'METHOD tag is missing.'
 if size(config.METHOD, /type) ne 7 then $
   message, base_err+'METHOD tag must be of type string.'
 if size(config.METHOD, /n_dim) ne 0 then $
   message, base_err+'METHOD tag must be a scalar.'
 if total(strupcase(config.METHOD) eq strupcase(fit_algorithm_options)) eq 0 then $
   message, base_err+"METHOD tag must be set to one of the following: '"+$
            strjoin(fit_algorithm_options, "', '")+"'."


     case 1 of
     ; ===================    Grid Fitting Algorithm Parameter Check   ==========================================
       ;strupcase(config.METHOD) eq 'GRID' : begin
       ;    if n_elements(where(strupcase(tags) eq 'GRID_ELEMENTS', /null)) eq 0 then $
       ;      message, base_err+'GRID_ELEMENTS tag is missing.'
       ;    if size(config.GRID_ELEMENTS, /type) lt 2 or size(config.GRID_ELEMENTS, /type) gt 3 then $
       ;      message, base_err+'GRID_ELEMENTS tag must be of type int.'
       ;    if size(config.GRID_ELEMENTS, /n_dim) ne 0 then $
       ;      message, base_err+'GRID_ELEMENTS tag must be a scalar.'
       ;    if config.GRID_ELEMENTS le 0 then $
       ;      message, base_err+'GRID_ELEMENTS tag must be a positive value.'
       ; end


     ; ===================    Inversion Fitting Algorithm Parameter Check   ===================================
       ;strupcase(config.METHOD) eq 'INVERSION' : begin
       ;    if n_elements(where(strupcase(tags) eq 'ALLOW_NEG_SFH', /null)) eq 0 then $
       ;      message, base_err+'ALLOW_NEG_SFH tag is missing.'
       ;    if size(config.ALLOW_NEG_SFH, /type) lt 2 or size(config.ALLOW_NEG_SFH, /type) gt 3 then $
       ;      message, base_err+'ALLOW_NEG_SFH tag must be of type int.'
       ;    if size(config.ALLOW_NEG_SFH, /n_dim) ne 0 then $
       ;      message, base_err+'ALLOW_NEG_SFH tag must be a scalar.'
       ;    if ~(config.ALLOW_NEG_SFH eq 0 or config.ALLOW_NEG_SFH eq 1) then $
       ;      message, base_err+'ALLOW_NEG_SFH tag must be a flag with a value of 0 or 1.'
       ;
       ;    if n_elements(where(strupcase(tags) eq 'INVERSION_GRID_ELEMENTS', /null)) eq 0 then $
       ;      message, base_err+'INVERSION_GRID_ELEMENTS tag is missing.'
       ;    if size(config.INVERSION_GRID_ELEMENTS, /type) lt 2 or size(config.INVERSION_GRID_ELEMENTS, /type) gt 3 then $
       ;      message, base_err+'INVERSION_GRID_ELEMENTS tag must be of type int.'
       ;    if size(config.INVERSION_GRID_ELEMENTS, /n_dim) ne 0 then $
       ;      message, base_err+'INVERSION_GRID_ELEMENTS tag must be a scalar.'
       ;    if config.INVERSION_GRID_ELEMENTS le 0 then $
       ;      message, base_err+'INVERSION_GRID_ELEMENTS tag must be a positive value.'
       ;  end


     ; ===================    MCMC Fitting Algorithm Parameter Check   ===================================
       strmatch(strupcase(config.METHOD), 'MCMC-*'): begin
           if n_elements(where(strupcase(tags) eq 'NTRIALS', /null)) eq 0 then $
             message, base_err+'NTRIALS tag is missing.'
           if size(config.NTRIALS, /type) lt 2 or size(config.NTRIALS, /type) gt 5 then $
             message, base_err+'NTRIALS tag must be of type int, float, or double.'
           if size(config.NTRIALS, /n_dim) ne 0 then $
             message, base_err+'NTRIALS tag must be a scalar.'
           if config.NTRIALS le 0 then $
             message, base_err+'NTRIALS tag must be a positive value.'

           if n_elements(where(strupcase(tags) eq 'NPARALLEL', /null)) eq 0 then $
             message, base_err+'NPARALLEL tag is missing.'
           if size(config.NPARALLEL, /type) lt 2 or size(config.NPARALLEL, /type) gt 5 then $
             message, base_err+'NPARALLEL tag must be of type int, float, or double.'
           if size(config.NPARALLEL, /n_dim) ne 0 then $
             message, base_err+'NPARALLEL tag must be a scalar.'
           if config.NPARALLEL le 0 then $
             message, base_err+'NPARALLEL tag must be a positive value.'

           if n_elements(where(strupcase(tags) eq 'C_STEP', /null)) eq 0 then $
             message, base_err+'C_STEP tag is missing.'
           if size(config.C_STEP, /type) lt 2 or size(config.C_STEP, /type) gt 5 then $
             message, base_err+'C_STEP tag must be of type int, float, or double.'
           if size(config.C_STEP, /n_dim) ne 0 then $
             message, base_err+'C_STEP tag must be a scalar.'
           if config.C_STEP le 0 then $
             message, base_err+'C_STEP tag must be a positive value.'

           if n_elements(where(strupcase(tags) eq 'TOLERANCE', /null)) eq 0 then $
             message, base_err+'TOLERANCE tag is missing.'
           if size(config.TOLERANCE, /type) lt 2 or size(config.TOLERANCE, /type) gt 5 then $
             message, base_err+'TOLERANCE tag must be of type int, float, or double.'
           if size(config.TOLERANCE, /n_dim) ne 0 then $
             message, base_err+'TOLERANCE tag must be a scalar.'
           if config.TOLERANCE le 0 then $
             message, base_err+'TOLERANCE tag must be a positive value.'

           case strupcase(config.METHOD) of
           ; ===================    Adaptive MCMC    ===================================
             'MCMC-ADAPTIVE' : begin
                 if n_elements(where(strupcase(tags) eq 'BETA_EXPONENT', /null)) eq 0 then $
                   message, base_err+'BETA_EXPONENT tag is missing.'
                 if size(config.BETA_EXPONENT, /type) lt 2 or size(config.BETA_EXPONENT, /type) gt 5 then $
                   message, base_err+'BETA_EXPONENT tag must be of type int, float, or double.'
                 if size(config.BETA_EXPONENT, /n_dim) ne 0 then $
                   message, base_err+'BETA_EXPONENT tag must be a scalar.'
                 if config.BETA_EXPONENT le 0 then $
                   message, base_err+'BETA_EXPONENT tag must be a positive value.'
               end

           ; ===================    Affine MCMC    ===================================
             'MCMC-AFFINE' : begin
                 if config.NPARALLEL lt 2 then $
                   message, base_err+'NPARALLEL tag must be a greater than 1 when using affine '+$
                                     'invariant MCMC. Ideally, it should be at least twice the '+$
                                     'number of free parameters.'

                 if n_elements(where(strupcase(tags) eq 'AFFINE_A', /null)) eq 0 then $
                   message, base_err+'AFFINE_A tag is missing.'
                 if size(config.AFFINE_A, /type) lt 2 or size(config.AFFINE_A, /type) gt 5 then $
                   message, base_err+'AFFINE_A tag must be of type int, float, or double.'
                 if size(config.AFFINE_A, /n_dim) ne 0 then $
                   message, base_err+'AFFINE_A tag must be a scalar.'
                 if config.AFFINE_A lt 1 then $
                   message, base_err+'AFFINE_A tag must be a value greater than or equal to 1.'
               end
           endcase

           unused_tags = ['NSOLVERS', 'FTOL', 'GTOL', 'XTOL', 'MAXITER']
           remove_tag = intarr(n_elements(unused_tags))
           for i=0, n_elements(unused_tags)-1 do remove_tag[i] = total(strmatch(strupcase(tags), unused_tags[i])) eq 1
           unused_tags = unused_tags[where(remove_tag, /null)]
           if n_elements(unused_tags) ne 0 then config_hash.remove, unused_tags

         end


     ; ===================    MPFIT Fitting Algorithm Parameter Check   ===================================
       strupcase(config.METHOD) eq 'MPFIT' : begin
           if n_elements(where(strupcase(tags) eq 'NSOLVERS', /null)) eq 0 then $
             message, base_err+'NSOLVERS tag is missing.'
           if size(config.NSOLVERS, /type) lt 2 or size(config.NSOLVERS, /type) gt 5 then $
             message, base_err+'NSOLVERS tag must be of type int, float, or double.'
           if size(config.NSOLVERS, /n_dim) ne 0 then $
             message, base_err+'NSOLVERS tag must be a scalar.'
           if config.NSOLVERS le 0 then $
             message, base_err+'NSOLVERS tag must be a positive value.'

           if n_elements(where(strupcase(tags) eq 'FTOL', /null)) eq 0 then $
             message, base_err+'FTOL tag is missing.'
           if size(config.FTOL, /type) lt 2 or size(config.FTOL, /type) gt 5 then $
             message, base_err+'FTOL tag must be of type int, float, or double.'
           if size(config.FTOL, /n_dim) ne 0 then $
             message, base_err+'FTOL tag must be a scalar.'
           if config.FTOL le 0 then $
             message, base_err+'FTOL tag must be a positive value.'

           if n_elements(where(strupcase(tags) eq 'GTOL', /null)) eq 0 then $
             message, base_err+'GTOL tag is missing.'
           if size(config.GTOL, /type) lt 2 or size(config.GTOL, /type) gt 5 then $
             message, base_err+'GTOL tag must be of type int, float, or double.'
           if size(config.GTOL, /n_dim) ne 0 then $
             message, base_err+'GTOL tag must be a scalar.'
           if config.GTOL le 0 then $
             message, base_err+'GTOL tag must be a positive value.'

           if n_elements(where(strupcase(tags) eq 'XTOL', /null)) eq 0 then $
             message, base_err+'XTOL tag is missing.'
           if size(config.XTOL, /type) lt 2 or size(config.XTOL, /type) gt 5 then $
             message, base_err+'XTOL tag must be of type int, float, or double.'
           if size(config.XTOL, /n_dim) ne 0 then $
             message, base_err+'XTOL tag must be a scalar.'
           if config.XTOL le 0 then $
             message, base_err+'XTOL tag must be a positive value.'

           if n_elements(where(strupcase(tags) eq 'MAXITER', /null)) eq 0 then $
             message, base_err+'MAXITER tag is missing.'
           if size(config.MAXITER, /type) lt 2 or size(config.MAXITER, /type) gt 5 then $
             message, base_err+'MAXITER tag must be of type int, float, or double.'
           if size(config.MAXITER, /n_dim) ne 0 then $
             message, base_err+'MAXITER tag must be a scalar.'
           if config.MAXITER le 0 then $
             message, base_err+'MAXITER tag must be a positive value.'

           unused_tags = ['NTRIALS', 'NPARALLEL', 'C_STEP', 'TOLERANCE', 'BETA_EXPONENT', 'AFFINE_A']
           remove_tag = intarr(n_elements(unused_tags))
           for i=0, n_elements(unused_tags)-1 do remove_tag[i] = total(strmatch(strupcase(tags), unused_tags[i])) eq 1
           unused_tags = unused_tags[where(remove_tag, /null)]
           if n_elements(unused_tags) ne 0 then config_hash.remove, unused_tags
         end
     endcase


;===========================   Post-processing Parameter Check   ===============================================
 if n_elements(where(strupcase(tags) eq 'KEEP_INTERMEDIATE_OUTPUT', /null)) eq 0 then $
   message, base_err+'KEEP_INTERMEDIATE_OUTPUT tag is missing.'
 if size(config.KEEP_INTERMEDIATE_OUTPUT, /type) lt 2 or size(config.KEEP_INTERMEDIATE_OUTPUT, /type) gt 3 then $
   message, base_err+'KEEP_INTERMEDIATE_OUTPUT tag must be of type int.'
 if size(config.KEEP_INTERMEDIATE_OUTPUT, /n_dim) ne 0 then $
   message, base_err+'KEEP_INTERMEDIATE_OUTPUT tag must be a scalar.'
 if ~(config.KEEP_INTERMEDIATE_OUTPUT eq 0 or config.KEEP_INTERMEDIATE_OUTPUT eq 1) then $
   message, base_err+'KEEP_INTERMEDIATE_OUTPUT tag must be a flag with a value of 0 or 1.'


     ; ===================    MCMC Post-processing Parameter Check    ===================================
     if (strupcase(config.METHOD) eq 'MCMC-ADAPTIVE' or strupcase(config.METHOD) eq 'MCMC-AFFINE') then begin
           if n_elements(where(strupcase(tags) eq 'BURN_IN', /null)) eq 0 then $
             message, base_err+'BURN_IN tag is missing.'
           if size(config.BURN_IN, /type) lt 2 or size(config.BURN_IN, /type) gt 5 then $
             message, base_err+'BURN_IN tag must be of type int, float, or double.'
           if size(config.BURN_IN, /n_dim) ne 0 then $
             message, base_err+'BURN_IN tag must be a scalar.'
           if config.BURN_IN lt 0 then $
             message, base_err+'BURN_IN tag must be a non-negative value.'
           if config.BURN_IN ge config.NTRIALS then $
             message, base_err+'BURN_IN tag must have a smaller value than NTRIALS tag.'

           if n_elements(where(strupcase(tags) eq 'THIN_FACTOR', /null)) eq 0 then $
             message, base_err+'THIN_FACTOR tag is missing.'
           if size(config.THIN_FACTOR, /type) lt 2 or size(config.THIN_FACTOR, /type) gt 5 then $
             message, base_err+'THIN_FACTOR tag must be of type int, float, or double.'
           if size(config.THIN_FACTOR, /n_dim) ne 0 then $
             message, base_err+'THIN_FACTOR tag must be a scalar.'
           if config.THIN_FACTOR lt 0 then $
             message, base_err+'THIN_FACTOR tag must be a non-negative value.'

           case strupcase(config.METHOD) of
             'MCMC-ADAPTIVE' : $
               if float(config.THIN_FACTOR) * config.FINAL_CHAIN_LENGTH gt float(config.NTRIALS) - config.BURN_IN then $
               message, base_err+'THIN_FACTOR tag is too large such that when the chain is truncated and thinned '+$
                        'there will not be enough trials left to have the specified value in FINAL_CHAIN_LENGTH tag. '+$
                        'Either increase NTRIALS, or decrease THIN_FACTOR or FINAL_CHAIN_LENGTH.'

             'MCMC-AFFINE' : $
               if (float(config.THIN_FACTOR) * config.FINAL_CHAIN_LENGTH)/config.NPARALLEL gt float(config.NTRIALS) - config.BURN_IN then $
               message, base_err+'THIN_FACTOR tag is too large such that when the chain is truncated, thinned, and merged '+$
                        'there will not be enough trials left to have the specified value in FINAL_CHAIN_LENGTH tag. '+$
                        'Either increase NTRIALS, or decrease THIN_FACTOR or FINAL_CHAIN_LENGTH.'
           endcase

           if n_elements(where(strupcase(tags) eq 'FINAL_CHAIN_LENGTH', /null)) eq 0 then $
             message, base_err+'FINAL_CHAIN_LENGTH tag is missing.'
           if size(config.FINAL_CHAIN_LENGTH, /type) lt 2 or size(config.FINAL_CHAIN_LENGTH, /type) gt 5 then $
             message, base_err+'FINAL_CHAIN_LENGTH tag must be of type int, float, or double.'
           if size(config.FINAL_CHAIN_LENGTH, /n_dim) ne 0 then $
             message, base_err+'FINAL_CHAIN_LENGTH tag must be a scalar.'
           if config.FINAL_CHAIN_LENGTH le 0 then $
             message, base_err+'FINAL_CHAIN_LENGTH tag must be a positive value.'
           if config.FINAL_CHAIN_LENGTH gt config.NTRIALS then $
             message, base_err+'FINAL_CHAIN_LENGTH tag must be less than or equal to the value in NTRIALS tag.'

           if n_elements(where(strupcase(tags) eq 'HIGH_RES_MODEL_FRACTION', /null)) eq 0 then $
             message, base_err+'HIGH_RES_MODEL_FRACTION tag is missing.'
           if size(config.HIGH_RES_MODEL_FRACTION, /type) lt 2 or size(config.HIGH_RES_MODEL_FRACTION, /type) gt 5 then $
             message, base_err+'HIGH_RES_MODEL_FRACTION tag must be of type int, float, or double.'
           if size(config.HIGH_RES_MODEL_FRACTION, /n_dim) ne 0 then $
             message, base_err+'HIGH_RES_MODEL_FRACTION tag must be a scalar.'
           if config.HIGH_RES_MODEL_FRACTION lt 0 or config.HIGH_RES_MODEL_FRACTION gt 1 then $
             message, base_err+'HIGH_RES_MODEL_FRACTION tag must be a value between 0 and 1.'

           if strupcase(config.METHOD) eq 'MCMC-AFFINE'then begin
             if n_elements(where(strupcase(tags) eq 'AFFINE_STRANDED_DEVIATION', /null)) eq 0 then $
               message, base_err+'AFFINE_STRANDED_DEVIATION tag is missing.'
             if size(config.AFFINE_STRANDED_DEVIATION, /type) lt 2 or size(config.AFFINE_STRANDED_DEVIATION, /type) gt 5 then $
               message, base_err+'AFFINE_STRANDED_DEVIATION tag must be of type int, float, or double.'
             if size(config.AFFINE_STRANDED_DEVIATION, /n_dim) ne 0 then $
               message, base_err+'AFFINE_STRANDED_DEVIATION tag must be a scalar.'
             if config.AFFINE_STRANDED_DEVIATION le 0 then $
               message, base_err+'AFFINE_STRANDED_DEVIATION tag must be a positive value.'
           endif
     endif else begin

           unused_tags = ['BURN_IN', 'THIN_FACTOR', 'FINAL_CHAIN_LENGTH', 'HIGH_RES_MODEL_FRACTION', 'AFFINE_STRANDED_DEVIATION']
           remove_tag = intarr(n_elements(unused_tags))
           for i=0, n_elements(unused_tags)-1 do remove_tag[i] = total(strmatch(strupcase(tags), unused_tags[i])) eq 1
           unused_tags = unused_tags[where(remove_tag, /null)]
           if n_elements(unused_tags) ne 0 then config_hash.remove, unused_tags

     endelse


 return, config_hash.tostruct(/recursive)

end
