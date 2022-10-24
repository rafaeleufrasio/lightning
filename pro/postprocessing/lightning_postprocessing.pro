pro lightning_postprocessing, input_dir, config, sed_id
;+
; Name
; ----
;   LIGHTNING_POSTPROCESSING
;
; Purpose
; -------
;   Post-processes the fitting results into more readily usable outputs.
;   The final post-processed output is saved to a table in the first
;   extension of a FITS file.
;
; Calling Sequence
; ----------------
;   ::
;
;       lightning_postprocessing, input_dir, config, sed_id
;
; Inputs
; ------
;   ``input_dir`` : string scalar
;       The path to the file containing the input SED data.
;   ``config`` : structure
;       A Lightning configuration structure. (See
;       ``lightning_configure_defaults.pro`` for details and contents.)
;   ``sed_id`` : string array(Nsed)
;       The ID of each SED.
;
; Output
; ------
;   A data table saved to the first extension of a FITS file with a file name (and path) of
;   ``<input_dir>/lightning_output/postprocessed_data_<date_in_UTC>.fits.gz``. The table contains the
;   post-processed data for all SEDs.
;   See :ref:`postprocessing-label` for full details and possible contents.
;
; Modification History
; --------------------
;   - 2022/07/12: Created (Keith Doore)
;   - 2022/07/29: Added check to remove occasional stranded walker in affine-MCMC from final distribution (Keith Doore)
;   - 2022/08/01: Fixed chain thinning bug (Keith Doore)
;   - 2022/08/01: Changed to element-wise flattening from segment-wise flattening for affine MCMC (Keith Doore)
;   - 2022/08/07: Fixed high resolution model indexing bug (Keith Doore)
;   - 2022/08/09: Fixed incorrect initialization size of high resolution xray models (Keith Doore)
;   - 2022/08/09: Fixed issue where ``config.HIGH_RES_MODEL_FRACTION=0`` was getting no models vs best fit (Keith Doore)
;   - 2022/08/11: Updated high resolution xray model function name (Keith Doore)
;   - 2022/08/17: Updated to include MPFIT outputs (Keith Doore)
;   - 2022/08/17: Included ``PARAMETER_NAMES`` tag to determine certain corresponding outputs (Keith Doore)
;   - 2022/08/18: Added progress printing (Keith Doore)
;   - 2022/08/18: Fixed burn-in and thinning issues with MCMC chain if ``Ntrials`` is small (Keith Doore)
;   - 2022/08/23: Added date to the end of output file name (Keith Doore)
;   - 2022/08/24: Updated to allow for different xray bandpasses, which were not checked at input (Keith Doore)
;   - 2022/08/24: Updated to make fixed parameters only be scalar vs array of same value (Keith Doore)
;   - 2022/09/01: Added handling for user-supplied X-ray count uncertainties (Erik B. Monson)
;   - 2022/09/01: Added ``chi2`` output in addition to ``lnprob`` (Erik B. Monson)
;   - 2022/09/01: Added X-ray emission to PPC (Erik B. Monson)
;   - 2022/09/14: Updates to allow fitting with X-ray fluxes (Erik B. Monson)
;   - 2022/09/21: Few bug fixes with X-ray fluxes (Keith Doore)
;   - 2022/09/22: Updated how we identify stranded walkers in affine MCMC (Keith Doore)
;   - 2022/09/23: Made ``autocorr_flag`` unique to affine MCMC (Keith Doore)
;   - 2022/09/26: Added ``DOF`` to MPFIT output (Keith Doore)
;   - 2022/10/24: Updated stranded walker search to use configuration input value (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(input_dir) eq 0 then message, 'Variable is undefined: INPUT_DIR.'
 if size(input_dir, /type) ne 7 then message, 'INPUT_DIR must be of type string.'
 if size(input_dir, /n_dim) ne 0 then message, 'INPUT_DIR must be a scalar.'
 if ~file_test(input_dir) then message, 'The input data directory does not exist. '+$
                                        'Please confirm specified directory is correct.'

 if n_elements(config) eq 0 then message, 'Variable is undefined: CONFIG.'
 if size(config, /type) ne 8 then message, 'CONFIG must be of type structure.'

 if n_elements(sed_id) eq 0 then message, 'Variable is not defined: SED_ID.'
 if size(sed_id, /type) ne 7 then message, 'SED_ID must be of type string.'
 if size(sed_id, /n_dim) gt 1 then message, 'SED_ID must be a scalar or 1D array.'



; Determine the number of filters and if using xray data determine max number of xray bandpasses
 Nsed = n_elements(sed_id)
 restore, input_dir+'lightning_output/input_sav_files/lightning_input_'+sed_id[0]+'.sav'
 Nfilters = n_elements(sed_data.FILTER_LABELS)
 if keyword_set(config.XRAY_EMISSION) then begin
   case strupcase(config.XRAY_UNIT) of
     'COUNTS': begin
         Nxray_max = !null
         for i=0, Nsed-1 do begin
           restore, input_dir+'lightning_output/input_sav_files/lightning_input_'+sed_id[i]+'.sav'
           Nxray_max = max([Nxray_max, n_elements(sed_data.NET_COUNTS)])
         endfor
     end
     'FLUX': Nxray_max = n_elements(sed_data.XRAY_LNU_OBS)
   endcase
 endif


; Extract the high resolution model wavelength
 void = lightning_model_lnu_highres(wave=wave, /get_wave)
 Nwave = n_elements(wave)
 void = lightning_model_lnuxray_highres(wave=xray_wave, /get_wave)
 Nxray_wave = n_elements(xray_wave)
 if strupcase(config.METHOD) eq 'MPFIT' then begin
   Nmodels = 1
   Nhighres_models = 1
 endif else begin
   Nmodels = config.FINAL_CHAIN_LENGTH
   Nhighres_models = long(config.FINAL_CHAIN_LENGTH * config.HIGH_RES_MODEL_FRACTION) > 1
 endelse


; Need to call generate_prior_struct once to get config_nopriors for generating the model data
;   Also gives the parameter names in priors without having to open a file
 priors = generate_prior_struct(config, sed_id[0], config_nopriors=config_nopriors)



; Create the blank output structure
 ; Prevents tags with arrays of length 1
 if Nmodels gt 1 then Nmodels_array = dblarr(Nmodels) else Nmodels_array = 0.d0

 ;Use hash to dynamically add tags
 out_hash = orderedhash()
 out_hash['SED_ID']          = ''
 out_hash['REDSHIFT']        = !values.D_NaN
 out_hash['LUMIN_DIST']      = !values.D_NaN
 out_hash['FILTER_LABELS']   = strarr(Nfilters)
 out_hash['WAVE_FILTERS']    = !values.D_NaN*dblarr(Nfilters)
 out_hash['LNU_OBS']         = !values.D_NaN*dblarr(Nfilters)
 out_hash['LNU_UNC']         = !values.D_NaN*dblarr(Nfilters)
 out_hash['LNU_MOD']         = !values.D_NaN*dblarr(Nfilters, Nmodels)
 out_hash['MODEL_UNC']       = !values.D_NaN
 out_hash['WAVE_HIRES']      = !values.D_NaN*dblarr(Nwave)
 out_hash['LNU_MOD_HIRES']   = !values.D_NaN*dblarr(Nwave, Nhighres_models)
 out_hash['LNPROB']          = !values.D_NaN * Nmodels_array
 out_hash['CHI2']            = !values.D_NaN * Nmodels_array
 out_hash['PARAMETER_NAMES'] = strarr(n_elements(priors.parameter_name))

 ; Add the free parameters
 ;  If a parameter has multiple priors (e.g., PSI), keep them as a single 2D array
 config_hash = orderedhash(config, /extract)
 foreach prior_substruct, config_hash, param_name do begin
   if size(prior_substruct, /type) eq 11 then begin
     if prior_substruct.haskey('PRIOR') then begin
       if n_elements(prior_substruct['PRIOR']) gt 1 then begin
         if total(strupcase(prior_substruct['PRIOR']) eq 'FIXED') eq n_elements(prior_substruct['PRIOR']) then begin
           out_hash[param_name] = !values.D_NaN*dblarr(n_elements(prior_substruct['PRIOR']))
         endif else begin
           out_hash[param_name] = !values.D_NaN*dblarr(n_elements(prior_substruct['PRIOR']), Nmodels)
         endelse
         if strupcase(config.METHOD) ne 'MPFIT' then begin
           out_hash[param_name+'_BESTFIT']     = !values.D_NaN*dblarr(n_elements(prior_substruct['PRIOR']))
           out_hash[param_name+'_PERCENTILES'] = !values.D_NaN*dblarr(n_elements(prior_substruct['PRIOR']), 3)
         endif else begin
           out_hash[param_name+'_UNC']     = !values.D_NaN*dblarr(n_elements(prior_substruct['PRIOR']))
         endelse
       endif else begin
         if strupcase(prior_substruct['PRIOR']) eq 'FIXED' then begin
           out_hash[param_name] = !values.D_NaN
         endif else begin
           out_hash[param_name] = !values.D_NaN * Nmodels_array
         endelse
         if strupcase(config.METHOD) ne 'MPFIT' then begin
           out_hash[param_name+'_BESTFIT']     = !values.D_NaN
           out_hash[param_name+'_PERCENTILES'] = !values.D_NaN*dblarr(3)
         endif else begin
           out_hash[param_name+'_UNC']     = !values.D_NaN
         endelse
       endelse
     endif
   endif
 endforeach


 case strupcase(config.SPS) of
   'PEGASE': begin
        case strupcase(config.SFH) of
          'NON-PARAMETRIC': begin
               Nsteps = n_elements(config.STEPS_BOUNDS) - 1
               out_hash['LNU_STARMOD']             = !values.D_NaN*dblarr(Nfilters, Nmodels)
               out_hash['LNU_STARMOD_UNRED']       = !values.D_NaN*dblarr(Nfilters, Nmodels)
               out_hash['LNU_STARMOD_HIRES']       = !values.D_NaN*dblarr(Nwave, Nhighres_models)
               out_hash['LNU_STARMOD_UNRED_HIRES'] = !values.D_NaN*dblarr(Nwave, Nhighres_models)
               out_hash['STEPS_BOUNDS']            = !values.D_NaN*dblarr(Nsteps + 1)
               out_hash['STEPS_MSTAR_COEFF']       = !values.D_NaN*dblarr(Nsteps)
             end
          'PARAMETRIC':
        endcase
      end
   'NONE':
 endcase

 case strupcase(config.DUST_MODEL) of
   'DL07': begin
        out_hash['LNU_DUSTMOD']       = !values.D_NaN*dblarr(Nfilters, Nmodels)
        out_hash['LNU_DUSTMOD_HIRES'] = !values.D_NaN*dblarr(Nwave, Nhighres_models)

        ; If LTIR is not a free parameter, still include its value in output
        if ~(out_hash.HasKey('LTIR')) then out_hash['LTIR'] = !values.D_NaN * Nmodels_array
      end
   'NONE':
 endcase

 if keyword_set(config.XRAY_EMISSION) then begin
   out_hash['XRAY_BANDPASS'] = !values.D_NaN*dblarr(2, Nxray_max)
   out_hash['GALACTIC_NH']   = !values.D_NaN
   out_hash['LNU_XRAYMOD']   = !values.D_NaN*dblarr(Nxray_max, Nmodels)

   case strupcase(config.XRAY_UNIT) of
     'COUNTS': begin
         out_hash['XRAY_EXPOSURE']     = !values.D_NaN*dblarr(Nxray_max)
         out_hash['NET_COUNTS']        = !values.D_NaN*dblarr(Nxray_max)
         out_hash['NET_COUNTS_UNC']    = !values.D_NaN*dblarr(Nxray_max)
         out_hash['XRAY_COUNTS_MOD']    = !values.D_NaN*dblarr(Nxray_max, Nmodels)
         out_hash['LNU_XRAY_OBS']    = !values.D_NaN*dblarr(Nxray_max, Nmodels)
         out_hash['LNU_XRAY_UNC']    = !values.D_NaN*dblarr(Nxray_max, Nmodels)
       end
     'FLUX': begin
         out_hash['LNU_XRAY_OBS']    = !values.D_NaN*dblarr(Nxray_max)
         out_hash['LNU_XRAY_UNC']    = !values.D_NaN*dblarr(Nxray_max)
       end
   endcase

   out_hash['WAVE_XRAYMOD_HIRES']            = !values.D_NaN*dblarr(Nxray_wave)
   out_hash['LNU_XRAYMOD_HIRES']             = !values.D_NaN*dblarr(Nxray_wave, Nhighres_models)
   out_hash['LNU_XRAYMOD_STAR_HIRES']        = !values.D_NaN*dblarr(Nxray_wave, Nhighres_models)
   out_hash['LNU_XRAYMOD_STAR_UNABS_HIRES']  = !values.D_NaN*dblarr(Nxray_wave, Nhighres_models)

   if strupcase(config.XRAY_AGN_MODEL) ne 'NONE' then begin
     out_hash['LNU_XRAYMOD_AGN_HIRES']         = !values.D_NaN*dblarr(Nxray_wave, Nhighres_models)
     out_hash['LNU_XRAYMOD_AGN_UNABS_HIRES']   = !values.D_NaN*dblarr(Nxray_wave, Nhighres_models)
   endif
 endif

 case strupcase(config.AGN_MODEL) of
   'SKIRTOR': begin
        out_hash['LNU_AGNMOD']             = !values.D_NaN*dblarr(Nfilters, Nmodels)
        out_hash['LNU_AGNMOD_UNRED']       = !values.D_NaN*dblarr(Nfilters, Nmodels)
        out_hash['LNU_AGNMOD_HIRES']       = !values.D_NaN*dblarr(Nwave, Nhighres_models)
        out_hash['LNU_AGNMOD_UNRED_HIRES'] = !values.D_NaN*dblarr(Nwave, Nhighres_models)
      end
   'NONE':
 endcase

 if strupcase(config.METHOD) ne 'MPFIT' then begin
   out_hash['LNPROB_BESTFIT']   = !values.D_NaN
   out_hash['CHI2_BESTFIT']     = !values.D_NaN
   out_hash['ACCEPTANCE_FRAC']  = !values.D_NaN*dblarr(config.NPARALLEL)
   out_hash['AUTOCORR_TIME']    = !values.D_NaN*dblarr(n_elements(priors.parameter_name))
   out_hash['BURN_IN_AUTOCORR'] = 0L
   out_hash['ACCEPTANCE_FLAG']  = intarr(config.NPARALLEL)
   out_hash['CONVERGENCE_FLAG'] = 0
   out_hash['SHORT_CHAIN_FLAG'] = 0
   out_hash['PVALUE']           = !values.D_NaN

   if strupcase(config.METHOD) eq 'MCMC-AFFINE' then begin
     out_hash['AUTOCORR_FLAG']    = intarr(n_elements(priors.parameter_name))
     out_hash['STRANDED_FLAG']    = intarr(config.NPARALLEL)
   endif

   if strupcase(config.METHOD) eq 'MCMC-ADAPTIVE' then begin
     out_hash['GELMAN_RUBIN_R_HAT']  = !values.D_NaN*dblarr(n_elements(priors.parameter_name))
     out_hash['BROOKS_GELMAN_R_HAT'] = !values.D_NaN
     out_hash['GELMAN_RUBIN_FLAG']   = intarr(n_elements(priors.parameter_name))
     out_hash['BROOKS_GELMAN_FLAG']  = 0
     out_hash['NCHAIN_R_HAT']        = 0L
   endif
 endif else begin
   out_hash['COVARIANCE']       = !values.D_NaN*dblarr(n_elements(priors.parameter_name), n_elements(priors.parameter_name))
   out_hash['STATUS']           = intarr(config.NSOLVERS)
   out_hash['STATUS_FLAG']      = intarr(config.NSOLVERS)
   out_hash['ERROR_MSG']        = strarr(config.NSOLVERS)
   out_hash['ITER_FRAC']        = !values.D_NaN*dblarr(config.NSOLVERS)
   out_hash['ITER_FLAG']        = intarr(config.NSOLVERS)
   out_hash['STUCK_FRAC']       = !values.D_NaN
   out_hash['STUCK_FLAG']       = 0
   out_hash['SIMILAR_FRAC']     = !values.D_NaN*dblarr(n_elements(priors.parameter_name))
   out_hash['SIMILAR_FLAG']     = intarr(n_elements(priors.parameter_name))
   out_hash['NFUNC_EVALS']      = intarr(config.NSOLVERS)
   out_hash['CONVERGENCE_FLAG'] = 0
   out_hash['PVALUE']           = !values.D_NaN*dblarr(config.NSOLVERS)
   out_hash['DOF']              = !values.D_NaN
 endelse


; Convert hash to structure as to replicate it and get array of structures
 out = out_hash.ToStruct()
 out = replicate(out, Nsed)



; Read in fitting algorithm results and compute commonly used properties from fit parameters
 t0 = systime(/sec)
 for i=0, Nsed-1 do begin
   restore, input_dir+'lightning_output/input_sav_files/lightning_input_'+sed_id[i]+'.sav'
   restore, input_dir+'lightning_output/output_sav_files/lightning_output_'+sed_id[i]+'.sav'

   case strupcase(config.METHOD) of
   ; If the affine MCMC algorithm, truncate, thin, flatten, then make chain to final size.
     'MCMC-AFFINE': begin
         if config.BURN_IN eq 0 then burn_in = convergence_metric.burn_in_autocorr else burn_in = config.BURN_IN
         ; If automatic, the burn-in may be longer than Ntrials.
         ;   If this happens, create a whole flattened chain of NaN and set a flag telling that this happened.
         if burn_in ge config.NTRIALS then begin
           short_chain_flag = 1
           parameters = replicate(!values.D_NaN, n_elements(parameter_name), config.FINAL_CHAIN_LENGTH)
           lnprob = replicate(!values.D_NaN, config.FINAL_CHAIN_LENGTH)
           chi2 = replicate(!values.D_NaN, config.FINAL_CHAIN_LENGTH)
         endif else begin
           chain = chain[*, burn_in:*, *]
           lnprob_chain = lnprob_chain[burn_in:*, *]
           chi2_chain = chi2_chain[burn_in:*, *]

           ; Check for stranded walkers and remove them.
           med_accept = median(convergence_metric.acceptance_frac)
           std_accept = stddev(convergence_metric.acceptance_frac)
           stranded = where(convergence_metric.acceptance_frac lt med_accept - config.AFFINE_STRANDED_DEVIATION*std_accept, $
                            Nstranded, comp=free, ncomp=Nfree, /null)
           stranded_flag = intarr(config.NPARALLEL)
           if Nstranded gt 0 then stranded_flag[stranded] = 1
           chain = chain[*, *, free]
           lnprob_chain = lnprob_chain[*, free]
           chi2_chain = chi2_chain[*, free]
           Nchain = Nfree

           if config.THIN_FACTOR eq 0 then thin_factor = max(convergence_metric.autocorr_time)*0.5d else $
             thin_factor = config.THIN_FACTOR
           thinning = [0:(config.NTRIALS-burn_in)-1:ceil(thin_factor, /L64)]
           chain = chain[*, thinning, *]
           lnprob_chain = lnprob_chain[thinning, *]
           chi2_chain = chi2_chain[thinning, *]

           flatten_length = long(Nchain) * n_elements(thinning)
           ; Need to transpose chains to flatten element-wise vs chain segment-wise
           chain = reform(transpose(chain, [0,2,1]), n_elements(parameter_name), flatten_length)
           lnprob_chain = reform(transpose(lnprob_chain), flatten_length)
           chi2_chain = reform(transpose(chi2_chain), flatten_length)

           ; If automatic, the burn-in and thinning may result in a final chain shorter than specified.
           ;   If this happens, keep what remains of the chain and set the rest to NaN.
           ;   Additionally, set a flag telling that this happened.
           if config.FINAL_CHAIN_LENGTH gt flatten_length then begin
             short_chain_flag = 1
             parameters = replicate(!values.D_NaN, n_elements(parameter_name), config.FINAL_CHAIN_LENGTH)
             lnprob = replicate(!values.D_NaN, config.FINAL_CHAIN_LENGTH)
             chi2 = replicate(!values.D_NaN, config.FINAL_CHAIN_LENGTH)
             parameters[*, 0:flatten_length-1] = chain
             lnprob[0:flatten_length-1] = lnprob_chain
             chi2[0:flatten_length-1] = chi2_chain
           endif else begin
             short_chain_flag = 0
             parameters = chain[*, -1*config.FINAL_CHAIN_LENGTH:*]
             lnprob = lnprob_chain[-1*config.FINAL_CHAIN_LENGTH:*]
             chi2 = chi2_chain[-1*config.FINAL_CHAIN_LENGTH:*]
           endelse
         endelse
       end

   ; If the adaptive MCMC algorithm, truncate, thin, select best fit chain, then make chain to final size.
     'MCMC-ADAPTIVE': begin
         if config.BURN_IN eq 0 then burn_in = convergence_metric.burn_in_autocorr else burn_in = config.BURN_IN
         ; If automatic, the burn-in may be longer than Ntrials.
         ;   If this happens, create a whole thinned single chain of NaN and set a flag telling that this happened.
         if burn_in ge config.NTRIALS then begin
           short_chain_flag = 1
           parameters = replicate(!values.D_NaN, n_elements(parameter_name), config.FINAL_CHAIN_LENGTH)
           lnprob = replicate(!values.D_NaN, config.FINAL_CHAIN_LENGTH)
           chi2 = replicate(!values.D_NaN, config.FINAL_CHAIN_LENGTH)
         endif else begin
           chain = chain[*, burn_in:*, *]
           lnprob_chain = lnprob_chain[burn_in:*, *]
           chi2_chain = chi2_chain[burn_in:*, *]

           if config.THIN_FACTOR eq 0 then thin_factor = max(convergence_metric.autocorr_time)*0.5d else $
             thin_factor = config.THIN_FACTOR
           thinning = [0:(config.NTRIALS-burn_in)-1:ceil(thin_factor, /L64)]
           chain = chain[*, thinning, *]
           lnprob_chain = lnprob_chain[thinning, *]
           chi2_chain = chi2_chain[thinning, *]
           thinned_length = n_elements(thinning)

           ; Do not flatten adaptive chain, rather select the parallel chain with the best fit model
           bestfit_chain = where(max(lnprob_chain, dim=1) eq max(lnprob_chain))
           chain = chain[*, *, bestfit_chain]
           lnprob_chain = lnprob_chain[*, bestfit_chain]
           chi2_chain = chi2_chain[*, bestfit_chain]

           ; If automatic, the burn-in and thinning may result in a final chain shorter than specified.
           ;   If this happens, keep what remains of the chain and set the rest to NaN.
           ;   Additionally, set a flag telling that this happened.
           if config.FINAL_CHAIN_LENGTH gt thinned_length then begin
             short_chain_flag = 1
             parameters = replicate(!values.D_NaN, n_elements(parameter_name), config.FINAL_CHAIN_LENGTH)
             lnprob = replicate(!values.D_NaN, config.FINAL_CHAIN_LENGTH)
             chi2 = replicate(!values.D_NaN, config.FINAL_CHAIN_LENGTH)
             parameters[*, 0:thinned_length-1] = chain
             lnprob[0:thinned_length-1] = lnprob_chain
             chi2[0:thinned_length-1] = chi2_chain
           endif else begin
             short_chain_flag = 0
             parameters = chain[*, -1*config.FINAL_CHAIN_LENGTH:*]
             lnprob = lnprob_chain[-1*config.FINAL_CHAIN_LENGTH:*]
             chi2 = chi2_chain[-1*config.FINAL_CHAIN_LENGTH:*]
           endelse
         endelse
       end

   ; If the MPFIT algorithm, select best fit solver.
     'MPFIT': begin
         bestfit_solver   = where(lnprob_mpfit eq max(lnprob_mpfit))
         parameters       = parameters_mpfit[*, bestfit_solver]
         parameters_error = parameters_error[*, bestfit_solver]
         covariance       = covariance[*, *, bestfit_solver]
         lnprob           = lnprob_mpfit[bestfit_solver]
         chi2             = -2.d * lnprob_mpfit[bestfit_solver]
       end
   endcase


   finite_model_idc = where(finite(lnprob), Nfinite_models, /null)
   ; Sort by negative lnprob to have best fit at front of array
   sorted_lnprob_idc = sort(-1.d*lnprob, /l64)
   bestfit_idc = sorted_lnprob_idc[0]
   highres_idc = sorted_lnprob_idc[0:(Nhighres_models - 1)]
   finite_highres_model_idc = highres_idc[where(finite(lnprob[highres_idc]), Nfinite_highres_models, /null)]


   Lnu_mod = replicate(!values.D_NaN, Nfilters, Nmodels)
   if (Nfinite_models gt 0) then begin
     ; Get the qsosed L2500 if need be, so that the AGN model luminosity is correct.
     L2500 = !null
     if config.XRAY_EMISSION and (strupcase(config.AGN_MODEL) eq 'SKIRTOR') then begin
       if strupcase(config.XRAY_AGN_MODEL) eq 'QSOSED' then begin
         agn_mass    = reform(parameters[(where(parameter_name eq 'AGN_MASS', /null))[0], finite_model_idc])
         agn_logmdot = reform(parameters[(where(parameter_name eq 'AGN_LOGMDOT', /null))[0], finite_model_idc])

         L2500 = qsosed_L2500(models.xray_models, agn_mass, agn_logmdot)
         L2500_qsosed_flag = 1
       endif
     endif
     Lnu_mod[*, finite_model_idc] = lightning_model_lnu(parameters[*, finite_model_idc], parameter_name, models, $
                                       Lnu_stellar=Lnu_stellar, Lnu_unred_stellar=Lnu_unred_stellar, Lnu_AGN=Lnu_AGN, $
                                       Lnu_unred_AGN=Lnu_unred_AGN, Lnu_dust=Lnu_dust, LTIR=LTIR, $
                                       Lbol_AGN_model=Lbol_AGN_model, L2500=L2500, _extra=config_nopriors)
   endif

   Lnu_mod_highres = replicate(!values.D_NaN, Nwave, Nhighres_models)
   if (Nfinite_highres_models gt 0) then begin
     ; Indexing L2500 by finite_highres_model_idc works since L2500 already only contains finite values.
     ;   Therefore, finite_highres_model_idc will either take the whole L2500 or just the bestfit ordered subset.
     if ~keyword_set(L2500_qsosed_flag) then L2500_highres = !null else L2500_highres = L2500[finite_highres_model_idc]

     Lnu_mod_highres[*, 0:Nfinite_highres_models-1] = lightning_model_lnu_highres(parameters[*, finite_highres_model_idc], $
                                  parameter_name, models, Lnu_stellar_highres=Lnu_stellar_highres, $
                                  Lnu_unred_stellar_highres=Lnu_unred_stellar_highres, $
                                  Lnu_AGN_highres=Lnu_AGN_highres, Lnu_unred_AGN_highres=Lnu_unred_AGN_highres, $
                                  Lnu_dust_highres=Lnu_dust_highres, Lbol_AGN_model=Lbol_AGN_model_highres, $
                                  L2500=L2500_highres, _extra=config_nopriors)
   endif


   if config.XRAY_EMISSION then begin
     Nxray = (size(sed_data.XRAY_BANDPASS, /dim))[1]
     xray_counts_mod = replicate(!values.D_NaN, Nxray, Nmodels)
     if strupcase(config.XRAY_UNIT eq 'COUNTS') then begin
       if (Nfinite_models gt 0) then $
         xray_counts_mod[*, finite_model_idc] = lightning_model_counts(parameters[*, finite_model_idc], parameter_name, models, $
                                                                     Lbol_AGN_model=Lbol_AGN_model, L2500=L2500, $
                                                                     _extra=config_nopriors)
     endif

     Lnu_xray_mod = replicate(!values.D_NaN, Nxray, Nmodels)
     if (Nfinite_models gt 0) then $
        Lnu_xray_mod[*, finite_model_idc] = lightning_model_lnuxray(parameters[*, finite_model_idc], parameter_name, models, $
                                                                     Lbol_AGN_model=Lbol_AGN_model, L2500=L2500, $
                                                                     Lnu_xray_stellar=Lnu_xray_stellar, Lnu_xray_unabs_stellar=Lnu_xray_unabs_stellar, $
                                                                     Lnu_xray_AGN=Lnu_xray_AGN, Lnu_xray_unabs_AGN=Lnu_xray_unabs_AGN, $
                                                                     _extra=config_nopriors)


     Lnu_xray_highres = replicate(!values.D_NaN, Nxray_wave, Nhighres_models)
     if (Nfinite_highres_models gt 0) then $
       Lnu_xray_highres[*, 0:Nfinite_highres_models-1] = lightning_model_lnuxray_highres(parameters[*, finite_highres_model_idc], $
                                  parameter_name, models, Lbol_AGN_model=Lbol_AGN_model_highres, L2500=L2500_highres, $
                                  Lnu_xray_stellar_highres=Lnu_xray_stellar_highres, $
                                  Lnu_xray_unabs_stellar_highres=Lnu_xray_unabs_stellar_highres, $
                                  Lnu_xray_AGN_highres=Lnu_xray_AGN_highres, $
                                  Lnu_xray_unabs_AGN_highres=Lnu_xray_unabs_AGN_highres, $
                                  _extra=config_nopriors)

   endif


  ; Generate goodness of fit statistic
   if strupcase(config.METHOD) ne 'MPFIT' then begin
     if (config.XRAY_EMISSION) then begin
       case strupcase(config.XRAY_UNIT) of
         'COUNTS': pvalue = ppc(2e3, sed_data.LNU_OBS, sed_data.LNU_UNC, Lnu_mod, lnprob, $
                                counts_obs=sed_data.NET_COUNTS, counts_unc=sed_data.NET_COUNTS_UNC, $
                                counts_predict=xray_counts_mod, _extra=config)
         'FLUX': pvalue = ppc(2e3, [sed_data.LNU_OBS, sed_data.XRAY_LNU_OBS], $
                                   [sed_data.LNU_UNC, sed_data.XRAY_LNU_UNC], $
                                   [Lnu_mod, Lnu_xray_mod], lnprob, _extra=config)
       endcase
     endif else begin
         pvalue = ppc(2e3, sed_data.LNU_OBS, sed_data.LNU_UNC, Lnu_mod, lnprob, _extra=config)
     endelse
   endif



 ; Place calculated and read-in values into output structure
   out[i].SED_ID          = sed_id[i]
   out[i].REDSHIFT        = sed_data.REDSHIFT
   out[i].LUMIN_DIST      = sed_data.LUMIN_DIST
   out[i].FILTER_LABELS   = sed_data.FILTER_LABELS
   out[i].WAVE_FILTERS    = models.stellar_models.WAVE_FILTERS
   out[i].LNU_OBS         = sed_data.LNU_OBS
   out[i].LNU_UNC         = sed_data.LNU_UNC
   out[i].LNU_MOD         = Lnu_mod
   out[i].MODEL_UNC       = config.MODEL_UNC
   out[i].WAVE_HIRES      = wave
   out[i].LNU_MOD_HIRES   = Lnu_mod_highres
   out[i].LNPROB          = lnprob
   out[i].CHI2            = chi2
   out[i].PARAMETER_NAMES[0:(n_elements(parameter_name)-1)] = parameter_name


   ; Set free parameters values by indexing the out structure
   out_tags = tag_names(out)
   for tag=0, n_elements(out_tags)-1 do begin
     ; Need to be specific in matching as we don't want to match tag_BESTFIT just tag_[0-9]* or tag
     parameter_name_mask = stregex(parameter_name, '^'+out_tags[tag]+'(_[0-9]*)?$', /bool)
     parameter_name_idc = where(parameter_name_mask, num, /null)
     if num gt 0 then begin
       if strupcase(config.METHOD) ne 'MPFIT' then begin
         out[i].(where(out_tags eq out_tags[tag]+'_BESTFIT'))[0:(num-1)] = parameters[parameter_name_idc, bestfit_idc]
         out_percent = !values.D_NaN*dblarr(3, num)
         for j=0, num-1 do begin
           out_percent[*, j] = percentile(reform(parameters[parameter_name_idc[j], finite_model_idc]), [0.16d, 0.5d, 0.84d])
         endfor
         if num gt 1 then begin
           out[i].(where(out_tags eq out_tags[tag]+'_PERCENTILES'))[0:(num-1), *] = transpose(out_percent)
         endif else out[i].(where(out_tags eq out_tags[tag]+'_PERCENTILES')) = transpose(out_percent)
       endif else begin
          out[i].(where(out_tags eq out_tags[tag]+'_UNC'))[0:(num-1)] = parameters_error[parameter_name_idc]
       endelse

       if num gt 1 then begin
         if total(finite((priors.FIXED)[parameter_name_idc])) eq n_elements(parameter_name_idc) then begin
           out[i].(tag)[0:(num-1)] = reform(parameters[parameter_name_idc, 0])
         endif else begin
           out[i].(tag)[0:(num-1), *] = reform(parameters[parameter_name_idc, *])
         endelse
       endif else begin
         if finite((priors.FIXED)[parameter_name_idc]) then begin
           out[i].(tag) = reform(parameters[parameter_name_idc, 0])
         endif else begin
           out[i].(tag) = reform(parameters[parameter_name_idc, *])
         endelse
       endelse
     endif
   endfor

   case strupcase(config.SPS) of
     'PEGASE': begin
          case strupcase(config.SFH) of
            'NON-PARAMETRIC': begin
                 out[i].LNU_STARMOD[*, finite_model_idc]       = Lnu_stellar
                 out[i].LNU_STARMOD_UNRED[*, finite_model_idc] = Lnu_unred_stellar
                 out[i].LNU_STARMOD_HIRES[*, 0:Nfinite_highres_models-1]       = Lnu_stellar_highres
                 out[i].LNU_STARMOD_UNRED_HIRES[*, 0:Nfinite_highres_models-1] = Lnu_unred_stellar_highres

                 ; Need to recalculate as true steps may have been truncated from input due to redshift
                 Nsteps_new = n_elements(models.stellar_models.BOUNDS) - 1
                 out[i].STEPS_BOUNDS[0:(Nsteps_new)]        = models.stellar_models.BOUNDS
                 out[i].STEPS_MSTAR_COEFF[0:(Nsteps_new-1)] = models.stellar_models.MSTAR
               end
            'PARAMETRIC':
          endcase
        end
     'NONE':
   endcase

   case strupcase(config.DUST_MODEL) of
     'DL07': begin
          out[i].LNU_DUSTMOD[*, finite_model_idc] = Lnu_dust
          out[i].LNU_DUSTMOD_HIRES[*, 0:Nfinite_highres_models-1] = Lnu_dust_highres
          ; Do this even if LTIR is a free parameter, since it will be the same values.
          out[i].LTIR[finite_model_idc] = LTIR
        end
     'NONE':
   endcase

   if keyword_set(config.XRAY_EMISSION) then begin
     out[i].XRAY_BANDPASS[*, 0:Nxray-1]   = sed_data.XRAY_BANDPASS
     out[i].GALACTIC_NH     = sed_data.GALACTIC_NH
     out[i].LNU_XRAYMOD[0:Nxray-1, *]  = Lnu_xray_mod

     case strupcase(config.XRAY_UNIT) of
         'COUNTS': begin
             out[i].XRAY_EXPOSURE[0:Nxray-1]   = sed_data.XRAY_EXPOSURE
             out[i].NET_COUNTS[0:Nxray-1]      = sed_data.NET_COUNTS
             out[i].NET_COUNTS_UNC[0:Nxray-1]  = sed_data.NET_COUNTS_UNC
             out[i].XRAY_COUNTS_MOD[0:Nxray-1, *] = xray_counts_mod
             out[i].LNU_XRAY_OBS[0:Nxray-1, finite_model_idc] = rebin(sed_data.NET_COUNTS, Nxray, Nfinite_models) / $
                                                                xray_counts_mod[*, finite_model_idc] * Lnu_xray_mod[*, finite_model_idc]
             out[i].LNU_XRAY_UNC[0:Nxray-1, finite_model_idc] = rebin(sed_data.NET_COUNTS_UNC, Nxray, Nfinite_models) / $
                                                                xray_counts_mod[*, finite_model_idc] * Lnu_xray_mod[*, finite_model_idc]
         end
         'FLUX': begin
            out[i].LNU_XRAY_OBS[0:Nxray-1] = sed_data.XRAY_LNU_OBS
            out[i].LNU_XRAY_UNC[0:Nxray-1] = sed_data.XRAY_LNU_UNC
        end
     endcase
     out[i].WAVE_XRAYMOD_HIRES  = xray_wave
     out[i].LNU_XRAYMOD_HIRES   = Lnu_xray_highres
     out[i].LNU_XRAYMOD_STAR_HIRES[*, 0:Nfinite_highres_models-1]        = Lnu_xray_stellar_highres
     out[i].LNU_XRAYMOD_STAR_UNABS_HIRES[*, 0:Nfinite_highres_models-1]  = Lnu_xray_unabs_stellar_highres

     if strupcase(config.XRAY_AGN_MODEL) ne 'NONE' then begin
       out[i].LNU_XRAYMOD_AGN_HIRES[*, 0:Nfinite_highres_models-1]         = Lnu_xray_AGN_highres
       out[i].LNU_XRAYMOD_AGN_UNABS_HIRES[*, 0:Nfinite_highres_models-1]   = Lnu_xray_unabs_AGN_highres
     endif
   endif

   case strupcase(config.AGN_MODEL) of
     'SKIRTOR': begin
          out[i].LNU_AGNMOD[*, finite_model_idc]       = Lnu_AGN
          out[i].LNU_AGNMOD_UNRED[*, finite_model_idc] = Lnu_unred_AGN
          out[i].LNU_AGNMOD_HIRES[*, 0:Nfinite_highres_models-1]       = Lnu_AGN_highres
          out[i].LNU_AGNMOD_UNRED_HIRES[*, 0:Nfinite_highres_models-1] = Lnu_unred_AGN_highres
        end
     'NONE':
   endcase

   if strupcase(config.METHOD) ne 'MPFIT' then begin
     out[i].LNPROB_BESTFIT   = lnprob[bestfit_idc]
     out[i].CHI2_BESTFIT     = chi2[bestfit_idc]
     out[i].ACCEPTANCE_FRAC  = convergence_metric.ACCEPTANCE_FRAC
     out[i].ACCEPTANCE_FLAG  = convergence_metric.ACCEPTANCE_FLAG
     out[i].AUTOCORR_TIME[0:(n_elements(parameter_name)-1)] = convergence_metric.AUTOCORR_TIME
     out[i].BURN_IN_AUTOCORR = convergence_metric.BURN_IN_AUTOCORR
     out[i].CONVERGENCE_FLAG = convergence_metric.CONVERGENCE_FLAG
     out[i].SHORT_CHAIN_FLAG = short_chain_flag
     out[i].PVALUE           = pvalue

     if strupcase(config.METHOD) eq 'MCMC-AFFINE' then begin
       out[i].AUTOCORR_FLAG[0:(n_elements(parameter_name)-1)] = convergence_metric.AUTOCORR_FLAG
       out[i].STRANDED_FLAG = stranded_flag
     endif

     if strupcase(config.METHOD) eq 'MCMC-ADAPTIVE' then begin
       out[i].GELMAN_RUBIN_R_HAT[0:(n_elements(parameter_name)-1)] = convergence_metric.GELMAN_RUBIN_R_HAT
       out[i].GELMAN_RUBIN_FLAG[0:(n_elements(parameter_name)-1)]  = convergence_metric.GELMAN_RUBIN_FLAG
       out[i].BROOKS_GELMAN_R_HAT = convergence_metric.BROOKS_GELMAN_R_HAT
       out[i].BROOKS_GELMAN_FLAG  = convergence_metric.BROOKS_GELMAN_FLAG
       out[i].NCHAIN_R_HAT        = convergence_metric.NCHAIN_R_HAT
     endif
   endif else begin
     out[i].COVARIANCE       = covariance
     out[i].STATUS           = convergence_metric.STATUS
     out[i].STATUS_FLAG      = convergence_metric.STATUS_FLAG
     out[i].ERROR_MSG        = convergence_metric.ERROR_MSG
     out[i].ITER_FRAC        = convergence_metric.ITER_FRAC
     out[i].ITER_FLAG        = convergence_metric.ITER_FLAG
     out[i].STUCK_FRAC       = convergence_metric.STUCK_FRAC
     out[i].STUCK_FLAG       = convergence_metric.STUCK_FLAG
     out[i].SIMILAR_FRAC[0:(n_elements(parameter_name)-1)] = convergence_metric.SIMILAR_FRAC
     out[i].SIMILAR_FLAG[0:(n_elements(parameter_name)-1)] = convergence_metric.SIMILAR_FLAG
     out[i].NFUNC_EVALS      = convergence_metric.NFUNC_EVALS
     out[i].CONVERGENCE_FLAG = convergence_metric.CONVERGENCE_FLAG
     out[i].PVALUE           = convergence_metric.PVALUE
     out[i].DOF              = convergence_metric.DOF
   endelse

   ; If progress printing, print progress bar.
   if config.PRINT_PROGRESS then lightning_print_progress, i, Nsed, t0, funcname='LIGHTNING_POSTPROCESSING'

 endfor

 date = (strmid(timestamp(/utc), 0, 19) + 'Z').replace(':','-')
 mwrfits, out, input_dir+'lightning_output/postprocessed_data_'+date+'.fits', /create
 file_gzip, input_dir+'lightning_output/postprocessed_data_'+date+'.fits', /delete

 if ~(config.KEEP_INTERMEDIATE_OUTPUT) then file_delete, input_dir+'lightning_output/output_sav_files', /recursive

end
