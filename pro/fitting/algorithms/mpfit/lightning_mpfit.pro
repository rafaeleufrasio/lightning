pro lightning_mpfit, input_dir, sed_data, config_nopriors, models, priors
;+
; Name
; ----
;   LIGHTNING_MPFIT
;
; Purpose
; -------
;   Fits an observed SED with the Lightning models using the IDL MPFIT package,
;   which utilizes the Levenberg-Marquardt gradient decent algorithm. Initial
;   parameter positions are randomly selected from the specified prior distribution.
;
; Calling Sequence
; ----------------
;   ::
;
;       lightning_mpfit, input_dir, sed_data, config_nopriors, models, priors
;
; Inputs
; ------
;   ``input_dir`` : string scalar
;       The path to the file containing the input SED data.
;   ``sed_data`` : structure
;       A structure containing the SED luminosities and uncertainties, filter
;       labels, distances, redshifts, and optional X-ray data. (See 
;       ``lightning_input.pro`` for details and contents.)
;   ``config_nopriors`` : structure
;       A Lightning configuration structure edited to remove the prior 
;       substructures for each parameter. (See ``lightning_configure_defaults.pro``
;       for details and contents.)
;   ``models`` : structure
;       A structure containing each model structure (stellar, dust, AGN, 
;       X-ray) as a substructure. (See ``lightning_models.pro`` for details
;       and contents.)
;   ``priors`` : structure
;        A structure containing the prior hyper-parameters. (See
;        ``generate_prior_struct.pro`` for details and contents.)
;
; Output
; ------
;   An IDL save file saved to ``<input_dir>/lightning_output/output_sav_files/`` named 
;   ``lightning_output_<galaxy_id>.sav'``, containing the resulting MPFIT best fit
;   parameters and log probability for each set of parameter starting points, and the
;   convergence metrics of the fits.
;
; References
; ----------
;   - `Markwardt, C. B. 2009, ASP Conf. Ser. 411, 251 <https://ui.adsabs.harvard.edu/abs/2009ASPC..411..251M/abstract>`_
;   - `Moré, J. 1978, Numerical Analysis, 630, 105 <https://ui.adsabs.harvard.edu/abs/1978LNM...630..105M/abstract>`_
;
; Modification History
; --------------------
;   - 2022/01/01: Created (Rafael Eufrasio)
;   - 2022/08/15: Major update to include new implementation (e.g., prior, config, etc.) (Keith Doore)
;   - 2022/08/18: Added progress printing (Keith Doore)
;   - 2023/01/23: Fixed bug where if ``status <= 0`` then the corresponding ``lnprob_mpfit = 0``.
;     Now have ``lnprob_mpfit = NaN`` if ``status <= 0`` (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(input_dir) eq 0 then message, 'Variable is undefined: INPUT_DIR.'
 if size(input_dir, /type) ne 7 then message, 'INPUT_DIR must be of type string.'
 if size(input_dir, /n_dim) ne 0 then message, 'INPUT_DIR must be a scalar.'
 if ~file_test(input_dir ,/dir) then message, 'INPUT_DIR is not a valid directory.'

 if n_elements(sed_data) eq 0 then message, 'Variable is undefined: SED_DATA.'
 if size(sed_data, /type) ne 8 then message, 'SED_DATA must be of type structure.'

 if n_elements(config_nopriors) eq 0 then message, 'Variable is undefined: CONFIG_NOPRIORS.'
 if size(config_nopriors, /type) ne 8 then message, 'CONFIG_NOPRIORS must be of type structure.'

 if n_elements(models) eq 0 then message, 'Variable is undefined: MODELS.'
 if size(models, /type) ne 8 then message, 'MODELS must be of type structure.'

 if n_elements(priors) eq 0 then message, 'Variable is undefined: PRIORS.'
 if size(priors, /type) ne 8 then message, 'PRIORS must be of type structure.'


; Generate starting parameter values and proposal distribution
 parameter_start = prior_sampled_initialization(priors, config_nopriors.NSOLVERS)


; Set up MPFIT parinfo
 Nparam = n_elements(priors.parameter_name)
 parinfo = replicate({fixed: 0, limited: [1, 1], limits: [0.d, 0.d]}, Nparam)
 ; Set fixed values flag
 ;     .FIXED - a boolean value, whether the parameter is to be held
 ;              fixed or not.  Fixed parameters are not varied by
 ;              MPFIT, but are passed on to MYFUNCT for evaluation.
 parinfo[priors.idcs.FIXED].fixed = 1
 ; Limited flags are always set as Lightning requires finite bounds on parameters
 ; Set limiting values
 ;     .LIMITS - a two-element float or double array.  Gives the
 ;               parameter limits on the lower and upper sides,
 ;               respectively.  Zero, one or two of these values can be
 ;               set, depending on the values of LIMITED.  Both LIMITED
 ;               and LIMITS must be given together.
 for i=0, Nparam-1 do parinfo[i].limits = [(priors.MIN_BOUND)[i], (priors.MAX_BOUND)[i]]


; Fit using MPFIT 
 status           = intarr(config_nopriors.NSOLVERS)
 error_msg        = strarr(config_nopriors.NSOLVERS)
 parameters_mpfit = dblarr(Nparam, config_nopriors.NSOLVERS)
 parameters_error = dblarr(Nparam, config_nopriors.NSOLVERS)
 covariance       = dblarr(Nparam, Nparam, config_nopriors.NSOLVERS)
 lnprob_mpfit     = replicate(!values.d_NaN, config_nopriors.NSOLVERS)
 Niterations      = lonarr(config_nopriors.NSOLVERS)
 Nfunc_evals      = lonarr(config_nopriors.NSOLVERS)
 DoF              = 0L


 t0 = systime(/sec)
 for i=0, (config_nopriors.NSOLVERS-1) do begin

   fcnargs = {sed_data:sed_data, config_nopriors:config_nopriors, models:models, priors:priors}

   params = mpfit('lightning_mpfit_function', parameter_start[*, i], functargs=fcnargs, /quiet, $
                  status=status_flag, errmsg=errmsg, covar=covar, parinfo=parinfo, perror=perror, $
                  Niter=Niter, bestnorm=chi2_bestfit, Nfev=Nfev, dof=deg_o_free, _extra=config_nopriors) 

   status[i]    = status_flag
   error_msg[i] = errmsg
   if status_flag gt 0 then begin
     parameters_mpfit[*, i]  = params
     parameters_error[*,  i] = perror
     covariance[*, *, i]     = covar
     lnprob_mpfit[i]         = -0.5 * chi2_bestfit
     Niterations[i]          = niter
     Nfunc_evals[i]          = Nfev
     DoF                     = deg_o_free
   endif else begin
     ; If run fails, then set parameters to the starting value to not cause issues in post-processing
     parameters_mpfit[*, i]  = parameter_start[*, i]
   endelse

   ; If progress printing, print progress bar.
   if config_nopriors.PRINT_PROGRESS then lightning_print_progress, i, config_nopriors.NSOLVERS, $
                                                                    t0, funcname='LIGHTNING_MPFIT'
 endfor

; Test for quality of fits
 convergence_metric = mpfit_convergence(parameters_mpfit, lnprob_mpfit, status, error_msg, Niterations, $
                                        Nfunc_evals, DoF, _extra=config_nopriors)

; Save output
 parameter_name = priors.parameter_name
 save, parameters_mpfit, parameters_error, covariance, lnprob_mpfit, models, parameter_name, convergence_metric, $
       /compress, filename=input_dir+'lightning_output/output_sav_files/lightning_output_'+sed_data.SED_ID+'.sav'

end