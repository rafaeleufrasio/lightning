pro lightning_mcmc, input_dir, sed_data, config_nopriors, models, priors
;+
; Name
; ----
;   LIGHTNING_MCMC
;
; Purpose
; -------
;   Fits an observed SED with the Lightning models using either the vanishing
;   adaptive Metropolis-Hastings algorithm (Algorithm 4 from Andrieu & Thoms 2008),
;   or the affine-invariant Goodman & Weare (2010) algorithm. Initial parameter
;   positions are randomly selected from the specified prior distribution.
;
; Calling Sequence
; ----------------
;   ::
;
;       lightning_mcmc, input_dir, sed_data, config_nopriors, models, priors
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
;   ``lightning_output_<galaxy_id>.sav'``, containing the resulting MCMC parameter
;   and log probability chain(s), and the convergence metrics of the chain(s).
;
; Notes
; -----
;   The Goodman & Weare affine invariant algorithm requires an ensemble of *at least* ``Ndim + 1`` walkers, which all
;   must be randomly initialized such that the ``[Ndim, Nparallel]`` matrix describing their positions is
;   nonsingular. These walkers are not independent and should not be treated as such in postprocessing: convergence
;   metrics such as Brooks-Gelman and Gelman-Rubin tests are not applicable, and the chains should be flattened into a
;   single vector to derive the posterior distributions. We recommend using the acceptance fraction (ideally 20-50%)
;   and the autocorrelation time (ideally much shorter than the length of the chains for all dimensions)
;   to assess the performance of this method.
;
; References
; ----------
;   - `Andrieu, C., & Thoms, J. 2008, Statistics and Computing, 18, 30 <https://link.springer.com/article/10.1007/s11222-008-9110-y>`_
;   - `Goodman, J., & Weare, J. 2010, CAMCS, 5, 65 <https://ui.adsabs.harvard.edu/abs/2010CAMCS...5...65G/abstract>`_
;
; Modification History
; --------------------
;   - 2020/04/27: Replaced if statements with ``n_elements`` on keywords to use ``keyword_set`` (Keith Doore)
;   - 2020/05/06: Set default starting attenuation parameters to ``0.0`` so if not wanted they would not be used (Keith Doore)
;   - 2020/05/06: Added ability to use pure Calzetti curve (no changes here due to ``_extra``) (Keith Doore)
;   - 2020/05/06: Added ``_ref_extra`` for keyword inheritance to cut down on list of keywords (Keith Doore)
;   - 2020/05/06: Changed all keywords that were the same to match across functions/procedures (Keith Doore)
;   - 2020/05/06: Removed any repetitive items at beginning that are set by other functions if not set in MCMC call (Keith Doore)
;   - 2020/05/06: Added needed items to run Tuffs attenuation as to match other attenuation (Keith Doore)
;   - 2021/03/17: Added UV-to-IR AGN fitting (Erik Monson)
;   - 2021/04/13: Modified Tuffs attenuation to have ``rdisk`` as intrinsic property and ``rbulge`` as B/D (Keith Doore)
;   - 2021/04/16: Added X-ray fitting (Erik Monson)
;   - 2021/04/16: Change sampling of AGN angle parameters: ``i`` -> ``cosi``; ``oa`` -> ``sin(oa)`` (Erik Monson)
;   - 2022/04/23: Changed ``lightning_mcmc.pro`` to use the new modular Lightning framework (Erik Monson)
;   - 2022/04/23: Added documentation (Erik Monson)
;   - 2022/06/07: Major update to include new implementation (e.g., prior, config, etc.) (Keith Doore)
;   - 2022/07/05: Renamed ``config`` to ``config_nopriors`` to reflect removed prior substructures (Keith Doore)
;   - 2022/07/06: Changed how ``sigma_new`` is initialized to eliminate for loop (Keith Doore)
;   - 2022/07/19: Added compression to output MCMC sav file to reduce memory usage (Keith Doore)
;   - 2022/07/19: Updated how ``mu_new`` is initialized to prevent bugs if first trial is rejected (Keith Doore)
;   - 2022/08/18: Added progress printing (Keith Doore)
;   - 2022/09/01: Added ``chi2_chain`` to output file in addition to ``lnprob_chain`` (Erik B. Monson)
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
 parameter_start = prior_sampled_initialization(priors, config_nopriors.Nparallel, sigma_start=sigma_start)



; Get inital log prob
 lnprob_old = lightning_model_lnprob(sed_data, parameter_start, config_nopriors, models, priors, UVIR_chi2=UVIR_chi2_old, xray_chi2=xray_chi2_old)
 chi2_old = UVIR_chi2_old + xray_chi2_old

; Initialize chains
 Nparam = n_elements(priors.parameter_name)
 chain = dblarr(Nparam, config_nopriors.Ntrials, config_nopriors.Nparallel)
 lnprob_chain = dblarr(config_nopriors.Ntrials, config_nopriors.Nparallel)
 chi2_chain = dblarr(config_nopriors.Ntrials, config_nopriors.Nparallel)
 accepted_trials = dblarr(config_nopriors.Nparallel)

 chain_old  = parameter_start

 chain[*,0,*] = chain_old
 lnprob_chain[0,*] = lnprob_old
 chi2_chain[0,*] = chi2_old

 variable_param = priors.idcs.variable
 nvar = n_elements(variable_param)

 if strupcase(config_nopriors.METHOD) eq 'MCMC-ADAPTIVE' then begin
   lambda_new = 2.38^2.d0 / replicate(double(nvar), config_nopriors.Nparallel)
   ; Need to offset mu from the chain, otherwise issues will occur if first iteration is rejected.
   ;   We randomly offset it, to make it a noisy estimate of the prior mean
   mu_new = chain_old + rebin(sigma_start, Nparam, config_nopriors.Nparallel) * randomn(seed, Nparam,config_nopriors.Nparallel)

   ; Sigma must only include variable parameters due to not needing to sample fixed parameters
   sigma_new=rebin(diag_matrix(sigma_start[variable_param]), nvar, nvar, config_nopriors.Nparallel)

   ; Optimal acceptance probabilities for tuning the proposal covariance.
   ; From Gelman et al. 1996 -- http://stat.columbia.edu/~gelman/research/published/baystat5.pdf
   p_jump = [0.441, 0.352, 0.316, 0.279, 0.275, 0.266, 0.261, 0.255, 0.261, 0.267]
   alpha_star = p_jump[(nvar - 1) < (n_elements(p_jump) - 1)]
 endif


; Start MCMC loop
 t0 = systime(/seconds)
 for trial = 0L, (config_nopriors.Ntrials - 2) do begin

   chain_new = chain_old

   case strupcase(config_nopriors.METHOD) of
     'MCMC-ADAPTIVE': begin
          gamma_new = 1.d0 / (trial + 1.d0)^config_nopriors.BETA_EXPONENT
          lambda_old = lambda_new
          mu_old = mu_new
          sigma_old = sigma_new
          covar = rebin(reform(lambda_old, 1, 1, config_nopriors.Nparallel), nvar, nvar, config_nopriors.Nparallel) * sigma_old

          for njump = 0, (config_nopriors.Nparallel - 1) do begin
            jump = mrandomn(seed, covar[*, *, njump], status=status)
            while (status gt 0) do begin
              ; Due to floating point rounding error, the covariance matrix constructed by the
              ; adpative MCMC is sometimes nonsymmetric. We try to remedy that here.
              covar[*, *, njump] = (covar[*, *, njump] + transpose(covar[*, *, njump])) / 2.d0
              ;evals = eigenql(covar[*,*,njump], EIGENVECTORS = evecs) ; Even after the above step covar is sometimes not symmetric enough for eigenql
              evals = la_eigenql(covar[*, *, njump], eigenvectors=evecs, status=status_ql)
              ; If the covariance is still too nonsymmetric to use a QL decomposition,
              ; we fall back on an eigenproblem solver for nonsymmetric matrices. It's about 20us slower.
              if (status_ql gt 0) then evals = la_eigenproblem(covar[*, *, njump], EIGENVECTORS=evecs)
              ; Pin vanishing or negative eigenvals to a small finite value
              evals[where(evals le 1.0d-16, /null)] = 1.0d-16
              ; Reconstruct covariance matrix from diagonal representation and force symmetry
              ; by averaging it with its transpose.
              covar[*, *, njump] = evecs # diag_matrix(evals) # invert(evecs)
              covar[*, *, njump] = (covar[*, *, njump] + transpose(covar[*, *, njump])) / 2.d0
              jump = mrandomn(seed, covar[*, *, njump], status=status)
            endwhile

            ; Apply initial jump to old chain
            chain_new[variable_param, njump] = chain_old[variable_param, njump] + jump

          endfor
        end

     'MCMC-AFFINE': begin
          ; Get a new ensemble state using the Goodman-Weare algorithm with the stretch move
          chain_new[variable_param, *] = gw_stretch_move(chain_old[variable_param, *], config_nopriors.AFFINE_A, z=z)
        end
   endcase


   ; Get new log prob
   lnprob_new = lightning_model_lnprob(sed_data, chain_new, config_nopriors, models, priors, UVIR_chi2=UVIR_chi2_new, xray_chi2=xray_chi2_new)
   chi2_new = UVIR_chi2_new + xray_chi2_new

   ; Calculate the acceptance probability
   ; a = p(new) / p(old)
   ; ln(a) = ln(p(new)) - ln(p(old))
   ; a = exp[ ln(p(new)) - ln(p(old)) ]
   metr_ratio = exp((lnprob_new - lnprob_old))

   ; In the GW10 algorithm with the stretch move, the acceptance probability
   ; depends on the scaling constant z that controls the size of the move as
   ; p propto z^(N - 1)
   if strupcase(config_nopriors.METHOD) eq 'MCMC-AFFINE' then $
     metr_ratio = metr_ratio * z^(double(nvar) - 1)

   accept = where(metr_ratio gt randomu(seed, config_nopriors.Nparallel), /null, compl=reject)

   if n_elements(accept) ne 0 then begin
     chain[*, trial + 1, accept] = chain_new[*, accept]
     chain_old[*, accept] = chain_new[*, accept]
     lnprob_chain[trial + 1,accept] = lnprob_new[accept]
     lnprob_old[accept] = lnprob_new[accept]
     chi2_chain[trial + 1, accept] = chi2_new[accept]
     chi2_old[accept] = chi2_new[accept]
     accepted_trials[accept]++
   endif

   if n_elements(reject) ne 0 then begin
     chain[*, trial + 1, reject] = chain_old[*, reject]
     chain_old[*, reject] = chain_old[*, reject]
     lnprob_chain[trial + 1, reject] = lnprob_old[reject]
     lnprob_old[reject] = lnprob_old[reject]
     chi2_chain[trial + 1, reject] = chi2_old[reject]
     chi2_old[reject] = chi2_old[reject]
   endif

   ; Adjust the covariance matrix using a vanishing adaptive algorithm.
   ; See Algorithm 4 in Andrieu 2008
   if strupcase(config_nopriors.METHOD) eq 'MCMC-ADAPTIVE' then begin
     lambda_new = 10.d0^(alog10(lambda_old) + gamma_new * (min([[replicate(1.d0, config_nopriors.Nparallel)], [metr_ratio]], dim=2) - alpha_star))
     mu_new = mu_old + gamma_new * (chain_old - mu_old)
     sigma_temp = dblarr(nvar, nvar, config_nopriors.Nparallel)
     for nsig = 0,(config_nopriors.Nparallel - 1) do sigma_temp[*, *, nsig] = (chain_old[variable_param, nsig] - mu_old[variable_param, nsig]) # $
                                                                     (chain_old[variable_param, nsig] - mu_old[variable_param, nsig])
     sigma_new = sigma_old + gamma_new * (sigma_temp - sigma_old)
   endif

   ; If progress printing, print progress bar.
   if config_nopriors.PRINT_PROGRESS then lightning_print_progress, trial+1, config_nopriors.Ntrials, $
                                                                    t0, funcname='LIGHTNING_MCMC'
 endfor


; Check final chain for convergence
 convergence_metric = mcmc_convergence(chain, lnprob_chain, accepted_trials, _extra=config_nopriors)


 parameter_name = priors.parameter_name
 save, chain, lnprob_chain, chi2_chain, models, parameter_name, convergence_metric, /compress, $
       filename=input_dir+'lightning_output/output_sav_files/lightning_output_'+sed_data.SED_ID+'.sav'

end
