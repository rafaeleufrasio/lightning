pro lightning_fit, input_dir, sed_data, config
;+
; Name
; ----
;   LIGHTNING_FIT
;
; Purpose
; -------
;   Generates the models and priors specified in the Lightning configuration
;   structure, and passes them and the SED data to the specified fitting
;   algorithm.
;
; Calling Sequence
; ----------------
;   ::
;
;       lightning_fit, input_dir, sed_data, config
;
; Inputs
; ------
;   ``input_dir`` : string scalar
;       The path to the file containing the input SED data.
;   ``sed_data`` : structure
;       A structure containing the SED luminosities and uncertainties, filter
;       labels, distances, redshifts, and optional X-ray data. (See
;       ``lightning_input.pro`` for details and contents.)
;   ``config`` : structure
;       A Lightning configuration structure. (See
;       ``lightning_configure_defaults.pro`` for details and contents.)
;
; Modification History
; --------------------
;   - 2022/05/13: Created (Keith Doore)
;   - 2022/07/11: Moved models generation to separate function (Keith Doore)
;   - 2022/07/27: Added non-parametric SFH bin cutting (Keith Doore)
;   - 2022/08/09: Added ``galactic_nh`` as input for xray emission from sed_data (Keith Doore)
;   - 2022/08/17: Fixed issue where ``config.PSI`` was being updated and returned to
;     ``lightning.pro`` if ``config.MAXCPU = 1`` (Keith Doore)
;   - 2022/09/15: Removed different calls to ``lightning_models.pro`` for X-ray emission and
;     no X-ray emission, since ``lightning_input.pro`` now handles including X-ray emission if
;     requested in configuration. (Keith Doore)
;   - 2022/10/25: Renamed SPS to SSP (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(input_dir) eq 0 then message, 'Variable is undefined: INPUT_DIR.'
 if size(input_dir, /type) ne 7 then message, 'INPUT_DIR is not of type string.'
 if size(input_dir, /n_dim) ne 0 then message, 'INPUT_DIR must be a scalar.'
 if ~file_test(input_dir ,/dir) then message, 'INPUT_DIR is not a valid directory.'

 if n_elements(sed_data) eq 0 then message, 'Variable is undefined: SED_DATA.'
 if size(sed_data, /type) ne 8 then message, 'SED_DATA is not of type structure.'

 if n_elements(config) eq 0 then message, 'Variable is undefined: CONFIG.'
 if size(config, /type) ne 8 then message, 'CONFIG is not of type structure.'



;====== Generate models =========
 models = lightning_models(config, input_dir=input_dir, _EXTRA=sed_data)


; Check if a stellar model is used and adjust PSI in configuration structure to reflect
;   potential changes in STEPS_BOUNDS due to redshift
 config_updated = config
 if strupcase(config.SSP) eq 'PEGASE' then begin
   if strupcase(config.SFH) eq 'NON-PARAMETRIC' then begin
     if n_elements(models.stellar_models.BOUNDS) ne n_elements(config.STEPS_BOUNDS) then begin
       Nsteps_new = n_elements(models.stellar_models.BOUNDS) - 1
       config_hash = orderedhash(config, /extract)
       config_hash['PSI', 'PRIOR'] = config_hash['PSI', 'PRIOR', 0:Nsteps_new-1]
       config_hash['PSI', 'PRIOR_ARG'] = config_hash['PSI', 'PRIOR_ARG', 0:Nsteps_new-1, *]
       config_updated = config_hash.ToStruct(/recursive)
     endif
   endif
 endif


;======= Generate Prior structure =======
 priors = generate_prior_struct(config_updated, sed_data.SED_ID, config_nopriors=config_nopriors)


;======= Fit the SED =========
 case strupcase(config.METHOD) of
   ;GRID':          lightning_grid,      input_dir, sed_data, config_nopriors, models, priors
   ;INVERSION':     lightning_inversion, input_dir, sed_data, config_nopriors, models, priors
   'MCMC-ADAPTIVE': lightning_mcmc,      input_dir, sed_data, config_nopriors, models, priors
   'MCMC-AFFINE':   lightning_mcmc,      input_dir, sed_data, config_nopriors, models, priors
   'MPFIT':         lightning_mpfit,     input_dir, sed_data, config_nopriors, models, priors
 endcase

end
