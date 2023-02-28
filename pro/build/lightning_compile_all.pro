;+
; Name
; ----
;   LIGHTNING_COMPILE_ALL
;
; Purpose
; -------
;   Compiles all of the functions and procedures required to run
;   Lightning.
;
; Calling Sequence
; ----------------
;   @lightning_compile_all
;
; Modification History
; --------------------
;   - 2022/03/17: Created (Erik B. Monson)
;   - 2022/06/27: Updated paths to files (Keith Doore)
;   - 2022/06/27: Added documentation (Keith Doore)
;   - 2022/08/11: Updated path locations and added new functions/procedures (Keith Doore)
;   - 2022/09/14: Updates to allow fitting with X-ray fluxes (Erik B. Monson)
;-


; Documentation
 .r pro/build/lightning_build_docs.pro

; Constants
 .r pro/miscellaneous/lightning_constants.pro

; Filter generator
 .r pro/filters/add_filters.pro
 .r pro/filters/get_filters.pro

; Miscellaneous
 .r pro/miscellaneous/mrandomn.pro
 .r pro/miscellaneous/percentile.pro
 .r pro/miscellaneous/print_to_width.pro
 .r pro/miscellaneous/trapez.pro
 .r pro/miscellaneous/trapex.pro
 .r pro/miscellaneous/trapow.pro
 .r pro/miscellaneous/trap_int.pro
 .r pro/miscellaneous/lightning_print_progress.pro

; Configuration
 .r pro/configuration/lightning_configure_defaults.pro
 .r pro/configuration/prior_interactive.pro
 .r pro/configuration/lightning_configure_interactive.pro
 .r pro/configuration/prior_check.pro
 .r pro/configuration/lightning_configure_check.pro
 .r pro/configuration/tabulated_prior_check.pro

; Input read
 .r pro/input/ascii_input_to_fits.pro
 .r pro/input/lightning_xray_input.pro
 .r pro/input/lightning_input.pro

; Dust emission
 .r pro/models/dust/DL07/dl07_models.pro
 .r pro/models/dust/DL07/dl07_fpdr.pro
 .r pro/models/dust/DL07/dl07_sed.pro
 .r pro/models/dust/DL07/dl07_spectrum.pro

; Calzetti00 attenuation
 .r pro/models/attenuation/calzetti00/calzetti00_atten.pro

; Doore21 attenuation
 .r pro/models/attenuation/doore21/doore21_atten.pro
 .r pro/models/attenuation/doore21/doore21_interp_lbol_abs_table.pro
 .r pro/models/attenuation/doore21/generate_doore21_lbol_abs_table.pro

; AGN emission
 .r pro/models/agn/stalevski2016/skirtor_models.pro
 .r pro/models/agn/stalevski2016/skirtor_sed.pro
 .r pro/models/agn/stalevski2016/skirtor_spectrum.pro

; Non-parametric Stellar emission
 .r pro/models/stellar/nonparametric_sfh/binned_stellar_models.pro
 .r pro/models/stellar/nonparametric_sfh/binned_stellar_sed.pro
 .r pro/models/stellar/nonparametric_sfh/binned_stellar_spectrum.pro

; Instantaneous Stellar emission
 .r pro/models/stellar/instantaneous/instantaneous_stellar_models.pro
 .r pro/models/stellar/instantaneous/instantaneous_stellar_sed.pro
 .r pro/models/stellar/instantaneous/instantaneous_stellar_spectrum.pro

; Xray emission
 .r pro/models/xray/L2keV_J07.pro
 .r pro/models/xray/L2keV_LR17.pro
 .r pro/models/xray/gilbertson22_LX_tau.pro
 .r pro/models/xray/interpolate_arf.pro
 .r pro/models/xray/xray_plaw_expcut.pro
 .r pro/models/xray/load_xray_abs.pro
 ; QSOSED emission
 .r pro/models/xray/qsosed/qsosed_models.pro
 .r pro/models/xray/qsosed/qsosed_L2500.pro
 .r pro/models/xray/qsosed/qsosed_spectrum.pro
 .r pro/models/xray/xrb_xagn_models.pro

; Model data
 .r pro/models/lightning_models.pro
 .r pro/models/lightning_model_counts.pro
 .r pro/models/lightning_model_lnu.pro
 .r pro/models/lightning_model_lnu_highres.pro
 .r pro/models/lightning_model_lnuxray.pro
 .r pro/models/lightning_model_lnuxray_highres.pro
 .r pro/models/lightning_model_lnprob.pro

; Priors
 .r pro/fitting/priors/generate_prior_struct.pro
 .r pro/fitting/priors/prior_sampled_initialization.pro
 .r pro/fitting/priors/lightning_priors.pro

; MCMC fitting algorithms
 .r pro/fitting/algorithms/mcmc/gw_proposal_z.pro
 .r pro/fitting/algorithms/mcmc/gw_stretch_move.pro
 .r pro/fitting/algorithms/mcmc/autocorr_time.pro
 .r pro/fitting/algorithms/mcmc/brooks_gelman_stat.pro
 .r pro/fitting/algorithms/mcmc/gelman_rubin_stat.pro
 .r pro/fitting/algorithms/mcmc/mcmc_convergence.pro
 .r pro/fitting/algorithms/mcmc/lightning_mcmc.pro

; MPFIT fitting algorithms
 .r pro/fitting/algorithms/mpfit/lightning_mpfit_function.pro
 .r pro/fitting/algorithms/mpfit/mpfit_convergence.pro
 .r pro/fitting/algorithms/mpfit/lightning_mpfit.pro

; Fit initialization
 .r pro/fitting/lightning_fit.pro

; Post-processing
 .r pro/postprocessing/ppc.pro
 .r pro/postprocessing/lightning_postprocessing.pro

; Main procedure
 .r pro/lightning.pro
