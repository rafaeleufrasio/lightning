function calc_integrated_AGN_luminosity, results, config, arf
;+
; Name
; ----
;   CALC_INTEGRATED_AGN_LUMINOSITY
;
; Purpose
; -------
;   Turn the X-ray model parameters into the integrated
;   luminosity of the UV-IR component of the AGN model.
;
; Input
; -----
;   results : structure
;       An IDL structure containing the fit to J033226
;   config : structure
;       An IDL structure containing the Lightning configuration
;       for the X-ray model fit to J033226
;   arf : structure
;       An IDL structure containing the auxiliary response function
;       (ARF)
;
; Output
; ------
;   Returns an IDL structure containing the parameters and luminosity
;
;-

    compile_opt IDL2

    ; Make sure our constants are in memory
    lightning_constants

    ; Construct an array for the parameters
    parameters = [[results.TAUV_DIFF],$
                  [results.DELTA],$
                  [replicate(results.TAUV_BC, 1000)],$
                  [results.UMIN],$
                  [replicate(results.UMAX, 1000)],$
                  [replicate(results.ALPHA, 1000)],$
                  [results.GAMMA],$
                  [replicate(results.QPAH, 1000)],$
                  [results.NH],$
                  [results.AGN_MASS],$
                  [results.AGN_LOGMDOT],$
                  [replicate(results.TAU97, 1000)],$
                  [results.AGN_COSI],$
                  [transpose(results.PSI)]]

    parameters = transpose(parameters)

    models = lightning_models(config,$
                              filter_labels=results.FILTER_LABELS,$
                              redshift=results.REDSHIFT,$
                              xray_bandpass=results.XRAY_BANDPASS,$
                              xray_exposure=results.XRAY_EXPOSURE,$
                              galactic_nh=results.GALACTIC_NH,$
                              arf_E_lo=arf.ENERG_LO, $
                              arf_E_hi=arf.ENERG_HI,$
                              arf_specresp=arf.SPECRESP)

    ; For each set of AGN parameters, we calculate the L2500
    L2500 = qsosed_L2500(models.XRAY_MODELS,$
                         results.AGN_MASS,$
                         results.AGN_LOGMDOT)


    ; We call SKIRTOR_SED to get the integrated luminosity and
    ; L2500 of each template, which differ slightly with cos i
    i_AGN = 180.d0 / !dpi * acos(results.AGN_COSI)
    mean_Lnu_AGN = skirtor_sed(models.AGN_MODELS,$
                               tau97=replicate(results.TAU97, 1000),$
                               i_agn=i_AGN,$
                               tauV_diff=results.TAUV_DIFF,$
                               delta=results.DELTA,$
                               Lbol_abs_AGN=Lbol_abs_SKIRTOR,$
                               Lbol_AGN=Lbol_SKIRTOR,$
                               L2500=L2500_SKIRTOR,$
                               mean_Lnu_unred_AGN=mean_Lnu_unred_AGN,$
                               _extra=config_nopriors)

    ; And then we normalize the integrated luminosity to the
    ; L2500 of the qsosed model.
    Lbol_AGN_model = L2500 * Lbol_SKIRTOR / L2500_SKIRTOR

    AGN_model_lum = {SED_ID: results.SED_ID,$
                     REDSHIFT: results.REDSHIFT,$
                     NH: results.NH,$
                     AGN_MASS: results.AGN_MASS,$
                     AGN_LOGMDOT: results.AGN_LOGMDOT,$
                     TAU97: replicate(results.TAU97, 1000),$
                     AGN_COSI: results.AGN_COSI,$
                     LBOL_AGN_MODEL: Lbol_AGN_model}

    return, AGN_model_lum

end
