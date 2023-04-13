function posterior_comparison, xray_res, noxray_res, agn_model_lum
;+
; Name
; ----
;   POSTERIOR_COMPARISON
;
; Purpose
; -------
;   Create a cornerplot figure comparing the posterior distributions
;   for the AGN parameters of J033226.49-274035.5 when we fit with
;   and without the X-ray model.
;
; Input
; -----
;   xray_res : structure
;       An IDL structure containing the fit to J033226 with
;       the X-ray model
;   noxray_res : structure
;       An IDL structure containing the fit to J033226 with
;       the X-ray model
;   agn_model_lum : structure
;       An IDL structure containing the integrated UV-IR AGN
;       luminosity for the fit with the X-ray model
;
; Output
; ------
;   Returns an IDL graphics object containing the plot
;
; Dependencies
; ------------
;   Keith's corner plot routine from lightning-visualization
;
;-

    compile_opt IDL2

    log_L_AGN_x = alog10(agn_model_lum.LBOL_AGN_MODEL)

    ; Build the array for the corner plot function
    x_arr = transpose([[reform(xray_res.AGN_MASS) / 1.d8],$
                       [reform(xray_res.AGN_LOGMDOT)],$
                       [reform(xray_res.NH) / 1e2],$
                       [log_L_AGN_x],$
                       [reform(xray_res.AGN_COSI)]])

    nox_arr = transpose([[randomn(seed, n_elements(xray_res.AGN_MASS))/1000],$
                         [randomn(seed, n_elements(xray_res.AGN_MASS))/1000],$
                         [-1 + randomn(seed, n_elements(xray_res.AGN_MASS))/1000],$
                         [reform(noxray_res.LOG_L_AGN)],$
                         [reform(noxray_res.AGN_COSI)]])

    full_arr = [[[nox_arr]], [[x_arr]]]

    labels = ['$M_{\rm SMBH}$ $[10^8$ $M_{\odot}]$',$
              'log !sm!r!A .!n!3 ',$
              '$N_H$ $[10^{22}$ $cm^{-2}]$',$
              '$log_{10}(L_{\rm AGN}$ $[L_{\odot}]$)',$
              'cos $i_{AGN}$']

    range_arr = [[9, 13],$
                 [-1.5, -1.38],$
                 [0.0, 0.25],$
                 [11.9, 12.25],$
                 [0.75, 0.99]]

    p = window(dimension=[860, 800])
    p = corner_plot(full_arr, labels,$
                    distribution_range=range_arr,$
                    /NORMALIZE,$
                    /SHOW_MEDIAN,$
                    CONTOUR_LEVELS=[0.68, 0.95],$
                    CONTOUR_SMOOTH=3,$
                    TICKINTERVAL=[2, 0.05, 0.1, 0.2, 0.1],$
                    DISTRIBUTION_COLOR=['orange', 'blue'],$
                    DISTRIBUTION_THICK=2,$
                    CONTOUR_THICK=2,$
                    font_size=10,$
                    position=[0.08, 0.08, 0.92, 0.98],$
                    /current)

    t1 = text(0.75, 0.82, 'J033226.49-274035.5', /CURRENT, /NORMAL, font_size=14, alignment=0.5)

    x = objarr(2)
    x[0] = plot([0, 0], [0, 0], color='blue', thick=2, name='with qsosed X-ray model', position=[2, 2, 3, 3], /current)
    x[1] = plot([0, 0], [0, 0], color='orange', thick=2, name='without X-ray model', position=[2, 2, 3, 3], /current)
    lgnd = legend(target=x, transparency=30, position=[0.75, 0.80], $
                  vertical_alignment=1, horizontal_alignment=0.5, font_size=14)

    return, p

end
