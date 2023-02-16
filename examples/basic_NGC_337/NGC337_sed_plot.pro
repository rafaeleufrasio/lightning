function NGC337_sed_plot, results
;+
; Name
; ----
;   NGC337_SED_PLOT
;
; Purpose
; -------
;   Create a plot showing the best-fit SED model for NGC 337
;
; Input
; -----
;   results : structure
;       The IDL structure containing the fit to NGC 337.
;
; Output
; ------
;   Returns the figure as an IDL graphics object.
;
;-


    compile_opt IDL2

    lightning_constants

    ; Low-res input and models
    wave_obs = results.WAVE_FILTERS
    nu_obs = !lightning_cgs.clight * 1e4 / wave_obs
    lnu_obs = results.LNU_OBS
    lnu_unc = results.LNU_UNC
    minchi2 = min(results.CHI2, bestfit)
    lnu_mod = results.LNU_MOD[*,bestfit]
    model_unc_vector = rebin(reform(results.MODEL_UNC, 1, n_elements(results)), n_elements(results[0].WAVE_FILTERS), n_elements(results))
    lnu_unc_total = sqrt(lnu_unc^2 + (model_unc_vector * lnu_mod)^2)
    delchi = (lnu_obs - lnu_mod) / lnu_unc_total

    ; High res models
    wave_hires = results.WAVE_HIRES
    nu_hires = !lightning_cgs.clight * 1e4 / wave_hires
    lnu_total_hires = results.LNU_MOD_HIRES

    ; We could easily plot the other components of the fit;
    ; we omit them here for brevity.
    ; lnu_stars_hires = results.LNU_STARMOD_HIRES
    ; lnu_AGN_hires = results.LNU_AGNMOD_HIRES
    ; lnu_dust_hires = results.LNU_DUSTMOD_HIRES

    lines = []

    p1 = plot(wave_hires, nu_hires * lnu_total_hires,$
              /XLOG,$
              /YLOG,$
              thick=2,$
              color='gray',$
              name='Total')
    lines = [lines, p1]

    p1 = errorplot(wave_obs, nu_obs * lnu_obs,$
                   nu_obs * lnu_unc_total,$
                   /OVERPLOT,$
                   /XLOG,$
                   /YLOG,$
                   linestyle='',$
                   sym='o',$
                   sym_color='k',$
                   /sym_filled,$
                   name='Data')
    lines = [lines, p1]

    p1.xrange = [0.09, 1000]
    p1.yrange = [1e7, 5e10]

    p1.xshowtext = 0
    p1.ytitle = '$\nuL_{\nu}$ $[L_{\odot}]$'
    ;p1.xtitle = 'Observed-Frame Wavelength $[\mu m]$'

    position = p1.position
    figwidth = position[2] - position[0]
    figheight = position[3] - position[1]
    residual_position = [position[0],$
                         position[1],$
                         position[2],$
                         position[1] + figheight * 4/15]
    position = [position[0],$
                position[1] + figheight * 4/15 + 0.02,$
                position[2],$
                position[3]]

    p1.position = position

    l1 = legend(target=lines, transp=30, POSITION=[position[0]+0.01, position[1]+0.01], $
                VERTICAL_ALIGNMENT=0, HORIZONTAL_ALIGNMENT=0,$
                font_size=8)

    t1 = text(position[2] - 0.02, position[3] - 0.03, results.SED_ID, VERTICAL_ALIGNMENT=1, ALIGNMENT=1)
    t2 = text(position[2] - 0.02, position[3] - 0.06, 'p = ' + strtrim(string(results.PVALUE, format='(F5.2)'),2), VERTICAL_ALIGNMENT=1, ALIGNMENT=1)

    p2 = plot(p1.xrange, [0, 0],$
              /CURRENT,$
              /XLOG,$
              xtickunits='scientific',$
              position=residual_position, $
              xrange=p1.xrange,$
              color='gray',$
              yrange=[-4.5, 4.5],$
              thick=2)
    p2 = errorplot(wave_obs, delchi, 1 + dblarr(n_elements(delchi)),$
                   /OVERPLOT,$
                   /XLOG,$
                   linestyle='',$
                   sym='o',$
                   sym_color='k',$
                   /sym_filled)

    p2.xtitle = 'Observed-Frame Wavelength $[\mu m]$'
    p2.ytitle = 'Residual [$\sigma$]'

    return, p2

end
