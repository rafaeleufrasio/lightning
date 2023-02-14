function sed_plots, results

    compile_opt IDL2

    lightning_constants

    ; Low-res input and models
    wave_obs = results.WAVE_FILTERS
    nu_obs = !lightning_cgs.clight * 1e4 / wave_obs
    lnu_obs = results.LNU_OBS
    lnu_unc = results.LNU_UNC
    help, lnu_unc
    minchi2 = min(results.CHI2, bestfit)
    lnu_mod = results.LNU_MOD[*,bestfit]
    model_unc_vector = rebin(reform(results.MODEL_UNC, 1, n_elements(results)), n_elements(results[0].WAVE_FILTERS), n_elements(results))
    help, model_unc_vector
    lnu_unc_total = sqrt(lnu_unc^2 + (model_unc_vector * lnu_mod)^2)
    help, lnu_unc_total
    delchi = (lnu_obs - lnu_mod) / lnu_unc_total

    ; High res models
    wave_hires = results.WAVE_HIRES
    nu_hires = !lightning_cgs.clight * 1e4 / wave_hires
    lnu_total_hires = results.LNU_MOD_HIRES
    ; lnu_stars_hires = results.LNU_STARMOD_HIRES
    ; lnu_AGN_hires = results.LNU_AGNMOD_HIRES
    ; lnu_dust_hires = results.LNU_DUSTMOD_HIRES

    lines = []

    p1 = plot(wave_hires[*,0], nu_hires[*,0] * lnu_total_hires[*,0],$
              /XLOG,$
              /YLOG,$
              thick=2,$
              color='green')
    lines = [lines, p1]

    p1 = plot(wave_hires[*,1], nu_hires[*,1] * lnu_total_hires[*,1],$
              /OVERPLOT,$
              /XLOG,$
              /YLOG,$
              thick=2,$
              color='orange')
    lines = [lines, p1]

    p1 = errorplot(wave_obs[*,0], nu_obs[*,0] * lnu_obs[*,0],$
                   nu_obs[*,0] * lnu_unc_total[*,0],$
                   /OVERPLOT,$
                   /XLOG,$
                   /YLOG,$
                   linestyle='',$
                   sym='o',$
                   sym_color='k',$
                   /sym_filled)
    lines = [lines, p1]

    p1 = errorplot(wave_obs[*,1], nu_obs[*,1] * lnu_obs[*,1],$
                   nu_obs[*,1] * lnu_unc_total[*,1],$
                   /OVERPLOT,$
                   /XLOG,$
                   /YLOG,$
                   linestyle='',$
                   sym='o',$
                   sym_color='gray',$
                   /sym_filled)
    lines = [lines, p1]

    p1.xrange = [0.09, 1000]
    ;p1.yrange = [1e10, 3e12]

    ;p1.xshowtext = 0
    p1.ytitle = '$\nuL_{\nu}$ $[L_{\odot}]$'
    p1.xtitle = 'Observed-Frame Wavelength $[\mu m]$'

    return, p1

end
