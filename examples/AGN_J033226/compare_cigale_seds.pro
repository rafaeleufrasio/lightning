function compare_cigale_seds, lghtng, cigale_sed, cigale_obs, cigale_results
;+
; Name
; ----
;   COMPARE_CIGALE_SEDS
;
; Purpose
; -------
;   Generates SED plots of J033226. Input tables
;   must be the post-processed outputs from Lightning and CIGALE.
;
; Calling Sequence
; ----------------
;   ::
;
;       fig = compare_cigale_seds(lghtng, cigale_sed, cigale_obs, cigale_results)
;
; Input
; -----
;   ``lghtng`` : structure
;       The Lightning post-processed output.
;   ``cigale_sed`` : structure
;       The CIGALE best-fit SED model output.
;   ``cigale_obs`` : structure
;       The CIGALE input observational data.
;   ``cigale_results`` : structure
;       The CIGALE results.
;
; Output
; ------
;   ``fig`` : object
;       The plot object containing the figure.
;
; Modification History
; --------------------
;   - 2023/04/12: Created (Keith Doore)
;-
 On_error, 2
 compile_opt IDL2

 ; Need to create unit conversion variable to convert from CIGALE units to Lightning units
 ; CIGALE in mJy
 conv = 4 * !dpi * (lghtng.lumin_dist * 1.d6 * !lightning_cgs.pc)^2 * !lightning_cgs.Jy/!lightning_cgs.Lsun
 ; Lsun/Hz/mJy   Lsun/Hz/Jy   Jy/mJy
 conv = conv * 1e-3

 ; Extract the CIGALE data from the structures and place them in arrays
 ; CIGALE wavelength is in units of nm
 cigale_wave = cigale_sed.wavelength * 1e-3
 cigale_wave_filters = ((lghtng.wave_filters)[sort(lghtng.wave_filters)])[where((lghtng.lnu_unc)[sort(lghtng.wave_filters)] ne 0)]

 ; Loop through and place input and model data into arrays, and convert units
 cigale_bands = strtrim(tag_names(cigale_obs), 2)
 cigale_bands = cigale_bands[(where(stregex(cigale_bands,'[^{_ERR}]$', /boolean) eq 1))[2:*]]
 cigale_phot = !null
 cigale_phot_err = !null
 cigale_phot_mod = !null
 for i=0,n_elements(cigale_bands)-1 do begin
   void = execute('cigale_phot = [cigale_phot, cigale_obs.'+strtrim(cigale_bands[i], 2)+']')
   void = execute('cigale_phot_err = [cigale_phot_err, cigale_obs.'+strtrim(cigale_bands[i], 2)+'_err]')
   void = execute('cigale_phot_mod = [cigale_phot_mod, cigale_results.best_'+strtrim(cigale_bands[i], 2)+']')
 endfor
 cigale_phot *= conv
 cigale_phot_err *= conv
 cigale_phot_mod *= conv

 ; Prep the Lightning X-ray data
 bestfit = where(lghtng.lnprob eq max(lghtng.lnprob))
 xray_bandpass = lghtng.xray_bandpass
 xray_obs = (lghtng.lnu_xray_obs)[*, bestfit]
 xray_unc = (lghtng.lnu_xray_unc)[*, bestfit]
 nxray = n_elements(xray_bandpass) / 2.d
 xray_obs[where(xray_obs le 0, /null)] = !values.d_nan
 xray_nu = xray_bandpass / (!lightning_cgs.hplanck / !lightning_cgs.kev)
 xray_wave = 1.d4*!lightning_cgs.clight/xray_nu
 xray_mean_wave = mean(xray_wave, dim=1)
 xray_mean_nu = 1.d4*!lightning_cgs.clight/xray_mean_wave


 ; Prepare the high resolution models. Select only the top 64% of best-fit models
 Lnu_mod_hires_lghtng = (lghtng.lnu_mod_hires)[*, 0:(0.64*n_elements((lghtng.lnu_mod_hires)[0, *]))]
 Lnu_mod_hires_star_lghtng = (lghtng.LNU_STARMOD_HIRES)[*, 0:(0.64*n_elements((lghtng.lnu_mod_hires)[0, *]))]
 Lnu_mod_hires_agn_lghtng = (lghtng.LNU_AGNMOD_HIRES)[*, 0:(0.64*n_elements((lghtng.lnu_mod_hires)[0, *]))]
 Lnu_mod_hires_dust_lghtng = (lghtng.LNU_DUSTMOD_HIRES)[*, 0:(0.64*n_elements((lghtng.lnu_mod_hires)[0, *]))]

 xray = lghtng.lnu_xraymod_agn_unabs_hires + lghtng.lnu_xraymod_star_unabs_hires
 Lnu_xmod_hires_lghtng = (xray)[*, 0:(0.64*n_elements((xray)[0, *]))]

 ; CIGALE components in units of W nm^-1
 cigale_nu = 1.d4 * !lightning_cgs.clight/cigale_wave

 ; Limit CIGALE high res model to above Lyman Break
 lyman_break = where(cigale_sed.wavelength * 1e-3 lt (0.0912 * (1 + lghtng.redshift)) and cigale_sed.wavelength * 1e-3 gt 1e-2)
 Lnu_mod_hires_cigale = cigale_nu * cigale_sed.fnu * conv
 Lnu_mod_hires_cigale[lyman_break] = 0


 ; Create SED plot window. We need a tall window to accommodate the SED and three residuals
 fig = window(dim=[600, 700])

 ; Plot the high resolution models
 nu_highres = 1.d4 * !lightning_cgs.clight/lghtng.wave_hires
 nu_highres_fill = rebin(reform(nu_highres, 1, n_elements(nu_highres)), 2, n_elements(nu_highres))

 xray_nu_highres = 1.d4 * !lightning_cgs.clight/lghtng.wave_xraymod_hires
 xray_nu_highres_fill = rebin(reform(xray_nu_highres, 1, n_elements(xray_nu_highres)), 2, n_elements(xray_nu_highres))


 ; Plot this uncertainty range of the Lightning models
 hires_fig = objarr(5)
 hires_fig_temp = fillplot(lghtng.wave_hires, nu_highres_fill*minmax(Lnu_mod_hires_star_lghtng, dim=2), /current, $
                           fill_color='red', linestyle='', transparency=70)
 hires_fig[4] = plot(lghtng.wave_hires, nu_highres*Lnu_mod_hires_star_lghtng[*, 0], /overplot, name='Stellar', $
                     color='red', thick=2)
 hires_fig_temp = fillplot(lghtng.wave_hires, nu_highres_fill*minmax(Lnu_mod_hires_dust_lghtng, dim=2), /over, $
                           fill_color='green', linestyle='', transparency=70)
 hires_fig[2] = plot(lghtng.wave_hires, nu_highres*Lnu_mod_hires_dust_lghtng[*, 0], /overplot, name='Dust', $
                     color='green', thick=2)
 hires_fig_temp = fillplot(lghtng.wave_hires, nu_highres_fill*minmax(Lnu_mod_hires_agn_lghtng, dim=2), /over, $
                           fill_color='orange', linestyle='', transparency=70)
 hires_fig[3] = plot(lghtng.wave_hires, nu_highres*Lnu_mod_hires_agn_lghtng[*, 0], /overplot, name='AGN', $
                     color='orange', thick=2)
 hires_fig_temp = fillplot(lghtng.wave_hires, nu_highres_fill*minmax(Lnu_mod_hires_lghtng, dim=2), /over, $
                           fill_color='gray', linestyle='', transparency=70)
 hires_fig_temp = fillplot(lghtng.wave_xraymod_hires, xray_nu_highres_fill*minmax(Lnu_xmod_hires_lghtng, dim=2), /over, $
                           fill_color='gray', linestyle='', transparency=70)
 hires_fig_temp = plot(lghtng.wave_xraymod_hires, xray_nu_highres*Lnu_xmod_hires_lghtng[*, 0], /overplot, $
                       color='gray', thick=2, name='Lightning')
 hires_fig[1] = plot(lghtng.wave_hires, nu_highres*Lnu_mod_hires_lghtng[*, 0], /overplot, name='Lightning', $
                     color='gray', thick=2)
 hires_fig[0] = plot(cigale_wave, cigale_nu * cigale_sed.fnu * conv, /overplot, name='CIGALE', $
                     color='magenta', thick=2)

 ; Plot the photometric data points
 nu = 1.d4 * !lightning_cgs.clight/lghtng.wave_filters
 fig = errorplot(lghtng.wave_filters, nu*lghtng.lnu_obs, nu*lghtng.lnu_unc, /overplot, /xlog, /ylog, $
                 linestyle='', name='Data', xrange=[3e-5, 1000], yrange=[4e10, 6e12], axis_style=1, $
                 symbol='o', /sym_filled, errorbar_thick=2, sym_size=1, xshowtext=0, $
                 ytitle='$\nu L_\nu$ [$L_\odot$]', position=[0.13, 0.38, 0.97, 0.93], ytickfont_size=13)

 fig_temp = errorplot(xray_mean_wave, xray_mean_nu*cigale_phot[0:14], xray_mean_nu*cigale_phot_err[0:14], /overplot, $
                      linestyle='', name='Data', symbol='o', /sym_filled, errorbar_thick=2, sym_size=1)
 fig_temp = plot(xray_mean_wave, xray_mean_nu*cigale_phot[0:14], /overplot, linestyle='', symbol='o', sym_size=1, sym_color='magenta', /sym_filled)
 fig_temp = plot(xray_mean_wave, xray_mean_nu*cigale_phot[0:14], /overplot, linestyle='', symbol='o', sym_size=1, sym_thick=2)

 fig_temp = errorplot(xray_mean_wave, xray_mean_nu*xray_obs, xray_mean_nu*xray_unc, /overplot, $
                      linestyle='', name='Data', symbol='o', /sym_filled, errorbar_thick=2, sym_size=1)
 fig_temp = plot(xray_mean_wave, xray_mean_nu*xray_obs, /overplot, linestyle='', symbol='o', sym_size=1, sym_color='gray', /sym_filled)
 fig_temp = plot(xray_mean_wave, xray_mean_nu*xray_obs, /overplot, linestyle='', symbol='o', sym_size=1, sym_thick=2)


 ; Add upper axis in energy
 yaxis = axis('Y', location='right', showtext=0, target=fig)
 xray_xrange=(1.d4 * !lightning_cgs.clight * (!lightning_cgs.hplanck / !lightning_cgs.keV)) / fig.xrange
 xpos = fig.position
 fig_temp = plot(xray_mean_wave, xray_mean_nu * xray_obs, /xlog, /ylog, xrange=xray_xrange, $
                 ytransparency=100, position=[xpos[0],xpos[3],xpos[2],xpos[3]+(xpos[3]-xpos[1])], $
                 xtickdir=1, xtextpos=1, /current, axis_style=1, yshowtext=0, xtickunits='scientific', $
                 /nodata, xtitle='Observed-Frame Energy [keV]', xtickfont_size=13)

 ; Include a legend
 lgnd = legend(target=[hires_fig, fig], transparency=30, position=[0.18, 0.65], $
               vertical_alignment=0, horizontal_alignment=0, font_size=14)

 ; Create residuals for each sed code. need to match model units to lightning units since using lightning observation
 lnu_unc_total = sqrt(lghtng.lnu_unc^2 + (lghtng.model_unc * lghtng.lnu_mod[*, bestfit])^2)
 lnu_unc_xray_total = sqrt(xray_unc^2 + (lghtng.model_unc * lghtng.lnu_xraymod[*, bestfit])^2)

 residuals_lghtng = (lghtng.lnu_obs - lghtng.lnu_mod[*, bestfit])/lnu_unc_total
 residuals_cigale = (cigale_phot - cigale_phot_mod)/cigale_phot_err

 residuals_xray_lghtng = (xray_obs - lghtng.lnu_xraymod[*, bestfit])/lnu_unc_xray_total
 residuals_xray_cigale = residuals_cigale[0:14]
 residuals_cigale = residuals_cigale[15:*]


 ; Determine range of residuals for y-axis range. Will be same for all codes for easier comparison
 residuals = [residuals_lghtng, residuals_xray_lghtng, residuals_cigale, residuals_xray_cigale]
 non_nan = where(finite(residuals), /null)
 yrange_resid = [-3.99, 3.99]

 ; Plot the residuals
 resid_fig = plot(fig.xrange, [0, 0], /current, ytitle='Residual [$\sigma$]', color='gray', ytickfont_size=13, $
                  /xlog, position=[0.13, 0.225, 0.97, 0.36], xrange=fig.xrange, yrange=yrange_resid, thick=2, $
                  xshowtext=0, ytickinterval=2)
 resid_fig = plot(fig.xrange, [1, 1], /over, color='light grey', thick=2, linestyle='--')
 resid_fig = plot(fig.xrange, [-1,-1], /over, color='light grey', thick=2, linestyle='--')
 fig = errorplot(lghtng.wave_filters, residuals_lghtng, replicate(1.d, n_elements(residuals_lghtng)), $
                 /overplot, linestyle='', symbol='o', /sym_filled, errorbar_thick=2, sym_size=1.2)
 fig = errorplot(xray_mean_wave, residuals_xray_lghtng, replicate(1.d, n_elements(residuals_xray_lghtng)), $
                 /overplot, linestyle='', symbol='o', /sym_filled, errorbar_thick=2, sym_size=1.2, sym_color='gray')
 fig = plot(xray_mean_wave, residuals_xray_lghtng, /overplot, linestyle='', symbol='o', sym_size=1.2, sym_color='gray', /sym_filled)
 fig = plot(xray_mean_wave, residuals_xray_lghtng, /overplot, linestyle='', symbol='o', sym_size=1.2, sym_thick=2)
 fig_text = text(0.15, 0.35, 'Lightning', vertical_alignment=1, font_size=14, color='gray')
 fig_text = text(0.40, 0.35, '$\chi^2_{best}$ = '+strtrim(string((lghtng.chi2)[bestfit], f='(F4.1)'), 2), $
                 vertical_alignment=1, font_size=14, color='gray')

 resid_fig = plot(fig.xrange, [0, 0], /current, ytitle='Residual [$\sigma$]', color='magenta', ytickfont_size=13, $
                  /xlog, position=[0.13, 0.07, 0.97, 0.205], xrange=fig.xrange, yrange=yrange_resid, thick=2, $
                  xtitle='Observed-Frame Wavelength [$\mu$m]', xtickfont_size=13, ytickinterval=2, xtickunits='scientific')
 resid_fig = plot(fig.xrange, [1, 1], /over, color='light grey', thick=2, linestyle='--')
 resid_fig = plot(fig.xrange, [-1,-1], /over, color='light grey', thick=2, linestyle='--')
 fig = errorplot(cigale_wave_filters, residuals_cigale, replicate(1.d, n_elements(residuals_cigale)), $
                 /overplot, linestyle='', symbol='o', /sym_filled, errorbar_thick=2, sym_size=1.2)
 fig = errorplot(xray_mean_wave, residuals_xray_cigale, replicate(1.d, n_elements(residuals_xray_cigale)), $
                 /overplot, linestyle='', symbol='o', /sym_filled, errorbar_thick=2, sym_size=1.2)
 fig = plot(xray_mean_wave, residuals_xray_cigale, /overplot, linestyle='', symbol='o', sym_size=1.2, sym_color='magenta', /sym_filled)
 fig = plot(xray_mean_wave, residuals_xray_cigale, /overplot, linestyle='', symbol='o', sym_size=1.2, sym_thick=2)
 fig_text = text(0.15, 0.195, 'CIGALE', vertical_alignment=1, font_size=14, color='magenta')
 fig_text = text(0.40, 0.195, '$\chi^2_{best}$ = '+strtrim(string(cigale_results.BEST_CHI_SQUARE, f='(F4.1)'), 2), $
                 vertical_alignment=1, font_size=14, color='magenta')

 return, fig

end