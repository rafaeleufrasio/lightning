function ngc628_sed_plots, lghtng, prospect, bags
;+
; Name
; ----
;   NGC628_SED_PLOTS
;
; Purpose
; -------
;   Generates SED plots of NGC 628. Input tables
;   must be the post-processed outputs from Lightning, Prospector, and BAGPIPES.
;
; Calling Sequence
; ----------------
;   ::
;
;       fig = ngc628_sed_plots(lghtng, prospect, bags)
;
; Input
; -----
;   ``lghtng`` : structure
;       The Lightning post-processed output.
;   ``prospect`` : structure
;       The Prospector post-processed output.
;   ``bags`` : structure
;       The BAGPIPES post-processed output.
;
; Output
; ------
;   ``fig`` : object
;       The plot object containing the figure.
;
; Modification History
; --------------------
;   - 2023/02/14: Created (Keith Doore)
;-
 On_error, 2
 compile_opt IDL2


; All units are to be converted to Lightning output units
; Prospector units in maggies
; Bagpipes units in erg/s/cm2/Angstrom
; Both have wavelengths in Angstroms

; conversion from Fnu in Jy to Lnu in Lsun/Hz
 conv = 4 * !dpi * (lghtng.lumin_dist * 1.d6 * !lightning_cgs.pc)^2 * !lightning_cgs.Jy/!lightning_cgs.Lsun

; Prospector luminosity conversion    1 maggie = Jy/3631
; (Lsun/Hz)/(maggies)   Jy/(maggie)   Lsun/Hz/Jy
 prospect_conv   =         3631.d   *    conv

; Bagpipes luminosity conversions
; Fnu = lambda^2/c * Flambda
; (Lsun/Hz)/(erg/s/cm2/Angstrom) Jy/(erg/s/cm2/Hz) Lsun/Hz/Jy     um      Angstrom/um        cm/s         Angstrom/cm 
 bags_phot_conv   =                  1d23          * conv * (lghtng.wave_filters * 1d4)^2/(!lightning_cgs.clight * 1d8)
 bags_spec_conv   =                  1d23          * conv *      (bags.spec_wave)^2      /(!lightning_cgs.clight * 1d8)


; Get sorting of models by bestfit for residual plots
 sorted_bestfit_idc_lghtng   = sort(-1.d*lghtng.lnprob, /l64)
 sorted_bestfit_idc_prospect = sort(-1.d*prospect.lnprob, /l64)
 sorted_bestfit_idc_bags     = sort(bags.chisq_phot, /l64)

 bestfit_lghtng   = sorted_bestfit_idc_lghtng[0]
 bestfit_prospect = sorted_bestfit_idc_prospect[0]
 bestfit_bags     = sorted_bestfit_idc_bags[0]

; Perform a PPC for the Prospector and BAGPIPES fits
 pvalue_prospect = ppc(2000, lghtng.lnu_obs, lghtng.lnu_unc, prospect.phot_mod * prospect_conv, prospect.lnprob)
 pvalue_bags     = ppc(2000, lghtng.lnu_obs, lghtng.lnu_unc, bags.photometry * rebin(bags_phot_conv, 30, 2000), $
                       -0.5 * bags.chisq_phot)


; Prepare the high resolution models. Do this by interpolating them all onto the same wavelength grid and units
; Limit wavelength range to that above the Lyman break. (Issues occur if we include these lower wavelengths due 
;    to low luminosities)
 wave_hires = lghtng.wave_hires[where(lghtng.wave_hires gt 0.0912, /null)]
 Lnu_mod_hires_lghtng = lghtng.lnu_mod_hires[where(lghtng.wave_hires gt 0.0912, /null)]
 Lnu_mod_hires_prospect = interpol(prospect.spec_mod * prospect_conv, prospect.wave_hires*1.d-4, wave_hires)
 Lnu_mod_hires_bags     = interpol(bags.spectrum_full * bags_spec_conv, bags.spec_wave*1.d-4, wave_hires)


; Create SED plot window. We need a tall window to accommodate the SED and three residuals
 sed_fig = window(dim=[600, 800])

; Plot the high resolution models
 nu_highres = 1.d4 * !lightning_cgs.clight/wave_hires
 nu_highres_fill = rebin(reform(nu_highres, 1, n_elements(nu_highres)), 2, n_elements(nu_highres))

; Plot the best fit high resolution models
 hires_fig = objarr(3)
 hires_fig[0] = plot(wave_hires, $
                     nu_highres*Lnu_mod_hires_lghtng, $
                     /overplot, $
                     name='Lightning', $
                     color='blue', $
                     thick=2)
 hires_fig[1] = plot(wave_hires, $
                     nu_highres*Lnu_mod_hires_prospect, $
                     /overplot, $
                     name='Prospector', $
                     color='green', $
                     thick=2, $
                     linestyle='--')
 hires_fig[2] = plot(wave_hires, $
                     nu_highres*Lnu_mod_hires_bags, $
                     /overplot, $
                     name='Bagpipes', $
                     color='orange', $
                     thick=2, $
                     linestyle='-:')

; Plot the photometric data points
 nu = 1.d4 * !lightning_cgs.clight/lghtng.wave_filters
 fig = errorplot(lghtng.wave_filters, $
                 nu*lghtng.lnu_obs, $
                 nu*lghtng.lnu_unc, $
                 /overplot, $
                 /xlog, $
                 /ylog, $
                 linestyle='', $
                 name='Data', $
                 xrange=[0.08, 1000], $
                 yrange=[1e7, 2e10], $
                 symbol='o', $
                 /sym_filled, $
                 errorbar_thick=2, $
                 sym_size=1.5, $
                 xshowtext=0, $
                 ytitle='$\nu L_\nu$ [$L_\odot$]', $
                 position=[0.12, 0.5, 0.97, 0.98], $
                 ytickfont_size=13)

; Include a legend
 lgnd = legend(target=[hires_fig, fig], $
               transparency=30, $
               position=[0.18, 0.53], $
               vertical_alignment=0, $
               horizontal_alignment=0, $
               font_size=14)

; Create residuals for each SED code. Need to match model units to lightning units since using lightning observation
 residuals_lghtng = (lghtng.lnu_obs - lghtng.lnu_mod[*, bestfit_lghtng])/lghtng.lnu_unc
 residuals_prospect = (lghtng.lnu_obs - prospect.phot_mod[*, bestfit_prospect] * prospect_conv)/lghtng.lnu_unc
 residuals_bags = (lghtng.lnu_obs - bags.photometry[*, bestfit_bags] * bags_phot_conv)/lghtng.lnu_unc

; Determine range of residuals for y-axis range. Will be same for all codes for easier comparison
 residuals = [residuals_lghtng, residuals_prospect, residuals_bags]
 non_nan = where(finite(residuals), /null)
 yrange_resid = [-3.99, 5.5]

; Plot the residuals
 ; Plot the Lightning residuals
 resid_fig = plot(fig.xrange, $
                  [0, 0], $
                  /current, $
                  ytitle='Residual [$\sigma$]', $
                  color='blue', $
                  ytickfont_size=13, $
                  /xlog, $
                  position=[0.12, 0.36, 0.97, 0.48], $
                  xrange=fig.xrange, $
                  yrange=yrange_resid, $
                  thick=2, $
                  xshowtext=0)
 resid_fig = plot(fig.xrange, $
                  [1, 1], $
                  /overplot, $
                  color='light grey', $
                  thick=2, $
                  linestyle='--')
 resid_fig = plot(fig.xrange, $
                  [-1,-1], $
                  /overplot, $
                  color='light grey', $
                  thick=2, $
                  linestyle='--')
 fig = errorplot(lghtng.wave_filters, $
                 residuals_lghtng, $
                 replicate(1.d, n_elements(residuals_lghtng)), $
                 /overplot, $
                 linestyle='', $
                 symbol='o', $
                 /sym_filled, $
                 errorbar_thick=2, $
                 sym_size=1.5)
 fig_text = text(0.945, 0.45, 'Lightning', $
                 alignment=1, $
                 vertical_alignment=0, $
                 font_size=14, $
                 color='blue')
 fig_text = text(0.14, 0.45, 'p-value = '+string(lghtng.pvalue, f='(F5.3)'), $
                 alignment=0, $
                 vertical_alignment=0, $
                 font_size=14, $
                 color='blue')

 ; Plot the Prospector residuals
 resid_fig = plot(fig.xrange, $
                  [0, 0], $
                  /current, $
                  ytitle='Residual [$\sigma$]', $
                  color='green', $
                  ytickfont_size=13, $
                  /xlog, $
                  position=[0.12, 0.22, 0.97, 0.34], $
                  xrange=fig.xrange, $
                  yrange=yrange_resid, $
                  thick=2, $
                  xshowtext=0)
 resid_fig = plot(fig.xrange, $
                  [1, 1], $
                  /overplot, $
                  color='light grey', $
                  thick=2, $
                  linestyle='--')
 resid_fig = plot(fig.xrange, $
                  [-1,-1], $
                  /overplot, $
                  color='light grey', $
                  thick=2, $
                  linestyle='--')
 fig = errorplot(lghtng.wave_filters, $
                 residuals_prospect, $
                 replicate(1.d, n_elements(residuals_prospect)), $
                 /overplot, $
                 linestyle='', $
                 symbol='o', $
                 /sym_filled, $
                 errorbar_thick=2, $
                 sym_size=1.5)
 fig_text = text(0.945, 0.31, 'Prospector', $
                 alignment=1, $
                 vertical_alignment=0, $
                 font_size=14, $
                 color='green')
 fig_text = text(0.14, 0.31, 'p-value = '+string(pvalue_prospect, f='(F5.3)'), $
                 alignment=0, $
                 vertical_alignment=0, $
                 font_size=14, $
                 color='green')

 ; Plot the BAGPIPES residuals
 resid_fig = plot(fig.xrange, $
                  [0, 0], $
                  /current, $
                  ytitle='Residual [$\sigma$]', $
                  color='orange', $
                  ytickfont_size=13, $
                  /xlog, $
                  position=[0.12, 0.08, 0.97, 0.20], $
                  xrange=fig.xrange, $
                  yrange=yrange_resid, $
                  thick=2, $
                  xtickfont_size=13, $
                  xtitle='Observed-Frame Wavelength [$\mu$m]', $
                  xtickunits='scientific')
 resid_fig = plot(fig.xrange, $
                  [1, 1], $
                  /overplot, $
                  color='light grey', $
                  thick=2, $
                  linestyle='--')
 resid_fig = plot(fig.xrange, $
                  [-1,-1], $
                  /overplot, $
                  color='light grey', $
                  thick=2, $
                  linestyle='--')
 fig = errorplot(lghtng.wave_filters, $
                 residuals_bags, $
                 replicate(1.d, n_elements(residuals_bags)), $
                 /overplot, $
                 linestyle='', $
                 symbol='o', $
                 /sym_filled, $
                 errorbar_thick=2, $
                 sym_size=1.5)
 fig_text = text(0.945, 0.17, 'Bagpipes', $
                 alignment=1, $
                 vertical_alignment=0, $
                 font_size=14, $
                 color='orange')
 fig_text = text(0.14, 0.17, 'p-value = '+string(pvalue_bags, f='(F5.3)'), $
                 alignment=0, $
                 vertical_alignment=0, $
                 font_size=14, $
                 color='orange')

 return, fig

end