; Block 1 - Create inclination prior
cd, !lightning_dir + 'examples/NGC_4631/'

; Axis ratio gathered from HyperLeda database.
logr25 = 0.82  ; logr25 is the axis ratio of the isophote (decimal logarithm of the ratio of the lengths of major to minor axes).
elogr25 = 0.05 ; error on logr25

; Convert from logr25 to q
q = 10.d^(-1.d * logr25) ; axis ratio is q = b/a, negative flips ratio from major/minor to minor/major.
q_err = q * alog(10.d) * elogr25

; Perform the Monte Carlo (MC) simulation
inclination_pdf, q, q_err, 'gamma_eps_distribution_tables/', inc_bins, inc_pdf, units='cosi'

; Now add the results from the MC simulation into the prior structure.
;   Also, let's smooth the pdf some since it is randomly generated.
prior = {cosi_values: inc_bins, cosi_pdf: smooth(inc_pdf, 5)}

; Save the structure to the correctly named FITS file.
mwrfits, prior, 'tabulated_prior_NGC_4631.fits', /create



; Block 2 - Plot inclination prior
plt = plot(1 - prior.cosi_values, $
           prior.cosi_pdf, $
           xtitle='1 - cos(i)', $
           ytitle='Probability', $
           thick=2)



; Block 3 - Create subset files
; Copy the photometry into the fitting directories
file_copy, 'ngc4631_dale17_photometry.fits', 'Calz/', /allow_same
file_copy, 'ngc4631_dale17_photometry.fits', 'Doore/', /allow_same
file_copy, 'ngc4631_dale17+xrays_photometry.fits', 'Calz+Xrays/', /allow_same
file_copy, 'ngc4631_dale17+xrays_photometry.fits', 'Doore+Xrays/', /allow_same



; Block 4 - Copy Lbol_abs_table
; Copy the Lbol_abs table directory into the Doore21 fitting directories
file_mkdir, ['Doore+Xrays/', 'Doore/'] + 'lightning_output/'
file_copy, replicate('doore21_Lbol_abs_table/', 2), $
           ['Doore+Xrays/', 'Doore/'] + 'lightning_output/', /allow_same, /recursive



; Block 5 - Fit the SED
restore, !lightning_dir + 'lightning.sav'

lightning, 'Calz/ngc4631_dale17_photometry.fits'
lightning, 'Calz+Xrays/ngc4631_dale17+xrays_photometry.fits'
lightning, 'Doore/ngc4631_dale17_photometry.fits'
lightning, 'Doore+Xrays/ngc4631_dale17+xrays_photometry.fits'



; Block 5 - Read in the fit data
; Read in each post-processed files
calz    = mrdfits('Calz/lightning_output/calz_no_xrays.fits.gz', 1)
calz_x  = mrdfits('Calz+Xrays/lightning_output/calz_with_xrays.fits.gz', 1)
doore   = mrdfits('Doore/lightning_output/doore_no_xrays.fits.gz', 1)
doore_x = mrdfits('Doore+Xrays/lightning_output/doore_with_xrays.fits.gz', 1)



; Block 6 - Initial convergence check
; Let's check each fit, if any, did not reach convergence by all our metrics.
print, '//Convergence Checks//'
print, '//Convergence for the Calzetti (no X-rays) fit//'
print, 'Mean acceptance fraction: ', strtrim(mean(calz.acceptance_frac), 2)
print, 'Convergence flag: ', strtrim(calz.convergence_flag, 2)
; Include the short chain flag since we had lightning automatically choose the burn-in and thinning
print, 'Short chain flag: ', strtrim(calz.short_chain_flag, 2)

print, ''
print, '//Convergence for the Calzetti with X-rays fit//'
print, 'Mean acceptance fraction: ', strtrim(mean(calz_x.acceptance_frac), 2)
print, 'Convergence flag: ', strtrim(calz_x.convergence_flag, 2)
print, 'Short chain flag: ', strtrim(calz_x.short_chain_flag, 2)

print, ''
print, '//Convergence for the Doore+21 (no X-rays) fit//'
print, 'Mean acceptance fraction: ', strtrim(mean(doore.acceptance_frac), 2)
print, 'Convergence flag: ', strtrim(doore.convergence_flag, 2)
print, 'Short chain flag: ', strtrim(doore.short_chain_flag, 2)

print, ''
print, '//Convergence for the Doore+21 with X-rays  fit//'
print, 'Mean acceptance fraction: ', strtrim(mean(doore_x.acceptance_frac), 2)
print, 'Convergence flag: ', strtrim(doore_x.convergence_flag, 2)
print, 'Short chain flag: ', strtrim(doore_x.short_chain_flag, 2)



; Block 7 - Detailed convergence check
print, '//Detailed Convergence Checks//'
print, 'Number of parameters with high autocorrelation times: ', strtrim(total(doore.autocorr_flag), 2)
print, 'Number of walkers with low acceptance fractions: ', strtrim(total(doore.acceptance_flag), 2)



; Block 8 - Autocorr_time check
print, '//Autocorrelation Time Checks//'
print, 'Parameter with high autocorrelation time: ', $
       strtrim((doore.parameter_names)[where(doore.autocorr_flag eq 1, /null)], 2)
print, 'Ratio of Ntrials to autocorr_time: ', $
       strtrim(5e4/((doore.autocorr_time)[where(doore.autocorr_flag eq 1, /null)]), 2)



; Block 9 - Make the plot
; Load the cgs constants, as the plotting script needs them
lightning_constants

; Plotting is simple with our figure script
fig = ngc4631_plots(calz, calz_x, doore, doore_x)


