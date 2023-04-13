; In this batch file we list the commands we use to analyze
; the fits to the AGN example. Commands can be copied and
; pasted into the IDL command line, or the batch file
; can be run as a whole with the @ syntax.

; Run Lightning
cd, !lightning_dir + 'examples/AGN_J033226/'
restore, !lightning_dir + 'lightning.sav'
lightning, 'noxray/J033226_lightning_input.fits'
lightning, 'xray/J033226_lightning_input.fits'

; Load the result tables
noxray_res = mrdfits('noxray/lightning_output/J033226_noxray_output.fits.gz', 1)
xray_res = mrdfits('xray/lightning_output/J033226_xray_output.fits.gz', 1)

; First we check for reasonable acceptance fractions and convergence in the MCMC chains
print, '//Convergence for the fit without the X-ray model//'
print, 'Mean acceptance fraction: ', strtrim(mean(noxray_res.ACCEPTANCE_FRAC), 2)
print, 'Convergence flag: ', strtrim(noxray_res.CONVERGENCE_FLAG, 2)
print, 'Short chain flag: ', strtrim(noxray_res.SHORT_CHAIN_FLAG, 2)
print, 'Number of "stranded" walkers: ', strtrim(total(noxray_res.STRANDED_FLAG), 2)

print,''

print, '//Convergence for the fit with the X-ray model//'
print, 'Mean acceptance fraction: ', strtrim(mean(xray_res.ACCEPTANCE_FRAC), 2)
print, 'Convergence flag: ', strtrim(xray_res.CONVERGENCE_FLAG, 2)
print, 'Short chain flag: ', strtrim(xray_res.SHORT_CHAIN_FLAG, 2)
print, 'Number of "stranded" walkers: ', strtrim(total(xray_res.STRANDED_FLAG), 2)
print, ''
print, 'Number of walkers with low acceptance fractions: ', strtrim(total(xray_res.ACCEPTANCE_FLAG), 2)
print, "Number of parameters with high autocorrelation times: ", strtrim(total(xray_res.AUTOCORR_FLAG), 2)

print, ''

print, '//Autocorrelation Time Checks//'
print, 'Parameter with high autocorrelation time: ', $
       strtrim((xray_res.parameter_names)[where(xray_res.autocorr_flag eq 1, /null)], 2)
print, 'Ratio of Ntrials to autocorr_time: ', $
       strtrim(4e4/((xray_res.autocorr_time)[where(xray_res.autocorr_flag eq 1, /null)]), 2)
print, ''


; Now we'll plot the SED and folded X-ray spectral fit for the fit with the X-ray
; model. Of course we could plot the other fit too, but it's qualitatively similar.
p = J033226_spectrum_plots(xray_res)

; p.save, 'images/figure_J033226_SED.eps' ; for the paper
; p.save, 'images/J033226_SED.png', dpi=400 ; for RTD
; p.close

; To compare the posteriors on equal footing, it'd be nice
; to have the integrate optical-IR luminosity of the AGN model
; for the fit with the X-ray model.
restore, 'xray/lightning_output/lightning_configure.sav'
; The ARF is not strictly necessary, but LIGHTNING_MODELS expects it
arf = mrdfits('xray/J033226_summed.arf', 1)
AGN_model_lum = calc_integrated_AGN_luminosity(xray_res, config, arf)

help, AGN_model_lum

; save, AGN_model_lum, 'AGN_model_luminosity_xray.sav'

; Now we'll make a cornerplot of all the AGN parameters
p = posterior_comparison(xray_res, noxray_res, AGN_model_lum)

; p.save, 'images/figure_J033226_corner_params.eps' ; for the paper
; p.save, 'images/J033226_corner_params.png', dpi=400 ; for RTD
; p.close
