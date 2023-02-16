; Block 1 - Run Lightning
restore, !lightning_dir + 'lightning.sav'

cd, !lightning_dir + 'examples/NGC_628/'
lightning, 'ngc628_dale17_photometry.fits'



; Block 2 - Read in the fit data
; Read in each post-processed files
lghtng   = mrdfits('lightning_output/ngc_628.fits.gz', 1)
prospect = mrdfits('prospector/ngc_628_prospector.fits.gz', 1)
bags     = mrdfits('bagpipes/ngc_628_bagpipes.fits.gz', 1)



; Block 3 - Initial convergence check
; Let's check the Lightning fit for convergence by all our metrics.
print, '//Convergence Checks//'
print, '//Convergence for Lightning fit//'
print, 'Mean acceptance fraction: ', strtrim(mean(lghtng.acceptance_frac), 2)
print, 'Convergence flag: ', strtrim(lghtng.convergence_flag, 2)
; Include the short chain flag since we had lightning automatically choose the burn-in and thinning
print, 'Short chain flag: ', strtrim(lghtng.short_chain_flag, 2)

; Let's check the Prospector fit for convergence by our autocorrelation and acceptance fraction metrics.
print, ''
print, '//Convergence for Prospector fit//'
print, 'Mean acceptance fraction: ', strtrim(mean(prospect.acceptance), 2)
print, 'Maximum autocorrelation time: ', strtrim(mean(prospect.autocorr_time), 2)



;Block 4 - Separate convergence checks
print, '//Detailed Lightning Convergence Checks//'
print, 'Number of parameters with high autocorrelation times: ', strtrim(total(lghtng.autocorr_flag), 2)
print, 'Number of walkers with low acceptance fractions: ', strtrim(total(lghtng.acceptance_flag), 2)



; Block 5 - Lightning acceptance frac check
print, '//Stranded Walkers Lightning fit//'
print, 'Number of walkers with low acceptance and not flagged as stranded: ', $
       strtrim(total((lghtng.acceptance_flag-lghtng.stranded_flag) > 0), 2)


; Block 6 - Make the SED plot
; Load the cgs constants, as the plotting script needs them
lightning_constants

; Plotting is simple with our figure script
fig = ngc628_sed_plots(lghtng, prospect, bags)



; Block 7 - Make the SED plot
; Plotting the distributions of the physical parameters
fig = ngc628_corner_sfh_plots(lghtng, prospect, bags)