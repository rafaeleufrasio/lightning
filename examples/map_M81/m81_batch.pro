; Block 1 - Create B-band mask
cd, !lightning_dir + 'examples/map_M81/'

; These are the original B-band ellipse values from HyperLeda.
ra = [09, 55, 33.15]   ; units in Hour, minute, and seconds.
dec = [69, 03, 55.2]   ; units in sexagesimal degrees.
logd25 = 2.34  ; logd25 is the decimal logarithm of the length the projected major axis (d25 in 0.1 arcmin).
logr25 = 0.28  ; logr25 is the axis ratio of the isophote (decimal logarithm of the ratio of the lengths of major to minor axes).
pa = 157.0	   ; Major axis position angle in degrees (North Eastwards).

; Let's convert these values to degrees.
ra = 15 * (ra[0] + ra[1]/60. + ra[2]/3600.) 
dec = (dec[0] + dec[1]/60. + dec[2]/3600.) 
a = ((10.d^logd25 * 0.1d) /2) / 60.d  ; divide by 2 for semi-major.
q = 10.d^(-1.d * logr25)   ; axis ratio is q = b/a, negative flips ratio from major/minor to minor/major.
b = a * q                  ; semi-minor axis.

; Read in the image containing the SED_ID of each pixel and its header.
m81_map  = mrdfits('m81_data.fits', 0, hd)
; Create the ellipse mask.
bband_mask = ellipse_region(ra, dec, a, b, pa, hd)



; Block 2 - B-band mask image
img_sed_id = image(m81_map, $
                   rgb_table=33, $
                   position=[0, 0, 0.8, 1])
cb = colorbar(target=img_sed_id, $
              orientation=1, $
              position=[0.83,0.05,0.86,0.95], $
              title='SED ID', $
              /border_on, $
              textpos=1, $
              tickdir=1, $
              font_size=12)
img_mask   = image(bband_mask, $
                   transparency=80, $
                   /overplot)



; Block 3 - Create subset files
; Read in the data table containing the SEDs.
m81_seds  = mrdfits('m81_data.fits', 1)

; We need the SED_IDs that correspond to the pixels within the mask.
; Note that we exclude any pixels that lack an SED, which is given by a NaN in m81_map.
bband_ids = m81_map[where(bband_mask eq 1 and finite(m81_map))]

; Since the SED_IDs in m81_seds are in ascending order based on value, we can just
;   index the table using the selected IDs.
fitting_data = m81_seds[bband_ids]

; Let's create individual files for the MCMC and MPFIT fittings.
mwrfits, fitting_data, 'MCMC/m81_input.fits', /create
mwrfits, fitting_data, 'MPFIT/m81_input.fits', /create



; Block 4 - Fit MCMC
;restore, !lightning_dir + 'lightning.sav'
;lightning, 'MCMC/m81_input.fits'



; Block 5 - Fit MPFIT
;lightning, 'MPFIT/m81_input.fits'



; Block 6 - Initial convergence check
; Read in the post-processed data table.
results = mrdfits('MCMC/lightning_output/m81_map_mcmc.fits.gz', 1)

; Let's check how many fits, if any, did not reach convergence by all our metrics.
print, '//Convergence Checks//'
print, 'Number with failed convergence flag: ', strtrim(total(results.convergence_flag), 2)



; Block 7 - Acceptance and autocorrelation check
; These totals count the number of SEDs with the corresponding flags set.
;   The < operator sets any value > 1 from the first total to 1.
;   This allows us to see the number of SED with an issue vs all flags.
print, '//Detailed Convergence Checks//'
print, 'Number with failed autocorrelation flag: ', strtrim(total(total(results.autocorr_flag, 1) < 1), 2)
print, 'Number with failed acceptance flag: ', strtrim(total(total(results.acceptance_flag, 1) < 1), 2)



; Block 8 - Stranded check
; This counts the number of SEDs with an acceptance flag set when the corresponding
;   stranded flag is not.
print, '//Stranded Walker Check//'
print, "Number of walkers with low acceptance fractions that /weren't/ flagged as stranged: ", $
       strtrim(total(total((results.acceptance_flag - results.stranded_flag) > 0, 1) < 1), 2)



; Block 9 - Acceptance plot
; Let's include the 1 standard deviation of the acceptance fraction, since
;   it is what is used to flag stranded walkers.
plt = errorplot(total(results.acceptance_flag, 1)/75, $
                median(results.acceptance_frac, dim=1), $
                stddev(results.acceptance_frac, dim=1), $
                linestyle='', $
                symbol='o', $
                xtitle='Fraction of walkers below 20% acceptance threshold', $
                ytitle='Acceptance Fraction')
plt = plot(plt.xrange, $
           [0.2, 0.2], $
           color='red', $
           linestyle='--', $
           thick=2, $
           /overplot)



; Block 10 - MPFIT convergence check
; Read in the MPFIT post-processed data table.
results_mpfit = mrdfits('MCMC/lightning_output/m81_map_mpfit.fits.gz', 1)

; Let's check how many fits, if any, did not reach convergence by all our metrics.
print, '//Convergence Checks//'
print, 'Number with failed convergence flag: ', strtrim(total(results_mpfit.convergence_flag), 2)



; Block 11 - MPFIT metric checks
; These totals count the number of SEDs with the corresponding flags set.
print, '//Detailed Convergence Checks//'
print, 'Number with failed status flag: ', strtrim(total(total(results_mpfit.status_flag, 1) < 1), 2)
print, 'Number that reach max iterations: ', strtrim(total(total(results_mpfit.iter_flag, 1) < 1), 2)
print, 'Number with < 50% of solvers near lowest chisqr: ', strtrim(total(results_mpfit.stuck_flag), 2)
print, 'Number with a parameter not within 1% of lowest chisqr solver: ', strtrim(total(total(results_mpfit.similar_flag, 1) < 1), 2)



; Block 12 - Iteration fraction check
print, '//Iteration Fraction Check//'
print, 'Maximum solvers to reach iteration limit: ', strtrim(max(total(results_mpfit.iter_flag, 1)), 2)



; Block 13 - Stuck fraction check
; Calculate the chi2 difference. chi2 values need to be recomputed from p-values.
;   This reconversion is very accurate as long as your chi2 is not a terrible fit (i.e., p-value = 0).
chi2 = dblarr(20, n_elements(results_mpfit))
for i=0, 19 do for j=0, n_elements(results_mpfit)-1 do $
    chi2[i, j] = CHISQR_CVF((results_mpfit[j].pvalue)[i], results_mpfit[j].dof)
chi2_diff = chi2 - rebin(reform(min(chi2, dim=1), 1, n_elements(results_mpfit)), 20, n_elements(results_mpfit))

; Determine how many solvers were below a given chi2 threshold difference from the minimum
reached_min = intarr(n_elements(results_mpfit))
for j=0, n_elements(results_mpfit)-1 do reached_min[j] = n_elements(where(chi2_diff[*, j] lt 4, /null))

; Only look at the fits with the stuck_flag set.
keep = where(results_mpfit.stuck_flag eq 1, nkeep, /null)

; Creat the distribution and plot
counts = histogram(reached_min[keep], nbins=max(reached_min)-min(reached_min[keep])+1, locations=binloc)
plt = plot(binloc, $
           counts, $
           xtitle='Number of solvers with $\chi^2 - 4 \leq \chi^2_{best}$', $
           ytitle='Counts', $
           thick=2, $
           /stair)



; Block 14 - Similar fraction check
; Calculate the number of non-stuck solvers with similar (1, 3 and 5% difference) solutions to the bestfit
;   (if the bestfit value is 0, just use difference)
Nsimilar = dblarr(3, n_elements(results_mpfit[0].parameter_names), n_elements(results_mpfit))
for j=0, n_elements(results_mpfit)-1 do begin $
    min_loc = (where(chi2_diff[*, j] eq 0, /null))[0] &$
    bestfit_params = (results_mpfit[j].parameter_values)[*, min_loc] &$
    zero_param = where(bestfit_params eq 0, comp=nonzero_param, ncomp=Nnonzero, /null) &$
    near_min = where(chi2_diff[*, j] lt 4, Nnear_min) &$
    param_per_diff = (results_mpfit[j].parameter_values)[*, near_min] - $
                     rebin(bestfit_params, n_elements(results_mpfit[0].parameter_names), Nnear_min) &$
    if Nnonzero gt 0 then $
      param_per_diff[nonzero_param, *] = param_per_diff[nonzero_param, *] / $
                                         rebin((results_mpfit[j].parameter_values)[nonzero_param, min_loc], Nnonzero, Nnear_min) &$
    Nsimilar[0, *, j] = fix(total(fix(param_per_diff lt 0.01d), 2)) &$
    Nsimilar[1, *, j] = fix(total(fix(param_per_diff lt 0.03d), 2)) &$
    Nsimilar[2, *, j] = fix(total(fix(param_per_diff lt 0.05d), 2))

; Determine the fraction of non-stuck solvers with solutions near the bestfit
similar_frac = Nsimilar / rebin(reform(reached_min, 1, 1, n_elements(results_mpfit)), 3, 13, n_elements(results_mpfit))

; Only look at the fits with the similar_flag set.
keep = where(total(results_mpfit.similar_flag, 1) gt 0, /null)

; Check if majority of non-stuck solvers reached the same solution
solver_frac = dblarr(3, 13)
for i=0, 2 do for j=0,12 do solver_frac[i, j] = n_elements(where(similar_frac[i, j, keep] gt 0.5,/null)) / float(n_elements(keep))

print, '// Similarity fractions per parameter //'
print, 'Parameter Name 1% diff         3% diff         5% diff'
forprint, results_mpfit[0].parameter_names, solver_frac[0,*], solver_frac[1,*], solver_frac[2,*]



; Block 15 - MCMC property maps
nsed = n_elements(results)
; Get the dimensions of the map so we can make blank ones to place our parameter values.
map_size = size(m81_map, /dim)

; We need to match each SED_ID from the result back to the map.
map_idc = lonarr(nsed)
for i=0, nsed-1 do map_idc[i] = where(m81_map eq float(strtrim(results[i].sed_id, 2)))

; Make a blank map to fill in with the results values.
sfr_map = replicate(!values.f_NaN, map_size)
; The SFR of the last 100 My is the weighted sum of the first two bins (0-10 and 10-100 Myr).
;   The bin index is the first dimension. The second dimension is the percentile values (16th, 50th, and 84th).
sfr_map[map_idc] = (results.psi_percentiles)[0, 1, *] * 0.1 + (results.psi_percentiles)[1, 1, *] * 0.9

mass_map = replicate(!values.f_NaN, map_size)
mass_map[map_idc] = median(results.mstar, dim=1)

tauv_diff_map = replicate(!values.f_NaN, map_size)
tauv_diff_map[map_idc] = median(results.tauv_diff, dim=1)

ltir_map = replicate(!values.f_NaN, map_size)
ltir_map[map_idc] = median(results.ltir, dim=1)

pvalue_map = replicate(!values.f_NaN, map_size)
pvalue_map[map_idc] = results.pvalue



; Block 16 - Crop MCMC maps
; Find the min and max values of where B-band mask eq 1 in each dimension.
mask_bounds_idc = minmax(array_indices(bband_mask, where(bband_mask eq 1)), dim=2)

; Crop each image to only include the B-band region.
sfr_map = sfr_map[mask_bounds_idc[0, 0]:mask_bounds_idc[1, 0], *]
sfr_map = sfr_map[*, mask_bounds_idc[0, 1]:mask_bounds_idc[1, 1]]

mass_map = mass_map[mask_bounds_idc[0, 0]:mask_bounds_idc[1, 0], *]
mass_map = mass_map[*, mask_bounds_idc[0, 1]:mask_bounds_idc[1, 1]]

tauv_diff_map = tauv_diff_map[mask_bounds_idc[0, 0]:mask_bounds_idc[1, 0], *]
tauv_diff_map = tauv_diff_map[*, mask_bounds_idc[0, 1]:mask_bounds_idc[1, 1]]

ltir_map = ltir_map[mask_bounds_idc[0, 0]:mask_bounds_idc[1, 0], *]
ltir_map = ltir_map[*, mask_bounds_idc[0, 1]:mask_bounds_idc[1, 1]]

pvalue_map = pvalue_map[mask_bounds_idc[0, 0]:mask_bounds_idc[1, 0], *]
pvalue_map = pvalue_map[*, mask_bounds_idc[0, 1]:mask_bounds_idc[1, 1]]



; Block 17 - Plot maps
; Plotting is simple with our figure script
fig = m81_maps(sfr_map, mass_map, sfr_map/mass_map, tauv_diff_map, ltir_map, pvalue_map)



; Block 18 - Ks-band mask
; Data gathered from Jarrett et al 2003.
ra = [09, 55, 33.1]     ; units in Hour, minute, and seconds.
dec = [69, 03, 54.9]    ; units in sexagesimal degrees.
r20 = 487.6   ; The length of the projected major axis at the isophotal level 20 mag/arcsec2 in the Ks-band (arcsec).
q = 0.51      ; The axis ratio of the isophote 20 mag/arcsec2 in the Ks-band (minor to major axes b/a).
pa = -31.0	  ; Major axis position angle in degrees (North Eastwards).

; Let's convert these values to degrees.
ra = 15 * (ra[0] + ra[1]/60. + ra[2]/3600.) 
dec = (dec[0] + dec[1]/60. + dec[2]/3600.) 
; Note that we use 1/2 the Ks-band 20 mag isophote
a = (r20/2.) / 3.6d3  ; semi-major axis in degrees
b = a * q             ; semi-minor axis in degrees

; Create the ellipse mask.
ksband_mask = ellipse_region(ra, dec, a, b, pa, hd)



; Block 19 - Ks-band mask image
img_sed_id = image(m81_map, $
                   rgb_table=33, $
                   position=[0, 0, 0.8, 1])
cb = colorbar(target=img_sed_id, $
              orientation=1, $
              position=[0.83,0.05,0.86,0.95], $
              title='SED ID', $
              /border_on, $
              textpos=1, $
              tickdir=1, $
              font_size=12)
img_mask   = image(ksband_mask, $
                   transparency=80, $
                   /overplot)



; Block 20 - Plot regions
; Plotting is simple with our figure script
fig = m81_regions(m81_map, hd, m81_seds, results, bband_mask, ksband_mask)



; Block 21 - MPFIT property maps
; We need to match each SED_ID from the result back to the map.
map_idc = lonarr(nsed)
for i=0, nsed-1 do map_idc[i] = where(m81_map eq float(strtrim(results_mpfit[i].sed_id, 2)))

; Make a blank map to fill in with the results values.
sfr_map_mpfit = replicate(!values.f_NaN, map_size)
; The SFR of the last 100 My is the weighted sum of the first two bins (0-10 and 10-100 Myr).
;   The bin index is the first dimension.
sfr_map_mpfit[map_idc] = (results_mpfit.psi)[0, *] * 0.1 + (results_mpfit.psi)[1, *] * 0.9

mass_map_mpfit = replicate(!values.f_NaN, map_size)
mass_map_mpfit[map_idc] = results_mpfit.mstar

tauv_diff_map_mpfit = replicate(!values.f_NaN, map_size)
tauv_diff_map_mpfit[map_idc] = results_mpfit.tauv_diff

ltir_map_mpfit = replicate(!values.f_NaN, map_size)
ltir_map_mpfit[map_idc] = results_mpfit.ltir

pvalue_map_mpfit = replicate(!values.f_NaN, map_size)
pvalue_map_mpfit[map_idc] = max(results_mpfit.pvalue, dim=1)



; Block 22 - Crop MPFIT maps
; Crop each image to only include the B-band region.
sfr_map_mpfit = sfr_map_mpfit[mask_bounds_idc[0, 0]:mask_bounds_idc[1, 0], *]
sfr_map_mpfit = sfr_map_mpfit[*, mask_bounds_idc[0, 1]:mask_bounds_idc[1, 1]]

mass_map_mpfit = mass_map_mpfit[mask_bounds_idc[0, 0]:mask_bounds_idc[1, 0], *]
mass_map_mpfit = mass_map_mpfit[*, mask_bounds_idc[0, 1]:mask_bounds_idc[1, 1]]

tauv_diff_map_mpfit = tauv_diff_map_mpfit[mask_bounds_idc[0, 0]:mask_bounds_idc[1, 0], *]
tauv_diff_map_mpfit = tauv_diff_map_mpfit[*, mask_bounds_idc[0, 1]:mask_bounds_idc[1, 1]]

ltir_map_mpfit = ltir_map_mpfit[mask_bounds_idc[0, 0]:mask_bounds_idc[1, 0], *]
ltir_map_mpfit = ltir_map_mpfit[*, mask_bounds_idc[0, 1]:mask_bounds_idc[1, 1]]

pvalue_map_mpfit = pvalue_map_mpfit[mask_bounds_idc[0, 0]:mask_bounds_idc[1, 0], *]
pvalue_map_mpfit = pvalue_map_mpfit[*, mask_bounds_idc[0, 1]:mask_bounds_idc[1, 1]]



; Block 23 - Scale the images to the MCMC images
keep_nan = where(~finite(sfr_map_mpfit))

; The > and < operators set values below/above that value to the specified value
pvalue_map_mpfit = 10^(((alog10(pvalue_map_mpfit) > min((alog10(pvalue_map))[where(finite(alog10(pvalue_map)))])) < $
                   max((alog10(pvalue_map))[where(finite(alog10(pvalue_map)))])))
sfr_map_mpfit  = ((sfr_map_mpfit > min(sfr_map)) < max(sfr_map))
mass_map_mpfit = ((mass_map_mpfit > min(mass_map)) < max(mass_map))
ssfr_map_mpfit = 10^((alog10(sfr_map_mpfit/mass_map_mpfit)   > (min(alog10(sfr_map/mass_map)))) < $
                  max(alog10(sfr_map/mass_map)))
tauv_diff_map_mpfit   = ((tauv_diff_map_mpfit > min(tauv_diff_map)) < max(tauv_diff_map))
ltir_map_mpfit   = ((ltir_map_mpfit > min(ltir_map)) < max(ltir_map))

; The > and < operators override NaNs, set the NaN pixels back to NaN
pvalue_map_mpfit[keep_nan] = !values.D_NaN
sfr_map_mpfit[keep_nan] = !values.D_NaN
mass_map_mpfit[keep_nan]   = !values.D_NaN
ssfr_map_mpfit[keep_nan]   = !values.D_NaN
tauv_diff_map_mpfit[keep_nan]   = !values.D_NaN
ltir_map_mpfit[keep_nan]   = !values.D_NaN



; Block 24 - Plot maps
; Plotting is simple with our figure script
fig = m81_maps(sfr_map_mpfit, mass_map_mpfit, ssfr_map_mpfit, tauv_diff_map_mpfit, ltir_map_mpfit, pvalue_map_mpfit)
