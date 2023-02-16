; Run the fit
restore, !lightning_dir + 'lightning.sav'

cd, !lightning_dir + 'examples/basic_NGC_337/'
lightning, 'ngc337_dale17_photometry.fits'

; Load the postprocessed results
res = mrdfits('lightning_output/ngc337_lightning_output.fits.gz',1)

; First we check on convergence and goodness of fit
print, '//Convergence for NGC 337//'
print, 'Mean acceptance fraction: ', strtrim(mean(res[0].acceptance_frac), 2)
print, 'Convergence flag: ', strtrim(res[0].convergence_flag, 2)
print, 'Short chain flag: ', strtrim(res[0].short_chain_flag, 2)
print, 'Number of "stranded" walkers: ', strtrim(total(res[0].stranded_flag), 2)
print, 'PPC p-value: ', strtrim(res[0].pvalue, 2)
print, ''

; Now we'll plot the best-fit model to make sure we're
; reaching a good solution
p = NGC337_sed_plot(res)

; p.save, 'images/NGC337_SED.png', dpi=400
; p.close

; Now, we investigate the stellar mass, SFR, and detailed SFH.
mstar_percentiles = percentile(res.mstar, [0.16, 0.50, 0.84])
SFR = 0.1 * res.psi[0,*] + 0.9 * res.psi[1,*]
SFR_percentiles = percentile(SFR, [0.16, 0.50, 0.84])

print, '//Basic properties for NGC 337//'
upstr = string(alog10(mstar_percentiles[2]) - alog10(mstar_percentiles[1]), format='(" (+",F5.2,")")')
downstr = string(alog10(mstar_percentiles[1]) - alog10(mstar_percentiles[0]), format='(" (-",F5.2,")")')
midstr = string(alog10(mstar_percentiles[1]), format='(F5.2)')
print, 'log(Mstar / Msun) = ', midstr, upstr, downstr

upstr = string(SFR_percentiles[2] - SFR_percentiles[1], format='(" (+",F5.2,")")')
downstr = string(SFR_percentiles[1] - SFR_percentiles[0], format='(" (-",F5.2,")")')
midstr = string(SFR_percentiles[1], format='(F5.2)')
print, 'SFR / [Msun yr-1] = ', midstr, upstr, downstr

p = NGC337_mstar_sfr_sfh_plots(res)

; p.save, 'images/NGC337_mstar_sfr_sfh.png', dpi=400
; p.close

; Finally, we can make a full corner plot for all the free parameters involved in the fit.
; Note that even for this 'simple' example our parameter space is 10-dimensional.
; Build an array of the posterior samples of our parameters
parameter_arr = transpose([[reform(res.psi[0,*])], $
                           [reform(res.psi[1,*])], $
                           [reform(res.psi[2,*])], $
                           [reform(res.psi[3,*])], $
                           [reform(res.psi[4,*])], $
                           [res.tauv_diff], $
                           [res.delta], $
                           [res.umin], $
                           [res.gamma], $
                           [res.qpah]])

; Labels for each parameter
labels = ['$\psi_1$', '$\psi_2$', '$\psi_3$', '$\psi_4$', '$\psi_5$',$
          '$\tau_V$', '$\delta$',$
          '$U_{rm min}$', '$\gamma$', '$q_{\rm PAH}$']

; Use the corner_plot function from the lightning-visualization repository
; to make our 10D corner plot
p = window(dimension=[900, 900])
p = corner_plot(parameter_arr, labels, $
                /normalize, $
                /show_median, $
                contour_levels=[0.68, 0.95], $
                contour_smooth=3, $
                tickinterval=[1,1,1,0.75,0.25,0.1,0.5,1.0,0.01,0.005], $
                distribution_thick=2, $
                contour_thick=2, $
                font_size=8, $
                position=[0.1, 0.1, 0.9, 0.9], $
                /current)

; p.save, 'images/NGC337_corner.png', dpi=400
; p.close
