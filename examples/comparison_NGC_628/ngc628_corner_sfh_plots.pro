function ngc628_corner_sfh_plots, lghtng, prospect, bags
;+
; Name
; ----
;   NGC628_CORNER_SFH_PLOTS
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
;       fig = ngc628_corner_sfh_plots(lghtng, prospect, bags)
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
 ;On_error, 2
 ;compile_opt IDL2


; Derive the total surviving stellar mass for each code
 tmass_lghtng   = lghtng.mstar
 tmass_prospect = total(prospect.chain[0:4, *], 1) * prospect.survive_mass_frac
 tmass_bags     = 10.d^bags.stellar_mass

; Derive the SFR of last 100 Myr for each code
 sfr_lghtng   = reform(0.1 * lghtng.psi[0, *] + 0.9 * lghtng.psi[1, *])
 sfr_prospect = reform(0.1 * prospect.sfh[0, *] + 0.9 * prospect.sfh[1, *])
 sfr_bags     = reform(0.1 * bags.sfh[0, *] + 0.9 * bags.sfh[1, *])

; Derive the sSFR of last 100 Myr for each code
 ssfr_lghtng   = sfr_lghtng/tmass_lghtng
 ssfr_prospect = sfr_prospect/tmass_prospect
 ssfr_bags     = sfr_bags/tmass_bags

; Derive the Av of each code
; tau_lam  = 0.4 * ln(10) * A_lam
; A_lam = 2.5/alog(10) * tau_lam
 av_lghtng   = 2.5d/alog(10.d) * lghtng.tauv
 av_prospect = 2.5d/alog(10.d) * reform(prospect.chain[where(prospect.THETA_LABELS eq 'dust2'), *])
 av_bags     = bags.dust_av

; Derive the LTIR (bolometric luminosity from 8um to 1000um) of Lightning
;  We did this already in the Prospector and BAGPIPES codes
 ltir_lghtng   = dblarr(2000)
 ltir_lghtng_idc   = where(lghtng.wave_hires ge 8 and lghtng.wave_hires le 1000, /null)
 ltir_prospect_idc = where(prospect.wave_hires/1e4 ge 8 and prospect.wave_hires/1e4 le 1000, /null)
 ltir_bags_idc     = where(bags.spec_wave/1e4 ge 8 and bags.spec_wave/1e4 le 1000, /null)
; Prospector and Bagpipes wave units are in Angstroms (i.e., 1e8 A = 1 cm)
 for i=0, 1999 do ltir_lghtng[i]   = -1.d * trapez(lghtng.lnu_mod_hires[ltir_lghtng_idc, i], $
                                    1.d4 * !lightning_cgs.clight / lghtng.wave_hires[ltir_lghtng_idc])
 ltir_prospect = prospect.ltir
 ltir_bags     = bags.ltir


; Compile the parameters into a single chain for plotting using homemade IDL corner plot function:
;    https://github.com/kjdoore/Useful_IDL_codes/blob/main/corner_plot.pro
 distribution_lghtng = transpose([[alog10(tmass_lghtng)], [sfr_lghtng], [alog10(ssfr_lghtng)], [av_lghtng], [alog10(ltir_lghtng)]])
 distribution_prospect = transpose([[alog10(tmass_prospect)], [sfr_prospect], [alog10(ssfr_prospect)], [av_prospect], [alog10(ltir_prospect)]])
 distribution_bags = transpose([[alog10(tmass_bags)], [sfr_bags], [alog10(ssfr_bags)], [av_bags], [alog10(ltir_bags)]])
 distribution = [[[distribution_bags]], [[distribution_prospect]], [[distribution_lghtng]]]

; Create the labels for each distribution
 distribution_labels = ['$log_{10}(M_\star$ [$M_\odot$])', 'SFR [$M_\odot$ $yr^{-1}$]', $
                        '$log_{10}$(sSFR [$yr^{-1}$])', '$A_V [mag]$', '$L_{TIR}$ [$L_\odot$]']

; Creat the plot window
 fig = window(dimension=[860, 800])

; Plot the corner plot using homemade function
 fig = corner_plot(distribution, $
                   distribution_labels, $
                   /normalize, $
                   /show_median, $
                   contour_levels=[0.68, 0.95], $
                   contour_smooth=5, $
                   distribution_color=['orange', 'green', 'blue'], $
                   distribution_thick=2, $
                   distribution_range=[[9.5, 10.25], [0, 2.5], [-11.5, -9.25], [0.2, 0.499], [9.75, 9.92]], $
                   contour_thick=2, $
                   tickinterval=[0.5, 1, 1, 0.1, 0.1], $
                   font_size=10, $
                   position=[0.08, 0.08, 0.92, 0.98], $
                   /current)

; Get age bin bounds and format for plotting. Also, truncate youngest bin, since it is at 0, and we ant log axes.
 Nsteps = n_elements(lghtng.steps_bounds) - 1
 age_bounds = lghtng.steps_bounds
 bounds = ([age_bounds, age_bounds])[sort([age_bounds,age_bounds])]
 bounds = bounds[1:-2]
 bounds[0] = (bounds[1])/10.d0

; Create the offset locations for the vertical bars indicating the error range on the SFH
 bounds_err = dblarr(3, 2, 5)
 for i=0, 4 do bounds_err[*, *, i] = rebin(10^(((alog10(bounds[(2*i)+1]) - alog10(bounds[(2*i)])) * $
                                                [0.25, 0.5, 0.75]) + alog10(bounds[(2*i)])), 3, 2)

; Create and format the SFH array and SFH uncertainty range array for plotting
 SFH_percent_lghtnng = dblarr(Nsteps, 3)
 SFH_percent_prospect = dblarr(Nsteps, 3)
 SFH_percent_bags = dblarr(Nsteps, 3)
 for i=0, Nsteps-1 do begin
   SFH_percent_lghtnng[i, *] = percentile(lghtng.psi[i, *], [0.5, 0.16, 0.84])
   SFH_percent_prospect[i, *] = percentile(prospect.sfh[i, *], [0.5, 0.16, 0.84])
   SFH_percent_bags[i, *] = percentile(bags.sfh[i, *], [0.5, 0.16, 0.84])
 endfor

 SFH_lghtng = dblarr(Nsteps * 2, 3)
 SFH_prospect = dblarr(Nsteps * 2, 3)
 SFH_bags = dblarr(Nsteps * 2, 3)
 for i=0, Nsteps-1 do begin
   SFH_lghtng[(2*i):(2*i)+1, *] = rebin(SFH_percent_lghtnng[i, *], 2, 3)
   SFH_prospect[(2*i):(2*i)+1, *] = rebin(SFH_percent_prospect[i, *], 2, 3)
   SFH_bags[(2*i):(2*i)+1, *] = rebin(SFH_percent_bags[i, *], 2, 3)
 endfor
 

; Plot the SFHs
; Plot the median SFH
 sfh_plt = objarr(3)
 sfh_plt[2] = plot(bounds, $
                   SFH_bags[*, 0], $
                   /current, $
                   color='orange', $
                   name='Bagpipes', $
                   thick=2, $
                   position=[0.615, 0.67, 0.98, 0.98], $
                   xrange=[1e6, max(bounds)], $
                   yrange=[7e-2, 10], $
                   xtitle='Stellar Age, t [yr]', $
                   ytitle='SFR(t) [$M_\odot$ $yr^{-1}$]', $
                   xtickfont_size=10, $
                   ytickfont_size=10, $
                   /xlog, $
                   /ylog, $
                   ytickunits='scientific')
 sfh_plt[1] = plot(bounds, $
                   SFH_prospect[*, 0], $
                   /overplot, $
                   color='green', $
                   name='Prospector', $
                   thick=2)
 sfh_plt[0] = plot(bounds, $
                   SFH_lghtng[*, 0], $
                   /overplot, $
                   color='blue', $
                   name='Lightning', $
                   thick=2)

; Plot the uncertainty range
 for i=0,4 do begin
   fig = plot(bounds_err[2, *, i], $
              SFH_bags[2*i, 1:2], $
              thick=2, $
              color='orange', $
              /overplot)
   fig = plot(bounds_err[0, *, i], $
              SFH_prospect[2*i, 1:2], $
              thick=2, $
              color='green', $
              /overplot)
   fig = plot(bounds_err[1, *, i], $
              SFH_lghtng[2*i, 1:2], $
              thick=2, $
              color='blue', $
              /overplot)
 endfor

; Include a legend
 lgnd = legend(target=sfh_plt, $
               transparency=30, $
               position=[0.28, 0.97], $
               vertical_alignment=1, $
               horizontal_alignment=0, $
               font_size=14)

 return, fig
 
end