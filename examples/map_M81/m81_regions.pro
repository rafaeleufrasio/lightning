function m81_regions, m81_map, hd, m81_seds, results, bband_mask, ksband_mask
;+
; Name
; ----
;   M81_REGIONS
;
; Purpose
; -------
;   Generates a figure showing the inner and outer regions of M81. Input must
;   be formatted and generated as show in Lightning's online documentation.
;
; Calling Sequence
; ----------------
;   ::
;
;       fig = m81_regions(m81_map, hd, m81_seds, results, bband_mask, ksband_mask)
;
; Input
; -----
;   ``m81_map`` : float or double array(N, M)
;       The map containing the ``SED_ID`` of each pixel.
;   ``hd`` : string array(I)
;       The FITS header associated with ``m81_map``.
;   ``m81_seds`` : structure
;       The original input data structure.
;   ``results`` : structure
;       The post-processed MCMC results structure.
;   ``bband_mask`` : int, float, or double array(N, M)
;       The B-band mask map.
;   ``ksband_mask`` : int, float, or double array(N, M)
;       The Ks-band mask map.
;
; Output
; ------
;   ``fig`` : object
;       The plot object containing the figure.
;
; Modification History
; --------------------
;   - 2023/02/08: Created (Keith Doore)
;-
 On_error, 2
 compile_opt IDL2

; First let's map our original SDSS SED data back to the map so we can make a g, r, i composite
 map_size = size(m81_map, /dimension)
 map_idc = lonarr(n_elements(m81_seds))
 for i=0, n_elements(m81_seds)-1 do map_idc[i] = where(m81_map eq float(strtrim(m81_seds[i].sed_id, 2)))

; IDL does not like mapping 2 dimensions with one index while fixing the other.
;   So, we use a temporary variable to do this.
 sdss_map = replicate(!values.f_NaN, [map_size, 3])
 temp = replicate(!values.f_NaN, map_size)
 temp[map_idc] = reform((m81_seds.fnu_obs)[where(strtrim(m81_seds[0].filter_labels, 2) eq 'SDSS_i'), *])
 sdss_map[*, *, 0] = temp
 temp[map_idc] = reform((m81_seds.fnu_obs)[where(strtrim(m81_seds[0].filter_labels, 2) eq 'SDSS_r'), *])
 sdss_map[*, *, 1] = temp
 temp[map_idc] = reform((m81_seds.fnu_obs)[where(strtrim(m81_seds[0].filter_labels, 2) eq 'SDSS_g'), *])
 sdss_map[*, *, 2] = temp


; Next, we need to match our results to the pixels within the associated region and sum the properties for all pixels
 bband_loc = m81_map[where(bband_mask eq 1 and finite(m81_map) and ksband_mask ne 1)]
 kband_loc = m81_map[where(ksband_mask eq 1 and finite(m81_map))]

 bband_idc = lonarr(n_elements(bband_loc))
 kband_idc = lonarr(n_elements(kband_loc))
 for i=0, n_elements(bband_loc)-1 do bband_idc[i] = where(bband_loc[i] eq float(strtrim(results.sed_id, 2)))
 for i=0, n_elements(kband_loc)-1 do kband_idc[i] = where(kband_loc[i] eq float(strtrim(results.sed_id, 2)))
  
; The observed luminosity
 outer_lnu_obs = total(results[bband_idc].lnu_obs, 2)
 inner_lnu_obs = total(results[kband_idc].lnu_obs, 2)

; The observed uncertainty on the luminosity
 outer_lnu_unc = total(sqrt((results[bband_idc].lnu_unc)^2.d), 2)
 inner_lnu_unc = total(sqrt((results[kband_idc].lnu_unc)^2.d), 2)

; The best-fit model luminosity
 lnu_mod_best = dblarr(n_elements(results[0].filter_labels), n_elements(results))
 for i=0, n_elements(results)-1 do lnu_mod_best[*, i] = (results[i].lnu_mod)[*, (where(results[i].lnprob eq max(results[i].lnprob)))[0]]
 outer_lnu_mod = total(lnu_mod_best[*, bband_idc], 2)
 inner_lnu_mod = total(lnu_mod_best[*, kband_idc], 2)

; The high resolution total model luminosity
 outer_lnu_hires = total(results[bband_idc].lnu_mod_hires, 2)
 inner_lnu_hires = total(results[kband_idc].lnu_mod_hires, 2)

; The high resolution stellar model luminosity
 outer_lnu_star_hires = total(results[bband_idc].lnu_starmod_hires, 2)
 inner_lnu_star_hires = total(results[kband_idc].lnu_starmod_hires, 2)

; The high resolution dust emission model luminosity
 outer_lnu_dust_hires = total(results[bband_idc].lnu_dustmod_hires, 2)
 inner_lnu_dust_hires = total(results[kband_idc].lnu_dustmod_hires, 2)

; The percentiles of the SFH of each region
 outer_sfh = total(results[bband_idc].psi, 3)
 inner_sfh = total(results[kband_idc].psi, 3)
 outer_sfh_percentiles = dblarr(5, 3)
 inner_sfh_percentiles = dblarr(5, 3)
 for i=0, 4 do outer_sfh_percentiles[i, *] = percentile(outer_sfh[i, *], [0.16, 0.5, 0.84])
 for i=0, 4 do inner_sfh_percentiles[i, *] = percentile(inner_sfh[i, *], [0.16, 0.5, 0.84])


; B-band ellipse defining values gathered from HyperLeda
 b_ra = [09, 55, 33.15]
 b_dec = [69, 03, 55.2]
 b_ra = 15 * (b_ra[0] + b_ra[1]/60. + b_ra[2]/3600.) 
 b_dec = (b_dec[0] + b_dec[1]/60. + b_dec[2]/3600.) 
 adxy, hd, b_ra, b_dec, xb, yb

 logd25 = 2.34  ; logd25 is the decimal logarithm of the length the projected major axis of a galaxy at the isophotal level 25 mag/arcsec2 in the B-band (d25 in 0.1 arcmin).
 logr25 = 0.28  ; logr25 is the axis ratio of the isophote 25 mag/arcsec2 in the  B-band for galaxies (decimal logarithm of the ratio of the lengths of major to minor axes).
 a = 10.d^logd25 * 0.1d/2 ; semi-major axis in arcmin
 q = 10.d^(-1.d*logr25)   ; axis ratio in q = b/a
 b = a * q                ; semi-minor axis in arcmin
 ab = a/60. / sxpar(hd, 'CD2_2')   ; convert to pixels using the header data
 bb = b/60. / sxpar(hd, 'CD2_2')   ; convert to pixels using the header data


; Ks-band ellipse defining values gathered from Jarrett et al 2003
 k_ra = [09, 55, 33.1]
 k_dec = [69, 03, 54.9]
 k_ra = 15 * (k_ra[0] + k_ra[1]/60. + k_ra[2]/3600.) 
 k_dec = (k_dec[0] + k_dec[1]/60. + k_dec[2]/3600.) 
 adxy, hd, k_ra, k_dec, xk, yk

 r20 = 487.6  ; The length of the projected major axis at the isophotal level 20 mag/arcsec2 in the Ks-band (arcsec).
 q = 0.51     ; The axis ratio of the isophote 20 mag/arcsec2 in the Ks-band (minor to major axes b/a).
 ; Note that we use 1/2 the Ks-band 20 mag isophote
 a = r20 / 2.           ; semi-major axis in arcsec
 b = a * q              ; semi-minor axis in arcsec
 ak = a/3.6d3 / sxpar(hd, 'CD2_2') ; convert to pixels using the header data
 bk = b/3.6d3 / sxpar(hd, 'CD2_2') ; convert to pixels using the header data



; Create a blank window with our desired dimensions
 fig_window = window(dimension=[1200, 600])

; Plot the SDSS image cropped to exclude the outer-most 35 pixels
 img = image(alog10(sdss_map[35:235,35:235, *]), $
             position=[0.01, 0.1, 0.31, 0.9], $
             aspect_ratio=1, $
             axis_style=4, $
             /current)
; Add the ellipses over the SDSS image
 ellps = ellipse(xB-35, $
                 yb-35, $
                 major=ab, $
                 minor=bb, $
                 theta=157-90, $
                 target=img, $
                 color='blue', $
                 /data, $
                 fill_background=0, $
                 thick=3)
 ellps = ellipse(xk-35, $
                 yk-35, $
                 major=ak, $
                 minor=bk, $
                 theta=-31-90, $
                 target=img, $
                 color='orange', $
                 /data, $
                 fill_background=0, $
                 thick=3)
; Add lines pointing the ellipses to the SED plots
 plyln = polyline([182, 405, 186, 405], $
                  [391, 588, 196, 348], $
                  color='blue', $
                  thick=3, $
                  /device, $
                  connectivity=[2, 0, 1, 2, 2, 3])
 plyln = polyline([210, 405, 167, 405], $
                  [337, 300, 261, 60], $
                  color='orange', $
                  thick=3, $
                  /device, $
                  connectivity=[2, 0, 1, 2, 2, 3])


; Convert the wavelengths in results to frequency for the filters and high res models
  nu = 1.d4*!lightning_cgs.clight/results[0].wave_filters
  nu_hires = 1.d4*!lightning_cgs.clight/results[0].wave_hires 

; Plot the outer region SED and high res models
 fig = objarr(4)
 fig[0] = plot(results[0].wave_hires, $
               nu_hires * outer_lnu_dust_hires, $
               position=[0.40, 0.68, 0.65, 0.98], $
               /current, $
               /ylog, $
               /xlog, $
               ytickfont_size=12, $
               xrange=[0.08, 500], $
               color='green', $
               thick=2, $
               name='Dust', $
               yrange=[1e6, 5e10], $
               xshowtext=0)
 fig[1] = plot(results[0].wave_hires, $
               nu_hires * outer_lnu_star_hires, $
               /overplot, $
               color='red', $
               thick=2, $
               name='Stellar')
 fig[2] = plot(results[0].wave_hires, $
               nu_hires * outer_lnu_hires, $
               /overplot, $
               color='gray', $
               thick=2, $
               name='Total')
 fig[3] = errorplot(results[0].wave_filters, $
                    nu * outer_lnu_obs, $
                    nu * outer_lnu_unc, $
                    /overplot, $
                    linestyle='', $
                    symbol='o', $
                    /sym_filled, $
                    name='Data')
 lgnd = legend(target=fig, $
               position=[0.41, 0.69], $
               vertical_alignment=0, $
               horizontal_alignment=0, $
               font_size=12, $
               vertical_spacing=0.01, $
               sample_width=0.07, $
               horizontal_spacing=0.03)
 fig_txt = text(0.63, 0.96, 'Outer Region', $
                font_size=14, $
                color='blue', $
                vertical_alignment=1, $
                alignment=1)
 fig_txt = text(0.35, 0.83, '$\nu L_\nu$ [$L_\odot$]', $
                font_size=12, $
                vertical_alignment=1, $
                alignment=0.5, $
                orientation=90)

; Plot the residuals for the outer region
 residuals = (outer_lnu_obs-outer_lnu_mod)/outer_lnu_unc
 fig = plot(fig[0].xrange, $
            [0, 0], $
            /current, $
            ytickinterval=3, $
            ytickfont_size=12, $
            /xlog, $
            xtickunits='scientific', $
            position=[0.40, 0.58, 0.65, 0.67], $
            xtickfont_size=12, $
            xrange=fig[0].xrange, $
            yrange=[-4.99, 4.99], $
            thick=2, $
            xtitle='Observed-Frame Wavelength [$\mu$m]')
 fig = plot(fig.xrange, $
            [1, 1], $
            /over, $
            color='light grey', $
            thick=2, $
            linestyle='--')
 fig = plot(fig.xrange, $
            [-1,-1], $
            /over, $
            color='light grey', $
            thick=2, $
            linestyle='--')
 fig = errorplot(results[0].wave_filters, $
                 residuals, $
                 replicate(1.d, n_elements(residuals)), $
                 /overplot, $
                 linestyle='', $
                 symbol='o', $
                 /sym_filled)
 fig_txt = text(0.35, 0.625, 'Residual [$\sigma$]', $
                font_size=12, $
                vertical_alignment=1, $
                alignment=0.5, $
                orientation=90)

; Plot the inner region SED and high res models
 fig = objarr(4)
 fig[0] = plot(results[0].wave_hires, $
               nu_hires * inner_lnu_dust_hires, $
               position=[0.40, 0.19, 0.65, 0.49], $
               /current, $
               /ylog, $
               /xlog, $
               ytickfont_size=12, $
               xrange=[0.08, 500], $
               color='green', $
               thick=2, $
               name='Dust', $
               yrange=[1e6, 5e10], $
               xshowtext=0)
 fig[1] = plot(results[0].wave_hires, $
               nu_hires * inner_lnu_star_hires, $
               /overplot, $
               color='red', $
               thick=2, $
               name='Stellar')
 fig[2] = plot(results[0].wave_hires, $
               nu_hires * inner_lnu_hires, $
               /overplot, $
               color='gray', $
               thick=2, $
               name='Total')
 fig[3] = errorplot(results[0].wave_filters, $
                    nu * inner_lnu_obs, $
                    nu * inner_lnu_unc, $
                    /overplot, $
                    linestyle='', $
                    symbol='o', $
                    /sym_filled, $
                    name='Data')
 lgnd = legend(target=fig, $
               position=[0.44, 0.20], $
               vertical_alignment=0, $
               horizontal_alignment=0, $
               font_size=12, $
               vertical_spacing=0.01, $
               sample_width=0.07, $
               horizontal_spacing=0.03)
 fig_txt = text(0.63, 0.47, 'Inner Region', $
                font_size=14, $
                color='orange', $
                vertical_alignment=1, $
                alignment=1)
 fig_txt = text(0.35, 0.34, '$\nu L_\nu$ [$L_\odot$]', $
                font_size=12, $
                vertical_alignment=1, $
                alignment=0.5, $
                orientation=90)

; Plot the residuals for the inner region
 residuals = (inner_lnu_obs-inner_lnu_mod)/inner_lnu_unc
 fig = plot(fig[0].xrange, $
            [0, 0], $
            /current, $
            ytickinterval=3, $
            ytickfont_size=12, $
            /xlog, $
            xtickunits='scientific', $
            position=[0.40, 0.09, 0.65, 0.18], $
            xtickfont_size=12, $
            xrange=fig[0].xrange, $
            yrange=[-4.99, 4.99], $
            thick=2, $
            xtitle='Observed-Frame Wavelength [$\mu$m]')
 fig = plot(fig.xrange, $
            [1, 1], $
            /overplot, $
            color='light grey', $
            thick=2, $
            linestyle='--')
 fig = plot(fig.xrange, $
            [-1,-1], $
            /overplot, $
            color='light grey', $
            thick=2, $
            linestyle='--')
 fig = errorplot(results[0].wave_filters, $
                 residuals, $
                 replicate(1.d, $
                 n_elements(residuals)), $
                 /overplot, $
                 linestyle='', $
                 symbol='o', $
                 /sym_filled)
 fig_txt = text(0.35, 0.135, 'Residual [$\sigma$]', $
                font_size=12, $
                vertical_alignment=1, $
                alignment=0.5, $
                orientation=90)


; Create age_bounds array for plotting
;   Allows for step function plotting
 bounds=([results[0].steps_bounds, results[0].steps_bounds])[sort([results[0].steps_bounds, results[0].steps_bounds])]
 bounds = bounds[1:-2]
 bounds[0] = (bounds[1])/10.d0

; Create SFH array for each region and uncertainty range for plotting
;   Allows for step function plotting
 SFH = dblarr(5 * 2, 2, 3)
 for k=0, 2 do for i=0, (5-1) do SFH[(2*i):(2*i)+1, 0, k] = outer_sfh_percentiles[i, k]  
 for k=0, 2 do for i=0, (5-1) do SFH[(2*i):(2*i)+1, 1, k] = inner_sfh_percentiles[i, k]  

; Plot the SFH for both regions on the same axis
 fig = FillPlot(bounds, $
                [[SFH[*, 0, 2]], [SFH[*, 0, 0]]], $
                xrange=[1e6, results[0].steps_bounds[-1]], $
                /current, $
                linestyle='', $
                fill_color='blue', $
                position=[0.73, 0.30, 0.98, 0.70], $
                /xlog, $
                /ylog, $
                xtitle='Stellar Age, t [yr]', $
                ytitle='SFR(t) [$M_\odot$ yr$^{-1}$]',$
                yrange=[3e-2, 1e1], $
                xtickfont_size=12, $
                ytickfont_size=12)
 fig = plot(bounds, $
            SFH[*, 0, 1], $
            /overplot, $
            color='blue', $
            thick=2)

 fig = FillPlot(bounds, $
                [[SFH[*, 1, 2]], $
                [SFH[*, 1, 0]]], $
                /overplot, $
                linestyle='', $
                fill_color='orange')
 fig = plot(bounds, $
           SFH[*, 1, 1], $
           /overplot, $
           color='orange', $
           thick=2)
 fig_txt = text(0.74, 0.68, 'Outer Region', $
                font_size=14, $
                color='blue', $
                vertical_alignment=1)
 fig_txt = text(0.74, 0.64, 'Inner Region', $
                font_size=14, $
                color='orange', $
                vertical_alignment=1)

 return, fig

end