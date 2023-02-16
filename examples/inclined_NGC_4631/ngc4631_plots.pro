function ngc4631_plots, calz, calz_x, doore, doore_x
;+
; Name
; ----
;   NGC4631_PLOTS
;
; Purpose
; -------
;   Generates SFR posterior comparison plots of NGC 4631. Input tables
;   must be the post-processed outputs from Lightning.
;
; Calling Sequence
; ----------------
;   ::
;
;       fig = ngc4631_plots(calz, calz_x, doore, doore_x)
;
; Input
; -----
;   ``calz`` : structure
;       The Lightning post-processed output for the Calzetti fit without X-rays.
;   ``calz_x`` : structure
;       The Lightning post-processed output for the Calzetti fit with X-rays.
;   ``doore`` : structure
;       The Lightning post-processed output for the Doore+21 fit without X-rays.
;   ``doore_x`` : structure
;       The Lightning post-processed output for the Doore+21 fit with X-rays.
;
; Output
; ------
;   ``fig`` : object
;       The plot object containing the figure.
;
; Modification History
; --------------------
;   - 2023/02/10: Created (Keith Doore)
;-
 On_error, 2
 compile_opt IDL2


; Let's combine all of the SFRs from each fit into a single array
;  We can do this with an execute statement since each structure has the same tags
 sfr = !null
 fit = ['calz', 'calz_x', 'doore','doore_x']
 for i=0,3 do begin
   void = execute('x = '+fit[i])
   ; SFR of last 100 Myr is weighted by each bin size (0-10 and 10-100 Myr)
   sfr = [[sfr], [0.1d * reform((x.psi)[0, *]) + 0.9d * reform((x.psi)[1, *])]]
 endfor

; Create plot window. We need a wide window to accommodate the histograms SED
 fig = window(dimension=[1200, 600])
; Set the colors for each of the fit models
 model_colors=['red', 'orange', 'green', 'blue']

; Determine histogram for each distribution and normalize
 for i=0, 3 do begin
   ; Binsize determined using Scott's normal reference rule
   binsize = 3.49d * stddev(sfr[*, i]) / (double(n_elements(sfr[*, i])))^(1/3.d)
   nbins = ceil((max(sfr[*, i]) - min(sfr[*, i]))/binsize)

   pdf = histogram(sfr[*, i], nbins=nbins, locations=binloc)
   ; Area normalize
   peak = total(pdf) * (binloc[1]-binloc[0])
   pdf = double(pdf)/peak

   if i eq 2 then xshowtext = !null else xshowtext = 0
   if i eq 0 then row = 0 else if i eq 2 then row = 1
   if i eq 0 then yrange = [0, 0.43] else if i eq 2 then yrange = [0, 0.19]
   if i eq 0 or i eq 2 then begin
     overplot = 0 & current = 1
   endif else begin
     overplot = 1 & current = 0
   endelse

   ; Plot the histogram
   fig = plot(binloc, $
              pdf, $
              /stair, $
              xtitle='SFR [$M_\odot$ $yr^{-1}$]', $
              ytitle='P(SFR)', $
              thick=2, $
              position=[0.08, 0.54 - row * 0.88/2, 0.3, 0.98 - row * 0.88/2], $
              color=model_colors[i], $
              xtickfont_size=14, $
              ytickfont_size=14, $
              xrange=[0, 13], $
              xshowtext=((i>2)-2), $
              yrange=yrange, $
              overplot=overplot, $
              current=current, $
              ytickinterval=0.1)
   fig['axis2'].showtext = 0
 endfor

 ; Plot the high resolution models
 nu_highres      = 1.d4 * !lightning_cgs.clight/doore_x.wave_hires
 nu_highres_xray = 1.d4 * !lightning_cgs.clight/doore_x.wave_xraymod_hires

 hires_fig = objarr(4)
 hires_fig[3] = plot(calz.wave_hires, $
                     nu_highres * (calz.lnu_mod_hires)[*, 0], $
                     /current, $
                     color='red', $
                     thick=2, $
                     name='Calzetti no X-rays')
 hires_fig[2] = plot(calz_x.wave_hires, $
                     nu_highres * (calz_x.lnu_mod_hires)[*, 0], $
                     /overplot, $
                     color='orange', $
                     thick=2, $
                     name='Calzetti with X-rays')
 hires_fig[2] = plot(calz_x.wave_xraymod_hires, $
                     nu_highres_xray * calz_x.lnu_xraymod_hires[*, 0], $
                     /overplot, $
                     color='orange', $
                     thick=2, $
                     name='Calzetti with X-rays')
 hires_fig[1] = plot(doore.wave_hires, $
                     nu_highres * (doore.lnu_mod_hires)[*, 0], $
                     /overplot, $
                     color='green', $
                     thick=2, $
                     name='Doore no X-rays')
 hires_fig[0] = plot(doore_x.wave_hires, $
                     nu_highres * (doore_x.lnu_mod_hires)[*, 0], $
                     /overplot, $
                     color='blue', $
                     thick=2, $
                     name='Doore with X-rays')
 hires_fig[0] = plot(doore_x.wave_xraymod_hires, $
                     nu_highres_xray * doore_x.lnu_xraymod_hires[*, 0], $
                     /overplot, $
                     color='blue', $
                     thick=2, $
                     name='Doore with X-rays')

; Create a legend using the high resolution models since we gave each plot a name
 lgnd = legend(target=hires_fig, $
               transparency=30, $
               position=[0.4, 0.83], $
               vertical_alignment=1, $
               horizontal_alignment=0, $
               font_size=14)

; Plot the observed UV-to-IR photometric data points
 nu = 1.d4 * !lightning_cgs.clight/doore_x.wave_filters
 fig = errorplot(doore_x.wave_filters, $
                 nu * doore_x.lnu_obs, $
                 nu * doore_x.lnu_unc, $
                 /overplot, $
                 /xlog, $
                 /ylog, $
                 linestyle='', $
                 name='Data', $
                 xrange=[1.001e-4, 1300], $
                 yrange=[1e4, 5e10], $
                 axis_style=1, $
                 symbol='o', $
                 /sym_filled, $
                 sym_size=1, $
                 xtitle='Observed-Frame Wavelength [$\mu$m]', $
                 ytitle='$\nu L_\nu$ [$L_\odot$]', $
                 position=[0.38, 0.1, 0.98, 0.9], $
                 ytickfont_size=14, $
                 xtickunits='scientific', $
                 xtickfont_size=14)

; Convert the X-ray bandpasses to wavelengths
 xray_nu = doore_x.xray_bandpass / (!lightning_cgs.hplanck / !lightning_cgs.keV)
 xray_wave = 1.d4*!lightning_cgs.clight/xray_nu
 xray_mean_wave = mean(xray_wave, dim=1)
 xray_mean_nu = 1.d4*!lightning_cgs.clight/xray_mean_wave
 xray_wave_unc = dblarr(size(doore_x.xray_bandpass, /dimension))
 xray_wave_unc[0, *] = xray_mean_wave - xray_wave[1, *]
 xray_wave_unc[1, *] = xray_wave[0, *] - xray_mean_wave 

; Plot the observed X-ray data
 fig_xray = errorplot(xray_mean_wave, $
                      xray_mean_nu * doore_x.lnu_xray_obs, $
                      xray_wave_unc, $
                      xray_mean_nu * doore_x.lnu_xray_unc, $
                      /overplot, $
                      linestyle='', $
                      symbol='o', $
                      /sym_filled, $
                      sym_size=1)

; Add upper axis in energy
 yaxis = axis('Y', $
              location='right', $
              showtext=0, $
              target=fig)
 xray_xrange=(1.d4 * !lightning_cgs.clight * (!lightning_cgs.hplanck / !lightning_cgs.keV)) / fig.xrange
 xpos = fig.position
 fig_temp = plot(xray_mean_wave, $
                 xray_mean_nu * doore_x.lnu_xray_obs, $
                 /xlog, $
                 /ylog, $
                 xrange=xray_xrange, $
                 ytransparency=100, $
                 position=[xpos[0],xpos[3],xpos[2],xpos[3]+(xpos[3]-xpos[1])], $
                 xtickdir=1, $
                 xtextpos=1, $
                 /current, $
                 axis_style=1, $
                 yshowtext=0, $
                 xtickunits='scientific', $
                 /nodata, $
                 xtitle='Observed-Frame Energy [keV]', $
                 xtickfont_size=14)


; Add histogram of inclination for the Doore+21 fits
 fit = ['doore','doore_x']
 for i=0, 1 do begin
   void = execute('x = '+fit[i])
   ; Binsize determined using Scott's normal reference rule
   binsize = 3.49d * stddev(1-x.cosi) / (double(n_elements(x.cosi)))^(1/3.d)
   nbins = ceil((max(1-x.cosi) - min(1-x.cosi))/binsize)

   pdf = histogram(1-x.cosi, nbins=nbins, locations=binloc)
   ; Area normalize
   peak = total(pdf) * (binloc[1]-binloc[0])
   pdf = double(pdf)/peak

   if i eq 0 then begin
     overplot = 0 & current = 1
   endif else begin
     overplot = 1 & current = 0
   endelse

   ; Plot the histogram
   fig = plot(binloc, $
              pdf, $
              /stair, $
              xtitle='1 - cos(i)', $
              ytitle='P(1 - cos i)', $
              thick=2, $
              position=[0.69, 0.21, 0.91, 0.63], $
              color=model_colors[i+2], $
              xtickfont_size=14, $
              ytickfont_size=14, $
              xrange=[0, 1], $
              ;yrange=[0, 0.43], $
              overplot=overplot, $
              current=current)
 endfor

 return, fig

end