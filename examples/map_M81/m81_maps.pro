function m81_maps, sfr_map, mass_map, ssfr_map, tauv_diff_map, ltir_map, pvalue_map
;+
; Name
; ----
;   M81_MAPS
;
; Purpose
; -------
;   Generates property maps of M81. Input maps must
;   be formatted and cropped as show in Lightning's online documentation.
;
; Calling Sequence
; ----------------
;   ::
;
;       fig = m81_maps(sfr_map, mass_map, ssfr_map, tauv_diff_map, ltir_map, pvalue_map)
;
; Input
; -----
;   ``sfr_map`` : float or double array(N, M)
;       The SFR map [:math:`M_\odot\ {\rm yr}^{-1}`].
;   ``mass_map`` : float or double array(N, M)
;       The stellar mass map [:math:`M_\odot`].
;   ``ssfr_map`` : float or double array(N, M)
;       The sSFR map [:math:`{\rm yr}`].
;   ``tauv_diff_map`` : float or double array(N, M)
;       The V-band optical depth map.
;   ``ltir_map`` : float or double array(N, M)
;       The LTIR map [:math:`L_\odot`].
;   ``pvalue_map`` : float or double array(N, M)
;       The pvalue map.
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


; Set the inputs that are common to each panel so we can feed them in as keyword inheritance
 img_keywords = {aspect_ratio: 1, $
                 axis_style: 4, $
                 background_color: 'black', $
                 current: 1}
 cb_keywords = {border_on: 1, $
                font_size: 16, $
                orientation: 1, $
                textpos: 1, $
                tickdir: 1}
 txt_keywords = {color: 'white', $
                 font_size: 20, $
                 vertical_alignment: 1}


; Create a blank window with our desired dimensions
 img = window(dimension=[1200, 864])


; Plot the SFR map in log scale with a colorbar within the created window 
 img = image(alog10(sfr_map), $
             position=[0.01, 0.53, 0.23, 0.98], $
             rgb_table=3, $
             _extra=img_keywords)
 cb = colorbar(target=img, $
               position=[0.235, 0.535, 0.25, 0.975], $
               title='$log_{10}(SFR$ [$M_\odot$ $yr^{-1}$])', $
               _extra=cb_keywords)
 img_txt = text(0.02, 0.97, '(a)', _extra=txt_keywords)

; Plot the stellar mass map in log scale with a colorbar within the created window 
 img = image(alog10(mass_map), $
             position=[0.345, 0.53, 0.565, 0.98], $
             rgb_table=3, $
             _extra=img_keywords)
 cb = colorbar(target=img, $
               position=[0.57, 0.535, 0.585, 0.975], $
               title='$log_{10}(M_\star$ [$M_\odot$])', $
               _extra=cb_keywords)
 img_txt = text(0.355, 0.97, '(b)', _extra=txt_keywords)

; Plot the sSFR map (SFR/mass) in log scale with a colorbar within the created window 
 img = image(alog10(ssfr_map), $
             position=[0.68, 0.53, 0.90, 0.98], $
             rgb_table=72, $
             _extra=img_keywords)
 cb = colorbar(target=img, $
               position=[0.905, 0.535, 0.92, 0.975], $
               title='$log_{10}(sSFR$ [$yr^{-1}$])', $
               tickinterval=1, $
               _extra=cb_keywords)
 img_txt = text(0.69, 0.97, '(c)', _extra=txt_keywords)

; Plot tauV_DIFF with a colorbar within the created window 
 img = image(tauv_diff_map, $
             position=[0.01, 0.02, 0.23, 0.47], $
             rgb_table=3, $
             _extra=img_keywords)
 cb = colorbar(target=img, $
               position=[0.235, 0.025, 0.25, 0.465], $
               title='$\tau_{DIFF, V}$', $
               _extra=cb_keywords)
 img_txt = text(0.02, 0.46, '(d)', _extra=txt_keywords)

; Plot LTIR in log scale with a colorbar within the created window 
 img = image(alog10(ltir_map), $
             position=[0.345, 0.02, 0.565, 0.47], $
             rgb_table=3, $
             _extra=img_keywords)
 cb = colorbar(target=img, $
               position=[0.57, 0.025, 0.585, 0.465], $
               title='$log_{10}(L_{TIR}$ [$L_\odot$])', $
               _extra=cb_keywords)
 img_txt = text(0.355, 0.46, '(e)', _extra=txt_keywords)

; Plot p-value in log scale with a colorbar within the created window 
 img = image(alog10(pvalue_map), $
             position=[0.68, 0.02, 0.90, 0.47], $
             rgb_table=3, $
             _extra=img_keywords)
 cb = colorbar(target=img, $
               position=[0.905, 0.025, 0.92, 0.465], $
               title='$log_{10}(p-value)$', $
               _extra=cb_keywords)
 img_txt = text(0.69, 0.46, '(f)', _extra=txt_keywords)

 return, img

end