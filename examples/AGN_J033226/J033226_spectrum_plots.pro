function J033226_spectrum_plots, xray_res
;+
; Name
; ----
;   J033226_SPECTRUM_PLOTS
;
; Purpose
; -------
;   Create a plot showing the SED fit for J033226.49-274035.5
;   (with the X-ray model; the other fit is not noticeably different)
;   and the folded X-ray spectral fit.
;
; Input
; -----
;   xray_res : structure
;       An IDL structure containing the fit to J033226 with
;       the X-ray model
;
; Output
; ------
;   Returns a IDL graphics object containing the plots
;
;-

    compile_opt IDL2

    lightning_constants

    ;--------------------------
    ; Construct UV-IR Models
    ;--------------------------
    ; Low-res
    wave_obs = xray_res.WAVE_FILTERS
    nu_obs = !lightning_cgs.clight * 1e4 / wave_obs
    lnu_obs = xray_res.LNU_OBS
    lnu_unc = xray_res.LNU_UNC
    minchi2 = min(xray_res.CHI2, bestfit)
    lnu_mod = xray_res.LNU_MOD[*,bestfit]
    lnu_unc_total = sqrt(lnu_unc^2 + (xray_res.MODEL_UNC * lnu_mod)^2)
    delchi = (lnu_obs - lnu_mod) / lnu_unc_total

    ; High res
    wave_hires = xray_res.WAVE_HIRES
    nu_hires = !lightning_cgs.clight * 1e4 / wave_hires
    lnu_total_hires = xray_res.LNU_MOD_HIRES[*,0]
    lnu_stars_hires = xray_res.LNU_STARMOD_HIRES[*,0]
    lnu_AGN_hires = xray_res.LNU_AGNMOD_HIRES[*,0]
    lnu_dust_hires = xray_res.LNU_DUSTMOD_HIRES[*,0]

    ;--------------------------
    ; Construct X-ray Models
    ;--------------------------
    ; Low-res
    xray_bandpass = xray_res.XRAY_BANDPASS
    xray_bandpass_mean = 0.5 * (xray_bandpass[0,*] + xray_bandpass[1,*])
    xray_bandpass_width = xray_bandpass[1,*] - xray_bandpass[0,*]
    xray_counts = xray_res.NET_COUNTS
    xray_counts_unc = xray_res.NET_COUNTS_UNC
    xray_counts_mod = xray_res.XRAY_COUNTS_MOD
    xray_exposure = xray_res.XRAY_EXPOSURE
    xray_countrate = xray_counts / xray_bandpass_width / xray_exposure
    xray_countrate_unc = xray_counts_unc / xray_bandpass_width / xray_exposure
    xray_countrate_mod = xray_counts_mod / xray_bandpass_width / xray_exposure
    xray_delchi = (xray_counts - xray_counts_mod) / xray_counts_unc

    ; For the legend
    lines = []

    p1 = plot(wave_hires, nu_hires * lnu_dust_hires,$
              /XLOG,$
              /YLOG,$
              thick=2,$
              color='green',$
              name='Dust')
    lines = [lines, p1]

    p1 = plot(wave_hires, nu_hires * lnu_AGN_hires,$
              /OVERPLOT,$
              /XLOG,$
              /YLOG,$
              thick=2,$
              color='orange',$
              name='AGN')
    lines = [lines, p1]

    p1 = plot(wave_hires, nu_hires * lnu_stars_hires,$
              /OVERPLOT,$
              /XLOG,$
              /YLOG,$
              thick=2,$
              color='red',$
              name='Stellar')
    lines = [lines, p1]

    p1 = plot(wave_hires, nu_hires * lnu_total_hires,$
              /OVERPLOT,$
              /XLOG,$
              /YLOG,$
              thick=2,$
              color='gray',$
              name='Total')
    lines = [lines, p1]

    p1 = errorplot(wave_obs, nu_obs * lnu_obs, nu_obs * lnu_unc_total,$
                   /OVERPLOT,$
                   /XLOG,$
                   /YLOG,$
                   linestyle='',$
                   sym='o',$
                   sym_color='k',$
                   /sym_filled,$
                   name='Data')
    lines = [lines, p1]

    p1['axis1'].showtext = 0
    p1['axis3'].showtext = 1

    p1.xrange = [0.09, 1000]
    p1.yrange = [1e10, 3e12]

    p1.xshowtext = 0
    p1.ytitle = '$\nuL_{\nu}$ $[L_{\odot}]$'

    position = p1.position
    figwidth = position[2] - position[0]
    figheight = position[3] - position[1]
    residual_position = [position[0] + figwidth * 0.48 + 0.01,$
                         position[1],$
                         position[0] + 2 * figwidth * 0.48 + 0.01,$
                         position[1] + figheight * 4/15]
    xray_position = [position[0],$
                     position[1] + figheight * 4/15 + 0.02,$
                     position[0] + figwidth * 0.48,$
                     position[3]]
    xray_residual_position = [position[0],$
                              position[1],$
                              position[0] + figwidth * 0.48,$
                              position[1] + figheight * 4/15]
    position = [position[0] + figwidth * 0.48 + 0.01,$
                position[1] + figheight * 4/15 + 0.02,$
                position[0] + 2 * figwidth * 0.48 + 0.01,$
                position[3]]

    p1.position = position
    l1 = legend(target=lines, transp=30, POSITION=[xray_position[2]-0.01, xray_position[1]+0.01], $
                VERTICAL_ALIGNMENT=0, HORIZONTAL_ALIGNMENT=1,$
                font_size=8)

    t1 = text(position[2] - 0.02, position[3] - 0.03, xray_res.SED_ID, VERTICAL_ALIGNMENT=1, ALIGNMENT=1)
    t2 = text(position[2] - 0.02, position[3] - 0.06, 'z = ' + strtrim(string(xray_res.REDSHIFT, format='(F5.2)'),2), VERTICAL_ALIGNMENT=1, ALIGNMENT=1)

    p2 = plot(p1.xrange, [0, 0], /CURRENT, ytitle='Residual [$\sigma$]', $
              /XLOG, xtickunits='scientific', position=residual_position, $
              xrange=p1.xrange, yrange=[-4.5, 4.5], thick=2)
    p2 = errorplot(wave_obs, delchi, 1 + dblarr(n_elements(delchi)),$
                   /OVERPLOT,$
                   /XLOG,$
                   linestyle='',$
                   sym='o',$
                   sym_color='k',$
                   yshowtext=0,$
                   /sym_filled)

    p2.xtitle = 'Observed-Frame Wavelength $[\mu m]$'

    p3 = plot(reverse(xray_bandpass_mean),$
              reverse(xray_countrate_mod),$
              /CURRENT,$
              /XLOG,$
              /YLOG,$
              xshowtext=0,$
              xrange=[8.0,0.5],$
              xtickvalues=[8,5,3,2,1],$
              ;xsubticklen=0,$
              xminor=1,$
              color='gray',$
              thick=2,$
              ytitle='Counts $s^{-1}$ $keV^{-1}$',$
              position=xray_position)

    p3 = errorplot(xray_bandpass_mean,$
                   xray_countrate,$
                   xray_countrate_unc,$
                   /OVERPLOT,$
                   /XLOG,$
                   /YLOG,$
                   ;xshowtext=0,$
                   ;xrange=[0.5,8.0],$
                   linestyle='',$
                   sym='o',$
                   sym_color='k',$
                   /sym_filled)

    p4 = plot(p3.xrange, [0, 0], /CURRENT, ytitle='Residual [$\sigma$]', $
              /XLOG, position=xray_residual_position, $
              xrange=p3.xrange, yrange=[-4.5, 4.5], xtickvalues=[8,5,3,2,1], xminor=1,$
              thick=2)

    p4 = errorplot(xray_bandpass_mean,$
                   xray_delchi,$
                   1 + dblarr(n_elements(xray_delchi)),$
                   /OVERPLOT,$
                   /XLOG,$
                   linestyle='',$
                   sym='o',$
                   sym_color='k',$
                   /sym_filled)

    p4.xtitle = 'Observed-Frame Energy [keV]'

    return, p1

end
