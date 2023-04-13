pro compare_cigale_properties, lghtng, cigale_results
;+
; Name
; ----
;   COMPARE_CIGALE_PROPERTIES
;
; Purpose
; -------
;   Generates and print AGN property comparison of J033226. Input tables
;   must be the post-processed outputs from Lightning and CIGALE.
;
; Calling Sequence
; ----------------
;   ::
;
;       compare_cigale_properties, lghtng, cigale_results)
;
; Input
; -----
;   ``lghtng`` : structure
;       The Lightning post-processed output.
;   ``cigale_results`` : structure
;       The CIGALE results.
;
; Modification History
; --------------------
;   - 2023/04/12: Created (Keith Doore)
;-
 On_error, 2
 compile_opt IDL2

 ; CIGALE gives luminosities in Watts, create variable to convert to Lsun
 ; Conversion from Watts to Lsun
 ;     (erg/s)/W   (ergs/sec)/Lsun
 conv = 1.d7 / !lightning_cgs.Lsun


 ; Calculate the LAGN with created function
 ; Read in the ARF and Lightning configuration file to calculate LAGN
 arf =  mrdfits('xray/J033226_summed.arf', 1, /silent)
 restore, 'xray/lightning_output/lightning_configure.sav'
 AGN_model_lum = calc_integrated_AGN_luminosity(lghtng, config, arf)


 ; Calculate alphaOX
 hc = (!lightning_cgs.hplanck / !lightning_cgs.keV) * (1e4 * !lightning_cgs.clight) ; h * c in keV * micron
 ; E_obs = hc / wave
 L2kev = dblarr(n_elements(AGN_model_lum.L2500))
 ; Need to calculate at rest-frame 2keV. So, shift model wavelengths blueward.
 for i=0, n_elements(AGN_model_lum.L2500)-1 do L2kev[i] = interpol((lghtng.LNU_XRAYMOD_AGN_UNABS_HIRES)[*,i], $
                                                                   lghtng.WAVE_XRAYMOD_HIRES/(1 + lghtng.redshift), hc/2.d0)
 alpha_ox = -0.3838 * alog10(AGN_model_lum.L2500/L2kev)


 ; Calculate fracAGN
 keep = where(lghtng.wave_hires ge 5 and lghtng.wave_hires le 1000)
 Lbol_ir = dblarr(n_elements(AGN_model_lum.L2500))
 Lagn_ir = dblarr(n_elements(AGN_model_lum.L2500))
 for i=0,n_elements(AGN_model_lum.L2500)-1 do Lbol_ir[i] = -1.d*trapez((lghtng.lnu_mod_hires)[keep, i], $
                                                                       1.d4 * !lightning_cgs.clight/(lghtng.wave_hires)[keep]) 
 for i=0,n_elements(AGN_model_lum.L2500)-1 do Lagn_ir[i] = -1.d*trapez((lghtng.LNU_AGNMOD_HIRES)[keep, i], $
                                                                       1.d4 * !lightning_cgs.clight/(lghtng.wave_hires)[keep])
 frac_AGN = Lagn_ir/lbol_ir

 ; Convert cosi to degrees
 i_AGN = acos(lghtng.AGN_COSI) / !dtor


 ; Calculate percentiles of each derived property
 inc_percent = percentile(i_agn, [0.5, 0.16, 0.84])
 lagn_percent = percentile(alog10(AGN_model_lum.LBOL_AGN_MODEL), [0.5, 0.16, 0.84])
 fracagn_percent = percentile(frac_AGN, [0.5, 0.16, 0.84])
 aox_percent = percentile(alpha_ox, [0.5, 0.16, 0.84])
 
 ; Adjust percentiles to be +/- the median
 inc_percent[1] = inc_percent[0]-inc_percent[1]
 inc_percent[2] = inc_percent[2]-inc_percent[0]
 lagn_percent[1] = lagn_percent[0]-lagn_percent[1]
 lagn_percent[2] = lagn_percent[2]-lagn_percent[0]
 fracagn_percent[1] = fracagn_percent[0]-fracagn_percent[1]
 fracagn_percent[2] = fracagn_percent[2]-fracagn_percent[0]
 aox_percent[1] = aox_percent[0]-aox_percent[1]
 aox_percent[2] = aox_percent[2]-aox_percent[0]

 ; Place the properties in to nice formatted strings: median (+84% / -16%)
 i_l = strtrim(string(inc_percent[0], f='(F4.1)'), 2)+' (+'+strtrim(string(inc_percent[2], f='(F3.1)'), 2)+$
       ' / -'+strtrim(string(inc_percent[1], f='(F3.1)'), 2)+')'
 lagn_l = strtrim(string(lagn_percent[0], f='(F5.2)'), 2)+' (+'+strtrim(string(lagn_percent[2], f='(F4.2)'), 2)+$
          ' / -'+strtrim(string(lagn_percent[1], f='(F4.2)'), 2)+')'
 frac_l = strtrim(string(fracagn_percent[0], f='(F4.2)'), 2)+' (+'+strtrim(string(fracagn_percent[2], f='(F4.2)'), 2)+$
          ' / -'+strtrim(string(fracagn_percent[1], f='(F4.2)'), 2)+')'
 a_ox_l = strtrim(string(aox_percent[0], f='(F5.2)'), 2)+' (+'+strtrim(string(aox_percent[2], f='(F4.2)'), 2)+$
          ' / -'+strtrim(string(aox_percent[1], f='(F4.2)'), 2)+')'

 ; Place into strings: median (+/- stddev)
 i_c = strtrim(string(cigale_results.BAYES_AGN_I, f='(F4.1)'), 2)+' (+/- '+strtrim(string(cigale_results.BAYES_AGN_I_err, f='(F4.1)'), 2)+')'
 ; CIGALE output LAGN in units of Watts
 lagn_c = strtrim(string(alog10(cigale_results.BAYES_AGN_LUMINOSITY*conv), f='(F5.2)'), 2)+' (+/- '+$
          strtrim(string(cigale_results.BAYES_AGN_LUMINOSITY_err/(cigale_results.BAYES_AGN_LUMINOSITY* alog(10)), f='(F4.2)'), 2)+')'
 frac_c = strtrim(string(cigale_results.BAYES_AGN_FRACAGN, f='(F4.2)'), 2)+' (+/- '+$
          strtrim(string(cigale_results.BAYES_AGN_FRACAGN_err, f='(F4.2)'), 2)+')'
 a_ox_c = strtrim(string(cigale_results.BAYES_XRAY_ALPHA_OX, f='(F5.2)'), 2)+' (+/- '+$
          strtrim(string(cigale_results.BAYES_XRAY_ALPHA_OX_err, f='(F4.2)'), 2)+')'

stop
 ; Print strings
 print, '//Inclination Comparison//'
 print, 'Lightning iAGN: '+i_l
 print, 'CIGALE iAGN:    '+i_c
 print, ''
 print, '//log10(LAGN) Comparison//'
 print, 'Lightning LAGN: '+lagn_l
 print, 'CIGALE LAGN:    '+lagn_c
 print, ''
 print, '//fracAGN Comparison//'
 print, 'Lightning fracAGN: '+frac_l
 print, 'CIGALE fracAGN:    '+frac_c
 print, ''
 print, '//alpha_OX Comparison//'
 print, 'Lightning alpha_OX: '+a_ox_l
 print, 'CIGALE alpha_OX:    '+a_ox_c
 print, ''

end