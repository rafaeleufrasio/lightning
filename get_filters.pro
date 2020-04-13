pro get_filters,filter_labels,nu,Filters,filters_dir=filters_dir,time=time,Lnu=Lnu,$
                mean_nu=mean_nu,mean_Lnu=mean_Lnu,sigma_nu=sigma_nu,             $
                mean_wave=mean_wave,mean_Lwave=mean_Lwave,sigma_wave=sigma_wave, $
                wave=wave,Lwave=Lwave,plot_filters=plot_filters,                 $
                plot_analysis=plot_analysis,read_filters=read_filters
;Input
;  Filter_labels
;  nu
;  filters_dir
;
;Optional inputs
;  time & nu & Lnu[time,nu]
;
;Outputs
;  Filters = Filters[filters,wave]
;  mean_nu = mean_nu[time,filters]
;  mean_Lnu = mean_Lnu[time,filters]
;  sigma_nu = mean_Lnu[time,filters]
;  mean_wave = mean_wave[time,filters]
;  mean_Lwave = mean_Lwave[time,filters]
;  sigma_wave = mean_Lwave[time,filters]
;  wave = wave[wave]
;  Lwave = Lwave[time,wave]

if (n_elements(filter_labels) eq 0) then begin print,'Please specify filter labels.' & end
if (n_elements(filters_dir) eq 0) then filters_dir='~/Filters/'
if (n_elements(read_filters) eq 0) then read_filters=1
if (n_elements(Lnu) eq 0) then analysis=0 else analysis=1
if (n_elements(plot_filters) eq 0) then plot_filters=0
if (n_elements(plot_analysis) eq 0) then plot_analysis=0

Nfilters = (size(filter_labels))[1]
Nwave = (size(nu))[1]

if read_filters then Filters=dblarr(Nfilters,Nwave) else Filters=double(Filters)

wave = 1.d4*!cv.clight/nu

if analysis then begin
  Lwave = 1.d4*!cv.clight*Lnu/wave^2       ;Lwave in Lsun/micron
  Ntime = (size(Lnu))[1]
  mean_nu = dblarr(Ntime,Nfilters)
  mean_Lnu = dblarr(Ntime,Nfilters)
  mean_nu2 = dblarr(Ntime,Nfilters)
  mean_Lwave=dblarr(Ntime,Nfilters)
  mean_wave = dblarr(Ntime,Nfilters)
  mean_wave2 = dblarr(Ntime,Nfilters)
endif

for f=0,(Nfilters-1) do begin
  if read_filters then begin
    case filter_labels[f] of
      'GALEX_FUV':     trans_curve = filters_dir+'GALEX/FUV.txt'
      'GALEX_NUV':     trans_curve = filters_dir+'GALEX/NUV.txt'
      'UVOT_W2':       trans_curve = filters_dir+'Swift/UVW2.txt'
      'UVOT_M2':       trans_curve = filters_dir+'Swift/UVM2.txt'
      'UVOT_W1':       trans_curve = filters_dir+'Swift/UVW1.txt'
      'UVOT_U':        trans_curve = filters_dir+'Swift/U.txt'
      'UVOT_B':        trans_curve = filters_dir+'Swift/B.txt'
      'UVOT_V':        trans_curve = filters_dir+'Swift/V.txt'
      'ACS_F435W':     trans_curve = filters_dir+'HST/ACS/wfc_F435W.txt'
      'ACS_F475W':     trans_curve = filters_dir+'HST/ACS/wfc_F475W.txt'
      'ACS_F555W':     trans_curve = filters_dir+'HST/ACS/wfc_F555W.txt'
      'ACS_F606W':     trans_curve = filters_dir+'HST/ACS/wfc_F606W.txt'
      'ACS_F658N':     trans_curve = filters_dir+'HST/ACS/wfc_F658N.txt'
      'ACS_F775W':     trans_curve = filters_dir+'HST/ACS/wfc_F775W.txt'
      'ACS_F814W':     trans_curve = filters_dir+'HST/ACS/wfc_F814W.txt'
      'ACS_F850LP':    trans_curve = filters_dir+'HST/ACS/wfc_F850LP.txt'
      'ACS_B':         trans_curve = filters_dir+'HST/ACS/ACS_B.txt'
      'ACS_V':         trans_curve = filters_dir+'HST/ACS/ACS_V.txt'
      'ACS_I':         trans_curve = filters_dir+'HST/ACS/ACS_I.txt'
      'ACS_Z':         trans_curve = filters_dir+'HST/ACS/ACS_Z.txt'
      'WFC3_F098M':    trans_curve = filters_dir+'HST/WFC3/WFC3_F098M.txt'
      'WFC3_F105W':    trans_curve = filters_dir+'HST/WFC3/WFC3_F105W.txt'
      'WFC3_F125W':    trans_curve = filters_dir+'HST/WFC3/WFC3_F125W.txt'
      'WFC3_F140W':    trans_curve = filters_dir+'HST/WFC3/WFC3_F140W.txt'
      'WFC3_F160W':    trans_curve = filters_dir+'HST/WFC3/WFC3_F160W.txt'
      'CTIO_U':        trans_curve = filters_dir+'CTIO/CTIO_U.txt'
      'VIMOS_U':       trans_curve = filters_dir+'VIMOS/VIMOS_U.txt'
      'VLT_U':         trans_curve = filters_dir+'VLT/FORS1_1033_U_BESS.txt'
      'VLT_B':         trans_curve = filters_dir+'VLT/FORS1_1034_B_BESS.txt'
      'VLT_V':         trans_curve = filters_dir+'VLT/FORS1_1035_V_BESS.txt'
      'VLT_R':         trans_curve = filters_dir+'VLT/FORS1_1036_R_BESS.txt'
      'VLT_I':         trans_curve = filters_dir+'VLT/FORS1_1037_I_BESS.txt'
      'KECK_g':        trans_curve = filters_dir+'Keck/Keck_g.txt'
      'KECK_Rs':       trans_curve = filters_dir+'Keck/Keck_Rs.txt'
      'KPNO_U':        trans_curve = filters_dir+'KITT_PEAK/KPNO_U.txt'
      'KP_U':          trans_curve = filters_dir+'KITT_PEAK/k1001bp_jul04.txt'
      'KP_B':          trans_curve = filters_dir+'KITT_PEAK/k1002bp_jul04.txt'
      'KP_V':          trans_curve = filters_dir+'KITT_PEAK/k1003bp_jul04.txt'
      'KP_R':          trans_curve = filters_dir+'KITT_PEAK/k1004bp_aug04.txt'
      'KP_I':          trans_curve = filters_dir+'KITT_PEAK/k1005bp_aug04.txt'
      'SDSS_u':        trans_curve = filters_dir+'SDSS/SDSS_u.txt'
      'SDSS_g':        trans_curve = filters_dir+'SDSS/SDSS_g.txt'
      'SDSS_r':        trans_curve = filters_dir+'SDSS/SDSS_r.txt'
      'SDSS_i':        trans_curve = filters_dir+'SDSS/SDSS_i.txt'
      'SDSS_z':        trans_curve = filters_dir+'SDSS/SDSS_z.txt'
      'SUBARU_B':      trans_curve = filters_dir+'SUBARU/SUBARU_B.txt'
      'SUBARU_V':      trans_curve = filters_dir+'SUBARU/SUBARU_V.txt'
      'SUBARU_Rc':     trans_curve = filters_dir+'SUBARU/SUBARU_Rc.txt'
      'SUBARU_I':      trans_curve = filters_dir+'SUBARU/SUBARU_I.txt'
      'SUBARU_zprime': trans_curve = filters_dir+'SUBARU/SUBARU_zprime.txt'
      'SUBARU_J':      trans_curve = filters_dir+'SUBARU/SUBARU_J.txt'
      'SUBARU_H':      trans_curve = filters_dir+'SUBARU/SUBARU_H.txt'
      'SUBARU_Ks':     trans_curve = filters_dir+'SUBARU/SUBARU_Ks.txt'
      '2MASS_J':       trans_curve = filters_dir+'2MASS/J_BAND_RSR.txt'
      '2MASS_H':       trans_curve = filters_dir+'2MASS/H_BAND_RSR.txt'
      '2MASS_Ks':      trans_curve = filters_dir+'2MASS/Ks_BAND_RSR.txt'
      'ISAAC_Ks':      trans_curve = filters_dir+'ISAAC/ISAAC_Ks.txt'
      'HAWKI_Ks':      trans_curve = filters_dir+'HAWKI/HAWKI_Ks.txt'
      'W1':            trans_curve = filters_dir+'WISE/RSR-W1.EE.txt'
      'W2':            trans_curve = filters_dir+'WISE/RSR-W2.EE.txt'
      'W3':            trans_curve = filters_dir+'WISE/RSR-W3.EE.txt'
      'W4':            trans_curve = filters_dir+'WISE/RSR-W4.EE.txt'
      'SPITZER_I1':    trans_curve = filters_dir+'SPITZER/IRAC/080924ch1trans_full.txt'
      'SPITZER_I2':    trans_curve = filters_dir+'SPITZER/IRAC/080924ch2trans_full.txt'
      'SPITZER_I3':    trans_curve = filters_dir+'SPITZER/IRAC/080924ch3trans_full.txt'
      'SPITZER_I4':    trans_curve = filters_dir+'SPITZER/IRAC/080924ch4trans_full.txt'
      'SPITZER_M1':    trans_curve = filters_dir+'SPITZER/MIPS/MIPSfilt_ch1.txt'
      'SPITZER_M2':    trans_curve = filters_dir+'SPITZER/MIPS/MIPSfilt_ch2.txt'
      'SPITZER_M3':    trans_curve = filters_dir+'SPITZER/MIPS/MIPSfilt_ch3.txt'
      'HAWC+A':        trans_curve = filters_dir+'SOFIA/HAWC+/bandA.txt'
      'HAWC+B':        trans_curve = filters_dir+'SOFIA/HAWC+/bandB.txt'
      'HAWC+C':        trans_curve = filters_dir+'SOFIA/HAWC+/bandC.txt'
      'HAWC+D':        trans_curve = filters_dir+'SOFIA/HAWC+/bandD.txt'
      'HAWC+E':        trans_curve = filters_dir+'SOFIA/HAWC+/bandE.txt'
      'HERSCHEL_P1':   trans_curve = filters_dir+'HERSCHEL/PACS/PACS1_EffectiveResponse_blue.txt'
      'HERSCHEL_P2':   trans_curve = filters_dir+'HERSCHEL/PACS/PACS2_EffectiveResponse_green.txt'
      'HERSCHEL_P3':   trans_curve = filters_dir+'HERSCHEL/PACS/PACS3_EffectiveResponse_red.txt'
      'HERSCHEL_S1':   trans_curve = filters_dir+'HERSCHEL/SPIRE/SVO_SPIRE250_ext.txt'
     ;'HERSCHEL_S1':   trans_curve = filters_dir+'HERSCHEL/SPIRE/SVO_SPIRE250.txt'
     ;'HERSCHEL_S1':   trans_curve = filters_dir+'HERSCHEL/SPIRE/HCSS_SPIRE250.txt'
      'HERSCHEL_S2':   trans_curve = filters_dir+'HERSCHEL/SPIRE/SVO_SPIRE350_ext.txt'
     ;'HERSCHEL_S2':   trans_curve = filters_dir+'HERSCHEL/SPIRE/SVO_SPIRE350.txt'
     ;'HERSCHEL_S2':   trans_curve = filters_dir+'HERSCHEL/SPIRE/HCSS_SPIRE350.txt'
      'HERSCHEL_S3':   trans_curve = filters_dir+'HERSCHEL/SPIRE/SVO_SPIRE500_ext.txt'
     ;'HERSCHEL_S3':   trans_curve = filters_dir+'HERSCHEL/SPIRE/SVO_SPIRE500.txt'
     ;'HERSCHEL_S3':   trans_curve = filters_dir+'HERSCHEL/SPIRE/HCSS_SPIRE500.txt'
    endcase
    readcol,trans_curve,wraw,Kraw,form='D,D',/silent              ;read wraw in microns and normalized trasmission
    wcomb = [wraw,wave]
      wcomb = wcomb[sort(wcomb)]
      wcomb = wcomb[unique(wcomb)]
    nu_wcomb = 1.d4*!cv.clight/wcomb
    offwraw       = where((wave  lt min(wraw)) or (wave  gt max(wraw)),/null,comp=inwraw)
    wcomb_offwraw = where((wcomb lt min(wraw)) or (wcomb gt max(wraw)),/null,comp=wcomb_inwraw)
    Kcomb = dblarr(wcomb.length)
      Kcomb[wcomb_inwraw] = interpol(Kraw,wraw,wcomb[wcomb_inwraw])
      Kcomb /= trap_int(wcomb,Kcomb)
    possible = ((min(wraw) ge min(wave)) and (max(wraw) le max(wave))) ; convolution is only possible if wraw is fully included in wave
    if possible then begin
      if (n_elements(inwraw) ne 0) then begin
        wcomb_inwave = []
          for kk=0,(wave.length-1) do wcomb_inwave = [wcomb_inwave,where(wcomb eq wave[kk])]
        Kwave = Kcomb[wcomb_inwave]
        if max(Kwave) ne 0 then begin
          Kwave /= trap_int(wave,Kwave)  ;normalizing filters in the common wavelength grid
                                         ;this normalization is just for plotting purposes.
                                         ;Denominator cancels out in all later calculations
          Filters[f,*] = (Kwave) 
        endif else begin
          Filters[f,*] = !values.D_NaN ;if all elements of Filter are zero
        endelse
      endif else begin
        ;if wave grid is too coarse, the full range of wraw may be between two consecutive points of wave.
        ;this happens, for instance with PEGASE2 models convolved with FIR filters.
        Filters[f,*] = !values.D_NaN
      endelse
    endif else begin
      Filters[f,*] = !values.D_NaN
    endelse
  endif
  
  if analysis then begin
    for t=0L,(Ntime-1) do begin
      Lnu_wcomb = interpol(reform(Lnu[t,*]),wave,wcomb)
      filter_Lnu_wcomb = Kcomb*Lnu_wcomb
      fnorm1 = trap_int(nu_wcomb,Kcomb)                                        ;fnorm1<0
      mean_Lnu[t,f] = trap_int(nu_wcomb,filter_Lnu_wcomb)/temporary(fnorm1)    ;both denom and numerators are negative
      fnorm2 = trap_int(nu_wcomb,filter_Lnu_wcomb)                             ;fnorm1<0
      mean_nu[t,f]  = trap_int(nu_wcomb,nu_wcomb  *filter_Lnu_wcomb)/          fnorm2
      mean_nu2[t,f] = trap_int(nu_wcomb,nu_wcomb^2*filter_Lnu_wcomb)/temporary(fnorm2)
      Lwave_wcomb = 1.d4*!cv.clight*Lnu_wcomb/wcomb^2
      filter_Lwave_wcomb = Kcomb*Lwave_wcomb
      fnorm3 = trap_int(wcomb,Kcomb)
      mean_Lwave[t,f]  = trap_int(wcomb,filter_Lwave_wcomb)/temporary(fnorm3)
      fnorm4 = trap_int(wcomb,filter_Lwave_wcomb)
      mean_wave [t,f] = trap_int(wcomb,wcomb  *filter_Lwave_wcomb)/          fnorm4
      mean_wave2[t,f] = trap_int(wcomb,wcomb^2*filter_Lwave_wcomb)/temporary(fnorm4)
    endfor
  endif

endfor

sigma_nu   = sqrt(mean_nu2 - mean_nu^2)
sigma_wave = sqrt(mean_wave2 - mean_wave^2)

if plot_filters then begin
  !p.charsize=2
  ;Plotting normalized Filters
  plot_oi,wave,wave*Filters[00,*],xr=minmax(wave),yr=[0,max((replicate(1.0,Nfilters)#wave)*Filters)]
    for i=0,(Nfilters-1) do oplot,wave,wave*Filters[i,*],color=54+200*i/(Nfilters-1)
    oplot,wave,wave*Filters[00,*]*0
endif

if plot_analysis then begin
  pmulti_o = !p.multi  
  window,1
  !p.multi=[0,5,3]
  for i=0,(Nfilters-1) do $
    plot_oi,time,(mean_wave*mean_nu)[*,i],psym=-4,xrange=[1.d6,1.d10],ystyle=1,title='!6<!7k!6><!7m!6> ('+Filter_labels[i]+')'
  window,2
  !p.multi=[0,5,3]
  for i=0,(Nfilters-1) do $
    plot_oi,time,mean_nu[*,i],psym=-4,xrange=[1.d6,1.d10],ystyle=1,title='!6<!7m!6> ('+Filter_labels[i]+')'
  window,3
  !p.multi=[0,5,3]
  for i=0,(Nfilters-1) do $
    plot_oi,time,mean_wave[*,i],psym=-4,xrange=[1.d6,1.d10],ystyle=1,title='!6<!7k!6> ('+Filter_labels[i]+')'
  window,4
  !p.multi=[0,5,3]
  for i=0,(Nfilters-1) do $
    plot_oi,time,sigma_wave[*,i],psym=-4,xrange=[1.d6,1.d10],ystyle=1,title='!6!7r!6!d<!7k!6>!n ('+Filter_labels[i]+')'
  window,5
  !p.multi=[0,5,3]
  for i=0,(Nfilters-1) do $
    plot_oi,time,(sigma_wave/mean_wave)[*,i],psym=-4,xrange=[1.d6,1.d10],ystyle=1,$
            title='!6!7r!6!d<!7k!6>!n/<!7k!6> ('+Filter_labels[i]+')'

  !p.multi = pmulti_o
endif

return

end
