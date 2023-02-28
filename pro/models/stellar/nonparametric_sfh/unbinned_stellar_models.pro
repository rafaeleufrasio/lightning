function unbinned_stellar_models, ages, filter_labels=filter_labels, redshift=redshift, $
                                dtime_SF=dtime_SF, Zmetal=Zmetal, IMF=IMF, $
                                cosmology=cosmology, no_emission_lines=no_emission_lines, $
                                no_nebular_extinction=no_nebular_extinction, error_check=error_check
;+
; Name
; ----
;   BINNED_STELLAR_MODELS
;
; Purpose
; -------
;   Generates the spectra, SEDs, and stellar parameters for given a set of
;   filters and steps boundaries. The spectra and SEDs can include or not
;   include nebular absorption and emission lines, and they depend on the
;   optionally input metallicity and initial mass function (IMF).
;
; Calling Sequence
; ----------------
;   ::
;
;       stellar_models = binned_stellar_models([filter_labels = , redshift = , $
;                               steps_bounds = , dtime_SF = , Zmetal = , IMF = , $
;                               cosmology = , /no_emission_lines, /no_nebular_extinction, $
;                               /error_check])
;
; Optional Inputs
; ---------------
;   ``ages`` : int, float, or double array(Nages)
;       The ages at which to evaluate the populations :math:`[{\rm yr}]`.
;   ``filter_labels`` : string array(Nfilters)
;       The filters labels for which models should be generated.
;       (Default = ``['GALEX_FUV', 'GALEX_NUV', 'SDSS_u', 'SDSS_g', 'SDSS_r',
;       'SDSS_i', 'SDSS_z', '2MASS_J', '2MASS_H', '2MASS_Ks', 'IRAC_CH1',
;       'IRAC_CH2', 'IRAC_CH3', 'IRAC_CH4', 'MIPS_CH1', 'PACS_green',
;       'PACS_red', 'SPIRE_250', 'SPIRE_350', 'SPIRE_500']``)
;   ``redshift`` : int, float, or double scalar
;       The redshift of the model. (Default = ``0.0``)
;   ``dtime_SF`` : int, float, or double scalar
;       The time step used when interpolating the stellar population to the age bins
;       defined by ``steps_bounds`` :math:`[{\rm yr}]`. (Default = ``5.d5``)
;   ``Zmetal`` : float or double scalar
;       The metallicity in terms of Z to be used when generating the stellar models.
;       (Default = ``0.02``)
;   ``IMF`` : string scalar
;       The IMF to be used when generating the stellar models. The only available
;       value is ``'Kroupa01'``. (Default = ``'Kroupa01'``)
;   ``cosmology`` : structure
;       A structure containing the cosmology parameters ``H0``, ``OMEGA_M``, ``LAMBDA0``,
;       ``Q0``, and ``K`` in individual tags.
;       (Default = ``{H0: 70, OMEGA_M: 0.3, LAMBDA0: 0.7, Q0: -0.55, K: 0}``)
;   ``no_emission_lines`` : flag
;       If set, no emission lines are included in the output spectra.
;   ``no_nebular_extinction`` : flag
;       If set, no nebular absorption is included in the output spectra.
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``stellar_models`` : structure
;       A structure including the full-resolution non-parametric stellar models
;       (in terms of :math:`L_\odot\ {\rm Hz}^{-1}`) and those convolved with the filters
;       specified by ``filter_labels``.
;       The full description of the structure is as follows:
;
;       =============     ========================     ========================================================================================
;       TAG               TYPE                         DESCRIPTION
;       =============     ========================     ========================================================================================
;       FILTER_LABELS     string(Nfilters)             Labels for filters, same as input
;       WAVE_FILTERS      double(Nfilters)             Mean wavelength of the filters :math:`[\mu \rm m]`
;       MEAN_LNU          double(Nfilters, Nages)     Model spectra convolved with the filters :math:`[L_\odot\ {\rm Hz}^{-1}]`
;       WAVE_REST         double(Nwave)                Rest-frame wavelength of the high-res models :math:`[\mu \rm m]`
;       WAVE_OBS          double(Nwave)                Observed-frame wavelength of the high-res models :math:`[\mu \rm m]`
;       LNU               double(Nwave, Nages)        Observed-frame Lnu spectrum emitted within each age bin :math:`[L_\odot\ {\rm Hz}^{-1}]`
;       Q0                double(Nages)               Ionizing flux rate emitted within each age bin :math:`[{\rm photons\ s^{-1}}]`
;       MSTAR             double(Nages)               Surviving stellar mass within each age bin :math:`[M_\odot]`
;       LBOL              double(Nages)               Bolometric luminosity emitted within each age bin :math:`[L_\odot]`
;       AGES              double(Nages+1)             Age bounds to separate each star formation step :math:`[{\rm yr}]`
;       FILTERS           double(Nfilters, Nwave)      Tabulated filter transmission functions, frequency normalized
;       REDSHIFT          float/double                 Redshift of model, same as input
;       =============     ========================     ========================================================================================
;
; Modification History
; --------------------
;   - 2016/05/01: Created (Rafael T. Eufrasio)
;   - 2020/03/09: Corrected (1+z) factor (Rafael T. Eufrasio)
;   - 2020/04/27: Replaced if statements using ``n_elements`` on keywords to use ``keyword_set`` (Keith Doore)
;   - 2022/03/15: Added proper error handling (Keith Doore)
;   - 2022/03/15: Renamed variables to standard format (Keith Doore)
;   - 2022/03/18: Updated documentation (Keith Doore)
;   - 2022/04/12: Allowed for ``filter_labels`` to be scalars (Keith Doore)
;   - 2022/04/12: Allowed inputs to be integers (Keith Doore)
;   - 2022/04/12: Replaced ``.length`` with ``n_elements`` (Keith Doore)
;   - 2022/04/13: Normalized Filters to the frequency grid (Keith Doore)
;   - 2022/05/16: Turned ``lightning_dir`` into system variable call (Keith Doore)
;   - 2022/05/16: Changed ``!cv`` constants call to ``!lightning_cgs`` call (Keith Doore)
;   - 2022/05/16: Removed ``Filters`` output, since it is redundant (Keith Doore)
;   - 2022/05/16: Changed ``no_lines`` to ``emission_lines`` to match naming scheme and inverted logic accordingly (Keith Doore)
;   - 2022/05/16: Changed ``no_nebular`` to ``nebular_extinction`` to match naming scheme and inverted logic accordingly (Keith Doore)
;   - 2022/05/16: Changed ``Zmet`` to ``Zmetal`` to match naming scheme (Keith Doore)
;   - 2022/05/16: Added bin truncation based on redshift (Keith Doore)
;   - 2022/05/16: Changed model save file from ``.idl`` to ``.sav`` (Keith Doore)
;   - 2022/05/16: Removed ``get_filters`` input of ``lightning_dir``, since folded into ``get_filters`` (Keith Doore)
;   - 2022/05/17: Added ``error_check`` keyword to do error handling (Keith Doore)
;   - 2022/05/17: Changed function name from ``steps_stellar`` to ``binned_stellar_models`` (Keith Doore)
;   - 2022/05/18: Added optional ``cosmology`` input to allow for unique cosmologies in ``galage`` call (Keith Doore)
;   - 2022/05/18: Updated path to model files (Keith Doore)
;   - 2022/06/29: Changed ``emission_lines`` to back to ``no_emission_lines`` as to be default if not set (Keith Doore)
;   - 2022/06/29: Changed ``nebular_extinction`` to back to ``no_nebular_extinction`` as to be default if not set (Keith Doore)
;-
 ;On_error, 2
 compile_opt idl2

 ; Error handling
 if keyword_set(error_check) then begin
   if n_elements(filter_labels) ne 0 then begin
     if size(filter_labels, /type) ne 7 then message, 'FILTER_LABELS must be of type string.'
     if size(reform(filter_labels), /n_dim) ne 1 then message, 'FILTER_LABELS must be a scalar or 1-D array.'
   endif

   if n_elements(ages) ne 0 then begin
     if size(ages, /type) lt 2 or size(ages,/type) gt 5 then $
              message, 'AGES must be of type int, float, or double.'
     if size(ages, /n_dim) ne 1 then message, 'AGES must be a 1-D array.'
     if n_elements(ages) eq 1 then message, 'AGES must have more than one element.'
     if total(ages lt 0) gt 0 then message, 'AGES must only contain non-negative values.'
   endif

   if n_elements(redshift) ne 0 then begin
     if size(redshift, /type) lt 2 or size(redshift,/type) gt 5 then $
              message, 'REDSHIFT must be of type int, float, or double.'
     if size(redshift, /n_dim) ne 0 then message, 'REDSHIFT must be a scalar.'
     if redshift lt 0 then message, 'REDSHIFT must be a non-negative value.'
   endif

   if n_elements(dtime_SF) ne 0 then begin
     if size(dtime_SF, /type) lt 2 or size(dtime_SF,/type) gt 5 then $
              message, 'DTIME_SF must be of type int, float, or double.'
     if size(dtime_SF, /n_dim) ne 0 then message, 'DTIME_SF must be a scalar.'
     if dtime_SF le 0 then message, 'DTIME_SF must be a positive value.'
   endif

   if n_elements(Zmetal) ne 0 then begin
     if size(Zmetal, /type) lt 4 or size(Zmetal,/type) gt 5 then $
              message, 'ZMETAL must be of type float or double.'
     if size(Zmetal, /n_dim) ne 0 then message, 'ZMETAL must be a scalar.'
     if Zmetal le 0 then message, 'ZMETAL must be a non-negative value.'
   endif

   if n_elements(IMF) ne 0 then begin
     if size(IMF, /type) ne 7 then message, 'IMF must be of type string.'
     if size(IMF, /n_dim) ne 0 then message, 'IMF must be a scalar.'
     if IMF ne 'Kroupa01' then message, "Only current IMF option is 'Kroupa01'."
   endif
 endif

 if n_elements(filter_labels) eq 0 then $
   filter_labels=['GALEX_FUV', 'GALEX_NUV', 'SDSS_u', 'SDSS_g', 'SDSS_r', 'SDSS_i', 'SDSS_z', $
                  '2MASS_J', '2MASS_H', '2MASS_Ks', 'IRAC_CH1', 'IRAC_CH2', 'IRAC_CH3', $
                  'IRAC_CH4', 'MIPS_CH1', 'PACS_green', 'PACS_red', 'SPIRE_250', 'SPIRE_350', 'SPIRE_500']
 if n_elements(steps_bounds) eq 0 then steps_bounds = [0.d0, 1.d7, 1.d8, 1.d9, 5.d9, 13.6d9]
 if n_elements(redshift) eq 0 then redshift = 0.0
 if n_elements(dtime_SF) eq 0 then dtime_SF = 5.d5
 if n_elements(Zmetal) eq 0 then Zmetal = 0.02
 if n_elements(IMF) eq 0 then IMF='kroupa01'
 Ztag = string(Zmetal,form='(F5.3)')


; Read in models from file
 tform = galage(redshift, 1.d3, /silent, _extra=cosmology)
 too_old = where(ages gt tform, Ntoo_old)
 if Ntoo_old gt 0 then begin
     print, 'WARNING -- At least one of the supplied AGES is older than the Universe at REDSHIFT. These ages will be set to the maximum allowed age in the output.'
     ages[too_old] = tform
 endif

 too_young = where(ages lt 1.d6, Ntoo_young)
 if Ntoo_young gt 0 then begin
     print, 'WARNING -- At least one of the supplied AGES is younger than 1 Myr. These ages will be set to 1 Myr in the output.'
     ages[too_young] = 1.d6
 endif

 Nages = n_elements(ages)


 if ~keyword_set(no_nebular_extinction) then nebtag='nebular_' else nebtag=''
 bursts_folder = !lightning_dir+'models/stellar/PEGASE/single_burst/'
 filename = bursts_folder+strlowcase(IMF)+'/'+strlowcase(IMF)+'_Z'+Ztag+'_'+nebtag+'spec.sav'

 save_obj = OBJ_NEW('IDL_Savefile', filename)
 save_obj->RESTORE, 'time'
 save_obj->RESTORE, 'wave'
 save_obj->RESTORE, 'nu'
 save_obj->RESTORE, 'Lnu'
 save_obj->RESTORE, 'Lbol'
 save_obj->RESTORE, 'Nlyc'
 save_obj->RESTORE, 'Mstars'
 Ntime = n_elements(time)
 Nwave_burst = n_elements(wave)

 if ~keyword_set(no_emission_lines) then begin
   save_obj->RESTORE, 'L_lines'
   save_obj->RESTORE, 'wlines'
   Nwlines = n_elements(wlines)
   vd = 5d6 ; 50km/s = 5,000,000 cm/s
   z = vd/!lightning_cgs.clight
   for tt=0, (Ntime-1) do begin
     Lnu_lines = wave*0
     for ll=0, (Nwlines-1) do begin
       x2 = ((wave/(wlines[ll])[0] - 1.d0)/z)^2
       Lnu_line = exp(-x2)
       Lnu_line /= abs(trap_int(nu, Lnu_line))
       Lnu_lines += (L_lines[tt, ll])[0]*Lnu_line
     endfor
     Lnu[tt, *] += Lnu_lines
   endfor
 endif

 age_burst   = temporary(time)
 nu_rest     = temporary(nu)                 ; restframe frequency
 nu_obs      = nu_rest / (1.d0 + redshift)   ; observed frequency
 wave_rest   = temporary(wave)               ; restframe wavelength
 wave_obs    = wave_rest * (1.d0 + redshift) ; observed wavelength
 Lnu_burst_rest = temporary(Lnu)             ; restframe Lnu spectrum
 Lnu_burst      = Lnu_burst_rest * (1.d0 + redshift) ; observed Lnu spectrum, where Lnu = 4*!dpi*DL^2*Fnu
 ;Lbol_burst = Lbol/(1.d0+redshift)  & Lbol   = !null
 ;Q0_burst   = Nlyc/(1.d0+redshift)  & Nlyc   = !null
 ;Lnu_burst   = Lnu                 & Lnu    = !null
 Lbol_burst  = Lbol                 & Lbol   = !null ; observed = restframe bolometric luminosity
 Q0_burst    = Nlyc                 & Nlyc   = !null ; restframe ionizing flux rate
 ;Q0_burst    = (1.d0+redshift)*Nlyc & Nlyc   = !null ; observed ionizing flux rate
 Mstar_burst = Mstars               & Mstars = !null
 get_filters, filter_labels, nu_obs, filters, Lnu=[1.d0] # [0*(nu_obs) + 3631],$
              mean_wave=mean_wave

 wave_filters = reform(mean_wave)

 Q0_ages    = dblarr(Nages)
 Lnu_ages   = dblarr(Nwave_burst, Nages)
 Lbol_ages  = dblarr(Nages)
 Mstar_ages = dblarr(Nages)

 for j=0, (Nages-1) do begin

   t = ages[j]
   Q0_ages[j]    = interpol(Q0_burst,    age_burst, t)
   Lbol_ages[j]  = interpol(Lbol_burst,  age_burst, t)
   Mstar_ages[j] = interpol(Mstar_burst, age_burst, t)

   for i=0,(n_elements(nu_obs) - 1) do begin
     Lnu_ages[i, j] = interpol(reform(Lnu_burst[*, i]), age_burst, t)
   endfor

 endfor

 get_filters, filter_labels, nu_obs, filters, Lnu=transpose(Lnu_ages),$
              mean_Lnu=mean_Lnu

 mean_Lnu_ages = transpose(temporary(mean_Lnu))

 stellar_models = {filter_labels : filter_labels,  $
                   wave_filters  : wave_filters,   $
                   mean_Lnu      : mean_Lnu_ages, $
                   wave_rest     : wave_rest,      $
                   wave_obs      : wave_obs,       $
                   Lnu           : Lnu_ages,      $
                   Q0            : Q0_ages,       $
                   Mstar         : Mstar_ages,    $
                   Lbol          : Lbol_ages,     $
                   ages          : ages,   $
                   filters       : filters,        $
                   redshift      : redshift}

 return, stellar_models

end
