function dl07_models, filter_labels=filter_labels, redshift=redshift, error_check=error_check
;+
; Name
; ----
;   DL07_MODELS
;
; Purpose
; -------
;   Generates the spectra, SEDs, and dust parameters for the Draine & Li (2007) 
;   dust models given a set of filters. The spectra and SEDs are tabulated as a 
;   function of dust model parameters.
;
; Calling Sequence
; ----------------
;   ::
;
;       dl07 = dl07_models([filter_labels = , redshift = , /error_check])
;
; Optional Inputs
; ---------------
;   ``filter_labels`` : string array(Nfilters)
;       The filters labels for which models should be generated.
;       (Default = ``['GALEX_FUV', 'GALEX_NUV', 'SDSS_u', 'SDSS_g', 'SDSS_r',
;       'SDSS_i', 'SDSS_z', '2MASS_J', '2MASS_H', '2MASS_Ks', 'IRAC_CH1',
;       'IRAC_CH2', 'IRAC_CH3', 'IRAC_CH4', 'MIPS_CH1', 'PACS_green',
;       'PACS_red', 'SPIRE_250', 'SPIRE_350', 'SPIRE_500']``)
;   ``redshift`` : int, float, or double scalar
;       The redshift of the model. (Default = ``0.0``)
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``dl07`` : structure
;       A structure including the full-resolution DL07 models grid (in terms
;       of Lsun Hz-1) and those convolved with the filters specified
;       by ``filter_labels``.
;       The full description of the structure is as follows:
;
;       =============     ===========================     ==================================================================================
;       TAG               TYPE                            DESCRIPTION
;       =============     ===========================     ==================================================================================
;       FILTER_LABELS     string(Nfilters)                Filters labels, same as input
;       WAVE_FILTERS      double(Nfilters)                Mean wavelength of the filters :math:`[\mu \rm m]`
;       MEAN_LNU          double(Nfilters, NU, Nqpah)     Model spectra convolved with the filters :math:`[L_\odot\ {\rm Hz}^{-1}]`
;       WAVE_REST         double(Nwave)                   Rest-frame wavelength of the high-res models :math:`[\mu \rm m]`
;       WAVE_OBS          double(Nwave)                   Observed-frame wavelength of the high-res models :math:`[\mu \rm m]`
;       LNU               double(Nwave, NU, Nqpah)        Observed-frame Lnu spectrum emitted by each model :math:`[L_\odot\ {\rm Hz}^{-1}]`
;       U                 string(NU)                      Tabulated radiation field intensities
;       MODEL             string(Nqpah)                   Milky Way models corresponding to each ``QPAH``
;       QPAH              double(Nqpah)                   Tabulated PAH index values
;       LBOL              double(NU, Nqpah)               Bolometric luminosity emitted by each model :math:`[L_\odot]`
;       FILTERS           double(Nfilters, Nwave)         Tabulated filter transmission functions, frequency normalized
;       REDSHIFT          float/double                    Redshift of model, same as input
;       =============     ===========================     ==================================================================================
;
; Reference
; ---------
;   `Draine, B. T., & Li, A. 2007, ApJ, 657, 810 <https://ui.adsabs.harvard.edu/abs/2007ApJ...657..810D/abstract>`_
;
; Modification History
; --------------------
;   - 2020/03/01: Created (Rafael T. Eufrasio)
;   - 2020/04/27: Made folder location ``lightning_dir`` vs fix location in code (Keith Doore)
;   - 2020/04/27: Changed ``folder`` keyword to ``lightning_dir`` for consistence between functions/procedures for ``_REF_EXTRA`` (Keith Doore)
;   - 2020/05/11: Corrected (1+z) term (Keith Doore)
;   - 2022/03/15: Added proper error handling (Keith Doore)
;   - 2022/03/15: Renamed variables to standard format (Keith Doore)
;   - 2022/03/18: Updated documentation (Keith Doore)
;   - 2022/04/12: Allowed for ``filter_labels`` to be scalars (Keith Doore)
;   - 2022/04/12: Changed ``readcol`` to skip extra beginning lines (Keith Doore)
;   - 2022/04/12: Replaced ``.length`` with ``n_elements`` (Keith Doore)
;   - 2022/05/16: Turned ``lightning_dir`` into system variable call (Keith Doore)
;   - 2022/05/16: Changed ``!cv`` constants call to ``!lightning_cgs`` call (Keith Doore)
;   - 2022/05/16: Removed ``get_filters`` input of ``lightning_dir``, since folded into ``get_filters`` (Keith Doore)
;   - 2022/05/17: Added ``error_check`` keyword to do error handling (Keith Doore)
;   - 2022/05/17: Changed function name from ``dl07_templates`` to ``dl07_models`` (Keith Doore)
;   - 2022/05/18: Updated path to model files (Keith Doore)
;   - 2022/06/29: Converted ``U`` in output to double from string (Keith Doore)
;   - 2022/06/29: Replaced for loops with ``reform`` for ``Lnu`` and ``mean_Lnu`` arrays (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(filter_labels) ne 0 then begin
     if size(filter_labels, /type) ne 7 then message, 'FILTER_LABELS must be of type string.'
     if size(reform(filter_labels), /n_dim) ne 1 then message, 'FILTER_LABELS must be a scalar or 1-D array.'
   endif
  
   if n_elements(redshift) ne 0 then begin
     if size(redshift, /type) lt 2 or size(redshift,/type) gt 5 then $
              message, 'REDSHIFT must be of type int, float, or double.'
     if size(redshift, /n_dim) ne 0 then message, 'REDSHIFT must be a scalar.'
     if redshift lt 0 then message, 'REDSHIFT must only contain non-negative values.'
   endif
 endif

 if n_elements(filter_labels) eq 0 then $
   filter_labels=['GALEX_FUV', 'GALEX_NUV', 'SDSS_u', 'SDSS_g', 'SDSS_r', 'SDSS_i', 'SDSS_z', $
                  '2MASS_J', '2MASS_H', '2MASS_Ks', 'IRAC_CH1', 'IRAC_CH2', 'IRAC_CH3', $
                  'IRAC_CH4', 'MIPS_CH1', 'PACS_green', 'PACS_red', 'SPIRE_250', 'SPIRE_350', 'SPIRE_500']
 if n_elements(redshift) eq 0 then redshift = 0.0


; Generate dust emission models
 Ugrid=['0.10','0.15','0.20','0.30','0.40','0.50','0.70','0.80','1.00',$
        '1.20','1.50','2.50','3.00','4.00','5.00','7.00','8.00','12.0',$
        '15.0','20.0','25.0','1e2','3e2','1e3','3e3','1e4','3e4','1e5','3e5']

 mod_grid = []                      &  qPAH_grid = []
 mod_grid = [mod_grid, 'MW3.1_00']  &  qPAH_grid = [qPAH_grid, 0.47d-2]  ;qPAH = 0.47 %
 mod_grid = [mod_grid, 'MW3.1_10']  &  qPAH_grid = [qPAH_grid, 1.12d-2]  ;qPAH = 1.12 %
 mod_grid = [mod_grid, 'MW3.1_20']  &  qPAH_grid = [qPAH_grid, 1.77d-2]  ;qPAH = 1.77 %
 mod_grid = [mod_grid, 'MW3.1_30']  &  qPAH_grid = [qPAH_grid, 2.50d-2]  ;qPAH = 2.50 %
 mod_grid = [mod_grid, 'MW3.1_40']  &  qPAH_grid = [qPAH_grid, 3.19d-2]  ;qPAH = 3.19 %
 mod_grid = [mod_grid, 'MW3.1_50']  &  qPAH_grid = [qPAH_grid, 3.90d-2]  ;qPAH = 3.90 %
 mod_grid = [mod_grid, 'MW3.1_60']  &  qPAH_grid = [qPAH_grid, 4.58d-2]  ;qPAH = 4.58 %
 ;mod_grid = [mod_grid,'LMC2_00']   &  qPAH_grid = [qPAH_grid,0.75d-2]   ;qPAH = 0.75 %
 ;mod_grid = [mod_grid,'LMC2_05']   &  qPAH_grid = [qPAH_grid,1.49d-2]   ;qPAH = 1.49 %
 ;mod_grid = [mod_grid,'LMC2_10']   &  qPAH_grid = [qPAH_grid,2.37d-2]   ;qPAH = 2.37 %
 ;mod_grid = [mod_grid,'smc']       &  qPAH_grid = [qPAH_grid,0.10d-2]   ;qPAH = 0.10 %

 Nfilters = n_elements(filter_labels)
 n_U = n_elements(Ugrid)
 n_mod = n_elements(mod_grid)

 ; => lambdas are all exactly the same for all delta functions
 Lnu  = dblarr(1001, n_U, n_mod)
 Lbol = dblarr(n_U, n_mod)
 for uu=0,(n_U-1) do begin
   U = Ugrid[uu]
   for mm=0,(n_mod-1) do begin
     modl = mod_grid[mm]
     readcol, !lightning_dir+'models/dust/DL07/DL07spec/U'+U+'/U'+U+'_'+U+'_'+modl+'.txt',$
              wave_uumm, nuLnu_uumm, jnu_min_uumm, SKIPLINE=61, /silent
     wave = reverse(wave_uumm) ; wave is increasing, in microns
     nu = 1.d4* !lightning_cgs.clight/wave
     Lnu[*, uu, mm] = reverse(nuLnu_uumm)/nu
     ; Integral is negative since nu is decreasing (big to small)
     Lbol[uu, mm] = abs(trap_int(nu, Lnu[*, uu, mm]))
   endfor
 endfor

 nu_rest   = temporary(nu)              ; restframe frequency
 nu_obs    = nu_rest/(1.d0+redshift)    ; observed frequency
 wave_rest = temporary(wave)            ; restframe wavelength
 wave_obs  = wave_rest*(1.d0+redshift)  ; observed wavelength
 Lnu_rest  = temporary(Lnu)             ; restframe Lnu spectrum
 Lnu       = Lnu_rest * (1.d0+redshift) ; observed Lnu spectrum, where Lnu = 4*!dpi*DL^2*Fnu

 get_filters, filter_labels, nu_obs, filters, Lnu=[1.d0] # [0*(nu_obs) + 3631], $ ; flat AB=0 source
              mean_wave=mean_wave

 wave_filters = reform(mean_wave)

 Lnu_reformed = reform(Lnu, n_elements(nu_obs), n_u*n_mod)
 get_filters, filter_labels, nu_obs, filters, Lnu=transpose(Lnu_reformed),$
              mean_Lnu=mean_Lnu_reformed

 mean_Lnu = reform(transpose(mean_Lnu_reformed), Nfilters, n_u, n_mod)
 ;filters not fully included in DL07 wave range produce mean_Lnu = NaN
 mean_Lnu[where(finite(mean_Lnu,/naN), /null)] = 0.0

 dl07 = {filter_labels : filter_labels, $
         wave_filters  : wave_filters,  $
         mean_Lnu      : mean_Lnu,      $
         wave_rest     : wave_rest,     $
         wave_obs      : wave_obs,      $
         Lnu           : Lnu,           $
         U             : double(Ugrid), $
         model         : mod_grid,      $
         qPAH          : qPAH_grid,     $
         Lbol          : Lbol,          $
         Filters       : Filters,       $
         redshift      : redshift}

 return, dl07

end
