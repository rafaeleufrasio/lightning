function skirtor_models, filter_labels=filter_labels, redshift=redshift, error_check=error_check
;+
; Name
; ----
;   SKIRTOR_MODELS
;
; Purpose
; -------
;   Generates the spectra, SEDs, and AGN parameters for the SKIRTOR model 
;   (Stalevski et al. 2016) given a set of filters. The model is a slice 
;   of :math:`\tau`-:math:`\delta`-:math:`i` space, with constant values of
;   :math:`R`, :math:`p`, and :math:`q`.
;
; Calling Sequence
; ----------------
;   ::
;
;       skirtor = skirtor_models([filter_labels = , redshift = , /error_check])
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
;   ``skirtor`` : structure
;       A structure including the full-resolution SKIRTOR models grid (in terms
;       of :math:`L_\odot\ {\rm Hz}^{-1}`) and those convolved with the filters specified
;       by ``filter_labels``.
;       The full description of the structure is as follows:
;
;       ===============     ============================    =================================================================================================
;       TAG                 TYPE                            DESCRIPTION
;       ===============     ============================    =================================================================================================
;       FILTER_LABELS       string(Nfilters)                Filters labels, same as input
;       WAVE_FILTERS        double(Nfilters)                Mean wavelength of the filters :math:`[\mu \rm m]`
;       MEAN_LNU            double(Nfilters, Ninc, Ntau)    Model spectra convolved with the filters :math:`[L_\odot\ {\rm Hz}^{-1}]`
;       WAVE_REST           double(Nwave)                   Rest-frame wavelength of the high-res models :math:`[\mu \rm m]`
;       WAVE_OBS            double(Nwave)                   Observed-frame wavelength of the high-res models :math:`[\mu \rm m]`
;       LNU_TOTAL           double(Nwave, Ninc, Ntau)       Total model spectra after reprocessing of disk emission by torus :math:`[L_\odot\ {\rm Hz}^{-1}]`
;       LNU_TRANSPARENT     double(Nwave, Ninc, Ntau)       Model accretion disk spectra prior to reprocessing :math:`[L_\odot\ {\rm Hz}^{-1}]`
;       LNU_INTEGRATED      double(Nwave, Ninc, Ntau)       Total inclination integrated model spectra :math:`[L_\odot\ {\rm Hz}^{-1}]`
;       DUST_MASS_TOTAL     double(Ninc, Ntau)              Total dust mass of each model :math:`[M_\odot]`
;       L2500A_ZERO         double                          Intrinsic accretion disk L2500 at ``i=0`` :math:`[L_\odot\ {\rm Hz}^{-1}]`
;       FILTERS             double(Nfilters, Nwave)         Tabulated filter transmission functions, frequency normalized
;       REDSHIFT            double                          Redshift of model, same as input
;       ===============     ============================    =================================================================================================
;
; Reference
; ---------
;   `Stalevski, M., Ricci, C., Ueda, Y., et al. 2016, MNRAS, 458, 2288 <https://ui.adsabs.harvard.edu/abs/2016MNRAS.458.2288S/abstract>`_
;
; Modification History
; --------------------
;   - 2021/03/03: Created (Erik B. Monson)
;   - 2022/03/15: Documentation update (Erik B. Monson)
;   - 2022/05/17: Turned ``lightning_dir`` into system variable call (Keith Doore)
;   - 2022/05/17: Removed ``get_filters`` input of ``lightning_dir``, since folded into ``get_filters`` (Keith Doore)
;   - 2022/05/17: Added proper error handling (Keith Doore)
;   - 2022/05/17: Added ``error_check`` keyword to do error handling (Keith Doore)
;   - 2022/05/17: Changed function name from ``load_SKIRTOR_grid`` to ``skirtor_models`` (Keith Doore)
;   - 2022/05/18: Updated path to model files (Keith Doore)
;   - 2022/05/18: Changed model save file from ``.idl`` to ``.sav`` (Keith Doore)
;   - 2022/06/01: Updated to fix covering factor and changed read in file format (Erik B. Monson)
;   - 2022/06/09: Merged updates between versions with and without ``CF`` (Keith Doore)
;   - 2022/06/30: Updated documentation (Keith Doore)
;   - 2022/06/30: Updated default values (Keith Doore)
;   - 2022/08/01: Added ``/silent`` to ``mrdfits`` (Keith Doore)
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
     if redshift lt 0 then message, 'REDSHIFT must be a non-negative value.'
   endif
 endif

 if n_elements(filter_labels) eq 0 then $
   filter_labels=['GALEX_FUV', 'GALEX_NUV', 'SDSS_u', 'SDSS_g', 'SDSS_r', 'SDSS_i', 'SDSS_z', $
                  '2MASS_J', '2MASS_H', '2MASS_Ks', 'IRAC_CH1', 'IRAC_CH2', 'IRAC_CH3', $
                  'IRAC_CH4', 'MIPS_CH1', 'PACS_green', 'PACS_red', 'SPIRE_250', 'SPIRE_350', 'SPIRE_500']
 if n_elements(redshift) eq 0 then redshift = 0.0


; Define the grids for these three parameters
 tau_grid = [3, 5, 7, 9, 11]
 ;delta_grid = [10, 20, 30, 40, 50, 60, 70, 80]
 i_grid = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
 Ninc = n_elements(i_grid)
 Ntau = n_elements(tau_grid)
 Nfilters = n_elements(filter_labels)
 Nwave = 132

  ; constant parameters: p = 1, q = 0, R = 20
 p_str = '1'
 q_str = '0'
 R_str = '20'


 SKIRTOR_models = mrdfits(!lightning_dir + 'models/agn/stalevski2016/SKIRTOR/SKIRTOR.fits.gz', 1, /silent)
 wave_rest = (SKIRTOR_models.WAVE_REST)[*,0] ; All the models exist on the same grid.
 nu_rest = (SKIRTOR_models.NU_REST)[*,0]
 wave_obs = (1 + redshift) * wave_rest
 nu_obs = nu_rest / (1 + redshift)
 Lnu_total_grid = SKIRTOR_models.LNU_TOTAL * (1 + redshift)
 Lnu_transparent_grid = SKIRTOR_models.LNU_TRANSPARENT * (1 + redshift)
 Lnu_integrated_grid = SKIRTOR_models.LNU_INTEGRATED * (1 + redshift)
 dust_mass = SKIRTOR_models.DUST_MASS_TOTAL

 ; Observe the total spectrum for all models at once
 ; Note that GET_FILTERS wants Lnu to have dimensions [Nmodels, Nwave]
 ; so we need to pass it the transpose to get a [Nmodels, Nwave] shaped array.
 get_filters, filter_labels, nu_obs, Filters, Lnu=transpose(Lnu_total_grid),$
              mean_Lnu=mean_Lnu, mean_wave=mean_wave

 mean_Lnu = transpose(reform(mean_Lnu, Ntau, Ninc, Nfilters), [2,1,0])
 Lnu_total_grid = transpose(reform(transpose(Lnu_total_grid), Ntau, Ninc, Nwave), [2,1,0])
 Lnu_transparent_grid = transpose(reform(transpose(Lnu_transparent_grid), Ntau, Ninc, Nwave), [2,1,0])
 Lnu_integrated_grid = transpose(reform(transpose(Lnu_integrated_grid), Ntau, Ninc, Nwave), [2,1,0])
 dust_mass = transpose(reform(transpose(dust_mass), Ntau, Ninc))

 ; Integrate in a box filter to get L2500A (Lnu at rest frame 2500 Angstroms)
 ; for the i=0 viewing angle (L_2500A_zero = L2500A(i = 0)). We'll use this to get
 ; the intrinsic L2500A for an arbitrary viewing angle later on.
 ; Note that the transparent flux is independent of all model parameters except viewing angle
 ; which it depends on as a quadratic in cos(i).
 lam0 = 0.25 ; 2500 Angstroms = 0.25 microns
 hw = 100.d-4 ; Filter half width is 100 Angstroms = 0.01 microns

 filter = 0 * wave_rest
 mask = (wave_rest ge (lam0 - hw)) and (wave_rest le (lam0 + hw))

 filter[where(mask, /NULL)] = 1

 filter = filter / trap_int(wave_rest, filter)

 L_2500A_zero = trap_int(wave_rest, reform(Lnu_transparent_grid[*,0,0]) * filter)

 skirtor = {FILTER_LABELS: filter_labels,          $
            WAVE_FILTERS: reform(mean_wave[0,*]),  $
            MEAN_LNU: mean_Lnu,                    $
            WAVE_REST: wave_rest,                  $
            WAVE_OBS: wave_obs,                    $
            LNU_TOTAL: Lnu_total_grid,             $
            LNU_TRANSPARENT: Lnu_transparent_grid, $
            LNU_INTEGRATED: Lnu_integrated_grid,   $
            DUST_MASS_TOTAL: dust_mass,            $
            L2500A_ZERO: L_2500A_zero,             $
            FILTERS: Filters,                      $
            REDSHIFT: redshift}

 return, skirtor

end