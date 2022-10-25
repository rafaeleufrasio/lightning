function lightning_models, config, filter_labels=filter_labels, redshift=redshift, lumin_dist=lumin_dist, $
                           xray_bandpass=xray_bandpass, xray_exposure=xray_exposure, arf_E_lo=arf_E_lo, $
                           arf_E_hi=arf_E_hi, arf_specresp=arf_specresp, galactic_nh=galactic_nh, $
                           input_dir=input_dir, error_check=error_check
;+
; Name
; ----
;   LIGHTNING_MODELS
;
; Purpose
; -------
;   Calls the functions to generate each model structure (stellar, dust, AGN,
;   and/or X-ray) specified in the Lightning configuration structure. The
;   structures are then placed in a parent structure for ease in passing to
;   other Lightning functions.
;
; Calling Sequence
; ----------------
;   ::
;
;       models = lightning_models(config [, filter_labels = , redshift = , lumin_dist = , $
;                                 xray_bandpass = , xray_exposure = , arf_E_lo = , arf_E_hi = , $
;                                 arf_specresp = , galactic_nh = , input_dir = , /error_check])
;
; Input
; -----
;   ``config`` : structure
;       A Lightning configuration structure. (See
;       ``lightning_configure_defaults.pro`` for details and contents.)
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
;   ``lumin_dist`` : int, float, double scalar
;       The luminosity distance of the model :math:`[\rm Mpc]`. (Default = ``10``)
;   ``xray_bandpass`` : int, float or double array(2, Nxray)
;       The bandpasses for the X-ray spectrum. The first column should be the lower
;       energy bound, and the second should be the upper energy bound :math:`[\rm keV]`.
;   ``xray_exposure`` : int, float or double array(Nxray)
;       The exposure time of the observations, one per band :math:`[\rm s]`.
;   ``arf_E_lo`` : float or double array(Nchannels)
;       Lower energy bounds of each channel in the ARF :math:`[\rm keV]`.
;   ``arf_E_hi`` : float or double array(Nchannels)
;       Upper energy bounds of each channel in the ARF :math:`[\rm keV]`.
;   ``arf_specresp`` : float or double array(Nchannels)
;       The spectral response of the ARF at each channel :math:`[\rm cm^2]`.
;   ``galactic_nH`` : int, float, or double scalar
;       Galactic, i.e. Milky Way, neutral Hydrogen column density along the line of
;       sight :math:`[10^{20}\ \rm{cm}^{-2}]`.
;   ``input_dir`` : string scalar
;       The path to the file containing the input SED data.
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``models`` : structure
;       A structure containing the filter labels, and each model structure
;       (stellar, dust, AGN, and/or X-ray) as a substructures. If a model
;       component is not specified in the Lighting configuration structure,
;       then the value of the substructure will be set to ``NaN``.
;       The full description of the structure is as follows:
;
;       ==============     ================     =====================================================================================
;       TAG                TYPE                 DESCRIPTION
;       ==============     ================     =====================================================================================
;       FILTER_LABELS      string(Nfilters)     Labels for filters, same as input
;       STELLAR_MODELS     structure            Stellar emission models (See ``binned_stellar_models.pro`` for details and contents.)
;       ATTEN_MODELS       structure            The preloaded files for the Doore+21 attenuation.
;       DUST_MODELS        structure            Dust emission models (See ``dl07_models.pro`` for details and contents.)
;       XRAY_MODELS        structure            X-ray emission models (See ``xrb_xagn_models.pro`` for details and contents.)
;       AGN_MODELS         structure            AGN emission models (See ``skirtor_models.pro`` for details and contents.)
;       ==============     ================     =====================================================================================
;
; Notes
; -----
;   - When using an X-ray emission model with ``XRAY_UNIT='COUNTS'``, the optional inputs ``xray_bandpass``,
;     ``xray_exposure``, ``arf_E_lo``, ``arf_E_hi``, and ``arf_specresp`` become
;     required inputs.
;   - When using an X-ray emission model with ``XRAY_UNIT='FLUX'``, the optional inputs ``xray_bandpass``,
;     becomes a required input.
;   - When using the ``DOORE21`` attenuation curves, the optional input ``input_dir``
;     becomes a required input.
;
; Modification History
; --------------------
;   - 2022/07/11: Created (Keith Doore)
;   - 2022/08/01: Added ``/silent`` to ``mrdfits`` (Keith Doore)
;   - 2022/08/09: Added ``GALACTIC_NH`` as input rather than keyword inheritance from 
;     ``config``, since now in input file (Keith Doore)
;   - 2022/10/25: Renamed SPS to SSP (Keith Doore)
;-
 ;On_error, 2
 compile_opt idl2

 ; Error handling
 if keyword_set(error_check) then begin
   if n_elements(config) eq 0 then message, 'Variable is undefined: CONFIG.'
   if size(config, /type) ne 8 then message, 'CONFIG is not of type structure.'

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

   if n_elements(lumin_dist) ne 0 then begin
     if size(lumin_dist, /type) lt 2 or size(lumin_dist,/type) gt 5 then $
              message, 'LUMIN_DIST must be of type int, float, or double.'
     if size(lumin_dist, /n_dim) ne 0 then message, 'LUMIN_DIST must be a scalar.'
     if lumin_dist le 0 then message, 'LUMIN_DIST must be a positive value.'
   endif

   if n_elements(xray_bandpass) ne 0 then begin
     if size(xray_bandpass, /type) lt 2 or size(xray_bandpass, /type) gt 5 then $
       message, 'XRAY_BANDPASS must be of type int, float, or double.'
     if size(reform(xray_bandpass), /n_dim) lt 1 or size(reform(xray_bandpass), /n_dim) gt 2 then $
       message, 'XRAY_BANDPASS must be a 1-D or 2-D array.'
     if min(xray_bandpass) le 0 then message, 'XRAY_BANDPASS must only contain positive values.'
     if size(xray_bandpass, /n_dim) le 1 then Nxray = 1 else Nxray = (size(xray_bandpass, /dim))[1]
   endif

   if n_elements(xray_exposure) ne 0 then begin
     if size(xray_exposure, /type) lt 2 or size(xray_exposure, /type) gt 5 then $
       message, 'XRAY_EXPOSURE must be of type int, float, or double.'
     if size(reform(xray_exposure), /n_dim) ne 1 then $
       message, 'XRAY_EXPOSURE must be a 1-D array.'
     if min(xray_exposure) le 0 then message, 'XRAY_EXPOSURE must only contain positive values.'
     if size(xray_bandpass, /n_dim) le 1 then if n_elements(xray_exposure) ne 1 then $
       message, 'XRAY_EXPOSURE must have the same number of elements as XRAY_BANDPASS'
     if size(xray_bandpass, /n_dim) gt 1 then if n_elements(xray_exposure) ne (size(xray_bandpass, /dim))[1] then $
       message, 'XRAY_EXPOSURE must have the same number of elements as the second dimension of XRAY_BANDPASS'
   endif

   if n_elements(arf_E_lo) ne 0 then begin
     if size(arf_E_lo, /type) lt 2 or size(arf_E_lo, /type) gt 5 then $
       message, 'ARF_E_LO must be of type int, float, or double.'
     if size(reform(arf_E_lo), /n_dim) ne 1 then $
       message, 'ARF_E_LO must be a 1-D array.'
     if min(arf_E_lo) le 0 then message, 'ARF_E_LO must only contain positive values.'
   endif

   if n_elements(arf_E_hi) ne 0 then begin
     if size(arf_E_hi, /type) lt 2 or size(arf_E_hi, /type) gt 5 then $
       message, 'ARF_E_HI must be of type int, float, or double.'
     if size(reform(arf_E_hi), /n_dim) ne 1 then $
       message, 'ARF_E_HI must be a 1-D array.'
     if min(arf_E_hi) le 0 then message, 'ARF_E_HI must only contain positive values.'
     if n_elements(arf_E_lo) ne 0 then if n_elements(arf_E_hi) ne n_elements(arf_E_lo) then $
       message, 'ARF_E_HI must have the same number of elements as ARF_E_LO.'
   endif

   if n_elements(arf_specresp) ne 0 then begin
     if size(arf_specresp, /type) lt 2 or size(arf_specresp, /type) gt 5 then $
       message, 'ARF_SPECRESP must be of type int, float, or double.'
     if size(reform(arf_specresp), /n_dim) ne 1 then $
       message, 'ARF_SPECRESP must be a 1-D array.'
     if min(arf_specresp) le 0 then message, 'ARF_SPECRESP must only contain positive values.'
     if n_elements(arf_E_lo) ne 0 then if n_elements(arf_specresp) ne n_elements(arf_E_lo) then $
       message, 'ARF_SPECRESP must have the same number of elements as ARF_E_LO.'
   endif

   if n_elements(galactic_nh) ne 0 then begin
     if size(galactic_nh, /type) lt 2 or size(galactic_nh, /type) gt 5 then $
       message, 'GALACTIC_NH must be of type int, float, or double.'
     if size(galactic_nh, /n_dim) ne 0 then $
       message, 'GALACTIC_NH must be a scalar.'
     if min(galactic_nh) lt 0 then message, 'GALACTIC_NH must only contain non-negative values.'
   endif

   if n_elements(input_dir) ne 0 then begin
     if size(input_dir, /type) ne 7 then $
       message, 'INPUT_DIR must be of type string.'
     if size(input_dir, /n_dim) ne 0 then $
       message, 'INPUT_DIR must be a scalar.'
     if ~file_test(input_dir, /directory) then message, 'INPUT_DIR is not a valid directory.'
   endif
 endif

 if n_elements(filter_labels) eq 0 then $
   filter_labels=['GALEX_FUV', 'GALEX_NUV', 'SDSS_u', 'SDSS_g', 'SDSS_r', 'SDSS_i', 'SDSS_z', $
                  '2MASS_J', '2MASS_H', '2MASS_Ks', 'IRAC_CH1', 'IRAC_CH2', 'IRAC_CH3', $
                  'IRAC_CH4', 'MIPS_CH1', 'PACS_green', 'PACS_red', 'SPIRE_250', 'SPIRE_350', 'SPIRE_500']
 if n_elements(redshift) eq 0 then redshift = 0.0
 if n_elements(lumin_dist) eq 0 then lumin_dist = 10.d
 if strupcase(config.ATTEN_CURVE) eq 'DOORE21' then $
   if n_elements(input_dir) eq 0 then $
     message, 'INPUT_DIR must be specified if using the DOORE21 attenuation model.'


;====== Generating models =========
 ; Generate stellar emission models
 case strupcase(config.SSP) of
   'PEGASE': begin
        case strupcase(config.SFH) of
          'NON-PARAMETRIC': stellar_models = binned_stellar_models(filter_labels=filter_labels, $
                                                                   redshift=redshift, _EXTRA=config)
          'PARAMETRIC': stellar_models = !values.D_NaN
        endcase
      end
   'NONE': stellar_models = !values.D_NaN
 endcase


 ; Load files for dust attenuation
 case strupcase(config.ATTEN_CURVE) of
   'DOORE21': begin
        restore, !lightning_dir+'models/attenuation/doore21/tuffs_coeff.sav'
        doore21_Lbol_abs_table = mrdfits(input_dir+'lightning_output/doore21_Lbol_abs_table/redshift_'+$
                                         string(redshift,f='(f0.6)')+'.fits', 1, /silent)
        atten_models = {tuffs_coeff : {abulge : temporary(abulge),       $
                                       adisk  : temporary(adisk),        $
                                       atdisk : temporary(atdisk)},      $
                        doore21_Lbol_abs_table  : doore21_Lbol_abs_table $
                        }
      end
   else: atten_models = !values.D_NaN
 endcase


 ; Generate dust emission models
 case strupcase(config.DUST_MODEL) of
   'DL07': dust_models = dl07_models(filter_labels=filter_labels, redshift=redshift)
   'NONE': dust_models = !values.D_NaN
 endcase


 ; Generate X-ray emission models
 if keyword_set(config.XRAY_EMISSION) then begin

   case strupcase(config.XRAY_UNIT) of

     'COUNTS': xray_models = xrb_xagn_models(xray_bandpass, xray_exposure=xray_exposure, arf_e_lo=arf_e_lo, $
                                             arf_e_hi=arf_e_hi, arf_specresp=arf_specresp, redshift=redshift, $
                                             lumin_dist=lumin_dist, galactic_nH=galactic_nH, $
                                             _EXTRA=config)

     'FLUX': xray_models = xrb_xagn_models(xray_bandpass, redshift=redshift, $
                                           lumin_dist=lumin_dist, galactic_nH=galactic_nH, $
                                           _EXTRA=config)
   endcase

 endif else xray_models = !values.D_NaN


 ; Generate AGN emission models
 case strupcase(config.AGN_MODEL) of
   'SKIRTOR': agn_models = skirtor_models(filter_labels=filter_labels, redshift=redshift)
   'NONE': agn_models = !values.D_NaN
 endcase


 models = {filter_labels  : filter_labels,   $
           stellar_models : stellar_models,  $
           atten_models   : atten_models,    $
           dust_models    : dust_models,     $
           xray_models    : xray_models,     $
           agn_models     : agn_models       $
           }

 return, models

end
