pro lightning_input, input_file_fits, cosmology=cosmology, xray_emission=xray_emission, xray_unit=xray_unit, xray_unc=xray_unc
;+
; Name
; ----
;   LIGHTNING_INPUT
;
; Purpose
; -------
;   Reads in the SED data that is to be fit by Lightning from a FITS file. The
;   input fluxes are converted to luminosities using the user supplied redshift
;   or luminosity distance. These luminosities are then saved to a ``.sav`` file
;   for each SED to be fit by Lightning.
;
; Calling Sequence
; ----------------
;   ::
;
;       lightning_input, input_file_fits [, cosmology = , /xray_emission, xray_unit = , xray_unc = ]
;
; Input
; -----
;   ``input_file_fits`` : string scalar
;       The name (including path) to the FITS file containing the SED fluxes
;       and distances (or redshifts) in a data table. The table must be in
;       the first extension. (See :ref:`fits-format-label` for full details,
;       required contents, and format of the data table.)
;
; Optional Inputs
; ---------------
;   ``cosmology`` : structure
;       A structure containing the cosmology parameters ``H0``, ``OMEGA_M``, ``LAMBDA0``,
;       ``Q0``, and ``K`` in individual tags.
;       (Default = ``{H0: 70, OMEGA_M: 0.3, LAMBDA0: 0.7, Q0: -0.55, K: 0}``)
;   ``xray_emission`` : flag
;       If set, Lightning will search for and load the additional data products necessary to fit
;       an X-ray emission model.
;   ``xray_unit`` : string scalar
;       The type of X-ray data to use in fitting. Current options are ``'COUNTS'`` and ``'FLUX'``.
;       (Default = ``'COUNTS'``)
;   ``xray_unc`` : string scalar
;       The errors to assume if using X-ray count data. Current options are ``'SQRT'``, ``'GEHRELS'``,
;       and ``'USER'``. (Default = ``'GEHRELS'``)
;
; Output
; ------
;   A save file (``<input_file_dir>/lightning_output/input_sav_files/lightning_input_<sed_id>.sav``)
;   for each SED that contains the data converted to luminosities.
;   The data for each SED are stored in a structure with the following format:
;
;   =======================     =================     ====================================================================================================================================
;   TAG                         TYPE                  DESCRIPTION
;   =======================     =================     ====================================================================================================================================
;   SED_ID                      string                Unique SED identifier
;   LNU_OBS                     double(Nfilters)      Luminosities of the SED for each set of filters in terms of :math:`L_\nu` :math:`[\rm{L_{\odot}\ Hz^{-1}}]`
;   LNU_UNC                     double(Nfilters)      Uncertainties associated with the luminosities :math:`[\rm{L_{\odot}\ Hz^{-1}}]`
;   FILTER_LABELS               string(Nfilters)      Names of the filters associated with each luminosity
;   REDSHIFT                    double                Redshift of the SED (if ``LUMIN_DIST`` was specified this is set to ``0``)
;   LUMIN_DIST                  double                Luminosity distance of the SED as input or converted from the redshift :math:`[\rm{Mpc}]`
;   XRAY_LNU_OBS [1]_           double(Nxray)         X-ray luminosities for each bandpass in terms of :math:`L_\nu` :math:`[\rm{L_{\odot}\ Hz^{-1}}]`
;   XRAY_LNU_UNC [1]_           double(Nxray)         Uncertainties associated with the X-ray luminosities :math:`[\rm{L_{\odot}\ Hz^{-1}}]`
;   XRAY_BANDPASS [1]_ [2]_     double(2, Nxray)      Bandpasses of X-ray observations: first column contains the lower energy bound, second column contains the upper. :math:`[\rm{keV}]`
;   XRAY_EXPOSURE [2]_          double(Nxray)         Exposure times of X-ray observations, one per band :math:`[\rm{s}]`
;   NET_COUNTS [2]_             double(Nxray)         Net counts in each X-ray band :math:`[\rm{counts}]`
;   NET_COUNTS_UNC [2]_         double(Nxray)         User-supplied uncertainty in net counts, if any :math:`[\rm{counts}]`
;   GALACTIC_NH [1]_ [2]_       double                Galactic (i.e., Milky Way) HI column density along the line of sight :math:`[10^{20}\ \rm{cm}^{-2}]`
;   ARF_E_LO [2]_               double(Nchannels)     Lower energy bounds of each channel in the ARF :math:`[\rm{keV}]`
;   ARF_E_HI [2]_               double(Nchannels)     Upper energy bounds of each channel in the ARF :math:`[\rm{keV}]`
;   ARF_SPECRESP [2]_           double(Nchannels)     Spectral response of the ARF at each channel :math:`[\rm{cm^2}]`
;   =======================     =================     ====================================================================================================================================
;
;   .. [1] These tags are present if the X-ray module is used with ``XRAY_UNIT = 'FLUX'``
;   .. [2] These tags are present if the X-ray module is used with ``XRAY_UNIT = 'COUNTS'``
;
; Modification History
; --------------------
;   - 2022/03/23: Created (Keith Doore)
;   - 2022/04/25: Updated to allow for Xray ARF file path (Keith Doore)
;   - 2022/04/28: Added unique ``SED_ID`` enforcement (Keith Doore)
;   - 2022/05/18: Added optional cosmology input to allow for unique cosmologies in ``lumdist`` (Keith Doore)
;   - 2022/05/18: Added check to make sure xray files existed (Keith Doore)
;   - 2022/05/19: Replace ``n_elements(where())`` statements in logicals with ``total()`` (Keith Doore)
;   - 2022/06/27: Removed missing data from SED before saving to file (Keith Doore)
;   - 2022/06/27: Fixed issue where xray_spec/arf file strings need to be trimmed of padded blanks (Keith Doore)
;   - 2022/08/01: Added ``/silent`` to ``mrdfits`` (Keith Doore)
;   - 2022/08/02: Undid removal of missing data from SED before saving to file to make post-processing simpler (Keith Doore)
;   - 2022/08/02: Added ability to use ``redshift`` if ``lumin_dist`` is ``0`` for some SEDs (Keith Doore)
;   - 2022/08/02: Added check to require both ``XRAY_SPEC_FILE`` and ``XRAY_ARF_FILE`` tags if including X-ray data (Keith Doore)
;   - 2022/08/09: Added ``GALACTIC_NH`` as input for xray emission (Keith Doore)
;   - 2022/08/10: Fixed bug with redshift distance being computed even with no redshifts (Keith Doore)
;   - 2022/09/01: Added handling for user-supplied X-ray count uncertainties (Erik B. Monson)
;   - 2022/09/14: Updated to allow fitting with X-ray fluxes (Erik B. Monson)
;   - 2022/09/15: Updated documentation (Keith Doore)
;   - 2022/10/24: Allowed for negative flux inputs (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(input_file_fits) eq 0 then message, 'Variable is undefined: INPUT_FILE_FITS.'
 if size(input_file_fits, /type) ne 7 then message, 'INPUT_FILE_FITS must be of type string.'
 if size(input_file_fits, /n_dim) ne 0 then message, 'INPUT_FILE_FITS must be a scalar.'
 if ~file_test(input_file_fits) then message, 'No input data file exists. '+$
                                        'Please confirm specified file name leads to the correct file.'

; Read in data
 data = mrdfits(input_file_fits, 1, /silent)
 tags = tag_names(data)


; Check for tags, any errors, and format
 if total(tags eq 'FNU_OBS') eq 0 then $
   message, 'FNU_OBS tag must be included in the input file table.'
 if size(data.FNU_OBS, /type) lt 4 or size(data.FNU_OBS, /type) gt 5 then $
   message, 'FNU_OBS tag in input file table must be of type float or double.'
 if size(data.FNU_OBS, /n_dim) lt 1 or size(data.FNU_OBS, /n_dim) gt 2 then $
   message, 'FNU_OBS tag in input file table must be a 1-D or 2-D array.'
 Fnu_obs = data.FNU_OBS
 Nfilters = (size(data.FNU_OBS, /dim))[0]
 if size(data.FNU_OBS, /n_dim) eq 1 then Nsed = 1 else Nsed = (size(data.FNU_OBS, /dim))[1]

 if total(tags eq 'SED_ID') eq 0 then begin
   sed_id = strtrim(sindgen(Nsed), 2)
 endif else begin
   if size(data.SED_ID, /type) ne 7 then $
     message, 'SED_ID tag in input file table must be of type string.'
   if size(data.SED_ID, /n_dim) gt 1 then message, 'SED_ID tag in input file table must be a 1-D array.'
   if size(reform(data.SED_ID), /dim) ne Nsed then $
     message, 'SED_ID tag in input file table must be present for each SED.'
   if n_elements(uniq(data.sed_id, sort(data.sed_id))) ne Nsed then $
     message, 'Each SED_ID tag in input file table must be unique.'
   sed_id = strtrim(data.SED_ID, 2)
 endelse

 if total(tags eq 'FNU_UNC') eq 0 then $
   message, 'FNU_UNC tag must be included in the input file table.'
 if size(data.FNU_UNC, /type) lt 4 or size(data.FNU_UNC, /type) gt 5 then $
   message, 'FNU_UNC tag in input file table must be of type float or double.'
 if size(data.FNU_UNC, /n_dim) lt 1 or size(data.FNU_UNC, /n_dim) gt 2 then $
   message, 'FNU_UNC tag in input file table must be a 1-D or 2-D array.'
 if total(size(data.FNU_UNC, /dim) ne size(data.FNU_OBS, /dim)) gt 0 then $
   message, 'FNU_UNC tag in input file table must have the same size as FNU_OBS.'
 if min(data.FNU_UNC) lt 0 then message, 'FNU_UNC tag in input file table must only contain non-negative values.'
 Fnu_unc = data.FNU_UNC

 if total(tags eq 'FILTER_LABELS') eq 0 then $
   message, 'FILTER_LABELS tag must be included in the input file table.'
 if size(data.FILTER_LABELS, /type) ne 7 then $
   message, 'FILTER_LABELS tag in input file table must be of type string.'
 if size(data.FILTER_LABELS, /n_dim) lt 1 or size(data.FILTER_LABELS, /n_dim) gt 2 then $
   message, 'FILTER_LABELS tag in input file table must be a 1-D or 2-D array.'
 if total(size(data.FILTER_LABELS, /dim) ne size(data.FNU_OBS, /dim)) gt 0 then $
   message, 'FILTER_LABELS in input file table must have the same size as FNU_OBS.'
 filter_labels = data.FILTER_LABELS
; Check if the labels are correct names for Lightning
 json = JSON_PARSE(!lightning_dir+'filters/filters.json')
 avail_filters = (json.keys()).toarray()
 uniq_filters = strtrim(filter_labels[uniq(filter_labels, sort(filter_labels))], 2)
 for i=0, (n_elements(uniq_filters)-1) do begin
   if total(uniq_filters[i] eq avail_filters) ne 1 then $
     message, "FILTER_LABELS value of '"+uniq_filters[i]+"' is not a recognized "+$
              "filter label in Lightning."
 endfor

 if total(tags eq 'LUMIN_DIST') eq 0 then begin
   if total(tags eq 'REDSHIFT') eq 0 then $
     message, 'REDSHIFT or LUMIN_DIST tag must be included in the input file table.'
   if size(data.REDSHIFT, /type) lt 2 or size(data.REDSHIFT, /type) gt 5 then $
     message, 'REDSHIFT tag in input file table must be of type int, float, or double.'
   if size(data.REDSHIFT, /n_dim) gt 1 then message, 'REDSHIFT tag in input file table must be a 1-D array.'
   if size(reform(data.REDSHIFT), /dim) ne Nsed then $
     message, 'REDSHIFT tag in input file table must be present for each SED.'
   if min(data.REDSHIFT) le 0 then message, 'REDSHIFT tag in input file table must only contain positive values if '+$
                                            'LUMIN_DIST tag is not in input file table.'
   if min(lumdist(data.REDSHIFT, /silent, _extra=cosmology)) le 10 then $
     message, 'REDSHIFT tag in input file table has very small redshifts, which result in nearby luminosity '+$
              'distances (<10 Mpc). Please input the luminosity distance instead of the redshift, as the '+$
              'luminosity distances from these redshifts are likely inaccurate.'
   redshift = data.REDSHIFT
   lumin_dist = dblarr(Nsed)
 endif else begin
   if size(data.LUMIN_DIST, /type) lt 2 or size(data.LUMIN_DIST, /type) gt 5 then $
     message, 'LUMIN_DIST tag in input file table must be of type int, float, or double.'
   if size(data.LUMIN_DIST, /n_dim) gt 1 then message, 'LUMIN_DIST tag in input file table must be a 1-D array.'
   if size(reform(data.LUMIN_DIST), /dim) ne Nsed then $
     message, 'LUMIN_DIST tag in input file table must be present for each SED.'
   lumin_dist = data.LUMIN_DIST
   redshift = dblarr(Nsed)

   if total(tags eq 'REDSHIFT') ne 0 then begin
     if min(data.LUMIN_DIST) lt 0 then message, 'LUMIN_DIST tag in input file table must only contain non-negative values '+$
                                                'if REDSHIFT tag is in input file table.'
     if size(data.REDSHIFT, /type) lt 2 or size(data.REDSHIFT, /type) gt 5 then $
       message, 'REDSHIFT tag in input file table must be of type int, float, or double.'
     if size(data.REDSHIFT, /n_dim) gt 1 then message, 'REDSHIFT tag in input file table must be a 1-D array.'
     if size(reform(data.REDSHIFT), /dim) ne Nsed then $
       message, 'REDSHIFT tag in input file table must be present for each SED.'
     if min(data.REDSHIFT) lt 0 then message, 'REDSHIFT tag in input file table must only contain non-negative values '+$
                                              'if LUMIN_DIST tag is in input file table.'
     ; If a luminosity is 0, then a redshift must be given instead
     red_idc = where(lumin_dist eq 0, Nred, /null)
     if Nred gt 0 then $
       if min(lumdist((data.REDSHIFT)[red_idc], /silent, _extra=cosmology)) le 10 then $
         message, 'REDSHIFT tag in input file table has very small redshifts for some SEDs without luminosity distances, '+$
                  'which result in nearby luminosity distances (<10 Mpc). Please input the luminosity distance '+$
                  'for these SEDs instead of the redshift, as the luminosity distances from these redshifts are likely inaccurate.'
     redshift[red_idc] = (data.REDSHIFT)[red_idc]
   endif else begin
     if min(data.LUMIN_DIST) le 0 then message, 'LUMIN_DIST tag in input file table must only contain positive values '+$
                                                'if REDSHIFT tag is not in input file table.'
   endelse
 endelse

 ; We'll only check the X-ray tags for correctness if X-ray emission is turned on
 if keyword_set(xray_emission) then begin
 
   if n_elements(xray_unit) ne 0 then begin
     if size(xray_unit, /type) ne 7 then message, 'XRAY_UNIT must be of type string.'
     if size(xray_unit, /dim) ne 0 then message, 'XRAY_UNIT must be a scalar.'
     if total(strupcase(xray_unit) eq ['COUNTS', 'FLUX']) ne 1 then $
       message, "XRAY_UNIT must be set to either 'COUNTS' or 'FLUX'."
   endif else xray_unit = 'COUNTS'

   if n_elements(xray_unc) ne 0 then begin
     if size(xray_unc, /type) ne 7 then message, 'XRAY_UNC must be of type string.'
     if size(xray_unc, /dim) ne 0 then message, 'XRAY_UNC must be a scalar.'
     if total(strupcase(xray_unc) eq ['SQRT', 'GEHRELS', 'USER']) ne 1 then $
       message, "XRAY_UNC must be set to either 'SQRT', 'GEHRELS', or 'USER'."
   endif else xray_unc = 'GEHRELS'

   if total(tags eq 'GALACTIC_NH') eq 0 then begin
     message, 'GALACTIC_NH tag must be included in the input file table if X-ray data is included.'
   endif else begin
     if size(data.GALACTIC_NH, /type) lt 2 or size(data.GALACTIC_NH, /type) gt 5 then $
       message, 'GALACTIC_NH tag in input file table must be of type int, float, or double.'
     if size(data.GALACTIC_NH, /n_dim) gt 1 then message, 'GALACTIC_NH tag in input file table must be a 1-D array.'
     if size(reform(data.GALACTIC_NH), /dim) ne Nsed then $
       message, 'GALACTIC_NH tag in input file table must be present for each SED.'
     if min(data.GALACTIC_NH) lt 0 then $
       message, base_err+'GALACTIC_NH tag in input file table must be a non-negative value.'
     galactic_nh = data.GALACTIC_NH
   endelse

   case strupcase(xray_unit) of

     'COUNTS': begin
         if total(tags eq 'XRAY_SPEC_FILE') eq 0 then begin
           message, "XRAY_SPEC_FILE tag must be included in the input file table if XRAY_UNIT is set to 'COUNTS'."
         endif else begin
           if size(data.XRAY_SPEC_FILE, /type) ne 7 then $
             message, 'XRAY_SPEC_FILE tag in input file table must be of type string.'
           if size(data.XRAY_SPEC_FILE, /n_dim) gt 1 then message, 'XRAY_SPEC_FILE tag in input file table must be a 1-D array.'
           if size(reform(data.XRAY_SPEC_FILE), /dim) ne Nsed then $
             message, 'XRAY_SPEC_FILE tag in input file table must be present for each SED.'
           foreach xray_file, data.XRAY_SPEC_FILE do begin
             if ~file_test(xray_file) then $
               path_message = 'XRAY_SPEC_FILE tag in input file table lists a file that does not exist. Please '+$
                              'specify the correct file. (Incorrect file: '+xray_file+')'
           endforeach
           xray_spec_file = data.XRAY_SPEC_FILE
         endelse

         if total(tags eq 'XRAY_ARF_FILE') eq 0 then begin
           message, "XRAY_ARF_FILE tag must be included in the input file table if XRAY_UNIT is set to 'COUNTS'."
         endif else begin
           if size(data.XRAY_ARF_FILE, /type) ne 7 then $
             message, 'XRAY_ARF_FILE tag in input file table must be of type string.'
           if size(data.XRAY_ARF_FILE, /n_dim) gt 1 then message, 'XRAY_ARF_FILE tag in input file table must be a 1-D array.'
           if size(reform(data.XRAY_ARF_FILE), /dim) ne Nsed then $
             message, 'XRAY_ARF_FILE tag in input file table must be present for each SED.'
           foreach arf_file, data.XRAY_ARF_FILE do begin
             if ~file_test(arf_file) then $
               path_message = 'XRAY_ARF_FILE tag in input file table lists a file that does not exist. Please '+$
                              'specify the correct file. (Incorrect file: '+xray_file+')'
           endforeach
           xray_arf_file = data.XRAY_ARF_FILE
         endelse
       end

     'FLUX': begin
         if total(tags eq 'XRAY_FLUX') eq 0 then begin
           message, "XRAY_FLUX tag must be included in the input file table if XRAY_UNIT is set to 'FLUX'."
         endif else begin
           if size(data.XRAY_FLUX, /type) lt 2 or size(data.XRAY_FLUX, /type) gt 5 then $
             message, 'XRAY_FLUX tag in input file table must be of type int, float, or double.'
           if size(data.XRAY_FLUX, /n_dim) lt 1 or size(data.XRAY_FLUX, /n_dim) gt 2 then $
             message, 'XRAY_FLUX tag in input file table must be a 1-D or 2-D array.'
           Nxray = (size(data.XRAY_FLUX, /dim))[0]
           xray_flux = data.XRAY_FLUX
         endelse
 
         if total(tags eq 'XRAY_FLUX_UNC') eq 0 then begin
           message, "XRAY_FLUX_UNC tag must be included in the input file table if XRAY_UNIT is set to 'FLUX'."
         endif else begin
           if size(data.XRAY_FLUX_UNC, /type) lt 2 or size(data.XRAY_FLUX_UNC, /type) gt 5 then $
             message, 'XRAY_FLUX_UNC tag in input file table must be of type int, float, or double.'
           if size(data.XRAY_FLUX_UNC, /n_dim) lt 1 or size(data.XRAY_FLUX_UNC, /n_dim) gt 2 then $
             message, 'XRAY_FLUX_UNC tag in input file table must be a 1-D or 2-D array.'
           if min(data.XRAY_FLUX_UNC) lt 0 then message, 'XRAY_FLUX_UNC tag in input file table must only contain non-negative values.'
           if (size(data.XRAY_FLUX_UNC, /dim))[0] ne Nxray then $
             message, 'XRAY_FLUX_UNC must contain the same number of observations as XRAY_FLUX.'
           xray_flux_unc = data.XRAY_FLUX_UNC
         endelse
 
         if total(tags eq 'XRAY_BANDPASS') eq 0 then begin
           message, "XRAY_BANDPASS tag must be included in the input file table if XRAY_UNIT is set to 'FLUX'."
         endif else begin
           if size(data.XRAY_BANDPASS, /type) lt 2 or size(data.XRAY_BANDPASS, /type) gt 5 then $
             message, 'XRAY_BANDPASS tag in input file table must be of type int, float, or double.'
           if size(data.XRAY_BANDPASS, /n_dim) lt 1 or size(data.XRAY_BANDPASS, /n_dim) gt 3 then $
             message, 'XRAY_BANDPASS tag in input file table must be a 1-D, 2-D, or 3-D array.'
           if min(data.XRAY_BANDPASS) lt 0 then message, 'XRAY_BANDPASS tag in input file table must only contain non-negative values.'
           if (size(data.XRAY_BANDPASS, /dim))[0] ne 2 then $
             message, 'XRAY_BANDPASS must have its first dimension equal to 2.'
           if (size(data.XRAY_BANDPASS, /dim))[1] ne Nxray then $
             message, 'XRAY_BANDPASS must have its second dimension correspond to the same number of observations as XRAY_FLUX.'
           xray_bandpass = data.XRAY_BANDPASS
         endelse
       end
   endcase
 endif


; Extract input file path and make a directory in it to save data
 save_dir = file_dirname(input_file_fits, /mark)+'lightning_output/input_sav_files/'
 file_mkdir, save_dir


; Convert the fluxes to luminosities
 ; If a luminosity is 0, then a redshift was given instead
 red_idc = where(lumin_dist eq 0, Nred, comp=dist_idc, /null)
 DL = dblarr(Nsed)
 DL[dist_idc] = lumin_dist[dist_idc] ; luminosity distance [Mpc]
 if Nred gt 0 then DL[red_idc] = lumdist(redshift[red_idc], /silent, _extra=cosmology) ; luminosity distance [Mpc]
 ; conversion from Fnu in Jy to Lnu in Lsun/Hz
 conv = 4*!dpi*(DL*1.d6*!lightning_cgs.pc)^2*!lightning_cgs.Jy/!lightning_cgs.Lsun
 conv = rebin(reform(conv, 1, Nsed), Nfilters, Nsed)

 Lnu_obs = Fnu_obs*conv
 Lnu_unc = Fnu_unc*conv


; Place the data into a new structure to output
 for i=0, Nsed-1 do begin

  ; Read in Xray data for the SED
   if (not keyword_set(xray_emission)) then begin
     sed_data = {SED_ID        : '',                                 $
                 REDSHIFT      : !values.D_NaN,                      $
                 LUMIN_DIST    : !values.D_NaN,                      $
                 FILTER_LABELS : strarr(Nfilters),                   $
                 LNU_OBS       : !values.D_NaN*dblarr(Nfilters),     $
                 LNU_UNC       : dblarr(Nfilters)                    $
                 }
   endif else begin

     case strupcase(xray_unit) of

       'COUNTS': begin

           xray_data = lightning_xray_input(sed_id[i], strtrim(xray_spec_file[i], 2), strtrim(xray_arf_file[i], 2),$
                                            xray_unc=xray_unc)

           sed_data = {SED_ID         : '',                                 $
                       REDSHIFT       : !values.D_NaN,                      $
                       LUMIN_DIST     : !values.D_NaN,                      $
                       FILTER_LABELS  : strarr(Nfilters),                   $
                       LNU_OBS        : !values.D_NaN*dblarr(Nfilters),     $
                       LNU_UNC        : dblarr(Nfilters),                   $
                       XRAY_BANDPASS  : xray_data.XRAY_BANDPASS,            $
                       XRAY_EXPOSURE  : xray_data.XRAY_EXPOSURE,            $
                       NET_COUNTS     : xray_data.NET_COUNTS,               $
                       NET_COUNTS_UNC : xray_data.NET_COUNTS_UNC,           $
                       GALACTIC_NH    : !values.D_NaN,                      $
                       ARF_E_LO       : xray_data.ARF_E_LO,                 $
                       ARF_E_HI       : xray_data.ARF_E_HI,                 $
                       ARF_SPECRESP   : xray_data.ARF_SPECRESP              $
                       }
         end
       'FLUX': begin

           conv = 4*!dpi*(DL[i]*1.d6*!lightning_cgs.pc)^2/!lightning_cgs.Lsun
           ;conv = rebin(conv, Nxray)
           ;conv = rebin(reform(conv, 1, Nsed), Nxray, Nsed)
           xray_bandpass_nu = xray_bandpass[*, *, i] / (!lightning_cgs.hplanck / !lightning_cgs.keV)
           xray_bandpass_dnu = abs(xray_bandpass_nu[1, *] - xray_bandpass_nu[0, *])

           xray_Lnu_obs = xray_flux[*, i] * conv / xray_bandpass_dnu
           xray_Lnu_unc = xray_flux_unc[*, i] * conv / xray_bandpass_dnu

           sed_data = {SED_ID         : '',                                 $
                       REDSHIFT       : !values.D_NaN,                      $
                       LUMIN_DIST     : !values.D_NaN,                      $
                       FILTER_LABELS  : strarr(Nfilters),                   $
                       LNU_OBS        : !values.D_NaN*dblarr(Nfilters),     $
                       LNU_UNC        : dblarr(Nfilters),                   $
                       XRAY_BANDPASS  : xray_bandpass[*, *, i],             $
                       XRAY_LNU_OBS   : xray_Lnu_obs,                       $
                       XRAY_LNU_UNC   : xray_Lnu_unc,                       $
                       GALACTIC_NH    : !values.D_NaN                       $
                       }
         end

     endcase

     sed_data.GALACTIC_NH  = galactic_nh[i]

   endelse

   sed_data.SED_ID         = sed_id[i]
   sed_data.REDSHIFT       = redshift[i]
   sed_data.LUMIN_DIST     = DL[i]
   sed_data.FILTER_LABELS  = filter_labels[*, i]
   sed_data.LNU_OBS        = lnu_obs[*, i]
   sed_data.LNU_UNC        = lnu_unc[*, i]

   save, sed_data, filename=save_dir+'lightning_input_'+sed_id[i]+'.sav'

 endfor

end
