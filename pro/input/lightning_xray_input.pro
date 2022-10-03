function lightning_xray_input, sed_id, xray_spec_file, xray_arf_file, xray_unc=xray_unc
;+
; Name
; ----
;   LIGHTNING_XRAY_INPUT
;
; Purpose
; -------
;   Reads in the necessary X-ray data from a spectral file (e.g., the outputs of
;   ACISExtract) and Auxiliary Response Function file.
;
; Calling Sequence
; ----------------
;   ::
;
;       xray_data = lightning_xray_input(sed_id, xray_spec_file, xray_arf_file [, xray_unc =])
;
; Inputs
; ------
;   ``sed_id`` : string scalar
;       A unique SED identifier
;   ``xray_spec_file`` : string scalar
;       The name (including path) to the FITS file containing the X-ray spectral
;       data table, which must be in the first extension.
;       The data table must have the following structure:
;
;       ==============     ===================     ==============================================================
;       TAG                TYPE                    DESCRIPTION
;       ==============     ===================     ==============================================================
;       ENERG_LO           float/double(Nxray)     Lower energy bound of each observation band :math:`[\rm{keV}]`
;       ENERG_HI           float/double(Nxray)     Upper energy bound of each observation band :math:`[\rm{keV}]`
;       NET_COUNTS         float/double(Nxray)     Net counts in each band :math:`[\rm{counts}]`
;       NET_COUNTS_UNC     float/double(Nxray)     Optional, uncertainty on net counts :math:`[\rm{counts}]`
;       EXPOSURE           float/double(Nxray)     Exposure time of each band :math:`[\rm{s}]`
;       ==============     ===================     ==============================================================
;
;   ``xray_arf_file`` : string scalar
;       The name (including path) to the FITS file containing the ARF
;       data table, which must be in the first extension.
;       The data table must have the following structure:
;
;       ========     =======================     ======================================================
;       TAG          TYPE                        DESCRIPTION
;       ========     =======================     ======================================================
;       ENERG_LO     float/double(Nchannels)     Lower energy bounds of each channel :math:`[\rm{keV}]`
;       ENERG_HI     float/double(Nchannels)     Upper energy bounds of each channel :math:`[\rm{keV}]`
;       SPECRESP     float/double(Nchannels)     Spectral response at each channel :math:`[\rm{cm}^2]`
;       ========     =======================     ======================================================
;
; Optional Input
; --------------
;   ``xray_unc`` : string scalar
;       The errors to assume if using X-ray count data. Current options are ``'SQRT'``, ``'GEHRELS'``,
;       and ``'USER'``. (Default = ``'GEHRELS'``)
;
; Output
; ------
;   ``xray_data`` : structure
;       This structure includes the X-ray data and associated ARF data.
;       The full description of the structure is as follows:
;
;       ==============     =================     ====================================================================================================================================
;       TAGS               TYPE                  DESCRIPTION
;       ==============     =================     ====================================================================================================================================
;       XRAY_BANDPASS      double(2, Nxray)      Bandpasses of X-ray observations: first column contains the lower energy bound, second column contains the upper. :math:`[\rm{keV}]`
;       XRAY_EXPOSURE      double(Nxray)         Exposure times of X-ray observations, one per band :math:`[\rm{s}]`
;       NET_COUNTS         double(Nxray)         Net counts in each X-ray band :math:`[\rm{counts}]`
;       NET_COUNTS_UNC     double(Nxray)         Uncertainty on the net counts in each X-ray band :math:`[\rm{counts}]`
;       ARF_E_LO           double(Nchannels)     Lower energy bounds of each channel in the ARF :math:`[\rm{keV}]`
;       ARF_E_HI           double(Nchannels)     Upper energy bounds of each channel in the ARF :math:`[\rm{keV}]`
;       ARF_SPECRESP       double(Nchannels)     Spectral response of the ARF at each channel :math:`[\rm{cm}^2]`
;       ==============     =================     ====================================================================================================================================
;
; Modification History
; --------------------
;   - 2022/06/08: Created (Erik B. Monson)
;   - 2022/06/08: Updated error handling (Keith Doore)
;   - 2022/09/01: Added handling for user-supplied X-ray count uncertainties (Erik B. Monson)
;   - 2022/09/01: Updated documentation (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(sed_id) eq 0 then message, 'Variable is undefined: SED_ID.'
 if size(sed_id, /type) ne 7 then message, 'SED_ID must be of type string.'
 if size(sed_id, /n_dim) ne 0 then message, 'SED_ID must be a scalar.'

 if n_elements(xray_spec_file) eq 0 then message, 'Variable is undefined: XRAY_SPEC_FILE.'
 if size(xray_spec_file, /type) ne 7 then message, 'XRAY_SPEC_FILE must be of type string.'
 if size(xray_spec_file, /n_dim) ne 0 then message, 'XRAY_SPEC_FILE must be a scalar.'
 if ~file_test(xray_spec_file) then message, 'No XRAY_SPEC_FILE file exists for SED_ID '+sed_id+'. '+$
                             'Please confirm specified file name leads to the correct file.'

 if n_elements(xray_arf_file) eq 0 then message, 'Variable is undefined: XRAY_ARF_FILE.'
 if size(xray_arf_file, /type) ne 7 then message, 'XRAY_ARF_FILE must be of type string.'
 if size(xray_arf_file, /n_dim) ne 0 then message, 'XRAY_ARF_FILE must be a scalar.'
 if ~file_test(xray_arf_file) then message, 'No XRAY_ARF_FILE file exists for SED_ID '+sed_id+'. '+$
                             'Please confirm specified file name leads to the correct file.'

 if n_elements(xray_unc) ne 0 then begin
   if size(xray_unc, /type) ne 7 then message, 'XRAY_UNC must be of type string.'
   if size(xray_unc, /dim) ne 0 then message, 'XRAY_UNC must be a scalar.'
   if total(strupcase(xray_unc) eq ['SQRT', 'GEHRELS', 'USER']) ne 1 then $
     message, "XRAY_UNC must be set to either 'SQRT', 'GEHRELS', or 'USER'."
 endif else xray_unc = 'GEHRELS'


; Read in data
 xray_spec = mrdfits(xray_spec_file, 1, /silent)
 xray_arf = mrdfits(xray_arf_file, 1, /silent)
 spec_tags = tag_names(xray_spec)
 arf_tags = tag_names(xray_arf)


; Check for tags, any errors, and format
 if total(spec_tags eq 'NET_COUNTS') eq 0 then $
   message, 'NET_COUNTS tag must be included in the XRAY_SPEC_FILE file for SED_ID '+sed_id+'.'
 if size(xray_spec.NET_COUNTS, /type) lt 2 or size(xray_spec.NET_COUNTS, /type) gt 5 then $
   message, 'NET_COUNTS tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must be of type int, float or double.'
 if size(xray_spec.NET_COUNTS, /n_dim) gt 1 then $
   message, 'NET_COUNTS tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must be a scalar or 1D array.'
 if min(xray_spec.NET_COUNTS) lt 0 then $
   message, 'NET_COUNTS tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must only contain positive values.'

 if strupcase(xray_unc) eq 'USER' then begin
   if (total(spec_tags eq 'NET_COUNTS_UNC') eq 1) then begin
     if size(xray_spec.NET_COUNTS_UNC, /type) lt 2 or size(xray_spec.NET_COUNTS, /type) gt 5 then $
        message, 'NET_COUNTS_UNC tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must be of type int, float or double.'
     if size(xray_spec.NET_COUNTS_UNC, /n_dim) gt 1 then $
        message, 'NET_COUNTS_UNC tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must be a scalar or 1D array.'
     if min(xray_spec.NET_COUNTS_UNC) lt 0 then $
        message, 'NET_COUNTS_UNC tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must only contain positive values.'
   endif else begin
     message, "XRAY_UNC is set to 'USER', but no uncertainties were found in the X-ray spectral files. "+$
              "Please specify the uncertainties in the NET_COUNTS_UNC column of the X-ray spectral files."
   endelse
 endif

 if total(spec_tags eq 'ENERG_LO') eq 0 then $
   message, 'ENERG_LO tag must be included in the XRAY_SPEC_FILE file for SED_ID '+sed_id+'.'
 if size(xray_spec.ENERG_LO, /type) lt 2 or size(xray_spec.ENERG_LO, /type) gt 5 then $
   message, 'ENERG_LO tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must be of type int, float or double.'
 if size(xray_spec.ENERG_LO, /n_dim) gt 1 then $
   message, 'ENERG_LO tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must be a scalar or 1D array.'
 if total(size(xray_spec.ENERG_LO, /dim) ne size(xray_spec.NET_COUNTS, /dim)) gt 0 then $
   message, 'ENERG_LO tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must have the same size as the NET_COUNTS tag.'
 if min(xray_spec.ENERG_LO) le 0 then $
   message, 'ENERG_LO tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must only contain positive values.'

 if total(spec_tags eq 'ENERG_HI') eq 0 then $
   message, 'ENERG_HI tag must be included in the XRAY_SPEC_FILE file for SED_ID '+sed_id+'.'
 if size(xray_spec.ENERG_HI, /type) lt 2 or size(xray_spec.ENERG_HI, /type) gt 5 then $
   message, 'ENERG_HI tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must be of type int, float or double.'
 if size(xray_spec.ENERG_HI, /n_dim) gt 1 then $
   message, 'ENERG_HI tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must be a scalar or 1D array.'
 if total(size(xray_spec.ENERG_HI, /dim) ne size(xray_spec.NET_COUNTS, /dim)) gt 0 then $
   message, 'ENERG_HI tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must have the same size as the NET_COUNTS tag.'
 if min(xray_spec.ENERG_HI) le 0 then $
   message, 'ENERG_HI tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must only contain positive values.'

 if total(spec_tags eq 'EXPOSURE') eq 0 then $
   message, 'EXPOSURE tag must be included in the XRAY_SPEC_FILE file for SED_ID '+sed_id+'.'
 if size(xray_spec.EXPOSURE, /type) lt 2 or size(xray_spec.EXPOSURE, /type) gt 5 then $
   message, 'EXPOSURE tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must be of type int, float or double.'
 if size(xray_spec.EXPOSURE, /n_dim) gt 1 then $
   message, 'EXPOSURE tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must be a scalar or 1D array.'
 if total(size(xray_spec.EXPOSURE, /dim) ne size(xray_spec.NET_COUNTS, /dim)) gt 0 then $
   message, 'EXPOSURE tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must have the same size as the NET_COUNTS tag.'
 if min(xray_spec.EXPOSURE) le 0 then $
   message, 'EXPOSURE tag in the XRAY_SPEC_FILE file for SED_ID '+sed_id+' must only contain positive values.'


 if total(arf_tags eq 'SPECRESP') eq 0 then $
   message, 'SPECRESP tag must be included in the XRAY_ARF_FILE file for SED_ID '+sed_id+'.'
 if size(xray_arf.SPECRESP, /type) lt 2 or size(xray_arf.SPECRESP, /type) gt 5 then $
   message, 'SPECRESP tag in the XRAY_ARF_FILE file for SED_ID '+sed_id+' must be of type int, float or double.'
 if size(xray_arf.SPECRESP, /n_dim) gt 1 then $
   message, 'SPECRESP tag in the XRAY_ARF_FILE file for SED_ID '+sed_id+' must be a scalar or 1D array.'
 if min(xray_arf.SPECRESP) le 0 then $
   message, 'SPECRESP tag in the XRAY_ARF_FILE file for SED_ID '+sed_id+' must only contain positive values.'

 if total(arf_tags eq 'ENERG_LO') eq 0 then $
   message, 'ENERG_LO tag must be included in the XRAY_ARF_FILE file for SED_ID '+sed_id+'.'
 if size(xray_arf.ENERG_LO, /type) lt 2 or size(xray_arf.ENERG_LO, /type) gt 5 then $
   message, 'ENERG_LO tag in the XRAY_ARF_FILE file for SED_ID '+sed_id+' must be of type int, float or double.'
 if size(xray_arf.ENERG_LO, /n_dim) gt 1 then $
   message, 'ENERG_LO tag in the XRAY_ARF_FILE file for SED_ID '+sed_id+' must be a scalar or 1D array.'
 if total(size(xray_arf.ENERG_LO, /dim) ne size(xray_arf.SPECRESP, /dim)) gt 0 then $
   message, 'ENERG_LO tag in the XRAY_ARF_FILE file for SED_ID '+sed_id+' must have the same size as the SPECRESP tag.'
 if min(xray_arf.ENERG_LO) le 0 then $
   message, 'ENERG_LO tag in the XRAY_ARF_FILE file for SED_ID '+sed_id+' must only contain positive values.'

 if total(arf_tags eq 'ENERG_HI') eq 0 then $
   message, 'ENERG_HI tag must be included in the XRAY_ARF_FILE file for SED_ID '+sed_id+'.'
 if size(xray_arf.ENERG_HI, /type) lt 2 or size(xray_arf.ENERG_HI, /type) gt 5 then $
   message, 'ENERG_HI tag in the XRAY_ARF_FILE file for SED_ID '+sed_id+' must be of type int, float or double.'
 if size(xray_arf.ENERG_HI, /n_dim) gt 1 then $
   message, 'ENERG_HI tag in the XRAY_ARF_FILE file for SED_ID '+sed_id+' must be a scalar or 1D array.'
 if total(size(xray_arf.ENERG_HI, /dim) ne size(xray_arf.SPECRESP, /dim)) gt 0 then $
   message, 'ENERG_HI tag in the XRAY_ARF_FILE file for SED_ID '+sed_id+' must have the same size as the SPECRESP tag.'
 if min(xray_arf.ENERG_HI) le 0 then $
   message, 'ENERG_HI tag in the XRAY_ARF_FILE file for SED_ID '+sed_id+' must only contain positive values.'



; Compile data
 Nxray = n_elements(xray_spec.ENERG_LO)
 xray_bandpass = dblarr(2, Nxray)
 xray_bandpass[0,*] = xray_spec.ENERG_LO
 xray_bandpass[1,*] = xray_spec.ENERG_HI

 xray_exposure = xray_spec.EXPOSURE
 net_counts = xray_spec.NET_COUNTS
 case strupcase(xray_unc) of
   'SQRT'    : net_counts_unc = sqrt(net_counts)
   'GEHRELS' : net_counts_unc = 1 + sqrt(0.75 + net_counts) ; The upper error from the Gehrels prescription
   'USER'    : net_counts_unc = xray_spec.NET_COUNTS_UNC
 endcase

 Nchannels = n_elements(xray_arf.ENERG_LO)
 arf_E_lo = xray_arf.ENERG_LO
 arf_E_hi = xray_arf.ENERG_HI
 arf_specresp = xray_arf.SPECRESP

; Pre-generating structure guarantees data type when input data is inserted
 xray_data = {XRAY_BANDPASS  : !values.D_NaN*dblarr(2, Nxray),  $
              XRAY_EXPOSURE  : !values.D_NaN*dblarr(Nxray),     $
              NET_COUNTS     : !values.D_NaN*dblarr(Nxray),     $
              NET_COUNTS_UNC : !values.D_NaN*dblarr(Nxray),     $
              ARF_E_LO       : !values.D_NaN*dblarr(Nchannels), $
              ARF_E_HI       : !values.D_NaN*dblarr(Nchannels), $
              ARF_SPECRESP   : !values.D_NaN*dblarr(Nchannels)  $
              }

 xray_data.XRAY_BANDPASS  = xray_bandpass
 xray_data.XRAY_EXPOSURE  = xray_exposure
 xray_data.NET_COUNTS     = net_counts
 xray_data.NET_COUNTS_UNC = net_counts_unc
 xray_data.ARF_E_LO       = arf_E_lo
 xray_data.ARF_E_HI       = arf_E_hi
 xray_data.ARF_SPECRESP   = arf_specresp

 return, xray_data

end
