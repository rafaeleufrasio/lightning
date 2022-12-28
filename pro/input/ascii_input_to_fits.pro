pro ascii_input_to_fits, input_file, input_file_fits=input_file_fits
;+
; Name
; ----
;   ASCII_INPUT_TO_FITS
;
; Purpose
; -------
;   Reads in the SED data that is to be fit by Lightning from an ASCII file. The 
;   read-in data is then placed in a formatted structure and saved to the first
;   extension of a FITS file with the same name as the input file.
;
; Calling Sequence
; ----------------
;   ::
;
;       ascii_input_to_fits, input_file [, input_file_fits=input_file_fits]
;
; Input
; -----
;   ``input_file`` : string scalar
;       The name (including path) to the ASCII file containing the SED fluxes 
;       and distances (or redshifts) in a data table. (See :ref:`ascii-format-label`
;       for full details, required contents, and format of the ASCII data table.)
;
; Output
; ------
;   A FITS file (``<input_file_dir>/<input_file_name>.fits``) containing the
;   SED fluxes and distances (or redshifts) as read in from the ASCII data
;   table. (See :ref:`fits-format-label` for full details and format of the
;   FITS data table.)
;
; Optional Output
; ---------------
;   ``input_file_fits`` : string scalar
;       The name (including path) of the output FITS file.
;
; Note
; ----
;   No error handling is performed on the data in the table as this is done
;   with the converted FITS file data in ``lightning_input.pro``.
;
; Modification History
; --------------------
;   - 2022/08/02: Created (Keith Doore)
;   - 2022/09/19: Allowed for X-ray fluxes to be input (Keith Doore)
;   - 2022/12/28: Fixed bug if only reading in one SED (Keith Doore)
;   - 2022/12/28: Fixed bug if no x-ray data (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(input_file) eq 0 then message, 'Variable is undefined: INPUT_FILE.'
 if size(input_file, /type) ne 7 then message, 'INPUT_FILE must be of type string.'
 if size(input_file, /n_dim) ne 0 then message, 'INPUT_FILE must be a scalar.'
 if ~file_test(input_file) then message, 'No input data file exists. '+$
                                        'Please confirm specified file name leads to the correct file.'


; Read in ASCII data file
; Determine the number of lines in the file.
 Nlines = file_lines(input_file)
 line_start = lonarr(Nlines + 1)

; Determine the line where the header ends.
 openr, lun, input_file, /Get_Lun
 in_header = 1
 current_line = 0L
 while in_header do begin
   line = ''
   readf, lun, line
   ; Check for delimiters of space, tab, or comma.
   ; Extract comment character even if no space between comment character and text
   comment_character = strmid((strsplit(line, ' '+string(9b)+',', /extract))[0], 0, 1)
   if comment_character ne '#' then in_header = 0
   point_lun, -lun, currentloc
   ; Add one to current_line index so, line_start[0] = 0 (i.e., beginning of file)
   line_start[current_line + 1] = currentloc
   current_line++
 endwhile

; Create a variable to hold the data and one to hold data column tags.
 ; Subtract one, since current line includes 1st data line.
 data = strarr(Nlines - (current_line - 1))
 tags = ''

; Reset the data file to the last line of header to read column tags, then read in data
 ; Subtract 2, since current line index is start of 2nd data line.
 point_lun, lun, line_start[current_line - 2]
 readf, lun, tags
 readf, lun, data
 free_lun, lun


; Split strings into each column name and data columns
 ; Remove potential comment character tag if a space between it and SED_ID column name
 ; Extract comment character again as it was overridden by first character of 1st data line
 tags = strsplit(tags, ' '+string(9b)+',', /extract)
 comment_character = strmid(tags[0], 0, 1)
 if tags[0] eq comment_character then tags = tags[1:*]

 data = strsplit(data, ' '+string(9b)+',', /extract)
 ; If data is a list (case of more than one data line), convert to array
 if size(data, /type) eq 11 then begin
   data = data.ToArray()
   ; Tranpose to have [Ncol, Nsed] (also makes indexing simple if Nsed = 1)
   data = transpose(data)
   Nsed = (size(data, /dim))[1]
 endif else Nsed = 1


; Extract data from each column using column names
 ; Check for SED ID
 id_idc = where(strupcase(tags) eq 'SED_ID', /null)
 if n_elements(id_idc) ne 0 then sed_id = reform(data[id_idc, *])

 ; Check for luminosity distance and redshift
 dist_idc = where(strupcase(tags) eq 'LUMIN_DIST', Ndist, /null)
 red_idc = where(strupcase(tags) eq 'REDSHIFT', Nred, /null)
 if Ndist ne 0 then lumin_dist = double(reform(data[dist_idc, *]))
 if Nred ne 0 then redshift = double(reform(data[red_idc, *]))

 ; Check for Xray-data
 ; GALACTIC_NH
 gal_nh_idc  = where(strupcase(tags) eq 'GALACTIC_NH', Ngal_nh, /null)
 if Ngal_nh ne 0 then galactic_nh = double(reform(data[gal_nh_idc, *]))

 ; XRAY_SPEC_FILE and XRAY_ARF_FILE
 xray_spec_idc = where(strupcase(tags) eq 'XRAY_SPEC_FILE', Nxray_spec, /null)
 xray_arf_idc  = where(strupcase(tags) eq 'XRAY_ARF_FILE', Nxray_arf, /null)
 if Nxray_spec ne 0 then xray_spec_file = reform(data[xray_spec_idc, *])
 if Nxray_arf ne 0 then xray_arf_file = reform(data[xray_arf_idc, *])

 ; XRAY_BANDPASS, XRAY_FLUX, and XRAY_FLUX_UNC
 xray_band_l_idc   = where(strmatch(strupcase(tags), 'XRAY_BANDPASS_L_*'), Nxray_band_l, /null)
 xray_band_u_idc   = where(strmatch(strupcase(tags), 'XRAY_BANDPASS_U_*'), Nxray_band_u, /null)
 xray_flux_idc     = where(strmatch(strupcase(tags), 'XRAY_FLUX_[!UNC]*'), Nxray_flux, /null)
 xray_flux_unc_idc = where(strmatch(strupcase(tags), 'XRAY_FLUX_UNC_*'), Nxray_flux_unc, /null)
 if Nxray_band_l ne 0 then xray_band_l = reform(data[xray_band_l_idc, *])
 if Nxray_band_u ne 0 then xray_band_u = reform(data[xray_band_u_idc, *])
 if Nxray_flux ne 0 then xray_flux = reform(data[xray_flux_idc, *])
 if Nxray_flux_unc ne 0 then xray_flux_unc = reform(data[xray_flux_unc_idc, *])
 ; Combine xray bandpass lower and upper bounds into one array
 Nxray_band = max([Nxray_band_l, Nxray_band_u])
 if Nxray_band gt 0 then begin
   xray_bandpass = !values.D_NaN*dblarr(2, Nxray_band, Nsed)
   xray_bandpass[0, 0:Nxray_band_l-1, *] = xray_band_l
   xray_bandpass[1, 0:Nxray_band_u-1, *] = xray_band_u
 endif

 ; Extract the filter indices
 nonfilter_idc = where(strupcase(tags) eq 'SED_ID'        or strupcase(tags) eq 'LUMIN_DIST'  or $
                       strupcase(tags) eq 'REDSHIFT'      or strupcase(tags) eq 'GALACTIC_NH' or $
                       strmatch(strupcase(tags), 'XRAY_*'), comp=filters_idc, ncomp=Nfilters, /null)
 filter_tags = tags[filters_idc]
 ; Sort filter indices so that they have order of label1 -> label1_err -> label2 -> label2_err -> ...
 filters_idc = filters_idc[sort(filter_tags)]
 ; Now only take every other index starting with the first to have label1 -> label2 -> ...
 ;   and every other starting with the second to have label1_err -> label2_err -> ...
 obs_idc = filters_idc[0:Nfilters-1:2]
 unc_idc = filters_idc[1:Nfilters-1:2]

 ; Divide Nfilters by two as it includes both observation and uncertainty
 Nfilters /= 2
 filter_labels = tags[obs_idc]
 Fnu_obs = double(data[obs_idc, *])
 Fnu_unc = double(data[unc_idc, *])


; Place extracted data into structure
 ; Make dynamically with hash, then convert to structure
 out = orderedhash()
 if n_elements(sed_id) ne 0 then out['SED_ID'] = ''

 out['FNU_OBS'] = !values.D_NaN*dblarr(Nfilters)
 out['FNU_UNC'] = dblarr(Nfilters)
 out['FILTER_LABELS'] = strarr(Nfilters)

 if Ndist ne 0 then out['LUMIN_DIST'] = !values.D_NaN
 if Nred ne 0 then out['REDSHIFT'] = !values.D_NaN

 if Ngal_nh ne 0 then out['GALACTIC_NH'] = !values.D_NaN

 if Nxray_spec ne 0 then out['XRAY_SPEC_FILE'] = ''
 if Nxray_arf ne 0 then out['XRAY_ARF_FILE'] = ''

 if Nxray_band ne 0 then out['XRAY_BANDPASS'] = !values.D_NaN*dblarr(2, Nxray_band)
 if Nxray_flux ne 0 then out['XRAY_FLUX'] = !values.D_NaN*dblarr(Nxray_flux)
 if Nxray_flux_unc ne 0 then out['XRAY_FLUX_UNC'] = !values.D_NaN*dblarr(Nxray_flux_unc)


 out = out.ToStruct(/recursive)
 out = replicate(out, Nsed)

 out_tags = tag_names(out)

 if total(strupcase(out_tags) eq 'SED_ID') eq 1 then out.SED_ID = sed_id

 out.FNU_OBS = Fnu_obs
 out.FNU_UNC = Fnu_unc
 for i=0, Nsed-1 do out[i].FILTER_LABELS = filter_labels

 if total(strupcase(out_tags) eq 'LUMIN_DIST') eq 1 then out.LUMIN_DIST = lumin_dist
 if total(strupcase(out_tags) eq 'REDSHIFT') eq 1 then out.REDSHIFT = redshift

 if total(strupcase(out_tags) eq 'GALACTIC_NH') eq 1 then out.GALACTIC_NH = galactic_nh

 if total(strupcase(out_tags) eq 'XRAY_SPEC_FILE') eq 1 then out.XRAY_SPEC_FILE = xray_spec_file
 if total(strupcase(out_tags) eq 'XRAY_ARF_FILE') eq 1 then out.XRAY_ARF_FILE = xray_arf_file

 if total(strupcase(out_tags) eq 'XRAY_BANDPASS') eq 1 then out.XRAY_BANDPASS = xray_bandpass
 if total(strupcase(out_tags) eq 'XRAY_FLUX') eq 1 then out.XRAY_FLUX = xray_flux
 if total(strupcase(out_tags) eq 'XRAY_FLUX_UNC') eq 1 then out.XRAY_FLUX_UNC = xray_flux_unc

; Save structure to FITS file
 input_file_fits = file_dirname(input_file , /mark_dir)+file_basename(input_file, '.txt')+'.fits'
 mwrfits, out, input_file_fits, /create

end