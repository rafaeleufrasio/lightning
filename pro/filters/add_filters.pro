pro add_filters, filter_profile, bandpass, instrument_survey, observatory=observatory, $
                 unit=unit, frequency=frequency, _extra=_extra
;+
; Name
; ----
;   ADD_FILTERS
;
; Purpose
; -------
;   Adds a user provided filter profile to Lightning. The filter profile can
;   either be input as a two column text file or array. It is then formatted
;   for use in Lightning. Additionally, the bandpass name, instrument and
;   observatory, or survey will need to be specified.
;
; Calling Sequence
; ----------------
;   ::
;
;       add_filters, filter_profile, bandpass, instrument_survey [, observatory = , $
;                    unit = , /frequency, _extra=_extra]
;
; Inputs
; ------
;   ``filter_profile`` : int, float, or double array(2, Nwave) or string scalar
;       If ``filter_profile`` is a string, then it contains the path and file name
;       to the tabulated filter profile that is to be read. For both the tabulated
;       filter profile and the directly input array, the first column must contain
;       the grid of wavelengths (or frequency, see below), and the second column
;       must contain the corresponding transmission profile.
;   ``bandpass`` : string scalar
;       The bandpass name to give the filter profile.
;   ``instrument_survey`` : string scalar
;       The instrument or survey associated with the filter profile.
;
; Optional Inputs
; ---------------
;   ``observatory`` : string scalar
;       The observatory associated with the filter profile. Omit if
;       the filter is intended to be associated with a survey. (See
;       Note below.)
;   ``unit`` : int, float, or double scalar
;       The unit conversion factor needed to convert the input wavelengths
;       (frequency) associated with the filter profile to microns (Hertz).
;       For example, if the input wavelengths are in nanometers, ``unit``
;       should be set to ``1d-3``. (Default = ``1``)
;   ``frequency`` : flag
;       If set, the input or read-in filter profile is gridded in frequency
;       rather than wavelength.
;   ``_extra`` : structure
;       Additional optional inputs that are passed to ``readcol.pro`` from the
;       NASA Astro library (https://idlastro.gsfc.nasa.gov/ftp/pro/misc/readcol.pro).
;
; Output
; ------
;   A text file containing the formatted filter profile that is saved to the users
;   local Lightning Filters directory. Additionally, the Lightning filter label
;   is printed to the screen.
;
; Notes
; -----
;   The output file and its path within the Lightning Filters directory is named
;   using the scheme: ``Filters/observatory/instrument/instrument_band.txt``. In the
;   case of survey specific filters (e.g., 2MASS, SDSS, etc.), the pattern is
;   changed to ``Filters/survey/survey_band.txt``. 
;
; Modification History
; --------------------
;   - 2022/07/28: Created (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(filter_profile) eq 0 then message, 'Variable is undefined: FILTER_PROFILE.'
 if size(filter_profile, /type) eq 7 then begin
   if size(filter_profile, /n_dim) ne 0 then message, 'FILTER_PROFILE must be a scalar if type string.'
   if ~file_test(filter_profile) then message, 'FILTER_PROFILE does not lead to an existing file. '+$
                                               'Please check that file exists.'
   file = 1
 endif else if (size(filter_profile, /type) ge 2 and size(filter_profile, /type) le 5) then begin
   if size(filter_profile, /n_dim) ne 2 then message, $
     'FILTER_PROFILE must be a 2-D array if type int, float, or double.'
   if min(filter_profile[0, *]) lt 0 then message, 'FILTER_PROFILE wavelength must only contain positive values.'
 endif else message, 'FILTER_PROFILE must be of type int, float, double, or string.'
 
 if n_elements(bandpass) eq 0 then message, 'Variable is undefined: BANDPASS.'
 if size(bandpass, /type) ne 7 then message, 'BANDPASS must be of type string.'
 if size(bandpass, /n_dim) ne 0 then message, 'BANDPASS must be a scalar.'

 if n_elements(instrument_survey) eq 0 then message, 'Variable is undefined: INSTRUMENT_SURVEY.'
 if size(instrument_survey, /type) ne 7 then message, 'INSTRUMENT_SURVEY must be of type string.'
 if size(instrument_survey, /n_dim) ne 0 then message, 'INSTRUMENT_SURVEY must be a scalar.'

 if n_elements(unit) ne 0 then begin
   if size(unit, /type) lt 2 or size(unit, /type) gt 5 then message, 'UNIT must be of type int, float, or double.'
   if size(unit, /n_dim) ne 0 then message, 'UNIT must be a scalar.'
 endif else unit=1.d0

 if n_elements(observatory) ne 0 then begin
   if size(observatory, /type) ne 7 then message, 'OBSERVATORY must be of type string.'
   if size(observatory, /n_dim) ne 0 then message, 'OBSERVATORY must be a scalar.'
   folder = observatory+'/'+instrument_survey+'/'
 endif else folder = instrument_survey+'/'


; Read-in or separate wavelength and transmission
 if keyword_set(file) then begin
   readcol, filter_profile, wave, transmission, _extra=_extra
 endif else begin
   wave = reform(filter_profile[0, *])
   transmission = reform(filter_profile[1, *])
 endelse

; Compile constants in case needed for frequency conversion
 lightning_constants

; Convert units  
 wave = wave * unit 
 if keyword_set(frequency) then wave = 1.d4*!lightning_cgs.clight/wave

; Normalize and sort
 transmission /= max(transmission)
 transmission = transmission[sort(wave)]
 transmission[where(transmission le 0, /null)] = 0
 wave = wave[sort(wave)]

; Print to file
 wave=string(wave, f='(E13.7)')
 transmission = string(transmission, f='(E13.7)')
 wave_trans_str='  '+wave+'    '+transmission

 file_mkdir, !lightning_dir+'filters/'+folder
 openw, tab, !lightning_dir+'filters/'+folder+instrument_survey+'_'+bandpass+'.txt', /get_lun

 printf, tab, '# wave[microns]    norm_trans'
 for j=0, (n_elements(wave)-1) do printf, tab, wave_trans_str[j]
 free_lun, tab
 close

 print_to_width, "Filter profile successfully added to Lightning. "+$
                 "Its Lightning filter label is: '"+instrument_survey+"_"+bandpass+"'."

; Update filters.json
 ; Extract all filter names and corresponding path
 filter_full_path = file_search(!lightning_dir+'filters', '*.txt')
 Nfilters_json = n_elements(filter_full_path)
 filter_path = strmid(filter_full_path, strlen(!lightning_dir+'filters/'))
 filter_name = file_basename(filter_path, '.txt')
 max_len = max(strlen(filter_name))

 openw, json, !lightning_dir+'filters/filters.json', /get_lun

 printf, json, '{'
 for i=0, Nfilters_json-1 do begin
   if i eq Nfilters_json-1 then comma = '' else comma = ','
   printf, json, '  "'+filter_name[i]+'":'+strjoin(replicate(' ', max_len-strlen(filter_name[i])+2))+'"'+filter_path[i]+'"'+comma
 endfor
 printf, json, '}'
 free_lun, json
 close

end