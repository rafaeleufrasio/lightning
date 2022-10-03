function percentile, array, percentile
;+
; Name
; ----
;   PERCENTILE
;
; Purpose
; -------
;   Determines the percentiles values for an input array by sorting the
;   input array and returning an interpolated percentile value.
;   It assumes the lowest value of the array to be the 0th percentile
;   and the highest to be the 100th percentile.
;
; Calling Sequence
; ----------------
;   ::
;
;       percentile_value = percentile(array, percentile)
;
; Inputs
; ------
;   ``array`` : int, float, or double array(N)
;       The data that will have its percentile values determined.
;   ``percentile`` : int, float, or double scalar or array(M, ...)
;       The percentiles in decimal form (i.e., 50th percentile = ``0.5``).
;
; Output
; ------
;   ``percentile_value`` : double scalar or array(M, ...)
;       The percentile values corresponding to each input percentile.
;
; Modification History
; --------------------
;   - 2016/05/01: Created (Rafael T. Eufrasio)
;   - 2022/03/15: Added error handling (Keith Doore)
;   - 2022/03/15: Updated documentation (Keith Doore)
;   - 2022/04/07: Updated error handling to allow for degenerate dimensions (Keith Doore)
;   - 2022/06/24: Replaced ``.length`` with ``n_elements`` (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error Handling
 if n_elements(array) eq 0 then message, 'Variable is undefined: ARRAY.'
 if size(array, /type) lt 2 or size(array, /type) gt 5 then $
   message, 'ARRAY must be of type int, float, or double.'
 if size(reform(array), /n_dim) ne 1 then message, 'ARRAY must be a 1-D array.'
 if n_elements(array) lt 2 then message, 'ARRAY must have more than one element.'

 if n_elements(percentile) eq 0 then message, 'Variable is undefined: PERCENTILE.'
 if size(percentile, /type) lt 2 or size(percentile, /type) gt 5 then $
   message, 'PERCENTILE must be of type int, float, or double.'
 if min(percentile) lt 0 or max(percentile) gt 1 then message, 'PERCENTILE must contain values between 0 and 1.'


; Determine percentile values
 Nsample = n_elements(array)
 percentile_grid = dindgen(Nsample)/(Nsample - 1)

 return, interpol(array[sort(array)], percentile_grid, percentile)

end