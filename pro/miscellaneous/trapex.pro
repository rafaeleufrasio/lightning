function trapex, f, x, lratmin=lratmin, vector=vector, delta_vector=delta_vector, $
                 error_check=error_check
;+
; Name
; ----
;   TRAPEX
;
; Purpose
; -------
;   Integrates an exponentially varying function given values on
;   increasing grid of x-points, using the trapeziodal rule in log space.
;   More accurate than linear trapezoidal integration for such a case.
;
; Calling Sequence
; ----------------
;   ::
;
;       fint = trapex(f, x [, lratmin = , /vector, /delta_vector, /error_check])
;
; Inputs
; ------
;   ``f`` : int, float, or double array(Nf)
;       Function values at points ``x``, must be all positive.
;   ``x`` : int, float, or double array(Nf)
;       Independent variable values.
;
; Optional Inputs
; ---------------
;   ``lratmin`` : int, float, or double scalar
;       If the absolute value of the log of adjacent ``f`` ratios
;       is less than ``lratmin`` then regular trapeziodal rule is used,
;       otherwise assume exponential type of variation. (Default = ``0.1``)
;   ``vector`` : flag
;       If set, returns a vector of cumulative definite integrals.
;   ``delta_vector`` : flag
;       If set, returns a vector of non-cumulative definite integrals
;       over the intervals between the ``x`` values.
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``fint`` : float or double scalar or array(Nf)
;       Approximate integral(s). The size of ``fint`` is determined if one of
;       the optional inputs, ``vector`` or ``delta_vector``, is set. If one is
;       set, then ``fint`` is an array. Otherwise, it is a scalar.
;
; Modification History
; --------------------
;   - 1988/01/01: Created (Frank Varosi)
;   - 2021/03/21: Standardized parameter names (Keith Doore)
;   - 2021/03/21: Updated documentation (Keith Doore)
;   - 2021/03/21: Corrected indexing to brackets from parentheses (Keith Doore)
;   - 2021/03/21: Added error handling (Keith Doore)
;   - 2022/04/08: Allowed for inputs to have degenerate dimensions (Keith Doore)
;   - 2022/04/08: Allowed integer inputs (Keith Doore)
;   - 2022/05/17: Added ``error_check`` keyword to do error handling (Keith Doore)
;-
 On_error, 2
 compile_opt IDL2

; Error Handling
 if keyword_set(error_check) then begin
   if n_elements(f) eq 0 then message, 'Variable is undefined: F.'
   if size(f, /type) lt 2 or size(f, /type) gt 5 then message, 'F must be of type int, float, or double.'
   if size(reform(f), /n_dim) ne 1 then message, 'F must be a 1-D array.'
   if n_elements(f) lt 2 then message, 'F must have more than one element.'
   if min(f) lt 0 then message, 'F must only contain positive values.'
  
   if n_elements(x) eq 0 then message, 'Variable is undefined: X.'
   if size(x, /type) lt 2 or size(x, /type) gt 5 then message, 'X must be of type int, float, or double.'
   if size(reform(x), /n_dim) ne 1 then message, 'X must be a 1-D array.'
   if n_elements(x) ne n_elements(f) then message, 'X must have the same number of elements as F.'
  
   if n_elements(Lratmin) ne 0 then begin
     if size(Lratmin, /type) lt 2 or size(Lratmin, /type) gt 5 then $
              message, 'LRATMIN must be of type int, float, or double.'
     if size(Lratmin, /n_dim) ne 0 then message, 'LRATMIN must be a scalar.'
     if Lratmin le 0 then message, 'LRATMIN must be a be a positive value.'
   endif
  
   if keyword_set(vector) and keyword_set(delta_vector) then $
     message, 'Both VECTOR and DELTA_VECTOR cannot be set at the same time.'
 endif

 if n_elements(Lratmin) eq 0 then Lratmin = 0.1

; integral calculation
 Lrat = aLog(f)
 Lrat = Lrat[1:*] - temporary(Lrat)

 df = make_array(SIZE=size(reform(f)))
 w = where(abs(Lrat) LT Lratmin, nz)

 if (nz GT 0) then begin

   w1 = w+1
   df[w1] = 0.5 * (f[w1] + f[w]) * (x[w1] - x[w])
   w = where(abs(Lrat) GE Lratmin, ng)

   if (ng GT 0) then begin
     w1 = w+1
     df[w1] = (f[w1] - f[w]) * (x[w1] - x[w]) / Lrat[w]
   endif

 endif else df[1] = (f[1:*] - f) * (x[1:*] - x) / Lrat

 if keyword_set(delta_vector) then return, df

 if keyword_set(vector) then begin

   for i=1L, n_elements(df)-1 do df[i] = df[i] + df[i-1]
   return, df

 endif else return, total(df)
end
