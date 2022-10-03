function trapez, f, x, vector=vector, delta_vector=delta_vector, error_check=error_check
;+
; Name
; ----
;   TRAPEZ
;
; Purpose
; -------
;   Integrates a real function given values on increasing grid, using
;   the trapezoidal rule.
;
; Calling Sequence
; ----------------
;   ::
;
;       fint = trapez(f, x [, /vector, /delta_vector, /error_check])
;
; Inputs
; ------
;   ``f`` : int, float, or double array(Nf)
;       Function values at points ``x``.
;   ``x`` : int, float, or double array(Nf)
;       Independent variable values.
;
; Optional Inputs
; ---------------
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
;   - 2021/03/21: Updated documentation (Keith Doore)
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
  
   if n_elements(x) eq 0 then message, 'Variable is undefined: X.'
   if size(x, /type) lt 2 or size(x, /type) gt 5 then message, 'X must be of type int, float, or double.'
   if size(reform(x), /n_dim) ne 1 then message, 'X must be a 1-D array.'
   if n_elements(x) ne n_elements(f) then message, 'X must have the same number of elements as F.'
  
   if keyword_set(vector) and keyword_set(delta_vector) then $
     message, 'Both VECTOR and DELTA_VECTOR cannot be set at the same time.'
 endif

; integral calculation
 df = make_array(SIZE=size(reform(f)))
 df[1] = (f[1:*] + f) * (x[1:*] - x)

 if keyword_set(delta_vector) then return, 0.5 * df

 if keyword_set(vector) then begin

   for i=1L, n_elements(df)-1 do df[i] = df[i] + df[i-1]
   return, 0.5 * df

 endif else return, 0.5 * total(df)
end
