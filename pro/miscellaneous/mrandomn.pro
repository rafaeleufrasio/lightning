function mrandomn, seed, covar, nrand=nrand, error_check=error_check, status=status
;+
; Name
; ----
;   MRANDOMN
;
; Purpose
; -------
;   Draws random deviates from a multivariate normal distribution
;   with zero mean and given covariance matrix.
;
; Calling Sequence
; ----------------
;   ::
;
;       rand_val = mrandomn(seed, covar [, nrand = , /error_check, status=status])
;
; Inputs
; ------
;   ``seed`` : any
;       The random number generator seed, the default is IDL's default in ``randomn()``.
;   ``covar`` : int, float, or double array(N, N)
;       The covariance matrix of the multivariate normal distribution.
;
; Optional Inputs
; ---------------
;   ``nrand`` : int, float, or double scalar
;       The number of randomn deviates to draw for each parameter, ``N``. (Default = ``1``)
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``rand_val`` : double array(nrand, N)
;       The random deviates.
;
; Optional Output
; ---------------
;   ``status`` : int scalar
;       Status of the Cholesky decomposition. If ``status = 0`` then the
;       computation was successful. If ``status > 0`` then the input covariance
;       matrix is not positive definite (see ``LA_CHOLDC.pro``), and ``rand_val = NaN``.
;       If not requested for output, then an error message will be printed
;       if ``status > 0``.
;
; Modification History
; --------------------
;   - 2004/09/01: Created (Brandon C. Kelly)
;   - 2013/10/01: Use ``LA_CHOLDC`` instead of ``CHOLDC`` to enable use of ``STATUS``
;     keyword. (W. Landsman)
;   - 2021/09/01: If ``status`` keyword is supplied then no error message will be printed,
;     and the function returns ``-1``. (Erik B. Monson)
;   - 2021/03/21: Made ``nrand`` an actual optional input (Keith Doore)
;   - 2021/03/21: Updated documentation (Keith Doore)
;   - 2021/03/21: Updated error handling (Keith Doore)
;   - 2021/04/07: Removed seed error handling and updated ``nrand`` error handling (Keith Doore)
;   - 2022/05/17: Added ``error_check`` keyword to do error handling (Keith Doore)
;   - 2022/06/24: Made ``rand_val = NaN`` if Cholesky decomposition failed (Keith Doore)
;   - 2022/06/24: Fixed issue of ``np`` creation being inside of ``error_check`` (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error Handling
 if keyword_set(error_check) then begin
   if n_elements(covar) eq 0 then message, 'Variable is undefined: COVAR.'
   if size(covar, /type) lt 2 or size(covar,/type) gt 5 then message, 'COVAR must be of type int, float, or double.'
   if size(covar, /n_dim) ne 2 then message, 'COVAR must be a 2-D array.'
   if (size(covar,/dim))[1] ne (size(covar,/dim))[0] then message, 'COVAR must be a square matrix.'

   if n_elements(nrand) ne 0 then begin
     if size(nrand, /type) lt 2 or size(nrand, /type) gt 5 then $
              message, 'NRAND must be of type int, float, or double.'
     if size(nrand, /n_dim) ne 0 then message, 'NRAND must be a scalar.'
     if nrand lt 1 then message, 'NRAND must be a positive value.'
   endif
 endif
 if n_elements(nrand) eq 0 then nrand = 1
 np = (size(covar,/dim))[0]


; Compute random deviates
 epsilon = randomn(seed, nrand, np) ;standard normal random deviates (NP x NRAND matrix)

 A = covar  ;store covariance into dummy variable for input into la_choldc

 la_choldc, A, /double, status=status        ;do Cholesky decomposition
 if status ne 0 then begin
   if ~arg_present(status) then message, 'COVAR is not positive definite, STATUS = ' + strtrim(status,2), /CON
   return, !values.D_NaN
 endif

 for i = 0, np - 2  do A[i+1:*,i] = 0d        ;Zero out upper triangular portion

;transform standard normal deviates so they have covariance matrix COVAR
 epsilon = A ## epsilon

 return, epsilon

end