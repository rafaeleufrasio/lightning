function trap_int, x, f, xrange=xrange, log_axes=log_axes, exponential=exponential, $
                   power_law=power_law, vector=vector, delta_vector=delta_vector, $
                   error_check=error_check
;+
; Name
; ----
;   TRAP_INT
;
; Purpose
; -------
;   Integrates a tabulated function using the trapeziodal rule, with options
;   to increase accuracy in the case of functions having exponential
;   or power-law variation by recasting trapeziodal rule in linear-log or
;   log-log (x,y) space, respectively. A user can also specify the limits
;   of integration, and the function will be automatically extrapolated
;   if requested range exceeds the tabulated range of function.
;
; Calling Sequence
; ----------------
;   ::
;
;      fint = trap_int(x, f [, xrange = , log_axes = , /exponential, $
;                      /power_law, /vector, /delta_vector, /error_check])
;
; Inputs
; ------
;   ``x`` : int, float, or double array(Nf)
;       Independent variable values.
;   ``f`` : int, float, or double array(Nf)
;       Function values at points ``x``.
;
; Optional Inputs
; ---------------
;   ``xrange`` : int, float, or double scalar or array(2)
;       Specifies the limits of integration. If just one number is given, it is
;       assumed to be the lower limit. If requested range exceeds tabulated range
;       of function, the values are extrapolated. (Default = ``minmax(x)``)
;   ``log_axes`` : int, float, or double scalar
;       Alternative indicator of integration method:
;       ``0`` = Linear-Linear space (Default), 
;       ``1`` = Linear-Log (equivalent to using ``/exponential``), 
;       ``2`` = Log-Log (equivalent to using ``/power_law``).
;   ``exponential`` : flag
;       If set, integrate (and interpolate) in linear-log space, which increases
;       accuracy for exponentially varying functions. Takes precedence over 
;       ``log_axes``. In this case, ``f`` array is thresholded > ``1.d-200``.
;   ``power_law`` : flag
;       If set, integrate (and interpolate) in log-log space, which increases
;       accuracy for functions with power-law variation. Takes precedence over
;       ``log_axes``. In this case, both ``x`` and ``f`` arrays are 
;       thresholded > ``1.d-200``.
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
;   ``fint`` : double scalar or array(Nf)
;       Approximate integral(s). The size of ``fint`` is determined if one of
;       the optional inputs, ``vector`` or ``delta_vector``, is set. If one is
;       set, then ``fint`` is an array. Otherwise, it is a scalar.
;
; Modification History
; --------------------
;   - 1988/01/01: Created (Frank Varosi)
;   - 2020/02/18: Truncation value changed from ``1.e-37`` to ``1.d-200`` (Rafael Eufrasio)
;   - 2021/03/21: Standardized parameter names (Keith Doore)
;   - 2021/03/21: Updated documentation (Keith Doore)
;   - 2021/03/21: Corrected indexing to brackets from parentheses (Keith Doore)
;   - 2021/03/21: Added error handling (Keith Doore)
;   - 2022/04/12: Allowed for inputs to have degenerate dimensions (Keith Doore)
;   - 2022/04/12: Allowed integer inputs (Keith Doore)
;   - 2022/05/17: Added ``error_check`` keyword to do error handling (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error Handling
 if keyword_set(error_check) then begin
   if n_elements(x) eq 0 then message, 'Variable is undefined: X.'
   if size(x, /type) lt 2 or size(x, /type) gt 5 then message, 'X must be of type int, float, or double.'
   if size(reform(x), /n_dim) ne 1 then message, 'X must be a 1-D array.'
   if n_elements(x) lt 2 then message, 'X must have more than one element.'
  
   if n_elements(f) eq 0 then message, 'Variable is undefined: F.'
   if size(f, /type) lt 2 or size(f, /type) gt 5 then message, 'F must be of type int, float, or double.'
   if size(reform(f), /n_dim) ne 1 then message, 'F must be a 1-D array.'
   if n_elements(f) ne n_elements(x) then message, 'F must have the same number of elements as X.'
  
   if n_elements(xrange) ne 0 then begin
     if size(xrange, /type) lt 2 or size(xrange,/type) gt 5 then $
              message, 'XRANGE must be of type int, float, or double.'
     if size(reform(xrange), /n_dim) gt 1 then message, 'XRANGE must be a scalar or 1D array.'
     if n_elements(xrange) lt 1 or n_elements(xrange) gt 2 then message, 'XRANGE must have one or two elements.'
   endif
  
   if keyword_set(vector) and keyword_set(delta_vector) then $
     message, 'Both VECTOR and DELTA_VECTOR cannot be set at the same time.'
  
   if keyword_set(exponential) and keyword_set(power_law) then $
     message, 'Both EXPONENTIAL and POWER_LAW cannot be set at the same time.'
  
   if n_elements(log_axes) ne 0 then begin
     if size(log_axes, /type) lt 2 or size(log_axes,/type) gt 5 then $
              message, 'LOG_AXES must be of type int, float, or double.'
     if size(log_axes, /n_dim) ne 0 then message, 'LOG_AXES must be a scalar.'
     if total(log_axes eq [0, 1, 2]) ne 1 then message,'LOG_AXES must only contain values of 0, 1, or 2.'
   endif
 endif

 if n_elements(log_axes) eq 0 then log_axes = 0
 if keyword_set(exponential) then log_axes = 1
 if keyword_set(power_law) then log_axes = 2


; Integral calculation
; Check for specified integration Limits:

 if N_elements(xrange) GT 0 then begin ;setup interpolation of f(xrange):

   maxx = max(x, MIN=minx)
   xmin = xrange[0]
   if N_elements(xrange) GT 1 then xmax = xrange[1] else xmax = maxx
   xmin = xmin < xmax
   xmax = xmax > xmin

   if (xmin GT maxx) OR (xmax LT minx) then begin
     message,"Cannot integrate if XRANGE is outside of definition"
   endif else if (xmin EQ xmax) then begin
      message,"Specified XRANGE has zero length"
   endif

   w = where((x GE xmin) AND (x LE xmax), nw)

   if (nw LE 0) then begin   ;xrange is between x-grid points.

     wmin = where(x LE xmin, nmin) > 0
     wmax = where(x GE xmax) > 0
     wmin = [wmin[nmin-1], wmax[0]]
     wmax = wmin

   endif else begin

     if (xmin LE minx) then begin
       if (xmin LT minx) then $
         message, /INFO, "extrapolating Xmin"
       wmin = [0,1]
     endif else wmin = w[0] - [1,0]

     if (xmax GE maxx) then begin
       if (xmax GT maxx) then $
         message, /INFO, "extrapolating Xmax"
       wmax = N_elements(x) - [2,1]
     endif else wmax = w[nw-1] + [0,1]

   endelse

; interpolate f(xmin) & f(xmax) (result -> fi) depending on log_axes type:

   wi = transpose([[wmin], [wmax]])
   xi = [xmin, xmax]
   if (log_axes GT 0) then fw = aLog(f[wi] > 1d-200) else fw = f[wi]

   if (log_axes EQ 2) then begin
     xw = aLog(x[wi] > 1d-200)
     xi = aLog(xi > 1d-200)
   endif else xw = x[wi]

   fi = fw[*,0] + (xi-xw[*,0])*(fw[*,1]-fw[*,0])/(xw[*,1]-xw[*,0])
   if (log_axes GT 0) then fi = exp(fi)

   if (nw GT 0) then begin   ;check for duplicates:

     if (xmin LT x[w[0]]) then begin
       if (xmax GT x[w[nw-1]]) then begin
         fi = [fi[0], f[w], fi[1]]
         xi = [xmin, x[w], xmax]
       endif else begin
         fi = [fi[0], f[w]]
         xi = [xmin, x[w]]
       endelse
     endif else if (xmax GT x[w[nw-1]]) then begin
       fi = [f[w], fi[1]]
       xi = [x[w], xmax]
     endif else begin
       fi = f[w]
       xi = x[w]
     endelse

   endif else xi = [xmin, xmax]

   return, trap_int(xi, fi, LOG=log_axes, VEC=vector, DEL=delta_vector)
 endif

;Perform the integration:

 CASE log_axes OF

   2: BEGIN
        if min(x) LE 0 then begin
          message, /INFO, $
            "Restricting integral to positive values of X"
          w = where(x GE 1d-200, nw)
          if (nw LT 2) then return, 0 $
          else return, Trapow(f[w]>1d-200, x[w], $
                              VEC=vector, DEL=delta_vector)
        endif

        return, Trapow(f > 1d-200, x, VEC=vector, DEL=delta_vector)
      END

   1: BEGIN
        return, Trapex(f > 1d-200, x, VEC=vector, DEL=delta_vector)
      END

   else: BEGIN
           return, Trapez(f, x, VEC=vector, DEL=delta_vector)
         END

 ENDCASE
end
