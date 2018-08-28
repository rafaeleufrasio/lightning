;+
; NAME:
;	Trapez
; PURPOSE:
;	Integrate  real function given values on increasing grid
;	using the trapeziodal rule.
;	To integrate rapidly varying functions try using function trapex.
; CALLING:
;	fint = Trapez( f, x )
; INPUTS:
;	f = array of function values at points x.
;	x = array of independent variable values.
; KEYWORDS:
;	/VECTOR : let the upper limit of integral vary over the x-values
;		and return a vector of definite integrals.
;	/DELTA : return a vector of definite integrals over
;		the intervals between the x-values (i.e. do not accumulate).
; OUTPUTS:
;	Function returns approximate integral(s).
; HISTORY:
;	Written: Frank Varosi, U.Md, 1988.
;-

function Trapez, f, x, VECTOR=vector, DELTA_VECTOR=delta_vec

	df = make_array( SIZE=size( f ) )
	df[1] = ( f[1:*] + f ) * ( x[1:*] - x )

	if keyword_set( delta_vec ) then return, 0.5 * df

	if keyword_set( vector ) then begin

		for i=1L,N_elements( df )-1 do df(i) = df(i) + df(i-1)
		return, 0.5 * df

	 endif else return, 0.5 * total( df )
end
