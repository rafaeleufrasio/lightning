;+
; NAME:
;	Trapow
; PURPOSE:
;	Integrate a function having power-law variation, given values on
;	increasing grid of x-points, using trapeziodal rule in Log-Log space.
;	More accurate than straight trapez. integration (exact for power-law).
; CALLING:
;	fint = Trapow( f, x )
; INPUTS:
;	f = array of function values at points x, must be all positive.
;	x = array of independent variable values,
;		must be all positive and strictly increasing.
; KEYWORDS:
;	/VECTOR : let the upper limit of integral vary over the x-values
;		and return a vector of definite integrals.
;	/DELTA : return a vector of definite integrals over
;		the intervals between the x-values (i.e. do not accumulate).
; OUTPUTS:
;	Returns approximate integral(s) (exact for power-law).
; PROCEDURE:
;	Assume power-law variation and cast in Log-Log space.
;	If power variation is close to -1, use Log integ. case.
; HISTORY:
;	Written: Frank Varosi, HSTX @ NASA/GSFC, 1995.
; MODIFICATION:
;      i was changed to be a LONG integer in the loop 
;      by Eli Dwek , March 7, 2001 
;-

function Trapow, f, x, VECTOR=vector, DELTA_VECTOR=delta_vec

	f = double( f )
	x = double( x )
	tmpf = aLog( f )
	pow = tmpf(1:*) - tmpf
	tmpf = x * f
	Logx = aLog( x )
	Logx = Logx(1:*) - Logx
	pow = 1 + temporary( pow )/Logx
	if N_elements( zero ) NE 1 then zero = 1.e-7

	df = make_array( SIZE=size( f ) )
	w = where( abs( pow ) LT zero, nz )

	if (nz GT 0) then begin

		w1 = w+1
		df(w1) = 0.5 * ( tmpf(w1) + tmpf(w) ) * Logx(w)
		w = where( abs( pow ) GE zero, nw )

		if (nw GT 0) then begin
			w1 = w+1
			df(w1) = ( tmpf(w1) - tmpf(w) )/ pow(w)
		   endif

	 endif else df(1) = ( tmpf(1:*) - tmpf ) / pow

	if keyword_set( delta_vec ) then return, df

	if keyword_set( vector ) then begin

		for i=1L,N_elements( df )-1 do df(i) = df(i) + df(i-1)
		return, df

	 endif else return, total( df )
end
