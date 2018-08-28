;+
; NAME:
;	Trapex
; PURPOSE:
;	Integrate an exponentially varying function given values on
;	increasing grid of x-points, using trapeziodal rule in Log space.
;	More accurate than straight trapez. integration for such a case.
; CALLING:
;	fint = Trapex( f, x )
; INPUTS:
;	f = array of function values at points x, must be all positive.
;	x = array of independent variable values.
; KEYWORDS:
;	/VECTOR : let the upper limit of integral vary over the x-values
;		and return a vector of definite integrals.
;	/DELTA : return a vector of definite integrals over
;		the intervals between the x-values (i.e. do not accumulate).
;	LRATMIN = floating, default = 0.1.
;		If the absolute value of the Log of adjacent f(x) ratios
;		is less than LRATMIN then regular trapeziodal rule is used,
;		otherwise assume exponential type of variation.
; OUTPUTS:
;	Returns approximate integral(s) (exact or exponential function).
; PROCEDURE:
;	Use Log of function variation to from sum which is exact for exp. func.
; HISTORY:
;	Written: Frank Varosi, U.Md, 1988.
;-

function Trapex, f, x, LRATMIN=Lratmin, VECTOR=vector, DELTA_VECTOR=delta_vec

	Lrat = aLog( f )
	Lrat = Lrat(1:*) - temporary( Lrat )
	if N_elements( Lratmin ) NE 1 then Lratmin = 0.1

	df = make_array( SIZE=size( f ) )
	w = where( abs( Lrat ) LT Lratmin, nz )

	if (nz GT 0) then begin

		w1 = w+1
		df(w1) = 0.5 * ( f(w1) + f(w) ) * ( x(w1) - x(w) )
		w = where( abs( Lrat ) GE Lratmin, ng )

		if (ng GT 0) then begin
			w1 = w+1
			df(w1) = ( f(w1) - f(w) ) * ( x(w1) - x(w) ) / Lrat(w)
		   endif

	 endif else df(1) = ( f(1:*) - f ) * ( x(1:*) - x ) / Lrat

	if keyword_set( delta_vec ) then return, df

	if keyword_set( vector ) then begin

		for i=1L,N_elements( df )-1 do df(i) = df(i) + df(i-1)
		return, df

	 endif else return, total( df )
end
