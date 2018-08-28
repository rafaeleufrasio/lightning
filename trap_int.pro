;+
; NAME:
;	Trap_Int
; PURPOSE:
;	Integrate a tabulated function using the trapeziodal rule, with options
;	to increase accuracy in the case of functions having exponential
;	or power-law variation by recasting trapeziodal rule in Linear-Log or
;	Log-Log (x,y) space, respectively. User can also specify the limits
;	of integration, and the function will be automatically extrapolated
;	if requested range exceeds the tabulated range of function. 
; CALLING:
;	fint = Trap_Int( x, f )
; INPUTS:
;	x = array of independent variable values, strictly increasing.
;	f = array of function values at points x, must be all positive.
; KEYWORDS:
;	XRANGE = two element vector specifying other Limits of integration.
;		If just one number it is assumed to be the Lower Limit.
;		If requested range exceeds tabulated range of function
;		the values are extrapolated. Default is range of x array.
;	/EXP : integrate (and interpolate) in Linear-Log space,
;		increases accuracy for exponentially varying functions.
;		In this case f array is thresholded > 1.e-37.
;	/POW : integrate (and interpolate) in Log-Log space,
;		increases accuracy for functions with power-law variation.
;		In this case both x and f arrays are thresholded > 1.e-37.
;	LOG_AXES = alternative indicator of integration method:
;			0 = Linear-Linear space (default, usualy trapezoidal).
;			1 = Linear-Log	(equiv to /EXP)
;			2 = Log-Log	(equiv to /POW)
;	/VECTOR : let the upper limit of integral vary over the x-values
;		and return a vector of cumulative definite integrals.
;	/DELTA : return a vector of definite integrals over
;		the intervals between the x-values (i.e. do not accumulate).
; OUTPUTS:
;	Returns approximate integral(s), scalar or vector depending on keywords.
; EXAMPLES:
;	fint = Trap_Int( x, f, XRAN=[.2,.9] )	;integrate over sub-range.
;	fint = Trap_Int( x, f, /EXP )		;exponentially varying function.
;	fint = Trap_Int( x, f, /POW )		;power-law varying function.
;	fintx = Trap_Int( x, f, /VEC )		;vector of integrals a func(x).
;	fintd = Trap_Int( x, f, /DELTA )	;vector of interval integrals.
; EXTERNAL CALLS:
;	function Trapez
;	function Trapex
;	function Trapow
; PROCEDURE:
;	Interpolate f-values at XRANGE points if necessary,
;	using Linear-Linear, Linear-Log, or Log-Log interpolation depending
;	on users specification of function variation.
;	If interpolation is performed, then create new xi & fi arrays
;	and call function Trap_Int recursively.
;	If /EXP or /POW recast trapeziodal rule into Lin-Log or Log-Log space,
;	by calling Trapex or Trapow respectively.
; HISTORY:
;	Written: Frank Varosi, HSTX @ NASA/GSFC, 1988-1995.
;-

function Trap_Int, x, f, XRANGE=xran, LOG_AXES=Log_axes, EXPONENTIAL=varexp, $
			POWER_LAW=varpow, VECTOR=vector, DELTA_VECTOR=delta

	if keyword_set( varexp ) then Log_axes = 1
	if keyword_set( varpow ) then Log_axes = 2
	if N_elements( Log_axes ) NE 1 then Log_axes = 0
	Log_axes = Log_axes < 2

;Check for specified integration Limits:

	if N_elements( xran ) GT 0 then begin	;setup interpolation of f(xran):

		maxx = max( x, MIN=minx )
		xmin = xran(0)
		if N_elements( xran ) GT 1 then xmax = xran(1) else xmax = maxx
		xmin = xmin < xmax
		xmax = xmax > xmin

		if (xmin GT maxx) OR (xmax LT minx) then begin
			message,/INFO, $
				"cannot integrate outside range of definition"
			return,0
		 endif else if (xmin EQ xmax) then begin
			message,/INFO,"specified X-range has zero length"
			return,0
		  endif

		w = where( (x GE xmin) AND (x LE xmax), nw )

		if (nw LE 0) then begin		;xran is between x-grid points.

			wmin = where( x LE xmin, nmin ) > 0
			wmax = where( x GE xmax ) > 0
			wmin = [ wmin(nmin-1), wmax(0) ]
			wmax = wmin
			
		 endif else begin

			if (xmin LE minx) then begin
				if (xmin LT minx) then $
					message,/INFO,"extrapolating Xmin"
				wmin = [0,1]
			  endif else wmin = w(0) - [1,0]

			if (xmax GE maxx) then begin
				if (xmax GT maxx) then $
					message,/INFO,"extrapolating Xmax"
				wmax = N_elements( x ) - [2,1]
			  endif else wmax = w(nw-1) + [0,1]
		  endelse

; interpolate f(xmin) & f(xmax) (result -> fi) depending on Log_axes type:

		wi = transpose( [ [wmin], [wmax] ] )
		xi = [ xmin, xmax ]
		if (Log_axes GT 0) then fw = aLog( f(wi)>1e-37 ) else fw = f(wi)

		if (Log_axes EQ 2) then begin
			xw = aLog( x(wi) > 1e-37 )
			xi = aLog( xi > 1e-37 )
		  endif else xw = x(wi)

		fi = fw(*,0) + (xi-xw(*,0))*(fw(*,1)-fw(*,0))/(xw(*,1)-xw(*,0))
		if (Log_axes GT 0) then fi = exp( fi )

		if (nw GT 0) then begin		;check for duplicates:

			if (xmin LT x(w(0))) then begin
				if (xmax GT x(w(nw-1))) then begin
					fi = [ fi(0), f(w), fi(1) ]
					xi = [ xmin, x(w), xmax ]
				 endif else begin
					fi = [ fi(0), f(w) ]
					xi = [ xmin, x(w) ]
				  endelse
			 endif else if (xmax GT x(w(nw-1))) then begin
				fi = [ f(w), fi(1) ]
				xi = [ x(w), xmax ]
			  endif else begin
				fi = f(w)
				xi = x(w)
			   endelse

		  endif else xi = [ xmin, xmax ]

		return, Trap_Int( xi,fi,LOG=Log_axes,VEC=vector,DEL=delta )
	   endif

;Perform the integration:

	CASE Log_axes OF

	     2: BEGIN
			if min( x ) LE 0 then begin
				message,/INFO,$
				  "restricting integral to positive values of X"
				w = where( x GE 1e-37, nw )
				if (nw LT 2) then return,0 $
				 else return,Trapow( f(w)>1e-37, x(w), $
							VEC=vector, DEL=delta )
			   endif

			return, Trapow( f > 1e-37, x, VEC=vector, DEL=delta )
		END

	     1: BEGIN
			return, Trapex( f > 1e-37, x, VEC=vector, DEL=delta )
		END

	  else: BEGIN
			return, Trapez( f, x, VEC=vector, DEL=delta )
		END

	 ENDCASE
end
