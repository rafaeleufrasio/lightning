function unique, array, Nuniq, RESOLUTION=resol, VERBOSE=verb, $
			SORT_FIRST=sortf, ASORT=asort, ORIGINAL_ORDER=orig
;+
; NAME:
;	unique
; PURPOSE:
;	Return subscripts of unique (first occurence) elements in array,
;	(thereby eliminating adjacent duplicates),
;	or optionally sort array before checking for duplicates,
;	and/or optionally return unique subscripts in original order.
;	If array is numerical: approximate duplicates (differing by epsilon)
;	can be found and eliminated by keyword RESOLUTION=epsilon.
;
; CALLING EXAMPLES:
;
;	u = unique( array, RESOL=.77, /VERB, /SORT )
;	x = array( unique( array, /SORT, /ORIGINAL ) )
;
; INPUTS:
;	array = array of numbers or strings.
;
; KEYWORDS:
;	RESOLUTION = the accuracy by which to determine equivalence,
;		if elements differ by RESOLUTION or less they are considered 
;		equal, (array must then be numerical, not strings).
;		Default is RESOL=0 (duplicates are defined by exact equality).
;
;	ASORT = optional output of the sorted array, if /SORT specified.
;	/VERBOSE to get brief message about # equal values.
;	/SORT to sort array first, and return sorted unique subscripts.
;	/ORIGINAL to return the unique subscripts in original order.
;
; OUTPUTS:
;	Nuniq = # of unique occurences found (optional output).
;
;	Function returns subscripts of unique elements in array.
;
; EXTERNAL CALLS:
;	function Fsort	(to sort and maintain original order when duplicates).
; PROCEDURE:
;	Create sorted array first (if requested) then
;	use the IDL intrinsic function "where" to find differences.
;	If original order is requested, just sort the unique subscripts at end.
; HISTORY:
;	Written: Frank Varosi NASA/GSFC 1990.
;	F.V. 1992, added /ORIGINAL_ORDER option.
;-
	sa = size( array )
	Nuniq = 1
	if sa(sa(0)+2) LE 1 then return,[0]

if sa(sa(0)+1) EQ 8 then begin
	message,"does not work on structures",/INFO
	Nuniq = N_elements( array )
	return,indgen( Nuniq )
   endif

if N_elements( resol ) EQ 1 then begin
	if (resol GT 0) then begin
		if sa(sa(0)+1) EQ 7 then begin
			message,"only exact equality for strings",/INFO
			resol=0
		   endif
	   endif
   endif else resol=0

if (resol GT 0) then begin

   if keyword_set( sortf ) then begin

	ais = Fsort( array, asort )
	wuniq = ais( where( abs( shift(asort,1)-asort ) GT resol, Nuniq ) > 0 )

     endif else wuniq = where( abs( shift(array,1)-array ) GT resol, Nuniq ) >0

   Nuniq = Nuniq > 1
   if keyword_set( verb ) then message, "Located " + strtrim( Nuniq, 2 ) + $
					" values unique by + or - " + $
					strtrim( resol, 2 ) + " out of " + $
					strtrim( sa(sa(0)+2),2 ),/INFO
 endif else begin

   if keyword_set( sortf ) then begin

		ais = Fsort( array, asort )
		wuniq = ais( where( shift(asort,1) NE asort, Nuniq ) > 0 )

     endif else wuniq = where( shift(array,1) NE array, Nuniq ) > 0

   Nuniq = Nuniq > 1
   if keyword_set( verb ) then message, "Located " + strtrim( Nuniq, 2 ) + $
					" unique values out of " + $
					strtrim( sa(sa(0)+2),2 ),/INFO
  endelse

if keyword_set( orig ) AND (Nuniq GT 1) then return, wuniq( sort( wuniq ) ) $
					else return, wuniq
end
