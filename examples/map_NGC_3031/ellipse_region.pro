function ellipse_region, ra, dec, a, b, pa, hd
;+
; Name
; ----
;   ELLIPSE_REGION
;
; Purpose
; -------
;   Generates an ellipse region mask defined by the input RA, DEC, and
;   ellipse shape arguments.
;
; Calling Sequence
; ----------------
;   ::
;
;       mask = ellipse_region(ra, dec, a, b, pa, hd)
;
; Inputs
; ------
;   ``ra`` : int, float, or double scalar
;       The Right Ascension of the center of the ellipse [degrees].
;   ``dec`` : int, float, or double scalar
;       The Declination of the center of the ellipse [degrees].
;   ``a`` : int, float, or double scalar
;       The semi-major axis size of the ellipse [degrees].
;   ``b`` : int, float, or double scalar
;       The semi-minor axis size of the ellipse [degrees].
;   ``pa`` : int, float, or double scalar
;       The North Eastwards position angle of the ellipse [degrees].
;   ``hd`` : string array(Nheader)
;       The FITS header of the image for which the mask is to be created. The
;       header must contain the non-standard keyword of ``PLTSCALE``, which
;       defines the pixel scale in arcsec.
;
; Output
; ------
;   ``mask`` : int array(Naxis1, Naxis2)
;       The ellipse region mask. Values of 1 indicate pixels within the ellipse, 
;       and values of 0 indicated pixels outside of the ellipse.
;
; Modification History
; --------------------
;   - 2023/02/07: Created (Keith Doore)
;-
 On_error, 2
 compile_opt IDL2

; Error Handling
 if n_elements(ra) eq 0 then message, 'Variable is undefined: RA.'
 if size(ra, /type) lt 2 or size(ra, /type) gt 5 then message, 'RA must be of type int, float, or double.'
 if size(ra, /n_dim) ne 0 then message, 'RA must be a scalar.'
 if min(ra) lt 0 then message, 'RA must only contain positive values.'
 
 if n_elements(dec) eq 0 then message, 'Variable is undefined: DEC.'
 if size(dec, /type) lt 2 or size(dec, /type) gt 5 then message, 'DEC must be of type int, float, or double.'
 if size(dec, /n_dim) ne 0 then message, 'DEC must be a scalar.'
 if min(dec) lt -90 or max(dec) gt 90 then message, 'DEC must only contain values between -90 and 90.'

 if n_elements(a) eq 0 then message, 'Variable is undefined: A.'
 if size(a, /type) lt 2 or size(a, /type) gt 5 then message, 'A must be of type int, float, or double.'
 if size(a, /n_dim) ne 0 then message, 'A must be a scalar.'
 if min(a) lt 0 then message, 'A must only contain positive values.'

 if n_elements(b) eq 0 then message, 'Variable is undefined: B.'
 if size(b, /type) lt 2 or size(b, /type) gt 5 then message, 'B must be of type int, float, or double.'
 if size(b, /n_dim) ne 0 then message, 'B must be a scalar.'
 if min(b) lt 0 then message, 'B must only contain positive values.'

 if n_elements(pa) eq 0 then message, 'Variable is undefined: PA.'
 if size(pa, /type) lt 2 or size(pa, /type) gt 5 then message, 'PA must be of type int, float, or double.'
 if size(pa, /n_dim) ne 0 then message, 'PA must be a scalar.'

 if n_elements(hd) eq 0 then message, 'Variable is undefined: HD.'
 if size(hd, /type) ne 7 then message, 'HD must be of type string.'
 if size(reform(hd), /n_dim) ne 1 then message, 'HD must be a 1-D array.'
 if sxpar(hd, 'PLTSCALE') le 0 then message, "HD must contain the 'PLTSCALE' keyword."


; Create mask
; Set the dimensions of the image
 dim = [sxpar(hd, 'NAXIS1'), sxpar(hd, 'NAXIS2')]
; Create and X and Y grid of locations
 dim_grid_x = indgen(dim[0]) # (intarr(1, dim[1]) + 1)
 dim_grid_y = (intarr(dim[0]) + 1) # indgen(dim[1])
 mask = intarr(dim)

; Convert ra and dec of input to image pixel locations
 extast, hd, astr
 ad2xy, ra, dec, astr, x_center, y_center

; Most headers have the plate scale in arcsec, so we assume that here
 pltscale = sxpar(hd, 'PLTSCALE') / 3.6d3
; Convert position angle to radians
 pa_rad = pa * !pi/1.8d2

 ellipse_values=(((dim_grid_x - x_center) * cos(pa_rad) + (dim_grid_y-y_center) * sin(pa_rad))^2.d0 / $
                 ((b/pltscale)^2.d0) + $
                 ((dim_grid_x - x_center) * sin(pa_rad) - (dim_grid_y-y_center) * cos(pa_rad))^2.d0 / $
                 ((a/pltscale)^2.d0))

 mask[where(ellipse_values le 1, /null)] = 1

 return, mask

end
