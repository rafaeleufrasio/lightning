;+
; This file contains some deprecated functionality of the X-ray model,
; which is no longer used in the course of Lightning fitting.
;
; Functions in this file are subject to deletion at any time.
;-

function flux_in_box, wave_obs, spectrum, box_bounds
    ; Compute mean flux of a spectrum inside a box or boxes,
    ; given the spectrum/a and the bounds of the box(es).
    ; The bounds of the boxes are assumed to be given in the same units
    ; as wave_obs; if this is not the case the result (if any) is not well
    ; behaved.
    ; The shape of box_bounds is assumed to be (2, Nboxes)

    compile_opt IDL2

    ; In all the other functions we have spectra/sed arrays with shape
    ; (N_points, N_parallel). We will preserve that here and expect the
    ; input to have this shape.
    dims = size(spectrum, /DIMENSIONS)
    N_wave = dims[0]
    N_spectra = dims[1]
    N_boxes = n_elements(box_bounds) / 2

    if (N_wave ne n_elements(wave_obs)) then message, 'ERROR -- Number of elements in provided wavelength array does not match number of elements in provided spectra'

    out_flux = dblarr(N_boxes, N_spectra)

    ; For each spectrum...
    for i=0, N_spectra - 1 do begin
        ; convolve it with each filter
        for j=0, N_boxes - 1 do begin
            if ((min(box_bounds[*,j]) lt min(wave_obs)) or (max(box_bounds[*,j]) gt max(wave_obs))) then $
                message, "ERROR -- Provided box filter bounds extend beyond bounds of spectral model"
            filter = dblarr(n_elements(wave_obs))
            filter[where((wave_obs ge min(box_bounds[*,j])) and (wave_obs le max(box_bounds[*,j])))] = 1.d0
            filter = filter / trap_int(wave_obs, filter)
            out_flux[j,i] = trap_int(wave_obs, filter * spectrum[*,i])
        endfor
    endfor

    return, out_flux
end
