function get_SKIRTOR_grid_cell, tau, delta, inc, idcs=idcs
;+
; Name
; ----
;   GET_SKIRTOR_GRID_CELL
;
; Purpose
; -------
;   This function finds the corners of the cube in the model space of the SKIRTOR grid
;   that encloses a given point (tau, delta, inc).
;
; Calling Sequence
; ----------------
;   ``Result = GET_SKIRTOR_GRID_CELL(tau, delta, inc)``
;
; Inputs
; ------
;   tau : scalar
;       Optical depth at 9.7 um.
;   delta : scalar
;       Opening angle of the dusty cone of the torus. For practical purposes we fix this at 40 degrees.
;   inc : scalar
;       Inclination from the polar axis to the line of sight.
;
; Outputs
; -------
;   This function returns a shape (2,3) array with the lower and upper bounds of each dimension,
;   one dimension per row.
;
; Optional Outputs
; ----------------
;   idcs : int array(2,3)
;       This is an array the same shape as the output, where the values correspond to the indices of the bounds
;       on the grids for each dimension.
;
; Modification History
; --------------------
;   - 2021-03-03: Created (Erik B. Monson)
;   - 2022-03-15: Documentation update (Erik B. Monson)
;-

    compile_opt IDL2

    tau_grid = [3, 5, 7, 9, 11]
    delta_grid = [10, 20, 30, 40, 50, 60, 70, 80]
    inc_grid = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]

    if(tau eq 3) then begin
        tau_lower_idx = 0
        tau_upper_idx = 1
    endif else begin
        if(tau eq 11) then begin
            tau_lower_idx = 3
            tau_upper_idx = 4
        endif else begin
            tau_lower_idx = where(tau_grid lt tau, /NULL)
            if(tau_lower_idx eq !null) then message, "Tau value is off-grid!"
            tau_lower_idx = tau_lower_idx[-1]
            tau_upper_idx = tau_lower_idx + 1
        endelse
    endelse

    if(delta eq 10) then begin
        delta_lower_idx = 0
        delta_upper_idx = 1
    endif else begin
        if(delta eq 80) then begin
            delta_lower_idx = 6
            delta_upper_idx = 7
        endif else begin
            delta_lower_idx = where(delta_grid lt delta, /NULL)
            if(delta_lower_idx eq !null) then message, "Opening angle value is off-grid!"
            delta_lower_idx = delta_lower_idx[-1]
            delta_upper_idx = delta_lower_idx + 1
        endelse
    endelse

    if(inc eq 0) then begin
        inc_lower_idx = 0
        inc_upper_idx = 1
    endif else begin
        if(inc eq 90) then begin
            inc_lower_idx = 8
            inc_upper_idx = 9
        endif else begin
            inc_lower_idx = where(inc_grid lt inc, /NULL)
            if(inc_lower_idx eq !null) then message, "Inclination value is off-grid!"
            inc_lower_idx = inc_lower_idx[-1]
            inc_upper_idx = inc_lower_idx + 1
        endelse
    endelse

    tau_lower = tau_grid[tau_lower_idx] & tau_upper = tau_grid[tau_upper_idx]
    delta_lower = delta_grid[delta_lower_idx] & delta_upper = delta_grid[delta_upper_idx]
    inc_lower = inc_grid[inc_lower_idx] & inc_upper = inc_grid[inc_upper_idx]

    bounds = [[tau_lower, tau_upper], [delta_lower, delta_upper], [inc_lower, inc_upper]]
    idcs = [[tau_lower_idx, tau_upper_idx], [delta_lower_idx, delta_upper_idx], [inc_lower_idx, inc_upper_idx]]

    return, bounds
end
