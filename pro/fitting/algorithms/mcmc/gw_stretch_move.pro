function gw_stretch_move, ensemble, a, error_check=error_check, z=z
;+
; Name
; ----
;   GW_STRETCH_MOVE
;
; Purpose
; -------
;   Proposes a new step for an MCMC ensemble following the "stretch move"
;   prescription in Goodman & Weare (2010).
;
; Calling Sequence
; ----------------
;   ::
;
;       ensemble_new = gw_stretch_move(ensemble, a [, /error_check, z=z])
;
; Inputs
; ------
;   ``ensemble`` : float or double array(Nparam, Nparallel)
;       The current state of all the chains in the MCMC ensemble at a given trial.
;   ``a`` : int, float, or double scalar
;       A real constant :math:`\geq 1` which controls the size of the proposal 
;       distribution (i.e., how close to 0 ``z`` is allowed to be and how 
;       large ``z`` is allowed to be). In practice, we can only sample 
;       values of ``z`` in ``[1/a, a]``.
;
; Optional Input
; --------------
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``ensemble_new`` : double array(Nparam, Nparallel)
;       The new proposed positions of the ensemble.
;
; Optional Output
; ---------------
;   ``z`` : double array(Nparallel)
;       The proposal scaling constant.
;
; Reference
; ---------
;   `Goodman, J., & Weare, J. 2010, CAMCS, 5, 65 <https://ui.adsabs.harvard.edu/abs/2010CAMCS...5...65G/abstract>`_
;
; Modification History
; --------------------
;   - 2021/12/02: Created (Erik B. Monson)
;   - 2022/04/21: Documentation (Erik B. Monson)
;   - 2022/06/17: Renamed variables to match naming scheme (Keith Doore)
;   - 2022/06/17: Updated documentation (Keith Doore)
;   - 2022/06/17: Added proper error handling (Keith Doore)
;   - 2022/06/17: Added ``error_check`` keyword to do error handling (Keith Doore)
;   - 2022/06/17: Removed ``bounds_arr`` input and subsequent code as this is now done in prior (Keith Doore)
;   - 2022/06/17: Removed ``seed`` input as we do not ever specify the random seed (Keith Doore)
;   - 2022/06/17: Made ``z`` a separate output from ``ensemble_new`` to prevent any issues or confusion (Keith Doore)
;   - 2022/07/05: Changed from ``round`` to ``floor`` in ``kprime`` as index could exceed array size (Keith Doore)
;-
 On_error, 2
 Compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(ensemble) eq 0 then message, 'Variable is undefined: ENSEMBLE.'
   if size(ensemble, /type) lt 4 or size(ensemble, /type) gt 5 then $
     message, 'ENSEMBLE is not of type float or double.'
   if size(ensemble, /n_dim) ne 2 then message, 'ENSEMBLE must be a 2-D array.'

   if n_elements(a) eq 0 then message, 'Variable is undefined: A.'
   if size(a, /type) lt 2 or size(a, /type) gt 5 then $
     message, 'A is not of type int, float, or double.'
   if size(a, /n_dim) ne 0 then message, 'A must be a scalar.'
   if a lt 1 then message, 'A must be a value greater than or equal to 1.'
 endif


; Let Nparallel be the number of walkers in the ensemble
; and let Nparam be the number of variable parameters.
 dims = size(ensemble, /dim)
 Nparam = dims[0]
 Nparallel = dims[1]

; Draw the scaling constant from the proposal distribution
 z = gw_proposal_z(Nparallel, a)

 ensemble_new = dblarr(Nparam, Nparallel)

 for k = 0, Nparallel - 1 do begin
   ; Pick a complementary ensemble
   ; k' from {0 .. Nparallel - 1} \ {k}
   kcompl = where(lindgen(Nparallel) ne k)
   kprime = kcompl[floor(randomu(seed) * (Nparallel - 1), /l64)]

  ; Propose a new point on the line through the current state of chains k and kprime
   proposal = ensemble[*, kprime] + replicate(z[k], Nparam) * (ensemble[*, k] - ensemble[*, kprime])
   ensemble_new[*, k] = proposal
 endfor

 return, ensemble_new

end
