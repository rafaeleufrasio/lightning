function gw_proposal_z, Nparallel, a, error_check=error_check
;+
; Name
; ----
;   GW_PROPOSAL_Z
;
; Purpose
; -------
;   Samples a specified number of values from the stretch move proposal distribution
;   from the Goodman & Weare (2010) affine invariant MCMC sampling algorithm.
;
; Calling Sequence
; ----------------
;   ::
;
;       z = gw_proposal_z(Nparallel, a [, /error_check])
;
; Inputs
; ------
;   ``Nparallel`` : int, float, or double scalar
;       Number of samples to draw (i.e., one sample for each ensemble walker).
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
;   ``z`` : double array(Nparallel)
;       The proposal scaling constant drawn from the proposal distribution.
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
;-
 On_error, 2
 Compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(Nparallel) eq 0 then message, 'Variable is undefined: NPARALLEL.'
   if size(Nparallel, /type) lt 2 or size(Nparallel, /type) gt 5 then $
     message, 'NPARALLEL is not of type int, float, or double.'
   if size(Nparallel, /n_dim) ne 0 then message, 'NPARALLEL must be a scalar.'
   if Nparallel le 0 then message, 'NPARALLEL must be a positive value.'

   if n_elements(a) eq 0 then message, 'Variable is undefined: A.'
   if size(a, /type) lt 2 or size(a, /type) gt 5 then $
     message, 'A is not of type int, float, or double.'
   if size(a, /n_dim) ne 0 then message, 'A must be a scalar.'
   if a lt 1 then message, 'A must be a value greater than or equal to 1.'
 endif

 u = randomu(seed, Nparallel, /double)

 ; This is a nice compact version of the
 ; inverse CDF from [1/a, a] where the CDF is
 ; monotonic
 z = (1.0d / a) * (u * (a - 1.0d) + 1.0d)^2.d

 return, z

end