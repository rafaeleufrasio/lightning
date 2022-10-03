function doore21_interp_lbol_abs_table, tauB_f, F_clump, b_to_d, rold0_ages, doore21_Lbol_abs_table, $
                               error_check=error_check
;+
; Name
; ----
;   DOORE21_INTERP_LBOL_ABS_TABLE
;
; Purpose
; -------
;   Linearly Interpolates the bolometric luminosity of the attenuated stellar
;   emission from precomputed tables for the Doore et al. (2021) attenuations curves.
;
; Calling Sequence
; ----------------
;   ::
;
;       steps_Lbol_abs = doore21_interp_lbol_abs_table(tauB_f, F_clump, b_to_d, $
;                             rold0_ages, doore21_Lbol_abs_table [, /error_check])
;
; Inputs
; ------
;   ``tauB_f`` : int, float, or double array(Nmodels)
;       The face-on optical depth in the B-band.
;   ``F_clump`` : int, float, or double array(Nmodels)
;       The clumpiness factor F.
;   ``b_to_d`` : int, float, or double array(Nmodels)
;       The bulge-to-disk ratio.
;   ``rold0_ages`` : int, float, or double array(Nsteps)
;       The binary parameter ``rold0``, designating each SFH age bin as part of
;       the young or old population when using the Doore et al. (2021) attenuation
;       curves. A value of ``0`` for the corresponding age bin considers it to be part of
;       the young population, and a value of ``1`` considers it to be part of the old
;       populations (see Doore et al. 2021 for more details).
;   ``doore21_Lbol_abs_table`` : structure
;       A structure giving the pre-computed table of ``steps_Lbol_abs`` for a given redshift
;       to have energy conservation with the Doore et al. (2021) attenuation curves.
;       (See ``generate_doore21_lbol_abs_table.pro`` for details and contents.)
;
; Optional Input
; --------------
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``steps_Lbol_abs`` : double array(Nsteps, Nmodels)
;       The interpolated bolometric luminosity of the attenuated stellar emission for
;       each SFH step :math:`[L_\odot]`.
;
; Reference
; ---------
;   `Doore, K., Eufrasio, R. T., Lehmer, B. D., et al. 2021, ApJ, 923, 26 <https://ui.adsabs.harvard.edu/abs/2021ApJ...923...26D/abstract>`_
;
; Modification History
; --------------------
;   - 2022/04/15: Created (Keith Doore)
;   - 2022/04/18: Added error handling (Keith Doore)
;   - 2022/06/09: Added ``error_check`` keyword to do error handling (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(tauB_f) eq 0 then message, 'Variable is undefined: TAUB_F.'
   if size(tauB_f, /type) lt 2 or size(tauB_f,/type) gt 5 then message, 'TAUB_F must be of type int, float, or double.'
   if size(reform(tauB_f), /n_dim) ne 1 then message, 'TAUB_F must be a scalar or 1-D array.'
   if min(tauB_f) lt 0 or max(tauB_f) gt 8 then message, 'TAUB_F must contain values between 0 and 8.'
  
   if n_elements(F_clump) eq 0 then message, 'Variable is undefined: F_CLUMP.'
   if size(F_clump, /type) lt 2 or size(F_clump,/type) gt 5 then message, 'F_CLUMP must be of type int, float, or double.'
   if size(reform(F_clump), /n_dim) ne 1 then message, 'F_CLUMP must be a scalar or 1-D array.'
   if min(F_clump) lt 0.d0 or max(F_clump) gt 0.61d0 then message, 'F_CLUMP must contain values between 0 and 0.61.'
  
   if n_elements(b_to_d) eq 0 then message, 'Variable is undefined: B_TO_D.'
   if size(b_to_d, /type) lt 2 or size(b_to_d,/type) gt 5 then message, 'B_TO_D must be of type int, float, or double.'
   if size(reform(b_to_d), /n_dim) ne 1 then message, 'B_TO_D must be a scalar or 1-D array.'
   if min(b_to_d) lt 0 then message, 'B_TO_D must only contain non-negative values.'
  
   if n_elements(tauB_f) ne n_elements(b_to_d) or n_elements(tauB_f) ne n_elements(F_clump) then $
      message, 'TAUB_F, F_CLUMP, and B_TO_D must have the same size.'
  
   if n_elements(rold0_ages) eq 0 then message, 'Variable is undefined: ROLD0_AGES.'
   if size(rold0_ages, /type) lt 2 or size(rold0_ages,/type) gt 5 then message, 'ROLD0_AGES must be of type int, float, or double.'
   if size(rold0_ages, /n_dim) ne 1 then message, 'ROLD0_AGES must be a 1-D array.'
   if total(rold0_ages eq 0 or rold0_ages eq 1) ne n_elements(rold0_ages) then message, 'ROLD0_AGES must only contain values of 0 or 1.'
  
   if n_elements(doore21_Lbol_abs_table) eq 0 then message, 'Variable is undefined: DOORE21_LBOL_ABS_TABLE.'
   if size(doore21_Lbol_abs_table, /type) ne 8 then message, 'DOORE21_LBOL_ABS_TABLE must be of type structure.'
   if doore21_Lbol_abs_table.nsteps ne n_elementS(rold0_ages) then $
      message, 'ROLD0_AGES must have the same number of elements as the number of steps in DOORE21_LBOL_ABS_TABLE'
 endif

; Interpolate table
 Nmodels = n_elements(taub_f)
 rdisk0 = 1.d0/(1+reform(b_to_d))

 Lbol_abs_model = doore21_Lbol_abs_table.Lbol_abs_model
 ntauB_f = n_elements(doore21_Lbol_abs_table.tauB_f_grid)
 nrdisk0 = n_elements(doore21_Lbol_abs_table.rdisk0_grid)
 nF = n_elements(doore21_Lbol_abs_table.F_grid)

 taub_loc = interpol(lindgen(ntaub_f), doore21_Lbol_abs_table.tauB_f_grid, reform(tauB_f)) > 0
 F_loc    = interpol(lindgen(nF), doore21_Lbol_abs_table.F_grid, reform(F_clump)) > 0
 rdisk0_loc  = interpol(lindgen(nrdisk0), doore21_Lbol_abs_table.rdisk0_grid, rdisk0) > 0

 steps_Lbol_abs = dblarr(doore21_Lbol_abs_table.nsteps, Nmodels)

 young_steps = where(rold0_ages eq 0, num_young, /null)
 old_steps = where(rold0_ages eq 1, num_old, /null)
 if num_young gt 0 then $
   steps_Lbol_abs[young_steps,*] = interpolate(Lbol_abs_model[young_steps,*,*,0], taub_loc, F_loc, /double)
 if num_old gt 0 then $
   steps_Lbol_abs[old_steps,*] = interpolate(Lbol_abs_model[old_steps,*,*,1], taub_loc, rdisk0_loc, /double)

 return, steps_Lbol_abs

end