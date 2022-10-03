function mcmc_convergence, chain, lnprob_chain, accepted_trials, C_step=C_step, tolerance=tolerance, $
                           burn_in=burn_in, thin_factor=thin_factor, final_chain_length=final_chain_length, $
                           method=method
;+
; Name
; ----
;   MCMC_CONVERGENCE
;
; Purpose
; -------
;   Tests the resulting MCMC chain for convergence. This includes the
;   autocorrelation time for both the affine-invariant and adaptive 
;   methods. Additionally, the Gelman-Rubin statistic as updated by
;   Brooks & Gelman (1998) is computed for the adaptive method.
;
; Calling Sequence
; ----------------
;   ::
;
;       convergence_metric = mcmc_convergence(chain, lnprob_chain, accepted_trials [, C_step = , $
;                                             tolerance = , burn_in = , thin_factor = , $
;                                             final_chain_length = , method = ])
;
; Input
; -----
;   ``chain`` : int, float, or double array(Nparam, Ntrials, Nparallel)
;       MCMC chain or ensemble of MCMC chains.
;   ``lnprob_chain`` : int, float, or double array(Ntrials, Nparallel)
;       Log-probability of each element in the MCMC chain or ensemble
;       of MCMC chains.
;   ``accepted_trials`` : int, float, or double array(Nparallel)
;       The number of accepted MCMC trials for each MCMC chain.
;
; Optional Inputs
; ---------------
;   ``C_step`` : int, float, or double scalar
;       Defines how many trials of the chain are used to calculate tau, where
;       we integrate tau to the smallest index ``M`` such that ``M > C_step * tau``.
;       (Default = ``5``)
;   ``tolerance`` : int, float, or double scalar
;       Defines how many taus the length of the chain should be for us to 
;       believe the estimate. (Default = ``50``)
;   ``burn_in`` : int, float, or double scalar
;       The number of initial MCMC trials to truncate as the burn-in phase. If set to ``0``,
;       then the number will chosen automatically from the autocorrelation time.
;       (Default = ``0``)
;   ``thin_factor`` : int, float, or double scalar
;       The factor to thin the MCMC chain after removing the burn-in trials. If set to ``0``,
;       then the number will be chosen automatically from the autocorrelation time.
;       (Default = ``0``)
;   ``final_chain_length`` : int, float, or double scalar
;       The number of MCMC trials to include for the final distributions as taken from 
;       end of the chain (only for adaptive MCMC method). (Default = ``1000``)
;   ``method`` : string scalar
;     The fitting algorithm used to fit the SED(s). Current options are: 
;     ``'MCMC-ADAPTIVE'`` and ``'MCMC-AFFINE'``. (Default = ``'MCMC-AFFINE'``)
;
; Output
; ------
;   ``convergence_metric`` : structure
;       A structure containing the convergence metrics and flags used to 
;       indicate if convergence was reached. (A flag of ``1`` indicates failure
;       to converge.)
;       The full description of the structure is as follows:
;
;       ===================     =================     ==================================================================================
;       TAG                     TYPE                  DESCRIPTION
;       ===================     =================     ==================================================================================
;       ACCEPTANCE_FRAC         double(Nparallel)     Fraction of accepted trials for each MCMC chain
;       AUTOCORR_TIME           double(Nparam)        Autocorrelation time of each parameter for the ensemble of MCMC chains
;       BURN_IN_AUTOCORR        int                   The number of burn-in trials as determined from the autocorrelation time
;       GELMAN_RUBIN_R_HAT      double(Nparam)        The Gelman-Rubin statistic (only for adaptive MCMC, ``NaN`` otherwise)
;       BROOKS_GELMAN_R_HAT     double                The Brooks-Gelman statistic (only for adaptive MCMC, ``NaN`` otherwise)
;       NCHAIN_R_HAT            int                   Number of parallel chains with ``lnprob`` < 2 above best fit chain element
;       ACCEPTANCE_FLAG         int(Nparallel)        Flag indicating if acceptance fraction is outside optimal range (``0.2-0.5``)
;       AUTOCORR_FLAG           int(Nparam)           Flag indicating if autocorrelation time is above specified tolerance
;       GELMAN_RUBIN_FLAG       int(Nparam)           Flag indicating if Gelman-Rubin stat is too large (``r_hat >= 1.2``)
;       BROOKS_GELMAN_FLAG      int                   Flag indicating if Brooks-Gelman stat is too large (``r_hat >= 1.2``)
;       CONVERGENCE_FLAG        int                   Flag indicating if any other flag was issued for a convergence metric
;       ===================     =================     ==================================================================================
;
; References
; ----------
;   - `Gelman, A., & Rubin, D. B. 1992, StaSc, 7, 457 <https://ui.adsabs.harvard.edu/abs/1992StaSc...7..457G/abstract>`_
;   - `Brooks, S. P., & Gelman, A. 1998, Journal of Computational and Graphical Statistics, 7, 434 <https://www.tandfonline.com/doi/abs/10.1080/10618600.1998.10474787>`_
;   - `Goodman, J., & Weare, J. 2010, CAMCS, 5, 65 <https://ui.adsabs.harvard.edu/abs/2010CAMCS...5...65G/abstract>`_
;
; Modification History
; --------------------
;   - 2022/07/14: Created (Keith Doore)
;   - 2022/08/03: Updated adaptive ``r_hat`` metric calculations (Keith Doore)
;   - 2022/08/03: Added optional ``burn_in`` input, and renamed automatic ``burn-in`` from autocorr to ``burn_in_autocorr`` (Keith Doore)
;   - 2022/08/03: Added optional ``thin_factor`` input to include full pre-thinned chain for MCMC-ADAPTIVE ``r_hat`` stat (Keith Doore)
;   - 2022/08/17: Fixed issue with Gelman-Rubin ``r_hat`` ``NaN`` array, when ``r_hat`` not calculated (Keith Doore)
;   - 2022/08/18: Fixed issue where pre-thinned adaptive chain length may be greater than ``Ntrials`` (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(chain) eq 0 then message, 'Variable is undefined: CHAIN.'
 if size(chain, /type) lt 2 or size(chain, /type) gt 5 then message, 'CHAIN must be of type int, float, or double.'
 if size(chain, /n_dim) lt 2 or size(chain, /n_dim) gt 3 then message, 'CHAIN must be a 2-D or 3-D array.'

 if n_elements(accepted_trials) eq 0 then message, 'Variable is undefined: ACCEPTED_TRIALS.'
 if size(accepted_trials, /type) lt 2 or size(accepted_trials, /type) gt 5 then $
   message, 'ACCEPTED_TRIALS must be of type int, float, or double.'
 if size(accepted_trials, /n_dim) gt 1 then message, 'ACCEPTED_TRIALS must be a scalar or 1-D array.'
 if max(accepted_trials) lt 0 then message, 'ACCEPTED_TRIALS must only contain non-negative values.'

 if n_elements(C_step) ne 0 then begin
   if size(C_step, /type) lt 2 or size(C_step, /type) gt 5 then $
     message, 'C_STEP must be of type int, float, or double.'
   if size(C_step, /n_dim) ne 0 then message, 'C_STEP must be a scalar.'
   if C_step le 0 then message, 'C_STEP must be a positive value.'
 endif else C_step = 5.d

 if n_elements(tolerance) ne 0 then begin
   if size(tolerance, /type) lt 2 or size(tolerance, /type) gt 5 then $
     message, 'TOLERANCE must be of type int, float, or double.'
   if size(tolerance, /n_dim) ne 0 then message, 'TOLERANCE must be a scalar.'
   if tolerance le 0 then message, 'TOLERANCE must be a positive value.'
 endif else tolerance = 50.d

 if n_elements(burn_in) ne 0 then begin
   if size(burn_in, /type) lt 2 or size(burn_in, /type) gt 5 then $
     message, 'BURN_IN must be of type int, float, or double.'
   if size(burn_in, /n_dim) ne 0 then message, 'BURN_IN must be a scalar.'
   if burn_in lt 0 then message, 'BURN_IN must be a non-negative value.'
 endif else burn_in = 0.d

 if n_elements(thin_factor) ne 0 then begin
   if size(thin_factor, /type) lt 2 or size(thin_factor, /type) gt 5 then $
     message, 'THIN_FACTOR must be of type int, float, or double.'
   if size(thin_factor, /n_dim) ne 0 then message, 'THIN_FACTOR must be a scalar.'
   if thin_factor lt 0 then message, 'THIN_FACTOR must be a non-negative value.'
 endif else thin_factor = 0.d

 if n_elements(final_chain_length) ne 0 then begin
   if size(final_chain_length, /type) lt 2 or size(final_chain_length, /type) gt 5 then $
     message, 'FINAL_CHAIN_LENGTH must be of type int, float, or double.'
   if size(final_chain_length, /n_dim) ne 0 then message, 'FINAL_CHAIN_LENGTH must be a scalar.'
   if final_chain_length le 0 then message, 'FINAL_CHAIN_LENGTH must be a positive value.'
 endif else final_chain_length = 1000.d

 if n_elements(method) ne 0 then begin
   if size(method, /type) ne 7 then message, 'METHOD must be of type string.'
   if size(method, /n_dim) ne 0 then message, 'METHOD must be a scalar.'
   if total(strupcase(method) eq ['MCMC-ADAPTIVE', 'MCMC-AFFINE']) ne 1 then $
     message, "METHOD must be set to either 'MCMC-ADAPTIVE' or 'MCMC-AFFINE'."
 endif else method = 'MCMC-AFFINE'


 dim = size(chain, /dimensions)
 Nparam = dim[0]
 Ntrials = dim[1]
 if n_elements(dim) eq 2 then Nparallel = 1 else Nparallel = dim[2]


; Check both algorithm's acceptance fraction
 acceptance_frac = accepted_trials/double(Ntrials)
 acceptance_flag = intarr(Nparallel)
 acceptance_flag[where(acceptance_frac lt 0.2 or acceptance_frac gt 0.5, /null)] = 1


; Compute autocorrelation time for affine invariant MCMC
; Compute Gelman-Rubin stat and autocorrelation time for the adaptive MCMC
 case strupcase(method) of
   'MCMC-AFFINE': begin
       tau = autocorr_time(chain, C_step=C_step)

       ; Determine burn-in from autocorrelation time (2x tau). Use this if the burn-in
       ;   is to be automatically determined.
       burn_in_autocorr = ceil(2*max(tau))
       if burn_in eq 0 then burn_in = burn_in_autocorr
       ; Remove burn-in and recompute autocorrelation time to check for convergence
       ;   If burn-in is longer than chain, don't recompute autocorrelation
       if burn_in lt Ntrials then tau = autocorr_time(chain[*, burn_in:*, *], C_step=C_step)

       ; Check if new autocorrelation time is below tolerance for each parameter
       ;   Fixed parameters have "converged", since they do not change
       autocorr_flag = intarr(Nparam)
       tol_idcs = where(Ntrials lt tolerance * tau[where(finite(tau), /null)], /null)
       var_idcs = (indgen(Nparam))[where(finite(tau), /null)]
       autocorr_flag[var_idcs[tol_idcs]] = 1

       ; Gelman-Rubin and Brooks-Gelman stats can't be computed for affine invariant method
       ;   So, set to NaN and flag to 0.
       r_hat_gr = replicate(!values.D_NaN, Nparam)
       r_hat_bg = !values.D_NaN
       nchain_r_hat = !values.D_NaN
       gelman_rubin_flag = intarr(Nparam)
       brooks_gelman_flag = 0
     end

   'MCMC-ADAPTIVE': begin
       tau = autocorr_time(chain, C_step=C_step)

       ; Determine burn-in from autocorrelation time (2x tau). This is likely
       ;   not enough for most adaptive chains, but calculate as such regardless
       ;   in case user specified automatic burn-in, which they should not for the 
       ;   adaptive MCMC.
       burn_in_autocorr = ceil(2*max(tau))
       if burn_in eq 0 then burn_in = burn_in_autocorr
       ; Remove burn-in and recompute autocorrelation time to check for convergence
       ;   If burn-in is longer than chain, don't recompute autocorrelation
       if burn_in lt Ntrials then tau = autocorr_time(chain[*, burn_in:*, *], C_step=C_step)

       ; Check if new autocorrelation time is below tolerance for each parameter
       ;   Fixed parameters have "converged", since they do not change
       autocorr_flag = intarr(Nparam)
       ;tol_idcs = where(Ntrials lt tolerance * tau[where(finite(tau), /null)], /null)
       ;var_idcs = (indgen(Nparam))[where(finite(tau), /null)]
       ;autocorr_flag[var_idcs[tol_idcs]] = 1

       ; Compute Gelman-Rubin and Brooks-Gelman statistics for parallel chains
       if Nparallel gt 1 then begin
         ; Only use chains whose most likely element is within 2 of the overall most likely element.
         ;   Chains that did not get within 2 likely did not reach the global (or best local)
         ;   minimum in chisqr space. The value of 2 is arbitrarily chosen but equates to only 
         ;   having a 13.5% chance of being accepted by the MCMC sampler, which is relatively low.
         max_lnprob_diff = max(lnprob_chain) - max(lnprob_chain, dim=1) 
         reached_min = where(max_lnprob_diff lt 2, nchain_r_hat)

         ; Determine length of final portion of pre-thinned chain
         if thin_factor eq 0 then thin_factor = max(tau)*0.5d
         prethinned_length = ceil(thin_factor, /L64) * final_chain_length
         ; Check if this length is less than the number of trials, otherwise set to number of trials
         if prethinned_length gt Ntrials then prethinned_length = Ntrials

         ; Use only final portion of pre-thinned chain to check for convergence by the end of sampling
         if nchain_r_hat gt 1 then begin
           r_hat_gr = gelman_rubin_stat(chain[*, -1*prethinned_length:*, reached_min])
           r_hat_bg = brooks_gelman_stat(chain[*, -1*prethinned_length:*, reached_min])
           gelman_rubin_flag = intarr(Nparam)
           gelman_rubin_flag[where(sqrt(r_hat_gr) ge 1.2, /null)] = 1
           brooks_gelman_flag = fix(sqrt(r_hat_bg) ge 1.2)
         endif else begin
           r_hat_gr = replicate(!values.D_NaN, Nparam)
           r_hat_bg = !values.D_NaN
           gelman_rubin_flag = intarr(Nparam) + 1
           brooks_gelman_flag = 1
         endelse

       endif else begin
         ; Gelman-Rubin and Brooks-Gelman stats can't be computed for a single chain
         ;   So, set to NaN and flag to 0.
         r_hat_gr = replicate(!values.D_NaN, Nparam)
         r_hat_bg = !values.D_NaN
         nchain_r_hat = !values.D_NaN
         gelman_rubin_flag = intarr(Nparam)
         brooks_gelman_flag = 0
       endelse
     end
 endcase


; Compute overall convergence flag that checks if any other flag was issued for a convergence metric
 convergence_flag = fix(total([acceptance_flag, autocorr_flag, gelman_rubin_flag, brooks_gelman_flag]) gt 0)


; Compile convergence tests and flags into structure
 convergence_metric = {acceptance_frac: acceptance_frac,       $
                       autocorr_time: tau,                     $
                       burn_in_autocorr: burn_in_autocorr,     $
                       gelman_rubin_r_hat: r_hat_gr,           $
                       brooks_gelman_r_hat: r_hat_bg,          $
                       nchain_r_hat: nchain_r_hat,             $
                       acceptance_flag: acceptance_flag,       $
                       autocorr_flag: autocorr_flag,           $
                       gelman_rubin_flag: gelman_rubin_flag,   $
                       brooks_gelman_flag: brooks_gelman_flag, $
                       convergence_flag: convergence_flag      $
                       }


 return, convergence_metric

end