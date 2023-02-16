function mpfit_convergence, parameters, lnprob, status, error_msg, Niterations, $
                            Nfunc_evals, DoF, maxiter=maxiter
;+
; Name
; ----
;   MPFIT_CONVERGENCE
;
; Purpose
; -------
;   Tests the resulting MPFIT solvers for convergence. This includes
;   checking the MPFIT status, the fraction of max iterations used by
;   each solver to reach its solution, the p-value of each solver 
;   from a simple :math:`\chi^2` test, how many solvers got stuck in
;   low probability regions, and how similar are the solver solutions
;   to each other.
;
; Calling Sequence
; ----------------
;   ::
;
;       convergence_metric = mpfit_convergence(parameters, lnprob, status, error_msg, $
;                                              Niterations, Nfunc_evals, $
;                                              DoF [, maxiter = ])
;
; Inputs
; ------
;   ``parameters`` : int, float, or double array(Nparam, Nsolvers)
;       The best fit parameters of the model as determined by MPFIT for each
;       solver. The actual parameters contained in this array depend on the
;       chosen model during configuration.
;   ``lnprob`` : int, float, or double array(Nsolvers)
;       Log-probability for each set of parameters.
;   ``status`` : int array(Nsolvers)
;       The status code given by MPFIT.
;   ``error_msg`` : string array(Nsolvers)
;       An error or warning message given by MPFIT.
;   ``Niterations`` : int array(Nsolvers)
;       The number of MPFIT iterations completed.
;   ``Nfunc_evals`` : int array(Nsolvers)
;       The number of ``lightning_mpfit_function.pro`` evaluations performed.
;   ``DoF`` : int
;       The number of degrees of freedom in the fit.
;
; Optional Input
; --------------
;   ``maxiter`` : int, float, or double scalar
;       The maximum number of MPFIT iterations that could have been perform.
;       (Default = ``200``)
;
; Output
; ------
;   ``convergence_metric`` : structure
;       A structure containing the convergence metrics and flags used to 
;       indicate if convergence was reached. (A flag of ``1`` indicates failure
;       to converge.)
;       The full description of the structure is as follows:
;
;       ================     ================     ==========================================================================================
;       TAG                  TYPE                 DESCRIPTION
;       ================     ================     ==========================================================================================
;       STATUS               int(Nsolvers)        The status code given by MPFIT
;       STATUS_FLAG          int(Nsolvers)        Flag indicating if the MPFIT status indicated failure of the algorithm (``<= 0``)
;       ERROR_MSG            string(Nsolvers)     Error or warning message given by MPFIT. Will be blank if no message
;       ITER_FRAC            double(Nsolvers)     Fraction of max iterations used by MPFIT to reach solution
;       ITER_FLAG            int(Nsolvers)        Flag indicating if the maximum number of iteration were used by MPFIT
;       PVALUE               double(Nsolvers)     P-value of fit as determined from a :math:`\chi^2` test using ``lnprob`` and ``DoF``
;       DOF                  int                  Degrees of freedom (same for each solver)
;       STUCK_FRAC           double               Fraction of solvers with ``lnprob`` > 2 above best fit solver
;       STUCK_FLAG           int                  Flag indicating if the majority of solvers had ``lnprob`` > 2 above best fit solver
;       SIMILAR_FLAG         int(Nparam)          Flag indicating if any non-stuck solvers had different solutions (>1% difference)
;       NFUNC_EVALS          int(Nsolvers)        Number of ``lightning_mpfit_function.pro`` evaluations performed by MPFIT
;       CONVERGENCE_FLAG     int                  Flag indicating if any other flag was issued for a convergence metric
;       ================     ================     ==========================================================================================
;
; Modification History
; --------------------
;   - 2022/08/16: Created (Keith Doore)
;   - 2022/09/26: Added degrees of freedom to output structure (Keith Doore)
;   - 2023/01/23: Adjusted ``pvalue`` calculation if ``lnprob = NaN`` (Keith Doore)
;   - 2023/01/23: Removed ``similar_frac`` output (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(parameters) eq 0 then message, 'Variable is undefined: PARAMETERS.'
 if size(parameters, /type) lt 2 or size(parameters, /type) gt 5 then $
   message, 'PARAMETERS must be of type int, float, or double.'
 if size(parameters, /n_dim) lt 1 or size(parameters, /n_dim) gt 2 then $
   message, 'PARAMETERS must be a 1-D or 2-D array.'
 dim = size(parameters, /dimensions)
 Nparam = dim[0]
 if size(parameters, /n_dim) eq 1 then Nsolvers = 1 else Nsolvers = dim[1]

 if n_elements(lnprob) eq 0 then message, 'Variable is undefined: LNPROB.'
 if size(lnprob, /type) lt 2 or size(lnprob, /type) gt 5 then $
   message, 'LNPROB must be of type int, float, or double.'
 if size(reform(lnprob), /n_dim) ne 1 then message, 'LNPROB must be a scalar or 1-D array.'
 if n_elements(lnprob) ne Nsolvers then message, 'LNPROB must have Nsolvers number of elements.'

 if n_elements(status) eq 0 then message, 'Variable is undefined: STATUS.'
 if size(status, /type) lt 2 or size(status, /type) gt 5 then $
   message, 'STATUS must be of type int, float, or double.'
 if size(reform(status), /n_dim) ne 1 then message, 'STATUS must be a scalar or 1-D array.'
 if n_elements(status) ne Nsolvers then message, 'STATUS must have Nsolvers number of elements.'

 if n_elements(error_msg) eq 0 then message, 'Variable is undefined: ERROR_MSG.'
 if size(error_msg, /type) ne 7 then message, 'ERROR_MSG must be of type string.'
 if size(reform(error_msg), /n_dim) ne 1 then message, 'ERROR_MSG must be a scalar or 1-D array.'
 if n_elements(error_msg) ne Nsolvers then message, 'ERROR_MSG must have Nsolvers number of elements.'

 if n_elements(Niterations) eq 0 then message, 'Variable is undefined: NITERATIONS.'
 if size(Niterations, /type) lt 2 or size(Niterations, /type) gt 5 then $
   message, 'NITERATIONS must be of type int, float, or double.'
 if size(reform(Niterations), /n_dim) ne 1 then message, 'NITERATIONS must be a scalar or 1-D array.'
 if n_elements(Niterations) ne Nsolvers then message, 'NITERATIONS must have Nsolvers number of elements.'

 if n_elements(Nfunc_evals) eq 0 then message, 'Variable is undefined: NFUNC_EVALS.'
 if size(Nfunc_evals, /type) lt 2 or size(Nfunc_evals, /type) gt 5 then $
   message, 'NFUNC_EVALS must be of type int, float, or double.'
 if size(reform(Nfunc_evals), /n_dim) ne 1 then message, 'NFUNC_EVALS must be a scalar or 1-D array.'
 if n_elements(Nfunc_evals) ne Nsolvers then message, 'NFUNC_EVALS must have Nsolvers number of elements.'

 if n_elements(DoF) eq 0 then message, 'Variable is undefined: DOF.'
 if size(DoF, /type) lt 2 or size(DoF, /type) gt 5 then $
   message, 'DOF must be of type int, float, or double.'
 if size(DoF, /n_dim) ne 0 then message, 'DOF must be a scalar.'

 if n_elements(maxiter) ne 0 then begin
   if size(maxiter, /type) lt 2 or size(maxiter, /type) gt 5 then $
     message, 'MAXITER must be of type int, float, or double.'
   if size(maxiter, /n_dim) ne 0 then message, 'MAXITER must be a scalar.'
   if maxiter le 0 then message, 'MAXITER must be a positive value.'
 endif else maxiter = 200


; Check the status and create flag. All status values greater
;   than zero can represent successful MPFIT run. However,
;   status=5 likely failed to converge. We do not include this
;   in the status flag, rather the iteration flag.
 status_flag = intarr(Nsolvers)
 status_flag[where(status le 0, /null)] = 1

; Check the fraction of iterations used and create flag if max was used.
 iter_frac = Niterations/double(maxiter)
 iter_flag = intarr(Nsolvers)
 iter_flag[where(iter_frac eq 1, /null)] = 1

; Check the quality of the fits based on chisqr test
 pvalue = replicate(!values.d_NaN, Nsolvers)
 if n_elements(where(finite(lnprob), /null)) gt 0 then $
   pvalue[where(finite(lnprob), /null)] = 1 - chisqr_pdf(-2.d * lnprob[where(finite(lnprob), /null)], DoF)
 
; Check for convergence between solvers
 if Nsolvers gt 1 then begin
   ; Find solvers whose chisqr is within 4 of the overall most likely solver.
   ;   Solvers that did not get within 4 likely did not reach the global (or 
   ;   best local) minimum in chisqr space and got stuck in a local minimum.
   ;   The value of 4 is arbitrarily chosen.
   chi2_diff = (-2.d * lnprob) - min(-2.d * lnprob, min_loc)
   reached_min = where(chi2_diff lt 4, Natmin)

   ; Check the fraction of solvers that were near max lnprob and set flag if less than half were not
   stuck_frac = 1 - Natmin/double(Nsolvers)
   if stuck_frac gt 0.5 then stuck_flag = 1 else stuck_flag = 0

   ; Check for solvers that were near min chisqr, did they have similar solutions
   ;   We assume similar solutions are those with parameter differences less than 1%
   ;   from the best fit solver. It is possible that the best fit solver could have a 
   ;   parameter with a value of 0. If that occurs just use the difference instead.
   if Natmin gt 1 then begin
     zero_param = where(parameters[*, min_loc] eq 0, comp=nonzero_param, ncomp=Nnonzero, /null)
     param_per_diff = parameters[*, reached_min] - rebin(parameters[*, min_loc], Nparam, Natmin)
     if Nnonzero gt 0 then $
       param_per_diff[nonzero_param, *] = param_per_diff[nonzero_param, *] / $
                                          rebin(parameters[nonzero_param, min_loc[0]], Nnonzero, Natmin)
     Nsimilar = fix(total(fix(param_per_diff lt 0.01d), 2))
   endif else Nsimilar = replicate(0, Nparam)
   similar_flag = intarr(Nparam)
   similar_flag[where(Nsimilar/double(Natmin) ne 1, /null)] = 1

 endif else begin
   ; Can't compute convergence between solvers if only using one solver
   ;   So, set values to NaN and flags to 0.
   stuck_frac = !values.D_NaN
   stuck_flag = 0
   similar_flag = replicate(0, Nparam)
 endelse


; Compute overall convergence flag that checks if any other flag was issued for a convergence metric
 convergence_flag = fix(total([status_flag, iter_flag, stuck_flag, similar_flag]) gt 0)


; Compile convergence tests and flags into structure
 convergence_metric = {status: status,                     $
                       status_flag: status_flag,           $
                       error_msg: error_msg,               $
                       iter_frac: iter_frac,               $
                       iter_flag: iter_flag,               $
                       pvalue: pvalue,                     $
                       dof: DoF,                           $
                       stuck_frac: stuck_frac,             $
                       stuck_flag: stuck_flag,             $
                       similar_flag: similar_flag,         $
                       Nfunc_evals: Nfunc_evals,           $
                       convergence_flag: convergence_flag  $
                       }


 return, convergence_metric

end