MPFIT_CONVERGENCE
=================

Name
----
MPFIT_CONVERGENCE

Purpose
-------
Tests the resulting MPFIT solvers for convergence. This includes
checking the MPFIT status, the fraction of max iterations used by
each solver to reach its solution, the p-value of each solver 
from a simple :math:`\chi^2` test, how many solvers got stuck in
low probability regions, and how similar are the solver solutions
to each other.

Calling Sequence
----------------
::

    convergence_metric = mpfit_convergence(parameters, lnprob, status, error_msg, $
                                           Niterations, Nfunc_evals, $
                                           DoF [, maxiter = ])

Inputs
------
``parameters`` : int, float, or double array(Nparam, Nsolvers)
    The best fit parameters of the model as determined by MPFIT for each
    solver. The actual parameters contained in this array depend on the
    chosen model during configuration.
``lnprob`` : int, float, or double array(Nsolvers)
    Log-probability for each set of parameters.
``status`` : int array(Nsolvers)
    The status code given by MPFIT.
``error_msg`` : string array(Nsolvers)
    An error or warning message given by MPFIT.
``Niterations`` : int array(Nsolvers)
    The number of MPFIT iterations completed.
``Nfunc_evals`` : int array(Nsolvers)
    The number of ``lightning_mpfit_function.pro`` evaluations performed.
``DoF`` : int
    The number of degrees of freedom in the fit.

Optional Input
--------------
``maxiter`` : int, float, or double scalar
    The maximum number of MPFIT iterations that could have been perform.
    (Default = ``200``)

Output
------
``convergence_metric`` : structure
    A structure containing the convergence metrics and flags used to 
    indicate if convergence was reached. (A flag of ``1`` indicates failure
    to converge.)
    The full description of the structure is as follows:

    ================     ================     ==========================================================================================
    TAG                  TYPE                 DESCRIPTION
    ================     ================     ==========================================================================================
    STATUS               int(Nsolvers)        The status code given by MPFIT
    STATUS_FLAG          int(Nsolvers)        Flag indicating if the MPFIT status indicated failure of the algorithm (``<= 0``)
    ERROR_MSG            string(Nsolvers)     Error or warning message given by MPFIT. Will be blank if no message
    ITER_FRAC            double(Nsolvers)     Fraction of max iterations used by MPFIT to reach solution
    ITER_FLAG            int(Nsolvers)        Flag indicating if the maximum number of iteration were used by MPFIT
    PVALUE               double(Nsolvers)     P-value of fit as determined from a :math:`\chi^2` test using ``lnprob`` and ``DoF``
    DOF                  int                  Degrees of freedom (same for each solver)
    STUCK_FRAC           double               Fraction of solvers with ``lnprob`` > 2 above best fit solver
    STUCK_FLAG           int                  Flag indicating if the majority of solvers had ``lnprob`` > 2 above best fit solver
    SIMILAR_FLAG         int(Nparam)          Flag indicating if any non-stuck solvers had different solutions (>1% difference)
    NFUNC_EVALS          int(Nsolvers)        Number of ``lightning_mpfit_function.pro`` evaluations performed by MPFIT
    CONVERGENCE_FLAG     int                  Flag indicating if any other flag was issued for a convergence metric
    ================     ================     ==========================================================================================

Modification History
--------------------
- 2022/08/16: Created (Keith Doore)
- 2022/09/26: Added degrees of freedom to output structure (Keith Doore)
- 2023/01/23: Adjusted ``pvalue`` calculation if ``lnprob = NaN`` (Keith Doore)
- 2023/01/23: Removed ``similar_frac`` output (Keith Doore)

