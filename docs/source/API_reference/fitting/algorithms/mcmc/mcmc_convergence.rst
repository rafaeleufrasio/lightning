MCMC_CONVERGENCE
================

Name
----
MCMC_CONVERGENCE

Purpose
-------
Tests the resulting MCMC chain for convergence. This includes the
autocorrelation time for both the affine-invariant and adaptive 
methods. Additionally, the Gelman-Rubin statistic as updated by
Brooks & Gelman (1998) is computed for the adaptive method.

Calling Sequence
----------------
::

    convergence_metric = mcmc_convergence(chain, lnprob_chain, accepted_trials [, C_step = , $
                                          tolerance = , burn_in = , thin_factor = , $
                                          final_chain_length = , method = ])

Input
-----
``chain`` : int, float, or double array(Nparam, Ntrials, Nparallel)
    MCMC chain or ensemble of MCMC chains.
``lnprob_chain`` : int, float, or double array(Ntrials, Nparallel)
    Log-probability of each element in the MCMC chain or ensemble
    of MCMC chains.
``accepted_trials`` : int, float, or double array(Nparallel)
    The number of accepted MCMC trials for each MCMC chain.

Optional Inputs
---------------
``C_step`` : int, float, or double scalar
    Defines how many trials of the chain are used to calculate tau, where
    we integrate tau to the smallest index ``M`` such that ``M > C_step * tau``.
    (Default = ``5``)
``tolerance`` : int, float, or double scalar
    Defines how many taus the length of the chain should be for us to 
    believe the estimate. (Default = ``50``)
``burn_in`` : int, float, or double scalar
    The number of initial MCMC trials to truncate as the burn-in phase. If set to ``0``,
    then the number will chosen automatically from the autocorrelation time.
    (Default = ``0``)
``thin_factor`` : int, float, or double scalar
    The factor to thin the MCMC chain after removing the burn-in trials. If set to ``0``,
    then the number will be chosen automatically from the autocorrelation time.
    (Default = ``0``)
``final_chain_length`` : int, float, or double scalar
    The number of MCMC trials to include for the final distributions as taken from 
    end of the chain (only for adaptive MCMC method). (Default = ``1000``)
``method`` : string scalar
  The fitting algorithm used to fit the SED(s). Current options are: 
  ``'MCMC-ADAPTIVE'`` and ``'MCMC-AFFINE'``. (Default = ``'MCMC-AFFINE'``)

Output
------
``convergence_metric`` : structure
    A structure containing the convergence metrics and flags used to 
    indicate if convergence was reached. (A flag of ``1`` indicates failure
    to converge.)
    The full description of the structure is as follows:

    ===================     =================     ==================================================================================
    TAG                     TYPE                  DESCRIPTION
    ===================     =================     ==================================================================================
    ACCEPTANCE_FRAC         double(Nparallel)     Fraction of accepted trials for each MCMC chain
    AUTOCORR_TIME           double(Nparam)        Autocorrelation time of each parameter for the ensemble of MCMC chains
    BURN_IN_AUTOCORR        int                   The number of burn-in trials as determined from the autocorrelation time
    GELMAN_RUBIN_R_HAT      double(Nparam)        The Gelman-Rubin statistic (only for adaptive MCMC, ``NaN`` otherwise)
    BROOKS_GELMAN_R_HAT     double                The Brooks-Gelman statistic (only for adaptive MCMC, ``NaN`` otherwise)
    NCHAIN_R_HAT            int                   Number of parallel chains with ``lnprob`` < 2 above best fit chain element
    ACCEPTANCE_FLAG         int(Nparallel)        Flag indicating if acceptance fraction is outside optimal range (``0.2-0.5``)
    AUTOCORR_FLAG           int(Nparam)           Flag indicating if autocorrelation time is above specified tolerance
    GELMAN_RUBIN_FLAG       int(Nparam)           Flag indicating if Gelman-Rubin stat is too large (``r_hat >= 1.2``)
    BROOKS_GELMAN_FLAG      int                   Flag indicating if Brooks-Gelman stat is too large (``r_hat >= 1.2``)
    CONVERGENCE_FLAG        int                   Flag indicating if any other flag was issued for a convergence metric
    ===================     =================     ==================================================================================

References
----------
- `Gelman, A., & Rubin, D. B. 1992, StaSc, 7, 457 <https://ui.adsabs.harvard.edu/abs/1992StaSc...7..457G/abstract>`_
- `Brooks, S. P., & Gelman, A. 1998, Journal of Computational and Graphical Statistics, 7, 434 <https://www.tandfonline.com/doi/abs/10.1080/10618600.1998.10474787>`_
- `Goodman, J., & Weare, J. 2010, CAMCS, 5, 65 <https://ui.adsabs.harvard.edu/abs/2010CAMCS...5...65G/abstract>`_

Modification History
--------------------
- 2022/07/14: Created (Keith Doore)
- 2022/08/03: Updated adaptive ``r_hat`` metric calculations (Keith Doore)
- 2022/08/03: Added optional ``burn_in`` input, and renamed automatic ``burn-in`` from autocorr to ``burn_in_autocorr`` (Keith Doore)
- 2022/08/03: Added optional ``thin_factor`` input to include full pre-thinned chain for MCMC-ADAPTIVE ``r_hat`` stat (Keith Doore)
- 2022/08/17: Fixed issue with Gelman-Rubin ``r_hat`` ``NaN`` array, when ``r_hat`` not calculated (Keith Doore)
- 2022/08/18: Fixed issue where pre-thinned adaptive chain length may be greater than ``Ntrials`` (Keith Doore)

