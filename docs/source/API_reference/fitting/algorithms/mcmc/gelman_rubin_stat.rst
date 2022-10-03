GELMAN_RUBIN_STAT
=================

Name
----
GELMAN_RUBIN_STAT

Purpose
-------
Computes the Gelman-Rubin statistic (:math:`\hat{R}`) of multiple Markov chains.
Uses the within-chain variance and between-chain variance to compare
convergence of the chains for each individual parameter. The use of 
degrees of freedom has been updated to the corrected version from 
Brooks & Gelman (1998).

Calling Sequence
----------------
::

    r_hat = gelman_rubin_stat(chain)

Input
-----
``chain`` : int, float, or double array(Nparam, Ntrials, Nparallel)
    An ensemble of Markov chains.

Output
------
``r_hat`` : double array(Nparam)
    The Gelman-Rubin statistic of the chains. The :math:`\sqrt{\hat{R}}`
    should be less than or equal to 1.2 if convergence of the
    chain was sufficiently reached. If a parameter is constant, the
    corresponding Gelman-Rubin stat is ``NaN``.

References
----------
- `Gelman, A., & Rubin, D. B. 1992, StaSc, 7, 457 <https://ui.adsabs.harvard.edu/abs/1992StaSc...7..457G/abstract>`_
- `Brooks, S. P., & Gelman, A. 1998, Journal of Computational and Graphical Statistics, 7, 434 <https://www.tandfonline.com/doi/abs/10.1080/10618600.1998.10474787>`_

Modification History
--------------------
- 2020/06/10: Created (Keith Doore)
- 2022/07/14: Updated documentation and error handling (Keith Doore)

