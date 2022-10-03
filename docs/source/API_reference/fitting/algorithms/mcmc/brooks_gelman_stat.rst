BROOKS_GELMAN_STAT
==================

Name
----
BROOKS_GELMAN_STAT

Purpose
-------
Computes the Brooks-Gelman statistic (i.e., the multivariate
Gelman-Rubin statistic as given in Brooks & Gelman 1998) of 
multiple Markov chains. Uses the within-chain variance and 
between-chain variance for all parameters to compare convergence
of the chains.

Calling Sequence
----------------
::

    r_hat = brooks_gelman_stat(chain)

Input
-----
``chain`` : int, float, or double array(Nparam, Ntrials, Nparallel)
    An ensemble of Markov chains.

Output
------
``r_hat`` : double scalar
    The Brooks-Gelman statistic of the chains. The :math:`\sqrt{\hat{R}}`
    should be less than or equal to 1.2 if convergence of the
    chain was sufficiently reached.

References
----------
- `Gelman, A., & Rubin, D. B. 1992, StaSc, 7, 457 <https://ui.adsabs.harvard.edu/abs/1992StaSc...7..457G/abstract>`_
- `Brooks, S. P., & Gelman, A. 1998, Journal of Computational and Graphical Statistics, 7, 434 <https://www.tandfonline.com/doi/abs/10.1080/10618600.1998.10474787>`_

Modification History
--------------------
- 2020/06/10: Created (Keith Doore)
- 2022/07/14: Updated documentation and error handling (Keith Doore)

