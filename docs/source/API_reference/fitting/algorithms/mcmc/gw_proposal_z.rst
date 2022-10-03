GW_PROPOSAL_Z
=============

Name
----
GW_PROPOSAL_Z

Purpose
-------
Samples a specified number of values from the stretch move proposal distribution
from the Goodman & Weare (2010) affine invariant MCMC sampling algorithm.

Calling Sequence
----------------
::

    z = gw_proposal_z(Nparallel, a [, /error_check])

Inputs
------
``Nparallel`` : int, float, or double scalar
    Number of samples to draw (i.e., one sample for each ensemble walker).
``a`` : int, float, or double scalar
    A real constant :math:`\geq 1` which controls the size of the proposal 
    distribution (i.e., how close to 0 ``z`` is allowed to be and how 
    large ``z`` is allowed to be). In practice, we can only sample 
    values of ``z`` in ``[1/a, a]``.

Optional Input
--------------
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Output
------
``z`` : double array(Nparallel)
    The proposal scaling constant drawn from the proposal distribution.

Reference
---------
`Goodman, J., & Weare, J. 2010, CAMCS, 5, 65 <https://ui.adsabs.harvard.edu/abs/2010CAMCS...5...65G/abstract>`_

Modification History
--------------------
- 2021/12/02: Created (Erik B. Monson)
- 2022/04/21: Documentation (Erik B. Monson)
- 2022/06/17: Renamed variables to match naming scheme (Keith Doore)
- 2022/06/17: Updated documentation (Keith Doore)
- 2022/06/17: Added proper error handling (Keith Doore)
- 2022/06/17: Added ``error_check`` keyword to do error handling (Keith Doore)

