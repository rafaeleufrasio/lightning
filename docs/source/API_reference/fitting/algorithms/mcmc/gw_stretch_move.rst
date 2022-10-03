GW_STRETCH_MOVE
===============

Name
----
GW_STRETCH_MOVE

Purpose
-------
Proposes a new step for an MCMC ensemble following the "stretch move"
prescription in Goodman & Weare (2010).

Calling Sequence
----------------
::

    ensemble_new = gw_stretch_move(ensemble, a [, /error_check, z=z])

Inputs
------
``ensemble`` : float or double array(Nparam, Nparallel)
    The current state of all the chains in the MCMC ensemble at a given trial.
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
``ensemble_new`` : double array(Nparam, Nparallel)
    The new proposed positions of the ensemble.

Optional Output
---------------
``z`` : double array(Nparallel)
    The proposal scaling constant.

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
- 2022/06/17: Removed ``bounds_arr`` input and subsequent code as this is now done in prior (Keith Doore)
- 2022/06/17: Removed ``seed`` input as we do not ever specify the random seed (Keith Doore)
- 2022/06/17: Made ``z`` a separate output from ``ensemble_new`` to prevent any issues or confusion (Keith Doore)
- 2022/07/05: Changed from ``round`` to ``floor`` in ``kprime`` as index could exceed array size (Keith Doore)

