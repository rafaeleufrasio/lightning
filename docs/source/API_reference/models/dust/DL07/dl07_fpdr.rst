DL07_FPDR
=========

Name
----
DL07_FPDR

Purpose
-------
Computes fPDR (fraction of total luminosity from the photodissociation
regions) according to Eq. 4 of Aniano et al. (2020). This calculation
assumes the parameters are defined as in Draine & Li (2007).

Calling Sequence
----------------
::

    fPDR = dl07_fpdr(alpha, umin, umax, gam [, /error_check])

Inputs
------
``alpha`` : int, float, or double array(Nmodels)
    The exponent of the power-law distribution of heating starlight
    intensities between ``umin`` and ``umax``.
``umin`` : int, float, or double array(Nmodels)
    The minimum radiation field intensity of the diffuse ISM radiation
    field heating the dust.
``umax`` : int, float, or double array(Nmodels)
    The maximum radiation field intensity of the power-law distribution
    of heating starlight intensities.
``gam`` : int, float, or double array(Nmodels)
    The fraction of the dust mass exposed to the power-law distribution of
    starlight intensities.

Optional Input
--------------
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Output
------
``fPDR`` : double array(Nmodel, ...)
    The fraction of the total dust luminosity that is radiated by
    dust in regions where ``U > 1e2`` (i.e., fraction of total luminosity
    from the photodissociation regions).

References
----------
- `Draine, B. T., & Li, A. 2007, ApJ, 657, 810 <https://ui.adsabs.harvard.edu/abs/2007ApJ...657..810D/abstract>`_
- `Aniano, G., Draine, B. T., Hunt, L. K., et al. 2020, ApJ, 889, 150 <https://ui.adsabs.harvard.edu/abs/2020ApJ...889..150A/abstract>`_

Modification History
--------------------
- 2020/08/19: Created (Rafael T. Eufrasio)
- 2022/03/15: Added error handling (Keith Doore)
- 2022/03/18: Updated documentation (Keith Doore)
- 2022/04/12: Allowed for inputs to have degenerate dimensions (Keith Doore)
- 2022/04/12: Allowed integer inputs (Keith Doore)
- 2022/07/12: Renamed reformed ``gam`` input as to not change input value (Keith Doore)

