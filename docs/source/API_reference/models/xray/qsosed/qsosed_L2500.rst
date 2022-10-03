QSOSED_L2500
============

Name
----
QSOSED_L2500

Purpose
-------
Generates the qsosed :math:`L_{2500}` (monochromatic luminosity at 2500 Angstroms)
corresponding to given values of :math:`M` and :math:`\log(\dot{M})` using the
models from Kubota & Done (2018).

Calling Sequence
----------------
::

    L2500 = qsosed_L2500(xray_models, agn_mass, agn_logmdot [, /error_check])

Inputs
------
``xray_models`` : structure
    A structure containing the spectra, counts, and X-ray model parameters.
    (See ``xrb_xagn_models.pro`` for details and contents.)
``agn_mass`` : int, float, or double array(Nmodels)
    Supermassive black hole mass :math:`[M_\odot]`.
``agn_logmdot`` : int, float, or double array(Nmodels)
    Log10 of SMBH accretion rate, normalized by the Eddington rate.

Optional Input
--------------
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Output
------
``L2500`` : double array(Nmodels)
    The interpolated rest-frame 2500 Angstrom monochromatic luminosity of 
    shifted to the observed frame :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Reference
---------
`Kubota, A., & Done, C. 2018, MNRAS, 480, 1247 <https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.1247K/abstract>`_

Modification History
--------------------
- 2022/04/18: Created (Erik B. Monson)
- 2022/06/10: Added error handling (Keith Doore)
- 2022/06/10: Added ``error_check`` keyword to do error handling (Keith Doore)
- 2022/06/10: Made ``mass`` and ``logmdot`` required inputs (Keith Doore)

