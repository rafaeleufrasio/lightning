L2KEV_J07
=========

Name
----
L2KEV_J07

Purpose
-------
Calculates the 2 keV monochromatic luminosity of an AGN 
given its intrinsic monochromatic luminosity at 2500
Angstroms using the Just et al. (2007) relationship:

    :math:`\alpha_{ox} = -0.137 \log(L_{2500}) + 2.638`

for :math:`L_{2500}` in :math:`\rm ergs\ s^{-1}\ Hz^{-1}`, and the definition of
:math:`\alpha_{ox}`, the UV-to-X-ray slope:

    :math:`\alpha_{ox} = -0.3838 \log(L_{2500} / L_{\rm 2keV})`.

Calling Sequence
----------------
::

    L2keV = L2keV_J07(L2500 [, dalpha = , /error_check])

Input
-----
``L2500`` : float or double array(...)
    The rest-frame 2500 Angstrom monochromatic luminosity
    :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Optional Inputs
---------------
``dalpha`` : float or double array(...)
    Deviation from Just et al. (2007) relationship. A "fudge factor" added
    to the expected value of :math:`\alpha_{ox}`. (Default = ``0.0``)
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Outputs
-------
``L2keV`` : double array(...)
    The 2 keV monochromatic luminosity :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Reference
---------
`Just, D. W., Brandt, W. N., Shemmer, O., et al. 2007, ApJ, 665, 2 <https://ui.adsabs.harvard.edu/abs/2007ApJ...665.1004J/abstract>`_

Modification History
--------------------
- 2021/03/19: Created (Erik B. Monson)
- 2022/03/16: Moved to separate file and documentation improved (Erik B. Monson).
- 2022/06/20: Replaced ``!cv`` with ``!lightning_cgs`` (Keith Doore)
- 2022/06/20: Updated documentation (Keith Doore)
- 2022/06/20: Added error handling (Keith Doore)
- 2022/06/20: Added ``error_check`` keyword to do error handling (Keith Doore)
- 2022/06/20: Removed ``alpha_ox_L2500_J07`` function, since only one line of code (Keith Doore)

