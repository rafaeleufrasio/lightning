CALZETTI00_ATTEN
================

Name
----
CALZETTI00_ATTEN

Purpose
-------
Generates attenuation values at the input wavelengths using the 
Calzetti et al. (2000) attenuation curve or its modified version.
Besides the normalization in the standard Calzetti curve, the 
modified version includes a variable slope as described in 
Noll et al. (2009), an optional 2175 Angstrom bump feature specified
in Kriek & Conroy (2013), and birth cloud attenuation as described in
Eufrasio et al. (2017).

Calling Sequence
----------------
::

    exp_neg_tau = calzetti00_atten(wave [, tauV_diff = , delta = , $
                                   tauV_BC = , /uv_bump, /error_check])

Input
-----
``wave`` : int, float, or double array(Nwave)
    The wavelength at which to determine the attenuation :math:`[\mu \rm m]`.

Optional Inputs
---------------
``tauV_diff`` : int, float, or double array(Nmodels)
    The V-band optical depth of the diffuse dust. (Default = ``1.d0``)
``delta`` : int, float, or double array(Nmodels)
    The power law value to change the attenuation curve slope.
    (Default = ``0.d0``)
``tauV_BC`` : int, float, or double array(Nmodels)
    The V-band optical depth of the birth cloud. (Default = ``0.d0``)
``uv_bump`` : flag
    If set, then a 2175 Angstrom UV bump feature will be added to the 
    attenuation curve.
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Output
------
``exp_neg_tau`` : double array(Nwave, Nmodels)
    The attenuation in terms of :math:`e^{-\tau}` (i.e., :math:`e^{-\tau} = 10^{-0.4 A_{\lambda}}`).

Note
----
Defaults will only be set if the optional ``error_check`` input is set.

References
----------
- `Calzetti, D., Armus, L., Bohlin, R. C., et al. 2000, ApJ, 533, 682 <https://ui.adsabs.harvard.edu/abs/2000ApJ...533..682C/abstract>`_
- `Noll, S., Burgarella, D., Giovannoli, E., et al. 2009, A&A, 507, 1793 <https://ui.adsabs.harvard.edu/abs/2009A%26A...507.1793N/abstract>`_
- `Kriek, M., & Conroy, C. 2013, ApJ, 775, L16 <https://ui.adsabs.harvard.edu/abs/2013ApJ...775L..16K/abstract>`_
- `Eufrasio, R. T., Lehmer, B. D., Zezas, A., et al. 2017, ApJ, 851, 10 <https://ui.adsabs.harvard.edu/abs/2017ApJ...851...10E/abstract>`_

Modification History
--------------------
- 2016/05/01: Created (Rafael T. Eufrasio)
- 2020/05/06: Added ability to run with pure Calzetti curve (Keith Doore)
- 2022/03/15: Added proper error handling (Keith Doore)
- 2022/03/15: Renamed variables to standard format (Keith Doore)
- 2022/03/15: Cleaned up method for separating pure and modified Calzetti (Keith Doore)
- 2022/03/18: Updated documentation (Keith Doore)
- 2022/04/07: Allowed for inputs to have degenerate dimensions (Keith Doore)
- 2022/04/07: Allowed for inputs to be scalars (Keith Doore)
- 2022/04/07: Allowed integer inputs (Keith Doore)
- 2022/04/07: Made optical depths positive input values (Keith Doore)
- 2022/04/07: Allowed for array of only one optional input to be given (Keith Doore)
- 2022/06/09: Changed ``no_bump`` keyword to ``uv_bump`` keyword and adjusted logic accordingly (Keith Doore)
- 2022/06/09: Added ``error_check`` keyword to do error handling (Keith Doore)
- 2022/11/04: Corrected extrapolation below 0.12um to use method from Noll+2009 (Keith Doore)

