DOORE21_ATTEN
=============

Name
----
DOORE21_ATTEN

Purpose
-------
Generates attenuation values at the input wavelengths using the attenuation
curves describe in Doore et al. (2021) given the attenuation curve
parameters. The attenuation curves are based on the Tuffs et al. (2004)
attenuation curves as updated by Popescu et al. (2011), which include the
parameters of inclination, face-on optical depth in the B-band, the
clumpiness factor, the fraction of intrinsic flux density from the old
stellar components compared to the total intrinsic flux density, and the
bulge-to-disk ratio.

Calling Sequence
----------------
::

    exp_neg_tau = doore21_atten(wave [, tauB_f = , F_clump = , cosi = , $
                                rold0 = , b_to_d = , tuffs_coeff = , /error_check])

Input
-----
``wave`` : int, float, or double array(Nwave)
    The wavelength at which to determine the attenuation :math:`[\mu \rm m]`.

Optional Inputs
---------------
``tauB_f`` : int, float, or double array(Nmodels)
    The face-on optical depth in the B-band. (Default = ``1.d0``)
``F_clump`` : int, float, or double array(Nmodels)
    The clumpiness factor F. (Default = ``0.d0``)
``cosi`` : int, float, or double array(Nmodels)
    The inclination of the galactic disk in terms of cos(i). (Default = ``1.d0``)
``rold0`` : int, float, or double array(Nmodels)
    The fraction of intrinsic flux density from the old stellar components
    compared to the total intrinsic flux density. (Default = ``0.d0``)
``b_to_d`` : int, float, or double array(Nmodels)
    The bulge-to-disk ratio. (Default = ``0.d0``)
``tuffs_coeff`` : structure
    A structure containing the polynomial coefficients needed for the 
    Tuffs model.
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
- `Tuffs, R. J., Popescu, C. C., Völk, H. J., Kylafis, N. D., & Dopita, M. A. 2004, A&A, 419, 821 <https://ui.adsabs.harvard.edu/abs/2004A%26A...419..821T/abstract>`_
- `Popescu, C. C., Tuffs, R. J., Dopita, M. A., et al. 2011, A&A, 527, A109 <https://ui.adsabs.harvard.edu/abs/2011A%26A...527A.109P/abstract>`_
- `Doore, K., Eufrasio, R. T., Lehmer, B. D., et al. 2021, ApJ, 923, 26 <https://ui.adsabs.harvard.edu/abs/2021ApJ...923...26D/abstract>`_

Modification History
--------------------
- 2019/09/22: Created (Keith Doore)
- 2021/04/13: Modified to have ``rdisk`` as intrinsic property (``rold0``) and ``rbulge`` as B/D (Keith Doore)
- 2022/03/15: Added proper error handling (Keith Doore)
- 2022/03/15: Renamed variables to standard format (Keith Doore)
- 2022/03/18: Updated Documentation (Keith Doore)
- 2022/03/31: Updated for efficiency, minimal arrays in memory (Keith Doore)
- 2022/04/07: Allowed for inputs to have degenerate dimensions (Keith Doore)
- 2022/04/07: Allowed for inputs to be scalars (Keith Doore)
- 2022/04/07: Allowed integer inputs (Keith Doore)
- 2022/04/07: Allowed for array of only one optional input to be given (Keith Doore)
- 2022/04/14: Allowed for ``b_to_d`` to be ``Inifinty`` and not cause ``NaNs`` in ``rbulge0`` (Keith Doore)
- 2022/04/15: Updated interpolation function for improved speed (Keith Doore)
- 2022/04/28: Removed for loop on 1-cosi polynomial for improved speed (Keith Doore)
- 2022/04/28: Removed unnecessary zeroing of arrays (Keith Doore)
- 2022/05/18: Turned ``lightning_dir`` into system variable call (Keith Doore)
- 2022/05/18: Added ``error_check`` keyword to do error handling (Keith Doore)
- 2022/05/18: Updated path to model files (Keith Doore)
- 2022/06/09: Made the coefficients an optional input for improved speed (Keith Doore)

