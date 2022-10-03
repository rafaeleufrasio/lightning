DOORE21_INTERP_LBOL_ABS_TABLE
=============================

Name
----
DOORE21_INTERP_LBOL_ABS_TABLE

Purpose
-------
Linearly Interpolates the bolometric luminosity of the attenuated stellar
emission from precomputed tables for the Doore et al. (2021) attenuations curves.

Calling Sequence
----------------
::

    steps_Lbol_abs = doore21_interp_lbol_abs_table(tauB_f, F_clump, b_to_d, $
                          rold0_ages, doore21_Lbol_abs_table [, /error_check])

Inputs
------
``tauB_f`` : int, float, or double array(Nmodels)
    The face-on optical depth in the B-band.
``F_clump`` : int, float, or double array(Nmodels)
    The clumpiness factor F.
``b_to_d`` : int, float, or double array(Nmodels)
    The bulge-to-disk ratio.
``rold0_ages`` : int, float, or double array(Nsteps)
    The binary parameter ``rold0``, designating each SFH age bin as part of
    the young or old population when using the Doore et al. (2021) attenuation
    curves. A value of ``0`` for the corresponding age bin considers it to be part of
    the young population, and a value of ``1`` considers it to be part of the old
    populations (see Doore et al. 2021 for more details).
``doore21_Lbol_abs_table`` : structure
    A structure giving the pre-computed table of ``steps_Lbol_abs`` for a given redshift
    to have energy conservation with the Doore et al. (2021) attenuation curves.
    (See ``generate_doore21_lbol_abs_table.pro`` for details and contents.)

Optional Input
--------------
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Output
------
``steps_Lbol_abs`` : double array(Nsteps, Nmodels)
    The interpolated bolometric luminosity of the attenuated stellar emission for
    each SFH step :math:`[L_\odot]`.

Reference
---------
`Doore, K., Eufrasio, R. T., Lehmer, B. D., et al. 2021, ApJ, 923, 26 <https://ui.adsabs.harvard.edu/abs/2021ApJ...923...26D/abstract>`_

Modification History
--------------------
- 2022/04/15: Created (Keith Doore)
- 2022/04/18: Added error handling (Keith Doore)
- 2022/06/09: Added ``error_check`` keyword to do error handling (Keith Doore)

