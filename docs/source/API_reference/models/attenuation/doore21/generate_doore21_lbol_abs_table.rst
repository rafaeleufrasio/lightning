GENERATE_DOORE21_LBOL_ABS_TABLE
===============================

Name
----
GENERATE_DOORE21_LBOL_ABS_TABLE

Purpose
-------
Generates and saves pre-computed table(s) containing the bolometric luminosity
of the stellar emission absorbed by dust (:math:`L_{\rm bol}^{\rm abs}``) for
the Doore et al. (2021) attenuation curves to have energy conservation. Since
the Doore et al. (2021) attenuation curves are inclination dependent, 
:math:`L_{\rm bol}^{\rm abs}`` cannot be simply computed as the difference
between the attenuated and unattenuated bolometric luminosities. Rather, the
attenuated bolometric luminosities must first be integrated over inclination.
Therefore, to save computational time, this process is performed before fitting
rather than during (see Section 4.4 of Doore et al. 2021 for further details).

Calling Sequence
----------------
::

    generate_doore21_lbol_abs_table, input_dir , redshifts, config

Inputs
------
``input_dir`` : string scalar
    The path to the file containing the input SED data.
``redshifts`` : int, float, or double array(Nred)
    The redshifts at which to compute the tables.
``config`` : structure
    A Lightning configuration structure. (See
    ``lightning_configure_defaults.pro`` for details and contents.)

Output
------
A FITS file per redshift containing a structure that has the pre-computed
table of :math:`L_{\rm bol}^{\rm abs}`` to have energy conservation with
the Doore et al. (2021) attenuation curves. Additionally, the parameter
grids and ``step_bounds`` are included in the structure. Files are saved
in the directory ``<input_dir>/lightning_output/doore21_Lbol_abs_table/``.
The full description of the structure is as follows:

==============     =============================     ========================================================
TAG                TYPE                              DESCRIPTION
==============     =============================     ========================================================
LBOL_ABS_MODEL     double(Nsteps, Ntaubf, Nn, 2)     Absorbed bolometric stellar luminosity :math:`[L_\odot]`
STEPS_BOUNDS       double(Nsteps+1)                  Age bounds :math:`[\rm yr]`
NSTEPS             long                              Number of age bins (steps)
TAUB_F_GRID        double(Ntaubf)                    Gridded values of ``tauB_f``
RDISK0_GRID        double(Nn)                        Gridded values of ``rdisk0``
F_GRID             double(Nn)                        Gridded values of ``F_clump``
==============     =============================     ========================================================

Notes
-----
The last dimension in the ``LBOL_ABS_MODEL`` tag indicates the young (``0``) and
old populations (``1``), respectively. Since the young population ignores ``rdisk0``
and the old population ignores ``F_clump``, this minimizes memory used to save
the files, while still encapsulating all parameter space.

Reference
---------
`Doore, K., Eufrasio, R. T., Lehmer, B. D., et al. 2021, ApJ, 923, 26 <https://ui.adsabs.harvard.edu/abs/2021ApJ...923...26D/abstract>`_

Modification History
--------------------
- 2022/04/19: Created (Keith Doore)
- 2022/05/17: Changed to only allow for specified redshift array (Keith Doore)
- 2022/05/17: Changed ``redshift`` and ``steps_bounds`` to be required inputs (Keith Doore)
- 2022/05/17: Replaced ``steps_bounds`` with configuration structure (Keith Doore)
- 2022/05/18: Allowed for unique cosmologies (Keith Doore)
- 2022/06/29: Replaced ``!cv`` with ``!lightning_cgs`` (Keith Doore)
- 2022/06/29: Fixed error where ``steps_bounds`` was not retrieved from ``config`` (Keith Doore)
- 2022/07/08: Fixed issue where age bins were not truncated to age of universe at ``z=0`` (Keith Doore)
- 2022/08/18: Added progress printing (Keith Doore)

