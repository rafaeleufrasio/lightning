QSOSED_MODELS
=============

Name
----
QSOSED_MODELS

Purpose
-------
Generates the spectra and X-ray parameters for the QSOSED model as
described in Kubota & Done (2018). The spectra can include or not
include Galactic absorption.

Calling Sequence
----------------
::

    model_agn_count_rate = qsosed_models(redshift, wave_rest, arf_interp, exp_neg_tau_xray_MW, $
                                galactic_nH, lumin_dist [, Lnu_x_agn=Lnu_x_agn, $
                                L2500=L2500, agn_mass=agn_mass, agn_logmdot=agn_logmdot])

Inputs
------
``redshift`` : int, float, or double scalar
    The redshift of the model.
``wave_rest`` : int, float or double array(Nwave)
    The rest-frame wavelengths at which to interpolate the qsosed models
    :math:`[\mu \rm m]`.
``arf_interp`` : double array(Nwave)
    A grid of ARF values interpolated from the input ARF file :math:`[{\rm cm}^2]`.
``exp_neg_tau_xray_MW`` : float or double array(Nwave)
    The attenuation from the Milky Way to be applied to the X-ray data in terms
    of :math:`e^{-\tau}`.
``galactic_nH`` : int, float, or double scalar
    Galactic, i.e. Milky Way, neutral Hydrogen column density along the line
    of sight :math:`[10^{20}\ {\rm cm}^2]`.
``lumin_dist`` : int, float, double scalar
    The luminosity distance of the model :math:`[{\rm Mpc}]`.

Output
------
``model_agn_count_rate`` : double array(Nwave, Nmass, Nmdot)
    Instrumental count rate density produced by qsosed model grid
    :math:`[{\rm counts\ s^{-1}\ Hz^{-1}}]`.

Optional Outputs
----------------
``Lnu_x_agn`` : double array(Nwave, Nmass, Nmdot)
    The qsosed model luminosity :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``L2500`` : double array(Nmass, Nmdot)
    The rest-frame 2500 Angstrom monochromatic luminosity shifted to the
    observed frame :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``agn_mass`` : double array(Nmass, Nmdot)
    The SMBH mass grid :math:`[M_\odot]`.
``agn_logmdot`` : double array(Nmass, Nmdot)
    The Log10 of SMBH accretion rate grid, normalized by the Eddington rate.

Reference
---------
`Kubota, A., & Done, C. 2018, MNRAS, 480, 1247 <https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.1247K/abstract>`_

Modification History
--------------------
- 2021/09/21: Created (Erik B. Monson).
- 2022/06/22: Moved to separate file and documentation improved (Keith Doore).
- 2022/06/22: Major update to include new implementation (Keith Doore)
- 2022/11/02: Galactic NH is now in units of 1e20 cm-2 (Erik B. Monson)

