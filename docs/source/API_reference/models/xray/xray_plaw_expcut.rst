XRAY_PLAW_EXPCUT
================

Name
----
XRAY_PLAW_EXPCUT

Purpose
-------
Generates an X-ray power law spectrum with a high energy exponential 
cutoff (i.e., :math:`L_{\nu} \propto F_{\nu} \propto E^{-\gamma + 1}
e^{-E / E_{\rm cut}}`). The spectrum is normalized by default to the
monochromatic flux at 2 keV, computed in a 0.5 keV (~100 Angstrom) wide
box filter.

Calling Sequence
----------------
::

    plaw_spec = xray_plaw_expcut(wave [, plaw_gamma = , E_cut = , normbox = ])

Input
-----
``wave`` : int, float or double array(Nwave)
    A grid of wavelengths at which evaluate the power law spectrum
    :math:`[\mu \rm m]`.

Optional Inputs
---------------
``plaw_gamma`` : int, float or double scalar
    The photon index of the power law. (Default = ``1.8``)
``E_cut`` : int, float or double scalar
    The high energy exponential cutoff value :math:`[{\rm keV}]`. (Default = ``300``)
``normbox`` : int, float or double array(2)
    The lower and upper energies of a box in which to normalize
    the spectrum :math:`[{\rm keV}]`. (Default = ``[1.75, 2.25]``)

Output
------
``plaw_spec`` : double array(Nwave)
    The normalized power law spectrum.

Modification History
--------------------
- 2021/03/19: Created (Erik B. Monson).
- 2022/06/22: Added documentation (Keith Doore)
- 2022/06/22: Replaced ``!cv`` with ``!lightning_cgs`` (Keith Doore)
- 2022/06/22: Updated documentation and renamed variables (Keith Doore)

