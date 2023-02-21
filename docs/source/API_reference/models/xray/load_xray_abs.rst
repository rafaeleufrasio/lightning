LOAD_XRAY_ABS
=============

Name
----
LOAD_XRAY_ABS

Purpose
-------
Loads the specified X-ray absorption model,and interpolates
it to the input wavelength grid.

Calling Sequence
----------------
::

    exp_neg_tau_xray = load_xray_abs(wave, xray_abs_model)

Inputs
------
``wave`` : int, float or double array(Nwave)
    A grid of wavelengths to be interpolated from the read in
    X-ray absorption model :math:`[\mu \rm m]`.
``xray_abs_model`` : string scalar
    The name of the X-ray absorption model to apply to the X-ray emission.
    Current options are ``'TBABS-WILM'``, ``'ATTEN'``,
    and ``'NONE'``.

Output
------
``exp_neg_tau_xray`` : float or double array(Nwave)
    The absorption in terms of :math:`e^{-\tau}` interpolated from the read
    in X-ray absorption model.

Modification History
--------------------
- 2021/09/21: Created (Erik B. Monson).
- 2022/06/21: Added documentation (Keith Doore)
- 2022/06/21: Converted to using ``config`` from ``curve`` keyword (Keith Doore)
- 2022/07/01: Converted from ``config`` to ``xray_abs_model`` keyword since it was the only accessed value from ``config`` (Keith Doore)

