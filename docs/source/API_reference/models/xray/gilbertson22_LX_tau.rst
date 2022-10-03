GILBERTSON22_LX_TAU
===================

Name
----
GILBERTSON22_LX_TAU

Purpose
-------
Produces the 2-10 keV LMXB and HMXB luminosities appropriate for the stellar 
population(s), according to the :math:`L_X/M` parametrizations with stellar age 
(that is, :math:`\tau = \log(\rm age)`) from Gilbertson et al. (2022).

Calling Sequence
----------------
::

    Lnu_XRB = gilbertson22_LX_tau(stellar_models, psi [, /error_check, $
                     Lnu_LMXB_steps=Lnu_LMXB_steps, Lnu_HMXB_steps=Lnu_HMXB_steps])

Inputs
------
``stellar_models`` : structure
    A structure containing the spectra, SEDs, and stellar parameters for the 
    non-parametric stellar model. (See ``binned_stellar_models.pro`` for
    details and contents.)
``psi`` : int, float, or double array(Nsteps, Nmodels)
    The non-parametric SFH coefficients :math:`[M_\odot\ {\rm yr}^{-1}]`.

Optional Input
--------------
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Output
------
``Lnu_XRB`` : double array(2, Nmodels) 
    The total LMXB (first column) and HMXB (second column) 2-10 keV
    luminosities :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Optional Outputs
----------------
``Lnu_LMXB_steps`` : double array(Nsteps, Nmodels)
    2-10 keV Lnu from LMXBs for each stellar age bin :math:`[L_\odot\ {\rm Hz}^{-1}]`.
``Lnu_HMXB_steps`` : double array(Nsteps, Nmodels)
    2-10 keV Lnu from HMXBs for each stellar age bin :math:`[L_\odot\ {\rm Hz}^{-1}]`.

Notes
-----
The parametrizations of :math:`\gamma = L_X / M` as a function of stellar age are as follows:

    :math:`\gamma(t_{\rm age}) = \gamma_{\rm HMXB}(t_{\rm age}) + \gamma_{\rm LMXB}(t_{\rm age})
    = [a_0(\tau - b_0)^2 + c_0] + [a_1(\tau - b_1)^2 + c_1]`

where :math:`\tau = \log_{10}(t_{\rm age})`. Gilbertson et al. (2022) has:

    - :math:`a_0 = -0.24 \ [{\rm erg s^{-1}\ M_\odot^{-1}}\ \log({\rm yr})^{-1}]`
    - :math:`b_0 = 5.23 \ [\log({\rm yr})^{-1}]`
    - :math:`c_0 = 32.54 \ [{\rm erg s^{-1}\ M_\odot^{-1}}]`
    - :math:`a_1 = -1.21 \ [{\rm erg s^{-1}\ M_\odot^{-1}}\ \log({\rm yr})^{-1}]`
    - :math:`b_1 = 9.32 \ [\log({\rm yr})^{-1}]`
    - :math:`c_1 = 29.09 \ [{\rm erg s^{-1}\ M_\odot^{-1}}]`

Reference
---------
`Gilbertson, W., Lehmer, B. D., Doore, K., et al. 2022, ApJ, 926, 28 <https://ui.adsabs.harvard.edu/abs/2022ApJ...926...28G/abstract>`_

Modification History
--------------------
- 2022/04/21: Created (Erik B. Monson)
- 2022/06/20: Renamed variables to common naming scheme (Keith Doore)
- 2022/06/20: Replaced ``!cv`` with ``!lightning_cgs`` (Keith Doore)
- 2022/06/20: Updated documentation (Keith Doore)
- 2022/06/20: Added error handling (Keith Doore)
- 2022/06/20: Added ``error_check`` keyword to do error handling (Keith Doore)
- 2022/07/12: Removed reforming of ``psi`` if ``Nmodels=1`` as it is not necessary (Keith Doore)

