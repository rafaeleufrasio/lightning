X-ray Model Outputs
===================

When fitting an X-ray model, the post-processed outputs will vary depending on your
input X-ray data type (i.e., counts or flux) and if you are fitting an X-ray AGN model.
Below, we describe the outputs that are always included regardless of of your
configuration choices. Then, we describe the outputs unique to each X-ray data
type and the X-ray AGN model.

.. note::

    We define additional array size variables here for convenience. 
    (See :ref:`here <array-size-label>` for the already defined array size variables.)

    - ``Nxray`` : the maximum number of X-ray bandpasses for an SED included in the input


``XRAY_BANDPASS`` : double array(2, Nxray)
    The X-ray bandpass(es) specified in the input. The first index of first dimension contains 
    the lower energy bound, and the second index of first dimension contains the upper bound :math:`[\rm{keV}]`.

``GALACTIC_NH`` : double
    The Galactic (i.e. Milky Way) HI column density along the line of sight specified in the input
    :math:`[10^{20}\ \rm{cm}^{-2}]`.

``LNU_XRAYMOD`` : double(Nxray, Nmodels)
    The mean :math:`L_\nu` produced by the model in each bandpass :math:`[L_\odot\ {\rm Hz}^{-1}]`.
    In the case of the MPFIT algorithm, this is the best-fitting mean :math:`L_\nu` in each bandpass.
    In the case of the MCMC algorithms, this is the posterior distribution on the mean :math:`L_\nu` in each bandpass.

``WAVE_XRAYMOD_HIRES`` : double(200)
    The wavelength grid for the high resolution X-ray model :math:`[\mu \rm m]`.

``LNU_XRAYMOD_HIRES`` : double(200, Nhighres_models)
    The total high resolution X-ray luminosities produced by the model :math:`[L_\odot\ {\rm Hz}^{-1}]`.

``LNU_XRAYMOD_STAR_HIRES`` : double(200, Nhighres_models)
    The stellar population (i.e., LMXB and HMXB) contribution to the total high resolution X-ray luminosities,
    after absorption :math:`[L_\odot\ {\rm Hz}^{-1}]`.

``LNU_XRAYMOD_STAR_UNABS_HIRES`` : double(200, Nhighres_models)
    The intrinsic stellar population (i.e., LMXB and HMXB) contribution to the total high resolution 
    X-ray luminosities, without absorption :math:`[L_\odot\ {\rm Hz}^{-1}]`.


Counts
------

If your input X-ray data was in counts, the following outputs will be included in the
post-processed file.

``XRAY_EXPOSURE`` : double array(Nxray)
    The exposure time of each bandpass specified in the input :math:`[\rm{s}]`.

``NET_COUNTS`` : double array(Nxray)
    The net counts in each bandpass specified in the input :math:`[\rm{counts}]`.

``NET_COUNTS_UNC`` : double array(Nxray)
    The uncertainty on the net counts in each bandpass as calculated from the net counts
    (if ``XRAY_UNC = 'SQRT'`` or ``XRAY_UNC = 'GEHRELS'``) or as specified in the input 
    (if ``XRAY_UNC = 'USER'``) :math:`[\rm{counts}]`.

``XRAY_COUNTS_MOD`` : double(Nxray, Nmodels)
    The net counts produced by the model in each bandpass :math:`[\rm{counts}]`.
    In the case of the MPFIT algorithm, this is the best-fitting net counts in each bandpass.
    In the case of the MCMC algorithms, this is the posterior distribution on the net counts in each bandpass.

``LNU_XRAY_OBS`` : double(Nxray, Nmodels)
    The model-dependent "observed" luminosities in each bandpass as derived from the net counts
    in terms of :math:`L_\nu` :math:`[L_\odot\ {\rm Hz}^{-1}]`.

``LNU_XRAY_UNC`` : double(Nxray, Nmodels)
    The uncertainty on the model-dependent "observed" luminosities in each bandpass as derived
    from the net counts uncertainty :math:`[L_\odot\ {\rm Hz}^{-1}]`.


Flux
----

If your input X-ray data was in terms of flux, the following outputs will be included in the
post-processed file.

``LNU_XRAY_OBS`` : double(Nxray)
    The observed X-ray luminosities converted from input X-ray flux data for each set of bandpasses:

    .. math::

    	L_\nu = \frac{4 \pi C (D_L)^2 F_{\nu}}{E_{\rm upper} - E_{\rm lower}},

    where :math:`L_{\nu}` is the observed X-ray luminosities :math:`[L_\odot\ {\rm Hz}^{-1}]`,
    :math:`C` is the unit conversion constant (:math:`C = 1.0247 \times 10^{-2}`),
    :math:`D_L` is the luminosity distance :math:`[{\rm Mpc}]`, :math:`F_{\nu}` is the
    input flux data :math:`[{\rm erg\ cm^{-2}\ s^{-1}}]`, :math:`E_{\rm upper}` is the upper 
    energy bound of the bandpass :math:`[{\rm keV}]`, and :math:`E_{\rm lower}` is the lower
    energy bound of the bandpass :math:`[{\rm keV}]`.


``LNU_XRAY_UNC`` : double(Nxray)
    The uncertainties on the observed X-ray luminosities converted from input 
    X-ray flux uncertainties :math:`[L_\odot\ {\rm Hz}^{-1}]`.



X-ray AGN
---------

If fitting an X-ray AGN model, the following outputs will be included in the
post-processed file.

``LNU_XRAYMOD_AGN_HIRES`` : double(200, Nhighres_models)
    The AGN contribution to the total high resolution X-ray luminosities,
    after absorption :math:`[L_\odot\ {\rm Hz}^{-1}]`.

``LNU_XRAYMOD_AGN_UNABS_HIRES`` : double(200, Nhighres_models)
    The intrinsic AGN contribution to the total high resolution 
    X-ray luminosities, without absorption :math:`[L_\odot\ {\rm Hz}^{-1}]`.

