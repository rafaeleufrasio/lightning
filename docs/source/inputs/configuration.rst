.. _configure-setting-label:

Configuration Settings
======================

The configuration has eight sections: :ref:`core-config-label`, :ref:`stellar-config-label`,
:ref:`dust-atten-config-label`, :ref:`dust-config-label`, :ref:`xray-config-label`,
:ref:`agn-config-label`, :ref:`fit-algor-config-label`, and :ref:`postpro-config-label`.
Of these sections, :ref:`core-config-label` gives the basic settings and core assumptions
within Lightning (e.g., cosmology). The next five, :ref:`stellar-config-label`,
:ref:`dust-atten-config-label`, :ref:`dust-config-label`, :ref:`xray-config-label`,
and :ref:`agn-config-label`, make up the model selection sections. Lastly,
:ref:`fit-algor-config-label` gives the choice of fitting algorithm and its corresponding
hyper-parameters, and :ref:`postpro-config-label` gives the choice of how to post-process
the resulting data.

.. note::

    If you are unsure on what model to select for fitting your SED(s), see our guide on
    :ref:`model-select-label`. The same goes for the fitting algorithm. See our guide on
    :ref:`algorithm-select-label` to help you decide on a fitting algorithm or if you
    are having issues with your fits.


.. _core-config-label:

Core
----

``PRINT_PROGRESS`` : flag
    A flag that indicates if the progress of Lightning should be printed to the
    terminal. This progress includes the current elapsed time, completed processes,
    and expected time remaining.

``MAX_CPUS`` : int scalar
    The maximum number of CPUs on the machine to utilize. If this value exceeds the
    actual number of CPUs on the machine, then all CPUs will be used. This setting
    allows for Lightning to fit SEDs in parallel, where one SED is fit per CPU
    (i.e., a number of SEDs equal to ``MAX_CPUS`` will be fit simultaneously in
    batches until all SEDs have been fit).

``ENERGY_BALANCE`` : flag
    A flag indicating if energy balance should be assumed in the SED fits. Energy balance
    is the assumption that the total integrated IR luminosity of the dust emission is equal
    to the total absorbed stellar (and, if set, AGN) emission.

    .. note::

        This is a key assumption in most SED fitting codes as it attempts to preserve conservation
        of energy. See our guide on :ref:`model-select-label` if you are unsure if you
        want energy balance in your model.


``MODEL_UNC`` : int, float, or double scalar
    The fractional model uncertainty to use in all filters when computing :math:`\chi^2`
    during the SED fitting. This form of uncertainty accounts for systematic effects
    in the models and is computed as

    .. math::

    	\sigma_{{\rm mod},\ i}^2 = \big({\tt MODEL\_UNC} * L_{\nu,\ i}^{\rm mod} \big)^2,

    where :math:`\sigma_{{\rm mod},\ i}` is the model uncertainty of filter :math:`i`, and
    :math:`L_{\nu,\ i}^{\rm mod}` is the model luminosity of filter :math:`i`. The total
    uncertainty used in the :math:`\chi^2` calculation is then given as

    .. math::

    	\sigma_{{\rm total},\ i}^2 = \sigma_{{\rm obs},\ i}^2 + \sigma_{{\rm mod},\ i}^2,

    where :math:`\sigma_{{\rm total},\ i}` is the total uncertainty of filter :math:`i`, and
    :math:`\sigma_{{\rm obs},\ i}` is the observed uncertainty of filter :math:`i` as given
    in the input.

    .. note::

        It is common in the literature to assume a fractional model uncertainty of 5-10%, regardless
        of SED fitting code. Therefore, we recommend using a fractional model uncertainty of 5%
        when fitting any SED for the first time.


Cosmology
^^^^^^^^^
The next five settings are the cosmology parameters to use in the SED fitting.
These parameters determine the assumed cosmology, which set the age of the universe
and the distance to objects if their distance was specified by redshift.

``H0`` : int, float, or double scalar
    The Hubble constant, :math:`H_0` :math:`[{\rm km\ s^{-1}\ Mpc^{-1}}]`.

``OMEGA_M`` : int, float, or double scalar
    The matter density normalized to the closure density, :math:`\Omega_m`.

``LAMBDA0`` : int, float, or double scalar
    The cosmological constant normalized to the closure density, :math:`\Lambda_0`.

``Q0`` : int, float, or double scalar
    The deceleration parameter, :math:`q_0`.

``K`` : int, float, or double scalar
    The curvature constant normalized to the closure density, :math:`k`.




.. _stellar-config-label:

Stellar Emission
----------------

``SSP`` : string scalar
    The simple stellar population (SSP)
    models to use for the stellar population. The only SSP models currently available in Lightning are the
    `PEGASE <http://www2.iap.fr/pegase/>`_ models. These models are selected by setting ``SSP`` to ``'PEGASE'``.
    To fit the SEDs without any stellar emission, set ``SSP`` to ``'NONE'``.

.. note::

    If no stellar emission model is chosen, all stellar emission model settings below
    can be skipped.


``IMF`` : string scalar
    The initial mass function (IMF)
    to use in the SSP models. The only IMF currently available in Lightning is that from
    `Kroupa (2001) <https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K/abstract>`_.
    This IMF is selected by setting ``IMF`` to ``'KROUPA01'``.

``ZMETAL`` : float or double scalar
    The metallicity to use in the SSP models
    in terms of Z, normalized to the solar metallicity.
    The current available metallicities in Lightning are 0.001, 0.004, 0.008, 0.02, 0.05, and 0.1
    in terms of :math:`Z`.

    .. note::

        Lightning currently assumes the chosen metallicity is constant for at all ages, and
        does not allow for metallicity evolution. To minimize any systematic effects caused
        by ignoring metallicity evolution, we recommend selecting a metallicity closest to
        current average metallicity of your input SEDs.


``EMISSION_LINES`` : flag
    A flag indicating if nebular emission lines should be included in the SSP models.

``NEBULAR_EXTINCTION`` : flag
    A flag indicating if nebular extinction should be included in the SSP models.

``SFH`` : string scalar
    The type of star formation history (SFH) to assume when fitting the SEDs. The only
    SFH type currently available in Lightning is the binned or "non-parametric" SFH.
    This SFH assumes a piece-wise constant SFH, where the SFR is a constant value within
    a set of age bins. This SFH type is selected by setting ``SFH`` to ``'NON-PARAMETRIC'``.


``STEPS_BOUNDS`` : int, float, or double array(Nsteps+1)
    The age bin (or step) boundaries to use in the "non-parametric" SFH in units of
    :math:`{\rm yr}`. Values must be in ascending order.

    .. note::

        If an age bin contains ages older than the universe at an input
        SED's redshift, the age bin upper bound will be automatically
        adjusted to the age of the universe at that redshift. If an entire
        age bin is older than universe at that redshift, then the entire
        age bin will be omitted and the next younger bin will be adjusted
        accordingly.

``DTIME_SF`` : int, float, or double scalar
    The time step used for interpolating the SSP models into the age
    bins in units of :math:`{\rm yr}`.

    .. warning::

        We do not recommend changing this value from its default. The only case
        in which it should be changed is if you specified age bins with differences less than
        the default value. However, in that case, your age bins are likely too small.


``PSI`` : structure
    The free parameter :math:`\psi_i`, the SFR for of the SFH age bin :math:`i` in :math:`M_\odot\ {\rm yr}^{-1}`.
    This structure contains the priors to assume for each :math:`\psi_i`.
    Values of :math:`\psi_i` are limited to being non-negative numbers.

    .. note::

        Check out the :ref:`priors-label` for details on what a prior structure contains
        and various examples.


.. _dust-atten-config-label:

Dust Attenuation
----------------

``ATTEN_CURVE`` : string scalar
    The assumed attenuation curve to apply to the stellar and/or AGN models. There are three attenuation
    curve options currently available in Lightning. They are the `Calzetti et al. (2000)
    <https://ui.adsabs.harvard.edu/abs/2000ApJ...533..682C/abstract>`_ attenuation curve,
    modified Calzetti et al. (2000) attenuation curve, and `Doore et al. (2021)
    <https://ui.adsabs.harvard.edu/abs/2021ApJ...923...26D/abstract>`_ attenuation curve. The
    modified Calzetti curve can include a variable slope as described in
    `Noll et al. (2009) <https://ui.adsabs.harvard.edu/abs/2009A%26A...507.1793N/abstract>`_,
    an optional 2175 Angstrom bump feature specified in `Kriek & Conroy (2013)
    <https://ui.adsabs.harvard.edu/abs/2013ApJ...775L..16K/abstract>`_, and birth cloud attenuation
    as described in `Eufrasio et al. (2017) <https://ui.adsabs.harvard.edu/abs/2017ApJ...851...10E/abstract>`_.
    The Doore et al (2021) attenuation curve is based on the `Tuffs et al. (2004)
    <https://ui.adsabs.harvard.edu/abs/2004A%26A...419..821T/abstract>`_ attenuation curves as updated
    by `Popescu et al. (2011) <https://ui.adsabs.harvard.edu/abs/2011A%26A...527A.109P/abstract>`_.
    These attenuation curves are selected by setting ``ATTEN_CURVE`` to ``'CALZETTI00'``, ``'CALZETTI_MOD'``,
    or ``'DOORE21'``, respectively.

    .. note::

        Attenuation of AGN can only use the ``'CALZETTI00'`` or ``'CALZETTI_MOD'``
        attenuation curves. Compatibility of the AGN models with the ``'DOORE21'``
        curve is currently not supported.


Calzetti+00
^^^^^^^^^^^

``TAUV`` : structure
    The free parameter :math:`\tau_V`, the V-band optical depth used for normalization
    in the Calzetti et al. (2000) attenuation curve.
    This structure contains the prior to assume for :math:`\tau_V`.
    Values of :math:`\tau_V` are limited to being non-negative numbers.


Modified Calzetti+00
^^^^^^^^^^^^^^^^^^^^

``TAUV_DIFF`` : structure
    The free parameter :math:`\tau_V^{\rm diff}`, the V-band optical depth of diffuse dust
    used for normalization in the Calzetti et al. (2000) attenuation curve.
    This structure contains the prior to assume for :math:`\tau_V^{\rm diff}`.
    Values of :math:`\tau_V^{\rm diff}` are limited to being non-negative numbers.

``DELTA`` : structure
    The free parameter :math:`\delta`, the power law value used to create a variable attenuation
    curve slope as described in Noll et al. (2009).
    This structure contains the prior to assume for :math:`\delta`.
    Values of :math:`\delta` can be any real numbers. A value of ``0`` indicates the same
    slope as the original Calzetti et al. (2000) attenuation curve.

``TAUV_BC`` : structure
    The free parameter :math:`\tau_V^{\rm BC}`, the V-band optical depth of the birth cloud component
    as described in Eufrasio et al. (2017).
    This structure contains the prior to assume for :math:`\tau_V^{\rm BC}`.
    Values of :math:`\tau_V^{\rm BC}` are limited to being non-negative numbers. A value of ``0``
    indicates no birth cloud attenuation.

``UV_BUMP`` : flag
    A flag indicating if a 2175 Angstrom UV bump feature as specified in Kriek & Conroy (2013)
    should be added to the attenuation curve.


Doore+21
^^^^^^^^
``TAUB_F`` : structure
    The free parameter :math:`\tau_B^{f}`, the face-on optical depth in the B-band.
    This structure contains the prior to assume for :math:`\tau_B^{f}`.
    Values of :math:`\tau_B^{f}` are limited to being between ``0`` and ``8``.

``F_CLUMP`` : structure
    The free parameter :math:`F`, the birth cloud clumpiness factor.
    This structure contains the prior to assume for :math:`F`.
    Values of :math:`F` are limited to being between ``0`` and ``0.61``.

``COSI`` : structure
    The free parameter :math:`\cos i`, the inclination of the galactic disk in terms of :math:`\cos i`.
    This structure contains the prior to assume for :math:`\cos i`.
    Values of :math:`\cos i` are limited to being between ``0`` and ``1``.

``B_TO_D`` : structure
    The free parameter :math:`B/D`, the bulge-to-disk ratio.
    This structure contains the prior to assume for :math:`B/D`.
    Values of :math:`B/D` are limited to being non-negative numbers.

``ROLD0_AGES`` : int, float, or double array(Nsteps)
    The binary parameter :math:`r^{0,\ {\rm old}}`, that designates each SFH age bin
    as part of the young or old population. A value of ``0`` for the corresponding age
    bin considers it to be part of the young population, and a value of ``1`` considers
    it to be part of the old populations (see section 4.3 of `Doore et al. 2021
    <https://ui.adsabs.harvard.edu/abs/2021ApJ...923...26D/abstract>`_ for more details).
    The number of elements must be one less than the number of elements in ``STEPS_BOUNDS``.

    .. note::

        We recommend setting age bins that contain ages :math:`< 500\ {\rm Myr}` to be part
        the young population as they can contain significant UV emission. If you choose
        to set age bins with ages :math:`< 500\ {\rm Myr}` to the old population, the SFR may
        be underestimated due to under-attenuation of the UV-emitting population.


.. _dust-config-label:

Dust Emission
-------------

``DUST_MODEL`` : string scalar
    The dust emission model to use. The only dust emission model currently available in Lightning is the
    `Draine & Li (2007) <https://ui.adsabs.harvard.edu/abs/2007ApJ...657..810D/abstract>`_ (DL07) model.
    This model is selected by setting ``DUST_MODEL`` to ``'DL07'``.
    To fit the SEDs without any dust emission, set ``DUST_MODEL`` to ``'NONE'``.

.. note::

    If no dust emission model is chosen, all dust emission model settings below
    can be skipped.


DL07
^^^^

``UMIN`` : structure
    The free parameter :math:`U_{\rm min}`, the minimum radiation field intensity
    of the diffuse ISM radiation field from the heated dust.
    This structure contains the prior to assume for :math:`U_{\rm min}`.
    Values of :math:`U_{\rm min}` are limited to being between ``0.1`` and ``25``.

``UMAX`` : structure
    The free parameter :math:`U_{\rm max}`, the maximum radiation field intensity
    of the power-law distribution of heating starlight intensities.
    This structure contains the prior to assume for :math:`U_{\rm max}`.
    Values of :math:`U_{\rm max}`` are limited to being between ``1e3`` and ``3e5``.

    .. note::

        The parameter range of :math:`U_{\rm max}` is slightly less than the quoted full
        range of the DL07 models (:math:`10^6`). This slightly limited range originates
        from the format of the `publicly available data
        <https://www.astro.princeton.edu/~draine/dust/irem.html>`_. The publicly available
        :math:`\delta`-functions of :math:`U`, from which :math:`U_{\rm max}`` can be
        calculated for any given :math:`\alpha`, have a maximum value of :math:`3 \times 10^5`.
        However, rather than extrapolating these :math:`\delta`-functions to
        :math:`U = 10^6`, we limit :math:`U_{\rm max}`` to the largest available value.

``ALPHA`` : structure
    The free parameter :math:`\alpha`, the exponent of the power-law distribution of
    heating starlight intensities between :math:`U_{\rm min}` and :math:`U_{\rm max}`.
    This structure contains the prior to assume for :math:`\alpha`.
    Values of :math:`\alpha` are limited to being between ``-10`` and ``4``.

``GAMMA`` : structure
    The free parameter :math:`\gamma`, the fraction of the dust mass exposed to
    the power-law distribution of radiation field intensities.
    This structure contains the prior to assume for :math:`\gamma`.
    Values of :math:`\gamma` are limited to being between ``0`` and ``1``.

``QPAH`` : structure
    The free parameter :math:`q_{\rm PAH}`, the fraction of the total grain mass
    corresponding to PAHs containing less than 1000 carbon atoms (PAH index).
    This structure contains the prior to assume for :math:`q_{\rm PAH}`.
    Values of :math:`q_{\rm PAH}` are limited to being between ``4.7e-3`` and ``4.58e-2``.

``LTIR`` : structure
    The free parameter :math:`L_{\rm TIR}`, the total integrated IR luminosity in :math:`L_\odot`.
    This structure contains the prior to assume for :math:`L_{\rm TIR}`.
    Values of :math:`L_{\rm TIR}` are limited to being non-negative numbers.

    .. note::

        ``LTIR`` is only a free parameter if ``ENERGY_BALANCE`` not is set. If ``ENERGY_BALANCE``
        is set then ``LTIR`` is determined instead by the absorbed the stellar (and, if set, AGN)
        emission.


.. _xray-config-label:

X-ray Emission
--------------

``XRAY_EMISSION`` : flag
    A flag indicating if an X-ray emission model will be used. This always includes
    stellar X-ray emission, but can optionally include AGN X-ray emission
    (:ref:`see below <xray-agn-config-label>`). The stellar X-ray emission is normalized
    according to the :math:`L_X/M` parametrizations with stellar age from `Gilbertson et
    al. (2022) <https://ui.adsabs.harvard.edu/abs/2022ApJ...926...28G/abstract>`_.

.. note::

    If no X-ray emission model is used, all X-ray emission model settings below
    can be skipped.


``XRAY_UNIT`` : string scalar
    The form (or unit type) of X-ray data within the input catalog.
    Currently, there are two types of X-ray data that can be input into Lightning.
    These are instrumental counts or fluxes (in :math:`{\rm erg\ cm^{-2}\ s^{-1}}`), which are
    selected by setting ``XRAY_UNIT`` to ``'COUNTS'`` or ``'FLUX'``, respectively.
    See the discussion on :ref:`input-formats-label` for more details on how to
    format the different X-ray data types.

    .. note::

        If set to ``'FLUX'``, the ``XRAY_UNC`` setting below is ignored.
        Uncertainties on the X-ray flux must always be provided in the input catalog.


``XRAY_UNC`` : string scalar
    The type of uncertainties to assume for the X-ray counts.
    In Lightning, the contribution to :math:`\chi^2` from the X-ray model is
    calculated as

    .. math::

        \chi^2_X = \sum_i \frac{(n^{\rm obs}_i - n^{\rm mod}_i)^2}{\sigma_{n,\ i}^2},

    where :math:`n^{\rm obs}_i` is the number of net (background-subtracted)
    counts in energy bin :math:`i`, :math:`n^{\rm mod}_i` is the number of model counts,
    and :math:`\sigma_{n,\ i}` is the uncertainty on the observed counts. There are
    three types of X-ray count uncertainties currently available in Lightning. They are
    the square root of the counts, the upper uncertainty from the `Gehrels (1986)
    <https://ui.adsabs.harvard.edu/abs/1986ApJ...303..336G/abstract>`_ approximation,
    and user input uncertainties. For the square root of the counts,
    :math:`\sigma_{n,\ i}` is assumed to be :math:`\sqrt{n^{\rm obs}_i}`.
    This is most appropriate for cases where the number of counts is large enough
    that the errors are approximately Gaussian. For the Gehrels (1986) approximation,
    :math:`\sigma_{n,\ i}` is assumed to be

    .. math::

        1 + \sqrt{0.75 + n^{\rm obs}_i}.

    This is most appropriate for data in the low-count regime.
    Finally, for the user input uncertainties, Lightning searches each X-ray spectral file
    for a column labeled ``NET_COUNTS_UNC`` and adopts this as the
    uncertainty on the net counts.
    These uncertainties types are selected by setting ``XRAY_UNC`` to ``'SQRT'``, ``'GEHRELS'``,
    or ``'USER'``, respectively.


``XRAY_ABS_MODEL`` : string scalar
    The X-ray absorption model to apply to the X-ray emission. There are
    three X-ray absorption models currently available in Lightning. They are
    the `"tbabs" absorption model <https://ui.adsabs.harvard.edu/abs/2000ApJ...542..914W/abstract>`_
    with the default `Wilms et al. (2000) <https://ui.adsabs.harvard.edu/abs/2000ApJ...542..914W/abstract>`_
    abundances, the "tbabs" model with `Anders & Grevesse (1986)
    <https://ui.adsabs.harvard.edu/abs/1989GeCoA..53..197A/abstract>`_ abundances, and the
    *Sherpa* "atten" model from `Rumph et al. (1994)
    <https://ui.adsabs.harvard.edu/abs/1994AJ....107.2108R/abstract>`_.
    These X-ray absorption models are selected by setting ``XRAY_ABS_MODEL`` to ``'TBABS-WILM'``,
    ``'TBABS-ANGR'``, or ``'ATTEN'``, respectively.

``NH`` : structure
    The free parameter :math:`N_H`, the intrinsic HI column density along
    the line of sight in :math:`10^{20}\ {\rm cm}^{-2}`.
    This structure contains the prior to assume for :math:`N_H`.
    Values of :math:`N_H` are limited to being between ``1e-4`` and ``1e5``.

    .. note::

    	While the value of :math:`N_H` is allowed to be larger than :math:`10^{24}\ {\rm cm^{-2}}`,
        we caution that our emission models are not suitable for the Compton-thick case.


.. _xray-agn-config-label:

``XRAY_AGN_MODEL`` : string scalar
    The AGN X-ray emission model to use. There are two AGN X-ray emission models currently available
    in Lightning, the `qsosed <https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node132.html>`_
    models from `Kubota & Done (2018) <https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.1247K/abstract>`_ and a power law model
    with and exponential cut off. The power law model has a photon index of :math:`\Gamma = 1.8` and an exponential
    cut off at 300 :math:`{\rm keV}`. This power law model is tied to the 2500 Angstrom emission using the
    relationship from `Lusso & Risaliti (2017) <https://ui.adsabs.harvard.edu/abs/2017A%26A...602A..79L/abstract>`_.
    These models are selected by setting ``XRAY_AGN_MODEL`` to ``'QSOSED'`` and ``'PLAW'``, respectively.
    To fit the SEDs without any AGN X-ray emission models, set ``XRAY_AGN_MODEL`` to ``'NONE'``.

.. note::

    If the ``'QSOSED'`` AGN X-ray emission model is not chosen, its corresponding settings
    below can be skipped.


QSOSED
^^^^^^

``AGN_MASS`` : structure
    The free parameter :math:`M_{\rm AGN}`, the supermassive black hole mass in
    :math:`M_\odot`.
    This structure contains the prior to assume for :math:`M_{\rm AGN}`.
    Values of :math:`M_{\rm AGN}` are limited to being between ``1e5`` and ``1e10``.

``AGN_LOGMDOT`` : structure
    The free parameter :math:`\log(\dot m)`, the :math:`\log_{10}` of :math:`\dot m`,
    the supermassive black hole accretion rate normalized by the Eddington rate.
    This structure contains the prior to assume for :math:`\log(\dot m)`.
    Values of :math:`\log(\dot m)` are limited to being between ``-1.5`` and ``0.3``.


.. _agn-config-label:

AGN Emission
------------

``AGN_MODEL`` : string scalar
    The UV-to-IR AGN emission model to use. The only AGN emission model currently available in Lightning is the
    `SKIRTOR <http://sites.google.com/site/skirtorus/home>`_ model from `Stalevski et al. (2016)
    <https://ui.adsabs.harvard.edu/abs/2016MNRAS.458.2288S/abstract>`_.
    This model is selected by setting ``AGN_MODEL`` to ``'SKIRTOR'``.
    To fit the SEDs without any AGN emission, set ``AGN_MODEL`` to ``'NONE'``.

.. note::

    If no AGN emission model is chosen, all AGN emission model settings below
    can be skipped.


SKIRTOR
^^^^^^^

``LOG_L_AGN`` : structure
    The free parameter :math:`\log(L_{\rm AGN})`, the total integrated luminosity of
    AGN model in :math:`\log_{10}(L_\odot)`, which is used for normalization.
    This structure contains the prior to assume for :math:`\log(L_{\rm AGN})`.
    Values of :math:`\log(L_{\rm AGN})` are limited to being between ``0`` and ``20``.

    .. note::

        :math:`\log(L_{\rm AGN})` will not be a free parameter if fitting using
        a ``'QSOSED'`` X-ray AGN model. Instead the normalization of
        UV-to-IR AGN model is tied to the rest-frame 2500 Angstrom monochromatic
        luminosity of the qsosed model.

``TAU97`` : structure
    The free parameter :math:`\tau_{9.7}`, the edge-on optical depth of AGN dust
    torus at 9.7 :math:`\mu \rm m`.
    This structure contains the prior to assume for :math:`\tau_{9.7}`.
    Values of :math:`\tau_{9.7}` are limited to being between ``3`` and ``11``.

``AGN_COSI`` : structure
    The free parameter :math:`\cos i_{\rm AGN}`, the inclination of the AGN
    disk in terms of :math:`\cos i`.
    This structure contains the prior to assume for :math:`\cos i_{\rm AGN}`.
    Values of :math:`\cos i_{\rm AGN}` are limited to being between ``0`` and ``1``.


.. _fit-algor-config-label:

Fitting Algorithm
-----------------

``METHOD`` : string scalar
    The fitting algorithm used to fit the SED(s). Lightning currently has three fitting algorithms
    that can be used: an adaptive MCMC, an affine-invariant MCMC, and a Levenberg–Marquardt algorithm.
    The adaptive MCMC algorithm is Algorithm 4 from `Andrieu & Thoms (2008)
    <https://link.springer.com/article/10.1007/s11222-008-9110-y>`_, the affine-invariant MCMC
    algorithm is the algorithm from `Goodman & Weare (2010)
    <https://ui.adsabs.harvard.edu/abs/2010CAMCS...5...65G/abstract>`_, and the Levenberg–Marquardt
    algorithm is `Craig Markwardt’s MPFIT <http://purl.com/net/mpfit>`_
    implementation.
    These fitting algorithms are selected by setting ``METHOD`` to ``'MCMC-ADAPTIVE'``, ``'MCMC-AFFINE'``,
    or ``'MPFIT'``, respectively.


.. note::

    See our guide on :ref:`algorithm-select-label` for more details on each algorithm and
    their corresponding hyper-parameters below. Additionally, the guide can help you decide
    on the best algorithm to fit your research needs.



MCMC
^^^^

``NTRIALS`` : int, float, or double scalar
    The number of MCMC trials to run for each parallel walker/chain.

``NPARALLEL`` : int, float, or double scalar
    The number of parallel walkers/chains.

    .. note::

        If using the affine-invariant algorithm, ``NPARALLEL`` must be greater than
        the number of free parameters plus one and ideally at least twice the number
        of free parameters for optimal sampling.


``C_STEP`` : int, float, or double scalar
    When calculating the autocorrelation time (:math:`\tau`) of the MCMC chain, this value
    defines how many trials of the chain are used to calculate :math:`\tau`, where
    we integrate :math:`\tau` to the smallest index :math:`M` such that :math:`M > C_{\rm step} \tau`.

``TOLERANCE`` : int, float, or double scalar
    When calculating the autocorrelation time (:math:`\tau`) of the MCMC chain, this value
    defines how many multiples of :math:`\tau` the length of the chain should be for us to believe
    the estimated value of :math:`\tau`.

.. note::

    We recommend using the default values for both ``C_STEP`` and ``TOLERANCE``. More details on these parameters
    can be found in the `emcee Autocorrelation Analysis documentation
    <https://emcee.readthedocs.io/en/stable/tutorials/autocorr/#autocorr>`_.


``BETA_EXPONENT`` : float or double scalar
    The factor controlling how fast the adaptiveness of the adaptive MCMC algorithm vanishes.
    Larger values stop the adaptiveness in fewer trials.

    .. note::

        This is a setting only for the adaptive MCMC algorithm.

``AFFINE_A`` : int, float, or double scalar
    The move scaling constant defining the maximum and
    minimum step size of the affine-invariant stretch move.

    .. note::

        This is a setting only for the affine-invariant MCMC algorithm.


MPFIT
^^^^^

``NSOLVERS`` : int, float, or double scalar
  The number of times to solve for the best fit SED using different
  starting locations in parameters space.

``FTOL`` : float or double scalar
  The relative error desired in the sum of squares. Termination
  of the MPFIT algorithm occurs when both the actual and predicted
  relative reductions in the sum of squares are at most ``FTOL``.

``GTOL`` : float or double scalar
  The orthogonality desired between the function vector and the
  columns of the Jacobian matrix. Termination of the MPFIT algorithm
  occurs when the cosine of the angle between function vector and any
  column of the Jacobian matrix is at most ``GTOL`` in absolute value.

``XTOL`` : float or double scalar
  The relative error desired in the approximate solution. Termination
  of the MPFIT algorithm occurs when the relative error between two
  consecutive iterates is at most ``XTOL``.

``MAXITER`` : int, float, or double scalar
  The maximum number of MPFIT iterations to perform.


.. _postpro-config-label:

Post-processing
---------------

``KEEP_INTERMEDIATE_OUTPUT`` : flag
    A flag indicating that the intermediate ``.sav`` files produced by the fitting algorithm
    should not be deleted.

    .. note::

        This is useful if needing to inspect the original fits before post-processing.
        Typically this will not be necessary, but if you are having trouble getting
        quality fits, inspecting the original fits can help determine the issue.


MCMC Post-processing
^^^^^^^^^^^^^^^^^^^^

The next four settings are the MCMC post-processing settings. These are only
used if fitting with an MCMC algorithm, and they determine the how the MCMC
chains are handled during post-processing for conversion to the posterior
distributions.


``BURN_IN`` : int, float, or double scalar
    The number of initial MCMC trials to discard as the burn-in phase. If set to ``0``,
    then the number will be chosen automatically from the autocorrelation time as

    .. math::

        {\tt BURN\_IN} = {\rm ceiling}(2\ {\rm max}(\tau)),

    where :math:`\tau` is the autocorrelation time and ``ceiling`` is the ceiling function
    that rounds values up to the nearest integer.

    .. note::
        We highly recommend specifying a value rather than using the automatic
        calculation when using the adaptive MCMC algorithm as the chains can vary widely
        in the number of autocorrelation times needed for burn-in.


``THIN_FACTOR`` : int, float, or double scalar
    The factor to thin the MCMC chain after removing the burn-in trials. Thinning of an
    MCMC chain is common practice, and it helps reduce the correlation between trials in
    the chain. To clarify what a value of ``THIN_FACTOR`` means, here are a few examples.
    A value of ``10`` will only keep every 10th trial in the chain, and a value of ``1``
    will keep every trial (i.e., no thinning). Finally, if set to ``0``, then the value
    will be chosen automatically from the autocorrelation time as

    .. math::

        {\tt THIN\_FACTOR} = {\rm ceiling}(0.5\tau),

    where :math:`\tau` is the autocorrelation time and ``ceiling`` is the ceiling function
    that rounds values up to the nearest integer.

    .. note::

        We recommend specifying a ``THIN_FACTOR`` of ``4`` and ``0`` (the automatic
        calculation) when using the adaptive and affine-invariant MCMC algorithms,
        respectively. The reason for the ``4`` value with the adaptive MCMC algorithm
        is that unique elements within the chains are minimally correlated. However,
        by design of the algorithm, a new unique element is only accepted into the
        chain every four or so trials. Therefore, by thinning by a factor of four,
        each element in the final chain will typically be unique.


``FINAL_CHAIN_LENGTH`` : int, float, or double scalar
    The number of MCMC trials to include for the final distributions as taken from
    the truncated, thinned, and if necessary, merged chains. In other words,
    ``FINAL_CHAIN_LENGTH`` specifies then number of samples to include in the
    posterior distributions. To get the posterior distributions, the raw chains output
    by the MCMC algorithm have their burn-in discarded and are thinned and merged if
    necessary. Then, a number of samples equal to ``FINAL_CHAIN_LENGTH`` will be taken
    from the end of the remaining chain to serve as the posterior distribution.

    .. note::

        We recommend specifying a nice round value for ``FINAL_CHAIN_LENGTH`` such
        as ``250``, ``500``, ``1000``, ``2000``, etc. Larger values will increase the fine detail of
        the posterior distribution at the cost of increased post-processed file size.


``HIGH_RES_MODEL_FRACTION`` : int, float, or double scalar
    The fraction of samples from ``FINAL_CHAIN_LENGTH``, sorted by quality of fit,
    from which to generate high resolution models. If set to ``0``, then only the
    best fit high resolution model will be generated. This setting dictates how
    many high resolution models per SED will be in the post-processed file, which
    are useful for plotting purposes. Having a value of ``HIGH_RES_MODEL_FRACTION``
    greater than zero would allow for having pointwise uncertainties on the best fit
    high resolution model, which can show, for example, areas of the model which are
    not well constrained by the data. It is important to stress that this setting
    gives the fraction of the **best-fitting** models in the posterior distribution
    for which high-resolution SEDs will be computed and saved. As an example,
    setting ``HIGH_RES_MODEL_FRACTION`` to ``0.68`` will return the high resolution
    models for the best 68% of fits in the posterior distribution.

    .. warning::

        Including more than the best fit high resolution model can cause the file size of
        the post-processed file to balloon dramatically. Be careful when
        increasing this value above ``0``. Doing so will increase the file size by **at least**
        ``FINAL_CHAIN_LENGTH`` * ``HIGH_RES_MODEL_FRACTION`` * ``8`` kB per SED per model component.


``AFFINE_STRANDED_DEVIATION`` : int, float, or double scalar
    The number of standard deviations a walker must be below the median
    acceptance fraction of the ensemble to be considered a stranded walker.
    (See the :ref:`affine-mcmc-label` description for more details on stranded walkers.)
