.. _outputs-label:

Outputs
=======

Below, we describe the contents of the post-processed files in detail. Note that
that the size of each output is only for one SED. If multiple SEDs were in the input
catalogue, all outputs will have an additional final dimension the length of the number of SEDs.

.. _array-size-label:
.. note::

    We define the array size variables here for convenience.

    - ``Nfilters`` : the number of filters included in the input.
    - ``Nmodels``: the number of models. If using the MPFIT algorithm, 
      this value will be ``1``. If using an MCMC algorithm, this value is equal to
      ``FINAL_CHAIN_LENGTH`` in the :ref:`configuration <configure-setting-label>`.
    - ``Nhighres_models``: the number of high resolution models. If using the MPFIT algorithm, 
      this value will be ``1``. If using an MCMC algorithm, this value is determined by
      ``HIGH_RES_MODEL_FRACTION`` in the :ref:`configuration <configure-setting-label>`.
    - ``Nparam``: the number of parameters used in the SED fit (including fixed parameters).


Basic Outputs
-------------

By default, the post-processed outputs include:

``SED_ID`` : string
    The unique identifier for each SED. If not specified in the input, this value
    will be automatically generated as integer values.

``REDSHIFT`` : double
    The redshift of each SED. If not specified in the input, this value will be ``0``.

``LUMIN_DIST`` : double
    The luminosity distance of each SED. If not specified in the input, this value will
    be determined from the redshift using the cosmology chosen during :ref:`configuration
    <configure-setting-label>` :math:`[{\rm Mpc}]`.

``FILTER_LABELS`` : string array(Nfilters)
    The list of filter labels specified in the input.

``WAVE_FILTERS`` : double array(Nfilters)
    The mean wavelength of each filter :math:`[{\mu \rm m}]`:

    .. math::

    	\bar\lambda = \frac{\int \lambda T(\lambda) d\lambda}{\int T(\lambda) d\lambda},

    where :math:`T(\lambda)` is the filter transmission function.

``LNU_OBS`` : double array(Nfilters)
    The observed luminosities converted from input flux data:

    .. math::

    	L_\nu = 4 \pi C (D_L)^2 F_{\nu},

    where :math:`L_{\nu}` is the observed luminosities :math:`[L_\odot\ {\rm Hz}^{-1}]`,
    :math:`C` is the unit conversion constant (:math:`C = 2.4778 \times 10^{-8}`),
    :math:`D_L` is the luminosity distance given by ``LUMIN_DIST`` :math:`[{\rm Mpc}]`,
    and :math:`F_{\nu}` is the input flux data :math:`[{\rm Jy}]`.

``LNU_UNC`` : double array(Nfilters)
    The uncertainties on the observed luminosities converted from input flux uncertainties :math:`[L_\odot\ {\rm Hz}^{-1}]`.

``LNU_MOD`` : double array(Nfilters, Nmodels)
    The mean :math:`L_\nu` produced by the model in each filter :math:`[L_\odot\ {\rm Hz}^{-1}]`:

    .. math::

    	\bar L_\nu = \frac{\int T(\lambda) L_\nu d\lambda}{\int T(\lambda) d\lambda}.

    In the case of the MPFIT algorithm, this is the best-fitting mean :math:`L_\nu` in each filter.
    In the case of the MCMC algorithms, this is the posterior distribution on the mean :math:`L_\nu` in each filter.

``MODEL_UNC`` : double
    The fractional uncertainties on the model chosen during :ref:`configuration
    <configure-setting-label>`.

``WAVE_HIRES`` : double array(1000)
    The wavelength grid for the high resolution UV-to-FIR model :math:`[\mu \rm m]`.

``LNU_MOD_HIRES`` : double array(1000, Nhighres_models)
    The total high resolution UV-to-FIR luminosities produced by the model :math:`[L_\odot\ {\rm Hz}^{-1}]`.

``LNPROB`` : double array(Nmodels)
    The natural log probability of each model. In the case of the MPFIT algorithm, this is the best-fitting
    log probability.  In the case of an MCMC algorithm, this is the sampled posterior log probability.

``LNPROB_BESTFIT`` : double
    The best-fitting log probability value.

    .. note::
   
        Only appears in the output if using an MCMC algorithm. 

``CHI2`` : double array(Nmodels)
    The :math:`\chi^2` of each model calculated as

    .. math::

        \chi^2 = \sum_i \frac{(L_{\nu,\ i}^{\rm obs} - L_{\nu,\ i}^{\rm mod})^2}{\sigma_{{\rm total},\ i}^2},

    where :math:`L_{\nu,\ i}^{\rm obs}` is ``LNU_OBS`` in filter :math:`i`, :math:`L_{\nu,\ i}^{\rm mod}` is
    ``LNU_MOD``  in filter :math:`i`, and
    :math:`\sigma_{{\rm total},\ i}` is the total uncertainty  in filter :math:`i` from the combined 
    observational and model uncertainty (see ``MODEL_UNC`` in the :ref:`configuration <configure-setting-label>`
    for details). In the case of the MPFIT algorithm, this is the best-fitting
    :math:`\chi^2`.  In the case of an MCMC algorithm, this is the sampled posterior :math:`\chi^2`.

``CHI2_BESTFIT`` : double
    The best-fitting :math:`\chi^2` value.

    .. note::
   
        Only appears in the output if using an MCMC algorithm. 


``PARAMETER_NAMES`` : string array(Nparam)
    The names of the parameters used in the SED fitting (including fixed parameters).

    .. _parameter-names-label:
    .. note::
   
        If multiple SEDs were input, it is possible that they may have a different number of parameters.
        This can happen, for example, with the parameters ``PSI_1``, ``PSI_2``, ``PSI_3``, etc. (the SFH 
        coefficients for each age bin) if a given age bin is older than the estimated age of the universe.
        In this case, the ``PSI_*`` parameters associated with these age bins will not be included in
        ``PARAMETER_NAMES``. If this happens, the last entries in ``PARAMETER_NAMES`` will be left blank.

``COVARIANCE`` : double array(Nparam, Nparam)
    The covariance matrix for the model parameters. The square root of the diagonal elements
    gives the estimated :math:`1\sigma` uncertainty for each parameter.

    .. note::
   
        Only appears in the output if using the MPFIT algorithm.


.. _parameter-output-label:

Parameter Outputs
-----------------

The post-processed files also include the model parameters for each component of the model as follows, where
``<PARAM-NAME>`` is a proxy for any of the parameter names given in ``PARAMETER_NAMES``:

.. note::

    While ``PARAMETER_NAMES`` has each parameter as individual entries, ``<PARAM-NAME>`` will compress
    like parameters into a single entry. This will occur, for example, with the ``PSI`` parameters.
    Rather than having ``PSI_1``, ``PSI_2``, ``PSI_3``, etc., they will all be compressed into a single
    ``PSI`` parameter that has an additional **leading** dimension with a length of the number of like parameters.
    If an entry in ``PARAMETER_NAMES`` is :ref:`blank <parameter-names-label>`, the corresponding ``<PARAM-NAME>``
    value will be ``NaN``.

``<PARAM-NAME>`` : double array(Nmodels)
    In the case of the MPFIT algorithm, this is the best-fitting value for each parameter. In the case of an MCMC algorithm,
    this is the sampled posterior for each parameter. 

    .. note::
   
        If the parameter was fixed, regardless of algorithm, this will be the fixed value. In the case of
        an MCMC algorithm, ``Nmodels`` will then be set to ``1`` for the fixed parameter to conserve memory.


``<PARAM-NAME>_PERCENTILES`` : double array(3, Nmodels)
    The 16th, 50th, and 84th percentiles of the sampled posterior distribution for each parameter.

    .. note::
   
        Only appears in the output if using an MCMC algorithm. 

``<PARAM-NAME>_BESTFIT`` : double
    The best-fitting value for each parameter.

    .. note::
   
        Only appears in the output if using an MCMC algorithm. 

``<PARAM-NAME>_UNC`` : double
    The estimated :math:`1\sigma` uncertainty for each parameter.

    .. note::
   
        Only appears in the output if using the MPFIT algorithm.


Other Outputs
-------------

The post-processed output will always contain the basic and parameter outputs. It will additionally contain
other outputs depending on the model and fitting algorithm chosen during configuration:

.. toctree::
    :maxdepth: 1

    stellar_population_outputs.rst
    dust_model_outputs.rst
    xray_model_outputs.rst
    agn_model_outputs.rst
    convergence_outputs.rst
