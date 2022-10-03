Dust Model Outputs
==================

``LNU_DUSTMOD`` : double array(Nfilters, Nmodels)
    The dust emission model contribution to the mean :math:`L_\nu` in each filter
    :math:`[L_\odot\ {\rm Hz}^{-1}]`.

``LNU_DUSTMOD_HIRES`` : double array(1000, Nhighres_models)
    The dust emission model contribution to the total high resolution UV-to-FIR luminosities
    :math:`[L_\odot\ {\rm Hz}^{-1}]`.

``LTIR`` : double array(Nmodels)
    The total infrared (TIR) luminosity (i.e., total bolometric luminosity) produced by
    the dust emission model :math:`[L_\odot]`.

    .. note::

        Since the dust model has wavelengths from 1-10000 :math:`\mu m`, :math:`L_{\rm TIR}`
        is integrated in this wavelength range. We note that this differs from the typically
        reported :math:`L_{\rm TIR}`, which is integrated from 8-1000 :math:`\mu m`. However,
        the dust emission model has minimal emission above 1000 :math:`\mu m` and below 5 :math:`\mu m`.
        Therefore, the differences in :math:`L_{\rm TIR}` values using these different wavelength ranges
        should be minimal.

        If :math:`L_{\rm TIR}` is a free parameter in the dust emission model (i.e., ``ENERGY_BALANCE`` is not 
        set in the :ref:`configuration <configure-setting-label>`), then this output will be
        considered part of the :ref:`parameter-output-label` and will not be duplicated here.
