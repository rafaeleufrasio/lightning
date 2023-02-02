AGN Model Outputs
=================

``LNU_AGNMOD`` : double array(Nfilters[, Nsamples])
    The AGN contribution to the mean :math:`L_\nu` in each filter,
    after reddening by the ISM dust attenuation model :math:`[L_\odot\ {\rm Hz}^{-1}]`.

``LNU_AGNMOD_UNRED`` : double array(Nfilters[, Nsamples])
    The AGN contribution to the mean :math:`L_\nu` in each filter,
    with no reddening by the ISM dust attenuation model :math:`[L_\odot\ {\rm Hz}^{-1}]`.

``LNU_AGNMOD_HIRES`` : double array(1000, Nhighres_models)
    The AGN contribution to the total high resolution UV-to-FIR luminosities,
    after reddening by the ISM dust attenuation model :math:`[L_\odot\ {\rm Hz}^{-1}]`.

``LNU_AGNMOD_UNRED_HIRES`` : double array(1000, Nhighres_models)
    The AGN contribution to the total high resolution UV-to-FIR luminosities,
    with no reddening by the ISM dust attenuation model :math:`[L_\odot\ {\rm Hz}^{-1}]`.
