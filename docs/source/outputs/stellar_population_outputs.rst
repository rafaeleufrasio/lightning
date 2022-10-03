Stellar Population Outputs
==========================

.. note::

    We define additional array size variables here for convenience. 
    (See :ref:`here <array-size-label>` for the already defined array size variables.)

    - ``Nsteps`` : the number of chosen SFH age bins (or steps).


``LNU_STARMOD`` : double array(Nfilters, Nmodels)
    The stellar population contribution to the mean :math:`L_\nu` in each filter,
    after reddening by the dust attenuation model :math:`[L_\odot\ {\rm Hz}^{-1}]`.

``LNU_STARMOD_UNRED`` : double array(Nfilters, Nmodels)
    The intrinsic stellar population contribution to the mean :math:`L_\nu` in each
    filter, with no reddening by the dust attenuation model :math:`[L_\odot\ {\rm Hz}^{-1}]`.

``LNU_STARMOD_HIRES`` : double array(1000, Nhighres_models)
    The stellar population contribution to the total high resolution UV-to-FIR luminosities,
    after reddening by the dust attenuation model :math:`[L_\odot\ {\rm Hz}^{-1}]`.

``LNU_STARMOD_UNRED_HIRES`` : double array(1000, Nhighres_models)
    The intrinsic stellar population contribution to the total high resolution UV-to-FIR luminosities,
    with no reddening by the dust attenuation model :math:`[L_\odot\ {\rm Hz}^{-1}]`.

``STEPS_BOUNDS`` : double array(Nsteps + 1)
    The bounds on the age bins of the star formation history :math:`[{\rm yr}]`.

``STEPS_MSTAR_COEFF`` : double array(Nsteps)
    The surviving stellar mass :math:`m_{\star,i}` in :math:`\rm M_\odot` produced by SFR of :math:`1\ \rm M_\odot\ yr^{-1}` in
    each age bin of the star formation history. The total stellar mass :math:`M_{\star}` is then

    .. math::

        M_{\star} = \sum_i m_{\star,i} \psi_i

    where :math:`\psi_i` are the SFH coefficients.
