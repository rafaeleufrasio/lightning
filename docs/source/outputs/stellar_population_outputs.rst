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

``MSTAR`` : double array(Nmodels)
    The total surviving stellar mass :math:`M_{\star}` in :math:`\rm M_\odot` produced by the resulting
    SFH.

``MSTAR_UNC`` : double array(Nmodels)
    The estimated :math:`1\sigma` uncertainty for the total stellar mass in :math:`\rm M_\odot`.

    .. note::
   
        Only appears in the output if using the MPFIT algorithm.

``STEPS_MSTAR`` : double array(Nsteps, Nmodels)
    The surviving stellar mass components :math:`M_{\star,i}` in :math:`\rm M_\odot` produced by each
    age bin :math:`i` of the star formation history.
    The total stellar mass :math:`M_{\star}` is then

    .. math::

        M_{\star} = \sum_i M_{\star,i}.

``STEPS_MSTAR_UNC`` : double array(Nsteps, Nmodels)
    The estimated :math:`1\sigma` uncertainty for the surviving stellar mass components in :math:`\rm M_\odot`.

    .. note::
   
        Only appears in the output if using the MPFIT algorithm.

