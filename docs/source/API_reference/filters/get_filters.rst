GET_FILTERS
===========

Name
----
GET_FILTERS

Purpose
-------
Generates a filter transmission array for the input filter labels.
The specified filters are read in from each tabulated filter profile
and interpolated onto the input frequency array.

Calling Sequence
----------------
::

    get_filters, filter_labels, nu, filters [, Lnu = ,   $
                 /error_check, /plot_filters, mean_nu=mean_nu, $
                 sigma_nu=sigma_nu, mean_Lnu=mean_Lnu, $
                 wave=wave, Llambda=Llambda, mean_wave=mean_wave, $
                 sigma_wave=sigma_wave, mean_Llambda=mean_Llambda]

Inputs
------
``filter_labels`` : string array(Nfilters)
    The names of the filters.
``nu`` : int, float, or double array(Nnu)
    The frequencies for interpolation :math:`[{\rm Hz}]`.

Optional Inputs
---------------
``Lnu`` : int, float, or double array(NLnu, Nnu)
    The luminosity densities at each element of ``nu`` to be used in
    the calculation of the optional outputs. The luminosity units can
    be user chosen and will dictate the optional output luminosity units.
    However, overall unit must be per Hertz (i.e., :math:`[L\ {\rm Hz}^{-1}]`).
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.
``plot_filters`` : flag
    If set, a plot will be generated showing the normalized filter
    transmissions as a function of wavelength.

Output
------
``filters`` : double array(Nfilters, Nnu)
    The tabulated filter transmission functions interpolated and normalized
    onto the input frequency array.

Optional Outputs
----------------
``mean_nu`` : double array(NLnu, Nfilters)
    The mean frequency of each filter as computed from each set of ``Lnu``
    :math:`[\rm Hz]`.
``sigma_nu`` : double array(NLnu, Nfilters)
    The standard deviation width of each filter as computed
    from each set of ``Lnu`` :math:`[\rm Hz]`.
``mean_Lnu`` : double array(NLnu, Nfilters)
    The mean luminosity density of each filter in the same units as
    ``Lnu`` (:math:`[L\ {\rm Hz}^{-1}]`) as computed from each set of ``Lnu``.
``wave`` : double array(Nnu)
    The wavelengths converted from ``nu`` :math:`[\mu \rm m]`.
``Llambda`` : double array(NLnu, Nnu)
    The luminosity density converted from ``Lnu`` in the same luminosity
    units as ``Lnu`` per micron (i.e. :math:`[L\ {\mu \rm m}^{-1}]`).
``mean_wave`` : double array(NLnu, Nfilters)
    The mean wavelength of each filter as computed from each set of ``Lnu``
    :math:`[\mu \rm m]`.
``sigma_wave`` : double array(NLnu, Nfilters)
    The standard deviation width of each filter as computed from each
    set of ``Lnu`` :math:`[\mu \rm m]`.
``mean_Llambda`` : double array(NLnu, Nfilters)
    The mean luminosity density of each filter in the same units as ``Llambda``
    (:math:`[L\ {\mu \rm m}^{-1}]`) as computed from each set of ``Llambda``.

Notes
-----
The optional outputs will only be computed if the optional input ``Lnu`` is
specified.

Modification History
--------------------
- 2016/05/01: Created (Rafael T. Eufrasio)
- 2019/04/28: Allowed for tabulated filter profile to exceed input ``nu`` range (Rafael T. Eufrasio)
- 2022/03/02: Added error handling (Keith Doore)
- 2022/03/02: Changed list of filters and directories from case statement to json file (Keith Doore)
- 2022/03/02: Minor bug fixes (Keith Doore)
- 2022/03/18: Updated documentation (Keith Doore)
- 2022/03/22: Updated to used plot function vs plot procedure (Keith Doore)
- 2022/04/12: Allowed for ``nu`` to have degenerate dimensions (Keith Doore)
- 2022/04/12: Allowed for ``filter_labels`` to be scalars (Keith Doore)
- 2022/04/12: Allowed integer inputs (Keith Doore)
- 2022/04/12: Made the output ``filters`` frequency normalized (Keith Doore)
- 2022/05/16: Removed ``filters_dir`` and changed to use ``!lightning_dir`` system variable (Keith Doore)
- 2022/05/16: Changed ``!cv`` constants call to ``!lightning_cgs`` call (Keith Doore)
- 2022/05/17: Added ``error_check`` keyword to do error handling (Keith Doore)
- 2022/06/08: Added ``strtrim`` to ``filter_labels`` to remove any extra blank padding (Keith Doore)
- 2022/06/29: Updated documentation and changed ``L_wave`` to ``Llambda`` (Keith Doore)

