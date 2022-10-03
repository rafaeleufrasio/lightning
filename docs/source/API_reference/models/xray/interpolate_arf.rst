INTERPOLATE_ARF
===============

Name
----
INTERPOLATE_ARF

Purpose
-------
Gets the values of the specified Auxiliary Response Function (ARF)
at the specified wavelengths.

Calling Sequence
----------------
::

    specresp_interp = interpolate_arf(E_lo, E_hi, specresp, wave)

Inputs
------
``E_lo`` : int, float or double array(Nchannels)
    Lower energy bounds of each channel in the ARF :math:`[{\rm keV}]`
``E_hi`` : int, float or double array(Nchannels)
    Upper energy bounds of each channel in the ARF :math:`[{\rm keV}]`
``specresp`` : int, float or double array(Nchannels)
    The spectral response of the ARF at each channel :math:`[{\rm cm}^2]`
``wave`` : int, float or double array(Nwave)
    A grid of wavelengths at which to interpolate the ARF :math:`[\mu \rm m]`.

Output
------
``specresp_interp`` : double array(Nwave)
    A grid of ARF values interpolated from the input ARF file at each
    specified wavelength :math:`[{\rm cm}^2]`.

Modification History
--------------------
- 2021/08/26: Created (Erik B. Monson).
- 2022/06/22: Added documentation (Keith Doore)
- 2022/06/22: Replaced ``!cv`` with ``!lightning_cgs`` (Keith Doore)

