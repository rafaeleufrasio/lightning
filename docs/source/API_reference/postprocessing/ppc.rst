PPC
===

Name
----
PPC

Purpose
-------
Computes a Posterior Predictive Check (PPC) on a MCMC Lightning output.
Uses the methods described in Rubin (1984) and Gelman et al. (1996).

Calling Sequence
----------------
::

    pvalue = ppc(Nrep, Lobs, Lunc, Lpredict, lnprob_chain [, $
                 counts_obs = , counts_unc = , counts_predict = , $
                 chisqr_obs=chisqr_obs, chisqr_rep=chisqr_rep])

Inputs
------
``Nrep`` : int, float or double scalar
    The number of times to replicate the predicted data.
``Lobs`` : int, float or double array(Nfilters)
    The observed luminosities [arbitrary units].
``Lunc`` : int, float or double array(Nfilters)
    The uncertainties on the observed luminosities in the same units as ``Lobs``.
``Lpredict`` : int, float, or double array(Nfilters, Nchain)
    The model luminosities predicted by the model given a set of parameters in
    same units as ``Lobs``.
``lnprob_chain`` : int, float, or double array(Nchain)
    The log posterior probability values of each SED predicted
    by the model given a set of parameters. Used to select the replicated data.

Optional Inputs
---------------
``counts_obs`` : int, float, or double array(Nxray)
    The observed X-ray counts. If specified, ``counts_unc`` and ``counts_predict``
    must also be given.
``counts_unc`` : int, float, or double array(Nxray)
    The uncertainties on the observed X-ray counts. If specified, ``counts_obs`` and ``counts_predict``
    must also be given.
``counts_predict``: int, float, or double array(Nxray, Nchain)
    The X-ray counts predicted by the model. If specified, ``counts_obs`` and ``counts_unc``
    must also be given.
``model_unc`` : int, float, or double scalar
    The fractional model uncertainty to use in all bands. (Default = ``0.d0``)

Output
------
``pvalue`` : double scalar
    The p-value associated with the PPC. Determined as the fraction
    of ``chisqr_rep`` that is greater than ``chisqr_obs``.

Optional Outputs
----------------
``chisqr_obs`` : double array(Nrep)
    The resulting :math:`\chi^2` values by comparing the predicted data with the
    observational data.
``chisqr_rep`` : double array(Nrep)
    The resulting :math:`\chi^2` values by comparing the predicted data with the
    replicated data.

Notes
-----
- If the X-ray model was fit using fluxes rather than counts, the X-ray luminosities can
  simply be appended to ``Lobs``, ``Lunc``, and ``Lpredict``.

References
----------
- `Rubin D. B., 1984, Ann. Stat., 12, 1151 <https://www.jstor.org/stable/2240995>`_
- `Gelman A., Meng X.-L., Stern H., 1996, Stat. Sin., 6, 733 <https://www.jstor.org/stable/24306036>`_

Modification History
--------------------
- 2022/03/15: Created (Keith Doore)
- 2022/04/18: Added documentation (Keith Doore)
- 2022/04/18: Added error handling (Keith Doore)
- 2022/07/21: Added ability to handle ``NaNs`` in ``lnprob_chain`` (Keith Doore)
- 2022/07/27: Fixed bug with missing keyword ``nrand`` in ``mrandomn`` (Keith Doore)
- 2022/07/27: Fixed bug where ``value_locate`` maybe selecting ``-1`` when we want ``0``,
  force ``-1`` to be ``0`` (Keith Doore)
- 2022/09/01: Added X-ray emission (Erik B. Monson)
- 2022/09/15: Changed from ``chi2_chain`` input to ``lnprob_chain`` (Keith Doore)
- 2022/09/30: Added ``model_unc`` input to include in uncertainties (Keith Doore)

