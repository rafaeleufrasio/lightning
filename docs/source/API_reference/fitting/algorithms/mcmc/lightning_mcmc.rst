LIGHTNING_MCMC
==============

Name
----
LIGHTNING_MCMC

Purpose
-------
Fits an observed SED with the Lightning models using either the vanishing
adaptive Metropolis-Hastings algorithm (Algorithm 4 from Andrieu & Thoms 2008),
or the affine-invariant Goodman & Weare (2010) algorithm. Initial parameter
positions are randomly selected from the specified prior distribution.

Calling Sequence
----------------
::

    lightning_mcmc, input_dir, sed_data, config_nopriors, models, priors

Inputs
------
``input_dir`` : string scalar
    The path to the file containing the input SED data.
``sed_data`` : structure
    A structure containing the SED luminosities and uncertainties, filter
    labels, distances, redshifts, and optional X-ray data. (See
    ``lightning_input.pro`` for details and contents.)
``config_nopriors`` : structure
    A Lightning configuration structure edited to remove the prior
    substructures for each parameter. (See ``lightning_configure_defaults.pro``
    for details and contents.)
``models`` : structure
    A structure containing each model structure (stellar, dust, AGN,
    X-ray) as a substructure. (See ``lightning_models.pro`` for details
    and contents.)
``priors`` : structure
     A structure containing the prior hyper-parameters. (See
     ``generate_prior_struct.pro`` for details and contents.)

Output
------
An IDL save file saved to ``<input_dir>/lightning_output/output_sav_files/`` named
``lightning_output_<galaxy_id>.sav'``, containing the resulting MCMC parameter
and log probability chain(s), and the convergence metrics of the chain(s).

Notes
-----
The Goodman & Weare affine invariant algorithm requires an ensemble of *at least* ``Ndim + 1`` walkers, which all
must be randomly initialized such that the ``[Ndim, Nparallel]`` matrix describing their positions is
nonsingular. These walkers are not independent and should not be treated as such in postprocessing: convergence
metrics such as Brooks-Gelman and Gelman-Rubin tests are not applicable, and the chains should be flattened into a
single vector to derive the posterior distributions. We recommend using the acceptance fraction (ideally 20-50%)
and the autocorrelation time (ideally much shorter than the length of the chains for all dimensions)
to assess the performance of this method.

References
----------
- `Andrieu, C., & Thoms, J. 2008, Statistics and Computing, 18, 30 <https://link.springer.com/article/10.1007/s11222-008-9110-y>`_
- `Goodman, J., & Weare, J. 2010, CAMCS, 5, 65 <https://ui.adsabs.harvard.edu/abs/2010CAMCS...5...65G/abstract>`_

Modification History
--------------------
- 2020/04/27: Replaced if statements with ``n_elements`` on keywords to use ``keyword_set`` (Keith Doore)
- 2020/05/06: Set default starting attenuation parameters to ``0.0`` so if not wanted they would not be used (Keith Doore)
- 2020/05/06: Added ability to use pure Calzetti curve (no changes here due to ``_extra``) (Keith Doore)
- 2020/05/06: Added ``_ref_extra`` for keyword inheritance to cut down on list of keywords (Keith Doore)
- 2020/05/06: Changed all keywords that were the same to match across functions/procedures (Keith Doore)
- 2020/05/06: Removed any repetitive items at beginning that are set by other functions if not set in MCMC call (Keith Doore)
- 2020/05/06: Added needed items to run Tuffs attenuation as to match other attenuation (Keith Doore)
- 2021/03/17: Added UV-to-IR AGN fitting (Erik Monson)
- 2021/04/13: Modified Tuffs attenuation to have ``rdisk`` as intrinsic property and ``rbulge`` as B/D (Keith Doore)
- 2021/04/16: Added X-ray fitting (Erik Monson)
- 2021/04/16: Change sampling of AGN angle parameters: ``i`` -> ``cosi``; ``oa`` -> ``sin(oa)`` (Erik Monson)
- 2022/04/23: Changed ``lightning_mcmc.pro`` to use the new modular Lightning framework (Erik Monson)
- 2022/04/23: Added documentation (Erik Monson)
- 2022/06/07: Major update to include new implementation (e.g., prior, config, etc.) (Keith Doore)
- 2022/07/05: Renamed ``config`` to ``config_nopriors`` to reflect removed prior substructures (Keith Doore)
- 2022/07/06: Changed how ``sigma_new`` is initialized to eliminate for loop (Keith Doore)
- 2022/07/19: Added compression to output MCMC sav file to reduce memory usage (Keith Doore)
- 2022/07/19: Updated how ``mu_new`` is initialized to prevent bugs if first trial is rejected (Keith Doore)
- 2022/08/18: Added progress printing (Keith Doore)
- 2022/09/01: Added ``chi2_chain`` to output file in addition to ``lnprob_chain`` (Erik B. Monson)

