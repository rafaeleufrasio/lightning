.. Lightning API Documentation documentation master file, created by
   sphinx-quickstart on Wed Mar 16 17:40:45 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Lightning
=======================================================

Lightning is a spectral energy distribution (SED) fitting code originally designed to rapidly fit
non-parametric star formation history (SFH) models to photometric galaxy data.

Now, Lightning employs a variety of methods, including an adaptive Markov-Chain Monte-Carlo,
affine-invariant Markov-Chain Monte-Carlo, and Levenberg-Marquardt gradient decent (MPFIT)
algorithms to fit X-ray to far-infrared (FIR) SEDs. The models fit to the SEDs can have
contributions from stellar populations, ISM dust emission, and AGN.

Documentation
-------------
.. toctree::
   :maxdepth: 1

   installation
   getting_started
   postprocessing
   examples/index
   filters
   FAQs
   API_reference/index

Attribution
-----------
If you make use of Lightning, please cite the following articles:

.. code-block:: latex

        @MISC{2017ascl.soft11009E,
            author = {{Eufrasio}, Rafael T.},
            title = "{Lightning: SED Fitting Package}",
            keywords = {Software},
            year = 2017,
            month = nov,
            eid = {ascl:1711.009},
            pages = {ascl:1711.009},
            archivePrefix = {ascl},
            eprint = {1711.009},
            adsurl = {https://ui.adsabs.harvard.edu/abs/2017ascl.soft11009E},
            adsnote = {Provided by the SAO/NASA Astrophysics Data System}
        }

        @ARTICLE{2017ApJ...851...10E,
            author = {{Eufrasio}, R.~T. and {Lehmer}, B.~D. and {Zezas}, A. and {Dwek}, E. and {Arendt}, R.~G. and {Basu-Zych}, A. and {Wiklind}, T. and {Yukita}, M. and {Fragos}, T. and {Hornschemeier}, A.~E. and {Markwardt}, L. and {Ptak}, A. and {Tzanavaris}, P.},
            title = "{On the Spatially Resolved Star Formation History in M51. I. Hybrid UV+IR Star Formation Laws and IR Emission from Dust Heated by Old Stars}",
            journal = {\apj},
            keywords = {galaxies: individual: NGC 5194, NGC 5195, galaxies: interactions, galaxies: spiral, galaxies: star formation, galaxies: stellar content, Astrophysics - Astrophysics of Galaxies, Astrophysics - Cosmology and Nongalactic Astrophysics},
            year = 2017,
            month = dec,
            volume = {851},
            number = {1},
            eid = {10},
            pages = {10},
            doi = {10.3847/1538-4357/aa9569},
            archivePrefix = {arXiv},
            eprint = {1710.09401},
            primaryClass = {astro-ph.GA},
            adsurl = {https://ui.adsabs.harvard.edu/abs/2017ApJ...851...10E},
            adsnote = {Provided by the SAO/NASA Astrophysics Data System}
        }

An incomplete list of other papers using ``Lightning`` are:

- `Lehmer et al. (2017) <https://ui.adsabs.harvard.edu/abs/2017ApJ...851...11L/abstract>`_
- `Lehmer et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020ApJS..248...31L/abstract>`_
- `Doore et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021ApJ...923...26D/abstract>`_
- `Monson et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021ApJ...919...51M/abstract>`_
- `Motiño Flores et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021ApJ...921..130M/abstract>`_
- `Abdeen et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022MNRAS.512..366A/abstract>`_
- `Doore et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022ApJ...931...53D/abstract>`_
- `Gilbertson et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022ApJ...926...28G/abstract>`_
- `Lehmer et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022ApJ...930..135L/abstract>`_

License
-------
Lightning is available under the terms of the MIT license.


.. Indices and tables
.. ==================
..
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
