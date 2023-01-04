# Lightning SED Fitting Package
[![Documentation Status](https://readthedocs.org/projects/lightning-sed/badge/?version=latest)](https://lightning-sed.readthedocs.io/en/latest/?badge=latest)

## Description

Lightning is a spectral energy distribution (SED) fitting code originally designed to rapidly fit non-parametric star formation history (SFH) models to photometric galaxy data.

Now, Lightning employs a variety of methods, including an adaptive Markov-Chain Monte-Carlo, affine-invariant Markov-Chain Monte-Carlo, and Levenberg-Marquardt gradient decent (MPFIT) algorithms to fit X-ray to far-infrared (FIR) SEDs. The models fit to the SEDs can have contributions from stellar populations, ISM dust emission, and AGN.


## Installation

To install the Lightning package, you can use `git` as follows:

```
cd <install_dir>
git clone https://github.com/rafaeleufrasio/lightning
```

> NOTE: Lightning is written in the Interactive Data Language (IDL) and requires IDL version 8.3 or later. It additionally requires the [IDL Astronomy User's Library](https://idlastro.gsfc.nasa.gov), [IDL Coyote](http://www.idlcoyote.com), and
[Craig Markwardt's MPFIT library](http://purl.com/net/mpfit). All three will need to be installed and [added to your IDL path](https://www.l3harrisgeospatial.com/Support/Self-Help-Tools/Help-Articles/Help-Articles-Detail/ArtMID/10220/ArticleID/16156/Quick-tips-for-customizing-your-IDL-program-search-path) before using Lightning.

Lightning will then need to be built by running:

```
cd <install_dir>/lightning
idl lightning_build
```


## How to Use & Documentation

To learn how to run Lightning and for a helpful user guide, please read our [Online Documentation](https://lightning-sed.readthedocs.io/en/latest/).

## Citation

If you make use of Lightning, please cite the following articles:


```
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
```

## License

Lightning is available under the terms of the MIT license.
