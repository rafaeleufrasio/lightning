#### Lightning SED Fitting Package ####



#### Introduction ####
Lightning is a spectral energy distribution (SED) fitting procedure that quickly and reliably recovers
star formation history (SFH) and extinction parameters. The SFH is modeled as discrete steps in time.
The Version 1 (v1) code consists of a fully vectorized inversion algorithm to determine SFH step intensities
and combines this with a grid-based approach to determine three attenuation parameters. The Version 2 (v2)
code consists of an adaptive Markov chain Monte Carlo (MCMC) algorithm to determine SFH step intensities, 
attenuation parameters, and dust emission parameters. Version 2 adds the dust emission models from
Draine & Li (2007) to model NIR to FIR emission, and the inclination-dependent attenuation curves from
Tuffs et al. (2004).

NOTE: We are currently working on compiling the fitting methods (and adding more)
to one streamlined package and writing thorough documentation to be easy to use for new users.
This basic documentation is for the most recent version of Lightning that uses the MCMC algorithm (v2).
However, the package contains both versions.




#### Installation ####
To install the Lightning package, enter the following:

    cd <install_dir>
    git clone https://github.com/rafaeleufrasio/lightning

Lightning is pure IDL and requires the IDL Astronomy User's Library (https://idlastro.gsfc.nasa.gov).




#### Post-Installation, Pre-Usage IDL Set-Up ####
Before running Lightning v2, the Lightning package must be added to you IDL path. This can be done 
permanently within an IDL session so that it does not have to be added each session with:

    PREF_SET, 'IDL_PATH', PREF_GET('IDL_PATH')+':+/PATH_TO_LIGHTNING_DIRECTORY/lightning', /COMMIT

where /PATH_TO_LIGHTNING_DIRECTORY/ is the path to the lightning directory, or only for 
the current session with:

    !PATH = EXPAND_PATH(PREF_GET('IDL_PATH')+':+/PATH_TO_LIGHTNING_DIRECTORY/lightning')

Then the lightning startup file must be run:

    @ startup_lightning.pro

To simplify this process, this line can be added to the users IDL_STARTUP file. If
one does not have a startup file, startup_RTE.pro can be made the startup file and 
added to the users preferences with:

    PREF_SET, 'IDL_STARTUP', /PATH_TO_LIGHTNING_DIRECTORY/startup_RTE.pro, /commit

Lightning v2 must then be compiled before running using:

    .r dl07_spec_vector.pro
    .r lightning_MCMC.pro




#### General Usage and Parameter Documentation ####
CALLING SEQUENCE:
  lightning_MCMC_vector, LNU_OBS, LNU_UNC, LTIR_OBS, LTIR_UNC, GALAXY_ID, FILTER_LABELS, /DUST_EMISSION,$
         [Z_SHIFT=Z_SHIFT, STEPS_BOUNDS=STEPS_BOUNDS, NPARALLEL=NPARALLEL, NTRIALS=NTRIALS, $
          OUTFOLDER=OUTFOLDER, LIGHTNING_FOLDER=LIGHTNING_FOLDER, ZMET=ZMET, IMF=IMF, $
          COEFF_START=COEFF_START, COEFF_SIMGA=COEFF_SIMGA, PARAMETER_START=PARAMETER_START, $
          PARAMETER_SIMGA=PARAMETER_SIGMA, BETA_EXPONENT=BETA_EXPONENT, PRIOR_DIST=PRIOR_DIST, $
          L_STAR_ABS_TABLE=L_STAR_ABS_TABLE, /PRINT_TIME, /NOLINES, /NONEBULAR, /PAR1_CONSTANT $
          /PAR2_CONSTANT, /PAR3_CONSTANT, /PAR4_CONSTANT, /PAR5_CONSTANT, /PAR6_CONSTANT, $
          /PAR7_CONSTANT, /PAR8_CONSTANT, /PAR9_CONSTANT, /PAR10_CONSTANT, /ADAPTIVE, $
          /CALZETTI_ATTENUATION, /TUFFS_ATTENUATION, /USE_PRIORS, /L_STAR_ABS_MODEL_TABLE]

REQUIRED INPUTS:
  LNU_OBS   - The observed photometry to be fit as Lnu with units of L_sun Hz^-1.
                This requires the observed flux to be converted to luminosity.
  LNU_UNC   - The observed photometry's uncertainty as Lnu with units of L_sun Hz^-1
  LTIR_OBS  - The observed total infrared luminosity. With the usage of the dust model, 
                this should be set to 0.d0
  LTIR_UNC  - The uncertainty on the total infrared luminosity. With the usage of the dust
                model, this should be set to 0.d0
  GALAXY_ID - The ID string label for the SED being fit (e.g., 'NGC224' or 'J004244.30+411609.00')
  FILTER_LABELS - The string labels for the filters corresponding to the observed photometry 
                    (LNU_OBS) in the same order. Filters' transmission functions must be stored 
                    within the Filters subdirectory of the Lightning package. The labels used 
                    MUST match those found listed in the case statement of get_filters.pro. 
                    If the filter is not in the case statement of get_filters.pro, then the 
                    name of the filter must be added to the procedure and the transmission 
                    function must be added to the Filters subdirectory following the layout 
                    of all other filter transmission functions.

OPTIONAL INPUTS:
  UNIVERSAL OPTIONAL INPUTS:
    Z_SHIFT          - Redshift of SED.
                        (Default is 0.0) 
    STEPS_BOUNDS     - Vector giving the boundaries of the age bins in years.
                        (Default is [0.d0,1.d7,1.d8,1.d9,5.d9,13.6d9])
    NPARALLEL        - Number of chains to run in one call of Lightning. This will only
                         fit the input SED multiple times, not multiple SEDs at once. Mainly
                         used to check for convergence of the MCMC chain.
                         (Default is 1, which means fit the SED once)
    NTRIALS          - Number of iterations to create the MCMC chain
                         (Default is 1e5)
    OUTFOLDER        - String containing the directory for the output files
                         (Default is '~/MCMC_Lightning_runs/')
    LIGHTNING_FOLDER - String containing the path to and including the installed Lightning 
                         package. (i.e., '/PATH_TO_LIGHTNING_DIRECTORY/lightning/')
                         (Default is '~/lightning/')
    ZMET             - Metallicity of to fit the SED in terms of Z. Currently available 
                         metallicites are 0.001, 0.004, 0.008, 0.02, 0.05, 0.1.
                         (Default is 0.02, which is solar metallicity)
    IMF              - String containing the IMF to fit the SED. Currently available IMFs
                         are Kroupa (2001).
                         (Default is Kroupa (2001))
    COEFF_START      - Vector of length n_elements(STEPS_BOUNDS)-1 containing the starting
                         MCMC values for the SFH coefficients. If NPARALLEL > 1, then 
                         COEFF_START can be an array with dimensions 
                         [n_elements(STEPS_BOUNDS)-1,NPARALLEL], or it can remain a vector,
                         then each parallel chain will have the same starting point.
                         (Default is a vector of all 1.0's)
    COEFF_SIGMA      - Vector of length n_elements(STEPS_BOUNDS)-1 containing the starting
                         standard deviations to use for the SFH coefficients in the proposal 
                         covariance matrix. If NPARALLEL > 1, then COEFF_SIMGA can be an array 
                         with dimensions [n_elements(STEPS_BOUNDS)-1,NPARALLEL], or it can 
                         remain a vector, then each parallel chain will have the same starting 
                         covariance matrix.
                         (Default is a vector of all 1.0's)
                         We recommend NOT using this optional input if using the Adaptive MCMC
                         algorithm.
    PARAMETER_START  - A 10 element vector containing the starting attenuation and dust emission
                         parameters. If NPARALLEL > 1, then PARAMETER_START can be an array 
                         with dimensions [10,NPARALLEL], or it can remain a vector,
                         then each parallel chain will have the same starting parameters.
                         
                         The first 5 elements of the vector are the dust attenuation parameters,
                         and the final 5 are the dust emission parameters.
                         
                         If the Calzetti or modified Calzetti curve is used, then the first 
                         parameter is the diffuse V-band optical depth (tauV_diff, range:[0.d0,3.d0]),
                         the second is the variable slope (delta, range:[-2.3d,0.4d]), the 
                         third is the birth cloud optical depth in the V-band 
                         (tauV_BC, range:[0.0d,4.0d]), and the final two attenuation parameters are 
                         void and can be any value (preferably 0.0).
                                               
                         If the Tuffs attenuation curves are used, then the first parameter is 
                         the face-on B-band optical depth (taub_f, range:[0.d0,8.d0]),
                         the second is fraction of attenuation coming from the disk (rdisk, range:[0.0d,1.0d]), 
                         the third is the birth cloud clumpiness factor (F, range:[0.0d,0.61d]), 
                         the fourth is the inclination in terms of cos(i) (cosi, range:[0.0d,1.0d]),
                         and the fifth is the fraction of attenuation coming from the bulge
                         (rbulge, range:[0.0d,1.0d]). 
                         
                         If the dust emission is used, then the sixth parameter (first dust emission 
                         parameter) is the slope of the power law (alpha, range:[-10.d0,4d0]),
                         the seventh is the minimum radiation field intensity (Umin, range:[0.1d0,25]),
                         the eighth is the maximum radiation field intensity (Umax, range:[1.d3,3.d5]),
                         the ninth is the fraction of dust mass exposed to the power-law radiation field
                         (gamma, range:[0.0d,1.0d]), and the tenth is the PAH index (qPAH, range:[0.0047d,0.0458d]).
                        
                         (Default is [0.2,0.0,0.2,0,0,-2.d0,1.0,1.d4,0.1,0.020] for the Calzetti
                         curve and [1.0,0.2,0.2,0.8,0.2,-2.d0,1.0,1.d4,0.1,0.020] for the Tuffs) 
    PARAMETER_SIGMA  - A 10 element vector containing the starting standard deviations to use for 
                         the attenuation and dust emission parameters in the proposal covariance 
                         matrix. If NPARALLEL > 1, then PARAMETER_START can be an array 
                         with dimensions [10,NPARALLEL], or it can remain a vector,
                         then each parallel chain will have the same starting parameters.
                         
                         The first 5 elements of the vector are the dust attenuation parameters,
                         and the final 5 are the dust emission parameters.

                         (Default is [0.1,0.1,0.2,0,0,0.5d0,0.1,1.d4,0.1,0.005] for the Calzetti
                         curve and [0.5,0.1,0.1,0.1,0.1,0.5d0,0.1,1.d4,0.1,0.005] for the Tuffs) 
                         We recommend NOT using this optional input if using the Adaptive MCMC
                         algorithm.

  ADAPTIVE MCMC OPTIONAL INPUTS:
   (ADAPTIVE keyword must be set for these to be used) 
    BETA_EXPONENT - Value which dictates how fast the adaptive MCMC stops adapting. Larger
                      values decrease the adaptiveness (stop adapting sooner), smaller values
                      increase the adaptiveness (adapt for longer). Must be between 0 and 1.
                      (Default is 0.35)
  
  TUFFS ATTENUATION OPTIONAL INPUTS:         
   (TUFFS_ATTENUATION keyword must be set for these to be used) 
    PRIOR_DIST       - Structure giving the tabulated prior distribution and corresponding 
                         bins from which to interpolate a prior. Current version of Lightning
                         only supports a prior on inclination. PRIOR_DIST must have the following
                         tags: inc_pdf and inc_bins, which are the probability of the inclination
                         (inc_pdf) at the corresponding inclination (inc_bins). USE_PROIR keyword
                         must be set to use prior.
                         (Default assumes a flat prior)
    L_STAR_ABS_TABLE - Structure giving the pre-computed table of L_STAR_ABS for a given redshift
                         to have energy conservation with the Tuffs attenuation curves. Not required
                         as Lightning will automatically read the table if the keyword 
                         L_STAR_ABS_MODEL_TABLE is set and the default STEPS_BOUNDS is used. 
                         However, pre-reading the table and inputting it will vastly increase 
                         speed of the code as Lightning re-reads the table for each MCMC iteration.
                         (This will be remedied in the next version)
                         (Structure for default STEPS_BOUNDS can be found in
                         ./lightning/L_star_abs_model_table/)

OPTIONAL KEYWORDS:
  UNIVESAL KEYWORDS:
    PRINT_TIME     - if set, then Lightning will print the time taken to fit the SED   
    NOLINES        - if set, then no emission lines will be included in the stellar models  
    NONEBULAR      - if set, then no nebular emission will be included in the stellar models
    PAR1_CONSTANT  - if set, then parameter 1 (tauV_diff for Calzetti models, taub_f for 
                       Tuffs model) will be held constant      
    PAR2_CONSTANT  - if set, then parameter 2 (delta for Calzetti models, rdisk for 
                       Tuffs model) will be held constant       
    PAR3_CONSTANT  - if set, then parameter 3 (tauV_BC for Calzetti models, F for 
                       Tuffs model) will be held constant       
    PAR4_CONSTANT  - if set, then parameter 4 (cosi for Tuffs model) will be held constant      
    PAR5_CONSTANT  - if set, then parameter 5 (rbulge for Tuffs model) will be held constant      
    PAR6_CONSTANT  - if set, then parameter 6 (alpha for dust model) will be held constant      
    PAR7_CONSTANT  - if set, then parameter 7 (Umin for dust model) will be held constant      
    PAR8_CONSTANT  - if set, then parameter 8 (Umax for dust model) will be held constant      
    PAR9_CONSTANT  - if set, then parameter 9 (gamma for dust model) will be held constant      
    PAR10_CONSTANT - if set, then parameter 10 (qPAH for dust model) will be held constant      
  
  DUST EMISSION KEYWORDS
    DUST_EMISSION  - if set, then the Draine & Li (2007) dust model will be used to fit the 
                       NIR to FIR emission. 
                       NOTE: This is a REQUIRED keyword in this version of
                       Lightning, without setting it, the code will fail! 

  ADAPTIVE MCMC KEYWORDS:
    ADAPTIVE  - if set, then the MCMC algorithm will be adaptive, automatically finding the
                  best proposal multivariate distribution to optimize the acceptance rate

  CALZETTI ATTENUATION KEYWORDS:         
    CALZETTI_ATTENUATION - if set, then the attenuation curve will be a pure Calzetti law.
                             If this or TUFFS_ATTENUATION are not set, then the attenuation
                             curve will be the modified Calzetti curve as described in Noll
                             et al. (2009).       

  TUFFS ATTENUATION KEYWORDS:         
    TUFFS_ATTENUATION      - if set, then the attenuation cuve will be the inclination-
                               dependent Tuffs attenuation curves. If this or CALZETTI_ATTENUATION
                               are not set, then the attenuation curve will be the modified 
                               Calzetti curve as described in Noll et al. (2009).      
    USE_PRIORS             - if set, then a prior will be used on inclination. This must be
                               set in conjunction with PRIOR_DIST.      
    L_STAR_ABS_MODEL_TABLE - if set, then energy balance will be determined from the pre-computed
                               tables of L_STAR_ABS_MODEL_TABLE. If not set, then L_STAR_ABS will
                               be calculated with each MCMC iteration, drastically increasing
                               computational time. 

OUTPUTS (NOTE: Output is a save file (.idl) in the OUTFOLDER directory with the name
         GALAXY_ID+'_chain.idl' containing the following variables):
  CHAIN            - An array with dimensions of [n_elements(STEPS_BOUNDS)-1+10,NTRIALS,NPARALLEL].
                       This is the resulting MCMC chain for each parameter and parallel chain.
                       The first n_elements(STEPS_BOUNDS)-1 in the first dimension are the SFH
                       coefficients, the next five are the attenuation parameters, and the final
                       five are the dust emission parameters.
  CHI2_CHAIN       - An array with dimensions of [NTRIALS,NPARALLEL] giving the chisqr of the fit
                       of the model to the data for each chain iteration and each parallel chain. 
  GALAXY_ID        - The ID of the galaxy (Same as input GALAXY_ID)
  SIGMA_NEW        - Final covariance matrix
  LAMBDA_NEW       - Final factor to multiply the covariance matrix by
  MU_NEW           - Final average of the chain
  STEPS_MSTAR      - A n_elements(STEPS_BOUNDS)-1 vector giving the factors needed to convert 
                       the SFH coefficients into stellar masses. To get a MCMC chain of masses
                       for each time bin do the following:
                       
                          steps_mstar_chain=chain[0:n_elements(STEPS_BOUNDS)-2,*,*]*rebin(steps_mstar,5,NTRIALS,NPARALLEL)
                       
                       To get the total mass, sum in the SFH coefficient dimension or:
                       
                          mstar_chain=total(steps_mstar_chain,1,/NAN)
                       
  STEPS_BOUNDS     - Vector giving the boundaries of the age bins in years. (Same as input STEPS_BOUNDS)
  ACCEPTED_TRIALS  - Number of accepted MCMC trials.
  CHI2_PRIOR_CHAIN - An array with dimensions of [NTRIALS,NPARALLEL] giving the chisqr 
                       from the prior for each chain iteration and each parallel chain. 
  
NOTES:
  Lightning does not save the output MCMC data into a well structured .fits file. It instead
  saves it as the raw MCMC chain. To compute most parameters of interest, the chain will have
  to be input into the MCMC_savefits_X.pro file, where X is the attenuation curve used. We do
  not document this process here, but will document and streamline this in the upcoming version
  of Lightning.
  



#### Example ####
An example code for fitting the SED of galaxy J123548.94+621144.7 can be found in the Example
subdirectory of the Lightning package.

To run the example, enter the following into the IDL command line after changing the
current directory to the Lightning directory (cd ./lightning/):

    @ ./Example/lightning_example.pro

This example results in two output files. One that contains the fitting results to the SED of J123548.94+621144.7
using the Calzetti et al. (2000) attenuation curve, and the other that contains the results using 
the Tuffs et al. (2004) inclination-dependent attenuation curves. Fits files for the parameters of interest
have been computed. The corresponding fits files (MCMC_lightning_calz.fits and MCMC_lightning_tuffs.fits, 
respectively) contain the final 5,000 iterations of the MCMC chain for various parameters of interest, 
such as the SFH, stellar mass, and AV; and the attenuation and dust emission parameters.




#### Citation ####
If you use this code please cite

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

along with

@ARTICLE{2017ApJ...851...10E,
       author = {{Eufrasio}, R.~T. and {Lehmer}, B.~D. and {Zezas}, A. and {Dwek}, E. and
         {Arendt}, R.~G. and {Basu-Zych}, A. and {Wiklind}, T. and {Yukita}, M. and
         {Fragos}, T. and {Hornschemeier}, A.~E. and {Markwardt}, L. and
         {Ptak}, A. and {Tzanavaris}, P.},
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
