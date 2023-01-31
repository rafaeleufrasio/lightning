function lightning_configure
;+
; Name
; ----
;   LIGHTNING_CONFIGURE
;
; Purpose
; -------
;   Generates the Lightning configuration structure. This structure
;   contains all of the possible selections users can make when
;   configuring Lightning (e.g., models, fitting algorithms, etc.).
;   The user can change the values for each tag to their desired values.
;   Small guiding comments are included. These comments give basic details
;   and indicate the expected values and type. Further details can be
;   found at :ref:`configure-setting-label`.
;
; Calling Sequence
; ----------------
;   ::
;
;       config = lightning_configure()
;
; Output
; ------
;   ``config`` : structure
;       A Lightning configuration structure.
;
; Notes
; -----
;   - Parameters that take a flag (value of ``0`` or ``1``) are set if the value is ``1``
;     and not set if the value is ``0``.
;   - Options that take arrays are indicated by brackets (i.e., ``[]``). These options
;     can have multiple values. Options without brackets must contain a single value.
;     Removing brackets may result in errors when the configuration structure is
;     checked for errors.
;   - All free parameters of the models have an associated prior structure, where
;     the distribution type, prior distribution shape arguments, and initialization
;     range are given. The prior distribution type options are: ``'fixed'``, ``'uniform'``,
;     ``'normal'``, and ``'tabulated'``. The number of values in the distribution
;     shape argument array depend on the chosen distribution type.
;
;     - ``'fixed'``: takes a single value (``Narg = 1``), the value at which to fix the
;       parameter. (If the prior is ``'fixed'``, then the initialization range is ignored,
;       since a fixed value does not need initialization.)
;     - ``'uniform'``: takes two values (``Narg = 2``), the minimum and maximum bounds of the
;       distribution in that order.
;     - ``'normal'``: takes four values (``Narg = 4``), the minimum bound, maximum bound,
;       distribution peak, and distribution standard deviation in that order.
;     - ``'tabulated'``: takes one value (``Narg = 1``), a string containing the path to the
;       directory containing the user tabulated prior file (see :ref:`tabulated-prior-label`
;       for more details).
;
;     The initialization range always has two values per prior,
;     indicating the minimum and maximum bounds for the random initialization of
;     the fitting algorithm.
;-

config = {                             $

;========================================    CORE    =========================================================
; Core options to Lightning


          ; OUTPUT_FILENAME`` : string scalar
          ;   The name (without the file extension suffix) to give to the FITS file containing
          ;   the output post-processed data.
          ; NOTE - A UTC timestamp can be automatically included in the filename so that you
          ;        can have a unique filename for multiple repeat runs to prevent accidentally
          ;        overwriting old runs. This is done by including a single ``%`` character in
          ;        the filename where you want the timestamp to appear. 
          OUTPUT_FILENAME: 'postprocessed_data_%'     ,$

          ; PRINT_PROGRESS : flag (0 or 1)
          ;   If set, the SED fitting progress and expected fitting time
          ;   remaining will be printed to the terminal.
          PRINT_PROGRESS: 1                           ,$

          ; MAX_CPUS : int scalar
          ;   The maximum number of CPUs to utilize. If this value exceeds the
          ;   number of CPUs on the machine, then all CPUs will be used.
          MAX_CPUS: 1                                 ,$

          ; ENERGY_BALANCE : flag (0 or 1)
          ;   If set, the total integrated IR luminosity (normalization) of
          ;   the dust emission is tied to the total absorbed stellar
          ;   (and, if set, AGN) emission.
          ENERGY_BALANCE: 1                           ,$

          ; MODEL_UNC : int, float, or double scalar
          ;   The fractional model uncertainty to use in all bands.
          MODEL_UNC: 0.d0                             ,$

          ; COSMOLOGY : H0 :      int, float, or double scalar
          ;             OMEGA_M : int, float, or double scalar
          ;             LAMBDA0 : int, float, or double scalar
          ;             Q0 :      int, float, or double scalar
          ;             K :       int, float, or double scalar
          ;   The cosmology parameters to use in the SED fitting. These parameters
          ;   set the age of the universe (and distance to the object if the object
          ;   has its distance specified by redshift). The parameters are H0
          ;   (Hubble constant in km/s/Mpc), Omega_m (Matter density, normalized
          ;   to the closure density), Lambda0 (Cosmological constant, normalized to
          ;   the closure density), q0 (Deceleration parameter), and k (Curvature
          ;   constant, normalized to the closure density).
          COSMOLOGY: {H0:      70.                    ,$
                      OMEGA_M: 0.3                    ,$
                      LAMBDA0: 0.7                    ,$
                      Q0:     -0.55                   ,$
                      K:       0.0                   },$




;==================================    STELLAR EMISSION    ===================================================
; Stellar emission module options

          ; SSP : string scalar
          ;   The simple stellar population (SSP) models to use for
          ;   the stellar population. For no stellar emission model, set
          ;   to 'NONE'.
          ;   Current options: 'PEGASE' and 'NONE'
          ; NOTE - If no stellar emission model is chosen, all stellar
          ;        emission model parameters will be ignored.
          SSP: 'PEGASE'                               ,$


          ;=========== STELLAR POPULATION ===========
          ; IMF : string scalar
          ;   The initial mass function (IMF) to use in the SSP models.
          ;   Current options: 'Kroupa01'
          IMF: 'Kroupa01'                             ,$

          ; ZMETAL : float or double scalar
          ;   The metallicity to use in the SSP models in terms of Z [Zsun = 0.02].
          ;   Current options: 0.001, 0.004, 0.008, 0.02, 0.05, 0.1
          ZMETAL: 0.02                                ,$

          ; EMISSION_LINES : flag (0 or 1)
          ;   If set, nebular emission lines are included in the SSP models.
          EMISSION_LINES: 1                           ,$

          ; NEBULAR_EXTINCTION : flag (0 or 1)
          ;   If set, nebular extinction is included in the SSP models.
          NEBULAR_EXTINCTION: 1                       ,$

          ; SFH : string scalar
          ;   The type of SFH to assume when fitting the SEDs.
          ;   Current options: 'NON-PARAMETRIC'
          SFH: 'NON-PARAMETRIC'                       ,$


          ;===================== NON-PARAMETRIC SFH =======================================
              ; STEPS_BOUNDS : int, float, or double array(Nsteps+1)
              ;   The age bin (or step) boundaries to use in the SFH in yr. Values must
              ;   be in ascending order.
              ; NOTE - If an age bin contains ages older than the universe at an input
              ;        SED's redshift, the age bin upper bound will be automatically
              ;        adjusted to the age of the universe at that redshift. If an entire
              ;        age bin is older than universe at that redshift, then the entire
              ;        age bin will be omitted and the next younger bin will be adjusted
              ;        accordingly.
              STEPS_BOUNDS: [0.d0, 1.d7, 1.d8, 1.d9, 5.d9, 13.6d9]  ,$

              ; DTIME_SF : int, float, or double scalar
              ;   The time step used for interpolating the SSP models into the age
              ;   bins in yr.
              ; NOTE - We do not recommend changing this value from its default,
              ;        unless you specified age bins with differences less than
              ;        the default value.
              DTIME_SF: 5.d5                                        ,$

              ; PSI : Prior : string array(Nsteps)
              ;       Prior_arg : int, float, or double array(Nsteps, Narg)
              ;                   or string array(Nsteps)
              ;       Initialization_range : int, float, or double array(Nsteps, 2)
              ;   The SFR for all SFH age bins in Msun yr-1. The number of elements
              ;   must be one less than the number of elements in STEPS_BOUNDS.
              ;   Allowed parameter range: 0 to Infinity.
              ; NOTE - PSI contains the priors for all SFH age bins. The prior
              ;        distribution types, arguments, and initialization range
              ;        can be unique for each age bin. This excludes the prior type
              ;        of 'tabulated', which all or none have to be.
              ;      - If the chosen priors types have different corresponding values
              ;        of Narg, then the array size uses the largest Narg. For priors with
              ;        smaller Narg values, they will need padded values (which
              ;        will be ignored) included where the additional arguments of
              ;        the more complex priors are required (see the documentation
              ;        for more details and examples).
              PSI :      {Prior: ['uniform', 'uniform', 'uniform', 'uniform', 'uniform'] ,$
                          Prior_arg: [[0.0d, 0.0d, 0.0d, 0.0d, 0.0d]                     ,$
                                      [1.d3, 1.d3, 1.d3, 1.d3, 1.d3]]                    ,$
                          Initialization_range: [[0.0d, 0.0d, 0.0d, 0.0d, 0.0d]          ,$
                                                 [1.d1, 1.d1, 1.d1, 1.d1, 1.d1]]          $
                          }                                                              ,$




;==================================    DUST ATTENUATION    ===================================================
; Dust attenuation module options

          ; ATTEN_CURVE : string scalar
          ;   The assumed attenuation curve applied to the stellar and/or AGN models.
          ;   Current options: 'CALZETTI00', 'CALZETTI_MOD', and 'DOORE21'
          ; NOTE - Only configurations of the chosen attenuation curve should be
          ;        updated from the default. Changes to the configurations
          ;        of the other attenuation curves will be ignored.
          ; NOTE - Attenuation of AGN can only use the 'CALZETTI00' or 'CALZETTI_MOD'
          ;        attenuation curves. Compatibility with the 'DOORE21' curve is
          ;        currently not supported.
          ATTEN_CURVE: 'CALZETTI00'                   ,$


          ;=============== CALZETTI00 PARAMETERS ===========================
              ; TAUV : Prior : string scalar
              ;        Prior_arg : int, float, double, or string array(Narg)
              ;        Initialization_range : int, float, or double array(2)
              ;   The V-band optical depth.
              ;   Allowed parameter range: 0 to Infinity.
              TAUV:      {Prior: 'uniform'                          ,$
                          Prior_arg: [0.0d, 10.0d]                  ,$
                          Initialization_range: [0.0d, 3.0d]         $
                          }                                         ,$


          ;============== CALZETTI_MOD PARAMETERS =========================
              ; TAUV_DIFF : Prior : string scalar
              ;             Prior_arg : int, float, double, or string array(Narg)
              ;             Initialization_range : int, float, or double array(2)
              ;   The V-band optical depth of diffuse dust.
              ;   Allowed parameter range: 0 to Infinity.
              TAUV_DIFF: {Prior: 'uniform'                          ,$
                          Prior_arg: [0.0d, 10.0d]                  ,$
                          Initialization_range: [0.0d, 3.0d]         $
                          }                                         ,$

              ; DELTA : Prior : string scalar
              ;         Prior_arg : int, float, double, or string array(Narg)
              ;         Initialization_range : int, float, or double array(2)
              ;   The power law value to change the attenuation curve slope.
              ;   Allowed parameter range: -Infinity to Infinity.
              DELTA:     {Prior: 'uniform'                          ,$
                          Prior_arg: [-2.3d, 0.4d]                  ,$
                          Initialization_range: [-1.0d, 0.0d]        $
                          }                                         ,$

              ; TAUV_BC : Prior : string scalar
              ;           Prior_arg : int, float, double, or string array(Narg)
              ;           Initialization_range : int, float, or double array(2)
              ;   The V-band optical depth of the birth cloud component.
              ;   Allowed parameter range: 0 to Infinity.
              TAUV_BC:   {Prior: 'fixed'                            ,$
                          Prior_arg: [0.0d]                         ,$
                          Initialization_range: [0.0d, 3.0d]         $
                          }                                         ,$

              ; UV_BUMP : flag (0 or 1)
              ;   If set, a 2175 Angstrom UV bump feature will be added
              ;   to the attenuation curve.
              UV_BUMP:   1                                          ,$


          ;================  DOORE21 PARAMETERS =============================
              ; TAUB_F : Prior : string scalar
              ;          Prior_arg : int, float, double, or string array(Narg)
              ;          Initialization_range : int, float, or double array(2)
              ;   The face-on optical depth in the B-band.
              ;   Allowed parameter range: 0 to 8.
              TAUB_F:    {Prior: 'uniform'                          ,$
                          Prior_arg: [0.0d, 8.0d]                   ,$
                          Initialization_range: [2.0d, 6.0d]         $
                          }                                         ,$

              ; F_CLUMP : Prior : string scalar
              ;           Prior_arg : int, float, double, or string array(Narg)
              ;           Initialization_range : int, float, or double array(2)
              ;   The birth cloud clumpiness factor, F.
              ;   Allowed parameter range: 0 to 0.61.
              F_CLUMP:   {Prior: 'uniform'                          ,$
                          Prior_arg: [0.0d, 0.61d]                  ,$
                          Initialization_range: [0.0d, 0.61d]        $
                          }                                         ,$

              ; COSI : Prior : string scalar
              ;        Prior_arg : int, float, double, or string array(Narg)
              ;        Initialization_range : int, float, or double array(2)
              ;   The inclination of the galactic disk in terms of cos(i).
              ;   Allowed parameter range: 0 to 1.
              COSI:      {Prior: 'uniform'                          ,$
                          Prior_arg: [0.0d, 1.0d]                   ,$
                          Initialization_range: [0.0d, 1.0d]         $
                          }                                         ,$

              ; B_TO_D : Prior : string scalar
              ;          Prior_arg : int, float, double, or string array(Narg)
              ;          Initialization_range : int, float, or double array(2)
              ;   The bulge-to-disk ratio.
              ;   Allowed parameter range: 0 to Infinity.
              B_TO_D:    {Prior: 'uniform'                          ,$
                          Prior_arg: [0.0d, 1.d3]                   ,$
                          Initialization_range: [0.0d, 5.d1]         $
                          }                                         ,$

              ; ROLD0_AGES : int, float, or double array(Nsteps)
              ;   The binary parameter rold0, designating each SFH age bin as part of
              ;   the young or old population. A value of 0 for the corresponding age
              ;   bin considers it to be part of the young population, and a value of
              ;   1 considers it to be part of the old populations (see Doore et al.
              ;   2021 for more details). The number of elements must be one less than
              ;   the number of elements in STEPS_BOUNDS.
              ; NOTE - Step bins that contain ages < 500 Myr should be
              ;        considered part of the young population as they
              ;        can contain significant UV emission. Set step bins
              ;        with ages < 500 Myr to the old population at your
              ;        own risk.
              ROLD0_AGES:  [0, 0, 0, 1, 1]                          ,$





;===================================    DUST EMISSION    =====================================================
; Dust emission module options

          ; DUST_MODEL : string scalar
          ;   The dust emission model to use. For no dust emission model,
          ;   set to 'NONE'.
          ;   Current options: 'DL07' and 'NONE'
          ; NOTE - If no dust emission model is chosen, all dust emission
          ;        model parameters will be ignored.
          DUST_MODEL: 'DL07'                          ,$


          ;============== DL07 PARAMETERS ===================================
              ; UMIN : Prior : string scalar
              ;        Prior_arg : int, float, double, or string array(Narg)
              ;        Initialization_range : int, float, or double array(2)
              ;   The minimum radiation field intensity of the diffuse ISM radiation
              ;   field heating the dust.
              ;   Allowed parameter range: 0.1 to 25.
              UMIN:    {Prior: 'uniform'                            ,$
                        Prior_arg: [0.1d, 25.d]                     ,$
                        Initialization_range: [0.1d, 10.d]           $
                        }                                           ,$

              ; UMAX : Prior : string scalar
              ;        Prior_arg : int, float, double, or string array(Narg)
              ;        Initialization_range : int, float, or double array(2)
              ;   The maximum radiation field intensity of the power-law distribution
              ;   of heating starlight intensities.
              ;   Allowed parameter range: 10^3 to 3x10^5.
              UMAX:    {Prior: 'fixed'                              ,$
                        Prior_arg: [3.d5]                           ,$
                        Initialization_range: [1.d5, 3.d5]           $
                        }                                           ,$

              ; ALPHA : Prior : string scalar
              ;         Prior_arg : int, float, double, or string array(Narg)
              ;         Initialization_range : int, float, or double array(2)
              ;       The exponent of the power-law distribution of heating starlight
              ;       intensities between Umin and Umax.
              ;   Allowed parameter range: -10 to 4.
              ALPHA:   {Prior: 'fixed'                              ,$
                        Prior_arg: [2.0d]                           ,$
                        Initialization_range: [1.0d, 3.0d]           $
                        }                                           ,$

              ; GAMMA : Prior : string scalar
              ;         Prior_arg : int, float, double, or string array(Narg)
              ;         Initialization_range : int, float, or double array(2)
              ;   The fraction of the dust mass exposed to the power-law
              ;   distribution of radiation field intensities.
              ;   Allowed parameter range: 0 to 1.
              GAMMA:   {Prior: 'uniform'                            ,$
                        Prior_arg: [0.0d, 1.0d]                     ,$
                        Initialization_range: [0.0d, 0.5d]           $
                        }                                           ,$

              ; QPAH : Prior : string scalar
              ;        Prior_arg : int, float, double, or string array(Narg)
              ;        Initialization_range : int, float, or double array(2)
              ;   The fraction of the total grain mass corresponding to
              ;   PAHs containing less than 1000 carbon atoms (PAH index).
              ;   Allowed parameter range: 4.7x10^(-3) to 4.58x10^(-2).
              QPAH:    {Prior: 'uniform'                            ,$
                        Prior_arg: [4.7d-3, 4.58d-2]                ,$
                        Initialization_range: [4.7d-3, 4.58d-2]      $
                        }                                           ,$

              ; LTIR : Prior : string scalar
              ;        Prior_arg : int, float, double, or string array(Narg)
              ;        Initialization_range : int, float, or double array(2)
              ;   The total integrated IR luminosity in Lsun.
              ; NOTE - If ENERGY_BALANCE is set, then this parameter is ignored
              ;        and determined instead by the absorbed the stellar (and,
              ;        if set, AGN) emission.
              ;   Allowed parameter range: 0 to Infinity.
              LTIR:    {Prior: 'uniform'                            ,$
                        Prior_arg: [0.0d, 1.d11]                    ,$
                        Initialization_range: [1.d7, 1.d11]          $
                        }                                           ,$




;===================================    X-RAY EMISSION    ====================================================
; X-ray emission module options

          ; XRAY_EMISSION : flag (0 or 1)
          ;   If set, an X-ray emission model will be used. This always includes
          ;   stellar X-ray emission, but can optionally include AGN X-ray
          ;   emission (see below).
          ; NOTE - If not set, all following X-ray emission model parameters will
          ;        be ignored.
          XRAY_EMISSION: 0                            ,$

          ; XRAY_UNIT : string scalar
          ;   The type of X-ray data to use for fitting.
          ; NOTE - If set to 'FLUX', the XRAY_UNC option below is ignored;
          ;        uncertainties on the X-ray flux must always be provided in the input catalog.
          ; Current options: 'COUNTS' and 'FLUX'
          XRAY_UNIT: 'COUNTS'                         ,$

          ; XRAY_UNC : string scalar
          ;   The uncertainties to assume for the X-ray counts.
          ;   Current options: 'SQRT', 'GEHRELS', and 'USER'
          XRAY_UNC: 'GEHRELS'                         ,$

          ; XRAY_ABS_MODEL : string scalar
          ;   The X-ray absorption model to apply to the X-ray emission.
          ;   Current options: 'TBABS-WILM', 'TBABS-ANGR', and 'ATTEN'
          XRAY_ABS_MODEL: 'TBABS-WILM'                ,$

          ; NH : Prior : string scalar
          ;      Prior_arg : int, float, double, or string array(Narg)
          ;      Initialization_range : int, float, or double array(2)
          ;   The intrinsic HI column density along the line of sight in 1e20 cm-2.
          ;   Allowed parameter range: 10^(-4) to 10^5.
          NH:          {Prior: 'uniform'                            ,$
                        Prior_arg: [1.d-4, 1.d5]                    ,$
                        Initialization_range: [1.d-1, 1.d2]          $
                        }                                           ,$

          ; XRAY_AGN_MODEL : string scalar
          ;   The AGN X-ray emission model to use. For no AGN X-ray emission model,
          ;   set to 'NONE'.
          ;   Current options: 'QSOSED', 'PLAW', and 'NONE'
          ; NOTE - If no X-ray AGN emission model is chosen, all AGN X-ray emission
          ;        model parameters will be ignored.
          XRAY_AGN_MODEL: 'QSOSED'                    ,$


          ;================= QSOSED PARAMETERS =========================================
              ; AGN_MASS : Prior : string scalar
              ;            Prior_arg : int, float, double, or string array(Narg)
              ;            Initialization_range : int, float, or double array(2)
              ;   The SMBH mass in Msun.
              ;   Allowed parameter range: 10^5 to 10^10.
              AGN_MASS:    {Prior: 'uniform'                        ,$
                            Prior_arg: [1.d5, 1.d10]                ,$
                            Initialization_range: [1.d6, 1.d9]       $
                            }                                       ,$

              ; AGN_LOGMDOT : Prior : string scalar
              ;               Prior_arg : int, float, double, or string array(Narg)
              ;               Initialization_range : int, float, or double array(2)
              ;   The log10 of the SMBH accretion rate normalized by the Eddington rate.
              ;   Allowed parameter range: -1.5 to 0.3.
              AGN_LOGMDOT: {Prior: 'uniform'                        ,$
                            Prior_arg: [-1.5d, 0.3d]                ,$
                            Initialization_range: [-1.d, 0.0d]       $
                            }                                       ,$




;====================================    AGN EMISSION    =====================================================
; AGN emission module options

          ; AGN_MODEL : string scalar
          ;   The UV-to-IR AGN emission model to use. For no UV-to-IR AGN
          ;   emission model, set to 'NONE'.
          ;   Current options: 'SKIRTOR' and 'NONE'
          ; NOTE - If no AGN emission model is chosen, all AGN emission
          ;        model parameters will be ignored.
          AGN_MODEL: 'NONE'                           ,$


          ;================= SKIRTOR PARAMETERS ==================================
              ; LOG_L_AGN : Prior : string scalar
              ;             Prior_arg : int, float, double, or string array(Narg)
              ;             Initialization_range : int, float, or double array(2)
              ;   The total integrated luminosity of AGN model in log10(Lsun).
              ; NOTE - This parameter is ignored if fitting an QSOSED X-ray AGN model.
              ;   Allowed parameter range: 0 to 20.
              LOG_L_AGN: {Prior: 'uniform'                          ,$
                          Prior_arg: [0.0d, 20.d0]                  ,$
                          Initialization_range: [6.d, 12.d]          $
                          }                                         ,$

              ; TAU97 : Prior : string scalar
              ;         Prior_arg : int, float, double, or string array(Narg)
              ;         Initialization_range : int, float, or double array(2)
              ;   The edge-on optical depth of AGN dust torus at 9.7 um.
              ;   Allowed parameter range: 3 to 11.
              TAU97:     {Prior: 'fixed'                            ,$
                          Prior_arg: [7.d0]                         ,$
                          Initialization_range: [3.0d, 11.0d]        $
                          }                                         ,$

              ; AGN_COSI : Prior : string scalar
              ;            Prior_arg : int, float, double, or string array(Narg)
              ;            Initialization_range : int, float, or double array(2)
              ;   The inclination of the AGN disk in terms of cos(i).
              ;   Allowed parameter range: 0 to 1.
              AGN_COSI:  {Prior: 'uniform'                          ,$
                          Prior_arg: [0.0d, 1.0d]                   ,$
                          Initialization_range: [0.0d, 1.0d]         $
                          }                                         ,$




;=================================    FITTING ALGORITHM    ===================================================
; Fitting algorithm module options

          ; METHOD : string scalar
          ;   The fitting algorithm used to fit the SED(s).
          ;   Current options: 'MCMC-ADAPTIVE', 'MCMC-AFFINE', and 'MPFIT'
          ; NOTE - Only configurations of the chosen fitting algorithm should be
          ;        updated from the default. Changes to the configurations
          ;        of the other fitting algorithms will be ignored.
          METHOD: 'MCMC-AFFINE'                       ,$


          ;===========  MCMC ===================================================
              ; NTRIALS : int, float, or double scalar
              ;   The number of MCMC trials to run.
              NTRIALS: 2e4                            ,$

              ; NPARALLEL : int, float, or double scalar
              ;   The number of parallel walkers/chains to run for each SED.
              ; NOTE - If using the affine-invariant algorithm, NPARALLEL must
              ;        be greater than the number of free parameters plus one
              ;        (ideally at least twice the number of free parameters)
              ;        for optimal sampling.
              NPARALLEL: 75                           ,$

              ; C_STEP : int, float, or double scalar
              ;   When calculating the autocorrelation time (tau) of the MCMC chain, this value
              ;   defines how many trials of the chain are used to calculate tau, where
              ;   we integrate tau to the smallest index ``m`` such that ``m > C_step * tau``.
              C_STEP: 5                               ,$

              ; TOLERANCE : int, float, or double scalar
              ;   When calculating the autocorrelation time (tau) of the MCMC chain, this value
              ;   defines how many taus the length of the chain should be for us to believe
              ;   the estimated value of tau.
              TOLERANCE: 50                           ,$

              ;======== ADAPTIVE ===============================================
                  ; BETA_EXPONENT : float or double scalar
                  ;   The factor controlling how fast the adaptiveness of the
                  ;   algorithm vanishes. Larger values stop the adaptiveness
                  ;   in fewer trials.
                  BETA_EXPONENT: 0.35                 ,$

              ;======== AFFINE =================================================
                  ; AFFINE_A : int, float, or double scalar
                  ;   The move scaling constant defining the maximum and
                  ;   minimum step size of the affine-invariant stretch move.
                  AFFINE_A: 2.0d                      ,$


          ;=========== MPFIT ===================================================
              ; NSOLVERS : int, float, or double scalar
              ;   The number of times to solve for the best fit SED using different
              ;   starting locations in parameters space.
              NSOLVERS: 100                           ,$

              ; FTOL : float or double scalar
              ;   The relative error desired in the sum of squares. Termination
              ;   of the MPFIT algorithm occurs when both the actual and predicted
              ;   relative reductions in the sum of squares are at most FTOL.
              FTOL: 1.d-10                            ,$

              ; GTOL : float or double scalar
              ;   The orthogonality desired between the function vector and the
              ;   columns of the Jacobian matrix. Termination of the MPFIT algorithm
              ;   occurs when the cosine of the angle between function vector and any
              ;   column of the Jacobian matrix is at most GTOL in absolute value.
              GTOL: 1.d-10                            ,$

              ; XTOL : float or double scalar
              ;   The relative error desired in the approximate solution. Termination
              ;   of the MPFIT algorithm occurs when the relative error between two
              ;   consecutive iterates is at most XTOL.
              XTOL: 1.d-10                            ,$

              ; MAXITER : int, float, or double scalar
              ;   The maximum number of MPFIT iterations to perform.
              MAXITER: 200                            ,$



;======================================    POST-PROCESSING   =========================================================
; Post-processing options

          ; KEEP_INTERMEDIATE_OUTPUT : flag (0 or 1)
          ;   If set, the intermediate sav files produced by the fitting algorithm
          ;   will not be deleted.
          ; NOTE - This is useful if needing to inspect the original fits before
          ;        post-processing.
          KEEP_INTERMEDIATE_OUTPUT: 0                 ,$


          ;================    MCMC Post-processing    ===================================
          ; NOTE - The following post-processing options only apply to the MCMC fitting algorithms.
          ;        If solely using the MPFIT algorithm, then a best fit model is determined for
          ;        each SED, and these options are not needed and will be ignored.

              ; BURN_IN : int, float, or double scalar
              ;   The number of initial MCMC trials to truncate as the burn-in phase. If set to 0,
              ;   then the number will be chosen automatically from the autocorrelation time.
              ; NOTE - We highly recommend specifying a value rather than using the automatic
              ;        calculation when using the MCMC-ADAPTIVE method as chains can vary widely
              ;        in the number of autocorrelation times needed for burn-in.
              BURN_IN: 0                              ,$

              ; THIN_FACTOR : int, float, or double scalar
              ;   The factor to thin the MCMC chain after removing the burn-in trials. For example,
              ;   a factor of 3 will include only every 3rd iteration in the final chain. A value of
              ;   1 means no thinning. If set to 0, then the number will be chosen automatically
              ;   from the autocorrelation time.
              ; NOTE - We recommend specifying a THIN_FACTOR of 1 when using the MCMC-ADAPTIVE method
              ;        as the final distribution will be from a single chain and each element in the
              ;        chain will be minimally correlated.
              THIN_FACTOR: 0                          ,$

              ; FINAL_CHAIN_LENGTH : int, float, or double scalar
              ;   The number of MCMC trials to include for the final distributions as taken from
              ;   the truncated and thinned chain.
              FINAL_CHAIN_LENGTH: 1000                ,$

              ; HIGH_RES_MODEL_FRACTION : int, float, or double scalar
              ;   The fraction of trials from FINAL_CHAIN_LENGTH sorted by quality of fit,
              ;   from which to generate high resolution models. If set to 0, then only the
              ;   best fit high resolution model will be generated.
              ; NOTE - Including more than the best fit high resolution model can drastically
              ;        increase file size. Be careful when increasing this value above 0. Doing
              ;        so will increase the post-processed file size by
              ;        FINAL_CHAIN_LENGTH * HIGH_RES_MODEL_FRACTION * 1e3 * 8 bytes per SED per
              ;        model component.
              HIGH_RES_MODEL_FRACTION: 0              ,$

              ;======== AFFINE =================================================
                  ; AFFINE_STRANDED_DEVIATION : int, float, or double scalar
                  ;  The number of standard deviations a walker must be below the median
                  ;  acceptance fraction of the ensemble to be considered a stranded walker.
                  AFFINE_STRANDED_DEVIATION: 2.0d      $

          }

 return, config

end
