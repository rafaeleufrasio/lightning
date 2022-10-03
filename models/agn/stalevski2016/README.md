#### SKIRTOR AGN Templates
This directory contains serialized (IDL savefile formatted) versions of the SKIRTOR modelset (Stalevski et al. 2012a; 2016; see https://sites.google.com/site/skirtorus/home).

The file naming convention is the same as the original release, described below. The files are sorted into folders `iXX/` where `XX` is the inclination of the models. The models are gridded over 6 parameters plus inclination. The grid values are listed in the table below.

|Parameter |Values                                  |Num. Values|
|----------|---------------------------------------|-----------:|
| t        | [3,5,7,9,11]                           |          5|
| p        | [0.0, 0.5, 1.0, 1.5]                   |          4|
| q        | [0.0, 0.5, 1.0, 1.5]                   |          4|
| oa       | [10, 20, 30, 40, 50, 60, 70, 80]       |          8|
| R        | [10, 20, 30]                           |          3|
| Mcl      | [0.97]                                 |          1|
| i        | [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]|         10|
|**Total** |                    ...                 |      19200|

The parameters are described in detail in the original documentation, part of which is reproduced below.

When loaded with the `RESTORE` command each file defines an structure named `SKIRTOR_model` which is arranged as follows:

|Attribute            |Type   |Dimension  |
|---------------------|-------|-----------|
|`MODEL_NAME`         |STRING |Scalar     |
|`INCLINATION`        |STRING |Scalar     |
|`WAVE_REST`          |FLOAT  |Array[132] |
|`NU_REST`            |DOUBLE |Array[132] |
|`LNU_TOTAL`          |DOUBLE |Array[132] |
|`LNU_DISK_DIRECT`    |DOUBLE |Array[132] |
|`LNU_DISK_SCATTERED` |DOUBLE |Array[132] |
|`LNU_DUST_TOTAL`     |DOUBLE |Array[132] |
|`LNU_DUST_SCATTERED` |DOUBLE |Array[132] |
|`LNU_TRANSPARENT`    |DOUBLE |Array[132] |
|`DUST_MASS_TOTAL`    |FLOAT  |Array[1]   |


`MODEL_NAME` refers to a string like `t5_p1_q0_oa50_R20_Mcl0.97` encoding the physical model parameters; `INCLINATION` refers to a string containing the inclination (e.g. `'30'`). The attribute `WAVE_REST` is rest-frame wavelength in microns; `NU_REST` is rest-frame frequency in Hz. We have converted the original models to luminosity density L_nu: the `LNU_*` attributes are rest-frame luminosity density, in L_sol / Hz. The order of the `LNU_*` attrubutes above is the same as in the flux columns in the original SKIRTOR text files from the 2016 release; the relevant portion of the documentation for those files is reproduced below. The final attribute, `DUST_MASS_TOTAL`, contains the total dust mass of the model in M_sol.

##### Original SKIRTOR Documentation
SED files naming convention and format
File name example: `t5_p1_q0_oa50_R20_Mcl0.97_i30_sed.dat`

- t: tau9.7, average edge-on optical depth at 9.7 micron; the actual one along the line of sight may vary depending on
the clumps distribution.
- p: power-law exponent that sets radial gradient of dust density
- q: index that sets dust density gradient with polar angle
- oa: angle measured between the equatorial plan and edge of the torus. Half-opening angle of the dust-free cone is 90-oa.
- R: ratio of outer to inner radius, R_out/R_in
- Mcl: fraction of total dust mass inside clumps. 0.97 means 97% of total mass is inside the clumps and 3% in the
interclump dust.
- i: inclination, i.e. viewing angle, i.e. position of the instrument w.r.t. the AGN axis. i=0: face-on, type 1 view; i=90: edge-on, type 2 view.

Flux in the SED files is calculated for a source at a distance of 10 Mpc and luminosity L_AGN=10^11 L_sol, where L_sol = 3.839e26 W (= 3.839e33 erg/s) (without solar neutrino radiation).

It is trivial to rescale the flux so that model SEDs can be applied to any given source luminosity and distance, just keeping in mind the scaling of dust emission, size and mass (see Ivezic & Elitzur 1997; equation 14 in Fritz et al. 2006; Honig & Kishimoto 2010).

Included in the SED library is the file with total dust masses for all the models; again, keep in mind these values (as well as size of the torus) need to be rescaled when applied to the sources with different L_AGN (see references above).

Primary source of radiation, the accretion disk, is represented by a point-like central source with anisotropic emission pattern, as suggested by Netzer (1987).

The adopted R_in is 0.5 pc (note that in Stalevski et al. 2016 the different value is reported; that is a typo). However, the actual R_in value depends on the polar angle, following anisotropic emission pattern of the accretion disk.

Header of the SED files contains the required description. A couple of clarifications:

###### column 3: direct stellar flux:
###### column 4: scattered stellar flux
"stellar" is a remnant from the original version of the code which provided only stars as primary source of radiation. In AGN models presented here, it stands for accretion disk emission. "direct" means what reached the instrument at the given viewing angle, after being absorbed or scattered away by the dust. If there is no dust along the line of sight for the given viewing angle it will be identical to the "transparent" flux.

###### column 7: transparent flux:
"transparent" stands for the original primary source flux (accretion disk), unaffected by the dust. "transparent" flux is the same for all the models, but different for each inclination, as anisotropic emission for primary source is assumed.