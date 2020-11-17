             Dust Emission Spectra from Draine & Li (2007) 

This directory contains results of theoretical modeling of infrared
emission from dust grains heated by starlight, as described in Draine
& Li (2007) -- please consult that paper for details.

In brief, the dust consists of a mixture of amorphous silicate grains
and carbonaceous grains.  The carbonaceous grains have the properties
of PAH molecules and ions when the effective radius a < 5.0nm, the
properties of graphite spheres when a >> 10 nm, and optical properties
intermediate between the properties of PAH particles and graphite
particles for 5 < a < 10 nm, as described by Draine & Li (2007).

The ionization fraction x_ion(a) of the PAH particles is assumed to be
the average for the diffuse ISM [see Figure 8 of Draine & Li (2007)].

Temperature distribution functions are calculated for all particles
small enough for quantized heating to be important; large grains are
treated as having steady-state temperatures determined by starlight
heating = radiative cooling.

Spectra are given for dust heated by starlight with the spectrum of
Mathis, Mezger, & Panagia (1983), scaled by a factor U.  U=1
corresponds to the starlight intensity estimate for the local ISM.
Band-averaged emissivities are also given for photometric bands used
by IRAS, COBE-DIRBE, Spitzer Space Telescope, Akari, and Herschel.

Spectra are given for 11 dust models:

  jm     name      q_PAH
  --   --------    ------
   1   MW3.1_00    0.47 %
   2   MW3.1_10    1.12 %
   3   MW3.1_20    1.77 %
   4   MW3.1_30    2.50 %
   5   MW3.1_40    3.19 %
   6   MW3.1_50    3.90 %
   7   MW3.1_60    4.58 %
   8   LMC2_00     0.75 %
   9   LMC2_05     1.49 %
  10   LMC2_10     2.37 %
  11   smc         0.10 %

q_PAH is the fraction of the total dust mass that is contributed by
PAH particles containing < 1000 C atoms.  All dust models use the dust
size distributions from Weingartner & Draine (2001), except for
adjustment of the parameters (a_01,a_02,sigma_1,sigma_2)
characterizing the PAH size distribution (see Draine & Li 2007, Table
2).

For each emission model, the dust is assumed to be exposed to a range
of radiation intensities, with the mass dM of dust exposed to
starlight intensities in [U,U+dU] being given by a power-law
distribution

   dM = const * U^{-2} * dU   for umin < U < umax

Results are available here for the following values of umin:

   umin= 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.70, 0.80, 1.00, 1.20,
   1.50, 2.00, 2.50, 3.00, 4.00, 5.00, 7.00, 8.00, 10.0, 12.0, 15.0,
   20.0, 25.0

and for the following values of umax:

   umax= 1e3, 1e4, 1e5, 1e6, 1e7

The model spectra are given in files

Uumin/Uumin_umax_modl.txt [e.g., U1.00/U1.00_1e6_MW3.1_60.txt ]

We also consider dust exposed to a single radiation intensity.  This
is simply the case where umax = umin.  These spectra are given in files

Uumin/Uumin_umin_modl.txt [e.g., U1.00/U1.00_1.00_MW3.1_60.txt ]

The files are ascii, and should be self-explanatory.  

Draine & Li (2007) propose fitting emission from galaxies, or large
regions within a galaxy, by a linear sum of emission from dust heated
by a single radiation intensity U=umin plus emission from dust heated
by a distribution of starlight intensities ranging from umin to umax,
with umax=1e6 appearing to often work well for star-forming galaxies.
The emission spectrum, expressed as emission per H nucleon, is then
simply

   j_nu = (1-gamma)*j_nu[umin,umin] + gamma*j_nu[umin,umax]

where gamma is the fraction of the total dust mass that is heated by
the distribution of starlight intensities [the remaining fraction
(1-gamma) is heated by starlight with U=umin], and j_nu[umin,umax] is
the emissivity from the file U$umin_umax_dustmod where dustmod is one
of {MW3.1_00, MW3.1_10, MW3.1_20, MW3.1_30, MW3.1_40, MW3.1_50,
MW3.1_60, LMC2_00, LMC2_05, LMC2_10, smc}.

Portions of the file U1.00/U1.00_1.00_MW3.1_60.txt are appended below.

References:

Draine, B.T., & Li, A. (2001), "Infrared Emission from Interstellar
   Dust, I.  Stochastic Heating of Small Grains", ApJ, 551, 807-824

Draine, B.T., & Li, A. (2007), "Infrared Emission from Interstellar
   Dust, IV. The Silicate-Graphite-PAH Model in the Post-Spitzer Era",
   ApJ, 657, 810-837

Li, A, & Draine, B.T. (2001), "Infrared Emission from Interstellar
   Dust, II.  The Diffuse Interstellar Medium", ApJ, 554, 778-802

Weingartner, J.C., & Draine, B.T. (2001), "Dust Grain Size
   Distributions and Extinction in the Milky Way, LMC, and SMC", ApJ,
   548, 296-309


         Appendix: The File  U1.00/U1.00_1.00_MW3.1_60.txt
         =================================================

The first part of the file gives some model properties, including the
spectrum convolved with 23 photometric bands used by COBE-DIRBE, IRAS,
Spitzer Space Telescope, Akari, and Herschel

%%%%%%%%%%%%%%%%%%%%% first 53 lines follow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dust emission from Draine & Li (2007), ApJ 657, 810-837
==========================================================
         7 = grain model
 4.000E-08 4.000E-01 4.500E-05 = a_01 (cm) , sigma_1 , b_C1/0.92
 2.000E-07 5.500E-01 1.500E-05 = a_02 (cm) , sigma_2 , b_C2/0.92
 1.000E+00 1.000E+00 1.000E+00 = Umin , Umax, beta (dN/dU propto U^{-beta})
 1.000E+00 = <U> 
     mmp83 = radfield
 4.558E-24 = power/H (erg s-1)
 band   nu*P_nu(band) j_nu(band)
 (um)   (erg s-1 H-1)(Jy cm2 sr-1 H-1)
1.270E+00  1.841E-26  6.202E-19 COBE-DIRBE
2.220E+00  2.349E-26  1.383E-18 COBE-DIRBE
3.530E+00  1.526E-25  1.429E-17 COBE-DIRBE
4.880E+00  6.843E-26  8.858E-18 COBE-DIRBE
1.229E+01  5.872E-25  1.914E-16 COBE-DIRBE
2.079E+01  3.442E-25  1.898E-16 COBE-DIRBE
5.599E+01  7.648E-25  1.136E-15 COBE-DIRBE
9.770E+01  2.311E-24  5.988E-15 COBE-DIRBE
1.479E+02  2.610E-24  1.024E-14 COBE-DIRBE
2.479E+02  1.127E-24  7.407E-15 COBE-DIRBE
1.200E+01  6.488E-25  2.065E-16 IRAS
2.500E+01  3.008E-25  1.995E-16 IRAS
6.000E+01  9.037E-25  1.438E-15 IRAS
1.000E+02  2.380E-24  6.312E-15 IRAS
3.550E+00  1.775E-25  1.671E-17 SST-IRAC
4.493E+00  4.980E-26  5.935E-18 SST-IRAC
5.731E+00  3.961E-25  6.021E-17 SST-IRAC
7.872E+00  9.580E-25  2.000E-16 SST-IRAC
2.368E+01  2.646E-25  1.662E-16 SST-MIPS
7.142E+01  1.243E-24  2.355E-15 SST-MIPS
1.559E+02  2.443E-24  1.010E-14 SST-MIPS
1.600E+01  4.360E-25  1.850E-16 SST-IRSPU
2.200E+01  2.913E-25  1.700E-16 SST-IRSPU
7.500E+01  1.045E-24  2.078E-15 Herschel-PACS
1.100E+02  1.757E-24  5.127E-15 Herschel-PACS
1.700E+02  1.866E-24  8.416E-15 Herschel-PACS
2.500E+02  9.866E-25  6.542E-15 Herschel-SPIRE
3.600E+02  3.264E-25  3.116E-15 Herschel-SPIRE
5.200E+02  8.803E-26  1.214E-15 Herschel-SPIRE
2.400E+00  2.633E-26  1.676E-18 Akari IRC NIR N2
3.200E+00  1.174E-25  9.966E-18 Akari IRC NIR N3
4.100E+00  5.402E-26  5.875E-18 Akari IRC NIR N4
7.000E+00  1.088E-24  2.020E-16 Akari IRC MIR-S S7
9.000E+00  7.666E-25  1.830E-16 Akari IRC MIR-S S9W
1.100E+01  5.932E-25  1.731E-16 Akari IRC MIR-S S11
1.500E+01  4.522E-25  1.799E-16 Akari IRC MIR-L L15
1.800E+01  3.594E-25  1.716E-16 Akari IRC MIR-L L18W
2.400E+01  2.800E-25  1.782E-16 Akari IRC MIR-L L24
6.500E+01  1.127E-24  1.942E-15 Akari FIS N60
8.000E+01  1.945E-24  4.127E-15 Akari FIS WIDE-S
1.400E+02  2.632E-24  9.775E-15 Akari FIS WIDE-L
1.600E+02  2.431E-24  1.032E-14 Akari FIS N160

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Below this are some blank lines, followed by the emission spectrum
between 1 cm and 1 micron (lines 62 - 1062, in frequency steps of
10**0.004).  The quantities tabulated are nu* dP/d\nu (where \nu is
frequency, and P = power radiated per H nucleon) and j_\nu = power
radiated per H nucleon per unit frequency per steradian.  Obviously,
j_\nu = (1/4*\pi)*(dP/d\nu)

%%%%%%%%%%%%%%%%%%%%%%%%%%% lines 60 - 70 follow %%%%%%%%%%%%%%%%%%%%%%
lambda    nu*dP/dnu     j_nu
 (um)   (erg s-1 H-1)(Jy cm2 sr-1 H-1)
1.000E+04  5.303E-31  1.407E-19
9.908E+03  5.497E-31  1.445E-19
9.817E+03  5.697E-31  1.484E-19
9.727E+03  5.905E-31  1.524E-19
9.638E+03  6.120E-31  1.565E-19
9.550E+03  6.344E-31  1.607E-19
9.462E+03  6.575E-31  1.650E-19
9.376E+03  6.815E-31  1.695E-19
9.290E+03  7.063E-31  1.741E-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
