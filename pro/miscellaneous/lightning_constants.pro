pro lightning_constants
;+
; Name
; ----
;   LIGHTNING_CONSTANTS
;
; Purpose
; -------
;   Generates a system variable that contains all of the constants needed
;   for unit conversions in Lightning. All constants are in CGS units.
;
; Calling Sequence
; ----------------
;   ::
;
;       lightning_constants
;
; Output
; ------
;   The created system variable, ``!lightning_cgs``.
;
; Modification History
; --------------------
;   - 2016/05/01: Created (Rafael T. Eufrasio)
;   - 2022/03/22: Changed file name and made into a procedure (Keith Doore)
;   - 2022/03/22: Added documentation (Keith Doore)
;   - 2022/04/18: Removed structure name to prevent multiple call issues (Keith Doore)
;   - 2022/05/09: Changed system variable name from ``!cv`` to ``!lightning_cgs`` (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

 cgs = {hplanck       : 6.62606957d-27       ,$ ;Planck's constant [erg s]
        hbar          : 1.05457172d-27       ,$ ;hbar = !cgs.hplanck/2.d0/!dpi [erg s]
        stefan        : 5.670373d-05         ,$ ;Stefan-Boltzmann constant [erg/cm^2/K^4]
        kboltz        : 1.3806488d-16        ,$ ;Boltzmann constant [erg/K]
        wtmax         : 2.8973d-01           ,$ ;Wien's displacement constant of BB spectrum [cm K]
        clight        : 2.99792458d+10       ,$ ;vacuum speed of light [cm/s]
        eV            : 1.602176565d-12      ,$ ;1  eV [erg]
        keV           : 1.602176565d-09      ,$ ;1 keV [erg]
        MeV           : 1.602176565d-06      ,$ ;1 MeV [erg]
        GeV           : 1.602176565d-03      ,$ ;1 GeV [erg]
        TeV           : 1.602176565d-00      ,$ ;1 TeV [erg]
        Rydberg       : 2.1798741d-11        ,$ ;Rydberg unit of energy = m_e*e^4/(2*h_bar^2) [ergs]
        G             : 6.67384d-08          ,$ ;Gravitation constant, dynes g^-2 cm^2 [cm^3 g^-1 s^-2]
        mH            : 1.673532666d-24      ,$ ;mass of Hydrogen [g]
        mel           : 9.10938291d-28       ,$ ;electron mass [g]
 ;      statC         : 3.33564095d-10       ,$ ;statC = !cgs.clight/1.d1
        elcharge      : 4.8032045057d-10     ,$ ;elementary charge [statC] = 1.602176565d-19 Coulombs = 1.602176565d-19/(10.d0/!cgs.clight) statC
                                                ;                          = sqrt(sqrt(2.d0*(!cgs.hplanck/2d0/!dpi)^2*!cgs.Rydberg/!cgs.mel))
        Watt          : 1.d7                 ,$ ;1 W [erg/s]
        Jy            : 1.d-23               ,$ ;1 Jy [erg s-1 cm-2 Hz-1]
        mJy           : 1.d-26               ,$ ;1 mJy [erg s-1 cm-2 Hz-1]
        Lsun          : 3.8427d+33           ,$ ;solar luminosity [ergs/sec]
        Msun          : 1.9891d+33           ,$ ;solar mass [g]
        Rsun          : 6.955d+10            ,$ ;solar radius [cm]
        a0            : 5.29177d-9           ,$ ;Bohr radius [cm]
        txsec         : 6.6524586d-25        ,$ ;Thomson cross section [cm^2]
        pc            : 3.08567758149d+18    ,$ ;parsec [cm] = !cgs.AU/tan(!cgs.arcsec)
        year          : 3.15576d+07          ,$ ;year [s]
        AU            : 1.49597870700d+13    ,$ ;sun to earth distance [cm]
        deg           : 1.7453292519943d-02  ,$ ; 1 degree [rad] = !dpi/1.8d2
        min           : 2.9088821866572d-04  ,$ ; 1 minute [rad] = !dpi/1.8d2/6.0d1
        arcsec        : 4.8481368110954d-06  ,$ ; 1 arcsec [rad] = !dpi/1.8d2/3.6d3
        sec           : 4.8481368110954d-06  ,$ ; 1 arcsec [rad]
        deg2          : 3.0461741978671d-04  ,$ ; 1 square degree [sr]
        min2          : 8.4615949940752d-08  ,$ ; 1 square arcmin [sr]
        sec2          : 2.3504430539098d-11  ,$ ; 1 square arcsec [sr]
        sr_to_arcsec2 : (1.8d2*3.6d3/!dpi)^2 ,$ ; 1 sr [arcsec^2]
        sr_to_arcmin2 : (1.8d2*6.0d1/!dpi)^2 ,$ ; 1 sr [arcmin^2]
        sr_to_deg2    : (1.8d2/!dpi)^2       ,$ ; 1 sr [deg^2]
        beam_to_FWHM2 : (!dpi^3/4.d0/alog(2.d0)/(1.8d2*3.6d3)^2) $ ; beam to FWHM/1", in sr
                                                                   ; i.e., BEAM/sr = !cgs.beam_to_FWHM2*FWHM^2
        }

 defsysv, "!lightning_cgs", cgs, 1 ;define as read-only sys-var.

end