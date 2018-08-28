;pro startup_RTE

;!prompt = "Chewbacca's IDL ----> "
;!quiet=1

;!EXCEPT=0
!P.title='!6'
!P.thick=2.
!x.thick=2
!y.thick=2
!P.charthick=1.5
!P.charsize=2

cv = { CONSTANTS_CGS, $
       hplanck       : 6.62606957D-27       ,$ ;Planck's constant [erg s]
       hbar          : 1.05457172d-27       ,$ ;hbar = !cv.hplanck/2.d0/!dpi [erg s]
       stefan        : 5.670373D-05         ,$ ;Stefan-Boltzmann constant [erg/cm^2/K^4]
       stephan       : 5.670373D-05         ,$ ;Stefan-Boltzmann constant [erg/cm^2/K^4] Eli's spelling
       kboltz        : 1.3806488D-16        ,$ ;Boltzmann constant [erg/K]
       wtmax         : 2.8973D+03           ,$ ;wavelength*max(T) of BB spectrum [um K]
       clight        : 2.99792458D+10       ,$ ;vacuum speed of light [cm/sec]
       eV            : 1.602176565D-12      ,$ ;1  eV [erg]
       keV           : 1.602176565D-09      ,$ ;1 keV [erg]
       MeV           : 1.602176565D-06      ,$ ;1 MeV [erg]
       GeV           : 1.602176565D-03      ,$ ;1 GeV [erg]
       TeV           : 1.602176565D-00      ,$ ;1 TeV [erg]
       Rydberg       : 2.1798741D-11        ,$ ;Rydberg in ergs = m_e*e^4/(2*h_bar^2)
       G             : 6.67384D-08          ,$ ;Gravitation constant, dynes/gm/cm^2 [cm^3 g^-1 s^-2]
       mH            : 1.673532666D-24      ,$ ;mass of Hydrogen, g.
       mel           : 9.10938291D-28       ,$ ;electron mass [g]
;      statC         : 3.33564095d-10       ,$ ;statC = !cv.clight/1.d1
       elcharge      : 4.8032045057d-10     ,$ ;elementary charge [statC] = 1.602176565d-19 Coulombs = 1.602176565d-19/(10.d0/!cv.clight) statC
;      elcharge      : 4.8032055659d-10     ,$ ;elementary charge [statC] = sqrt(sqrt(2.d0*(!cv.hplanck/2d0/!dpi)^2*!cv.Rydberg/!cv.mel))
       Watt          : 1.d7                 ,$ ;1 W in erg s-1
       Jy            : 1.d-23               ,$ ;1 Jy in erg s-1 cm-2 Hz-1
       mJy           : 1.d-26               ,$ ;1 Jy in erg s-1 cm-2 Hz-1
       Lsun          : 3.8427D+33           ,$ ;solar luminosity, ergs/sec
       Msun          : 1.9891D+33           ,$ ;solar mass, g.
       Rsun          : 6.955D+10            ,$ ;solar radius, cm.
       a0            : 5.29177D-9           ,$ ;Bohr radius, cm.
       txsec         : 6.6524586D-25        ,$ ;Thomson cross section, cm^2.
       pc            : 3.08567758149D+18    ,$ ;parsec [cm] = !cv.AU/tan(!cv.arcsec)
       year          : 3.15576D+07          ,$ ;year in seconds
       AU            : 1.49597870700D+13    ,$ ;sun to earth distance [cm].
       deg           : 1.7453292519943D-02  ,$ ; 1 degree [rad] = !dpi/1.8d2
       min           : 2.9088821866572D-04  ,$ ; 1 minute [rad] = !dpi/1.8d2/6.0d1
       arcsec        : 4.8481368110954D-06  ,$ ; 1 arcsec [rad] = !dpi/1.8d2/3.6d3
       sec           : 4.8481368110954D-06  ,$ ; 1 arcsec [rad]
       deg2          : 3.0461741978671D-04  ,$ ; 1 square degree [sr]
       min2          : 8.4615949940752D-08  ,$ ; 1 square arcmin [sr
       sec2          : 2.3504430539098D-11  ,$ ; 1 square arcsec [sr]
       sr_to_arcsec2 : (1.8d2*3.6d3/!dpi)^2 ,$ ; 1 sr in square arcsec
       sr_to_arcmin2 : (1.8d2*6.0d1/!dpi)^2 ,$ ; 1 sr in square arcmin
       sr_to_deg2    : (1.8d2/!dpi)^2       ,$ ; 1 sr in square degree
       beam_to_FWHM2 : (!dpi^3/4.d0/alog(2.d0)/(1.8d2*3.6d3)^2) $ ; beam to FWHM/1", in sr
                                                                  ; i.e., BEAM/sr = !cv.beam_to_FWHM2*FWHM^2
     }

defsysv, "!CV", cv, 1 ;define as read-only sys-var.

;device,decompose=0    ; allows color table changes at 24-bit color depth

;get_colors

;end
