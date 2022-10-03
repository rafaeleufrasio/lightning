function doore21_atten, wave, tauB_f=tauB_f, F_clump=F_clump, cosi=cosi, $
                        rold0=rold0, b_to_d=b_to_d, tuffs_coeff=tuffs_coeff, $
                        error_check=error_check
;+
; Name
; ----
;   DOORE21_ATTEN
;
; Purpose
; -------
;   Generates attenuation values at the input wavelengths using the attenuation
;   curves describe in Doore et al. (2021) given the attenuation curve
;   parameters. The attenuation curves are based on the Tuffs et al. (2004)
;   attenuation curves as updated by Popescu et al. (2011), which include the
;   parameters of inclination, face-on optical depth in the B-band, the
;   clumpiness factor, the fraction of intrinsic flux density from the old
;   stellar components compared to the total intrinsic flux density, and the
;   bulge-to-disk ratio.
;
; Calling Sequence
; ----------------
;   ::
;
;       exp_neg_tau = doore21_atten(wave [, tauB_f = , F_clump = , cosi = , $
;                                   rold0 = , b_to_d = , tuffs_coeff = , /error_check])
;
; Input
; -----
;   ``wave`` : int, float, or double array(Nwave)
;       The wavelength at which to determine the attenuation :math:`[\mu \rm m]`.
;
; Optional Inputs
; ---------------
;   ``tauB_f`` : int, float, or double array(Nmodels)
;       The face-on optical depth in the B-band. (Default = ``1.d0``)
;   ``F_clump`` : int, float, or double array(Nmodels)
;       The clumpiness factor F. (Default = ``0.d0``)
;   ``cosi`` : int, float, or double array(Nmodels)
;       The inclination of the galactic disk in terms of cos(i). (Default = ``1.d0``)
;   ``rold0`` : int, float, or double array(Nmodels)
;       The fraction of intrinsic flux density from the old stellar components
;       compared to the total intrinsic flux density. (Default = ``0.d0``)
;   ``b_to_d`` : int, float, or double array(Nmodels)
;       The bulge-to-disk ratio. (Default = ``0.d0``)
;   ``tuffs_coeff`` : structure
;       A structure containing the polynomial coefficients needed for the 
;       Tuffs model.
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``exp_neg_tau`` : double array(Nwave, Nmodels)
;       The attenuation in terms of :math:`e^{-\tau}` (i.e., :math:`e^{-\tau} = 10^{-0.4 A_{\lambda}}`).
;
; Note
; ----
;   Defaults will only be set if the optional ``error_check`` input is set.
;
; References
; ----------
;   - `Tuffs, R. J., Popescu, C. C., Völk, H. J., Kylafis, N. D., & Dopita, M. A. 2004, A&A, 419, 821 <https://ui.adsabs.harvard.edu/abs/2004A%26A...419..821T/abstract>`_
;   - `Popescu, C. C., Tuffs, R. J., Dopita, M. A., et al. 2011, A&A, 527, A109 <https://ui.adsabs.harvard.edu/abs/2011A%26A...527A.109P/abstract>`_
;   - `Doore, K., Eufrasio, R. T., Lehmer, B. D., et al. 2021, ApJ, 923, 26 <https://ui.adsabs.harvard.edu/abs/2021ApJ...923...26D/abstract>`_
;
; Modification History
; --------------------
;   - 2019/09/22: Created (Keith Doore)
;   - 2021/04/13: Modified to have ``rdisk`` as intrinsic property (``rold0``) and ``rbulge`` as B/D (Keith Doore)
;   - 2022/03/15: Added proper error handling (Keith Doore)
;   - 2022/03/15: Renamed variables to standard format (Keith Doore)
;   - 2022/03/18: Updated Documentation (Keith Doore)
;   - 2022/03/31: Updated for efficiency, minimal arrays in memory (Keith Doore)
;   - 2022/04/07: Allowed for inputs to have degenerate dimensions (Keith Doore)
;   - 2022/04/07: Allowed for inputs to be scalars (Keith Doore)
;   - 2022/04/07: Allowed integer inputs (Keith Doore)
;   - 2022/04/07: Allowed for array of only one optional input to be given (Keith Doore)
;   - 2022/04/14: Allowed for ``b_to_d`` to be ``Inifinty`` and not cause ``NaNs`` in ``rbulge0`` (Keith Doore)
;   - 2022/04/15: Updated interpolation function for improved speed (Keith Doore)
;   - 2022/04/28: Removed for loop on 1-cosi polynomial for improved speed (Keith Doore)
;   - 2022/04/28: Removed unnecessary zeroing of arrays (Keith Doore)
;   - 2022/05/18: Turned ``lightning_dir`` into system variable call (Keith Doore)
;   - 2022/05/18: Added ``error_check`` keyword to do error handling (Keith Doore)
;   - 2022/05/18: Updated path to model files (Keith Doore)
;   - 2022/06/09: Made the coefficients an optional input for improved speed (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(wave) eq 0 then message, 'Variable is undefined: WAVE.'
   if size(wave, /type) lt 2 or size(wave,/type) gt 5 then message, 'WAVE must be of type int, float, or double.'
   if size(reform(wave), /n_dim) ne 1 then message, 'WAVE must be a scalar or 1-D array.'
   if min(wave) le 0 then message, 'WAVE must only contain positive values.'
  
   if n_elements(tauB_f) ne 0 then begin
     if size(tauB_f, /type) lt 2 or size(tauB_f,/type) gt 5 then $
              message, 'TAUB_F must be of type int, float, or double.'
     if size(reform(tauB_f), /n_dim) ne 1 then message, 'TAUB_F must be a scalar or 1-D array.'
     if min(tauB_f) lt 0 or max(tauB_f) gt 8 then message, 'TAUB_F must contain values between 0 and 8.'
     Nmodels = n_elements(taub_f)
   endif
  
   if n_elements(cosi) ne 0 then begin
     if size(cosi, /type) lt 2 or size(cosi,/type) gt 5 then $
              message, 'COSI must be of type int, float, or double.'
     if size(reform(cosi), /n_dim) ne 1 then message, 'COSI must be a scalar or 1-D array.'
     if min(cosi) lt 0 or max(cosi) gt 1 then message, 'COSI must contain values between 0 and 1.'
     if n_elements(Nmodels) eq 0 then Nmodels = n_elements(cosi)
   endif
  
   if n_elements(F_clump) ne 0 then begin
     if size(F_clump, /type) lt 2 or size(F_clump,/type) gt 5 then $
              message, 'F_CLUMP must be of type int, float, or double.'
     if size(reform(F_clump), /n_dim) ne 1 then message, 'F_CLUMP must be a scalar or 1-D array.'
     if min(F_clump) lt 0.d0 or max(F_clump) gt 0.61d0 then message, 'F_CLUMP must contain values between 0 and 0.61.'
     if n_elements(Nmodels) eq 0 then Nmodels = n_elements(F_clump)
   endif
  
   if n_elements(rold0) ne 0 then begin
     if size(rold0, /type) lt 2 or size(rold0,/type) gt 5 then $
              message, 'ROLD0 must be of type int, float, or double.'
     if size(reform(rold0), /n_dim) ne 1 then message, 'ROLD0 must be a scalar or 1-D array.'
     if total(rold0 eq 0 or rold0 eq 1) ne n_elements(rold0) then $
       message ,'ROLD0 must only contain values of 0 or 1.'
     if n_elements(Nmodels) eq 0 then Nmodels = n_elements(rold0)
   endif
  
   if n_elements(b_to_d) ne 0 then begin
     if size(b_to_d, /type) lt 2 or size(b_to_d,/type) gt 5 then $
              message, 'B_TO_D must be of type int, float, or double.'
     if size(reform(b_to_d), /n_dim) ne 1 then message, 'B_TO_D must be a scalar or 1-D array.'
     if min(b_to_d) lt 0 then message, 'B_TO_D must only contain non-negative values.'
     if n_elements(Nmodels) eq 0 then Nmodels = n_elements(b_to_d)
   endif
  
   if n_elements(Nmodels) ne 0 then begin
     if n_elements(tauB_f ) eq 0 then tauB_f  = replicate(1.d0, Nmodels)
     if n_elements(rold0  ) eq 0 then rold0   = replicate(0.d0, Nmodels)
     if n_elements(b_to_d ) eq 0 then b_to_d  = replicate(0.d0, Nmodels)
     if n_elements(cosi   ) eq 0 then cosi    = replicate(1.d0, Nmodels)
     if n_elements(F_clump) eq 0 then F_clump = replicate(0.d0, Nmodels)
   endif else begin
     Nmodels = 1
     tauB_f  = 1.d0
     rold0   = 0.d0
     b_to_d  = 0.d0
     cosi    = 1.d0
     F_clump = 0.d0
   endelse

   if n_elements(tauB_f) ne n_elements(rold0) or n_elements(tauB_f) ne n_elements(b_to_d)  or $
      n_elements(tauB_f) ne n_elements(cosi)  or n_elements(tauB_f) ne n_elements(F_clump) then $
       message,'TAUB_F, ROLD0, B_TO_D, COSI, and F_CLUMP must have the same size.'

   if n_elements(tuffs_coeff) ne 0 then begin
     if size(tuffs_coeff, /type) ne 8 then $
              message, 'TUFFS_COEFF must be of type structure.'
   endif  
 endif
 Nmodels = n_elements(taub_f) 


; Derive the appropriate attenuation curve
 ; coeff.dat file from ftp site (http://cdsarc.u-strasbg.fr/ftp/J/A+A/527/A109/coeff.dat) in Popescu et al (2011)
 ; Converted to convenient 3d matrix where [tau,wavebands,coeff] are the dimensions and values
 if keyword_set(tuffs_coeff) then begin
   abulge = tuffs_coeff.abulge
   adisk  = tuffs_coeff.adisk
   atdisk = tuffs_coeff.atdisk
 endif else restore, !lightning_dir+'models/attenuation/doore21/tuffs_coeff.sav'

 ; tau values used in Table 3 in Popescu et al (2011)
 tau = [0.1d, 0.3, 0.5, 1.0, 2.0, 4.0, 8.0]

 ; Include waveband at 50000 Angstroms where mlambda is 0 for interpolation smoothness
 ; Wavelengths of the bands used (Angstroms)
 wavebands = [912.d, 1350, 1500, 1650, 2000, 2200, 2500, 2800, 3650, 4430, 5640, 8090, 12590, 22000, 50000]

 ; Values from Table E.4 in Popescu et al (2011)
 Fcal = 0.35d0
 flamb_fcal_column = [0.427d, 0.484, 0.527, 0.545, 0.628, 0.676, 0.739, 0.794, 0.892, 0.932, 0.962, 0.984, 0.991, 0.994, 0.999]

 ; [value=(1−Fcal*fλ)]->[(1-value)/Fcal=fλ]
 ;fλ values
 flambda = (1-flamb_fcal_column)/Fcal

 ; Convert b_to_d ratio into rbulge0 and rdisk0
 ; rdisk0+rbulge0=rold0 and rdisk0*B/D=rbulge0
 rdisk0 = reform(rold0/(1+b_to_d))
 rbulge0 = reform(rold0-rdisk0)

 ; Add a tau at 0 for interpolation smoothness
 ; Set the diffuse attenuation to 0 for this case
 tau = [0.0, tau]
 Ntau = n_elements(tau)
 Nwave = n_elements(wave)
 Nwavebands = n_elements(wavebands)


 mlambda_d  = dblarr(ntau, nwavebands, Nmodels)
 mlambda_td = dblarr(ntau, nwavebands, Nmodels)
 mlambda_b  = dblarr(ntau, nwavebands, Nmodels)

 ; Calculate delta_mlambda for each component of galaxy
 ; Using Equation 6 from Tuffs et al (2004)
 ; First tau index is tau=0 and last wavelength index is wavelength=50000 Angstroms. Set both to 0 for mlambda
 one_minus_cosi_mat = rebin(reform((1-cosi), 1, 1, Nmodels), ntau, nwavebands, Nmodels)
 a_disk  = rebin(reform(adisk,  ntau-1, nwavebands-1, 1, 6), ntau-1, nwavebands-1, Nmodels, 6)
 mlambda_d[1:*,0:-2,*]= (one_minus_cosi_mat[1:*,0:-2,*])^0.d0*a_disk[*,*,*,0]+ (one_minus_cosi_mat[1:*,0:-2,*])^1.d0*a_disk[*,*,*,1]+$
                        (one_minus_cosi_mat[1:*,0:-2,*])^2.d0*a_disk[*,*,*,2]+ (one_minus_cosi_mat[1:*,0:-2,*])^3.d0*a_disk[*,*,*,3]+$
                        (one_minus_cosi_mat[1:*,0:-2,*])^4.d0*a_disk[*,*,*,4]+ (one_minus_cosi_mat[1:*,0:-2,*])^5.d0*a_disk[*,*,*,5]
 a_disk  = !null
 a_tdisk = rebin(reform(atdisk, ntau-1, nwavebands-1, 1, 6), ntau-1, nwavebands-1, Nmodels, 6)
 mlambda_td[1:*,0:-2,*]=(one_minus_cosi_mat[1:*,0:-2,*])^0.d0*a_tdisk[*,*,*,0]+(one_minus_cosi_mat[1:*,0:-2,*])^1.d0*a_tdisk[*,*,*,1]+$
                        (one_minus_cosi_mat[1:*,0:-2,*])^2.d0*a_tdisk[*,*,*,2]+(one_minus_cosi_mat[1:*,0:-2,*])^3.d0*a_tdisk[*,*,*,3]+$
                        (one_minus_cosi_mat[1:*,0:-2,*])^4.d0*a_tdisk[*,*,*,4]+(one_minus_cosi_mat[1:*,0:-2,*])^5.d0*a_tdisk[*,*,*,5]
 a_tdisk = !null
 a_bulge = rebin(reform(abulge, ntau-1, nwavebands-1, 1, 6), ntau-1, nwavebands-1, Nmodels, 6)
 mlambda_b[1:*,0:-2,*]= (one_minus_cosi_mat[1:*,0:-2,*])^0.d0*a_bulge[*,*,*,0]+(one_minus_cosi_mat[1:*,0:-2,*])^1.d0*a_bulge[*,*,*,1]+$
                        (one_minus_cosi_mat[1:*,0:-2,*])^2.d0*a_bulge[*,*,*,2]+(one_minus_cosi_mat[1:*,0:-2,*])^3.d0*a_bulge[*,*,*,3]+$
                        (one_minus_cosi_mat[1:*,0:-2,*])^4.d0*a_bulge[*,*,*,4]+(one_minus_cosi_mat[1:*,0:-2,*])^5.d0*a_bulge[*,*,*,5]
 a_bulge = !null
 one_minus_cosi_mat = !null

 ; Calculate delta_mlambda for entire galaxy using Equation 4 given in Doore et al. (2021)
 ; Equation is the intrinsic version of Equation 16 from Tuffs et al (2004)
 rbulge0_mat = rebin(reform(rbulge0, 1, 1, Nmodels), ntau, nwavebands, Nmodels)
 rdisk0_mat  = rebin(reform(rdisk0,  1, 1, Nmodels), ntau, nwavebands, Nmodels)
 att_bulge = rbulge0_mat*10.d0^(-1.d0*temporary(mlambda_b)/2.5d0)
 att_disk = rdisk0_mat*10.d0^(-1.d0*temporary(mlambda_d)/2.5d0)
 att_disk_bulge = temporary(att_disk) + temporary(att_bulge)

 att_tdisk_1 = (1-temporary(rdisk0_mat)-temporary(rbulge0_mat))*10.d0^(-1.d0*temporary(mlambda_td)/2.5d0)
 F_mat = rebin(reform(F_clump, 1, 1, Nmodels), ntau, nwavebands, Nmodels)
 flam_mat = rebin(reform(flambda, 1, nwavebands, 1), ntau, nwavebands, Nmodels)
 att_tdisk = (1-(temporary(F_mat)*temporary(flam_mat)))*temporary(att_tdisk_1)

 mlambda=-2.5d0*alog10(temporary(att_tdisk)+temporary(att_disk_bulge))
 ;There is no attenuation taub_f=0, non-zeros occur from F_clump in att_tdisk
 mlambda[0,*,*]=0.d0

 ;interpolate data for input tauB_f's and input wavelength's
 delta_m=dblarr(nwave, Nmodels)
 
 taub_loc = interpol(lindgen(Ntau), tau, tauB_f) > 0
 wave_loc = interpol(lindgen(nwavebands), wavebands, (reform(wave)*10000.d0)) > 0

 for i=0,(Nmodels-1) do $
   delta_m[*,i] = INTERPOLATE(mlambda[*,*,i], taub_loc[i], wave_loc, /double, /grid)

 ;Set attenuation to 0 at wavelengths greater than 5.0 microns. Do this in case interpolation has error.
 delta_m[where(reform(wave) gt 5.d0 or reform(wave) lt 9.12d-2,/null),*]=0.0  

 exp_neg_tau=10.d0^(-0.4d0*delta_m)

 if Nmodels eq 1 and Nwave eq 1 then return, exp_neg_tau[0] else $
   return, exp_neg_tau

end
