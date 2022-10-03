pro generate_doore21_lbol_abs_table, input_dir, redshifts, config
;+
; Name
; ----
;   GENERATE_DOORE21_LBOL_ABS_TABLE
;
; Purpose
; -------
;   Generates and saves pre-computed table(s) containing the bolometric luminosity
;   of the stellar emission absorbed by dust (:math:`L_{\rm bol}^{\rm abs}``) for
;   the Doore et al. (2021) attenuation curves to have energy conservation. Since
;   the Doore et al. (2021) attenuation curves are inclination dependent, 
;   :math:`L_{\rm bol}^{\rm abs}`` cannot be simply computed as the difference
;   between the attenuated and unattenuated bolometric luminosities. Rather, the
;   attenuated bolometric luminosities must first be integrated over inclination.
;   Therefore, to save computational time, this process is performed before fitting
;   rather than during (see Section 4.4 of Doore et al. 2021 for further details).
;
; Calling Sequence
; ----------------
;   ::
;
;       generate_doore21_lbol_abs_table, input_dir , redshifts, config
;
; Inputs
; ------
;   ``input_dir`` : string scalar
;       The path to the file containing the input SED data.
;   ``redshifts`` : int, float, or double array(Nred)
;       The redshifts at which to compute the tables.
;   ``config`` : structure
;       A Lightning configuration structure. (See
;       ``lightning_configure_defaults.pro`` for details and contents.)
;
; Output
; ------
;   A FITS file per redshift containing a structure that has the pre-computed
;   table of :math:`L_{\rm bol}^{\rm abs}`` to have energy conservation with
;   the Doore et al. (2021) attenuation curves. Additionally, the parameter
;   grids and ``step_bounds`` are included in the structure. Files are saved
;   in the directory ``<input_dir>/lightning_output/doore21_Lbol_abs_table/``.
;   The full description of the structure is as follows:
;
;   ==============     =============================     ========================================================
;   TAG                TYPE                              DESCRIPTION
;   ==============     =============================     ========================================================
;   LBOL_ABS_MODEL     double(Nsteps, Ntaubf, Nn, 2)     Absorbed bolometric stellar luminosity :math:`[L_\odot]`
;   STEPS_BOUNDS       double(Nsteps+1)                  Age bounds :math:`[\rm yr]`
;   NSTEPS             long                              Number of age bins (steps)
;   TAUB_F_GRID        double(Ntaubf)                    Gridded values of ``tauB_f``
;   RDISK0_GRID        double(Nn)                        Gridded values of ``rdisk0``
;   F_GRID             double(Nn)                        Gridded values of ``F_clump``
;   ==============     =============================     ========================================================
;
; Notes
; -----
;   The last dimension in the ``LBOL_ABS_MODEL`` tag indicates the young (``0``) and
;   old populations (``1``), respectively. Since the young population ignores ``rdisk0``
;   and the old population ignores ``F_clump``, this minimizes memory used to save
;   the files, while still encapsulating all parameter space.
;
; Reference
; ---------
;   `Doore, K., Eufrasio, R. T., Lehmer, B. D., et al. 2021, ApJ, 923, 26 <https://ui.adsabs.harvard.edu/abs/2021ApJ...923...26D/abstract>`_
;
; Modification History
; --------------------
;   - 2022/04/19: Created (Keith Doore)
;   - 2022/05/17: Changed to only allow for specified redshift array (Keith Doore)
;   - 2022/05/17: Changed ``redshift`` and ``steps_bounds`` to be required inputs (Keith Doore)
;   - 2022/05/17: Replaced ``steps_bounds`` with configuration structure (Keith Doore)
;   - 2022/05/18: Allowed for unique cosmologies (Keith Doore)
;   - 2022/06/29: Replaced ``!cv`` with ``!lightning_cgs`` (Keith Doore)
;   - 2022/06/29: Fixed error where ``steps_bounds`` was not retrieved from ``config`` (Keith Doore)
;   - 2022/07/08: Fixed issue where age bins were not truncated to age of universe at ``z=0`` (Keith Doore)
;   - 2022/08/18: Added progress printing (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(input_dir) eq 0 then message, 'Variable is undefined: INPUT_DIR'
 if size(input_dir, /type) ne 7 then message, 'INPUT_DIR must be of type string.'
 if size(input_dir, /n_dim) ne 0 then message, 'INPUT_DIR must be a scalar.'

 if n_elements(redshifts) eq 0 then message, 'Variable is undefined: REDSHIFTS'
 if size(redshifts, /type) lt 2 or size(redshifts, /type) gt 5 then $
          message, 'REDSHIFTS must be of type int, float, or double.'
 if size(reform(redshifts), /n_dim) ne 1 then message, 'REDSHIFTS must be a 1-D array.'
 if total(redshifts lt 0) gt 0 then message, 'REDSHIFTS must only contain non-negative values.'

 if n_elements(config) eq 0 then message, 'Variable is not defined: CONFIG.'
 if size(config, /type) ne 8 then message, 'CONFIG is not of type structure.'



 ; Add redshift zero since it is the base redshift when steps_bounds is unaltered.
 redshifts = [0.d0, reform(redshifts)]
 ; Remove duplicate redshifts
 redshifts = redshifts[uniq(redshifts, sort(redshifts))]
 steps_bounds = config.STEPS_BOUNDS

 ; Remove old Lbol_abs_tables and generate new directory
 lbol_abs_table_dir = input_dir+'lightning_output/doore21_Lbol_abs_table/'
 file_delete, lbol_abs_table_dir, /recursive, /allow_nonexistent
 temp_dir = lbol_abs_table_dir+'temp_files/'
 file_mkdir, temp_dir


; The parameter space density that allows for minimal error when interpolating between
;  gridded values (see Section 4.3 of Doore+2021 for details)
 ntauB_f = 51
 nrdisk0 = 51
 nF      = 51
 ncosi   = 70
 taub_f_extra = [0.0d:0.1d:0.01d]
 taub_f_extra2 = [0.0d:0.32d:0.02d]

 ran_tauB_f  = [ 0.d, 8.0d]  & tauB_f  = ran_tauB_f[0] +dindgen(ntauB_f)*(ran_tauB_f[1] -ran_tauB_f[0] )/(ntauB_f-1)
 ran_rdisk0  = [ 0.d, 1.0d]  & rdisk0  = ran_rdisk0[0] +dindgen(nrdisk0)*(ran_rdisk0[1] -ran_rdisk0[0] )/(nrdisk0-1)
 ran_F_clump = [0.d, 0.61d]  & F_clump = ran_F_clump[0]+dindgen(nF)     *(ran_F_clump[1]-ran_F_clump[0])/(nF-1)
 ran_cosi    = [ 0.d, 1.0d]  & cosi    = ran_cosi[0]   +dindgen(ncosi)  *(ran_cosi[1]   -ran_cosi[0]   )/(ncosi-1)
 ; rbulge0 is tied to rdisk0 due the the binary parameter, rold0 (rold0 = rdisk0 + rbulge0).
 ;  Therefore it does not need to be computed.

 tauB_f = [tauB_f,tauB_f_extra,tauB_f_extra2]
 taub_f = taub_f[uniq(taub_f, sort(taub_f))]

 ntauB_f  = n_elements(tauB_f)
 nrdisk0  = n_elements(rdisk0)
 nF       = n_elements(F_clump)
 ncosi    = n_elements(cosi)

 ; create vectors of all possible models
 Nmodel = ntaub_f*ncosi*nF*2L
 taub_f_vector  = reform(rebin(reform(taub_f,  ntaub_f, 1, 1, 1), ntaub_f, ncosi, nF, 2), Nmodel)
 cosi_vector    = reform(rebin(reform(cosi,      1, ncosi, 1, 1), ntaub_f, ncosi, nF, 2), Nmodel)
 F_clump_vector = reform(rebin(reform(F_clump,   1,    1, nF, 1), ntaub_f, ncosi, nF, 2), Nmodel)
 ; Convert rdisk0 to b_to_d ratio
 b_to_d_vector  = reform(rebin(reform(1.0/rdisk0-1, 1, 1, nF, 1), ntaub_f, ncosi, nF, 2), Nmodel)
 ; Rebin and reform replaces Inifinity with -NaN. Replace with Infinity
 b_to_d_vector[where(finite(b_to_d_vector,/nan) eq 1)]=!values.D_infinity
 rold0_vector   = reform(rebin(reform([0.0, 1.0], 1,  1,  1, 2),  ntaub_f, ncosi, nF, 2), Nmodel)

 ; Need stellar models grid of wavelength
 stellar_models = binned_stellar_models(filter_labels='WFC3_F125W', $
                                        redshift=0.0, _extra=config)

 exp_neg_tau    = doore21_atten(stellar_models.wave_rest, $
                                tauB_f=taub_f_vector,     $
                                rold0=rold0_vector,       $
                                b_to_d=b_to_d_vector,     $
                                F_clump=F_clump_vector,   $
                                cosi=cosi_vector)


; Loop through redshift and compute table of Lbol_abs.
;  Nonzero redshifts only have final adjusted age bin computed.
 t0 = systime(/sec)
 for k=0, (n_elements(redshifts)-1) do begin

   ; If redshift is not 0, only use the last age 2 age bins as they are the only ones that change.
   steps_bounds_temp = steps_bounds < galage(redshifts[k], 1e3, /silent, _extra=config.COSMOLOGY)
   steps_bounds_temp = steps_bounds_temp[uniq(steps_bounds_temp, sort(steps_bounds_temp))]
   if redshifts[k] ne 0 then $
     steps_bounds_temp = steps_bounds_temp[-2:-1]

   Nsteps = n_elements(steps_bounds_temp) - 1

   ; Can't pass config by keyword inheritance as it will override steps_bounds_temp
   stellar_models = binned_stellar_models(filter_labels='WFC3_F125W', steps_bounds=steps_bounds_temp, $
                                          redshift=redshifts[k], dtime_SF=config.dtime_SF, Zmetal=config.Zmetal, $
                                          IMF=config.IMF, cosmology=config.cosmology, $
                                          no_emission_lines=config.no_emission_lines, $
                                          no_nebular_extinction=config.no_nebular_extinction)

   steps_Lbol_abs = dblarr(Nsteps, ntaub_f, nF, 2)

   for st=0,(Nsteps-1) do begin

     steps_Lbol_abs_inc = dblarr(Nmodel)
     steps_Lnu=rebin(stellar_models.lnu[*,st], n_elements(stellar_models.wave_obs), Nmodel)
     steps_Lnu_red = exp_neg_tau*steps_Lnu
     Lbol = (stellar_models.Lbol[st])[0]

     for j=0,(Nmodel-1) do begin
       ; gives negative value from integral bc nu_obs is reversed (big to small)
        steps_Lbol_abs_inc[j] = trapez(steps_Lnu_red[*,j], 1.d4*!lightning_cgs.clight/stellar_models.wave_obs)
     endfor

     steps_Lbol_abs_inc_reformed = reform(steps_Lbol_abs_inc, ntaub_f, ncosi, nF, 2)
     for j=0,(ntaub_f-1) do begin
       for i=0,(nF-1) do begin
         for k1=0,(2-1) do begin
            ; Need negative since acos(cosi) is reversed (big to small)
   	        steps_Lbol_abs[st,j,i,k1]=-1.d0*trapez((Lbol+steps_Lbol_abs_inc_reformed[j,*,i,k1])*sin(acos(cosi)),$
                                                   acos(cosi))
         endfor
       endfor
     endfor

   endfor

   Lstar_abs = temporary(steps_Lbol_abs)
   ; Differences below 0.000001 in redshift result in insignificant changes in steps_bounds.
   ;   Therefore, use this as the precision in the save file name.
   save, Lstar_abs, filename=temp_dir+'temp_z_'+string(redshifts[k],f='(f0.6)')+'.sav'

   ; If progress printing, print progress bar.
   if config.PRINT_PROGRESS then lightning_print_progress, k, n_elements(redshifts), t0, funcname='GENERATE_DOORE21_LBOL_ABS_TABLE'

 endfor

; Compile nonzero redshift tables with beginning ages of zero redshift table and save
 for i=0,(n_elements(redshifts)-1) do begin

   steps_bounds_temp = steps_bounds < galage(redshifts[i], 1e3, /silent, _extra=config.COSMOLOGY)
   steps_bounds_temp = steps_bounds_temp[uniq(steps_bounds_temp, sort(steps_bounds_temp))]

   Nsteps = n_elements(steps_bounds_temp) - 1

   steps_Lbol_abs = dblarr(Nsteps, ntaub_f, nF, 2)
   steps_Lbol_abs = dblarr(Nsteps, ntaub_f, nF, 2)

   restore, temp_dir+'temp_z_'+string(redshifts[i],f='(f0.6)')+'.sav'
   steps_Lbol_abs[-1,*,*,*] = Lstar_abs[-1,*,*,*]

   restore,filename=temp_dir+'temp_z_'+string(redshifts[0],f='(f0.6)')+'.sav'
   steps_Lbol_abs[0:Nsteps-2,*,*,*] = Lstar_abs[0:Nsteps-2,*,*,*]

   Lbol_abs_table={Lbol_abs_model:steps_Lbol_abs,  $
                   steps_bounds:steps_bounds_temp, $
                   nsteps:nsteps,                  $
                   tauB_f_grid:tauB_f,             $
                   rdisk0_grid:rdisk0,             $
                   F_grid:F_clump                  $
                   }


   mwrfits, Lbol_abs_table, lbol_abs_table_dir+'redshift_'+string(redshifts[i],f='(f0.6)')+'.fits',/create

 endfor

 file_delete, temp_dir, /recursive

end