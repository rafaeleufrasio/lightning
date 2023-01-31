function lightning_configure_interactive, config_edit=config_edit
;+
; Name
; ----
;   LIGHTNING_CONFIGURE_INTERACTIVE
;
; Purpose
; -------
;   Generates the Lightning configuration structure via interactive
;   terminal prompts. Inputs are checked for errors to make sure all
;   inputs are correctly formatted and can be understood by Lightning.
;   Additionally, already existing configuration structures can be edited
;   interactively if it is optionally input. Further details can be
;   found at :ref:`configure-setting-label`.
;
; Calling Sequence
; ----------------
;   ::
;
;       config = lightning_configure_interactive([config_edit = ])
;
; Optional Input
; --------------
;   ``config_edit`` : structure
;       A Lightning configuration structure that is to be
;       interactively edited.
;
; Output
; ------
;   ``config`` : structure
;       A Lightning configuration structure. (See
;       ``lightning_configure_defaults.pro`` for details and contents.)
;
; Modification History
; --------------------
;   - 2022/04/28: Created (Keith Doore)
;   - 2022/05/05: Allowed for interactively editing structure (Keith Doore)
;   - 2022/05/16: Added ability to choose no stellar model (Keith Doore)
;   - 2022/05/17: Added ability to choose cosmology (Keith Doore)
;   - 2022/05/20: Removed loguniform and lognormal prior options (Keith Doore)
;   - 2022/08/30: Placed cosmology parameters at end of Core section for online documentation purposes (Keith Doore)
;   - 2022/09/01: Replaced ``XRAY_STAT`` with ``XRAY_UNC`` (Erik B. Monson)
;   - 2022/09/14: Updates to allow fitting with X-ray fluxes (Erik B. Monson)
;   - 2022/10/24: Added option to choose stranded walker deviation value for affine MCMC (Keith Doore)
;   - 2022/10/25: Renamed SPS to SSP (Keith Doore)
;   - 2022/12/13: Prevented ``XRAY_UNC`` from begin set if ``XRAY_UNIT='FLUX'`` (Keith Doore)
;   - 2023/01/31: Added ``OUTPUT_FILENAME`` option to allow for setting of post-processed filename (Keith Doore)
;-
 On_error, 2
 Compile_opt idl2

; Error handling
 if n_elements(config_edit) ne 0 then $
   if size(config_edit, /type) ne 8 then message, 'config_edit must be of type structure.'

; Generate a orderedhash that is used to dynamically make configuration structure
 config = orderedhash()

; Check if existing file is to be edited
 if n_elements(config_edit) eq 0 then edit = 0 else edit = 1


 if edit then begin
   instruct = 'Welcome to the interactive version of the Lightning SED fitting code. '+$
              'The following prompts will allow you to edit your previous configuration '+$
              'of Lightning. These prompts allow you to change your desired model '+$
              'configuration and corresponding parameter priors. Additionally, '+$
              'you can choose a different fitting method and corresponding '+$
              'hyperparameters. The previous configuration choices will be given '+$
              'at each prompt. For ease of use, these choices can be reselected '+$
              'by simply leaving the prompt blank and pressing <ENTER>. '+$
              'For a more detailed description of all configuration settings than those '+$
              'given at each prompt, please see the online documentation.'
 endif else begin
   instruct = 'Welcome to the interactive version of the Lightning SED fitting code. '+$
              'The following prompts will lead you through configuring Lightning. '+$
              'These prompts allow you to choose your desired model '+$
              'configuration and corresponding parameter priors. Additionally, '+$
              'you will specify your desired fitting method and corresponding '+$
              'hyperparameters. If you are unsure on what values to choose for each '+$
              'input, default/recommended values for most configuration options are given '+$
              'at each prompt. For ease of use, the default/recommended value can be '+$
              'chosen by simply leaving the prompt blank and pressing <ENTER>. '+$
              'For a more detailed description of all configuration settings than those '+$
              'given at each prompt, please see the online documentation.'
 endelse
 for i=0,(terminal_size())[1] do print_to_width, ''
 print_to_width, 'INSTRUCTIONS'
 print_to_width, '============'
 print_to_width, instruct
 for i=0,2 do print_to_width, ''
 print_to_width, 'BEGINNING CONFIGURATION'

 prior_options = ['fixed', 'uniform', 'normal', 'tabulated']
;========================================    CORE    =========================================================
; Core options to Lightning
 ;============ Output Filename ============
 error = 1
 default_val = 'postprocessed_data_%'
 default = '(Default: '+default_val+')'
 if edit then begin
   default_val = config_edit.OUTPUT_FILENAME
   default = '(Previous configuration: '+default_val+')'
 endif
 filename_message = 'What name (without the file extension suffix) would you like to give '+$
                    'to the FITS file containing the output post-processed data. NOTE: A UTC '+$
                    'timestamp can be automatically included in the filename so that you '+$
                    'can have a unique filename for multiple repeat runs to prevent accidentally '+$
                    'overwriting old runs. This is done by including a single % character in '+$
                    'the filename where you want the timestamp to appear. '+default
 print_to_width, '======================='
 print_to_width, ''
 print_to_width, filename_message

 temp = ''
 read,temp, prompt='OUTPUT_FILENAME: ',form='(A0)'
 temp = temp
 if temp eq '' then temp = default_val
 config['OUTPUT_FILENAME'] = temp


 ;============ Print Progress ============
 error = 1
 default_val = 'YES'
 default = '(Default: '+strlowcase(default_val)+')'
 if edit then begin
   if config_edit.PRINT_PROGRESS eq 1 then default_val = 'YES' else default_val = 'NO'
   default = '(Previous configuration: '+strlowcase(default_val)+')'
 endif
 print_message = 'Would you like the SED fitting progress and expected fitting time '+$
                 'remaining to be printed to the terminal? Please indicate yes or no. '+$
                 default
 print_to_width, '======================='
 while error do begin
   print_to_width, ''
   print_to_width, print_message

   temp = ''
   read,temp, prompt='PRINT_PROGRESS: ',form='(A0)'
   temp = strtrim(temp, 2)
   if temp eq '' then temp = default_val

   if ~(strupcase(temp) eq 'YES' or strupcase(temp) eq 'NO') then begin
     print_message = 'Please type yes or no.'
   endif else error = 0
 endwhile
 if strupcase(temp) eq 'YES' then temp = 1 else temp = 0
 config['PRINT_PROGRESS'] = temp


 ;============ Max CPUs ============
 error = 1
 default_val = strtrim(string(!CPU.hw_ncpu, f='(I0)'), 2)
 default = '(Recommended: '+default_val+')'
 if edit then begin
   default_val = strtrim(string(config_edit.MAX_CPUS,f='(I0)'), 2)
   default = '(Previous configuration: '+default_val+')'
 endif
 cpu_message = 'Please specify the maximum number of CPUs to utilize. If this value '+$
               'exceeds the number of CPUs on the machine, then all CPUs will be used. '+$
               'Current number of CPUs available on this machine: '+$
               strtrim(string(!CPU.hw_ncpu,f='(I0)'), 2)+'. '+$
               default
 print_to_width, '======================='
 while error do begin
   print_to_width, ''
   print_to_width, cpu_message

   temp = ''
   read,temp, prompt='MAX_CPUS: ',form='(A0)'
   if temp eq '' then temp = default_val
   temp = fix(temp)

   if temp le 0 then begin
     cpu_message = 'Please specify a number of CPUs of at least 1.'
   endif else error = 0
 endwhile
 config['MAX_CPUS'] = temp


 ;============ Energy Balance ============
 error = 1
 default_val = 'YES'
 default = '(Default: '+strlowcase(default_val)+')'
 if edit then begin
   if config_edit.ENERGY_BALANCE eq 1 then default_val = 'YES' else default_val = 'NO'
   default = '(Previous configuration: '+strlowcase(default_val)+')'
 endif
 energy_message = 'Would you like the total integrated IR luminosity (normalization) of '+$
                  'the dust emission to be tied to the total absorbed stellar (and, if set, '+$
                  'AGN) emission? Please indicate yes or no. '+$
                  default
 print_to_width, '======================='
 while error do begin
   print_to_width, ''
   print_to_width, energy_message

   temp = ''
   read, temp, prompt='ENERGY_BALANCE: ',form='(A0)'
   temp = strtrim(temp, 2)
   if temp eq '' then temp = default_val

   if ~(strupcase(temp) eq 'YES' or strupcase(temp) eq 'NO') then begin
     energy_message = 'Please type yes or no.'
   endif else error = 0
 endwhile
 if strupcase(temp) eq 'YES' then temp = 1 else temp = 0
 config['ENERGY_BALANCE'] = temp


 ;============ Model Uncertainty ============
 error = 1
 default_val = '0.0'
 default = '(Default: '+default_val+')'
 if edit then begin
   default_val = strtrim(string(config_edit.MODEL_UNC,f='(F0)'), 2)
   default = '(Previous configuration: '+default_val+')'
 endif
 mod_unc_message = 'Please specify the fractional model uncertainty to use. This '+$
                   'uncertainty is assumed for all bandpasses. If no model uncertainty '+$
                   'is to be included, set this value to 0. '+$
                   default
 print_to_width, '======================='
 while error do begin
   print_to_width, ''
   print_to_width, mod_unc_message

   temp = ''
   read,temp, prompt='MODEL_UNC: ',form='(A0)'
   if temp eq '' then temp = default_val
   temp = double(temp)

   if temp lt 0 then begin
     mod_unc_message = 'The model uncertainty must be zero or a positive value.'
   endif else error = 0
 endwhile
 config['MODEL_UNC'] = temp


 ;============ Cosmology ============
 cosmo_message = 'Please specify the cosmology parameters to use in the SED fitting. '+$
                 'These parameters set the age of the universe (and distance to the '+$
                 'object if the object has its distance specified by redshift).'
 print_to_width, '======================='
 print_to_width, ''
 print_to_width, cosmo_message
    ;========= H0 =============
     error = 1
     default_val = string(70, f='(F0.3)')
     default = '(Default: '+default_val+')'
     if edit then begin
       default_val = string(config_edit.COSMOLOGY.H0, f='(F0.3)')
       default = '(Previous configuration: '+default_val+')'
     endif
     cosmo_message = 'Please specify H0, the Hubble constant in km/s/Mpc. '+default

     while error do begin
       print_to_width, ''
       print_to_width, cosmo_message

       temp = ''
       read, temp, prompt='H0: ',form='(A0)'
       if temp eq '' then temp = default_val
       temp = double(temp)

       if temp le 0 then begin
         cosmo_message = 'H0 must be a postive value.'
       endif else error = 0
     endwhile
     config['COSMOLOGY', 'H0'] = temp

    ;========= OMEGA_M =============
     error = 1
     default_val = string(0.3, f='(F0.3)')
     default = '(Default: '+default_val+')'
     if edit then begin
       default_val = string(config_edit.COSMOLOGY.OMEGA_M, f='(F0.3)')
       default = '(Previous configuration: '+default_val+')'
     endif
     cosmo_message = 'Please specify OMEGA_M, the Matter density, normalized to the '+$
                     'closure density. '+default

     while error do begin
       print_to_width, ''
       print_to_width, cosmo_message

       temp = ''
       read, temp, prompt='OMEGA_M: ',form='(A0)'
       if temp eq '' then temp = default_val
       temp = double(temp)

       if temp lt 0 then begin
         cosmo_message = 'OMEGA_M must be zero or positive.'
       endif else error = 0
     endwhile
     config['COSMOLOGY', 'OMEGA_M'] = temp

    ;========= LAMBDA0 =============
     error = 1
     default_val = string(0.7, f='(F0.3)')
     default = '(Default: '+default_val+')'
     if edit then begin
       default_val = string(config_edit.COSMOLOGY.LAMBDA0, f='(F0.3)')
       default = '(Previous configuration: '+default_val+')'
     endif
     cosmo_message = 'Please specify LAMBDA0, the cosmological constant, normalized to the '+$
                     'closure density. '+default

     print_to_width, ''
     print_to_width, cosmo_message

     temp = ''
     read, temp, prompt='LAMBDA0: ',form='(A0)'
     if temp eq '' then temp = default_val
     temp = double(temp)

     config['COSMOLOGY', 'LAMBDA0'] = temp

    ;========= Q0 =============
     error = 1
     default_val = string(-0.55, f='(F0.3)')
     default = '(Default: '+default_val+')'
     if edit then begin
       default_val = string(config_edit.COSMOLOGY.Q0, f='(F0.3)')
       default = '(Previous configuration: '+default_val+')'
     endif
     cosmo_message = 'Please specify Q0, the deceleration parameter. '+default

     print_to_width, ''
     print_to_width, cosmo_message

     temp = ''
     read, temp, prompt='Q0: ',form='(A0)'
     if temp eq '' then temp = default_val
     temp = double(temp)

     config['COSMOLOGY', 'Q0'] = temp

    ;========= K =============
     error = 1
     default_val = string(0.0, f='(F0.3)')
     default = '(Default: '+default_val+')'
     if edit then begin
       default_val = string(config_edit.COSMOLOGY.K, f='(F0.3)')
       default = '(Previous configuration: '+default_val+')'
     endif
     cosmo_message = 'Please specify K, the curvature constant, normalized to the '+$
                     'closure density. '+default

     print_to_width, ''
     print_to_width, cosmo_message

     temp = ''
     read, temp, prompt='K: ',form='(A0)'
     if temp eq '' then temp = default_val
     temp = double(temp)

     config['COSMOLOGY', 'K'] = temp



;==================================    STELLAR EMISSION    ===================================================
; Stellar emission module options
 ;============ SSP ============
 ssp_options = ['PEGASE', 'none']
 error = 1
 default_val = 'PEGASE'
 default = '(Default: '+default_val+')'
 if edit then begin
   default_val = strupcase(config_edit.SSP)
   default = '(Previous configuration: '+default_val+')'
 endif
 ssp_message = 'Please specify the stellar population synthesis (SSP) models to use '+$
               'for the stellar population. For no stellar emission model, set '+$
               "to 'NONE'. Current options: "+strjoin(ssp_options, ', ')+'. '+$
               default
 print_to_width, '======================='
 while error do begin
   print_to_width, ''
   print_to_width, ssp_message

   temp = ''
   read,temp, prompt='SSP: ',form='(A0)'
   if temp eq '' then temp = default_val
   temp = strtrim(temp, 2)

   if total(strupcase(temp) eq strupcase(ssp_options)) eq 0 then begin
     ssp_message = 'Please specify one of the current options: '+strjoin(ssp_options, ', ')+'.'
   endif else error = 0
 endwhile
 config['SSP'] = strupcase(temp)


 case strupcase(config['SSP']) of
   'PEGASE': begin
      ;============ IMF ============
      imf_options = ['Kroupa01']
      error = 1
      default_val = 'Kroupa01'
      default = '(Default: '+default_val+')'
      if edit then begin
        if strupcase(config_edit.SSP) eq 'PEGASE' then begin
          default_val = strupcase(config_edit.IMF)
          default = '(Previous configuration: '+default_val+')'
        endif
      endif
      imf_message = 'Please specify the initial mass function (IMF) to use in the SSP models. '+$
                    'Current options: '+strjoin(imf_options, ', ')+'. '+$
                    default
      print_to_width, '======================='
      while error do begin
        print_to_width, ''
        print_to_width, imf_message

        temp = ''
        read,temp, prompt='IMF: ',form='(A0)'
        if temp eq '' then temp = default_val
        temp = strtrim(temp, 2)

        if total(strupcase(temp) eq strupcase(imf_options)) eq 0 then begin
          imf_message = 'Please specify one of the current options: '+strjoin(imf_options, ', ')+'.'
        endif else error = 0
      endwhile
      config['IMF'] = strupcase(temp)


      ;============ ZMETAL ============
      zmetal_options = [0.001d, 0.004d, 0.008d, 0.02d, 0.05d, 0.1d]
      error = 1
      default_val = '0.02'
      default = '(Default: '+default_val+')'
      if edit then begin
        if strupcase(config_edit.SSP) eq 'PEGASE' then begin
          default_val = strtrim(string(config_edit.ZMETAL, f='(F0.3)'), 2)
          default = '(Previous configuration: '+default_val+')'
        endif
      endif
      metal_message = 'Please specify the metallicity to use in the SSP models in terms of Z [Zsun = 0.02]. '+$
                      'Current options: '+strjoin(string(zmetal_options, f='(F0.3)'), ', ')+'. '+$
                      default
      print_to_width, '======================='
      while error do begin
        print_to_width, ''
        print_to_width, metal_message

        temp = ''
        read,temp, prompt='ZMETAL: ',form='(A0)'
        if temp eq '' then temp = default_val
        temp = double(temp)

        if total(temp eq zmetal_options) eq 0 then begin
          metal_message = 'Please specify one of the current options: '+$
                          strjoin(string(zmetal_options, f='(F0.3)'), ', ')+'.'
        endif else error = 0
      endwhile
      config['ZMETAL'] = temp


      ;============ EMISSION_LINES ============
      error = 1
      default_val = 'YES'
      default = '(Default: '+strlowcase(default_val)+')'
      if edit then begin
        if strupcase(config_edit.SSP) eq 'PEGASE' then begin
          if config_edit.EMISSION_LINES eq 1 then default_val = 'YES' else default_val = 'NO'
          default = '(Previous configuration: '+strlowcase(default_val)+')'
        endif
      endif
      emis_line_message = 'Would you like the SSP models to include nebular emission lines? '+$
                          'Please indicate yes or no. '+$
                          default
      print_to_width, '======================='
      while error do begin
        print_to_width, ''
        print_to_width, emis_line_message

        temp = ''
        read,temp, prompt='EMISSION_LINES: ',form='(A0)'
        if temp eq '' then temp = default_val
        temp = strtrim(temp, 2)

        if ~(strupcase(temp) eq 'YES' or strupcase(temp) eq 'NO') then begin
          emis_line_message = 'Please type yes or no.'
        endif else error = 0
      endwhile
      if strupcase(temp) eq 'YES' then temp = 1 else temp = 0
      config['EMISSION_LINES'] = temp


      ;============ NEBULAR_EXTINCTION ============
      error = 1
      default_val = 'YES'
      default = '(Default: '+strlowcase(default_val)+')'
      if edit then begin
        if strupcase(config_edit.SSP) eq 'PEGASE' then begin
          if config_edit.NEBULAR_EXTINCTION eq 1 then default_val = 'YES' else default_val = 'NO'
          default = '(Previous configuration: '+strlowcase(default_val)+')'
        endif
      endif
      neb_ext_message = 'Would you like the SSP models to include nebular extinction? '+$
                        'Please indicate yes or no. '+$
                        default
      print_to_width, '======================='
      while error do begin
        print_to_width, ''
        print_to_width, neb_ext_message

        temp = ''
        read,temp, prompt='NEBULAR_EXTINCTION: ',form='(A0)'
        if temp eq '' then temp = default_val
        temp = strtrim(temp, 2)

        if ~(strupcase(temp) eq 'YES' or strupcase(temp) eq 'NO') then begin
          neb_ext_message = 'Please type yes or no.'
        endif else error = 0
      endwhile
      if strupcase(temp) eq 'YES' then temp = 1 else temp = 0
      config['NEBULAR_EXTINCTION'] = temp


      ;============ SFH ============
      sfh_options = ['Non-Parametric']
      error = 1
      default_val = 'Non-Parametric'
      default = '(Default: '+default_val+')'
      if edit then begin
        if strupcase(config_edit.SSP) eq 'PEGASE' then begin
          default_val = config_edit.SFH
          default = '(Previous configuration: '+config_edit.SFH+')'
        endif
      endif
      sfh_message = 'Please specify the type of SFH to assume when fitting the SEDs. '+$
                    'Current options: '+strjoin(sfh_options, ', ')+'. '+$
                    default
      print_to_width, '======================='
      while error do begin
        print_to_width, ''
        print_to_width, sfh_message

        temp = ''
        read,temp, prompt='SFH: ',form='(A0)'
        if temp eq '' then temp = strupcase(default_val)
        temp = strtrim(temp, 2)

        if total(strupcase(temp) eq strupcase(sfh_options)) eq 0 then begin
          sfh_message = 'Please specify one of the current options: '+strjoin(sfh_options, ', ')+'.'
        endif else error = 0
      endwhile
      config['SFH'] = strupcase(temp)


      case strupcase(config['SFH']) of
      ; ===================    Non-Parametric SFH   =============================================
        'NON-PARAMETRIC' : begin
          ;============ Number of SFH steps ============
           error = 1
           default_val = '5'
           default = '(Default: '+default_val+')'
           if edit then begin
             if strupcase(config_edit.SSP) eq 'PEGASE' then begin
               if strupcase(config_edit.SFH) eq 'NON-PARAMETRIC' then begin
                 default_val = strtrim(string(n_elements(config_edit.STEPS_BOUNDS)-1, f='(I0)'), 2)
                 default = '(Previous configuration: '+default_val+')'
               endif
             endif
           endif
           nsteps_message = 'Please specify the number of non-parametric SFH age bins (or steps) '+$
                            'to use. '+default
           print_to_width, '======================='
           while error do begin
             print_to_width, ''
             print_to_width, nsteps_message

             temp = ''
             read,temp, prompt='NSTEPS: ',form='(A0)'
             if temp eq '' then temp = default_val
             temp = long(double(temp))

             if temp le 0 then begin
               nsteps_message = 'The number of steps must be 1 or more.'
             endif else error = 0
           endwhile
           Nsteps = temp


          ;============ Steps Bounds ============
           error = 1
           if nsteps eq 1 then default_val = [0, 13.6d9] else $
             default_val = [0, 10.d^(7+dindgen(nsteps)*(alog10(13.6d9)-7)/(nsteps-1))]
           default_val = strjoin(string(default_val, f='(E0.3)'), ', ')
           default = '(Default: '+default_val+')'
           if edit then begin
             if strupcase(config_edit.SSP) eq 'PEGASE' then begin
                if strupcase(config_edit.SFH) eq 'NON-PARAMETRIC' then begin
                  if n_elements(config_edit.STEPS_BOUNDS)-1 eq nsteps then begin
                    default_val = strjoin(string(config_edit.STEPS_BOUNDS, f='(E0.3)'), ', ')
                    default = '(Previous configuration: '+default_val+')'
                    temp = ''
                  endif
                endif
             endif
           endif
           steps_message = 'Please specify the age bin (or step) boundaries to use in the '+$
                           'non-parametric SFH. If an age bin contains ages older than the '+$
                           "universe at an input SED's redshift, the age bin upper bound will "+$
                           'be automatically adjusted to the age of the universe at that '+$
                           'redshift. If an entire age bin is older than universe at that '+$
                           'redshift, then the entire age bin will be omitted and the next '+$
                           'younger bin will be adjusted accordingly. Values must be given in '+$
                           'ascending order and units of years. Please separate values with '+$
                           'a space or comma. Number of expected values: '+$
                           strtrim(string(Nsteps+1, f='(I0)'), 2)+'. '+default
           print_to_width, '======================='
           while error do begin
             print_to_width, ''
             print_to_width, steps_message

             temp = ''
             read,temp, prompt='STEPS_BOUNDS: ',form='(A0)'
             if temp eq '' then temp = default_val
             temp = double(STRSPLIT(temp, ' ,', /EXTRACT))

             case 1 of
               total(temp lt 0) gt 0: $
                 steps_message = 'Please specify only zero or positive values.'

               n_elements(temp) ne Nsteps+1: $
                 steps_message = 'Please specify '+strtrim(string(Nsteps+1, f='(I0)'), 2)+ ' values.'

               total(sort(temp) ne indgen(n_elements(temp))) ne 0: $
                 steps_message = 'Please specify values in ascending order.'

               n_elements(uniq(temp, sort(temp))) ne n_elements(temp): $
                 steps_message = 'Please specify only unique values.'

               else: error = 0
             endcase
           endwhile
           config['STEPS_BOUNDS'] = temp


          ;============ Dtime_SF ============
           error = 1
           default_val = '5.d5'
           default = '(Default: '+default_val+')'
           if edit then begin
             if strupcase(config_edit.SSP) eq 'PEGASE' then begin
               if strupcase(config_edit.SFH) eq 'NON-PARAMETRIC' then begin
                 default_val = strtrim(string(config_edit.DTIME_SF, f='(E0.3)'), 2)
                 default = '(Previous configuration: '+default_val+')'
               endif
             endif
           endif
           dtime_message = 'Please specify the time step used for interpolating the SSP '+$
                           'models into the age bins in units of years. NOTE: We do not '+$
                           'recommend changing this value from its default, unless you '+$
                           'specified age bins with differences less than the default '+$
                           'value. '+default
           print_to_width, '======================='
           while error do begin
             print_to_width, ''
             print_to_width, dtime_message

             temp = ''
             read,temp, prompt='DTIME_SF: ',form='(A0)'
             if temp eq '' then temp = default_val
             temp = double(temp)

             if temp le 0 then begin
               dtime_message = 'The time step must be a positive value.'
             endif else error = 0
           endwhile
           config['DTIME_SF'] = temp


          ;============ PSI Priors ============
           default_priors = replicate('uniform', Nsteps)
           edit_temp = 0
           if edit then begin
             if strupcase(config_edit.SSP) eq 'PEGASE' then begin
               if strupcase(config_edit.SFH) eq 'NON-PARAMETRIC' then begin
                 if total(config_edit.STEPS_BOUNDS eq config['STEPS_BOUNDS']) eq Nsteps+1 then begin
                   default_priors = config_edit.PSI.PRIOR
                   edit_temp = 1
                 endif
               endif
             endif
           endif

           default_prior_args = [0.0d, 0.0, 1000, 1.0, 10.0]
           default_initialization = [0.0d, 1.d1]
           narg_max = 0
           prior = strarr(Nsteps)
           prior_arg = strarr(Nsteps, 4)
           initialization_range = dblarr(Nsteps, 2)
           for i=0, Nsteps-1 do begin
             if edit then begin
               if strupcase(config_edit.SSP) eq 'PEGASE' then begin
                 if strupcase(config_edit.SFH) eq 'NON-PARAMETRIC' then begin
                   if total(config_edit.STEPS_BOUNDS eq config['STEPS_BOUNDS']) eq Nsteps+1 then begin
                     if Nsteps gt 1 then prior_args = reform((config_edit.PSI.PRIOR_ARG)[i, *]) else $
                                         prior_args = config_edit.PSI.PRIOR_ARG
                     case strupcase(default_priors[i]) of
                       'FIXED':       default_prior_args = [prior_args, 0.0d, 1000, 1.0, 10.0]
                       'UNIFORM':     default_prior_args = [0.0d, prior_args, 1.0, 10.0]
                       'NORMAL':      default_prior_args = [0.0d, prior_args]
                       'TABULATED':
                     endcase
                     if Nsteps gt 1 then default_initialization = reform((config_edit.PSI.INITIALIZATION_RANGE)[i, *]) else $
                                         default_initialization = config_edit.PSI.INITIALIZATION_RANGE
                   endif
                 endif
               endif
             endif

             steps = [string(config['STEPS_BOUNDS', i], f='(E0.3)'), $
                      string(config['STEPS_BOUNDS', i+1], f='(E0.3)')]
             prior_arg_temp = prior_interactive('psi_'+strtrim(string(i+1, f='(I0)'), 2), $
                                       'the SFR of age bin '+strtrim(string(i+1, f='(I0)'), 2)+$
                                       ' ('+steps[0]+' - '+steps[1]+' Myr) in Msun yr-1. NOTE: All '+$
                                       ' or no priors must be tabulated', $
                                       [0.0d, !values.d_infinity], $
                                       default_priors[i], $
                                       default_prior_args, $
                                       default_initialization, $
                                       prior_options, prior=prior_temp, $
                                       initialization_range=initialization_range_temp, $
                                       edit=edit_temp)
             prior[i] = prior_temp
             initialization_range[i, *] = initialization_range_temp

             while total(total(strupcase(prior) eq 'TABULATED') eq [0, i+1]) eq 0 do begin
               print_to_width, ''
               print_to_width, '======================='
               print_to_width, '======================='
               print_to_width, 'WARNING! All or no psi priors must be tabulated. Previous prior choices: '+$
                      strjoin(prior, ', ')+'. Please redo the previous prior to meet '+$
                      'this requirement.'
               print_to_width, '======================='
               prior_arg_temp = prior_interactive('psi_'+strtrim(string(i+1, f='(I0)'), 2), $
                                    'the SFR of age bin '+strtrim(string(i+1, f='(I0)'), 2)+$
                                    ' ('+steps[0]+' - '+steps[1]+' Myr) in Msun yr-1. NOTE: All '+$
                                    ' or no priors must be tabulated.', $
                                    [0.0d, !values.d_infinity], $
                                    default_priors[i], $
                                    default_prior_args, $
                                    default_initialization, $
                                    prior_options, prior=prior_temp, $
                                    initialization_range=initialization_range_temp, $
                                    edit=edit_temp)
               prior[i] = prior_temp
               initialization_range[i, *] = initialization_range_temp
             endwhile

             if total(strupcase(prior[i]) eq ['FIXED', 'TABULATED']) eq 1 then narg_max = narg_max > 1
             if total(strupcase(prior[i]) eq ['UNIFORM']) eq 1 then narg_max = narg_max > 2
             if total(strupcase(prior[i]) eq ['NORMAL']) eq 1 then narg_max = narg_max > 4
             prior_arg[i, 0:(n_elements(prior_arg_temp)-1)] = prior_arg_temp
           endfor

           prior_arg = prior_arg[*, 0:(narg_max-1)]
           if total(strupcase(prior) eq 'TABULATED') eq 0 then prior_arg = double(prior_arg)
           if Nsteps eq 1 then prior = prior[0]
           config['PSI','PRIOR'] = strupcase(prior)
           config['PSI','PRIOR_ARG'] = reform(prior_arg)
           config['PSI','INITIALIZATION_RANGE'] = reform(initialization_range)
        end
      endcase
   end

   'NONE':

 endcase



;==================================    DUST ATTENUATION    ===================================================
; Dust attenuation module options
 ;============ ATTENUATION CURVE ============
 atten_curve_options = ['CALZETTI00', 'CALZETTI_MOD', 'DOORE21']
 error = 1
 default_val = 'CALZETTI00'
 default = '(Default: '+default_val+')'
 if edit then begin
   default_val = strupcase(config_edit.ATTEN_CURVE)
   default = '(Previous configuration: '+default_val+')'
 endif
 ; Check if stellar model is used, since DOORE21 requires it
 if strupcase(config['SSP']) eq 'NONE' then $
   atten_curve_options = ['CALZETTI00', 'CALZETTI_MOD']
 atten_message = 'Please specify the attenuation curve to apply to the stellar and/or AGN '+$
                 'models. Current options: '+strjoin(atten_curve_options, ', ')+'. '+$
                 default
 print_to_width, '======================='
 while error do begin
   print_to_width, ''
   print_to_width, atten_message

   temp = ''
   read,temp, prompt='ATTEN_CURVE: ',form='(A0)'
   if temp eq '' then temp = default_val
   temp = strtrim(temp, 2)

   if total(strupcase(temp) eq strupcase(atten_curve_options)) eq 0 then begin
     atten_message = 'Please specify one of the current options: '+strjoin(atten_curve_options, ', ')+'.'
   endif else error = 0
 endwhile
 config['ATTEN_CURVE'] = strupcase(temp)


 case strupcase(config['ATTEN_CURVE']) of
 ; ===================    Calzetti00   =============================================
   'CALZETTI00' : begin
      ;============ TauV prior ============
       default_priors = 'uniform'
       default_prior_args = [0.0d, 0.0, 10., 1.0, 1.0]
       default_initialization = [0.0d, 3]
       edit_temp = 0
       if edit then begin
         if strupcase(config_edit.ATTEN_CURVE) eq 'CALZETTI00' then begin
           default_priors = config_edit.TAUV.PRIOR
           prior_args = config_edit.TAUV.PRIOR_ARG
           case strupcase(default_priors) of
             'FIXED':       default_prior_args = [prior_args, 0.0d, 10., 1.0, 1.0]
             'UNIFORM':     default_prior_args = [0.0d, prior_args, 1.0, 1.0]
             'NORMAL':      default_prior_args = [0.0d, prior_args]
             'TABULATED':
           endcase
           default_initialization = config_edit.TAUV.INITIALIZATION_RANGE
           edit_temp = 1
         endif
       endif
       prior_arg = prior_interactive('tauV', 'the V-band optical depth', $
                                     [0.0d, !values.d_infinity], $
                                     default_priors, $
                                     default_prior_args, $
                                     default_initialization, $
                                     prior_options, prior=prior, $
                                     initialization_range=initialization_range, $
                                     edit=edit_temp)
       config['TAUV','PRIOR'] = strupcase(prior)
       config['TAUV','PRIOR_ARG'] = prior_arg
       config['TAUV','INITIALIZATION_RANGE'] = initialization_range
     end


 ; ===================    Modified Calzetti   =============================================
   'CALZETTI_MOD' : begin
      ;============ TauV_DIFF prior ============
       default_priors = 'uniform'
       default_prior_args = [0.0d, 0.0, 10., 1.0, 1.0]
       default_initialization = [0.0d, 3]
       edit_temp = 0
       if edit then begin
         if strupcase(config_edit.ATTEN_CURVE) eq 'CALZETTI_MOD' then begin
           default_priors = config_edit.TAUV_DIFF.PRIOR
           prior_args = config_edit.TAUV_DIFF.PRIOR_ARG
           case strupcase(default_priors) of
             'FIXED':       default_prior_args = [prior_args, 0.0d, 10., 1.0, 1.0]
             'UNIFORM':     default_prior_args = [0.0d, prior_args, 1.0, 1.0]
             'NORMAL':      default_prior_args = [0.0d, prior_args]
             'TABULATED':
           endcase
           default_initialization = config_edit.TAUV_DIFF.INITIALIZATION_RANGE
           edit_temp = 1
         endif
       endif
       prior_arg = prior_interactive('tauV_DIFF', 'the V-band optical depth of the diffuse dust', $
                                     [0.0d, !values.d_infinity], $
                                     default_priors, $
                                     default_prior_args, $
                                     default_initialization, $
                                     prior_options, prior=prior, $
                                     initialization_range=initialization_range, $
                                     edit=edit_temp)
       config['TAUV_DIFF','PRIOR'] = strupcase(prior)
       config['TAUV_DIFF','PRIOR_ARG'] = prior_arg
       config['TAUV_DIFF','INITIALIZATION_RANGE'] = initialization_range

      ;============ Delta prior ============
       default_priors = 'uniform'
       default_prior_args = [0.0, -2.3d, 0.4d, 0.0, 1.0]
       default_initialization = [-1.0d, 0.0d]
       edit_temp = 0
       if edit then begin
         if strupcase(config_edit.ATTEN_CURVE) eq 'CALZETTI_MOD' then begin
           default_priors = config_edit.DELTA.PRIOR
           prior_args = config_edit.DELTA.PRIOR_ARG
           case strupcase(default_priors) of
             'FIXED':       default_prior_args = [prior_args, -2.3d, 0.4d, 0.0, 1.0]
             'UNIFORM':     default_prior_args = [0.0d, prior_args, 0.0, 1.0]
             'NORMAL':      default_prior_args = [0.0d, prior_args]
             'TABULATED':
           endcase
           default_initialization = config_edit.DELTA.INITIALIZATION_RANGE
           edit_temp = 1
         endif
       endif
       prior_arg = prior_interactive('delta', 'the power law value to change the attenuation curve slope', $
                                     [-1.0d*!values.d_infinity, !values.d_infinity], $
                                     default_priors, $
                                     default_prior_args, $
                                     default_initialization, $
                                     prior_options, prior=prior, $
                                     initialization_range=initialization_range, $
                                     edit=edit_temp)
       config['DELTA','PRIOR'] = strupcase(prior)
       config['DELTA','PRIOR_ARG'] = prior_arg
       config['DELTA','INITIALIZATION_RANGE'] = initialization_range

      ;============ TauV_BC prior ============
       default_priors = 'fixed'
       default_prior_args = [0.0d, 0.0, 10., 1.0, 1.0]
       default_initialization = [0.0d, 3]
       edit_temp = 0
       if edit then begin
         if strupcase(config_edit.ATTEN_CURVE) eq 'CALZETTI_MOD' then begin
           default_priors = config_edit.TAUV_BC.PRIOR
           prior_args = config_edit.TAUV_BC.PRIOR_ARG
           case strupcase(default_priors) of
             'FIXED':       default_prior_args = [prior_args, 0.0d, 10., 1.0, 1.0]
             'UNIFORM':     default_prior_args = [0.0d, prior_args, 1.0, 1.0]
             'NORMAL':      default_prior_args = [0.0d, prior_args]
             'TABULATED':
           endcase
           default_initialization = config_edit.TAUV_BC.INITIALIZATION_RANGE
           edit_temp = 1
         endif
       endif
       prior_arg = prior_interactive('tauV_BC', 'the V-band optical depth of the birth cloud component', $
                                     [0.0d, !values.d_infinity], $
                                     default_priors, $
                                     default_prior_args, $
                                     default_initialization, $
                                     prior_options, prior=prior, $
                                     initialization_range=initialization_range, $
                                     edit=edit_temp)
       config['TAUV_BC','PRIOR'] = strupcase(prior)
       config['TAUV_BC','PRIOR_ARG'] = prior_arg
       config['TAUV_BC','INITIALIZATION_RANGE'] = initialization_range

      ;============ UV Bump ============
       error = 1
       default_val = 'YES'
       default = '(Default: '+strlowcase(default_val)+')'
       if edit then begin
         if strupcase(config_edit.ATTEN_CURVE) eq 'CALZETTI_MOD' then begin
           if config_edit.UV_BUMP eq 1 then default_val = 'YES' else default_val = 'NO'
           default = '(Previous configuration: '+strlowcase(default_val)+')'
         endif
       endif
       bump_message = 'Would you like a 2175 Angstrom UV bump feature added to '+$
                      'the attenuation curve? Please indicate yes or no. '+$
                      default
       print_to_width, '======================='
       while error do begin
         print_to_width, ''
         print_to_width, bump_message

         temp = ''
         read,temp, prompt='UV_BUMP: ',form='(A0)'
         if temp eq '' then temp = default_val
         temp = strtrim(temp, 2)

         if ~(strupcase(temp) eq 'YES' or strupcase(temp) eq 'NO') then begin
           bump_message = 'Please type yes or no.'
         endif else error = 0
       endwhile
       if strupcase(temp) eq 'YES' then temp = 1 else temp = 0
       config['UV_BUMP'] = temp
     end


 ; ===================    Doore21   =============================================
   'DOORE21' : begin
      ;============ TauB_F prior ============
       default_priors = 'uniform'
       default_prior_args = [0.0d, 0.0, 8.0, 4.0, 2.0]
       default_initialization = [2.0d, 6.0d]
       edit_temp = 0
       if edit then begin
         if strupcase(config_edit.ATTEN_CURVE) eq 'DOORE21' then begin
           default_priors = config_edit.TAUB_F.PRIOR
           prior_args = config_edit.TAUB_F.PRIOR_ARG
           case strupcase(default_priors) of
             'FIXED':       default_prior_args = [prior_args, 0.0d, 8.0, 4.0, 2.0]
             'UNIFORM':     default_prior_args = [0.0d, prior_args, 4.0, 2.0]
             'NORMAL':      default_prior_args = [0.0d, prior_args]
             'TABULATED':
           endcase
           default_initialization = config_edit.TAUB_F.INITIALIZATION_RANGE
           edit_temp = 1
         endif
       endif
       prior_arg = prior_interactive('tauB_F', 'the face-on optical depth in the B-band', $
                                     [0.0d, 8.0d], $
                                     default_priors, $
                                     default_prior_args, $
                                     default_initialization, $
                                     prior_options, prior=prior, $
                                     initialization_range=initialization_range, $
                                     edit=edit_temp)
       config['TAUB_F','PRIOR'] = strupcase(prior)
       config['TAUB_F','PRIOR_ARG'] = prior_arg
       config['TAUB_F','INITIALIZATION_RANGE'] = initialization_range

      ;============ F_clump prior ============
       default_priors = 'uniform'
       default_prior_args = [0.0, 0.0, 0.61d, 0.3d, 0.1d]
       default_initialization = [0.0, 0.61d]
       edit_temp = 0
       if edit then begin
         if strupcase(config_edit.ATTEN_CURVE) eq 'DOORE21' then begin
           default_priors = config_edit.F_CLUMP.PRIOR
           prior_args = config_edit.F_CLUMP.PRIOR_ARG
           case strupcase(default_priors) of
             'FIXED':       default_prior_args = [prior_args, 0.0d, 0.61d, 0.3d, 0.1d]
             'UNIFORM':     default_prior_args = [0.0d, prior_args, 0.3d, 0.1d]
             'NORMAL':      default_prior_args = [0.0d, prior_args]
             'TABULATED':
           endcase
           default_initialization = config_edit.F_CLUMP.INITIALIZATION_RANGE
           edit_temp = 1
         endif
       endif
       prior_arg = prior_interactive('F_clump', 'the birth cloud clumpiness factor', $
                                     [0.0d, 0.61d], $
                                     default_priors, $
                                     default_prior_args, $
                                     default_initialization, $
                                     prior_options, prior=prior, $
                                     initialization_range=initialization_range, $
                                     edit=edit_temp)
       config['F_CLUMP','PRIOR'] = strupcase(prior)
       config['F_CLUMP','PRIOR_ARG'] = prior_arg
       config['F_CLUMP','INITIALIZATION_RANGE'] = initialization_range

      ;============ Cos(i) prior ============
       default_priors = 'uniform'
       default_prior_args = [1.0, 0.0, 1.0, 0.5d, 0.2d]
       default_initialization = [0.0d, 1]
       edit_temp = 0
       if edit then begin
         if strupcase(config_edit.ATTEN_CURVE) eq 'DOORE21' then begin
           default_priors = config_edit.COSI.PRIOR
           prior_args = config_edit.COSI.PRIOR_ARG
           case strupcase(default_priors) of
             'FIXED':       default_prior_args = [prior_args, 0.0, 1.0, 0.5d, 0.2d]
             'UNIFORM':     default_prior_args = [1.0d, prior_args, 0.5d, 0.2d]
             'NORMAL':      default_prior_args = [1.0d, prior_args]
             'TABULATED':
           endcase
           default_initialization = config_edit.COSI.INITIALIZATION_RANGE
           edit_temp = 1
         endif
       endif
       prior_arg = prior_interactive('cosi', 'the inclination of the galactic disk in terms of cos(i).', $
                                     [0.0d, 1.0d], $
                                     default_priors, $
                                     default_prior_args, $
                                     default_initialization, $
                                     prior_options, prior=prior, $
                                     initialization_range=initialization_range, $
                                     edit=edit_temp)
       config['COSI','PRIOR'] = strupcase(prior)
       config['COSI','PRIOR_ARG'] = prior_arg
       config['COSI','INITIALIZATION_RANGE'] = initialization_range

      ;============ Bulge-to-disk prior ============
       default_priors = 'uniform'
       default_prior_args = [0.0, 0.0, 1.d3, 10.0, 10.0]
       default_initialization = [0.0d, 5.d1]
       edit_temp = 0
       if edit then begin
         if strupcase(config_edit.ATTEN_CURVE) eq 'DOORE21' then begin
           default_priors = config_edit.B_TO_D.PRIOR
           prior_args = config_edit.B_TO_D.PRIOR_ARG
           case strupcase(default_priors) of
             'FIXED':       default_prior_args = [prior_args, 0.0, 1.d3, 10.0, 10.0]
             'UNIFORM':     default_prior_args = [0.0d, prior_args, 10.0, 10.0]
             'NORMAL':      default_prior_args = [0.0d, prior_args]
             'TABULATED':
           endcase
           default_initialization = config_edit.B_TO_D.INITIALIZATION_RANGE
           edit_temp = 1
         endif
       endif
       prior_arg = prior_interactive('B_TO_D', 'the bulge-to-disk ratio', $
                                     [0.0d, !values.d_infinity], $
                                     default_priors, $
                                     default_prior_args, $
                                     default_initialization, $
                                     prior_options, prior=prior, $
                                     initialization_range=initialization_range, $
                                     edit=edit_temp)
       config['B_TO_D','PRIOR'] = strupcase(prior)
       config['B_TO_D','PRIOR_ARG'] = prior_arg
       config['B_TO_D','INITIALIZATION_RANGE'] = initialization_range

      ;============ Rold0_ages ============
       error = 1
       default_val = dblarr(Nsteps)
       steps_bounds_temp = config['STEPS_BOUNDS', 1:*]
       default_val[where(steps_bounds_temp gt 5e8, /null)] = 1
       default_val = strjoin(string(default_val, f='(I0)'), ', ')
       default = '(Default: '+default_val+')'
       if edit then begin
         if strupcase(config_edit.ATTEN_CURVE) eq 'DOORE21' then begin
           if total(config_edit.STEPS_BOUNDS eq config['STEPS_BOUNDS']) eq Nsteps+1 then begin
             default_val = config_edit.ROLD0_AGES
             default_val = strjoin(string(default_val, f='(I0)'), ', ')
             default = '(Previous configuration: '+default_val+')'
           endif
         endif
       endif
       rold0_message = 'Please specify the binary parameter, rold0, designating each SFH age '+$
                       'bin as part of the young or old population. '+$
                       'A value of 0 for the corresponding age bin considers it to be part '+$
                       'of the young population, and a value of 1 considers it to be part '+$
                       'of the old population. NOTE: Step bins that contain ages < 500 Myr '+$
                       'should be considered part of the young population as they can '+$
                       'contain significant UV emission. Set step bins with ages < 500 Myr '+$
                       'to the old population at your own risk. Please separate values with '+$
                       'a space or comma. Number of expected values: '+$
                       strtrim(string(Nsteps, f='(I0)'), 2)+'. '+default
       print_to_width, '======================='
       while error do begin
         print_to_width, ''
         print_to_width, rold0_message

         temp = ''
         read,temp, prompt='ROLD0_AGES: ',form='(A0)'
         if temp eq '' then temp = default_val
         temp = fix(STRSPLIT(temp, ' ,', /EXTRACT))

         if total(temp eq 0 or temp eq 1) ne Nsteps then begin
           rold0_message = 'Please only specify values of 0 or 1.'
         endif else error = 0
       endwhile
       config['ROLD0_AGES'] = temp
     end
 endcase




;===================================    DUST EMISSION    =====================================================
; Dust emission module options
 ;============ DUST MODEL ============
 dust_emission_options = ['DL07', 'NONE']
 error = 1
 default_val = 'DL07'
 default = '(Default: '+default_val+')'
 if edit then begin
   default_val = strupcase(config_edit.DUST_MODEL)
   default = '(Previous configuration: '+default_val+')'
 endif
 dust_message = 'Please specify the dust emission model to use. '+$
                'Current options: '+strjoin(dust_emission_options, ', ')+'. '+$
                default
 print_to_width, '======================='
 while error do begin
   print_to_width, ''
   print_to_width, dust_message

   temp = ''
   read,temp, prompt='DUST_MODEL: ',form='(A0)'
   if temp eq '' then temp = default_val
   temp = strtrim(temp, 2)

   if total(strupcase(temp) eq strupcase(dust_emission_options)) eq 0 then begin
     dust_message = 'Please specify one of the current options: '+strjoin(dust_emission_options, ', ')+'.'
   endif else error = 0
 endwhile
 config['DUST_MODEL'] = strupcase(temp)


 case strupcase(config['DUST_MODEL']) of
 ; ===================    Draine and Li 2007   =============================================
   'DL07': begin
      ;============ Umin prior ============
       default_priors = 'uniform'
       default_prior_args = [1.0, 0.1d, 25., 1.0, 3.0]
       default_initialization = [0.1d, 10.d]
       edit_temp = 0
       if edit then begin
         if strupcase(config_edit.DUST_MODEL) eq 'DL07' then begin
           default_priors = config_edit.UMIN.PRIOR
           prior_args =     config_edit.UMIN.PRIOR_ARG
           case strupcase(default_priors) of
             'FIXED':       default_prior_args = [prior_args, 0.1d, 25., 1.0, 3.0]
             'UNIFORM':     default_prior_args = [1.0d, prior_args, 10.0, 10.0]
             'NORMAL':      default_prior_args = [1.0d, prior_args]
             'TABULATED':
           endcase
           default_initialization = config_edit.UMIN.INITIALIZATION_RANGE
           edit_temp = 1
         endif
       endif
       prior_arg = prior_interactive('Umin', 'the minimum radiation field intensity of the diffuse ISM '+$
                                     'radiation field heating the dust.', $
                                     [0.10d, 25.d], $
                                     default_priors, $
                                     default_prior_args, $
                                     default_initialization, $
                                     prior_options, prior=prior, $
                                     initialization_range=initialization_range, $
                                     edit=edit_temp)
       config['UMIN','PRIOR'] = strupcase(prior)
       config['UMIN','PRIOR_ARG'] = prior_arg
       config['UMIN','INITIALIZATION_RANGE'] = initialization_range

      ;============ Umax prior ============
       default_priors = 'fixed'
       default_prior_args = [3.d5, 1.d3, 3.d5, 3.d5, 1.d4]
       default_initialization = [1.d5, 3.d5]
       edit_temp = 0
       if edit then begin
         if strupcase(config_edit.DUST_MODEL) eq 'DL07' then begin
           default_priors = config_edit.UMAX.PRIOR
           prior_args =     config_edit.UMAX.PRIOR_ARG
           case strupcase(default_priors) of
             'FIXED':       default_prior_args = [prior_args, 1.d3, 3.d5, 3.d5, 1.d4]
             'UNIFORM':     default_prior_args = [3.d5, prior_args, 3.d5, 1.d4]
             'NORMAL':      default_prior_args = [3.d5, prior_args]
             'TABULATED':
           endcase
           default_initialization = config_edit.UMAX.INITIALIZATION_RANGE
           edit_temp = 1
         endif
       endif
       prior_arg = prior_interactive('Umax', 'the maximum radiation field intensity of the '+$
                                     'power-law distribution of heating starlight intensities', $
                                     [1.d3, 3.d5], $
                                     default_priors, $
                                     default_prior_args, $
                                     default_initialization, $
                                     prior_options, prior=prior, $
                                     initialization_range=initialization_range, $
                                     edit=edit_temp)
       config['UMAX','PRIOR'] = strupcase(prior)
       config['UMAX','PRIOR_ARG'] = prior_arg
       config['UMAX','INITIALIZATION_RANGE'] = initialization_range

      ;============ Alpha prior ============
       default_priors = 'fixed'
       default_prior_args = [2.0, -10.d, 4.0d, 2.0, 1.0]
       default_initialization = [1.0d, 3.0d]
       edit_temp = 0
       if edit then begin
         if strupcase(config_edit.DUST_MODEL) eq 'DL07' then begin
           default_priors = config_edit.ALPHA.PRIOR
           prior_args =     config_edit.ALPHA.PRIOR_ARG
           case strupcase(default_priors) of
             'FIXED':       default_prior_args = [prior_args, -10.d, 4.0d, 2.0, 1.0]
             'UNIFORM':     default_prior_args = [2.0d, prior_args, 2.0, 1.0]
             'NORMAL':      default_prior_args = [2.0d, prior_args]
             'TABULATED':
           endcase
           default_initialization = config_edit.ALPHA.INITIALIZATION_RANGE
           edit_temp = 1
         endif
       endif
       prior_arg = prior_interactive('Alpha', 'the exponent of the power-law distribution of heating starlight '+$
                                     'intensities between Umin and Umax.', $
                                     [-10.d, 4.0d], $
                                     default_priors, $
                                     default_prior_args, $
                                     default_initialization, $
                                     prior_options, prior=prior, $
                                     initialization_range=initialization_range, $
                                     edit=edit_temp)
       config['ALPHA','PRIOR'] = strupcase(prior)
       config['ALPHA','PRIOR_ARG'] = prior_arg
       config['ALPHA','INITIALIZATION_RANGE'] = initialization_range

      ;============ Gamma prior ============
       default_priors = 'uniform'
       default_prior_args = [0.05d, 0.0, 1.0, 0.2d, 0.3d]
       default_initialization = [0.0d, 0.5d]
       edit_temp = 0
       if edit then begin
         if strupcase(config_edit.DUST_MODEL) eq 'DL07' then begin
           default_priors = config_edit.GAMMA.PRIOR
           prior_args =     config_edit.GAMMA.PRIOR_ARG
           case strupcase(default_priors) of
             'FIXED':       default_prior_args = [prior_args, 0.0, 1.0, 0.2d, 0.3d]
             'UNIFORM':     default_prior_args = [0.05d, prior_args, 0.2d, 0.3d]
             'NORMAL':      default_prior_args = [0.05d, prior_args]
             'TABULATED':
           endcase
           default_initialization = config_edit.GAMMA.INITIALIZATION_RANGE
           edit_temp = 1
         endif
       endif
       prior_arg = prior_interactive('Gamma', 'the fraction of the dust mass exposed to the power-law '+$
                                     'distribution of radiation field intensities.', $
                                     [0.0d, 1.0d], $
                                     default_priors, $
                                     default_prior_args, $
                                     default_initialization, $
                                     prior_options, prior=prior, $
                                     initialization_range=initialization_range, $
                                     edit=edit_temp)
       config['GAMMA','PRIOR'] = strupcase(prior)
       config['GAMMA','PRIOR_ARG'] = prior_arg
       config['GAMMA','INITIALIZATION_RANGE'] = initialization_range

      ;============ QPAH prior ============
       default_priors = 'uniform'
       default_prior_args = [1d-2, 4.7d-3, 4.58d-2, 1d-2, 1d-2]
       default_initialization = [4.7d-3, 4.58d-2]
       edit_temp = 0
       if edit then begin
         if strupcase(config_edit.DUST_MODEL) eq 'DL07' then begin
           default_priors = config_edit.QPAH.PRIOR
           prior_args =     config_edit.QPAH.PRIOR_ARG
           case strupcase(default_priors) of
             'FIXED':       default_prior_args = [prior_args, 4.7d-3, 4.58d-2, 1d-2, 1d-2]
             'UNIFORM':     default_prior_args = [1d-2, prior_args, 1d-2, 1d-2]
             'NORMAL':      default_prior_args = [1d-2, prior_args]
             'TABULATED':
           endcase
           default_initialization = config_edit.QPAH.INITIALIZATION_RANGE
           edit_temp = 1
         endif
       endif
       prior_arg = prior_interactive('qPAH', 'the fraction of the total grain mass corresponding to '+$
                                     'PAHs containing less than 1000 carbon atoms (PAH index)', $
                                     [4.7d-3, 4.58d-2], $
                                     default_priors, $
                                     default_prior_args, $
                                     default_initialization, $
                                     prior_options, prior=prior, $
                                     initialization_range=initialization_range, $
                                     edit=edit_temp)
       config['QPAH','PRIOR'] = strupcase(prior)
       config['QPAH','PRIOR_ARG'] = prior_arg
       config['QPAH','INITIALIZATION_RANGE'] = initialization_range

       if ~(config['ENERGY_BALANCE']) then begin
        ;============ LTIR prior ============
         default_priors = 'uniform'
         default_prior_args = [1.d10, 0.0, 1.d11, 1.d10, 3.d10]
         default_initialization = [1.d7, 1.d11]
         edit_temp = 0
         if edit then begin
           if strupcase(config_edit.DUST_MODEL) eq 'DL07' and ~(config_edit.ENERGY_BALANCE) then begin
             default_priors = config_edit.LTIR.PRIOR
             prior_args =     config_edit.LTIR.PRIOR_ARG
             case strupcase(default_priors) of
               'FIXED':       default_prior_args = [prior_args,  0.0, 1d11, 1d10, 3d10]
               'UNIFORM':     default_prior_args = [1d10, prior_args, 1d10, 3d10]
               'NORMAL':      default_prior_args = [1d10, prior_args]
             'TABULATED':
             endcase
             default_initialization = config_edit.LTIR.INITIALIZATION_RANGE
             edit_temp = 1
           endif
         endif
         prior_arg = prior_interactive('LTIR', 'the total integrated IR luminosity in Lsun', $
                                       [0.0d, !values.d_infinity], $
                                       default_priors, $
                                       default_prior_args, $
                                       default_initialization, $
                                       prior_options, prior=prior, $
                                       initialization_range=initialization_range, $
                                       edit=edit_temp)
         config['LTIR','PRIOR'] = strupcase(prior)
         config['LTIR','PRIOR_ARG'] = prior_arg
         config['LTIR','INITIALIZATION_RANGE'] = initialization_range
       endif
     end

    'NONE':
 endcase




;===================================    X-RAY EMISSION    ====================================================
; X-ray emission module options
 ;============ Xray Emission ============
 error = 1
 default_val = 'NO'
 default = '(Default: '+strlowcase(default_val)+')'
 if edit then begin
   if config_edit.XRAY_EMISSION eq 1 then default_val = 'YES' else default_val = 'NO'
   default = '(Previous configuration: '+strlowcase(default_val)+')'
 endif
 ; Check if SSP is not NONE, if it is NONE, no xray emission can be used.
 if strupcase(config['SSP']) eq 'NONE' then begin
   xray_message = "You can not have X-ray emission if not using a stellar population model (i.e., SSP='NONE')."+$
                  ' Please indicate no, for no X-ray emission. (Required value: no)'
   default_val = 'NO'
 endif else begin
   xray_message = 'Would you like to include an X-ray emission model? This always includes '+$
                  'stellar X-ray emission, but can optionally include AGN X-ray emission. '+$
                  'Please indicate yes or no. '+default
 endelse
 print_to_width, '======================='
 while error do begin
   print_to_width, ''
   print_to_width, xray_message

   temp = ''
   read,temp, prompt='XRAY_EMISSION: ',form='(A0)'
   if temp eq '' then temp = default_val
   temp = strtrim(temp, 2)

   if strupcase(config['SSP']) eq 'NONE' then begin
     if strupcase(temp) ne 'NO' then begin
       xray_message = 'Please type no.'
     endif else error = 0
   endif else begin
     if ~(strupcase(temp) eq 'NO' or strupcase(temp) eq 'YES') then begin
       xray_message = 'Please type yes or no.'
     endif else error = 0
   endelse
 endwhile
 if strupcase(temp) eq 'YES' then temp = 1 else temp = 0
 config['XRAY_EMISSION'] = temp


 if (config['XRAY_EMISSION']) then begin
  ;============ Xray units ============
   xray_unit_options = ['COUNTS', 'FLUX']
   error = 1
   default_val = 'COUNTS'
   default = '(Default: '+default_val+')'
   if edit then begin
     if config_edit.XRAY_EMISSION then begin
       default_val = strupcase(config_edit.XRAY_UNIT)
       default = '(Previous configuration: '+default_val+')'
     endif
   endif
   xray_unit_message = 'Please specify the type of X-ray data to use for fitting. '+$
                       'Current options: '+strjoin(xray_unit_options, ', ')+$
                       '. '+default
   print_to_width, '======================='
   while error do begin
     print_to_width, ''
     print_to_width, xray_unit_message

     temp = ''
     read,temp, prompt='XRAY_UNIT: ',form='(A0)'
     if temp eq '' then temp = default_val
     temp = strtrim(temp, 2)

     if total(strupcase(temp) eq strupcase(xray_unit_options)) eq 0 then begin
       xray_unit_message = 'Please specify one of the current options: '+strjoin(xray_unit_options, ', ')+'.'
     endif else error = 0
   endwhile
   config['XRAY_UNIT'] = strupcase(temp)

  ;============ Xray uncertainties ============
   if config['XRAY_UNIT'] eq 'COUNTS' then begin
     xray_unc_options = ['SQRT', 'GEHRELS', 'USER']
     error = 1
     default_val = 'GEHRELS'
     default = '(Default: '+default_val+')'
     if edit then begin
       if config_edit.XRAY_EMISSION then begin
         default_val = strupcase(config_edit.XRAY_UNC)
         default = '(Previous configuration: '+default_val+')'
       endif
     endif
     xray_unc_message = 'Please specify the uncertainties to assume for the X-ray counts. '+$
                         'Current options: '+strjoin(xray_unc_options, ', ')+$
                         '. '+default
     print_to_width, '======================='
     while error do begin
       print_to_width, ''
       print_to_width, xray_unc_message

       temp = ''
       read,temp, prompt='XRAY_UNC: ',form='(A0)'
       if temp eq '' then temp = default_val
       temp = strtrim(temp, 2)

       if total(strupcase(temp) eq strupcase(xray_unc_options)) eq 0 then begin
         xray_unc_message = 'Please specify one of the current options: '+strjoin(xray_unc_options, ', ')+'.'
       endif else error = 0
     endwhile
     config['XRAY_UNC'] = strupcase(temp)
   endif

  ;============ Xray Absorption ============
   xray_abs_options = ['TBABS-WILM', 'TBABS-ANGR', 'ATTEN']
   error = 1
   default_val = 'TBABS-WILM'
   default = '(Default: '+default_val+')'
   if edit then begin
     if config_edit.XRAY_EMISSION then begin
       default_val = strupcase(config_edit.XRAY_ABS_MODEL)
       default = '(Previous configuration: '+default_val+')'
     endif
   endif
   xray_abs_message = 'Please specify the X-ray absorption model to apply to the X-ray emission. '+$
                        'Current options: '+strjoin(xray_abs_options, ', ')+$
                        '. '+default
   print_to_width, '======================='
   while error do begin
     print_to_width, ''
     print_to_width, xray_abs_message

     temp = ''
     read,temp, prompt='XRAY_ABS_MODEL: ',form='(A0)'
     if temp eq '' then temp = default_val
     temp = strtrim(temp, 2)

     if total(strupcase(temp) eq strupcase(xray_abs_options)) eq 0 then begin
       xray_abs_message = 'Please specify one of the current options: '+strjoin(xray_abs_options, ', ')+'.'
     endif else error = 0
   endwhile
   config['XRAY_ABS_MODEL'] = strupcase(temp)


  ;============ nH prior ============
   default_priors = 'uniform'
   default_prior_args = [1.d2, 1.d-4, 1.d5, 1d2, 1d3]
   default_initialization = [1.d-1, 1.d2]
   edit_temp = 0
   if edit then begin
     if config_edit.XRAY_EMISSION then begin
       default_priors = config_edit.NH.PRIOR
       prior_args =     config_edit.NH.PRIOR_ARG
       case strupcase(default_priors) of
         'FIXED':       default_prior_args = [prior_args, 1.d-4, 1.d5, 1d2, 1d3]
         'UNIFORM':     default_prior_args = [1d2, prior_args, 1d2, 1d3]
         'NORMAL':      default_prior_args = [1d2, prior_args]
         'TABULATED':
       endcase
       default_initialization = config_edit.NH.INITIALIZATION_RANGE
       edit_temp = 1
     endif
   endif
   prior_arg = prior_interactive('nH', 'the intrinsic HI column density along the line of sight in 1e20 cm-2', $
                                 [1.d-4, 1.d5], $
                                 default_priors, $
                                 default_prior_args, $
                                 default_initialization, $
                                 prior_options, prior=prior, $
                                 initialization_range=initialization_range, $
                                 edit=edit_temp)
   config['NH','PRIOR'] = strupcase(prior)
   config['NH','PRIOR_ARG'] = prior_arg
   config['NH','INITIALIZATION_RANGE'] = initialization_range


  ;============ Xray AGN model ============
   xray_agn_options = ['QSOSED', 'PLAW', 'NONE']
   error = 1
   default_val = 'QSOSED'
   default = '(Default: '+default_val+')'
   if edit then begin
     if config_edit.XRAY_EMISSION then begin
       default_val = strupcase(config_edit.XRAY_AGN_MODEL)
       default = '(Previous configuration: '+default_val+')'
     endif
   endif
   xray_agn_message = 'Please specify the AGN X-ray emission model to use. '+$
                      'Current options: '+strjoin(xray_agn_options, ', ')+$
                      '. '+default
   print_to_width, '======================='
   while error do begin
     print_to_width, ''
     print_to_width, xray_agn_message

     temp = ''
     read,temp, prompt='XRAY_AGN_MODEL: ',form='(A0)'
     if temp eq '' then temp = default_val
     temp = strtrim(temp, 2)

     if total(strupcase(temp) eq strupcase(xray_agn_options)) eq 0 then begin
       xray_agn_message = 'Please specify one of the current options: '+strjoin(xray_agn_options, ', ')+'.'
     endif else error = 0
   endwhile
   config['XRAY_AGN_MODEL'] = strupcase(temp)


   case strupcase(config['XRAY_AGN_MODEL']) of
    ;================= QSOSED PARAMETERS =========================================
     'QSOSED': begin
       ;============ AGN SMBH prior ============
         default_priors = 'uniform'
         default_prior_args = [1.d8, 1.d5, 1.d10, 1.d8, 1.d9]
         default_initialization = [1.d6, 1.d9]
         edit_temp = 0
         if edit then begin
           if config_edit.XRAY_EMISSION then begin
             if strupcase(config_edit.XRAY_AGN_MODEL) eq 'QSOSED' then begin
               default_priors = config_edit.AGN_MASS.PRIOR
               prior_args =     config_edit.AGN_MASS.PRIOR_ARG
               case strupcase(default_priors) of
                 'FIXED':       default_prior_args = [prior_args, 1.d5, 1.d10, 1.d8, 1.d9]
                 'UNIFORM':     default_prior_args = [1.d8, prior_args, 1.d8, 1.d9]
                 'NORMAL':      default_prior_args = [1.d8, prior_args]
                 'TABULATED':
               endcase
               default_initialization = config_edit.AGN_MASS.INITIALIZATION_RANGE
               edit_temp = 1
             endif
           endif
         endif
         prior_arg = prior_interactive('AGN_mass', 'the SMBH mass in Msun', $
                                       [1.d5, 1.d10], $
                                       default_priors, $
                                       default_prior_args, $
                                       default_initialization, $
                                       prior_options, prior=prior, $
                                       initialization_range=initialization_range, $
                                       edit=edit_temp)
         config['AGN_MASS','PRIOR'] = strupcase(prior)
         config['AGN_MASS','PRIOR_ARG'] = prior_arg
         config['AGN_MASS','INITIALIZATION_RANGE'] = initialization_range

       ;============ AGN mdot prior ============
         default_priors = 'uniform'
         default_prior_args = [-0.2d, -1.5d, 0.3d, -0.2d, 0.3d]
         default_initialization = [-1.d, 0.0d]
         edit_temp = 0
         if edit then begin
           if config_edit.XRAY_EMISSION then begin
             if strupcase(config_edit.XRAY_AGN_MODEL) eq 'QSOSED' then begin
               default_priors = config_edit.AGN_LOGMDOT.PRIOR
               prior_args =     config_edit.AGN_LOGMDOT.PRIOR_ARG
               case strupcase(default_priors) of
                 'FIXED':       default_prior_args = [prior_args, -1.5d, 0.3d, -0.2d, 0.3d]
                 'UNIFORM':     default_prior_args = [-0.2d, prior_args, -0.2d, 0.3d]
                 'NORMAL':      default_prior_args = [-0.2d, prior_args]
                 'TABULATED':
               endcase
               default_initialization = config_edit.AGN_LOGMDOT.INITIALIZATION_RANGE
               edit_temp = 1
             endif
           endif
         endif
         prior_arg = prior_interactive('AGN_log_mdot', 'the log10 of the SMBH accretion rate '+$
                                       'normalized by the Eddington rate', $
                                       [-1.5d, 0.3d], $
                                       default_priors, $
                                       default_prior_args, $
                                       default_initialization, $
                                       prior_options, prior=prior, $
                                       initialization_range=initialization_range, $
                                       edit=edit_temp)
         config['AGN_LOGMDOT','PRIOR'] = strupcase(prior)
         config['AGN_LOGMDOT','PRIOR_ARG'] = prior_arg
         config['AGN_LOGMDOT','INITIALIZATION_RANGE'] = initialization_range
       end

     'PLAW':

     'NONE':

   endcase
 endif




;====================================    AGN EMISSION    =====================================================
; AGN emission module options
 ;============ AGN MODEL ============
 agn_emission_options = ['SKIRTOR', 'NONE']
 error = 1
 default_val = 'NONE'
 default = '(Default: '+default_val+')'
 if edit then begin
   default_val = strupcase(config_edit.AGN_MODEL)
   default = '(Previous configuration: '+default_val+')'
 endif
 ; Check if DOORE21 attenuation is used, since it is not compatible
 if strupcase(config['ATTEN_CURVE']) eq 'DOORE21' then begin
   agn_emission_options = ['none']
   agn_message = 'The AGN emission model is currently incompatible with the DOORE21 '+$
                 'attenuation curve. Please specify: '+strjoin(agn_emission_options, ', ')+$
                 '. '
   default_val = 'NONE'
 endif else begin
   agn_message = 'Please specify the UV-to-IR AGN emission model to use. '+$
                 'Current options: '+strjoin(agn_emission_options, ', ')+$
                 '. '+default
 endelse
 print_to_width, '======================='
 while error do begin
   print_to_width, ''
   print_to_width, agn_message

   temp = ''
   read,temp, prompt='AGN_MODEL: ',form='(A0)'
   if temp eq '' then temp = default_val
   temp = strtrim(temp, 2)

   if total(strupcase(temp) eq strupcase(agn_emission_options)) eq 0 then begin
     agn_message = 'Please specify one of the current options: '+strjoin(agn_emission_options, ', ')+'.'
   endif else error = 0
 endwhile
 config['AGN_MODEL'] = strupcase(temp)


 case strupcase(config['AGN_MODEL']) of
 ; ===================    AGN MODEL PARAMETERS   =============================================
   'SKIRTOR': begin
     ;============ Tau97 prior ============
       default_priors = 'fixed'
       default_prior_args = [7.0d, 3.0, 11.0, 5.0, 1.0]
       default_initialization = [3.0d, 11.0d]
       edit_temp = 0
       if edit then begin
         if strupcase(config_edit.AGN_MODEL) eq 'SKIRTOR' then begin
           default_priors = config_edit.TAU97.PRIOR
           prior_args =     config_edit.TAU97.PRIOR_ARG
           case strupcase(default_priors) of
             'FIXED':       default_prior_args = [prior_args, 3.0d, 11.0, 5.0, 1.0]
             'UNIFORM':     default_prior_args = [7.0d, prior_args, 5.0, 1.0]
             'NORMAL':      default_prior_args = [7.0d, prior_args]
             'TABULATED':
           endcase
           default_initialization = config_edit.TAU97.INITIALIZATION_RANGE
           edit_temp = 1
         endif
       endif
       prior_arg = prior_interactive('tau97', 'the edge-on optical depth of AGN dust torus at 9.7 um', $
                                     [3.0d, 11.d0], $
                                     default_priors, $
                                     default_prior_args, $
                                     default_initialization, $
                                     prior_options, prior=prior, $
                                     initialization_range=initialization_range, $
                                     edit=edit_temp)
       config['TAU97','PRIOR'] = strupcase(prior)
       config['TAU97','PRIOR_ARG'] = prior_arg
       config['TAU97','INITIALIZATION_RANGE'] = initialization_range

     ;============ AGN cos(i) prior ============
       default_priors = 'uniform'
       default_prior_args = [1.0, 0.0, 1.0, 0.5d, 0.2d]
       default_initialization = [0.0d, 1.0d]
       edit_temp = 0
       if edit then begin
         if strupcase(config_edit.AGN_MODEL) eq 'SKIRTOR' then begin
           default_priors = config_edit.AGN_COSI.PRIOR
           prior_args =     config_edit.AGN_COSI.PRIOR_ARG
           case strupcase(default_priors) of
             'FIXED':       default_prior_args = [prior_args, 0.0, 1.0, 0.5d, 0.2d]
             'UNIFORM':     default_prior_args = [1.0, prior_args, 0.5d, 0.2d]
             'NORMAL':      default_prior_args = [1.0, prior_args]
             'TABULATED':
           endcase
           default_initialization = config_edit.AGN_COSI.INITIALIZATION_RANGE
           edit_temp = 1
         endif
       endif
       prior_arg = prior_interactive('AGN_COSI', 'inclination of the AGN disk in terms of cos(i)', $
                                     [0.0d, 1.0d], $
                                     default_priors, $
                                     default_prior_args, $
                                     default_initialization, $
                                     prior_options, prior=prior, $
                                     initialization_range=initialization_range, $
                                     edit=edit_temp)
       config['AGN_COSI','PRIOR'] = strupcase(prior)
       config['AGN_COSI','PRIOR_ARG'] = prior_arg
       config['AGN_COSI','INITIALIZATION_RANGE'] = initialization_range


       ; Complicated logic due to 'XRAY_AGN_MODEL' key not existing if 'XRAY_EMISSION' tag is 0
       if ~(config['XRAY_EMISSION']) then begin
         default_priors = 'uniform'
         default_prior_args = [10.d, 0.0, 20., 10., 3.]
         default_initialization = [6.d, 12.d]
         edit_temp = 0
         if edit then begin
           if ~(config_edit.XRAY_EMISSION) and strupcase(config_edit.AGN_MODEL) eq 'SKIRTOR' then begin
             default_priors = config_edit.LOG_L_AGN.PRIOR
             prior_args =     config_edit.LOG_L_AGN.PRIOR_ARG
             case strupcase(default_priors) of
               'FIXED':       default_prior_args = [prior_args, 0.0d, 20., 10, 3.]
               'UNIFORM':     default_prior_args = [10.d, prior_args, 10, 3.]
               'NORMAL':      default_prior_args = [10.d, prior_args]
               'TABULATED':
             endcase
             default_initialization = config_edit.LOG_L_AGN.INITIALIZATION_RANGE
             edit_temp = 1
           endif
         endif
         prior_arg = prior_interactive('LOG_L_AGN', 'the total integrated luminosity of AGN model in log10(Lsun)', $
                                       [0.0d, 20.d0], $
                                       default_priors, $
                                       default_prior_args, $
                                       default_initialization, $
                                       prior_options, prior=prior, $
                                       initialization_range=initialization_range, $
                                       edit=edit_temp)
         config['LOG_L_AGN','PRIOR'] = strupcase(prior)
         config['LOG_L_AGN','PRIOR_ARG'] = prior_arg
         config['LOG_L_AGN','INITIALIZATION_RANGE'] = initialization_range
       endif else begin
         default_priors = 'uniform'
         default_prior_args = [10.d, 0.0, 20., 10., 3.]
         default_initialization = [6.d, 12.d]
         edit_temp = 0
         if strupcase(config['XRAY_AGN_MODEL']) ne 'QSOSED' then begin
           if edit then begin
             if config_edit.XRAY_EMISSION and strupcase(config_edit.AGN_MODEL) eq 'SKIRTOR' then begin
               if strupcase(config_edit.XRAY_AGN_MODEL) ne 'QSOSED'  then begin
                 default_priors = config_edit.LOG_L_AGN.PRIOR
                 prior_args =     config_edit.LOG_L_AGN.PRIOR_ARG
                 case strupcase(default_priors) of
                   'FIXED':       default_prior_args = [prior_args, 0.0d, 20., 10, 3.]
                   'UNIFORM':     default_prior_args = [10.d, prior_args, 10, 3.]
                   'NORMAL':      default_prior_args = [10.d, prior_args]
                   'TABULATED':
                 endcase
                 default_initialization = config_edit.LOG_L_AGN.INITIALIZATION_RANGE
                 edit_temp = 1
               endif
             endif
           endif
           prior_arg = prior_interactive('LOG_L_AGN', 'the total integrated luminosity of AGN model in log10(Lsun)', $
                                         [0.0d, 20.d0], $
                                         default_priors, $
                                         default_prior_args, $
                                         default_initialization, $
                                         prior_options, prior=prior, $
                                         initialization_range=initialization_range, $
                                         edit=edit_temp)
           config['LOG_L_AGN','PRIOR'] = strupcase(prior)
           config['LOG_L_AGN','PRIOR_ARG'] = prior_arg
           config['LOG_L_AGN','INITIALIZATION_RANGE'] = initialization_range
         endif
       endelse
     end

    'NONE':
 endcase




;=================================    FITTING ALGORITHM    ===================================================
; Fitting algorithm module options
 ;============ Fitting algorithm ============
 ;fit_algorithm_options = ['GRID', 'INVERSION', 'MCMC-ADAPTIVE', 'MCMC-AFFINE', 'MPFIT']
 fit_algorithm_options = ['MCMC-ADAPTIVE', 'MCMC-AFFINE', 'MPFIT']
 error = 1
 if edit then begin
   default_val = strupcase(config_edit.METHOD)
   default = '(Previous configuration: '+default_val+')'
 endif else begin
   default_val = 'MCMC-AFFINE'
   default = '(Default: '+default_val+')'
 endelse
 fit_message = 'Please specify the fitting algorithm used to fit the SED(s). '+$
               'Current options: '+strjoin(fit_algorithm_options, ', ')+'. '+$
               default
 print_to_width, '======================='
 while error do begin
   print_to_width, ''
   print_to_width, fit_message

   temp = ''
   read,temp, prompt='FITTING_ALGORITHM: ',form='(A0)'
   if temp eq '' then temp = default_val
   temp = strtrim(temp, 2)

   if total(strupcase(temp) eq strupcase(fit_algorithm_options)) eq 0 then begin
     fit_message = 'Please specify one of the current options: '+strjoin(fit_algorithm_options, ', ')+'.'
   endif else error = 0
 endwhile
 config['METHOD'] = strupcase(temp)


 case strupcase(config['METHOD']) of
 ;============ Grid algorithm ============
   ;'GRID': begin
   ;;============ Grid elements ============
   ;    error = 1
   ;    default_val = '10'
   ;    default = '(Default: '+default_val+')'
   ;    if edit then begin
   ;      if strupcase(config_edit.METHOD) eq 'GRID' then begin
   ;        default_val = strtrim(string(config_edit.GRID_ELEMENTS, f='(I0)'), 2)
   ;        default = '(Previous configuration: '+default_val+')'
   ;      endif
   ;    endif
   ;    grid_message = 'Please specify the number of elements to evenly divide each '+$
   ;                   'non-fixed parameter to generate the grid based solution. '+$
   ;                   default
   ;    print_to_width, '======================='
   ;    while error do begin
   ;      print_to_width, ''
   ;      print_to_width, grid_message
   ;
   ;      temp = ''
   ;      read,temp, prompt='GRID_ELEMENTS: ',form='(A0)'
   ;      if temp eq '' then temp = default_val
   ;      temp = long(temp)
   ;
   ;      if temp le 0 then begin
   ;        grid_message = 'Please specify a positive value.'
   ;      endif else error = 0
   ;    endwhile
   ;    config['GRID_ELEMENTS'] = temp
   ;  end


 ;============ Inversion algorithm ============
   ;'INVERSION': begin
   ;;============ Negative SFR ============
   ;    error = 1
   ;    default_val = 'NO'
   ;    default = '(Default: '+strlowcase(default_val)+')'
   ;    if edit then begin
   ;      if strupcase(config_edit.METHOD) eq 'INVERSION' then begin
   ;        if config_edit.ALLOW_NEG_SFH eq 1 then default_val = 'YES' else default_val = 'NO'
   ;        default = '(Previous configuration: '+strlowcase(default_val)+')'
   ;      endif
   ;    endif
   ;    neg_sfh_message = 'Would you like to allow for negative SFH bin values when '+$
   ;                      'computing the SFH. Otherwise, zero or positive values will be '+$
   ;                      'enforced. Please indicate yes or no. '+default
   ;    print_to_width, '======================='
   ;    while error do begin
   ;      print_to_width, ''
   ;      print_to_width, neg_sfh_message
   ;
   ;      temp = ''
   ;      read,temp, prompt='ALLOW_NEG_SFH: ',form='(A0)'
   ;      if temp eq '' then temp = default_val
   ;      temp = strtrim(temp, 2)
   ;
   ;      if ~(strupcase(temp) eq 'YES' or strupcase(temp) eq 'NO') then begin
   ;        neg_sfh_message = 'Please type yes or no.'
   ;      endif else error = 0
   ;    endwhile
   ;    if strupcase(temp) eq 'YES' then temp=1 else temp=0
   ;    config['ALLOW_NEG_SFH']=temp
   ;
   ;;============ Inversion Grid elements ============
   ;    error = 1
   ;    default_val = '10'
   ;    default = '(Default: '+default_val+')'
   ;    if edit then begin
   ;      if strupcase(config_edit.METHOD) eq 'INVERSION' then begin
   ;        default_val = strtrim(string(config_edit.INVERSION_GRID_ELEMENTS, f='(I0)'), 2)
   ;        default = '(Previous configuration: '+default_val+')'
   ;      endif
   ;    endif
   ;    grid_message = 'Please specify the number of elements to evenly divide each '+$
   ;                   'non-fixed, non-SFH age bin parameter to generate the matrix '+$
   ;                   'inversion based solution. '+default
   ;    print_to_width, '======================='
   ;    while error do begin
   ;      print_to_width, ''
   ;      print_to_width, grid_message
   ;
   ;      temp = ''
   ;      read,temp, prompt='INVERSION_GRID_ELEMENTS: ',form='(A0)'
   ;      if temp eq '' then temp = default_val
   ;      temp = long(temp)
   ;
   ;      if temp le 0 then begin
   ;        grid_message = 'Please specify a positive value.'
   ;      endif else error = 0
   ;    endwhile
   ;    config['INVERSION_GRID_ELEMENTS']=temp
   ;  end


 ;============ Adaptive MCMC algorithm ============
   'MCMC-ADAPTIVE': begin
    ;============ Number of trials ============
       error = 1
       default_val = '1.5e5'
       default = '(Default: '+default_val+')'
       if edit then begin
         if strupcase(config_edit.METHOD) eq 'MCMC-ADAPTIVE' then begin
           default_val = strtrim(string(config_edit.NTRIALS, f='(I0)'), 2)
           default = '(Previous configuration: '+default_val+')'
         endif
       endif
       ntrials_message = 'Please specify the number of MCMC trials to run. '+$
                         default
       print_to_width, '======================='
       while error do begin
         print_to_width, ''
         print_to_width, ntrials_message

         temp = ''
         read,temp, prompt='NTRIALS: ',form='(A0)'
         if temp eq '' then temp = default_val
         temp = long(double(temp))

         if temp le 0 then begin
           ntrials_message = 'Please specify a positive value.'
         endif else error = 0
       endwhile
       config['NTRIALS'] = temp

   ;============ Number of parallel chains ============
       error = 1
       default_val = '1'
       default = '(Default: '+default_val+')'
       if edit then begin
         if strupcase(config_edit.METHOD) eq 'MCMC-ADAPTIVE' then begin
           default_val = strtrim(string(config_edit.NPARALLEL, f='(I0)'), 2)
           default = '(Previous configuration: '+default_val+')'
         endif
       endif
       npar_message = 'Please specify the number of parallel chains to run for each SED. '+$
                      default
       print_to_width, '======================='
       while error do begin
         print_to_width, ''
         print_to_width, npar_message

         temp = ''
         read,temp, prompt='NPARALLEL: ',form='(A0)'
         if temp eq '' then temp = default_val
         temp = long(double(temp))

         if temp le 0 then begin
           npar_message = 'Please specify a positive value.'
         endif else error = 0
       endwhile
       config['NPARALLEL'] = temp

   ;============ Autocorrelation step ============
       error = 1
       default_val = '5'
       default = '(Default: '+default_val+')'
       if edit then begin
         if strupcase(config_edit.METHOD) eq 'MCMC-ADAPTIVE' then begin
           default_val = strtrim(string(config_edit.C_STEP, f='(I0)'), 2)
           default = '(Previous configuration: '+default_val+')'
         endif
       endif
       cstep_message = 'Please specify the step size defining how many trials of the chain are used to calculate '+$
                      'the autocorrelation time (tau), where we integrate tau to the smallest index m such that '+$
                      'm > C_step * tau. We recommend not changing this value from the default. '+$
                      default
       print_to_width, '======================='
       while error do begin
         print_to_width, ''
         print_to_width, cstep_message

         temp = ''
         read,temp, prompt='C_STEP: ',form='(A0)'
         if temp eq '' then temp = default_val
         temp = double(temp)

         if temp le 0 then begin
           cstep_message = 'Please specify a positive value.'
         endif else error = 0
       endwhile
       config['C_STEP'] = temp

   ;============ Tolerance of chain length to autocorrelation time ============
       error = 1
       default_val = '50'
       default = '(Default: '+default_val+')'
       if edit then begin
         if strupcase(config_edit.METHOD) eq 'MCMC-ADAPTIVE' then begin
           default_val = strtrim(string(config_edit.TOLERANCE, f='(I0)'), 2)
           default = '(Previous configuration: '+default_val+')'
         endif
       endif
       toler_message = 'Please specify the tolerance to use that defines how many taus '+$
                       'the length of the chain should be for us to believe the estimated value of tau. '+$
                       default
       print_to_width, '======================='
       while error do begin
         print_to_width, ''
         print_to_width, toler_message

         temp = ''
         read,temp, prompt='TOLERANCE: ',form='(A0)'
         if temp eq '' then temp = default_val
         temp = double(temp)

         if temp le 0 then begin
           toler_message = 'Please specify a positive value.'
         endif else error = 0
       endwhile
       config['TOLERANCE'] = temp

    ;============ Adaptiveness factor ============
       error = 1
       default_val = '0.35'
       default = '(Default: '+default_val+')'
       if edit then begin
         if strupcase(config_edit.METHOD) eq 'MCMC-ADAPTIVE' then begin
           default_val = strtrim(string(config_edit.BETA_EXPONENT, f='(F0.3)'), 2)
           default = '(Previous configuration: '+default_val+')'
         endif
       endif
       beta_message = 'Please specify the factor controlling how fast the adaptiveness of the '+$
                      'algorithm vanishes. Larger values stop the adaptiveness '+$
                      'in fewer trials. '+default
       print_to_width, '======================='
       while error do begin
         print_to_width, ''
         print_to_width, beta_message

         temp = 0.d0
         read,temp, prompt='BETA_EXPONENT: ',form='(A0)'
         if temp eq '' then temp = default_val
         temp = double(temp)

         if temp lt 0 then begin
           beta_message = 'Please specify a positive value.'
         endif else error = 0
       endwhile
       config['BETA_EXPONENT'] = temp
     end


 ;============ Affine invariant MCMC algorithm ============
   'MCMC-AFFINE': begin
   ;============ Number of trials ============
       error = 1
       default_val = '2e4'
       default = '(Default: '+default_val+')'
       if edit then begin
         if strupcase(config_edit.METHOD) eq 'MCMC-AFFINE' then begin
           default_val = strtrim(string(config_edit.NTRIALS, f='(I0)'), 2)
           default = '(Previous configuration: '+default_val+')'
         endif
       endif
       ntrials_message = 'Please specify the number of MCMC trials to run. '+default
       print_to_width, '======================='
       while error do begin
         print_to_width, ''
         print_to_width, ntrials_message

         temp = ''
         read,temp, prompt='NTRIALS: ',form='(A0)'
         if temp eq '' then temp = default_val
         temp = long(double(temp))

         if temp le 0 then begin
           ntrials_message = 'Please specify a positive value.'
         endif else error = 0
       endwhile
       config['NTRIALS'] = temp

   ;============ Number of walkers ============
       error = 1
       default_val = '75'
       default = '(Default: '+default_val+')'
       if edit then begin
         if strupcase(config_edit.METHOD) eq 'MCMC-AFFINE' then begin
           default_val = strtrim(string(config_edit.NPARALLEL, f='(I0)'), 2)
           default = '(Previous configuration: '+default_val+')'
         endif
       endif
       npar_message = 'Please specify the number of parallel walkers to run for each SED. '+$
                      'This value must be greater than the number of free parameters plus one '+$
                      '(ideally at least twice the number of free parameters) or optimal '+$
                      'sampling. '+default
       print_to_width, '======================='
       while error do begin
         print_to_width, ''
         print_to_width, npar_message

         temp = ''
         read,temp, prompt='NPARALLEL: ',form='(A0)'
         if temp eq '' then temp = default_val
         temp = long(double(temp))

         if temp le 0 then begin
           npar_message = 'Please specify a positive value.'
         endif else error = 0
       endwhile
       config['NPARALLEL'] = temp

   ;============ Autocorrelation step ============
       error = 1
       default_val = '5'
       default = '(Default: '+default_val+')'
       if edit then begin
         if strupcase(config_edit.METHOD) eq 'MCMC-AFFINE' then begin
           default_val = strtrim(string(config_edit.C_STEP, f='(I0)'), 2)
           default = '(Previous configuration: '+default_val+')'
         endif
       endif
       cstep_message = 'Please specify the step size defining how many trials of the chain are used to calculate '+$
                      'the autocorrelation time (tau), where we integrate tau to the smallest index m such that '+$
                      'm > C_step * tau. We recommend not changing this value from the default. '+$
                      default
       print_to_width, '======================='
       while error do begin
         print_to_width, ''
         print_to_width, cstep_message

         temp = ''
         read,temp, prompt='C_STEP: ',form='(A0)'
         if temp eq '' then temp = default_val
         temp = double(temp)

         if temp le 0 then begin
           cstep_message = 'Please specify a positive value.'
         endif else error = 0
       endwhile
       config['C_STEP'] = temp

   ;============ Tolerance of chain length to autocorrelation time ============
       error = 1
       default_val = '50'
       default = '(Default: '+default_val+')'
       if edit then begin
         if strupcase(config_edit.METHOD) eq 'MCMC-AFFINE' then begin
           default_val = strtrim(string(config_edit.TOLERANCE, f='(I0)'), 2)
           default = '(Previous configuration: '+default_val+')'
         endif
       endif
       toler_message = 'Please specify the tolerance to use that defines how many taus '+$
                       'the length of the chain should be for us to believe the estimated value of tau. '+$
                       default
       print_to_width, '======================='
       while error do begin
         print_to_width, ''
         print_to_width, toler_message

         temp = ''
         read,temp, prompt='TOLERANCE: ',form='(A0)'
         if temp eq '' then temp = default_val
         temp = double(temp)

         if temp le 0 then begin
           toler_message = 'Please specify a positive value.'
         endif else error = 0
       endwhile
       config['TOLERANCE'] = temp

    ;============ Move scaling constant ============
       error = 1
       default_val = '2.0'
       default = '(Default: '+default_val+')'
       if edit then begin
         if strupcase(config_edit.METHOD) eq 'MCMC-AFFINE' then begin
           default_val = strtrim(string(config_edit.AFFINE_A, f='(I0)'), 2)
           default = '(Previous configuration: '+default_val+')'
         endif
       endif
       a_message = 'Please specify the move scaling constant defining the maximum and '+$
                   'minimum step size of the affine-invariant stretch move. '+default
       print_to_width, '======================='
       while error do begin
         print_to_width, ''
         print_to_width, a_message

         temp = ''
         read,temp, prompt='AFFINE_A: ',form='(A0)'
         if temp eq '' then temp = default_val
         temp = double(temp)

         if temp lt 1 then begin
           a_message = 'Please specify a value greater than or equal to 1.'
         endif else error = 0
       endwhile
       config['AFFINE_A'] = temp
     end


 ;============ MPFIT algorithm ============
   'MPFIT': begin
   ;============ Number of solvers ============
       error = 1
       default_val = '100'
       default = '(Default: '+default_val+')'
       if edit then begin
         if strupcase(config_edit.METHOD) eq 'MPFIT' then begin
           default_val = strtrim(string(config_edit.NSOLVERS, f='(I0)'), 2)
           default = '(Previous configuration: '+default_val+')'
         endif
       endif
       nsolvers_message = 'Please specify the number of times to solve for the best fit SED '+$
                          'using different starting locations in parameters space. '+default
       print_to_width, '======================='
       while error do begin
         print_to_width, ''
         print_to_width, nsolvers_message

         temp = ''
         read,temp, prompt='NSOLVERS: ',form='(A0)'
         if temp eq '' then temp = default_val
         temp = long(double(temp))

         if temp le 0 then begin
           nsolvers_message = 'Please specify a positive value.'
         endif else error = 0
       endwhile
       config['NSOLVERS'] = temp

   ;============ Relative error in sum of squares tolerance ============
       error = 1
       default_val = '1.d-10'
       default = '(Default: '+default_val+')'
       if edit then begin
         if strupcase(config_edit.METHOD) eq 'MPFIT' then begin
           default_val = strtrim(string(config_edit.FTOL, f='(E0)'), 2)
           default = '(Previous configuration: '+default_val+')'
         endif
       endif
       ftol_message = 'Please specify the relative error desired in the sum of squares. '+$
                      'Termination of the MPFIT algorithm occurs when both the actual and '+$
                      'predicted relative reductions in the sum of squares are at most FTOL. '+default
       print_to_width, '======================='
       while error do begin
         print_to_width, ''
         print_to_width, ftol_message

         temp = ''
         read,temp, prompt='FTOL: ',form='(A0)'
         if temp eq '' then temp = default_val
         temp = double(temp)

         if temp le 0 then begin
           ftol_message = 'The relative error desired in the sum of squares must be a positive value.'
         endif else error = 0
       endwhile
       config['FTOL'] = temp

   ;============ Orthogonality tolerance ============
       error = 1
       default_val = '1.d-10'
       default = '(Default: '+default_val+')'
       if edit then begin
         if strupcase(config_edit.METHOD) eq 'MPFIT' then begin
           default_val = strtrim(string(config_edit.GTOL, f='(E0)'), 2)
           default = '(Previous configuration: '+default_val+')'
         endif
       endif
       gtol_message = 'Please specify the orthogonality desired between the function vector and the '+$
                      'columns of the Jacobian matrix. Termination of the MPFIT algorithm '+$
                      'occurs when the cosine of the angle between function vector and any '+$
                      'column of the Jacobian matrix is at most GTOL in absolute value. '+default
       print_to_width, '======================='
       while error do begin
         print_to_width, ''
         print_to_width, gtol_message

         temp = ''
         read,temp, prompt='GTOL: ',form='(A0)'
         if temp eq '' then temp = default_val
         temp = double(temp)

         if temp le 0 then begin
           gtol_message = 'The orthogonality tolerance must be a positive value.'
         endif else error = 0
       endwhile
       config['GTOL'] = temp

   ;============ Relative error in solution tolerance ============
       error = 1
       default_val = '1.d-10'
       default = '(Default: '+default_val+')'
       if edit then begin
         if strupcase(config_edit.METHOD) eq 'MPFIT' then begin
           default_val = strtrim(string(config_edit.XTOL, f='(E0)'), 2)
           default = '(Previous configuration: '+default_val+')'
         endif
       endif
       xtol_message = 'Please specify the relative error desired in the approximate solution. '+$
                      'Termination of the MPFIT algorithm occurs when the relative error between two '+$
                      'consecutive iterates is at most XTOL. '+default
       print_to_width, '======================='
       while error do begin
         print_to_width, ''
         print_to_width, xtol_message

         temp = ''
         read,temp, prompt='XTOL: ',form='(A0)'
         if temp eq '' then temp = default_val
         temp = double(temp)

         if temp le 0 then begin
           xtol_message = 'The relative error desired in the approximate solution must be a positive value.'
         endif else error = 0
       endwhile
       config['XTOL'] = temp

   ;============ Max iterations ============
       error = 1
       default_val = '200'
       default = '(Default: '+default_val+')'
       if edit then begin
         if strupcase(config_edit.METHOD) eq 'MPFIT' then begin
           default_val = strtrim(string(config_edit.MAXITER, f='(I0)'), 2)
           default = '(Previous configuration: '+default_val+')'
         endif
       endif
       maxiter_message = 'Please specify the maximum number of MPFIT iterations to perform. '+default
       print_to_width, '======================='
       while error do begin
         print_to_width, ''
         print_to_width, maxiter_message

         temp = ''
         read,temp, prompt='MAXITER: ',form='(A0)'
         if temp eq '' then temp = default_val
         temp = long(double(temp))

         if temp le 0 then begin
           maxiter_message = 'Please specify a positive value.'
         endif else error = 0
       endwhile
       config['MAXITER'] = temp
     end
 endcase



;======================================    POST-PROCESSING   =========================================================
; Post-processing options
 ;============ Keep Intermediate outputs ============
 error = 1
 default_val = 'NO'
 default = '(Default: '+strlowcase(default_val)+')'
 if edit then begin
   if config_edit.KEEP_INTERMEDIATE_OUTPUT eq 1 then default_val = 'YES' else default_val = 'NO'
   default = '(Previous configuration: '+strlowcase(default_val)+')'
 endif
 keep_message = 'Would you like to keep the intermediate .sav files produced by '+$
                'the fitting algorithm? Please indicate yes or no. '+$
                default
 print_to_width, '======================='
 while error do begin
   print_to_width, ''
   print_to_width, keep_message

   temp = ''
   read,temp, prompt='KEEP_INTERMEDIATE_OUTPUT: ',form='(A0)'
   if temp eq '' then temp = default_val
   temp = strtrim(temp, 2)

   if ~(strupcase(temp) eq 'YES' or strupcase(temp) eq 'NO') then begin
     keep_message = 'Please type yes or no.'
   endif else error = 0
 endwhile
 if strupcase(temp) eq 'YES' then temp = 1 else temp = 0
 config['KEEP_INTERMEDIATE_OUTPUT'] = temp


 ;================    MCMC Post-processing    ===================================
 if (strupcase(config['METHOD']) eq 'MCMC-ADAPTIVE' or strupcase(config['METHOD']) eq 'MCMC-AFFINE') then begin
   ;============ Burn-in iterations ============
   error = 1
   default_val = '0'
   default = '(Default: '+default_val+')'
   if edit then begin
     if strupcase(config_edit.METHOD) eq 'MCMC-ADAPTIVE' or $
        strupcase(config_edit.METHOD) eq 'MCMC-AFFINE' then begin
       default_val = strtrim(string(config_edit.BURN_IN, f='(I0)'), 2)
       default = '(Previous configuration: '+default_val+')'
     endif
   endif
   if strupcase(config['METHOD']) eq 'MCMC-ADAPTIVE' then $
     adaptive_note = 'NOTE: We highly recommend specifying a value rather than using the automatic '+$
                     'calculation as chains can vary widely in the number of autocorrelation times '+$
                     'needed for burn-in when using the MCMC-ADAPTIVE method. ' else adaptive_note = ''

   burnin_message = 'Please specify the number of initial MCMC trials to '+$
                    'truncate as the burn-in phase. If set to 0, then the number '+$
                    'will be chosen automatically from the autocorrelation time. '+$
                    adaptive_note+default
   print_to_width, '======================='
   while error do begin
     print_to_width, ''
     print_to_width, burnin_message

     temp = ''
     read,temp, prompt='BURN_IN: ',form='(A0)'
     if temp eq '' then temp = default_val
     temp = long(double(temp))

     if temp lt 0 then begin
       burnin_message = 'Please specify a non-negative value.'
     endif else error = 0
   endwhile
   config['BURN_IN'] = temp

   ;============ Thinning of chain ============
   error = 1
   default_val = '0'
   default = '(Default: '+default_val+')'
   if edit then begin
     if strupcase(config_edit.METHOD) eq 'MCMC-ADAPTIVE' or $
        strupcase(config_edit.METHOD) eq 'MCMC-AFFINE' then begin
       default_val = strtrim(string(config_edit.THIN_FACTOR, f='(I0)'), 2)
       default = '(Previous configuration: '+default_val+')'
     endif
   endif
   if strupcase(config['METHOD']) eq 'MCMC-ADAPTIVE' then $
     adaptive_note = 'NOTE: We recommend specifying a THIN_FACTOR of 1 as the final distribution will be '+$
                     'from a single chain and each element in the chain will be minimally correlated, '+$
                     'when using the MCMC-ADAPTIVE method. ' else adaptive_note = ''

   thin_message = 'Please specify the factor to thin the MCMC chain after removing '+$
                  'the burn-in trials. For example, a factor of 3 will include only '+$
                  'every 3rd iteration in the final chain. A value of 1 means no thinning. '+$
                  'If set to 0, then the number will chosen automatically from the '+$
                  'autocorrelation time. '+adaptive_note+default
   print_to_width, '======================='
   while error do begin
     print_to_width, ''
     print_to_width, thin_message

     temp = ''
     read,temp, prompt='THIN_FACTOR: ',form='(A0)'
     if temp eq '' then temp = default_val
     temp = long(double(temp))

     if temp lt 0 then begin
       thin_message = 'Please specify a non-negative value.'
     endif else error = 0
   endwhile
   config['THIN_FACTOR'] = temp

   ;============ Number of final iterations ============
   error = 1
   default_val = '1000'
   default = '(Default: '+default_val+')'
   if edit then begin
     if strupcase(config_edit.METHOD) eq 'MCMC-ADAPTIVE' or $
        strupcase(config_edit.METHOD) eq 'MCMC-AFFINE' then begin
       default_val = strtrim(string(config_edit.FINAL_CHAIN_LENGTH, f='(I0)'), 2)
       default = '(Previous configuration: '+default_val+')'
     endif
   endif
   final_message = 'Please specify the number of MCMC trials to include for the final '+$
                   'distributions as taken from the truncated and thinned chain. '+default
   print_to_width, '======================='
   while error do begin
     print_to_width, ''
     print_to_width, final_message

     temp = ''
     read,temp, prompt='FINAL_CHAIN_LENGTH: ',form='(A0)'
     if temp eq '' then temp = default_val
     temp = long(double(temp))

     if temp le 0 then begin
       final_message = 'Please specify a positive value.'
     endif else error = 0
   endwhile
   config['FINAL_CHAIN_LENGTH'] = temp

   ;============ Fraction to generate high res models ============
   error = 1
   default_val = '0'
   default = '(Default: '+default_val+')'
   if edit then begin
     if strupcase(config_edit.METHOD) eq 'MCMC-ADAPTIVE' or $
        strupcase(config_edit.METHOD) eq 'MCMC-AFFINE' then begin
       default_val = strtrim(string(config_edit.HIGH_RES_MODEL_FRACTION, f='(I0)'), 2)
       default = '(Previous configuration: '+default_val+')'
     endif
   endif
   hires_message = 'Please specify the fraction of trials from the final chain, sorted '+$
                   'by quality of fit, from which to generate high resolution models. '+$
                   'If set to 0, then only the best fit high resolution model will be generated. '+$
                   'NOTE: Including more than the best fit high resolution model can drastically '+$
                   'increase file size. Be careful when increasing this value above 0. '+default
   print_to_width, '======================='
   while error do begin
     print_to_width, ''
     print_to_width, hires_message

     temp = ''
     read,temp, prompt='HIGH_RES_MODEL_FRACTION: ',form='(A0)'
     if temp eq '' then temp = default_val
     temp = double(temp)

     if temp lt 0 or temp gt 1 then begin
       hires_message = 'Please specify a value between 0 and 1.'
     endif else error = 0
   endwhile
   config['HIGH_RES_MODEL_FRACTION'] = temp

   ;============ Stranded walker acceptance fraction deviation ============
   if strupcase(config['METHOD']) eq 'MCMC-AFFINE' then begin
     error = 1
     default_val = '2'
     default = '(Default: '+default_val+')'
     if edit then begin
       if strupcase(config_edit.METHOD) eq 'MCMC-AFFINE' then begin
         default_val = strtrim(string(config_edit.AFFINE_STRANDED_DEVIATION, f='(I0)'), 2)
         default = '(Previous configuration: '+default_val+')'
       endif
     endif
     stranded_message = 'Please specify the number of standard deviations a walker must be below the median '+$
                        'acceptance fraction of the ensemble to be considered a stranded walker. '+default
     print_to_width, '======================='
     while error do begin
       print_to_width, ''
       print_to_width, stranded_message

       temp = ''
       read, temp, prompt='AFFINE_STRANDED_DEVIATION: ', form='(A0)'
       if temp eq '' then temp = default_val
       temp = double(temp)

       if temp le 0 then begin
         stranded_message = 'Please specify a positive value.'
       endif else error = 0
     endwhile
     config['AFFINE_STRANDED_DEVIATION'] = temp
   endif
 endif


 return, config.tostruct(/rec)


end
