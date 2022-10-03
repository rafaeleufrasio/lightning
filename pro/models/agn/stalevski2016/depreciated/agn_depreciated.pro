;+
; This file contains some deprecated functionality of the AGN model,
; which is no longer used in the course of Lightning fitting.
;
; Functions in this file are subject to deletion at any time.
;-

function verify_SKIRTOR_parameters, tau97=tau97, p_cl=p_cl, q_cl=q_cl, delta_tor=delta_tor, R_tor=R_tor, i_AGN=i_AGN
    ; Check that the parameters provided are on the grid for the SKIRTOR models.
    ; The grids are hard coded here.
    ; This function does not check that the parameter vectors are non-empty or if they're the same length.
    ; HISTORY:
    ; 2021-02-24: Created (Erik B. Monson)

    compile_opt IDL2

    ; Define the grids for the models.
    tau_97_grid = [3, 5, 7, 9, 11]
    p_grid = [0.0, 0.5, 1.0, 1.5]
    q_grid = [0.0, 0.5, 1.0, 1.5]
    delta_grid = [10, 20, 30, 40, 50, 60, 70, 80]
    R_grid = [10, 20, 30]
    Mcl_grid = [0.97]
    i_grid = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]

    N_models = n_elements(tau97)

    all_good = 1

    for j=0, N_models - 1 do begin

        t_idx = where(tau_97_grid eq tau97[j])
        p_idx = where(p_grid eq p_cl[j])
        q_idx = where(q_grid eq q_cl[j])
        delta_idx = where(delta_grid eq delta_tor[j])
        R_idx = where(R_grid eq R_tor[j])
        i_idx = where(i_grid eq i_AGN[j])

        ; If parameters are not on grid, fail loudly.
        if ((t_idx lt 0) OR (p_idx lt 0) OR (q_idx lt 0) OR $
            (delta_idx lt 0) OR (R_idx lt 0) OR (i_idx lt 0)) then begin
            all_good = 0
            text = string([tau97[j], p_cl[j], q_cl[j], delta_tor[j], R_tor[j], i_AGN[j]], $
                   format='("Supplied model parameters (",(I5 ,D5.1 ,D5.1 ,I5 ,I5 ,I5 ),") are not on available SKIRTOR grid.")')
            message, text
        endif

    endfor

    return, all_good
end

function parse_SKIRTOR_filename, tau97=tau97, p_cl=p_cl, q_cl=q_cl, delta_tor=delta_tor, R_tor=R_tor, i_AGN=i_AGN,$
         lightning_folder=lightning_folder
    ; Construct paths to SKIRTOR files given model parameters.
    ; HISTORY:
    ; 2021-02-24: Created (Erik B. Monson)

    compile_opt IDL2

    ; For now we'll assume that the relevant input-checking
    ; has already been done in the functions calling this one, so
    ; we won't check that the parameter vectors are non-empty and
    ; the same length, or that the parameters are all on-grid.
    ; I might move that error handling into this
    ; function later on.

    ; Not a fan of making arrays that grow over a loop
    ; but it's hard not to when the entries are strings
    ; of a different length.
    paths = !null

    N_files = n_elements(tau97)

    for j=0, N_files - 1 do begin

        ; Assemble the file path for the model
        ; File name example: t5_p1_q0_oa50_R20_Mcl0.97_i30_sed.idl
        file_path = lightning_folder + 'agn_emission/stalevski2016/SKIRTOR/' + 'i' + strtrim(i_AGN[j],2) + '/'

        ; Format the model parameters
        t_str = strtrim(tau97[j],2)
        ; They made this as annoying as possible
        ; and I'm probably making it worse.
        p_fmt = !null
        q_fmt = !null
        if ((p_cl[j] eq 0.0) || (p_cl[j] eq 1.0)) then p_fmt = '(I)' else p_fmt='(D3.1)'
        if ((q_cl[j] eq 0.0) || (q_cl[j] eq 1.0)) then q_fmt = '(I)' else q_fmt='(D3.1)'
        p_str = strtrim(string(p_cl[j], format=p_fmt), 2)
        q_str = strtrim(string(q_cl[j], format=q_fmt), 2)
        oa_str = strtrim(delta_tor[j], 2)
        R_str = strtrim(R_tor[j], 2)
        Mcl_str = '0.97' ; Hard coded for now; the models are only gridded with Mcl = 0.97
        i_str = strtrim(i_AGN[j], 2)

        file_path = file_path + 't' + t_str + '_p' + p_str + '_q' + q_str + '_oa' + oa_str + $
                    '_R' + R_str + '_Mcl' + Mcl_str + '_i' + i_str + '_sed.idl'

        paths = [paths, file_path]
    endfor

    return, paths
end

function SKIRTOR_template, tau97=tau97, p_cl=p_cl, q_cl=q_cl, delta_tor=delta_tor, R_tor=R_tor, i_AGN=i_AGN,$
         filter_labels=filter_labels, z_shift=z_shift, lightning_folder=lightning_folder
    ; Load SKIRTOR model for given params, convert to Lnu, and compute the SED in the given filters.
    ; For the DL07 dust models we load the entire modelset at the start of fitting, but that's not feasible
    ; for the SKIRTOR models (there are 19200 of them and they take up about 150 MB even in binary form)
    ; So, we load them one at a time. Unless that turns out to be a speed bottleneck in testing.
    ; HISTORY:
    ; 2021-02-17: Created (Erik B. Monson)

    compile_opt IDL2

    ; If the parameters aren't defined we'll default to some values for consistency with existing Lightning code
    if (n_elements(tau97) eq 0) then tau97 = 3 else tau97 = (tau97)[0]
    if (n_elements(p_cl)  eq 0) then p_cl = 0.0 else p_cl = (p_cl)[0]
    if (n_elements(q_cl)  eq 0) then q_cl = 0.0 else q_cl = (q_cl)[0]
    if (n_elements(delta_tor) eq 0) then delta_tor = 10 else delta_tor = (delta_tor)[0]
    if (n_elements(R_tor)  eq 0) then R_tor = 10 else R_tor = (R_tor)[0]
    if (n_elements(i_AGN) eq 0) then i_AGN = 0 else i_AGN = (i_AGN)[0]

    allgood = verify_SKIRTOR_parameters(tau97=tau97,  p_cl=p_cl, q_cl=q_cl,$
              delta_tor=delta_tor, R_tor=R_tor, i_AGN=i_AGN)

    if (n_elements(filter_labels) eq 0) then $
      filter_labels = ['GALEX_FUV', 'UVOT_W2', 'UVOT_M2', 'GALEX_NUV', 'UVOT_W1', 'SDSS_u', 'SDSS_g', 'SDSS_r',$
                       'SDSS_i', 'SDSS_z', '2MASS_J', '2MASS_H', '2MASS_Ks', 'W1', 'SPITZER_I1', 'SPITZER_I2', 'W2',$
                       'W3', 'W4', 'HAWC+A', 'HAWC+B', 'HAWC+C', 'HAWC+D', 'HAWC+E']
    if (n_elements(z_shift) eq 0) then z_shift = 0.0
    if (n_elements(lightning_folder) eq 0) then lightning_folder = '~/Lightning/'
    if (not lightning_folder.EndsWith('/')) then lightning_folder = lightning_folder + '/'

    ; Assemble the file path for the model
    file_path = parse_SKIRTOR_filename(tau97=tau97, p_cl=p_cl, q_cl=q_cl, delta_tor=delta_tor, R_tor=R_tor,$
                                        i_AGN=i_AGN, lightning_folder=lightning_folder)

    ; Load the file -- this defines the structure `SKIRTOR_model` saved in the file
    restore, file_path

    ; Observed-frame frequency and wavelength
    nu_obs = SKIRTOR_model.nu_rest / (1.d0 + z_shift)
    wave_obs = SKIRTOR_model.wave_rest * (1.d0 + z_shift)

    ; Compute "observed-frame Lnu"
    ; Lnu_obs = Lnu_rest * (1+z_shift)

    Lnu_obs_total = SKIRTOR_model.lnu_total * (1 + z_shift)
    Lnu_obs_disk_direct = SKIRTOR_model.lnu_disk_direct * (1 + z_shift)
    Lnu_obs_disk_scattered = SKIRTOR_model.lnu_disk_scattered * (1 + z_shift)
    Lnu_obs_dust_total = SKIRTOR_model.lnu_dust_total * (1 + z_shift)
    Lnu_obs_dust_scattered = SKIRTOR_model.lnu_dust_scattered * (1 + z_shift)
    Lnu_obs_transparent = SKIRTOR_model.lnu_transparent * (1 + z_shift)

    ; Observe the total spectrum
    ; Note that GET_FILTERS wants Lnu to have dimensions [n_models, n_wave]
    ; so we need to pass it the transpose to get a [1, n_wave] shaped array.
    get_filters, filter_labels, nu_obs, Filters, Lnu=transpose(Lnu_obs_total),$
                 filters_dir=lightning_folder+'Filters/', mean_Lnu=mean_Lnu,$
                 mean_wave=mean_wave

    ; This structure probably holds more information than necessary
    template = {MODEL_NAME: SKIRTOR_model.model_name,$
                TAU_97: tau97,$
                P: p_cl,$
                Q: q_cl,$
                OA: delta_tor,$
                R: R_tor,$
                INCLINATION: i_AGN,$
                Z_SHIFT: z_shift,$
                LNU_TOTAL: Lnu_obs_total,$ ; The Lnu in this output structure will be the
                LNU_DISK_DIRECT: Lnu_obs_disk_direct,$ ; observed frame Lnu. I'm dropping the
                LNU_DISK_SCATTER: Lnu_obs_disk_scattered,$ ; _obs tag to avoid confusion with the
                LNU_DUST_TOTAL: Lnu_obs_dust_total,$ ; 'observed' SED in the provided filters.
                LNU_DUST_SCATTER: Lnu_obs_dust_scattered,$
                LNU_TRANSPARENT: Lnu_obs_transparent,$
                WAVE_REST: SKIRTOR_model.wave_rest,$ ; rest-frame wavelength in microns
                WAVE_OBS: wave_obs,$ ; observed-frame wavelength in microns
                NU_REST: SKIRTOR_model.nu_rest,$ ; rest-frame frequency in Hz
                NU_OBS: nu_obs,$ ; observed-frame frequency in Hz
                LNU_TOTAL_OBS: reform(mean_Lnu),$ ; Lnu observed in filters
                WAVE_FILTERS: reform(mean_wave),$ ; Mean wavelength of filters
                FILTER_LABELS: filter_labels,$ ; Names of filters
                FILTERS: Filters} ; Structure containing filter info

    return, template
end

function SKIRTOR_sed_vector, tau97=tau97, p_cl=p_cl, q_cl=q_cl, delta_tor=delta_tor, R_tor=R_tor, i_AGN=i_AGN,$
         filter_labels=filter_labels, z_shift=z_shift, lightning_folder=lightning_folder, mean_wave=mean_wave
    ; Load SKIRTOR models for given params and compute the SEDs in the given filters.
    ; Vectorized version: this function expects N values of each parameter and returns
    ; N-many SEDs.
    ; Following the example of the DL07 implementation, this returns /just/ the SEDs, with no
    ; wavelength information or anything.
    ; HISTORY:
    ; 2021-02-24: Created (Erik B. Monson)

    compile_opt IDL2

    ; If the parameters aren't defined we'll default to some values for consistency with existing Lightning code

    n_tau = n_elements(tau97)
    n_p = n_elements(p_cl)
    n_q = n_elements(q_cl)
    n_delta = n_elements(delta_tor)
    n_R = n_elements(R_tor)
    n_i = n_elements(i_AGN)

    if (n_tau eq 0) then tau97 = 3
    if (n_p  eq 0) then p_cl = 0.0
    if (n_q  eq 0) then q_cl = 0.0
    if (n_delta eq 0) then delta_tor = 10
    if (n_R  eq 0) then R_tor = 10
    if (n_i eq 0) then i_AGN = 0

    if ((n_tau ne n_p) OR (n_tau ne n_q) OR (n_tau ne n_delta) OR (n_tau ne n_R) OR (n_tau ne n_i) OR $
    (n_p ne n_q) OR (n_p ne n_delta) OR (n_p ne n_R) OR (n_p ne n_i) OR $
    (n_q ne n_delta) OR (n_q ne n_R) OR (n_q ne n_i) OR (n_delta ne n_q) OR $
    (n_delta ne n_R) OR (n_delta ne n_i) OR (n_R ne n_i)) then message, 'Input parameter vectors must have same length.'

    N_models = n_tau

    ; Check if the model parameters are on-grid
    allgood = verify_SKIRTOR_parameters(tau97=tau97, p_cl=p_cl, q_cl=q_cl,$
              delta_tor=delta_tor, R_tor=R_tor, i_AGN=i_AGN)

    if (n_elements(filter_labels) eq 0) then $
      filter_labels = ['GALEX_FUV', 'UVOT_W2', 'UVOT_M2', 'GALEX_NUV', 'UVOT_W1', 'SDSS_u', 'SDSS_g', 'SDSS_r',$
                       'SDSS_i', 'SDSS_z', '2MASS_J', '2MASS_H', '2MASS_Ks', 'W1', 'SPITZER_I1', 'SPITZER_I2', 'W2',$
                       'W3', 'W4', 'HAWC+A', 'HAWC+B', 'HAWC+C', 'HAWC+D', 'HAWC+E']
    N_filters = n_elements(filter_labels)
    if (n_elements(z_shift) eq 0) then z_shift = 0.0
    if (n_elements(lightning_folder) eq 0) then lightning_folder = '~/Lightning/'
    if (not lightning_folder.EndsWith('/')) then lightning_folder = lightning_folder + '/'

    file_paths = parse_SKIRTOR_filename(tau97=tau97, p_cl=p_cl, q_cl=q_cl, delta_tor=delta_tor, R_tor=R_tor,$
                                        i_AGN=i_AGN, lightning_folder=lightning_folder)

    seds = dblarr(N_filters, N_models)

    for j=0, N_models - 1 do begin

        restore, file_paths[j]

        nu_obs = SKIRTOR_model.nu_rest / (1.d0 + z_shift)
        Lnu_obs_total = SKIRTOR_model.lnu_total * (1 + z_shift)

        ; Observe the total spectrum
        ; Note that GET_FILTERS wants Lnu to have dimensions [n_models, n_wave]
        ; so we need to pass it the transpose to get a [1, n_wave] shaped array.
        get_filters, filter_labels, nu_obs, Filters, Lnu=transpose(Lnu_obs_total),$
                     filters_dir=lightning_folder+'Filters/', mean_Lnu=mean_Lnu,$
                     mean_wave=mean_wave

        seds[*,j] = mean_Lnu

    endfor

    return, seds
end


function SKIRTOR_spec_vector, tau97=tau97, p_cl=p_cl, q_cl=q_cl, delta_tor=delta_tor, R_tor=R_tor, i_AGN=i_AGN,$
         filter_labels=filter_labels, z_shift=z_shift, lightning_folder=lightning_folder, nu_obs=nu_obs
    ; Load SKIRTOR models for given params and return the total AGN emission as Lnu_obs.
    ; Vectorized version: this function expects N values of each parameter and returns
    ; N-many SEDs.
    ; Following the example of the DL07 implementation, this returns /just/ the SEDs, with no
    ; wavelength information or anything.
    ; HISTORY:
    ; 2021-02-24: Created (Erik B. Monson)

    compile_opt IDL2

    ; If the parameters aren't defined we'll default to some values for consistency with existing Lightning code
    n_tau = n_elements(tau97)
    n_p = n_elements(p_cl)
    n_q = n_elements(q_cl)
    n_delta = n_elements(delta_tor)
    n_R = n_elements(R_tor)
    n_i = n_elements(i_AGN)

    if (n_tau eq 0) then tau97 = 3
    if (n_p  eq 0) then p_cl = 0.0
    if (n_q  eq 0) then q_cl = 0.0
    if (n_delta eq 0) then delta_tor = 10
    if (n_R  eq 0) then R_tor = 10
    if (n_i eq 0) then i_AGN = 0

    if ((n_tau ne n_p) OR (n_tau ne n_q) OR (n_tau ne n_delta) OR (n_tau ne n_R) OR (n_tau ne n_i) OR $
    (n_p ne n_q) OR (n_p ne n_delta) OR (n_p ne n_R) OR (n_p ne n_i) OR $
    (n_q ne n_delta) OR (n_q ne n_R) OR (n_q ne n_i) OR (n_delta ne n_q) OR $
    (n_delta ne n_R) OR (n_delta ne n_i) OR (n_R ne n_i)) then message, 'Input parameter vectors must have same length.'

    N_models = n_tau

    ; Check if the model parameters are on-grid
    allgood = verify_SKIRTOR_parameters(tau97=tau97, p_cl=p_cl, q_cl=q_cl,$
              delta_tor=delta_tor, R_tor=R_tor, i_AGN=i_AGN)

    if (n_elements(z_shift) eq 0) then z_shift = 0.0
    if (n_elements(lightning_folder) eq 0) then lightning_folder = '~/Lightning/'
    if (not lightning_folder.EndsWith('/')) then lightning_folder = lightning_folder + '/'

    file_paths = parse_SKIRTOR_filename(tau97=tau97, p_cl=p_cl, q_cl=q_cl, delta_tor=delta_tor, R_tor=R_tor,$
                                        i_AGN=i_AGN, lightning_folder=lightning_folder)

    ; Need a way to know how many wavelength elements are in the
    ; files; or maybe not, since there's alway 132 for every one of the
    ; SKIRTOR files.
    spec = dblarr(132, N_models)

    for j=0, N_models - 1 do begin

        restore, file_paths[j]

        nu_obs = SKIRTOR_model.nu_rest / (1.d0 + z_shift)
        Lnu_obs_total = SKIRTOR_model.lnu_total * (1.d0 + z_shift)

        spec[*,j] = Lnu_obs_total

    endfor

    return, spec
end
