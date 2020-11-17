; Example IDL code for fitting the SED of J123548.94+621144.7
; Need the current directory to be the lightning directory
; To run, type the following into the IDL command line: @ ./Example/lightning_example.pro

; Creates system variable for physical constants
 @ ./startup_lightning.pro

; Compile Lightning functions and programs
.r ./dust_emission/DL07/dl07_spec_vector.pro
.r ./lightning_MCMC.pro
.r ./MCMC_savefits_Tuffs.pro
.r ./MCMC_savefits_Calz.pro

; Read in the data to be fit
; Luminosity data are in L_nu (L_sun Hz^-1) and have been recalibrated as stated in section 2.1 
; of Doore et al. (2020, in prep).
data=mrdfits('./Example/J123548.94+621144.7.fits',1)

; Define output folder and lightning folder
outfolder='./Example/'
lightning_folder='./'

; Read in precomputed file of L_star_abs for the energy balance assumption
L_star_abs_table=mrdfits('./L_star_abs_model_table/L_star_abs_model_table_z_'+string(data.redshift,f='(F4.2)')+'.fits',1)


; Limit parameters to set ranges or values for Tuffs Model
rant=dblarr(2,10)
rant[*,0] = [ 0.d,8.0d]        ;tau_b_f
rant[*,1] = [ 0.d,1.0d]        ;r_disk
rant[*,2] = [0.d,0.61d]        ;F
rant[*,3] = [ 0.d,1.0d]        ;cosi
rant[*,4] = [ 0.d,0.0d]        ;r_bulge
rant[*,5] = [ 2.d0, 2.d0]      ;alpha
rant[*,6] = [ 0.7d0,25]        ;Umin
rant[*,7] = [ 3.d5,  3.d5]     ;Umax
rant[*,8] = [ 0.d0,  1.d0]     ;gamma_dust
rant[*,9]= [ 0.0047,0.0458d0]  ;q_PAH

; Randomly generate starting point for each parameter and SFH coefficient
parameter_start_t=randomu(seed,10)*rebin(reform(rant[1,*]-rant[0,*]),10)+rebin(reform(rant[0,*]),10)
coeff_start_t=randomu(seed,5)*10

; Limit parameters to set ranges or values for Calzetti Model
ranc=dblarr(2,10)
ranc[*,0] = [ 0.d,3.0d]       ;tau_v_diff
ranc[*,1] = [ 0.d,0.0d]       ;
ranc[*,2] = [ 0.d,0.0d]       ;
ranc[*,3] = [ 0.d,0.0d]       ;
ranc[*,4] = [ 0.d,0.0d]       ;
ranc[*,5] = [ 2.d0, 2.d0]     ;alpha
ranc[*,6] = [ 0.7d0,25]       ;Umin
ranc[*,7] = [ 3.d5,  3.d5]    ;Umax
ranc[*,8] = [ 0.d0,  1.d0]    ;gamma_dust
ranc[*,9]= [ 0.0047,0.0458d0] ;q_PAH

; Randomly generate starting point for each parameter and SFH coefficient
parameter_start_c=randomu(seed,10)*rebin(reform(ranc[1,*]-ranc[0,*]),10)+rebin(reform(ranc[0,*]),10)
coeff_start_c=randomu(seed,5)*10


; Run lightning for the Tuffs model with image-based prior on inclination and the Calzetti model for comparison
; Each run is saved separately in the output folder as a .idl file
lightning_MCMC_vector,data.Lnu_obs,data.Lnu_unc,0.d0,0.d0,data.galaxy_id+'_tuffs',strtrim(data.filter_labels,2),$
         z_shift=data.redshift,Ntrials=2e5,lightning_folder=lightning_folder,$
         /par5_constant,outfolder=outfolder,/adaptive,/dust_emission,$
         /Tuffs,L_star_abs_table=L_star_abs_table,/L_star_abs_model_table,parameter_start=parameter_start_t,coeff_start=coeff_start_t,$
         /par6_constant,/par8_constant,prior_dist=data,/use_priors

lightning_MCMC_vector,data.Lnu_obs,data.Lnu_unc,0.d0,0.d0,data.galaxy_id+'_calz',strtrim(data.filter_labels,2),$
         z_shift=data.redshift,lightning_folder=lightning_folder,/par2_constant,/par3_constant,outfolder=outfolder,$
         /adaptive,/dust_emission,/calzetti,/par6_constant,/par8_constant,parameter_start=parameter_start_c,coeff_start=coeff_start_c
         
         
; Save the lightning output into a more convenient fits file for analysis
MCMC_savefits_tuffs,data,strtrim(data.filter_labels,2),data.lnu_obs,data.lnu_unc,outfolder,outfolder,5000,5,$
         file_name='MCMC_lightning_tuffs',/prior,LTIR_table_loc='./L_star_abs_model_table/',lightning_folder=lightning_folder

MCMC_savefits_calz,data,strtrim(data.filter_labels,2),data.lnu_obs,data.lnu_unc,outfolder,outfolder,5000,5,$
         file_name='MCMC_lightning_calz',lightning_folder=lightning_folder,/calzetti

