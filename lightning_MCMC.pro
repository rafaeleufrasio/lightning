function percentile,array,percentile
;routine sorts an input array and returns an interpolated required percentile
;  it assumes the lowest value of the array to be the 0-th percentile
;  and the highest to be the 100-th percentile

Nsample = array.length
percentile_grid = dindgen(Nsample)/(Nsample-1)

return,interpol(array[sort(array)],percentile_grid,percentile)

end


function steps_stellar,filter_labels=filter_labels,steps_bounds=steps_bounds,$
                       z_shift=z_shift,dtime_SF=dtime_SF,nolines=nolines,$
                       nonebular=nonebular,Zmet=Zmet,IMF=IMF,Filters=Filters,$
                       lightning_folder=lightning_folder
; this function generates the spectra, SEDs, and stellar parameters, given a set of filters and steps boundaries
;
;Modification History
; Rafael Eufrasio - Corrected (1+z) factor - March 9th 2020
; Keith Doore     - Replaced if statements with n_elements() on keywords to use keyword_set() - April 27th 2020

if (n_elements(filter_labels) eq 0) then $
  filter_labels=['GALEX_FUV','UVOT_W2','UVOT_M2','GALEX_NUV','UVOT_W1','SDSS_u','SDSS_g','SDSS_r','SDSS_i','SDSS_z',$
                 '2MASS_J','2MASS_H','2MASS_Ks','W1','SPITZER_I1','SPITZER_I2','W2']
if (n_elements(steps_bounds) eq 0) then steps_bounds = [0.d0,1.d7,1.d8,1.d9,5.d9,13.6d9]
if (n_elements(z_shift) eq 0) then z_shift = 0.0
if (n_elements(dtime_SF) eq 0) then dtime_SF = 5.d5
if ~keyword_set(nolines) then nolines = 0
if ~keyword_set(nonebular) then nebtag='nebular_' else nebtag=''
if (n_elements(Zmet) eq 0) then Ztag='0.020' else Ztag=string(Zmet,form='(F5.3)')
if (n_elements(IMF) eq 0) then IMFtag='Kroupa01' else IMFtag=IMF
if (n_elements(lightning_folder) eq 0) then Lightning='~/Lightning/' else Lightning=lightning_folder
;lightning_folder is the folder that contains both Filters/ and models/ subfolders

Nsteps = N_elements(steps_bounds) - 1

bursts_folder = Lightning+'models/04-single_burst/'
filename = bursts_folder+IMFtag+'/'+IMFtag+'_Z'+Ztag+'_'+nebtag+'spec.idl'

sObj = OBJ_NEW('IDL_Savefile',filename)
  sObj->RESTORE,'time'
  sObj->RESTORE,'wave'
  sObj->RESTORE,'nu'
  sObj->RESTORE,'Lnu'
  sObj->RESTORE,'Lbol'
  sObj->RESTORE,'Nlyc'
  sObj->RESTORE,'Mstars'
  Ntime = time.length
  Nwave_burst = wave.length

    sObj->RESTORE,'L_lines'
    sObj->RESTORE,'wlines'
    Nwlines=wlines.length
if nolines eq 0 then begin
  vd = 5d6 ; 50km/s = 5,000,000 cm/s
  z = vd/!cv.clight
  for tt=0,(Ntime-1) do begin
    Lnu_lines = wave*0
    for ll=0,(Nwlines-1) do begin
      x2 = ((wave/(wlines[ll])[0] - 1.d0)/z)^2
      Lnu_line = exp(-x2)
      Lnu_line /= abs(trap_int(nu,Lnu_line))
      Lnu_lines += (L_lines[tt,ll])[0]*Lnu_line
    endfor
    Lnu[tt,*] += Lnu_lines
  endfor
endif

age_burst   = temporary(time)
nu_rest     = temporary(nu)            ; restframe frequency
nu_obs      = nu_rest/(1.d0+z_shift)   ; observed frequency
wave_rest   = temporary(wave)          ; restframe wavelength
wave_obs    = wave_rest*(1.d0+z_shift) ; observed wavelength
Lnu_burst_rest = temporary(Lnu)                  ; restframe Lnu spectrum
Lnu_burst      = Lnu_burst_rest * (1.d0+z_shift) ; observed Lnu spectrum, where Lnu = 4*!dpi*DL^2*Fnu
;Lbol_burst = Lbol/(1.d0+z_shift)  & Lbol   = !null
;Q0_burst   = Nlyc/(1.d0+z_shift)  & Nlyc   = !null
;Lnu_burst   = Lnu                 & Lnu    = !null
Lbol_burst  = Lbol                 & Lbol   = !null ; observed = restframe bolometric luminosity 
Q0_burst    = Nlyc                 & Nlyc   = !null ; restframe ionizing flux rate
;Q0_burst    = (1.d0+z_shift)*Nlyc & Nlyc   = !null ; observed ionizing flux rate
Mstar_burst = Mstars               & Mstars = !null
GET_FILTERS,filter_labels,nu_obs,Filters,Lnu=[1.d0] # [0*(nu_obs) + 3631],$
            filters_dir=Lightning+'Filters/',mean_wave=mean_wave;,$
            ;/plot_filters;,mean_nu=mean_nu,mean_Lnu=mean_Lnu,mean_Lnu=mean_Lnu,sigma_wave=sigma_wave

wave_filters = reform(mean_wave)

Q0_steps    = dblarr(Nsteps)
Lnu_steps   = dblarr(Nwave_burst,Nsteps)
Lbol_steps  = dblarr(Nsteps)
Mstar_steps = dblarr(Nsteps)

dtime_SF = 5.d5
;dtime_SF = 1.d5
for st=0,(Nsteps-1) do begin
  ti = steps_bounds[st]
  tf = steps_bounds[st+1]
  dt_step = tf - ti
  Ntime_SF = tf/dtime_SF + 1
  time_SF = dtime_SF*dindgen(Ntime_SF) ;time from the onset of SF, progressing linearly
  Q0_steps[st]    = trap_int(time_SF,interpol(Q0_burst,   age_burst,tf-time_SF),xrange=[0,dt_step])
  Lbol_steps[st]  = trap_int(time_SF,interpol(Lbol_burst, age_burst,tf-time_SF),xrange=[0,dt_step])
  Mstar_steps[st] = trap_int(time_SF,interpol(Mstar_burst,age_burst,tf-time_SF),xrange=[0,dt_step])
  for kk=0,(nu_obs.length-1) do begin
    Lnu_steps[kk,st] = trap_int(time_SF,interpol(reform(Lnu_burst[*,kk]),age_burst,tf-time_SF),xrange=[0,dt_step])
  endfor
endfor

GET_FILTERS,filter_labels,nu_obs,Filters,Lnu=transpose(Lnu_steps),$
            filters_dir=Lightning+'Filters/',mean_Lnu=mean_Lnu;,$
            ;/plot_filters;,mean_wave=mean_wave,mean_nu=mean_nu,sigma_wave=sigma_wave

mean_Lnu_steps = transpose(temporary(mean_Lnu))

steps = {filter_labels : filter_labels, $
         wave_filters  : wave_filters,  $
         mean_Lnu      : mean_Lnu_steps,$
         wave_rest     : wave_rest,     $ ; restframe wavelength
         wave_obs      : wave_obs,      $
         Lnu           : Lnu_steps,     $ ; observed Lnu spectrum, where Lnu = 4*!dpi*DL^2*Fnu = (1+z)*Lnu_rest 
         Q0            : Q0_steps,      $
         Mstar         : Mstar_steps,   $
         Lbol          : Lbol_steps,    $
         bounds        : steps_bounds,  $
         Filters       : Filters,       $
         z_shift       : z_shift}

return,steps

end




function matrix_unred_Keith_vector,wave,tauV_DIFF_arr=tauV_DIFF_arr,delta_arr=delta_arr,tauV_BC_arr=tauV_BC_arr,$
                            Calzetti_attenuation=Calzetti_attenuation


;Modification History
; Keith Doore     - Added ability to run with pure Calzetti curve - May 6th 2020


;wave[Nwave] in microns
if (n_elements(tauV_DIFF_arr) eq 0) then tauV_DIFF_arr = 0.d0
if (n_elements(delta_arr) eq 0)  then delta_arr  = 0.d0
if (n_elements(tauV_BC_arr) eq 0) then tauV_BC_arr = 0.d0

if min(tauv_diff_arr) gt 0 then message,'Error -- TauV_diff must be negative'

Nwave = wave.length
n1 = tauV_DIFF_arr.length
n2 = delta_arr.length
n3 = tauV_BC_arr.length

if n1 ne n2 or n1 ne n3 or n2 ne n3 then message,'Input parameters must have same length for vectorization'
nn=n1

Dlam = dblarr(Nwave,nn)
  FWHM_bump = 0.0350d0 ; 350 Angs
  lambda0_bump = 0.2175d0   ; 2175 Angs
  ;Ebump = 0.85 - 1.9*delta ;
  ;Dlam = Ebump / ( ((wave^2 - lambda0_bump^2)/(wave*FWHM_bump))^2 + 1.d0 )
  Dlam = (0.85d0-1.9d0*rebin(reform(delta_arr,1,nn),nwave,nn)) / $
         (rebin((((wave^2 - lambda0_bump^2)/(wave*FWHM_bump))^2 + 1.d0),nwave,nn))

;Calzetti extinction curve
klam = wave*0.d0 
  RV = 4.05d0 ;Calzetti et al. 2000
  ;AV = 1.0
  ;EBV = AV / RV
  ;
  ; Calzettic curve extrapolated from 0.01 to 5.0 microns
  ; originally derived from 0.0912 to 2.2 microns
  w1 = where((wave GE 0.6300) AND (wave LE 2.2000),nw1,/null)
  w2 = where((wave GE 0.0912) AND (wave LT 0.6300),nw2,/null)
  w3 = where((wave LT 0.0912) AND (wave GT 2.2000),nw3,/null)
  x  = 1.d0/wave                      ;Wavelength in inverse microns
  if (nw1 ne 0) then klam[w1] = 2.659d0*(-1.857d0 + 1.040d0*x[w1]) + RV
  if (nw2 ne 0) then klam[w2] = 2.659d0*(poly(x[w2], [-2.156, 1.509d0, -0.198d0, 0.011d0])) + RV
  if (nw3 ne 0) then klam[w3] = 0.

flam1 = dblarr(Nwave,nn)
if ~keyword_set(Calzetti_attenuation) then begin
  flam1 = (rebin(klam,nwave,nn) + Dlam)/4.05 * $
          (rebin(wave,nwave,nn)/0.55)^(rebin(reform(delta_arr,1,nn),nwave,nn))
          ;flam1(lambda,delta) = tau_DIFF(lambda,delta) / tauV_DIFF
endif else begin
  flam1 = rebin(klam,nwave,nn)/4.05
endelse

; birth cloud attenuation
;flam2 = tau_BC(lambda) / tauV_BC
flam2 = 1.d0 / (wave/0.55d0)

tau_lam_diff_4D = rebin(flam1,nwave,nn)*rebin(reform(tauv_diff_arr,1,nn),nwave,nn)

tau_lam_BC_4D = flam2 # tauV_BC_arr

if ~keyword_set(Calzetti_attenuation) then begin
  exp_tau = exp(temporary(tau_lam_DIFF_4D) + temporary(tau_lam_BC_4D))
endif else begin
  exp_tau = exp(temporary(tau_lam_DIFF_4D))
endelse

;exp_tau[k,k1,k2,k3] = exp(tauV_DIFF[k1]*flam1[k,k2] + tauV_BC[k3]*flam2[k])
;Lnu[k]*exp_tau[k,k1,k2,k3]
;Lnu_unred[k,k1,k2,k3] = Lnu[k]*exp(tauV_DIFF[k1]*flam1[k,k2] + tauV_BC[k3]*flam2[k])

return,exp_tau


end



function tuffs_matrix_unred_vector,wave,tauB_f_arr=tauB_f_arr,F_arr=F_arr,cosi_arr=cosi_arr, $
			rold0_arr=rold0_arr,b_to_d_arr=b_to_d_arr,lightning_folder=lightning_folder

;Modification History
; Keith Doore     - modified to have rdisk as intrinsic property and rbulge as B/D - April 13th 2021

;input wavelengths must be in microns
if (n_elements(tauB_f_arr) eq 0) then tauB_f_arr = 1.d0
if (n_elements(rold0_arr) eq 0)  then rold0_arr  = 0.d0
if (n_elements(b_to_d_arr) eq 0) then b_to_d_arr = 0.d0
if (n_elements(F_arr) eq 0) 	 then F_arr      = 0.d0
if (n_elements(cosi_arr) eq 0)   then cosi_arr   = 1.d0
if (n_elements(lightning_folder) eq 0) then lightning_folder='~/Lightning/'


if max(cosi_arr) gt 1 or min(cosi_arr) lt 0 then message,'Error -- cosi must be between 0 and 1'
if max(F_arr) ge 0.61 or min(F_arr) lt 0 then message,'Error -- F must be between 0 and 0.61'
if min(b_to_d_arr) lt 0 then message,'Error -- b_to_d must be greater than or equal to 0'
if min(rold0_arr) lt 0 or max(rold0_arr) gt 1 then message,'Error -- rold0 must be between 0 and 1'
if min(tauB_f_arr) lt 0 then message,'Error -- TauB_f must be greater than or equal to 0'


;coeff.dat file from ftp site (http://cdsarc.u-strasbg.fr/ftp/J/A+A/527/A109/coeff.dat) in Popescu et al (2011)
;converted to convenient 3d matrix where [tau,wavebands,coeff] are the dimensions and values
restore,lightning_folder+'tuffs_coeff.idl'

tau=[0.1d,0.3,0.5,1.0,2.0,4.0,8.0]  ;tau values used in Table 3 in Popescu et al (2011)

;Include waveband at 50000 Angstroms where mlambda is 0 for interpolation smoothness
wavebands=[912.d,1350,1500,1650,2000,2200,2500,2800,3650,4430,5640,8090,12590,22000,50000]   ;wavelengths of the bands used (Angstroms)

;Values from Table E.4 in Popescu et al (2011)
Fcal=0.35d0
flamb_fcal_column=[0.427d,0.484,0.527,0.545,0.628,0.676,0.739,0.794,0.892,0.932,0.962,0.984,0.991,0.994,0.999]

;[value=(1−Fcal*fλ)]->[(1-value)/Fcal=fλ]
flambda=(1-flamb_fcal_column)/Fcal ;fλ values

;Convert b_to_d ratio into rbulge and rdisk
;rdisk+rbulge=rold and rdisk*B/D=rbulge
rdisk0_arr=rold0_arr/(1+b_to_d_arr)
rbulge0_arr=rold0_arr*b_to_d_arr/(1+b_to_d_arr)

;Add a tau at 0 for interpolation smoothness
;Set the diffuse attenuation to 0 for this case
tau=[0.0,tau]
ntau=n_elements(tau)
nwave=n_elements(wave)
ncosi=n_elements(cosi_arr)
ntauB_f=n_elements(tauB_f_arr)
nF=n_elements(f_arr)
nrdisk0=n_elements(rdisk0_arr)
nrbulge0=n_elements(rbulge0_arr)
nwavebands=n_elements(wavebands)

if ncosi ne ntaub_f or ncosi ne nf or ncosi ne nrdisk0 or ncosi ne nrbulge0 or ntauB_f ne nf or ntauB_f ne nrdisk0 or ntauB_f ne nrbulge0 or $
nF ne nrdisk0 or nF ne nrbulge0 or nrdisk0 ne nrbulge0 then message,'Input parameters must have same length for vectorization'
nn=ncosi

mlambda_d=dblarr(ntau,nwavebands,nn)
mlambda_td=dblarr(ntau,nwavebands,nn)
mlambda_b=dblarr(ntau,nwavebands,nn)

;Calculate delta_mlambda for each component of galaxy  
;Using Equation 6 from Tuffs et al (2004)
;First tau index is tau=0 and last wavelength index is wavelength=50000 Angstroms. Set both to 0 for mlambda
One_minus_cosi_arr=rebin(reform((1-cosi_arr),1,1,nn),ntau,nwavebands,nn)
for i=0,(nn-1) do begin
	mlambda_d[1:*,0:-2,i]=(One_minus_cosi_arr[1:*,0:-2,i])^0.d0*adisk[*,*,0]+(One_minus_cosi_arr[1:*,0:-2,i])^1.d0*adisk[*,*,1]+$
                          (One_minus_cosi_arr[1:*,0:-2,i])^2.d0*adisk[*,*,2]+(One_minus_cosi_arr[1:*,0:-2,i])^3.d0*adisk[*,*,3]+$
                          (One_minus_cosi_arr[1:*,0:-2,i])^4.d0*adisk[*,*,4]+(One_minus_cosi_arr[1:*,0:-2,i])^5.d0*adisk[*,*,5]
	mlambda_td[1:*,0:-2,i]=(One_minus_cosi_arr[1:*,0:-2,i])^0.d0*atdisk[*,*,0]+(One_minus_cosi_arr[1:*,0:-2,i])^1.d0*atdisk[*,*,1]+$
                           (One_minus_cosi_arr[1:*,0:-2,i])^2.d0*atdisk[*,*,2]+(One_minus_cosi_arr[1:*,0:-2,i])^3.d0*atdisk[*,*,3]+$
                           (One_minus_cosi_arr[1:*,0:-2,i])^4.d0*atdisk[*,*,4]+(One_minus_cosi_arr[1:*,0:-2,i])^5.d0*atdisk[*,*,5]
	mlambda_b[1:*,0:-2,i]=(One_minus_cosi_arr[1:*,0:-2,i])^0.d0*abulge[*,*,0]+(One_minus_cosi_arr[1:*,0:-2,i])^1.d0*abulge[*,*,1]+$
                          (One_minus_cosi_arr[1:*,0:-2,i])^2.d0*abulge[*,*,2]+(One_minus_cosi_arr[1:*,0:-2,i])^3.d0*abulge[*,*,3]+$
                          (One_minus_cosi_arr[1:*,0:-2,i])^4.d0*abulge[*,*,4]+(One_minus_cosi_arr[1:*,0:-2,i])^5.d0*abulge[*,*,5]
endfor

;calculate delta_mlambda for entire galaxy using Equation 4 given in Doore et al. (2021)
; Equation is the intrinsic version of Equation 16 from Tuffs et al (2004)
rbulge0_mat=rebin(reform(rbulge0_arr,1,1,nn),ntau,nwavebands,nn)
rdisk0_mat=rebin(reform(rdisk0_arr,1,1,nn),ntau,nwavebands,nn)

mlambda_b_mat=mlambda_b
att_bulge=rbulge0_mat*10.d0^(-1.d0*temporary(mlambda_b_mat)/2.5d0)
;There is no attenuation from the bulge in the UV
att_bulge[*,0:8]=0.d0


mlambda_d_mat=mlambda_d
att_disk=rdisk0_mat*10.d0^(-1.d0*temporary(mlambda_d_mat)/2.5d0)
;There is no attenuation from the disk in the UV
att_disk[*,0:8]=0.d0

att_disk_bulge=temporary(att_disk)+temporary(att_bulge)


mlambda_td_mat=mlambda_td
att_tdisk_1=(1-temporary(rdisk0_mat)-temporary(rbulge0_mat))*10.d0^(-1.d0*temporary(mlambda_td_mat)/2.5d0)

F_mat=rebin(reform(F_arr,1,1,nn),ntau,nwavebands,nn)
flam_mat=rebin(reform(flambda,1,nwavebands,1),ntau,nwavebands,nn)
att_tdisk=(1-(temporary(F_mat)*temporary(flam_mat)))*temporary(att_tdisk_1)


mlambda=-2.5d0*alog10(temporary(att_tdisk)+temporary(att_disk_bulge))
;There is no attenuation at 5um or taub_f=0
mlambda[0,*]=0.d0
mlambda[*,-1]=0.d0
;Infinities occur if rdisk0 or rbulge0=1, due to no attenuation in UV. Set to 0
mlambda[where(finite(mlambda) eq 0,/null)]=0.d0


;interpolate data for input tauB_f's and input wavelength's
delta_m_tau=dblarr(nwavebands)
delta_m=dblarr(nwave,nn)
for i=0,(nn-1) do begin $
  for m=0,(nwavebands-1) do delta_m_tau[m]=interpol(mlambda[*,m,i],tau,tauB_f_arr[i])
  delta_m[*,i]=interpol(delta_m_tau,wavebands,(wave*10000.d0))
  ;Set attenuation to 0 at wavelengths greater than 5.0 microns. Do this in case interpolation has error.
  delta_m[where(wave gt 5.d0 or wave lt 9.12d-2,/null),i]=0.0  
endfor

;F = F_0 * exp(-tau) = F_0 * 10^(-0.4*delta_m)
;e^(-tau) = 10^(-0.4*delta_m)

exp_tau=10.d0^(-0.4d0*delta_m)

;output form is [wavelength,tauB_f,rdisk,F,cosi,rbulge]
return,exp_tau

end



function lightning_models_vector,vectors=vectors,steps=steps,Tuffs_attenuation=Tuffs_attenuation,exp_tau=exp_tau, $
                                Filters=Filters,print_time=print_time,time=time,L_star_abs_model_table=L_star_abs_model_table, $
			                    L_star_abs_table=L_star_abs_table,BC_exp_tau=BC_exp_tau,rold0_ages=rold0_ages,$
                                _EXTRA=_e_models


;Modification History
; Keith Doore     - Replaced if statements with n_elements() on keywords to use keyword_set() - April 27th 2020
; Keith Doore     - Added use of birth cloud component on modified Calzetti if grid not set - April 27th 2020
; Keith Doore     - Added ability to input birth cloud component exp_tau in modified Calzetti - April 27th 2020
; Keith Doore     - Fixed issue with L_star_abs_table if a data point was at the final point in the table - April 27th 2020
; Keith Doore     - Added ability to use pure Calzetti curve - May 6th 2020
; Keith Doore     - Added _REF_EXTRA for keyword inheritance to cut down on list of keywords - May 6th 2020
; Keith Doore     - modified Tuffs attenuation to have rdisk as intrinsic property and rbulge as B/D - April 13th 2021
; Keith Doore     - added rold0_ages, is 1 or 0 for the corresponding age bins - April 13th 2021



if ~keyword_set(Tuffs_attenuation) and keyword_set(L_star_abs_model_table) then message,'L_star_abs_model_table keyword can only be called with Tuffs_attenuation'
if keyword_set(Tuffs_attenuation) and ~keyword_set(L_star_abs_model_table) then message,'Need LTIR table to run vectorized'

t0 = SYSTIME(/seconds)

if (N_ELEMENTS(steps) eq 0) then begin
  steps = STEPS_STELLAR(filter_labels=filter_labels,steps_bounds=steps_bounds,$
                        _EXTRA=_e_models)
endif
Nfilters = n_elements(steps.FILTER_LABELS)
Nsteps = steps.bounds.length - 1
wave_rest = steps.wave_rest ; restframe wavelength
wave_obs  = steps.wave_obs  ; observed wavelength
Nwave = wave_rest.length
nu_rest = 1.d4 * !cv.clight / wave_rest
nu_obs  = 1.d4 * !cv.clight / wave_obs
z_shift = steps.z_shift
;normalizing filters
if (n_elements(Filters) eq 0) then begin
  Filters = steps.filters
  for k=0,(Nfilters-1) do Filters[k,*] /= (trapez(Filters[k,*],nu_obs))[0]
endif

; check that rold0_ages is of correct format
if n_elements(rold0_ages) ne 0 then begin
  if n_elements(rold0_ages) ne nsteps then $
    message,'rold0_ages must have same number of elements as number of age bins'
  if n_elements(where(rold0_ages eq 1.d0 or rold0_ages eq 0.d0,/null)) ne nsteps then $
    message,'rold0_ages must have values of either 1 or 0'
endif else begin
  rold0_ages=replicate(-1.0d0,nsteps)
endelse


if keyword_set(Tuffs_attenuation) then begin
  tauB_f_vectors = vectors.tauB_f_vectors & n1 = n_elements(tauB_f_vectors)
  rold0_vectors  = vectors.rold0_vectors  & n2 = n_elements(rold0_vectors)
  F_vectors      = vectors.F_vectors      & n3 = n_elements(F_vectors)
  cosi_vectors   = vectors.cosi_vectors   & n4 = n_elements(cosi_vectors)
  b_to_d_vectors = vectors.b_to_d_vectors & n5 = n_elements(b_to_d_vectors)
  if n_elements(exp_tau) eq 0 and total(rold0_ages) lt 0 then begin
    exp_tau=Tuffs_matrix_unred_vector(wave_rest,               $
                               tauB_f_arr = tauB_f_vectors,$
                               rold0_arr = rold0_vectors, $
                               b_to_d_arr = b_to_d_vectors,$
                               F_arr=F_vectors,            $
                               cosi_arr=cosi_vectors,      $
                               _EXTRA=_e_models)
  endif else if n_elements(exp_tau) eq 0 and total(rold0_ages) ge 0 then begin
    rold0_y_vectors = replicate(0.d0,n1)
    exp_tau_young=Tuffs_matrix_unred_vector(wave_rest,               $
                               tauB_f_arr = tauB_f_vectors,$
                               rold0_arr = rold0_y_vectors, $
                               b_to_d_arr = b_to_d_vectors,$
                               F_arr=F_vectors,            $
                               cosi_arr=cosi_vectors,      $
                               _EXTRA=_e_models)
    rold0_o_vectors = replicate(1.d0,n1)
    exp_tau_old=Tuffs_matrix_unred_vector(wave_rest,               $
                               tauB_f_arr = tauB_f_vectors,$
                               rold0_arr = rold0_o_vectors, $
                               b_to_d_arr = b_to_d_vectors,$
                               F_arr=F_vectors,            $
                               cosi_arr=cosi_vectors,      $
                               _EXTRA=_e_models)
  endif
endif else begin
  tauV_DIFF_vectors = vectors.tauV_DIFF_vectors & n1 = n_elements(tauV_DIFF_vectors)
  delta_vectors     = vectors.delta_vectors     & n2 = n_elements(delta_vectors)
  tauV_BC_vectors   = vectors.tauV_BC_vectors   & n3 = n_elements(tauV_BC_vectors)
  n4 = n1
  n5 = n1
  if (n_elements(BC_exp_tau) eq 0) then $      
    BC_exp_tau = matrix_unred_Keith_vector(wave_rest,tauV_DIFF_arr=-1.d0*tauV_DIFF_vectors,delta_arr=delta_vectors,$
                 tauV_BC_arr=-1.d0*tauV_BC_vectors,_EXTRA=_e_models)
  if (n_elements(exp_tau) eq 0) then $     
    exp_tau    = matrix_unred_Keith_vector(wave_rest,tauV_DIFF_arr=-1.d0*tauV_DIFF_vectors,delta_arr=delta_vectors,$
                 tauV_BC_arr=dblarr(n3),_EXTRA=_e_models)
endelse

if n1 ne n2 or n1 ne n3 or n1 ne n4 or n1 ne n5 or n2 ne n3 or n2 ne n4 or n2 ne n5 or $
n3 ne n4 or n3 ne n5 or n4 ne n5 then message,'Input parameters must have same length for vectorization'
nn=n1


steps_mean_Lnu_red_grid = dblarr(Nfilters+1,Nsteps,nn)
transmission = transpose(Filters)
Lnu_inc = dblarr(Nsteps,nn)

if keyword_set(L_star_abs_model_table) then begin
  if n_elements(L_star_abs_table) eq 0 then $
    L_star_abs_table=mrdfits(_e_models.lightning_folder+'L_star_abs_model_table/L_star_abs_model_table_z_'+string(z_shift,'(F4.2)')+'.fits',1)

  n1_L_star_abs_table = n_elements(L_star_abs_table.tauB_f_grid)
  n2_L_star_abs_table = n_elements(L_star_abs_table.rdisk0_grid)
  n3_L_star_abs_table = n_elements(L_star_abs_table.F_grid)
  n5_L_star_abs_table = n_elements(L_star_abs_table.rbulge0_grid)

  ;Needed to allow for below process to work if n5_L_star_abs_table=1, which it currently is
  if n5_L_star_abs_table eq 1 then LABS_STAR_MODEL=rebin(L_star_abs_table.LABS_STAR_MODEL,L_star_abs_table.nsteps,n1_L_star_abs_table,n2_L_star_abs_table,n3_L_star_abs_table,2)

  taub_loc    = interpol(lindgen(n1_L_star_abs_table),L_star_abs_table.tauB_f_grid,tauB_f_vectors)>0
  F_loc       = interpol(lindgen(n3_L_star_abs_table),L_star_abs_table.F_grid,F_vectors)>0
  tau_weights=taub_loc-fix(taub_loc)
  f_weights=f_loc-fix(f_loc)

  if total(rold0_ages) lt 0 then begin
    rdisk0_vectors=rold0_vectors/(1+b_to_d_vectors)
    rbulge0_vectors=rold0_vectors*b_to_d_vectors/(1+b_to_d_vectors)

    rdisk0_loc  = interpol(lindgen(n2_L_star_abs_table),L_star_abs_table.rdisk0_grid,rdisk0_vectors)>0
    rbulge0_loc = interpol(lindgen(n5_L_star_abs_table),L_star_abs_table.rbulge0_grid,rbulge0_vectors)>0
    
    rdisk0_weights=rdisk0_loc-fix(rdisk0_loc)
    rbulge0_weights=rbulge0_loc-fix(rbulge0_loc)
    
    for k1=0,(nn-1) do begin
      if taub_loc[k1] eq (n1_L_star_abs_table-1) then taub_loc_temp=[taub_loc[k1],taub_loc[k1]] $
        else taub_loc_temp=[taub_loc[k1]:(taub_loc[k1]+1)]
      if rdisk0_loc[k1] eq (n2_L_star_abs_table-1) then rdisk0_loc_temp=[rdisk0_loc[k1],rdisk0_loc[k1]] $
        else rdisk0_loc_temp=[rdisk0_loc[k1]:(rdisk0_loc[k1]+1)]
      if f_loc[k1] eq (n3_L_star_abs_table-1) then f_loc_temp=[f_loc[k1],f_loc[k1]] $
        else f_loc_temp=[f_loc[k1]:(f_loc[k1]+1)] 
      if rbulge0_loc[k1] eq (n5_L_star_abs_table-1) then rbulge0_loc_temp=[rbulge0_loc[k1],rbulge0_loc[k1]] $
        else rbulge0_loc_temp=[rbulge0_loc[k1]:(rbulge0_loc[k1]+1)] 
        
      averaging_region=reform(LABS_STAR_MODEL[*,taub_loc_temp,rdisk0_loc_temp,f_loc_temp,rbulge0_loc_temp])
      temp1=reform(averaging_region[*,1,*,*,*]*tau_weights[k1]+averaging_region[*,0,*,*,*]*(1-tau_weights[k1]))
      temp2=reform(temp1[*,1,*,*]*rdisk0_weights[k1]+temp1[*,0,*,*]*(1-rdisk0_weights[k1]))
      temp3=reform(temp2[*,1,*]*f_weights[k1]+temp2[*,0,*]*(1-f_weights[k1]))
      Lnu_inc[*,k1]=reform(temp3[*,1]*rbulge0_weights[k1]+temp3[*,0]*(1-rbulge0_weights[k1]))
    endfor
  endif else begin
    rdisk0_vectors=rold0_o_vectors/(1+b_to_d_vectors)
    rbulge0_vectors=rold0_o_vectors*b_to_d_vectors/(1+b_to_d_vectors)

    rdisk0_o_loc  = interpol(lindgen(n2_L_star_abs_table),L_star_abs_table.rdisk0_grid,rdisk0_vectors)>0
    rbulge0_o_loc = interpol(lindgen(n5_L_star_abs_table),L_star_abs_table.rbulge0_grid,rbulge0_vectors)>0
    
    rdisk0_o_weights=rdisk0_o_loc-fix(rdisk0_o_loc)
    rbulge0_o_weights=rbulge0_o_loc-fix(rbulge0_o_loc)
   
    for k1=0,(nn-1) do begin
      if taub_loc[k1] eq (n1_L_star_abs_table-1) then taub_loc_temp=[taub_loc[k1],taub_loc[k1]] $
        else taub_loc_temp=[taub_loc[k1]:(taub_loc[k1]+1)]
      if rdisk0_o_loc[k1] eq (n2_L_star_abs_table-1) then rdisk0_o_loc_temp=[rdisk0_o_loc[k1],rdisk0_o_loc[k1]] $
        else rdisk0_o_loc_temp=[rdisk0_o_loc[k1]:(rdisk0_o_loc[k1]+1)]
      if f_loc[k1] eq (n3_L_star_abs_table-1) then f_loc_temp=[f_loc[k1],f_loc[k1]] $
        else f_loc_temp=[f_loc[k1]:(f_loc[k1]+1)] 
      if rbulge0_o_loc[k1] eq (n5_L_star_abs_table-1) then rbulge0_o_loc_temp=[rbulge0_o_loc[k1],rbulge0_o_loc[k1]] $
        else rbulge0_o_loc_temp=[rbulge0_o_loc[k1]:(rbulge0_o_loc[k1]+1)] 
        
      averaging_region=LABS_STAR_MODEL[where(rold0_ages eq 1,/null),*,*,*,*]
      averaging_region=averaging_region[*,taub_loc_temp,rdisk0_o_loc_temp,f_loc_temp,rbulge0_o_loc_temp]
      temp1=averaging_region[*,1,*,*,*]*tau_weights[k1]+averaging_region[*,0,*,*,*]*(1-tau_weights[k1])
      temp2=temp1[*,*,1,*,*]*rdisk0_o_weights[k1]+temp1[*,*,0,*,*]*(1-rdisk0_o_weights[k1])
      temp3=temp2[*,*,*,1,*]*f_weights[k1]+temp2[*,*,*,0,*]*(1-f_weights[k1])
      Lnu_inc[where(rold0_ages eq 1,/null),k1]=reform(temp3[*,*,*,*,1]*rbulge0_o_weights[k1]+temp3[*,*,*,*,0]*(1-rbulge0_o_weights[k1]))

      averaging_region=LABS_STAR_MODEL[where(rold0_ages eq 0,/null),*,0,*,0]
      averaging_region=averaging_region[*,taub_loc_temp,*,f_loc_temp,*]
      temp1=averaging_region[*,1,*,*,*]*tau_weights[k1]+averaging_region[*,0,*,*,*]*(1-tau_weights[k1])
      Lnu_inc[where(rold0_ages eq 0,/null),k1]=reform(temp1[*,*,*,1,*]*f_weights[k1]+temp1[*,*,*,0,*]*(1-f_weights[k1]))
    endfor
  endelse
endif

if ~keyword_set(L_star_abs_model_table) then Lnu_att_total = dblarr(Nsteps,nn)
for st=0,(Nsteps-1) do begin
  Lnu_5D=rebin(steps.lnu[*,st],nwave,nn)
  if ~keyword_set(Tuffs_attenuation) and st eq 0 then begin
    steps_Lnu_red = BC_exp_tau*Lnu_5D
  endif else if keyword_set(Tuffs_attenuation) and rold0_ages[st] eq 0 then begin
    steps_Lnu_red = exp_tau_young*Lnu_5D
  endif else if keyword_set(Tuffs_attenuation) and rold0_ages[st] eq 1 then begin
    steps_Lnu_red = exp_tau_old*Lnu_5D
  endif else begin
    steps_Lnu_red = exp_tau*Lnu_5D
  endelse
  
  ;Lbol = (1+z_shift)*trapez(lnu_5d,nu_rest)*(-1.d0) WRONG 1+z factor
  ;Lbol = trapez(lnu_5d,nu_rest)/(1+z_shift)*(-1.d0) ; correct 1+z factor, does not include nebular emission absorption
  ;Lbol = trapez(lnu_5d,nu_obs)*(-1.d0)              ; same as previous line
  Lbol = (steps.Lbol[st])[0]
    for k1=0,(nn-1) do begin
      for k=0,(Nfilters-1) do begin
        steps_mean_Lnu_red_grid[k,st,k1] = trapez(steps_Lnu_red[*,k1]*transmission[*,k],nu_obs)
      endfor
      if ~keyword_set(L_star_abs_model_table) then $
    	Lnu_att_total[st,k1] = trapez(steps_Lnu_red[*,k1],nu_obs)
      ; gives negative value from integral bc nu_rest is inverted           
      if keyword_set(Tuffs_attenuation) then begin
        steps_mean_Lnu_red_grid[Nfilters,st,k1] = Lnu_inc[st,k1]
      endif
    endfor 
  if ~keyword_set(Tuffs_attenuation) then begin
    steps_mean_Lnu_red_grid[Nfilters,st,*] = Lbol + Lnu_att_total[st,*]
  endif
endfor

steps_Lnu_red = !null

t2 = systime(/seconds)
if keyword_set(print_time) then print,'Total models running time: ',t2-t0,' seconds'
time=t2-t0
models = temporary(steps_mean_Lnu_red_grid)
return,models

end



function lightning_priors,alpha=alpha,Umin=Umin,Umax=Umax,gamma=gamma,qPAH=qPAH,tauv_diff=tauv_diff,$
                     delta=delta,tauv_bc=tauv_bc,taub_f=taub_f,rdisk=rdisk,F=F,cosi=cosi,rbulge=rbulge,$
                     Tuffs_attenuation=Tuffs_attenuation,prior_dist=prior_dist
                    
; alpha, umin, umax, etc are the input parameters from the MCMC.
; Prior_dist is a structure that has the tabulated pdf and parameter sample space (ss) spanning max parameter range
;   i.e., prior_dist.cosi_pdf and prior_dist.cosi_ss
;   In general would be prior_dist.x_pdf and prior_dist.x_ss where x is one of our parameters
; If the PDF was a function this could be added or just made into tabulated version and added to prior_dist
;   I think making the function tabulated would be the best method for consistency
; Assume we just multiply all the prior density of each parameter together

;MCMC questions
; Should we have a prior keyword that initiates the priors, and if not set then the prior ratio is 1

prior_density_cosi=interpol(prior_dist.inc_pdf,prior_dist.inc_bins,acos(cosi)*180/!pi)

prior_density=prior_density_cosi > 1d-140

return,prior_density

end







pro lightning_MCMC_vector,Lnu_obs,Lnu_unc,LTIR_obs,LTIR_unc,galaxy_id,filter_labels,z_shift=z_shift,steps_bounds=steps_bounds,$
                   print_time=print_time,Tuffs_attenuation=Tuffs_attenuation,Nparallel=Nparallel,$
                   par1_constant=par1_constant,par2_constant=par2_constant,par3_constant=par3_constant,$
                   par4_constant=par4_constant,par5_constant=par5_constant,par6_constant=par6_constant,$
                   par7_constant=par7_constant,par8_constant=par8_constant,par9_constant=par9_constant,$
                   par10_constant=par10_constant,Ntrials=Ntrials,outfolder=outfolder,beta_exponent=beta_exponent,$
                   mu_input=mu_input,sigma_input=sigma_input,lambda_input=lambda_input,$
                   adaptive=adaptive,coeff_start=coeff_start,coeff_sigma=coeff_sigma,use_priors=use_priors,$
                   parameter_start=parameter_start,parameter_sigma=parameter_sigma,dust_emission=dust_emission,_REF_EXTRA=_extra_MCMC
                   
;Required inputs
;  Lnu_obs [Nfilters] or Lnu_obs [Nfilters,Ngal]
;  Lnu_unc [Nfilters] or Lnu_unc [Nfilters,Ngal]
;  galaxy_ID[Ngal] type string 
;  filter_labels[Nfilters] type string preespecified from list
;
;Optional Input
;  upper_limit_flag[Nfilters,Ngal] either zero or 1    => Lmod > Lnu_obs

;Modification History
; Keith Doore     - Replaced if statements with n_elements() on keywords to use keyword_set() - April 27th 2020
; Keith Doore     - Set default starting attenuation parameters to 0.0 so if not wanted they would not be used - May 6th 2020
; Keith Doore     - Added ability to use pure Calzetti curve (no changes here due to _extra) - May 6th 2020
; Keith Doore     - Added _REF_EXTRA for keyword inheritance to cut down on list of keywords - May 6th 2020
; Keith Doore     - Changed all keywords that were the same to match across functions/procedures - May 6th 2020
; Keith Doore     - Removed any repetitive items at beginning that are set by other functions if not set in MCMC call - May 6th 2020
; Keith Doore     - Added needed items to run Tuffs attenuation as to match other attenuation - May 6th 2020
; Keith Doore     - modified Tuffs attenuation to have rdisk as intrinsic property and rbulge as B/D - April 13th 2021




  t0=systime(/sec)
  
  if typename(galaxy_id) ne 'STRING' then message,'ERROR -- galaxy_id input must be of type STRING'
  if (Lnu_obs.length ne filter_labels.length) then message,'ERROR -- SED and Filter_labels must have same sizes'
  if (Lnu_obs.length ne Lnu_unc.length) then message,'ERROR -- SED and uncertainties must have same sizes'
  if (LTIR_obs.length ne LTIR_unc.length) then message,'ERROR -- TIR and uncertainties must have same sizes'
  
  if (n_elements(z_shift) eq 0) then z_shift = 0.0
  if (n_elements(steps_bounds) eq 0) then steps_bounds = [0.d0,1.d7,1.d8,1.d9,5.d9,13.6d9]
  if (n_elements(Ntrials) eq 0) then Ntrials = 1e5
  if (n_elements(outfolder) eq 0) then outfolder_plots='~/MCMC_Lightning_runs/'
  if (n_elements(beta_exponent) eq 0) then beta_exponent = 0.35
  if (n_elements(Nparallel) eq 0) then Nparallel = 1
  
  if (n_elements(parameter_start) eq 0) and keyword_set(dust_emission) and ~keyword_set(Tuffs_attenuation) then begin
    parameter_start=[0.2,0.0,0.2,0,0,-2.d0,1.0,1.d4,0.1,0.020]
    parameter_start=rebin(parameter_start,n_elements(parameter_start),Nparallel)
  endif
  if (n_elements(parameter_sigma) eq 0) and keyword_set(dust_emission) and ~keyword_set(Tuffs_attenuation) then begin
     parameter_sigma=[0.1,0.1,0.2,0,0,0.5d0,0.1,1.d4,0.1,0.005]
    parameter_sigma=rebin(parameter_sigma,n_elements(parameter_sigma),Nparallel)
  endif
  if (n_elements(parameter_start) eq 0) and keyword_set(dust_emission) and keyword_set(Tuffs_attenuation) then begin
    parameter_start=[1.0,0.0,0.2,0.8,0.0,-2.d0,1.0,1.d4,0.1,0.020]
    parameter_start=rebin(parameter_start,n_elements(parameter_start),Nparallel)
  endif
  if (n_elements(parameter_sigma) eq 0) and keyword_set(dust_emission) and keyword_set(Tuffs_attenuation) then begin
    parameter_sigma=[0.5,0.1,0.1,0.1,0.1,0.5d0,0.1,1.d4,0.1,0.005]
    parameter_sigma=rebin(parameter_sigma,n_elements(parameter_sigma),Nparallel)
  endif

  if Nparallel eq 1 and (size(parameter_start))[0] ne 1 and (size(parameter_sigma))[0] ne 1 then $
    message,'ERROR -- parameter_start and parameter_sigma must have one dimension'
  if Nparallel gt 1 and (size(parameter_start))[2] ne Nparallel and (size(parameter_sigma))[2] ne Nparallel then $
    message,'ERROR -- parameter_start and parameter_sigma must have second dimension equal to Nparallel'
 
  
  steps_bounds = steps_bounds < galage(z_shift,25,/silent)
  steps_bounds = steps_bounds[unique(steps_bounds)]
  Nsteps = n_elements(steps_bounds) - 1
  steps = STEPS_STELLAR(filter_labels=filter_labels,steps_bounds=steps_bounds,$
                        z_shift=z_shift,_EXTRA=_extra_MCMC)
                        

  Nfilters = filter_labels.length
  Filters = steps.filters
  for k=0,(Nfilters-1) do Filters[k,*] /= (trapez(Filters[k,*],(1.d4*!cv.clight/steps.wave_obs)))[0]

  Lmod_obs = [Lnu_obs,LTIR_obs]
  Lmod_unc = [Lnu_unc,LTIR_unc]
  Lmod_obs = rebin(Lmod_obs,Nfilters+1,Nparallel)
  Lmod_unc = rebin(Lmod_unc,Nfilters+1,Nparallel)

  if keyword_set(Tuffs_attenuation) then begin
    ran1=[ 0.d,8.0d]
    ran2=[ 0.d,1.0d]
    ran3=[0.d,0.61d]
    ran4=[ 0.d,1.0d]
    ran5=[ 0.d,!values.d_infinity]
  
    par1_guess = parameter_start[0,*]  &  par1_sigma = parameter_sigma[0,*]
    par2_guess = parameter_start[1,*]  &  par2_sigma = parameter_sigma[1,*]
    par3_guess = parameter_start[2,*]  &  par3_sigma = parameter_sigma[2,*]
    par4_guess = parameter_start[3,*]  &  par4_sigma = parameter_sigma[3,*]
    par5_guess = parameter_start[4,*]  &  par5_sigma = parameter_sigma[4,*]
    if keyword_set(par1_constant) then par1_sigma=rebin(reform(0.d0,1,1),1,Nparallel)
    if keyword_set(par2_constant) then par2_sigma=rebin(reform(0.d0,1,1),1,Nparallel)
    if keyword_set(par3_constant) then par3_sigma=rebin(reform(0.d0,1,1),1,Nparallel)
    if keyword_set(par4_constant) then par4_sigma=rebin(reform(0.d0,1,1),1,Nparallel)
    if keyword_set(par5_constant) then par5_sigma=rebin(reform(0.d0,1,1),1,Nparallel)

    vectors = {tauB_f_vectors: reform(par1_guess), rold0_vectors: reform(par2_guess), F_vectors: reform(par3_guess), cosi_vectors: reform(par4_guess), b_to_d_vectors: reform(par5_guess)}
    models=lightning_models_vector(steps=steps,vectors=vectors,Filters=Filters,/tuffs,_EXTRA=_extra_MCMC)
    models[where(finite(models,/NaN),/null)] = 0.0
  endif else begin
    ran1=[ 0.0d,3.0d]
    ran2=[-2.3d,0.4d]
    ran3=[ 0.0d,4.0d]
    ran4=[ 0.0d,0.0d]
    ran5=[ 0.0d,0.0d]

    par1_guess = parameter_start[0,*] & par1_sigma = parameter_sigma[0,*]
    par2_guess = parameter_start[1,*] & par2_sigma = parameter_sigma[1,*]
    par3_guess = parameter_start[2,*] & par3_sigma = parameter_sigma[2,*]
    if keyword_set(par1_constant) then par1_sigma=rebin(reform(0.d0,1,1),1,Nparallel)
    if keyword_set(par2_constant) then par2_sigma=rebin(reform(0.d0,1,1),1,Nparallel)
    if keyword_set(par3_constant) then par3_sigma=rebin(reform(0.d0,1,1),1,Nparallel)
    par4_guess = rebin(reform(0.d0,1,1),1,Nparallel)  &  par4_sigma = rebin(reform(0.d0,1,1),1,Nparallel)  &  par4_constant=1                             
    par5_guess = rebin(reform(0.d0,1,1),1,Nparallel)  &  par5_sigma = rebin(reform(0.d0,1,1),1,Nparallel)  &  par5_constant=1   

    vectors   = {tauV_DIFF_vectors: reform(par1_guess), delta_vectors: reform(par2_guess), tauV_BC_vectors: reform(par3_guess)}
    models = lightning_models_vector(steps=steps,vectors=vectors,Filters=Filters,_EXTRA=_extra_MCMC)
    models[where(finite(models,/NaN),/null)] = 0.0
  endelse


  if keyword_set(dust_emission) then begin
    ran6 = [-10.d0,  4d0]       ; alpha
    ran7 = [ 0.7d0,25]          ; Umin
    ran8 = [ 1.d3,  3.d5]       ; Umax
    ran9 = [ 0.d0,  1d0]        ; gamma
    ran10= [ 0.0047,0.0458d0]  ; qPAH

    par6_guess = parameter_start[5,*] & par6_sigma = parameter_sigma[5,*] & if keyword_set(par6_constant)  then par6_sigma=rebin(reform(0.d0,1,1),1,Nparallel)
    par7_guess = parameter_start[6,*] & par7_sigma = parameter_sigma[6,*] & if keyword_set(par7_constant)  then par7_sigma=rebin(reform(0.d0,1,1),1,Nparallel)
    par8_guess = parameter_start[7,*] & par8_sigma = parameter_sigma[7,*] & if keyword_set(par8_constant)  then par8_sigma=rebin(reform(0.d0,1,1),1,Nparallel)
    par9_guess = parameter_start[8,*] & par9_sigma = parameter_sigma[8,*] & if keyword_set(par9_constant)  then par9_sigma=rebin(reform(0.d0,1,1),1,Nparallel)
    par10_guess= parameter_start[9,*] & par10_sigma= parameter_sigma[9,*] & if keyword_set(par10_constant) then par10_sigma=rebin(reform(0.d0,1,1),1,Nparallel)
    
    dl07     = dl07_templates(filter_labels=filter_labels,z_shift=z_shift,_EXTRA=_extra_MCMC)
    de_model = dl07_sed_vector(dl07,alpha=reform(par6_guess),umin=reform(par7_guess),umax=reform(par8_guess),$
               gam=reform(par9_guess),qPAH=reform(par10_guess),filter_labels=filter_labels,Lbol=de_LTIR)
  endif


  if (n_elements(coeff_start) eq 0) then coeff_start = replicate(1.d0,nsteps,Nparallel)
  if ((size(coeff_start))[1] ne nsteps) then message, 'ERROR -- Number of starting SFH coefficients must equal number of steps'
  if (n_elements(coeff_sigma) eq 0) then coeff_sigma = replicate(1.d0,nsteps,Nparallel)
  if ((size(coeff_sigma))[1] ne nsteps) then message, 'ERROR -- Number of starting SFH coefficients must equal number of steps'
  if Nparallel eq 1 and (size(coeff_start))[0] ne 1 and (size(coeff_sigma))[0] ne 1 then $
    message,'ERROR -- coeff_start and coeff_sigma must have one dimension'
  if Nparallel gt 1 and (size(coeff_start))[2] ne Nparallel and (size(coeff_sigma))[2] ne Nparallel then $
    message,'ERROR -- coeff_start and coeff_sigma must have second dimension equal to Nparallel'

  if keyword_set(par1_constant) then random2  = 0       else random2  = 1
  if keyword_set(par2_constant) then random3  = random2 else random3  = 1 + random2
  if keyword_set(par3_constant) then random4  = random3 else random4  = 1 + random3
  if keyword_set(par4_constant) then random5  = random4 else random5  = 1 + random4
  if keyword_set(par5_constant) then random6  = random5 else random6  = 1 + random5
  if keyword_set(par6_constant) then random7  = random6 else random7  = 1 + random6
  if keyword_set(par7_constant) then random8  = random7 else random8  = 1 + random7
  if keyword_set(par8_constant) then random9  = random8 else random9  = 1 + random8
  if keyword_set(par9_constant) then random10 = random9 else random10 = 1 + random9


  coeff_guess = coeff_start
  
  coeff_old = coeff_guess                     
  par1_old = par1_guess
  par2_old = par2_guess
  par3_old = par3_guess
  par4_old = par4_guess
  par5_old = par5_guess
  par6_old = par6_guess
  par7_old = par7_guess
  par8_old = par8_guess
  par9_old = par9_guess
  par10_old= par10_guess
  Lmod_old  = total(models * (rebin(reform(coeff_old,1,Nsteps,Nparallel),Nfilters+1,Nsteps,Nparallel)),2) ;models(Nfilters+1,Nsteps,Nparallel)
  Lmod_old[0:(Nfilters-1),*] = Lmod_old[0:(Nfilters-1),*] + $
                               rebin(reform((reform(Lmod_old[Nfilters,*])/de_LTIR),1,Nparallel),Nfilters,Nparallel)*de_model
  chi2_old = total(((Lmod_old - Lmod_obs)/Lmod_unc)^2,/NaN,1)

  chain_old  = [coeff_old,par1_old,par2_old,par3_old,par4_old,par5_old,par6_old,par7_old,par8_old,par9_old,par10_old]
  sigma_vals = [coeff_sigma,par1_sigma,par2_sigma,par3_sigma,par4_sigma,par5_sigma,par6_sigma,par7_sigma,par8_sigma,par9_sigma,par10_sigma]

  if keyword_set(use_priors) then chi2_prior_old=-2.d0*alog(lightning_priors(cosi=reform(chain_old[Nsteps+3,*]),_EXTRA=_extra_MCMC)) else chi2_prior_old=dblarr(Nparallel)
  
  chain = dblarr(nsteps+10,Ntrials,Nparallel)
  chi2_chain = dblarr(Ntrials,Nparallel)
  chi2_prior_chain = dblarr(Ntrials,Nparallel)
  chain[*,0,*] = chain_old
  chi2_chain[0,*] = chi2_old
  chi2_prior_chain[0,*] = chi2_prior_old
  
  constant_param = where(sigma_vals[*,0] eq 0,ncomp=Npar,comp=variable_param,/null)

  if (n_elements(lambda_input) eq 0) then lambda_new=2.38^2.d0/replicate(Npar,Nparallel) else lambda_new=lambda_input
  ;if (n_elements(mu_input) eq 0) then mu_new=chain else mu_new=mu_input
  if (n_elements(mu_input) eq 0) then mu_new=chain_old+sigma_vals else mu_new=mu_input
  if (n_elements(sigma_input) eq 0) then begin
    sigma_new=dblarr(Npar,Npar,Nparallel)
    for i=0,(nparallel-1) do sigma_new[*,*,i]=diag_matrix(sigma_vals[variable_param,i]) 
  endif else begin
    sigma_new=sigma_input
  endelse

  p_jump=[0.441,0.352,0.316,0.279,0.275,0.266,0.261,0.255,0.261,0.267]
  alpha_star=p_jump[(Npar-1) < (p_jump.length-1)]
  accepted_trials=fltarr(Nparallel)
  t0=systime(/sec)

  for trial=0,(Ntrials-2) do begin
    gamma_new  = 1.d0/(trial+1.d0)^beta_exponent
    lambda_old = lambda_new
    mu_old     = mu_new
    sigma_old  = sigma_new
    covar      = rebin(reform(lambda_old,1,1,Nparallel),Npar,Npar,Nparallel)*sigma_old
    chain_new = chain_old
    for njump=0,(Nparallel-1) do begin
      jump       = mrandomn(seed,covar[*,*,njump],status=status)
      while (status gt 0) do begin
        evals = eigenql(covar[*,*,njump],EIGENVECTORS = evecs)
        evals[where(evals le 1.0d-16,/null)] = 1.0d-16
        covar[*,*,njump] = evecs # diag_matrix(evals) # invert(evecs)
        covar[*,*,njump] = (covar[*,*,njump] + transpose(covar[*,*,njump]))/2.d0
        jump  = mrandomn(seed,covar[*,*,njump],status=status)
      endwhile
      chain_new[variable_param,njump] = chain_old[variable_param,njump] + jump
      in = 1
      coeff_new = chain_new[0:(Nsteps-1),njump] & for st=0,(nsteps-1) do in = in and (coeff_new[st] ge 0)
      par1_new  = chain_new[Nsteps  ,njump] & in = in and ((par1_new  ge min(ran1 )) and (par1_new  le max(ran1 )))
      par2_new  = chain_new[Nsteps+1,njump] & in = in and ((par2_new  ge min(ran2 )) and (par2_new  le max(ran2 )))
      par3_new  = chain_new[Nsteps+2,njump] & in = in and ((par3_new  ge min(ran3 )) and (par3_new  le max(ran3 )))
      par4_new  = chain_new[Nsteps+3,njump] & in = in and ((par4_new  ge min(ran4 )) and (par4_new  le max(ran4 )))
      par5_new  = chain_new[Nsteps+4,njump] & in = in and ((par5_new  ge min(ran5 )) and (par5_new  le max(ran5 )))
      par6_new  = chain_new[Nsteps+5,njump] & in = in and ((par6_new  ge min(ran6 )) and (par6_new  le max(ran6 )))
      par7_new  = chain_new[Nsteps+6,njump] & in = in and ((par7_new  ge min(ran7 )) and (par7_new  le max(ran7 )))
      par8_new  = chain_new[Nsteps+7,njump] & in = in and ((par8_new  ge min(ran8 )) and (par8_new  le max(ran8 )))
      par9_new  = chain_new[Nsteps+8,njump] & in = in and ((par9_new  ge min(ran9 )) and (par9_new  le max(ran9 )))
      par10_new = chain_new[Nsteps+9,njump] & in = in and ((par10_new ge min(ran10)) and (par10_new le max(ran10)))
      while in eq 0 do begin
        jump  = mrandomn(seed,covar[*,*,njump],status=status)
        chain_new[variable_param,njump] = chain_old[variable_param,njump] + jump
        in = 1
        coeff_new = chain_new[0:(Nsteps-1),njump] & for st=0,(nsteps-1) do in = in and (coeff_new[st] ge 0)
        par1_new  = chain_new[Nsteps  ,njump] & in = in and ((par1_new  ge min(ran1 )) and (par1_new  le max(ran1 )))
        par2_new  = chain_new[Nsteps+1,njump] & in = in and ((par2_new  ge min(ran2 )) and (par2_new  le max(ran2 )))
        par3_new  = chain_new[Nsteps+2,njump] & in = in and ((par3_new  ge min(ran3 )) and (par3_new  le max(ran3 )))
        par4_new  = chain_new[Nsteps+3,njump] & in = in and ((par4_new  ge min(ran4 )) and (par4_new  le max(ran4 )))
        par5_new  = chain_new[Nsteps+4,njump] & in = in and ((par5_new  ge min(ran5 )) and (par5_new  le max(ran5 )))
        par6_new  = chain_new[Nsteps+5,njump] & in = in and ((par6_new  ge min(ran6 )) and (par6_new  le max(ran6 )))
        par7_new  = chain_new[Nsteps+6,njump] & in = in and ((par7_new  ge min(ran7 )) and (par7_new  le max(ran7 )))
        par8_new  = chain_new[Nsteps+7,njump] & in = in and ((par8_new  ge min(ran8 )) and (par8_new  le max(ran8 )))
        par9_new  = chain_new[Nsteps+8,njump] & in = in and ((par9_new  ge min(ran9 )) and (par9_new  le max(ran9 )))
        par10_new = chain_new[Nsteps+9,njump] & in = in and ((par10_new ge min(ran10)) and (par10_new le max(ran10)))
      endwhile
    endfor
;    coeff_new=!null
;    for i=0,(nsteps-1) do begin                   
;      coeff_buffer = chain_new[i] & while (coeff_buffer lt 0) do coeff_buffer = chain_old[i] + (mrandomn(seed,covar))[i]
;      coeff_new = [coeff_new,coeff_buffer]
;    endfor 
;    par1_new = chain_new[nsteps]   & while ((par1_new lt min(ran1)) or (par1_new gt max(ran1))) do par1_new = chain_old[nsteps]   + (mrandomn(seed,covar))[nsteps]
;    par2_new = chain_new[nsteps+1] & while ((par2_new lt min(ran2)) or (par2_new gt max(ran2))) do par2_new = chain_old[nsteps+1] + (mrandomn(seed,covar))[nsteps+random2]
;    par3_new = chain_new[nsteps+2] & while ((par3_new lt min(ran3)) or (par3_new gt max(ran3))) do par3_new = chain_old[nsteps+2] + (mrandomn(seed,covar))[nsteps+random3]
;    par4_new = chain_new[nsteps+3] & while ((par4_new lt min(ran4)) or (par4_new gt max(ran4))) do par4_new = chain_old[nsteps+3] + (mrandomn(seed,covar))[nsteps+random4]
;    par5_new = chain_new[nsteps+4] & while ((par5_new lt min(ran5)) or (par5_new gt max(ran5))) do par5_new = chain_old[nsteps+4] + (mrandomn(seed,covar))[nsteps+random5]
;    if (n_elements(dust_emission) eq 1) then begin
;      par6_new = chain_new[nsteps+5]  & while ((par6_new  lt min(ran6))  or (par6_new  gt max(ran6)))  do par6_new  = chain_old[nsteps+5] + (mrandomn(seed,covar))[nsteps+random6]
;      par7_new = chain_new[nsteps+6]  & while ((par7_new  lt min(ran7))  or (par7_new  gt max(ran7)))  do par7_new  = chain_old[nsteps+6] + (mrandomn(seed,covar))[nsteps+random7]
;      par8_new = chain_new[nsteps+7]  & while ((par8_new  lt min(ran8))  or (par8_new  gt max(ran8)))  do par8_new  = chain_old[nsteps+7] + (mrandomn(seed,covar))[nsteps+random8]
;      par9_new = chain_new[nsteps+8]  & while ((par9_new  lt min(ran9))  or (par9_new  gt max(ran9)))  do par9_new  = chain_old[nsteps+8] + (mrandomn(seed,covar))[nsteps+random9]
;      par10_new = chain_new[nsteps+9] & while ((par10_new lt min(ran10)) or (par10_new gt max(ran10))) do par10_new = chain_old[nsteps+9] + (mrandomn(seed,covar))[nsteps+random10]
;      de_model = dl07_sed(dl07,alpha=par6_new,umin=par7_new,umax=par8_new,qPAH=par10_new,filter_labels=filter_labels,Lbol=de_LTIR)
;    endif
;stop
    if keyword_set(Tuffs_attenuation) then begin
      vectors = {tauB_f_vectors: reform(chain_new[Nsteps,*]), rold0_vectors: reform(chain_new[Nsteps+1,*]), $
                 F_vectors: reform(chain_new[Nsteps+2,*]), cosi_vectors: reform(chain_new[Nsteps+3,*]), $
                 b_to_d_vectors: reform(chain_new[Nsteps+4,*])}
      models=lightning_models_vector(steps=steps,vectors=vectors,Filters=Filters,/tuffs,_EXTRA=_extra_MCMC)
      models[where(finite(models,/NaN),/null)] = 0.0
    endif else begin
      vectors  = {tauV_DIFF_vectors: reform(chain_new[Nsteps,*]), delta_vectors: reform(chain_new[Nsteps+1,*]), $
                  tauV_BC_vectors: reform(chain_new[Nsteps+2,*])}
      models = lightning_models_vector(steps=steps,vectors=vectors,Filters=Filters,_EXTRA=_extra_MCMC)
      models[where(finite(models,/NaN),/null)] = 0.0
    endelse
    if keyword_set(dust_emission) then begin
      de_model = dl07_sed_vector(dl07,alpha=reform(chain_new[Nsteps+5,*]),umin=reform(chain_new[Nsteps+6,*]),$
        umax=reform(chain_new[Nsteps+7,*]),gam=reform(chain_new[Nsteps+8,*]),qPAH=reform(chain_new[Nsteps+9,*]),$
        filter_labels=filter_labels,Lbol=de_LTIR)
    endif
        
    Lmod_new  = total(models * (rebin(reform(chain_new[0:(Nsteps-1),*],1,Nsteps,Nparallel),Nfilters+1,Nsteps,Nparallel)),2) ;models(Nfilters+1,Nsteps,Nparallel)
    Lmod_new[0:(Nfilters-1),*] = Lmod_new[0:(Nfilters-1),*] + $
                                 rebin(reform((reform(Lmod_new[Nfilters,*])/de_LTIR),1,Nparallel),Nfilters,Nparallel)*de_model
    chi2_new = total(((Lmod_new - Lmod_obs)/Lmod_unc)^2,/NaN,1)
    
    if keyword_set(use_priors) then chi2_prior_new=-2.d0*alog(lightning_priors(cosi=reform(chain_new[Nsteps+3,*]),_EXTRA=_extra_MCMC)) else chi2_prior_new=dblarr(Nparallel)
           
    ratio = exp((-(chi2_new-chi2_old)-(chi2_prior_new-chi2_prior_old))/2.d0)
    
    accept=where(ratio gt randomu(seed,Nparallel),/null,compl=reject)
    if n_elements(accept) ne 0 then begin
      chain[*,trial+1,accept] = chain_new[*,accept]
      chain_old[*,accept] = chain_new[*,accept]            
      chi2_chain[trial+1,accept] = chi2_new[accept]
      chi2_old[accept] = chi2_new[accept]
      chi2_prior_chain[trial+1,accept] = chi2_prior_new[accept]
      chi2_prior_old[accept] = chi2_prior_new[accept]
      accepted_trials[accept]++
    endif
    if n_elements(reject) ne 0 then begin
      chain[*,trial+1,reject] = chain_old[*,reject]
      chain_old[*,reject] = chain_old[*,reject]
      chi2_chain[trial+1,reject] = chi2_old[reject]
      chi2_old[reject] = chi2_old[reject]
      chi2_prior_chain[trial+1,reject] = chi2_prior_old[reject]
      chi2_prior_old[reject] = chi2_prior_old[reject]
    endif

    if keyword_set(adaptive) then begin
      lambda_new = 10.d0^(alog10(lambda_old)+gamma_new*(min([[replicate(1,Nparallel)],[ratio]],dim=2)-alpha_star))
      mu_new=mu_old+gamma_new*(chain_old-mu_old)
      sigma_temp = dblarr(Npar,Npar,nparallel)
      for nsig=0,(Nparallel-1) do sigma_temp[*,*,nsig]=(chain_old[variable_param,nsig]-mu_old[variable_param,nsig])#(chain_old[variable_param,nsig]-mu_old[variable_param,nsig])
      sigma_new=sigma_old+gamma_new*(sigma_temp-sigma_old)
    endif
  endfor

  steps_Mstar = steps.Mstar

  ;if n_elements(Tuffs_attenuation) eq 0 then chain=chain[0:(nsteps+2),*]
  
  save,chain,chi2_chain,galaxy_id,sigma_new,lambda_new,mu_new,steps_Mstar,steps_bounds, $
       accepted_trials,chi2_prior_chain,filename=outfolder+galaxy_id+'_chain.idl'
     
  if keyword_set(print_time) then begin
    t1 = systime(/sec)
    dt = t1-t0
    nT=string(ceil(alog10(dt)+2) > 6)
    nC=string(ceil(alog10(accepted_trials[0]+1)))
    
    print,''
      print,'Done in '+string(dt,form='(F'+nT+'.1)')+' seconds'
      print,'Number of accepted trials: '+string(accepted_trials[0],form='(I'+nC+')')
      print,'Percentage of total trials: '+string(accepted_trials[0]/ntrials*100.d0,form='(F5.1)')+' %'
      print,''
  endif

end
