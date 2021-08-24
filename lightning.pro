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


function steps_derivative_chisqr_minimize,Lnu_obs,Lnu_unc,Lnu_steps,coeff1=coeff1,first_status=status1,status=status

; Function returns non-negative coefficients, possibly negative coefficients can be returned with coeff1
; status of the first inversion is status1 ; 0 is good, 1 & 2 means problems with the matrix inversion
; status of the second inversion is status2

; Lnu_obs[Nfilters]
; Lnu_unc[Nfilters]
; Lnu_steps[Nfilters,Nsteps]

;Nsteps=(size(Lnu_steps))[2]
Nsteps=(Lnu_steps.dim)[1]

coeff1 = dblarr(Nsteps,/nozero) ; first inversion coefficients, possibly negative
coeff2 = dblarr(Nsteps,/nozero) ; posterior inversions, after all coefficients are non-negative

A = dblarr(Nsteps,Nsteps,/nozero)
for j=0,(Nsteps-1) do for k=0,(Nsteps-1) do A[k,j] = total(Lnu_steps[*,j] * Lnu_steps[*,k] / Lnu_unc^2,/NaN)
;this array only need to be calculated once...

B = dblarr(Nsteps,/nozero)
for k=0,(Nsteps-1) do  B[k] = total(Lnu_steps[*,k] * Lnu_obs / Lnu_unc^2,/NaN)

coeff1 = invert(A,status1,/double) # B
;print,coeff1_sim[*,ss],form='(F10.2)'
coeff2 = coeff1
while total(coeff2 lt 0) do begin
  good = where(coeff2 gt 0,complement=bad,/null)
  Bgood = B[good]
  coeff2[bad] = 0
  if n_elements(good) then begin
    Agood = (A[good,*])[*,good]
    coeff2[good] = invert(Agood,status,/double) # Bgood
  endif
endwhile

return,coeff2

end


function calz_plus_bump_unred,Lnu,wave,tauV_DIFF=tauV_DIFF,delta=delta,tauV_BC=tauV_BC,RV=RV,exp_tau=exp_tau

;wave[Nwave] in microns
;Lwave[Nwave] or Lnu[Nwave] in any units

if (n_elements(delta) eq 0) then delta = 0d0 ; Calzetti curve
if (n_elements(tauV_BC) eq 0) then tauV_BC = 0d0
if (n_elements(RV) eq 0) then RV = 4.05 ; RV for Calzetti curve

klam = wave*0.d0
flam1 = wave*0.d0
if (tauV_DIFF ne 0) then begin
  FWHM_bump = 0.0350d0 ; 350 Angs
  lambda0_bump = 0.2175d0   ; 2175 Angs
  ;Ebump = 0.85 - 1.9*delta ;
  ;Dlam = Ebump / ( ((wave^2 - lambda0_bump^2)/(wave*FWHM_bump))^2 + 1.d0 )
  Dlam = (0.85d0-1.9d0*delta)/(((wave^2 - lambda0_bump^2)/(wave*FWHM_bump))^2 + 1.d0)

  ;RV = 4.05 ;Calzetti et al. 2000
  ;AV = 1.0
  ;EBV = AV / RV
  ;
  ; Calzettic curve extrapolated from 0.01 to 5.0 microns
  ; originally derived from 0.0912 to 2.2 microns
  w1 = where((wave GE 0.6300) AND (wave LE 2.2000),n1,/null)
  w2 = where((wave GE 0.0912) AND (wave LT 0.6300),n2,/null)
  w3 = where((wave LT 0.0912) AND (wave GT 2.2000),n3,/null)
  x  = 1.d0/wave                      ;Wavelength in inverse microns
  if (n1 ne 0) then klam[w1] = 2.659d0*(-1.857d0 + 1.040d0*x[w1]) + RV
  if (n2 ne 0) then klam[w2] = 2.659d0*(poly(x[w2], [-2.156, 1.509d0, -0.198d0, 0.011d0])) + RV
  if (n3 ne 0) then klam[w3] = 0.

  flam1 = (klam + Dlam)/4.05 * (wave/0.55)^delta ;flam1 = tau_DIFF(lambda) / tauV_DIFF
endif

flam2 = wave*0.d0              ;flam2 = tau_BC(lambda) / tauV_BC
if (tauV_BC ne 0) then begin
  flam2  = 1.d0 / (wave/0.55d0) ; birth cloud attenuation
endif

exp_tau = exp(tauV_DIFF*flam1 + tauV_BC*flam2)

return,(Lnu*exp_tau)

end


function matrix_unred,Lnu,wave,tauV_DIFF_arr=tauV_DIFF_arr,delta_arr=delta_arr,tauV_BC_arr=tauV_BC_arr

;wave[Nwave] in microns
;Lwave[Nwave] or Lnu[Nwave] in any units

Nwave = wave.length
n1 = tauV_DIFF_arr.length
n2 = delta_arr.length
n3 = tauV_BC_arr.length

Dlam = dblarr(Nwave,n2)
  FWHM_bump = 0.0350d0 ; 350 Angs
  lambda0_bump = 0.2175d0   ; 2175 Angs
  ;Ebump = 0.85 - 1.9*delta ;
  ;Dlam = Ebump / ( ((wave^2 - lambda0_bump^2)/(wave*FWHM_bump))^2 + 1.d0 )
  Dlam = (replicate(1.d0,Nwave) # (0.85d0-1.9d0*delta_arr)) / $
         ((((wave^2 - lambda0_bump^2)/(wave*FWHM_bump))^2 + 1.d0) # replicate(1.d0,n2))

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

flam1 = dblarr(Nwave,n2)
flam1 = ((klam # replicate(1.d0,n2)) + Dlam)/4.05 * $
        (wave # replicate(1.d0,n2)/0.55)^(replicate(1.d0,Nwave) # delta_arr)
        ;flam1(lambda,delta) = tau_DIFF(lambda,delta) / tauV_DIFF

; birth cloud attenuation
;flam2 = tau_BC(lambda) / tauV_BC
flam2 = 1.d0 / (wave/0.55d0)

tau_lam_DIFF = dblarr(Nwave,n1,n2)
for k1=0,(n1-1) do tau_lam_DIFF[*,k1,*] = (tauV_DIFF_arr[k1])[0]*flam1
flam=!null

tau_lam_DIFF_4D = dblarr(Nwave,n1,n2,n3)
for k3=0,(n3-1) do tau_lam_DIFF_4D[*,*,*,k3] = tau_lam_DIFF & tau_lam_DIFF = !null

tau_lam_BC = flam2 # tauV_BC_arr
tau_lam_BC_3D = dblarr(nwave,n1,n3)
for k1=0,(n1-1) do tau_lam_BC_3D[*,k1,*] = tau_lam_BC       & tau_lam_BC = !null
tau_lam_BC_4D = dblarr(nwave,n1,n2,n3)
for k2=0,(n2-1) do tau_lam_BC_4D[*,*,k2,*] = tau_lam_BC_3D  & tau_lam_BC_3D = !null

exp_tau_4D = exp(temporary(tau_lam_DIFF_4D) + temporary(tau_lam_BC_4D))

Lnu_4D = dblarr(nwave,n1,n2,n3)
for k=0,(nwave-1) do Lnu_4D[k,*,*,*] = (Lnu[k])[0]

;exp_tau[k,k1,k2,k3] = exp(tauV_DIFF[k1]*flam1[k,k2] + tauV_BC[k3]*flam2[k])
;Lnu[k]*exp_tau[k,k1,k2,k3]
;Lnu_unred[k,k1,k2,k3] = Lnu[k]*exp(tauV_DIFF[k1]*flam1[k,k2] + tauV_BC[k3]*flam2[k])

return,(Lnu_4D*exp_tau_4D)

end


function matrix_unred_Keith,wave,tauV_DIFF_arr=tauV_DIFF_arr,delta_arr=delta_arr,tauV_BC_arr=tauV_BC_arr,$
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

Dlam = dblarr(Nwave,n2)
  FWHM_bump = 0.0350d0 ; 350 Angs
  lambda0_bump = 0.2175d0   ; 2175 Angs
  ;Ebump = 0.85 - 1.9*delta ;
  ;Dlam = Ebump / ( ((wave^2 - lambda0_bump^2)/(wave*FWHM_bump))^2 + 1.d0 )
  Dlam = (0.85d0-1.9d0*rebin(reform(delta_arr,1,n2),nwave,n2)) / $
         (rebin((((wave^2 - lambda0_bump^2)/(wave*FWHM_bump))^2 + 1.d0),nwave,n2))

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

flam1 = dblarr(Nwave,n2)
if ~keyword_set(Calzetti_attenuation) then begin
  flam1 = (rebin(klam,nwave,n2) + Dlam)/4.05 * $
          (rebin(wave,nwave,n2)/0.55)^(rebin(reform(delta_arr,1,n2),nwave,n2))
          ;flam1(lambda,delta) = tau_DIFF(lambda,delta) / tauV_DIFF
endif else begin
  flam1 = rebin(klam,nwave,n2)/4.05
endelse

; birth cloud attenuation
;flam2 = tau_BC(lambda) / tauV_BC
flam2 = 1.d0 / (wave/0.55d0)

tau_lam_diff_4D = rebin(reform(flam1,nwave,1,n2,1),nwave,n1,n2,n3)*rebin(reform(tauv_diff_arr,1,n1,1,1),nwave,n1,n2,n3)

tau_lam_BC = flam2 # tauV_BC_arr
tau_lam_BC_4D=rebin(reform(tau_lam_BC,nwave,1,1,n3),nwave,n1,n2,n3)

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


function calz_plus_bump_unred_Keith,wave,tauV_DIFF=tauV_DIFF,delta=delta,tauV_BC=tauV_BC,$
                                    Calzetti_attenuation=Calzetti_attenuation

;wave[Nwave] in microns


if (n_elements(tauV_DIFF) eq 0) then tauV_DIFF = 0.d0
if (n_elements(delta) eq 0) then delta = 0d0 ; Calzetti curve
if (n_elements(tauV_BC) eq 0) then tauV_BC = 0d0
if (n_elements(RV) eq 0) then RV = 4.05 ; RV for Calzetti curve

if tauv_diff gt 0 then message,'Error -- TauB_f must be negative'


klam = wave*0.d0
flam1 = wave*0.d0
if (tauV_DIFF ne 0) then begin
  FWHM_bump = 0.0350d0 ; 350 Angs
  lambda0_bump = 0.2175d0   ; 2175 Angs
  ;Ebump = 0.85 - 1.9*delta ;
  ;Dlam = Ebump / ( ((wave^2 - lambda0_bump^2)/(wave*FWHM_bump))^2 + 1.d0 )
  Dlam = (0.85d0-1.9d0*delta)/(((wave^2 - lambda0_bump^2)/(wave*FWHM_bump))^2 + 1.d0)

  ;RV = 4.05 ;Calzetti et al. 2000
  ;AV = 1.0
  ;EBV = AV / RV
  ;
  ; Calzettic curve extrapolated from 0.01 to 5.0 microns
  ; originally derived from 0.0912 to 2.2 microns
  w1 = where((wave GE 0.6300) AND (wave LE 2.2000),n1,/null)
  w2 = where((wave GE 0.0912) AND (wave LT 0.6300),n2,/null)
  w3 = where((wave LT 0.0912) AND (wave GT 2.2000),n3,/null)
  x  = 1.d0/wave                      ;Wavelength in inverse microns
  if (n1 ne 0) then klam[w1] = 2.659d0*(-1.857d0 + 1.040d0*x[w1]) + RV
  if (n2 ne 0) then klam[w2] = 2.659d0*(poly(x[w2], [-2.156, 1.509d0, -0.198d0, 0.011d0])) + RV
  if (n3 ne 0) then klam[w3] = 0.

  if ~keyword_set(Calzetti_attenuation) then begin
    flam1 = (klam + Dlam)/4.05 * (wave/0.55)^delta ;flam1 = tau_DIFF(lambda) / tauV_DIFF
  endif else begin
    flam1 = klam/4.05
  endelse

endif

flam2 = wave*0.d0              ;flam2 = tau_BC(lambda) / tauV_BC
if (tauV_BC ne 0) then begin
  flam2  = 1.d0*0.55d0 /(wave) ; birth cloud attenuation
endif

if ~keyword_set(Calzetti_attenuation) then begin
  exp_tau = exp(tauV_DIFF*flam1 + tauV_BC*flam2)
endif else begin
  exp_tau = exp(tauV_DIFF*flam1)
endelse


return,exp_tau

end


function tuffs_unred,wave,tauB_f=tauB_f,rdisk=rdisk,rbulge=rbulge,F=F,cosi=cosi,lightning_folder=lightning_folder

;input wavelengths must be in microns
if (n_elements(tauB_f) eq 0) then tauB_f = 0.d0
if (n_elements(rdisk) eq 0)  then rdisk  = 0.d0
if (n_elements(rbulge) eq 0) then rbulge = 0.d0
if (n_elements(F) eq 0) 	 then F      = 0.d0
if (n_elements(cosi) eq 0)   then cosi   = 1.d0
if (n_elements(lightning_folder) eq 0) then lightning_folder='~/Lightning/'


if cosi gt 1 or cosi lt 0 then message,'Error -- cosi must be between 0 and 1'
if F ge 0.61 or F lt 0 then message,'Error -- F must be between 0 and 0.61'
if rbulge gt 1 or rbulge lt 0 then message,'Error -- rbulge must be between 0 and 1'
if rdisk gt 1 or rdisk lt 0 then message,'Error -- rdisk must be between 0 and 1'
if tauB_f lt 0 then message,'Error -- TauB_f must be greater than 0'
;if max() gt 1 then message,'Error -- (rdisk + rbulge) must be less than or equal to 1'


;coeff.dat file from ftp site (http://cdsarc.u-strasbg.fr/ftp/J/A+A/527/A109/coeff.dat) in Popescu et al (2011)
;converted to convenient 3d matrix where dimensions are [tau,wavebands,coeff_number] for each component
restore,lightning_folder+'tuffs_coeff.idl'

tau=[0.1,0.3,0.5,1.0,2.0,4.0,8.0]  ;tau values used in Table 3 in Popescu et al (2011)

;Include waveband at 50000 Angstroms where mlambda is 0 for interpolation smoothness
wavebands=[912.,1350,1500,1650,2000,2200,2500,2800,3650,4430,5640,8090,12590,22000,50000]   ;wavelengths of the bands used (Angstroms)

;Values from Table E.4 in Popescu et al (2011)
Fcal=0.35d0
wave_flambda=[912.,1350,1500,1650,2000,2200,2500,2800,3650,4430,5640,8090,12590,22000,50000] ;wavelengths are in Angstroms
flamb_fcal_column=[0.427,0.484,0.527,0.545,0.628,0.676,0.739,0.794,0.892,0.932,0.962,0.984,0.991,0.994,0.999]

;[value=(1−Fcal*fλ)]->[(1-value)/Fcal=fλ]
flambda=(1-flamb_fcal_column)/Fcal ;fλ values

;Add a tau at 0 for interpolation smoothness
;Set the diffuse attenuation to 0 for this case
tau=[0.0,tau]
ntau=n_elements(tau)
nwave=n_elements(wave)
nwavebands=n_elements(wavebands)

;calculate delta_mlambda for each component of galaxy  
;Using Equation 6 from Tuffs et al (2004)
;First tau index is tau=0 and last wavelength index is wavelength=50000 Angstroms. Set both to 0 for mlambda
One_minus_cosi=1.0-cosi
mlambda_d=dblarr(ntau,nwavebands)
mlambda_td=dblarr(ntau,nwavebands)
mlambda_b=dblarr(ntau,nwavebands)

mlambda_d[1:*,0:-2]=(One_minus_cosi)[0]^0.d0*adisk[*,*,0]+(One_minus_cosi)[0]^1.d0*adisk[*,*,1]+(One_minus_cosi)[0]^2.d0*adisk[*,*,2]+$
			(One_minus_cosi)[0]^3.d0*adisk[*,*,3]+(One_minus_cosi)[0]^4.d0*adisk[*,*,4]+(One_minus_cosi)[0]^5.d0*adisk[*,*,5]
mlambda_td[1:*,0:-2]=(One_minus_cosi)[0]^0.d0*atdisk[*,*,0]+(One_minus_cosi)[0]^1.d0*atdisk[*,*,1]+(One_minus_cosi)[0]^2.d0*atdisk[*,*,2]+$
			(One_minus_cosi)[0]^3.d0*atdisk[*,*,3]+(One_minus_cosi)[0]^4.d0*atdisk[*,*,4]+(One_minus_cosi)[0]^5.d0*atdisk[*,*,5]
mlambda_b[1:*,0:-2]=(One_minus_cosi)[0]^0.d0*abulge[*,*,0]+(One_minus_cosi)[0]^1.d0*abulge[*,*,1]+(One_minus_cosi)[0]^2.d0*abulge[*,*,2]+$
			(One_minus_cosi)[0]^3.d0*abulge[*,*,3]+(One_minus_cosi)[0]^4.d0*abulge[*,*,4]+(One_minus_cosi)[0]^5.d0*abulge[*,*,5]
			

;calculate delta_mlambda for entire galaxy using Equation 16 from Tuffs et al (2004)
flambda_mat=rebin(reform(flambda,1,nwavebands),ntau,nwavebands)
att_tdisk_temp=((1.d0-(rdisk)[0]-(rbulge)[0]))/((1-((F)[0]*temporary(flambda_mat))))
;NaNs occur due to 0/0. If this occurs there is no tdisk contribution to emission therefore there can be no attenuation
att_tdisk_temp[where(finite(att_tdisk_temp,/nan) eq 1,/null)]=0.d0
mlambda=2.5d0*alog10(((rdisk)[0]*10.d0^(temporary(mlambda_d)/2.5d0))+(temporary(att_tdisk_temp)*$
		10.d0^(temporary(mlambda_td)/2.5d0))+((rbulge)[0]*10.d0^(temporary(mlambda_b)/2.5d0)))

;interpolate data for input tauB_f's and input wavelength's
delta_m_tau=dblarr(nwavebands)
for i=0,(nwavebands-1) do begin $
  delta_m_tau[i]=interpol(mlambda[*,i],tau,tauB_f)
endfor
delta_m=interpol(delta_m_tau,wavebands,(wave*10000.d0))
delta_m[where(wave gt 5.0d0,/null)]=0.0d0
delta_m[where(wave lt 9.12d-2,/null)]=0.0d0

;F = F_0 * exp(-tau) = F_0 * 10^(-0.4*delta_m)
;e^(-tau) = 10^(-0.4*delta_m)

exp_tau=10.d0^(-0.4d0*delta_m)

;output form is [wavelength,tauB_f,rdisk,F,cosi,rbulge]
return,exp_tau

end


function tuffs_matrix_unred,wave,tauB_f_arr=tauB_f_arr,rdisk0_arr=rdisk0_arr,rbulge0_arr=rbulge0_arr, $
			F_arr=F_arr,cosi_arr=cosi_arr,lightning_folder=lightning_folder

;input wavelengths must be in microns
if (n_elements(tauB_f_arr) eq 0) then tauB_f_arr = 0.d0
if (n_elements(rdisk0_arr) eq 0)  then rdisk0_arr  = 0.d0
if (n_elements(rbulge0_arr) eq 0) then rbulge0_arr = 0.d0
if (n_elements(F_arr) eq 0) 	 then F_arr      = 0.d0
if (n_elements(cosi_arr) eq 0)   then cosi_arr   = 1.d0
if (n_elements(lightning_folder) eq 0) then lightning_folder='~/Lightning/'


if max(cosi_arr) gt 1 or min(cosi_arr) lt 0 then message,'Error -- cosi must be between 0 and 1'
if max(F_arr) ge 0.61 or min(F_arr) lt 0 then message,'Error -- F must be between 0 and 0.61'
if max(rbulge0_arr) gt 1 or min(rbulge0_arr) lt 0 then message,'Error -- rbulge0 must be between 0 and 1'
if max(rdisk0_arr) gt 1 or min(rdisk0_arr) lt 0 then message,'Error -- rdisk0 must be between 0 and 1'
if min(tauB_f_arr) lt 0 then message,'Error -- TauB_f must be greater than 0'


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


mlambda_d=dblarr(ntau,nwavebands,ncosi)
mlambda_td=dblarr(ntau,nwavebands,ncosi)
mlambda_b=dblarr(ntau,nwavebands,ncosi)

;Calculate delta_mlambda for each component of galaxy  
;Using Equation 6 from Tuffs et al (2004)
;First tau index is tau=0 and last wavelength index is wavelength=50000 Angstroms. Set both to 0 for mlambda
One_minus_cosi_arr=rebin(reform((1-cosi_arr),1,1,ncosi),ntau,nwavebands,ncosi)
for i=0,(ncosi-1) do begin
	mlambda_d[1:*,0:-2,i]= (One_minus_cosi_arr[1:*,0:-2,i])^0.d0*adisk[*,*,0]+(One_minus_cosi_arr[1:*,0:-2,i])^1.d0*adisk[*,*,1]+$
                           (One_minus_cosi_arr[1:*,0:-2,i])^2.d0*adisk[*,*,2]+(One_minus_cosi_arr[1:*,0:-2,i])^3.d0*adisk[*,*,3]+$
                           (One_minus_cosi_arr[1:*,0:-2,i])^4.d0*adisk[*,*,4]+(One_minus_cosi_arr[1:*,0:-2,i])^5.d0*adisk[*,*,5]
	mlambda_td[1:*,0:-2,i]=(One_minus_cosi_arr[1:*,0:-2,i])^0.d0*atdisk[*,*,0]+(One_minus_cosi_arr[1:*,0:-2,i])^1.d0*atdisk[*,*,1]+$
                           (One_minus_cosi_arr[1:*,0:-2,i])^2.d0*atdisk[*,*,2]+(One_minus_cosi_arr[1:*,0:-2,i])^3.d0*atdisk[*,*,3]+$
                           (One_minus_cosi_arr[1:*,0:-2,i])^4.d0*atdisk[*,*,4]+(One_minus_cosi_arr[1:*,0:-2,i])^5.d0*atdisk[*,*,5]
	mlambda_b[1:*,0:-2,i]= (One_minus_cosi_arr[1:*,0:-2,i])^0.d0*abulge[*,*,0]+(One_minus_cosi_arr[1:*,0:-2,i])^1.d0*abulge[*,*,1]+$
                           (One_minus_cosi_arr[1:*,0:-2,i])^2.d0*abulge[*,*,2]+(One_minus_cosi_arr[1:*,0:-2,i])^3.d0*abulge[*,*,3]+$
                           (One_minus_cosi_arr[1:*,0:-2,i])^4.d0*abulge[*,*,4]+(One_minus_cosi_arr[1:*,0:-2,i])^5.d0*abulge[*,*,5]
endfor

;calculate delta_mlambda for entire galaxy using Equation 16 from Tuffs et al (2004)
rdisk0_mat=rebin(reform(rdisk0_arr,1,1,nrdisk0,1,1,1),ntau,nwavebands,nrdisk0,nF,ncosi,nrbulge0)
rbulge0_mat=rebin(reform(rbulge0_arr,1,1,1,1,1,nrbulge0),ntau,nwavebands,nrdisk0,nF,ncosi,nrbulge0)
att_tdisk_1=(1.d0-temporary(rdisk0_mat)-temporary(rbulge0_mat))

F_mat=rebin(reform(F_arr,1,1,1,nF,1,1),ntau,nwavebands,nrdisk0,nF,ncosi,nrbulge0)
flam_mat=rebin(reform(flambda,1,nwavebands,1,1,1,1),ntau,nwavebands,nrdisk0,nF,ncosi,nrbulge0)
att_tdisk_temp=temporary(att_tdisk_1)*(1-(temporary(F_mat)*temporary(flam_mat)))

mlambda_td_mat=rebin(reform(mlambda_td,ntau,nwavebands,1,1,ncosi,1),ntau,nwavebands,nrdisk0,nF,ncosi,nrbulge0)
att_tdisk=temporary(att_tdisk_temp)*10.d0^(-1.d0*temporary(mlambda_td_mat)/2.5d0)


rdisk0_mat=rebin(reform(rdisk0_arr,1,1,nrdisk0,1,1,1),ntau,nwavebands,nrdisk0,nF,ncosi,nrbulge0)
mlambda_d_mat=rebin(reform(mlambda_d,ntau,nwavebands,1,1,ncosi,1),ntau,nwavebands,nrdisk0,nF,ncosi,nrbulge0)
att_disk=temporary(rdisk0_mat)*10.d0^(-1.d0*temporary(mlambda_d_mat)/2.5d0)
;There is no attenuation from the disk in the UV
att_disk[*,0:8,*,*,*,*]=0.d0

att_tdisk_disk=temporary(att_tdisk)+temporary(att_disk)


rbulge0_mat=rebin(reform(rbulge0_arr,1,1,1,1,1,nrbulge0),ntau,nwavebands,nrdisk0,nF,ncosi,nrbulge0)
mlambda_b_mat=rebin(reform(mlambda_b,ntau,nwavebands,1,1,ncosi,1),ntau,nwavebands,nrdisk0,nF,ncosi,nrbulge0)
att_bulge=temporary(rbulge0_mat)*10.d0^(-1.d0*temporary(mlambda_b_mat)/2.5d0)
;There is no attenuation from the bulge in the UV
att_bulge[*,0:8,*,*,*,*]=0.d0


mlambda=-2.5d0*alog10(temporary(att_tdisk_disk)+temporary(att_bulge))
;There is no attenuation at 5um or taub_f=0
mlambda[0,*,*,*,*,*]=0.d0
mlambda[*,-1,*,*,*,*]=0.d0
;Infinities occur if rdisk or rbulge=1, due to no attenuation in UV. Set to 0
mlambda[where(finite(mlambda) eq 0,/null)]=0.d0



;interpolate data for input tauB_f's and input wavelength's
delta_m_tau=dblarr(nwavebands,ntauB_f,nrdisk0,nF,ncosi,nrbulge0)
delta_m=dblarr(nwave,ntauB_f,nrdisk0,nF,ncosi,nrbulge0)
for i=0,(nrbulge0-1) do begin
  for j=0,(ncosi-1) do begin 
    for k=0,(nF-1) do begin 
      for l=0,(nrdisk0-1) do begin 
        for m=0,(nwavebands-1) do begin 
          delta_m_tau[m,*,l,k,j,i]=interpol(mlambda[*,m,l,k,j,i],tau,tauB_f_arr)
        endfor
        for n=0,(ntauB_F-1) do begin 
          delta_m[*,n,l,k,j,i]=interpol(delta_m_tau[*,n,l,k,j,i],wavebands,(wave*10000.d0))
          ;Set attenuation to 0 at wavelengths greater than 5.0 microns. Do this in case interpolation has error.
          delta_m[where(wave gt 5.d0 or wave lt 9.12d-2,/null),n,l,k,j,i]=0.0
        endfor 
      endfor 
    endfor 
  endfor
endfor

DELVAR,mlambda
DELVAR,delta_m_tau

;F = F_0 * exp(-tau) = F_0 * 10^(-0.4*delta_m)
;e^(-tau) = 10^(-0.4*delta_m)

exp_tau=10.d0^(-0.4d0*delta_m)

;output form is [wavelength,tauB_f,rdisk,F,cosi,rbulge]
return,exp_tau

end


function steps_fitting,Lnu_obs,Lnu_unc,Nsim=Nsim,steps=steps,grid=grid
t0 = SYSTIME(/seconds)

Nfilters = (Lnu_obs.DIM)[0]
if (SIZE(Lnu_obs))[0] eq 2 then Nreg = (Lnu_obs.DIM)[1]
if (N_ELEMENTS(Nsim) eq 0) then Nsim=1

if (((Lnu_unc.Dim)[0] ne Nfilters) or ((Lnu_unc.DIM)[1] ne Nreg)) then begin
  message,'ERROR - SED and uncertainties must have same sizes'
  return,!null
endif

Lnu_sim = dblarr(Nfilters,Nreg,Nsim)
for ks=1,(Nsim-1) do Lnu_sim[*,*,ks] = Lnu_obs[*,*,ks] + Lnu_unc[*,*] * randomn(seed,Nfilters,Nreg)
Lnu_sim[*,*,0] = Lnu_obs ; first simulation is the observed data

if (N_ELEMENTS(steps) eq 0) then begin
  steps = STEPS_STELLAR(filter_labels=filter_labels,steps_bounds=steps_bounds,z_shift=0.0,dtime_SF=5.d5)
endif
Nsteps = steps.bounds.length - 1

wave = steps.wave
Nwave = wave.length

nu = 1.d4 * !cv.clight / wave

;normalizing filters
for k=0,(Nfilters-1) do steps.filters[k,*] = steps.filters[k,*]/trap_int(nu,steps.filters[k,*])

;Lnu_true = total((replicate(1.d0,Nfilters) # coeff_input) * steps_mean_Lnu_red,2)
;Lnu_unc = 0.05 * Lnu_true
;Lnu_obs = Lnu_true; + Lnu_unc * randomn(seed,Nfilters)
;Lnu_sim = (Lnu_true # replicate(1.d0,Nsim)) + (Lnu_unc # replicate(1.d0,Nsim))*randomn(seed,Nfilters,Nsim)
;Lnu_sim[*,0] = Lnu_obs ; first simulation is the observed data

if n_elements(grid) eq 0 then begin
  n1 = 41
  n2 = 41
  n3 = 41

  ran1 = [   0d,4d0] & tauV_DIFF_grid = ran1[0] + dindgen(n1)*(ran1[1]-ran1[0])/(n1-1)
  ran2 = [-2.2d,.4d] & delta_grid     = ran2[0] + dindgen(n2)*(ran2[1]-ran2[0])/(n2-1)
  ran3 = [   0d,4d0] & tauV_BC_grid   = ran3[0] + dindgen(n3)*(ran3[1]-ran3[0])/(n3-1)
endif else begin
  tauV_DIFF_grid = grid.tauV_DIFF_grid & n1 = n_elements(tauV_DIFF_grid)
  delta_grid     = grid.delta_grid     & n2 = n_elements(delta_grid)
  tauV_BC_grid   = grid.tauV_BC_grid   & n3 = n_elements(tauV_BC_grid)
endelse

;as in Leja et al. 2016:
;ran1 = [   0d,4d0]  & tauV_DIFF_grid = ran1[0] + dindgen(n1)*(ran1[1]-ran1[0])/(n1-1)
;ran2 = [-2.2d,.4d]  & delta_grid     = ran2[0] + dindgen(n2)*(ran2[1]-ran2[0])/(n2-1)
;ran3 = [   0d,4d0]  & tauV_BC_grid   = ran3[0] + dindgen(n3)*(ran3[1]-ran3[0])/(n3-1)

;coeff_grid = dblarr(Nsteps,n1,n2,n3,Nreg,Nsim)
coeff_grid = dblarr(Nsteps,n1,n2,n3,Nreg,Nsim,/nozero)     & replicate_inplace,coeff_grid,0.d0
;coeff_grid = fltarr(Nsteps,n1,n2,n3,Nreg,Nsim,/nozero)     & replicate_inplace,coeff_grid,1.e0
;status_grid = dblarr(n1,n2,n3)
;chisqr_grid = dblarr(n1,n2,n3,Nreg,Nsim)
chisqr_grid = dblarr(n1,n2,n3,Nreg,Nsim,/nozero)           & replicate_inplace,chisqr_grid,1.e0
;chisqr_grid = fltarr(n1,n2,n3,Nreg,Nsim,/nozero)           & replicate_inplace,chisqr_grid,1.e0
;Lnu_mod_grid = dblarr(Nfilters,n1,n2,n3,Nreg,Nsim)
;Lnu_mod_grid = dblarr(Nfilters,n1,n2,n3,Nreg,Nsim,/nozero) & replicate_inplace,Lnu_mod_grid,0.d0
Lnu_mod_grid = dblarr(Nfilters,n1,n2,n3,Nreg,Nsim,/nozero) & replicate_inplace,Lnu_mod_grid,0.d0
;Lnu_mod_grid = fltarr(Nfilters,n1,n2,n3,Nreg,Nsim,/nozero) & replicate_inplace,Lnu_mod_grid,0.e0
;steps_mean_Lnu_red_grid = dblarr(Nfilters,n1,n2,n3,Nsteps)

steps_Lnu_red = steps.Lnu * 0.d0
steps_mean_Lnu_red = dblarr(Nfilters,Nsteps,/nozero)       & replicate_inplace,steps_mean_Lnu_red,0.d0
;steps_mean_Lnu_red = fltarr(Nfilters,Nsteps,/nozero)       & replicate_inplace,steps_mean_Lnu_red,0.e0
steps_mean_Lnu_red_grid = dblarr(Nfilters,Nsteps,n1,n2,n3,/nozero) & replicate_inplace,steps_mean_Lnu_red_grid,0.d0
;steps_mean_Lnu_red_grid = fltarr(Nfilters,Nsteps,n1,n2,n3,/nozero) & replicate_inplace,steps_mean_Lnu_red_grid,0.e0
t1 = systime(/seconds)
for kk=0,(n3-1) do begin
  tauV_BC = tauV_BC_grid[kk]
  for jj=0,(n2-1) do begin
    delta = delta_grid[jj]
    for ii=0,(n1-1) do begin
      tauV_DIFF = tauV_DIFF_grid[ii]
      for st=0,(Nsteps-1) do begin
        steps_Lnu_red[*,st] = calz_plus_bump_unred(steps.Lnu[*,st],steps.wave,$
          tauV_DIFF=-1.d0*tauV_DIFF,delta=delta,tauV_BC=(tauV_BC*[-1d,0d,0d,0d,0d])[st]);,exp_tau=exp_tau)
        ;flam[*,st] = alog(exp_tau)/alog(interpol(exp_tau,wave,.55d0))
      endfor
      for k=0,(Nfilters-1) do begin
        filter_transmission = reform(steps.Filters[k,*])
        for st=0,(Nsteps-1) do begin
          steps_mean_Lnu_red[k,st] = trap_int(nu,steps_Lnu_red[*,st]*filter_transmission);/trap_int(nu,filter_transmission)
        endfor
      endfor
      steps_mean_Lnu_red_grid[*,*,ii,jj,kk] = steps_mean_Lnu_red
      ;could save models to save time on next runs... File would not be too big ~ few megabytes
    endfor
  endfor
endfor
t2=systime(/sec)
print,'Generate models: ',t2-t1,' seconds'

t3=systime(/sec)
for ss=0,(Nsim-1) do begin
  Lnu_tmp = Lnu_sim[*,*,ss]
  for sr=0,(Nreg-1) do begin
    Lnu_fit = Lnu_tmp[*,sr]
    Lnu_sig = Lnu_unc[*,sr]
    for kk=0,(n3-1) do begin
      for jj=0,(n2-1) do begin
        for ii=0,(n1-1) do begin
          steps_mean_Lnu_red = steps_mean_Lnu_red_grid[*,*,ii,jj,kk]
          coeff_temp = $
            steps_derivative_chisqr_minimize(Lnu_fit,Lnu_sig,steps_mean_Lnu_red)
          ;status_grid[ii,jj,kk] = temporary(status_temp)
          Lnu_mod_temp = total((replicate(1.d0,Nfilters) # coeff_temp)*steps_mean_Lnu_red,2)
          coeff_grid[*,ii,jj,kk,sr,ss] = temporary(coeff_temp)
          chisqr_grid[ii,jj,kk,sr,ss] = total(((Lnu_mod_temp - Lnu_fit)/Lnu_sig)^2)
          Lnu_mod_grid[*,ii,jj,kk,sr,ss] = temporary(Lnu_mod_temp)
        endfor
      endfor
    endfor
    ti=systime(/seconds)
    print,'Finished outer loop #'+string(sr+1,form='(I3)')+' out of '+strtrim(Nreg,2)+'. '+$
          string((sr+1.)*100./Nreg,form='(F5.1)')+'% completed. '+string(ti-t3,form='(F10.2)')+$
          ' seconds elapsed. Min chisqr: '+string(min(chisqr_grid[*,*,*,sr,ss]),form='(F10.2)')
  endfor
endfor
t4 = systime(/seconds)
print,'Main loop running time: ',t4-t3,' seconds'
print,'Total running time: ',t4-t0,' seconds'

output = {chisqr_grid:temporary(chisqr_grid),coeff_grid:temporary(coeff_grid),Lnu_mod_grid:temporary(Lnu_mod_grid)}

return,output

end


pro plot_lightning,Lnu_obs,Lnu_unc,LIR_obs,LIR_unc,steps=steps,grid=grid,fit=fit,bestfit=bestfit,$
                   maxpost=maxpost,outfolder=outfolder,Regions=regions,noprior=noprior,$
                   Tuffs_attenuation=Tuffs_attenuation
;Modification History
; Keith Doore     - added Tuffs extinction option - 2019
; Keith Doore     - Replaced if statements with n_elements() on keywords to use keyword_set() - April 27th 2020
 
if keyword_set(Tuffs_attenuation) then begin
  tauB_f_grid = grid.tauB_f_grid & n1 = N_ELEMENTS(tauB_f_grid)
  rdisk_grid  = grid.rdisk_grid  & n2 = N_ELEMENTS(rdisk_grid)
  F_grid      = grid.F_grid      & n3 = N_ELEMENTS(F_grid)
  cosi_grid   = grid.cosi_grid   & n4 = N_ELEMENTS(cosi_grid)
  rbulge_grid = grid.rbulge_grid & n5 = N_ELEMENTS(rbulge_grid)
endif else begin
  tauV_DIFF_grid = grid.tauV_DIFF_grid & n1 = n_elements(tauV_DIFF_grid)
  delta_grid     = grid.delta_grid     & n2 = n_elements(delta_grid)
  tauV_BC_grid   = grid.tauV_BC_grid   & n3 = n_elements(tauV_BC_grid)
  n4 = 1
  n5 = 1
endelse

Nsteps   = N_ELEMENTS(steps.Mstar)
Nfilters = (Lnu_obs.dim)[0]
if (SIZE(Lnu_obs))[0] eq 2 then Nreg=(Lnu_obs.DIM)[1] else Nreg=1
if (N_ELEMENTS(regions) eq 0) then begin
  Nfit = Nreg
  regions = [0,Nreg-1]
endif else begin
  Nfit = Regions[1] - Regions[0] + 1
endelse
N0 = regions[0]
Nf = regions[1]

if N_ELEMENTS(fit.CHISQR_GRID.dim) lt 7 then Nsim=1 else Nsim = (fit.CHISQR_GRID.dim)[6]

nu_obs = 1.d4*!cv.clight/steps.wave_filters

if ~keyword_set(noprior) then begin
  if keyword_set(Tuffs_attenuation) then begin
    P = DOUBLE((tauB_f_grid # REPLICATE(1.d0,n2)) le (2.0*(REPLICATE(1.d0,n1) # rdisk_grid)))
  endif else begin
    P = double((((tauV_DIFF_grid # (1.d0/tauV_BC_grid)) lt 2.0) and $
                 (tauV_DIFF_grid # (1.d0/tauV_BC_grid)) gt 0.5))
  endelse
  P3D = dblarr(n1,n2,n3,/nozero)
  for k2=0,(n2-1) do P3D[*,k2,*] = P
  P4D = dblarr(n1,n2,n3,n4,/nozero)
    for k4=0,(n4-1) do P4D[*,*,*,k4] = P3D
  P5D = dblarr(n1,n2,n3,n4,n5,/nozero)
    for k5=0,(n5-1) do P5D[*,*,*,*,k5] = P4D
    P5D = P5D/total(P5D)
  P6D = dblarr(n1,n2,n3,n4,n5,1,nsim,/nozero)
    for ss=0,(nsim-1) do P6D[*,*,*,*,*,0,ss] = P5D
endif else begin
  P5D = dblarr(n1,n2,n3,n4,n5)+1
  P5D = P5D/total(P5D)
  P6D = rebin(reform(P5D,n1,n2,n3,n4,n5,1,1),n1,n2,n3,n4,n5,1,nsim)
endelse
    
Post = dblarr(n1,n2,n3,n4,n5,/nozero)

minchisqr           = dblarr(Nfit,Nsim)
minchisqr_loc       = lonarr(Nfit,Nsim)
minchisqr_ind       = lonarr(5,Nfit,Nsim)
coeff_minchisqr     = dblarr(Nsteps,Nfit,Nsim)
minchisqr_Lnu       = dblarr(Nfilters,Nfit,Nsim)
minchisqr_Lnu_unred = dblarr(Nfilters,Nfit,Nsim)
minchisqr_LIR       = dblarr(Nfit,Nsim)

maxpost_ind     = lonarr(5,Nfit,Nsim)
chisqr_maxpost  = dblarr(Nfit,Nsim)
coeff_maxpost   = dblarr(Nsteps,Nfit,Nsim)
Lnu_mod         = dblarr(Nfilters,Nfit,Nsim)
Lnu_mod_unred   = dblarr(Nfilters,Nfit,Nsim)
LIR_mod         = dblarr(Nfit,Nsim)

rep_filt = replicate(1.d0,Nfilters)

for ss=0,(Nsim-1) do begin
  for h=0,(Nfit-1) do begin
    minchisqr[h,ss] = min(fit.chisqr_grid[*,*,*,*,*,h,ss],minchi2_loc)
    minchisqr_loc[h,ss] = minchi2_loc
  endfor
endfor

reg_form = '(I0'+strtrim(ceil(alog10(Nreg)),2)+')'
mydevice = !d.name
set_plot,'ps'
for h=0,(Nfit-1) do begin
    coeff_reg = fit.coeff_grid[*,*,*,*,*,*,h,*]
    Lmod_reg = fit.Lmod_grid[*,*,*,*,*,*,h,*]
    chisqr_reg = fit.chisqr_grid[*,*,*,*,*,h,*]
    minchisqr_reg = 0*chisqr_reg
    for ss=0,(Nsim-1) do minchisqr_reg[*,*,*,*,*,0,ss] = minchisqr[h,ss]
    L6D = exp(-(chisqr_reg-minchisqr_reg)/2.d0)
    Post6D = P6D * L6D
  for ss=0,(Nsim-1) do begin
    L = L6D[*,*,*,*,*,0,ss]       & L    /= (total(L))[0]
    Post = Post6D[*,*,*,*,*,0,ss] & Post /= (total(Post))[0]

    minchi2_ind = lonarr(5,/nozero)
    minchi2_loc = (minchisqr_loc[h,ss])[0]
    
    if (n5 eq 1) then begin
      minchi2_ind[4] = 0
      if (n4 eq 1) then begin
        minchi2_ind[3] = 0
        if (n3 eq 1) then begin
          minchi2_ind[2] = 0
          if (n2 eq 1) then begin
            minchi2_ind[1] = 0
            minchi2_ind[0] = minchi2_loc
          endif else begin
            minchi2_ind[0:1] = array_indices(dblarr(n1,n2,/nozero),minchi2_loc)
          endelse
        endif else if n2 eq 1 then begin
          minchi2_ind[1] = 0
          minchi2_ind[[0,2]] = array_indices(dblarr(n1,n3,/nozero),minchi2_loc)
        endif else begin
          minchi2_ind[0:2] = array_indices(dblarr(n1,n2,n3,/nozero),minchi2_loc)
        endelse
      endif else if n3 eq 1 then begin
        minchi2_ind[2] = 0
        if (n2 eq 1) then begin
          minchi2_ind[1] = 0
          minchi2_ind[[0,3]] = array_indices(dblarr(n1,n4,/nozero),minchi2_loc)
        endif else begin
          minchi2_ind[[0,1,3]] = array_indices(dblarr(n1,n2,n4,/nozero),minchi2_loc)
        endelse
      endif else if n2 eq 1 then begin
          minchi2_ind[1] = 0
          minchi2_ind[[0,2,3]] = array_indices(dblarr(n1,n3,n4,/nozero),minchi2_loc)
      endif else begin
          minchi2_ind[0:3] = array_indices(dblarr(n1,n2,n3,n4,/nozero),minchi2_loc)
      endelse
    endif else if n4 eq 1 then begin
      minchi2_ind[3] = 0
      if (n3 eq 1) then begin
        minchi2_ind[2] = 0
        if (n2 eq 1) then begin
          minchi2_ind[1] = 0
          minchi2_ind[[0,4]] = array_indices(dblarr(n1,n5,/nozero),minchi2_loc)
        endif else begin
          minchi2_ind[[0,1,4]] = array_indices(dblarr(n1,n2,n5,/nozero),minchi2_loc)
        endelse
      endif else if n2 eq 1 then begin
        minchi2_ind[1] = 0
        minchi2_ind[[0,2,4]] = array_indices(dblarr(n1,n3,n5,/nozero),minchi2_loc)
      endif else begin
        minchi2_ind[[0,1,2,4]] = array_indices(dblarr(n1,n2,n3,n5,/nozero),minchi2_loc)
      endelse
    endif else if n3 eq 1 then begin
      minchi2_ind[2] = 0
      if (n2 eq 1) then begin
        minchi2_ind[1] = 0
        minchi2_ind[[0,3,4]] = array_indices(dblarr(n1,n4,n5,/nozero),minchi2_loc)
      endif else begin
        minchi2_ind[[0,1,3,4]] = array_indices(dblarr(n1,n2,n4,n5,/nozero),minchi2_loc)
      endelse
    endif else if n2 eq 1 then begin
        minchi2_ind[1] = 0
        minchi2_ind[[0,2,3,4]] = array_indices(dblarr(n1,n3,n4,n5,/nozero),minchi2_loc)
    endif else begin
        minchi2_ind[0:4] = array_indices(dblarr(n1,n2,n3,n4,n5,/nozero),minchi2_loc)
    endelse
 
    
    minchisqr_ind[*,h,ss] = minchi2_ind
    coeff = coeff_reg[*,minchi2_ind[0],minchi2_ind[1],minchi2_ind[2],minchi2_ind[3],minchi2_ind[4],0,ss]
    coeff_minchisqr[*,h,ss] = coeff
    Lmod = Lmod_reg[*,minchi2_ind[0],minchi2_ind[1],minchi2_ind[2],minchi2_ind[3],minchi2_ind[4],0,ss]
    minchisqr_Lnu_unred[*,h,ss] = total((rep_filt # coeff)*steps.mean_Lnu,2)
    minchisqr_Lnu[*,h,ss] = Lmod[0:-2]
    minchisqr_LIR[h,ss] = Lmod[-1]

    maxpost = max(Post,loc)
    ind = lonarr(5,/nozero)
    if (n5 eq 1) then begin
      ind[4] = 0
      if (n4 eq 1) then begin
        ind[3] = 0
        if (n3 eq 1) then begin
          ind[2] = 0
          if (n2 eq 1) then begin
            ind[1] = 0
            ind[0] = loc
          endif else begin
            ind[0:1] = array_indices(dblarr(n1,n2,/nozero),loc)
          endelse
        endif else if n2 eq 1 then begin
          ind[1] = 0
          ind[[0,2]] = array_indices(dblarr(n1,n3,/nozero),loc)
        endif else begin
          ind[0:2] = array_indices(dblarr(n1,n2,n3,/nozero),loc)
        endelse
      endif else if n3 eq 1 then begin
        ind[2] = 0
        if (n2 eq 1) then begin
          ind[1] = 0
          ind[[0,3]] = array_indices(dblarr(n1,n4,/nozero),loc)
        endif else begin
          ind[[0,1,3]] = array_indices(dblarr(n1,n2,n4,/nozero),loc)
        endelse
      endif else if n2 eq 1 then begin
          ind[1] = 0
          ind[[0,2,3]] = array_indices(dblarr(n1,n3,n4,/nozero),loc)
      endif else begin
          ind[0:3] = array_indices(dblarr(n1,n2,n3,n4,/nozero),loc)
      endelse
    endif else if n4 eq 1 then begin
      ind[3] = 0
      if (n3 eq 1) then begin
        ind[2] = 0
        if (n2 eq 1) then begin
          ind[1] = 0
          ind[[0,4]] = array_indices(dblarr(n1,n5,/nozero),loc)
        endif else begin
          ind[[0,1,4]] = array_indices(dblarr(n1,n2,n5,/nozero),loc)
        endelse
      endif else if n2 eq 1 then begin
        ind[1] = 0
        ind[[0,2,4]] = array_indices(dblarr(n1,n3,n5,/nozero),loc)
      endif else begin
        ind[[0,1,2,4]] = array_indices(dblarr(n1,n2,n3,n5,/nozero),loc)
      endelse
    endif else if n3 eq 1 then begin
      ind[2] = 0
      if (n2 eq 1) then begin
        ind[1] = 0
        ind[[0,3,4]] = array_indices(dblarr(n1,n4,n5,/nozero),loc)
      endif else begin
        ind[[0,1,3,4]] = array_indices(dblarr(n1,n2,n4,n5,/nozero),loc)
      endelse
    endif else if n2 eq 1 then begin
        ind[1] = 0
        ind[[0,2,3,4]] = array_indices(dblarr(n1,n3,n4,n5,/nozero),loc)
    endif else begin
        ind[0:4] = array_indices(dblarr(n1,n2,n3,n4,n5,/nozero),loc)
    endelse
    
    maxpost_ind[*,h,ss] = ind
    chisqr_maxpost[h,ss] = ((chisqr_reg)[*,*,*,*,*,0,ss])[loc]
    coeff = coeff_reg[*,ind[0],ind[1],ind[2],ind[3],ind[4],0,ss]
    coeff_maxpost[*,h,ss] = coeff
    Lnu_mod_unred[*,h,ss] = total((rep_filt # coeff)*steps.mean_Lnu,2)
    Lmod = Lmod_reg[*,ind[0],ind[1],ind[2],ind[3],ind[4],0,ss]
    Lnu_mod[*,h,ss] = Lmod[0:-2]
    LIR_mod[h,ss] = Lmod[-1]
    
    
    if ss eq 0 then begin
      Nh = N0+h
      if keyword_set(Tuffs_attenuation) then begin
    print,Nh,minchisqr[h,ss],chisqr_maxpost[h,ss],coeff_maxpost[*,h,ss],tauB_f_grid[ind[0]],$
            rdisk_grid[ind[1]],F_grid[ind[2]],cosi_grid[ind[3]],rbulge_grid[ind[4]],form='(I,2F10.2,'+strtrim(Nsteps,2)+'E10.2,5F12.3)'
      endif else begin
        print,Nh,minchisqr[h,ss],chisqr_maxpost[h,ss],coeff_maxpost[*,h,ss],tauV_DIFF_grid[ind[0]],$
            delta_grid[ind[1]],tauV_BC_grid[ind[2]],form='(I,2F10.2,'+strtrim(Nsteps,2)+'E10.2,3F12.3)'
      endelse

      prefix = 'Reg'+string(Nh,form=reg_form)
      
      fileout=outfolder+prefix+'_SED.eps'
      device,filename=fileout,/encaps,bits=8,/color
        plot_oo,minmax(steps.wave_filters),$
                minmax([nu_obs*Lnu_obs[*,Nh],nu_obs*Lnu_mod_unred[*,h],$
                        nu_obs*Lnu_mod[*,h],LIR_obs[Nh],LIR_mod[h]]),$
                xtitle='!6Wavelength, !7k!6 [!7l!6m]',$
                ytitle='!7m!6L!d!7m!6!n [L'+sunsymbol()+']',title='Reg '+strtrim(Nh,2),/nodata
          oploterror,steps.wave_filters,nu_obs*Lnu_obs[*,Nh],nu_obs*Lnu_unc[*,Nh],psym=4
          oploterror,7.*[1,1],LIR_obs[Nh]*[1,1],LIR_unc[Nh]*[1,1],psym=4
          oplot,steps.wave_filters,nu_obs*Lnu_mod_unred[*,h],psym=-2,color=54,/linestyle
          oplot,steps.wave_filters,nu_obs*Lnu_mod[*,h],psym=-2,color=254,/linestyle
          oplot,7.*[1,1],LIR_mod[h]*[1,1],psym=-2,color=254,/linestyle

      device,/close

      if n1 gt 1 then begin
        if n5 gt 1 then begin
          Post1 = total(Post,5) & L1 = total(L,5) & P1 = total(P5D,5)
          if n4 gt 1 then begin
            Post1 = total(Post1,4) & L1 = total(L1,4) & P1 = total(P1,4)
            if n3 gt 1 then begin
              Post1 = total(Post1,3) & L1 = total(L1,3) & P1 = total(P1,3)
              if n2 gt 1 then begin
                Post1 = total(Post1,2) & L1 = total(L1,2) & P1 = total(P1,2)
              endif
            endif
          endif
        endif else begin
          if n4 gt 1 then begin
            Post1 = total(Post,4) & L1 = total(L,4) & P1 = total(P5D,4)
            if n3 gt 1 then begin
              Post1 = total(Post1,3) & L1 = total(L1,3) & P1 = total(P1,3)
              if n2 gt 1 then begin
                Post1 = total(Post1,2) & L1 = total(L1,2) & P1 = total(P1,2)
              endif
            endif 
          endif else begin
            if n3 gt 1 then begin
              Post1 = total(Post,3) & L1 = total(L,3) & P1 = total(P5D,3)
              if n2 gt 1 then begin
                Post1 = total(Post1,2) & L1 = total(L1,2) & P1 = total(P1,2)
              endif
            endif else begin
              if n2 gt 1 then begin
                Post1 = total(Post,2) & L1 = total(L,2) & P1 = total(P5D,2)
              endif else begin
                Post1 = Post  &  L1 = L  & P1 = P5D
              endelse
            endelse
          endelse
        endelse
        fileout=outfolder+prefix+'_par1.eps'
        device,filename=fileout,/encaps,bits=8,/color
          if keyword_set(Tuffs_attenuation) then begin
            plot,minmax(tauB_f_grid),minmax([Post1,L1,P1]),$
               xtitle='Face on attenuation, tauB_f_grid',ytitle='Probalibity Densities',/nodata
            oplot,tauB_f_grid,Post1,color=54
            oplot,[1,1]*(tauB_f_grid[ind[0]])[0],!y.crange,color=54
            oplot,tauB_f_grid,L1,linestyle=2,color=254
            oplot,[1,1]*(tauB_f_grid[minchi2_ind[0]])[0],!y.crange,linestyle=2,color=254
            oplot,tauB_f_grid,P1,linestyle=1,color=54
          endif else begin
            plot,minmax(tauV_DIFF_grid),minmax([Post1,L1,P1]),$
               xtitle='Diffuse attenuation, tau_V_DIFF',ytitle='Probalibity Densities',/nodata
            oplot,tauV_DIFF_grid,Post1,color=54
            oplot,[1,1]*(tauV_DIFF_grid[ind[0]])[0],!y.crange,color=54
            oplot,tauV_DIFF_grid,L1,linestyle=2,color=254
            oplot,[1,1]*(tauV_DIFF_grid[minchi2_ind[0]])[0],!y.crange,linestyle=2,color=254
            oplot,tauV_DIFF_grid,P1,linestyle=1,color=54
          endelse
        device,/close
      endif

      if n2 gt 1 then begin
        if n5 gt 1 then begin
          Post2 = total(Post,5) & L2 = total(L,5) & P2 = total(P5D,5)
          if n4 gt 1 then begin
            Post2 = total(Post2,4) & L2 = total(L2,4) & P2 = total(P2,4)
            if n3 gt 1 then begin
              Post2 = total(total(Post2,3),1) & L2 = total(total(L2,3),1) & P2 = total(total(P2,3),1)
            endif
          endif
        endif else begin
          if n4 gt 1 then begin
            Post2 = total(Post,4) & L2 = total(L,4) & P2 = total(P5D,4)
            if n3 gt 1 then begin
              Post2 = total(total(Post2,3),1) & L2 = total(total(L2,3),1) & P2 = total(total(P2,3),1)
            endif 
          endif else begin
            if n3 gt 1 then begin
              Post2 = total(total(Post,3),1) & L2 = total(total(L,3),1) & P2 = total(total(P5D,3),1)
            endif else begin
              Post2 = total(Post,1) & L2 = total(L,1) & P2 = total(P5D,1)
            endelse
          endelse
        endelse
        fileout=outfolder+prefix+'_par2.eps'
        device,filename=fileout,/encaps,bits=8,/color
          if keyword_set(Tuffs_attenuation) then begin
            plot,minmax(rdisk_grid),minmax([Post2,L2,P2]),$
               xtitle='!8L!6!ddisk!n/!8L!6!dtotal!n, rdisk_grid',ytitle='Probalibity Densities',/nodata
            oplot,rdisk_grid,Post2,color=54
            oplot,[1,1]*(rdisk_grid[ind[1]])[0],!y.crange,color=54
            oplot,rdisk_grid,L2,linestyle=2,color=254
            oplot,[1,1]*(rdisk_grid[minchi2_ind[1]])[0],!y.crange,linestyle=2,color=254
            oplot,rdisk_grid,P2,linestyle=1,color=54
          endif else begin
            plot,minmax(delta_grid),minmax([Post2,L2,P2]),$
               xtitle='Slope, delta',ytitle='Probalibity Densities',/nodata
            oplot,delta_grid,Post2,color=54
            oplot,[1,1]*(delta_grid[ind[1]])[0],!y.crange,color=54
            oplot,delta_grid,L2,linestyle=2,color=254
            oplot,[1,1]*(delta_grid[minchi2_ind[1]])[0],!y.crange,linestyle=2,color=254
            oplot,delta_grid,P2,linestyle=1,color=54
          endelse
        device,/close
      endif

      if n3 gt 1 then begin 
        if n5 gt 1 then begin
          Post3 = total(Post,5) & L3 = total(L,5) & P3 = total(P5D,5)
          if n4 gt 1 then begin
            Post3 = total(total(total(Post3,4),2),1) & L3 = total(total(total(L3,4),2),1) & P3 = total(total(total(P3,4),2),1)
          endif
        endif else begin
          if n4 gt 1 then begin
            Post3 = total(total(total(Post,4),2),1) & L3 = total(total(total(L,4),2),1) & P3 = total(total(total(P5D,4),2),1)
          endif else begin
            Post3 = total(total(Post,2),1) & L3 = total(total(L,2),1) & P3 = total(total(P5D,2),1)
          endelse
        endelse        
        fileout=outfolder+prefix+'_par3.eps'
        device,filename=fileout,/encaps,bits=8,/color
          if keyword_set(Tuffs_attenuation) then begin
            plot,minmax(F_grid),minmax([Post3,L3,P3]),$
               xtitle='Light absorbed in SF regions, F_grid',ytitle='Probalibity Densities',/nodata
            oplot,F_grid,Post3,color=54
            oplot,[1,1]*(F_grid[ind[2]])[0],!y.crange,color=54
            oplot,F_grid,L3,linestyle=2,color=254
            oplot,[1,1]*(F_grid[minchi2_ind[2]])[0],!y.crange,linestyle=2,color=254
            oplot,F_grid,P3,linestyle=1,color=54
          endif else begin
            plot,minmax(tauV_BC_grid),minmax([Post3,L3,P3]),$
               xtitle='Birth could attenuation, tau_V_BC',ytitle='Probalibity Densities',/nodata
            oplot,tauV_BC_grid,Post3,color=54
            oplot,[1,1]*(tauV_BC_grid[ind[2]])[0],!y.crange,color=54
            oplot,tauV_BC_grid,L3,linestyle=2,color=254
            oplot,[1,1]*(tauV_BC_grid[minchi2_ind[2]])[0],!y.crange,linestyle=2,color=254
            oplot,tauV_BC_grid,P3,linestyle=1,color=54
          endelse
        device,/close
      endif

      if n4 gt 1 then begin
        if n5 gt 1 then begin
          Post4 = total(total(total(total(Post,5),3),2),1)
          L4    = total(total(total(total(L,5),3),2),1)
          P4    = total(total(total(total(P5D,5),3),2),1)
        endif else begin
          Post4 = total(total(total(Post,3),2),1)
          L4    = total(total(total(L,3),2),1)
          P4    = total(total(total(P5D,3),2),1)          
        endelse
        fileout=outfolder+prefix+'_par4.eps'
        device,filename=fileout,/encaps,bits=8,/color
          plot,minmax(cosi_grid),minmax([Post4,L4,P4]),$
               xtitle='Cosine of inclination, cosi_grid',ytitle='Probalibity Densities',/nodata
            oplot,cosi_grid,Post4,color=54
            oplot,[1,1]*(cosi_grid[ind[3]])[0],!y.crange,color=54
            oplot,cosi_grid,L4,linestyle=2,color=254
            oplot,[1,1]*(cosi_grid[minchi2_ind[3]])[0],!y.crange,linestyle=2,color=254
            oplot,cosi_grid,P4,linestyle=1,color=54
        device,/close
      endif
      
      if n5 gt 1 then begin
        Post5 = total(total(total(total(Post,4),3),2),1)
        L5    = total(total(total(total(L,4),3),2),1)
        P5    = total(total(total(total(P5D,4),3),2),1)
        fileout=outfolder+prefix+'_par5.eps'
        device,filename=fileout,/encaps,bits=8,/color
          plot,minmax(rbulge_grid),minmax([Post5,L5,P5]),$
               xtitle='!8L!6!dbulge!n/!8L!6!dtotal!n, rbulge_grid',ytitle='Probalibity Densities',/nodata
            oplot,rbulge_grid,Post5,color=54
            oplot,[1,1]*(rbulge_grid[ind[4]])[0],!y.crange,color=54
            oplot,rbulge_grid,L5,linestyle=2,color=254
            oplot,[1,1]*(rbulge_grid[minchi2_ind[4]])[0],!y.crange,linestyle=2,color=254
            oplot,rbulge_grid,P5,linestyle=1,color=54
        device,/close
      endif
    endif
  endfor
endfor
set_plot,mydevice

bestfit = {Lnu:minchisqr_Lnu,Lnu_unred:minchisqr_Lnu_unred,LIR:minchisqr_LIR,$
           chisqr:minchisqr,ind:minchisqr_ind,coeff:coeff_minchisqr}
maxpost = {Lnu:Lnu_mod,Lnu_unred:Lnu_mod_unred,LIR:LIR_mod,$
           chisqr:chisqr_maxpost,ind:maxpost_ind,coeff:coeff_maxpost}

return

end



function lightning_models_Keith,steps=steps,grid=grid,Tuffs_attenuation=Tuffs_attenuation,exp_tau=exp_tau, $
                                Filters=Filters,print_time=print_time,time=time,LTIR_model_table=LTIR_model_table, $
			                    LTIR_table=LTIR_table,BC_exp_tau=BC_exp_tau,_REF_EXTRA=_e_models


;Modification History
; Keith Doore     - Replaced if statements with n_elements() on keywords to use keyword_set() - April 27th 2020
; Keith Doore     - Added use of birth cloud component on modified Calzetti if grid not set - April 27th 2020
; Keith Doore     - Added ability to input birth cloud component exp_tau in modified Calzetti - April 27th 2020
; Keith Doore     - Fixed issue with LTIR_table if a data point was at the final point in the table - April 27th 2020
; Keith Doore     - Added ability to use pure Calzetti curve - May 6th 2020
; Keith Doore     - Added _REF_EXTRA for keyword inheritance to cut down on list of keywords - May 6th 2020



if ~keyword_set(Tuffs_attenuation) and keyword_set(LTIR_model_table) then message,'LTIR_model_table keyword can only be called with Tuffs_attenuation'

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

if n_elements(grid) eq 0 then begin
  if keyword_set(Tuffs_attenuation) then begin
    n1 = 7 
    n2 = 7
    n3 = 7
    n4 = 7
    n5 = 1

    ran1=[ 0.d,8.0d]  & tauB_f_grid =ran1[0]+dindgen(n1)*(ran1[1]-ran1[0])/(n1-1) & if n1 eq 1 then tauB_f_grid=[0.d0]
    ran2=[ 0.d,1.0d]  & rdisk0_grid =ran2[0]+dindgen(n2)*(ran2[1]-ran2[0])/(n2-1) & if n2 eq 1 then rdisk0_grid=[0.d0]
    ran3=[0.d,0.61d]  & F_grid      =ran3[0]+dindgen(n3)*(ran3[1]-ran3[0])/(n3-1) & if n3 eq 1 then F_grid=[0.d0]
    ran4=[ 0.d,1.0d]  & cosi_grid   =ran4[0]+dindgen(n4)*(ran4[1]-ran4[0])/(n4-1) & if n4 eq 1 then cosi_grid=[0.d0]
    ran5=[ 0.d,1.0d]  & rbulge0_grid=ran5[0]+dindgen(n5)*(ran5[1]-ran5[0])/(n5-1) & if n5 eq 1 then rbulge0_grid=[0.d0]
    if (n_elements(exp_tau) eq 0) then                        $
      exp_tau = Tuffs_matrix_unred(wave_rest,                 $
                                   tauB_f_arr   = tauB_f_grid,$
                                   rdisk0_arr   = rdisk0_grid, $
                                   rbulge0_arr  = rbulge0_grid,$
                                   F_arr        = F_grid,     $
                                   cosi_arr     = cosi_grid,  $
                                   _EXTRA=_e_models)
  endif else begin
    n1 = 21
    n2 = 21
    n3 = 21
    n4 = 1
    n5 = 1

    ran1 = [   0d,4d0] & tauV_DIFF_grid = ran1[0] + dindgen(n1)*(ran1[1]-ran1[0])/(n1-1) ;as in Leja et al. 2016
    ran2 = [-2.2d,.4d] & delta_grid     = ran2[0] + dindgen(n2)*(ran2[1]-ran2[0])/(n2-1) ;as in Leja et al. 2016
    ran3 = [   0d,4d0] & tauV_BC_grid   = ran3[0] + dindgen(n3)*(ran3[1]-ran3[0])/(n3-1) ;as in Leja et al. 2016
    if (n_elements(BC_exp_tau)) eq 0 then                               $
      BC_exp_tau = matrix_unred_Keith(wave_rest,                        $
                                   tauV_DIFF_arr = -1.d0*tauV_DIFF_grid,$
                                   delta_arr     = delta_grid,          $
                                   tauV_BC_arr   = -1.d0*tauV_BC_grid,  $
                                   _EXTRA=_e_models)
    if (n_elements(exp_tau)) eq 0 then                                  $
      exp_tau    = matrix_unred_Keith(wave_rest,                        $
                                   tauV_DIFF_arr = -1.d0*tauV_DIFF_grid,$
                                   delta_arr     = delta_grid,          $
                                   tauV_BC_arr   = dblarr(n3),          $
                                   _EXTRA=_e_models)
  endelse
endif else begin
  if keyword_set(Tuffs_attenuation) then begin
    tauB_f_grid = grid.tauB_f_grid & n1 = n_elements(tauB_f_grid)
    rdisk0_grid = grid.rdisk0_grid & n2 = n_elements(rdisk0_grid)
    F_grid      = grid.F_grid      & n3 = n_elements(F_grid)
    cosi_grid   = grid.cosi_grid   & n4 = n_elements(cosi_grid)
      if n4 eq 1 and ~keyword_set(LTIR_model_table) then begin 
        n4=13
        ran4=[ 0.d,1.0d]  & cosi_grid  =ran4[0]+dindgen(n4)*(ran4[1]-ran4[0])/(n4-1)
        cosi_grid=([cosi_grid,grid.cosi_grid])[sort([cosi_grid,grid.cosi_grid])]
        n4 = n_elements(cosi_grid)
        for i=0,n4-1 do begin
          if cosi_grid[i] eq grid.cosi_grid then input_cosi=i
        endfor
      endif
    rbulge0_grid = grid.rbulge0_grid & n5 = n_elements(rbulge0_grid)
    if n_elements(exp_tau) eq 0 then $
      exp_tau=Tuffs_matrix_unred(wave_rest,               $
                                 tauB_f_arr = tauB_f_grid,$
                                 rdisk0_arr= rdisk0_grid, $
                                 rbulge0_arr=rbulge0_grid,$
                                 F_arr=F_grid,            $
                                 cosi_arr=cosi_grid,      $
                                 _EXTRA=_e_models)
  endif else begin
    tauV_DIFF_grid = grid.tauV_DIFF_grid & n1 = n_elements(tauV_DIFF_grid)
    delta_grid     = grid.delta_grid     & n2 = n_elements(delta_grid)
    tauV_BC_grid   = grid.tauV_BC_grid   & n3 = n_elements(tauV_BC_grid)
    n4 = 1
    n5 = 1
    if (n_elements(BC_exp_tau) eq 0) then $      
      BC_exp_tau = matrix_unred_Keith(wave_rest,tauV_DIFF_arr=-1.d0*tauV_DIFF_grid,delta_arr=delta_grid,$
                   tauV_BC_arr=-1.d0*tauV_BC_grid,_EXTRA=_e_models)
    if (n_elements(exp_tau) eq 0) then $     
      exp_tau    = matrix_unred_Keith(wave_rest,tauV_DIFF_arr=-1.d0*tauV_DIFF_grid,delta_arr=delta_grid,$
                   tauV_BC_arr=dblarr(n3),_EXTRA=_e_models)
  endelse
endelse

steps_mean_Lnu_red_grid = dblarr(Nfilters+1,Nsteps,n1,n2,n3,n4,n5)
transmission = transpose(Filters)
Lnu_inc = dblarr(Nsteps,n1,n2,n3,n5)

if keyword_set(LTIR_model_table) then begin
  if n_elements(LTIR_table) eq 0 then LTIR_table=mrdfits('~/Lightning/LTIR_model_table_z_'+string(z_shift,'(F4.2)')+'.fits',1)

  n1_LTIR_table = n_elements(LTIR_table.tauB_f_grid)
  n2_LTIR_table = n_elements(LTIR_table.rdisk0_grid)
  n3_LTIR_table = n_elements(LTIR_table.F_grid)
  n5_LTIR_table = n_elements(LTIR_table.rbulge0_grid)
  taub_loc   = interpol(lindgen(n1_LTIR_table),LTIR_table.tauB_f_grid,tauB_f_grid)>0
  rdisk0_loc  = interpol(lindgen(n2_LTIR_table),LTIR_table.rdisk0_grid,rdisk0_grid)>0
  F_loc      = interpol(lindgen(n3_LTIR_table),LTIR_table.F_grid,F_grid)>0
  rbulge0_loc = interpol(lindgen(n5_LTIR_table),LTIR_table.rbulge0_grid,rbulge0_grid)>0
  
  tau_weights=taub_loc-fix(taub_loc)
  rdisk0_weights=rdisk0_loc-fix(rdisk0_loc)
  f_weights=f_loc-fix(f_loc)
  rbulge0_weights=rbulge0_loc-fix(rbulge0_loc)
  
  ;Needed to allow for below process to work if n5_LTIR_table=1, which it currently is
  if n5_LTIR_table eq 1 then LTIR_model=rebin(LTIR_table.LABS_STAR_MODEL,LTIR_table.nsteps,n1_LTIR_table,n2_LTIR_table,n3_LTIR_table,2)

  
  for k1=0,(n1-1) do begin
    for k2=0,(n2-1) do begin
      for k3=0,(n3-1) do begin
        for k5=0,(n5-1) do begin
          
          if taub_loc[k1] eq (n1_LTIR_table-1) then taub_loc_temp=[taub_loc[k1],taub_loc[k1]] $
            else taub_loc_temp=[taub_loc[k1]:(taub_loc[k1]+1)]
          if rdisk0_loc[k2] eq (n2_LTIR_table-1) then rdisk0_loc_temp=[rdisk0_loc[k2],rdisk0_loc[k2]] $
            else rdisk0_loc_temp=[rdisk0_loc[k2]:(rdisk0_loc[k2]+1)]
          if f_loc[k3] eq (n3_LTIR_table-1) then f_loc_temp=[f_loc[k3],f_loc[k3]] $
            else f_loc_temp=[f_loc[k3]:(f_loc[k3]+1)] 
          if rbulge0_loc[k5] eq (n5_LTIR_table-1) then rbulge0_loc_temp=[rbulge0_loc[k5],rbulge0_loc[k5]] $
            else rbulge0_loc_temp=[rbulge0_loc[k5]:(rbulge0_loc[k5]+1)] 
            
          averaging_region=reform(LTIR_model[*,taub_loc_temp,rdisk0_loc_temp,f_loc_temp,rbulge0_loc_temp])
          temp1=reform(averaging_region[*,1,*,*,*]*tau_weights[k1]+averaging_region[*,0,*,*,*]*(1-tau_weights[k1]))
          temp2=reform(temp1[*,1,*,*]*rdisk0_weights[k2]+temp1[*,0,*,*]*(1-rdisk0_weights[k2]))
          temp3=reform(temp2[*,1,*]*f_weights[k3]+temp2[*,0,*]*(1-f_weights[k3]))
          Lnu_inc[*,k1,k2,k3,k5]=reform(temp3[*,1]*rbulge0_weights[k5]+temp3[*,0]*(1-rbulge0_weights[k5]))
        endfor
      endfor
    endfor
  endfor
endif

if ~keyword_set(LTIR_model_table) then Lnu_att_total = dblarr(Nsteps,n1,n2,n3,n4,n5)
for st=0,(Nsteps-1) do begin
  Lnu_5D=rebin(steps.lnu[*,st],nwave,n1,n2,n3,n4,n5)
  if ~keyword_set(Tuffs_attenuation) and st eq 0 then begin
    steps_Lnu_red = BC_exp_tau*Lnu_5D
  endif else begin
    steps_Lnu_red = exp_tau*Lnu_5D
  endelse
  ;Lbol = (1+z_shift)*trapez(lnu_5d,nu_rest)*(-1.d0) WRONG 1+z factor
  ;Lbol = trapez(lnu_5d,nu_rest)/(1+z_shift)*(-1.d0) ; correct 1+z factor, does not include nebular emission absorption
  ;Lbol = trapez(lnu_5d,nu_obs)*(-1.d0)              ; same as previous line
  Lbol = (steps.Lbol[st])[0]
  for k3=0,(n3-1) do begin
    for k2=0,(n2-1) do begin
      for k1=0,(n1-1) do begin
        for k5=0,(n5-1) do begin
          for k4=0,(n4-1) do begin
            for k=0,(Nfilters-1) do begin
              steps_mean_Lnu_red_grid[k,st,k1,k2,k3,k4,k5] = trapez(steps_Lnu_red[*,k1,k2,k3,k4,k5]*transmission[*,k],nu_obs)
            endfor
        	if ~keyword_set(LTIR_model_table) then $
        		Lnu_att_total[st,k1,k2,k3,k4,k5] = trapez(steps_Lnu_red[*,k1,k2,k3,k4,k5],nu_obs)
            ; gives negative value from integral bc nu_rest is inverted           
          endfor
      	  if keyword_set(Tuffs_attenuation) then begin
      	    if ~keyword_set(LTIR_model_table) then $
      	      Lnu_inc[st,k1,k2,k3,k5]=-1.d0*trapez((Lbol+Lnu_att_total[st,k1,k2,k3,*,k5])*sin(acos(cosi_grid)),acos(cosi_grid))
 		   	  ;integral is negative bc cosi grid is inverted (π/2 to 0) and Lnu_att_total is negative from above but added 
 		   	  ;to Lbol actually to subtract
      	    steps_mean_Lnu_red_grid[Nfilters,st,k1,k2,k3,*,k5] = Lnu_inc[st,k1,k2,k3,k5]
      	  endif
        endfor 
      endfor
    endfor
  endfor
  if ~keyword_set(Tuffs_attenuation) then begin
    steps_mean_Lnu_red_grid[Nfilters,st,*,*,*,*,*] = Lbol + Lnu_att_total[st,*,*,*,*,*]
  endif
endfor

steps_Lnu_red = !null

t2 = systime(/seconds)
if keyword_set(print_time) then print,'Total models running time: ',t2-t0,' seconds'
time=t2-t0
if n_elements(input_cosi) ne 0 then steps_mean_Lnu_red_grid=steps_mean_Lnu_red_grid[*,*,*,*,*,input_cosi[0],*]
models = temporary(steps_mean_Lnu_red_grid)
return,models

end



function lightning,Lnu_obs,Lnu_unc,LIR_obs,LIR_unc,Nsim=Nsim,steps=steps,grid=grid,models=models,exp_tau=exp_tau,$
                   float=float,regions=regions,Lnu_sim=Lnu_sim,LIR_sim=LIR_sim,Tuffs_attenuation=Tuffs_attenuation,$
                   lightning_folder=lightning_folder
;Modification History
; Keith Doore     - Created models function & added Tuffs extinction option - 2019
; Rafael Eufrasio - Corrected (1+z) factor - March 9th 2020
; Keith Doore     - Replaced if statements with n_elements() on keywords to use keyword_set() - April 27th 2020

t0 = SYSTIME(/seconds)

Nfilters = (Lnu_obs.DIM)[0]
if (SIZE(Lnu_obs))[0] eq 2 then Nreg=(Lnu_obs.DIM)[1] else Nreg=1
if (n_elements(regions) eq 0) then begin
  Nfit = Nreg
  regions = [0,Nreg-1]
endif else begin
  Nfit = Regions[1] - Regions[0] + 1
endelse
N0 = regions[0]
Nf = regions[1]

if (N_ELEMENTS(Nsim) eq 0) then Nsim=1
if Nsim eq 0 then Nsim=1

if ((Lnu_obs.Dim)[0] ne Nfilters) then begin
  message,'ERROR - SED and uncertainties must have same sizes'
  return,!null
endif

if ((Lnu_unc.Dim)[0] ne Nfilters) then begin
  message,'ERROR - SED and uncertainties must have same sizes'
  return,!null
endif

if (N_ELEMENTS(Lnu_sim) eq 0) then begin
  Lnu_sim = dblarr(Nfilters,Nreg,Nsim,/nozero)
  Lnu_sim[*,*,0] = Lnu_obs   ; first simulation is the observed data
  for ss=1,(Nsim-1) do Lnu_sim[*,*,ss] = Lnu_obs + Lnu_unc*randomn(seed,Nfilters,Nreg)
endif

if (N_ELEMENTS(LIR_sim) eq 0) then begin
  LIR_sim = dblarr(Nreg,Nsim)
  LIR_sim[*,0]   = LIR_obs   ; first simulation is the observed data
  if (Nsim gt 1) then begin
    rep_sim = replicate(1.d0,Nsim-1)
    LIR_sim[*,1:*] = (LIR_obs # rep_sim) + (LIR_unc # rep_sim)*randomn(seed,Nreg,Nsim)
  endif
endif

Lsim = dblarr(Nfilters+1,Nfit,Nsim,/nozero)
Lunc = dblarr(Nfilters+1,Nfit,/nozero)  
  Lsim[0:(Nfilters-1),*,*] = Lnu_sim[*,N0:Nf,*] ; First entries are broadband filters
  Lunc[0:(Nfilters-1),*]   = Lnu_unc[*,N0:Nf]   ; First entries are broadband filters
  Lsim[Nfilters,*,*]       = LIR_sim[N0:Nf,*]   ; Last entry is TIR
  Lunc[Nfilters,*]         = LIR_unc[N0:Nf]     ; Last entry is TIR

if (N_ELEMENTS(steps) eq 0) then begin
  steps = STEPS_STELLAR(filter_labels=filter_labels,steps_bounds=steps_bounds,z_shift=0.0,dtime_SF=5.d5)
endif
Nsteps = steps.bounds.length - 1

wave_rest = steps.wave_rest ; restframe wavelength
wave_obs  = steps.wave_obs  ; observed wavelength
Nwave = wave_rest.length

nu_rest = 1.d4 * !cv.clight / wave_rest
nu_obs  = 1.d4 * !cv.clight / wave_obs
z_shift = steps.z_shift

;normalizing filters
Filters = steps.filters
for k=0,(Nfilters-1) do Filters[k,*] /= (trapez(Filters[k,*],nu_obs))[0]

if n_elements(grid) eq 0 then begin
  if keyword_set(Tuffs_attenuation) then begin
      n1 = 7 
      n2 = 7
      n3 = 7
      n4 = 7
      n5 = 1
    
      ran1=[ 0.d,8.0d]  & tauB_f_grid=ran1[0]+dindgen(n1)*(ran1[1]-ran1[0])/(n1-1) & if n1 eq 1 then tauB_f_grid=[0.d0]
      ran2=[ 0.d,1.0d]  & rdisk_grid =ran2[0]+dindgen(n2)*(ran2[1]-ran2[0])/(n2-1) & if n2 eq 1 then rdisk_grid=[0.d0]
      ran3=[0.d,0.61d]  & F_grid     =ran3[0]+dindgen(n3)*(ran3[1]-ran3[0])/(n3-1) & if n3 eq 1 then F_grid=[0.d0]
      ran4=[ 0.d,1.0d]  & cosi_grid  =ran4[0]+dindgen(n4)*(ran4[1]-ran4[0])/(n4-1) & if n4 eq 1 then cosi_grid=[0.d0]
      ran5=[ 0.d,1.0d]  & rbulge_grid=ran5[0]+dindgen(n5)*(ran5[1]-ran5[0])/(n5-1) & if n5 eq 1 then rbulge_grid=[0.d0]
  endif else begin
    n1 = 21
      n2 = 21
      n3 = 21
      n4 = 1
      n5 = 1
      
      ran1 = [   0d,4d0] & tauV_DIFF_grid = ran1[0] + dindgen(n1)*(ran1[1]-ran1[0])/(n1-1) ;as in Leja et al. 2016
      ran2 = [-2.2d,.4d] & delta_grid     = ran2[0] + dindgen(n2)*(ran2[1]-ran2[0])/(n2-1) ;as in Leja et al. 2016
      ran3 = [   0d,4d0] & tauV_BC_grid   = ran3[0] + dindgen(n3)*(ran3[1]-ran3[0])/(n3-1) ;as in Leja et al. 2016
  endelse
endif else begin
  if keyword_set(Tuffs_attenuation) then begin
    tauB_f_grid = grid.tauB_f_grid & n1 = n_elements(tauB_f_grid)
    rdisk_grid  = grid.rdisk_grid  & n2 = n_elements(rdisk_grid)
    F_grid      = grid.F_grid      & n3 = n_elements(F_grid)
    cosi_grid   = grid.cosi_grid   & n4 = n_elements(cosi_grid)
    rbulge_grid = grid.rbulge_grid & n5 = n_elements(rbulge_grid)
  endif else begin
    tauV_DIFF_grid = grid.tauV_DIFF_grid & n1 = n_elements(tauV_DIFF_grid)
      delta_grid     = grid.delta_grid     & n2 = n_elements(delta_grid)
      tauV_BC_grid   = grid.tauV_BC_grid   & n3 = n_elements(tauV_BC_grid)
      n4 = 1
      n5 = 1
  endelse
endelse

if n_elements(models) eq 0 then begin

  if keyword_set(Tuffs_attenuation) then begin
    steps_mean_Lnu_red_grid=lightning_models_Keith(steps=steps,grid=grid,/Tuffs,exp_tau=exp_tau, $
      lightning_folder=lightning_folder,Filters=Filters,time=time)
  endif else begin
    steps_mean_Lnu_red_grid=lightning_models_Keith(steps=steps,grid=grid,exp_tau=exp_tau,Filters=Filters,time=time)
  endelse
  
  print,'Time to generate models : ',time,' seconds'
endif else begin
  steps_mean_Lnu_red_grid = models
endelse


if keyword_set(float) then begin
  coeff_grid = fltarr(Nsteps,n1,n2,n3,n4,n5,Nfit,Nsim,/nozero)    & replicate_inplace,coeff_grid,1.e0
  chisqr_grid = fltarr(n1,n2,n3,n4,n5,Nfit,Nsim,/nozero)          & replicate_inplace,chisqr_grid,1.e0
  Lmod_grid = fltarr(Nfilters+1,n1,n2,n3,n4,n5,Nfit,Nsim,/nozero) & replicate_inplace,Lmod_grid,0.e0
  coeff_temp = fltarr(Nfilters+1,Nsteps,n1,n2,n3,n4,n5,/nozero)   & replicate_inplace,coeff_temp,1.e0
  Lfit = fltarr(Nfilters+1,n1,n2,n3,n4,n5,/nozero)                & replicate_inplace,Lfit,1.e0
  Lsig = fltarr(Nfilters+1,n1,n2,n3,n4,n5,/nozero)                & replicate_inplace,Lsig,1.e0
endif else begin
  coeff_grid = dblarr(Nsteps,n1,n2,n3,n4,n5,Nfit,Nsim,/nozero)    & replicate_inplace,coeff_grid,0.d0
  chisqr_grid = dblarr(n1,n2,n3,n4,n5,Nfit,Nsim,/nozero)          & replicate_inplace,chisqr_grid,1.d0
  Lmod_grid = dblarr(Nfilters+1,n1,n2,n3,n4,n5,Nfit,Nsim,/nozero) & replicate_inplace,Lmod_grid,0.d0
  coeff_temp = dblarr(Nfilters+1,Nsteps,n1,n2,n3,n4,n5,/nozero)   & replicate_inplace,coeff_temp,1.d0
  Lfit = dblarr(Nfilters+1,n1,n2,n3,n4,n5,/nozero)                & replicate_inplace,Lfit,1.d0
  Lsig = dblarr(Nfilters+1,n1,n2,n3,n4,n5,/nozero)                & replicate_inplace,Lsig,1.d0
endelse

t1=systime(/sec)
for sr=0,(Nfit-1) do begin
  Lreg_tmp = Lsim[*,sr,*]
  Lsig_tmp = Lunc[*,sr]
  for ss=0,(Nsim-1) do begin
    Lsim_tmp = Lreg_tmp[*,0,ss]
    for k=0,Nfilters do begin
      Lfit[k,*,*,*,*,*] = Lsim_tmp[k]
      Lsig[k,*,*,*,*,*] = Lsig_tmp[k]
    endfor 
    for k5=0,(n5-1) do begin 
      for k4=0,(n4-1) do begin
        for k3=0,(n3-1) do begin
          for k2=0,(n2-1) do begin
            for k1=0,(n1-1) do begin
              coeff_grid[*,k1,k2,k3,k4,k5,sr,ss] = $
                steps_derivative_chisqr_minimize(Lsim_tmp,Lsig_tmp,steps_mean_Lnu_red_grid[*,*,k1,k2,k3,k4,k5])
            endfor
          endfor
        endfor
      endfor
    endfor
    coeff = coeff_grid[*,*,*,*,*,*,sr,ss]
    for k=0,Nfilters do coeff_temp[k,*,*,*,*,*,*] = coeff
    Lmod_grid[*,*,*,*,*,*,sr,ss] = total(coeff_temp * steps_mean_Lnu_red_grid,2)
    chisqr_grid[*,*,*,*,*,sr,ss] = total(((Lmod_grid[*,*,*,*,*,*,sr,ss] - Lfit)/Lsig)^2,1,/NaN)
  endfor
  ti=systime(/seconds)
                                             nR=string(ceil(alog10(Nreg+1)))
  dt = ti-t1                               & nT=string(ceil(alog10(dt)+2) > 6)
  minchisqr = min(chisqr_grid[*,*,*,*,*,sr,0]) & nC=string(ceil(alog10(minchisqr)+3) > 7)
  print,string((sr+N0+1.)*100./Nreg,form='(F5.1)')+'% completed --- Done with region # '+$
        string(sr+N0+1,form='(I'+nR+')')+' out of '+strtrim(Nreg,2)+' @ '+$
        string(dt,form='(F'+nT+'.1)')+' seconds --- Min chisqr: '+string(minchisqr,form='(F'+nC+'.2)')  
endfor

t2 = systime(/seconds)
print,'Main loop running time: ',t2-t1,' seconds'
print,'Total running time: ',t2-t0,' seconds'

models = temporary(steps_mean_Lnu_red_grid)

output = {chisqr_grid:temporary(chisqr_grid),$
          coeff_grid:temporary(coeff_grid),$
          Lmod_grid:temporary(Lmod_grid)}

return,output

end














pro lightning_MCMC_RTE,Lnu_obs,Lnu_unc,LTIR_obs,LTIR_unc,galaxy_id,filter_labels,z_shift=z_shift,steps_bounds=steps_bounds,$
                   print_time=print_time,Tuffs_attenuation=Tuffs_attenuation,$
                   par1_constant=par1_constant,par2_constant=par2_constant,par3_constant=par3_constant,$
                   par4_constant=par4_constant,par5_constant=par5_constant,par6_constant=par6_constant,$
                   par7_constant=par7_constant,par8_constant=par8_constant,par9_constant=par9_constant,$
                   par10_constant=par10_constant,Ntrials=Ntrials,outfolder=outfolder,beta_exponent=beta_exponent,$
                   mu_input=mu_input,sigma_input=sigma_input,lambda_input=lambda_input,$
                   adaptive=adaptive,coeff_start=coeff_start,coeff_sigma=coeff_sigma,$
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
  
  if (n_elements(parameter_start) eq 0) and keyword_set(dust_emission) and ~keyword_set(Tuffs_attenuation) then parameter_start=[0.2,0.0,0.2,0,0,-2.d0,1.0,1.d4,0.1,0.020]
  if (n_elements(parameter_sigma) eq 0) and keyword_set(dust_emission) and ~keyword_set(Tuffs_attenuation) then parameter_sigma=[0.1,0.1,0.2,0,0,0.5d0,0.1,1.d4,0.1,0.005]
  if (n_elements(parameter_start) eq 0) and keyword_set(dust_emission) and keyword_set(Tuffs_attenuation) then parameter_start=[1.0,0.2,0.2,0.8,0.2,-2.d0,1.0,1.d4,0.1,0.020]
  if (n_elements(parameter_sigma) eq 0) and keyword_set(dust_emission) and keyword_set(Tuffs_attenuation) then parameter_sigma=[0.5,0.1,0.1,0.1,0.1,0.5d0,0.1,1.d4,0.1,0.005]

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

  if keyword_set(Tuffs_attenuation) then begin
    ran1=[ 0.d,8.0d]
    ran2=[ 0.d,1.0d]
    ran3=[0.d,0.61d]
    ran4=[ 0.d,1.0d]
    ran5=[ 0.d,1.0d]
  
    par1_guess = parameter_start[0]  &  par1_sigma = parameter_sigma[0]  &  if keyword_set(par1_constant) then par1_sigma=0.d0
    par2_guess = parameter_start[1]  &  par2_sigma = parameter_sigma[1]  &  if keyword_set(par2_constant) then par2_sigma=0.d0
    par3_guess = parameter_start[2]  &  par3_sigma = parameter_sigma[2]  &  if keyword_set(par3_constant) then par3_sigma=0.d0
    par4_guess = parameter_start[3]  &  par4_sigma = parameter_sigma[3]  &  if keyword_set(par4_constant) then par4_sigma=0.d0
    par5_guess = parameter_start[4]  &  par5_sigma = parameter_sigma[4]  &  if keyword_set(par5_constant) then par5_sigma=0.d0
    
    grid = {tauB_f_grid:par1_guess, rdisk_grid:par2_guess, F_grid:par3_guess, cosi_grid:par4_guess, rbulge_grid:par5_guess}
    models=lightning_models_Keith(steps=steps,grid=grid,Filters=Filters,/tuffs,_EXTRA=_extra_MCMC)
    models[where(finite(models,/NaN),/null)] = 0.0
  endif else begin
    ran1=[ 0.0d,3.0d]
    ran2=[-2.3d,0.4d]
    ran3=[ 0.0d,4.0d]
    ran4=[ 0.0d,0.0d]
    ran5=[ 0.0d,0.0d]

    par1_guess = parameter_start[0] & par1_sigma = parameter_sigma[0] & if keyword_set(par1_constant) then par1_sigma=0.d0
    par2_guess = parameter_start[1] & par2_sigma = parameter_sigma[1] & if keyword_set(par2_constant) then par2_sigma=0.d0
    par3_guess = parameter_start[2] & par3_sigma = parameter_sigma[2] & if keyword_set(par3_constant) then par3_sigma=0.d0
    par4_guess = 0.0  &  par4_sigma = 0.0d0  &  par4_constant=1                             
    par5_guess = 0.0  &  par5_sigma = 0.0d0  &  par5_constant=1   

    grid   = {tauV_DIFF_grid:par1_guess, delta_grid:par2_guess, tauV_BC_grid:par3_guess}
    models = lightning_models_Keith(steps=steps,grid=grid,Filters=Filters,_EXTRA=_extra_MCMC)
    models[where(finite(models,/NaN),/null)] = 0.0
  endelse


  if keyword_set(dust_emission) then begin
    ran6 = [-10.d0,  4d0]       ; alpha
    ran7 = [ 0.1d0,25]        ; Umin
    ran8 = [ 1.d3,  3.d5]       ; Umax
    ran9 = [ 0.d0,  1d0]        ; gamma
    ran10= [ 0.0047,0.0458d0]  ; qPAH

    par6_guess = parameter_start[5] & par6_sigma = parameter_sigma[5] & if keyword_set(par6_constant)  then par6_sigma=0.d0
    par7_guess = parameter_start[6] & par7_sigma = parameter_sigma[6] & if keyword_set(par7_constant)  then par7_sigma=0.d0
    par8_guess = parameter_start[7] & par8_sigma = parameter_sigma[7] & if keyword_set(par8_constant)  then par8_sigma=0.d0
    par9_guess = parameter_start[8] & par9_sigma = parameter_sigma[8] & if keyword_set(par9_constant)  then par9_sigma=0.d0
    par10_guess= parameter_start[9] & par10_sigma= parameter_sigma[9] & if keyword_set(par10_constant) then par10_sigma=0.d0
    
    dl07     = dl07_templates(filter_labels=filter_labels,z_shift=z_shift,_EXTRA=_extra_MCMC)
    de_model = dl07_sed(dl07,alpha=par6_guess,umin=par7_guess,umax=par8_guess,gam=par9_guess,qPAH=par10_guess,filter_labels=filter_labels,Lbol=de_LTIR)
  endif


  if (n_elements(coeff_start) eq 0) then coeff_start = replicate(1.d0,nsteps)
  if (n_elements(coeff_start) ne nsteps) then message, 'ERROR -- Number of starting SFH coefficients must equal number of steps'
  if (n_elements(coeff_sigma) eq 0) then coeff_sigma = replicate(1.d0,nsteps)
  if (n_elements(coeff_sigma) ne nsteps) then message, 'ERROR -- Number of starting SFH coefficients must equal number of steps'

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
  Lmod_old  = total(models * (rebin(reform(coeff_old,1,Nsteps),Nfilters+1,Nsteps)),2) ;models(Nfilters+1,Nsteps)
  ;Lmod_old = total(models * (rebin(reform(coeff_old,1,n_elements(coeff_old)),Nfilters+1,n_elements(coeff_old))),2)
  Lmod_old[0:(Nfilters-1)] = Lmod_old[0:(Nfilters-1)] + (Lmod_old[Nfilters])[0]*de_model/de_LTIR
  chi2_old = total(((Lmod_old - Lmod_obs)/Lmod_unc)^2,/NaN)

  chain_old  = [coeff_old,par1_old,par2_old,par3_old,par4_old,par5_old,par6_old,par7_old,par8_old,par9_old,par10_old]
  sigma_vals = [coeff_sigma,par1_sigma,par2_sigma,par3_sigma,par4_sigma,par5_sigma,par6_sigma,par7_sigma,par8_sigma,par9_sigma,par10_sigma]
  chain = [chain_old]
  chi2_chain = [chi2_old]
  constant_param = where(sigma_vals eq 0,ncomp=Npar,comp=variable_param,/null)

  if (n_elements(lambda_input) eq 0) then lambda_new=2.38^2.d0/Npar else lambda_new=lambda_input
  ;if (n_elements(mu_input) eq 0) then mu_new=chain else mu_new=mu_input
  if (n_elements(mu_input) eq 0) then mu_new=chain+sigma_vals else mu_new=mu_input
  if (n_elements(sigma_input) eq 0) then sigma_new=diag_matrix(sigma_vals[variable_param]) else sigma_new=sigma_input

  p_jump=[0.441,0.352,0.316,0.279,0.275,0.266,0.261,0.255,0.261,0.267]
  alpha_star=p_jump[(Npar-1) < (p_jump.length-1)]


  for trial=0,(Ntrials-1) do begin
    gamma_new  = 1.d0/(trial+1.d0)^beta_exponent
    lambda_old = lambda_new
    mu_old     = mu_new
    sigma_old  = sigma_new
    covar      = lambda_old*sigma_old
    jump       = mrandomn(seed,covar,status=status)
    while (status gt 0) do begin
      ;evals[where(evals le 1.0d-200,/null)] = 1.0d-200
;      stop
      evals = eigenql(covar,EIGENVECTORS = evecs)
      evals[where(evals le 1.0d-16,/null)] = 1.0d-16
      ;evals[where(evals le 1.0d-6,/null)] = 1.0d-6
      ;evals[where(evals le 0.D0,/null)] = 1.0d-6
      covar = evecs # diag_matrix(evals) # invert(evecs)
      covar = (covar + transpose(covar))/2.d0
      jump  = mrandomn(seed,covar,status=status)
    endwhile
    ;stop
    chain_new = chain_old
    chain_new[variable_param] = chain_old[variable_param] + jump
    in = 1
    coeff_new = chain_new[0:(Nsteps-1)] & for st=0,(nsteps-1) do in = in and (coeff_new[st] ge 0)
    par1_new = chain_new[Nsteps  ] & in = in and ((par1_new  ge min(ran1 )) and (par1_new  le max(ran1 )))
    par2_new = chain_new[Nsteps+1] & in = in and ((par2_new  ge min(ran2 )) and (par2_new  le max(ran2 )))
    par3_new = chain_new[Nsteps+2] & in = in and ((par3_new  ge min(ran3 )) and (par3_new  le max(ran3 )))
    par4_new = chain_new[Nsteps+3] & in = in and ((par4_new  ge min(ran4 )) and (par4_new  le max(ran4 )))
    par5_new = chain_new[Nsteps+4] & in = in and ((par5_new  ge min(ran5 )) and (par5_new  le max(ran5 )))
    par6_new = chain_new[Nsteps+5] & in = in and ((par6_new  ge min(ran6 )) and (par6_new  le max(ran6 )))
    par7_new = chain_new[Nsteps+6] & in = in and ((par7_new  ge min(ran7 )) and (par7_new  le max(ran7 )))
    par8_new = chain_new[Nsteps+7] & in = in and ((par8_new  ge min(ran8 )) and (par8_new  le max(ran8 )))
    par9_new = chain_new[Nsteps+8] & in = in and ((par9_new  ge min(ran9 )) and (par9_new  le max(ran9 )))
    par10_new= chain_new[Nsteps+9] & in = in and ((par10_new ge min(ran10)) and (par10_new le max(ran10)))
    ;draw=1
    while in eq 0 do begin
    ;  draw++
      jump  = mrandomn(seed,covar,status=status)
      chain_new = chain_old
      chain_new[variable_param] = chain_old[variable_param] + jump
      in = 1
      coeff_new = chain_new[0:(Nsteps-1)] & for st=0,(nsteps-1) do in = in and (coeff_new[st] ge 0)
      par1_new = chain_new[Nsteps  ] & in = in and ((par1_new  ge min(ran1 )) and (par1_new  le max(ran1 )))
      par2_new = chain_new[Nsteps+1] & in = in and ((par2_new  ge min(ran2 )) and (par2_new  le max(ran2 )))
      par3_new = chain_new[Nsteps+2] & in = in and ((par3_new  ge min(ran3 )) and (par3_new  le max(ran3 )))
      par4_new = chain_new[Nsteps+3] & in = in and ((par4_new  ge min(ran4 )) and (par4_new  le max(ran4 )))
      par5_new = chain_new[Nsteps+4] & in = in and ((par5_new  ge min(ran5 )) and (par5_new  le max(ran5 )))
      par6_new = chain_new[Nsteps+5] & in = in and ((par6_new  ge min(ran6 )) and (par6_new  le max(ran6 )))
      par7_new = chain_new[Nsteps+6] & in = in and ((par7_new  ge min(ran7 )) and (par7_new  le max(ran7 )))
      par8_new = chain_new[Nsteps+7] & in = in and ((par8_new  ge min(ran8 )) and (par8_new  le max(ran8 )))
      par9_new = chain_new[Nsteps+8] & in = in and ((par9_new  ge min(ran9 )) and (par9_new  le max(ran9 )))
      par10_new= chain_new[Nsteps+9] & in = in and ((par10_new ge min(ran10)) and (par10_new le max(ran10)))
    ;  print,'draw = ',draw,chain_new
    endwhile
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
      grid = {tauB_f_grid:par1_new, rdisk_grid:par2_new, F_grid:par3_new, cosi_grid:par4_new, rbulge_grid:par5_new}
      models=lightning_models_Keith(steps=steps,grid=grid,Filters=Filters,/tuffs,_EXTRA=_extra_MCMC)
      models[where(finite(models,/NaN),/null)] = 0.0
    endif else begin
      grid = {tauV_DIFF_grid:par1_new, delta_grid:par2_new, tauV_BC_grid:par3_new}
      models = lightning_models_Keith(steps=steps,grid=grid,Filters=Filters,_EXTRA=_extra_MCMC)
      models[where(finite(models,/NaN),/null)] = 0.0
    endelse
    if keyword_set(dust_emission) then begin
      de_model = dl07_sed(dl07,alpha=par6_new,umin=par7_new,umax=par8_new,gam=par9_new,qPAH=par10_new,filter_labels=filter_labels,Lbol=de_LTIR)
    endif
    chain_new = [coeff_new,par1_new,par2_new,par3_new,par4_new,par5_new,par6_new,par7_new,par8_new,par9_new,par10_new]          
    Lmod_new  = total(models * (rebin(reform(coeff_new,1,Nsteps),Nfilters+1,Nsteps)),2) ;models(Nfilters+1,Nsteps)
    ;Lmod_new = total(models * (rebin(reform(coeff_new,1,n_elements(coeff_new)),Nfilters+1,n_elements(coeff_new))),2)
    Lmod_new[0:(Nfilters-1)] = Lmod_new[0:(Nfilters-1)] + (Lmod_new[-1])[0]*de_model/de_LTIR
    chi2_new = total(((Lmod_new - Lmod_obs)/Lmod_unc)^2,/NaN)       
    ratio = exp(-(chi2_new-chi2_old)/2.d0)
    if ratio gt randomu(seed,1) then begin                   
      chain = [[chain],[chain_new]]                        
      chi2_chain = [chi2_chain,chi2_new]    
      chain_old = chain_new            
      chi2_old = chi2_new
    endif

    ;print,trial,chi2_new,chain_old

    if keyword_set(adaptive) then begin
;stop
      lambda_new = 10.d0^(alog10(lambda_old)+gamma_new*(min([1,ratio])-alpha_star))
      mu_new=mu_old+gamma_new*(chain_old-mu_old)
      sigma_new=sigma_old+$
                gamma_new*((chain_old[variable_param]-mu_old[variable_param])#(chain_old[variable_param]-mu_old[variable_param])-sigma_old)
;stop
    endif
;    stop
  endfor

  Nchain = (chain.dim)[1]
  steps_Mstar = steps.Mstar
  
  ;if n_elements(Tuffs_attenuation) eq 0 then chain=chain[0:(nsteps+2),*]
  
  save,chain,chi2_chain,galaxy_id,sigma_new,lambda_new,mu_new,steps_Mstar,steps_bounds, $
       filename=outfolder+galaxy_id+'_chain.idl'
     
  if keyword_set(print_time) then begin
    t1 = systime(/sec)
    dt = t1-t0
    nT=string(ceil(alog10(dt)+2) > 6)
    nC=string(ceil(alog10(nchain+1)))
    
    print,''
      print,'Done in '+string(dt,form='(F'+nT+'.1)')+' seconds'
      print,'Number of accepted trials: '+string(nchain,form='(I'+nC+')')
      print,'Percentage of total trials: '+string(nchain/ntrials*100.d0,form='(F5.1)')+' %'
      print,''
  endif

end









pro lightning_MCMC_Keith,Lnu_obs,Lnu_unc,LTIR_obs,LTIR_unc,galaxy_id,filter_labels,z_shift=z_shift,steps_bounds=steps_bounds,$
                   print_time=print_time,Tuffs_attenuation=Tuffs_attenuation,lightning_folder=lightning_folder,NGC=NGC,$
                   folder_Tuffs_attenuation=folder_Tuffs_attenuation,par1_constant=par1_constant,par2_constant=par2_constant,$
                   par3_constant=par3_constant, par4_constant=par4_constant, par5_constant=par5_constant, Ntrials=Ntrials,$
                   outfolder_plots=outfolder_plots,beta_exponent=beta_exponent,mu_input=mu_input,sigma_input=sigma_input,$
                   lambda_input=lambda_input,dtime_SF=dtime_SF,nolines=nolines,nonebular=nonebular,Zmet=Zmet,IMF_type=IMF_type,$
                   adaptive=adaptive,coeff_start=coeff_start,parameter_start=parameter_start


  t0=systime(/sec)
  
  if typename(galaxy_id) ne 'STRING' then message,'ERROR -- galaxy_id input must be of type STRING'
  if (Lnu_obs.length ne filter_labels.length) then message,'ERROR -- SED and Filter_labels must have same sizes'
  if (Lnu_obs.length ne Lnu_unc.length) then message,'ERROR -- SED and uncertainties must have same sizes'
  if (LTIR_obs.length ne LTIR_unc.length) then message,'ERROR -- TIR and uncertainties must have same sizes'
  
  if (n_elements(z_shift) eq 0) then z_shift = 0.0
  if (n_elements(steps_bounds) eq 0) then steps_bounds = [0.d0,1.d7,1.d8,1.d9,5.d9,13.3d9]
  ;if (n_elements(lightning_folder) eq 0) then lightning_folder = '~/Box/Keith/Lightning/'
  if (n_elements(lightning_folder) eq 0) then lightning_folder = '~/Lightning/'
  if (n_elements(folder_Tuffs_attenuation) eq 0) then folder_Tuffs_attenuation='~/Box/Keith/Ellipticity/Data/'
  if (n_elements(Ntrials) eq 0) then Ntrials = 1e5
  ;if (n_elements(outfolder_plots) eq 0) then outfolder_plots='~/Box/Keith/MCMC_Lightning_runs/'
  if (n_elements(outfolder_plots) eq 0) then outfolder_plots='~/MCMC_Lightning_runs/'
  if (n_elements(beta_exponent) eq 0) then beta_exponent = 0.35
  if (n_elements(dtime_SF) eq 0) then dtime_SF=!null
  if ~keyword_set(nolines) then nolines=!null
  if ~keyword_set(nonebular) then nonebular=!null
  if (n_elements(Zmet) eq 0) then Zmet=!null
  if (n_elements(IMF_type) eq 0) then IMF_type=!null
  if (n_elements(parameter_start) eq 0) and keyword_set(Tuffs_attenuation)  then parameter_start=fltarr(5)
  if (n_elements(parameter_start) eq 0) and ~keyword_set(Tuffs_attenuation) then parameter_start=fltarr(3)

  if (n_elements(parameter_start) ne 5) and keyword_set(Tuffs_attenuation)  then message, 'ERROR -- Must have five starting parameters for Tuffs extinction'
  if (n_elements(parameter_start) ne 3) and ~keyword_set(Tuffs_attenuation) then message, 'ERROR -- Must have three starting parameters for modified Calzetti extinction'
  
  steps_bounds = steps_bounds < galage(z_shift,25,/silent)
    steps_bounds = steps_bounds[unique(steps_bounds)]
  Nsteps = n_elements(steps_bounds) - 1
  steps = STEPS_STELLAR(filter_labels=filter_labels,steps_bounds=steps_bounds,$
                        dtime_SF=dtime_SF,nolines=nolines,nonebular=nonebular,$
                        Zmet=Zmet,IMF_type=IMF_type,z_shift=z_shift,lightning_folder=lightning_folder)


  Nfilters = filter_labels.length
  Filters = steps.filters
  for k=0,(Nfilters-1) do Filters[k,*] /= (trapez(Filters[k,*],(1.d4*!cv.clight/steps.wave_obs)))[0]

  Lmod_obs = [Lnu_obs,LTIR_obs]
  Lmod_unc = [Lnu_unc,LTIR_unc]

  if keyword_set(Tuffs_attenuation) then begin
    ran1=[ 0.d,8.0d]
    ran2=[ 0.d,1.0d]
    ran3=[0.d,0.61d]
    ran4=[ 0.d,1.0d]
    ran5=[ 0.d,1.0d]
  
    par1_guess = parameter_start[0]  &  par1_sigma = 0.5d0  &   if keyword_set(par1_constant) then par1_sigma=0.d0
    par2_guess = parameter_start[1]  &  par2_sigma = 0.1d0  &   if keyword_set(par2_constant) then par2_sigma=0.d0
    par3_guess = parameter_start[2]  &  par3_sigma = 0.1d0  &   if keyword_set(par3_constant) then par3_sigma=0.d0
    par4_guess = parameter_start[3]  &  par4_sigma = 0.1d0  &   if keyword_set(par4_constant) then par4_sigma=0.d0
    par5_guess = parameter_start[4]  &  par5_sigma = 0.1d0  &   if keyword_set(par5_constant) then par5_sigma=0.d0
    
    grid = {tauB_f_grid:par1_guess, rdisk_grid:par2_guess, F_grid:par3_guess, cosi_grid:par4_guess, rbulge_grid:par5_guess}
    models=lightning_models_Keith(steps=steps,grid=grid,Filters=Filters,/Tuffs,coeff=folder_Tuffs_attenuation)
  endif else begin
    ran1=[ 0.0d,3.0d]                     
    ran2=[-2.3d,0.4d]                     
    ran3=[ 0.0d,4.0d]                     
    ran4=[ 0.0d,0.0d]                     
    ran5=[ 0.0d,0.0d]                     
    
    par1_guess = parameter_start[0]  &  par1_sigma = 0.1d0  &   if keyword_set(par1_constant) then par1_sigma=0.d0
    par2_guess = parameter_start[1]  &  par2_sigma = 0.1d0  &   if keyword_set(par2_constant) then par2_sigma=0.d0
    par3_guess = parameter_start[2]  &  par3_sigma = 0.2d0  &   if keyword_set(par3_constant) then par3_sigma=0.d0
    par4_guess = 0.0  &  par4_sigma = 0.0d0  &  par4_constant=1                             
    par5_guess = 0.0  &  par5_sigma = 0.0d0  &  par5_constant=1   

    grid = {tauV_DIFF_grid:par1_guess, delta_grid:par2_guess, tauV_BC_grid:par3_guess}
    models=lightning_models_Keith(steps=steps,grid=grid,Filters=Filters)
  endelse

  if keyword_set(par1_constant) then random2=0 else random2=1
  if keyword_set(par2_constant) then random3=random2 else random3=1+random2
  if keyword_set(par3_constant) then random4=random3 else random4=1+random3
  if keyword_set(par4_constant) then random5=random4 else random5=1+random4
  
  if (n_elements(coeff_start) eq 0) then coeff_start=replicate(1.d0,nsteps)
  if (n_elements(coeff_start) ne nsteps) then message, 'ERROR -- Number of starting SFH coefficients must equal number of steps'

  coeff_guess = coeff_start & coeff_sigma = replicate(1.d0,nsteps)
  
  coeff_old = coeff_guess                     
  par1_old = par1_guess
  par2_old = par2_guess
  par3_old = par3_guess
  par4_old = par4_guess
  par5_old = par5_guess     
  Lmod_old = total(models * (rebin(reform(coeff_old,1,n_elements(coeff_old)),Nfilters+1,n_elements(coeff_old))),2)
  chi2_old = total(((Lmod_old - Lmod_obs)/Lmod_unc)^2,/NaN)       
                                            
  chain_old = [coeff_old,par1_old,par2_old,par3_old,par4_old,par5_old]
  sigma_vals = [coeff_sigma,par1_sigma,par2_sigma,par3_sigma,par4_sigma,par5_sigma]
  chain = [chain_old]
  chi2_chain = [chi2_old]
  Constant_param=where(sigma_vals eq 0,ncomp=Npar,comp=variable_param,/null)
  if (n_elements(lambda_input) eq 0) then lambda_new=2.38^2.d0/Npar else lambda_new=lambda_input
  if (n_elements(mu_input) eq 0) then mu_new=chain else mu_new=mu_input
  if (n_elements(sigma_input) eq 0) then sigma_new=diag_matrix(sigma_vals[variable_param]) else sigma_new=sigma_input

  p_jump=[0.441,0.352,0.316,0.279,0.275,0.266,0.261,0.255,0.261,0.267]
  alpha_star=p_jump[(Npar-1)]

  for trial=0,(Ntrials-1) do begin
    gamma_new=1/(trial+1)^beta_exponent
    lambda_old=lambda_new
    mu_old=mu_new
    sigma_old=sigma_new
    covar=lambda_old*sigma_old
    jump = mrandomn(seed,covar,status=status)
    while status gt 0 do begin
      evals=eigenql(covar,EIGENVECTORS = evecs)
      evals[where(evals le 1.0d-6)]=1.0d-6
      covar=evecs#diag_matrix(evals)#invert(evecs)
      covar=(covar+transpose(covar))/2
      jump=mrandomn(seed,covar,status=status)
    endwhile
    chain_new = chain_old
    chain_new[variable_param] = chain_old[variable_param] + jump
    coeff_new=!null
    for i=0,(nsteps-1) do begin                   
      coeff_buffer = chain_new[i] & while (coeff_buffer lt 0) do coeff_buffer = chain_old[i] + (mrandomn(seed,covar))[i]
      coeff_new = [coeff_new,coeff_buffer]
    endfor 
    par1_new = chain_new[nsteps]   & while ((par1_new lt min(ran1)) or (par1_new gt max(ran1))) do par1_new = chain_old[nsteps]   + (mrandomn(seed,covar))[nsteps]
    par2_new = chain_new[nsteps+1] & while ((par2_new lt min(ran2)) or (par2_new gt max(ran2))) do par2_new = chain_old[nsteps+1] + (mrandomn(seed,covar))[nsteps+random2]
    par3_new = chain_new[nsteps+2] & while ((par3_new lt min(ran3)) or (par3_new gt max(ran3))) do par3_new = chain_old[nsteps+2] + (mrandomn(seed,covar))[nsteps+random3]
    par4_new = chain_new[nsteps+3] & while ((par4_new lt min(ran4)) or (par4_new gt max(ran4))) do par4_new = chain_old[nsteps+3] + (mrandomn(seed,covar))[nsteps+random4]
    par5_new = chain_new[nsteps+4] & while ((par5_new lt min(ran5)) or (par5_new gt max(ran5))) do par5_new = chain_old[nsteps+4] + (mrandomn(seed,covar))[nsteps+random5]
    if keyword_set(Tuffs_attenuation) then begin
      grid = {tauB_f_grid:par1_new, rdisk_grid:par2_new, F_grid:par3_new, cosi_grid:par4_new, rbulge_grid:par5_new}
      models=lightning_models_Keith(steps=steps,grid=grid,Filters=Filters,/Tuffs,coeff=folder_Tuffs_attenuation)
    endif else begin
      grid = {tauV_DIFF_grid:par1_new, delta_grid:par2_new, tauV_BC_grid:par3_new}
      models=lightning_models_Keith(steps=steps,grid=grid,Filters=Filters)                                            
    endelse
    chain_new = [coeff_new,par1_new,par2_new,par3_new,par4_new,par5_new]          
    Lmod_new = total(models * (rebin(reform(coeff_new,1,n_elements(coeff_new)),Nfilters+1,n_elements(coeff_new))),2) 
    chi2_new = total(((Lmod_new - Lmod_obs)/Lmod_unc)^2,/NaN)       
    ratio = exp(-(chi2_new-chi2_old)/2.d0)
    if ratio gt randomu(seed,1) then begin  $                   
      chain = [[chain],[chain_new]]    &$                        
      chi2_chain = [chi2_chain,chi2_new]   &$    
      chain_old = chain_new       &$            
      chi2_old = chi2_new
    endif
    if keyword_set(adaptive) then begin
      lambda_new = 10.d0^(alog10(lambda_old)+gamma_new*(min([1,ratio])-alpha_star))
      mu_new=mu_old+gamma_new*(chain_old-mu_old)
      sigma_new=sigma_old+gamma_new*((chain_old[variable_param]-mu_old[variable_param])#$
            (chain_old[variable_param]-mu_old[variable_param])-sigma_old)
    endif
  endfor

  Nchain = (chain.dim)[1]
  steps_Mstar = steps.Mstar
  
  if ~keyword_set(Tuffs_attenuation) then chain=chain[0:(nsteps+2),*]
  
  if ~keyword_set(NGC) then begin
    save,chain,chi2_chain,galaxy_id,sigma_new,lambda_new,mu_new,steps_Mstar,steps_bounds, $
       filename=outfolder_plots+'ID_'+galaxy_id+'_chain.idl'
    endif else begin
    save,chain,chi2_chain,galaxy_id,sigma_new,lambda_new,mu_new,steps_Mstar,steps_bounds, $
       filename=outfolder_plots+'NGC_'+galaxy_id+'_chain.idl'
  endelse     
     
  if keyword_set(print_time) then begin
    t1 = systime(/sec)
    dt = t1-t0
    nT=string(ceil(alog10(dt)+2) > 6)
    nC=string(ceil(alog10(nchain+1)))
    
    print,''
      print,'Done in '+string(dt,form='(F'+nT+'.1)')+' seconds'
      print,'Number of accepted trials: '+string(nchain,form='(I'+nC+')')
      print,'Percentage of total trials: '+string(nchain/ntrials*100.d0,form='(F5.1)')+' %'
      print,''
  endif

end

