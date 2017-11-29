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
                       nonebular=nonebular,Zmet=Zmet,IMF=IMF,Filters=Filters

; this function generates the spectra, SEDs, and stellar parameters, given a set of filters and steps boundaries

if (n_elements(filter_labels) eq 0) then $
  filter_labels=['GALEX_FUV','UVOT_W2','UVOT_M2','GALEX_NUV','UVOT_W1','SDSS_u','SDSS_g','SDSS_r','SDSS_i','SDSS_z',$
                 '2MASS_J','2MASS_H','2MASS_Ks','W1','SPITZER_I1','SPITZER_I2','W2']
if (n_elements(steps_bounds) eq 0) then steps_bounds = [0d0,1.d7,1.d8,1.d9,5.d9,13.8d9]
if (n_elements(z_shift) eq 0) then z_shift = 0.0
if (n_elements(dtime_SF) eq 0) then dtime_SF = 5.d5
if (n_elements(nolines) eq 0) then nolines = 0
if (n_elements(nonebular) eq 0) then nebtag='nebular_' else nebtag=''
if (n_elements(Zmet) eq 0) then Ztag='0.020' else Ztag=string(Zmet,form='(F5.3)')
if (n_elements(IMF) eq 0) then IMFtag='Kroupa01' else IMFtag=IMF

Nsteps = N_elements(steps_bounds) - 1

Dropbox = '~/Dropbox/'
bursts_folder = Dropbox+'PEGASE2/PEGASE2/01-models/04-single_burst/'
filename = bursts_folder+IMFtag+'/'+IMFtag+'_Z'+Ztag+'_'+nebtag+'spec.idl'

sObj = OBJ_NEW('IDL_Savefile',filename)
;sObj = OBJ_NEW('IDL_Savefile',Dropbox+'PEGASE2/PEGASE2/01-models/04-single_burst/Kroupa01/Kroupa01_Z0.020_nebular_spec.idl')
;sObj = OBJ_NEW('IDL_Savefile',Dropbox+'PEGASE2/PEGASE2/01-models/04-single_burst/Kroupa01/Kroupa01_Z0.020_spec.idl')
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
Lnu_burst      = Lnu_burst_rest / (1.d0+z_shift) ; observed Lnu spectrum, where Lnu = 4*!dpi*DL^2*Fnu
;Lbol_burst = Lbol/(1.d0+z_shift) & Lbol   = !null
;Q0_burst   = Nlyc/(1.d0+z_shift) & Nlyc   = !null
;Lnu_burst   = Lnu                 & Lnu    = !null
Lbol_burst  = Lbol                & Lbol   = !null ; restframe bolometric luminosity 
Q0_burst    = Nlyc                & Nlyc   = !null
Mstar_burst = Mstars              & Mstars = !null
GET_FILTERS,filter_labels,nu_obs,Filters,Lnu=[1.d0] # [0*(nu_obs) + 3631],$
            filters_dir=Dropbox+'Filters/',mean_wave=mean_wave,$
            /plot_filters;,mean_nu=mean_nu,mean_Lnu=mean_Lnu,mean_Lnu=mean_Lnu,sigma_wave=sigma_wave

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
            filters_dir=Dropbox+'Filters/',mean_Lnu=mean_Lnu;,$
            ;/plot_filters;,mean_wave=mean_wave,mean_nu=mean_nu,sigma_wave=sigma_wave

mean_Lnu_steps = transpose(temporary(mean_Lnu))

steps = {filter_labels : filter_labels, $
         wave_filters  : wave_filters,  $
         mean_Lnu      : mean_Lnu_steps,$
         wave_rest     : wave_rest,     $ ; restframe wavelength
         wave_obs      : wave_obs,      $
         Lnu           : Lnu_steps,     $ ; observed Lnu spectrum, where Lnu = 4*!dpi*DL^2*Fnu = Lnu_rest / (1+z) 
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
                   maxpost=maxpost,outfolder=outfolder,Regions=regions,noprior=noprior

tauV_DIFF_grid = grid.tauV_DIFF_grid & n1 = N_ELEMENTS(tauV_DIFF_grid)
delta_grid     = grid.delta_grid     & n2 = N_ELEMENTS(delta_grid)
tauV_BC_grid   = grid.tauV_BC_grid   & n3 = N_ELEMENTS(tauV_BC_grid)

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

if N_ELEMENTS(fit.CHISQR_GRID.dim) lt 5 then Nsim=1 else Nsim = (fit.CHISQR_GRID.dim)[4]

nu_obs = 1.d4*!cv.clight/steps.wave_filters

if (N_ELEMENTS(noprior) eq 0) then begin
  P = DOUBLE((tauV_DIFF_grid # REPLICATE(1.d0,n3)) le (2.0*(REPLICATE(1.d0,n1) # tauV_BC_grid)))
  ;P = double((((tauV_DIFF_grid # (1.d0/tauV_BC_grid)) lt 2.0) and $
  ;             (tauV_DIFF_grid # (1.d0/tauV_BC_grid)) gt 0.5))
endif else begin
  P = replicate(1.d0,n1) # replicate(1.d0,n3)
endelse
    
P3D = dblarr(n1,n2,n3,/nozero)
  for k2=0,(n2-1) do P3D[*,k2,*] = P
  P3D = P3D/total(P3D)
P4D = dblarr(n1,n2,n3,1,nsim,/nozero)
  for ss=0,(nsim-1) do P4D[*,*,*,0,ss] = P3D

Post = dblarr(n1,n2,n3,/nozero)

minchisqr           = dblarr(Nfit,Nsim)
minchisqr_loc       = lonarr(Nfit,Nsim)
minchisqr_ind       = lonarr(3,Nfit,Nsim)
coeff_minchisqr     = dblarr(Nsteps,Nfit,Nsim)
minchisqr_Lnu       = dblarr(Nfilters,Nfit,Nsim)
minchisqr_Lnu_unred = dblarr(Nfilters,Nfit,Nsim)
minchisqr_LIR       = dblarr(Nfit,Nsim)

maxpost_ind     = lonarr(3,Nfit,Nsim)
chisqr_maxpost  = dblarr(Nfit,Nsim)
coeff_maxpost   = dblarr(Nsteps,Nfit,Nsim)
Lnu_mod         = dblarr(Nfilters,Nfit,Nsim)
Lnu_mod_unred   = dblarr(Nfilters,Nfit,Nsim)
LIR_mod         = dblarr(Nfit,Nsim)

rep_filt = replicate(1.d0,Nfilters)

for ss=0,(Nsim-1) do begin
  for h=0,(Nfit-1) do begin
    minchisqr[h,ss] = min(fit.chisqr_grid[*,*,*,h,ss],minchi2_loc)
    minchisqr_loc[h,ss] = minchi2_loc
  endfor
endfor

reg_form = '(I0'+strtrim(ceil(alog10(Nreg)),2)+')'
mydevice = !d.name
set_plot,'ps'
for h=0,(Nfit-1) do begin
    ;min(fit.chisqr_grid[*,*,*,0,*],dim=5,loc)  
    coeff_reg = fit.coeff_grid[*,*,*,*,h,*]
    Lmod_reg = fit.Lmod_grid[*,*,*,*,h,*]
    chisqr_reg = fit.chisqr_grid[*,*,*,h,*]
    minchisqr_reg = 0*chisqr_reg
    for ss=0,(Nsim-1) do minchisqr_reg[*,*,*,0,ss] = minchisqr[h,ss]
    L4D = exp(-(chisqr_reg-minchisqr_reg)/2.d0) ;& L4D /= (total(total(total(L,1),1),1))[0]
    Post4D = P4D * L4D ;                             & Post4D /= (total(total(total(Post4d,1),1),1))[0]
  for ss=0,(Nsim-1) do begin
    L = L4D[*,*,*,0,ss]       & L    /= (total(L))[0]
    Post = Post4D[*,*,*,0,ss] & Post /= (total(Post))[0]

    minchi2_ind = lonarr(3,/nozero)
    minchi2_loc = (minchisqr_loc[h,ss])[0]
    if (n3 eq 1) then begin
      minchi2_ind[2] = 0
      if n2 eq 1 then begin
        minchi2_ind[1] = 0
        minchi2_ind[0] = minchi2_loc
      endif else begin
        minchi2_ind[0:1] = array_indices(dblarr(n1,n2,/nozero),minchi2_loc)
      endelse
    endif else if n2 eq 1 then begin
      minchi2_ind[1] = 0
      minchi2_ind[[0,2]] = array_indices(dblarr(n1,n3,/nozero),minchi2_loc)
    endif else begin
      minchi2_ind = array_indices(dblarr(n1,n2,n3,/nozero),minchi2_loc)
    endelse
    minchisqr_ind[*,h,ss] = minchi2_ind
    coeff = coeff_reg[*,minchi2_ind[0],minchi2_ind[1],minchi2_ind[2],0,ss]
    coeff_minchisqr[*,h,ss] = coeff
    Lmod = Lmod_reg[*,minchi2_ind[0],minchi2_ind[1],minchi2_ind[2],0,ss]
    minchisqr_Lnu_unred[*,h,ss] = total((rep_filt # coeff)*steps.mean_Lnu,2)
    minchisqr_Lnu[*,h,ss] = Lmod[0:-2]
    minchisqr_LIR[h,ss] = Lmod[-1]

    maxpost = max(Post,loc)
    ind = array_indices(dblarr(n1,n2,n3,/nozero),loc)
    ind = lonarr(3,/nozero)
    if n3 eq 1 then begin
      ind[2] = 0
      if n2 eq 1 then begin
        ind[1] = 0
        ind[0] = loc
      endif else begin
        ind[0:1] = array_indices(dblarr(n1,n2,/nozero),loc)
      endelse
    endif else if n2 eq 1 then begin
      ind[1] = 0
      ind[[0,2]] = array_indices(dblarr(n1,n3,/nozero),loc)
    endif else begin
      ind = array_indices(dblarr(n1,n2,n3,/nozero),loc)
    endelse
    maxpost_ind[*,h,ss] = ind
    chisqr_maxpost[h,ss] = ((chisqr_reg)[*,*,*,0,ss])[loc]
    coeff = coeff_reg[*,ind[0],ind[1],ind[2],0,ss]
    coeff_maxpost[*,h,ss] = coeff
    Lnu_mod_unred[*,h,ss] = total((rep_filt # coeff)*steps.mean_Lnu,2)
    Lmod = Lmod_reg[*,ind[0],ind[1],ind[2],0,ss]
    Lnu_mod[*,h,ss] = Lmod[0:-2]
    LIR_mod[h,ss] = Lmod[-1]
    
    if ss eq 0 then begin
      Nh = N0+h
      print,Nh,minchisqr[h,ss],chisqr_maxpost[h,ss],coeff_maxpost[*,h,ss],tauV_DIFF_grid[ind[0]],$
            delta_grid[ind[1]],tauV_BC_grid[ind[2]],form='(I,2F10.2,'+strtrim(Nsteps,2)+'E10.2,3F12.3)'

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
      ;spawn,'open '+fileout+' &'
    
      if n1 gt 1 then begin     ; n1 > 1
        if n3 gt 1 then begin   ; n1 > 1,         n3 > 1
          Post1 = total(Post,3) & L1 = total(L,3) & P1 = total(P3D,3)
          ;if n2 eq 1 idl collapses second dimension after this sum
          if n2 gt 1 then begin ; n1 > 1, n2 > 1, n3 > 1
            Post1 = total(Post1,2) & L1 = total(L1,2) & P1 = total(P1,2)
          endif
        endif else begin   ; n1 > 1,         n3 = 1
          if n2 gt 1 then begin ; n1 > 1, n2 > 1, n3 = 1
            Post1 = total(Post,2) & L1 = total(L,2) & P1 = total(P3D,2)
          endif else begin      ; n1 > 1, n2 = 1, n3 = 1
            Post1 = Post & L1 = L & P1 = P3D
          endelse
        endelse
        fileout=outfolder+prefix+'_par1.eps'
        device,filename=fileout,/encaps,bits=8,/color
          plot,minmax(tauV_DIFF_grid),minmax([Post1,L1,P1]),$
               xtitle='Diffuse attenuation, tau_V_DIFF',ytitle='Probalibity Densities',/nodata
            oplot,tauV_DIFF_grid,Post1,color=54
            oplot,[1,1]*(tauV_DIFF_grid[ind[0]])[0],!y.crange,color=54
            oplot,tauV_DIFF_grid,L1,linestyle=2,color=254
            oplot,[1,1]*(tauV_DIFF_grid[minchi2_ind[0]])[0],!y.crange,linestyle=2,color=254
            oplot,tauV_DIFF_grid,P1,linestyle=1,color=54
        device,/close
        ;spawn,'open '+fileout+' &'
      endif
    
      if n2 gt 1 then begin
        if n3 gt 1 then begin   ;         n2 > 1, n3 > 1
          Post2 = total(total(Post,3),1)
          L2    = total(total(L,3),1)
          P2    = total(total(P3D,3),1)
        endif else begin        ; n1 > 1, n2 > 1, n3 = 1
          Post2 = total(Post,1)
          L2    = total(L,1)
          P2    = total(P3D,1)
        endelse
        fileout=outfolder+prefix+'_par2.eps'
        device,filename=fileout,/encaps,bits=8,/color
          plot,minmax(delta_grid),minmax([Post2,L2,P2]),$
               xtitle='Slope, delta',ytitle='Probalibity Densities',/nodata
            oplot,delta_grid,Post2,color=54
            oplot,[1,1]*(delta_grid[ind[1]])[0],!y.crange,color=54
            oplot,delta_grid,L2,linestyle=2,color=254
            oplot,[1,1]*(delta_grid[minchi2_ind[1]])[0],!y.crange,linestyle=2,color=254
            oplot,delta_grid,P2,linestyle=1,color=54
        device,/close
        ;spawn,'open '+fileout+' &'
      endif
      
      if n3 gt 1 then begin
        Post3 = total(total(Post,2),1)
        L3    = total(total(L,2),1)
        P3    = total(total(P3D,2),1)
        fileout=outfolder+prefix+'_par3.eps'
        device,filename=fileout,/encaps,bits=8,/color
          plot,minmax(tauV_BC_grid),minmax([Post3,L3,P3]),$
               xtitle='Birth could attenuation, tau_V_BC',ytitle='Probalibity Densities',/nodata
            oplot,tauV_BC_grid,Post3,color=54
            oplot,[1,1]*(tauV_BC_grid[ind[2]])[0],!y.crange,color=54
            oplot,tauV_BC_grid,L3,linestyle=2,color=254
            oplot,[1,1]*(tauV_BC_grid[minchi2_ind[2]])[0],!y.crange,linestyle=2,color=254
            oplot,tauV_BC_grid,P3,linestyle=1,color=54
        device,/close
        ;spawn,'open '+fileout+' &'
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


function lightning,Lnu_obs,Lnu_unc,LIR_obs,LIR_unc,Nsim=Nsim,steps=steps,grid=grid,models=models,$
                   float=float,regions=regions,Lnu_sim=Lnu_sim,LIR_sim=LIR_sim

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
  n1 = 41
  n2 = 41
  n3 = 41
  ran1 = [   0d,4d0] & tauV_DIFF_grid = ran1[0] + dindgen(n1)*(ran1[1]-ran1[0])/(n1-1) ;as in Leja et al. 2016
  ran2 = [-2.2d,.4d] & delta_grid     = ran2[0] + dindgen(n2)*(ran2[1]-ran2[0])/(n2-1) ;as in Leja et al. 2016
  ran3 = [   0d,4d0] & tauV_BC_grid   = ran3[0] + dindgen(n3)*(ran3[1]-ran3[0])/(n3-1) ;as in Leja et al. 2016
endif else begin
  tauV_DIFF_grid = grid.tauV_DIFF_grid & n1 = n_elements(tauV_DIFF_grid)
  delta_grid     = grid.delta_grid     & n2 = n_elements(delta_grid)
  tauV_BC_grid   = grid.tauV_BC_grid   & n3 = n_elements(tauV_BC_grid)
endelse

if n_elements(models) eq 0 then begin
  ttt = systime(/seconds)
  steps_LIR_grid = dblarr(Nsteps,n1,n2,n3)
  steps_mean_Lnu_red_grid = dblarr(Nfilters+1,Nsteps,n1,n2,n3)
  transmission = transpose(Filters)
  
  st=0
  steps_Lnu_red = matrix_unred(steps.Lnu[*,st],wave_rest,tauV_DIFF_arr=-1.d0*tauV_DIFF_grid,$
                               delta_arr=delta_grid,tauV_BC_arr=-1.d0*tauV_BC_grid) ; observed attenuated Lnu 
  Lbol = (steps.Lbol[st])[0]
  for k3=0,(n3-1) do begin
    for k2=0,(n2-1) do begin
      for k1=0,(n1-1) do begin
        for k=0,(Nfilters-1) do begin
          steps_mean_Lnu_red_grid[k,st,k1,k2,k3] = $
            trapez(steps_Lnu_red[*,k1,k2,k3]*transmission[*,k],nu_obs)
          ;integral will be positive because transmission is negative and nu is inverted 
          ;if steps_mean_Lnu_red_grid[k,st,k1,k2,k3] lt 0 then stop
        endfor
        steps_mean_Lnu_red_grid[Nfilters,st,k1,k2,k3] = steps.Lbol[st] + $
                                                        (1+z_shift)*trapez(steps_Lnu_red[*,k1,k2,k3],nu_rest)
        ;integral is negative because nu is inverted, so negative sign becames positive
        ;if trapez(steps_Lnu_red[*,k1,k2,k3],nu) gt 0 then stop
        ;IR Luminosity is measured in restframe units
      endfor  
    endfor
  endfor
  
  for st=1,(Nsteps-1) do begin
    ;steps_Lnu_red = dblarr(Nwave,n1,n2)
    steps_Lnu_red = matrix_unred(steps.Lnu[*,st],wave_rest,tauV_DIFF_arr=-1.d0*tauV_DIFF_grid,$
                                 delta_arr=delta_grid,tauV_BC_arr= 0.d0)
    for k2=0,(n2-1) do begin
      for k1=0,(n1-1) do begin
        for k=0,(Nfilters-1) do begin
          steps_mean_Lnu_red_grid[k,st,k1,k2,*] = $
            trapez(steps_Lnu_red[*,k1,k2]*transmission[*,k],nu_obs)
          ;integral will be positive because transmission is negative and nu is inverted 
          ;if steps_mean_Lnu_red_grid[k,st,k1,k2] lt 0 then stop
        endfor  
        steps_mean_Lnu_red_grid[Nfilters,st,k1,k2,*] = steps.Lbol[st] + $
                                                       (1+z_shift)*trapez(steps_Lnu_red[*,k1,k2],nu_rest)
        ;integral is negative because nu is inverted, so negative sign becames positive
        ;if trapez(steps_Lnu_red[*,k1,k2],nu) gt 0 then stop
      endfor
    endfor
  endfor
  steps_Lnu_red = !null
  tttt = systime(/seconds)
  print,'Time to generate models : ',tttt-ttt,' seconds'
endif else begin
  steps_mean_Lnu_red_grid = models
endelse

if n_elements(float) eq 0 then float=0
if (float eq 1) then begin
  coeff_grid = fltarr(Nsteps,n1,n2,n3,Nfit,Nsim,/nozero)    & replicate_inplace,coeff_grid,1.e0
  chisqr_grid = fltarr(n1,n2,n3,Nfit,Nsim,/nozero)          & replicate_inplace,chisqr_grid,1.e0
  Lmod_grid = fltarr(Nfilters+1,n1,n2,n3,Nfit,Nsim,/nozero) & replicate_inplace,Lmod_grid,0.e0
  coeff_temp = fltarr(Nfilters+1,Nsteps,n1,n2,n3,/nozero)   & replicate_inplace,coeff_temp,1.e0
  Lfit = fltarr(Nfilters+1,n1,n2,n3,/nozero)                & replicate_inplace,Lfit,1.e0
  Lsig = fltarr(Nfilters+1,n1,n2,n3,/nozero)                & replicate_inplace,Lfit,1.e0
endif else begin
  coeff_grid = dblarr(Nsteps,n1,n2,n3,Nfit,Nsim,/nozero)    & replicate_inplace,coeff_grid,0.d0
  chisqr_grid = dblarr(n1,n2,n3,Nfit,Nsim,/nozero)          & replicate_inplace,chisqr_grid,1.d0
  Lmod_grid = dblarr(Nfilters+1,n1,n2,n3,Nfit,Nsim,/nozero) & replicate_inplace,Lmod_grid,0.d0
  coeff_temp = dblarr(Nfilters+1,Nsteps,n1,n2,n3,/nozero)   & replicate_inplace,coeff_temp,1.d0
  Lfit = dblarr(Nfilters+1,n1,n2,n3,/nozero)                & replicate_inplace,Lfit,1.d0
  Lsig = dblarr(Nfilters+1,n1,n2,n3,/nozero)                & replicate_inplace,Lfit,1.d0
endelse

t1=systime(/sec)
for sr=0,(Nfit-1) do begin
  Lreg_tmp = Lsim[*,sr,*]
  Lsig_tmp = Lunc[*,sr]
  for ss=0,(Nsim-1) do begin
    Lsim_tmp = Lreg_tmp[*,0,ss]
    for k=0,Nfilters do begin
      Lfit[k,*,*,*] = Lsim_tmp[k]
      Lsig[k,*,*,*] = Lsig_tmp[k]
    endfor  
    for k3=0,(n3-1) do begin
      ;tauV_BC = tauV_BC_grid[k3]
      for k2=0,(n2-1) do begin
        ;delta = delta_grid[k2]
        for k1=0,(n1-1) do begin
          ;tauV_DIFF = tauV_DIFF_grid[k1]
          coeff_grid[*,k1,k2,k3,sr,ss] = $
            steps_derivative_chisqr_minimize(Lsim_tmp,Lsig_tmp,steps_mean_Lnu_red_grid[*,*,k1,k2,k3])
          ;status_grid[k1,k2,k3] = temporary(status_temp)
        endfor
      endfor
    endfor
    coeff = coeff_grid[*,*,*,*,sr,ss]
    for k=0,Nfilters do coeff_temp[k,*,*,*,*] = coeff
    Lmod_grid[*,*,*,*,sr,ss] = total(coeff_temp * steps_mean_Lnu_red_grid,2)
    chisqr_grid[*,*,*,sr,ss] = total(((Lmod_grid[*,*,*,*,sr,ss] - Lfit)/Lsig)^2,1,/NaN)
  endfor
  ti=systime(/seconds)
                                             nR=string(ceil(alog10(Nreg+1)))
  dt = ti-t1                               & nT=string(ceil(alog10(dt)+2) > 6)
  minchisqr = min(chisqr_grid[*,*,*,sr,0]) & nC=string(ceil(alog10(minchisqr)+3) > 7)
  print,string((sr+N0+1.)*100./Nreg,form='(F5.1)')+'% completed --- Done with region # '+$
        string(sr+N0+1,form='(I'+nR+')')+' out of '+strtrim(Nreg,2)+' @ '+$
        string(dt,form='(F'+nT+'.1)')+' seconds --- Min chisqr: '+string(minchisqr,form='(F'+nC+'.2)')
  
;  help,min(chisqr_grid[*,*,*,sr,0],loc)
;  ind = array_indices(chisqr_grid[*,*,*,sr,0],loc)
;  st=2
;  steps_Lnu_red = matrix_unred(steps.Lnu[*,st],wave_rest,tauV_DIFF_arr=-1.d0*tauV_DIFF_grid,$
;                               delta_arr=delta_grid,tauV_BC_arr= 0.d0)
;  help,trapez(steps_Lnu_red[*,ind[0],ind[1]],nu)
;  coeff = coeff_grid[2,ind[0],ind[1],ind[2],sr,0]
;  print,Lmod_grid[12,ind[0],ind[1],ind[2],0,0],steps.Lbol[2]*coeff-trapez(steps_Lnu_red[*,ind[0],ind[1]],nu)*coeff
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
