.compile ~/Dropbox/Libraries/GitHub/lightning/lightning.pro
@~/Dropbox/Libraries/GitHub/lightning/startup_RTE.pro

device,decompose=0
loadct,39

filter_labels=['GALEX_FUV',$
               'GALEX_NUV',$
               'SDSS_u',$;'KP_B',$
               'SDSS_g',$;'KP_V',$
               'SDSS_r',$;'KP_R',$
               'SDSS_i',$;'KP_I',$
               'SDSS_z',$
               '2MASS_J',$
               '2MASS_H',$
               '2MASS_Ks',$
               'W1',$
               'SPITZER_I1',$
               'SPITZER_I2',$
               'W2']

folder = './Libraries/GitHub/lightning/'

steps_bounds = [0d0,1.d7,1.d8,1.d9,5.d9,13.6d9] & Nsteps = n_elements(steps_bounds) - 1
steps = STEPS_STELLAR(filter_labels=filter_labels,steps_bounds=steps_bounds,z_shift=0.0,dtime_SF=5.d5,folder=folder)

n1 = 21 ; run with 81 to save
n2 = 27
n3 = 21 ; run with 81 to save

ran1=[   0d,4.0d] & tauV_DIFF_grid=ran1[0]+dindgen(n1)*(ran1[1]-ran1[0])/(n1-1) & if n1 eq 1 then tauV_DIFF_grid=[0.d0]
ran2=[-2.2d,0.4d] & delta_grid    =ran2[0]+dindgen(n2)*(ran2[1]-ran2[0])/(n2-1) & if n2 eq 1 then delta_grid=[0.d0]
ran3=[   0d,4.0d] & tauV_BC_grid  =ran3[0]+dindgen(n3)*(ran3[1]-ran3[0])/(n3-1) & if n3 eq 1 then tauV_BC_grid=[0.d0]

grid = {tauV_DIFF_grid:tauV_DIFF_grid, delta_grid:delta_grid, tauV_BC_grid:tauV_BC_grid}

Lnu_obs = total(steps.mean_Lnu,2)   &   Lnu_unc = 0.05*Lnu_obs
LIR_obs = 0                         &   LIR_unc = 0

Nsim=1

models=!null
fit = lightning(Lnu_obs,Lnu_unc,LIR_obs,LIR_unc,steps=steps,grid=grid,models=models,Nsim=Nsim)


;Mstar_bestfit_sim = total(reform(bestfit.coeff) * (steps.mstar # replicate(1.d0,Nsim)),1)
;Mstar_maxpost_sim = total(reform(maxpost.coeff) * (steps.mstar # replicate(1.d0,Nsim)),1)
;
;print,percentile(Mstar_bestfit_sim,[0.16,0.5,0.84])
;print,percentile(Mstar_maxpost_sim,[0.16,0.5,0.84])




Nfilters = filter_labels.length

coeff_tru = [1.d,1.d,10.d,10.d,1.d]
par1_tru = 0.4d  &  k1_tru = interpol(lindgen(n1),tauV_DIFF_grid,par1_tru)
par2_tru = 0.0d  &  k2_tru = interpol(lindgen(n2),delta_grid,    par2_tru)
par3_tru = 0.4d  &  k3_tru = interpol(lindgen(n3),tauV_BC_grid,  par3_tru)
Lmod_tru = total(models[*,*,k1_tru,k2_tru,k3_tru] * (replicate(1.d0,Nfilters+1)#coeff_tru),2)
Lmod_unc = 0.05*Lmod_tru
Lmod_obs = Lmod_tru + Lmod_unc*randomn(seed,Nfilters+1)

Lnu_obs = Lmod_obs[0:-2] & Lnu_unc = Lmod_unc[0:-2]
LIR_obs = Lmod_obs[-1]   & LIR_unc = Lmod_unc[-1]

Nsim=101
fit = lightning(Lnu_obs,Lnu_unc,LIR_obs,LIR_unc,steps=steps,grid=grid,models=models,Nsim=Nsim)
outfolder_plots=folder+'plots/'
plot_lightning,Lnu_obs,Lnu_unc,LIR_obs,LIR_unc,steps=steps,grid=grid,$
               fit=fit,bestfit=bestfit,maxpost=maxpost,outfolder=outfolder_plots




;coeff_guess = [2.0d,2.0d,2.0d,2.0d,2.0d] & logpsi_scale = 0.05d0*coeff_guess/coeff_guess
coeff_guess = coeff_tru & logpsi_scale = 0.1d0*coeff_guess/coeff_guess
;par1_guess = 1.0d  &  par1_scale = 0.01d  &  k1_guess = interpol(lindgen(n1),tauV_DIFF_grid,par1_guess)
;par2_guess = -.4d  &  par2_scale = 0.01d  &  k2_guess = interpol(lindgen(n2),delta_grid,    par2_guess)
;par3_guess = 1.0d  &  par3_scale = 0.01d  &  k3_guess = interpol(lindgen(n3),tauV_BC_grid,  par3_guess)
par1_guess = 0.4d  &  par1_scale = 0.01d  &  k1_guess = interpol(lindgen(n1),tauV_DIFF_grid,par1_guess)
par2_guess = 0.0d  &  par2_scale = 0.01d  &  k2_guess = interpol(lindgen(n2),delta_grid,    par2_guess)
par3_guess = 0.4d  &  par3_scale = 0.01d  &  k3_guess = interpol(lindgen(n3),tauV_BC_grid,  par3_guess)


coeff_old = coeff_guess
par1_old = par1_guess
par2_old = par2_guess
par3_old = par3_guess
;models[k,st,k1,k2,k3]
Lmod_old = total(models[*,*,k1_guess,k2_guess,k3_guess] * (replicate(1.d0,Nfilters+1)#coeff_old),2)
chi2_old = total(((Lmod_old - Lmod_obs)/Lmod_unc)^2,/NaN)

chain_old = [alog10(coeff_old),par1_old,par2_old,par3_old]
scale = [logpsi_scale,par1_scale,par2_scale,par3_scale]
chain = [chain_old]
chi2_chain = [chi2_old]

t0=systime(/sec)
Npar = chain_old.length
Ntrials = 1e5
for trial=0,Ntrials-1 do begin                                $
  jump = randomn(seed,Npar)                                  &$
  chain_new = chain_old + scale*jump                         &$
  logcoeff_new = chain_new[0:4]                              &$
  ;logpsi1_new = chain_new[0]                                 &$
  ;logpsi2_new = chain_new[1]                                 &$
  ;logpsi3_new = chain_new[2]                                 &$
  ;logpsi4_new = chain_new[3]                                 &$
  ;logpsi5_new = chain_new[4]                                 &$
  ;logcoeff_new = [logpsi1_new,logpsi2_new,logpsi3_new,logpsi4_new,logpsi5_new] &$
  par1_new = chain_new[5] & while ((par1_new lt min(tauV_DIFF_grid)) or (par1_new gt max(tauV_DIFF_grid))) do par1_new = chain_old[5] + scale[5]*randomn(seed)  &$
  par2_new = chain_new[6] & while ((par2_new lt min(delta_grid))     or (par2_new gt max(delta_grid)))     do par2_new = chain_old[6] + scale[6]*randomn(seed)  &$
  par3_new = chain_new[7] & while ((par3_new lt min(tauV_BC_grid))   or (par3_new gt max(tauV_BC_grid)))   do par3_new = chain_old[7] + scale[7]*randomn(seed)  &$
  k1_new = interpol(lindgen(n1),tauV_DIFF_grid,par1_new)     &$
  k2_new = interpol(lindgen(n2),delta_grid,    par2_new)     &$
  k3_new = interpol(lindgen(n3),tauV_BC_grid,  par3_new)     &$
  chain_new = [logcoeff_new,par1_new,par2_new,par3_new]      &$
  Lmod_new = total(models[*,*,k1_new,k2_new,k3_new] * (replicate(1.d0,Nfilters+1)#(1.d1^logcoeff_new)),2) &$
  chi2_new = total(((Lmod_new - Lmod_obs)/Lmod_unc)^2)       &$
  ratio = exp(-(chi2_new-chi2_old)/2.d0)                     &$ ;= likenew/likeold
  if ratio gt randomu(seed,1) then begin                      $
    chain = [[chain],[chain_new]]                            &$
    chi2_chain = [chi2_chain,chi2_new]                       &$
    chain_old = chain_new                                    &$
    chi2_old = chi2_new
t1 = systime(/sec)
dt = t1-t0
print,'####################'
print,dt,'seconds'
print,'####################'

Nchain = (chain.dim)[1] 
print, Nchain

Mstar_chain = total((1.d1^chain[0:4,*]) * (steps.Mstar # replicate(1.d0,Nchain)),1)
SFR100_chain = 0.1*reform(1.d1^chain[0,*]) + 0.9*reform(1.d1^chain[1,*])
print,percentile(SFR100_chain,[0.025,0.16,0.50,0.84,0.975])
print,percentile(Mstar_chain, [0.025,0.16,0.50,0.84,0.975])

Mstar_tru = total(steps.mstar*coeff_tru)
SFR100_tru = 0.1*coeff_tru[0] + 0.9*coeff_tru[1]


erase
!p.multi=[64,8,8]    & plothist,chain[0,*],color=255
!p.multi=[64-8,8,8]  & plot,chain[0,*],chain[1,*]
!p.multi=[64-16,8,8] & plot,chain[0,*],chain[2,*]
!p.multi=[64-24,8,8] & plot,chain[0,*],chain[3,*]
!p.multi=[64-32,8,8] & plot,chain[0,*],chain[4,*]
!p.multi=[64-40,8,8] & plot,chain[0,*],chain[5,*]
!p.multi=[64-48,8,8] & plot,chain[0,*],chain[6,*]
!p.multi=[64-56,8,8] & plot,chain[0,*],chain[7,*]

!p.multi=[63-8,8,8]  & plothist,chain[1,*],color=255
!p.multi=[63-16,8,8] & plot,chain[1,*],chain[2,*]
!p.multi=[63-24,8,8] & plot,chain[1,*],chain[3,*]
!p.multi=[63-32,8,8] & plot,chain[1,*],chain[4,*]
!p.multi=[63-40,8,8] & plot,chain[1,*],chain[5,*]
!p.multi=[63-48,8,8] & plot,chain[1,*],chain[6,*]
!p.multi=[63-56,8,8] & plot,chain[1,*],chain[7,*]

!p.multi=[62-16,8,8] & plothist,chain[2,*],color=255
!p.multi=[62-24,8,8] & plot,chain[2,*],chain[3,*]
!p.multi=[62-32,8,8] & plot,chain[2,*],chain[4,*]
!p.multi=[62-40,8,8] & plot,chain[2,*],chain[5,*]
!p.multi=[62-48,8,8] & plot,chain[2,*],chain[6,*]
!p.multi=[62-56,8,8] & plot,chain[2,*],chain[7,*]

!p.multi=[61-24,8,8] & plothist,chain[3,*],color=255
!p.multi=[61-32,8,8] & plot,chain[3,*],chain[4,*] 
!p.multi=[61-40,8,8] & plot,chain[3,*],chain[5,*]
!p.multi=[61-48,8,8] & plot,chain[3,*],chain[6,*]
!p.multi=[61-56,8,8] & plot,chain[3,*],chain[7,*]

!p.multi=[60-32,8,8] & plothist,chain[4,*],color=255
!p.multi=[60-40,8,8] & plot,chain[4,*],chain[5,*]
!p.multi=[60-48,8,8] & plot,chain[4,*],chain[6,*]
!p.multi=[60-56,8,8] & plot,chain[4,*],chain[7,*]

!p.multi=[59-40,8,8] & plothist,chain[5,*],color=255
!p.multi=[59-48,8,8] & plot,chain[5,*],chain[6,*]
!p.multi=[59-56,8,8] & plot,chain[5,*],chain[7,*]

!p.multi=[58-48,8,8] & plothist,chain[6,*],color=255
!p.multi=[58-56,8,8] & plot,chain[6,*],chain[7,*]

!p.multi=[57-56,8,8] & plothist,chain[7,*],color=255



erase
!p.multi=[64,8,8]    & plothist,chain[0,*],color=255
!p.multi=[64-8,8,8]  & density,chain[0,*],chain[1,*]
!p.multi=[64-16,8,8] & density,chain[0,*],chain[2,*]
!p.multi=[64-24,8,8] & density,chain[0,*],chain[3,*]
!p.multi=[64-32,8,8] & density,chain[0,*],chain[4,*]
!p.multi=[64-40,8,8] & density,chain[0,*],chain[5,*]
!p.multi=[64-48,8,8] & density,chain[0,*],chain[6,*]
!p.multi=[64-56,8,8] & density,chain[0,*],chain[7,*]

!p.multi=[63-8,8,8]  & plothist,chain[1,*],color=255
!p.multi=[63-16,8,8] & density,chain[1,*],chain[2,*]
!p.multi=[63-24,8,8] & density,chain[1,*],chain[3,*]
!p.multi=[63-32,8,8] & density,chain[1,*],chain[4,*]
!p.multi=[63-40,8,8] & density,chain[1,*],chain[5,*]
!p.multi=[63-48,8,8] & density,chain[1,*],chain[6,*]
!p.multi=[63-56,8,8] & density,chain[1,*],chain[7,*]

!p.multi=[62-16,8,8] & plothist,chain[2,*],color=255
!p.multi=[62-24,8,8] & density,chain[2,*],chain[3,*]
!p.multi=[62-32,8,8] & density,chain[2,*],chain[4,*]
!p.multi=[62-40,8,8] & density,chain[2,*],chain[5,*]
!p.multi=[62-48,8,8] & density,chain[2,*],chain[6,*]
!p.multi=[62-56,8,8] & density,chain[2,*],chain[7,*]

!p.multi=[61-24,8,8] & plothist,chain[3,*],color=255
!p.multi=[61-32,8,8] & density,chain[3,*],chain[4,*] 
!p.multi=[61-40,8,8] & density,chain[3,*],chain[5,*]
!p.multi=[61-48,8,8] & density,chain[3,*],chain[6,*]
!p.multi=[61-56,8,8] & density,chain[3,*],chain[7,*]

!p.multi=[60-32,8,8] & plothist,chain[4,*],color=255
!p.multi=[60-40,8,8] & density,chain[4,*],chain[5,*]
!p.multi=[60-48,8,8] & density,chain[4,*],chain[6,*]
!p.multi=[60-56,8,8] & density,chain[4,*],chain[7,*]

!p.multi=[59-40,8,8] & plothist,chain[5,*],color=255
!p.multi=[59-48,8,8] & density,chain[5,*],chain[6,*]
!p.multi=[59-56,8,8] & density,chain[5,*],chain[7,*]

!p.multi=[58-48,8,8] & plothist,chain[6,*],color=255
!p.multi=[58-56,8,8] & density,chain[6,*],chain[7,*]

!p.multi=[57-56,8,8] & plothist,chain[7,*],color=255





dist = (chain - (median(chain,dim=2)#replicate(1.d0,nchain)))/(stddev(chain,dim=2)#replicate(1.d0,nchain))
dist *= dist


