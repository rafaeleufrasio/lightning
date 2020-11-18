pro MCMC_savefits_tuffs,galaxy,filter_labels,lnu_obs,lnu_unc,Chain_folder,outfolder,Nchain,Nsteps_max,file_name=file_name,$
                      L_star_abs_table_loc=L_star_abs_table_loc,prior=prior,_REF_EXTRA=_extra_savefits

Nfilters = (size(lnu_obs))[1]

x = {   Galaxy_ID			    : '',										$
        RAJ2000				    : !values.D_NaN,							$
        DECJ2000			    : !values.D_NaN,							$
        redshift			    : !values.D_NaN,							$ 
        metallicity			    : !values.D_NaN,							$ 
        N_J                     : !values.D_NaN,                            $
        N_J_ERR                 : !values.D_NaN,                            $
        Q_J                     : !values.D_NaN,                            $
        Q_J_ERR                 : !values.D_NaN,                            $
        filter_labels		    : strarr(Nfilters),			                $
        WAVE_FILTERS            : !values.D_NaN*dblarr(Nfilters),			$ 
        Lnu_obs				    : !values.D_NaN*dblarr(Nfilters),			$ 
        Lnu_unc				    : !values.D_NaN*dblarr(Nfilters),			$ 
        Lnu_mod                 : !values.D_NaN*dblarr(Nfilters,Nchain),    $
        Lnu_mod_unred           : !values.D_NaN*dblarr(Nfilters,Nchain),    $
        LTIR_mod                : !values.D_NaN*dblarr(Nchain),             $
        LFUV_mod                : !values.D_NaN*dblarr(Nchain),             $
        wave_hires_dustmod      : !values.D_NaN*dblarr(1001),               $
        wave_hires_starmod      : !values.D_NaN*dblarr(1221),               $
        wave_hires_totalmod     : !values.D_NaN*dblarr(1971),               $
        lnu_hires_dustmod       : !values.D_NaN*dblarr(1001),               $
        lnu_hires_starmod       : !values.D_NaN*dblarr(1221),               $
        lnu_hires_starmod_unred : !values.D_NaN*dblarr(1221),               $
        lnu_hires_totalmod      : !values.D_NaN*dblarr(1971),               $
        Afuv                    : !values.D_NaN*dblarr(Nchain),             $
        Ab                      : !values.D_NaN*dblarr(Nchain),             $
        Av                      : !values.D_NaN*dblarr(Nchain),             $
        Anir                    : !values.D_NaN*dblarr(Nchain),             $
        Nsteps				    : 0L,										$ 
        Steps_bounds		    : !values.D_NaN*dblarr(Nsteps_max+1),		$
        chisqr_lightning 	    : !values.D_NaN*dblarr(Nchain), 			$
        chisqr_prior            : !values.D_NaN*dblarr(Nchain), 			$
        SFH					    : !values.D_NaN*dblarr(Nsteps_max,Nchain),	$
        steps_Mstar_coeff       : !values.D_NaN*dblarr(Nsteps_max),        	$
        Mstar                   : !values.D_NaN*dblarr(Nchain),             $
        steps_Mstar             : !values.D_NaN*dblarr(Nsteps_max,Nchain),  $
        tauB_F 			        : !values.D_NaN*dblarr(Nchain), 			$
        rdisk				    : !values.D_NaN*dblarr(Nchain),				$
        F				        : !values.D_NaN*dblarr(Nchain),				$
        cosi			        : !values.D_NaN*dblarr(Nchain),				$
        rbulge			        : !values.D_NaN*dblarr(Nchain),				$
        alpha 			        : !values.D_NaN*dblarr(Nchain), 			$
        u_min				    : !values.D_NaN*dblarr(Nchain),				$
        u_max				    : !values.D_NaN*dblarr(Nchain),				$
        gamma			        : !values.D_NaN*dblarr(Nchain),				$
        q_pah			        : !values.D_NaN*dblarr(Nchain)				}


ngal=size(galaxy,/str)
ngal=ngal.n_elements 

out=replicate(x,ngal)



for i=0,(ngal-1) do begin

  radec,galaxy[i].raj2000,galaxy[i].decj2000,ihr,imin,xsec,ideg,imn,xsc
  if galaxy[i].decj2000 ge 0 then begin
    gal_id = 'J'+string(ihr,form='(I02)')+''+string(imin,form='(I02)')+''+string(xsec,form='(F05.2)')+ $
      '+'+string(ideg,form='(I02)')+''+string(imn,form='(I02)')+''+string(xsc,form='(F04.1)')
    field='N'
  endif else begin
    gal_id = 'J'+string(ihr,form='(I02)')+''+string(imin,form='(I02)')+''+string(xsec,form='(F05.2)')+ $
      ''+string(ideg,form='(I03)')+''+string(imn,form='(I02)')+''+string(xsc,form='(F04.1)')
    field='S'
  endelse

  restore,Chain_folder+gal_id+'_tuffs_chain.idl'
  Nsteps = steps_bounds.length - 1
  
  steps = STEPS_STELLAR(filter_labels=filter_labels,steps_bounds=steps_bounds,$
                        z_shift=galaxy[i].redshift,_EXTRA=_extra_savefits)
  steps_alambda = STEPS_STELLAR(filter_labels=['GALEX_FUV','ACS_F435W','ACS_F555W','SUBARU_Ks'],steps_bounds=steps_bounds,$
	                      z_shift=0.d0,_EXTRA=_extra_savefits)
  dl07     = dl07_templates(filter_labels=filter_labels,z_shift=galaxy[i].redshift,_EXTRA=_extra_savefits)
	                      
  L_star_abs_table=mrdfits(L_star_abs_table_loc+'L_star_abs_model_table_z_'+string(galaxy[i].redshift,f='(f4.2)')+'.fits',1)

  exp_tau=!null
  vectors = {tauB_f_vectors: reform(chain[nsteps,-Nchain:-1]), rdisk_vectors: reform(chain[nsteps+1,-Nchain:-1]), F_vectors: reform(chain[nsteps+2,-Nchain:-1]), $
          cosi_vectors: reform(chain[nsteps+3,-Nchain:-1]), rbulge_vectors: reform(chain[nsteps+4,-Nchain:-1])}
  models = lightning_models_vector(steps=steps,vectors=vectors,/tuffs,/L_star_abs_model_table,L_star_abs_table=L_star_abs_table,exp_tau=exp_tau,_EXTRA=_extra_savefits)
  models[where(finite(models,/NaN),/null)] = 0.0
  models_alambda = lightning_models_vector(steps=steps_alambda,vectors=vectors,/tuffs,/L_star_abs_model_table,L_star_abs_table=L_star_abs_table,_EXTRA=_extra_savefits)
  models_alambda[where(finite(models_alambda,/NaN),/null)] = 0.0
  de_models = dl07_sed_vector(dl07,alpha=reform(chain[nsteps+5,-Nchain:-1]),umin=reform(chain[nsteps+6,-Nchain:-1]),umax=reform(chain[nsteps+7,-Nchain:-1]),$
                           gam=reform(chain[nsteps+8,-Nchain:-1]),qPAH=reform(chain[nsteps+9,-Nchain:-1]),filter_labels=filter_labels,Lbol=dust_LTIR)
  de_LTIR=dust_LTIR

  
  Lnu_mod_star = total(models * (rebin(reform(chain[0:(nsteps-1),-Nchain:-1],1,Nsteps,Nchain),Nfilters+1,Nsteps,Nchain)),2)
  Lnu_mod = Lnu_mod_star[0:(Nfilters-1),*] + rebin(Lnu_mod_star[Nfilters,*]/reform(de_LTIR,1,nchain),Nfilters,Nchain)*de_models
  Lnu_unred = total((rebin(reform(chain[0:(nsteps-1),-Nchain:-1],1,Nsteps,Nchain),Nfilters,Nsteps,Nchain))*$
              (rebin(reform(steps.mean_Lnu,Nfilters,Nsteps,1),Nfilters,Nsteps,Nchain)),2)

  Lnu_mod_star_alambda = total(models_alambda * (rebin(reform(chain[0:(nsteps-1),-Nchain:-1],1,Nsteps,Nchain),5,Nsteps,Nchain)),2)
  Lnu_unred_alambda = total((rebin(reform(chain[0:(nsteps-1),-Nchain:-1],1,Nsteps,Nchain),4,Nsteps,Nchain))*$
                 (rebin(reform(steps_alambda.mean_Lnu,4,Nsteps,1),4,Nsteps,Nchain)),2)
  Alambda=-2.5*alog10(Lnu_mod_star_alambda[0:3,*]/Lnu_unred_alambda)


  min_chi2=(where(chi2_chain[-Nchain:-1] eq min(chi2_chain[-Nchain:-1])))[0]
  nsteps=n_elements(steps.bounds)-1
  min_chi_chain=chain[*,-Nchain:-1]
  min_chi_chain=min_chi_chain[*,min_chi2]

  Lhi_res_unred=total(steps.lnu * rebin(reform(min_chi_chain[0:(nsteps-1)],1,nsteps),n_elements(steps.wave_rest),nsteps),2)
  steps_lnu_red=steps.lnu*rebin(reform(exp_tau[*,min_chi2],n_elements(steps.wave_rest),1),n_elements(steps.wave_rest),nsteps)
  Lhi_res_red=total(steps_lnu_red * rebin(reform(min_chi_chain[0:(nsteps-1)],1,nsteps),n_elements(steps.wave_rest),nsteps),2)

  dust_Lnu_pow_v03=dl07_spec(dl07,alpha=(reform(chain[nsteps+5,-Nchain:-1]))[min_chi2],umin=(reform(chain[nsteps+6,-Nchain:-1]))[min_chi2],$
     umax=(reform(chain[nsteps+7,-Nchain:-1]))[min_chi2],gam=(reform(chain[nsteps+8,-Nchain:-1]))[min_chi2],$
     qPAH=(reform(chain[nsteps+9,-Nchain:-1]))[min_chi2],Lbol=Lbol_v03)
  dust_lnu_v03=dust_Lnu_pow_v03*Lnu_mod_star[Nfilters,min_chi2]/Lbol_v03

  total_wave=[steps.wave_obs,(dl07.wave_obs)[where(dl07.wave_obs gt (steps.wave_obs)[-9])]]
  total_wave=total_wave[sort(total_wave)]

  total_lnu_star_v03=total_wave*0.d0
  total_lnu_dust_v03=total_wave*0.d0
  total_lnu_star_v03[where(total_wave le max(steps.wave_obs))]=interpol(Lhi_res_red,steps.wave_obs,$
     total_wave[where(total_wave le max(steps.wave_obs))])
  total_lnu_dust_v03[where(total_wave ge min(dl07.wave_obs))]=interpol(dust_lnu_v03,dl07.wave_obs,$
     total_wave[where(total_wave ge min(dl07.wave_obs))])
  total_lnu_v03=total_lnu_star_v03+total_lnu_dust_v03


  out[i].Galaxy_ID	                      = gal_id
  out[i].RAJ2000	                      = galaxy[i].raj2000
  out[i].DECJ2000	                      = galaxy[i].decj2000
  out[i].redshift	                      = galaxy[i].redshift
  out[i].metallicity                      = 0.020
  out[i].N_J                              = galaxy[i].N_J    
  out[i].N_J_ERR                          = galaxy[i].N_J_ERR
  out[i].Q_J                              = galaxy[i].Q_J    
  out[i].Q_J_ERR                          = galaxy[i].Q_J_ERR
  out[i].filter_labels                    = filter_labels
  out[i].WAVE_FILTERS                     = steps.WAVE_FILTERS
  out[i].Lnu_obs	                      = lnu_obs[*,i]
  out[i].Lnu_unc	                      = lnu_unc[*,i]
  out[i].Lnu_mod                          = Lnu_mod
  out[i].Lnu_mod_unred                    = Lnu_unred
  out[i].LTIR_mod                         = reform(lnu_mod_star[Nfilters,*])
  out[i].LFUV_mod                         = reform(lnu_mod_star_alambda[0,*])
  out[i].wave_hires_dustmod               = dl07.wave_obs
  out[i].wave_hires_starmod               = steps.wave_obs
  out[i].wave_hires_totalmod              = total_wave
  out[i].lnu_hires_dustmod                = dust_lnu_v03
  out[i].lnu_hires_starmod                = Lhi_res_red
  out[i].lnu_hires_starmod_unred          = Lhi_res_unred
  out[i].lnu_hires_totalmod               = total_lnu_v03
  out[i].Av		                          = reform(Alambda[2,*])
  out[i].Ab		                          = reform(Alambda[1,*])
  out[i].Afuv		                      = reform(Alambda[0,*])
  out[i].Anir		                      = reform(Alambda[3,*])
  out[i].Nsteps		                      = Nsteps
  out[i].Steps_bounds	                  = steps_bounds
  out[i].chisqr_lightning                 = chi2_chain[-Nchain:-1]
  out[i].SFH[0:(Nsteps-1),*]              = chain[0:(nsteps-1),-Nchain:-1]
  out[i].steps_Mstar_coeff[0:(Nsteps-1)]  = steps_mstar
  out[i].Mstar                            = total(chain[0:(nsteps-1),-Nchain:-1]*rebin(steps_mstar,nsteps,Nchain),1)
  out[i].steps_Mstar[0:(Nsteps-1),*]      = chain[0:(nsteps-1),-Nchain:-1]*rebin(steps_mstar,nsteps,Nchain)
  out[i].tauB_F 		                  = chain[nsteps,-Nchain:-1] 
  out[i].rdisk		                      = chain[nsteps+1,-Nchain:-1] 
  out[i].F			                      = chain[nsteps+2,-Nchain:-1] 
  out[i].cosi		                      = chain[nsteps+3,-Nchain:-1] 
  out[i].rbulge		                      = chain[nsteps+4,-Nchain:-1] 
  out[i].alpha 			                  = chain[nsteps+5,-Nchain:-1]   
  out[i].u_min			                  = chain[nsteps+6,-Nchain:-1] 
  out[i].u_max			                  = chain[nsteps+7,-Nchain:-1] 
  out[i].gamma			                  = chain[nsteps+8,-Nchain:-1]   
  out[i].q_pah			                  = chain[nsteps+9,-Nchain:-1]   
 
  if keyword_set(prior) then out[i].chisqr_prior = chi2_prior_chain[-Nchain:-1]
 
  print,i
  
endfor

spawn,'rm -r '+outfolder+file_name+'.fits'
mwrfits,out, outfolder+file_name+'.fits',/create

end