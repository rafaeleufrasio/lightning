function dl07_templates,filter_labels=filter_labels,z_shift=z_shift,lightning_folder=lightning_folder

; MODIFICATION HISTORY:
;   Written, Rafael Eufrasio, March 2020
;   Made folder location lightning_folder vs fix location in code, Keith Doore, April 27th 2020
;   Changed "folder" keyword to "lightning_folder" for consistence between functions/procedures for _REF_EXTRA, Keith Doore, May 6th 2020
;   Corrected (1+z) term, Keith Doore, May 11th 2020

if (n_elements(filter_labels) eq 0) then $
  filter_labels=['GALEX_FUV','UVOT_W2','UVOT_M2','GALEX_NUV','UVOT_W1','SDSS_u','SDSS_g','SDSS_r','SDSS_i','SDSS_z',$
                 '2MASS_J','2MASS_H','2MASS_Ks','W1','SPITZER_I1','SPITZER_I2','W2','W3','W4','HAWC+A','HAWC+B','HAWC+C','HAWC+D','HAWC+E']
if (n_elements(z_shift) eq 0) then z_shift = 0.0
if (n_elements(lightning_folder) eq 0) then Lightning='~/Lightning/' else Lightning=lightning_folder

Ugrid=['0.10','0.15','0.20','0.30','0.40','0.50','0.70','0.80','1.00','1.20','1.50','2.50','3.00','4.00',$
         '5.00','7.00','8.00','12.0','15.0','20.0','25.0','1e2','3e2','1e3','3e3','1e4','3e4','1e5','3e5']

mod_grid = []                     &  qPAH_grid = []
mod_grid = [mod_grid,'MW3.1_00']  &  qPAH_grid = [qPAH_grid,0.47d-2]  ;qPAH = 0.47 %
mod_grid = [mod_grid,'MW3.1_10']  &  qPAH_grid = [qPAH_grid,1.12d-2]  ;qPAH = 1.12 %
mod_grid = [mod_grid,'MW3.1_20']  &  qPAH_grid = [qPAH_grid,1.77d-2]  ;qPAH = 1.77 %
mod_grid = [mod_grid,'MW3.1_30']  &  qPAH_grid = [qPAH_grid,2.50d-2]  ;qPAH = 2.50 %
mod_grid = [mod_grid,'MW3.1_40']  &  qPAH_grid = [qPAH_grid,3.19d-2]  ;qPAH = 3.19 %
mod_grid = [mod_grid,'MW3.1_50']  &  qPAH_grid = [qPAH_grid,3.90d-2]  ;qPAH = 3.90 %
mod_grid = [mod_grid,'MW3.1_60']  &  qPAH_grid = [qPAH_grid,4.58d-2]  ;qPAH = 4.58 %
;mod_grid = [mod_grid,'LMC2_00']   &  qPAH_grid = [qPAH_grid,0.75d-2]  ;qPAH = 0.75 %
;mod_grid = [mod_grid,'LMC2_05']   &  qPAH_grid = [qPAH_grid,1.49d-2]  ;qPAH = 1.49 %
;mod_grid = [mod_grid,'LMC2_10']   &  qPAH_grid = [qPAH_grid,2.37d-2]  ;qPAH = 2.37 %
;mod_grid = [mod_grid,'smc']       &  qPAH_grid = [qPAH_grid,0.10d-2]  ;qPAH = 0.10 %

Nfilters=filter_labels.length
Udbl = double(Ugrid)

n_U = Ugrid.length
n_mod = mod_grid.length

wave = dblarr(1001)
; => lambdas are all exactly the same for all delta functions
Lnu  = dblarr(1001,n_U,n_mod)
Lbol = dblarr(n_U,n_mod)

for uu=0,(n_U-1) do begin
  U = Ugrid[uu]
  for mm=0,(n_mod-1) do begin
    modl = mod_grid[mm]
    readcol,Lightning+'dust_emission/DL07/DL07spec/U'+U+'/U'+U+'_'+U+'_'+modl+'.txt',wave_uumm,nuLnu_uumm,jnu_min_uumm,/silent
    ;readcol,'~/Dropbox/DL07/DL07spec/U'+U+'/U'+U+'_'+U+'_'+modl+'.txt',wave_uumm,nuLnu_uumm
    ;readcol,'~/Dropbox/DL07/DL07spec/U'+U+'/U'+U+'_'+U+'_'+modl+'.txt',wave_uumm,nuLnu_uumm,/silent
    l0 = ((where(wave_uumm eq 1e4,/null))[-1])[0]
    ;lambda = wave_uumm[l0:*]
    ;print,'  ',U,modl,lambda.length
    wave = reverse(wave_uumm[l0:*]) ; wave is increasing, in microns
    nu = 1.d4* !cv.clight/wave
    Lnu[*,uu,mm] = reverse(nuLnu_uumm)/nu
    Lbol[ uu,mm] = abs(trap_int(nu,Lnu[*,uu,mm]))
  endfor
endfor


nu_rest   = temporary(nu)             ; restframe frequency
nu_obs    = nu_rest/(1.d0+z_shift)    ; observed frequency
wave_rest = temporary(wave)           ; restframe wavelength
wave_obs  = wave_rest*(1.d0+z_shift)  ; observed wavelength
Lnu_rest  = temporary(Lnu)            ; restframe Lnu spectrum
Lnu       = Lnu_rest * (1.d0+z_shift) ; observed Lnu spectrum, where Lnu = 4*!dpi*DL^2*Fnu
GET_FILTERS,filter_labels,nu_obs,Filters,Lnu=[1.d0] # [0*(nu_obs) + 3631],     $ ; flat AB=0 source
            filters_dir=Lightning+'Filters/',mean_wave=mean_wave;,/plot_filters;$
            ;,mean_nu=mean_nu,mean_Lnu=mean_Lnu,mean_Lnu=mean_Lnu,sigma_wave=sigma_wave

wave_filters = reform(mean_wave)

Lnu_kkii = dblarr(nu_obs.length,n_u*n_mod)
for kk=0,(nu_obs.length-1L) do Lnu_kkii[kk,*] = (Lnu[kk,*,*])[lindgen(n_u*n_mod)] 

GET_FILTERS,filter_labels,nu_obs,Filters,Lnu=transpose(Lnu_kkii),$
            filters_dir=Lightning+'Filters/',mean_Lnu=mean_Lnu;,$
            ;/plot_filters;,mean_wave=mean_wave,mean_nu=mean_nu,sigma_wave=sigma_wave

mean_Lnu[where(finite(mean_Lnu,/naN),/null)] = 0.0 ;filters not fully included in DL07 wave range produce mean_Lnu = NaN

mean_Lnu_kii = transpose(temporary(mean_Lnu)); mean_Lnukii[Nfilter,n_u*n_mod]
mean_Lnu = dblarr(Nfilters,n_u,n_mod)
for k=0,(Nfilters-1L) do mean_Lnu[k,*,*] = mean_Lnu_kii[k,lindgen(n_u,n_mod)]

dl07 = {filter_labels : filter_labels,$
        wave_filters  : wave_filters, $
        U             : Ugrid,        $
        model         : mod_grid,     $
        qPAH          : qPAH_grid,    $
        mean_Lnu      : mean_Lnu,     $
        wave_rest     : wave_rest,    $ ; restframe wavelength
        wave_obs      : wave_obs,     $
        Lnu           : Lnu,          $ ; observed Lnu spectrum, where Lnu = 4*!dpi*DL^2*Fnu = Lnu_rest * (1+z) 
        Lbol          : Lbol,         $
        Filters       : Filters,      $
        z_shift       : z_shift}

return,dl07

end




function dl07_sed_vector,dl07,alpha=alpha,umin=umin,umax=umax,gam=gam,modeltype=modeltype,qPAH=qPAH,filter_labels=filter_labels,$
                  Lbol=Lbol,Lpow=Lpow,Ldel=Ldel,mean_Lnu_pow=mean_Lnu_pow,mean_Lnu_del=mean_Lnu_del

; output:
;   mean_Lnu
;   
; optional output:
;   Lbol : integrated luminosity
;
; MODIFICATION HISTORY:
;   Written, Rafael Eufrasio, March 2020
;   added gamma and del component, Sept 30 2020, Keith Doore and Rafael Eufrasio


if n_elements(alpha) eq 0 then alpha = 2.0  
if n_elements(Umin)  eq 0 then Umin  = 0.1  
if n_elements(Umax)  eq 0 then Umax  = 1.e5 
if n_elements(gam)   eq 0 then gam   = 1.0  
if n_elements(qPAH)  eq 0 then qPAH  = 0.025

if n_elements(alpha) ne n_elements(Umin) or n_elements(alpha) ne n_elements(Umax) or n_elements(alpha) ne n_elements(qPAH) or $
  n_elements(alpha) ne n_elements(gam) or n_elements(Umin) ne n_elements(Umax) or n_elements(Umin) ne n_elements(qPAH) or $
  n_elements(Umin) ne n_elements(gam) or n_elements(Umax) ne n_elements(qPAH) or n_elements(Umax) ne n_elements(gam) or $
  n_elements(qPAH) ne n_elements(gam) then message,'Input parameters must have same length for vectorization'
nn =  n_elements(alpha) 


Nfilters = n_elements(dl07.filter_labels)
N_u      = n_elements(dl07.U)
Nmod     = n_elements(dl07.model)

qPAH_vec = dl07.qPAH

;(dl07.Lbol)[uu,mm]       --- N_u x Nmod
;(dl07.mean_Lnu)[k,uu,mm] --- Nfilters x N_u x Nmod
Lbol_temp = dblarr(N_u,nn)
mean_Lnu = dblarr(N_u,Nfilters,nn)
for uu=0,(N_u-1) do begin 
  Lbol_temp[uu,*] = interpol(reform((dl07.Lbol)[uu,*]),qPAH_vec,qPAH) ;Lbol[uu]
  mean_Lnu_uu = reform((dl07.mean_Lnu)[*,uu,*])                 ;mean_Lnu_uu[k,mm]
  for k=0,(Nfilters-1) do begin 
    mean_Lnu_kuu = reform(mean_Lnu_uu[k,*])                     ;mean_Lnu_kuu[mm]
    mean_Lnu[uu,k,*] = interpol(mean_Lnu_kuu,qPAH_vec,qPAH) 
    ;mean_Lnu[uu,k] = finterpol(mean_Lnu_kuu,qPAH_vec,qPAH)
  endfor
endfor

U = rebin(double(dl07.U),n_elements(dl07.U),nn)

;power-law component
fU = U^(rebin(reform(-1*alpha,1,nn),n_elements(dl07.U),nn))  
Lpow=dblarr(nn)
mean_Lnu_pow = dblarr(Nfilters,nn)
for i=0,(nn-1) do begin
  fU[*,i] /= trap_int(U[*,i],fU[*,i],xrange=[Umin[i],Umax[i]],/pow)
  Lpow[i] = gam[i]*trap_int(U[*,i],Lbol_temp[*,i]*fU[*,i],xrange=[Umin[i],Umax[i]],/pow)
  for k=0,(Nfilters-1) do mean_Lnu_pow[k,i] = gam[i]*trap_int(U[*,i],mean_Lnu[*,k,i]*fU[*,i],xrange=[Umin[i],Umax[i]],/pow)
endfor
mean_Lnu_pow[where(mean_Lnu_pow lt 1.d-150,/null)] = 0.d0


;delta function at Umin
Ldel=dblarr(nn)
mean_Lnu_del = dblarr(Nfilters,nn)
for i=0,(nn-1) do begin
  Ldel[i] = (1.d0-gam[i])*interpol(Lbol_temp[*,i],U[*,i],Umin[i])
  for k=0,(Nfilters-1) do mean_Lnu_del[k,i] = (1.d0-gam[i])*interpol(mean_Lnu[*,k,i],U[*,i],Umin[i])
endfor
mean_Lnu_del[where(mean_Lnu_del lt 1.d-150,/null)] = 0.d0

;combine delta function and power-law
Lbol     = Ldel         + Lpow
mean_Lnu = mean_Lnu_del + mean_Lnu_pow


return,mean_Lnu

end






function dl07_spec_vector,dl07,alpha=alpha,umin=umin,umax=umax,gam=gam,qPAH=qPAH,modeltype=modeltype,$
                   Lbol=Lbol,Lpow=Lpow,Ldel=Ldel,Lnu_pow=Lnu_pow,Lnu_del=Lnu_del
;output:
;  Lnu
;  
;optional output:
;  Lbol : integrated luminosity
;  Lpow: power law component bolbometric luminosity
;  Ldel:
;
; MODIFICATION HISTORY:
;   Written, Rafael Eufrasio, March 2020
;   added gamma and del component, Sept 30 2020, Keith Doore and Rafael Eufrasio

  if n_elements(alpha) eq 0 then alpha = 2.0  
  if n_elements(Umin)  eq 0 then Umin  = 0.1  
  if n_elements(Umax)  eq 0 then Umax  = 1.e5 
  if n_elements(gam)   eq 0 then gam   = 1.0  
  if n_elements(qPAH)  eq 0 then qPAH  = 0.025


  if n_elements(alpha) ne n_elements(Umin) or n_elements(alpha) ne n_elements(Umax) or n_elements(alpha) ne n_elements(qPAH) or $
    n_elements(alpha) ne n_elements(gam) or n_elements(Umin) ne n_elements(Umax) or n_elements(Umin) ne n_elements(qPAH) or $
    n_elements(Umin) ne n_elements(gam) or n_elements(Umax) ne n_elements(qPAH) or n_elements(Umax) ne n_elements(gam) or $
    n_elements(qPAH) ne n_elements(gam) then message,'Input parameters must have same length for vectorization'
  nn =  n_elements(alpha) 


  Nwave = n_elements(dl07.wave_obs)
  N_u   = n_elements(dl07.U)
  Nmod  = n_elements(dl07.model)

  qPAH_vec = dl07.qPAH

  ;(dl07.Lbol)[uu,mm]        --- N_u x Nmod
  ;(dl07.mean_Lnu)[kk,uu,mm] --- Nwave x N_u x Nmod
  Lbol_temp = dblarr(N_u,nn)
  Lnu = dblarr(N_u,Nwave,nn)
  for uu=0,(N_u-1) do begin
    Lbol_temp[uu,*] = interpol(reform((dl07.Lbol)[uu,*]),qPAH_vec,qPAH) ;Lbol[uu]
    Lnu_uu = reform((dl07.Lnu)[*,uu,*])                          ;Lnu_uu[k,mm]
    for kk=0,(Nwave-1) do begin
      Lnu_kuu = reform(Lnu_uu[kk,*])                             ;Lnu_kuu[mm]
      Lnu[uu,kk,*] = interpol(Lnu_kuu,qPAH_vec,qPAH)
      ;mean_Lnu[uu,k] = finterpol(mean_Lnu_kuu,qPAH_vec,qPAH)
    endfor
  endfor

  U = rebin(double(dl07.U),n_elements(dl07.U),nn)
  
  ;power-law component
  fU = U^(rebin(reform(-1*alpha,1,nn),n_elements(dl07.U),nn))  
  Lpow=dblarr(nn)
  Lnu_pow = dblarr(Nwave,nn)
  for i=0,(nn-1) do begin
    fU[*,i] /= trap_int(U[*,i],fU[*,i],xrange=[Umin[i],Umax[i]],/pow)
    Lpow[i] = gam[i]*trap_int(U[*,i],Lbol_temp[*,i]*fU[*,i],xrange=[Umin[i],Umax[i]],/pow)
    for kk=0,(Nwave-1) do Lnu_pow[kk,i] = gam[i]*trap_int(U[*,i],Lnu[*,kk,i]*fU[*,i],xrange=[Umin[i],Umax[i]],/pow)
  endfor


  ;delta function at Umin
  Ldel=dblarr(nn)
  Lnu_del = dblarr(Nwave,nn)
  for i=0,(nn-1) do begin
    Ldel[i] = (1.d0-gam[i])*interpol(Lbol_temp[*,i],U[*,i],Umin[i])
    for kk=0,(Nwave-1) do Lnu_del[kk,i] = (1.d0-gam[i])*interpol(Lnu[*,kk,i],U[*,i],Umin[i])
  endfor

  ;combine delta function and power-law
  Lbol = Ldel    + Lpow
  Lnu  = Lnu_del + Lnu_pow

  return,Lnu

end
