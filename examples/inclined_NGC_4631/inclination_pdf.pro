pro inclination_pdf,q,q_err,folder,inc_bins,inc_pdf,NTRIALS=NTRIALS,$
           NQGAUSSIAN=NQGAUSSIAN,NBINS_INC=NBINS_INC,RAN_Q=RAN_Q,RAN_INC=RAN_INC,$
           GEDIST=GEDIST,UNITS=UNITS,GAM_EPS=GAM_EPS
;+
; NAME:
;   INCLINATION_PDF
; PURPOSE:
;   To create the probability distribution function (PDF) of inclination as described 
;     in Doore et al. (2021)
; EXPLANATION: 
;   Creates a distribution of inclination from the input axis ratio q and error,
;     assuming the input q has a Gaussian distribution, using a Monte Carlo method. 
;     Can use the distributions of asymmetry and intrinsic thickness of spiral galaxies 
;     from Rodriguez & Padilla (2013).
;
; CALLING SEQUENCE:
;   inclination_pdf, q, q_err, folder, inc_bins, inc_pdf, [NTRIALS=NTRIALS, $
;           NQGAUSSIAN=NQGAUSSIAN, NBINS_INC=NBINS_INC, RAN_Q=RAN_Q, RAN_INC=RAN_INC,$
;           GEDIST=GEDIST,UNITS=UNITS,GAM_EPS=GAM_EPS]
;
; INPUTS:
;   q      - a N element vector containing the values of axis ratio to be considered
;   q_err  - a N element vector containing the corresponding axis ratio error values
;   folder - a string containing the directory location of the files of the gamma and epsilon 
;              distributions from Rodriguez & Padilla (2013)
;
; OPTIONAL INPUTS:
;   NTRIALS       - a single value containing the number of trials to run the Monte Carlo 
;                     simulation (default = 1e7)
;   NQGAUSSIAN    - a single value containing the number of elements to use in the Gaussian 
;                     distribution of q (default = 1e6)
;   NBINS_INC     - a single value containing the number of bins for which to determine
;                     the distribution of inclination (default = 180)
;   GEDIST        - a string containing the distribution to use for gamma and epsilon from 
;                     Figure 11 of Rodriguez & Padilla (2013). Options are the color of the 
;                     distributions: red (default), blue, or green
;   UNITS         - a string containing the desired units of the output inc_bins. Options
;                     are degrees (default), radians, cosi (cosine of the inclination)
;   GAM_EPS       - a 2 element vector containing the fixed values of gamma and epsilon, 
;                     respectively. Values for both must be between 0 and 1. If set, then 
;                     GEDIST input will be ignored.
;
; OUTPUTS:
;   inc_bins - an NBINS_INC element vector containing the inclinations in the units given 
;                in the optional input UNITS corresponding to the PDF values. Points are 
;                at the center of each bin
;   inc_pdf  - an NBINS_INC x N element array containing the normalized PDF values for 
;                the inclination of each input galaxies axis ratio. 
;
; OPTIONAL OUTPUTS:
;   RAN_Q   - a NTRIALS element vector containing the q values from the Monte Carlo simulation
;   RAN_INC - a NTRIALS element vector containing the corresponding inclination values to
;               q from the Monte Carlo simulation
;
; EXAMPLE USAGE:
;   IDL> inclination_pdf,0.3,0.01,'~/YOUR_FOLDER/',inc_bins,inc_pdf
;   IDL> help, inc_bins,inc_pdf
;   IDL> plot(inc_bins,inc_pdf)
;
; REVISON HISTORY:
;   Written by K. Doore, 9/01/2020
;   Updated by K. Doore, 3/03/2021
;     - Replaced the optional keyword RADIANS with optional input of UNITS to allow for
;         specifying the units as DEGREES, RADIANS, or COSI. Changed INC_BINWIDTH to be
;         NBINS_INC to be general for the variation in the possible units
;     - Added the optional input of GAM_EPS to allow for the values of gamma and epsilon
;         to be fixed to a specific value
;     - Updated the descriptions of each input/output to be more detailed
;-

  Compile_opt idl2
  On_error,2

; Check for incorrect inputs
  if size(q,/type) lt 2 or size(q,/type) gt 5 then begin
    print,'q is incorrect data type'
    return
  endif
  if size(q,/n_dim) gt 1 then begin
    print,'q must be 1 dimensional'
    return
  endif

  if min(q) lt 0 or max(q) gt 1 then begin
    print,'q must be between 0 and 1'
    return
  endif

  if size(q_err,/type) lt 2 or size(q_err,/type) gt 5 then begin
    print,'q_err is incorrect data type'
    return
  endif 
  if size(q_err,/n_dim) gt 1 then begin
    print,'q_err must be 1 dimensional'
    return
  endif
  if n_elements(q_err) ne n_elements(q) then begin
    print,'q and q_err need same number of elements'
    return
  endif
  if min(q_err) le 0 then begin
    print,'q_err must be greater than 0'
    return
  endif

  if size(folder,/type) ne 7 then begin
    print,'folder is not of type STRING'
    return
  endif 

  if n_elements(NTRIALS) gt 0 then begin
    if size(NTRIALS,/type) lt 2 or size(NTRIALS,/type) gt 5 then begin
      print,'NTRIALS is incorrect data type'
      return
    endif
    if min(NTRIALS) lt 0 then begin
      print,'NTRIALS must be greater than 0'
      return
    endif
  endif else begin
    NTRIALS = 1e7
  endelse

  if n_elements(NQGAUSSIAN) gt 0 then begin
    if size(NQGAUSSIAN,/type) lt 2 or size(NQGAUSSIAN,/type) gt 5 then begin
      print,'NQGAUSSIAN is incorrect data type'
      return
    endif
    if min(NQGAUSSIAN) lt 0 then begin
      print,'NQGAUSSIAN must be greater than 0'
      return
    endif
  endif else begin
    NQGAUSSIAN = 1e6
  endelse
  
  if n_elements(NTRIALS) lt n_elements(NQGAUSSIAN) then begin
    print,'NTRIALS needs to be larger than NQGAUSSIAN'
    return
  endif

  if n_elements(NBINS_INC) gt 0 then begin
    if size(NBINS_INC,/type) lt 2 or size(NBINS_INC,/type) gt 5 then begin
      print,'NBINS_INC is incorrect data type'
      return
    endif
    if min(NBINS_INC) le 0 then begin
      print,'NBINS_INC must be positive'
      return
    endif
  endif else begin
    NBINS_INC = 180
  endelse
  
  if n_elements(GEDIST) gt 0 then begin
    if size(GEDIST,/type) ne 7 then begin
      print,'GEDIST must be of type string'
      return
    endif
    if ~STRCMP(GEDIST,'red',/fold_case) and ~STRCMP(GEDIST,'blue',/fold_case) and $
       ~STRCMP(GEDIST,'green',/fold_case) then begin
      print,'GEDIST must be either red, blue, or green'
      return
    endif
  endif else begin
    if n_elements(GAM_EPS) eq 0 then GEDIST = 'red'
  endelse

  if n_elements(UNITS) gt 0 then begin
    if size(UNITS,/type) ne 7 then begin
      print,'UNITS must be of type string'
      return
    endif
    if ~STRCMP(UNITS,'degrees',/fold_case) and ~STRCMP(UNITS,'radians',/fold_case) and $
       ~STRCMP(UNITS,'cosi',/fold_case) then begin
      print,'UNITS must be either degrees, radians, or cosi'
      return
    endif
  endif else begin
    UNITS = 'degrees'
  endelse

  if n_elements(GAM_EPS) gt 0 then begin
    if size(GAM_EPS,/type) lt 2 or size(GAM_EPS,/type) gt 5 then begin
      print,'GAM_EPS is incorrect data type'
      return
    endif
    if n_elements(GAM_EPS) ne 2 then begin
      print,'GAM_EPS must be a two element vector'
      return
    endif
    if GAM_EPS[0] lt 0 or GAM_EPS[0] gt 1 then begin
      print,'The first element of GAM_EPS (gamma) must be between 0 and 1'
      return
    endif
    if GAM_EPS[1] lt 0 or GAM_EPS[1] gt 1 then begin
      print,'The second element of GAM_EPS (epsilon) must be between 0 and 1'
      return
    endif
  endif


; compute the cumulative PDFs of gamma and epsilon
  if n_elements(GEDIST) gt 0 then begin
    if STRCMP(GEDIST,'red',/fold_case) then begin
      readcol,folder+'gamma_spirals_red.txt',gam_val,gam_PDF,f='f,f',/silent
      readcol,folder+'epsilon_spirals_red.txt',eps_val,eps_PDF,f='f,f',/silent
    endif else if STRCMP(GEDIST,'blue',/fold_case) then begin
      readcol,folder+'gamma_spiralsfracDeV_blue.txt',gam_val,gam_PDF,f='f,f',/silent
      readcol,folder+'epsilon_spiralsfracDeV_blue.txt',eps_val,eps_PDF,f='f,f',/silent
    endif else if STRCMP(GEDIST,'green',/fold_case) then begin
      readcol,folder+'gamma_spiralsPS08_green.txt',gam_val,gam_PDF,f='f,f',/silent
      readcol,folder+'epsilon_spiralsPS08_green.txt',eps_val,eps_PDF,f='f,f',/silent
    endif
    
    gam=[0:1:0.0001]
    eps=[0:1:0.0001]
    PDF_gam=interpol(gam_PDF,gam_val,gam)
    PDF_eps=interpol(eps_PDF,eps_val,eps)
    CUMPDF_gam=total(pdf_gam,/cum)/max(total(pdf_gam,/cum))
    CUMPDF_eps=total(pdf_eps,/cum)/max(total(pdf_eps,/cum))
  

; run the Monte Carlo simulation
    ran_gam=interpol(gam,cumpdf_gam,randomu(seed,NTRIALS))
    ran_eps=interpol(eps,cumpdf_eps,randomu(seed,NTRIALS))
  endif else begin
    ran_gam=replicate(GAM_EPS[0],NTRIALS)
    ran_eps=replicate(GAM_EPS[1],NTRIALS)
  endelse

  ran_phi=randomu(seed,NTRIALS)*2*!pi
  ran_cosi=randomu(seed,NTRIALS)
  ran_q=sqrt(ran_cosi^2*(1-ran_gam^2)+ran_gam^2)*(1-ran_eps*cos(2*ran_phi))
  
  keep=where(ran_q le 1 and ran_q ge 0)

  if STRCMP(UNITS,'degrees',/fold_case) then begin
    ran_inc=acos(ran_cosi)/!dtor
    binvals=[0:90:(90/double(NBINS_INC))]
    inc_bins=binvals[1:*]-((90/double(NBINS_INC))/2.d)
  endif else if STRCMP(UNITS,'radians',/fold_case) then begin
    ran_inc=acos(ran_cosi)
    binvals=[0:(!dpi/2.d):((!dpi/2.d)/double(NBINS_INC))]
    inc_bins=binvals[1:*]-(((!dpi/2.d)/double(NBINS_INC))/2.d)
  endif else if STRCMP(UNITS,'cosi',/fold_case) then begin
    ran_inc=ran_cosi
    binvals=[0:1:(1/double(NBINS_INC))]
    inc_bins=binvals[1:*]-((1/double(NBINS_INC))/2.d)
  endif

  ran_inc=ran_inc[keep]  
  ran_q=ran_q[keep]  
  ran_inc=ran_inc[sort(ran_q)]
  ran_q=ran_q[sort(ran_q)]


; compute the inclination PDF
  nbins=n_elements(inc_bins)

  inc_PDF=dblarr(nbins,n_elements(q))
  for i=0,(n_elements(q)-1) do begin
    
    q_dist=randomn(seed,NQGAUSSIAN)*q_err[i]+q[i]
    q_dist=q_dist[where(q_dist ge 0 and q_dist le 1)]
    values=ran_inc[value_locate(ran_q,q_dist,/l64)]
    
    for j=0,(n_elements(binvals)-2) do begin
      temp_num=where(values ge binvals[j] and values lt binvals[j+1],num)
      inc_PDF[j,i]=num
    endfor
    
    inc_PDF[*,i]=inc_pdf[*,i]/INT_TABULATED(inc_bins,inc_pdf[*,i])
  endfor

return

end