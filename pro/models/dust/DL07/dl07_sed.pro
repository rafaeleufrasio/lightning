function dl07_sed, dl07, alpha=alpha, umin=umin, umax=umax, gam=gam, qPAH=qPAH, error_check=error_check, $
                   LTIR=LTIR, pow_LTIR=pow_LTIR, del_LTIR=del_LTIR, pow_mean_Lnu=pow_mean_Lnu, $
                   del_mean_Lnu=del_mean_Lnu
;+
; Name
; ----
;   DL07_SED
;
; Purpose
; -------
;   Generates the dust model SED for a given set of dust model parameters and
;   filters. The dust model is based on the Draine & Li (2007) dust model, which
;   include the parameters of :math:`\alpha`, :math:`U_{\rm min}`, :math:`U_{\rm max}`,
;   :math:`\gamma`, and :math:`q_{\rm PAH}`.
;
; Calling Sequence
; ----------------
;   ::
;
;       mean_Lnu_dust = dl07_sed(dl07 [, alpha = , umin = , umax = , gam = , qPAH = , $
;                                /error_check, LTIR=LTIR, pow_LTIR=pow_LTIR, $
;                                del_LTIR=del_LTIR, pow_mean_Lnu=pow_mean_Lnu, $
;                                del_mean_Lnu=del_mean_Lnu])
;
; Input
; -----
;   ``dl07`` : structure
;       A structure containing the spectra, SEDs, and dust parameters for 
;       the DL07 model. (See ``dl07_models.pro`` for details and contents.)
;
; Optional Inputs
; ---------------
;   ``alpha`` : int, float, or double array(Nmodels)
;       The exponent of the power-law distribution of heating starlight
;       intensities between ``umin`` and ``umax``. (Default = ``2.0``)
;   ``umin`` : int, float, or double array(Nmodels)
;       The minimum radiation field intensity of the diffuse ISM radiation
;       field heating the dust. (Default = ``1.0``)
;   ``umax`` : int, float, or double array(Nmodels)
;       The maximum radiation field intensity of the power-law distribution
;       of heating starlight intensities. (Default = ``3e5``)
;   ``gam`` : int, float, or double array(Nmodels)
;       The fraction of the dust mass exposed to the power-law distribution of
;       starlight intensities.  (Default = ``0.1``)
;   ``qPAH`` : float or double array(Nmodels)
;       The PAH index.  (Default = ``0.025``)
;   ``error_check`` : flag
;       If set, all inputs are checked for errors. Otherwise, all inputs are
;       assumed to be of correct format.
;
; Output
; ------
;   ``mean_Lnu_dust`` : double array(Nfilters, Nmodels)
;       The mean luminosity of each filter and set of dust model parameters
;       :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;
; Optional Outputs
; ----------------
;   ``LTIR`` : double array(Nmodels)
;       The total integrated IR luminosity (i.e., the bolometric luminosity)
;       of the dust model for each set of model parameters :math:`[L_\odot]`.
;   ``pow_LTIR`` : double array(Nmodels)
;       The power-law component of ``LTIR`` :math:`[L_\odot]`.
;   ``del_LTIR`` : double array(Nmodels)
;       The delta-function component of ``LTIR`` :math:`[L_\odot]`.
;   ``pow_mean_Lnu`` : double array(Nfilters, Nmodels)
;       The power-law component of ``mean_Lnu_dust`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;   ``del_mean_Lnu`` : double array(Nfilters, Nmodels)
;       The delta-function component of ``mean_Lnu_dust`` :math:`[L_\odot\ {\rm Hz}^{-1}]`.
;
; Note
; ----
;   Defaults will only be set if the optional ``error_check`` input is set.
;
; Reference
; ---------
;   `Draine, B. T., & Li, A. 2007, ApJ, 657, 810 <https://ui.adsabs.harvard.edu/abs/2007ApJ...657..810D/abstract>`_
;
; Modification History
; --------------------
;   - 2020/03/01: Created (Rafael T. Eufrasio)
;   - 2020/09/30: Added ``gamma`` and ``del`` component (Rafael T. Eufrasio and Keith Doore)
;   - 2022/03/15: Added proper error handling (Keith Doore)
;   - 2022/03/18: Updated documentation (Keith Doore)
;   - 2022/04/12: Allowed for inputs to have degenerate dimensions (Keith Doore)
;   - 2022/04/12: Allowed for inputs to be scalars (Keith Doore)
;   - 2022/04/12: Allowed integer inputs (Keith Doore)
;   - 2022/04/12: Allowed for array of only one optional input to be given (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if keyword_set(error_check) then begin
   if n_elements(dl07) eq 0 then message, 'Variable is undefined: DL07.'
   if size(dl07, /type) ne 8 then message, 'DL07 must be a of type structure.'
  
   if n_elements(alpha) ne 0 then begin
     if size(alpha, /type) lt 2 or size(alpha,/type) gt 5 then $
              message, 'ALPHA must be a of type int, float, or double.'
     if size(reform(alpha), /n_dim) ne 1 then message, 'ALPHA must be a scalar or 1-D array.'
     if min(alpha) lt -10 or max(alpha) gt 4 then message,'ALPHA must contain values between -10 and 4.'
     Nmodels = n_elements(alpha)
   endif
  
   if n_elements(Umin) ne 0 then begin
     if size(Umin, /type) lt 2 or size(Umin,/type) gt 5 then $
              message, 'UMIN  must be a of type int, float, or double.'
     if size(reform(Umin), /n_dim) ne 1 then message, 'UMIN must be a scalar or 1-D array.'
     if min(Umin) lt 0.1d0 or max(Umin) gt 25 then message,'UMIN must contain values between 0.1 and 25.'
     if n_elements(Nmodels) eq 0 then Nmodels = n_elements(Umin)
   endif
  
   if n_elements(Umax) ne 0 then begin
     if size(Umax, /type) lt 2 or size(Umax,/type) gt 5 then $
              message, 'UMAX must be a of type int, float, or double.'
     if size(reform(Umax), /n_dim) ne 1 then message, 'UMAX must be a scalar or 1-D array.'
     if min(Umax) lt 1.d3 or max(Umax) gt 3.d5 then message,'UMAX must contain values between 1e3 and 3e5.'
     if n_elements(Nmodels) eq 0 then Nmodels = n_elements(Umax)
   endif
  
   if n_elements(gam) ne 0 then begin
     if size(gam, /type) lt 2 or size(gam,/type) gt 5 then $
              message, 'GAM must be a of type int, float, or double.'
     if size(reform(gam), /n_dim) ne 1 then message, 'GAM must be a scalar or 1-D array.'
     if min(gam) lt 0 or max(gam) gt 1 then message,'GAM must contain values between 0 and 1.'
     if n_elements(Nmodels) eq 0 then Nmodels = n_elements(gam)
   endif
  
   if n_elements(qPAH) ne 0 then begin
     if size(qPAH, /type) lt 4 or size(qPAH,/type) gt 5 then $
              message, 'QPAH must be a of type float or double.'
     if size(reform(qPAH), /n_dim) ne 1 then message, 'QPAH must be a scalar or 1-D array.'
     if min(qPAH) lt 0.0047d or max(qPAH) gt 0.0458d0 then message,'QPAH must contain values between 0.0047 and 0.0458.'
     if n_elements(Nmodels) eq 0 then Nmodels = n_elements(qPAH)
   endif
  
   if n_elements(Nmodels) ne 0 then begin
     if n_elements(alpha) eq 0 then alpha = replicate(2.d   ,Nmodels)
     if n_elements(Umin ) eq 0 then Umin  = replicate(1.d   ,Nmodels)
     if n_elements(Umax ) eq 0 then Umax  = replicate(3.d5  ,Nmodels)
     if n_elements(gam  ) eq 0 then gam   = replicate(0.1d  ,Nmodels)
     if n_elements(qPAH ) eq 0 then qPAH  = replicate(2.5d-2,Nmodels)
   endif else begin
     Nmodels = 1
     alpha = 2.d
     Umin  = 1.d
     Umax  = 3.d5
     gam   = 0.1d
     qPAH  = 2.5d-2
   endelse
  
   if n_elements(alpha) ne n_elements(Umin)   or n_elements(alpha) ne n_elements(Umax)  or $
      n_elements(alpha) ne n_elements(gam)    or n_elements(alpha) ne n_elements(qPAH) then $
      message,'ALPHA, UMIN, UMAX, GAM, and QPAH must have the same size.'
 endif

; Generate model SED
 Nmodels   = n_elements(alpha)
 Nfilters = n_elements(dl07.filter_labels)
 N_u      = n_elements(dl07.U)
 Nmod     = n_elements(dl07.model)

 qPAH_vec = dl07.qPAH

 ;(dl07.Lbol)[uu,mm]       --- N_u x Nmod
 ;(dl07.mean_Lnu)[k,uu,mm] --- Nfilters x N_u x Nmod
 Lbol = dblarr(N_u, Nmodels)
 mean_Lnu = dblarr(N_u, Nfilters, Nmodels)
 for uu=0, (N_u-1) do begin
   Lbol[uu, *] = interpol(reform((dl07.Lbol)[uu, *]), qPAH_vec, reform(qPAH)) ;Lbol[uu]
   for k=0,(Nfilters-1) do begin
     mean_Lnu[uu, k, *] = interpol((dl07.mean_Lnu)[k, uu, *], qPAH_vec, reform(qPAH))
   endfor
 endfor

 U = rebin(dl07.U, n_elements(dl07.U), Nmodels)

 ;power-law component
 fU = U^(rebin(reform(-1*alpha, 1, Nmodels), n_elements(dl07.U), Nmodels))
 pow_LTIR = dblarr(Nmodels)
 pow_mean_Lnu = dblarr(Nfilters, Nmodels)
 for i=0, (Nmodels-1) do begin
   fU[*, i] /= trap_int(U[*, i], fU[*, i], xrange=[Umin[i], Umax[i]], /pow)
   pow_LTIR[i] = gam[i]*trap_int(U[*, i], Lbol[*, i]*fU[*, i], xrange=[Umin[i], Umax[i]], /pow)
   for k=0, (Nfilters-1) do pow_mean_Lnu[k, i] = gam[i]*trap_int(U[*, i], mean_Lnu[*, k, i]*fU[*, i], xrange=[Umin[i], Umax[i]], /pow)
 endfor
 pow_mean_Lnu[where(pow_mean_Lnu lt 1.d-150,/null)] = 0.d0


 ;delta function at Umin
 del_LTIR=dblarr(Nmodels)
 del_mean_Lnu = dblarr(Nfilters, Nmodels)
 for i=0, (Nmodels-1) do begin
   del_LTIR[i] = (1.d0-gam[i])*interpol(Lbol[*, i], U[*, i], Umin[i])
   for k=0, (Nfilters-1) do del_mean_Lnu[k, i] = (1.d0-gam[i])*interpol(mean_Lnu[*, k, i], U[*, i], Umin[i])
 endfor
 del_mean_Lnu[where(del_mean_Lnu lt 1.d-150,/null)] = 0.d0

 ;combine delta function and power-law
 LTIR     = del_LTIR         + pow_LTIR
 mean_Lnu_dust = del_mean_Lnu + pow_mean_Lnu

 if Nmodels eq 1 and Nfilters eq 1 then return, mean_Lnu_dust[0] else $
   return, mean_Lnu_dust

end