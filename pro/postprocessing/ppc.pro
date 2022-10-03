function ppc, Nrep, Lobs, Lunc, Lpredict, lnprob_chain, $
              counts_obs=counts_obs, counts_unc=counts_unc, counts_predict=counts_predict, $
              chisqr_obs=chisqr_obs, chisqr_rep=chisqr_rep, model_unc=model_unc
;+
; Name
; ----
;   PPC
;
; Purpose
; -------
;   Computes a Posterior Predictive Check (PPC) on a MCMC Lightning output.
;   Uses the methods described in Rubin (1984) and Gelman et al. (1996).
;
; Calling Sequence
; ----------------
;   ::
;
;       pvalue = ppc(Nrep, Lobs, Lunc, Lpredict, lnprob_chain [, $
;                    counts_obs = , counts_unc = , counts_predict = , $
;                    chisqr_obs=chisqr_obs, chisqr_rep=chisqr_rep])
;
; Inputs
; ------
;   ``Nrep`` : int, float or double scalar
;       The number of times to replicate the predicted data.
;   ``Lobs`` : int, float or double array(Nfilters)
;       The observed luminosities [arbitrary units].
;   ``Lunc`` : int, float or double array(Nfilters)
;       The uncertainties on the observed luminosities in the same units as ``Lobs``.
;   ``Lpredict`` : int, float, or double array(Nfilters, Nchain)
;       The model luminosities predicted by the model given a set of parameters in
;       same units as ``Lobs``.
;   ``lnprob_chain`` : int, float, or double array(Nchain)
;       The log posterior probability values of each SED predicted
;       by the model given a set of parameters. Used to select the replicated data.
;
; Optional Inputs
; ---------------
;   ``counts_obs`` : int, float, or double array(Nxray)
;       The observed X-ray counts. If specified, ``counts_unc`` and ``counts_predict``
;       must also be given.
;   ``counts_unc`` : int, float, or double array(Nxray)
;       The uncertainties on the observed X-ray counts. If specified, ``counts_obs`` and ``counts_predict``
;       must also be given.
;   ``counts_predict``: int, float, or double array(Nxray, Nchain)
;       The X-ray counts predicted by the model. If specified, ``counts_obs`` and ``counts_unc``
;       must also be given.
;   ``model_unc`` : int, float, or double scalar
;       The fractional model uncertainty to use in all bands. (Default = ``0.d0``)
;
; Output
; ------
;   ``pvalue`` : double scalar
;       The p-value associated with the PPC. Determined as the fraction
;       of ``chisqr_rep`` that is greater than ``chisqr_obs``.
;
; Optional Outputs
; ----------------
;   ``chisqr_obs`` : double array(Nrep)
;       The resulting :math:`\chi^2` values by comparing the predicted data with the
;       observational data.
;   ``chisqr_rep`` : double array(Nrep)
;       The resulting :math:`\chi^2` values by comparing the predicted data with the
;       replicated data.
;
; Notes
; -----
;   - If the X-ray model was fit using fluxes rather than counts, the X-ray luminosities can
;     simply be appended to ``Lobs``, ``Lunc``, and ``Lpredict``.
;
; References
; ----------
;   - `Rubin D. B., 1984, Ann. Stat., 12, 1151 <https://www.jstor.org/stable/2240995>`_
;   - `Gelman A., Meng X.-L., Stern H., 1996, Stat. Sin., 6, 733 <https://www.jstor.org/stable/24306036>`_
;
; Modification History
; --------------------
;   - 2022/03/15: Created (Keith Doore)
;   - 2022/04/18: Added documentation (Keith Doore)
;   - 2022/04/18: Added error handling (Keith Doore)
;   - 2022/07/21: Added ability to handle ``NaNs`` in ``lnprob_chain`` (Keith Doore)
;   - 2022/07/27: Fixed bug with missing keyword ``nrand`` in ``mrandomn`` (Keith Doore)
;   - 2022/07/27: Fixed bug where ``value_locate`` maybe selecting ``-1`` when we want ``0``,
;     force ``-1`` to be ``0`` (Keith Doore)
;   - 2022/09/01: Added X-ray emission (Erik B. Monson)
;   - 2022/09/15: Changed from ``chi2_chain`` input to ``lnprob_chain`` (Keith Doore)
;   - 2022/09/30: Added ``model_unc`` input to include in uncertainties (Keith Doore)
;-
 On_error,2
 compile_opt idl2

; Error handling
 if n_elements(Nrep) eq 0 then message, 'Variable is undefined: Nrep.'
 if size(Nrep, /type) lt 2 or size(Nrep,/type) gt 5 then message, 'Nrep must be of type int, float, or double.'
 if size(Nrep, /n_dim) ne 0 then message, 'Nrep must be a scalar.'
 if min(Nrep) le 0 then message,'Nrep must be a positive value.'

 if n_elements(Lobs) eq 0 then message, 'Variable is undefined: LOBS.'
 if size(Lobs, /type) lt 2 or size(Lobs,/type) gt 5 then message, 'LOBS must be of type int, float, or double.'
 if size(reform(Lobs), /n_dim) ne 1 then message, 'LOBS must be a 1-D array.'

 if n_elements(Lunc) eq 0 then message, 'Variable is undefined: LUNC.'
 if size(Lunc, /type) lt 2 or size(Lunc,/type) gt 5 then message, 'LUNC must be of type int, float, or double.'
 if size(reform(Lunc), /n_dim) ne 1 then message, 'LUNC must be a 1-D array.'
 if n_elements(Lunc) ne n_elements(Lobs) then message, 'LUNC must have the same number of elements as LOBS.'

 if n_elements(Lpredict) eq 0 then message, 'Variable is undefined: LPREDICT.'
 if size(Lpredict, /type) lt 2 or size(Lpredict,/type) gt 5 then message, 'LPREDICT must be of type int, float, or double.'
 if size(reform(Lpredict), /n_dim) ne 2 then message, 'LPREDICT must be a 2-D array.'
 if (size(reform(Lpredict), /dim))[0] ne n_elements(Lobs) then $
    message, 'LPREDICT must have a first dimension with the same number of elements as LOBS.'

 if n_elements(lnprob_chain) eq 0 then message, 'Variable is undefined: LNPROB_CHAIN.'
 if size(lnprob_chain, /type) lt 2 or size(lnprob_chain,/type) gt 5 then $
   message, 'LNPROB_CHAIN must be of type int, float, or double.'
 if size(reform(lnprob_chain), /n_dim) ne 1 then message, 'LNPROB_CHAIN must be a 1-D array.'
 if n_elements(lnprob_chain) ne (size(reform(Lpredict), /dim))[1] then $
    message, 'LNPROB_CHAIN must have the same number of elements as the second dimension of LPREDICT.'

; X-ray input handling
 xray_emission = 0
 if n_elements(counts_obs) ne 0 then begin
   if size(counts_obs, /type) lt 2 or size(counts_obs,/type) gt 5 then $
            message, 'COUNTS_OBS must be of type int, float, or double.'
   if size(reform(counts_obs), /n_dim) ne 1 then message, 'COUNTS_OBS must be a scalar or 1-D array.'
   if min(counts_obs) lt 0 then message, 'COUNTS_OBS must only contain non-negative values.'
   xray_emission += 1
   Nxray = n_elements(counts_obs)
 endif

 if n_elements(counts_unc) ne 0 then begin
   if size(counts_unc, /type) lt 2 or size(counts_unc,/type) gt 5 then $
            message, 'COUNTS_UNC must be of type int, float, or double.'
   if size(reform(counts_unc), /n_dim) ne 1 then message, 'COUNTS_UNC must be a scalar or 1-D array.'
   if min(counts_unc) lt 0 then message, 'COUNTS_UNC must only contain non-negative values.'
   xray_emission += 1
 endif

 if n_elements(counts_predict) ne 0 then begin
   if size(counts_predict, /type) lt 2 or size(counts_predict,/type) gt 5 then $
            message, 'COUNTS_PREDICT must be of type int, float, or double.'
   if size(reform(counts_predict), /n_dim) ne 2 then message, 'COUNTS_PREDICT must be a 2-D array.'
   if min(counts_predict) lt 0 then message, 'COUNTS_PREDICT must only contain non-negative values.'
   if (size(reform(counts_predict), /dim))[1] ne n_elements(lnprob_chain) then $
       message, 'The second dimension of COUNTS_PREDICT must have the same number of elements as LNPROB_CHAIN.'
   xray_emission += 1
 endif

 if keyword_set(xray_emission) then if xray_emission ne 3 then $
     message, 'COUNTS_OBS, COUNTS_UNC, and the first dimension of COUNTS_PREDICT must all have the same length.'

 if n_elements(model_unc) ne 0 then begin
   if size(model_unc, /type) lt 2 or size(model_unc, /type) gt 5 then $
     message, 'MODEL_UNC must be of type int, float, or double.'
   if size(model_unc, /n_dim) ne 0 then $
     message, 'MODEL_UNC must be a scalar.'
   if model_unc lt 0 then $
     message, 'MODEL_UNC must be a be a non-negative value.'
 endif else model_unc = 0.d0


; Select the indices to use from the lnprob_chain, with replacement
 pdf = exp(reform(lnprob_chain))
 sorted = sort(pdf)
 cdf = total(pdf[sorted], /cum)/max(total(pdf[sorted], /cum))
 finite_cdf = cdf[where(finite(cdf), /null)]
 ; If all elements in cdf are NaN, then end calculation and return pvalue of 0
 if n_elements(finite_cdf) eq 0 then return, 0.d0
 sorted_index = value_locate(finite_cdf, randomu(seed, Nrep))
 ; value_locate will output -1 index if a random value is smaller than the smallest value in finite_cdf
 ;   In that case we want 0 instead.
 sorted_index = sorted_index > 0
 index = sorted[sorted_index]
 Lpredict_rep = (reform(Lpredict))[*, index]
 if keyword_set(xray_emission) then counts_predict_rep = (reform(counts_predict))[*, index]

; Compute the replicated and observed chisqr
 ; Priors are not required in chisqr computation as they will be the same
 ;  in both the replicated and actual data. Therefore, they would just
 ;  be an additive factor in both parts, resulting in no change.
 Nfilters = n_elements(Lobs)
 Lrep = dblarr(Nfilters, Nrep)

 ; Find zeros that may be contained in Lunc
 nonzero = where(Lunc ne 0, comp=zero, /null)

 ; Determine total uncertainty, including model uncertainty
 Lunc_total = sqrt(rebin(reform(Lunc^2), Nfilters, Nrep) + (model_unc * Lpredict_rep)^2)
 Lunc_total[zero, *] = 0.d0

 ; Generate replicated data from predicted data perturbed by the uncertainties
 ;Lrep[nonzero, *] = transpose(mrandomn(seed, diag_matrix(reform(Lunc[nonzero]^2)), nrand=Nrep)) + Lpredict_rep[nonzero, *]
 for i=0, Nrep-1 do $
   Lrep[nonzero, i] = reform(mrandomn(seed, diag_matrix(reform(Lunc_total[nonzero, i]^2)))) + Lpredict_rep[nonzero, i]

 ; Compute chisqrs
 chisqr_obs = total(((rebin(reform(Lobs), Nfilters, Nrep) - Lpredict_rep)/Lunc_total)^2, 1, /NaN)
 chisqr_rep = total(((Lrep - Lpredict_rep)/Lunc_total)^2, 1, /NaN)

 if keyword_set(xray_emission) then begin
    counts_rep = poidev(counts_predict_rep)
    xray_sigma = rebin(reform(counts_unc), Nxray, Nrep)
    xray_chisqr_obs = total((rebin(reform(counts_obs), Nxray, Nrep) - counts_predict_rep)^2 / xray_sigma^2, 1, /NaN)
    xray_chisqr_rep = total((counts_predict_rep - counts_rep)^2 / xray_sigma^2, 1, /NaN)
    chisqr_obs += xray_chisqr_obs
    chisqr_rep += xray_chisqr_rep
 endif

 pvalue = n_elements(where(chisqr_rep gt chisqr_obs, /null))/double(Nrep)

 return, pvalue

end
