function gelman_rubin_stat, chain
;+
; Name
; ----
;   GELMAN_RUBIN_STAT
;
; Purpose
; -------
;   Computes the Gelman-Rubin statistic (:math:`\hat{R}`) of multiple Markov chains.
;   Uses the within-chain variance and between-chain variance to compare
;   convergence of the chains for each individual parameter. The use of 
;   degrees of freedom has been updated to the corrected version from 
;   Brooks & Gelman (1998).
;
; Calling Sequence
; ----------------
;   ::
;
;       r_hat = gelman_rubin_stat(chain)
;
; Input
; -----
;   ``chain`` : int, float, or double array(Nparam, Ntrials, Nparallel)
;       An ensemble of Markov chains.
;
; Output
; ------
;   ``r_hat`` : double array(Nparam)
;       The Gelman-Rubin statistic of the chains. The :math:`\sqrt{\hat{R}}`
;       should be less than or equal to 1.2 if convergence of the
;       chain was sufficiently reached. If a parameter is constant, the
;       corresponding Gelman-Rubin stat is ``NaN``.
;
; References
; ----------
;   - `Gelman, A., & Rubin, D. B. 1992, StaSc, 7, 457 <https://ui.adsabs.harvard.edu/abs/1992StaSc...7..457G/abstract>`_
;   - `Brooks, S. P., & Gelman, A. 1998, Journal of Computational and Graphical Statistics, 7, 434 <https://www.tandfonline.com/doi/abs/10.1080/10618600.1998.10474787>`_
;
; Modification History
; --------------------
;   - 2020/06/10: Created (Keith Doore)
;   - 2022/07/14: Updated documentation and error handling (Keith Doore)
;-
 Compile_opt idl2
 On_error,2

; Error handling
 if n_elements(chain) eq 0 then message, 'Variable is undefined: CHAIN.'
 if size(chain, /type) lt 2 or size(chain, /type) gt 5 then message, 'CHAIN must be of type int, float, or double.'
 if size(chain, /n_dim) lt 2 or size(chain, /n_dim) gt 3 then message, 'CHAIN must be a 2-D or 3-D array.'


 n_dim = size(chain, /n_dim)
 dim = size(chain, /dim)
 if n_dim eq 2 then Nparam = 1 else Nparam = dim[0]
 ; Ntrials is N in Gelman & Rubin (1992)
 if n_dim eq 2 then Ntrials = dim[0] else Ntrials = dim[1]
 ; Nparallel is M in Gelman & Rubin (1992)
 if n_dim eq 2 then Nparallel = dim[1] else Nparallel = dim[2]

 if Ntrials lt 10 then message, 'CHAIN must have more trials to compute Gelman-Rubin statistic.'


; Compute Gelman-Rubin statistic
  mean_i = mean(chain, dim=n_dim-1)
  var_i = variance(chain, dim=n_dim-1)
  xbar = mean(mean_i, dim=n_dim-1)

  B = Ntrials * variance(mean_i, dim=n_dim-1)
  W = mean(var_i, dim=n_dim-1)

  var_hat = W * (Ntrials-1.d) / Ntrials + B / Ntrials
  v_hat = var_hat + B / (Nparallel * Ntrials)
  
  covar1 = dblarr(Nparam)
  covar2 = dblarr(Nparam)
  for i=0, Nparam-1 do begin
    covar1[i] = CORRELATE(var_i[i, *], mean_i[i, *]^2.d, /covariance)
    covar2[i] = CORRELATE(var_i[i, *], mean_i[i, *]    , /covariance)
  endfor
  
  var_hat_v_hat = ((Ntrials-1.d) / Ntrials)^2.d / Nparallel * variance(var_i, dim=n_dim-1) + $
                  ((Nparallel+1.d) / (Nparallel * Ntrials))^2.d * 2.d/(Nparallel-1.d) * B^2.d + $
                  2.d*((Nparallel+1.d) * (Ntrials-1.d)) / (Nparallel * Ntrials^2.d) * $
                  Ntrials / Nparallel * (covar1 - 2.d * xbar * covar2)
                  
  df=2.d*v_hat / var_hat_v_hat
  
  r_hat = (df+3.d)/(df+1.d) * v_hat / W

  return, abs(r_hat)

end