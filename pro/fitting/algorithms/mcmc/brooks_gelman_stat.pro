function brooks_gelman_stat, chain
;+
; Name
; ----
;   BROOKS_GELMAN_STAT
;
; Purpose
; -------
;   Computes the Brooks-Gelman statistic (i.e., the multivariate
;   Gelman-Rubin statistic as given in Brooks & Gelman 1998) of 
;   multiple Markov chains. Uses the within-chain variance and 
;   between-chain variance for all parameters to compare convergence
;   of the chains.
;
; Calling Sequence
; ----------------
;   ::
;
;       r_hat = brooks_gelman_stat(chain)
;
; Input
; -----
;   ``chain`` : int, float, or double array(Nparam, Ntrials, Nparallel)
;       An ensemble of Markov chains.
;
; Output
; ------
;   ``r_hat`` : double scalar
;       The Brooks-Gelman statistic of the chains. The :math:`\sqrt{\hat{R}}`
;       should be less than or equal to 1.2 if convergence of the
;       chain was sufficiently reached.
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
 if size(chain, /n_dim) ne 3 then message, 'CHAIN must be a 3-D array.'


 n_dim = size(chain, /n_dim)
 dim = size(chain, /dim)
 ; Nparam is P in Brooks & Gelman (1998)
 Nparam = dim[0]
 ; Ntrials is N in Brooks & Gelman (1998)
 Ntrials = dim[1]
 ; Nparallel is M in Brooks & Gelman (1998)
 Nparallel = dim[2]

 if Ntrials lt 10 then message, 'CHAIN must have more trials to compute Gelman-Rubin statistic.'


; Compute Gelman-Rubin statistic
 psi_j_t = chain
 psi_j_dot = mean(chain, dim=2)
 psi_dot_dot = mean(psi_j_dot, dim=2)


 w_sum = 0.d0
 bn_sum = 0.d0
 for j = 0, Nparallel-1 do begin
   for t = 0, Ntrials-1 do begin
     w_sum += (psi_j_t[*, t, j] - psi_j_dot[*, j]) ## (psi_j_t[*, t, j] - psi_j_dot[*, j])
   endfor
   bn_sum += (psi_j_dot[*, j] - psi_dot_dot) ## (psi_j_dot[*, j] - psi_dot_dot)
 endfor

 W = 1.d/(Nparallel * (Ntrials-1.d)) * w_sum
 Bn = 1.d/(Nparallel-1.d) * bn_sum

 temp = invert(W) ## Bn
 ; Confirm matrix is symmetric, as it may not be from machine precision.
 temp = (temp+transpose(temp))/2.d
 lam1 = max(eigenql(temp))

 r_hat = (Ntrials-1.d) / Ntrials + (Nparallel+1.d) / Nparallel * lam1

 return, r_hat

end