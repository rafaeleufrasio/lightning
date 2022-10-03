function autocorr_func, x
; Name
; ----
;   AUTOCORR_FUNC
;
; Purpose
; -------
;   Computes the autocorrelation function of a time series x, computed with a
;   fast Fourier transform. This is quite noisy, and can potentially be quite
;   wrong when tested with analytic functions with few numbers of samples.
;   This noise is mediated in practice by averaging over the many and longer
;   Markov chains.
;
; Calling Sequence
; ----------------
;   ::
;
;       acf = autocorr_func(x)
;
; Input
; -----
;   ``x`` : int, float, or double array(Ntrials)
;       A time series.
;
; Output
; ------
;   ``acf`` : array(Nparam)
;       The autocorrelation function of the time series. If the time series is
;       constant, the corresponding autocorrelation function is NaN.
;
; Modificiation History
; ---------------------
;   - 2022/02/03: Created (E.B. Monson)
;   - 2022/07/14: Added documentation and improved error handling (Keith Doore)
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(x) eq 0 then message, 'Variable is undefined: X.'
 if size(x, /type) lt 2 or size(x, /type) gt 5 then message, 'X must be of type int, float, or double.'
 if size(x, /n_dim) ne 1 then message, 'X must be a 1-D array.'


 n = n_elements(x)
 ; Ignore constant dimensions.
 ; Sets the autocorrelation function to NaN, since the time series never changes.
 if (variance(x, /double) eq 0) then return, replicate(!values.D_NaN, n)


 ; FFT is far more efficient on arrays with length 2^n for integer n.
 ;   So, we pad our input out to the next power of 2 with zeros.
 n2 = long(1)
 while (n2 lt n) do n2 = ishft(n2, 1)
 f = dblarr(n2)
 f[0:n-1] = x - mean(x, /double)

 ; Compute the autocorrelation function with an FFT
 ff = fft(f)
 acf = real_part(fft(ff * conj(ff), /inverse))
 acf = acf[0:n-1]
 acf = acf / acf[0]

 return, acf

end


function autocorr_time, chain, C_step=C_step
;+
; Name
; ----
;   AUTOCORR_TIME
;
; Purpose
; -------
;   Computes the integrated autocorrelation time of an ensemble of Markov chains,
;   a measure of how many steps it takes for the chain to forget where it started.
;   In other words, the autocorrelation time is how many samples are needed, on
;   average, to get an independent sample. Integrating the autocorrelation time
;   over the entire chain is noisy. So, instead, it is integrated to the smallest
;   index :math:`M` such that :math:`M > C_{\rm step} \tau`.
;
; Calling Sequence
; ----------------
;   ::
;
;       tau = autocorr_time(chain [, C_step =])
;
; Input
; -----
;   ``chain`` : int, float, or double array(Nparam, Ntrials, Nparallel)
;       An ensemble of Markov chains.
;
; Optional Input
; --------------
;   ``C_step`` : int, float, or double scalar
;       Defines how many trials of the chain are used to calculate ``tau``, where
;       ``tau`` is integrated to the smallest index ``M`` such that ``M > C_step * tau``.
;       (Default = ``5``)
;
; Output
; ------
;   ``tau`` : array(Nparam)
;       Integrated autocorrelation time for each parameter, averaged over the number
;       of chains. If a parameter is constant, the corresponding autocorrelation
;       time is ``NaN``.
;
; Notes
; -----
;   The length of the chain divided by the autocorrelation time is a measurement of how many
;   independent MCMC samples were generated. If that number is smaller than, say, a few thousand
;   for any parameter, it may be necessary to restrict the model or run a longer chain. The
;   autocorrelation time can also be viewed as giving us the scales for burn-in and thinning.
;   The examples in the `emcee <https://emcee.readthedocs.io/en/stable/tutorials/autocorr/#autocorr>`_
;   documentation suggest discarding samples up to a few times the autocorrelation time and
;   thinning by one-half to one times the autocorrelation time.
;
;   Ported from the autocorrelation analysis in the python package `emcee`. See:
;   - `<https://emcee.readthedocs.io/en/stable/tutorials/autocorr/#autocorr>`_
;   - `<https://github.com/dfm/emcee>`_
;
; Modificiation History
; ---------------------
;   - 2022/02/03: Created (E.B. Monson)
;   - 2022/06/01: Allow ``tolerance`` to be set to 0 for silent operation. (E.B. Monson)
;   - 2022/07/14: Removed ``tolerance`` input and moved outside of function, since we just want the autocorrelation time (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Error handling
 if n_elements(chain) eq 0 then message, 'Variable is undefined: CHAIN.'
 if size(chain, /type) lt 2 or size(chain, /type) gt 5 then message, 'CHAIN must be of type int, float, or double.'
 if size(chain, /n_dim) lt 2 or size(chain, /n_dim) gt 3 then message, 'CHAIN must be a 2-D or 3-D array.'

 if n_elements(C_step) ne 0 then begin
   if size(C_step, /type) lt 2 or size(C_step, /type) gt 5 then message, 'C_STEP must be of type int, float, or double.'
   if size(C_step, /n_dim) ne 0 then message, 'C_STEP must be a scalar.'
   if C_step le 0 then message, 'C_STEP must be a positive value.'
 endif else C_step = 5.d


; Compute autocorrelation time
 dim = size(chain, /dimensions)
 Nparam = dim[0]
 Ntrials = dim[1]
 if n_elements(dim) eq 2 then Nparallel = 1 else Nparallel = dim[2]

 tau_est = dblarr(Nparam)

 for p=0, Nparam-1 do begin
   ac = dblarr(Ntrials)
   for c=0, Nparallel-1 do begin
       ac = ac + autocorr_func(reform(chain[p, *, c]))
   endfor
   ; Average the autocorrelation function across the ensemble of chains
   ac = ac / double(Nparallel)
   ; Subtract one rather than add, due to rho(0) being included in the sum and rho(0) = 1
   ;   tau = 1 + 2 * sum(rho(tau), from=1, to=N) = 2 * sum(rho(tau), from=0, to=N) - 1
   tau_series = 2.0 * total(ac, /cumulative) - 1

   ; Get the window value M
   idcs = where(findgen(n_elements(tau_series)) gt (C_step * tau_series), count)
   if (count gt 0) then m = min(idcs) else m = n_elements(tau_series) - 1

   tau_est[p] = tau_series[m]
 endfor

 return, tau_est

end
