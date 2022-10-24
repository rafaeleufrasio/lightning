pro lightning_print_progress, iteration, Niteration, t0, funcname=funcname
;+
; Name
; ----
;   LIGHTNING_PRINT_PROGRESS
;
; Purpose
; -------
;   Prints a progress bar to the screen to track the speed and
;   progress of an iterative function in Lightning.
;
; Calling Sequence
; ----------------
;   ::
;
;       lightning_print_progress, iteration, Niteration, t0 [, funcname = ]
;
; Inputs
; ------
;   ``iteration`` : int, float, or double scalar
;       Number of current iteration.
;   ``Niteration`` : int, float, or double scalar
;       Total number of iterations.
;   ``t0`` : int, float, or double scalar
;       The overall starting system time :math:`[\rm{s}]`. Should be 
;       the same for each iteration.
;
; Optional Input
; --------------
;   ``funcname`` : string scalar
;       Name of calling function that is printed in front of progress bar.
;       (Default = ``'LIGHTNING_MCMC'`` )
;
; Output
; ------
;   A progress bar that is printed to the screen, formatted as follows::
;
;       FUNCNAME 22%|####      | 22000/100000 | xxx.x steps/s | xxxxx.x s remaining
;
;   The elements are, in order, the argument of ``funcname``, percentage comlete,
;   a progress bar, ``iteration/Niteration``, average speed, time remaining in seconds.
;
; Notes
; -----
;   Overhead on this function is about 17 microseconds per iteration.
;
; Modification History
; --------------------
;   - 2022/02/07: Created (Erik B. Monson)
;   - 2022/03/26: Documentation update (Erik B. Monson)
;   - 2022/06/27: Converted to procedure since no variable output (Keith Doore)
;   - 2022/06/27: Updated variable names to be more generic (Keith Doore)
;   - 2022/06/27: Updated from ``(iteration+2)`` to ``(iteration+1)`` (Keith Doore)
;   - 2022/06/28: Simplified number of chunks using percentage (Keith Doore)
;   - 2022/06/28: Removed remainder for double the speed (Keith Doore)
;   - 2022/08/18: Made it so that if ``iteration+1 == Niteration`` then end with new line (Keith Doore)
;   - 2022/10/24: Update to make changing values to take up a constant number of characters (Keith Doore)
;-
 On_error,2
 compile_opt idl2

; We do not perform error handling on this procedure due to its simplicity and doing
;   so would decreased speed.

 if ~keyword_set(funcname) then funcname = 'LIGHTNING_MCMC'

 tnow = systime(/sec)
 t_elapsed =  tnow - t0
 avg_speed = (iteration + 1) / t_elapsed
 t_remaining = (Niteration - (iteration + 1)) / avg_speed

 percent = floor(100. * (iteration + 1) / Niteration)

 bar = ''
 N_chunks = 10.
 chunks_done = floor(percent / 100. * N_chunks)
 ;remainder = percent mod (100./N_chunks)
 for chunk=0, chunks_done-1 do bar = bar + '#' 
 ;bar = bar + string(remainder, f='(I0)')
 for chunk=0, N_chunks-chunks_done do bar = bar + ' '

 if iteration + 1 eq Niteration then begin
   print, string(13b), funcname, percent, bar, (iteration + 1), Niteration, avg_speed, t_remaining, $
          format='(A,A,": ",(I3),"%|",(A10),"| ",(I7),"/",(I7)," | ",(F5.1)," steps/s | ",(F9.1)," s remaining")'
 endif else begin
   print, string(13b), funcname, percent, bar, (iteration + 1), Niteration, avg_speed, t_remaining, $
          format='(A,A,": ",(I3),"%|",(A10),"| ",(I7),"/",(I7)," | ",(F5.1)," steps/s | ",(F9.1)," s remaining",$)'
 endelse

end