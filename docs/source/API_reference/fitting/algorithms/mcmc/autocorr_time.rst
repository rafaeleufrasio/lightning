AUTOCORR_TIME
=============

Name
----
AUTOCORR_TIME

Purpose
-------
Computes the integrated autocorrelation time of an ensemble of Markov chains,
a measure of how many steps it takes for the chain to forget where it started.
In other words, the autocorrelation time is how many samples are needed, on
average, to get an independent sample. Integrating the autocorrelation time
over the entire chain is noisy. So, instead, it is integrated to the smallest
index :math:`M` such that :math:`M > C_{\rm step} \tau`.

Calling Sequence
----------------
::

    tau = autocorr_time(chain [, C_step =])

Input
-----
``chain`` : int, float, or double array(Nparam, Ntrials, Nparallel)
    An ensemble of Markov chains.

Optional Input
--------------
``C_step`` : int, float, or double scalar
    Defines how many trials of the chain are used to calculate ``tau``, where
    ``tau`` is integrated to the smallest index ``M`` such that ``M > C_step * tau``.
    (Default = ``5``)

Output
------
``tau`` : array(Nparam)
    Integrated autocorrelation time for each parameter, averaged over the number
    of chains. If a parameter is constant, the corresponding autocorrelation
    time is ``NaN``.

Notes
-----
The length of the chain divided by the autocorrelation time is a measurement of how many
independent MCMC samples were generated. If that number is smaller than, say, a few thousand
for any parameter, it may be necessary to restrict the model or run a longer chain. The
autocorrelation time can also be viewed as giving us the scales for burn-in and thinning.
The examples in the `emcee <https://emcee.readthedocs.io/en/stable/tutorials/autocorr/#autocorr>`_
documentation suggest discarding samples up to a few times the autocorrelation time and
thinning by one-half to one times the autocorrelation time.

Ported from the autocorrelation analysis in the python package `emcee`. See:
- `<https://emcee.readthedocs.io/en/stable/tutorials/autocorr/#autocorr>`_
- `<https://github.com/dfm/emcee>`_

Modificiation History
---------------------
- 2022/02/03: Created (E.B. Monson)
- 2022/06/01: Allow ``tolerance`` to be set to 0 for silent operation. (E.B. Monson)
- 2022/07/14: Removed ``tolerance`` input and moved outside of function, since we just want the autocorrelation time (Keith Doore)

