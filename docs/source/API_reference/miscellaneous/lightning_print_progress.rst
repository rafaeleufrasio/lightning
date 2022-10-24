LIGHTNING_PRINT_PROGRESS
========================

Name
----
LIGHTNING_PRINT_PROGRESS

Purpose
-------
Prints a progress bar to the screen to track the speed and
progress of an iterative function in Lightning.

Calling Sequence
----------------
::

    lightning_print_progress, iteration, Niteration, t0 [, funcname = ]

Inputs
------
``iteration`` : int, float, or double scalar
    Number of current iteration.
``Niteration`` : int, float, or double scalar
    Total number of iterations.
``t0`` : int, float, or double scalar
    The overall starting system time :math:`[\rm{s}]`. Should be 
    the same for each iteration.

Optional Input
--------------
``funcname`` : string scalar
    Name of calling function that is printed in front of progress bar.
    (Default = ``'LIGHTNING_MCMC'`` )

Output
------
A progress bar that is printed to the screen, formatted as follows::

    FUNCNAME 22%|####      | 22000/100000 | xxx.x steps/s | xxxxx.x s remaining

The elements are, in order, the argument of ``funcname``, percentage comlete,
a progress bar, ``iteration/Niteration``, average speed, time remaining in seconds.

Notes
-----
Overhead on this function is about 17 microseconds per iteration.

Modification History
--------------------
- 2022/02/07: Created (Erik B. Monson)
- 2022/03/26: Documentation update (Erik B. Monson)
- 2022/06/27: Converted to procedure since no variable output (Keith Doore)
- 2022/06/27: Updated variable names to be more generic (Keith Doore)
- 2022/06/27: Updated from ``(iteration+2)`` to ``(iteration+1)`` (Keith Doore)
- 2022/06/28: Simplified number of chunks using percentage (Keith Doore)
- 2022/06/28: Removed remainder for double the speed (Keith Doore)
- 2022/08/18: Made it so that if ``iteration+1 == Niteration`` then end with new line (Keith Doore)
- 2022/10/24: Update to make changing values to take up a constant number of characters (Keith Doore)

