TRAP_INT
========

Name
----
TRAP_INT

Purpose
-------
Integrates a tabulated function using the trapeziodal rule, with options
to increase accuracy in the case of functions having exponential
or power-law variation by recasting trapeziodal rule in linear-log or
log-log (x,y) space, respectively. A user can also specify the limits
of integration, and the function will be automatically extrapolated
if requested range exceeds the tabulated range of function.

Calling Sequence
----------------
::

   fint = trap_int(x, f [, xrange = , log_axes = , /exponential, $
                   /power_law, /vector, /delta_vector, /error_check])

Inputs
------
``x`` : int, float, or double array(Nf)
    Independent variable values.
``f`` : int, float, or double array(Nf)
    Function values at points ``x``.

Optional Inputs
---------------
``xrange`` : int, float, or double scalar or array(2)
    Specifies the limits of integration. If just one number is given, it is
    assumed to be the lower limit. If requested range exceeds tabulated range
    of function, the values are extrapolated. (Default = ``minmax(x)``)
``log_axes`` : int, float, or double scalar
    Alternative indicator of integration method:
    ``0`` = Linear-Linear space (Default), 
    ``1`` = Linear-Log (equivalent to using ``/exponential``), 
    ``2`` = Log-Log (equivalent to using ``/power_law``).
``exponential`` : flag
    If set, integrate (and interpolate) in linear-log space, which increases
    accuracy for exponentially varying functions. Takes precedence over 
    ``log_axes``. In this case, ``f`` array is thresholded > ``1.d-200``.
``power_law`` : flag
    If set, integrate (and interpolate) in log-log space, which increases
    accuracy for functions with power-law variation. Takes precedence over
    ``log_axes``. In this case, both ``x`` and ``f`` arrays are 
    thresholded > ``1.d-200``.
``vector`` : flag
    If set, returns a vector of cumulative definite integrals.
``delta_vector`` : flag
    If set, returns a vector of non-cumulative definite integrals
    over the intervals between the ``x`` values.
``error_check`` : flag
    If set, all inputs are checked for errors. Otherwise, all inputs are
    assumed to be of correct format.

Output
------
``fint`` : double scalar or array(Nf)
    Approximate integral(s). The size of ``fint`` is determined if one of
    the optional inputs, ``vector`` or ``delta_vector``, is set. If one is
    set, then ``fint`` is an array. Otherwise, it is a scalar.

Modification History
--------------------
- 1988/01/01: Created (Frank Varosi)
- 2020/02/18: Truncation value changed from ``1.e-37`` to ``1.d-200`` (Rafael Eufrasio)
- 2021/03/21: Standardized parameter names (Keith Doore)
- 2021/03/21: Updated documentation (Keith Doore)
- 2021/03/21: Corrected indexing to brackets from parentheses (Keith Doore)
- 2021/03/21: Added error handling (Keith Doore)
- 2022/04/12: Allowed for inputs to have degenerate dimensions (Keith Doore)
- 2022/04/12: Allowed integer inputs (Keith Doore)
- 2022/05/17: Added ``error_check`` keyword to do error handling (Keith Doore)

