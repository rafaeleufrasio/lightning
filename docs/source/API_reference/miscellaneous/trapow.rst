TRAPOW
======

Name
----
TRAPOW

Purpose
-------
Integrates a function having power-law variation given values on an
increasing grid of x-points, using the trapeziodal rule in log-log space.

Calling Sequence
----------------
::

    fint = trapow(f, x [, /vector, /delta_vector, /error_check])

Inputs
------
``f`` : int, float, or double array(Nf)
    Function values at points ``x``, must be all positive.
``x`` : int, float, or double array(Nf)
    Independent variable values, must be all positive.

Optional Inputs
---------------
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
- 1995/01/01: Created (Frank Varosi)
- 2001/03/07: ``i`` was changed to be a LONG integer in the loop (Eli Dwek)
- 2021/03/21: Standardized parameter names (Keith Doore)
- 2021/03/21: Updated documentation (Keith Doore)
- 2021/03/21: Corrected indexing to brackets from parentheses (Keith Doore)
- 2021/03/21: Added error handling (Keith Doore)
- 2022/04/08: Allowed for inputs to have degenerate dimensions (Keith Doore)
- 2022/04/08: Allowed integer inputs (Keith Doore)
- 2022/05/17: Added ``error_check`` keyword to do error handling (Keith Doore)

