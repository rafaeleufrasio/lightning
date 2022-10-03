PERCENTILE
==========

Name
----
PERCENTILE

Purpose
-------
Determines the percentiles values for an input array by sorting the
input array and returning an interpolated percentile value.
It assumes the lowest value of the array to be the 0th percentile
and the highest to be the 100th percentile.

Calling Sequence
----------------
::

    percentile_value = percentile(array, percentile)

Inputs
------
``array`` : int, float, or double array(N)
    The data that will have its percentile values determined.
``percentile`` : int, float, or double scalar or array(M, ...)
    The percentiles in decimal form (i.e., 50th percentile = ``0.5``).

Output
------
``percentile_value`` : double scalar or array(M, ...)
    The percentile values corresponding to each input percentile.

Modification History
--------------------
- 2016/05/01: Created (Rafael T. Eufrasio)
- 2022/03/15: Added error handling (Keith Doore)
- 2022/03/15: Updated documentation (Keith Doore)
- 2022/04/07: Updated error handling to allow for degenerate dimensions (Keith Doore)
- 2022/06/24: Replaced ``.length`` with ``n_elements`` (Keith Doore)

