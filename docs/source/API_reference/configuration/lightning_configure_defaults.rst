LIGHTNING_CONFIGURE_DEFAULTS
============================

Name
----
LIGHTNING_CONFIGURE

Purpose
-------
Generates the Lightning configuration structure. This structure
contains all of the possible selections users can make when
configuring Lightning (e.g., models, fitting algorithms, etc.).
The user can change the values for each tag to their desired values.
Small guiding comments are included. These comments give basic details
and indicate the expected values and type. Further details can be
found at :ref:`configure-setting-label`.

Calling Sequence
----------------
::

    config = lightning_configure()

Output
------
``config`` : structure
    A Lightning configuration structure.

Notes
-----
- Parameters that take a flag (value of ``0`` or ``1``) are set if the value is ``1``
  and not set if the value is ``0``.
- Options that take arrays are indicated by brackets (i.e., ``[]``). These options
  can have multiple values. Options without brackets must contain a single value.
  Removing brackets may result in errors when the configuration structure is
  checked for errors.
- All free parameters of the models have an associated prior structure, where
  the distribution type, prior distribution shape arguments, and initialization
  range are given. The prior distribution type options are: ``'fixed'``, ``'uniform'``,
  ``'normal'``, and ``'tabulated'``. The number of values in the distribution
  shape argument array depend on the chosen distribution type.

  - ``'fixed'``: takes a single value (``Narg = 1``), the value at which to fix the
    parameter. (If the prior is ``'fixed'``, then the initialization range is ignored,
    since a fixed value does not need initialization.)
  - ``'uniform'``: takes two values (``Narg = 2``), the minimum and maximum bounds of the
    distribution in that order.
  - ``'normal'``: takes four values (``Narg = 4``), the minimum bound, maximum bound,
    distribution peak, and distribution standard deviation in that order.
  - ``'tabulated'``: takes one value (``Narg = 1``), a string containing the path to the
    directory containing the user tabulated prior file (see :ref:`tabulated-prior-label`
    for more details).

  The initialization range always has two values per prior,
  indicating the minimum and maximum bounds for the random initialization of
  the fitting algorithm.

