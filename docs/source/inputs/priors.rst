.. _priors-label:

Prior Examples
==============

Lightning only has two options for analytical priors: uniform and truncated normal (Gaussian) priors.
Additionally, there is the option to input a prior of any shape in tabulated form.
Finally, while not actually a prior, we consider fixing a free parameter in a model a prior for implementation purposes.

For general cases, we recommend using using the default uniform or fixed priors specified in the configuration.
However, the choice of prior strongly depends on your research goals, which means different kinds of SEDs will
likely need different priors. Below, we give examples on how to input the different prior options into
the configuration using the file method, along with an example of how to
format the tabulated priors.

.. note::

    Analytical priors must be the same for all SEDs of a given run of Lightning. While this
    is typically not an issue for uniform priors, this can be a problem if you want unique
    normal priors for each SED in the batch. However, tabulated priors can be unique for each
    SED, and unique normal priors can be implemented this way by inputting them in tabulated form.
    Alternatively, you can run a unique configuration of Lightning for each unique normal prior.


File Configuration Method
-------------------------

There are three components specifying each prior when configuring Lightning with the file method. They
are the ``PRIOR``, which is the prior distribution type; ``PRIOR_ARG``, which is the prior distribution shape arguments;
and ``INITIALIZATION_RANGE``, which is the range to initialization the fitting algorithm.
Setting the prior distribution (``PRIOR``) for a parameter is simple: it is set to a string containing one of
the four prior distribution types: ``'fixed'``, ``'uniform'``, ``'normal'``, or ``'tabulated'``.
If the parameter has multiple parts (e.g., ``PSI``), each part can be a unique prior type, excluding
``'tabulated'``, which can only apply to an entire parameter. Here is and example for ``TAUV`` and ``PSI`` parameters:

.. code-block:: idl
    :emphasize-lines: 1,5

    TAUV:      {Prior: 'uniform'                                           ,$
                Prior_arg: [0.0d, 10.d]                                    ,$
                Initialization_range: [0.0d, 3.0d]                          }

    PSI:       {Prior:['uniform', 'normal', 'fixed', 'uniform', 'normal'] ,$
                Prior_arg: [[0.d0, 0.d0, 5.0d, 0.d0, 0.d0]                 ,$
                            [1.d3, 1.d3, 0.d0, 1.d3, 1.d3]                 ,$
                            [0.d0, 1.d1, 0.d0, 0.d0, 5.d0]                 ,$
                            [0.d0, 2.d0, 0.d0, 0.d0, 2.d0]]                ,$
                Initialization_range: [[0.0d, 5.0d, 0.0d, 0.0d, 1.0d]      ,$
                                       [1.d1, 15.d, 0.0d, 1.d1, 9.0d]]      }

As for the prior distribution shape arguments (``PRIOR_ARG``), these give the values
which define the shape of the prior distribution. In the case of a fixed prior, this is
the fixed value. For a uniform prior, this is the minimum and maximum bounds of the
distribution. For the normal prior, this is the minimum bound, maximum bound,
distribution peak, and distribution standard deviation. Finally, for the tabulated
prior, this is a string containing the path to the directory containing the user
defined :ref:`tabulated prior file <tabulated-prior-label>`. If the parameter has multiple
parts (e.g., ``PSI``) and the chosen priors types are different, then the array size used to define
the shape of the prior will use the largest number of shape arguments. For priors with
less shape arguments, they will need padded values (we recommend zeros) included
where the additional arguments of the more complex priors are required. The shape arguments
in our examples are highlighted below:

.. code-block:: idl
    :emphasize-lines: 2,6-9

    TAUV:      {Prior: 'uniform'                                           ,$
                Prior_arg: [0.0d, 10.d]                                    ,$
                Initialization_range: [0.0d, 3.0d]                          }

    PSI:       {Prior:['uniform', 'normal', 'fixed', 'uniform', 'normal'] ,$
                Prior_arg: [[0.d0, 0.d0, 5.0d, 0.d0, 0.d0]                 ,$
                            [1.d3, 1.d3, 0.d0, 1.d3, 1.d3]                 ,$
                            [0.d0, 1.d1, 0.d0, 0.d0, 5.d0]                 ,$
                            [0.d0, 2.d0, 0.d0, 0.d0, 2.d0]]                ,$
                Initialization_range: [[0.0d, 5.0d, 0.0d, 0.0d, 1.0d]      ,$
                                       [1.d1, 15.d, 0.0d, 1.d1, 9.0d]]      }

For the ``PSI`` parameter, the first column of prior shape arguments corresponds to the
first uniform prior, and the second column to the first normal prior and so on. Notice that
there are four rows, since the normal priors require four values. However, the uniform and
fixed priors only require the first two and first row(s), respectively, and the remaining row
values are padded with zeros.

Finally, for the initialization range (``INITIALIZATION_RANGE``),
this gives the minimum and maximum bounds to :ref:`randomly initialize the fitting algorithm
<random-initialize-label>`. This will be a two element array regardless of prior type.
The initialization ranges in our example are highlighted below:

.. code-block:: idl
    :emphasize-lines: 3,10,11

    TAUV:      {Prior: 'uniform'                                           ,$
                Prior_arg: [0.0d, 10.d]                                    ,$
                Initialization_range: [0.0d, 3.0d]                          }

    PSI:       {Prior:['uniform', 'normal', 'fixed', 'uniform', 'normal'] ,$
                Prior_arg: [[0.d0, 0.d0, 5.0d, 0.d0, 0.d0]                 ,$
                            [1.d3, 1.d3, 0.d0, 1.d3, 1.d3]                 ,$
                            [0.d0, 1.d1, 0.d0, 0.d0, 5.d0]                 ,$
                            [0.d0, 2.d0, 0.d0, 0.d0, 2.d0]]                ,$
                Initialization_range: [[0.0d, 5.0d, 0.0d, 0.0d, 1.0d]      ,$
                                       [1.d1, 15.d, 0.0d, 1.d1, 9.0d]]      }

Notice that the fixed prior on ``PSI`` has both values set to zero (third column) rather than
the fixed value. This is because Lightning will automatically ignore the initialization range
for fixed priors and will use the fixed value instead.

.. _tabulated-prior-label:

Tabulated Priors
----------------

As stated above, the prior argument for a tabulated prior is a string containing
the path to the directory containing the user defined tabulated prior file. This directory
must contain a tabulated prior file for each SED within the input data set, since all
SEDs must have the same prior distribution type. These tabulated
prior files must be FITS files with the tabulated prior data within a data table in the first
extension of the FITS file. The file must be named following the scheme of:
``tabulated_prior_<sed_id>.fits``. All tabulated priors for this SED must be within the single
data table and have the following vector columns for each parameter: ``<PARAM-NAME>_values`` and
``<PARAM-NAME>_pdf``, where ``<PARAM-NAME>`` should be replaced with the actual parameter name
(e.g., ``PSI`` or ``TAUV``). ``<PARAM-NAME>_values`` are the gridded values of the parameter
corresponding to ``<PARAM-NAME>_pdf``, which is the probability of the prior at each gridded value.
The only restrictions on these columns are:

1) ``<PARAM-NAME>_values`` must be in ascending order,
2) ``<PARAM-NAME>_values`` must all be unique,
3) ``<PARAM-NAME>_pdf`` must contain non-negative values, and
4) ``<PARAM-NAME>_pdf`` must be area normalized (i.e., total area sums to 1).

.. highlight:: idl

Otherwise, any distribution shape is allowed. To help clarify this,
we give an example of two overlapping Gaussians for the prior of the parameter ``TAUV``::

    ; Set tauV range from 0 to 3
    tauv_values = [0:3:0.001d]
    ; First Gaussian of mean 1 and width 0.2
    first_gaussian = exp(-0.5d*((tauv_values - 1)/0.2d)^2)
    ; Second Gaussian of mean 2 and width 0.5
    second_gaussian = exp(-0.5d*((tauv_values - 2)/0.5d)^2)
    tauv_pdf = first_gaussian + second_gaussian
    ; Normalize area to 1
    tauv_pdf /= int_tabulated(tauv_values, tauv_pdf, /double)
    data_table = {tauv_values: tauv_values, tauv_pdf: tauv_pdf}
