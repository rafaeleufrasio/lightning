Getting Started
===============

.. toctree::
    :hidden:

    inputs/input_format.rst
    inputs/configuration.rst
    inputs/models.rst
    inputs/fitting_algorithm.rst
    inputs/priors.rst


The first step in getting started with Lightning requires you to format your :ref:`input-label`
such that they can be understood by Lightning. In the next step, :ref:`configure-label`, you select
a model, fitting algorithm, and more. Finally, once the input and configuration steps are complete,
you can start :ref:`run-lightning-label`.


.. _input-label:

Input Data
----------

Lightning requires you to input your data as either an ASCII or FITS data table. The data must contain:

1) the fluxes for the SED(s) in terms of :math:`F_\nu` and units of :math:`\rm Jy`,
2) the uncertainties on the fluxes in terms of :math:`F_\nu` and units of :math:`\rm Jy`,
3) the names (labels) of the filters corresponding to each flux (a list of filter labels understood
   by Lightning can be found :ref:`here <filter_table>`),
4) a distance indicator, either a luminosity distance (in units of :math:`{\rm Mpc}`) or redshift.

Additional inputs, such as an SED ID, can be included in the tables. See :ref:`input-formats-label`
for the full details on the input data and file formats.


.. _configure-label:

Configuring Lightning
---------------------

Before running, Lightning requires you to configure it (e.g., select a model, priors, fitting algorithm, etc.).
We have included two different methods for configuring Lightning: the
first is more complex, and requires you to edit a file directly, while the other gives interactive
prompts to the terminal for you to answer. Both methods are described below, and further details on
each configuration setting are discussed in the :ref:`configure-setting-label`.

.. note::

    We recommend that first time users to Lightning (and those newer to SED fitting)
    use the interactive method to get a better understanding of Lightning's configuration settings
    and minimize frustration.


.. _file-method-label:

File Configuration
^^^^^^^^^^^^^^^^^^

To use the default file method, you will first need to copy the ``lightning_configure.pro`` file
within the top level of the Lightning installation into your input data directory without changing
its name. You can then open the copied file and edit the values associated with each structure
tag to your desired configuration. The values already given for each tag indicate the default
Lightning configuration. Guiding comments are given in the notes of the header and
for each configuration setting above the corresponding tag. These comments give
basic details and indicate the expected values and type. Further details are discussed in
:ref:`configure-setting-label`.

.. note::

    We recommend making a new copy of ``lightning_configure.pro`` for each project or new set of fits.
    Additionally, if you want to check your edited file for errors before running Lightning (which will
    check the file for errors) check out the ``lightning_configure_check.pro`` function in the
    :ref:`api-label`.


Interactive Configuration
^^^^^^^^^^^^^^^^^^^^^^^^^

With the interactive method, prompts are displayed to the terminal, which lead you through configuring Lightning.
This method prevents any user error in the configuration, as responses to the prompts are required to be
allowed values. If a value is not allowed, a new prompt will be displayed stating your error and asking you
to update your value appropriately. Additionally, if you are unsure on what values to choose for each
input, we give default values for most configuration options. To use the non-default interactive method,
run Lightning in :ref:`interactive-label`.


.. highlight:: idl
.. _run-lightning-label:

Running Lightning
-----------------

Running Lightning only requires you to make a call to the main ``lightning`` procedure::

    IDL> lightning, '<input_dir>/<your_photometric_catalog>'

All Lightning processes take place within this procedure. They consist of:

1) Searching for a :ref:`configuration script <file-method-label>` (i.e., an edited copy of ``lightning_configure.pro``)
   within ``<input_dir>``, compiling the configuration script, checking it for errors, and saving the resulting
   configuration structure as ``<input_dir>/lightning_output/lightning_configure.sav``.

2) Loading the :ref:`input photometry <input-label>`, checking it for errors, and saving the photometry for
   each SED in the input catalog to individual IDL ``.sav`` files as
   ``<input_dir>/lightning_output/input_sav_files/lightning_input_<sed_id>.sav``.

3) Generating the models specified in the configuration and then fitting the models to the input data using
   the specified fitting algorithm. (If ``MAX_CPUS`` in the configuration is set to a number larger than ``1``,
   Lightning will parallelize the fitting where one SED is fit per CPU.)

4) Post-processing of the the resulting SED fits into more readily usable outputs.


Resuming
^^^^^^^^

If any issue caused Lightning to be interrupted mid-run, you can resume running where it left off
by passing the ``/resume`` flag to ``lightning``::

    IDL> lightning, '<input_dir>/<your_photometric_catalog>', /resume

.. warning::

    Resuming causes Lightning to search and use the saved configuration structure
    (``<input_dir>/lightning_output/lightning_configure.sav``). Without this ``.sav``
    file, Lightning cannot be resumed.

.. _interactive-label:

Interactive Mode
^^^^^^^^^^^^^^^^

You can run Lightning in interactive mode by passing the ``/interactive`` flag to ``lightning``::

    IDL> lightning, '<input_dir>/<your_photometric_catalog>', /interactive

In interactive mode, Lightning prints a series of prompts to the screen, allowing you to select a configuration without
any pre-existing configuration script.
