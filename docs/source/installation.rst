============
Installation
============

The latest Lightning release can be installed as source from our `GitHub <https://github.com/rafaeleufrasio/lightning>`_ page via ``git``::

    cd <install_dir>
    git clone --depth 1 https://github.com/rafaeleufrasio/lightning

.. note::

    Lightning is written in the Interactive Data Language (IDL) and requires IDL version 8.3 or later. It additionally
    requires the `IDL Astronomy User's Library <https://idlastro.gsfc.nasa.gov>`_, the `IDL Coyote Library <http://www.idlcoyote.com>`_, and
    `Craig Markwardt's MPFIT library <http://purl.com/net/mpfit>`_. All three will need to be installed and `added to your IDL
    path <https://www.l3harrisgeospatial.com/Support/Self-Help-Tools/Help-Articles/Help-Articles-Detail/ArtMID/10220/ArticleID/16156/Quick-tips-for-customizing-your-IDL-program-search-path>`_
    before using Lightning.


Before running Lightning, you will need to build the package using::

    cd <install_dir>/lightning
    idl lightning_build

.. highlight:: idl

This creates a ``.sav`` file in the Lightning base directory containing all the compiled routines you need to run Lightning.
Then, at the beginning of an IDL session where you want to run Lightning, you should do::

    IDL> restore, !lightning_dir+'lightning.sav'

Or, you can add the above line of code to your `IDL startup file <https://www.l3harrisgeospatial.com/Support/Self-Help-Tools/Help-Articles/Help-Articles-Detail/ArtMID/10220/ArticleID/18093/How-do-I-specify-a-program-to-automatically-run-when-my-IDL-session-starts-up>`_.

.. note::

    We recommend adding this line of code to your startup file. If you do not already have one,
    a startup file will be automatically generated when building Lightning. The file will be
    named ``startup.pro``, and it will be located in the current user's ``.idl`` directory in the
    home directory.
