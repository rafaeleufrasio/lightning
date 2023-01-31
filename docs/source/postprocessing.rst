.. _postprocessing-label:

Results and Post-Processing
===========================

.. toctree::
    :hidden:

    outputs/convergence_metrics.rst
    outputs/outputs.rst


Lightning automatically post-processes the resulting intermediate SED fit data into more
readily usable outputs. The final post-processed output is a single data table saved in the first
extension of a FITS file as ``<input_dir>/lightning_output/<OUTPUT_FILENAME>.fits.gz``, where
``<OUTPUT_FILENAME>`` is the value set in the :ref:`configure-setting-label`.
This data table contains the resulting outputs for all SEDs included in the input. 
The data within the table contains a variety of outputs depending on your chosen configuration.
This always includes the modified input data, such as the luminosities to which the model was fit
(converted from the input fluxes) and the distance indicator. It will additionally include the 
basic results from the fitting, such as the estimated free parameters and model luminosities. 


Also included in the table are the :ref:`converge-good-label`. These outputs can be used
to evaluate if the fitting algorithm converged to a solution and if the resulting solution
reasonably fits the data. Checking these outputs is essential to ensure that the other
outputs can be believed. For the convergence metrics,
each fitting algorithm has several of their own metrics, which can be used to check for convergence.
Therefore, a detailed discussion is dedicated to these :ref:`convergence-describe-label`.
For the goodness of fit, this is simply an estimated `p`-value that
indicates how well the fit model compares to the data. If the fitting algorithm converged,
these values can indicate if an SED was poorly fit by the model. 

Full details about all the outputs, and which ones are included in the post-processed file, can be found in
the :ref:`outputs-label`.
