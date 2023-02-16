.. _model-select-label:

Selecting a Model
=================

This guide will help you with selecting a model to fit to your SED(s). Typically,
the choice of model strongly depends on your data and research goals. For example,
you would obviously not use an X-ray model if you did not have any X-ray data, or if
you are studying a sample of AGN, you would likely want to include an AGN model component.
To help you select the model that best fits your research goals, below we discuss why
you would or would not want to include a certain component.

.. note::

    Most model components in Lightning are independent, meaning that a component
    can be excluded from the model without impacting the other components. The exceptions
    to this are:

    1) the X-ray component, which requires a stellar component;
    2) the X-ray AGN component, which requires an UV-to-IR AGN component;
    3) when energy balance is set and a dust emission component is used, a stellar or AGN component is required.


.. _stellar-emission-model:

Stellar Emission
----------------

Typically, the purpose for fitting a galaxy's SED is to derive its stellar mass and SFR by modeling its SFH.
Therefore, it is most likely that you will want to include a stellar component. One of the few reasons you may not
want a stellar component is if you are fitting a quasar or like object whose emission
is dominated by the AGN. Another reason may be if you only have IR data and are wanting
to just fit a dust emission model to your IR SED. However,
in general, we recommend using a stellar emission component unless you know that one
is not needed.

As for the parameters in the stellar model, the only two we discuss further here are
the metallicity and SFH age bins (``STEPS_BOUNDS``). As noted in the
:ref:`configure-setting-label`, Lightning currently assumes the chosen metallicity
is constant for all ages. This assumption ignores any underlying metallicity evolution
(however strong or weak that may be), and it could cause systematic variation in the resulting
SFHs estimates. For example, as metallicity decreases for our models, the intrinsic UV-optical
emission increases for a fixed SFR. This can lead to slightly decreased SFR estimates
for the younger populations of the SFH, assuming fixed attenuation, due to the younger
populations dominating the UV-optical emission. Since this effect is typically limited to
the younger population, we recommend selecting a metallicity closest to current average
metallicity of your input SEDs to minimize this systematic effect. This can be from
measured metallicities or from a redshift-metallicity relation (e.g., `Madau & Fragos 2017
<https://ui.adsabs.harvard.edu/abs/2017ApJ...840...39M/abstract>`_).

Selecting the appropriate age bins for the stellar model is not simple task, as there are several
ways to space the bins and select the number of bins. We find that equally spacing the bins in log
space (besides the first bin, which must be manually set as it should begin at 0 yr)
works well for most general uses. As for the number of bins, this will depend on the number of
photometric data points that you have per SED. If you have 20+ data points, you could easily
use 7 or more bins if you are using a relatively simple model. However, as
you increase the number of bins, the correlation between bins increases and can increase each
bin's uncertainty. Therefore, **as a rule of thumb we recommend using 5-7 log spaced stellar age
bins for your stellar models.**


.. _dust-attenuation-model:

Dust Attenuation
----------------

The dust attenuation curve determines how the stellar and/or AGN emission is attenuated.
The most commonly used attenuation curve to fit large samples of SEDs is the Calzetti+00
attenuation curve due to its simple parametric form. However, the curve is inflexible and can
have problems fitting a variety of SEDs. This is why the modified Calzetti attenuation
curve was introduced and has since taken over as the primarily used attenuation curve in
SED fitting codes. Unless your research needs involve you wanting to use the Doore+21 attenuation
curve, we highly recommend using the modified Calzetti curve with a UV bump, variable slope,
and no birth cloud attenuation.


.. _dust-emission-model:

Dust Emission
-------------

Since the majority of galaxies in the universe contain dust, using a dust emission model is typically
recommended unless you know that one is not needed. The only time we recommend not using a dust
emission model is if you do not have the data
to constrain it. For example, in our models, dust emission does not begin to dominate the stellar
emission until rest-frame wavelengths > 5 :math:`\mu {\rm m}`. Therefore, if you do not have any
data beyond 5 :math:`\mu {\rm m}` rest-frame, we recommend not using a dust emission model.
For best results constraining the dust model, we recommend having several IR observations ranging
from 5 to 500 :math:`\mu {\rm m}` rest-frame, with one around 100 :math:`\mu {\rm m}` to constrain the
peak of the dust emission. However, this many observations are not required, but having fewer will
result in larger uncertainties on the dust emission.

Finally, we recommend having energy balance set if you are using the dust emission model and a stellar and/or
AGN component. While, energy balance will not influence the dust emission model, it can better constrain the
attenuated portion of the stellar and/or AGN models. From our tests, using energy balance typically results in
approximately the same solution, but it reduces
the uncertainty on the attenuation parameters and younger SFH coefficients by a factor or 3 compared to not using
energy balance.


.. _xray-emission-model:

X-ray Emission
--------------

The X-ray emission model has two components, a stellar X-ray binary component and an optional AGN X-ray
emission component. As noted above, including an X-ray emission model requires there be a stellar
emission model to tie to the stellar X-ray binary component. Additionally, the optional AGN X-ray
emission component requires an AGN emission model. In general, we recommend including
an X-ray emission component, if you have X-ray data and a stellar emission model, as it can help
constrain the stellar emission. The same goes for the AGN X-ray emission component. If you have a UV-to-IR
AGN emission model, we recommend including the AGN X-ray emission component.

As for the parameters in the X-ray emission model, the only two we discuss further are the
X-ray uncertainty type and absorption model. For the uncertainty type, the choice depends on
your X-ray data. For data with large numbers of counts, the ``'SQRT'``
uncertainties may be appropriate. For smaller numbers of counts, we recommend using the ``'GEHRELS'`` uncertainties.
If you have your own uncertainties on your counts (if, e.g. the background contribution to the uncertainties is significant,
or you wish to adjust the weighting of the X-ray data relative to the UV-IR data), you can include them in the input X-ray
data and select ``'USER'`` uncertainties. As for the X-ray absorption model, we recommend using the
``'TBABS-WILM'`` model for most cases.


.. _agn-emission-model:

AGN Emission (UV-to-IR)
-----------------------

The AGN emission model only includes AGN emission in the UV-to-IR. This includes both the broken power law UV-to-NIR
emission from the accretion disk and the reprocessed IR dust emission from the surrounding torus. In general,
you would only want to include an AGN emission model if you know that your SED contains AGN emission.
Including one otherwise could lead to degeneracies between free parameters in the the other emission
model components. Therefore, we recommend only using the AGN emission model if you know your SED
contains AGN emission or if including one is relevant to your research goals. In cases where
quality NIR-MIR data are available, fitting an AGN model can also be used to select AGN by comparing fit quality with
and without an AGN component.
