.. _algorithm-select-label:

Selecting an Algorithm
======================

This guide will help you with selecting an algorithm and its corresponding hyper-parameters
needed to fit your SED(s).
Lightning currently has three fitting algorithms
to select from: the adaptive MCMC algorithm (Algorithm 4) from `Andrieu & Thoms (2008)
<https://link.springer.com/article/10.1007/s11222-008-9110-y>`_, the affine-invariant MCMC
algorithm from `Goodman & Weare (2010)
<https://ui.adsabs.harvard.edu/abs/2010CAMCS...5...65G/abstract>`_, and a Levenberg–Marquardt
algorithm implemented via `Craig Markwardt’s MPFIT <http://purl.com/net/mpfit>`_.
Below, we first compare the MPFIT and MCMC algorithms and give advice on how to select between
them. Then, we explain how the algorithms are initialized. Finally, we describe each algorithm
and its corresponding hyper-parameters in more detail.


MCMC vs MPFIT
-------------

Like all other settings in Lightning, the choice of algorithm will depend on your research goals.
The main differences between the MCMC algorithms and the MPFIT algorithm in practical terms are
the speed of the algorithms and the form of the uncertainties on the free parameters.
Since the MPFIT algorithm just searches for the best-fit solution rather than probing the shape
of the posterior distribution like the MCMC algorithms, it is requires significantly fewer loss
function evaluations and hence is significantly faster than the MCMC algorithms. However,
this speed comes at a loss in terms of the uncertainty estimations. The MPFIT algorithm simply
outputs a :math:`1\sigma` uncertainty estimation on each free parameter. Comparatively,
the MCMC algorithms give a full Bayesian posterior distribution for each free parameter.
Therefore, if you are looking for a quick best-fit solution to your SEDs with uncertainties
you will not heavily rely on for other analyses, we recommend using the MPFIT algorithm.
However, if you are looking for a full Bayesian posterior distribution with well estimated uncertainties,
we recommend using the MCMC algorithms.


.. _random-initialize-label:

Random Initialization
---------------------

All three algorithms in Lightning require an initial starting value for each free parameter.
Currently, Lightning automatically selects the starting values by randomly sampling the prior
distribution of each parameter independently. Since uniform priors can have a much larger
range than we would want to randomly sample, we allow the user to input an initialization range.
This limits the random sampling of the prior to only occur within the initialization range.

To give an example of this, say we have a uniform prior for a parameter with a range of 0 to 100.
However, we typically expect this parameter's solution to be less than 10. If
we were to randomly initialize using the full prior and randomly select a value near the maximum
limit, our algorithms may not reach the best solution within its allotted iteration limit.
Therefore, by setting the initialization range to be between 0 and 10, we allow for a more appropriate
initial value and increase our chances to reach the best solution.


Adaptive MCMC
-------------

Description
^^^^^^^^^^^

The adaptive MCMC algorithm is an adaptive version of the standard `Metropolis–Hastings algorithm
<https://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm>`_ created by `Andrieu & Thoms (2008)
<https://link.springer.com/article/10.1007/s11222-008-9110-y>`_. The algorithm simply adjusts the
proposal density distribution to achieve an optimal acceptance ratio (~25%). Additionally, this
adjustment of the proposal density is vanishing, meaning the adaptiveness decreases with
each subsequent iteration. Therefore, after many iterations the adaptiveness is insignificant
and the algorithm is practically equivalent to the standard Metropolis–Hastings algorithm.


Hyper-parameters
^^^^^^^^^^^^^^^^

The adaptive MCMC algorithm has three hyper-parameters ``NTRIALS``, ``NPARALLEL``, and ``BETA_EXPONENT``.
The value of ``NTRIALS`` gives the number of MCMC trials (iterations) to run for each parallel chain,
the number of which is given by ``NPARALLEL``. Having ``NPARALLEL = 1`` will only run a single MCMC chain,
while having ``NPARALLEL > 1`` will run that many **independent** parallel MCMC chains. Running a single MCMC chain
presents the issue of us not being able to determine if we have reached the best solution or are stuck in a local
minimum in :math:`\chi^2` space, since we have no way of comparing the chain's solution to the best solution.
By running multiple parallel chains from different starting locations, we can compare the ending
segment of each parallel chain. If all or most of these segments have reached the same solution,
then we can be confident that this is likely the best solution (see :ref:`convergence-describe-label`
for more details on what qualifies as the same solution). Therefore, we recommend setting ``NPARALLEL``
to be a value between ``5`` and ``20`` to determine if the chains have reached the best solution.
A larger value will allow you to be more confident that you have reached the best solution, but will
require more computational time. Additionally, we recommend setting ``NTRIALS`` to an integer multiple
of ``5e4`` for each :ref:`model component <model-select-label>` you are using. If your parallel chains
are not converging (see the :ref:`convergence-describe-label` description for more details),
we recommend increasing this value by a few extra integer multiples to give more
iterations for convergence to occur.

``BETA_EXPONENT`` is the parameter unique to the adaptive MCMC algorithm. It determines how fast
the adaptiveness decreases with each subsequent iteration. Values must be positive, and larger values
stop the adaptiveness in fewer trials. We recommend using the default value for ``BETA_EXPONENT``
in most cases. However, if your parallel chains are not converging even after increasing the number
of trials, we recommend decreasing this value slightly (e.g., from ``0.35`` to ``0.3``).


Post-Processed Chain
^^^^^^^^^^^^^^^^^^^^

During post-processing, only a portion of the raw MCMC chain from fitting is kept for each parameter.
For the adaptive MCMC, this post-processed chain portion is determined as follows:

1) the parallel chains have ``BURN_IN`` number of initial trials discarded as the burn-in phase,
2) these truncated chains are then thinned by keeping only every ``THIN_FACTOR`` elements,
3) the thinned and truncated chain with the highest probability is selected,
4) the final ``FINAL_CHAIN_LENGTH`` elements from this highest probability chain are kept as the
   post-processed chain.

Here, we explain our reason for selecting the highest probability chain in step 3. Assuming
convergence of the parallel chains has been reached, all parallel chains will have very similar
distributions, and it does not matter which one we select to use. Therefore, we find selecting
the one with the highest probability guarantees it is at the best solution even when convergence
is not reached.


.. _affine-mcmc-label:

Affine-Invariant MCMC
---------------------

Description
^^^^^^^^^^^

The affine-invariant MCMC uses an ensemble of samplers to adjust the proposal density distribution
and sample the posterior distribution. This ensemble consists of multiple chains (or walkers) that
are run in parallel and allowed to interact with one another so that they can adapt their proposal
densities. Partial re-sampling, which is a generalized version of the `Gibbs sampling
<https://en.wikipedia.org/wiki/Gibbs_sampling>`_, is implemented to achieve this interaction. In
our implementation, we use the affine-invariant stretch move method as presented in `Goodman & Weare
(2010) <https://ui.adsabs.harvard.edu/abs/2010CAMCS...5...65G/abstract>`_.


Hyper-parameters
^^^^^^^^^^^^^^^^

The affine-invariant MCMC algorithm has three hyper-parameters ``NTRIALS``, ``NPARALLEL``, and ``AFFINE_A``.
The value of ``NPARALLEL`` gives the number of walkers to include in the ensemble, and ``NTRIALS`` gives
the number of MCMC trials (iterations) to run for each walker. Unlike the adaptive MCMC, the affine-invariant
MCMC must have ``NPARALLEL > 1``. Specifically, ``NPARALLEL`` must be greater than the number of free
parameters plus one, and ideally, it should be at least twice the number of free parameters for optimal
sampling. Therefore, we recommend setting ``NPARALLEL`` to be a value between ``50`` and ``100``, which is
3 to 5 times the maximum number of free parameters that is expected from Lightning's most complex models.
Additionally, we recommend setting ``NTRIALS`` to an integer multiple of ``1e4`` for each
:ref:`model component <model-select-label>` you are using. If your ensemble is not converging
(see the :ref:`convergence-describe-label` description for more details), we recommend increasing
this value by a few extra integer multiples to give more iterations for convergence to occur.

``AFFINE_A`` is the parameter unique to the affine-invariant MCMC algorithm. It specifies
the move scaling constant, which defines the maximum and minimum step size of the stretch move.
Values must be greater than or equal to 1, and larger values allow for larger stretch moves in
parameter space. We recommend using the default value for ``AFFINE_A``
in most cases. However, if your ensemble is not converging even after increasing the number
of trials or has a low overall acceptance rate (< 20%), we recommend decreasing this value slightly
(e.g., from ``2`` to ``1.8``).


Post-Processed Chain
^^^^^^^^^^^^^^^^^^^^

During post-processing, only a portion of the raw MCMC ensemble from fitting is kept for each parameter.
For the affine-invariant MCMC, the post-processed chain portion is determined as follows:

1) each walker in the ensemble has ``BURN_IN`` number of initial trials discarded as the burn-in phase,
2) if a walker has an acceptance fraction less than ``AFFINE_STRANDED_DEVIATION`` standard deviations
   below the median acceptance fraction, we consider them stranded walkers and remove them from the ensemble,
3) the non-stranded truncated ensemble is then thinned by keeping only every ``THIN_FACTOR`` elements,
4) the thinned and truncated ensemble is flattened element-wise into a single chain,
5) the final ``FINAL_CHAIN_LENGTH`` elements from this flattened chain are kept as the post-processed chain.

Here, we explain our reason for removing stranded walkers in step 2. Due to the boundaries of the
free parameters, the affine-invariant MCMC can have trouble accepting moves of walkers separated from
the ensemble when the ensemble is near a boundary. This results in the walkers becoming stranded and
having a very low acceptance rates, since they are failing to have any proposal jumps accepted. With
enough iterations, these walkers will get lucky and have a jump that rejoins them with the ensemble.
However, we do not have an infinite amount of iterations to allow for this to occur. Therefore, once
our iteration limit has been reached, we want to remove any stranded walkers that may remain. We have
found that the most effective method for correctly selecting stranded walkers is to compare each walker's
acceptance fraction with that of the median of the ensemble. Those that have an abnormally low
acceptance fractions compared to the rest of the ensemble are usually stranded.

.. note::

    We find that only a few walkers within the ensemble become stranded when using a standard amount
    of iterations. Therefore, having ``AFFINE_STRANDED_DEVIATION = 2`` effectively removes these walkers without
    removing non-stranded ones. However, when using a smaller amount of iteration for quick sampling, more
    walkers may end up remaining stranded. Therefore, we recommend setting ``AFFINE_STRANDED_DEVIATION = 1`` to
    account for the increase in the ensemble's standard deviation and better classify stranded walkers.



Adaptive vs Affine-Invariant MCMC
---------------------------------

The main differences between the affine-invariant MCMC and the adaptive MCMC algorithms is their
speed and consistency for reaching the best solution. From some general tests, we find that
both algorithms result in very similar posterior distributions if the best solution is reached,
as should be expected. However, the adaptive MCMC algorithm is less effective at searching parameter
space for the best solution. It can spend a significant portion, if not all, of its trials
stuck in local minima if it does not start near the best solution, especially with more complex models.
In comparison, we find the affine-invariant MCMC algorithm regularly reaches the best solution rapidly
without getting stuck in local minima. Therefore, we recommend using the affine-invariant MCMC over
the adaptive MCMC algorithm as it more consistently reaches the best solution and requires
less cost function evaluations to get the needed posterior distribution.


.. _mpfit-algorithm-label:

MPFIT
-----

Description
^^^^^^^^^^^

The MPFIT algorithm is `Craig Markwardt’s implementation <http://purl.com/net/mpfit>`_ 
of the gradient-descent `Levenberg–Marquardt algorithm <https://en.wikipedia.org/wiki/Levenberg–Marquardt_algorithm>`_,
which is used to solve non-linear least squares problems. The MPFIT implementation allows for several
necessary constraints in Lightning, such as fixing parameters and setting parameter bounds. Additionally,
the algorithm calculates the parameter covariance matrix to give estimated parameter uncertainties.


Hyper-parameters
^^^^^^^^^^^^^^^^

The MPFIT algorithm has five hyper-parameters ``NSOLVERS``, ``FTOL``, ``GTOL``, ``XTOL``,
and ``MAXITER``. The value of ``NSOLVERS`` gives the number of "solvers" to run in parallel,
where each solver is a fit to the SED using different starting locations in parameters space.
Numerous solvers are necessary, since like the adaptive MCMC algorithm, running a single solver
presents the issue of us not being able to determine if we have reached the best solution or are
stuck in a local minimum. By running multiple solvers from different starting locations, we can
compare each solver's solution. If the majority of the solvers have reached the same solution,
then we can be confident that this is likely the best solution. Therefore, we recommend setting
``NSOLVERS`` to an integer multiple of ``50`` for each :ref:`model component <model-select-label>`
you are using to determine if the solvers have reached the best solution. A larger value will allow
you to be more confident that you have reached the best solution, but will require more computational
time.

``FTOL``, ``GTOL``, and ``XTOL`` give the tolerances indicating when the MPFIT algorithm should
terminate. Smaller values of the tolerances mean MPFIT will continue to run until smaller
differences are produced in the relative error. We recommend using the default values for each
of the tolerances. However, if you find that the majority of your solvers are not finding the
same best solution, then we recommend decreasing ``FTOL`` or ``XTOL`` slightly (e.g., from
``1d-10`` to ``1d-12``) if the MPFIT ``STATUS`` is ``1`` or ``2``, respectively. This may allow
for some solvers to escape local minima and reach the best solution.

Finally, ``MAXITER`` gives the maximum number of MPFIT iterations to perform per solver. If the
MPFIT algorithm has not terminated already from reaching one of the tolerances, then it will
terminate after performing this maximum number of iterations. This number is to prevent MPFIT
from potentially running indefinitely if the tolerances cannot be reached. We recommend using
the default value for ``MAXITER``. However, if a reasonable portion of your set of solvers is
reaching the maximum iterations, then you can increase this value to allow for more iteration
for the tolerances to be reached.


Post-Processed Fits
^^^^^^^^^^^^^^^^^^^

The post-processing for the MPFIT algorithm is simple. The solver with the best-fit solution
(i.e., lowest :math:`\chi^2`) is kept, and its parameter and uncertainty estimations are used
as the best fit. The rest are discarded as they were only needed to test convergence to the
best solution.
