import time, sys, os
import gzip
import h5py
import numpy as np
import scipy
from matplotlib.pyplot import *
from astropy.table import Table

import fsps
import sedpy
import prospect

import emcee
import dynesty

from prospect.fitting import fit_model
from prospect.fitting import lnprobfn
import prospect.io.read_results as reader



def read_lightning_data():
    """
    This will need to be directly edited to set the FITS file path and
    update the filters labels from Lightning to sedpy.

    :returns lightning_data:
        A FITS table containing the input data in Lightning format.

    :returns sedpy_filters:
        A list of `Filter()` objects generated from `sedpy`.
    """
    from astropy.io import fits
    from astropy.io import ascii
    import sedpy

    # Read in Lightning data from FITS table
    hdul = fits.open('../ngc628_dale17_photometry.fits')
    lightning_data = hdul[1].data
    hdul.close()


    # Update the names of the filters from the FITS table, 
    # in the same order as the photometric data (see below)
    galex = ['galex_FUV', 'galex_NUV']
    sdss = ['sdss_{0}0'.format(b) for b in ['u','g','r','i','z']]
    bessell = ['bessell_'+n for n in ['B','V','R','I']]
    twomass = ['twomass_'+n for n in ['J','H','Ks']]
    irac = ['spitzer_irac_ch'+n for n in ['1','2','3','4']]
    mips = ['spitzer_mips_'+n for n in ['24','70','160']]
    wise = ['wise_w3']
    pacs = ['herschel_pacs_'+n for n in ['70','100','160']]
    spire = ['herschel_spire_'+n for n in ['250','350','500']]

    # Prospector does not have SCUBA or Planck HFI bands, so add them
    # We have these saved in the BAGPIPES example folder. So let's use those.
    hfi_data = ascii.read('../bagpipes/filters/HFI_F353.txt')
    scuba_data = ascii.read('../bagpipes/filters/SCUBA2_850.txt')

    # Wavelength units in bagpipes/filters is Angstrom, which is what sedpy wants
    hfi = sedpy.observate.Filter("hfi_353", data=(hfi_data["wave[Angstroms]"], hfi_data["norm_trans"]))
    scuba = sedpy.observate.Filter("scuba2_850", data=(scuba_data["wave[Angstroms]"], scuba_data["norm_trans"]))

    # Combine the filter names from above into a single list
    filternames = galex + [sdss[0]] + [bessell[0]] + [sdss[1]] + [bessell[1]] + [sdss[2]] + [bessell[2]] + \
                  [sdss[3]] + [bessell[3]] + [sdss[4]] + twomass + irac + wise + mips[0:2] + pacs[0:2] + \
                  [mips[2]] + [pacs[2]] + spire

    # Read the sedpy filters that are defaults in sedpy and add on the hfi and scuba filters we created
    sedpy_filters = sedpy.observate.load_filters(filternames) + [hfi] + [scuba]

    return lightning_data, sedpy_filters


# --------------
# OBSERVATIONS
# --------------


def build_obs(lightning_data, sedpy_filters):
    """
    :returns obs:
        A dictionary of observational data to use in the fit.

    :param lightning_data:
        The Lightning input FITS table that the required input observations.
        
    :param sedpy_filters: 
        A list `Filter()` objects generated from `sedpy` (Prospector may not contain
        the same filters as Lightning. So, we create this outside of this function).
    """
    from prospect.utils.obsutils import fix_obs

    # The obs dictionary, empty for now
    obs = {}

    # And here we instantiate the `Filter()` objects using methods in `sedpy`,
    # and put the resulting list of Filter objects in the "filters" key of the `obs` dictionary
    obs["filters"] = sedpy_filters

    # Now we store the measured fluxes for a single object, **in the same order as "filters"**
    # The units of the fluxes need to be maggies (Jy/3631) so we will do the conversion here too.
    obs["maggies"] = np.squeeze(lightning_data['fnu_obs']/3631)

    # And now we store the uncertainties (again in units of maggies)
    obs["maggies_unc"] = np.squeeze(lightning_data['fnu_unc']/3631)

    # Now we need a mask, which says which flux values to consider in the likelihood.
    # IMPORTANT: the mask is *True* for values that you *want* to fit, 
    # and *False* for values you want to ignore.
    obs["phot_mask"] = np.array(['scuba' not in f.name for f in obs["filters"]])

    # This is an array of effective wavelengths for each of the filters.  
    # It is not necessary, but it can be useful for plotting so we store it here as a convenience
    obs['phot_wave'] = np.array([filt.wave_effective for filt in obs['filters']])

    # We do not have a spectrum, so we set some required elements of the obs dictionary to None.
    # (this would be a vector of vacuum wavelengths in angstroms)
    obs["wavelength"] = None
    # (this would be the spectrum in units of maggies)
    obs["spectrum"] = None
    # (spectral uncertainties are given here)
    obs['unc'] = None
    # (again, to ignore a particular wavelength set the value of the 
    #  corresponding elemnt of the mask to *False*)
    obs['mask'] = None

    # This function ensures all required keys are present in the obs dictionary,
    # adding default values if necessary
    obs = fix_obs(obs)

    return obs


# --------------
# MODEL
# --------------


def build_model(lightning_data, **extras):
    """Build a prospect.models.SedModel object
    
    :returns obs:
        An instance of prospect.models.SedModel

    :param lightning_data:
        The Lightning input FITS table that the required input observations.
    """
    from prospect.models.sedmodel import SedModel
    from prospect.models.templates import TemplateLibrary
    from prospect.models import priors

    # Get the prepackaged non-parametric SFH model set dictionaries.
    # This is, somewhat confusingly, a dictionary of dictionaries, keyed by parameter name
    model_params = TemplateLibrary["logm_sfh"]
    # Add a dust emission model
    model_params.update(TemplateLibrary["dust_emission"])
    
    # Now add the lumdist parameter by hand as another entry in the dictionary.
    # This will control the distance since we are setting the redshift to zero.  
    model_params["lumdist"] = {"N": 1, "isfree": False, "init": lightning_data['lumin_dist'], "units":"Mpc"}

    # Let's fix the parameters to match the default Lightning setup.
    model_params["logzsol"]["isfree"] = False
    model_params["logzsol"]['init'] = 0.0
    model_params["zred"]["isfree"] = False
    model_params["zred"]['init'] = 0.0
    
    # Let's make some changes to the parameters settings
    # dust2 is the tauV parameter in Lightning, so let's give it the same prior
    model_params["dust2"]["init"] = 1.0
    model_params["dust2"]["prior"] = priors.TopHat(mini=0.0, maxi=10.0)

    # Set the umin dust emission prior
    model_params["duste_umin"]["init"] = 5.0
    model_params["duste_umin"]["isfree"] = True
    model_params["duste_umin"]["prior"] = priors.TopHat(mini=0.1, maxi=25.0)

    # Set the qpah dust emission prior
    # Note that qpah is in percentage, not fractional form like Lightning
    model_params["duste_qpah"]["init"] = 2.0
    model_params["duste_qpah"]["isfree"] = True
    model_params["duste_qpah"]["prior"] = priors.TopHat(mini=0.47, maxi=4.58)

    # Set the gamma dust emission prior
    model_params["duste_gamma"]["init"] = 0.1
    model_params["duste_gamma"]["isfree"] = True
    model_params["duste_gamma"]["prior"] = priors.TopHat(mini=0.0, maxi=1.0)

    # Update model parameters (i.e., dust attenuation and age bins)
    # Dust_type = 2 means base Calzetti attenuation
    model_params["dust_type"]["init"] = 2
    # Needs to be set to work correctly as stated in fsps manual
    model_params["dust1"] = {"N": 1, "isfree": False, "init": 0.0}

    # Now we adjust the age bins to match the age bins in Lightning
    model_params["agebins"]["N"] = 5
    # Note that the agelims are in log10 space. So we set the youngest age to 1 yr
    agelims = [0, 7, 8, 9, 9.69897, 10.1271]
    agebins = np.array([agelims[:-1], agelims[1:]])
    model_params["agebins"]["init"] = agebins.T
    model_params["mass"]["N"] = 5
    # Now we set the mass prior to match bagpipe's priors as much as possible
    model_params["mass"]["init"] = np.array([1.e+7, 1.e+7, 1.e+7, 1.e+7, 1.e+7])
    model_params["mass"]["prior"] = priors.TopHat(mini=np.array([0.0, 0.0, 0.0, 0.0, 0.0]), \
                                                  maxi=np.array([1.e+12, 1.e+12, 1.e+12, 1.e+12, 1.e+12]))

    # We are going to be using emcee. So, it is useful to provide a 
    # minimum scale for the cloud of walkers (the default is 0.1)
    # We set this to a reasonable range that is similar to Lightning's initialization range
    model_params["mass"]["disp_floor"] = 1.e+7
    model_params["dust2"]["disp_floor"] = 0.5
    model_params["duste_umin"]["disp_floor"] = 1.0
    model_params["duste_qpah"]["disp_floor"] = 1.0
    model_params["duste_gamma"]["disp_floor"] = 0.5
    

    # Now instantiate the model object using this dictionary of parameter specifications
    model = SedModel(model_params)

    return model


# --------------
# SPS OBJECT
# --------------

def build_sps():

    # Build the SPS models
    from prospect.sources import FastStepBasis
    sps = FastStepBasis()
    return sps



if __name__ == "__main__":

    from prospect.io import write_results as writer

    # See https://github.com/bd-j/prospector/blob/main/demo/InteractiveDemo.ipynb
    #   for helpful info.

    # Readin the data observations, create the models, and SPS
    lightning_data, sedpy_filters = read_lightning_data()
    obs = build_obs(lightning_data, sedpy_filters)
    model = build_model(lightning_data)
    sps = build_sps()
    
    # --- start minimization ----
    run_params = {}
    run_params["optimize"] = False
    run_params["emcee"] = True
    run_params["dynesty"] = False
    # Number of emcee walkers, matches Lightning
    run_params["nwalkers"] = 75
    # Number of iterations of the MCMC sampling, matches Lightning
    run_params["niter"] = 10000
    # Number of iterations in each round of burn-in
    # After each round, the walkers are reinitialized based on the 
    # locations of the highest probablity half of the walkers.
    # We would preferably not want to use this. However Prospector breaks if it does not have
    #   at least one burn in trial. So, that is what we do, and we will manually remove
    #   the burn-in phase on our own later.
    run_params["nburn"] = [1]

    hfile = "NGC0628_emcee_mcmc.h5"

    # Fit the model to the data and save the results
    if not os.path.isfile(hfile): 
        start = time.time()
        output = fit_model(obs, model, sps, lnprobfn=lnprobfn, **run_params)
        end = time.time()
        print('Prospector finished fitting in '+str(end-start)+'s.')

        print("Done optmization in {}s".format(output["sampling"][1]))

        writer.write_hdf5(hfile, run_params, model, obs,
                          output["sampling"][0], sps=sps,
                          tsample=output["sampling"][1])


    # grab results (dictionary), the obs dictionary, and our corresponding models
    result, obs, model = reader.results_from(hfile)

    # Compute the autocorrelation time to see if we can trust our parameter estimates
    autocorr_time = emcee.autocorr.integrated_time(np.transpose(result['chain'], (1, 0, 2)), quiet=True)
    # Remove 2x tau samples (i.e., burn-in) and recompute to get better estimate
    burn_in = 2 * np.ceil(np.max(autocorr_time))
    autocorr_time = emcee.autocorr.integrated_time(np.transpose((result['chain'])[:, burn_in.astype(int):, :], (1, 0, 2)), quiet=True)

    # Recreate the observation and models as the saved data we loaded is not the same
    obs = build_obs(lightning_data, sedpy_filters)
    model = build_model(lightning_data)

    # Truncate, thin, and merge MCMC results into final chain
    thin = 0.5 * np.ceil(np.max(autocorr_time))
    final_chain = result['chain'][:, burn_in.astype(int):run_params["niter"]:thin.astype(int), :]
    final_chain = np.reshape(final_chain, (run_params["nwalkers"]*final_chain.shape[1], 9), order='F')
    final_chain = final_chain[-2000:, :]
    lnprob = result['lnprobability'][:, burn_in.astype(int):run_params["niter"]:thin.astype(int)]
    lnprob = (lnprob.flatten(order='F'))[-2000:]
    lnprob = lnprob[-2000:]


    # Create SFH parameters for each mass bin
    agebins = model.params['agebins']
    dt = (10**agebins[:, 1] - 10**agebins[:, 0])
    masses = final_chain[:, 0:5]
    sfh = masses / np.broadcast_to(dt, masses.shape)

    # Create model photometry for each final chain element, and the spectra for the best-fit model
    wspec = sps.wavelengths
    mod_spec = np.zeros((2000, len(wspec)))
    mod_phot = np.zeros((2000, len(obs['phot_wave'])))
    survive_mass_frac = np.zeros(2000)
    for i in range(2000):
        theta = final_chain[i, :]
        mod_spec[i, :], mod_phot[i, :], survive_mass_frac[i] = model.predict(theta, obs, sps=sps)

    # Only keep the best-fit model spectra
    # Before we do that we need to calculate LTIR from it
    ltir_idc = ((wspec/1e4) >= 8) & ((wspec/1e4) <= 1000)

    # Converts Fnu in Maggies to Lnu in Lsun Hz-1
    conv = 3631 * 4 * np.pi * (7.2 * 1e6 * 3.08567758149e18)**2 * 2.6023368985348840e-57
    ltir_wave_grid = np.repeat((1e8 * 2.99792458e10 / wspec[ltir_idc])[np.newaxis, :], 2000, axis=0)
    ltir = -1 * np.trapz(mod_spec[:, ltir_idc] * conv, x=ltir_wave_grid, axis=1)

    mod_spec = mod_spec[np.argmax(lnprob), :]

    # Create FITS table to hold all necessary outputs
    t = Table()
    t['phot_mod'] = mod_phot[None, ...]
    t['wave_hires'] = wspec[None, ...]
    t['spec_mod'] = mod_spec[None, ...]
    t['lnprob'] = lnprob[None, ...]
    t['survive_mass_frac'] = survive_mass_frac[None, ...]
    t['theta_labels'] = (np.array(result['theta_labels']))[None, ...]
    t['sfh'] = sfh[None, ...]
    t['ltir'] = ltir[None, ...]
    t['chain'] = final_chain[None, ...]
    t['agebins'] = agebins[None, ...]
    t['burn_in'] = (burn_in.astype(int))[None, ...]

    t['ntrials'] = (np.array(result['run_params']['niter']))[None, ...]
    t['acceptance'] = result['acceptance'][None, ...]
    t['sampling_duration'] = (np.array(result['sampling_duration']))[None, ...]
    t['autocorr_time'] = autocorr_time[None, ...]

    # Save the table and gunzip it to reduce size
    t.write('ngc_628_prospector.fits', format='fits', overwrite=True)
    with open('ngc_628_prospector.fits', 'rb') as f_in, gzip.open('ngc_628_prospector.fits.gz', 'wb') as f_out:
        f_out.writelines(f_in)
    os.remove('ngc_628_prospector.fits')