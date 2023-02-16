import bagpipes as pipes
import numpy as np
from astropy.io import fits
from astropy.table import Table
import os
import gzip
import time


# load_data function needed to generate observational data
def load_data(ID):

    # Read in Lightning data from FITS table
    hdul = fits.open('../ngc628_dale17_photometry.fits')
    lightning_data = hdul[1].data
    hdul.close()
    
    # Extract the fluxes and uncertainties, converting from Jy to uJy
    fluxes = np.squeeze(lightning_data['fnu_obs']) * 1e6
    fluxerrs = np.squeeze(lightning_data['fnu_unc']) * 1e6

    # Turn these into a 2D array.
    photometry = np.c_[fluxes, fluxerrs]

    # blow up the errors associated with any missing fluxes.
    for i in range(len(photometry)):
        if (photometry[i, 0] == 0.) or (photometry[i, 1] <= 0):
            photometry[i,:] = [0., 9.9*10**99.]
            
    return photometry



# Read in Lightning data from FITS table to get filters
hdul = fits.open('../ngc628_dale17_photometry.fits')
lightning_data = hdul[1].data
hdul.close()

# Load filter list
# Can't directly use Lightning's filters since wavelength is in microns, while Bagpipes needs Angstroms
# Converted filters from Lightning to have wavelength units in Angstroms and saved to bagpipes dir
filt_list = np.squeeze(lightning_data['filter_labels']).tolist()
filt_list = ['filters/'+x.strip(' ')+'.txt' for x in filt_list]

# Generate observational data
galaxy = pipes.galaxy("NGC_0628", load_data, spectrum_exists=False, filt_list=filt_list)

# Set the model we want to use, as close to that of Lightning as we can get.
# Use Calzetti attenauation
# By default Bagpipes uses uniform priors, so no need to set specifically
# Just need to set the pior range of tauV = 0.4*ln(10)*Av    (tauV_range = (0, 10))
dust = {}                                # Dust component
dust["type"] = "Calzetti"                # Attenuation law: "Calzetti", "Cardelli", "CF00" or "Salim"
dust["Av"] = (0.0, 10/(0.4*np.log(10)))  # Absolute attenuation in the V band: magnitudes

# Dust emission parameters
# Bagpipes has Draine & Li (2007) model with alpha=2 and Umax=1e6 permanently fixed
# Umax = 1e6 is close enough to Lightning's 3e5
# Set the uniform prior range for qpah, umin, and gamma
# Note that qpah is in percentage, not fractional form like Lightning
dust["qpah"] = (0.47, 4.58)       # PAH mass fraction
dust["umin"] = (0.1,  25.0)       # Lower limit of starlight intensity distribution
dust["gamma"] = (0.0,  1.0)       # Fraction of stars at umin

# Create SFH
# We can get the same SFH as Lightning by using a collection of constant (top-hat) functions
# Lightning SFH age limits are STEPS_BOUNDS: [0.d0, 1.d7, 1.d8, 1.d9, 5.d9, 1.34e10] in yr
# Lightning has linear SFH vs log mass for normalization. However, this is not an option.
# So, we will have to set the log mass range to be something reasonable, and the same as Prospector.
# Also metallicity is fixed at Z=0.02
constant1 = {}                       # tophat function
constant1["age_max"] = 1e-2          # Time since SF switched on: Gyr
constant1["age_min"] = 0.0           # Time since SF switched off: Gyr
constant1["massformed"] = (0., 12.)  # Log_10 total stellar mass formed: M_Solar
# pow_10 prior gives a linear uniform prior since mass is by default in log10
constant1["massformed_prior"] = 'pow_10'  # Log_10 total stellar mass formed: M_Solar
constant1["metallicity"] = 0.02      # Metallicity: Z_sol = 0.02

constant2 = {}                       # tophat function
constant2["age_max"] = 1e-1          # Time since SF switched on: Gyr
constant2["age_min"] = 1e-2          # Time since SF switched off: Gyr
constant2["massformed"] = (0., 12.)  # Log_10 total stellar mass formed: M_Solar
constant2["massformed_prior"] = 'pow_10'  # Log_10 total stellar mass formed: M_Solar
constant2["metallicity"] = 0.02      # Metallicity: Z_sol = 0.02

constant3 = {}                       # tophat function
constant3["age_max"] = 1e0           # Time since SF switched on: Gyr
constant3["age_min"] = 1e-1          # Time since SF switched off: Gyr
constant3["massformed"] = (0., 12.)  # Log_10 total stellar mass formed: M_Solar
constant3["massformed_prior"] = 'pow_10'  # Log_10 total stellar mass formed: M_Solar
constant3["metallicity"] = 0.02      # Metallicity: Z_sol = 0.02

constant4 = {}                       # tophat function
constant4["age_max"] = 5e0           # Time since SF switched on: Gyr
constant4["age_min"] = 1e0           # Time since SF switched off: Gyr
constant4["massformed"] = (0., 12.)  # Log_10 total stellar mass formed: M_Solar
constant4["massformed_prior"] = 'pow_10'  # Log_10 total stellar mass formed: M_Solar
constant4["metallicity"] = 0.02      # Metallicity: Z_sol = 0.02

constant5 = {}                       # tophat function
constant5["age_max"] = 1.34e1        # Time since SF switched on: Gyr
constant5["age_min"] = 5e0           # Time since SF switched off: Gyr
constant5["massformed"] = (0., 12.)  # Log_10 total stellar mass formed: M_Solar
constant5["massformed_prior"] = 'pow_10'  # Log_10 total stellar mass formed: M_Solar
constant5["metallicity"] = 0.02      # Metallicity: Z_sol = 0.02


# Create fit instructions which contains the model we want to fit
fit_instructions = {}                   # The fit instructions dictionary

# Unsure how distance is determined if redshift is 0 and there is no way to specify distance...
# Instead we will set the redshift such that it gives the expected distance
# z = 0.001679 leads to distance of ~7.2 Mpc using Bagpipes cosmology functions
fit_instructions["redshift"] = 0.001679      # Observed redshift  
fit_instructions["dust"] = dust
fit_instructions["constant1"] = constant1
fit_instructions["constant2"] = constant2
fit_instructions["constant3"] = constant3
fit_instructions["constant4"] = constant4
fit_instructions["constant5"] = constant5

start = time.time()
# Prep the fit
fit = pipes.fit(galaxy, fit_instructions, n_posterior=2000)

# Run the fit and time
if not os.path.isfile('pipes/posterior/NGC_0628.h5'):
    fit.fit(verbose=True, n_live=1000)
    duration = time.time() - start
    print('Bagpipes finished fitting in '+str(duration)+'s.')

# Generate the advanced quantities that BAGPIPES can generate
fit.posterior.get_advanced_quantities()

# Generate a simple model to just get the high res model wavelengths
model_components = {}
model_components["redshift"] = 0.001679      # Observed redshift  

constant1 = {}                       # tophat function
constant1["age_max"] = 1e-2          # Time since SF switched on: Gyr
constant1["age_min"] = 0.0           # Time since SF switched off: Gyr
constant1["massformed"] = (fit.posterior.samples['constant1:massformed'])[0]   # Log_10 total stellar mass formed: M_Solar
constant1["metallicity"] = 0.02      # Metallicity: Z_sol = 0.02

model_components["constant1"] = constant1

mod = pipes.model_galaxy(model_components, filt_list=filt_list)
# Get the high res model wavelengths since they are not in default output
spec_wave = mod.wavelengths

# Create astropy table for saving output dictionary to FITS file
t = Table()
for key in fit.posterior.samples.keys():
    t[key] = fit.posterior.samples[key][None, ...]

# Remove the columns in the table that we do not need for our analysis
t.remove_columns(['sfr', 'ssfr', 'nsfr', 'mass_weighted_age', 'tform', 
                  'tquench', 'uvj', 'dust_curve', 'formed_mass'])

# The SFH output by BAGPIPES is very detailed in time and very redundant.
#   Let's remove the repeated data
_, bags_sfh_uniq_idc = np.unique(t['sfh'], axis=2, return_index=True)
bags_sfh_uniq_idc = np.sort(bags_sfh_uniq_idc)

# We don't want the last index as it is 0 from the last age bin being not quite the age of the universe in Bagpipes
bags_sfh_uniq_idc = bags_sfh_uniq_idc[:-1]
t['sfh'] = (t['sfh'])[:, :, bags_sfh_uniq_idc]

# We only need the best fit model spectrum, not one for each sample
#   So, let's select it and remove the rest
#   First though we need to calculate LTIR from this spectrum
ltir_idc = ((spec_wave/1e4) >= 8) & ((spec_wave/1e4) <= 1000)

# Converts Flambda in erg/s/cm2/Angstrom to Lnu in Lsun Hz-1
conv = 4 * np.pi * (7.2 * 1e6 * 3.08567758149e18)**2 * 2.6023368985348840e-57
spec_conv = np.repeat((1e23 * conv * (spec_wave)**2 /(2.99792458e18))[np.newaxis, :], 2000, axis=0)

ltir_wave_grid = np.repeat((1e8 * 2.99792458e10 /  spec_wave[ltir_idc])[np.newaxis, :], 2000, axis=0)
ltir     = -1 * np.trapz(((t['spectrum_full'])[:, :, ltir_idc])[0, :, :] * spec_conv[:, ltir_idc],
                         x=ltir_wave_grid, axis=1)
t['ltir'] = ltir[None, ...]


t['spectrum_full'] = (t['spectrum_full'])[:, np.argmin(t['chisq_phot'])]

# Add the high res wavelength grid and time to run the fit
t['spec_wave'] = spec_wave[None, ...]
t['time'] = (np.array(duration))[None, ...]

# Save the table to a FITS file
t.write('ngc_628_bagpipes.fits', format='fits', overwrite=True)

# Gunzip the FITS file to reduce file size
with open('ngc_628_bagpipes.fits', 'rb') as f_in, gzip.open('ngc_628_bagpipes.fits.gz', 'wb') as f_out:
    f_out.writelines(f_in)
os.remove('ngc_628_bagpipes.fits')
