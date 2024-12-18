from BP_PreLoad import *
import bagpipes as pipes
import pandas as pd
from astropy.io import fits
import numpy as np
import time
from astropy.cosmology import LambdaCDM

if '--survey' in args:
        
    indx = args.index("--survey")

    survey = sys.argv[indx + 1]
    
else:
    print('No survey detected, Defaulting to CANDELS')
    survey = 'CEERS'

print('Grabbing the necessary filter curves')
filters = make_filters(filter_set = survey)

def get_redshift(DF, ID):
    
    z_tesla = DF.loc[ID, 'Redshift']

    return z_tesla

def Galaxy_Model_Builder(ID, load_func, filters):
    '''
    This function will make a bagpipes galaxy model for one ID
    
    '''
    #need to take the extra step of converting ID to string for bagpipes galaxy function to work properly
    galaxy_bp = pipes.galaxy(str(ID), load_func, filt_list = filters, spectrum_exists=False)
    
    return galaxy_bp


def load_phot_CEERS(ID):
    
    ##########
    #WILL NEED TO CHANGE THIS BASED OFF OF THE INPUT CATALOG
    #WE NEED TO DO THIS 
    ##########
    Bagpipes_Phot_DF = pd.read_csv('../data/MARISSA_CEERS_PHOTOM_MATCHED.txt', 
                                   sep = ' ', index_col = 0)
    
    ID = int(ID)
    
    #defining the columns we will use in the photometry
    flux_cols = flux = [x for x in Bagpipes_Phot_DF.columns.values if 'FLUX_F' in x]

    flux_err_cols = [x for x in Bagpipes_Phot_DF.columns.values if 'FLUXERR_F' in x]

    
    
    #getting the full flux and flux error info
    photom_flux = Bagpipes_Phot_DF.loc[ID, flux_cols]
    photom_flux_err = Bagpipes_Phot_DF.loc[ID, flux_err_cols]

    #we are artificially inflating the flux errors for non_irac filters by 5%
    photom_flux_err[flux_err_cols] = np.sqrt(((photom_flux_err[flux_err_cols].values.astype(float))**2 + 
                                         (.05 * photom_flux[flux_cols].values.astype(float))**2))



    #############
    # We do not have IRAC for CEERS Data
    #############
    
#     flux_irac = ['IRAC_CH1_FLUX', 'IRAC_CH2_FLUX']
#     error_flux_irac = ['IRAC_CH1_FLUXERR', 'IRAC_CH2_FLUXERR']

#     photom_flux_err[error_flux_irac] = np.sqrt(((photom_flux_err[error_flux_irac].values.astype(float))**2 + 
#                                                (.2 *  photom_flux[flux_irac].values.astype(float))**2))

    #getting the snr of sources
    snr = photom_flux/photom_flux_err
    
    #if the snr is below -5 then we know it is bad we make the flux 0 and error really big
    bad_flux_idx = snr < -5

    #setting bad flux to a really small value and error to be really big
    photom_flux[bad_flux_idx] = 0
    photom_flux_err[bad_flux_idx] = 1e16

    photom = np.c_[photom_flux.astype(float), photom_flux_err.astype(float)]

    photom = photom/1000
    
    return photom

def load_phot_CANDELS(ID):
    
    ##########
    #WILL NEED TO CHANGE THIS BASED OFF OF THE INPUT CATALOG
    #WE NEED TO DO THIS 
    ##########
    Bagpipes_Phot_DF = pd.read_csv('MARISSA_CANDELS_RUBIES_PHOTOM.txt', 
                                   sep = ' ', index_col = 0)
    
    ID = int(ID)
    
    #defining the columns we will use in the photometry
    flux_cols = flux = [x for x in Bagpipes_Phot_DF.columns.values if 'F' in x and 'DF' not in x]

    flux_err_cols = [x for x in Bagpipes_Phot_DF.columns.values if 'DF' in x]

    
    non_irac_fluxes = flux_cols[:-2]
    non_irac_flux_err = flux_err_cols[:-2]
    
    #getting the full flux and flux error info
    photom_flux_non_irac = Bagpipes_Phot_DF.loc[ID, non_irac_fluxes]
    photom_flux_err_non_irac = Bagpipes_Phot_DF.loc[ID, non_irac_flux_err]

    #we are artificially inflating the flux errors for non_irac filters by 5%
    photom_flux_err_non_irac[non_irac_flux_err] = np.sqrt(((photom_flux_err_non_irac[non_irac_flux_err].values.astype(float))**2 + 
                                                 (.05 * photom_flux_non_irac[non_irac_fluxes].values.astype(float))**2))



    #############
    # We do not have IRAC for CEERS Data
    #############
    
    flux_irac = ['F36', 'F45']
    error_flux_irac = ['DF36', 'DF45']
    
    photom_flux_irac = Bagpipes_Phot_DF.loc[ID, flux_irac]
    photom_flux_err_irac = Bagpipes_Phot_DF.loc[ID, error_flux_irac]

    photom_flux_err_irac[error_flux_irac] = np.sqrt(((photom_flux_err_irac[error_flux_irac].values.astype(float))**2 + 
                                               (.2 *  photom_flux_irac[flux_irac].values.astype(float))**2))
    
    
    photom_flux = np.concatenate((photom_flux_non_irac.values, photom_flux_irac.values))
    photom_flux_err = np.concatenate((photom_flux_err_non_irac.values, photom_flux_err_irac.values))

    #getting the snr of sources
    snr = photom_flux/photom_flux_err
    
    #if the snr is below -5 then we know it is bad we make the flux 0 and error really big
    bad_flux_idx = snr < -5

    #setting bad flux to a really small value and error to be really big
    photom_flux[bad_flux_idx] = 0
    photom_flux_err[bad_flux_idx] = 1e16

    photom = np.c_[photom_flux.astype(float), photom_flux_err.astype(float)]
    
    photom = photom/1000

    return photom

def load_phot_UDS(ID):

    ##########
    #WILL NEED TO CHANGE THIS BASED OFF OF THE INPUT CATALOG
    #WE NEED TO DO THIS 
    ##########
    Bagpipes_Phot_DF = pd.read_csv('MARISSA_UDS_PHOTOM.txt', 
                                   sep = ' ', index_col = 0)
    
    ID = int(ID)

    #defining the columns we will use in the photometry
    flux_cols = flux = [x for x in Bagpipes_Phot_DF.columns.values if 'FLUX_F' in x]

    flux_err_cols = [x for x in Bagpipes_Phot_DF.columns.values if 'FLUXERR_F' in x]

    
    
    #getting the full flux and flux error info
    photom_flux = Bagpipes_Phot_DF.loc[ID, flux_cols]
    photom_flux_err = Bagpipes_Phot_DF.loc[ID, flux_err_cols]

    #we are artificially inflating the flux errors for non_irac filters by 5%
    photom_flux_err[flux_err_cols] = np.sqrt(((photom_flux_err[flux_err_cols].values.astype(float))**2 + 
                                         (.05 * photom_flux[flux_cols].values.astype(float))**2))
    

    #getting the snr of sources
    snr = photom_flux/photom_flux_err
    
    #if the snr is below -5 then we know it is bad we make the flux 0 and error really big
    bad_flux_idx = snr < -5

    #setting bad flux to a really small value and error to be really big
    photom_flux[bad_flux_idx] = 0
    photom_flux_err[bad_flux_idx] = 1e16

    photom = np.c_[photom_flux.astype(float), photom_flux_err.astype(float)]

    photom = photom/1000
    
    return photom



def fit_instruction_nebular_fixedz(z, model = 'delayed_tau'):
    
    '''
    Function that creates the bagpipes fit model, this is what bagpipes will attempt to fit with.
    
    Input:
    z (float): redshift of the source, but feel free to change this if you do not have redshift
    model (str): 
    returns:
    fit_instruction (dictionary): the intructions BP will use to fit the galaxy
    '''
    
    if model == 'delayed_tau':
        print("Making Fit Instructions for Delayed-Tau SFH Model")
        
        #Model Building 
        model = {}
        model['age'] = (.01, 13)                  # Age of the galaxy in Gyr

        model['tau'] = (.02, 14)                  # Delayed decayed time
        model["metallicity"] = (0., 2)          # Metallicity in terms of Z_sol
        model["massformed"] = (4., 13.)           # log10 of mass formed
        
        dust = {}                                 # Dust component
        dust["type"] = "Calzetti"                 # Define the shape of the attenuation curve
        dust["Av"] = (0., 3.)                     # Vary Av between 0 and 3 magnitudes

        #############
        #NOTE: Before next run need to talk to steve about fixing this at 1.44 or letting it be a param
        #      BP will fit
        #############
        dust["eta"] = 1

        #will need to include this, this includes SF regions and their emission to the spectrum
        nebular = {}

        #changed upper limits from -4 to -2 as it seems I keep getting an error with -1
        nebular["logU"] = (-4, -1)

        fit_instructions = {}
        fit_instructions['delayed'] = model
        
        fit_instructions['redshift'] = z
        
        fit_instructions['dust'] = dust
        fit_instructions['nebular'] = nebular
    
    
        return fit_instructions
    
    elif model == 'nonparam':
        #values taken from the PLANCK CMB 2018 paper
        Om0 = .315
        Ode0 = 1 - Om0
        cosmo = LambdaCDM(H0 = 67.4, 
                          Om0 = .315, 
                          Ode0 = Ode0)
        
        age_Gyr = cosmo.age(z).value
        age_Myr = age_Gyr * 1e3
        
        starting_bin = np.array([0])
        bin_end = np.log10(age_Myr) - .05
        
        bins = np.logspace(np.log10(5), bin_end, 9)
        
        age_bins = np.append(starting_bin, bins)
        
        
        print("Making Fit Instructions for Non-Parametric SFH Model")
        dust = {}
        dust["type"] = "Calzetti"
        dust["eta"] = 1.
        dust["Av"] = (0., 3.)

        nebular = {}
        nebular["logU"] = (-4, -1)

        fit_instructions = {}
        fit_instructions["dust"] = dust
        fit_instructions["nebular"] = nebular
        fit_instructions["redshift"] = z

        print(age_bins)
        
        continuity = {}
        continuity["massformed"] = (0.0001, 13.)
        continuity["metallicity"] = (0.01, 3.)
        continuity["metallicity_prior"] = "log_10"
        continuity["bin_edges"] = list(age_bins)

        for i in range(1, len(continuity["bin_edges"])-1):
            continuity["dsfr" + str(i)] = (-10., 10.)
            continuity["dsfr" + str(i) + "_prior"] = "student_t"

        fit_instructions["continuity"] = continuity
        
        return fit_instructions
        
    elif model == "bursty":
        
        #values taken from the PLANCK CMB 2018 paper
        Om0 = .315
        Ode0 = 1 - Om0
        cosmo = LambdaCDM(H0 = 67.4, 
                          Om0 = .315, 
                          Ode0 = Ode0)
        
        age_Gyr = cosmo.age(z).value
        age_Myr = age_Gyr * 1e3
        
        starting_bin = np.array([0])
        bin_end = np.log10(age_Myr) - .01
        
        bins = np.logspace(np.log10(5), bin_end, 9)
        
        age_bins = np.append(starting_bin, bins)
        
        print("Making Fit Instructions for Bursty Non-Parametric SFH Model")
        dust = {}
        dust["type"] = "Calzetti"
        dust["eta"] = 1.
        dust["Av"] = (0., 3.)

        nebular = {}
        nebular["logU"] = (-4, -1)

        fit_instructions = {}
        fit_instructions["dust"] = dust
        fit_instructions["nebular"] = nebular
        fit_instructions["redshift"] = z


        continuity = {}
        continuity["massformed"] = (0.0001, 13.)
        continuity["metallicity"] = (0.01, 5.)
        continuity["metallicity_prior"] = "log_10"
        continuity["bin_edges"] = list(age_bins)

        for i in range(1, len(continuity["bin_edges"])-1):
            continuity["dsfr" + str(i)] = (-10., 10.)
            continuity["dsfr" + str(i) + "_prior"] = "student_t"
            
            #adding this prior scale to make it bursty
            continuity["dsfr" + str(i) + "_prior_scale"] =2.0
            
        fit_instructions["continuity"] = continuity
        
        print(fit_instructions)
        
        return fit_instructions

def fit_BP(index, filters, load_func, z, run, only_fit = True, model = 'delayed_tau'):

    print('Making the BP Galaxy Model')
    BP_Galaxy = Galaxy_Model_Builder(index, load_func, filters)
    
    print('Getting the BP Fit Instructions')
    print(f'Redshift is: {z: .5f}')
    fit_instructions = fit_instruction_nebular_fixedz(z, model = model)
    
    if only_fit:
        
        start = time.time()
        fit = pipes.fit(BP_Galaxy, fit_instructions, run = run)
    
        fit.fit(verbose=True)
        end = time.time()
        
        duration = end - start
        print(f'Full Time of the Fit is: {duration:.2f} seconds, {duration/60:.2f} Minutes')
        
    else:
        
        fit = pipes.fit(BP_Galaxy, fit_instructions, run = run)
    
        fit.fit(verbose=True)
    
        return fit
    
def fit_serial_bp(DF, IDs, run,
                  load_func, 
                  filters = filters,
                  only_fit = True, 
                  test = False, model = 'nonparam'):
    
    if test:
        print('Testing the Code on the First 10 Sources')
        
        for idx in IDs[:10]:
            print(f'Fitting Galaxy ID: {idx}')
            z_tesla = get_redshift(DF, idx)
            fit_BP(idx, filters, load_func, only_fit = only_fit, model = model, z = z_tesla, run = run)
        
    else:
        
        print(f'Running on the Full Sample of: {DF.shape[0]} Sources')
        
        for idx in IDs:
            
            z_tesla = get_redshift(DF, idx)
            
            fit_BP(idx, filters, load_func, only_fit = only_fit, z = z_tesla)
        
def get_ID():
    
    if '--id' in args:
        
        indx = args.index("--id")
        
        ID = int(sys.argv[indx + 1])
    
    else:
        print('No ID detected.')
        sys.exit(-1) 
    
    return ID


if __name__ == '__main__':
    
    print('Grabbing Bagpipes Run Name')
    run = get_run_name()
    print(f'Run-Name: {run} Acquired')

    print('Attempting to read in Bagpipes Photometric Catalog')
    Bagpipes_Phot_DF = read_input_phot_cat()
    print('Read in Bagpipes Catalog')
    
    if '--test' in args:
        
        indx = args.index("--test")
        
        test = sys.argv[indx + 1]
        
        test = bool(test)
    
    else:
        print('No test detected. Defaulting to no Test')
        test = False 
        
    if survey == "CEERS":
        load_phot = load_phot_CEERS

    elif survey == "CANDELS":
        load_phot = load_phot_CANDELS

    elif survey == 'UDS':
        load_phot = load_phot_UDS
    else:
        print('No Survey detected, aborting')
        sys.exit(-1) 
    
    #ID = get_ID()
    for ID in Bagpipes_Phot_DF.index.values:
        z_tesla = get_redshift(Bagpipes_Phot_DF, ID)
        try:
            
            fit_BP(ID, filters, load_phot, z_tesla, run, only_fit = True, model = 'bursty')
            print(f'Successfully fitted ID: {ID}')
            print()
        except e:
            print('ERROR IN FITTING!!!')
        #    print(f'Check on ID: {ID}')
        #    print()
        
    