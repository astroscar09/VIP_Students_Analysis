from Bagpipes_Fitting import *

Bagpipes_Phot_DF = read_input_phot_cat()

#ID = get_ID()
run = get_run_name()

if survey == "CEERS":
    load_phot = load_phot_CEERS

elif survey == "CANDELS":
    load_phot = load_phot_CANDELS

elif survey == 'UDS':
    load_phot = load_phot_UDS


bp_fits = {}

for ID in Bagpipes_Phot_DF.index.values:
    
    z_tesla = get_redshift(Bagpipes_Phot_DF, ID)
    fit = fit_BP(ID, filters, load_phot, z_tesla, run, only_fit = False, model = 'bursty')
    
    bp_fits[ID] = fit
    
    #print(f'Successfully fitted ID: {ID}')
    #print()
    
stellar_mass = {}
chi2 = {}

for key, fit in bp_fits.items():
    
    fit.posterior.get_advanced_quantities()
    mass = fit.posterior.samples['stellar_mass']
    chi2_source = fit.posterior.samples['chisq_phot']
    
    stellar_mass[key] = mass
    chi2[key] = chi2_source
    
    
Mass = pd.DataFrame(stellar_mass).T
Chi2 = pd.DataFrame(chi2).T

Mass.to_csv('Stellar_Mass_Marissa_RUBIES_CEERS.txt', sep = ' ')
Chi2.to_csv('Chi2_Marissa_RUBIES_CEERS.txt', sep = ' ')
