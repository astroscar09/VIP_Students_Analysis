import sys
args = list(map(str.lower,sys.argv))
import pandas as pd

def GetIndex():
    
    if '--index' in args:
        
        indx = args.index("--index")
        
        gal_index = sys.argv[indx + 1]
        
        return gal_index
    
    else:
        print('No Galaxy Index detected. Use --index <index>')
        sys.exit(-1) 
        
def read_input_phot_cat():
    
    if '--bp_input_cat' in args:
        
        indx = args.index("--bp_input_cat")
        
        bp_phot_path = sys.argv[indx + 1]
        
        bp_phot_cat = pd.read_csv(bp_phot_path, 
                                  sep = ' ', 
                                  index_col = 0)
        
        return bp_phot_cat
    
    else:
        print('No Bagpipes input catalog. Use --bp_input_cat <path to BP phot Cat>')
        sys.exit(-1) 
        
def get_run_name():
    
    if '--run_name' in args:
        
        indx = args.index("--run_name")
        
        run_name = sys.argv[indx + 1]
        
        return run_name
    
    else:
        print('No Bagpipes input catalog. Use --bp_input_cat <path to BP phot Cat>')
        sys.exit(-1) 

def make_filters(filter_set = 'TESLA'):
    
    if filter_set == 'TESLA':
        
        print('TESLA Filter Set to be used \n')
        
        base_TESLA = '/Users/oac466/Desktop/filters/'

        TESLA_Filters = ['CFHT/H20_CFHT_Megaprime.u.dat', 
                         'HSC_Filters/HSC_g_filt.txt',
                         'HSC_Filters/HSC_r_filt.txt',
                         'HSC_Filters/HSC_i_filt.txt',
                         'HSC_Filters/HSC_z_filt.txt',
                         'HSC_Filters/HSC_y_filt.txt',
                         'HSC_Filters/IRAC_36_filt.txt',
                         'HSC_Filters/IRAC_45_filt.txt']    

        TESLA_filts = [base_TESLA+x for x in TESLA_Filters]
        
        print('READING IN TESLA FILTERS: ')
        for x in TESLA_filts:
            print(x)
        print()
        return TESLA_filts
    
    elif filter_set == 'SHELA':
        
        print('SHELA Filter Set to be used \n')
        
        base_SHELA = '/Users/oac466/Desktop/filters/'
        
#         TESLA_Filters = ['CFHT/H20_CFHT_Megaprime.u.dat', 
#                          'HSC_Filters/HSC_g_filt.txt',
#                          'HSC_Filters/HSC_r_filt.txt',
#                          'HSC_Filters/HSC_i_filt.txt',
#                          'HSC_Filters/HSC_z_filt.txt',
#                          'HSC_Filters/HSC_y_filt.txt',
#                          'HSC_Filters/IRAC_36_filt.txt',
#                          'HSC_Filters/IRAC_45_filt.txt']    

#         TESLA_filts = [base_TESLA+x for x in TESLA_Filters]

#         return TESLA_filts
        pass

    elif filter_set == 'CEERS':
        '''
        The Filters used in the CEERS program observations are listed below, got this form the UNICORN Catalog
        #          'FLUX_F606W',
        #          'FLUX_F814W',
        #          'FLUX_F105W',
        #          'FLUX_F125W',
        #          'FLUX_F140W',
        #          'FLUX_F160W',
        #          'FLUX_F115W',
        #          'FLUX_F150W',
        #          'FLUX_F200W',
        #          'FLUX_F277W',
        #          'FLUX_F356W',
        #          'FLUX_F410M',
        #          'FLUX_F444W'
        '''
    
        base_ceers = '/Users/oac466/Desktop/filters/'
        
        filter_files =  ['HST/ACS/ACS_F435W.txt',
                         'HST/ACS/ACS_F606W.txt', 
                         'HST/ACS/ACS_F814W.txt',
                         'HST/WFC3/WFC3_F105W.txt',
                         'HST/WFC3/WFC3_F125W.txt',
                         'HST/WFC3/WFC3_F140W.txt',
                         'HST/WFC3/WFC3_F160W.txt',
                         'JWST/F090W.txt',
                         'JWST/F115W.txt',
                         'JWST/F150W.txt', 
                         'JWST/F200W.txt', 
                         'JWST/F277W.txt', 
                         'JWST/F356W.txt', 
                         'JWST/F410M.txt', 
                         'JWST/F444W.txt', 
                         'JWST/F470N.txt']
        
        CEERS_filts = [base_ceers+x for x in filter_files]
        
        print('READING IN the CEERS FILTERS: ')
        for x in CEERS_filts:
            print(x)
        print('-------------------------------------------------------------------------------------')
        print()
        
    
        return CEERS_filts
    
    elif filter_set == 'CANDELS':
        
        #'F606W', 'F814W', 'F105W', 'F125W', 'F140W', 'F160W', 'F36', 'F45'
        
        base_candels = '/Users/oac466/Desktop/filters/'
        filter_files =  ['HST/ACS/ACS_F606W.txt', 
                         'HST/ACS/ACS_F814W.txt',
                         'HST/WFC3/WFC3_F105W.txt',
                         'HST/WFC3/WFC3_F125W.txt',
                         'HST/WFC3/WFC3_F140W.txt',
                         'HST/WFC3/WFC3_F160W.txt',
                         'Spitzer/IRAC/IRAC_36.txt', 
                         'Spitzer/IRAC/IRAC_45.txt']
        
        CANDELS_filts = [base_candels+x for x in filter_files]
        print('READING IN the CANDELS FILTERS: ')
        for x in CANDELS_filts:
            print(x)
        print('-------------------------------------------------------------------------------------')
        print()
        
        return CANDELS_filts
    
    elif filter_set == 'UDS':

        # 'FLUX_F435W',
        # 'FLUX_F606W',
        # 'FLUX_F814W',
        # 'FLUX_F125W',
        # 'FLUX_F140W',
        # 'FLUX_F160W',
        # 'FLUX_F090W',
        # 'FLUX_F115W',
        # 'FLUX_F150W',
        # 'FLUX_F200W',
        # 'FLUX_F277W',
        # 'FLUX_F356W',
        # 'FLUX_F410M',
        # 'FLUX_F444W'

        base_uds = '/Users/oac466/Desktop/filters/'
        filter_files =  ['HST/ACS/ACS_F435W.txt',
                        'HST/ACS/ACS_F606W.txt', 
                        'HST/ACS/ACS_F814W.txt',
                        'HST/WFC3/WFC3_F125W.txt',
                        'HST/WFC3/WFC3_F140W.txt',
                        'HST/WFC3/WFC3_F160W.txt',
                        'JWST/F090W.txt',
                        'JWST/F115W.txt',
                        'JWST/F150W.txt',
                        'JWST/F200W.txt',
                        'JWST/F277W.txt',
                        'JWST/F356W.txt',
                        'JWST/F410M.txt',
                        'JWST/F444W.txt']
        
        UDS_filts = [base_uds+x for x in filter_files]

        print('Reading in the UDS Filters: ')
        for x in UDS_filts:
            print(x)
        print('-------------------------------------------------------------------------------------')
        print()

        return UDS_filts

    
    else:
        print('No Valid Filter Set Given')
        print('Only have TESLA, SHELA, CEERS')
        sys.exit(-1)