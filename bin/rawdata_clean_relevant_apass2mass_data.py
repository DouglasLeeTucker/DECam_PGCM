#!/usr/bin/env python
"""
    rawdata_clean_relevant_apass2mass_data.py

    Example:
    
    rawdata_clean_relevant_apass2mass_data.py --help

    rawdata_clean_relevant_apass2mass_data.py --inputFile apass2mass_new_rawdata_rawdata.csv --outputFile apass2mass_new_y2a1_rawdata.u.csv.tmp --verbose 2

    """

##################################

def main():

    import argparse
    import time

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputFile', help='name of the input CSV file', default='input.csv')
    parser.add_argument('--outputFile', help='name of the output CSV file', default='output.csv')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    status = clean_relevant_apass2mass_data(args)

    return status


##################################
# clean_relevant_apass2mass_data
#

def clean_relevant_apass2mass_data(args):

    import numpy as np 
    import os
    import sys
    import datetime
    import fitsio
    import pandas as pd

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'clean_relevant_apass2mass_data'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    inputFile = args.inputFile
    outputFile = args.outputFile
    
    # Read selected columns from inputFile...
    columns = ['RA_WRAP','RAJ2000_APASS','DEJ2000_APASS','GMAG_APASS','RMAG_APASS','IMAG_APASS','JMAG_2MASS','HMAG_2MASS','KMAG_2MASS','g_des','r_des','i_des','z_des','Y_des']
    print datetime.datetime.now()
    print """Reading in selected columns from %s...""" % (inputFile)
    df = pd.read_csv(inputFile, usecols=columns)
    print datetime.datetime.now()

    # Rename a couple columns...
    df.rename(columns={'RAJ2000_APASS':'RA', 'DEJ2000_APASS':'DEC'}, inplace=True)

    # Add APASS "g-r" column...
    df.loc[:,'gr_apass'] = df.loc[:,'GMAG_APASS'] - df.loc[:,'RMAG_APASS']

    # Transformation equation (20161123, from run of y2a1_new_u_from_apass2massGR.py):
    #    u_des = g_apass + 1.699*(g_apass-r_apass)**2 - 0.1106*(g_apass-r_apass) + 0.6307, 
    #  which is appropriate for "0.2 <= (g-r)_apass <= 0.8".
    #  Use a signal alue of -9999 for stars with colors outside this range....
    df.loc[:,'UMAG_DES'] = -9999.
    mask1 = ( (df.GMAG_APASS > 0.0) & (df.RMAG_APASS > 0.0) )
    mask2 = ( (df.gr_apass >=  0.2) & (df.gr_apass <= 0.8) )
    mask_all = ( mask1 & mask2 )
    df.loc[mask_all,'UMAG_DES'] = df.loc[mask_all,'GMAG_APASS'] + 1.699*df.loc[mask_all,'gr_apass']*df.loc[mask_all,'gr_apass']  - 0.1106*df.loc[mask_all,'gr_apass'] + 0.6307

    # For the time being (until we've updated the transformation equations for these 
    #  other filters), we'll keep the current values of GMAG_DES, RMAG_DES, IMAG_DES, 
    #  ZMAG_DES, and YMAG_DES
    df.loc[:,'GMAG_DES'] = df.loc[:,'g_des']
    df.loc[:,'RMAG_DES'] = df.loc[:,'r_des']
    df.loc[:,'IMAG_DES'] = df.loc[:,'i_des']
    df.loc[:,'ZMAG_DES'] = df.loc[:,'z_des']
    df.loc[:,'YMAG_DES'] = df.loc[:,'Y_des']

    # Output results...
    outcolumns = ['RA_WRAP','RA','DEC','GMAG_APASS','RMAG_APASS','IMAG_APASS','JMAG_2MASS','HMAG_2MASS','KMAG_2MASS','UMAG_DES','GMAG_DES','RMAG_DES','IMAG_DES','ZMAG_DES','YMAG_DES']
    df.to_csv(outputFile,  columns=outcolumns, index=False, float_format='%.6f')

    return 0


##################################

if __name__ == "__main__":
    main()

##################################

