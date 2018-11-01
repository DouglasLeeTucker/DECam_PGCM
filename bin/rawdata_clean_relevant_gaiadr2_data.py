#!/usr/bin/env python
"""
    rawdata_clean_relevant_gaiadr2_data.py

    Example:
    
    rawdata_clean_relevant_gaiadr2_data.py --help

    rawdata_clean_relevant_gaiadr2_data.py --inputFile gaiadr2_new_rawdata_rawdata.csv --outputFile gaiadr2_new_y2a1_rawdata.u.csv.tmp --verbose 2

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

    status = clean_relevant_gaiadr2_data(args)

    return status


##################################
# clean_relevant_gaiadr2_data
#

def clean_relevant_gaiadr2_data(args):

    import numpy as np 
    import os
    import sys
    import datetime
    import fitsio
    import pandas as pd

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'clean_relevant_gaiadr2_data'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    inputFile = args.inputFile
    outputFile = args.outputFile
    
    # Read selected columns from inputFile...
    columns = ['RA_WRAP','RA','DEC',
               'PHOT_G_MEAN_MAG','PHOT_G_MEAN_FLUX_OVER_ERROR',
               'PHOT_BP_MEAN_MAG','PHOT_BP_MEAN_FLUX_OVER_ERROR',
               'PHOT_RP_MEAN_MAG','PHOT_RP_MEAN_FLUX_OVER_ERROR',
               'BP_RP','BP_G','G_RP','PHOT_BP_RP_EXCESS_FACTOR']

    print datetime.datetime.now()
    print """Reading in selected columns from %s...""" % (inputFile)
    df = pd.read_csv(inputFile, usecols=columns)
    print datetime.datetime.now()

    # Includes masks for mag, magerr, color, main stellar locus outliers,
    #  and BP_RP photometric excess...
    mask =  ( (1.086/df.PHOT_G_MEAN_FLUX_OVER_ERROR < 0.3) & 
              (1.086/df.PHOT_BP_MEAN_FLUX_OVER_ERROR < 0.3) & 
              (1.086/df.PHOT_RP_MEAN_FLUX_OVER_ERROR < 0.3) & 
              (df.PHOT_G_MEAN_MAG < 19.0) &  
              (df.BP_G > 0.2) & (df.BP_G < 1.6) & 
              (np.abs(df.G_RP - 0.45*(df.BP_RP + 0.2)) < 0.2) &
              (df.PHOT_BP_RP_EXCESS_FACTOR > (1.0 + 0.015*df.BP_RP*df.BP_RP)) &
              (df.PHOT_BP_RP_EXCESS_FACTOR < (1.3 + 0.060*df.BP_RP*df.BP_RP)) ) 


    # Steve Kent's Gaia DR2 -> DES transformations, of the format:
    #   des_mag = Gaia_G + intercept + slope*( (Gaia_BP-Gaia_G) - color0 ),
    #             one relation for (Gaia_BP-Gaia_G) < color0 [blue], 
    #             and another for (Gaia_BP-Gaia_G) > color0 [red].
    #
    # See S Kent's e-mail from 31 August 2018...

    skent1 = {}

    skent1['g.color0'] = 0.899
    skent1['g.intercept'] = 1.339
    skent1['g.blue.slope'] = 1.682
    skent1['g.red.slope'] = 1.015

    skent1['r.color0'] = 0.78
    skent1['r.intercept'] = -0.124
    skent1['r.blue.slope'] = -0.174
    skent1['r.red.slope'] = 0.767

    skent1['i.color0'] = 0.90
    skent1['i.intercept'] = -0.674
    skent1['i.blue.slope'] = -0.879
    skent1['i.red.slope'] = -0.437

    skent1['z.color0'] = 1.12
    skent1['z.intercept'] = -1.216
    skent1['z.blue.slope'] = -1.247
    skent1['z.red.slope'] = -0.706

    skent1['Y.color0'] = 0.91
    skent1['Y.intercept'] = -1.052
    skent1['Y.blue.slope'] = -1.441
    skent1['Y.red.slope'] = -1.028


    skent2 = {}

    skent2['g.color0'] = 0.899
    skent2['g.intercept'] = 1.349
    skent2['g.blue.slope'] = 1.702
    skent2['g.red.slope'] = 0.907

    skent2['r.color0'] = 0.78
    skent2['r.intercept'] = -0.116
    skent2['r.blue.slope'] = -0.151
    skent2['r.red.slope'] = 0.747

    skent2['i.color0'] = 0.90
    skent2['i.intercept'] = -0.691
    skent2['i.blue.slope'] = -0.925
    skent2['i.red.slope'] = -0.410

    skent2['z.color0'] = 1.12
    skent2['z.intercept'] = -1.217
    skent2['z.blue.slope'] = -1.282
    skent2['z.red.slope'] = -0.637

    skent2['Y.color0'] = 0.91
    skent2['Y.intercept'] = -1.055
    skent2['Y.blue.slope'] = -1.514
    skent2['Y.red.slope'] = -0.992


    skent3 = {}

    skent3['g.color0'] = 0.899
    skent3['g.intercept'] = 1.306
    skent3['g.blue.slope'] = 1.634
    skent3['g.red.slope'] = 0.939

    skent3['r.color0'] = 0.78
    skent3['r.intercept'] = -0.136
    skent3['r.blue.slope'] = -0.179
    skent3['r.red.slope'] = 0.747

    skent3['i.color0'] = 0.90
    skent3['i.intercept'] = -0.678
    skent3['i.blue.slope'] = -0.905
    skent3['i.red.slope'] = -0.444

    skent3['z.color0'] = 1.12
    skent3['z.intercept'] = -1.193
    skent3['z.blue.slope'] = -1.256
    skent3['z.red.slope'] = -0.873

    skent3['Y.color0'] = 0.91
    skent3['Y.intercept'] = -1.034
    skent3['Y.blue.slope'] = -1.464
    skent3['Y.red.slope'] = -1.094

    for band in ['g', 'r', 'i', 'z', 'Y']:

        # S Kent #1:

        desMagColName1 = """%sMAG_DES_1""" % (band.upper())
        color0     = """%s.color0"""     % (band)
        intercept  = """%s.intercept"""  % (band)
        blue_slope = """%s.blue.slope""" % (band)
        red_slope  = """%s.red.slope"""  % (band)

        df.loc[:,desMagColName1] = -9999.

        blueMask = (mask & (df.BP_G <= skent1[color0]))
        redMask  = (mask & (df.BP_G >  skent1[color0]))

        df.loc[blueMask,desMagColName1] = df.loc[blueMask,'PHOT_G_MEAN_MAG'] + \
             + skent1[intercept] + skent1[blue_slope]*(df.loc[blueMask,'BP_G'] - skent1[color0])

        df.loc[redMask,desMagColName1] = df.loc[redMask,'PHOT_G_MEAN_MAG'] + \
             + skent1[intercept] + skent1[red_slope]*(df.loc[redMask,'BP_G'] - skent1[color0])


        # S Kent #2:

        desMagColName2 = """%sMAG_DES_2""" % (band.upper())
        color0     = """%s.color0"""     % (band)
        intercept  = """%s.intercept"""  % (band)
        blue_slope = """%s.blue.slope""" % (band)
        red_slope  = """%s.red.slope"""  % (band)

        df.loc[:,desMagColName2] = -9999.

        blueMask = (mask & (df.BP_G <= skent2[color0]))
        redMask  = (mask & (df.BP_G >  skent2[color0]))

        df.loc[blueMask,desMagColName2] = df.loc[blueMask,'PHOT_G_MEAN_MAG'] + \
             + skent2[intercept] + skent2[blue_slope]*(df.loc[blueMask,'BP_G'] - skent1[color0])

        df.loc[redMask,desMagColName2] = df.loc[redMask,'PHOT_G_MEAN_MAG'] + \
             + skent2[intercept] + skent2[red_slope]*(df.loc[redMask,'BP_G'] - skent1[color0])

        # S Kent #3:

        desMagColName3 = """%sMAG_DES_3""" % (band.upper())
        color0     = """%s.color0"""     % (band)
        intercept  = """%s.intercept"""  % (band)
        blue_slope = """%s.blue.slope""" % (band)
        red_slope  = """%s.red.slope"""  % (band)

        df.loc[:,desMagColName3] = -9999.

        blueMask = (mask & (df.BP_G <= skent3[color0]))
        redMask  = (mask & (df.BP_G >  skent3[color0]))

        df.loc[blueMask,desMagColName3] = df.loc[blueMask,'PHOT_G_MEAN_MAG'] + \
             + skent3[intercept] + skent3[blue_slope]*(df.loc[blueMask,'BP_G'] - skent1[color0])

        df.loc[redMask,desMagColName3] = df.loc[redMask,'PHOT_G_MEAN_MAG'] + \
             + skent3[intercept] + skent3[red_slope]*(df.loc[redMask,'BP_G'] - skent1[color0])


        # S Kent average...
        desMagColName = """%sMAG_DES""" % (band.upper())
        df.loc[:,desMagColName] = ( df.loc[:,desMagColName1] + \
                                        df.loc[:,desMagColName2] + \
                                        df.loc[:,desMagColName3] ) / 3.
        

    # Output results...
    outcolumns = columns.extend(['GMAG_DES','RMAG_DES','IMAG_DES','ZMAG_DES','YMAG_DES'])
    df.to_csv(outputFile,  columns=outcolumns, index=False, float_format='%.6f')

    return 0


##################################

if __name__ == "__main__":
    main()

##################################

