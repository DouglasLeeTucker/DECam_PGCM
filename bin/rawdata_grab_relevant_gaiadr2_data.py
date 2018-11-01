#!/usr/bin/env python
"""
    rawdata_grab_relevant_gaiadr2_data.py

    Example:
    
    rawdata_grab_relevant_gaiadr2_data.py --help

    rawdata_grab_relevant_gaiadr2_data.py --inputFile inputFile.csv --outputFile outputFile.csv --verbose 2

    WARNING:  avoid using on ginormous regions.  It seems to handle SDSSDR13 coverage area well, though.
    
    """

##################################

def main():

    import argparse
    import time

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputFile', help='name of the input CSV file', default='input.csv')
    parser.add_argument('--outputFile', help='name of the output CSV file', default='output.csv')
    parser.add_argument('--inputFileRAColName', help='(case-sensitive) name of the RA column in the input CSV file', default='RA')
    parser.add_argument('--inputFileDecColName', help='(case-sensitive) name of the DEC column in the input CSV file', default='DEC')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    status = grab_relevant_gaiadr2_data(args)

    return status


##################################
# 

def grab_relevant_gaiadr2_data(args):

    import numpy as np 
    import os
    import sys
    import datetime
    import pandas as pd
    import healpixTools

    if args.verbose > 0:
        print "grab_relevant_gaiadr2_data"

    inputFile = args.inputFile
    outputFile = args.outputFile
    inputFileRAColName = args.inputFileRAColName
    inputFileDecColName = args.inputFileDecColName


    # Does input file exist for this band?
    if os.path.isfile(inputFile)==False:
        print """input file %s does not exist.  Exiting...""" % (inputFile)
        return 1


    # Read selected columns from inputFile into a pandas DataFrame...
    if args.verbose > 0:
        print """Reading in RA,DEC columns from %s as a pandas DataFrame...""" % (inputFile)
        print datetime.datetime.now()
    try:
        df = pd.read_csv(inputFile, usecols=[inputFileRAColName, inputFileDecColName])
    except:
        print """Cannot read columns %s, %s from input file %s""" % (inputFileRAColName, inputFileDecColName)
        print """Don't forget that the pandas read_csv usecols option is case-sensitive"""
        print "Exiting..."
        return 1

    if args.verbose > 0:
        print datetime.datetime.now()
        print

    #  Third, add healpix (nside=8) column and find list of unique healpix values...
    df['IPIX8'] = healpixTools.getipix(8, df[inputFileRAColName], df[inputFileDecColName])
    uniqIPIX8Array = df['IPIX8'].unique()
    
    # We no longer need df...
    del df


    # Temporary fix while working on new BLISS star cluster data...
    if True:

        myfile='/data/des40.a/data/dtucker/BLISS/StarClusterPhotom/gaiadr2_around_new_bliss_starcluster.csv'

        print 'Temporary fix while working on new BLISS star cluster data...'
        print 'Only reading in this file:'
        print '    '+myfile

        df = pd.read_csv(myfile, dtype={ 'datalink_url' : np.object_ , 'epoch_photometry_url' : np.object_ })

        # Capitalize all column names (for consistency)...
        df.columns = [x.upper() for x in df.columns]

    else:

        # This section needs updating for the Gaia DR2 catalog on disk...

        all_files = []
        for ipix8 in uniqIPIX8Array:
            myfile="""/data/des20.b/data/sallam/pyPSM_Year2/TWOMASS/ALL-2MASS/2arcsec/apass_TWO_MASS_%d.csv""" % (ipix8)
            if os.path.isfile(myfile)==False:
                print """%s does not exist.  Skipping...""" % (myfile)
            else:
                all_files.append(myfile)

        if args.verbose > 0:
            print 'Reading in columns from the various apass/2mass healpix files...'
            print datetime.datetime.now()
        # Trick from 
        #  http://stackoverflow.com/questions/20906474/import-multiple-csv-files-into-pandas-and-concatenate-into-one-dataframe
        #df = pd.concat((pd.read_csv(f) for f in all_files))
        df = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)
        if args.verbose > 0:
            print datetime.datetime.now()
            print


    df.loc[:, 'RA_WRAP'] = df.loc[:, 'RA']
    
    # If RA_WRAP > 180., subtract 360...
    mask = ( df.RA_WRAP > 180. )
    df.loc[mask,'RA_WRAP'] = df.loc[mask,'RA_WRAP'] - 360.

    # Sort by RA_WRAP...
    df.sort_values(by='RA_WRAP', ascending=True, inplace=True)

    # Write to outputFile...
    if args.verbose > 0:
        print """Writing output file %s...""" % (outputFile)
        print datetime.datetime.now()
    df.to_csv(outputFile, index=False, float_format='%.6f')
    if args.verbose > 0:
        print datetime.datetime.now()
        print

    # Delete dataframe...
    del df

    return 0



##################################

if __name__ == "__main__":
    main()

##################################
