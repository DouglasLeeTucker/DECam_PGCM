#!/usr/bin/env python
"""
    tie_to_fgcm_stds.py

    Example:
    
    tie_to_fgcm_stds.py --help

    tie_to_fgcm_stds.py --inputFile inputFile.csv --outputFile outputFile.csv --verbose 2
    
    """

##################################

def main():

    import argparse
    import time

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputFile', help='name of the input CSV file', default='input.csv')
    parser.add_argument('--outputFileBaseName', help='base name of the output CSV file for those methods with multiple similarly named output files', default='outputBaseName')
    parser.add_argument('--magStdBaseName', help='base name of the FGCM std mag columns (e.g., for Y3A2, it was MAG_PSF; for Y6A1, it was MAG_STD)', default='MAG_PSF')
    parser.add_argument('--magerrStdBaseName', help='base name of the FGCM std magerr columns (e.g., for Y3A2, it was MAG_PSF_ERR; for Y6A1, it was MAGERR_STD)', default='MAG_PSF_ERR')
    parser.add_argument('--band', help='name of filter band', default='g')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    status = tie_to_fgcm_stds(args)


##################################
# 

def tie_to_fgcm_stds(args):

    import numpy as np 
    import os
    import sys
    import datetime
    import pandas as pd

    inputFile = args.inputFile
    outputFileBaseName = args.outputFileBaseName
    outputFilePSF = """%s.magpsf.csv""" % (outputFileBaseName)
    outputFileAper8 = """%s.magaper8.csv""" % (outputFileBaseName)
    magStdBaseName = args.magStdBaseName
    magerrStdBaseName = args.magerrStdBaseName
    band = args.band

    # Does input file exist for this band?
    if os.path.isfile(inputFile)==False:
        print """tie_to_stds input file %s does not exist.  Exiting...""" % (inputFile)
        return 1

    # Read selected columns from inputFile into a pandas DataFrame...
    print datetime.datetime.now()
    magStd = """%s_%s_1""" % (magStdBaseName.upper(), band.upper())
    magerrStd = """%s_%s_1""" % (magerrStdBaseName.upper(), band.upper())
    print """Reading in selected columns from %s as a pandas DataFrame...""" % (inputFile)
    dataFrame = pd.read_csv(inputFile, usecols=[magStd, magerrStd, "BAND_2", "EXPNUM_2", "FILENAME_2", "PFW_ATTEMPT_ID_2", "CCDNUM_2", "RA_WRAP_2", "DEC_2", "FLUX_PSF_2", "FLUXERR_PSF_2", "FLUX_APER_8_2", "FLUXERR_APER_8_2", "AIRMASS_2", "EXPTIME_2", "SKYTILT_2", "T_EFF_2", "CALNAC_2", "MJD_OBS_2"])
    print datetime.datetime.now()
    print

    # Add a MAG_PSF_2 column and a MAG_PSF_DIFF column to the pandas DataFrame...
    dataFrame['MAG_PSF_2'] = -2.5*np.log10(dataFrame['FLUX_PSF_2'])
    dataFrame['MAG_PSF_DIFF'] = dataFrame[magStd]-dataFrame['MAG_PSF_2']

    # Add a MAG_APER_8_2 column and a MAG_APER_8_DIFF column to the pandas DataFrame...
    dataFrame['MAG_APER_8_2'] = -2.5*np.log10(dataFrame['FLUX_APER_8_2'])
    dataFrame['MAG_APER_8_DIFF'] = dataFrame[magStd]-dataFrame['MAG_APER_8_2']


    ###############################################
    # MAG_PSF:  Aggregate by FILENAME_2...
    ###############################################

    # Make a copy of original dataFrame...
    df = dataFrame.copy()

    # Create an initial mask...
    mask1 = ( (df[magStd] >= 0.) & (df[magStd] <= 25.) & (df[magerrStd] < 0.1) )
    magPsfDiffGlobalMedian = df[mask1]['MAG_PSF_DIFF'].median()
    magPsfDiffMin = magPsfDiffGlobalMedian - 5.0
    magPsfDiffMax = magPsfDiffGlobalMedian + 5.0
    mask2 = (df['MAG_PSF_DIFF'] > magPsfDiffMin) & (df['MAG_PSF_DIFF'] < magPsfDiffMax)
    mask = mask1 & mask2

    # Iterate over the copy of dataFrame 3 times, removing outliers...
    #  We are using "Method 2/Group by item" from
    #  http://nbviewer.jupyter.org/urls/bitbucket.org/hrojas/learn-pandas/raw/master/lessons/07%20-%20Lesson.ipynb
    print "Sigma-clipping..."
    niter = 0
    for i in range(3):

        niter = i + 1
        print """   iter%d...""" % ( niter )

        # make a copy of original df, and then delete the old one...
        newdf = df[mask].copy()
        del df

        # group by FILENAME_2...
        grpnewdf = newdf.groupby(['FILENAME_2'])
        
        # add/update new columns to newdf
        print datetime.datetime.now()
        newdf['Outlier']  = grpnewdf['MAG_PSF_DIFF'].transform( lambda x: abs(x-x.mean()) > 3.00*x.std() )
        print datetime.datetime.now()
        del grpnewdf
        print datetime.datetime.now()
        #print newdf

        nrows = newdf['MAG_PSF_DIFF'].size
        print """  Number of rows remaining:  %d""" % ( nrows )

        df = newdf
        mask = ( df['Outlier'] == False )  


    # Perform pandas grouping/aggregating functions on sigma-clipped Data Frame...
    print datetime.datetime.now()
    print 'Performing grouping/aggregating functions on sigma-clipped pandas DataFrame...'
    groupedDataFrame = df.groupby(['FILENAME_2'])
    magPsfZeroMedian = groupedDataFrame['MAG_PSF_DIFF'].median()
    magPsfZeroMean = groupedDataFrame['MAG_PSF_DIFF'].mean()
    magPsfZeroStd = groupedDataFrame['MAG_PSF_DIFF'].std()
    magPsfZeroNum = groupedDataFrame['MAG_PSF_DIFF'].count()
    magPsfZeroErr = magPsfZeroStd/np.sqrt(magPsfZeroNum-1)
    raMedian  = groupedDataFrame['RA_WRAP_2'].median()
    decMedian = groupedDataFrame['DEC_2'].median()
    skyTilt = groupedDataFrame['SKYTILT_2'].first()
    teff = groupedDataFrame['T_EFF_2'].first()
    calnac = groupedDataFrame['CALNAC_2'].first()
    airmass = groupedDataFrame['AIRMASS_2'].first()
    exptime = groupedDataFrame['EXPTIME_2'].first()
    expnum = groupedDataFrame['EXPNUM_2'].first()
    ccdnum = groupedDataFrame['CCDNUM_2'].first()
    pfw = groupedDataFrame['PFW_ATTEMPT_ID_2'].first()
    mjdobs = groupedDataFrame['MJD_OBS_2'].first()
    print datetime.datetime.now()
    print

    # Rename some of the pandas series...
    magPsfZeroMedian.name = 'MAG_ZERO_MEDIAN'
    magPsfZeroMean.name = 'MAG_ZERO_MEAN'
    magPsfZeroStd.name = 'MAG_ZERO_STD'
    magPsfZeroNum.name = 'MAG_ZERO_NUM'
    magPsfZeroErr.name = 'MAG_ZERO_MEAN_ERR'
    raMedian.name = 'RA_CCD'
    decMedian.name = 'DEC_CCD'
    skyTilt.name = 'SKYTILT_EXP'
    teff.name = 'T_EFF_EXP'
    calnac.name = 'CALNAC'
    airmass.name = 'AIRMASS'
    exptime.name = 'EXPTIME'
    expnum.name = 'EXPNUM'
    ccdnum.name = 'CCDNUM'
    pfw.name = 'PFW_ATTEMPT_ID'
    mjdobs.name = 'MJD_OBS'

    # Create new data frame containing all the relevant aggregate quantities
    newDataFrame = pd.concat( [expnum, ccdnum, pfw, exptime, airmass, \
                               magPsfZeroMedian, magPsfZeroMean, magPsfZeroStd, \
                               magPsfZeroErr, magPsfZeroNum, \
                               raMedian, decMedian, skyTilt, teff, calnac, mjdobs], \
                               join='outer', axis=1 )

    # Saving catname-based results to output files...
    print datetime.datetime.now()
    print """Writing %s output file (using pandas to_csv method)...""" % (outputFilePSF)
    newDataFrame.to_csv(outputFilePSF, float_format='%.4f')
    print datetime.datetime.now()
    print


    ###############################################
    # MAG_APER_8:  Aggregate by FILENAME_2...
    ###############################################

    # Make a copy of original dataFrame...
    df = dataFrame.copy()

    # Create an initial mask...
    mask1 = ( (df[magStd] >= 0.) & (df[magStd] <= 25.) & (df[magerrStd] < 0.1) )
    magAper8DiffGlobalMedian = df[mask1]['MAG_APER_8_DIFF'].median()
    magAper8DiffMin = magAper8DiffGlobalMedian - 5.0
    magAper8DiffMax = magAper8DiffGlobalMedian + 5.0
    mask2 = (df['MAG_APER_8_DIFF'] > magAper8DiffMin) & (df['MAG_APER_8_DIFF'] < magAper8DiffMax)
    mask = mask1 & mask2

    # Iterate over the copy of dataFrame 3 times, removing outliers...
    #  We are using "Method 2/Group by item" from
    #  http://nbviewer.jupyter.org/urls/bitbucket.org/hrojas/learn-pandas/raw/master/lessons/07%20-%20Lesson.ipynb
    print "Sigma-clipping..."
    niter = 0
    for i in range(3):

        niter = i + 1
        print """   iter%d...""" % ( niter )

        # make a copy of original df, and then delete the old one...
        newdf = df[mask].copy()
        del df

        # group by FILENAME_2...
        grpnewdf = newdf.groupby(['FILENAME_2'])
        
        # add/update new columns to newdf
        print datetime.datetime.now()
        newdf['Outlier']  = grpnewdf['MAG_APER_8_DIFF'].transform( lambda x: abs(x-x.mean()) > 3.00*x.std() )
        print datetime.datetime.now()
        del grpnewdf
        print datetime.datetime.now()

        nrows = newdf['MAG_APER_8_DIFF'].size
        print """  Number of rows remaining:  %d""" % ( nrows )

        df = newdf
        mask = ( df['Outlier'] == False )  


    # Perform pandas grouping/aggregating functions on sigma-clipped Data Frame...
    print datetime.datetime.now()
    print 'Performing grouping/aggregating functions on sigma-clipped pandas DataFrame...'
    groupedDataFrame = df.groupby(['FILENAME_2'])
    magAper8ZeroMedian = groupedDataFrame['MAG_APER_8_DIFF'].median()
    magAper8ZeroMean = groupedDataFrame['MAG_APER_8_DIFF'].mean()
    magAper8ZeroStd = groupedDataFrame['MAG_APER_8_DIFF'].std()
    magAper8ZeroNum = groupedDataFrame['MAG_APER_8_DIFF'].count()
    magAper8ZeroErr = magAper8ZeroStd/np.sqrt(magAper8ZeroNum-1)
    raMedian  = groupedDataFrame['RA_WRAP_2'].median()
    decMedian = groupedDataFrame['DEC_2'].median()
    skyTilt = groupedDataFrame['SKYTILT_2'].first()
    teff = groupedDataFrame['T_EFF_2'].first()
    calnac = groupedDataFrame['CALNAC_2'].first()
    airmass = groupedDataFrame['AIRMASS_2'].first()
    exptime = groupedDataFrame['EXPTIME_2'].first()
    expnum = groupedDataFrame['EXPNUM_2'].first()
    ccdnum = groupedDataFrame['CCDNUM_2'].first()
    pfw = groupedDataFrame['PFW_ATTEMPT_ID_2'].first()
    mjdobs = groupedDataFrame['MJD_OBS_2'].first()
    print datetime.datetime.now()
    print

    # Rename some of the pandas series...
    magAper8ZeroMedian.name = 'MAG_ZERO_MEDIAN'
    magAper8ZeroMean.name = 'MAG_ZERO_MEAN'
    magAper8ZeroStd.name = 'MAG_ZERO_STD'
    magAper8ZeroNum.name = 'MAG_ZERO_NUM'
    magAper8ZeroErr.name = 'MAG_ZERO_MEAN_ERR'
    raMedian.name = 'RA_CCD'
    decMedian.name = 'DEC_CCD'
    skyTilt.name = 'SKYTILT_EXP'
    teff.name = 'T_EFF_EXP'
    calnac.name = 'CALNAC'
    airmass.name = 'AIRMASS'
    exptime.name = 'EXPTIME'
    expnum.name = 'EXPNUM'
    ccdnum.name = 'CCDNUM'
    pfw.name = 'PFW_ATTEMPT_ID'
    mjdobs.name = 'MJD_OBS'

    # Create new data frame containing all the relevant aggregate quantities
    newDataFrame = pd.concat( [expnum, ccdnum, pfw, exptime, airmass, \
                               magAper8ZeroMedian, magAper8ZeroMean, magAper8ZeroStd, \
                               magAper8ZeroErr, magAper8ZeroNum, \
                               raMedian, decMedian, skyTilt, teff, calnac, mjdobs], \
                               join='outer', axis=1 )

    # Saving catname-based results to output files...
    print datetime.datetime.now()
    print """Writing %s output file (using pandas to_csv method)...""" % (outputFileAper8)
    newDataFrame.to_csv(outputFileAper8, float_format='%.4f')
    print datetime.datetime.now()
    print


    # Cleaning up to save memory...
    del dataFrame
    del newDataFrame
    del magPsfZeroMedian
    del magPsfZeroMean
    del magPsfZeroStd
    del magPsfZeroNum
    del magPsfZeroErr
    del magAper8ZeroMedian
    del magAper8ZeroMean
    del magAper8ZeroStd
    del magAper8ZeroNum
    del magAper8ZeroErr
    del raMedian
    del decMedian
    del skyTilt
    del teff
    del calnac
    del airmass
    del exptime
    del expnum
    del ccdnum
    del pfw
    del mjdobs

    return 0


##################################

if __name__ == "__main__":
    main()

##################################

