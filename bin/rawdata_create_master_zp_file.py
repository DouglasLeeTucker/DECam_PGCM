#!/usr/bin/env python
"""
    NO UPDATES SINCE Y3A1 u-band calibrations (December 2016)!!!

    y2a1_create_master_zp_file_2.py

    Example:
    
    y2a1_create_master_zp_file_2.py --help

    y2a1_create_master_zp_file_2.py 
         --inputRawDataFile y2a1_rawdata.sorted.ccd_radec.u.csv
         --inputPSMCCDZPFile zps.psm_y2a1_rawdata.sorted.ccd_radec.u.magpsf.good_a.csv
         --inputDR13CCDZPFile zps.pgcm_forced.new_sec.y2a1_rawdata.ccd_radec.u.magpsf.ccd.csv
         --inputDR13EXPZPFile zps.pgcm_forced.new_sec.y2a1_rawdata.ccd_radec.u.magpsf.exp.csv 
         --inputAPASSCCDZPFile zps.pgcm_forced.apass2mass_new_fixed_match_y2a1_rawdata.sorted.ccd_radec.u.magpsf.ccd.csv
         --inputAPASSEXPZPFile zps.pgcm_forced.apass2mass_new_fixed_match_y2a1_rawdata.sorted.ccd_radec.u.magpsf.exp.csv
         --outputFile zps.pgcm_sdssd13_apass_pypsm.u.magpsf.master.csv 
         --verbose 2
    
    """

##################################

def main():

    import argparse
    import time

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputRawDataFile', help='name of the input rawdata CSV file', default='inputRawData.csv')
    parser.add_argument('--inputPSMCCDZPFile', help='name of the input ccd-level pyPSM-based zp CSV file', default='inputPSMCCDZP.csv')
    parser.add_argument('--inputDR13CCDZPFile', help='name of the input ccd-level DES-transformed SDSSDR13-based zp CSV file', default='inputDR13CCDZP.csv')
    parser.add_argument('--inputDR13EXPZPFile', help='name of the input exposure-level DES-transformed SDSSDR13-based zp CSV file', default='inputDR13EXPZP.csv')
    parser.add_argument('--inputAPASSCCDZPFile', help='name of the input ccd-level DES-transformed APASS/2MASS-based zp CSV file', default='inputAPASSCCDZP.csv')
    parser.add_argument('--inputAPASSEXPZPFile', help='name of the input exposure-level DES-transformed APASS/2MASS-based zp CSV file', default='inputAPASSEXPZP.csv')
    parser.add_argument('--outputFile', help='name of the output CSV file', default='output.csv')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    status = create_master_zp_file_2(args)
        

##################################
# 

def create_master_zp_file_2(args):

    import numpy as np
    import pandas as pd
    import math
    import datetime
    import decimal

    print 'Start of program: ',
    print datetime.datetime.now()

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'create_master_zp_file_2'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    inputRawDataFile    = args.inputRawDataFile
    inputPSMCCDZPFile   = args.inputPSMCCDZPFile
    inputDR13CCDZPFile  = args.inputDR13CCDZPFile
    inputDR13EXPZPFile  = args.inputDR13EXPZPFile
    inputAPASSCCDZPFile = args.inputAPASSCCDZPFile
    inputAPASSEXPZPFile = args.inputAPASSEXPZPFile
    outputFile          = args.outputFile


    # Read in input rawdata file and create a master list of 
    #  zeropoints for each filename (expnum+ccdnum)...
    print """Reading %s into a data frame...""" % (inputRawDataFile)
    df = pd.read_csv(inputRawDataFile, usecols=['FILENAME','EXPNUM','CCDNUM','PFW_ATTEMPT_ID','EXPTIME','AIRMASS','EXPRA','EXPDEC','RA_CCD','DEC_CCD','SKYTILT','T_EFF','CALNAC','MJD_OBS'])
    
    # Group by filename...
    df_grouped = df.groupby(['FILENAME'])

    # Extract relevant quantities...
    expnum = df_grouped['EXPNUM'].first()
    ccdnum = df_grouped['CCDNUM'].first()
    pfw = df_grouped['PFW_ATTEMPT_ID'].first()
    exptime = df_grouped['EXPTIME'].first()
    airmass = df_grouped['AIRMASS'].first()
    expra = df_grouped['EXPRA'].first()
    expdec = df_grouped['EXPDEC'].first()
    ccdra = df_grouped['RA_CCD'].first()
    ccddec = df_grouped['DEC_CCD'].first()
    skyTilt = df_grouped['SKYTILT'].first()
    teff = df_grouped['T_EFF'].first()
    calnac = df_grouped['CALNAC'].first()
    mjdobs = df_grouped['MJD_OBS'].first()

    # Ensure that the columns have the appropriate names...
    expnum.name = 'EXPRA'
    expnum.name = 'EXPNUM'
    ccdnum.name = 'CCDNUM'
    pfw.name = 'PFW_ATTEMPT_ID'
    exptime.name = 'EXPTIME'
    airmass.name = 'AIRMASS'
    expra.name = 'EXPRA'
    expdec.name = 'EXPDEC'
    ccdra.name = 'RA_CCD'
    ccddec.name = 'DEC_CCD'
    skyTilt.name = 'SKYTILT_EXP'
    teff.name = 'T_EFF_EXP'
    calnac.name = 'CALNAC'
    mjdobs.name = 'MJD_OBS'

    # Create master data frame of the appropriate quantities...
    df_master = pd.concat( [expnum, ccdnum, pfw, exptime, airmass, ccdra, ccddec, skyTilt, teff, calnac, mjdobs], join='outer', axis=1 )

    # Sort master data frame by EXPNUM, CCDNUM...
    df_master.sort(['EXPNUM','CCDNUM'], inplace=True)

    # Delete old data frame and grouped data frame, to save memory...
    del df_grouped
    del df

    # Add and initialize values for columns for MAG_ZERO_MEDIAN, MAG_ZERO_MEAN, 
    #  MAG_ZERO_STD, MAG_ZERO_MEAN_ERR, and MAG_ZERO_NUM to master data frame...
    magZeroColList = ['MAG_ZERO_MEDIAN','MAG_ZERO_MEAN','MAG_ZERO_STD','MAG_ZERO_MEAN_ERR','MAG_ZERO_NUM']
    for magZeroCol in magZeroColList:
        df_master[magZeroCol] = np.nan
    # Reset the MAG_ZER_NUM column to 0...
    df_master.loc[:, 'MAG_ZERO_NUM'] = 0
    # Add a calibration source flag, intializing to a default of "not calibrated"
    # (0=not calibrated, 1=PGM_FORCED_EXP, 2=PGCM_FORCED_CCD, ...)...
    df_master.loc[:, 'CALIB_SRC_FLAG'] = 0

    # Read input PSM-based CCD ZP file into a data frame, using the filename as the index
    print """Reading %s into a data frame...""" % (inputPSMCCDZPFile)
    df_psmzpccd = pd.read_csv(inputPSMCCDZPFile, index_col='FILENAME')
    # ... and add new columns based on old columns...
    df_psmzpccd.loc[:,'MAG_ZERO_MEDIAN']   = df_psmzpccd.loc[:,'MAG_ZERO_C26202']
    df_psmzpccd.loc[:,'MAG_ZERO_MEAN']     = df_psmzpccd.loc[:,'MAG_ZERO_C26202']
    df_psmzpccd.loc[:,'MAG_ZERO_STD']      = df_psmzpccd.loc[:,'MAG_ZERO_C26202_ERR']
    df_psmzpccd.loc[:,'MAG_ZERO_MEAN_ERR'] = df_psmzpccd.loc[:,'MAG_ZERO_C26202_ERR']
    df_psmzpccd.loc[:,'MAG_ZERO_NUM']      = -1
    # ... and add the median MAG_ZERO_C26202 per exposure as a column (for later comparisons)...
    grpdf_psmzpccd = df_psmzpccd.groupby(['EXPNUM'])
    df_psmzpccd['MAG_ZERO_C26202_MEDIAN_EXP'] = grpdf_psmzpccd['MAG_ZERO_C26202'].transform( lambda x: x.median() )
    # ... and drop any rows for which MAG_ZERO_C26202 or MAG_ZERO_C26202_MEDIAN_EXP are NaN...
    df_psmzpccd.dropna(subset=['MAG_ZERO_C26202','MAG_ZERO_C26202_MEDIAN_EXP'],inplace=True)

    # Read input SDSSDR13-based CCD ZP file into a data frame, using the filename as the index
    print """Reading %s into a data frame...""" % (inputDR13CCDZPFile)
    df_dr13zpccd = pd.read_csv(inputDR13CCDZPFile, index_col='FILENAME_2')
    # ... and add 0.05mag RMS in quadrature to MAG_ZERO_MEAN_ERR...
    df_dr13zpccd.loc[:,'MAG_ZERO_MEAN_ERR'] = \
        np.sqrt(df_dr13zpccd.loc[:,'MAG_ZERO_MEAN_ERR']*df_dr13zpccd.loc[:,'MAG_ZERO_MEAN_ERR'] + 0.05*0.05)
    # ... and drop any rows for which MAG_ZERO_MEDIAN is NaN...
    df_dr13zpccd.dropna(subset=['MAG_ZERO_MEDIAN'],inplace=True)

    # Read input SDSSDR13-based EXP ZP file into a data frame, using the expnum as the index
    print """Reading %s into a data frame...""" % (inputDR13EXPZPFile)
    df_dr13zpexp = pd.read_csv(inputDR13EXPZPFile, index_col='EXPNUM_2')
    # ... and add 0.05mag RMS in quadrature to MAG_ZERO_MEAN_ERR...
    df_dr13zpexp.loc[:,'MAG_ZERO_MEAN_ERR'] = \
        np.sqrt(df_dr13zpexp.loc[:,'MAG_ZERO_MEAN_ERR']*df_dr13zpexp.loc[:,'MAG_ZERO_MEAN_ERR'] + 0.05*0.05)
    # ... and drop any rows for which MAG_ZERO_MEDIAN is NaN...
    df_dr13zpexp.dropna(subset=['MAG_ZERO_MEDIAN'],inplace=True)

    # Read input APASS/2MASS-based CCD ZP file into a data frame, using the filename as the index
    print """Reading %s into a data frame...""" % (inputAPASSCCDZPFile)
    df_apasszpccd = pd.read_csv(inputAPASSCCDZPFile, index_col='FILENAME_2')
    # ... and add 0.10mag RMS in quadrature to MAG_ZERO_MEAN_ERR...
    df_apasszpccd.loc[:,'MAG_ZERO_MEAN_ERR'] = \
        np.sqrt(df_apasszpccd.loc[:,'MAG_ZERO_MEAN_ERR']*df_apasszpccd.loc[:,'MAG_ZERO_MEAN_ERR'] + 0.10*0.10)
    # ... and drop any rows for which MAG_ZERO_MEDIAN is NaN...
    df_apasszpccd.dropna(subset=['MAG_ZERO_MEDIAN'],inplace=True)

    # Read input APASS/2MASS-based EXP ZP file into a data frame, using the expnum as the index
    print """Reading %s into a data frame...""" % (inputAPASSEXPZPFile)
    df_apasszpexp = pd.read_csv(inputAPASSEXPZPFile, index_col='EXPNUM_2')
    # ... and add 0.10mag RMS in quadrature to MAG_ZERO_MEAN_ERR...
    df_apasszpexp.loc[:,'MAG_ZERO_MEAN_ERR'] = \
        np.sqrt(df_apasszpexp.loc[:,'MAG_ZERO_MEAN_ERR']*df_apasszpexp.loc[:,'MAG_ZERO_MEAN_ERR'] + 0.10*0.10)
    # ... and drop any rows for which MAG_ZERO_MEDIAN is NaN...
    df_apasszpexp.dropna(subset=['MAG_ZERO_MEDIAN'],inplace=True)

    # Create lists of indices for the 3 currently existing data frames...
    masterfilenameList   = df_master.index.tolist()
    psmccdfilenameList   = df_psmzpccd.index.tolist()
    dr13ccdfilenameList  = df_dr13zpccd.index.tolist()
    dr13expexpnumList    = df_dr13zpexp.index.tolist()
    apassccdfilenameList = df_apasszpccd.index.tolist()
    apassexpexpnumList   = df_apasszpexp.index.tolist()

    # Loop through the master data frame, updating 
    #  the zeropoints as appropriate...
    print 'Updating master data frame...'
    for filename in masterfilenameList:

        expnum = df_master.at[filename,'EXPNUM']

        if args.verbose > 1:
            print filename, expnum

        # Determine which -- if any -- calibration to use...
        #  (working our way from worst to best...)

        # Default flag is 0 ("no calibration")...
        flag = 0
        # First, check if there is an apass/2mass exposure-level calibration...
        #  This serves as a basic calibration (flag=1)...
        if expnum in apassexpexpnumList:
            if df_apasszpexp.at[expnum,'MAG_ZERO_NUM'] > 2:
                flag = 1
        # Second, check if there is an apass/2mass ccd-level calibration...
        #  This serves as an improvement over an apass/2mass exposure-level 
        #  calibration (flag=2)...
        if filename in apassccdfilenameList:
            if df_apasszpccd.at[filename,'MAG_ZERO_NUM'] > 2:
                flag = 2 
        # Third, check if there is an SDSSDR13-based exposure-level calibration...
        #  This serves as an improvement over an apass/2mass calibration (flag=3)...
        if expnum in dr13expexpnumList:
            if df_dr13zpexp.at[expnum,'MAG_ZERO_NUM'] > 2:
                flag = 3
        # Fourth, check if there is an SDSSDR13-based ccd-level calibration...
        #  This serves as an improvement over an SDSSDR13-based exposure-level 
        #  calibration (flag=4)...
        if filename in dr13ccdfilenameList:
            if df_dr13zpccd.at[filename,'MAG_ZERO_NUM'] > 2:
                flag = 4 
        # Fifth, check if there is a pyPSM ccd-level calibration...
        #  This serves as an improvement over an SDSSDR13-based calibration, 
        #  BUT ONLY IF THE EXPOSURE WAS OBSERVED UNDER PHOTOMETRIC CONDITIONS 
        #  (INCLUDING NO DOME OCCLUSIONS) (flag=5)...
        if filename in psmccdfilenameList:
            psm_mag_zero_exp = df_psmzpccd.at[filename,'MAG_ZERO_C26202_MEDIAN_EXP']
            if ( (flag==1) or (flag==2) ):
                apass_mag_zero_exp = df_apasszpexp.at[expnum,'MAG_ZERO_MEDIAN']
                dmag_zero_exp = psm_mag_zero_exp - apass_mag_zero_exp 
                if (abs(dmag_zero_exp) <= 0.10):  flag = 5
            elif flag==3:
                dr13_mag_zero_exp = df_dr13zpexp.at[expnum,'MAG_ZERO_MEDIAN']
                dmag_zero_exp = psm_mag_zero_exp - dr13_mag_zero_exp 
                if (abs(dmag_zero_exp) <= 0.05):  flag = 5



        # Now assign the chosen MAG_ZERO values for this filename in the
        #  df_master dataframe...

        #  ...  set to APASS/2MASS-based exposure-level MAG_ZERO...
        if flag==1:
            df_master.at[filename,'CALIB_SRC_FLAG'] = flag
            for magZeroCol in magZeroColList:
                df_master.at[filename,magZeroCol] = df_apasszpexp.at[expnum,magZeroCol]

        #  ...  set to APASS/2MASS-based ccd-level MAG_ZERO...
        elif flag==2:
            df_master.at[filename,'CALIB_SRC_FLAG'] = flag
            for magZeroCol in magZeroColList:
                df_master.at[filename,magZeroCol] = df_apasszpccd.at[filename,magZeroCol]

        #  ...  set to SDSSDR13-based exposure-level MAG_ZERO...
        elif flag==3:
            df_master.at[filename,'CALIB_SRC_FLAG'] = flag
            for magZeroCol in magZeroColList:
                df_master.at[filename,magZeroCol] = df_dr13zpexp.at[expnum,magZeroCol]

        #  ...  set to SDSSDR13-based ccd-level MAG_ZERO...
        elif flag==4:
            df_master.at[filename,'CALIB_SRC_FLAG'] = flag
            for magZeroCol in magZeroColList:
                df_master.at[filename,magZeroCol] = df_dr13zpccd.at[filename,magZeroCol]

        #  ...  set to pyPSM-based exposure-level MAG_ZERO...
        elif flag==5:
            df_master.at[filename,'CALIB_SRC_FLAG'] = flag
            for magZeroCol in magZeroColList:
                df_master.at[filename,magZeroCol] = df_psmzpccd.at[filename,magZeroCol]

        # otherwise, leave "as is" and continue to the next 
        #  entry in the master data frame...
        else:
            continue


    # Writing master data frame to output file...
    print """Writing master data frame to %s""" % (outputFile)
    outputCols = ['EXPNUM','CCDNUM','PFW_ATTEMPT_ID','EXPTIME','AIRMASS','MAG_ZERO_MEDIAN','MAG_ZERO_MEAN','MAG_ZERO_STD','MAG_ZERO_MEAN_ERR','MAG_ZERO_NUM','RA_CCD','DEC_CCD','SKYTILT_EXP','T_EFF_EXP','CALNAC','MJD_OBS','CALIB_SRC_FLAG']
    df_master.to_csv(outputFile, columns=outputCols, float_format='%.4f')


    return 0


##################################

if __name__ == "__main__":
    main()

##################################





