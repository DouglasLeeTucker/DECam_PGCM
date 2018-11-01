#!/usr/bin/env python
"""
    sdss_galex_transform_u.py

    Example:
    
    sdss_galex_transform_u --help

    sdss_galex_transform_u.py --matchFile sdssdr13GUVCat_297.csv --verbose 2

    """

##################################

import numpy as np 
import pandas as pd
import math
from scipy import interpolate
from scipy.optimize import leastsq
import os
import sys
import glob
import datetime
import extinction
import matplotlib.pyplot as plt
import plotly
from plotly.offline import download_plotlyjs, plot, iplot
import plotly.graph_objs as go


import healpy as hp
import healpixTools

import paramFile
    

def main():

    import argparse
    import time

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--matchFile', help='name of the input CSV match file', default='default')
    parser.add_argument('--resultsFile', help='name of the fitting results output CSV file', default='default')
    parser.add_argument('--planFile', help='name of the input plan file', default='sdss_galex_transform_u.par')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    status = sdss_galex_transform_u(args)

    return status


##################################
# sdss_galex_transform_u
#

def sdss_galex_transform_u(args):

    import os
    import sys
    import paramFile

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'sdss_galex_transform_u'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 


    # Read contents of planFile into python dictionary...
    planDict = paramFile.readParamFile(args.planFile, args.verbose)

    # Check that all the required keywords were found...
    requiredKeywordList = ['norder',
                           'nsigma',
                           'niter',
                           'matchDir',
                           'matchFile',
                           'resultsFile']

    flag = 0
    for requiredKeyword in requiredKeywordList:
        if requiredKeyword not in planDict:
            print """Required keyword '%s' not in planFile %s""" % (requiredKeyword, args.planFile)
            flag = flag + 1
    if flag > 0:
        print 'Returning with error code 1 now...'
        return 1


    # Extract the relevant info from the planDict...
    
    #  Extract norder...
    #  (Default is 2)
    if planDict['norder'].lower() == 'default':
        norder = 2
    else:
        norder = int(planDict['norder'])    
    if args.verbose > 0:
        print 'norder: ', norder

    #  Extract nsigma...
    #  (Default is 3.0)
    if planDict['nsigma'].lower() == 'default':
        nsigma = 3.0
    else:
        nsigma = float(planDict['nsigma'])    
    if args.verbose > 0:
        print 'nsigma: ', nsigma

    #  Extract niter...
    #  (Default is 3)
    if planDict['niter'].lower() == 'default':
        niter = 3
    else:
        niter = int(planDict['niter'])    
    if args.verbose > 0:
        print 'niter: ', norder

    #  Grab matchFile name from command-line argument list....
    #  If it is not set, then extract the full path 
    #  (matchDir+matchFile) from the paramFile...
    matchFile = args.matchFile
    if matchFile == 'default':
        #  Extract the name of the matchDir...
        matchDir = planDict['matchDir']
        if matchDir.lower() == 'default':
            matchDir = '/data/des20.b/data/sallam/GALEX/GUVCat_AIS_Bianchietal2017/DLT-03-07-18/sdssdr13GUVCat'
        #  Extract the name of the matchFile...
        matchFile = planDict['matchFile']
        if matchFile.lower() == 'default':
            matchFile = 'sdssdr13GUVCat_297.csv'
        #  Create full patch of the matchFile
        matchFile = os.path.join(matchDir,matchFile)

    # Check to make sure matchFile exists...
    if os.path.isfile(matchFile)==False:
        print """matchFile %s does not exist...""" % (matchFile)
        print 'Returning with error code 1 now...'
        return 1
    if args.verbose > 0:
        print 'matchFile: ', matchFile

    #  Grab resultsFile name from command-line argument list....
    #  If it is not set, then extract the full path
    #  (resultsDir+resultsFile) from the paramFile...
    resultsFile = args.resultsFile
    if resultsFile == 'default':
        #  Extract the name of the resultsDir...
        resultsDir = planDict['resultsDir']
        if resultsDir.lower() == 'default':
            resultsDir = os.getcwd()
        #  Extract the name of the matchFile...
        resultsFile = planDict['resultsFile']
        if resultsFile.lower() == 'default':
            resultsFile = 'sdss_galex_transform_results.csv'
        #  Create full patch of the matchFile
        resultsFile = os.path.join(resultsDir,resultsFile)


    # Create a full path base name for QA files based on the value of resultsFile...
    qaFileBaseName =  os.path.splitext(resultsFile)[0]


    # Read in selected columns from file into Pandas dataframe:
    columns = ['psfMag_u','psfMag_g','psfMag_r','psfMag_i','psfMag_z','FUV_MAG','NUV_MAG']
    df = pd.read_csv(matchFile, usecols=columns)
    df.head(10)


    # Rename columns...
    df.rename(columns={'psfMag_u':'u_sdss',
                       'psfMag_g':'g_sdss',
                       'psfMag_r':'r_sdss',
                       'psfMag_i':'i_sdss',
                       'psfMag_z':'z_sdss',
                       'FUV_MAG':'FUV_galex',
                       'NUV_MAG':'NUV_galex'
                       },inplace=True)


    # Add color columns...
    df.loc[:,'ug_sdss'] = df.loc[:,'u_sdss'] - df.loc[:,'g_sdss']
    df.loc[:,'gr_sdss'] = df.loc[:,'g_sdss'] - df.loc[:,'r_sdss']
    df.loc[:,'ri_sdss'] = df.loc[:,'r_sdss'] - df.loc[:,'i_sdss']
    df.loc[:,'iz_sdss'] = df.loc[:,'i_sdss'] - df.loc[:,'z_sdss']
    df.loc[:,'gi_sdss'] = df.loc[:,'g_sdss'] - df.loc[:,'i_sdss']
    df.loc[:,'uNUV']    = df.loc[:,'u_sdss'] - df.loc[:,'NUV_galex']
    df.loc[:,'NUVg']    = df.loc[:,'NUV_galex'] - df.loc[:,'g_sdss']
    df.loc[:,'FUVNUV']  = df.loc[:,'FUV_galex'] - df.loc[:,'NUV_galex']

    # Create initial (and generous)  mask...
    mask = ( ( df['ug_sdss'] > 0.5 ) & ( df['ug_sdss'] < 2.0 ) &
             ( df['gr_sdss'] > 0.2 ) & ( df['gr_sdss'] < 0.8 ) &
             ( df['gi_sdss'] > 0.3 ) & ( df['gi_sdss'] < 1.0 ) &  
             ( df['uNUV'] > -10 ) & ( df['uNUV'] < 10.0 )      &
             ( df['NUVg'] > -10 ) & ( df['NUVg'] < 10.0 )        )

    # Make a backup copy of original df...
    df_orig = df.copy()

    # Make a backup copy of original mask...
    mask_orig = mask.copy()


    # Open fit results output file...
    try:
        fout = open(resultsFile, 'w')
    except IOError:
        sys.exit('Unable to write to file ' + resultsFile)

    # Write header to fit results output file...
    hdr = createFitResultsHeaderOutputLine(norder)
    fout.write(hdr+'\n')


    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Fit (u_sdss-NUV) vs. (NUV-g_sdss)...                            #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    if args.verbose > 0:
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print """# Fit (u_sdss-NUV) vs. (NUV-g_sdss)                               #"""
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print

    # Create names for use in QA plots...
    dmagName = '$u_{sdss} - NUV_{galex}$'
    colorName1 = '$(NUV_{galex} - g_{sdss})$'

    # Grab the original version of df from the backup copy...
    df = df_orig.copy()

    # Grab the original version of mask from the backup copy...
    mask = mask_orig.copy()

    # Iterate, with sigma-clipping...
    for i in range(niter):

        iiter = i + 1
        if args.verbose > 0:
            print """   iter%d...""" % ( iiter )

        # make a copy of original df, overwriting the old one...
        df = df[mask].copy()

        # Identify dmag and color1 series...
        dmag =  df.loc[:,'uNUV']
        color1 = df.loc[:,'NUVg']

        # Perform fit...
        p,perr,rms = transformFit1(color1, dmag, norder, args.verbose)
        df.loc[:,'res'] = residuals1(p, color1, dmag)

        # Identify outliers...
        stddev = df['res'].std()
        mask = (np.abs(df.res)< nsigma*stddev)

    outputLine = createFitResultsOutputLine(2, p, perr, rms, 'uNUV', 'NUVg')
    fout.write(outputLine+'\n')
    if args.verbose > 0:
        print outputLine
        print 

    # Create QA plots...
    res =  df.loc[:,'res']
    outputFileName = """%s.%s.%s.qa1.png""" % (qaFileBaseName, 'uNUV', 'NUVg')
    status = sdssGalexTransform1ColorQAPlots1(dmag, color1, res, norder, dmagName, colorName1, p, rms, outputFileName)


    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Fit (u_sdss-g_sdss) vs. (g_sdss-r_sdss)...                      #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    if args.verbose > 0:
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print """# Fit (u_sdss-g_sdss) vs. (g_sdss-r_sdss)                         #"""
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print

    # Create names for use in QA plots...
    dmagName = '$u_{sdss} - g_{sdss}$'
    colorName1 = '$(g-r)_{sdss}$'

    # Grab the original version of df from the backup copy...
    df = df_orig.copy()

    # Grab the original version of mask from the backup copy...
    mask = mask_orig.copy()

    # Iterate, with sigma-clipping...
    for i in range(niter):

        iiter = i + 1
        if args.verbose > 0:
            print """   iter%d...""" % ( iiter )

        # make a copy of original df, overwriting the old one...
        df = df[mask].copy()

        # Identify dmag and color1 series...
        dmag =  df.loc[:,'ug_sdss']
        color1 = df.loc[:,'gr_sdss']

        # Perform fit...
        p,perr,rms = transformFit1(color1, dmag, norder, args.verbose)
        df.loc[:,'res'] = residuals1(p, color1, dmag)

        # Identify outliers...
        stddev = df['res'].std()
        mask = (np.abs(df.res)< nsigma*stddev)

    outputLine = createFitResultsOutputLine(2, p, perr, rms, 'ug_sdss', 'gr_sdss')
    fout.write(outputLine+'\n')
    if args.verbose > 0:
        print outputLine
        print 

    # Create QA plots...
    res =  df.loc[:,'res']
    outputFileName = """%s.%s.%s.qa1.png""" % (qaFileBaseName, 'ug_sdss', 'gr_sdss')
    status = sdssGalexTransform1ColorQAPlots1(dmag, color1, res, norder, dmagName, colorName1, p, rms, outputFileName)


    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Fit (u_sdss-g_sdss) vs. (g_sdss-i_sdss)...                      #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    if args.verbose > 0:
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print """# Fit (u_sdss-g_sdss) vs. (g_sdss-i_sdss)                         #"""
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print

    # Create names for use in QA plots...
    dmagName = '$u_{sdss} - g_{sdss}$'
    colorName1 = '$(g-i)_{sdss}$'

    # Grab the original version of df from the backup copy...
    df = df_orig.copy()

    # Grab the original version of mask from the backup copy...
    mask = mask_orig.copy()

    # Iterate, with sigma-clipping...
    for i in range(niter):

        iiter = i + 1
        if args.verbose > 0:
            print """   iter%d...""" % ( iiter )

        # make a copy of original df, overwriting the old one...
        df = df[mask].copy()

        # Identify dmag and color1 series...
        dmag =  df.loc[:,'ug_sdss']
        color1 = df.loc[:,'gi_sdss']

        # Perform fit...
        p,perr,rms = transformFit1(color1, dmag, norder, args.verbose)
        df.loc[:,'res'] = residuals1(p, color1, dmag)

        # Identify outliers...
        stddev = df['res'].std()
        mask = (np.abs(df.res)< nsigma*stddev)

    outputLine = createFitResultsOutputLine(2, p, perr, rms, 'ug_sdss', 'gi_sdss')
    fout.write(outputLine+'\n')
    if args.verbose > 0:
        print outputLine
        print 

    # Create QA plots...
    res =  df.loc[:,'res']
    outputFileName = """%s.%s.%s.qa1.png""" % (qaFileBaseName, 'ug_sdss', 'gi_sdss')
    status = sdssGalexTransform1ColorQAPlots1(dmag, color1, res, norder, dmagName, colorName1, p, rms, outputFileName)


    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Fit (u_sdss-NUV) vs. (NUV-g_sdss) and (g_sdss-r_sdss)...        #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    if args.verbose > 0:
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print """# Fit (u_sdss-NUV) vs. (NUV-g_sdss) and (g_sdss-r_sdss)...        #"""
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print

    # Create names for use in QA plots...
    dmagName = '$u_{sdss} - NUV_{galex}$'
    colorName1 = '$(NUV_{galex} - g_{sdss})$'
    colorName2 = '$(g-r)_{sdss}$'

    # Grab the original version of df from the backup copy...
    df = df_orig.copy()

    # Grab the original version of mask from the backup copy...
    mask = mask_orig.copy()

    # Iterate, with sigma-clipping...
    for i in range(niter):

        iiter = i + 1
        if args.verbose > 0:
            print """   iter%d...""" % ( iiter )

        # make a copy of original df, overwriting the old one...
        df = df[mask].copy()

        # Identify dmag, color1, and color2 series...
        dmag =  df.loc[:,'uNUV']
        color1 = df.loc[:,'NUVg']
        color2 = df.loc[:,'gr_sdss']

        # Perform fit...
        p,perr,rms = transformFit2(color1, color2, dmag, norder, args.verbose)
        df.loc[:,'res'] = residuals2(p, color1, color2, dmag)

        # Identify outliers...
        stddev = df['res'].std()
        mask = (np.abs(df.res)< nsigma*stddev)

    outputLine = createFitResultsOutputLine(norder, p, perr, rms, 'uNUV', 'NUVg', 'gr_sdss')
    fout.write(outputLine+'\n')
    if args.verbose > 0:
        print outputLine
        print 

    # Create QA plots...
    res =  df.loc[:,'res']
    outputFileName = """%s.%s.%s.%s.qa1.png""" % (qaFileBaseName, 'uNUV', 'NUVg', 'gr_sdss')
    status = sdssGalexTransform2ColorQAPlots1(dmag, color1, color2, res, norder, dmagName, colorName1, colorName2, p, rms, outputFileName)

    outputFileName = """%s.%s.%s.%s.qa2.html""" % (qaFileBaseName, 'uNUV', 'NUVg', 'gr_sdss')
    status = sdssGalexTransform2ColorQAPlots2(df, 'uNUV', 'NUVg', 'gr_sdss', 'res', dmagName, colorName1, colorName2, norder, p, rms, outputFileName)

    outputFileName = """%s.%s.%s.%s.qa3_res.html""" % (qaFileBaseName, 'uNUV', 'NUVg', 'gr_sdss')
    status = sdssGalexTransform2ColorQAPlots3(df, 'uNUV', 'NUVg', 'gr_sdss', 'res', dmagName, colorName1, colorName2, norder, p, rms, outputFileName)





    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # Fit (u_sdss-NUV) vs. (NUV-g_sdss) and (g_sdss-i_sdss)...        #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    if args.verbose > 0:
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print """# Fit (u_sdss-NUV) vs. (NUV-g_sdss) and (g_sdss-i_sdss)...        #"""
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print

    # Create names for use in QA plots...
    dmagName = '$u_{sdss} - NUV_{galex}$'
    colorName1 = '$(NUV_{galex} - g_{sdss})$'
    colorName2 = '$(g-i)_{sdss}$'

    # Grab the original version of df from the backup copy...
    df = df_orig.copy()

    # Grab the original version of mask from the backup copy...
    mask = mask_orig.copy()

    # Iterate, with sigma-clipping...
    for i in range(niter):

        iiter = i + 1
        if args.verbose > 0:
            print """   iter%d...""" % ( iiter )

        # make a copy of original df, overwriting the old one...
        df = df[mask].copy()

        # Identify dmag, color1, and color2 series...
        dmag =  df.loc[:,'uNUV']
        color1 = df.loc[:,'NUVg']
        color2 = df.loc[:,'gi_sdss']

        # Perform fit...
        p,perr,rms = transformFit2(color1, color2, dmag, norder, args.verbose)
        df.loc[:,'res'] = residuals2(p, color1, color2, dmag)

        # Identify outliers...
        stddev = df['res'].std()
        mask = (np.abs(df.res)< nsigma*stddev)
    
    outputLine = createFitResultsOutputLine(norder, p, perr, rms, 'uNUV', 'NUVg', 'gi_sdss')
    fout.write(outputLine+'\n')
    if args.verbose > 0:
        print outputLine
        print 

    # Create QA plots...
    res =  df.loc[:,'res']
    outputFileName = """%s.%s.%s.%s.qa1.png""" % (qaFileBaseName, 'uNUV', 'NUVg', 'gi_sdss')
    status = sdssGalexTransform2ColorQAPlots1(dmag, color1, color2, res, norder, dmagName, colorName1, colorName2, p, rms, outputFileName)

    outputFileName = """%s.%s.%s.%s.qa2.html""" % (qaFileBaseName, 'uNUV', 'NUVg', 'gi_sdss')
    status = sdssGalexTransform2ColorQAPlots2(df, 'uNUV', 'NUVg', 'gi_sdss', 'res', dmagName, colorName1, colorName2, norder, p, rms, outputFileName)

    outputFileName = """%s.%s.%s.%s.qa3_res.html""" % (qaFileBaseName, 'uNUV', 'NUVg', 'gi_sdss')
    status = sdssGalexTransform2ColorQAPlots3(df, 'uNUV', 'NUVg', 'gi_sdss', 'res', dmagName, colorName1, colorName2, norder, p, rms, outputFileName)



    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # Fit (u_sdss-g_sdss) vs. (g_sdss-r_sdss) and (r_sdss-i_sdss)...  #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    if args.verbose > 0:
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print """# Fit (u_sdss-g_sdss) vs. (g_sdss-r_sdss) and (r_sdss-i_sdss)...  #"""
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print

    # Create names for use in QA plots...
    dmagName = '$u_{sdss} - g_{sdss}$'
    colorName1 = '$(g-r)_{sdss}$'
    colorName2 = '$(r-i)_{sdss}$'

    # Grab the original version of df from the backup copy...
    df = df_orig.copy()

    # Grab the original version of mask from the backup copy...
    mask = mask_orig.copy()

    # Iterate, with sigma-clipping...
    for i in range(niter):

        iiter = i + 1
        if args.verbose > 0:
            print """   iter%d...""" % ( iiter )

        # make a copy of original df, overwriting the old one...
        df = df[mask].copy()

        # Identify dmag, color1, and color2 series...
        dmag =  df.loc[:,'ug_sdss']
        color1 = df.loc[:,'gr_sdss']
        color2 = df.loc[:,'ri_sdss']

        # Perform fit...
        p,perr,rms = transformFit2(color1, color2, dmag, norder, args.verbose)
        df.loc[:,'res'] = residuals2(p, color1, color2, dmag)

        # Identify outliers...
        stddev = df['res'].std()
        mask = (np.abs(df.res)< nsigma*stddev)
    
    outputLine = createFitResultsOutputLine(norder, p, perr, rms, 'ug_sdss', 'gr_sdss', 'ri_sdss')
    fout.write(outputLine+'\n')
    if args.verbose > 0:
        print outputLine
        print 

    # Create QA plots...
    res =  df.loc[:,'res']
    outputFileName = """%s.%s.%s.%s.qa1.png""" % (qaFileBaseName, 'ug_sdss', 'gr_sdss', 'ri_sdss')
    status = sdssGalexTransform2ColorQAPlots1(dmag, color1, color2, res, norder, dmagName, colorName1, colorName2, p, rms, outputFileName)

    outputFileName = """%s.%s.%s.%s.qa2.html""" % (qaFileBaseName, 'ug_sdss', 'gr_sdss', 'ri_sdss')
    status = sdssGalexTransform2ColorQAPlots2(df, 'ug_sdss', 'gr_sdss', 'ri_sdss', 'res', dmagName, colorName1, colorName2, norder, p, rms, outputFileName)

    outputFileName = """%s.%s.%s.%s.qa3_res.html""" % (qaFileBaseName, 'ug_sdss', 'gr_sdss', 'ri_sdss')
    status = sdssGalexTransform2ColorQAPlots3(df, 'ug_sdss', 'gr_sdss', 'ri_sdss', 'res', dmagName, colorName1, colorName2, norder, p, rms, outputFileName)



    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # Fit (u_sdss-g_sdss) vs. (g_sdss-i_sdss) and (r_sdss-i_sdss)...  #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    if args.verbose > 0:
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print """# Fit (u_sdss-g_sdss) vs. (g_sdss-i_sdss) and (r_sdss-i_sdss)...  #"""
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print

    # Create names for use in QA plots...
    dmagName = '$u_{sdss} - g_{sdss}$'
    colorName1 = '$(g-i)_{sdss}$'
    colorName2 = '$(r-i)_{sdss}$'

    # Grab the original version of df from the backup copy...
    df = df_orig.copy()

    # Grab the original version of mask from the backup copy...
    mask = mask_orig.copy()

    # Iterate, with sigma-clipping...
    for i in range(niter):

        iiter = i + 1
        if args.verbose > 0:
            print """   iter%d...""" % ( iiter )

        # make a copy of original df, overwriting the old one...
        df = df[mask].copy()

        # Identify dmag, color1, and color2 series...
        dmag =  df.loc[:,'ug_sdss']
        color1 = df.loc[:,'gi_sdss']
        color2 = df.loc[:,'ri_sdss']

        # Perform fit...
        p,perr,rms = transformFit2(color1, color2, dmag, norder, args.verbose)
        df.loc[:,'res'] = residuals2(p, color1, color2, dmag)

        # Identify outliers...
        stddev = df['res'].std()
        mask = (np.abs(df.res)< nsigma*stddev)

    outputLine = createFitResultsOutputLine(norder, p, perr, rms, 'ug_sdss', 'gi_sdss', 'ri_sdss')
    fout.write(outputLine+'\n')
    if args.verbose > 0:
        print outputLine
        print 

    # Create QA plots...
    res =  df.loc[:,'res']
    outputFileName = """%s.%s.%s.%s.qa1.png""" % (qaFileBaseName, 'ug_sdss', 'gi_sdss', 'ri_sdss')
    status = sdssGalexTransform2ColorQAPlots1(dmag, color1, color2, res, norder, dmagName, colorName1, colorName2, p, rms, outputFileName)

    outputFileName = """%s.%s.%s.%s.qa2.html""" % (qaFileBaseName, 'ug_sdss', 'gi_sdss', 'ri_sdss')
    status = sdssGalexTransform2ColorQAPlots2(df, 'ug_sdss', 'gi_sdss', 'ri_sdss', 'res', dmagName, colorName1, colorName2, norder, p, rms, outputFileName)

    outputFileName = """%s.%s.%s.%s.qa3_res.html""" % (qaFileBaseName, 'ug_sdss', 'gi_sdss', 'ri_sdss')
    status = sdssGalexTransform2ColorQAPlots3(df, 'ug_sdss', 'gi_sdss', 'ri_sdss', 'res', dmagName, colorName1, colorName2, norder, p, rms, outputFileName)

    # Close outputFile...
    fout.close()

    return 0



##################################
#
# Define some functions for fitting dmag vs. color...
#
# These functions are based on a scripts found at 
# http://linuxgazette.net/115/andreasen.html (by Anders Andreasen)
# and at
# http://www.phy.uct.ac.za/courses/python/examples/fitresonance.py (University of Cape Town)


##################################

# Parametric function:  
#  p is the parameter vector; 
#  For fp1, we assume a polynomial function in one color...
def fp1(p,color1_array):
    #retValue = p[0] + p[1]*color1_array + p[2]*color1_array*color1_array
    norder = p.size-1
    retValue = p[0]
    for i in range(norder):
        retValue = retValue + p[i+1]*color1_array**(i+1)
    return retValue


##################################

# Error function:
def residuals1(p,color1_array,dmag_array):
    err = (dmag_array-fp1(p,color1_array))
    return err


##################################

# Fitting code:
def transformFit1(color1_array, dmag_array, norder=2, verbose=0):

    # Calculate the median of dmag for use as an initial guess
    # for the overall zeropoint offset..
    mdn = np.median( dmag_array, None )

    # Parameter names
    #pname = (['c_0', 'c_1', 'c_2'])
    pname = []
    for i in range(0,norder+1):
        pname.append("""c_%d""" % i)

    # Initial parameter values
    #p0 = [mdn, 0.0, 0.0]
    p0 = (1+norder)*[0.0]
    p0[0] = mdn

    if verbose > 0:
        print 
        print 'Initial parameter values:  ', p0

    # Perform fit

    p,cov,infodict,mesg,ier = leastsq(residuals1, p0, 
                                      args=(color1_array, dmag_array), 
                                      maxfev=10000, full_output=1)

    if ( ier>=1 and ier <=4):
        if verbose > 0:  print "Converged"
    else:
        # Add an exception error or a non-zero return value?
        print "Not converged"
        print mesg


    # Calculate some descriptors of the fit 
    # (similar to the output from gnuplot 2d fits)

    chisq=sum(infodict['fvec']*infodict['fvec'])
    dof=len(dmag_array)-len(p)
    rms=math.sqrt(chisq/dof)
    
    if verbose > 0:
        print "Converged with chi squared ",chisq
        print "degrees of freedom, dof ", dof
        print "RMS of residuals (i.e. sqrt(chisq/dof)) ", rms
        print "Reduced chisq (i.e. variance of residuals) ", chisq/dof
        print


    # uncertainties are calculated as per gnuplot, "fixing" the result
    # for non unit values of the reduced chisq.
    # values at min match gnuplot
    perr = []
    if verbose > 0:  
        print "Fitted parameters at minimum, with 68% C.I.:"
    for i,pmin in enumerate(p):
        if verbose > 0:  
            print "%-10s %13g +/- %13g   (%5f percent)" % (pname[i],pmin,math.sqrt(cov[i,i])*math.sqrt(chisq/dof),
                                                           100.*math.sqrt(cov[i,i])*math.sqrt(chisq/dof)/abs(pmin))
        perr.append(math.sqrt(cov[i,i])*math.sqrt(chisq/dof))
    if verbose > 0: print

    if verbose > 0:
        print "Correlation matrix:"
        # correlation matrix close to gnuplot
        print "               ",
        for i in range(len(pname)): print "%-10s" % (pname[i],),
        print
        for i in range(len(p)):
            print "%-10s" % pname[i],
            for j in range(i+1):
                print "%10f" % (cov[i,j]/math.sqrt(cov[i,i]*cov[j,j]),),
            #endfor
            print
        #endfor
        print
        print
        print
    
    return p, perr, rms


##################################
#
# Define some functions for fitting dmag vs. color1 and color2...
#
# These functions are based on a scripts found at 
# http://linuxgazette.net/115/andreasen.html (by Anders Andreasen)
# and at
# http://www.phy.uct.ac.za/courses/python/examples/fitresonance.py (University of Cape Town)


##################################

# Parametric function:  
#  p is the parameter vector; 
#  For fp2, we assume a polynomial in each of the 2 colors 
#   but with no cross terms... 
def fp2(p,color1_array,color2_array):
    #retValue = p[0] + \
    #    p[1]*color1_array + p[2]*color1_array*color1_array + \
    #    p[3]*color2_array + p[4]*color2_array*color2_array
    norder = (p.size-1)/2
    retValue = p[0]
    for i in range(norder):
        retValue = retValue + p[i+1]*color1_array**(i+1)
        retValue = retValue + p[i+norder+1]*color2_array**(i+1)
    return retValue


##################################

# Error function:
def residuals2(p,color1_array,color2_array,dmag_array):
    err = (dmag_array-fp2(p,color1_array,color2_array))
    return err


##################################

# Fitting code:
def transformFit2(color1_array, color2_array, dmag_array, norder=2, verbose=0):
    
    # Calculate the median of dmag for use as an initial guess
    # for the overall zeropoint offset..
    mdn = np.median( dmag_array, None )

    # Parameter names
    #pname = (['c_0', 'c_1', 'c_2', 'c_3', 'c_4'])
    pname = []
    for i in range(0,2*norder+1):
        pname.append("""c_%d""" % i)

    # Initial parameter values
    #p0 = [mdn, 0.0, 0.0, 0.0, 0.0]
    p0 = (1+2*norder)*[0.0]
    p0[0] = mdn

    if verbose > 0:
        print 
        print 'Initial parameter values:  ', p0


    # Perform fit
    p,cov,infodict,mesg,ier = leastsq(residuals2, p0, 
                                      args=(color1_array, color2_array, dmag_array), 
                                      maxfev=10000, full_output=1)

    if ( ier>=1 and ier <=4):
        if verbose > 0: print "Converged"
    else:
        # Add an exception error or a non-zero return value?
        print "Not converged"
        print mesg


    # Calculate some descriptors of the fit 
    # (similar to the output from gnuplot 2d fits)

    chisq=sum(infodict['fvec']*infodict['fvec'])
    dof=len(dmag_array)-len(p)
    rms=math.sqrt(chisq/dof)
    
    if verbose > 0:
        print "Converged with chi squared ",chisq
        print "degrees of freedom, dof ", dof
        print "RMS of residuals (i.e. sqrt(chisq/dof)) ", rms
        print "Reduced chisq (i.e. variance of residuals) ", chisq/dof
        print


    # uncertainties are calculated as per gnuplot, "fixing" the result
    # for non unit values of the reduced chisq.
    # values at min match gnuplot
    perr = []
    if verbose > 0:
        print "Fitted parameters at minimum, with 68% C.I.:"
    for i,pmin in enumerate(p):
        if verbose > 0:
            print "%-10s %13g +/- %13g   (%5f percent)" % (pname[i],pmin,math.sqrt(cov[i,i])*math.sqrt(chisq/dof),
                                                           100.*math.sqrt(cov[i,i])*math.sqrt(chisq/dof)/abs(pmin))
        perr.append(math.sqrt(cov[i,i])*math.sqrt(chisq/dof))
    if verbose > 0: print

    if verbose > 0:
        print "Correlation matrix:"
        # correlation matrix close to gnuplot
        print "               ",
        for i in range(len(pname)): print "%-10s" % (pname[i],),
        print
        for i in range(len(p)):
            print "%-10s" % pname[i],
            for j in range(i+1):
                print "%10f" % (cov[i,j]/math.sqrt(cov[i,i]*cov[j,j]),),
            #endfor
            print
        #endfor
        print
        print
        print
    
    return p, perr, rms


##################################

def createFitResultsOutputLine(norder, p, perr, rms, dmag_name, color1_name, color2_name=''):

    outputList = (2*(2*norder+1)+4)*[-9999.]
    outputList[0] = dmag_name
    outputList[1] = color1_name
    outputList[2] = color2_name
    for j in range(p.size):
        outputList[2*j+3] = p[j]
        outputList[2*j+4] = perr[j]
    outputList[2*(2*norder+1)+3] = rms
    outputLine = ','.join(map(str, outputList))
    return outputLine


##################################

def createFitResultsHeaderOutputLine(norder):

    outputList = (2*(2*norder+1)+4)*['c_']
    outputList[0] = 'dmag_name'
    outputList[1] = 'color1_name'
    outputList[2] = 'color2_name'
    for j in range(2*norder+1):
        outputList[2*j+3] = ("""c_%d""" % j)
        outputList[2*j+4] = ("""cerr_%d""" % j)
    outputList[2*(2*norder+1)+3] = 'rms'
    outputLine = ','.join(map(str, outputList))
    return outputLine


##################################

def sdssGalexTransform1ColorQAPlots1(dmag, color1, res, norder, dmagName, colorName1, p, rms, outputFileName):

    # Prepare QA plots...
    fig = plt.figure(figsize=(10,5))
    fig.subplots_adjust(hspace=0.3)
    #fig.suptitle("This is a supertitle!")

    # We will exclude the lowest and highets 1% of color1, color2, 
    #  dmag, and residuals when plotting the QA figures...
    color1_desc = color1.describe(percentiles=[0.01, 0.99])
    dmag_desc = dmag.describe(percentiles=[0.01, 0.99])
    #res_desc = df.res.describe(percentiles=[0.01, 0.99])
    res_desc = res.describe(percentiles=[0.01, 0.99])
    color1_min = color1_desc['1%']
    color1_max = color1_desc['99%']
    dmag_min = dmag_desc['1%']
    dmag_max = dmag_desc['99%']
    res_min = res_desc['1%']
    res_max = res_desc['99%']

    # Plot 1:  Descriptive text...
    plt.subplot(231)
    if norder == 1:
        plot1Text = """%s = \n %.3f + \n %.3f*%s \n\n [rms: %.3f]""" % \
            (dmagName, p[0], p[1], colorName1, rms)
    elif norder == 2:
        plot1Text = """%s = \n %.3f + \n %.3f*%s + \n %.3f*%s^2 \n\n [rms: %.3f]""" % \
            (dmagName, p[0], p[1], colorName1, p[2], colorName1, rms)
    else:
        plot1Text = ''
    plt.text(0.1,0.25,plot1Text)
    plt.axis('off')

    # Plot 2:  2D hexbin histogram of dmag vs. color1...
    plt.subplot(232) 
    hb=plt.hexbin(color1, dmag, gridsize=100, cmap='inferno')
    plt.axis([color1_min, color1_max, dmag_min, dmag_max])
    plt.xlabel(colorName1)
    plt.ylabel(dmagName)
    cb = fig.colorbar(hb)
    cb.set_label('Number')
    plt.grid(color='white')
    plt.grid(True)

    # Plot 3:  N/A

    # Plot 4:  1d histogram of residuals...
    plt.subplot(234) 
    #plt.hist(df.loc[:,'res'],bins=100)
    plt.hist(res,bins=100)
    plt.xlabel('residuals [mag]')
    plt.ylabel('Number')
    plt.grid(True)
    plt.grid(color='black')

    # Plot 5:  2d hexbin histogram of residuals vs. color1...
    plt.subplot(235) 
    #hb = plt.hexbin(color1, df.loc[:,'res'], gridsize=100, cmap='inferno')
    hb = plt.hexbin(color1, res, gridsize=100, cmap='inferno')
    plt.axis([color1_min, color1_max, res_min, res_max])
    plt.xlabel(colorName1)
    plt.ylabel('residuals [mag]')
    cb = plt.colorbar(hb)
    cb.set_label('Number')
    plt.grid(True)
    plt.grid(color='white')

    # Plot 6:  N/A

    # Plot...
    plt.tight_layout()
    #plt.show()
    plt.savefig(outputFileName)

    return 0


##################################

def sdssGalexTransform2ColorQAPlots1(dmag, color1, color2, res, norder, dmagName, colorName1, colorName2, p, rms, outputFileName):

    # Prepare QA plots...
    fig = plt.figure(figsize=(10,5))
    fig.subplots_adjust(hspace=0.3)
    #fig.suptitle("This is a supertitle!")

    # We will exclude the lowest and highets 1% of color1, color2, 
    #  dmag, and residuals when plotting the QA figures...
    color1_desc = color1.describe(percentiles=[0.01, 0.99])
    color2_desc = color2.describe(percentiles=[0.01, 0.99])
    dmag_desc = dmag.describe(percentiles=[0.01, 0.99])
    #res_desc = df.res.describe(percentiles=[0.01, 0.99])
    res_desc = res.describe(percentiles=[0.01, 0.99])
    color1_min = color1_desc['1%']
    color1_max = color1_desc['99%']
    color2_min = color2_desc['1%']
    color2_max = color2_desc['99%']
    dmag_min = dmag_desc['1%']
    dmag_max = dmag_desc['99%']
    res_min = res_desc['1%']
    res_max = res_desc['99%']

    # Plot 1:  Descriptive text...
    plt.subplot(231)
    if norder == 1:
        plot1Text = """%s = \n %.3f + \n %.3f*%s + \n %.3f*%s \n\n [rms: %.3f]""" % \
            (dmagName, p[0], p[1], colorName1, p[2], colorName2, rms)
    elif norder == 2:
        plot1Text = """%s = \n %.3f + \n %.3f*%s + \n %.3f*%s^2 + \n %.3f*%s + \n %.3f*%s^2 \n\n [rms: %.3f]""" % \
            (dmagName, p[0], p[1], colorName1, p[2], colorName1, p[3], colorName2, p[4], colorName2, rms)
    else:
        plot1Text = ''
    plt.text(0.1,0.25,plot1Text)
    plt.axis('off')

    # Plot 2:  2D hexbin histogram of dmag vs. color1...
    plt.subplot(232) 
    hb=plt.hexbin(color1, dmag, gridsize=100, cmap='inferno')
    plt.axis([color1_min, color1_max, dmag_min, dmag_max])
    plt.xlabel(colorName1)
    plt.ylabel(dmagName)
    cb = fig.colorbar(hb)
    cb.set_label('Number')
    plt.grid(color='white')
    plt.grid(True)

    # Plot 3:  2D hexbin histogram of dmag vs. color2...
    plt.subplot(233) 
    hb=plt.hexbin(color2, dmag, gridsize=100, cmap='inferno')
    plt.axis([color2_min, color2_max, dmag_min, dmag_max])
    plt.xlabel(colorName2)
    plt.ylabel(dmagName)
    cb = plt.colorbar(hb)
    cb.set_label('Number')
    plt.grid(color='white')
    plt.grid(True)

    # Plot 4:  1d histogram of residuals...
    plt.subplot(234) 
    #plt.hist(df.loc[:,'res'],bins=100)
    plt.hist(res,bins=100)
    plt.xlabel('residuals [mag]')
    plt.ylabel('Number')
    plt.grid(True)
    plt.grid(color='black')

    # Plot 5:  2d hexbin histogram of residuals vs. color1...
    plt.subplot(235) 
    #hb = plt.hexbin(color1, df.loc[:,'res'], gridsize=100, cmap='inferno')
    hb = plt.hexbin(color1, res, gridsize=100, cmap='inferno')
    plt.axis([color1_min, color1_max, res_min, res_max])
    plt.xlabel(colorName1)
    plt.ylabel('residuals [mag]')
    cb = plt.colorbar(hb)
    cb.set_label('Number')
    plt.grid(True)
    plt.grid(color='white')

    # Plot 6:  2d hexbin histogram of residuals vs. color2...
    plt.subplot(236) 
    #hb = plt.hexbin(color2, df.loc[:,'res'], gridsize=100, cmap='inferno')
    hb = plt.hexbin(color2, res, gridsize=100, cmap='inferno')
    plt.axis([color2_min, color2_max, res_min, res_max])
    plt.xlabel(colorName2)
    plt.ylabel('residuals [mag]')
    cb = plt.colorbar(hb)
    cb.set_label('Number')
    plt.grid(True)
    plt.grid(color='white')

    # Plot...
    plt.tight_layout()
    #plt.show()
    plt.savefig(outputFileName)

    return 0


##################################

# Based on plotly code at http://inversionlabs.com/2016/03/21/best-fit-surface.html
def sdssGalexTransform2ColorQAPlots2(df, colName_dmag, colName_color1, colName_color2, colName_res, dmagName, colorName1, colorName2, norder, p, rms, outputFileName):

    # An initial sanity check...
    if norder > 2:
        print 'sdssGalexTransform2ColorQAPlots2 can not deal with norder > 2...  skipping...'
        return 1

    # Create data from color1, color2, and dmag...
    #  If the sample size is larger than 1000, 
    #  take a random sample of 1000 elements...
    n_elements = df[colName_res].size
    if n_elements <= 1000:
        x = df.loc[:,colName_color1].values
        y = df.loc[:,colName_color2].values
        z = df.loc[:,colName_dmag].values
    else:
        df1000 = df.sample(n=1000,axis=0)
        n_elements = df1000[colName_res].size
        x = df1000.loc[:,colName_color1].values
        y = df1000.loc[:,colName_color2].values
        z = df1000.loc[:,colName_dmag].values    
    data = np.c_[x,y,z]

    
    # Regular grid covering the domain of the data...
    mn = np.min(data, axis=0)
    mx = np.max(data, axis=0)
    X,Y = np.meshgrid(np.linspace(mn[0], mx[0], 20), np.linspace(mn[1], mx[1], 20))
    XX = X.flatten()
    YY = Y.flatten()
    
    # Evaluate it on grid...
    if norder == 1:
        Z = p[0] + p[1]*X + p[2]*Y
    elif norder == 2:
        Z = p[0] + p[1]*X + p[2]*X**2 + p[3]*Y + p[4]*Y**2        
    
    # Create data_fit from color1, color2, and the fit parameters p[*]...
    x_fit = x
    y_fit = y
    if norder == 1:
        z_fit = p[0] + p[1]*x_fit + p[2]*y_fit
    elif norder == 2:
        z_fit = p[0] + p[1]*x_fit + p[2]*x_fit*x_fit + p[3]*y_fit + p[4]*y_fit*y_fit
    data_fit = np.c_[x_fit,y_fit,z_fit]
    #delta_z = z - z_fit

    # trace1 is the scatter plot of the original points...
    trace1 = go.Scatter3d(
        x=data[:,0],
        y=data[:,1],
        z=data[:,2], 
        #z=delta_z, 
        mode='markers',
        marker=dict(size=1, color='red', line=dict(color='black', width=0.5), opacity=0.5)
        )

    # trace2 is the scatter plot of the fit values at the x,y positions of the original points...
    trace2 = go.Scatter3d(
        x=data_fit[:,0],
        y=data_fit[:,1],
        z=data_fit[:,2],
        mode='markers',
        marker=dict(size=2, color='yellow', line=dict(color='black', width=0.5), opacity=0.8)
        )

    # trace3 is the 2D surface of the fit equation...
    trace3 = go.Surface(z=Z, x=X, y=Y, colorscale='RdBu', opacity=0.333)

    # Package the trace dictionaries into a data object
    data_go = go.Data([trace1, trace2, trace3])

    # Dictionary of style options for all axes
    axis = dict(
        showbackground=True, # show axis background
        backgroundcolor="rgb(204, 204, 204)", # set background color to grey
        gridcolor="rgb(255, 255, 255)",       # set grid line color
        zerolinecolor="rgb(255, 255, 255)",   # set zero grid line color
        )

    # Create a title...
    if norder == 1:
        titleText = """%s = %.3f + %.3f*%s + %.3f*%s [rms: %.3f]""" % \
            (dmagName, p[0], p[1], colorName1, p[2], colorName2, rms)
    elif norder == 2:
        titleText = """%s = %.3f + %.3f*%s + %.3f*%s^2 + %.3f*%s + %.3f*%s^2\n[npts=%d, rms: %.3f]""" % \
            (dmagName, p[0], p[1], colorName1, p[2], colorName1, p[3], colorName2, p[4], colorName2, n_elements, rms)
    else:
        titleText = ''
    titleText = titleText.replace('$','')

    # Make a layout object
    layout = go.Layout(
        title=titleText, # set plot title
        scene=go.Scene(  # axes are part of a 'scene' in 3d plots
            xaxis=go.XAxis(axis), # set x-axis style
            yaxis=go.YAxis(axis), # set y-axis style
            zaxis=go.ZAxis(axis)),  # set z-axis style
        )


    # Make a figure object
    fig = go.Figure(data=data_go, layout=layout)

    # Create interactive plot and save as javascript to html file...
    #plotly.offline.iplot(fig, filename=outputFileName)
    plotly.offline.plot(fig, filename=outputFileName, auto_open=False)

    return 0


##################################

# Based on plotly code at http://inversionlabs.com/2016/03/21/best-fit-surface.html
def sdssGalexTransform2ColorQAPlots3(df, colName_dmag, colName_color1, colName_color2, colName_res, dmagName, colorName1, colorName2, norder, p, rms, outputFileName):

    # An initial sanity check...
    if norder > 2:
        print 'sdssGalexTransform2ColorQAPlots3 can not deal with norder > 2...  skipping...'
        return 1

    # Create data from color1, color2, and res...
    #  If the sample size is larger than 1000, 
    #  take a random sample of 1000 elements...
    n_elements = df[colName_res].size
    if n_elements <= 1000:
        x = df.loc[:,colName_color1].values
        y = df.loc[:,colName_color2].values
        z = df.loc[:,colName_res].values
    else:
        df1000 = df.sample(n=1000,axis=0)
        n_elements = df1000[colName_res].size
        x = df1000.loc[:,colName_color1].values
        y = df1000.loc[:,colName_color2].values
        z = df1000.loc[:,colName_res].values    
    data = np.c_[x,y,z]

    
    # Regular grid covering the domain of the data...
    mn = np.min(data, axis=0)
    mx = np.max(data, axis=0)
    X,Y = np.meshgrid(np.linspace(mn[0], mx[0], 20), np.linspace(mn[1], mx[1], 20))
    XX = X.flatten()
    YY = Y.flatten()
    
    # Evaluate it on grid...
    Z = 0.00 + 0.00*X + 0.00*Y
    
    # Create data_fit from color1, color2, and the fit parameters p[*]...
    x_fit = x
    y_fit = y
    z_fit = 0.00 + 0.00*x_fit + 0.00*y_fit
    data_fit = np.c_[x_fit,y_fit,z_fit]

    # trace1 is the scatter plot of the original points...
    trace1 = go.Scatter3d(
        x=data[:,0],
        y=data[:,1],
        z=data[:,2], 
        mode='markers',
        marker=dict(size=1, color='red', line=dict(color='black', width=0.5), opacity=0.5)
        )

    # trace2 is the scatter plot of the fit values at the x,y positions of the original points...
    trace2 = go.Scatter3d(
        x=data_fit[:,0],
        y=data_fit[:,1],
        z=data_fit[:,2],
        mode='markers',
        marker=dict(size=2, color='yellow', line=dict(color='black', width=0.5), opacity=0.8)
        )

    # trace3 is the 2D surface of the fit equation...
    #trace3 = go.Surface(z=Z, x=X, y=Y, colorscale='RdBu', opacity=0.667)
    trace3 = go.Surface(z=Z, x=X, y=Y, colorscale='Greys', opacity=0.667)

    # Package the trace dictionaries into a data object
    data_go = go.Data([trace1, trace2, trace3])

    # Dictionary of style options for all axes
    axis = dict(
        showbackground=True, # show axis background
        backgroundcolor="rgb(204, 204, 204)", # set background color to grey
        gridcolor="rgb(255, 255, 255)",       # set grid line color
        zerolinecolor="rgb(255, 255, 255)",   # set zero grid line color
        )

    # Create a title...
    if norder == 1:
        titleText = """%s = %.3f + %.3f*%s + %.3f*%s [rms: %.3f]""" % \
            (dmagName, p[0], p[1], colorName1, p[2], colorName2, rms)
    elif norder == 2:
        titleText = """%s = %.3f + %.3f*%s + %.3f*%s^2 + %.3f*%s + %.3f*%s^2\n[npts=%d, rms: %.3f]""" % \
            (dmagName, p[0], p[1], colorName1, p[2], colorName1, p[3], colorName2, p[4], colorName2, n_elements, rms)
    else:
        titleText = ''
    titleText = titleText.replace('$','')

    # Make a layout object
    layout = go.Layout(
        title=titleText, # set plot title
        scene=go.Scene(  # axes are part of a 'scene' in 3d plots
            xaxis=go.XAxis(axis), # set x-axis style
            yaxis=go.YAxis(axis), # set y-axis style
            zaxis=go.ZAxis(axis)),  # set z-axis style
        )


    # Make a figure object
    fig = go.Figure(data=data_go, layout=layout)

    # Create interactive plot and save as javascript to html file...
    #plotly.offline.iplot(fig, filename=outputFileName)
    plotly.offline.plot(fig, filename=outputFileName, auto_open=False)

    return 0


##################################

if __name__ == "__main__":
    main()

##################################
