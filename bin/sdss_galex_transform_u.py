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
             ( df['gi_sdss'] > 0.3 ) & ( df['gr_sdss'] < 1.0 ) &  
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


    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Fit (u_sdss-g_sdss) vs. (g_sdss-r_sdss)...                      #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    if args.verbose > 0:
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print """# Fit (u_sdss-g_sdss) vs. (g_sdss-r_sdss)                         #"""
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print

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


    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Fit (u_sdss-g_sdss) vs. (g_sdss-i_sdss)...                      #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    if args.verbose > 0:
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print """# Fit (u_sdss-g_sdss) vs. (g_sdss-i_sdss)                         #"""
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print

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


    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # Fit (u_sdss-NUV) vs. (NUV-g_sdss) and (g_sdss-r_sdss)...        #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    if args.verbose > 0:
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print """# Fit (u_sdss-NUV) vs. (NUV-g_sdss) and (g_sdss-r_sdss)...        #"""
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print

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


    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # Fit (u_sdss-NUV) vs. (NUV-g_sdss) and (g_sdss-i_sdss)...        #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    if args.verbose > 0:
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print """# Fit (u_sdss-NUV) vs. (NUV-g_sdss) and (g_sdss-i_sdss)...        #"""
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print

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


    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # Fit (u_sdss-g_sdss) vs. (g_sdss-r_sdss) and (r_sdss-i_sdss)...  #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    if args.verbose > 0:
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print """# Fit (u_sdss-g_sdss) vs. (g_sdss-r_sdss) and (r_sdss-i_sdss)...  #"""
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print

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


    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # Fit (u_sdss-g_sdss) vs. (g_sdss-i_sdss) and (r_sdss-i_sdss)...  #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    if args.verbose > 0:
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print """# Fit (u_sdss-g_sdss) vs. (g_sdss-i_sdss) and (r_sdss-i_sdss)...  #"""
        print """# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"""
        print

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


    fout.close()


    x=color1
    #y=dmag
    y = df.loc[:,'res']
    #xmin = x.min()
    #xmax = x.max()
    #ymin = y.min()
    #ymax = y.max()
    xdesc = x.describe(percentiles=[0.01, 0.05, 0.50, 0.95, 0.99])
    ydesc = y.describe(percentiles=[0.01, 0.05, 0.50, 0.95, 0.99])
    xmin = xdesc['1%']
    xmax = xdesc['99%']
    ymin = ydesc['1%']
    ymax = ydesc['99%']
    fig, axs = plt.subplots(ncols=1)
    ax=axs
    #hb = ax.hexbin(x, y, gridsize=100, bins='log', cmap='inferno')
    hb = ax.hexbin(x, y, gridsize=100, cmap='inferno')
    ax.axis([xmin, xmax, ymin, ymax])
    ax.set_title('Test')
    ax.set_xlabel('color')
    #ax.set_ylabel('dmag')
    ax.set_ylabel('residuals')
    cb = fig.colorbar(hb, ax=ax)
    #cb.set_label('log(N)')
    cb.set_label('Number')
    plt.grid(True)
    ax.grid(color='white')
    plt.show()

    x=color1
    y=color2
    z = df.loc[:,'res']
    xdesc = x.describe(percentiles=[0.01, 0.05, 0.50, 0.95, 0.99])
    ydesc = y.describe(percentiles=[0.01, 0.05, 0.50, 0.95, 0.99])
    zdesc = y.describe(percentiles=[0.01, 0.05, 0.50, 0.95, 0.99])
    xmin = xdesc['1%']
    xmax = xdesc['99%']
    ymin = ydesc['1%']
    ymax = ydesc['99%']
    zmin = zdesc['1%']
    zmax = zdesc['99%']
    fig, axs = plt.subplots(ncols=1)
    ax=axs
    hb = ax.hexbin(x, y, C=z, gridsize=100, cmap='rainbow', reduce_C_function=np.median)
    ax.axis([xmin, xmax, ymin, ymax])
    ax.set_title('Test')
    ax.set_xlabel('color1')
    ax.set_ylabel('color2')
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('residuals')
    plt.grid(True)
    ax.grid(color='white')
    plt.show()



    return 0



##################################
#
# Define some functions for fitting dmag vs. color...
#
# These functions are based on a scripts found at 
# http://linuxgazette.net/115/andreasen.html (by Anders Andreasen)
# and at
# http://www.phy.uct.ac.za/courses/python/examples/fitresonance.py (University of Cape Town)

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

# Error function:
def residuals1(p,color1_array,dmag_array):
    err = (dmag_array-fp1(p,color1_array))
    return err

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

# Error function:
def residuals2(p,color1_array,color2_array,dmag_array):
    err = (dmag_array-fp2(p,color1_array,color2_array))
    return err

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

if __name__ == "__main__":
    main()

##################################
