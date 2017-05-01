#!/usr/bin/env python

# Initial setup...
import numpy as np
import pandas as pd
from astropy.io import fits
import fitsio
from scipy import interpolate
import glob
import math
import os
import shutil
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

#--------------------------------------------------------------------------
# Main code.

def main():

    import argparse

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--planFile', help='name of the input plan file', default='Zeropoints_for_QA_plan.par')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    status = zeropoints_for_qa(args)

    return status

#--------------------------------------------------------------------------

def zeropoints_for_qa(args):

    import paramFile

    # Read contents of planFile into python dictionary...
    planDict = paramFile.readParamFile(args.planFile, args.verbose)

    # Check that all the required keywords were found...
    requiredKeywordList = ['matchFileListFile', 
                           'qaDir', 
                           'blancoOpticsFile', 
                           'procEpochsFile',
                           'combEpochsFile']

    flag = 0
    for requiredKeyword in requiredKeywordList:
        if requiredKeyword not in planDict:
            print """Required keyword '%s' not in planFile %s""" % (requiredKeyword, args.planFile)
            flag = flag + 1
    if flag > 0:
        print 'Returning with error code 1 now...'
        return 1

    # Read in the matchFileListFile...
    matchFileListFile = planDict['matchFileListFile']
    if os.path.isfile(matchFileListFile)==False:
        print """matchFileListFile %s does not exist...""" % (matchFileListFile)
        print 'Returning with error code 1 now...'
        return 1
    inputFileArray = np.genfromtxt(matchFileListFile,dtype='str')

    # Create qaDir if it does not already exist...
    # See 
    # http://stackoverflow.com/questions/273192/how-to-check-if-a-directory-exists-and-create-it-if-necessary/14364249#14364249 
    # for avoiding the race condition beteween checking  a directory exists 
    # and creating the directory... 
    qaDir = planDict['qaDir']
    try: 
        os.makedirs(qaDir)
    except OSError:
        if not os.path.isdir(qaDir):
            raise

    # Read in the Blanco optics cleaning history file...
    blancoOpticsFile = planDict['blancoOpticsFile']
    if blancoOpticsFile.lower() == 'default':
        # Grab path and name of blanco_optics_cleaning_events.csv file 
        #  in the DECam_PGCM data directory...
        #  Is there a better way to do this?
        #  See also:
        #  http://stackoverflow.com/questions/779495/python-access-data-in-package-subdirectory
        # Absolute path for the directory containing this module:
        moduleAbsPathName = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        blancoOpticsFile = os.path.join(moduleAbsPathName, "data", "blanco_optics_cleaning_events.csv")
    if os.path.isfile(blancoOpticsFile)==False:
        print """blancoOpticsFile %s does not exist...""" % (blancoOpticsFile)
        print 'Returning with error code 1 now...'
        return 1
    blancoOpticsDF = pd.read_csv(blancoOpticsFile, comment='#')
    if args.verbose > 1:
        print blancoOpticsDF


    # Read in the data processing epochs history file...
    procEpochsFile = planDict['procEpochsFile']
    if procEpochsFile.lower() == 'default':
        # Grab path and name of processing_epochs.csv file 
        #  in the DECam_PGCM data directory...
        #  Is there a better way to do this?
        #  See also:
        #  http://stackoverflow.com/questions/779495/python-access-data-in-package-subdirectory
        # Absolute path for the directory containing this module:
        moduleAbsPathName = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        procEpochsFile = os.path.join(moduleAbsPathName, "data", "processing_epochs.csv")
    if os.path.isfile(procEpochsFile)==False:
        print """procEpochsFile %s does not exist...""" % (procEpochsFile)
        print 'Returning with error code 1 now...'
        return 1
    procEpochsDF = pd.read_csv(procEpochsFile, comment='#')
    if args.verbose > 1:
        print procEpochsDF

    # Read in the combined epochs history file...
    combEpochsFile = planDict['combEpochsFile']
    if combEpochsFile.lower() == 'default':
        # Grab path and name of combined_epochs.csv file 
        #  in the DECam_PGCM data directory...
        #  Is there a better way to do this?
        #  See also:
        #  http://stackoverflow.com/questions/779495/python-access-data-in-package-subdirectory
        # Absolute path for the directory containing this module:
        moduleAbsPathName = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        combEpochsFile = os.path.join(moduleAbsPathName, "data", "combined_epochs.csv")
    if os.path.isfile(combEpochsFile)==False:
        print """combEpochsFile %s does not exist...""" % (combEpochsFile)
        print 'Returning with error code 1 now...'
        return 1
    combEpochsDF = pd.read_csv(combEpochsFile, comment='#')
    if args.verbose > 1:
        print combEpochsDF


    # Loop through the input files in inpuFileArray
    for inputFile in inputFileArray:

        print inputFile

        baseName = os.path.basename(inputFile)
        baseNameNoExt = os.path.splitext(baseName)[0]

        df = pd.read_csv(inputFile, usecols=['EXPNUM_2',
                                             'BAND_2',
                                             'MJD_OBS_2',
                                             'AIRMASS_2',
                                             'FLUX_PSF_2', 
                                             'EXPTIME_2', 
                                             'MAG_PSF_MEDIAN_1'])

        df.rename(columns={'EXPNUM_2':'expnum', 
                           'BAND_2':'band', 
                           'MJD_OBS_2':'mjd', 
                           'AIRMASS_2':'airmass', 
                           'FLUX_PSF_2':'flux_psf',  
                           'EXPTIME_2':'exptime', 
                           'MAG_PSF_MEDIAN_1':'mag_std'
                           }, inplace=True)
    
        # Check to ensure that there is only one type of filter band in this data file...
        bandsList = df.band.unique() 
        if bandsList.size > 1:
            print 
            print """inputFile %s appears to contain multiple bands (%s)...  Skipping""" % \
                (inputFile, bandsList)
            print 
            continue
        else:
            band = bandsList[0]
            print """inputFile %s contains data for band (%s)...""" % \
                (inputFile, band)

        # Add dmag [= mag_std - mag_inst] column
        df['dmag'] = df['mag_std'] + 2.5*np.log10(df['flux_psf']/df['exptime'])

        # Let's plot a 2D histogram of log(Nobs), binned by dmag and airmass...
        x=df.mjd
        y=df.dmag
        fig, axs = plt.subplots(ncols=1)
        ax=axs
        hb = ax.hexbin(x, y, gridsize=100, bins='log', cmap='inferno')
        ax.set_title(baseName)
        ax.set_xlabel("MJD")
        ax.set_ylabel("dmag")
        cb = fig.colorbar(hb, ax=ax)
        cb.set_label('log10(N)')
        plt.grid(True)
        ax.grid(color='white')
        #plt.show()
        outputFile = """%s/%s.dmag_vs_mjd.png"""% (qaDir, baseNameNoExt)
        fig.savefig(outputFile)
        plt.close()
        del x
        del y

        # Save df in df_orig for future reference...
        df_orig = df.copy()

        ## Split df into two epochs (split at Dec 6 mirror washing)...
        #df1 = df_orig[df.mjd<56268.0].copy()
        #df2 = df_orig[df.mjd>=56268.0].copy()
        #
        #for epoch in [1, 2]:

        # Loop through combEpochsDF, peforming fits for
        #  epochs that have data in this inputFile... 
        for index, row in combEpochsDF.iterrows():

            epoch = row['name']
            mjd_start = row['mjd_start']
            mjd_end = row['mjd_end']

            print index, epoch, mjd_start, mjd_end
        
            # Extract data from inputFile that lie within this epoch...
            df = df_orig[ ( (df_orig.mjd >= mjd_start) & (df_orig.mjd < mjd_end) ) ].copy()
            # Are there at least 10 entries left?  If not, skip...
            if df.mjd.size < 10: 
                if args.verbose > 1:
                    print """WARNING:  Only %d entries in epoch %s...  skipping...""" % \
                          (df.mjd.size, epoch)
                del df
                continue
            # Are there at least 2 unique values of the airmass left?  If not, skip...
            if df.airmass.unique().size < 2: 
                if args.verbose > 1:
                    print """WARNING:  Fewer than 2 unique airmasses for entries in epoch %s...  skipping...""" % \
                          (epoch)
                del df
                continue
            # Maybe add other constraints, multiple unique time stamps or expnums, 
            #  to avoid problems with fits.


            #
            # Perform iterative sigma-clipping fit for this epoch...
            #

            # Create initial (and generous)  mask...
            #mjd0 = 1.0*int(df['mjd'].min())
            mjd0 = 1.0*int(mjd_start+0.5)
            df['dmjd'] = df['mjd'] - mjd0
            print "mjd0:", mjd0
            mask1 = ( np.abs( df['dmag'] - df['dmag'].median() ) < 1.0 )
            df['vigintile'] = pd.qcut( df['dmag'], 20, labels=False )
            mask2 = (df['vigintile']>0)
            mask = mask1 & mask2

            # Sigma-clipping parameters...
            dmag = 'dmag'
            nsigma = 3.0
            niter = 5

            for i in range(niter):

                iiter = i + 1
                print """   iter%d...""" % ( iiter )

                # make a copy of original df, overwritting the old one...
                df = df[mask].copy()

                # Perform fit...
                outputFile = """%s/%s.%s.fitlog.iter%d.txt"""% (qaDir, baseNameNoExt, epoch, iiter)
                p,rms = aTmCamTestFit(df.loc[:,'dmjd'], df.loc[:,'airmass'], df.loc[:,dmag], outputFile)

                # Append a line to outputFile...
                fout = open(outputFile, 'a')
                if p[1] < 0.0: 
                    sign1 = ''
                else:
                    sign1 = '+'
                if p[2] < 0.0: 
                    sign2 = ''
                else:
                    sign2 = '+'
                outputLine = """%s (MJD%.4f-%.4f) Fit Equation (%s-band):  dmag = %.4f %s %.4f*(MJD-%.4f) %s %.4f*AIRMASS    (RMS: %.4f)\n""" % \
                    (epoch, mjd_start, mjd_end, band, p[0], sign1, p[1], mjd0, sign2, p[2], rms)
                print outputLine
                fout.write(outputLine)
                fout.close()

                # If this is the final iteration, copy the outputFile to the "final" output file
                if iiter == niter:
                    finalOutputFile = """%s/%s.%s.fitlog.final.txt"""% (qaDir, baseNameNoExt, epoch)
                    shutil.copy (outputFile, finalOutputFile)
                
                # Calculate residuals and estimate new mask...
                df.loc[:,'res'] = residuals(p,df.loc[:,'dmjd'],df.loc[:,'airmass'],df.loc[:,dmag])
                stddev = df['res'].std()
                mask = (np.abs(df.res)< nsigma*stddev)
    
                # Output plot of resdiduals vs. airmass...
                x=df['airmass']
                y=df['res']
                fig, axs = plt.subplots(ncols=1)
                ax=axs
                hb = ax.hexbin(x, y, gridsize=100, bins='log', cmap='inferno')
                #ax.axis([xmin, xmax, ymin, ymax])
                title = """File=%s; Epoch:%s; RMS=%.3f (Iter%d)""" % (baseName, epoch, stddev, iiter)
                ax.set_title(title)
                ax.set_xlabel("airmass")
                ylabel="""%sres"""  % (dmag)
                ax.set_ylabel(ylabel)
                cb = fig.colorbar(hb, ax=ax)
                cb.set_label('log10(N)')
                plt.grid(True)
                ax.grid(color='white')
                outputFile = """%s/%s.%s.res_vs_airmass.iter%d.png"""% (qaDir, baseNameNoExt, epoch, iiter)
                fig.savefig(outputFile)
                plt.close()
                del x
                del y

                # If this is the final iteration, copy the outputFile to the "final" output file
                if iiter == niter:
                    finalOutputFile = """%s/%s.%s.res_vs_airmass.final.png"""% (qaDir, baseNameNoExt, epoch)
                    shutil.copy (outputFile, finalOutputFile)

                # Output plot of resdiduals vs. MJD...
                x=df['mjd']
                y=df['res']
                fig, axs = plt.subplots(ncols=1)
                ax=axs
                hb = ax.hexbin(x, y, gridsize=100, bins='log', cmap='inferno')
                title = """File=%s; Epoch:%s; RMS=%.3f (Iter%d)""" % (baseName, epoch, stddev, iiter)
                ax.set_title(title)
                ax.set_xlabel("MJD")
                ylabel="""%sres"""  % (dmag)
                ax.set_ylabel(ylabel)
                cb = fig.colorbar(hb, ax=ax)
                cb.set_label('log10(N)')
                plt.grid(True)
                ax.grid(color='white')
                outputFile = """%s/%s.%s.res_vs_mjd.iter%d.png"""% (qaDir, baseNameNoExt, epoch, iiter)
                fig.savefig(outputFile)
                plt.close()
                del x
                del y

                # If this is the final iteration, copy the outputFile to the "final" output file
                if iiter == niter:
                    finalOutputFile = """%s/%s.%s.res_vs_mjd.final.png"""% (qaDir, baseNameNoExt, epoch)
                    shutil.copy (outputFile, finalOutputFile)


                # Output histogram of resdiduals...
                ax=df.hist('res', grid=True, bins=100)
                ax=df['res'].hist(grid=True, bins=100)
                title = """File=%s; Epoch:%s; RMS=%.3f (Iter%d)""" % (baseName, epoch, stddev, iiter)
                ax.set_title(title)
                fig = ax.get_figure()
                outputFile = """%s/%s.%s.res_hist.iter%d.png"""% (qaDir, baseNameNoExt, epoch, iiter)
                fig.savefig(outputFile)
                plt.close()

                # If this is the final iteration, copy the outputFile to the "final" output file
                if iiter == niter:
                    finalOutputFile = """%s/%s.%s.res_hist.final.png"""% (qaDir, baseNameNoExt, epoch)
                    shutil.copy (outputFile, finalOutputFile)


            #endfor

            del df

        #endfor

    #endfor


#--------------------------------------------------------------------------
# Define some functions for fitting dmag vs. (dmjd,airmass)
#
# These functions are based on a scripts found at 
# http://linuxgazette.net/115/andreasen.html (by Anders Andreasen)
# and at
# http://www.phy.uct.ac.za/courses/python/examples/fitresonance.py (University of Cape Town)
#
#--------------------------------------------------------------------------
# Parametric function:  
#  p is the parameter vector; 
def fp(p,dmjd_array,airmass_array):
    return p[0] + p[1]*dmjd_array + p[2]*airmass_array

#--------------------------------------------------------------------------
# Error function:
def residuals(p,dmjd_array,airmass_array,dmag_array):
    err = (dmag_array-fp(p,dmjd_array,airmass_array))
    return err

#--------------------------------------------------------------------------
# Fitting code:
def aTmCamTestFit(dmjd_array, airmass_array, dmag_array, outputFile=None):

    if outputFile is not None:
        # Should really verify if outputFile is a valid path... 
        fout = open(outputFile, 'w')

    # Calculate the median of dmag for use as an initial guess
    # for the overall zeropoint offset..
    mdn = np.median( dmag_array, None )

    # Parameter names
    pname = (['a_0', 'a_1', 'k'])

    # Initial parameter values
    p0 = [mdn, 0.0, 0.0]

    outputLine = """\nInitial parameter values:  %s\n""" % (p0)
    print outputLine,
    if outputFile is not None:
        fout.write(outputLine)
    #print 
    #print 'Initial parameter values:  ', p0


    # Perform fit

    p,cov,infodict,mesg,ier = leastsq(residuals, p0, args=(dmjd_array, airmass_array, dmag_array), maxfev=10000, full_output=1)

    if ( ier>=1 and ier <=4):
        outputLine = 'Converged\n'
        print outputLine,
        if outputFile is not None:
            fout.write(outputLine)
        #print "Converged"
    else:
        outputLine = """Not converged\n%s\n""" (mesg)
        print outputLine,
        if outputFile is not None:
            fout.write(outputLine)
        #print "Not converged"
        #print mesg
        raise ValueError('scipy.optimize.leastsq fit did not converge.')   

    # Calculate some descriptors of the fit 
    # (similar to the output from gnuplot 2d fits)

    chisq=sum(infodict['fvec']*infodict['fvec'])
    dof=len(dmag_array)-len(p)
    rms=math.sqrt(chisq/dof)
    
    outputLine1 = """Converged with chi squared %g\n""" % (chisq)
    outputLine2 = """degrees of freedom, dof %d\n""" % (dof)
    outputLine3 = """RMS of residuals (i.e. sqrt(chisq/dof)) %g\n""" % (rms)
    outputLine4 = """Reduced chisq (i.e. variance of residuals) %g\n""" % (chisq/dof)
    outputLine = outputLine1+outputLine2+outputLine3+outputLine4+'\n'
    print outputLine,
    if outputFile is not None:
        fout.write(outputLine)
    #print "Converged with chi squared ",chisq
    #print "degrees of freedom, dof ", dof
    #print "RMS of residuals (i.e. sqrt(chisq/dof)) ", rms
    #print "Reduced chisq (i.e. variance of residuals) ", chisq/dof
    #print


    # uncertainties are calculated as per gnuplot, "fixing" the result
    # for non unit values of the reduced chisq.
    # values at min match gnuplot
    print "Fitted parameters at minimum, with 68% C.I.:"
    for i,pmin in enumerate(p):
        outputLine = """%-10s %13g +/- %13g   (%5f percent)\n""" % \
            (pname[i],pmin,math.sqrt(cov[i,i])*math.sqrt(chisq/dof),\
                 100.*math.sqrt(cov[i,i])*math.sqrt(chisq/dof)/abs(pmin))
        print outputLine,
        if outputFile is not None:
            fout.write(outputLine)
        #print "%-10s %13g +/- %13g   (%5f percent)" % (pname[i],pmin,math.sqrt(cov[i,i])*math.sqrt(chisq/dof),100.*math.sqrt(cov[i,i])*math.sqrt(chisq/dof)/abs(pmin))
    #print

    # correlation matrix close to gnuplot
    outputLine = """\nCorrelation matrix:\n"""
    outputLine = outputLine + "               "
    for i in range(len(pname)): 
        outputSnippet = "%-10s" % (pname[i])
        outputLine = outputLine + outputSnippet
    outputLine = outputLine + '\n'
    for i in range(len(p)):
        outputSnippet = "%-10s" % (pname[i])
        outputLine = outputLine + outputSnippet
        for j in range(i+1):
	    outputSnippet = "%10f" % (cov[i,j]/math.sqrt(cov[i,i]*cov[j,j]))
            outputLine = outputLine + outputSnippet
        #endfor
        outputLine = outputLine + '\n'
    #endfor
    outputLine = outputLine + '\n\n\n'
    print outputLine, 
    if outputFile is not None:
        fout.write(outputLine)
    

    #print "Correlation matrix:"
    #print "               ",
    #for i in range(len(pname)): print "%-10s" % (pname[i],),
    #print
    #for i in range(len(p)):
    #    print "%-10s" % pname[i],
    #    for j in range(i+1):
    #        print "%10f" % (cov[i,j]/math.sqrt(cov[i,i]*cov[j,j]),),
    #    #endfor
    #    print
    ##endfor
    #
    #print
    #print
    #print
    
    if outputFile is not None:
        fout.close()

    return p, rms


#--------------------------------------------------------------------------

if __name__ == "__main__":
    main()

#--------------------------------------------------------------------------
