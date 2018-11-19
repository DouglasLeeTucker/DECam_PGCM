#!/usr/bin/env python
"""
Calculate AB mags in a given set of bandpasses for a given spectrum.

Currently allowable spectrum types are:
     'flat' (Flam = constant)
     'synphot' (spectrum is contained within a synphot-like FITS file) 
     'blackbody' (blackbody spectrum, normalized a la Suzuki & Fukugita (2018)

Parameters are taken from a planFile.

Values of the parameters from the planFile are overwritten if these
parameters are included in the command-line.

"""

# Initial setup...
import numpy as np
import pandas as pd
from scipy import interpolate

#--------------------------------------------------------------------------
# Main code.

def main():

    import argparse

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--planFile', help='name of the input plan file', default='abmag.par')

    parser.add_argument('--bandList', help='comma-separated list of filter bands to consider', default='')
    parser.add_argument('--bandpassFile', help='full or relative path name of the bandpass file', default='')
    parser.add_argument('--spectrumType', help='type of spectrum to use (flat, synphot, or blackbody)', default='')
    parser.add_argument('--spectrumFile', help='full or relative path name of the spectrum file (for spectrumType=synphot)', default='')
    parser.add_argument('--Teff', help='Teff of spectrum (for spectrumType=blackbody)', default='')
    parser.add_argument('--a', help='Suzuki & Fukugita (2018) scale factor (for spectrumType=blackbody', default='')

    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    status = abmag(args)

    return status

#--------------------------------------------------------------------------

def abmag(args):

    import os

    import spectrum
    import paramFile

    # Read contents of planFile into python dictionary...
    # NOTE:  might wish to switch over to python configparser module
    #        in the future for paramFile-related input.
    planDict = paramFile.readParamFile(args.planFile, args.verbose)

    # Check that all the required keywords were found...
    requiredKeywordList = ['bandList',
                           'bandpassFile',
                           'spectrumType']

    flag = 0
    for requiredKeyword in requiredKeywordList:
        if requiredKeyword not in planDict:
            print """Required keyword '%s' not in planFile %s""" % (requiredKeyword, args.planFile)
            flag = flag + 1
    if flag > 0:
        print 'Returning with error code 1 now...'
        return 1


    # Extract the relevant info from the planDict...
    # Overwrite values from the command-line if provided...
    
    #  Extract the bandList...
    #  (Default is u,g,r,i,z,Y.)
    if  planDict['bandList'].lower() == 'default':
        bandList = ['u','g','r','i','z','Y']
    else:
        bandList = planDict['bandList']
    bandList = bandList.split(',')
    # Overwrite with value from command-line if provided...
    if args.bandList is not '':
        bandList = args.bandList.strip().split(',')
    if args.verbose > 0:
        print 'bandList: ', bandList

    #  Extract the spectrum type...
    #  (Default is flat, i.e., Flam=const..)
    if planDict['spectrumType'].lower() == 'default':
        spectrumType = 'flat'
    else:
        spectrumType = planDict['spectrumType']
    # Overwrite with value from command-line if provided...
    if args.spectrumType is not '':
        spectrumType = args.spectrumType
    if args.verbose > 0:
        print 'spectrumType: ', spectrumType

    #  Extract the name of the spectrum file...
    #  (if spectrumType=synphot)
    if spectrumType.lower() == 'synphot':
        spectrumFile = planDict['spectrumFile']
        # Overwrite with value from command-line if provided...
        if args.spectrumFile is not '':
            spectrumFile = args.spectrumFile
        if os.path.isfile(spectrumFile)==False:
            print """spectrumFile %s does not exist...""" % (spectrumFile)
            print 'Returning with error code 1 now...'
            return 1
        if args.verbose > 0:
            print 'spectrumFile: ', spectrumFile

    #  Extract the temperature and scale factor of the blackbody spectrum...
    #  (if spectrumType=blackbody)
    if spectrumType.lower() == 'blackbody':
        Teff = float(planDict['Teff'])
        a = float(planDict['a'])
        # Overwrite with values from command-line if provided...
        if args.Teff is not '':
            Teff = float(args.Teff)
        if args.a is not '':
            a = float(args.a)
        if args.verbose > 0:
            print 'Teff: ', Teff
            print 'a:    ', a


    #  Extract the name of the bandpassFile...
    #  (Default is the bandpass file in the DECam_PGCM/data directory.)
    bandpassFile = planDict['bandpassFile']
    if bandpassFile.lower() == 'default':
        # Grab path and name of the DES_STD_BANDPASSES_Y3A2_ugrizY.test.csv file
        #  in the DECam_PGCM data directory...
        #  Is there a better way to do this?
        #  See also:
        #  http://stackoverflow.com/questions/779495/python-access-data-in-package-subdirectory
        # Absolute path for the directory containing this module:
        moduleAbsPathName = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        bandpassFile = os.path.join(moduleAbsPathName, "data", "bandpasses", "DES_STD_BANDPASSES_Y3A2_ugrizY.test.csv")
    # Overwrite with value from command-line if provided...
    if args.bandpassFile is not '':
        bandpassFile = args.bandpassFile
    if os.path.isfile(bandpassFile)==False:
        print """bandpassFile %s does not exist...""" % (bandpassFile)
        print 'Returning with error code 1 now...'
        return 1
    if args.verbose > 0:
        print 'bandpassFile: ', bandpassFile


    # Read bandpass file into pandas dataframe...
    df_resp = pd.read_csv(bandpassFile, comment='#')

    # Get spectrum...
    if spectrumType.lower() == 'flat':
        flam = spectrum.getSpectrumFlat()
    elif spectrumType.lower() == 'synphot':
        flam = spectrum.getSpectrumSynphot(spectrumFile, fluxFactor=1.0)
    elif spectrumType.lower() == 'blackbody':
        flam = spectrum.getSpectrumBlackBody(Teff, a)


    #delta_wavelength = 2.5 # angstroms
    delta_wavelength = 1.0 # angstroms
    # Temporary kludge for GALEX fuv, nuv bandpasses:
    if ( df_resp.LAMBDA.min() > 1000.) and  (df_resp.LAMBDA.max() < 3600.):
        wavelength_array = np.arange(1000, 3600 + delta_wavelength, delta_wavelength)
    else:
        wavelength_array = np.arange(2600, 25000 + delta_wavelength, delta_wavelength)
    flam_array = flam(wavelength_array)

    for band in bandList:
        print """      %3s   """ % (band),
    print

    for band in bandList:
    
        response = interpolate.interp1d(df_resp['LAMBDA'], df_resp[band], 
                                            bounds_error=False, fill_value=0., 
                                            kind='linear')
        response_array = response(wavelength_array)


        # Calculate the abmag in this filter band for this spectrum... 
        abmag = calc_abmag(wavelength_array, response_array, flam_array)

        # Print out results...
        #print '%s, abmag = %.4f' % (band, abmag)
        print """  %10.4f""" % (abmag),

# Calculate the abmag in this filter band for this spectrum... 
def calc_abmag(wavelength_array, response_array, flam_array):

    c = 2.99792458e8  # m / s
    c_angstroms_per_sec = c * 1e10

    fnu_array = wavelength_array*wavelength_array*flam_array/c_angstroms_per_sec

    numerator = np.sum(fnu_array * response_array / wavelength_array)
    denominator = np.sum(response_array / wavelength_array )

    abmag = -2.5 * (np.log10(numerator) - np.log10(denominator)) - 48.60

    return abmag


#--------------------------------------------------------------------------

if __name__ == "__main__":
    main()

#--------------------------------------------------------------------------
