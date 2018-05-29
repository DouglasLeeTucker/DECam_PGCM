#!/usr/bin/env python

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
    parser.add_argument('--planFile', help='name of the input plan file', default='interstellar_Rband.par')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    status = interstellar_Rband(args)

    return status

#--------------------------------------------------------------------------

def interstellar_Rband(args):

    import os

    import reddening
    import spectrum
    import paramFile

    # Read contents of planFile into python dictionary...
    planDict = paramFile.readParamFile(args.planFile, args.verbose)

    # Check that all the required keywords were found...
    requiredKeywordList = ['bandList',
                           'bandpassFile',
                           'Rv',
                           'c',
                           'N',
                           'ebv']

    flag = 0
    for requiredKeyword in requiredKeywordList:
        if requiredKeyword not in planDict:
            print """Required keyword '%s' not in planFile %s""" % (requiredKeyword, args.planFile)
            flag = flag + 1
    if flag > 0:
        print 'Returning with error code 1 now...'
        return 1


    # Extract the relevant info from the planDict...
    
    #  Extract the bandList...
    #  (Default is u,g,r,i,z,Y.)
    if  planDict['bandList'].lower() == 'default':
        bandList = ['u','g','r','i','z','Y']
    else:
        bandList = planDict['bandList']
    bandList = bandList.split(',')
    if args.verbose > 0:
        print 'bandList: ', bandList

    #  Extract the reddening law...
    #  (Default is fitzpatrick99.)
    if planDict['law'].lower() == 'default':
        law = 'fitzpatrick99'
    else:
        law = planDict['law']
    if args.verbose > 0:
        print 'law: ', law

    #  Extract Rv...
    #  (Default is 3.1)
    if planDict['Rv'].lower() == 'default':
        Rv = 3.1
    else:
        Rv = float(planDict['Rv'])    
    if args.verbose > 0:
        print 'Rv: ', Rv

    #  Extract the value for "E(B-V)/A(1micron)"...
    #  (Default is 1.319; see Appendix of Schlafly et al. 2010.)
    if planDict['c'].lower() == 'default':
        c = 1.319 # Value for 
    else:
        c = float(planDict['c'])    
    if args.verbose > 0:
        print 'c: ', c
    
    #  Extract the normalization factor N...
    #  (Default is 0.78 for SFD98 maps and 1.00 for Planck13 and Lenz17 maps; 
    #   see Appendices of Schlafly et al. 2010 and of Schlafly & Finkbeiner 2011.)
    if planDict['N'].lower() == 'default':
        N = 0.78
    else:
        N = float(planDict['N'])    
    if args.verbose > 0:
        print 'N: ', N
    
    #  Extract the value for "E(B-V)" to be used in the full calculation...
    #  (Default is 1.0e-6, the very low-exnction limit.)
    if planDict['ebv'].lower() == 'default':
        ebv = 1.0e-6
    else:
        ebv = float(planDict['ebv'])    
    if args.verbose > 0:
        print 'ebv: ', ebv

    
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
    if os.path.isfile(bandpassFile)==False:
        print """bandpassFile %s does not exist...""" % (bandpassFile)
        print 'Returning with error code 1 now...'
        return 1
    if args.verbose > 0:
        print 'bandpassFile: ', bandpassFile


    df_resp = pd.read_csv(bandpassFile, comment='#')
    rlaw = reddening.getReddeningLaw(inv=True)
    flam = spectrum.getSpectrumFlat()

    # For purposes of debugging...
    if args.verbose > 0:
        w = [2600, 2700, 4110, 4670, 5470, 6000, 12200, 26500, 100000]
        print w
        #print r(w)
        print rlaw(1./np.array(w))
        print flam(w)


    delta_wavelength = 2.5 # angstroms
    # Temporary kludge for GALEX fuv, nuv bandpasses:
    if ( df_resp.LAMBDA.min() > 1000.) and  (df_resp.LAMBDA.max() < 3600.):
        wavelength_array = np.arange(1000, 3600 + delta_wavelength, delta_wavelength)
    else:
        wavelength_array = np.arange(2600, 25000 + delta_wavelength, delta_wavelength)
    flam_array = flam(wavelength_array)
    rlaw_array = rlaw(1. / wavelength_array)
    #rlaw_array = rlaw(wavelength_array)
    # Normalize the reddening law to its value at 1micron.
    rlaw_array /= rlaw(1. / 1.e4) 
    #rlaw_array /= rlaw(1.e4)

    for band in bandList:
    
        response = interpolate.interp1d(df_resp['LAMBDA'], df_resp[band], 
                                            bounds_error=False, fill_value=0., 
                                            kind='linear')
        response_array = response(wavelength_array)


        # Calculate R_band using the exact (non-linear) equation, 
        #  where "exact" is pretty loosely defined, given that this 
        #  is a numerical integration...
        #  (via Schlegel & Finkbeiner, Keith, and Eli)...
        exact_value = rb_exact(wavelength_array, response_array, flam_array, rlaw_array, c, N, ebv)

        # Calculate R_band using the approximate (linearized) equation (via Eli)...
        approx_value = rb_approx(wavelength_array, response_array, flam_array, rlaw_array, c, N)

        # Print out results...
        print '%s, R_f[\"exact\"] = %.4f, R_f[approx] = %.4f'%(band, exact_value, approx_value)


# Calculate R_band using the exact (non-linear) equation...
def rb_exact(wavelength_array, response_array, flam_array, rlaw_array, c, N, ebv):

    # See Equations 2&3 of Rykoff, Burke, et al. "Linearized Reddening Corrections",
    #   but see also the Appendixes of Schlafly et al. 2010 and Schlafly & Finkbeiner 2011.
    # This is essentially what is calculated in Table 6 of Schlafly et al. 2011
    #  (i.e., "Delta m_b / E(B-V)", except here we have used a flat Flam spectrum 
    #  rather than a Teff=7000K spectrum for the source.)
    # Recall that, here, rlaw_array is the reddening law normalized at 1micron.
    # Further, "c * N * ebv" is the N times the extinction @ 1micron [in mags] 
    #  according to E(B-V) map. See Schlafly et al. 2011 Eq. A1 and related material.
    a_array = c * N * ebv * rlaw_array 
    numerator = np.sum(wavelength_array * response_array * flam_array * 10.**(-0.4 * a_array))
    denominator = np.sum(wavelength_array * response_array * flam_array ) 
    extinction = -2.5 * (np.log10(numerator) - np.log10(denominator))
    exact_value = (extinction / ebv)  

    return exact_value


# Calculate R_band the using the linearized approximation...
def rb_approx(wavelength_array, response_array, flam_array, rlaw_array, c, N):

    # See Equation 10 of Rykoff, Burke, et al. "Linearized Reddening Corrections",
    #  which is the first-order [linear] approximation, in the limit of small 
    #  extinctions, of Equations 2&3 of the same document...
    # Recall that, here, rlaw_array is the reddening law normalized at 1micron.
    # The "exact" method and the "approximate" method result in the same values 
    #   for Rlam in the limit of small E(B-V) [hence the use of ebv=1.0e-6 above]...
    numerator = np.sum(wavelength_array * response_array * flam_array * rlaw_array)
    denominator = np.sum(wavelength_array * response_array * flam_array )
    approx_value = c * N * numerator / denominator

    return approx_value


#--------------------------------------------------------------------------

if __name__ == "__main__":
    main()

#--------------------------------------------------------------------------
