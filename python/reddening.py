# Wrapper around Kyle Barbary's "extinction" module
#  (https://github.com/kbarbary/extinction).
# Wrapper function based on Keith Becthol's
#  fiducial_r_f_values_ugrizy_r_v_scan.py code.
def getReddeningLaw(law='fitzpatrick99',Rv=3.1,inv=False):

    import numpy as np
    from scipy import interpolate
    import extinction

    # Wavelength ranges (lambda_min - lambda_max) of the various reddening laws 
    #  (in Angstroms)...
    lambda_min = {'ccm89':          1250., 
                  'odonnell94':     1250., 
                  'calzetti00':     1200., 
                  'fitzpatrick99':   910., 
                  'fm07':            910.}
    lambda_max = {'ccm89':         33000.,
                  'odonnell94':    33000.,
                  'calzetti00':    22000.,
                  'fitzpatrick99': 60000.,
                  'fm07':          60000.}
    # We can extract the list of supported reddening laws by
    #  grabbing those that are keys within the lambda_min dictionary...
    supported_laws = lambda_min.keys()

    # If reddening law not in the the list of supported reddening laws,
    #  return an Exception...
    if law not in supported_laws: 
        print """Un-supported reddening law:  %s""" % (law)
        print 'Supported reddening laws are: ', supported_laws 
        print 'Returning exception'
        return Exception

    # Calculate and return the reddening law in either
    #  inverse wavelength form (inv=True) or in wavelength
    #  form (inv=False)...
    if inv==True:

        # Use inverse microns to call to "extinction" module
        #  and return reddening law in inverse Angstroms...

        # Calculate inverse wavelengths...
        x_lambda_min = 1.0e4/lambda_max[law]
        x_lambda_max = 1.0e4/lambda_min[law]
        x_micron = np.linspace(x_lambda_min, x_lambda_max, 2000) # microns^-1
        x_angstrom = x_micron * 1.0e-4 # Convert from microns^-1 to Anstroms^-1

        # Call appropriate reddening law function...
        if law == 'ccm89':
            r_array = Rv*extinction.ccm89(x_micron, 1.0, Rv, unit='invum')
        elif law == 'odonnell94':
            r_array = Rv*extinction.odonnell94(x_micron, 1.0, Rv, unit='invum')
        elif law == 'calzetti00':
            r_array = Rv*extinction.calzetti00(x_micron, 1.0, Rv, unit='invum')
        elif law == 'fitzpatrick99':
            r_array = Rv*extinction.fitzpatrick99(x_micron, 1.0, Rv, unit='invum')
        elif law == 'fm07':
            r_array = Rv*extinction.fm07(x_micron, 1.0, unit='invum')

        # Create interpolation function for reddening law...
        r = interpolate.interp1d(x_angstrom, r_array, 
                                 bounds_error=False, fill_value=0., kind=3)

    else:

        # Use Angstroms to call to "extinction" module
        #  and return reddening law in Angstroms...

        # Create wavelength array...
        angstrom = np.logspace(np.log10(lambda_min[law]), np.log10(lambda_max[law]), 2000)

        # Call appropriate reddening law function...
        if law == 'ccm89':
            r_array = Rv*extinction.ccm89(angstrom, 1.0, Rv, unit='aa')
        elif law == 'odonnell94':
            r_array = Rv*extinction.odonnell94(angstrom, 1.0, Rv, unit='aa')
        elif law == 'calzetti00':
            r_array = Rv*extinction.calzetti00(angstrom, 1.0, Rv, unit='aa')
        elif law == 'fitzpatrick99':
            r_array = Rv*extinction.fitzpatrick99(angstrom, 1.0, Rv, unit='aa')
        elif law == 'fm07':
            r_array = Rv*extinction.fm07(angstrom, 1.0, unit='aa')

        # Create interpolation function for reddening law...
        r = interpolate.interp1d(angstrom, r_array, 
                                 bounds_error=False, fill_value=0., kind='linear')

    # Return interpolation fucntion...
    return r


# Based on code from Keith Bechtol's synthesize_locus.py code...
#  Calculates and returns the Fitzpatrick 1999 reddening law
#  (for Rv=3.1) in inverse Angstroms...
def getReddeningOrig():

    import numpy as np
    from scipy import interpolate

    # Fitzpatrick 1999
    wavelength, a = zip(*[[2600, 6.591],
                          [2700, 6.265],
                          [4110, 4.315],
                          [4670, 3.806],
                          [5470, 3.055],
                          [6000, 2.688],
                          [12200, 0.829],
                          [26500, 0.265],
                          [1000000, 0.]])
    
    r = interpolate.interp1d(1. / np.array(wavelength), a, 
                             bounds_error=False, fill_value=0., kind=3) # NORMAL

    return r
