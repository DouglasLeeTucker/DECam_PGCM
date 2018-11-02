# Unless otherwise noted, fluxes are assumed to be Flam vs. 
#  wavelength in Angstroms...

# Return a spectrum that is constant vs. wavelength...
#  (Based on code from Keith Bechtol's synthesize_locus.py.)
def getSpectrumFlat(flux=1.0):

    import numpy as np
    from scipy import interpolate    

    data = {'wavelength': [100., 100000.],
            'flux': [flux, flux]}

    f = interpolate.interp1d(data['wavelength'], data['flux'], 
                             bounds_error=False, fill_value=0., 
                             kind='linear')

    return f

def getSpectrumSynphot(synphotFileName, fluxFactor=1.0):

    import numpy as np
    from scipy import interpolate 
    from astropy.io import fits
    from astropy.table import Table
    import sys

    try:

        hdulist = fits.open(synphotFileName)
        t = Table.read(hdulist[1])
        hdulist.close()

    except IOError:

        print """Could not read %s""" % synphotFileName
        sys.exit(1)


    wave = t['WAVELENGTH'].data.tolist()
    #flam = fluxFactor*t['FLUX'].data.tolist()   
    flam = t['FLUX'].data.tolist()   
    data = {'wavelength': wave, 'flux': flam}

    f = interpolate.interp1d(data['wavelength'], data['flux'], 
                             bounds_error=False, fill_value=0., 
                             kind='linear')

    return f

def getSpectrumBlackBody(Teff, a=1.00e-23):

    # Here, "a" is the flux normalization factor from Eq. 1 of 
    #  Suzuki & Fukugita (2018, AJ, 156, 219).  The default 
    #  value we use here (a=1.00e-23) is equivalent to 
    #  a blackbody at 100 pc with a radius of 5506 km (i.e.,  
    #  about the size of relatively small -- but not inconceivably 
    #  small -- white dwarf).
    # 
    # Note:  pysynphot's BlackBody Spectrum object is normalized to
    #  a spherical blackbody with the same radius as the Sun 
    #  (Rsun = 695700 km) at a distance of 1 kpc (d = 1 kpc = 3.086e+16 km).
    #  This, according to Eq. 2 of Suzuki & Fukugita (2018), if the 
    #  equivalent of a Suzuki & Fukugita flux normalization factor of 
    #     a = pi*(Rsun/d)**2 = 1.5966e-21
    #       = math.pi*(695700./3.086e+16)*(695700./3.086e+16)
    #       = 1.5966e-21
    # 
    # Thus, to get pysynphot to yield the same blackbody spectrum
    #  as the current function, getSpectrumBlackBody, one would 
    #  need to do the following:
    #     a_syn = math.pi*(695700./3.086e+16)*(695700./3.086e+16)
    #     bb = (a/a_syn)*S.BlackBody(Teff)
    #     bb.convert('flam')
    #  where a is the Suzuki & Fukugita flux normalization 
    #  factor for the blackbody star in question (see their Table 1).


    import numpy as np
    from scipy import interpolate 

    # Planck's constant:
    h = 6.62607004e-34  # m2 kg/s
    # Speed of light:
    c = 2.99792458e8    # m/s
    # Boltzmann constant:
    k = 1.38064852e-23  # m2 kg s-2 K-1

    # Wavelength range...
    wave_array = np.arange(100.,100000.)
    wave_meters_array = (1.00e-10)*wave_array
    # Factor of (1.00e-7) converts from MKS to ergs/s/cm**2/Angstrom:
    flam_array = (1.00e-7) * a * (2.*h*c*c/wave_meters_array**5) / ( np.exp(h*c/(wave_meters_array*k*Teff)) - 1. )

    # Combine wave_array and flam_array...
    wave = wave_array.tolist()
    flam = flam_array.tolist()
    data = {'wavelength': wave, 'flux': flam}

    # Create interpolation function...
    f = interpolate.interp1d(data['wavelength'], data['flux'], 
                             bounds_error=False, fill_value=0., 
                             kind='linear')

    return f

