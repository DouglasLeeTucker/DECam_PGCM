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
