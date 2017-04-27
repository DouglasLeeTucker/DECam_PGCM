#!/usr/bin/env python
"""
    decamInfo.py

    
    """

##################################
# getDecamFocalPlaneInfoAsNDArray
#

def getDecamFocalPlaneInfoAsNDArray():

    import numpy as np
    import os

    # Grab path and name of decam.txt file in the PGCM data directory...
    # Is there a better way to do this?
    # See also http://stackoverflow.com/questions/779495/python-access-data-in-package-subdirectory
    # Absolute path for the directory containing this module:
    moduleAbsPathName = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    #decamInfoFile = moduleAbsPathName+'/data/decam.txt'
    decamInfoFile = os.path.join(moduleAbsPathName, "data", "decam.txt")

    ndarray = np.genfromtxt(decamInfoFile,dtype=None,names=['ccdnum','detpos','xpos','ypos','xsize','ysize','DAboard'])

    return ndarray



##################################
# getDecamFocalPlaneInfoAsDataFrame
#

# NOT YET TESTED!:
def getDecamFocalPlaneInfoAsDataFrame():

    import numpy as np
    import pandas as pd
    import os

    ndarray = getDecamFocalPlaneInfoAsNDArray()

    # Convert whole ndarray to pandas DataFrame...
    df_decam = pd.DataFrame(ndarray)
    # Remove those annoying ccdnum==0 ccds...
    df_decam = df_decam[df_decam.ccdnum != 0]
    # Sort by ccdnum...
    df_decam.sort('ccdnum', ascending=True, inplace=True)
    # Replace original index with the ccdnum...
    df_decam.set_index('ccdnum',inplace=True)    

    return df


##################################


