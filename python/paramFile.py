# Reads first 2 columns from a whitespace-delimited text file 
#  into a python dictionary, taking the 1st column as the 
#  keyword and the 2nd column as the value for that keyword.
#  Lines having '#' as the first character and lines having
#  fewer than 2 whitespace-delimited columns are ignored.
#
# WARNING:  all values are returned as strings.
#
def readParamFile(fileName, verbose=0):

    import numpy as np
    import os
    import sys
    
    # Read in fileName...
    if os.path.isfile(fileName)==False:
        print """File %s does not exist...""" % (fileName)
        return 1

    paramDict = {}
    with open(fileName,'r') as fin:

        for lin in fin:

            if lin[0] == '#': continue
            lins=lin.strip().split()
            if len(lins) < 2: continue
            keyword = lins[0].strip()
            value = lins[1].strip()
            if verbose > 0: print keyword+": "+value
            paramDict[keyword] = value            
            
    return paramDict
