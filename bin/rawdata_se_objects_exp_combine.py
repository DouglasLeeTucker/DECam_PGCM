#!/usr/bin/env python
"""
    rawdata_se_objects_exp_combine.py

    Example:
    
    rawdata_se_objects_exp_combine.py --help

    rawdata_se_objects_exp_combine.py --inputSEObjFile seobjfile.csv
                                      --inputExpFile expfile.csv
                                      --outputFile outputfile.csv
                                      --verbose 2
    
    """

##################################

def main():

    import os
    import argparse
    import time

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputSEObjFile', help='name of CSV file containing the se_object info', default='seobjfile.csv')
    parser.add_argument('--inputExpFile', help='name of CSV file containing the exposure info', default='expfile.csv')
    parser.add_argument('--outputFile', help='name of output CSV file', default='outputfile.csv')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    # Execute method...
    status = rawdata_se_objects_exp_combine(args)


##################################
# rawdata_se_objects_exp_combine
#

def rawdata_se_objects_exp_combine(args):

    import os
    import datetime
    import numpy as np
    import pandas as pd

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'rawdata_se_objects_exp_combine'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    inputSEObjFile = args.inputSEObjFile
    inputExpFile = args.inputExpFile
    outputFile = args.outputFile

    # Read in inputExpFile...
    if os.path.isfile(inputExpFile)==False:
        print """Input file %s does not exist...""" % (inputExpFile)
        print """Exiting rawdata_se_objects_exp_combine method with return code 1"""
        return 1
    if args.verbose > 1:
        print """Reading in input file %s...""" % (inputExpFile)
        print datetime.datetime.now()
    inputExpFileDF = pd.read_csv(inputExpFile) 

    # Read in contents of inputSEObjFile in chunks of 1 million lines
    #  and combine with contents of inputExpFile...
    if os.path.isfile(inputSEObjFile)==False:
        print """Input file %s does not exist...""" % (inputSEObjFile)
        print """Exiting rawdata_se_objects_exp_combine method with return code 1"""
        return 1
    if args.verbose > 1:
        print """Reading in input file %s...""" % (inputSEObjFile)
        print datetime.datetime.now()

    # If output file of the same name already exists,
    #  rename the original...
    if os.path.isfile(outputFile):
        newFileName = outputFile+'~'
        if args.verbose > 1:
            print """Renaming original %s to %s...""" % \
                (outputFile, newFileName)
        os.rename(outputFile,newFileName)

    chunksize = 10 ** 6
    ichunk = 1
    for chunk in pd.read_csv(inputSEObjFile, chunksize=chunksize, low_memory=False):

        if args.verbose > 1:
            print """Working on chunk %d...""" % (ichunk)

        inputSEObjFileDF = chunk
    
        # Combine the two dataframes...
        # To drop columns with duplicate names, we make use of 
        #  recommendations from an e-mail from S. Allam (27 July 2016)
        #  and from these two stackoverflow questions:
        #  http://stackoverflow.com/questions/27313647/merging-two-pandas-dataframes-results-in-duplicate-columns 
        #  and http://stackoverflow.com/questions/19125091/pandas-merge-how-to-avoid-duplicating-columns
        if args.verbose > 1:
            print """Combining the two dataframes"""
            print datetime.datetime.now()

        combinedDF = pd.merge(inputExpFileDF, inputSEObjFileDF, \
                                  on=['EXPNUM'], \
                                  how='inner', \
                                  suffixes=('','_y')).reset_index(drop=True)
        to_drop = [colname for colname in combinedDF if colname.endswith('_y')]
        combinedDF.drop(to_drop, axis=1, inplace=True)   

        # If outputFile does not exist, create it; 
        #  otherwise append to it...
        if args.verbose > 1:
            print """Writing out combined dataframe to %s""" % (outputFile)
            print datetime.datetime.now()
        if not os.path.isfile(outputFile):
            combinedDF.to_csv(outputFile, index=False)
        else: 
            combinedDF.to_csv(outputFile, mode='a', index=False, header=False)

        # Delete chunk dataframe...
        del inputSEObjFileDF

        #increment ichunk
        ichunk += 1


    # Final clean up...
    del combinedDF
    del inputExpFileDF

    if args.verbose > 1:
        print """Finished with rawdata_se_objects_exp_combine..."""
        print datetime.datetime.now()

    if args.verbose>0: print

    return 0


##################################

if __name__ == "__main__":
    main()

##################################


