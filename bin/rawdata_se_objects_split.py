#!/usr/bin/env python
"""
    rawdata_se_objects_split.py

    Example:
    
    rawdata_se_objects_split.py --help

    rawdata_se_objects_split.py --inputFileListFile inputfilelist.csv
                                --outputFileListFile outputfilelist.csv
                                --verbose 2
    
    """

##################################

def main():

    import os
    import argparse
    import time

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputFileListFile', help='name of CSV file containing list of input files', default='inputfilelist.csv')
    parser.add_argument('--outputFileListFile', help='name of CSV file containing list of filter bands and output files', default='outputfilelist.csv')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    # Execute method...
    status = rawdata_se_objects_split(args)


##################################
# rawdata_se_objects_split
#
# Based on sepBands from y2a1_tertiaries.py
#

def rawdata_se_objects_split(args):

    import os
    import datetime
    import numpy as np
    import pandas as pd
    from astropy.table import Table

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'rawdata_se_objects_split'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    # Read in inputFileListFile...
    inputFileListFile = args.inputFileListFile
    if os.path.isfile(inputFileListFile)==False:
        print """Input filelist file %s does not exist...""" % (inputFileListFile)
        print """Exiting rawdata_se_objects_split method with return code 1"""
        return 1
    inputFileListDF = pd.read_csv(inputFileListFile, 
                                  header=None, 
                                  names=['FILENAME'], 
                                  comment='#')
    inputFileListDF['FILENAME'] = inputFileListDF['FILENAME'].str.strip()
    inputFileListSeries = inputFileListDF['FILENAME']
    inputFileList = inputFileListSeries.values.tolist()

    if args.verbose>1:  
        print 'Input file list:'
        print inputFileList


    # Read in outputFileListFile and convert it to a python 
    #  dictionary, ensuring there are no extraneous white spaces
    #  in the file names listed...
    outputFileListFile = args.outputFileListFile
    if os.path.isfile(outputFileListFile)==False:
        print """Output filelist file %s does not exist...""" % (outputFileListFile)
        print """Exiting rawdata_se_objects_split method with return code 1"""
        return 1
    outputFileListDF = pd.read_csv(outputFileListFile, 
                                  header=None, 
                                  names=['BAND','FILENAME'], 
                                  index_col='BAND', 
                                  comment='#')
    outputFileListDF['FILENAME'] = outputFileListDF['FILENAME'].str.strip()
    outputFileListSeries = outputFileListDF['FILENAME']
    outputFileListDict = outputFileListSeries.to_dict()

    # Also, grab the band list from the outputFileListFile series...
    bandList = outputFileListSeries.index.values.tolist()

    if args.verbose>1:  
        print 'Output file list dictionary:'
        print outputFileListDict
        print 'Band list:'
        print bandList


    # Loop through inputFileList...
    firstFile=True
    for inputFile in inputFileList:

        if args.verbose > 1:
            print """Working on input file %s...""" % inputFile
            print datetime.datetime.now()

        # Read in file...
        t = Table.read(inputFile)

        # Convert astropy Table to pandas dataframe...
        df = t.to_pandas()

        # Verify that 'BAND' is a column in the dataframe; 
        #  otherwise, skip...
        if 'BAND' not in df.columns:
            print """Could not find 'BAND' in header of %s...  Skipping""" \
                % (inputFile)
            del df
            continue

        # Verify that 'FILENAME' is a column in the dataframe; 
        #  otherwise, skip...
        if 'FILENAME' not in df.columns:
            print """Could not find 'FILENAME' in header of %s...  Skipping""" \
                % (inputFile)
            del df
            continue

        # Trim leading and trailing white space from the FILENAME column...
        df['FILENAME'] = df['FILENAME'].str.strip()

        # If this is the first (valid) file, create initial 
        #   output files (empty except for the CSV header)...
        if firstFile is True:
            for band in bandList:
                outputFile = outputFileListDict[band]
                # Create a mask with all entries equal to False...
                mask = pd.Series(np.zeros(df.BAND.size, dtype=bool))        
                df[mask].to_csv(outputFile,index=False)
            firstFile = False

        # Loop through band list, appending the rows from
        #  each band to the appropriate output file...
        for band in bandList:
            outputFile = outputFileListDict[band]
            mask = (df.BAND == band)
            # Faster if we move the "open" to outside the loop?:
            with open(outputFile, 'a') as f:
                df[mask].to_csv(f, index=False, header=False)
            f.close()

        # Clean up some space before moving to next file...
        del df
        #del t

        if args.verbose > 1:
            print datetime.datetime.now()


    if args.verbose>0: print

    return 0


##################################

if __name__ == "__main__":
    main()

##################################


