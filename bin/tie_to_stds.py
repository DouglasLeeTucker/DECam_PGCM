#!/usr/bin/env python
"""
    tie_to_stds.py

    Generic version of tie_to_fgcm_stds.py.

    Example:
    
    tie_to_stds.py --help

    tie_to_stds.py --inputFile inputFile.csv --outputFile outputFile.csv --verbose 2
    
    """

##################################

def main():

    import argparse
    import time

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputFile', 
                        help='name of the input CSV file', 
                        default='input.csv')
    parser.add_argument('--outputFile', 
                        help='name of the output CSV file', 
                        default='output.csv')
    parser.add_argument('--magStdColName', 
                        help='name of the column containing the standard magnitude', 
                        default='GMAG_DES_1')
    parser.add_argument('--magerrStdColName', 
                        help='name of the column in inputFile containing the standard magnitude error', 
                        default='None')
    parser.add_argument('--fluxObsColName', 
                        help='name of the column in inputFile containing the observed flux', 
                        default='FLUX_PSF_2')
    parser.add_argument('--fluxerrObsColName', 
                        help='name of the column in inputFile containing the observed fluxerr', 
                        default='FLUXERR_PSF_2')
    parser.add_argument('--aggFieldColName', 
                        help='name of the column in inputFile containing the field to aggregate', 
                        default='FILENAME_2')
    parser.add_argument('--verbose', 
                        help='verbosity level of output to screen (0,1,2,...)', 
                        default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    status = tie_to_stds(args)


##################################
# 

def tie_to_stds(args):

    import numpy as np 
    import os
    import sys
    import datetime
    import pandas as pd

    inputFile = args.inputFile
    outputFile = args.outputFile
    magStdColName = args.magStdColName
    magerrStdColName = args.magerrStdColName
    fluxObsColName = args.fluxObsColName
    fluxerrObsColName = args.fluxerrObsColName
    aggFieldColName = args.aggFieldColName

    reqColList = [magStdColName,magerrStdColName,fluxObsColName,fluxerrObsColName,aggFieldColName]

    # If magerrStdColName is 'none', remove it from reqColList...
    if magerrStdColName.lower() == 'none':  
        magerrStdColName = 'None'
        del reqColList[1]

    # Does the input file exist?
    if os.path.isfile(inputFile)==False:
        print """tie_to_stds input file %s does not exist.  Exiting...""" % (inputFile)
        return 1

    # Read inputFile into a pandas DataFrame...
    print datetime.datetime.now()
    print """Reading in %s as a pandas DataFrame...""" % (inputFile)
    dataFrame = pd.read_csv(inputFile)
    print datetime.datetime.now()
    print

    reqColFlag = 0
    colList = dataFrame.columns.tolist()
    for reqCol in reqColList:
        if reqCol not in colList:
            print """ERROR:  Required column %s is not in the header""" % (reqCol)
            reqColFlag = 1
    if reqColFlag == 1:
        print """Missing required columns in header of %s...  Exiting now!""" (inputFile)
        return 1

    # Add a 'MAG_OBS' column and a 'MAG_DIFF' column to the pandas DataFrame...
    dataFrame['MAG_OBS'] = -2.5*np.log10(dataFrame[fluxObsColName])
    dataFrame['MAG_DIFF'] = dataFrame[magStdColName]-dataFrame['MAG_OBS']


    ###############################################
    # Aggregate by aggFieldColName
    ###############################################

    # Make a copy of original dataFrame...
    df = dataFrame.copy()

    # Create an initial mask...
    mask1 = ( (df[magStdColName] >= 0.) & (df[magStdColName] <= 25.) )
    if magerrStdColName != 'None':  
        mask1 =( mask1 & (df[magerrStdColName] < 0.1) )
    magPsfDiffGlobalMedian = df[mask1]['MAG_DIFF'].median()
    magPsfDiffMin = magPsfDiffGlobalMedian - 5.0
    magPsfDiffMax = magPsfDiffGlobalMedian + 5.0
    mask2 = ( (df['MAG_DIFF'] > magPsfDiffMin) & (df['MAG_DIFF'] < magPsfDiffMax) )
    mask = mask1 & mask2

    # Iterate over the copy of dataFrame 3 times, removing outliers...
    #  We are using "Method 2/Group by item" from
    #  http://nbviewer.jupyter.org/urls/bitbucket.org/hrojas/learn-pandas/raw/master/lessons/07%20-%20Lesson.ipynb
    print "Sigma-clipping..."
    niter = 0
    for i in range(3):

        niter = i + 1
        print """   iter%d...""" % ( niter )

        # make a copy of original df, and then delete the old one...
        newdf = df[mask].copy()
        del df

        # group by aggFieldColName...
        grpnewdf = newdf.groupby([aggFieldColName])
        
        # add/update new columns to newdf
        print datetime.datetime.now()
        newdf['Outlier']  = grpnewdf['MAG_DIFF'].transform( lambda x: abs(x-x.mean()) > 3.00*x.std() )
        print datetime.datetime.now()
        del grpnewdf
        print datetime.datetime.now()
        #print newdf

        nrows = newdf['MAG_DIFF'].size
        print """  Number of rows remaining:  %d""" % ( nrows )

        df = newdf
        mask = ( df['Outlier'] == False )  


    # Perform pandas grouping/aggregating functions on sigma-clipped Data Frame...
    print datetime.datetime.now()
    print 'Performing grouping/aggregating functions on sigma-clipped pandas DataFrame...'
    groupedDataFrame = df.groupby([aggFieldColName])
    magPsfZeroMedian = groupedDataFrame['MAG_DIFF'].median()
    magPsfZeroMean = groupedDataFrame['MAG_DIFF'].mean()
    magPsfZeroStd = groupedDataFrame['MAG_DIFF'].std()
    magPsfZeroNum = groupedDataFrame['MAG_DIFF'].count()
    magPsfZeroErr = magPsfZeroStd/np.sqrt(magPsfZeroNum-1)
    print datetime.datetime.now()
    print

    # Rename these pandas series...
    magPsfZeroMedian.name = 'MAG_ZERO_MEDIAN'
    magPsfZeroMean.name = 'MAG_ZERO_MEAN'
    magPsfZeroStd.name = 'MAG_ZERO_STD'
    magPsfZeroNum.name = 'MAG_ZERO_NUM'
    magPsfZeroErr.name = 'MAG_ZERO_MEAN_ERR'

    # Also, calculate group medians for all columns in df that have a numerical data type...
    numericalColList = df.select_dtypes(include=[np.number]).columns.tolist()
    groupedDataMedian = {}
    for numericalCol in numericalColList:
        groupedDataMedian[numericalCol] = groupedDataFrame[numericalCol].median()
        groupedDataMedian[numericalCol].name = """%s_MEDIAN""" % (numericalCol)

    # Create new data frame containing all the relevant aggregate quantities
    #newDataFrame = pd.concat( [magPsfZeroMedian, magPsfZeroMean, magPsfZeroStd, \
    #                           magPsfZeroErr, magPsfZeroNum], \
    #                           join='outer', axis=1 )
    seriesList = []
    for numericalCol in numericalColList:
        seriesList.append(groupedDataMedian[numericalCol])
    seriesList.extend([magPsfZeroMedian, magPsfZeroMean, magPsfZeroStd, \
                               magPsfZeroErr, magPsfZeroNum])
    print seriesList
    newDataFrame = pd.concat( seriesList, join='outer', axis=1 )
                               

    #newDataFrame.index.rename('FILENAME', inplace=True)

    # Saving catname-based results to output files...
    print datetime.datetime.now()
    print """Writing %s output file (using pandas to_csv method)...""" % (outputFile)
    newDataFrame.to_csv(outputFile, float_format='%.4f')
    print datetime.datetime.now()
    print

    return 0


##################################

if __name__ == "__main__":
    main()

##################################

