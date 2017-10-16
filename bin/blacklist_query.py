#!/usr/bin/env python
"""
    blacklist_query.py

    Example:
    
    blacklist_query.py --help

    blacklist_query.py --tag Y4A1_FINALCUT \
                         --outputFile y4a1_rawdata.y4.blacklist.csv \
                         --verbose 2
    
    """

##################################

def main():

    import argparse
    import time

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--outputFile', help='name of an output file', default='output.csv')
    parser.add_argument('--expnumMin', help='Minimum value for EXPNUM', default=1, type=int)
    parser.add_argument('--expnumMax', help='Maximum value for EXPNUM', default=99999999, type=int)
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    # Run exposure query...
    status = blacklist_query(args)


##################################
# 

def blacklist_query(args):

    import easyaccess as ea
    import pandas as pd
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'blacklist_query'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    outputFile=args.outputFile    

    query="""
        select * 
        from prod.blacklist 
        where expnum between %d and %d
        order by expnum, ccdnum
        """ % (args.expnumMin, args.expnumMax)
    
    if args.verbose>0: print query
    
    # Open connection to db...
    connection=ea.connect('desoper') 

    # Make query and save to output file...
    connection.query_and_save(query,outputFile)
    
    # Close connection to db...
    connection.close()
    
    if args.verbose>0: print

    return 0


##################################

if __name__ == "__main__":
    main()

##################################
