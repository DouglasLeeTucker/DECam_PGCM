#!/usr/bin/env python
"""
    rawdata_exp_query.py

    Example:
    
    rawdata_exp_query.py --help

    rawdata_exp_query.py --tag Y4A1_FINALCUT \
                         --outputFile y4a1_rawdata.expinfo.csv \
                         --verbose 2
    
    """

##################################

def main():

    import argparse
    import time

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--outputFile', help='name of an output file', default='output.csv')
    parser.add_argument('--tag', help='PROCTAG tag name', default='Y4A1_FINALCUT')
    parser.add_argument('--expnumMin', help='Minimum value for EXPNUM', default=1, type=int)
    parser.add_argument('--expnumMax', help='Maximum value for EXPNUM', default=99999999, type=int)
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    # Run exposure query...
    status = rawdata_image_query(args)


##################################
# 


def rawdata_image_query(args):

    import easyaccess as ea 
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'rawdata_image_query'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    connection=ea.connect('desoper') 

    query="""
    with x as (
        select /*+materialize */ t.pfw_attempt_id, qa.expnum 
            from prod.proctag t, prod.qa_summary qa 
            where t.pfw_attempt_id=qa.pfw_attempt_id and 
                  tag='%s' and
                  qa.expnum between %d and %d
        ) 
    select i.pfw_attempt_id, i.band, i.expnum, i.ccdnum, 
           i.rac1, i.decc1, i.rac2, i.decc2, i.rac3, i.decc3, i.rac4, i.decc4, 
           i.ra_cent, i.dec_cent, i.racmax, i.racmin, i.crossra0 
    from prod.image i, x 
    where x.pfw_attempt_id=i.pfw_attempt_id and 
          i.filetype='red_starflat' and 
          i.expnum=x.expnum 
    order by i.expnum, i.ccdnum
    """ % (args.tag, args.expnumMin, args.expnumMax)


    if args.verbose>0: print query

    outputFile = args.outputFile

    connection.query_and_save(query,outputFile)

    connection.close()

    if args.verbose>0: print

    return 0


##################################

if __name__ == "__main__":
    main()

##################################
