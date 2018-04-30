#!/usr/bin/env python
"""
    rawdata_exp_firstcut_query.py

    Example:
    
    rawdata_exp_firstcut_query.py --help

    rawdata_exp_firstcut_query.py --tag Y5N_FIRSTCUT \
                                  --outputFile y5n_rawdata.expinfo.csv \
                                  --verbose 2
    
    """

##################################

def main():

    import argparse
    import time

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--outputFile', help='name of an output file', default='output.csv')
    parser.add_argument('--tag', help='PROCTAG tag name', default='Y5N_FIRSTCUT')
    parser.add_argument('--expnumMin', help='Minimum value for EXPNUM', default=1, type=int)
    parser.add_argument('--expnumMax', help='Maximum value for EXPNUM', default=99999999, type=int)
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    # Run exposure query...
    status = rawdata_exp_firstcut_query(args)


##################################
# 

def rawdata_exp_firstcut_query(args):

    import easyaccess as ea
    import pandas as pd
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'rawdata_exp_firstcut_query'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    query="""
    with x as (
          select /*+materialize */ pfw_attempt_id, tag
          from prod.proctag
          where tag='%s'
          ),
    y as (
          select /*+materialize */ a.expnum as expnum,
                                   max(a.lastchanged_time) as evaltime
          from prod.firstcut_eval a
          where a.analyst!='SNQUALITY'
          group by a.expnum
          ), 
    z as (
          select /*+materialize */ b.*
          from prod.firstcut_eval b, y
          where b.expnum=y.expnum and 
                b.lastchanged_time=y.evaltime and 
                b.expnum between %d and %d
          )
    select e.expnum, 
           e.radeg as EXPRA, e.decdeg as EXPDEC, e.exptime, e.airmass, 
           e.band, e.nite, e.mjd_obs, e.field, e.object, e.program, 
           x.pfw_attempt_id, x.tag,
           z.analyst as SOURCE, z.t_eff, z.fwhm_asec as PSF_FWHM, z.f_eff,
           z.skybrightness, z.b_eff, z.cloud_apass, z.cloud_nomad,
           z.n_apass, z.n_nomad, z.c_eff, z.cloud_des, z.n_des, z.accepted
    from prod.exposure e, x, z
    where e.expnum=z.expnum and
          x.pfw_attempt_id=z.pfw_attempt_id and
         (z.accepted='True' or e.program='photom-std-field') 
    order by e.expnum
    """ % \
        (args.tag, args.expnumMin, args.expnumMax)

    
    if args.verbose>0: print query
    
    outputFile=args.outputFile
    
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
