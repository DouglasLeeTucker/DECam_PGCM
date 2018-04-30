#!/usr/bin/env python
"""
    rawdata_se_object_firstcut_query.py

    Example:
    
    rawdata_se_object_firstcut_query.py --help

    rawdata_se_object_firstcut_query.py --tag Y5N_FIRSTCUT 
                                        --outputFiles y5n_rawdata.Y5N_FIRSTCUT.objinfo.fits 
                                        --verbose 2
    """

##################################

def main():

    import argparse
    import time

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--outputFiles', help='name of output files', default='output.csv')
    parser.add_argument('--tag', help='PROCTAG tag name', default='Y5N_FIRSTCUT')
    parser.add_argument('--expnumMin', help='Minimum value for EXPNUM', default=1, type=int)
    parser.add_argument('--expnumMax', help='Maximum value for EXPNUM', default=99999999, type=int)
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    status = rawdata_se_object_firstcut_query(args)


##################################
# 
# Based on a query from Robert Gruendl, with some 
#  modifications suggested by Eli Rykoff.

def rawdata_se_object_firstcut_query(args):

    import easyaccess as ea 
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'rawdata_se_object_firstcut_query'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    connection=ea.connect('desoper') 

    # Based on query from Robert Gruendl from early Y3A1 tests,
    #  with some more modern (Y4A1 vintage) suggestions from 
    #  Eli Rykoff...
    # 
    # (Originally, included "b.accepted='True'" in WHERE
    #  statement for z, but realized that this omitted 
    #  lots of standard star field stars.  Now, this 
    #  WHERE clause is included as part of a  
    #  "(z.accepted='True' or e.program='photom-std-field')"
    #  WHERE clause for the exposure query.
    #  Thus, be careful where you want to include standard
    #  star fields or not downstream...)

    query="""
    with x as (
          select /*+materialize */ pfw_attempt_id
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
          select /*+materialize */ b.expnum
          from prod.firstcut_eval b, y
          where b.expnum=y.expnum and
                b.lastchanged_time=y.evaltime and
                b.expnum between %d and %d
          )
    select cast(c.expnum as NUMBER(8)) as EXPNUM,
           c.filename, c.band, 
           cast(c.ccdnum as NUMBER(4)) as CCDNUM, 
           o.object_number, o.x_image, o.y_image, o.ra, o.dec,
           o.flux_psf, o.fluxerr_psf, o.flux_aper_8, fluxerr_aper_8,
           o.class_star, o.spread_model, o.spreaderr_model, o.flags
    from prod.se_object o, prod.catalog c, x, z
    where x.pfw_attempt_id=c.pfw_attempt_id and
          c.filetype='cat_firstcut' and
          c.expnum=z.expnum and
          o.filename=c.filename and
         (o.x_image between 100 and 1900) and
         (o.y_image between 100 and 3950) and
          o.flags = 0 and o.imaflags_iso = 0 and
         (o.spread_model between -0.01 and 0.01) and o.class_star > 0.80 and
        ((1.086*o.fluxerr_psf/o.flux_psf) between 0.0001 and 0.1000)
        """ % \
        (args.tag, args.expnumMin, args.expnumMax)


    if args.verbose>0: print query

    outputFiles=args.outputFiles

    connection.query_and_save(query,outputFiles)

    connection.close()

    if args.verbose>0: print

    return 0


##################################

if __name__ == "__main__":
    main()

##################################
