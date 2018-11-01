#!/usr/bin/env python
"""
    rawdata_se_object_query.py

    Example:
    
    rawdata_se_object_query.py --help

    rawdata_se_object_query.py --tag Y4A1_FINALCUT 
                               --outputFiles y4a1_rawdata.Y4A1_FINALCUT.objinfo.fits 
                               --verbose 2
    """

##################################

def main():

    import argparse
    import time

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--outputFiles', help='name of output files', default='output.csv')
    parser.add_argument('--tag', help='PROCTAG tag name', default='Y2A1_FINALCUT')
    parser.add_argument('--reqnum', help='REQNUM (used by se_object_query_by_reqnum)', default=2285, type=int)
    parser.add_argument('--expnumMin', help='Minimum value for EXPNUM', default=1, type=int)
    parser.add_argument('--expnumMax', help='Maximum value for EXPNUM', default=99999999, type=int)
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    status = rawdata_se_object_query(args)


##################################
# 
# Based on query from Eli Rykoff and Dave Burke, which 
#  is in turn is based on a query from Robert Gruendl.

def rawdata_se_object_query(args):

    import easyaccess as ea 
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'rawdata_se_object_query'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    connection=ea.connect('desoper') 

    # Based on query from Eli Rykoff (9 Oct 2017), 
    #  and input from days of old from Robert Gruendl...
    #query="""
    #with x as (
    #    select /*+materialize */ pfw_attempt_id
    # 	from prod.proctag
	#where tag='%s'
	#),
    #y as (
    #    select /*+materialize */ f.expnum, f.pfw_attempt_id
	#from prod.finalcut_eval f, x
	#where f.pfw_attempt_id = x.pfw_attempt_id
	#),
    #z as (
    #    select /*+materialize */ e.expnum,
    #        cast(e.band as VARCHAR(1)) as BAND,
    #        c.filename,c.ccdnum
    # 	from prod.catalog c, prod.exposure e, y
	#where c.expnum=y.expnum
	#    and e.expnum=c.expnum
	#    and c.pfw_attempt_id=y.pfw_attempt_id
	#    and c.filetype='cat_finalcut'
        #    and e.expnum between %d and %d
    	#)
    #select cast(z.expnum as NUMBER(8)) as EXPNUM,
    #       z.filename, z.band, 
    #       cast(z.ccdnum as NUMBER(4)) as CCDNUM,
    #       o.object_number, o.x_image, o.y_image, o.ra, o.dec, 
    #       o.flux_psf, o.fluxerr_psf, o.flux_aper_8, fluxerr_aper_8,
    #       o.class_star, o.spread_model, o.spreaderr_model, o.flags
    #from prod.se_object o, z
    #where o.filename=z.filename
    #    and (o.x_image between 100 and 1900)
    #    and (o.y_image between 100 and 3950)
    #    and o.flags = 0
    #    and o.imaflags_iso = 0
    #    and (o.spread_model between -0.01 and 0.01)
    #    and o.class_star > 0.80
    #    and ((1.086*o.fluxerr_psf/o.flux_psf) between 0.0001 and 0.1000)
    #    """ % \
    #    (args.tag, args.expnumMin, args.expnumMax)

    # Based on query from Eli Rykoff (9 Oct 2017), 
    #  and input from days of old from Robert Gruendl...
    # 20180629:  changed "o.flags = 0" and "o.imaflags_iso = 0"
    #            to "o.flags < 2" and "o.imaflags_iso < 2"

    query="""
    with x as (
        select /*+materialize */ pfw_attempt_id
     	from prod.proctag
	where tag='%s'
	),
    y as (
        select /*+materialize */ qa.expnum, qa.pfw_attempt_id,
            qa.flag,qa.calnac
	from prod.qa_summary qa, x
	where qa.pfw_attempt_id = x.pfw_attempt_id
	    and qa.flag < 20000
	    and qa.calnac < 20000
	),
    z as (
        select /*+materialize */ e.expnum,
            cast(e.band as VARCHAR(1)) as BAND,
            c.filename,c.ccdnum
     	from prod.catalog c, prod.exposure e, y
	where c.expnum=y.expnum
	    and e.expnum=c.expnum
	    and c.pfw_attempt_id=y.pfw_attempt_id
	    and c.filetype='cat_finalcut'
            and e.expnum between %d and %d
    	)
    select cast(z.expnum as NUMBER(8)) as EXPNUM,
           z.filename, z.band, 
           cast(z.ccdnum as NUMBER(4)) as CCDNUM,
           o.object_number, o.x_image, o.y_image, o.ra, o.dec, 
           o.flux_psf, o.fluxerr_psf, o.flux_aper_8, fluxerr_aper_8,
           o.class_star, o.spread_model, o.spreaderr_model, o.flags
    from prod.se_object o, z
    where o.filename=z.filename
        and (o.x_image between 100 and 1900)
        and (o.y_image between 100 and 3950)
        and o.flags < 2
        and o.imaflags_iso < 2
        and (o.spread_model between -0.01 and 0.01)
        and o.class_star > 0.80
        and ((1.086*o.fluxerr_psf/o.flux_psf) between 0.0001 and 0.1000)
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
