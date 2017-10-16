#!/usr/bin/env python
"""
    y2a1_rawdata.py

    Example:
    
    y2a1_rawdata.py --help

    y2a1_rawdata.py --tag Y2A1_FINALCUT --baseName y2a1_rawdata --verbose 2
    
    """

##################################

def main():

    import argparse
    import time

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--baseName', help='base name of the output CSV file', default='rawdata')
    parser.add_argument('--fileNamePattern', help='pattern of filenames to be globbed (used by combine_stilts_internal_match_files and reformat_stilts_internal_match_files_for_gcm_exp)', default='y2a1_rawdata.svy1y2y3.sorted.tmp.z.csv.??.inmatch')
    parser.add_argument('--outputFileName', help='name of an output file', default='output.csv')
    parser.add_argument('--tag', help='PROCTAG tag name', default='Y2A1_FINALCUT')
    parser.add_argument('--reqnum', help='REQNUM (used by se_object_query_by_reqnum)', default=2285, type=int)
    parser.add_argument('--expnumMin', help='Minimum value for EXPNUM', default=1, type=int)
    parser.add_argument('--expnumMax', help='Maximum value for EXPNUM', default=99999999, type=int)
    parser.add_argument('--nite', help='nite to query (used by se_object_query_by_nite_and_band)', default='20130113')
    parser.add_argument('--bands', help='comma-separated list of filter bands to consider', default='g,r,i,z,Y')
    parser.add_argument('--band', help='single filter band to consider (used by se_object_query_by_nite_and_band)', default='u')
    parser.add_argument('--cleanup', help='remove less interesting intermediate files', default=False, action='store_true')
    parser.add_argument('--do_exposure_query',help='do exposure query', default=False, action='store_true')
    parser.add_argument('--do_reqnum_query',help='do reqnum query', default=False, action='store_true')
    parser.add_argument('--do_se_object_query_by_reqnum',help='do object query', default=False, action='store_true')
    parser.add_argument('--do_se_object_query_by_expnumRange',help='do object query', default=False, action='store_true')
    parser.add_argument('--do_reqnum_query_results_to_csv_band_files',help='do convert reqnum_query results to a set of CSV files (one per filter band)', default=False, action='store_true')
    parser.add_argument('--do_nite_band_query_results_to_csv_band_files',help='do convert night_band_query results to a set of CSV files (one per filter band)', default=False, action='store_true')
    parser.add_argument('--do_add_exp_info',help='do add exposure info to the se_object files', default=False, action='store_true')
    parser.add_argument('--do_combine_stilts_internal_match_files',help='do combine STILTS internally matched files', default=False, action='store_true')
    parser.add_argument('--do_reformat_stilts_internal_match_files_for_gcm_exp',help='do reformat STILTS internally matched filesfor GCM (exposure-by-exposure fit)', default=False, action='store_true')
    parser.add_argument('--do_blacklist_query',help='do blacklist query', default=False, action='store_true')
    parser.add_argument('--do_y3a1_blacklist_query',help='do y3a1_blacklist query', default=False, action='store_true')
    parser.add_argument('--do_image_query',help='do image query', default=False, action='store_true')
    parser.add_argument('--do_exp_image_cat_y3a1blacklist_query',help='do exposure/image/catalog/y3a1_blacklist query and merge', default=False, action='store_true')
    parser.add_argument('--do_se_object_query_by_nite_and_band',help='do nite/band objects query', default=False, action='store_true')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    # Checking bandList...
    supportedBandList = ['u','g','r','i','z','Y']
    bandList = args.bands.strip().split(',')
    badBandCount = 0
    for band in bandList:
        try:
            supportedBandIndex = supportedBandList.index(band)
        except:
            print """Filter band %s not found in list of supported bands...""" % (band)
            badBandCount = badBandCount + 1
    if badBandCount > 0:
        print """The list of supported bands is %s""" % (supportedBandList)
        print """Note that the filter band names in this list are case-sensitive."""
        print """Exiting now!..."""
        print 
        return 1


    # Run exposure query...
    if args.do_exposure_query==True:
        status = exposure_query(args)

    # Run object query by reqnum...
    if args.do_se_object_query_by_reqnum==True:
        status = se_object_query_by_reqnum(args)

    # Run object query by expnum range...
    if args.do_se_object_query_by_expnumRange==True:
        status = se_object_query_by_expnumRange(args)

    # Run object query by nite and band...
    if args.do_se_object_query_by_nite_and_band==True:
        status = se_object_query_by_nite_and_band(args)

    # Run query to find unique list of reqnums for
    #  a given processing tag and within a given expnum range...
    if args.do_reqnum_query==True:
        reqnum_query(args)

    # Convert results from reqnum_query to a set of 
    #  CSV files (one per filter band)...
    if args.do_reqnum_query_results_to_csv_band_files==True:
        reqnum_query_results_to_csv_band_files(args)

    # Convert results from nite_band_query to a set of 
    #  CSV files (one per filter band)...
    if args.do_nite_band_query_results_to_csv_band_files==True:
        nite_band_query_results_to_csv_band_files(args)

    # Add exposure info to the se_object files...
    if args.do_add_exp_info==True:
        add_exp_info(args)

    # Combine files from STILTS internal match into a
    #  single large file...
    if args.do_combine_stilts_internal_match_files==True:
        combine_stilts_internal_match_files(args)

    # Reformat files from STILTS internal match into a
    #  single large file that is ingestible by GCM for
    #  the exposure-by-exposure fit...
    if args.do_reformat_stilts_internal_match_files_for_gcm_exp==True:
        reformat_stilts_internal_match_files_for_gcm_exp(args)

    # Run blacklist query...
    if args.do_blacklist_query==True:
        status = blacklist_query(args)

    # Run y3a1_blacklist query...
    #  (this includes Eli Rykoff's spread_model blacklist just for y3a1...)
    if args.do_y3a1_blacklist_query==True:
        status = y3a1_blacklist_query(args)

    # Run image query...
    if args.do_image_query==True:
        status = image_query(args)

    if args.do_exp_image_cat_y3a1blacklist_query==True:
        status = exp_image_cat_y3a1blacklist_query(args)


##################################
# 

def exposure_query(args):

    import easyaccess as ea 
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'exposure_query'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    print 'TBD...'

    connection=ea.connect('desoper') 
    
    query="""
        select e.expnum, t.unitname, t.reqnum, t.attnum, 
               e.radeg as EXPRA, e.decdeg as EXPDEC, e.exptime, e.airmass, e.band, 
               e.nite, e.mjd_obs, e.field, e.object, e.program, 
               qa.*
        from prod.exposure e, prod.proctag t, prod.qa_summary qa
        where t.pfw_attempt_id=qa.pfw_attempt_id and
              e.expnum=qa.expnum and 
              t.tag='%s' and
              qa.expnum between %d and %d
        order by qa.expnum
        """ % (args.tag, args.expnumMin, args.expnumMax)
    
    if args.verbose>0: print query
    
    outputFile="""%s.expinfo.csv""" % (args.baseName)
    
    connection.query_and_save(query,outputFile)
    
    connection.close()
    
    if args.verbose>0: print

    return 0


##################################
#
# Based on query from Eli Rykoff and Dave Burke, which 
#  is in turn is based on a query from Robert Gruendl.

def se_object_query_orig(args):

    import easyaccess as ea 
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'y2a1_se_object_query'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    connection=ea.connect('desoper') 

    query="""
    with x as (
        select /*+materialize */ pfw_attempt_id
     	    from prod.proctag
	    where tag='Y2A1_FINALCUT'
	),
    y as (
        select /*+materialize */ a.expnum as expnum, max(a.lastchanged_time) as evaltime
source     	    from prod.finalcut_eval a
	    where a.analyst!='SNQUALITY'
	    and (upper(a.program)='SURVEY' or upper(a.program)='SUPERNOVA')
	    group by a.expnum
        ),
    z as (
        select /*+materialize */ b.expnum
     	    from prod.finalcut_eval b, y
	    where b.expnum=y.expnum
	    and b.lastchanged_time=y.evaltime
	    and b.accepted='True'
	    )
    select o.ra, o.dec, o.flux_psf, o.fluxerr_psf, o.class_star,
       o.spread_model, cast(o.band as VARCHAR(1)) as band,
       cast(c.expnum as NUMBER(8)) as expnum,
       cast(c.ccdnum as NUMBER(4)) as ccdnum, e.mjd_obs, e.exptime,
       e.tradeg, e.tdecdeg, cast(e.ha as varchar(12)) as ha
    from prod.se_object o, prod.catalog c, prod.exposure e, x, z
    where x.pfw_attempt_id=c.pfw_attempt_id
        and c.filetype='cat_finalcut'
        and c.expnum=z.expnum
        and e.expnum=c.expnum
        and o.filename=c.filename
        and (o.x_image between 100 and 1900)
        and (o.y_image between 100 and 3950)
        and o.flags = 0
        and o.class_star > 0.50
        and o.spread_model < 0.01
        and ((1.086*o.fluxerr_psf/o.flux_psf) between 0.001 and 0.100)
        and (c.expnum between 383995 and 410586)""" 
     
    if args.verbose>0: print query

    outputFiles="""%s.objinfo.fits""" % (args.baseName)

    connection.query_and_save(query,outputFiles)

    connection.close()

    if args.verbose>0: print

    return 0


##################################
#
# Based on query from Eli Rykoff and Dave Burke, which 
#  is in turn is based on a query from Robert Gruendl.

def se_object_query_by_reqnum_old(args):

    import easyaccess as ea 
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'y2a1_se_object_query_by_reqnum_old'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    connection=ea.connect('desoper') 

    query="""
    with x as (
        select /*+materialize */ pfw_attempt_id
     	    from prod.proctag
	    where tag='%s' and reqnum=%d
	),
    y as (
        select /*+materialize */ a.expnum as expnum, max(a.lastchanged_time) as evaltime
     	    from prod.finalcut_eval a
	    where a.analyst!='SNQUALITY'
	    group by a.expnum
        ),
    z as (
        select /*+materialize */ b.expnum
     	    from prod.finalcut_eval b, y
	    where b.expnum=y.expnum
	    and b.lastchanged_time=y.evaltime
	    and b.accepted='True'
	    )
    select c.expnum, c.filename,
           cast(o.band as VARCHAR(1)) as BAND,
           c.ccdnum as CCDNUM,
           o.object_number, o.x_image, o.y_image, o.ra, o.dec, 
           o.flux_psf, o.fluxerr_psf, o.flux_aper_8, fluxerr_aper_8, 
           o.class_star, o.spread_model, o.spreaderr_model, o.flags
    from prod.se_object o, prod.catalog c, x, z
    where x.pfw_attempt_id=c.pfw_attempt_id
        and c.filetype='cat_finalcut'
        and c.expnum=z.expnum
        and o.filename=c.filename
        and o.flags = 0
        and o.class_star > 0.80
        and o.spread_model < 0.01
        and ((1.086*o.fluxerr_psf/o.flux_psf) between 0.001 and 0.01)
        and (c.expnum between %d and %d)""" % \
        (args.tag, args.reqnum, args.expnumMin, args.expnumMax)
     
    if args.verbose>0: print query

    outputFiles="""%s.%s.%d-%d.%d.objinfo.fits""" % \
        (args.baseName, args.tag, args.expnumMin, args.expnumMax, args.reqnum)

    connection.query_and_save(query,outputFiles)

    connection.close()

    if args.verbose>0: print

    return 0


##################################
#
# Based on query from Eli Rykoff and Dave Burke, which 
#  is in turn is based on a query from Robert Gruendl.

def se_object_query_by_reqnum(args):

    import easyaccess as ea 
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'y2a1_se_object_query_by_reqnum'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    connection=ea.connect('desoper') 

    query="""
    with x as (
        select /*+materialize */ t.pfw_attempt_id, qa.expnum
     	    from prod.proctag t, prod.qa_summary qa
	    where t.pfw_attempt_id=qa.pfw_attempt_id and 
                  tag='%s' and reqnum=%d and 
                  qa.expnum between %d and %d
        )
    select c.expnum, c.filename,
           cast(o.band as VARCHAR(1)) as BAND,
           c.ccdnum as CCDNUM,
           o.object_number, o.x_image, o.y_image, o.ra, o.dec, 
           o.flux_psf, o.fluxerr_psf, o.flux_aper_8, fluxerr_aper_8, 
           o.class_star, o.spread_model, o.spreaderr_model, o.flags
    from prod.se_object o, prod.catalog c, x
    where x.pfw_attempt_id=c.pfw_attempt_id
        and c.filetype='cat_finalcut'
        and c.expnum=x.expnum
        and o.filename=c.filename
        and o.flags = 0
        and o.class_star > 0.80
        and o.spread_model < 0.01
        and ((1.086*o.fluxerr_psf/o.flux_psf) between 0.001 and 0.01)
        """ % \
        (args.tag, args.reqnum, args.expnumMin, args.expnumMax)
     
    if args.verbose>0: print query

    outputFiles="""%s.%s.%d-%d.%d.objinfo.fits""" % \
        (args.baseName, args.tag, args.expnumMin, args.expnumMax, args.reqnum)

    connection.query_and_save(query,outputFiles)

    connection.close()

    if args.verbose>0: print

    return 0


##################################
#
# Based on query from Eli Rykoff and Dave Burke, which 
#  is in turn is based on a query from Robert Gruendl.

def se_object_query_by_nite_and_band(args):

    import easyaccess as ea 
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'y2a1_se_object_query_by_nite_and_band'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    connection=ea.connect('desoper') 


    # Note looser cuts for o.flags and "magerr_psf"...
    query="""
    with x as (
        select /*+materialize */ t.pfw_attempt_id, qa.expnum, e.band 
     	    from prod.proctag t, prod.qa_summary qa, prod.exposure e 
	    where t.pfw_attempt_id=qa.pfw_attempt_id and 
                  qa.expnum=e.expnum and 
                  t.tag='%s' and 
                  e.nite='%s' and 
                  e.band='%s'
        )
    select c.expnum, c.filename,
           x.band,
           c.ccdnum as CCDNUM,
           o.object_number, o.x_image, o.y_image, o.ra, o.dec, 
           o.flux_psf, o.fluxerr_psf, o.flux_aper_8, o.fluxerr_aper_8, 
           o.class_star, o.spread_model, o.spreaderr_model, o.flags
    from prod.se_object o, prod.catalog c, x
    where x.pfw_attempt_id=c.pfw_attempt_id
        and c.filetype='cat_finalcut'
        and c.expnum=x.expnum
        and o.filename=c.filename
        and o.flags < 3 
        and o.class_star > 0.80
        and o.spread_model < 0.01
        and ((1.086*o.fluxerr_psf/o.flux_psf) between 0.001 and 0.03)
        """ % \
        (args.tag, args.nite, args.band)
     
    if args.verbose>0: print query

    outputFiles="""%s.%s.%s.%s.objinfo.fits""" % \
        (args.baseName, args.tag, args.nite, args.band)

    connection.query_and_save(query,outputFiles)

    connection.close()

    if args.verbose>0: print

    return 0


##################################
#
# Based on query from Eli Rykoff and Dave Burke, which 
#  is in turn is based on a query from Robert Gruendl.

def se_object_query_by_expnumRange_old(args):

    import easyaccess as ea 
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'y2a1_se_object_query'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    connection=ea.connect('desoper') 

    query="""
    with x as (
        select /*+materialize */ pfw_attempt_id
     	    from prod.proctag
	    where tag='%s'
	),
    y as (
        select /*+materialize */ a.expnum as expnum, max(a.lastchanged_time) as evaltime
     	    from prod.finalcut_eval a
	    where a.analyst!='SNQUALITY'
	    group by a.expnum
        ),
    z as (
        select /*+materialize */ b.expnum
     	    from prod.finalcut_eval b, y
	    where b.expnum=y.expnum
	    and b.lastchanged_time=y.evaltime
	    and b.accepted='True'
	    )
    select c.expnum, c.filename,
           cast(o.band as VARCHAR(1)) as BAND,
           c.ccdnum as CCDNUM,
           o.object_number, o.x_image, o.y_image, o.ra, o.dec, 
           o.flux_psf, o.fluxerr_psf, o.flux_aper_8, fluxerr_aper_8, 
           o.class_star, o.spread_model, o.spreaderr_model, o.flags
    from prod.se_object o, prod.catalog c, x, z
    where x.pfw_attempt_id=c.pfw_attempt_id
        and c.filetype='cat_finalcut'
        and c.expnum=z.expnum
        and o.filename=c.filename
        and o.flags = 0
        and o.class_star > 0.80
        and o.spread_model < 0.01
        and ((1.086*o.fluxerr_psf/o.flux_psf) between 0.001 and 0.01)
        and (c.expnum between %d and %d)""" % \
        (args.tag, args.expnumMin, args.expnumMax)
     
    if args.verbose>0: print query

    outputFiles="""%s.%s.exp%d-%d.objinfo.fits""" % \
        (args.baseName, args.tag, args.expnumMin, args.expnumMax)

    connection.query_and_save(query,outputFiles)

    connection.close()

    if args.verbose>0: print

    return 0


##################################
##################################
#
# Based on query from Eli Rykoff and Dave Burke, which 
#  is in turn is based on a query from Robert Gruendl.

def se_object_query_by_expnumRange(args):

    import easyaccess as ea 
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'y2a1_se_object_query'
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
    select c.expnum, c.filename,
           cast(o.band as VARCHAR(1)) as BAND,
           c.ccdnum as CCDNUM,
           o.object_number, o.x_image, o.y_image, o.ra, o.dec, 
           o.flux_psf, o.fluxerr_psf, o.flux_aper_8, fluxerr_aper_8, 
           o.class_star, o.spread_model, o.spreaderr_model, o.flags
    from prod.se_object o, prod.catalog c, x
    where x.pfw_attempt_id=c.pfw_attempt_id
        and c.filetype='cat_finalcut'
        and c.expnum=x.expnum
        and o.filename=c.filename
        and o.flags = 0
        and o.class_star > 0.80
        and o.spread_model < 0.01
        and ((1.086*o.fluxerr_psf/o.flux_psf) between 0.001 and 0.01)
        """ % \
        (args.tag, args.expnumMin, args.expnumMax)
     
    if args.verbose>0: print query

    outputFiles="""%s.%s.exp%d-%d.objinfo.fits""" % \
        (args.baseName, args.tag, args.expnumMin, args.expnumMax)

    connection.query_and_save(query,outputFiles)

    connection.close()

    if args.verbose>0: print

    return 0


##################################
#
# Based on query from Robert Gruendl.

def reqnum_query_old(args):

    import easyaccess as ea 
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'reqnum_query_old'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    connection=ea.connect('desoper') 
    
    query="""
        with x as (
           select /*+materialize */ pfw_attempt_id, unitname, reqnum, attnum
           from prod.proctag
           where tag='%s'
        ),
        y as (
           select /*+materialize */ a.expnum as expnum, max(a.lastchanged_time) as evaltime
           from prod.finalcut_eval a
           where a.analyst!='SNQUALITY'
           group by a.expnum
        ),
        z as (
           select /*+materialize */ b.expnum, b.unitname, b.reqnum, b.attnum
           from prod.finalcut_eval b, y
           where b.expnum=y.expnum and 
                 b.lastchanged_time=y.evaltime and 
                 b.accepted='True'
        )
        select distinct x.reqnum
        from x, z
        where z.unitname=x.unitname and 
              z.reqnum=x.reqnum and
              z.attnum=x.attnum and 
              z.expnum between %d and %d
        order by x.reqnum
        """ % (args.tag, args.expnumMin, args.expnumMax)
    

    if args.verbose>0: print query
    
    outputFile="""%s.reqnums.csv""" % (args.baseName)
    
    connection.query_and_save(query,outputFile)
    
    connection.close()
    
    if args.verbose>0: print

    return 0


##################################
#
# Based on query from Robert Gruendl.

def reqnum_query(args):

    import easyaccess as ea 
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'reqnum_query'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    connection=ea.connect('desoper') 
    
    query="""
        select distinct (t.reqnum)  
        from prod.proctag t,  prod.qa_summary qa
        where t.pfw_attempt_id=qa.pfw_attempt_id and 
              t.tag='%s' and
              qa.expnum between %d and %d
        order by t.reqnum
        """ % (args.tag, args.expnumMin, args.expnumMax)
    

    if args.verbose>0: print query
    
    outputFile="""%s.reqnums.csv""" % (args.baseName)
    
    connection.query_and_save(query,outputFile)
    
    connection.close()
    
    if args.verbose>0: print

    return 0


##################################
# reqnum_query_results_to_csv_band_files_old
#
# Based on sepBands from y2a1_tertiaries.py
#

def reqnum_query_results_to_csv_band_files_old(args):

    import os
    import sys
    import glob
    from astropy.io import fits

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'reqnum_query_results_to_csv_band_files_old'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    # This pattern matching will need to be made more general
    #  (maybe instead pass a list of files instead of trying to 
    #   pattern match?)...
    pattern = """%s.%s.??????-??????.????.objinfo.fits""" % \
        (args.baseName, args.tag)
    fileNameList = sorted(glob.glob(pattern))

    if args.verbose>1:  print fileNameList

    # Open up a CSV output file for each filter band...
    bandList = args.bands.strip().split(',')
    fout = {}
    for band in bandList:
        outputFile = """%s.objinfo.%s.csv""" % (args.baseName, band)
        fout[band] = open(outputFile,'w')

    # Read header from the "zeroth" input file from fileNameList
    hdulist = fits.open(fileNameList[0])
    columnNamesList = hdulist[1].columns.names
    hdulist.close()    

    # Check if there is column called 'BAND'...
    try:
        bandCol = columnNamesList.index('BAND')
    except:
        print """Could not find 'BAND' in header of %s...""" % (fileNameList[0])
        print """Exiting 'reqnum_query_results_to_csv_band_files' method with return code 1 now..."""
        return 1

    # Create a CSV header line from columnNamesList...
    hout = ''
    for columnName in columnNamesList:
        hout = hout+columnName+','
    hout = hout[:-1]
    hout = hout+'\n'

    # Write header to each of the output files (one per filter band)...    
    for band in bandList:
        fout[band].write(hout)

    # Loop through each of the input files from fileNameList...
    invalidBandList = []
    for fileName in fileNameList:

        print fileName

        # Open the input file...
        hdulist = fits.open(fileName)
        tbdata = hdulist[1].data

        for linecnt in range(tbdata.size):

            if ( (linecnt/1000.0 == int(linecnt/1000.0)) and (args.verbose > 1) ):
                print '\r'+'Progress (lines read from '+fileName+'):  ',linecnt,
                sys.stdout.flush()
                    
            lins = tbdata[linecnt]
            #band = tbdata['BAND'][linecnt]
            band = lins[bandCol]

            if band in bandList:

                # Create a CSV data line from lins...
                lin = ''
                for colcnt in range(len(lins)):
                    lin = lin+str(lins[colcnt])+','
                lin = lin[:-1]
                lin = lin+'\n'

                # Output CSV data line to appropriate output file...
                fout[band].write(lin)

            elif band not in invalidBandList: 

                invalidBandList.append(band)
                print 'Unrecognized band:  '+band


        # Close this input file...
        hdulist.close()

        if args.verbose > 1:  
            print


    # Close the output file for each filter band...
    for band in bandList:
        fout[band].close
                          
    if args.verbose>0: print

    return 0


##################################
# reqnum_query_results_to_csv_band_files
#
# Based on sepBands from y2a1_tertiaries.py
# Uses STILTS for quick conversion of FITS to CSV of input files.
#

def reqnum_query_results_to_csv_band_files(args):

    import os
    import sys
    import datetime
    import glob
    from astropy.io import fits

    stiltsDir='/usrdevel/dp0/dtucker/STILTS3.0/latest'

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'reqnum_query_results_to_csv_band_files'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    # This pattern matching will need to be made more general
    #  (maybe instead pass a list of files instead of trying to 
    #   pattern match?)...
    #pattern = """%s.%s.??????-??????.????.objinfo.fits""" % \
    pattern = """%s.%s.??????-??????.????.objinfo*.fits""" % \
        (args.baseName, args.tag)
    fileNameList = sorted(glob.glob(pattern))

    if args.verbose>1:  print fileNameList

    # Convert each FITS file to CSV format via the STILTS tcopy command...
    print 'Converting files from FITS to CSV format...'
    for fileName in fileNameList:
        print fileName
        inputFile = fileName
        outputFile = """%s.csv""" % (inputFile)
        cmd = """%s/stilts -Xmx8000m tcopy in=%s ifmt=fits out=%s ofmt=csv""" % (stiltsDir,inputFile,outputFile)
        print datetime.datetime.now()
        print """Running:  %s""" % (cmd)
        status = os.system(cmd)
        print datetime.datetime.now()
        print 
        if status !=0:
            print """%s failed.""" % (cmd)
            return 1

    # Open up a CSV output file for each filter band...
    bandList = args.bands.strip().split(',')
    fout = {}
    for band in bandList:
        outputFile = """%s.objinfo.%s.csv""" % (args.baseName, band)
        fout[band] = open(outputFile,'w')

    # Read header from the "zeroth" input file from fileNameList
    #  (well, actually from its equivalent CSV file)...
    csvFileName = """%s.csv""" % (fileNameList[0])
    fin = open(csvFileName)
    hin = fin.readline()
    fin.close()

    # Identify which column contains the filter band...
    hins=hin.upper().strip().split(',')
    try:
        bandCol = hins.index('BAND')
    except:
        print """Could not find 'BAND' in header of %s...""" % (fileNameList[0])
        print """Exiting 'reqnum_query_results_to_csv_band_files' method with return code 1 now..."""
        return 1

    # Write header to each of the output files (one per filter band)...
    for band in bandList:
        fout[band].write(hin)

    # Loop through each of the input files from fileNameList...
    invalidBandList = []
    for fileName in fileNameList:

        csvFileName = """%s.csv""" % (fileName)
        print csvFileName

        # Open the input file...
        with open(csvFileName,'r') as fin:

            # Skip header line...
            try:
                next(fin)
            except:
                print """%s is empty.  Skipping...""" % (csvFileName)
                continue
            
            print 
            linecnt = 0
            for lin in fin:
            
                # Increment line count...
                linecnt += 1
                if ( (linecnt/1000.0 == int(linecnt/1000.0)) and (args.verbose > 1) ):
                    print '\r'+'Progress (lines read from '+csvFileName+'):  ',linecnt,
                    sys.stdout.flush()
                    
                lins=lin.strip().split(',')
                band = lins[bandCol]

                if band in bandList:
                    fout[band].write(lin)
                else:
                    print 'Unrecognized band:  '+band

        if args.verbose > 1:  
            print


    # Close the output file for each filter band...
    for band in bandList:
        fout[band].close

    # If the "--cleanup" option was used, remove the csv versions
    #  of the original files...
    if args.cleanup:
        for fileName in fileNameList:
            csvFileName = """%s.csv""" % (fileName)
            os.remove(csvFileName)

    if args.verbose>0: print

    return 0

##################################
# nite_band_query_results_to_csv_band_files
#
# Almost identical to reqnum_quer_results_to_csv_bands_files above.
# Uses STILTS for quick conversion of FITS to CSV of input files.
#

def nite_band_query_results_to_csv_band_files(args):

    import os
    import sys
    import datetime
    import glob
    from astropy.io import fits

    stiltsDir='/usrdevel/dp0/dtucker/STILTS3.0/latest'

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'nite_band_query_results_to_csv_band_files'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    # This pattern matching will need to be made more general
    #  (maybe instead pass a list of files instead of trying to 
    #   pattern match?)...
    #pattern = """%s.%s.??????-??????.????.objinfo.fits""" % \
    #pattern = """%s.%s.??????-??????.????.objinfo*.fits""" % \
    #    (args.baseName, args.tag)

    # The different pattern is basically the only difference
    #  between nite_band_query_results_to_csv_band_files and
    #  reqnum_query_results_to_csv_band_files...
    pattern = """%s.%s.????????.?.objinfo*.fits""" % \
        (args.baseName, args.tag)
    fileNameList = sorted(glob.glob(pattern))

    if args.verbose>1:  print fileNameList

    # Convert each FITS file to CSV format via the STILTS tcopy command...
    print 'Converting files from FITS to CSV format...'
    for fileName in fileNameList:
        print fileName
        inputFile = fileName
        outputFile = """%s.csv""" % (inputFile)
        cmd = """%s/stilts -Xmx8000m tcopy in=%s ifmt=fits out=%s ofmt=csv""" % (stiltsDir,inputFile,outputFile)
        print datetime.datetime.now()
        print """Running:  %s""" % (cmd)
        status = os.system(cmd)
        print datetime.datetime.now()
        print 
        if status !=0:
            print """%s failed.""" % (cmd)
            return 1

    # Open up a CSV output file for each filter band...
    bandList = args.bands.strip().split(',')
    fout = {}
    for band in bandList:
        outputFile = """%s.objinfo.%s.csv""" % (args.baseName, band)
        fout[band] = open(outputFile,'w')

    # Read header from the "zeroth" input file from fileNameList
    #  (well, actually from its equivalent CSV file)...
    csvFileName = """%s.csv""" % (fileNameList[0])
    fin = open(csvFileName)
    hin = fin.readline()
    fin.close()

    # Identify which column contains the filter band...
    hins=hin.upper().strip().split(',')
    try:
        bandCol = hins.index('BAND')
    except:
        print """Could not find 'BAND' in header of %s...""" % (fileNameList[0])
        print """Exiting 'reqnum_query_results_to_csv_band_files' method with return code 1 now..."""
        return 1

    # Write header to each of the output files (one per filter band)...
    for band in bandList:
        fout[band].write(hin)

    # Loop through each of the input files from fileNameList...
    invalidBandList = []
    for fileName in fileNameList:

        csvFileName = """%s.csv""" % (fileName)
        print csvFileName

        # Open the input file...
        with open(csvFileName,'r') as fin:

            # Skip header line...
            try:
                next(fin)
            except:
                print """%s is empty.  Skipping...""" % (csvFileName)
                continue
            
            print 
            linecnt = 0
            for lin in fin:
            
                # Increment line count...
                linecnt += 1
                if ( (linecnt/1000.0 == int(linecnt/1000.0)) and (args.verbose > 1) ):
                    print '\r'+'Progress (lines read from '+csvFileName+'):  ',linecnt,
                    sys.stdout.flush()
                    
                lins=lin.strip().split(',')
                band = lins[bandCol]

                if band in bandList:
                    fout[band].write(lin)
                else:
                    print 'Unrecognized band:  '+band

        if args.verbose > 1:  
            print


    # Close the output file for each filter band...
    for band in bandList:
        fout[band].close

    # If the "--cleanup" option was used, remove the csv versions
    #  of the original files...
    if args.cleanup:
        for fileName in fileNameList:
            csvFileName = """%s.csv""" % (fileName)
            os.remove(csvFileName)

    if args.verbose>0: print

    return 0

##################################
# add_exp_info
#

def add_exp_info(args):

    import numpy as np 
    import os
    import sys

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'add_exp_info'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    # Read in expInfoFile...
    expInfoFile="""%s.expinfo.csv""" % (args.baseName)
    if os.path.isfile(expInfoFile)==False:
        print """Exposure info file %s does not exist...""" % (expInfoFile)
        print """Exiting addExpInfo method with return code 1"""
        return 1

    expInfoDict = {}
    with open(expInfoFile,'r') as fin:

        # Read header line...
        try:
            hin = fin.readline()
        except:
            print """%s is empty."""  % (expInfoFile)
            print """Exiting 'add_exp_info' method with return code 1 now..."""
            return 1

        # Identify column in the expInfoFile CSV file corresponding to EXPNUM...
        hins=hin.upper().strip().split(',')
        try:
            expnumCol = hins.index('EXPNUM')
        except:
            print """Could not find 'EXPNUM' in header of %s...""" % (expInfoFile)
            print """Exiting 'add_exp_info' method with return code 1 now..."""
            return 1

        hin_orig = hin

        print
        linecnt = 0
        for lin in fin:
            
            # Increment line count...
            linecnt += 1
            if ( (linecnt/1000.0 == int(linecnt/1000.0)) and (args.verbose > 1) ):
                print '\r'+'Progress (lines read from '+expInfoFile+'):  ',linecnt,
                sys.stdout.flush()
                    
            lins=lin.strip().split(',')
            expnum = int(lins[expnumCol])
            expnumString = str(expnum)
            expInfoDict[expnumString] = lin.strip()

    if args.verbose > 1:  
        print

    # Loop through each filter band in the official bandList...
    bandList = args.bands.strip().split(',')
    for band in bandList:

        # Check for existence of objInfoFile for this filter band...
        objInfoFile = """%s.objinfo.%s.csv""" % (args.baseName, band)
        if os.path.isfile(objInfoFile)==False:
            print """%s does not exist... Skipping""" % (objInfoFile)
            continue

        # Open the outputFile for this filter band...
        outputFile = """%s.%s.csv""" % (args.baseName, band)
        fout = open(outputFile,'w')
        
        # Open the objInfoFile for this filter band...
        with open(objInfoFile,'r') as fin:

            # Read header line...
            try:
                hin = fin.readline()
            except:
                print """%s is empty.  Skipping...""" % (objInfoFile)
                continue
            
            # Identify column in the objInfoFile CSV file corresponding to EXPNUM...
            hins=hin.upper().strip().split(',')
            try:
                expnumCol = hins.index('EXPNUM')
            except:
                print """Could not find 'EXPNUM' in header of %s...  Skipping""" % (objInfoFile)
                continue

            # Write header line...
            hout = hin_orig.upper().strip()+','+hin.upper().strip()+'\n'
            fout.write(hout)
            
            linecnt = 0
            for lin in fin:
            
                # Increment line count...
                linecnt += 1
                if ( (linecnt/1000.0 == int(linecnt/1000.0)) and (args.verbose > 1) ):
                    print '\r'+'Progress (lines read from '+objInfoFile+'):  ',linecnt,
                    sys.stdout.flush()
                    
                lins=lin.strip().split(',')
                expnum = int(lins[expnumCol])
                expnumString = str(expnum)

                try:
                    outputLine=expInfoDict[expnumString]
                except:
                    if args.verbose > 2:
                        print """Expnum %d not in %s.  Skipping...""" % (expnum, expInfoFile)
                else:
                    outputLine = outputLine+','+lin.strip()+'\n'
                    fout.write(outputLine)

        if args.verbose > 1:  
            print

    if args.verbose>0: print

    return 0


##################################
# 

def combine_stilts_internal_match_files(args):

    import os
    import sys
    import glob

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'combine_stilts_internal_match_files'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 


    fileNamePattern = 'y2a1_rawdata.svy1y2y3.sorted.tmp.g.csv.??.inmatch'
    fileNamePattern = args.fileNamePattern
    fileNameList = sorted(glob.glob(fileNamePattern))

    if args.verbose>2:  print fileNameList

    # Read header line from the "zeroth" input file from fileNameList
    fin = open(fileNameList[0])
    hin = fin.readline()
    fin.close()

    # Identify which column contains the GroupID...
    hins=hin.upper().strip().split(',')
    try:
        groupIDCol = hins.index('GROUPID')
    except:
        print """Could not find 'GROUPID' in header of %s...""" % (fileNameList[0])
        print """Exiting now..."""
        return 1

    # Open up the output file and write header to it...
    outputFileName = args.outputFileName
    fout = open(outputFileName,'w')
    fout.write(hin)

    groupIDOffset = 0

    # Loop through each of the input files from fileNameList...
    for fileName in fileNameList:

        print fileName

        # Create an array of newGroupIDs contained within this fileName...
        newGroupIDList = []

        # Open the input file...
        with open(fileName,'r') as fin:

            # Skip header line...
            try:
                next(fin)
            except:
                print """%s is empty.  Skipping...""" % (fileName)
                continue
            
            #print 
            linecnt = 0
            for lin in fin:
            
                # Increment line count...
                linecnt += 1
                if ( (linecnt/10000.0 == int(linecnt/10000.0)) and (args.verbose > 1) ):
                    print '\r'+'Progress (lines read from '+fileName+'):  ',linecnt,
                    sys.stdout.flush()

                # Update GroupID...
                lins=lin.strip().split(',')
                groupID = int(lins[groupIDCol])
                newGroupID = groupID + groupIDOffset
                newGroupIDList.append(newGroupID)
                lins[groupIDCol] = str(newGroupID)

                # Output updated line...
                outputLine = ''
                for value in lins:
                    outputLine = outputLine+value+','
                outputLine = outputLine[:-1]
                outputLine = outputLine+'\n'
                fout.write(outputLine)


        # After we finish with a given fileName,
        #  we update the groupIDOffset to be
        #  the value of the largest newGroupID
        #  in that fileName...
        groupIDOffset = max(newGroupIDList)
        print groupIDOffset

        if args.verbose > 1:  
            print
            print 

    # Close the output file...
    fout.close()
                          
    if args.verbose>0: print

    return 0


##################################
# 

def reformat_stilts_internal_match_files_for_gcm_exp_old(args):

    import os
    import sys
    import glob
    import numpy as np
    import pandas as pd
    import math
    import gcmPythonTools
    import datetime

    print 'Start of program: ',
    print datetime.datetime.now()

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'reformat_stilts_internal_match_files_for_gcm_exp'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 


    # Add blacklist, y3a1_blacklist, y3a1_graylist (non-photometric exposures), 
    #  and updated qa_summary table to "prune" non-photometric or otherwise
    #  questionable data from the fit.

    fileNamePattern = 'y2a1_rawdata.svy1y2y3.sorted.tmp.g.csv.ac.inmatch.first1000'
    #fileNamePattern = 'y2a1_rawdata.svy1y2y3.sorted.tmp.g.csv.??.inmatch'
    fileNamePattern = args.fileNamePattern
    fileNameList = sorted(glob.glob(fileNamePattern))

    if args.verbose>2:  print fileNameList

    # Open up the output file and write header to it...
    outputFileName = args.outputFileName
    fout = open(outputFileName,'w')
    hdr = gcmPythonTools.gcmInputFileHeader("stdstars")
    #if args.verbose > 1:
    #    print hdr
    fout.write(hdr)

    groupIDOffset = 0

    # Loop through each of the input files from fileNameList...
    for fileName in fileNameList:

        #print fileName

        # Create an array of newGroupIDs contained within this fileName...
        newGroupIDList = []

        # Open the input file as a pandas dataframe...
        print 'Start reading '+fileName+': ',
        print datetime.datetime.now()
        dataFrame = pd.pandas.read_csv(fileName)
        print 'Finished reading '+fileName+': ',
        print datetime.datetime.now()

        dataFrame['MAG'] = -2.5*np.log10(dataFrame['FLUX_PSF']/dataFrame['EXPTIME'])
        dataFrame['MAGERR'] = 1.086*dataFrame['FLUXERR_PSF']/dataFrame['FLUX_PSF']


        #print dataFrame

        # Update GroupID and groupIDOffset
        dataFrame['GroupID'] = dataFrame['GroupID'] + groupIDOffset
        groupIDOffset = max(dataFrame['GroupID'])
        
        # Identify unique set of groupIDs...
        uniqGroupIDArray = np.sort(dataFrame['GroupID'].unique())

        print 'Starting loop through individual groups: ',
        print datetime.datetime.now()

        for groupID in uniqGroupIDArray:

            #if args.verbose > 1:
                #print """Working on GroupID %d (max GroupID=%d)""" % (groupID, groupIDOffset)

            # Create a temporary data frame just containing 
            #  members of this particular group...
            mask = (dataFrame['GroupID'] == groupID)
            tmpDataFrame = dataFrame[mask]

            # How many members are in this particular group
            #  (especially after masking for bad expnums, etc.)?
            nmem = tmpDataFrame['GroupID'].size

            # Loop over all unique pairs in this group, outputting
            # info in a format digestible by the Global Calibrations
            # Module code GlobalZPSolverDC6.java...
        
            for i in range(0,nmem-1):

                regionid1        = tmpDataFrame['EXPNUM_old'].iloc[i]
                regionRaCenDeg1  = tmpDataFrame['EXPRA'].iloc[i]
                regionDecCenDeg1 = tmpDataFrame['EXPDEC'].iloc[i]
                regionQuality1   = 0
                starid1          = tmpDataFrame['GLOBJID'].iloc[i]
                ximage1          = tmpDataFrame['X_IMAGE'].iloc[i]
                yimage1          = tmpDataFrame['Y_IMAGE'].iloc[i]
                raDeg1           = tmpDataFrame['RA'].iloc[i]
                decDeg1          = tmpDataFrame['DEC'].iloc[i]
                mag1             = tmpDataFrame['MAG'].iloc[i]
                magErr1          = tmpDataFrame['MAGERR'].iloc[i]
                #flux1            = tmpDataFrame['FLUX_PSF'].iloc[i]
                #fluxErr1         = tmpDataFrame['FLUXERR_PSF'].iloc[i]
                #exptime1         = tmpDataFrame['EXPTIME'].iloc[i]

                #if ( (flux1 > 0.) and (exptime1 > 0.) ) :
                #    mag1 = -2.5*math.log10(flux1/exptime1)
                #else:
                #    mag1 = -9999.

                #if ( (flux1 > 0.) and (fluxErr1 >= 0.) ):
                #    magErr1 = 1.086*fluxErr1/flux1
                #else:
                #    magErr1 = -9999.

                for j in range(i+1,nmem):

                    regionid2        = tmpDataFrame['EXPNUM_old'].iloc[j]
                    regionRaCenDeg2  = tmpDataFrame['EXPRA'].iloc[j]
                    regionDecCenDeg2 = tmpDataFrame['EXPDEC'].iloc[j]
                    regionQuality2   = 0
                    starid2          = tmpDataFrame['GLOBJID'].iloc[j]
                    ximage2          = tmpDataFrame['X_IMAGE'].iloc[j]
                    yimage2          = tmpDataFrame['Y_IMAGE'].iloc[j]
                    raDeg2           = tmpDataFrame['RA'].iloc[j]
                    decDeg2          = tmpDataFrame['DEC'].iloc[j]
                    mag2             = tmpDataFrame['MAG'].iloc[j]
                    magErr2          = tmpDataFrame['MAGERR'].iloc[j]
                    #flux2            = tmpDataFrame['FLUX_PSF'].iloc[j]
                    #fluxErr2         = tmpDataFrame['FLUXERR_PSF'].iloc[j]
                    #exptime2         = tmpDataFrame['EXPTIME'].iloc[j]
                    
                    #if ( (flux2 > 0.) and (exptime2 > 0.) ) :
                    #    mag2 = -2.5*math.log10(flux2/exptime2)
                    #else:
                    #    mag2 = -9999.
                        
                    #if ( (flux2 > 0.) and (fluxErr2 >= 0.) ):
                    #    magErr2 = 1.086*fluxErr2/flux2
                    #else:
                    #    magErr2 = -9999.

                    #try:
                    #    sepArcsec = 3600.*gcmPythonTools.sepDegGet(raDeg1,decDeg1,raDeg2,decDeg2)
                    #except:
                    #    sepArcsec = -9999.
                    # Newly added - need to replace:
                    sepArcsec = 0.00

                    #if ( mag1 > -99. and magErr1 > 0. and mag2 > -99. and magErr2 > 0\
                    #         and sepArcsec >= 0. and regionid1 != regionid2 ):

                    outputLine = '%d %f %f %d %d %f %f %f %f    %d %f %f %d %d %f %f %f %f   %f\n' % \
                        (    regionid1, regionRaCenDeg1, regionDecCenDeg1, regionQuality1, \
                                 starid1, raDeg1, decDeg1, mag1, magErr1,\
                                 regionid2, regionRaCenDeg2, regionDecCenDeg2, regionQuality2, \
                                 starid2, raDeg2, decDeg2, mag2, magErr2,\
                                 sepArcsec)

                    fout.write(outputLine)

                    #endif

                #endfor_j

            #endfor_i           

        #print groupIDOffset

        #endfor_groupID

        del dataFrame

        #if args.verbose > 1:  
        #    print
        #    print 

    #endfor_filename

    # Close the output file...
    fout.close()
    
    print 'Finished: ',
    print datetime.datetime.now()
    if args.verbose>0: print

    return 0


##################################
# 

def reformat_stilts_internal_match_files_for_gcm_exp(args):

    import os
    import sys
    import glob
    import numpy as np
    import pandas as pd
    import math
    import gcmPythonTools
    import datetime

    print 'Start of program: ',
    print datetime.datetime.now()

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'reformat_stilts_internal_match_files_for_gcm_exp'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 


    # Add blacklist, y3a1_blacklist, y3a1_graylist (non-photometric exposures), 
    #  and updated qa_summary table to "prune" non-photometric or otherwise
    #  questionable data from the fit.

    #fileNamePattern = 'y2a1_rawdata.svy1y2y3.sorted.tmp.g.csv.ac.inmatch.first1000'
    fileNamePattern = 'y2a1_rawdata.svy1y2y3.sorted.tmp.g.csv.??.inmatch'
    fileNamePattern = args.fileNamePattern
    fileNameList = sorted(glob.glob(fileNamePattern))

    if args.verbose>2:  print fileNameList

    # Open up the output file and write header to it...
    outputFileName = args.outputFileName
    fout = open(outputFileName,'w')
    hdr = gcmPythonTools.gcmInputFileHeader("stdstars")
    #if args.verbose > 1:
    #    print hdr
    fout.write(hdr)

    # Loop through each of the input files from fileNameList...
    for fileName in fileNameList:

        #print fileName

        # Open the input file as a pandas dataframe...
        print 'Start reading '+fileName+': ',
        print datetime.datetime.now()
        dataFrame = pd.pandas.read_csv(fileName)
        print 'Finished reading '+fileName+': ',
        print datetime.datetime.now()

        dataFrame['MAG'] = -2.5*np.log10(dataFrame['FLUX_PSF']/dataFrame['EXPTIME'])
        dataFrame['MAGERR'] = 1.086*dataFrame['FLUXERR_PSF']/dataFrame['FLUX_PSF']

        maxGroupID = max(dataFrame['GroupID'])

        grouped = dataFrame.groupby('GroupID')


        print 'Starting loop through individual groups: ',
        print datetime.datetime.now()

        ntot = 0

        for name, group in grouped:

            if args.verbose > 1:
                if (ntot % 100 == 0):
                    print datetime.datetime.now(), 
                    print """Working on GroupID %d (max GroupID=%d)""" % (name, maxGroupID)
            ntot = ntot + 1

            # How many members are in this particular group
            #  (especially after masking for bad expnums, etc.)?
            nmem = group['GroupID'].size

            # Loop over all unique pairs in this group, outputting
            # info in a format digestible by the Global Calibrations
            # Module code GlobalZPSolverDC6.java...
        
            for i in range(0,nmem-1):

                regionid1        = group['EXPNUM_old'].iloc[i]
                regionRaCenDeg1  = group['EXPRA'].iloc[i]
                regionDecCenDeg1 = group['EXPDEC'].iloc[i]
                regionQuality1   = 0
                starid1          = group['GLOBJID'].iloc[i]
                ximage1          = group['X_IMAGE'].iloc[i]
                yimage1          = group['Y_IMAGE'].iloc[i]
                raDeg1           = group['RA'].iloc[i]
                decDeg1          = group['DEC'].iloc[i]
                mag1             = group['MAG'].iloc[i]
                magErr1          = group['MAGERR'].iloc[i]

                for j in range(i+1,nmem):

                    regionid2        = group['EXPNUM_old'].iloc[j]
                    regionRaCenDeg2  = group['EXPRA'].iloc[j]
                    regionDecCenDeg2 = group['EXPDEC'].iloc[j]
                    regionQuality2   = 0
                    starid2          = group['GLOBJID'].iloc[j]
                    ximage2          = group['X_IMAGE'].iloc[j]
                    yimage2          = group['Y_IMAGE'].iloc[j]
                    raDeg2           = group['RA'].iloc[j]
                    decDeg2          = group['DEC'].iloc[j]
                    mag2             = group['MAG'].iloc[j]
                    magErr2          = group['MAGERR'].iloc[j]

                    #try:
                    #    sepArcsec = 3600.*gcmPythonTools.sepDegGet(raDeg1,decDeg1,raDeg2,decDeg2)
                    #except:
                    #    sepArcsec = -9999.
                    # Newly added - need to replace:
                    sepArcsec = 0.00

                    #if ( mag1 > -99. and magErr1 > 0. and mag2 > -99. and magErr2 > 0\
                    #         and sepArcsec >= 0. and regionid1 != regionid2 ):

                    outputLine = '%d %f %f %d %d %f %f %f %f    %d %f %f %d %d %f %f %f %f   %f\n' % \
                        (    regionid1, regionRaCenDeg1, regionDecCenDeg1, regionQuality1, \
                                 starid1, raDeg1, decDeg1, mag1, magErr1,\
                                 regionid2, regionRaCenDeg2, regionDecCenDeg2, regionQuality2, \
                                 starid2, raDeg2, decDeg2, mag2, magErr2,\
                                 sepArcsec)

                    fout.write(outputLine)

                    #endif

                #endfor_j

            #endfor_i           

        #print groupIDOffset

        #endfor_groupID

        del dataFrame

        #if args.verbose > 1:  
        #    print
        #    print 

    #endfor_filename

    # Close the output file...
    fout.close()
    
    print 'Finished: ',
    print datetime.datetime.now()
    if args.verbose>0: print

    return 0


##################################
# 

def blacklist_query(args):

    import easyaccess as ea 
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'blacklist_query'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    connection=ea.connect('desoper') 
    
    query="""
        select * 
        from prod.blacklist 
        where expnum between %d and %d
        order by expnum, ccdnum
        """ % (args.expnumMin, args.expnumMax)
    
    if args.verbose>0: print query
    
    outputFile="""%s.blacklist.csv""" % (args.baseName)
    
    connection.query_and_save(query,outputFile)
    
    connection.close()
    
    if args.verbose>0: print

    return 0


##################################
# PROD.y3a1_blacklist contains the merger of the
#  original PROD.blacklist table and Eli Rykoff's
#  special, just-for-Y3A1 spread_model-based
#  y3a1_blacklist...

def y3a1_blacklist_query(args):

    import easyaccess as ea 
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'y3a1_blacklist_query'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    connection=ea.connect('desoper') 
    
    query="""
        select * 
        from prod.y3a1_blacklist 
        where expnum between %d and %d
        order by expnum, ccdnum
        """ % (args.expnumMin, args.expnumMax)
    
    if args.verbose>0: print query
    
    outputFile="""%s.y3a1_blacklist.csv""" % (args.baseName)
    
    connection.query_and_save(query,outputFile)
    
    connection.close()
    
    if args.verbose>0: print

    return 0


##################################
#
# Based on se_object_query_by_reqnum...

def image_query(args):

    import easyaccess as ea 
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'image_query'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    connection=ea.connect('desoper') 

    query="""
    with x as (
        select /*+materialize */ t.pfw_attempt_id, qa.expnum 
            from prod.proctag t, prod.qa_summary qa 
            where t.pfw_attempt_id=qa.pfw_attempt_id and 
                  tag='%s'
        ) 
    select i.pfw_attempt_id, i.band, i.expnum, i.ccdnum, 
           i.rac1, i.decc1, i.rac2, i.decc2, i.rac3, i.decc3, i.rac4, i.decc4, 
           i.ra_cent, i.dec_cent, i.racmax, i.racmin, i.crossra0 
    from prod.image i, x 
    where x.pfw_attempt_id=i.pfw_attempt_id and 
          i.filetype='red_starflat' and 
          i.expnum=x.expnum 
    order by i.expnum, i.ccdnum
    """ % (args.tag)


    if args.verbose>0: print query

    outputFile = args.outputFileName

    connection.query_and_save(query,outputFile)

    connection.close()

    if args.verbose>0: print

    return 0


##################################
#
# Based on se_object_query_by_reqnum...

def exp_image_cat_y3a1blacklist_query(args):

    import easyaccess as ea 
    import numpy as np 
    import pandas as pd
    import datetime

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'exp_image_cat_y3a1blacklist_query'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 


    # Open database connection...
    connection=ea.connect('desoper') 


    # Prepare and run query0...

    if args.verbose>0:
        print 
        print "Running query0..."
        print datetime.datetime.now()
        print 

    query0="""
        select e.expnum, t.unitname, t.reqnum, t.attnum, 
               e.radeg as EXPRA, e.decdeg as EXPDEC, e.exptime, e.airmass, e.band, 
               e.nite, e.mjd_obs, e.field, e.object, e.program, 
               qa.pfw_attempt_id, qa.source, qa.t_eff, qa.psf_fwhm, 
               qa.f_eff, qa.skybrightness, qa.b_eff, qa.cloud_apass, 
               qa.cloud_nomad, qa.n_apass, qa.n_nomad, qa.c_eff, 
               qa.skytilt, qa.astrom_sigma, qa.astrom_offset, 
               qa.astrom_ndets, qa.astrom_chi2, qa.nobjects, 
               qa.flag, qa.calnac
        from prod.exposure e, prod.proctag t, prod.qa_summary qa
        where t.pfw_attempt_id=qa.pfw_attempt_id and
              e.expnum=qa.expnum and 
              t.tag='%s' 
        """ % (args.tag)
    
    if args.verbose>0: print query0

    df0 = connection.query_to_pandas(query0).sort('EXPNUM')

    if args.verbose > 0:
        print datetime.datetime.now()

    if args.verbose > 0:
        print 
        print "Outputting query0 to test_df0.csv..."
        print datetime.datetime.now()

    df0.to_csv('test_df0.csv',index=False,float_format='%.8f')

    if args.verbose > 0:
        print datetime.datetime.now()


    # Prepare and run query1...

    if args.verbose>0:
        print 
        print "Running query1..."
        print datetime.datetime.now()
        print 

    query1="""
    with x as (
        select /*+materialize */ t.pfw_attempt_id, qa.expnum 
            from prod.proctag t, prod.qa_summary qa 
            where t.pfw_attempt_id=qa.pfw_attempt_id and 
                  tag='%s'
        ) 
    select i.pfw_attempt_id, i.band, i.expnum, i.ccdnum,
           i.gaina, i.gainb, i.rdnoisea, i.rdnoiseb,  
           i.rac1, i.decc1, i.rac2, i.decc2, i.rac3, i.decc3, i.rac4, i.decc4, 
           i.ra_cent, i.dec_cent, i.racmax, i.racmin, i.crossra0 
    from prod.image i, x 
    where x.pfw_attempt_id=i.pfw_attempt_id and 
          i.filetype='red_starflat' and 
          i.expnum=x.expnum 
    """ % (args.tag)

    if args.verbose>0: print query1

    df1 = connection.query_to_pandas(query1).sort(['EXPNUM','CCDNUM'])
    #df1 = pd.read_csv('test_df1.csv')

    if args.verbose > 0:
        print datetime.datetime.now()

    if args.verbose > 0:
        print 
        print "Outputting query1 to test_df1.csv..."
        print datetime.datetime.now()

    df1.to_csv('test_df1.csv',index=False,float_format='%.8f')

    if args.verbose > 0:
        print datetime.datetime.now()


    # Prepare and run query2...

    if args.verbose>0:
        print 
        print "Running query2..."
        print datetime.datetime.now()
        print 

    query2="""
    with x as (
        select /*+materialize */ t.pfw_attempt_id, qa.expnum 
            from prod.proctag t, prod.qa_summary qa 
            where t.pfw_attempt_id=qa.pfw_attempt_id and 
                  tag='%s'
        ) 
    select c.filename, c.band, c.expnum, c.ccdnum, c.objects 
    from prod.catalog c, x 
    where x.pfw_attempt_id=c.pfw_attempt_id and 
          c.filetype='cat_finalcut' and 
          c.expnum=x.expnum 
    """ % (args.tag)

    if args.verbose>0: print query2

    df2 = connection.query_to_pandas(query2).sort(['EXPNUM','CCDNUM'])

    if args.verbose > 0:
        print datetime.datetime.now()

    if args.verbose > 0:
        print 
        print "Outputting query2 to test_df2.csv..."
        print datetime.datetime.now()

    df2.to_csv('test_df2.csv',index=False,float_format='%.8f')

    if args.verbose > 0:
        print datetime.datetime.now()


    # Prepare and run query3...

    if args.verbose>0:
        print 
        print "Running query3..."
        print datetime.datetime.now()
        print 

    query3="""
        select * 
        from prod.y3a1_blacklist 
        """
    
    if args.verbose>0: print query3
    
    df3 = connection.query_to_pandas(query3).sort(['EXPNUM','CCDNUM'])

    if args.verbose > 0:
        print datetime.datetime.now()

    # Add EXPNUMCCDNUM column...
    df3.loc[:, 'EXPNUMCCDNUM'] = 100*df3['EXPNUM'] + df3['CCDNUM']

    if args.verbose > 0:
        print 
        print "Outputting query3 to test_df3.csv..."
        print datetime.datetime.now()

    df3.to_csv('test_df3.csv',index=False,float_format='%.8f')

    if args.verbose > 0:
        print datetime.datetime.now()

    # Create blacklist array based on EXPNUMCCDNUM column...
    blacklistArray = df3['EXPNUMCCDNUM'].values


    # Merge queries...
    # 
    # To drop columns with duplicate names, we make use of 
    #  recommendations from an e-mail from S. Allam (27 July 2016)
    #  and from these two stackoverflow questions:
    #  http://stackoverflow.com/questions/27313647/merging-two-pandas-dataframes-results-in-duplicate-columns 
    #  and http://stackoverflow.com/questions/19125091/pandas-merge-how-to-avoid-duplicating-columns

    # Merge df0 and df1...

    if args.verbose > 0:
        print 
        print "Merging data frames df0 and df1..."
        print datetime.datetime.now()

    df01 = df0.merge(df1, on=['EXPNUM'], suffixes=('', '_y')).sort(['EXPNUM'], ascending=True).reset_index(drop=True)
    to_drop = [colname for colname in df01 if colname.endswith('_y')]
    df01.drop(to_drop, axis=1, inplace=True)    

    if args.verbose > 0:
        print datetime.datetime.now()

    if args.verbose > 0:
        print 
        print "Outputting data frame df01 to test_df01.csv..."
        print datetime.datetime.now()

    df01.to_csv('test_df01.csv',index=False,float_format='%.8f')

    if args.verbose > 0:
        print datetime.datetime.now()


    # Merge df01 and df2...

    if args.verbose > 0:
        print 
        print "Merging data frames df01 and df2..."
        print datetime.datetime.now()

    df012 = df01.merge(df2, on=['EXPNUM','CCDNUM'], suffixes=('', '_y')).sort(['EXPNUM','CCDNUM'], ascending=True).reset_index(drop=True)
    to_drop = [colname for colname in df012 if colname.endswith('_y')]
    df012.drop(to_drop, axis=1, inplace=True)    

    if args.verbose > 0:
        print datetime.datetime.now()

    # Add EXPNUMCCDNUM column...
    df012.loc[:, 'EXPNUMCCDNUM'] = 100*df012['EXPNUM'] + df012['CCDNUM']

    if args.verbose > 0:
        print 
        print "Outputting data frame df012 to test_df012.csv..."
        print datetime.datetime.now()

    df012.to_csv('test_df012.csv',index=False,float_format='%.8f')

    if args.verbose > 0:
        print datetime.datetime.now()


    # Create blacklist mask...
    #  Note that this makes use of the inverse operator ("~")...
    mask = ~df012.EXPNUMCCDNUM.isin(blacklistArray)

    # Output the mask==True entries...
    outputFile = args.outputFileName
    if args.verbose > 0:
        print 
        print """Outputting blacklist-cleaned data frame df012 to %s...""" % (outputFile)
        print datetime.datetime.now()

    df012[mask].to_csv(outputFile, index=False)


    # Close database connection...
    connection.close()

    if args.verbose>0: print

    return 0


##################################

if __name__ == "__main__":
    main()

##################################


