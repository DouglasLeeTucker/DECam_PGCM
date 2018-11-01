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
    status = rawdata_exp_query(args)


##################################
# 

def rawdata_exp_query(args):

    import easyaccess as ea
    import pandas as pd
    import numpy as np 

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'rawdata_exp_query'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    # Updated query to get rid of any commas in the field, object, or program name
    #  (which could play havoc downstream)...
    query="""
        select e.expnum, 
               e.radeg as EXPRA, e.decdeg as EXPDEC, e.exptime, e.airmass, 
               e.band, e.nite, e.mjd_obs, 
               REPLACE(e.field,   ',', ' ') as field, 
               REPLACE(e.object,  ',', ' ') as object, 
               REPLACE(e.program, ',', ' ') as program, 
               t.pfw_attempt_id, t.tag, 
               qa.source, qa.t_eff, qa.psf_fwhm, qa.f_eff, 
               qa.skybrightness, qa.b_eff, qa.cloud_apass, qa.cloud_nomad,
               qa.n_apass, qa.n_nomad, qa.c_eff, qa.skytilt, qa.astrom_sigma,
               qa.astrom_offset, qa.astrom_ndets, qa.astrom_chi2, qa.nobjects,
               qa.flag, qa.calnac, qa.cloud_des, qa.n_des
        from prod.exposure e, prod.proctag t, prod.qa_summary qa
        where t.pfw_attempt_id=qa.pfw_attempt_id and
              e.expnum=qa.expnum and 
              t.tag='%s' and
              qa.expnum between %d and %d
        order by qa.expnum
        """ % (args.tag, args.expnumMin, args.expnumMax)
    

    #select e.expnum,
    #       e.radeg as EXPRA, e.decdeg as EXPDEC, e.exptime, e.airmass, 
    #       e.band, e.nite, e.mjd_obs, e.field, e.object, e.program,  
    #       q.reqnum,q.unitname,q.attnum,
    #       qa.t_eff,qa.c_eff,qa.b_eff,qa.skybrightness,qa.fwhm_asec 
    #from exposure e,finalcut_eval qa,proctag t,pfw_attempt q  
    #where t.pfw_attempt_id=qa.pfw_attempt_id and 
    #      e.expnum=qa.expnum and 
    #      q.id=t.pfw_attempt_id and 
    #      tag like 'Y3A2_MISC' 
    #order by e.expnum;

    #query="""
    #    select e.expnum, 
    #           e.radeg as EXPRA, e.decdeg as EXPDEC, e.exptime, e.airmass, 
    #           e.band, e.nite, e.mjd_obs, e.field, e.object, e.program, 
    #           t.pfw_attempt_id, t.tag, 
    #           f.t_eff, f.f_eff, f.b_eff, f.c_eff, 
    #           f.fwhm_asec, f.ellipticity, f.skybrightness, 
    #           f.cloud_apass, f.cloud_nomad, f.cloud_catalog, 
    #           f.n_apass, f.n_nomad,
    #           f.cloud_des, f.n_des
    #    from prod.exposure e, prod.proctag t, prod.finalcut_eval f
    #    where t.pfw_attempt_id=f.pfw_attempt_id and
    #          e.expnum=f.expnum and 
    #          t.tag='%s' and
    #          e.expnum between %d and %d
    #    order by e.expnum
    #    """ % (args.tag, args.expnumMin, args.expnumMax)
    
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
