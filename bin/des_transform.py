#!/usr/bin/env python
"""
    des_transform.py

    Example:
    
    des_transform --help

    des_transform.py --inputFile match_fgcm_in_SDSSNorthSouthGalCapDR13_dtucker_des.csv--verbose 2

    """

##################################

def main():

    import argparse
    import time

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputFile', help='name of the input CSV file', default='input.csv')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    status = des_transform(args)

    return status


##################################
# des_transform
#

def des_transform(args):

    import numpy as np 
    import os
    import sys
    import datetime
    import fitsio
    import pandas as pd
    import matplotlib.pyplot as plt

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'des_transform'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    inputFile = args.inputFile
    
    # Read selected columns from inputFile
    columns = ['MAG_PSF_G_1','MAG_PSF_R_1','MAG_PSF_I_1','MAG_PSF_Z_1','MAG_PSF_Y_1','G_SDSS_2','R_SDSS_2','I_SDSS_2','Z_SDSS_2']
    print datetime.datetime.now()
    print """Reading in selected columns from %s...""" % (inputFile)
    df = pd.read_csv(inputFile, usecols=columns)
    print datetime.datetime.now()

    # Rename columns...
    df.rename(columns={'MAG_PSF_G_1':'g_des','MAG_PSF_R_1':'r_des','MAG_PSF_I_1':'i_des','MAG_PSF_Z_1':'z_des','MAG_PSF_Y_1':'Y_des','G_SDSS_2':'g_sdss','R_SDSS_2':'r_sdss','I_SDSS_2':'i_sdss','Z_SDSS_2':'z_sdss'},inplace=True)

    # Add color columns...
    df.loc[:,'gr_des'] = df.loc[:,'g_des'] - df.loc[:,'r_des']
    df.loc[:,'iz_des'] = df.loc[:,'i_des'] - df.loc[:,'z_des']
    df.loc[:,'gi_des'] = df.loc[:,'g_des'] - df.loc[:,'i_des']
    df.loc[:,'gr_sdss'] = df.loc[:,'g_sdss'] - df.loc[:,'r_sdss']
    df.loc[:,'iz_sdss'] = df.loc[:,'i_sdss'] - df.loc[:,'z_sdss']
    df.loc[:,'gi_sdss'] = df.loc[:,'g_sdss'] - df.loc[:,'i_sdss']


    # Convert to numpy ndarrays...
    g_des = df['g_des'].as_matrix()
    r_des = df['r_des'].as_matrix()
    i_des = df['i_des'].as_matrix()
    z_des = df['z_des'].as_matrix()
    Y_des = df['Y_des'].as_matrix()
    gr_des = df['gr_des'].as_matrix()
    iz_des = df['iz_des'].as_matrix()
    gi_des = df['gi_des'].as_matrix()
    g_sdss = df['g_sdss'].as_matrix()
    r_sdss = df['r_sdss'].as_matrix()
    i_sdss = df['i_sdss'].as_matrix()
    z_sdss = df['z_sdss'].as_matrix()
    gr_sdss = df['gr_sdss'].as_matrix()
    iz_sdss = df['iz_sdss'].as_matrix()
    gi_sdss = df['gi_sdss'].as_matrix()

    # No longer need the data frame...
    del df


    bandList = ['g','r','i','z','Y']

    for band in bandList:

        print band

        if band=='g':

            mag_des = g_des
            mag_sdss = g_sdss
            color_des = gr_des
            color_sdss = gr_sdss
            color_lo = -0.6
            color_hi =  2.0
        
        elif band=='r':

            mag_des = r_des
            mag_sdss = r_sdss
            color_des = gr_des
            color_sdss = gr_sdss
            color_lo = -0.6
            color_hi =  2.0
        
        elif band=='i':

            mag_des = i_des
            mag_sdss = i_sdss
            color_des = iz_des
            color_sdss = iz_sdss
            color_lo = -0.5
            color_hi =  1.1
        
        elif band=='z':

            mag_des = z_des
            mag_sdss = z_sdss
            color_des = iz_des
            color_sdss = iz_sdss
            color_lo = -0.5
            color_hi =  1.1
        
        elif band=='Y':

            mag_des = Y_des
            mag_sdss = z_sdss
            color_des = iz_des
            color_sdss = iz_sdss
            color_lo = -0.5
            color_hi =  1.1
        
        else:
            
            continue


        # First, fit the SDSS->DES tranformation and the "normal" color...
        # We want to consider 1st, 2nd, and 3rd order polynomials...

        for norder in range(1,4):

            print norder

            color = color_sdss
            dmag =  mag_des - mag_sdss
            mask1 = ( (mag_des > 0.0) & (mag_sdss > 0.0) )
            mask2 = ( (color >= color_lo) & (color < color_hi) )
            mask3 = ( mag_des > -1000. ) # initially a dummy value

            niter = 0
            for i in range(10):

                niter = i + 1
                mask = mask1 & mask2 & mask3
                dmag_masked = dmag[np.where(mask)]
                color_masked = color[np.where(mask)]
        
                if norder == 1:
                    p1,p0 = np.polyfit(color_masked,dmag_masked,1)
                    p = np.poly1d([p1,p0])
                    resid = dmag - (p1*color + p0)

                elif norder == 2:
                    p2,p1,p0 = np.polyfit(color_masked,dmag_masked,2)
                    p = np.poly1d([p2,p1,p0])
                    resid = dmag - (p2*color*color + p1*color + p0)

                elif norder == 3:
                    p3,p2,p1,p0 = np.polyfit(color_masked,dmag_masked,3)
                    p = np.poly1d([p3,p2,p1,p0])
                    resid = dmag - (p3*color*color*color + p2*color*color + p1*color + p0)

                else:
                    continue


                resid_masked = resid[np.where(mask)]
                mean_offset = np.mean(resid[np.where(mask)])
                median_offset = np.median(resid[np.where(mask)])
                stddev = np.std(resid[np.where(mask)])
                ntot = resid[np.where(mask)].size
                mask3 = (resid > -3.0*stddev) & (resid < 3.0*stddev)


                # Prepare QA plots...
                if norder == 1:
                    title = """Iter%d: (p1,p0)=(%.4f,%.4f)\n  median offset: %.3f, mean offset: %.3f, stddev: %.3f, ntot: %d""" % \
                        (niter, p1, p0, median_offset, mean_offset, stddev, ntot)
                elif norder == 2:
                    title = """Iter%d: (p2,p1,p0)=(%.4f,%.4f,%.4f)\n  median offset: %.3f, mean offset: %.3f, stddev: %.3f, ntot: %d""" % \
                        (niter, p2, p1, p0, median_offset, mean_offset, stddev, ntot)
                elif norder == 3:
                    title = """Iter%d: (p3,p2,p1,p0)=(%.4f,%.4f,%.4f,%.4f)\n  median offset: %.3f, mean offset: %.3f, stddev: %.3f, ntot: %d""" % \
                        (niter, p3, p2, p1, p0, median_offset, mean_offset, stddev, ntot)
                else:
                    continue

                if args.verbose>0:
                    print title
                    print np.poly1d(p)

                xlinlo = int(100*color_lo)
                xlinhi = int(100*color_hi)

                # Plot1:  1D histogram...
                plotName1 = """%s_sdss_to_des_hist1d.norder%d.niter%d.png""" % (band, norder, niter)
                plt.hist(resid_masked,bins=100)
                plt.title(title)
                plt.xlabel('Delta_mag Residuals to Fit')
                plt.ylabel('Number')
                plt.grid(True)
                plt.savefig(plotName1,format='png')
                plt.clf()

                # Plot2:  2D histogram...
                plotName2 = """%s_sdss_to_des_hist2d_1.norder%d.niter%d.png""" % (band, norder, niter)
                plt.hist2d(color_masked, dmag_masked, bins=100)
                plt.title(title)
                plt.xlabel('color')
                plt.ylabel('Delta_mag')
                cb = plt.colorbar()
                cb.set_label('Number')
                x = 0.01*np.linspace(xlinlo, xlinhi)
                plt.plot(x,p(x),'m-')
                plt.savefig(plotName2,format='png')
                plt.clf()

                # Plot3:  2D histogram (residuals)...
                plotName3 = """%s_sdss_to_des_hist2d_2.norder%d.niter%d.png""" % (band, norder, niter)
                plt.hist2d(color_masked, resid_masked, bins=100)
                plt.title(title)
                plt.xlabel('color')
                plt.ylabel('Delta_mag Residuals to Fit')
                cb = plt.colorbar()
                cb.set_label('Number')
                x = 0.01*np.linspace(xlinlo, xlinhi)
                y = 0.00*np.linspace(xlinlo, xlinhi)
                plt.plot(x,y,'m-')
                plt.savefig(plotName3,format='png')
                plt.clf()

                print

            print
            print


        # Second, fit the SDSS->DES tranformation versus the (g-i) color...
        # We want to consider 1st, 2nd, and 3rd order polynomials...

        color_des = gi_des
        color_sdss = gi_sdss
        color_lo = -0.5
        color_hi =  3.5
        for norder in range(1,4):

            print norder

            color = color_sdss
            dmag =  mag_des - mag_sdss
            mask1 = ( (mag_des > 0.0) & (mag_sdss > 0.0) )
            mask2 = ( (color >= color_lo) & (color < color_hi) )
            mask3 = ( mag_des > -1000. ) # initially a dummy value

            niter = 0
            for i in range(10):

                niter = i + 1
                mask = mask1 & mask2 & mask3
                dmag_masked = dmag[np.where(mask)]
                color_masked = color[np.where(mask)]
        
                if norder == 1:
                    p1,p0 = np.polyfit(color_masked,dmag_masked,1)
                    p = np.poly1d([p1,p0])
                    resid = dmag - (p1*color + p0)

                elif norder == 2:
                    p2,p1,p0 = np.polyfit(color_masked,dmag_masked,2)
                    p = np.poly1d([p2,p1,p0])
                    resid = dmag - (p2*color*color + p1*color + p0)

                elif norder == 3:
                    p3,p2,p1,p0 = np.polyfit(color_masked,dmag_masked,3)
                    p = np.poly1d([p3,p2,p1,p0])
                    resid = dmag - (p3*color*color*color + p2*color*color + p1*color + p0)

                else:
                    continue


                resid_masked = resid[np.where(mask)]
                mean_offset = np.mean(resid[np.where(mask)])
                median_offset = np.median(resid[np.where(mask)])
                stddev = np.std(resid[np.where(mask)])
                ntot = resid[np.where(mask)].size
                mask3 = (resid > -3.0*stddev) & (resid < 3.0*stddev)


                # Prepare QA plots...
                if norder == 1:
                    title = """Iter%d: (p1,p0)=(%.4f,%.4f)\n  median offset: %.3f, mean offset: %.3f, stddev: %.3f, ntot: %d""" % \
                        (niter, p1, p0, median_offset, mean_offset, stddev, ntot)
                elif norder == 2:
                    title = """Iter%d: (p2,p1,p0)=(%.4f,%.4f,%.4f)\n  median offset: %.3f, mean offset: %.3f, stddev: %.3f, ntot: %d""" % \
                        (niter, p2, p1, p0, median_offset, mean_offset, stddev, ntot)
                elif norder == 3:
                    title = """Iter%d: (p3,p2,p1,p0)=(%.4f,%.4f,%.4f,%.4f)\n  median offset: %.3f, mean offset: %.3f, stddev: %.3f, ntot: %d""" % \
                        (niter, p3, p2, p1, p0, median_offset, mean_offset, stddev, ntot)
                else:
                    continue

                if args.verbose>0:
                    print title
                    print np.poly1d(p)

                xlinlo = int(100*color_lo)
                xlinhi = int(100*color_hi)

                # Plot1:  1D histogram...
                plotName1 = """%s_sdss_to_des_hist1d.norder%d.niter%d.g-i.png""" % (band, norder, niter)
                plt.hist(resid_masked,bins=100)
                plt.title(title)
                plt.xlabel('Delta_mag Residuals to Fit')
                plt.ylabel('Number')
                plt.grid(True)
                plt.savefig(plotName1,format='png')
                plt.clf()

                # Plot2:  2D histogram...
                plotName2 = """%s_sdss_to_des_hist2d_1.norder%d.niter%d.g-i.png""" % (band, norder, niter)
                plt.hist2d(color_masked, dmag_masked, bins=100)
                plt.title(title)
                plt.xlabel('color')
                plt.ylabel('Delta_mag')
                cb = plt.colorbar()
                cb.set_label('Number')
                x = 0.01*np.linspace(xlinlo, xlinhi)
                plt.plot(x,p(x),'m-')
                plt.savefig(plotName2,format='png')
                plt.clf()

                # Plot3:  2D histogram (residuals)...
                plotName3 = """%s_sdss_to_des_hist2d_2.norder%d.niter%d.g-i.png""" % (band, norder, niter)
                plt.hist2d(color_masked, resid_masked, bins=100)
                plt.title(title)
                plt.xlabel('color')
                plt.ylabel('Delta_mag Residuals to Fit')
                cb = plt.colorbar()
                cb.set_label('Number')
                x = 0.01*np.linspace(xlinlo, xlinhi)
                y = 0.00*np.linspace(xlinlo, xlinhi)
                plt.plot(x,y,'m-')
                plt.savefig(plotName3,format='png')
                plt.clf()

                print 

            print
            print

        # Third, fit the DES->SDSS tranformation and the "normal" color...
        # We want to consider 1st, 2nd, and 3rd order polynomials...

        for norder in range(1,4):

            print norder

            color = color_des
            dmag =  mag_sdss - mag_des
            mask1 = ( (mag_des > 0.0) & (mag_sdss > 0.0) )
            mask2 = ( (color >= color_lo) & (color < color_hi) )
            mask3 = ( mag_sdss > -1000. ) # initially a dummy value

            niter = 0
            for i in range(10):

                niter = i + 1
                mask = mask1 & mask2 & mask3
                dmag_masked = dmag[np.where(mask)]
                color_masked = color[np.where(mask)]
        
                if norder == 1:
                    p1,p0 = np.polyfit(color_masked,dmag_masked,1)
                    p = np.poly1d([p1,p0])
                    resid = dmag - (p1*color + p0)

                elif norder == 2:
                    p2,p1,p0 = np.polyfit(color_masked,dmag_masked,2)
                    p = np.poly1d([p2,p1,p0])
                    resid = dmag - (p2*color*color + p1*color + p0)

                elif norder == 3:
                    p3,p2,p1,p0 = np.polyfit(color_masked,dmag_masked,3)
                    p = np.poly1d([p3,p2,p1,p0])
                    resid = dmag - (p3*color*color*color + p2*color*color + p1*color + p0)

                else:
                    continue


                resid_masked = resid[np.where(mask)]
                mean_offset = np.mean(resid[np.where(mask)])
                median_offset = np.median(resid[np.where(mask)])
                stddev = np.std(resid[np.where(mask)])
                ntot = resid[np.where(mask)].size
                mask3 = (resid > -3.0*stddev) & (resid < 3.0*stddev)


                # Prepare QA plots...
                if norder == 1:
                    title = """Iter%d: (p1,p0)=(%.4f,%.4f)\n  median offset: %.3f, mean offset: %.3f, stddev: %.3f, ntot: %d""" % \
                        (niter, p1, p0, median_offset, mean_offset, stddev, ntot)
                elif norder == 2:
                    title = """Iter%d: (p2,p1,p0)=(%.4f,%.4f,%.4f)\n  median offset: %.3f, mean offset: %.3f, stddev: %.3f, ntot: %d""" % \
                        (niter, p2, p1, p0, median_offset, mean_offset, stddev, ntot)
                elif norder == 3:
                    title = """Iter%d: (p3,p2,p1,p0)=(%.4f,%.4f,%.4f,%.4f)\n  median offset: %.3f, mean offset: %.3f, stddev: %.3f, ntot: %d""" % \
                        (niter, p3, p2, p1, p0, median_offset, mean_offset, stddev, ntot)
                else:
                    continue

                if args.verbose>0:
                    print title
                    print np.poly1d(p)

                xlinlo = int(100*color_lo)
                xlinhi = int(100*color_hi)

                # Plot1:  1D histogram...
                plotName1 = """%s_des_to_sdss_hist1d.norder%d.niter%d.png""" % (band, norder, niter)
                plt.hist(resid_masked,bins=100)
                plt.title(title)
                plt.xlabel('Delta_mag Residuals to Fit')
                plt.ylabel('Number')
                plt.grid(True)
                plt.savefig(plotName1,format='png')
                plt.clf()

                # Plot2:  2D histogram...
                plotName2 = """%s_des_to_sdss_hist2d_1.norder%d.niter%d.png""" % (band, norder, niter)
                plt.hist2d(color_masked, dmag_masked, bins=100)
                plt.title(title)
                plt.xlabel('color')
                plt.ylabel('Delta_mag')
                cb = plt.colorbar()
                cb.set_label('Number')
                x = 0.01*np.linspace(xlinlo, xlinhi)
                plt.plot(x,p(x),'m-')
                plt.savefig(plotName2,format='png')
                plt.clf()

                # Plot3:  2D histogram (residuals)...
                plotName3 = """%s_des_to_sdss_hist2d_2.norder%d.niter%d.png""" % (band, norder, niter)
                plt.hist2d(color_masked, resid_masked, bins=100)
                plt.title(title)
                plt.xlabel('color')
                plt.ylabel('Delta_mag Residuals to Fit')
                cb = plt.colorbar()
                cb.set_label('Number')
                x = 0.01*np.linspace(xlinlo, xlinhi)
                y = 0.00*np.linspace(xlinlo, xlinhi)
                plt.plot(x,y,'m-')
                plt.savefig(plotName3,format='png')
                plt.clf()

                print

            print
            print


        # Fourth, fit the DES->SDSS tranformation versus the (g-i) color...
        # We want to consider 1st, 2nd, and 3rd order polynomials...

        color_des = gi_des
        color_sdss = gi_sdss
        color_lo = -0.5
        color_hi =  3.5
        for norder in range(1,4):

            print norder

            color = color_des
            dmag =  mag_sdss - mag_des
            mask1 = ( (mag_des > 0.0) & (mag_sdss > 0.0) )
            mask2 = ( (color >= color_lo) & (color < color_hi) )
            mask3 = ( mag_sdss > -1000. ) # initially a dummy value

            niter = 0
            for i in range(10):

                niter = i + 1
                mask = mask1 & mask2 & mask3
                dmag_masked = dmag[np.where(mask)]
                color_masked = color[np.where(mask)]
        
                if norder == 1:
                    p1,p0 = np.polyfit(color_masked,dmag_masked,1)
                    p = np.poly1d([p1,p0])
                    resid = dmag - (p1*color + p0)

                elif norder == 2:
                    p2,p1,p0 = np.polyfit(color_masked,dmag_masked,2)
                    p = np.poly1d([p2,p1,p0])
                    resid = dmag - (p2*color*color + p1*color + p0)

                elif norder == 3:
                    p3,p2,p1,p0 = np.polyfit(color_masked,dmag_masked,3)
                    p = np.poly1d([p3,p2,p1,p0])
                    resid = dmag - (p3*color*color*color + p2*color*color + p1*color + p0)

                else:
                    continue


                resid_masked = resid[np.where(mask)]
                mean_offset = np.mean(resid[np.where(mask)])
                median_offset = np.median(resid[np.where(mask)])
                stddev = np.std(resid[np.where(mask)])
                ntot = resid[np.where(mask)].size
                mask3 = (resid > -3.0*stddev) & (resid < 3.0*stddev)


                # Prepare QA plots...
                if norder == 1:
                    title = """Iter%d: (p1,p0)=(%.4f,%.4f)\n  median offset: %.3f, mean offset: %.3f, stddev: %.3f, ntot: %d""" % \
                        (niter, p1, p0, median_offset, mean_offset, stddev, ntot)
                elif norder == 2:
                    title = """Iter%d: (p2,p1,p0)=(%.4f,%.4f,%.4f)\n  median offset: %.3f, mean offset: %.3f, stddev: %.3f, ntot: %d""" % \
                        (niter, p2, p1, p0, median_offset, mean_offset, stddev, ntot)
                elif norder == 3:
                    title = """Iter%d: (p3,p2,p1,p0)=(%.4f,%.4f,%.4f,%.4f)\n  median offset: %.3f, mean offset: %.3f, stddev: %.3f, ntot: %d""" % \
                        (niter, p3, p2, p1, p0, median_offset, mean_offset, stddev, ntot)
                else:
                    continue

                if args.verbose>0:
                    print title
                    print np.poly1d(p)

                xlinlo = int(100*color_lo)
                xlinhi = int(100*color_hi)

                # Plot1:  1D histogram...
                plotName1 = """%s_des_to_sdss_hist1d.norder%d.niter%d.g-i.png""" % (band, norder, niter)
                plt.hist(resid_masked,bins=100)
                plt.title(title)
                plt.xlabel('Delta_mag Residuals to Fit')
                plt.ylabel('Number')
                plt.grid(True)
                plt.savefig(plotName1,format='png')
                plt.clf()

                # Plot2:  2D histogram...
                plotName2 = """%s_des_to_sdss_hist2d_1.norder%d.niter%d.g-i.png""" % (band, norder, niter)
                plt.hist2d(color_masked, dmag_masked, bins=100)
                plt.title(title)
                plt.xlabel('color')
                plt.ylabel('Delta_mag')
                cb = plt.colorbar()
                cb.set_label('Number')
                x = 0.01*np.linspace(xlinlo, xlinhi)
                plt.plot(x,p(x),'m-')
                plt.savefig(plotName2,format='png')
                plt.clf()

                # Plot3:  2D histogram (residuals)...
                plotName3 = """%s_des_to_sdss_hist2d_2.norder%d.niter%d.g-i.png""" % (band, norder, niter)
                plt.hist2d(color_masked, resid_masked, bins=100)
                plt.title(title)
                plt.xlabel('color')
                plt.ylabel('Delta_mag Residuals to Fit')
                cb = plt.colorbar()
                cb.set_label('Number')
                x = 0.01*np.linspace(xlinlo, xlinhi)
                y = 0.00*np.linspace(xlinlo, xlinhi)
                plt.plot(x,y,'m-')
                plt.savefig(plotName3,format='png')
                plt.clf()

                print 

            print
            print



        print
        print
        print 

    return 0


##################################

if __name__ == "__main__":
    main()

##################################

