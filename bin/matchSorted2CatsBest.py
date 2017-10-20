#!/usr/bin/env python
"""
matchSorted2CatsBest.py

Take two star catalogs (pre-sorted in ascending order by RA), read
star catalog 1 in line-by-line, and find the best match (by RA,DEC
coordinates, within a given tolerance)for the star catalog 1 entry 
from star catalog 2.

Examples:

matchSorted2CatsBest.py --help

matchSorted2CatsBest.py --inputCatFile1 inputCatFile1.csv --inputCatFile2 inputCatFile2.csv --racolCatFile1 0 --deccolCatFile1 1 --racolCatFile2 5 --deccolCatFile2 6 --outputMatchFile matched.csv --matchTolerance 2.0 --verbose 2 

"""

import math
import sys

##################################

def main():
    import argparse
    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputCatFile1', help='CSV file containing the star catalog 1', default='inputCatFile1.csv')
    parser.add_argument('--inputCatFile2', help='CSV file containing the star catalog 2', default='inputCatFile2.csv')
    parser.add_argument('--outputMatchFile', help='CSV file containing the output of the match', default='matched.csv')
    parser.add_argument('--racolCatFile1', help='zero-indexed column number of RA in inputCatFile1 (0, 1, 2, ...)', type=int, default=0)
    parser.add_argument('--deccolCatFile1', help='zero-indexed column number of DEC in inputCatFile1 (0, 1, 2, ...)', type=int, default=1)
    parser.add_argument('--racolCatFile2', help='zero-indexed column number of RA in inputCatFile2 (0, 1, 2, ...)', type=int, default=0)
    parser.add_argument('--deccolCatFile2', help='zero-indexed column number of DEC in inputCatFile2 (0, 1, 2, ...)', type=int, default=1)
    parser.add_argument('--matchTolerance', help='match tolerance (in arcsec)', type=float, default=1.0)
    parser.add_argument('--verbose', help='verbosity level of output to screen (0, 1, 2, ...)', type=int, default=0)
                        
    args = parser.parse_args()

    matchSorted2CatsBest(args)


##################################


# The upper-level match method...
def matchSorted2CatsBest(args):
    
    inputFile1=args.inputCatFile1
    inputFile2=args.inputCatFile2
    outputFile=args.outputMatchFile
    racol1=args.racolCatFile1
    deccol1=args.deccolCatFile1
    racol2=args.racolCatFile2
    deccol2=args.deccolCatFile2
    matchTolArcsec=args.matchTolerance
    verbose=args.verbose

    # Initialize lists of catalog 2 stars within the sliding window of RA...
    ra2_win=[]
    dec2_win=[]
    line2_win=[]
    
    # Open the output file for the catalog 1/catalog 2 star matches...
    fout=open(outputFile,'w')

    # Initialize running match id
    matchid=0
    
    # Open the CSV file for catalog 1...
    fin1=open(inputFile1)
    
    # Read header line of the catalog 1 CSV file...
    h1=fin1.readline()
    h1n=h1.strip().split(',')
    
    # Open the CSV file for catalog 2...
    fin2=open(inputFile2)
    
    # Read header line of catalog 2 CSV file...
    h2=fin2.readline()
    h2n=h2.strip().split(',')

    # Create and output header for the output CSV file...
    #  Note that the column names from the catalog 1 file
    #  now have a suffix of "_1", and that column names from
    #  the catalog 2 file now have a suffix of "_2".
    outputHeader='MATCHID'
    for colhead in h1n:
        outputHeader=outputHeader+','+colhead.upper()+'_1'
    for colhead in h2n:
        outputHeader=outputHeader+','+colhead.upper()+'_2'
    outputHeader=outputHeader+',SEPARCSEC\n'
    fout.write(outputHeader)
    

    # initialize some variables
    #  done1 = are we done reading catalog 1?
    #  done2 = are we done reading catalog 2?
    #  ra1, dec1 = current values of RA,DEC for star 1"
    #  ra2, dec2 = current values of RA,DEC for star 2"
    #  tol = sky angular separation tolerance (in degrees)
    #  tol2 = square of tol
    #  tolrawin = half-range of RA window (in degrees)
    #  linecnt1 = catolog 1 line count
    done1=0
    done2=0
    ra1=-9999.
    dec1=-9999.
    ra2=-9999.
    dec2=-9999.
    tol=matchTolArcsec/3600.0
    tol2=tol*tol
    tolrawin=3.*tol

    linecnt1=0
    
    # Loop through catalog 1 line-by-line...
    while (done1 == 0):

        # Increment line count from the catalog 1 file...
        linecnt1 += 1
        if ( (linecnt1/1000.0 == int(linecnt1/1000.0)) and (verbose > 1) ):
            print '\r'+'Progress (lines read from catalog 1):  ',linecnt1,
            sys.stdout.flush()

        # Read line from catalog 1 data file...
        l1=fin1.readline()

        # Are we done reading through the catalog 1 file yet?
        # If so, set done1=1 and ignore the rest of the loop;
        # otherwise, process the data line and continue with the
        # rest of the loop...
        if l1 == "":
            done1 = 1
            continue
        else:
            #line1 holds the whole line of information for this entry for
            # future use...
            line1=l1.strip()
            l1s=l1.strip().split(',')
            ra1=float(l1s[racol1])
            dec1=float(l1s[deccol1])
        #endif


        # Update the sliding RA window of catalog 2 stars...
        #  ... but only if ra2-ra1 <= tolrawin, 
        #  ... and only if we haven't previously finished
        #      reading the catalog 2 file...
        while ( (ra2-ra1 <= tolrawin) and (done2 == 0) ):

            # Read the next line from the catalog 2 data file...
            l2=fin2.readline()
        
            # if we have reached the end of the standard star file,
            # set done2=1 and skip the rest of this code block; 
            # otherwise, process the new line...
            if l2 ==  "":

                done2=1

            else:

                l2s=l2.strip().split(',')
                ra2_new=float(l2s[racol2])
                dec2_new=float(l2s[deccol2])

                # if the new star 2 RA (ra2_new) is at or above the 
                #  lower bound, add this catalog 2 star to the sliding RA 
                #  window...
                if ((ra2_new-ra1) >= -tolrawin):

                    # update values of ra2, dec2...
                    ra2 = ra2_new
                    dec2 = dec2_new

                    # add the catalog 2 star info to lists of ra, dec, and general
                    #  data for this sliding window...
                    ra2_win.append(ra2)
                    dec2_win.append(dec2)
                    line2_win.append(l2.strip())
            
                #endif
            
            #endif

        #endwhile -- inner "while" loop

        
        # Find the best "within radial tolerance" match between this catalog 1 star and 
        #  the set of catalog 2 stars within the sliding RA window... 
        
        cosd=math.cos(math.radians(dec1))
        best_flag=0
        best_delta2=tol2
        
        # Loop through all catalog 2 stars i in the sliding RA window for this
        #  catalog 1 star...
        for i in range(0,len(ra2_win)):

            delta2=(ra1-ra2_win[i])*(ra1-ra2_win[i])*cosd*cosd+(dec1-dec2_win[i])*(dec1-dec2_win[i])

            # Is the sky separation between the current catalog 1 star and catalog 2 star i 
            #  better than the previous "best" separation in this sliding RA window?
            #  (Note, initially the "best" separation is the radial tolerance for a good match.)
            if float(delta2) < float(best_delta2):
                best_flag = 1
                best_delta2 = delta2
                best_line2 = line2_win[i]
            #endif 

        #endfor

        if best_flag == 1:
            #output match

            # increment the running star id
            matchid += 1

            sepArcsec = 3600.*math.sqrt(best_delta2)
            
            # output line to the match file...
            outputLine = """%d,%s,%s,%.2f\n""" % (matchid,line1,best_line2,sepArcsec)
            fout.write(outputLine)

        #endif
            

        # Do some cleanup of the lists associated with the 
        #  sliding RA window...
        #  For each iteration of this while loop, we look at the "zeroth"
        #  catalog 2 star in the sliding RA window and remove it if it
        #  turns out to be now outside the RA tolerance window.
        while ( (len(ra2_win) > 1) and (ra1-ra2_win[0] > tolrawin) ):

            # Delete the lists associated with standard star
            #  (star "0" in the sliding RA window)...
            del ra2_win[0]
            del dec2_win[0]
            del line2_win[0]

        #endwhile -- inner "while" loop

    #endwhile -- outer "while" loop
        

    # close the input and output files...
    fin1.close()
    fin2.close()
    fout.close()

    return 0


##################################

if __name__ == "__main__":
    main()

##################################
