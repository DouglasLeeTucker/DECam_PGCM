#!/usr/bin/env python
"""
matchSortedStdObsCats.py

Take a standard star catalog (pre-sorted in ascending order by RA) 
and an observed catalog (also pre-sorted in ascending order by RA) 
and match them by their RA,DEC coordinates.

Note that it is possible for multiple entries from the observed 
catalog to be matched to a single entry from the standard star 
catalog.

CAVEAT:  current algorithm has issues if the separation between 
         individual standard stars is less than the match 
         tolerance radius.  (See section in the code on finding 
         the first good match.)

Examples:

matchSortedStdObsCats.py --help

matchSortedStdObsCats.py --inputStdStarCatFile standard_stars_all_id6_pypsmFormat.csv --inputObsCatFile obsquery-20141011-g.list.csv --racolStdStarCatFile 0 --deccolStdStarCatFile 1 --racolObsCatFile 5 --deccolObsCatFile 6 --outputMatchFile matched-20141011-g.list.csv --matchTolerance 2.0 --verbose 2 

"""

import math
import sys

##################################

def main():
    import argparse
    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputStdStarCatFile', help='CSV file containing the standard star catalog to be used in the match', default='standardstarcat.csv')
    parser.add_argument('--inputObsCatFile', help='CSV file containing the catalog of observations to be used in the match', default='obsquery.csv')
    parser.add_argument('--outputMatchFile', help='CSV file containing the output of the match', default='matched.csv')
    parser.add_argument('--racolStdStarCatFile', help='zero-indexed column number of RA in inputStdStarCatFile (0, 1, 2, ...)', type=int, default=0)
    parser.add_argument('--deccolStdStarCatFile', help='zero-indexed column number of DEC in inputStdStarCatFile (0, 1, 2, ...)', type=int, default=1)
    parser.add_argument('--racolObsCatFile', help='zero-indexed column number of RA in inputObsCatFile (0, 1, 2, ...)', type=int, default=0)
    parser.add_argument('--deccolObsCatFile', help='zero-indexed column number of DEC in inputObsCatFile (0, 1, 2, ...)', type=int, default=1)
    parser.add_argument('--matchTolerance', help='match tolerance (in arcsec)', type=float, default=1.0)
    parser.add_argument('--verbose', help='verbosity level of output to screen (0, 1, 2, ...)', type=int, default=0)
                        
    args = parser.parse_args()

    matchSortedStdObsCats(args)


##################################


# The upper-level match method...
def matchSortedStdObsCats(args):
    
    f1=args.inputStdStarCatFile
    f2=args.inputObsCatFile
    outfile=args.outputMatchFile
    stdracol=args.racolStdStarCatFile
    stddeccol=args.deccolStdStarCatFile
    obsracol=args.racolObsCatFile
    obsdeccol=args.deccolObsCatFile
    matchTolArcsec=args.matchTolerance
    verbose=args.verbose

    # Initialize "dictionaries"...
    # Each element of a "dictionary" is associated with a standard star.
    # Each element is a list of information from the potential matches from the
    #  observed data.
    raDict=[]
    decDict=[]
    obslineDict=[]
    
    # Initialize lists of standard stars...
    # These are actually lists of standards within a given sliding window of RA.
    stdra_win=[]
    stddec_win=[]
    stdline_win=[]
    
    # Open the output file for the standard star/observed star matches...
    ofd=open(outfile,'w')

    # Initialize running match id
    matchid=0
    
    # Open the standard star CSV file...
    fd1=open(f1)
    
    # Read header line of standard star CSV file...
    h1=fd1.readline()
    h1n=h1.strip().split(',')
    
    # Open CSV file of observed data...
    fd2=open(f2)
    
    # Read header line of observed data CSV file...
    h2=fd2.readline()
    h2n=h2.strip().split(',')

    # Create and output header for the output CSV file...
    #  Note that the column names from the standard star file
    #  now have a suffix of "_1", and that column names from
    #  the observed star file now have a suffix of "_2".
    outputHeader='MATCHID'
    for colhead in h1n:
        outputHeader=outputHeader+','+colhead.upper()+'_1'
    for colhead in h2n:
        outputHeader=outputHeader+','+colhead.upper()+'_2'
    outputHeader=outputHeader+'\n'
    ofd.write(outputHeader)
    

    # initialize some variables
    #  done_std = "are we done reading the standard stars file?"
    #  done_obs = "are we done reading the observations file?"
    #  stdra, stddec = "current values for standard star RA,DEC"
    #  obsra, obsdec = "current values for observed star RA,DEC"
    #  tol = sky angular separation tolerance (in degrees)
    #  tol2 = square of tol
    #  tolrawin = half-range of RA window (in degrees)
    #  linecnt = "line count"
    done_std=0
    done_obs=0
    stdra=-999
    stddec=-999
    obsra=-999
    obsdec=-999
    tol=matchTolArcsec/3600.0
    tol2=tol*tol
    tolrawin=3.*tol

    linecnt=0
    
    # Loop through file of observed data...
    while (done_obs == 0):

        # Increment line count from observed data file...
        linecnt += 1
        if ( (linecnt/1000.0 == int(linecnt/1000.0)) and (verbose > 1) ):
            print '\r'+'Progress (lines read from observed catalog):  ',linecnt,
            sys.stdout.flush()

        # Read line from observed data file...
        l2=fd2.readline()

        # Are we done reading through the file of observed data yet?
        # If so, set done_obs=1 and ignore the rest of the loop;
        # otherwise, process the data line and continue with the
        # rest of the loop...
        if l2 == "":
            done_obs = 1
            continue
        else:
            #obsline2 holds the whole line of information for this entry for
            # future use...
            obsline2=l2.strip()
            l2s=l2.strip().split(',')
            obsra=float(l2s[obsracol])
            obsdec=float(l2s[obsdeccol])
        #endif


        # Update the sliding RA window of standard stars...
        #  ... but only if stdra-obsra <= tolrawin, 
        #  ... and only if we haven't previously finished
        #      reading the standard star file...
        while ( (stdra-obsra <= tolrawin) and (done_std == 0) ):

            # Read the next line from the standard star file...
            l1=fd1.readline()
        
            # if we have reached the end of the standard star file,
            # set done_std=1 and skip the rest of this code block; 
            # otherwise, process the new line...
            if l1 ==  "":

                done_std=1

            else:

                l1s=l1.strip().split(',')
                stdra_new=float(l1s[stdracol])
                stddec_new=float(l1s[stddeccol])

                # if the new standard star RA (stdra_new) is at or above the 
                #  lower bound, add this standard star to the sliding RA 
                #  window...
                if ((stdra_new-obsra) >= -tolrawin):

                    # update values of stdra, stddec...
                    stdra = stdra_new
                    stddec = stddec_new

                    # add the standard star info to lists of ra, dec, and general
                    #  data for this sliding window...
                    stdra_win.append(stdra)
                    stddec_win.append(stddec)
                    stdline_win.append(l1.strip())
            
                    # initialize lists for possible observed star/standard star 
                    #  matches and add these (still empty) lists to "dictionaries" 
                    #  associated with this sliding window of standard stars...
                    raDict.append([])
                    decDict.append([])
                    obslineDict.append([])
            
                #endif
            
            #endif

        #endwhile -- inner "while" loop

        
        # Find the first good match (not necessarily the best match) between this
        # observed star and the set of standard stars within the sliding RA
        # window... 
        # (We might want to revisit this choice -- i.e., of first match vs. best
        #  match -- in the future.)
        
        cosd=math.cos(math.radians(obsdec))

        # Loop through all standards stars i in the sliding RA window for that
        #  observed star...
        for i in range(0,len(stdra_win)):

            delta2=(obsra-stdra_win[i])*(obsra-stdra_win[i])*cosd*cosd+(obsdec-stddec_win[i])*(obsdec-stddec_win[i])

            # Is the sky position of standard star i (in the sliding RA window)
            #  within the given radial tolerance of the observed star?  If so, 
            #  add the observed info to that standard star's dictionaries...
            if float(delta2) < float(tol2):
                raDict[i].append(obsra)
                decDict[i].append(obsdec)
                obslineDict[i].append(obsline2)
                # if we found one match, we take it and break out of this "for"
                #  loop...
                break
            #endif 

        #endfor


        # Do some cleanup of the lists and "dictionaries" associated with the 
        #  sliding RA window and output matches to output file...
        #  For each iteration of this while loop, we look at the "zeroth"
        #  standard star in the sliding RA window and remove it if it
        #  turns out to be now outside the RA tolerance window.
        while ( (len(stdra_win) > 1) and (obsra-stdra_win[0] > tolrawin) ):

            # Loop through all the observations matched with this standard 
            # star...
            # (Note that many standard stars may have zero matches...)
            for j in range(0,len(raDict[0])):

                # increment the running star id
                matchid += 1
                        
                # output line to the match file...
                outputLine = """%d,%s,%s\n""" % (matchid,stdline_win[0],obslineDict[0][j])
                ofd.write(outputLine)

            #endfor

            # Delete the dictionaries associated with this standard star
            #  (star "0" in the sliding RA window)...
            del raDict[0]
            del decDict[0]
            del obslineDict[0]
    
            # Delete the lists associated with standard star
            #  (star "0" in the sliding RA window)...
            del stdra_win[0]
            del stddec_win[0]
            del stdline_win[0]

        #endwhile -- inner "while" loop

    #endwhile -- outer "while" loop
        

    # Do some cleanup of the lists and "dictionaries" associated with the sliding 
    #  RA window after reading last line of observed data file and output matches
    #  to output file...
    while (len(stdra_win) > 0):

        # Loop through all the observations matched with this standard star...
        # (Note that many standard stars may have zero matches...)
        for j in range(0,len(raDict[0])):

            # increment the running star id
            matchid += 1

            # output line to the match file...
            outputLine = """%d,%s,%s\n""" % (matchid,stdline_win[0],obslineDict[0][j])
            ofd.write(outputLine)

        #endfor

        # Delete the dictionaries associated with this standard star
        #  (star "0" in the sliding RA window)...
        del raDict[0]
        del decDict[0]
        del obslineDict[0]
    
        # Delete the lists associated with standard star
        #  (star "0" in the sliding RA window)...
        del stdra_win[0]
        del stddec_win[0]
        del stdline_win[0]
            
    #endwhile
        

    # close the input and output files...
    fd1.close()
    fd2.close()
    ofd.close()

    return 0


##################################

if __name__ == "__main__":
    main()

##################################
