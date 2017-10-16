#!/usr/bin/env python
"""
    uniq_RAsorted_rows.py

    Example:
    
    uniq_RAsorted_rows.py --help

    uniq_RAsorted_rows.py --inputFile input.csv --outputFile output.cvs --cleanup --verbose 2
    
    """

##################################

def main():

    import argparse
    import time

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputFile', help='name of the input CSV file', default='input.csv')
    parser.add_argument('--outputFile', help='name of the output CSV file', default='output.csv')
    parser.add_argument('--cleanup', help='remove intermediate files', default=False, action='store_true')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    args = parser.parse_args()

    if args.verbose > 0: print args

    status = uniq_RAsorted_rows(args)

    return status


##################################
# Sort by RA and keep only unique rows...
# Assumes that the column to be sorted has column name "RA"...

def uniq_RAsorted_rows(args):

    import numpy as np 
    import os
    import sys
    import datetime

    if args.verbose>0: 
        print 
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 'uniq_RAsorted_rows'
        print '* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *'
        print 

    inputFile = args.inputFile
    tmpFile = """%s.tmp.csv""" % (inputFile)
    outputFile = args.outputFile

    # Does input file exist?
    if os.path.isfile(inputFile)==False:
        print """uniq_RASorted_rows input file %s does not exist...  Exiting...""" % (inputFile)
        return 1

    # Copy inputFile to tmpFile, but adding an "RA_WRAP" column, which 
    #  steps around the RA=360/0deg wrap problems, as the first column
    #  in the tmpFile...
    # (assumes there is a column with name "RA")...
    with open(inputFile,'r') as fin:

        # Try to read header line...
        try:
            hin = fin.readline()
        except:
            print """%s is empty.  Exiting...""" % (inputFile)
            return 1

        # Try to identify which column from input file contains object RA...
        hins=hin.upper().strip().split(',')
        try:
            # Python is zero-indexed; unix sort & awk are one-indexed...
            raColumnID = hins.index('RA')
            unixRaColumnID = raColumnID + 1
        except:
            print """Could not find 'RA' in header of %s.  Exiting...""" % (inputFile)
            return 1

        hin = 'RA_WRAP,'+hin
        ftmp = open(tmpFile, 'w')
        ftmp.write(hin)

        linecnt = 0
        for lin in fin:
            
            # Increment line count...
            linecnt += 1
            if ( (linecnt/1000.0 == int(linecnt/1000.0)) and (args.verbose > 1) ):
                print '\r'+'Progress (lines read from '+inputFile+'):  ',linecnt,
                sys.stdout.flush()
                    
            lins=lin.strip().split(',')
            ra = float(lins[raColumnID])
            rawrap = ra
            if ra > 180.0:
                rawrap = ra - 360.

            outputLine = """%.8f,%s\n""" % (rawrap,lin.strip())
            ftmp.write(outputLine)

        ftmp.close()

    if args.verbose > 1:  
        print

    # Open outputFile and write header...
    fout = open(outputFile, 'w')
    fout.write(hin)
    fout.close()
        
    # Hard to beat Unix sort and awk for this sort of thing...
    #cmd = """awk 'NR>1' %s | sort -t, -g -k%d | uniq >> %s""" % (tmpFile,unixRaColumnID,outputFile)
    #cmd = """awk 'NR>1' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)

    totErrors = 0

    # Split into 10 parts, based on RA_WRAP, to help avoid memory issues:
    cmd = """awk -F, 'NR>1 && $1<-60.' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    cmd = """awk -F, 'NR>1 && ($1>=-60. && $1<-50.)' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    cmd = """awk -F, 'NR>1 && ($1>=-50. && $1<-40.)' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    cmd = """awk -F, 'NR>1 && ($1>=-40. && $1<-30.)' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    cmd = """awk -F, 'NR>1 && ($1>=-30. && $1<-20.)' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    cmd = """awk -F, 'NR>1 && ($1>=-20. && $1<-10.)' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    cmd = """awk -F, 'NR>1 && ($1>=-10. && $1<0.)' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    cmd = """awk -F, 'NR>1 && ($1>=0. && $1<10.)' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    cmd = """awk -F, 'NR>1 && ($1>=10. && $1<20.)' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    cmd = """awk -F, 'NR>1 && ($1>=20. && $1<30.)' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    cmd = """awk -F, 'NR>1 && ($1>=30. && $1<40.)' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    cmd = """awk -F, 'NR>1 && ($1>=40. && $1<50.)' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    cmd = """awk -F, 'NR>1 && ($1>=50. && $1<60.)' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    cmd = """awk -F, 'NR>1 && ($1>=60. && $1<70.)' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    cmd = """awk -F, 'NR>1 && ($1>=70. && $1<80.)' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    cmd = """awk -F, 'NR>1 && ($1>=80. && $1<90.)' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    cmd = """awk -F, 'NR>1 && ($1>=90. && $1<100.)' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    cmd = """awk -F, 'NR>1 && $1>=100.' %s | sort -t, -g -k1 | uniq >> %s""" % (tmpFile,outputFile)
    print datetime.datetime.now()
    print """Running:  %s""" % (cmd)
    status = os.system(cmd)
    print datetime.datetime.now()
    if status !=0:
        print """%s failed.""" % (cmd)
        totErrors += 1

    if args.verbose>0: print

    #  Finally:  cleanup of intermediate files...
    if args.cleanup==True:
        cmd = """rm -f %s""" % (tmpFile)
        status = os.system(cmd)
        if status !=0:
            print """%s failed.""" % (cmd)
            totErrors += 1

    if totErrors>0: 
        return 1
    else:
        return 0


##################################

if __name__ == "__main__":
    main()

##################################


