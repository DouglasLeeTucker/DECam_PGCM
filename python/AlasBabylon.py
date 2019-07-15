# Convert between sexagesimal units and decimal degs
#     and between sexagesimal units and decimal radians
#Created on 18.06.2009

import math
import string


#---------------------------------------------------------------------------

def HMSToDeg(hms):
    hmslist = hms.split(':')
    if hmslist[0][0] == '-' :
        hmsSign = -1
    else:
        hmsSign = 1
    hh = abs(int(hmslist[0]))
    mm = abs(int(hmslist[1]))
    ss = abs(float(hmslist[2]))
    hours = hh + (mm + ss/60.)/60.
    deg = 15.*hours
    deg = hmsSign*deg
    return deg

#---------------------------------------------------------------------------

def degToHMS(deg):
    if deg < 0:
        degSign = -1
        deg = abs(deg)
    else:
        degSign = 1
    hours = deg/15.
    # round to a certain number of digits to avoid int->float conversion difficulties
    hh = int( round( hours, 5) )
    mm = int( round( (hours - hh)*60., 5) )
    ss = abs( round( (hours - hh - mm/60.)*3600., 3) )
    hms = "%2.2d:%2.2d:%06.3f" % (hh, mm, ss)
    if degSign == -1:
        hms = '-' + hms
    return hms

#---------------------------------------------------------------------------

def DMSToDeg(dms):
    dmslist = dms.split(':')
    if dmslist[0][0] == '-' :
        dmsSign = -1
    else:
        dmsSign = 1
    dd = abs(int(dmslist[0]))
    mm = abs(int(dmslist[1]))
    ss = abs(float(dmslist[2]))
    deg = dd + (mm + ss/60.)/60.
    deg = dmsSign*deg
    return deg

#---------------------------------------------------------------------------

def degToDMS(deg):
    if deg < 0:
        degSign = -1
        deg = abs(deg)
    else:
        degSign = 1
    # round to a certain number of digits to avoid int->float conversion difficulties
    dd = int( round( deg, 5) )
    mm = int( round( (deg - dd)*60., 5) )
    ss = abs( round( (deg - dd - mm/60.)*3600., 3) )
    dms = "%2.2d:%2.2d:%06.3f" % (dd, mm, ss)
    if degSign == -1:
        dms = '-' + dms
    else:
        dms = '+' + dms
    return dms

#---------------------------------------------------------------------------

def HMSToRad(hms):
    deg = HMSToDeg(hms)
    rad = math.radians(deg)
    return rad

#---------------------------------------------------------------------------

def radToHMS(rad):
    deg = math.degrees(rad)
    hms = degToHMS(deg)
    return hms

#---------------------------------------------------------------------------

def DMSToRad(dms):
    deg = DMSToDeg(dms)
    rad = math.radians(deg)
    return rad

#---------------------------------------------------------------------------

def radToDMS(rad):
    deg = math.degrees(rad)
    dms = degToDMS(deg)
    return dms

#---------------------------------------------------------------------------
