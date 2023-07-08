# class TLE:
#     def __init__(somename, name, company):
#         somename.name = name
#         somename.company = company
 
#     def show(somename):
#         print("Hello my name is " + somename.name +
#               " and I work in "+somename.company+".")
        
#     def displayError():
#         pass

#     def parseTLE():
#         pass    

import numpy as np
from datetime import datetime, timedelta
import pytz
import ephem

pdt = pytz.timezone('US/Pacific')

sqrt = np.sqrt
pi = np.pi
sin = np.sin
cos = np.cos

# Standard Gravitational parameter in km^3 / s^2 of Earth
GM = 398600.4418

def splitElem(tle):
    "Splits a two line element set into title and it's two lines with stripped lines"
    return map(lambda x: x.strip(), tle.split('\n'))

def checkValid(tle):
    "Checks with checksum to make sure element is valid"
    title, line1, line2 =  splitElem(tle)

    return line1[0] == '1' and line2[0] == '2' and \
           line1[2:7] == line2[2:7] and \
           int(line1[-1]) == doChecksum(line1) and int(line2[-1]) == doChecksum(line2)

def stringScientificNotationToFloat(sn):
    "Specific format is 5 digits, a + or -, and 1 digit, ex: 01234-5 which is 0.01234e-5"
    return 0.00001*float(sn[5]) * 10**int(sn[6:])

def eccentricAnomalyFromMean(mean_anomaly, eccentricity, initValue,
                             maxIter=500, maxAccuracy = 0.0001):
    """Approximates Eccentric Anomaly from Mean Anomaly
       All input and outputs are in radians"""
    mean_anomaly = mean_anomaly
    e0 = initValue
    for x in range(maxIter):
        e1 = e0 - (e0 - eccentricity * sin(e0) - mean_anomaly) / (1.0 - eccentricity * cos(e0))
        if (abs(e1-e0) < maxAccuracy):
            break
    return e1

def pretty_print(tle, printInfo = True, labels = True):
    "Returns commented information on a two line element"
    title, line1, line2 =  splitElem(tle)
    if not checkValid(tle):
        print("Invalid element.")
        return

    satellite_number                                        = int(line1[2:7])
    classification                                          = line1[7:8]
    international_designator_year                           = int(line1[9:11])
    international_designator_launch_number                  = int(line1[11:14])
    international_designator_piece_of_launch                = line1[14:17]
    epoch_year                                              = int(line1[18:20])
    epoch                                                   = float(line1[20:32])
    first_time_derivative_of_the_mean_motion_divided_by_two = float(line1[33:43])
    second_time_derivative_of_mean_motion_divided_by_six    = stringScientificNotationToFloat(line1[44:52])
    bstar_drag_term                                         = stringScientificNotationToFloat(line1[53:61])
    the_number_0                                            = float(line1[62:63])
    element_number                                          = float(line1[64:68])
    checksum1                                               = float(line1[68:69])

    satellite        = int(line2[2:7])
    inclination      = float(line2[8:16])
    right_ascension  = float(line2[17:25])
    eccentricity     = float(line2[26:33]) * 0.0000001
    argument_perigee = float(line2[34:42])
    mean_anomaly     = float(line2[43:51])
    mean_motion      = float(line2[52:63])
    revolution       = float(line2[63:68])
    checksum2        = float(line2[68:69])

    # Inferred Epoch date
    year = 2000 + epoch_year if epoch_year < 70 else 1900 + epoch_year
    epoch_date = datetime(year=year, month=1, day=1, tzinfo=pytz.utc) + timedelta(days=epoch-1) # Have to subtract one day to get correct midnight

    # Time difference of now from epoch, offset in radians
    diff = datetime.now().replace(tzinfo=pytz.utc) + timedelta(hours=8) - epoch_date # Offset for PDT
    diff_seconds = 24*60*60*diff.days + diff.seconds + 1e-6*diff.microseconds # sec
    print ("Time offset: %s" % diff)
    motion_per_sec = mean_motion * 2*pi / (24*60*60) # rad/sec
    print ("Radians per second: %g" % motion_per_sec)
    offset = diff_seconds * motion_per_sec #rad
    print ("Offset to apply: %g" % offset)
    mean_anomaly += offset * 180/pi % 360

    # Inferred period
    day_seconds = 24*60*60
    period = day_seconds * 1./mean_motion

    # Inferred semi-major axis (in km)
    semi_major_axis = ((period/(2*pi))**2 * GM)**(1./3)

    # Inferred true anomaly
    eccentric_anomaly = eccentricAnomalyFromMean(mean_anomaly * pi/180, eccentricity, mean_anomaly * pi/180)
    true_anomaly = 2*np.arctan2(sqrt(1+eccentricity) * sin(eccentric_anomaly/2.0), sqrt(1-eccentricity) * cos(eccentric_anomaly/2.0))
    # Convert to degrees
    eccentric_anomaly *= 180/pi
    true_anomaly *= 180/pi

    tle_rec = ephem.readtle(title, line1, line2)
    tle_rec.compute()

    if (printInfo):
        print ("----------------------------------------------------------------------------------------")
        print ("Satellite Name                                            = %s" % title)
        print ("Satellite number                                          = %g (%s)" % (satellite_number, "Unclassified" if classification == 'U' else "Classified"))
        print ("International Designator                                  = YR: %02d, LAUNCH #%d, PIECE: %s" % (international_designator_year, international_designator_launch_number, international_designator_piece_of_launch))
        print ("Epoch Date                                                = %s  (YR:%02d DAY:%.11g)" % (epoch_date.strftime("%Y-%m-%d %H:%M:%S.%f %Z"), epoch_year, epoch))
        print ("First Time Derivative of the Mean Motion divided by two   = %g" % first_time_derivative_of_the_mean_motion_divided_by_two)
        print ("Second Time Derivative of Mean Motion divided by six      = %g" % second_time_derivative_of_mean_motion_divided_by_six)
        print ("BSTAR drag term                                           = %g" % bstar_drag_term)
        print ("The number 0                                              = %g" % the_number_0)
        print ("Element number                                            = %g" % element_number)
        print
        print ("Inclination [Degrees]                                     = %g" % inclination)
        print ("Right Ascension of the Ascending Node [Degrees]           = %g" % right_ascension)
        print ("Eccentricity                                              = %g" % eccentricity)
        print ("Argument of Perigee [Degrees]                             = %g" % argument_perigee)
        print ("Mean Anomaly [Degrees] Anomaly                            = %g" % mean_anomaly)
        print ("Eccentric Anomaly                                         = %g" % eccentric_anomaly)
        print ("True Anomaly                                              = %g" % true_anomaly)
        print ("Mean Motion [Revs per day] Motion                         = %g" % mean_motion)
        print ("Period                                                    = %s" % timedelta(seconds=period))
        print ("Revolution number at epoch [Revs]                         = %g" % revolution)

        print("----------------------------------------------------------------------------------------")
        print ("Semi-Major Axis (a)                                       = %gkm" % semi_major_axis)
        print ("Eccentricity    (e)                                       = %g" % eccentricity)
        print ("Inclination     (i)                                       = %g" % inclination)
        print ("Argument of Periapsis (w)                                 = %g" % argument_perigee)
        print ("Right Ascension of the Ascending Node (Ω)                 = %g" % right_ascension)
        print ("True Anomaly (v)                                          = %g" % true_anomaly)
        print ("----------------------------------------------------------------------------------------")
        print ("Longitude                                                 = %s" % tle_rec.sublong)
        print ("Latitude                                                  = %s" % tle_rec.sublat)
        print("----------------------------------------------------------------------------------------")


def doChecksum(line):
    """The checksums for each line are calculated by adding the all numerical digits on that line, including the 
       line number. One is added to the checksum for each negative sign (-) on that line. All other non-digit 
       characters are ignored.
       @note this excludes last char for the checksum thats already there."""
    return sum(map(int, filter(lambda c: c >= '0' and c <= '9', line[:-1].replace('-','1')))) % 10

elem1 = """ISS (ZARYA)
1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927
2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"""

name = input("Paste in Name of Satellite: ")
lineOne = input("Paste in TLE Line One: ")
lineTwo = input("Paste in TLE Line Two: ")

constructedTLE = name + "\n" + lineOne + "\n" + lineTwo

# dragon = """DRAGON CRS-2            
# 1 39115U 13010A   13062.62492353  .00008823  00000-0  14845-3 0   188
# 2 39115  51.6441 272.5899 0012056 334.2535  68.5574 15.52501943   306"""

pretty_print(constructedTLE)