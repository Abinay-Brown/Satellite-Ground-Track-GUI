import numpy as np
from numpy import pi

R_e = 6378.137  # Earth Eq radius in [km]
f_par = 0.00335281066  # Earth flatness parameter (f) used for conversion from ECF to Latitude and Longitude


# JulianDate function, converts UTC time vector Time = [ YY , MM , DD , Hour , Min , Sec ] into Julian Date
# Input parameter is a list holding year, month, day, hour, minute, second (floats)
# Output is a float number giving the current Julian Date (float)
def JulianDate(Time):
    status = 0
    ls = [[1959, 12, 31], [1960, 12, 31], [1961, 7, 31], [1961, 12, 31], [1963, 10, 31], [1963, 12, 31], [1964, 3, 31],
          [1964, 8, 31], [1964, 12, 31], [1965, 2, 28], [1965, 6, 30], [1965, 8, 31], [1965, 12, 31], [1968, 1, 31],
          [1971, 12, 31], [1972, 6, 30], [1972, 12, 31], [1973, 12, 31], [1974, 12, 31], [1975, 12, 31],
          [1976, 12, 31], [1977, 12, 31], [1978, 12, 31], [1979, 12, 31], [1981, 6, 30], [1982, 6, 30], [1983, 6, 30],
          [1985, 6, 30], [1987, 12, 31], [1989, 12, 31], [1990, 12, 31], [1992, 6, 30], [1993, 6, 30], [1994, 6, 30],
          [1995, 12, 31], [1997, 6, 30], [1998, 12, 31], [2005, 12, 31], [2008, 12, 31], [2012, 6, 30], [2015, 6, 30],
          [2016, 12, 31]]
    # from https://cdf.gsfc.nasa.gov/html/CDFLeapSeconds.txt

    if (len(Time) == 6):
        LS = 60
        for i in ls:
            if (Time[0] == i[0] and Time[1] == i[1] and Time[2] == i[2] and Time[3] == 23 and Time[4] == 59):
                LS = 61
                break

        JD = 1721013.5 + 367.0 * Time[0] - float(int((7.0 / 4.0) * (Time[0] + float(int(((Time[1] + 9.0) / 12.0)))))) + \
             float(int((275.0 * Time[1] / 9.0))) + Time[2] + (60.0 * Time[3] + Time[4] + Time[5] / LS) / 1440.0

    else:
        print("Invalid array size in JD function - insert Time in UTC format [ YY , MM , DD , Hour , Min , Sec ]")
        JD = 0.0

    return JD
