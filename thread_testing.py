import multiprocessing 
import threading
from Satellites import Satellite
from Ground_Track import *
from datetime import datetime
from datetime import timedelta
import datetime
import time
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
  
if __name__ == "__main__": 
	
	'''
	line1= "CartoSat-f" 
	line2= "1 43111U 18004A   20361.84719535  .00000716  00000-0  37053-4 0  9993"
	line3= "2 43111  97.3618  56.9921 0005899 102.5229  19.5777 15.19302019163928"
	'''
	'''
	line1= "COSMOS 2361 [-]" 
	line2= "1 25590U 98076A   20361.42762958  .00000033  00000-0  19320-4 0  9991"
	line3= "2 25590  82.9345 243.8287 0030228 264.1302 269.5101 13.72988143102408"
	'''

	'''
	line1 = "MOLNIYA 2-9"
	line2 = "1 07276U 74026A   20361.24330652 -.00000118  00000-0  00000-0 0  9999"
	line3 = "2 07276  64.1628 308.3540 6738500 286.6354  13.5045  2.45096145236169"
	'''
	
	line1= "ISS (ZARYA)" 
	line2= "1 25544U 98067A   20364.69439405  .00000748  00000-0  21559-4 0  9990"
	line3= "2 25544  51.6460 101.0967 0001179 165.2966 341.7300 15.49234108262321"
	
	'''
	line1= "NOAA 15" 
	line2= "1 25338U 98030A   20361.69032464  .00000041  00000-0  35661-4 0  9997"
	line3= "2 25338  98.7005  25.5997 0011323  66.3979 293.8388 14.26005702176629"
	'''
	ISS = Satellite(line1,line2,line3)
	ECI, LLA, TS = ISS.run_solver()
	track2(LLA)


    