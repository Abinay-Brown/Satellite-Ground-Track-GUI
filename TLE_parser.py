'''
"TLE_parser"
Author: Abinay Brown
Email: abinay.joel@gmail.com
Description:
This software contains the read_TLE() function that takes
list of 3 strings of TLE data as input. It returns the 
initial conditions of position and velocity (cartesian). 
Ballistic coefficient, and epoch date and time in UTC.
'''
import numpy as np
from scipy import optimize
from datetime import datetime
from datetime import timedelta


def read_TLE(TLE_lines):
	
	name = TLE_lines[0];
	line1 = TLE_lines[1];
	line2 = TLE_lines[2];

	# Parsing line 1
	line1 = line1.split()
	catNum=(line1[1][0:-1])
	clss = line1[1][-1]
	if (clss == 'U'):
		clss = "Unclassified"
	elif (clss =='C'):
		clss = "Classified"
	elif (clss == 'S'):
		clss = "Secret"


	year = int(line1[3][0:2])
	day_frac = float(line1[3][2:-1])
	if year >=0 and year<=22:
		year = year + 2000
	elif year>=90 and year<=99:
		year = year + 1900
	date = datetime(year, 1, 1, 0, 0, 0);
	date = date + timedelta(days= day_frac)
	year = int(date.year)
	month = int(date.month)
	day = int(date.day)-1
	if (day == 0):
		day = 1
	hour = int(date.hour)
	minute = int(date.minute)
	second = int(date.second)
	if day == 0 and month == 1:
		day = 31
		month = 12

	
	if line1[6][0] == '-' or line1[6][0] == '+':
		power = int(line1[6][-2:])
		if line1[6][0] == '-':
			sign = -1
		if line1[6][0] == '+':
			sign = 1
		Bstar = sign*float('.' +line1[6][1:-2])*(10**power)
	else:
		power = int(line1[6][-2:])
		Bstar = float('.'+line1[6][:-2])*(10**power)
	if Bstar != 0:
		bc = 1/(12.741621*Bstar)
	else:
		bc = 1

	# Parsing line 2
	line2 = line2.split()
	inc = np.deg2rad(float(line2[2]))
	raan = np.deg2rad(float(line2[3]))
	ecc  = float('.'+ line2[4])
	argp = np.deg2rad(float(line2[5]))
	M = np.deg2rad(float(line2[6]))
	n = float(line2[7]) # rev/day

	#print(np.rad2deg(inc))
	#print(np.rad2deg(raan))
	#print(ecc)
	#print(np.rad2deg(argp))
	#print(np.rad2deg(M))
	# Calculating true anomaly
	def func(E):
		return E-ecc*np.sin(E)-M;
	def func_prime(E):
		return 1-ecc*np.cos(E);
	E = optimize.newton(func, M, func_prime, maxiter=50)
	v = 2*np.arctan((np.sqrt((1 + ecc)/(1 - ecc)))* np.tan(E/2))

	earthMu = 3.986004418 * (10**(14))
	t = 24*60*60/n
	a = ((earthMu*t**2)/(4.0*np.pi**2))**(1/3)
	r = a * (1- (ecc*np.cos(E)))

	# Perifocal coordinate system
	
	x   =  r*np.cos(v)
	y   =  r*np.sin(v)        
	z   = 0
	
	p   = a*(1 - ecc**2)
	vx  = -np.sqrt(earthMu/p) * np.sin(v)             
	vy  = np.sqrt(earthMu/p) * (ecc + np.cos(v))        
	vz  = 0

	# Perifocal to ECI coordinate system

	position = np.array([x, y, z])
	velocity = np.array([vx, vy, vz])
	
	# Creating rotation matrices for multiplication
	rightAscenMatrix = np.array([[np.cos(raan), -np.sin(raan), 0], [np.sin(raan), np.cos(raan), 0], [0, 0, 1]])
	inclinationMatrix = np.array([[1, 0, 0], [0, np.cos(inc), -np.sin(inc)], [0,np.sin(inc),np.cos(inc)]])
	argPeriapsisMatrix = np.array([[np.cos(argp), -np.sin(argp), 0], [np.sin(argp), np.cos(argp), 0], [0, 0, 1]])
	
	# Matrix multiplication of all the rotation matrices
	matrixConv = rightAscenMatrix.dot(inclinationMatrix.dot(argPeriapsisMatrix))
	
	# Multiplying the rotation matrix to the position and velocity vectors to get in ECI system
	pos = matrixConv.dot(position.transpose())
	vel = matrixConv.dot(velocity.transpose()) 
	return pos[0], pos[1], pos[2], vel[0], vel[1], vel[2], bc, year, month, day, hour, minute, second, t, name, catNum, clss



# run the code to test
TLE_lines=[]
TLE_lines.append("ISS (ZARYA)");
TLE_lines.append("1 25544U 98067A   21153.33610538  .00000648  00000-0  19959-4 0  9996");
TLE_lines.append("2 25544  51.6454  56.1809 0003452  55.5733  83.1214 15.48941137286297");
t0=0
#print(read_TLE(TLE_lines));
