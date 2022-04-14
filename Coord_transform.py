'''
This python file contains functions for coordinate system 
transformations.
-> Orbital Elements to Inertial frame
-> Inertial frame to Orbital Elements
-> Orbital Elements to Perifocal frame
-> Inertial fram to Perifocal frame

Author: Abinay Brown
Date: 26-Sep-2020
Email: abinayb@vt.edu

measuring units:
->length : meters
->velocity: meters/sec
->gravParam: m^3/s^2
->angle: deg (default)
'''
import numpy as np
import julian
from datetime import datetime as dt
import datetime

def ECItoLLA(state,yy,mm,dd,hh,mn,ss):
	omega =  7.292115486*(10**(-5))
	dt = datetime.datetime(yy,mm,dd,hh,mn,ss)
	jd = julian.to_jd(dt, fmt='jd')
	
	ut = (jd+0.5)- int(jd+0.5)
	jd = jd - ut
	du = jd - 2451545.0
	tu = du/36525

	GMST = 24110.54841 + (tu * (8640184.812866 + tu * (0.093104 - tu * 6.2*10**(-6))))
	GMST = np.mod(GMST + 86400.0*1.00273790934*ut,86400.0)
	gmst = 2*np.pi * GMST/86400.0
	
	'''
	#Solution 1
	d0 = datetime.datetime(2017,1,1,0,0,0)
	d1 = datetime.datetime(yy,mm,dd,hh,mn,ss)
	
	n0 = 366 + d0.toordinal() + (d0 - dt.fromordinal(d0.toordinal())).total_seconds()/(24*60*60)
	n1 = 366 + d1.toordinal() + (d1 - dt.fromordinal(d1.toordinal())).total_seconds()/(24*60*60)
	D = n1-n0
	g0 =100.833
	gmst = np.mod(g0 + (1.002737909349889*2*np.pi*D),2*np.pi) 
	print(gmst)
	'''
	
	rotMatrix = np.array([[np.cos(gmst), np.sin(gmst), 0],\
		[-np.sin(gmst), np.cos(gmst), 0]\
		,[0,0,1]])
	pos = np.array([state[0], state[1], state[2]])
	#vel = np.array([state[3], state[4], state[5]])
	
	pos_ECEF = rotMatrix.dot(pos.transpose())
	#vel_ECEF = rotMatrix.dot(vel.transpose())

	r_equ= np.linalg.norm(pos_ECEF[0:2])
	sinA = pos_ECEF[1]/r_equ
	cosA = pos_ECEF[0]/r_equ
	Long = np.rad2deg(np.arctan2(sinA,cosA));
	Lat  = np.rad2deg(np.arcsin(pos_ECEF[2]/np.linalg.norm(pos_ECEF)));
	Alt  = np.linalg.norm(pos_ECEF)-6378.137;
	return Lat, Long, Alt
	
	

def OEtoECI(a, e, i, argp, W, v, gravParam, deg = True):
	'''
	a - Semimajor axis (m)
	e - Eccentricity
	i - Inclination    (deg)
	argp - Argument of Periapsis (deg)
	W - Longitude of Ascending Node (deg)
	v - True Anomaly (deg)
	'''
	# conversion to radians
	if deg == True:
		i = np.deg2rad(i)
		argp = np.deg2rad(argp)
		W = np.deg2rad(W)
		v = np.deg2rad(v)

	marsMu = gravParam
	# Semi-Latus Rectum
	p = a*(1-(e**2))
	radius = p/(1 + (e*np.cos(v)))

	# Perifocal Coordinate System
	x = radius * np.cos(v)
	y = radius * np.sin(v)
	z = 0
	velx = -((marsMu/p)**(0.5)) * np.sin(v)             
	vely = ((marsMu/p)**(0.5)) * (e + np.cos(v))        
	velz = 0

	position = np.array([x, y, z])
	velocity = np.array([velx, vely, velz])

	# Creating rotation matrices for multiplication
	rightAscenMatrix = np.array([[np.cos(W), -np.sin(W), 0], [np.sin(W), np.cos(W), 0], [0, 0, 1]])
	inclinationMatrix = np.array([[1, 0, 0], [0, np.cos(i), -np.sin(i)], [0,np.sin(i),np.cos(i)]])
	argPeriapsisMatrix = np.array([[np.cos(argp), -np.sin(argp), 0], [np.sin(argp), np.cos(argp), 0], [0, 0, 1]])
	
	# Matrix multiplication of all the rotation matrices
	matrixConv = rightAscenMatrix.dot(inclinationMatrix.dot(argPeriapsisMatrix))
	
	# Multiplying the rotation matrix to the position and velocity vectors to get in ECI system
	pos = matrixConv.dot(position.transpose())
	vel = matrixConv.dot(velocity.transpose())        
	state = np.concatenate((pos, vel)) 
	return state

def ECItoOE(state, gravParam, deg = True):

	x= state[0]
	y= state[1]
	z= state[2]
	velx = state[3]
	vely = state[4]
	velz = state[5]
	marsMu = gravParam
	pos = np.array([x, y, z])
	vel = np.array([velx, vely, velz])
	
	posMag = np.linalg.norm(pos)
	velMag = np.linalg.norm(vel)
	# Calculating angular momentum
	h = np.cross(pos, vel)


	# Calculating perpendicular node vector
	K = np.array([0, 0, 1])
	n = np.cross(K, h)

	# Calculating eccentricity 
	eVector = ((((velMag**2)-(marsMu/posMag))*pos) - (np.linalg.norm((pos.dot(vel)))*vel))/marsMu
	e = np.linalg.norm(eVector)

	# Calculating Semi-Major axis
	# Specific Mechanical Energy (KE - PE)/m
	E = ((velMag**2)/2) - (marsMu/posMag)
	if (e-1) > 0.000001:
		a = -(marsMu/(2*E))
	else:
		a = float('inf')

	h_mag = np.linalg.norm(h)
	a= (h_mag**2)/(gravParam*(1-e**2))
	
	# Calculating Inclination
	i = np.arccos(h[2]/np.linalg.norm(h))
	i = np.rad2deg(i)

	# Calculating Argument of Periapsis
	argp = np.arccos(np.linalg.norm(n.dot(eVector))/(np.linalg.norm(n)*e))
	argp = np.rad2deg(argp)
	if (eVector[2] < 0):
		argp = 360 - argp

	# Calculating Longitude of Ascending Node
	W = np.arccos(n[0]/np.linalg.norm(n))
	W = np.rad2deg(W)
	if (eVector[1] < 0):
		W = 360 - W

	# Calculating True Anomaly
	v = np.arccos(np.linalg.norm(eVector.dot(pos))/(e*posMag))
	v = np.rad2deg(v)
	if (np.linalg.norm(pos.dot(vel)) < 0):
		v = 360 - v

	return a, e, i, argp, W, v

def ECItoPF(x, y, z, velx, vely, velz, gravParam, deg = True):

	marsMu = gravParam
	a, e, i, argp, W, v = ECItoOE(x, y, z, velx, vely, velz, marsMu)

	# conversion to radians
	if deg == True:
		i = np.deg2rad(i)
		argp = np.deg2rad(argp)
		W = np.deg2rad(W)
		v = np.deg2rad(v)
	# Semi-Latus Rectum
	p = a*(1-(e**2))
	radius = p/(1 + (e*np.cos(v)))

	# Perifocal Coordinate System
	x = radius * np.cos(v)
	y = radius * np.sin(v)
	z = 0
	velx = -((marsMu/p)**(0.5)) * np.sin(v)             
	vely = ((marsMu/p)**(0.5)) * (e + np.cos(v))        
	velz = 0

	return x, y, z, velx, vely, velz

def OEtoPF(a, e, i, argp, W, v, gravParam, deg = True):
	if deg == True:
		i = np.deg2rad(i)
		argp = np.deg2rad(argp)
		W = np.deg2rad(W)
		v = np.deg2rad(v)

	marsMu = gravParam
	# Semi-Latus Rectum
	p = a*(1-(e**2))
	radius = p/(1 + (e*np.cos(v)))

	# Perifocal Coordinate System
	x = radius * np.cos(v)
	y = radius * np.sin(v)
	z = 0
	velx = -((marsMu/p)**(0.5)) * np.sin(v)             
	vely = ((marsMu/p)**(0.5)) * (e + np.cos(v))        
	velz = 0

	return x, y, z, velx, vely, velz


print(ECItoOE( [-5529.203 * 1000, -2217.254 * 1000, 3399.353 * 1000, 3.049295 * 1000, 2.478617 * 1000, 6.576514 * 1000],5.972 * 6.674 * (10 ** -11) * (10 ** 24)));


