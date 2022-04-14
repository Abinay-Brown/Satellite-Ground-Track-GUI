import numpy as np
from Atm_Model import *


def diffEqn(state, time, mu, bc):
	earthRadius = 6378137	# meters
	earthOmega =  7.292115486*(10**(-5)) # rad/s 
	J2 = 1.0826*(10**(-3))
	# unpacking position data
	x, y, z = state[0:3]
	r = state[0:3]
	
	#unpacking velocity data
	velx, vely, velz = state[3:6]
	v = state[3:6]

	# Newtonian Gravity model

	accelNewtonX = -(mu/((np.linalg.norm(r))**3))*x
	accelNewtonY = -(mu/((np.linalg.norm(r))**3))*y
	accelNewtonZ = -(mu/((np.linalg.norm(r))**3))*z
	#print(accelNewtonX,accelNewtonY,accelNewtonZ)

	# J2 Perturbations 
	accelJ2X = -(3*mu*(earthRadius**2)*J2*(x)*((x**2)+(y**2)-(4*(z)**2)))
	accelJ2X = accelJ2X / (2*(np.linalg.norm(r))**7)
	accelJ2Y = -(3*mu*(earthRadius**2)*J2*(y)*((x**2)+(y**2)-(4*(z)**2)))
	accelJ2Y = accelJ2Y / (2*(np.linalg.norm(r))**7)
	accelJ2Z = -(3*mu*(earthRadius**2)*J2*(z)*((3*(x)**2)+(3*(y)**2)-(2*(z)**2)))
	accelJ2Z = accelJ2Z / (2*(np.linalg.norm(r))**7)
	#print(accelJ2X,accelJ2Y,accelJ2Z)    
	# Atmospheric drag

	velRel = []
	velRel.append(v[0] + (earthOmega*x))
	velRel.append(v[1] - (earthOmega*y))
	velRel.append(v[2])
	#************************************************
	# find atmospheric drag
	#************************************************
	rho = atmModel(np.linalg.norm(r))
	accelDragX = -(rho*np.linalg.norm(velRel)*(velRel[0])/bc)/2 
	accelDragY = -(rho*np.linalg.norm(velRel)*(velRel[1])/bc)/2
	accelDragZ = -(rho*np.linalg.norm(velRel)*(velRel[2])/bc)/2
	#print(accelDragX,accelDragY,accelDragZ)
	# Net Acceleration
	accelNetX = accelNewtonX + accelJ2X + accelDragX
	accelNetY = accelNewtonY + accelJ2Y + accelDragY
	accelNetZ = accelNewtonZ + accelJ2Z + accelDragZ

	stateDerivative = [velx, vely, velz, accelNetX, accelNetY, accelNetZ]
	print(stateDerivative);
	return stateDerivative
