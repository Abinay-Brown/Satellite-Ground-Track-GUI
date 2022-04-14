

from scipy.integrate import ode
from scipy.integrate import odeint
from TLE_parser import *
from Ground_Track import track, animated_track
from Coord_transform import ECItoLLA
from simulator import *
from Atm_Model import *
from ODE import *
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from datetime import timedelta
import datetime
from tools import JulianDate
class Satellite:

	# Satellite Name
	name = "" 
	# initial conditions
	state = []
	# Date and time of epoch
	epoch = []
	# orbit period
	period = 0 
	earthMu = 0
	bc = 0

	# Solution Variables
	time_index = 0
	timestamp = 0
	ECI =0
	LLA=0

	def __init__(self, TLE_Line1, TLE_Line2, TLE_Line3):

		# setting up initial conditions
		TLE_lines=[]
		TLE_lines.append(TLE_Line1);
		TLE_lines.append(TLE_Line2);
		TLE_lines.append(TLE_Line3);

		x,y,z,velx,vely,velz,bc,yy,mm,dd,hh,mn,sc,t,name,cn,clss= read_TLE(TLE_lines);
		self.name = name
		self.state = [x, y, z, velx, vely, velz]
		self.epoch = [yy, mm, dd, hh, mn, sc]
		self.period = t
		self.earthMu = 3.986004418 * (10**(14)) 
		self.bc = bc


		
		
		#print(self.epoch)

	def run_solver(self, numOrbits, stepSize = 5):
		earthMu = self.earthMu
		bc = self.bc
		state = self.state
		t0=0
		yy = self.epoch[0]
		mm = self.epoch[1]
		dd = self.epoch[2]
		hh = self.epoch[3]
		mn = self.epoch[4]
		sc = self.epoch[5]

		solver = ode(diffEqn)
		solver.set_integrator('dop853', atol = 10**(-15), rtol = 10**(-15)) 
		solver.set_f_params(earthMu, bc)
		solver.set_initial_value(state, t0)
		
		
		stepSize = 3 # seconds
		now = datetime.datetime.utcnow().replace(microsecond=0)
		end = now + timedelta(seconds = numOrbits*self.period)
		epoch = datetime.datetime(yy,mm,dd,hh,mn,sc).replace(microsecond=0)
		diff = (end-epoch)
		t1 = int(diff.total_seconds())  
		
		# Creating solution variables
		t = np.arange(t0,t1, stepSize)
		timestamp = np.empty((len(t),6))
		sol = np.empty((len(t), 6))
		LLA_sol = np.empty((len(t), 3))
		
		# initializing solution variables
		lat,lon,alt =ECItoLLA(state,yy,mm,dd,hh,mn,sc)
		sol[0] = state
		LLA_sol[0] = [lat, lon, alt]
		timestamp[0]= [yy,mm,dd,hh,mn,sc]
		
		dt = datetime.datetime(yy,mm,dd,hh,mn,sc)
		k = 1
		ind = 0


		while solver.successful() and k < len(t):
			solver.integrate(t[k])
			sol[k]= solver.y
			dt = dt + timedelta(seconds = stepSize)
			yy = int(dt.year)
			mm = int(dt.month)
			dd = int(dt.day)
			hh = int(dt.hour)
			mn = int(dt.minute)
			sc = int(dt.second)
			lat,lon,alt =ECItoLLA(solver.y,yy,mm,dd,hh,mn,sc)
			LLA_sol[k]= [lat,lon,alt]
			timestamp[k]= [yy,mm,dd,hh,mn,sc]
			if (int(now.year) == yy and int(now.month) == mm and int(now.day) == dd and \
				int(now.hour)==hh and int(now.minute) == mn and abs(now.second-sc)<=stepSize):
				ind = k
			#print(solver.y)
			print(lat,lon,alt)
			#print(timestamp[k])
			k=k+1

		self.time_index = ind
		self.ECI = sol[ind:,:]
		self.LLA = LLA_sol[ind:,:]
		self.timestamp = timestamp[ind:,:]
		#track(LLA_sol)
		#animated_track(LLA_sol)
		return self.ECI, self.LLA, self.timestamp

	def run_solver2(self, numOrbits, stepSize=5):
		earthMu = self.earthMu
		bc = self.bc
		state = self.state
		t0 = 0
		yy = self.epoch[0]
		mm = self.epoch[1]
		dd = self.epoch[2]
		hh = self.epoch[3]
		mn = self.epoch[4]
		sc = self.epoch[5]

		stepSize = 3  # seconds
		now = datetime.datetime.utcnow().replace(microsecond=0)
		end = now + timedelta(seconds=numOrbits * self.period)
		epoch = datetime.datetime(yy, mm, dd, hh, mn, sc).replace(microsecond=0)
		diff = (end - epoch)
		t1 = int(diff.total_seconds())

		# Creating solution variables
		t = np.arange(t0, t1, stepSize)
		sol = odeint(diffEqn, state, t, rtol=10**(-11), atol=10**(-11), args=(earthMu,bc))

		dt = datetime.datetime(yy, mm, dd, hh, mn, sc)
		timestamp = np.empty((len(sol), 6))
		timestamp[0] = [yy, mm, dd, hh, mn, sc]

		LLA_sol = np.empty((len(sol), 3))
		lat, lon, alt = ECItoLLA(state, yy, mm, dd, hh, mn, sc)
		LLA_sol[0] = [lat, lon, alt]

		k = 1
		ind = 0
		######################### Delete Later#############################
		fine = open("ISS_ORBIT_2021_7_4.txt", 'w')
		fine.write("|-------Julian Day Number----------|---------ECI Position Vector (km)--------|---------ECI Velocity Vector (km/s)------|\n")
		fine.write("|----------------------------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
		JDN = JulianDate(timestamp[0])
		fine.write(f"     {JDN:.6f}                ,")
		fine.write(f"  {state[0]/1000:.3f},    {state[1]/1000:.3f},    {state[2]/1000:.3f}")
		fine.write(f",         {state[3]/1000:.3f},        {state[4]/1000:.3f},        {state[5]/1000:.3f}\n")
		####################################################################
		while k < len(sol):
			dt = dt + timedelta(seconds=stepSize)
			yy = int(dt.year)
			mm = int(dt.month)
			dd = int(dt.day)
			hh = int(dt.hour)
			mn = int(dt.minute)
			sc = int(dt.second)
			########################## Delete Later###################################
			JDN = JulianDate([yy,mm,dd,hh,mn,sc])
			fine.write(f"     {JDN:.6f}                ,")
			fine.write(f"  {sol[k,0] / 1000:.3f},    {sol[k,1] / 1000:.3f},    {sol[k,2] / 1000:.3f}")
			fine.write(f",         {sol[k,3] / 1000:.3f},        {sol[k,4] / 1000:.3f},        {sol[k,5] / 1000:.3f}\n")
			##########################################################################
			lat, lon, alt = ECItoLLA(sol[k, :], yy, mm, dd, hh, mn, sc)
			LLA_sol[k] = [lat, lon, alt]
			timestamp[k] = [yy, mm, dd, hh, mn, sc]
			if (int(now.year) == yy and int(now.month) == mm and int(now.day) == dd and \
					int(now.hour) == hh and int(now.minute) == mn and abs(now.second - sc) <= stepSize):
				ind = k
			#print(lat, lon, alt)
			k = k + 1
		fine.close()
		self.time_index = ind
		self.ECI = sol[ind:, :]
		self.LLA = LLA_sol[ind:, :]
		self.timestamp = timestamp[ind:, :]

		return self.ECI, self.LLA, self.timestamp

