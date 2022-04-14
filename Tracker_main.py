
from scipy.integrate import ode
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

# Earth Constants 
earthMu = 3.986004418 * (10**(14)) 	# m^3/s^2

if __name__ == "__main__":
	
	# setting up initial conditions
	TLE_lines=[]
	TLE_lines.append("ISS (ZARYA)");
	TLE_lines.append("1 25544U 98067A   20360.82210118  .00001438  00000-0  33988-4 0  9995");
	TLE_lines.append("2 25544  51.6447 120.2583 0001352 157.5694 338.5516 15.49226385261722");
	t0=0
	x,y,z,velx,vely,velz,bc,yy,mm,dd,hh,mn,sc,t,name= read_TLE(TLE_lines);
	state = [x, y, z, velx, vely, velz]
	#print(read_TLE(TLE_lines))
	# setting up the ODE solver
	solver = ode(diffEqn)
	solver.set_integrator('dop853', atol = 10**(-9), rtol = 10**(-9)) 
	solver.set_f_params(earthMu, bc)
	solver.set_initial_value(state, t0)
	
	
	stepSize = 5 # seconds
	now = datetime.datetime.utcnow() #+ timedelta(seconds=300)
	epoch = datetime.datetime(yy,mm,dd,hh,mn,sc)
	diff = now-epoch
	t1 = diff.total_seconds()  # 300 minutes in seconds
	
	t = np.arange(t0,t1, stepSize)
	sol = np.empty((len(t), 6))
	sol[0]= state
	k = 1
	LLA_sol = np.empty((len(t)-1, 3))
	dt = datetime.datetime(yy,mm,dd,hh,mn,sc)
	l=0
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
		#print(solver.y)
		LLA_sol[l]= [lat,lon,alt] 
		print(lat,lon,alt)
		l = l+1
		k=k+1
	#plt.plot(LLA_sol[:,1],LLA_sol[:,0])
	#plt.plot(t,sol[:,1])
	#plt.show()
	#RunSimulation(sol)

	track(LLA_sol)
	#animated_track(LLA_sol)
	