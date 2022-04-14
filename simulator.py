'''
vpython trajectory simulation for Aerobraking.
'''
from vpython import *

def RunSimulation(state):

	# Mars 
	mars = sphere(pos = vector(0,0,0), radius = 1000, color = color.red)

	# Spacecraft
	x = state[0, 0]
	y = state[0, 1]
	z = state[0, 2]
	

	spacecraft = sphere(pos = vector(x/1000,y/1000,z/1000), radius = 75, color = color.white, make_trail=True, trail_type="points",
              )
	
	i = 0
	while (i < len(state)):
		x = state[i, 0] /1000
		y = state[i, 1] /1000
		z = state[i, 2] /1000
		spacecraft.pos =  vector(x, y, z)
		rate(1000)
		i = i+1
	return
