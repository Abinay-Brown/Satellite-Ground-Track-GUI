from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import time
import datetime
i=-1
def track(sol):
	# setup Lambert Conformal basemap.
	# set resolution=None to skip processing of boundary datasets.
	fig = plt.figure(figsize = (12, 9))
	m = Basemap(projection='gall')
	m.etopo(scale=0.5, alpha=0.8)

	#date = datetime.datetime.utcnow()
	#CS = m.nightshade(date)

	m.drawparallels(np.arange(-90.,91.,30.), labels=[True, False, False, False])
	m.drawmeridians(np.arange(-180.,181.,30.), labels=[0,0,0,1])
	
	m.scatter(sol[:,1], sol[:,0], latlon= True, s= 1, c='b')
	m.scatter(sol[-1,1], sol[-1,0], latlon= True, s= 50, c='r')
	#plt.title('Ground Trace for %s (UTC)' % date.strftime("%d %b %Y %H:%M:%S"))
	plt.show()
	return 0;
def track2(sol):
	# setup Lambert Conformal basemap.
	# set resolution=None to skip processing of boundary datasets.
	fig = plt.figure(figsize = (12, 9))
	m = Basemap(projection='gall')
	m.etopo(scale=0.5, alpha=0.8)

	#date = datetime.datetime.utcnow()
	#CS = m.nightshade(date)

	m.drawparallels(np.arange(-90.,91.,30.), labels=[True, False, False, False])
	m.drawmeridians(np.arange(-180.,181.,30.), labels=[0,0,0,1])
	
	m.scatter(sol[:,1],sol[:,0], latlon= True, s= 1, c='b')
	m.scatter(sol[0,1], sol[0,0], latlon= True, s= 50, c='r')
	#plt.title('Ground Trace for %s (UTC)' % date.strftime("%d %b %Y %H:%M:%S"))
	plt.show()
	return 0;

def animated_track(sol,my_map,fig):
	#my_map = Basemap(projection='gall')
	#my_map.etopo(scale=0.5, alpha=0.8)
	#my_map.drawcoastlines()
	#my_map.drawcountries()
	
	#my_map.drawmapboundary()
	#my_map.drawmeridians(np.arange(0, 360, 30))
	#my_map.drawparallels(np.arange(-90, 90, 30))
	my_map.scatter(sol[:,1], sol[:,0], latlon= True, s= 1, c='b')
	x,y = my_map(sol[0,1], sol[0,0])
	point = my_map.plot(x, y, 'ro', markersize=5)[0]

	def init():
	    point.set_data([], [])
	    return point,

	# animation function.  This is called sequentially
	
	def animate(frame):
		global i
		i=i+1
		lon = sol[i][1]
		lat = sol[i][0]
		print(i)
		x, y = my_map(lon, lat)
		point.set_data(x, y)
		return point,
	    

	# call the animator.  blit=True means only re-draw the parts that have changed.
	anim = animation.FuncAnimation(plt.gcf(), animate, init_func=init,
	                               frames=len(sol[:][0]), interval=1, blit=True)
	#video = anim.to_html5_video()
	#html=display.HTML(video)
	#display.display(html)
	#fig.show()
	plt.show()