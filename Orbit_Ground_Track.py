from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from tkinter import *
import tkinter as tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
import threading
from Satellites import Satellite
from Ground_Track import *
from TLE_parser import *
from Coord_transform import *
import datetime


class mainWindow(Tk):

  LLA_sol = []
  ECI_sol = []
  TS_sol= []
  TLE_lines = []
  a = 6378137.0
  b = 6356752.314
  n = 0
  earthMu = 3.986004418 * (10**(14))
  yy=0
  mm=0
  dd=0
  hh=0
  mn=0
  sc=0
  def __init__(self):
    
    self.i=-1
    super(mainWindow,self).__init__()
    #self.title("Orbit Ground Track");

    self.matplotCanvas()
    self.logo = tk.Label(self,text="Orbit Ground Track", font = "Courier 18 bold")
    self.logo.place(x=420,y=5)
    
    # TLE group box
    self.LB = LabelFrame(self,text="TLE", font ="Courier 15 bold",width =620, height=230)
    self.LB.place(x=0,y=570)

    self.line1_label = Label(self,text="Line 1: ", font ="Courier 12 bold")
    self.line1_label.place(x=20, y = 600)

    self.line1_entry = Entry(self, font ="Courier 12", width = 40)
    self.line1_entry.place(x=120, y = 600)

    self.line2_label = Label(self,text="Line 2: ", font ="Courier 12 bold")
    self.line2_label.place(x=20, y = 650)

    self.line2_entry = Entry(self, font ="Courier 12", width = 40)
    self.line2_entry.place(x=120, y = 650)

    self.line3_label = Label(self,text="Line 3: ", font ="Courier 12 bold")
    self.line3_label.place(x=20, y = 700)
    
    self.line3_entry = Entry(self, font ="Courier 12", width = 40)
    self.line3_entry.place(x=120, y = 700)

    self.read_Button = Button(self, text = "READ", font ="Courier 12 bold", height = 1, width =5, command=lambda:threading.Thread(target=self.click_readButton).start())
    self.read_Button.place(x = 530, y = 750)
    
    # Satellite group box
    self.RB = LabelFrame(self,text="Satellite", font ="Courier 15 bold",width =530, height=230)
    self.RB.place(x=620,y=570)

    self.name_label = Label(self,text="Name: ", font ="Courier 12 bold")
    self.name_label.place(x=650, y = 620)

    self.catno_label = Label(self,text="Catalog no. ", font ="Courier 12 bold")
    self.catno_label.place(x=650, y = 700)

    self.class_label = Label(self,text="Classification: ", font ="Courier 12 bold")
    self.class_label.place(x=650, y = 770)

    self.ed_label = Label(self,text="Epoch Date: ", font ="Courier 12 bold")
    self.ed_label.place(x=900, y = 620)

    self.ep_label = Label(self,text="Epoch Time: ", font ="Courier 12 bold")
    self.ep_label.place(x=900, y = 700)


    # Live group box
    self.RV = LabelFrame(self,text="Live", font ="Courier 15 bold",width =390, height=800)
    self.RV.place(x=1150,y=0)

    factor_y = -150
    factor_x = -1970
    self.Num_Orbits_label = Label(self,text="Number of Orbits: ", font ="Courier 12 bold")
    self.Num_Orbits_label.place(x=3170+factor_x, y = 50)

    self.Orbit_entry_label = Entry(self, font ="Courier 12", width =2)
    self.Orbit_entry_label.place(x=3350+factor_x, y = 50)

    self.run_button = Button(self, text ="RUN", font ="Courier 12", width =5, command = self.click_runButton)
    #lambda:threading.Thread(target=self.click_runButton).start()
    self.run_button.place(x=3400+factor_x, y= 45)

    self.lat_label = Label(self,text=" Lat: ", font ="Courier 12 bold")
    self.lat_label.place(x=3170+factor_x, y = 250+factor_y)
    self.latres_label = Label(self,text="", font ="Courier 12 bold")
    self.latres_label.place(x=3250+factor_x, y = 250+factor_y)

    self.lon_label = Label(self,text=" Lon: ", font ="Courier 12 bold")
    self.lon_label.place(x=3170+factor_x, y = 300+factor_y)
    self.lonres_label = Label(self,text="", font ="Courier 12 bold")
    self.lonres_label.place(x=3250+factor_x, y = 300+factor_y)

    self.alt_label = Label(self,text=" Alt: ", font ="Courier 12 bold")
    self.alt_label.place(x=3170+factor_x, y = 350+factor_y)
    self.altres_label = Label(self,text="", font ="Courier 12 bold")
    self.altres_label.place(x=3250+factor_x, y = 350+factor_y)

    self.vel_label = Label(self,text=" Vel: ", font ="Courier 12 bold")
    self.vel_label.place(x=3170+factor_x, y = 400+factor_y)
    self.velres_label = Label(self,text="", font ="Courier 12 bold")
    self.velres_label.place(x=3250+factor_x, y = 400+factor_y)


    self.period_label = Label(self,text=" Period: ", font ="Courier 12 bold")
    self.period_label.place(x=3170+factor_x, y = 450+factor_y)
    self.periodres_label = Label(self,text="", font ="Courier 12 bold")
    self.periodres_label.place(x=3250+factor_x, y = 450+factor_y)

    self.inc_label = Label(self,text=" Inc: ", font ="Courier 12 bold")
    self.inc_label.place(x=3170+factor_x, y = 500+factor_y)
    self.incres_label = Label(self,text="", font ="Courier 12 bold")
    self.incres_label.place(x=3250+factor_x, y = 500+factor_y)


    self.ecc_label = Label(self,text=" Ecc: ", font ="Courier 12 bold")
    self.ecc_label.place(x=3170+factor_x, y = 550+factor_y)
    self.eccres_label = Label(self,text="", font ="Courier 12 bold")
    self.eccres_label.place(x=3250+factor_x, y = 550+factor_y)

    self.ARGP_label = Label(self,text=" ARGP: ", font ="Courier 12 bold")
    self.ARGP_label.place(x=3170+factor_x, y = 600+factor_y)
    self.ARGPres_label = Label(self,text="", font ="Courier 12 bold")
    self.ARGPres_label.place(x=3250+factor_x, y = 600+factor_y)


    self.RAAN_label = Label(self,text=" RAAN: ", font ="Courier 12 bold")
    self.RAAN_label.place(x=3170+factor_x, y = 650+factor_y)
    self.RAANres_label = Label(self,text="", font ="Courier 12 bold")
    self.RAANres_label.place(x=3250+factor_x, y = 650+factor_y)


    self.ta_label = Label(self,text=" True Anomaly: ", font ="Courier 12 bold")
    self.ta_label.place(x=3170+factor_x, y = 700+factor_y)
    self.tares_label = Label(self,text="", font ="Courier 12 bold")
    self.tares_label.place(x=3350+factor_x, y = 700+factor_y)

    self.save_button = Button(self, text ="SAVE", font ="Courier 12 bold", width =5, command=lambda:threading.Thread(target=self.click_saveButton).start())
    self.save_button.place(x=3250+factor_x, y= 800+factor_y)

    self.run_button = Label(self, text ="Created by Abinay Brown", font ="Courier 12 bold")
    self.run_button.place(x=3250+factor_x, y= 900+factor_y)
    


    
    


  
  def matplotCanvas(self):
    self.fig=plt.figure(figsize = (12, 7))
    self.m = Basemap(projection='gall')

    #self.m.etopo(scale=0.1, alpha=0.8)

    self.m.drawcoastlines()
    self.m.fillcontinents(color='lightgreen',lake_color='lightblue')
    self.m.drawmapboundary(fill_color='lightblue')
    self.m.drawparallels(np.arange(-90.,91.,30.), labels=[True, False, False, False])
    self.m.drawmeridians(np.arange(-180.,181.,30.), labels=[0,0,0,1])
    
    
    self.canvas = FigureCanvasTkAgg(self.fig,self)
    self.canvas.get_tk_widget().place(x=-50,y=-70)
    self.canvas.draw()



  def click_runButton(self):
    
    ISS = Satellite(self.TLE_lines[0],self.TLE_lines[1],self.TLE_lines[2])
    num = int(self.Orbit_entry_label.get())
    ECI, LLA, TS = ISS.run_solver2(num)
    self.LLA_sol = LLA
    self.ECI_sol = ECI
    self.TS_sol= TS
    self.n =0
    self.n = int(self.n)
    self.yy = TS[self.n,0]
    self.mm = TS[self.n,1]
    self.dd = TS[self.n,2]
    self.hh = TS[self.n,3]
    self.mn = TS[self.n,4]
    self.sc = TS[self.n,5]
    #x,y = self.m(0, 0)
    #self.point = self.m.plot(x, y, 'ro', markersize=5)[0]
    #self.point.set_data([], [])
    self.m.scatter(self.LLA_sol[:,1], self.LLA_sol[:,0], latlon= True, s= 2, c='b', zorder=7)
    self.canvas.draw()
    x,y = self.m(self.LLA_sol[0,1], self.LLA_sol[0,0])
    point = self.m.plot(x, y, 'ro', markersize=5,zorder=8)[0]

    def init():
        point.set_data([], [])
        return point,

    # animation function.  This is called sequentially
    
    def animate(frame):
      now = datetime.datetime.utcnow().replace(microsecond=0)
      
      while not (int(now.year) == self.yy and int(now.month) == self.mm and int(now.day) == self.dd and \
          int(now.hour)==self.hh and int(now.minute) == self.mn and abs(now.second-self.sc)<=2):
          self.yy = TS[self.n,0]
          self.mm = TS[self.n,1]
          self.dd = TS[self.n,2]
          self.hh = TS[self.n,3]
          self.mn = TS[self.n,4]
          self.sc = TS[self.n,5]
          self.n=self.n+1
      
      self.i=self.i+1
      lon = self.LLA_sol[self.n][1]
      lat = self.LLA_sol[self.n][0]
      print(self.i)
      self.update_labels(self.n)
      x, y = self.m(lon, lat)
      point.set_data(x, y)

      return point,

    a = animation.FuncAnimation(self.fig, animate, frames=len(self.LLA_sol[:,0])-1, repeat=False, interval=3000)
    self.canvas.draw()

    return 0

  def click_readButton(self):
    self.TLE_lines = []
    self.TLE_lines.append(self.line1_entry.get())
    self.TLE_lines.append(self.line2_entry.get())
    self.TLE_lines.append(self.line3_entry.get())
    x,y,z,velx,vely,velz,bc,yy,mm,dd,hh,mn,sc,t,name,cn,clss= read_TLE(self.TLE_lines);
    self.name_label['text'] = "Name: " + name
    self.catno_label['text'] = "Catalog no. " + cn
    self.class_label['text'] = "Classification: " + clss
    self.ed_label['text'] = "Epoch Date: " + str(mm)+'/'+str(dd)+'/'+str(yy) 
    self.ep_label['text'] = "Epoch Time: " + str(hh)+':'+str(mn)+':'+str(sc)

    return 0

  def click_saveButton(self):
    return 0

  def update_labels(self,ind):
      lat = np.deg2rad(round(self.LLA_sol[ind,0],6))
      f1 = ((self.a**2)*(np.cos(lat)))**2
      f2 = ((self.b**2)*(np.sin(lat)))**2
      f3 = (self.a*np.cos(lat))**2
      f4 = (self.b*np.sin(lat))**2
      radius = np.sqrt((f1+f2)/(f3+f4))
      alt =round((self.LLA_sol[ind,2]-radius)/1000,3)
      vel = round(np.linalg.norm(self.ECI_sol[ind,3:])/1000,4)
      a, e, i, argp, W, v = ECItoOE(self.ECI_sol[ind,:],self.earthMu)
      period = round(2*np.pi*np.sqrt((a**3)/self.earthMu)/60,3)

      self.latres_label['text'] = str(round(self.LLA_sol[ind,0],4))
      self.lonres_label['text'] = str(round(self.LLA_sol[ind,1],4))
      self.altres_label['text'] = str(alt)
      self.velres_label['text'] = str(vel)
      self.periodres_label['text'] = str(period)
      self.incres_label['text'] = str(round(i,2))
      self.eccres_label['text'] = str(round(e,6))
      self.ARGPres_label['text'] = str(round(argp,2))
      self.RAANres_label['text'] = str(round(W,2))
      self.tares_label['text'] = str(round(v,2))
      return 0
    
  
  
    

    

if __name__=='__main__':
  root =mainWindow()
  root.mainloop()