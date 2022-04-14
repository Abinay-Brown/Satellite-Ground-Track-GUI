import numpy as np
import matplotlib.pyplot as plt

x=np.linspace(0,2*np.pi,30)
y=np.sin(x)
y_prime=np.cos(x)
plt.plot(x,y,x,y_prime,color='green',marker='*')
plt.xlabel('bobby is too slow')
plt.ylabel('Monica is also too slow')
plt.grid()
plt.show()

