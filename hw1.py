import math
import matplotlib.pyplot as plt
import numpy as np

'''
# section a
x = np.linspace(0, 2*math.pi, 100)
y = np.sin(x)
y_t = x-x**3/6+x**5/120
plt.plot(x, y)
plt.plot(x, y_t)
plt.show()
'''
'''
# section b
x = [0, 0.5*math.pi, math.pi, 3/2*math.pi, 2*math.pi]
y = np.cos(x)
y_t = 1-np.power(x, 2)/2+np.power(x, 4)/24
print(y)
print(y_t)

x = np.linspace(0, 2*math.pi, 100)
y = np.cos(x)
y_t = 1-x**2+x**4/24
plt.plot(x, y)
plt.plot(x, y_t)
plt.show()
'''
# section d
x = [-2*math.pi, -1.5*math.pi, -1*math.pi, -0.5*math.pi, 0, 0.5*math.pi, math.pi, 1.5*math.pi, 2*math.pi]
y = 1/np.cos(x)
y_t = 1+np.power(x, 2)/2+5*np.power(x, 4)/24
print(y)
print(y_t)

x = np.linspace(-2*math.pi, 2*math.pi, 100)
y = 1/np.cos(x)
y_t = 1+np.power(x, 2)/2+5*np.power(x, 4)/24
plt.plot(x, y)
plt.plot(x, y_t)
plt.show()