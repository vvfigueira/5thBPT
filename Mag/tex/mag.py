import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative
import matplotlib.animation as animation

theta = 0
phi = 0
def theta1(temp): 
    return np.sin(temp)
def phi1(temp):
    return temp*np.pi/6
omega = np.pi/2

r = np.linspace(0.1,1,1000)
t = np.linspace(0.1,1,100)
M = 1
g = 10

mu0 = 1

m = 1
m1 = 1

def U1(dist):
    return mu0*m*m1*np.sin(theta1)*np.cos(phi1)/(4*np.pi*(dist**3)) - dist*M*g

def U(dist, temp):
    return (-mu0*m*m1*(3*np.sin(theta1(temp))*((np.sin(theta))**2)*np.cos(phi1(temp)-phi)*np.cos(phi-omega*temp)
                    +3*np.cos(theta1(temp))*np.cos(theta)*np.sin(theta1(temp))*np.cos(phi-omega*temp)
                    -np.sin(theta1(temp))*np.cos(phi1(temp)-omega*temp))/(4*np.pi*(dist**3)) 
            -dist*M*g*np.cos(theta) + dist**2 *(derivative(theta1,temp))**2
            +dist**2 *np.sin(theta1(temp))*(derivative(phi1, temp))**2)            

# plt.plot(r,U(r, 1))
# plt.ylim(-10,10)
# fig, ax = plt.subplots()

# ax.plot(r,derivative(U,r,0.01))

fig2, ax2 = plt.subplots()

line, = ax2.plot([], [], lw = 2, color = 'red')

ax2.set_xlim(0,1)
ax2.set_ylim(-100,100)

def init():
    line.set_data([], [])
    return line,

def animate(i):
    line.set_data(r, U(r, i/100))
    return line,

ani = animation.FuncAnimation(fig2, animate, init_func = init, frames= 7000, interval = 1, blit = True)

plt.show()