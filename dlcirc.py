import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.special as sc
import numpy as np
from special_functions import besselj

tdel = 50
tmax = 10000
step = 0.01
ndiv = 50
lado = 1.
vel = 1.
pos = [lado/2., 2*np.pi/3]

tempo = np.arange(0, tmax, step)

t = 0

rvec = np.linspace(0, lado, ndiv)
phivec = np.linspace(0, 2*np.pi, ndiv)

theta = np.pi / 4

def aharmonic(A, B, C, D, m, n):
    return (A*np.cos(vel*sc.jn_zeros(m,n)[n-1]*tempo/lado)+B*np.sin(vel*sc.jn_zeros(m,n)[n-1]*tempo/lado))*(besselj(m, sc.jn_zeros(m,n)[n-1]*pos[0]/lado, 1)*(C*np.cos(m*pos[1])+D*np.sin(m*pos[1]))*sc.jn_zeros(m,n)[n-1]*np.cos(pos[1])/lado
        +sc.jn(m, sc.jn_zeros(m,n)[n-1]*pos[0]/lado)*(C*m*np.sin(pos[1])*np.sin(m*pos[1])/pos[0] - D*m*np.sin(pos[1])*np.cos(m*pos[1])/pos[0]))

def bharmonic(A, B, C, D, m, n):
    return (A*np.cos(vel*sc.jn_zeros(m,n)[n-1]*tempo/lado)+B*np.sin(vel*sc.jn_zeros(m,n)[n-1]*tempo/lado))*(besselj(m, sc.jn_zeros(m,n)[n-1]*pos[0]/lado, 1)*(C*np.cos(m*pos[1])+D*np.sin(m*pos[1]))*sc.jn_zeros(m,n)[n-1]*np.sin(pos[1])/lado
        +sc.jn(m, sc.jn_zeros(m,n)[n-1]*pos[0]/lado)*(-C*m*np.cos(pos[1])*np.sin(m*pos[1])/pos[0] + D*m*np.cos(pos[1])*np.cos(m*pos[1])/pos[0]))

def y(a, b, th):
    return 2*np.cos(b)*np.cos(b)*np.cos(th+a)*np.cos(a)/(1-np.sin(a)*np.sin(a)*np.sin(b)*np.sin(b)) - np.cos(th)

def x(a, b, th):
    return -2*np.sin(b)*np.cos(b)*np.cos(th+a)*np.cos(a)/(1-np.sin(a)*np.sin(a)*np.sin(b)*np.sin(b))

alpha = np.arctan((aharmonic(1., 1., 1., 1., 0, 1)) + aharmonic(1., 1., 1., 1., 1, 1))
beta = np.arctan((bharmonic(1., 1., 1., 1., 0, 1)) + bharmonic(1., 1., 1., 1., 1, 1))

xplot = x(alpha, beta, theta)
yplot = y(alpha, beta, theta)

# def f(a, b, c, m, n, g, h, s):
#     return np.array([aharmonic(a, b, m, n, c, i, j, lado, s)/np.pi for i in g for j in h]).reshape(ndiv,ndiv)

# plt.xlim(-1.2,1.2)
# plt.ylim(-1.2,1.2)

# plt.scatter(xplot, yplot, color = 'red', s = 20)

fig2, ax2 = plt.subplots()

# ims = []
# for i in range(5000):
#     im = ax2.scatter(xplot[i], yplot[i], color = 'red', s = 20)
#     ims.append([im])

line, = ax2.plot([], [], lw = 2, color = 'red')

ax2.set_xlim(-1.2,1.2)
ax2.set_ylim(-1.2,1.2)

def init():
    line.set_data([], [])
    return line,

def animate(i):
    if(i<tdel):
        line.set_data(xplot[:i], yplot[:i])
    else :
        line.set_data(xplot[i-tdel:i], yplot[i-tdel:i])
    return line,

ani = animation.FuncAnimation(fig2, animate, init_func = init, frames= 7000, interval = 1, blit = True)
# ani.save('growingCoil.gif', writer = 'ffmpeg', fps = 1000)

# fig, ax = plt.subplots(2,2)

# ax[0,0].scatter(tempo, alpha, color = 'black', s = 20)
# ax[0,0].set_title('alpha')
# ax[0,1].scatter(tempo, beta, color = 'blue', s = 20)
# ax[0,1].set_title('beta')
# ax[1,0].scatter(tempo, xplot, color = 'black', s = 20)
# ax[1,0].set_title('x')
# ax[1,1].scatter(tempo, yplot, color = 'blue', s = 20)
# ax[1,1].set_title('y')

plt.show()