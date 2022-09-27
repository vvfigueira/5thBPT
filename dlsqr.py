import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

tmax = 10000
step = 0.01
ndiv = 50
lado = 1.
vel = 1.
pos = [lado/3., lado/2. + 0.02]

tempo = np.arange(0, tmax, step)

t = 0

xvec = np.linspace(0, lado, ndiv)
yvec = np.linspace(0, lado, ndiv)

theta = np.pi / 4

def aharmonic(A, B, m, n):
    return (A*np.cos(vel*np.pi*np.sqrt((m*m) + (n*n))*tempo/lado)+B*np.sin(vel*np.pi*np.sqrt((m*m) + (n*n))*tempo/lado))*m*np.pi*np.cos(m*np.pi*pos[0]/lado)*np.sin(n*np.pi*pos[1]/lado)/lado

def bharmonic(A, B, m, n):
    return (A*np.cos(vel*np.pi*np.sqrt((m*m) + (n*n))*tempo/lado)+B*np.sin(vel*np.pi*np.sqrt((m*m) + (n*n))*tempo/lado))*n*np.pi*np.sin(m*np.pi*pos[0]/lado)*np.cos(n*np.pi*pos[1]/lado)/lado

def y(a, b, th):
    return 2*np.cos(b)*np.cos(b)*np.cos(th+a)*np.cos(a)/(1-np.sin(a)*np.sin(a)*np.sin(b)*np.sin(b)) - np.cos(th)

def x(a, b, th):
    return -2*np.sin(b)*np.cos(b)*np.cos(th+a)*np.cos(a)/(1-np.sin(a)*np.sin(a)*np.sin(b)*np.sin(b))

alpha = np.arctan((aharmonic(1., 1., 1, 1)) + aharmonic(1., 1., 3, 2))
beta = np.arctan((bharmonic(1., 1., 1, 1)) + bharmonic(1., 1., 3, 2))

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
    if(i<50):
        line.set_data(xplot[:i], yplot[:i])
    else :
        line.set_data(xplot[i-50:i], yplot[i-50:i])
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