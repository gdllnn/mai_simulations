import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def EquationsOfLodka(X,Y,Phi,Vx,Vy,Omega,F_dv,alpha):
    global m,J,a,Cl,Cb,Cvr,Sl,Sb,Rho

    Vl = Vx*np.cos(Phi) + Vy*np.sin(Phi)
    Vb = - Vx * np.sin(Phi) + Vy * np.cos(Phi)

    Fl = Rho*Cl*Sl*Vl*abs(Vl)/2
    Fb = Rho*Cb * Sb * Vb * abs(Vb) / 2
    Ms = Rho*Cvr * Sb * Omega * abs(Omega) / 2
    print('Ms=', Ms,'M_dv=',a*F_dv*np.sin(alpha))

    Wx = (-Fl*np.cos(Phi) + Fb*np.sin(Phi) + F_dv*np.cos(Phi-alpha))/m
    Wy = (-Fl * np.sin(Phi) - Fb * np.cos(Phi) + F_dv * np.sin(Phi - alpha)) / m
    Epsilon = (-Ms + a*F_dv*np.sin(alpha))/J
    return [Vx,Vy,Omega,Wx,Wy,Epsilon]

def Rot2D(X, Y, Alpha):
    RX = X * np.cos(Alpha) - Y * np.sin(Alpha)
    RY = X * np.sin(Alpha) + Y * np.cos(Alpha)
    return RX, RY

global m,J,a,Cl,Cb,Cvr,Sl,Sb,Rho
m = 400

a = 1.5
J = m*(8/3*a)**2/12
h = 0.3
Cl = 0.2
Cb = 0.9
Cvr = 2
Sl = a/3*h
Sb = 2*a*h*2
Rho = 1000

X = 0
Y = 0
Phi = 0
Vx = 0
Vy = 0
Omega = 0

fig = plt.figure(figsize=[15,7])
ax = fig.add_subplot(1,1,1)
ax.axis('equal')
ax.axes.set(xlim=[-70*a, 70*a],ylim=[-70*a,70*a])
L_X = np.array([-a, -a/3, a/3, a, 5/3*a, a, a/3, -a/3, -a, -a])
L_Y = np.array([0.4*a, a/2, a/2, 0.4*a, 0, -0.4*a, -a/2, -a/2, -0.4*a, 0.4*a])
LodkaX,LodkaY = Rot2D(L_X,L_Y,Phi)

#DrawedMore = ax.plot([-1000*a, 1000],LodkaY,color=[1,1,1])[0]
DrawedLodka = ax.plot(X+LodkaX,Y+LodkaY,color=[0,0,1])[0]
plt.draw()

F_dv = 200
alpha = 0.1
dt = 0.05
t=0


def kadr(j):
    global X,Y,Phi,Vx,Vy,Omega,F_dv,alpha,t
    alpha = 0.5*np.sin(0.1*t)
    Vx, Vy, Omega, Wx, Wy, Epsilon = EquationsOfLodka(X,Y,Phi,Vx,Vy,Omega,F_dv,alpha)
    X+=Vx*dt
    Y+=Vy*dt
    Phi+=Omega*dt
    Vx+=Wx*dt
    Vy+=Wy*dt
    Omega+=Epsilon*dt
    t+=dt
    LodkaX, LodkaY = Rot2D(L_X, L_Y, Phi)
    DrawedLodka.set_data(X+LodkaX,Y+LodkaY)
    print('t=',t,'omega=',Omega)

    return [DrawedLodka]


kino = FuncAnimation(fig, kadr, interval=dt*1000, blit=True)

plt.show()



