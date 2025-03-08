import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pickle

class Planet:
    def __init__(self,X0,Y0,Vx0,Vy0,m,R,Co,Cz):
        self.X0=float(X0)
        self.Y0=float(Y0)
        self.Vx0=float(Vx0)
        self.Vy0=float(Vy0)
        self.X = float(X0)
        self.Y = float(Y0)
        self.Vx = float(Vx0)
        self.Vy = float(Vy0)
        self.m=m
        self.R=R
        self.Co=Co
        self.Cz=Cz
        self.TraceX=np.array([self.X])
        self.TraceY=np.array([self.Y])

    def DrawPlanet(self,ax):
        self.Dp = ax.plot(self.X0,self.Y0,'o',color=self.Co,markersize=self.R)[0]
        self.Dt = ax.plot(self.TraceX,self.TraceY,':',color=self.Co)[0]

    def RefreshPlanet(self,X,Y,Vx,Vy):
        self.X = X
        self.Y = Y
        self.Vx = Vx
        self.Vy = Vy
        self.TraceX = np.append(self.TraceX,X)
        self.TraceY = np.append(self.TraceY,Y)
        self.Dp.set_data(X,Y)
        self.Dt.set_data(self.TraceX, self.TraceY)

class PlanetSystem:
    def __init__(self,Planets,Xlim,Ylim,G=1):
        self.Planets=Planets
        self.N = len(self.Planets)
        self.X = np.array([self.Planets[i].X for i in range(self.N)],'float')
        self.Y = np.array([self.Planets[i].Y for i in range(self.N)],'float')
        self.VX = np.array([self.Planets[i].Vx for i in range(self.N)],'float')
        self.VY = np.array([self.Planets[i].Vy for i in range(self.N)],'float')
        self.Xlim = Xlim
        self.Ylim = Ylim
        self.G=G

    def DrawPlanets(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(1, 1, 1)
        self.ax.axis('equal')
        self.ax.set(xlim=self.Xlim, ylim=self.Ylim)
        for planet in self.Planets:
            planet.DrawPlanet(self.ax)

    def RefreshSystem(self):
        for i in range(self.N):
            self.Planets[i].RefreshPlanet(self.X[i],self.Y[i],self.VX[i],self.VY[i])

    def GetMoveEquations(self):
        X = sp.symbols('x:'+str(self.N))
        Y = sp.symbols('y:' + str(self.N))
        VX = sp.symbols('Vx:' + str(self.N))
        VY = sp.symbols('Vy:' + str(self.N))
        DX = VX
        DY = VY
        DVX = [sum([self.Planets[j].m * (X[j]-X[i]) * self.G / ((X[j]-X[i])**2 + (Y[j]-Y[i])**2)**1.5
                    for j in range(self.N) if j!=i])
               for i in range(self.N)]
        DVY = [sum([self.Planets[j].m * (Y[j] - Y[i]) * self.G / ((X[j] - X[i]) ** 2 + (Y[j] - Y[i]) ** 2) ** 1.5
                    for j in range(self.N) if j != i])
               for i in range(self.N)]

        self.MoveEquations = sp.lambdify([X,Y,VX,VY],[DX,DY,DVX,DVY],'numpy')

f = 2

if f == 1:
    Planet1 = Planet(10,0,0,1.5,40,5,[1,0,0],[0,1,0])
    Planet2 = Planet(-10,0,0,-2.1,30,5,[0,1,0],[0,0,1])
    Planet3 = Planet(0,0,0,-0.33,30,5,[0,0,1],[1,0,0])
    PS = PlanetSystem([Planet1,Planet2,Planet3],[-12,12],[-12,12],10)

    fimeName = "Universe1.universe"
    dict = {'PS': PS}
    with open(fimeName,'wb') as file:
        pickle.dump(dict, file)

    print("Сделано")
if f == 2:
    Filename = "Universe1.universe"
    with open(Filename, 'rb') as file:
        PS = pickle.load(file)['PS']



    PS.DrawPlanets()
    PS.GetMoveEquations()
    dt = 0.001
    def kadr(i):
        DX,DY,DVX,DVY = PS.MoveEquations(PS.X,PS.Y,PS.VX,PS.VY)
        XK1, YK1, VXK1, VYK1 = np.array(DX), np.array(DY), np.array(DVX), np.array(DVY)
        DX, DY, DVX, DVY = PS.MoveEquations(PS.X+dt/2*XK1, PS.Y+dt/2*YK1, PS.VX+dt/2*VXK1, PS.VY+dt/2*VYK1)
        XK2, YK2, VXK2, VYK2 = np.array(DX), np.array(DY), np.array(DVX), np.array(DVY)
        DX, DY, DVX, DVY = PS.MoveEquations(PS.X+dt/2*XK2, PS.Y+dt/2*YK2, PS.VX+dt/2*VXK2,PS.VY+dt/2*VYK2)
        XK3, YK3, VXK3, VYK3 = np.array(DX), np.array(DY), np.array(DVX), np.array(DVY)
        DX, DY, DVX, DVY = PS.MoveEquations(PS.X + dt*XK3, PS.Y + dt*YK3, PS.VX + dt*VXK3, PS.VY + dt*VYK3)
        XK4, YK4, VXK4, VYK4 = np.array(DX), np.array(DY), np.array(DVX), np.array(DVY)

        PS.X += (XK1+2*XK2+2*XK3+XK4) * dt/6
        PS.Y += (YK1+2*YK2+2*YK3+YK4) * dt/6
        PS.VX += (VXK1+2*VXK2+2*VXK3+VXK4) * dt/6
        PS.VY += (VYK1+2*VYK2+2*VYK3+VYK4) * dt/6

        PS.RefreshSystem()

        return [PS.Planets[i].Dp for i in range(PS.N)] + \
               [PS.Planets[i].Dt for i in range(PS.N)]

    multik = FuncAnimation(PS.fig, kadr, interval=dt*100, blit=True)
    plt.show()

