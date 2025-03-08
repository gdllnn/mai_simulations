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


class SpaceShip:
    def __init__(self,X0,Y0,Vx0,Vy0,Fmax,R,Co,Cz):
        self.X0=float(X0)
        self.Y0=float(Y0)
        self.Vx0=float(Vx0)
        self.Vy0=float(Vy0)
        self.X = float(X0)
        self.Y = float(Y0)
        self.Vx = float(Vx0)
        self.Vy = float(Vy0)
        self.R=R
        self.Co=Co
        self.Cz=Cz
        self.Fmax = float(Fmax)
        self.TraceX=np.array([self.X])
        self.TraceY=np.array([self.Y])
        self.Ship_X = self.R*np.array([-0.4, 0, 0.5, 1, 0.5, 0, -0.4, -0.7, -0.15, 0, -0.4, -0.4, -0.7, -0.15, 0])
        self.Ship_Y = self.R*np.array([0.2, 0.3, 0.25, 0, -0.25, -0.3, -0.2, -0.5, -0.5, -0.3, -0.2, 0.2, 0.5, 0.5, 0.3])
        self.Flame_X = self.R*np.array([0, -0.4, -0.3, -0.8, -0.7, -1, -0.7, -0.8, -0.3, -0.4, 0])
        self.Flame_Y = self.R*np.array([0.2, 0.25, 0.2, 0.18, 0.1, 0, -0.1, -0.18, -0.2, -0.25, -0.2])

    def DrawShip(self,ax):
        self.Dfl = ax.plot(self.X0 + self.Ship_X[0] * self.R + self.Flame_X * self.R, self.Y0 + self.Flame_Y * self.R,
                           color=self.Cz)[0]
        self.Dsh = ax.plot(self.X0+self.Ship_X*self.R,self.Y0+self.Ship_Y*self.R,color=self.Co)[0]
        self.Dt = ax.plot(self.TraceX,self.TraceY,':',color=self.Co)[0]

    def RefreshShip(self,X,Y,Vx,Vy,Phi,Fdv):
        self.X = X
        self.Y = Y
        self.Vx = Vx
        self.Vy = Vy
        self.TraceX = np.append(self.TraceX,X)
        self.TraceY = np.append(self.TraceY,Y)
        RShip_X, RShip_Y = Rot2D(self.Ship_X, self.Ship_Y, np.pi/2 - Phi)
        RFlame_X, RFlame_Y = Rot2D(self.Flame_X*Fdv+self.Ship_X[0], self.Flame_Y, np.pi/2 - Phi)

        self.Dsh.set_data(self.X + RShip_X, self.Y + RShip_Y)
        self.Dfl.set_data(self.X + RFlame_X, self.Y + RFlame_Y)
        self.Dt.set_data(self.TraceX, self.TraceY)

class PlanetSystem:
    def __init__(self,Planets,Ship='Null'):
        self.Planets=Planets
        self.N = len(self.Planets)
        self.X = np.array([self.Planets[i].X for i in range(self.N)],'float')
        self.Y = np.array([self.Planets[i].Y for i in range(self.N)],'float')
        self.VX = np.array([self.Planets[i].Vx for i in range(self.N)],'float')
        self.VY = np.array([self.Planets[i].Vy for i in range(self.N)],'float')
        self.Ship = Ship
        if self.Ship!='Null':
            self.ShX = float(Ship.X)
            self.ShY = float(Ship.Y)
            self.ShVX = float(Ship.Vx)
            self.ShVY = float(Ship.Vy)
            self.ShPhi = 0.0
            self.ShFdv = 0.0
        self.Camera = 0

    def SetGravity(self,G):
        self.G = G
    def DrawSystem(self,ax):
        self.ax = ax
        self.ax.cla()
        #self.ax.plot()
        for planet in self.Planets:
            planet.DrawPlanet(self.ax)
        if self.Ship!='Null':
            self.Ship.DrawShip(self.ax)


    def RefreshSystem(self):
        if self.Camera == 1:
            self.ax.set(xlim=self.Xlim+self.ShX, ylim=self.Ylim+self.ShY)
        for i in range(self.N):
            self.Planets[i].RefreshPlanet(self.X[i],self.Y[i],self.VX[i],self.VY[i])
        if self.Ship != 'Null':
            self.Ship.RefreshShip(self.ShX, self.ShY, self.ShVX, self.ShVY, self.ShPhi, self.ShFdv)
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
        if self.Ship != 'Null':
            ShX = sp.symbols('Shx')
            ShY = sp.symbols('Shy')
            ShVX = sp.symbols('ShVx')
            ShVY = sp.symbols('ShVy')
            Phi = sp.symbols('phi')
            Fdv = sp.symbols('fdv')
            ShDX = ShVX
            ShDY = ShVY
            ShDVX = sum([self.Planets[j].m * (X[j] - ShX) * self.G / ((X[j] - ShX) ** 2 + (Y[j] - ShY) ** 2) ** 1.5
                        for j in range(self.N)]) + self.Ship.Fmax*Fdv*sp.sin(Phi)
            ShDVY = sum([self.Planets[j].m * (Y[j] - ShY) * self.G / ((X[j] - ShX) ** 2 + (Y[j] - ShY) ** 2) ** 1.5
                        for j in range(self.N)]) + self.Ship.Fmax*Fdv*sp.cos(Phi)

            self.ShipMoveEquations = sp.lambdify([X, Y, VX, VY, ShX, ShY, ShVX, ShVY, Phi, Fdv], [ShDX, ShDY, ShDVX, ShDVY], 'numpy')

Planet1 = Planet(10,0,0,1.5,40,5,[1,0,0],[0,1,0])
Planet2 = Planet(-10,0,0,-2.1,30,5,[0,1,0],[0,0,1])
Planet3 = Planet(0,0,0,-0.33,30,5,[0,0,1],[1,0,0])
Ship = SpaceShip(5,0,0,1,100,1,[1,1,1],[1,0,0])
PS = PlanetSystem([Planet1,Planet2,Planet3],Ship)
PS.SetGravity(10)
fimeName = "Universe2.universe"
dict = {'PS': PS}
with open(fimeName,'wb') as file:
    pickle.dump(dict, file)

print("Сделано")