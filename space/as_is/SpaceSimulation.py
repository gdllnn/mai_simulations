from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
import SpaceWindow
import space_widget
import math
import numpy as np
import random as rd
from matplotlib.animation import FuncAnimation
import sympy as sp
import pprint
import time
import scipy.io as io
import pickle
#from matplotlib.backends.backend_qt5agg import FigureCanvas

from matplotlib.figure import Figure


def Rot2D(X, Y, Alpha):
    RX = X * np.cos(Alpha) - Y * np.sin(Alpha)
    RY = X * np.sin(Alpha) + Y * np.cos(Alpha)
    return RX, RY


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
        print('--------------------------')
        print(self.TraceX, len(self.TraceX))
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
        #self.ax.cla()
        #self.ax.plot()
        for planet in self.Planets:
            planet.DrawPlanet(self.ax)
        if self.Ship!='Null':
            self.Ship.DrawShip(self.ax)


    def RefreshSystem(self):
        if self.Camera == 1:
            self.ax.axis('scaled')
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

def DrawTheSpace(axes):
    print(2)
    axes.fill([-100, 100, 100, - 100], [-100, - 100, 100, 100], 'black')
    nstars = 5000
    xstars = 200 * np.random.random(nstars) - 100
    ystars = 200 * np.random.random(nstars) - 100
    brstars = 0.5 + np.random.random(nstars) / 2
    sizestars = np.random.random(nstars)
    for i in range(nstars):
        axes.plot(xstars[i], ystars[i], marker='o', markersize=sizestars[i], color=[brstars[i], brstars[i], brstars[i]])


class SpaceWidget(QMainWindow, SpaceWindow.Ui_MainWindow):

    def __init__(self):
        QMainWindow.__init__(self)

        self.setupUi(self)
        global Simulation
        Simulation = 0
        self.setWindowTitle("Просторы галактики")

        self.StartButton.clicked.connect(self.HereAreWeGo)


    def HereAreWeGo(self):
        global Simulation, PS, multik
        #PS = 0
        #multik = 0
        if Simulation == 0:
            #Simulation = 1
            Filename = self.File_Edit.text()+'.universe'

            dt = float(self.Step_Edit.text())
            Xmin = float(self.Xmin_Edit.text())
            Xmax = float(self.Xmax_Edit.text())
            Ymin = float(self.Ymin_Edit.text())
            Ymax = float(self.Ymax_Edit.text())


            #self.spacewidget.canvas.axes.axis('manual')
            self.spacewidget.canvas.axes.clear()

            self.spacewidget.canvas.axes.axis('scaled')
            self.spacewidget.canvas.axes.set(xlim=[Xmin, Xmax], ylim=[Xmin*12/22.3, Xmax*12/22.3])
            #self.spacewidget.canvas.axes.axis('equal')
            print(1)
            DrawTheSpace(self.spacewidget.canvas.axes)
            print(3)
            #print(self.spacewidget.canvas.axes.lims)

            with open(Filename, 'rb') as file:
                PS = pickle.load(file)['PS']
            PS.Xlim = [Xmin, Xmax]
            PS.Ylim = [Xmin*12/22.3, Xmax*12/22.3]
            if self.Camera_Ship.isChecked():
                PS.Camera = 1
            PS.DrawSystem(self.spacewidget.canvas.axes)
            self.spacewidget.canvas.show()
            PS.GetMoveEquations()

            def kadr(i):

                DX, DY, DVX, DVY = PS.MoveEquations(PS.X, PS.Y, PS.VX, PS.VY)
                XK1, YK1, VXK1, VYK1 = np.array(DX), np.array(DY), np.array(DVX), np.array(DVY)
                DX, DY, DVX, DVY = PS.MoveEquations(PS.X + dt / 2 * XK1, PS.Y + dt / 2 * YK1, PS.VX + dt / 2 * VXK1,
                                                    PS.VY + dt / 2 * VYK1)
                XK2, YK2, VXK2, VYK2 = np.array(DX), np.array(DY), np.array(DVX), np.array(DVY)
                DX, DY, DVX, DVY = PS.MoveEquations(PS.X + dt / 2 * XK2, PS.Y + dt / 2 * YK2, PS.VX + dt / 2 * VXK2,
                                                    PS.VY + dt / 2 * VYK2)
                XK3, YK3, VXK3, VYK3 = np.array(DX), np.array(DY), np.array(DVX), np.array(DVY)
                DX, DY, DVX, DVY = PS.MoveEquations(PS.X + dt * XK3, PS.Y + dt * YK3, PS.VX + dt * VXK3, PS.VY + dt * VYK3)
                XK4, YK4, VXK4, VYK4 = np.array(DX), np.array(DY), np.array(DVX), np.array(DVY)

                if PS.Ship != 'Null':
                    PS.ShPhi = self.Course_dial.value()/100
                    PS.ShFdv = self.Force_Slider.value()/100
                    DShX, DShY, DShVX, DShVY = PS.ShipMoveEquations(PS.X, PS.Y, PS.VX, PS.VY,
                                                                    PS.ShX, PS.ShY, PS.ShVX, PS.ShVY,
                                                                    PS.ShPhi, PS.ShFdv)
                    ShXK1, ShYK1, ShVXK1, ShVYK1 = np.array(DShX), np.array(DShY), np.array(DShVX), np.array(DShVY)
                    DShX, DShY, DShVX, DShVY = PS.ShipMoveEquations(PS.X + dt / 2 * XK1, PS.Y + dt / 2 * YK1,
                                                                    PS.VX + dt / 2 * VXK1, PS.VY + dt / 2 * VYK1,
                                                                    PS.ShX + dt / 2 * ShXK1, PS.ShY + dt / 2 * ShYK1,
                                                                    PS.ShVX + dt / 2 * ShVXK1, PS.ShVY + dt / 2 * ShVYK1,
                                                                    PS.ShPhi, PS.ShFdv)
                    ShXK2, ShYK2, ShVXK2, ShVYK2 = np.array(DShX), np.array(DShY), np.array(DShVX), np.array(DShVY)
                    DShX, DShY, DShVX, DShVY = PS.ShipMoveEquations(PS.X + dt / 2 * XK2, PS.Y + dt / 2 * YK2,
                                                                    PS.VX + dt / 2 * VXK2, PS.VY + dt / 2 * VYK2,
                                                                    PS.ShX + dt / 2 * ShXK2, PS.ShY + dt / 2 * ShYK2,
                                                                    PS.ShVX + dt / 2 * ShVXK2, PS.ShVY + dt / 2 * ShVYK2,
                                                                    PS.ShPhi, PS.ShFdv)
                    ShXK3, ShYK3, ShVXK3, ShVYK3 = np.array(DShX), np.array(DShY), np.array(DShVX), np.array(DShVY)
                    DShX, DShY, DShVX, DShVY = PS.ShipMoveEquations(PS.X + dt * XK3, PS.Y + dt * YK3,
                                                                    PS.VX + dt * VXK3, PS.VY + dt * VYK3,
                                                                    PS.ShX + dt * ShXK3, PS.ShY + dt * ShYK3,
                                                                    PS.ShVX + dt * ShVXK3, PS.ShVY + dt * ShVYK3,
                                                                    PS.ShPhi, PS.ShFdv)
                    ShXK4, ShYK4, ShVXK4, ShVYK4 = np.array(DShX), np.array(DShY), np.array(DShVX), np.array(DShVY)

                    PS.ShX += (ShXK1 + 2 * ShXK2 + 2 * ShXK3 + ShXK4) * dt / 6
                    PS.ShY += (ShYK1 + 2 * ShYK2 + 2 * ShYK3 + ShYK4) * dt / 6
                    PS.ShVX += (ShVXK1 + 2 * ShVXK2 + 2 * ShVXK3 + ShVXK4) * dt / 6
                    PS.ShVY += (ShVYK1 + 2 * ShVYK2 + 2 * ShVYK3 + ShVYK4) * dt / 6

                PS.X += (XK1 + 2 * XK2 + 2 * XK3 + XK4) * dt / 6
                PS.Y += (YK1 + 2 * YK2 + 2 * YK3 + YK4) * dt / 6
                PS.VX += (VXK1 + 2 * VXK2 + 2 * VXK3 + VXK4) * dt / 6
                PS.VY += (VYK1 + 2 * VYK2 + 2 * VYK3 + VYK4) * dt / 6

                PS.RefreshSystem()
                if PS.Ship != 'Null':
                    return [PS.Planets[i].Dp for i in range(PS.N)] + \
                           [PS.Planets[i].Dt for i in range(PS.N)] + \
                           [PS.Ship.Dsh, PS.Ship.Dfl, PS.Ship.Dt]
                return [PS.Planets[i].Dp for i in range(PS.N)] + \
                       [PS.Planets[i].Dt for i in range(PS.N)]

            multik = FuncAnimation(self.spacewidget.canvas.figure, kadr, interval=dt * 100, blit=True)

            self.spacewidget.canvas.draw()
        else:
            Simulation = 0
            multik.repeat = False
            self.spacewidget.canvas.axes.clear()
            self.spacewidget.canvas.figure.clf()
            self.spacewidget.canvas.draw()
        return






app = QApplication([])
window = SpaceWidget()
window.show()
app.exec_()