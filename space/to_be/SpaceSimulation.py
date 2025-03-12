import sys
import math
import numpy as np
import sympy as sp
import pickle
from math import pi, sqrt
from matplotlib.animation import FuncAnimation
from matplotlib.figure import Figure
from PyQt5.QtWidgets import QApplication, QMainWindow
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

# Импорт UI
import SpaceWindow
import space_widget
from space_widget import SpaceWidget

# ------------------------
# Дополнительные функции
# ------------------------
def Rot2D(X, Y, Alpha):
    RX = X * np.cos(Alpha) - Y * np.sin(Alpha)
    RY = X * np.sin(Alpha) + Y * np.cos(Alpha)
    return RX, RY

def DrawTheSpace(axes):
    axes.fill([-100, 100, 100, -100], [-100, -100, 100, 100], 'black')
    nstars = 5000
    xstars = 200 * np.random.random(nstars) - 100
    ystars = 200 * np.random.random(nstars) - 100
    brstars = 0.5 + np.random.random(nstars)/2
    sizestars = np.random.random(nstars)
    for i in range(nstars):
        axes.plot(xstars[i], ystars[i], marker='o', markersize=sizestars[i],
                  color=[brstars[i], brstars[i], brstars[i]])

# ------------------------
# Классы симуляции
# ------------------------
class Planet:
    def __init__(self, X0, Y0, Vx0, Vy0, m, R, Co, Cz):
        self.X0 = float(X0)
        self.Y0 = float(Y0)
        self.Vx0 = float(Vx0)
        self.Vy0 = float(Vy0)
        self.X = float(X0)
        self.Y = float(Y0)
        self.Vx = float(Vx0)
        self.Vy = float(Vy0)
        self.m = m
        self.R = R
        self.Co = Co
        self.Cz = Cz
        self.TraceX = np.array([self.X])
        self.TraceY = np.array([self.Y])
    def DrawPlanet(self, ax):
        self.Dp = ax.plot(self.X0, self.Y0, 'o', color=self.Co, markersize=self.R)[0]
        self.Dt = ax.plot(self.TraceX, self.TraceY, ':', color=self.Co)[0]
    def RefreshPlanet(self, X, Y, Vx, Vy):
        self.X = X; self.Y = Y; self.Vx = Vx; self.Vy = Vy
        self.TraceX = np.append(self.TraceX, X)
        self.TraceY = np.append(self.TraceY, Y)
        self.Dp.set_data(X, Y)
        self.Dt.set_data(self.TraceX, self.TraceY)

class Satellite:
    def __init__(self, parent, angle, r_orbit, G):
        self.parent = parent
        self.R = parent.R / 1000.0
        self.X = parent.X + r_orbit * np.cos(angle)
        self.Y = parent.Y + r_orbit * np.sin(angle)
        v_orbit = 1.1 * sqrt(G * parent.m / r_orbit)
        self.Vx = parent.Vx - v_orbit * np.sin(angle)
        self.Vy = parent.Vy + v_orbit * np.cos(angle)
        self.TraceX = np.array([self.X])
        self.TraceY = np.array([self.Y])
        self.destroyed = False
    def DrawSatellite(self, ax):
        self.Dp = ax.plot(self.X, self.Y, 'o', color='yellow', markersize=3)[0]
        self.Dt = ax.plot(self.TraceX, self.TraceY, '-', color='yellow', linewidth=1)[0]
    def RefreshSatellite(self, X, Y, Vx, Vy):
        self.X = X; self.Y = Y; self.Vx = Vx; self.Vy = Vy
        self.TraceX = np.append(self.TraceX, X)
        self.TraceY = np.append(self.TraceY, Y)
        self.Dp.set_data(X, Y)
        self.Dt.set_data(self.TraceX, self.TraceY)

class SpaceShip:
    def __init__(self, X0, Y0, Vx0, Vy0, Fmax, R, Co, Cz):
        self.X0 = float(X0)
        self.Y0 = float(Y0)
        self.Vx0 = float(Vx0)
        self.Vy0 = float(Vy0)
        self.X = float(X0)
        self.Y = float(Y0)
        self.Vx = float(Vx0)
        self.Vy = float(Vy0)
        self.m = 1.0
        self.R = R
        self.Co = Co
        self.Cz = Cz
        self.Fmax = float(Fmax)
        self.TraceX = np.array([self.X])
        self.TraceY = np.array([self.Y])
        self.Ship_X = self.R * np.array([-0.4, 0, 0.5, 1, 0.5, 0,
                                         -0.4, -0.7, -0.15, 0,
                                         -0.4, -0.4, -0.7, -0.15, 0])
        self.Ship_Y = self.R * np.array([0.2, 0.3, 0.25, 0, -0.25, -0.3,
                                         -0.2, -0.5, -0.5, -0.3,
                                         -0.2, 0.2, 0.5, 0.5, 0.3])
        self.Flame_X = self.R * np.array([0, -0.4, -0.3, -0.8, -0.7, -1,
                                          -0.7, -0.8, -0.3, -0.4, 0])
        self.Flame_Y = self.R * np.array([0.2, 0.25, 0.2, 0.18, 0.1, 0,
                                          -0.1, -0.18, -0.2, -0.25, -0.2])
    def DrawShip(self, ax):
        self.Dfl = ax.plot(self.X0 + self.Ship_X[0]*self.R + self.Flame_X*self.R,
                           self.Y0 + self.Flame_Y*self.R, color=self.Cz)[0]
        self.Dsh = ax.plot(self.X0 + self.Ship_X*self.R, self.Y0 + self.Ship_Y*self.R, color=self.Co)[0]
        self.Dt = ax.plot(self.TraceX, self.TraceY, ':', color=self.Co)[0]
    def RefreshShip(self, X, Y, Vx, Vy, Phi, Fdv):
        self.X = X; self.Y = Y; self.Vx = Vx; self.Vy = Vy
        self.TraceX = np.append(self.TraceX, X)
        self.TraceY = np.append(self.TraceY, Y)
        RShip_X, RShip_Y = Rot2D(self.Ship_X, self.Ship_Y, np.pi/2 - Phi)
        RFlame_X, RFlame_Y = Rot2D(self.Flame_X*Fdv + self.Ship_X[0], self.Flame_Y, np.pi/2 - Phi)
        self.Dsh.set_data(self.X + RShip_X, self.Y + RShip_Y)
        self.Dfl.set_data(self.X + RFlame_X, self.Y + RFlame_Y)
        self.Dt.set_data(self.TraceX, self.TraceY)

class PlanetSystem:
    def __init__(self, Planets, Ship='Null', Xlim=[-12,12], Ylim=[-12,12], G=10):
        self.Planets = Planets
        self.N = len(self.Planets)
        self.X = np.array([p.X for p in self.Planets], 'float')
        self.Y = np.array([p.Y for p in self.Planets], 'float')
        self.VX = np.array([p.Vx for p in self.Planets], 'float')
        self.VY = np.array([p.Vy for p in self.Planets], 'float')
        self.Ship = Ship
        self.Satellites = []   # Поддержка спутников
        if self.Ship != 'Null':
            self.ShX = float(Ship.X)
            self.ShY = float(Ship.Y)
            self.ShVX = float(Ship.Vx)
            self.ShVY = float(Ship.Vy)
            self.ShPhi = 0.0
            self.ShFdv = 0.0
        self.Xlim = Xlim
        self.Ylim = Ylim
        self.G = G
        self.Camera = 0

    def SetGravity(self, G):
        self.G = G

    def DrawSystem(self, ax):
        self.ax = ax
        for planet in self.Planets:
            planet.DrawPlanet(self.ax)
        if self.Ship != 'Null':
            self.Ship.DrawShip(self.ax)
        for sat in self.Satellites:
            sat.DrawSatellite(self.ax)

    def RefreshSystem(self):
        if getattr(self, 'Camera', 0) == 1:
            self.ax.axis('scaled')
            self.ax.set(xlim=[self.Xlim[0] + self.ShX, self.Xlim[1] + self.ShX],
                        ylim=[self.Ylim[0] + self.ShY, self.Ylim[1] + self.ShY])
        for i in range(self.N):
            self.Planets[i].RefreshPlanet(self.X[i], self.Y[i], self.VX[i], self.VY[i])
        if self.Ship != 'Null':
            self.Ship.RefreshShip(self.ShX, self.ShY, self.ShVX, self.ShVY, self.ShPhi, self.ShFdv)
        for sat in self.Satellites:
            if not sat.destroyed:
                sat.RefreshSatellite(sat.X, sat.Y, sat.Vx, sat.Vy)

    def GetMoveEquations(self):
        X = sp.symbols('x:' + str(self.N))
        Y = sp.symbols('y:' + str(self.N))
        VX = sp.symbols('Vx:' + str(self.N))
        VY = sp.symbols('Vy:' + str(self.N))
        DX = VX
        DY = VY
        gravity_factor = 0.01
        DVX = [sum([self.Planets[j].m * (X[j] - X[i]) * self.G * gravity_factor /
                    ((X[j] - X[i])**2 + (Y[j] - Y[i])**2)**1.5
                    for j in range(self.N) if j != i])
               for i in range(self.N)]
        DVY = [sum([self.Planets[j].m * (Y[j] - Y[i]) * self.G * gravity_factor /
                    ((X[j] - X[i])**2 + (Y[j] - Y[i])**2)**1.5
                    for j in range(self.N) if j != i])
               for i in range(self.N)]
        self.MoveEquations = sp.lambdify([X, Y, VX, VY], [DX, DY, DVX, DVY], 'numpy')
        if self.Ship != 'Null':
            ShX = sp.symbols('Shx')
            ShY = sp.symbols('Shy')
            ShVX = sp.symbols('ShVx')
            ShVY = sp.symbols('ShVy')
            Phi = sp.symbols('phi')
            Fdv = sp.symbols('fdv')
            ShDX = ShVX
            ShDY = ShVY
            ShDVX = sum([self.Planets[j].m * (X[j] - ShX) * self.G /
                         ((X[j] - ShX)**2 + (Y[j] - ShY)**2)**1.5
                         for j in range(self.N)]) + self.Ship.Fmax * Fdv * sp.sin(Phi)
            ShDVY = sum([self.Planets[j].m * (Y[j] - ShY) * self.G /
                         ((X[j] - ShX)**2 + (Y[j] - ShY)**2)**1.5
                         for j in range(self.N)]) + self.Ship.Fmax * Fdv * sp.cos(Phi)
            self.ShipMoveEquations = sp.lambdify([X, Y, VX, VY, ShX, ShY, ShVX, ShVY, Phi, Fdv],
                                                 [ShDX, ShDY, ShDVX, ShDVY], 'numpy')

    def SatelliteAcceleration(self, x, y, parent):
        ax_val = 0.0
        ay_val = 0.0
        gravity_factor = 0.5  # для чужих планет
        epsilon = 0.1
        for planet in self.Planets:
            dx = planet.X - x
            dy = planet.Y - y
            r_sq = dx**2 + dy**2 + epsilon**2
            r3 = r_sq**1.5
            if planet == parent:
                ax_val += self.G * planet.m * dx / r3
                ay_val += self.G * planet.m * dy / r3
            else:
                ax_val += gravity_factor * self.G * planet.m * dx / r3
                ay_val += gravity_factor * self.G * planet.m * dy / r3
        ship_threshold = 2.0
        if self.Ship != 'Null':
            dx = self.ShX - x
            dy = self.ShY - y
            r = np.sqrt(dx**2 + dy**2)
            if r < ship_threshold:
                ax_val += self.G * self.Ship.m * dx / ((dx**2+dy**2+epsilon**2)**1.5)
                ay_val += self.G * self.Ship.m * dy / ((dx**2+dy**2+epsilon**2)**1.5)
        return ax_val, ay_val

# ------------------------
# Устанавливаем __module__ для pickle
# ------------------------
Planet.__module__ = "__main__"
Satellite.__module__ = "__main__"
SpaceShip.__module__ = "__main__"
PlanetSystem.__module__ = "__main__"

# ------------------------
# Приложение с UI
# ------------------------
class SpaceSimulationWidget(QMainWindow, SpaceWindow.Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.setWindowTitle("Просторы галактики")
        self.resize(1866, 945)
        self.StartButton.clicked.connect(self.HereAreWeGo)

    def HereAreWeGo(self):
        global Simulation, PS, multik
        from matplotlib.animation import FuncAnimation
        if not hasattr(self, 'Simulation'):
            self.Simulation = 0
        if self.Simulation == 0:
            Filename = self.File_Edit.text() + '.universe'
            dt = float(self.Step_Edit.text())
            Xmin = float(self.Xmin_Edit.text())
            Xmax = float(self.Xmax_Edit.text())
            Ymin = float(self.Ymin_Edit.text())
            Ymax = float(self.Ymax_Edit.text())
            self.spacewidget.canvas.axes.clear()
            self.spacewidget.canvas.axes.axis('scaled')
            self.spacewidget.canvas.axes.set(xlim=[Xmin, Xmax], ylim=[Ymin, Ymax])
            DrawTheSpace(self.spacewidget.canvas.axes)
            with open(Filename, 'rb') as file:
                data = pickle.load(file)
                PS = data['PS']
            if self.Camera_Ship.isChecked():
                PS.Camera = 1
            PS.DrawSystem(self.spacewidget.canvas.axes)
            self.spacewidget.canvas.show()
            PS.SetGravity(10)
            PS.GetMoveEquations()
            def kadr(i):
                DX, DY, DVX, DVY = PS.MoveEquations(PS.X, PS.Y, PS.VX, PS.VY)
                XK1, YK1, VXK1, VYK1 = np.array(DX), np.array(DY), np.array(DVX), np.array(DVY)
                DX, DY, DVX, DVY = PS.MoveEquations(PS.X + dt/2*XK1, PS.Y + dt/2*YK1,
                                                    PS.VX + dt/2*VXK1, PS.VY + dt/2*VYK1)
                XK2, YK2, VXK2, VYK2 = np.array(DX), np.array(DY), np.array(DVX), np.array(DVY)
                DX, DY, DVX, DVY = PS.MoveEquations(PS.X + dt/2*XK2, PS.Y + dt/2*YK2,
                                                    PS.VX + dt/2*VXK2, PS.VY + dt/2*VYK2)
                XK3, YK3, VXK3, VYK3 = np.array(DX), np.array(DY), np.array(DVX), np.array(DVY)
                DX, DY, DVX, DVY = PS.MoveEquations(PS.X + dt*XK3, PS.Y + dt*YK3,
                                                    PS.VX + dt*VXK3, PS.VY + dt*VYK3)
                XK4, YK4, VXK4, VYK4 = np.array(DX), np.array(DY), np.array(DVX), np.array(DVY)
                PS.X += (XK1 + 2*XK2 + 2*XK3 + XK4)*dt/6
                PS.Y += (YK1 + 2*YK2 + 2*YK3 + YK4)*dt/6
                PS.VX += (VXK1 + 2*VXK2 + 2*VXK3 + VXK4)*dt/6
                PS.VY += (VYK1 + 2*VYK2 + 2*VYK3 + VYK4)*dt/6
                if PS.Ship != 'Null':
                    PS.ShPhi = self.Course_dial.value()/100
                    PS.ShFdv = self.Force_Slider.value()/100
                    # RK4 для корабля
                    DShX, DShY, DShVX, DShVY = PS.ShipMoveEquations(PS.X, PS.Y, PS.VX, PS.VY,
                                                                       PS.ShX, PS.ShY, PS.ShVX, PS.ShVY,
                                                                       PS.ShPhi, PS.ShFdv)
                    ShXK1, ShYK1, ShVXK1, ShVYK1 = np.array(DShX), np.array(DShY), np.array(DShVX), np.array(DShVY)
                    DShX, DShY, DShVX, DShVY = PS.ShipMoveEquations(PS.X + dt/2*XK1, PS.Y + dt/2*YK1,
                                                                       PS.VX + dt/2*VXK1, PS.VY + dt/2*VYK1,
                                                                       PS.ShX + dt/2*ShXK1, PS.ShY + dt/2*ShYK1,
                                                                       PS.ShVX + dt/2*ShVXK1, PS.ShVY + dt/2*ShVYK1,
                                                                       PS.ShPhi, PS.ShFdv)
                    ShXK2, ShYK2, ShVXK2, ShVYK2 = np.array(DShX), np.array(DShY), np.array(DShVX), np.array(DShVY)
                    DShX, DShY, DShVX, DShVY = PS.ShipMoveEquations(PS.X + dt/2*XK2, PS.Y + dt/2*YK2,
                                                                       PS.VX + dt/2*VXK2, PS.VY + dt/2*VYK2,
                                                                       PS.ShX + dt/2*ShXK2, PS.ShY + dt/2*ShYK2,
                                                                       PS.ShVX + dt/2*ShVXK2, PS.ShVY + dt/2*ShVYK2,
                                                                       PS.ShPhi, PS.ShFdv)
                    ShXK3, ShYK3, ShVXK3, ShVYK3 = np.array(DShX), np.array(DShY), np.array(DShVX), np.array(DShVY)
                    DShX, DShY, DShVX, DShVY = PS.ShipMoveEquations(PS.X + dt*XK3, PS.Y + dt*YK3,
                                                                       PS.VX + dt*VXK3, PS.VY + dt*VYK3,
                                                                       PS.ShX + dt*ShXK3, PS.ShY + dt*ShYK3,
                                                                       PS.ShVX + dt*ShVXK3, PS.ShVY + dt*ShVYK3,
                                                                       PS.ShPhi, PS.ShFdv)
                    ShXK4, ShYK4, ShVXK4, ShVYK4 = np.array(DShX), np.array(DShY), np.array(DShVX), np.array(DShVY)
                    PS.ShX += (ShXK1 + 2*ShXK2 + 2*ShXK3 + ShXK4)*dt/6
                    PS.ShY += (ShYK1 + 2*ShYK2 + 2*ShYK3 + ShYK4)*dt/6
                    PS.ShVX += (ShVXK1 + 2*ShVXK2 + 2*ShVXK3 + ShVXK4)*dt/6
                    PS.ShVY += (ShVYK1 + 2*ShVYK2 + 2*ShVYK3 + ShVYK4)*dt/6

                for sat in PS.Satellites:
                    if sat.destroyed:
                        continue
                    def dSat(s):
                        x, y, vx, vy = s
                        ax_sat, ay_sat = PS.SatelliteAcceleration(x, y, sat.parent)
                        return np.array([vx, vy, ax_sat, ay_sat])
                    state0 = np.array([sat.X, sat.Y, sat.Vx, sat.Vy])
                    k1_sat = dSat(state0)
                    k2_sat = dSat(state0 + dt/2 * k1_sat)
                    k3_sat = dSat(state0 + dt/2 * k2_sat)
                    k4_sat = dSat(state0 + dt * k3_sat)
                    new_state = state0 + dt/6 * (k1_sat + 2*k2_sat + 2*k3_sat + k4_sat)
                    sat.X, sat.Y, sat.Vx, sat.Vy = new_state
                    # При сближении со своей планетой – уничтожаем спутник
                    if np.hypot(sat.X - sat.parent.X, sat.Y - sat.parent.Y) < 0.01:
                        PS.ax.plot(sat.X, sat.Y, '*', color='red', markersize=15)
                        sat.destroyed = True

                PS.RefreshSystem()
                objs = []
                for planet in PS.Planets:
                    objs.append(planet.Dp)
                    objs.append(planet.Dt)
                if PS.Ship != 'Null':
                    objs.append(PS.Ship.Dsh)
                    objs.append(PS.Ship.Dfl)
                    objs.append(PS.Ship.Dt)
                for sat in getattr(PS, 'Satellites', []):
                    if not sat.destroyed:
                        objs.append(sat.Dp)
                        objs.append(sat.Dt)
                return objs
            multik = FuncAnimation(self.spacewidget.canvas.figure, kadr, interval=dt*100, blit=True)
            self.spacewidget.canvas.draw()
            self.Simulation = 1
        else:
            self.Simulation = 0
            multik.repeat = False
            self.spacewidget.canvas.axes.clear()
            self.spacewidget.canvas.figure.clf()
            self.spacewidget.canvas.draw()

# ------------------------
# Запуск приложения
# ------------------------
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = SpaceSimulationWidget()
    window.show()
    sys.exit(app.exec_())
