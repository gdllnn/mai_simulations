import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pickle
from math import pi, sqrt

# ------------------------
# Класс Planet
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

# ------------------------
# Класс Satellite
# ------------------------
class Satellite:
    def __init__(self, parent, angle, r_orbit, G):
        self.parent = parent
        # Спутники будут очень малы по сравнению с планетами:
        self.R = parent.R / 1000.0
        self.X = parent.X + r_orbit * np.cos(angle)
        self.Y = parent.Y + r_orbit * np.sin(angle)
        # Устойчивая орбитальная скорость с небольшим увеличением (1.1)
        v_orbit = 1.1 * sqrt(G * parent.m / r_orbit)
        self.Vx = parent.Vx - v_orbit * np.sin(angle)
        self.Vy = parent.Vy + v_orbit * np.cos(angle)
        self.TraceX = np.array([self.X])
        self.TraceY = np.array([self.Y])
        self.destroyed = False
    def DrawSatellite(self, ax):
        # Спутники отрисовываются жёлтым цветом
        self.Dp = ax.plot(self.X, self.Y, 'o', color='yellow', markersize=3)[0]
        self.Dt = ax.plot(self.TraceX, self.TraceY, '-', color='yellow', linewidth=1)[0]
    def RefreshSatellite(self, X, Y, Vx, Vy):
        self.X = X; self.Y = Y; self.Vx = Vx; self.Vy = Vy
        self.TraceX = np.append(self.TraceX, X)
        self.TraceY = np.append(self.TraceY, Y)
        self.Dp.set_data(X, Y)
        self.Dt.set_data(self.TraceX, self.TraceY)

# ------------------------
# Класс SpaceShip
# ------------------------
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
        self.m = 1.0  # Масса корабля
        self.R = R
        self.Co = Co
        self.Cz = Cz
        self.Fmax = float(Fmax)
        self.TraceX = np.array([self.X])
        self.TraceY = np.array([self.Y])
        # Геометрия корабля (можно расширить при необходимости)
        self.Ship_X = self.R * np.array([-0.4, 0, 0.5, 1, 0.5, 0, -0.4, -0.7, -0.15, 0, -0.4, -0.4, -0.7, -0.15, 0])
        self.Ship_Y = self.R * np.array([0.2, 0.3, 0.25, 0, -0.25, -0.3, -0.2, -0.5, -0.5, -0.3, -0.2, 0.2, 0.5, 0.5, 0.3])
        self.Flame_X = self.R * np.array([0, -0.4, -0.3, -0.8, -0.7, -1, -0.7, -0.8, -0.3, -0.4, 0])
        self.Flame_Y = self.R * np.array([0.2, 0.25, 0.2, 0.18, 0.1, 0, -0.1, -0.18, -0.2, -0.25, -0.2])
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

# ------------------------
# Функция Rot2D
# ------------------------
def Rot2D(X, Y, Alpha):
    RX = X * np.cos(Alpha) - Y * np.sin(Alpha)
    RY = X * np.sin(Alpha) + Y * np.cos(Alpha)
    return RX, RY

# ------------------------
# Класс PlanetSystem
# ------------------------
class PlanetSystem:
    def __init__(self, Planets, Ship='Null', Xlim=[-12,12], Ylim=[-12,12], G=10):
        self.Planets = Planets
        self.N = len(self.Planets)
        self.X = np.array([p.X for p in self.Planets], 'float')
        self.Y = np.array([p.Y for p in self.Planets], 'float')
        self.VX = np.array([p.Vx for p in self.Planets], 'float')
        self.VY = np.array([p.Vy for p in self.Planets], 'float')
        self.Ship = Ship
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
        self.Satellites = []
    def DrawPlanets(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(1,1,1)
        self.ax.axis('equal')
        self.ax.set(xlim=self.Xlim, ylim=self.Ylim)
        self.ax.set_facecolor('black')
        for planet in self.Planets:
            planet.DrawPlanet(self.ax)
        if self.Ship != 'Null':
            self.Ship.DrawShip(self.ax)
        for sat in self.Satellites:
            sat.DrawSatellite(self.ax)
    def RefreshSystem(self):
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
        gravity_factor = 0.01  # уменьшенный коэффициент гравитации
        DVX = [sum([self.Planets[j].m * (X[j] - X[i]) * self.G * gravity_factor /
                    ((X[j] - X[i])**2 + (Y[j] - Y[i])**2)**1.5
                    for j in range(self.N) if j != i])
               for i in range(self.N)]
        DVY = [sum([self.Planets[j].m * (Y[j] - Y[i]) * self.G * gravity_factor /
                    ((X[j] - X[i])**2 + (Y[j] - Y[i])**2)**1.5
                    for j in range(self.N) if j != i])
               for i in range(self.N)]
        self.MoveEquations = sp.lambdify([X,Y,VX,VY], [DX, DY, DVX, DVY], 'numpy')
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
# Устанавливаем __module__ для корректной сериализации
Planet.__module__ = "__main__"
Satellite.__module__ = "__main__"
SpaceShip.__module__ = "__main__"
PlanetSystem.__module__ = "__main__"

# ------------------------
# Режим работы: Создание и сохранение системы
# ------------------------
f = 2  # Запустите с f = 1 для создания и сохранения

if f == 1:
    # Инициализируем планеты, корабль и спутники
    Planet1 = Planet(10, 0, 0, 0, 30, 5, [0,0,1], [1,0,0])
    Planet2 = Planet(-10, 0, 0, 0, 30, 5, [0,0,1], [1,0,0])
    Planet3 = Planet(0, 0, 0, 0, 30, 5, [0,0,1], [1,0,0])
    Ship = SpaceShip(5, 0, 0, 1, 100, 0.3, [1,1,1], [1,0,0])
    PS = PlanetSystem([Planet1, Planet2, Planet3], Ship)
    # Добавляем спутники с орбитальным радиусом = 2
    sat1 = Satellite(Planet1, angle=0, r_orbit=2, G=PS.G)
    sat2 = Satellite(Planet2, angle=pi/2, r_orbit=2, G=PS.G)
    sat3 = Satellite(Planet3, angle=pi/4, r_orbit=2, G=PS.G)
    PS.Satellites.append(sat1)
    PS.Satellites.append(sat2)
    PS.Satellites.append(sat3)
    fimeName = "Universe1.universe"
    with open(fimeName, 'wb') as file:
        pickle.dump({'PS': PS}, file, protocol=pickle.HIGHEST_PROTOCOL)
    print("Сделано")
    
elif f == 2:
    Filename = "Universe1.universe"
    with open(Filename, 'rb') as file:
        PS = pickle.load(file)['PS']
    PS.DrawPlanets()
    PS.GetMoveEquations()
    dt = 0.001
    def kadr(i):
        DX, DY, DVX, DVY = PS.MoveEquations(PS.X, PS.Y, PS.VX, PS.VY)
        XK1, YK1, VXK1, VYK1 = np.array(DX), np.array(DY), np.array(DVX), np.array(DVY)
        DX, DY, DVX, DVY = PS.MoveEquations(PS.X+dt/2*XK1, PS.Y+dt/2*YK1,
                                            PS.VX+dt/2*VXK1, PS.VY+dt/2*VYK1)
        XK2, YK2, VXK2, VYK2 = np.array(DX), np.array(DY), np.array(DVX), np.array(DVY)
        DX, DY, DVX, DVY = PS.MoveEquations(PS.X+dt/2*XK2, PS.Y+dt/2*YK2,
                                            PS.VX+dt/2*VXK2, PS.VY+dt/2*VYK2)
        XK3, YK3, VXK3, VYK3 = np.array(DX), np.array(DY), np.array(DVX), np.array(DVY)
        DX, DY, DVX, DVY = PS.MoveEquations(PS.X+dt*XK3, PS.Y+dt*YK3,
                                            PS.VX+dt*VXK3, PS.VY+dt*VYK3)
        XK4, YK4, VXK4, VYK4 = np.array(DX), np.array(DY), np.array(DVX), np.array(DVY)
        PS.X += (XK1+2*XK2+2*XK3+XK4)*dt/6
        PS.Y += (YK1+2*YK2+2*YK3+YK4)*dt/6
        PS.VX += (VXK1+2*VXK2+2*VXK3+VXK4)*dt/6
        PS.VY += (VYK1+2*VYK2+2*VYK3+VYK4)*dt/6

        # Обновляем спутники
        remove_list = []
        for sat in PS.Satellites:
            if sat.destroyed:
                continue
            def dSat(s, parent):
                x, y, vx, vy = s
                ax_sat, ay_sat = PS.SatelliteAcceleration(x, y, parent)
                return np.array([vx, vy, ax_sat, ay_sat])
            state0 = np.array([sat.X, sat.Y, sat.Vx, sat.Vy])
            k1_sat = dSat(state0, sat.parent)
            k2_sat = dSat(state0 + dt/2*k1_sat, sat.parent)
            k3_sat = dSat(state0 + dt/2*k2_sat, sat.parent)
            k4_sat = dSat(state0 + dt*k3_sat, sat.parent)
            new_state = state0 + dt/6*(k1_sat + 2*k2_sat + 2*k3_sat + k4_sat)
            sat.X, sat.Y, sat.Vx, sat.Vy = new_state
            for planet in PS.Planets:
                threshold = 0.1
                if np.hypot(sat.X-planet.X, sat.Y-planet.Y) < threshold:
                    PS.ax.plot(sat.X, sat.Y, '*', color='red', markersize=15)
                    sat.destroyed = True
                    remove_list.append(sat)
                    break
        for sat in remove_list:
            if sat in PS.Satellites:
                sat.Dp.remove()
                sat.Dt.remove()
                PS.Satellites.remove(sat)
        PS.RefreshSystem()
        objs = []
        for planet in PS.Planets:
            objs.append(planet.Dp)
            objs.append(planet.Dt)
        if PS.Ship != 'Null':
            objs.append(PS.Ship.Dsh)
            objs.append(PS.Ship.Dfl)
            objs.append(PS.Ship.Dt)
        for sat in PS.Satellites:
            if not sat.destroyed:
                objs.append(sat.Dp)
                objs.append(sat.Dt)
        return objs

        from matplotlib.animation import FuncAnimation
    from matplotlib.animation import FuncAnimation
    multik = FuncAnimation(PS.fig, kadr, interval=dt*100, blit=True)
    plt.show()