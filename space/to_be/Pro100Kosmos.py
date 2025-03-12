import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pickle

# ------------------------
# Класс планеты (без изменений)
# ------------------------
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


# ------------------------
# Класс спутника (новый)
# ------------------------
class Satellite:
    def __init__(self, parent, angle, r_orbit, G):
        """
        parent : Planet (планета, вокруг которой вращается)
        angle : начальный угол орбиты (радианы)
        r_orbit : радиус орбиты
        G : гравитационная постоянная
        """
        self.parent = parent
        # Берём визуальный радиус спутника поменьше, чтобы отличался
        self.R = parent.R / 2.0
        # Начальные координаты: смещение от планеты
        self.X = parent.X + r_orbit*np.cos(angle)
        self.Y = parent.Y + r_orbit*np.sin(angle)
        # Устойчивое круговое движение: v = sqrt(G*M / r)
        v_orbit = np.sqrt(G * parent.m / r_orbit)
        # Начальная скорость = скорость планеты + орбитальная компонента
        self.Vx = parent.Vx - v_orbit*np.sin(angle)
        self.Vy = parent.Vy + v_orbit*np.cos(angle)
        self.TraceX = np.array([self.X])
        self.TraceY = np.array([self.Y])
        self.destroyed = False  # Флаг "уничтожен"

    def DrawSatellite(self, ax):
        # 's' - квадрат, markersize = 2*self.R
        self.Dp = ax.plot(self.X, self.Y, 's', color='yellow', markersize=self.R*2)[0]
        self.Dt = ax.plot(self.TraceX, self.TraceY, '--', color='yellow')[0]

    def RefreshSatellite(self, X, Y, Vx, Vy):
        self.X = X
        self.Y = Y
        self.Vx = Vx
        self.Vy = Vy
        self.TraceX = np.append(self.TraceX, X)
        self.TraceY = np.append(self.TraceY, Y)
        self.Dp.set_data(X, Y)
        self.Dt.set_data(self.TraceX, self.TraceY)


# ------------------------
# Класс системы планет
# ------------------------
class PlanetSystem:
    def __init__(self,Planets,Xlim,Ylim,G=1):
        self.Planets=Planets
        self.N = len(self.Planets)
        self.X = np.array([p.X for p in self.Planets],'float')
        self.Y = np.array([p.Y for p in self.Planets],'float')
        self.VX = np.array([p.Vx for p in self.Planets],'float')
        self.VY = np.array([p.Vy for p in self.Planets],'float')
        self.Xlim = Xlim
        self.Ylim = Ylim
        self.G=G
        self.Satellites = []  # Список спутников

    def DrawPlanets(self):
        # Создаём фигуру и ось
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(1, 1, 1)
        self.ax.axis('equal')  # Чтобы не искажались планеты
        self.ax.set(xlim=self.Xlim, ylim=self.Ylim)
        # Заливаем чёрный фон
        self.ax.set_facecolor('black')
        # Рисуем планеты
        for planet in self.Planets:
            planet.DrawPlanet(self.ax)
        # Рисуем спутники
        for sat in self.Satellites:
            sat.DrawSatellite(self.ax)

    def RefreshSystem(self):
        # Обновляем положение планет
        for i in range(self.N):
            self.Planets[i].RefreshPlanet(self.X[i], self.Y[i], self.VX[i], self.VY[i])
        # Обновляем положение спутников
        for sat in self.Satellites:
            if not sat.destroyed:
                sat.RefreshSatellite(sat.X, sat.Y, sat.Vx, sat.Vy)

    def GetMoveEquations(self):
        X = sp.symbols('x:'+str(self.N))
        Y = sp.symbols('y:'+str(self.N))
        VX = sp.symbols('Vx:'+str(self.N))
        VY = sp.symbols('Vy:'+str(self.N))

        DX = VX
        DY = VY
        # Сумма гравитационных сил для каждой планеты
        DVX = [sum([self.Planets[j].m * (X[j]-X[i]) * self.G /
                    ((X[j]-X[i])**2 + (Y[j]-Y[i])**2)**1.5
                    for j in range(self.N) if j!=i])
               for i in range(self.N)]
        DVY = [sum([self.Planets[j].m * (Y[j]-Y[i]) * self.G /
                    ((X[j]-X[i])**2 + (Y[j]-Y[i])**2)**1.5
                    for j in range(self.N) if j!=i])
               for i in range(self.N)]

        # Сохраняем функцию для вычисления dX, dY, dVx, dVy (для планет)
        self.MoveEquations = sp.lambdify([X,Y,VX,VY],[DX,DY,DVX,DVY],'numpy')

    def SatelliteAcceleration(self, x, y):
        # Спутники не влияют ни на кого, но все планеты (и их G, m) влияют на спутник
        ax_val = 0.0
        ay_val = 0.0
        for planet in self.Planets:
            dx = planet.X - x
            dy = planet.Y - y
            rr = (dx**2 + dy**2)**1.5
            if rr < 1e-9:
                continue
            ax_val += self.G * planet.m * dx / rr
            ay_val += self.G * planet.m * dy / rr
        return ax_val, ay_val


# ------------------------
# Параметр f
# ------------------------
f = 2  # 1 - создать систему, 2 - загрузить и запустить

if f == 1:
    # Создаём 3 планеты
    Planet1 = Planet(10, 0, 0, 1.5, 40, 5, [1, 0, 0], [0, 1, 0])
    Planet2 = Planet(-10, 0, 0, -2.1, 30, 5, [0, 1, 0], [0, 0, 1])
    Planet3 = Planet(0, 0, 0, -0.33, 30, 5, [0, 0, 1], [1, 0, 0])

    # Создаём систему с границами [-12,12], G=10
    PS = PlanetSystem([Planet1, Planet2, Planet3], [-12,12], [-12,12], 10)

    # Добавляем спутники (пример: по одному на каждую планету)
    # r_orbit = 4 (для наглядности)
    # Начальный угол берём 0, π/2,  π и т.д. по желанию
    from math import pi, sqrt, sin, cos
    sat1 = Satellite(Planet1, angle=0, r_orbit=4, G=PS.G)
    sat2 = Satellite(Planet2, angle=np.pi/2, r_orbit=4, G=PS.G)
    sat3 = Satellite(Planet3, angle=np.pi/4, r_orbit=4, G=PS.G)

    PS.Satellites.append(sat1)
    PS.Satellites.append(sat2)
    PS.Satellites.append(sat3)

    # Сохраняем в файл Universe1.universe
    fimeName = "Universe1.universe"
    with open(fimeName,'wb') as file:
        pickle.dump({'PS': PS}, file)

    print("Система создана и сохранена в", fimeName)

elif f == 2:
    # Загрузка системы из файла
    Filename = "Universe1.universe"
    with open(Filename, 'rb') as file:
        PS = pickle.load(file)['PS']

    # Рисуем систему
    PS.DrawPlanets()
    PS.GetMoveEquations()

    dt = 0.001

    def kadr(i):
        # 1) Считаем приращения для планет методом Рунге–Кутты 4
        DX,DY,DVX,DVY = PS.MoveEquations(PS.X,PS.Y,PS.VX,PS.VY)
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

        # Обновляем координаты планет
        PS.X += (XK1+2*XK2+2*XK3+XK4) * dt/6
        PS.Y += (YK1+2*YK2+2*YK3+YK4) * dt/6
        PS.VX += (VXK1+2*VXK2+2*VXK3+VXK4) * dt/6
        PS.VY += (VYK1+2*VYK2+2*VYK3+VYK4) * dt/6

        # 2) Обновляем спутники
        #    Т.к. спутники "не влияют" на планеты, мы просто считаем их ускорение
        #    от каждой планеты (методом Рунге–Кутты 4)
        to_remove = []
        for sat in PS.Satellites:
            if sat.destroyed:
                continue

            def dSatellite(s):
                x, y, vx, vy = s
                ax, ay = PS.SatelliteAcceleration(x,y)
                return np.array([vx, vy, ax, ay])

            state0 = np.array([sat.X, sat.Y, sat.Vx, sat.Vy])
            k1 = dSatellite(state0)
            k2 = dSatellite(state0 + dt/2*k1)
            k3 = dSatellite(state0 + dt/2*k2)
            k4 = dSatellite(state0 + dt*k3)
            state_new = state0 + dt/6*(k1 + 2*k2 + 2*k3 + k4)

            sat.X, sat.Y, sat.Vx, sat.Vy = state_new

            # Проверка сближения со "своей" планетой или любой другой
            for planet in PS.Planets:
                dist = np.hypot(sat.X - planet.X, sat.Y - planet.Y)
                if dist < planet.R * 1.2:
                    # Уничтожаем спутник, рисуем взрыв
                    PS.ax.plot(sat.X, sat.Y, '*', color='red', markersize=15)
                    sat.destroyed = True
                    to_remove.append(sat)
                    break

            # Проверка сближения с другими планетами (в том же цикле выше),
            # при желании — корабль и т.д.

        # Удаляем уничтоженные спутники из списка
        for sat in to_remove:
            if sat in PS.Satellites:
                sat.Dp.remove()
                sat.Dt.remove()
                PS.Satellites.remove(sat)

        # 3) Обновляем отрисовку
        PS.RefreshSystem()

        # 4) Возвращаем список графических объектов для blit
        #    (планеты + их траектории + спутники)
        objs = []
        for planet in PS.Planets:
            objs.append(planet.Dp)
            objs.append(planet.Dt)
        for sat in PS.Satellites:
            if not sat.destroyed:
                objs.append(sat.Dp)
                objs.append(sat.Dt)
        return objs

    # Запуск анимации
    from matplotlib.animation import FuncAnimation
    multik = FuncAnimation(PS.fig, kadr, interval=0.001*100, blit=True)
    plt.show()
