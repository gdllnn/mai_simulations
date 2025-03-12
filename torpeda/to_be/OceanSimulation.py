from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
import OceanWindow
import ocean
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

# Функция для расчёта динамики лодки
def EquationsOfLodka(X, Y, Phi, Vx, Vy, Omega, F_dv, alpha):
    global m, J, a, Cl, Cb, Cvr, Sl, Sb, Rho
    Vl = Vx * np.cos(Phi) + Vy * np.sin(Phi)
    Vb = -Vx * np.sin(Phi) + Vy * np.cos(Phi)
    Fl = Rho * Cl * Sl * Vl * abs(Vl) / 2
    Fb = Rho * Cb * Sb * Vb * abs(Vb) / 2
    Ms = Rho * Cvr * Sb * Omega * abs(Omega) / 2
    Wx = (-Fl * np.cos(Phi) + Fb * np.sin(Phi) + F_dv * np.cos(Phi - alpha)) / m
    Wy = (-Fl * np.sin(Phi) - Fb * np.cos(Phi) + F_dv * np.sin(Phi - alpha)) / m
    Epsilon = (-Ms + a * F_dv * np.sin(alpha)) / J
    return [Vx, Vy, Omega, Wx, Wy, Epsilon]

# Функция для расчёта динамики торпеды
def EquationsOfTorpeda(X, Y, Phi, Vx, Vy, Omega, alpha):
    global mT, JT, aT, ClT, CbT, CvrT, SlT, SbT, Rho, VoT, R
    Vl = Vx * np.cos(Phi) + Vy * np.sin(Phi)
    Vb = -Vx * np.sin(Phi) + Vy * np.cos(Phi)
    Fl = Rho * ClT * SlT * Vl * abs(Vl) / 2
    Fb = Rho * CbT * SbT * Vb * abs(Vb) / 2
    Ms = Rho * CvrT * SbT * Omega * abs(Omega) / 2
    F_v = (VoT - Vl) * 3.14 * R**2 * (Vl + VoT) * Rho
    F_r = (Vl + VoT) * np.sin(alpha)**3 * 3.14 * R**2 * (Vl + VoT) * Rho / 2
    Wx = (-Fl * np.cos(Phi) + Fb * np.sin(Phi) + F_v * np.cos(Phi) + F_r * np.sin(Phi - alpha)) / mT
    Wy = (-Fl * np.sin(Phi) - Fb * np.cos(Phi) + F_v * np.sin(Phi) - F_r * np.cos(Phi - alpha)) / mT
    Epsilon = (-Ms + aT * F_r * np.cos(alpha)) / JT
    return [Vx, Vy, Omega, Wx, Wy, Epsilon]

# Функция для вычисления управляющего сигнала по курсу (резервный вариант)
def GetCourse(X, Y, XT, YT, PhiT):
    Beta = np.arctan2(Y - YT, X - XT) - PhiT
    while abs(Beta) > np.pi:
        Beta -= 2 * np.pi * np.sign(Beta)
    if abs(Beta) < 0.7 and (Y - YT)**2 + (X - XT)**2 < 30**2:
        V = 1
        alpha = Beta
        if alpha > 0.6:
            alpha = 0.6
        elif alpha < -0.6:
            alpha = -0.6
    else:
        V = 0
        alpha = 0
    return alpha, V

# Функция поворота координат на заданный угол
def Rot2D(X, Y, Alpha):
    RX = X * np.cos(Alpha) - Y * np.sin(Alpha)
    RY = X * np.sin(Alpha) + Y * np.cos(Alpha)
    return RX, RY

# Основной класс приложения
class OceanWidget(QMainWindow, OceanWindow.Ui_MainWindow):

    def __init__(self):
        QMainWindow.__init__(self)
        self.setupUi(self)
        self.setWindowTitle("Водная гладь")
        self.Start_Button.clicked.connect(self.HereAreWeGo)

    def HereAreWeGo(self):
        global m, J, a, Cl, Cb, Cvr, Sl, Sb, Rho, \
               X, Y, Phi, Vx, Vy, Omega, F_dv, alpha, t, F_dv_max, \
               mT, JT, aT, ClT, CbT, CvrT, SlT, SbT, VoT, R, \
               XT, YT, PhiT, VxT, VyT, OmegaT, Explode, r_E

        # Считываем параметры лодки из интерфейса
        m = float(self.m_Edit.text())
        a = float(self.a_Edit.text())
        J = float(self.J_Edit.text())
        h = float(self.h_Edit.text())
        Cl = float(self.Cl_Edit.text())
        Cb = float(self.Cb_Edit.text())
        Cvr = float(self.Cvr_Edit.text())
        Sl = float(self.Sl_Edit.text())
        Sb = float(self.Sb_Edit.text())
        Rho = float(self.Rho_Edit.text())

        # Считываем параметры торпеды
        mT = float(self.mT_Edit.text())
        aT = float(self.aT_Edit.text())
        JT = float(self.JT_Edit.text())
        R = float(self.RT_Edit.text())
        ClT = float(self.ClT_Edit.text())
        CbT = float(self.CbT_Edit.text())
        CvrT = float(self.CvrT_Edit.text())
        SlT = float(self.SlT_Edit.text())
        SbT = float(self.SbT_Edit.text())
        VoT = float(self.VoT_Edit.text())

        # Считываем начальные данные лодки
        X = float(self.X_Edit.text())
        Y = float(self.Y_Edit.text())
        Phi = float(self.Phi_Edit.text())
        Vx = float(self.Vx_Edit.text())
        Vy = float(self.Vy_Edit.text())
        Omega = float(self.Omega_Edit.text())

        # Считываем начальные данные торпеды
        XT = float(self.XT_Edit.text())
        YT = float(self.YT_Edit.text())
        PhiT = float(self.PhiT_Edit.text())
        VxT = float(self.VxT_Edit.text())
        VyT = float(self.VyT_Edit.text())
        OmegaT = float(self.OmegaT_Edit.text())

        # Настройка графики
        self.morewidget.canvas.axes.axis('equal')
        self.morewidget.canvas.axes.set(xlim=[-40 * a, 40 * a], ylim=[-40 * a, 40 * a])

        # Вычисляем контуры лодки и торпеды для отрисовки
        L_X = np.array([-a, -a/3, a/3, a, 5/3*a, a, a/3, -a/3, -a, -a])
        L_Y = np.array([0.4*a, a/2, a/2, 0.4*a, 0, -0.4*a, -a/2, -a/2, -0.4*a, 0.4*a])
        T_X = aT * np.array([-1.1, -0.9, -0.8, 0.75, 0.8, 0.75, -0.8, -0.9, -1.1, -1.1])
        T_Y = aT * np.array([0.1, 0.1, -0.1, -0.1, 0, 0.1, 0.1, -0.1, -0.1, 0.1])
        LodkaX, LodkaY = Rot2D(L_X, L_Y, Phi)
        TorpedaX, TorpedaY = Rot2D(T_X, T_Y, PhiT)
        DrawedLodka = self.morewidget.canvas.axes.plot(X + LodkaX, Y + LodkaY, color=[0, 0, 1])[0]
        DrawedTorpeda = self.morewidget.canvas.axes.plot(XT + TorpedaX, YT + TorpedaY, color=[1, 0, 0])[0]

        # Инициализация переменных симуляции
        F_dv_max = 200
        alpha = 0.1
        dt = 0.1
        t = 0
        Explode = 0
        r_E = 0.5

        # Рассчитываем эффективную скорость торпеды (VT) по её характеристикам
        VT = np.sqrt(np.pi * 2) * VoT * R / np.sqrt(ClT * SlT + 2 * np.pi * R**2)
        print("VT:", VT)
        # Начальное вычисление интерцепционного времени и корректировка PhiT (однократно)
        T0 = (np.sqrt( -(Y - YT)**2 * Vx**2 + 2*Vy*(Y - YT)*(X - XT)*Vx - (X - XT)**2 * Vy**2 + VT**2*(X**2 - 2*X*XT + XT**2 + (Y - YT)**2)) + (X - XT)*Vx + (Y - YT)*Vy) / (-Vy**2 - Vx**2 + VT**2)
        print("T0:", T0)
        PhiT = np.arctan2((Y + Vy * T0 - YT) / (VT * T0), (X + Vx * T0 - XT) / (VT * T0))
        print("Initial PhiT:", PhiT)

        # Основная функция, вызываемая на каждом шаге анимации
        def kadr(j):
            global X, Y, Phi, Vx, Vy, Omega, F_dv, alpha, t, F_dv_max, \
                   XT, YT, PhiT, VxT, VyT, OmegaT, Explode, r_E

            # Проверка на столкновение: если объекты слишком близко – запускаем эффект взрыва
            if ((X - XT)**2 + (Y - YT)**2 < (a + aT)**2):
                Explode = 1

            if Explode == 0:
                # Получаем управляющие воздействия для лодки
                alpha = self.alpha_dial.value() / 100
                F_dv = self.ForceSlider.value() * F_dv_max / 100

                # Вычисляем текущее расстояние между лодкой и торпедой
                current_distance = np.sqrt((X - XT)**2 + (Y - YT)**2)

                # *** Интерцепционное наведение ***
                # Решаем квадратное уравнение для времени перехвата:
                # (dx + Vx*T)^2 + (dy + Vy*T)^2 = (VT*T)^2,
                # где dx = X - XT, dy = Y - YT.
                dx = X - XT
                dy = Y - YT
                A = Vx**2 + Vy**2 - VT**2
                B = 2 * (dx * Vx + dy * Vy)
                C = dx**2 + dy**2
                discriminant = B**2 - 4 * A * C

                if A != 0 and discriminant >= 0:
                    T_candidate = (-B - np.sqrt(discriminant)) / (2 * A)
                    if T_candidate <= 0:
                        T_candidate = (-B + np.sqrt(discriminant)) / (2 * A)
                    if T_candidate > 0:
                        # Вычисляем точку перехвата для корабля
                        X_int = X + Vx * T_candidate
                        Y_int = Y + Vy * T_candidate
                        # Определяем желаемое направление от торпеды к точке перехвата
                        desired_heading = np.arctan2(Y_int - YT, X_int - XT)
                        alphaT = desired_heading - PhiT
                        # Ограничиваем управляющее воздействие
                        max_alpha = 0.6
                        if alphaT > max_alpha:
                            alphaT = max_alpha
                        elif alphaT < -max_alpha:
                            alphaT = -max_alpha
                        V = 1
                    else:
                        # Если положительного решения нет, используем резервное управление
                        alphaT, V = GetCourse(X, Y, XT, YT, PhiT)
                else:
                    # Если нет решения (например, отрицательный дискриминант), используем резервное управление
                    alphaT, V = GetCourse(X, Y, XT, YT, PhiT)

                # Расчет динамики для лодки
                Vx, Vy, Omega, Wx, Wy, Epsilon = EquationsOfLodka(X, Y, Phi, Vx, Vy, Omega, F_dv, alpha)
                # Расчет динамики для торпеды с использованием интерцепционного управления
                VxT, VyT, OmegaT, WxT, WyT, EpsilonT = EquationsOfTorpeda(XT, YT, PhiT, VxT, VyT, OmegaT, alphaT)

                # Обновляем состояния объектов методом Эйлера
                X += Vx * dt
                Y += Vy * dt
                Phi += Omega * dt
                Vx += Wx * dt
                Vy += Wy * dt
                Omega += Epsilon * dt

                XT += VxT * dt
                YT += VyT * dt
                PhiT += OmegaT * dt
                VxT += WxT * dt
                VyT += WyT * dt
                OmegaT += EpsilonT * dt

                t += dt

                # Обновляем отрисовку: поворачиваем контуры объектов согласно их текущим углам
                LodkaX, LodkaY = Rot2D(L_X, L_Y, Phi)
                DrawedLodka.set_data(X + LodkaX, Y + LodkaY)
                TorpedaX, TorpedaY = Rot2D(T_X, T_Y, PhiT)
                DrawedTorpeda.set_data(XT + TorpedaX, YT + TorpedaY)

                print("PhiT:", PhiT)
                if V == 0:
                    DrawedTorpeda.set_color([0, 0.7, 0])
                else:
                    DrawedTorpeda.set_color([1, 0, 0])
            else:
                # Эффект взрыва: визуальное отображение
                Betas = np.linspace(0, 2, 41) * np.pi
                r_Explosion = (1 + np.sin(5 * Betas)) * 5 * a * (1 + 1 / (-1 - r_E))
                X_Explosion = np.concatenate([r_Explosion * np.cos(Betas), r_Explosion * np.sin(Betas)])
                Y_Explosion = np.concatenate([r_Explosion * np.sin(Betas), r_Explosion * np.cos(Betas)])
                DrawedTorpeda.set_data(XT + X_Explosion, YT + Y_Explosion)
                DrawedTorpeda.set_marker('o')
                r_E += 0.01

            return [DrawedLodka, DrawedTorpeda]

        # Запускаем анимацию с заданным интервалом
        multik = FuncAnimation(self.morewidget.canvas.figure, kadr, interval=dt * 200, blit=True)
        self.morewidget.canvas.draw()
        return

# Запуск приложения
app = QApplication([])
window = OceanWidget()
window.show()
app.exec_()
