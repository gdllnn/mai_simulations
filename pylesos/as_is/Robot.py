import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from PIL import Image
import pickle

# ------------------------------------------------
# 1) ПРОВЕРКА ЦВЕТА
# ------------------------------------------------
def IsColor(C1, C2):
    """
    Проверяет, похож ли цвет C1 на C2 (евклидово расстояние < 0.5 в пространстве RGB/255).
    """
    R = (C1 - C2)/255
    return (R[0]**2 + R[1]**2 + R[2]**2) < 0.5

# ------------------------------------------------
# 2) ЗАГРУЗКА КАРТЫ
# ------------------------------------------------
jpg = Image.open("Map1.jpg")
fig = plt.figure(figsize=[17,9])
ax = fig.add_subplot(1,1,1)
ax.imshow(jpg)

im2arr = np.array(jpg)
print(im2arr.shape)
print(im2arr[4])

Razmer = 63
N_x = im2arr.shape[1] // Razmer
N_y = im2arr.shape[0] // Razmer

Map = np.zeros((N_x, N_y))

ax.plot([N_x*Razmer, N_x*Razmer],[0, N_y*Razmer],color=[0,1,0])
ax.plot([0, N_x*Razmer],[N_y*Razmer, N_y*Razmer],color=[0,1,0])

f=2  # 1->создать карту, 2->читать

if f==1:
    # Создаём карту и сохраняем
    for i in range(N_x):
        ax.plot([i*Razmer,i*Razmer],[0,N_y*Razmer], color=[0,1,0])
        print(i,'/',N_x)
        for j in range(N_y):
            if i==0:
                ax.plot([0,N_x*Razmer],[j*Razmer,j*Razmer], color=[0,1,0])
            N_p=0
            N_s=0
            N_c=0
            # Перебор пикселей в ячейке
            for iz in range(Razmer):
                for jz in range(Razmer):
                    pix = im2arr[j*Razmer+jz][i*Razmer+iz]
                    if IsColor(pix,[0,0,0]):
                        N_p+=1
                    elif IsColor(pix,[0,0,255]):
                        N_s+=1
                    elif IsColor(pix,[255,0,0]):
                        N_c+=1
            if N_s>=10:
                Map[i,j] = -2  # Старт
            if N_c>=10:
                Map[i,j] = -3  # Цель
            if N_p>=10:
                Map[i,j] = -1  # Препятствие

    with open("Map1.map",'wb') as file:
        pickle.dump({"Map":Map}, file)
    print("Карта создана и сохранена")

elif f==2:
    for i in range(N_x):
        ax.plot([i*Razmer,i*Razmer],[0,N_y*Razmer],color=[0,1,0])
        print(i,'/',N_x)
        for j in range(N_y):
            if i==0:
                ax.plot([0,N_x*Razmer],[j*Razmer,j*Razmer],color=[0,1,0])
    with open("Map1.map",'rb') as file:
        Map = pickle.load(file)['Map']

# ------------------------------------------------
# 3) ОПРЕДЕЛЯЕМ СТАРТ/ЦЕЛЬ
# ------------------------------------------------
I_s=0
J_s=0
I_c=0
J_c=0
N_s=0
N_c=0
for i in range(N_x):
    for j in range(N_y):
        if Map[i,j]==-2:
            I_s+=i
            J_s+=j
            N_s+=1
        elif Map[i,j]==-3:
            I_c+=i
            J_c+=j
            N_c+=1

if N_s>0:
    I_s = int(round(I_s/N_s))
    J_s = int(round(J_s/N_s))
if N_c>0:
    I_c = int(round(I_c/N_c))
    J_c = int(round(J_c/N_c))

for i in range(N_x):
    for j in range(N_y):
        # Оставляем только одну стартовую и одну целевую
        if Map[i,j]==-2 and not (i==I_s and j==J_s):
            Map[i,j]=0
        if Map[i,j]==-3 and not (i==I_c and j==J_c):
            Map[i,j]=0

        # Рисуем
        if Map[i,j]==0:
            ax.plot((i+0.5)*Razmer,(j+0.5)*Razmer,'o',color=[1,1,1])
        elif Map[i,j]==-1:
            ax.plot((i+0.5)*Razmer,(j+0.5)*Razmer,'o',color=[0,0,0])
        elif Map[i,j]==-2:
            ax.plot((i+0.5)*Razmer,(j+0.5)*Razmer,'o',color=[0,0,1])
        elif Map[i,j]==-3:
            ax.plot((i+0.5)*Razmer,(j+0.5)*Razmer,'o',color=[1,0,0])

# ------------------------------------------------
# 4) ВОЛНОВОЙ АЛГОРИТМ
# ------------------------------------------------
def wave_algorithm(start, Map):
    """
    Ищет путь к ближайшей -3, объезжая -1.
    Возвращает список (i,j) или [] если не найдено.
    """
    steps=[[ (start[0], start[1]) ]]
    InverseRoute=[[start[0],start[1]]]
    i=1
    fin=0

    while len(steps[-1]) and fin==0:
        wave=[]
        for (px,py) in steps[-1]:
            # вправо
            if px+1<Map.shape[0]:
                if Map[px+1,py]==0:
                    wave.append((px+1,py))
                    Map[px+1,py] = i
                elif Map[px+1,py]==-3:
                    wave.append((px+1,py))
                    Map[px+1,py]= i
                    fin=1
                    InverseRoute=[[px+1,py]]
                    break
            # влево
            if px-1>=0 and fin==0:
                if Map[px-1,py]==0:
                    wave.append((px-1,py))
                    Map[px-1,py]= i
                elif Map[px-1,py]==-3:
                    wave.append((px-1,py))
                    Map[px-1,py]= i
                    fin=1
                    InverseRoute=[[px-1,py]]
                    break
            # вверх
            if py+1<Map.shape[1] and fin==0:
                if Map[px,py+1]==0:
                    wave.append((px,py+1))
                    Map[px,py+1]=i
                elif Map[px,py+1]==-3:
                    wave.append((px,py+1))
                    Map[px,py+1]= i
                    fin=1
                    InverseRoute=[[px,py+1]]
                    break
            # вниз
            if py-1>=0 and fin==0:
                if Map[px,py-1]==0:
                    wave.append((px,py-1))
                    Map[px,py-1]= i
                elif Map[px,py-1]==-3:
                    wave.append((px,py-1))
                    Map[px,py-1]= i
                    fin=1
                    InverseRoute=[[px,py-1]]
                    break
        steps.append(wave)
        i+=1

    if fin==0:
        print("Цель не найдена!")
        return []

    # Восстанавливаем
    N_steps = len(steps)
    for _ in range(N_steps):
        pLast = InverseRoute[-1]
        val = Map[pLast[0], pLast[1]]
        added=False
        for dx,dy in [(1,0),(-1,0),(0,1),(0,-1)]:
            nx,ny = pLast[0]+dx, pLast[1]+dy
            if 0<=nx<Map.shape[0] and 0<=ny<Map.shape[1]:
                if Map[nx,ny]==val-1:
                    InverseRoute.append([nx,ny])
                    added=True
                    break
                elif Map[nx,ny]==-2:  # старт
                    InverseRoute.append([nx,ny])
                    added=True
                    break
        if not added:
            break

    return InverseRoute[::-1]


# ------------------------------------------------
# 5) ПРОВЕРКА line-of-sight
# ------------------------------------------------
def line_of_sight(Map, x1,y1, x2,y2):
    """
    Проверяем, нет ли -1 между (x1,y1) и (x2,y2).
    Дробим отрезок на мелкие шаги ~0.25 клетки.
    """
    steps = int(max(abs(x2-x1),abs(y2-y1))*4 +1)
    for s in range(steps+1):
        t= s/steps
        ix = int(round(x1 + (x2-x1)*t))
        iy = int(round(y1 + (y2-y1)*t))
        if ix<0 or ix>=Map.shape[0] or iy<0 or iy>=Map.shape[1]:
            return False
        if Map[ix,iy] == -1:
            return False
    return True

# ------------------------------------------------
# 6) СГЛАЖИВАНИЕ (LOOKAHEAD=3)
# ------------------------------------------------
def simplify_path_los(route, Map, LOOKAHEAD=3):
    """
    route - список [(i0,j0), (i1,j1), ...]
    Пытаемся пропускать до 3 шагов вперёд (или меньше, если не хватает).
    """
    if len(route)<2:
        return route
    
    new_route=[route[0]]
    k=0
    while k<len(route)-1:
        best = k+1
        # идём от (k+1) до (k+LOOKAHEAD+1), ищем самую дальнюю точку,
        # до которой есть видимость
        for test_i in range(k+1, min(len(route), k+1+LOOKAHEAD+1)):
            x0,y0 = route[k]
            xT,yT = route[test_i]
            if line_of_sight(Map, x0,y0, xT,yT):
                best = test_i
        new_route.append(route[best])
        k=best
    return new_route

# ------------------------------------------------
# 7) ИСПОЛЬЗУЕМ ВОЛНОВОЙ АЛГОРИТМ, СГЛАЖИВАЕМ
# ------------------------------------------------
RouteOriginal = wave_algorithm((I_s,J_s), Map)

if not RouteOriginal:  
    print("Путь не найден, выходим.")
    plt.show()
    exit()

# Рисуем исходный (зелёные точки)
for pt in RouteOriginal:
    if pt != [I_s,J_s]:
        ax.plot((pt[0]+0.5)*Razmer, (pt[1]+0.5)*Razmer,
                'o', color=[0,0.5,0], markersize=4)

# Сглаживаем (LOOKAHEAD=3)
RouteSmoothed = simplify_path_los(RouteOriginal, Map, LOOKAHEAD=3)

# Рисуем сглаженный (красные точки)
for pt in RouteSmoothed:
    ax.plot((pt[0]+0.5)*Razmer, (pt[1]+0.5)*Razmer,
            'o', color=[1,0,0], markersize=4)

# Переводим в формат (x+0.5, y+0.5), добавляем [1,1] в начало
RouteSmoothed.insert(0, [1,1])
Route = [[p[0]+0.5, p[1]+0.5] for p in RouteSmoothed]

near_index=0

# ------------------------------------------------
# 8) ПАРАМЕТРЫ РОБОТА И ФИЗИЧЕСКАЯ МОДЕЛЬ
# ------------------------------------------------
m=1
R=10
J= m*(R**2)/2
mu=0.1
g=9.81

X=(I_s+0.5)
Y=(J_s+0.5)
Phi= -np.arctan2(Route[1][1]-Route[0][1],
                 Route[1][0]-Route[0][0])
VX=0
VY=0
Omega=0

dt=0.015
F0 = 0.7 * mu*m*g

# Окружность + луч (последняя точка= (0,0))
alpha_R = np.linspace(0,2*np.pi,51)
X_R = np.append(Razmer/2*np.cos(alpha_R), 0)
Y_R = np.append(Razmer/2*np.sin(alpha_R), 0)
Robo = ax.plot(X_R, Y_R, color=[1,0,1])[0]

beta=0
beta_prev=0
betas=0
rho=0
rho_prev=0
rhos=0
K1,K2,K3,K4,K5= 2,7,6,10,0.5

t=[0]
N_sred=50
y1,y2,y3,y4 = (np.zeros(N_sred) for _ in range(4))
I_sred=0

stop_distance=0.3  # порог остановки возле последней точки

def Rot2D(Xv, Yv, Alpha):
    RX = Xv*np.cos(Alpha) - Yv*np.sin(Alpha)
    RY = Xv*np.sin(Alpha) + Yv*np.cos(Alpha)
    return RX, RY

def get_nearest_points(index):
    global near_index
    if index>=len(Route)-1:
        return Route[-1],Route[-1]

    dist_curr = (Route[index][0]-X)**2 + (Route[index][1]-Y)**2
    dist_next = (Route[index+1][0]-X)**2 + (Route[index+1][1]-Y)**2
    if dist_curr>dist_next:
        index+=1
        near_index+=1
    if index>=len(Route)-1:
        return Route[-1],Route[-1]
    return Route[index], Route[index+1]

def get_distance(nA, nB, X, Y):
    A = nA[1]-nB[1]
    B = nB[0]-nA[0]
    C = -(nB[0]*nA[1]) + nB[0]*nB[1] \
        + nA[0]*nB[1] - nB[0]*nB[1]
    denom = np.sqrt(A**2 + B**2)
    if denom<1e-12:
        return 0.
    return (A*X + B*Y + C)/denom

def RobotSystemOfEquations(X, Y, Phi, Vx, Vy, Omega, F1, F2):
    V = Vx*np.cos(Phi) - Vy*np.sin(Phi)
    Ftr = 2/np.pi*np.arctan(10000*V)*mu*m*g
    Vst = (F1+F2 - Ftr)/m
    Phitt= (F1-F2)*R/J
    Xtt = Vst*np.cos(Phi) - V*Omega*np.sin(Phi)
    Ytt = -Vst*np.sin(Phi) - V*Omega*np.cos(Phi)
    return [Vx, Vy, Omega, Xtt, Ytt, Phitt]

from matplotlib.animation import FuncAnimation

def kadr(frame):
    global X,Y,Phi,VX,VY,Omega
    global beta, betas, beta_prev, rho, rho_prev, rhos, I_sred

    t.append(t[-1]+1)

    ## Если уже на последней
    #if near_index>=len(Route)-1:
    #    dx = Route[-1][0]-X
    #    dy = Route[-1][1]-Y
    #    if np.hypot(dx,dy)<stop_distance:
    #        print("Достигли цели, останавливаемся.")
    #        multik.event_source.stop()
    #        return [Robo]

    pA, pB = get_nearest_points(near_index)
    PhiT = -np.arctan2( pB[1]-pA[1], pB[0]-pA[0] )

    beta_prev=beta
    beta= np.arctan2( np.sin(PhiT-Phi), np.cos(PhiT-Phi))
    rho_prev=rho
    rho= get_distance(pA, pB, X, Y)

    betas=(beta-beta_prev)/dt
    rhos=(rho-rho_prev)/dt

    s_beta= sum(y1)/N_sred
    s_betas= sum(y2)/N_sred
    s_rho= sum(y3)/N_sred
    s_rhos= sum(y4)/N_sred

    dFb= 2*s_beta + 7*s_betas
    dFr= 6*s_rho + 10*s_rhos
    speed_now= np.hypot(VX,VY)
    F0m= F0 - 0.5*speed_now
    F1= F0m*(1 + dFb + dFr)
    F2= F0m*(1 - dFb - dFr)

    XK1, YK1, PhiK1, VXK1, VYK1, OmegaK1 = RobotSystemOfEquations(X,Y,Phi,VX,VY,Omega,F1,F2)
    XK2, YK2, PhiK2, VXK2, VYK2, OmegaK2 = RobotSystemOfEquations(X+dt/2*XK1,
                                                                 Y+dt/2*YK1,
                                                                 Phi+dt/2*PhiK1,
                                                                 VX+dt/2*VXK1,
                                                                 VY+dt/2*VYK1,
                                                                 Omega+dt/2*OmegaK1,
                                                                 F1,F2)
    XK3, YK3, PhiK3, VXK3, VYK3, OmegaK3 = RobotSystemOfEquations(X+dt/2*XK2,
                                                                 Y+dt/2*YK2,
                                                                 Phi+dt/2*PhiK2,
                                                                 VX+dt/2*VXK2,
                                                                 VY+dt/2*VYK2,
                                                                 Omega+dt/2*OmegaK2,
                                                                 F1,F2)
    XK4, YK4, PhiK4, VXK4, VYK4, OmegaK4 = RobotSystemOfEquations(X+dt*XK3,
                                                                 Y+dt*YK3,
                                                                 Phi+dt*PhiK3,
                                                                 VX+dt*VXK3,
                                                                 VY+dt*VYK3,
                                                                 Omega+dt*OmegaK3,
                                                                 F1,F2)

    X  += (XK1+2*XK2+2*XK3+XK4)*dt/6
    Y  += (YK1+2*YK2+2*YK3+YK4)*dt/6
    Phi+= (PhiK1+2*PhiK2+2*PhiK3+PhiK4)*dt/6
    VX += (VXK1+2*VXK2+2*VXK3+VXK4)*dt/6
    VY += (VYK1+2*VYK2+2*VYK3+VYK4)*dt/6
    Omega += (OmegaK1+2*OmegaK2+2*OmegaK3+OmegaK4)*dt/6

    y1[I_sred]=beta
    y2[I_sred]=betas
    y3[I_sred]=rho
    y4[I_sred]=rhos
    I_sred+=1
    if I_sred>=N_sred:
        I_sred=0
    # Окружность + луч => последнее значение = (0,0)
    # Поворачиваем на -Phi
    RX_rot, RY_rot = Rot2D(X_R, Y_R, -Phi)
    # Сдвигаем
    Robo.set_data( X*Razmer + RX_rot, Y*Razmer + RY_rot )

    if near_index >= len(Route) - 1:
        stop_distance = 1.0
    # Смотрим, на каком расстоянии робот от конечной точки
        dx = Route[-1][0] - X
        dy = Route[-1][1] - Y
        dist_to_goal = np.hypot(dx, dy)
        if dist_to_goal < stop_distance:
            print("Цель достигнута, останавливаем анимацию")
            multik.event_source.stop()
            return [Robo]
    return [Robo]

multik = FuncAnimation(fig, kadr, interval=dt*1, blit=True)
plt.show()
