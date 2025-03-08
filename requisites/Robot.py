import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from PIL import Image
import pickle


def IsColor(C1,C2):
    R = (C1-C2)/255
    if R[0]**2+R[1]**2+R[2]**2<0.5:
        return True
    else:
        return False


jpg = Image.open("Map1.jpg")
fig = plt.figure(figsize=[17,9])
ax = fig.add_subplot(1, 1, 1)
ax.imshow(jpg)

im2arr = np.array(jpg)
print(im2arr.shape)
print(im2arr[4])
Razmer = 63
N_x = im2arr.shape[1] // Razmer
N_y = im2arr.shape[0] // Razmer
Map = np.zeros((N_x, N_y))
ax.plot([N_x * Razmer, N_x * Razmer], [0, N_y * Razmer], color=[0, 1, 0])
ax.plot([0, N_x * Razmer], [N_y * Razmer, N_y * Razmer], color=[0, 1, 0])

f = 2

if f == 1:
    for i in range(N_x):
        ax.plot([i * Razmer, i * Razmer], [0, N_y * Razmer], color=[0, 1, 0])
        print(i, '/', N_x)
        for j in range(N_y):
            if i == 0:
                ax.plot([0, N_x * Razmer], [j * Razmer, j * Razmer], color=[0, 1, 0])
            N_p = 0
            N_s = 0
            N_c = 0
            for iz in range(Razmer):
                for jz in range(Razmer):
                    if IsColor(im2arr[j * Razmer + jz][i * Razmer + iz], [0, 0, 0]):
                        N_p += 1
                    elif IsColor(im2arr[j * Razmer + jz][i * Razmer + iz], [0, 0, 255]):
                        N_s += 1
                    elif IsColor(im2arr[j * Razmer + jz][i * Razmer + iz], [255, 0, 0]):
                        N_c += 1
            # print('ГђВЇГ‘вЂЎГђВµГђВ№ГђВєГђВ° [', i, j,']: ГђЕёГ‘в‚¬ГђВµГђВїГ‘ВЏГ‘вЂљГ‘ВЃГ‘вЂљГђВІГђВёГђВ№: ', N_p, 'ГђВЎГ‘вЂљГђВ°Г‘в‚¬Г‘вЂљ: ', N_s, 'ГђВ¦ГђВµГђВ»Г‘Е’: ', N_c)
            if N_s >= 10:
                Map[i, j] = -2
            if N_c >= 10:
                Map[i, j] = -3
            if N_p >= 10:
                Map[i, j] = -1
    fimeName = "Map1.map"
    dict = {'Map': Map}
    with open(fimeName,'wb') as file:
        pickle.dump(dict, file)

    print("РЎРґРµР»Р°РЅРѕ")
if f == 2:
    for i in range(N_x):
        ax.plot([i * Razmer, i * Razmer], [0, N_y * Razmer], color=[0, 1, 0])
        print(i, '/', N_x)
        for j in range(N_y):
            if i == 0:
                ax.plot([0, N_x * Razmer], [j * Razmer, j * Razmer], color=[0, 1, 0])
    Filename = "Map1.map"
    with open(Filename, 'rb') as file:
        Map = pickle.load(file)['Map']


I_s = 0
J_s = 0
I_c = 0
J_c = 0
N_s = 0
N_c = 0
for i in range(N_x):
    for j in range(N_y):
        if Map[i, j] == -2:
            N_s += 1
            I_s += i
            J_s += j
        if Map[i, j] == -3:
            N_c += 1
            I_c += i
            J_c += j

I_s = int(np.round(I_s/N_s))
J_s = int(np.round(J_s/N_s))

I_c = int(np.round(I_c/N_c))
J_c = int(np.round(J_c/N_c))

for i in range(N_x):
    for j in range(N_y):
        if Map[i, j] == -2 and (i!=I_s or j!=J_s):
            Map[i, j] = 0
        if Map[i, j] == -3 and (i!=I_c or j!=J_c):
            Map[i, j] = 0
        if Map[i,j] == 0:
            ax.plot((i+0.5)*Razmer,(j+0.5)*Razmer,'o',color=[1,1,1])
        if Map[i,j] == -1:
            ax.plot((i+0.5)*Razmer,(j+0.5)*Razmer,'o',color=[0,0,0])
        if Map[i,j] == -2:
            ax.plot((i+0.5)*Razmer,(j+0.5)*Razmer,'o',color=[0,0,1])
        if Map[i,j] == -3:
            ax.plot((i+0.5)*Razmer,(j+0.5)*Razmer,'o',color=[1,0,0])


def wave_algorithm(start, Map):
    steps=[ [( start[0],start[1] )]]
    InverseRoute = [[(start[0], start[1])]]
    i = 1
    fin = 0
    while len(steps[-1]) and fin==0:
        wave = []
        for point in steps[-1]:
            if point[0]+1 < Map.shape[0] and point[0]+1 >= 0:
                if Map[point[0]+1, point[1]] == 0:
                    wave.append((point[0]+1, point[1]))
                    Map[point[0] + 1, point[1]] = i
                elif Map[point[0]+1, point[1]] == -3:
                    wave.append((point[0] + 1, point[1]))
                    Map[point[0] + 1, point[1]] = i
                    fin = 1
                    InverseRoute = [[point[0] + 1, point[1]]]
                    break
            if point[0]-1 < Map.shape[0] and point[0]-1 >= 0:
                if Map[point[0] - 1, point[1]] == 0:
                    wave.append((point[0] - 1, point[1]))
                    Map[point[0] - 1, point[1]] = i
                elif Map[point[0] - 1, point[1]] == -3:
                    wave.append((point[0] - 1, point[1]))
                    Map[point[0] - 1, point[1]] = i
                    fin = 1
                    InverseRoute = [[point[0] - 1, point[1]]]
                    break
            if point[1] + 1 < Map.shape[1] and point[1] + 1 >= 0:
                if Map[point[0], point[1]+1] == 0:
                    wave.append((point[0], point[1]+1))
                    Map[point[0], point[1]+1] = i
                elif Map[point[0], point[1]+1] == -3:
                    wave.append((point[0], point[1]+1))
                    Map[point[0], point[1]+1] = i
                    fin = 1
                    InverseRoute = [[point[0], point[1]+1]]
                    break
            if point[1] - 1 < Map.shape[1] and point[1] - 1 >= 0:
                if Map[point[0], point[1]-1] == 0:
                    wave.append((point[0], point[1]-1))
                    Map[point[0], point[1]-1] = i
                elif Map[point[0], point[1]-1] == -3:
                    wave.append((point[0], point[1]-1))
                    Map[point[0], point[1]-1] = i
                    fin = 1
                    InverseRoute = [[point[0], point[1]-1]]
                    break
        steps.append(wave)
        i+=1
    print(Map)
    N_Steps = len(steps)
    print('ГђЕЎГђВѕГђВ»ГђВёГ‘вЂЎГђВµГ‘ВЃГ‘вЂљГђВІГђВѕ Г‘Л†ГђВ°ГђВіГђВѕГђВІ = ',N_Steps)
    if fin == 0:
        print('ГђВ¦ГђВµГђВ»Г‘Е’ ГђВЅГђВµГђВґГђВѕГ‘ВЃГ‘вЂљГђВёГђВ¶ГђВёГђВјГђВ°')
    else:
        for i in range(N_Steps):
            point = InverseRoute[i]
            print(point)
            print('ГђЛњГђВЅГђВґГђВµГђВєГ‘ВЃ ГђВєГђВ»ГђВµГ‘вЂљГђВєГђВё = ',Map[point[0],point[1]])
            Add = 0
            if point[0]+1 < Map.shape[0] and point[0]+1 >= 0 and Add == 0:
                if Map[point[0], point[1]] - Map[point[0]+1, point[1]] == 1:
                    InverseRoute.append([point[0]+1, point[1]])
                    Add = 1
                elif Map[point[0]+1, point[1]] == -2:
                    InverseRoute.append([point[0]+1, point[1]])
                    Add = 1
                    break
            if point[0]-1 < Map.shape[0] and point[0]-1 >= 0 and Add == 0:
                if Map[point[0], point[1]] - Map[point[0] - 1, point[1]] == 1:
                    InverseRoute.append([point[0] - 1, point[1]])
                    Add = 1
                elif Map[point[0] - 1, point[1]] == -2:
                    InverseRoute.append([point[0] - 1, point[1]])
                    Add = 1
                    break
            if point[1] + 1 < Map.shape[1] and point[1] + 1 >= 0 and Add == 0:
                if Map[point[0], point[1]] - Map[point[0], point[1]+1] == 1:
                    InverseRoute.append([point[0], point[1]+1])
                    Add = 1
                elif Map[point[0], point[1]+1] == -2:
                    InverseRoute.append([point[0], point[1]+1])
                    Add = 1
                    break
            if point[1] - 1 < Map.shape[1] and point[1] - 1 >= 0 and Add == 0:
                if Map[point[0], point[1]] - Map[point[0], point[1] - 1] == 1:
                    InverseRoute.append([point[0], point[1] - 1])
                    Add = 1
                elif Map[point[0], point[1] - 1] == -2:
                    InverseRoute.append([point[0], point[1] - 1])
                    Add = 1
                    break
            if Add == 0:
                print('ГђвЂњГђВѕГ‘вЂљГђВѕГђВІГђВѕ')
                break
        print('ГђЕѕГђВ±Г‘в‚¬ГђВ°Г‘вЂљГђВЅГ‘вЂ№ГђВ№ ГђВїГ‘Ж’Г‘вЂљГ‘Е’ = ', InverseRoute)
        Route = InverseRoute[::-1]
        print('ГђЕёГ‘в‚¬Г‘ВЏГђВјГђВѕГђВ№ ГђВїГ‘Ж’Г‘вЂљГ‘Е’ = ', Route)
    return Route
Route = wave_algorithm(np.array([I_s, J_s]), Map)
for point in Route:
    if Map[point[0],point[1]] != -2:
        ax.plot((point[0] + 0.5) * Razmer, (point[1] + 0.5) * Razmer, 'o', color=[1*Map[point[0],point[1]]/(len(Route)), 0, 1*(len(Route)-Map[point[0],point[1]])/(len(Route))])

Route.insert(0, [1, 1])
Route = [[i[0] + 0.5, i[1] + 0.5] for i in Route]

near_index = 0

m = 1
R = 10
J = m * R**2 / 2
mu = 0.2
g = 9.81

X = (I_s+0.5)
Y = (J_s+0.5)
Phi = -np.arctan2(Route[1][1]-Route[0][1],Route[1][0]-Route[0][0])
VX = 0
VY = 0
Omega = 0

epsilon = 0.05


def get_nearest_points(index):
    global near_index
    if not index == len(Route)-1:
        if (Route[index][0] - X)**2+(Route[index][1] - Y)**2 > (Route[index+1][0] - X)**2 + (Route[index+1][1] - Y)**2:
            index += 1
            near_index += 1
    return Route[index], Route[index + 1]


def get_distance(near_point, next_point, X, Y):
    A = near_point[1] - next_point[1]
    B = next_point[0] - near_point[0]
    C = -next_point[0]*near_point[1] + next_point[0]*next_point[1] + near_point[0]*next_point[1] - next_point[0]*next_point[1]
    d = (A*X + B*Y + C) / np.sqrt(A**2 + B**2)
    return d


def RobotSystemOfEquations(X, Y, Phi, Vx, Vy, Omega, F1, F2):

    V = Vx * np.cos(Phi) - Vy * np.sin(Phi)
    Ftr = 2 / np.pi * np.arctan(10000 * V) * mu * m * g
    Vst = (F1 + F2 - Ftr) / m
    Phitt = (F1 - F2) * R / J

    Xtt = Vst * np.cos(Phi) - V * Omega * np.sin(Phi)
    Ytt = -Vst * np.sin(Phi) - V * Omega * np.cos(Phi)

    return [Vx, Vy, Omega, Xtt, Ytt, Phitt]


# fig = plt.figure()
dt = 0.01
F0 = 0.7*mu*m*g
alpha_R = np.linspace(0,2*np.pi,51)
X_R = np.append(Razmer/2*np.cos(alpha_R),0)
Y_R = np.append(Razmer/2*np.sin(alpha_R),0)
Robo = ax.plot(X_R, Y_R, color=[1, 0, 1])[0]

beta = 0
beta_prev = 0
betas = 0
rho = 0
rho_prev = 0
rhos = 0
K1, K2, K3, K4, K5 = 2, 7, 6, 10, 0.5

# figg = plt.figure(figsize=(16, 9))
#ax1 = fig.add_subplot(2,1,2)
# ax2 = figg.add_subplot(2,2,2)
# ax3 = figg.add_subplot(2,2,3)
# ax4 = figg.add_subplot(2,2,4)
t = [0]
N_sred = 50
y1, y2, y3, y4 = np.zeros(N_sred), np.zeros(N_sred), np.zeros(N_sred), np.zeros(N_sred)
# rhog = ax1.plot(t, y1, color=[1, 0, 0])[0]
# betag = ax1.plot(t, y2, color=[0,1,0])[0]
# rhosg = ax1.plot(t, y3, color=[0,0,1])[0]
# betasg = ax1.plot(t, y4,color=[0,1,1])[0]


def Rot2D(X, Y, Alpha):
    RX = X * np.cos(Alpha) - Y * np.sin(Alpha)
    RY = X * np.sin(Alpha) + Y * np.cos(Alpha)
    return RX, RY

I_sred = 0
def kadr(i):
    global X, Y, Phi, VX, VY, Omega, beta, betas, beta_prev, rho, rho_prev, rhos, x, y, I_sred
    t.append(t[-1]+1)


    near_point, next_point = get_nearest_points(near_index)
    PhiT = -np.arctan2(next_point[1] - near_point[1], next_point[0] - near_point[0])
    # print(near_point, next_point)
    # print(PhiT)
    beta_prev = beta
    beta = np.arctan2(np.sin(PhiT - Phi), np.cos(PhiT - Phi))
    print(beta)
    rho_prev = rho
    rho = get_distance(near_point, next_point, X, Y)
    print(rho)

    betas = (beta - beta_prev) / dt
    rhos = (rho - rho_prev) / dt

    s_beta = (sum(y1)/N_sred)
    s_betas = (sum(y2)/N_sred)
    s_rho = (sum(y3) / N_sred)
    s_rhos = (sum(y4) / N_sred)

    print(s_beta)
    print(s_rho)
    dFb = K1 * s_beta + K2 * s_betas
    dFr = K3 * s_rho + K4 * s_rhos
    F0m = F0 - K5*(VX**2+VY**2)**0.5
    F1 = F0m*(1 + dFb + dFr)
    F2 = F0m*(1 - dFb - dFr)
    print('____________________________')

    XK1, YK1, PhiK1, VXK1, VYK1, OmegaK1 = RobotSystemOfEquations(X, Y, Phi, VX, VY, Omega, F1, F2)
    # XK1, YK1, PhiK1, VXK1, VYK1, OmegaK1 = np.array(DX), np.array(DY), np.array(DPhi), np.array(DVX), np.array(DVY), np.array(DOmega)
    XK2, YK2, PhiK2, VXK2, VYK2, OmegaK2 = RobotSystemOfEquations(X+dt/2*XK1, Y+dt/2*YK1, Phi+dt/2*PhiK1, VX+dt/2*VXK1, VY+dt/2*VYK1, Omega+dt/2*OmegaK1, F1, F2)
    # XK2, YK2, PhiK2, VXK2, VYK2, OmegaK2 = np.array(DX), np.array(DY), np.array(DPhi), np.array(DVX), np.array(DVY), np.array(DOmega)
    XK3, YK3, PhiK3, VXK3, VYK3, OmegaK3 = RobotSystemOfEquations(X+dt/2*XK2, Y+dt/2*YK2, Phi+dt/2*PhiK2, VX+dt/2*VXK2, VY+dt/2*VYK2, Omega+dt/2*OmegaK2, F1, F2)
    # XK3, YK3, PhiK3, VXK3, VYK3, OmegaK3 = np.array(DX), np.array(DY), np.array(DPhi), np.array(DVX), np.array(DVY), np.array(DOmega)
    XK4, YK4, PhiK4, VXK4, VYK4, OmegaK4 = RobotSystemOfEquations(X+dt/2*XK3, Y+dt/2*YK3, Phi+dt/2*PhiK3, VX+dt/2*VXK3, VY+dt/2*VYK3, Omega+dt/2*OmegaK3, F1, F2)
    # XK4, YK4, PhiK4, VXK4, VYK4, OmegaK4 = np.array(DX), np.array(DY), np.array(DPhi), np.array(DVX), np.array(DVY), np.array(DOmega)

    X += (XK1+2*XK2+2*XK3+XK4) * dt/6
    Y += (YK1+2*YK2+2*YK3+YK4) * dt/6
    Phi += (PhiK1+2*PhiK2+2*PhiK3+PhiK4) * dt/6
    VX += (VXK1+2*VXK2+2*VXK3+VXK4) * dt/6
    VY += (VYK1+2*VYK2+2*VYK3+VYK4) * dt/6
    Omega += (OmegaK1+2*OmegaK2+2*OmegaK3+OmegaK4) * dt/6
    RX_R, RY_R = Rot2D(X_R, Y_R, -Phi)

    y3[I_sred] = rho
    y4[I_sred] = rhos
    y1[I_sred] = beta
    y2[I_sred] = betas
    I_sred+=1
    if I_sred>=N_sred:
        I_sred -= N_sred
    # rhog.set_data(t, y1)
    #
    # y2.append(beta)
    # betag.set_data(t, y2)
    #
    # y3.append(rhos)
    # rhosg.set_data(t, y3)
    #
    # y4.append(betas)
    # betasg.set_data(t, y4)
    # ax1.set_xlim(0, t[-1])
    # ax1.set_ylim(min(min(y1), min(y2), min(y3), min(y4)), max(max(y1), max(y2), max(y3), max(y4)))
    # ax2.set_ydata(y)

    Robo.set_data(X*Razmer + RX_R, Y*Razmer + RY_R)
    # print(i,X,Y,Phi)
    return [Robo]


multik = FuncAnimation(fig, kadr, interval=dt*1, blit=True)

plt.show()