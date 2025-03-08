import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.image as mpimg
from PIL import Image

m = 1
L = 0.25
k_A = 0.1
k_B = 0.2
k_C = 0.4
A = m*L**2*k_A
B = m*L**2*k_B
C = m*L**2*k_C
g = 9.81

X = 0
Y = 0
Z = 0

q0 = 1
q1 = 0
q2 = 0
q3 = 0

VX = 0
VY = 0
VZ = 0

Omegax = 0
Omegay = 0
Omegaz = 0

k_t = 1
k_m = 0.2
R = 0.1
S_v = np.pi*R**2
rho = 1.2
Omega0 = np.sqrt(m*g / (2 * k_t * S_v * rho * R ** 2))
print(Omega0)

Cx = 0.3
Cy = 0.3
Cz = 0.9
Cvx = 0.6
Cvy = 0.6
Cvz = 0.3
Sx = 2*L*L/2*0.4
Sy = 2*L*L/2*0.4
Sz = 2*L*2*L*0.4
Svz = np.sqrt(Sx**2 + Sy**2)

k_Z = 2
k_VZ = 2

x_VC = np.array([L, -L, -L, L])
y_VC = np.array([L, L, -L, -L])
z_VC = np.array([0, 0, 0, 0])
rC = np.array([0,1,0])
TC = np.array([0,1,0])

k_XY = 0.5
k_dXY = 0.7

k_N = 5
k_dN = 1

k_R = 0.14
k_dR = 0.2

O_max = 4
O_min = -1


# Р’Р°СЂРёР°РЅС‚ 1: РџРѕСЃС‚РѕСЏРЅРЅС‹Р№ РІРµС‚РµСЂ
def wind_constant(X, Y, Z):
    VwindX = 10
    VwindY = 0
    VwindZ = 0
    return VwindX, VwindY, VwindZ


# Р’Р°СЂРёР°РЅС‚ 2: Р Р°РІРЅРѕРјРµСЂРЅРѕ РјРµРЅСЏСЋС‰РёР№СЃСЏ РІРµС‚РµСЂ
def wind_linear(X, Y, Z):
    ax, ay, az = 0.1, 0.1, 0.05
    VwindX = ax * X
    VwindY = ay * Y
    VwindZ = az * Z
    return VwindX, VwindY, VwindZ


# Р’Р°СЂРёР°РЅС‚ 3: РЎРєР°С‡РєРѕРѕР±СЂР°Р·РЅРѕ РјРµРЅСЏСЋС‰РёР№СЃСЏ РІРµС‚РµСЂ
def wind_stepwise(X, Y, Z):
    x0, y0, z0 = 4, 4, 4  # Р“СЂР°РЅРёС†С‹ РїРµСЂРµРєР»СЋС‡РµРЅРёР№
    VwindX = 2 if X < x0 else -2
    VwindY = 1 if Y < y0 else -1
    VwindZ = 0.5 if Z < z0 else -0.5
    return VwindX, VwindY, VwindZ



def CopterSystemOfEquations(X, Y, Z, q0, q1, q2, q3, VX, VY, VZ, Omegax, Omegay, Omegaz):
    global Omega1, Omega2, Omega3, Omega4, PhiZ, \
        CelX, CelY, CelZ, CelVX, CelVY, CelVZ
    dZ = - k_Z*np.arctan(Z - CelZ) - k_VZ*(VZ - CelVZ)

    X_VC, Y_VC, Z_VC = Rot3D(x_VC, y_VC, z_VC, q0, q1, q2, q3)
    Sm = np.arctan(((X - CelX)**2+(Y - CelY)**2)**0.5)/((X - CelX)**2+(Y - CelY)**2)**0.5

    dXY1 = k_XY * ((X - CelX) * X_VC[0] + (Y - CelY) * Y_VC[0]) * Sm + k_dXY * ((VX - CelVX) * X_VC[0] + (VY - CelVY) * Y_VC[0])
    dXY2 = k_XY * ((X - CelX) * X_VC[1] + (Y - CelY) * Y_VC[1]) * Sm + k_dXY * ((VX - CelVX) * X_VC[1] + (VY - CelVY) * Y_VC[1])
    dXY3 = k_XY * ((X - CelX) * X_VC[2] + (Y - CelY) * Y_VC[2]) * Sm + k_dXY * ((VX - CelVX) * X_VC[2] + (VY - CelVY) * Y_VC[2])
    dXY4 = k_XY * ((X - CelX) * X_VC[3] + (Y - CelY) * Y_VC[3]) * Sm + k_dXY * ((VX - CelVX) * X_VC[3] + (VY - CelVY) * Y_VC[3])

    OmX_VC, OmY_VC, OmZ_VC = Rot3D(Omegay*z_VC-Omegaz*y_VC, Omegaz*x_VC-Omegax*z_VC, Omegax*y_VC-Omegay*x_VC, q0, q1, q2, q3)
    dX_VC, dY_VC, dZ_VC = OmX_VC, OmY_VC, OmZ_VC
    dN1 = - k_N * Z_VC[0] - k_dN * dZ_VC[0]
    dN2 = - k_N * Z_VC[1] - k_dN * dZ_VC[1]
    dN3 = - k_N * Z_VC[2] - k_dN * dZ_VC[2]
    dN4 = - k_N * Z_VC[3] - k_dN * dZ_VC[3]

    TC = np.array([CelVX, CelVY, 0])/(CelVX**2 + CelVY**2)**0.5
    RC = Rot3D(rC[0], rC[1], rC[2], q0, q1, q2, q3)
    Cos = RC[0]*TC[0] + RC[1]*TC[1]
    Sin = RC[1]*TC[0] - RC[0]*TC[1]
    PhiZ = np.arctan2(Sin,Cos)
    dR1 = -k_R*PhiZ - k_dR*Omegaz
    dR2 = k_R*PhiZ + k_dR*Omegaz
    dR3 = -k_R*PhiZ - k_dR*Omegaz
    dR4 = k_R*PhiZ + k_dR*Omegaz
    print('PhiZ = ', PhiZ, 'dRs = ', dR1,dR2,dR3,dR4)

    ROmega1 = Omega0*(1 + dZ + dXY1 + dN1 + dR1)
    ROmega2 = Omega0*(1 + dZ + dXY2 + dN2 + dR2)
    ROmega3 = Omega0*(1 + dZ + dXY3 + dN3 + dR3)
    ROmega4 = Omega0*(1 + dZ + dXY4 + dN4 + dR4)

    KO1B,KO1M,KO2B,KO2M,KO3B,KO3M,KO4B,KO4M = 1,1,1,1,1,1,1,1
    if ROmega1 > Omega0*O_max:
        KO1B = ROmega1 / (Omega0 * O_max)
    elif ROmega1 < Omega0*O_min:
        KO1M = -ROmega1 / (Omega0 * O_min)
    if ROmega2 > Omega0*O_max:
        KO2B = ROmega2 / (Omega0 * O_max)
    elif ROmega2 < Omega0*O_min:
        KO2M = -ROmega2 / (Omega0 * O_min)
    if ROmega3 > Omega0*O_max:
        KO3B = ROmega3 / (Omega0 * O_max)
    elif ROmega3 < Omega0*O_min:
        KO3M = -ROmega3 / (Omega0 * O_min)
    if ROmega4 > Omega0*O_max:
        KO4B = ROmega4 / (Omega0 * O_max)
    elif ROmega4 < Omega0*O_min:
        KO4M = -ROmega4 / (Omega0 * O_min)

    KO = max(KO1B,KO1M,KO2B,KO2M,KO3B,KO3M,KO4B,KO4M)
    Omega1 = ROmega1 / KO
    Omega2 = ROmega2 / KO
    Omega3 = ROmega3 / KO
    Omega4 = ROmega4 / KO

    print('Sm = ', Sm, ' Sz = ', np.arctan(Z) / Z, ' KO = ', KO)

    print('Omega1 = ', Omega1, 'Omega2 = ', Omega2, 'Omega3 = ', Omega3, 'Omega4 = ', Omega4)

    F1 = k_t * S_v * rho * R ** 2 * Omega1 * abs(Omega1) / 2
    M1 = k_m * S_v * rho * R ** 3 * Omega1 * abs(Omega1) / 2
    F2 = k_t * S_v * rho * R ** 2 * Omega2 * abs(Omega2) / 2
    M2 = k_m * S_v * rho * R ** 3 * Omega2 * abs(Omega2) / 2
    F3 = k_t * S_v * rho * R ** 2 * Omega3 * abs(Omega3) / 2
    M3 = k_m * S_v * rho * R ** 3 * Omega3 * abs(Omega3) / 2
    F4 = k_t * S_v * rho * R ** 2 * Omega4 * abs(Omega4) / 2
    M4 = k_m * S_v * rho * R ** 3 * Omega4 * abs(Omega4) / 2

    VwindX, VwindY, VwindZ = wind_constant(X, Y, Z)
    #VwindX, VwindY, VwindZ = wind_linear(X, Y, Z)
    #VwindX, VwindY, VwindZ = wind_stepwise(X, Y, Z)

    F_VX, F_VY, F_VZ = Rot3D(0, 0, F1+F2+F3+F4, q0, q1, q2, q3)
    Vx, Vy, Vz = TRot3D(VX - VwindX, VY - VwindY, VZ - VwindZ, q0, q1, q2, q3)
    FSx = rho * Cx * Sx * Vx * abs(Vx) / 2
    FSy = rho * Cy * Sy * Vy * abs(Vy) / 2
    FSz = rho * Cz * Sz * Vz * abs(Vz) / 2
    FSX, FSY, FSZ = Rot3D(FSx, FSy, FSz, q0, q1, q2, q3)

    MSx = rho * Cvx * Sz * Omegax * abs(Omegax) / 2
    MSy = rho * Cvy * Sz * Omegay * abs(Omegay) / 2
    MSz = rho * Cvz * Svz * Omegaz * abs(Omegaz) / 2

    FX = F_VX - FSX
    FY = F_VY - FSY
    FZ = - m * g + F_VZ - FSZ

    MCx = L*F1 + L*F2 - L*F3 - L*F4 - MSx
    MCy = -L*F1 + L*F2 + L*F3 - L*F4 - MSy
    MCz = M1 - M2 + M3 - M4 - MSz

    dX = VX
    dY = VY
    dZ = VZ
    dq0 = 0.5*(- q1*Omegax - q2*Omegay - q3*Omegaz)
    dq1 = 0.5*(q0*Omegax + q2*Omegaz - q3*Omegay)
    dq2 = 0.5*(q0*Omegay - q1*Omegaz + q3*Omegax)
    dq3 = 0.5*(q0*Omegaz + q1*Omegay - q2*Omegax)
    dVX = FX/m
    dVY = FY/m
    dVZ = FZ/m
    dOmegax = (MCx - (C-B)*Omegay*Omegaz)/A
    dOmegay = (MCy - (A-C)*Omegaz*Omegax)/B
    dOmegaz = (MCz - (B-A)*Omegax*Omegay)/C

    return [dX, dY, dZ, dq0, dq1, dq2, dq3, dVX, dVY, dVZ, dOmegax, dOmegay, dOmegaz]


dt = 0.008

fig = plt.figure(figsize=(19, 9))
ax = fig.add_subplot(1,2,1,projection='3d')
axg = fig.add_subplot(2,4,3)
axo = fig.add_subplot(2,4,4)
axv = fig.add_subplot(2,2,4)
#ax.axis('equal')
ax.set(xlim=[-20*L, 20*L], ylim=[-20*L, 20*L], zlim=[-20*L, 20*L])
axg.set(xlim=[0, 10], ylim=[-600, 600])
axo.set(xlim=[0, 10], ylim=[-1, 1])
axv.axis('equal')
axv.set(xlim=[-0.65, 0.65], ylim=[-0.65, 0.65])

X_Copter = np.array([0, L, 0, -L, 0, -L, 0, L, 0])
Y_Copter = np.array([0, L, 0, L, 0, -L, 0, -L, 0])
Z_Copter = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0])
K_Copter = ax.plot(X + X_Copter, Y + Y_Copter, Z + Z_Copter, color=[0, 0, 1])[0]

C_Copter = ax.plot(X, Y, Z, 'o', color=[1, 0, 0])[0]
A_Copter = ax.plot(X + X_Copter[1], Y + Y_Copter[1], Z + Z_Copter[1],'o', color=[0, 1, 0])[0]
B_Copter = ax.plot(X + X_Copter[3], Y + Y_Copter[3], Z + Z_Copter[3],'o', color=[1, 0, 1])[0]
Ce_Copter = ax.plot(X + X_Copter[5], Y + Y_Copter[5], Z + Z_Copter[5],'o', color=[0, 1, 1])[0]
D_Copter = ax.plot(X + X_Copter[7], Y + Y_Copter[7], Z + Z_Copter[7],'o', color=[1, 1, 0])[0]


def Rot3D(x,y,z, q0,q1,q2,q3):
    RX = x * (q0**2 + q1**2 - q2**2 - q3**2) + y * (- 2*q0*q3 + 2*q1*q2) + z * (2*q1*q3 + 2*q0*q2)
    RY = x * (2*q0*q3 + 2*q1*q2) + y * (q0**2 - q1**2 + q2**2 - q3**2) + z * (- 2*q0*q1 + 2*q2*q3)
    RZ = x * (2*q1*q3 - 2*q0*q2) + y * (2*q0*q1 + 2*q2*q3) + z * (q0**2 - q1**2 - q2**2 + q3**2)
    return RX, RY, RZ

def TRot3D(X,Y,Z, q0,q1,q2,q3):
    Rx = X * (q0**2 + q1**2 - q2**2 - q3**2) + Y * (2*q0*q3 + 2*q1*q2) + Z * (2*q1*q3 - 2*q0*q2)
    Ry = X * (- 2*q0*q3 + 2*q1*q2) + Y * (q0**2 - q1**2 + q2**2 - q3**2) + Z * (2*q0*q1 + 2*q2*q3)
    Rz = X * (2*q1*q3 + 2*q0*q2) + Y * (- 2*q0*q1 + 2*q2*q3) + Z * (q0**2 - q1**2 - q2**2 + q3**2)
    return Rx, Ry, Rz

MOmega1 = np.array([])
MOmega2 = np.array([])
MOmega3 = np.array([])
MOmega4 = np.array([])
MT = np.array([])
Omega1, Omega2, Omega3, Omega4 = 0, 0, 0, 0
Pomega1 = axg.plot(MT,MOmega1)[0]
Pomega2 = axg.plot(MT,MOmega2)[0]
Pomega3 = axg.plot(MT,MOmega3)[0]
Pomega4 = axg.plot(MT,MOmega4)[0]

MOX = np.array([])
MOY = np.array([])
MOZ = np.array([])
MOVX = np.array([])
MOVY = np.array([])
MOVZ = np.array([])
MOmegax = np.array([])
MOmegay = np.array([])
MOmegaz = np.array([])
MPhiZ = np.array([])

OX = axo.plot(MT,MOX)[0]
OY = axo.plot(MT,MOY)[0]
OZ = axo.plot(MT,MOZ)[0]
OVX = axo.plot(MT, MOVX)[0]
OVY = axo.plot(MT, MOVY)[0]
OVZ = axo.plot(MT, MOVZ)[0]
Pomegax = axo.plot(MT,MOmegax)[0]
Pomegay = axo.plot(MT,MOmegay)[0]
Pomegaz = axo.plot(MT,MOmegaz)[0]
PPhiZ = axo.plot(MT,MPhiZ)[0]
"""
Nlines = 51
ShLines = 2
ZLines = -30
LPoints = np.linspace(-(Nlines-1)/2*ShLines,(Nlines-1)/2*ShLines,Nlines)
XPoints,YPoints = np.meshgrid(LPoints,LPoints)
XPoints,YPoints = XPoints.ravel(),YPoints.ravel()
ZPoints = XPoints*0 + ZLines
print('XPts = ', XPoints, 'YPts = ', YPoints, 'ZPts = ', ZPoints)

image = Image.open('padro.png')
image_array = np.array(image)  # РџСЂРµРѕР±СЂР°Р·СѓРµРј РІ РјР°СЃСЃРёРІ
image_height, image_width, _ = image_array.shape  # СѓР·РЅР°РµРј СЂР°Р·РјРµСЂС‹

def GetAngles(X, Y, Z, XPoints,YPoints,ZPoints, q0, q1, q2, q3):
    Norms = ((XPoints - X) ** 2 + (YPoints - Y) ** 2 + (ZPoints - Z) ** 2) ** 0.5
    VXP, VYP, VZP = (XPoints - X)/Norms, (YPoints - Y)/Norms, (ZPoints - Z)/Norms
    TVXP, TVYP, TVZP = TRot3D(VXP, VYP, VZP, q0, q1, q2, q3)
    #Psi = np.arcsin(VXP)
    #Phi = np.arcsin(VYP/np.cos(Psi))
    return TVXP,TVYP#Psi,Phi
VXP,VYP = GetAngles(X, Y, Z, XPoints,YPoints,ZPoints, q0, q1, q2, q3)
PPoints = axv.plot(VXP,VYP,'o')[0]
axv.plot([-0.65, 0.65, 0.65, -0.65, -0.65],[-0.65, -0.65, 0.65, 0.65, -0.65],'r')"""

Nlines = 101
ShLines = 0.5
ZLines = -30

LPoints = np.linspace(-(Nlines - 1) / 2 * ShLines, (Nlines - 1) / 2 * ShLines, Nlines)
XPoints, YPoints = np.meshgrid(LPoints, LPoints)
XPoints, YPoints = XPoints.ravel(), YPoints.ravel()
ZPoints = XPoints * 0 + ZLines
print('XPts = ', XPoints, 'YPts = ', YPoints, 'ZPts = ', ZPoints)

#image = Image.open('padro.png')
image = Image.open('karta.jpg')
image_array = np.array(image)
image_height, image_width, _ = image_array.shape

def GetAngles(X, Y, Z, XPoints, YPoints, ZPoints, q0, q1, q2, q3, q0_cam, q1_cam, q2_cam, q3_cam):
    Norms = ((XPoints - X) ** 2 + (YPoints - Y) ** 2 + (ZPoints - Z) ** 2) ** 0.5
    VXP, VYP, VZP = (XPoints - X) / Norms, (YPoints - Y) / Norms, (ZPoints - Z) / Norms
    TVXP, TVYP, TVZP = TRot3D(VXP, VYP, VZP, q0, q1, q2, q3)
    TVXP, TVYP, TVZP = TRot3D(TVXP, TVYP, TVZP, q0_cam, q1_cam, q2_cam, q3_cam)
    return TVXP, VYP
Alpha_cam = 1
q0_cam = np.cos(Alpha_cam/2)
q1_cam = np.sin(Alpha_cam/2)
q2_cam = 0
q3_cam = 0

VXP,VYP = GetAngles(X, Y, Z, XPoints,YPoints,ZPoints, q0, q1, q2, q3, q0_cam, q1_cam, q2_cam, q3_cam)

VXP_norm = (VXP + 1) / 2
VYP_norm = (VYP + 1) / 2

u_pix = (VXP * image_width).astype(int)
v_pix = (VYP * image_height).astype(int)

u_pix = np.clip(u_pix, 0, image_width - 1)
v_pix = np.clip(v_pix, 0, image_height - 1)

colors = image_array[v_pix, u_pix]

PPoints = axv.scatter(VXP, VYP, c=colors/255, s=10)

axv.plot([-0.65, 0.65, 0.65, -0.65, -0.65], [-0.65, -0.65, 0.65, 0.65, -0.65], 'r')

CelX = 0
CelY = 5
CelZ = 0
CelVX = 5*0.4
CelVY = 0
CelVZ = 2.5*0.8

X, Y, Z = CelX, CelY+0.1, CelZ+0.1

PCel = ax.plot(CelX, CelY, CelZ, 'o', color=[0.7, 0, 0.5])[0]


def kadr(i):
    global X, Y, Z, q0, q1, q2, q3, VX, VY, VZ, Omegax, Omegay, Omegaz, \
        Omega1, Omega2, Omega3, Omega4, PhiZ, \
        MOmega1, MOmega2, MOmega3, MOmega4, MOX, MOY, MOZ, \
        MOVX, MOVY, MOVZ, MOmegax, MOmegay, MOmegaz, MPhiZ, MT, \
        CelX, CelY, CelZ, CelVX, CelVY, CelVZ, q0_cam, q1_cam, q2_cam, q3_cam

    CelX = 5 * np.sin(0.4 * i * dt)
    CelY = 5 * np.cos(0.4 * i * dt)
    CelZ = 2.5 * np.sin(0.8 * i * dt)
    CelVX = 5 * 0.4 * np.cos(0.4 * i * dt)
    CelVY = -5 * 0.4 * np.sin(0.4 * i * dt)
    CelVZ = 2.5 * 0.8 * np.cos(0.8 * i * dt)

    XK1, YK1, ZK1, q0K1, q1K1, q2K1, q3K1, VXK1, VYK1, VZK1, OmegaxK1, OmegayK1, OmegazK1 = \
        CopterSystemOfEquations(
            X, Y, Z, q0, q1, q2, q3,
            VX, VY, VZ, Omegax, Omegay, Omegaz)
    Omega1K1, Omega2K1, Omega3K1, Omega4K1 = Omega1, Omega2, Omega3, Omega4

    XK2, YK2, ZK2, q0K2, q1K2, q2K2, q3K2, VXK2, VYK2, VZK2, OmegaxK2, OmegayK2, OmegazK2 = \
        CopterSystemOfEquations(
            X + dt / 2 * XK1, Y + dt / 2 * YK1, Z + dt / 2 * ZK1,
            q0 + dt / 2 * q0K1, q1 + dt / 2 * q1K1, q2 + dt / 2 * q2K1, q3 + dt / 2 * q3K1,
            VX + dt / 2 * VXK1, VY + dt / 2 * VYK1, VZ + dt / 2 * VZK1,
            Omegax + dt / 2 * OmegaxK1, Omegay + dt / 2 * OmegayK1, Omegaz + dt / 2 * OmegazK1)
    Omega1K2, Omega2K2, Omega3K2, Omega4K2 = Omega1, Omega2, Omega3, Omega4

    XK3, YK3, ZK3, q0K3, q1K3, q2K3, q3K3, VXK3, VYK3, VZK3, OmegaxK3, OmegayK3, OmegazK3 = \
        CopterSystemOfEquations(
            X + dt / 2 * XK2, Y + dt / 2 * YK2, Z + dt / 2 * ZK2,
            q0 + dt / 2 * q0K2, q1 + dt / 2 * q1K2, q2 + dt / 2 * q2K2, q3 + dt / 2 * q3K2,
            VX + dt / 2 * VXK2, VY + dt / 2 * VYK2, VZ + dt / 2 * VZK2,
            Omegax + dt / 2 * OmegaxK2, Omegay + dt / 2 * OmegayK2, Omegaz + dt / 2 * OmegazK2)
    Omega1K3, Omega2K3, Omega3K3, Omega4K3 = Omega1, Omega2, Omega3, Omega4

    XK4, YK4, ZK4, q0K4, q1K4, q2K4, q3K4, VXK4, VYK4, VZK4, OmegaxK4, OmegayK4, OmegazK4 = \
        CopterSystemOfEquations(
            X + dt * XK3, Y + dt * YK3, Z + dt * ZK3,
            q0 + dt * q0K3, q1 + dt * q1K3, q2 + dt * q2K3, q3 + dt * q3K3,
            VX + dt * VXK3, VY + dt * VYK3, VZ + dt * VZK3,
            Omegax + dt * OmegaxK3, Omegay + dt * OmegayK3, Omegaz + dt * OmegazK3)
    Omega1K4, Omega2K4, Omega3K4, Omega4K4 = Omega1, Omega2, Omega3, Omega4

    X += (XK1 + 2 * XK2 + 2 * XK3 + XK4) * dt / 6
    Y += (YK1 + 2 * YK2 + 2 * YK3 + YK4) * dt / 6
    Z += (ZK1 + 2 * ZK2 + 2 * ZK3 + ZK4) * dt / 6

    q0 += (q0K1 + 2 * q0K2 + 2 * q0K3 + q0K4) * dt / 6
    q1 += (q1K1 + 2 * q1K2 + 2 * q1K3 + q1K4) * dt / 6
    q2 += (q2K1 + 2 * q2K2 + 2 * q2K3 + q2K4) * dt / 6
    q3 += (q3K1 + 2 * q3K2 + 2 * q3K3 + q3K4) * dt / 6

    VX += (VXK1 + 2 * VXK2 + 2 * VXK3 + VXK4) * dt / 6
    VY += (VYK1 + 2 * VYK2 + 2 * VYK3 + VYK4) * dt / 6
    VZ += (VZK1 + 2 * VZK2 + 2 * VZK3 + VZK4) * dt / 6

    Omegax += (OmegaxK1 + 2 * OmegaxK2 + 2 * OmegaxK3 + OmegaxK4) * dt / 6
    Omegay += (OmegayK1 + 2 * OmegayK2 + 2 * OmegayK3 + OmegayK4) * dt / 6
    Omegaz += (OmegazK1 + 2 * OmegazK2 + 2 * OmegazK3 + OmegazK4) * dt / 6

    Omega1 = (Omega1K1 + 2 * Omega1K2 + 2 * Omega1K3 + Omega1K4) / 6
    Omega2 = (Omega2K1 + 2 * Omega2K2 + 2 * Omega2K3 + Omega2K4) / 6
    Omega3 = (Omega3K1 + 2 * Omega3K2 + 2 * Omega3K3 + Omega3K4) / 6
    Omega4 = (Omega4K1 + 2 * Omega4K2 + 2 * Omega4K3 + Omega4K4) / 6

    MOmega1 = np.append(MOmega1, Omega1)
    MOmega2 = np.append(MOmega2, Omega2)
    MOmega3 = np.append(MOmega3, Omega3)
    MOmega4 = np.append(MOmega4, Omega4)

    MOX = np.append(MOX, X - CelX)
    MOY = np.append(MOY, Y - CelY)
    MOZ = np.append(MOZ, Z - CelZ)
    MOVX = np.append(MOVX, VX - CelVX)
    MOVY = np.append(MOVY, VY - CelVY)
    MOVZ = np.append(MOVZ, VZ - CelVZ)
    MOmegax = np.append(MOmegax, Omegax)
    MOmegay = np.append(MOmegay, Omegay)
    MOmegaz = np.append(MOmegaz, Omegaz)
    MPhiZ = np.append(MPhiZ, PhiZ)

    MT = np.append(MT, i * dt)

    print('-----------------  ', i, '  -----------------')
    print('X = ', X, 'Y = ', Y, 'Z = ', Z)
    print('VX = ', VX, 'VY = ', VY, 'VZ = ', VZ)
    print('Omegax = ', Omegax, ', Omegay = ', Omegay, ', Omegaz = ', Omegaz)



    C_Copter.set_data_3d([X], [Y], [Z])
    RX_Copter, RY_Copter, RZ_Copter = Rot3D(X_Copter, Y_Copter, Z_Copter, q0, q1, q2, q3)

    K_Copter.set_data_3d(X + RX_Copter, Y + RY_Copter, Z + RZ_Copter)
    A_Copter.set_data_3d([X + RX_Copter[1]], [Y + RY_Copter[1]], [Z + RZ_Copter[1]])
    B_Copter.set_data_3d([X + RX_Copter[3]], [Y + RY_Copter[3]], [Z + RZ_Copter[3]])
    Ce_Copter.set_data_3d([X + RX_Copter[5]], [Y + RY_Copter[5]], [Z + RZ_Copter[5]])
    D_Copter.set_data_3d([X + RX_Copter[7]], [Y + RY_Copter[7]], [Z + RZ_Copter[7]])
    PCel.set_data_3d([CelX], [CelY], [CelZ])

    Pomega1.set_data(MT, MOmega1)
    Pomega2.set_data(MT, MOmega2)
    Pomega3.set_data(MT, MOmega3)
    Pomega4.set_data(MT, MOmega4)

    OX.set_data(MT, MOX)
    OY.set_data(MT, MOY)
    OZ.set_data(MT, MOZ)
    OVX.set_data(MT, MOVX)
    OVY.set_data(MT, MOVY)
    OVZ.set_data(MT, MOVZ)
    Pomegax.set_data(MT, MOmegax)
    Pomegay.set_data(MT, MOmegay)
    Pomegaz.set_data(MT, MOmegaz)
    PPhiZ.set_data(MT, MPhiZ)
    VXP, VYP = GetAngles(X, Y, Z, XPoints, YPoints, ZPoints, q0, q1, q2, q3, q0_cam, q1_cam, q2_cam, q3_cam)

    VXP_norm = (VXP + 1) / 2
    VYP_norm = (VYP + 1) / 2
    u_pix = (VXP_norm * image_width).astype(int)
    v_pix = (VYP_norm * image_height).astype(int)
    u_pix = np.clip(u_pix, 0, image_width - 1)
    v_pix = np.clip(v_pix, 0, image_height - 1)

    # -------------------- РџРѕР»СѓС‡РµРЅРёРµ РЅРѕРІС‹С… С†РІРµС‚РѕРІ --------------------
    colors = image_array[v_pix, u_pix] / 255.0

    # -------------------- РћР±РЅРѕРІР»РµРЅРёРµ С†РІРµС‚РѕРІ С‚РѕС‡РµРє --------------------
    PPoints.set_facecolors(colors)
    #PPoints.set_data(VXP, VYP)
    #PPoints.set_offsets(np.c_[VXP, VYP])

    return [C_Copter, K_Copter, A_Copter, B_Copter, Ce_Copter, D_Copter,
            Pomega1, Pomega2, Pomega3, Pomega4, OX, OY, OZ, OVX, OVY, OVZ,
            Pomegax, Pomegay, Pomegaz, PPhiZ, PPoints, PCel]


multik = FuncAnimation(fig, kadr, interval=dt*10, blit=True)

plt.show()