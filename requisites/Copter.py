import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import scipy.optimize as sc
#import cv2 as cv



m = 1
L = 0.25
k_A = 0.1
k_B = 0.2
k_C = 0.4
A = m*L**2*k_A
B = m*L**2*k_B
C = m*L**2*k_C
g = 9.81

X = 5
Y = 3
Z = 0

q0 = 1
q1 = 0
q2 = 0
q3 = 0

VX = 2
VY = 0
VZ = 1

Omegax = 1
Omegay = 0
Omegaz = 2

k_t = 1
k_cr = np.array([1.0, 1.0, 1.0, 1.0])
k_cp = np.array([1.0, 1.0, 1.0, 1.0])
k_m = 0.2
R = 0.1
S_v = np.pi*R**2
rho = 1.2
Omega0 = np.sqrt(m*g / (2 * k_t * S_v * rho * R ** 2))
print(Omega0)
DR = 50

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

k_N = 1
k_dN = 0.2

k_R = 0.14
k_dR = 0.2

O_max = 4
O_min = -1

def GetOmegas(X, Y, Z, q0, q1, q2, q3, VX, VY, VZ, Omegax, Omegay, Omegaz):
    global Omega1, Omega2, Omega3, Omega4, PhiZ, \
        CelX, CelY, CelZ, CelVX, CelVY, CelVZ, CelWX, CelWY, CelWZ

    CelVx, CelVy, CelVz = TRot3D(CelVX, CelVY, CelVZ, q0, q1, q2, q3)

    CelFSx = rho * Cx * Sx * CelVx * abs(CelVx) / 2
    CelFSy = rho * Cy * Sy * CelVy * abs(CelVy) / 2
    CelFSz = rho * Cz * Sz * CelVz * abs(CelVz) / 2

    CelFSX, CelFSY, CelFSZ = Rot3D(CelFSx, CelFSy, CelFSz, q0, q1, q2, q3)

    CelFX = - CelFSX - m * CelWX
    CelFY = - CelFSY - m * CelWY
    CelFZ = - m * g - CelFSZ - m * CelWZ

    F_0 = (CelFX**2 + CelFY**2 + CelFZ**2)**0.5
    Omega0 = np.sqrt(F_0 / (2 * k_t * S_v * rho * R ** 2))

    n = - np.array([CelFX, CelFY, CelFZ]) / F_0

    dZ = - k_Z * np.arctan(Z - CelZ) - k_VZ * (VZ - CelVZ)

    X_VC, Y_VC, Z_VC = Rot3D(x_VC, y_VC, z_VC, q0, q1, q2, q3)
    Sm = np.arctan(((X - CelX) ** 2 + (Y - CelY) ** 2) ** 0.5) / ((X - CelX) ** 2 + (Y - CelY) ** 2) ** 0.5

    dXY1 = k_XY * ((X - CelX) * X_VC[0] + (Y - CelY) * Y_VC[0]) * Sm + k_dXY * (
                (VX - CelVX) * X_VC[0] + (VY - CelVY) * Y_VC[0])
    dXY2 = k_XY * ((X - CelX) * X_VC[1] + (Y - CelY) * Y_VC[1]) * Sm + k_dXY * (
                (VX - CelVX) * X_VC[1] + (VY - CelVY) * Y_VC[1])
    dXY3 = k_XY * ((X - CelX) * X_VC[2] + (Y - CelY) * Y_VC[2]) * Sm + k_dXY * (
                (VX - CelVX) * X_VC[2] + (VY - CelVY) * Y_VC[2])
    dXY4 = k_XY * ((X - CelX) * X_VC[3] + (Y - CelY) * Y_VC[3]) * Sm + k_dXY * (
                (VX - CelVX) * X_VC[3] + (VY - CelVY) * Y_VC[3])

    OmX_VC, OmY_VC, OmZ_VC = Rot3D(Omegay * z_VC - Omegaz * y_VC, Omegaz * x_VC - Omegax * z_VC,
                                   Omegax * y_VC - Omegay * x_VC, q0, q1, q2, q3)
    dX_VC, dY_VC, dZ_VC = OmX_VC, OmY_VC, OmZ_VC

    dN1 = - k_N * (X_VC[0]*n[0] + Y_VC[0]*n[1] + Z_VC[0]*n[2]) - k_dN * (dX_VC[0]*n[0] + dY_VC[0]*n[1] + dZ_VC[0]*n[2])
    dN2 = - k_N * (X_VC[1]*n[0] + Y_VC[1]*n[1] + Z_VC[1]*n[2]) - k_dN * (dX_VC[1]*n[0] + dY_VC[1]*n[1] + dZ_VC[1]*n[2])
    dN3 = - k_N * (X_VC[2]*n[0] + Y_VC[2]*n[1] + Z_VC[2]*n[2]) - k_dN * (dX_VC[2]*n[0] + dY_VC[2]*n[1] + dZ_VC[2]*n[2])
    dN4 = - k_N * (X_VC[3]*n[0] + Y_VC[3]*n[1] + Z_VC[3]*n[2]) - k_dN * (dX_VC[3]*n[0] + dY_VC[3]*n[1] + dZ_VC[3]*n[2])

    TC = np.array([CelVX, CelVY, 0]) / (CelVX ** 2 + CelVY ** 2) ** 0.5
    RC = Rot3D(rC[0], rC[1], rC[2], q0, q1, q2, q3)
    Cos = RC[0] * TC[0] + RC[1] * TC[1]
    Sin = RC[1] * TC[0] - RC[0] * TC[1]
    PhiZ = np.arctan2(Sin, Cos)
    dR1 = -k_R * PhiZ - k_dR * Omegaz
    dR2 = k_R * PhiZ + k_dR * Omegaz
    dR3 = -k_R * PhiZ - k_dR * Omegaz
    dR4 = k_R * PhiZ + k_dR * Omegaz
    #print('PhiZ = ', PhiZ, 'dRs = ', dR1, dR2, dR3, dR4)

    ROmega1 = Omega0 * (1 + dZ + dXY1 + dN1 + dR1)
    ROmega2 = Omega0 * (1 + dZ + dXY2 + dN2 + dR2)
    ROmega3 = Omega0 * (1 + dZ + dXY3 + dN3 + dR3)
    ROmega4 = Omega0 * (1 + dZ + dXY4 + dN4 + dR4)
    print('dXY = (', round(dXY1,3),round(dXY2,3),round(dXY3,3),round(dXY4,3), ')')
    print('dN = (', round(dN1, 3), round(dN2, 3), round(dN3, 3), round(dN4, 3), ')')
    print('dZ = ', dZ, ', dXY1 = ', dXY1, 'dN1 = ', dN1, 'dR1 = ', dR1)

    KO1B, KO1M, KO2B, KO2M, KO3B, KO3M, KO4B, KO4M = 1, 1, 1, 1, 1, 1, 1, 1
    if ROmega1 > Omega0 * O_max:
        KO1B = ROmega1 / (Omega0 * O_max)
    elif ROmega1 < Omega0 * O_min:
        KO1M = -ROmega1 / (Omega0 * O_min)
    if ROmega2 > Omega0 * O_max:
        KO2B = ROmega2 / (Omega0 * O_max)
    elif ROmega2 < Omega0 * O_min:
        KO2M = -ROmega2 / (Omega0 * O_min)
    if ROmega3 > Omega0 * O_max:
        KO3B = ROmega3 / (Omega0 * O_max)
    elif ROmega3 < Omega0 * O_min:
        KO3M = -ROmega3 / (Omega0 * O_min)
    if ROmega4 > Omega0 * O_max:
        KO4B = ROmega4 / (Omega0 * O_max)
    elif ROmega4 < Omega0 * O_min:
        KO4M = -ROmega4 / (Omega0 * O_min)

    KO = max(KO1B, KO1M, KO2B, KO2M, KO3B, KO3M, KO4B, KO4M)
    Omega1 = ROmega1 / KO
    Omega2 = ROmega2 / KO
    Omega3 = ROmega3 / KO
    Omega4 = ROmega4 / KO

    #print('Sm = ', Sm, ' Sz = ', np.arctan(Z) / Z, ' KO = ', KO)

    print('Omega1 = ', Omega1, 'Omega2 = ', Omega2, 'Omega3 = ', Omega3, 'Omega4 = ', Omega4)
    return Omega1, Omega2, Omega3, Omega4

def CopterSystemOfEquations(X, Y, Z, q0, q1, q2, q3, VX, VY, VZ, Omegax, Omegay, Omegaz, k_c):
    global Omega1, Omega2, Omega3, Omega4, PhiZ, \
        CelX, CelY, CelZ, CelVX, CelVY, CelVZ

    F1 = k_c[0] * k_t * S_v * rho * R ** 2 * Omega1 * abs(Omega1) / 2
    M1 = k_c[0] * k_m * S_v * rho * R ** 3 * Omega1 * abs(Omega1) / 2
    F2 = k_c[1] * k_t * S_v * rho * R ** 2 * Omega2 * abs(Omega2) / 2
    M2 = k_c[1] * k_m * S_v * rho * R ** 3 * Omega2 * abs(Omega2) / 2
    F3 = k_c[2] * k_t * S_v * rho * R ** 2 * Omega3 * abs(Omega3) / 2
    M3 = k_c[2] * k_m * S_v * rho * R ** 3 * Omega3 * abs(Omega3) / 2
    F4 = k_c[3] * k_t * S_v * rho * R ** 2 * Omega4 * abs(Omega4) / 2
    M4 = k_c[3] * k_m * S_v * rho * R ** 3 * Omega4 * abs(Omega4) / 2

    F_VX, F_VY, F_VZ = Rot3D(0, 0, F1+F2+F3+F4, q0, q1, q2, q3)
    Vx, Vy, Vz = TRot3D(VX, VY, VZ, q0, q1, q2, q3)
    #print('+++++++++++++++++++++++++')
    #print('Vx = ', Vx, ', Vy = ', Vy, ', Vz = ', Vz)
    #print('VX = ', VX, ', VY = ', VY, ', VZ = ', VZ)
    #print('+++++++++++++++++++++++++')
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

def GetAccuracy(k_cp,X, Y, Z, q0, q1, q2, q3, VX, VY, VZ,
                Omegax, Omegay, Omegaz,
                VXK1, VYK1, VZK1, OmegaxK1, OmegayK1, OmegazK1):
    XKo, YKo, ZKo, q0Ko, q1Ko, q2Ko, q3Ko, VXKo, VYKo, VZKo, OmegaxKo, OmegayKo, OmegazKo = \
        CopterSystemOfEquations(
            X, Y, Z, q0, q1, q2, q3,
            VX, VY, VZ, Omegax, Omegay, Omegaz, k_cp)
    return (VXKo - VXK1) ** 2 + (VYKo - VYK1) ** 2 + (VZKo - VZK1) ** 2 + (OmegaxKo - OmegaxK1) ** 2 + \
          (OmegayKo - OmegayK1) ** 2 + (OmegazKo - OmegazK1) ** 2


dt = 0.01

fig = plt.figure(figsize=(19, 9))
ax = fig.add_subplot(1,2,1,projection='3d')
axg = fig.add_subplot(2,4,3)
axo = fig.add_subplot(2,4,4)
axv = fig.add_subplot(2,2,4)
#ax.axis('equal')
ax.set(xlim=[-15*L, 15*L], ylim=[-15*L, 15*L], zlim=[-15*L, 15*L])
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

Nlines = 11
ShLines = 10
ZLines = -30
LPoints = np.linspace(-(Nlines-1)/2*ShLines,(Nlines-1)/2*ShLines,Nlines)
XPoints,YPoints = np.meshgrid(LPoints,LPoints)
XPoints,YPoints = XPoints.ravel(),YPoints.ravel()
ZPoints = XPoints*0 + ZLines
print('XPts = ', XPoints, 'YPts = ', YPoints, 'ZPts = ', ZPoints)
def GetAngles(X, Y, Z, XPoints,YPoints,ZPoints, q0, q1, q2, q3,
              q0_cam, q1_cam, q2_cam, q3_cam):
    Norms = ((XPoints - X) ** 2 + (YPoints - Y) ** 2 + (ZPoints - Z) ** 2) ** 0.5
    VXP, VYP, VZP = (XPoints - X)/Norms, (YPoints - Y)/Norms, (ZPoints - Z)/Norms
    TVXP, TVYP, TVZP = TRot3D(VXP, VYP, VZP, q0, q1, q2, q3)
    TVXP, TVYP, TVZP = TRot3D(TVXP, TVYP, TVZP, q0_cam, q1_cam, q2_cam, q3_cam)
    #Psi = np.arcsin(VXP)
    #Phi = np.arcsin(VYP/np.cos(Psi))
    return TVXP,TVYP#Psi,Phi

Alpha_cam = 1.2
q0_cam = np.cos(Alpha_cam/2)
q1_cam = np.sin(Alpha_cam/2)
q2_cam = 0
q3_cam = 0
VXP,VYP = GetAngles(X, Y, Z, XPoints,YPoints,ZPoints,
                    q0, q1, q2, q3, q0_cam, q1_cam, q2_cam, q3_cam)
PPoints = axv.plot(VXP,VYP,'o')[0]
axv.plot([-0.65, 0.65, 0.65, -0.65, -0.65],[-0.65, -0.65, 0.65, 0.65, -0.65],'r')

CelX = 0
CelY = 5
CelZ = 0
CelVX = 5*0.4
CelVY = 0
CelVZ = 2.5*0.8

CelX = 1
CelY = 0
CelZ = 1
CelVX = 0
CelVY = 0.1
CelVZ = 0

X, Y, Z = CelX, CelY+1, CelZ+1
VX, VY, VZ = CelVX, CelVY+1, CelVZ+1

PCel = ax.plot(CelX, CelY, CelZ, 'o', color=[0.7, 0, 0.5])[0]
Delta0 = np.sqrt((X - CelX)**2 + (Y - CelY)**2 +(Z - CelZ)**2)
MDelta = np.array([])
Delta = axg.plot(MT,MOX)[0]

def kadr(i):
    global X, Y, Z, q0, q1, q2, q3, VX, VY, VZ, Omegax, Omegay, Omegaz,\
        Omega1, Omega2, Omega3, Omega4, PhiZ, \
        MOmega1, MOmega2, MOmega3, MOmega4, MDelta, MOX, MOY, MOZ, \
        MOVX, MOVY, MOVZ, MOmegax, MOmegay, MOmegaz, MPhiZ, MT, \
        CelX, CelY, CelZ, CelVX, CelVY, CelVZ, CelWX, CelWY, CelWZ, \
        q0_cam, q1_cam, q2_cam, q3_cam, k_cp

    CelPPX = 5 * np.sin(0.4 * (i - 2) * dt)
    CelPPY = 5 * np.cos(0.4 * (i - 2) * dt)
    CelPPZ = 2.5 * np.sin(1.6 * (i - 2) * dt)

    CelPX = 5 * np.sin(0.4 * (i - 1) * dt)
    CelPY = 5 * np.cos(0.4 * (i - 1) * dt)
    CelPZ = 2.5 * np.sin(1.6 * (i - 1) * dt)

    CelX = 5 * np.sin(0.4 * i*dt)
    CelY = 5 * np.cos(0.4 * i*dt)
    CelZ = 2.5 * np.sin(1.6 * i*dt)

    CelVX = (CelX - CelPX) / dt
    CelVY = (CelY - CelPY) / dt
    CelVZ = (CelZ - CelPZ) / dt

    CelWX = (CelX - 2*CelPX + CelPPX) / dt**2
    CelWY = (CelY - 2*CelPY + CelPPY) / dt**2
    CelWZ = (CelZ - 2*CelPZ + CelPPZ) / dt**2

    # CelX = 1
    # CelY = 0
    # CelZ = 1
    # CelVX = 0
    # CelVY = 0.1
    # CelVZ = 0
    if i>0:
        k_cr = np.array([1.0, 1.0, 1.0, 1.0])
    else:
        k_cr = np.array([0.6, 0.5, 1.0, 1.0])
    Omega1, Omega2, Omega3, Omega4 = GetOmegas(X, Y, Z, q0, q1, q2, q3, VX, VY, VZ, Omegax, Omegay, Omegaz)

    XK1, YK1, ZK1, q0K1, q1K1, q2K1, q3K1, VXK1, VYK1, VZK1, OmegaxK1, OmegayK1, OmegazK1 = \
        CopterSystemOfEquations(
            X, Y, Z, q0, q1, q2, q3,
            VX, VY, VZ, Omegax, Omegay, Omegaz, k_cr)
    XKo, YKo, ZKo, q0Ko, q1Ko, q2Ko, q3Ko, VXKo, VYKo, VZKo, OmegaxKo, OmegayKo, OmegazKo = \
        CopterSystemOfEquations(
            X, Y, Z, q0, q1, q2, q3,
            VX, VY, VZ, Omegax, Omegay, Omegaz, k_cp)
    Otl = (VXKo - VXK1)**2 + (VYKo - VYK1)**2 + (VZKo - VZK1)**2 + (OmegaxKo - OmegaxK1)**2 + \
          (OmegayKo - OmegayK1)**2 + (OmegazKo - OmegazK1)**2
    if Otl > DR:
        print("Khhryas")
        Opt = sc.minimize(GetAccuracy,k_cp,(X, Y, Z, q0, q1, q2, q3,
            VX, VY, VZ, Omegax, Omegay, Omegaz,VXK1, VYK1, VZK1, OmegaxK1, OmegayK1, OmegazK1))
        k_cp = Opt.x
    else:
        print('Norm')
    print('k_cp = ',k_cp)
    print('Otl = ',Otl)
    #Omega1K1, Omega2K1, Omega3K1, Omega4K1 = Omega1, Omega2, Omega3, Omega4

    XK2, YK2, ZK2, q0K2, q1K2, q2K2, q3K2, VXK2, VYK2, VZK2, OmegaxK2, OmegayK2, OmegazK2 = \
        CopterSystemOfEquations(
            X+dt/2*XK1, Y+dt/2*YK1, Z+dt/2*ZK1,
            q0+dt/2*q0K1, q1+dt/2*q1K1, q2+dt/2*q2K1, q3+dt/2*q3K1,
            VX+dt/2*VXK1, VY+dt/2*VYK1, VZ+dt/2*VZK1,
            Omegax+dt/2*OmegaxK1, Omegay+dt/2*OmegayK1, Omegaz+dt/2*OmegazK1, k_cr)
    #Omega1K2, Omega2K2, Omega3K2, Omega4K2 = Omega1, Omega2, Omega3, Omega4

    XK3, YK3, ZK3, q0K3, q1K3, q2K3, q3K3, VXK3, VYK3, VZK3, OmegaxK3, OmegayK3, OmegazK3 = \
        CopterSystemOfEquations(
            X + dt/2 * XK2, Y + dt / 2 * YK2, Z + dt / 2 * ZK2,
            q0 + dt/2*q0K2, q1 + dt / 2 * q1K2, q2 + dt / 2 * q2K2, q3 + dt / 2 * q3K2,
            VX + dt/2*VXK2, VY + dt / 2 * VYK2, VZ + dt / 2 * VZK2,
            Omegax + dt/2*OmegaxK2, Omegay + dt/2*OmegayK2, Omegaz + dt / 2 * OmegazK2, k_cr)
    #Omega1K3, Omega2K3, Omega3K3, Omega4K3 = Omega1, Omega2, Omega3, Omega4

    XK4, YK4, ZK4, q0K4, q1K4, q2K4, q3K4, VXK4, VYK4, VZK4, OmegaxK4, OmegayK4, OmegazK4 = \
        CopterSystemOfEquations(
            X + dt * XK3, Y + dt * YK3, Z + dt * ZK3,
            q0 + dt * q0K3, q1 + dt * q1K3, q2 + dt * q2K3, q3 + dt * q3K3,
            VX + dt * VXK3, VY + dt * VYK3, VZ + dt * VZK3,
            Omegax + dt * OmegaxK3, Omegay + dt * OmegayK3, Omegaz + dt * OmegazK3, k_cr)
    #Omega1K4, Omega2K4, Omega3K4, Omega4K4 = Omega1, Omega2, Omega3, Omega4

    X += (XK1 + 2 * XK2 + 2 * XK3 + XK4) * dt / 6
    Y += (YK1 + 2 * YK2 + 2 * YK3 + YK4) * dt / 6
    Z += (ZK1 + 2 * ZK2 + 2 * ZK3 + ZK4) * dt / 6

    q0 += (q0K1 + 2 * q0K2 + 2 * q0K3 + q0K4) * dt / 6
    q1 += (q1K1 + 2 * q1K2 + 2 * q1K3 + q1K4) * dt / 6
    q2 += (q2K1 + 2 * q2K2 + 2 * q2K3 + q2K4) * dt / 6
    q3 += (q3K1 + 2 * q3K2 + 2 * q3K3 + q3K4) * dt / 6
    m_q = np.sqrt(q0**2+q1**2+q2**2+q3**2)
    q0 = q0 / m_q
    q1 = q1 / m_q
    q2 = q2 / m_q
    q3 = q3 / m_q

    VX += (VXK1 + 2 * VXK2 + 2 * VXK3 + VXK4) * dt / 6
    VY += (VYK1 + 2 * VYK2 + 2 * VYK3 + VYK4) * dt / 6
    VZ += (VZK1 + 2 * VZK2 + 2 * VZK3 + VZK4) * dt / 6

    Omegax += (OmegaxK1 + 2 * OmegaxK2 + 2 * OmegaxK3 + OmegaxK4) * dt / 6
    Omegay += (OmegayK1 + 2 * OmegayK2 + 2 * OmegayK3 + OmegayK4) * dt / 6
    Omegaz += (OmegazK1 + 2 * OmegazK2 + 2 * OmegazK3 + OmegazK4) * dt / 6

    #Omega1 = (Omega1K1 + 2 * Omega1K2 + 2 * Omega1K3 + Omega1K4) / 6
    #Omega2 = (Omega2K1 + 2 * Omega2K2 + 2 * Omega2K3 + Omega2K4) / 6
    #Omega3 = (Omega3K1 + 2 * Omega3K2 + 2 * Omega3K3 + Omega3K4) / 6
    #Omega4 = (Omega4K1 + 2 * Omega4K2 + 2 * Omega4K3 + Omega4K4) / 6

    MOmega1 = np.append(MOmega1,Omega1)
    MOmega2 = np.append(MOmega2,Omega2)
    MOmega3 = np.append(MOmega3,Omega3)
    MOmega4 = np.append(MOmega4,Omega4)

    MOX = np.append(MOX,X-CelX)
    MOY = np.append(MOY,Y-CelY)
    MOZ = np.append(MOZ,Z-CelZ)
    MOVX = np.append(MOVX,VX-CelVX)
    MOVY = np.append(MOVY,VY-CelVY)
    MOVZ = np.append(MOVZ,VZ-CelVZ)
    MOmegax = np.append(MOmegax,Omegax)
    MOmegay = np.append(MOmegay,Omegay)
    MOmegaz = np.append(MOmegaz,Omegaz)
    MPhiZ = np.append(MPhiZ,PhiZ)

    MT = np.append(MT,i*dt)

    MDelta = np.append(MDelta,500*np.sqrt((X - CelX) ** 2 + (Y - CelY) ** 2 + (Z - CelZ) ** 2)/Delta0)


    print('-----------------  ', i,'  -----------------')
    print('X = ', X,'Y = ', Y,'Z = ',Z)
    print('VX = ', VX, 'VY = ', VY, 'VZ = ', VZ)
    print('Omegax = ', Omegax, ', Omegay = ', Omegay, ', Omegaz = ', Omegaz)
    print('q = (', q0,q1,q2,q3, '), |q| = ', np.sqrt(q0**2+q1**2+q2**2+q3**2), '.')

    C_Copter.set_data_3d(X, Y, Z)
    RX_Copter, RY_Copter, RZ_Copter = Rot3D(X_Copter, Y_Copter, Z_Copter, q0,q1,q2,q3)
    K_Copter.set_data_3d(X + RX_Copter, Y + RY_Copter, Z + RZ_Copter)
    A_Copter.set_data_3d(X + RX_Copter[1], Y + RY_Copter[1], Z + RZ_Copter[1])
    B_Copter.set_data_3d(X + RX_Copter[3], Y + RY_Copter[3], Z + RZ_Copter[3])
    Ce_Copter.set_data_3d(X + RX_Copter[5], Y + RY_Copter[5], Z + RZ_Copter[5])
    D_Copter.set_data_3d(X + RX_Copter[7], Y + RY_Copter[7], Z + RZ_Copter[7])
    PCel.set_data_3d(CelX,CelY,CelZ)

    Pomega1.set_data(MT, MOmega1)
    Pomega2.set_data(MT, MOmega2)
    Pomega3.set_data(MT, MOmega3)
    Pomega4.set_data(MT, MOmega4)
    Delta.set_data(MT, MDelta)

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
    VXP, VYP = GetAngles(X, Y, Z, XPoints, YPoints, ZPoints,
                         q0, q1, q2, q3, q0_cam, q1_cam, q2_cam, q3_cam)
    PPoints.set_data(VXP, VYP)

    return [C_Copter,K_Copter,A_Copter,B_Copter,Ce_Copter,D_Copter,
            Pomega1,Pomega2,Pomega3,Pomega4, OX,OY,OZ,OVX,OVY,OVZ,
            Pomegax,Pomegay,Pomegaz,PPhiZ,PPoints, PCel,Delta]


multik = FuncAnimation(fig, kadr, interval=dt*10, blit=True)

plt.show()





def PositionProba(X, Y, Z, q0, q1, q2, q3, XPoints,YPoints,ZPoints, CTVXP,CTVYP):
    Norms = ((XPoints - X) ** 2 + (YPoints - Y) ** 2 + (ZPoints - Z) ** 2) ** 0.5
    VXP, VYP, VZP = (XPoints - X)/Norms, (YPoints - Y)/Norms, (ZPoints - Z)/Norms
    TVXP, TVYP, TVZP = TRot3D(VXP, VYP, VZP, q0, q1, q2, q3)
    Different = sum((CTVXP-TVXP)**2) + sum((CTVYP-TVYP)**2)
    return Different

#PositionProba->min
#X, Y, Z, q0, q1, q2, q3