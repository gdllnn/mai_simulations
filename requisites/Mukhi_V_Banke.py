import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

n1 = 20
n2 = 20
n = n1+n2
R = 1
r=0.05
V = 2
dt=0.01
k = 1
k_pr = 0.98
k_m = 1
k_m_pr = 1

Rs = R*np.random.rand(n)
Phis1 = np.pi*np.random.rand(n1)
Phis2 = np.pi*(1+np.random.rand(n2))
Phis = np.concatenate([Phis1,Phis2])
X = Rs*np.cos(Phis)
Y = Rs*np.sin(Phis)
VX = V*(2*np.random.rand(n)-1)
VY = V*(2*np.random.rand(n)-1)

fig = plt.figure(figsize=[15,7])
ax = fig.add_subplot(1,1,1)
ax.axis('equal')
Mukhi1 = ax.plot(X[0:n1],Y[0:n1],'bo')[0]
Mukhi2 = ax.plot(X[n1:n],Y[n1:n],'ro')[0]
Trace1X = np.array(X[0])
Trace1Y = np.array(Y[0])
Trace1 = ax.plot(Trace1X,Trace1Y,'--')[0]
Trace2X = np.array(X[1])
Trace2Y = np.array(Y[1])
Trace2 = ax.plot(Trace2X,Trace2Y,'--')[0]
Trace3X = np.array(X[2])
Trace3Y = np.array(Y[2])
Trace3 = ax.plot(Trace3X,Trace3Y,'--')[0]
PhiB = np.linspace(0,6.29,1001)
Banka = ax.plot(R*np.cos(PhiB),R*np.sin(PhiB),color=[0,0,0])[0]
plt.draw()


def kadr(j):
    global X,Y,VX,VY,dt,Trace1X,Trace1Y,Trace2X,Trace2Y,Trace3X,Trace3Y
    Xn = X + VX * dt
    Yn = Y + VY * dt
    for i in range(n):
        if Xn[i]**2+Yn[i]**2>=R**2:
            Vper_minus = (VX[i]*X[i]+VY[i]*Y[i]) / (X[i]**2+Y[i]**2)**0.5
            Vpar_minus = (VX[i]*Y[i]-VY[i]*X[i]) / (X[i]**2+Y[i]**2)**0.5
            Vpar_plus = Vpar_minus*k_pr
            Vper_plus = -Vper_minus*k
            VX[i] = (Vper_plus*X[i]+Vpar_plus*Y[i]) / (X[i]**2+Y[i]**2)**0.5
            VY[i] = (Vper_plus*Y[i]-Vpar_plus*X[i]) / (X[i]**2+Y[i]**2)**0.5

    for i in range(n-1):
        for j in range(i+1,n):
            if (Xn[i]-Xn[j])**2+(Yn[i]-Yn[j])**2<4*r**2:
                Vper_i_minus = (VX[i] * (Xn[i]-Xn[j]) + VY[i] * (Yn[i]-Yn[j])) / ((Xn[i]-Xn[j]) ** 2 + (Yn[i]-Yn[j]) ** 2) ** 0.5
                Vpar_i_minus = (VX[i] * (Yn[i]-Yn[j]) - VY[i] * (Xn[i]-Xn[j])) / ((Xn[i]-Xn[j]) ** 2 + (Yn[i]-Yn[j]) ** 2) ** 0.5
                Vper_j_minus = (VX[j] * (Xn[i] - Xn[j]) + VY[j] * (Yn[i] - Yn[j])) / ((Xn[i] - Xn[j]) ** 2 + (Yn[i] - Yn[j]) ** 2) ** 0.5
                Vpar_j_minus = (VX[j] * (Yn[i] - Yn[j]) - VY[j] * (Xn[i] - Xn[j])) / ((Xn[i] - Xn[j]) ** 2 + (Yn[i] - Yn[j]) ** 2) ** 0.5
                Vper_i_plus = Vper_j_minus*k_m
                Vpar_i_plus = Vpar_i_minus*k_m_pr
                Vper_j_plus = Vper_i_minus*k_m
                Vpar_j_plus = Vpar_j_minus*k_m_pr
                VX[i] = (Vper_i_plus * (Xn[i]-Xn[j]) + Vpar_i_plus * (Yn[i]-Yn[j])) / ((Xn[i]-Xn[j]) ** 2 + (Yn[i]-Yn[j]) ** 2) ** 0.5
                VY[i] = (Vper_i_plus * (Yn[i]-Yn[j]) - Vpar_i_plus * (Xn[i]-Xn[j])) / ((Xn[i]-Xn[j]) ** 2 + (Yn[i]-Yn[j]) ** 2) ** 0.5
                VX[j] = (Vper_j_plus * (Xn[i] - Xn[j]) + Vpar_j_plus * (Yn[i] - Yn[j])) / ((Xn[i] - Xn[j]) ** 2 + (Yn[i] - Yn[j]) ** 2) ** 0.5
                VY[j] = (Vper_j_plus * (Yn[i] - Yn[j]) - Vpar_j_plus * (Xn[i] - Xn[j])) / ((Xn[i] - Xn[j]) ** 2 + (Yn[i] - Yn[j]) ** 2) ** 0.5
                print('Vzzhh',i,j)
    X+=VX*dt
    Y+=VY*dt
    VX = VX
    VY = VY
    Trace1X = np.append(Trace1X,X[0])
    Trace1Y = np.append(Trace1Y,Y[0])
    Trace2X = np.append(Trace2X, X[1])
    Trace2Y = np.append(Trace2Y, Y[1])
    Trace3X = np.append(Trace3X, X[2])
    Trace3Y = np.append(Trace3Y, Y[2])
    Mukhi1.set_data(X[0:n1], Y[0:n1])
    Mukhi2.set_data(X[n1:n], Y[n1:n])
    Trace1.set_data(Trace1X, Trace1Y)
    Trace2.set_data(Trace2X, Trace2Y)
    Trace3.set_data(Trace3X, Trace3Y)
    return [Mukhi1,Mukhi2,Trace1,Trace2,Trace3]


kino = FuncAnimation(fig, kadr, interval=dt*1000, blit=True)

plt.show()



