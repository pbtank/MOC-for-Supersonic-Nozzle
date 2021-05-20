#############################################################################
#                          Original Author                                  #
#                          Shubham Maurya                                   #
#         Indian Institute of Space Science and Technology                  #
#                           2017 | MATLAB                                   #
#                                                                           #
#                            Modified by                                    #
#                           Priyansu Tank                                   #
#                     priyansutank3044@gmail.com                            #
#                Vishwakarma Govt. Engg. College, GTU                       #
#                           2021 | python                                   #
#                                                                           #
# This program generates contour for a minimum length nozzle. The flow      #
# expands through sharp corner at throat.                                   #
# The input can be of two types:                                            #
# [1] You have chamber properties. Then specify T_c, P_c, g, FT, m_dot, D,  #
#     ALT and n.                                                            #
# [2] If you already have Me, then only Me, g, D and n are required.        #
# Also you can export nozzle contour points as x,y,z format in a .csv file. #
#############################################################################

import numpy as np
import scipy as sp
import math
from matplotlib import pyplot as plt

############# PROBLEM PAERAMETERS ################
T_c = 3673          # CHAMBER TEMPERATURE    [K]
P_c = 9.72e6        # CHAMBER PRESSURE       [Pa]
g = 1.4             # SP. HEAT RATIO(GAMMA)  
FT = 654333.397     # THRUST                 [N]
m_dot = 236.6       # MASS FLOW RATE         [kg/s]
D = 0.1156          # THROAT RADIUS          [m]
ALT = 0             # ALTITUDE               [m]

n = 10              #Incident expansion wave conditions

# ATMOSPHERIC MODEL BY NASA GRC
if ((11000<ALT) and (ALT<25000)):
    T_amb = -56.46
    p_amb = 1000*(22.65*np.exp(1.73-0.000157*ALT))
elif ALT>= 25000:
    T_amb = -131.21 + 0.00299*ALT
    p_amb = 1000*(2.488*((T_amb+273.1)/216.6)**-11.388)
else:
    T_amb = 15.04 - 0.00649*ALT
    p_amb = 1000*(101.29*((T_amb+273.1)/288.08)**5.256)

# EXIT MACH NUMBER 
# (IF YOU HAVE Me THEN ONLY REQUIRED PARAMETERS 
# ARE g, Me, n, D)
# Me = math.sqrt(((P_c/p_amb)**((g-1)/g)-1)*(2/(g-1)))
Me = 3

RTOD = 180/math.pi
DTOR = math.pi/180

def PrandtlMeyer(M, g):
    A = math.sqrt((g+1)/(g-1))
    B = (g-1)/(g+1)
    v = A*math.atan(math.sqrt(B*((M**2)-1)))-math.atan(math.sqrt((M**2)-1))
    return(v)

## characteristic parameter solver
# v - Prandtl-Meyer function
# KL - (k+)Left running characteristic constant
# KR - (k-)Right running characteristic constant
# theta - Flow angle relative to horizontal
def moc2d(theta_max, theta_0, n):
    dtheta = (theta_max - theta_0)/(n-1)

    node = int(0.5*n*(4+n-1))
    theta = np.zeros(node)
    v = np.zeros(node)
    KL = np.zeros(node);
    KR = np.zeros(node);
    for i in range(0,n):
        theta[i] = theta_0 + (i)*dtheta
        v[i] = theta[i]
        KL[i] = theta[i]-v[i]  #k+
        KR[i] = theta[i]+v[i]  #k-
    i = n
    theta[i] = theta[i-1]
    v[i] = v[i-1]
    KL[i] = KL[i-1]
    KR[i] = KR[i-1]

    st = np.zeros(n)
    en = np.zeros(n)
    for i in range(0,n):
        if i==0:
            st[i] = 0
        else:
            st[i] = en[i-1]+1
        en[i] = st[i]+n-i

    l = 1
    for i in range(1,n):
        k = n
        for j in range(int(en[i]), int(st[i]-1), -1):
            if j==st[i]:
                theta[j] = 0
            KR[j] = KR[k]
            KL[j] = -KR[l]
            theta[j] = 0.5*(KR[j]+KL[j])
            v[j] = 0.5*(KR[j]-KL[j])
            k=k-1
        l=l+1

    return (np.round(v,2), np.round(KL,2), np.round(KR,2), np.round(theta,2))

# Based on Hall, I. M. "Inversion of the prandtl-meyer relation."
# Aeronautical Journal 79 (1975): 417.
def InvPrandtlMeyer(v):
    v = v*math.pi/180
    A=1.3604
    B=0.0962
    C=-0.5127
    D=-0.6722
    E=-0.3278
    v_0 = 0.5*math.pi*(math.sqrt(6)-1)
    y = (v/v_0)**(2/3)
    Mach = (1 + A*y + B*y**2 + C*y**3)/(1 + D*y + E*y**2)
    return(Mach)

def Mu(MM):
    mu = math.asin(1/MM)
    return(mu*RTOD)
    
############################################################################

exp = input("Export export nozzle contour point? [y/n] : ")

if exp=="y":
    f = open("points.csv", "w")

theta_max = PrandtlMeyer(Me, g)*RTOD*0.5

theta_0 = round(theta_max - math.floor(theta_max), 2)

v, KL, KR, theta = moc2d(theta_max, theta_0, int(n))

# Mach number and Mach angle at each node
node = int(0.5*n*(4+n-1))
M = np.zeros(node)
mu = np.zeros(node)
for i in range(0,node):
    M[i] = InvPrandtlMeyer(v[i])
    mu[i] = Mu(M[i])
np.round(M, 4)
np.round(mu, 4)

# Grid generator
i = 0
x = np.zeros(node)
y = np.zeros(node)
wall = theta_max
while (i<=n):
    if i==0:
        x[i] = round(-D/(math.tan((theta[i]-mu[i])*DTOR)), 5)
        y[i] = 0
        plt.plot([0, x[i]], [D, 0], 'r')
    elif i==n:
        x[i] = round((y[i-1]-D-x[i-1]*math.tan(DTOR*(theta[i-1]+theta[i]+mu[i-1]+mu[i])*0.5))/(math.tan(DTOR*0.5*(wall+theta[i]))-math.tan(DTOR*(theta[i-1]+theta[i]+mu[i-1]+mu[i])*0.5)), 5)
        y[i] = round(D + (x[i]*math.tan(DTOR*0.5*(wall+theta[i]))), 5)
        plt.plot([x[i-1], x[i]], [y[i-1], y[i]], 'b')
        plt.plot([0, x[i]], [D, y[i]], 'g', marker='o')
        if exp=="y":
            f.write(str(0) + "," + str(D) + "," + str(0) + "\n")
            f.write(str(x[i]) + "," + str(y[i]) + "," + str(0) + "\n")
    else:
        x[i] = round((D-y[i-1]+x[i-1]*math.tan(DTOR*(theta[i-1]+theta[i]+mu[i-1]+mu[i])*0.5))/(math.tan(DTOR*0.5*(mu[i-1]+theta[i-1]+mu[i]+theta[i]))-math.tan(DTOR*(theta[i]-mu[i]))), 5)
        y[i] = round(D + x[i]*math.tan(DTOR*(theta[i]-mu[i])), 5)
        plt.plot([x[i-1], x[i]], [y[i-1], y[i]], 'b')
        plt.plot([0, x[i]], [D, y[i]], 'r')
    i=i+1

h=i
k=0
i=h
for j in range(0,n-1):
    while (i<=h+n-k-1):
        if i==h:
            x[i] = round(x[i-n+k]-y[i-n+k]/(math.tan(DTOR*0.5*(theta[i-n+k]+theta[i]-mu[i-n+k]-mu[i]))), 5)
            y[i] = 0

            plt.plot([x[i-n+k], x[i]], [y[i-n+k], y[i]], 'r')
        elif i==h+n-k-1:
            x[i] = round((x[i-n+k]*math.tan(DTOR*0.5*(theta[i-n+k]+theta[i]))-y[i-n+k]+y[i-1]-x[i-1]*math.tan(DTOR*(theta[i-1]+theta[i]+mu[i-1]+mu[i])*0.5))/(math.tan(DTOR*0.5*(theta[i-n+k]+theta[i]))-math.tan(DTOR*(theta[i-1]+theta[i]+mu[i-1]+mu[i])*0.5)), 5)
            y[i] = round(y[i-n+k]+(x[i]-x[i-n+k])*math.tan(DTOR*0.5*(theta[i-n+k]+theta[i])), 5)
            plt.plot([x[i-1], x[i]], [y[i-1], y[i]], 'b')
            plt.plot([x[i-n+k], x[i]], [y[i-n+k], y[i]], 'g', marker='o')
            if exp=="y":
                f.write(str(x[i]) + "," + str(y[i]) + "," + str(0) + "\n")
        else:
            s1 = math.tan(DTOR*0.5*(theta[i]+theta[i-1]+mu[i]+mu[i-1]))
            s2 = math.tan(DTOR*0.5*(theta[i]+theta[i-n+k]-mu[i]-mu[i-n+k]))
            x[i] = round((y[i-n+k]-y[i-1]+s1*x[i-1]-s2*x[i-n+k])/(s1-s2), 5)
            y[i] = round(y[i-1]+(x[i]-x[i-1])*s1, 5)
            plt.plot([x[i-1], x[i]], [y[i-1], y[i]], 'b')
            plt.plot([x[i-n+k], x[i]], [y[i-n+k], y[i]], 'r')
        i=i+1
    k=k+1
    h=i
    i=h

if exp=="y":
    f.close()
    print("Exported successfully in file points.csv")
plt.axis([0,x[i-1]+0.5, 0, y[i-1]+0.5])
plt.show()



























