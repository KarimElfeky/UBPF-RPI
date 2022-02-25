# Importing COM Interface and Creating the OpenDSS engine object.
from cmath import cos, pi, sqrt,sin
from timeit import default_timer, timeit
from GandBbus_V3 import *
from Loadextract import *
import numpy as nmp
from guppy import hpy
from math import *
import subprocess
Ybus = nmp.loadtxt("Ybus.csv",delimiter=",")
Ybus = tuple(Ybus)
BusNo = 76

Busnamesfile = open("Busnames.txt","r")
file_lines = Busnamesfile.read()
Busnames = file_lines.split("\n")
G,B= GandBbus(BusNo,Ybus)


#Input to the Algorithm

P = nmp.loadtxt("P.csv",delimiter=",")
Q = nmp.loadtxt("Q.csv",delimiter=",")
G1 = nmp.array(G)
B1 = nmp.array(B)

Vbase = 230.94
Sbase = 10e3
Zbase = Vbase**2/Sbase
#G = G*Zbase
#B = B*Zbase
Th_sl = nmp.array([0,-120,120,0])
V_sl = nmp.array([1,1,1])
Th_sl = Th_sl*pi/180
PC = -1*P*1000/Sbase
QC = -1*Q*1000/Sbase

#Theta = nmp.zeros((BusNo,4))
Theta = [[0]*4 for i in range(BusNo)]
for i in range(BusNo):
    Theta[i][:] = Th_sl
#V = nmp.ones((BusNo,4))
#V[:,3] = nmp.ones((BusNo,))/Vbase

V = [[1,1,1,1/Vbase] for i in range(BusNo)]

sl_busNo = nmp.array([1]) #input the slack bus number

PQbuses = list(Busnames[:])
PQbusind = list(range(BusNo))
for i in sl_busNo:
    sl_ind = PQbuses.index(str(i))
    PQbuses.remove(str(i))
    V[sl_ind][3] = 0
    PQbusind.remove(sl_ind)

PQbusesNos = nmp.array([int(i) for i in PQbuses])
NoPQBuses = len(PQbusesNos)

iter_reduce = iter_reduction(G1,B1)
#iter_reduce1 = nmp.array(iter_reduce)
#iter_reduce = nmp.array(iter_reduce)
#iter_reduce = iter_reduce.astype(int)
#nmp.savetxt("iter_reduce_py.csv",iter_reduce,delimiter=",")
iter =0
Threshold = 0.0001
#Jtime=list()
startt = timeit,default_timer()
#h = hpy()
while True:
    xV = nmp.array(V[1:])
    xTheta = nmp.array(Theta[1:])
    xV = nmp.reshape(xV,(NoPQBuses*4,),'C')
    xTheta = nmp.reshape(xTheta,(NoPQBuses*4,),'C')
    #xV = nmp.reshape(V[1:,:],(NoPQBuses*4,),'C')
    #xTheta = nmp.reshape(Theta[1:,:],(NoPQBuses*4,),'C')
    x = nmp.concatenate([xTheta,xV])
    iter = iter +1
    #mismatch functions
    P = nmp.zeros((NoPQBuses,3))
    Q = nmp.zeros((NoPQBuses,3))
    fCbr = [0]*NoPQBuses
    fCbi = [0]*NoPQBuses
    for n in iter_reduce:
        if (n[0]!=0):
            f_PQind = n[0]
            f_phase = n[2]
            f_PQ = f_PQind-1  
            t_bus = n[1]
            t_phase = n[3]
            ind1,ind2 = Index_calc(f_PQind,t_bus,f_phase,t_phase)
            the_V1 = V[f_PQind][f_phase]
            the_V2 = V[t_bus][t_phase]
            the_G = G[ind1][ind2]
            the_B = B[ind1][ind2]
            the_theta1 = Theta[f_PQind][f_phase]
            the_theta2 = Theta[t_bus][t_phase]
            d_theta = the_theta1 - the_theta2
            fCbr[f_PQ] = fCbr[f_PQ]+the_G*(the_V1*cos(the_theta1)-the_V2*cos(the_theta2))+the_B*(-the_V1*sin(the_theta1)+the_V2*sin(the_theta2))
            fCbi[f_PQ] = fCbi[f_PQ]+the_G*(the_V1*sin(the_theta1)-the_V2*sin(the_theta2))+the_B*(+the_V1*cos(the_theta1)-the_V2*cos(the_theta2))        
            if (n[2]!=3):
                P[f_PQ,f_phase] = P[f_PQ,f_phase]+the_V1*the_V2*((the_G*cos(d_theta))+(the_B*sin((d_theta))))
                Q[f_PQ,f_phase] = Q[f_PQ,f_phase]+the_V1*the_V2*((the_G*sin(d_theta))-(the_B*cos((d_theta))))       
    Pn = nmp.zeros((NoPQBuses,3))
    Qn = nmp.zeros((NoPQBuses,3))
    for f_PQ in range(NoPQBuses):
        f_PQind = PQbusind[f_PQ]
        for f_phase in range(3):
            the_Vn = V[f_PQind][3]
            the_V1 = V[f_PQind][f_phase]
            the_P = P[f_PQ,f_phase]
            the_Q = Q[f_PQ,f_phase]
            d_theta2 = Theta[f_PQind][3]-Theta[f_PQind][f_phase]
            Pn[f_PQ,f_phase]=the_Vn/the_V1*(-the_P*cos(d_theta2)+the_Q*sin(d_theta2))
            Qn[f_PQ,f_phase]=the_Vn/the_V1*(-the_P*sin(d_theta2)-the_Q*cos(d_theta2))
    #P = nmp.array(P)
    #Q = nmp.array(Q)
    #Pn = nmp.array(Pn)
    #Qn = nmp.array(Qn)
    fP = P+Pn-PC[1:,]
    fQ = Q+Qn-QC[1:,]
    #P = P.tolist()
    #Q = Q.tolist()
    fCb = nmp.concatenate([fCbr,fCbi])
    fP = nmp.reshape(fP,(NoPQBuses*3,),'F')
    fQ = nmp.reshape(fQ,(NoPQBuses*3,),'F')
    fmismatch = nmp.concatenate([fP,fQ,fCb])
    if max(nmp.absolute(fmismatch))<Threshold:
        break
    # Jacobian Matrix Calculation
    # The change of active and reactive power with respect to the phase angle (Theta)
    dP_dTh = nmp.zeros((NoPQBuses*3,NoPQBuses*4))
    dQ_dTh = nmp.zeros((NoPQBuses*3,NoPQBuses*4))
    dP_dV = nmp.zeros((NoPQBuses*3,NoPQBuses*4))
    dQ_dV = nmp.zeros((NoPQBuses*3,NoPQBuses*4))
    dPn_dTh = nmp.zeros((NoPQBuses*3,NoPQBuses*4))
    dQn_dTh = nmp.zeros((NoPQBuses*3,NoPQBuses*4))
    dPn_dV = nmp.zeros((NoPQBuses*3,NoPQBuses*4))
    dQn_dV = nmp.zeros((NoPQBuses*3,NoPQBuses*4))
    dfCbr_dTh = nmp.zeros((NoPQBuses,NoPQBuses*4))
    dfCbi_dTh = nmp.zeros((NoPQBuses,NoPQBuses*4))
    dfCbr_dV = nmp.zeros((NoPQBuses,NoPQBuses*4))
    dfCbi_dV = nmp.zeros((NoPQBuses,NoPQBuses*4))
    for n in iter_reduce:
        t_bus = n[1]
        t_phase = n[3]
        if t_bus !=0 and n[0]!=0:
            t_PQind = t_bus
            t_bus = t_bus - 1
            f_PQind = n[0]
            f_phase = n[2]
            f_PQ = f_PQind-1 
            ind1,ind2 = Index_calc(f_PQind,t_PQind,f_phase,t_phase)
            ind1J,ind2J = Index_calcJ(f_PQ,t_bus,f_phase,t_phase)
            ind1Jo,ind2Jo = Index_calcJ(f_PQ,t_bus,f_phase,t_phase)
            indJ = Index_calcJ2(t_bus,t_phase)
            the_V1 = V[f_PQind][f_phase]
            the_V2 = V[t_PQind][t_phase]
            the_Vn = V[f_PQind][3]
            the_G = G[ind1][ind2]
            the_B = B[ind1][ind2]
            the_theta1 = Theta[f_PQind][f_phase]
            the_theta2 = Theta[t_PQind][t_phase]
            d_theta = the_theta1 - the_theta2
            d_theta2 = Theta[f_PQind][3] - the_theta1
            dfCbr_dTh[f_PQ,indJ] = dfCbr_dTh[f_PQ,indJ] +  the_G*the_V2*sin(the_theta2)+the_B*the_V2*cos(the_theta2)
            dfCbi_dTh[f_PQ,indJ] = dfCbi_dTh[f_PQ,indJ] -  the_G*the_V2*cos(the_theta2)+the_B*the_V2*sin(the_theta2)
            dfCbr_dV[f_PQ,indJ] = dfCbr_dV[f_PQ,indJ] -  the_G*cos(the_theta2)+the_B*sin(the_theta2)
            dfCbi_dV[f_PQ,indJ] = dfCbi_dV[f_PQ,indJ] -  the_G*sin(the_theta2)-the_B*cos(the_theta2)
            if n[2]!=3:
                if f_PQ != t_bus or f_phase != t_phase:
                    dP_dTh[ind1J,ind2J] = the_V1*the_V2*(+the_G*sin(d_theta)-the_B*cos(d_theta))
                    dQ_dTh[ind1J,ind2J] = the_V1*the_V2*(-the_G*cos(d_theta)-the_B*sin(d_theta))
                    dP_dV[ind1J,ind2J] = the_V1*(the_G*cos(d_theta)+the_B*sin(d_theta))
                    dQ_dV[ind1J,ind2J] = the_V1*(the_G*sin(d_theta)-the_B*cos(d_theta))
                ind11J,ind22J = Index_calcJ(f_PQ,f_PQ,f_phase,f_phase)
                #or dQ_dTh[ind11J,ind22J]==0 or dP_dV[ind11J,ind22J]==0 or dQ_dV[ind11J,ind22J]==0
                if dP_dTh[ind11J,ind22J]==0:
                    ind11,ind22 = Index_calc(f_PQind,f_PQind,f_phase,f_phase)
                    the_G2 = G[ind11][ind22]
                    the_B2 = B[ind11][ind22]
                    dP_dTh[ind11J,ind22J]=-Q[f_PQ,f_phase]-the_B2*the_V1**2
                    dQ_dTh[ind11J,ind22J]= P[f_PQ,f_phase]-the_G2*the_V1**2
                    dP_dV[ind11J,ind22J]=  (P[f_PQ,f_phase]/the_V1) + the_V1*the_G2
                    dQ_dV[ind11J,ind22J]=  (Q[f_PQ,f_phase]/the_V1) - the_V1*the_B2
                if t_PQind == f_PQind and t_phase==3:
                    dPn_dTh[ind1J,ind2J] = the_Vn/the_V1*(-dP_dTh[ind1Jo,ind2Jo]*cos(d_theta2)+P[f_PQ,f_phase]*sin(d_theta2)+dQ_dTh[ind1Jo,ind2Jo]*sin(d_theta2)+Q[f_PQ,f_phase]*cos(d_theta2))
                    dQn_dTh[ind1J,ind2J] = the_Vn/the_V1*(-dP_dTh[ind1Jo,ind2Jo]*sin(d_theta2)-P[f_PQ,f_phase]*cos(d_theta2)-dQ_dTh[ind1Jo,ind2Jo]*cos(d_theta2)+Q[f_PQ,f_phase]*sin(d_theta2))
                    dPn_dV[ind1J,ind2J] = 1/the_V1*(-P[f_PQ,f_phase]*cos(d_theta2)+Q[f_PQ,f_phase]*sin(d_theta2))+the_Vn/the_V1*(-dP_dV[ind1Jo,ind2Jo]*cos(d_theta2)+dQ_dV[ind1Jo,ind2Jo]*sin(d_theta2))
                    dQn_dV[ind1J,ind2J] = 1/the_V1*(-P[f_PQ,f_phase]*sin(d_theta2)-Q[f_PQ,f_phase]*cos(d_theta2))+the_Vn/the_V1*(-dP_dV[ind1Jo,ind2Jo]*sin(d_theta2)-dQ_dV[ind1Jo,ind2Jo]*cos(d_theta2))
                elif t_PQind == f_PQind and t_phase==f_phase:
                    dPn_dTh[ind1J,ind2J] = the_Vn/the_V1*(-dP_dTh[ind1Jo,ind2Jo]*cos(d_theta2)-P[f_PQ,f_phase]*sin(d_theta2)+dQ_dTh[ind1Jo,ind2Jo]*sin(d_theta2)-Q[f_PQ,f_phase]*cos(d_theta2))
                    dQn_dTh[ind1J,ind2J] = the_Vn/the_V1*(-dP_dTh[ind1Jo,ind2Jo]*sin(d_theta2)+P[f_PQ,f_phase]*cos(d_theta2)-dQ_dTh[ind1Jo,ind2Jo]*cos(d_theta2)-Q[f_PQ,f_phase]*sin(d_theta2))
                    dPn_dV[ind1J,ind2J] = the_Vn/the_V1**2*(-P[f_PQ,f_phase]*cos(d_theta2)+Q[f_PQ,f_phase]*sin(d_theta2))+the_Vn/the_V1*(-dP_dV[ind1Jo,ind2Jo]*cos(d_theta2)+dQ_dV[ind1Jo,ind2Jo]*sin(d_theta2))
                    dQn_dV[ind1J,ind2J] = the_Vn/the_V1**2*(-P[f_PQ,f_phase]*sin(d_theta2)-Q[f_PQ,f_phase]*cos(d_theta2))+the_Vn/the_V1*(-dP_dV[ind1Jo,ind2Jo]*sin(d_theta2)-dQ_dV[ind1Jo,ind2Jo]*cos(d_theta2))
                else:
                    dPn_dTh[ind1J,ind2J] = the_Vn/the_V1*(-dP_dTh[ind1Jo,ind2Jo]*cos(d_theta2)+dQ_dTh[ind1Jo,ind2Jo]*sin(d_theta2))
                    dQn_dTh[ind1J,ind2J] = the_Vn/the_V1*(-dP_dTh[ind1Jo,ind2Jo]*sin(d_theta2)-dQ_dTh[ind1Jo,ind2Jo]*cos(d_theta2))
                    dPn_dV[ind1J,ind2J] = the_Vn/the_V1*(-dP_dV[ind1Jo,ind2Jo]*cos(d_theta2)+dQ_dV[ind1Jo,ind2Jo]*sin(d_theta2))
                    dQn_dV[ind1J,ind2J] = the_Vn/the_V1*(-dP_dV[ind1Jo,ind2Jo]*sin(d_theta2)-dQ_dV[ind1Jo,ind2Jo]*cos(d_theta2))
    dfP_dTh = dP_dTh+dPn_dTh
    dfQ_dTh = dQ_dTh+dQn_dTh
    dfP_dV = dP_dV+dPn_dV
    dfQ_dV = dQ_dV+dQn_dV
    dfCb_dTh = nmp.concatenate([dfCbr_dTh,dfCbi_dTh])
    dfCb_dV = nmp.concatenate([dfCbr_dV,dfCbi_dV])
    J1 = nmp.concatenate([dfP_dTh,dfP_dV],1)
    J2 = nmp.concatenate([dfQ_dTh,dfQ_dV],1)
    J3 = nmp.concatenate([dfCb_dTh,dfCb_dV],1)
    J = nmp.concatenate([J1,J2,J3])
    #start = timeit,default_timer()
    dx = nmp.linalg.solve(J,-1*fmismatch)
    #stop = timeit,default_timer()
    #Jtime.append(stop[1]-start[1])
    x = x+dx
    x = nmp.reshape(x,(2*NoPQBuses,4))
    V = nmp.array(V)
    Theta = nmp.array(Theta)
    Theta[1:,:] = x[0:NoPQBuses,]
    V[1:,:] = x[NoPQBuses:,]
    V = V.tolist()
    Theta = Theta.tolist()
subprocess.call("/home/pi/Documents/Karim/Unbalanced_PowerFlow-main/memmon_v1")
print("\n")
Vreal = nmp.multiply(V,nmp.cos(Theta))
Vimag = nmp.multiply(V,nmp.sin(Theta))
Vcomplex = Vreal + 1j*Vimag
Theta = nmp.angle(Vcomplex)
V = nmp.absolute(Vcomplex)
V = V*Vbase
Theta = Theta*180/pi
stopt = timeit,default_timer()
#print("Total Jacobian Time: ",sum(Jtime))
print('Total Time: ',stopt[1]-startt[1])
debug = 1