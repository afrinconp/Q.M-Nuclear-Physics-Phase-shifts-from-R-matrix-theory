import matplotlib.pyplot as plt
import mpmath as mp
from scipy.special import eval_legendre 
from scipy.special import roots_legendre  
from scipy.special import erf  
import numpy as np
import pandas as pd


#-----------------------------------------Parameters---------------------------------------#

e2= 1.44
Ap= 1 
At= 12 
zp= 1
zt= 6
h2mu=20.736*(Ap+At)/(Ap*At)   
a=8 #leer
V0=73.8
R0=2.70
l=0

#-------------------------------Legendre roots and function--------------------------------#

def root(x):  #zeros legendre
    zero=roots_legendre(x)[0]
    xi=(zero+1)/2
    return xi

def L(N,x): #eval legendre function
    b=eval_legendre (N,x)
    return b

def f(N): #Regulazed Legendre functions evaluated in a
    M=[]
    zerosN=root(N)
    k=0
    for xi in zerosN:
        p=(-1)**(N+k+1)*(1/xi)*(a*xi*(1.-xi))**(1./2)*1/(a-a*xi)
        k=k+1
        M.append(p)
    MM=np.array(M)
    return MM



#----------------------------------------Potential--------------------------------------------#

def V(x,l):
    aa=Vn(x)
    bb=Vc(x)
    cc=Vl(x,l)
    return aa+bb +cc

def Vl(x,l):
    cent=h2mu*l*(l+1)/(x*x)
    return cent

def Vn(x):
    xx=(x/R0)**2
    VN=-V0*np.exp(-xx)
    return VN

"""
    if x>R0:
        a=4*e/x
    else:
        a=4*(e/(2*R0))*(3-np.power(x/R0,2))
    return a
"""
def Vc(x):
    a=zp*zt*e2/x
    return a
#--------------------------------------Matrix-------------------------------------------------#


def Cmatrix(zeros1,zeros2,a,N,E): #Matrix Elements
	A = np.matrix(np.zeros((N, N), dtype = np.float))
	N=float(N)
	i=0
	for xi in zeros2:
		j=0
		for xj in zeros1:
			if xi!=xj:
				element=((-1)**(i+j)/(a**2*(xi*xj*(1.-xi)*(1.-xj))**(0.5)))*(N**2+N+1+(xi+xj-2*xi*xj)/(xi-xj)**2-1./(1.-xi)-1./(1.-xj))
				A[i,j]=element*h2mu
			elif xi==xj:
				A[i,j]=((4*N**2+4*N+3)*xi*(1-xi)-6*xi+1)/(3*a**2*xi**2*(1-xi)**2)*h2mu+V(a*xi,l)-E
			j=j+1
		i=i+1	
	return A

#--------------------------------------R-Matrix-------------------------------------------------#
def RE(x,N):
	zeros1=root(N)
	zeros2=root(N)
	Regulazed_Legendre_functions=f(N)
	C=Cmatrix(zeros1,zeros2,a,N,x)                              #C matrix
	Cinv = np.linalg.inv(C)                                     #C inverse matrix
	RE1=np.dot(Cinv,Regulazed_Legendre_functions)               #Dot matrix vector
	p=RE1.tolist()[0]
	RE2=np.array(p)
	REa=(h2mu/a)*Regulazed_Legendre_functions.dot(RE2)          #R-matrix
	return REa
#------------------------------------Scattering matrix------------------------------------------#
def F(l,eta,x):                                   #F Coulumb function
    Fc=mp.coulombf(l,eta,x)
    return Fc

def G(l,eta,x):                                   #G Coulumb function
    Gc=mp.coulombg(l,eta,x)
    return Gc

def derivateH(m,l,eta,x):                         #Derivate of a function
    dx=1e-10
    h=(H(m,l,eta,x+dx/2)-H(m,l,eta,x-dx/2))/dx
    return h

def H(m,l,eta,x):                                 #Coulumb Hankel (+,-) function
    if m==1:
        h=G(l,eta,x)+1j*F(l,eta,x)
    else:
        h=G(l,eta,x)-1j*F(l,eta,x)
    return h
	
def U(l,eta,x,E,N):                                 #Scattering matrix
    ul=(H(0,l,eta,x)-x*RE(E,N)*derivateH(0,l,eta,x))/(H(1,l,eta,x)-x*RE(E,N)*derivateH(1,l,eta,x))
    return ul 

#---------------------------------------Phase shifts (in degrees)------------------------------------------#

def angle(i,l,N):
	k=(i/h2mu)**(1./2)
	eta=zp*zt*e2/(2*k*h2mu)
	scatterM=U(l,eta,k*a,i,N)
	y=float(scatterM.imag)
	x=float(scatterM.real)
	alpha=np.arctan2(y,x)/2
	if alpha<0:
		alpha=alpha+np.pi
	deg=np.degrees(alpha)
	return deg	


#---------------------------------------------------Plot-------------------------------------------------#
E=np.arange(0.01,2,0.02)
Angle=[]
for j in 7,10:
	M=[]
	for i in E:
		shift=angle(i,l,j)
		print(j,i,shift)
		M.append(shift)
	Angle.append(M)

plt.plot(E,Angle[0],label="N=7")
plt.plot(E,Angle[1],label="N=10")
plt.legend()
plt.show()

"""
def DF(l,eta,x):                         #Fs Derivate
	dx=0.000000001
	h=(F(l,eta,x+dx/2)-F(l,eta,x-dx/2))/dx	
	return h

def DG(l,eta,x):                         #Gs Derivate
	dx=0.000000001
	h=(G(l,eta,x+dx/2)-G(l,eta,x-dx/2))/dx	
	return h

"""

"""

def Pl(l,eta,x,E):
	f=F(l,eta,x)
	g=G(l,eta,x)
	P=(x)/(f**2+g**2)
	return P

def Sl(l,eta,x,E):
	f=F(l,eta,x)
	g=G(l,eta,x)
	ff=DF(l,eta,x)
	gg=DG(l,eta,x)
	S=(x)/(f**2+g**2)*(f*ff+g*gg)
	return S
	
def angle(i,N):
	k=(i/h2mu)**(1./2)
	B=0
	eta=zp*zt*e2/(2*k*h2mu)
	uu=-1*np.arctan(0.005)
	u1=float(F(l,eta,k*a)/G(l,eta,k*a))
	u11=np.arctan(u1)
	u2=float((Pl(l,eta,k*a,i)*RE(i,N))/(1-(Sl(l,eta,k*a,i)-B)*RE(i,N)))
	phi=-np.arctan(u1)+np.arctan(u2)
	return phi
	
E=np.arange(0.01,2,0.02)	
ANGLE=[]
for i in E:
	shift=angle(i,7)
	if shift<0:
		shift=shift+np.pi
	shift1=np.degrees(shift)
	print(i,shift1)
	ANGLE.append(shift1)

ANGLE2=[]
for i in E:
	shift=angle(i,10)
	if shift<0:
		shift=shift+np.pi
	shift1=np.degrees(shift)
	print(i,shift1)
	ANGLE2.append(shift1)

plt.title("Alpha shift")
plt.xlabel("E(MeV)")
plt.plot(E,ANGLE,label="N=7")
plt.plot(E,ANGLE2,label="N=10")
plt.legend()
plt.show()

"""



	

