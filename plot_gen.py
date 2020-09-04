from sympy import *
from scipy import *
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
import gen_T as gT
from numpy import linalg as LA
from statistics import *
from scipy.linalg import eigh
# Note: last 2 functions are still under construction, can only use up to make_isolated

def create_TZ(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype):
    M = gT.rand_mat
    A=array([v,n,p11,p12,p22,mini,maxi,meshpoint,mattype])
    M.setValue(M,A)
    T = M.make_T(M)
    Z=np.linspace(mini,maxi,meshpoint)
    return T,Z

def make_rho(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype):
    T=create_TZ(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype)[0]
    N1=int(v*n)
    N2=n-N1
    NIT11=N1*T[:, array([1]) ]
    NIT12=N2*T[:, array([3]) ]
    NIT1=concatenate((NIT11,NIT12),axis=1)
    rho=(-1/(n*pi))*NIT1.sum(axis=1)
    return rho


def make_plot_rho(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype,col):
    Rho = make_rho(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype)
    Z=create_TZ(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype)[1]
    plt.plot(Z,abs(Rho),col,linewidth=3.0)
    
def make_plot_t(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype,plottype):
    T,Z = create_TZ(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype)
    if plottype=='real_t1':
        k=0
    elif plottype=='im_t1':
        k=1
    elif plottype=='real_t2':
        k=2
    elif plottype=='im_t2':
        k=3
    plt.plot(Z,T[:,k])

def make_plot_all(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype):
    T,Z=create_TZ(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype)
    color = ['r','g','b','m']
    typ=['$a_{1}(z)$','$b_{1}(z)$','$a_{2}(z)$','$b_{2}(z)$']
    for i in arange(4):
        plt.plot(Z,T[:,i],color[i],label=typ[i])
        plt.title('s = ' + str(v))
        plt.xlabel('z')
        plt.legend(loc='upper left')  

def make_figure(A,n,p11,p12,p22,mini,maxi,meshpoint,mattype):
    for a in A:
        make_plot_all(a,n,p11,p12,p22,mini,maxi,meshpoint,mattype)
        plt.show()

def make_empirical_eig(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype):
    M = gT.rand_mat
    A=array([v,n,p11,p12,p22,mini,maxi,meshpoint,mattype])
    M.setValue(M,A)
    E,V = eigh( M.make_empirical(M) )
    return E,V

def make_emp_fig(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype):
    Z=create_TZ(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype)[1]
    E=make_empirical_eig(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype)[0]
    #if (mattype=='adjacency') or (mattype=='laplacian'):
    #    numbin=maxi-mini
    #elif mattype=='Nlaplacian':
    #    numbin=200
    hist, edges = np.histogram(E, bins=int(len(Z)), range=(mini,maxi), normed=None, weights=None, density=True)
    plt.scatter(Z,hist,marker='v',color='g')
    make_plot_rho(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype,'r')
    
def make_all_rhos(A,n,p11,p12,p22,mini,maxi,meshpoint,mattype):
    color = ['r','g','b']
    typ=['s='+str(A[0]),'s='+str(A[1]),'s='+str(A[2])]
    i=0
    Z=create_TZ(A[0],n,p11,p12,p22,mini,maxi,meshpoint,mattype)[1]
    plt.title('Spectrum ' +  r"$\rho(z)$" + ' of M')
    for a in A:
        Rho=make_rho(a,n,p11,p12,p22,mini,maxi,meshpoint,mattype)
        plt.plot(Z,abs(Rho),color[i],label=typ[i])
        plt.xlabel('z')
        plt.legend(loc='upper left')
        i=i+1
        
def make_isolated(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype,numbin,xrange,yrange):
    N1=int(v*n)
    N2=n-N1
    d1=N1*p11+N2*p12
    d2=N1*p12+N2*p22
    T,Z = create_TZ(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype)
    E,V=make_empirical_eig(v,n,p11,p12,p22,mini,maxi,meshpoint,mattype)
    #if (mattype=='adjacency') or (mattype=='laplacian'):
    #    numbin=maxi-mini
    #elif mattype=='Nlaplacian':
    #    numbin=200
    a1=T[:,0]
    a2=T[:,2]
    if mattype=='adjacency':
        D=(1-N1*a1*p11)*(1-N2*a2*p22)-N1*N2*a1*a2*(p12**2)
    elif mattype=='laplacian':
        D=(1+N1*a1*p11)*(1+N2*a2*p22)-N1*N2*a1*a2*(p12**2)
    elif mattype=='Nlaplacian':
        D=(1+(N1*a1*p11)/d1)*(1+(N2*a2*p22)/d2)-(N1*N2*a1*a2*(p12**2))/(d1*d2)
    C=np.zeros(len(Z))
    plt.hist(E, bins=numbin, range=( mini,maxi ),fc='g',ec='black',density=None)
    plt.plot(Z,D,'r',linewidth=3)
    plt.plot(Z,C,'b')
    plt.xlim(xrange[0],xrange[1])
    plt.ylim(yrange[0],yrange[1])

def make_iso_eigvec(v,n,p11,p12,p22,mini,maxi,mattype):
    N1=int(v*n)
    N2=n-N1
    d1=N1*p11+N2*p12
    d2=N1*p12+N2*p22
    T,Z = create_TZ(v,n,p11,p12,p22,mini,maxi,mattype)
    E,V=make_empirical_eig(v,n,p11,p12,p22,mini,maxi,mattype)
    sorted_E_ids = np.argsort(E)
    E = E[sorted_E_ids]
    V = V[:,sorted_E_ids]
    if (mattype=='adjacency') or (mattype=='laplacian'):
        numbin=maxi-mini
    elif mattype=='Nlaplacian':
        numbin=200
    hist, edges = np.histogram(E, bins=numbin, range=( mini,maxi ), normed=None, weights=None, density=True)
    index_iso=[]
    a1=T[:,0]
    a2=T[:,2]
    if (mattype=='adjacency') or (mattype=='laplacian'):
        for i in arange( len(hist) ):
            if hist[i] == 1:
                index_iso.append(i)
        index_iso.sort()
    elif mattype=='Nlaplacian':
        indl = list(np.nonzero(hist))[0]
        ind2=indl[1]    
    if mattype=='adjacency':
        TB = np.array([[ a1[ index_iso[-2:][0] ], 0 ],[ 0, a2[ index_iso[-2:][0] ] ]])
        MB=np.array([ [p11, p12], [ p12, p22 ] ])
        ind=-2
    elif mattype=='laplacian':
        TB = np.array([[ a1[ index_iso[:2][1] ], 0 ],[ 0, a2[ index_iso[:2][1] ] ]])
        MB=np.array([ [-p11, -p12], [ -p12, -p22 ] ])
        ind=1
    elif mattype=='Nlaplacian':
        TB = np.array([[ a1[ ind2 ], 0 ],[ 0, a2[ ind2 ] ]])
        MB=np.array([ [-p11/d1, -p12/np.sqrt(d1*d2)], [ -p12/np.sqrt(d1*d2), -p22/d2 ] ])
        ind=1
    NB=np.array([ [ N1,0 ],[ 0,N2 ] ])
    K=np.dot(TB, np.dot( MB,NB ) )
    e1,v1=LA.eig(K)
    if abs(e1[0]-1)>0.1 and abs(e1[1]-1)>0.1:
        ind2=indl[0]
    TB = np.array([[ a1[ ind2 ], 0 ],[ 0, a2[ ind2 ] ]])
    K=np.dot(TB, np.dot( MB,NB ) )
    e1,v1=LA.eig(K)
    for i in np.arange(len(e1)):
        if (v1[0,i]*v1[1,i]<0):
            ind3=i
    vec=v1[:,ind3]
    w1=vec[0]*np.ones(N1)/np.sqrt(N1)
    w2=vec[1]*np.ones(N2)/np.sqrt(N2)
    w=np.concatenate((w1,w2),axis=0)
    return w,V[:,ind],abs(np.dot(V[:,ind],w))**2

def plot_iso_eigvec(v,n,p11,p12,p22,mini,maxi,mattype):
    w,V=make_iso_eigvec(v,n,p11,p12,p22,mini,maxi,mattype)[:2]
    if np.dot(w,V)<0:
        w=-w
    plt.figure(figsize=(6,6))
    plt.plot(V,'+',label='emprically observed')
    plt.plot(w,'r',label='asymptotic',linewidth=3)
    plt.title('ASYMPTOTIC AND EMPRICALLY OBSERVED EIGENVECTORS')
    plt.xlabel('Node index, i')
    plt.ylabel('Eigenvalue entry  '+ r'$ v_{i}$')
    plt.legend(loc='lower left')

def detect(v,n,dp,pout,number,mini,maxi,mattype):
    pin=pout+dp
    Dprod=np.zeros(number)
    for i in np.arange(number):
        Dprod[i]=make_iso_eigvec(v,n,pin[i],pout[i],pin[i],mini,maxi,mattype)[2]
    plt.plot(dp,Dprod)