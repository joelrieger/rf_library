from spar_class import *
from scipy import optimize as opt
from scipy.interpolate import interp1d
import numpy as np
import copy
import sys, traceback #for error handling

def SSS(start,stop,step):
    return [start+step*i for i in range(int((stop-start)/step)+1)]

def passbyval(func):
    def new(*args,**kargs):
        cargs = [copy.copy(arg) for arg in args]
        ckargs = copy.copy(kargs)
        return func(*cargs,**ckargs)
    return new

class ErrorCatch(LookupError):
    def __init__(self,ErrCode,string):
        print "Error Code:", ErrCode
        print string
        sys.exit()


@passbyval
def stoz(S,Z0=50.0):

    Z0=complex(Z0)
    
    #Error Check
    if not (isinstance(S,np.matrix)):
        tb=sys.exc_info()[-1]
        raise ErrorCatch(0,'stoz requires numpy matrix, got '+str(type(a)))

    n,m=S.shape
    if n!=m:
        raise ErrorCatch(0,'square matrix required, got '+str(S.shape))

    E=np.identity(n)
    sqrtzn=np.sqrt(Z0)*E
    Z=sqrtzn*(E+S)*(E-S)**(-1)*sqrtzn

    return Z


def zC(C,freq):
    return 1.0/(2*np.pi*freq*C*1j)


def zL(L,freq):
    return (2*np.pi*freq*L*1j)


@passbyval
def zin(Z,InputPort=0,ZTerm=[]):
    #Check Zterm length and Z size
    n,m=Z.shape

    V=np.zeros((n,1))
    V[InputPort-1,0]=1 #Apply a 1V source to determine input impedance

    ZTerm[InputPort-1]=0
    for i,zval in enumerate(ZTerm):
        Z[i,i]+=zval

    I=Z**(-1)*V

    return (I[InputPort-1,InputPort-1])**(-1)

#@passbyval
def TruncateZ(Z,NewPorts=[1,2],ZTerm=[]):
    """Given a Znm matrix, use ZTerm to apply impedances to ports not in NewPorts"""
    N,M=Z.shape

    NewZ=np.zeros((len(NewPorts),len(NewPorts)),np.complex)

    for m in [tmp+1 for tmp in range(M)]:
        Zx=copy.copy(Z)
        Vx=np.zeros((N,1),np.complex)
        for n in [tmp+1 for tmp in range(N)]:
            if n==m:
                Vx[n-1,0]=1
            elif n not in NewPorts:
                Zx[n-1,n-1]+=ZTerm[n-1]
            else:
                Zx[n-1,n-1]+=10e6
                
        I=Zx**(-1)*Vx
        V=Z*I
        
        for n in NewPorts:
            Newn=NewPorts.index(n)
            if m in NewPorts:
                Newm=NewPorts.index(m)
                NewZ[Newn,Newm]=(V[n-1]/I[m-1])[0,0]

    return NewZ

#@passbyval
def NWLoss(Z,Z0=50.0):

    Z11=Z[0,0]
    Z12=Z[0,1]
    Z21=Z[1,0]
    Z22=Z[1,1]

    def convert(Z11,Z12,Z21,Z22):
        dZ=(Z11+Z0)*(Z22+Z0)-Z12*Z21

        S11=((Z11-Z0)*(Z22+Z0)-Z12*Z21)/dZ
        S12=(2*Z12*Z0)/dZ
        S21=(2*Z21*Z0)/dZ
        S22=((Z11+Z0)*(Z22-Z0)-Z12*Z21)/dZ

        return [S11,S12,S21,S22]

    S11,S12,S21,S22=convert(Z11,Z12,Z21,Z22)

    return 10*np.log10(abs(S21)**2/(1-abs(S11)**2))

if __name__=='__main__':  

    sdata=spar("TEST8PORT.s8p")
    S=sdata.getdata(1000e6)
    S=np.matrix(S)
    Zarr=[]
    Z=stoz(S)
