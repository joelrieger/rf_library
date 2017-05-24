import re
import numpy as np

class spar:

    global unitsf
    unitsf={'THZ':1e12,'GHZ':1e9,'MHZ':1e6,'KHZ':1e3,'HZ':1.0}

    def __init__(self,filename):

        self.filename=filename

        self.desc={'unit':'GHz','parameter':'S','format':'MA','R':50.0}
        self.data=[]
        self.freq=[]

        self.readfile(filename)

    def MA2RI(self,mag,ang):
        R=mag*np.cos(ang*np.pi/180.0)
        I=mag*np.sin(ang*np.pi/180.0)
        return R+1j*I
    
    def readfile(self,filename):
        fh=open(filename)

        #test=[]
        freq=[]

        Started=False
        n=0

        for rawline in fh.readlines():
            line=rawline.split('!')[0]

            if line.strip()=='':
                continue
            else:
                if Started==True:
                    numeric=r'[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
                    numarr=[float(tmp) for tmp in re.findall(numeric,line)]

                    if len(numarr)%2==0:
                        NewLine=False
                    else:
                        NewLine=True
                
                    if NewLine==False:
                        if self.desc['format']=='RI':
                            cmplxnumarr=[numarr[2*i]+1j*numarr[2*i+1] for i in range((len(numarr))/2)]
                        elif self.desc['format']=='MA':
                            cmplxnumarr=[self.MA2RI(numarr[2*i],numarr[2*i+1]) for i in range((len(numarr))/2)]

                        self.data[-1].extend(cmplxnumarr)
                    else:
                        if n>0:
                            NumPorts=np.sqrt((len(self.data[-1])))
                            self.data[-1]=np.matrix(self.data[-1]).reshape((NumPorts,NumPorts))

                        self.freq.append(numarr[0]*unitsf[self.desc['unit']]/unitsf['HZ'])
                        if self.desc['format']=='RI':
                            cmplxnumar=[numarr[2*i+1]+1j*numarr[2*i+2] for i in range((len(numarr)-1)/2)]
                        elif self.desc['format']=='MA':
                            cmplxnumar=[self.MA2RI(numarr[2*i+1],numarr[2*i+2]) for i in range((len(numarr)-1)/2)]
                        self.data.append(cmplxnumar)
                        NewLine=False

            if line.strip().startswith('#'):
                self.desc['unit'],self.desc['parameter'],self.desc['format'],nul,self.desc['R']=[tmp.upper() for tmp in line.split()[1::]]
                self.desc['R']=float(self.desc['R'])
                Started=True
                NewLine=True

        #Missed the last one; computation performed when new frequency is read.
        NumPorts=np.sqrt((len(self.data[-1])))
        self.data[-1]=np.matrix(self.data[-1]).reshape((NumPorts,NumPorts))

    def getdata(self,frequency):
        try:
            i=self.freq.index(frequency)
        except ValueError:
            print "ERROR: Interpolation is not yet supported"
            return -1
        else:
            return self.data[i]
