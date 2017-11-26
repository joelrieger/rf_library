import numpy as np

text="""R:R1  2 1 50 
R:R2  0 2 25 
R:R3  3 2 12 
R:R4  0 3 12
Vsrc:V1 3 0 2"""


text="""Vsrc:V1 0 1 1
R:R1 1 2 5
R:R2 0 2 10
Isrc:I1 0 2 1"""

text="""Vsrc:V1 0 1 1
Vsrc:V2 3 4 1
R:R1 1 2 10
R:R4 1 4 10
R:R3 2 3 10
R:R2 0 2 10"""


ElemDict={'R','Isrc','Vsrc'}

#Create array from nets
RawNets=[]
RawNetsV=[]
RawNetsI=[]
nets=set([])
for i,line in enumerate(text.split('\n')):
    element=list(line.split())
    if element[0].startswith('R'):
        RawNets.append(([element[0],float(element[1]),float(element[2]),float(element[3])]))
        nets.add(float(element[1]))
        nets.add(float(element[2]))
    elif element[0].startswith('Vsrc'):
        RawNetsV.append(([element[0],float(element[1]),float(element[2]),float(element[3])]))
    elif element[0].startswith('Isrc'):
        RawNetsI.append(([element[0],float(element[1]),float(element[2]),float(element[3])]))
    if len(element)>4:
        print "ERROR: Line "+str(i+1)+", too many arguments"
        break

elems=zip(*RawNets)

#Count Nets
M=len(RawNets)
M=len(nets)
N=len(RawNetsV) #Sources

A=np.matrix(np.zeros([M+N,M+N]))

for typ,n1,n2,val in RawNets:
    if n1!=0:
        A[n1-1,n1-1]+=1.0/val
    if n2!=0:
        A[n2-1,n2-1]+=1.0/val
    if ((n2!=0)&(n1!=0)):
        A[n1-1,n2-1]+=-1.0/val
        A[n2-1,n1-1]+=-1.0/val

z=np.matrix(np.zeros([M+N,1]))
for i,(typ,n1,n2,val) in enumerate(RawNetsI):
    if n1!=0:
        z[n1-1]-=val
    if n2!=0:
        z[n2-1]+=val

for i,(typ,n1,n2,val) in enumerate(RawNetsV):
    if n2!=0:
        A[n2-1,M]=1
        A[M,n2-1]=1
    if n1!=0:
        A[n2-1,M]=1
        A[M,n2-1]=1
    z[M+i]=val

print np.linalg.inv(A)*z



