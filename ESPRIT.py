import socket
import ipaddress
import numpy as np
import scipy.linalg as LA
import scipy.signal as ss
import math
import re
import csv

TIMEOUT = 1
L = 4  # number of sources
N = 2  # number of ULA elements
M=4
array = np.linspace(0,(N-1)/2,N)
Angles = np.linspace(-np.pi/2,np.pi/2,360)
output_file = open("ESPRIT.csv", "w")
writer=csv.writer(output_file)
output2 = open("Coordinates.csv", "w")
writer2=csv.writer(output2)
def array_response_vector(array,theta):
    N = array.shape
    v = np.exp(1j*2*np.pi*array*np.sin(theta))
    return v/np.sqrt(N)

def music(CovMat,L,N,array,Angles):
    # CovMat is the signal covariance matrix, L is the number of sources, N is the number of antennas
    # array holds the positions of antenna elements
    # Angles are the grid of directions in the azimuth angular domain
    _,V = LA.eig(CovMat)
    Qn  = V[:,L:N]
    numAngles = Angles.size
    pspectrum = np.zeros(numAngles)
    for i in range(numAngles):
        av = array_response_vector(array,Angles[i])
        pspectrum[i] = 1/LA.norm((Qn.conj().transpose()@av))
    psindB    = np.log10(10*pspectrum/pspectrum.min())
    DoAsMUSIC,_= ss.find_peaks(psindB,height=1.35, distance=1.5)
    return DoAsMUSIC,pspectrum

def esprit(CovMat,L,N):
    # CovMat is the signal covariance matrix, L is the number of sources, N is the number of antennas
    _,U = LA.eig(CovMat)
    S = U[:,0:L]
    Phi = LA.pinv(S[0:N-1]) @ S[1:N] # the original array is divided into two subarrays [0,1,...,N-2] and [1,2,...,N-1]
    eigs,_ = LA.eig(Phi)
    DoAsESPRIT = np.degrees(np.arcsin(np.angle(eigs)/np.pi))
    #DoAsESPRIT = DoAsESPRIT*180
    return DoAsESPRIT

def cordinate(x1,x2,x3,x4,y1,y2,y3,y4,angle1,angle2,angle3,angle4):
    X=x1*math.tan(angle1)-x2*math.tan(angle2)-x3*math.tan(angle3)-x4*math.tan(angle4)-y1+y2-y3+y4
    Y=y1+(X-x1)*math.tan(angle1)
    return X,Y

def doa(data1,data2):
    angle=0
    tt=0
    pr=0
    amp2=0
    for k in range(len(data1)):
        amp=math.sqrt(data1[k]*data1[k]+data2[k]+data2[k])
        faza=math.degrees(math.atan2(data1[k],data2[k]))
        time=(1/data1[k])*(faza/360)
        x=time*(331*0.606*21)
        tt+=time*amp
        pr+=x*amp
        amp2+=amp
    r0=pr/amp2
    angle=math.asin(r0)/0.2
    x = math.cos(angle) * r0
    y = math.sin(angle) * r0

    print(f"angle: {angle}, (x, y): ({x}, {y})")

    return (x, y)
    #return angle

HOST = '192.168.3.11'
# print(str(wifi.radio.ipv4_address))
PORT = 80  # Port to listen on

print("Creating socket")
sock = socket.socket()
buf_size=256*9*2
buf = bytearray(buf_size)

sock.bind((HOST, PORT))
sock.listen(4)
A=[]
B=[]
C=[]
D=[]

def calc_coord(A, B, C, D):
    midd = len(A) // 2
    A_coord = doa(A[:midd], A[midd:])

    midd = len(B) // 2
    B_coord = list(doa(B[:midd], B[midd:]))
    B_coord[0] = 0.4 - B_coord[0]

    midd = len(C) // 2
    C_coord = list(doa(C[:midd], C[midd:]))
    C_coord[1] = 0.4 - C_coord[1]

    midd = len(D) // 2
    D_coord = list(doa(D[:midd], D[midd:]))
    D_coord[0] = 0.4 - D_coord[0]
    D_coord[1] = 0.4 - D_coord[1]

    return ((A_coord[0] + B_coord[0] + C_coord[0] + D_coord[0]) / 4,
            (A_coord[1] + B_coord[1] + C_coord[1] + D_coord[1]) / 4)

x1=0
x2=0
x3=50
x4=50
y1=0
y2=50
y3=50
y4=0

while True:
    conn, addr = sock.accept()
    conn.settimeout(TIMEOUT)
    print("Accepted from", addr)

    size = conn.recv_into(buf)
    byteObject = buf[:size].decode()
    conn.close()
    #print(byteObject)
    if byteObject.find('A'):
        A=re.findall(r'\d+', byteObject)
        A = [int(row) for row in A]
        if (len(A)%2!=0):
            A.pop(-1)
    if byteObject.find('B'):
        B=re.findall(r'\d+', byteObject)
        B = [int(row) for row in B]
        if (len(B)%2!=0):
            B.pop(-1)
    if byteObject.find('C'):
        C=re.findall(r'\d+', byteObject)
        C = [int(row) for row in C]
        if (len(C)%2!=0):
            C.pop(-1)
    if byteObject.find('D'):
        D = re.findall(r'\d+', byteObject)
        D = [int(row) for row in D]
        if (len(D)%2!=0):
            D.pop(-1)

    #byteObject=re.findall(r'\d+', byteObject)
    #byteObject = [int(row) for row in byteObject]
    #if (len(byteObject) % 2 != 0):
        #byteObject.pop(-1)
    #byteObject=np.split(byteObject, 2)
    #print(len(byteObject))
    if len(A)!=0 and len(B)!=0:
        midd = len(A) // 2
        midd2 = len(B) // 2
        midd3=len(C) // 2
        midd4=len(D) // 2
        SignalMAtrix = np.array([A[:midd], A[midd:], B[:midd2], B[midd2:],C[:midd3],C[midd3:],D[:midd4],D[midd4:]])
        K = math.floor(SignalMAtrix.shape[1] / N)  # number of snapshots
        NS = K * N
        SignalMAtrix = SignalMAtrix[:, :NS]
        SignalMAtrix = np.flipud(SignalMAtrix)
        X = np.zeros(SignalMAtrix.shape, 'complex')
        A=[]
        B=[]
        C=[]
        D=[]
        for k in range(0, K):
            for m in range(0, M):
                s_k = SignalMAtrix[m][k * N:k * N + N]
                S = np.fft.fft(s_k)
                X[m][k:K * N:K] = S

        CovMat = X @ X.conj().transpose()

        DoAsESPRIT = esprit(CovMat, L, N)

    #DoAsMUSIC, psindB = music(CovMat, L, N, array, Angles)
    #DOA1 = doa(SignalMAtrix[:midd], SignalMAtrix[midd:])
    #DOA1=doa(byteObject[:midd],byteObject[midd:])
        #DOA1 = calc_coord(A, B, C, D)
    #print('MUSIC DoAs:', DoAsMUSIC, '\n')
        print('ESPRIT DoAs:', DoAsESPRIT, '\n')
        #print(list(DoAsESPRIT))
        coordinates=cordinate(x1,x2,x3,x4,y1,y2,y3,y4,DoAsESPRIT[0],DoAsESPRIT[1],DoAsESPRIT[2],DoAsESPRIT[3])
        writer.writerow([DoAsESPRIT])
        print(coordinates)
        writer2.writerow([coordinates])
    #print('(x, y): ', DOA1, '\n')

