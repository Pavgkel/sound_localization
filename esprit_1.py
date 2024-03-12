import socket
import ipaddress
import numpy as np
import scipy.linalg as LA
import scipy.signal as ss
import math
import re
import csv

TIMEOUT = 1
L = 1  # number of sources
N = 2  # number of ULA elements
M=2
array = np.linspace(0,(N-1)/2,N)
Angles = np.linspace(-np.pi/2,np.pi/2,360)
output_file = open("angle.csv", "w")
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
sock.listen(1)

while True:
    conn, addr = sock.accept()
    conn.settimeout(TIMEOUT)
    print("Accepted from", addr)

    size = conn.recv_into(buf)
    print(buf[:size])
    byteObject = buf[:size].decode()
    conn.close()

    A = re.findall(r'\d+', byteObject)
    A = [int(row) for row in A]
    if (len(A) % 2 != 0):
        A.pop(-1)
    midd = len(A) // 2
    SignalMAtrix = np.array([A[:midd], A[midd:]])
    K = math.floor(SignalMAtrix.shape[1] / N)  # number of snapshots
    NS = K * N
    SignalMAtrix = SignalMAtrix[:, :NS]
    SignalMAtrix = np.flipud(SignalMAtrix)
    X = np.zeros(SignalMAtrix.shape, 'complex')
    A = []

    for k in range(0, K):
        for m in range(0, M):
            s_k = SignalMAtrix[m][k * N:k * N + N]
            S = np.fft.fft(s_k)
            X[m][k:K * N:K] = S

    CovMat = X @ X.conj().transpose()

    DoAsESPRIT = esprit(CovMat, L, N)

    # DoAsMUSIC, psindB = music(CovMat, L, N, array, Angles)
    # DOA1 = doa(SignalMAtrix[:midd], SignalMAtrix[midd:])
    # DOA1=doa(byteObject[:midd],byteObject[midd:])
    # DOA1 = calc_coord(A, B, C, D)
    # print('MUSIC DoAs:', DoAsMUSIC, '\n')
    print('ESPRIT DoAs:', DoAsESPRIT, '\n')
    # print(list(DoAsESPRIT))
    #coordinates = cordinate(x1, x2, x3, x4, y1, y2, y3, y4, DoAsESPRIT[0], DoAsESPRIT[1], DoAsESPRIT[2], DoAsESPRIT[3])
    writer.writerow([DoAsESPRIT])


