import re
import numpy as np
import serial
import math
import scipy.linalg as LA
import scipy.signal as ss
L = 1  # number of sources
N = 2  # number of ULA elements
M=2

arduino = serial.Serial('COM3', 115200)
def doa(data1,data2):
    angle=0
    tt=0
    pr=0
    amp2=0
    for k in range(len(data1)):
        amp=math.sqrt(data1[k]*data1[k]+data2[k]+data2[k])
        faza=math.degrees(math.atan2(data1[k],data2[k]))
        time=(1/data1[k])*(faza/360)
        x=time*(331*0.606*20)
        tt+=time*amp
        pr+=x*amp
        amp2+=amp
    r0=pr/amp2
    angle=math.asin(r0)/0.2
    x = math.cos(angle) * r0
    y = math.sin(angle) * r0

    print(f"angle: {angle}, (x, y): ({x}, {y})")

    return (x, y)
def esprit(CovMat,L,N):
    # CovMat is the signal covariance matrix, L is the number of sources, N is the number of antennas
    _,U = LA.eig(CovMat)
    S = U[:,0:L]
    Phi = LA.pinv(S[0:N-1]) @ S[1:N] # the original array is divided into two subarrays [0,1,...,N-2] and [1,2,...,N-1]
    eigs,_ = LA.eig(Phi)
    DoAsESPRIT = np.degrees(np.arcsin(np.angle(eigs)/np.pi))
    #DoAsESPRIT = DoAsESPRIT*180
    return DoAsESPRIT
try:
    while True:
        for i in range(10):
            print(arduino.readline())
            arduinoString = arduino.readline().decode()

            dataArray = arduinoString.split()
        byteObject = re.findall(r'\d+', arduinoString)

        mic0 = [int(row[0]) for row in byteObject]
        mic1 = [int(row[1]) for row in byteObject]
        s = []
        midd = len(byteObject) // 2
        SignalMAtrix = np.array([mic0, mic1])
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
        print('ESPRIT DoAs:', DoAsESPRIT, '\n')
        #for k in range(len(byteObject)):
            #s.append(np.fft.fft(byteObject[k]))

        DOA1 = doa(mic0, mic1)

except KeyboardInterrupt:
    arduino.close()

print("serial connection closed")