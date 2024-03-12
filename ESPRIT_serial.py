import numpy as np
import scipy.linalg as LA
import scipy.signal as ss
import operator
import serial
import math
import csv

def array_response_vector(array,theta):
    N = array.shape
    v = np.exp(1j*2*np.pi*array*np.sin(theta))
    return v/np.sqrt(N)

def esprit(CovMat,L,N):
    # CovMat is the signal covariance matrix, L is the number of sources, N is the number of antennas
    _,U = LA.eig(CovMat)
    S = U[:,0:L]
    Phi = LA.pinv(S[0:N-1]) @ S[1:N] # the original array is divided into two subarrays [0,1,...,N-2] and [1,2,...,N-1]
    eigs,_ = LA.eig(Phi)
    DoAsESPRIT = np.arcsin(np.angle(eigs)/np.pi)
    DoAsESPRIT = DoAsESPRIT*180
    return DoAsESPRIT

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

np.random.seed(6)

lamda = 1 # wavelength
kappa = np.pi/lamda # wave number
L = 1  # number of sources
N = 1  # number of ULA elements
M=1
snr = 10 # signal to noise ratio
arduino = serial.Serial('COM9', 115200)
output_file = open("ESPRIT.csv", "w")
writer=csv.writer(output_file)
array = np.linspace(0,(N-1)/2,N)
Angles = np.linspace(-np.pi/2,np.pi/2,360)

try:
    while True:
        Q = []
        for i in range(10):
            arduinoString = arduino.readline().decode()
            dataArray = arduinoString.split()
            # print()
            Q.append(dataArray)
        mic0 = [int(row[0]) for row in Q]
        mic1 = [int(row[1]) for row in Q]
        SignalMAtrix = np.array([mic0, mic1], 'float')
        K = math.floor(SignalMAtrix.shape[1] / N)  # number of snapshots
        NS = K * N
        SignalMAtrix = SignalMAtrix[:, :NS]
        SignalMAtrix = np.flipud(SignalMAtrix)
        X = np.zeros(SignalMAtrix.shape, 'complex')

        for k in range(0, K):
            for m in range(0, M):
                s_k = SignalMAtrix[m][k * N:k * N + N]
                S = np.fft.fft(s_k)
                X[m][k:K * N:K] = S

        CovMat = X @ X.conj().transpose()
        DoAsESPRIT = esprit(CovMat, L, N)
        #DoAsMUSIC, psindB = music(CovMat, L, N, array, Angles)
        aa=np.max(DoAsESPRIT)
        if aa.size>0:
            writer.writerow([aa, 1])
        print('ESPRIT DoAs:', aa, '\n')
        #print('DoAsMUSIC DoAs:', DoAsMUSIC, '\n')
except KeyboardInterrupt:
    arduino.close()

print("serial connection closed")