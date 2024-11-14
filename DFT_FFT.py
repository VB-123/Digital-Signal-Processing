import time 
import numpy as np 
import matplotlib.pyplot as plt 

def dft_matrix(N):
    w_n = np.exp(-2j*np.pi/N)
    return [[w_n ** (j*k) for j in range (N)] for k in range (N)]

def compute_dft(x, N):
    M = dft_matrix(N)
    return np.dot(M,x)

def dft_time(x, N):
    start = time.time()
    _ = compute_dft(x, N)
    end = time.time()
    return end - start

def fft_time(x, N):
    start = time.time()
    _ = np.fft.fft(x,N)
    end = time.time()
    return end - start

dft_time_values, fft_time_values = [],[]
values = [2**j for j in range(2,12)]
for val in values:
    x = np.random.rand(val)
    dft_time_values.append(dft_time(x,val))
    fft_time_values.append(fft_time(x,val))

plt.plot(values,dft_time_values,label = "DFT")
plt.xlabel("Values of N")
plt.plot(values,fft_time_values,label = "FFT")
plt.ylabel("Times")
plt.grid()
plt.figure(figsize = (12,5))
plt.show()

