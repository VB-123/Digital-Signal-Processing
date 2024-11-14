#Overlap Add method
import numpy as np 

def circular_convolution(x1, x2):
    N = max(len(x1), len(x2))
    x1_padded = np.pad(x1, (0, N - len(x1)), mode = "constant")
    x2_padded = np.pad(x2, (0, N - len(x2)), mode = "constant")
    y = np.fft.ifft(np.fft.fft(x1_padded)*np.fft.fft(x2_padded))
    return y

def Overlap_Add(x, h, N):
    M = len(h)
    L = N - M + 1 

    x_padded = np.pad(x, (0, (L - len(x)%L)), mode = "constant")
    #Split X into blocks
    x_split = np.split(x_padded, len(x_padded) // L)

    #Pad the impulse sequence
    h_padded = np.pad(h, (0, N - M), mode = "constant")

    #Declare the output array
    y = np.zeros(len(x_padded) + M - 1, dtype = "complex")

    for i, x_i in enumerate(x_split):
        x_i_padded = np.pad(x_i, (0, N - L), mode = "constant")
        y_i = circular_convolution(x_i_padded, h_padded)
        y[i*L : i*L+N] += y_i

    return y 

x = np.array([1+1j, 2+2j, 3+3j, 4+4j])
h = np.array([1+1j, 1+1j])
print(Overlap_Add(x, h, 4))
