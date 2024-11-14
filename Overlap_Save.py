#Overlap save
import numpy as np

def circular_convolution(x1, x2):
    N = max(len(x1), len(x2))
    x1_padded = np.pad(x1, (0, N -len(x1)), mode = "constant")
    x2_padded = np.pad(x2, (0, N -len(x2)), mode = "constant")

    y = np.fft.ifft(np.fft.fft(x1_padded) * np.fft.fft(x2_padded))
    return y

def overlap_save(x, h, N):
    M = len(h)
    L = N - M + 1 

    x_padded = np.pad(x, (M-1,0), mode = "constant")
    h_padded = np.pad(h, (0, N-M), mode = "constant")
    y = np.array([], dtype = "complex")
    for i in range(0, len(x_padded), L):
        block = x_padded[i : i+N]
        if (len(block) < N):
            block = np.pad(block, (0, N - len(block)), mode = "constant")
        ym = circular_convolution(block, h_padded)
        ym = ym[M-1:]
        y = np.append(y,ym)
    return y
x = np.array([1+1j, 2+2j, 3+3j, 4+4j])
h = np.array([1+1j, 1+1j])
print(overlap_save(x, h, 4))
