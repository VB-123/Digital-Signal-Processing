#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
  float real;
  float imag;
}complex;

complex add(complex a, complex b) {
    complex c;
    c.real = a.real + b.real;
    c.imag = a.imag + b.imag;
    return c;
}

complex sub(complex a, complex b) {
    complex c;
    c.real = a.real - b.real;
    c.imag = a.imag - b.imag;
    return c;
}

complex mul(complex a, complex b) {
    complex c;
    c.real = (a.real * b.real) - (a.imag * b.imag);
    c.imag = (a.real * b.imag) + (a.imag * b.real);
    return c;
}

void print_complex(complex c) {
    if (c.imag >= 0) {
        printf("%.6f + %.6fi", c.real, c.imag);
    } else {
        printf("%.6f - %.6fi", c.real, -c.imag);
    }
}

int bit_reverse(int num, int N) {
    int reversed = 0;
    for (int i = N; i > 0; i--) {
        reversed = (reversed << 1) | (num & 1);
        num >>= 1;
    }
    return reversed;

}

complex twiddle_factor(int k,int N){
   complex w = {cos(2*M_PI*k/N), -sin(2*M_PI*k/N)};
   return w;
}

complex twiddle_factor_ifft(int k, int N){
   complex w = {cos(2*M_PI*k/N), sin(2*M_PI*k/N)};
  return w;
}

void butterfly(complex* a, complex* b, complex w) {
    complex temp = *b;
    *b = sub(*a, mul(w, temp));
    *a = add(*a, mul(w, temp));
}

int main() {
    int N;
    printf("Enter N, for N-point fft (must be a power of 2): ");
    scanf("%d", &N);

    complex* arr = (complex*)malloc(N * sizeof(complex));
    complex* result = (complex*)malloc(N * sizeof(complex));
    complex* ifft = (complex*)malloc(N * sizeof(complex));

    // Initialize input array (example values)
    for(int i = 0; i < N; i++) {
        arr[i].real = i + 1;
        arr[i].imag = 0;
    }

    printf("Input array:\n");
    for(int i = 0; i < N; i++) {
        print_complex(arr[i]);
        printf("\n");
    }

    // Copy input to result array
    for(int i = 0; i < N; i++) {
        result[i] = arr[i];
    }

    // FFT Computation
    complex w;
    for (int s = N; s > 1; s /= 2) { // The stage numbers
        int m = s / 2;
        for (int k = 0; k < N ; k += s) { // The groups in each stage
            for (int j = 0; j < m; j++) { // The inner Butterfly
                w = twiddle_factor(j*N/s,N);
                complex t = result[k + j];
                complex u = result[k + j + m];
                result[k + j] = add(t, u);
                complex temp = sub(t, u);
                result[k + j + m] = mul(w, temp);
            }
        }
    }

    // Generate bit-reverse array
    int bits = (int)log2(N);
    for(int i = 0; i < N; i++) {
        int rev = bit_reverse(i, bits);
        if (i < rev) {
            complex temp = result[i];
            result[i] = result[rev];
            result[rev] = temp;
        }
    }

    printf("FFT Result:\n");
    for (int i = 0; i < N; i++) {    
      print_complex(result[i]);
      printf("\n");
    }
    //IFFT

    // Copy result to ifft array
    for(int i = 0; i < N; i++) {        
        ifft[i] = result[i];
    }

    // Generate bit-reverse array
    //int bits = (int)log2(N);
    for(int i = 0; i < N; i++) {
        int rev = bit_reverse(i, bits);
        if (i < rev) {
            complex temp = ifft[i];
            ifft[i] = ifft[rev];
            ifft[rev] = temp;
        }
    }

    // IFFT Computation
    for(int i = 2; i <= N; i *= 2) { // i = 2^s, s = 1,2... is the stage number
        for(int k = 0; k < N; k += i) { // k for groups in each stage
            for(int j = 0; j < i/2; j++) { // j for butterflies in each group
                complex wm = twiddle_factor_ifft(j * N / i, N);
                complex t = mul(wm, ifft[k + j + i/2]);
                complex u = ifft[k + j];
                ifft[k + j] = add(u, t);
                ifft[k + j + i/2] = sub(u, t);
            }
        }
    }

    // Scale the results by 1/N
    for(int i = 0; i < N; i++) {
       ifft[i].real /= N;
       ifft[i].imag /= N;
    }

    // Print results
    printf("IFFT Result:\n");
    for (int i = 0; i < N; i++) {
        print_complex(ifft[i]);
        printf("\n");
    }

    // Free allocated memory


    free(arr);
    free(result);
    return 0;
}
