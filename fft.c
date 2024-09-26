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

complex twiddle_factor(int N){
   complex w = {cos(2*M_PI/N), -sin(2*M_PI/N)};
  return w;
}

int main() {
    int N;
    printf("Enter N, for N-point fft (must be a power of 2): ");
    scanf("%d", &N);

    complex* arr = (complex*)malloc(N * sizeof(complex));
    complex* result = (complex*)malloc(N * sizeof(complex));

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

    // FFT Computation
    complex w, wm;
    for (int i = 2; i <= N; i = i * 2) { // i is 2^s; s is the stage number
        wm = twiddle_factor(i);
        for (int k = 0; k < N; k += i) { // The groups in each stage 
            w = (complex){1.0, 0.0};
            for (int j = 0; j < i/2; j++) { // The inner Butterfly
                complex t = mul(w, result[k + j + i/2]);
                complex u = result[k + j];
                result[k + j] = add(u, t);
                result[k + j + i/2] = sub(u, t);
                w = mul(w, wm);
            }
        }
    }

    // Print results
    printf("FFT Result:\n");
    for (int i = 0; i < N; i++) {
        printf("%f + %fi\n", result[i].real, result[i].imag);
    }

    // Free allocated memory
    free(arr);
    free(result);

    return 0;
}
