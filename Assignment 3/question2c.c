
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

#define PI 3.141592653589793
#define M 11

// Compute Hermite polynomial using recurrence relation
double hermiteH(int n, double x) {
    if (n == 0) return 1.0;
    if (n == 1) return 2.0 * x;
    
    double Hn_2 = 1.0, Hn_1 = 2.0 * x, Hn;
    for (int i = 2; i <= n; i++) {
        Hn = 2.0 * x * Hn_1 - 2.0 * (i - 1) * Hn_2;
        Hn_2 = Hn_1;
        Hn_1 = Hn;
    }
    return Hn;
}

// Function to evaluate the integrand f(x) = e^(-x^2) * H_j(x) * H_k(x)
double integrand(int j, int k, double x) {
    return hermiteH(j, x) * hermiteH(k, x); // Remove exp(-x*x)
}

// Gaussian Quadrature using Hermite weights and roots (for m=10)
double gaussian_quadrature(int j, int k) {
    
    // Define the roots (x_i) and weights (w_i) for Hermite polynomials (for m=10)
    double x[M] = {-3.43616, -2.53273, -1.75668, -1.03661, -0.342901, 0.342901, 1.03661, 1.75668, 2.53273, 3.43616};
    
    double w[M] = {7.64043286e-6,1.34364575e-3,3.38743945e-2,0.24013861,0.61086263,0.61086263,0.24013861,3.38743945e-2,1.34364575e-3,7.64043286e-6};
    // Compute the integral using Gaussian quadrature
    double integral = 0.0;
    for (int i = 0; i < M; i++) {
        integral += w[i] * integrand(j, k, x[i]);
    }
    
    return integral;
}

int main() {
    printf(" j | k |  Numerical | Exact  | Error\n");
    printf("---------------------------------------\n");
    
    for (int j = 0; j <= 10; j++) {
        for (int k = 0; k <= 10; k++) {
            double exact = (j == k) ? pow(2, j) * tgamma(j + 1) * sqrt(PI) : 0.0; 
            double numerical = gaussian_quadrature(j, k);
            double error = fabs(numerical - exact);

            //changes from printing doubles to scientific notation at j=5
            if (j<5){
                printf("%2d |%2d |%10.5f  |%7.5f |%8.2e \n",
                    j, k, numerical, exact, error);
            } else {
                printf("%2d |%2d |%10.5e |%7.5e |%8.2e \n",
                    j, k, numerical, exact, error); }
        }
    }
    return 0;
}
