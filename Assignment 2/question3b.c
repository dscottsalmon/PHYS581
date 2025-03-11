#include <stdio.h>
#include <math.h>

// Function to compute grid points: x_k = cos(pi/m * (k - 0.5))
// Note that we need to be very careful because C is a 0-index language but equation is 1-index. Hence why x-index is k-1.
void computeGridPoints(int m, double x[]) {
    for (int k = 1; k <= m; k++) {
        x[k-1] = cos(M_PI / m * (k - 0.5));
    }
}

// Chebyshev polynomial of degree n
double ChebyshevT(int n, double x) {
    return cos(n * acos(x));
}

// Function to calculate barycentric weights (Equation 1.31)
void computeBarycentricWeights(int N, double x[], double w[]) {
    for (int j = 0; j < N; j++) {
        
        // Initialize the product for w[j]
        double product = 1.0; 
        
        for (int i = 0; i < N; i++) {
            if (i != j) {
                product = product * (x[j] - x[i]); // Compute the product of (x[j] - x[i])
            }
        }
        // Take the reciprocal of product to get w[j]
        w[j] = 1.0 / product; 
    }
}

// Barycentric interpolation function (Equation 1.36)
double BarycentricInterpolate(int m, double x[], double w[], double f[], double x_eval) {
    double numerator = 0.0, denominator = 0.0;
    
    for (int j = 0; j < m; j++) {
        double term = w[j] / (x_eval - x[j]);
        numerator += term * f[j];
        denominator += term;
    }
    return numerator / denominator;
}

int main() {
    
    int m = 10;
    
    //initializing variables
    double x[10], w[10], f[10];

    // Compute grid points using cos(pi/m * (k - 0.5))
    computeGridPoints(m, x);

    // Compute barycentric weights
    computeBarycentricWeights(m, x, w);

    // Test for different Chebyshev polynomials
    for (int n = 1; n <= m; n++) {
        
        // Compute exact values of T_n(x) at grid points
        
        for (int j = 0; j < m; j++) {
            f[j] = ChebyshevT(n, x[j]);
        }

        printf("\nChebyshev Polynomial T_%d(x) Interpolation:\n", n);
        printf("x_eval\t     Exact T_%d(x)\tInterpolated T_%d(x)\n", n, n);
        
        //comparing interpolated and exact data
        for (double x_eval = -1.0; x_eval <= 1.0; x_eval += 0.2) {
            double exact = ChebyshevT(n, x_eval);
            double interpolated = BarycentricInterpolate(m, x, w, f, x_eval);
            printf("%0.4f\t\t%0.4f\t\t    %0.4f\n", x_eval, exact, interpolated);
        }
    }
    return 0;
}

