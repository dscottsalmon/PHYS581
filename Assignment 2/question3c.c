#include <stdio.h>
#include <math.h>

// Function to compute grid points: x_k = cos(pi/m * (k - 0.5))
// Note that we need to be very careful because C is a 0-index language but equation is 1-index. Hence why x-index is k-1.
void computeGridPoints(int m, double x[]) {
    for (int k = 1; k <= m; k++) {
        x[k-1] = cos(M_PI / m * (k - 0.5));
    }
}

// Chebyshev polynomials of degree n
double ChebyshevT(int n, double x) {
    return cos(n * acos(x));
}

// derivative of Chebyshev polynomials of degree n
double ChebyshevTPrime(int n, double x) {
    return n*sin(n*acos(x))/(sqrt(1-x*x));
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

//Uses quotient rule and mathematica results to find derivative of Equation 1.36
double BarycentricDerivative(int m, double x[], double w[], double f[], double x_eval) {

    double numerator = 0.0, denominator = 0.0;
    double sum_N_prime = 0.0, sum_D_prime = 0.0;

    for (int j = 0; j < m; j++) {
        
        //this block is the same as regular interpolation equation
        double diff = x_eval - x[j];
        double term = w[j] / diff;
        numerator += term * f[j];
        denominator += term;

        //results from mathematica on what dNx and dDx are
        double term_deriv = w[j] / (diff * diff);
        sum_N_prime += term_deriv * f[j];
        sum_D_prime += term_deriv;
    }

    //adding negative signs
    double N_prime = -1 * sum_N_prime;
    double D_prime = -1 * sum_D_prime;

    //line below is applying quotient rule to find derivative
    double derivative = (N_prime * denominator - numerator * D_prime) / (denominator * denominator);

    return derivative;
}

int main() {
    int m = 10;
    double x[10], w[10], f[10];

    computeGridPoints(m, x);
    computeBarycentricWeights(m, x, w);

    for (int n = 1; n <= m; n++) {
        for (int j = 0; j < m; j++) {
            f[j] = ChebyshevT(n, x[j]);
        }

        printf("\nChebyshev Derivative T'_%d(x):\n", n);
        printf("x_eval\t     Exact Deriv\tInterpolated Deriv\n");

        for (double x_eval = -1.0; x_eval <= 1.0; x_eval += 0.2) {
            double exact_deriv = ChebyshevTPrime(n, x_eval);
            double interpolated_deriv = BarycentricDerivative(m, x, w, f, x_eval);
            printf("%.4f\t%12.6f\t%16.6f\n", x_eval, exact_deriv, interpolated_deriv);
        }
    }

    return 0;
}
