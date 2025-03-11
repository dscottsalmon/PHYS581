#include <stdio.h>
#include <math.h>
#include <float.h>

// Compute physicist's Hermite polynomials using recurrence relation
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

// Function to determine xmax for given j, k
double find_xmax(int j, int k) {
    double x = 5.0;  // Initial guess for xmax
    double threshold = 2.22e-16; //machine precision for double is ~2.22e-16
    double step = 0.1;  // Step size for increasing x

    while (1) {
        double H_j = hermiteH(j, x);
        double H_k = hermiteH(k, x);
        double integrand = exp(-x * x) * H_j * H_k;
        
        if (fabs(integrand) < threshold) break;
        
        x += step;
        if (x > 50) break; // Safety limit
    }
    return x;
}

int main() {
    printf("j \\ k |");
    for (int k = 0; k <= 10; k++) printf(" %2d  |", k);
    printf("\n-----------------------------------------------------------------------------\n");

    for (int j = 0; j <= 10; j++) {
        printf(" %2d   |", j);
        for (int k = 0; k <= 10; k++) {
            double xmax = find_xmax(j, k);
            printf(" %.1f |", xmax);
        }
        printf("\n");
    }
    
    return 0;
}

