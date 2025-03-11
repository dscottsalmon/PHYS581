
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

#define PI 3.141592653589793

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

void compute_weights(int n, double w[], double h) {
    if (n == 2) { // Trapezoidal Rule
        w[0] = h/2;
        w[1] = h/2;
    } else if (n == 3) { // Simpson's 1/3 Rule
        w[0] = h/3;
        w[1] = 4*h/3;
        w[2] = h/3;
    } else if (n == 4) { // Simpson's 3/8 Rule
        w[0] = 3*h/8;
        w[1] = 9*h/8;
        w[2] = 9*h/8;
        w[3] = 3*h/8;
    } else if (n == 5) { // n=5 from q1 mathematica
        w[0] = 14*h/45;
        w[1] = 64*h/45;
        w[2] = 8*h/15;
        w[3] = 64*h/45;
        w[4] = 14*h/45; 
    } else if (n == 6) { // n=6 from q1 mathematica
        w[0] = 95*h/288;
        w[1] = 125*h/96;
        w[2] = 125*h/144;
        w[3] = 125*h/144;
        w[4] = 125*h/96; 
        w[5] = 95*h/288; 
    } else if (n == 7) { // n=7 from q1 mathematica
        w[0] = 41*h/140;
        w[1] = 54*h/35;
        w[2] = 27*h/140;
        w[3] = 68*h/35;
        w[4] = 27*h/140; 
        w[5] = 54*h/35; 
        w[6] = 41*h/140; 
    } else if (n == 8) {  //n=8 from q1 mathematica
        w[0] = 5257*h/17280;
        w[1] = 25039*h/17280;
        w[2] = 343*h/640;
        w[3] = 20923*h/17280;
        w[4] = 20923*h/17280; 
        w[5] = 343*h/640; 
        w[6] = 25039*h/17280;
        w[7] = 5257*h/17280;
    } else if (n == 9) { //n=9 from q1 mathematica
        w[0] = 3956*h/14175;
        w[1] = 23552*h/14175;
        w[2] = -3712*h/14175;
        w[3] = 41984*h/14175;
        w[4] = -3632*h/2835; 
        w[5] = 41984*h/14175; 
        w[6] = -3712*h/14175;
        w[7] = 23552*h/14175; 
        w[8] = 3956*h/14175; 
    } else if (n == 10) { //n=10 from q1 mathematica
        w[0] = 25713*h/89600;
        w[1] = 141669*h/89600;
        w[2] = 243*h/2240;
        w[3] = 10881*h/5600;
        w[4] = 26001*h/44800; 
        w[5] = 26001*h/44800; 
        w[6] = 10881*h/5600;
        w[7] = 243*h/2240; 
        w[8] = 141669*h/89600;
        w[9] = 25713*h/89600; 
    } else {
        printf("Error: n-point rule not implemented!\n");
        exit(1);
    }
}

// Function to compute the integral using specified n-point integration rule
double compute_integral(int j, int k, int m, int n, double xmax) {
    
    // Step size for the grid
    double dx = (2.0 * xmax) / (m - 1);  
    
    //initializing integral and weights variables for specified n rule
    double integral = 0.0;
    double weights[n];
    compute_weights(n, weights, dx);

    int step = n - 1; //stepsize for n-sized subset of the m grid
    
    //looping for all the n-sized subsets in m
    for (int i = 0; i <= m - n; i += step) {
        double sub_integral = 0.0;
        double x_center = -xmax + i * dx;

        for (int l = 0; l < n; l++) {
            double x_l = x_center + l * dx;
            double H_j_val = hermiteH(j, x_l);
            double H_k_val = hermiteH(k, x_l);
            sub_integral += weights[l] * exp(-x_l * x_l) * H_j_val * H_k_val;
        }

        integral += sub_integral;
    }

    return integral;
}

int main() {
    int m_values[] = {100, 1000, 10000, 100000}; // Test different m values
    int n_values[] = {2, 3, 4, 5, 6, 7, 8, 9, 10}; // Test different n-point rules
    printf(" j | k |  n  |    m     |  Numerical  |  Exact  |  Error\n");
    printf("------------------------------------------------------------\n");

    //looping over all j and k values 0-10
    for (int j = 0; j <= 10; j++) {
        for (int k = 0; k <= 10; k++) {
            
            double xmax = find_xmax(j, k);
            //if j==k, do the equation, else, exact = 0.0. tgamma function is eqal to (n-1)! for positive integers
            double exact = (j == k) ? pow(2, j) * tgamma(j + 1) * sqrt(PI) : 0.0; 

            //looping for all n rules from 2-10
            for (int ni = 0; ni < 9; ni++) {
                int n = n_values[ni];
                
                //looping for all m values
                for (int mi = 0; mi < 4; mi++) {
                    int m = m_values[mi];
                    
                    //calls function to do numerical solution, and then finds difference between exact and numerical solution
                    double numerical = compute_integral(j, k, m, n, xmax);
                    double error = fabs(numerical - exact);

                    //changes from printing doubles to scientific notation at j=5
                    if (j<5){
                    printf("%2d |%2d |%3d  |%9d |%10.3f |%7.3f |%8.2e \n",
                           j, k, n, m, numerical, exact, error);
                    } else {
                    printf("%2d |%2d |%3d  |%9d |%10.3e |%7.3e |%8.2e \n",
                           j, k, n, m, numerical, exact, error);
                    }
                }
            }
        }
    }
    return 0;
}
