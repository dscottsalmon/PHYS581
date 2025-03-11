#include <stdio.h>
#include <math.h>

//Function to compute T_10(x)
double T10(double x) {
    if (x > 1.0) x = 1.0;
    if (x < -1.0) x = -1.0;
    return cos(10 * acos(x));
}

//Function to compute the first derivative of T_10(x)
double T10_prime(double x) {
    return 10*sin(10*acos(x))/(sqrt(1-x*x));
}

//Function to compute the second derivative of T_10(x)
double T10_double_prime(double x) {
    return -100 * cos(10 * acos(x))/(1-x*x) + 10 * x * sin(10 * acos(x))/pow((1-x*x), (1.5));
}

int main() {
    
    //label for differentiation method
    printf("\nDerivatives to Chebyshev n=10 polynomial, using analytic solutions given by Wolfram Mathematica.\n");
    printf("I then applied the analytic solution with 101 different points (each seperated by an equal distance from -1 to 1) to the function.\n\n");
    
    //Number of points
    int m = 101;
    
    //distance between each point
    double a = 2.0/(m-1);
    
    //initial value for x
    double x = -1.0;
    
    //initializing 3 double arrays with size m for storing values
    double T10_values[m], T10_prime_values[m], T10_double_prime_values[m];
    
    //Loop through the equally spaced grid points between -1 and 1
    for (int i = 0; i < m; i++) {
        
        // Store the values in arrays
        T10_values[i] = T10(x);
        T10_prime_values[i] = T10_prime(x);
        T10_double_prime_values[i] = T10_double_prime(x);
        
        //iterate x with distance a
        x = x + a;
        
        //safety check to ensure -1.0<=x<=1.0 (else trig functions will be undefined)
        if (x > 1.0) x = 1.0;
        
    }

    //reinitializing x
    x = -1.0;
    // Print the results
    printf("x\t\tT10(x)\t\tT10'(x)\t\tT10''(x)\n");
    printf("-------------------------------------------------------------\n");
    for (int i = 0; i < m; i++) {
        printf("%f\t%f\t%f\t%f\n", x, T10_values[i], T10_prime_values[i], T10_double_prime_values[i]);
        x = x + a;
        if (x > 1.0) x = 1.0;
    }

    return 0;
}
