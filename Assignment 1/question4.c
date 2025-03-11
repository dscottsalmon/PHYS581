
#include <stdio.h>
#include <math.h>

//Function to compute T_10(x)
double T10(double x) {
    
    if (x > 1.0) x = 1.0;
    if (x < -1.0) x = -1.0;
    
    return cos(10 * acos(x));
}

//Function to compute the first derivative of T_10(x) using finite difference
double T10_prime_fd(int i, int m, double x, double h) {
    
    //left boundary, uses forwards difference (since we cant do x-h)
    if (i==0) {
        return (T10(x+h) - T10(x)) / h;
        }
    
    //right boundary, uses backwards difference (since we cant do x+h)
    else if(i == m-1) {
        return (T10(x) - T10(x-h)) / h;
        }
    
    //else, central difference (because we can do x+h and x-h)
    else{return (T10(x+h) - T10(x-h)) / (2*h);}
}


// Function to compute the second derivative of T_10(x) using finite difference
double T10_double_prime_fd(int i, int m, double x, double h) {
    
    //left boundary (i == 0): use forward difference
    if (i==0) {
        return (T10(x + 2 * h) - 2 * T10(x + h) + T10(x)) / (h * h);
    }
    // Right boundary (i == m - 1): use backward difference
    else if (i == m - 1) {
        return (T10(x) - 2 * T10(x - h) + T10(x - 2 * h)) / (h * h);
    }
    // Central difference for the interior points
    else {
        return (T10(x + h) - 2 * T10(x) + T10(x - h)) / (h * h);
    }
}

int main() {
    
    //label for differentiation method
    printf("\nDerivatives to Chebyshev n=10 polynomial, using finite difference numerical methods.\n");
    printf("Primarily central finite difference (outside of boundaries, which used forwards/backwards difference).\n\n");
    
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
        T10_prime_values[i] = T10_prime_fd(i, m, x, a);
        T10_double_prime_values[i] = T10_double_prime_fd(i, m, x, a);
        
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
