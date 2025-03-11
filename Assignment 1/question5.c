
#include <stdio.h>
#include <math.h>

//Function to compute T_10(x)
double T10(double x) {
    
    if (x > 1.0) x = 1.0;
    if (x < -1.0) x = -1.0;
    
    return cos(10 * acos(x));
}

//Function to compute the first derivative of T_10(x), analytically
double T10_prime(double x) {
    return 10*sin(10*acos(x))/(sqrt(1-x*x));
}

//Function to compute the second derivative of T_10(x), analytically
double T10_double_prime(double x) {
    return -100 * cos(10 * acos(x))/(1-x*x) + 10 * x * sin(10 * acos(x))/pow((1-x*x), (1.5));
}

//Function to compute the first derivative of T_10(x) using finite difference method
double T10_prime_fd(int i, int m, double x, double h) {
    
    //left boundary, uses forwards difference (since we cant do x-h)
    if (i==0) {
        return (T10(x+h) - T10(x)) / h;
    }
    
    //right boundary, uses backwards difference (since we cant do x+h)
    else if (i == m-1) {
        return (T10(x) - T10(x-h)) / h;
    }
    
    //else, central difference (because we can do x+h and x-h)
    else{
        return (T10(x+h) - T10(x-h)) / (2*h);
    }
}


// Function to compute the second derivative of T_10(x) using finite difference method
double T10_double_prime_fd(int i, int m, double x, double h) {
    
    //left boundary, uses forwards difference (since we cant do x-h)
    if (i == 0) {
        return (T10(x + 2 * h) - 2 * T10(x + h) + T10(x)) / (h * h);
    }
    
    //right boundary, uses backwards difference (since we cant do x+h)
    else if (i == m - 1) { 
        return (T10(x) - 2 * T10(x - h) + T10(x - 2 * h)) / (h * h);
    }
    
    //else, central difference (because we can do x+h and x-h)
    else {
        return (T10(x + h) - 2 * T10(x) + T10(x - h)) / (h * h);
    }
}

int main() {
    
    printf("\nComparing finite difference methods with analytical solutions.\n");
    printf("Error is the absolute value of the average difference between the two solutions.\n\n");
    
    
    //Number of points to test
    int pointValues[] = {10, 15, 20, 25, 30, 35, 40, 45, 50, 75, 100, 125, 150, 175, 200, 250, 275, 300, 350, 500, 999};
    int length = sizeof(pointValues)/sizeof(pointValues[0]);
    
    //header print statement for table
    printf("# of Points    Point Spacing (a)    Average Error T10'(x)   Average Error T10''(x)\n");
    printf("----------------------------------------------------------------------------------\n");
    
    
    //main loop that will iterate through each of the point values
    for (int j=0; j<length; j++){
        
        //assigning m to the number of points we're testing this iteration
        int m = pointValues[j];
        
        //distance between each point
        double a = 2.0/(m-1);
        
        //initial value for x
        double x = -1.0;
        
        //initializing 3 double arrays with size m for storing analytic values
        double T10_values_analytic[m], T10_prime_values_analytic[m], T10_double_prime_values_analytic[m];
        
        //initializing 3 double arrays with size m for storing finite difference values
        double T10_values[m], T10_prime_values[m], T10_double_prime_values[m];
        
        //initializing average tracking variable
        double averageprime = 0;
        double averagedoubleprime = 0;
        
        //Loop through the equally spaced grid points between -1 and 1
        for (int i = 0; i < m; i++) {
            
            //Store the analytic values in their arrays
            T10_prime_values_analytic[i] = T10_prime(x);
            T10_double_prime_values_analytic[i] = T10_double_prime(x);
            
            //Store the fd values in their arrays
            T10_prime_values[i] = T10_prime_fd(i, m, x, a);
            T10_double_prime_values[i] = T10_double_prime_fd(i, m, x, a);
            
            //finds differences between fd and analytic solutions
            double primedifference = fabs(T10_prime_values_analytic[i]-T10_prime_values[i]);
            double doubleprimedifference = fabs(T10_double_prime_values_analytic[i]-T10_double_prime_values[i]);
            
            //if condition to ignore boundaries (nan for analytic solution), adding to sum
            if (i!= 0 && i!= m-1) {
                averageprime = averageprime + primedifference;
                averagedoubleprime = averagedoubleprime + doubleprimedifference;
            }
            
            //iterate x with distance a
            x = x + a;
            
            //safety check to ensure -1.0<=x<=1.0 (else trig functions will be undefined)
            if (x > 1.0) x = 1.0;
        }
        
        //divides sums by total number of points to find average
        averageprime = averageprime/m;
        averagedoubleprime = averagedoubleprime/m;
        
        //prints results in table
        printf("    %i\t\t   %f\t\t%f\t\t%f\n", m, a, averageprime, averagedoubleprime);
        }
    return 0;
}
