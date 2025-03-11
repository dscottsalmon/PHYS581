#include <stdio.h>

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

int main() {
    
    //arbitrary grid of 10 evenly spaced points from -1 to 1
    double x[10] = {-1.0, -0.7777, -0.5556, -0.3333, -0.1111, 0.1111, 0.3333, 0.5556, 0.7777, 1.0};
    double w[10];

    computeBarycentricWeights(10, x, w);

    // Print results
    printf("Barycentric Weights:\n");
    for (int i = 0; i < 10; i++) {
        printf("w[%d] = %f\n", i, w[i]);
    }

    return 0;
}
