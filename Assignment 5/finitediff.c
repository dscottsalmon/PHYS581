
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <direct.h>
#include <errno.h>

//Constants
#define L 1.0         // Box length
#define a 0.1         // Initial Gaussian width
#define c 1.0         // Speed of sound
#define Nx 100        // Number of spatial points
#define dx (L / (Nx - 1))   // Spatial step size (fixed)

//Initialize the wave with Gaussian profile and zero initial velocity
void initialize_wave(double *f, double *v, double *x) {
    for (int i = 0; i < Nx; i++) {
        x[i] = -L / 2 + i * dx;
        f[i] = exp(-x[i] * x[i] / (a * a));  // Gaussian initial condition
        v[i] = 0.0;
    }
}

//Compute second spatial derivative using fourth-order finite differencing
void compute_spatial_derivatives(double *f, double *d2f_dx2) {
    
    for (int i = 2; i < Nx - 2; i++) {
        d2f_dx2[i] = (-f[i - 2] + 16 * f[i - 1] - 30 * f[i] + 16 * f[i + 1] - f[i + 2]) / (12 * dx * dx);
    }

    // Enforce boundary conditions
    d2f_dx2[0] = d2f_dx2[1] = d2f_dx2[Nx - 2] = d2f_dx2[Nx - 1] = 0.0;
}

// Perform one RK45 step to update f and v
void rk45_step(double *f, double *v, double *d2f_dx2, double dt, double *f_next, double *v_next) {
    double f_stage[Nx], d2f_stage[Nx];
    double kf1[Nx], kv1[Nx], kf2[Nx], kv2[Nx], kf3[Nx], kv3[Nx];
    double kf4[Nx], kv4[Nx], kf5[Nx], kv5[Nx], kf6[Nx], kv6[Nx];

    for (int i = 0; i < Nx; i++) {
        kf1[i] = dt * v[i];
        kv1[i] = dt * c * c * d2f_dx2[i];
    }

    for (int i = 0; i < Nx; i++) f_stage[i] = f[i] + 0.25 * kf1[i];
    compute_spatial_derivatives(f_stage, d2f_stage);
    for (int i = 0; i < Nx; i++) {
        kf2[i] = dt * (v[i] + 0.25 * kv1[i]);
        kv2[i] = dt * c * c * d2f_stage[i];
    }

    for (int i = 0; i < Nx; i++)
        f_stage[i] = f[i] + (3.0/32.0)*kf1[i] + (9.0/32.0)*kf2[i];
    compute_spatial_derivatives(f_stage, d2f_stage);
    for (int i = 0; i < Nx; i++) {
        kf3[i] = dt * (v[i] + (3.0/32.0)*kv1[i] + (9.0/32.0)*kv2[i]);
        kv3[i] = dt * c * c * d2f_stage[i];
    }

    for (int i = 0; i < Nx; i++)
        f_stage[i] = f[i] + (1932.0/2197.0)*kf1[i] - (7200.0/2197.0)*kf2[i] + (7296.0/2197.0)*kf3[i];
    compute_spatial_derivatives(f_stage, d2f_stage);
    for (int i = 0; i < Nx; i++) {
        kf4[i] = dt * (v[i] + (1932.0/2197.0)*kv1[i] - (7200.0/2197.0)*kv2[i] + (7296.0/2197.0)*kv3[i]);
        kv4[i] = dt * c * c * d2f_stage[i];
    }

    for (int i = 0; i < Nx; i++)
        f_stage[i] = f[i] + (439.0/216.0)*kf1[i] - 8.0*kf2[i] + (3680.0/513.0)*kf3[i] - (845.0/4104.0)*kf4[i];
    compute_spatial_derivatives(f_stage, d2f_stage);
    for (int i = 0; i < Nx; i++) {
        kf5[i] = dt * (v[i] + (439.0/216.0)*kv1[i] - 8.0*kv2[i] + (3680.0/513.0)*kv3[i] - (845.0/4104.0)*kv4[i]);
        kv5[i] = dt * c * c * d2f_stage[i];
    }

    for (int i = 0; i < Nx; i++)
        f_stage[i] = f[i] - (8.0/27.0)*kf1[i] + 2.0*kf2[i] - (3544.0/2565.0)*kf3[i]
                     + (1859.0/4104.0)*kf4[i] - (11.0/40.0)*kf5[i];
    compute_spatial_derivatives(f_stage, d2f_stage);
    for (int i = 0; i < Nx; i++) {
        kf6[i] = dt * (v[i] - (8.0/27.0)*kv1[i] + 2.0*kv2[i] - (3544.0/2565.0)*kv3[i]
                        + (1859.0/4104.0)*kv4[i] - (11.0/40.0)*kv5[i]);
        kv6[i] = dt * c * c * d2f_stage[i];
    }

    for (int i = 0; i < Nx; i++) {
        f_next[i] = f[i] + (16.0/135.0)*kf1[i] + (6656.0/12825.0)*kf3[i]
                          + (28561.0/56430.0)*kf4[i] - (9.0/50.0)*kf5[i] + (2.0/55.0)*kf6[i];

        v_next[i] = v[i] + (16.0/135.0)*kv1[i] + (6656.0/12825.0)*kv3[i]
                          + (28561.0/56430.0)*kv4[i] - (9.0/50.0)*kv5[i] + (2.0/55.0)*kv6[i];
    }

    // Enforce Dirichlet boundary conditions
    f_next[0] = f_next[Nx - 1] = 0.0;
    v_next[0] = v_next[Nx - 1] = 0.0;
}

// Save snapshot to file
void save_snapshot(double *f, int timestep) {
    char filename[200];
    snprintf(filename, sizeof(filename), "Finite_Diff_Snapshots/wave_timestep_%d.txt", timestep);

    FILE* file = fopen(filename, "w");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < Nx; i++) {
        fprintf(file, "%f\n", f[i]);
    }

    fclose(file);
}

int main() {
    
    // Arrays for spatial grid and wave amplitudes
    double x[Nx];          // Spatial grid points
    double f[Nx];          // Current wave amplitude
    double v[Nx];          // Current velocity
    double f_next[Nx];     // Next wave amplitude
    double v_next[Nx];     // Next velocity
    double d2f_dx2[Nx];    // Second spatial derivative

    // Create the directory for saving snapshots (ignore if it already exists)
    if (mkdir("Finite_Diff_Snapshots") != 0 && errno != EEXIST) {
        perror("mkdir failed");
        exit(EXIT_FAILURE);
    }

    // Initialize the wave and spatial grid
    initialize_wave(f, v, x);

    // Time propagation loop
    double dt = 0.0005;     // Time step size
    int Nt = 4000;          // Number of time steps

    for (int n = 0; n < Nt; n++) {
        compute_spatial_derivatives(f, d2f_dx2);      // Compute spatial derivatives
        rk45_step(f, v, d2f_dx2, dt, f_next, v_next); // RK45 time integration

        if (n % 20 == 0) {
        save_snapshot(f_next, n);
    }

        // Update arrays for the next iteration
        for (int i = 0; i < Nx; i++) {
            f[i] = f_next[i];
            v[i] = v_next[i];
        }
    }
    return 0;
}
