
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>

#ifdef _WIN32
#include <direct.h>
#define mkdir_if_needed(dir) _mkdir(dir)
#else
#define mkdir_if_needed(dir) mkdir(dir, 0777)
#endif

#define L 1.0
#define a 0.1
#define c 1.0
#define Nx 64
#define PI 3.141592653589793
#define dx (L / (Nx + 1))

// Function prototypes
void initialize_wave(double *f, double *v);
void dst(double *in, double *out);
void idst(double *in, double *out);
void compute_d2f_dx2_sine(double *f, double *d2f_dx2);
void rk45_step(double *f, double *v, double *d2f_dx2, double dt, double *f_next, double *v_next);
void save_snapshot(double *f, int timestep);

int main() {
    double f[Nx], v[Nx], f_next[Nx], v_next[Nx], d2f_dx2[Nx];

    if (mkdir_if_needed("SineFourier_Snapshots") != 0 && errno != EEXIST) {
        perror("mkdir failed");
        exit(EXIT_FAILURE);
    }

    initialize_wave(f, v);

    double dt = 0.0005;
    int Nt = 4000;

    for (int n = 0; n < Nt; n++) {
        compute_d2f_dx2_sine(f, d2f_dx2);
        rk45_step(f, v, d2f_dx2, dt, f_next, v_next);

        if (n % 20 == 0)
            save_snapshot(f_next, n);

        for (int i = 0; i < Nx; i++) {
            f[i] = f_next[i];
            v[i] = v_next[i];
        }
    }

    return 0;
}

// Gaussian initial condition, 0 boundaries implied
void initialize_wave(double *f, double *v) {
    for (int i = 0; i < Nx; i++) {
        double x = -L/2 + (i + 1) * dx;
        f[i] = exp(-x * x / (a * a));
        v[i] = 0.0;
    }
}

// Naive DST-I
void dst(double *in, double *out) {
    for (int k = 0; k < Nx; k++) {
        out[k] = 0.0;
        for (int j = 0; j < Nx; j++) {
            out[k] += in[j] * sin(PI * (j + 1) * (k + 1) / (Nx + 1));
        }
    }
}

// Inverse DST-I (same as forward, scaled)
void idst(double *in, double *out) {
    for (int j = 0; j < Nx; j++) {
        out[j] = 0.0;
        for (int k = 0; k < Nx; k++) {
            out[j] += in[k] * sin(PI * (j + 1) * (k + 1) / (Nx + 1));
        }
        out[j] *= 2.0 / (Nx + 1);
    }
}

void compute_d2f_dx2_sine(double *f, double *d2f_dx2) {
    double f_hat[Nx], d2f_hat[Nx];
    dst(f, f_hat);

    for (int k = 0; k < Nx; k++) {
        double lambda_k = PI * (k + 1) / L;
        d2f_hat[k] = -lambda_k * lambda_k * f_hat[k];
    }

    idst(d2f_hat, d2f_dx2);
}

void rk45_step(double *f, double *v, double *d2f_dx2, double dt, double *f_next, double *v_next) {
    double f_stage[Nx], d2f_stage[Nx];
    double kf1[Nx], kv1[Nx], kf2[Nx], kv2[Nx], kf3[Nx], kv3[Nx];
    double kf4[Nx], kv4[Nx], kf5[Nx], kv5[Nx], kf6[Nx], kv6[Nx];

    for (int i = 0; i < Nx; i++) {
        kf1[i] = dt * v[i];
        kv1[i] = dt * c * c * d2f_dx2[i];
    }

    for (int i = 0; i < Nx; i++) f_stage[i] = f[i] + 0.25 * kf1[i];
    compute_d2f_dx2_sine(f_stage, d2f_stage);
    for (int i = 0; i < Nx; i++) {
        kf2[i] = dt * (v[i] + 0.25 * kv1[i]);
        kv2[i] = dt * c * c * d2f_stage[i];
    }

    for (int i = 0; i < Nx; i++)
        f_stage[i] = f[i] + (3.0/32.0)*kf1[i] + (9.0/32.0)*kf2[i];
    compute_d2f_dx2_sine(f_stage, d2f_stage);
    for (int i = 0; i < Nx; i++) {
        kf3[i] = dt * (v[i] + (3.0/32.0)*kv1[i] + (9.0/32.0)*kv2[i]);
        kv3[i] = dt * c * c * d2f_stage[i];
    }

    for (int i = 0; i < Nx; i++)
        f_stage[i] = f[i] + (1932.0/2197.0)*kf1[i] - (7200.0/2197.0)*kf2[i] + (7296.0/2197.0)*kf3[i];
    compute_d2f_dx2_sine(f_stage, d2f_stage);
    for (int i = 0; i < Nx; i++) {
        kf4[i] = dt * (v[i] + (1932.0/2197.0)*kv1[i] - (7200.0/2197.0)*kv2[i] + (7296.0/2197.0)*kv3[i]);
        kv4[i] = dt * c * c * d2f_stage[i];
    }

    for (int i = 0; i < Nx; i++)
        f_stage[i] = f[i] + (439.0/216.0)*kf1[i] - 8.0*kf2[i] + (3680.0/513.0)*kf3[i] - (845.0/4104.0)*kf4[i];
    compute_d2f_dx2_sine(f_stage, d2f_stage);
    for (int i = 0; i < Nx; i++) {
        kf5[i] = dt * (v[i] + (439.0/216.0)*kv1[i] - 8.0*kv2[i] + (3680.0/513.0)*kv3[i] - (845.0/4104.0)*kv4[i]);
        kv5[i] = dt * c * c * d2f_stage[i];
    }

    for (int i = 0; i < Nx; i++)
        f_stage[i] = f[i] - (8.0/27.0)*kf1[i] + 2.0*kf2[i] - (3544.0/2565.0)*kf3[i]
                     + (1859.0/4104.0)*kf4[i] - (11.0/40.0)*kf5[i];
    compute_d2f_dx2_sine(f_stage, d2f_stage);
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
}

void save_snapshot(double *f, int timestep) {
    char filename[200];
    snprintf(filename, sizeof(filename), "SineFourier_Snapshots/wave_timestep_%d.txt", timestep);
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Error writing snapshot");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < Nx; i++)
        fprintf(file, "%f\n", f[i]);
    fclose(file);
}
