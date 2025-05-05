#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M_PI 3.14159265358979323846
#define G 4 * M_PI * M_PI  // Constante gravitacional en UA^3 / (MS * año^2)
#define N 8                // Número de planetas (Mercurio a Neptuno)
#define SOL_INDEX 0

double dt = 0.001;         // Paso temporal en años
int steps = 10000;         // Número de pasos

// Arrays con el Sol + planetas (índice 0 = Sol, 1–8 = planetas)
double mass[N + 1];
double pos[N + 1][3];
double vel[N + 1][3];
double acc[N + 1][3];
double ecc[N + 1];  // Solo útil para planetas

void rescale_system() {
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j < 3; j++) {
            pos[i][j] /= 1.0;
            vel[i][j] /= 1.0;
        }
    }
}

void unscale_system() {
    // Nada que hacer en este ejemplo
}

void compute_accelerations() {
    for (int i = 1; i <= N; i++) {  // Solo los planetas
        double dx = pos[SOL_INDEX][0] - pos[i][0];
        double dy = pos[SOL_INDEX][1] - pos[i][1];
        double dz = pos[SOL_INDEX][2] - pos[i][2];
        double r3 = pow(dx*dx + dy*dy + dz*dz, 1.5);

        acc[i][0] = G * mass[SOL_INDEX] * dx / r3;
        acc[i][1] = G * mass[SOL_INDEX] * dy / r3;
        acc[i][2] = G * mass[SOL_INDEX] * dz / r3;
    }
}

void update_positions() {
    for (int i = 1; i <= N; i++) {
        for (int j = 0; j < 3; j++) {
            pos[i][j] += vel[i][j] * dt + 0.5 * acc[i][j] * dt * dt;
        }
    }
}

void update_velocities(double acc_new[N + 1][3]) {
    for (int i = 1; i <= N; i++) {
        for (int j = 0; j < 3; j++) {
            vel[i][j] += 0.5 * (acc[i][j] + acc_new[i][j]) * dt;
            acc[i][j] = acc_new[i][j];
        }
    }
}

void write_output(FILE* files[N + 1]) {
    for (int i = 1; i <= N; i++) {
        fprintf(files[i], "%lf %lf %lf\n", pos[i][0], pos[i][1], pos[i][2]);
    }
}

int main() {
    // Inicializa el Sol
    mass[SOL_INDEX] = 1.0;
    pos[SOL_INDEX][0] = pos[SOL_INDEX][1] = pos[SOL_INDEX][2] = 0.0;
    vel[SOL_INDEX][0] = vel[SOL_INDEX][1] = vel[SOL_INDEX][2] = 0.0;

    // Masas, posiciones y velocidades iniciales simplificadas
    double initial_mass[N] = {1.0e-6, 1.66e-7, 2.45e-6, 3.0e-6, 9.54e-4, 2.86e-4, 4.4e-5, 5.15e-5};
    double initial_pos[N][3] = {
        {0.39, 0.0, 0.0}, {0.72, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.52, 0.0, 0.0},
        {5.2, 0.0, 0.0}, {9.58, 0.0, 0.0}, {19.2, 0.0, 0.0}, {30.1, 0.0, 0.0}
    };
    double initial_vel[N][3] = {
        {0.0, 9.91, 0.0}, {0.0, 7.36, 0.0}, {0.0, 6.28, 0.0}, {0.0, 5.08, 0.0},
        {0.0, 2.63, 0.0}, {0.0, 2.03, 0.0}, {0.0, 1.43, 0.0}, {0.0, 1.14, 0.0}
    };
    double initial_ecc[N] = {0.206, 0.007, 0.017, 0.093, 0.049, 0.057, 0.046, 0.010};

    for (int i = 1; i <= N; i++) {
        mass[i] = initial_mass[i - 1];
        ecc[i] = initial_ecc[i - 1];
        for (int j = 0; j < 3; j++) {
            pos[i][j] = initial_pos[i - 1][j];
            vel[i][j] = initial_vel[i - 1][j];
        }
    }

    rescale_system();

    FILE* files[N + 1];
    char filename[20];
    for (int i = 1; i <= N; i++) {
        sprintf(filename, "planet_%d.txt", i);
        files[i] = fopen(filename, "w");
        if (!files[i]) {
            perror("Error abriendo fichero");
            exit(1);
        }
    }

    compute_accelerations();

    for (int step = 0; step < steps; step++) {
        update_positions();
        double acc_new[N + 1][3];
        compute_accelerations();
        for (int i = 1; i <= N; i++)
            for (int j = 0; j < 3; j++)
                acc_new[i][j] = acc[i][j];
        update_velocities(acc_new);
        write_output(files);
    }

    for (int i = 1; i <= N; i++) fclose(files[i]);

    unscale_system();

    return 0;
}
