#include <GL/glut.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cmath>

// Структура для описания частицы (тела)
struct particle {
    double x, y, z;
};

int nBodies = 9;

struct particle *p = NULL;       // Положения тел
struct particle *v = NULL;       // Скорости тел
double *m = NULL;                // Гравитационные параметры GM (в а.е.^3/дней^2)
double *body_radius = NULL;      // Радиусы для отрисовки
double (*body_color)[3] = NULL;   // Цвета тел

const double dt = 0.1;
const double DRAW_SCALE = 300.0;

// Скорость света в а.е./день и её квадрат:
const double c = 173.144632674240;
const double c2 = c * c;

int winWidth = 1920, winHeight = 1080;

struct particle g_barycenter = {0,0,0};
struct particle g_momentum = {0,0,0};
struct particle g_angularMomentum = {0,0,0};
double g_totalEnergy = 0.0;

double camX = 0.0, camY = 0.0, camZ = 2000.0;

std::vector<double> time_points;
std::vector<double> cm_magnitude;
std::vector<double> energy_rel;
std::vector<double> momentum_rel;
std::vector<double> angular_momentum_rel;

double initial_energy;
static bool firstCall = true;
particle initial_momentum;
particle initial_angularMomentum;
double simulation_time = 0.0;

//DOPRI-8
const int s = 12;
static const double c_dopri[s] = {
    0.0,
    1.0/18.0,
    1.0/12.0,
    1.0/8.0,
    5.0/16.0,
    3.0/8.0,
    59.0/400.0,
    93.0/200.0,
    5490023248.0/9719169821.0,
    13.0/20.0,
    1.0,
    1.0/2.0
};
static const double A[s][s] = {
    { 0.0 },
    { 1.0/18.0,                0.0 },
    { 1.0/48.0, 1.0/16.0,      0.0 },
    { 1.0/32.0, 0.0, 3.0/32.0, 0.0 },
    { 5.0/16.0, 0.0, -75.0/64.0, 75.0/64.0, 0.0 },
    { 3.0/80.0, 0.0, 0.0, 3.0/16.0, 3.0/20.0, 0.0 },
    { 29443841.0/614563906.0, 0.0, 0.0,
      77736538.0/692538347.0, -28693883.0/1125000000.0, 23124283.0/1800000000.0, 0.0 },
    { 16016141.0/946692911.0, 0.0, 0.0,
      61564180.0/158732637.0, 22789713.0/633445777.0, 545815736.0/2771057229.0,
      -180193667.0/1043307555.0, 0.0 },
    { 39632708.0/573591083.0, 0.0, 0.0,
      -433636366.0/683701615.0, -421739975.0/2616292301.0,
      100302831.0/723423059.0, 790204164.0/839813087.0,
      800635310.0/3783071287.0, 0.0 },
    { 246121993.0/1340847787.0, 0.0, 0.0,
      -37695042795.0/15268766246.0, -309121744.0/1061227803.0,
      -12992083.0/490766935.0, 6005943493.0/2108947869.0,
      393006217.0/1396673457.0, 123872331.0/1001029789.0, 0.0 },
    { 1028468189.0/846180014.0, 0.0, 0.0,
      8478235783.0/508512852.0, 1311729495.0/1432422823.0,
      -10304129995.0/1701304382.0, -48777925059.0/3047939560.0,
      15336726248.0/1032824649.0, -45442868181.0/3398467696.0,
      3065993473.0/597172653.0, 0.0 },
    { 185892177.0/718116043.0, 0.0, 0.0,
      -3185094517.0/667107341.0, -477755414.0/1098053517.0,
      -703635378.0/230739211.0, 5731566787.0/1027545527.0,
      5232866602.0/850066563.0, -4093664535.0/808688257.0,
      3962137247.0/1805957418.0, 65686358.0/487910083.0, 0.0 }
};
static const double b_main[s] = {
    14005451.0/335480064.0,
    0.0,
    0.0,
    0.0,
    0.0,
    -59238493.0/1068277825.0,
    181606767.0/758867731.0,
    561292985.0/797845732.0,
    -1041891430.0/1371343529.0,
    760417239.0/1151165299.0,
    118820643.0/751138087.0,
    -528747749.0/2220607170.0
};
static const double b_err[s] = {
    13451932.0/455176623.0,
    0.0,
    0.0,
    0.0,
    0.0,
    -808719846.0/976000145.0,
    1757004468.0/5645159321.0,
    656045339.0/265891186.0,
    -3867574721.0/1518517206.0,
    465885868.0/322736535.0,
    53011238.0/667516719.0,
    2.0/45.0
};


void compute_accelerations(struct particle *p_array, struct particle *v_array, struct particle *a_array) {
    for (int a = 0; a < nBodies; a++) {
        a_array[a].x = 0.0;
        a_array[a].y = 0.0;
        a_array[a].z = 0.0;
        for (int b = 0; b < nBodies; b++) {
            if (b == a) continue;
            double dx = p_array[a].x - p_array[b].x;
            double dy = p_array[a].y - p_array[b].y;
            double dz = p_array[a].z - p_array[b].z;
            double r_sq = dx*dx + dy*dy + dz*dz;
            double r = sqrt(r_sq);
            if (r == 0.0) continue;
            double r3 = r_sq * r;
            // Ньютоновский вклад
            a_array[a].x += - (m[b] / r3) * dx;
            a_array[a].y += - (m[b] / r3) * dy;
            a_array[a].z += - (m[b] / r3) * dz;
        }
    }
    double sum_a[ nBodies ];
    double v_sq[ nBodies ];
    for (int a = 0; a < nBodies; a++) {
        sum_a[a] = 0.0;
        v_sq[a] = v_array[a].x*v_array[a].x + v_array[a].y*v_array[a].y + v_array[a].z*v_array[a].z;
        for (int c = 0; c < nBodies; c++) {
            if (c == a) continue;
            double dx = p_array[a].x - p_array[c].x;
            double dy = p_array[a].y - p_array[c].y;
            double dz = p_array[a].z - p_array[c].z;
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            if (r == 0.0) continue;
            sum_a[a] += m[c] / r;
        }
    }
    // Релятивистские поправки
    for (int a = 0; a < nBodies; a++) {
        for (int b = 0; b < nBodies; b++) {
            if (b == a) continue;
            double dx = p_array[a].x - p_array[b].x;
            double dy = p_array[a].y - p_array[b].y;
            double dz = p_array[a].z - p_array[b].z;
            double r_sq = dx*dx + dy*dy + dz*dz;
            double r = sqrt(r_sq);
            if (r == 0.0) continue;
            double r3 = r_sq * r;

            double v_b_sq = v_array[b].x*v_array[b].x + v_array[b].y*v_array[b].y + v_array[b].z*v_array[b].z;
            double dot_v = v_array[a].x*v_array[b].x + v_array[a].y*v_array[b].y + v_array[a].z*v_array[b].z;
            double dot_x_vb = dx*v_array[b].x + dy*v_array[b].y + dz*v_array[b].z;
            double sum_b = 0.0;
            for (int c = 0; c < nBodies; c++) {
                if (c == b) continue;
                double dx_b = p_array[b].x - p_array[c].x;
                double dy_b = p_array[b].y - p_array[c].y;
                double dz_b = p_array[b].z - p_array[c].z;
                double r_bc = sqrt(dx_b*dx_b + dy_b*dy_b + dz_b*dz_b);
                if (r_bc == 0.0) continue;
                sum_b += m[c] / r_bc;
            }
            double bracket = 1.0
                - (4.0/c2)*sum_a[a]
                - (1.0/c2)*sum_b
                + (v_sq[a] / c2)
                + (2.0 * v_b_sq / c2)
                - (4.0 * dot_v / c2)
                - (3.0/(2.0*c2)) * ((dot_x_vb * dot_x_vb) / r_sq);
            double corr = bracket - 1.0;
            double T1_corr = - (m[b] / r3) * corr;

            double dot_term = dx*(3*v_array[b].x - 4*v_array[a].x)
                            + dy*(3*v_array[b].y - 4*v_array[a].y)
                            + dz*(3*v_array[b].z - 4*v_array[a].z);
            double T2 = (m[b] / r3) * (dot_term / c2);

            a_array[a].x += T1_corr * dx + T2 * (v_array[a].x - v_array[b].x);
            a_array[a].y += T1_corr * dy + T2 * (v_array[a].y - v_array[b].y);
            a_array[a].z += T1_corr * dz + T2 * (v_array[a].z - v_array[b].z);
        }
    }
    for (int a = 0; a < nBodies; a++) {
        double T3_x = 0.0, T3_y = 0.0, T3_z = 0.0;
        for (int b = 0; b < nBodies; b++) {
            if (b == a) continue;
            double dx_ab = p_array[a].x - p_array[b].x;
            double dy_ab = p_array[a].y - p_array[b].y;
            double dz_ab = p_array[a].z - p_array[b].z;
            double r_ab = sqrt(dx_ab*dx_ab + dy_ab*dy_ab + dz_ab*dz_ab);
            if (r_ab == 0.0) continue;
            double factor_b = m[b] / r_ab;
            for (int c = 0; c < nBodies; c++) {
                if (c == a || c == b) continue;
                double dx_ac = p_array[a].x - p_array[c].x;
                double dy_ac = p_array[a].y - p_array[c].y;
                double dz_ac = p_array[a].z - p_array[c].z;
                double r_ac_sq = dx_ac*dx_ac + dy_ac*dy_ac + dz_ac*dz_ac;
                double r_ac = sqrt(r_ac_sq);
                if (r_ac == 0.0) continue;
                double factor_c = m[c] / (r_ac_sq * r_ac);
                T3_x += factor_b * factor_c * dx_ac;
                T3_y += factor_b * factor_c * dy_ac;
                T3_z += factor_b * factor_c * dz_ac;
            }
        }
        double T3_coeff = 7.0 / (2.0 * c2);
        a_array[a].x += T3_coeff * T3_x;
        a_array[a].y += T3_coeff * T3_y;
        a_array[a].z += T3_coeff * T3_z;
    }
}


void rk4_step() {
    static struct particle *buffer = NULL;
    static struct particle *k1_p, *k1_v, *k2_p, *k2_v, *k3_p, *k3_v, *k4_p, *k4_v;
    static struct particle *temp_p, *temp_v, *a;
    static size_t prev_nBodies = 0;
    if (buffer == NULL || prev_nBodies != nBodies) {
        free(buffer);
        size_t total = 11 * nBodies;
        buffer = (struct particle*)calloc(total, sizeof(struct particle));
        if (!buffer) {
            exit(EXIT_FAILURE);
        }
        k1_p = buffer;
        k1_v = k1_p + nBodies;
        k2_p = k1_v + nBodies;
        k2_v = k2_p + nBodies;
        k3_p = k2_v + nBodies;
        k3_v = k3_p + nBodies;
        k4_p = k3_v + nBodies;
        k4_v = k4_p + nBodies;
        temp_p = k4_v + nBodies;
        temp_v = temp_p + nBodies;
        a = temp_v + nBodies;

        prev_nBodies = nBodies;
    }
    compute_accelerations(p, v, a);
    for (int i = 0; i < nBodies; i++) {
        k1_p[i] = v[i];
        k1_v[i] = a[i];
    }
    for (int i = 0; i < nBodies; i++) {
        temp_p[i].x = p[i].x + 0.5 * dt * k1_p[i].x;
        temp_p[i].y = p[i].y + 0.5 * dt * k1_p[i].y;
        temp_p[i].z = p[i].z + 0.5 * dt * k1_p[i].z;
        temp_v[i].x = v[i].x + 0.5 * dt * k1_v[i].x;
        temp_v[i].y = v[i].y + 0.5 * dt * k1_v[i].y;
        temp_v[i].z = v[i].z + 0.5 * dt * k1_v[i].z;
    }
    compute_accelerations(temp_p, temp_v, a);
    for (int i = 0; i < nBodies; i++) {
        k2_p[i] = temp_v[i];
        k2_v[i] = a[i];
    }
    for (int i = 0; i < nBodies; i++) {
        temp_p[i].x = p[i].x + 0.5 * dt * k2_p[i].x;
        temp_p[i].y = p[i].y + 0.5 * dt * k2_p[i].y;
        temp_p[i].z = p[i].z + 0.5 * dt * k2_p[i].z;
        temp_v[i].x = v[i].x + 0.5 * dt * k2_v[i].x;
        temp_v[i].y = v[i].y + 0.5 * dt * k2_v[i].y;
        temp_v[i].z = v[i].z + 0.5 * dt * k2_v[i].z;
    }
    compute_accelerations(temp_p, temp_v, a);
    for (int i = 0; i < nBodies; i++) {
        k3_p[i] = temp_v[i];
        k3_v[i] = a[i];
    }
    for (int i = 0; i < nBodies; i++) {
        temp_p[i].x = p[i].x + dt * k3_p[i].x;
        temp_p[i].y = p[i].y + dt * k3_p[i].y;
        temp_p[i].z = p[i].z + dt * k3_p[i].z;
        temp_v[i].x = v[i].x + dt * k3_v[i].x;
        temp_v[i].y = v[i].y + dt * k3_v[i].y;
        temp_v[i].z = v[i].z + dt * k3_v[i].z;
    }
    compute_accelerations(temp_p, temp_v, a);
    for (int i = 0; i < nBodies; i++) {
        k4_p[i] = temp_v[i];
        k4_v[i] = a[i];
    }
    for (int i = 0; i < nBodies; i++) {
        p[i].x += (dt / 6.0) * (k1_p[i].x + 2.0 * k2_p[i].x + 2.0 * k3_p[i].x + k4_p[i].x);
        p[i].y += (dt / 6.0) * (k1_p[i].y + 2.0 * k2_p[i].y + 2.0 * k3_p[i].y + k4_p[i].y);
        p[i].z += (dt / 6.0) * (k1_p[i].z + 2.0 * k2_p[i].z + 2.0 * k3_p[i].z + k4_p[i].z);
        
        v[i].x += (dt / 6.0) * (k1_v[i].x + 2.0 * k2_v[i].x + 2.0 * k3_v[i].x + k4_v[i].x);
        v[i].y += (dt / 6.0) * (k1_v[i].y + 2.0 * k2_v[i].y + 2.0 * k3_v[i].y + k4_v[i].y);
        v[i].z += (dt / 6.0) * (k1_v[i].z + 2.0 * k2_v[i].z + 2.0 * k3_v[i].z + k4_v[i].z);
    }
    simulation_time += dt;
}


void predictor_corrector_ab8_step() {
    const int iterations = 10;

    static struct particle *p_pred = NULL;
    static struct particle *v_pred = NULL;
    static struct particle *a_old = NULL;
    static struct particle *a_new = NULL;
    
    static struct particle *a_hist[7] = {NULL};
    static struct particle *v_hist[7] = {NULL};
    static size_t prev_nBodies = 0;
    static bool initialized = false;

    if (!initialized || prev_nBodies != nBodies) {
        free(p_pred); free(v_pred); free(a_old); free(a_new);
        
        for (int i = 0; i < 7; i++) {
            free(a_hist[i]);
            free(v_hist[i]);
            a_hist[i] = (struct particle*)calloc(nBodies, sizeof(struct particle));
            v_hist[i] = (struct particle*)calloc(nBodies, sizeof(struct particle));
        }

        p_pred = (struct particle*)calloc(nBodies, sizeof(struct particle));
        v_pred = (struct particle*)calloc(nBodies, sizeof(struct particle));
        a_old = (struct particle*)calloc(nBodies, sizeof(struct particle));
        a_new = (struct particle*)calloc(nBodies, sizeof(struct particle));

        compute_accelerations(p, v, a_old);
        for (int i = 0; i < 7; i++) {
            memcpy(a_hist[i], a_old, nBodies * sizeof(struct particle));
            memcpy(v_hist[i], v, nBodies * sizeof(struct particle));
        }

        prev_nBodies = nBodies;
        initialized = true;
    }

    compute_accelerations(p, v, a_old);
    
    double coeffs[8] = { 4277.0/1440.0, -7923.0/1440.0, 9982.0/1440.0, -7298.0/1440.0, 2877.0/1440.0, -475.0/1440.0, 27.0/1440.0, -1.0/1440.0 };

    for (int i = 0; i < nBodies; i++) {
        v_pred[i] = v[i];
        p_pred[i] = p[i];
        for (int j = 0; j < 8; j++) {
            v_pred[i].x += dt * coeffs[j] * (j < 7 ? a_hist[j][i].x : a_old[i].x);
            v_pred[i].y += dt * coeffs[j] * (j < 7 ? a_hist[j][i].y : a_old[i].y);
            v_pred[i].z += dt * coeffs[j] * (j < 7 ? a_hist[j][i].z : a_old[i].z);
            p_pred[i].x += dt * coeffs[j] * (j < 7 ? v_hist[j][i].x : v[i].x);
            p_pred[i].y += dt * coeffs[j] * (j < 7 ? v_hist[j][i].y : v[i].y);
            p_pred[i].z += dt * coeffs[j] * (j < 7 ? v_hist[j][i].z : v[i].z);
        }
    }

    for (int iter = 0; iter < iterations; iter++) {
        compute_accelerations(p_pred, v_pred, a_new);
        for (int i = 0; i < nBodies; i++) {
            v_pred[i].x = v[i].x + (dt / 2.0) * (a_old[i].x + a_new[i].x);
            v_pred[i].y = v[i].y + (dt / 2.0) * (a_old[i].y + a_new[i].y);
            v_pred[i].z = v[i].z + (dt / 2.0) * (a_old[i].z + a_new[i].z);
            p_pred[i].x = p[i].x + (dt / 2.0) * (v[i].x + v_pred[i].x);
            p_pred[i].y = p[i].y + (dt / 2.0) * (v[i].y + v_pred[i].y);
            p_pred[i].z = p[i].z + (dt / 2.0) * (v[i].z + v_pred[i].z);
        }
    }

    for (int i = 0; i < nBodies; i++) {
        p[i] = p_pred[i];
        v[i] = v_pred[i];
    }
    simulation_time += dt;

    for (int i = 6; i > 0; i--) {
        memmove(a_hist[i], a_hist[i - 1], nBodies * sizeof(struct particle));
        memmove(v_hist[i], v_hist[i - 1], nBodies * sizeof(struct particle));
    }
    memcpy(a_hist[0], a_old, nBodies * sizeof(struct particle));
    memcpy(v_hist[0], v, nBodies * sizeof(struct particle));
}


void predictor_corrector_am8_step() {
    const int iterations = 3;

    static struct particle *p_pred = NULL;
    static struct particle *v_pred = NULL;
    static struct particle *a_old = NULL;
    static struct particle *a_new = NULL;
    
    static struct particle *a_hist[7] = {NULL};
    static struct particle *v_hist[7] = {NULL};
    static size_t prev_nBodies = 0;
    static bool initialized = false;

    if (!initialized || prev_nBodies != nBodies) {
        free(p_pred); free(v_pred); free(a_old); free(a_new);
        for (int i = 0; i < 7; i++) {
            free(a_hist[i]);
            free(v_hist[i]);
            a_hist[i] = (struct particle*)calloc(nBodies, sizeof(struct particle));
            v_hist[i] = (struct particle*)calloc(nBodies, sizeof(struct particle));
        }

        p_pred = (struct particle*)calloc(nBodies, sizeof(struct particle));
        v_pred = (struct particle*)calloc(nBodies, sizeof(struct particle));
        a_old = (struct particle*)calloc(nBodies, sizeof(struct particle));
        a_new = (struct particle*)calloc(nBodies, sizeof(struct particle));

        compute_accelerations(p, v, a_old);
        for (int i = 0; i < 7; i++) {
            memcpy(a_hist[i], a_old, nBodies * sizeof(struct particle));
            memcpy(v_hist[i], v, nBodies * sizeof(struct particle));
        }

        prev_nBodies = nBodies;
        initialized = true;
    }

    compute_accelerations(p, v, a_old);
    
    double am_coeffs[8] = { 251.0/720.0, 646.0/720.0, -264.0/720.0, 106.0/720.0, -19.0/720.0, 0, 0, 0 };
    
    for (int i = 0; i < nBodies; i++) {
        v_pred[i] = v[i];
        p_pred[i] = p[i];
        for (int j = 0; j < 8; j++) {
            if (j < 7) {
                v_pred[i].x += dt * am_coeffs[j] * a_hist[j][i].x;
                v_pred[i].y += dt * am_coeffs[j] * a_hist[j][i].y;
                v_pred[i].z += dt * am_coeffs[j] * a_hist[j][i].z;
                p_pred[i].x += dt * am_coeffs[j] * v_hist[j][i].x;
                p_pred[i].y += dt * am_coeffs[j] * v_hist[j][i].y;
                p_pred[i].z += dt * am_coeffs[j] * v_hist[j][i].z;
            } else {
                v_pred[i].x += dt * am_coeffs[j] * a_old[i].x;
                v_pred[i].y += dt * am_coeffs[j] * a_old[i].y;
                v_pred[i].z += dt * am_coeffs[j] * a_old[i].z;
                p_pred[i].x += dt * am_coeffs[j] * v[i].x;
                p_pred[i].y += dt * am_coeffs[j] * v[i].y;
                p_pred[i].z += dt * am_coeffs[j] * v[i].z;
            }
        }
    }
    

    for (int iter = 0; iter < iterations; iter++) {
        compute_accelerations(p_pred, v_pred, a_new);
        for (int i = 0; i < nBodies; i++) {
            v_pred[i].x = v[i].x + (dt / 2.0) * (a_old[i].x + a_new[i].x);
            v_pred[i].y = v[i].y + (dt / 2.0) * (a_old[i].y + a_new[i].y);
            v_pred[i].z = v[i].z + (dt / 2.0) * (a_old[i].z + a_new[i].z);
            p_pred[i].x = p[i].x + (dt / 2.0) * (v[i].x + v_pred[i].x);
            p_pred[i].y = p[i].y + (dt / 2.0) * (v[i].y + v_pred[i].y);
            p_pred[i].z = p[i].z + (dt / 2.0) * (v[i].z + v_pred[i].z);
        }
    }

    for (int i = 0; i < nBodies; i++) {
        p[i] = p_pred[i];
        v[i] = v_pred[i];
    }
    simulation_time += dt;

    for (int i = 6; i > 0; i--) {
        memmove(a_hist[i], a_hist[i - 1], nBodies * sizeof(struct particle));
        memmove(v_hist[i], v_hist[i - 1], nBodies * sizeof(struct particle));
    }
    memcpy(a_hist[0], a_old, nBodies * sizeof(struct particle));
    memcpy(v_hist[0], v, nBodies * sizeof(struct particle));
}


void dopri8_step() {
    static particle *k_p = nullptr;
    static particle *k_v = nullptr;
    static size_t prev_nBodies = 0;
    if (k_p == nullptr || prev_nBodies != (size_t)nBodies) {
        free(k_p); free(k_v);
        k_p = (particle*)calloc(s * nBodies, sizeof(particle));
        k_v = (particle*)calloc(s * nBodies, sizeof(particle));
        prev_nBodies = nBodies;
    }
    std::vector<particle> a_old(nBodies);
    compute_accelerations(p, v, a_old.data());
    for (int i = 0; i < nBodies; i++) {
        k_p[0 * nBodies + i] = v[i];
        k_v[0 * nBodies + i] = a_old[i];
    }
    std::vector<particle> p_temp(nBodies), v_temp(nBodies), a_stage(nBodies);
    for (int stage = 1; stage < s; stage++) {
        for (int i = 0; i < nBodies; i++) {
            p_temp[i] = p[i];
            v_temp[i] = v[i];
            for (int j = 0; j < stage; j++) {
                p_temp[i].x += dt * A[stage][j] * k_p[j * nBodies + i].x;
                p_temp[i].y += dt * A[stage][j] * k_p[j * nBodies + i].y;
                p_temp[i].z += dt * A[stage][j] * k_p[j * nBodies + i].z;
                
                v_temp[i].x += dt * A[stage][j] * k_v[j * nBodies + i].x;
                v_temp[i].y += dt * A[stage][j] * k_v[j * nBodies + i].y;
                v_temp[i].z += dt * A[stage][j] * k_v[j * nBodies + i].z;
            }
        }
        compute_accelerations(p_temp.data(), v_temp.data(), a_stage.data());
        for (int i = 0; i < nBodies; i++) {
            k_p[stage * nBodies + i] = v_temp[i];
            k_v[stage * nBodies + i] = a_stage[i];
        }
    }
    std::vector<particle> p_pred(nBodies), v_pred(nBodies);
    for (int i = 0; i < nBodies; i++) {
        p_pred[i] = p[i];
        v_pred[i] = v[i];
        for (int j = 0; j < s; j++) {
            p_pred[i].x += dt * b_main[j] * k_p[j * nBodies + i].x;
            p_pred[i].y += dt * b_main[j] * k_p[j * nBodies + i].y;
            p_pred[i].z += dt * b_main[j] * k_p[j * nBodies + i].z;
            
            v_pred[i].x += dt * b_main[j] * k_v[j * nBodies + i].x;
            v_pred[i].y += dt * b_main[j] * k_v[j * nBodies + i].y;
            v_pred[i].z += dt * b_main[j] * k_v[j * nBodies + i].z;
        }
    }
    std::vector<particle> p_err(nBodies), v_err(nBodies);
    for (int i = 0; i < nBodies; i++) {
        p_err[i].x = 0.0; p_err[i].y = 0.0; p_err[i].z = 0.0;
        v_err[i].x = 0.0; v_err[i].y = 0.0; v_err[i].z = 0.0;
        for (int j = 0; j < s; j++) {
            p_err[i].x += dt * b_err[j] * k_p[j * nBodies + i].x;
            p_err[i].y += dt * b_err[j] * k_p[j * nBodies + i].y;
            p_err[i].z += dt * b_err[j] * k_p[j * nBodies + i].z;
            
            v_err[i].x += dt * b_err[j] * k_v[j * nBodies + i].x;
            v_err[i].y += dt * b_err[j] * k_v[j * nBodies + i].y;
            v_err[i].z += dt * b_err[j] * k_v[j * nBodies + i].z;
        }
    }
    for (int i = 0; i < nBodies; i++) {
        p[i] = p_pred[i];
        v[i] = v_pred[i];
    }
    simulation_time += dt;
}


void calculate_energy_momentum() {
    double E_kin_lab = 0.0, E_pot_lab = 0.0;
    struct particle L_lab = {0, 0, 0};
    struct particle P_lab = {0, 0, 0};
    for (int i = 0; i < nBodies; i++) {
        double v2 = v[i].x*v[i].x + v[i].y*v[i].y + v[i].z*v[i].z;
        E_kin_lab += 0.5 * m[i] * v2;
        double px = m[i] * v[i].x;
        double py = m[i] * v[i].y;
        double pz = m[i] * v[i].z;
        L_lab.x += p[i].y * pz - p[i].z * py;
        L_lab.y += p[i].z * px - p[i].x * pz;
        L_lab.z += p[i].x * py - p[i].y * px;
        P_lab.x += m[i] * v[i].x;
        P_lab.y += m[i] * v[i].y;
        P_lab.z += m[i] * v[i].z;
    }
    for (int i = 0; i < nBodies; i++) {
        for (int j = i + 1; j < nBodies; j++) {
            double dx = p[j].x - p[i].x;
            double dy = p[j].y - p[i].y;
            double dz = p[j].z - p[i].z;
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            E_pot_lab += - (m[i] * m[j]) / r;
        }
    }
    double E_lab = E_kin_lab + E_pot_lab;
    struct particle cm = {0, 0, 0};
    double total_mass = 0.0;
    for (int i = 0; i < nBodies; i++) {
        cm.x += m[i] * p[i].x;
        cm.y += m[i] * p[i].y;
        cm.z += m[i] * p[i].z;
        total_mass += m[i];
    }
    cm.x /= total_mass;
    cm.y /= total_mass;
    cm.z /= total_mass;
    double E_new = E_lab;
    double P_new_mag = sqrt(P_lab.x*P_lab.x + P_lab.y*P_lab.y + P_lab.z*P_lab.z);
    double L_new_mag = sqrt(L_lab.x*L_lab.x + L_lab.y*L_lab.y + L_lab.z*L_lab.z);
    
    g_totalEnergy = E_new;
    g_momentum = P_lab;
    g_angularMomentum = L_lab;
}

void initBodies() {
    p = (struct particle*)malloc(nBodies * sizeof(struct particle));
    v = (struct particle*)malloc(nBodies * sizeof(struct particle));
    m = (double*)malloc(nBodies * sizeof(double));
    body_radius = (double*)malloc(nBodies * sizeof(double));
    body_color = (double (*)[3])malloc(nBodies * 3 * sizeof(double));

    double radii[9] = {
        25.0,  // Солнце
         2.0,  // Меркурий
         4.0,  // Венера
         4.5,  // Земля
         3.0,  // Марс
        10.0,  // Юпитер
         9.0,  // Сатурн
         5.0,  // Уран
         5.0   // Нептун
    };

    double colors[9][3] = {
        {1.0, 1.0, 0.0},
        {0.7, 0.7, 0.7},
        {1.0, 0.8, 0.6},
        {0.0, 0.0, 1.0},
        {1.0, 0.0, 0.0},
        {1.0, 0.5, 0.0},
        {0.9, 0.8, 0.5},
        {0.5, 1.0, 1.0},
        {0.0, 0.0, 0.5}
    };

    // Инициализация Солнца
    p[0].x = 0.0; p[0].y = 0.0; p[0].z = 0.0;
    v[0].x = 0.0; v[0].y = 0.0; v[0].z = 0.0;
    m[0] = 0.000295912208287202;
    body_radius[0] = radii[0];
    body_color[0][0] = colors[0][0];
    body_color[0][1] = colors[0][1];
    body_color[0][2] = colors[0][2];

    // Данные для 8 планет: {GM, x, y, z, vx, vy, vz}
    double init_data[8][7] = {
        {4.91248045036476e-11, -0.180936397763299, -0.385147512043214, -0.186946626365338, 0.0202536935932389, -0.00773490114534313, -0.00623293831111668},
        {7.24345233264412e-10,  0.353683857895012, -0.571732114793733, -0.279571521887592, 0.0175395112483742,  0.00930788808619176,  0.00307674893669492},
        {8.88769246302608e-10,  0.824137420123391,  0.509462200362359,  0.220903233023091, -0.00990108870586615, 0.0130341335464125,   0.00565223913525109},
        {9.54954869555077e-11,  1.18050055345426,  -0.641235699022141, -0.326064357170515, 0.00780844617790429,  0.0120318494722293,   0.00530714335760668},
        {2.82534582597155e-07,  1.61789448364208,  -4.5078997844356,   -1.97171823063617,  0.00707630602337104,  0.00255314883766558,  0.000921969043533649},
        {8.45970607324503e-08, -6.51160967219589,  -6.98339235308506,  -2.60414973032736,  0.00390650223666086, -0.00334762436325148, -0.00155039042288},
        {1.29202657962908e-08, -5.39442972255829, -16.7542834486181,   -7.26135557045396,  0.00375053973389636, -0.00116844828942368, -0.000564972666531907},
        {1.52435734788511e-08,  0.539531581225625, -27.9945229379515,  -11.4718890860437,  0.00312669573863622,  9.79150978551659e-05, -3.77167939161073e-05}
    };

    for (int i = 1; i < nBodies; i++) {
        p[i].x = init_data[i-1][1];
        p[i].y = init_data[i-1][2];
        p[i].z = init_data[i-1][3];
        v[i].x = init_data[i-1][4];
        v[i].y = init_data[i-1][5];
        v[i].z = init_data[i-1][6];
        m[i] = init_data[i-1][0];

        body_radius[i] = radii[i];
        body_color[i][0] = colors[i][0];
        body_color[i][1] = colors[i][1];
        body_color[i][2] = colors[i][2];
    }
}

void renderBitmapString(float x, float y, void *font, const char *string) {
    while(*string) {
        glutBitmapCharacter(font, *string);
        string++;
    }
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(camX, camY, camZ,  0.0, 0.0, 0.0,  0.0, 1.0, 0.0);
    for (int i = 0; i < nBodies; i++) {
        glPushMatrix();
        glTranslatef(p[i].x * DRAW_SCALE, p[i].y * DRAW_SCALE, p[i].z * DRAW_SCALE);
        glColor3dv(body_color[i]);
        glutSolidSphere(body_radius[i], 20, 20);
        glPopMatrix();
    }
    
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, winWidth, 0, winHeight);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    char buffer[256];
    glColor3f(1.0, 1.0, 1.0);
    sprintf(buffer, "Barycenter: (%.3e, %.3e, %.3e)", g_barycenter.x, g_barycenter.y, g_barycenter.z);
    renderBitmapString(10, winHeight - 20, GLUT_BITMAP_HELVETICA_18, buffer);
    
    sprintf(buffer, "Total Energy: %.3e", g_totalEnergy);
    renderBitmapString(10, winHeight - 40, GLUT_BITMAP_HELVETICA_18, buffer);
    
    sprintf(buffer, "Total Momentum: (%.3e, %.3e, %.3e)", g_momentum.x, g_momentum.y, g_momentum.z);
    renderBitmapString(10, winHeight - 60, GLUT_BITMAP_HELVETICA_18, buffer);
    
    sprintf(buffer, "Angular Momentum: (%.3e, %.3e, %.3e)", g_angularMomentum.x, g_angularMomentum.y, g_angularMomentum.z);
    renderBitmapString(10, winHeight - 80, GLUT_BITMAP_HELVETICA_18, buffer);
    
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    
    glutSwapBuffers();
}

void reshape(int w, int h) {
    winWidth = w;
    winHeight = h;
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, (double)w / h, 1.0, 10000.0);
    glMatrixMode(GL_MODELVIEW);
}

void timer(int value) {
    predictor_corrector_am8_step();
    calculate_energy_momentum();
    glutPostRedisplay();
    glutTimerFunc(16, timer, 0);
}

void keyboard(unsigned char key, int x, int y) {
    switch(key) {
        case 'w': camZ -= 100; break;
        case 's': camZ += 100; break;
        case 'a': camX -= 100; break;
        case 'd': camX += 100; break;
        case 'q': camY += 100; break;
        case 'e': camY -= 100; break;
        case 27: exit(0); break;
    }
    glutPostRedisplay();
}

int main(int argc, char *argv[]) {
    const size_t initial_capacity = 20000000;
    time_points.reserve(initial_capacity);
    cm_magnitude.reserve(initial_capacity);
    energy_rel.reserve(initial_capacity);
    momentum_rel.reserve(initial_capacity);
    angular_momentum_rel.reserve(initial_capacity);
    srand((unsigned)time(NULL));
    initBodies();

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(winWidth, winHeight);
    glutCreateWindow("Solar System Simulation (3D)");
    
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutTimerFunc(16, timer, 0);
    glutMainLoop();

    
    free(p);
    free(v);
    free(m);
    free(body_radius);
    free(body_color);
    
    return 0;
}
