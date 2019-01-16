#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* #include <mkl.h> */

#ifdef WATCHTIME
#include <time.h>
#endif

int calc_inverse_matrix_of_3x3(double (*mat)[3], double (*out)[3]);

void wwz_c(double t_sample_arr[],
           double omega_sample_arr[],
           double t_data_arr[],
           double x_data_arr[],
           double c_coef,
           double **out_wwz,
           double **out_wwa,
           int    N_t_sample_arr,
           int    N_omega_sample_arr,
           int    N_t_data_arr) {
    int i, j, k, l;

    double omega, tau;
    double sum_weights = 0.;
    double Neff1, Vx1, Vy1;

    double phi_1_arr[N_t_data_arr];
    double phi_2_arr[N_t_data_arr];
    double weights[N_t_data_arr];
    double S_matrix_inv1[3][3] = 
    {{0., 0., 0.},
     {0., 0., 0.},
     {0., 0., 0.}};
    double y_coef_arr1[3] = {0., 0., 0.};
    double y_model_arr1[N_t_data_arr];

    int size = 3;
    int lda = 3;
    /* double matrix_arr[9]; */
    int ipiv;

    double bufd;
    double bufds[5] = {0., 0., 0., 0., 0.};
#ifdef WATCHTIME
    FILE *fptime;
    char tempchar[1024];
    double times[10];
    double tstart;
    int i_time = 0;

    sprintf(tempchar, "time_log_from_c.txt");
    fptime = fopen(tempchar, "w");
    if(fptime==NULL) {
        printf("Can't open file %s \n", tempchar);
        exit(1);
    }

    for (i = 0; i < 10; i++) {
        times[i] = 0.;
    }
#endif

    for (i = 0; i < N_t_sample_arr; i++) {
        tau = t_sample_arr[i];
        for (j = 0; j < N_omega_sample_arr; j++) {
            omega = omega_sample_arr[j];
#ifdef WATCHTIME
            tstart = clock();
#endif
            sum_weights = 0.;
            for (k = 0; k < N_t_data_arr; k++) {
                weights[k] = exp(- c_coef * pow((omega * (t_data_arr[k] - tau)), 2));
                sum_weights += weights[k];
            }
#ifdef DEBUG
            printf ("\n");
            printf ("omega %+e   tau %+e\n", omega, tau);
            printf ("weights\n");
            for (k = 0; k < N_t_data_arr; k++) {
                printf("%+.2lf ", weights[k]);
            }
            printf ("\n");
#endif

#ifdef WATCHTIME
            times[i_time] += (double)(clock() - tstart);
            i_time++;
            tstart = clock();
#endif

            S_matrix_inv1[0][1] = 0.;
            S_matrix_inv1[0][2] = 0.;
            S_matrix_inv1[1][1] = 0.;
            S_matrix_inv1[1][2] = 0.;
            S_matrix_inv1[2][2] = 0.;

            for (k = 0; k < N_t_data_arr; k++) {
                phi_1_arr[k] = cos(omega * (t_data_arr[k] - tau));
                phi_2_arr[k] = sin(omega * (t_data_arr[k] - tau));

                S_matrix_inv1[0][1] += weights[k] * phi_1_arr[k];
                S_matrix_inv1[0][2] += weights[k] * phi_2_arr[k];
                S_matrix_inv1[1][1] += weights[k] * phi_1_arr[k] * phi_1_arr[k];
                S_matrix_inv1[1][2] += weights[k] * phi_1_arr[k] * phi_2_arr[k];
                S_matrix_inv1[2][2] += weights[k] * phi_2_arr[k] * phi_2_arr[k];
            }

            S_matrix_inv1[0][1] /= sum_weights;
            S_matrix_inv1[0][2] /= sum_weights;
            S_matrix_inv1[1][1] /= sum_weights;
            S_matrix_inv1[1][2] /= sum_weights;
            S_matrix_inv1[2][2] /= sum_weights;

            S_matrix_inv1[0][0] = 1.;
            S_matrix_inv1[1][0] = S_matrix_inv1[0][1];
            S_matrix_inv1[2][0] = S_matrix_inv1[0][2];
            S_matrix_inv1[2][1] = S_matrix_inv1[1][2];

            // calc inverse of S_matrix
#ifndef DEBUG
            calc_inverse_matrix_of_3x3(S_matrix_inv1, S_matrix_inv1);
#else
            /* for (k = 0; k < 3; k++) { */
            /*     for (l = 0; l < 3; l++) { */
            /*         S_matrix_inv1[k][l] = 1. + 3.*k + l; */
            /*     } */
            /* } */
            /* S_matrix_inv1[0][0] = 0.; */
            printf("S_matrix\n");
            printf("%+.2lf %+.2lf %+.2lf\n", S_matrix_inv1[0][0], S_matrix_inv1[0][1], S_matrix_inv1[0][2]);
            printf("%+.2lf %+.2lf %+.2lf\n", S_matrix_inv1[1][0], S_matrix_inv1[1][1], S_matrix_inv1[1][2]);
            printf("%+.2lf %+.2lf %+.2lf\n", S_matrix_inv1[2][0], S_matrix_inv1[2][1], S_matrix_inv1[2][2]);
            /* double out_matrix[3][3]; */
            calc_inverse_matrix_of_3x3(S_matrix_inv1, S_matrix_inv1);
            printf("S_matrix_inv1\n");
            printf("%+.2lf %+.2lf %+.2lf\n", S_matrix_inv1[0][0], S_matrix_inv1[0][1], S_matrix_inv1[0][2]);
            printf("%+.2lf %+.2lf %+.2lf\n", S_matrix_inv1[1][0], S_matrix_inv1[1][1], S_matrix_inv1[1][2]);
            printf("%+.2lf %+.2lf %+.2lf\n", S_matrix_inv1[2][0], S_matrix_inv1[2][1], S_matrix_inv1[2][2]);
#endif

#ifdef WATCHTIME
            times[i_time] += (double)(clock() - tstart);
            i_time++;
            tstart = clock();
#endif

            for (k = 0; k < 3; k++) {
                bufds[k] = 0.;
            }
            for (k = 0; k < N_t_data_arr; k++) {
                bufd = weights[k] * x_data_arr[k];
                bufds[0] += bufd;
                bufds[1] += bufd * phi_1_arr[k];
                bufds[2] += bufd * phi_2_arr[k];
            }
            for (k = 0; k < 3; k++) {
                y_coef_arr1[k] = 0.;
                for (l = 0; l < 3; l++) {
                    y_coef_arr1[k] += S_matrix_inv1[k][l]  * bufds[l];
                }
                y_coef_arr1[k] /= sum_weights;
            }

#ifdef WATCHTIME
            times[i_time] += (double)(clock() - tstart);
            i_time++;
            tstart = clock();
#endif

            for (k = 0; k < N_t_data_arr; k++) {
                y_model_arr1[k] = y_coef_arr1[0] + y_coef_arr1[1] * phi_1_arr[k] + y_coef_arr1[2] * phi_2_arr[k];
            }

#ifdef WATCHTIME
            times[i_time] += (double)(clock() - tstart);
            i_time++;
            tstart = clock();
#endif

            for (k = 0; k < 5; k++) {
                bufds[k] = 0.;
            }
            for (k = 0; k < N_t_data_arr; k++) {
                bufds[0] += exp(- c_coef * pow((sqrt(2) * omega * (t_data_arr[k] - tau)), 2));
                bufds[1] += weights[k] * x_data_arr[k];
                bufds[2] += weights[k] * x_data_arr[k] * x_data_arr[k];
                bufds[3] += weights[k] * y_model_arr1[k];
                bufds[4] += weights[k] * y_model_arr1[k] * y_model_arr1[k];
            }

            Neff1 = pow(sum_weights, 2) / bufds[0];
            Vx1   = bufds[2] / sum_weights - pow(bufds[1] / sum_weights, 2);
            Vy1   = bufds[4] / sum_weights - pow(bufds[3] / sum_weights, 2);

            out_wwz[i][j] = 0.5 * (Neff1 - 3.) * Vy1 /  (Vx1 - Vy1);
            out_wwa[i][j] = sqrt(y_coef_arr1[1] * y_coef_arr1[1] + y_coef_arr1[2] * y_coef_arr1[2]);

#ifdef DEBUG
            printf("Neff %+e\n", Neff1);
            printf("Vx1  %+e\n", Vx1);
            printf("Vy1  %+e\n", Vy1);
            printf("wwz  %+e\n", out_wwz[i][j]);
            printf("wwa  %+e\n", out_wwa[i][j]);
#endif

#ifdef WATCHTIME
            times[i_time] += (double)(clock() - tstart);
            i_time++;
            tstart = clock();
            i_time = 0;
#endif
        }
    }

#ifdef WATCHTIME
    for (i = 0; i < 10; i++) {
        fprintf(fptime, "%+e\n", times[i] / CLOCKS_PER_SEC);
    }
    fclose(fptime);
#endif


}


int calc_inverse_matrix_of_3x3(double (*mat)[3], double (*out)[3]){
    int i, j;
    double denominator;
    double buf[3][3];
    denominator =   mat[0][0] * mat[1][1] * mat[2][2]
                  + mat[0][1] * mat[1][2] * mat[2][0]
                  + mat[0][2] * mat[1][0] * mat[2][1]
                  - mat[0][0] * mat[1][2] * mat[2][1]
                  - mat[0][1] * mat[1][0] * mat[2][2]
                  - mat[0][2] * mat[1][1] * mat[2][0];

    if (denominator != 0.) {
        buf[0][0] =   mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
        buf[0][1] = - mat[0][1] * mat[2][2] + mat[0][2] * mat[2][1];
        buf[0][2] =   mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];

        buf[1][0] = - mat[1][0] * mat[2][2] + mat[1][2] * mat[2][0];
        buf[1][1] =   mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0];
        buf[1][2] = - mat[0][0] * mat[1][2] + mat[0][2] * mat[1][0];

        buf[2][0] =   mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];
        buf[2][1] = - mat[0][0] * mat[2][1] + mat[0][1] * mat[2][0];
        buf[2][2] =   mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];

        for(i = 0; i < 3 ; i++) {
            for(j = 0; j < 3 ; j++) {
                out[i][j] = buf[i][j] / denominator;
            }
        }

        return 0;

    } else {
        printf("cannot calc inverse matrix!: denominator == 0\n");
        return -1;

    }

}


#ifdef DEBUG
int main() {
    int i;
    int N_t_sample_arr = 100;
    int N_omega_sample_arr = 20;
    int N_t_data_arr = 60;
    double t_sample_arr[N_t_sample_arr];
    double omega_sample_arr[N_omega_sample_arr];
    double t_data_arr[N_t_data_arr];
    double x_data_arr[N_t_data_arr];
    double c_coef = 0.0125;

    double **wwz, **wwa;
    double *base_wwz, *base_wwa;
    wwz = malloc(sizeof(double *) * N_t_sample_arr);
    wwa = malloc(sizeof(double *) * N_t_sample_arr);
    base_wwz = malloc(sizeof(double) * N_t_sample_arr * N_omega_sample_arr);
    base_wwa = malloc(sizeof(double) * N_t_sample_arr * N_omega_sample_arr);
    for (i = 0; i < N_t_sample_arr; i++) {
        wwz[i] = base_wwz + i * N_omega_sample_arr;
        wwa[i] = base_wwa + i * N_omega_sample_arr;
    }
    /* double wwa[N_t_sample_arr][N_omega_sample_arr]; */
    /* double (*wwz)[N_omega_sample_arr]; */
    /* double (*wwa)[N_omega_sample_arr]; */
    /* wwz = malloc(sizeof(double) * N_t_sample_arr * N_omega_sample_arr); */
    /* wwa = malloc(sizeof(double) * N_t_sample_arr * N_omega_sample_arr); */

    srand((unsigned int)time(NULL));
    for (i = 0; i < N_t_data_arr; i++) {
        t_data_arr[i] = i * 0.208 + 0.001 * (double)rand() / RAND_MAX;
        x_data_arr[i] = 3 * sin(1. * t_data_arr[i]) + 0.001 * (double)rand() / RAND_MAX;
    }
    for (i = 0; i < N_t_sample_arr; i++) {
        t_sample_arr[i] = i * 0.125;
    }
    for (i = 0; i < N_omega_sample_arr; i++) {
        omega_sample_arr[i] = (i+1) *0.25;
    }

    /* t_sample_arr[0] = 5.; */
    /* omega_sample_arr[0] = 1.; */

    printf("\n");
    for (i = 0; i < N_t_data_arr; i++) {
        printf("%+e %+e\n", t_data_arr[i], x_data_arr[i]);
    }
    printf("\n");

    wwz_c(t_sample_arr,
          omega_sample_arr,
          t_data_arr,
          x_data_arr,
          c_coef,
          wwz,
          wwa,
          N_t_sample_arr,
          N_omega_sample_arr,
          N_t_data_arr);

    free(wwa);
    free(wwz);

    printf("end\n");
}
#endif
