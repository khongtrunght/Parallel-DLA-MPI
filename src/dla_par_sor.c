#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define iter 800
#define N 200
#define tol 0.001
#define omega 1.9
#define c0 1
#define cN 0

int i, j, k, r, rank, size, Np;

float *C_current, *C_old, *alpha;
float *C_above, *C_below;
int *O, *O_below, *candidates, *O_above;
float *nutri;

char snum[20];
char fpath[40] = "output/";
FILE *fp;
MPI_Status status;

//=========================
float r2()
{
    return (float)rand() / (float)RAND_MAX;
}
void KhoiTao(float *C_old, float *C_current, int *O, int *candidates)
{
    // C  Initialization
    for (i = 0; i < Np; ++i)
        for (j = 0; j < N; ++j)
        {
            *(C_current + i * N + j) = 0;
            *(C_old + i * N + j) = 0;
            *(O + i * N + j) = 0;
        }

    // Object Initialization
    if (rank == size - 1)
        *(O + (Np - 1) * N + N / 2) = 1;
}
//=========================
void SOR()
{
    float right, left, up, down;
    float global_alpha;
    int dieu_chinh;
    MPI_Datatype everytwice;
    MPI_Type_vector((N + 1) / 2, 1, 2, MPI_FLOAT, &everytwice);
    MPI_Type_commit(&everytwice);

    do
    {
        *alpha = 0;
        for (r = 0; r < 2; r++)
        {
            // He so cho tung Np khac nhau
            dieu_chinh = (Np * rank + 1 + r) % 2;

            // Send down
            if (rank != size - 1)
                MPI_Send(C_current + (Np - 1) * N + (dieu_chinh + Np) % 2, 1, everytwice, rank + 1, 0, MPI_COMM_WORLD);

            // Send up
            if (rank != 0)
                MPI_Send(C_current + 1 - dieu_chinh, 1, everytwice, rank - 1, 1, MPI_COMM_WORLD);

            // Receive from below
            if (rank != size - 1)
                MPI_Recv(C_below, (N + 1) / 2, MPI_FLOAT, rank + 1, 1, MPI_COMM_WORLD, &status);
            else
            {
                for (i = 0; i < (N + 1) / 2; ++i)
                    C_below[i] = cN;
            }

            // Receive from above
            if (rank != 0)
                MPI_Recv(C_above, (N + 1) / 2, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &status);
            else
            {
                for (i = 0; i < (N + 1) / 2; ++i)
                    C_above[i] = c0;
            }

            for (i = 0; i < Np; i++)
            {
                for (j = 0; j < N; j++)
                {
                    // ignore object position
                    if (*(O + i * N + j) == 1)
                        continue;

                    // ignore half of the grid
                    if ((rank * Np + i + j) % 2 == r)
                        continue;

                    left = (j == 0) ? *(C_current + i * N + (N - 1)) : *(C_current + i * N + j - 1);
                    right = (j == N - 1) ? *(C_current + i * N) : *(C_current + i * N + j + 1);
                    up = (i == 0) ? *(C_above + j / 2) : *(C_current + (i - 1) * N + j);
                    down = (i == Np - 1) ? *(C_below + j / 2) : *(C_current + (i + 1) * N + j);

                    *(C_current + i * N + j) = (1 - omega) * *(C_current + i * N + j) + omega * (left + right + up + down) / 4;

                    *alpha = fmax(*alpha, fabs(*(C_current + i * N + j) - *(C_old + i * N + j)));
                }
            }
        }

        for (i = 0; i < Np; ++i)
            for (j = 0; j < N; ++j)
                *(C_old + i * N + j) = *(C_current + i * N + j);

        MPI_Allreduce(alpha, &global_alpha, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    } while (global_alpha > tol);

    MPI_Type_free(&everytwice);
}
//=========================
void growth()
{
    int left, right, up, down;
    float global_nutri;

    for (i = 0; i < N; ++i)
    {
        *(O_below + i) = 0;
        *(O_above + i) = 0;
    }

    // send the first row of object to above
    if (rank != 0)
        MPI_Send(O, N, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);

    // send the last row of object to below
    if (rank != size - 1)
        MPI_Send(O + (Np - 1) * N, N, MPI_INT, rank + 1, 1, MPI_COMM_WORLD);

    // receive object from below
    if (rank != size - 1)
        MPI_Recv(O_below, N, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, &status);

    // receive object from above
    if (rank != 0)
        MPI_Recv(O_above, N, MPI_INT, rank - 1, 1, MPI_COMM_WORLD, &status);

    // reset candidates
    for (i = 0; i < Np; ++i)
        for (j = 0; j < N; ++j)
            *(candidates + i * N + j) = 0;

    // reset nutri
    *nutri = 0;

    // calculate candidates and nutri
    for (i = 0; i < Np; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (*(O + i * N + j) == 1)
            {
                continue;
            }

            left = (j == 0) ? 0 : *(O + i * N + j - 1);
            right = (j == N - 1) ? 0 : *(O + i * N + j + 1);
            up = (i == 0) ? *(O_above + j) : *(O + (i - 1) * N + j);
            down = (i == Np - 1) ? *(O_below + j) : *(O + (i + 1) * N + j);

            if ((left == 1 || right == 1 || up == 1 || down == 1))
            {
                *(candidates + i * N + j) = 1;
                *nutri += *(C_current + i * N + j);
            }
            else
                *(candidates + i * N + j) = 0;
        }
    }

    // Reduce
    MPI_Allreduce(nutri, &global_nutri, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    for (i = 0; i < Np; ++i)
        for (j = 0; j < N; ++j)
            if (*(candidates + i * N + j) == 1 && r2() <= (*(C_current + i * N + j) / global_nutri))
            {
                *(O + i * N + j) = 1;
                *(C_current + i * N + j) = 0;
            }
}
//=========================
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Np = N / size;

    C_current = (float *)malloc(Np * N * sizeof(float));
    C_old = (float *)malloc(Np * N * sizeof(float));
    alpha = (float *)malloc(sizeof(float));

    C_above = (float *)malloc((N + 1) / 2 * sizeof(float));
    C_below = (float *)malloc((N + 1) / 2 * sizeof(float));

    O = (int *)malloc(Np * N * sizeof(int));
    O_below = (int *)malloc(N * sizeof(int));
    O_above = (int *)malloc(N * sizeof(int));

    candidates = (int *)malloc(Np * N * sizeof(int));

    nutri = (float *)malloc(sizeof(float));

    KhoiTao(C_old, C_current, O, candidates);

    for (k = 0; k < iter; k++)
    {
        SOR();
        for (i = 0; i < Np; ++i)
            for (j = 0; j < N; ++j)
                *(C_old + i * N + j) = *(C_current + i * N + j);
        growth();
    }

    float *C_current_global;
    int *O_global;
    // print
    if (rank == 0)
    {
        C_current_global = (float *)malloc(N * N * sizeof(float));
        O_global = (int *)malloc(N * N * sizeof(int));
    }
    MPI_Gather(C_current, Np * N, MPI_FLOAT, C_current_global, Np * N, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(O, Np * N, MPI_INT, O_global, Np * N, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("Printing");
        sprintf(snum, "%d", N);
        strcat(snum, "par.txt");
        strcat(fpath, snum);
        fp = fopen(fpath, "w");
        for (i = 0; i < N; ++i)
        {
            for (j = 0; j < N; ++j)
            {
                if (*(O_global + i * N + j) == 1)
                    *(C_current_global + i * N + j) = 1.0;
                fprintf(fp, "%f ", *(C_current_global + i * N + j));
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }

    MPI_Finalize();
    return 0;
}