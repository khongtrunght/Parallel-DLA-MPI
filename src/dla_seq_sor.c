#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define iter 800
#define N 200
#define tol 0.001
#define omega 1.9
#define c0 1
#define cN 0

int i, j, k, r;

float *C_current, *C_old, *alpha;
int *O, *candidates;
float *nutri;

char snum[10];
char fpath[40] = "output/";
FILE *fp;

//=========================
float r2()
{
    return (float)rand() / (float)RAND_MAX;
}
//=========================

void DisplayMatrix(float *A, int row, int col)
{
    int i, j;
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
            printf("%lf\t", *(A + i * col + j));
        printf("\n");
    }
}
//=========================
void KhoiTao(float *C_old, float *C_current, int *O, int *candidates)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            if (i == 0)
                *(C_old + i * N + j) = 1;
            else
                *(C_old + i * N + j) = 0;
            *(C_current + i * N + j) = *(C_old + i * N + j);
            *(O + i * N + j) = 0;
            *(candidates + i * N + j) = 0;
        }

    *(O + (N - 1) * N + N / 2) = 1;
}
//=========================
void SOR(float *C_old, float *C_current, int *O)
{
    float delta;
    float right, left, up, down;

    do
    {
        delta = 0;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                if (*(O + i * N + j) == 1)
                    continue;

                left = (j == 0) ? *(C_current + i * N + (N - 1)) : *(C_current + i * N + j - 1);
                right = (j == N - 1) ? *(C_current + i * N + 0) : *(C_current + i * N + j + 1);
                up = (i == 0) ? c0 : *(C_current + (i - 1) * N + j);
                down = (i == N - 1) ? cN : *(C_current + (i + 1) * N + j);
                *(C_current + i * N + j) = 0.25 * omega * (left + right + up + down) + (1 - omega) * *(C_current + i * N + j);
                delta = fmax(delta, fabs(*(C_current + i * N + j) - *(C_old + i * N + j)));
            }
        // Update C_old
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                *(C_old + i * N + j) = *(C_current + i * N + j);

    } while (delta > tol);
}
//=========================
void growth(float *C_current, int *O, int *candidates)
{
    int left, right, up, down;

    *nutri = 0;

    // reset candidates
    for (i = 0; i < N; ++i)
        for (j = 0; j < N; ++j)
            *(candidates + i * N + j) = 0;

    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
        {
            if (*(O + i * N + j) == 1)
                continue;

            left = (j == 0) ? 0 : *(O + i * N + j - 1);
            right = (j == N - 1) ? 0 : *(O + i * N + j + 1);
            up = (i == 0) ? 0 : *(O + (i - 1) * N + j);
            down = (i == N - 1) ? 0 : *(O + (i + 1) * N + j);
            if ((left == 1 || right == 1 || up == 1 || down == 1))
            {
                *(candidates + i * N + j) = 1;
                *nutri += *(C_current + i * N + j);
            }
            else
                *(candidates + i * N + j) = 0;
        }

    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            if ((*(candidates + i * N + j) == 1) && (r2() <= (*(C_current + i * N + j) / *nutri)))
            {
                *(O + i * N + j) = 1;
                *(C_current + i * N + j) = 0;
            }
}

//=========================

int main()
{
    srand(time(NULL));
    sprintf(snum, "%d", N);
    strcat(snum, ".txt");
    strcat(fpath, snum);
    fp = fopen(fpath, "w");

    C_current = (float *)malloc(sizeof(float) * N * N);
    C_old = (float *)malloc(sizeof(float) * N * N);
    O = (int *)malloc(sizeof(int) * N * N);
    candidates = (int *)malloc(sizeof(int) * N * N);
    nutri = (float *)malloc(sizeof(float));

    KhoiTao(C_old, C_current, O, candidates);

    for (int k = 0; k < iter; k++)
    {

        SOR(C_old, C_current, O);
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                *(C_old + i * N + j) = *(C_current + i * N + j);
        growth(C_current, O, candidates);
    }

    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (*(O + i * N + j) == 1)
                *(C_current + i * N + j) = 1.0;
            fprintf(fp, "%lf\t", *(C_current + i * N + j));
        }
        fprintf(fp, "\n");
    }
    return 0;
}
