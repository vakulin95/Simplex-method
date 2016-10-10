#include <stdio.h>
#include <stdlib.h>

#define FILENAME "matrix4_max.txt"
#define MAXIMIZATION 1

int read(double**, int*, FILE*);
int write(double**, int*);
int def_nbas(int*, int*);
int s_lead_column(double**, int*);
int s_lead_line(double**, int*, int);
int simpl_t_conv(double**, int*, int, int);
int answer(double**, int*);

int exN, exM, exB;

int main()
{
    int i;
    int *basis, *not_basis;
    int lc, ll;
    double **mass;
    FILE *in;

    if((in = fopen(FILENAME, "r")) == NULL)
    {
        printf("Opening file error!\n");
        return 1;
    }

    fscanf(in, "%d %d %d", &exN, &exM, &exB);

    basis = (int*)malloc(exB * sizeof(int));
    not_basis = (int*)malloc((exM - exB - 1) * sizeof(int));
    mass = (double**)malloc(exN * sizeof(double*));
    for (i = 0; i < exN; i++)
        mass[i] = (double*)malloc(exM * sizeof(double));

    read(mass, basis, in);
    fclose(in);
    write(mass, basis);

    if(!def_nbas(basis, not_basis))
        return 0;

    while((lc = s_lead_column(mass, basis)) != -1)
    {
        printf("leading column: %d\nleading line: %d\n", lc, ll = s_lead_line(mass, basis, lc));
        simpl_t_conv(mass, basis, ll, lc);
        write(mass, basis);
    }
    answer(mass, not_basis);

    for(i = 0; i < exN; i++)
        free(mass[i]);
    free(mass);
    free(basis);

    return 0;
}

int read(double **x, int *basis, FILE *in)
{
    int i, j;
    char c;

    for(i = 0; i < exB; i++)
        basis[i] = exM - exB - 1 + i;

    for(i = 0; i < exN; i++)
        for(j = 0; j < exM; j++)
            fscanf(in, "%lf", &x[i][j]);

    return 0;
}

int write(double **x, int *basis)
{
    int i, j;

    printf("\nBasis: ");
    for(i = 0; i < exB; i++)
        printf("%d ", basis[i]);
    printf("\n");

    printf("Matrix:\n");
    for(i = 0; i < exN; ++i)
    {
        for(j = 0; j < exM; j++)
            printf("%6.2lf", x[i][j]);
        printf("\n");
    }

    return 0;
}

int def_nbas(int *basis, int *not_basis)
{
    int i, j, k;
    int cond = 0;

    for(i = 0, k = 0; i < exM - 1; i++)
    {
        for(j = 0; j < exB; j++)
        {
            if(i == basis[j])
            {
                cond = 0;
                break;
            }
            cond = 1;
        }
        if(cond)
        {
            if(k >= exM - exB - 1)
            {
                printf("ERROR in def_nbas()\n");
                return 0;
            }
            not_basis[k] = i;
            k++;
        }
    }

    return 1;
}

int s_lead_column(double **mass, int *basis)
{
    int i, j, cond = 0;
    int Y = 0;
    double minmax = mass[exN - 1][0];

    for(i = 1; i < exM - 1; i++)
    {
        for(j = 0; j < exB; j++)
        {
            if(i == basis[j])
            {
                cond = 0;
                break;
            }
            cond = 1;
        }
        if(cond)
        {
            if(MAXIMIZATION && mass[exN - 1][i] < 0 && mass[exN - 1][i] < minmax)
                Y = i;
            else if(!MAXIMIZATION && (int)mass[exN - 1][i] > 0 && mass[exN - 1][i] > minmax)
                Y = i;
        }
    }

    if((MAXIMIZATION && Y == 0 && minmax >= 0) || ( !MAXIMIZATION && Y == 0 && minmax <= 0 ))
    {
        printf("No leading column!\n");
        return -1;
    }

    return Y;
}

int s_lead_line(double **mass, int *basis, int lc)
{
    int i;
    int Y = -1;
    double min = 1000; 

    for(i = 0; i < exN - 1; i++)
    {
        if(mass[i][lc] <= 0)
            continue;
        if(mass[i][exM - 1] / mass[i][lc] < min)
        {
            min = mass[i][exM - 1] / mass[i][lc];
            //printf("%lf\n", min);
            Y = i;
        }
    }

    if(Y == -1)
        printf("No leading line!\n");

    return Y;
}

int simpl_t_conv(double **mass, int *basis, int ll, int lc)
{
    int i, j;
    double coeff = mass[ll][lc];

    //new basis[]
    basis[ll] = lc;

    //new leading line calculation
    for(i = 0; i < exM; i++)
        mass[ll][i] /= coeff;

    //other lines calculation
    for(i = 0; i < exN; i++)
    {
        if(i == ll)
            continue;

        coeff = mass[i][lc];

        for(j = 0; j < exM; j++)
        {
            mass[i][j] -= coeff*mass[ll][j];
        }
    }
}

int answer(double **mass, int *not_basis)
{
    int i;

    printf("\nAnswer:\n");

    printf("Z* = %.2lf\nX* = ( ", mass[exN - 1][exM - 1]);

    for(i = 0; i < exM - exB - 1; i++)
        printf("%.2lf; ", mass[not_basis[i]][exM - 1]);
    printf("\b\b )\n");

    return 0;
}
