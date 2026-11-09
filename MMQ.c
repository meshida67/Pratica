#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//-----------------------------------------------------------------------------
#define degrau 1000
#define DEGRAU 100000

typedef struct
{
    double *coeficientes;
    double *incertezas;
    int n_parametros;
} ResultadoAjuste;
//-----------------------------------------------------------------------------
double** multMatriz(double **A, double **B, int l, int m, int n)
{
    double **AB = (double**)malloc(l * sizeof(double*));
    for (int i = 0; i < l; i++)
        AB[i] = (double*)malloc(n * sizeof(double));

    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < n; j++)
        {
             AB[i][j] = 0;
             for (int k = 0; k < m; k++)
                AB[i][j] += A[i][k] * B[k][j];
        }
    }

    return AB;
}
//-----------------------------------------------------------------------------
double **Transposta(double **A, int m, int n)
{
    double **At = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++)
        At[i] = (double*)malloc( m * sizeof(double));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
            At[i][j] = A[j][i];
    }
    return At;
}
//-----------------------------------------------------------------------------
double *elimGauss(double **A, double **b, int n)
{
    double *x = (double*)malloc(n * sizeof(double)); 
    for (int i = 0; i < n; i++)
    {
        int linha_max = i;
        for (int j = i + 1; j < n; j++)
        {
            if (fabs(A[j][i]) > fabs(A[i][i]))
                linha_max = j;
        }

        if (linha_max != i)
        {
            for (int j = 0; j < n; j++)
            {
                double tmp = A[i][j];
                A[i][j] = A[linha_max][j];
                A[linha_max][j] = tmp;
            }
            double tmp = b[i][0];
            b[i][0] = b[linha_max][0];
            b[linha_max][0] = tmp;
        }

        for (int j = i + 1; j < n; j++)
        {
            double fator = A[j][i] / A[i][i];
            for (int k = i; k < n; k++)
                A[j][k] -= fator * A[i][k];
            b[j][0] -= fator * b[i][0];
        }
    }

    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = b[i][0];
        for (int j = i + 1; j < n; j++)
            x[i] -= x[j] * A[i][j];
        x[i] /= A[i][i];
    }

    return x;
}
//-----------------------------------------------------------------------------
double **InvMatriz(double **A, int n)
{
    double **Inv = (double**)malloc(n * sizeof (double*));
    for (int i = 0; i < n; i++)
    {
        Inv[i] = (double*)malloc(n * sizeof(double));
        double **A_copia = (double**)malloc(n * sizeof(double*));
        double **Id = (double**)malloc(n * sizeof(double*));
        for (int j = 0; j < n; j++)
        {
            A_copia[j] = (double*)malloc(n * sizeof(double));
            Id[j] = (double*)malloc(n * sizeof(double));
            for (int k = 0; k < n; k++)
            {
                A_copia[j][k] = A[j][k];
                if (i == j)
                    Id[j][0] = 1.0;
                else
                    Id[j][0] = 0.0;
            }
        }
        double *col_Inv = elimGauss(A_copia, Id, n);
        for (int j = 0; j < n; j++)
        {
            Inv[i][j] = col_Inv[i];
            free(A_copia[j]);
            free(Id[j]);
        }
        free(A_copia);
        free(Id);
        free(col_Inv);
    }
    return Inv;
}
//-----------------------------------------------------------------------------
ResultadoAjuste ajusteQuad(double **A, int m, int n)
{
    int mn = n * m;
    double **X = (double**)malloc(mn * sizeof(double*));
    double **y = (double**)malloc(mn * sizeof(double*));
    for (int i = 0; i < mn; i++)
    {
        X[i] = (double*)malloc(3 * sizeof(double));
        y[i] = (double*)malloc(sizeof(double));
    }

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
            y[n * i + j][0] = A[i][j];
    }

    for (int i = 0; i < m; i++)
    {
        double normal = (double)(i + 1);
        for (int j = 0; j < n; j++)
        {
            int indice = n * i + j;
            X[indice][0] = 1.0;
            X[indice][1] = normal;
            X[indice][2] = normal * normal;
        }
    }

    double **Xt = Transposta(X, n * m, 3);
    double **cov = multMatriz(Xt, X, 3, n * m, 3);
    double **y_ = multMatriz(Xt, y, 3, n * m, 1);

    double **cov_copia = (double**)malloc(3 * sizeof(double*));
    for(int i = 0; i < 3; i++) {
        cov_copia[i] = (double*)malloc(3 * sizeof(double));
        for(int j = 0; j < 3; j++)
            cov_copia[i][j] = cov[i][j];
    }

    double *b = elimGauss(cov, y_, 3);

    b[1] = b[1] / degrau;
    b[2] = b[2] / (degrau * degrau);

    double SQe = 0.0;
    for (int i = 0; i < m; i++)
    {
        double x_val = (double)(i + 1) * degrau;
        double pred = b[0] + b[1] * x_val + b[2] * x_val * x_val;
        for (int j = 0; j < n; j++)
        {
            double dif = A[i][j] - pred;
            SQe += dif * dif;
        }
    }

    double sigma2 = SQe / (mn - 3);

    double **invCov = InvMatriz(cov_copia, 3);




    ResultadoAjuste ajuste;
    ajuste.coeficientes = b;
    ajuste.incertezas = (double*)malloc(3 * sizeof(double));
    ajuste.incertezas[0] = sqrt(sigma2 * invCov[0][0]);
    ajuste.incertezas[1] = sqrt(sigma2 * invCov[1][1]) / degrau;
    ajuste.incertezas[2] = sqrt(sigma2 * invCov[2][2]) / (degrau * degrau);

    for (int i = 0; i < mn; i++)
    {
        free(X[i]);
        free(y[i]);
    }
    for (int i = 0; i < 3; i++)
    {
        free(cov_copia[i]);
        free(cov[i]);
        free(y_[i]);
    }
    free(X);
    free(y);
    free(cov);
    free(cov_copia);
    free(y_);

    return ajuste;
}
//-----------------------------------------------------------------------------
ResultadoAjuste ajusteNlogN(double **A, int m, int n)
{
    int mn = n * m;
    double **X = (double**)malloc(mn * sizeof(double*));
    double **y = (double**)malloc(mn * sizeof(double*));
    for (int i = 0; i < mn; i++)
    {
        X[i] = (double*)malloc(2 * sizeof(double));
        y[i] = (double*)malloc(sizeof(double));
    }

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
            y[n * i + j][0] = A[i][j];
    }

    for (int i = 0; i < m; i++)
    {
        double normal = (double)(i + 1) * DEGRAU;
        double n_logn = normal * (log(normal) / log(2.0));
        for (int j = 0; j < n; j++)
        {
            int indice = n * i + j;
            X[indice][0] = 1.0;
            X[indice][1] = n_logn / DEGRAU;
        }
    }


    double **Xt = Transposta(X, n * m, 2);
    double **cov = multMatriz(Xt, X, 2, n * m, 2);
    double **y_ = multMatriz(Xt, y, 2, n * m, 1);

    double **cov_copia = (double**)malloc(2 * sizeof(double*));
    for(int i = 0; i < 2; i++)
    {
        cov_copia[i] = (double*)malloc(2 * sizeof(double));
        for(int j = 0; j < 2; j++)
            cov_copia[i][j] = cov[i][j];
    }

    double *b = elimGauss(cov, y_, 2);

    b[1] = b[1] / DEGRAU;

    double SQe = 0.0;
    for (int i = 0; i < m; i++)
    {
        double x_val = (double)(i + 1) * DEGRAU;
        double pred = b[0] + b[1] * x_val * (log(x_val) / log(2.0));
        for (int j = 0; j < n; j++)
        {
            double dif = A[i][j] - pred;
            SQe += dif * dif;
        }
    }

    double sigma2 = SQe / (mn - 2);

    double **invCov = InvMatriz(cov_copia, 2);

    ResultadoAjuste ajuste;
    ajuste.coeficientes = b;
    ajuste.incertezas = (double*)malloc(2 * sizeof(double));
    ajuste.incertezas[0] = sqrt(sigma2 * invCov[0][0]);
    ajuste.incertezas[1] = sqrt(sigma2 * invCov[1][1]) / DEGRAU;


    for (int i = 0; i < mn; i++)
    {
        free(X[i]);
        free(y[i]);
    }
    free(X);
    free(y);

    return ajuste;
}
//-----------------------------------------------------------------------------
double modeloQuad(double x, double a, double b, double c)
{
    x *= degrau;
    return a + b * x + c * x * x;
}
//-----------------------------------------------------------------------------
double modeloNlogN(double x, double a, double b)
{
    x *= DEGRAU;
    return a + b * x * (log(x) / log(2.0));
}
//-----------------------------------------------------------------------------
int main()
{
    FILE *dados = fopen("dados.txt", "r");

    double **Bubble;
    double **Quick;
    double **Merge;

    double **comparaBubble;
    double **comparaQuick;
    double **comparaMerge;

    Bubble = (double**)malloc(20 * sizeof(double*));
    Quick = (double**)malloc(20 * sizeof(double*));
    Merge = (double**)malloc(20 * sizeof(double*));
    comparaBubble = (double**)malloc(20 * sizeof(double*));
    comparaQuick = (double**)malloc(20 * sizeof(double*));
    comparaMerge = (double**)malloc(20 * sizeof(double*));



    for (int i = 0; i < 20; i++)
    {
        Bubble[i] = (double*)malloc(100 * sizeof(double));
        comparaBubble[i] = (double*)malloc(sizeof(double));
        for (int j = 0; j < 100; j++)
            fscanf(dados, "%lf %lf ", &Bubble[i][j], &comparaBubble[i][0]);
    }
    ResultadoAjuste coefBubble = ajusteQuad(Bubble, 20 ,100);

    for (int i = 0; i < 20; i++)
    {
        Quick[i] = (double*)malloc(100 * sizeof(double));
        comparaQuick[i] = (double*)malloc(10 * sizeof(double));
        for (int j = 0; j < 100; j++)
            fscanf(dados, "%lf  %lf ", &Quick[i][j], &comparaQuick[i][j / 10]);
    }
    ResultadoAjuste coefQuick = ajusteNlogN(Quick, 20, 100);

    for (int i = 0; i < 20; i++)
    {
        Merge[i] = (double*)malloc(100 * sizeof(double));
        comparaMerge[i] = (double*)malloc(10 * sizeof(double));
        for (int j = 0; j < 100; j++)
            fscanf(dados, "%lf %lf ", &Merge[i][j], &comparaMerge[i][j / 10]);
    }
    ResultadoAjuste coefMerge = ajusteNlogN(Merge, 20, 100);

    printf("Bubble: %e + %e * n + %e * n * n\n", coefBubble.coeficientes[0], coefBubble.coeficientes[1], coefBubble.coeficientes[2]);
    printf("Incertezas Bubble: %e, %e, %e\n",coefBubble.incertezas[0], coefBubble.incertezas[1], coefBubble.incertezas[2]);
    printf("Quick: %e + %e * n * log(n)\n", coefQuick.coeficientes[0], coefQuick.coeficientes[1]);
    printf("Incertezas Quick: %e, %e\n", coefQuick.incertezas[0], coefQuick.incertezas[1]);
    printf("Merge: %e + %e * n * log(n)\n", coefMerge.coeficientes[0], coefMerge.coeficientes[1]);
    printf("Incertezas Merge: %e, %e\n", coefMerge.incertezas[0], coefMerge.incertezas[1]);
//-----------------------------------------------------------------------------
double SQR_bubble = 0.0;
double SQe_bubble = 0.0;
double SQT_bubble = 0.0;

double media_bubble[20];
double media_global_bubble = 0.0;
for (int i = 0; i < 20; i++)
{
    media_bubble[i] = 0.0;
    for (int j = 0; j < 100; j++)
        media_bubble[i] += Bubble[i][j];
    media_global_bubble += media_bubble[i];
    media_bubble[i] /= 100;
}
media_global_bubble /= 2000;
for (int i = 0; i < 20; i++)
{
    for (int j = 0; j < 100; j++)
    {
        double tmp = Bubble[i][j] - media_global_bubble;
        SQT_bubble += tmp * tmp;

        tmp = modeloQuad(i + 1, coefBubble.coeficientes[0], coefBubble.coeficientes[1], coefBubble.coeficientes[2]);
        SQR_bubble += (tmp - media_global_bubble) * (tmp - media_global_bubble);
        tmp = Bubble[i][j] - media_bubble[i];
        SQe_bubble += tmp * tmp;
    }
}
double R2 = SQR_bubble / SQT_bubble;
double R2max = 1 - SQe_bubble / SQT_bubble;
R2 *= 100;
R2max *= 100;
printf("R2 Bubble: %lf %%\nR2max Bubble: %lf %%\n", R2, R2max);
//-----------------------------------------------------------------------------
double SQR_quick = 0.0;
double SQe_quick = 0.0;
double SQT_quick = 0.0;

double media_quick[20];
double media_global_quick = 0.0;
for (int i = 0; i < 20; i++)
{
    media_quick[i] = 0.0;
    for (int j = 0; j < 100; j++)
        media_quick[i] += Quick[i][j];
    media_global_quick += media_quick[i];
    media_quick[i] /= 100;
}
media_global_quick /= 2000;
for (int i = 0; i < 20; i++)
{
    for (int j = 0; j < 100; j++)
    {
        double tmp = Quick[i][j] - media_global_quick;
        SQT_quick += tmp * tmp;

        tmp = modeloNlogN(i + 1, coefQuick.coeficientes[0], coefQuick.coeficientes[1]);
        SQR_quick += (tmp - media_global_quick) * (tmp - media_global_quick);
        tmp = Quick[i][j] - media_quick[i];
        SQe_quick += tmp * tmp;
    }
}
R2 = SQR_quick / SQT_quick;
R2max = 1 - SQe_quick / SQT_quick;
R2 *= 100;
R2max *= 100;
printf("R2 Quick: %lf %%\nR2max Quick: %lf %%\n", R2, R2max);
//-----------------------------------------------------------------------------
double SQR_merge = 0.0;
double SQe_merge = 0.0;
double SQT_merge = 0.0;

double media_merge[20];
double media_global_merge = 0.0;
for (int i = 0; i < 20; i++)
{
    media_merge[i] = 0.0;
    for (int j = 0; j < 100; j++)
        media_merge[i] += Merge[i][j];
    media_global_merge += media_merge[i];
    media_merge[i] /= 100;
}
media_global_merge /= 2000;
for (int i = 0; i < 20; i++)
{
    for (int j = 0; j < 100; j++)
    {
        double tmp = Merge[i][j] - media_global_merge;
        SQT_merge += tmp * tmp;

        tmp = modeloNlogN(i + 1, coefMerge.coeficientes[0], coefMerge.coeficientes[1]);
        SQR_merge += (tmp - media_global_merge) * (tmp - media_global_merge);
        tmp = Merge[i][j] - media_merge[i];
        SQe_merge += tmp * tmp;
    }
}
R2 = SQR_merge / SQT_merge;
R2max = 1 - SQe_merge / SQT_merge;
R2 *= 100;
R2max *= 100;
printf("R2 Merge: %lf %%\nR2max Merge: %lf %%\n", R2, R2max);
//-----------------------------------------------------------------------------
double media_compara_quick[20];
double media_compara_merge[20];
double sigma_compara_quick[20];
double sigma_compara_merge[20];
printf("Comparacoes Bubble:\n");
for (int i = 0; i < 20; i++)
{
    printf("%lf\n", comparaBubble[i][0]);
    media_compara_quick[i] = 0.0;
    media_compara_merge[i] = 0.0;
    for (int j = 0; j < 10; j++)
    {
        media_compara_quick[i] += (double)comparaQuick[i][j];
        media_compara_merge[i] += (double)comparaMerge[i][j];
    }
    media_compara_quick[i] /= 10.0;
    media_compara_merge[i] /= 10.0;

    sigma_compara_quick[i] = 0.0;
    sigma_compara_merge[i] = 0.0;
    for (int j = 0; j < 10; j++)
    {
        double tmp = (double)comparaQuick[i][j] - media_compara_quick[i];
        sigma_compara_quick[i] += tmp * tmp;
        tmp = comparaMerge[i][j] - media_compara_merge[i];
        sigma_compara_merge[i] += tmp * tmp;
    }
    sigma_compara_quick[i] /= 9;
    sigma_compara_merge[i] /= 9;
    sigma_compara_quick[i] = sqrt((float)sigma_compara_quick[i]);
    sigma_compara_merge[i] = sqrt((float)sigma_compara_merge[i]);

}
printf("Comparacoes Quick:\n");
for (int i = 0; i < 20; i++)
{
    int a = (int)media_compara_quick[i];
    int b = (int)sigma_compara_quick[i];
    printf("%d +- %d\n", a, b);
}
printf("Comparacoes Merge:\n");
for (int i = 0; i < 20; i++)
{
        int a = (int)media_compara_merge[i];
        int b = (int)sigma_compara_merge[i];
    printf("%d +- %d\n", a, b);
}
//-----------------------------------------------------------------------------
ResultadoAjuste coef_compara_bubble = ajusteQuad(comparaBubble, 20, 1);
ResultadoAjuste coef_compara_quick = ajusteNlogN(comparaQuick, 20, 10);
ResultadoAjuste coef_compara_merge = ajusteNlogN(comparaMerge, 20, 10);

printf("Compara Bubble: %e + %e * n + %e * n * n\n", coef_compara_bubble.coeficientes[0], coef_compara_bubble.coeficientes[1], coef_compara_bubble.coeficientes[2]);
printf("Incertezas Compara Bubble: %e, %e, %e\n", coef_compara_bubble.incertezas[0], coef_compara_bubble.incertezas[1], coef_compara_bubble.incertezas[2]);
printf("Compara Quick: %e + %e * n * log(n)\n", coef_compara_quick.coeficientes[0], coef_compara_quick.coeficientes[1]);
printf("Incertezas Compara Quick: %e, %e\n", coef_compara_quick.incertezas[0], coef_compara_quick.incertezas[1]);
printf("Compara Merge: %e + %e * n * log(n)\n", coef_compara_merge.coeficientes[0], coef_compara_merge.coeficientes[1]);
printf("Incertezas Compara Merge: %e, %e\n", coef_compara_merge.incertezas[0], coef_compara_merge.incertezas[1]);
//-----------------------------------------------------------------------------
double SQR_compara = 0.0;
double SQe_compara = 0.0;
double SQT_compara = 0.0;

double media_compara[20];
double media_global_compara = 0.0;
for (int i = 0; i < 20; i++)
{
    media_compara[i] = 0.0;
    for (int j = 0; j < 10; j++)
         media_compara[i] += comparaQuick[i][j];
    media_global_compara += media_compara[i];
    media_compara[i] /= 10;
}
media_global_compara /= 200;
for (int i = 0; i < 20; i++)
{
    for (int j = 0; j < 10; j++)
    {
        double tmp = comparaQuick[i][j] - media_global_compara;
        SQT_compara += tmp * tmp;

        tmp = modeloNlogN(i + 1, coef_compara_quick.coeficientes[0], coef_compara_quick.coeficientes[1]);
        SQR_compara += (tmp - media_global_compara) * (tmp - media_global_compara);
        tmp = comparaQuick[i][j] - media_compara[i];
        SQe_compara += tmp * tmp;
    }
}

R2 = SQR_compara / SQT_compara;
R2max = 1 - SQe_compara / SQT_compara;
R2 *= 100;
R2max *= 100;
printf("R2 Compara Quick: %lf %%\nR2max Compara Quick: %lf %%\n", R2, R2max);


SQR_compara = 0.0;
SQe_compara = 0.0;
SQT_compara = 0.0;

media_global_compara = 0.0;
for (int i = 0; i < 20; i++)
{
    media_compara[i] = 0.0;
    for (int j = 0; j < 10; j++)
         media_compara[i] += comparaMerge[i][j];
    media_global_compara += media_compara[i];
    media_compara[i] /= 10;
}
media_global_compara /= 200;
for (int i = 0; i < 20; i++)
{
    for (int j = 0; j < 10; j++)
    {
        double tmp = comparaMerge[i][j] - media_global_compara;
        SQT_compara += tmp * tmp;

        tmp = modeloNlogN(i + 1, coef_compara_merge.coeficientes[0], coef_compara_merge.coeficientes[1]);
        SQR_compara += (tmp - media_global_compara) * (tmp - media_global_compara);
        tmp = comparaMerge[i][j] - media_compara[i];
        SQe_compara += tmp * tmp;
    }
}

R2 = SQR_compara / SQT_compara;
R2max = 1 - SQe_compara / SQT_compara;
R2 *= 100;
R2max *= 100;
printf("R2 Compara Merge: %lf %%\nR2max Compara Merge: %lf %%\n", R2, R2max);
//-----------------------------------------------------------------------------
    for (int i = 0; i < 20; i++)
    {
        free(Bubble[i]);
        free(Quick[i]);
        free(Merge[i]);
        free(comparaQuick[i]);
        free(comparaMerge[i]);
    }
    free(comparaBubble[0]);
    free(Bubble);
    free(Quick);
    free(Merge);
    free(comparaBubble);
    free(comparaQuick);
    free(comparaMerge);

    fclose(dados);
    return 0;
}
