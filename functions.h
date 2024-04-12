#include <stdio.h>
#include <stdlib.h>
#include <math.h>//

double max_norm(double *y, const int n);

void printvec(double *y, const int n);

int power10(const int n);

void obnul(double *y, int n);

double f1(double x, double y);

double f2(double x, double y);

double f3(double x);

double f4(double x);

double f5(double x);

double f6(double x);

double f7(double x);

double ortF1(double x, double y, int m, int n);

double Simpson(double a, double b, double (*f)(double));

double Gauss(double a, double b, double (*f)(double));

double SimpsonP(double a, double b, int n);

double GaussP(double a, double b, int n);

double MPow(int k);

double Norm(double **c, double **vertices, double *F, double(*fOrt)(double, double, int, int), int m, int n,
            double verticesCount);

double PowerX(double x, const int n);

double IntegralPowerX(double a, double b, int n);

double ErrorRateGaussP(double a, double b, int n);

double ErrorRateSimpsonP(double a, double b, int n);

double DividedSimpson(double a, double b, double(*f)(double), int N);

double DivededSimpsonP(double a, double b, int n, int N);

double DividedErrorRateSimpson(double a, double b, double derevativeMax, int N);

double DividedErrorRateGauss(double a, double b, double derevativeMax, int N);

double cos100(double x);

double exp1000(double x);

double sqrtdiv(double x);


double f1(double x, double y) { return 1 + x + 2 * y + x * y; }

double f2(double x, double y) { return cos(x) + cos(y); }

double f3(double x) { return fabs(x); };

double f4(double x) { return 1. / (1. + x * x); }

double f5(double x) {
    if (x > 0) return 1.;
    return -1.;
}

double cos100(double x) { return cos(100. * x); }

double exp1000(double x) { return exp(-1000. * x); }

double sqrtdiv(double x) { return 1. / sqrt(1. - x * x); }

double f6(double x) { return sin(x - 5 * exp(x)) * (x * x + sin(x)); }

double f7(double x) { return PowerX(x, 6); }

int WriteToFile(FILE *file, double *x, double *y, const int n) {
    for (int i = 0; i < n; ++i) {
        fprintf(file, "%lf %.12lf\n", x[i], y[i]);
    }
    return 1;
}

void RavUzl(double *x, const double a, const double b, const int n) {
    for (int i = 0; i < n; ++i) {
        x[i] = (b - a) / (n - 1.) * i + a;
        //printf("x[%d]=%lf   ", i, x[i]);
    }
}

double max_norm(double *y, const int n) {
    double max = 0;
    double k = 0;
    for (int i = 0; i < n; ++i) {
        if (y[i] < 0) {
            k = -y[i];
        } else {
            k = y[i];
        }
        if (k > max) { max = k; }
    }
    return max;
}


void printvec(double *y, const int n) {
    int k = 0;
    for (int i = 0; i < n / 10; ++i) {
        for (int j = 0; j < 10; ++j) {
            printf("%.12lf  ", y[k]);
            ++k;
        }
        printf("\n");
    }
    for (int i = 0; i < n % 10; ++i) {
        printf("%.12lf  ", y[k]);
        ++k;
    }
    printf("\n\n");
}

int power10(const int n) {
    int var = 1;
    for (int i = 0; i < n; ++i) var *= 10;
    return var;
}


double PowerX(double x, const int n) {
    double p;
    p = 1.;
    for (int i = 0; i < n; ++i) {
        p *= x;
    }
    return p;
}

double ErrorRateSimpsonP(double a, double b, int n) {
    double mp;
    double k = 1.;
    if (n < 4) return 0.;
    for (int i = 0; i < 4; ++i) {
        k *= (n - i);
    }
    if (fabs(a) > fabs(b)) {
        mp = a;
    } else mp = b;
    return k * PowerX(mp, n - 4) * PowerX(b - a, 5) / 2880.;
}

double ErrorRateGaussP(double a, double b, int n) {
    double mp;
    int k = 1;
    if (n < 6) return 0.;
    for (int i = 0; i < 6; ++i) {
        k *= n - i;
    }
    if (fabs(a) > fabs(b)) {
        mp = fabs(a);
    } else mp = fabs(b);
    return k * PowerX(mp, n - 6) * PowerX(b - a, 7) / 2016000.;
}

double DividedErrorRateSimpson(double a, double b, double derevativeMax, int N) {
    return derevativeMax * pow(b - a, 5.) / (2880. * pow(N, 4.));
}

double DividedErrorRateGauss(double a, double b, double derevativeMax, int N) {
    return derevativeMax * pow(b - a, 7.) / 2016000. / pow(N, 6);
}

double IntegralPowerX(double a, double b, int n) {
    return 1. / (n + 1) * (PowerX(b, n + 1) - PowerX(a, n + 1));
}

double MPow(int k) {
    if (k % 2 == 1) return -1.;
    return 1.;
}


void Find2DMap(double **vertices, double *F, double (*f)(double, double), int verticesCount) {
    for (int i = 0; i < verticesCount; ++i) {
        F[i] = f(vertices[i][0], vertices[i][1]);
        //printf("i: %d, x: %lf, y: %lf\n", i, vertices[i][0], vertices[i][1]);
        //printf("F: %lf\n", F[i]);
    }
}

double ortF1(double x, double y, int m, int n) {
    return sin(M_PI * (m + 1) * x) * cos(M_PI * (n + 1) * y);
}

double ortF2(double x, double y, int m, int n) {
    return pow(x, m) * pow(y, n);
}

/*
void FindGrad(double x, double y, double **grad, double **c, double **vertices, double *F,
              double(*Norm)(double **, double **, double *, double(*)(double, double, int, int), int, int,
                            double), double(*fOrt)(double, double, int, int), int m, int n, int verticesCount) {
    int maxiter = 100;
    double delta = 0.1;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            grad[i][j] = Norm(c, vertices, F, )
        }
    }
}
*/

void
FindCVar(double **c, double **vertices, double *F, double(*fOrt)(double, double, int, int), int m, int n,
         int verticesCount) {
    double ch = 0., zn = 0.;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            ch = 0.;
            zn = 0.;
            for (int k = 0; k < verticesCount; ++k) {
                ch += fOrt(vertices[k][0], vertices[k][1], i, j) * F[k];
                //printf("g: %lf, %lf\n", fOrt(vertices[k][0], vertices[k][1], i, j), F[k]);
                zn += fOrt(vertices[k][0], vertices[k][1], i, j) * fOrt(vertices[k][0], vertices[k][1], i, j);
                //printf("x: %lf, y: %lf\n", vertices[k][0], vertices[k][1]);
            }

            printf("ch: %.1e, zn: %.1e\n", ch, zn);
            if (zn != 0.) {
                c[i][j] = ch / zn;
            } else
                c[i][j] = 0.;
            printf("c: %lf\n", c[i][j]);
        }
    }
}

double QuasiF(double x, double y, double **c, double(*fOrt)(double, double, int, int), int m, int n) {
    double sum = 0;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            sum += c[i][j] * fOrt(x, y, i, j);
        }
    }
    return sum;
}

double Norm(double **c, double **vertices, double *F, double(*fOrt)(double, double, int, int), int m, int n,
            double verticesCount) {
    double sum = 0;
    double summn;
    for (int i = 0; i < verticesCount; ++i) {
        summn = 0;
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < n; ++k) {
                summn += c[j][k] * fOrt(vertices[i][0], vertices[i][1], j, k);
                //printf("c: %lf\n", c[j][k]);
            }
        }
        //printf("vert: %lf %lf\n", vertices[i][0], vertices[i][1]);
        //printf("%lf\n", sum);
        sum += (F[i] - summn) * (F[i] - summn);
        //printf("Phi: %lf, F: %lf\n", summn, F[i]);
    }
    if (sum < 2.2 * 1e-16) sum = 0;
    return sqrt(sum);
}

void NormalizeVector(double **v, int m, int n) {
    double norm = 0;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            //norm += v[i][j] * v[i][j];
            if (norm < fabs(v[i][j])) norm = fabs(v[i][j]);
        }
    }
    //printf("норма: %lf\n", norm);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            v[i][j] = v[i][j] / norm;
            //printf("v: %lf\n", v[i][j]);
        }
    }
}

void GaussSolve(double **a, double *x, double *y, int n) {
    double max, temp;
    int k, index;
    const double eps = 3e-16;  // точность
    k = 0;
    while (k < n) {
        // Поиск строки с максимальным a[i][k]
        max = fabs(a[k][k]);
        index = k;
        for (int i = k + 1; i < n; ++i) {
            if (fabs(a[i][k]) > max) {
                max = fabs(a[i][k]);
                index = i;
            }
        }
        // Перестановка строк
        if (max < eps) {
            // нет ненулевых диагональных элементов
            printf("Решение получить невозможно из-за нулевого столбца ");
        }
        for (int j = 0; j < n; ++j) {
            temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;
        }
        temp = y[k];
        y[k] = y[index];
        y[index] = temp;
        // Нормализация уравнений
        for (int i = k; i < n; ++i) {
            temp = a[i][k];
            if (fabs(temp) < eps) continue; // для нулевого коэффициента пропустить
            for (int j = k; j < n; j++)
                a[i][j] = a[i][j] / temp;
            y[i] = y[i] / temp;
            if (i == k) continue; // уравнение не вычитать само из себя
            for (int j = 0; j < n; ++j)
                a[i][j] = a[i][j] - a[k][j];
            y[i] = y[i] - y[k];
        }
        k++;
    }
    // обратная подстановка
    for (k = n - 1; k >= 0; --k) {
        x[k] = y[k];
        for (int i = 0; i < k; i++)
            y[i] = y[i] - a[i][k] * x[k];
    }
}

void PrintMatrix(double **a, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%lf  ", a[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
}

double GaussSolve1(double **a, double *x, double *b, int n) { //a[i][j] i строка, j столбец
    double max, tmp;
    int iRMax;
    for (int i = 0; i < n; ++i) { //по строкам
        max = 0;
        for (int j = i; j < n; ++j) { //по строкам находим строку с максимальным элементом
            if (fabs(a[j][i]) > max) {
                max = fabs(a[j][i]);
                iRMax = j;
            }
        }
        //printf("iRmax: %d\n", iRMax);
        //PrintMatrix(a, n);

        if (i != iRMax) { //Переставляем строки
            for (int j = 0; j < n; ++j) {
                tmp = a[i][j];
                a[i][j] = a[iRMax][j];
                a[iRMax][j] = tmp;
            }
            tmp = b[i];
            b[i] = b[iRMax];
            b[iRMax] = tmp;
        }

        //printf("Переставили строки\n");
        //PrintMatrix(a, n);

        for (int j = i; j < n; ++j) { //приведение к еденице диаганального элемента
            a[i][j] /= a[i][i];
        }
        b[i] /= a[i][i];

        //printf("Привели к единице диаганльный элемент\n");
        //PrintMatrix(a, n);

        for (int k = i + 1; k < n; ++k) { //Вычитание из строк
            for (int j = 0; j < n; ++j) { // Лучше взять j = i
                a[k][j] -= a[i][j] * a[k][i];
            }
            b[k] -= b[i] * a[k][i];
        }
        //printf("Пошли нахуй\n");
        //PrintMatrix(a, n);
    }

    for (int i = n - 1; i >= 0; --i) {
        for (int k = i - 1; k >= 0; --k) {
            a[k][i] -= a[i][i] * a[k][i];
            b[k] -= b[i] * a[k][i];
        }
    }
    for(int i = 0; i < n; ++i){
        x[i] = b[i];
    }
    //printf("Пошли нахуй\n");
    //PrintMatrix(a, n);
}


void MakeMatrix(double **A, double **v, double(*fOrt)(double, double, int, int), int m, int n, int verticesCount) {
    double sum = 0;
    for (int i = 0; i < m * n; ++i) {
        for (int j = 0; j < m * n; ++j) {
            sum = 0;
            for (int k = 0; k < verticesCount; ++k) {
                sum += fOrt(v[k][0], v[k][1], i / m, i % m) * fOrt(v[k][0], v[k][1], j / m, j % m);
            }
            A[i][j] = sum;
        }
    }
}

void
MakeY(double *y, double **v, double *F, double (*fOrt)(double, double, int, int), int m, int n, int verticesCount) {
    double sum = 0;
    for (int i = 0; i < m * n; ++i) {
        sum = 0;
        for (int j = 0; j < verticesCount; ++j) {
            sum += F[j] * fOrt(v[j][0], v[j][1], i / m, i % m);
            //printf("f* phi: %lf", F[j] * fOrt(v[j][0], v[j][1], i / m, i % m));
            //printf("sum: %lf\n", sum);
        }
        //printf("\n");
        y[i] = sum;
    }
}