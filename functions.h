#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void PrintVec(double *y, int n);

double f1(double x, double y);

double f2(double x, double y);

double f3(double x, double y);

double f4(double x, double y);

double f5(double x, double y);

double f6(double x, double y);

double ortF1(double x, double y, int m, int n);

double ortF2(double x, double y, int m, int n);

double Norm(double **c, double **vertices, double *F, double(*fOrt)(double, double, int, int), int m, int n,
            double verticesCount);

double NormMax(double **c, double **vertices, double *F, double(*fOrt)(double, double, int, int), int m, int n,
               double verticesCount);

void RavUzl(double *x, double a, double b, int n);

void RandomUzl(double *x, double a, double b, int n);

void RandRemoveVertices(double **v, double *F, int percent, int *count);

double PolCheb(double x, int n);

void PrintVertices(double **v, int n);

void Find2DMap(double **vertices, double *F, double (*f)(double, double), int verticesCount);

void PrintMatrixWithb(double **a, double *b, int n);

void PrintMatrix(double **a, int n);

void
FindCVar(double **c, double **vertices, double *F, double(*fOrt)(double, double, int, int), int m, int n,
         int verticesCount);

double QuasiF(double x, double y, double **c, double(*fOrt)(double, double, int, int), int m, int n);

double GaussSolve(double **a, double *x, double *b, int n);

void MakeMatrix(double **A, double **v, double(*fOrt)(double, double, int, int), int m, int n, int verticesCount);

void
MakeY(double *y, double **v, double *F, double (*fOrt)(double, double, int, int), int m, int n, int verticesCount);

double f1(double x, double y) { return 1 + x + 2 * y + x * y; }

double f2(double x, double y) { return sin(2 * M_PI * x) * cos(2 * M_PI * y); }

double f3(double x, double y) { return cos(x); }

double f4(double x, double y) { return sin(x); }

double f5(double x, double y) { return (exp(x) - 1) * (y + 1); }

double f6(double x, double y) { return y * x * sin(5 * x + y) / exp(8 * x * y); }

void RavUzl(double *x, double a, double b, int n) {
    for (int i = 0; i < n; ++i) {
        x[i] = (b - a) / (n - 1.) * i + a;
        //printf("x[%d]=%lf   ", i, x[i]);
    }
}

void RandomUzl(double *x, double a, double b, int n) {
    for (int i = 0; i < n; ++i) {
        x[i] = rand() / (double) RAND_MAX * (b - a) + a;
    }
}

void RandRemoveVertices(double **v, double *F, int percent, int *count) {
    int r;
    int verticesCount = *count;
    int deleteCount = verticesCount * percent / 100;
    for (int i = 0; i < deleteCount; ++i) {
        r = rand() % verticesCount;
        //printf("r: %d\n", r);
        for (int j = r; j < verticesCount - 1; ++j) {
            v[j][0] = v[j + 1][0];
            v[j][1] = v[j + 1][1];
            F[j] = F[j + 1];
        }
        --verticesCount;
    }
    *count = verticesCount;
}

void PrintVertices(double **v, int n) {
    int k = 0;
    for (int i = 0; i < n / 10; ++i) {
        for (int j = 0; j < 10; ++j) {
            printf("(%lf, %lf)  ", v[k][0], v[k][1]);
            ++k;
        }
        printf("\n");
    }
    for (int i = 0; i < n % 10; ++i) {
        printf("(%lf, %lf)  ", v[k][0], v[k][1]);
        ++k;
    }
    printf("\n\n");
}

double PolCheb(double x, int n) {
    if (n == 0) return 1;
    if (n == 1) return x;
    return 2 * x * PolCheb(x, n - 1) - PolCheb(x, n - 2);
}

void PrintVec(double *y, int n) {
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
    return PolCheb(x, m) * PolCheb(y, n);
}

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

double NormMax(double **c, double **vertices, double *F, double(*fOrt)(double, double, int, int), int m, int n,
               double verticesCount) {
    double max = 0;
    double summn = 0;
    double xMax, yMax, fMax, phiMax;
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
        if (fabs(F[i] - summn) > max) {
            max = fabs(F[i] - summn);
            fMax = F[i];
            phiMax = summn;
            xMax = vertices[i][0];
            yMax = vertices[i][1];
        }
        //printf("Phi: %lf, F: %lf\n", summn, F[i]);
    }
    printf("Максимальное отклонение в точке (%lf, %lf)\n", xMax, yMax);
    printf("F: %lf   Phi: %lf\n", fMax, phiMax);
    return max;
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

void PrintMatrixWithb(double **a, double *b, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%lf  ", a[i][j]);
        }
        printf("|%lf  ", b[i]);
        printf("\n");
    }
    printf("\n\n");
}

double GaussSolve(double **a, double *x, double *b, int n) { //a[i][j] i строка, j столбец
    double max, tmp;
    int iRMax;

    //PrintVec(b, n);
    for (int i = 0; i < n; ++i) { //по строкам
        max = 0;
        for (int j = i; j < n; ++j) { //по строкам находим строку с максимальным элементом
            if (fabs(a[j][i]) > max) {
                max = fabs(a[j][i]);
                iRMax = j;
            }
        }
        //printf("iRmax: %d\n", iRMax);
        //PrintMatrixWithb(a, b, n);

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
        //PrintMatrixWithb(a, b, n);

        for (int j = i + 1; j < n; ++j) { //приведение к еденице диаганального элемента
            a[i][j] /= a[i][i];
        }
        b[i] /= a[i][i];
        a[i][i] = 1;

        //printf("Привели к единице диаганльный элемент\n");
        //PrintMatrixWithb(a, b, n);

        for (int k = i + 1; k < n; ++k) { //Вычитание из строк
            tmp = a[k][i];
            b[k] -= b[i] * tmp;
            for (int j = 0; j < n; ++j) { // Лучше взять j = i
                a[k][j] -= a[i][j] * tmp;
            }
        }
        //printf("Пошли нахуй\n");
        //PrintMatrixWithb(a, b, n);
    }
    //PrintVec(b, n);
    for (int i = n - 1; i >= 0; --i) {
        for (int k = i - 1; k >= 0; --k) {
            b[k] -= b[i] * a[k][i];
            a[k][i] -= a[i][i] * a[k][i];
        }
        //printf("Сильно идём нахуй\n");
        //PrintMatrixWithb(a, b, n);
    }
    for (int i = 0; i < n; ++i) {
        x[i] = b[i];
    }
    //PrintVec(b, n);
    //printf("Пошли нахуй\n");
    //PrintMatrixWithb(a, b, n);
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