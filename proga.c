#include "functions.h"


/*
 c - двумерный массив клэффециентов по m, n
 cLine - одномерный массив коэффециентов
 T - кол-во функций
 */

int main() {
    FILE *file;
    double **vertices, **c, **grad, **A; //m, n
    double *xUzl, *yUzl, *F, *y, *cLine;
    double (*fOrt)(double, double, int, int);
    double x_0, x_1, y_0, y_1;
    int xSplitCount, ySplitCount, verticesCount;
    int m, n, T, kVar;
    int outSplitCount;
    double deviation;
    double delta = 1e-8, h = 1e-4, stepDelta = 1e-3;
    int maxIter = 1000000;

    //Объявление переменных
    m = 11;
    n = 12;
    xSplitCount = 12;
    ySplitCount = 10;
    x_0 = 1;
    x_1 = 5.2;
    y_0 = 1;
    y_1 = 5;
    T = m * n;
    outSplitCount = 100;
    fOrt = ortF1;
    file = fopen("data.txt", "w");
    verticesCount = xSplitCount * ySplitCount;
    printf("verticesCount: %d\n", verticesCount);
    //Нормировка
    ++m;
    ++n;

    //Выделение памяти
    xUzl = (double *) malloc((xSplitCount) * sizeof(double));
    yUzl = (double *) malloc((ySplitCount) * sizeof(double));
    cLine = (double *) malloc(m * n * sizeof(double));
    y = (double *) malloc(m * n * sizeof(double));
    F = (double *) malloc((xSplitCount * ySplitCount) * sizeof(double));

    c = (double **) malloc(m * sizeof(double *));
    A = (double **) malloc(m * n * sizeof(double *));
    grad = (double **) malloc(m * sizeof(double *));
    vertices = (double **) malloc(verticesCount * sizeof(double *));

    //Выделение памяти двумерным массивам
    for (int i = 0; i < m; ++i) {
        c[i] = (double *) malloc(n * sizeof(double));
        grad[i] = (double *) malloc(n * sizeof(double));
    }
    for (int i = 0; i < xSplitCount * ySplitCount; ++i) {
        vertices[i] = (double *) malloc(2 * sizeof(double));
    }
    for (int i = 0; i < m * n; ++i) {
        A[i] = (double *) malloc(m * n * sizeof(double));
    }
    //Разбиение плоскости на квадраты
    RavUzl(xUzl, x_0, x_1, xSplitCount);
    //printvec(xUzl, xSplitCount);
    RavUzl(yUzl, y_0, y_1, ySplitCount);
    //printvec(yUzl, ySplitCount);
    kVar = 0;
    for (int i = 0; i < xSplitCount; ++i) {
        for (int j = 0; j < ySplitCount; ++j) {
            vertices[kVar][0] = xUzl[i];
            vertices[kVar][1] = yUzl[j];
            ++kVar;
            //printf("i: %d, x: %lf, y: %lf\n", kVar - 1, vertices[kVar - 1][0], vertices[kVar - 1][1]);
        }
    } // Тут всё зашибись, не проверяй
    Find2DMap(vertices, F, f1, verticesCount);
    MakeY(y, vertices, F, fOrt, m, n, verticesCount);

    MakeMatrix(A, vertices, fOrt, m, n, verticesCount);
    //printf("Выберете метод:\n1 - Градиентный спуск\n2 - Решение СЛУ\n3 - Градиентный спуск V2\n");
    //scanf("%d", &kVar);

    GaussSolve1(A, cLine, y, m * n);
    for (int i = 0; i < m * n; ++i) {
        //printf("Cline: ");
        //printvec(cLine, m * n);
        c[i / n][i % n] = cLine[i];
    }
    /*
    } else if(kVar == 1) {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                c[i][j] = 0;
            }
        }
        double norm1 = 0;
        double norm2 = 0;
        double lastNorm = 0;
        for (int t = 0; t < maxIter; ++t) {
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    c[i][j] -= delta;
                    norm1 = Norm(c, vertices, F, ortF2, m, n, verticesCount);
                    c[i][j] += 2 * delta;
                    norm2 = Norm(c, vertices, F, ortF2, m, n, verticesCount);
                    grad[i][j] = (norm2 - norm1) / delta / 2.;
                    c[i][j] -= delta;
                    if (t % 40000 == 0)
                        printf("%lf\n", grad[i][j]);
                }
            }
            //norm1 = Norm(c, vertices, F, ortF2, m, n, verticesCount); // счиитаем норму на данном моменте
            if (t % 40000 == 0) {
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < m; ++j) {
                        printf("%lf  ", c[i][j]);
                    }
                    printf("\n");
                }
                printf("h: %.1e\n", h);
                printf("|| ||: %lf   %.1e\n", norm1, norm1);
                printf("\n");
            }
            double max = 0;
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (max < fabs(grad[i][j])) max = fabs(grad[i][j]);
                }
            }
            if (max < 1e-12) break;
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    grad[i][j] = grad[i][j] / max;
                }
            }

            //NormalizeVector(grad, m, n);
            //for (int i = 0; i < m; ++i) {
            //    for (int j = 0; j < n; ++j) {
            //        c[i][j] -= grad[i][j] * h;
            //    }
            //}

            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    c[i][j] -= grad[i][j] * h;
                }
            }

            if (max > lastNorm) {
                h /= 2.;
            } else {
                h *= 2.;
            }


            if (fabs(lastNorm - max) < stepDelta) {
                h *= 2.;
            } else if (max > lastNorm && h > 1e-12) {
                h /= 2.;
            }

            lastNorm = max;
            //printf("\n");
        }
    }
    */
    //Нахождение градиента
    /*
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            c[i][j] = 0;
        }
    }
    double norm1 = 0;
    double norm2 = 0;
    double lastNorm = 0;
    for (int t = 0; t < maxIter; ++t) {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                c[i][j] -= delta;
                norm1 = Norm(c, vertices, F, ortF2, m, n, verticesCount);
                c[i][j] += 2 * delta;
                norm2 = Norm(c, vertices, F, ortF2, m, n, verticesCount);
                grad[i][j] = (norm2 - norm1) / delta / 2.;
                c[i][j] -= delta;
                if (t % 40000 == 0)
                    printf("%lf\n", grad[i][j]);
            }
        }
        //norm1 = Norm(c, vertices, F, ortF2, m, n, verticesCount); // счиитаем норму на данном моменте
        if (t % 40000 == 0) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < m; ++j) {
                    printf("%lf  ", c[i][j]);
                }
                printf("\n");
            }
            printf("h: %.1e\n", h);
            printf("|| ||: %lf\n", norm1);
            printf("\n");
        }
        double max = 0;
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                if (max < fabs(grad[i][j])) max = fabs(grad[i][j]);
            }
        }
        if (max < 1e-12) break;
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                grad[i][j] = grad[i][j] / max;
            }
        }

        //NormalizeVector(grad, m, n);
        //for (int i = 0; i < m; ++i) {
        //    for (int j = 0; j < n; ++j) {
        //        c[i][j] -= grad[i][j] * h;
        //    }
        //}

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            c[i][j] -= grad[i][j] * h;
        }
    }

    //if (max > lastNorm) h /= 2.;
    if (fabs(lastNorm - max) < stepDelta) {
        h *= 2.;
    } else if (max > lastNorm && h > 1e-12) {
        h /= 2.;
    }
    lastNorm = max;
    //printf("\n");
}    grad = (double **) malloc(m * sizeof(double *));
    vertices = (double **) malloc(verticesCount * sizeof(double *));

*/

//Нахождение Cmn
//
//FindCVar(c, vertices, F, ortF2, m ,n, verticesCount);

//Вывод данных в файл
/*
for(int i = 0; i < m; ++i){
    for(int j = 0; i < n; ++j){
        double xDot = (x_1 - x_0) / (outSplitCount - 1.) * i + x_0;
        double yDot = (y_1 - y_0) / (outSplitCount - 1.) * j + y_0;
        fprintf()
    }
}
*/


    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%lf  ", c[i][j]);
        }
        printf("\n");
    }
    deviation = Norm(c, vertices, F, fOrt, m, n, verticesCount);
/*
for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
        printf("%lf  ", c[i][j]);
    }
    printf("\n");
}

    for (int i = 0; i < verticesCount; ++i) {
        printf("F %lf\n", F[i]);
    }
    */
    printf("Погрешность:\n%lf\n%.1e\n", deviation, deviation);

//Очистка памяти
    free(vertices);
    free(c);
    free(y);
    free(A);
    free(cLine);
    free(xUzl);
    free(yUzl);
    free(F);
    free(grad);
    fclose(file);
    return 0;
}
