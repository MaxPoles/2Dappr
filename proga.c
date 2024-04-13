#include "functions.h"


/*
 c - двумерный массив клэффециентов по m, n
 cLine - одномерный массив коэффециентов
 T - кол-во функций
 */

int main() {
    FILE *fileX, *fileY, *fileZ;
    double **vertices, **c, **grad, **A; //m, n
    double *xUzl, *yUzl, *F, *y, *cLine;
    double (*fOrt)(double, double, int, int);
    double (*f)(double, double);
    void (*UzlDivide)(double *, double, double, int);
    double x_0, x_1, y_0, y_1;
    int percent;
    int xSplitCount, ySplitCount, verticesCount;
    int m, n, kVar;
    int outSplitCount;
    double deviation;
    //m 12,n 12,x 14,y 14, f6, 15
    //Объявление переменных
    m = 10;
    n = 10;
    xSplitCount = 15;
    ySplitCount = 15;
    x_0 = -1;
    x_1 = 1;
    y_0 = -1;
    y_1 = 1;
    outSplitCount = 100;
    percent = 0;

    fOrt = ortF2;
    f = f7;
    UzlDivide = RandomUzl;

    fileX = fopen("x.txt", "w");
    fileY = fopen("y.txt", "w");
    fileZ = fopen("z.txt", "w");
    verticesCount = xSplitCount * ySplitCount;
    printf("Кол-во точек: %d\n", verticesCount);
    //Нормировка
    ++m;
    ++n;

    //Выделение памяти
    srand(time(NULL));

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
    UzlDivide(xUzl, x_0, x_1, xSplitCount);
    //PrintVec(xUzl, xSplitCount);
    UzlDivide(yUzl, y_0, y_1, ySplitCount);
    //PrintVec(yUzl, ySplitCount);
    kVar = 0;
    for (int i = 0; i < xSplitCount; ++i) {
        for (int j = 0; j < ySplitCount; ++j) {
            vertices[kVar][0] = xUzl[i];
            vertices[kVar][1] = yUzl[j];
            ++kVar;
            //printf("i: %d, x: %lf, y: %lf\n", kVar - 1, vertices[kVar - 1][0], vertices[kVar - 1][1]);
        }
    } // Тут всё зашибись, не проверяй

    //Строим F и выкидываем часть точек
    Find2DMap(vertices, F, f, verticesCount);
    RandRemoveVertices(vertices, F, percent, &verticesCount);

    //Создаём матрицу и вектора
    MakeY(y, vertices, F, fOrt, m, n, verticesCount);
    MakeMatrix(A, vertices, fOrt, m, n, verticesCount);

    //PrintVertices(vertices, verticesCount);
    //PrintVertices(vertices, verticesCount);
    //printf("verticesCount: %d\n", verticesCount);
    printf("Кол-во точек после выброса: %d\n", verticesCount);


    //Решаем СЛУ
    GaussSolve(A, cLine, y, m * n);
    for (int i = 0; i < m * n; ++i) {
        c[i / n][i % n] = cLine[i];
    }

    //Выводим матрицу коэффециентов c
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%lf  ", c[i][j]);
        }
        printf("\n");
    }
    deviation = Norm(c, vertices, F, fOrt, m, n, verticesCount);
    printf("Норма:\n%lf\n%.1e\n", deviation, deviation);
    deviation = NormMax(c, vertices, F, fOrt, m, n, verticesCount);
    printf("Максимальное отклонение:\n%lf\n%.1e\n", deviation, deviation);
    //SaveCToFile(file, c, m, n);

    free(xUzl);
    free(yUzl);
    xUzl = (double *) malloc((outSplitCount) * sizeof(double));
    yUzl = (double *) malloc((outSplitCount) * sizeof(double));

    RavUzl(xUzl, x_0, x_1, outSplitCount);
    RavUzl(yUzl, y_0, y_1, outSplitCount);

    for (int i = 0; i < outSplitCount; ++i) {
        for (int j = 0; j < outSplitCount; ++j) {
            fprintf(fileX, "%.12lf ", xUzl[j]);
            fprintf(fileY, "%.12lf ", yUzl[i]);

            fprintf(fileZ, "%.12lf ", QuasiF(xUzl[j], yUzl[i], c, fOrt, m, n));
        }
        fprintf(fileX, "\n");
        fprintf(fileY, "\n");
        fprintf(fileZ, "\n");
    }


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
    fclose(fileX);
    fclose(fileY);
    fclose(fileZ);
    return 0;
}
