#include "functions.h"

int main() {
    FILE *file;
    double **vertices, **c; //m, n
    double *xUzl, *yUzl, *F;
    double x_0, x_1, y_0, y_1;
    int xSplitCount, ySplitCount, verticesCount;
    int m, n;
    int outSplitCount;
    double deviation;

    //Объявление переменных
    m = 5;
    n = 6;
    xSplitCount = 20;
    ySplitCount = 21;
    x_0 = 0;
    x_1 = 1;
    y_0 = 0;
    y_1 = 1;
    outSplitCount = 100;
    file = fopen("data.txt", "w");
    verticesCount = xSplitCount * ySplitCount;

    //Выделение памяти
    xUzl = (double *) malloc((xSplitCount) * sizeof(double));
    yUzl = (double *) malloc((ySplitCount) * sizeof(double));
    c = (double **) malloc(m * sizeof(double *));
    vertices = (double **) malloc(verticesCount * sizeof(double *));
    F = (double *) malloc((xSplitCount * ySplitCount) * sizeof(double));

    //Выделение памяти двумерным массивам
    for (int i = 0; i < m; ++i) {
        c[i] = (double *) malloc(n * sizeof(double));
    }
    for (int i = 0; i < xSplitCount * ySplitCount; ++i) {
        vertices[i] = (double *) malloc(2 * sizeof(double));
    }

    //Разбиение плоскости на квадраты
    RavUzl(xUzl, x_0, x_1, xSplitCount);
    RavUzl(yUzl, y_0, y_1, ySplitCount);
    for (int i = 0; i < xSplitCount; ++i) {
        for (int j = 0; j < ySplitCount; ++j) {
            vertices[i * xSplitCount + j][0] = xUzl[i];
            vertices[i * xSplitCount + j][1] = yUzl[j];
        }
    }

    //Нахождение Fij, Cmn
    Find2DMap(vertices, F, f1, verticesCount);
    FindCVar(c, vertices, F, ortF1, m ,n, verticesCount);

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


    deviation = Norm(c, vertices, F, ortF1, m, n, verticesCount);
    printf("Погрешность:\n%lf\n%.1e\n", deviation, deviation);

    //Очистка памяти
    free(vertices);
    free(c);
    free(xUzl);
    free(yUzl);
    free(F);
    fclose(file);
    return 0;
}
