#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

#define FILAS 500
#define COLUMNAS 500

void multiplicarMatrizPorDos(float matriz[FILAS][COLUMNAS], int filas, int columnas, float resultado[FILAS][COLUMNAS]) {
    for (int i = 0; i < filas; i++) {
        for (int j = 0; j < columnas; j++) {
            resultado[i][j] = matriz[i][j] * 2;
        }
    }
}

int main() {
    clock_t inicio, fin;
    double tiempo;
    float matriz[FILAS][COLUMNAS];
    float resultado[FILAS][COLUMNAS];


    //Inicializamos la matriz original
    for(int i = 0; i<FILAS; i++){
        for(int j = 0; j<COLUMNAS; j++){
            matriz[i][j] = i*j;
        }
    }
    /*
    printf("Matriz original:\n");
    for (int i = 0; i < FILAS; i++) {
        for (int j = 0; j < COLUMNAS; j++) {
            printf("%f ", matriz[i][j]);
        }
        printf("\n");
    }
    */
    //Multiplicamos la matriz por dos y guardamos el resultado en otra matriz llamada resultado
    //Sin paralelizar
    inicio = clock();
    multiplicarMatrizPorDos(matriz, FILAS, COLUMNAS, resultado);
    fin = clock();
    tiempo = ((double)(fin - inicio))/CLOCKS_PER_SEC;
    
    /*
    printf("\nMatriz multiplicada por 2:\n");
    for (int i = 0; i < FILAS; i++) {
        for (int j = 0; j < COLUMNAS; j++) {
            printf("%f ", resultado[i][j]);
        }
        printf("\n");
    }
    */
    printf("\nEl tiempo que tarda en ejecutarse sin paralelizar es %f s\n", tiempo);
    return 0;
}