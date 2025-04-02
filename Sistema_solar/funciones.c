# include "funciones.h"
# include <stdlib.h>
# include <math.h>

typedef double Vector2D[2];

Vector2D *reescalar_posiciones(Vector2D *posiciones, int cuerpos){
    if (posiciones == NULL || cuerpos <= 0) {
        return NULL;
    }

    int i, j;
    const double C = 1.496E11;

    for(i = 0; i < cuerpos; i++){
        for (j = 0; j < 2; j++)
        {
            posiciones[i][j] = posiciones[i][j] / C;
        }
        
    }
    return posiciones;
}

double reescalar_tiempo(double tiempo){
    const double C = 1.496E11;
    const double Ms = 1.989E30;
    const double G = 6.672E-11;

    tiempo = tiempo * sqrt(G * Ms / pow(C, 3));
    return tiempo;
}

double *reescalar_masas(double *masas, int size){
    if(masas == NULL || size <= 0){
        return NULL;
    }
    const double Ms = 1.989E30;

    for (int i = 0; i < size; i++)
    {
        masas[i] = masas[i] / Ms;
    }
    
    return masas;
}

Vector2D *reescalar_posiciones_inv(Vector2D *posiciones, int cuerpos){
    if (posiciones == NULL || cuerpos <= 0) {
        return NULL;
    }

    int i, j;
    const double C = 1.496E11;

    for(i = 0; i < cuerpos; i++){
        for (j = 0; j < 2; j++)
        {
            posiciones[i][j] = posiciones[i][j] * C;
        }
        
    }
    return posiciones;
}

double reescalar_tiempo_inv(double tiempo){
    const double C = 1.496E11;
    const double Ms = 1.989E30;
    const double G = 6.672E-11;

    tiempo = tiempo / sqrt(G * Ms / pow(C, 3));
    return tiempo;
}

double *reescalar_masas_inv(double *masas, int size){
    if(masas == NULL || size <= 0){
        return NULL;
    }
    const double Ms = 1.989E30;

    for (int i = 0; i < size; i++)
    {
        masas[i] = masas[i] * Ms;
    }
    
    return masas;
}