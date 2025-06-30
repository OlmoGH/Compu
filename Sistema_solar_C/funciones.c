# include "funciones.h"
# include <stdlib.h>
# include <stdio.h>
# include <math.h>

typedef double Vector2D[2];

double norma_cuad(Vector2D posicion){
    double d = posicion[0] * posicion[0] + posicion[1] * posicion[1];
    return d;
}

double momento_angular(Vector2D posicion, double masa, Vector2D velocidad){
    return masa * (posicion[0] * velocidad[1] - posicion[1] * velocidad[0]);
}

double energia_mecanica(Vector2D posicion, double masa, Vector2D velocidad, double masa_origen){
    return 0.5 * masa * norma_cuad(velocidad) - masa * masa_origen / sqrt(norma_cuad(posicion));
}

double *reescalar_distancias(double *r, int cuerpos){
    if(r == NULL || cuerpos <= 0){
        return NULL;
    }
    const double C = 1.496E11;
    
    for (int i = 0; i < cuerpos; i++)
    {
        r[i] = r[i] / C;
    }
    
    return r;
}

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

double *reescalar_distancias_inv(double *r, int cuerpos){
    if(r == NULL || cuerpos <= 0){
        return NULL;
    }
    const double C = 1.496E11;
    
    for (int i = 0; i < cuerpos; i++)
    {
        r[i] = r[i] * C;
    }
    
    return r;
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