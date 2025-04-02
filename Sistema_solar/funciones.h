# ifndef FUNCIONES_H
# define FUNCIONES_H

//Declaración de funciones
typedef double Vector2D[2];

Vector2D *reescalar_posiciones_inv(Vector2D *posiciones, int cuerpos);

double reescalar_tiempo(double tiempo);

double *reescalar_masas(double *masas, int size);

Vector2D *reescalar_posiciones_inv(Vector2D *posiciones, int cuerpos);

double reescalar_tiempo_inv(double tiempo);

double *reescalar_masas_inv(double *masas, int size);

# endif