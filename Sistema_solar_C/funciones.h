# ifndef FUNCIONES_H
# define FUNCIONES_H

//Declaración de funciones
typedef double Vector2D[2];

double norma_cuad(Vector2D posicion);

double momento_angular(Vector2D posicion, double masa, Vector2D velocidad);

double energia_mecanica(Vector2D posicion, double masa, Vector2D velocidad, double masa_origen);

double *reescalar_distancias(double *r, int cuerpos);

Vector2D *reescalar_posiciones(Vector2D *posiciones, int cuerpos);

double reescalar_tiempo(double tiempo);

double *reescalar_masas(double *masas, int size);

double *reescalar_distancias_inv(double *r, int cuerpos);

Vector2D *reescalar_posiciones_inv(Vector2D *posiciones, int cuerpos);

double reescalar_tiempo_inv(double tiempo);

double *reescalar_masas_inv(double *masas, int size);

# endif