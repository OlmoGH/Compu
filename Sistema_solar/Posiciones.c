# include "funciones.h"
# include <stdlib.h>
# include <math.h>
# include <stdio.h>
# define PI 3.14159265
# define CUERPOS 1 //Cuantos cuerpos, planetas, a parte del sol estamos tomando

typedef double Vector2D[2];
int main(){
    //Guardamos espacio en memoria para las masas, las posiciones, velocidades y aceleraciones
    double *masas = (double *) malloc(CUERPOS * sizeof(double));
    Vector2D *posiciones = (Vector2D *) malloc(CUERPOS * 2 * sizeof(double)); 
    Vector2D *velocidades = (Vector2D *) malloc(CUERPOS * 2 * sizeof(double));
    Vector2D *aceleraciones = (Vector2D *) malloc(CUERPOS * 2 * sizeof(double));
    double posiciones_iniciales[CUERPOS][2] = {{1.496E11, 0}};
    double velocidades_iniciales[CUERPOS][2] = {{0, 0}};
    double excentricidades[CUERPOS] = {0.01671};

    for (int i = 0; i < CUERPOS; i++)
    {
        reescalar_posiciones(&posiciones_iniciales[i], CUERPOS);
    }
    
    //Conociendo el valor de las excentricidades y la distancia respecto al sol calculamos
    //la velocidad inicial
    for (int i = 0; i < CUERPOS; i++)
    {
        velocidades_iniciales[i][0] =  0;
        velocidades_iniciales[i][1] = sqrt((1 - excentricidades[i] * excentricidades[i])/(2 * posiciones_iniciales[i][0]));
    }
    
    printf("Magnitud de la velocidad inicial de la Tierra %f\n", velocidades_iniciales[0][1]);
    free(masas);
    free(posiciones);
    free(velocidades);
    free(aceleraciones);
    return 0;
}
