# include "funciones.h"
# include <stdlib.h>

typedef double Vector2D[2];
int main(){
    const int cuerpos = 2;
    double *masas;
    Vector2D *posiciones = (Vector2D *) malloc(cuerpos * sizeof(Vector2D)); 
    Vector2D *velocidades = (Vector2D *) malloc(cuerpos * sizeof(Vector2D));
    Vector2D *aceleraciones = (Vector2D *) malloc(cuerpos * sizeof(Vector2D));
    return 0;
}
