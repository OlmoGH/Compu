# include "funciones.h"
# include <stdlib.h>
# include <math.h>
# include <stdio.h>
# define PI 3.14159265
# define CUERPOS 2 //Cuantos cuerpos, planetas, a parte del sol estamos tomando

typedef double Vector2D[2];
int main(){
    double paso = 0.1;
    int iteraciones = 100;

    //Guardamos espacio en memoria para las masas, las posiciones, velocidades y aceleraciones
    double vector_masas[] = {1.989E30, 5.972E24};
    double *masas = vector_masas;
    Vector2D *posiciones = (Vector2D *) malloc(CUERPOS * 2 * sizeof(double)); 
    Vector2D *velocidades = (Vector2D *) malloc(CUERPOS * 2 * sizeof(double));
    Vector2D *aceleraciones = (Vector2D *) malloc(CUERPOS * 2 * sizeof(double));

    //Definimos las posiciones iniciales y las excentricidades de las órbitas
    double posiciones_iniciales[CUERPOS][2] = {{0, 0}, {1.496E11, 0}};
    double velocidades_iniciales[CUERPOS][2] = {{0, 0}, {0, 0}};
    double excentricidades[CUERPOS] = {0, 0.01671};

    //Reescalamos las posiciones y las masas
    for (int i = 0; i < CUERPOS; i++)
    {
        reescalar_posiciones(&posiciones_iniciales[i], CUERPOS);
    }
    for (int i = 0; i < CUERPOS; i++)
    {
        reescalar_masas(masas, CUERPOS);
    }
    
    /**********************************************************************************************/
    /**********************  CÁLCULO DE LAS CONDICIONES INICIALES *********************************/
    /**********************************************************************************************/

    //Conociendo el valor de las excentricidades y la distancia respecto al sol calculamos
    //la velocidad inicial (la posición inicial es (x, 0)
    for (int i = 1; i < CUERPOS; i++)
    {
        velocidades_iniciales[i][0] =  0;
        velocidades_iniciales[i][1] = sqrt((1 - excentricidades[i] * excentricidades[i])/(2 * posiciones_iniciales[i][0]));
    }
    
    //Inicializamos las posiciones de los planetas en el instante 0 como posiciones_iniciales
    for (int i = 0; i < CUERPOS; i++)
    {
        posiciones[i][0] = posiciones_iniciales[i][0];
        posiciones[i][1] = posiciones_iniciales[i][1];
    }
    //Inicializamos las velocidades de los planetas en el instante 0 como velocidades_iniciales
    for (int i = 0; i < CUERPOS; i++)
    {
        velocidades[i][0] = velocidades_iniciales[i][0];
        velocidades[i][1] = velocidades_iniciales[i][1];
    }
    //Inicializamos las velocidades de los planetas en el instante 0 teniendo en cuenta las interacciones entre ellos
    for (int i = 0; i < CUERPOS; i++)//Recorremos todos los cuerpos
    {
        for (int a = 0; i < 2; a++)//Aplicamos para la coordenada x e y
        {
            aceleraciones[i][a] = 0;
            for (int j = 0; j < CUERPOS; j++)//Calculamos el efecto del cuerpo j sobre el i
            {
                if (i != j)
                {
                    double R[] = posiciones[i] - posiciones[j];//Definimos R como la distancia entre los cuerpos i y j
                    aceleraciones[i][a] += masas[j] * R[a] / pow(norma_cuad(R), 1.5);
                }
            }
        }
            
            
    }

    /************************************************************************************************/
    /********************************* APERTURA DEL FICHERO DE ESCRITURA ****************************/
    /************************************************************************************************/
    //Abrimos el fichero donde queremos escribir los datos
    FILE *archivo = fopen("datos.txt", "w");
    if(archivo == NULL){
        printf("Error al abrir el archivo.\n");
        return 1;

    }

    //Escribimos las condiciones iniciales en el fichero
    fprintf(archivo, "POSICION_X\tPOSICION_Y\n");
    for (int i = 0; i < CUERPOS; i++)
    {
        
        fprintf(archivo, "%f\t%f\n", posiciones_iniciales[i][0], posiciones_iniciales[i][1]);
        
        
    }
    
    /***********************************************************************************************/
    /************* CALCULAMOS Y ESCRIBIMOS LA POSICION DE LOS PLANETAS EN TIEMPO t > O *************/
    /***********************************************************************************************/
    for (int n = 0; n < iteraciones; n++)//Iteramos el código
    {
        for (int i = 0; i < CUERPOS; i++)//Calculamos los vectores para el cuerpo i
        {
            for (int a = 0; a < 2; a++)//Iteramos sobre las dos componentes de los vectores
            {
                aceleraciones[i][a] = 0;
                for (int j = 0; j < CUERPOS; j++)//Calculamos el efecto del cuerpo j != i con i
                {
                    if (i != j)
                    {
                        double R[] = posiciones[i] - posiciones[j];//Definimos R como la distancia entre los cuerpos i y j
                        aceleraciones[i][a] += masas[j] * R[a] / pow(norma_cuad(R), 1.5);
                    }
                }
                posiciones[i][a] = posiciones[i][a] + paso * velocidades[i][a] + paso * paso / 2. * aceleraciones[i][a];
                velocidades[i][a] = velocidades[i][a] + paso * aceleraciones[i][a];
            }
        }
    }
    
    fclose(archivo);//Cerramos el fichero 
    free(posiciones);
    free(velocidades);
    free(aceleraciones);
    return 0;
}
