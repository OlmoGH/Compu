# include "funciones.h"
# include <stdlib.h>
# include <math.h>
# include <stdio.h>
# include <string.h>
# define PI 3.14159265
# define CUERPOS 6 //Cuantos cuerpos estamos tomando

typedef double Vector2D[2];
int main(){
    double paso = 0.001;
    int iteraciones = 100000;

    //Guardamos espacio en memoria para las masas, las posiciones, velocidades y aceleraciones
    double vector_masas[] = {1.989E30, 5.972E24, 0.642E24, 1898E24, 568E24, 86.8E24};
    double *masas = vector_masas;
    Vector2D *posiciones = (Vector2D *) malloc(CUERPOS * 2 * sizeof(double)); 
    Vector2D *velocidades = (Vector2D *) malloc(CUERPOS * 2 * sizeof(double));
    Vector2D *aceleraciones = (Vector2D *) malloc(CUERPOS * 2 * sizeof(double));

    //Definimos las posiciones iniciales y las excentricidades de las órbitas
    double excentricidades[CUERPOS] = {0, 0.01671, 0.094, 0.049, 0.052, 0.047};
    Vector2D posiciones_iniciales[CUERPOS] = {{0, 0}, {1.498E11, 0}, {2.493E11, 0}, {8.164E11, 0}, {15.065E11, 0}, {30.014E11, 0}};
    Vector2D velocidades_iniciales[CUERPOS];

    //Reescalamos las posiciones y las masas
    reescalar_posiciones(posiciones_iniciales, CUERPOS);
    
    reescalar_masas(masas, CUERPOS);

    printf("Las masas reescaladas son %f para el sol y %f para la tierra\n", masas[0], masas[1]);
    
    /**********************************************************************************************/
    /**********************  CÁLCULO DE LAS CONDICIONES INICIALES *********************************/
    /**********************************************************************************************/

    //Conociendo el valor de las excentricidades y la distancia respecto al sol calculamos
    //la velocidad inicial (la posición inicial es (x, 0)
    for (int i = 1; i < CUERPOS; i++)
    {
        velocidades_iniciales[i][0] =  0;
        velocidades_iniciales[i][1] = 0.7 * sqrt(masas[0] * (1 + excentricidades[i]) / (1 - excentricidades[i]) / posiciones_iniciales[i][0]);
    }
    printf("Posicion inicial cuerpo 1: rx = %f, ry = %f\n", posiciones_iniciales[1][0], posiciones_iniciales[1][1]);
    printf("Velocidad inicial cuerpo 1: vx = %f, vy = %f\n", velocidades_iniciales[1][0], velocidades_iniciales[1][1]);
    
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
    //Inicializamos las aceleraciones de los planetas en el instante 0 teniendo en cuenta las interacciones entre ellos
    for (int i = 0; i < CUERPOS; i++)//Recorremos todos los cuerpos
    {
        for (int a = 0; a < 2; a++)//Aplicamos para la coordenada x e y
        {
            aceleraciones[i][a] = 0;
            for (int j = 0; j < CUERPOS; j++)//Calculamos el efecto del cuerpo j sobre el i
            {
                if (i != j)
                {
                    Vector2D R;//Definimos R como la distancia entre los cuerpos i y j
                    R[0] = posiciones[i][0] - posiciones[j][0];
                    R[1] = posiciones[i][1] - posiciones[j][1];
                    aceleraciones[i][a] -= masas[j] * R[a] / pow(norma_cuad(R), 1.5);
                }
            }
        }     
            
    }

    /************************************************************************************************/
    /********************************* APERTURA DEL FICHERO DE ESCRITURA ****************************/
    /************************************************************************************************/
    //Abrimos el fichero donde queremos escribir los datos
    FILE *archivos[CUERPOS];
    char filename[20];
    for (int i = 0; i < CUERPOS; i++)
    {
        sprintf(filename, "cuerpo_%d.txt", i);
        archivos[i] = fopen(filename, "w");
        if (archivos[i] == NULL)
        {
            perror("Error al abrir el archivo");
            return 1;
        }
        
    }

    //Escribimos las condiciones iniciales en el fichero
    
    
    for (int i = 0; i < CUERPOS; i++)
    {
        fprintf(archivos[i], "POSICION_X\tPOSICION_Y\n");
        fprintf(archivos[i], "%f\t%f\n", posiciones_iniciales[i][0], posiciones_iniciales[i][1]);
    }
    
    /***********************************************************************************************/
    /************* CALCULAMOS Y ESCRIBIMOS LA POSICIÓN, VELOCIDAD Y ACELERACIÓN ********************/
    /************* DE LOS PLANETAS EN TIEMPO t > O *************************************************/
    /***********************************************************************************************/
    
    //Iteramos el código el número de iteraciones elegidas
    for (int n = 0; n < iteraciones; n++)
    {
        //Calculamos la posición del cuerpo i
        for (int i = 0; i < CUERPOS; i++)
        {
            //Iteramos sobre las dos componentes de la posición
            for (int a = 0; a < 2; a++)
            {
                posiciones[i][a] +=  paso * velocidades[i][a] + 0.5 * paso * paso * aceleraciones[i][a];
            }

            //Escribimos la posición del cuerpo i en su archivo
            fprintf(archivos[i], "%f\t%f\n", posiciones[i][0], posiciones[i][1]);
        }

        //Guardo los valores previos de las aceleraciones para el cálculo de la velocidad
        Vector2D aceleraciones_previa[CUERPOS];
        memcpy(aceleraciones_previa, aceleraciones, sizeof(aceleraciones));

        //Recorremos todos los cuerpos calculando sus aceleraciones
        for (int i = 0; i < CUERPOS; i++)
        {
            //Reseteamos el valor de las aceleraciones a 0
            for (int a = 0; a < 2; a++)
            {
                aceleraciones[i][a] = 0.0;
            }
            
            //Calculamos la aceleración que genera j sobre i
            for (int j = 0; j < CUERPOS; j++)
            {
                if (j != i)
                {
                    //Calculamos la distancia entre los cuerpos i y j
                    Vector2D R;
                    for (int a = 0; a < 2; a++)
                    {
                        R[a] = posiciones[i][a] - posiciones[j][a];
                    }

                    //Calculamos la aceleración mediante la ley de gravitación universal
                    for (int a = 0; a < 2; a++)
                    {
                        aceleraciones[i][a] -=  masas[j] * R[a] / pow(norma_cuad(R), 1.5);
                    }
                }
            }
        }

        //Calculamos la velocidad de los cuerpos
        for (int i = 0; i < CUERPOS; i++)
        {
            for (int a = 0; a < 2; a++)
            {
                velocidades[i][a] += 0.5 * paso * (aceleraciones[i][a] + aceleraciones_previa[i][a]);
            }
        }
    }
    
    for (int i = 0; i < CUERPOS; i++)
    {
        //Cerramos los ficheros
        fclose(archivos[i]); 
    }

    //Liberamos el espacio de los punteros
    free(posiciones);
    free(velocidades);
    free(aceleraciones);
    return 0;
}
