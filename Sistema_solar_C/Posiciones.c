# include "funciones.h"
# include <stdlib.h>
# include <math.h>
# include <stdio.h>
# include <string.h>
# include <omp.h>

# define PI 3.14159265
# define G 6.6743E-11
# define Ms 1.989E30
# define UA 1.496E11
# define CUERPOS 10 //Cuantos cuerpos estamos tomando

typedef double Vector2D[2];
int main(){
    double paso = 0.001;
    int iteraciones = 100000;
    int calcular_energias = 1;
    int calcular_momentos = 1;

    //Guardamos espacio en memoria para las masas, las posiciones, velocidades y aceleraciones
    Vector2D *posiciones = (Vector2D *) malloc(CUERPOS * 2 * sizeof(double)); 
    Vector2D *velocidades = (Vector2D *) malloc(CUERPOS * 2 * sizeof(double));
    Vector2D *aceleraciones = (Vector2D *) malloc(CUERPOS * 2 * sizeof(double));

    //Definimos las posiciones iniciales (afelios) y las excentricidades de las órbitas
    double vector_masas[] = {1.989E30, 3.303E23, 4.869E24, 5.972E24, 0.642E24, 1898E24, 568E24, 86.8E24, 102.43E24, 1.303E22};
    double excentricidades[] = {0, 0.205630690, 0.0068, 0.01671, 0.0934, 0.0484, 0.0542, 0.0444, 0.00859, 0.25};
    Vector2D posiciones_iniciales[] = {{0, 0}, {0.69817445E11, 0}, {1.08942780E11, 0}, {1.52098233E11, 0}, {2.49232432E11, 0}, {8.16001807E11, 0}, {15.03509229E11, 0}, {30.06318143E11, 0}, {45.37039826E11, 0}, {73.76124302E11, 0}};
    Vector2D velocidades_iniciales[CUERPOS];
    double *masas = vector_masas;
    double energia_inicial[CUERPOS];
    double momento_angular_inicial[CUERPOS];
    double energia_total;
    double momento_total;


    //Reescalamos las posiciones y las masas
    reescalar_posiciones(posiciones_iniciales, CUERPOS);
    
    reescalar_masas(masas, CUERPOS);
    
    /**********************************************************************************************/
    /**********************  CÁLCULO DE LAS CONDICIONES INICIALES *********************************/
    /**********************************************************************************************/

    //Conociendo el valor de las excentricidades y la distancia respecto al sol calculamos
    //la velocidad inicial (la posición inicial es (x, 0)
    for (int i = 1; i < CUERPOS; i++)
    {
        velocidades_iniciales[i][0] =  0;
        velocidades_iniciales[i][1] = sqrt(masas[0] * (1 - excentricidades[i]) / (1 + excentricidades[i]) / posiciones_iniciales[i][0]);
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

    //Calculamos la energía mecánica inicial de cada cuerpo
    for (int i = 0; i < CUERPOS; i++)
    {
        energia_inicial[i] = 0;
        for (int j = 0; j < CUERPOS; j++)
        {
            if (i!=j)
            {
                Vector2D R;
                R[0] = posiciones_iniciales[i][0] - posiciones_iniciales[j][0];
                R[1] = posiciones_iniciales[i][1] - posiciones_iniciales[j][1];
                energia_inicial[i] += (R, masas[i], velocidades_iniciales[i], masas[j]);
            }
        }
    }

    //Calculamos los mmomentos angulares iniciales
    for (int i = 0; i < CUERPOS; i++)
    {
        momento_angular_inicial[i] = momento_angular(posiciones_iniciales[i], masas[i], velocidades_iniciales[i]);
    }


    /************************************************************************************************/
    /********************************* APERTURA DEL FICHERO DE ESCRITURA ****************************/
    /************************************************************************************************/
    //Abrimos los ficheros donde queremos escribir los datos de la posición
    FILE *archivos[CUERPOS];
    char filename[20];
    for (int i = 0; i < CUERPOS; i++)
    {
        sprintf(filename, "Datos_cuerpos/cuerpo_%d.txt", i);
        archivos[i] = fopen(filename, "w");
        if (archivos[i] == NULL)
        {
            perror("Error al abrir el archivo");
            return 1;
        }
        
    }

    //Abrimos los ficheros donde queremos escribir los datos de las energías
    FILE *archivos_energias[CUERPOS];
    char filename_energias[20];
    for (int i = 0; i < CUERPOS; i++)
    {
        sprintf(filename_energias, "Datos_cuerpos/energia_%d.txt", i);
        archivos_energias[i] = fopen(filename_energias, "w");
        if (archivos_energias[i] == NULL)
        {
            perror("Error al abrir el archivo");
            return 1;
        }
        
    }
    
    //Abrimos los ficheros donde queremos escribir los datos de los momentos angulares
    FILE *archivos_momentos[CUERPOS];
    char filename_momentos[20];
    for (int i = 0; i < CUERPOS; i++)
    {
        sprintf(filename_momentos, "Datos_cuerpos/momentos_%d.txt", i);
        archivos_momentos[i] = fopen(filename_momentos, "w");
        if (archivos_momentos[i] == NULL)
        {
            perror("Error al abrir el archivo");
            return 1;
        }
        
    }

    FILE *archivo_energia_total = fopen("Datos_sistema/energia_total.txt", "w");
    if (archivo_energia_total == NULL)
    {
        perror("Error al abrir el archivo");
        return 1;
    }

    FILE *archivo_momento_total = fopen("Datos_sistema/momento_total.txt", "w");
    if (archivo_momento_total == NULL)
    {
        perror("Error al abrir el archivo");
        return 1;
    }


    //Escribimos las condiciones iniciales en el fichero
    energia_total = 0.0;
    momento_total = 0.0;
    for (int i = 0; i < CUERPOS; i++)
    {
        fprintf(archivos[i], "POSICION_X\tPOSICION_Y\n");
        fprintf(archivos[i], "%f\t%f\n", posiciones_iniciales[i][0], posiciones_iniciales[i][1]);
        if (calcular_energias)
        {
            fprintf(archivos_energias[i], "ENERGIA\n");
            fprintf(archivos_energias[i], "%.17lf\n", energia_inicial[i]);
        
        }
        if (calcular_momentos)
        {
            fprintf(archivos_momentos[i], "MOMENTO_ANGULAR\n");
            fprintf(archivos_momentos[i], "%.17f\n", momento_angular(posiciones[i], masas[i], velocidades[i]));

        }
        energia_total += energia_inicial[i];
        momento_total += momento_angular(posiciones[i], masas[i], velocidades[i]);
    }

    fprintf(archivo_energia_total, "ENERGIA_TOTAL\n");
    fprintf(archivo_energia_total, "%.17f\n", energia_total);
    fprintf(archivo_momento_total, "MOMENTO_TOTAL\n");
    fprintf(archivo_momento_total, "%.17f\n", momento_total);
    
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
        for (int a = 0; a < 2; a++)
        {
            for (int i = 0; i < CUERPOS; i++)
            {
                aceleraciones_previa[i][a] = aceleraciones[i][a];
            }
        }

        //Recorremos todos los cuerpos calculando sus aceleraciones y velocidades
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

        //Escribimos en el fichero la energia mecanica de los cuerpos
        energia_total = 0.0;
        for (int i = 0; i < CUERPOS; i++)
        {
            energia_inicial[i] = 0;
            for (int j = 0; j < CUERPOS; j++)
            {
                if (i!=j)
                {
                    Vector2D R;
                    R[0] = posiciones_iniciales[i][0] - posiciones_iniciales[j][0];
                    R[1] = posiciones_iniciales[i][1] - posiciones_iniciales[j][1];
                    energia_inicial[i] += energia_mecanica(R, masas[i], velocidades_iniciales[i], masas[j]);
                }
            }
            fprintf(archivos_energias[i], "%.17lf\n", energia_inicial[i]);
            energia_total += energia_inicial[i];
        }
        fprintf(archivo_energia_total, "%.17f\n", energia_total);

        //Escribimos en el fichero el momento angular de los cuerpos
        momento_total = 0.0;
        for (int i = 0; i < CUERPOS; i++)
        {
            fprintf(archivos_momentos[i], "%.17f\n", momento_angular(posiciones[i], masas[i], velocidades[i]));
            momento_total += momento_angular(posiciones[i], masas[i], velocidades[i]);
        }
        fprintf(archivo_momento_total, "%.17f\n", momento_total);
    } 

    for (int i = 0; i < CUERPOS; i++)
    {
        //Cerramos los ficheros
        fclose(archivos[i]);
        fclose(archivos_energias[i]);
        fclose(archivos_momentos[i]); 
    }
    fclose(archivo_energia_total);
    fclose(archivo_momento_total);

    //Liberamos el espacio de los punteros
    free(posiciones);
    free(velocidades);
    free(aceleraciones);
    return 0;
}
