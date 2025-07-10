// Lennard_Jones_separado.cpp

# include "lennard-jones-lib/include/Simulador.h"
# include "lennard-jones-lib/include/Archivo.h"
# include <string>
# include <iostream>
# include <cmath>
# include <string>
# include <omp.h>

int main(){
    /**
     * @brief Magnitudes del sistema que pueden ser analizadas
     * 
     */
    enum class Magnitudes {Posiciones, Velocidades, Energias, VarPosicion, Separacion};
    int N = 16;
    double Lx = 4;
    double Ly = 4;
    int iteraciones = 3000000;
    double vmax = 0;
    double empuje = 0.0;
    double paso = 0.0001;
    int n_threads = 4;
    int incremento = 300000;
    int indice = 0;

    int id_1 = 1; // Número de la primera partícula para la que se calcula la separación entre partículas
    int id_2 = 7; // Número de la segunda partícula para la que se calcula la separación entre partículas

    Simulador simulacion(N, iteraciones, paso, Lx, Ly, vmax, true, n_threads, empuje);

    // Inicializamos las partículas
    simulacion.inicializar(Estructura::Cuadricula, Modo::Compensado);

    // Escribimos el estado del sistema
    simulacion.datos.escribir_posiciones(simulacion.particulas);
    simulacion.velocidades.escribir_velocidades(simulacion.particulas);

    // Calculamos la variación de la posición con respecto a la posición inicial
    double variacion_media = simulacion.calcular_var_posicion();

    // Calculamos la separación entre dos partículas
    double separacion_value = simulacion.calcular_separacion(id_1, id_2);
        
    // Escribimos la variación de la posición media de todas las partículas en el archivo
    simulacion.var_posicion.escribir(std::to_string(variacion_media) + "\n");    
    
    // Escribimos la separación entre dos partículas centrales
    simulacion.separacion.escribir(std::to_string(separacion_value) + "\n");

    // Calculamos la energía cinética y potencial inicial
    simulacion.calcular_energia_cinetica(0);
    simulacion.calcular_energia_potencial(0);

    // Iniciamos el bucle de simulación
    for (int n = 1; n < iteraciones; n++)
    {
        if (n % 10000 == 0)
        std::cout << "Ejecutando iteracion " << n << " de " << iteraciones << std::endl;

        if(n == incremento)
        {
            for (int i = 0; i < N; i++)
            {
                auto v = simulacion.particulas[i].get_vxvy();
                simulacion.particulas[i].set_vxvy(v[0] * 1.1, v[1] * 1.1);
            }
            incremento += 300000;
            std::cout << "Holiii" << std::endl;
        }

        // Aplicamos integración mediante el método Verlet-Velocity
        simulacion.aplicar_Verlet(n);

        // Guardamos el estado del sistema
        simulacion.datos.escribir_posiciones(simulacion.particulas);
        simulacion.velocidades.escribir_velocidades(simulacion.particulas);

        simulacion.calcular_energia_cinetica(n);
        simulacion.calcular_energia_potencial(n);
        
        // Calculamos la variación de la posición con respecto a la posición inicial
        double variacion_media = simulacion.calcular_var_posicion();
        
        // Escribimos la variación de la posición media de todas las partículas en el archivo
        simulacion.var_posicion.escribir(std::to_string(variacion_media) + "\n");

        // Calculamos la separación entre dos partículas
        double separacion_value = simulacion.calcular_separacion(id_1, id_2);
        
        // Escribimos la separación entre dos partículs centrales
        simulacion.separacion.escribir(std::to_string(separacion_value) + "\n");
    }

    simulacion.guardar_energias();

    return 0;
}
