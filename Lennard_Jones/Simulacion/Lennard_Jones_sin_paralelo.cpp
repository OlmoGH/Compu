// Lennard_Jones_separado.cpp

# include "lennard-jones-lib/include/Simulador.h"
# include <string>
# include <iostream>

int main(){
    int N = 16;
    double L = 4.0;
    int iteraciones = 5000000;
    double vmax = 0;
    double paso = 0.001;

    // Se le pasan como argumentos a la simulación
    // - int N: número de partículas
    // - int iteraciones: número de pasos del programa
    // - double paso: paso de tiempo entre dos pasos del programa
    // - double L: tamaño de la caja
    // - int n_threads: número de hilos para la paralelización
    // - double vmax: módulo de la velocidad de las partículas
    // - bool guardar_datos: elección si se quieren escribir los datos del sistema, velocidad y posición, en archivos

    Simulador simulacion(N, iteraciones, paso, L, vmax, true);

    // Inicializamos las partículas
    simulacion.inicializar();
    
    // Escribimos el estado del sistema
    //simulacion.datos.escribir_posiciones(simulacion.particulas);
    // simulacion.velocidades.escribir_velocidades(simulacion.particulas);

    // Calculamos la variación de la posición con respecto a la posición inicial
    //double variacion_media = simulacion.calcular_var_posicion();

    // Calculamos la separación entre dos partículas
    double separacion_value = simulacion.calcular_separacion();
        
    // Escribimos la variación de la posición media de todas las partículas en el archivo
    //simulacion.var_posicion.escribir(std::to_string(variacion_media) + "\n");    
    
    // Escribimos la separación entre dos partículs centrales
    simulacion.separacion.escribir(std::to_string(separacion_value) + "\n");

    // Calculamos la energía cinética y potencial inicial
    //simulacion.calcular_energia_cinetica(0);
    //simulacion.calcular_energia_potencial(0);

    int incremento = 300000;

    // Iniciamos el bucle de simulación
    for (int n = 1; n < iteraciones; n++)
    {
        if (n == incremento)
        {
            for (int i = 0; i < N; i++)
            {
                auto v = simulacion.particulas[i].get_vxvy();
                simulacion.particulas[i].set_vxvy(1.1 * v[0], 1.1 * v[1]);
            }
            incremento += incremento;
        }
        // Aplicamos integración mediante el método Verlet-Velocity
        simulacion.aplicar_Verlet();

        // Guardamos el estado del sistema
        //simulacion.datos.escribir_posiciones(simulacion.particulas);
        // simulacion.velocidades.escribir_velocidades(simulacion.particulas);

        //simulacion.calcular_energia_cinetica(n);
        //simulacion.calcular_energia_potencial(n);
        
        // Calculamos la variación de la posición con respecto a la posición inicial
        //double variacion_media = simulacion.calcular_var_posicion();
        
        // Escribimos la variación de la posición media de todas las partículas en el archivo
        //simulacion.var_posicion.escribir(std::to_string(variacion_media) + "\n");

        // Calculamos la separación entre dos partículas
        double separacion_value = simulacion.calcular_separacion();
        
        // Escribimos la separación entre dos partículs centrales
        simulacion.separacion.escribir(std::to_string(separacion_value) + "\n");

        if (n % 10000 == 0)
        std::cout << "Simulando iteracion " << n << " de " << iteraciones << std::endl;

    }
    //simulacion.guardar_energias();

    /*
    int n_pruebas = 10;
    vector<double> T(n_pruebas, 0.0);
    double T_media = 0.0;
    int inicio = 10000;
    ofstream presiones("Datos_presion.txt");

    presiones << "v\tT\tp" << endl;

    for (int m = 0; m < 100; m++)
    {
        Simulador simulacion(N, pasos, 0.002, L, 1, sqrt(m));
        simulacion.guardar_archivo(false);

        for (int n = 0; n < n_pruebas; n++)
        {
            simulacion.inicializar();
            simulacion.integrar();

            auto T_array = simulacion.get_energia_cinetica();
            for (int i = inicio; i < pasos; i++)
            {
                T_media += T_array[i];
            }

            T_media /= (pasos - 10000) * N;
            presiones << sqrt(m) << "\t" << T_media << "\t" << simulacion.get_presion() / (pasos) << endl;
            cout << "Completado " << fixed << setprecision(4) << (n_pruebas * m + n) * 1.0 / (n_pruebas * 100) << "%" << "\r" ;
        }
    }

    */
    
    return 0;
}