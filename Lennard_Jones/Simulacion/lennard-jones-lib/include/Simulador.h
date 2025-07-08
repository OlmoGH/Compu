# ifndef SIMULADOR_H
# define SIMULADOR_H

# include <vector>
# include <array>
# include "Particula.h"
# include "Archivo.h"
# include <omp.h>

/**
 * @brief Modos de inicialización del sistema
 * 
 */
enum class Modo {Azar, Cuadricula, CuadriculaAzar, Panal};

///*Clase Simulador*///
class Simulador
{
public:
    const int N;
    const int iteraciones;
    const double paso;
    const double L;
    double v_max;
    std::vector<Particula> particulas;
    std::vector<std::array<double, 2>> aceleraciones, nuevas_aceleraciones;
    double presion;
    std::vector<double> energia_cinetica;
    std::vector<double> energia_potencial;
    std::vector<Particula> posicion_inicial;
    bool save_to_file;
    Archivo datos;
    Archivo velocidades;
    Archivo energias;
    Archivo var_posicion;
    Archivo separacion;
    Simulador(const int N_, const int iteraciones_, const double paso_, const double L_, double v_max_, bool guardar_datos, int n_threads);
    ~Simulador();
    void inicializar(Modo modo);
    void aplicar_Verlet(int n);
    double calcular_var_posicion();
    double calcular_separacion(int id_1, int id_2);
    std::vector<std::array<double, 2>> calcular_aceleraciones() const;
    void calcular_energia_cinetica(int n);
    void calcular_energia_potencial(int n);
    void guardar_energias();
    double get_presion();
    std::vector<double> get_energia_cinetica();
};

# endif