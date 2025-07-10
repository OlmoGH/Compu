# ifndef SIMULADOR_H
# define SIMULADOR_H

# include <vector>
# include <array>
# include "Particula.h"
# include "Archivo.h"
# include <omp.h>

/**
 * @brief Modos de inicialización de la estructura del sistema
 * 
 */
enum class Estructura {Azar, Cuadricula, CuadriculaAzar, Panal};


/**
 * @brief Modos de inicialización de la velocidad de las partículas
 * 
 */
enum class Modo {Compensado, Horizontal, Descompensado};

///*Clase Simulador*///
class Simulador
{
public:
    const int N;
    const int iteraciones;
    const double paso;
    const double Lx;
    const double Ly;
    double v_max;
    double empuje;
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
    Simulador(const int N_, const int iteraciones_, const double paso_, const double Lx_, const double Ly_, double v_max_, bool guardar_datos, int n_threads, double empuje);
    ~Simulador();
    void inicializar(Estructura estructura, Modo modo);
    void aplicar_Verlet(int n);
    double calcular_var_posicion();
    double calcular_separacion(int id_1, int id_2);
    std::vector<std::array<double, 2>> calcular_aceleraciones() const;
    void calcular_energia_cinetica(int n);
    void calcular_energia_potencial(int n);
    void guardar_energias();
    double get_presion();
    std::vector<double> get_energia_cinetica();
    double get_temperatura();
};

# endif