// Lennard_Jones_separado.cpp

#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <random>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <sstream>
#include <iomanip>

#define PI 3.1415926535897932

using namespace std;

//***Clase Particula***//
class Particula
{
private:
    double x;
    double y;
    double vx;
    double vy;
public:
    Particula(double x_, double y_);
    Particula();
    void set_xy(double x_, double y_);
    void set_vxvy(double vx_, double vy_);
    std::array<double, 2> get_xy() const;
    std::array<double, 2> get_vxvy() const;
    void actualizar_posicion(const array<double, 2>& aceleracion, double paso);
    void actualizar_velocidad(const array<double, 2>& aceleracion, const array<double, 2>& aceleracion_nueva, double paso);
};

void Particula::set_xy(double x_, double y_){
    x = x_;
    y = y_;
}

void Particula::set_vxvy(double vx_, double vy_){
    vx = vx_;
    vy = vy_;
}

std::array<double, 2> Particula::get_xy() const {
    return {x, y};
}

std::array<double, 2> Particula::get_vxvy() const {
    return {vx, vy};
}

void Particula::actualizar_posicion(const array<double, 2>& aceleracion, double paso){
    x += paso * vx + 0.5 * paso * paso * aceleracion[0];
    y += paso * vy + 0.5 * paso * paso * aceleracion[1];
}

void Particula::actualizar_velocidad(const std::array<double, 2>& aceleracion, const array<double, 2>& aceleracion_nueva, double paso){
    vx += 0.5 * paso * (aceleracion[0] + aceleracion_nueva[0]);
    vy += 0.5 * paso * (aceleracion[1] + aceleracion_nueva[1]);
}

Particula:: Particula(double x_, double y_) : x(x_), y(y_), vx(0.0), vy(0.0)
{
}

Particula:: Particula() : x(0.0), y(0.0), vx(0.0), vy(0.0)
{
}
///*********************///

//***Clase Archicvo***/
class Archivo
{
private:
    string nombre;
    ofstream archivo;
    bool habilitado;
public:
    Archivo(const string& nombre_, bool habilitado_);
    ~Archivo();
    void escribir_cabecera(const string& cabecera);
    void escribir_posiciones(const vector<Particula>& p);
    void escribir_velocidades(const vector<Particula>& p);
    void escribir_energia(const vector<double>& K, const vector<double>& V);
    void escribir(const string& texto);
};

Archivo::Archivo(const string& nombre_, bool habilitado_) : habilitado(habilitado_), nombre(nombre_)
{
    if (habilitado)
    {
        archivo.open(nombre);
        if(!archivo.is_open())
        {
            throw std::runtime_error("No se pudo abrir el archivo: " + nombre);
        }
    }
}

Archivo::~Archivo()
{
    if(archivo.is_open()){
        archivo.close();
    }
}

void Archivo::escribir_cabecera(const string& cabecera)
{
    if (habilitado)
    {
        archivo << cabecera << endl;
    }
    
}

void Archivo::escribir_posiciones(const vector<Particula>& p)
{
    if (habilitado)
    {
        for (int i = 0; i < p.size()-1; i++)
        {
            archivo << p[i].get_xy()[0] << "," << p[i].get_xy()[1] << ",";
        }

        archivo << p.back().get_xy()[0] << "," << p.back().get_xy()[1] << endl;
    }
}

void Archivo::escribir_velocidades(const vector<Particula>& p)
{
    if (habilitado)
    {
        for (int i = 0; i < p.size()-1; i++)
        {
            archivo << p[i].get_vxvy()[0] << "," << p[i].get_vxvy()[1] << ",";
        }

        archivo << p.back().get_vxvy()[0] << "," << p.back().get_vxvy()[1] << endl;
    }
}

void Archivo::escribir_energia(const vector<double>& K, const vector<double>& V)
{
    if (habilitado)
    {
        for (int i = 0; i < K.size(); i++)
        {
            archivo << K[i] << "," << V[i] << endl;
        }
        
    }
    
}

void Archivo::escribir(const string& texto)
{
    archivo << texto;
}


///*Clase Simulador*///
class Simulador
{
private:
    vector<Particula> particulas;
    vector<array<double, 2>> aceleraciones, nuevas_aceleraciones;
    const int N;
    const int iteraciones;
    const double paso;
    const double L;
    double presion;
    int num_threads;
    double v_max;
    vector<double> energia_cinetica;
    vector<double> energia_potencial;
    bool save_to_file;
    Archivo datos;
    Archivo velocidades;
    Archivo energias;
    Archivo var_posicion;
    vector<double> variacion_posicion;
public:
    Simulador(const int N_, const int iteraciones_, const double paso_, const double L_, int num_threads, double v_max_, bool guardar_datos);
    ~Simulador();
    void inicializar();
    void integrar();
    vector<array<double, 2>> calcular_aceleraciones(int n) const;
    void calcular_energia_cinetica(int n);
    void calcular_energia_potencial(int n);
    void guardar_energias();
    double get_presion();
    vector<double> get_energia_cinetica();
};

Simulador::Simulador(const int N_, const int iteraciones_, const double paso_, const double L_, int num_threads_, double v_max_, bool guardar_datos) : 
N(N_), 
iteraciones(iteraciones_), 
paso(paso_), L(L_),
particulas(N_), 
num_threads(num_threads_),
aceleraciones(N_, {0.0, 0.0}), 
nuevas_aceleraciones(N_, {0.0, 0.0}), 
v_max(v_max_), 
energia_cinetica(iteraciones_), 
energia_potencial(iteraciones_),
variacion_posicion(N_),
save_to_file(guardar_datos),
datos("Datos/Posiciones/Lennard_Jones_" + to_string(N_) + "_" + to_string(iteraciones_) + "_" + to_string((int)L_) + ".csv", save_to_file),
velocidades("Datos/Velocidades/Velocidades_" + to_string(N) + "_" + to_string(iteraciones) + "_" + to_string((int)L) + ".csv", save_to_file),
energias("Datos/Energias/Energias_" + to_string(N) + "_" + to_string(iteraciones) + "_" + to_string((int)L) + ".csv", save_to_file),
var_posicion("Datos/Variacion_posicion/Var_posicion_" + to_string(N_) + "_" + to_string(iteraciones_) + "_" + to_string((int)L_) + ".csv", save_to_file)
{   
    string cabecera_datos, cabecera_velocidades;

    // Escribimos la cabecera del archivo de velocidades y posiciones
    for (size_t i = 0; i < particulas.size()-1; i++)
    {
        cabecera_datos += "x" + to_string(i) + ",y" + to_string(i) + ",";
        cabecera_velocidades += "vx" + to_string(i) + ",vy" + to_string(i) + ",";
    }
    cabecera_datos += "x" + to_string(N-1) + ",y" + to_string(N-1);
    cabecera_velocidades += "vx" + to_string(N-1) + ",vy" + to_string(N-1);

    datos.escribir_cabecera(cabecera_datos);
    velocidades.escribir_cabecera(cabecera_velocidades);

    // Escribimos la cabecera de la variación de las posiciones
    var_posicion.escribir_cabecera("<r>^2");

    omp_set_num_threads(num_threads);

    // Escribimos los parámetros de la simulación 
    // para no tener que volver a escribirlos en el programa de python
    Archivo parametros("Analisis/Parametros.txt", save_to_file);
    parametros.escribir("N\t" + to_string(N) + "\n");
    parametros.escribir("L\t" + to_string(L) + "\n");
    parametros.escribir("Pasos\t" + to_string(iteraciones) + "\n");
    parametros.escribir("vmax\t" + to_string(v_max) + "\n");
}

Simulador::~Simulador()
{
}

void Simulador::inicializar()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> real(0.0, 2 * PI); 

    for (size_t i = 0; i < particulas.size(); i++)
    {
        double randnum = real(gen);

        // Inicio al azar
        //particulas[i].set_xy(((double) rand() / RAND_MAX) * L, ((double) rand() / RAND_MAX) * L); 

        // Inicio dentro de cuadriculas en una posicion al azar 
        //particulas[i][0].set_xy(L / sqrt(N) * (posx + i % ((int) sqrt(N))), L / sqrt(N) * (posy + i / (int) sqrt(N)));

        // Inicio cuadriculado
        particulas[i].set_xy(L / sqrt(N) * (0.5 + i % ((int) sqrt(N))), L / sqrt(N) * (0.5 + i / (int) sqrt(N)));

        // Inicio en panal
        //particulas[i].set_xy(L / sqrt(N) * (0.5 + i % ((int) sqrt(N))), L / sqrt(N) * (0.5 + i / (int) sqrt(N)));
        //if (i % (2 * (int)L) > (int)L)
        //particulas[i].set_xy(L / sqrt(N) * (1 + i % ((int) sqrt(N))), L / sqrt(N) * (1 + i / (int) sqrt(N)));

        // Velocidades radiales
        //particulas[i].set_vxvy(2 * (particulas[i].get_xy()[0] - L/2), 2 * (particulas[i].get_xy()[1] - L/2));

        // Velocidades aleatorias
        particulas[i].set_vxvy(v_max * cos(randnum), v_max * sin(randnum));
    }

}

void Simulador::integrar()
{
    vector<array<double, 2>> posicion_inicial(N);
    for (int i = 0; i < N; i++)
    {
        posicion_inicial[i][0] = particulas[i].get_xy()[0];
        posicion_inicial[i][1] = particulas[i].get_xy()[1];
    }
    
    presion = 0.0;
    datos.escribir_posiciones(particulas);
    velocidades.escribir_velocidades(particulas);
    calcular_energia_cinetica(0);
    calcular_energia_potencial(0);
    for (int n = 1; n < iteraciones; n++)
    {
        // Calculamos la aceleración de todas las partículas
        aceleraciones = calcular_aceleraciones(n);

        // Calculamos la posición y la velocidad acumulada de todas las partículas
        for (size_t i = 0; i < N; i++)
        {
            //calcular_presion(particulas[i], n);

            particulas[i].actualizar_posicion(aceleraciones[i], paso);
            auto r = particulas[i].get_xy();
            auto v = particulas[i].get_vxvy();
            if (r[0] >= L)     
            {
                r[0] -= L;
                presion += 2 * fabs(v[0]);
            }
            else if (r[0] < 0) 
            {
                r[0] += L;
                presion += 2 * fabs(v[0]);
            }

            if (r[1] >= L)     
            {
                r[1] -= L;
                presion += 2 * fabs(v[1]);
            }
            else if (r[1] < 0) 
            {
                r[1] += L;
                presion += 2 * fabs(v[1]);
            }

            particulas[i].set_xy(r[0], r[1]);
        }

        // Guardamos el estado del sistema
        datos.escribir_posiciones(particulas);
        velocidades.escribir_velocidades(particulas);

        // Calculamos las nuevas aceleraciones de las partículas
        nuevas_aceleraciones = calcular_aceleraciones(n);

        // Calculamos las nuevas velocidades de las partículas
        for (size_t i = 0; i < particulas.size(); i++)
        {
            particulas[i].actualizar_velocidad(aceleraciones[i], nuevas_aceleraciones[i], paso);
        }

        calcular_energia_cinetica(n);
        calcular_energia_potencial(n);

        // Calculamos la variación en la posición
        for (int i = 0; i < N; i++)
        {
            double x0 = posicion_inicial[i][0];
            double y0 = posicion_inicial[i][1];
            double x = particulas[i].get_xy()[0];
            double y = particulas[i].get_xy()[1];
            variacion_posicion[i] = (x - x0) * (x - x0) + (y - y0) * (y - y0);
        }

        // Hacemos la media de la variación de la posición de todas las partículas
        double variacion_media = 0.0;
        for (int i = 0; i < N; i++)
        {
            variacion_media += variacion_posicion[i];
        }
        
        variacion_media = variacion_media / N;
        
        // Escribimos la variación de la posición media de todas las partículas en el archivo
        var_posicion.escribir(to_string(variacion_media) + "\n");

    }

}

vector<array<double, 2>> Simulador::calcular_aceleraciones(int n) const{

    vector<array<double, 2>> aceleraciones(N, {0.0, 0.0});

    #pragma omp parallel for
    for (size_t i = 0; i < particulas.size(); i++)
    {
        auto ri = particulas[i].get_xy();
        // Calculamos la aceleración que sufre el cuerpo i debido a todos los cuerpos j != i
        for (size_t j = 0; j < particulas.size(); j++)
        {
            if(i != j)
            {
                auto rj = particulas[j].get_xy();
                std::array<double, 2> rij = {ri[0] - rj[0], ri[1] - rj[1]};
                // aplicar condición de mínima imagen
                for (size_t a = 0; a < 2; ++a) {
                    if (rij[a] >  0.5 * L) rij[a] -= L;
                    if (rij[a] < -0.5 * L) rij[a] += L;
                }
                if(rij[0] * rij[0] + rij[1] * rij[1] <= 9.0)
                {
                    double rij2 = rij[0] * rij[0] + rij[1] * rij[1];
                    double rij4 = rij2 * rij2;
                    double rij8 = rij4 * rij4;
                    double rij14 = rij8 * rij4 * rij2;
                    double fuerza = 48 / rij14 - 24 / rij8;
                    for (int a = 0; a < 2; ++a) 
                    {
                        aceleraciones[i][a] += fuerza * rij[a];
                    }
                }
                

            }
        }
    }
    return aceleraciones;

}

void Simulador::calcular_energia_cinetica(int n)
{
    double energia_local = 0.0;

    #pragma omp parallel for reduction(+:energia_local)
    for (int i = 0; i < N; i++)
    {
        auto v = particulas[i].get_vxvy();
        energia_local += 0.5 * (v[0] * v[0] + v[1] * v[1]);
    }

    energia_cinetica[n] = energia_local;
}

void Simulador::calcular_energia_potencial(int n)
{
    double energia_local = 0.0;
    double rc = 3.0;

    #pragma omp parallel for reduction(+:energia_local) schedule(dynamic)
    for (size_t i = 0; i < N; i++)
    {
        auto ri = particulas[i].get_xy();
        for (size_t j = i+1; j < N; j++)
        {
            {
                auto rj = particulas[j].get_xy();
                std::array<double, 2> rij = {ri[0] - rj[0], ri[1] - rj[1]};
                // aplicar condición de mínima imagen
                for (size_t a = 0; a < 2; ++a) {
                    if (rij[a] >  0.5 * L) rij[a] -= L;
                    if (rij[a] < -0.5 * L) rij[a] += L;
                }

                if(rij[0] * rij[0] + rij[1] * rij[1] <= 9.0)
                {
                    double rij2 = rij[0]*rij[0] + rij[1] * rij[1];
                    double rij4 = rij2 * rij2;
                    double rij6 = rij2 * rij4;
                    double rij12 = rij6 * rij6;
                    energia_local += 4 * (1 / rij12 - 1 / rij6);
                }
                
            }
            
        }
        
    }

    energia_potencial[n] = energia_local;

}

double Simulador::get_presion()
{
    return presion;
}

vector<double> Simulador::get_energia_cinetica()
{
    return energia_cinetica;
}

void Simulador::guardar_energias()
{
    energias.escribir_energia(energia_cinetica, energia_potencial);
}

int main(){
    int N = 16;
    double L = 10.0;
    int iteraciones = 25000;
    double vmax = 1.0;

    // Se le pasan como argumentos a la simulación
    // - int N: número de partículas
    // - int iteraciones: número de pasos del programa
    // - double paso: paso de tiempo entre dos pasos del programa
    // - double L: tamaño de la caja
    // - int n_threads: número de hilos para la paralelización
    // - double vmax: módulo de la velocidad de las partículas
    // - bool guardar_datos: elección si se quieren escribir los datos del sistema, velocidad y posición, en archivos


    Simulador simulacion(N, iteraciones, 0.002, L, 1, vmax, true);
    simulacion.inicializar();
    simulacion.integrar();
    simulacion.guardar_energias();


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