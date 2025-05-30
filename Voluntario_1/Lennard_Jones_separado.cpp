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
    ofstream datos;
    ofstream velocidades;
public:
    Simulador(const int N_, const int iteraciones_, const double paso_, const double L_, int num_threads, double v_max_);
    ~Simulador();
    void inicializar();
    void integrar();
    void guardar_estado(bool save_to_file, int n);
    void guardar_velocidades(bool save_to_file, int n);
    vector<array<double, 2>> calcular_aceleraciones(int n) const;
    void calcular_energia_cinetica(int n);
    void calcular_energia_potencial(int n);
    void guardar_energias();
    int analizar_simulacion(int paso_particulas);
    double get_presion();
    void guardar_archivo(bool eleccion)
    {
        save_to_file = eleccion;
    }
    vector<double> get_energia_cinetica();
};

Simulador::Simulador(const int N_, const int iteraciones_, const double paso_, const double L_, int num_threads_, double v_max_) : 
N(N_), 
iteraciones(iteraciones_), 
paso(paso_), L(L_),
particulas(N_), 
num_threads(num_threads_),
aceleraciones(N_, {0.0, 0.0}), 
nuevas_aceleraciones(N_, {0.0, 0.0}), 
v_max(v_max_), 
energia_cinetica(iteraciones_), 
energia_potencial(iteraciones_)
{   
    string nombre_datos = "Datos/Lennard_Jones_" + to_string(N) + "_" + to_string(iteraciones) + "_" + to_string((int)L) + ".csv";
    string nombre_vel = "Datos/Velocidades_" + to_string(N) + "_" + to_string(iteraciones) + "_" + to_string((int)L) + ".csv";
    datos.open(nombre_datos);
    velocidades.open(nombre_vel);
    if(!datos || !velocidades)
    {
        cout << "No se pudo abrir el archivo de datos o velocidades\n";
        return;
    }

    // Escribimos la cabecera del archivo
    for (size_t i = 0; i < particulas.size()-1; i++)
    {
        datos << "x" << i << ",y" << i << ",";
        velocidades << "vx" << i << ",vy" << i << ",";
    }
    datos << "x" << N-1 << ",y" << N-1 << "\n";
    velocidades << "vx" << N-1 << ",vy" << N-1 << "\n";

    omp_set_num_threads(num_threads);

    // Escribimos los parámetros de la simulación 
    // para no tener que volver a escribirlos en el programa de python
    ofstream archivo_parametros("Parametros.txt");
    if (!archivo_parametros)
    {
        cerr << "No se pudo abrir el archivo de los parametros" << endl;
    }
    archivo_parametros << "N\t" << N << endl;
    archivo_parametros << "L\t" << L << endl;
    archivo_parametros << "Pasos\t" << iteraciones << endl;
    archivo_parametros << "vmax\t" << v_max << endl;

    archivo_parametros.close();
}

Simulador::~Simulador()
{
    datos.close();
    velocidades.close();
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
    presion = 0.0;
    guardar_estado(save_to_file, 0);
    guardar_velocidades(save_to_file, 0);
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
        guardar_estado(save_to_file, n);
        guardar_velocidades(save_to_file, n);

        // Calculamos las nuevas aceleraciones de las partículas
        nuevas_aceleraciones = calcular_aceleraciones(n);

        // Calculamos las nuevas velocidades de las partículas
        for (size_t i = 0; i < particulas.size(); i++)
        {
            particulas[i].actualizar_velocidad(aceleraciones[i], nuevas_aceleraciones[i], paso);
        }

        calcular_energia_cinetica(n);
        calcular_energia_potencial(n);
    }

}

void Simulador::guardar_estado(bool save_to_file, int n)
{

    // Comprobamos si debemos escribir los datos
    if(!save_to_file)
    {
        return;
    }
    // Escribimos los datos de la posición de las partículas en la iteración n
    for (int i = 0; i < N-1; i++)
    {
        datos << particulas[i].get_xy()[0] << "," << particulas[i].get_xy()[1] << ",";
    }
    datos << particulas[N-1].get_xy()[0] << "," << particulas[N-1].get_xy()[1];

    if (n != iteraciones-1)
    datos << "\n";
    
}

void Simulador::guardar_velocidades(bool save_to_file, int n)
{

    // Comprobamos si debemos escribir los datos
    if(!save_to_file)
    {
        return;
    }
    // Escribimos los datos de la posición de las partículas en la iteración n
    for (int i = 0; i < N-1; i++)
    {
        velocidades << particulas[i].get_vxvy()[0] << "," << particulas[i].get_vxvy()[1] << ",";
    }
    velocidades << particulas[N-1].get_vxvy()[0] << "," << particulas[N-1].get_vxvy()[1];

    if (n != iteraciones-1)
    velocidades << "\n";
    
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

void Simulador::guardar_energias()
{
    string nombre_energias = "Datos/Energias_" + to_string(N) + "_" + to_string(iteraciones) + "_" + to_string((int)L) + ".csv";
    ofstream archivo_energias(nombre_energias);

    if (!archivo_energias)
    {
        cout << "No se pudo abrir el archivo de energias\n";
        return;
    }

    archivo_energias << "K,V\n";
    for (int n = 0; n < iteraciones-1; n++)
    {
        archivo_energias << energia_cinetica[n] << "," << energia_potencial[n] << endl;
    }
    archivo_energias << energia_cinetica[iteraciones-1] << "," << energia_potencial[iteraciones-1];

    archivo_energias.close();
}

double Simulador::get_presion()
{
    return presion;
}

vector<double> Simulador::get_energia_cinetica()
{
    return energia_cinetica;
}

int main(){
    int N = 16;
    double L = 10.0;
    int pasos = 100000;
    double vmax = 1.0;
    int n_pruebas = 100;
    vector<double> T(n_pruebas, 0.0);
    double T_media = 0.0;
    int inicio = 10000;
    ofstream presiones("Datos_presion.txt");

    for (int m = 1; m < 10; m++)
    {
        Simulador simulacion(N, pasos, 0.002, L, 1, m);
        simulacion.guardar_archivo(false);

        presiones << "v\tT\tp" << endl;

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
        presiones << m << "\t" << T_media << "\t" << simulacion.get_presion() << endl;
    }
    
    }
    


    



    return 0;
}