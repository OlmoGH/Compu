# include "../include/Simulador.h"
# include "Particula.h"
# include "Archivo.h"
# include <vector>
# include <array>
# include <string>
# include <iostream>
# include <random>
# include <omp.h> 

#define PI 3.1415926535897932

using namespace std;
/**
 * @brief Construct a new Simulador:: Simulador object
 * 
 * @param N_ Número de partículas
 * @param iteraciones_ Número de pasos del programa
 * @param paso_ Paso de tiempo entre dos pasos del programa
 * @param L_ Tamaño de la caja
 * @param v_max_ Módulo de la velocidad de las partículas
 * @param guardar_datos Elección si se quieren escribir los datos del sistema, velocidad y posición, en archivos
 * @param n_threads Número de hilos para la paralelización
 * @param empuje_-Empuje Velocidad inicial del sistema en el eje x
 */
Simulador::Simulador(const int N_, const int iteraciones_, const double paso_, const double Lx_, const double Ly_,  double v_max_, bool guardar_datos, int n_threads, double empuje_) : 
N(N_), 
iteraciones(iteraciones_), 
paso(paso_), 
Lx(Lx_),
Ly(Ly_),
v_max(v_max_), 
empuje(empuje_),
particulas(N_), 
aceleraciones(N_, {0.0, 0.0}), 
nuevas_aceleraciones(N_, {0.0, 0.0}),
presion(0.0), 
energia_cinetica(iteraciones_), 
energia_potencial(iteraciones_),
posicion_inicial(N_),
save_to_file(guardar_datos),
datos("Datos/Posiciones/Lennard_Jones_" + to_string(N_) + "_" + to_string(iteraciones_) + "_" + to_string((int)Lx_) + "_" + to_string((int)v_max_) + ".csv", guardar_datos),
velocidades("Datos/Velocidades/Velocidades_" + to_string(N) + "_" + to_string(iteraciones) + "_" + to_string((int)Lx_) + "_" + to_string((int)v_max_) + "_" + to_string((int)empuje_) + ".csv", guardar_datos),
energias("Datos/Energias/Energias_" + to_string(N) + "_" + to_string(iteraciones) + "_" + to_string((int)Lx_) + "_" + to_string((int)v_max_) + ".csv", guardar_datos),
var_posicion("Datos/Variacion_posicion/Var_posicion_" + to_string(N_) + "_" + to_string(iteraciones_) + "_" + to_string((int)Lx_) + "_" + to_string((int)v_max_) + ".csv", guardar_datos),
separacion("Datos/Variacion_posicion/Separacion_" + to_string(N_) + "_" + to_string(iteraciones_) + "_" + to_string((int)Lx_) + "_" + to_string((int)v_max_) + ".csv", guardar_datos)
{   
    // Elegimos el número de threads de la simulación
    omp_set_num_threads(n_threads);
    
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
    
    // Escribimos la cabecera de la separacion entre partículas
    separacion.escribir_cabecera("r2 - r1");

    // Escribimos los parámetros de la simulación 
    // para no tener que volver a escribirlos en el programa de python
    Archivo parametros("Datos/Parametros.txt", save_to_file);
    parametros.escribir("N\t" + to_string(N) + "\n");
    parametros.escribir("Lx\t" + to_string((int)Lx) + "\n");
    parametros.escribir("Ly\t" + to_string((int)Ly) + "\n");
    parametros.escribir("Pasos\t" + to_string(iteraciones) + "\n");
    parametros.escribir("vmax\t" + to_string(v_max) + "\n");
}

Simulador::~Simulador()
{
}


/**
 * @brief Inicializa el sistema con una configuración concreta
 * @param Modo Configuración de posición en la que se inicializa el sistema
 * - Modo::Azar
 * - Modo::Cuadricula
 * - Modo::CuadriculaAzar
 * - Modo::Panal
 * 
 */
void Simulador::inicializar(Estructura estructura, Modo modo)
{
    double x, y;
    double vx, vy;
    #pragma omp parallel
    {
        // Generador de números aleatorios por hilo (thread-safe)
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis_ang(0.0, 2 * PI);
        std::uniform_real_distribution<> dis_pos(0.0, 1.0);
        std::uniform_real_distribution<> dis_vx(0.0, v_max);

        #pragma omp for
        for (size_t i = 0; i < particulas.size(); i++)
        {
            double randnum = dis_ang(gen);
            double randvx = dis_vx(gen);
            double posx = dis_pos(gen) - 0.5;
            double posy = dis_pos(gen) - 0.5;

            switch (estructura)
            {
            case Estructura::Azar:
                x = ((double) rand() / RAND_MAX) * Lx;
                y = ((double) rand() / RAND_MAX) * Ly;
                break;

            case Estructura::Cuadricula:
                x = Lx / sqrt(N) * (0.5 + i % (int)sqrt(N));
                y = Ly / sqrt(N) * (0.5 + i / (int)sqrt(N));
                break;

            case Estructura::CuadriculaAzar:
                x = Lx / sqrt(N) * (0.5 + posx * 0.5 + i % (int)sqrt(N));
                y = Ly / sqrt(N) * (0.5 + posy * 0.5 + i / (int)sqrt(N));
                break;

            case Estructura::Panal:
                int n;
                n = int(1 + sqrt(1 + 2 * N)) / 2;
                double factor;
                factor = sqrt(3);
                y = (i / n / 2.0 + 0.25);
                x = (i % n + 1) + (i % n + 1) / 2;
                if (i / n % 2 == 1)
                x = (i % n) + (i % n) / 2 + 1.5;
                x = x / factor;
                break;

            default:
                break;
            }

            switch (modo)
            {
            case Modo::Compensado:
                // Módulo constante, ángulo aleatorio
                vx = v_max * cos(randnum);
                vy = v_max * sin(randnum);
                break;
                
            case Modo::Descompensado:
                vx = v_max * cos(randnum) + empuje;
                vy = v_max * sin(randnum);
                break;

            case Modo::Horizontal:
                vx =  randvx;
                vy = 0;
                break;
            
            default:
                break;
            }

            particulas[i].set_xy(x, y);
            posicion_inicial[i].set_xy(x, y);
            particulas[i].set_vxvy(vx, vy);
        }

    }
    
}

void Simulador::aplicar_Verlet(int n)
{
        // Calculamos la aceleración de todas las partículas para la iteracion inicial
        if (n == 1)
        {
            aceleraciones = calcular_aceleraciones();
        }

        double presion_por_paso = 0.0;

        // Actualizamos las posiciones de todas las partículas usando la imagen mínima y calculamos la presión
        #pragma omp parallel for reduction(+:presion_por_paso)
        for (size_t i = 0; i < N; i++)
        {
            particulas[i].actualizar_posicion(aceleraciones[i], paso);
            auto r = particulas[i].get_xy();
            auto v = particulas[i].get_vxvy();
            if (r[0] >= Lx)     { r[0] -= Lx; presion += 2 * fabs(v[0]); }
            else if (r[0] < 0) { r[0] += Lx; presion += 2 * fabs(v[0]); }
            if (r[1] >= Ly)     { r[1] -= Ly; presion += 2 * fabs(v[1]); }
            else if (r[1] < 0) { r[1] += Ly; presion += 2 * fabs(v[1]); }

            particulas[i].set_xy(r[0], r[1]);
        }
        presion += presion_por_paso;

        // Calculamos las nuevas aceleraciones de las partículas
        nuevas_aceleraciones = calcular_aceleraciones();

        // Calculamos las nuevas velocidades de las partículas
        #pragma omp for
        for (size_t i = 0; i < particulas.size(); i++)
        {
            particulas[i].actualizar_velocidad(aceleraciones[i], nuevas_aceleraciones[i], paso);
        }

        // Sobreescribimos las aceleraciones con las nuevas calculadas
        aceleraciones = nuevas_aceleraciones;

}

double Simulador::calcular_var_posicion()
{
    double variacion_media = 0.0;

    // Calculamos la variación en la posición
    #pragma omp parallel for reduction(+:variacion_media)
    for (int i = 0; i < N; i++)
    {
        double x0 = posicion_inicial[i].x;
        double y0 = posicion_inicial[i].y;
        double x = particulas[i].x;
        double y = particulas[i].y;
        double dx = x - x0;
        double dy = y - y0;

        // Aplicar condición de mínima imagen
        if (dx >  0.5 * Lx) dx -= Lx;
        if (dx < -0.5 * Lx) dx += Lx;        
        if (dy >  0.5 * Ly) dy -= Ly;
        if (dy < -0.5 * Ly) dy += Ly;
        
        variacion_media += dx * dx + dy * dy;
    }
    
    variacion_media = variacion_media / N;
    return variacion_media;
}

double Simulador::calcular_separacion(int id_1, int id_2)
{
    double separacion = 0.0;

            double xi = particulas[id_1].x;
            double yi = particulas[id_1].y;
            double xj = particulas[id_2].x;
            double yj = particulas[id_2].y;
            double dx = xi - xj;
            double dy = yi - yj;

            // Aplicar condición de mínima imagen
            if (dx >  0.5 * Lx) dx -= Lx;
            if (dx < -0.5 * Lx) dx += Lx;        
            if (dy >  0.5 * Ly) dy -= Ly;
            if (dy < -0.5 * Ly) dy += Ly;
            
            separacion += dx * dx + dy * dy;
    ;
    return separacion;
}

vector<array<double, 2>> Simulador::calcular_aceleraciones() const{

    vector<array<double, 2>> aceleraciones(N, {0.0, 0.0});

    #pragma omp parallel
    {
        #pragma omp for
        for (size_t i = 0; i < particulas.size(); i++)
        {
            auto ri = particulas[i].get_xy();
            // Calculamos la aceleración que sufre el cuerpo i debido a todos los cuerpos j != i a partir de i+1
            // ya que la aceleración que sufre i debido a j es igual en módulo y de sentido opuesto a la que sufre j
            // debido a i 
            for (size_t j = i+1; j < particulas.size(); j++)
            {
                auto rj = particulas[j].get_xy();
                array<double, 2> rij = {ri[0] - rj[0], ri[1] - rj[1]};

                // aplicar condición de mínima imagen
                    array<double, 2> L = {Lx, Ly};
                    for (size_t a = 0; a < 2; ++a) {
                        if (rij[a] >  0.5 * L[a]) rij[a] -= L[a];
                        if (rij[a] < -0.5 * L[a]) rij[a] += L[a];
                    }

                double rij2 = rij[0] * rij[0] + rij[1] * rij[1];

                // Comprobamos si la distancia entre partículas es pequeña
                // Si r^2 < 6.0 se calcula la aceleración
                // Si r^2 > 6.0 no se suma nada
                if (rij2 < 6.0)
                {
                    double rij4 = rij2 * rij2;
                    double rij8 = rij4 * rij4;
                    double rij14 = rij8 * rij4 * rij2;
                    double fuerza = 48 / rij14 - 24 / rij8;
                    for (int a = 0; a < 2; ++a) 
                    {
                        aceleraciones[i][a] += fuerza * rij[a];
                        aceleraciones[j][a] -= fuerza * rij[a];
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

    #pragma omp for
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

    #pragma omp for
    for (size_t i = 0; i < N; i++)
    {
        auto ri = particulas[i].get_xy();
        for (size_t j = i+1; j < N; j++)
        {
            {
                auto rj = particulas[j].get_xy();
                array<double, 2> rij = {ri[0] - rj[0], ri[1] - rj[1]};
                // aplicar condición de mínima imagen
                array<double, 2> L = {Lx, Ly};
                for (size_t a = 0; a < 2; ++a) {
                    if (rij[a] >  0.5 * L[a]) rij[a] -= L[a];
                    if (rij[a] < -0.5 * L[a]) rij[a] += L[a];
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

double Simulador::get_temperatura()
{
    double v_cuadrado = 0.0;
    double T;
    for (int i = 0; i < N; i++)
    {
        auto v = particulas[i].get_vxvy();
        v_cuadrado += v[0] * v[0] + v[1] * v[1];
    }

    T = v_cuadrado / (2 * N);
    return T;
}