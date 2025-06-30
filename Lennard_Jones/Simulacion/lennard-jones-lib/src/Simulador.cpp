# include "Simulador.h"
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

Simulador::Simulador(const int N_, const int iteraciones_, const double paso_, const double L_, double v_max_, bool guardar_datos) : 
N(N_), 
iteraciones(iteraciones_), 
paso(paso_), 
L(L_),
v_max(v_max_), 
particulas(N_), 
aceleraciones(N_, {0.0, 0.0}), 
nuevas_aceleraciones(N_, {0.0, 0.0}),
presion(0.0), 
energia_cinetica(iteraciones_), 
energia_potencial(iteraciones_),
posicion_inicial(N_),
save_to_file(guardar_datos),
datos("Datos/Posiciones/Lennard_Jones_" + to_string(N_) + "_" + to_string(iteraciones_) + "_" + to_string((int)L_) + ".csv", guardar_datos),
velocidades("Datos/Velocidades/Velocidades_" + to_string(N) + "_" + to_string(iteraciones) + "_" + to_string((int)L) + ".csv", guardar_datos),
energias("Datos/Energias/Energias_" + to_string(N) + "_" + to_string(iteraciones) + "_" + to_string((int)L) + ".csv", guardar_datos),
var_posicion("Datos/Variacion_posicion/Var_posicion_" + to_string(N_) + "_" + to_string(iteraciones_) + "_" + to_string((int)L_) + ".csv", guardar_datos),
separacion("Datos/Variacion_posicion/Separacion_" + to_string(N_) + "_" + to_string(iteraciones_) + "_" + to_string((int)L_) + ".csv", guardar_datos)
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
    
    // Escribimos la cabecera de la separacion entre partículas
    separacion.escribir_cabecera("r2 - r1");

    // Escribimos los parámetros de la simulación 
    // para no tener que volver a escribirlos en el programa de python
    Archivo parametros("Datos/Parametros.txt", save_to_file);
    parametros.escribir("N\t" + to_string(N) + "\n");
    parametros.escribir("L\t" + to_string((int)L) + "\n");
    parametros.escribir("Pasos\t" + to_string(iteraciones) + "\n");
    parametros.escribir("vmax\t" + to_string(v_max) + "\n");
}

Simulador::~Simulador()
{
}

void Simulador::inicializar()
{
    #pragma omp parallel
    {
        // Generador de números aleatorios por hilo (thread-safe)
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis_ang(0.0, 2 * PI);

        #pragma omp for
        for (size_t i = 0; i < particulas.size(); i++)
        {
            double randnum = dis_ang(gen);

            // Inicio al azar
            // double x = ((double) rand() / RAND_MAX) * L;
            // double y = ((double) rand() / RAND_MAX) * L;

            // Inicio dentro de cuadriculas en una posicion al azar 
            //posicion_inicial[i][0].set_xy(L / sqrt(N) * (posx + i % ((int) sqrt(N))), L / sqrt(N) * (posy + i / (int) sqrt(N)));

            // Inicio cuadriculado
            double x = L / sqrt(N) * (0.5 + i % (int)sqrt(N));
            double y = L / sqrt(N) * (0.5 + i / (int)sqrt(N));

            // Inicio en panal
            //posicion_inicial[i].set_xy(L / sqrt(N) * (0.5 + i % ((int) sqrt(N))), L / sqrt(N) * (0.5 + i / (int) sqrt(N)));
            //if (i % (2 * (int)L) > (int)L)
            //posicion_inicial[i].set_xy(L / sqrt(N) * (1 + i % ((int) sqrt(N))), L / sqrt(N) * (1 + i / (int) sqrt(N)));

            // Velocidades radiales
            //posicion_inicial[i].set_vxvy(2 * (particulas[i].get_xy()[0] - L/2), 2 * (particulas[i].get_xy()[1] - L/2));

            // Velocidades aleatorias
            double vx = v_max * cos(randnum);
            double vy = v_max * sin(randnum);

            particulas[i].set_xy(x, y);
            posicion_inicial[i].set_xy(x, y);
            particulas[i].set_vxvy(vx, vy);
        }

    }
    
}

void Simulador::aplicar_Verlet()
{
    #pragma omp parallel
    {
        // Calculamos la aceleración de todas las partículas
        #pragma omp single
        aceleraciones = calcular_aceleraciones();
        double presion_por_paso = 0.0;

        // Actualizamos las posiciones de todas las partículas usando la imagen mínima y calculamos la presión
        #pragma omp for reduction(+:presion_por_paso)
        for (size_t i = 0; i < N; i++)
        {
            particulas[i].actualizar_posicion(aceleraciones[i], paso);
            auto r = particulas[i].get_xy();
            auto v = particulas[i].get_vxvy();
            if (r[0] >= L)     { r[0] -= L; presion += 2 * fabs(v[0]); }
            else if (r[0] < 0) { r[0] += L; presion += 2 * fabs(v[0]); }
            if (r[1] >= L)     { r[1] -= L; presion += 2 * fabs(v[1]); }
            else if (r[1] < 0) { r[1] += L; presion += 2 * fabs(v[1]); }

            particulas[i].set_xy(r[0], r[1]);
        }
        presion += presion_por_paso;

        // Calculamos las nuevas aceleraciones de las partículas
        #pragma omp single
        nuevas_aceleraciones = calcular_aceleraciones();

        // Calculamos las nuevas velocidades de las partículas
        #pragma omp for
        for (size_t i = 0; i < particulas.size(); i++)
        {
            particulas[i].actualizar_velocidad(aceleraciones[i], nuevas_aceleraciones[i], paso);
        }
    }

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
        if (dx >  0.5 * L) dx -= L;
        if (dx < -0.5 * L) dx += L;        
        if (dy >  0.5 * L) dy -= L;
        if (dy < -0.5 * L) dy += L;
        
        variacion_media += dx * dx + dy * dy;
    }
    
    variacion_media = variacion_media / N;
    return variacion_media;
}

double Simulador::calcular_separacion()
{
    double separacion = 0.0;
    // #pragma omp parallel for reduction(+:separacion)
    // for (size_t i = 0; i < particulas.size(); i++)
    // {
    //     for (size_t j = i+1; j < particulas.size(); j++)
    //     {
    //         double xi = particulas[i].x;
    //         double yi = particulas[i].y;
    //         double xj = particulas[j].x;
    //         double yj = particulas[j].y;
    //         double dx = xi - xj;
    //         double dy = yi - yj;

    //         // Aplicar condición de mínima imagen
    //         if (dx >  0.5 * L) dx -= L;
    //         if (dx < -0.5 * L) dx += L;        
    //         if (dy >  0.5 * L) dy -= L;
    //         if (dy < -0.5 * L) dy += L;
            
    //         separacion += dx * dx + dy * dy;
    //     }
        
    // }

            double xi = particulas[N / 2].x;
            double yi = particulas[N / 2].y;
            double xj = particulas[N / 2 + 1].x;
            double yj = particulas[N / 2 + 1].y;
            double dx = xi - xj;
            double dy = yi - yj;

            // Aplicar condición de mínima imagen
            if (dx >  0.5 * L) dx -= L;
            if (dx < -0.5 * L) dx += L;        
            if (dy >  0.5 * L) dy -= L;
            if (dy < -0.5 * L) dy += L;
            
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
                for (size_t a = 0; a < 2; ++a) {
                    if (rij[a] >  0.5 * L) rij[a] -= L;
                    if (rij[a] < -0.5 * L) rij[a] += L;
                }
                double rij2 = rij[0] * rij[0] + rij[1] * rij[1];
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