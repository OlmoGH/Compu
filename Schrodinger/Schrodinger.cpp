# include <iostream>
# include <vector>
# include <complex>
# include <cmath>
# include <string>
# include <fstream>
# include <omp.h>

using namespace std;
using Vcomplex = vector<complex<double>>;
constexpr double PI = 3.14159265358979323846;

/**
 * @brief Clase que lleva a cabo la simulación
 * 
 */
class Simulacion
{
private:
    int N;
    int pasos;
    Vcomplex FuncionDeOnda;
    vector<double> FuncionNormaCuadrado;
    int n_ciclos;
    double s;
    double h;
    double lambda;
    double mu;
    double sigma;
    double norma_derecha;
    Vcomplex alpha;
    Vcomplex beta;
    vector<double> V;
    const string nombre_archivo;
    ofstream norma_file;
    int inicio_potencial;
    int ancho_potencial;

public:
    Simulacion(int N_, int n_ciclos_, double lambda_, double L_, double sigma_, double mu_, double s_, int pasos_); 
    ~Simulacion();
    void IniciarFuncion();
    void CalcularV(int inicio, int ancho);
    void CalcularAlpha();
    void CalcularBeta();
    void AplicarPaso();
    double get_norma(Vcomplex& funcion);
    void CalcularNorma();
    void CalcularNormaDerecha(int inicio_detector, int ancho_detector);
    void ExportarDatos(ofstream& archivo);

};

    /**
     * @brief Construct a new Simulacion object
     * 
     * @param N_ Número de particiones del espacio
     * @param n_ciclos_ Número de onda
     * @param lambda_ Altura del potencial
     * @param L_ Ancho de la caja
     * @param sigma_ Anchura de la gaussiana
     * @param mu_ Media de la gaussiana
     * @param s_ Paso temporal
     * @param pasos_ Número de pasos
     * @param nombre_archivo_ Nombre del archivo donde se escriben los datos
     */
    
    Simulacion::Simulacion(int N_, int n_ciclos_, double lambda_, double L_, double sigma_, double mu_, double s_, int pasos_) 
        : N(N_), 
        n_ciclos(n_ciclos_),
        pasos(pasos_), 
        alpha(N_-1), 
        beta(N_-1), 
        FuncionDeOnda(N_, {0.0, 0.0}), 
        FuncionNormaCuadrado(N_, 0.0),
        h(L_ / N_),
        s(s_ * N_ * N_ /(L_ * L_)), 
        V(N_, 0.0), 
        lambda(lambda_),
        sigma(sigma_),
        mu(mu_), 
        norma_file()
        {
            // cout << "N = " << N << endl;
            // cout << "n_ciclos = " << n_ciclos << endl;
            // cout << "lambda = " << lambda << endl;
            // cout << "h = " << h << endl;
            // cout << "sigma = " << sigma << endl;
            // cout << "mu = " << mu << endl;
            // cout << "s = " << s << endl;
            // cout << "pasos = " << pasos << endl;
            norma_file.open("Norma.txt");
            if (!norma_file)
            {
                cout << "No se pudo abrir el archivo\n";
            }
        }

    Simulacion::~Simulacion() 
    {
        norma_file.close();
    }
    
    /**
     * @brief Inicializa la función de onda a una gaussiana con bordes nulos
     * 
     */
    void Simulacion::IniciarFuncion()
    {
        // Inicializamos la función de onda compleja
        #pragma omp parallel for
        for (int j = 0; j < N; j++)
        {
            FuncionDeOnda[j] = polar(1.0, 2 * PI * n_ciclos * j / N) * exp(-(j * h - mu) * (j * h - mu) / (2 * sigma * sigma));
        }
        FuncionDeOnda[0] = {0.0, 0.0};
        FuncionDeOnda[N-1] = {0.0, 0.0};

        // Dividimos entre la norma para normalizar la función
        double norma = get_norma(FuncionDeOnda);
        for (int j = 0; j < N; j++)
        {
            FuncionDeOnda[j] = FuncionDeOnda[j] / sqrt(norma);
        }
    }
    
    /**
     * @brief Lee y guarda el valor del potencial
     * 
     */
    void Simulacion::CalcularV(int inicio, int ancho)
    {
        inicio_potencial = inicio;
        ancho_potencial = ancho;
        ofstream V_file;
        V_file.open("V.txt");
        double k0 = (2.0 * PI * n_ciclos) / N;

        #pragma omp parallel for
        for (int j = inicio_potencial; j < inicio_potencial +  ancho_potencial + 1; j++)
        {
            V[j] = lambda * k0 * k0;
        }


        for (int j = 0; j < N; j++)
        {
            V_file << V[j] / (k0 * k0);
            if (j < N-1)
            V_file << endl;
        }
        V_file.close();
    }

    /**
     * @brief Función que calcula el valor de alpha al inicio de la simuación
     * 
     */
    void Simulacion::CalcularAlpha()
    {
        alpha[N-2] = {0.0, 0.0};
        for (int j = N - 2; j > 0; j--)
        {
            complex<double> denominador(-2 - V[j], 2.0 / s);
            denominador += alpha[j];
            alpha[j-1] = -1.0 / denominador;
        }        
    }

    /**
     * @brief Función que calcula beta en cada iteración
     * 
     */
    void Simulacion::CalcularBeta()
    {
        beta[N-2] = {0.0, 0.0};
        for (int j = N-2; j > 0; j--)
        {
            complex<double> numerador(0.0, 4);
            numerador *= FuncionDeOnda[j] / s;
            numerador += -beta[j];
            complex<double> denominador(-2 - V[j], 2.0 / s);
            denominador += alpha[j];
            beta[j-1] = numerador / denominador;
        }
    }

    /**
     * @brief Aplica un paso temporal a la función de onda
     * 
     * @param n número del paso
     */
    void Simulacion::AplicarPaso()
    {
        Vcomplex chi(N, {0.0, 0.0});
        for (int j = 0; j < N-1; j++)
        {
            chi[j+1] = alpha[j] * chi[j] + beta[j];
        }

        #pragma omp parallel for
        for (int j = 1; j < N-1; j++)
        {
            FuncionDeOnda[j] = chi[j] - FuncionDeOnda[j];
            FuncionNormaCuadrado[j] = norm(FuncionDeOnda[j]);
        }
        
    }
    
    /**
     * @brief Get the norma object
     * 
     * @param funcion 
     * @return double 
     */
    double Simulacion::get_norma(Vcomplex& funcion)
    {
        int l = funcion.size();
        double norma = 0.0;
        #pragma omp parallel for reduction(+:norma)
        for (int i = 0; i < l; i++)
        {
            double mod_cuadrado = funcion[i].real() * funcion[i].real() + funcion[i].imag() * funcion[i].imag();
            norma += mod_cuadrado * h;
        }
        return norma;
    }

    /**
     * @brief Calcula la norma de la función para un paso y la guarda en un archivo
     * 
     * @param n Número del paso
     */
    void Simulacion::CalcularNorma()
     {
        double norma = 0.0;
        for (int j = 0; j < N; j++)
        {
            norma += FuncionNormaCuadrado[j] * h;
        }

        norma_file << norma;
     }

    /**
     * @brief Calcula y escribe la norma de la función a la derecha del potencial para un cierto paso
     * 
     * @param n Número del paso
     */
    void Simulacion::CalcularNormaDerecha(int inicio_detector, int ancho_detector)
     {
        norma_derecha = 0.0;

        for (int j = inicio_detector; j < inicio_detector + ancho_detector; j++)
        {
            norma_derecha += FuncionNormaCuadrado[j] * h;
        }

        norma_file << norma_derecha <<  "\t";
     }
    
     /**
     * @brief Guarda los datos del módulo de la función de onda en un archivo
     * 
     * @param matriz 
     */
    void Simulacion::ExportarDatos(ofstream& archivo) 
    {
        
        for (int j = 0; j < N; j++)
        archivo << FuncionNormaCuadrado[j] << '\t';
        archivo <<  endl;
    }

int main(int argc, char* argv[])
{
    ofstream archivo("Estados.txt");
    if (!archivo) {
        cout << "No se pudo abrir el archivo Estados.txt\n";
        return 1;
    }
    
    int N = 1000;
    int n_ciclos = N / 10;
    double lambda = 0.8;
    double L = 1;
    double sigma = L / 16;
    double mu = L / 4;
    double s = L * L / (8 * PI * N * n_ciclos);
    int pasos = 2000;
    int inicio = 4 * N / 5;
    int ancho = N / 5;

        for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-N" && i + 1 < argc) {
            N = std::stoi(argv[++i]);
        }
        else if ((arg == "-steps" || arg == "-s") && i + 1 < argc) {
            pasos = std::stoi(argv[++i]);
        }
    }

    Simulacion simulacion(N, n_ciclos, lambda, L, sigma, mu, s, pasos);

    simulacion.IniciarFuncion();
    simulacion.CalcularV(2 * N / 5, N / 5);

    simulacion.CalcularAlpha();
    for (int n = 0; n < pasos; n++)
    {
        simulacion.CalcularBeta();
        simulacion.AplicarPaso();
        // simulacion.CalcularNormaDerecha(n, inicio, ancho);
        simulacion.CalcularNorma();
        simulacion.ExportarDatos(archivo);
    }
    
    archivo.close();


    // Archivo con los datos de la simulación empleados para crear la animación
    // ofstream output("datos.txt");
    // if(!output)
    // {
    //     cerr << "No se pudo abrir el archivo\n";
    //     return 1;
    // }
    // output << N << endl;
    // output << n_ciclos << endl;
    // output << lambda << endl;
    // output << L << endl;
    // output << sigma << endl;
    // output << mu << endl;
    // output << s << endl;
    // output << pasos << endl;
    // output << inicio << endl;
    // output << ancho << endl;

    // output.close();
    
    return 0;
}