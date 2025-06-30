# include <iostream>
# include <vector>
# include <complex>
# include <cmath>
# include <string>
# include <fstream>

using namespace std;
using Vcomplex = vector<complex<double>>;
constexpr double PI = 3.14159265358979323846;


class Simulacion
{
private:
    int N;
    int pasos;
    Vcomplex FuncionDeOnda;
    vector<vector<double>> FuncionNormaCuadrado;
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

public:
    Simulacion(int N_, int n_ciclos_, double lambda_, double L_, double sigma_, double mu_, double s_, int pasos_, const string nombre_archivo_) 
        : N(N_), 
        n_ciclos(n_ciclos_),
        pasos(pasos_), 
        alpha(N_-1), 
        beta(N_-1), 
        FuncionDeOnda(N_), 
        FuncionNormaCuadrado(N_, vector<double>(pasos_, 0.0)),
        h(L_ / N_),
        s(s_ * N_ * N_ /(L_ * L_)), 
        V(N_, 0.0), 
        lambda(lambda_),
        sigma(sigma_),
        mu(mu_), 
        nombre_archivo(nombre_archivo_) 
        {
            cout << "N = " << N << endl;
            cout << "n_ciclos = " << n_ciclos << endl;
            cout << "lambda = " << lambda << endl;
            cout << "h = " << h << endl;
            cout << "sigma = " << sigma << endl;
            cout << "mu = " << mu << endl;
            cout << "s = " << s << endl;
            cout << "pasos = " << pasos << endl;
            norma_file.open("Norma.txt");
            if (!norma_file)
            cout << "No se pudo abrir el archivo\n";
        }

    ~Simulacion() 
    {
        norma_file.close();
    }
    
    void IniciarFuncion()
    {
        for (int j = 0; j < N; j++)
        {
            FuncionDeOnda[j] = polar(1.0, 2 * PI * n_ciclos * j / N) * exp(-(j * h - mu) * (j * h - mu) / (2 * sigma * sigma));
        }
        FuncionDeOnda[0] = {0.0, 0.0};
        FuncionDeOnda[N-1] = {0.0, 0.0};
    }
    
    void CalcularV()
    {
        ofstream V_file;
        V_file.open("V.txt");
        double k0 = (2.0 * PI * n_ciclos) / N;
        int inicio = 2 * N / 5;
        int ancho = 1;
        for (int j = inicio; j < inicio +  ancho + 1; j++)
        {
            V[j] = lambda * k0 * k0;
        }
        for (int j = 0; j < N; j++)
        {
            V_file << V[j];
            if (j < N-1)
            V_file << endl;
        }
        V_file.close();
    }

    void CalcularAlpha()
    {
        alpha[N-2] = {0.0, 0.0};
        for (int j = N - 2; j > 0; j--)
        {
            complex<double> denominador(-2 - V[j], 2.0 / s);
            denominador += alpha[j];
            alpha[j-1] = -1.0 / denominador;
        }        
    }

    void CalcularBeta()
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

    void AplicarPaso(int n)
    {
        Vcomplex chi(N, {0.0, 0.0});
        for (int j = 0; j < N-1; j++)
        {
            chi[j+1] = alpha[j] * chi[j] + beta[j];
        }
        for (int j = 1; j < N-1; j++)
        {
            FuncionDeOnda[j] = chi[j] - FuncionDeOnda[j];
            FuncionNormaCuadrado[j][n] = norm(FuncionDeOnda[j]);
        }
        
    }

    void Simular()
    {
        IniciarFuncion();
        CalcularV();
        CalcularAlpha();
        for (int n = 0; n < pasos; n++)
        {
            CalcularBeta();
            AplicarPaso(n);
            CalcularNorma(n);
        }
        ExportarDatos(FuncionNormaCuadrado);
        
    }
    
    void CalcularNorma(int n)
     {
        double norma = 0.0;
        for (int j = 0; j < N; j++)
        {
            norma += FuncionNormaCuadrado[j][n] * h;
        }

        norma_file << norma;
        if (n < pasos-1)
        norma_file << endl;
     }

    void CalcularNormaDerecha(int n)
     {
        norma_derecha = 0.0;
        int inicio = 2 * N / 5;
        int ancho = 1;
        for (int j = inicio + ancho; j < N; j++)
        {
            norma_derecha += FuncionNormaCuadrado[j][n] * h;
        }
     }
    void ExportarDatos(vector<vector<double>> &matriz) {
    ofstream archivo(nombre_archivo);
    if (!archivo)
    {
        cerr << "No se pudo abrir el archivo\n";
        return;
    }
    for (int n = 0; n < pasos; n++) {
        for (int j = 0; j < N-1; j++) {
            archivo << matriz[j][n] << "\t";
        }
        archivo << matriz[N-1][n] << "\n";
    }
    archivo.close();
    }

};



int main()
{
    int N = 1000;
    int n_ciclos = N / 10;
    double lambda = 3;
    double L = 1;
    double sigma = L / 16;
    double mu = L / 4;
    double s = L * L / (8 * PI * N * n_ciclos);
    int pasos = 3000;
    string nombre = "Estados.txt";

    Simulacion simulacion(N, n_ciclos, lambda, L, sigma, mu, s, pasos, nombre);
    simulacion.Simular();

    ofstream output("datos.txt");
    if(!output)
    {
        cerr << "No se pudo abrir el archivo\n";
        return 1;
    }
    output << N << endl;
    output << n_ciclos << endl;
    output << lambda << endl;
    output << L << endl;
    output << sigma << endl;
    output << mu << endl;
    output << s << endl;
    output << pasos << endl;

    output.close();
    
    return 0;
}