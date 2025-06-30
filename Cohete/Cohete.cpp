// Cohete.cpp
#include <iostream>
#include <array>
#include <math.h>
#include <fstream>
#include <algorithm>

using namespace std;

class Cuerpo
{
public:
    double masa;
    double r;
    double phi;
    double pr;
    double pPhi;
    Cuerpo(double r_, double phi_, double pr_, double pPhi_);
    ~Cuerpo();
    void set_posicion(double r_, double phi_);
    array<double, 2> get_velocidad();
};

Cuerpo::Cuerpo(double r_, double phi_, double pr_, double pPhi_) : r(r_), phi(phi_), pr(pr_), pPhi(pPhi_)
{
}

Cuerpo::~Cuerpo()
{
}

void Cuerpo::set_posicion(double r_, double phi_)
{
    r = r_;
    phi = phi_;
}

array<double, 2> Cuerpo::get_velocidad()
{
    double vr = pr;
    double vTan = pPhi / r;
    return {vr, vTan};
}

class Simulador
{
public:
    Cuerpo cohete;
    double delta;
    double w;
    double mu;
    ofstream archivo;
    ofstream momentos;
    ofstream energia;
    Simulador(double r0_, double phi0_, double pr0_, double pPhi0_);
    ~Simulador();
    double rL(double r_, double phi_, double h_, double n_);
    double rPunto(double pr_);
    double prPunto(double r_, double phi_, double pPhi_, double h_, double n_);
    double phiPunto(double r_, double pPhi_);
    double pPhiPunto(double r_, double phi_, double h_, double n_);
    void aplicarRungeKutta(double n, double h);
    array<double, 3> calcularEnergia(double n, double h);
    double pasoVariable(double h0, double n, double errMax);
    
};

Simulador::Simulador(double r0_, double phi0_, double pr0_, double pPhi0_) : 
cohete(r0_, phi0_, pr0_, pPhi0_),
archivo("Posiciones.txt"),
momentos("Momentos.txt"),
energia("Energia.txt")
{
    mu = 0.01230246418;
    delta = 6.67 * 5.9736 / pow(3.844, 3) * 1.0E-11;
    w = 0.0000026617;
    archivo << "r\tphi" << endl;
    momentos << "pr\tPphi" << endl;
    energia << "T\tV\tH'" << endl;
}

Simulador::~Simulador()
{
}

double Simulador::rL(double r_, double phi_, double h_, double n_)
{
    double resultado = sqrt(1 + r_ * r_ - 2 * r_ * cos(phi_ - w * n_ * h_));
    return resultado;
}

double Simulador::rPunto(double pr_)
{
    return pr_;
}

double Simulador::prPunto(double r_, double phi_, double pPhi_, double h_, double n_)
{
    double rPrima = rL(r_, phi_, h_, n_);
    double termino1 = pPhi_ * pPhi_ / (r_ * r_ * r_);
    double termino2 = delta * (1.0 / (r_ * r_) + mu / (rPrima * rPrima * rPrima) * (r_ - cos(phi_ - w * n_ * h_)));
    return termino1 - termino2;
}

double Simulador::phiPunto(double r_, double pPhi_)
{
    return pPhi_ / (r_ * r_);
}

double Simulador::pPhiPunto(double r_, double phi_, double h_, double n_)
{
    double rPrima = rL(r_, phi_, h_, n_);
    double resultado = - delta * mu * r_ * sin(phi_ - w * h_ * n_) / (rPrima * rPrima * rPrima);
    return resultado;
}

void Simulador::aplicarRungeKutta(double n, double h)
{
    double k1[4], k2[4], k3[4], k4[4];

    // Asignamos los valores iniciales de las funciones
    double r0 = cohete.r;
    double phi0 = cohete.phi;
    double pr0 = cohete.pr;
    double pPhi0 = cohete.pPhi;

    // Calculamos el primer array de coeficientes 
    // k1[i] = h * f_i(r, phi, pr, pPhi, h)
    k1[0] = h * rPunto(pr0);
    k1[1] = h * phiPunto(r0, pPhi0);
    k1[2] = h * prPunto(r0, phi0, pPhi0, h, n);
    k1[3] = h * pPhiPunto(r0, phi0, h, n);    
    
    // Calculamos el segundo array de coeficientes 
    // k2[i] = h * f_i(r + k1[0]/2, phi + k1[1]/2, pr + k1[2]/2, pPhi + k1[3]/2, h/2)
    k2[0] = h * rPunto(pr0 + k1[2] / 2);
    k2[1] = h * phiPunto(r0 + k1[0] / 2, pPhi0 + k1[3] / 2);
    k2[2] = h * prPunto(r0 + k1[0] / 2, phi0 + k1[1] / 2, pPhi0 + k1[3] / 2, h / 2, n);
    k2[3] = h * pPhiPunto(r0 + k1[0] / 2, phi0 + k1[1] / 2, h / 2, n);       
    
    // Calculamos el tercer array de coeficientes 
    // k3[i] = h * f_i(r + k2[0]/2, phi + k2[1]/2, pr + k2[2]/2, pPhi + k2[3]/2, h/2)
    k3[0] = h * rPunto(pr0 + k2[2] / 2);
    k3[1] = h * phiPunto(r0 + k2[0] / 2, pPhi0 + k2[3] / 2);
    k3[2] = h * prPunto(r0 + k2[0] / 2, phi0 + k2[1] / 2, pPhi0 + k2[3] / 2, h / 2, n);
    k3[3] = h * pPhiPunto(r0 + k2[0] / 2, phi0 + k2[1] / 2, h / 2, n);    
    
    // Calculamos el cuarto array de coeficientes 
    // k4[i] = h * f_i(r + k4[0], phi + k4[1], pr + k4[2], pPhi + k4[4], h)
    k4[0] = h * rPunto(pr0 + k3[2]);
    k4[1] = h * phiPunto(r0 + k3[0], pPhi0 + k3[3]);
    k4[2] = h * prPunto(r0 + k3[0], phi0 + k3[1], pPhi0 + k3[3], h, n);
    k4[3] = h * pPhiPunto(r0 + k3[0], phi0 + k3[1], h, n); 
    
    
    // Calculamos el nuevo coeficiente de la aproximación
    double k[4];
    for (int i = 0; i < 4; i++)
    {
        k[i] = (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
    }

    // Calculamos los valores finales de las funciones
    cohete.r += k[0];
    cohete.phi += k[1];
    cohete.pr += k[2];
    cohete.pPhi += k[3];
    
}

array<double, 3> Simulador::calcularEnergia(double n, double h)
{
    auto v = cohete.get_velocidad();
    double v_mod_cuadrado = v[0] * v[0] + v[1] * v[1];
    double T = 0.5 * v_mod_cuadrado;
    double V = -delta * (1.0 / cohete.r + mu / rL(cohete.r, cohete.phi, h, n));

    return {T, V, T + V - w * cohete.pPhi};
}

double Simulador::pasoVariable(double h0, double n, double tolerancia)
{
    Cuerpo estado_original = cohete;

    // Aplicamos una iteración sobre cohete con paso h0
    aplicarRungeKutta(n, h0);
    Cuerpo y1 = cohete;

    // Devolvemos el sistema al estado original
    cohete = estado_original;

    // Aplicamos dos pasos con paso h0/2
    aplicarRungeKutta(n, h0/2);
    aplicarRungeKutta(n + 0.5, h0/2);
    Cuerpo y2 = cohete;

    // Calculamos el error entre los pasos
    double error_r = 15.0 * std::abs(y2.r - y1.r) / 16.0;
    double error_phi = 15.0 * std::abs(y2.phi - y1.phi) / 16.0;
    double error_pr = 15.0 * std::abs(y2.pr - y1.pr) / 16.0;
    double error_pPhi = 15.0 * std::abs(y2.pPhi - y1.pPhi) / 16.0;

    double max_error = max({error_r, error_phi, error_pr, error_pPhi});
    if ((int)n % 100 == 0)
    cout << n << " --> " << max_error << endl;

    double a = pow(max_error / tolerancia, 0.2);
    double b = 1E-8;
    double s = max(a, b);
    double h_max = h0 / s;

    // Si el paso es demasiado preciso se aumenta el tamaño del paso
    if (s < 1.0)
    return pasoVariable(2 * h0, n, tolerancia);

    // Si el paso no es lo suficientemente preciso se disminuye hasta el tamaño máximo del paso
    else if (s > 2.0)
    return pasoVariable(h_max, n, tolerancia);

    // Si el paso es correcto se acepta el segundo método en dos pasos
    else 
    return h0;

}
int main()
{
    int iteraciones = 30000;
    double h_inicial = 0.1;
    double tolerancia = 0.00001;

    // Introducimos las condiciones iniciales del cohete en unidades de SI
    double r0 = 6378160.0; // Radio de la tierra en metros
    double phi0 = -.95;
    double v0 = 11186.0 * 0.991; // Velocidad inicial del cohete, similar a la de escape
    double parteTangencial = 0.9; // Fracción de la velocidad en la dirección tangencial
    double parteNormal = sqrt(1 - parteTangencial * parteTangencial);

    // Calculamos los momentos conjugados asociados a dicha velocidad
    double pPhi0 = r0 * parteTangencial * v0;
    double pr0 = parteNormal * v0;

    // Reescalamos las unidades con la masa del cohete (hemos tomado como 1) y la distancia Tierra-Luna
    double dTL = 3.844E8;
    r0 = r0 / dTL;
    pr0 = pr0 / dTL;
    pPhi0 = pPhi0 / (dTL * dTL);

    // Declaramos la variable simulación
    Simulador simulacion(r0, phi0, pr0, pPhi0);

    double h_actual = h_inicial;
    auto energias = simulacion.calcularEnergia(0, h_actual);

    for (int n = 0; n < iteraciones; n++)
    {
        // Escribimos los datos de la simulación en archivos de salida
        simulacion.archivo << simulacion.cohete.r << "\t" << simulacion.cohete.phi << endl;
        // simulacion.momentos << simulacion.cohete.pr << "\t" << simulacion.cohete.pPhi << endl;
        simulacion.energia << energias[0] << "\t" << energias[1] << "\t" << energias[2] << endl;

        h_actual = simulacion.pasoVariable(h_actual, n, tolerancia);
    }

    
    return 0;
}