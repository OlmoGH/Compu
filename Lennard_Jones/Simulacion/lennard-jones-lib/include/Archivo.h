# ifndef ARCHIVO_H
# define ARCHIVO_H
# include <string>
# include <fstream>
# include <vector> 
# include "Particula.h"

//***Clase Archicvo***/
class Archivo
{
private:
    std::string nombre;
    std::ofstream archivo;
    bool habilitado;
public:
    Archivo(const std::string& nombre_, bool habilitado_);
    ~Archivo();
    void escribir_cabecera(const std::string& cabecera);
    void escribir_posiciones(const std::vector<Particula>& p);
    void escribir_velocidades(const std::vector<Particula>& p);
    void escribir_energia(const std::vector<double>& K, const std::vector<double>& V);
    void escribir(const std::string& texto);
};

# endif