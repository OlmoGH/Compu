# include "Archivo.h"
# include <string>
# include <iomanip>
# include <vector>

Archivo::Archivo(const std::string& nombre_, bool habilitado_) : nombre(nombre_), archivo(), habilitado(habilitado_)
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

void Archivo::escribir_cabecera(const std::string& cabecera)
{
    if (habilitado)
    {
        archivo << cabecera << std::endl;
    }
    
}

void Archivo::escribir_posiciones(const std::vector<Particula>& p)
{
    if (habilitado)
    {
        for (int i = 0; i < p.size()-1; i++)
        {
            archivo << p[i].get_xy()[0] << "," << p[i].get_xy()[1] << ",";
        }

        archivo << p.back().get_xy()[0] << "," << p.back().get_xy()[1] << std::endl;
    }
}

void Archivo::escribir_velocidades(const std::vector<Particula>& p)
{
    if (habilitado)
    {
        for (int i = 0; i < p.size()-1; i++)
        {
            archivo << p[i].get_vxvy()[0] << "," << p[i].get_vxvy()[1] << ",";
        }

        archivo << p.back().get_vxvy()[0] << "," << p.back().get_vxvy()[1] << std::endl;
    }
}

void Archivo::escribir_energia(const std::vector<double>& K, const std::vector<double>& V)
{
    if (habilitado)
    {
        for (int i = 0; i < K.size(); i++)
        {
            archivo << K[i] << "," << V[i] << std::endl;
        }
        
    }
    
}

void Archivo::escribir(const std::string& texto)
{
    archivo << texto;
}