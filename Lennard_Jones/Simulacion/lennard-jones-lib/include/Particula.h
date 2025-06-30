# ifndef PARTICULA_H
# define PARTICULA_H
# include <array>

//***Clase Particula***//
class Particula
{
public:
    double x;
    double y;
    double vx;
    double vy;
    Particula(double x_, double y_);
    Particula();
    void set_xy(double x_, double y_);
    void set_vxvy(double vx_, double vy_);
    std::array<double, 2> get_xy() const;
    std::array<double, 2> get_vxvy() const;
    void actualizar_posicion(const std::array<double, 2>& aceleracion, double paso);
    void actualizar_velocidad(const std::array<double, 2>& aceleracion, const std::array<double, 2>& aceleracion_nueva, double paso);
};

# endif