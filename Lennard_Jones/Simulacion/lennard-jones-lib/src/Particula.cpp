# include "Particula.h"
# include <array>

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

void Particula::actualizar_posicion(const std::array<double, 2>& aceleracion, double paso){
    x += paso * vx + 0.5 * paso * paso * aceleracion[0];
    y += paso * vy + 0.5 * paso * paso * aceleracion[1];
}

void Particula::actualizar_velocidad(const std::array<double, 2>& aceleracion, const std::array<double, 2>& aceleracion_nueva, double paso){
    vx += 0.5 * paso * (aceleracion[0] + aceleracion_nueva[0]);
    vy += 0.5 * paso * (aceleracion[1] + aceleracion_nueva[1]);
}

Particula:: Particula(double x_, double y_) : x(x_), y(y_), vx(0.0), vy(0.0)
{
}

Particula:: Particula() : x(0.0), y(0.0), vx(0.0), vy(0.0)
{
}