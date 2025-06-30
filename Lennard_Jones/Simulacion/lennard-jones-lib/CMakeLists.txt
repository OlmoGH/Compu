cmake_minimum_required(VERSION 3.10)
project(lennardjones)

# Incluye los headers
include_directories(include)

# Crea la librería estática
add_library(lennardjones STATIC
    src/Archivo.cpp
    src/Particula.cpp
    src/Simulador.cpp
)

# Opcional: si quieres instalar los headers
install(DIRECTORY include/ DESTINATION include)