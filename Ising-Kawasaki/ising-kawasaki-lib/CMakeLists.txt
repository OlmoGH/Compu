cmake_minimum_required(VERSION 3.10)
project(isingkawasaki)

# Incluye los headers
include_directories(include)

# Crea la librería estática
add_library(isingkawasaki STATIC

)

# Opcional: si quieres instalar los headers
install(DIRECTORY include/ DESTINATION include)