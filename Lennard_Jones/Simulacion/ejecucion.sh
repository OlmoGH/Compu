g++ -I./lennard-jones-lib/include -L./lennard-jones-lib/build -fopenmp Lennard_Jones_sin_paralelo.cpp -llennardjones -o simulacion.exe
time ./simulacion.exe