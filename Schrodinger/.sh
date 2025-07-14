g++ -fopenmp Schrodinger.cpp -o SchrodingerO0.exe -static -lm -O0
g++ -fopenmp Schrodinger.cpp -o SchrodingerO1.exe -static -lm -O1
g++ -fopenmp Schrodinger.cpp -o SchrodingerO2.exe -static -lm -O2
g++ -fopenmp Schrodinger.cpp -o SchrodingerO3.exe -static -lm -O3

hyperfine \
  --prepare "del /f /q Norma.txt Estados.txt V.txt 2> nul" \
  --warmup 2 \
  --runs 3 \
  --export-json results2.json \
  -L opt_level 0,1,2,3 \
  -L threads 1,2,4,8 \
  -L N 1000,10000,100000 \
  -L steps 1000,10000,100000 \
  -n "O{opt_level}-Threads{threads}-N{N}-Steps{steps}" \
  "set OMP_NUM_THREADS={threads} && SchrodingerO{opt_level}.exe -N {N} -steps {steps}"