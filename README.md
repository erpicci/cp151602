# Progetto Calcolo Parallelo 2015/2016
Calcolo degli autovalori di una matrice simmetrica tridiagonale, ponendo
a confronto diversi algoritmi seriali e paralleli (basati su MPI).

## Compilazione
Compila i sorgenti con impostazioni di default (openMPI)

    make

Compila i sorgenti su Power7

    make cc=mpcc

Sposta gli eseguibili generati sotto la cartella bin/

    make install

Genera la documentazione con Doxygen in doc/html/

    make doc


## Utilizzo
Genera (su standard output) una matrice di dimensione 10

    ./generator -n 10


Calcola gli autovalori della matrice memorizzata in file.dat (esempio usando il risolutore di lapack)

    ./lapack < file.dat
