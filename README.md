# dnapunctuation

The DNA Punctuation program searches for stem-loop structures in a DNA sequence

Danila Zaev, Anton Zaikin, Maria Poptsova

Faculty of Computer Science, National Research University Higher School of Economics, Moscow, Russia

------------------------------------------------------------------------------------------------------

COMPILLATION:

You can compile the program as easy as
>g++ dnapunctuation.cpp -o dnapunctuation 

To accelerate processing time one can try these parameters (essential when annotating large chromosomes):
>g++ -O2 -march=native dnapunctuation.cpp -o dnapunctuation

Here, the executable version was compiled under Mac OS High Sierra, version 10.13.6 

------------------------------------------------------------------------------------------------------

RUNNING:

The program takes 7 arguments:

1. imput sequence file in fasta format
2. Minimum size of a stem
3. Maximum size of a stem
4. Minimum size of a loop
5. Maximum size of a loop
6. Number of mismatches allowed

One should run the program as
>./dnapuncutation L1.fna 6 15 3 10 1

where 
L1.fna is a sequence file containing L1 human transposon
6 is the minimum size of a stem
10 is the maximum size of a stem
4 is the minimum size of a loop
5 is the maximum size of a loop
1 is the number of gaps/mismatches allowed

The promram output is the text tab-delimitted file named according to the stem-loop parameters specified as input: L1.fna.S6-15_L3-10_M1.pal

L1.fna and L1.fna.S6-15_L3-10_M1.pal can be found in the main directory





