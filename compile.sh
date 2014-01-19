rm cmtool.o coactivation.o Analyze.o
gcc -Wall -g -c cmtool.c coactivation.c Analyze.c
gcc -Wall cmtool.o coactivation.o Analyze.o -o cmtool
