
cmapdir=~/Applications/brainbits/f.CoactivationMap/src
ramonesdir=~/Applications/brainbits/RAMONES

rm cmtool.o coactivation.o Analyze.o
gcc -Wall -g -c cmtool.c $cmapdir/coactivation.c $ramonesdir/Analyze.c -I $cmapdir -I $ramonesdir
gcc -Wall cmtool.o coactivation.o Analyze.o -o cmtool
