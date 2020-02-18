#! /bin/bash
MAIN=liftforcemap_lagint.c
SUB=`ls *.c | fgrep -v $MAIN`
echo $MAIN $SUB
gcc -lm $MAIN $SUB
./a.out
./contours.plt
ps2pdf14 output.ps
rm output.ps
