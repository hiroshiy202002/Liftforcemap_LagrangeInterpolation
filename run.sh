#! /bin/bash
MAIN=liftforcemap_lagint.c
SUB=`ls *.c | fgrep -v $MAIN`
echo $MAIN $SUB
gcc -lm $MAIN $SUB
./a.out
./colormap.py IntpltdFx.dat ; mv output.ps output1.ps
./colormap.py IntpltdFr.dat ; mv output.ps output2.ps
./colormap.py IntpltdFt.dat ; mv output.ps output3.ps
./nullclines.py ; mv output.ps output4.ps
cat output*.ps > output.ps
ps2pdf14 output.ps
rm output*.ps
