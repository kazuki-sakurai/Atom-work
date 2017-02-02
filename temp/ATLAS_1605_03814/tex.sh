#!/bin/bash

runname=$1

if [[ runname == "" ]]; then
    echo "[runname]"
    exit
fi

cd tex
pdflatex $runname.tex
open $runname.pdf
cd ../

exit