#!/bin/bash

echo " ===========> RUNNING CMAKE <=========== "
cmake -Bbuild
echo " ===========> RUNNING  MAKE <=========== "
make -Cbuild --no-print-directory
echo " ===========> RUNNING  MAIN <=========== "
./build/main
#echo " ===========> RUNNING  PLOT <=========== "
#python plot.py
