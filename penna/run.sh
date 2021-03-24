#!/bin/bash

cmake -S./ -Bbuild
make --no-print-directory -Cbuild
make --no-print-directory -Cbuild test
./build/main
python plot.py